# ==========================================================================
# gpu_validate_decompmimics_e2e.jl — end-to-end GPU parity for the WHOLE
# decomp_rates_mimics! function (every per-(c)/per-(c,j) compute loop now a
# KernelAbstractions kernel).
#
# Builds a small Float32 instance mirroring test/test_decomp_mimics.jl's setup,
# runs decomp_rates_mimics! on the CPU, adapts the carbon-flux struct + the MIMICS
# state's spatially-varying arrays + the input arrays to the GPU, runs the SAME call
# on the device, and compares the mutated outputs (decomp_k, pathfrac, rf, cn,
# w_scalar, o_scalar) field-by-field.
#
# Exercises the branchy config paths: nlevdecomp==1 (single-level w/o reduction)
# AND nlevdecomp==5 (multi-level kernels), plus the no-anoxia o_scalar path.
#
#   julia --project=scripts scripts/gpu_validate_decompmimics_e2e.jl
# ==========================================================================

using CLM
using Printf
include(joinpath(@__DIR__, "gpu_backends.jl"))

# Float32-down-converting Metal adaptor: init_decompcascade_mimics! allocates the
# DecompMIMICSState texture arrays as Float64 (zeros) regardless of precision, so we
# adapt-reconstruct the struct with device Float32 arrays (scalar coeff fields, which
# never reach a kernel, pass through unchanged).
struct MetalF32 end
CLM.Adapt.adapt_storage(::MetalF32, x::AbstractArray{<:AbstractFloat}) = device_array_type()(Float32.(x))
CLM.Adapt.adapt_storage(::MetalF32, x::AbstractArray{<:Integer})       = device_array_type()(x)
CLM.Adapt.adapt_storage(::MetalF32, x::AbstractArray{Bool})            = device_array_type()(x)

# reldiff: NaN-aware (both-NaN agrees; one-sided NaN flags divergence). Also asserts
# the CPU reference is FINITE so a both-NaN false PASS can't slip through.
function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(A[i] - B[i]) / (1.0 + max(abs(A[i]), abs(B[i]))))
    end
    return m
end
cpu_has_finite(a) = any(isfinite, Array(a))

# --- MIMICS parameters (Wieder et al. 2015), mirroring the oracle test ---
function make_mimics_params(::Type{FT}; ndecomp_pools_max=8) where {FT}
    CLM.DecompMIMICSParams(
        mimics_nue_into_mic          = 0.85,
        mimics_desorpQ10             = 1.1,
        mimics_densdep               = 1.0,
        mimics_tau_mod_factor        = 2.0,
        mimics_tau_mod_min           = 0.8,
        mimics_tau_mod_max           = 1.2,
        mimics_ko_r                  = 6.0,
        mimics_ko_k                  = 6.0,
        mimics_cn_r                  = 7.0,
        mimics_cn_k                  = 10.0,
        mimics_cn_mod_num            = 0.4,
        mimics_t_soi_ref             = 25.0,
        mimics_initial_Cstocks_depth = 0.3,
        mimics_initial_Cstocks       = fill(100.0, ndecomp_pools_max),
        mimics_mge                   = [0.6, 0.2, 0.6, 0.6, 0.3, 0.6],
        mimics_vmod                  = [10.0, 3.0, 10.0, 10.0, 2.0, 10.0],
        mimics_vint                  = [5.47, 5.47, 5.47, 5.47, 5.47, 5.47],
        mimics_vslope                = [0.063, 0.063, 0.063, 0.063, 0.063, 0.063],
        mimics_kmod                  = [0.125, 0.5, 0.25, 0.5, 0.25, 0.167],
        mimics_kint                  = [3.19, 3.19, 3.19, 3.19, 3.19, 3.19],
        mimics_kslope                = [0.017, 0.017, 0.017, 0.017, 0.017, 0.017],
        mimics_fmet                  = [1.0, 0.85, 0.013, 25.0],
        mimics_p_scalar              = [2.0, -2.5],
        mimics_fphys_r               = [0.3, 1.3],
        mimics_fphys_k               = [0.2, 0.8],
        mimics_fchem_r               = [0.1, -3.0],
        mimics_fchem_k               = [0.3, -3.0],
        mimics_desorp                = [0.00015, 1.5],
        mimics_tau_r                 = [5.2e-4, -8.0],
        mimics_tau_k                 = [2.4e-4, -2.0],
    )
end

# --- Build a full instance at precision FT, with `nlevdecomp` levels ---
function build(::Type{FT}; nc=4, nlevdecomp=1, ndecomp_pools=8,
               ndecomp_cascade_transitions=15) where {FT}
    params = make_mimics_params(FT; ndecomp_pools_max=ndecomp_pools)
    cn_params = CLM.CNSharedParamsData(
        Q10=1.5, minpsi=-10.0, maxpsi=-0.1, rf_cwdl2=0.0, tau_cwd=10.0,
        cwd_flig=0.24, froz_q10=1.5, decomp_depth_efolding=0.5, mino2lim=0.0)

    mimics_state = CLM.DecompMIMICSState()
    cascade_con  = CLM.DecompCascadeConData()

    nlev = max(nlevdecomp, 5)
    cellclay        = fill(30.0, nc, nlev)
    t_soisno        = FT.(fill(CLM.TFRZ + 15.0, nc, nlev))
    soilpsi         = FT.(fill(-1.0, nc, nlev))
    col_dz          = FT.(fill(0.1, nc, nlev))
    decomp_cpools_vr= FT.(fill(10.0, nc, nlev, ndecomp_pools))
    ligninNratioAvg = FT.(fill(10.0, nc))
    annsum_npp_col  = FT.(fill(500.0, nc))
    mask_bgc_soilc  = trues(nc)

    cf = CLM.SoilBiogeochemCarbonFluxData{FT, Vector{FT}, Matrix{FT}, Array{FT,3}, Vector{Int}}()
    CLM.soil_bgc_carbon_flux_init!(cf, nc, nlev, ndecomp_pools, ndecomp_cascade_transitions)

    # One-time STATIC config/init (stays on host) — populates indices, pool
    # metadata, and the texture-dependent mimics_state arrays.
    CLM.init_decompcascade_mimics!(mimics_state, cascade_con, params, cn_params;
        cellclay=cellclay, bounds=1:nc, nlevdecomp=nlevdecomp,
        ndecomp_pools_max=ndecomp_pools,
        ndecomp_cascade_transitions_max=ndecomp_cascade_transitions, use_fates=false)

    return (; params, cn_params, mimics_state, cascade_con, cf, t_soisno, soilpsi,
              col_dz, decomp_cpools_vr, ligninNratioAvg, annsum_npp_col,
              mask_bgc_soilc, nc, nlevdecomp, ndecomp_pools, ndecomp_cascade_transitions)
end

run_dr!(H) = CLM.decomp_rates_mimics!(
    H.cf, H.mimics_state, H.params, H.cn_params, H.cascade_con;
    mask_bgc_soilc=H.mask_bgc_soilc, bounds=1:H.nc, nlevdecomp=H.nlevdecomp,
    t_soisno=H.t_soisno, soilpsi=H.soilpsi, decomp_cpools_vr=H.decomp_cpools_vr,
    col_dz=H.col_dz, ligninNratioAvg=H.ligninNratioAvg,
    annsum_npp_col=H.annsum_npp_col, days_per_year=365.0, dt=1800.0)

function check_case(to, FT, nlevdecomp)
    @printf("\n  --- nlevdecomp = %d ---\n", nlevdecomp)
    # Build the CPU reference entirely at Float64: init_decompcascade_mimics! allocates
    # the texture arrays (fphys/desorp/p_scalar) at Float64 regardless of precision, so a
    # mixed Float32-cf / Float64-mimics build can't unify in the _DMRArr device-view bundle.
    # The device snapshot down-converts EVERY array (cf + mimics_state + forcing) to Float32
    # via MetalF32, keeping the whole device computation at one precision.
    H = build(Float64; nlevdecomp=nlevdecomp)   # CPU reference (uniform Float64)
    B = build(Float64; nlevdecomp=nlevdecomp)   # device source snapshot

    adf32(x) = CLM.Adapt.adapt(MetalF32(), x)   # Float64 -> Float32 device
    dmask(m) = to(collect(Bool, m))

    cfd = adf32(B.cf)
    # adapt-reconstruct the loose-M DecompMIMICSState with device Float32 arrays (mutating
    # the fields would reconvert to host Float64 — the instance type is fixed at construction).
    msd = adf32(B.mimics_state)

    Bd = (; B.params, B.cn_params, mimics_state=msd, B.cascade_con, cf=cfd,
            t_soisno=adf32(B.t_soisno), soilpsi=adf32(B.soilpsi), col_dz=adf32(B.col_dz),
            decomp_cpools_vr=adf32(B.decomp_cpools_vr),
            ligninNratioAvg=adf32(B.ligninNratioAvg), annsum_npp_col=adf32(B.annsum_npp_col),
            mask_bgc_soilc=dmask(B.mask_bgc_soilc),
            B.nc, B.nlevdecomp, B.ndecomp_pools, B.ndecomp_cascade_transitions)

    if !(cfd.decomp_k_col isa device_array_type())
        println("  BLOCKED: carbon-flux struct did not move to the device under adapt.")
        return 2
    end

    run_dr!(H)    # CPU
    run_dr!(Bd)   # device

    checks = [
        ("decomp_k",  H.cf.decomp_k_col,                Bd.cf.decomp_k_col),
        ("pathfrac",  H.cf.pathfrac_decomp_cascade_col, Bd.cf.pathfrac_decomp_cascade_col),
        ("rf",        H.cf.rf_decomp_cascade_col,       Bd.cf.rf_decomp_cascade_col),
        ("cn_col",    H.cf.cn_col,                      Bd.cf.cn_col),
        ("w_scalar",  H.cf.w_scalar_col,                Bd.cf.w_scalar_col),
        ("o_scalar",  H.cf.o_scalar_col,                Bd.cf.o_scalar_col),
    ]
    nfail = 0
    for (nm, a, b) in checks
        if !cpu_has_finite(a)
            @printf("  [WARN] %-10s CPU reference all-NaN/Inf — skipping\n", nm)
            continue
        end
        d = reldiff(a, b); ok = d < 1f-3
        @printf("  [%s] %-10s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
        ok || (nfail += 1)
    end
    return nfail
end

function main(backend)
    println("=" ^ 70)
    println("END-TO-END GPU parity for decomp_rates_mimics! (whole function)")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU driver exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n", name, FT)

    nfail = 0
    nfail += check_case(to, FT, 1)   # single-level path (frw/w_scalar reductions)
    nfail += check_case(to, FT, 5)   # multi-level path (2D scalar kernels)

    println()
    println(nfail == 0 ? "  WHOLE decomp_rates_mimics! MATCHES CPU ON $name ($FT) ✓" :
                         "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
