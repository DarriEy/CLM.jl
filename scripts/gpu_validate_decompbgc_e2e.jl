# ==========================================================================
# gpu_validate_decompbgc_e2e.jl — end-to-end Metal parity for the WHOLE
# decomp_rate_constants_bgc! BGC decomposition rate-constant driver.
#
# Builds a small multi-column instance mirroring test/test_decomp_bgc.jl, runs
# decomp_rate_constants_bgc! on the CPU, converts the carbon-flux + bgc_state
# structs to Float32 + adapts to Metal, runs the SAME call on the device, and
# compares the mutated outputs (t/w/o scalars, decomp_k, pathfrac, rf) field by
# field. Several branchy config paths are exercised in separate scenarios:
#   - nlevdecomp == 1, Q10 temperature function, no anoxia
#   - nlevdecomp == 1, Q10, with anoxia (use_lch4 && anoxia O2 scalar)
#   - nlevdecomp == 5 (multi-level), Q10, with spinup_state=1 (geog-term kernel)
#   - nlevdecomp == 1, CENTURY arctangent temperature function
#
#   julia --project=scripts scripts/gpu_validate_decompbgc_e2e.jl
#
# The *_init! routines allocate Float64 regardless of the struct type param, so we
# build at Float64 and down-convert array fields to Float32 for the device snapshot
# (Metal has no Float64).
# ==========================================================================

using CLM
using Printf
import Metal   # MtlArray is the Adapt adaptor type for the device-view structs
include(joinpath(@__DIR__, "gpu_backends.jl"))

# NaN-aware mixed abs/rel diff; asserts the CPU reference is finite so a both-NaN
# false PASS can't slip through.
function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(A[i] - B[i]) / (1.0 + max(abs(A[i]), abs(B[i]))))
    end
    return m
end
cpu_has_finite(a) = any(isfinite, Array(a))

# Float32-down-converting Metal adaptor (same pattern as cstateupdate1 harness).
struct MetalF32 end
CLM.Adapt.adapt_storage(::MetalF32, x::AbstractArray{<:AbstractFloat}) = Metal.MtlArray(Float32.(x))
CLM.Adapt.adapt_storage(::MetalF32, x::AbstractArray{<:Integer})       = Metal.MtlArray(x)
CLM.Adapt.adapt_storage(::MetalF32, x::AbstractArray{Bool})            = Metal.MtlArray(x)
# DecompBGCState mixes scalar FT fields with array fields under one {FT,M<:AbstractMatrix{FT}}
# struct — down-convert scalar floats too so the whole struct lands as Float32 (FT==Float32,
# M==MtlMatrix{Float32}) and no Float64 scalar reaches a Metal kernel.
CLM.Adapt.adapt_storage(::MetalF32, x::AbstractFloat) = Float32(x)

const NC = 4
const NDECOMP_POOLS = 7
const NCASCADE = 10

# --- BGC params (mirror test/test_decomp_bgc.jl make_decomp_bgc_data) ---
function make_params()
    CLM.DecompBGCParams(
        cn_s1_bgc     = 12.0,
        cn_s2_bgc     = 12.0,
        cn_s3_bgc     = 10.0,
        rf_l1s1_bgc   = 0.39,
        rf_l2s1_bgc   = 0.55,
        rf_l3s2_bgc   = 0.29,
        rf_s2s1_bgc   = 0.55,
        rf_s2s3_bgc   = 0.55,
        rf_s3s1_bgc   = 0.55,
        rf_cwdl3_bgc  = 0.0,
        tau_l1_bgc    = 1.0 / 18.5,
        tau_l2_l3_bgc = 1.0 / 4.9,
        tau_s1_bgc    = 1.0 / 7.3,
        tau_s2_bgc    = 1.0 / 0.2,
        tau_s3_bgc    = 1.0 / 0.0045,
        cwd_fcel_bgc  = 0.45,
        bgc_initial_Cstocks       = fill(200.0, NDECOMP_POOLS),
        bgc_initial_Cstocks_depth = 0.3,
    )
end

make_cn_params() = CLM.CNSharedParamsData(
    Q10                   = 1.5,
    minpsi                = -10.0,
    maxpsi                = -0.1,
    rf_cwdl2              = 0.0,
    tau_cwd               = 10.0,
    cwd_flig              = 0.24,
    froz_q10              = 1.5,
    decomp_depth_efolding = 0.5,
    mino2lim              = 0.0,
)

# Build an initialized (cf, bgc_state) pair plus the inputs for one scenario.
function build(; nlevdecomp, use_century_tfunc=false, spinup_state=0)
    nlevfull = max(nlevdecomp, 5)
    params    = make_params()
    cn_params = make_cn_params()
    bgc_state = CLM.DecompBGCState()
    bgc_state.use_century_tfunc = use_century_tfunc
    bgc_state.normalize_q10_to_century_tfunc = !use_century_tfunc
    cascade_con = CLM.DecompCascadeConData()

    cellsand = fill(50.0, NC, nlevfull)
    CLM.init_decomp_cascade_bgc!(bgc_state, cascade_con, params, cn_params;
        cellsand=cellsand, bounds=1:NC, nlevdecomp=nlevdecomp,
        ndecomp_pools_max=NDECOMP_POOLS,
        ndecomp_cascade_transitions_max=NCASCADE, use_fates=false)

    cf = CLM.SoilBiogeochemCarbonFluxData()
    CLM.soil_bgc_carbon_flux_init!(cf, NC, nlevfull, NDECOMP_POOLS, NCASCADE)

    t_soisno = fill(CLM.TFRZ + 15.0, NC, nlevfull)
    soilpsi  = fill(-1.0, NC, nlevfull)
    col_dz   = fill(0.1, NC, nlevfull)
    o2stress = fill(0.7, NC, nlevfull)
    zsoi_vals = [0.01, 0.04, 0.09, 0.16, 0.26, 0.40, 0.58, 0.80, 1.06, 1.36]
    if nlevdecomp > length(zsoi_vals)
        append!(zsoi_vals, range(1.7, stop=3.0, length=nlevdecomp - length(zsoi_vals)))
    end

    # Spinup geog-term inputs (only used when spinup_state >= 1).
    col_gridcell = [1, 1, 2, 2]
    latdeg       = [65.0, 10.0]

    mask = trues(NC)

    return (; params, cn_params, bgc_state, cascade_con, cf,
              t_soisno, soilpsi, col_dz, o2stress, zsoi_vals,
              col_gridcell, latdeg, mask, nlevdecomp, spinup_state)
end

function run!(cf, bgc_state, params, cn_params, cascade_con, S, mask;
              use_lch4, anoxia, t_soisno, soilpsi, col_dz, o2stress, zsoi_vals,
              col_gridcell, latdeg)
    CLM.decomp_rate_constants_bgc!(cf, bgc_state, params, cn_params, cascade_con;
        mask_bgc_soilc=mask, bounds=1:NC, nlevdecomp=S.nlevdecomp,
        t_soisno=t_soisno, soilpsi=soilpsi,
        days_per_year=365.0, dt=1800.0, zsoi_vals=zsoi_vals,
        spinup_state=S.spinup_state, use_lch4=use_lch4, anoxia=anoxia,
        use_fates=false, o2stress_unsat=o2stress, col_dz=col_dz,
        col_gridcell=col_gridcell, latdeg=latdeg)
end

function scenario(name, backend; nlevdecomp, use_century_tfunc=false,
                  spinup_state=0, use_lch4=false, anoxia=false)
    nm, to, FT = backend
    H = build(; nlevdecomp, use_century_tfunc, spinup_state)
    B = build(; nlevdecomp, use_century_tfunc, spinup_state)

    # CPU reference (Float64).
    run!(H.cf, H.bgc_state, H.params, H.cn_params, H.cascade_con, H, H.mask;
         use_lch4=use_lch4, anoxia=anoxia, t_soisno=H.t_soisno, soilpsi=H.soilpsi,
         col_dz=H.col_dz, o2stress=H.o2stress, zsoi_vals=H.zsoi_vals,
         col_gridcell=H.col_gridcell, latdeg=H.latdeg)

    # Device snapshot: adapt the mutated/read state structs to Metal (Float32).
    ad(x) = CLM.Adapt.adapt(MetalF32(), x)
    cf_d  = ad(B.cf)
    bgc_d = ad(B.bgc_state)
    dmask(m) = to(collect(Bool, m))

    if !(cf_d.decomp_k_col isa Metal.MtlArray)
        @printf("  [%s] BLOCKED: cf did not move to device under adapt.\n", name)
        return 2
    end

    run!(cf_d, bgc_d, B.params, B.cn_params, B.cascade_con, B, dmask(B.mask);
         use_lch4=use_lch4, anoxia=anoxia,
         t_soisno=to(Float32.(B.t_soisno)), soilpsi=to(Float32.(B.soilpsi)),
         col_dz=to(Float32.(B.col_dz)), o2stress=to(Float32.(B.o2stress)),
         zsoi_vals=to(Float32.(B.zsoi_vals)),
         col_gridcell=to(B.col_gridcell), latdeg=to(Float32.(B.latdeg)))

    pairs = [
        ("t_scalar",  H.cf.t_scalar_col,                cf_d.t_scalar_col),
        ("w_scalar",  H.cf.w_scalar_col,                cf_d.w_scalar_col),
        ("o_scalar",  H.cf.o_scalar_col,                cf_d.o_scalar_col),
        ("decomp_k",  H.cf.decomp_k_col,                cf_d.decomp_k_col),
        ("pathfrac",  H.cf.pathfrac_decomp_cascade_col, cf_d.pathfrac_decomp_cascade_col),
        ("rf",        H.cf.rf_decomp_cascade_col,       cf_d.rf_decomp_cascade_col),
    ]
    nfail = 0
    @printf("  -- scenario: %s --\n", name)
    for (cn, a, b) in pairs
        if !cpu_has_finite(a)
            @printf("    [WARN] %-10s CPU reference all-NaN — skipping\n", cn)
            continue
        end
        d = reldiff(a, b); ok = d < 1f-3
        @printf("    [%s] %-10s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", cn, d)
        ok || (nfail += 1)
    end
    return nfail
end

function main(backend)
    println("=" ^ 70)
    println("END-TO-END Metal parity for decomp_rate_constants_bgc!")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU path exercised by the suite).")
        return 0
    end
    nm, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", nm, FT)

    nfail = 0
    nfail += scenario("nlevdecomp=1, Q10, no anoxia",    backend; nlevdecomp=1)
    nfail += scenario("nlevdecomp=1, Q10, anoxia",       backend; nlevdecomp=1, use_lch4=true, anoxia=true)
    nfail += scenario("nlevdecomp=5, Q10, spinup=1",     backend; nlevdecomp=5, spinup_state=1)
    nfail += scenario("nlevdecomp=1, CENTURY tfunc",     backend; nlevdecomp=1, use_century_tfunc=true)

    println()
    println(nfail == 0 ? "  WHOLE decomp_rate_constants_bgc! MATCHES CPU ON $nm ($FT) ✓" :
                         "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
