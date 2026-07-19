# ==========================================================================
# gpu_validate_ch4aere_e2e.jl — end-to-end GPU parity for the WHOLE ch4_aere!
# BGC driver (aerenchyma transport: site_ox_aere! math inlined into a per-PATCH
# kernel that scatters into ch4_aere_depth / ch4_tran_depth / o2_aere_depth).
#
# Builds a small multi-column / multi-patch CH4Data mirroring test/test_methane.jl
# (there is no ch4_aere! unit test, so the call is assembled from the ch4! call
# site), runs ch4_aere! on the CPU, converts the CH4Data + arg arrays to Float32
# and adapts to the GPU, runs the SAME call on the device, and compares the mutated
# column outputs field-by-field. The scenario exercises both the unsaturated
# (sat=0) and saturated (sat=1) branches, the transpirationloss + prognostic
# aerenchyma-oxidation flags, and a FATES column (forced-vegetated path).
#
#   julia --project=scripts scripts/gpu_validate_ch4aere_e2e.jl
#
# CH4Data is {FT,…}-parametric + @adapt_structure'd, so it adapts whole. ch4_aere!
# also takes loose arg arrays (watsat, rootfr, …) + Int/Bool index/flag vectors
# which we move to the device preserving their eltype (floats→Float32, Int/Bool
# as-is). Metal has no Float64, so the device snapshot is Float32 throughout.
# ==========================================================================

using CLM
using Printf
include(joinpath(@__DIR__, "gpu_backends.jl"))

function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(A[i] - B[i]) / (1.0 + max(abs(A[i]), abs(B[i]))))
    end
    return m
end
cpu_has_finite(a) = any(isfinite, Array(a))

# Float32-down-converting Metal adaptor (CH4Data is concretely typed; build at
# Float64 then down-convert float arrays as adapt rebuilds the struct).
struct MetalF32 end
CLM.Adapt.adapt_storage(::MetalF32, x::AbstractArray{<:AbstractFloat}) = device_array_type()(Float32.(x))
CLM.Adapt.adapt_storage(::MetalF32, x::AbstractArray{<:Integer})       = device_array_type()(x)
CLM.Adapt.adapt_storage(::MetalF32, x::AbstractArray{Bool})            = device_array_type()(x)

const NC = 4
const NP = 6
const NG = 2
const NLEVSOI = 5
const NOVEG = 0

# Build the CPU reference state (Float64) + the loose ch4_aere! args.
function build()
    params = CLM.CH4Params()
    ch4vc  = CLM.CH4VarCon()
    ch4vc.transpirationloss = true
    ch4vc.use_aereoxid_prog = true   # exercise the O2-diffusion branch

    z  = zeros(NC, NLEVSOI); dz = fill(0.1, NC, NLEVSOI); zi = zeros(NC, NLEVSOI)
    for j in 1:NLEVSOI, c in 1:NC
        z[c, j]  = (j - 0.5) * 0.1
        zi[c, j] = j * 0.1
    end

    ch4 = CLM.CH4Data{Float64}(
        # scatter targets + the prod input ch4_aere! reads
        ch4_aere_depth_sat_col   = zeros(NC, NLEVSOI),
        ch4_aere_depth_unsat_col = zeros(NC, NLEVSOI),
        ch4_tran_depth_sat_col   = zeros(NC, NLEVSOI),
        ch4_tran_depth_unsat_col = zeros(NC, NLEVSOI),
        o2_aere_depth_sat_col    = zeros(NC, NLEVSOI),
        o2_aere_depth_unsat_col  = zeros(NC, NLEVSOI),
        ch4_prod_depth_sat_col   = fill(1.0e-8, NC, NLEVSOI),
        ch4_prod_depth_unsat_col = fill(2.0e-8, NC, NLEVSOI),
        conc_ch4_sat_col         = fill(1.0e-4, NC, NLEVSOI),
        conc_ch4_unsat_col       = fill(1.0e-5, NC, NLEVSOI),
        conc_o2_sat_col          = fill(0.01,  NC, NLEVSOI),
        conc_o2_unsat_col        = fill(0.05,  NC, NLEVSOI),
        # c_atm: CH4 (col 1) low, O2 (col 2) high so the prognostic O2-aerenchyma
        # uptake branch produces nonzero o2_aere_depth (else it smooth_max's to 0).
        c_atm_grc                = [0.03 8.0 0.0; 0.03 8.0 0.0],
        annavg_agnpp_patch       = fill(1.0e-5, NP),
        annavg_bgnpp_patch       = fill(1.5e-5, NP),
        grnd_ch4_cond_patch      = fill(0.01,  NP),
    )

    # Soil / patch driver arrays
    watsat     = fill(0.45, NC, NLEVSOI)
    h2osoi_vol = fill(0.30, NC, NLEVSOI)
    t_soisno   = fill(CLM.TFRZ + 15.0, NC, NLEVSOI)
    rootr      = fill(0.1,  NP, NLEVSOI)
    rootfr     = fill(0.2,  NP, NLEVSOI)
    elai          = fill(2.0,  NP)
    qflx_tran_veg = fill(1.0e-3, NP)
    annsum_npp    = fill(1.0e-4, NP)

    # Topology + branch selectors. patch 1 is bare ground (itype==NOVEG → not
    # vegetated); patch 5/6 sit on column 4 which is FATES (forced vegetated).
    patch_column = [1, 2, 2, 3, 4, 4]
    col_gridcell = [1, 1, 2, 2]
    patch_itype  = [NOVEG, 1, 2, 3, 4, 5]
    patch_wtcol  = [1.0, 0.6, 0.4, 1.0, 0.5, 0.5]
    is_fates     = [false, false, false, true]   # column 4 is FATES
    jwt          = [2, 1, 0, 3]                   # vary water-table layer
    mask_soil    = trues(NC)
    mask_soilp   = trues(NP)

    return (; params, ch4vc, ch4, mask_soil, mask_soilp, patch_column, patch_itype,
              patch_wtcol, col_gridcell, is_fates, watsat, h2osoi_vol, t_soisno,
              rootr, rootfr, elai, qflx_tran_veg, annsum_npp, z, dz, jwt)
end

# Run ch4_aere! on host (Float64) or device (Float32) — the masks/index/flag
# vectors and loose arrays are pre-converted by the caller.
run_aere!(S, ch4, m_soil, m_soilp, pcol, pitype, pwt, cgrid, fates,
          watsat, h2osoi_vol, t_soisno, rootr, rootfr, elai, qtran, anpp,
          z, dz, jwt, sat, dtime) =
    CLM.ch4_aere!(ch4, S.params, S.ch4vc, m_soil, m_soilp, pcol, pitype, pwt,
                  cgrid, fates, watsat, h2osoi_vol, t_soisno, rootr, rootfr,
                  elai, qtran, anpp, z, dz, jwt, sat, false, NLEVSOI, dtime, NOVEG)

function main(backend)
    println("=" ^ 70)
    println("END-TO-END GPU parity for ch4_aere! (BGC aerenchyma scatter)")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU path exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    H = build()    # CPU reference (Float64)
    B = build()    # source for the device snapshot

    ad(x) = CLM.Adapt.adapt(MetalF32(), x)
    ch4_d = ad(B.ch4)
    if !(ch4_d.ch4_aere_depth_unsat_col isa device_array_type())
        println("  BLOCKED: CH4Data did not move to the device under adapt.")
        return 2
    end

    # Device copies of the loose args (floats→Float32, Int/Bool eltype preserved).
    dF(x)  = to(Float32.(x))
    dI(x)  = to(x)
    dB(x)  = to(collect(Bool, x))
    d_watsat = dF(B.watsat); d_h2v = dF(B.h2osoi_vol); d_tsoi = dF(B.t_soisno)
    d_rootr = dF(B.rootr); d_rootfr = dF(B.rootfr); d_elai = dF(B.elai)
    d_qtran = dF(B.qflx_tran_veg); d_anpp = dF(B.annsum_npp)
    d_z = dF(B.z); d_dz = dF(B.dz)
    d_pcol = dI(B.patch_column); d_pit = dI(B.patch_itype); d_pwt = dF(B.patch_wtcol)
    d_cgrid = dI(B.col_gridcell); d_jwt = dI(B.jwt)
    d_fates = dB(B.is_fates); d_msoil = dB(B.mask_soil); d_msoilp = dB(B.mask_soilp)

    dt_h = 1800.0
    dt_d = FT(1800.0)

    nfail = 0
    for sat in (0, 1)
        tag = sat == 0 ? "unsat" : "sat"

        # CPU
        run_aere!(H, H.ch4, H.mask_soil, H.mask_soilp, H.patch_column, H.patch_itype,
                  H.patch_wtcol, H.col_gridcell, H.is_fates, H.watsat, H.h2osoi_vol,
                  H.t_soisno, H.rootr, H.rootfr, H.elai, H.qflx_tran_veg, H.annsum_npp,
                  H.z, H.dz, H.jwt, sat, dt_h)
        # Device
        run_aere!(B, ch4_d, d_msoil, d_msoilp, d_pcol, d_pit, d_pwt, d_cgrid, d_fates,
                  d_watsat, d_h2v, d_tsoi, d_rootr, d_rootfr, d_elai, d_qtran, d_anpp,
                  d_z, d_dz, d_jwt, sat, dt_d)

        if sat == 0
            checks = [
                ("ch4_aere_depth_unsat", H.ch4.ch4_aere_depth_unsat_col, ch4_d.ch4_aere_depth_unsat_col),
                ("ch4_tran_depth_unsat", H.ch4.ch4_tran_depth_unsat_col, ch4_d.ch4_tran_depth_unsat_col),
                ("o2_aere_depth_unsat",  H.ch4.o2_aere_depth_unsat_col,  ch4_d.o2_aere_depth_unsat_col),
            ]
        else
            checks = [
                ("ch4_aere_depth_sat", H.ch4.ch4_aere_depth_sat_col, ch4_d.ch4_aere_depth_sat_col),
                ("ch4_tran_depth_sat", H.ch4.ch4_tran_depth_sat_col, ch4_d.ch4_tran_depth_sat_col),
                ("o2_aere_depth_sat",  H.ch4.o2_aere_depth_sat_col,  ch4_d.o2_aere_depth_sat_col),
            ]
        end
        for (nm, a, b) in checks
            if !cpu_has_finite(a)
                @printf("  [WARN] %-24s CPU reference is all-NaN/Inf — skipping\n", nm)
                continue
            end
            d = reldiff(a, b); ok = d < 1f-3
            @printf("  [%s] %-24s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
            ok || (nfail += 1)
        end
    end
    println()
    println(nfail == 0 ? "  WHOLE ch4_aere! MATCHES CPU ON $name ($FT) ✓" :
                         "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
