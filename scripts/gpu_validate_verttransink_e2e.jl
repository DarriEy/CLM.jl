# ==========================================================================
# gpu_validate_verttransink_e2e.jl — end-to-end Metal parity for the WHOLE
# compute_effec_rootfrac_and_vert_tran_sink! wrapper (and its children:
# _Default, _HydStress_Roads, _HydStress).
#
# Builds a small Float32 instance with three column groups (pervious road,
# non-istsoil "other", natural vegetation), each with multiple patches, runs
# the wrapper on the CPU, adapts every state struct (+ the hydrology mask) to
# Metal, runs the SAME wrapper on the device, and compares the mutated outputs
# field-by-field. Both the default path (use_hydrstress=false) and the
# hydraulic-stress path (use_hydrstress=true) are exercised, covering:
#   * the per-column transpiration-weighted accumulation/normalize/sink kernel
#   * the per-column hydraulic-stress kernel (per-patch hydr-redist + phs_neg)
#   * the host filter construction + device-filter movement in the wrapper.
#
#   julia --project=scripts scripts/gpu_validate_verttransink_e2e.jl
# ==========================================================================

using CLM
using Printf
import Metal   # Metal-specific; MtlArray is the Adapt adaptor type for the structs
include(joinpath(@__DIR__, "gpu_backends.jl"))

# reldiff: NaN-aware (both-NaN agrees; one-sided NaN flags divergence).
function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(A[i] - B[i]) / (1.0 + max(abs(A[i]), abs(B[i]))))
    end
    return m
end
cpu_has_finite(a) = any(isfinite, Array(a))

# Builds one Float-precision instance. Three columns:
#   c=1 pervious road (urban landunit)     -> _HydStress_Roads / _Default
#   c=2 crop (not road, not istsoil)        -> _HydStress / _Default
#   c=3 natural vegetation (istsoil)        -> _HydStress / _Default
# Columns 1,3 carry 2 patches; column 2 carries 1 patch.
function build(::Type{FT}) where {FT}
    nlevsno = 5; nlevgrnd = 10; nlevsoi = 4; nlevmaxurbgrnd = 10
    CLM.varpar.nlevsno        = nlevsno
    CLM.varpar.nlevgrnd       = nlevgrnd
    CLM.varpar.nlevurb        = 5
    CLM.varpar.nlevmaxurbgrnd = nlevmaxurbgrnd
    CLM.varpar.nlevsoi        = nlevsoi
    CLM.varpar.nlevlak        = 10

    nc = 3; nl = 3
    patches_per_col = [2, 1, 2]
    np = sum(patches_per_col)

    # --- Column ---
    col = CLM.ColumnData{FT}(); CLM.column_init!(col, nc)
    col.itype[1] = CLM.ICOL_ROAD_PERV
    col.itype[2] = CLM.ISTCROP            # arbitrary non-road col itype
    col.itype[3] = 1
    col.landunit[1] = 1; col.landunit[2] = 2; col.landunit[3] = 3
    col.active[1:nc] .= true
    pidx = 1
    for c in 1:nc
        col.patchi[c]   = pidx
        col.npatches[c] = patches_per_col[c]
        col.patchf[c]   = pidx + patches_per_col[c] - 1
        pidx += patches_per_col[c]
    end
    for c in 1:nc, j in 1:nlevsoi
        col.z[c, j + nlevsno] = 0.1 * j
    end

    # --- Landunit ---
    lun = CLM.LandunitData{FT}(); CLM.landunit_init!(lun, nl)
    lun.itype[1] = CLM.ISTURB_MD          # pervious-road landunit
    lun.itype[2] = CLM.ISTCROP            # NOT istsoil -> "other" group
    lun.itype[3] = CLM.ISTSOIL            # natural vegetation group

    # --- Patch ---
    patch_data = CLM.PatchData{FT}(); CLM.patch_init!(patch_data, np)
    patch_data.active[1:np] .= true
    # column membership (defaults to ISPVAL otherwise)
    for c in 1:nc, p in col.patchi[c]:col.patchf[c]
        patch_data.column[p] = c
    end
    # mix of weights; include a >0 weight per column so hydstress contributes
    patch_data.wtcol[1] = 0.6; patch_data.wtcol[2] = 0.4   # col 1
    patch_data.wtcol[3] = 1.0                              # col 2
    patch_data.wtcol[4] = 0.7; patch_data.wtcol[5] = 0.3   # col 3

    # --- Soil state ---
    ss = CLM.SoilStateData{FT}(); CLM.soilstate_init!(ss, np, nc)
    for p in 1:np, j in 1:nlevsoi
        ss.rootr_patch[p, j] = 0.1 * j + 0.05 * p   # distinct, finite
    end
    for c in 1:nc, j in 1:nlevgrnd
        ss.rootr_col[c, j] = 0.0
        ss.smp_l_col[c, j] = -300.0 - 50.0 * c      # finite matric potential
    end
    for p in 1:np, j in 1:nlevsoi
        ss.k_soil_root_patch[p, j] = 0.005 + 0.001 * j
    end

    # --- Water flux bulk ---
    wfb = CLM.WaterFluxBulkData{FT}(); CLM.waterfluxbulk_init!(wfb, nc, np, nl, 1)
    tran_patch = [0.001, 0.002, 0.0015, 0.0025, 0.0012]
    for p in 1:np
        wfb.wf.qflx_tran_veg_patch[p] = tran_patch[p]
    end
    for c in 1:nc
        wfb.wf.qflx_tran_veg_col[c] = 0.001 * c
    end
    for c in 1:nc, j in 1:nlevsoi
        wfb.qflx_rootsoi_col[c, j] = 0.0
    end
    for c in 1:nc
        wfb.qflx_phs_neg_col[c] = 0.0
    end
    for p in 1:np
        wfb.qflx_hydr_redist_patch[p] = 0.0
    end

    # --- Canopy state ---
    cs = CLM.CanopyStateData{FT}(); CLM.canopystate_init!(cs, np)
    cs.frac_veg_nosno_patch[1:np] .= 1
    for p in 1:np
        cs.vegwp_patch[p, 4] = -700.0 - 30.0 * p     # root segment (root_idx=4)
    end

    # --- Energy flux (unused by the math, required by signature) ---
    ef = CLM.EnergyFluxData{FT}(); CLM.energyflux_init!(ef, np, nc, nl, 1)

    S = (; col, lun, patch_data, ss, wfb, cs, ef)
    return (; nc, np, nl, nlevsoi, S)
end

run_wrapper!(S, nlevsoi, mask_hydrology, nc; use_hydrstress) =
    CLM.compute_effec_rootfrac_and_vert_tran_sink!(
        1:nc, nlevsoi, mask_hydrology,
        S.ss, S.cs, S.wfb, S.ef, S.col, S.lun, S.patch_data;
        use_hydrstress = use_hydrstress)

function run_one(to, FT, use_hydrstress::Bool)
    H = build(FT)               # CPU reference
    B = build(FT)               # source for the device snapshot

    ad(x) = CLM.Adapt.adapt(Metal.MtlArray, x)
    Sd = map(ad, B.S)

    if !(Sd.ss.rootr_col isa Metal.MtlArray)
        println("  BLOCKED: a state struct did not move to the device under adapt.")
        return (2, [])
    end

    # mask_hydrology stays a BitVector on the host (the wrapper builds host
    # filters from it + col/lun metadata, then moves the int filters to device).
    mask_h_cpu = trues(H.nc)
    mask_h_dev = trues(H.nc)   # device path reads the same host BitVector

    run_wrapper!(H.S, H.nlevsoi, mask_h_cpu, H.nc; use_hydrstress = use_hydrstress)
    run_wrapper!(Sd,   B.nlevsoi, mask_h_dev, B.nc; use_hydrstress = use_hydrstress)

    checks = [
        ("rootr_col",            H.S.ss.rootr_col,                Sd.ss.rootr_col),
        ("qflx_rootsoi_col",     H.S.wfb.qflx_rootsoi_col,        Sd.wfb.qflx_rootsoi_col),
        ("qflx_phs_neg_col",     H.S.wfb.qflx_phs_neg_col,        Sd.wfb.qflx_phs_neg_col),
        ("qflx_hydr_redist_pat", H.S.wfb.qflx_hydr_redist_patch,  Sd.wfb.qflx_hydr_redist_patch),
    ]
    return (0, checks)
end

function main(backend)
    println("=" ^ 70)
    println("END-TO-END Metal parity for compute_effec_rootfrac_and_vert_tran_sink!")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU wrapper exercised by suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    nfail = 0
    for use_hs in (false, true)
        @printf("  --- use_hydrstress = %s ---\n", use_hs)
        status, checks = run_one(to, FT, use_hs)
        if status != 0
            return status
        end
        for (nm, a, b) in checks
            if !cpu_has_finite(a)
                @printf("  [WARN] %-22s CPU reference all-NaN/Inf — skipping (no signal)\n", nm)
                continue
            end
            d = reldiff(a, b); ok = d < 1f-3
            @printf("  [%s] %-22s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
            ok || (nfail += 1)
        end
        println()
    end

    println(nfail == 0 ?
        "  WHOLE compute_effec_rootfrac_and_vert_tran_sink! MATCHES CPU ON $name ($FT) OK" :
        "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
