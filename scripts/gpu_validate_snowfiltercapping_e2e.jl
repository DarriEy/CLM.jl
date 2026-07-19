# ==========================================================================
# gpu_validate_snowfiltercapping_e2e.jl — end-to-end GPU parity for the A7
# snow-filter + snow-capping group, run WHOLE-FUNCTION on the device:
#   build_snow_filter!                       (snow / no-snow column split by snl)
#   init_flux_snow_capping!                  (zero the 4 capping flux fields)
#   bulk_flux_snow_capping_fluxes!           (excess + reset override + flux split)
#     → snow_capping_excess!                 (standard cap + reset/glc override)
#   update_state_remove_snow_capping_fluxes! (subtract capped mass from bottom lyr)
#   snow_capping_update_dz_and_aerosols!     (rescale dz + 8 aerosol species)
#
# Columns span the snl branches AND the capping branches so every path runs on
# BOTH backends:
#   col 1: snl=-1, h2osno  > H2OSNO_MAX           -> snow filter + capping fires
#   col 2: snl=-3, h2osno  < H2OSNO_MAX           -> snow filter, NO capping
#   col 3: snl=-3, h2osno >> H2OSNO_MAX           -> snow filter + capping (deep)
#   col 4: snl= 0 (no resolved snow)              -> NO-SNOW filter branch
#
# CPU reference fields are asserted finite so a both-NaN false PASS cannot slip
# through. snl-branch coverage + per-column capping flags are printed explicitly.
#
#   julia --project=scripts scripts/gpu_validate_snowfiltercapping_e2e.jl
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

const NC = 4
const NLEVSNO = 5
const JJ_BOTTOM = NLEVSNO              # Fortran layer 0 -> Julia nlevsno (bottom snow layer)

# Build a self-contained snow-capping instance at precision FT.
function build(::Type{FT}) where {FT}
    CLM.varpar.nlevsno = NLEVSNO   # aerosol_init! reads varpar.nlevsno
    nc = NC
    H2OSNO_MAX = FT(CLM.H2OSNO_MAX)

    snl        = Int[-1, -3, -3, 0]
    mask_nolake = trues(nc)

    # h2osno per column (mm SWE): cols 1 & 3 exceed the cap.
    h2osno_total = FT[H2OSNO_MAX + FT(500.0),   # col1: small excess
                      H2OSNO_MAX - FT(2000.0),  # col2: under cap
                      H2OSNO_MAX + FT(8000.0),   # col3: large excess
                      FT(0.0)]                   # col4: no snow

    topo = fill(FT(1000.0), nc)
    col_landunit = Int[1, 2, 3, 4]
    lun_itype    = Int[CLM.ISTSOIL, CLM.ISTSOIL, CLM.ISTSOIL, CLM.ISTSOIL]

    # bottom-layer state (one value per column at jj_bottom)
    dz_bottom  = FT[FT(0.10), FT(0.20), FT(0.25), FT(0.0)]
    ice_bottom = FT[FT(80.0), FT(120.0), FT(200.0), FT(0.0)]
    liq_bottom = FT[FT(5.0),  FT(8.0),  FT(15.0),  FT(0.0)]

    # aerosol masses in the bottom snow layer (nc x nlevsno)
    aer = CLM.AerosolData{FT}()
    CLM.aerosol_init!(aer, nc)
    for fld in (:mss_bcphi_col, :mss_bcpho_col, :mss_ocphi_col, :mss_ocpho_col,
                :mss_dst1_col, :mss_dst2_col, :mss_dst3_col, :mss_dst4_col)
        M = getfield(aer, fld)
        for c in 1:nc, j in 1:NLEVSNO
            M[c, j] = FT(0.0)
        end
        # put a nonzero mass in the bottom layer so frac_adjust scaling is visible
        for c in 1:nc
            M[c, JJ_BOTTOM] = FT(c) * FT(1.0e-6)
        end
    end

    # capping flux outputs + scratch
    qflx_snwcp_ice           = fill(FT(NaN), nc)
    qflx_snwcp_liq           = fill(FT(NaN), nc)
    qflx_snwcp_discarded_ice = fill(FT(NaN), nc)
    qflx_snwcp_discarded_liq = fill(FT(NaN), nc)
    rho_orig_bottom = zeros(FT, nc)
    frac_adjust     = zeros(FT, nc)

    # snow-filter outputs
    mask_snow   = falses(nc)
    mask_nosnow = falses(nc)
    mask_capping = falses(nc)

    return (; nc, snl, mask_nolake, h2osno_total, topo, col_landunit, lun_itype,
            dz_bottom, ice_bottom, liq_bottom, aer,
            qflx_snwcp_ice, qflx_snwcp_liq, qflx_snwcp_discarded_ice,
            qflx_snwcp_discarded_liq, rho_orig_bottom, frac_adjust,
            mask_snow, mask_nosnow, mask_capping, dtime = FT(1800.0))
end

# Run the whole filter+capping pipeline on the supplied (host or device) arrays.
function run_pipeline!(S, nstep::Int)
    bounds = 1:S.nc

    CLM.build_snow_filter!(S.mask_snow, S.mask_nosnow, S.snl, S.mask_nolake, bounds)

    CLM.init_flux_snow_capping!(
        S.qflx_snwcp_ice, S.qflx_snwcp_liq,
        S.qflx_snwcp_discarded_ice, S.qflx_snwcp_discarded_liq,
        S.mask_snow, bounds)

    CLM.bulk_flux_snow_capping_fluxes!(
        S.mask_capping, S.rho_orig_bottom, S.frac_adjust,
        S.qflx_snwcp_ice, S.qflx_snwcp_liq,
        S.qflx_snwcp_discarded_ice, S.qflx_snwcp_discarded_liq,
        S.dtime, S.dz_bottom, S.topo, S.h2osno_total,
        S.ice_bottom, S.liq_bottom,
        S.col_landunit, S.lun_itype,
        S.mask_snow, bounds, NLEVSNO, nstep)

    CLM.update_state_remove_snow_capping_fluxes!(
        S.ice_bottom, S.liq_bottom, S.dtime,
        S.qflx_snwcp_ice, S.qflx_snwcp_liq,
        S.qflx_snwcp_discarded_ice, S.qflx_snwcp_discarded_liq,
        S.mask_capping, bounds)

    CLM.snow_capping_update_dz_and_aerosols!(
        S.dz_bottom, S.aer, JJ_BOTTOM,
        S.rho_orig_bottom, S.ice_bottom, S.frac_adjust,
        S.mask_capping, bounds)
    return nothing
end

function main(backend)
    println("=" ^ 72)
    println("END-TO-END GPU parity for build_snow_filter! + snow-capping group (A7)")
    println("=" ^ 72)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU path exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    # Snow-resetting OFF so the standard cap branch drives the result.
    CLM.RESET_SNOW[] = false
    CLM.RESET_SNOW_GLC[] = false

    H = build(FT)     # CPU reference
    B = build(FT)     # device snapshot source

    # device-resident copies (Bool masks, not BitArrays)
    db(x) = to(collect(Bool, x))
    Sd = (; nc = B.nc,
            snl = to(B.snl), mask_nolake = db(B.mask_nolake),
            h2osno_total = to(B.h2osno_total), topo = to(B.topo),
            col_landunit = to(B.col_landunit), lun_itype = to(B.lun_itype),
            dz_bottom = to(B.dz_bottom), ice_bottom = to(B.ice_bottom),
            liq_bottom = to(B.liq_bottom), aer = CLM.Adapt.adapt(device_array_type(), B.aer),
            qflx_snwcp_ice = to(B.qflx_snwcp_ice), qflx_snwcp_liq = to(B.qflx_snwcp_liq),
            qflx_snwcp_discarded_ice = to(B.qflx_snwcp_discarded_ice),
            qflx_snwcp_discarded_liq = to(B.qflx_snwcp_discarded_liq),
            rho_orig_bottom = to(B.rho_orig_bottom), frac_adjust = to(B.frac_adjust),
            mask_snow = db(B.mask_snow), mask_nosnow = db(B.mask_nosnow),
            mask_capping = db(B.mask_capping), dtime = B.dtime)

    if !(Sd.dz_bottom isa device_array_type())
        println("  BLOCKED: device arrays did not move under adapt/to.")
        return 2
    end

    nstep = 1
    run_pipeline!(H,  nstep)   # CPU
    run_pipeline!(Sd, nstep)   # device

    aer_fields = (:mss_bcphi_col, :mss_bcpho_col, :mss_ocphi_col, :mss_ocpho_col,
                  :mss_dst1_col, :mss_dst2_col, :mss_dst3_col, :mss_dst4_col)

    checks = Any[
        ("mask_snow",         H.mask_snow,                Sd.mask_snow),
        ("mask_nosnow",       H.mask_nosnow,              Sd.mask_nosnow),
        ("mask_capping",      H.mask_capping,             Sd.mask_capping),
        ("qflx_snwcp_ice",    H.qflx_snwcp_ice,           Sd.qflx_snwcp_ice),
        ("qflx_snwcp_liq",    H.qflx_snwcp_liq,           Sd.qflx_snwcp_liq),
        ("qflx_disc_ice",     H.qflx_snwcp_discarded_ice, Sd.qflx_snwcp_discarded_ice),
        ("qflx_disc_liq",     H.qflx_snwcp_discarded_liq, Sd.qflx_snwcp_discarded_liq),
        ("rho_orig_bottom",   H.rho_orig_bottom,          Sd.rho_orig_bottom),
        ("frac_adjust",       H.frac_adjust,              Sd.frac_adjust),
        ("ice_bottom",        H.ice_bottom,               Sd.ice_bottom),
        ("liq_bottom",        H.liq_bottom,               Sd.liq_bottom),
        ("dz_bottom",         H.dz_bottom,                Sd.dz_bottom),
    ]
    for fld in aer_fields
        push!(checks, (String(fld), getfield(H.aer, fld), getfield(Sd.aer, fld)))
    end

    # Convert Bool masks to 0/1 for the numeric diff
    tonum(x) = Array(x) isa AbstractArray{Bool} ? Float64.(Array(x)) : Array(x)

    nfail = 0
    for (nm, a, b) in checks
        an = tonum(a); bn = tonum(b)
        if !any(isfinite, an)
            @printf("  [WARN] %-18s CPU reference all-NaN/Inf — skipping (no signal)\n", nm)
            continue
        end
        d = reldiff(an, bn); ok = d < 1f-3
        @printf("  [%s] %-18s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
        ok || (nfail += 1)
    end

    println()
    sc = Array(H.mask_snow); ns = Array(H.mask_nosnow); cp = Array(H.mask_capping)
    for c in 1:H.nc
        @printf("    col %d: snl=%2d  snow=%-5s nosnow=%-5s capping=%-5s\n",
                c, H.snl[c], sc[c], ns[c], cp[c])
    end

    # --- second pass with snow-RESETTING ON so the reset-override kernel in
    # snow_capping_excess! also runs on both backends (nstep within the reset
    # window; RESET_SNOW_H2OSNO=35mm << H2OSNO_MAX so every snow column resets). ---
    CLM.RESET_SNOW[] = true
    CLM.RESET_SNOW_GLC[] = false
    Hr = build(FT); Br = build(FT)
    Sdr = (; nc = Br.nc,
            snl = to(Br.snl), mask_nolake = db(Br.mask_nolake),
            h2osno_total = to(Br.h2osno_total), topo = to(Br.topo),
            col_landunit = to(Br.col_landunit), lun_itype = to(Br.lun_itype),
            dz_bottom = to(Br.dz_bottom), ice_bottom = to(Br.ice_bottom),
            liq_bottom = to(Br.liq_bottom), aer = CLM.Adapt.adapt(device_array_type(), Br.aer),
            qflx_snwcp_ice = to(Br.qflx_snwcp_ice), qflx_snwcp_liq = to(Br.qflx_snwcp_liq),
            qflx_snwcp_discarded_ice = to(Br.qflx_snwcp_discarded_ice),
            qflx_snwcp_discarded_liq = to(Br.qflx_snwcp_discarded_liq),
            rho_orig_bottom = to(Br.rho_orig_bottom), frac_adjust = to(Br.frac_adjust),
            mask_snow = db(Br.mask_snow), mask_nosnow = db(Br.mask_nosnow),
            mask_capping = db(Br.mask_capping), dtime = Br.dtime)
    run_pipeline!(Hr,  1)
    run_pipeline!(Sdr, 1)
    CLM.RESET_SNOW[] = false
    dr = max(reldiff(Hr.qflx_snwcp_discarded_ice, Sdr.qflx_snwcp_discarded_ice),
             reldiff(Hr.frac_adjust, Sdr.frac_adjust),
             reldiff(Hr.dz_bottom, Sdr.dz_bottom))
    okr = dr < 1f-3 && any(Array(Hr.mask_capping))
    @printf("\n  [%s] RESET-override pass     rel|dev-cpu| = %.3e (capping cols=%d, apply_runoff=false)\n",
            okr ? "PASS" : "FAIL", dr, count(Array(Hr.mask_capping)))
    okr || (nfail += 1)

    println()
    println(nfail == 0 ? "  WHOLE filter+capping group MATCHES CPU ON $name ($FT)" :
                         "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
