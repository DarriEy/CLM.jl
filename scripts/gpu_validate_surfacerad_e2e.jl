# ==========================================================================
# gpu_validate_surfacerad_e2e.jl — end-to-end GPU parity for the WHOLE
# radiation pair: canopy_sun_shade_fracs! + surface_radiation!.
#
# Builds a small Float32 instance with TWO patches that exercise both sides of
# the major branches:
#   patch 1: non-urban SOIL column WITH snow layers (snl=-2) — drives the
#            snow-layer absorption pass, the SNICAR no-aerosol sums, the local-
#            noon diagnostics and snow diagnostics (snl<0).
#   patch 2: URBAN column (snl=0, no snow) — drives the urban pass and the
#            no-snow / fsun-zeroing paths.
# Runs BOTH whole functions on the CPU, adapts every state struct (+ masks /
# forcing) to the GPU, runs the SAME calls on the device, and compares every
# mutated field with reldiff.
#
# CRITICAL: every required param/input is set to a real CLM default so the CPU
# reference is FINITE — reldiff silently PASSES when both sides are NaN, so the
# harness ASSERTS the CPU reference fields are finite before trusting parity.
#
#   julia --project=scripts scripts/gpu_validate_surfacerad_e2e.jl
# ==========================================================================

using CLM
using Printf
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

# use_snicar_frc + use_SSRE are exercised so the SNICAR / snow-free passes run.
const USE_SNICAR = true
const USE_SSRE   = true

function build(::Type{FT}) where {FT}
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
    np = 2; nc = 2; ng = 1; nl = 2
    numrad  = CLM.NUMRAD
    nlevsno = CLM.varpar.nlevsno
    nlev_sno1 = nlevsno + 1

    surfalb = CLM.SurfaceAlbedoData{FT}(); CLM.surfalb_init!(surfalb, np, nc, ng)
    CLM.surfalb_init_cold!(surfalb, 1:nc, 1:np)
    canopystate = CLM.CanopyStateData{FT}(); CLM.canopystate_init!(canopystate, np)
    solarabs = CLM.SolarAbsorbedData{FT}(); CLM.solarabs_init!(solarabs, np, nl)
    surfrad  = CLM.SurfaceRadiationData{FT}(); CLM.surfrad_init!(surfrad, np)
    waterdiag = CLM.WaterDiagnosticBulkData{FT}(); CLM.waterdiagnosticbulk_init!(waterdiag, nc, np, nl, ng)

    col = CLM.ColumnData{FT}(); CLM.column_init!(col, nc)
    lun = CLM.LandunitData{FT}(); CLM.landunit_init!(lun, nl)
    grc = CLM.GridcellData{FT}(); CLM.gridcell_init!(grc, ng)
    pch = CLM.PatchData{FT}(); CLM.patch_init!(pch, np)

    # --- topology: patch1 -> col1/lun1 (soil, snow), patch2 -> col2/lun2 (urban) ---
    pch.column   .= [1, 2]
    pch.gridcell .= [1, 1]
    pch.landunit .= [1, 2]
    col.snl      .= [-2, 0]
    col.landunit .= [1, 2]
    lun.itype    .= [CLM.ISTSOIL, CLM.ISTURB_TBD]
    grc.londeg   .= [0.0]          # local noon at current_tod=43200 → exercises ln branch

    # --- canopy: patch1 has a 1-layer canopy; patch2 (urban) bare ---
    surfalb.nrad_patch .= [1, 1]
    surfalb.tlai_z_patch[1, 1] = 2.0; surfalb.tlai_z_patch[2, 1] = 0.0
    surfalb.fsun_z_patch[1, 1] = 0.6; surfalb.fsun_z_patch[2, 1] = 0.0
    surfalb.fabd_sun_z_patch[1, 1] = 0.30; surfalb.fabd_sha_z_patch[1, 1] = 0.10
    surfalb.fabi_sun_z_patch[1, 1] = 0.25; surfalb.fabi_sha_z_patch[1, 1] = 0.08
    canopystate.elai_patch .= [2.0, 0.0]
    canopystate.esai_patch .= [0.5, 0.0]
    canopystate.tlai_patch .= [2.0, 0.0]
    canopystate.fsun_patch .= [0.5, 0.5]

    # --- albedos / canopy transfer (both columns, both bands) ---
    for c in 1:nc, ib in 1:numrad
        surfalb.albgrd_col[c, ib]     = 0.20 + 0.05 * (ib - 1)
        surfalb.albgri_col[c, ib]     = 0.20 + 0.05 * (ib - 1)
        surfalb.albsod_col[c, ib]     = 0.12 + 0.06 * (ib - 1)
        surfalb.albsoi_col[c, ib]     = 0.12 + 0.06 * (ib - 1)
        surfalb.albsnd_hst_col[c, ib] = 0.70 - 0.20 * (ib - 1)
        surfalb.albsni_hst_col[c, ib] = 0.70 - 0.20 * (ib - 1)
        # SNICAR no-aerosol ground albedos (slightly lower → positive forcing)
        surfalb.albgrd_bc_col[c, ib]  = 0.18 + 0.05 * (ib - 1)
        surfalb.albgri_bc_col[c, ib]  = 0.18 + 0.05 * (ib - 1)
        surfalb.albgrd_oc_col[c, ib]  = 0.19 + 0.05 * (ib - 1)
        surfalb.albgri_oc_col[c, ib]  = 0.19 + 0.05 * (ib - 1)
        surfalb.albgrd_dst_col[c, ib] = 0.17 + 0.05 * (ib - 1)
        surfalb.albgri_dst_col[c, ib] = 0.17 + 0.05 * (ib - 1)
        surfalb.albgrd_pur_col[c, ib] = 0.10 + 0.05 * (ib - 1)
        surfalb.albgri_pur_col[c, ib] = 0.10 + 0.05 * (ib - 1)
    end
    for p in 1:np, ib in 1:numrad
        surfalb.fabd_patch[p, ib]   = 0.15
        surfalb.fabi_patch[p, ib]   = 0.10
        surfalb.ftdd_patch[p, ib]   = 0.70
        surfalb.ftid_patch[p, ib]   = 0.05
        surfalb.ftii_patch[p, ib]   = 0.80
        surfalb.albd_patch[p, ib]   = 0.20 + 0.05 * (ib - 1)
        surfalb.albi_patch[p, ib]   = 0.20 + 0.05 * (ib - 1)
        surfalb.albdSF_patch[p, ib] = 0.18 + 0.05 * (ib - 1)
        surfalb.albiSF_patch[p, ib] = 0.18 + 0.05 * (ib - 1)
    end

    # --- SNICAR per-layer flux absorption factors (col1 has 2 snow layers) ---
    surfalb.flx_absdv_col = fill(FT(0), nc, nlev_sno1)
    surfalb.flx_absdn_col = fill(FT(0), nc, nlev_sno1)
    surfalb.flx_absiv_col = fill(FT(0), nc, nlev_sno1)
    surfalb.flx_absin_col = fill(FT(0), nc, nlev_sno1)
    # snl=-2: fortran active layers -1,0,1 → julia nlevsno-1, nlevsno, nlevsno+1
    surfalb.flx_absdv_col[1, nlevsno - 1] = 0.10
    surfalb.flx_absdv_col[1, nlevsno]     = 0.15
    surfalb.flx_absdv_col[1, nlevsno + 1] = 0.05
    surfalb.flx_absdn_col[1, nlevsno - 1] = 0.08
    surfalb.flx_absdn_col[1, nlevsno]     = 0.20
    surfalb.flx_absdn_col[1, nlevsno + 1] = 0.22
    surfalb.flx_absiv_col[1, nlevsno - 1] = 0.12
    surfalb.flx_absiv_col[1, nlevsno]     = 0.13
    surfalb.flx_absiv_col[1, nlevsno + 1] = 0.05
    surfalb.flx_absin_col[1, nlevsno - 1] = 0.10
    surfalb.flx_absin_col[1, nlevsno]     = 0.18
    surfalb.flx_absin_col[1, nlevsno + 1] = 0.22

    # --- water diagnostics (col1 snowy, col2 dry) ---
    waterdiag.snow_depth_col .= [0.5, 0.0]
    waterdiag.frac_sno_col   .= [0.9, 0.0]
    waterdiag.frac_sno_eff_col .= [0.9, 0.0]

    # --- forcing ---
    forc_solad_col = fill(FT(0), nc, numrad)
    forc_solad_col[1, 1] = 400.0; forc_solad_col[1, 2] = 200.0
    forc_solad_col[2, 1] = 300.0; forc_solad_col[2, 2] = 150.0
    forc_solai = fill(FT(0), ng, numrad)
    forc_solai[1, 1] = 100.0; forc_solai[1, 2] = 50.0

    S = (; surfalb, canopystate, solarabs, surfrad, waterdiag, col, lun, grc, pch)
    return (; np, nc, ng, nl, S,
              forc_solad_col, forc_solai, dtime = FT(3600.0), current_tod = FT(43200.0))
end

function run_pair!(S, m_np, m_up, fad, fai, dtime, current_tod, np)
    CLM.canopy_sun_shade_fracs!(S.surfalb, S.canopystate, S.solarabs,
                                fad, fai, S.pch, m_np, 1:np)
    CLM.surface_radiation!(S.surfalb, S.canopystate, S.solarabs, S.surfrad,
                           S.waterdiag, S.col, S.lun, S.grc, S.pch,
                           fad, fai, m_np, m_up, 1:np;
                           dtime = dtime, current_tod = current_tod,
                           use_subgrid_fluxes = true,
                           use_snicar_frc = USE_SNICAR, use_SSRE = USE_SSRE,
                           do_sno_oc = false)
    return nothing
end

function main(backend)
    println("=" ^ 70)
    println("END-TO-END GPU parity for canopy_sun_shade_fracs! + surface_radiation!")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU driver exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    H = build(FT)   # CPU reference
    B = build(FT)   # source for the device snapshot

    ad(x) = CLM.Adapt.adapt(device_array_type(), x)
    Sd = map(ad, B.S)
    dmask(m) = to(collect(Bool, m))

    if !(Sd.solarabs.sabg_patch isa device_array_type())
        println("  BLOCKED: a state struct did not move to the device under adapt.")
        return 2
    end

    np = H.np
    m_np = [true, false]   # patch1 non-urban
    m_up = [false, true]   # patch2 urban

    run_pair!(H.S, m_np, m_up, H.forc_solad_col, H.forc_solai,
              H.dtime, H.current_tod, np)                                  # CPU
    run_pair!(Sd, dmask(m_np), dmask(m_up), to(B.forc_solad_col), to(B.forc_solai),
              B.dtime, B.current_tod, np)                                  # device

    # Guard against a false PASS: CPU reference fields must be FINITE.
    refchecks = [
        ("sabg",        H.S.solarabs.sabg_patch),
        ("fsa",         H.S.solarabs.fsa_patch),
        ("fsr",         H.S.solarabs.fsr_patch),
        ("sabg_lyr",    H.S.solarabs.sabg_lyr_patch),
        ("parsun_z",    H.S.solarabs.parsun_z_patch),
        ("fsds_vis_d",  H.S.surfrad.fsds_vis_d_patch),
        ("sfc_frc_aer", H.S.surfrad.sfc_frc_aer_patch),
    ]
    for (nm, a) in refchecks
        if !cpu_has_finite(a)
            @printf("  BLOCKED: CPU reference field %s has no finite entry — parity would be a false PASS.\n", nm)
            return 2
        end
    end

    checks = [
        # canopy_sun_shade_fracs! outputs
        ("parsun_z",      H.S.solarabs.parsun_z_patch,        Sd.solarabs.parsun_z_patch),
        ("parsha_z",      H.S.solarabs.parsha_z_patch,        Sd.solarabs.parsha_z_patch),
        ("laisun_z",      H.S.canopystate.laisun_z_patch,     Sd.canopystate.laisun_z_patch),
        ("laisha_z",      H.S.canopystate.laisha_z_patch,     Sd.canopystate.laisha_z_patch),
        ("laisun",        H.S.canopystate.laisun_patch,       Sd.canopystate.laisun_patch),
        ("laisha",        H.S.canopystate.laisha_patch,       Sd.canopystate.laisha_patch),
        ("fsun",          H.S.canopystate.fsun_patch,         Sd.canopystate.fsun_patch),
        # surface_radiation! solarabs outputs
        ("sabv",          H.S.solarabs.sabv_patch,            Sd.solarabs.sabv_patch),
        ("sabg",          H.S.solarabs.sabg_patch,            Sd.solarabs.sabg_patch),
        ("sabg_soil",     H.S.solarabs.sabg_soil_patch,       Sd.solarabs.sabg_soil_patch),
        ("sabg_snow",     H.S.solarabs.sabg_snow_patch,       Sd.solarabs.sabg_snow_patch),
        ("sabg_pen",      H.S.solarabs.sabg_pen_patch,        Sd.solarabs.sabg_pen_patch),
        ("sabg_lyr",      H.S.solarabs.sabg_lyr_patch,        Sd.solarabs.sabg_lyr_patch),
        ("sub_surf_abs",  H.S.solarabs.sub_surf_abs_SW_patch, Sd.solarabs.sub_surf_abs_SW_patch),
        ("fsa",           H.S.solarabs.fsa_patch,             Sd.solarabs.fsa_patch),
        ("fsa_r",         H.S.solarabs.fsa_r_patch,           Sd.solarabs.fsa_r_patch),
        ("fsr",           H.S.solarabs.fsr_patch,             Sd.solarabs.fsr_patch),
        ("fsrSF",         H.S.solarabs.fsrSF_patch,           Sd.solarabs.fsrSF_patch),
        ("ssre_fsr",      H.S.solarabs.ssre_fsr_patch,        Sd.solarabs.ssre_fsr_patch),
        ("fsds_nir_d",    H.S.solarabs.fsds_nir_d_patch,      Sd.solarabs.fsds_nir_d_patch),
        ("fsds_nir_i",    H.S.solarabs.fsds_nir_i_patch,      Sd.solarabs.fsds_nir_i_patch),
        ("fsr_nir_d",     H.S.solarabs.fsr_nir_d_patch,       Sd.solarabs.fsr_nir_d_patch),
        ("fsr_nir_i",     H.S.solarabs.fsr_nir_i_patch,       Sd.solarabs.fsr_nir_i_patch),
        ("fsr_nir_d_ln",  H.S.solarabs.fsr_nir_d_ln_patch,    Sd.solarabs.fsr_nir_d_ln_patch),
        # surface_radiation! surfrad outputs
        ("fsds_vis_d",    H.S.surfrad.fsds_vis_d_patch,       Sd.surfrad.fsds_vis_d_patch),
        ("fsds_vis_i",    H.S.surfrad.fsds_vis_i_patch,       Sd.surfrad.fsds_vis_i_patch),
        ("fsds_vis_d_ln", H.S.surfrad.fsds_vis_d_ln_patch,    Sd.surfrad.fsds_vis_d_ln_patch),
        ("fsds_vis_i_ln", H.S.surfrad.fsds_vis_i_ln_patch,    Sd.surfrad.fsds_vis_i_ln_patch),
        ("fsr_vis_d",     H.S.surfrad.fsr_vis_d_patch,        Sd.surfrad.fsr_vis_d_patch),
        ("fsr_vis_i",     H.S.surfrad.fsr_vis_i_patch,        Sd.surfrad.fsr_vis_i_patch),
        ("fsr_vis_d_ln",  H.S.surfrad.fsr_vis_d_ln_patch,     Sd.surfrad.fsr_vis_d_ln_patch),
        ("parveg_ln",     H.S.surfrad.parveg_ln_patch,        Sd.surfrad.parveg_ln_patch),
        ("fsds_sno_vd",   H.S.surfrad.fsds_sno_vd_patch,      Sd.surfrad.fsds_sno_vd_patch),
        ("fsds_sno_nd",   H.S.surfrad.fsds_sno_nd_patch,      Sd.surfrad.fsds_sno_nd_patch),
        ("fsr_sno_vd",    H.S.surfrad.fsr_sno_vd_patch,       Sd.surfrad.fsr_sno_vd_patch),
        ("fsr_sno_ni",    H.S.surfrad.fsr_sno_ni_patch,       Sd.surfrad.fsr_sno_ni_patch),
        # SNICAR forcings
        ("sfc_frc_bc",    H.S.surfrad.sfc_frc_bc_patch,       Sd.surfrad.sfc_frc_bc_patch),
        ("sfc_frc_oc",    H.S.surfrad.sfc_frc_oc_patch,       Sd.surfrad.sfc_frc_oc_patch),
        ("sfc_frc_dst",   H.S.surfrad.sfc_frc_dst_patch,      Sd.surfrad.sfc_frc_dst_patch),
        ("sfc_frc_aer",   H.S.surfrad.sfc_frc_aer_patch,      Sd.surfrad.sfc_frc_aer_patch),
        ("sfc_frc_bc_sno",H.S.surfrad.sfc_frc_bc_sno_patch,   Sd.surfrad.sfc_frc_bc_sno_patch),
        ("sfc_frc_aer_sno",H.S.surfrad.sfc_frc_aer_sno_patch, Sd.surfrad.sfc_frc_aer_sno_patch),
        # SSRE diagnostics
        ("ssre_fsr_vis_d",H.S.surfrad.ssre_fsr_vis_d_patch,   Sd.surfrad.ssre_fsr_vis_d_patch),
        # fsun zeroed for urban
        ("fsun_zero_urb", H.S.canopystate.fsun_patch,         Sd.canopystate.fsun_patch),
    ]
    nfail = 0; worst = 0.0
    for (nm, a, b) in checks
        d = reldiff(a, b); ok = d < 1f-3
        worst = max(worst, d)
        @printf("  [%s] %-17s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
        ok || (nfail += 1)
    end
    @printf("\n  worst reldiff = %.3e\n", worst)
    println(nfail == 0 ? "  WHOLE radiation pair MATCHES CPU ON $name ($FT)" :
                         "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
