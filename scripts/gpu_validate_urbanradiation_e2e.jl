# ==========================================================================
# gpu_validate_urbanradiation_e2e.jl — end-to-end GPU parity for the WHOLE
# urban canyon radiation function urban_radiation! (which internally drives the
# kernelized net_longwave! canyon longwave solve).
#
# Builds a small Float32 instance with TWO landunits:
#   lun 1: URBAN landunit with 5 columns / 5 patches — one per canyon facet
#          (roof, sunwall, shadewall, impervious road, pervious road). This
#          exercises the snow-effect input-forcing pass (roof/road snow branch),
#          the full multiple-reflection canyon longwave solve, the roof longwave,
#          and every branch of the per-patch history/atm-output pass.
#   lun 2: NON-URBAN soil landunit (1 column, 0 patches) — drives the SPVAL
#          restart-field fill pass for non-urban landunits.
#
# Two scenes are validated:
#   DAY:   nonzero direct+diffuse solar (drives sabg / fsa / fsa_u).
#   NIGHT: zero solar (sabg -> 0; still finite) — exercises the zero-sun branch
#          while the longwave canyon solve still runs.
#
# Runs urban_radiation! WHOLE on the CPU, adapts every state struct (+ masks /
# forcing) to the GPU, runs the SAME call on the device, and compares every
# mutated field with reldiff.
#
# CRITICAL: every required param/input is a real CLM default so the CPU
# reference is FINITE — reldiff silently PASSES when both sides are NaN, so the
# harness ASSERTS the CPU reference fields are finite before trusting parity.
#
#   julia --project=scripts scripts/gpu_validate_urbanradiation_e2e.jl
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

function build(::Type{FT}; daytime::Bool) where {FT}
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
    numrad = CLM.NUMRAD
    nl = 2          # lun1 urban, lun2 non-urban soil
    nc = 6          # cols 1..5 urban facets, col 6 non-urban soil
    np = 5          # 5 urban patches (one per facet)
    ng = 1

    hwr = 1.0
    wtroad_perv = 0.4

    # --- LandunitData ---
    lun = CLM.LandunitData{FT}(); CLM.landunit_init!(lun, nl)
    lun.gridcell .= [1, 1]
    lun.coli     .= [1, 6]
    lun.colf     .= [5, 6]
    lun.canyon_hwr  .= [FT(hwr), FT(0)]
    lun.wtroad_perv .= [FT(wtroad_perv), FT(0)]
    lun.urbpoi   .= [true, false]
    lun.itype    .= [CLM.ISTURB_MIN, CLM.ISTSOIL]

    # --- ColumnData ---
    col = CLM.ColumnData{FT}(); CLM.column_init!(col, nc)
    coltypes = [CLM.ICOL_ROOF, CLM.ICOL_SUNWALL, CLM.ICOL_SHADEWALL,
                CLM.ICOL_ROAD_IMPERV, CLM.ICOL_ROAD_PERV, 0]
    for c in 1:nc
        col.itype[c]    = coltypes[c]
        col.landunit[c] = c <= 5 ? 1 : 2
        col.snl[c]      = 0
    end

    # --- PatchData (5 urban patches, one per facet column) ---
    pch = CLM.PatchData{FT}(); CLM.patch_init!(pch, np)
    for p in 1:np
        pch.column[p]   = p
        pch.landunit[p] = 1
        pch.gridcell[p] = 1
    end

    # --- UrbanParamsData ---
    up = CLM.UrbanParamsData{FT}(); CLM.urbanparams_init!(up, nl; numrad = numrad)
    up.em_roof[1]    = 0.90
    up.em_improad[1] = 0.95
    up.em_perroad[1] = 0.95
    up.em_wall[1]    = 0.85
    up.vf_sr[1] = sqrt(hwr^2 + 1.0) - hwr
    up.vf_wr[1] = 0.5 * (1.0 - up.vf_sr[1])
    up.vf_sw[1] = 0.5 * (hwr + 1.0 - sqrt(hwr^2 + 1.0)) / hwr
    up.vf_rw[1] = up.vf_sw[1]
    up.vf_ww[1] = 1.0 - up.vf_sw[1] - up.vf_rw[1]

    # --- TemperatureData ---
    temperature = CLM.TemperatureData{FT}(); CLM.temperature_init!(temperature, np, nc, nl, ng)
    # distinct facet ground temperatures so each branch is exercised
    temperature.t_grnd_col[1] = 300.0   # roof
    temperature.t_grnd_col[2] = 295.0   # sunwall
    temperature.t_grnd_col[3] = 292.0   # shadewall
    temperature.t_grnd_col[4] = 288.0   # imperv road
    temperature.t_grnd_col[5] = 287.0   # perv road
    temperature.t_grnd_col[6] = 285.0   # non-urban soil

    # --- SolarAbsorbedData ---
    solarabs = CLM.SolarAbsorbedData{FT}(); CLM.solarabs_init!(solarabs, np, nl)
    CLM.solarabs_init_cold!(solarabs, 1:nl)
    for ib in 1:numrad
        solarabs.sabs_roof_dir_lun[1, ib]      = 0.80
        solarabs.sabs_roof_dif_lun[1, ib]      = 0.80
        solarabs.sabs_sunwall_dir_lun[1, ib]   = 0.30
        solarabs.sabs_sunwall_dif_lun[1, ib]   = 0.30
        solarabs.sabs_shadewall_dir_lun[1, ib] = 0.20
        solarabs.sabs_shadewall_dif_lun[1, ib] = 0.20
        solarabs.sabs_improad_dir_lun[1, ib]   = 0.70
        solarabs.sabs_improad_dif_lun[1, ib]   = 0.70
        solarabs.sabs_perroad_dir_lun[1, ib]   = 0.65
        solarabs.sabs_perroad_dif_lun[1, ib]   = 0.65
    end

    # --- EnergyFluxData ---
    energyflux = CLM.EnergyFluxData{FT}(); CLM.energyflux_init!(energyflux, np, nc, nl, ng)

    # --- WaterDiagnosticBulkData (snow on roof + imperv road) ---
    waterdiag = CLM.WaterDiagnosticBulkData{FT}(); CLM.waterdiagnosticbulk_init!(waterdiag, nc, np, nl, ng)
    waterdiag.frac_sno_col .= FT(0)
    waterdiag.frac_sno_col[1] = 0.5   # roof snow  -> em_roof_s snow branch
    waterdiag.frac_sno_col[4] = 0.3   # imperv road snow -> em_improad_s snow branch
    waterdiag.frac_sno_col[5] = 0.2   # perv road snow   -> em_perroad_s snow branch

    # --- Forcing ---
    forc_lwrad = fill(FT(350.0), ng)
    forc_solad = fill(FT(0), ng, numrad)
    forc_solai = fill(FT(0), ng, numrad)
    if daytime
        forc_solad[1, 1] = 200.0; forc_solad[1, 2] = 100.0
        forc_solai[1, 1] = 80.0;  forc_solai[1, 2] = 40.0
    end

    # --- Masks ---
    mask_nourbanl = falses(nl); mask_nourbanl[2] = true
    mask_urbanl   = falses(nl); mask_urbanl[1] = true
    mask_urbanc   = falses(nc); mask_urbanc[1:5] .= true
    mask_urbanp   = falses(np); mask_urbanp[1:np] .= true

    S = (; solarabs, energyflux, col, lun, pch, urbanparams = up,
          temperature, waterdiag)
    return (; np, nc, nl, ng, S, forc_lwrad, forc_solad, forc_solai,
              mask_nourbanl, mask_urbanl, mask_urbanc, mask_urbanp)
end

function run!(S, fl, fad, fai, m_nourbanl, m_urbanl, m_urbanc, m_urbanp, nl, np)
    CLM.urban_radiation!(S.solarabs, S.energyflux, S.col, S.lun, S.pch,
                         S.urbanparams, S.temperature, S.waterdiag,
                         fl, fad, fai,
                         m_nourbanl, m_urbanl, m_urbanc, m_urbanp,
                         1:nl, 1:np)
    return nothing
end

function validate(name, to, FT, daytime)
    H = build(FT; daytime = daytime)   # CPU reference
    B = build(FT; daytime = daytime)   # source for the device snapshot

    ad(x) = CLM.Adapt.adapt(device_array_type(), x)
    Sd = map(ad, B.S)
    dmask(m) = to(collect(Bool, m))

    if !(Sd.solarabs.sabg_patch isa device_array_type())
        println("  BLOCKED: a state struct did not move to the device under adapt.")
        return 2
    end

    run!(H.S, H.forc_lwrad, H.forc_solad, H.forc_solai,
         H.mask_nourbanl, H.mask_urbanl, H.mask_urbanc, H.mask_urbanp, H.nl, H.np)
    run!(Sd, to(B.forc_lwrad), to(B.forc_solad), to(B.forc_solai),
         dmask(B.mask_nourbanl), dmask(B.mask_urbanl), dmask(B.mask_urbanc),
         dmask(B.mask_urbanp), B.nl, B.np)

    # Guard against false PASS: CPU reference fields must be FINITE.
    refchecks = [
        ("eflx_lwrad_out", H.S.energyflux.eflx_lwrad_out_patch),
        ("eflx_lwrad_net", H.S.energyflux.eflx_lwrad_net_patch),
        ("sabg",           H.S.solarabs.sabg_patch),
        ("fsa",            H.S.solarabs.fsa_patch),
        ("sabs_roof_dir(non-urb SPVAL)", view(H.S.solarabs.sabs_roof_dir_lun, 2, :)),
    ]
    for (nm, a) in refchecks
        if !cpu_has_finite(a)
            @printf("  BLOCKED [%s]: CPU reference field %s has no finite entry — parity would be a false PASS.\n",
                    daytime ? "DAY" : "NIGHT", nm)
            return 2
        end
    end

    checks = [
        ("eflx_lwrad_out",   H.S.energyflux.eflx_lwrad_out_patch,   Sd.energyflux.eflx_lwrad_out_patch),
        ("eflx_lwrad_net",   H.S.energyflux.eflx_lwrad_net_patch,   Sd.energyflux.eflx_lwrad_net_patch),
        ("eflx_lwrad_net_u", H.S.energyflux.eflx_lwrad_net_u_patch, Sd.energyflux.eflx_lwrad_net_u_patch),
        ("sabg",             H.S.solarabs.sabg_patch,               Sd.solarabs.sabg_patch),
        ("sabv",             H.S.solarabs.sabv_patch,               Sd.solarabs.sabv_patch),
        ("fsa",              H.S.solarabs.fsa_patch,                Sd.solarabs.fsa_patch),
        ("fsa_u",            H.S.solarabs.fsa_u_patch,              Sd.solarabs.fsa_u_patch),
        ("sabs_roof_dir",    H.S.solarabs.sabs_roof_dir_lun,        Sd.solarabs.sabs_roof_dir_lun),
        ("sabs_roof_dif",    H.S.solarabs.sabs_roof_dif_lun,        Sd.solarabs.sabs_roof_dif_lun),
        ("sabs_perroad_dif", H.S.solarabs.sabs_perroad_dif_lun,     Sd.solarabs.sabs_perroad_dif_lun),
    ]
    nfail = 0; worst = 0.0
    @printf("  --- %s scene ---\n", daytime ? "DAY (nonzero solar)" : "NIGHT (zero solar)")
    for (nm, a, b) in checks
        d = reldiff(a, b); ok = d < 1f-3
        worst = max(worst, d)
        @printf("    [%s] %-18s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
        ok || (nfail += 1)
    end
    @printf("    worst reldiff (%s) = %.3e\n", daytime ? "DAY" : "NIGHT", worst)
    return nfail
end

function main(backend)
    println("=" ^ 70)
    println("END-TO-END GPU parity for urban_radiation! (canyon longwave + solar)")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU driver exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    r_day   = validate(name, to, FT, true)
    (r_day == 2) && return 2
    r_night = validate(name, to, FT, false)
    (r_night == 2) && return 2

    nfail = r_day + r_night
    println()
    println(nfail == 0 ? "  WHOLE urban_radiation! MATCHES CPU ON $name ($FT) — day + night" :
                         "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
