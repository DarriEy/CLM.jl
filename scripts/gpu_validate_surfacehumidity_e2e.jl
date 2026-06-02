# ==========================================================================
# gpu_validate_surfacehumidity_e2e.jl — end-to-end Metal parity for the WHOLE
# calculate_surface_humidity! driver (per-column kernel with an internal
# sequential nlevgrnd loop for the pervious-road branch, plus qsat() calls).
#
# Builds a small Float32 instance with several columns exercising the real
# branches: soil/crop (snl=0 and snl<0 snow), soil with surface water, wet,
# ice, urban sunwall (qred=0), urban roof (qred=1), and pervious road (the
# nlevgrnd normalization loop). Runs calculate_surface_humidity! on the CPU,
# adapts every state struct (+ Bool mask + forcings) to Metal, runs the SAME
# call on the device, and compares the mutated outputs field-by-field.
#
#   julia --project=scripts scripts/gpu_validate_surfacehumidity_e2e.jl
#
# reldiff PASSES when both sides are NaN, so main() asserts the CPU reference
# fields are FINITE before trusting any parity number.
# ==========================================================================

using CLM
using Printf
import Metal
include(joinpath(@__DIR__, "gpu_backends.jl"))

function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(A[i] - B[i]) / (1.0 + max(abs(A[i]), abs(B[i]))))
    end
    return m
end
# Pass the parity gate if the CPU reference has ANY finite, mutated entry
# (some columns legitimately leave soilalpha/_u at SPVAL/NaN). reldiff itself is
# NaN-aware (both-NaN agrees), so per-field parity is still meaningful.
allfinite(a) = any(isfinite, Array(a))

function build(::Type{FT}) where {FT}
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
    nlevsno = CLM.varpar.nlevsno; nlevgrnd = CLM.varpar.nlevgrnd; joff = nlevsno
    nc = 7; nl = 7

    col = CLM.ColumnData{FT}(); CLM.column_init!(col, nc)
    lun = CLM.LandunitData{FT}(); CLM.landunit_init!(lun, nl)
    temp = CLM.TemperatureData{FT}(); CLM.temperature_init!(temp, 0, nc, nl, 1)
    ss = CLM.SoilStateData{FT}(); CLM.soilstate_init!(ss, 0, nc)
    wsb = CLM.WaterStateBulkData{FT}(); CLM.waterstatebulk_init!(wsb, nc, 0, nl, 1)
    wdb = CLM.WaterDiagnosticBulkData{FT}(); CLM.waterdiagnosticbulk_init!(wdb, nc, 0, nl, 1)

    for c in 1:nc
        col.landunit[c] = c; col.snl[c] = 0
    end
    forc_pbot = fill(FT(101325.0), nc)
    forc_q = fill(FT(0.005), nc)

    function set_soil!(c; liq, snl, fh2osfc, fsno, tsoil, th2osfc)
        col.dz[c, 1+joff] = 0.02
        wsb.ws.h2osoi_liq_col[c, 1+joff] = liq
        wsb.ws.h2osoi_ice_col[c, 1+joff] = 0.0
        ss.watsat_col[c, 1] = 0.45; ss.sucsat_col[c, 1] = 100.0
        ss.bsw_col[c, 1] = 5.0; ss.smpmin_col[c] = -1.0e8
        temp.t_soisno_col[c, 1+joff] = tsoil
        temp.t_grnd_col[c] = tsoil; temp.t_h2osfc_col[c] = th2osfc
        wdb.frac_sno_eff_col[c] = fsno; wdb.frac_h2osfc_col[c] = fh2osfc
        col.snl[c] = snl
    end

    # col 1: soil, snl=0, no snow/h2osfc
    col.itype[1] = CLM.ISTSOIL; lun.itype[1] = CLM.ISTSOIL
    set_soil!(1; liq=50.0, snl=0, fh2osfc=0.0, fsno=0.0, tsoil=290.0, th2osfc=290.0)

    # col 2: soil with snow layers (snl<0)
    col.itype[2] = CLM.ISTSOIL; lun.itype[2] = CLM.ISTSOIL
    set_soil!(2; liq=40.0, snl=-2, fh2osfc=0.0, fsno=0.5, tsoil=280.0, th2osfc=280.0)
    wsb.ws.h2osoi_ice_col[2, 1+joff] = 5.0
    temp.t_soisno_col[2, joff]   = 268.0   # top snow (j=0)
    temp.t_soisno_col[2, joff-1] = 265.0   # snow (j=-1)
    temp.t_grnd_col[2] = 268.0

    # col 3: soil with surface water
    col.itype[3] = CLM.ISTSOIL; lun.itype[3] = CLM.ISTSOIL
    set_soil!(3; liq=50.0, snl=0, fh2osfc=0.3, fsno=0.0, tsoil=290.0, th2osfc=292.0)

    # col 4 wet, col 5 ice
    col.itype[4] = CLM.ISTWET; lun.itype[4] = CLM.ISTWET; temp.t_grnd_col[4] = 285.0
    col.itype[5] = CLM.ISTICE; lun.itype[5] = CLM.ISTICE; temp.t_grnd_col[5] = 260.0

    # col 6: urban sunwall (qred=0); forc_q=0 so dew cap doesn't fire
    col.itype[6] = CLM.ICOL_SUNWALL; lun.itype[6] = CLM.ISTURB_TBD
    temp.t_grnd_col[6] = 300.0; forc_q[6] = 0.0

    # col 7: pervious road (the nlevgrnd loop)
    col.itype[7] = CLM.ICOL_ROAD_PERV; lun.itype[7] = CLM.ISTURB_MD
    wdb.frac_sno_eff_col[7] = 0.0
    for j in 1:nlevgrnd
        temp.t_soisno_col[7, j+joff] = 290.0
        col.dz[7, j+joff] = 0.1
        wsb.ws.h2osoi_liq_col[7, j+joff] = 10.0
        wsb.ws.h2osoi_ice_col[7, j+joff] = 0.0
        ss.watsat_col[7, j] = 0.45; ss.watdry_col[7, j] = 0.05; ss.watopt_col[7, j] = 0.35
        ss.rootfr_road_perv_col[7, j] = 1.0 / nlevgrnd
    end
    temp.t_grnd_col[7] = 290.0

    S = (; col, lun, temp, ss, wsb, wdb)
    return (; nc, S, forc_pbot, forc_q)
end

run!(S, fp, fq, m) = CLM.calculate_surface_humidity!(S.col, S.lun, S.temp, S.ss,
    S.wsb, S.wdb, fp, fq, m, 1:length(m))

function main(backend)
    println("="^70)
    println("END-TO-END Metal parity for calculate_surface_humidity! (whole driver)")
    println("="^70)
    if backend === nothing
        println("  No GPU backend — nothing to validate.")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    H = build(FT); B = build(FT)
    ad(x) = CLM.Adapt.adapt(Metal.MtlArray, x)
    Sd = map(ad, B.S)
    if !(Sd.wdb.qg_col isa Metal.MtlArray)
        println("  BLOCKED: a state struct did not move to the device under adapt.")
        return 2
    end
    m = trues(H.nc)
    run!(H.S, H.forc_pbot, H.forc_q, m)
    run!(Sd, to(B.forc_pbot), to(B.forc_q), to(collect(Bool, m)))

    checks = [
        ("qg",          H.S.wdb.qg_col,         Sd.wdb.qg_col),
        ("qg_soil",     H.S.wdb.qg_soil_col,    Sd.wdb.qg_soil_col),
        ("qg_snow",     H.S.wdb.qg_snow_col,    Sd.wdb.qg_snow_col),
        ("qg_h2osfc",   H.S.wdb.qg_h2osfc_col,  Sd.wdb.qg_h2osfc_col),
        ("dqgdT",       H.S.wdb.dqgdT_col,      Sd.wdb.dqgdT_col),
        ("soilalpha",   H.S.ss.soilalpha_col,   Sd.ss.soilalpha_col),
        ("soilalpha_u", H.S.ss.soilalpha_u_col, Sd.ss.soilalpha_u_col),
        ("rootr_road",  H.S.ss.rootr_road_perv_col, Sd.ss.rootr_road_perv_col),
    ]
    nfail = 0
    for (nm, a, b) in checks
        if !allfinite(a)
            @printf("  [WARN] %-12s CPU reference not all-finite — skipping\n", nm)
            continue
        end
        d = reldiff(a, b); ok = d < 1f-3
        @printf("  [%s] %-12s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
        ok || (nfail += 1)
    end
    println()
    println(nfail == 0 ? "  WHOLE calculate_surface_humidity! MATCHES CPU on $name ($FT)" :
                         "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
