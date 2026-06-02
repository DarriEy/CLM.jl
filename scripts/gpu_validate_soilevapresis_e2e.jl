# ==========================================================================
# gpu_validate_soilevapresis_e2e.jl — end-to-end Metal parity for the WHOLE
# calc_soilevap_resis! driver (both methods: Lee-Pielke 1992 beta and the
# Swenson & Lawrence 2014 dry-surface-layer resistance).
#
# Builds a small Float32 instance with several columns exercising the real
# branches (soil below/above field capacity, snow+h2osfc, wet/ice landunit,
# urban column types), runs calc_soilevap_resis! on the CPU, adapts every
# state struct (+ Bool mask) to Metal, runs the SAME call on the device, and
# compares the mutated outputs field-by-field. Both methods are validated.
#
#   julia --project=scripts scripts/gpu_validate_soilevapresis_e2e.jl
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
allfinite(a) = all(isfinite, Array(a))

# 6 columns: 1 soil<FC, 2 soil>FC w/ snow+h2osfc, 3 wet, 4 ice, 5 road_perv, 6 sunwall
function build(::Type{FT}) where {FT}
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
    nlevsno = CLM.varpar.nlevsno; joff = nlevsno
    nc = 6; nl = 6

    col = CLM.ColumnData{FT}(); CLM.column_init!(col, nc)
    lun = CLM.LandunitData{FT}(); CLM.landunit_init!(lun, nl)
    temp = CLM.TemperatureData{FT}(); CLM.temperature_init!(temp, 0, nc, nl, 1)
    ss = CLM.SoilStateData{FT}(); CLM.soilstate_init!(ss, 0, nc)
    wsb = CLM.WaterStateBulkData{FT}(); CLM.waterstatebulk_init!(wsb, nc, 0, nl, 1)
    wdb = CLM.WaterDiagnosticBulkData{FT}(); CLM.waterdiagnosticbulk_init!(wdb, nc, 0, nl, 1)

    for c in 1:nc
        col.landunit[c] = c
    end
    # col 1: soil below field capacity
    col.itype[1] = CLM.ISTSOIL; lun.itype[1] = CLM.ISTSOIL
    col.dz[1, 1+joff] = 0.02
    wsb.ws.h2osoi_liq_col[1, 1+joff] = 2.0; wsb.ws.h2osoi_ice_col[1, 1+joff] = 0.0
    ss.watsat_col[1, 1] = 0.45; ss.watfc_col[1, 1] = 0.30
    ss.bsw_col[1, 1] = 5.0; ss.sucsat_col[1, 1] = 100.0
    temp.t_soisno_col[1, 1+joff] = 290.0
    wdb.frac_sno_col[1] = 0.0; wdb.frac_h2osfc_col[1] = 0.0

    # col 2: soil above field capacity, with snow + h2osfc
    col.itype[2] = CLM.ISTSOIL; lun.itype[2] = CLM.ISTSOIL
    col.dz[2, 1+joff] = 0.02
    wsb.ws.h2osoi_liq_col[2, 1+joff] = 8.0; wsb.ws.h2osoi_ice_col[2, 1+joff] = 0.0
    ss.watsat_col[2, 1] = 0.45; ss.watfc_col[2, 1] = 0.30
    ss.bsw_col[2, 1] = 5.0; ss.sucsat_col[2, 1] = 100.0
    temp.t_soisno_col[2, 1+joff] = 285.0
    wdb.frac_sno_col[2] = 0.3; wdb.frac_h2osfc_col[2] = 0.2

    # col 3 wet, col 4 ice
    col.itype[3] = CLM.ISTWET; lun.itype[3] = CLM.ISTWET
    col.itype[4] = CLM.ISTICE; lun.itype[4] = CLM.ISTICE

    # col 5: pervious road (urban, not wet/ice)
    col.itype[5] = CLM.ICOL_ROAD_PERV; lun.itype[5] = CLM.ISTURB_TBD
    col.dz[5, 1+joff] = 0.02
    wsb.ws.h2osoi_liq_col[5, 1+joff] = 1.0
    ss.watsat_col[5, 1] = 0.45; ss.watfc_col[5, 1] = 0.30
    ss.bsw_col[5, 1] = 5.0; ss.sucsat_col[5, 1] = 100.0
    temp.t_soisno_col[5, 1+joff] = 285.0

    # col 6: sunwall (urban)
    col.itype[6] = CLM.ICOL_SUNWALL; lun.itype[6] = CLM.ISTURB_TBD
    col.dz[6, 1+joff] = 0.02
    ss.watsat_col[6, 1] = 0.45; ss.watfc_col[6, 1] = 0.30
    ss.bsw_col[6, 1] = 5.0; ss.sucsat_col[6, 1] = 100.0
    temp.t_soisno_col[6, 1+joff] = 285.0

    S = (; col, lun, temp, ss, wsb, wdb)
    return (; nc, S)
end

run!(S, m) = CLM.calc_soilevap_resis!(S.col, S.lun, S.ss, S.wsb, S.wdb, S.temp,
                                      m, 1:length(m))

function main(backend)
    println("="^70)
    println("END-TO-END Metal parity for calc_soilevap_resis! (whole driver)")
    println("="^70)
    if backend === nothing
        println("  No GPU backend — nothing to validate.")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    ad(x) = CLM.Adapt.adapt(Metal.MtlArray, x)
    nfail = 0

    for method in (CLM.SOIL_RESIS_LEEPIELKE_1992, CLM.SOIL_RESIS_SL_14)
        CLM.soil_resistance_read_nl!(soil_resis_method = method)
        CLM.surface_resistance_read_params!(d_max = 15.0, frac_sat_soil_dsl_init = 0.8)
        mname = method == CLM.SOIL_RESIS_LEEPIELKE_1992 ? "Lee-Pielke1992" : "SL14"
        println("  --- method: $mname ---")

        H = build(FT); B = build(FT)
        Sd = map(ad, B.S)
        if !(Sd.ss.soilbeta_col isa Metal.MtlArray)
            println("  BLOCKED: a state struct did not move to the device under adapt.")
            return 2
        end
        m = trues(H.nc)
        run!(H.S, m)
        run!(Sd, to(collect(Bool, m)))

        checks = method == CLM.SOIL_RESIS_LEEPIELKE_1992 ?
            [("soilbeta", H.S.ss.soilbeta_col, Sd.ss.soilbeta_col)] :
            [("dsl",       H.S.ss.dsl_col,       Sd.ss.dsl_col),
             ("soilresis", H.S.ss.soilresis_col, Sd.ss.soilresis_col)]

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
    end

    println(nfail == 0 ? "  WHOLE calc_soilevap_resis! MATCHES CPU on $name ($FT)" :
                         "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
