# ==========================================================================
# gpu_validate_soilwater_mf_e2e.jl — end-to-end GPU parity for the WHOLE
# moisture-form Richards solver (soilwater_moisture_form!, the non-default
# soil_water! path used when use_aquifer_layer=false → MOISTURE_FORM + BC_ZERO_FLUX).
#
# Exercises the full per-column adaptive-substep kernel (hydraulic properties,
# boundary fluxes, RHS/LHS, in-thread tridiagonal solve, substep accept/reject)
# end-to-end: builds a Float32 instance, runs soil_water! on CPU, adapts every
# struct to the GPU, runs the SAME call on device, compares mutated outputs.
#
#   julia --project=scripts scripts/gpu_validate_soilwater_mf_e2e.jl
# ==========================================================================

using CLM
using Printf
include(joinpath(@__DIR__, "gpu_backends.jl"))

const NS = 5; const NG = 10; const NSOI = 10; const NMAX = 10
const JOFF = NS

function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(A[i] - B[i]) / (1.0 + max(abs(A[i]), abs(B[i]))))
    end
    return m
end

function build(::Type{FT}) where {FT}
    nc = 4; np = 4; nl = 4; ng = 1
    col = CLM.ColumnData{FT}();         CLM.column_init!(col, nc)
    temp = CLM.TemperatureData{FT}();   CLM.temperature_init!(temp, np, nc, nl, ng)
    ef = CLM.EnergyFluxData{FT}();      CLM.energyflux_init!(ef, np, nc, nl, ng)
    ss = CLM.SoilStateData{FT}();       CLM.soilstate_init!(ss, np, nc)
    wsb = CLM.WaterStateBulkData{FT}(); CLM.waterstatebulk_init!(wsb, nc, np, nl, ng)
    wfb = CLM.WaterFluxBulkData{FT}();  CLM.waterfluxbulk_init!(wfb, nc, np, nl, ng)
    cs = CLM.CanopyStateData{FT}();     CLM.canopystate_init!(cs, np)
    sh = CLM.SoilHydrologyData{FT}();   CLM.soilhydrology_init!(sh, nc)
    return (; nc, np, nl, ng, col, temp, ef, ss, wsb, wfb, cs, sh)
end

function populate!(B, mask_hydrology, active)
    col, temp, ss, wsb, wfb, sh = B.col, B.temp, B.ss, B.wsb, B.wfb, B.sh
    DENH2O = CLM.DENH2O
    for c in active
        mask_hydrology[c] = true
        col.landunit[c] = c; col.itype[c] = CLM.ISTSOIL
        col.snl[c] = 0; col.nbedrock[c] = NSOI
        for j in 1:NMAX
            jj = j + JOFF
            col.dz[c, jj] = 0.1 * j
            col.z[c, jj] = sum(0.1 * k for k in 1:j) - 0.5 * 0.1 * j
        end
        col.zi[c, JOFF + 1] = 0.0
        for j in 1:NMAX; col.zi[c, j + JOFF + 1] = col.zi[c, j + JOFF] + col.dz[c, j + JOFF]; end
        for j in 1:NG
            ss.watsat_col[c, j] = 0.45; ss.bsw_col[c, j] = 5.0
            ss.sucsat_col[c, j] = 100.0; ss.hksat_col[c, j] = 1.0e-4
            ss.eff_porosity_col[c, j] = 0.45; ss.smp_l_col[c, j] = 0.0; ss.hk_l_col[c, j] = 0.0
            sh.icefrac_col[c, j] = 0.0
        end
        ss.smpmin_col[c] = -1.0e8
        sh.zwt_col[c] = 6.0; sh.qcharge_col[c] = 0.0; sh.num_substeps_col[c] = 0.0
        # Wet column with a strong infiltration so the substep error control engages.
        for j in 1:NMAX
            jj = j + JOFF
            wsb.ws.h2osoi_liq_col[c, jj] = 20.0 + 2.0 * j
            wsb.ws.h2osoi_ice_col[c, jj] = 0.0
            wsb.ws.h2osoi_vol_col[c, j]  = wsb.ws.h2osoi_liq_col[c, jj] / (col.dz[c, jj] * DENH2O)
        end
        wsb.ws.wa_col[c] = 4000.0
        wfb.wf.qflx_infl_col[c] = 5.0e-2
        for j in 1:NSOI; wfb.qflx_rootsoi_col[c, j] = 1.0e-3; end
    end
    return nothing
end

function run_sw!(B, m_h, m_u, cfg, dtime)
    CLM.soil_water!(B.col, m_h, m_u, B.sh, B.ss, B.wfb, B.wsb, B.temp, B.cs, B.ef,
        CLM.SoilWaterRetentionCurveClappHornberg1978(), cfg, dtime)
end

function main(backend)
    println("=" ^ 70)
    println("END-TO-END GPU parity for soilwater_moisture_form! (whole solver)")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU solver exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    vp = CLM.varpar
    saved = (vp.nlevsno, vp.nlevgrnd, vp.nlevsoi, vp.nlevmaxurbgrnd)
    vp.nlevsno, vp.nlevgrnd, vp.nlevsoi, vp.nlevmaxurbgrnd = NS, NG, NSOI, NMAX
    try
        B = build(FT)
        m_h = falses(B.nc); m_u = falses(B.nc)
        populate!(B, m_h, [1, 2])
        # Non-default path: moisture form with zero-flux lower boundary.
        cfg = CLM.SoilWaterMovementConfig(soilwater_movement_method = CLM.MOISTURE_FORM,
                                          lower_boundary_condition = CLM.BC_ZERO_FLUX)
        dt = FT(1800)

        ad(x) = CLM.Adapt.adapt(device_array_type(), x)
        D = (; nc = B.nc, col = ad(B.col), temp = ad(B.temp), ef = ad(B.ef), ss = ad(B.ss),
             wsb = ad(B.wsb), wfb = ad(B.wfb), cs = ad(B.cs), sh = ad(B.sh))
        dmask(m) = to(collect(Bool, m))

        if !(D.sh.zwt_col isa device_array_type() && D.wsb.ws.h2osoi_liq_col isa device_array_type())
            println("  BLOCKED: a state struct did not move to the device under adapt.")
            return 2
        end

        run_sw!(B, m_h, m_u, cfg, dt)                       # CPU reference
        run_sw!(D, dmask(m_h), dmask(m_u), cfg, dt)         # device

        checks = [
            ("h2osoi_liq",  B.wsb.ws.h2osoi_liq_col, D.wsb.ws.h2osoi_liq_col),
            ("smp_l",       B.ss.smp_l_col,          D.ss.smp_l_col),
            ("hk_l",        B.ss.hk_l_col,           D.ss.hk_l_col),
            ("qcharge",     B.sh.qcharge_col,        D.sh.qcharge_col),
            ("num_substeps",B.sh.num_substeps_col,   D.sh.num_substeps_col),
        ]
        nfail = 0
        for (nm, a, b) in checks
            d = reldiff(a, b)
            ok = d < 1f-2
            @printf("  [%s] %-13s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
            ok || (nfail += 1)
        end
        ns = Array(B.sh.num_substeps_col)
        @printf("\n  (CPU substeps for active cols: %s)\n",
                string([isnan(x) ? -1 : Int(round(x)) for x in ns]))
        println()
        println(nfail == 0 ? "  WHOLE soilwater_moisture_form! MATCHES CPU ON $name ($FT) ✓" :
                             "  DIVERGENCE — investigate.")
        return nfail == 0 ? 0 : 1
    finally
        vp.nlevsno, vp.nlevgrnd, vp.nlevsoi, vp.nlevmaxurbgrnd = saved
    end
end

const BACKEND = detect_backend()
exit(main(BACKEND))
