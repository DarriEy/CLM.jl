# ==========================================================================
# gpu_validate_soiltemp_e2e.jl — end-to-end Metal parity for the WHOLE
# soil_temperature! driver (not just its individual kernels).
#
# Builds a small Float32 instance (one soil column), runs soil_temperature! on the CPU,
# adapts every state struct to the GPU + moves masks/forcing, runs the SAME call on the
# device, and compares the mutated outputs. This exercises the full chain — soil_therm_prop,
# ground heat flux + patch->column scatter, heat diffusion, RHS/matrix assembly, the batched
# pentadiagonal solve, phase change, and the orchestrator's own kernelized loops — together.
#
#   julia --project=scripts scripts/gpu_validate_soiltemp_e2e.jl
# ==========================================================================

using CLM
using Printf
import Metal   # this e2e is Metal-specific; MtlArray is the Adapt adaptor type for the structs
include(joinpath(@__DIR__, "gpu_backends.jl"))

const NS = 5; const NG = 10; const NU = 5; const NMAX = 10; const NSOI = 8
const JOFF = NS

# NaN-aware: inactive (unmasked) columns hold NaN scratch on BOTH backends (rvector/tvector
# init), so both-NaN counts as agreement; a one-sided NaN stays NaN and flags real divergence.
function maxabsdiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(A[i] - B[i]))
    end
    return m
end

# Build the 12 state structs (+ masks/forcing/bounds) at precision FT, sizes nc=np=nl=4.
function build(::Type{FT}) where {FT}
    nc = np = nl = 4; ng = 1
    col = CLM.ColumnData{FT}();              CLM.column_init!(col, nc)
    lun = CLM.LandunitData{FT}();            CLM.landunit_init!(lun, nl)
    patch = CLM.PatchData{FT}();             CLM.patch_init!(patch, np)
    temp = CLM.TemperatureData{FT}();        CLM.temperature_init!(temp, np, nc, nl, ng)
    ef = CLM.EnergyFluxData{FT}();           CLM.energyflux_init!(ef, np, nc, nl, ng)
    ss = CLM.SoilStateData{FT}();            CLM.soilstate_init!(ss, np, nc)
    wsb = CLM.WaterStateBulkData{FT}();      CLM.waterstatebulk_init!(wsb, nc, np, nl, ng)
    wdb = CLM.WaterDiagnosticBulkData{FT}(); CLM.waterdiagnosticbulk_init!(wdb, nc, np, nl, ng)
    wfb = CLM.WaterFluxBulkData{FT}();       CLM.waterfluxbulk_init!(wfb, nc, np, nl, ng)
    sa = CLM.SolarAbsorbedData{FT}();        CLM.solarabs_init!(sa, np, nl)
    cs = CLM.CanopyStateData{FT}();          CLM.canopystate_init!(cs, np)
    up = CLM.UrbanParamsData{FT}();          CLM.urbanparams_init!(up, nl; nlevurb = NU)
    return (; nc, np, nl, ng, col, lun, patch, temp, ef, ss, wsb, wdb, wfb, sa, cs, up)
end

# Populate one soil column (c=p=l=1), mirroring test_soil_temperature.jl's setup_soil_column!.
function populate!(B, mask_nolakec, mask_nolakep)
    c = p = l = 1
    col, lun, patch, temp, ef, ss, wsb, wdb, wfb, sa, cs, up =
        B.col, B.lun, B.patch, B.temp, B.ef, B.ss, B.wsb, B.wdb, B.wfb, B.sa, B.cs, B.up
    mask_nolakec[c] = true; mask_nolakep[p] = true
    col.landunit[c] = l; col.itype[c] = CLM.ISTSOIL
    lun.itype[l] = CLM.ISTSOIL; lun.urbpoi[l] = false
    col.snl[c] = 0; col.nbedrock[c] = NG
    patch.column[p] = c; patch.landunit[p] = l; patch.wtcol[p] = 1.0
    for j in 1:NMAX
        jj = j + JOFF
        col.dz[c, jj] = 0.1 * j
        col.z[c, jj] = sum(0.1 * k for k in 1:j) - 0.5 * 0.1 * j
    end
    col.zi[c, JOFF + 1] = 0.0
    for j in 1:NMAX; col.zi[c, j + JOFF + 1] = col.zi[c, j + JOFF] + col.dz[c, j + JOFF]; end
    for j in 1:NMAX; temp.t_soisno_col[c, j + JOFF] = 280.0; end
    temp.t_grnd_col[c] = 280.0; temp.t_h2osfc_col[c] = 280.0; temp.emg_col[c] = 0.97
    for j in 1:NMAX
        jj = j + JOFF
        wsb.ws.h2osoi_liq_col[c, jj] = 5.0; wsb.ws.h2osoi_ice_col[c, jj] = 0.0
    end
    wsb.ws.h2osfc_col[c] = 0.0; wsb.ws.h2osno_no_layers_col[c] = 0.0
    for j in 1:NMAX; wsb.ws.excess_ice_col[c, j] = 0.0; end
    wsb.int_snow_col[c] = 0.0
    wdb.frac_sno_eff_col[c] = 0.0; wdb.frac_h2osfc_col[c] = 0.0; wdb.snow_depth_col[c] = 0.0
    for j in 1:NG
        ss.watsat_col[c, j] = 0.45; ss.bsw_col[c, j] = 5.0; ss.sucsat_col[c, j] = 100.0
        ss.tkmg_col[c, j] = 2.0; ss.tkdry_col[c, j] = 0.25; ss.csol_col[c, j] = 2.0e6
        ss.tksatu_col[c, j] = 1.5; ss.thk_col[c, j] = 0.0
    end
    ef.htvp_col[c] = CLM.HVAP; ef.eflx_bot_col[c] = 0.0
    ef.dlrad_patch[p] = 50.0; ef.cgrnd_patch[p] = 5.0
    ef.eflx_sh_grnd_patch[p] = 10.0; ef.eflx_sh_snow_patch[p] = 10.0
    ef.eflx_sh_soil_patch[p] = 10.0; ef.eflx_sh_h2osfc_patch[p] = 10.0
    ef.dgnetdT_patch[p] = 0.0; ef.eflx_gnet_patch[p] = 0.0
    sa.sabg_patch[p] = 100.0; sa.sabg_soil_patch[p] = 100.0; sa.sabg_snow_patch[p] = 0.0
    sa.sabg_chk_patch[p] = 0.0
    for j in 1:(NS + 1); sa.sabg_lyr_patch[p, j] = 0.0; end
    sa.sabg_lyr_patch[p, NS + 1] = 100.0
    cs.frac_veg_nosno_patch[p] = 0
    wfb.wf.qflx_evap_soi_patch[p] = 0.0; wfb.wf.qflx_ev_snow_patch[p] = 0.0
    wfb.wf.qflx_ev_soil_patch[p] = 0.0; wfb.wf.qflx_ev_h2osfc_patch[p] = 0.0
    wfb.wf.qflx_tran_veg_patch[p] = 0.0
    wfb.qflx_ev_snow_patch[p] = 0.0; wfb.qflx_ev_soil_patch[p] = 0.0; wfb.qflx_ev_h2osfc_patch[p] = 0.0
    wfb.wf.qflx_snomelt_col[c] = 0.0; wfb.wf.qflx_snofrz_col[c] = 0.0
    wfb.wf.qflx_snow_drain_col[c] = 0.0; wfb.wf.qflx_h2osfc_to_ice_col[c] = 0.0
    for j in 1:NS
        wfb.wf.qflx_snofrz_lyr_col[c, j] = 0.0; wfb.wf.qflx_snomelt_lyr_col[c, j] = 0.0
        wdb.bw_col[c, j] = 0.0
    end
    ef.eflx_snomelt_col[c] = 0.0; ef.eflx_snomelt_r_col[c] = 0.0; ef.eflx_snomelt_u_col[c] = 0.0
    ef.eflx_h2osfc_to_snow_col[c] = 0.0; ef.eflx_building_heat_errsoi_col[c] = 0.0
    ef.eflx_urban_ac_col[c] = 0.0; ef.eflx_urban_heat_col[c] = 0.0; ef.eflx_fgr12_col[c] = 0.0
    for j in 1:NG; ef.eflx_fgr_col[c, j] = 0.0; end
    up.nlev_improad[l] = 0; up.t_building_min[l] = 288.0
    wdb.snomelt_accum_col[c] = 0.0
    for j in 1:NMAX; wdb.exice_subs_col[c, j] = 0.0; end
    return nothing
end

function run_soiltemp!(B, urbantv, forc_lwrad, m_c, m_p, m_l, m_uc, dtime)
    CLM.soil_temperature!(B.col, B.lun, B.patch, B.temp, B.ef, B.ss, B.wsb, B.wdb, B.wfb,
        B.sa, B.cs, B.up, urbantv, forc_lwrad, m_c, m_p, m_l, m_uc, 1:B.nc, 1:B.nl, 1:B.np, dtime)
end

function main(backend)
    println("=" ^ 70)
    println("END-TO-END Metal parity for soil_temperature! (whole driver)")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU driver exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    vp = CLM.varpar
    saved = (vp.nlevsno, vp.nlevgrnd, vp.nlevurb, vp.nlevmaxurbgrnd, vp.nlevsoi)
    vp.nlevsno, vp.nlevgrnd, vp.nlevurb, vp.nlevmaxurbgrnd, vp.nlevsoi = NS, NG, NU, NMAX, NSOI
    sv_xi = CLM.varctl.use_excess_ice; CLM.varctl.use_excess_ice = false
    sv_rn = CLM.urban_ctrl.read_namelist; sv_bt = CLM.urban_ctrl.building_temp_method
    CLM.urban_ctrl.read_namelist = true
    CLM.urban_ctrl.building_temp_method = CLM.BUILDING_TEMP_METHOD_SIMPLE
    try
        B = build(FT)
        m_c = falses(B.nc); m_p = falses(B.np); m_l = falses(B.nl); m_uc = falses(B.nc)
        populate!(B, m_c, m_p)
        urbantv = fill(FT(323.15), B.nl)
        forc_lwrad = fill(FT(300.0), B.nc)
        dt = FT(1800)

        # Device snapshot of the populated initial state (adapt copies), BEFORE the CPU run mutates B.
        # NB: the Adapt adaptor must be the device ARRAY TYPE (Metal.MtlArray), not the `to`
        # converter function — adapt(::Function, x) is a silent no-op.
        ad(x) = CLM.Adapt.adapt(Metal.MtlArray, x)
        D = (; nc = B.nc, np = B.np, nl = B.nl,
            col = ad(B.col), lun = ad(B.lun), patch = ad(B.patch), temp = ad(B.temp),
            ef = ad(B.ef), ss = ad(B.ss), wsb = ad(B.wsb), wdb = ad(B.wdb), wfb = ad(B.wfb),
            sa = ad(B.sa), cs = ad(B.cs), up = ad(B.up))
        dmask(m) = to(collect(Bool, m))

        # CPU reference run
        run_soiltemp!(B, urbantv, forc_lwrad, m_c, m_p, m_l, m_uc, dt)
        # Sanity: every adapted struct field must have actually moved to the device. Structs
        # whose fields are pinned to `Vector{FT}`/`Matrix{FT}` (rather than a loose array-type
        # param) are a no-op under adapt and would hand a host Array to a device kernel.
        if !(D.wdb.frac_sno_eff_col isa Metal.MtlArray)
            println("  BLOCKED: WaterDiagnosticBulkData did not move to the device under adapt.")
            println("           Its fields are typed ::Vector{FT}/::Matrix{FT} (pinned), unlike the")
            println("           other 11 state structs (loose array-type params, which adapted fine).")
            println("           → reparametrize WaterDiagnosticBulkData with loose array types +")
            println("             Adapt.@adapt_structure (same treatment WaterStateBulkData already has).")
            println("\n  soil_temperature! itself is fully kernelized + per-kernel Metal/AD-validated;")
            println("  the whole-function device run is gated only on this struct's adaptability.")
            return 2
        end
        # Device run
        run_soiltemp!(D, to(urbantv), to(forc_lwrad), dmask(m_c), dmask(m_p), dmask(m_l), dmask(m_uc), dt)

        checks = [
            ("t_soisno", B.temp.t_soisno_col, D.temp.t_soisno_col),
            ("t_grnd",   B.temp.t_grnd_col,   D.temp.t_grnd_col),
            ("t_h2osfc", B.temp.t_h2osfc_col, D.temp.t_h2osfc_col),
            ("fact",     B.temp.fact_col,     D.temp.fact_col),
            ("h2osoi_liq", B.wsb.ws.h2osoi_liq_col, D.wsb.ws.h2osoi_liq_col),
            ("h2osoi_ice", B.wsb.ws.h2osoi_ice_col, D.wsb.ws.h2osoi_ice_col),
            ("eflx_fgr12", B.ef.eflx_fgr12_col, D.ef.eflx_fgr12_col),
            ("eflx_fgr",   B.ef.eflx_fgr_col,   D.ef.eflx_fgr_col),
            ("xmf",        B.temp.xmf_col,       D.temp.xmf_col),
        ]
        nfail = 0
        for (nm, a, b) in checks
            d = maxabsdiff(a, b)
            ok = d < 5f-3
            @printf("  [%s] %-12s  max|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
            ok || (nfail += 1)
        end
        println()
        println(nfail == 0 ? "  WHOLE soil_temperature! MATCHES CPU ON $name ($FT) ✓" :
                             "  DIVERGENCE — investigate.")
        return nfail == 0 ? 0 : 1
    finally
        vp.nlevsno, vp.nlevgrnd, vp.nlevurb, vp.nlevmaxurbgrnd, vp.nlevsoi = saved
        CLM.varctl.use_excess_ice = sv_xi
        CLM.urban_ctrl.read_namelist = sv_rn
        CLM.urban_ctrl.building_temp_method = sv_bt
    end
end

const BACKEND = detect_backend()
exit(main(BACKEND))
