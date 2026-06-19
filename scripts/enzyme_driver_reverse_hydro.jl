# =============================================================================
# EXTENSION of scripts/enzyme_driver_reverse_full.jl by more clm_drv! HydrologyNoDrainage
# phases: water_table! ([G]/[H]) and hydrology_no_drainage! ([I]/[J], the namesake
# end-of-step diagnostic). All through the production CLM.compositional_reverse! engine
# (src/biogeophys/canopy_fluxes_reverse.jl); the phase wrappers are productionized in
# src/driver/driver_reverse.jl.
#
#   [I] hydrology_no_drainage! STANDALONE reverse (machine precision 1.1e-11): perturb
#       h2osoi_liq, L = sum(h2osoi_vol²) (the diagnostic recomputes vol from liq/ice).
#   [J] canopy+soil_temp+soil_water+water_table+hydrology_no_drainage FIVE-module chain
#       (4.1e-7, FD-bracketed): perturb t_grnd → gradient flows across five module
#       boundaries to L = sum(h2osoi_vol²).
#
#   [G] water_table! STANDALONE reverse on the real Bow cold-start inst. Perturb
#       qcharge_col (the live aquifer-recharge INPUT the phase integrates: zwt -=
#       qcharge·dtime/…) and FD-validate dL/d(qcharge), L = sum(abs2, zwt_col). This
#       is the guaranteed-live test that Enzyme reverse-differentiates water_table!
#       (machine precision: rel ≈ 1.3e-16).
#
#   [H] canopy + soil_temperature! + soil_water! + water_table! FOUR-module chain in
#       ONE compositional_reverse! call — extends [F] by one module. L = sum(abs2,
#       h2osoi_liq_col) (the proven-live [F] output) drives the reverse END-TO-END back
#       through water_table!, soil_water!, soil_temp!, canopy!, ACROSS four module
#       boundaries; FD-validates dL/d(t_grnd) (~1.4e-3, FD-bracketed). NB: seeding zwt
#       gives a true-zero here — at this cold-start state the water-table outputs are a
#       flat branch w.r.t. the upstream chain (correct AD); water_table!'s own reverse is
#       validated live + exactly by [G].
#
#   julia +1.10 --project=/tmp/clm_jl10_<id> scripts/enzyme_driver_reverse_hydro.jl
# =============================================================================
using CLM, Enzyme, Printf
const C = CLM
const FT = Float64

const FSURDAT = get(ENV, "CLM_FSURDAT",
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc")
const PARAMFILE = get(ENV, "CLM_PARAMFILE",
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc")

# ---- real Bow cold-start builder (same as enzyme_driver_reverse_full.jl) ----
function setup_forcing!(a2l, T, ng)
    for g in 1:ng
        a2l.forc_t_not_downscaled_grc[g]=T; a2l.forc_pbot_not_downscaled_grc[g]=85000.0
        a2l.forc_th_not_downscaled_grc[g]=T*(100000.0/85000.0)^(C.RAIR/C.CPAIR)
        a2l.forc_rho_not_downscaled_grc[g]=85000.0/(C.RAIR*T); a2l.forc_lwrad_not_downscaled_grc[g]=300.0
        a2l.forc_vp_grc[g]=800.0; a2l.forc_hgt_grc[g]=30.0; a2l.forc_topo_grc[g]=0.0; a2l.forc_wind_grc[g]=3.0
        a2l.forc_hgt_u_grc[g]=30.0; a2l.forc_hgt_t_grc[g]=30.0; a2l.forc_hgt_q_grc[g]=30.0
        for b in 1:C.NUMRAD; a2l.forc_solad_not_downscaled_grc[g,b]=200.0; a2l.forc_solai_grc[g,b]=80.0; end
        a2l.forc_solar_not_downscaled_grc[g]=560.0; a2l.forc_rain_not_downscaled_grc[g]=0.0001
        a2l.forc_snow_not_downscaled_grc[g]=0.0
    end
end

function build_real()
    (inst, bounds, filt, tm) = C.clm_initialize!(; fsurdat=FSURDAT, paramfile=PARAMFILE)
    setup_forcing!(inst.atm2lnd, 285.0, bounds.endg)
    C.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
    config = C.CLMDriverConfig(); filt_ia = C.clump_filter_inactive_and_active
    (declin, eccf) = C.compute_orbital(120.0); nextsw_cday = 120.0 + 1800.0/C.SECSPDAY
    runstep!(n; first=false) = C.clm_drv!(config, inst, filt, filt_ia, bounds,
        true, nextsw_cday, declin, declin, C.ORB_OBLIQR_DEFAULT, false, false, "", false;
        nstep=n, is_first_step=first, is_beg_curr_day=first, dtime=1800.0, mon=1, day=1,
        photosyns=inst.photosyns)
    for n in 1:3; runstep!(n; first=(n==1)); end
    return inst, bounds, filt, config
end

# ---- soil_temperature! phase (same as base) ----
soil_aux(bounds, filt) = (;
    urbantv = fill(323.15, bounds.endl), dtime = 1800.0,
    bc_col = bounds.begc:bounds.endc, bc_lun = bounds.begl:bounds.endl,
    bc_patch = bounds.begp:bounds.endp,
    nolakec = filt.nolakec, nolakep = filt.nolakep, urbanl = filt.urbanl, urbanc = filt.urbanc)
function soil_phase!(b, aux)
    C.soil_temperature!(b.column, b.landunit, b.patch, b.temperature, b.energyflux,
        b.soilstate, b.water.waterstatebulk_inst, b.water.waterdiagnosticbulk_inst,
        b.water.waterfluxbulk_inst, b.solarabs, b.canopystate, b.urbanparams,
        aux.urbantv, b.atm2lnd.forc_lwrad_downscaled_col, aux.nolakec, aux.nolakep,
        aux.urbanl, aux.urbanc, aux.bc_col, aux.bc_lun, aux.bc_patch, aux.dtime)
    return nothing
end

# ---- soil_water! phase (same as base) ----
function water_aux(filt, config)
    swm_cfg = config.use_aquifer_layer ? C.SoilWaterMovementConfig() :
        C.SoilWaterMovementConfig(soilwater_movement_method=C.MOISTURE_FORM,
                                  lower_boundary_condition=C.BC_ZERO_FLUX)
    return (; hydrologyc = filt.hydrologyc, urbanc = filt.urbanc,
              swrc = C.SoilWaterRetentionCurveClappHornberg1978(), cfg = swm_cfg, dtime = 1800.0)
end
function water_phase!(b, aux)
    C.soil_water!(b.column, aux.hydrologyc, aux.urbanc,
        b.soilhydrology, b.soilstate, b.water.waterfluxbulk_inst, b.water.waterstatebulk_inst,
        b.temperature, b.canopystate, b.energyflux, aux.swrc, aux.cfg, aux.dtime)
    return nothing
end

# ---- NEW: water_table! phase. The single Bool kwarg (recompute_frost_table) is
#      resolved INSIDE the wrapper (primal-time) — not an Enzyme arg — so Enzyme
#      differentiates the wrapper body without a kwarg-NamedTuple thunk. ----
wt_aux(bounds, filt, config) = (; hydrologyc = filt.hydrologyc,
    bc_col = bounds.begc:bounds.endc, nlevsoi = C.varpar.nlevsoi, dtime = 1800.0,
    recompute_frost_table = config.use_aquifer_layer)
function wt_phase!(b, aux)
    C.water_table!(b.soilhydrology, b.soilstate,
        b.temperature.t_soisno_col, b.water.waterstatebulk_inst.ws, b.water.waterfluxbulk_inst,
        b.column.dz, b.column.z, b.column.zi,
        aux.hydrologyc, aux.bc_col, aux.nlevsoi, aux.dtime;
        recompute_frost_table = aux.recompute_frost_table)
    return nothing
end

# ---- canopy chain helpers (same as base [F]) ----
chain_bundle(s) = (; inst = s, scratch = C.cf_rev_scratch(FT, length(s.patch.column)))
chain_cf_view(b) = C.cf_rev_bundle(b.inst.canopystate, b.inst.energyflux, b.inst.frictionvel,
    b.inst.temperature, b.inst.solarabs, b.inst.soilstate, b.inst.water.waterfluxbulk_inst,
    b.inst.water.waterstatebulk_inst, b.inst.water.waterdiagnosticbulk_inst, b.inst.photosyns,
    b.scratch)
function chain_canopy_phase!(b, aux, Ncanopy)
    cv = chain_cf_view(b)
    for (f, cargs) in C.cf_rev_phases(aux, Ncanopy); f(cv, cargs...); end
    return nothing
end
chain_soil_phase!(b, saux)  = soil_phase!(b.inst, saux)
chain_water_phase!(b, waux) = water_phase!(b.inst, waux)
chain_wt_phase!(b, wtaux)   = wt_phase!(b.inst, wtaux)

# hydrology_no_drainage! phase (end-of-step diagnostics: recompute h2osoi_vol etc.)
hnd_aux(bounds, filt) = (; nolakec=filt.nolakec, hydrologyc=filt.hydrologyc, urbanc=filt.urbanc,
    snowc=filt.snowc, nosnowc=filt.nosnowc, bc_col=bounds.begc:bounds.endc, dtime=1800.0,
    nlevsno=C.varpar.nlevsno, nlevsoi=C.varpar.nlevsoi, nlevgrnd=C.varpar.nlevgrnd, nlevurb=C.varpar.nlevurb)
function hnd_phase!(b, aux)
    C.hydrology_no_drainage!(b.temperature, b.soilstate, b.water.waterstatebulk_inst,
        b.water.waterdiagnosticbulk_inst, b.column, b.landunit,
        aux.nolakec, aux.hydrologyc, aux.urbanc, aux.snowc, aux.nosnowc,
        aux.bc_col, aux.dtime, aux.nlevsno, aux.nlevsoi, aux.nlevgrnd, aux.nlevurb)
    return nothing
end
chain_hnd_phase!(b, hndaux) = hnd_phase!(b.inst, hndaux)

function canopy_aux(inst, bounds, filt)
    NP = bounds.endp; ev = filt.exposedvegp
    fp = Int[p for p in bounds.begp:bounds.endp if ev[p]]; fn = length(fp)
    a2l = inst.atm2lnd
    forc_q_col = fill!(similar(a2l.forc_pbot_downscaled_col), zero(FT))
    C.compute_forc_q!(forc_q_col, inst.column.gridcell, a2l.forc_vp_grc, a2l.forc_pbot_downscaled_col)
    forc = (; lwrad = a2l.forc_lwrad_downscaled_col, q = forc_q_col,
              pbot = a2l.forc_pbot_downscaled_col, th = a2l.forc_th_downscaled_col,
              rho = a2l.forc_rho_downscaled_col, t = a2l.forc_t_downscaled_col,
              u_grc = a2l.forc_u_grc, v_grc = a2l.forc_v_grc, pco2 = a2l.forc_pco2_grc, po2 = a2l.forc_po2_grc,
              hgt_t = a2l.forc_hgt_t_grc, hgt_u = a2l.forc_hgt_u_grc, hgt_q = a2l.forc_hgt_q_grc,
              dayl = inst.gridcell.dayl, max_dayl = inst.gridcell.max_dayl,
              downreg = fill(FT(1.0), NP), leafn = fill(FT(1.0), NP))
    MP = C.MXPFT + 1
    pft = (; dleaf = fill(FT(0.04), MP), z0v_Cr = fill(FT(0.35), MP), z0v_Cs = fill(FT(0.003), MP),
             z0v_c = fill(FT(0.25), MP), z0v_cw = fill(FT(2.0), MP), z0v_LAImax = fill(FT(8.0), MP),
             grnd_ch4 = fill(FT(0.0), NP))
    psn = (; c3psn = fill(FT(1.0), MP), leafcn = fill(FT(25.0), MP), flnr = fill(FT(0.1), MP),
             fnitr = fill(FT(1.0), MP), slatop = fill(FT(0.01), MP), mbbopt = fill(FT(9.0), MP),
             medlynintercept = fill(FT(100.0), MP), medlynslope = fill(FT(6.0), MP),
             nrad = fill(1, NP), tlai_z = fill(FT(1.0), NP, C.NLEVCAN),
             parsun_z = fill(FT(0.0), NP, C.NLEVCAN), parsha_z = fill(FT(0.0), NP, C.NLEVCAN),
             laisun_z = fill(FT(0.0), NP, C.NLEVCAN), laisha_z = fill(FT(0.0), NP, C.NLEVCAN),
             vcmaxcint_sun = fill(FT(1.0), NP), vcmaxcint_sha = fill(FT(0.6), NP),
             o3coefv = fill(FT(1.0), NP), o3coefg = fill(FT(1.0), NP), t10 = inst.temperature.t_a10_patch)
    return (; patch = inst.patch, col = inst.column, grid = inst.gridcell, forc = forc, pft = pft, psn = psn,
        filterp = fp, fn = fn, active = Bool[ev[p] for p in 1:NP], mask = ev, ivt = inst.patch.itype .+ 1,
        forc_pbot_patch = FT[a2l.forc_pbot_downscaled_col[inst.patch.column[p]] for p in 1:NP],
        soilevap_beta = C.do_soilevap_beta(), soil_resis_sl14 = C.do_soil_resistance_sl14(),
        nlevsno = C.varpar.nlevsno, dtime = FT(1800.0), use_psn = false)
end

# shared FD-bracket verdict (same as base report_chain)
function report_chain(tag, label, L_pert, g_rev, seedname, idx)
    @printf("perturb %s[%s]; L = sum(h2osoi_liq^2) AFTER %s\n", seedname, string(idx), label)
    hs = (2e-2, 1e-2, 5e-3, 2.5e-3); cfd(h) = (L_pert(h) - L_pert(-h)) / (2h)
    fds = [cfd(h) for h in hs]
    for (h, f) in zip(hs, fds); @printf("  FD(h=%.2e) = % .8e\n", h, f); end
    g_rich = (4 * fds[end] - fds[end-1]) / 3
    lo, hi = minimum(fds), maximum(fds)
    bracketed = (lo - abs(lo)*1e-9) <= g_rev <= (hi + abs(hi)*1e-9)
    relerr = abs(g_rev - g_rich) / max(abs(g_rich), 1e-10)
    @printf("  Richardson = % .8e\n  rev = % .8e\n", g_rich, g_rev)
    pass = relerr < 5e-5 || bracketed
    @printf("[%s] rel err vs Richardson = %.3e  bracketed=%s  %s\n",
        tag, relerr, string(bracketed), pass ? "PASS ✓" : "FAIL ✗")
    return relerr
end

# =============================================================================
# [G] water_table! STANDALONE reverse — perturb qcharge_col (live recharge input).
# =============================================================================
function section_G(inst, bounds, filt, config)
    println("\n", "#"^70, "\n# [G] water_table! REVERSE (real cold-start inst) vs FD\n", "#"^70)
    wtaux = wt_aux(bounds, filt, config)
    hc = filt.hydrologyc; cs = [c for c in bounds.begc:bounds.endc if hc[c]]
    @assert !isempty(cs) "no hydrology columns"; c0 = cs[1]

    # zwt is the live state water_table! evolves here: zwt -= qcharge·dtime/… (the
    # aquifer-recharge integration). wa stays at baseline in this regime (a flat branch),
    # so seed zwt and perturb qcharge — d(zwt)/d(qcharge) is smooth and nonzero.
    L(s) = sum(abs2, s.soilhydrology.zwt_col)
    let s = deepcopy(inst); wt_phase!(s, wtaux)
        @printf("standalone primal: zwt[%d]=%.6f qcharge[%d]=%.3e (finite: %s)\n", c0,
            s.soilhydrology.zwt_col[c0], c0, inst.soilhydrology.qcharge_col[c0],
            string(isfinite(s.soilhydrology.zwt_col[c0])))
    end
    function L_pert(δ)
        s = deepcopy(inst); s.soilhydrology.qcharge_col[c0] += δ
        wt_phase!(s, wtaux); return L(s)
    end
    cfd(h) = (L_pert(h) - L_pert(-h)) / (2h)
    g_fd = (4*cfd(5e-3) - cfd(1e-2)) / 3
    seed!(db, b) = (db.soilhydrology.zwt_col .= 2 .* b.soilhydrology.zwt_col)
    db = C.compositional_reverse!(Any[(wt_phase!, (wtaux,))], deepcopy(inst), seed!)
    g_rev = db.soilhydrology.qcharge_col[c0]
    rl = abs(g_rev - g_fd) / max(abs(g_fd), 1e-10)
    @printf("[G] perturb qcharge[%d]; L=sum(zwt^2); FD = % .8e  rev = % .8e  rel = %.3e  %s\n",
        c0, g_fd, g_rev, rl, rl < 1e-5 ? "PASS ✓" : "FAIL ✗")
    return rl
end

# =============================================================================
# [H] canopy + soil_temp + soil_water + water_table FOUR-module chain.
# =============================================================================
function section_H(inst, bounds, filt, config)
    println("\n", "#"^70, "\n# [H] canopy+soil_temp+soil_water+water_table FOUR-module CHAIN vs FD\n", "#"^70)
    aux  = canopy_aux(inst, bounds, filt); saux = soil_aux(bounds, filt)
    waux = water_aux(filt, config); wtaux = wt_aux(bounds, filt, config)
    Ncanopy = parse(Int, get(ENV, "CLM_NCANOPY", "12"))
    hc = filt.hydrologyc; cs = [c for c in bounds.begc:bounds.endc if hc[c]]; c0 = cs[1]

    chain_phases = Any[(chain_canopy_phase!, (aux, Ncanopy)), (chain_soil_phase!, (saux,)),
                       (chain_water_phase!, (waux,)), (chain_wt_phase!, (wtaux,))]
    # Seed h2osoi_liq (the proven-live chain output, as in [F]): this drives the reverse
    # END-TO-END through all FOUR phases — the gradient flows back through water_table!'s
    # reverse, then soil_water!, soil_temp!, canopy! — demonstrating the 4-module chain is
    # differentiable. (Seeding zwt instead gives a TRUE-ZERO here: at this cold-start state
    # the water-table outputs zwt/qcharge are a flat branch w.r.t. the upstream — correct AD,
    # but uninformative. water_table!'s own reverse is validated live + exactly by [G].)
    let bs = chain_bundle(deepcopy(inst))
        for (f, ca) in chain_phases; f(bs, ca...); end
        @printf("chained primal: t_veg ok=%s  h2osoi_liq[%d,%d]=%.6f (finite: %s)\n",
            string(all(isfinite, bs.inst.temperature.t_veg_patch)), c0, C.varpar.nlevsno+4,
            bs.inst.water.waterstatebulk_inst.ws.h2osoi_liq_col[c0, C.varpar.nlevsno+4],
            string(isfinite(bs.inst.water.waterstatebulk_inst.ws.h2osoi_liq_col[c0, C.varpar.nlevsno+4])))
    end
    Lchain(s) = sum(abs2, s.water.waterstatebulk_inst.ws.h2osoi_liq_col)
    function L_pert(δ)
        s = deepcopy(inst); s.temperature.t_grnd_col[c0] += δ; bs = chain_bundle(s)
        for (f, ca) in chain_phases; f(bs, ca...); end
        return Lchain(bs.inst)
    end
    seed_chain!(db, b) = (db.inst.water.waterstatebulk_inst.ws.h2osoi_liq_col .=
                          2 .* b.inst.water.waterstatebulk_inst.ws.h2osoi_liq_col)
    db = C.compositional_reverse!(chain_phases, chain_bundle(deepcopy(inst)), seed_chain!)
    g_rev = db.inst.temperature.t_grnd_col[c0]
    return report_chain("H", "canopy→soiltemp→soilwater→watertable", L_pert, g_rev, "t_grnd", c0)
end

# =============================================================================
# [I] hydrology_no_drainage! STANDALONE reverse — perturb h2osoi_liq → h2osoi_vol.
# =============================================================================
function section_I(inst, bounds, filt, config)
    println("\n", "#"^70, "\n# [I] hydrology_no_drainage! REVERSE (real cold-start inst) vs FD\n", "#"^70)
    hndaux = hnd_aux(bounds, filt)
    hc = filt.hydrologyc; c0 = [c for c in bounds.begc:bounds.endc if hc[c]][1]; j0 = C.varpar.nlevsno + 4
    vol(s) = s.water.waterstatebulk_inst.ws.h2osoi_vol_col
    L(s) = sum(abs2, vol(s))
    let s = deepcopy(inst); hnd_phase!(s, hndaux)
        @printf("standalone primal: h2osoi_vol[%d,4]=%.6f (finite: %s)\n", c0, vol(s)[c0,4],
            string(isfinite(vol(s)[c0,4])))
    end
    function L_pert(δ)
        s = deepcopy(inst); s.water.waterstatebulk_inst.ws.h2osoi_liq_col[c0,j0] += δ
        hnd_phase!(s, hndaux); return L(s)
    end
    cfd(h) = (L_pert(h) - L_pert(-h)) / (2h)
    g_fd = (4*cfd(5e-3) - cfd(1e-2)) / 3
    seed!(db, b) = (db.water.waterstatebulk_inst.ws.h2osoi_vol_col .=
                    2 .* b.water.waterstatebulk_inst.ws.h2osoi_vol_col)
    db = C.compositional_reverse!(Any[(hnd_phase!, (hndaux,))], deepcopy(inst), seed!)
    g_rev = db.water.waterstatebulk_inst.ws.h2osoi_liq_col[c0,j0]
    rl = abs(g_rev - g_fd) / max(abs(g_fd), 1e-10)
    @printf("[I] perturb h2osoi_liq[%d,%d]; L=sum(vol^2); FD = % .8e  rev = % .8e  rel = %.3e  %s\n",
        c0, j0, g_fd, g_rev, rl, rl < 1e-5 ? "PASS ✓" : "FAIL ✗")
    return rl
end

# =============================================================================
# [J] canopy + soil_temp + soil_water + water_table + hydrology_no_drainage
#     FIVE-module chain — extends [H] by the namesake diagnostic phase.
# =============================================================================
function section_J(inst, bounds, filt, config)
    println("\n", "#"^70, "\n# [J] FIVE-module CHAIN (+ hydrology_no_drainage!) vs FD\n", "#"^70)
    aux  = canopy_aux(inst, bounds, filt); saux = soil_aux(bounds, filt)
    waux = water_aux(filt, config); wtaux = wt_aux(bounds, filt, config); hndaux = hnd_aux(bounds, filt)
    Ncanopy = parse(Int, get(ENV, "CLM_NCANOPY", "12"))
    hc = filt.hydrologyc; c0 = [c for c in bounds.begc:bounds.endc if hc[c]][1]
    chain_phases = Any[(chain_canopy_phase!, (aux, Ncanopy)), (chain_soil_phase!, (saux,)),
                       (chain_water_phase!, (waux,)), (chain_wt_phase!, (wtaux,)), (chain_hnd_phase!, (hndaux,))]
    let bs = chain_bundle(deepcopy(inst))
        for (f, ca) in chain_phases; f(bs, ca...); end
        @printf("chained primal: t_veg ok=%s  h2osoi_vol[%d,4]=%.6f (finite: %s)\n",
            string(all(isfinite, bs.inst.temperature.t_veg_patch)), c0,
            bs.inst.water.waterstatebulk_inst.ws.h2osoi_vol_col[c0,4],
            string(isfinite(bs.inst.water.waterstatebulk_inst.ws.h2osoi_vol_col[c0,4])))
    end
    Lchain(s) = sum(abs2, s.water.waterstatebulk_inst.ws.h2osoi_vol_col)
    function L_pert(δ)
        s = deepcopy(inst); s.temperature.t_grnd_col[c0] += δ; bs = chain_bundle(s)
        for (f, ca) in chain_phases; f(bs, ca...); end
        return Lchain(bs.inst)
    end
    seed_chain!(db, b) = (db.inst.water.waterstatebulk_inst.ws.h2osoi_vol_col .=
                          2 .* b.inst.water.waterstatebulk_inst.ws.h2osoi_vol_col)
    db = C.compositional_reverse!(chain_phases, chain_bundle(deepcopy(inst)), seed_chain!)
    g_rev = db.inst.temperature.t_grnd_col[c0]
    return report_chain("J", "canopy→…→watertable→hydnodrain", L_pert, g_rev, "t_grnd", c0)
end

# ---- driver ----
inst, bounds, filt, config = build_real()
rG = try section_G(inst, bounds, filt, config) catch e; @printf("\n[G] ERRORED: %s\n", sprint(showerror, e)); NaN end
rH = try section_H(inst, bounds, filt, config) catch e; @printf("\n[H] ERRORED: %s\n", sprint(showerror, e)); NaN end
rI = try section_I(inst, bounds, filt, config) catch e; @printf("\n[I] ERRORED: %s\n", sprint(showerror, e)); NaN end
rJ = try section_J(inst, bounds, filt, config) catch e; @printf("\n[J] ERRORED: %s\n", sprint(showerror, e)); NaN end
println("\n", "="^70)
@printf("SUMMARY  [G] water_table standalone           = %.3e\n", rG)
@printf("         [H] canopy+soiltemp+soilwater+wt      = %.3e\n", rH)
@printf("         [I] hydrology_no_drainage standalone  = %.3e\n", rI)
@printf("         [J] FIVE-module chain (+hydnodrain)   = %.3e\n", rJ)
println("="^70)
