# =============================================================================
# PRE-soil_water! surface-hydrology reverse phases (HydrologyNoDrainage), through the
# production CLM.compositional_reverse! engine. Extends the reverse coverage to the
# surface water partitioning that feeds soil_water!'s top boundary:
#   [K] saturated_excess_runoff!  (fsat = wtfact·exp(-0.5·fff·zwt) → saturated runoff)
#   [L] infiltration_excess_runoff! (fsat → qinmax, Hortonian excess)
#   [M] infiltration!              (qflx_infl = qflx_in_soil_limited − qflx_h2osfc_drain)
#   [N] saturated_excess_runoff! → infiltration_excess_runoff! sub-chain (the fsat link)
# Wrappers are productionized in src/driver/driver_reverse.jl (satexcess/inflexcess/
# infil _rev_phase!). julia +1.10 --project=/tmp/clm_jl10_<id> this_script.jl
# =============================================================================
using CLM, Enzyme, Printf
const C = CLM; const FT = Float64
const FSURDAT  = get(ENV, "CLM_FSURDAT",  "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc")
const PARAMFILE = get(ENV, "CLM_PARAMFILE", "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc")

function setup_forcing!(a2l, T, ng)
    for g in 1:ng
        a2l.forc_t_not_downscaled_grc[g]=T; a2l.forc_pbot_not_downscaled_grc[g]=85000.0
        a2l.forc_th_not_downscaled_grc[g]=T*(100000.0/85000.0)^(C.RAIR/C.CPAIR)
        a2l.forc_rho_not_downscaled_grc[g]=85000.0/(C.RAIR*T); a2l.forc_lwrad_not_downscaled_grc[g]=300.0
        a2l.forc_vp_grc[g]=800.0; a2l.forc_hgt_grc[g]=30.0; a2l.forc_topo_grc[g]=0.0; a2l.forc_wind_grc[g]=3.0
        a2l.forc_hgt_u_grc[g]=30.0; a2l.forc_hgt_t_grc[g]=30.0; a2l.forc_hgt_q_grc[g]=30.0
        for b in 1:C.NUMRAD; a2l.forc_solad_not_downscaled_grc[g,b]=200.0; a2l.forc_solai_grc[g,b]=80.0; end
        a2l.forc_solar_not_downscaled_grc[g]=560.0; a2l.forc_rain_not_downscaled_grc[g]=0.001  # wet → live infil
        a2l.forc_snow_not_downscaled_grc[g]=0.0
    end
end
function build_real()
    (inst, bounds, filt, tm) = C.clm_initialize!(; fsurdat=FSURDAT, paramfile=PARAMFILE)
    setup_forcing!(inst.atm2lnd, 285.0, bounds.endg)
    C.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
    config = C.CLMDriverConfig(); fia = C.clump_filter_inactive_and_active
    (declin, _) = C.compute_orbital(120.0); nextsw = 120.0 + 1800.0/C.SECSPDAY
    for n in 1:3
        C.clm_drv!(config, inst, filt, fia, bounds, true, nextsw, declin, declin,
            C.ORB_OBLIQR_DEFAULT, false, false, "", false; nstep=n, is_first_step=(n==1),
            is_beg_curr_day=(n==1), dtime=1800.0, mon=1, day=1, photosyns=inst.photosyns)
    end
    return inst, bounds, filt, config
end

# ---- phase wrappers (operate on inst i; Bool kwargs resolved inside) ----
saux(bounds, filt) = (; hydrologyc=filt.hydrologyc, bc_col=bounds.begc:bounds.endc)
function satx_phase!(i, aux)
    C.saturated_excess_runoff!(i.sat_excess_runoff, aux.hydrologyc, aux.bc_col,
        i.column, i.landunit, i.soilhydrology, i.soilstate, i.water.waterfluxbulk_inst)
    return nothing
end
function inflx_phase!(i, aux)
    C.infiltration_excess_runoff!(i.infilt_excess_runoff, i.soilhydrology, i.soilstate,
        i.sat_excess_runoff.fsat_col, i.water.waterfluxbulk_inst, i.water.waterdiagnosticbulk_inst,
        aux.hydrologyc, aux.bc_col)
    return nothing
end
function infil_phase!(i, aux)
    C.infiltration!(i.water.waterfluxbulk_inst, aux.hydrologyc, aux.bc_col)
    return nothing
end

richardson_fd(L_pert) = begin
    cfd(h) = (L_pert(h) - L_pert(-h)) / (2h)
    (4*cfd(5e-3) - cfd(1e-2)) / 3
end
function verdict(tag, g_fd, g_rev)
    rl = abs(g_rev - g_fd) / max(abs(g_fd), 1e-12)
    @printf("[%s] FD = % .8e  rev = % .8e  rel = %.3e  %s\n", tag, g_fd, g_rev, rl,
            rl < 1e-4 ? "PASS ✓" : "FAIL ✗")
    return rl
end

inst, bounds, filt, config = build_real()
aux = saux(bounds, filt)
hc = filt.hydrologyc; c0 = [c for c in bounds.begc:bounds.endc if hc[c]][1]
j0 = C.varpar.nlevsno + 4

# [K] saturated_excess_runoff!: perturb zwt → fsat, L = sum(fsat²)
rK = try
    println("\n", "#"^70, "\n# [K] saturated_excess_runoff! REVERSE vs FD\n", "#"^70)
    Lf(s) = sum(abs2, s.sat_excess_runoff.fsat_col)
    Lp(δ) = (s=deepcopy(inst); s.soilhydrology.zwt_col[c0]+=δ; satx_phase!(s, aux); Lf(s))
    gfd = richardson_fd(Lp)
    db = C.compositional_reverse!(Any[(satx_phase!, (aux,))], deepcopy(inst),
        (db,b)->(db.sat_excess_runoff.fsat_col .= 2 .* b.sat_excess_runoff.fsat_col))
    verdict("K", gfd, db.soilhydrology.zwt_col[c0])
catch e; @printf("[K] ERRORED: %s\n", sprint(showerror,e)); NaN end

# [L] infiltration_excess_runoff!: perturb fsat_col → qinmax, L = sum(qinmax²)
rL = try
    println("\n", "#"^70, "\n# [L] infiltration_excess_runoff! REVERSE vs FD\n", "#"^70)
    # need fsat populated first
    let s=deepcopy(inst); satx_phase!(s, aux); inst.sat_excess_runoff.fsat_col .= s.sat_excess_runoff.fsat_col; end
    Lf(s) = sum(abs2, s.infilt_excess_runoff.qinmax_col)
    Lp(δ) = (s=deepcopy(inst); s.sat_excess_runoff.fsat_col[c0]+=δ; inflx_phase!(s, aux); Lf(s))
    gfd = richardson_fd(Lp)
    db = C.compositional_reverse!(Any[(inflx_phase!, (aux,))], deepcopy(inst),
        (db,b)->(db.infilt_excess_runoff.qinmax_col .= 2 .* b.infilt_excess_runoff.qinmax_col))
    verdict("L", gfd, db.sat_excess_runoff.fsat_col[c0])
catch e; @printf("[L] ERRORED: %s\n", sprint(showerror,e)); NaN end

# [M] infiltration!: perturb qflx_in_soil_limited → qflx_infl, L = sum(qflx_infl²)
rM = try
    println("\n", "#"^70, "\n# [M] infiltration! REVERSE vs FD\n", "#"^70)
    Lf(s) = sum(abs2, s.water.waterfluxbulk_inst.wf.qflx_infl_col)
    Lp(δ) = (s=deepcopy(inst); s.water.waterfluxbulk_inst.qflx_in_soil_limited_col[c0]+=δ; infil_phase!(s, aux); Lf(s))
    gfd = richardson_fd(Lp)
    db = C.compositional_reverse!(Any[(infil_phase!, (aux,))], deepcopy(inst),
        (db,b)->(db.water.waterfluxbulk_inst.wf.qflx_infl_col .= 2 .* b.water.waterfluxbulk_inst.wf.qflx_infl_col))
    verdict("M", gfd, db.water.waterfluxbulk_inst.qflx_in_soil_limited_col[c0])
catch e; @printf("[M] ERRORED: %s\n", sprint(showerror,e)); NaN end

# [N] sat_excess → infl_excess SUB-CHAIN: perturb zwt → fsat → qinmax, L = sum(qinmax²)
rN = try
    println("\n", "#"^70, "\n# [N] saturated_excess → infiltration_excess SUB-CHAIN vs FD\n", "#"^70)
    phases = Any[(satx_phase!, (aux,)), (inflx_phase!, (aux,))]
    Lf(s) = sum(abs2, s.infilt_excess_runoff.qinmax_col)
    Lp(δ) = (s=deepcopy(inst); s.soilhydrology.zwt_col[c0]+=δ; for (f,ca) in phases; f(s,ca...); end; Lf(s))
    gfd = richardson_fd(Lp)
    db = C.compositional_reverse!(phases, deepcopy(inst),
        (db,b)->(db.infilt_excess_runoff.qinmax_col .= 2 .* b.infilt_excess_runoff.qinmax_col))
    verdict("N", gfd, db.soilhydrology.zwt_col[c0])
catch e; @printf("[N] ERRORED: %s\n", sprint(showerror,e)); NaN end

println("\n", "="^70)
@printf("SUMMARY  [K] sat_excess=%.3e  [L] infl_excess=%.3e  [M] infiltration=%.3e  [N] subchain=%.3e\n", rK, rL, rM, rN)
println("="^70)
