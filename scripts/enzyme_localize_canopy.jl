# Localization probe: WHICH sub-call of canopy_fluxes_core! triggers the
# Enzyme sret_ty / propagate_returned! bug?
#
# canopy_fluxes_core! fails to reverse-differentiate standalone (see
# enzyme_phase_probe_canopy.jl). soil_temperature! — which ALSO uses KA kernels —
# differentiates cleanly (enzyme_phase_probe.jl). So the wall is NOT a blanket
# "Enzyme can't do KA kernels"; one specific sub-function inside canopy_fluxes_core!
# triggers it. This probe differentiates the big leaf calls in ISOLATION to localize
# the culprit. The sret bug is a COMPILE-TIME type-analysis failure, so crude-but-
# finite scratch inputs are sufficient to reproduce (or clear) it.
#
#   julia +1.10 --project=/tmp/clm_jl10_env scripts/enzyme_localize_canopy.jl
#
using CLM, Enzyme, Printf
const FSURDAT = get(ENV, "CLM_FSURDAT",
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc")
const PARAMFILE = get(ENV, "CLM_PARAMFILE",
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc")

(inst, bounds, filt, tm) = CLM.clm_initialize!(; fsurdat=FSURDAT, paramfile=PARAMFILE)
ng = bounds.endg; nc = bounds.endc; np = bounds.endp
config = CLM.CLMDriverConfig()
filt_ia = CLM.clump_filter_inactive_and_active
(declin, eccf) = CLM.compute_orbital(120.0)
nextsw_cday = 120.0 + 1800.0 / CLM.SECSPDAY
T0 = 285.0
a2l = inst.atm2lnd
for g in 1:ng
    a2l.forc_t_not_downscaled_grc[g] = T0
    a2l.forc_pbot_not_downscaled_grc[g] = 85000.0
    a2l.forc_th_not_downscaled_grc[g] = T0 * (100000.0/85000.0)^(CLM.RAIR/CLM.CPAIR)
    a2l.forc_rho_not_downscaled_grc[g] = 85000.0/(CLM.RAIR*T0)
    a2l.forc_lwrad_not_downscaled_grc[g] = 300.0
    a2l.forc_vp_grc[g] = 800.0; a2l.forc_hgt_grc[g] = 30.0; a2l.forc_topo_grc[g] = 0.0
    a2l.forc_wind_grc[g] = 3.0
    for b in 1:CLM.NUMRAD
        a2l.forc_solad_not_downscaled_grc[g, b] = 200.0; a2l.forc_solai_grc[g, b] = 80.0
    end
    a2l.forc_solar_not_downscaled_grc[g] = 560.0
    a2l.forc_rain_not_downscaled_grc[g] = 0.0001; a2l.forc_snow_not_downscaled_grc[g] = 0.0
end
CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
step!(n) = CLM.clm_drv!(config, inst, filt, filt_ia, bounds,
    true, nextsw_cday, declin, declin, CLM.ORB_OBLIQR_DEFAULT, false, false, "", false;
    nstep=n, is_first_step=(n==1), is_beg_curr_day=(n==1),
    dtime=1800.0, mon=1, day=1, photosyns=inst.photosyns)
println("Warmup (3 steps)..."); for n in 1:3; step!(n); end
flush(stdout)

# ---------------------------------------------------------------------------
# Shared setup mirroring the canopy_fluxes_core! call sites.
# ---------------------------------------------------------------------------
bounds_patch = bounds.begp:bounds.endp
begp = bounds.begp; endp = bounds.endp
dtime = 1800.0
pc = CLM.pftcon
mask_ev = filt.exposedvegp

# Crude-but-finite scratch (compile-time probe — numeric values irrelevant).
mk(v) = fill(v, endp)
svpts   = mk(1200.0)   # leaf saturation vapor pressure [Pa]
eah     = mk(900.0)    # canopy air vapor pressure [Pa]
o2_arr  = mk(21000.0)  # O2 partial pressure [Pa]
co2_arr = mk(40.0)     # CO2 partial pressure [Pa]
rb      = mk(40.0)     # leaf boundary layer resistance [s/m]
dayl_factor = mk(1.0)
downreg = zeros(np); leafn = zeros(np); o3v = zeros(np); o3g = zeros(np)
t10_patch   = inst.temperature.t_a10_patch

# forc_pbot patch-gathered from column (as canopy_fluxes_core! does)
forc_pbot_patch = zeros(endp)
for p in bounds_patch
    forc_pbot_patch[p] = a2l.forc_pbot_downscaled_col[inst.patch.column[p]]
end
ivt_vec = inst.patch.itype .+ 1

# ---------------------------------------------------------------------------
# Sub-call closures: each mutates inst and returns a scalar objective.
# ---------------------------------------------------------------------------
function call_photosynthesis!(inst)
    CLM.photosynthesis!(inst.photosyns,
        svpts, eah, o2_arr, co2_arr, rb,
        inst.energyflux.btran_patch, dayl_factor, leafn,
        forc_pbot_patch, inst.temperature.t_veg_patch, t10_patch,
        inst.temperature.thm_patch, inst.surfalb.nrad_patch,
        inst.surfalb.tlai_z_patch, inst.canopystate.tlai_patch,
        inst.solarabs.parsun_z_patch, inst.canopystate.laisun_z_patch,
        inst.surfalb.vcmaxcintsun_patch, o3v, o3g,
        pc.c3psn, pc.leafcn, pc.flnr, pc.fnitr, pc.slatop,
        pc.mbbopt, pc.medlynintercept, pc.medlynslope,
        ivt_vec, inst.patch.column,
        mask_ev, bounds_patch, "sun")
    return sum(abs2, inst.photosyns.psnsun_patch)
end

function call_friction_velocity!(inst)
    fv = inst.frictionvel
    fn = count(mask_ev[p] for p in bounds_patch)
    filt_arr = [p for p in bounds_patch if mask_ev[p]]
    ur    = mk(3.0)
    temp1 = mk(0.0); temp2 = mk(0.0); temp12m = mk(0.0); temp22m = mk(0.0); fm = mk(0.0)
    active = [mask_ev[p] for p in 1:endp]
    CLM.friction_velocity!(fv, fn, filt_arr,
        inst.canopystate.displa_patch, fv.z0mv_patch, fv.z0hv_patch, fv.z0qv_patch,
        fv.obu_patch, 1, ur, fv.um_patch, fv.ustar_patch,
        temp1, temp2, temp12m, temp22m, fm; active = active)
    return sum(abs2, fv.ustar_patch)
end

# Shared filter constants, built ONCE outside the differentiated closures (integer
# filter construction is not a differentiation target — Enzyme must not see it).
const FN_EV      = count(mask_ev[p] for p in bounds_patch)
const FILTERP_EV = Int[p for p in bounds_patch if mask_ev[p]]
const ACTIVE_EV  = Bool[mask_ev[p] for p in 1:endp]
params_cf = CLM.canopy_fluxes_params

function call_cf_resist!(inst)
    fn = FN_EV; filterp = FILTERP_EV; active = ACTIVE_EV
    z(args...) = zeros(args...)
    temp1=z(endp); temp2=z(endp); tlbef=z(endp); del2=z(endp); del_arr=z(endp)
    rah=z(endp,2); raw=z(endp,2); uuc=z(endp); rb=mk(40.0); svpts=mk(1200.0)
    eah=mk(900.0); el=mk(1000.0); grnd_ch4=z(endp)
    CLM.cf_resist_update!(inst.frictionvel, inst.canopystate, inst.temperature,
        inst.water.waterdiagnosticbulk_inst, inst.patch, filterp, fn, active,
        temp1, temp2, tlbef, del2, del_arr, rah, raw, uuc, rb, svpts, eah, el,
        pc.dleaf, grnd_ch4, a2l.forc_pbot_downscaled_col,
        params_cf.csoilc, false, false, false, params_cf)
    return sum(abs2, inst.frictionvel.ram1_patch)
end

function call_cf_energy!(inst)
    fn = FN_EV; filterp = FILTERP_EV; active = ACTIVE_EV
    z(args...) = zeros(args...)
    rah=z(endp,2); raw=z(endp,2); rb=mk(40.0); rstem=z(endp); sa_leaf=z(endp)
    sa_stem=z(endp); sa_internal=z(endp); frac=z(endp); air=z(endp); bir=z(endp)
    cir=z(endp); cp_leaf=z(endp); tl_ini=mk(285.0); tlbef=mk(285.0); zldis=mk(25.0)
    temp1=z(endp); temp2=z(endp); ur=mk(3.0); efeb=z(endp); wtg=z(endp); wtl0=z(endp)
    wta0=z(endp); wtstem0=z(endp); wtga=z(endp); wtal=z(endp); lw_stem=z(endp)
    lw_leaf=z(endp); wtgq=z(endp); wtlq0=z(endp); wtaq0=z(endp); wtalq=z(endp)
    efe=z(endp); dt_veg=z(endp); del_arr=z(endp); err_arr=z(endp); qsatl=mk(0.01)
    el=mk(1000.0); qsatldT=z(endp); dth=z(endp); dqh=z(endp); delq=z(endp)
    obuold=z(endp); nmozsgn=zeros(Int, endp)
    CLM.cf_energy_update!(inst.canopystate, inst.energyflux, inst.frictionvel,
        inst.temperature, inst.solarabs, inst.soilstate, inst.water.waterfluxbulk_inst,
        inst.water.waterstatebulk_inst, inst.water.waterdiagnosticbulk_inst,
        inst.photosyns, inst.patch, inst.column, filterp, fn, active,
        rah, raw, rb, rstem, sa_leaf, sa_stem, sa_internal, frac, air, bir, cir,
        cp_leaf, tl_ini, tlbef, zldis, temp1, temp2, ur, efeb, wtg, wtl0, wta0,
        wtstem0, wtga, wtal, lw_stem, lw_leaf, wtgq, wtlq0, wtaq0, wtalq, efe,
        dt_veg, del_arr, err_arr, qsatl, el, qsatldT, dth, dqh, delq, obuold, nmozsgn,
        a2l.forc_q_downscaled_col, a2l.forc_th_downscaled_col,
        a2l.forc_pbot_downscaled_col, a2l.forc_rho_downscaled_col,
        false, false, false, false, CLM.varpar.nlevsno, dtime, params_cf)
    return sum(abs2, inst.temperature.t_veg_patch)
end

# ---------------------------------------------------------------------------
# Drive each isolated probe.
# ---------------------------------------------------------------------------
Enzyme.API.strictAliasing!(false)

function probe(name, f, grad_field)
    println("\n", "="^70, "\n", name, "\n", "="^70)
    try
        v = f(inst); @printf("  primal value = %.6g\n", v)
    catch e
        println("  PRIMAL FAILED:"); showerror(stdout, e); println(); return
    end
    dinst = Enzyme.make_zero(inst)
    print("  Enzyme.Reverse ... "); flush(stdout)
    try
        Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.Reverse), f, Enzyme.Active,
                        Enzyme.Duplicated(inst, dinst))
        g = grad_field(dinst)
        @printf("SUCCESS — %d/%d nonzero shadow entries\n", count(!=(0.0), g), length(g))
    catch e
        println("BLOCKED:")
        bt = catch_backtrace()
        showerror(stdout, e, bt)
        println()
        # one-line signature for the summary
        msg = sprint(showerror, e)
        sig = first(split(msg, '\n'))
        println("  >>> SIGNATURE: ", sig)
    end
    flush(stdout)
end

# Full monolith call (the known failure) — same process, for an airtight A/B.
bc_col = bounds.begc:bounds.endc
forc_q_col = zeros(nc)
CLM.compute_forc_q!(forc_q_col, inst.column.gridcell, a2l.forc_vp_grc, a2l.forc_pbot_downscaled_col)
function call_monolith!(inst)
    CLM.canopy_fluxes_core!(inst.canopystate, inst.energyflux, inst.frictionvel,
        inst.temperature, inst.solarabs, inst.soilstate,
        inst.water.waterfluxbulk_inst, inst.water.waterstatebulk_inst,
        inst.water.waterdiagnosticbulk_inst, inst.photosyns,
        inst.patch, inst.column, inst.gridcell,
        mask_ev, bounds_patch, bc_col,
        a2l.forc_lwrad_downscaled_col, forc_q_col, a2l.forc_pbot_downscaled_col,
        a2l.forc_th_downscaled_col, a2l.forc_rho_downscaled_col, a2l.forc_t_downscaled_col,
        a2l.forc_u_grc, a2l.forc_v_grc, a2l.forc_pco2_grc, a2l.forc_po2_grc,
        a2l.forc_hgt_t_grc, a2l.forc_hgt_u_grc, a2l.forc_hgt_q_grc,
        inst.gridcell.dayl, inst.gridcell.max_dayl, downreg, leafn, dtime,
        inst.temperature.t_a10_patch, inst.surfalb.nrad_patch, inst.surfalb.tlai_z_patch,
        inst.surfalb.vcmaxcintsun_patch, inst.surfalb.vcmaxcintsha_patch,
        inst.solarabs.parsun_z_patch, inst.solarabs.parsha_z_patch,
        inst.canopystate.laisun_z_patch, inst.canopystate.laisha_z_patch,
        o3v, o3g, pc.dleaf, pc.slatop, pc.leafcn, pc.flnr, pc.fnitr, pc.mbbopt,
        pc.c3psn, pc.woody, inst.overrides)
    return sum(abs2, inst.temperature.t_veg_patch)
end

probe("photosynthesis! (sun)", call_photosynthesis!,
      d -> d.photosyns.psnsun_patch)
probe("friction_velocity!", call_friction_velocity!,
      d -> d.frictionvel.ustar_patch)
probe("cf_resist_update!  [struct-arg kernel]", call_cf_resist!,
      d -> d.frictionvel.ram1_patch)
probe("cf_energy_update!  [struct-arg kernel]", call_cf_energy!,
      d -> d.temperature.t_veg_patch)
probe("canopy_fluxes_core! [MONOLITH — known wall]", call_monolith!,
      d -> d.temperature.t_veg_patch)

println("\nDone.")
