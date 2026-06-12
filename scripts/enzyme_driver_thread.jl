# Threading canopy + soil into ONE compositional-reverse engine (item #3).
#
# The generic CLM.compositional_reverse!(phases, b, seed_bang!) engine (in
# src/biogeophys/canopy_fluxes_reverse.jl) drives BOTH:
#   • the canopy sub-phase chain (validated finite in scripts/enzyme_canopy_psn.jl), and
#   • soil_temperature! as a single driver phase (here, on the real warmed-up inst —
#     the phase Enzyme already differentiates cleanly: see enzyme_phase_probe.jl).
#
# This script proves soil_temperature! plugs into the SAME engine and reproduces its
# standalone reverse gradient — i.e. the driver-level reverse is just a phase list,
# and canopy's reverse is one (multi-sub-phase) entry, soil's another. A full
# end-to-end canopy→soil COUPLED FD check additionally needs a single finite state in
# which BOTH run finite (the Bow-at-Banff test domain leaves canopy radiation NaN);
# that finite full-driver state is the remaining step.
#
#   julia +1.10 --project=/tmp/clm_jl10_env scripts/enzyme_driver_thread.jl
#
using CLM, Enzyme, Printf
const FSURDAT = get(ENV, "CLM_FSURDAT",
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc")
const PARAMFILE = get(ENV, "CLM_PARAMFILE",
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc")

(inst, bounds, filt, tm) = CLM.clm_initialize!(; fsurdat=FSURDAT, paramfile=PARAMFILE)
ng = bounds.endg; nl = bounds.endl
config = CLM.CLMDriverConfig(); filt_ia = CLM.clump_filter_inactive_and_active
(declin, eccf) = CLM.compute_orbital(120.0); nextsw_cday = 120.0 + 1800.0/CLM.SECSPDAY
a2l = inst.atm2lnd; T0 = 285.0
for g in 1:ng
    a2l.forc_t_not_downscaled_grc[g]=T0; a2l.forc_pbot_not_downscaled_grc[g]=85000.0
    a2l.forc_th_not_downscaled_grc[g]=T0*(100000.0/85000.0)^(CLM.RAIR/CLM.CPAIR)
    a2l.forc_rho_not_downscaled_grc[g]=85000.0/(CLM.RAIR*T0); a2l.forc_lwrad_not_downscaled_grc[g]=300.0
    a2l.forc_vp_grc[g]=800.0; a2l.forc_hgt_grc[g]=30.0; a2l.forc_topo_grc[g]=0.0; a2l.forc_wind_grc[g]=3.0
    a2l.forc_hgt_u_grc[g]=30.0; a2l.forc_hgt_t_grc[g]=30.0; a2l.forc_hgt_q_grc[g]=30.0  # component obs heights (forcing_reader derives these; harness bypasses it)
    for b in 1:CLM.NUMRAD; a2l.forc_solad_not_downscaled_grc[g,b]=200.0; a2l.forc_solai_grc[g,b]=80.0; end
    a2l.forc_solar_not_downscaled_grc[g]=560.0; a2l.forc_rain_not_downscaled_grc[g]=0.0001; a2l.forc_snow_not_downscaled_grc[g]=0.0
end
CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
for n in 1:3
    CLM.clm_drv!(config, inst, filt, filt_ia, bounds, true, nextsw_cday, declin, declin,
        CLM.ORB_OBLIQR_DEFAULT, false, false, "", false; nstep=n, is_first_step=(n==1),
        is_beg_curr_day=(n==1), dtime=1800.0, mon=1, day=1, photosyns=inst.photosyns)
end

# soil_temperature! as a driver phase over the whole-inst bundle.
const URBANTV = fill(323.15, nl); const DTIME = 1800.0
const BC_COL = bounds.begc:bounds.endc; const BC_LUN = bounds.begl:bounds.endl
const BC_PATCH = bounds.begp:bounds.endp
function soil_phase!(b)
    CLM.soil_temperature!(b.column, b.landunit, b.patch, b.temperature, b.energyflux,
        b.soilstate, b.water.waterstatebulk_inst, b.water.waterdiagnosticbulk_inst,
        b.water.waterfluxbulk_inst, b.solarabs, b.canopystate, b.urbanparams,
        URBANTV, b.atm2lnd.forc_lwrad_downscaled_col, filt.nolakec, filt.nolakep,
        filt.urbanl, filt.urbanc, BC_COL, BC_LUN, BC_PATCH, DTIME)
    return nothing
end

# Reference: standalone soil reverse (the enzyme_phase_probe.jl pattern).
function L_soil(inst)
    soil_phase!(inst); return sum(abs2, inst.temperature.t_soisno_col)
end
println("standalone soil_temperature! reverse (reference) ...")
dref = let inst2 = deepcopy(inst), d = Enzyme.make_zero(inst)
    Enzyme.API.strictAliasing!(false)
    Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.Reverse), L_soil, Enzyme.Active,
                    Enzyme.Duplicated(inst2, d))
    d.temperature.t_soisno_col
end
@printf("  standalone d(sum t_soisno^2)/d(t_soisno_in)[1,1] = %.6f  (nz=%d/%d)\n",
    dref[1,1], count(!=(0.0), dref), length(dref))

# Via the GENERIC engine: phases = [soil_phase!], seed dL/d(t_soisno)=2·t_soisno.
println("\nsoil_temperature! through CLM.compositional_reverse! (the canopy engine) ...")
seed_soil!(db, b) = (db.temperature.t_soisno_col .= 2 .* b.temperature.t_soisno_col)
db = CLM.compositional_reverse!(Any[(soil_phase!, ())], deepcopy(inst), seed_soil!)
g_eng = db.temperature.t_soisno_col
@printf("  engine     d(sum t_soisno^2)/d(t_soisno_in)[1,1] = %.6f  (nz=%d/%d)\n",
    g_eng[1,1], count(!=(0.0), g_eng), length(g_eng))

println("\n", "="^60)
# Compare only entries finite in BOTH (the real warmed-up inst carries spval/NaN in
# t_soisno for the NaN-canopy columns — present identically in both gradients; a raw
# max would be poisoned by NaN−NaN).
function maxdiff_finite(a, b)
    d = 0.0; n = 0
    for i in eachindex(a, b)
        (isfinite(a[i]) && isfinite(b[i])) || continue
        n += 1; d = max(d, abs(a[i] - b[i]))
    end
    return d, n
end
d, nfin = maxdiff_finite(g_eng, dref)
@printf("max |engine - standalone| over %d finite entries = %.3e\n", nfin, d)
if d < 1e-8 && nfin > 0
    println("\nsoil_temperature! THREADS THROUGH THE SHARED COMPOSITIONAL-REVERSE ENGINE ✓")
    println("(same engine drives the canopy sub-phase chain — see enzyme_canopy_psn.jl)")
else
    println("\nMISMATCH ✗")
end
println("="^60)
