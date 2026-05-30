# Granularity de-risk for compositional Enzyme: can a single large CLM phase be
# reverse-mode compiled in isolation without segfaulting?
#
# Target: soil_temperature! — a standalone phase with a clean positional
# signature that also exercises the custom band_solve! adjoint rule in context.
# We wrap it as an inst-closure and run Enzyme.Reverse over Duplicated(inst).
#
#   julia --project=. scripts/enzyme_phase_probe.jl

using CLM
using Enzyme
using Printf

const FSURDAT = get(ENV, "CLM_FSURDAT",
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc")
const PARAMFILE = get(ENV, "CLM_PARAMFILE",
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc")

function setup_forcing!(a2l, T, ng)
    for g in 1:ng
        a2l.forc_t_not_downscaled_grc[g] = T
        a2l.forc_pbot_not_downscaled_grc[g] = 85000.0
        a2l.forc_th_not_downscaled_grc[g] = T * (100000.0 / 85000.0)^(CLM.RAIR / CLM.CPAIR)
        a2l.forc_rho_not_downscaled_grc[g] = 85000.0 / (CLM.RAIR * T)
        a2l.forc_lwrad_not_downscaled_grc[g] = 300.0
        a2l.forc_vp_grc[g] = 800.0
        a2l.forc_hgt_grc[g] = 30.0
        a2l.forc_topo_grc[g] = 0.0
        a2l.forc_wind_grc[g] = 3.0
        for b in 1:CLM.NUMRAD
            a2l.forc_solad_not_downscaled_grc[g, b] = 200.0
            a2l.forc_solai_grc[g, b] = 80.0
        end
        a2l.forc_solar_not_downscaled_grc[g] = 560.0
        a2l.forc_rain_not_downscaled_grc[g] = 0.0001
        a2l.forc_snow_not_downscaled_grc[g] = 0.0
    end
end

println("Initializing + warmup...")
(inst, bounds, filt, tm) = CLM.clm_initialize!(; fsurdat=FSURDAT, paramfile=PARAMFILE)
ng = bounds.endg; nl = bounds.endl
config = CLM.CLMDriverConfig()
filt_ia = CLM.clump_filter_inactive_and_active
(declin, eccf) = CLM.compute_orbital(120.0)
nextsw_cday = 120.0 + 1800.0 / CLM.SECSPDAY
setup_forcing!(inst.atm2lnd, 285.0, ng)
CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
for n in 1:3
    CLM.clm_drv!(config, inst, filt, filt_ia, bounds, true, nextsw_cday, declin, declin,
        CLM.ORB_OBLIQR_DEFAULT, false, false, "", false;
        nstep=n, is_first_step=(n == 1), is_beg_curr_day=(n == 1),
        dtime=1800.0, mon=1, day=1, photosyns=inst.photosyns)
end

# Driver-local constants for the soil_temperature! call (Const for Enzyme).
bc_col = bounds.begc:bounds.endc
bc_lun = bounds.begl:bounds.endl
bc_patch = bounds.begp:bounds.endp
urbantv = fill(323.15, nl)
dtime = 1800.0

# Phase as an inst-closure: differentiate soil_temperature! alone.
function phase!(inst)
    CLM.soil_temperature!(inst.column, inst.landunit, inst.patch, inst.temperature,
        inst.energyflux, inst.soilstate, inst.water.waterstatebulk_inst,
        inst.water.waterdiagnosticbulk_inst, inst.water.waterfluxbulk_inst,
        inst.solarabs, inst.canopystate, inst.urbanparams,
        urbantv, inst.atm2lnd.forc_lwrad_downscaled_col,
        filt.nolakec, filt.nolakep, filt.urbanl, filt.urbanc,
        bc_col, bc_lun, bc_patch, dtime)
    return sum(abs2, inst.temperature.t_soisno_col)
end

println("phase! primal value = ", phase!(inst))

dinst = Enzyme.make_zero(inst)
println("Attempting Enzyme.Reverse on soil_temperature! alone...")
Enzyme.API.strictAliasing!(false)
if get(ENV, "CLM_ENZYME_VERBOSE", "0") == "1"
    Enzyme.Compiler.VERBOSE_ERRORS[] = true
end
flush(stdout)
try
    Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.Reverse), phase!, Enzyme.Active,
                    Enzyme.Duplicated(inst, dinst))
    g = dinst.temperature.t_soisno_col
    @printf("\nSUCCESS: soil_temperature! reverse-mode compiled.\n")
    @printf("  d(sum t_soisno^2)/d(t_soisno_in)[1,1] = %.6f\n", g[1, 1])
    @printf("  nonzero shadow entries = %d / %d\n", count(!=(0.0), g), length(g))
catch e
    println("\nPHASE BLOCKER:")
    showerror(stdout, e, catch_backtrace())
    println()
end
