# Exploratory probe: attempt Enzyme reverse-mode through the full clm_drv!.
#
# Differentiates t_grnd (one timestep, after warmup) w.r.t. the forcing
# temperature, using Duplicated(inst) reverse-mode. The custom band-solve rule
# (enzyme_rules.jl) handles the soil-temperature LAPACK solve.
#
# This is HIGH-RISK / exploratory: Enzyme compiling reverse-mode through ~89K
# LOC may take very long, OOM, or error on an unsupported construct. The point
# is to capture the FIRST real blocker (or success), so run with a timeout.
#
#   julia --project=. scripts/enzyme_fulldriver_probe.jl

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

println("Initializing...")
(inst, bounds, filt, tm) = CLM.clm_initialize!(; fsurdat=FSURDAT, paramfile=PARAMFILE)
ng = bounds.endg
config = CLM.CLMDriverConfig()
filt_ia = CLM.clump_filter_inactive_and_active
(declin, eccf) = CLM.compute_orbital(120.0)
nextsw_cday = 120.0 + 1800.0 / CLM.SECSPDAY
T0 = 285.0
setup_forcing!(inst.atm2lnd, T0, ng)
CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)

step!(inst, n) = CLM.clm_drv!(config, inst, filt, filt_ia, bounds,
    true, nextsw_cday, declin, declin, CLM.ORB_OBLIQR_DEFAULT,
    false, false, "", false;
    nstep=n, is_first_step=(n == 1), is_beg_curr_day=(n == 1),
    dtime=1800.0, mon=1, day=1, photosyns=inst.photosyns)

println("Warmup (3 steps)...")
for n in 1:3
    step!(inst, n)
end

# Scalar objective for reverse-mode: one more timestep, read t_grnd.
# Calls the fully-positional clm_drv_core! (no kwargs) so Enzyme need not build
# an augmented kwarg NamedTuple containing the active photosyns struct.
function objective!(inst)
    CLM.clm_drv_core!(config, inst, filt, filt_ia, bounds,
        true, nextsw_cday, declin, declin, CLM.ORB_OBLIQR_DEFAULT,
        false, false, "", false,
        4, false, false, false, false, 1800.0, 1, 1, inst.photosyns)
    return inst.temperature.t_grnd_col[1]
end

println("Building zero shadow (make_zero)...")
dinst = Enzyme.make_zero(inst)

println("Attempting Enzyme.Reverse through clm_drv!  (this may take a while)...")
# strictAliasing(false): CLM kernels hit small Union types in inference that
# Enzyme's strict type analysis rejects (IllegalTypeAnalysisException). This is
# the documented workaround.
Enzyme.API.strictAliasing!(false)
flush(stdout)
try
    # set_runtime_activity: many CLM kernels bundle constant PFT params with
    # active state in the same NamedTuple/struct, which static activity analysis
    # cannot separate. Runtime activity handles this (slight perf cost).
    Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.Reverse), objective!, Enzyme.Active,
                    Enzyme.Duplicated(inst, dinst))
    g = dinst.atm2lnd.forc_t_not_downscaled_grc[1]
    @printf("\nSUCCESS: d(t_grnd)/d(forc_t) [reverse] = %.6f\n", g)
catch e
    println("\nFIRST BLOCKER:")
    showerror(stdout, e, catch_backtrace())
    println()
end
