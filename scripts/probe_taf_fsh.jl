# probe_taf_fsh.jl — per-timestep canopy-air-temp / resistance dump for the FSH_V
# (vegetation sensible heat) coupled-taf Monin-Obukhov chase.
#
# Dumps, for every exposed-veg patch each timestep, the converged canopy quantities
# that drive eflx_sh_veg: taf, t_veg, t_grnd, thm, rah_above, rah_below, rb,
# elai, esai, dt_veg, eflx_sh_veg, eflx_sh_grnd. Match to the Fortran PARITYTAF
# stdout dump by date+sec+itype to find which term diverges.
#
#   DOMAIN=Savanna julia +1.12 --project=. scripts/probe_taf_fsh.jl
#
# Output: paper/data/probe_taf_<domain>.csv

using NCDatasets, Dates, Printf, CLM
const run_clm! = getfield(CLM, Symbol("clm_run!"))
const DATA = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data"
_snicar(f) = let d1 = "$DATA/installs/cesm-inputdata/lnd/clm2/snicardata",
                 d2 = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata"
    isfile(joinpath(d1, f)) ? joinpath(d1, f) : joinpath(d2, f)
end

const DOMAINS = Dict(
    "Savanna" => (year = 2017, dtime = 3600, baseflow = 0.001, int_snow = 2000.0,
        caldir  = "$DATA/domain_Savanna_Donga_Benin/settings/CLM/parameters",
        forcing = "$DATA/domain_Savanna_Donga_Benin/data/forcing/CLM_input/clmforc.2017.nc",
        restart = "$DATA/domain_Savanna_Donga_Benin/simulations/clm_savanna/CLM/Savanna_Donga_Benin.clm2.r.2017-01-01-00000.nc"),
    "Tagus" => (year = 2013, dtime = 3600, baseflow = 0.001, int_snow = 2000.0,
        caldir  = "$DATA/domain_Mediterranean_Tagus_Spain/settings/CLM/parameters",
        forcing = "$DATA/domain_Mediterranean_Tagus_Spain/data/forcing/CLM_input/clmforc.2013.nc",
        restart = "$DATA/domain_Mediterranean_Tagus_Spain/simulations/clm_mediterranean/CLM/Mediterranean_Tagus_Spain.clm2.r.2013-01-01-00000.nc"),
    "Stillwater" => (year = 2003, dtime = 3600, baseflow = 0.016035343416707443, int_snow = 2000.0,
        caldir  = "$DATA/domain_Stillwater_Oklahoma/optimization/CLM/dds_clm_dds_calibration/final_evaluation/settings/CLM/parameters",
        forcing = "$DATA/domain_Stillwater_Oklahoma/data/forcing/CLM_input/clmforc.2003.nc",
        restart = "$DATA/domain_Stillwater_Oklahoma/optimization/CLM/dds_clm_dds_calibration/final_evaluation/Stillwater_Oklahoma.clm2.r.2003-01-01-00000.nc"),
)

const DOM = get(ENV, "DOMAIN", "Savanna")
cfg = DOMAINS[DOM]
yr  = cfg.year
fsurdat   = joinpath(cfg.caldir, "surfdata_clm.nc")
paramfile = joinpath(cfg.caldir, "clm5_params.nc")
outdir    = abspath(joinpath(@__DIR__, "..", "paper", "data"))
fhistory  = joinpath(outdir, "probe_taf_hist_$(lowercase(DOM)).nc")
csvpath   = joinpath(outdir, "probe_taf_$(lowercase(DOM)).csv")

CLM.rooting_profile_config.rooting_profile_method_water  = CLM.JACKSON_1996_ROOT
CLM.rooting_profile_config.rooting_profile_method_carbon = CLM.JACKSON_1996_ROOT
CLM.snow_hydrology_set_control_for_testing!(;
    wind_dep_snow_density = true,
    overburden_compaction_method = CLM.OVERBURDEN_COMPACTION_VIONNET2012)

io = open(csvpath, "w")
println(io, "yr,mon,day,sec,itype,wt,taf,t_veg,t_grnd,thm,wta0,wtl0,wtg0,wtga,rah1,rah2,rb,elai,esai,dt_veg,eflx_sh_veg,eflx_sh_grnd,uaf,sabv,ram1,ustar")

function probe(inst, tm)
    d = tm.current_date
    sec = Dates.hour(d) * 3600 + Dates.minute(d) * 60 + Dates.second(d)
    fv = inst.frictionvel; ef = inst.energyflux; cs = inst.canopystate; sa = inst.solarabs
    tp = inst.temperature; pch = inst.patch
    gp(a, p) = (p <= length(a) ? Float64(a[p]) : NaN)
    for p in eachindex(pch.itype)
        it = pch.itype[p]
        it == 0 && continue                       # skip bare ground
        pch.wtgcell[p] > 0.0 || continue
        elai = gp(cs.elai_patch, p); esai = gp(cs.esai_patch, p)
        (isfinite(elai) && (elai + esai) > 0.0) || continue
        c = Int(pch.column[p])
        rah1 = gp(fv.rah1_patch, p); rah2 = gp(fv.rah2_patch, p); rb = gp(fv.rb1_patch, p)
        # Reconstruct the normalized heat-conductance weights from the resistances,
        # exactly as canopy_fluxes does: wta=1/rah_above, wtg=1/rah_below,
        # wtl=(elai+esai)/rb, wtstem=0 (biomass heat storage off). wtshi normalizes.
        wta = 1.0 / rah1; wtg = 1.0 / rah2; wtl = (elai + esai) / rb
        wtshi = 1.0 / (wta + wtl + wtg)
        wta0 = wta * wtshi; wtg0 = wtg * wtshi; wtl0 = wtl * wtshi
        wtga = wta0 + wtg0
        @printf(io, "%d,%d,%d,%d,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.8f,%.8f,%.8f,%.8f,%.6f,%.6f,%.6f,%.6f,%.6f,%.8f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
            Dates.year(d), Dates.month(d), Dates.day(d), sec, Int(it), gp(pch.wtgcell, p),
            gp(fv.taf_patch, p), gp(tp.t_veg_patch, p), gp(tp.t_grnd_col, c), gp(tp.thm_patch, p),
            wta0, wtl0, wtg0, wtga, rah1, rah2, rb, elai, esai,
            gp(tp.dt_veg_patch, p), gp(ef.eflx_sh_veg_patch, p), gp(ef.eflx_sh_grnd_patch, p),
            gp(fv.uaf_patch, p), gp(sa.sabv_patch, p), gp(fv.ram1_patch, p), gp(fv.ustar_patch, p))
    end
end

t0 = time()
run_clm!(;
    fsurdat = fsurdat, paramfile = paramfile, fforcing = cfg.forcing, fhistory = fhistory,
    start_date = DateTime(yr, 1, 1), end_date = DateTime(yr + 1, 1, 1),
    dtime = cfg.dtime, use_cn = false, verbose = false,
    use_aquifer_layer = false,
    use_hydrstress = true, use_luna = true, h2osfcflag = 1,
    baseflow_scalar = cfg.baseflow, int_snow_max = cfg.int_snow,
    ffortran_restart = cfg.restart, interp_forcing = true,
    fsnowoptics = isfile(_snicar("snicar_optics_5bnd_c013122.nc")) ? _snicar("snicar_optics_5bnd_c013122.nc") : "",
    fsnowaging  = isfile(_snicar("snicar_drdt_bst_fit_60_c070416.nc")) ? _snicar("snicar_drdt_bst_fit_60_c070416.nc") : "",
    step_probe = probe)
close(io)
@printf("Done: %s in %.1f s -> %s\n", DOM, time() - t0, csvpath)
