# probe_luna.jl — dump Julia's persisted LUNA 10-day inputs + converged vcmx25_z
# once per day, to compare against the Fortran LUNADUMP (LunaMod.F90 SourceMod).
#
#   DOMAIN=Stillwater julia +1.12 --project=. scripts/probe_luna.jl
# Output: paper/data/probe_luna_<domain>.csv

using NCDatasets, Dates, Printf, CLM
const run_clm! = getfield(CLM, Symbol("clm_run!"))
const DATA = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data"
_snicar(f) = let d1="$DATA/installs/cesm-inputdata/lnd/clm2/snicardata", d2="/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata"
    isfile(joinpath(d1,f)) ? joinpath(d1,f) : joinpath(d2,f) end

const DOMAINS = Dict(
    "Stillwater" => (year=2003, dtime=3600, baseflow=0.016035343416707443, int_snow=2000.0,
        caldir="$DATA/domain_Stillwater_Oklahoma/optimization/CLM/dds_clm_dds_calibration/final_evaluation/settings/CLM/parameters",
        forcing="$DATA/domain_Stillwater_Oklahoma/data/forcing/CLM_input/clmforc.2003.nc",
        restart="$DATA/domain_Stillwater_Oklahoma/optimization/CLM/dds_clm_dds_calibration/final_evaluation/Stillwater_Oklahoma.clm2.r.2003-01-01-00000.nc"),
)
const DOM = get(ENV, "DOMAIN", "Stillwater"); cfg = DOMAINS[DOM]; yr = cfg.year
outdir = abspath(joinpath(@__DIR__, "..", "paper", "data"))
csvpath = joinpath(outdir, "probe_luna_$(lowercase(DOM)).csv")

CLM.rooting_profile_config.rooting_profile_method_water  = CLM.JACKSON_1996_ROOT
CLM.rooting_profile_config.rooting_profile_method_carbon = CLM.JACKSON_1996_ROOT
CLM.snow_hydrology_set_control_for_testing!(; wind_dep_snow_density=true,
    overburden_compaction_method=CLM.OVERBURDEN_COMPACTION_VIONNET2012)

const TFRZ = 273.15
io = open(csvpath, "w")
println(io, "yr,mon,day,itype,vcmx25_z,tleafd10,tleafn10,tair10,relh10,par240d,CO2a10,forc_pbot10,rb10")
lastday = Ref(-1)

function probe(inst, tm)
    d = tm.current_date
    Dates.day(d) == lastday[] && return
    lastday[] = Dates.day(d)
    ps=inst.photosyns; tp=inst.temperature; wd=inst.water.waterdiagnosticbulk_inst
    sa=inst.solarabs; a2l=inst.atm2lnd; fv=inst.frictionvel; pch=inst.patch
    gp(a,p)=(p<=length(a) ? Float64(a[p]) : NaN)
    for p in eachindex(pch.itype)
        it=pch.itype[p]; it==0 && continue
        pch.wtgcell[p]>0.0 || continue
        vz = (size(ps.vcmx25_z_patch,1)>=p) ? Float64(ps.vcmx25_z_patch[p,1]) : NaN
        (isfinite(vz) && vz>0) || continue
        @printf(io, "%d,%d,%d,%d,%.4f,%.4f,%.4f,%.4f,%.5f,%.4f,%.5f,%.2f,%.4f\n",
            Dates.year(d),Dates.month(d),Dates.day(d),Int(it), vz,
            gp(tp.t_veg10_day_patch,p)-TFRZ, gp(tp.t_veg10_night_patch,p)-TFRZ,
            gp(tp.t_a10_patch,p)-TFRZ, min(1.0,gp(wd.rh10_af_patch,p)),
            (size(sa.par240d_z_patch,1)>=p ? Float64(sa.par240d_z_patch[p,1]) : NaN),
            gp(a2l.forc_pco2_240_patch,p), gp(a2l.forc_pbot240_downscaled_patch,p), gp(fv.rb10_patch,p))
    end
end

run_clm!(; fsurdat=joinpath(cfg.caldir,"surfdata_clm.nc"), paramfile=joinpath(cfg.caldir,"clm5_params.nc"),
    fforcing=cfg.forcing, fhistory=joinpath(outdir,"probe_luna_hist_$(lowercase(DOM)).nc"),
    start_date=DateTime(yr,1,1), end_date=DateTime(yr+1,1,1), dtime=cfg.dtime, use_cn=false, verbose=false,
    use_hydrstress=true, use_luna=true, h2osfcflag=1,
    baseflow_scalar=cfg.baseflow, int_snow_max=cfg.int_snow, ffortran_restart=cfg.restart, interp_forcing=true,
    fsnowoptics=isfile(_snicar("snicar_optics_5bnd_c013122.nc")) ? _snicar("snicar_optics_5bnd_c013122.nc") : "",
    fsnowaging =isfile(_snicar("snicar_drdt_bst_fit_60_c070416.nc")) ? _snicar("snicar_drdt_bst_fit_60_c070416.nc") : "",
    step_probe=probe)
close(io); @printf("Done LUNA probe: %s -> %s\n", DOM, csvpath)
