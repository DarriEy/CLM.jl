# probe_h2osfc_subdaily.jl — dump per-timestep surface-water state for a domain.
#
# Confirms the FH2OSFC diurnal-cycle hypothesis: records h2osfc, frac_h2osfc and
# the surface-water in/out fluxes every timestep (not daily-averaged) for the
# soil column, so we can see whether Julia's within-day h2osfc is flat (steady
# ponding) vs cycles diurnally like Fortran.
#
#   DOMAIN=Krycklan julia +1.12 --project=. scripts/probe_h2osfc_subdaily.jl
#
# Output: scripts/../paper/data/probe_h2osfc_<domain>.csv
# Reuses the parity_run_domain.jl per-domain registry + config.

using NCDatasets, Dates, Printf, CLM
const run_clm! = getfield(CLM, Symbol("clm_run!"))
const DATA = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data"

_snicar(f) = let d1 = "$DATA/installs/cesm-inputdata/lnd/clm2/snicardata",
                 d2 = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata"
    isfile(joinpath(d1, f)) ? joinpath(d1, f) : joinpath(d2, f)
end

# Minimal per-domain registry (subset of parity_run_domain.jl fields we need).
const DOMAINS = Dict(
    "Krycklan" => (year = 2013, dtime = 3600, baseflow = 0.001, int_snow = 2000.0,
        caldir  = "$DATA/domain_Boreal_Krycklan_Sweden/settings/CLM/parameters",
        forcing = "$DATA/domain_Boreal_Krycklan_Sweden/data/forcing/CLM_input/clmforc.2013.nc",
        restart = "$DATA/domain_Boreal_Krycklan_Sweden/simulations/clm_boreal/CLM/Boreal_Krycklan_Sweden.clm2.r.2013-01-01-00000.nc"),
    "Abisko" => (year = 2013, dtime = 3600, baseflow = 0.001, int_snow = 2000.0,
        caldir  = "$DATA/domain_Arctic_Abisko_Sweden/settings/CLM/parameters",
        forcing = "$DATA/domain_Arctic_Abisko_Sweden/data/forcing/CLM_input/clmforc.2013.nc",
        restart = "$DATA/domain_Arctic_Abisko_Sweden/simulations/clm_arctic/CLM/Arctic_Abisko_Sweden.clm2.r.2013-01-01-00000.nc"),
)

const DOM = get(ENV, "DOMAIN", "Krycklan")
cfg = DOMAINS[DOM]
yr  = cfg.year
fsurdat   = joinpath(cfg.caldir, "surfdata_clm.nc")
paramfile = joinpath(cfg.caldir, "clm5_params.nc")
outdir    = abspath(joinpath(@__DIR__, "..", "paper", "data"))
fhistory  = joinpath(outdir, "probe_hist_$(lowercase(DOM)).nc")

CLM.rooting_profile_config.rooting_profile_method_water  = CLM.JACKSON_1996_ROOT
CLM.rooting_profile_config.rooting_profile_method_carbon = CLM.JACKSON_1996_ROOT
CLM.snow_hydrology_set_control_for_testing!(;
    wind_dep_snow_density = true,
    overburden_compaction_method = CLM.OVERBURDEN_COMPACTION_VIONNET2012)

# Per-timestep records for the (first) soil column.
recs = NamedTuple[]
csoil = Ref(0)

function probe(inst, tm)
    # Locate a soil column once (istsoil landunit).
    if csoil[] == 0
        for c in eachindex(inst.column.snl)
            csoil[] = c; break
        end
    end
    c = csoil[]
    ws  = inst.water.waterstatebulk_inst.ws
    wd  = inst.water.waterdiagnosticbulk_inst
    wfb = inst.water.waterfluxbulk_inst
    wf  = wfb.wf
    getc(a) = (c <= length(a) ? Float64(a[c]) : NaN)
    push!(recs, (
        date   = tm.current_date,
        h2osfc = getc(ws.h2osfc_col),
        frac   = getc(wd.frac_h2osfc_col),
        q_in   = getc(wfb.qflx_top_soil_to_h2osfc_col),  # input to h2osfc (mm/s)
        q_drain= getc(wfb.qflx_h2osfc_drain_col),        # h2osfc -> soil (mm/s)
        q_spill= getc(wfb.qflx_h2osfc_surf_col),         # h2osfc -> runoff (mm/s)
        q_melt = getc(wf.qflx_snomelt_col),              # snowmelt (mm/s)
        q_infl = getc(wf.qflx_infl_col),                 # total infiltration (mm/s)
    ))
end

println("Probe run: $DOM $yr (per-timestep surface water)")
t0 = time()
run_clm!(;
    fsurdat = fsurdat, paramfile = paramfile, fforcing = cfg.forcing,
    fhistory = fhistory,
    start_date = DateTime(yr, 1, 1), end_date = DateTime(yr + 1, 1, 1),
    dtime = cfg.dtime, use_cn = false, verbose = false,
    use_aquifer_layer = false, use_hydrstress = true, use_luna = true,
    h2osfcflag = 1,
    baseflow_scalar = cfg.baseflow, int_snow_max = cfg.int_snow,
    ffortran_restart = cfg.restart, interp_forcing = true,
    fsnowoptics = isfile(_snicar("snicar_optics_5bnd_c013122.nc")) ? _snicar("snicar_optics_5bnd_c013122.nc") : "",
    fsnowaging  = isfile(_snicar("snicar_drdt_bst_fit_60_c070416.nc")) ? _snicar("snicar_drdt_bst_fit_60_c070416.nc") : "",
    step_probe = probe)

csv = joinpath(outdir, "probe_h2osfc_$(lowercase(DOM)).csv")
open(csv, "w") do io
    println(io, "date,h2osfc,frac,q_in,q_drain,q_spill,q_melt,q_infl")
    for r in recs
        @printf(io, "%s,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n",
            r.date, r.h2osfc, r.frac, r.q_in, r.q_drain, r.q_spill, r.q_melt, r.q_infl)
    end
end
@printf("Done %s in %.1fs -> %s (%d steps)\n", DOM, time() - t0, csv, length(recs))
