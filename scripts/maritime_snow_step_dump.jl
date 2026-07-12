#!/usr/bin/env julia

using CLM
using Dates
using NCDatasets
using Printf

const run_clm! = getfield(CLM, Symbol("clm_run!"))
const DATA = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data"
const DOMAIN = "$DATA/domain_Maritime_HJAndrews_USA"
const CALDIR = "$DOMAIN/settings/CLM/parameters"
const FSURDAT = joinpath(CALDIR, "surfdata_clm.nc")
const FPARAM = joinpath(CALDIR, "clm5_params.nc")
const FFORCING = "$DOMAIN/data/forcing/CLM_input/clmforc.2017.nc"
const FRESTART = "$DOMAIN/simulations/clm_maritime/CLM/Maritime_HJAndrews_USA.clm2.r.2017-01-01-00000.nc"

_snicar(f) = let d1 = "$DATA/installs/cesm-inputdata/lnd/clm2/snicardata",
                 d2 = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata"
    isfile(joinpath(d1, f)) ? joinpath(d1, f) : joinpath(d2, f)
end

const OUTDIR = get(ENV, "JULIA_SNOW_DUMPDIR", "/tmp/hj_julia_snow")
const FORTRAN_OFFSET = parse(Int, get(ENV, "FORTRAN_NSTEP_OFFSET", "8760"))
const GLOBAL_LO = parse(Int, get(ENV, "JULIA_DUMP_GLOBAL_LO", "9900"))
const GLOBAL_HI = parse(Int, get(ENV, "JULIA_DUMP_GLOBAL_HI", "10960"))
const LOCAL_LO = GLOBAL_LO - FORTRAN_OFFSET
const LOCAL_HI = GLOBAL_HI - FORTRAN_OFFSET

LOCAL_LO >= 1 || error("LOCAL_LO=$LOCAL_LO; check FORTRAN_NSTEP_OFFSET/GLOBAL_LO")
isdir(OUTDIR) || mkpath(OUTDIR)

function put1d!(ds, name, value)
    v = defVar(ds, name, Float64, ("column",))
    v[:] = [Float64(value)]
end

function put2d!(ds, name, values, dim)
    v = defVar(ds, name, Float64, (dim, "column"))
    v[:, :] = reshape(Float64.(values), :, 1)
end

function write_snow_dump(inst, tm)
    gstep = tm.nstep + FORTRAN_OFFSET
    (GLOBAL_LO <= gstep <= GLOBAL_HI) || return nothing

    nlevsno = CLM.varpar.nlevsno
    nlevgrnd = CLM.varpar.nlevgrnd
    nlevtot = nlevsno + nlevgrnd

    col = inst.column
    temp = inst.temperature
    ws = inst.water.waterstatebulk_inst.ws
    wd = inst.water.waterdiagnosticbulk_inst
    wfb = inst.water.waterfluxbulk_inst
    wf = wfb.wf
    ef = inst.energyflux

    f = joinpath(OUTDIR, @sprintf("julia_snow_step_n%d.nc", gstep))
    ds = NCDataset(f, "c")
    defDim(ds, "column", 1)
    defDim(ds, "levtot", nlevtot)
    defDim(ds, "levsno", nlevsno)

    put2d!(ds, "T_SOISNO", temp.t_soisno_col[1, 1:nlevtot], "levtot")
    put2d!(ds, "H2OSOI_LIQ", ws.h2osoi_liq_col[1, 1:nlevtot], "levtot")
    put2d!(ds, "H2OSOI_ICE", ws.h2osoi_ice_col[1, 1:nlevtot], "levtot")
    put2d!(ds, "DZSNO", col.dz[1, 1:nlevsno], "levsno")
    put2d!(ds, "QFLX_SNOMELT_LYR", wfb.qflx_snomelt_lyr_col[1, 1:nlevsno], "levsno")
    put2d!(ds, "QFLX_SNOFRZ_LYR", wf.qflx_snofrz_lyr_col[1, 1:nlevsno], "levsno")
    put2d!(ds, "QFLX_SNOW_PERC", wf.qflx_snow_percolation_col[1, 1:nlevsno], "levsno")

    put1d!(ds, "SNLSNO", col.snl[1])
    put1d!(ds, "T_GRND", temp.t_grnd_col[1])
    put1d!(ds, "H2OSFC", ws.h2osfc_col[1])
    put1d!(ds, "INT_SNOW", inst.water.waterstatebulk_inst.int_snow_col[1])
    put1d!(ds, "H2OSNO_NO_LAYERS", ws.h2osno_no_layers_col[1])
    put1d!(ds, "SNOW_DEPTH", wd.snow_depth_col[1])
    put1d!(ds, "frac_sno", wd.frac_sno_col[1])
    put1d!(ds, "frac_sno_eff", wd.frac_sno_eff_col[1])
    put1d!(ds, "SNOMELT_ACCUM", wd.snomelt_accum_col[1])
    put1d!(ds, "QFLX_SNOMELT", wf.qflx_snomelt_col[1])
    put1d!(ds, "QFLX_SNOFRZ", wf.qflx_snofrz_col[1])
    put1d!(ds, "QFLX_SNOW_DRAIN", wf.qflx_snow_drain_col[1])
    put1d!(ds, "QFLX_RAIN_PLUS_SNOMELT", wf.qflx_rain_plus_snomelt_col[1])
    put1d!(ds, "QFLX_LIQ_GRND", wf.qflx_liq_grnd_col[1])
    put1d!(ds, "QFLX_SNOW_GRND", wf.qflx_snow_grnd_col[1])
    put1d!(ds, "QFLX_SOLIDDEW_TOP", wf.qflx_soliddew_to_top_layer_col[1])
    put1d!(ds, "QFLX_LIQDEW_TOP", wf.qflx_liqdew_to_top_layer_col[1])
    put1d!(ds, "QFLX_SOLIDEVAP_TOP", wf.qflx_solidevap_from_top_layer_col[1])
    put1d!(ds, "QFLX_LIQEVAP_TOP", wf.qflx_liqevap_from_top_layer_col[1])
    put1d!(ds, "EFLX_SNOMELT", ef.eflx_snomelt_col[1])

    defVar(ds, "nstep", Int32, ())[] = Int32(gstep)
    defVar(ds, "julia_nstep", Int32, ())[] = Int32(tm.nstep)
    ymd = year(tm.current_date) * 10000 + month(tm.current_date) * 100 + day(tm.current_date)
    tod = hour(tm.current_date) * 3600 + minute(tm.current_date) * 60 + second(tm.current_date)
    defVar(ds, "timemgr_rst_curr_ymd", Int32, ())[] = Int32(ymd)
    defVar(ds, "timemgr_rst_curr_tod", Int32, ())[] = Int32(tod)
    close(ds)
    return nothing
end

println("="^64)
println("  Maritime HJ Andrews Julia snow per-step dump")
println("  global nstep $GLOBAL_LO:$GLOBAL_HI -> local $LOCAL_LO:$LOCAL_HI")
println("  outdir: $OUTDIR")
println("="^64)

CLM.rooting_profile_config.rooting_profile_method_water = CLM.JACKSON_1996_ROOT
CLM.rooting_profile_config.rooting_profile_method_carbon = CLM.JACKSON_1996_ROOT
CLM.snow_hydrology_set_control_for_testing!(;
    wind_dep_snow_density = true,
    overburden_compaction_method = CLM.OVERBURDEN_COMPACTION_VIONNET2012)

t0 = time()
run_clm!(;
    fsurdat = FSURDAT,
    paramfile = FPARAM,
    fforcing = FFORCING,
    fhistory = joinpath(OUTDIR, "julia_snow_probe_history.nc"),
    start_date = DateTime(2017, 1, 1),
    end_date = DateTime(2017, 1, 1) + Hour(LOCAL_HI),
    dtime = 3600,
    use_cn = false,
    verbose = true,
    use_aquifer_layer = false,
    use_hydrstress = true,
    use_luna = true,
    h2osfcflag = 1,
    baseflow_scalar = 0.001,
    int_snow_max = 2000.0,
    ffortran_restart = FRESTART,
    interp_forcing = true,
    fsnowoptics = isfile(_snicar("snicar_optics_5bnd_c013122.nc")) ? _snicar("snicar_optics_5bnd_c013122.nc") : "",
    fsnowaging = isfile(_snicar("snicar_drdt_bst_fit_60_c070416.nc")) ? _snicar("snicar_drdt_bst_fit_60_c070416.nc") : "",
    step_probe = write_snow_dump)

@printf("Done in %.1f s. Dumps in %s\n", time() - t0, OUTDIR)
