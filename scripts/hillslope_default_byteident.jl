#!/usr/bin/env julia
# DEFAULT (use_hillslope=false) run — dumps the water-balance trajectory so the
# branch output can be diffed byte-for-byte against origin/main. Passes NO
# hillslope keywords, so it runs identically on main and the branch.
using CLM
using Dates, Printf

const FSURDAT = "/private/tmp/claude-501/hillslope_fortran_run/surfdata_aripuana_hillslope.nc"
const PARAMF  = "/Users/darri.eythorsson/projects/scratch/aripuana_data/params/clm5_params.nc"
const FORC    = "/Users/darri.eythorsson/projects/scratch/aripuana_data/forcing/clmforc.2002.nc"
const SNOWOPT = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_optics_5bnd_c013122.nc"
const SNOWAGE = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_drdt_bst_fit_60_c070416.nc"
outfile = length(ARGS) >= 1 ? ARGS[1] : "/tmp/default_run.txt"

nsteps = 24; dtime = 1800
start_date = DateTime(2002,1,1,0,0,0)
end_date   = start_date + Second(nsteps*dtime)

io = open(outfile, "w")
probe = function(inst, tm)
    wb = inst.water.waterbalancebulk_inst
    ws = inst.water.waterstatebulk_inst.ws
    liqsum = sum(ws.h2osoi_liq_col)
    @printf(io, "%.17e %.17e %.17e\n", wb.begwb_grc[1], wb.endwb_grc[1], liqsum)
end

CLM.clm_run!(; fsurdat=FSURDAT, paramfile=PARAMF, fforcing=FORC,
    fhistory=tempname()*".nc", start_date=start_date, end_date=end_date,
    dtime=dtime, use_cn=false, use_bedrock=false, use_aquifer_layer=false,
    fsnowoptics=SNOWOPT, fsnowaging=SNOWAGE, verbose=false, step_probe=probe)
close(io)
println("wrote ", outfile)
