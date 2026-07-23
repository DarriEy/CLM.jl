#!/usr/bin/env julia
# End-to-end hillslope driver run: full clm_run! timestep on the 4-column
# Aripuanã catena (SP + hillslope + routing), matching the Fortran reference
# config in /private/tmp/claude-501/hillslope_fortran_run/, then diffed against
# the Fortran per-column dump.
using CLM
using Dates, Printf, NCDatasets, Statistics

const FSURDAT = "/private/tmp/claude-501/hillslope_fortran_run/surfdata_aripuana_hillslope.nc"
const HILLF   = "/private/tmp/claude-501/hillslope_aripuana_synthetic_c.nc"
const PARAMF  = "/Users/darri.eythorsson/projects/scratch/aripuana_data/params/clm5_params.nc"
const FORC    = "/Users/darri.eythorsson/projects/scratch/aripuana_data/forcing/clmforc.2002.nc"
const SNOWOPT = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_optics_5bnd_c013122.nc"
const SNOWAGE = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_drdt_bst_fit_60_c070416.nc"
const FDUMP   = "/private/tmp/claude-501/hillslope_fortran_run/hillslope_aripuana_SP_fortran_dump.h0.nc"
const AREA_KM2 = 12259.269   # Fortran grc%area (from ESMF mesh)

nsteps = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 48
dtime  = 1800
start_date = DateTime(2002,1,1,0,0,0)
end_date   = start_date + Second(nsteps*dtime)

rec = NamedTuple[]
probe = function(inst, tm)
    wb = inst.water.waterbalancebulk_inst
    ws = inst.water.waterstatebulk_inst.ws
    sh = inst.soilhydrology
    col = inst.column
    hcols = findall(col.is_hillslope_column)
    nsno = CLM.varpar.nlevsno; nsoi = CLM.varpar.nlevsoi
    # h2osoi_vol_col is (nc, nlevgrnd); take soil layers 1:nsoi for hillslope cols
    h2ov = isempty(hcols) ? zeros(0,0) : copy(ws.h2osoi_vol_col[hcols, 1:nsoi])
    zwt  = isempty(hcols) ? Float64[] : copy(sh.zwt_col[hcols])
    errg = isempty(wb.errh2o_grc) ? NaN : maximum(abs.(wb.errh2o_grc))
    swv  = copy(ws.stream_water_volume_lun)
    push!(rec, (errh2o_grc=errg, h2osoi_vol=h2ov, zwt=zwt, stream_vol=copy(swv),
                allfinite = all(isfinite, ws.h2osoi_liq_col) &&
                            all(isfinite, sh.zwt_col) && all(isfinite, swv)))
end

println("=== HILLSLOPE E2E: $nsteps steps (area=$AREA_KM2 km2) ===")
inst = CLM.clm_run!(;
    fsurdat=FSURDAT, paramfile=PARAMF, fforcing=FORC, fhistory=tempname()*".nc",
    start_date=start_date, end_date=end_date, dtime=dtime, use_cn=false,
    use_hillslope=true, use_hillslope_routing=true, hillslope_file=HILLF,
    gridcell_area_km2=AREA_KM2, use_bedrock=false, use_aquifer_layer=false,
    fsnowoptics=SNOWOPT, fsnowaging=SNOWAGE, verbose=false, step_probe=probe)

col = inst.column
hcols = findall(col.is_hillslope_column)
println("\nhillslope columns: ", hcols)
println("cold (downhill): ", col.cold[hcols], "  [Fortran: -9999,1,2,3]")
println("colu (uphill):   ", col.colu[hcols], "  [Fortran: 2,3,4,-9999]")

maxerr = maximum(r.errh2o_grc for r in rec if isfinite(r.errh2o_grc); init=0.0)
allfin = all(r.allfinite for r in rec)
@printf("\nMAX |errh2o_grc| = %.3e mm ; all-finite = %s ; steps=%d\n", maxerr, allfin, length(rec))

# ---- Fortran diff ----
# Dump records: rec 1 = cold-start IC (time 0), recs 2..49 = nstep 1..48.
# Port probe records: rec[k] = state AFTER driver step k (k=1..nsteps).
fd = NCDataset(FDUMP,"r")
fh2o = Array(fd["H2OSOI"])          # (col=4, levsoi=20, time=49)
fzwt = Array(fd["ZWT"])             # (col=4, time=49)
fswv = Array(fd["STREAM_WATER_VOLUME"]) # (lun=1, time=49)
close(fd)

println("\nstep | errh2o(mm) | H2OSOI maxAbsΔ | ZWT maxAbsΔ(m) | port_streamV | fort_streamV")
for k in 1:min(nsteps, 48)
    r = rec[k]
    frec = k + 1                       # Fortran record for nstep k
    # port h2osoi_vol[col, layer] vs fortran fh2o[col, layer, frec]
    dh = 0.0; dz = 0.0
    for c in 1:4, j in 1:20
        dh = max(dh, abs(r.h2osoi_vol[c,j] - Float64(fh2o[c,j,frec])))
    end
    for c in 1:4
        dz = max(dz, abs(r.zwt[c] - Float64(fzwt[c,frec])))
    end
    pv = isempty(r.stream_vol) ? NaN : r.stream_vol[1]
    fv = Float64(fswv[1,frec])
    if k <= 6 || k == nsteps || k % 12 == 0
        @printf("%4d | %10.2e | %13.3e | %13.3e | %11.3e | %11.3e\n",
                k, r.errh2o_grc, dh, dz, pv, fv)
    end
end

# overall H2OSOI/ZWT diff stats over all steps
allh = Float64[]; allz = Float64[]
for k in 1:min(nsteps,48)
    r = rec[k]; frec = k+1
    for c in 1:4, j in 1:20; push!(allh, abs(r.h2osoi_vol[c,j]-Float64(fh2o[c,j,frec]))); end
    for c in 1:4; push!(allz, abs(r.zwt[c]-Float64(fzwt[c,frec]))); end
end
@printf("\nH2OSOI vs Fortran: maxAbsΔ=%.3e  meanAbsΔ=%.3e (vol frac)\n", maximum(allh), mean(allh))
@printf("ZWT    vs Fortran: maxAbsΔ=%.3e  meanAbsΔ=%.3e (m)\n", maximum(allz), mean(allz))
