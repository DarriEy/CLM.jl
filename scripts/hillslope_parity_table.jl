#!/usr/bin/env julia
# Full port-vs-Fortran parity table on the 48-step Aripuanã hillslope catena.
using CLM
using Dates, Printf, NCDatasets, Statistics

const FSURDAT = "/private/tmp/claude-501/hillslope_fortran_run/surfdata_aripuana_hillslope.nc"
const HILLF   = "/private/tmp/claude-501/hillslope_aripuana_synthetic_c.nc"
const PARAMF  = "/Users/darri.eythorsson/projects/scratch/aripuana_data/params/clm5_params.nc"
const FORC    = "/Users/darri.eythorsson/projects/scratch/aripuana_data/forcing/clmforc.2002.nc"
const SNOWOPT = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_optics_5bnd_c013122.nc"
const SNOWAGE = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_drdt_bst_fit_60_c070416.nc"
const FDUMP   = "/private/tmp/claude-501/hillslope_fortran_run/hillslope_aripuana_SP_fortran_dump.h0.nc"
const AREA_KM2 = 12259.269

nsteps = 48; dtime = 1800
start_date = DateTime(2002,1,1,0,0,0); end_date = start_date + Second(nsteps*dtime)

rec = NamedTuple[]
probe = function(inst, tm)
    wb = inst.water.waterbalancebulk_inst
    ws = inst.water.waterstatebulk_inst.ws
    wf = inst.water.waterfluxbulk_inst.wf
    sh = inst.soilhydrology
    col = inst.column
    hc = findall(col.is_hillslope_column)
    nsoi = CLM.varpar.nlevsoi
    push!(rec, (
        errg = isempty(wb.errh2o_grc) ? NaN : maximum(abs.(wb.errh2o_grc)),
        errc = isempty(wb.errh2o_col) ? NaN : maximum(abs.(filter(isfinite, wb.errh2o_col))),
        h2ov = copy(ws.h2osoi_vol_col[hc, 1:nsoi]),
        zwt  = copy(sh.zwt_col[hc]),
        swv  = copy(ws.stream_water_volume_lun),
        lat  = copy(wf.qflx_latflow_out_col[hc]),
        qdr  = copy(wf.qflx_drain_col[hc]),
        qdrp = copy(wf.qflx_drain_perched_col[hc]),
        finite = all(isfinite, ws.h2osoi_liq_col) && all(isfinite, sh.zwt_col)))
end

inst = CLM.clm_run!(; fsurdat=FSURDAT, paramfile=PARAMF, fforcing=FORC, fhistory=tempname()*".nc",
    start_date=start_date, end_date=end_date, dtime=dtime, use_cn=false,
    use_hillslope=true, use_hillslope_routing=true, hillslope_file=HILLF,
    gridcell_area_km2=AREA_KM2, use_bedrock=false, use_aquifer_layer=false,
    fsnowoptics=SNOWOPT, fsnowaging=SNOWAGE, verbose=false, step_probe=probe)

fd = NCDataset(FDUMP,"r")
fh2o = Array(fd["H2OSOI"]); fzwt = Array(fd["ZWT"]); fswv = Array(fd["STREAM_WATER_VOLUME"])
flat = Array(fd["QLATFLOWOUT"]); fqdr = Array(fd["QDRAI"]); fqdrp = Array(fd["QDRAI_PERCH"])
close(fd)

# accumulate abs + rel diffs over steps 1..48, Fortran record = k+1
function stats(getport, getfort, ncol=4, nlay=1)
    ad = Float64[]; rd = Float64[]; pmax = 0.0
    for k in 1:48
        r = rec[k]; frec = k+1
        for c in 1:ncol, j in 1:nlay
            pv = nlay==1 ? getport(r)[c] : getport(r)[c,j]
            fv = nlay==1 ? Float64(getfort(frec,c)) : Float64(getfort(frec,c,j))
            push!(ad, abs(pv-fv)); pmax = max(pmax, abs(fv))
            den = max(abs(fv), 1e-12); push!(rd, abs(pv-fv)/den)
        end
    end
    (maxabs=maximum(ad), meanabs=mean(ad), maxrel=maximum(rd), scale=pmax)
end

println("MAX |errh2o_grc| = ", @sprintf("%.3e", maximum(r.errg for r in rec)), " mm")
println("MAX |errh2o_col| = ", @sprintf("%.3e", maximum(r.errc for r in rec if isfinite(r.errc))), " mm")
println("all-finite over 48 steps = ", all(r.finite for r in rec))
println()
@printf("%-22s | %-11s | %-11s | %-11s | %-11s\n","field","maxAbsΔ","meanAbsΔ","maxRelΔ","fort|scale|")
sH  = stats(r->r.h2ov, (f,c,j)->fh2o[c,j,f], 4, 20); @printf("%-22s | %.3e | %.3e | %.3e | %.3e\n","H2OSOI (vol frac)",sH.maxabs,sH.meanabs,sH.maxrel,sH.scale)
sZ  = stats(r->r.zwt, (f,c)->fzwt[c,f]); @printf("%-22s | %.3e | %.3e | %.3e | %.3e\n","ZWT (m)",sZ.maxabs,sZ.meanabs,sZ.maxrel,sZ.scale)
sS  = stats(r->r.swv, (f,c)->fswv[1,f], 1); @printf("%-22s | %.3e | %.3e | %.3e | %.3e\n","STREAM_WATER_VOL (m3)",sS.maxabs,sS.meanabs,sS.maxrel,sS.scale)
sL  = stats(r->r.lat, (f,c)->flat[c,f]); @printf("%-22s | %.3e | %.3e | %.3e | %.3e\n","QLATFLOWOUT (mm/s)",sL.maxabs,sL.meanabs,sL.maxrel,sL.scale)
sD  = stats(r->r.qdr, (f,c)->fqdr[c,f]); @printf("%-22s | %.3e | %.3e | %.3e | %.3e\n","QDRAI (mm/s)",sD.maxabs,sD.meanabs,sD.maxrel,sD.scale)
sDP = stats(r->r.qdrp,(f,c)->fqdrp[c,f]); @printf("%-22s | %.3e | %.3e | %.3e | %.3e\n","QDRAI_PERCH (mm/s)",sDP.maxabs,sDP.meanabs,sDP.maxrel,sDP.scale)
println("\nDiffed against dump records 2..49 (nstep 1..48); record 1 = cold-start IC.")
