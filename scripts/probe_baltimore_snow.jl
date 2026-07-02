# probe_baltimore_snow.jl — dump per-step snow-column state for Baltimore to
# localize the thin-warm-snow density (+12.6% depth) residual. Records mass-
# weighted snow temperature, liquid fraction, bulk density, layer count, and the
# compaction inputs (frac_sno, imelt, frac_iceold, forcing T/wind).
#   julia +1.12 --project=. scripts/probe_baltimore_snow.jl
# Output: paper/data/probe_baltimore_snow.csv

using NCDatasets, Dates, Printf, CLM
const run_clm! = getfield(CLM, Symbol("clm_run!"))
const DATA = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data"
_snicar(f) = let d1 = "$DATA/installs/cesm-inputdata/lnd/clm2/snicardata",
                 d2 = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata"
    isfile(joinpath(d1, f)) ? joinpath(d1, f) : joinpath(d2, f)
end
const CAL = "$DATA/domain_Urban_DeadRun_Baltimore/settings/CLM/parameters"
cfg = (year=2013, dtime=3600, baseflow=0.001, int_snow=2000.0,
    forcing="$DATA/domain_Urban_DeadRun_Baltimore/data/forcing/CLM_input/clmforc.2013.nc",
    restart="$DATA/domain_Urban_DeadRun_Baltimore/simulations/clm_urban/CLM/Urban_DeadRun_Baltimore.clm2.r.2013-01-01-00000.nc")

CLM.rooting_profile_config.rooting_profile_method_water  = CLM.JACKSON_1996_ROOT
CLM.rooting_profile_config.rooting_profile_method_carbon = CLM.JACKSON_1996_ROOT
CLM.snow_hydrology_set_control_for_testing!(; wind_dep_snow_density=true,
    overburden_compaction_method=CLM.OVERBURDEN_COMPACTION_VIONNET2012)

recs = NamedTuple[]
csnow = Ref(0)

function probe(inst, tm)
    ws  = inst.water.waterstatebulk_inst.ws
    wd  = inst.water.waterdiagnosticbulk_inst
    col = inst.column
    temp = inst.temperature
    a2l = inst.atm2lnd
    nsno = CLM.varpar.nlevsno
    # pick the column with the most SWE (the snowy one) each step
    swe_of(c) = ws.h2osno_no_layers_col[c]
    # find max-snow column
    cbest = 0; sbest = -1.0
    for c in eachindex(col.snl)
        s = 0.0
        for j in 1:nsno
            vi = ws.h2osoi_ice_col[c,j]; vl = ws.h2osoi_liq_col[c,j]
            isfinite(vi) && (s += vi); isfinite(vl) && (s += vl)
        end
        s > sbest && (sbest = s; cbest = c)
    end
    cbest == 0 && return
    c = cbest
    snl = Int(col.snl[c])
    snl == 0 && return   # no explicit snow layers
    # mass-weighted snow temp, liquid fraction, total ice+liq, total dz
    txmass=0.0; mtot=0.0; liq=0.0; ice=0.0; dztot=0.0
    for j in (nsno+snl+1):nsno
        vi = ws.h2osoi_ice_col[c,j]; vl = ws.h2osoi_liq_col[c,j]
        m = vi+vl
        txmass += temp.t_soisno_col[c,j]*m; mtot += m
        liq += vl; ice += vi; dztot += col.dz[c,j]
    end
    g = Int(col.gridcell[c])
    fsno = c <= length(wd.frac_sno_col) ? wd.frac_sno_col[c] : NaN
    # top snow layer (topmost = index nsno+snl+1) compaction inputs
    jt = nsno + snl + 1
    ficeold = wd.frac_iceold_col[c, jt]
    imelt = inst.temperature.imelt_col[c, jt]
    wx_t = ws.h2osoi_ice_col[c,jt] + ws.h2osoi_liq_col[c,jt]
    fi_t = wx_t > 0 ? ws.h2osoi_ice_col[c,jt]/wx_t : NaN
    bi_t = ws.h2osoi_ice_col[c,jt] / (fsno * col.dz[c,jt])
    push!(recs, (date=tm.current_date, snl=snl,
        swe = ice+liq,
        snowdp = c <= length(wd.snow_depth_col) ? wd.snow_depth_col[c] : NaN,
        rho = dztot>0 ? (ice+liq)/dztot : NaN,
        tsno = mtot>0 ? txmass/mtot - 273.15 : NaN,
        liqfrac = (ice+liq)>0 ? liq/(ice+liq) : NaN,
        fsno = fsno,
        forc_t = a2l.forc_t_downscaled_col[c]-273.15,
        forc_wind = g<=length(a2l.forc_wind_grc) ? a2l.forc_wind_grc[g] : NaN,
        # top-layer compaction inputs
        t_top = inst.temperature.t_soisno_col[c,jt]-273.15,
        bi_top = bi_t, ficeold=ficeold, fi=fi_t, imelt=Int(imelt),
        liq_top = ws.h2osoi_liq_col[c,jt], dz_top = col.dz[c,jt]))
end

println("Baltimore snow probe")
run_clm!(; fsurdat=joinpath(CAL,"surfdata_clm.nc"), paramfile=joinpath(CAL,"clm5_params.nc"),
    fforcing=cfg.forcing, fhistory=joinpath(@__DIR__,"..","paper","data","probe_balt_snow_hist.nc"),
    start_date=DateTime(cfg.year,1,1), end_date=DateTime(cfg.year+1,1,1),
    dtime=cfg.dtime, use_cn=false, verbose=false, use_aquifer_layer=false,
    use_hydrstress=true, use_luna=true, h2osfcflag=1,
    baseflow_scalar=cfg.baseflow, int_snow_max=cfg.int_snow,
    ffortran_restart=cfg.restart, interp_forcing=true,
    fsnowoptics=isfile(_snicar("snicar_optics_5bnd_c013122.nc")) ? _snicar("snicar_optics_5bnd_c013122.nc") : "",
    fsnowaging=isfile(_snicar("snicar_drdt_bst_fit_60_c070416.nc")) ? _snicar("snicar_drdt_bst_fit_60_c070416.nc") : "",
    step_probe=probe)

csv = joinpath(@__DIR__,"..","paper","data","probe_baltimore_snow.csv")
open(csv, "w") do io
    println(io, "date,snl,swe,snowdp,rho,tsno,liqfrac,fsno,forc_t,forc_wind,t_top,bi_top,ficeold,fi,imelt,liq_top,dz_top")
    for r in recs
        @printf(io, "%s,%d,%.5g,%.5g,%.5g,%.5g,%.5g,%.5g,%.5g,%.5g,%.5g,%.5g,%.5g,%.5g,%d,%.5g,%.5g\n",
            r.date, r.snl, r.swe, r.snowdp, r.rho, r.tsno, r.liqfrac, r.fsno, r.forc_t, r.forc_wind,
            r.t_top, r.bi_top, r.ficeold, r.fi, r.imelt, r.liq_top, r.dz_top)
    end
end
println("-> $csv ($(length(recs)) snow steps)")
