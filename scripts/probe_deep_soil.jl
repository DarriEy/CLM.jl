# probe_deep_soil.jl — dump per-timestep deep-layer soil moisture + water table for
# Krycklan, to localize the QDRAI −10% residual (Julia water table ~2.3 cm deeper
# during draining → less baseflow head near the 2 m bedrock). Records h2osoi_vol
# (volumetric, = Fortran H2OSOI) and saturation for soil layers 8–14 (nodes
# 0.80–2.99 m), plus zwt and qflx_drain, for the (first) soil column.
#
#   DOMAIN=Krycklan julia +1.12 --project=. scripts/probe_deep_soil.jl
# Output: paper/data/probe_deepsoil_<domain>.csv

using NCDatasets, Dates, Printf, CLM
const run_clm! = getfield(CLM, Symbol("clm_run!"))
const DATA = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data"
_snicar(f) = let d1 = "$DATA/installs/cesm-inputdata/lnd/clm2/snicardata",
                 d2 = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata"
    isfile(joinpath(d1, f)) ? joinpath(d1, f) : joinpath(d2, f)
end
const DOM = get(ENV, "DOMAIN", "Krycklan")
const CAL = "$DATA/domain_Boreal_Krycklan_Sweden/settings/CLM/parameters"
cfg = (year=2013, dtime=3600, baseflow=0.001, int_snow=2000.0,
    forcing="$DATA/domain_Boreal_Krycklan_Sweden/data/forcing/CLM_input/clmforc.2013.nc",
    restart="$DATA/domain_Boreal_Krycklan_Sweden/simulations/clm_boreal/CLM/Boreal_Krycklan_Sweden.clm2.r.2013-01-01-00000.nc")

CLM.rooting_profile_config.rooting_profile_method_water  = CLM.JACKSON_1996_ROOT
CLM.rooting_profile_config.rooting_profile_method_carbon = CLM.JACKSON_1996_ROOT
CLM.snow_hydrology_set_control_for_testing!(; wind_dep_snow_density=true,
    overburden_compaction_method=CLM.OVERBURDEN_COMPACTION_VIONNET2012)

const LAYERS = 8:14   # soil layers to record (nodes ~0.80–2.99 m)
recs = NamedTuple[]
csoil = Ref(0)

function probe(inst, tm)
    csoil[] == 0 && (csoil[] = 1)
    c = csoil[]
    ws  = inst.water.waterstatebulk_inst.ws
    sh  = inst.soilhydrology
    col = inst.column
    ss  = inst.soilstate
    nsno = CLM.varpar.nlevsno
    DH2O = 1000.0; DICE = 917.0
    vol = Float64[]
    for j in LAYERS
        jj = j + nsno
        v = ws.h2osoi_liq_col[c, jj] / (col.dz[c, jj] * DH2O) +
            ws.h2osoi_ice_col[c, jj] / (col.dz[c, jj] * DICE)
        push!(vol, v)
    end
    push!(recs, (date=tm.current_date,
                 zwt = c <= length(sh.zwt_col) ? sh.zwt_col[c] : NaN,
                 qdrai = inst.water.waterfluxbulk_inst.wf.qflx_drain_col[c],
                 vol = copy(vol)))
end

println("Deep-soil probe: $DOM")
run_clm!(; fsurdat=joinpath(CAL,"surfdata_clm.nc"), paramfile=joinpath(CAL,"clm5_params.nc"),
    fforcing=cfg.forcing, fhistory=joinpath(@__DIR__,"..","paper","data","probe_deep_hist.nc"),
    start_date=DateTime(cfg.year,1,1), end_date=DateTime(cfg.year+1,1,1),
    dtime=cfg.dtime, use_cn=false, verbose=false, use_aquifer_layer=false,
    use_hydrstress=true, use_luna=true, h2osfcflag=1,
    baseflow_scalar=cfg.baseflow, int_snow_max=cfg.int_snow,
    ffortran_restart=cfg.restart, interp_forcing=true,
    fsnowoptics=isfile(_snicar("snicar_optics_5bnd_c013122.nc")) ? _snicar("snicar_optics_5bnd_c013122.nc") : "",
    fsnowaging=isfile(_snicar("snicar_drdt_bst_fit_60_c070416.nc")) ? _snicar("snicar_drdt_bst_fit_60_c070416.nc") : "",
    step_probe=probe)

csv = joinpath(@__DIR__,"..","paper","data","probe_deepsoil_$(lowercase(DOM)).csv")
open(csv, "w") do io
    println(io, "date,zwt,qdrai," * join(["vol$j" for j in LAYERS], ","))
    for r in recs
        @printf(io, "%s,%.6g,%.6g,%s\n", r.date, r.zwt, r.qdrai,
                join((@sprintf("%.6g", v) for v in r.vol), ","))
    end
end
println("-> $csv ($(length(recs)) steps)")
