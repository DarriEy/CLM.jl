using NCDatasets, Dates, Statistics, CLM
const run_clm! = getfield(CLM, Symbol("clm_run!"))
const basedir = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped"

fh = tempname() * "_runoff.nc"
# Use calibrated parameters (from DDS optimization) for parity with Fortran reference
const caldir = joinpath(basedir, "optimization/CLM/dds_run_1/final_evaluation/settings/CLM/parameters")
inst = run_clm!(;
    fsurdat=joinpath(caldir, "surfdata_clm.nc"),
    paramfile=joinpath(caldir, "clm5_params.nc"),
    fforcing=joinpath(basedir, "data/forcing/CLM_input/clmforc.2002.nc"),
    fhistory=fh,
    start_date=DateTime(2002,1,1), end_date=DateTime(2003,1,1),
    fsnowoptics="/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_optics_5bnd_c013122.nc",
    fsnowaging="/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_drdt_bst_fit_60_c070416.nc",
    use_aquifer_layer=false,
    baseflow_scalar=0.0022119554,
    int_snow_max=3113.2227,
    verbose=false)

# Read Julia output — annual summary
ds_j = NCDataset(fh)
println("=== Julia annual summary ===")
for v in ["QRUNOFF","QOVER","QFLX_DRAIN","QFLX_INFL","QFLX_TRAN_VEG","H2OSNO","FSAT",
          "QFLX_RAIN_PLUS_SNOMELT","QFLX_SNOMELT"]
    if haskey(ds_j, v)
        d = ds_j[v][:] |> skipmissing |> collect
        println("  $v: min=$(round(minimum(d),sigdigits=4)), max=$(round(maximum(d),sigdigits=4)), mean=$(round(mean(d),sigdigits=4))")
    end
end

# Monthly comparison
ndays_per_month = [31,28,31,30,31,30,31,31,30,31,30,31]
for vname in ["QRUNOFF", "QOVER", "QFLX_DRAIN"]
    if haskey(ds_j, vname)
        vd = ds_j[vname][:] |> skipmissing |> collect
        local idx = 1
        println("\nMonthly Julia $vname (mm/day):")
        for (m, nd) in enumerate(ndays_per_month)
            if idx+nd-1 <= length(vd)
                chunk = vd[idx:idx+nd-1]
                println("  $(lpad(m,2)): $(round(mean(chunk)*86400, sigdigits=4))")
                idx += nd
            end
        end
    end
end
close(ds_j)

# Fortran reference
f_h0 = joinpath(basedir, "optimization/CLM/dds_run_1/final_evaluation/Bow_at_Banff_lumped.clm2.h0.2002-01-01-00000.nc")
ds_f = NCDataset(f_h0)
println("\n=== Fortran annual summary ===")
for v in ["QRUNOFF","QOVER","QDRAI","QINFL"]
    if haskey(ds_f, v)
        d = ds_f[v][:] |> skipmissing |> collect
        println("  $v: min=$(round(minimum(d),sigdigits=4)), max=$(round(maximum(d),sigdigits=4)), mean=$(round(mean(d),sigdigits=4))")
    end
end

for vf in ["QRUNOFF","QOVER","QDRAI"]
    if haskey(ds_f, vf)
        d = ds_f[vf][:] |> skipmissing |> collect
        local idx = 1
        println("\nMonthly Fortran $vf (mm/day):")
        for (m, nd) in enumerate(ndays_per_month)
            if idx+nd-1 <= length(d)
                chunk = d[idx:idx+nd-1]
                println("  $(lpad(m,2)): $(round(mean(chunk)*86400, sigdigits=4))")
                idx += nd
            end
        end
    end
end
close(ds_f)

# Final instantaneous state
wf = inst.water.waterfluxbulk_inst.wf
ws = inst.water.waterstatebulk_inst.ws
sh = inst.soilhydrology
println("\n=== Julia final state ===")
println("  zwt = $(sh.zwt_col[1])")
println("  wa = $(ws.wa_col[1])")
println("  qflx_runoff = $(wf.qflx_runoff_col[1])")
println("  qflx_surf = $(wf.qflx_surf_col[1])")
println("  qflx_drain = $(wf.qflx_drain_col[1])")
println("  qflx_drain_perched = $(wf.qflx_drain_perched_col[1])")
println("  qflx_qrgwl = $(wf.qflx_qrgwl_col[1])")
println("  qflx_infl = $(wf.qflx_infl_col[1])")

rm(fh, force=true)
