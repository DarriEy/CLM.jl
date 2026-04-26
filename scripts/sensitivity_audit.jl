using CLM, NCDatasets, Dates, Statistics, Printf

const BASE_PF = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc"
const BASE_SF = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc"
const FORCING = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/data/forcing/CLM_input/clmforc.2000_2004_spinup.nc"

const DDS = Dict{String,Float64}(
    "baseflow_scalar"=>0.002212,"fff"=>0.10915,"wimp"=>0.01033,
    "ksatdecay"=>0.33006,"n_baseflow"=>2.6479,"e_ice"=>3.4173,
    "perched_baseflow_scalar"=>1.676e-6,"interception_fraction"=>0.6455,
    "max_leaf_wetted_frac"=>0.07266,"fmax"=>0.9809,"bsw_mult"=>1.201,
    "sucsat_mult"=>0.6879,"watsat_mult"=>1.1608,"hksat_mult"=>4.1194,
    "organic_max"=>53.226,"fresh_snw_rds_max"=>70.232,
    "snw_aging_bst"=>90.153,"SNO_Z0MV"=>0.009816,"accum_factor"=>0.001556,
    "SNOW_DENSITY_MAX"=>505.73,"SNOW_DENSITY_MIN"=>141.67,
    "n_melt_coef"=>440.13,"int_snow_max"=>3113.2,
    "medlynslope"=>11.153,"slatop"=>0.00667,"flnr"=>0.07278,
    "froot_leaf"=>2.8143,"stem_leaf"=>2.6721,"route_k"=>31.006,
)

function do_run(pd)
    cp(BASE_PF*".bak", BASE_PF, force=true)
    cp(BASE_SF*".bak", BASE_SF, force=true)
    CLM.apply_all_params!(BASE_PF, BASE_SF, pd)
    ov = CLM.CalibrationOverrides()
    haskey(pd,"baseflow_scalar") && (ov.baseflow_scalar = pd["baseflow_scalar"])
    haskey(pd,"fff") && (ov.fff = pd["fff"])
    haskey(pd,"hksat_mult") && (ov.ksat_scale = pd["hksat_mult"])
    haskey(pd,"medlynslope") && (ov.medlyn_slope = pd["medlynslope"])
    haskey(pd,"bsw_mult") && (ov.bsw_mult = pd["bsw_mult"])
    haskey(pd,"watsat_mult") && (ov.watsat_mult = pd["watsat_mult"])
    haskey(pd,"sucsat_mult") && (ov.sucsat_mult = pd["sucsat_mult"])
    h = tempname()*".nc"
    CLM.clm_run!(;fsurdat=BASE_SF,paramfile=BASE_PF,fforcing=FORCING,fhistory=h,
        start_date=DateTime(2000,1,1),end_date=DateTime(2005,1,1),dtime=1800,verbose=false,
        overrides=ov,int_snow_max=get(pd,"int_snow_max",2000.0))
    ds = NCDataset(h,"r")
    q = haskey(ds,"QRUNOFF") ? Float64.(ds["QRUNOFF"][:]) : Float64[]
    lh = haskey(ds,"EFLX_LH_TOT") ? Float64.(ds["EFLX_LH_TOT"][:]) : Float64[]
    swe = haskey(ds,"H2OSNO") ? Float64.(ds["H2OSNO"][:]) : Float64[]
    close(ds); rm(h, force=true)
    return mean(q[1462:end]), mean(lh[1462:end]), mean(swe[1462:end])
end

function main()
    cp(BASE_PF, BASE_PF*".bak", force=true)
    cp(BASE_SF, BASE_SF*".bak", force=true)

    q0, lh0, swe0 = do_run(DDS)
    println("Baseline: Q=$(round(q0*86400,digits=3))mm/d LH=$(round(lh0,digits=2))W/m² SWE=$(round(swe0,digits=1))mm\n")

    n_active = 0
    n_dead = 0
    @printf("%-26s %8s %8s %9s %9s %9s %s\n","Parameter","TestVal","Default","ΔQ%","ΔLH%","ΔSWE%","Status")
    println("-"^90)

    for pname in sort(collect(keys(DDS)))
        pname == "route_k" && continue
        haskey(CLM.PARAM_BOUNDS, pname) || continue
        lo, hi = CLM.PARAM_BOUNDS[pname]
        test_val = (DDS[pname] > 0.7*hi) ? lo : hi
        d = copy(DDS)
        d[pname] = test_val
        try
            qp, lhp, swep = do_run(d)
            dq = (qp-q0)/(q0+1e-30)*100
            dlh = (lhp-lh0)/(abs(lh0)+1e-30)*100
            dswe = (swep-swe0)/(abs(swe0)+1e-30)*100
            is_act = abs(dq)>0.1 || abs(dlh)>0.5 || abs(dswe)>0.5
            st = is_act ? "ACTIVE" : "DEAD"
            if is_act; n_active += 1; else; n_dead += 1; end
            @printf("%-26s %8.4g %8.4g %+8.2f%% %+8.2f%% %+8.2f%% %s\n",
                pname, test_val, DDS[pname], dq, dlh, dswe, st)
        catch e
            n_dead += 1
            println("  $pname: ERROR — $(sprint(showerror,e)[1:min(80,end)])")
        end
    end

    println("\nTotal: $n_active ACTIVE, $n_dead DEAD out of 28")

    cp(BASE_PF*".bak", BASE_PF, force=true)
    cp(BASE_SF*".bak", BASE_SF, force=true)
    rm(BASE_PF*".bak", force=true)
    rm(BASE_SF*".bak", force=true)
end

main()
