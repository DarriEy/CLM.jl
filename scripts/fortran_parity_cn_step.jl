# ==========================================================================
# fortran_parity_cn_step.jl — CN/BGC single-step parity (Track 2).
# Inject the before_step BGC dump as the shared IC, run ONE use_cn clm_drv!
# step, and diff the post-step CN pools against the end-of-step Fortran dump
# (after_hydrologydrainage), which captures the post-CNDriver pools.
#
#   julia +1.12 --project=. scripts/fortran_parity_cn_step.jl
# ==========================================================================
include(joinpath(@__DIR__, "fortran_parity_common.jl"))

const BGC_DUMPDIR = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_dumps_2202"
const NSTEP = 1753153

function main()
    bdump = joinpath(BGC_DUMPDIR, "pdump_before_step_n$(NSTEP).nc")
    edump = joinpath(BGC_DUMPDIR, "pdump_after_hydrologydrainage_n$(NSTEP).nc")
    (isfile(bdump) && isfile(edump)) || (println("dumps missing"); return 1)

    println("Running ONE use_cn=true clm_drv! step from BGC IC (nstep=$NSTEP) ...")
    local inst, bounds
    try
        inst, bounds = run_one_parity_step!(NSTEP; use_cn=true, dumpdir=BGC_DUMPDIR,
                                            step_date=DateTime(2003,1,1))
    catch e
        println("\n!!! use_cn clm_drv! step CRASHED:")
        Base.showerror(stdout, e, catch_backtrace())
        println()
        return 2
    end
    println("step completed without error.\n")

    ccs = inst.bgc_vegetation.cnveg_carbonstate_inst
    cns = inst.bgc_vegetation.cnveg_nitrogenstate_inst
    scs = inst.soilbiogeochem_carbonstate
    sns = inst.soilbiogeochem_nitrogenstate

    ds = NCDataset(edump, "r")
    jl_pfts = Int.(inst.patch.itype); f_pfts = Int.(ds["pfts1d_itypveg"][:])
    jl2f = Dict{Int,Int}()
    for pj in 1:length(jl_pfts)
        for pf in 1:length(f_pfts)
            if f_pfts[pf] == jl_pfts[pj]; jl2f[pj] = pf; break; end
        end
    end

    gmax = 0.0; rows = Tuple{String,Float64,Float64}[]
    function patchcmp(name, jlarr)
        d = _dumpvar(ds, name); d === nothing && return
        mabs=0.0; mrel=0.0
        for pj in 1:min(length(jlarr), length(f_pfts))
            haskey(jl2f,pj) || continue
            fv=d[jl2f[pj]]; jv=Float64(jlarr[pj])
            (isnan(fv)&&isnan(jv)) && continue
            a=abs(fv-jv); mabs=max(mabs,a); mrel=max(mrel, a/(1+max(abs(fv),abs(jv))))
        end
        gmax=max(gmax,mrel); push!(rows,(name,mabs,mrel))
    end
    function poolcmp(name, arr3, pidx)
        d=_dumpvar(ds,name); d===nothing && return
        nlev=min(length(d),size(arr3,2)); mabs=0.0; mrel=0.0
        for j in 1:nlev
            fv=d[j]; jv=Float64(arr3[1,j,pidx]); (isnan(fv)&&isnan(jv)) && continue
            a=abs(fv-jv); mabs=max(mabs,a); mrel=max(mrel,a/(1+max(abs(fv),abs(jv))))
        end
        gmax=max(gmax,mrel); push!(rows,(name,mabs,mrel))
    end
    function col2cmp(name, arr)
        d=_dumpvar(ds,name); d===nothing && return
        nlev=min(length(d),size(arr,2)); mabs=0.0; mrel=0.0
        for j in 1:nlev
            fv=d[j]; jv=Float64(arr[1,j]); (isnan(fv)&&isnan(jv)) && continue
            a=abs(fv-jv); mabs=max(mabs,a); mrel=max(mrel,a/(1+max(abs(fv),abs(jv))))
        end
        gmax=max(gmax,mrel); push!(rows,(name,mabs,mrel))
    end

    for s in ("leafc","frootc","livestemc","deadstemc","livecrootc","deadcrootc",
              "cpool","xsmrpool","gresp_storage")
        f=Symbol(s*"_patch"); hasproperty(ccs,f) && patchcmp(s, getfield(ccs,f))
    end
    for s in ("leafn","frootn","deadstemn","retransn","npool")
        f=Symbol(s*"_patch"); hasproperty(cns,f) && patchcmp(s, getfield(cns,f))
    end
    for (vn,p) in ("litr1c_vr"=>1,"litr2c_vr"=>2,"litr3c_vr"=>3,"soil1c_vr"=>4,
                   "soil2c_vr"=>5,"soil3c_vr"=>6,"cwdc_vr"=>7)
        poolcmp(vn, scs.decomp_cpools_vr_col, p)
    end
    for (vn,p) in ("litr1n_vr"=>1,"soil1n_vr"=>4,"cwdn_vr"=>7)
        poolcmp(vn, sns.decomp_npools_vr_col, p)
    end
    col2cmp("sminn_vr", sns.sminn_vr_col)
    col2cmp("smin_no3_vr", sns.smin_no3_vr_col)
    col2cmp("smin_nh4_vr", sns.smin_nh4_vr_col)
    close(ds)

    @printf("\n  %-16s %12s %12s\n", "CN field", "max|abs|", "max|rel|")
    @printf("  %s\n", "-"^42)
    for (n,a,r) in sort(rows, by=x->-x[3])
        @printf("  %-16s %12.3e %12.3e %s\n", n, a, r, r>1e-6 ? "DIFF" : "ok")
    end
    @printf("  %s\n  global CN max|rel| = %.3e\n", "-"^42, gmax)
    return 0
end

exit(main())
