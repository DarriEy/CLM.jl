# ==========================================================================
# fortran_parity_cn_read.jl — validate that read_fortran_restart! correctly
# injects the CN/BGC pools from a spun-up Fortran BGC dump into a use_cn=true
# Bow instance. This is the FIRST step of CN parity (Track 2): confirm the
# reader plumbing before exercising the use_cn clm_drv! path.
#
#   julia +1.12 --project=. scripts/fortran_parity_cn_read.jl
# ==========================================================================
include(joinpath(@__DIR__, "fortran_parity_common.jl"))

const BGC_DUMPDIR = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_dumps_2202"
const NSTEP = 1753153

function main()
    dumpfile = joinpath(BGC_DUMPDIR, "pdump_before_step_n$(NSTEP).nc")
    isfile(dumpfile) || (println("dump missing: $dumpfile"); return 1)

    println("Building use_cn=true Bow instance ...")
    (inst, bounds, filt, tm) = build_bow_inst(; dtime=3600, use_cn=true,
                                                start_date=DateTime(2003,1,1))
    ccs = inst.bgc_vegetation.cnveg_carbonstate_inst
    cns = inst.bgc_vegetation.cnveg_nitrogenstate_inst
    scs = inst.soilbiogeochem_carbonstate
    sns = inst.soilbiogeochem_nitrogenstate
    @printf("  veg C alloc: leafc_patch len=%d ; soil C alloc: decomp_cpools_vr_col=%s\n",
            length(ccs.leafc_patch), string(size(scs.decomp_cpools_vr_col)))

    println("Injecting BGC dump (read_fortran_restart!) ...")
    inject_dump!(inst, bounds, dumpfile)

    ds = NCDataset(dumpfile, "r")
    jl_pfts = Int.(inst.patch.itype)
    f_pfts  = haskey(ds, "pfts1d_itypveg") ? Int.(ds["pfts1d_itypveg"][:]) : Int[]
    println("\nJulia patch PFT types: ", jl_pfts, "   Fortran: ", f_pfts)

    # map Julia patch -> Fortran patch by PFT type (same logic as set_patch_1d!)
    jl2f = Dict{Int,Int}()
    for pj in 1:length(jl_pfts)
        for pf in 1:length(f_pfts)
            if f_pfts[pf] == jl_pfts[pj]; jl2f[pj] = pf; break; end
        end
    end

    fail = 0
    function check_patch(name, jlarr)
        d = _dumpvar(ds, name); d === nothing && return
        @printf("  %-16s", name)
        # compare col-1 patches only (Julia patch 4 is col-2, not in the 1-col ref)
        for pj in 1:min(length(jlarr), length(f_pfts))
            fv = haskey(jl2f, pj) ? d[jl2f[pj]] : NaN
            jv = Float64(jlarr[pj])
            ok = (isnan(fv) && isnan(jv)) || abs(fv-jv) <= 1e-9*(1+abs(fv))
            ok || (fail += 1)
            @printf("  p%d:%.4g%s", pj, jv, ok ? "" : "(F=$(round(fv,sigdigits=4))!)")
        end
        println()
    end
    function check_pool_vr(name, arr3, pidx)
        d = _dumpvar(ds, name); d === nothing && return
        nlev = min(length(d), size(arr3,2))
        m = 0.0
        for j in 1:nlev
            fv = d[j]; jv = Float64(arr3[1,j,pidx])
            (isnan(fv)&&isnan(jv)) && continue
            m = max(m, abs(fv-jv))
        end
        m > 1e-9 && (fail += 1)
        @printf("  %-16s pool%d max|abs|=%.3e (top3 Jl=%.4g,%.4g,%.4g)\n",
                name, pidx, m, arr3[1,1,pidx], arr3[1,2,pidx], arr3[1,3,pidx])
    end
    function check_col2d(name, arr)
        d = _dumpvar(ds, name); d === nothing && return
        nlev = min(length(d), size(arr,2)); m = 0.0
        for j in 1:nlev
            fv=d[j]; jv=Float64(arr[1,j]); (isnan(fv)&&isnan(jv)) && continue
            m=max(m, abs(fv-jv))
        end
        m > 1e-9 && (fail += 1)
        @printf("  %-16s max|abs|=%.3e (top3 Jl=%.4g,%.4g,%.4g)\n",
                name, m, arr[1,1], arr[1,2], arr[1,3])
    end

    println("\n--- vegetation carbon (Julia value; F=dump if mismatch) ---")
    check_patch("leafc", ccs.leafc_patch)
    check_patch("frootc", ccs.frootc_patch)
    check_patch("livestemc", ccs.livestemc_patch)
    check_patch("deadstemc", ccs.deadstemc_patch)
    check_patch("cpool", ccs.cpool_patch)
    check_patch("xsmrpool", ccs.xsmrpool_patch)
    println("--- vegetation nitrogen ---")
    check_patch("leafn", cns.leafn_patch)
    check_patch("deadstemn", cns.deadstemn_patch)
    check_patch("retransn", cns.retransn_patch)
    check_patch("npool", cns.npool_patch)

    println("--- soil decomposition carbon pools (decomp_cpools_vr_col) ---")
    check_pool_vr("litr1c_vr", scs.decomp_cpools_vr_col, 1)
    check_pool_vr("litr2c_vr", scs.decomp_cpools_vr_col, 2)
    check_pool_vr("litr3c_vr", scs.decomp_cpools_vr_col, 3)
    check_pool_vr("soil1c_vr", scs.decomp_cpools_vr_col, 4)
    check_pool_vr("soil2c_vr", scs.decomp_cpools_vr_col, 5)
    check_pool_vr("soil3c_vr", scs.decomp_cpools_vr_col, 6)
    check_pool_vr("cwdc_vr",   scs.decomp_cpools_vr_col, 7)
    println("--- soil decomposition nitrogen pools (decomp_npools_vr_col) ---")
    check_pool_vr("litr1n_vr", sns.decomp_npools_vr_col, 1)
    check_pool_vr("soil1n_vr", sns.decomp_npools_vr_col, 4)
    check_pool_vr("cwdn_vr",   sns.decomp_npools_vr_col, 7)
    println("--- mineral N ---")
    check_col2d("sminn_vr",    sns.sminn_vr_col)
    check_col2d("smin_no3_vr", sns.smin_no3_vr_col)
    check_col2d("smin_nh4_vr", sns.smin_nh4_vr_col)
    close(ds)

    @printf("\nCN read validation: %s (%d mismatches)\n", fail==0 ? "PASS" : "FAIL", fail)
    return fail == 0 ? 0 : 1
end

exit(main())
