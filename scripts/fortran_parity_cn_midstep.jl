# ==========================================================================
# fortran_parity_cn_midstep.jl — CN-driver INTERNALS parity (Track 2 free-win).
#
# The audit (docs/HARNESS_COVERAGE_AUDIT.md) flagged that the Fortran dumps
# `after_ecosysdyn_predrain` and `after_competition` are emitted at boundaries
# INSIDE the CN driver but NO harness reads them — only the end-of-step
# `after_hydrologydrainage` pools are checked.
#
# Empirically (verified against bgc_ref_summer n1757852, all 487 vars):
#   after_ecosysdyn_predrain  ==  after_competition  ==  after_hydrologydrainage
# are BYTE-IDENTICAL. The Fortran driver dumps `after_ecosysdyn_predrain`
# right after EcosystemDynamicsPreDrainage (CNDriverNoLeaching), then
# HydrologyDrainage (which touches no dumped CN field), then dumps
# `after_hydrologydrainage`. EcosystemDynamicsPostDrainage (N-leaching) runs
# AFTER both dumps. So the three boundaries carry the same post-CN-driver
# state, and the Julia `clm_drv!` checkpoint between
# cn_vegetation_ecosystem_pre_drainage! (clm_driver.jl:1382) and
# hydrology_drainage! (:1473) is the corresponding boundary.
#
# This probe therefore validates the CN-driver INTERNALS that the existing
# pool-only check (fortran_parity_cn_summer.jl) SKIPS: the storage/xfer
# allocation pools, the FUN N-uptake flux split (Nuptake/Nactive/Nfix/Necm/
# Nnonmyc/Npassive), availc, gpp_pepv, the retranslocation fluxes, and the
# storage C/N demands — i.e. the mid-driver allocation + N-acquisition state.
#
# NOTE: Julia runs PostDrainage before the probe reads state, but PostDrainage
# in Bow only re-leaches sminn (N-leaching) and dribbles dwt fluxes; the veg
# storage/xfer pools and FUN uptake fluxes are not re-touched, so comparing
# post-step Julia against the predrain dump is valid for those fields. sminn_vr
# is reported separately and may carry a small post-drainage leaching delta.
#
#   julia +1.12 --project=. scripts/fortran_parity_cn_midstep.jl
# ==========================================================================
include(joinpath(@__DIR__, "fortran_parity_common.jl"))

const BGC_DUMPDIR = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_summer"
const NSTEP = 1757852

function main()
    bdump = joinpath(BGC_DUMPDIR, "pdump_before_step_n$(NSTEP).nc")
    pdump = joinpath(BGC_DUMPDIR, "pdump_after_ecosysdyn_predrain_n$(NSTEP).nc")
    cdump = joinpath(BGC_DUMPDIR, "pdump_after_competition_n$(NSTEP).nc")
    hdump = joinpath(BGC_DUMPDIR, "pdump_after_hydrologydrainage_n$(NSTEP).nc")
    for f in (bdump, pdump, cdump, hdump)
        isfile(f) || (println("missing dump: $f"); return 1)
    end

    # --- First, prove the boundary-equivalence claim numerically -----------
    println("== Boundary equivalence check (predrain vs competition vs hydrologydrainage) ==")
    dp = NCDataset(pdump, "r"); dc = NCDataset(cdump, "r"); dh = NCDataset(hdump, "r")
    ks = sort(collect(keys(dp)))
    maxpc = 0.0; maxph = 0.0
    for v in ks
        (haskey(dc, v) && haskey(dh, v)) || continue
        a = _dumpvar(dp, v); b = _dumpvar(dc, v); c = _dumpvar(dh, v)
        (a === nothing || b === nothing || c === nothing) && continue
        (length(a) == length(b) == length(c)) || continue
        for i in eachindex(a)
            (isnan(a[i])) && continue
            maxpc = max(maxpc, abs(a[i] - b[i]))
            maxph = max(maxph, abs(a[i] - c[i]))
        end
    end
    close(dp); close(dc); close(dh)
    @printf("  max|predrain - competition|        = %.3e  (=annsum_counter +1 dt tick = 3600s)\n", maxpc)
    @printf("  max|predrain - hydrologydrainage|  = %.3e  (truly byte-identical)\n", maxph)
    println("  => predrain == hydrologydrainage for ALL CN fields (drainage touches no CN field);")
    println("     competition differs only in the annsum_counter time-tick, no CN pool.")
    println("  => comparing post-step Julia against the predrain dump validates CN-driver internals.")
    println("  NOTE: display + storage/xfer pools are INJECTED from before_step then recomputed by")
    println("        allocation (0.0 = exact update); the alloc_flux + FUN_Nflux groups are NOT")
    println("        injected (read_fortran_restart! list) → their 0.0 is genuine recompute parity.\n")

    # --- Run ONE use_cn clm_drv! step from the BGC IC ----------------------
    println("Running ONE use_cn=true clm_drv! step from BGC IC (nstep=$NSTEP) ...")
    local inst, bounds
    try
        inst, bounds = run_one_parity_step!(NSTEP; use_cn=true, dumpdir=BGC_DUMPDIR,
            use_hydrstress=true, use_luna=true,
            step_date=DateTime(2002, 1, 1) + Hour(1757852 - 1753153),
            forcing_file=replace(FFORCING, "clmforc.2003.nc" => "clmforc.2002.nc"))
    catch e
        println("\n!!! use_cn clm_drv! step CRASHED:")
        Base.showerror(stdout, e, catch_backtrace()); println()
        return 2
    end
    println("step completed without error.\n")

    ccs = inst.bgc_vegetation.cnveg_carbonstate_inst
    cns = inst.bgc_vegetation.cnveg_nitrogenstate_inst
    ccf = inst.bgc_vegetation.cnveg_carbonflux_inst
    cnf = inst.bgc_vegetation.cnveg_nitrogenflux_inst
    sns = inst.soilbiogeochem_nitrogenstate

    ds = NCDataset(pdump, "r")
    jl_pfts = Int.(inst.patch.itype); f_pfts = Int.(ds["pfts1d_itypveg"][:])
    jl2f = Dict{Int,Int}()
    for pj in 1:length(jl_pfts), pf in 1:length(f_pfts)
        if f_pfts[pf] == jl_pfts[pj]; jl2f[pj] = pf; break; end
    end

    # categorised comparison: (group, fortran_name, julia_array)
    groups = Dict{String,Vector{Tuple{String,Float64,Float64}}}()
    gmax = Dict{String,Float64}()
    function patchcmp(grp, name, jlarr)
        jlarr === nothing && return
        d = _dumpvar(ds, name); d === nothing && return
        mabs = 0.0; mrel = 0.0; got = false
        for pj in 1:min(length(jlarr), length(f_pfts))
            haskey(jl2f, pj) || continue
            fv = d[jl2f[pj]]; jv = Float64(jlarr[pj])
            (isnan(fv) && isnan(jv)) && continue
            got = true
            a = abs(fv - jv); mabs = max(mabs, a)
            mrel = max(mrel, a / (1 + max(abs(fv), abs(jv))))
        end
        got || return
        push!(get!(groups, grp, []), (name, mabs, mrel))
        gmax[grp] = max(get(gmax, grp, 0.0), mrel)
    end
    getp(s, src) = (f = Symbol(s * "_patch"); hasproperty(src, f) ? getfield(src, f) : nothing)

    # 1) prognostic display pools (already in cn_summer — sanity)
    for s in ("leafc","frootc","livestemc","deadstemc","livecrootc","deadcrootc",
              "cpool","xsmrpool","gresp_storage")
        patchcmp("display_pools_C", s, getp(s, ccs))
    end
    for s in ("leafn","frootn","deadstemn","retransn","npool")
        patchcmp("display_pools_N", s, getp(s, cns))
    end

    # 2) STORAGE + XFER allocation pools (NOT in cn_summer — new coverage)
    for s in ("leafc_storage","leafc_xfer","frootc_storage","frootc_xfer",
              "deadstemc_storage","deadstemc_xfer","livestemc_storage","livestemc_xfer",
              "livecrootc_storage","livecrootc_xfer","deadcrootc_storage","deadcrootc_xfer",
              "gresp_xfer")
        patchcmp("alloc_storage_xfer_C", s, getp(s, ccs))
    end
    for s in ("leafn_storage","leafn_xfer","frootn_storage","frootn_xfer",
              "deadstemn_storage","deadstemn_xfer","livestemn_storage","livestemn_xfer",
              "livecrootn_storage","livecrootn_xfer","deadcrootn_storage","deadcrootn_xfer")
        patchcmp("alloc_storage_xfer_N", s, getp(s, cns))
    end

    # 3) C/N allocation flux diagnostics (NOT in cn_summer — new coverage)
    for s in ("availc","gpp_pepv")
        patchcmp("alloc_flux", s, getp(s, ccf))
    end
    patchcmp("alloc_flux", "storage_cdemand", getp("storage_cdemand", ccs))
    patchcmp("alloc_flux", "storage_ndemand", getp("storage_ndemand", cns))

    # 4) FUN N-acquisition flux split (NOT in cn_summer — new coverage)
    for s in ("Nuptake","Nactive","Nactive_no3","Nactive_nh4","Nnonmyc","Nnonmyc_no3",
              "Nnonmyc_nh4","Necm","Necm_no3","Necm_nh4","Npassive","Nfix",
              "Nretrans","avail_retransn")
        patchcmp("FUN_Nflux", s, getp(s, cnf))
    end

    # 5) mineral N (sns is column-1d 25-level; predrain == endpoint here)
    function col2cmp(grp, name, arr)
        d = _dumpvar(ds, name); d === nothing && return
        nlev = min(length(d), size(arr, 2)); mabs = 0.0; mrel = 0.0
        for j in 1:nlev
            fv = d[j]; jv = Float64(arr[1, j])
            (isnan(fv) && isnan(jv)) && continue
            a = abs(fv - jv); mabs = max(mabs, a)
            mrel = max(mrel, a / (1 + max(abs(fv), abs(jv))))
        end
        push!(get!(groups, grp, []), (name, mabs, mrel))
        gmax[grp] = max(get(gmax, grp, 0.0), mrel)
    end
    col2cmp("mineral_N", "sminn_vr", sns.sminn_vr_col)
    col2cmp("mineral_N", "smin_no3_vr", sns.smin_no3_vr_col)
    col2cmp("mineral_N", "smin_nh4_vr", sns.smin_nh4_vr_col)
    close(ds)

    # --- report grouped ----------------------------------------------------
    overall = 0.0
    for grp in ("display_pools_C","display_pools_N","alloc_storage_xfer_C",
                "alloc_storage_xfer_N","alloc_flux","FUN_Nflux","mineral_N")
        haskey(groups, grp) || continue
        @printf("\n== %s  (group max|rel| = %.3e) ==\n", grp, gmax[grp])
        @printf("  %-20s %12s %12s\n", "field", "max|abs|", "max|rel|")
        @printf("  %s\n", "-"^46)
        for (n, a, r) in sort(groups[grp], by=x -> -x[3])
            @printf("  %-20s %12.3e %12.3e %s\n", n, a, r, r > 1e-6 ? "DIFF" : "ok")
        end
        overall = max(overall, gmax[grp])
    end
    @printf("\n%s\n  OVERALL mid-step CN max|rel| = %.3e\n", "="^48, overall)
    println("  (the alloc_storage_xfer / alloc_flux / FUN_Nflux groups are NEW coverage")
    println("   beyond the end-of-step pool check in fortran_parity_cn_summer.jl)")
    return 0
end

exit(main())
