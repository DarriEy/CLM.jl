# =============================================================================
# fortran_parity_cn_coldstart.jl — validate the Julia use_cn COLD START against a
# Fortran use_cn cold start.
#
# This was impossible before PR #212: nothing on the live use_cn init path wrote
# leafc / deadstemc / decomp_cpools_vr, so a Julia CN cold start was entirely NaN
# and only a dump injection produced usable state (which is precisely why the
# single-step harness never noticed).
#
# Fortran reference: /Users/.../clm_cn_coldstart/pdump_before_step_n0.nc
#   start_type = startup (cold), use_cn = .true., spinup_state = 0, start_ymd = 20020101.
#   `before_step` at nstep 0 is the state at the top of the FIRST driver call, i.e.
#   immediately after InitCold -- the cold-start state itself.
#
# NOTE ON EXPECTATIONS. A CN cold start is NOT a parity target in the same sense as
# a warm mid-run step: CTSM's cold CN state is a crude seeding (a nominal leafc, zero
# woody pools, a prescribed decomposition-pool profile), and any difference in how the
# two codes seed it is a legitimate finding, not necessarily a bug. What this harness
# asserts is (a) FINITENESS -- no NaN anywhere, the #212 regression guard -- and
# (b) that the seeded values agree with Fortran where both codes claim to seed them.
#
#   julia +1.12 --project=. scripts/fortran_parity_cn_coldstart.jl
# =============================================================================
include(joinpath(@__DIR__, "fortran_parity_common.jl"))

const COLD_DUMP = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_cn_coldstart/pdump_before_step_n0.nc"

_rel(f, j) = abs(f - j) / (1.0 + max(abs(f), abs(j)))

function main()
    isfile(COLD_DUMP) || (println("Fortran cold-start dump missing: $COLD_DUMP"); return 1)

    # A REAL cold start: default COLDSTART_MATCH_FORTRAN, no dump injection.
    # Bow's lnd_in sets &cnvegcarbonstate initial_vegc = 100 (the CTSM default is 20).
    CLM.cnvegcstate_const.initial_vegC = 100.0
    CLM.rooting_profile_config.rooting_profile_method_water  = CLM.JACKSON_1996_ROOT
    CLM.rooting_profile_config.rooting_profile_method_carbon = CLM.JACKSON_1996_ROOT
    (inst, bounds, filt, tm) = CLM.clm_initialize!(;
        fsurdat=FSURDAT, paramfile=FPARAM,
        start_date=DateTime(2002, 1, 1), dtime=3600,
        use_cn=true, use_luna=true, use_bedrock=true, use_aquifer_layer=false,
        h2osfcflag=0, fsnowoptics=FSNOWOPT, fsnowaging=FSNOWAGE,
        int_snow_max=INT_SNOW_MAX)

    ccs = inst.bgc_vegetation.cnveg_carbonstate_inst
    cns = inst.bgc_vegetation.cnveg_nitrogenstate_inst
    cvs = inst.bgc_vegetation.cnveg_state_inst
    scs = inst.soilbiogeochem_carbonstate
    sns = inst.soilbiogeochem_nitrogenstate

    # ---- (a) FINITENESS: the #212 regression guard ------------------------------
    println("="^80)
    println("JULIA use_cn COLD START — finiteness (the #212 guard)")
    println("="^80)
    # Only the SOIL/CROP columns carry soil-BGC state. Bow's column 2 is a deep-lake
    # landunit (lun itype 5): Fortran fills its decomp/mineral-N arrays with spval and
    # Julia with NaN — the same "inactive" convention, and the CN filters exclude it.
    # Checking it would be a false alarm, so restrict to the columns CN actually runs on.
    soilc = [c for c in 1:bounds.endc
             if inst.landunit.itype[inst.column.landunit[c]] in (CLM.ISTSOIL, CLM.ISTCROP)]
    println("  soil-BGC columns: $soilc  (of 1:$(bounds.endc); others are non-soil landunits)")
    nan_report = String[]
    function chk(name, arr)
        A = Array(arr)
        isempty(A) && (push!(nan_report, "$name: EMPTY"); return)
        nn = count(isnan, A)
        nn > 0 && push!(nan_report, "$name: $nn/$(length(A)) NaN")
    end
    function chkc(name, arr)   # column-dimensioned: soil columns only
        A = Array(arr)
        isempty(A) && (push!(nan_report, "$name: EMPTY"); return)
        S = ndims(A) == 3 ? A[soilc, :, :] : A[soilc, :]
        nn = count(isnan, S)
        nn > 0 && push!(nan_report, "$name: $nn/$(length(S)) NaN on soil columns")
    end
    for s in ("leafc","leafc_storage","leafc_xfer","frootc","frootc_storage",
              "livestemc","deadstemc","livecrootc","deadcrootc","cpool","xsmrpool","gresp_storage")
        f = Symbol(s*"_patch"); hasproperty(ccs, f) && chk(s, getfield(ccs, f))
    end
    for s in ("leafn","frootn","livestemn","deadstemn","retransn","npool")
        f = Symbol(s*"_patch"); hasproperty(cns, f) && chk(s, getfield(cns, f))
    end
    chkc("decomp_cpools_vr", scs.decomp_cpools_vr_col)
    chkc("decomp_npools_vr", sns.decomp_npools_vr_col)
    chkc("sminn_vr", sns.sminn_vr_col)
    chkc("smin_no3_vr", sns.smin_no3_vr_col)
    chkc("smin_nh4_vr", sns.smin_nh4_vr_col)
    if isempty(nan_report)
        println("  PASS — every CN state array is fully finite (no NaN).")
    else
        println("  *** NaN / EMPTY found:")
        for r in nan_report; println("    ", r); end
    end

    # ---- (b) values vs the Fortran cold start -----------------------------------
    ds = NCDataset(COLD_DUMP, "r")
    jl_pfts = Int.(inst.patch.itype); jl_cols = Int.(inst.patch.column)
    f_pfts = Int.(ds["pfts1d_itypveg"][:])
    jl2f = Dict{Int,Int}()
    for pj in eachindex(jl_pfts)
        jl_cols[pj] == 1 || continue
        for pf in eachindex(f_pfts)
            if f_pfts[pf] == jl_pfts[pj] && !haskey(jl2f, pj); jl2f[pj] = pf; end
        end
    end

    println("\n" * "="^80)
    println("JULIA vs FORTRAN cold-start CN state (Bow, 2002-01-01, use_cn, spinup_state=0)")
    println("="^80)
    @printf("  %-20s %14s %14s %11s\n", "field (grass pft12)", "Fortran", "Julia", "rel")
    println("  " * "-"^62)

    gj = findfirst(p -> jl_pfts[p] == 12 && jl_cols[p] == 1, eachindex(jl_pfts))
    gf = jl2f[gj]
    worst = 0.0; rows = Tuple{String,Float64,Float64,Float64}[]
    function pcmp(name, jlarr)
        haskey(ds, name) || return
        d = _dumpvar(ds, name)
        fv = d[gf]; jv = Float64(jlarr[gj])
        (isnan(fv) || isnan(jv)) && return
        r = _rel(fv, jv); worst = max(worst, r)
        push!(rows, (name, fv, jv, r))
    end
    for s in ("leafc","leafc_storage","leafc_xfer","frootc","frootc_storage",
              "livestemc","deadstemc","livecrootc","deadcrootc","cpool","xsmrpool","gresp_storage")
        f = Symbol(s*"_patch"); hasproperty(ccs, f) && pcmp(s, getfield(ccs, f))
    end
    for s in ("leafn","frootn","livestemn","deadstemn","retransn","npool")
        f = Symbol(s*"_patch"); hasproperty(cns, f) && pcmp(s, getfield(cns, f))
    end
    for s in ("dormant_flag","offset_flag","onset_flag","annavg_t2m")
        f = Symbol(s*"_patch"); hasproperty(cvs, f) && pcmp(s, getfield(cvs, f))
    end
    for (n, fv, jv, r) in rows
        @printf("  %-20s %14.6g %14.6g %11.3e %s\n", n, fv, jv, r, r > 1e-6 ? "DIFF" : "ok")
    end

    # column-level pools (top layer + column total)
    println("\n  %-20s" |> x -> "")
    @printf("  %-20s %14s %14s %11s\n", "soil pool (lev1)", "Fortran", "Julia", "rel")
    println("  " * "-"^62)
    function ccmp(name, arr3, pidx)
        haskey(ds, name) || return
        d = _dumpvar(ds, name)
        fv = d[1]; jv = Float64(arr3[1, 1, pidx])
        (isnan(fv) || isnan(jv)) && return
        r = _rel(fv, jv); worst = max(worst, r)
        @printf("  %-20s %14.6g %14.6g %11.3e %s\n", name, fv, jv, r, r > 1e-6 ? "DIFF" : "ok")
    end
    for (vn, p) in ("litr1c_vr"=>1, "litr2c_vr"=>2, "litr3c_vr"=>3, "soil1c_vr"=>4,
                    "soil2c_vr"=>5, "soil3c_vr"=>6, "cwdc_vr"=>7)
        ccmp(vn, scs.decomp_cpools_vr_col, p)
    end
    for (vn, p) in ("litr1n_vr"=>1, "soil1n_vr"=>4, "cwdn_vr"=>7)
        ccmp(vn, sns.decomp_npools_vr_col, p)
    end
    function c2cmp(name, arr)
        haskey(ds, name) || return
        d = _dumpvar(ds, name)
        fv = d[1]; jv = Float64(arr[1, 1])
        (isnan(fv) || isnan(jv)) && return
        r = _rel(fv, jv); worst = max(worst, r)
        @printf("  %-20s %14.6g %14.6g %11.3e %s\n", name, fv, jv, r, r > 1e-6 ? "DIFF" : "ok")
    end
    c2cmp("sminn_vr", sns.sminn_vr_col)
    c2cmp("smin_no3_vr", sns.smin_no3_vr_col)
    c2cmp("smin_nh4_vr", sns.smin_nh4_vr_col)
    close(ds)

    @printf("\n  cold-start worst rel = %.3e\n", worst)
    println("\n  Interpretation: a cold CN start is a crude SEEDING, not a converged state;")
    println("  what matters is that it is finite and that the seeded values agree with CTSM.")
    return isempty(nan_report) ? 0 : 2
end

exit(main())
