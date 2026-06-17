# ==========================================================================
# fortran_parity_daylength.jl — daylength (dayl / prev_dayl) probe.
#
# Tier-0 candidate from docs/FULL_PLATFORM_COVERAGE_PLAN.md. init_daylength! /
# update_daylength! compute the gridcell daylength `dayl` and `prev_dayl` from
# the solar declination and latitude (clm_driver.jl:734).
#
# FINDING: the Fortran pdump files do NOT carry `dayl` or `prev_dayl` (verified
# — absent from every pdump_before_step / boundary dump in the BGC summer
# window; only `coszen` is present). So there is NO Fortran ground-truth field
# to diff against on this machine. This probe therefore cannot establish a
# Fortran-dump parity for daylength, and daylength is LEFT OUT of the asserted
# free-wins regression test.
#
# As the best available substitute, this probe verifies that Julia's
# init_daylength! is INTERNALLY consistent with the closed-form astronomical
# daylength dayl = 2 * (SECSPDAY/2pi) * acos(-tan(lat)*tan(decl)) at the dump's
# declination — i.e. that the Julia routine matches the textbook formula the
# Fortran shares. This is a self-consistency check, NOT a Fortran-dump parity.
#
#   julia +1.12 --project=. scripts/fortran_parity_daylength.jl
# ==========================================================================
include(joinpath(@__DIR__, "fortran_parity_common.jl"))

const BGC_DUMPDIR = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_summer"
const NSTEP = 1757852

# Closed-form astronomical daylength (CLM DaylengthMod formula), in seconds.
function _ref_daylength(lat, decl)
    secs_per_rad = CLM.SECSPDAY / (2.0 * π)   # 13750.9871 s/rad
    # clamp latitude away from the poles exactly as DaylengthMod does
    lat_eps = lat
    my_lat  = min(max(lat_eps, -(π/2) + 1e-12), (π/2) - 1e-12)
    temp = -(sin(my_lat) * sin(decl)) / (cos(my_lat) * cos(decl))
    temp = min(1.0, max(-1.0, temp))
    return 2.0 * secs_per_rad * acos(temp)
end

function main()
    bdump = joinpath(BGC_DUMPDIR, "pdump_before_step_n$(NSTEP).nc")
    isfile(bdump) || (println("dump missing"); return 1)

    # Is dayl actually in the dump? Document the absence.
    ds = NCDataset(bdump, "r")
    has_dayl = haskey(ds, "dayl") || haskey(ds, "prev_dayl")
    @printf("  Fortran dump carries dayl/prev_dayl? %s\n", has_dayl ? "YES" : "NO")
    if !has_dayl
        println("  => no Fortran ground-truth daylength field in the dumps; running a")
        println("     SELF-CONSISTENCY check of init_daylength! vs the closed-form formula.\n")
    end
    close(ds)

    step_date = DateTime(2002, 1, 1) + Hour(NSTEP - 1753153)
    forcing   = replace(FFORCING, "clmforc.2003.nc" => "clmforc.2002.nc")
    (inst, bounds, filt, tm) = build_bow_inst(; dtime=3600, start_date=step_date, use_cn=true)

    calday = CLM.get_curr_calday(tm)
    (declin, _)   = CLM.compute_orbital(calday)
    (declinm1, _) = CLM.compute_orbital(calday - 3600.0 / CLM.SECSPDAY)
    obliqr = CLM.ORB_OBLIQR_DEFAULT
    CLM.init_daylength!(inst.gridcell, declin, declinm1, obliqr, 1:bounds.endg)

    grc = inst.gridcell
    gmax = 0.0
    @printf("  %-12s %16s %16s %12s\n", "field", "Julia", "ref formula", "rel")
    @printf("  %s\n", "-"^60)
    for g in 1:bounds.endg
        lat = grc.lat[g]
        jl_dayl = Float64(grc.dayl[g]);     ref_dayl = _ref_daylength(lat, declin)
        jl_prev = Float64(grc.prev_dayl[g]); ref_prev = _ref_daylength(lat, declinm1)
        r1 = abs(jl_dayl - ref_dayl) / (1.0 + max(abs(jl_dayl), abs(ref_dayl)))
        r2 = abs(jl_prev - ref_prev) / (1.0 + max(abs(jl_prev), abs(ref_prev)))
        gmax = max(gmax, r1, r2)
        @printf("  %-12s %16.4f %16.4f %12.3e\n", "dayl[$g]",      jl_dayl, ref_dayl, r1)
        @printf("  %-12s %16.4f %16.4f %12.3e\n", "prev_dayl[$g]", jl_prev, ref_prev, r2)
    end
    @printf("  %s\n", "-"^60)
    @printf("  daylength self-consistency max|rel| = %.3e  (NOT a Fortran-dump parity)\n", gmax)
    return 0
end

exit(main())
