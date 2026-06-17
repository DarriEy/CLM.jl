# ==========================================================================
# fortran_parity_activelayer.jl — active-layer-thickness (ALT) parity.
#
# Tier-0 "free win" from docs/FULL_PLATFORM_COVERAGE_PLAN.md: the dumps carry
# altmax/altmax_indx/altmax_lastyear/altmax_lastyear_indx (col) but no harness
# diffs them. alt_calc! runs automatically inside clm_drv! (clm_driver.jl:637),
# so we let Julia COMPUTE the annual-max thaw depth from the injected t_soisno
# and diff the post-step values against the after_hydrologydrainage dump.
#
# read_fortran_restart! does NOT inject altmax/altmax_indx (verified — absent
# from the injection list), so altmax_indx (the current-step ALT) is a genuine
# recompute from the injected soil-temperature column. altmax/altmax_lastyear
# are RUNNING annual maxima; with the IC not seeding them, Julia takes the max
# of its cold-start value and this step's ALT. The summer Bow column is fully
# thawed (altmax = 41.998 m, altmax_indx = 25 = bottom soil layer), so the
# current-step ALT pins to the bottom and altmax/indx match exactly once Julia
# also reaches the bottom. We compare altmax_indx (the per-step recompute) as
# the strong field and report altmax/lastyear as context.
#
#   julia +1.12 --project=. scripts/fortran_parity_activelayer.jl
# ==========================================================================
include(joinpath(@__DIR__, "fortran_parity_common.jl"))

const BGC_DUMPDIR = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_summer"
const NSTEP = 1757852

function main()
    bdump = joinpath(BGC_DUMPDIR, "pdump_before_step_n$(NSTEP).nc")
    edump = joinpath(BGC_DUMPDIR, "pdump_after_hydrologydrainage_n$(NSTEP).nc")
    (isfile(bdump) && isfile(edump)) || (println("dumps missing"); return 1)

    println("Running ONE use_cn=true clm_drv! step from BGC IC (nstep=$NSTEP) ...")
    local inst, bounds
    try
        inst, bounds = run_one_parity_step!(NSTEP; use_cn=true, dumpdir=BGC_DUMPDIR,
            use_hydrstress=true, use_luna=true,
            step_date=DateTime(2002, 1, 1) + Hour(NSTEP - 1753153),
            forcing_file=replace(FFORCING, "clmforc.2003.nc" => "clmforc.2002.nc"))
    catch e
        println("\n!!! clm_drv! step CRASHED:")
        Base.showerror(stdout, e, catch_backtrace()); println()
        return 2
    end
    println("step completed.\n")

    al = inst.active_layer
    ds = NCDataset(edump, "r")

    fields = [
        ("altmax",               () -> Float64(al.altmax_col[1])),
        ("altmax_indx",          () -> Float64(al.altmax_indx_col[1])),
        ("altmax_lastyear",      () -> Float64(al.altmax_lastyear_col[1])),
        ("altmax_lastyear_indx", () -> Float64(al.altmax_lastyear_indx_col[1])),
    ]
    @printf("  %-22s %14s %14s %12s\n", "field", "Julia", "Fortran", "rel")
    @printf("  %s\n", "-"^66)
    gmax = 0.0; idx_rel = NaN
    for (fname, getter) in fields
        haskey(ds, fname) || continue
        f = Float64(ds[fname][1]); jl = getter()
        rel = abs(jl - f) / (1.0 + max(abs(jl), abs(f)))
        fname == "altmax_indx" && (idx_rel = rel)
        gmax = max(gmax, rel)
        @printf("  %-22s %14.6f %14.6f %12.3e %s\n", fname, jl, f, rel, rel > 1e-6 ? "DIFF" : "ok")
    end
    close(ds)
    @printf("  %s\n", "-"^66)
    @printf("  active-layer global max|rel| = %.3e   (altmax_indx rel = %.3e, the per-step recompute)\n",
            gmax, idx_rel)
    return 0
end

exit(main())
