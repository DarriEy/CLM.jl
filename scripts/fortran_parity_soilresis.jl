# ==========================================================================
# fortran_parity_soilresis.jl — bare-soil evaporation resistance parity.
#
# Tier-0 "free win" from docs/FULL_PLATFORM_COVERAGE_PLAN.md: the dump carries
# `SOILRESIS` (col, S&L14 soil evaporative resistance, s/m) but no harness
# diffs it. calc_soilevap_resis! runs automatically inside clm_drv!
# (clm_driver.jl:832), so we let Julia COMPUTE it from the injected IC and diff
# the post-step value against the after_hydrologydrainage dump.
#
# NOTE on masking: read_fortran_restart! DOES inject SOILRESIS (the before_step
# value) — but the step OVERWRITES it (calc_soilevap_resis! recomputes from the
# evolved dsl / t_soisno / vapor diffusivity). The BGC summer window confirms
# SOILRESIS changes every step (e.g. 2712.6 -> 2671.6 over n1757852), so the
# post-step Julia value is a genuine recompute, not the injected number.
#
#   julia +1.12 --project=. scripts/fortran_parity_soilresis.jl
# ==========================================================================
include(joinpath(@__DIR__, "fortran_parity_common.jl"))

const BGC_DUMPDIR = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_summer"
const NSTEP = 1757852

function main()
    bdump = joinpath(BGC_DUMPDIR, "pdump_before_step_n$(NSTEP).nc")
    edump = joinpath(BGC_DUMPDIR, "pdump_after_hydrologydrainage_n$(NSTEP).nc")
    (isfile(bdump) && isfile(edump)) || (println("dumps missing"); return 1)

    # injected SOILRESIS (before_step) — recorded so we prove the step overwrote it
    ds0 = NCDataset(bdump, "r"); inj = Float64(ds0["SOILRESIS"][1]); close(ds0)

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

    jl = Float64(inst.soilstate.soilresis_col[1])
    ds = NCDataset(edump, "r"); f = Float64(ds["SOILRESIS"][1]); close(ds)

    rel = abs(jl - f) / (1.0 + max(abs(jl), abs(f)))
    relpct = abs(jl - f) / abs(f)
    @printf("  injected (before_step) SOILRESIS = %.4f\n", inj)
    @printf("  Julia recomputed differs from injected? %s\n",
            abs(jl - inj) > 1e-6 ? "YES (genuine recompute)" : "NO — WARNING masking")
    @printf("  %-12s %14s %14s\n", "field", "Julia", "Fortran")
    @printf("  %-12s %14.6f %14.6f\n", "SOILRESIS", jl, f)
    @printf("  max|abs| = %.3e   rel = %.3e   (%.4f%% of |Fortran|)\n", abs(jl - f), rel, 100 * relpct)
    return 0
end

exit(main())
