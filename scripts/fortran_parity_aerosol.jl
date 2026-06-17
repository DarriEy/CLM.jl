# ==========================================================================
# fortran_parity_aerosol.jl — snow aerosol mass parity.
#
# Tier-0 "free win" from docs/FULL_PLATFORM_COVERAGE_PLAN.md: the dumps carry
# the snow aerosol masses mss_bcphi/mss_bcpho/mss_dst1..4/mss_ocphi/mss_ocpho
# (col x levsno, kg) but no harness diffs them. aerosol_masses! runs
# automatically inside clm_drv! (clm_driver.jl:1324/1342), so we let Julia
# COMPUTE the masses from the injected IC and diff the post-step values against
# the after_hydrologydrainage dump.
#
# IMPORTANT CAVEAT: the only available pdump_before_step dumps are the BGC
# SUMMER window (n1757845..n1757897), which is entirely SNOW-FREE — every
# mss_* field is 0.0 in every dump. read_fortran_restart! DOES inject the
# mss_* fields (all 0.0); aerosol_masses! then recomputes them, and with no
# snow the answer is 0.0. So this probe verifies that Julia's aerosol_masses!
# produces the same (zero) masses as Fortran on a snow-free column — a
# structural / no-crash / mass-conservation sanity check, NOT a stress test of
# the aerosol redistribution logic (which would need a snowy dump that does not
# exist on this machine). We report it honestly as such.
#
#   julia +1.12 --project=. scripts/fortran_parity_aerosol.jl
# ==========================================================================
include(joinpath(@__DIR__, "fortran_parity_common.jl"))

const BGC_DUMPDIR = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_summer"
const NSTEP = 1757852

const AER_FIELDS = [
    ("mss_bcphi", :mss_bcphi_col), ("mss_bcpho", :mss_bcpho_col),
    ("mss_ocphi", :mss_ocphi_col), ("mss_ocpho", :mss_ocpho_col),
    ("mss_dst1",  :mss_dst1_col),  ("mss_dst2",  :mss_dst2_col),
    ("mss_dst3",  :mss_dst3_col),  ("mss_dst4",  :mss_dst4_col),
]

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

    nlevsno = CLM.varpar.nlevsno
    ds = NCDataset(edump, "r")
    @printf("  %-12s %12s %12s %6s  %s\n", "field", "max|abs|", "max|rel|", "npts", "status")
    @printf("  %s\n", "-"^56)
    gmax = 0.0; gabs = 0.0; fdumptot = 0.0
    rows = Tuple{String,Float64,Float64,Int}[]
    for (fname, jfield) in AER_FIELDS
        haskey(ds, fname) || continue
        hasproperty(inst.aerosol, jfield) || continue
        f = ds[fname][:, 1]                       # (levsno_dump,)
        jl = getfield(inst.aerosol, jfield)        # (col, nlevsno)
        nlev = min(nlevsno, length(f), size(jl, 2))
        mabs = 0.0; mrel = 0.0
        for k in 1:nlev
            fv = ismissing(f[k]) ? NaN : Float64(f[k])
            jv = Float64(jl[1, k])
            (isnan(fv) && isnan(jv)) && continue
            fdumptot += abs(fv)
            a = abs(fv - jv); mabs = max(mabs, a)
            mrel = max(mrel, a / (1.0 + max(abs(fv), abs(jv))))
        end
        gmax = max(gmax, mrel); gabs = max(gabs, mabs)
        push!(rows, (fname, mabs, mrel, nlev))
    end
    close(ds)
    for (n, a, r, np) in rows
        @printf("  %-12s %12.3e %12.3e %6d  %s\n", n, a, r, np, r > 1e-9 ? "DIFF" : "ok")
    end
    @printf("  %s\n", "-"^56)
    @printf("  aerosol-mass global max|rel| = %.3e  max|abs| = %.3e\n", gmax, gabs)
    @printf("  (NOTE: dump total |mss| = %.3e — SNOW-FREE summer window, all masses 0;\n", fdumptot)
    @printf("   this is a structural / no-crash parity check, not an aerosol-physics stress test)\n")
    return 0
end

exit(main())
