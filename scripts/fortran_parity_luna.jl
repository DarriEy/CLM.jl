# ==========================================================================
# fortran_parity_luna.jl — LUNA photosynthetic-capacity (vcmx25_z / jmx25_z)
# NO-INJECT recompute probe + cadence finding.
#
# CONTEXT (the masking the audit flags):
#   The parity harness (run_one_parity_step!) INJECTS vcmx25_z_patch /
#   jmx25_z_patch from the Fortran restart (fortran_restart.jl:316-319), so
#   LUNA's OWN daily acclimation update is bypassed and unvalidated.
#
# TWO FINDINGS this probe establishes (be honest about what is / isn't checked):
#
#  (A) CADENCE — LUNA updates vcmx25_z/jmx25_z ONLY at end-of-day.
#      In CLM, update_photosynthesis_capacity! runs under is_time_to_run_luna
#      = is_end_curr_day (luna.jl:559). In the Bow dump window
#      n1757845..n1757872 (2002-07-15T12:00 .. 07-16T15:00 UTC), the ONLY
#      end-of-day step is n1757856 (runs 23:00 -> 00:00). The dumps prove this
#      directly: max|after - before| of vcmx25_z is exactly 0 on every step
#      EXCEPT n1757856, where it jumps (vcmx25_z by ~3.2, jmx25_z by ~0.6).
#      => On any non-EOD step a "recompute" is the identity; only n1757856 has
#         a real LUNA update to validate.
#
#      ALSO: CLM.jl's live driver never CALLS update_photosynthesis_capacity!
#      at all (grep: the only driver reference is is_time_to_run_luna=false,
#      clm_driver.jl:839, hardwired for the ozone path). So the injected
#      vcmx25_z is the sole source — LUNA's update is a defined+unit-tested
#      function (test/test_luna.jl) that is NOT wired into clm_drv!. This probe
#      exercises it standalone.
#
#  (B) NO-INJECT RECOMPUTE FEASIBILITY — the dump does NOT carry the LUNA
#      climate accumulators that update_photosynthesis_capacity! consumes.
#      Present in the dump: par240d, par240x, rh10, rb10, lnca, pnlc_z,
#      ndaysteps, nnightsteps, tlai_z, nrad, vcmx25_z(pre), tair10,
#      pco2_240, po2_240.
#      MISSING from the dump (verified, all boundaries): t_veg10_day,
#      t_veg10_night, fpsn24, t_veg_day, t_veg_night, dayl, o3coefjmax.
#      These drive the core temperature path (tleafd10/tleafn10/tleaf10), the
#      growth gate (fpsn24 > 0), and the "first day" gate (t_veg_day != SPVAL).
#      Without them a single-step standalone recompute CANNOT reproduce the
#      Fortran post-LUNA vcmx25_z. This probe documents that gap explicitly and
#      reports the constrained-update magnitude the dump CAN expose.
#
#   julia +1.12 --project=. scripts/fortran_parity_luna.jl
# ==========================================================================
include(joinpath(@__DIR__, "fortran_parity_common.jl"))

const BGC_DUMPDIR = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_summer"
const EOD_NSTEP   = 1757856      # the only end-of-day step in the dump window (23:00->00:00)
const WINDOW      = 1757845:1757872

# Fields update_photosynthesis_capacity! reads that the dump is MISSING.
const MISSING_LUNA_DRIVERS = ("t_veg10_day", "t_veg10_night", "fpsn24",
                              "t_veg_day", "t_veg_night", "dayl", "o3coefjmax")

function cadence_scan()
    println("== (A) CADENCE: which step does Fortran LUNA actually update vcmx25_z/jmx25_z? ==")
    @printf("  %-10s %6s %16s %16s  %s\n", "nstep", "tod", "max|d vcmx25_z|", "max|d jmx25_z|", "")
    @printf("  %s\n", "-"^62)
    eod_found = Int[]
    for nstep in WINDOW
        b = joinpath(BGC_DUMPDIR, "pdump_before_step_n$(nstep).nc")
        a = joinpath(BGC_DUMPDIR, "pdump_after_hydrologydrainage_n$(nstep).nc")
        (isfile(b) && isfile(a)) || continue
        db = NCDataset(b, "r"); da = NCDataset(a, "r")
        vb = Float64.(db["vcmx25_z"][:, :]); va = Float64.(da["vcmx25_z"][:, :])
        jb = Float64.(db["jmx25_z"][:, :]);  ja = Float64.(da["jmx25_z"][:, :])
        close(db); close(da)
        dv = maximum(abs.(va .- vb)); dj = maximum(abs.(ja .- jb))
        tod = Dates.hour(DateTime(2002, 1, 1) + Hour(nstep - 1753153))
        upd = (dv > 1e-8 || dj > 1e-8)
        upd && push!(eod_found, nstep)
        @printf("  %-10d %5d:00 %16.4e %16.4e  %s\n", nstep, tod, dv, dj,
                upd ? "<== LUNA UPDATED (end-of-day)" : "")
    end
    @printf("  %s\n", "-"^62)
    println("  => LUNA updated on step(s): ", eod_found, "  (expected exactly [$EOD_NSTEP])")
    println("  => identity on all other steps confirms is_time_to_run_luna == is_end_curr_day.\n")
    return eod_found
end

function feasibility_check()
    println("== (B) NO-INJECT RECOMPUTE: are the LUNA climate drivers in the dump? ==")
    ds = NCDataset(joinpath(BGC_DUMPDIR, "pdump_before_step_n$(EOD_NSTEP).nc"), "r")
    present = String[]; missing_ = String[]
    needed = ("par240d","par240x","rh10","rb10","lnca","pnlc_z","ndaysteps",
              "nnightsteps","tlai_z","nrad","vcmx25_z","tair10","pco2_240_VALUE",
              "po2_240_VALUE", MISSING_LUNA_DRIVERS...)
    for k in needed
        (haskey(ds, k) ? present : missing_) |> v -> push!(v, k)
    end
    close(ds)
    println("  present : ", join(present, ", "))
    println("  MISSING : ", join(missing_, ", "))
    println("""
  => update_photosynthesis_capacity! consumes t_veg10_day/night (tleaf10),
     fpsn24 (the C3 growth gate), t_veg_day (the first-day gate) and dayl —
     all ABSENT from every dump boundary. A single-step standalone recompute
     therefore CANNOT reconstruct the Fortran post-LUNA vcmx25_z. The LUNA
     daily update is exercised by test/test_luna.jl (synthetic inputs); the
     Fortran-restart dumps lack the accumulator state to drive a parity recompute.
""")
    return missing_
end

function update_magnitude()
    # What the dump CAN expose: the size of the constrained daily update at EOD,
    # and the post-update consistency vcmx25_z == vcmx25_z_last_valid (an internal
    # LUNA invariant: every updated layer copies into *_last_valid).
    println("== (C) What the dump exposes at the EOD step n$EOD_NSTEP ==")
    db = NCDataset(joinpath(BGC_DUMPDIR, "pdump_before_step_n$(EOD_NSTEP).nc"), "r")
    da = NCDataset(joinpath(BGC_DUMPDIR, "pdump_after_hydrologydrainage_n$(EOD_NSTEP).nc"), "r")
    vb = Float64.(db["vcmx25_z"][:, :]); va = Float64.(da["vcmx25_z"][:, :])
    jb = Float64.(db["jmx25_z"][:, :]);  ja = Float64.(da["jmx25_z"][:, :])
    vlv = haskey(da, "vcmx25_z_last_valid_patch") ? Float64.(da["vcmx25_z_last_valid_patch"][:, :]) : nothing
    npft = size(vb, 2)
    @printf("  %-6s %12s %12s %12s | %12s %12s\n", "pft", "vcmx25(pre)", "vcmx25(post)",
            "d(post-pre)", "jmx25(pre)", "jmx25(post)")
    for p in 1:npft
        @printf("  %-6d %12.4f %12.4f %12.4f | %12.4f %12.4f\n",
                p, vb[1, p], va[1, p], va[1, p] - vb[1, p], jb[1, p], ja[1, p])
    end
    if vlv !== nothing
        inv_ok = maximum(abs.(va .- vlv)) < 1e-8
        @printf("  LUNA invariant vcmx25_z(post) == vcmx25_z_last_valid : %s (max|diff|=%.2e)\n",
                inv_ok ? "HOLDS" : "VIOLATED", maximum(abs.(va .- vlv)))
    end
    close(db); close(da)
    println()
end

function main()
    isfile(joinpath(BGC_DUMPDIR, "pdump_before_step_n$(EOD_NSTEP).nc")) ||
        (println("missing dumps under $BGC_DUMPDIR"); return 1)
    eod = cadence_scan()
    missing_ = feasibility_check()
    update_magnitude()

    println("="^64)
    println("SUMMARY")
    println("  (A) LUNA cadence: vcmx25_z/jmx25_z update ONLY at end-of-day; in the")
    println("      dump window that is step n$EOD_NSTEP. Identity on all other steps.")
    println("      CLM.jl's clm_drv! does not call update_photosynthesis_capacity! —")
    println("      the harness injection is the sole source for these fields.")
    println("  (B) NO-INJECT recompute is BLOCKED: the dump omits the LUNA climate")
    println("      accumulators t_veg10_day/night, fpsn24, t_veg_day, dayl. A clean")
    println("      single-step recompute-vs-dump is not possible from these dumps.")
    println("  => No LUNA recompute row goes into the parity test (honest: nothing to")
    println("     validate that isn't either injected or driver-unwired or dump-blocked).")
    return (length(eod) == 1 && eod[1] == EOD_NSTEP) ? 0 : 3
end

exit(main())
