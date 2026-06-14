# ==========================================================================
# fortran_parity_smoke.jl — validate the parity infrastructure end to end.
#
# 1. Build a Bow single-column Julia CLM instance.
# 2. Inject the Fortran `before_step` dump as the shared initial condition.
# 3. Compare the live Julia state against that SAME dump:
#      - fields the reader populates  -> should match to ~0 (round-trip OK),
#      - fields the reader does NOT set -> show the Julia-vs-Fortran gap.
#
#   julia +1.12 --project=. scripts/fortran_parity_smoke.jl
# ==========================================================================

include(joinpath(@__DIR__, "fortran_parity_common.jl"))

const NSTEP = length(ARGS) >= 1 ? ARGS[1] : "8761"
const DUMP_BEFORE = joinpath(DUMPDIR, "pdump_before_step_n$(NSTEP).nc")

println("="^64)
println("  Julia ↔ Fortran parity — INFRASTRUCTURE SMOKE TEST")
println("  Dump: ", DUMP_BEFORE)
println("="^64)

isfile(DUMP_BEFORE) || error("dump not found: $DUMP_BEFORE")
isfile(FSURDAT) || error("surfdata not found: $FSURDAT (check Bow paths)")

println("\n[1] Building Bow single-column instance ...")
(inst, bounds, filt, tm) = build_bow_inst()
println("    cols=", bounds.endc, " patches=", bounds.endp,
        " nlevsno=", CLM.varpar.nlevsno, " nlevgrnd=", CLM.varpar.nlevgrnd)

println("\n[2] Injecting Fortran before_step dump ...")
n_set = inject_dump!(inst, bounds, DUMP_BEFORE)
println("    read_fortran_restart! set $n_set variable groups")

println("\n[3] Round-trip comparison (inst vs the SAME dump):")
results, gmax = compare_inst_to_dump(inst, DUMP_BEFORE; label="round-trip", tol=1e-9)

println()
println("="^64)
# Reader-populated fields should be ~0; report how many round-tripped clean.
clean = count(r -> r[3] <= 1e-9, results)
println("  $clean / $(length(results)) fields round-tripped to <= 1e-9")
println("  (fields with larger diffs are ones read_fortran_restart! does not set,")
println("   i.e. they hold Julia's cold-start value, not the Fortran state.)")
println("="^64)
