# ==========================================================================
# test_accum_wiring.jl — guard the "ported accumulator, never driven" bug CLASS.
#
# Sibling of test_init_cold_wiring.jl. Same disease, different phase.
#
# Fortran splits accumulator handling three ways (accumulMod.F90):
#     InitAccBuffer   — at init: REGISTER each field (period, acctype, init_value)
#     InitAccVars     — after restart read: EXTRACT buffers into type members
#     UpdateAccVars   — EVERY STEP, from clm_driver: accumulate + extract
#
# CLM.jl ported all three families and then:
#   * never called ANY `*_init_acc_buffer!` or `*_init_acc_vars!`;
#   * called the `*_update_acc_vars!` family from the driver, but several were
#     NO-OP STUBS (`water_update_acc_vars!` dispatched to two empty functions and
#     a comment) or PASS-THROUGHS that assigned the INSTANTANEOUS forcing value to
#     a field named "240hr average";
#   * never called `crop_update_acc_vars!` at all.
#
# The damage was silent and everywhere: LUNA acclimated to instantaneous CO2/O2/
# pressure instead of their 10-day means; `soila10`/`t_a5min`/`snow_5day`/
# `t_ref2m_min`/`t_ref2m_max`/`AnnET` stayed NaN for the whole run and poisoned
# every comparison they fed; `prec365` was a hardcoded zero buffer so CNDV could
# never establish a PFT; crop HUI sat at SPVAL, i.e. "instantly mature".
#
# Unit tests for these routines all PASSED throughout — each called its own
# routine directly with hand-made buffers. That is exactly why the ports looked
# tested and were dead. So, as with InitCold, this file asserts the WIRING:
#
#   1. STATIC — every accumulator routine defined in src/ has a real call site on
#      a LIVE path (init for the Init* family, the driver for Update*).
#   2. STATIC — no `*_update_acc_vars!` is a no-op body. A routine that is called
#      but does nothing is a dead port wearing a live call site as a disguise.
#   3. RUNTIME — after a short run the accumulated fields are FINITE on active
#      points and, for the running means, actually AVERAGE (a 10-day mean of a
#      swinging signal must not equal the instantaneous value).
#
# See src/infrastructure/accumul.jl (`accum_runmean`, `accum_window_steps`) for
# the in-array accumulator primitives and why they reproduce accumulMod exactly.
# ==========================================================================
using Test, CLM, Dates

const ACC_SRC_DIR = normpath(joinpath(@__DIR__, "..", "src"))

# Files that constitute the LIVE path for accumulators.
const LIVE_ACC_FILES = [
    joinpath(ACC_SRC_DIR, "infrastructure", "init_cold.jl"),   # init_acc_buffer! / init_acc_vars!
    joinpath(ACC_SRC_DIR, "driver", "clm_initialize.jl"),      # calls both, once
    joinpath(ACC_SRC_DIR, "driver", "clm_driver.jl"),          # the per-step UpdateAccVars block
    joinpath(ACC_SRC_DIR, "types", "water.jl"),                # WaterType's sub-dispatcher
]

# Accumulator routines that are genuinely dead and STAYING dead, with the reason.
# Not an escape hatch: adding an entry requires justifying it.
const ACC_EXPECTED_DEAD = Dict(
    "accum_rest!" =>
    """
    Restart I/O stub. CLM.jl has no accumulator-buffer restart path (the whole
    restart layer is REPLACE-light). Fortran's accumulRest read/writes the raw
    accum(:) buffers so a restarted run resumes its running means mid-window;
    CLM.jl's in-array accumulators would instead re-seed from the first post-restart
    sample. That is a real (documented) restart-fidelity gap, but it is a restart
    feature, not a dead accumulator — wiring an empty function would not close it.
    """,
)

"Names of every accumulator routine DEFINED under src/, split by family."
function defined_accum_routines()
    buffers = Set{String}()
    initvars = Set{String}()
    updates = Set{String}()
    for (root, _, files) in walkdir(ACC_SRC_DIR), f in files
        endswith(f, ".jl") || continue
        for line in eachline(joinpath(root, f))
            for (re, set) in (
                (r"^function\s+([A-Za-z0-9_]+_init_acc_buffer!)\s*\(", buffers),
                (r"^function\s+([A-Za-z0-9_]+_init_acc_vars!)\s*\(", initvars),
                (r"^function\s+([A-Za-z0-9_]+_update_acc_vars!)\s*\(", updates),
            )
                m = match(re, line)
                m === nothing || push!(set, m.captures[1])
            end
        end
    end
    return (buffers, initvars, updates)
end

"Concatenated source of the live accumulator path, comment lines stripped so a
mention inside a doc/comment cannot masquerade as a call site."
function live_accum_source()
    io = IOBuffer()
    for f in LIVE_ACC_FILES
        isfile(f) || continue
        for line in eachline(f)
            startswith(lstrip(line), "#") && continue
            println(io, line)
        end
    end
    return String(take!(io))
end

"""
Return the body lines of `fname` as defined in src/, with comments and blank
lines stripped. Used to catch a routine that is CALLED but does NOTHING.

Every accumulator routine is a TOP-LEVEL `function ... end`, so the body runs
from the signature line to the first column-0 `end`.
"""
function function_body(fname::AbstractString)
    pat = Regex("^function\\s+" * fname * "\\s*\\(")
    for (root, _, files) in walkdir(ACC_SRC_DIR), f in files
        endswith(f, ".jl") || continue
        lines = readlines(joinpath(root, f))
        for (i, line) in enumerate(lines)
            match(pat, line) === nothing && continue
            body = String[]
            for j in (i + 1):length(lines)
                lines[j] == "end" && break          # column-0 end closes the function
                push!(body, lines[j])
            end
            return filter(body) do l
                s = strip(l)
                !isempty(s) && !startswith(s, "#")
            end
        end
    end
    return String[]
end

@testset "Accumulator wiring (the dead-accumulator guard)" begin

    buffers, initvars, updates = defined_accum_routines()
    src = live_accum_source()

    @testset "the scan finds the accumulator families at all" begin
        @test !isempty(buffers)
        @test !isempty(initvars)
        @test !isempty(updates)
    end

    @testset "every accumulator routine is CALLED on a live path" begin
        dead = String[]
        for fn in sort(collect(union(buffers, initvars, updates)))
            haskey(ACC_EXPECTED_DEAD, fn) && continue
            occursin(fn * "(", src) || push!(dead, fn)
        end
        if !isempty(dead)
            @error """
            DEAD ACCUMULATOR PORT(S) — defined in src/ but never called on the live
            path. Fortran calls InitAccBuffer + InitAccVars once at init and
            UpdateAccVars EVERY STEP. A running mean that is never updated is stuck at
            its init value (or NaN) and silently poisons phenology / LUNA / fire / GDD.
            Wire it into src/infrastructure/init_cold.jl (init_acc_buffer! /
            init_acc_vars!) or the accumulator block of clm_drv!, NOT into a unit test.
            """ dead
        end
        @test isempty(dead)
    end

    # A call site is necessary but NOT sufficient: `water_update_acc_vars!` had one
    # and its body was three no-ops. Assert the Update family actually does work.
    @testset "no *_update_acc_vars! is a no-op body" begin
        hollow = String[]
        for fn in sort(collect(updates))
            body = function_body(fn)
            # A body that is nothing but `return nothing` / `nothing` / `end` is dead.
            meat = filter(body) do l
                s = strip(l)
                s != "end" && s != "nothing" && s != "return nothing"
            end
            isempty(meat) && push!(hollow, fn)
        end
        if !isempty(hollow)
            @error """
            HOLLOW ACCUMULATOR(S) — called every step, but the body does nothing. This
            is the dead-port class wearing a live call site as a disguise: the driver
            "wires" the routine, the wiring test passes, and the accumulator still never
            accumulates. Implement the body or delete the routine.
            """ hollow
        end
        @test isempty(hollow)
    end
end

# --------------------------------------------------------------------------
# RUNTIME: the accumulators must actually AVERAGE, not pass values through.
#
# This is the check that would have caught the original bug. Every atm2lnd
# accumulator used to assign the instantaneous forcing straight into the field
# named "240hr average", so the field was always exactly equal to this step's
# forcing. Drive a swinging signal and assert the running mean DIVERGES from it.
# --------------------------------------------------------------------------
@testset "running means actually average (not pass-through)" begin
    np = 4
    a2l = CLM.Atm2LndData{Float64}()
    temp = CLM.TemperatureData{Float64}()

    # Minimal hand-built accumulator state: one gridcell, one column, np patches.
    a2l.fsd24_patch  = fill(NaN, np)
    a2l.fsd240_patch = fill(NaN, np)
    a2l.fsi24_patch  = fill(NaN, np)
    a2l.fsi240_patch = fill(NaN, np)
    a2l.forc_solad_not_downscaled_grc = zeros(1, 2)
    a2l.forc_solai_grc                = zeros(1, 2)

    pg = ones(Int, np)
    pc = ones(Int, np)
    dtime = 1800

    # A day/night-like swing. With a correct 1-day / 10-day running mean the
    # accumulated value must NOT track the instantaneous value.
    nsteps = 96   # 2 days at dtime=1800
    local last_inst = 0.0
    for s in 1:nsteps
        inst = 400.0 * max(0.0, sin(2π * (s - 12) / 48))
        a2l.forc_solad_not_downscaled_grc[1, 1] = inst
        a2l.forc_solai_grc[1, 1]                = inst
        CLM.atm2lnd_update_acc_vars!(a2l, 1:np, pg, pc; nstep = s, dtime = dtime)
        last_inst = inst
    end

    # The 1-day mean of a half-rectified sine is ~ its mean (400/π ≈ 127), and the
    # instantaneous value at the final step is far from that.
    @test all(isfinite, a2l.fsd24_patch)
    @test all(isfinite, a2l.fsd240_patch)
    # Pass-through would make these EXACTLY equal. A real mean cannot be.
    @test !isapprox(a2l.fsd24_patch[1], last_inst; atol = 1e-8)
    @test !isapprox(a2l.fsd240_patch[1], last_inst; atol = 1e-8)
    # And the 24h and 240h means must differ from each other (different windows).
    @test !isapprox(a2l.fsd24_patch[1], a2l.fsd240_patch[1]; atol = 1e-8)
    # Sanity: a running mean of a non-negative signal stays in range.
    @test 0.0 <= a2l.fsd24_patch[1] <= 400.0
    @test 0.0 <= a2l.fsd240_patch[1] <= 400.0
end

@testset "accum_runmean reproduces accumulMod's runmean recurrence" begin
    # Fortran update_accum_field_runmean:
    #   nsteps = min(nsteps+1, period); val = ((nsteps-1)*val + field)/nsteps
    # Drive it directly against a reference implementation.
    period = 48
    ref_val = 0.0
    ref_n = 0
    cur = NaN            # cold-start missing-value flag, as allocated
    for s in 1:200
        field = 3.0 + 2.0 * sin(s / 7)
        # reference (Fortran)
        ref_n = min(ref_n + 1, period)
        ref_val = ((ref_n - 1) * ref_val + field) / ref_n
        # port
        cur = CLM.accum_runmean(cur, field, s, period)
        @test isapprox(cur, ref_val; rtol = 1e-14)
    end
end

@testset "accum_window_steps converts DAYS to timesteps like init_accum_field" begin
    # Fortran: period = (-accum_period) * SHR_CONST_CDAY / get_step_size()
    @test CLM.accum_window_steps(1, 1800) == 48
    @test CLM.accum_window_steps(10, 1800) == 480
    @test CLM.accum_window_steps(10, 3600) == 240     # must SCALE with dtime
    @test CLM.accum_window_steps(5, 1800) == 240
    @test CLM.accum_window_steps(365, 1800) == 17520
    @test CLM.accum_window_steps(1, 86400) == 1       # never zero
end
