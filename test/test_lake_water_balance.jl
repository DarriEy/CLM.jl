# ==========================================================================
# test_lake_water_balance.jl — lake column water balance closes (gated).
#
# Runs the full CLM.jl driver on a pure-LAKE column (Bow surfdata with
# PCT_LAKE=100; test_inputs/lake/surfdata_lake100.nc) from a cold start on real
# hourly forcing and asserts the lake column's water balance CLOSES:
#   |errh2o_col| < 1e-5 mm  (CLM's BALANCE_ERROR_THRESH), AND
#   endwb != begwb on some step (endwb is genuinely computed, not a stub).
#
# Before the ComputeWaterMassLake port, endwb was stubbed as endwb = begwb, so
# the lake balance was unfalsifiable; and the patch-indexed lake kernels were
# launched with a COLUMN-sized ndrange, skipping the lake patch entirely (its
# qflx_qrgwl was never written and errh2o_col was NaN — the check passed on NaN).
#
# GATED: inputs are machine-local; absent -> skip. Runs in an isolated subprocess
# (the harness mutates module globals).
# ==========================================================================
using Test, CLM

@testset "Lake column water balance closes (gated)" begin
    script = joinpath(@__DIR__, "..", "scripts", "lake_water_balance.jl")
    if !isfile(script)
        @info "lake water balance: harness absent, skipping" script
        @test_skip isfile(script)
    else
        mod = Module(:LAKE_WB)
        Core.eval(mod, :(using Test, CLM, NCDatasets, Dates, Printf))
        Base.include(mod, script)   # no auto-run: PROGRAM_FILE guard
        # 168 steps (1 week) — long enough to hit a snow event that actually moves
        # the lake column's storage, so the closure test is not vacuous.
        # getfield itself must be invokelatest'd: under Julia 1.12's stricter world-age
        # rules, reading a binding defined by the include above otherwise warns.
        report_fn = Base.invokelatest(getfield, mod, :lake_water_balance_report)
        rep = Base.invokelatest(report_fn; nsteps=168)
        if rep === missing
            @info "lake water balance: inputs absent, skipping"
            @test_skip rep === true
        else
            @test rep.all_finite                       # errh2o is a number, not NaN
            @test rep.endwb_ne_begwb                   # endwb genuinely computed
            @test rep.max_abs_dwb > 0.0                # storage really moved
            @test rep.max_abs_errh2o < 1.0e-5          # and the balance CLOSES
        end
    end
end
