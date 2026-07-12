# ==========================================================================
# test_glacier_robustness.jl — glacier (istice) driver robustness (gated).
#
# Runs the full CLM.jl driver on a pure-GLACIER column (Bow surfdata with
# PCT_GLACIER=100; test_inputs/glacier/surfdata_glacier100.nc) from a cold start
# under a cold, SNOWY alpine forcing, and asserts:
#
#   (a) snow layers actually FORM (snl goes 0 -> negative), and
#   (b) no field goes NaN/Inf — including across the layer-formation step and
#       the steps after it.
#
# (a) is not a nicety: the previous version of this gate used a light snowfall
# that never pushed snl below 0, so the entire snow-layer path (layer creation,
# the snow band of the soil-temperature solve, SNICAR on a resolved snow layer,
# compaction / combine / divide) was structurally unreachable and the test could
# not fail on it. It didn't — a glacier column NaN'd on the step AFTER the first
# snow layer formed (dead `aerosol_init_cold!` -> allocation-time NaN in the
# snow-layer aerosol masses -> NaN SNICAR albedo/flx_abs -> NaN sabg_lyr -> NaN
# soil-temperature solve) while this test stayed green. `layers_formed` is now
# asserted so the gate cannot go blind again.
#
# The soil-only multisite robustness test never exercises the istice land-unit
# path; this also guards the earlier glacier cold-start fixes (cv bedrock floor,
# cold-start albedo seed, canopy-fall flux mask) against regression. GATED:
# inputs are machine-local; absent -> skip. Runs in an isolated subprocess (the
# harness mutates module globals).
# ==========================================================================
using Test, CLM

@testset "Glacier (istice) driver robustness (gated)" begin
    script = joinpath(@__DIR__, "..", "scripts", "glacier_smoke.jl")
    if !isfile(script)
        @info "glacier robustness: harness absent, skipping" script
        @test_skip isfile(script)
    else
        mod = Module(:GLAC_SMOKE)
        Core.eval(mod, :(using Test, CLM, NCDatasets, Dates, Printf))
        Base.include(mod, script)   # no auto-run: PROGRAM_FILE guard
        r = Base.invokelatest(getfield(mod, :glacier_smoke_result); nsteps=24)
        if r === missing
            @info "glacier robustness: inputs absent, skipping"
            @test_skip r !== missing
        else
            # The gate is only meaningful if the snow-layer path was exercised.
            @test r.layers_formed          # snl went negative at least once
            @test r.min_snl <= -1
            # ...and with it exercised, every field stayed finite for all 24 steps.
            @test r.nbad == 0
        end
    end
end
