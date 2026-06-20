# ==========================================================================
# test_glacier_robustness.jl — glacier (istice) driver robustness (gated).
#
# Runs the full CLM.jl driver on a pure-GLACIER column (Bow surfdata with
# PCT_GLACIER=100; test_inputs/glacier/surfdata_glacier100.nc) from a cold start
# with a cold, snowy alpine forcing snapshot, and asserts no field goes NaN/Inf.
#
# The soil-only multisite robustness test never exercises the istice land-unit
# path; this guards the glacier cold-start fixes (cv bedrock floor, cold-start
# albedo seed, canopy-fall flux mask) against regression. GATED: inputs are
# machine-local; absent -> skip. Runs in an isolated subprocess (the harness
# mutates module globals).
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
        ok = Base.invokelatest(getfield(mod, :glacier_smoke_ok); nsteps=6)
        if ok === missing
            @info "glacier robustness: inputs absent, skipping"
            @test_skip ok === true
        else
            @test ok === true   # driver ran 6 steps on a pure-glacier column, all finite
        end
    end
end
