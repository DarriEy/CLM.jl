# ==========================================================================
# test_urban_robustness.jl — urban (isturb) driver robustness (gated).
#
# Runs the full CLM.jl driver on a surfdata that actually contains urban
# land units with real morphology (the CTSM Mexico-City single-point test
# surfdata, surfdata_1x1_mexicocityMEX_*.nc), from a cold start with a warm
# dry daytime forcing snapshot, and asserts every audited urban column/patch
# field stays finite.
#
# The soil-only / glacier / snow robustness tests never exercise the urban
# roof / wall / road land-unit path. This guards the urban-input reader
# (UrbanInputMod port -> urbinp morphology), the urban cold-start init
# (building temperatures, albedo, AC/heat) and the urban building-thermal
# band solve (incl. the standing-water "layer 0" decouple for roof/wall and
# the band-solver partial-pivoting fix the roof columns first exposed).
# GATED: inputs are machine-local; absent -> skip. Runs in an isolated
# subprocess (the harness mutates module globals).
# ==========================================================================
using Test, CLM

@testset "Urban (isturb) driver robustness (gated)" begin
    script = joinpath(@__DIR__, "..", "scripts", "urban_smoke.jl")
    if !isfile(script)
        @info "urban robustness: harness absent, skipping" script
        @test_skip isfile(script)
    else
        mod = Module(:URB_SMOKE)
        Core.eval(mod, :(using Test, CLM, NCDatasets, Dates, Printf))
        Base.include(mod, script)   # no auto-run: PROGRAM_FILE guard
        ok = Base.invokelatest(getfield(mod, :urban_smoke_ok); nsteps=6)
        if ok === missing
            @info "urban robustness: inputs absent, skipping"
            @test_skip ok === true
        else
            @test ok === true   # driver ran 6 steps on urban land units, all finite
        end
    end
end
