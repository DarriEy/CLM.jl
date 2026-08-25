# ==========================================================================
# test_snow_robustness.jl — deep-snowpack / cold-soil driver robustness (gated).
#
# Runs the CLM.jl driver on the Bow natural-veg column (needleleaf-evergreen +
# grass = boreal-like) from a cold start with cold, heavy-snow forcing for enough
# steps to build a MULTI-LAYER snowpack, and asserts no field goes NaN/Inf and that
# snow layers actually form (snl < 0). Exercises the snow-layer build/combine/divide,
# snow compaction, and frozen-soil paths that the boreal and Arctic eval domains hit
# — which the short glacier/multisite smokes never reach. GATED: inputs are
# machine-local; absent -> skip. Isolated subprocess (the harness mutates globals).
# ==========================================================================
using Test, CLM

include(joinpath(@__DIR__, "testdata.jl"))

@testset "Snow / cold-soil driver robustness (gated)" begin
    script = joinpath(@__DIR__, "..", "scripts", "snow_smoke.jl")
    if !isfile(script)
        testdata_missing("snow robustness harness", script)
    else
        mod = Module(:SNOW_SMOKE)
        Core.eval(mod, :(using Test, CLM, NCDatasets, Dates, Printf))
        Base.include(mod, script)
        smoke_fn = Base.invokelatest(getfield, mod, :snow_smoke_ok)
        ok = Base.invokelatest(smoke_fn; nsteps=240)
        if ok === missing
            testdata_unavailable("snow robustness")
        else
            @test ok === true   # 240 steps finite AND a multi-layer snowpack formed
        end
    end
end
