# ==========================================================================
# test_multisite_robustness.jl — off-Bow driver robustness (gated).
#
# Runs the full CLM.jl driver at CONTRASTING sites from a cold start with a
# synthetic forcing snapshot, for both the biogeophysics-only (use_cn=false) and
# full-BGC (use_cn=true) paths, and asserts no field goes NaN/Inf. This catches
# Bow-specific assumptions (alpine, snow-dominated, 3 PFTs). Robustness check,
# NOT Fortran parity — no reference dumps exist off-Bow.
#
#   - Aripuanã Amazon  (tropical: broadleaf evergreen/deciduous tree, no snow)
#   - Stillwater, OK   (semi-arid: temperate shrub + C3/C4 grass, hot summer)
#
# Each site harness defines its own `main`/forcing helpers, so each is included
# into a FRESH module to avoid name collisions. GATED: site inputs are
# machine-local; absent → skip (CI/other machines stay green). Runs in an
# isolated subprocess (the harnesses mutate module globals).
# ==========================================================================
using Test, CLM

@testset "Multi-site driver robustness (gated)" begin
    sites = [
        ("Aripuanã (tropical)",  "multisite_smoke.jl",            :multisite_smoke_ok),
        ("Stillwater (semiarid)", "multisite_smoke_stillwater.jl", :multisite_stillwater_ok),
    ]
    for (label, fname, okfn) in sites
        script = joinpath(@__DIR__, "..", "scripts", fname)
        @testset "$label" begin
            if !isfile(script)
                @info "multisite robustness: harness absent, skipping" label script
                @test_skip isfile(script)
                continue
            end
            # fresh module: each harness has its own `main`/`set_*_forcing!`
            mod = Module(Symbol("MS_", okfn))
            Core.eval(mod, :(using Test, CLM, NCDatasets, Dates, Printf))
            Base.include(mod, script)   # no auto-run: PROGRAM_FILE guard
            ok = Base.invokelatest(getfield(mod, okfn); nsteps=6)
            if ok === missing
                @info "multisite robustness: site data absent, skipping" label
                @test_skip ok === true
            else
                @test ok === true   # both use_cn=false and use_cn=true ran with no NaN/Inf
            end
        end
    end
end
