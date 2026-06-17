# ==========================================================================
# test_multisite_robustness.jl — off-Bow driver robustness (gated).
#
# Runs the full CLM.jl driver at a CONTRASTING site (tropical Aripuanã Amazon:
# broadleaf evergreen/deciduous tropical tree + C3/C4 grass, no snow, ~99 kPa)
# from a cold start with a synthetic wet-season forcing snapshot, for both the
# biogeophysics-only (use_cn=false) and full-BGC (use_cn=true) paths, and
# asserts no field goes NaN/Inf. This catches Bow-specific assumptions (alpine,
# snow-dominated, 3 PFTs). It is a ROBUSTNESS check, not Fortran parity — no
# reference dumps exist off-Bow.
#
# GATED: the Aripuanã domain inputs are machine-local; when absent the test is
# skipped (not failed) so CI and other machines stay green. Runs in an isolated
# subprocess (it mutates module globals like rooting_profile_config).
# ==========================================================================
using Test, CLM

@testset "Multi-site driver robustness (gated)" begin
    script = joinpath(@__DIR__, "..", "scripts", "multisite_smoke.jl")
    if !isfile(script)
        @info "multisite robustness: harness script absent, skipping"
        @test_skip isfile(script)
    else
        Base.include(@__MODULE__, script)   # defines multisite_smoke_ok (no auto-run: PROGRAM_FILE guard)
        ok = Base.invokelatest(multisite_smoke_ok; nsteps=6)
        if ok === missing
            @info "multisite robustness: Aripuanã site data absent, skipping"
            @test_skip ok === true
        else
            @test ok === true   # both use_cn=false and use_cn=true ran 6 steps with no NaN/Inf
        end
    end
end
