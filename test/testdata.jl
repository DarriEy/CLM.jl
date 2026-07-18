# ==========================================================================
# testdata.jl — ONE place that knows where the external test data lives.
#
# WHY THIS EXISTS
# ---------------
# ~30 test files used to hardcode the Calgary-laptop data root
# `/Users/darri.eythorsson/compHydro/SYMFLUENCE_data` and then do:
#
#     if !isfile(fsurdat); @warn "skipping"; @test true; return; end
#
# A `@test true` is VACUOUSLY GREEN. On CI (and on any box without that
# laptop) every one of those testsets reported "passed" while testing
# nothing — that is how #238 merged green while BOTH AD suites were
# erroring (see #252). When the laptop is returned the root disappears
# and the blindness becomes permanent and total.
#
# The data has been migrated to a Drive copy, so the roots are now
# env-overridable. The defaults are unchanged, so nothing changes on the
# original box; set the env var to point at any other copy:
#
#     SYMFLUENCE_DATA=/path/to/SYMFLUENCE_data julia --project=. test/runtests.jl
#     CESM_INPUTDATA=/path/to/cesm-inputdata   julia --project=. test/runtests.jl
#
# STRICT MODE
# -----------
# Set `CLM_REQUIRE_TESTDATA=1` to turn every data-gated skip into a hard
# FAILURE. That is the only way to assert "these tests really ran" — use
# it on a box that is supposed to have the data. Default behaviour
# (skip + warn) is unchanged.
#
# Usage in a test file:
#
#     include(joinpath(@__DIR__, "testdata.jl"))
#     fsurdat, paramfile = bow_params()
#     if !isfile(fsurdat) || !isfile(paramfile)
#         testdata_missing("my testset", fsurdat, paramfile) && return
#     end
# ==========================================================================

# Idempotent: several test files include this, and runtests.jl includes them
# all into one scope. Re-including must be a no-op, not a redefinition error.
if !isdefined(@__MODULE__, :SYMFLUENCE_DATA_ROOT)

using Test

"""
    symfluence_data_root() -> String

Root of the SYMFLUENCE data tree (domains, parity dumps, `installs/`).
Override with `ENV["SYMFLUENCE_DATA"]`.
"""
symfluence_data_root() =
    get(ENV, "SYMFLUENCE_DATA", "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data")

"""
    cesm_inputdata_root() -> String

Root of the CESM input-data tree (SNICAR optics, etc).
Override with `ENV["CESM_INPUTDATA"]`.
"""
cesm_inputdata_root() =
    get(ENV, "CESM_INPUTDATA", "/Users/darri.eythorsson/projects/cesm-inputdata")

const SYMFLUENCE_DATA_ROOT = symfluence_data_root()

"""
    symfluence_path(parts...) -> String

Join `parts` onto the SYMFLUENCE data root.
"""
symfluence_path(parts...) = joinpath(symfluence_data_root(), parts...)

"""
    domain_params_dir(domain) -> String

`<root>/<domain>/settings/CLM/parameters` — where `surfdata_clm.nc` and
`clm5_params.nc` live for a domain.
"""
domain_params_dir(domain::AbstractString) =
    symfluence_path(domain, "settings", "CLM", "parameters")

"""
    domain_params(domain) -> (fsurdat, paramfile)

The surface-dataset and parameter-file pair for `domain`.
"""
function domain_params(domain::AbstractString)
    d = domain_params_dir(domain)
    return (joinpath(d, "surfdata_clm.nc"), joinpath(d, "clm5_params.nc"))
end

"""
    bow_domain_dir() -> String

The Bow-at-Banff lumped domain directory (the default test domain).
"""
bow_domain_dir() = symfluence_path("domain_Bow_at_Banff_lumped")

"""
    bow_params() -> (fsurdat, paramfile)

Shorthand for `domain_params("domain_Bow_at_Banff_lumped")`.
"""
bow_params() = domain_params("domain_Bow_at_Banff_lumped")

"""
    bow_forcing(file) -> String

A forcing file under the Bow domain's `data/forcing/CLM_input`.
"""
bow_forcing(file::AbstractString) =
    joinpath(bow_domain_dir(), "data", "forcing", "CLM_input", file)

"""
    snicar_optics() -> String
    snicar_aging()  -> String

SNICAR lookup tables from the CESM input-data tree.
"""
snicar_optics() = joinpath(cesm_inputdata_root(), "lnd", "clm2", "snicardata",
                           "snicar_optics_5bnd_c013122.nc")
snicar_aging()  = joinpath(cesm_inputdata_root(), "lnd", "clm2", "snicardata",
                           "snicar_drdt_bst_fit_60_c070416.nc")

"""
    require_testdata() -> Bool

True when `CLM_REQUIRE_TESTDATA` is set truthy, meaning a missing data
asset must FAIL rather than skip.
"""
require_testdata() = get(ENV, "CLM_REQUIRE_TESTDATA", "0") in ("1", "true", "TRUE", "yes")

"""
    testdata_missing(what, paths...) -> Bool

Report that `what` cannot run because `paths` are absent, and return
`true` so the caller can `return`.

Under `CLM_REQUIRE_TESTDATA=1` this records a real `@test` FAILURE naming
the missing paths. Otherwise it warns and records a *skipped* marker
(`@test_skip`) — deliberately NOT `@test true`, so a skip is visible as a
skip in the summary instead of masquerading as a pass.
"""
function testdata_missing(what::AbstractString, paths::AbstractString...)
    missing_paths = [p for p in paths if !ispath(p)]
    isempty(missing_paths) && (missing_paths = collect(paths))
    if require_testdata()
        @error "CLM_REQUIRE_TESTDATA=1 but test data for $what is missing" missing_paths
        @test "test data present for $what" == "MISSING: $(join(missing_paths, ", "))"
    else
        @warn "Skipping $what: input files not found (set SYMFLUENCE_DATA / CESM_INPUTDATA, or CLM_REQUIRE_TESTDATA=1 to make this fatal)" missing_paths
        @test_skip "test data present for $what"
    end
    return true
end

end # isdefined guard
