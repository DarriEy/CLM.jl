# test_fates_landuse.jl
# FATES land use, host side (fates_driver_landuse.jl) and its configuration:
#
#   A. fates_harvest_mode -> use_logging / use_lu_harvest / num_lu_harvest_cats
#      (CLMFatesGlobals1) and the clm_initialize! land-use option checks.
#   B. The LUH2 reader (dynFatesLandUseInit / dynFatesLandUseInterp) on a synthetic
#      three-year file: year selection (END_OF_TIMESTEP, clamped after the series),
#      zeroing before the series except at init, harvest variable sets per mode,
#      potential vegetation (no file).
#   C. GetLandusePFTData and the bc_in pack on the committed single-point subsets
#      (data/fates/luh_1pt, nearest 4x5 cell of the CESM inputdata LUH2 files).
#   D. (data-gated) logging and LUH2 runs on the fates_parity_1pt site through the
#      first logging event, scored against the committed Fortran CTSM extracts
#      (scripts/validation/fates_fortran_parity/fates_lu_<mode>_fortran.csv).

using Test
using CLM
using NCDatasets

include("testdata.jl")

const _LU1PT = joinpath(@__DIR__, "..", "data", "fates", "luh_1pt")

function _synthetic_luh(path; years = [2000, 2001, 2002])
    NCDataset(path, "c") do o
        nt = length(years)
        defDim(o, "time", nt); defDim(o, "lsmlon", 1); defDim(o, "lsmlat", 1)
        defVar(o, "YEAR", Float64.(years), ("time",))
        names = vcat(CLM.landuse_state_varnames, CLM.landuse_transition_varnames,
                     CLM.landuse_harvest_area_varnames, CLM.landuse_harvest_mass_varnames)
        for (k, v) in enumerate(names)
            # value encodes (variable k, year index t)
            x = Float32[k * 1.0f-3 + t for t in 1:nt]
            defVar(o, v, reshape(x, 1, 1, nt), ("lsmlon", "lsmlat", "time"))
        end
    end
    return path
end

@testset "FATES land use (host LUH2 / harvest)" begin

    @testset "A. configuration" begin
        @test CLM.fates_harvest_ctrl("no_harvest") == (CLM.ifalse, CLM.ifalse, 0)
        @test CLM.fates_harvest_ctrl("event_code") == (CLM.itrue, CLM.ifalse, 0)
        @test CLM.fates_harvest_ctrl("luhdata_area") == (CLM.itrue, CLM.itrue, 5)
        @test CLM.fates_harvest_ctrl("luhdata_mass") == (CLM.itrue, CLM.itrue, 5)
        @test_throws ErrorException CLM.fates_harvest_ctrl("landuse_timeseries")
        @test_throws ErrorException CLM.fates_harvest_ctrl("bogus")

        bad = (; fsurdat = "nonexistent.nc", paramfile = "nonexistent.nc")
        # land-use options need FATES
        @test_throws ErrorException CLM.clm_initialize!(; bad..., use_fates_luh = true)
        @test_throws ErrorException CLM.clm_initialize!(; bad..., fates_harvest_mode = "event_code")
        # lupft needs the landuse x pft file; potential veg forbids an LUH2 file
        @test_throws ErrorException CLM.clm_initialize!(; bad..., use_fates = true,
                                                        use_fates_lupft = true)
        @test_throws ErrorException CLM.clm_initialize!(; bad..., use_fates = true,
            use_fates_potentialveg = true, fluh_timeseries = "x.nc")
    end

    @testset "B. LUH2 reader" begin
        f = _synthetic_luh(joinpath(mktempdir(), "luh.nc"))
        nst, ntr = CLM.num_landuse_state_vars, CLM.num_landuse_transition_vars
        val(k, t) = Float64(Float32(k * 1.0f-3 + t))

        luh = CLM.dyn_fates_landuse_init(1, f; current_year = 2001)
        @test size(luh.landuse_states) == (12, 1)
        @test size(luh.landuse_transitions) == (108, 1)
        @test luh.landuse_states[1, 1] == val(1, 2)               # 2001 -> record 2
        @test luh.landuse_transitions[3, 1] == val(nst + 3, 2)
        @test all(iszero, luh.landuse_harvest)                     # no LUH harvest mode
        CLM.dyn_fates_landuse_interp!(luh, 2002)
        @test luh.landuse_states[12, 1] == val(12, 3)
        CLM.dyn_fates_landuse_interp!(luh, 2050)                   # after: last record
        @test luh.landuse_states[2, 1] == val(2, 3)
        CLM.dyn_fates_landuse_interp!(luh, 1990)                   # before: zeroed
        @test all(iszero, luh.landuse_states) && all(iszero, luh.landuse_transitions)
        # ... but the init read before the series takes the first record
        luh0 = CLM.dyn_fates_landuse_init(1, f; current_year = 1990)
        @test luh0.landuse_states[1, 1] == val(1, 1)

        la = CLM.dyn_fates_landuse_init(1, f; current_year = 2002, harvest_mode = "luhdata_area")
        @test la.landuse_harvest_varnames == CLM.landuse_harvest_area_varnames
        @test la.landuse_harvest_units == CLM.landuse_harvest_area_units
        @test la.landuse_harvest[1, 1] == val(nst + ntr + 1, 3)
        lm = CLM.dyn_fates_landuse_init(1, f; current_year = 2002, harvest_mode = "luhdata_mass")
        @test lm.landuse_harvest_units == CLM.landuse_harvest_mass_units
        @test lm.landuse_harvest[5, 1] == val(nst + ntr + 10, 3)

        pv = CLM.dyn_fates_landuse_init(1, ""; use_potentialveg = true, current_year = 2002)
        @test pv.file === nothing && all(iszero, pv.landuse_states)
        CLM.dyn_fates_landuse_interp!(pv, 2002)
        @test all(iszero, pv.landuse_states)
        @test_throws ErrorException CLM.dyn_fates_landuse_init(1, ""; current_year = 2002)
    end

    @testset "C. single-point subsets + bc_in pack" begin
        luhf = joinpath(_LU1PT, "LUH2_1pt_Bow_hist_simyr1850-2015.nc")
        pftf = joinpath(_LU1PT, "fates_landuse_pft_map_1pt_Bow.nc")
        luh = CLM.dyn_fates_landuse_init(1, luhf; current_year = 2002,
                                         harvest_mode = "luhdata_area")
        # LUH2 2002 states of the 50N 245E cell (primf, secdf; fractions of gridcell)
        @test luh.landuse_states[1, 1] ≈ 0.35397348 rtol = 1e-6
        @test luh.landuse_states[3, 1] ≈ 0.2827094 rtol = 1e-6
        @test luh.landuse_harvest[1, 1] ≈ 0.0031246035 rtol = 1e-6   # primf_harv
        CLM.get_landuse_pft_data!(luh, pftf, 1)
        @test size(luh.landuse_pft_map) == (1, 14, 4)
        @test luh.landuse_bareground[1] ≈ 0.03202889043911796
        @test sum(luh.landuse_pft_map[1, :, 1]) ≈ 1.0 atol = 1e-12  # frac_primr sums to 1

        bc = CLM.bc_in_type()
        bc.hlm_harvest_rates = zeros(5); bc.hlm_harvest_catnames = fill("", 5)
        bc.hlm_luh_states = zeros(12); bc.hlm_luh_state_names = fill("", 12)
        bc.hlm_luh_transitions = zeros(108); bc.hlm_luh_transition_names = fill("", 108)
        old = CLM.hlm_use_luh[]
        try
            CLM.hlm_use_luh[] = CLM.itrue
            CLM.fates_pack_bcin_landuse!(bc, luh, 1)
        finally
            CLM.hlm_use_luh[] = old
        end
        @test bc.hlm_luh_states == luh.landuse_states[:, 1]
        @test bc.hlm_luh_state_names == CLM.landuse_state_varnames
        @test bc.hlm_luh_transition_names[end] == "c3nfx_to_urban"
        @test bc.hlm_harvest_catnames == CLM.landuse_harvest_area_varnames
        @test bc.hlm_harvest_units == CLM.landuse_harvest_area_units
        # the FATES state vector built from the pack sums to one, urban excluded
        bc2 = deepcopy(bc)
        oldn = CLM.hlm_num_luh2_states[]; oldpv = CLM.hlm_use_potentialveg[]
        try
            CLM.hlm_num_luh2_states[] = 12; CLM.hlm_use_potentialveg[] = CLM.ifalse
            sv = CLM.GetLUHStatedata(bc2)
            @test sum(sv) ≈ 1.0
            @test sv[CLM.primaryland] ≈ (0.35397348 + 0.039089218) / (1 - 0.00135295) rtol = 1e-6
        finally
            CLM.hlm_num_luh2_states[] = oldn; CLM.hlm_use_potentialveg[] = oldpv
        end
    end

    @testset "D. live runs vs Fortran (data-gated)" begin
        sitedir = symfluence_path("domain_Bow_at_Banff_lumped_era5", "settings", "CLM",
                                  "parameters")
        fsurdat = joinpath(sitedir, "surfdata_clm.nc")
        fparam  = joinpath(sitedir, "clm5_params.nc")
        forcing = symfluence_path("domain_Bow_at_Banff_lumped_era5", "data", "forcing",
                                  "CLM_input", "clmforc.2002.nc")
        if !(isfile(fsurdat) && isfile(fparam) && isfile(forcing))
            testdata_missing("FATES land-use live runs", fsurdat, fparam, forcing)
            return
        end
        harness = joinpath(@__DIR__, "..", "scripts", "validation", "fates_fortran_parity",
                           "fates_landuse_parity.jl")
        ENV["SYMFLUENCE_DATA"] = symfluence_data_root()
        mod = Module(:_FatesLandUseHarness)
        Core.eval(mod, :(include(p::AbstractString) = Base.include($mod, p)))
        Base.include(mod, harness)
        # 4 days + 1 step: through the 2002-07-04 logging event (nstep 72)
        nsteps = 73
        for (mode, checks) in (
                ("logging", ("FATES_SECONDARY_FOREST_FRACTION" => 5e-3,
                             "FATES_DISTURBANCE_RATE_LOGGING" => 5e-3,
                             "FATES_VEGC" => 5e-3)),
                ("luh", ("FATES_PATCHAREA_LU_1" => 1e-3, "FATES_PATCHAREA_LU_2" => 1e-3,
                         "FATES_PATCHAREA_LU_3" => 1e-3, "FATES_PATCHAREA_LU_5" => 1e-3,
                         "FATES_VEGC" => 5e-3)))
            ref = Base.invokelatest(getglobal(mod, :fortran_csv), mode)
            isfile(ref) || continue
            ft = Base.invokelatest(getglobal(mod, :read_fortran), ref)
            jl = Base.invokelatest(getglobal(mod, :run_julia), mode, nsteps)
            sc = Base.invokelatest(getglobal(mod, :score), jl, ft; io = devnull)
            @testset "$mode" begin
                for (f, tol) in checks
                    @test haskey(sc, f)
                    haskey(sc, f) && @test sc[f][2] < tol
                end
                if mode == "logging"
                    @test jl[nsteps + 1]["FATES_NPATCHES_SECONDARY"] ==
                          ft["FATES_NPATCHES_SECONDARY"][nsteps + 1]
                    # the logging-day litter pulse (no direct logging in ctsm5.3.012
                    # event-code mode; with it the port had 2.4x this pulse)
                    @test jl[nsteps + 1]["FATES_LITTER_IN"] ≈
                          ft["FATES_LITTER_IN"][nsteps + 1] rtol = 5e-3
                end
            end
        end
    end
end
