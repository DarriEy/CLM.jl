# test_fates_reduced_modes.jl
# FATES reduced-complexity modes: satellite phenology (use_fates_sp), no competition
# (use_fates_nocomp) and fixed biogeography (use_fates_fixed_biogeog).
#
#   A. Configuration: the CTSM/FATES dependencies (SP => nocomp .and. fixed_biogeog;
#      use_fates_bgc = .not. use_fates_sp) are enforced, and the HLM patch reservation
#      follows SetFatesGlobalElements1 (fates_maxPatchesPerSite) per mode.
#   B. (data-gated) Each mode cold-starts and runs through the first daily FATES step on
#      the fates_parity_1pt site (Bow at Banff, 2002-07-01) and is scored against the
#      committed hourly Fortran CTSM extracts of the same case
#      (scripts/validation/fates_fortran_parity/fates_<mode>_fortran.csv, produced by
#      stage_fates_modes_run.sh; see the README "REDUCED-COMPLEXITY MODES").

using Test
using CLM

include("testdata.jl")

@testset "FATES reduced-complexity modes" begin

    @testset "A. configuration" begin
        # fates_maxPatchesPerSite per mode (default CDL: 14 PFTs, maxpatches 9,4,1,1,1)
        mp(; kw...) = CLM.fates_max_patches_per_site(; surf_numpft = 15, surf_numcft = 2, kw...)
        @test mp() == CLM.fates_maxpatch_total() + 1 == 17
        @test mp(use_nocomp = true) == 14 + 4 + 1 + 1 + 1 + 1     # primary raised to numpft
        @test mp(use_sp = true, use_nocomp = true) == max(15 + 2, 14 + 1)

        # SP needs nocomp + fixed biogeography; SP forbids soil BGC; the modes need FATES.
        bad = (; fsurdat = "nonexistent.nc", paramfile = "nonexistent.nc")
        @test_throws ErrorException CLM.clm_initialize!(; bad..., use_fates = true,
                                                        use_fates_sp = true)
        @test_throws ErrorException CLM.clm_initialize!(; bad..., use_fates = true,
                                                        use_fates_sp = true,
                                                        use_fates_nocomp = true)
        @test_throws ErrorException CLM.clm_initialize!(; bad..., use_fates = true,
            use_fates_sp = true, use_fates_nocomp = true, use_fates_fixed_biogeog = true,
            use_fates_bgc = true)
        @test_throws ErrorException CLM.clm_initialize!(; bad..., use_fates = false,
                                                        use_fates_nocomp = true)
        @test_throws ErrorException CLM.clm_initialize!(; bad..., use_fates = true,
            use_fates_fixed_biogeog = true, fates_biogeog_screen = :drop_cold_deciduous)

        # Driver: SP turns the FATES soil BGC off (controlMod); non-SP keeps it.
        @test !CLM.CLMDriverConfig(use_fates = true, use_fates_sp = true).use_fates_bgc
        @test CLM.CLMDriverConfig(use_fates = true).use_fates_bgc
    end

    @testset "B. live runs vs Fortran (data-gated)" begin
        sitedir = symfluence_path("domain_Bow_at_Banff_lumped_era5", "settings", "CLM",
                                  "parameters")
        fsurdat = joinpath(sitedir, "surfdata_clm.nc")
        fparam  = joinpath(sitedir, "clm5_params.nc")
        forcing = symfluence_path("domain_Bow_at_Banff_lumped_era5", "data", "forcing",
                                  "CLM_input", "clmforc.2002.nc")
        if !(isfile(fsurdat) && isfile(fparam) && isfile(forcing))
            testdata_missing("FATES reduced-mode live runs", fsurdat, fparam, forcing)
            return
        end
        harness = joinpath(@__DIR__, "..", "scripts", "validation", "fates_fortran_parity",
                           "fates_modes_parity.jl")
        ENV["SYMFLUENCE_DATA"] = symfluence_data_root()
        m = Module(:_FatesModesHarness)
        Core.eval(m, :(include(p::AbstractString) = Base.include($m, p)))
        Base.include(m, harness)
        run_julia   = getglobal(m, :run_julia)
        read_ft     = getglobal(m, :read_fortran)
        fortran_csv = getglobal(m, :fortran_csv)
        NSTEPS = 30                        # crosses the first daily FATES step (nstep 25)
        for mode in ("sp", "nocomp", "fixedbg")
            @testset "$mode" begin
                recs = Base.invokelatest(run_julia, mode, NSTEPS)
                ft = Base.invokelatest(read_ft, Base.invokelatest(fortran_csv, mode))
                n = min(length(recs), length(ft["FATES_GPP"]))
                J(f) = [recs[k][f] for k in 2:n]
                F(f) = ft[f][2:n]
                for f in ("FATES_GPP", "TLAI", "ELAI", "FSH", "EFLX_LH_TOT", "FATES_VEGC")
                    @test all(isfinite, J(f))
                end
                # Census: identical patch/cohort counts every step.
                @test J("FATES_NPATCHES") == F("FATES_NPATCHES")
                @test J("FATES_NCOHORTS") == F("FATES_NCOHORTS")
                # Vegetation state and canopy structure track Fortran closely (SP: to
                # ~1e-8, the FATES side only reads the prescribed LAI); the first-day
                # fluxes carry the known host top-soil thaw residual (see the README
                # "REDUCED-COMPLEXITY MODES" for the 10-day numbers), so they are held
                # to a loose window-mean bound here.
                relmean(f) = abs(sum(J(f)) - sum(F(f))) / abs(sum(F(f)))
                @test relmean("FATES_LAI") < 0.01
                @test relmean("ELAI") < 0.02
                @test relmean("TLAI") < 0.02
                @test relmean("FATES_VEGC") < 0.01
                mode == "sp" && @test relmean("FATES_VEGC") < 1e-6
                @test relmean("FATES_GPP") < 0.3
                @test relmean("FSH") < 0.5
                @test relmean("EFLX_LH_TOT") < 0.15
                # Patch area by PFT (nocomp) and plant density by PFT, end of window.
                for k in 1:14
                    if mode != "fixedbg"
                        f = "FATES_NOCOMP_PATCHAREA_PF_$k"
                        @test recs[n][f] ≈ ft[f][n] atol = 1e-6
                    end
                    f = "FATES_NPLANT_PF_$k"
                    @test recs[n][f] ≈ ft[f][n] rtol = 0.02 atol = 1e-8
                end
            end
        end
    end
end
