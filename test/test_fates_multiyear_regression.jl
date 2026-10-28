# ==========================================================================
# test_fates_multiyear_regression.jl — multi-year FATES demography vs Fortran.
#
# Scores a multi-year Julia run of the fates_parity_1pt case (Bow at Banff,
# 2002-07-01, dtime 3600 s; scripts/fates_parity_multiyear.jl) against the
# committed Fortran CTSM-FATES daily-history reference
# (scripts/validation/fates_fortran_parity/fates_multiyear_fortran_daily.csv.gz,
# 4 model years). See docs/FATES_MULTIYEAR_ENERGY.md "Defect 2".
#
# A 4-year run takes ~40 min, so it is OPT-IN:
#   FATES_MULTIYEAR_CSV=<harness csv>   score an existing harness output, or
#   FATES_MULTIYEAR_RUN=1               run the harness here (FATES_MY_DAYS, default 730).
# Otherwise the scoring block is reported as SKIPPED (`@test_skip`, never a vacuous
# pass). The reference self-checks below always run.
#
# What it catches (both measured on the pre-fix tree, PR body):
#   * recruitment collapsing ~85x below Fortran in year 2 (the never-advanced FATES
#     clock: GDD reset daily, no canopy trimming);
#   * unphysical WOODY growth: Fortran keeps every woody PFT in the 0-5 cm class
#     for all 4 years, so a woody max-dbh increment above 5 cm/yr is a defect.
#   Grass "dbh" is FATES's allometric proxy (structc = 0, leaf targets capped at
#   allom_dbh_maxheight = 20-30 cm): Fortran's C4 grass reaches the >= 100 cm class
#   within 2 years on this case, so large grass dbh is FAITHFUL; only the timing of
#   that class is checked.
# ==========================================================================
using Test
include(joinpath(@__DIR__, "testdata.jl"))
Base.include(@__MODULE__, joinpath(@__DIR__, "..", "scripts", "fates_multiyear_compare.jl"))

const _FMY_REF = joinpath(@__DIR__, "..", "scripts", "validation", "fates_fortran_parity",
                          "fates_multiyear_fortran_daily.csv.gz")

@testset "FATES multi-year vs Fortran reference" begin
    @testset "reference self-check" begin
        @test isfile(_FMY_REF)
        F, nf = fmy_read_csv(_FMY_REF)
        @test nf >= 4 * 365
        rec = F("recruit_d")
        @test all(isfinite, rec) && all(>=(0), rec)
        # the Fortran signature this test leans on: woody stays < 5 cm, grass >= 100 cm by yr 3
        @test maximum(maximum(F("maxsz_$p")) for p in FMY_WOODY) == 0
        @test maximum(maximum(F("maxsz_$p")[1:3*365]) for p in FMY_GRASS) == 100
    end

    csv = get(ENV, "FATES_MULTIYEAR_CSV", "")
    if isempty(csv) && get(ENV, "FATES_MULTIYEAR_RUN", "") == "1"
        fsur = joinpath(symfluence_data_root(), "domain_Bow_at_Banff_lumped_era5", "settings",
                        "CLM", "parameters", "surfdata_clm.nc")
        if !isfile(fsur)
            testdata_missing("FATES multi-year regression", fsur) && return
        end
        csv = joinpath(mktempdir(), "fates_multiyear.csv")
        ENV["FATES_MY_OUT"] = csv
        haskey(ENV, "FATES_MY_DAYS") || (ENV["FATES_MY_DAYS"] = "730")
        m = Module(:FatesMultiyearRun)
        Base.include(m, joinpath(@__DIR__, "..", "scripts", "fates_parity_multiyear.jl"))
        Base.invokelatest(m.main_multiyear)
    end

    if isempty(csv)
        @info "FATES multi-year scoring skipped: set FATES_MULTIYEAR_CSV or FATES_MULTIYEAR_RUN=1"
        @test_skip "multi-year Julia run vs Fortran"
    else
        s = multiyear_summary(csv, _FMY_REF)
        print_summary(s)
        @test length(s) >= 1
        for r in s
            # recruitment within a factor 3 of Fortran every model year
            # (pre-fix year 2: 0.012; post-fix 1.0 / 1.4 / 1.6 / 2.2 for years 1-4)
            @test 1/3 <= r.rec_j / r.rec_f <= 3
            # woody growth stays physical: Fortran never leaves the 0-5 cm class
            @test r.woody_ddbh_j <= 5.0
            @test r.woody_maxbin_j == r.woody_maxbin_f
            # grass reaches the >= 100 cm class in the same model year (+-1) as Fortran
            if r.grass_maxbin_f == 100
                @test r.grass_maxbin_j >= 60
            end
        end
    end
end
