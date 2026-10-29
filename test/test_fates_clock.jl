# ==========================================================================
# test_fates_clock.jl — the FATES clock must advance with the host date.
#
# Fortran: clmfates_interfaceMod.F90 `GetAndSetTime` (3730-3786) runs at FATES
# init and at the top of EVERY daily `dynamics_driv` (1068). The port called
# SetFatesTime once at init with a hardcoded day_of_year = 1 and never again, so
# for the whole run FATES believed it was 2 January:
#   * phenology's `hlm_day_of_year == gddstart` (NH gddstart = 1) reset the growing
#     degree days EVERY day -> cold-deciduous PFTs could never leaf back on;
#   * the annual `trim_canopy` (doy == days_per_year-1) never fired;
#   * the cold-start leaf-on/off dates (`cleafon - hlm_day_of_year`) were offset.
# On the fates_parity_1pt case this collapsed year-2 recruitment ~50x below the
# Fortran reference (docs/FATES_MULTIYEAR_ENERGY.md, "Defect 2").
# ==========================================================================
using Test
using CLM

@testset "FATES clock (GetAndSetTime)" begin
    saved = (CLM.hlm_current_year[], CLM.hlm_current_month[], CLM.hlm_current_day[],
             CLM.hlm_current_tod[], CLM.hlm_current_date[], CLM.hlm_reference_date[],
             CLM.hlm_model_day[], CLM.hlm_day_of_year[], CLM.hlm_days_per_year[],
             CLM.hlm_freq_day[])
    try
        CLM.hlm_reference_date[] = 20020701

        # the run start: 0-based NOLEAP day of year, model_day 0
        CLM.fates_get_and_set_time!(2002, 7, 1, 0)
        @test CLM.hlm_day_of_year[] == 181        # floor(datediff(jan01, jul01))
        @test CLM.hlm_model_day[] == 0.0
        @test CLM.hlm_current_date[] == 20020701
        @test CLM.hlm_days_per_year[] == 365
        @test CLM.hlm_freq_day[] ≈ 1 / 365

        # one day later, an hour into the day
        CLM.fates_get_and_set_time!(2002, 7, 2, 3600)
        @test CLM.hlm_day_of_year[] == 182
        @test CLM.hlm_model_day[] ≈ 1 + 1 / 24

        # the trim_canopy gate: second-to-last day of the year
        CLM.fates_get_and_set_time!(2002, 12, 31, 0)
        @test CLM.hlm_day_of_year[] == CLM.hlm_days_per_year[] - 1

        # year = 0 ("host does not track it"): carried forward, advanced on wrap
        CLM.fates_get_and_set_time!(0, 1, 1, 0)
        @test CLM.hlm_current_year[] == 2003
        @test CLM.hlm_day_of_year[] == 0          # NH gddstart (1) is NOT day 0
        @test CLM.hlm_model_day[] == 184.0
        CLM.fates_get_and_set_time!(0, 3, 1, 0)
        @test CLM.hlm_current_year[] == 2003
        @test CLM.hlm_day_of_year[] == 59         # NOLEAP: Mar 1 = day 59 (0-based)

        # an explicit year wins
        CLM.fates_get_and_set_time!(2005, 7, 1, 0)
        @test CLM.hlm_current_year[] == 2005
        @test CLM.hlm_model_day[] == 3 * 365.0

        # NH growing-degree-day reset fires once a year, not every day
        CLM.hlm_reference_date[] = 20020101
        nreset = 0
        for d in 0:729
            # walk a NOLEAP calendar day by day
            yr = 2002 + d ÷ 365; doy0 = d % 365
            mon = findlast(<=(doy0), collect(CLM._FATES_NOLEAP_CUMDAYS))
            dd = doy0 - CLM._FATES_NOLEAP_CUMDAYS[mon] + 1
            CLM.fates_get_and_set_time!(yr, mon, dd, 0)
            CLM.hlm_day_of_year[] == 1 && (nreset += 1)   # NH gddstart
        end
        @test nreset == 2

        # GREGORIAN host calendar (the fates_parity_1pt case): leap years count
        CLM.FATES_HLM_CALENDAR[] = :gregorian
        CLM.hlm_reference_date[] = 20020701
        CLM.fates_get_and_set_time!(2004, 3, 1, 0)
        @test CLM.hlm_day_of_year[] == 60        # Feb 29 exists in 2004
        @test CLM.hlm_days_per_year[] == 366
        @test CLM.hlm_freq_day[] ≈ 1 / 366
        @test CLM.hlm_model_day[] == 609.0      # 2002-07-01 -> 2004-03-01
        CLM.fates_get_and_set_time!(2004, 12, 31, 0)
        @test CLM.hlm_day_of_year[] == CLM.hlm_days_per_year[] - 1   # trim_canopy gate
        CLM.fates_get_and_set_time!(0, 1, 1, 0)                      # year wrap
        @test CLM.hlm_current_year[] == 2005 && CLM.hlm_days_per_year[] == 365
    finally
        CLM.FATES_HLM_CALENDAR[] = :noleap
        (CLM.hlm_current_year[], CLM.hlm_current_month[], CLM.hlm_current_day[],
         CLM.hlm_current_tod[], CLM.hlm_current_date[], CLM.hlm_reference_date[],
         CLM.hlm_model_day[], CLM.hlm_day_of_year[], CLM.hlm_days_per_year[],
         CLM.hlm_freq_day[]) = saved
    end
end
