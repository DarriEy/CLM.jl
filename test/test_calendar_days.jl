# Calendar-aware average-days-per-year (mirrors clm_time_manager.F90's
# get_average_days_per_year). Feeds the CN/MIMICS decomposition rate constants;
# the per-pool rate-constant SHIFT under GREGORIAN is proven in test_decomp_bgc.jl
# (k_l1) and test_decomp_mimics.jl (k_frag). Here we pin the helper + wiring.
@testset "Calendar days-per-year" begin
    @testset "get_average_days_per_year(calendar)" begin
        @test CLM.get_average_days_per_year("NO_LEAP")   == 365.0
        @test CLM.get_average_days_per_year("GREGORIAN") == 365.2425
        # case-insensitive (Fortran to_upper)
        @test CLM.get_average_days_per_year("no_leap")   == 365.0
        @test CLM.get_average_days_per_year("Gregorian") == 365.2425
        # unrecognized calendar aborts (Fortran shr_sys_abort)
        @test_throws ErrorException CLM.get_average_days_per_year("JULIAN")
    end

    @testset "module-level CLM_CALENDAR default (NO_LEAP → 365.0)" begin
        # CTSM default is NO_LEAP → the decomposition rate constants stay
        # byte-identical to the historical port.
        @test CLM.CLM_CALENDAR[] == "NO_LEAP"
        @test CLM.get_average_days_per_year() == 365.0

        # Flipping the selector to GREGORIAN flows through the no-arg method,
        # then restore so no other test sees a changed global.
        old = CLM.CLM_CALENDAR[]
        try
            CLM.CLM_CALENDAR[] = "GREGORIAN"
            @test CLM.get_average_days_per_year() == 365.2425
        finally
            CLM.CLM_CALENDAR[] = old
        end
        @test CLM.CLM_CALENDAR[] == "NO_LEAP"
        @test CLM.get_average_days_per_year() == 365.0
    end

    @testset "TimeManager overload" begin
        tm = CLM.TimeManager()                       # defaults to NO_LEAP
        @test CLM.get_average_days_per_year(tm) == 365.0
        tm.calendar = "GREGORIAN"
        @test CLM.get_average_days_per_year(tm) == 365.2425
    end
end
