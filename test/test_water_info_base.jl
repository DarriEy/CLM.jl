@testset "WaterInfoBase" begin

    # Define a concrete subtype for testing
    struct TestWaterInfo <: CLM.WaterInfoBaseType
        ratio::Float64
    end

    # Implement required interface
    CLM.get_name(::TestWaterInfo) = "test_tracer"
    CLM.fname(info::TestWaterInfo, basename::String) = basename * "_test"
    CLM.lname(info::TestWaterInfo, basename::String) = basename * " test tracer"
    CLM.is_communicated_with_coupler(::TestWaterInfo) = true
    CLM.is_included_in_consistency_check(::TestWaterInfo) = false

    @testset "concrete subtype instantiation" begin
        info = TestWaterInfo(1.5)
        @test info isa CLM.WaterInfoBaseType
        @test info.ratio == 1.5
    end

    @testset "get_name" begin
        info = TestWaterInfo(1.0)
        @test CLM.get_name(info) == "test_tracer"
    end

    @testset "fname" begin
        info = TestWaterInfo(1.0)
        @test CLM.fname(info, "H2OSOI") == "H2OSOI_test"
    end

    @testset "lname" begin
        info = TestWaterInfo(1.0)
        @test CLM.lname(info, "soil water") == "soil water test tracer"
    end

    @testset "is_communicated_with_coupler" begin
        info = TestWaterInfo(1.0)
        @test CLM.is_communicated_with_coupler(info) == true
    end

    @testset "is_included_in_consistency_check" begin
        info = TestWaterInfo(1.0)
        @test CLM.is_included_in_consistency_check(info) == false
    end

    @testset "get_ratio" begin
        info = TestWaterInfo(2.5)
        @test CLM.get_ratio(info) == 2.5
    end

    @testset "set_metadata!" begin
        # Need a mutable subtype for set_metadata!
        mutable struct MutableTestWaterInfo <: CLM.WaterInfoBaseType
            ratio::Float64
        end

        info = MutableTestWaterInfo(0.0)
        CLM.set_metadata!(info, 3.14)
        @test info.ratio == 3.14
        @test CLM.get_ratio(info) == 3.14
    end

end
