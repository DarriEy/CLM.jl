@testset "Accumulation" begin

    # Helper: create a simple manager with one field
    function make_mgr(; accum_type="timeavg", period=4, numlev=1,
                        npts=3, init_value=0.0, active=nothing)
        mgr = CLM.AccumManager()
        if active === nothing
            active = fill(true, npts)
        end
        CLM.init_accum_field!(mgr;
            name = "TEST",
            units = "K",
            desc = "test field",
            accum_type = accum_type,
            accum_period = period,
            numlev = numlev,
            subgrid_type = "column",
            init_value = init_value,
            active = active,
            npts = npts)
        return mgr
    end

    # ------------------------------------------------------------------
    @testset "init_accum_field!" begin
        mgr = make_mgr(accum_type="timeavg", period=4, npts=5)
        @test length(mgr.fields) == 1
        af = mgr.fields[1]
        @test af.name == "TEST"
        @test af.units == "K"
        @test af.acctype == CLM.ACCTYPE_TIMEAVG
        @test af.period == 4
        @test af.numlev == 1
        @test af.beg1d == 1
        @test af.end1d == 5
        @test size(af.val) == (5, 1)
        @test all(af.val .== 0.0)
        @test all(.!af.reset)
        @test all(af.nsteps .== 0)

        # Invalid accum_type
        mgr2 = CLM.AccumManager()
        @test_throws ErrorException CLM.init_accum_field!(mgr2;
            name="BAD", units="", desc="", accum_type="bogus",
            accum_period=1, numlev=1, subgrid_type="column",
            init_value=0.0, active=Bool[true], npts=1)

        # timeavg requires init_value == 0
        mgr3 = CLM.AccumManager()
        @test_throws ErrorException CLM.init_accum_field!(mgr3;
            name="BAD", units="", desc="", accum_type="timeavg",
            accum_period=1, numlev=1, subgrid_type="column",
            init_value=5.0, active=Bool[true], npts=1)

        # runaccum requires init_value == 0
        mgr4 = CLM.AccumManager()
        @test_throws ErrorException CLM.init_accum_field!(mgr4;
            name="BAD", units="", desc="", accum_type="runaccum",
            accum_period=1, numlev=1, subgrid_type="column",
            init_value=1.0, active=Bool[true], npts=1)

        # runmean allows non-zero init_value
        mgr5 = CLM.AccumManager()
        CLM.init_accum_field!(mgr5;
            name="RM", units="K", desc="runmean test",
            accum_type="runmean", accum_period=10, numlev=1,
            subgrid_type="column", init_value=300.0,
            active=Bool[true], npts=1)
        @test mgr5.fields[1].initval == 300.0

        # Negative period (days → timesteps)
        mgr6 = CLM.AccumManager()
        CLM.init_accum_field!(mgr6;
            name="DY", units="K", desc="day period",
            accum_type="runmean", accum_period=-1, numlev=1,
            subgrid_type="column", init_value=0.0,
            active=Bool[true], npts=1, step_size=1800)
        # -1 day → 86400/1800 = 48 timesteps
        @test mgr6.fields[1].period == 48
    end

    # ------------------------------------------------------------------
    @testset "find_field" begin
        mgr = make_mgr()
        @test CLM.find_field(mgr, "TEST") == 1
        @test_throws ErrorException CLM.find_field(mgr, "NONEXISTENT")
    end

    # ------------------------------------------------------------------
    @testset "acctype_to_string" begin
        @test CLM.acctype_to_string(CLM.ACCTYPE_TIMEAVG)  == "timeavg"
        @test CLM.acctype_to_string(CLM.ACCTYPE_RUNMEAN)  == "runmean"
        @test CLM.acctype_to_string(CLM.ACCTYPE_RUNACCUM) == "runaccum"
        @test CLM.acctype_to_string(999) == "unknown"
    end

    # ------------------------------------------------------------------
    @testset "timeavg — basic accumulation" begin
        # Period = 4, 1 point, all active
        mgr = make_mgr(accum_type="timeavg", period=4, npts=1)

        # Update steps 1–4 with values 1,2,3,4
        for (step, v) in [(1, 1.0), (2, 2.0), (3, 3.0), (4, 4.0)]
            CLM.update_accum_field!(mgr, "TEST", [v], step)
        end

        # Extract at step 4 (end of period): should be average
        out = zeros(1)
        CLM.extract_accum_field!(mgr, "TEST", out, 4)
        @test out[1] ≈ 2.5  # (1+2+3+4)/4

        # Extract at step 3 (mid-period): should be SPVAL
        out2 = zeros(1)
        CLM.extract_accum_field!(mgr, "TEST", out2, 3)
        @test out2[1] == CLM.SPVAL
    end

    # ------------------------------------------------------------------
    @testset "timeavg — period=1 (every step)" begin
        mgr = make_mgr(accum_type="timeavg", period=1, npts=1)

        CLM.update_accum_field!(mgr, "TEST", [7.0], 1)
        out = zeros(1)
        CLM.extract_accum_field!(mgr, "TEST", out, 1)
        @test out[1] ≈ 7.0

        CLM.update_accum_field!(mgr, "TEST", [3.0], 2)
        CLM.extract_accum_field!(mgr, "TEST", out, 2)
        @test out[1] ≈ 3.0
    end

    # ------------------------------------------------------------------
    @testset "timeavg — multiple points" begin
        mgr = make_mgr(accum_type="timeavg", period=2, npts=3)

        CLM.update_accum_field!(mgr, "TEST", [10.0, 20.0, 30.0], 1)
        CLM.update_accum_field!(mgr, "TEST", [30.0, 40.0, 50.0], 2)

        out = zeros(3)
        CLM.extract_accum_field!(mgr, "TEST", out, 2)
        @test out[1] ≈ 20.0  # (10+30)/2
        @test out[2] ≈ 30.0  # (20+40)/2
        @test out[3] ≈ 40.0  # (30+50)/2
    end

    # ------------------------------------------------------------------
    @testset "runmean — basic" begin
        # Period = 3, 1 point
        mgr = make_mgr(accum_type="runmean", period=3, npts=1, init_value=0.0)

        # Step 1: val = (0*0 + 6)/1 = 6.0
        CLM.update_accum_field!(mgr, "TEST", [6.0], 1)
        out = zeros(1)
        CLM.extract_accum_field!(mgr, "TEST", out, 1)
        @test out[1] ≈ 6.0

        # Step 2: val = (1*6 + 3)/2 = 4.5
        CLM.update_accum_field!(mgr, "TEST", [3.0], 2)
        CLM.extract_accum_field!(mgr, "TEST", out, 2)
        @test out[1] ≈ 4.5

        # Step 3: val = (2*4.5 + 9)/3 = 6.0
        CLM.update_accum_field!(mgr, "TEST", [9.0], 3)
        CLM.extract_accum_field!(mgr, "TEST", out, 3)
        @test out[1] ≈ 6.0

        # Step 4: nsteps capped at period=3, accper=3
        # val = (2*6.0 + 12)/3 = 8.0
        CLM.update_accum_field!(mgr, "TEST", [12.0], 4)
        CLM.extract_accum_field!(mgr, "TEST", out, 4)
        @test out[1] ≈ 8.0
    end

    # ------------------------------------------------------------------
    @testset "runmean — with non-zero init_value" begin
        mgr = CLM.AccumManager()
        CLM.init_accum_field!(mgr;
            name="TINIT", units="K", desc="test init",
            accum_type="runmean", accum_period=3, numlev=1,
            subgrid_type="column", init_value=300.0,
            active=Bool[true], npts=1)

        # Before any update, extract should return initval
        out = zeros(1)
        CLM.extract_accum_field!(mgr, "TINIT", out, 0)
        @test out[1] ≈ 300.0

        # Step 1: val = (0*300 + 280)/1 = 280
        CLM.update_accum_field!(mgr, "TINIT", [280.0], 1)
        CLM.extract_accum_field!(mgr, "TINIT", out, 1)
        @test out[1] ≈ 280.0
    end

    # ------------------------------------------------------------------
    @testset "runmean — reset" begin
        mgr = make_mgr(accum_type="runmean", period=10, npts=1, init_value=0.0)

        # Accumulate a few steps
        CLM.update_accum_field!(mgr, "TEST", [10.0], 1)
        CLM.update_accum_field!(mgr, "TEST", [20.0], 2)

        # Mark reset
        CLM.markreset_accum_field!(mgr, "TEST")

        # Next update should reset then accumulate
        CLM.update_accum_field!(mgr, "TEST", [5.0], 3)

        out = zeros(1)
        CLM.extract_accum_field!(mgr, "TEST", out, 3)
        # After reset: val=initval=0, nsteps=0
        # Then: nsteps=1, accper=1, val=(0*0+5)/1=5
        @test out[1] ≈ 5.0
    end

    # ------------------------------------------------------------------
    @testset "runaccum — basic" begin
        mgr = make_mgr(accum_type="runaccum", period=1, npts=1)

        # Step 1: val = min(max(0+5, 0), 99999) = 5
        CLM.update_accum_field!(mgr, "TEST", [5.0], 1)
        out = zeros(1)
        CLM.extract_accum_field!(mgr, "TEST", out, 1)
        @test out[1] ≈ 5.0

        # Step 2: val = min(max(5+3, 0), 99999) = 8
        CLM.update_accum_field!(mgr, "TEST", [3.0], 2)
        CLM.extract_accum_field!(mgr, "TEST", out, 2)
        @test out[1] ≈ 8.0

        # Step 3: val = min(max(8-2, 0), 99999) = 6
        CLM.update_accum_field!(mgr, "TEST", [-2.0], 3)
        CLM.extract_accum_field!(mgr, "TEST", out, 3)
        @test out[1] ≈ 6.0
    end

    # ------------------------------------------------------------------
    @testset "runaccum — clamping" begin
        mgr = make_mgr(accum_type="runaccum", period=1, npts=1)

        # Should clamp to 0 on negative
        CLM.update_accum_field!(mgr, "TEST", [-100.0], 1)
        out = zeros(1)
        CLM.extract_accum_field!(mgr, "TEST", out, 1)
        @test out[1] ≈ 0.0

        # Should clamp to 99999
        CLM.update_accum_field!(mgr, "TEST", [200000.0], 2)
        CLM.extract_accum_field!(mgr, "TEST", out, 2)
        @test out[1] ≈ 99999.0
    end

    # ------------------------------------------------------------------
    @testset "runaccum — reset" begin
        mgr = make_mgr(accum_type="runaccum", period=1, npts=1)

        CLM.update_accum_field!(mgr, "TEST", [10.0], 1)
        CLM.update_accum_field!(mgr, "TEST", [5.0], 2)

        # Mark reset
        CLM.markreset_accum_field!(mgr, "TEST")

        # Step 3: reset happens (val=0), no update in same call
        CLM.update_accum_field!(mgr, "TEST", [99.0], 3)
        out = zeros(1)
        CLM.extract_accum_field!(mgr, "TEST", out, 3)
        @test out[1] ≈ 0.0  # reset, no accumulate in same step

        # Step 4: normal accumulate
        CLM.update_accum_field!(mgr, "TEST", [7.0], 4)
        CLM.extract_accum_field!(mgr, "TEST", out, 4)
        @test out[1] ≈ 7.0
    end

    # ------------------------------------------------------------------
    @testset "inactive points" begin
        active = [true, false, true]
        mgr = make_mgr(accum_type="runmean", period=10, npts=3,
                        init_value=0.0, active=active)

        CLM.update_accum_field!(mgr, "TEST", [10.0, 999.0, 30.0], 1)
        out = zeros(3)
        CLM.extract_accum_field!(mgr, "TEST", out, 1)
        @test out[1] ≈ 10.0
        @test out[2] ≈ 0.0   # inactive: stays at init_value
        @test out[3] ≈ 30.0
    end

    # ------------------------------------------------------------------
    @testset "multi-level field" begin
        mgr = CLM.AccumManager()
        CLM.init_accum_field!(mgr;
            name="ML", units="K", desc="multi-level test",
            accum_type="runmean", accum_period=5, numlev=2,
            subgrid_type="column", init_value=0.0,
            active=Bool[true, true], npts=2,
            type2d="levsoi", scale_by_thickness=false)

        af = mgr.fields[1]
        @test af.numlev == 2
        @test size(af.val) == (2, 2)

        # Update with 2×2 matrix [npts × numlev]
        field = [1.0 2.0; 3.0 4.0]
        CLM.update_accum_field!(mgr, "ML", field, 1)

        out = zeros(2, 2)
        CLM.extract_accum_field!(mgr, "ML", out, 1)
        @test out[1, 1] ≈ 1.0
        @test out[1, 2] ≈ 2.0
        @test out[2, 1] ≈ 3.0
        @test out[2, 2] ≈ 4.0
    end

    # ------------------------------------------------------------------
    @testset "markreset_accum_field!" begin
        mgr = make_mgr(accum_type="runmean", period=10, npts=3, init_value=0.0)

        # Reset all
        CLM.markreset_accum_field!(mgr, "TEST")
        r = CLM.get_accum_reset(mgr, "TEST")
        @test all(r)

        # Clear by updating
        CLM.update_accum_field!(mgr, "TEST", [1.0, 1.0, 1.0], 1)

        # Reset single point
        CLM.markreset_accum_field!(mgr, "TEST"; kf=2)
        r2 = CLM.get_accum_reset(mgr, "TEST")
        @test !r2[1, 1]
        @test r2[2, 1]
        @test !r2[3, 1]

        # Multi-level reset by level
        mgr2 = CLM.AccumManager()
        CLM.init_accum_field!(mgr2;
            name="ML2", units="K", desc="test",
            accum_type="runmean", accum_period=5, numlev=3,
            subgrid_type="column", init_value=0.0,
            active=Bool[true, true], npts=2,
            type2d="levsoi", scale_by_thickness=false)

        CLM.markreset_accum_field!(mgr2, "ML2"; level=2)
        r3 = CLM.get_accum_reset(mgr2, "ML2")
        @test !r3[1, 1]
        @test r3[1, 2]
        @test !r3[1, 3]
        @test !r3[2, 1]
        @test r3[2, 2]
        @test !r3[2, 3]
    end

    # ------------------------------------------------------------------
    @testset "get_accum_reset returns a copy" begin
        mgr = make_mgr(accum_type="runmean", period=10, npts=2, init_value=0.0)
        CLM.markreset_accum_field!(mgr, "TEST")
        r = CLM.get_accum_reset(mgr, "TEST")
        r .= false  # modify the copy
        # Original should still be true
        @test all(mgr.fields[1].reset)
    end

    # ------------------------------------------------------------------
    @testset "clean_accum_fields!" begin
        mgr = make_mgr()
        @test length(mgr.fields) == 1
        CLM.clean_accum_fields!(mgr)
        @test length(mgr.fields) == 0
    end

    # ------------------------------------------------------------------
    @testset "print_accum_fields (no error)" begin
        mgr = make_mgr()
        # Just verify it doesn't error
        CLM.print_accum_fields(mgr)
        @test true
    end

    # ------------------------------------------------------------------
    @testset "accum_rest! (stub)" begin
        mgr = make_mgr()
        # Just verify it doesn't error
        @test CLM.accum_rest!(mgr) === nothing
    end

    # ------------------------------------------------------------------
    @testset "timeavg — two consecutive periods" begin
        mgr = make_mgr(accum_type="timeavg", period=3, npts=1)

        # Period 1: steps 1-3, values 3,6,9 → avg = 6
        CLM.update_accum_field!(mgr, "TEST", [3.0], 1)
        CLM.update_accum_field!(mgr, "TEST", [6.0], 2)
        CLM.update_accum_field!(mgr, "TEST", [9.0], 3)

        out = zeros(1)
        CLM.extract_accum_field!(mgr, "TEST", out, 3)
        @test out[1] ≈ 6.0

        # Period 2: steps 4-6, values 12,15,18 → avg = 15
        CLM.update_accum_field!(mgr, "TEST", [12.0], 4)
        CLM.update_accum_field!(mgr, "TEST", [15.0], 5)
        CLM.update_accum_field!(mgr, "TEST", [18.0], 6)

        CLM.extract_accum_field!(mgr, "TEST", out, 6)
        @test out[1] ≈ 15.0
    end

    # ------------------------------------------------------------------
    @testset "multiple fields in one manager" begin
        mgr = CLM.AccumManager()
        active = Bool[true, true]

        CLM.init_accum_field!(mgr; name="F1", units="K", desc="field 1",
            accum_type="runmean", accum_period=5, numlev=1,
            subgrid_type="column", init_value=0.0,
            active=active, npts=2)

        CLM.init_accum_field!(mgr; name="F2", units="mm", desc="field 2",
            accum_type="runaccum", accum_period=1, numlev=1,
            subgrid_type="column", init_value=0.0,
            active=active, npts=2)

        @test length(mgr.fields) == 2
        @test CLM.find_field(mgr, "F1") == 1
        @test CLM.find_field(mgr, "F2") == 2

        # Update independently
        CLM.update_accum_field!(mgr, "F1", [10.0, 20.0], 1)
        CLM.update_accum_field!(mgr, "F2", [5.0, 3.0], 1)

        out1 = zeros(2)
        out2 = zeros(2)
        CLM.extract_accum_field!(mgr, "F1", out1, 1)
        CLM.extract_accum_field!(mgr, "F2", out2, 1)

        @test out1[1] ≈ 10.0
        @test out1[2] ≈ 20.0
        @test out2[1] ≈ 5.0
        @test out2[2] ≈ 3.0
    end

end
