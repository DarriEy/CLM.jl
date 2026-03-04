@testset "DecompPrecisionControl" begin

    # ----------------------------------------------------------------
    # Helper: create minimal C and N state structures for testing
    # ----------------------------------------------------------------
    function make_precision_control_data(; nc=4, nlevdecomp=3, ndecomp_pools=7)
        cs = CLM.SoilBiogeochemCarbonStateData()
        ns = CLM.SoilBiogeochemNitrogenStateData()

        # Allocate key arrays
        cs.decomp_cpools_vr_col = zeros(nc, nlevdecomp, ndecomp_pools)
        cs.ctrunc_vr_col        = zeros(nc, nlevdecomp)

        ns.decomp_npools_vr_col = zeros(nc, nlevdecomp, ndecomp_pools)
        ns.ntrunc_vr_col        = zeros(nc, nlevdecomp)
        ns.smin_no3_vr_col      = zeros(nc, nlevdecomp)
        ns.smin_nh4_vr_col      = zeros(nc, nlevdecomp)

        mask = trues(nc)

        return cs, ns, mask
    end

    # ----------------------------------------------------------------
    # Constants
    # ----------------------------------------------------------------
    @testset "default thresholds" begin
        @test CLM.CCRIT_DEFAULT == 1.0e-8
        @test CLM.NCRIT_DEFAULT == 1.0e-8
    end

    # ----------------------------------------------------------------
    # Initialization
    # ----------------------------------------------------------------
    @testset "precision_control_init" begin
        cs = CLM.SoilBiogeochemCarbonStateData(totvegcthresh=NaN)
        ns = CLM.SoilBiogeochemNitrogenStateData(totvegcthresh=NaN)

        ccrit, ncrit = CLM.soil_bgc_precision_control_init!(cs, ns)
        @test ccrit == 1.0e-8
        @test ncrit == 1.0e-8
        @test cs.totvegcthresh == 1.0
        @test ns.totvegcthresh == 1.0
    end

    @testset "precision_control_init with isotopes" begin
        cs    = CLM.SoilBiogeochemCarbonStateData(totvegcthresh=NaN)
        c13cs = CLM.SoilBiogeochemCarbonStateData(totvegcthresh=NaN)
        c14cs = CLM.SoilBiogeochemCarbonStateData(totvegcthresh=NaN)
        ns    = CLM.SoilBiogeochemNitrogenStateData(totvegcthresh=NaN)

        ccrit, ncrit = CLM.soil_bgc_precision_control_init!(
            cs, ns;
            c13cs=c13cs, c14cs=c14cs,
            use_c13=true, use_c14=true,
            totvegcthresh=2.0)

        @test ccrit == 1.0e-8
        @test ncrit == 1.0e-8
        @test cs.totvegcthresh    == 2.0
        @test c13cs.totvegcthresh == 2.0
        @test c14cs.totvegcthresh == 2.0
        @test ns.totvegcthresh    == 2.0
    end

    @testset "precision_control_init skips isotopes when flags false" begin
        cs    = CLM.SoilBiogeochemCarbonStateData(totvegcthresh=NaN)
        c13cs = CLM.SoilBiogeochemCarbonStateData(totvegcthresh=NaN)
        c14cs = CLM.SoilBiogeochemCarbonStateData(totvegcthresh=NaN)
        ns    = CLM.SoilBiogeochemNitrogenStateData(totvegcthresh=NaN)

        CLM.soil_bgc_precision_control_init!(
            cs, ns;
            c13cs=c13cs, c14cs=c14cs,
            use_c13=false, use_c14=false)

        @test cs.totvegcthresh == 1.0
        @test isnan(c13cs.totvegcthresh)  # not set because use_c13=false
        @test isnan(c14cs.totvegcthresh)  # not set because use_c14=false
        @test ns.totvegcthresh == 1.0
    end

    # ----------------------------------------------------------------
    # Main precision control: basic truncation
    # ----------------------------------------------------------------
    @testset "truncation of small C pools" begin
        nc = 3; nlevdecomp = 2; ndecomp_pools = 4
        cs, ns, mask = make_precision_control_data(nc=nc, nlevdecomp=nlevdecomp,
                                                    ndecomp_pools=ndecomp_pools)

        # Set a few pools to tiny values below the threshold
        tiny = 5.0e-9  # below ccrit=1e-8
        cs.decomp_cpools_vr_col[1, 1, 1] = tiny
        cs.decomp_cpools_vr_col[1, 1, 2] = -tiny  # negative tiny
        cs.decomp_cpools_vr_col[1, 1, 3] = 1.0    # above threshold, should remain

        # Corresponding N pools
        ns.decomp_npools_vr_col[1, 1, 1] = tiny * 0.1
        ns.decomp_npools_vr_col[1, 1, 2] = -tiny * 0.1
        ns.decomp_npools_vr_col[1, 1, 3] = 0.1   # should remain

        CLM.soil_bgc_precision_control!(cs, ns;
            mask_bgc_soilc=mask,
            nlevdecomp=nlevdecomp,
            ndecomp_pools=ndecomp_pools)

        # Tiny pools should be zeroed
        @test cs.decomp_cpools_vr_col[1, 1, 1] == 0.0
        @test cs.decomp_cpools_vr_col[1, 1, 2] == 0.0
        @test ns.decomp_npools_vr_col[1, 1, 1] == 0.0
        @test ns.decomp_npools_vr_col[1, 1, 2] == 0.0

        # Large pool should be untouched
        @test cs.decomp_cpools_vr_col[1, 1, 3] == 1.0
        @test ns.decomp_npools_vr_col[1, 1, 3] == 0.1

        # Truncation sinks should accumulate the removed mass
        @test cs.ctrunc_vr_col[1, 1] == tiny + (-tiny)  # net zero for these two
        @test cs.ctrunc_vr_col[1, 1] == 0.0
        @test ns.ntrunc_vr_col[1, 1] == tiny * 0.1 + (-tiny * 0.1)
        @test ns.ntrunc_vr_col[1, 1] == 0.0
    end

    @testset "truncation accumulates into existing ctrunc" begin
        nc = 2; nlevdecomp = 1; ndecomp_pools = 2
        cs, ns, mask = make_precision_control_data(nc=nc, nlevdecomp=nlevdecomp,
                                                    ndecomp_pools=ndecomp_pools)

        # Pre-existing truncation value
        cs.ctrunc_vr_col[1, 1] = 1.0e-5
        ns.ntrunc_vr_col[1, 1] = 2.0e-5

        tiny = 3.0e-9
        cs.decomp_cpools_vr_col[1, 1, 1] = tiny
        ns.decomp_npools_vr_col[1, 1, 1] = tiny * 0.5

        CLM.soil_bgc_precision_control!(cs, ns;
            mask_bgc_soilc=mask,
            nlevdecomp=nlevdecomp,
            ndecomp_pools=ndecomp_pools)

        @test cs.decomp_cpools_vr_col[1, 1, 1] == 0.0
        @test cs.ctrunc_vr_col[1, 1] ≈ 1.0e-5 + tiny
        @test ns.ntrunc_vr_col[1, 1] ≈ 2.0e-5 + tiny * 0.5
    end

    @testset "values at exactly the threshold are NOT truncated" begin
        nc = 1; nlevdecomp = 1; ndecomp_pools = 1
        cs, ns, mask = make_precision_control_data(nc=nc, nlevdecomp=nlevdecomp,
                                                    ndecomp_pools=ndecomp_pools)

        # Set value exactly at ccrit -- abs(val) is NOT < ccrit, so it stays
        cs.decomp_cpools_vr_col[1, 1, 1] = 1.0e-8
        ns.decomp_npools_vr_col[1, 1, 1] = 1.0e-8

        CLM.soil_bgc_precision_control!(cs, ns;
            mask_bgc_soilc=mask,
            nlevdecomp=nlevdecomp,
            ndecomp_pools=ndecomp_pools)

        @test cs.decomp_cpools_vr_col[1, 1, 1] == 1.0e-8
        @test ns.decomp_npools_vr_col[1, 1, 1] == 1.0e-8
        @test cs.ctrunc_vr_col[1, 1] == 0.0
        @test ns.ntrunc_vr_col[1, 1] == 0.0
    end

    @testset "values just below threshold ARE truncated" begin
        nc = 1; nlevdecomp = 1; ndecomp_pools = 1
        cs, ns, mask = make_precision_control_data(nc=nc, nlevdecomp=nlevdecomp,
                                                    ndecomp_pools=ndecomp_pools)

        val = prevfloat(1.0e-8)  # just below 1e-8
        cs.decomp_cpools_vr_col[1, 1, 1] = val
        ns.decomp_npools_vr_col[1, 1, 1] = val * 0.1

        CLM.soil_bgc_precision_control!(cs, ns;
            mask_bgc_soilc=mask,
            nlevdecomp=nlevdecomp,
            ndecomp_pools=ndecomp_pools)

        @test cs.decomp_cpools_vr_col[1, 1, 1] == 0.0
        @test cs.ctrunc_vr_col[1, 1] == val
    end

    # ----------------------------------------------------------------
    # Mask filtering
    # ----------------------------------------------------------------
    @testset "mask excludes columns" begin
        nc = 3; nlevdecomp = 1; ndecomp_pools = 2
        cs, ns, mask = make_precision_control_data(nc=nc, nlevdecomp=nlevdecomp,
                                                    ndecomp_pools=ndecomp_pools)

        tiny = 5.0e-9
        # Set all columns to tiny values
        cs.decomp_cpools_vr_col .= tiny
        ns.decomp_npools_vr_col .= tiny * 0.1

        # Only column 2 is active
        mask .= false
        mask[2] = true

        CLM.soil_bgc_precision_control!(cs, ns;
            mask_bgc_soilc=mask,
            nlevdecomp=nlevdecomp,
            ndecomp_pools=ndecomp_pools)

        # Column 2 should be zeroed
        @test cs.decomp_cpools_vr_col[2, 1, 1] == 0.0
        @test cs.decomp_cpools_vr_col[2, 1, 2] == 0.0

        # Columns 1 and 3 should still have the tiny values
        @test cs.decomp_cpools_vr_col[1, 1, 1] == tiny
        @test cs.decomp_cpools_vr_col[3, 1, 1] == tiny
    end

    # ----------------------------------------------------------------
    # C13/C14 isotope truncation
    # ----------------------------------------------------------------
    @testset "C13 isotope truncation" begin
        nc = 2; nlevdecomp = 1; ndecomp_pools = 2
        cs, ns, mask = make_precision_control_data(nc=nc, nlevdecomp=nlevdecomp,
                                                    ndecomp_pools=ndecomp_pools)

        c13cs = CLM.SoilBiogeochemCarbonStateData()
        c13cs.decomp_cpools_vr_col = zeros(nc, nlevdecomp, ndecomp_pools)
        c13cs.ctrunc_vr_col        = zeros(nc, nlevdecomp)

        tiny = 2.0e-9
        cs.decomp_cpools_vr_col[1, 1, 1]    = tiny   # triggers truncation
        c13cs.decomp_cpools_vr_col[1, 1, 1] = tiny * 0.011  # C13 ratio
        ns.decomp_npools_vr_col[1, 1, 1]    = tiny * 0.1

        CLM.soil_bgc_precision_control!(cs, ns;
            mask_bgc_soilc=mask,
            nlevdecomp=nlevdecomp,
            ndecomp_pools=ndecomp_pools,
            c13cs=c13cs,
            use_c13=true)

        @test cs.decomp_cpools_vr_col[1, 1, 1]    == 0.0
        @test c13cs.decomp_cpools_vr_col[1, 1, 1] == 0.0
        @test c13cs.ctrunc_vr_col[1, 1] ≈ tiny * 0.011
    end

    @testset "C14 isotope truncation" begin
        nc = 2; nlevdecomp = 1; ndecomp_pools = 2
        cs, ns, mask = make_precision_control_data(nc=nc, nlevdecomp=nlevdecomp,
                                                    ndecomp_pools=ndecomp_pools)

        c14cs = CLM.SoilBiogeochemCarbonStateData()
        c14cs.decomp_cpools_vr_col = zeros(nc, nlevdecomp, ndecomp_pools)
        c14cs.ctrunc_vr_col        = zeros(nc, nlevdecomp)

        tiny = 4.0e-9
        cs.decomp_cpools_vr_col[1, 1, 2]    = tiny
        c14cs.decomp_cpools_vr_col[1, 1, 2] = tiny * 1.2e-12
        ns.decomp_npools_vr_col[1, 1, 2]    = tiny * 0.05

        CLM.soil_bgc_precision_control!(cs, ns;
            mask_bgc_soilc=mask,
            nlevdecomp=nlevdecomp,
            ndecomp_pools=ndecomp_pools,
            c14cs=c14cs,
            use_c14=true)

        @test cs.decomp_cpools_vr_col[1, 1, 2]    == 0.0
        @test c14cs.decomp_cpools_vr_col[1, 1, 2] == 0.0
        @test c14cs.ctrunc_vr_col[1, 1] ≈ tiny * 1.2e-12
    end

    @testset "C13 not truncated when use_c13 is false" begin
        nc = 1; nlevdecomp = 1; ndecomp_pools = 1
        cs, ns, mask = make_precision_control_data(nc=nc, nlevdecomp=nlevdecomp,
                                                    ndecomp_pools=ndecomp_pools)

        c13cs = CLM.SoilBiogeochemCarbonStateData()
        c13cs.decomp_cpools_vr_col = zeros(nc, nlevdecomp, ndecomp_pools)
        c13cs.ctrunc_vr_col        = zeros(nc, nlevdecomp)

        tiny = 3.0e-9
        cs.decomp_cpools_vr_col[1, 1, 1]    = tiny
        c13cs.decomp_cpools_vr_col[1, 1, 1] = tiny * 0.011
        ns.decomp_npools_vr_col[1, 1, 1]    = tiny * 0.1

        CLM.soil_bgc_precision_control!(cs, ns;
            mask_bgc_soilc=mask,
            nlevdecomp=nlevdecomp,
            ndecomp_pools=ndecomp_pools,
            c13cs=c13cs,
            use_c13=false)

        # C12 pool should still be truncated
        @test cs.decomp_cpools_vr_col[1, 1, 1] == 0.0
        # C13 should NOT be touched
        @test c13cs.decomp_cpools_vr_col[1, 1, 1] == tiny * 0.011
        @test c13cs.ctrunc_vr_col[1, 1] == 0.0
    end

    # ----------------------------------------------------------------
    # Nitrification/denitrification mineral N stability
    # ----------------------------------------------------------------
    @testset "small negative NO3 reset to zero" begin
        nc = 2; nlevdecomp = 2; ndecomp_pools = 1
        cs, ns, mask = make_precision_control_data(nc=nc, nlevdecomp=nlevdecomp,
                                                    ndecomp_pools=ndecomp_pools)

        # Set small negative NO3 value (below ncrit/1e4 = 1e-12)
        ns.smin_no3_vr_col[1, 1] = -5.0e-13

        CLM.soil_bgc_precision_control!(cs, ns;
            mask_bgc_soilc=mask,
            nlevdecomp=nlevdecomp,
            ndecomp_pools=ndecomp_pools,
            use_nitrif_denitrif=true,
            use_fun=false)

        @test ns.smin_no3_vr_col[1, 1] == 0.0
    end

    @testset "small negative NH4 reset to zero" begin
        nc = 2; nlevdecomp = 2; ndecomp_pools = 1
        cs, ns, mask = make_precision_control_data(nc=nc, nlevdecomp=nlevdecomp,
                                                    ndecomp_pools=ndecomp_pools)

        ns.smin_nh4_vr_col[1, 2] = -3.0e-13

        CLM.soil_bgc_precision_control!(cs, ns;
            mask_bgc_soilc=mask,
            nlevdecomp=nlevdecomp,
            ndecomp_pools=ndecomp_pools,
            use_nitrif_denitrif=true,
            use_fun=false)

        @test ns.smin_nh4_vr_col[1, 2] == 0.0
    end

    @testset "small positive NO3/NH4 NOT reset" begin
        nc = 1; nlevdecomp = 1; ndecomp_pools = 1
        cs, ns, mask = make_precision_control_data(nc=nc, nlevdecomp=nlevdecomp,
                                                    ndecomp_pools=ndecomp_pools)

        # Small positive values should remain (they are small but not negative)
        ns.smin_no3_vr_col[1, 1] = 5.0e-13
        ns.smin_nh4_vr_col[1, 1] = 3.0e-13

        CLM.soil_bgc_precision_control!(cs, ns;
            mask_bgc_soilc=mask,
            nlevdecomp=nlevdecomp,
            ndecomp_pools=ndecomp_pools,
            use_nitrif_denitrif=true,
            use_fun=false)

        @test ns.smin_no3_vr_col[1, 1] == 5.0e-13
        @test ns.smin_nh4_vr_col[1, 1] == 3.0e-13
    end

    @testset "NO3/NH4 not touched when use_fun is true" begin
        nc = 1; nlevdecomp = 1; ndecomp_pools = 1
        cs, ns, mask = make_precision_control_data(nc=nc, nlevdecomp=nlevdecomp,
                                                    ndecomp_pools=ndecomp_pools)

        ns.smin_no3_vr_col[1, 1] = -5.0e-13
        ns.smin_nh4_vr_col[1, 1] = -3.0e-13

        CLM.soil_bgc_precision_control!(cs, ns;
            mask_bgc_soilc=mask,
            nlevdecomp=nlevdecomp,
            ndecomp_pools=ndecomp_pools,
            use_nitrif_denitrif=true,
            use_fun=true)

        # Should remain negative because use_fun=true skips this cleanup
        @test ns.smin_no3_vr_col[1, 1] == -5.0e-13
        @test ns.smin_nh4_vr_col[1, 1] == -3.0e-13
    end

    @testset "NO3/NH4 not touched when use_nitrif_denitrif is false" begin
        nc = 1; nlevdecomp = 1; ndecomp_pools = 1
        cs, ns, mask = make_precision_control_data(nc=nc, nlevdecomp=nlevdecomp,
                                                    ndecomp_pools=ndecomp_pools)

        ns.smin_no3_vr_col[1, 1] = -5.0e-13
        ns.smin_nh4_vr_col[1, 1] = -3.0e-13

        CLM.soil_bgc_precision_control!(cs, ns;
            mask_bgc_soilc=mask,
            nlevdecomp=nlevdecomp,
            ndecomp_pools=ndecomp_pools,
            use_nitrif_denitrif=false,
            use_fun=false)

        @test ns.smin_no3_vr_col[1, 1] == -5.0e-13
        @test ns.smin_nh4_vr_col[1, 1] == -3.0e-13
    end

    @testset "larger negative NO3/NH4 NOT reset (above ncrit/1e4 threshold)" begin
        nc = 1; nlevdecomp = 1; ndecomp_pools = 1
        cs, ns, mask = make_precision_control_data(nc=nc, nlevdecomp=nlevdecomp,
                                                    ndecomp_pools=ndecomp_pools)

        # These are negative but abs > ncrit/1e4 = 1e-12, so not reset
        ns.smin_no3_vr_col[1, 1] = -5.0e-8
        ns.smin_nh4_vr_col[1, 1] = -3.0e-8

        CLM.soil_bgc_precision_control!(cs, ns;
            mask_bgc_soilc=mask,
            nlevdecomp=nlevdecomp,
            ndecomp_pools=ndecomp_pools,
            use_nitrif_denitrif=true,
            use_fun=false)

        @test ns.smin_no3_vr_col[1, 1] == -5.0e-8
        @test ns.smin_nh4_vr_col[1, 1] == -3.0e-8
    end

    # ----------------------------------------------------------------
    # Multi-column, multi-level, multi-pool integration test
    # ----------------------------------------------------------------
    @testset "multi-column multi-level integration" begin
        nc = 3; nlevdecomp = 4; ndecomp_pools = 5
        cs, ns, mask = make_precision_control_data(nc=nc, nlevdecomp=nlevdecomp,
                                                    ndecomp_pools=ndecomp_pools)

        # Set various values -- some tiny, some large
        cs.decomp_cpools_vr_col[1, 1, 1] = 1.0e-10   # tiny, truncate
        cs.decomp_cpools_vr_col[1, 2, 3] = -2.0e-9    # tiny negative, truncate
        cs.decomp_cpools_vr_col[2, 3, 5] = 5.0e-9     # tiny, truncate
        cs.decomp_cpools_vr_col[3, 4, 2] = 1.0        # large, keep
        cs.decomp_cpools_vr_col[1, 1, 2] = 0.5        # large, keep

        ns.decomp_npools_vr_col[1, 1, 1] = 1.0e-11
        ns.decomp_npools_vr_col[1, 2, 3] = -3.0e-10
        ns.decomp_npools_vr_col[2, 3, 5] = 7.0e-10
        ns.decomp_npools_vr_col[3, 4, 2] = 0.1
        ns.decomp_npools_vr_col[1, 1, 2] = 0.05

        CLM.soil_bgc_precision_control!(cs, ns;
            mask_bgc_soilc=mask,
            nlevdecomp=nlevdecomp,
            ndecomp_pools=ndecomp_pools)

        # Truncated pools
        @test cs.decomp_cpools_vr_col[1, 1, 1] == 0.0
        @test cs.decomp_cpools_vr_col[1, 2, 3] == 0.0
        @test cs.decomp_cpools_vr_col[2, 3, 5] == 0.0

        # Kept pools
        @test cs.decomp_cpools_vr_col[3, 4, 2] == 1.0
        @test cs.decomp_cpools_vr_col[1, 1, 2] == 0.5

        # Truncation sinks
        @test cs.ctrunc_vr_col[1, 1] ≈ 1.0e-10
        @test cs.ctrunc_vr_col[1, 2] ≈ -2.0e-9
        @test cs.ctrunc_vr_col[2, 3] ≈ 5.0e-9

        @test ns.ntrunc_vr_col[1, 1] ≈ 1.0e-11
        @test ns.ntrunc_vr_col[1, 2] ≈ -3.0e-10
        @test ns.ntrunc_vr_col[2, 3] ≈ 7.0e-10

        # N pools for kept C should also be kept
        @test ns.decomp_npools_vr_col[3, 4, 2] == 0.1
        @test ns.decomp_npools_vr_col[1, 1, 2] == 0.05
    end

    # ----------------------------------------------------------------
    # Mass conservation test
    # ----------------------------------------------------------------
    @testset "mass conservation (C pool + ctrunc)" begin
        nc = 2; nlevdecomp = 3; ndecomp_pools = 4
        cs, ns, mask = make_precision_control_data(nc=nc, nlevdecomp=nlevdecomp,
                                                    ndecomp_pools=ndecomp_pools)

        # Fill with random small values (some above, some below threshold)
        for c in 1:nc, j in 1:nlevdecomp, k in 1:ndecomp_pools
            cs.decomp_cpools_vr_col[c, j, k] = (rand() - 0.5) * 2.0e-8
            ns.decomp_npools_vr_col[c, j, k] = (rand() - 0.5) * 2.0e-8
        end

        # Compute total mass before (pools + truncation sinks)
        total_c_before = sum(cs.decomp_cpools_vr_col) + sum(cs.ctrunc_vr_col)
        total_n_before = sum(ns.decomp_npools_vr_col) + sum(ns.ntrunc_vr_col)

        CLM.soil_bgc_precision_control!(cs, ns;
            mask_bgc_soilc=mask,
            nlevdecomp=nlevdecomp,
            ndecomp_pools=ndecomp_pools)

        total_c_after = sum(cs.decomp_cpools_vr_col) + sum(cs.ctrunc_vr_col)
        total_n_after = sum(ns.decomp_npools_vr_col) + sum(ns.ntrunc_vr_col)

        @test total_c_after ≈ total_c_before atol=1.0e-20
        @test total_n_after ≈ total_n_before atol=1.0e-20
    end

    # ----------------------------------------------------------------
    # Custom ccrit/ncrit
    # ----------------------------------------------------------------
    @testset "custom ccrit threshold" begin
        nc = 1; nlevdecomp = 1; ndecomp_pools = 1
        cs, ns, mask = make_precision_control_data(nc=nc, nlevdecomp=nlevdecomp,
                                                    ndecomp_pools=ndecomp_pools)

        # Use a larger threshold
        cs.decomp_cpools_vr_col[1, 1, 1] = 0.5e-4  # below custom ccrit of 1e-4
        ns.decomp_npools_vr_col[1, 1, 1] = 0.5e-5

        CLM.soil_bgc_precision_control!(cs, ns;
            mask_bgc_soilc=mask,
            nlevdecomp=nlevdecomp,
            ndecomp_pools=ndecomp_pools,
            ccrit=1.0e-4)

        @test cs.decomp_cpools_vr_col[1, 1, 1] == 0.0
        @test cs.ctrunc_vr_col[1, 1] ≈ 0.5e-4
    end

    # ----------------------------------------------------------------
    # All zeros -- no changes expected
    # ----------------------------------------------------------------
    @testset "all-zero pools unchanged" begin
        nc = 2; nlevdecomp = 2; ndecomp_pools = 3
        cs, ns, mask = make_precision_control_data(nc=nc, nlevdecomp=nlevdecomp,
                                                    ndecomp_pools=ndecomp_pools)

        # Everything starts at zero
        CLM.soil_bgc_precision_control!(cs, ns;
            mask_bgc_soilc=mask,
            nlevdecomp=nlevdecomp,
            ndecomp_pools=ndecomp_pools)

        @test all(cs.decomp_cpools_vr_col .== 0.0)
        @test all(ns.decomp_npools_vr_col .== 0.0)
        @test all(cs.ctrunc_vr_col .== 0.0)
        @test all(ns.ntrunc_vr_col .== 0.0)
    end

    # ----------------------------------------------------------------
    # Empty mask -- nothing should change
    # ----------------------------------------------------------------
    @testset "empty mask changes nothing" begin
        nc = 2; nlevdecomp = 1; ndecomp_pools = 2
        cs, ns, mask = make_precision_control_data(nc=nc, nlevdecomp=nlevdecomp,
                                                    ndecomp_pools=ndecomp_pools)

        tiny = 3.0e-9
        cs.decomp_cpools_vr_col .= tiny
        ns.decomp_npools_vr_col .= tiny

        mask .= false  # no active columns

        CLM.soil_bgc_precision_control!(cs, ns;
            mask_bgc_soilc=mask,
            nlevdecomp=nlevdecomp,
            ndecomp_pools=ndecomp_pools)

        # All values should remain unchanged
        @test all(cs.decomp_cpools_vr_col .== tiny)
        @test all(ns.decomp_npools_vr_col .== tiny)
    end

    # ----------------------------------------------------------------
    # N pool truncation tied to C decision
    # ----------------------------------------------------------------
    @testset "N pool zeroed only when C pool triggers truncation" begin
        nc = 1; nlevdecomp = 1; ndecomp_pools = 2
        cs, ns, mask = make_precision_control_data(nc=nc, nlevdecomp=nlevdecomp,
                                                    ndecomp_pools=ndecomp_pools)

        # Pool 1: C is tiny (triggers truncation), N is large
        cs.decomp_cpools_vr_col[1, 1, 1] = 5.0e-9
        ns.decomp_npools_vr_col[1, 1, 1] = 100.0  # large N value

        # Pool 2: C is large (no truncation), N is tiny
        cs.decomp_cpools_vr_col[1, 1, 2] = 100.0
        ns.decomp_npools_vr_col[1, 1, 2] = 5.0e-9  # tiny N

        CLM.soil_bgc_precision_control!(cs, ns;
            mask_bgc_soilc=mask,
            nlevdecomp=nlevdecomp,
            ndecomp_pools=ndecomp_pools)

        # Pool 1: both C and N zeroed because C triggered
        @test cs.decomp_cpools_vr_col[1, 1, 1] == 0.0
        @test ns.decomp_npools_vr_col[1, 1, 1] == 0.0
        @test ns.ntrunc_vr_col[1, 1] ≈ 100.0

        # Pool 2: neither changed because C did not trigger
        @test cs.decomp_cpools_vr_col[1, 1, 2] == 100.0
        @test ns.decomp_npools_vr_col[1, 1, 2] == 5.0e-9
    end

end
