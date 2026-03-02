@testset "N Dynamics" begin

    # ----------------------------------------------------------------
    # Helper: create minimal data structures for testing
    # ----------------------------------------------------------------
    function make_n_dynamics_data(; nc=4, np=8, ng=2)
        ndecomp_pools = 7
        ndecomp_cascade_transitions = 10
        nlevdecomp = 1
        nlevsoi = 3

        # --- Soil biogeochem N flux ---
        nf = CLM.SoilBiogeochemNitrogenFluxData()
        CLM.soil_bgc_nitrogen_flux_init!(nf, nc, nlevsoi, ndecomp_pools, ndecomp_cascade_transitions)

        # --- Soil biogeochem N state ---
        ns = CLM.SoilBiogeochemNitrogenStateData()
        CLM.soil_bgc_nitrogen_state_init!(ns, nc, 1, nlevsoi, ndecomp_pools)
        ns.sminn_col .= 20.0  # gN/m2

        # --- Soil biogeochem state ---
        soilbgc_st = CLM.SoilBiogeochemStateData()
        CLM.soil_bgc_state_init!(soilbgc_st, nc, np, nlevsoi, ndecomp_cascade_transitions)
        soilbgc_st.fpg_col .= 0.5

        # --- CNVeg carbon flux ---
        cf = CLM.CNVegCarbonFluxData()
        cf.annsum_npp_col = fill(200.0, nc)    # gC/m2/yr
        cf.lag_npp_col    = fill(5.0e-6, nc)   # gC/m2/s

        # --- CNVeg nitrogen flux ---
        cnveg_nf = CLM.CNVegNitrogenFluxData()
        cnveg_nf.fert_patch          = fill(1.0e-7, np)   # gN/m2/s
        cnveg_nf.soyfixn_patch       = zeros(np)
        cnveg_nf.plant_ndemand_patch = fill(2.0e-7, np)   # gN/m2/s

        # --- CNVeg state ---
        cnveg_state = CLM.CNVegStateData()
        cnveg_state.gddmaturity_patch = fill(1500.0, np)

        # --- Crop ---
        crop = CLM.CropData()
        crop.hui_patch      = fill(600.0, np)    # growing degree-days
        crop.croplive_patch = fill(true, np)

        # --- Water diagnostic bulk ---
        wdiag = CLM.WaterDiagnosticBulkData()
        wdiag.wf_col = fill(0.6, nc)

        # --- Water flux bulk ---
        wfb = CLM.WaterFluxBulkData()
        wfb.AnnET = fill(2.0e-5, nc)  # mm H2O/s

        # --- Column data (set only fields needed for N dynamics) ---
        col = CLM.ColumnData()
        col.gridcell = [1, 1, 2, 2]
        col.is_fates = fill(false, nc)

        # --- Patch data (set only fields needed for N dynamics) ---
        patch = CLM.PatchData()
        patch.column = [1, 1, 2, 2, 3, 3, 4, 4]
        patch.wtcol  = fill(0.5, np)   # equal weighting
        patch.itype  = fill(1, np)     # default non-soybean type

        # --- Masks and bounds ---
        mask_soilc = trues(nc)
        mask_soilp = trues(np)

        # --- N dynamics params ---
        params = CLM.NDynamicsParams()

        return (nf=nf, ns=ns, soilbgc_st=soilbgc_st,
                cf=cf, cnveg_nf=cnveg_nf, cnveg_state=cnveg_state,
                crop=crop, wdiag=wdiag, wfb=wfb,
                col=col, patch=patch, params=params,
                mask_soilc=mask_soilc, mask_soilp=mask_soilp,
                nc=nc, np=np, ng=ng)
    end

    # ================================================================
    # Test 1: NDynamicsParams construction and read
    # ================================================================
    @testset "NDynamicsParams construction" begin
        params = CLM.NDynamicsParams()
        @test params.freelivfix_intercept == 0.0117
        @test params.freelivfix_slope_wET == 0.0006
    end

    @testset "n_dynamics_read_nml!" begin
        params = CLM.NDynamicsParams()
        CLM.n_dynamics_read_nml!(params;
            freelivfix_intercept=0.02,
            freelivfix_slope_wET=0.001)
        @test params.freelivfix_intercept == 0.02
        @test params.freelivfix_slope_wET == 0.001
    end

    # ================================================================
    # Test 2: N deposition
    # ================================================================
    @testset "n_deposition!" begin
        d = make_n_dynamics_data()
        forc_ndep = [1.0e-8, 2.0e-8]  # gN/m2/s per gridcell

        CLM.n_deposition!(d.nf;
            forc_ndep=forc_ndep,
            col_gridcell=d.col.gridcell,
            bounds=1:d.nc)

        # Columns 1,2 → gridcell 1; columns 3,4 → gridcell 2
        @test d.nf.ndep_to_sminn_col[1] == 1.0e-8
        @test d.nf.ndep_to_sminn_col[2] == 1.0e-8
        @test d.nf.ndep_to_sminn_col[3] == 2.0e-8
        @test d.nf.ndep_to_sminn_col[4] == 2.0e-8
    end

    # ================================================================
    # Test 3: N free living fixation
    # ================================================================
    @testset "n_free_living_fixation!" begin
        d = make_n_dynamics_data()
        dayspyr = 365.0

        CLM.n_free_living_fixation!(d.nf, d.params;
            mask_soilc=d.mask_soilc,
            bounds=1:d.nc,
            AnnET=d.wfb.AnnET,
            dayspyr=dayspyr)

        secs_per_year = dayspyr * 24.0 * 3600.0
        expected = (d.params.freelivfix_slope_wET * (max(0.0, d.wfb.AnnET[1]) * secs_per_year) +
                    d.params.freelivfix_intercept) / secs_per_year

        for c in 1:d.nc
            @test d.nf.ffix_to_sminn_col[c] ≈ expected
            @test d.nf.ffix_to_sminn_col[c] > 0.0
        end
    end

    @testset "n_free_living_fixation! with negative ET" begin
        d = make_n_dynamics_data()
        d.wfb.AnnET .= -1.0  # negative ET

        CLM.n_free_living_fixation!(d.nf, d.params;
            mask_soilc=d.mask_soilc,
            bounds=1:d.nc,
            AnnET=d.wfb.AnnET,
            dayspyr=365.0)

        secs_per_year = 365.0 * 24.0 * 3600.0
        # With max(0,AnnET)=0, only intercept contributes
        expected = d.params.freelivfix_intercept / secs_per_year

        for c in 1:d.nc
            @test d.nf.ffix_to_sminn_col[c] ≈ expected
        end
    end

    # ================================================================
    # Test 4: N fixation (exponential relaxation mode)
    # ================================================================
    @testset "n_fixation! (relaxation mode)" begin
        d = make_n_dynamics_data()
        dayspyr = 365.0

        CLM.n_fixation!(d.nf, d.cf;
            mask_soilc=d.mask_soilc,
            bounds=1:d.nc,
            col_is_fates=d.col.is_fates,
            dayspyr=dayspyr,
            nfix_timeconst=10.0,  # triggers relaxation mode
            use_fun=false)

        for c in 1:d.nc
            @test d.nf.nfix_to_sminn_col[c] >= 0.0
            @test isfinite(d.nf.nfix_to_sminn_col[c])
        end

        # Verify against hand-calculation
        npp = d.cf.lag_npp_col[1]
        t_expected = (1.8 * (1.0 - exp(-0.003 * npp * (CLM.SECSPDAY * dayspyr)))) /
                     (CLM.SECSPDAY * dayspyr)
        @test d.nf.nfix_to_sminn_col[1] ≈ max(0.0, t_expected)
    end

    # ================================================================
    # Test 5: N fixation (annual mean mode)
    # ================================================================
    @testset "n_fixation! (annual mean mode)" begin
        d = make_n_dynamics_data()
        dayspyr = 365.0

        CLM.n_fixation!(d.nf, d.cf;
            mask_soilc=d.mask_soilc,
            bounds=1:d.nc,
            col_is_fates=d.col.is_fates,
            dayspyr=dayspyr,
            nfix_timeconst=-1.0,  # triggers annual mean mode
            use_fun=false)

        for c in 1:d.nc
            @test d.nf.nfix_to_sminn_col[c] >= 0.0
            @test isfinite(d.nf.nfix_to_sminn_col[c])
        end

        # Verify against hand-calculation
        npp = d.cf.annsum_npp_col[1]
        t_expected = (1.8 * (1.0 - exp(-0.003 * npp))) / (CLM.SECSPDAY * dayspyr)
        @test d.nf.nfix_to_sminn_col[1] ≈ max(0.0, t_expected)
    end

    # ================================================================
    # Test 6: N fixation with use_fun=true sets to zero
    # ================================================================
    @testset "n_fixation! (use_fun=true)" begin
        d = make_n_dynamics_data()

        CLM.n_fixation!(d.nf, d.cf;
            mask_soilc=d.mask_soilc,
            bounds=1:d.nc,
            col_is_fates=d.col.is_fates,
            dayspyr=365.0,
            nfix_timeconst=-1.0,
            use_fun=true)

        for c in 1:d.nc
            @test d.nf.nfix_to_sminn_col[c] == 0.0
        end
    end

    # ================================================================
    # Test 7: N fixation with FATES columns
    # ================================================================
    @testset "n_fixation! (FATES columns)" begin
        d = make_n_dynamics_data()
        d.col.is_fates .= true  # all columns are FATES

        CLM.n_fixation!(d.nf, d.cf;
            mask_soilc=d.mask_soilc,
            bounds=1:d.nc,
            col_is_fates=d.col.is_fates,
            dayspyr=365.0,
            nfix_timeconst=-1.0,
            use_fun=false)

        # With FATES, npp=0 → t = 1.8*(1 - exp(0))/(sec/yr) = 0
        for c in 1:d.nc
            @test d.nf.nfix_to_sminn_col[c] == 0.0
        end
    end

    # ================================================================
    # Test 8: N fertilizer (p2c)
    # ================================================================
    @testset "n_fert!" begin
        d = make_n_dynamics_data()

        # Set distinct fert values per patch
        d.cnveg_nf.fert_patch .= [1.0e-7, 3.0e-7,   # col 1 patches
                                   2.0e-7, 4.0e-7,   # col 2 patches
                                   0.0, 0.0,          # col 3 patches
                                   5.0e-7, 1.0e-7]    # col 4 patches

        CLM.n_fert!(d.nf, d.cnveg_nf;
            mask_soilc=d.mask_soilc,
            bounds=1:d.nc,
            patch=d.patch,
            mask_soilp=d.mask_soilp,
            bounds_p=1:d.np)

        # Column 1: (1e-7 * 0.5) + (3e-7 * 0.5) = 2e-7
        @test d.nf.fert_to_sminn_col[1] ≈ 2.0e-7
        # Column 2: (2e-7 * 0.5) + (4e-7 * 0.5) = 3e-7
        @test d.nf.fert_to_sminn_col[2] ≈ 3.0e-7
        # Column 3: 0
        @test d.nf.fert_to_sminn_col[3] ≈ 0.0
        # Column 4: (5e-7 * 0.5) + (1e-7 * 0.5) = 3e-7
        @test d.nf.fert_to_sminn_col[4] ≈ 3.0e-7
    end

    # ================================================================
    # Test 9: Soy fixation — active soybean
    # ================================================================
    @testset "n_soyfix! (active soybean)" begin
        d = make_n_dynamics_data()

        ntmp_soybean       = 23
        nirrig_tmp_soybean = 58
        ntrp_soybean       = 24
        nirrig_trp_soybean = 59

        # Set patches 1,2 (col 1) to soybean type
        d.patch.itype[1] = ntmp_soybean
        d.patch.itype[2] = ntmp_soybean
        d.crop.croplive_patch[1] = true
        d.crop.croplive_patch[2] = true

        # Set GDD fraction to be in the active range (0.30-0.55 → fxg=1.0)
        d.crop.hui_patch[1] = 600.0
        d.crop.hui_patch[2] = 600.0
        d.cnveg_state.gddmaturity_patch[1] = 1500.0
        d.cnveg_state.gddmaturity_patch[2] = 1500.0
        # GDDfrac = 600/1500 = 0.4 (in [0.30, 0.55] → fxg = 1.0)

        # fpg < 1 means there's unmet N demand
        d.soilbgc_st.fpg_col[1] = 0.5

        # sminn below threshold2 → fxn = 1.0
        d.ns.sminn_col[1] = 5.0

        # wf/0.85 = 0.6/0.85 ≈ 0.706 → fxw ≈ 0.706
        d.wdiag.wf_col[1] = 0.6

        CLM.n_soyfix!(d.nf, d.cnveg_nf, d.soilbgc_st, d.ns,
            d.cnveg_state, d.crop, d.wdiag;
            mask_soilc=d.mask_soilc,
            bounds=1:d.nc,
            mask_soilp=d.mask_soilp,
            bounds_p=1:d.np,
            patch=d.patch,
            ntmp_soybean=ntmp_soybean,
            nirrig_tmp_soybean=nirrig_tmp_soybean,
            ntrp_soybean=ntrp_soybean,
            nirrig_trp_soybean=nirrig_trp_soybean)

        # Soybean patches should have positive fixation
        @test d.cnveg_nf.soyfixn_patch[1] > 0.0
        @test d.cnveg_nf.soyfixn_patch[2] > 0.0

        # Non-soybean patches should have zero
        @test d.cnveg_nf.soyfixn_patch[3] == 0.0

        # Column-level should be positive
        @test d.nf.soyfixn_to_sminn_col[1] > 0.0
    end

    # ================================================================
    # Test 10: Soy fixation — non-soybean patches
    # ================================================================
    @testset "n_soyfix! (non-soybean)" begin
        d = make_n_dynamics_data()
        # All patches are non-soybean (default itype=1)

        CLM.n_soyfix!(d.nf, d.cnveg_nf, d.soilbgc_st, d.ns,
            d.cnveg_state, d.crop, d.wdiag;
            mask_soilc=d.mask_soilc,
            bounds=1:d.nc,
            mask_soilp=d.mask_soilp,
            bounds_p=1:d.np,
            patch=d.patch,
            ntmp_soybean=23,
            nirrig_tmp_soybean=58,
            ntrp_soybean=24,
            nirrig_trp_soybean=59)

        for p in 1:d.np
            @test d.cnveg_nf.soyfixn_patch[p] == 0.0
        end
        for c in 1:d.nc
            @test d.nf.soyfixn_to_sminn_col[c] == 0.0
        end
    end

    # ================================================================
    # Test 11: Soy fixation — high soil N suppresses fixation
    # ================================================================
    @testset "n_soyfix! (high soil N → fxn=0)" begin
        d = make_n_dynamics_data()

        ntmp_soybean = 23
        d.patch.itype[1] = ntmp_soybean
        d.crop.croplive_patch[1] = true
        d.soilbgc_st.fpg_col[1] = 0.5

        # sminn > 30 → fxn = 0, which zeroes out fixation
        d.ns.sminn_col[1] = 50.0

        # GDDfrac in active range
        d.crop.hui_patch[1] = 600.0
        d.cnveg_state.gddmaturity_patch[1] = 1500.0

        CLM.n_soyfix!(d.nf, d.cnveg_nf, d.soilbgc_st, d.ns,
            d.cnveg_state, d.crop, d.wdiag;
            mask_soilc=d.mask_soilc,
            bounds=1:d.nc,
            mask_soilp=d.mask_soilp,
            bounds_p=1:d.np,
            patch=d.patch,
            ntmp_soybean=ntmp_soybean,
            nirrig_tmp_soybean=58,
            ntrp_soybean=24,
            nirrig_trp_soybean=59)

        # fxn=0 → fxr=0 → soyfixn=0
        @test d.cnveg_nf.soyfixn_patch[1] == 0.0
    end

    # ================================================================
    # Test 12: Soy fixation — fpg>=1 means demand met, no fixation
    # ================================================================
    @testset "n_soyfix! (demand met, fpg>=1)" begin
        d = make_n_dynamics_data()

        ntmp_soybean = 23
        d.patch.itype[1] = ntmp_soybean
        d.crop.croplive_patch[1] = true
        d.soilbgc_st.fpg_col[1] = 1.0  # demand fully met

        d.crop.hui_patch[1] = 600.0
        d.cnveg_state.gddmaturity_patch[1] = 1500.0
        d.ns.sminn_col[1] = 5.0

        CLM.n_soyfix!(d.nf, d.cnveg_nf, d.soilbgc_st, d.ns,
            d.cnveg_state, d.crop, d.wdiag;
            mask_soilc=d.mask_soilc,
            bounds=1:d.nc,
            mask_soilp=d.mask_soilp,
            bounds_p=1:d.np,
            patch=d.patch,
            ntmp_soybean=ntmp_soybean,
            nirrig_tmp_soybean=58,
            ntrp_soybean=24,
            nirrig_trp_soybean=59)

        @test d.cnveg_nf.soyfixn_patch[1] == 0.0
    end

    # ================================================================
    # Test 13: Mask filtering
    # ================================================================
    @testset "Mask filtering" begin
        d = make_n_dynamics_data()
        d.mask_soilc[2] = false
        d.mask_soilc[4] = false

        CLM.n_free_living_fixation!(d.nf, d.params;
            mask_soilc=d.mask_soilc,
            bounds=1:d.nc,
            AnnET=d.wfb.AnnET,
            dayspyr=365.0)

        # Active columns should have computed values
        @test d.nf.ffix_to_sminn_col[1] > 0.0
        @test d.nf.ffix_to_sminn_col[3] > 0.0

        # Masked columns should still have NaN (not updated)
        @test isnan(d.nf.ffix_to_sminn_col[2])
        @test isnan(d.nf.ffix_to_sminn_col[4])
    end

end
