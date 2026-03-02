@testset "Phenology" begin

    # =====================================================================
    # Constants
    # =====================================================================
    @testset "Phenology Constants" begin
        @test CLM.NOT_Planted == 999
        @test CLM.NOT_Harvested == 999
        @test CLM.inNH == 1
        @test CLM.inSH == 2

        @test CLM.critical_daylight_constant == 1
        @test CLM.critical_daylight_depends_on_lat == 2
        @test CLM.critical_daylight_depends_on_veg == 3
        @test CLM.critical_daylight_depends_on_latnveg == 4

        @test CLM.critical_offset_high_lat == 65.0

        @test CLM.HARVEST_REASON_MATURE == 1.0
        @test CLM.HARVEST_REASON_MAXSEASLENGTH == 2.0
        @test CLM.HARVEST_REASON_SOWNBADDEC31 == 3.0
        @test CLM.HARVEST_REASON_SOWTODAY == 4.0
        @test CLM.HARVEST_REASON_SOWTOMORROW == 5.0
        @test CLM.HARVEST_REASON_IDOPTOMORROW == 6.0
        @test CLM.HARVEST_REASON_VERNFREEZEKILL == 7.0
    end

    # =====================================================================
    # PhenologyParams struct
    # =====================================================================
    @testset "PhenologyParams defaults" begin
        pp = CLM.PhenologyParams()
        @test pp.crit_dayl ≈ 39200.0
        @test pp.crit_dayl_at_high_lat ≈ 54000.0
        @test pp.crit_dayl_lat_slope ≈ 720.0
        @test pp.ndays_off ≈ 30.0
        @test pp.fstor2tran ≈ 0.5
        @test pp.crit_onset_fdd ≈ 15.0
        @test pp.crit_onset_swi ≈ 15.0
        @test pp.soilpsi_on ≈ -0.6
        @test pp.crit_offset_fdd ≈ 15.0
        @test pp.crit_offset_swi ≈ 15.0
        @test pp.soilpsi_off ≈ -0.8
        @test pp.lwtop ≈ 0.7
        @test pp.phenology_soil_depth ≈ 0.08
        @test pp.snow5d_thresh_for_onset ≈ 0.2
    end

    # =====================================================================
    # PhenologyState struct
    # =====================================================================
    @testset "PhenologyState defaults" begin
        ps = CLM.PhenologyState()
        @test ps.dt == 0.0
        @test ps.fracday == 0.0
        @test ps.phenology_soil_layer == 1
        @test ps.p1d ≈ 0.004
        @test ps.p1v ≈ 0.003
        @test ps.hti ≈ 1.0
        @test ps.tbase ≈ 0.0
        @test ps.jdayyrstart == [1, 182]
    end

    # =====================================================================
    # PftConPhenology struct
    # =====================================================================
    @testset "PftConPhenology defaults" begin
        pfc = CLM.PftConPhenology()
        @test isempty(pfc.evergreen)
        @test isempty(pfc.season_decid)
        @test isempty(pfc.stress_decid)
        @test isempty(pfc.woody)
        @test isempty(pfc.leaf_long)
        @test isempty(pfc.leafcn)
    end

    # =====================================================================
    # cn_phenology_set_params!
    # =====================================================================
    @testset "cn_phenology_set_params!" begin
        pp = CLM.PhenologyParams(crit_dayl=0.0, ndays_off=0.0)
        CLM.cn_phenology_set_params!(pp)
        @test pp.crit_dayl ≈ 39200.0
        @test pp.ndays_off ≈ 30.0
        @test pp.fstor2tran ≈ 0.5
        @test pp.lwtop ≈ 0.7
    end

    # =====================================================================
    # cn_phenology_set_nml!
    # =====================================================================
    @testset "cn_phenology_set_nml!" begin
        CLM.cn_phenology_set_nml!(onset_thresh_depends_on_veg=false,
                                   critical_daylight_method_in=CLM.critical_daylight_constant)
        @test CLM._onset_thresh_depends_on_veg[] == false
        @test CLM._critical_daylight_method[] == CLM.critical_daylight_constant

        CLM.cn_phenology_set_nml!(onset_thresh_depends_on_veg=true,
                                   critical_daylight_method_in=CLM.critical_daylight_depends_on_lat)
        @test CLM._onset_thresh_depends_on_veg[] == true
        @test CLM._critical_daylight_method[] == CLM.critical_daylight_depends_on_lat

        # Out-of-range method
        @test_throws ErrorException CLM.cn_phenology_set_nml!(critical_daylight_method_in=0)
        @test_throws ErrorException CLM.cn_phenology_set_nml!(critical_daylight_method_in=5)

        # Reset
        CLM.cn_phenology_set_nml!(onset_thresh_depends_on_veg=false,
                                   critical_daylight_method_in=CLM.critical_daylight_constant)
    end

    # =====================================================================
    # cn_phenology_init!
    # =====================================================================
    @testset "cn_phenology_init!" begin
        ps = CLM.PhenologyState()
        pp = CLM.PhenologyParams()
        dt = 1800.0  # 30-min timestep

        CLM.cn_phenology_init!(ps, pp, dt)

        @test ps.dt ≈ 1800.0
        @test ps.fracday ≈ 1800.0 / CLM.SECSPDAY
        @test ps.crit_dayl ≈ 39200.0
        @test ps.ndays_off ≈ 30.0
        @test ps.fstor2tran ≈ 0.5
        @test ps.phenology_soil_layer == 1
        @test ps.crit_onset_fdd ≈ 15.0
        @test ps.soilpsi_on ≈ -0.6
        @test ps.crit_offset_fdd ≈ 15.0
        @test ps.soilpsi_off ≈ -0.8
        # live wood turnover: annual to per-second
        @test ps.lwtop ≈ 0.7 / 31536000.0
    end

    @testset "cn_phenology_init! error checks" begin
        ps = CLM.PhenologyState()
        pp = CLM.PhenologyParams(crit_dayl=39200.0, crit_dayl_at_high_lat=30000.0)

        # Set a method that requires validation
        CLM.cn_phenology_set_nml!(critical_daylight_method_in=CLM.critical_daylight_depends_on_lat)
        @test_throws ErrorException CLM.cn_phenology_init!(ps, pp, 1800.0)

        # Reset
        CLM.cn_phenology_set_nml!(critical_daylight_method_in=CLM.critical_daylight_constant)
    end

    # =====================================================================
    # Helper: make_test_data — sets up minimal data structures for 1 patch
    # =====================================================================
    function make_test_data(;np::Int=1, nc::Int=1, ng::Int=1, nlevdecomp::Int=1, nlitr::Int=3)
        pstate = CLM.PhenologyState()
        params = CLM.PhenologyParams()
        CLM.cn_phenology_init!(pstate, params, 1800.0)

        # PFT constants (indexed by PFT type)
        npft = 20
        pftcon = CLM.PftConPhenology(
            evergreen             = zeros(npft),
            season_decid          = zeros(npft),
            season_decid_temperate= zeros(npft),
            stress_decid          = zeros(npft),
            woody                 = zeros(npft),
            leaf_long             = fill(1.0, npft),
            leafcn                = fill(25.0, npft),
            frootcn               = fill(42.0, npft),
            lflitcn               = fill(50.0, npft),
            livewdcn              = fill(50.0, npft),
            deadwdcn              = fill(500.0, npft),
            ndays_on              = fill(30.0, npft),
            crit_onset_gdd_sf     = fill(1.0, npft),
            lf_f                  = fill(1.0/nlitr, npft, nlitr),
            fr_f                  = fill(1.0/nlitr, npft, nlitr),
            biofuel_harvfrac      = zeros(npft),
            repr_structure_harvfrac = zeros(npft, 1),
            minplanttemp          = fill(280.0, npft),
            planttemp             = fill(283.0, npft),
            gddmin                = fill(50.0, npft),
            lfemerg               = fill(0.04, npft),
            grnfill               = fill(0.65, npft),
            hybgdd                = fill(1700.0, npft),
            mxmat                 = fill(200, npft),
            manunitro             = zeros(npft),
            is_pft_known_to_model = fill(true, npft),
            mnNHplantdate         = fill(100.0, npft),
            mxNHplantdate         = fill(170.0, npft),
            mnSHplantdate         = fill(280.0, npft),
            mxSHplantdate         = fill(350.0, npft),
        )

        patch_data = CLM.PatchData()
        patch_data.itype    = fill(5, np)   # PFT type 5
        patch_data.column   = fill(1, np)
        patch_data.gridcell = fill(1, np)
        patch_data.wtcol    = fill(1.0, np)

        gridcell = CLM.GridcellData()
        gridcell.latdeg   = fill(45.0, ng)
        gridcell.dayl     = fill(45000.0, ng)
        gridcell.prev_dayl = fill(44000.0, ng)

        temperature = CLM.TemperatureData()
        temperature.t_ref2m_patch     = fill(285.0, np)
        temperature.t_ref2m_min_patch = fill(280.0, np)
        temperature.t_ref2m_max_patch = fill(290.0, np)
        temperature.t_soisno_col      = fill(285.0, nc, 10)
        temperature.soila10_col       = fill(280.0, nc)
        temperature.t_a5min_patch     = fill(278.0, np)

        water_diag = CLM.WaterDiagnosticBulkData()
        water_diag.snow_5day_col  = fill(0.0, nc)
        water_diag.snow_depth_col = fill(0.0, nc)

        soil_state = CLM.SoilStateData()
        soil_state.soilpsi_col = fill(-0.5, nc, 10)

        canopy_state = CLM.CanopyStateData()
        canopy_state.tlai_patch = fill(2.0, np)

        cnveg_state = CLM.CNVegStateData()
        cnveg_state.tempavg_t2m_patch      = fill(280.0, np)
        cnveg_state.annavg_t2m_patch       = fill(280.0, np)
        cnveg_state.onset_flag_patch       = fill(0.0, np)
        cnveg_state.onset_counter_patch    = fill(0.0, np)
        cnveg_state.onset_gddflag_patch    = fill(0.0, np)
        cnveg_state.onset_gdd_patch        = fill(0.0, np)
        cnveg_state.onset_fdd_patch        = fill(0.0, np)
        cnveg_state.onset_swi_patch        = fill(0.0, np)
        cnveg_state.offset_flag_patch      = fill(0.0, np)
        cnveg_state.offset_counter_patch   = fill(0.0, np)
        cnveg_state.offset_fdd_patch       = fill(0.0, np)
        cnveg_state.offset_swi_patch       = fill(0.0, np)
        cnveg_state.dormant_flag_patch     = fill(1.0, np)
        cnveg_state.days_active_patch      = fill(0.0, np)
        cnveg_state.bglfr_patch            = fill(0.0, np)
        cnveg_state.bgtr_patch             = fill(0.0, np)
        cnveg_state.lgsf_patch             = fill(0.0, np)
        cnveg_state.huigrain_patch         = fill(0.8, np)
        cnveg_state.huileaf_patch          = fill(0.1, np)
        cnveg_state.idop_patch             = fill(0, np)
        cnveg_state.iyop_patch             = fill(0, np)
        cnveg_state.aleaf_patch            = fill(0.0, np)
        cnveg_state.aleafi_patch           = fill(0.0, np)
        cnveg_state.astem_patch            = fill(0.0, np)
        cnveg_state.astemi_patch           = fill(0.0, np)
        cnveg_state.aroot_patch            = fill(0.0, np)
        cnveg_state.cumvd_patch            = fill(0.0, np)
        cnveg_state.hdidx_patch            = fill(0.0, np)

        cnveg_cs = CLM.CNVegCarbonStateData()
        cnveg_cs.leafc_storage_patch        = fill(10.0, np)
        cnveg_cs.frootc_storage_patch       = fill(5.0, np)
        cnveg_cs.livestemc_storage_patch    = fill(3.0, np)
        cnveg_cs.deadstemc_storage_patch    = fill(2.0, np)
        cnveg_cs.livecrootc_storage_patch   = fill(1.5, np)
        cnveg_cs.deadcrootc_storage_patch   = fill(1.0, np)
        cnveg_cs.gresp_storage_patch        = fill(0.5, np)
        cnveg_cs.leafc_xfer_patch           = fill(0.0, np)
        cnveg_cs.frootc_xfer_patch          = fill(0.0, np)
        cnveg_cs.livestemc_xfer_patch       = fill(0.0, np)
        cnveg_cs.deadstemc_xfer_patch       = fill(0.0, np)
        cnveg_cs.livecrootc_xfer_patch      = fill(0.0, np)
        cnveg_cs.deadcrootc_xfer_patch      = fill(0.0, np)
        cnveg_cs.leafc_patch                = fill(50.0, np)
        cnveg_cs.frootc_patch               = fill(20.0, np)
        cnveg_cs.livestemc_patch            = fill(100.0, np)
        cnveg_cs.livecrootc_patch           = fill(40.0, np)

        cnveg_cf = CLM.CNVegCarbonFluxData()
        cnveg_cf.leafc_storage_to_xfer_patch        = fill(0.0, np)
        cnveg_cf.frootc_storage_to_xfer_patch       = fill(0.0, np)
        cnveg_cf.livestemc_storage_to_xfer_patch    = fill(0.0, np)
        cnveg_cf.deadstemc_storage_to_xfer_patch    = fill(0.0, np)
        cnveg_cf.livecrootc_storage_to_xfer_patch   = fill(0.0, np)
        cnveg_cf.deadcrootc_storage_to_xfer_patch   = fill(0.0, np)
        cnveg_cf.gresp_storage_to_xfer_patch        = fill(0.0, np)
        cnveg_cf.leafc_xfer_to_leafc_patch          = fill(0.0, np)
        cnveg_cf.frootc_xfer_to_frootc_patch        = fill(0.0, np)
        cnveg_cf.livestemc_xfer_to_livestemc_patch  = fill(0.0, np)
        cnveg_cf.deadstemc_xfer_to_deadstemc_patch  = fill(0.0, np)
        cnveg_cf.livecrootc_xfer_to_livecrootc_patch= fill(0.0, np)
        cnveg_cf.deadcrootc_xfer_to_deadcrootc_patch= fill(0.0, np)
        cnveg_cf.leafc_to_litter_patch              = fill(0.0, np)
        cnveg_cf.frootc_to_litter_patch             = fill(0.0, np)
        cnveg_cf.prev_leafc_to_litter_patch         = fill(0.0, np)
        cnveg_cf.prev_frootc_to_litter_patch        = fill(0.0, np)
        cnveg_cf.crop_seedc_to_leaf_patch           = fill(0.0, np)
        cnveg_cf.crop_harvestc_to_cropprodc_patch   = fill(0.0, np)
        cnveg_cf.leafc_to_biofuelc_patch            = fill(0.0, np)
        cnveg_cf.livestemc_to_biofuelc_patch        = fill(0.0, np)
        cnveg_cf.leafc_to_removedresiduec_patch     = fill(0.0, np)
        cnveg_cf.livestemc_to_removedresiduec_patch = fill(0.0, np)
        cnveg_cf.livestemc_to_deadstemc_patch       = fill(0.0, np)
        cnveg_cf.livecrootc_to_deadcrootc_patch     = fill(0.0, np)
        cnveg_cf.livestemc_to_litter_patch          = fill(0.0, np)
        cnveg_cf.phenology_c_to_litr_c_col          = zeros(nc, nlevdecomp, nlitr)

        cnveg_ns = CLM.CNVegNitrogenStateData()
        cnveg_ns.leafn_storage_patch        = fill(0.4, np)
        cnveg_ns.frootn_storage_patch       = fill(0.12, np)
        cnveg_ns.livestemn_storage_patch    = fill(0.06, np)
        cnveg_ns.deadstemn_storage_patch    = fill(0.004, np)
        cnveg_ns.livecrootn_storage_patch   = fill(0.03, np)
        cnveg_ns.deadcrootn_storage_patch   = fill(0.002, np)
        cnveg_ns.leafn_xfer_patch           = fill(0.0, np)
        cnveg_ns.frootn_xfer_patch          = fill(0.0, np)
        cnveg_ns.livestemn_xfer_patch       = fill(0.0, np)
        cnveg_ns.deadstemn_xfer_patch       = fill(0.0, np)
        cnveg_ns.livecrootn_xfer_patch      = fill(0.0, np)
        cnveg_ns.deadcrootn_xfer_patch      = fill(0.0, np)
        cnveg_ns.leafn_patch                = fill(2.0, np)
        cnveg_ns.frootn_patch               = fill(0.5, np)
        cnveg_ns.livestemn_patch            = fill(2.0, np)
        cnveg_ns.livecrootn_patch           = fill(0.8, np)

        cnveg_nf = CLM.CNVegNitrogenFluxData()
        cnveg_nf.leafn_storage_to_xfer_patch        = fill(0.0, np)
        cnveg_nf.frootn_storage_to_xfer_patch       = fill(0.0, np)
        cnveg_nf.livestemn_storage_to_xfer_patch    = fill(0.0, np)
        cnveg_nf.deadstemn_storage_to_xfer_patch    = fill(0.0, np)
        cnveg_nf.livecrootn_storage_to_xfer_patch   = fill(0.0, np)
        cnveg_nf.deadcrootn_storage_to_xfer_patch   = fill(0.0, np)
        cnveg_nf.leafn_xfer_to_leafn_patch          = fill(0.0, np)
        cnveg_nf.frootn_xfer_to_frootn_patch        = fill(0.0, np)
        cnveg_nf.livestemn_xfer_to_livestemn_patch  = fill(0.0, np)
        cnveg_nf.deadstemn_xfer_to_deadstemn_patch  = fill(0.0, np)
        cnveg_nf.livecrootn_xfer_to_livecrootn_patch= fill(0.0, np)
        cnveg_nf.deadcrootn_xfer_to_deadcrootn_patch= fill(0.0, np)
        cnveg_nf.leafn_to_litter_patch              = fill(0.0, np)
        cnveg_nf.frootn_to_litter_patch             = fill(0.0, np)
        cnveg_nf.leafn_to_retransn_patch            = fill(0.0, np)
        cnveg_nf.crop_seedn_to_leaf_patch           = fill(0.0, np)
        cnveg_nf.crop_harvestn_to_cropprodn_patch   = fill(0.0, np)
        cnveg_nf.leafn_to_biofueln_patch            = fill(0.0, np)
        cnveg_nf.livestemn_to_biofueln_patch        = fill(0.0, np)
        cnveg_nf.leafn_to_removedresiduen_patch     = fill(0.0, np)
        cnveg_nf.livestemn_to_removedresiduen_patch = fill(0.0, np)
        cnveg_nf.livestemn_to_deadstemn_patch       = fill(0.0, np)
        cnveg_nf.livecrootn_to_deadcrootn_patch     = fill(0.0, np)
        cnveg_nf.livestemn_to_retransn_patch        = fill(0.0, np)
        cnveg_nf.livecrootn_to_retransn_patch       = fill(0.0, np)
        cnveg_nf.livestemn_to_litter_patch          = fill(0.0, np)
        cnveg_nf.phenology_n_to_litr_n_col          = zeros(nc, nlevdecomp, nlitr)

        crop = CLM.CropData()
        CLM.crop_init!(crop, np)

        cn_params = CLM.CNSharedParamsData()

        mask_soilp  = BitVector(fill(true, np))
        mask_pcropp = BitVector(fill(false, np))
        mask_soilc  = BitVector(fill(true, nc))

        return (pstate=pstate, params=params, pftcon=pftcon,
                patch_data=patch_data, gridcell=gridcell,
                temperature=temperature, water_diag=water_diag,
                soil_state=soil_state, canopy_state=canopy_state,
                cnveg_state=cnveg_state, cnveg_cs=cnveg_cs,
                cnveg_cf=cnveg_cf, cnveg_ns=cnveg_ns, cnveg_nf=cnveg_nf,
                crop=crop, cn_params=cn_params,
                mask_soilp=mask_soilp, mask_pcropp=mask_pcropp,
                mask_soilc=mask_soilc)
    end

    # =====================================================================
    # seasonal_decid_onset
    # =====================================================================
    @testset "seasonal_decid_onset" begin
        fracday = 1800.0 / CLM.SECSPDAY

        # Reset NML state
        CLM.cn_phenology_set_nml!(onset_thresh_depends_on_veg=false,
                                   critical_daylight_method_in=CLM.critical_daylight_constant)

        # Case 1: GDD flag not set, winter solstice → should start accumulation
        result = CLM.seasonal_decid_onset(
            0.0,   # onset_gdd
            0.0,   # onset_gddflag
            CLM.TFRZ + 10.0,  # soilt (warm)
            280.0,  # soila10
            278.0,  # t_a5min
            50000.0, # dayl (long)
            0.0,    # snow_5day
            1.0,    # ws_flag (increasing daylight)
            100.0,  # crit_onset_gdd
            1.0,    # season_decid_temperate
            fracday,
            0.2     # snow5d_thresh
        )
        @test result.onset_gddflag == 1.0
        @test result.onset_gdd > 0.0
        @test result.result == false  # not enough GDD yet

        # Case 2: GDD already above threshold → onset should trigger
        result = CLM.seasonal_decid_onset(
            200.0,  # onset_gdd (above crit)
            1.0,    # onset_gddflag
            CLM.TFRZ + 10.0,
            280.0, 278.0, 50000.0, 0.0,
            1.0,    # ws_flag
            100.0,  # crit_onset_gdd
            1.0, fracday, 0.2
        )
        @test result.result == true

        # Case 3: GDD below threshold, still accumulating
        result = CLM.seasonal_decid_onset(
            50.0,   # onset_gdd (below 100)
            1.0,    # onset_gddflag
            CLM.TFRZ + 5.0,  # soilt
            280.0, 278.0, 50000.0, 0.0,
            1.0, 100.0, 1.0, fracday, 0.2
        )
        @test result.result == false
        @test result.onset_gdd > 50.0  # accumulated some

        # Case 4: Past summer solstice, GDD flag reset
        result = CLM.seasonal_decid_onset(
            50.0, 1.0,
            CLM.TFRZ + 10.0, 280.0, 278.0, 50000.0, 0.0,
            0.0,    # ws_flag=0 (decreasing daylight)
            100.0, 1.0, fracday, 0.2
        )
        @test result.onset_gddflag == 0.0
        @test result.onset_gdd == 0.0
        @test result.result == false
    end

    # =====================================================================
    # seasonal_critical_daylength
    # =====================================================================
    @testset "seasonal_critical_daylength" begin
        d = make_test_data()

        # Method = constant
        CLM.cn_phenology_set_nml!(critical_daylight_method_in=CLM.critical_daylight_constant)
        cd = CLM.seasonal_critical_daylength(1, 1, d.pstate.crit_dayl, d.params, d.pftcon,
                                              d.gridcell, d.patch_data)
        @test cd ≈ d.pstate.crit_dayl

        # Method = depends_on_veg, temperate → should return crit_dayl
        CLM.cn_phenology_set_nml!(critical_daylight_method_in=CLM.critical_daylight_depends_on_veg)
        d.pftcon.season_decid_temperate[5] = 1.0
        cd = CLM.seasonal_critical_daylength(1, 1, d.pstate.crit_dayl, d.params, d.pftcon,
                                              d.gridcell, d.patch_data)
        @test cd ≈ d.pstate.crit_dayl

        # Method = depends_on_veg, not temperate → should return crit_dayl_at_high_lat
        d.pftcon.season_decid_temperate[5] = 0.0
        cd = CLM.seasonal_critical_daylength(1, 1, d.pstate.crit_dayl, d.params, d.pftcon,
                                              d.gridcell, d.patch_data)
        @test cd ≈ d.params.crit_dayl_at_high_lat

        # Method = depends_on_lat
        CLM.cn_phenology_set_nml!(critical_daylight_method_in=CLM.critical_daylight_depends_on_lat)
        cd = CLM.seasonal_critical_daylength(1, 1, d.pstate.crit_dayl, d.params, d.pftcon,
                                              d.gridcell, d.patch_data)
        expected = d.params.crit_dayl_at_high_lat -
            d.params.crit_dayl_lat_slope * (CLM.critical_offset_high_lat - abs(d.gridcell.latdeg[1]))
        expected = max(expected, d.pstate.crit_dayl)
        @test cd ≈ expected

        # Reset
        CLM.cn_phenology_set_nml!(critical_daylight_method_in=CLM.critical_daylight_constant)
    end

    # =====================================================================
    # _is_doy_in_interval
    # =====================================================================
    @testset "_is_doy_in_interval" begin
        # Normal interval (non-wrapping)
        @test CLM._is_doy_in_interval(100, 200, 150) == true
        @test CLM._is_doy_in_interval(100, 200, 50) == false
        @test CLM._is_doy_in_interval(100, 200, 250) == false
        @test CLM._is_doy_in_interval(100, 200, 100) == true
        @test CLM._is_doy_in_interval(100, 200, 200) == true

        # Wrapping interval (start > end, e.g. Nov to Feb)
        @test CLM._is_doy_in_interval(300, 60, 350) == true
        @test CLM._is_doy_in_interval(300, 60, 30) == true
        @test CLM._is_doy_in_interval(300, 60, 150) == false
        @test CLM._is_doy_in_interval(300, 60, 300) == true
        @test CLM._is_doy_in_interval(300, 60, 60) == true
    end

    # =====================================================================
    # get_swindow
    # =====================================================================
    @testset "get_swindow" begin
        # No prescribed windows (all negative) → use params
        rx_starts = [-1, -1]
        rx_ends   = [-1, -1]
        result = CLM.get_swindow(100, rx_starts, rx_ends, 90, 180)
        @test result.w == 1
        @test result.start_w == 90
        @test result.end_w == 180

        # Prescribed windows: find the right one
        rx_starts = [80, 200]
        rx_ends   = [120, 260]

        result = CLM.get_swindow(100, rx_starts, rx_ends, 80, 120)
        @test result.w == 1
        @test result.start_w == 80
        @test result.end_w == 120

        result = CLM.get_swindow(210, rx_starts, rx_ends, 80, 120)
        @test result.w == 2
        @test result.start_w == 200
        @test result.end_w == 260

        # Past last window → wraps to first
        result = CLM.get_swindow(300, rx_starts, rx_ends, 80, 120)
        @test result.w == 1
        @test result.start_w == 80
        @test result.end_w == 120
    end

    # =====================================================================
    # was_sown_in_this_window
    # =====================================================================
    @testset "was_sown_in_this_window" begin
        # Non-wrapping window, sown inside, before today → true
        @test CLM.was_sown_in_this_window(100, 200, 150, 110, true) == true

        # Non-wrapping window, sown after today → false (same window, future)
        @test CLM.was_sown_in_this_window(100, 200, 120, 150, true) == false

        # Not in window → false
        @test CLM.was_sown_in_this_window(100, 200, 50, 110, true) == false
        @test CLM.was_sown_in_this_window(100, 200, 250, 110, true) == false
    end

    # =====================================================================
    # days_past_planting
    # =====================================================================
    @testset "days_past_planting" begin
        @test CLM.days_past_planting(100, 150) == 50
        @test CLM.days_past_planting(100, 100) == 0
        # Year wraparound
        @test CLM.days_past_planting(350, 10) == 25
        @test CLM.days_past_planting(350, 10; dayspyr=365) == 25
    end

    # =====================================================================
    # cn_evergreen_phenology!
    # =====================================================================
    @testset "cn_evergreen_phenology!" begin
        d = make_test_data()
        # Set PFT 5 as evergreen woody
        d.pftcon.evergreen[5] = 1.0
        d.pftcon.woody[5] = 1.0
        d.pftcon.leaf_long[5] = 3.0

        CLM.cn_evergreen_phenology!(d.pstate, d.mask_soilp, d.pftcon,
            d.cnveg_state, d.cnveg_cs, d.cnveg_ns, d.cnveg_cf, d.cnveg_nf,
            d.patch_data, 365.0)

        # Background litterfall rate
        expected_bglfr = 1.0 / (3.0 * 365.0 * CLM.SECSPDAY)
        @test d.cnveg_state.bglfr_patch[1] ≈ expected_bglfr

        # Storage → transfer fluxes (fstor2tran=0.5, dt=1800)
        fstor2tran = d.pstate.fstor2tran
        dt = d.pstate.dt

        @test d.cnveg_cf.leafc_storage_to_xfer_patch[1] ≈ fstor2tran * 10.0 / dt
        @test d.cnveg_cf.frootc_storage_to_xfer_patch[1] ≈ fstor2tran * 5.0 / dt

        # Woody: stem and coarse root storage transfers
        @test d.cnveg_cf.livestemc_storage_to_xfer_patch[1] ≈ fstor2tran * 3.0 / dt
        @test d.cnveg_cf.deadstemc_storage_to_xfer_patch[1] ≈ fstor2tran * 2.0 / dt
        @test d.cnveg_cf.livecrootc_storage_to_xfer_patch[1] ≈ fstor2tran * 1.5 / dt
        @test d.cnveg_cf.deadcrootc_storage_to_xfer_patch[1] ≈ fstor2tran * 1.0 / dt
        @test d.cnveg_cf.gresp_storage_to_xfer_patch[1] ≈ fstor2tran * 0.5 / dt

        # Nitrogen fluxes
        @test d.cnveg_nf.leafn_storage_to_xfer_patch[1] ≈ fstor2tran * 0.4 / dt
        @test d.cnveg_nf.frootn_storage_to_xfer_patch[1] ≈ fstor2tran * 0.12 / dt
    end

    @testset "cn_evergreen_phenology! non-woody" begin
        d = make_test_data()
        d.pftcon.evergreen[5] = 1.0
        d.pftcon.woody[5] = 0.0  # not woody

        CLM.cn_evergreen_phenology!(d.pstate, d.mask_soilp, d.pftcon,
            d.cnveg_state, d.cnveg_cs, d.cnveg_ns, d.cnveg_cf, d.cnveg_nf,
            d.patch_data, 365.0)

        # Leaf and froot transfers happen
        @test d.cnveg_cf.leafc_storage_to_xfer_patch[1] > 0.0
        @test d.cnveg_cf.frootc_storage_to_xfer_patch[1] > 0.0

        # Woody transfers should NOT happen
        @test d.cnveg_cf.livestemc_storage_to_xfer_patch[1] == 0.0
        @test d.cnveg_cf.deadstemc_storage_to_xfer_patch[1] == 0.0
    end

    @testset "cn_evergreen_phenology! skips non-evergreen" begin
        d = make_test_data()
        d.pftcon.evergreen[5] = 0.0  # not evergreen

        CLM.cn_evergreen_phenology!(d.pstate, d.mask_soilp, d.pftcon,
            d.cnveg_state, d.cnveg_cs, d.cnveg_ns, d.cnveg_cf, d.cnveg_nf,
            d.patch_data, 365.0)

        @test d.cnveg_cf.leafc_storage_to_xfer_patch[1] == 0.0
    end

    # =====================================================================
    # cn_onset_growth!
    # =====================================================================
    @testset "cn_onset_growth! during onset" begin
        d = make_test_data()
        d.cnveg_state.onset_flag_patch[1] = 1.0
        d.cnveg_state.onset_counter_patch[1] = 30.0 * CLM.SECSPDAY

        d.cnveg_cs.leafc_xfer_patch[1]  = 5.0
        d.cnveg_cs.frootc_xfer_patch[1] = 2.5

        CLM.cn_onset_growth!(d.pstate, d.mask_soilp, d.pftcon,
            d.cnveg_state, d.cnveg_cs, d.cnveg_ns, d.cnveg_cf, d.cnveg_nf,
            d.patch_data)

        # t1 = 2 / onset_counter
        t1 = 2.0 / (30.0 * CLM.SECSPDAY)
        @test d.cnveg_cf.leafc_xfer_to_leafc_patch[1] ≈ t1 * 5.0
        @test d.cnveg_cf.frootc_xfer_to_frootc_patch[1] ≈ t1 * 2.5
    end

    @testset "cn_onset_growth! last timestep" begin
        d = make_test_data()
        dt = d.pstate.dt
        d.cnveg_state.onset_flag_patch[1] = 1.0
        d.cnveg_state.onset_counter_patch[1] = dt  # exactly one timestep left

        d.cnveg_cs.leafc_xfer_patch[1] = 5.0

        CLM.cn_onset_growth!(d.pstate, d.mask_soilp, d.pftcon,
            d.cnveg_state, d.cnveg_cs, d.cnveg_ns, d.cnveg_cf, d.cnveg_nf,
            d.patch_data)

        # t1 = 1/dt on last timestep
        @test d.cnveg_cf.leafc_xfer_to_leafc_patch[1] ≈ 5.0 / dt
    end

    @testset "cn_onset_growth! background transfer" begin
        d = make_test_data()
        dt = d.pstate.dt
        d.cnveg_state.onset_flag_patch[1] = 0.0  # not in onset
        d.cnveg_state.bgtr_patch[1] = 1.0e-7     # positive background transfer

        d.cnveg_cs.leafc_xfer_patch[1] = 3.0
        d.cnveg_cs.frootc_xfer_patch[1] = 1.5

        CLM.cn_onset_growth!(d.pstate, d.mask_soilp, d.pftcon,
            d.cnveg_state, d.cnveg_cs, d.cnveg_ns, d.cnveg_cf, d.cnveg_nf,
            d.patch_data)

        @test d.cnveg_cf.leafc_xfer_to_leafc_patch[1] ≈ 3.0 / dt
        @test d.cnveg_cf.frootc_xfer_to_frootc_patch[1] ≈ 1.5 / dt
    end

    # =====================================================================
    # cn_offset_litterfall!
    # =====================================================================
    @testset "cn_offset_litterfall! during offset" begin
        d = make_test_data()
        dt = d.pstate.dt
        d.cnveg_state.offset_flag_patch[1] = 1.0
        d.cnveg_state.offset_counter_patch[1] = 20.0 * CLM.SECSPDAY

        CLM.cn_offset_litterfall!(d.pstate, d.mask_soilp, d.pftcon,
            d.cnveg_state, d.cnveg_cs, d.cnveg_ns, d.cnveg_cf, d.cnveg_nf,
            d.crop, d.patch_data)

        # Should produce positive litter fluxes
        @test d.cnveg_cf.leafc_to_litter_patch[1] > 0.0
        @test d.cnveg_cf.frootc_to_litter_patch[1] > 0.0

        # Nitrogen should also flow to litter
        @test d.cnveg_nf.leafn_to_litter_patch[1] > 0.0
        @test d.cnveg_nf.frootn_to_litter_patch[1] > 0.0
    end

    @testset "cn_offset_litterfall! last timestep" begin
        d = make_test_data()
        dt = d.pstate.dt
        d.cnveg_state.offset_flag_patch[1] = 1.0
        d.cnveg_state.offset_counter_patch[1] = dt  # last step

        CLM.cn_offset_litterfall!(d.pstate, d.mask_soilp, d.pftcon,
            d.cnveg_state, d.cnveg_cs, d.cnveg_ns, d.cnveg_cf, d.cnveg_nf,
            d.crop, d.patch_data)

        # Dump all remaining: leafc_to_litter = (1/dt) * leafc
        @test d.cnveg_cf.leafc_to_litter_patch[1] ≈ d.cnveg_cs.leafc_patch[1] / dt
        @test d.cnveg_cf.frootc_to_litter_patch[1] ≈ d.cnveg_cs.frootc_patch[1] / dt
    end

    @testset "cn_offset_litterfall! CNratio_floating" begin
        d = make_test_data()
        dt = d.pstate.dt
        d.cnveg_state.offset_flag_patch[1] = 1.0
        d.cnveg_state.offset_counter_patch[1] = dt

        CLM.cn_offset_litterfall!(d.pstate, d.mask_soilp, d.pftcon,
            d.cnveg_state, d.cnveg_cs, d.cnveg_ns, d.cnveg_cf, d.cnveg_nf,
            d.crop, d.patch_data; CNratio_floating=true)

        # With floating CN: N litter = 0.5 * leafc_to_litter * (leafn/leafc)
        leafc_to_lit = d.cnveg_cs.leafc_patch[1] / dt
        ntovr = leafc_to_lit * (d.cnveg_ns.leafn_patch[1] / d.cnveg_cs.leafc_patch[1])
        expected_n_litter = 0.5 * ntovr
        @test d.cnveg_nf.leafn_to_litter_patch[1] ≈ expected_n_litter
        @test d.cnveg_nf.leafn_to_retransn_patch[1] ≈ ntovr - expected_n_litter
    end

    @testset "cn_offset_litterfall! skips non-offset" begin
        d = make_test_data()
        d.cnveg_state.offset_flag_patch[1] = 0.0

        CLM.cn_offset_litterfall!(d.pstate, d.mask_soilp, d.pftcon,
            d.cnveg_state, d.cnveg_cs, d.cnveg_ns, d.cnveg_cf, d.cnveg_nf,
            d.crop, d.patch_data)

        @test d.cnveg_cf.leafc_to_litter_patch[1] == 0.0
    end

    # =====================================================================
    # cn_background_litterfall!
    # =====================================================================
    @testset "cn_background_litterfall!" begin
        d = make_test_data()
        bglfr = 1.0e-7
        d.cnveg_state.bglfr_patch[1] = bglfr

        CLM.cn_background_litterfall!(d.pstate, d.mask_soilp, d.pftcon,
            d.cnveg_state, d.cnveg_cs, d.cnveg_ns, d.cnveg_cf, d.cnveg_nf,
            d.patch_data)

        @test d.cnveg_cf.leafc_to_litter_patch[1] ≈ bglfr * d.cnveg_cs.leafc_patch[1]
        @test d.cnveg_cf.frootc_to_litter_patch[1] ≈ bglfr * d.cnveg_cs.frootc_patch[1]

        # Fixed CN: N litter = C litter / lflitcn
        @test d.cnveg_nf.leafn_to_litter_patch[1] ≈ d.cnveg_cf.leafc_to_litter_patch[1] / d.pftcon.lflitcn[5]
    end

    @testset "cn_background_litterfall! zero bglfr" begin
        d = make_test_data()
        d.cnveg_state.bglfr_patch[1] = 0.0

        CLM.cn_background_litterfall!(d.pstate, d.mask_soilp, d.pftcon,
            d.cnveg_state, d.cnveg_cs, d.cnveg_ns, d.cnveg_cf, d.cnveg_nf,
            d.patch_data)

        @test d.cnveg_cf.leafc_to_litter_patch[1] == 0.0
    end

    # =====================================================================
    # cn_livewood_turnover!
    # =====================================================================
    @testset "cn_livewood_turnover! woody" begin
        d = make_test_data()
        d.pftcon.woody[5] = 1.0

        CLM.cn_livewood_turnover!(d.pstate, d.mask_soilp, d.pftcon,
            d.cnveg_cs, d.cnveg_ns, d.cnveg_cf, d.cnveg_nf,
            d.patch_data)

        lwtop = d.pstate.lwtop
        # live stem → dead stem
        expected_stem_c = d.cnveg_cs.livestemc_patch[1] * lwtop
        @test d.cnveg_cf.livestemc_to_deadstemc_patch[1] ≈ expected_stem_c

        expected_stem_n_dead = expected_stem_c / d.pftcon.deadwdcn[5]
        @test d.cnveg_nf.livestemn_to_deadstemn_patch[1] ≈ expected_stem_n_dead

        # Retranslocation
        ntovr = expected_stem_c / d.pftcon.livewdcn[5]
        @test d.cnveg_nf.livestemn_to_retransn_patch[1] ≈ ntovr - expected_stem_n_dead

        # live coarse root → dead coarse root
        expected_croot_c = d.cnveg_cs.livecrootc_patch[1] * lwtop
        @test d.cnveg_cf.livecrootc_to_deadcrootc_patch[1] ≈ expected_croot_c
    end

    @testset "cn_livewood_turnover! non-woody skipped" begin
        d = make_test_data()
        d.pftcon.woody[5] = 0.0

        CLM.cn_livewood_turnover!(d.pstate, d.mask_soilp, d.pftcon,
            d.cnveg_cs, d.cnveg_ns, d.cnveg_cf, d.cnveg_nf,
            d.patch_data)

        @test d.cnveg_cf.livestemc_to_deadstemc_patch[1] == 0.0
        @test d.cnveg_cf.livecrootc_to_deadcrootc_patch[1] == 0.0
    end

    # =====================================================================
    # cn_crop_harvest_to_product_pools!
    # =====================================================================
    @testset "cn_crop_harvest_to_product_pools!" begin
        d = make_test_data()

        # Set biofuel and residue fluxes
        d.cnveg_cf.leafc_to_biofuelc_patch[1] = 1.0
        d.cnveg_cf.livestemc_to_biofuelc_patch[1] = 2.0
        d.cnveg_cf.leafc_to_removedresiduec_patch[1] = 0.5
        d.cnveg_cf.livestemc_to_removedresiduec_patch[1] = 0.3

        d.cnveg_nf.leafn_to_biofueln_patch[1] = 0.04
        d.cnveg_nf.livestemn_to_biofueln_patch[1] = 0.08
        d.cnveg_nf.leafn_to_removedresiduen_patch[1] = 0.02
        d.cnveg_nf.livestemn_to_removedresiduen_patch[1] = 0.01

        CLM.cn_crop_harvest_to_product_pools!(d.mask_soilp, d.mask_soilc,
            d.cnveg_cf, d.cnveg_nf, d.patch_data; use_crop=true)

        @test d.cnveg_cf.crop_harvestc_to_cropprodc_patch[1] ≈ 1.0 + 2.0 + 0.5 + 0.3
        @test d.cnveg_nf.crop_harvestn_to_cropprodn_patch[1] ≈ 0.04 + 0.08 + 0.02 + 0.01
    end

    @testset "cn_crop_harvest_to_product_pools! no-op without crop" begin
        d = make_test_data()
        d.cnveg_cf.leafc_to_biofuelc_patch[1] = 1.0

        CLM.cn_crop_harvest_to_product_pools!(d.mask_soilp, d.mask_soilc,
            d.cnveg_cf, d.cnveg_nf, d.patch_data; use_crop=false)

        @test d.cnveg_cf.crop_harvestc_to_cropprodc_patch[1] == 0.0
    end

    # =====================================================================
    # cn_litter_to_column!
    # =====================================================================
    @testset "cn_litter_to_column!" begin
        d = make_test_data()
        nlevdecomp = 1
        nlitr = 3

        leaf_prof  = ones(1, nlevdecomp)
        froot_prof = ones(1, nlevdecomp)

        # Set leaf and froot litter
        d.cnveg_cf.leafc_to_litter_patch[1]  = 10.0
        d.cnveg_cf.frootc_to_litter_patch[1] = 5.0
        d.cnveg_nf.leafn_to_litter_patch[1]  = 0.4
        d.cnveg_nf.frootn_to_litter_patch[1] = 0.12

        CLM.cn_litter_to_column!(d.mask_soilp, d.pftcon,
            d.cnveg_state, d.cnveg_cf, d.cnveg_nf,
            d.patch_data, leaf_prof, froot_prof;
            nlevdecomp=nlevdecomp, i_litr_min=1, i_litr_max=nlitr)

        # Each litter fraction gets lf_f = 1/3 of total, wtcol=1.0, leaf_prof=1.0
        for i in 1:nlitr
            expected_c = 10.0 * (1.0/nlitr) * 1.0 * 1.0 + 5.0 * (1.0/nlitr) * 1.0 * 1.0
            @test d.cnveg_cf.phenology_c_to_litr_c_col[1, 1, i] ≈ expected_c
        end

        # Total C across litter types
        total_c = sum(d.cnveg_cf.phenology_c_to_litr_c_col[1, 1, :])
        @test total_c ≈ 15.0  # 10 + 5
    end

    # =====================================================================
    # plant_crop!
    # =====================================================================
    @testset "plant_crop!" begin
        d = make_test_data()
        dt = d.pstate.dt
        leafcn = 25.0

        CLM.plant_crop!(1, leafcn, 120, 2025,
            d.crop, d.cnveg_state, d.cnveg_cs, d.cnveg_ns,
            d.cnveg_cf, d.cnveg_nf, dt)

        @test d.crop.croplive_patch[1] == true
        @test d.crop.sown_in_this_window[1] == true
        @test d.cnveg_state.idop_patch[1] == 120
        @test d.cnveg_state.iyop_patch[1] == 2025
        @test d.crop.harvdate_patch[1] == CLM.NOT_Harvested

        seed = CLM._initial_seed_at_planting[]
        @test d.cnveg_cs.leafc_xfer_patch[1] ≈ seed
        @test d.cnveg_ns.leafn_xfer_patch[1] ≈ seed / leafcn
        @test d.cnveg_cf.crop_seedc_to_leaf_patch[1] ≈ seed / dt
        @test d.cnveg_nf.crop_seedn_to_leaf_patch[1] ≈ (seed / leafcn) / dt

        @test d.cnveg_state.aleaf_patch[1] ≈ 1.0
        @test d.cnveg_state.aleafi_patch[1] ≈ 1.0
        @test d.cnveg_state.astem_patch[1] ≈ 0.0
    end

    # =====================================================================
    # vernalization!
    # =====================================================================
    @testset "vernalization! warm conditions" begin
        d = make_test_data()
        # Warm conditions: t_ref2m > TFRZ
        d.temperature.t_ref2m_patch[1]     = CLM.TFRZ + 10.0
        d.temperature.t_ref2m_min_patch[1] = CLM.TFRZ + 5.0
        d.temperature.t_ref2m_max_patch[1] = CLM.TFRZ + 15.0

        force_harvest = CLM.vernalization!(1, d.pstate, d.canopy_state,
            d.temperature, d.water_diag, d.cnveg_state, d.crop, d.patch_data)

        @test force_harvest == false
        # Vernalization days should have accumulated (warm but within range)
        @test d.cnveg_state.cumvd_patch[1] >= 0.0
        @test isfinite(d.crop.vf_patch[1])
    end

    @testset "vernalization! cold conditions" begin
        d = make_test_data()
        # Cold conditions: t_ref2m well below freezing
        d.temperature.t_ref2m_patch[1]     = CLM.TFRZ - 5.0
        d.temperature.t_ref2m_min_patch[1] = CLM.TFRZ - 10.0
        d.temperature.t_ref2m_max_patch[1] = CLM.TFRZ + 5.0
        d.canopy_state.tlai_patch[1] = 2.0

        force_harvest = CLM.vernalization!(1, d.pstate, d.canopy_state,
            d.temperature, d.water_diag, d.cnveg_state, d.crop, d.patch_data)

        @test typeof(force_harvest) == Bool
        @test isfinite(d.cnveg_state.cumvd_patch[1])
    end

    # =====================================================================
    # crop_phase!
    # =====================================================================
    @testset "crop_phase!" begin
        d = make_test_data()
        mask = BitVector([true])
        phase_out = fill(CLM.cphase_not_planted, 1)

        # Not live → no change
        d.crop.croplive_patch[1] = false
        CLM.crop_phase!(mask, d.crop, d.cnveg_state, phase_out)
        @test phase_out[1] == CLM.cphase_not_planted

        # Planted but below leaf emerge threshold
        d.crop.croplive_patch[1] = true
        d.crop.gddtsoi_patch[1] = 0.0
        d.crop.hui_patch[1] = 0.0
        d.cnveg_state.huileaf_patch[1] = 0.1
        d.cnveg_state.huigrain_patch[1] = 0.8
        CLM.crop_phase!(mask, d.crop, d.cnveg_state, phase_out)
        @test phase_out[1] == CLM.cphase_planted

        # Leaf emerged
        d.crop.gddtsoi_patch[1] = 100.0
        d.crop.hui_patch[1] = 0.5
        CLM.crop_phase!(mask, d.crop, d.cnveg_state, phase_out)
        @test phase_out[1] == CLM.cphase_leafemerge

        # Grain fill
        d.crop.hui_patch[1] = 0.9
        CLM.crop_phase!(mask, d.crop, d.cnveg_state, phase_out)
        @test phase_out[1] == CLM.cphase_grainfill
    end

    # =====================================================================
    # cn_phenology_climate!
    # =====================================================================
    @testset "cn_phenology_climate!" begin
        d = make_test_data()
        old_tempavg = d.cnveg_state.tempavg_t2m_patch[1]

        CLM.cn_phenology_climate!(d.pstate, d.mask_soilp, d.temperature,
            d.cnveg_state, d.crop, d.patch_data, d.pftcon)

        # tempavg should move towards t_ref2m
        @test d.cnveg_state.tempavg_t2m_patch[1] != old_tempavg
        # Should move in direction of t_ref2m
        t_ref = d.temperature.t_ref2m_patch[1]
        if t_ref > old_tempavg
            @test d.cnveg_state.tempavg_t2m_patch[1] > old_tempavg
        else
            @test d.cnveg_state.tempavg_t2m_patch[1] < old_tempavg
        end
    end

    # =====================================================================
    # Full cn_phenology! driver — phase error
    # =====================================================================
    @testset "cn_phenology! bad phase" begin
        d = make_test_data()
        varctl = CLM.VarCtl()
        leaf_prof = ones(1, 1)
        froot_prof = ones(1, 1)

        @test_throws ErrorException CLM.cn_phenology!(d.pstate, d.params, d.pftcon,
            d.mask_soilp, d.mask_pcropp, d.mask_soilc,
            d.temperature, d.water_diag, d.canopy_state, d.soil_state,
            d.cnveg_state, d.cnveg_cs, d.cnveg_cf, d.cnveg_ns, d.cnveg_nf,
            d.crop, d.patch_data, d.gridcell, d.cn_params,
            leaf_prof, froot_prof, 3; varctl=varctl)
    end

    # =====================================================================
    # Full cn_phenology! — phase 1 evergreen
    # =====================================================================
    @testset "cn_phenology! phase 1 evergreen" begin
        d = make_test_data()
        d.pftcon.evergreen[5] = 1.0
        d.pftcon.woody[5] = 1.0
        d.pftcon.leaf_long[5] = 3.0

        varctl = CLM.VarCtl()
        leaf_prof = ones(1, 1)
        froot_prof = ones(1, 1)

        CLM.cn_phenology!(d.pstate, d.params, d.pftcon,
            d.mask_soilp, d.mask_pcropp, d.mask_soilc,
            d.temperature, d.water_diag, d.canopy_state, d.soil_state,
            d.cnveg_state, d.cnveg_cs, d.cnveg_cf, d.cnveg_ns, d.cnveg_nf,
            d.crop, d.patch_data, d.gridcell, d.cn_params,
            leaf_prof, froot_prof, 1; varctl=varctl)

        # Evergreen should have set bglfr and storage→transfer
        @test d.cnveg_state.bglfr_patch[1] > 0.0
        @test d.cnveg_cf.leafc_storage_to_xfer_patch[1] > 0.0
    end

    # =====================================================================
    # Full cn_phenology! — phase 2 onset+offset+bg+turnover+litter
    # =====================================================================
    @testset "cn_phenology! phase 2 onset growth" begin
        d = make_test_data()
        d.pftcon.woody[5] = 1.0
        d.cnveg_state.onset_flag_patch[1] = 1.0
        d.cnveg_state.onset_counter_patch[1] = 20.0 * CLM.SECSPDAY
        d.cnveg_cs.leafc_xfer_patch[1] = 5.0
        d.cnveg_ns.leafn_xfer_patch[1] = 0.2

        varctl = CLM.VarCtl()
        leaf_prof = ones(1, 1)
        froot_prof = ones(1, 1)

        CLM.cn_phenology!(d.pstate, d.params, d.pftcon,
            d.mask_soilp, d.mask_pcropp, d.mask_soilc,
            d.temperature, d.water_diag, d.canopy_state, d.soil_state,
            d.cnveg_state, d.cnveg_cs, d.cnveg_cf, d.cnveg_ns, d.cnveg_nf,
            d.crop, d.patch_data, d.gridcell, d.cn_params,
            leaf_prof, froot_prof, 2; varctl=varctl)

        @test d.cnveg_cf.leafc_xfer_to_leafc_patch[1] > 0.0
    end

    # =====================================================================
    # crop_phenology_init!
    # =====================================================================
    @testset "crop_phenology_init!" begin
        ps = CLM.PhenologyState()
        pftcon = CLM.PftConPhenology(
            is_pft_known_to_model = fill(true, 20),
            mnNHplantdate = fill(100.0, 20),
            mxNHplantdate = fill(170.0, 20),
            mnSHplantdate = fill(280.0, 20),
            mxSHplantdate = fill(350.0, 20),
        )
        patch_data = CLM.PatchData()
        patch_data.itype = [5, 17]
        patch_data.gridcell = [1, 1]
        gridcell = CLM.GridcellData()
        gridcell.latdeg = [45.0]

        CLM.crop_phenology_init!(ps, pftcon, patch_data, gridcell, 17, 18, 20)

        @test length(ps.inhemi) == 2
        @test ps.inhemi[1] == CLM.inNH  # lat > 0
        @test ps.inhemi[2] == CLM.inNH
        @test ps.minplantjday[17, CLM.inNH] == 100
        @test ps.maxplantjday[17, CLM.inNH] == 170
        @test ps.minplantjday[17, CLM.inSH] == 280
        @test ps.maxplantjday[17, CLM.inSH] == 350
    end

    @testset "crop_phenology_init! southern hemisphere" begin
        ps = CLM.PhenologyState()
        pftcon = CLM.PftConPhenology(
            is_pft_known_to_model = fill(true, 20),
            mnNHplantdate = fill(100.0, 20),
            mxNHplantdate = fill(170.0, 20),
            mnSHplantdate = fill(280.0, 20),
            mxSHplantdate = fill(350.0, 20),
        )
        patch_data = CLM.PatchData()
        patch_data.itype = [17]
        patch_data.gridcell = [1]
        gridcell = CLM.GridcellData()
        gridcell.latdeg = [-30.0]

        CLM.crop_phenology_init!(ps, pftcon, patch_data, gridcell, 17, 18, 20)

        @test ps.inhemi[1] == CLM.inSH
    end

end
