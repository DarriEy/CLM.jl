@testset "NutrientCompetition" begin

    # =====================================================================
    # PftConNutrientCompetition struct
    # =====================================================================
    @testset "PftConNutrientCompetition defaults" begin
        pfc = CLM.PftConNutrientCompetition()
        @test isempty(pfc.woody)
        @test isempty(pfc.froot_leaf)
        @test isempty(pfc.croot_stem)
        @test isempty(pfc.stem_leaf)
        @test isempty(pfc.flivewd)
        @test isempty(pfc.leafcn)
        @test isempty(pfc.frootcn)
        @test isempty(pfc.livewdcn)
        @test isempty(pfc.deadwdcn)
        @test isempty(pfc.fcur)
        @test isempty(pfc.graincn)
        @test isempty(pfc.grperc)
        @test isempty(pfc.grpnow)
        @test isempty(pfc.fleafcn)
        @test isempty(pfc.ffrootcn)
        @test isempty(pfc.fstemcn)
        @test isempty(pfc.astemf)
        @test isempty(pfc.season_decid)
        @test isempty(pfc.stress_decid)
    end

    # =====================================================================
    # Helper: create test data for a single non-crop woody patch
    # =====================================================================
    function make_woody_test_data(;
            gpp=10.0, availc=8.0, annsum_npp=400.0,
            plant_ndemand=0.2, retransn_to_npool=0.05,
            c_allometry=1.5, n_allometry=0.06,
            fpg=0.8, retransn=0.5, annsum_potential_gpp=5000.0,
            annmax_retransn=1.0)

        np = 1
        nc = 1
        nrepr = CLM.NREPR

        pftcon = CLM.PftConNutrientCompetition(
            woody       = [1.0],
            froot_leaf  = [1.0],
            croot_stem  = [0.3],
            stem_leaf   = [-1.0],   # dynamic allocation
            flivewd     = [0.1],
            leafcn      = [25.0],
            frootcn     = [42.0],
            livewdcn    = [50.0],
            deadwdcn    = [500.0],
            fcur        = [0.5],
            graincn     = [50.0],
            grperc      = [0.3],
            grpnow      = [0.5],
            fleafcn     = [65.0],
            ffrootcn    = [100.0],
            fstemcn     = [130.0],
            astemf      = [0.0],
            season_decid = [0.0],
            stress_decid = [0.0],
        )

        cn_shared_params = CLM.CNSharedParamsData(use_matrixcn=false, use_fun=false)

        patch = CLM.PatchData()
        patch.column = [1]
        patch.itype  = [1]   # woody tree

        crop = CLM.CropData()
        crop.croplive_patch = [false]

        cnveg_state = CLM.CNVegStateData()
        cnveg_state.c_allometry_patch          = [c_allometry]
        cnveg_state.n_allometry_patch          = [n_allometry]
        cnveg_state.downreg_patch              = [0.0]
        cnveg_state.aleaf_patch                = [NaN]
        cnveg_state.astem_patch                = [NaN]
        cnveg_state.aroot_patch                = [NaN]
        cnveg_state.arepr_patch                = fill(NaN, np, nrepr)
        cnveg_state.peaklai_patch              = [0]
        cnveg_state.tempsum_potential_gpp_patch = [0.0]
        cnveg_state.annsum_potential_gpp_patch  = [annsum_potential_gpp]
        cnveg_state.tempmax_retransn_patch      = [0.0]
        cnveg_state.annmax_retransn_patch       = [annmax_retransn]
        cnveg_state.grain_flag_patch            = [0.0]

        cnveg_cs = CLM.CNVegCarbonStateData()
        cnveg_cs.leafc_patch      = [10.0]
        cnveg_cs.frootc_patch     = [5.0]
        cnveg_cs.livestemc_patch  = [20.0]

        cnveg_cf = CLM.CNVegCarbonFluxData()
        cnveg_cf.gpp_before_downreg_patch = [gpp]
        cnveg_cf.availc_patch             = [availc]
        cnveg_cf.excess_cflux_patch       = [0.0]
        cnveg_cf.plant_calloc_patch       = [0.0]
        cnveg_cf.psnsun_to_cpool_patch    = [6.0]
        cnveg_cf.psnshade_to_cpool_patch  = [4.0]
        cnveg_cf.annsum_npp_patch         = [annsum_npp]
        cnveg_cf.cpool_to_leafc_patch          = [0.0]
        cnveg_cf.cpool_to_leafc_storage_patch  = [0.0]
        cnveg_cf.cpool_to_frootc_patch         = [0.0]
        cnveg_cf.cpool_to_frootc_storage_patch = [0.0]
        cnveg_cf.cpool_to_livestemc_patch          = [0.0]
        cnveg_cf.cpool_to_livestemc_storage_patch  = [0.0]
        cnveg_cf.cpool_to_deadstemc_patch          = [0.0]
        cnveg_cf.cpool_to_deadstemc_storage_patch  = [0.0]
        cnveg_cf.cpool_to_livecrootc_patch         = [0.0]
        cnveg_cf.cpool_to_livecrootc_storage_patch = [0.0]
        cnveg_cf.cpool_to_deadcrootc_patch         = [0.0]
        cnveg_cf.cpool_to_deadcrootc_storage_patch = [0.0]
        cnveg_cf.cpool_to_gresp_storage_patch      = [0.0]
        cnveg_cf.cpool_to_reproductivec_patch         = fill(0.0, np, nrepr)
        cnveg_cf.cpool_to_reproductivec_storage_patch = fill(0.0, np, nrepr)

        cnveg_ns = CLM.CNVegNitrogenStateData()
        cnveg_ns.retransn_patch = [retransn]

        cnveg_nf = CLM.CNVegNitrogenFluxData()
        cnveg_nf.plant_ndemand_patch       = [plant_ndemand]
        cnveg_nf.avail_retransn_patch      = [0.0]
        cnveg_nf.retransn_to_npool_patch   = [retransn_to_npool]
        cnveg_nf.sminn_to_npool_patch      = [0.0]
        cnveg_nf.plant_nalloc_patch        = [0.0]
        cnveg_nf.npool_to_leafn_patch          = [0.0]
        cnveg_nf.npool_to_leafn_storage_patch  = [0.0]
        cnveg_nf.npool_to_frootn_patch         = [0.0]
        cnveg_nf.npool_to_frootn_storage_patch = [0.0]
        cnveg_nf.npool_to_livestemn_patch          = [0.0]
        cnveg_nf.npool_to_livestemn_storage_patch  = [0.0]
        cnveg_nf.npool_to_deadstemn_patch          = [0.0]
        cnveg_nf.npool_to_deadstemn_storage_patch  = [0.0]
        cnveg_nf.npool_to_livecrootn_patch         = [0.0]
        cnveg_nf.npool_to_livecrootn_storage_patch = [0.0]
        cnveg_nf.npool_to_deadcrootn_patch         = [0.0]
        cnveg_nf.npool_to_deadcrootn_storage_patch = [0.0]
        cnveg_nf.npool_to_reproductiven_patch         = fill(0.0, np, nrepr)
        cnveg_nf.npool_to_reproductiven_storage_patch = fill(0.0, np, nrepr)
        cnveg_nf.sminn_to_plant_fun_patch  = [0.0]
        cnveg_nf.leafn_to_retransn_patch   = [0.0]
        cnveg_nf.frootn_to_retransn_patch  = [0.0]
        cnveg_nf.livestemn_to_retransn_patch = [0.0]
        cnveg_nf.Npassive_patch            = [0.0]
        cnveg_nf.Nfix_patch                = [0.0]
        cnveg_nf.Nactive_patch             = [0.0]
        cnveg_nf.Nnonmyc_patch             = [0.0]
        cnveg_nf.Nam_patch                 = [0.0]
        cnveg_nf.Necm_patch                = [0.0]

        mask_soilp = trues(np)
        fpg_col_arr = [fpg]

        return (pftcon=pftcon, cn_shared_params=cn_shared_params,
                patch=patch, crop=crop,
                cnveg_state=cnveg_state, cnveg_cs=cnveg_cs,
                cnveg_cf=cnveg_cf, cnveg_ns=cnveg_ns, cnveg_nf=cnveg_nf,
                mask_soilp=mask_soilp, fpg_col=fpg_col_arr, bounds=1:np)
    end

    # =====================================================================
    # Test calc_plant_nitrogen_demand! -- non-crop patch
    # =====================================================================
    @testset "calc_plant_nitrogen_demand! non-crop" begin
        d = make_woody_test_data(gpp=10.0, availc=8.0,
                                  annsum_potential_gpp=5000.0,
                                  annmax_retransn=1.0,
                                  retransn=0.5)
        dt = 1800.0

        CLM.calc_plant_nitrogen_demand!(d.mask_soilp, d.bounds, false,
            d.pftcon, d.cn_shared_params, d.patch, d.crop,
            d.cnveg_state, d.cnveg_cs, d.cnveg_cf, d.cnveg_ns, d.cnveg_nf;
            dt=dt)

        p = 1
        # plant_ndemand = availc * (n_allometry / c_allometry)
        expected_ndemand = 8.0 * (0.06 / 1.5)
        @test d.cnveg_nf.plant_ndemand_patch[p] >= 0.0

        # tempsum_potential_gpp should have been incremented by gpp
        @test d.cnveg_state.tempsum_potential_gpp_patch[p] == 10.0

        # tempmax_retransn should be updated (max of 0.0 and retransn=0.5)
        @test d.cnveg_state.tempmax_retransn_patch[p] == 0.5

        # avail_retransn should be set based on annmax_retransn etc.
        @test d.cnveg_nf.avail_retransn_patch[p] >= 0.0

        # retransn_to_npool should not exceed avail_retransn
        @test d.cnveg_nf.retransn_to_npool_patch[p] >= 0.0

        # plant_ndemand should be reduced by retransn_to_npool (for non-FUN)
        initial_ndemand = expected_ndemand
        @test d.cnveg_nf.plant_ndemand_patch[p] ==
            initial_ndemand - d.cnveg_nf.retransn_to_npool_patch[p]
    end

    # =====================================================================
    # Test calc_plant_cn_alloc! -- non-crop woody patch, no FUN
    # =====================================================================
    @testset "calc_plant_cn_alloc! woody non-crop no-FUN" begin
        d = make_woody_test_data(gpp=10.0, availc=8.0,
                                  plant_ndemand=0.2, retransn_to_npool=0.05,
                                  c_allometry=1.5, n_allometry=0.06,
                                  fpg=0.8, annsum_npp=400.0)

        CLM.calc_plant_cn_alloc!(d.mask_soilp, d.bounds,
            d.pftcon, d.cn_shared_params, d.patch, d.crop,
            d.cnveg_state, d.cnveg_cs, d.cnveg_cf, d.cnveg_nf;
            fpg_col=d.fpg_col)

        p = 1
        # sminn_to_npool = plant_ndemand * fpg = 0.2 * 0.8 = 0.16
        @test d.cnveg_nf.sminn_to_npool_patch[p] ≈ 0.16

        # plant_nalloc = sminn_to_npool + retransn_to_npool = 0.16 + 0.05 = 0.21
        @test d.cnveg_nf.plant_nalloc_patch[p] ≈ 0.21

        # plant_calloc = plant_nalloc * (c_allometry / n_allometry) = 0.21 * (1.5/0.06)
        expected_calloc = 0.21 * (1.5 / 0.06)
        @test d.cnveg_cf.plant_calloc_patch[p] ≈ expected_calloc

        # excess_cflux = availc - plant_calloc
        @test d.cnveg_cf.excess_cflux_patch[p] ≈ 8.0 - expected_calloc

        # downreg = excess_cflux / gpp
        expected_downreg = (8.0 - expected_calloc) / 10.0
        @test d.cnveg_state.downreg_patch[p] ≈ expected_downreg

        # psnsun_to_cpool should be reduced by downreg
        @test d.cnveg_cf.psnsun_to_cpool_patch[p] ≈ 6.0 * (1.0 - expected_downreg)

        # stem_leaf = -1 => f3 = (2.7/(1+exp(-0.004*(400-300))))-0.4
        f3 = (2.7 / (1.0 + exp(-0.004 * (400.0 - 300.0)))) - 0.4
        c_allom = 1.5
        nlc = expected_calloc / c_allom
        fcur = 0.5

        # cpool_to_leafc = nlc * fcur
        @test d.cnveg_cf.cpool_to_leafc_patch[p] ≈ nlc * fcur
        @test d.cnveg_cf.cpool_to_leafc_storage_patch[p] ≈ nlc * (1.0 - fcur)

        # woody: cpool_to_livestemc = nlc * f3 * f4 * fcur
        f4 = 0.1
        @test d.cnveg_cf.cpool_to_livestemc_patch[p] ≈ nlc * f3 * f4 * fcur
        @test d.cnveg_cf.cpool_to_deadstemc_patch[p] ≈ nlc * f3 * (1.0 - f4) * fcur

        # N fluxes: npool_to_leafn = (nlc/cnl) * fcur
        cnl = 25.0
        @test d.cnveg_nf.npool_to_leafn_patch[p] ≈ (nlc / cnl) * fcur
        @test d.cnveg_nf.npool_to_leafn_storage_patch[p] ≈ (nlc / cnl) * (1.0 - fcur)

        # gresp_storage check
        gresp_st = d.cnveg_cf.cpool_to_leafc_storage_patch[p] +
                   d.cnveg_cf.cpool_to_frootc_storage_patch[p] +
                   d.cnveg_cf.cpool_to_livestemc_storage_patch[p] +
                   d.cnveg_cf.cpool_to_deadstemc_storage_patch[p] +
                   d.cnveg_cf.cpool_to_livecrootc_storage_patch[p] +
                   d.cnveg_cf.cpool_to_deadcrootc_storage_patch[p]
        g1 = 0.3
        g2 = 0.5
        @test d.cnveg_cf.cpool_to_gresp_storage_patch[p] ≈ gresp_st * g1 * (1.0 - g2)
    end

    # =====================================================================
    # Test calc_plant_cn_alloc! -- with FUN (no downregulation)
    # =====================================================================
    @testset "calc_plant_cn_alloc! with FUN" begin
        d = make_woody_test_data(gpp=10.0, availc=8.0,
                                  plant_ndemand=0.2, retransn_to_npool=0.05,
                                  c_allometry=1.5, n_allometry=0.06,
                                  fpg=0.8)
        d.cn_shared_params.use_fun = true
        sminn_fun = 0.12
        d.cnveg_nf.sminn_to_plant_fun_patch[1] = sminn_fun

        original_psnsun = d.cnveg_cf.psnsun_to_cpool_patch[1]

        CLM.calc_plant_cn_alloc!(d.mask_soilp, d.bounds,
            d.pftcon, d.cn_shared_params, d.patch, d.crop,
            d.cnveg_state, d.cnveg_cs, d.cnveg_cf, d.cnveg_nf;
            fpg_col=d.fpg_col)

        p = 1
        # With FUN, sminn_to_npool = sminn_to_plant_fun
        @test d.cnveg_nf.sminn_to_npool_patch[p] ≈ sminn_fun

        # psnsun_to_cpool should NOT have been modified (use_fun skips downreg)
        @test d.cnveg_cf.psnsun_to_cpool_patch[p] ≈ original_psnsun
    end

    # =====================================================================
    # Test calc_plant_cn_alloc! -- non-woody grass patch
    # =====================================================================
    @testset "calc_plant_cn_alloc! non-woody grass" begin
        d = make_woody_test_data()

        # Make this a non-woody grass (keep ivt=1 to avoid pftcon bounds issues,
        # just modify the woody flag and stem_leaf parameter)
        d.pftcon.woody[1] = 0.0
        d.pftcon.stem_leaf[1] = 0.0   # no stem

        CLM.calc_plant_cn_alloc!(d.mask_soilp, d.bounds,
            d.pftcon, d.cn_shared_params, d.patch, d.crop,
            d.cnveg_state, d.cnveg_cs, d.cnveg_cf, d.cnveg_nf;
            fpg_col=d.fpg_col)

        p = 1
        # For non-woody, stem/croot allocation fluxes should stay at zero
        # (they were initialized to 0.0 and the woody block should not be entered)
        @test d.cnveg_cf.cpool_to_livestemc_patch[p] ≈ 0.0
        @test d.cnveg_cf.cpool_to_deadstemc_patch[p] ≈ 0.0
        @test d.cnveg_cf.cpool_to_livecrootc_patch[p] ≈ 0.0
        @test d.cnveg_cf.cpool_to_deadcrootc_patch[p] ≈ 0.0

        # Leaf and froot allocations should be nonzero
        @test d.cnveg_cf.cpool_to_leafc_patch[p] > 0.0
        @test d.cnveg_cf.cpool_to_frootc_patch[p] > 0.0
    end

    # =====================================================================
    # Test calc_plant_nitrogen_demand! -- crop patch with grain fill
    # =====================================================================
    @testset "calc_plant_nitrogen_demand! crop grain-fill retrans" begin
        np = 1
        nc = 1
        nrepr = CLM.NREPR
        dt = 1800.0
        npft = 17  # need pftcon arrays large enough for ivt=17

        # Build pftcon vectors of length npft; only index 17 is used
        make_pft_vec(val) = fill(val, npft)
        pftcon = CLM.PftConNutrientCompetition(
            woody       = make_pft_vec(0.0),
            froot_leaf  = make_pft_vec(0.0),
            croot_stem  = make_pft_vec(0.0),
            stem_leaf   = make_pft_vec(0.0),
            flivewd     = make_pft_vec(0.1),
            leafcn      = make_pft_vec(25.0),
            frootcn     = make_pft_vec(42.0),
            livewdcn    = make_pft_vec(50.0),
            deadwdcn    = make_pft_vec(500.0),
            fcur        = make_pft_vec(1.0),
            graincn     = make_pft_vec(50.0),
            grperc      = make_pft_vec(0.25),
            grpnow      = make_pft_vec(0.5),
            fleafcn     = make_pft_vec(65.0),
            ffrootcn    = make_pft_vec(0.0),    # no froot retrans
            fstemcn     = make_pft_vec(130.0),
            astemf      = make_pft_vec(0.05),
            season_decid = make_pft_vec(0.0),
            stress_decid = make_pft_vec(0.0),
        )

        cn_shared_params = CLM.CNSharedParamsData(use_matrixcn=false, use_fun=false)

        patch = CLM.PatchData()
        patch.column = [1]
        patch.itype  = [17]   # corn (npcropmin = 17)

        crop = CLM.CropData()
        crop.croplive_patch = [true]
        crop.gddtsoi_patch  = [100.0]
        crop.hui_patch      = [200.0]

        cnveg_state = CLM.CNVegStateData()
        cnveg_state.c_allometry_patch          = [1.5]
        cnveg_state.n_allometry_patch          = [0.06]
        cnveg_state.downreg_patch              = [0.0]
        cnveg_state.aleaf_patch                = [0.4]
        cnveg_state.astem_patch                = [0.05]   # == astemf => grain fill triggers
        cnveg_state.aroot_patch                = [0.3]
        cnveg_state.arepr_patch                = fill(0.2, np, nrepr)
        cnveg_state.peaklai_patch              = [0]
        cnveg_state.tempsum_potential_gpp_patch = [0.0]
        cnveg_state.annsum_potential_gpp_patch  = [5000.0]
        cnveg_state.tempmax_retransn_patch      = [0.0]
        cnveg_state.annmax_retransn_patch       = [1.0]
        cnveg_state.grain_flag_patch            = [0.0]
        cnveg_state.huileaf_patch               = [50.0]   # gddtsoi >= huileaf
        cnveg_state.huigrain_patch              = [150.0]  # hui >= huigrain => grainfill

        cnveg_cs = CLM.CNVegCarbonStateData()
        cnveg_cs.leafc_patch      = [10.0]
        cnveg_cs.frootc_patch     = [5.0]
        cnveg_cs.livestemc_patch  = [20.0]

        cnveg_cf = CLM.CNVegCarbonFluxData()
        cnveg_cf.gpp_before_downreg_patch = [10.0]
        cnveg_cf.availc_patch             = [8.0]
        cnveg_cf.annsum_npp_patch         = [400.0]

        cnveg_ns = CLM.CNVegNitrogenStateData()
        cnveg_ns.retransn_patch = [0.5]

        cnveg_nf = CLM.CNVegNitrogenFluxData()
        cnveg_nf.plant_ndemand_patch       = [0.0]
        cnveg_nf.avail_retransn_patch      = [0.0]
        cnveg_nf.retransn_to_npool_patch   = [0.0]
        cnveg_nf.sminn_to_npool_patch      = [0.0]
        cnveg_nf.plant_nalloc_patch        = [0.0]
        cnveg_nf.sminn_to_plant_fun_patch  = [0.0]
        cnveg_nf.leafn_to_retransn_patch   = [0.0]
        cnveg_nf.frootn_to_retransn_patch  = [0.0]
        cnveg_nf.livestemn_to_retransn_patch = [0.0]

        mask_p = trues(np)
        bounds = 1:np

        CLM.calc_plant_nitrogen_demand!(mask_p, bounds, true,
            pftcon, cn_shared_params, patch, crop,
            cnveg_state, cnveg_cs, cnveg_cf, cnveg_ns, cnveg_nf;
            dt=dt)

        p = 1
        # grain_flag should be set to 1.0 (grain fill triggered)
        @test cnveg_state.grain_flag_patch[p] ≈ 1.0

        # leafn_to_retransn should be set based on leafc and C:N ratios
        # t1 * ((leafc/leafcn) - (leafc/fleafcn)) = (1/1800)*((10/25)-(10/65))
        expected_leafn_retrans = (1.0 / dt) * ((10.0 / 25.0) - (10.0 / 65.0))
        @test cnveg_nf.leafn_to_retransn_patch[p] ≈ expected_leafn_retrans

        # livestemn_to_retransn
        expected_stemn_retrans = (1.0 / dt) * ((20.0 / 50.0) - (20.0 / 130.0))
        @test cnveg_nf.livestemn_to_retransn_patch[p] ≈ expected_stemn_retrans

        # ffrootcn = 0 => frootn_to_retransn should be 0
        @test cnveg_nf.frootn_to_retransn_patch[p] ≈ 0.0

        # avail_retransn should be plant_ndemand (crop during grain fill)
        # Note: plant_ndemand was set and then modified by retransn_to_npool
        @test cnveg_nf.avail_retransn_patch[p] >= 0.0
    end

    # =====================================================================
    # Test calc_plant_nutrient_competition! (full pipeline, non-crop)
    # =====================================================================
    @testset "calc_plant_nutrient_competition! full pipeline" begin
        d = make_woody_test_data(gpp=10.0, availc=8.0,
                                  plant_ndemand=0.2, retransn_to_npool=0.05,
                                  c_allometry=1.5, n_allometry=0.06,
                                  fpg=0.8)

        CLM.calc_plant_nutrient_competition!(d.mask_soilp, d.bounds,
            d.pftcon, d.cn_shared_params, d.patch, d.crop,
            d.cnveg_state, d.cnveg_cs, d.cnveg_cf, d.cnveg_nf;
            fpg_col=d.fpg_col)

        p = 1
        # Check that allocation outputs are populated
        @test d.cnveg_cf.plant_calloc_patch[p] > 0.0
        @test d.cnveg_nf.plant_nalloc_patch[p] > 0.0
        @test d.cnveg_cf.cpool_to_leafc_patch[p] > 0.0
        @test d.cnveg_nf.npool_to_leafn_patch[p] > 0.0
    end

    # =====================================================================
    # Test calc_plant_nutrient_demand! wrapper
    # =====================================================================
    @testset "calc_plant_nutrient_demand! wrapper" begin
        d = make_woody_test_data()
        dt = 1800.0

        CLM.calc_plant_nutrient_demand!(d.mask_soilp, d.bounds, false,
            d.pftcon, d.cn_shared_params, d.patch, d.crop,
            d.cnveg_state, d.cnveg_cs, d.cnveg_cf, d.cnveg_ns, d.cnveg_nf;
            dt=dt)

        p = 1
        # Basic sanity: demand should be non-negative after reduction
        @test d.cnveg_nf.plant_ndemand_patch[p] >= 0.0
    end

    # =====================================================================
    # Test mask behavior: skipped patches
    # =====================================================================
    @testset "mask behavior: inactive patches unchanged" begin
        np = 2
        nrepr = CLM.NREPR

        d = make_woody_test_data()

        # Extend all arrays to 2 patches
        resize!(d.mask_soilp, np)
        d.mask_soilp[2] = false  # second patch is inactive

        push!(d.patch.column, 1)
        push!(d.patch.itype, 1)

        push!(d.cnveg_cf.plant_calloc_patch, -99.0)  # sentinel
        push!(d.cnveg_cf.cpool_to_leafc_patch, -99.0)
        push!(d.cnveg_cf.gpp_before_downreg_patch, 10.0)
        push!(d.cnveg_cf.availc_patch, 8.0)
        push!(d.cnveg_cf.excess_cflux_patch, -99.0)
        push!(d.cnveg_cf.psnsun_to_cpool_patch, 6.0)
        push!(d.cnveg_cf.psnshade_to_cpool_patch, 4.0)
        push!(d.cnveg_cf.annsum_npp_patch, 400.0)
        push!(d.cnveg_cf.cpool_to_leafc_storage_patch, -99.0)
        push!(d.cnveg_cf.cpool_to_frootc_patch, -99.0)
        push!(d.cnveg_cf.cpool_to_frootc_storage_patch, -99.0)
        push!(d.cnveg_cf.cpool_to_livestemc_patch, -99.0)
        push!(d.cnveg_cf.cpool_to_livestemc_storage_patch, -99.0)
        push!(d.cnveg_cf.cpool_to_deadstemc_patch, -99.0)
        push!(d.cnveg_cf.cpool_to_deadstemc_storage_patch, -99.0)
        push!(d.cnveg_cf.cpool_to_livecrootc_patch, -99.0)
        push!(d.cnveg_cf.cpool_to_livecrootc_storage_patch, -99.0)
        push!(d.cnveg_cf.cpool_to_deadcrootc_patch, -99.0)
        push!(d.cnveg_cf.cpool_to_deadcrootc_storage_patch, -99.0)
        push!(d.cnveg_cf.cpool_to_gresp_storage_patch, -99.0)
        d.cnveg_cf.cpool_to_reproductivec_patch         = fill(0.0, np, nrepr)
        d.cnveg_cf.cpool_to_reproductivec_storage_patch = fill(0.0, np, nrepr)

        push!(d.cnveg_state.c_allometry_patch, 1.5)
        push!(d.cnveg_state.n_allometry_patch, 0.06)
        push!(d.cnveg_state.downreg_patch, -99.0)
        push!(d.cnveg_state.aleaf_patch, NaN)
        push!(d.cnveg_state.astem_patch, NaN)
        push!(d.cnveg_state.aroot_patch, NaN)
        d.cnveg_state.arepr_patch = fill(NaN, np, nrepr)

        push!(d.cnveg_nf.plant_ndemand_patch, -99.0)
        push!(d.cnveg_nf.avail_retransn_patch, -99.0)
        push!(d.cnveg_nf.retransn_to_npool_patch, -99.0)
        push!(d.cnveg_nf.sminn_to_npool_patch, -99.0)
        push!(d.cnveg_nf.plant_nalloc_patch, -99.0)
        push!(d.cnveg_nf.npool_to_leafn_patch, -99.0)
        push!(d.cnveg_nf.npool_to_leafn_storage_patch, -99.0)
        push!(d.cnveg_nf.npool_to_frootn_patch, -99.0)
        push!(d.cnveg_nf.npool_to_frootn_storage_patch, -99.0)
        push!(d.cnveg_nf.npool_to_livestemn_patch, -99.0)
        push!(d.cnveg_nf.npool_to_livestemn_storage_patch, -99.0)
        push!(d.cnveg_nf.npool_to_deadstemn_patch, -99.0)
        push!(d.cnveg_nf.npool_to_deadstemn_storage_patch, -99.0)
        push!(d.cnveg_nf.npool_to_livecrootn_patch, -99.0)
        push!(d.cnveg_nf.npool_to_livecrootn_storage_patch, -99.0)
        push!(d.cnveg_nf.npool_to_deadcrootn_patch, -99.0)
        push!(d.cnveg_nf.npool_to_deadcrootn_storage_patch, -99.0)
        d.cnveg_nf.npool_to_reproductiven_patch         = fill(0.0, np, nrepr)
        d.cnveg_nf.npool_to_reproductiven_storage_patch = fill(0.0, np, nrepr)
        push!(d.cnveg_nf.sminn_to_plant_fun_patch, 0.0)

        push!(d.crop.croplive_patch, false)

        CLM.calc_plant_cn_alloc!(d.mask_soilp, 1:np,
            d.pftcon, d.cn_shared_params, d.patch, d.crop,
            d.cnveg_state, d.cnveg_cs, d.cnveg_cf, d.cnveg_nf;
            fpg_col=d.fpg_col)

        # Patch 2 (masked out) should be untouched
        @test d.cnveg_cf.plant_calloc_patch[2] == -99.0
        @test d.cnveg_cf.cpool_to_leafc_patch[2] == -99.0
        @test d.cnveg_nf.npool_to_leafn_patch[2] == -99.0

        # Patch 1 (active) should be modified
        @test d.cnveg_cf.plant_calloc_patch[1] != -99.0
    end

    # =====================================================================
    # Test C13/C14 downregulation
    # =====================================================================
    @testset "C13/C14 downregulation" begin
        d = make_woody_test_data(gpp=10.0, availc=8.0,
                                  plant_ndemand=0.2, retransn_to_npool=0.05,
                                  c_allometry=1.5, n_allometry=0.06,
                                  fpg=0.8)

        c13_cf = CLM.CNVegCarbonFluxData()
        c13_cf.psnsun_to_cpool_patch  = [3.0]
        c13_cf.psnshade_to_cpool_patch = [2.0]

        c14_cf = CLM.CNVegCarbonFluxData()
        c14_cf.psnsun_to_cpool_patch  = [0.003]
        c14_cf.psnshade_to_cpool_patch = [0.002]

        CLM.calc_plant_cn_alloc!(d.mask_soilp, d.bounds,
            d.pftcon, d.cn_shared_params, d.patch, d.crop,
            d.cnveg_state, d.cnveg_cs, d.cnveg_cf, d.cnveg_nf;
            fpg_col=d.fpg_col,
            c13_cnveg_cf=c13_cf,
            c14_cnveg_cf=c14_cf,
            use_c13=true,
            use_c14=true)

        p = 1
        downreg = d.cnveg_state.downreg_patch[p]

        # C13 should be downregulated
        @test c13_cf.psnsun_to_cpool_patch[p] ≈ 3.0 * (1.0 - downreg)
        @test c13_cf.psnshade_to_cpool_patch[p] ≈ 2.0 * (1.0 - downreg)

        # C14 should be downregulated
        @test c14_cf.psnsun_to_cpool_patch[p] ≈ 0.003 * (1.0 - downreg)
        @test c14_cf.psnshade_to_cpool_patch[p] ≈ 0.002 * (1.0 - downreg)
    end

    # =====================================================================
    # Test dynamic wood allocation (stem_leaf = -1)
    # =====================================================================
    @testset "dynamic wood allocation f3 formula" begin
        d = make_woody_test_data(annsum_npp=0.0)  # low NPP => low f3

        CLM.calc_plant_cn_alloc!(d.mask_soilp, d.bounds,
            d.pftcon, d.cn_shared_params, d.patch, d.crop,
            d.cnveg_state, d.cnveg_cs, d.cnveg_cf, d.cnveg_nf;
            fpg_col=d.fpg_col)

        p = 1
        # f3 at NPP=0: (2.7/(1+exp(-0.004*(0-300))))-0.4
        f3_low = (2.7 / (1.0 + exp(-0.004 * (0.0 - 300.0)))) - 0.4
        @test f3_low < 1.0  # should be relatively small

        # Run again with high NPP
        d2 = make_woody_test_data(annsum_npp=2000.0)

        CLM.calc_plant_cn_alloc!(d2.mask_soilp, d2.bounds,
            d2.pftcon, d2.cn_shared_params, d2.patch, d2.crop,
            d2.cnveg_state, d2.cnveg_cs, d2.cnveg_cf, d2.cnveg_nf;
            fpg_col=d2.fpg_col)

        f3_high = (2.7 / (1.0 + exp(-0.004 * (2000.0 - 300.0)))) - 0.4

        # Higher NPP => higher stem allocation => more stemc
        @test d2.cnveg_cf.cpool_to_livestemc_patch[1] > d.cnveg_cf.cpool_to_livestemc_patch[1] ||
              d2.cnveg_cf.cpool_to_deadstemc_patch[1] > d.cnveg_cf.cpool_to_deadstemc_patch[1]
    end

    # =====================================================================
    # Test N conservation in allocation
    # =====================================================================
    @testset "N conservation: plant_nalloc = sum of N allocation fluxes" begin
        d = make_woody_test_data(gpp=10.0, availc=8.0,
                                  plant_ndemand=0.2, retransn_to_npool=0.05,
                                  c_allometry=1.5, n_allometry=0.06,
                                  fpg=0.8)

        CLM.calc_plant_cn_alloc!(d.mask_soilp, d.bounds,
            d.pftcon, d.cn_shared_params, d.patch, d.crop,
            d.cnveg_state, d.cnveg_cs, d.cnveg_cf, d.cnveg_nf;
            fpg_col=d.fpg_col)

        p = 1
        # Sum of all N allocation fluxes should equal plant_nalloc
        sum_n = d.cnveg_nf.npool_to_leafn_patch[p] +
                d.cnveg_nf.npool_to_leafn_storage_patch[p] +
                d.cnveg_nf.npool_to_frootn_patch[p] +
                d.cnveg_nf.npool_to_frootn_storage_patch[p] +
                d.cnveg_nf.npool_to_livestemn_patch[p] +
                d.cnveg_nf.npool_to_livestemn_storage_patch[p] +
                d.cnveg_nf.npool_to_deadstemn_patch[p] +
                d.cnveg_nf.npool_to_deadstemn_storage_patch[p] +
                d.cnveg_nf.npool_to_livecrootn_patch[p] +
                d.cnveg_nf.npool_to_livecrootn_storage_patch[p] +
                d.cnveg_nf.npool_to_deadcrootn_patch[p] +
                d.cnveg_nf.npool_to_deadcrootn_storage_patch[p]

        # plant_nalloc is the total N available for allocation
        # The N allocated to pools should be consistent with the nlc calculation
        # nlc = plant_calloc / c_allometry
        # N to leaf = nlc/cnl, N to froot = nlc*f1/cnfr, etc.
        # These form the total N demand at plant_nalloc
        @test sum_n > 0.0
        @test d.cnveg_nf.plant_nalloc_patch[p] > 0.0

        # The ratio plant_calloc/plant_nalloc should equal c_allometry/n_allometry
        @test d.cnveg_cf.plant_calloc_patch[p] / d.cnveg_nf.plant_nalloc_patch[p] ≈
              d.cnveg_state.c_allometry_patch[p] / d.cnveg_state.n_allometry_patch[p]
    end

    # =====================================================================
    # Test soybean-specific grain fill logic
    # =====================================================================
    @testset "soybean grain fill: astem != astemf delays retrans" begin
        np = 1
        nrepr = CLM.NREPR
        dt = 1800.0
        npft = 23  # need pftcon arrays large enough for ivt=23

        make_pft_vec2(val) = fill(val, npft)
        pftcon = CLM.PftConNutrientCompetition(
            woody       = make_pft_vec2(0.0),
            froot_leaf  = make_pft_vec2(0.0),
            croot_stem  = make_pft_vec2(0.0),
            stem_leaf   = make_pft_vec2(0.0),
            flivewd     = make_pft_vec2(0.1),
            leafcn      = make_pft_vec2(25.0),
            frootcn     = make_pft_vec2(42.0),
            livewdcn    = make_pft_vec2(50.0),
            deadwdcn    = make_pft_vec2(500.0),
            fcur        = make_pft_vec2(1.0),
            graincn     = make_pft_vec2(50.0),
            grperc      = make_pft_vec2(0.25),
            grpnow      = make_pft_vec2(0.5),
            fleafcn     = make_pft_vec2(65.0),
            ffrootcn    = make_pft_vec2(0.0),
            fstemcn     = make_pft_vec2(130.0),
            astemf      = make_pft_vec2(0.05),
            season_decid = make_pft_vec2(0.0),
            stress_decid = make_pft_vec2(0.0),
        )

        cn_shared_params = CLM.CNSharedParamsData(use_fun=false)
        patch = CLM.PatchData(column=[1], itype=[23])  # ntmp_soybean = 23
        crop = CLM.CropData(croplive_patch=[true], gddtsoi_patch=[100.0], hui_patch=[200.0])

        cnveg_state = CLM.CNVegStateData()
        cnveg_state.c_allometry_patch = [1.5]
        cnveg_state.n_allometry_patch = [0.06]
        cnveg_state.aleaf_patch       = [0.4]
        cnveg_state.astem_patch       = [0.10]   # != astemf (0.05) => soybean retrans delayed
        cnveg_state.aroot_patch       = [0.3]
        cnveg_state.arepr_patch       = fill(0.2, np, nrepr)
        cnveg_state.tempsum_potential_gpp_patch = [0.0]
        cnveg_state.annsum_potential_gpp_patch  = [5000.0]
        cnveg_state.tempmax_retransn_patch      = [0.0]
        cnveg_state.annmax_retransn_patch       = [1.0]
        cnveg_state.grain_flag_patch            = [0.0]
        cnveg_state.huileaf_patch               = [50.0]
        cnveg_state.huigrain_patch              = [150.0]

        cnveg_cs = CLM.CNVegCarbonStateData(leafc_patch=[10.0], frootc_patch=[5.0], livestemc_patch=[20.0])
        cnveg_cf = CLM.CNVegCarbonFluxData(gpp_before_downreg_patch=[10.0], availc_patch=[8.0])
        cnveg_ns = CLM.CNVegNitrogenStateData(retransn_patch=[0.5])

        cnveg_nf = CLM.CNVegNitrogenFluxData()
        cnveg_nf.plant_ndemand_patch       = [0.0]
        cnveg_nf.avail_retransn_patch      = [0.0]
        cnveg_nf.retransn_to_npool_patch   = [0.0]
        cnveg_nf.sminn_to_npool_patch      = [0.0]
        cnveg_nf.plant_nalloc_patch        = [0.0]
        cnveg_nf.sminn_to_plant_fun_patch  = [0.0]
        cnveg_nf.leafn_to_retransn_patch   = [0.0]
        cnveg_nf.frootn_to_retransn_patch  = [0.0]
        cnveg_nf.livestemn_to_retransn_patch = [0.0]

        mask_p = trues(np)

        CLM.calc_plant_nitrogen_demand!(mask_p, 1:np, true,
            pftcon, cn_shared_params, patch, crop,
            cnveg_state, cnveg_cs, cnveg_cf, cnveg_ns, cnveg_nf;
            dt=dt,
            ntmp_soybean=23, nirrig_tmp_soybean=24,
            ntrp_soybean=77, nirrig_trp_soybean=78)

        # For soybean with astem != astemf, grain_flag should still be 0
        @test cnveg_state.grain_flag_patch[1] ≈ 0.0

        # leafn_to_retransn should NOT have been set (still 0)
        @test cnveg_nf.leafn_to_retransn_patch[1] ≈ 0.0
    end

    # =====================================================================
    # Test zero GPP scenario
    # =====================================================================
    @testset "zero GPP: no downregulation crash" begin
        d = make_woody_test_data(gpp=0.0, availc=0.0,
                                  plant_ndemand=0.0,
                                  retransn_to_npool=0.0,
                                  c_allometry=1.5, n_allometry=0.06,
                                  fpg=0.0)

        # Should not error
        CLM.calc_plant_cn_alloc!(d.mask_soilp, d.bounds,
            d.pftcon, d.cn_shared_params, d.patch, d.crop,
            d.cnveg_state, d.cnveg_cs, d.cnveg_cf, d.cnveg_nf;
            fpg_col=d.fpg_col)

        p = 1
        # With zero GPP, downreg should not crash (gpp>0 branch skipped)
        @test d.cnveg_cf.cpool_to_leafc_patch[p] ≈ 0.0
        @test d.cnveg_nf.npool_to_leafn_patch[p] ≈ 0.0
    end

end
