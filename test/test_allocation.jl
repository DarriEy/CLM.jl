@testset "Allocation" begin

    # =====================================================================
    # AllocationParams struct
    # =====================================================================
    @testset "AllocationParams defaults" begin
        ap = CLM.AllocationParams()
        @test ap.dayscrecover ≈ 30.0
    end

    # =====================================================================
    # PftConAllocation struct
    # =====================================================================
    @testset "PftConAllocation defaults" begin
        pfc = CLM.PftConAllocation()
        @test isempty(pfc.woody)
        @test isempty(pfc.froot_leaf)
        @test isempty(pfc.croot_stem)
        @test isempty(pfc.stem_leaf)
        @test isempty(pfc.flivewd)
        @test isempty(pfc.leafcn)
        @test isempty(pfc.frootcn)
        @test isempty(pfc.livewdcn)
        @test isempty(pfc.deadwdcn)
        @test isempty(pfc.graincn)
        @test isempty(pfc.grperc)
        @test isempty(pfc.arooti)
        @test isempty(pfc.arootf)
        @test isempty(pfc.bfact)
        @test isempty(pfc.fleafi)
        @test isempty(pfc.aleaff)
        @test isempty(pfc.astemf)
        @test isempty(pfc.allconss)
        @test isempty(pfc.allconsl)
        @test isempty(pfc.declfact)
    end

    # =====================================================================
    # Helper: create test data for a single woody patch
    # =====================================================================
    function make_woody_test_data(;
            psnsun=5.0, psnsha=2.0,
            laisun=2.0, laisha=3.0,
            leaf_mr=0.001, froot_mr=0.0005,
            livestem_mr=0.0003, livecroot_mr=0.0002,
            xsmrpool=0.0)

        np = 1
        nrepr = CLM.NREPR

        # Allocation params
        alloc_params = CLM.AllocationParams(dayscrecover=30.0)

        # PFT constants: woody tree (ivt=1)
        pftcon = CLM.PftConAllocation(
            woody      = [1.0],
            froot_leaf = [1.0],
            croot_stem = [0.3],
            stem_leaf  = [-1.0],  # dynamic allocation
            flivewd    = [0.1],
            leafcn     = [25.0],
            frootcn    = [42.0],
            livewdcn   = [50.0],
            deadwdcn   = [500.0],
            graincn    = [50.0],
            grperc     = [0.3],
            arooti     = [0.0],
            arootf     = [0.0],
            bfact      = [0.0],
            fleafi     = [0.0],
            aleaff     = [0.0],
            astemf     = [0.0],
            allconss   = [0.0],
            allconsl   = [0.0],
            declfact   = [0.0],
        )

        cn_shared_params = CLM.CNSharedParamsData(use_matrixcn=false, use_fun=false)

        patch = CLM.PatchData()
        CLM.patch_init!(patch, np)
        patch.itype[1] = 1

        crop = CLM.CropData()
        CLM.crop_init!(crop, np)
        crop.croplive_patch[1] = false

        photosyns = CLM.PhotosynthesisData()
        CLM.photosynthesis_data_init!(photosyns, np)
        photosyns.psnsun_patch[1] = psnsun
        photosyns.psnsha_patch[1] = psnsha
        photosyns.c13_psnsun_patch[1] = 0.0
        photosyns.c13_psnsha_patch[1] = 0.0
        photosyns.c14_psnsun_patch[1] = 0.0
        photosyns.c14_psnsha_patch[1] = 0.0

        canopystate = CLM.CanopyStateData()
        CLM.canopystate_init!(canopystate, np)
        canopystate.laisun_patch[1] = laisun
        canopystate.laisha_patch[1] = laisha

        cnveg_cs = CLM.CNVegCarbonStateData()
        CLM.cnveg_carbon_state_init!(cnveg_cs, np, 1, 1)
        cnveg_cs.xsmrpool_patch[1] = xsmrpool

        cnveg_cf = CLM.CNVegCarbonFluxData()
        CLM.cnveg_carbon_flux_init!(cnveg_cf, np, 1, 1)
        cnveg_cf.leaf_mr_patch[1]      = leaf_mr
        cnveg_cf.froot_mr_patch[1]     = froot_mr
        cnveg_cf.livestem_mr_patch[1]  = livestem_mr
        cnveg_cf.livecroot_mr_patch[1] = livecroot_mr
        for k in 1:nrepr
            cnveg_cf.reproductive_mr_patch[1, k] = 0.0
        end
        cnveg_cf.annsum_npp_patch[1] = 400.0
        cnveg_cf.xsmrpool_recover_patch[1] = 0.0
        cnveg_cf.cpool_to_xsmrpool_patch[1] = 0.0

        mask = BitVector([true])
        bounds = 1:np

        return alloc_params, pftcon, cn_shared_params, patch, crop, photosyns,
               canopystate, cnveg_cs, cnveg_cf, mask, bounds
    end

    # =====================================================================
    # calc_gpp_mr_availc! — basic woody patch
    # =====================================================================
    @testset "calc_gpp_mr_availc! woody patch basic" begin
        alloc_params, pftcon, cn_shared_params, patch, crop, photosyns,
            canopystate, cnveg_cs, cnveg_cf, mask, bounds = make_woody_test_data(
                psnsun=5.0, psnsha=2.0, laisun=2.0, laisha=3.0,
                leaf_mr=0.001, froot_mr=0.0005,
                livestem_mr=0.0003, livecroot_mr=0.0002)

        CLM.calc_gpp_mr_availc!(mask, bounds,
            alloc_params, pftcon, cn_shared_params,
            patch, crop, photosyns, canopystate,
            cnveg_cs, cnveg_cf)

        # psnsun_to_cpool = 5.0 * 2.0 * 12.011e-6
        expected_psnsun = 5.0 * 2.0 * 12.011e-6
        expected_psnsha = 2.0 * 3.0 * 12.011e-6
        @test cnveg_cf.psnsun_to_cpool_patch[1] ≈ expected_psnsun
        @test cnveg_cf.psnshade_to_cpool_patch[1] ≈ expected_psnsha

        expected_gpp = expected_psnsun + expected_psnsha
        @test cnveg_cf.gpp_before_downreg_patch[1] ≈ expected_gpp

        # mr = leaf + froot + livestem + livecroot (woody)
        expected_mr = 0.001 + 0.0005 + 0.0003 + 0.0002
        expected_availc = max(expected_gpp - expected_mr, 0.0)
        @test cnveg_cf.availc_patch[1] ≈ expected_availc

        # Check curmr partitioning
        if expected_gpp >= expected_mr
            # curmr_ratio = 1.0
            @test cnveg_cf.leaf_curmr_patch[1] ≈ 0.001
            @test cnveg_cf.leaf_xsmr_patch[1] ≈ 0.0
        else
            # curmr_ratio = gpp/mr
            curmr_ratio = expected_gpp / expected_mr
            @test cnveg_cf.leaf_curmr_patch[1] ≈ 0.001 * curmr_ratio
            @test cnveg_cf.leaf_xsmr_patch[1] ≈ 0.001 - 0.001 * curmr_ratio
        end
    end

    # =====================================================================
    # calc_gpp_mr_availc! — MR exceeds GPP
    # =====================================================================
    @testset "calc_gpp_mr_availc! MR exceeds GPP" begin
        alloc_params, pftcon, cn_shared_params, patch, crop, photosyns,
            canopystate, cnveg_cs, cnveg_cf, mask, bounds = make_woody_test_data(
                psnsun=0.01, psnsha=0.01, laisun=1.0, laisha=1.0,
                leaf_mr=0.01, froot_mr=0.005,
                livestem_mr=0.003, livecroot_mr=0.002)

        CLM.calc_gpp_mr_availc!(mask, bounds,
            alloc_params, pftcon, cn_shared_params,
            patch, crop, photosyns, canopystate,
            cnveg_cs, cnveg_cf)

        # GPP = (0.01 * 1 + 0.01 * 1) * 12.011e-6 = very small
        gpp = 0.01 * 1.0 * 12.011e-6 + 0.01 * 1.0 * 12.011e-6
        mr  = 0.01 + 0.005 + 0.003 + 0.002  # = 0.02

        # mr > gpp, so availc < 0 → curmr_ratio = gpp/mr
        @test cnveg_cf.gpp_before_downreg_patch[1] ≈ gpp
        @test cnveg_cf.availc_patch[1] ≈ 0.0  # clamped to 0

        curmr_ratio = gpp / mr
        @test cnveg_cf.leaf_curmr_patch[1] ≈ 0.01 * curmr_ratio
        @test cnveg_cf.leaf_xsmr_patch[1] ≈ 0.01 - 0.01 * curmr_ratio
    end

    # =====================================================================
    # calc_gpp_mr_availc! — xsmrpool deficit recovery
    # =====================================================================
    @testset "calc_gpp_mr_availc! xsmrpool deficit" begin
        alloc_params, pftcon, cn_shared_params, patch, crop, photosyns,
            canopystate, cnveg_cs, cnveg_cf, mask, bounds = make_woody_test_data(
                psnsun=5.0, psnsha=2.0, laisun=2.0, laisha=3.0,
                leaf_mr=0.0, froot_mr=0.0,
                livestem_mr=0.0, livecroot_mr=0.0,
                xsmrpool=-100.0)

        CLM.calc_gpp_mr_availc!(mask, bounds,
            alloc_params, pftcon, cn_shared_params,
            patch, crop, photosyns, canopystate,
            cnveg_cs, cnveg_cf)

        # xsmrpool_recover = 100 / (30 * 86400)
        expected_recover = 100.0 / (30.0 * CLM.SECSPDAY)
        gpp = 5.0 * 2.0 * 12.011e-6 + 2.0 * 3.0 * 12.011e-6

        if expected_recover < gpp
            @test cnveg_cf.xsmrpool_recover_patch[1] ≈ expected_recover
            @test cnveg_cf.availc_patch[1] ≈ gpp - expected_recover
        else
            @test cnveg_cf.xsmrpool_recover_patch[1] ≈ gpp
            @test cnveg_cf.availc_patch[1] ≈ 0.0
        end
        @test cnveg_cf.cpool_to_xsmrpool_patch[1] ≈ cnveg_cf.xsmrpool_recover_patch[1]
    end

    # =====================================================================
    # calc_allometry! — woody patch
    # =====================================================================
    @testset "calc_allometry! woody patch" begin
        np = 1
        nrepr = CLM.NREPR

        pftcon = CLM.PftConAllocation(
            woody      = [1.0],
            froot_leaf = [1.0],
            croot_stem = [0.3],
            stem_leaf  = [2.0],
            flivewd    = [0.1],
            leafcn     = [25.0],
            frootcn    = [42.0],
            livewdcn   = [50.0],
            deadwdcn   = [500.0],
            graincn    = [50.0],
            grperc     = [0.3],
            arooti     = [0.0],
            arootf     = [0.0],
            bfact      = [0.0],
            fleafi     = [0.0],
            aleaff     = [0.0],
            astemf     = [0.0],
            allconss   = [0.0],
            allconsl   = [0.0],
            declfact   = [0.0],
        )

        cn_shared_params = CLM.CNSharedParamsData(use_fun=false)

        patch = CLM.PatchData()
        CLM.patch_init!(patch, np)
        patch.itype[1] = 1

        cnveg_cf = CLM.CNVegCarbonFluxData()
        CLM.cnveg_carbon_flux_init!(cnveg_cf, np, 1, 1)
        cnveg_cf.annsum_npp_patch[1] = 400.0

        cnveg_state = CLM.CNVegStateData()
        CLM.cnveg_state_init!(cnveg_state, np, 1)
        cnveg_state.aleaf_patch[1] = 0.5
        cnveg_state.astem_patch[1] = 0.2
        cnveg_state.aroot_patch[1] = 0.3
        for k in 1:nrepr
            cnveg_state.arepr_patch[1, k] = 0.0
        end

        mask = BitVector([true])
        bounds = 1:np

        CLM.calc_allometry!(mask, bounds, pftcon, cn_shared_params, patch,
            cnveg_cf, cnveg_state)

        # For woody: f1=1.0, f2=0.3, f3=2.0, f4=0.1, g1=0.3
        f1 = 1.0; f2 = 0.3; f3 = 2.0; f4 = 0.1; g1 = 0.3
        expected_c = (1.0 + g1) * (1.0 + f1 + f3 * (1.0 + f2))
        expected_n = 1.0/25.0 + f1/42.0 + (f3*f4*(1.0+f2))/50.0 + (f3*(1.0-f4)*(1.0+f2))/500.0
        @test cnveg_state.c_allometry_patch[1] ≈ expected_c
        @test cnveg_state.n_allometry_patch[1] ≈ expected_n
    end

    # =====================================================================
    # calc_allometry! — non-woody, non-crop patch
    # =====================================================================
    @testset "calc_allometry! grass patch" begin
        np = 1
        nrepr = CLM.NREPR
        npft = 14  # C4 grass = PFT 14

        pftcon = CLM.PftConAllocation(
            woody      = fill(0.0, npft),
            froot_leaf = fill(1.5, npft),
            croot_stem = fill(0.0, npft),
            stem_leaf  = fill(0.0, npft),
            flivewd    = fill(0.0, npft),
            leafcn     = fill(25.0, npft),
            frootcn    = fill(42.0, npft),
            livewdcn   = fill(50.0, npft),
            deadwdcn   = fill(500.0, npft),
            graincn    = fill(50.0, npft),
            grperc     = fill(0.3, npft),
            arooti     = fill(0.0, npft),
            arootf     = fill(0.0, npft),
            bfact      = fill(0.0, npft),
            fleafi     = fill(0.0, npft),
            aleaff     = fill(0.0, npft),
            astemf     = fill(0.0, npft),
            allconss   = fill(0.0, npft),
            allconsl   = fill(0.0, npft),
            declfact   = fill(0.0, npft),
        )

        cn_shared_params = CLM.CNSharedParamsData(use_fun=false)

        patch = CLM.PatchData()
        CLM.patch_init!(patch, np)
        patch.itype[1] = 14  # C4 grass

        cnveg_cf = CLM.CNVegCarbonFluxData()
        CLM.cnveg_carbon_flux_init!(cnveg_cf, np, 1, 1)
        cnveg_cf.annsum_npp_patch[1] = 200.0

        cnveg_state = CLM.CNVegStateData()
        CLM.cnveg_state_init!(cnveg_state, np, 1)

        mask = BitVector([true])
        bounds = 1:np

        CLM.calc_allometry!(mask, bounds, pftcon, cn_shared_params, patch,
            cnveg_cf, cnveg_state)

        # Non-woody, non-crop: c_allometry = 1 + g1a + f1 + f1*g1a
        f1 = 1.5; g1 = 0.3
        expected_c = 1.0 + g1 + f1 + f1 * g1
        expected_n = 1.0/25.0 + f1/42.0
        @test cnveg_state.c_allometry_patch[1] ≈ expected_c
        @test cnveg_state.n_allometry_patch[1] ≈ expected_n
    end

    # =====================================================================
    # calc_allometry! — use_fun = true (g1a = 0)
    # =====================================================================
    @testset "calc_allometry! with use_fun" begin
        np = 1
        nrepr = CLM.NREPR

        pftcon = CLM.PftConAllocation(
            woody      = [1.0],
            froot_leaf = [1.0],
            croot_stem = [0.3],
            stem_leaf  = [2.0],
            flivewd    = [0.1],
            leafcn     = [25.0],
            frootcn    = [42.0],
            livewdcn   = [50.0],
            deadwdcn   = [500.0],
            graincn    = [50.0],
            grperc     = [0.3],
            arooti     = [0.0],
            arootf     = [0.0],
            bfact      = [0.0],
            fleafi     = [0.0],
            aleaff     = [0.0],
            astemf     = [0.0],
            allconss   = [0.0],
            allconsl   = [0.0],
            declfact   = [0.0],
        )

        cn_shared_params = CLM.CNSharedParamsData(use_fun=true)

        patch = CLM.PatchData()
        CLM.patch_init!(patch, np)
        patch.itype[1] = 1

        cnveg_cf = CLM.CNVegCarbonFluxData()
        CLM.cnveg_carbon_flux_init!(cnveg_cf, np, 1, 1)
        cnveg_cf.annsum_npp_patch[1] = 400.0

        cnveg_state = CLM.CNVegStateData()
        CLM.cnveg_state_init!(cnveg_state, np, 1)

        mask = BitVector([true])
        bounds = 1:np

        CLM.calc_allometry!(mask, bounds, pftcon, cn_shared_params, patch,
            cnveg_cf, cnveg_state)

        # use_fun=true => g1a=0
        f1 = 1.0; f2 = 0.3; f3 = 2.0; f4 = 0.1
        expected_c = (1.0 + 0.0) * (1.0 + f1 + f3 * (1.0 + f2))
        @test cnveg_state.c_allometry_patch[1] ≈ expected_c
    end

    # =====================================================================
    # calc_allometry! — dynamic stem_leaf (=-1)
    # =====================================================================
    @testset "calc_allometry! dynamic stem_leaf" begin
        np = 1
        nrepr = CLM.NREPR

        pftcon = CLM.PftConAllocation(
            woody      = [1.0],
            froot_leaf = [1.0],
            croot_stem = [0.3],
            stem_leaf  = [-1.0],
            flivewd    = [0.1],
            leafcn     = [25.0],
            frootcn    = [42.0],
            livewdcn   = [50.0],
            deadwdcn   = [500.0],
            graincn    = [50.0],
            grperc     = [0.3],
            arooti     = [0.0],
            arootf     = [0.0],
            bfact      = [0.0],
            fleafi     = [0.0],
            aleaff     = [0.0],
            astemf     = [0.0],
            allconss   = [0.0],
            allconsl   = [0.0],
            declfact   = [0.0],
        )

        cn_shared_params = CLM.CNSharedParamsData(use_fun=false)

        patch = CLM.PatchData()
        CLM.patch_init!(patch, np)
        patch.itype[1] = 1

        cnveg_cf = CLM.CNVegCarbonFluxData()
        CLM.cnveg_carbon_flux_init!(cnveg_cf, np, 1, 1)
        cnveg_cf.annsum_npp_patch[1] = 800.0

        cnveg_state = CLM.CNVegStateData()
        CLM.cnveg_state_init!(cnveg_state, np, 1)

        mask = BitVector([true])
        bounds = 1:np

        CLM.calc_allometry!(mask, bounds, pftcon, cn_shared_params, patch,
            cnveg_cf, cnveg_state)

        # f3 = (2.7/(1+exp(-0.004*(800-300)))) - 0.4
        f3 = (2.7 / (1.0 + exp(-0.004 * (800.0 - 300.0)))) - 0.4
        f1 = 1.0; f2 = 0.3; g1 = 0.3
        expected_c = (1.0 + g1) * (1.0 + f1 + f3 * (1.0 + f2))
        @test cnveg_state.c_allometry_patch[1] ≈ expected_c
    end

    # =====================================================================
    # calc_crop_allocation_fractions! — not live crop
    # =====================================================================
    @testset "calc_crop_allocation_fractions! not live" begin
        np = 1
        nrepr = CLM.NREPR
        npft = 17  # crop PFT

        pftcon = CLM.PftConAllocation(
            woody      = fill(0.0, npft),
            froot_leaf = fill(0.0, npft),
            croot_stem = fill(0.0, npft),
            stem_leaf  = fill(0.0, npft),
            flivewd    = fill(0.0, npft),
            leafcn     = fill(25.0, npft),
            frootcn    = fill(42.0, npft),
            livewdcn   = fill(50.0, npft),
            deadwdcn   = fill(500.0, npft),
            graincn    = fill(50.0, npft),
            grperc     = fill(0.25, npft),
            arooti     = fill(0.6, npft),
            arootf     = fill(0.1, npft),
            bfact      = fill(0.01, npft),
            fleafi     = fill(0.85, npft),
            aleaff     = fill(0.0, npft),
            astemf     = fill(0.0, npft),
            allconss   = fill(1.0, npft),
            allconsl   = fill(1.0, npft),
            declfact   = fill(1.0, npft),
        )

        patch = CLM.PatchData()
        CLM.patch_init!(patch, np)
        patch.itype[1] = 17  # crop

        crop = CLM.CropData()
        CLM.crop_init!(crop, np)
        crop.croplive_patch[1] = false

        cnveg_state = CLM.CNVegStateData()
        CLM.cnveg_state_init!(cnveg_state, np, 1)

        mask = BitVector([true])
        bounds = 1:np

        CLM.calc_crop_allocation_fractions!(mask, bounds, pftcon, patch, crop, cnveg_state)

        @test cnveg_state.aleaf_patch[1] ≈ 1.0
        @test cnveg_state.astem_patch[1] ≈ 0.0
        @test cnveg_state.aroot_patch[1] ≈ 0.0
        for k in 1:nrepr
            @test cnveg_state.arepr_patch[1, k] ≈ 0.0
        end
    end

    # =====================================================================
    # calc_crop_allocation_fractions! — planted phase (pre-emergence)
    # =====================================================================
    @testset "calc_crop_allocation_fractions! planted phase" begin
        np = 1
        nrepr = CLM.NREPR
        npft = 17

        pftcon = CLM.PftConAllocation(
            woody      = fill(0.0, npft),
            froot_leaf = fill(0.0, npft),
            croot_stem = fill(0.0, npft),
            stem_leaf  = fill(0.0, npft),
            flivewd    = fill(0.0, npft),
            leafcn     = fill(25.0, npft),
            frootcn    = fill(42.0, npft),
            livewdcn   = fill(50.0, npft),
            deadwdcn   = fill(500.0, npft),
            graincn    = fill(50.0, npft),
            grperc     = fill(0.25, npft),
            arooti     = fill(0.6, npft),
            arootf     = fill(0.1, npft),
            bfact      = fill(0.01, npft),
            fleafi     = fill(0.85, npft),
            aleaff     = fill(0.0, npft),
            astemf     = fill(0.0, npft),
            allconss   = fill(1.0, npft),
            allconsl   = fill(1.0, npft),
            declfact   = fill(1.0, npft),
        )

        patch = CLM.PatchData()
        CLM.patch_init!(patch, np)
        patch.itype[1] = 17

        crop = CLM.CropData()
        CLM.crop_init!(crop, np)
        crop.croplive_patch[1] = true
        crop.hui_patch[1]      = 10.0
        crop.gddtsoi_patch[1]  = 5.0  # below huileaf -> cphase_planted

        cnveg_state = CLM.CNVegStateData()
        CLM.cnveg_state_init!(cnveg_state, np, 1)
        cnveg_state.huileaf_patch[1]    = 100.0
        cnveg_state.huigrain_patch[1]   = 500.0
        cnveg_state.gddmaturity_patch[1] = 1000.0

        mask = BitVector([true])
        bounds = 1:np

        CLM.calc_crop_allocation_fractions!(mask, bounds, pftcon, patch, crop, cnveg_state)

        # Pre-emergence: aleaf=1, astem=0, aroot=0
        @test cnveg_state.aleaf_patch[1] ≈ 1.0
        @test cnveg_state.astem_patch[1] ≈ 0.0
        @test cnveg_state.aroot_patch[1] ≈ 0.0
    end

    # =====================================================================
    # calc_crop_allocation_fractions! — leaf emerge phase
    # =====================================================================
    @testset "calc_crop_allocation_fractions! leaf emerge phase" begin
        np = 1
        nrepr = CLM.NREPR
        npft = 17

        pftcon = CLM.PftConAllocation(
            woody      = fill(0.0, npft),
            froot_leaf = fill(0.0, npft),
            croot_stem = fill(0.0, npft),
            stem_leaf  = fill(0.0, npft),
            flivewd    = fill(0.0, npft),
            leafcn     = fill(25.0, npft),
            frootcn    = fill(42.0, npft),
            livewdcn   = fill(50.0, npft),
            deadwdcn   = fill(500.0, npft),
            graincn    = fill(50.0, npft),
            grperc     = fill(0.25, npft),
            arooti     = fill(0.6, npft),
            arootf     = fill(0.1, npft),
            bfact      = fill(0.01, npft),
            fleafi     = fill(0.85, npft),
            aleaff     = fill(0.0, npft),
            astemf     = fill(0.0, npft),
            allconss   = fill(1.0, npft),
            allconsl   = fill(1.0, npft),
            declfact   = fill(1.0, npft),
        )

        patch = CLM.PatchData()
        CLM.patch_init!(patch, np)
        patch.itype[1] = 17

        crop = CLM.CropData()
        CLM.crop_init!(crop, np)
        crop.croplive_patch[1] = true
        crop.hui_patch[1]      = 200.0
        crop.gddtsoi_patch[1]  = 150.0  # >= huileaf

        cnveg_state = CLM.CNVegStateData()
        CLM.cnveg_state_init!(cnveg_state, np, 1)
        cnveg_state.huileaf_patch[1]    = 100.0
        cnveg_state.huigrain_patch[1]   = 500.0
        cnveg_state.gddmaturity_patch[1] = 1000.0
        cnveg_state.peaklai_patch[1]    = 0

        mask = BitVector([true])
        bounds = 1:np

        CLM.calc_crop_allocation_fractions!(mask, bounds, pftcon, patch, crop, cnveg_state)

        # Check that allocations sum to ~1
        total = cnveg_state.aleaf_patch[1] + cnveg_state.astem_patch[1] + cnveg_state.aroot_patch[1]
        for k in 1:nrepr
            total += cnveg_state.arepr_patch[1, k]
        end
        @test total ≈ 1.0 atol=1e-10

        # Check aroot is between arooti and arootf
        @test cnveg_state.aroot_patch[1] >= 0.0
        @test cnveg_state.aroot_patch[1] <= 1.0
        @test cnveg_state.aleaf_patch[1] >= 1.0e-5

        # astemi and aleafi should be saved
        @test cnveg_state.astemi_patch[1] ≈ cnveg_state.astem_patch[1]
        @test cnveg_state.aleafi_patch[1] ≈ cnveg_state.aleaf_patch[1]
    end

    # =====================================================================
    # calc_gpp_mr_availc! — masked patches are skipped
    # =====================================================================
    @testset "calc_gpp_mr_availc! respects mask" begin
        alloc_params, pftcon, cn_shared_params, patch, crop, photosyns,
            canopystate, cnveg_cs, cnveg_cf, _, bounds = make_woody_test_data()

        # Set mask to false — patch should not be touched
        mask = BitVector([false])
        cnveg_cf.psnsun_to_cpool_patch[1] = -999.0  # sentinel

        CLM.calc_gpp_mr_availc!(mask, bounds,
            alloc_params, pftcon, cn_shared_params,
            patch, crop, photosyns, canopystate,
            cnveg_cs, cnveg_cf)

        @test cnveg_cf.psnsun_to_cpool_patch[1] ≈ -999.0  # unchanged
    end

end
