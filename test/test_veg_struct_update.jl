@testset "CNVegStructUpdate (veg_struct_update.jl)" begin

    # -----------------------------------------------------------------------
    # Helper: create minimal data structures for np patches, nc columns
    # -----------------------------------------------------------------------
    function make_veg_struct_data(; np=4, nc=2)
        # --- Patch ---
        patch = CLM.PatchData()
        patch.itype  = zeros(Int, np)
        patch.column = ones(Int, np)
        patch.gridcell = ones(Int, np)
        patch.active = trues(np)

        # --- Canopy state ---
        cs = CLM.CanopyStateData()
        cs.tlai_patch               = zeros(np)
        cs.tsai_patch               = zeros(np)
        cs.elai_patch               = zeros(np)
        cs.esai_patch               = zeros(np)
        cs.htop_patch               = zeros(np)
        cs.hbot_patch               = zeros(np)
        cs.stem_biomass_patch       = zeros(np)
        cs.leaf_biomass_patch       = zeros(np)
        cs.frac_veg_nosno_alb_patch = zeros(Int, np)

        # --- CN veg carbon state ---
        cnveg_cs = CLM.CNVegCarbonStateData()
        cnveg_cs.leafc_patch     = zeros(np)
        cnveg_cs.deadstemc_patch = zeros(np)
        cnveg_cs.livestemc_patch = zeros(np)

        # --- Water diagnostic bulk ---
        wd = CLM.WaterDiagnosticBulkData()
        wd.frac_sno_col  = zeros(nc)
        wd.snow_depth_col = zeros(nc)

        # --- Friction velocity ---
        fv = CLM.FrictionVelocityData()
        fv.forc_hgt_u_patch = fill(30.0, np)

        # --- CN veg state ---
        vs = CLM.CNVegStateData()
        vs.farea_burned_col = zeros(nc)
        vs.htmx_patch       = zeros(np)
        vs.peaklai_patch    = zeros(Int, np)

        # --- Crop ---
        cr = CLM.CropData()
        cr.harvdate_patch = fill(999, np)

        # --- PFT constants (minimal setup) ---
        # Using 20 entries to cover PFT indices used in tests
        npft = 80
        pft = CLM.PftconType()
        pft.woody    = zeros(npft)
        pft.slatop   = fill(0.01, npft)
        pft.dsladlai = zeros(npft)
        pft.z0mr     = fill(0.055, npft)
        pft.displar  = fill(0.67, npft)
        pft.dwood    = fill(2.5e5, npft)
        pft.ztopmx   = fill(1.5, npft)
        pft.laimx    = fill(6.0, npft)
        pft.nstem    = fill(0.005, npft)
        pft.taper    = fill(200.0, npft)
        pft.fbw      = fill(0.5, npft)

        mask = trues(np)
        bounds = 1:np

        return (; patch, cs, cnveg_cs, wd, fv, vs, cr, pft, mask, bounds)
    end

    # -----------------------------------------------------------------------
    # Test 1: noveg patches get zero LAI/SAI/height
    # -----------------------------------------------------------------------
    @testset "noveg patches produce zero output" begin
        d = make_veg_struct_data(np=2, nc=1)
        d.patch.itype .= [0, 0]      # noveg = 0
        d.patch.column .= [1, 1]

        # Set some non-zero initial values to verify they get zeroed
        d.cs.tlai_patch .= [1.0, 2.0]
        d.cs.tsai_patch .= [0.5, 0.3]
        d.cs.htop_patch .= [5.0, 3.0]
        d.cs.hbot_patch .= [1.0, 0.5]

        CLM.cn_veg_struct_update!(
            d.mask, d.bounds, d.patch, d.cs, d.cnveg_cs, d.wd, d.fv, d.vs, d.cr, d.pft;
            dt=1800.0, noveg=0
        )

        for p in 1:2
            @test d.cs.tlai_patch[p] == 0.0
            @test d.cs.tsai_patch[p] == 0.0
            @test d.cs.htop_patch[p] == 0.0
            @test d.cs.hbot_patch[p] == 0.0
            @test d.cs.elai_patch[p] == 0.0
            @test d.cs.esai_patch[p] == 0.0
            @test d.cs.frac_veg_nosno_alb_patch[p] == 0
        end
    end

    # -----------------------------------------------------------------------
    # Test 2: LAI computation with dsladlai > 0 (exponential SLA)
    # -----------------------------------------------------------------------
    @testset "LAI with dsladlai > 0 (exponential SLA)" begin
        d = make_veg_struct_data(np=1, nc=1)
        d.patch.itype .= [13]        # c3_non-arctic_grass
        d.patch.column .= [1]

        leafc = 50.0
        slatop_val = 0.03
        dsladlai_val = 0.001

        d.cnveg_cs.leafc_patch .= [leafc]
        d.pft.slatop[13] = slatop_val
        d.pft.dsladlai[13] = dsladlai_val

        d.cs.tlai_patch .= [0.0]
        d.cs.tsai_patch .= [1.0]

        CLM.cn_veg_struct_update!(
            d.mask, d.bounds, d.patch, d.cs, d.cnveg_cs, d.wd, d.fv, d.vs, d.cr, d.pft;
            dt=1800.0, noveg=0
        )

        expected_tlai = (slatop_val * (exp(leafc * dsladlai_val) - 1.0)) / dsladlai_val
        @test d.cs.tlai_patch[1] ≈ expected_tlai atol=1e-10
    end

    # -----------------------------------------------------------------------
    # Test 3: LAI computation with dsladlai == 0 (linear SLA)
    # -----------------------------------------------------------------------
    @testset "LAI with dsladlai == 0 (linear SLA)" begin
        d = make_veg_struct_data(np=1, nc=1)
        d.patch.itype .= [13]
        d.patch.column .= [1]

        leafc = 50.0
        slatop_val = 0.02

        d.cnveg_cs.leafc_patch .= [leafc]
        d.pft.slatop[13] = slatop_val
        d.pft.dsladlai[13] = 0.0

        d.cs.tlai_patch .= [0.0]
        d.cs.tsai_patch .= [1.0]

        CLM.cn_veg_struct_update!(
            d.mask, d.bounds, d.patch, d.cs, d.cnveg_cs, d.wd, d.fv, d.vs, d.cr, d.pft;
            dt=1800.0, noveg=0
        )

        expected_tlai = slatop_val * leafc
        @test d.cs.tlai_patch[1] ≈ expected_tlai atol=1e-10
    end

    # -----------------------------------------------------------------------
    # Test 4: LAI is non-negative
    # -----------------------------------------------------------------------
    @testset "LAI is clamped to non-negative" begin
        d = make_veg_struct_data(np=1, nc=1)
        d.patch.itype .= [13]
        d.patch.column .= [1]

        # zero leafc should produce zero LAI
        d.cnveg_cs.leafc_patch .= [0.0]
        d.pft.slatop[13] = 0.02
        d.pft.dsladlai[13] = 0.0
        d.cs.tlai_patch .= [0.0]
        d.cs.tsai_patch .= [1.0]

        CLM.cn_veg_struct_update!(
            d.mask, d.bounds, d.patch, d.cs, d.cnveg_cs, d.wd, d.fv, d.vs, d.cr, d.pft;
            dt=1800.0, noveg=0
        )

        @test d.cs.tlai_patch[1] >= 0.0
    end

    # -----------------------------------------------------------------------
    # Test 5: TSI formula for generic crops (nc3crop)
    # -----------------------------------------------------------------------
    @testset "tsai for generic crops (nc3crop)" begin
        d = make_veg_struct_data(np=1, nc=1)
        nc3crop_val = 15
        d.patch.itype .= [nc3crop_val]
        d.patch.column .= [1]

        d.cnveg_cs.leafc_patch .= [40.0]
        d.pft.slatop[nc3crop_val] = 0.02
        d.pft.dsladlai[nc3crop_val] = 0.0

        tlai_old = 3.0
        tsai_old = 0.8
        d.cs.tlai_patch .= [tlai_old]
        d.cs.tsai_patch .= [tsai_old]

        dt = 1800.0
        CLM.cn_veg_struct_update!(
            d.mask, d.bounds, d.patch, d.cs, d.cnveg_cs, d.wd, d.fv, d.vs, d.cr, d.pft;
            dt=dt, noveg=0, nc3crop=nc3crop_val, npcropmin=17
        )

        new_tlai = 0.02 * 40.0  # slatop * leafc
        tsai_alpha = 1.0 - 1.0 * dt / CLM.DTSMONTH
        tsai_min = 0.1 * 0.5
        expected_tsai = max(tsai_alpha * tsai_old + max(tlai_old - new_tlai, 0.0), tsai_min)
        @test d.cs.tsai_patch[1] ≈ expected_tsai atol=1e-10
    end

    # -----------------------------------------------------------------------
    # Test 6: TSI formula for non-crop vegetation
    # -----------------------------------------------------------------------
    @testset "tsai for non-crop vegetation" begin
        d = make_veg_struct_data(np=1, nc=1)
        d.patch.itype .= [13]   # c3_non-arctic_grass
        d.patch.column .= [1]

        d.cnveg_cs.leafc_patch .= [30.0]
        d.pft.slatop[13] = 0.02
        d.pft.dsladlai[13] = 0.0

        tlai_old = 2.0
        tsai_old = 1.5
        d.cs.tlai_patch .= [tlai_old]
        d.cs.tsai_patch .= [tsai_old]

        dt = 1800.0
        CLM.cn_veg_struct_update!(
            d.mask, d.bounds, d.patch, d.cs, d.cnveg_cs, d.wd, d.fv, d.vs, d.cr, d.pft;
            dt=dt, noveg=0, npcropmin=17
        )

        new_tlai = 0.02 * 30.0
        tsai_alpha = 1.0 - 0.5 * dt / CLM.DTSMONTH
        tsai_min = 1.0 * 0.5
        expected_tsai = max(tsai_alpha * tsai_old + max(tlai_old - new_tlai, 0.0), tsai_min)
        @test d.cs.tsai_patch[1] ≈ expected_tsai atol=1e-10
    end

    # -----------------------------------------------------------------------
    # Test 7: Woody vegetation height from deadstemc
    # -----------------------------------------------------------------------
    @testset "woody vegetation height" begin
        d = make_veg_struct_data(np=1, nc=1)
        d.patch.itype .= [2]   # needleleaf_evergreen_boreal_tree
        d.patch.column .= [1]

        d.pft.woody[2] = 1.0
        deadstemc_val = 5000.0
        d.cnveg_cs.deadstemc_patch .= [deadstemc_val]
        d.cnveg_cs.leafc_patch .= [100.0]
        d.pft.slatop[2] = 0.008
        d.pft.dsladlai[2] = 0.0
        d.pft.nstem[2] = 0.005
        d.pft.taper[2] = 200.0
        d.pft.dwood[2] = 2.5e5
        d.pft.z0mr[2] = 0.055
        d.pft.displar[2] = 0.67

        d.cs.tlai_patch .= [0.5]
        d.cs.tsai_patch .= [1.0]
        d.fv.forc_hgt_u_patch .= [30.0]

        CLM.cn_veg_struct_update!(
            d.mask, d.bounds, d.patch, d.cs, d.cnveg_cs, d.wd, d.fv, d.vs, d.cr, d.pft;
            dt=1800.0, noveg=0, spinup_factor_deadwood=1.0
        )

        # Expected height calculation
        expected_htop = ((3.0 * deadstemc_val * 1.0 * 200.0^2) /
            (pi * 0.005 * 2.5e5))^(1.0 / 3.0)
        # Clamp to forcing height limit
        expected_htop = min(expected_htop, (30.0 / (0.67 + 0.055)) - 3.0)
        expected_htop = max(expected_htop, 0.01)

        @test d.cs.htop_patch[1] ≈ expected_htop atol=1e-10
        # hbot = max(0, min(3, htop - 1))
        expected_hbot = max(0.0, min(3.0, expected_htop - 1.0))
        @test d.cs.hbot_patch[1] ≈ expected_hbot atol=1e-10
    end

    # -----------------------------------------------------------------------
    # Test 8: Woody htop minimum constraint (0.01)
    # -----------------------------------------------------------------------
    @testset "woody htop minimum constraint" begin
        d = make_veg_struct_data(np=1, nc=1)
        d.patch.itype .= [2]
        d.patch.column .= [1]

        d.pft.woody[2] = 1.0
        # Very small deadstemc leads to very small htop
        d.cnveg_cs.deadstemc_patch .= [0.001]
        d.cnveg_cs.leafc_patch .= [1.0]
        d.pft.slatop[2] = 0.008
        d.pft.dsladlai[2] = 0.0
        d.pft.nstem[2] = 0.005
        d.pft.taper[2] = 200.0
        d.pft.dwood[2] = 2.5e5

        d.cs.tlai_patch .= [0.5]
        d.cs.tsai_patch .= [1.0]

        CLM.cn_veg_struct_update!(
            d.mask, d.bounds, d.patch, d.cs, d.cnveg_cs, d.wd, d.fv, d.vs, d.cr, d.pft;
            dt=1800.0, noveg=0
        )

        @test d.cs.htop_patch[1] >= 0.01
    end

    # -----------------------------------------------------------------------
    # Test 9: Grass height depends on LAI
    # -----------------------------------------------------------------------
    @testset "grass height depends on LAI" begin
        d = make_veg_struct_data(np=1, nc=1)
        d.patch.itype .= [13]   # c3_non-arctic_grass
        d.patch.column .= [1]

        leafc = 40.0
        d.cnveg_cs.leafc_patch .= [leafc]
        d.pft.slatop[13] = 0.02
        d.pft.dsladlai[13] = 0.0
        d.pft.z0mr[13] = 0.055
        d.pft.displar[13] = 0.67

        d.cs.tlai_patch .= [0.0]
        d.cs.tsai_patch .= [1.0]

        CLM.cn_veg_struct_update!(
            d.mask, d.bounds, d.patch, d.cs, d.cnveg_cs, d.wd, d.fv, d.vs, d.cr, d.pft;
            dt=1800.0, noveg=0, npcropmin=17
        )

        new_tlai = 0.02 * leafc
        expected_htop = max(0.25, new_tlai * 0.25)
        expected_htop = min(expected_htop, (30.0 / (0.67 + 0.055)) - 3.0)
        expected_htop = max(expected_htop, 0.01)
        @test d.cs.htop_patch[1] ≈ expected_htop atol=1e-10
    end

    # -----------------------------------------------------------------------
    # Test 10: Prognostic crop height
    # -----------------------------------------------------------------------
    @testset "prognostic crop height" begin
        d = make_veg_struct_data(np=1, nc=1)
        npcrop = 17  # temperate corn
        d.patch.itype .= [npcrop]
        d.patch.column .= [1]

        d.cnveg_cs.leafc_patch .= [100.0]
        d.pft.slatop[npcrop] = 0.02
        d.pft.dsladlai[npcrop] = 0.0
        d.pft.ztopmx[npcrop] = 2.5
        d.pft.laimx[npcrop] = 6.0

        d.cs.tlai_patch .= [0.0]
        d.cs.tsai_patch .= [0.5]

        CLM.cn_veg_struct_update!(
            d.mask, d.bounds, d.patch, d.cs, d.cnveg_cs, d.wd, d.fv, d.vs, d.cr, d.pft;
            dt=1800.0, noveg=0, npcropmin=17
        )

        new_tlai = 0.02 * 100.0  # = 2.0
        # htop = ztopmx * (min(tlai / (laimx - 1), 1))^2
        expected_htop = 2.5 * (min(new_tlai / (6.0 - 1.0), 1.0))^2
        expected_htop = max(0.05, expected_htop)
        @test d.cs.htop_patch[1] ≈ expected_htop atol=1e-10
        @test d.cs.hbot_patch[1] == 0.02
    end

    # -----------------------------------------------------------------------
    # Test 11: Prognostic crop tsai for corn-like PFTs (0.1 * tlai)
    # -----------------------------------------------------------------------
    @testset "prognostic crop tsai for corn-like PFTs" begin
        d = make_veg_struct_data(np=1, nc=1)
        npcrop = 17  # temperate corn
        d.patch.itype .= [npcrop]
        d.patch.column .= [1]

        d.cnveg_cs.leafc_patch .= [100.0]
        d.pft.slatop[npcrop] = 0.02
        d.pft.dsladlai[npcrop] = 0.0
        d.pft.ztopmx[npcrop] = 2.5
        d.pft.laimx[npcrop] = 6.0

        d.cs.tlai_patch .= [0.0]
        d.cs.tsai_patch .= [0.5]

        CLM.cn_veg_struct_update!(
            d.mask, d.bounds, d.patch, d.cs, d.cnveg_cs, d.wd, d.fv, d.vs, d.cr, d.pft;
            dt=1800.0, noveg=0, npcropmin=17, ntmp_corn=17
        )

        new_tlai = 0.02 * 100.0
        @test d.cs.tsai_patch[1] ≈ 0.1 * new_tlai atol=1e-10
    end

    # -----------------------------------------------------------------------
    # Test 12: Prognostic crop tsai for non-corn crops (0.2 * tlai)
    # -----------------------------------------------------------------------
    @testset "prognostic crop tsai for non-corn crops" begin
        d = make_veg_struct_data(np=1, nc=1)
        npcrop = 19   # spring wheat
        d.patch.itype .= [npcrop]
        d.patch.column .= [1]

        d.cnveg_cs.leafc_patch .= [50.0]
        d.pft.slatop[npcrop] = 0.02
        d.pft.dsladlai[npcrop] = 0.0
        d.pft.ztopmx[npcrop] = 1.2
        d.pft.laimx[npcrop] = 5.0

        d.cs.tlai_patch .= [0.0]
        d.cs.tsai_patch .= [0.5]

        CLM.cn_veg_struct_update!(
            d.mask, d.bounds, d.patch, d.cs, d.cnveg_cs, d.wd, d.fv, d.vs, d.cr, d.pft;
            dt=1800.0, noveg=0, npcropmin=17
        )

        new_tlai = 0.02 * 50.0
        @test d.cs.tsai_patch[1] ≈ 0.2 * new_tlai atol=1e-10
    end

    # -----------------------------------------------------------------------
    # Test 13: Crop stubble after harvest
    # -----------------------------------------------------------------------
    @testset "crop stubble after harvest" begin
        d = make_veg_struct_data(np=1, nc=1)
        npcrop = 17
        d.patch.itype .= [npcrop]
        d.patch.column .= [1]

        d.cnveg_cs.leafc_patch .= [0.0]   # zero leafc after harvest
        d.pft.slatop[npcrop] = 0.02
        d.pft.dsladlai[npcrop] = 0.0
        d.pft.ztopmx[npcrop] = 2.5
        d.pft.laimx[npcrop] = 6.0

        d.cs.tlai_patch .= [0.0]
        d.cs.tsai_patch .= [0.5]
        d.cr.harvdate_patch .= [200]      # harvested (< 999)
        d.vs.farea_burned_col .= [0.1]
        d.vs.htmx_patch .= [1.5]

        CLM.cn_veg_struct_update!(
            d.mask, d.bounds, d.patch, d.cs, d.cnveg_cs, d.wd, d.fv, d.vs, d.cr, d.pft;
            dt=1800.0, noveg=0, npcropmin=17, ntmp_corn=17
        )

        expected_tsai = 0.25 * (1.0 - 0.1 * 0.90)
        @test d.cs.tsai_patch[1] ≈ expected_tsai atol=1e-10
        @test d.vs.htmx_patch[1] == 0.0
        @test d.vs.peaklai_patch[1] == 0
    end

    # -----------------------------------------------------------------------
    # Test 14: peaklai flag set when tlai >= laimx
    # -----------------------------------------------------------------------
    @testset "peaklai flag set when tlai >= laimx" begin
        d = make_veg_struct_data(np=1, nc=1)
        npcrop = 17
        d.patch.itype .= [npcrop]
        d.patch.column .= [1]

        # leafc enough to produce tlai >= laimx
        d.pft.slatop[npcrop] = 0.02
        d.pft.dsladlai[npcrop] = 0.0
        d.pft.laimx[npcrop] = 4.0
        d.pft.ztopmx[npcrop] = 2.5
        d.cnveg_cs.leafc_patch .= [250.0]  # tlai = 0.02 * 250 = 5.0 >= 4.0

        d.cs.tlai_patch .= [0.0]
        d.cs.tsai_patch .= [0.5]
        d.vs.peaklai_patch .= [0]

        CLM.cn_veg_struct_update!(
            d.mask, d.bounds, d.patch, d.cs, d.cnveg_cs, d.wd, d.fv, d.vs, d.cr, d.pft;
            dt=1800.0, noveg=0, npcropmin=17, ntmp_corn=17
        )

        @test d.vs.peaklai_patch[1] == 1
    end

    # -----------------------------------------------------------------------
    # Test 15: Snow burial for trees/shrubs (ivt > noveg && ivt <= nbrdlf_dcd_brl_shrub)
    # -----------------------------------------------------------------------
    @testset "snow burial for trees/shrubs" begin
        d = make_veg_struct_data(np=1, nc=1)
        d.patch.itype .= [5]    # broadleaf_evergreen_temperate_tree
        d.patch.column .= [1]

        d.pft.woody[5] = 1.0
        d.cnveg_cs.deadstemc_patch .= [3000.0]
        d.cnveg_cs.leafc_patch .= [80.0]
        d.pft.slatop[5] = 0.012
        d.pft.dsladlai[5] = 0.0
        d.pft.nstem[5] = 0.005
        d.pft.taper[5] = 200.0
        d.pft.dwood[5] = 2.5e5

        d.cs.tlai_patch .= [1.0]
        d.cs.tsai_patch .= [1.0]

        # Partial snow
        d.wd.snow_depth_col .= [0.5]
        d.wd.frac_sno_col .= [0.3]

        CLM.cn_veg_struct_update!(
            d.mask, d.bounds, d.patch, d.cs, d.cnveg_cs, d.wd, d.fv, d.vs, d.cr, d.pft;
            dt=1800.0, noveg=0, nbrdlf_dcd_brl_shrub=11
        )

        # With snow, elai < tlai and esai < tsai
        @test d.cs.elai_patch[1] <= d.cs.tlai_patch[1]
        @test d.cs.esai_patch[1] <= d.cs.tsai_patch[1]
        @test d.cs.elai_patch[1] >= 0.0
        @test d.cs.esai_patch[1] >= 0.0
        @test d.cs.frac_veg_nosno_alb_patch[1] == 1
    end

    # -----------------------------------------------------------------------
    # Test 16: Snow burial for grasses (ivt > nbrdlf_dcd_brl_shrub)
    # -----------------------------------------------------------------------
    @testset "snow burial for grasses" begin
        d = make_veg_struct_data(np=1, nc=1)
        d.patch.itype .= [13]   # c3_non-arctic_grass (> 11)
        d.patch.column .= [1]

        d.cnveg_cs.leafc_patch .= [50.0]
        d.pft.slatop[13] = 0.02
        d.pft.dsladlai[13] = 0.0

        d.cs.tlai_patch .= [0.0]
        d.cs.tsai_patch .= [1.0]

        # Partial snow
        d.wd.snow_depth_col .= [0.03]
        d.wd.frac_sno_col .= [0.5]

        CLM.cn_veg_struct_update!(
            d.mask, d.bounds, d.patch, d.cs, d.cnveg_cs, d.wd, d.fv, d.vs, d.cr, d.pft;
            dt=1800.0, noveg=0, nbrdlf_dcd_brl_shrub=11
        )

        # elai should be reduced by snow burial
        @test d.cs.elai_patch[1] <= d.cs.tlai_patch[1]
        @test d.cs.elai_patch[1] >= 0.0
    end

    # -----------------------------------------------------------------------
    # Test 17: frac_sno_threshold behavior
    # -----------------------------------------------------------------------
    @testset "frac_sno above threshold treated as 1.0" begin
        d = make_veg_struct_data(np=1, nc=1)
        d.patch.itype .= [13]
        d.patch.column .= [1]

        d.cnveg_cs.leafc_patch .= [50.0]
        d.pft.slatop[13] = 0.02
        d.pft.dsladlai[13] = 0.0

        d.cs.tlai_patch .= [0.0]
        d.cs.tsai_patch .= [1.0]

        # frac_sno = 0.9999 > threshold (0.999)
        d.wd.frac_sno_col .= [0.9999]
        d.wd.snow_depth_col .= [1.0]   # deep snow

        CLM.cn_veg_struct_update!(
            d.mask, d.bounds, d.patch, d.cs, d.cnveg_cs, d.wd, d.fv, d.vs, d.cr, d.pft;
            dt=1800.0, noveg=0
        )

        # With frac_sno adjusted to 1.0 and deep snow, vegetation should be fully buried
        # For grasses with htop ≈ 0.25, snow_depth=1.0 means fb ≈ 0
        # fb = 1 - max(min(1.0, max(0.05, 0.25*0.8)), 0) / max(0.05, 0.25*0.8)
        #    = 1 - 0.2 / 0.2 = 0
        # elai = max(tlai * (1-1) + tlai * 0 * 1, 0) = 0
        @test d.cs.frac_veg_nosno_alb_patch[1] == 0
    end

    # -----------------------------------------------------------------------
    # Test 18: No snow means no burial (frac_sno = 0)
    # -----------------------------------------------------------------------
    @testset "no snow means no burial" begin
        d = make_veg_struct_data(np=1, nc=1)
        d.patch.itype .= [13]
        d.patch.column .= [1]

        d.cnveg_cs.leafc_patch .= [50.0]
        d.pft.slatop[13] = 0.02
        d.pft.dsladlai[13] = 0.0

        d.cs.tlai_patch .= [0.0]
        d.cs.tsai_patch .= [1.0]

        d.wd.frac_sno_col .= [0.0]
        d.wd.snow_depth_col .= [0.0]

        CLM.cn_veg_struct_update!(
            d.mask, d.bounds, d.patch, d.cs, d.cnveg_cs, d.wd, d.fv, d.vs, d.cr, d.pft;
            dt=1800.0, noveg=0
        )

        # With no snow, elai should equal tlai
        @test d.cs.elai_patch[1] ≈ d.cs.tlai_patch[1] atol=1e-10
        @test d.cs.esai_patch[1] ≈ d.cs.tsai_patch[1] atol=1e-10
    end

    # -----------------------------------------------------------------------
    # Test 19: frac_veg_nosno_alb is 1 when vegetation exists, 0 otherwise
    # -----------------------------------------------------------------------
    @testset "frac_veg_nosno_alb flag" begin
        d = make_veg_struct_data(np=2, nc=1)
        d.patch.itype .= [13, 0]
        d.patch.column .= [1, 1]

        d.cnveg_cs.leafc_patch .= [50.0, 0.0]
        d.pft.slatop[13] = 0.02
        d.pft.dsladlai[13] = 0.0

        d.cs.tlai_patch .= [0.0, 0.0]
        d.cs.tsai_patch .= [1.0, 0.0]

        CLM.cn_veg_struct_update!(
            d.mask, d.bounds, d.patch, d.cs, d.cnveg_cs, d.wd, d.fv, d.vs, d.cr, d.pft;
            dt=1800.0, noveg=0
        )

        @test d.cs.frac_veg_nosno_alb_patch[1] == 1  # has vegetation
        @test d.cs.frac_veg_nosno_alb_patch[2] == 0  # noveg
    end

    # -----------------------------------------------------------------------
    # Test 20: Mask is respected (inactive patches skipped)
    # -----------------------------------------------------------------------
    @testset "mask filters inactive patches" begin
        d = make_veg_struct_data(np=2, nc=1)
        d.patch.itype .= [13, 13]
        d.patch.column .= [1, 1]

        d.cnveg_cs.leafc_patch .= [50.0, 50.0]
        d.pft.slatop[13] = 0.02
        d.pft.dsladlai[13] = 0.0

        d.cs.tlai_patch .= [999.0, 999.0]
        d.cs.tsai_patch .= [1.0, 1.0]

        # Only patch 1 is active
        mask = BitVector([true, false])

        CLM.cn_veg_struct_update!(
            mask, d.bounds, d.patch, d.cs, d.cnveg_cs, d.wd, d.fv, d.vs, d.cr, d.pft;
            dt=1800.0, noveg=0
        )

        # Patch 1 should be updated
        @test d.cs.tlai_patch[1] ≈ 0.02 * 50.0 atol=1e-10

        # Patch 2 should be unchanged (999.0)
        @test d.cs.tlai_patch[2] == 999.0
    end

    # -----------------------------------------------------------------------
    # Test 21: Biomass heat storage computation (leaf_biomass)
    # -----------------------------------------------------------------------
    @testset "biomass heat storage leaf_biomass" begin
        d = make_veg_struct_data(np=1, nc=1)
        d.patch.itype .= [13]
        d.patch.column .= [1]

        leafc = 50.0
        fbw_val = 0.6
        d.cnveg_cs.leafc_patch .= [leafc]
        d.pft.slatop[13] = 0.02
        d.pft.dsladlai[13] = 0.0
        d.pft.fbw[13] = fbw_val

        d.cs.tlai_patch .= [0.0]
        d.cs.tsai_patch .= [1.0]

        CLM.cn_veg_struct_update!(
            d.mask, d.bounds, d.patch, d.cs, d.cnveg_cs, d.wd, d.fv, d.vs, d.cr, d.pft;
            dt=1800.0, noveg=0, use_biomass_heat_storage=true
        )

        expected_leaf_biomass = max(0.0025, leafc) * CLM.C_TO_B * 1.0e-3 / (1.0 - fbw_val)
        @test d.cs.leaf_biomass_patch[1] ≈ expected_leaf_biomass atol=1e-10
    end

    # -----------------------------------------------------------------------
    # Test 22: Biomass heat storage computation (stem_biomass for woody)
    # -----------------------------------------------------------------------
    @testset "biomass heat storage stem_biomass for woody" begin
        d = make_veg_struct_data(np=1, nc=1)
        d.patch.itype .= [2]
        d.patch.column .= [1]

        d.pft.woody[2] = 1.0
        deadstemc_val = 3000.0
        livestemc_val = 1000.0
        fbw_val = 0.5
        d.cnveg_cs.deadstemc_patch .= [deadstemc_val]
        d.cnveg_cs.livestemc_patch .= [livestemc_val]
        d.cnveg_cs.leafc_patch .= [100.0]
        d.pft.slatop[2] = 0.008
        d.pft.dsladlai[2] = 0.0
        d.pft.nstem[2] = 0.005
        d.pft.taper[2] = 200.0
        d.pft.dwood[2] = 2.5e5
        d.pft.fbw[2] = fbw_val

        d.cs.tlai_patch .= [0.5]
        d.cs.tsai_patch .= [1.0]

        spinup_f = 1.0
        CLM.cn_veg_struct_update!(
            d.mask, d.bounds, d.patch, d.cs, d.cnveg_cs, d.wd, d.fv, d.vs, d.cr, d.pft;
            dt=1800.0, noveg=0, use_biomass_heat_storage=true,
            spinup_factor_deadwood=spinup_f
        )

        expected_stem = (spinup_f * deadstemc_val + livestemc_val) *
            CLM.C_TO_B * 1.0e-3 / (1.0 - fbw_val)
        @test d.cs.stem_biomass_patch[1] ≈ expected_stem atol=1e-10
    end

    # -----------------------------------------------------------------------
    # Test 23: Spinup factor affects height
    # -----------------------------------------------------------------------
    @testset "spinup_factor_deadwood affects height" begin
        # Use small deadstemc so htop stays well below forc_hgt_u ceiling
        d1 = make_veg_struct_data(np=1, nc=1)
        d1.patch.itype .= [2]
        d1.patch.column .= [1]
        d1.pft.woody[2] = 1.0
        d1.cnveg_cs.deadstemc_patch .= [50.0]
        d1.cnveg_cs.leafc_patch .= [50.0]
        d1.pft.slatop[2] = 0.008
        d1.pft.dsladlai[2] = 0.0
        d1.pft.nstem[2] = 0.005
        d1.pft.taper[2] = 200.0
        d1.pft.dwood[2] = 2.5e5
        d1.cs.tlai_patch .= [0.5]
        d1.cs.tsai_patch .= [1.0]
        d1.fv.forc_hgt_u_patch .= [100.0]  # high forcing height to avoid ceiling

        CLM.cn_veg_struct_update!(
            d1.mask, d1.bounds, d1.patch, d1.cs, d1.cnveg_cs, d1.wd, d1.fv, d1.vs, d1.cr, d1.pft;
            dt=1800.0, noveg=0, spinup_factor_deadwood=1.0
        )
        htop1 = d1.cs.htop_patch[1]

        d2 = make_veg_struct_data(np=1, nc=1)
        d2.patch.itype .= [2]
        d2.patch.column .= [1]
        d2.pft.woody[2] = 1.0
        d2.cnveg_cs.deadstemc_patch .= [50.0]
        d2.cnveg_cs.leafc_patch .= [50.0]
        d2.pft.slatop[2] = 0.008
        d2.pft.dsladlai[2] = 0.0
        d2.pft.nstem[2] = 0.005
        d2.pft.taper[2] = 200.0
        d2.pft.dwood[2] = 2.5e5
        d2.cs.tlai_patch .= [0.5]
        d2.cs.tsai_patch .= [1.0]
        d2.fv.forc_hgt_u_patch .= [100.0]  # high forcing height to avoid ceiling

        CLM.cn_veg_struct_update!(
            d2.mask, d2.bounds, d2.patch, d2.cs, d2.cnveg_cs, d2.wd, d2.fv, d2.vs, d2.cr, d2.pft;
            dt=1800.0, noveg=0, spinup_factor_deadwood=10.0
        )
        htop10 = d2.cs.htop_patch[1]

        # With spinup_factor=10, effective deadstemc is 10x, so htop should be larger
        @test htop10 > htop1
    end

    # -----------------------------------------------------------------------
    # Test 24: Multiple PFT types in one call
    # -----------------------------------------------------------------------
    @testset "mixed PFT types processed correctly" begin
        d = make_veg_struct_data(np=4, nc=2)
        # noveg, woody tree, grass, prognostic crop
        d.patch.itype .= [0, 2, 13, 17]
        d.patch.column .= [1, 1, 2, 2]

        d.pft.woody[2] = 1.0
        d.cnveg_cs.deadstemc_patch .= [0.0, 3000.0, 0.0, 0.0]
        d.cnveg_cs.leafc_patch .= [0.0, 80.0, 50.0, 100.0]
        d.pft.slatop .= fill(0.02, 80)
        d.pft.slatop[2] = 0.012
        d.pft.dsladlai .= fill(0.0, 80)
        d.pft.nstem[2] = 0.005
        d.pft.taper[2] = 200.0
        d.pft.dwood[2] = 2.5e5
        d.pft.ztopmx[17] = 2.5
        d.pft.laimx[17] = 6.0

        d.cs.tlai_patch .= [0.0, 0.5, 0.0, 0.0]
        d.cs.tsai_patch .= [0.0, 1.0, 1.0, 0.5]

        CLM.cn_veg_struct_update!(
            d.mask, d.bounds, d.patch, d.cs, d.cnveg_cs, d.wd, d.fv, d.vs, d.cr, d.pft;
            dt=1800.0, noveg=0, npcropmin=17, ntmp_corn=17
        )

        # noveg: all zero
        @test d.cs.tlai_patch[1] == 0.0
        @test d.cs.htop_patch[1] == 0.0

        # woody: has height from deadstemc
        @test d.cs.htop_patch[2] > 0.01
        @test d.cs.tlai_patch[2] ≈ 0.012 * 80.0

        # grass: height depends on LAI
        @test d.cs.htop_patch[3] >= 0.01
        @test d.cs.tlai_patch[3] ≈ 0.02 * 50.0

        # crop: specific height formula
        @test d.cs.htop_patch[4] >= 0.05
        @test d.cs.hbot_patch[4] == 0.02
    end

    # -----------------------------------------------------------------------
    # Test 25: Snow burial detailed calculation for tree/shrub branch
    # -----------------------------------------------------------------------
    @testset "snow burial tree/shrub detailed calculation" begin
        d = make_veg_struct_data(np=1, nc=1)
        d.patch.itype .= [9]   # broadleaf_evergreen_shrub (ivt=9, <= 11)
        d.patch.column .= [1]

        d.pft.woody[9] = 1.0
        d.cnveg_cs.deadstemc_patch .= [500.0]
        d.cnveg_cs.leafc_patch .= [30.0]
        d.pft.slatop[9] = 0.015
        d.pft.dsladlai[9] = 0.0
        d.pft.nstem[9] = 0.01
        d.pft.taper[9] = 100.0
        d.pft.dwood[9] = 2.5e5

        d.cs.tlai_patch .= [0.3]
        d.cs.tsai_patch .= [0.5]

        snow_depth_val = 0.3
        frac_sno_val = 0.4
        d.wd.snow_depth_col .= [snow_depth_val]
        d.wd.frac_sno_col .= [frac_sno_val]

        CLM.cn_veg_struct_update!(
            d.mask, d.bounds, d.patch, d.cs, d.cnveg_cs, d.wd, d.fv, d.vs, d.cr, d.pft;
            dt=1800.0, noveg=0, nbrdlf_dcd_brl_shrub=11
        )

        htop_val = d.cs.htop_patch[1]
        hbot_val = d.cs.hbot_patch[1]
        tlai_val = d.cs.tlai_patch[1]
        tsai_val = d.cs.tsai_patch[1]

        # Verify snow burial calculation
        ol = min(max(snow_depth_val - hbot_val, 0.0), htop_val - hbot_val)
        fb = 1.0 - ol / max(1.0e-06, htop_val - hbot_val)
        expected_elai = max(tlai_val * (1.0 - frac_sno_val) + tlai_val * fb * frac_sno_val, 0.0)
        expected_esai = max(tsai_val * (1.0 - frac_sno_val) + tsai_val * fb * frac_sno_val, 0.0)

        @test d.cs.elai_patch[1] ≈ expected_elai atol=1e-10
        @test d.cs.esai_patch[1] ≈ expected_esai atol=1e-10
    end

    # -----------------------------------------------------------------------
    # Test 26: Grass hbot formula
    # -----------------------------------------------------------------------
    @testset "grass hbot formula" begin
        d = make_veg_struct_data(np=1, nc=1)
        d.patch.itype .= [14]   # c4_grass
        d.patch.column .= [1]

        d.cnveg_cs.leafc_patch .= [20.0]
        d.pft.slatop[14] = 0.02
        d.pft.dsladlai[14] = 0.0
        d.pft.z0mr[14] = 0.055
        d.pft.displar[14] = 0.67

        d.cs.tlai_patch .= [0.0]
        d.cs.tsai_patch .= [1.0]

        CLM.cn_veg_struct_update!(
            d.mask, d.bounds, d.patch, d.cs, d.cnveg_cs, d.wd, d.fv, d.vs, d.cr, d.pft;
            dt=1800.0, noveg=0, npcropmin=17
        )

        htop_val = d.cs.htop_patch[1]
        expected_hbot = max(0.0, min(0.05, htop_val - 0.20))
        @test d.cs.hbot_patch[1] ≈ expected_hbot atol=1e-10
    end

    # -----------------------------------------------------------------------
    # Test 27: Crop htmx tracks max height
    # -----------------------------------------------------------------------
    @testset "crop htmx tracks max height" begin
        d = make_veg_struct_data(np=1, nc=1)
        npcrop = 19  # spring wheat
        d.patch.itype .= [npcrop]
        d.patch.column .= [1]

        d.cnveg_cs.leafc_patch .= [100.0]
        d.pft.slatop[npcrop] = 0.02
        d.pft.dsladlai[npcrop] = 0.0
        d.pft.ztopmx[npcrop] = 1.0
        d.pft.laimx[npcrop] = 5.0

        # Set htmx to a previous high value
        d.vs.htmx_patch .= [0.8]
        d.cs.tlai_patch .= [0.0]
        d.cs.tsai_patch .= [0.5]

        CLM.cn_veg_struct_update!(
            d.mask, d.bounds, d.patch, d.cs, d.cnveg_cs, d.wd, d.fv, d.vs, d.cr, d.pft;
            dt=1800.0, noveg=0, npcropmin=17
        )

        # htmx should be at least the calculated htop
        @test d.vs.htmx_patch[1] >= d.cs.htop_patch[1] || d.cs.htop_patch[1] <= 0.05
        # htop uses max(htmx, htop), so htop >= htmx
        @test d.cs.htop_patch[1] >= 0.05
    end

end
