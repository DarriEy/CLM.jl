@testset "CNDV (Dynamic Vegetation)" begin

    # -----------------------------------------------------------------------
    # Helper: create minimal PftconType for CNDV tests
    # -----------------------------------------------------------------------
    function make_cndv_pftcon(npft::Int)
        pft = CLM.PftconType()
        CLM.pftcon_allocate!(pft)
        n = length(pft.woody)

        # Set up PFT properties for first npft+1 entries (0-based noveg + PFTs)
        # PFT 0 = bare ground (noveg)
        # PFT 1 = needleleaf evr tmp tree
        # PFT 2 = needleleaf evr brl tree
        # PFT 3 = needleleaf dcd brl tree (deciduous)
        # PFT 12 = C3 arctic grass
        # PFT 13 = C3 non-arctic grass

        # Julia 1-based index = Fortran PFT + 1
        # Trees: indices 2,3,4 (PFTs 1,2,3)
        for i in [2, 3, 4]
            pft.woody[i] = 1.0
            pft.is_tree[i] = true
            pft.is_shrub[i] = false
            pft.dwood[i] = 2.5e5
            pft.slatop[i] = 0.012
            pft.dsladlai[i] = 0.0
        end

        # Grasses: indices 13,14 (PFTs 12,13)
        for i in [13, 14]
            pft.woody[i] = 0.0
            pft.is_tree[i] = false
            pft.is_shrub[i] = false
            pft.dwood[i] = 0.0
            pft.slatop[i] = 0.030
            pft.dsladlai[i] = 0.0
        end

        # DGVS ecophysiological constants
        pft.pftpar20 .= 15.0      # crownarea_max
        pft.pftpar28 .= -35.0     # tcmin
        pft.pftpar29 .= 28.0      # tcmax
        pft.pftpar30 .= 350.0     # gddmin
        pft.pftpar31 .= 1000.0    # twmax (>999 = no heat limit for most PFTs)

        return pft
    end

    # -----------------------------------------------------------------------
    # Helper: create minimal DGVS ecophys constants
    # -----------------------------------------------------------------------
    function make_cndv_eco(pft::CLM.PftconType)
        eco = CLM.DGVEcophysCon()
        CLM.dgv_ecophyscon_init!(eco, pft)
        return eco
    end

    # -----------------------------------------------------------------------
    # Helper: create a minimal CNDV test configuration
    # Patches: 1=bare, 2=tree(PFT1), 3=tree(PFT2), 4=grass(PFT12)
    # -----------------------------------------------------------------------
    function make_cndv_test_data()
        np = 4
        nc = 1
        ng = 1

        pft = make_cndv_pftcon(15)
        eco = make_cndv_eco(pft)

        dgvs = CLM.DGVSData()
        CLM.dgvs_init!(dgvs, np)

        # Patch data
        patch = CLM.PatchData()
        CLM.patch_init!(patch, np)
        patch.itype .= [CLM.noveg, 1, 2, 12]  # 0-based Fortran PFT indices
        patch.gridcell .= 1
        patch.column .= 1
        patch.landunit .= 1
        patch.wtcol .= [0.1, 0.3, 0.3, 0.3]

        # Landunit data
        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, 1)
        lun.itype[1] = CLM.ISTSOIL
        lun.wtgcell[1] = 1.0

        return (; dgvs, eco, pft, patch, lun, np, nc, ng)
    end

    # -----------------------------------------------------------------------
    # 1. DGVSData construction and init
    # -----------------------------------------------------------------------
    @testset "DGVSData construction and init" begin
        d = CLM.DGVSData()
        @test length(d.agdd_patch) == 0
        @test length(d.present_patch) == 0

        CLM.dgvs_init!(d, 5)
        @test length(d.agdd_patch) == 5
        @test length(d.present_patch) == 5
        @test length(d.nind_patch) == 5
        @test length(d.fpcgrid_patch) == 5
        @test all(iszero, d.nind_patch)
    end

    @testset "DGVSData cold-start init" begin
        d = CLM.DGVSData()
        CLM.dgvs_init_cold!(d, 3)
        @test length(d.present_patch) == 3
        @test all(x -> x == false, d.present_patch)
        @test all(iszero, d.crownarea_patch)
        @test all(iszero, d.nind_patch)
        @test d.tmomin20_patch[1] ≈ CLM.TFRZ - 5.0
    end

    # -----------------------------------------------------------------------
    # 2. DGVEcophysCon initialization
    # -----------------------------------------------------------------------
    @testset "DGVEcophysCon init from pftcon" begin
        pft = make_cndv_pftcon(15)
        eco = CLM.DGVEcophysCon()
        CLM.dgv_ecophyscon_init!(eco, pft)

        @test length(eco.crownarea_max) == length(pft.pftpar20)
        @test length(eco.tcmin) == length(pft.pftpar28)
        @test eco.crownarea_max[1] ≈ 15.0
        @test eco.tcmin[1] ≈ -35.0
        @test eco.tcmax[1] ≈ 28.0
        @test eco.gddmin[1] ≈ 350.0
        @test eco.twmax[1] ≈ 1000.0
        @test eco.reinickerp[1] ≈ CLM.REINICKERP
        @test eco.allom1[1] ≈ CLM.ALLOM1  # non-shrub
    end

    # -----------------------------------------------------------------------
    # 3. dyn_cndv_init! — fpcgrid = wtcol
    # -----------------------------------------------------------------------
    @testset "dyn_cndv_init!" begin
        d = make_cndv_test_data()
        CLM.dyn_cndv_init!(d.dgvs, d.patch, 1:d.np)

        for p in 1:d.np
            @test d.dgvs.fpcgrid_patch[p] ≈ d.patch.wtcol[p]
            @test d.dgvs.fpcgridold_patch[p] ≈ d.patch.wtcol[p]
        end
    end

    # -----------------------------------------------------------------------
    # 4. dyn_cndv_interp! — weight interpolation
    # -----------------------------------------------------------------------
    @testset "dyn_cndv_interp!" begin
        d = make_cndv_test_data()
        # Set fpcgrid and fpcgridold to known values
        d.dgvs.fpcgrid_patch .= [0.05, 0.35, 0.30, 0.30]
        d.dgvs.fpcgridold_patch .= [0.10, 0.30, 0.30, 0.30]

        # wt1=0.5 means halfway between old and new
        CLM.dyn_cndv_interp!(d.dgvs, d.patch, d.lun, 1:d.np;
                              wt1=0.5, is_beg_curr_year=false)

        # wtcol = fpcgrid + wt1*(fpcgridold - fpcgrid)
        @test d.patch.wtcol[1] ≈ 0.05 + 0.5 * (0.10 - 0.05)  # 0.075
        @test d.patch.wtcol[2] ≈ 0.35 + 0.5 * (0.30 - 0.35)  # 0.325

        # fpcgridold should NOT be updated (is_beg_curr_year=false)
        @test d.dgvs.fpcgridold_patch[1] ≈ 0.10

        # Now test with is_beg_curr_year=true
        CLM.dyn_cndv_interp!(d.dgvs, d.patch, d.lun, 1:d.np;
                              wt1=0.0, is_beg_curr_year=true)

        # With wt1=0, wtcol = fpcgrid
        @test d.patch.wtcol[1] ≈ 0.05
        # fpcgridold should be updated to fpcgrid
        @test d.dgvs.fpcgridold_patch[1] ≈ 0.05
    end

    # -----------------------------------------------------------------------
    # 5. cndv_update_acc_vars! — GDD accumulation
    # -----------------------------------------------------------------------
    @testset "cndv_update_acc_vars! GDD accumulation" begin
        d = make_cndv_test_data()
        dtime = 1800.0  # 30 min timestep

        # Set initial values
        d.dgvs.agddtw_patch .= 10.0
        d.dgvs.agdd_patch .= 20.0

        # Temperature above threshold: t_ref2m > TFRZ + 5
        t_ref2m = fill(CLM.TFRZ + 15.0, d.np)
        # t_a10 above twmax of boreal deciduous tree
        # twmax for ndllf_dcd_brl_tree (PFT 3, Julia index 4)
        twmax_brl = d.eco.twmax[CLM.ndllf_dcd_brl_tree + 1]
        t_a10 = fill(CLM.TFRZ + twmax_brl + 5.0, d.np)

        CLM.cndv_update_acc_vars!(d.dgvs, d.eco, t_a10, t_ref2m,
                                   1:d.np, dtime;
                                   month=6, day=15, secs=1800)

        # AGDD should increase: max(0, (15 - 5)) * dtime/SECSPDAY
        agdd_inc = (15.0 - 5.0) * dtime / CLM.SECSPDAY
        @test d.dgvs.agdd_patch[1] ≈ 20.0 + agdd_inc

        # AGDDTW should increase: t_a10 is above twmax
        @test d.dgvs.agddtw_patch[1] > 10.0
    end

    @testset "cndv_update_acc_vars! annual reset" begin
        d = make_cndv_test_data()
        dtime = 1800.0

        d.dgvs.agddtw_patch .= 100.0
        d.dgvs.agdd_patch .= 200.0

        t_ref2m = fill(CLM.TFRZ + 15.0, d.np)
        t_a10 = fill(CLM.TFRZ + 10.0, d.np)

        # Trigger annual reset: month=1, day=1, secs=dtime
        CLM.cndv_update_acc_vars!(d.dgvs, d.eco, t_a10, t_ref2m,
                                   1:d.np, dtime;
                                   month=1, day=1, secs=Int(dtime))

        # After reset, values should be only the current timestep contribution
        agdd_inc = max(0.0, (15.0 - 5.0)) * dtime / CLM.SECSPDAY
        @test d.dgvs.agdd_patch[1] ≈ agdd_inc  # reset to 0 then accumulated
    end

    @testset "cndv_update_acc_vars! cold temperature" begin
        d = make_cndv_test_data()
        dtime = 1800.0

        d.dgvs.agdd_patch .= 50.0
        d.dgvs.agddtw_patch .= 30.0

        # Temperature below thresholds
        t_ref2m = fill(CLM.TFRZ + 2.0, d.np)  # below TFRZ+5
        twmax_brl = d.eco.twmax[CLM.ndllf_dcd_brl_tree + 1]
        t_a10 = fill(CLM.TFRZ + twmax_brl - 5.0, d.np)  # below twmax

        CLM.cndv_update_acc_vars!(d.dgvs, d.eco, t_a10, t_ref2m,
                                   1:d.np, dtime;
                                   month=6, day=15, secs=1800)

        # AGDD should not increase (T < TFRZ+5)
        @test d.dgvs.agdd_patch[1] ≈ 50.0
        # AGDDTW should not increase (t_a10 < TFRZ + twmax)
        @test d.dgvs.agddtw_patch[1] ≈ 30.0
    end

    # -----------------------------------------------------------------------
    # 6. cndv_light! — light competition
    # -----------------------------------------------------------------------
    @testset "cndv_light! crown area and FPC computation" begin
        d = make_cndv_test_data()
        mask_natvegp = trues(d.np)

        # Set up present trees with known state
        d.dgvs.present_patch .= [false, true, true, true]
        d.dgvs.nind_patch .= [0.0, 0.05, 0.03, 1.0]
        d.dgvs.fpcgrid_patch .= [0.0, 0.3, 0.3, 0.3]
        d.dgvs.crownarea_patch .= [0.0, 8.0, 8.0, 1.0]

        # Carbon pools
        deadstemc = [0.0, 5000.0, 3000.0, 0.0]
        leafcmax = [0.0, 200.0, 150.0, 50.0]

        CLM.cndv_light!(d.dgvs, d.eco, deadstemc, leafcmax,
                         d.pft, d.patch, mask_natvegp, 1:d.np;
                         bounds_gridcell=1:d.ng)

        # Trees (patches 2,3) should have updated crown area
        @test d.dgvs.crownarea_patch[2] > 0.0
        @test d.dgvs.crownarea_patch[3] > 0.0

        # FPC should be computed: fpcgrid = crownarea * nind * fpc_ind
        @test d.dgvs.fpcgrid_patch[2] >= 0.0
        @test d.dgvs.fpcgrid_patch[3] >= 0.0

        # Grass (patch 4) fpc should be non-negative
        @test d.dgvs.fpcgrid_patch[4] >= 0.0

        # Sum of all fpcgrid should be reasonable (<=1)
        total_fpc = sum(d.dgvs.fpcgrid_patch)
        @test total_fpc <= 1.5  # may exceed 1 before establishment adjustment
    end

    @testset "cndv_light! tree competition when fpc_tree > 0.95" begin
        d = make_cndv_test_data()
        mask_natvegp = trues(d.np)

        # Set up trees dominating > 0.95
        d.dgvs.present_patch .= [false, true, true, true]
        d.dgvs.nind_patch .= [0.0, 0.10, 0.10, 1.0]
        d.dgvs.fpcgrid_patch .= [0.0, 0.50, 0.50, 0.05]
        d.dgvs.crownarea_patch .= [0.0, 10.0, 10.0, 1.0]

        deadstemc = [0.0, 10000.0, 10000.0, 0.0]
        leafcmax = [0.0, 500.0, 500.0, 50.0]

        CLM.cndv_light!(d.dgvs, d.eco, deadstemc, leafcmax,
                         d.pft, d.patch, mask_natvegp, 1:d.np;
                         bounds_gridcell=1:d.ng)

        # After light competition, tree fpc should be capped
        fpc_tree = d.dgvs.fpcgrid_patch[2] + d.dgvs.fpcgrid_patch[3]
        @test fpc_tree <= 1.0  # should be reduced

        # nind should be reduced for trees
        @test d.dgvs.nind_patch[2] >= 0.0
        @test d.dgvs.nind_patch[3] >= 0.0
    end

    # -----------------------------------------------------------------------
    # 7. cndv_establishment! — bioclim checks
    # -----------------------------------------------------------------------
    @testset "cndv_establishment! PFT survives warm climate" begin
        d = make_cndv_test_data()

        # Present tree with warm climate (survives)
        d.dgvs.present_patch .= [false, true, true, true]
        d.dgvs.nind_patch .= [0.0, 0.05, 0.03, 1.0]
        d.dgvs.fpcgrid_patch .= [0.0, 0.3, 0.3, 0.3]
        d.dgvs.crownarea_patch .= [0.0, 8.0, 8.0, 1.0]
        d.dgvs.pftmayexist_patch .= true

        # Set climate: warm enough to survive
        # tmomin20 >= tcmin + TFRZ for all PFTs
        d.dgvs.tmomin20_patch .= CLM.TFRZ + 0.0  # 0 C > tcmin=-35
        d.dgvs.agdd20_patch .= 500.0  # > gddmin=350

        prec365 = [200.0 / (365.0 * CLM.SECSPDAY)]  # above min precip
        annsum_npp = fill(100.0, d.np)
        annsum_litfall = fill(50.0, d.np)
        deadstemc = [0.0, 5000.0, 3000.0, 0.0]
        leafcmax = [0.0, 200.0, 150.0, 50.0]

        CLM.cndv_establishment!(d.dgvs, d.eco, prec365,
                                 annsum_npp, annsum_litfall,
                                 deadstemc, leafcmax,
                                 d.pft, d.patch, d.lun, 1:d.np;
                                 bounds_gridcell=1:d.ng)

        # Trees should survive
        @test d.dgvs.present_patch[2] == true
        @test d.dgvs.present_patch[3] == true
        # Grass should survive
        @test d.dgvs.present_patch[4] == true
    end

    @testset "cndv_establishment! PFT dies in cold climate" begin
        d = make_cndv_test_data()

        # Present tree but too cold
        d.dgvs.present_patch .= [false, true, false, false]
        d.dgvs.nind_patch .= [0.0, 0.05, 0.0, 0.0]
        d.dgvs.fpcgrid_patch .= [0.0, 0.3, 0.0, 0.0]
        d.dgvs.crownarea_patch .= [0.0, 8.0, 0.0, 0.0]
        d.dgvs.pftmayexist_patch .= true

        # Set climate: too cold (tmomin20 < tcmin + TFRZ)
        tcmin_pft2 = d.eco.tcmin[2]  # Julia index 2 = Fortran PFT 1
        d.dgvs.tmomin20_patch .= CLM.TFRZ + tcmin_pft2 - 10.0  # well below tcmin
        d.dgvs.agdd20_patch .= 500.0

        prec365 = [200.0 / (365.0 * CLM.SECSPDAY)]
        annsum_npp = fill(100.0, d.np)
        annsum_litfall = fill(50.0, d.np)
        deadstemc = [0.0, 5000.0, 0.0, 0.0]
        leafcmax = [0.0, 200.0, 0.0, 0.0]

        CLM.cndv_establishment!(d.dgvs, d.eco, prec365,
                                 annsum_npp, annsum_litfall,
                                 deadstemc, leafcmax,
                                 d.pft, d.patch, d.lun, 1:d.np;
                                 bounds_gridcell=1:d.ng)

        # Tree should be killed (too cold)
        @test d.dgvs.present_patch[2] == false
        @test d.dgvs.nind_patch[2] ≈ 0.0
        @test d.dgvs.fpcgrid_patch[2] ≈ 0.0
    end

    @testset "cndv_establishment! new PFT introduction" begin
        d = make_cndv_test_data()

        # No PFTs present initially
        d.dgvs.present_patch .= false
        d.dgvs.nind_patch .= 0.0
        d.dgvs.fpcgrid_patch .= 0.0
        d.dgvs.crownarea_patch .= 0.0
        d.dgvs.pftmayexist_patch .= true

        # Climate suitable for establishment
        d.dgvs.tmomin20_patch .= CLM.TFRZ + 0.0  # > tcmin + TFRZ
        d.dgvs.agdd20_patch .= 500.0  # > gddmin
        d.dgvs.agddtw_patch .= 0.0    # twmax > 999, so agddtw check passes

        prec365 = [200.0 / (365.0 * CLM.SECSPDAY)]  # above min precip
        annsum_npp = fill(0.0, d.np)
        annsum_litfall = fill(0.0, d.np)
        deadstemc = fill(0.1, d.np)
        leafcmax = fill(1.0, d.np)

        CLM.cndv_establishment!(d.dgvs, d.eco, prec365,
                                 annsum_npp, annsum_litfall,
                                 deadstemc, leafcmax,
                                 d.pft, d.patch, d.lun, 1:d.np;
                                 bounds_gridcell=1:d.ng)

        # Trees should be introduced (PFTs 1,2 have twmax > 999)
        @test d.dgvs.present_patch[2] == true
        @test d.dgvs.present_patch[3] == true
        # Initial fpcgrid for trees = 0.000844
        @test d.dgvs.fpcgrid_patch[2] > 0.0
        @test d.dgvs.fpcgrid_patch[3] > 0.0

        # Grass should be introduced with fpcgrid = 0.05
        @test d.dgvs.present_patch[4] == true
    end

    @testset "cndv_establishment! tree establishment rate" begin
        d = make_cndv_test_data()

        # One tree present, one absent but suitable
        d.dgvs.present_patch .= [false, true, false, true]
        d.dgvs.nind_patch .= [0.0, 0.01, 0.0, 1.0]
        d.dgvs.fpcgrid_patch .= [0.0, 0.2, 0.0, 0.3]
        d.dgvs.crownarea_patch .= [0.0, 8.0, 0.0, 1.0]
        d.dgvs.pftmayexist_patch .= true

        d.dgvs.tmomin20_patch .= CLM.TFRZ + 0.0
        d.dgvs.agdd20_patch .= 500.0
        d.dgvs.agddtw_patch .= 0.0

        prec365 = [200.0 / (365.0 * CLM.SECSPDAY)]
        annsum_npp = fill(100.0, d.np)
        annsum_litfall = fill(50.0, d.np)
        deadstemc = [0.0, 5000.0, 0.1, 0.0]
        leafcmax = [0.0, 200.0, 1.0, 50.0]

        nind_before = d.dgvs.nind_patch[2]

        CLM.cndv_establishment!(d.dgvs, d.eco, prec365,
                                 annsum_npp, annsum_litfall,
                                 deadstemc, leafcmax,
                                 d.pft, d.patch, d.lun, 1:d.np;
                                 bounds_gridcell=1:d.ng)

        # Existing tree (patch 2) should gain individuals from sapling establishment
        @test d.dgvs.nind_patch[2] > nind_before
        # New tree (patch 3) should be introduced
        @test d.dgvs.present_patch[3] == true
    end

    @testset "cndv_establishment! mortality diagnostics" begin
        d = make_cndv_test_data()

        d.dgvs.present_patch .= [false, true, true, true]
        d.dgvs.nind_patch .= [0.0, 0.05, 0.03, 1.0]
        d.dgvs.fpcgrid_patch .= [0.0, 0.3, 0.3, 0.3]
        d.dgvs.crownarea_patch .= [0.0, 8.0, 8.0, 1.0]
        d.dgvs.pftmayexist_patch .= true
        d.dgvs.tmomin20_patch .= CLM.TFRZ + 0.0
        d.dgvs.agdd20_patch .= 500.0
        d.dgvs.agddtw_patch .= 0.0

        prec365 = [200.0 / (365.0 * CLM.SECSPDAY)]
        annsum_npp = fill(500.0, d.np)
        annsum_litfall = fill(100.0, d.np)
        deadstemc = [0.0, 5000.0, 3000.0, 0.0]
        leafcmax = [0.0, 200.0, 150.0, 50.0]

        CLM.cndv_establishment!(d.dgvs, d.eco, prec365,
                                 annsum_npp, annsum_litfall,
                                 deadstemc, leafcmax,
                                 d.pft, d.patch, d.lun, 1:d.np;
                                 bounds_gridcell=1:d.ng)

        # Growth efficiency should be computed for trees
        @test d.dgvs.greffic_patch[2] >= 0.0
        @test d.dgvs.greffic_patch[3] >= 0.0

        # Heatstress should be 0 when twmax > 999 (no heat limit)
        @test d.dgvs.heatstress_patch[2] ≈ 0.0
        @test d.dgvs.heatstress_patch[3] ≈ 0.0
    end

    # -----------------------------------------------------------------------
    # 8. cndv_driver! — annual CNDV cycle
    # -----------------------------------------------------------------------
    @testset "cndv_driver! annual cycle" begin
        d = make_cndv_test_data()

        # Set up initial state
        d.dgvs.present_patch .= [false, true, true, true]
        d.dgvs.nind_patch .= [0.0, 0.05, 0.03, 1.0]
        d.dgvs.fpcgrid_patch .= [0.0, 0.3, 0.3, 0.3]
        d.dgvs.crownarea_patch .= [0.0, 8.0, 8.0, 1.0]
        d.dgvs.pftmayexist_patch .= true
        d.dgvs.tmomin20_patch .= CLM.TFRZ + 0.0
        d.dgvs.agdd20_patch .= 500.0
        d.dgvs.agdd_patch .= 600.0

        mask_natvegp = trues(d.np)
        t_mo_min = fill(CLM.TFRZ - 5.0, d.np)
        prec365 = [200.0 / (365.0 * CLM.SECSPDAY)]
        annsum_npp = fill(200.0, d.np)
        annsum_litfall = fill(80.0, d.np)
        deadstemc = [0.0, 5000.0, 3000.0, 0.0]
        leafcmax = [0.0, 200.0, 150.0, 50.0]

        CLM.cndv_driver!(d.dgvs, d.eco, t_mo_min, prec365,
                          annsum_npp, annsum_litfall,
                          deadstemc, leafcmax,
                          d.pft, d.patch, d.lun, mask_natvegp,
                          1:d.np, 3;  # kyr=3 (not 2, so running means updated)
                          bounds_gridcell=1:d.ng)

        # Climate20 running means should be updated
        # tmomin20 = (19 * tmomin20_old + t_mo_min) / 20
        expected_tmomin20 = (19.0 * CLM.TFRZ + (CLM.TFRZ - 5.0)) / 20.0
        @test d.dgvs.tmomin20_patch[1] ≈ expected_tmomin20 atol=0.1

        # Annual reset: leafcmax should be zeroed
        @test leafcmax[2] ≈ 0.0
        # t_mo_min should be reset to 1e36
        @test t_mo_min[1] ≈ 1.0e36
    end

    @testset "cndv_driver! year 2 initialization" begin
        d = make_cndv_test_data()

        d.dgvs.present_patch .= [false, true, true, true]
        d.dgvs.nind_patch .= [0.0, 0.05, 0.03, 1.0]
        d.dgvs.fpcgrid_patch .= [0.0, 0.3, 0.3, 0.3]
        d.dgvs.crownarea_patch .= [0.0, 8.0, 8.0, 1.0]
        d.dgvs.pftmayexist_patch .= true
        d.dgvs.tmomin20_patch .= 0.0  # will be overwritten at kyr=2
        d.dgvs.agdd20_patch .= 0.0    # will be overwritten at kyr=2
        d.dgvs.agdd_patch .= 600.0

        mask_natvegp = trues(d.np)
        t_mo_min = fill(CLM.TFRZ - 5.0, d.np)
        prec365 = [200.0 / (365.0 * CLM.SECSPDAY)]
        annsum_npp = fill(200.0, d.np)
        annsum_litfall = fill(80.0, d.np)
        deadstemc = [0.0, 5000.0, 3000.0, 0.0]
        leafcmax = [0.0, 200.0, 150.0, 50.0]

        CLM.cndv_driver!(d.dgvs, d.eco, t_mo_min, prec365,
                          annsum_npp, annsum_litfall,
                          deadstemc, leafcmax,
                          d.pft, d.patch, d.lun, mask_natvegp,
                          1:d.np, 2;  # kyr=2 triggers direct init
                          bounds_gridcell=1:d.ng)

        # At kyr=2, tmomin20 is set to t_mo_min before running-mean update
        # So tmomin20 = (19 * t_mo_min + t_mo_min) / 20 = t_mo_min
        @test d.dgvs.tmomin20_patch[1] ≈ CLM.TFRZ - 5.0 atol=0.01
        # agdd20 = (19 * agdd + agdd) / 20 = agdd
        @test d.dgvs.agdd20_patch[1] ≈ 600.0 atol=0.01
    end

end
