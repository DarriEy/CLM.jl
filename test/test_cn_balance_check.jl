@testset "CN Balance Check (CNBalanceCheckMod)" begin

    # -----------------------------------------------------------------------
    # Helper: create minimal test data for balance checks
    # -----------------------------------------------------------------------
    function make_balance_test_data(; nc=3, ng=1)
        bal = CLM.CNBalanceData()
        CLM.cn_balance_init!(bal, nc, ng)

        # Column data
        col_data = CLM.ColumnData()
        col_data.gridcell = fill(1, nc)
        col_data.wtgcell  = fill(1.0 / nc, nc)
        col_data.is_fates = fill(false, nc)

        # Gridcell data
        grc_data = CLM.GridcellData()
        grc_data.latdeg = [45.0]
        grc_data.londeg = [-90.0]

        # SoilBiogeochem carbon state
        soilbgc_cstate = CLM.SoilBiogeochemCarbonStateData()
        soilbgc_cstate.totc_col = fill(100.0, nc)
        soilbgc_cstate.totc_grc = fill(100.0, ng)

        # SoilBiogeochem nitrogen state
        soilbgc_nstate = CLM.SoilBiogeochemNitrogenStateData()
        soilbgc_nstate.totn_col = fill(10.0, nc)
        soilbgc_nstate.totn_grc = fill(10.0, ng)

        # SoilBiogeochem carbon flux
        soilbgc_cflux = CLM.SoilBiogeochemCarbonFluxData()
        soilbgc_cflux.hr_col = fill(0.0, nc)
        soilbgc_cflux.som_c_leached_col = fill(0.0, nc)
        soilbgc_cflux.fates_litter_flux = fill(0.0, nc)

        # SoilBiogeochem nitrogen flux
        soilbgc_nflux = CLM.SoilBiogeochemNitrogenFluxData()
        soilbgc_nflux.ndep_to_sminn_col = fill(0.0, nc)
        soilbgc_nflux.nfix_to_sminn_col = fill(0.0, nc)
        soilbgc_nflux.ffix_to_sminn_col = fill(0.0, nc)
        soilbgc_nflux.fert_to_sminn_col = fill(0.0, nc)
        soilbgc_nflux.soyfixn_to_sminn_col = fill(0.0, nc)
        soilbgc_nflux.supplement_to_sminn_col = fill(0.0, nc)
        soilbgc_nflux.denit_col = fill(0.0, nc)
        soilbgc_nflux.sminn_leached_col = fill(0.0, nc)
        soilbgc_nflux.smin_no3_leached_col = fill(0.0, nc)
        soilbgc_nflux.smin_no3_runoff_col = fill(0.0, nc)
        soilbgc_nflux.f_n2o_nit_col = fill(0.0, nc)
        soilbgc_nflux.som_n_leached_col = fill(0.0, nc)
        soilbgc_nflux.sminn_to_plant_col = fill(0.0, nc)
        soilbgc_nflux.fates_litter_flux = fill(0.0, nc)

        # CNVeg carbon flux
        cnveg_cflux = CLM.CNVegCarbonFluxData()
        cnveg_cflux.gpp_col = fill(0.0, nc)
        cnveg_cflux.er_col = fill(0.0, nc)
        cnveg_cflux.fire_closs_col = fill(0.0, nc)
        cnveg_cflux.hrv_xsmrpool_to_atm_col = fill(0.0, nc)
        cnveg_cflux.xsmrpool_to_atm_col = fill(0.0, nc)
        cnveg_cflux.wood_harvestc_col = fill(0.0, nc)
        cnveg_cflux.gru_conv_cflux_col = fill(0.0, nc)
        cnveg_cflux.gru_wood_productc_gain_col = fill(0.0, nc)
        cnveg_cflux.crop_harvestc_to_cropprodc_col = fill(0.0, nc)
        cnveg_cflux.nbp_grc = fill(0.0, ng)
        cnveg_cflux.dwt_seedc_to_leaf_grc = fill(0.0, ng)
        cnveg_cflux.dwt_seedc_to_deadstem_grc = fill(0.0, ng)

        # CNVeg nitrogen flux
        cnveg_nflux = CLM.CNVegNitrogenFluxData()
        cnveg_nflux.fire_nloss_col = fill(0.0, nc)
        cnveg_nflux.wood_harvestn_col = fill(0.0, nc)
        cnveg_nflux.gru_conv_nflux_col = fill(0.0, nc)
        cnveg_nflux.gru_conv_nflux_grc = fill(0.0, ng)
        cnveg_nflux.gru_wood_productn_gain_col = fill(0.0, nc)
        cnveg_nflux.gru_wood_productn_gain_grc = fill(0.0, ng)
        cnveg_nflux.crop_harvestn_to_cropprodn_col = fill(0.0, nc)
        cnveg_nflux.dwt_seedn_to_leaf_grc = fill(0.0, ng)
        cnveg_nflux.dwt_seedn_to_deadstem_grc = fill(0.0, ng)
        cnveg_nflux.dwt_conv_nflux_grc = fill(0.0, ng)

        # Products
        c_products = CLM.CNProductsData()
        CLM.cn_products_init!(c_products, ng)
        n_products = CLM.CNProductsData()
        CLM.cn_products_init!(n_products, ng)

        mask_soil = trues(nc)
        bounds_c = 1:nc
        bounds_g = 1:ng
        dt = 1800.0

        return (bal=bal, col_data=col_data, grc_data=grc_data,
                soilbgc_cstate=soilbgc_cstate, soilbgc_nstate=soilbgc_nstate,
                soilbgc_cflux=soilbgc_cflux, soilbgc_nflux=soilbgc_nflux,
                cnveg_cflux=cnveg_cflux, cnveg_nflux=cnveg_nflux,
                c_products=c_products, n_products=n_products,
                mask_soil=mask_soil, bounds_c=bounds_c, bounds_g=bounds_g, dt=dt)
    end

    # -----------------------------------------------------------------------
    # Test CNBalanceData initialization
    # -----------------------------------------------------------------------
    @testset "cn_balance_init!" begin
        bal = CLM.CNBalanceData()
        CLM.cn_balance_init!(bal, 5, 2)

        @test length(bal.begcb_col) == 5
        @test length(bal.endcb_col) == 5
        @test length(bal.begnb_col) == 5
        @test length(bal.endnb_col) == 5
        @test length(bal.begcb_grc) == 2
        @test length(bal.endcb_grc) == 2
        @test length(bal.begnb_grc) == 2
        @test length(bal.endnb_grc) == 2

        @test all(isnan, bal.begcb_col)
        @test all(isnan, bal.begcb_grc)

        @test bal.cwarning == 1.0e-8
        @test bal.nwarning == 1.0e-7
        @test bal.cerror == 1.0e-7
        @test bal.nerror == 1.0e-3
    end

    # -----------------------------------------------------------------------
    # Test CNProductsData initialization
    # -----------------------------------------------------------------------
    @testset "cn_products_init!" begin
        prod = CLM.CNProductsData()
        CLM.cn_products_init!(prod, 3)

        @test length(prod.cropprod1_grc) == 3
        @test length(prod.tot_woodprod_grc) == 3
        @test length(prod.product_loss_grc) == 3
        @test all(prod.cropprod1_grc .== 0.0)
    end

    # -----------------------------------------------------------------------
    # Test begin_cn_gridcell_balance!
    # -----------------------------------------------------------------------
    @testset "begin_cn_gridcell_balance! without FATES" begin
        d = make_balance_test_data()
        d.soilbgc_cstate.totc_grc[1] = 500.0
        d.soilbgc_nstate.totn_grc[1] = 50.0
        d.c_products.tot_woodprod_grc[1] = 10.0
        d.c_products.cropprod1_grc[1] = 5.0
        d.n_products.tot_woodprod_grc[1] = 1.0
        d.n_products.cropprod1_grc[1] = 0.5

        hrv = [2.0]
        gru = [3.0]
        dwt = [4.0]

        CLM.begin_cn_gridcell_balance!(d.bal, d.soilbgc_cstate, d.soilbgc_nstate,
            d.c_products, d.n_products, d.bounds_g;
            use_fates_bgc=false,
            hrv_xsmrpool_amount_left=hrv,
            gru_conv_cflux_amount_left=gru,
            dwt_conv_cflux_amount_left=dwt)

        @test d.bal.begcb_grc[1] ≈ 500.0 + 10.0 + 5.0 + 2.0 + 3.0 + 4.0
        @test d.bal.begnb_grc[1] ≈ 50.0 + 1.0 + 0.5
    end

    @testset "begin_cn_gridcell_balance! with FATES" begin
        d = make_balance_test_data()
        d.soilbgc_cstate.totc_grc[1] = 500.0
        d.soilbgc_nstate.totn_grc[1] = 50.0
        d.c_products.tot_woodprod_grc[1] = 10.0
        d.c_products.cropprod1_grc[1] = 5.0

        CLM.begin_cn_gridcell_balance!(d.bal, d.soilbgc_cstate, d.soilbgc_nstate,
            d.c_products, d.n_products, d.bounds_g;
            use_fates_bgc=true)

        # With FATES, dribbler amounts are not added
        @test d.bal.begcb_grc[1] ≈ 500.0 + 10.0 + 5.0
    end

    # -----------------------------------------------------------------------
    # Test begin_cn_column_balance!
    # -----------------------------------------------------------------------
    @testset "begin_cn_column_balance!" begin
        d = make_balance_test_data()
        d.soilbgc_cstate.totc_col .= [100.0, 200.0, 300.0]
        d.soilbgc_nstate.totn_col .= [10.0, 20.0, 30.0]

        CLM.begin_cn_column_balance!(d.bal, d.soilbgc_cstate, d.soilbgc_nstate,
            d.mask_soil, d.bounds_c)

        @test d.bal.begcb_col[1] ≈ 100.0
        @test d.bal.begcb_col[2] ≈ 200.0
        @test d.bal.begcb_col[3] ≈ 300.0
        @test d.bal.begnb_col[1] ≈ 10.0
        @test d.bal.begnb_col[2] ≈ 20.0
        @test d.bal.begnb_col[3] ≈ 30.0
    end

    @testset "begin_cn_column_balance! with mask" begin
        d = make_balance_test_data()
        d.soilbgc_cstate.totc_col .= [100.0, 200.0, 300.0]
        d.soilbgc_nstate.totn_col .= [10.0, 20.0, 30.0]
        d.mask_soil[2] = false

        CLM.begin_cn_column_balance!(d.bal, d.soilbgc_cstate, d.soilbgc_nstate,
            d.mask_soil, d.bounds_c)

        @test d.bal.begcb_col[1] ≈ 100.0
        @test isnan(d.bal.begcb_col[2])  # masked out, stays NaN
        @test d.bal.begcb_col[3] ≈ 300.0
    end

    # -----------------------------------------------------------------------
    # Test c_balance_check! — balanced case (no error)
    # -----------------------------------------------------------------------
    @testset "c_balance_check! balanced (no error)" begin
        d = make_balance_test_data()

        # Set up a balanced scenario:
        # begcb = totcolc_begin = 100 for all columns
        CLM.begin_cn_column_balance!(d.bal, d.soilbgc_cstate, d.soilbgc_nstate,
            d.mask_soil, d.bounds_c)

        # Set gridcell beginning balance (use FATES path to skip dribblers)
        CLM.begin_cn_gridcell_balance!(d.bal, d.soilbgc_cstate, d.soilbgc_nstate,
            d.c_products, d.n_products, d.bounds_g;
            use_fates_bgc=true)

        # GPP = 0.01 gC/m2/s, ER = 0.01 gC/m2/s → net zero change
        # totc_col stays the same (100.0)
        d.cnveg_cflux.gpp_col .= 0.01
        d.cnveg_cflux.er_col .= 0.01

        # Should not throw (balance is zero)
        CLM.c_balance_check!(d.bal, d.soilbgc_cflux, d.soilbgc_cstate,
            d.cnveg_cflux, d.c_products, d.col_data, d.grc_data,
            d.mask_soil, d.bounds_c, d.bounds_g, d.dt;
            is_fates_col=d.col_data.is_fates,
            use_fates_bgc=true)

        # endcb should equal totcolc
        @test d.bal.endcb_col[1] ≈ 100.0
        @test d.bal.endcb_col[2] ≈ 100.0
        @test d.bal.endcb_col[3] ≈ 100.0
    end

    # -----------------------------------------------------------------------
    # Test c_balance_check! — with GPP changing storage
    # -----------------------------------------------------------------------
    @testset "c_balance_check! GPP changes storage" begin
        d = make_balance_test_data()
        dt = d.dt

        # begin balance
        CLM.begin_cn_column_balance!(d.bal, d.soilbgc_cstate, d.soilbgc_nstate,
            d.mask_soil, d.bounds_c)
        CLM.begin_cn_gridcell_balance!(d.bal, d.soilbgc_cstate, d.soilbgc_nstate,
            d.c_products, d.n_products, d.bounds_g;
            use_fates_bgc=true)

        # GPP input, no output → storage should increase by gpp*dt
        gpp_rate = 0.05
        d.cnveg_cflux.gpp_col .= gpp_rate

        # After timestep, totc_col increases by gpp*dt
        d.soilbgc_cstate.totc_col .= 100.0 .+ gpp_rate * dt

        # Should not throw
        CLM.c_balance_check!(d.bal, d.soilbgc_cflux, d.soilbgc_cstate,
            d.cnveg_cflux, d.c_products, d.col_data, d.grc_data,
            d.mask_soil, d.bounds_c, d.bounds_g, d.dt;
            is_fates_col=d.col_data.is_fates,
            use_fates_bgc=true)

        @test d.bal.endcb_col[1] ≈ 100.0 + gpp_rate * dt
    end

    # -----------------------------------------------------------------------
    # Test c_balance_check! — error when balance is violated
    # -----------------------------------------------------------------------
    @testset "c_balance_check! error on large imbalance" begin
        d = make_balance_test_data()

        CLM.begin_cn_column_balance!(d.bal, d.soilbgc_cstate, d.soilbgc_nstate,
            d.mask_soil, d.bounds_c)
        CLM.begin_cn_gridcell_balance!(d.bal, d.soilbgc_cstate, d.soilbgc_nstate,
            d.c_products, d.n_products, d.bounds_g;
            use_fates_bgc=true)

        # Large GPP but no change in storage → imbalance
        d.cnveg_cflux.gpp_col .= 1.0  # large input
        # totc_col stays at 100.0 (no change)

        @test_throws ErrorException CLM.c_balance_check!(
            d.bal, d.soilbgc_cflux, d.soilbgc_cstate,
            d.cnveg_cflux, d.c_products, d.col_data, d.grc_data,
            d.mask_soil, d.bounds_c, d.bounds_g, d.dt;
            is_fates_col=d.col_data.is_fates,
            use_fates_bgc=true)
    end

    # -----------------------------------------------------------------------
    # Test c_balance_check! — FATES column
    # -----------------------------------------------------------------------
    @testset "c_balance_check! FATES column balanced" begin
        d = make_balance_test_data(nc=1, ng=1)
        d.col_data.is_fates = [true]

        CLM.begin_cn_column_balance!(d.bal, d.soilbgc_cstate, d.soilbgc_nstate,
            d.mask_soil, d.bounds_c)
        CLM.begin_cn_gridcell_balance!(d.bal, d.soilbgc_cstate, d.soilbgc_nstate,
            d.c_products, d.n_products, d.bounds_g;
            use_fates_bgc=true)

        # FATES: litter input matched by HR
        litter_rate = 0.02
        d.soilbgc_cflux.fates_litter_flux[1] = litter_rate
        d.soilbgc_cflux.hr_col[1] = litter_rate
        # Storage unchanged

        CLM.c_balance_check!(d.bal, d.soilbgc_cflux, d.soilbgc_cstate,
            d.cnveg_cflux, d.c_products, d.col_data, d.grc_data,
            d.mask_soil, d.bounds_c, d.bounds_g, d.dt;
            is_fates_col=d.col_data.is_fates,
            use_fates_bgc=true)

        @test d.bal.endcb_col[1] ≈ 100.0
    end

    # -----------------------------------------------------------------------
    # Test n_balance_check! — balanced case
    # -----------------------------------------------------------------------
    @testset "n_balance_check! balanced (no error)" begin
        d = make_balance_test_data()

        CLM.begin_cn_column_balance!(d.bal, d.soilbgc_cstate, d.soilbgc_nstate,
            d.mask_soil, d.bounds_c)
        CLM.begin_cn_gridcell_balance!(d.bal, d.soilbgc_cstate, d.soilbgc_nstate,
            d.c_products, d.n_products, d.bounds_g;
            use_fates_bgc=true)

        # ndep input matched by denit output
        ndep_rate = 0.001
        d.soilbgc_nflux.ndep_to_sminn_col .= ndep_rate
        d.soilbgc_nflux.denit_col .= ndep_rate
        # Storage unchanged

        CLM.n_balance_check!(d.bal, d.soilbgc_nflux, d.soilbgc_nstate,
            d.cnveg_nflux, d.n_products, d.col_data, d.grc_data,
            d.mask_soil, d.bounds_c, d.bounds_g, d.dt;
            is_fates_col=d.col_data.is_fates,
            use_fates_bgc=true)

        @test d.bal.endnb_col[1] ≈ 10.0
    end

    # -----------------------------------------------------------------------
    # Test n_balance_check! — N deposition increases storage
    # -----------------------------------------------------------------------
    @testset "n_balance_check! ndep increases storage" begin
        d = make_balance_test_data()
        dt = d.dt

        CLM.begin_cn_column_balance!(d.bal, d.soilbgc_cstate, d.soilbgc_nstate,
            d.mask_soil, d.bounds_c)
        CLM.begin_cn_gridcell_balance!(d.bal, d.soilbgc_cstate, d.soilbgc_nstate,
            d.c_products, d.n_products, d.bounds_g;
            use_fates_bgc=true)

        ndep_rate = 0.002
        d.soilbgc_nflux.ndep_to_sminn_col .= ndep_rate
        d.soilbgc_nstate.totn_col .= 10.0 .+ ndep_rate * dt

        CLM.n_balance_check!(d.bal, d.soilbgc_nflux, d.soilbgc_nstate,
            d.cnveg_nflux, d.n_products, d.col_data, d.grc_data,
            d.mask_soil, d.bounds_c, d.bounds_g, d.dt;
            is_fates_col=d.col_data.is_fates,
            use_fates_bgc=true)

        @test d.bal.endnb_col[1] ≈ 10.0 + ndep_rate * dt
    end

    # -----------------------------------------------------------------------
    # Test n_balance_check! — error on large imbalance
    # -----------------------------------------------------------------------
    @testset "n_balance_check! error on large imbalance" begin
        d = make_balance_test_data()

        CLM.begin_cn_column_balance!(d.bal, d.soilbgc_cstate, d.soilbgc_nstate,
            d.mask_soil, d.bounds_c)
        CLM.begin_cn_gridcell_balance!(d.bal, d.soilbgc_cstate, d.soilbgc_nstate,
            d.c_products, d.n_products, d.bounds_g;
            use_fates_bgc=true)

        # Large N input but no change in storage → imbalance > nerror
        d.soilbgc_nflux.ndep_to_sminn_col .= 1.0

        @test_throws ErrorException CLM.n_balance_check!(
            d.bal, d.soilbgc_nflux, d.soilbgc_nstate,
            d.cnveg_nflux, d.n_products, d.col_data, d.grc_data,
            d.mask_soil, d.bounds_c, d.bounds_g, d.dt;
            is_fates_col=d.col_data.is_fates,
            use_fates_bgc=true)
    end

    # -----------------------------------------------------------------------
    # Test n_balance_check! with nitrif_denitrif paths
    # -----------------------------------------------------------------------
    @testset "n_balance_check! nitrif_denitrif=false" begin
        d = make_balance_test_data()
        dt = d.dt

        CLM.begin_cn_column_balance!(d.bal, d.soilbgc_cstate, d.soilbgc_nstate,
            d.mask_soil, d.bounds_c)
        CLM.begin_cn_gridcell_balance!(d.bal, d.soilbgc_cstate, d.soilbgc_nstate,
            d.c_products, d.n_products, d.bounds_g;
            use_fates_bgc=true)

        # Use sminn_leached path instead of smin_no3_leached
        leach_rate = 0.001
        d.soilbgc_nflux.ndep_to_sminn_col .= leach_rate
        d.soilbgc_nflux.sminn_leached_col .= leach_rate
        # Storage unchanged (input = output)

        CLM.n_balance_check!(d.bal, d.soilbgc_nflux, d.soilbgc_nstate,
            d.cnveg_nflux, d.n_products, d.col_data, d.grc_data,
            d.mask_soil, d.bounds_c, d.bounds_g, d.dt;
            is_fates_col=d.col_data.is_fates,
            use_fates_bgc=true,
            use_nitrif_denitrif=false)

        @test d.bal.endnb_col[1] ≈ 10.0
    end

    # -----------------------------------------------------------------------
    # Test c2g_unity! aggregation
    # -----------------------------------------------------------------------
    @testset "c2g_unity! aggregation" begin
        nc = 4
        ng = 2
        garr = zeros(ng)
        carr = [10.0, 20.0, 30.0, 40.0]
        col_gridcell = [1, 1, 2, 2]
        col_wtgcell = [0.3, 0.7, 0.5, 0.5]

        CLM.c2g_unity!(garr, carr, col_gridcell, col_wtgcell, 1:nc, 1:ng)

        @test garr[1] ≈ 10.0 * 0.3 + 20.0 * 0.7
        @test garr[2] ≈ 30.0 * 0.5 + 40.0 * 0.5
    end

    # -----------------------------------------------------------------------
    # Test n_balance_check! — gridcell-level (non-FATES) path
    # -----------------------------------------------------------------------
    @testset "n_balance_check! gridcell-level non-FATES balanced" begin
        d = make_balance_test_data()
        dt = d.dt

        CLM.begin_cn_column_balance!(d.bal, d.soilbgc_cstate, d.soilbgc_nstate,
            d.mask_soil, d.bounds_c)
        CLM.begin_cn_gridcell_balance!(d.bal, d.soilbgc_cstate, d.soilbgc_nstate,
            d.c_products, d.n_products, d.bounds_g;
            use_fates_bgc=false,
            hrv_xsmrpool_amount_left=[0.0],
            gru_conv_cflux_amount_left=[0.0],
            dwt_conv_cflux_amount_left=[0.0])

        # ndep input matched by denit output at column level
        ndep_rate = 0.001
        d.soilbgc_nflux.ndep_to_sminn_col .= ndep_rate
        d.soilbgc_nflux.denit_col .= ndep_rate
        # Storage unchanged (input = output)

        # Set gridcell-level product fields to zero to keep it simple
        d.n_products.product_loss_grc .= 0.0

        CLM.n_balance_check!(d.bal, d.soilbgc_nflux, d.soilbgc_nstate,
            d.cnveg_nflux, d.n_products, d.col_data, d.grc_data,
            d.mask_soil, d.bounds_c, d.bounds_g, d.dt;
            is_fates_col=d.col_data.is_fates,
            use_fates_bgc=false)

        # endnb_grc should be set correctly
        @test d.bal.endnb_grc[1] ≈ d.soilbgc_nstate.totn_grc[1]
    end

    # -----------------------------------------------------------------------
    # Test masked columns are skipped properly
    # -----------------------------------------------------------------------
    @testset "masked columns skipped" begin
        d = make_balance_test_data()
        d.mask_soil[2] = false

        CLM.begin_cn_column_balance!(d.bal, d.soilbgc_cstate, d.soilbgc_nstate,
            d.mask_soil, d.bounds_c)

        @test d.bal.begcb_col[1] ≈ 100.0
        @test isnan(d.bal.begcb_col[2])  # masked
        @test d.bal.begcb_col[3] ≈ 100.0
    end

end
