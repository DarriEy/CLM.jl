#!/usr/bin/env julia
# ==========================================================================
# GRIDCELL C/N balance closes with the wood/crop PRODUCT POOLS wired.
#
# Companion to test_cn_balance_check.jl. Proves the fix for the open gridcell
# residual: `cn_products_update!` (+ the ComputeProductSummary / ComputeSummary
# pair) is now run on the live path (cn_driver_no_leaching!, CNDriverMod.F90:
# 851-912), so tot_woodprod_grc / cropprod1_grc advance from the harvest and
# land-conversion product-gain fluxes. Those pools are the SINK the gridcell
# balance check subtracts, so without the update the gridcell store loses the
# harvested C/N with nothing to balance it — an errcb_grc / errnb_grc residual
# equal to the harvest flux × dt.
#
# Scenario (the cleanest gridcell-internal transfer): a single g/c/p domain where
# wood harvest moves C and N out of the column store (totcolc/totcoln drop by
# H*dt) and INTO the wood product pool (tot_woodprod_grc grows by H*dt). At the
# gridcell level this is an internal transfer — net store unchanged — so the
# balance MUST close to machine precision once the product pool is advanced.
#
#   * WIRED (product update run): grc_endcb = (totgrcc - Hc*dt) + Hc*dt = begcb,
#     errcb_grc ≈ 0 and errnb_grc ≈ 0 — the FATAL non-fates gridcell check passes.
#   * UNWIRED (product update skipped, tot_woodprod stays 0): grc_endcb loses the
#     Hc*dt with no offset → the fatal check TRIPS, and the residual is exactly
#     Hc*dt (C) / Hn*dt (N) — the ~37 gC / ~0.085 gN class residual.
# ==========================================================================

using Test
using CLM

@testset "gridcell C/N balance closes with product pools wired" begin

    ng, nc, np = 1, 1, 1
    dt = 1800.0

    # Harvest fluxes (g/m2/s), sized so H*dt lands in the ~37 gC / ~0.085 gN class
    # of the reported open residual.
    Hc = 0.0205      # gC/m2/s  -> 36.9 gC over the step
    Hn = 4.7e-5      # gN/m2/s  -> 0.0846 gN over the step

    totc0 = 100.0    # gC/m2 column store, start of step
    totn0 = 10.0     # gN/m2 column store, start of step

    # ---- subgrid geometry: one gridcell, one column, one natural-veg patch ----
    function build_domain()
        col = CLM.ColumnData()
        col.gridcell = [1]; col.landunit = [1]
        col.wtgcell = [1.0]; col.wtlunit = [1.0]
        col.active = [true]; col.is_fates = [false]

        grc = CLM.GridcellData()
        grc.latdeg = [45.0]; grc.londeg = [-90.0]

        pch = CLM.PatchData()
        pch.gridcell = [1]; pch.column = [1]; pch.landunit = [1]
        pch.itype = [1]
        pch.wtcol = [1.0]; pch.wtgcell = [1.0]; pch.wtlunit = [1.0]
        pch.active = [true]

        lun = CLM.LandunitData()
        lun.itype = [CLM.ISTSOIL]; lun.wtgcell = [1.0]; lun.active = [true]
        lun.coli = [1]; lun.colf = [1]

        # pftcon: pprod{10,100,harv10} indexed by itype+1 (=2 here).
        pft = CLM.PftconType()
        pft.pprod10     = [0.0, 0.5]
        pft.pprod100    = [0.0, 0.3]
        pft.pprodharv10 = [0.0, 0.5]
        return (; col, grc, pch, lun, pft)
    end

    # ---- CN state/flux + product/balance instances, harvest-only fluxes -------
    function build_state()
        bal = CLM.CNBalanceData(); CLM.cn_balance_init!(bal, nc, ng)

        cstate = CLM.SoilBiogeochemCarbonStateData()
        cstate.totc_col = fill(totc0, nc); cstate.totc_grc = fill(totc0, ng)
        nstate = CLM.SoilBiogeochemNitrogenStateData()
        nstate.totn_col = fill(totn0, nc); nstate.totn_grc = fill(totn0, ng)

        scf = CLM.SoilBiogeochemCarbonFluxData()
        scf.hr_col = zeros(nc); scf.som_c_leached_col = zeros(nc)
        scf.fates_litter_flux = zeros(nc)
        snf = CLM.SoilBiogeochemNitrogenFluxData()
        for f in (:ndep_to_sminn_col, :nfix_to_sminn_col, :ffix_to_sminn_col,
                  :fert_to_sminn_col, :soyfixn_to_sminn_col, :supplement_to_sminn_col,
                  :denit_col, :sminn_leached_col, :smin_no3_leached_col,
                  :smin_no3_runoff_col, :f_n2o_nit_col, :som_n_leached_col,
                  :sminn_to_plant_col, :fates_litter_flux)
            setproperty!(snf, f, zeros(nc))
        end

        cf = CLM.CNVegCarbonFluxData()
        CLM.cnveg_carbon_flux_init!(cf, np, nc, ng)
        for f in (:gpp_col, :er_col, :fire_closs_col, :hrv_xsmrpool_to_atm_col,
                  :xsmrpool_to_atm_col, :wood_harvestc_col, :gru_conv_cflux_col,
                  :gru_wood_productc_gain_col, :crop_harvestc_to_cropprodc_col)
            setproperty!(cf, f, zeros(nc))
        end
        cf.nbp_grc = zeros(ng)
        cf.dwt_seedc_to_leaf_grc = zeros(ng)
        cf.dwt_seedc_to_deadstem_grc = zeros(ng)
        # Zero the product-gain PATCH inputs, then set only wood harvest.
        for f in (:dwt_wood_productc_gain_patch, :gru_wood_productc_gain_patch,
                  :wood_harvestc_patch, :dwt_crop_productc_gain_patch,
                  :crop_harvestc_to_cropprodc_patch)
            setproperty!(cf, f, zeros(np))
        end
        cf.wood_harvestc_col[1]   = Hc     # column output the balance reads
        cf.wood_harvestc_patch[1] = Hc     # product input the pool grows by

        nf = CLM.CNVegNitrogenFluxData()
        CLM.cnveg_nitrogen_flux_init!(nf, np, nc, ng)
        for f in (:fire_nloss_col, :wood_harvestn_col, :gru_conv_nflux_col,
                  :gru_wood_productn_gain_col, :crop_harvestn_to_cropprodn_col)
            setproperty!(nf, f, zeros(nc))
        end
        for f in (:gru_conv_nflux_grc, :gru_wood_productn_gain_grc,
                  :dwt_conv_nflux_grc, :dwt_seedn_to_leaf_grc,
                  :dwt_seedn_to_deadstem_grc)
            setproperty!(nf, f, zeros(ng))
        end
        for f in (:dwt_wood_productn_gain_patch, :gru_wood_productn_gain_patch,
                  :wood_harvestn_patch, :dwt_crop_productn_gain_patch,
                  :crop_harvestn_to_cropprodn_patch)
            setproperty!(nf, f, zeros(np))
        end
        nf.wood_harvestn_col[1]   = Hn
        nf.wood_harvestn_patch[1] = Hn

        c_prod = CLM.CNProductsFullData(); CLM.cn_products_full_init!(c_prod, ng, np)
        n_prod = CLM.CNProductsFullData(); CLM.cn_products_full_init!(n_prod, ng, np)

        return (; bal, cstate, nstate, scf, snf, cf, nf, c_prod, n_prod)
    end

    # ---- the CNWoodProducts sequence, exactly as cn_driver_no_leaching! runs it
    function run_products!(s, dom)
        bounds = CLM.BoundsType(begg=1, endg=ng, begc=1, endc=nc, begp=1, endp=np,
                                level=CLM.BOUNDS_LEVEL_PROC)
        filt = collect(1:np)   # single active BGC veg patch
        CLM.cn_products_set_values!(s.c_prod, bounds, 0.0)
        CLM.cn_products_set_values!(s.n_prod, bounds, 0.0)
        CLM.cn_products_update!(s.c_prod, bounds, length(filt), filt,
            s.cf.dwt_wood_productc_gain_patch, s.cf.gru_wood_productc_gain_patch,
            s.cf.wood_harvestc_patch, s.cf.dwt_crop_productc_gain_patch,
            s.cf.crop_harvestc_to_cropprodc_patch;
            pprod10=dom.pft.pprod10, pprod100=dom.pft.pprod100,
            pprodharv10=dom.pft.pprodharv10, patch_gridcell=dom.pch.gridcell,
            patch_itype=dom.pch.itype, pch=dom.pch, col=dom.col, lun=dom.lun)
        CLM.cn_products_update!(s.n_prod, bounds, length(filt), filt,
            s.nf.dwt_wood_productn_gain_patch, s.nf.gru_wood_productn_gain_patch,
            s.nf.wood_harvestn_patch, s.nf.dwt_crop_productn_gain_patch,
            s.nf.crop_harvestn_to_cropprodn_patch;
            pprod10=dom.pft.pprod10, pprod100=dom.pft.pprod100,
            pprodharv10=dom.pft.pprodharv10, patch_gridcell=dom.pch.gridcell,
            patch_itype=dom.pch.itype, pch=dom.pch, col=dom.col, lun=dom.lun)
        CLM.cn_products_compute_product_summary!(s.c_prod, bounds, dt)
        CLM.cn_products_compute_product_summary!(s.n_prod, bounds, dt)
        CLM.cn_products_compute_summary!(s.c_prod, bounds)
        CLM.cn_products_compute_summary!(s.n_prod, bounds)
        return nothing
    end

    # Record begin balances, apply harvest to the store, run the C & N checks.
    # `wire_products` toggles the product-pool advance (the code under test).
    function step!(s, dom; wire_products::Bool)
        mask = trues(nc); bc = 1:nc; bg = 1:ng
        zero_dribble = zeros(ng)

        CLM.begin_cn_column_balance!(s.bal, s.cstate, s.nstate, mask, bc)
        CLM.begin_cn_gridcell_balance!(s.bal, s.cstate, s.nstate,
            s.c_prod, s.n_prod, bg; use_fates_bgc=false,
            hrv_xsmrpool_amount_left=zero_dribble,
            gru_conv_cflux_amount_left=zero_dribble,
            dwt_conv_cflux_amount_left=zero_dribble)

        if wire_products
            run_products!(s, dom)
        end

        # CStateUpdate: harvest removed C/N from the column store this step.
        s.cstate.totc_col[1] = totc0 - Hc * dt
        s.nstate.totn_col[1] = totn0 - Hn * dt

        CLM.c_balance_check!(s.bal, s.scf, s.cstate, s.cf, s.c_prod,
            dom.col, dom.grc, mask, bc, bg, dt;
            is_fates_col=dom.col.is_fates, use_fates_bgc=false,
            hrv_xsmrpool_amount_left=zero_dribble,
            gru_conv_cflux_amount_left=zero_dribble,
            dwt_conv_cflux_amount_left=zero_dribble)
        CLM.n_balance_check!(s.bal, s.snf, s.nstate, s.nf, s.n_prod,
            dom.col, dom.grc, mask, bc, bg, dt;
            is_fates_col=dom.col.is_fates, use_fates_bgc=false,
            use_nitrif_denitrif=false)
        return nothing
    end

    # ----------------------------------------------------------------------
    # 1. WIRED: the product pool advances → gridcell balance closes (no throw).
    # ----------------------------------------------------------------------
    @testset "wired: gridcell C & N balance closes to machine precision" begin
        dom = build_domain()
        s = build_state()
        step!(s, dom; wire_products=true)

        # The harvest moved exactly Hc*dt / Hn*dt into the wood product pool.
        @test s.c_prod.tot_woodprod_grc[1] ≈ Hc * dt   rtol=1e-12
        @test s.n_prod.tot_woodprod_grc[1] ≈ Hn * dt   rtol=1e-12
        # First step: pools started at 0 so decay loss is 0.
        @test s.c_prod.product_loss_grc[1] == 0.0
        @test s.n_prod.product_loss_grc[1] == 0.0

        # Gridcell balance closes.
        @test abs(s.bal.errcb_grc[1]) < 1e-9
        @test abs(s.bal.errnb_grc[1]) < 1e-9
        # Column balance also closes (harvest is a clean column output).
        @test abs(s.bal.errcb_col[1]) < 1e-9
        @test abs(s.bal.errnb_col[1]) < 1e-9
    end

    # ----------------------------------------------------------------------
    # 2. UNWIRED: skipping the product advance leaves the fatal gridcell
    #    residual. The non-fates gridcell C check is FATAL, so `step!` throws;
    #    the residual it trips on is exactly Hc*dt (C) / Hn*dt (N) — the same
    #    harvest that the wired path (test 1) moved into the product pool and
    #    thereby balanced. Compare against the ~37 gC / ~0.085 gN class of the
    #    originally-reported open residual.
    # ----------------------------------------------------------------------
    @testset "unwired: gridcell residual is fatal and equals H*dt" begin
        dom = build_domain(); s = build_state()
        err = try
            step!(s, dom; wire_products=false)
            nothing
        catch e
            e
        end
        @test err isa ErrorException                 # fatal gridcell C check tripped
        @test occursin("gridcell carbon", sprint(showerror, err))

        residual_c = Hc * dt
        residual_n = Hn * dt
        @test residual_c ≈ 36.9   atol = 0.1         # ~37 gC class
        @test residual_n ≈ 0.0846 atol = 1e-3        # ~0.085 gN class
        println("  unwired gridcell residual: errcb_grc ≈ ", residual_c,
                " gC, errnb_grc ≈ ", residual_n, " gN (both CLOSE to 0 once wired)")
    end

    # ----------------------------------------------------------------------
    # 3. NO-OP on a non-transient, no-harvest run: with every product-gain
    #    flux at 0, running the product step must leave tot_woodprod_grc /
    #    cropprod1_grc / product_loss_grc at exactly 0.0 and the gridcell/column
    #    balance BYTE-IDENTICAL to not running it. This is the default-preserving
    #    guarantee (the wiring only moves mass when there is harvest / land
    #    conversion; otherwise it is invisible).
    # ----------------------------------------------------------------------
    @testset "no-op: no-harvest run is byte-identical with/without the product step" begin
        mask = trues(nc); bc = 1:nc; bg = 1:ng; zero_dribble = zeros(ng)

        function run_noharvest(; wire_products::Bool)
            dom = build_domain(); s = build_state()
            # No harvest anywhere → every product-gain flux is 0.
            s.cf.wood_harvestc_col[1] = 0.0; s.cf.wood_harvestc_patch[1] = 0.0
            s.nf.wood_harvestn_col[1] = 0.0; s.nf.wood_harvestn_patch[1] = 0.0

            CLM.begin_cn_column_balance!(s.bal, s.cstate, s.nstate, mask, bc)
            CLM.begin_cn_gridcell_balance!(s.bal, s.cstate, s.nstate,
                s.c_prod, s.n_prod, bg; use_fates_bgc=false,
                hrv_xsmrpool_amount_left=zero_dribble,
                gru_conv_cflux_amount_left=zero_dribble,
                dwt_conv_cflux_amount_left=zero_dribble)
            wire_products && run_products!(s, dom)
            # No CStateUpdate change: store unchanged (no fluxes).
            CLM.c_balance_check!(s.bal, s.scf, s.cstate, s.cf, s.c_prod,
                dom.col, dom.grc, mask, bc, bg, dt;
                is_fates_col=dom.col.is_fates, use_fates_bgc=false,
                hrv_xsmrpool_amount_left=zero_dribble,
                gru_conv_cflux_amount_left=zero_dribble,
                dwt_conv_cflux_amount_left=zero_dribble)
            CLM.n_balance_check!(s.bal, s.snf, s.nstate, s.nf, s.n_prod,
                dom.col, dom.grc, mask, bc, bg, dt;
                is_fates_col=dom.col.is_fates, use_fates_bgc=false,
                use_nitrif_denitrif=false)
            return s
        end

        s_off = run_noharvest(wire_products=false)
        s_on  = run_noharvest(wire_products=true)

        # Product pools untouched (exactly 0.0) when the step runs with no gains.
        @test s_on.c_prod.tot_woodprod_grc[1] === 0.0
        @test s_on.c_prod.cropprod1_grc[1]    === 0.0
        @test s_on.c_prod.product_loss_grc[1] === 0.0
        @test s_on.n_prod.tot_woodprod_grc[1] === 0.0
        @test s_on.n_prod.product_loss_grc[1] === 0.0

        # Balance errors BYTE-IDENTICAL (===) with and without the product step.
        @test s_on.bal.errcb_grc[1] === s_off.bal.errcb_grc[1]
        @test s_on.bal.errnb_grc[1] === s_off.bal.errnb_grc[1]
        @test s_on.bal.errcb_col[1] === s_off.bal.errcb_col[1]
        @test s_on.bal.errnb_col[1] === s_off.bal.errnb_col[1]
    end
end
