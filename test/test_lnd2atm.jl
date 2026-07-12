@testset "Lnd2AtmData" begin

    @testset "Lnd2AtmParamsData default construction" begin
        p = CLM.Lnd2AtmParamsData()
        @test p.melt_non_icesheet_ice_runoff == false
    end

    @testset "Lnd2AtmData default construction" begin
        l = CLM.Lnd2AtmData()
        @test length(l.t_rad_grc) == 0
        @test length(l.eflx_sh_ice_to_liq_col) == 0
        @test size(l.albd_grc) == (0, 0)
        @test size(l.albi_grc) == (0, 0)
        @test size(l.flxdst_grc) == (0, 0)
        @test size(l.ddvel_grc) == (0, 0)
        @test size(l.flxvoc_grc) == (0, 0)
        @test size(l.fireflx_grc) == (0, 0)
        @test length(l.fireztop_grc) == 0
    end

    @testset "lnd2atm_params_init!" begin
        p = CLM.Lnd2AtmParamsData()
        CLM.lnd2atm_params_init!(p; melt_non_icesheet_ice_runoff = true)
        @test p.melt_non_icesheet_ice_runoff == true

        CLM.lnd2atm_params_init!(p; melt_non_icesheet_ice_runoff = false)
        @test p.melt_non_icesheet_ice_runoff == false
    end

    @testset "lnd2atm_init! basic allocation" begin
        ng, nc = 4, 8
        l = CLM.Lnd2AtmData()
        CLM.lnd2atm_init!(l, ng, nc)

        # gridcell-level 1D fields
        @test length(l.t_rad_grc) == ng
        @test all(x -> x == 0.0, l.t_rad_grc)
        @test length(l.t_ref2m_grc) == ng
        @test length(l.u_ref10m_grc) == ng
        @test length(l.taux_grc) == ng
        @test length(l.tauy_grc) == ng
        @test length(l.eflx_lwrad_out_grc) == ng
        @test length(l.eflx_sh_tot_grc) == ng
        @test all(x -> x == 0.0, l.eflx_sh_tot_grc)
        @test length(l.eflx_sh_precip_conversion_grc) == ng
        @test length(l.eflx_lh_tot_grc) == ng
        @test length(l.fsa_grc) == ng
        @test length(l.z0m_grc) == ng
        @test length(l.net_carbon_exchange_grc) == ng
        @test length(l.nem_grc) == ng
        @test length(l.ram1_grc) == ng
        @test length(l.fv_grc) == ng
        @test length(l.ch4_surf_flux_tot_grc) == ng

        # gridcell-level 2D fields
        @test size(l.albd_grc) == (ng, CLM.NUMRAD)
        @test all(x -> x == 0.0, l.albd_grc)
        @test size(l.albi_grc) == (ng, CLM.NUMRAD)
        @test size(l.flxdst_grc) == (ng, CLM.NDST)

        # column-level
        @test length(l.eflx_sh_ice_to_liq_col) == nc
        @test all(x -> x == 0.0, l.eflx_sh_ice_to_liq_col)

        # conditional fields not allocated by default
        @test size(l.flxvoc_grc) == (0, 0)
        @test size(l.fireflx_grc) == (0, 0)
        @test length(l.fireztop_grc) == 0
        @test size(l.ddvel_grc) == (0, 0)
    end

    @testset "lnd2atm_init! with conditional allocations" begin
        ng, nc = 3, 6
        l = CLM.Lnd2AtmData()
        CLM.lnd2atm_init!(l, ng, nc;
            n_megan_comps = 5,
            n_fire_emis_comps = 3,
            n_drydep = 2)

        @test size(l.flxvoc_grc) == (ng, 5)
        @test all(x -> x == 0.0, l.flxvoc_grc)
        @test size(l.fireflx_grc) == (ng, 3)
        @test all(x -> x == 0.0, l.fireflx_grc)
        @test length(l.fireztop_grc) == ng
        @test all(x -> x == 0.0, l.fireztop_grc)
        @test size(l.ddvel_grc) == (ng, 2)
        @test all(x -> x == 0.0, l.ddvel_grc)
    end

    @testset "lnd2atm_read_namelist!" begin
        l = CLM.Lnd2AtmData()
        CLM.lnd2atm_init!(l, 3, 5)
        CLM.lnd2atm_read_namelist!(l; melt_non_icesheet_ice_runoff = true)
        @test l.params.melt_non_icesheet_ice_runoff == true

        CLM.lnd2atm_read_namelist!(l; melt_non_icesheet_ice_runoff = false)
        @test l.params.melt_non_icesheet_ice_runoff == false
    end

    @testset "lnd2atm_init_for_testing! default params" begin
        ng, nc = 3, 5
        l = CLM.Lnd2AtmData()
        CLM.lnd2atm_init_for_testing!(l, ng, nc)

        @test l.params.melt_non_icesheet_ice_runoff == false
        @test length(l.t_rad_grc) == ng
        @test length(l.eflx_sh_ice_to_liq_col) == nc
    end

    @testset "lnd2atm_init_for_testing! with custom params" begin
        ng, nc = 2, 4
        custom_params = CLM.Lnd2AtmParamsData(
            melt_non_icesheet_ice_runoff = true)

        l = CLM.Lnd2AtmData()
        CLM.lnd2atm_init_for_testing!(l, ng, nc; params = custom_params)

        @test l.params.melt_non_icesheet_ice_runoff == true
        @test length(l.t_rad_grc) == ng
    end

    @testset "lnd2atm_clean!" begin
        ng, nc = 3, 5
        l = CLM.Lnd2AtmData()
        CLM.lnd2atm_init!(l, ng, nc;
            n_megan_comps = 2, n_fire_emis_comps = 1, n_drydep = 3)

        @test length(l.t_rad_grc) == ng
        @test size(l.flxvoc_grc) == (ng, 2)

        CLM.lnd2atm_clean!(l)

        @test length(l.t_rad_grc) == 0
        @test length(l.t_ref2m_grc) == 0
        @test length(l.u_ref10m_grc) == 0
        @test size(l.albd_grc) == (0, 0)
        @test size(l.albi_grc) == (0, 0)
        @test length(l.taux_grc) == 0
        @test length(l.tauy_grc) == 0
        @test length(l.eflx_lh_tot_grc) == 0
        @test length(l.eflx_sh_tot_grc) == 0
        @test length(l.eflx_sh_precip_conversion_grc) == 0
        @test length(l.eflx_lwrad_out_grc) == 0
        @test length(l.fsa_grc) == 0
        @test length(l.z0m_grc) == 0
        @test length(l.net_carbon_exchange_grc) == 0
        @test length(l.nem_grc) == 0
        @test length(l.ram1_grc) == 0
        @test length(l.fv_grc) == 0
        @test size(l.flxdst_grc) == (0, 0)
        @test size(l.ddvel_grc) == (0, 0)
        @test size(l.flxvoc_grc) == (0, 0)
        @test size(l.fireflx_grc) == (0, 0)
        @test length(l.fireztop_grc) == 0
        @test length(l.ch4_surf_flux_tot_grc) == 0
        @test length(l.eflx_sh_ice_to_liq_col) == 0
    end

    @testset "stub functions run without error" begin
        l = CLM.Lnd2AtmData()
        CLM.lnd2atm_init!(l, 3, 5)
        @test CLM.lnd2atm_init_history!(l, 3, 5) === nothing
    end

    @testset "field mutability" begin
        l = CLM.Lnd2AtmData()
        CLM.lnd2atm_init!(l, 3, 5)

        l.t_rad_grc[1] = 300.0
        @test l.t_rad_grc[1] == 300.0

        l.eflx_sh_ice_to_liq_col[2] = 42.0
        @test l.eflx_sh_ice_to_liq_col[2] == 42.0

        l.albd_grc[1, 1] = 0.5
        @test l.albd_grc[1, 1] == 0.5

        l.flxdst_grc[2, 3] = 1.0e-6
        @test l.flxdst_grc[2, 3] == 1.0e-6
    end

    @testset "re-init overwrites previous state" begin
        l = CLM.Lnd2AtmData()
        CLM.lnd2atm_init!(l, 3, 5)
        l.t_rad_grc[1] = 999.0

        CLM.lnd2atm_init!(l, 10, 20)
        @test length(l.t_rad_grc) == 10
        @test all(x -> x == 0.0, l.t_rad_grc)
        @test length(l.eflx_sh_ice_to_liq_col) == 20
    end

    # ------------------------------------------------------------------
    # handle_ice_runoff! / add_liq_from_ice_to_runoff! (lnd2atmMod.F90)
    #
    # qflx_ice_runoff_col = qflx_ice_runoff_snwcp_col + qflx_ice_runoff_xs_col is the
    # ice-runoff sink the column water balance subtracts (BalanceCheckMod.F90:650).
    # ------------------------------------------------------------------
    function make_ice_runoff_data(; nc = 3, ng = 2, nl = 2)
        CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)

        l2a = CLM.Lnd2AtmData()
        CLM.lnd2atm_init!(l2a, ng, nc)

        col = CLM.ColumnData()
        col.gridcell = [1, 1, 2]
        col.landunit = [1, 2, 2]     # col 1 on a soil landunit, cols 2-3 on glacier
        col.active   = fill(true, nc)

        lun = CLM.LandunitData()
        lun.itype = [CLM.ISTSOIL, CLM.ISTICE]

        water = CLM.WaterData()
        CLM.water_init!(water, nc, nc, nl, ng)
        wfb = water.waterfluxbulk_inst
        wfb.wf.qflx_ice_runoff_snwcp_col .= [1.0e-3, 2.0e-3, 0.0]
        wfb.wf.qflx_ice_runoff_xs_col    .= [4.0e-4, 0.0, 0.0]
        wfb.wf.qflx_qrgwl_col            .= 0.0
        wfb.wf.qflx_runoff_col           .= 0.0

        return (l2a=l2a, col=col, lun=lun, wfb=wfb, nc=nc, ng=ng, bounds=1:nc)
    end

    @testset "handle_ice_runoff! sums snwcp + excess-ice runoff" begin
        d = make_ice_runoff_data()
        CLM.handle_ice_runoff!(d.l2a, d.wfb, d.col, d.lun, d.bounds)

        # Fortran: qflx_ice_runoff_col = qflx_ice_runoff_snwcp_col + qflx_ice_runoff_xs_col
        @test d.l2a.qflx_ice_runoff_col ≈ [1.4e-3, 2.0e-3, 0.0]
        # Non-zero where snow capping / excess ice produced solid runoff — this is the
        # term the driver must pass to balance_check! (it used to pass zeros).
        @test d.l2a.qflx_ice_runoff_col[1] > 0.0
        @test d.l2a.qflx_ice_runoff_col[2] > 0.0

        # melt_non_icesheet_ice_runoff defaults to false (Fortran namelist default)
        @test d.l2a.params.melt_non_icesheet_ice_runoff == false
        @test all(iszero, d.l2a.qflx_liq_from_ice_col)
        @test all(iszero, d.l2a.eflx_sh_ice_to_liq_col)
    end

    @testset "handle_ice_runoff! melt_non_icesheet_ice_runoff conversion" begin
        d = make_ice_runoff_data()
        CLM.lnd2atm_params_init!(d.l2a.params; melt_non_icesheet_ice_runoff = true)

        # No gridcell has its icesheet runoff melted en route (glc_behavior default).
        CLM.handle_ice_runoff!(d.l2a, d.wfb, d.col, d.lun, d.bounds)

        # col 1 (non-icesheet landunit) converts: ice runoff -> liquid runoff, and the
        # ice->liquid phase change absorbs energy (negative sensible heat to atm).
        @test d.l2a.qflx_ice_runoff_col[1] == 0.0
        @test d.l2a.qflx_liq_from_ice_col[1] ≈ 1.4e-3
        @test d.l2a.eflx_sh_ice_to_liq_col[1] ≈ -1.4e-3 * CLM.HFUS

        # col 2 (istice, gridcell not flagged) keeps its ice runoff.
        @test d.l2a.qflx_ice_runoff_col[2] ≈ 2.0e-3
        @test d.l2a.qflx_liq_from_ice_col[2] == 0.0
        @test d.l2a.eflx_sh_ice_to_liq_col[2] == 0.0

        # Flag the icesheet gridcell (col 2 is in gridcell 1) -> it converts too.
        d2 = make_ice_runoff_data()
        CLM.lnd2atm_params_init!(d2.l2a.params; melt_non_icesheet_ice_runoff = true)
        CLM.handle_ice_runoff!(d2.l2a, d2.wfb, d2.col, d2.lun, d2.bounds;
                               ice_runoff_melted_grc = [true, true])
        @test d2.l2a.qflx_ice_runoff_col[2] == 0.0
        @test d2.l2a.qflx_liq_from_ice_col[2] ≈ 2.0e-3
    end

    @testset "add_liq_from_ice_to_runoff! routes melted ice into qrgwl" begin
        d = make_ice_runoff_data()
        CLM.lnd2atm_params_init!(d.l2a.params; melt_non_icesheet_ice_runoff = true)
        CLM.handle_ice_runoff!(d.l2a, d.wfb, d.col, d.lun, d.bounds)
        CLM.add_liq_from_ice_to_runoff!(d.wfb, d.l2a, d.col, d.bounds)

        # The melted ice must leave via qflx_qrgwl (which the balance check accounts
        # for), otherwise the converted water would vanish from the balance.
        @test d.wfb.wf.qflx_qrgwl_col[1] ≈ 1.4e-3
        @test d.wfb.wf.qflx_runoff_col[1] ≈ 1.4e-3
        @test d.wfb.wf.qflx_qrgwl_col[2] == 0.0   # istice: still ice runoff, not liquid

        # Default (no conversion) => nothing added.
        d0 = make_ice_runoff_data()
        CLM.handle_ice_runoff!(d0.l2a, d0.wfb, d0.col, d0.lun, d0.bounds)
        CLM.add_liq_from_ice_to_runoff!(d0.wfb, d0.l2a, d0.col, d0.bounds)
        @test all(iszero, d0.wfb.wf.qflx_qrgwl_col)
        @test all(iszero, d0.wfb.wf.qflx_runoff_col)
    end

end
