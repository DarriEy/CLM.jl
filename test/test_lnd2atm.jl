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

end
