@testset "Instance factory (clm_instMod)" begin
    ng = 3   # gridcells
    nl = 5   # landunits
    nc = 10  # columns
    np = 20  # patches

    nlevdecomp_full = 10
    ndecomp_pools = 7
    ndecomp_cascade_transitions = 5

    @testset "Construction and init" begin
        inst = CLM.CLMInstances()
        CLM.clm_instInit!(inst; ng=ng, nl=nl, nc=nc, np=np,
                          nlevdecomp_full=nlevdecomp_full,
                          ndecomp_pools=ndecomp_pools,
                          ndecomp_cascade_transitions=ndecomp_cascade_transitions)

        # --- Grid hierarchy dimensions ---
        @test length(inst.gridcell.latdeg) == ng
        @test length(inst.landunit.itype) == nl
        @test length(inst.column.gridcell) == nc
        @test length(inst.patch.gridcell) == np

        # --- Topography ---
        @test length(inst.topo.topo_col) == nc

        # --- Urban ---
        @test length(inst.urbanparams.em_roof) == nl

        # --- Core physics ---
        @test size(inst.temperature.t_soisno_col, 1) == nc
        @test size(inst.energyflux.eflx_sh_tot_patch, 1) == np
        @test length(inst.canopystate.laisun_patch) == np
        @test size(inst.soilstate.watsat_col, 1) == nc
        @test length(inst.soilhydrology.zwt_col) == nc
        @test length(inst.frictionvel.u10_patch) == np
        @test length(inst.lakestate.lake_icefrac_col) > 0
        @test size(inst.solarabs.sabg_patch, 1) == np
        @test size(inst.surfalb.albgrd_col, 1) == nc

        # --- Water ---
        @test length(inst.water.waterfluxbulk_inst.qflx_ev_snow_patch) == np

        # --- Coupling ---
        @test length(inst.atm2lnd.forc_u_grc) == ng
        @test length(inst.lnd2atm.eflx_lh_tot_grc) == ng

        # --- Crop ---
        @test length(inst.crop.nyrs_crop_active_patch) == np

        # --- Soil biogeochemistry ---
        @test size(inst.soilbiogeochem_carbonstate.decomp_cpools_vr_col, 1) == nc
        @test size(inst.soilbiogeochem_carbonflux.decomp_cpools_sourcesink_col, 1) == nc
        @test size(inst.soilbiogeochem_nitrogenstate.decomp_npools_vr_col, 1) == nc
        @test size(inst.soilbiogeochem_nitrogenflux.decomp_npools_sourcesink_col, 1) == nc
        @test size(inst.soilbiogeochem_state.fpi_vr_col, 1) == nc

        # --- CN vegetation ---
        @test length(inst.bgc_vegetation.cnveg_state_inst.dormant_flag_patch) == np
        @test length(inst.cn_products.cropprod1_grc) == ng
    end

    @testset "Multiple independent instances" begin
        inst1 = CLM.CLMInstances()
        inst2 = CLM.CLMInstances()
        CLM.clm_instInit!(inst1; ng=2, nl=3, nc=5, np=8)
        CLM.clm_instInit!(inst2; ng=4, nl=6, nc=12, np=24)

        # Instances should be independent
        @test length(inst1.temperature.t_grnd_col) == 5
        @test length(inst2.temperature.t_grnd_col) == 12
        @test length(inst1.patch.gridcell) == 8
        @test length(inst2.patch.gridcell) == 24
    end
end
