@testset "WaterFluxData" begin

    # Ensure varpar is initialized for dimension parameters
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)

    nlevsno = CLM.varpar.nlevsno
    nlevsoi = CLM.varpar.nlevsoi

    @testset "default construction" begin
        wf = CLM.WaterFluxData()

        # Patch fields
        @test length(wf.qflx_through_snow_patch) == 0
        @test length(wf.qflx_through_liq_patch) == 0
        @test length(wf.qflx_intercepted_snow_patch) == 0
        @test length(wf.qflx_intercepted_liq_patch) == 0
        @test length(wf.qflx_snocanfall_patch) == 0
        @test length(wf.qflx_liqcanfall_patch) == 0
        @test length(wf.qflx_snow_unload_patch) == 0
        @test length(wf.qflx_tran_veg_patch) == 0
        @test length(wf.qflx_evap_tot_patch) == 0
        @test length(wf.qflx_irrig_drip_patch) == 0
        @test length(wf.qflx_irrig_sprinkler_patch) == 0

        # Column fields
        @test length(wf.qflx_liq_grnd_col) == 0
        @test length(wf.qflx_snow_grnd_col) == 0
        @test length(wf.qflx_infl_col) == 0
        @test length(wf.qflx_surf_col) == 0
        @test length(wf.qflx_drain_col) == 0
        @test length(wf.qflx_runoff_col) == 0
        @test length(wf.qflx_snomelt_col) == 0
        @test length(wf.qflx_snow_drain_col) == 0
        @test length(wf.qflx_sfc_irrig_col) == 0
        @test length(wf.qflx_gw_uncon_irrig_col) == 0
        @test length(wf.qflx_gw_con_irrig_col) == 0

        # Column 2D
        @test size(wf.qflx_snofrz_lyr_col) == (0, 0)
        @test size(wf.qflx_snow_percolation_col) == (0, 0)
        @test size(wf.qflx_gw_uncon_irrig_lyr_col) == (0, 0)

        # Landunit
        @test length(wf.volumetric_streamflow_lun) == 0

        # Gridcell
        @test length(wf.qflx_liq_dynbal_grc) == 0
        @test length(wf.qflx_ice_dynbal_grc) == 0
    end

    @testset "waterflux_init!" begin
        nc, np, nl, ng = 10, 15, 3, 2
        wf = CLM.WaterFluxData()
        CLM.waterflux_init!(wf, nc, np, nl, ng)

        # Patch 1D sizes
        @test length(wf.qflx_through_snow_patch) == np
        @test length(wf.qflx_through_liq_patch) == np
        @test length(wf.qflx_intercepted_snow_patch) == np
        @test length(wf.qflx_intercepted_liq_patch) == np
        @test length(wf.qflx_snocanfall_patch) == np
        @test length(wf.qflx_liqcanfall_patch) == np
        @test length(wf.qflx_snow_unload_patch) == np
        @test length(wf.qflx_solidevap_from_top_layer_patch) == np
        @test length(wf.qflx_tran_veg_patch) == np
        @test length(wf.qflx_liqdew_to_top_layer_patch) == np
        @test length(wf.qflx_soliddew_to_top_layer_patch) == np
        @test length(wf.qflx_evap_veg_patch) == np
        @test length(wf.qflx_evap_can_patch) == np
        @test length(wf.qflx_evap_soi_patch) == np
        @test length(wf.qflx_evap_tot_patch) == np
        @test length(wf.qflx_liqevap_from_top_layer_patch) == np
        @test length(wf.qflx_irrig_drip_patch) == np
        @test length(wf.qflx_irrig_sprinkler_patch) == np

        # Column 1D sizes
        @test length(wf.qflx_liq_grnd_col) == nc
        @test length(wf.qflx_snow_grnd_col) == nc
        @test length(wf.qflx_rain_plus_snomelt_col) == nc
        @test length(wf.qflx_infl_col) == nc
        @test length(wf.qflx_surf_col) == nc
        @test length(wf.qflx_drain_col) == nc
        @test length(wf.qflx_runoff_col) == nc
        @test length(wf.qflx_snomelt_col) == nc
        @test length(wf.qflx_snow_drain_col) == nc
        @test length(wf.qflx_sfc_irrig_col) == nc
        @test length(wf.qflx_gw_uncon_irrig_col) == nc
        @test length(wf.qflx_gw_con_irrig_col) == nc
        @test length(wf.qflx_ice_runoff_snwcp_col) == nc
        @test length(wf.qflx_ice_runoff_xs_col) == nc
        @test length(wf.qflx_glcice_dyn_water_flux_col) == nc

        # Column 2D sizes
        @test size(wf.qflx_snofrz_lyr_col) == (nc, nlevsno)
        @test size(wf.qflx_snow_percolation_col) == (nc, nlevsno)
        @test size(wf.qflx_gw_uncon_irrig_lyr_col) == (nc, nlevsoi)

        # Landunit 1D
        @test length(wf.volumetric_streamflow_lun) == nl

        # Gridcell 1D
        @test length(wf.qflx_liq_dynbal_grc) == ng
        @test length(wf.qflx_ice_dynbal_grc) == ng

        # qflx_ice_runoff_xs_col must be ZERO-initialized, not NaN: Fortran
        # (WaterFluxType.F90:897, InitCold) zeroes it because it "only gets set in the
        # hydrology filter" — columns outside that filter (lakes) would otherwise carry
        # NaN into qflx_ice_runoff_col (= snwcp + xs) and poison their errh2o_col.
        @test all(iszero, wf.qflx_ice_runoff_xs_col)

        # NaN initialization for standard fields
        @test all(isnan, wf.qflx_through_snow_patch)
        @test all(isnan, wf.qflx_through_liq_patch)
        @test all(isnan, wf.qflx_liq_grnd_col)
        @test all(isnan, wf.qflx_snow_grnd_col)
        @test all(wf.qflx_infl_col .== 0.0)
        @test all(wf.qflx_surf_col .== 0.0)
        @test all(wf.qflx_drain_col .== 0.0)
        @test all(wf.qflx_runoff_col .== 0.0)
        @test all(isnan, wf.qflx_snofrz_lyr_col)
        @test all(isnan, wf.qflx_snow_percolation_col)
        @test all(isnan, wf.qflx_gw_uncon_irrig_lyr_col)
        @test all(isnan, wf.qflx_liq_dynbal_grc)
        @test all(isnan, wf.qflx_ice_dynbal_grc)
        @test all(isnan, wf.volumetric_streamflow_lun)

        # Zero initialization for special fields (ival=0 in Fortran)
        @test all(wf.qflx_solidevap_from_top_layer_patch .== 0.0)
        @test all(wf.qflx_solidevap_from_top_layer_col .== 0.0)
    end

    @testset "waterflux_clean!" begin
        wf = CLM.WaterFluxData()
        CLM.waterflux_init!(wf, 5, 8, 2, 1)
        CLM.waterflux_clean!(wf)

        # Patch 1D
        @test length(wf.qflx_through_snow_patch) == 0
        @test length(wf.qflx_tran_veg_patch) == 0
        @test length(wf.qflx_irrig_drip_patch) == 0

        # Column 1D
        @test length(wf.qflx_liq_grnd_col) == 0
        @test length(wf.qflx_infl_col) == 0
        @test length(wf.qflx_runoff_col) == 0
        @test length(wf.qflx_sfc_irrig_col) == 0

        # Column 2D
        @test size(wf.qflx_snofrz_lyr_col) == (0, 0)
        @test size(wf.qflx_snow_percolation_col) == (0, 0)
        @test size(wf.qflx_gw_uncon_irrig_lyr_col) == (0, 0)

        # Landunit
        @test length(wf.volumetric_streamflow_lun) == 0

        # Gridcell
        @test length(wf.qflx_liq_dynbal_grc) == 0
        @test length(wf.qflx_ice_dynbal_grc) == 0
    end

    @testset "waterflux_init_cold! (soil columns)" begin
        nc, np, nl, ng = 3, 4, 1, 1
        wf = CLM.WaterFluxData()
        CLM.waterflux_init!(wf, nc, np, nl, ng)

        landunit_col = fill(1, nc)
        lun_itype = [CLM.ISTSOIL]

        CLM.waterflux_init_cold!(wf, 1:nc, 1:np, 1:nl;
            landunit_col=landunit_col,
            lun_itype=lun_itype,
            use_hillslope_routing=false)

        # Patch zeros
        for p in 1:np
            @test wf.qflx_snocanfall_patch[p] == 0.0
            @test wf.qflx_liqcanfall_patch[p] == 0.0
            @test wf.qflx_snow_unload_patch[p] == 0.0
            @test wf.qflx_liqevap_from_top_layer_patch[p] == 0.0
            @test wf.qflx_liqdew_to_top_layer_patch[p] == 0.0
            @test wf.qflx_soliddew_to_top_layer_patch[p] == 0.0
            @test wf.qflx_irrig_drip_patch[p] == 0.0
            @test wf.qflx_irrig_sprinkler_patch[p] == 0.0
            @test wf.qflx_tran_veg_patch[p] == 0.0
            @test wf.qflx_evap_veg_patch[p] == 0.0
        end

        # Column zeros
        for c in 1:nc
            @test wf.qflx_sfc_irrig_col[c] == 0.0
            @test wf.qflx_gw_uncon_irrig_col[c] == 0.0
            @test wf.qflx_gw_con_irrig_col[c] == 0.0
            @test wf.qflx_liqevap_from_top_layer_col[c] == 0.0
            @test wf.qflx_liqdew_to_top_layer_col[c] == 0.0
            @test wf.qflx_soliddew_to_top_layer_col[c] == 0.0
            @test wf.qflx_snow_drain_col[c] == 0.0
            @test wf.qflx_ice_runoff_xs_col[c] == 0.0
            @test wf.qflx_glcice_dyn_water_flux_col[c] == 0.0
        end

        # Soil/crop columns should have drain/surf/latflow zeroed
        for c in 1:nc
            @test wf.qflx_drain_col[c] == 0.0
            @test wf.qflx_surf_col[c] == 0.0
            @test wf.qflx_latflow_in_col[c] == 0.0
            @test wf.qflx_latflow_out_col[c] == 0.0
            @test wf.volumetric_discharge_col[c] == 0.0
        end

        # Irrigation layer columns should be zeroed
        for c in 1:nc
            for j in 1:nlevsoi
                @test wf.qflx_gw_uncon_irrig_lyr_col[c, j] == 0.0
            end
        end
    end

    @testset "waterflux_init_cold! (with hillslope routing)" begin
        nc, np, nl, ng = 2, 3, 2, 1
        wf = CLM.WaterFluxData()
        CLM.waterflux_init!(wf, nc, np, nl, ng)

        landunit_col = [1, 2]
        lun_itype = [CLM.ISTSOIL, CLM.ISTCROP]

        CLM.waterflux_init_cold!(wf, 1:nc, 1:np, 1:nl;
            landunit_col=landunit_col,
            lun_itype=lun_itype,
            use_hillslope_routing=true)

        # Both landunits are soil/crop → streamflow should be zeroed
        @test wf.volumetric_streamflow_lun[1] == 0.0
        @test wf.volumetric_streamflow_lun[2] == 0.0
    end

    @testset "waterflux_init_cold! (non-soil column skips drain/surf)" begin
        nc, np, nl, ng = 2, 2, 2, 1
        wf = CLM.WaterFluxData()
        CLM.waterflux_init!(wf, nc, np, nl, ng)

        # Column 1 → landunit 1 (ISTSOIL), column 2 → landunit 2 (ISTWET)
        landunit_col = [1, 2]
        lun_itype = [CLM.ISTSOIL, CLM.ISTWET]

        CLM.waterflux_init_cold!(wf, 1:nc, 1:np, 1:nl;
            landunit_col=landunit_col,
            lun_itype=lun_itype)

        # Soil column should have drain zeroed
        @test wf.qflx_drain_col[1] == 0.0
        @test wf.qflx_surf_col[1] == 0.0

        # Non-soil column drain/surf should still be 0.0 (zero-initialized)
        @test wf.qflx_drain_col[2] == 0.0
        @test wf.qflx_surf_col[2] == 0.0
    end

    @testset "stub functions run without error" begin
        wf = CLM.WaterFluxData()
        CLM.waterflux_init!(wf, 5, 8, 2, 1)

        @test CLM.waterflux_init_history!(wf, 1:5) === nothing
        @test CLM.waterflux_restart!(wf, 1:5) === nothing
    end

    @testset "field mutability" begin
        wf = CLM.WaterFluxData()
        CLM.waterflux_init!(wf, 3, 4, 2, 1)

        wf.qflx_through_snow_patch[1] = 0.5
        @test wf.qflx_through_snow_patch[1] == 0.5

        wf.qflx_infl_col[2] = 1.5
        @test wf.qflx_infl_col[2] == 1.5

        wf.qflx_snofrz_lyr_col[1, 1] = 0.01
        @test wf.qflx_snofrz_lyr_col[1, 1] == 0.01

        wf.qflx_gw_uncon_irrig_lyr_col[2, 3] = 0.02
        @test wf.qflx_gw_uncon_irrig_lyr_col[2, 3] == 0.02

        wf.volumetric_streamflow_lun[1] = 99.0
        @test wf.volumetric_streamflow_lun[1] == 99.0

        wf.qflx_liq_dynbal_grc[1] = 0.003
        @test wf.qflx_liq_dynbal_grc[1] == 0.003
    end

    @testset "re-init overwrites previous state" begin
        wf = CLM.WaterFluxData()
        CLM.waterflux_init!(wf, 3, 4, 2, 1)
        wf.qflx_infl_col[1] = 999.0

        CLM.waterflux_init!(wf, 7, 10, 3, 2)
        @test length(wf.qflx_infl_col) == 7
        @test all(wf.qflx_infl_col .== 0.0)
        @test length(wf.qflx_through_snow_patch) == 10
        @test size(wf.qflx_snofrz_lyr_col) == (7, nlevsno)
        @test size(wf.qflx_gw_uncon_irrig_lyr_col) == (7, nlevsoi)
        @test length(wf.volumetric_streamflow_lun) == 3
        @test length(wf.qflx_liq_dynbal_grc) == 2
    end

    # =====================================================================
    # Regression guard for the DEAD-InitCold class of bug.
    #
    # `waterflux_init_cold!` (waterflux_type%InitCold, WaterFluxType.F90:860) and
    # `waterfluxbulk_init_cold!` (InitBulkCold, WaterFluxBulkType.F90:255) were both
    # DEAD in the port: defined, unit-tested above, and never called from any run
    # path. Fortran zeroes these fields because several are ONLY written inside the
    # hydrology filter, so columns outside it (lake, glacier, urban roof/wall,
    # wetland) keep their allocation-time NaN forever. `qflx_ice_runoff_xs_col` was
    # one instance of this class (PR #207); the audit found many more, e.g.
    # `qflx_gw_uncon_irrig_col` NaN on every column, `qflx_tran_veg_patch` /
    # `qflx_evap_veg_patch` NaN on lake + urban patches, `qflx_phs_neg_col` NaN
    # everywhere with PHS off.
    #
    # Both routines are now wired into `cold_start_initialize!` via
    # `init_water_flux_cold!` (src/infrastructure/cold_start.jl). This guard asserts
    # that EVERY field Fortran's InitCold/InitBulkCold zeroes really is non-NaN on
    # EVERY column/patch — including a non-soil (lake) landunit that never runs the
    # hydrology filter. It fails if either routine is ever unwired again.
    # =====================================================================
    @testset "InitCold covers every Fortran-zeroed field on non-hydrology landunits" begin
        # Column 1 = ISTSOIL (runs the hydrology filter),
        # Column 2 = ISTDLAK (deep lake — outside it; this is the poisoned case).
        nc, np, nl, ng = 2, 2, 2, 1
        landunit_col = [1, 2]
        lun_itype    = [CLM.ISTSOIL, CLM.ISTDLAK]

        wfb = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(wfb, nc, np, nl, ng)
        wf = wfb.wf

        CLM.waterflux_init_cold!(wf, 1:nc, 1:np, 1:nl;
            landunit_col=landunit_col, lun_itype=lun_itype,
            use_hillslope_routing=false)
        CLM.waterfluxbulk_init_cold!(wfb, 1:nc, 1:np)

        # --- Fields Fortran's InitCold zeroes for ALL columns (WaterFluxType.F90:883-903) ---
        all_col = [:qflx_sfc_irrig_col, :qflx_gw_uncon_irrig_col, :qflx_gw_con_irrig_col,
                   :qflx_liqevap_from_top_layer_col, :qflx_liqdew_to_top_layer_col,
                   :qflx_soliddew_to_top_layer_col, :qflx_snow_drain_col,
                   :qflx_ice_runoff_xs_col, :qflx_glcice_dyn_water_flux_col]
        # --- ...and for ALL patches (WaterFluxType.F90:875-888, 906-908) ---
        all_pat = [:qflx_snocanfall_patch, :qflx_liqcanfall_patch, :qflx_snow_unload_patch,
                   :qflx_liqevap_from_top_layer_patch, :qflx_liqdew_to_top_layer_patch,
                   :qflx_soliddew_to_top_layer_patch, :qflx_irrig_drip_patch,
                   :qflx_irrig_sprinkler_patch, :qflx_tran_veg_patch, :qflx_evap_veg_patch]
        # --- InitBulkCold (WaterFluxBulkType.F90:264-268) ---
        bulk_col = [:qflx_phs_neg_col, :qflx_h2osfc_surf_col]
        bulk_pat = [:qflx_snowindunload_patch, :qflx_snotempunload_patch]

        for f in all_col, c in 1:nc
            v = getfield(wf, f)[c]
            @test !isnan(v)
            @test v == 0.0
        end
        for f in all_pat, p in 1:np
            v = getfield(wf, f)[p]
            @test !isnan(v)
            @test v == 0.0
        end
        for f in bulk_col, c in 1:nc
            v = getfield(wfb, f)[c]
            @test !isnan(v)
            @test v == 0.0
        end
        for f in bulk_pat, p in 1:np
            v = getfield(wfb, f)[p]
            @test !isnan(v)
            @test v == 0.0
        end
        # 2D irrigation-by-layer, all columns (WaterFluxType.F90:885)
        for c in 1:nc, j in 1:nlevsoi
            @test !isnan(wf.qflx_gw_uncon_irrig_lyr_col[c, j])
            @test wf.qflx_gw_uncon_irrig_lyr_col[c, j] == 0.0
        end

        # Fortran zeroes drain/surf/latflow/discharge ONLY on soil/crop — the lake
        # column is intentionally left alone here (it is not in the CNNLeaching loop).
        # Assert the soil column got them, i.e. the branch really ran.
        @test wf.qflx_drain_col[1] == 0.0
        @test wf.qflx_surf_col[1] == 0.0
        @test wf.qflx_latflow_in_col[1] == 0.0
        @test wf.qflx_latflow_out_col[1] == 0.0
        @test wf.volumetric_discharge_col[1] == 0.0
    end

    @testset "init_water_flux_cold! is wired into the cold-start path" begin
        # Guards the wiring itself: cold_start_initialize! must call InitCold. If the
        # call is ever dropped, `init_water_flux_cold!` vanishing (or no longer being
        # invoked from cold_start_initialize!) is what re-opens the whole bug class.
        @test isdefined(CLM, :init_water_flux_cold!)
        src = read(joinpath(@__DIR__, "..", "src", "infrastructure", "cold_start.jl"), String)
        body = match(r"function cold_start_initialize!.*?\n(.*?)\nend\n"s, src)
        @test body !== nothing
        @test occursin("init_water_flux_cold!", body.captures[1])
    end
end
