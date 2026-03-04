# Include the source file into the CLM module (will be added to CLM.jl later)
if !isdefined(CLM, :heat_base_temp)
    Base.include(CLM, joinpath(@__DIR__, "..", "src", "biogeophys", "total_water_heat.jl"))
end

@testset "TotalWaterAndHeatMod" begin

    # Ensure varpar is initialized for dimension parameters
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
    CLM.varcon_init!()

    nlevsno        = CLM.varpar.nlevsno
    nlevgrnd       = CLM.varpar.nlevgrnd
    nlevurb        = CLM.varpar.nlevurb
    nlevmaxurbgrnd = CLM.varpar.nlevmaxurbgrnd
    nlevlak        = CLM.varpar.nlevlak
    nlevtot        = nlevsno + nlevmaxurbgrnd

    # -----------------------------------------------------------------------
    # Helper: set up minimal column/landunit/water state data for nc columns
    # -----------------------------------------------------------------------
    function setup_basic(nc; col_types=nothing, lun_types=nothing,
                         hydro_active=nothing, snl_vals=nothing)
        col = CLM.ColumnData()
        CLM.column_init!(col, nc)

        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, nc)  # one landunit per column for simplicity

        # Link columns to landunits (1:1 mapping)
        for c in 1:nc
            col.landunit[c] = c
            lun.itype[c] = lun_types !== nothing ? lun_types[c] : CLM.ISTSOIL
            lun.urbpoi[c] = (CLM.ISTURB_MIN <= lun.itype[c] <= CLM.ISTURB_MAX)
            lun.lakpoi[c] = (lun.itype[c] == CLM.ISTDLAK)
        end

        if col_types !== nothing
            for c in 1:nc
                col.itype[c] = col_types[c]
            end
        else
            for c in 1:nc
                col.itype[c] = CLM.ISTSOIL
            end
        end

        if hydro_active !== nothing
            for c in 1:nc
                col.hydrologically_active[c] = hydro_active[c]
            end
        else
            for c in 1:nc
                col.hydrologically_active[c] = true
            end
        end

        if snl_vals !== nothing
            for c in 1:nc
                col.snl[c] = snl_vals[c]
            end
        else
            for c in 1:nc
                col.snl[c] = 0
            end
        end

        # Initialize layer thicknesses to 0.1m
        for c in 1:nc
            for j in 1:nlevtot
                col.dz[c, j] = 0.1
            end
            for j in 1:nlevlak
                col.dz_lake[c, j] = 1.0
            end
        end

        # Water state
        ws = CLM.WaterStateData()
        CLM.waterstate_init!(ws, nc, nc, nc, nc)

        for c in 1:nc
            ws.h2osno_no_layers_col[c] = 0.0
            ws.h2osfc_col[c] = 0.0
            ws.wa_col[c] = CLM.AQUIFER_WATER_BASELINE
            ws.dynbal_baseline_liq_col[c] = 0.0
            ws.dynbal_baseline_ice_col[c] = 0.0
            for j in 1:nlevtot
                ws.h2osoi_liq_col[c, j] = 0.0
                ws.h2osoi_ice_col[c, j] = 0.0
                ws.excess_ice_col[c, j] = 0.0
            end
        end
        ws.aquifer_water_baseline = CLM.AQUIFER_WATER_BASELINE

        # Water diagnostic bulk
        wdb = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(wdb, nc, nc, nc, nc)
        for c in 1:nc
            wdb.total_plant_stored_h2o_col[c] = 0.0
        end

        # Water state bulk
        # NB: waterstatebulk_init! calls waterstate_init!(wsb.ws, ...) internally,
        # which would reinitialize ws to NaN if wsb.ws already pointed to it.
        # So we call init first (on the default empty WaterStateData), then replace.
        wsb = CLM.WaterStateBulkData()
        CLM.waterstatebulk_init!(wsb, nc, nc, nc, nc)
        wsb.ws = ws

        # Temperature
        temp = CLM.TemperatureData()
        CLM.temperature_init!(temp, nc, nc, nc, nc)
        for c in 1:nc
            temp.t_h2osfc_col[c] = CLM.TFRZ
            temp.dynbal_baseline_heat_col[c] = 0.0
            for j in 1:nlevtot
                temp.t_soisno_col[c, j] = CLM.TFRZ
            end
            for j in 1:nlevlak
                temp.t_lake_col[c, j] = CLM.TFRZ + 5.0
            end
        end

        # Soil state
        ss = CLM.SoilStateData()
        CLM.soilstate_init!(ss, nc, nc)
        for c in 1:nc
            for j in 1:nlevmaxurbgrnd
                if j <= nlevgrnd
                    ss.watsat_col[c, j] = 0.4
                    ss.csol_col[c, j] = 2.0e6
                else
                    ss.watsat_col[c, j] = 0.0
                end
            end
        end

        # Lake state
        ls = CLM.LakeStateData()
        CLM.lakestate_init!(ls, nc, nc)
        for c in 1:nc
            for j in 1:nlevlak
                ls.lake_icefrac_col[c, j] = 0.0
            end
        end

        # Urban params
        up = CLM.UrbanParamsData()
        CLM.urbanparams_init!(up, nc; nlevurb=nlevurb)
        for l in 1:nc
            up.nlev_improad[l] = 2
            for j in 1:nlevurb
                up.cv_wall[l, j] = 1.5e6
                up.cv_roof[l, j] = 1.2e6
                up.cv_improad[l, j] = 1.8e6
            end
        end

        return col, lun, ws, wdb, wsb, temp, ss, ls, up
    end

    # ===================================================================
    @testset "Constants" begin
        @test CLM.heat_base_temp == CLM.TFRZ
        @test CLM.DeltaLiqMinTemp == CLM.TFRZ
        @test CLM.DeltaLiqMaxTemp == CLM.TFRZ + 35.0
    end

    # ===================================================================
    @testset "_temp_to_heat" begin
        # At freezing point, heat should be zero
        @test CLM._temp_to_heat(CLM.TFRZ, 100.0) == 0.0

        # Above freezing
        temp = CLM.TFRZ + 10.0
        cv = 200.0
        @test CLM._temp_to_heat(temp, cv) ≈ cv * 10.0

        # Below freezing
        temp = CLM.TFRZ - 5.0
        @test CLM._temp_to_heat(temp, cv) ≈ cv * (-5.0)
    end

    # ===================================================================
    @testset "_accumulate_liquid_water_heat (with cv_liquid)" begin
        temp = CLM.TFRZ + 10.0
        h2o = 5.0  # kg/m2
        heat_liq = 0.0
        latent_liq = 0.0
        cv_liq = 0.0

        hl, ll, cvl = CLM._accumulate_liquid_water_heat(temp, h2o, heat_liq, latent_liq, cv_liq)

        cv_expected = h2o * CLM.CPLIQ
        @test cvl ≈ cv_expected
        @test hl ≈ cv_expected * 10.0
        @test ll ≈ h2o * CLM.HFUS
    end

    # ===================================================================
    @testset "_accumulate_liquid_water_heat (without cv_liquid)" begin
        temp = CLM.TFRZ + 10.0
        h2o = 5.0
        heat_liq = 0.0
        latent_liq = 0.0

        hl, ll = CLM._accumulate_liquid_water_heat(temp, h2o, heat_liq, latent_liq)

        cv_expected = h2o * CLM.CPLIQ
        @test hl ≈ cv_expected * 10.0
        @test ll ≈ h2o * CLM.HFUS
    end

    # ===================================================================
    @testset "liquid_water_heat" begin
        # At freezing: heat = 0 + hfus * h2o
        h2o = 3.0
        @test CLM.liquid_water_heat(CLM.TFRZ, h2o) ≈ h2o * CLM.HFUS

        # Above freezing: sensible + latent
        temp = CLM.TFRZ + 20.0
        expected = h2o * CLM.CPLIQ * 20.0 + h2o * CLM.HFUS
        @test CLM.liquid_water_heat(temp, h2o) ≈ expected
    end

    # ===================================================================
    @testset "compute_water_mass_non_lake - zero state" begin
        nc = 3
        col, lun, ws, wdb, wsb, temp, ss, ls, up = setup_basic(nc)

        mask = BitVector([true, true, true])
        water_mass = zeros(nc)

        CLM.compute_water_mass_non_lake!(mask, col, ws, wdb, false, water_mass)

        # Everything is zero, so water mass should be zero
        for c in 1:nc
            @test water_mass[c] ≈ 0.0
        end
    end

    # ===================================================================
    @testset "compute_water_mass_non_lake - with soil water" begin
        nc = 2
        col, lun, ws, wdb, wsb, temp, ss, ls, up = setup_basic(nc)

        # Set soil liquid water in first ground layer (j=1, jj = 1 + nlevsno)
        jj = 1 + nlevsno
        ws.h2osoi_liq_col[1, jj] = 10.0  # kg/m2
        ws.h2osoi_ice_col[1, jj] = 5.0   # kg/m2

        mask = BitVector([true, true])
        water_mass = zeros(nc)

        CLM.compute_water_mass_non_lake!(mask, col, ws, wdb, false, water_mass)

        # Column 1 should have water, column 2 should be zero
        @test water_mass[1] ≈ 15.0
        @test water_mass[2] ≈ 0.0
    end

    # ===================================================================
    @testset "compute_water_mass_non_lake - with surface water" begin
        nc = 2
        col, lun, ws, wdb, wsb, temp, ss, ls, up = setup_basic(nc)

        ws.h2osfc_col[1] = 3.0  # mm H2O

        mask = BitVector([true, true])
        water_mass = zeros(nc)

        CLM.compute_water_mass_non_lake!(mask, col, ws, wdb, false, water_mass)

        @test water_mass[1] ≈ 3.0
        @test water_mass[2] ≈ 0.0
    end

    # ===================================================================
    @testset "compute_water_mass_non_lake - urban impervious road no surface water" begin
        nc = 2
        col, lun, ws, wdb, wsb, temp, ss, ls, up = setup_basic(nc;
            col_types=[CLM.ICOL_ROAD_IMPERV, CLM.ISTSOIL],
            lun_types=[CLM.ISTURB_MD, CLM.ISTSOIL])

        ws.h2osfc_col[1] = 3.0  # should NOT be counted
        ws.h2osfc_col[2] = 3.0  # should be counted

        mask = BitVector([true, true])
        water_mass = zeros(nc)

        CLM.compute_water_mass_non_lake!(mask, col, ws, wdb, false, water_mass)

        # Impervious road should not add surface water
        @test water_mass[1] ≈ 0.0
        # Normal soil should add it
        @test water_mass[2] ≈ 3.0
    end

    # ===================================================================
    @testset "compute_water_mass_non_lake - dynbal baselines" begin
        nc = 2
        col, lun, ws, wdb, wsb, temp, ss, ls, up = setup_basic(nc)

        ws.h2osfc_col[1] = 10.0
        ws.dynbal_baseline_liq_col[1] = 3.0
        ws.dynbal_baseline_ice_col[1] = 2.0

        mask = BitVector([true, true])
        water_mass = zeros(nc)

        # Without subtracting baselines
        CLM.compute_water_mass_non_lake!(mask, col, ws, wdb, false, water_mass)
        @test water_mass[1] ≈ 10.0

        # With subtracting baselines
        CLM.compute_water_mass_non_lake!(mask, col, ws, wdb, true, water_mass)
        @test water_mass[1] ≈ 10.0 - 3.0 - 2.0
    end

    # ===================================================================
    @testset "compute_water_mass_non_lake - snow layers" begin
        nc = 1
        col, lun, ws, wdb, wsb, temp, ss, ls, up = setup_basic(nc; snl_vals=[-2])

        # Snow layers: snl=-2, so active layers are j=-1 and j=0
        # j=-1 -> jj = -1 + nlevsno = nlevsno - 1
        # j=0  -> jj = 0 + nlevsno  = nlevsno
        ws.h2osoi_liq_col[1, nlevsno - 1] = 1.0
        ws.h2osoi_ice_col[1, nlevsno - 1] = 2.0
        ws.h2osoi_liq_col[1, nlevsno] = 0.5
        ws.h2osoi_ice_col[1, nlevsno] = 3.0

        mask = BitVector([true])
        water_mass = zeros(1)

        CLM.compute_water_mass_non_lake!(mask, col, ws, wdb, false, water_mass)

        @test water_mass[1] ≈ (1.0 + 0.5) + (2.0 + 3.0)
    end

    # ===================================================================
    @testset "compute_water_mass_non_lake - aquifer water" begin
        nc = 2
        col, lun, ws, wdb, wsb, temp, ss, ls, up = setup_basic(nc;
            hydro_active=[true, false])

        ws.wa_col[1] = CLM.AQUIFER_WATER_BASELINE + 100.0
        ws.wa_col[2] = CLM.AQUIFER_WATER_BASELINE + 100.0

        mask = BitVector([true, true])
        water_mass = zeros(nc)

        CLM.compute_water_mass_non_lake!(mask, col, ws, wdb, false, water_mass)

        # Column 1 is hydrologically active, should have +100 from aquifer
        @test water_mass[1] ≈ 100.0
        # Column 2 is not active, aquifer water excluded
        @test water_mass[2] ≈ 0.0
    end

    # ===================================================================
    @testset "accumulate_soil_liq_ice_mass_non_lake - sunwall/shadewall excluded" begin
        nc = 3
        col, lun, ws, wdb, wsb, temp, ss, ls, up = setup_basic(nc;
            col_types=[CLM.ICOL_SUNWALL, CLM.ICOL_SHADEWALL, CLM.ISTSOIL],
            lun_types=[CLM.ISTURB_MD, CLM.ISTURB_MD, CLM.ISTSOIL])

        # Set soil water in first ground layer
        jj = 1 + nlevsno
        for c in 1:nc
            ws.h2osoi_liq_col[c, jj] = 10.0
            ws.h2osoi_ice_col[c, jj] = 5.0
        end

        mask = BitVector([true, true, true])
        liquid_mass = zeros(nc)
        ice_mass = zeros(nc)

        CLM.accumulate_soil_liq_ice_mass_non_lake!(mask, col, ws, liquid_mass, ice_mass)

        # Sunwall and shadewall should have zero
        @test liquid_mass[1] ≈ 0.0
        @test ice_mass[1] ≈ 0.0
        @test liquid_mass[2] ≈ 0.0
        @test ice_mass[2] ≈ 0.0

        # Normal soil column should have water
        @test liquid_mass[3] ≈ 10.0
        @test ice_mass[3] ≈ 5.0
    end

    # ===================================================================
    @testset "accumulate_soil_liq_ice_mass_non_lake - roof limited to nlevurb" begin
        nc = 1
        col, lun, ws, wdb, wsb, temp, ss, ls, up = setup_basic(nc;
            col_types=[CLM.ICOL_ROOF],
            lun_types=[CLM.ISTURB_MD])

        # Set water in layer nlevurb (should be counted) and nlevurb+1 (should not)
        jj_urb = nlevurb + nlevsno
        jj_below = (nlevurb + 1) + nlevsno
        ws.h2osoi_liq_col[1, jj_urb] = 7.0
        if jj_below <= nlevtot
            ws.h2osoi_liq_col[1, jj_below] = 99.0
        end

        mask = BitVector([true])
        liquid_mass = zeros(1)
        ice_mass = zeros(1)

        CLM.accumulate_soil_liq_ice_mass_non_lake!(mask, col, ws, liquid_mass, ice_mass)

        @test liquid_mass[1] ≈ 7.0
    end

    # ===================================================================
    @testset "accumulate_soil_liq_ice_mass_non_lake - excess ice" begin
        nc = 1
        col, lun, ws, wdb, wsb, temp, ss, ls, up = setup_basic(nc)

        jj = 1 + nlevsno
        ws.h2osoi_ice_col[1, jj] = 5.0
        ws.excess_ice_col[1, jj] = 2.0

        mask = BitVector([true])
        liquid_mass = zeros(1)
        ice_mass = zeros(1)

        CLM.accumulate_soil_liq_ice_mass_non_lake!(mask, col, ws, liquid_mass, ice_mass)

        @test ice_mass[1] ≈ 7.0  # 5.0 + 2.0
    end

    # ===================================================================
    @testset "compute_liq_ice_mass_lake - basic" begin
        nc = 1
        col, lun, ws, wdb, wsb, temp, ss, ls, up = setup_basic(nc;
            lun_types=[CLM.ISTDLAK])

        # Set lake depth and ice fraction
        for j in 1:nlevlak
            col.dz_lake[1, j] = 2.0
            ls.lake_icefrac_col[1, j] = 0.3
        end

        # Set some soil water under the lake
        jj = 1 + nlevsno
        ws.h2osoi_liq_col[1, jj] = 4.0
        ws.h2osoi_ice_col[1, jj] = 1.0

        mask = BitVector([true])
        liquid_mass = zeros(1)
        ice_mass = zeros(1)

        CLM.compute_liq_ice_mass_lake!(mask, col, ws, ls, true, liquid_mass, ice_mass)

        # Lake contribution: for each layer, liq = dz * denh2o * (1 - icefrac) * ratio
        expected_lake_liq = nlevlak * 2.0 * CLM.DENH2O * 0.7
        expected_lake_ice = nlevlak * 2.0 * CLM.DENH2O * 0.3

        # Snow (none), soil water
        expected_soil_liq = 4.0
        expected_soil_ice = 1.0

        # With dynbal baselines at zero
        @test liquid_mass[1] ≈ expected_lake_liq + expected_soil_liq
        @test ice_mass[1] ≈ expected_lake_ice + expected_soil_ice
    end

    # ===================================================================
    @testset "compute_liq_ice_mass_lake - without lake water" begin
        nc = 1
        col, lun, ws, wdb, wsb, temp, ss, ls, up = setup_basic(nc;
            lun_types=[CLM.ISTDLAK])

        for j in 1:nlevlak
            col.dz_lake[1, j] = 2.0
            ls.lake_icefrac_col[1, j] = 0.5
        end

        mask = BitVector([true])
        liquid_mass = zeros(1)
        ice_mass = zeros(1)

        # Don't add lake water
        CLM.compute_liq_ice_mass_lake!(mask, col, ws, ls, false, liquid_mass, ice_mass)

        # Should have zero from lake, only snow and soil (both zero here)
        @test liquid_mass[1] ≈ 0.0
        @test ice_mass[1] ≈ 0.0
    end

    # ===================================================================
    @testset "accumulate_liq_ice_mass_lake" begin
        nc = 1
        col, lun, ws, wdb, wsb, temp, ss, ls, up = setup_basic(nc;
            lun_types=[CLM.ISTDLAK])

        for j in 1:nlevlak
            col.dz_lake[1, j] = 3.0
            ls.lake_icefrac_col[1, j] = 0.0  # all liquid
        end

        mask = BitVector([true])
        liquid_mass = zeros(1)
        ice_mass = zeros(1)

        CLM.accumulate_liq_ice_mass_lake!(mask, col, ls, 1.0, liquid_mass, ice_mass)

        expected_liq = nlevlak * 3.0 * CLM.DENH2O
        @test liquid_mass[1] ≈ expected_liq
        @test ice_mass[1] ≈ 0.0
    end

    # ===================================================================
    @testset "accumulate_liq_ice_mass_lake - all frozen" begin
        nc = 1
        col, lun, ws, wdb, wsb, temp, ss, ls, up = setup_basic(nc;
            lun_types=[CLM.ISTDLAK])

        for j in 1:nlevlak
            col.dz_lake[1, j] = 1.0
            ls.lake_icefrac_col[1, j] = 1.0
        end

        mask = BitVector([true])
        liquid_mass = zeros(1)
        ice_mass = zeros(1)

        CLM.accumulate_liq_ice_mass_lake!(mask, col, ls, 1.0, liquid_mass, ice_mass)

        expected_ice = nlevlak * 1.0 * CLM.DENH2O
        @test liquid_mass[1] ≈ 0.0
        @test ice_mass[1] ≈ expected_ice
    end

    # ===================================================================
    @testset "compute_heat_non_lake - zero state at freezing" begin
        nc = 1
        col, lun, ws, wdb, wsb, temp, ss, ls, up = setup_basic(nc)

        mask = BitVector([true])
        heat = zeros(1)
        heat_liquid = zeros(1)
        cv_liquid = zeros(1)

        # All temperatures at freezing, all water zero -> heat should be zero
        CLM.compute_heat_non_lake!(mask, col, lun, up, ss, temp, wsb, wdb,
                                    heat, heat_liquid, cv_liquid)

        @test heat[1] ≈ 0.0 atol=1e-10
    end

    # ===================================================================
    @testset "compute_heat_non_lake - with soil water above freezing" begin
        nc = 1
        col, lun, ws, wdb, wsb, temp, ss, ls, up = setup_basic(nc)

        jj = 1 + nlevsno
        temp.t_soisno_col[1, jj] = CLM.TFRZ + 10.0
        wsb.ws.h2osoi_liq_col[1, jj] = 5.0  # kg/m2

        mask = BitVector([true])
        heat = zeros(1)
        heat_liquid = zeros(1)
        cv_liquid = zeros(1)

        CLM.compute_heat_non_lake!(mask, col, lun, up, ss, temp, wsb, wdb,
                                    heat, heat_liquid, cv_liquid)

        # Expected: soil dry mass heat + liquid water heat + latent heat
        # Dry mass: csol * (1-watsat) * dz * (T - Tfrz)
        dry_heat = ss.csol_col[1, 1] * (1.0 - ss.watsat_col[1, 1]) * col.dz[1, jj] * 10.0
        # Liquid water: cpliq * h2o * (T - Tfrz) + hfus * h2o
        liq_heat = CLM.CPLIQ * 5.0 * 10.0 + CLM.HFUS * 5.0

        @test heat[1] ≈ dry_heat + liq_heat atol=1e-6
        @test cv_liquid[1] ≈ 5.0 * CLM.CPLIQ atol=1e-6
    end

    # ===================================================================
    @testset "compute_heat_non_lake - with frozen soil water" begin
        nc = 1
        col, lun, ws, wdb, wsb, temp, ss, ls, up = setup_basic(nc)

        jj = 1 + nlevsno
        temp.t_soisno_col[1, jj] = CLM.TFRZ - 5.0
        wsb.ws.h2osoi_ice_col[1, jj] = 8.0  # kg/m2

        mask = BitVector([true])
        heat = zeros(1)
        heat_liquid = zeros(1)
        cv_liquid = zeros(1)

        CLM.compute_heat_non_lake!(mask, col, lun, up, ss, temp, wsb, wdb,
                                    heat, heat_liquid, cv_liquid)

        # Expected: dry mass heat + ice heat
        dry_heat = ss.csol_col[1, 1] * (1.0 - ss.watsat_col[1, 1]) * col.dz[1, jj] * (-5.0)
        ice_heat = CLM.CPICE * 8.0 * (-5.0)

        @test heat[1] ≈ dry_heat + ice_heat atol=1e-6
        @test cv_liquid[1] ≈ 0.0 atol=1e-10
    end

    # ===================================================================
    @testset "compute_heat_non_lake - dynbal baseline subtraction" begin
        nc = 1
        col, lun, ws, wdb, wsb, temp, ss, ls, up = setup_basic(nc)

        temp.dynbal_baseline_heat_col[1] = 100.0

        mask = BitVector([true])
        heat = zeros(1)
        heat_liquid = zeros(1)
        cv_liquid = zeros(1)

        CLM.compute_heat_non_lake!(mask, col, lun, up, ss, temp, wsb, wdb,
                                    heat, heat_liquid, cv_liquid)

        # With zero water and freezing temps, heat = 0 - baseline = -100
        @test heat[1] ≈ -100.0 atol=1e-10
    end

    # ===================================================================
    @testset "compute_heat_lake - basic" begin
        nc = 1
        col, lun, ws, wdb, wsb, temp, ss, ls, up = setup_basic(nc;
            lun_types=[CLM.ISTDLAK])

        # Lake all liquid at TFRZ + 5
        for j in 1:nlevlak
            col.dz_lake[1, j] = 1.0
            ls.lake_icefrac_col[1, j] = 0.0
            temp.t_lake_col[1, j] = CLM.TFRZ + 5.0
        end

        mask = BitVector([true])
        heat = zeros(1)
        heat_liquid = zeros(1)
        cv_liquid = zeros(1)

        CLM.compute_heat_lake!(mask, col, ss, temp, wsb, ls,
                                heat, heat_liquid, cv_liquid)

        # Lake heat should be positive (above freezing)
        @test heat[1] > 0.0
    end

    # ===================================================================
    @testset "compute_heat_lake - frozen lake" begin
        nc = 1
        col, lun, ws, wdb, wsb, temp, ss, ls, up = setup_basic(nc;
            lun_types=[CLM.ISTDLAK])

        # Lake all frozen at TFRZ - 10
        for j in 1:nlevlak
            col.dz_lake[1, j] = 1.0
            ls.lake_icefrac_col[1, j] = 1.0
            temp.t_lake_col[1, j] = CLM.TFRZ - 10.0
        end

        mask = BitVector([true])
        heat = zeros(1)
        heat_liquid = zeros(1)
        cv_liquid = zeros(1)

        CLM.compute_heat_lake!(mask, col, ss, temp, wsb, ls,
                                heat, heat_liquid, cv_liquid)

        # All ice below freezing -> heat should be negative
        @test heat[1] < 0.0
        # No liquid water contribution to heat_liquid from lake
        # (heat_liquid should only come from soil if any)
        @test heat_liquid[1] ≈ 0.0 atol=1e-10
    end

    # ===================================================================
    @testset "accumulate_heat_lake" begin
        nc = 1
        col, lun, ws, wdb, wsb, temp, ss, ls, up = setup_basic(nc;
            lun_types=[CLM.ISTDLAK])

        # Single lake layer, all liquid at TFRZ + 10
        for j in 1:nlevlak
            col.dz_lake[1, j] = 0.0
            ls.lake_icefrac_col[1, j] = 0.0
            temp.t_lake_col[1, j] = CLM.TFRZ
        end
        col.dz_lake[1, 1] = 2.0
        temp.t_lake_col[1, 1] = CLM.TFRZ + 10.0

        mask = BitVector([true])
        heat = zeros(1)

        CLM.accumulate_heat_lake!(mask, col, temp, ls, heat)

        # Expected: cpliq * mass * 10 + hfus * mass
        mass_liq = 2.0 * CLM.DENH2O
        expected_heat = CLM.CPLIQ * mass_liq * 10.0 + CLM.HFUS * mass_liq

        @test heat[1] ≈ expected_heat atol=1e-6
    end

    # ===================================================================
    @testset "adjust_delta_heat_for_delta_liq - positive delta_liq" begin
        delta_liq = [10.0]
        temp1 = [CLM.TFRZ + 5.0]
        temp2 = [CLM.TFRZ + 15.0]
        delta_heat = [0.0]

        CLM.adjust_delta_heat_for_delta_liq!(1:1, delta_liq, temp1, temp2, delta_heat)

        # delta_liq > 0: use temp2
        expected = CLM.liquid_water_heat(CLM.TFRZ + 15.0, 10.0)
        @test delta_heat[1] ≈ -expected
    end

    # ===================================================================
    @testset "adjust_delta_heat_for_delta_liq - negative delta_liq" begin
        delta_liq = [-5.0]
        temp1 = [CLM.TFRZ + 10.0]
        temp2 = [CLM.TFRZ + 20.0]
        delta_heat = [0.0]

        CLM.adjust_delta_heat_for_delta_liq!(1:1, delta_liq, temp1, temp2, delta_heat)

        # delta_liq < 0: use temp1
        expected = CLM.liquid_water_heat(CLM.TFRZ + 10.0, -5.0)
        @test delta_heat[1] ≈ -expected
    end

    # ===================================================================
    @testset "adjust_delta_heat_for_delta_liq - zero delta_liq" begin
        delta_liq = [0.0]
        temp1 = [CLM.TFRZ + 10.0]
        temp2 = [CLM.TFRZ + 20.0]
        delta_heat = [42.0]

        CLM.adjust_delta_heat_for_delta_liq!(1:1, delta_liq, temp1, temp2, delta_heat)

        # No change
        @test delta_heat[1] ≈ 42.0
    end

    # ===================================================================
    @testset "adjust_delta_heat_for_delta_liq - temperature clamping" begin
        delta_liq = [1.0]
        # Extreme temperatures that should be clamped
        temp1 = [100.0]       # way below DeltaLiqMinTemp
        temp2 = [500.0]       # way above DeltaLiqMaxTemp
        delta_heat = [0.0]

        CLM.adjust_delta_heat_for_delta_liq!(1:1, delta_liq, temp1, temp2, delta_heat)

        # delta_liq > 0: use temp2 clamped to DeltaLiqMaxTemp
        expected = CLM.liquid_water_heat(CLM.DeltaLiqMaxTemp, 1.0)
        @test delta_heat[1] ≈ -expected
    end

    # ===================================================================
    @testset "mask filtering" begin
        nc = 3
        col, lun, ws, wdb, wsb, temp, ss, ls, up = setup_basic(nc)

        # Set water only on column 2
        ws.h2osfc_col[2] = 10.0

        # Mask only includes columns 1 and 3 (not 2)
        mask = BitVector([true, false, true])
        water_mass = fill(NaN, nc)

        CLM.compute_water_mass_non_lake!(mask, col, ws, wdb, false, water_mass)

        @test water_mass[1] ≈ 0.0
        @test isnan(water_mass[2])  # should be untouched
        @test water_mass[3] ≈ 0.0
    end

    # ===================================================================
    @testset "compute_water_mass_lake" begin
        nc = 1
        col, lun, ws, wdb, wsb, temp, ss, ls, up = setup_basic(nc;
            lun_types=[CLM.ISTDLAK])

        for j in 1:nlevlak
            col.dz_lake[1, j] = 1.0
            ls.lake_icefrac_col[1, j] = 0.5
        end

        mask = BitVector([true])
        water_mass = zeros(1)

        CLM.compute_water_mass_lake!(mask, col, ws, ls, true, water_mass)

        expected = nlevlak * 1.0 * CLM.DENH2O  # all water (liq + ice)
        @test water_mass[1] ≈ expected
    end

    # ===================================================================
    @testset "consistency: water mass = liq + ice" begin
        nc = 2
        col, lun, ws, wdb, wsb, temp, ss, ls, up = setup_basic(nc)

        # Set various water components
        jj = 1 + nlevsno
        ws.h2osoi_liq_col[1, jj] = 10.0
        ws.h2osoi_ice_col[1, jj] = 5.0
        ws.h2osfc_col[1] = 2.0
        ws.h2osno_no_layers_col[2] = 3.0

        mask = BitVector([true, true])
        water_mass = zeros(nc)
        liquid_mass = zeros(nc)
        ice_mass = zeros(nc)

        CLM.compute_water_mass_non_lake!(mask, col, ws, wdb, false, water_mass)
        CLM.compute_liq_ice_mass_non_lake!(mask, col, ws, wdb, false,
                                            liquid_mass, ice_mass)

        for c in 1:nc
            @test water_mass[c] ≈ liquid_mass[c] + ice_mass[c]
        end
    end

    # ===================================================================
    @testset "accumulate_soil_heat_non_lake - urban wall heat" begin
        nc = 1
        col, lun, ws, wdb, wsb, temp, ss, ls, up = setup_basic(nc;
            col_types=[CLM.ICOL_SUNWALL],
            lun_types=[CLM.ISTURB_MD])

        jj = 1 + nlevsno
        temp.t_soisno_col[1, jj] = CLM.TFRZ + 20.0

        mask = BitVector([true])
        heat = zeros(1)
        heat_liquid = zeros(1)
        cv_liquid = zeros(1)

        CLM.accumulate_soil_heat_non_lake!(mask, col, lun, up, ss, temp, wsb,
                                            heat, heat_liquid, cv_liquid)

        # Wall has dry mass heat from cv_wall but no water
        expected_dry = up.cv_wall[1, 1] * col.dz[1, jj] * 20.0
        @test heat[1] ≈ expected_dry atol=1e-6
        @test cv_liquid[1] ≈ 0.0 atol=1e-10
    end

end
