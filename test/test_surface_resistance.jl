@testset "SurfaceResistance" begin
    # Setup: initialize varpar
    vp = CLM.varpar
    saved_nlevsno = vp.nlevsno
    saved_nlevgrnd = vp.nlevgrnd
    saved_nlevurb = vp.nlevurb
    saved_nlevmaxurbgrnd = vp.nlevmaxurbgrnd
    saved_nlevsoi = vp.nlevsoi

    vp.nlevsno = 5
    vp.nlevgrnd = 10
    vp.nlevurb = 5
    vp.nlevmaxurbgrnd = 10
    vp.nlevsoi = 8

    nlevsno = vp.nlevsno
    nlevgrnd = vp.nlevgrnd
    nlevmaxurbgrnd = vp.nlevmaxurbgrnd
    joff = nlevsno  # offset for combined snow+soil indexing

    nc = 6  # number of columns
    nl = 6  # number of landunits

    # Save and restore surface resistance control/params state
    saved_method = CLM.surface_resistance_ctrl.soil_resis_method
    saved_d_max = CLM.surface_resistance_params.d_max
    saved_frac = CLM.surface_resistance_params.frac_sat_soil_dsl_init

    function setup_test_data()
        col = CLM.ColumnData()
        CLM.column_init!(col, nc)

        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, nl)

        temp = CLM.TemperatureData()
        CLM.temperature_init!(temp, 0, nc, nl, 1)

        soilstate = CLM.SoilStateData()
        CLM.soilstate_init!(soilstate, 0, nc)

        waterstatebulk = CLM.WaterStateBulkData()
        CLM.waterstatebulk_init!(waterstatebulk, nc, 0, nl, 1)

        waterdiagbulk = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(waterdiagbulk, nc, 0, nl, 1)

        mask_nolakec = falses(nc)
        bounds = 1:nc

        return col, lun, temp, soilstate, waterstatebulk, waterdiagbulk,
               mask_nolakec, bounds
    end

    @testset "Control functions" begin
        # Test do_soilevap_beta and do_soil_resistance_sl14
        CLM.soil_resistance_read_nl!(soil_resis_method=CLM.SOIL_RESIS_LEEPIELKE_1992)
        @test CLM.do_soilevap_beta() == true
        @test CLM.do_soil_resistance_sl14() == false

        CLM.soil_resistance_read_nl!(soil_resis_method=CLM.SOIL_RESIS_SL_14)
        @test CLM.do_soilevap_beta() == false
        @test CLM.do_soil_resistance_sl14() == true
    end

    @testset "Parameter setting" begin
        CLM.surface_resistance_read_params!(d_max=20.0, frac_sat_soil_dsl_init=0.7)
        @test CLM.surface_resistance_params.d_max == 20.0
        @test CLM.surface_resistance_params.frac_sat_soil_dsl_init == 0.7

        # Restore defaults for subsequent tests
        CLM.surface_resistance_read_params!(d_max=15.0, frac_sat_soil_dsl_init=0.8)
    end

    @testset "Lee-Pielke 1992 beta — soil column below field capacity" begin
        CLM.soil_resistance_read_nl!(soil_resis_method=CLM.SOIL_RESIS_LEEPIELKE_1992)
        col, lun, temp, ss, wsb, wdb, mask, bounds = setup_test_data()

        c = 1
        mask[c] = true
        col.landunit[c] = 1
        col.itype[c] = CLM.ISTSOIL
        lun.itype[1] = CLM.ISTSOIL

        # Set layer thickness
        col.dz[c, 1 + joff] = 0.02  # 2cm top layer

        # Set soil water — below field capacity
        wsb.ws.h2osoi_liq_col[c, 1 + joff] = 2.0   # low liquid
        wsb.ws.h2osoi_ice_col[c, 1 + joff] = 0.0

        # Soil properties
        ss.watsat_col[c, 1] = 0.45
        ss.watfc_col[c, 1] = 0.30

        # Fractions
        wdb.frac_sno_col[c] = 0.0
        wdb.frac_h2osfc_col[c] = 0.0

        CLM.calc_soilevap_resis!(col, lun, ss, wsb, wdb, temp, mask, bounds)

        # soilbeta should be between 0 and 1
        @test ss.soilbeta_col[c] > 0.0
        @test ss.soilbeta_col[c] < 1.0

        # Verify formula: 0.25*(1-cos(π*fac_fc))^2
        wx = (2.0 / CLM.DENH2O + 0.0 / CLM.DENICE) / 0.02
        fac_fc = min(1.0, wx / 0.30)
        fac_fc = max(fac_fc, 0.01)
        expected_beta = 0.25 * (1.0 - cos(π * fac_fc))^2.0
        @test ss.soilbeta_col[c] ≈ expected_beta atol=1e-10
    end

    @testset "Lee-Pielke 1992 beta — soil column above field capacity" begin
        CLM.soil_resistance_read_nl!(soil_resis_method=CLM.SOIL_RESIS_LEEPIELKE_1992)
        col, lun, temp, ss, wsb, wdb, mask, bounds = setup_test_data()

        c = 1
        mask[c] = true
        col.landunit[c] = 1
        col.itype[c] = CLM.ISTSOIL
        lun.itype[1] = CLM.ISTSOIL

        col.dz[c, 1 + joff] = 0.02

        # Water content above field capacity
        wsb.ws.h2osoi_liq_col[c, 1 + joff] = 8.0  # high liquid
        wsb.ws.h2osoi_ice_col[c, 1 + joff] = 0.0

        ss.watsat_col[c, 1] = 0.45
        ss.watfc_col[c, 1] = 0.30

        wdb.frac_sno_col[c] = 0.0
        wdb.frac_h2osfc_col[c] = 0.0

        CLM.calc_soilevap_resis!(col, lun, ss, wsb, wdb, temp, mask, bounds)

        # Above field capacity → soilbeta = 1.0
        @test ss.soilbeta_col[c] ≈ 1.0
    end

    @testset "Lee-Pielke 1992 beta — with snow and h2osfc" begin
        CLM.soil_resistance_read_nl!(soil_resis_method=CLM.SOIL_RESIS_LEEPIELKE_1992)
        col, lun, temp, ss, wsb, wdb, mask, bounds = setup_test_data()

        c = 1
        mask[c] = true
        col.landunit[c] = 1
        col.itype[c] = CLM.ISTSOIL
        lun.itype[1] = CLM.ISTSOIL

        col.dz[c, 1 + joff] = 0.02
        wsb.ws.h2osoi_liq_col[c, 1 + joff] = 2.0
        wsb.ws.h2osoi_ice_col[c, 1 + joff] = 0.0
        ss.watsat_col[c, 1] = 0.45
        ss.watfc_col[c, 1] = 0.30

        # 30% snow, 20% surface water
        wdb.frac_sno_col[c] = 0.3
        wdb.frac_h2osfc_col[c] = 0.2

        CLM.calc_soilevap_resis!(col, lun, ss, wsb, wdb, temp, mask, bounds)

        # Snow/h2osfc fractions should increase soilbeta towards 1
        wx = (2.0 / CLM.DENH2O + 0.0 / CLM.DENICE) / 0.02
        fac_fc = min(1.0, wx / 0.30)
        fac_fc = max(fac_fc, 0.01)
        bare_part = (1.0 - 0.3 - 0.2) * 0.25 * (1.0 - cos(π * fac_fc))^2.0
        expected = bare_part + 0.3 + 0.2
        @test ss.soilbeta_col[c] ≈ expected atol=1e-10
    end

    @testset "Lee-Pielke 1992 beta — wet/ice landunits" begin
        CLM.soil_resistance_read_nl!(soil_resis_method=CLM.SOIL_RESIS_LEEPIELKE_1992)
        col, lun, temp, ss, wsb, wdb, mask, bounds = setup_test_data()

        # Wet landunit
        c = 1
        mask[c] = true
        col.landunit[c] = 1
        col.itype[c] = CLM.ISTWET
        lun.itype[1] = CLM.ISTWET

        # Ice landunit
        c2 = 2
        mask[c2] = true
        col.landunit[c2] = 2
        col.itype[c2] = CLM.ISTICE
        lun.itype[2] = CLM.ISTICE

        CLM.calc_soilevap_resis!(col, lun, ss, wsb, wdb, temp, mask, bounds)

        @test ss.soilbeta_col[c] ≈ 1.0
        @test ss.soilbeta_col[c2] ≈ 1.0
    end

    @testset "Lee-Pielke 1992 beta — urban columns" begin
        CLM.soil_resistance_read_nl!(soil_resis_method=CLM.SOIL_RESIS_LEEPIELKE_1992)
        col, lun, temp, ss, wsb, wdb, mask, bounds = setup_test_data()

        # Set up 4 urban column types
        for (idx, ctype) in enumerate([CLM.ICOL_ROAD_PERV, CLM.ICOL_SUNWALL,
                                        CLM.ICOL_ROOF, CLM.ICOL_ROAD_IMPERV])
            mask[idx] = true
            col.landunit[idx] = idx
            col.itype[idx] = ctype
            lun.itype[idx] = CLM.ISTURB_TBD  # urban but not wet/ice
        end

        CLM.calc_soilevap_resis!(col, lun, ss, wsb, wdb, temp, mask, bounds)

        for idx in 1:4
            @test ss.soilbeta_col[idx] ≈ 0.0
        end
    end

    @testset "SL14 resistance — soil column" begin
        CLM.soil_resistance_read_nl!(soil_resis_method=CLM.SOIL_RESIS_SL_14)
        CLM.surface_resistance_read_params!(d_max=15.0, frac_sat_soil_dsl_init=0.8)
        col, lun, temp, ss, wsb, wdb, mask, bounds = setup_test_data()

        c = 1
        mask[c] = true
        col.landunit[c] = 1
        col.itype[c] = CLM.ISTSOIL
        lun.itype[1] = CLM.ISTSOIL

        col.dz[c, 1 + joff] = 0.02
        wsb.ws.h2osoi_liq_col[c, 1 + joff] = 5.0
        wsb.ws.h2osoi_ice_col[c, 1 + joff] = 0.0

        ss.watsat_col[c, 1] = 0.45
        ss.bsw_col[c, 1] = 5.0
        ss.sucsat_col[c, 1] = 100.0

        temp.t_soisno_col[c, 1 + joff] = 290.0

        CLM.calc_soilevap_resis!(col, lun, ss, wsb, wdb, temp, mask, bounds)

        # soilresis should be finite positive
        @test ss.soilresis_col[c] > 0.0
        @test ss.soilresis_col[c] < 1.0e6
        @test isfinite(ss.soilresis_col[c])

        # dsl should be non-negative and bounded
        @test ss.dsl_col[c] >= 0.0
        @test ss.dsl_col[c] <= 200.0
    end

    @testset "SL14 resistance — very wet soil gives small resistance" begin
        CLM.soil_resistance_read_nl!(soil_resis_method=CLM.SOIL_RESIS_SL_14)
        CLM.surface_resistance_read_params!(d_max=15.0, frac_sat_soil_dsl_init=0.8)
        col, lun, temp, ss, wsb, wdb, mask, bounds = setup_test_data()

        c = 1
        mask[c] = true
        col.landunit[c] = 1
        col.itype[c] = CLM.ISTSOIL
        lun.itype[1] = CLM.ISTSOIL

        col.dz[c, 1 + joff] = 0.02
        # Nearly saturated soil
        wsb.ws.h2osoi_liq_col[c, 1 + joff] = 0.02 * 1000.0 * 0.44  # close to saturation
        wsb.ws.h2osoi_ice_col[c, 1 + joff] = 0.0

        ss.watsat_col[c, 1] = 0.45
        ss.bsw_col[c, 1] = 5.0
        ss.sucsat_col[c, 1] = 100.0

        temp.t_soisno_col[c, 1 + joff] = 290.0

        CLM.calc_soilevap_resis!(col, lun, ss, wsb, wdb, temp, mask, bounds)

        # When wet: dsl should be small (clamped by max(0.001,...) in numerator),
        # resistance should be relatively low
        @test ss.dsl_col[c] < 1.0
        @test ss.soilresis_col[c] < 100.0
    end

    @testset "SL14 resistance — dry soil gives large resistance" begin
        CLM.soil_resistance_read_nl!(soil_resis_method=CLM.SOIL_RESIS_SL_14)
        CLM.surface_resistance_read_params!(d_max=15.0, frac_sat_soil_dsl_init=0.8)
        col, lun, temp, ss, wsb, wdb, mask, bounds = setup_test_data()

        c = 1
        mask[c] = true
        col.landunit[c] = 1
        col.itype[c] = CLM.ISTSOIL
        lun.itype[1] = CLM.ISTSOIL

        col.dz[c, 1 + joff] = 0.02
        # Very dry soil
        wsb.ws.h2osoi_liq_col[c, 1 + joff] = 0.001
        wsb.ws.h2osoi_ice_col[c, 1 + joff] = 0.0

        ss.watsat_col[c, 1] = 0.45
        ss.bsw_col[c, 1] = 5.0
        ss.sucsat_col[c, 1] = 100.0

        temp.t_soisno_col[c, 1 + joff] = 290.0

        CLM.calc_soilevap_resis!(col, lun, ss, wsb, wdb, temp, mask, bounds)

        # Dry soil should give high resistance
        @test ss.soilresis_col[c] > 100.0
        @test ss.dsl_col[c] > 0.0
    end

    @testset "SL14 resistance — wet/ice landunits" begin
        CLM.soil_resistance_read_nl!(soil_resis_method=CLM.SOIL_RESIS_SL_14)
        col, lun, temp, ss, wsb, wdb, mask, bounds = setup_test_data()

        c = 1
        mask[c] = true
        col.landunit[c] = 1
        col.itype[c] = CLM.ISTWET
        lun.itype[1] = CLM.ISTWET

        c2 = 2
        mask[c2] = true
        col.landunit[c2] = 2
        col.itype[c2] = CLM.ISTICE
        lun.itype[2] = CLM.ISTICE

        CLM.calc_soilevap_resis!(col, lun, ss, wsb, wdb, temp, mask, bounds)

        @test ss.soilresis_col[c] ≈ 0.0
        @test ss.soilresis_col[c2] ≈ 0.0
    end

    @testset "SL14 resistance — urban columns" begin
        CLM.soil_resistance_read_nl!(soil_resis_method=CLM.SOIL_RESIS_SL_14)
        col, lun, temp, ss, wsb, wdb, mask, bounds = setup_test_data()

        for (idx, ctype) in enumerate([CLM.ICOL_ROAD_PERV, CLM.ICOL_SUNWALL,
                                        CLM.ICOL_ROOF, CLM.ICOL_ROAD_IMPERV])
            mask[idx] = true
            col.landunit[idx] = idx
            col.itype[idx] = ctype
            lun.itype[idx] = CLM.ISTURB_TBD
        end

        CLM.calc_soilevap_resis!(col, lun, ss, wsb, wdb, temp, mask, bounds)

        for idx in 1:4
            @test ss.soilresis_col[idx] ≈ 1.0e6
        end
    end

    @testset "Dispatch — method selection" begin
        col, lun, temp, ss, wsb, wdb, mask, bounds = setup_test_data()

        c = 1
        mask[c] = true
        col.landunit[c] = 1
        col.itype[c] = CLM.ISTSOIL
        lun.itype[1] = CLM.ISTSOIL
        col.dz[c, 1 + joff] = 0.02
        wsb.ws.h2osoi_liq_col[c, 1 + joff] = 5.0
        wsb.ws.h2osoi_ice_col[c, 1 + joff] = 0.0
        ss.watsat_col[c, 1] = 0.45
        ss.watfc_col[c, 1] = 0.30
        ss.bsw_col[c, 1] = 5.0
        ss.sucsat_col[c, 1] = 100.0
        temp.t_soisno_col[c, 1 + joff] = 290.0
        wdb.frac_sno_col[c] = 0.0
        wdb.frac_h2osfc_col[c] = 0.0

        # Test with Lee-Pielke method
        CLM.soil_resistance_read_nl!(soil_resis_method=CLM.SOIL_RESIS_LEEPIELKE_1992)
        ss.soilbeta_col[c] = NaN
        CLM.calc_soilevap_resis!(col, lun, ss, wsb, wdb, temp, mask, bounds)
        @test !isnan(ss.soilbeta_col[c])

        # Test with SL14 method
        CLM.soil_resistance_read_nl!(soil_resis_method=CLM.SOIL_RESIS_SL_14)
        ss.soilresis_col[c] = NaN
        ss.dsl_col[c] = NaN
        CLM.calc_soilevap_resis!(col, lun, ss, wsb, wdb, temp, mask, bounds)
        @test !isnan(ss.soilresis_col[c])
        @test !isnan(ss.dsl_col[c])
    end

    @testset "Masked columns are skipped" begin
        CLM.soil_resistance_read_nl!(soil_resis_method=CLM.SOIL_RESIS_SL_14)
        col, lun, temp, ss, wsb, wdb, mask, bounds = setup_test_data()

        c = 1
        # Leave mask[c] = false
        col.landunit[c] = 1
        col.itype[c] = CLM.ISTSOIL
        lun.itype[1] = CLM.ISTSOIL

        ss.soilresis_col[c] = -999.0
        ss.dsl_col[c] = -999.0

        CLM.calc_soilevap_resis!(col, lun, ss, wsb, wdb, temp, mask, bounds)

        @test ss.soilresis_col[c] == -999.0
        @test ss.dsl_col[c] == -999.0
    end

    @testset "SL14 — frozen soil with ice" begin
        CLM.soil_resistance_read_nl!(soil_resis_method=CLM.SOIL_RESIS_SL_14)
        CLM.surface_resistance_read_params!(d_max=15.0, frac_sat_soil_dsl_init=0.8)
        col, lun, temp, ss, wsb, wdb, mask, bounds = setup_test_data()

        c = 1
        mask[c] = true
        col.landunit[c] = 1
        col.itype[c] = CLM.ISTSOIL
        lun.itype[1] = CLM.ISTSOIL

        col.dz[c, 1 + joff] = 0.02
        # Mostly ice, very little liquid
        wsb.ws.h2osoi_liq_col[c, 1 + joff] = 0.001
        wsb.ws.h2osoi_ice_col[c, 1 + joff] = 8.0  # significant ice

        ss.watsat_col[c, 1] = 0.45
        ss.bsw_col[c, 1] = 5.0
        ss.sucsat_col[c, 1] = 100.0

        temp.t_soisno_col[c, 1 + joff] = 260.0  # below freezing

        CLM.calc_soilevap_resis!(col, lun, ss, wsb, wdb, temp, mask, bounds)

        # With ice filling pores, eff_porosity is reduced → higher resistance
        @test ss.soilresis_col[c] > 0.0
        @test isfinite(ss.soilresis_col[c])
        @test ss.dsl_col[c] >= 0.0
    end

    # Restore module state
    CLM.surface_resistance_ctrl.soil_resis_method = saved_method
    CLM.surface_resistance_params.d_max = saved_d_max
    CLM.surface_resistance_params.frac_sat_soil_dsl_init = saved_frac

    # Restore varpar
    vp.nlevsno = saved_nlevsno
    vp.nlevgrnd = saved_nlevgrnd
    vp.nlevurb = saved_nlevurb
    vp.nlevmaxurbgrnd = saved_nlevmaxurbgrnd
    vp.nlevsoi = saved_nlevsoi
end
