@testset "SurfaceHumidity" begin
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
    nlev_soisno = nlevsno + nlevmaxurbgrnd
    joff = nlevsno  # offset for combined snow+soil indexing

    nc = 6  # number of columns
    nl = 6  # number of landunits (one per column for simplicity)

    # Helper to build minimal test data structures
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

        forc_pbot = fill(101325.0, nc)   # standard pressure
        forc_q = fill(0.005, nc)         # typical specific humidity
        mask_nolakec = falses(nc)        # default: all masked out; each test enables what it needs
        bounds = 1:nc

        return col, lun, temp, soilstate, waterstatebulk, waterdiagbulk,
               forc_pbot, forc_q, mask_nolakec, bounds
    end

    @testset "Soil/crop column basic humidity" begin
        col, lun, temp, ss, wsb, wdb, forc_pbot, forc_q, mask, bounds = setup_test_data()

        # Configure column 1 as soil column
        c = 1
        mask[c] = true
        col.landunit[c] = 1
        col.itype[c] = CLM.ISTSOIL
        lun.itype[1] = CLM.ISTSOIL
        col.snl[c] = 0  # no snow layers

        # Set soil temperature (layer 1)
        temp.t_soisno_col[c, 1 + joff] = 290.0  # warm soil
        temp.t_grnd_col[c] = 290.0
        temp.t_h2osfc_col[c] = 290.0

        # Set soil water (layer 1, combined indexing)
        wsb.ws.h2osoi_liq_col[c, 1 + joff] = 50.0   # kg/m2 liquid
        wsb.ws.h2osoi_ice_col[c, 1 + joff] = 0.0     # no ice

        # Set layer thickness
        col.dz[c, 1 + joff] = 0.02  # 2cm top layer

        # Set soil hydraulic properties (soil-only indexing)
        ss.watsat_col[c, 1] = 0.45
        ss.sucsat_col[c, 1] = 100.0
        ss.bsw_col[c, 1] = 5.0
        ss.smpmin_col[c] = -1.0e8

        # Set fractions
        wdb.frac_sno_eff_col[c] = 0.0
        wdb.frac_h2osfc_col[c] = 0.0

        CLM.calculate_surface_humidity!(col, lun, temp, ss, wsb, wdb,
            forc_pbot, forc_q, mask, bounds)

        # soilalpha should be set (qred = hr since no snow/h2osfc)
        @test !isnan(ss.soilalpha_col[c])
        @test ss.soilalpha_col[c] > 0.0
        @test ss.soilalpha_col[c] <= 1.0

        # qg_soil should be set
        @test !isnan(wdb.qg_soil_col[c])
        @test wdb.qg_soil_col[c] > 0.0

        # With snl=0, qg_snow = qg_soil
        @test wdb.qg_snow_col[c] == wdb.qg_soil_col[c]

        # With frac_h2osfc=0, qg_h2osfc = qg_soil
        @test wdb.qg_h2osfc_col[c] == wdb.qg_soil_col[c]

        # qg should equal qg_soil when no snow and no h2osfc
        @test wdb.qg_col[c] ≈ wdb.qg_soil_col[c]

        # dqgdT should be set
        @test !isnan(wdb.dqgdT_col[c])
    end

    @testset "Soil column with snow layers" begin
        col, lun, temp, ss, wsb, wdb, forc_pbot, forc_q, mask, bounds = setup_test_data()

        c = 1
        mask[c] = true
        col.landunit[c] = 1
        col.itype[c] = CLM.ISTSOIL
        lun.itype[1] = CLM.ISTSOIL
        col.snl[c] = -2  # 2 snow layers

        # Set temperatures
        temp.t_soisno_col[c, 1 + joff] = 280.0  # soil layer 1
        # Snow layer: snl+1 = -1, index = -1 + 1 + joff = joff
        temp.t_soisno_col[c, joff] = 268.0  # top snow layer (j=0)
        temp.t_soisno_col[c, joff - 1] = 265.0  # snow layer (j=-1)
        temp.t_grnd_col[c] = 268.0
        temp.t_h2osfc_col[c] = 280.0

        # Soil water
        wsb.ws.h2osoi_liq_col[c, 1 + joff] = 40.0
        wsb.ws.h2osoi_ice_col[c, 1 + joff] = 5.0
        col.dz[c, 1 + joff] = 0.02

        # Soil properties
        ss.watsat_col[c, 1] = 0.45
        ss.sucsat_col[c, 1] = 100.0
        ss.bsw_col[c, 1] = 5.0
        ss.smpmin_col[c] = -1.0e8

        # Fractions
        wdb.frac_sno_eff_col[c] = 0.5
        wdb.frac_h2osfc_col[c] = 0.0

        CLM.calculate_surface_humidity!(col, lun, temp, ss, wsb, wdb,
            forc_pbot, forc_q, mask, bounds)

        # With snow: qg_snow should be saturated humidity at snow temperature
        @test !isnan(wdb.qg_snow_col[c])
        @test wdb.qg_snow_col[c] > 0.0

        # qg_soil should differ from qg_snow (different temperatures)
        @test wdb.qg_soil_col[c] != wdb.qg_snow_col[c]

        # qg should be weighted average
        expected_qg = wdb.frac_sno_eff_col[c] * wdb.qg_snow_col[c] +
                      (1.0 - wdb.frac_sno_eff_col[c] - wdb.frac_h2osfc_col[c]) * wdb.qg_soil_col[c] +
                      wdb.frac_h2osfc_col[c] * wdb.qg_h2osfc_col[c]
        @test wdb.qg_col[c] ≈ expected_qg
    end

    @testset "Soil column with surface water" begin
        col, lun, temp, ss, wsb, wdb, forc_pbot, forc_q, mask, bounds = setup_test_data()

        c = 1
        mask[c] = true
        col.landunit[c] = 1
        col.itype[c] = CLM.ISTSOIL
        lun.itype[1] = CLM.ISTSOIL
        col.snl[c] = 0

        temp.t_soisno_col[c, 1 + joff] = 290.0
        temp.t_grnd_col[c] = 290.0
        temp.t_h2osfc_col[c] = 292.0  # slightly warmer surface water

        wsb.ws.h2osoi_liq_col[c, 1 + joff] = 50.0
        wsb.ws.h2osoi_ice_col[c, 1 + joff] = 0.0
        col.dz[c, 1 + joff] = 0.02

        ss.watsat_col[c, 1] = 0.45
        ss.sucsat_col[c, 1] = 100.0
        ss.bsw_col[c, 1] = 5.0
        ss.smpmin_col[c] = -1.0e8

        wdb.frac_sno_eff_col[c] = 0.0
        wdb.frac_h2osfc_col[c] = 0.3  # 30% surface water

        CLM.calculate_surface_humidity!(col, lun, temp, ss, wsb, wdb,
            forc_pbot, forc_q, mask, bounds)

        # qg_h2osfc should be saturated at t_h2osfc temperature
        @test !isnan(wdb.qg_h2osfc_col[c])
        @test wdb.qg_h2osfc_col[c] > 0.0
        # qg_h2osfc should differ from qg_soil due to different temperatures
        @test wdb.qg_h2osfc_col[c] != wdb.qg_soil_col[c]

        # qg should include h2osfc contribution
        expected_qg = (1.0 - wdb.frac_h2osfc_col[c]) * wdb.qg_soil_col[c] +
                      wdb.frac_h2osfc_col[c] * wdb.qg_h2osfc_col[c]
        @test wdb.qg_col[c] ≈ expected_qg
    end

    @testset "Wet landunit" begin
        col, lun, temp, ss, wsb, wdb, forc_pbot, forc_q, mask, bounds = setup_test_data()

        c = 1
        mask[c] = true
        col.landunit[c] = 1
        col.itype[c] = CLM.ISTWET
        lun.itype[1] = CLM.ISTWET
        col.snl[c] = 0

        temp.t_grnd_col[c] = 285.0

        CLM.calculate_surface_humidity!(col, lun, temp, ss, wsb, wdb,
            forc_pbot, forc_q, mask, bounds)

        # For wet landunit: soilalpha = spval
        @test ss.soilalpha_col[c] == CLM.SPVAL

        # qred = 1.0, so qg = qsatg
        (qs_expected, _, _, _) = CLM.qsat(285.0, 101325.0)
        @test wdb.qg_col[c] ≈ qs_expected

        # All humidity components should be equal
        @test wdb.qg_snow_col[c] == wdb.qg_col[c]
        @test wdb.qg_soil_col[c] == wdb.qg_col[c]
        @test wdb.qg_h2osfc_col[c] == wdb.qg_col[c]
    end

    @testset "Ice landunit" begin
        col, lun, temp, ss, wsb, wdb, forc_pbot, forc_q, mask, bounds = setup_test_data()

        c = 1
        mask[c] = true
        col.landunit[c] = 1
        col.itype[c] = CLM.ISTICE
        lun.itype[1] = CLM.ISTICE
        col.snl[c] = 0

        temp.t_grnd_col[c] = 260.0

        CLM.calculate_surface_humidity!(col, lun, temp, ss, wsb, wdb,
            forc_pbot, forc_q, mask, bounds)

        @test ss.soilalpha_col[c] == CLM.SPVAL

        (qs_expected, _, _, _) = CLM.qsat(260.0, 101325.0)
        @test wdb.qg_col[c] ≈ qs_expected
    end

    @testset "Urban sunwall (qred=0)" begin
        col, lun, temp, ss, wsb, wdb, forc_pbot, forc_q, mask, bounds = setup_test_data()

        c = 1
        mask[c] = true
        col.landunit[c] = 1
        col.itype[c] = CLM.ICOL_SUNWALL
        lun.itype[1] = CLM.ISTURB_TBD  # urban landunit (not wet or ice)
        col.snl[c] = 0

        temp.t_grnd_col[c] = 300.0
        # Set forc_q=0 so the dew cap (qsatg > forc_q && forc_q > qred*qsatg) doesn't trigger
        forc_q[c] = 0.0

        CLM.calculate_surface_humidity!(col, lun, temp, ss, wsb, wdb,
            forc_pbot, forc_q, mask, bounds)

        @test ss.soilalpha_u_col[c] == CLM.SPVAL

        # qred=0 → qg=0
        @test wdb.qg_col[c] ≈ 0.0 atol=1e-15
        @test wdb.dqgdT_col[c] ≈ 0.0 atol=1e-15
    end

    @testset "Urban roof (qred=1)" begin
        col, lun, temp, ss, wsb, wdb, forc_pbot, forc_q, mask, bounds = setup_test_data()

        c = 1
        mask[c] = true
        col.landunit[c] = 1
        col.itype[c] = CLM.ICOL_ROOF
        lun.itype[1] = CLM.ISTURB_HD  # urban
        col.snl[c] = 0

        temp.t_grnd_col[c] = 300.0

        CLM.calculate_surface_humidity!(col, lun, temp, ss, wsb, wdb,
            forc_pbot, forc_q, mask, bounds)

        @test ss.soilalpha_u_col[c] == CLM.SPVAL

        # qred=1 → qg = qsatg
        (qs_expected, _, _, _) = CLM.qsat(300.0, 101325.0)
        @test wdb.qg_col[c] ≈ qs_expected
    end

    @testset "Urban pervious road" begin
        col, lun, temp, ss, wsb, wdb, forc_pbot, forc_q, mask, bounds = setup_test_data()

        c = 1
        mask[c] = true
        col.landunit[c] = 1
        col.itype[c] = CLM.ICOL_ROAD_PERV
        lun.itype[1] = CLM.ISTURB_MD  # urban
        col.snl[c] = 0

        wdb.frac_sno_eff_col[c] = 0.0

        # Set soil temperatures and water for all levels
        for j in 1:nlevgrnd
            temp.t_soisno_col[c, j + joff] = 290.0
            col.dz[c, j + joff] = 0.1
            wsb.ws.h2osoi_liq_col[c, j + joff] = 10.0
            wsb.ws.h2osoi_ice_col[c, j + joff] = 0.0
            ss.watsat_col[c, j] = 0.45
            ss.watdry_col[c, j] = 0.05
            ss.watopt_col[c, j] = 0.35
            ss.rootfr_road_perv_col[c, j] = 1.0 / nlevgrnd
        end

        temp.t_grnd_col[c] = 290.0

        CLM.calculate_surface_humidity!(col, lun, temp, ss, wsb, wdb,
            forc_pbot, forc_q, mask, bounds)

        # soilalpha_u should be set for pervious road
        @test !isnan(ss.soilalpha_u_col[c])
        @test ss.soilalpha_u_col[c] > 0.0
        @test ss.soilalpha_u_col[c] <= 1.0

        # rootr_road_perv should be normalized (sum to 1.0)
        rootr_sum = sum(ss.rootr_road_perv_col[c, j] for j in 1:nlevgrnd)
        @test rootr_sum ≈ 1.0 atol=1e-10

        # qg should be set
        @test !isnan(wdb.qg_col[c])
        @test wdb.qg_col[c] > 0.0
    end

    @testset "Masked columns are skipped" begin
        col, lun, temp, ss, wsb, wdb, forc_pbot, forc_q, mask, bounds = setup_test_data()

        # Configure column 1 as soil but leave it masked out (default is falses)
        c = 1
        col.landunit[c] = 1
        col.itype[c] = CLM.ISTSOIL
        lun.itype[1] = CLM.ISTSOIL
        col.snl[c] = 0

        # Set outputs to known values
        wdb.qg_col[c] = -999.0
        wdb.qg_soil_col[c] = -999.0

        CLM.calculate_surface_humidity!(col, lun, temp, ss, wsb, wdb,
            forc_pbot, forc_q, mask, bounds)

        # Masked column should not be modified
        @test wdb.qg_col[c] == -999.0
        @test wdb.qg_soil_col[c] == -999.0
    end

    @testset "Dew formation capping (qsatg > forc_q > hr*qsatg)" begin
        col, lun, temp, ss, wsb, wdb, forc_pbot, forc_q, mask, bounds = setup_test_data()

        c = 1
        mask[c] = true
        col.landunit[c] = 1
        col.itype[c] = CLM.ISTSOIL
        lun.itype[1] = CLM.ISTSOIL
        col.snl[c] = 0

        # Use a dry soil to get a small hr value
        temp.t_soisno_col[c, 1 + joff] = 300.0
        temp.t_grnd_col[c] = 300.0
        temp.t_h2osfc_col[c] = 300.0

        # Very dry soil → small hr
        wsb.ws.h2osoi_liq_col[c, 1 + joff] = 0.5
        wsb.ws.h2osoi_ice_col[c, 1 + joff] = 0.0
        col.dz[c, 1 + joff] = 0.02

        ss.watsat_col[c, 1] = 0.45
        ss.sucsat_col[c, 1] = 200.0
        ss.bsw_col[c, 1] = 8.0
        ss.smpmin_col[c] = -1.0e8

        wdb.frac_sno_eff_col[c] = 0.0
        wdb.frac_h2osfc_col[c] = 0.0

        # Compute qsatg to find the right forc_q
        (qs, _, _, _) = CLM.qsat(300.0, 101325.0)

        # Set forc_q between hr*qsatg and qsatg to trigger the dew cap
        # We need a moderate forc_q that will be between hr*qsatg and qsatg
        forc_q[c] = qs * 0.5  # half of saturation

        CLM.calculate_surface_humidity!(col, lun, temp, ss, wsb, wdb,
            forc_pbot, forc_q, mask, bounds)

        # If dew cap was triggered, qg_soil = forc_q (note: this depends on hr value)
        # We can at least verify the output is reasonable
        @test wdb.qg_soil_col[c] > 0.0
        @test !isnan(wdb.qg_col[c])
    end

    @testset "Wet/ice else-branch dew capping" begin
        col, lun, temp, ss, wsb, wdb, forc_pbot, forc_q, mask, bounds = setup_test_data()

        c = 1
        mask[c] = true
        col.landunit[c] = 1
        col.itype[c] = CLM.ISTWET
        lun.itype[1] = CLM.ISTWET
        col.snl[c] = 0

        temp.t_grnd_col[c] = 300.0

        # Compute saturation humidity
        (qs, _, _, _) = CLM.qsat(300.0, 101325.0)

        # Set forc_q between qred*qsatg and qsatg (qred=1 for wet, so this won't trigger)
        # For wet, qred=1, so forc_q must be between 1*qs and qs — impossible
        # So qg should just be qs
        forc_q[c] = qs * 0.999

        CLM.calculate_surface_humidity!(col, lun, temp, ss, wsb, wdb,
            forc_pbot, forc_q, mask, bounds)

        @test wdb.qg_col[c] ≈ qs
    end

    @testset "Frozen pervious road (fac=0)" begin
        col, lun, temp, ss, wsb, wdb, forc_pbot, forc_q, mask, bounds = setup_test_data()

        c = 1
        mask[c] = true
        col.landunit[c] = 1
        col.itype[c] = CLM.ICOL_ROAD_PERV
        lun.itype[1] = CLM.ISTURB_TBD
        col.snl[c] = 0

        wdb.frac_sno_eff_col[c] = 0.0

        # Set frozen soil temperatures
        for j in 1:nlevgrnd
            temp.t_soisno_col[c, j + joff] = 260.0  # frozen
            col.dz[c, j + joff] = 0.1
            wsb.ws.h2osoi_liq_col[c, j + joff] = 5.0
            wsb.ws.h2osoi_ice_col[c, j + joff] = 50.0
            ss.watsat_col[c, j] = 0.45
            ss.watdry_col[c, j] = 0.05
            ss.watopt_col[c, j] = 0.35
            ss.rootfr_road_perv_col[c, j] = 1.0 / nlevgrnd
        end

        temp.t_grnd_col[c] = 260.0

        CLM.calculate_surface_humidity!(col, lun, temp, ss, wsb, wdb,
            forc_pbot, forc_q, mask, bounds)

        # All frozen → fac=0 for all layers → hr_road_perv=0, qred=frac_sno_eff=0
        @test ss.soilalpha_u_col[c] ≈ 0.0 atol=1e-15

        # qg should be 0 since qred=0
        @test wdb.qg_col[c] ≈ 0.0 atol=1e-15
    end

    # Restore varpar
    vp.nlevsno = saved_nlevsno
    vp.nlevgrnd = saved_nlevgrnd
    vp.nlevurb = saved_nlevurb
    vp.nlevmaxurbgrnd = saved_nlevmaxurbgrnd
    vp.nlevsoi = saved_nlevsoi
end
