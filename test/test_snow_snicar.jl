@testset "Snow SNICAR" begin

    # ---- Constants ----
    @testset "SNICAR constants" begin
        @test CLM.SNO_NBR_AER == 8
        @test CLM.DO_SNO_AER == true
        @test CLM.DEFAULT_NUMBER_BANDS == 5
        @test CLM.HIGH_NUMBER_BANDS == 480
        @test CLM.IDX_MIE_SNW_MX == 1471
        @test CLM.SNW_RDS_MAX_TBL == 1500
        @test CLM.SNW_RDS_MIN_TBL == 30
        @test CLM.SNW_RDS_MAX == 1500.0
        @test CLM.MIN_SNW ≈ 1.0e-30
        @test CLM.IDX_T_MAX == 11
        @test CLM.IDX_TGRD_MAX == 31
        @test CLM.IDX_RHOS_MAX == 8
    end

    # ---- SnicarParams ----
    @testset "SnicarParams default construction" begin
        p = CLM.SnicarParams()
        @test p.xdrdt == 1.0
        @test p.snw_rds_refrz == 1000.0
        @test p.snw_rds_min == 54.526
        @test p.fresh_snw_rds_max == 204.526
        @test p.C2_liq_Brun89 > 0.0
    end

    # ---- SnicarOpticsData ----
    @testset "SnicarOpticsData default construction" begin
        od = CLM.SnicarOpticsData()
        @test size(od.ss_alb_snw_drc) == (0, 0)
        @test isempty(od.ss_alb_bc_hphil)
    end

    # ---- SnicarAgingData ----
    @testset "SnicarAgingData default construction" begin
        ad = CLM.SnicarAgingData()
        @test size(ad.snowage_tau) == (0, 0, 0)
    end

    # ---- piecewise_linear_interp1d ----
    @testset "piecewise_linear_interp1d" begin
        # Single data point
        @test CLM.piecewise_linear_interp1d([1.0], [5.0], 2.0) == 5.0

        # Linear interpolation between two points
        xd = [0.0, 1.0]
        yd = [0.0, 2.0]
        @test CLM.piecewise_linear_interp1d(xd, yd, 0.5) ≈ 1.0
        @test CLM.piecewise_linear_interp1d(xd, yd, 0.0) ≈ 0.0
        @test CLM.piecewise_linear_interp1d(xd, yd, 1.0) ≈ 2.0

        # Extrapolation below range
        @test CLM.piecewise_linear_interp1d(xd, yd, -0.5) ≈ -1.0

        # Extrapolation above range
        @test CLM.piecewise_linear_interp1d(xd, yd, 1.5) ≈ 3.0

        # Multiple segments
        xd3 = [0.0, 1.0, 3.0]
        yd3 = [0.0, 2.0, 6.0]
        @test CLM.piecewise_linear_interp1d(xd3, yd3, 0.5) ≈ 1.0
        @test CLM.piecewise_linear_interp1d(xd3, yd3, 2.0) ≈ 4.0
        @test CLM.piecewise_linear_interp1d(xd3, yd3, 3.0) ≈ 6.0
    end

    # ---- fresh_snow_radius ----
    @testset "fresh_snow_radius" begin
        p = CLM.SnicarParams(snw_rds_min=54.526, fresh_snw_rds_max=204.526)

        # Very cold: should return minimum
        @test CLM.fresh_snow_radius(200.0; params=p) == p.snw_rds_min

        # At freezing: should return maximum
        @test CLM.fresh_snow_radius(CLM.TFRZ; params=p) == p.fresh_snw_rds_max

        # At tmin (TFRZ - 30): should return minimum
        @test CLM.fresh_snow_radius(CLM.TFRZ - 30.0; params=p) ≈ p.snw_rds_min

        # Mid-range: linear interpolation
        tmid = CLM.TFRZ - 15.0  # halfway
        expected = 0.5 * p.snw_rds_min + 0.5 * p.fresh_snw_rds_max
        @test CLM.fresh_snow_radius(tmid; params=p) ≈ expected

        # When fresh_snw_rds_max <= snw_rds_min, always return min
        p2 = CLM.SnicarParams(snw_rds_min=100.0, fresh_snw_rds_max=50.0)
        @test CLM.fresh_snow_radius(260.0; params=p2) == p2.snw_rds_min
    end

    # ---- snow_optics_init! ----
    @testset "snow_optics_init!" begin
        od = CLM.SnicarOpticsData()
        CLM.snow_optics_init!(od; numrad_snw=5)

        @test size(od.ss_alb_snw_drc) == (CLM.IDX_MIE_SNW_MX, 5)
        @test size(od.ext_cff_mss_snw_dfs) == (CLM.IDX_MIE_SNW_MX, 5)
        @test length(od.ss_alb_bc_hphil) == 5
        @test length(od.ss_alb_dst4) == 5
        @test length(od.flx_wgt_dir) == 5
        @test length(od.flx_wgt_dif) == 5
    end

    # ---- snowage_init! ----
    @testset "snowage_init!" begin
        ad = CLM.SnicarAgingData()
        CLM.snowage_init!(ad)

        @test size(ad.snowage_tau) == (CLM.IDX_RHOS_MAX, CLM.IDX_TGRD_MAX, CLM.IDX_T_MAX)
        @test size(ad.snowage_kappa) == (CLM.IDX_RHOS_MAX, CLM.IDX_TGRD_MAX, CLM.IDX_T_MAX)
        @test size(ad.snowage_drdt0) == (CLM.IDX_RHOS_MAX, CLM.IDX_TGRD_MAX, CLM.IDX_T_MAX)
    end

    # ---- snicar_rt! basic smoke test ----
    @testset "snicar_rt! smoke test (no snow)" begin
        nlevsno = 5
        ncols = 2
        numrad_snw = 5

        # Set up optics with dummy values
        od = CLM.SnicarOpticsData()
        CLM.snow_optics_init!(od; numrad_snw=numrad_snw)
        # Fill with reasonable placeholder values
        fill!(od.ss_alb_snw_drc, 0.9999)
        fill!(od.asm_prm_snw_drc, 0.85)
        fill!(od.ext_cff_mss_snw_drc, 10.0)
        fill!(od.ss_alb_snw_dfs, 0.9999)
        fill!(od.asm_prm_snw_dfs, 0.85)
        fill!(od.ext_cff_mss_snw_dfs, 10.0)
        # Set flux weights so they sum to something nonzero
        od.flx_wgt_dir .= [1.0, 0.5, 0.3, 0.15, 0.05]
        od.flx_wgt_dif .= [1.0, 0.5, 0.3, 0.15, 0.05]

        coszen = [0.5, 0.0]
        h2osno_liq = zeros(ncols, nlevsno)
        h2osno_ice = zeros(ncols, nlevsno)
        h2osno_total = [0.0, 0.0]  # no snow
        snw_rds = fill(100, ncols, nlevsno)
        mss_cnc_aer = zeros(ncols, nlevsno, CLM.SNO_NBR_AER)
        albsfc = fill(0.2, ncols, CLM.NUMRAD)
        snl_vec = [0, 0]  # no snow layers
        frac_sno = [0.0, 0.0]
        albout = zeros(ncols, CLM.NUMRAD)
        flx_abs = zeros(ncols, nlevsno + 1, CLM.NUMRAD)

        # Save current snicar_numrad_snw and restore after
        old_numrad = CLM.varctl.snicar_numrad_snw
        CLM.varctl.snicar_numrad_snw = numrad_snw
        old_shape = CLM.varctl.snicar_snw_shape
        CLM.varctl.snicar_snw_shape = "sphere"

        CLM.snicar_rt!(coszen, 1, h2osno_liq, h2osno_ice, h2osno_total,
                        snw_rds, mss_cnc_aer, albsfc, snl_vec, frac_sno,
                        albout, flx_abs, nlevsno;
                        optics=od)

        # No snow: albedo should be 0 (no sun col 2) or 0 (no snow col 1)
        @test albout[1, CLM.IVIS] == 0.0
        @test albout[1, CLM.INIR] == 0.0
        @test albout[2, CLM.IVIS] == 0.0
        @test albout[2, CLM.INIR] == 0.0

        # Restore
        CLM.varctl.snicar_numrad_snw = old_numrad
        CLM.varctl.snicar_snw_shape = old_shape
    end

    # ---- snicar_rt! with single snow layer ----
    @testset "snicar_rt! with single snow layer" begin
        nlevsno = 5
        ncols = 1
        numrad_snw = 5

        # Set up optics with physically meaningful values
        od = CLM.SnicarOpticsData()
        CLM.snow_optics_init!(od; numrad_snw=numrad_snw)
        # Pure ice grains: high single-scatter albedo, moderate asymmetry
        fill!(od.ss_alb_snw_drc, 0.9999)
        fill!(od.asm_prm_snw_drc, 0.85)
        fill!(od.ext_cff_mss_snw_drc, 10.0)
        fill!(od.ss_alb_snw_dfs, 0.9999)
        fill!(od.asm_prm_snw_dfs, 0.85)
        fill!(od.ext_cff_mss_snw_dfs, 10.0)
        # No aerosol contribution
        fill!(od.ss_alb_bc_hphil, 0.0)
        fill!(od.ss_alb_bc_hphob, 0.0)
        fill!(od.ss_alb_oc_hphil, 0.0)
        fill!(od.ss_alb_oc_hphob, 0.0)
        fill!(od.ss_alb_dst1, 0.0)
        fill!(od.ss_alb_dst2, 0.0)
        fill!(od.ss_alb_dst3, 0.0)
        fill!(od.ss_alb_dst4, 0.0)
        od.flx_wgt_dir .= [1.0, 0.5, 0.3, 0.15, 0.05]
        od.flx_wgt_dif .= [1.0, 0.5, 0.3, 0.15, 0.05]

        coszen = [0.5]
        h2osno_liq = zeros(ncols, nlevsno)
        h2osno_ice = zeros(ncols, nlevsno)
        # Single snow layer: layer index for j=0 is nlevsno=5
        h2osno_ice[1, nlevsno] = 50.0  # 50 kg/m2 ice in top layer
        h2osno_total = [50.0]
        snw_rds = fill(100, ncols, nlevsno)  # 100 micron grain size
        mss_cnc_aer = zeros(ncols, nlevsno, CLM.SNO_NBR_AER)
        albsfc = fill(0.2, ncols, CLM.NUMRAD)
        snl_vec = [-1]  # one snow layer
        frac_sno = [1.0]
        albout = zeros(ncols, CLM.NUMRAD)
        flx_abs = zeros(ncols, nlevsno + 1, CLM.NUMRAD)

        old_numrad = CLM.varctl.snicar_numrad_snw
        CLM.varctl.snicar_numrad_snw = numrad_snw
        old_shape = CLM.varctl.snicar_snw_shape
        CLM.varctl.snicar_snw_shape = "sphere"

        CLM.snicar_rt!(coszen, 1, h2osno_liq, h2osno_ice, h2osno_total,
                        snw_rds, mss_cnc_aer, albsfc, snl_vec, frac_sno,
                        albout, flx_abs, nlevsno;
                        optics=od)

        # With snow and sun, should get nonzero positive albedo between 0 and 1
        @test albout[1, CLM.IVIS] > 0.0
        @test albout[1, CLM.IVIS] <= 1.0
        @test albout[1, CLM.INIR] > 0.0
        @test albout[1, CLM.INIR] <= 1.0

        # VIS albedo should be higher than NIR for clean snow
        @test albout[1, CLM.IVIS] >= albout[1, CLM.INIR]

        # Restore
        CLM.varctl.snicar_numrad_snw = old_numrad
        CLM.varctl.snicar_snw_shape = old_shape
    end

    # ---- snowage_grain! smoke test ----
    @testset "snowage_grain! smoke test" begin
        nlevsno = 5
        ncols = 2
        joff = nlevsno

        # Set up aging tables with dummy values
        ag = CLM.SnicarAgingData()
        CLM.snowage_init!(ag)
        fill!(ag.snowage_tau, 100.0)    # hours
        fill!(ag.snowage_kappa, 1.0)    # unitless
        fill!(ag.snowage_drdt0, 1.0)    # um/hr

        params = CLM.SnicarParams()

        snl_vec = [-2, 0]  # col 1 has 2 snow layers, col 2 has none
        dz = fill(0.1, ncols, nlevsno + 25)  # need enough layers for soil too
        frac_sno = [1.0, 0.0]
        h2osoi_liq = zeros(ncols, nlevsno + 25)
        h2osoi_ice = zeros(ncols, nlevsno + 25)
        t_soisno = fill(260.0, ncols, nlevsno + 25)
        t_grnd = [260.0, 260.0]
        qflx_snow_grnd_col = [0.0, 0.0]
        qflx_snofrz_lyr = zeros(ncols, nlevsno + 25)
        h2osno_no_layers = [0.0, 5.0]
        forc_t = [260.0, 260.0]

        # Initialize snow grain radius and ice content
        snw_rds = fill(0.0, ncols, nlevsno + 25)
        # col 1: 2 snow layers at indices (snl+1+joff) to joff
        # snl=-2, so layers are at Fortran indices -1 and 0, Julia indices 4 and 5
        snw_rds[1, 4] = 100.0
        snw_rds[1, 5] = 100.0
        h2osoi_ice[1, 4] = 20.0
        h2osoi_ice[1, 5] = 30.0

        snw_rds_top = zeros(ncols)
        sno_liq_top = zeros(ncols)
        snot_top = zeros(ncols)
        dTdz_top = zeros(ncols)

        dtime = 1800.0  # 30 minute timestep

        CLM.snowage_grain!(snl_vec, dz, frac_sno, h2osoi_liq, h2osoi_ice,
                            t_soisno, t_grnd, qflx_snow_grnd_col, qflx_snofrz_lyr,
                            h2osno_no_layers, forc_t,
                            snw_rds, snw_rds_top, sno_liq_top,
                            snot_top, dTdz_top,
                            nlevsno, dtime;
                            params=params, aging=ag)

        # Column 1: grain radius should have grown from aging
        @test snw_rds[1, 4] > 100.0  # should have increased
        @test snw_rds[1, 5] > 100.0
        @test snw_rds[1, 4] <= CLM.SNW_RDS_MAX
        @test snw_rds[1, 5] <= CLM.SNW_RDS_MAX

        # Top layer variables should be set for col 1
        @test snot_top[1] ≈ 260.0
        @test snw_rds_top[1] > 0.0

        # Column 2: no snow layers, but h2osno_no_layers > 0
        # should set snw_rds at layer index nlevsno to snw_rds_min
        @test snw_rds[2, nlevsno] ≈ params.snw_rds_min
    end

    # ---- Gaussian quadrature constants ----
    @testset "Gaussian quadrature constants" begin
        @test length(CLM.SNICAR_DIFGAUSPT) == CLM.SNICAR_NGMAX
        @test length(CLM.SNICAR_DIFGAUSWT) == CLM.SNICAR_NGMAX
        @test all(CLM.SNICAR_DIFGAUSPT .> 0.0)
        @test all(CLM.SNICAR_DIFGAUSWT .> 0.0)
    end

    # ---- Aspherical snow grain constants ----
    @testset "Aspherical snow grain constants" begin
        @test length(CLM.SNICAR_G_WVL) == CLM.SNICAR_SEVEN_BANDS + 1
        @test length(CLM.SNICAR_G_B0) == CLM.SNICAR_SEVEN_BANDS
        @test length(CLM.SNICAR_G_B1) == CLM.SNICAR_SEVEN_BANDS
        @test length(CLM.SNICAR_G_B2) == CLM.SNICAR_SEVEN_BANDS
    end

end
