@testset "LUNA" begin

    # =====================================================================
    # Constants
    # =====================================================================
    @testset "LUNA Constants" begin
        @test CLM.LUNA_Cv ≈ 1.2e-5 * 3600.0
        @test CLM.LUNA_Fc25 == 294.2
        @test CLM.LUNA_Fj25 == 1257.0
        @test CLM.LUNA_NUEr25 == 33.69
        @test CLM.LUNA_Cb == 1.78
        @test CLM.LUNA_O2ref == 209460.0
        @test CLM.LUNA_CO2ref == 380.0
        @test CLM.LUNA_forc_pbot_ref == 101325.0
        @test CLM.LUNA_Q10Enz == 2.0
        @test CLM.LUNA_NMCp25 == 0.715
        @test CLM.LUNA_Trange1 == 5.0
        @test CLM.LUNA_Trange2 == 42.0
        @test CLM.LUNA_SNC == 0.004
        @test CLM.LUNA_mp == 9.0
        @test CLM.LUNA_PARLowLim == 200.0
    end

    # =====================================================================
    # LunaParamsData
    # =====================================================================
    @testset "LunaParamsData init/clean" begin
        lp = CLM.LunaParamsData()
        CLM.luna_params_init!(lp, CLM.MXPFT)
        @test length(lp.jmaxb0) == CLM.MXPFT + 1
        @test length(lp.jmaxb1) == CLM.MXPFT + 1
        @test length(lp.wc2wjb0) == CLM.MXPFT + 1
        @test all(isnan, lp.jmaxb0)
        @test all(isnan, lp.jmaxb1)

        CLM.luna_params_clean!(lp)
        @test isempty(lp.jmaxb0)
        @test isempty(lp.jmaxb1)
        @test isempty(lp.wc2wjb0)
    end

    # Helper to create typical LUNA params for testing
    function make_test_luna_params()
        lp = CLM.LunaParamsData()
        # Typical CLM parameter values (after unit conversion in readParams)
        lp.cp25_yr2000 = 4.275            # CO2 compensation point (~4.275 Pa)
        lp.kc25_coef = 40.49              # Michaelis-Menten for CO2
        lp.ko25_coef = 27840.0            # Michaelis-Menten for O2
        lp.luna_theta_cj = 0.9
        lp.enzyme_turnover_daily = 0.1
        lp.relhExp = 6.0
        lp.minrelh = 0.25
        return lp
    end

    # =====================================================================
    # Temperature response functions
    # =====================================================================
    @testset "vcmx_t_kattge" begin
        # At 25°C tgrow and tleaf, should be close to 1.0
        val = CLM.vcmx_t_kattge(25.0, 25.0)
        @test val ≈ 1.0 atol=1e-10
        @test isfinite(val)

        # Temperature increases -> value increases (up to optimum)
        val30 = CLM.vcmx_t_kattge(25.0, 30.0)
        @test val30 > val
        @test isfinite(val30)

        # Low temperature -> value decreases
        val10 = CLM.vcmx_t_kattge(25.0, 10.0)
        @test val10 < val
        @test isfinite(val10)
    end

    @testset "jmx_t_kattge" begin
        val = CLM.jmx_t_kattge(25.0, 25.0)
        @test val ≈ 1.0 atol=1e-10
        @test isfinite(val)

        val30 = CLM.jmx_t_kattge(25.0, 30.0)
        @test val30 > val
        @test isfinite(val30)
    end

    @testset "vcmx_t_leuning" begin
        val = CLM.vcmx_t_leuning(25.0, 25.0)
        @test val ≈ 1.0 atol=1e-10
        @test isfinite(val)
    end

    @testset "jmx_t_leuning" begin
        val = CLM.jmx_t_leuning(25.0, 25.0)
        @test val ≈ 1.0 atol=1e-10
        @test isfinite(val)
    end

    @testset "resp_t_bernacchi" begin
        val25 = CLM.resp_t_bernacchi(25.0)
        @test val25 > 0.0
        @test isfinite(val25)

        # Higher temperature → higher respiration
        val35 = CLM.resp_t_bernacchi(35.0)
        @test val35 > val25
    end

    # =====================================================================
    # Quadratic solver
    # =====================================================================
    @testset "quadratic_luna" begin
        # x^2 - 3x + 2 = 0 => roots 1, 2
        r1, r2 = CLM.quadratic_luna(1.0, -3.0, 2.0)
        @test r1 ≈ 2.0 atol=1e-10
        @test r2 ≈ 1.0 atol=1e-10

        # a=0 => special case, returns 1e36
        r1, r2 = CLM.quadratic_luna(0.0, -3.0, 2.0)
        @test r1 == 1.0e36
        @test r2 == 1.0e36

        # x^2 - 5x + 6 = 0 => roots 2, 3
        r1, r2 = CLM.quadratic_luna(1.0, -5.0, 6.0)
        @test r1 ≈ 3.0 atol=1e-10
        @test r2 ≈ 2.0 atol=1e-10
    end

    # =====================================================================
    # NUE reference
    # =====================================================================
    @testset "nue_ref" begin
        lp = make_test_luna_params()
        NUEjref, NUEcref, Kj2Kcref = CLM.nue_ref(lp)
        @test NUEjref > 0.0
        @test NUEcref > 0.0
        @test Kj2Kcref > 0.0
        @test isfinite(NUEjref)
        @test isfinite(NUEcref)
        @test isfinite(Kj2Kcref)
    end

    # =====================================================================
    # NUE calculation
    # =====================================================================
    @testset "nue_calc" begin
        lp = make_test_luna_params()
        O2a = 21232.0  # Pa
        ci = 27.0      # Pa
        tgrow = 25.0
        tleaf = 25.0

        NUEj, NUEc, Kj2Kc = CLM.nue_calc(O2a, ci, tgrow, tleaf, lp)
        @test NUEj > 0.0
        @test NUEc > 0.0
        @test Kj2Kc > 0.0
        @test isfinite(NUEj)
        @test isfinite(NUEc)
    end

    # =====================================================================
    # Photosynthesis LUNA
    # =====================================================================
    @testset "photosynthesis_luna!" begin
        lp = make_test_luna_params()

        forc_pbot = 101325.0
        tleafd = 25.0
        relh = 0.6
        CO2a = 38.0    # Pa
        O2a = 21232.0  # Pa
        rb = 25.0
        Vcmax = 50.0
        JmeanL = 40.0

        ci, Kc, Kj, A = CLM.photosynthesis_luna!(forc_pbot, tleafd, relh, CO2a, O2a, rb, Vcmax, JmeanL, lp)
        @test ci > 0.0
        @test Kc >= 0.0
        @test Kj >= 0.0
        @test A >= 0.0
        @test isfinite(ci)
        @test isfinite(A)
    end

    # =====================================================================
    # is_time_to_run_luna
    # =====================================================================
    @testset "is_time_to_run_luna" begin
        @test CLM.is_time_to_run_luna(true) == true
        @test CLM.is_time_to_run_luna(false) == false
    end

    # =====================================================================
    # Nitrogen allocation
    # =====================================================================
    @testset "nitrogen_allocation!" begin
        lp = make_test_luna_params()

        FNCa = 2.0
        forc_pbot10 = 101325.0
        relh10 = 0.6
        CO2a10 = 38.0
        O2a10 = 21232.0
        PARi10 = 500.0
        PARimx10 = 800.0
        rb10 = 25.0
        hourpd = 12.0
        tair10 = 25.0
        tleafd10 = 28.0
        tleafn10 = 18.0
        jmaxb0_val = 0.03
        jmaxb1_val = 0.22
        wc2wjb0_val = 0.8
        PNlcold = 0.15
        PNetold = 0.0
        PNrespold = 0.0
        PNcbold = 0.0
        dayl_factor = 1.0
        o3coefjmax = 1.0

        PNstoreopt, PNlcopt, PNetopt, PNrespopt, PNcbopt =
            CLM.nitrogen_allocation!(FNCa, forc_pbot10, relh10, CO2a10, O2a10,
                PARi10, PARimx10, rb10, hourpd, tair10, tleafd10, tleafn10,
                jmaxb0_val, jmaxb1_val, wc2wjb0_val, PNlcold, PNetold,
                PNrespold, PNcbold, dayl_factor, o3coefjmax, lp)

        # All fractions should be finite and non-negative
        @test isfinite(PNstoreopt)
        @test isfinite(PNlcopt)
        @test isfinite(PNetopt)
        @test isfinite(PNrespopt)
        @test isfinite(PNcbopt)
        @test PNlcopt >= 0.0
        @test PNetopt >= 0.0
        @test PNrespopt >= 0.0
        @test PNcbopt >= 0.0

        # Fractions should roughly sum to 1
        total = PNstoreopt + PNlcopt + PNetopt + PNrespopt + PNcbopt
        @test total ≈ 1.0 atol=0.1
    end

    # =====================================================================
    # Clear 24hr climate
    # =====================================================================
    @testset "clear24_climate_luna!" begin
        CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
        np = 3
        nc = 5
        nl = 2
        ng = 1

        solarabs = CLM.SolarAbsorbedData()
        CLM.solarabs_init!(solarabs, np, nl; use_luna=true)
        photosyns = CLM.PhotosynthesisData()
        CLM.photosynthesis_data_init!(photosyns, np; use_luna=true)
        temperature = CLM.TemperatureData()
        CLM.temperature_init!(temperature, np, nc, nl, ng)
        patchdata = CLM.PatchData()
        CLM.patch_init!(patchdata, np)

        # Set some non-zero values
        temperature.t_veg_day_patch[1] = 300.0
        solarabs.par24d_z_patch[1, 1] = 100.0
        photosyns.fpsn24_patch[1] = 5.0
        temperature.ndaysteps_patch[1] = 10

        mask = trues(np)
        bounds = 1:np

        CLM.clear24_climate_luna!(solarabs, photosyns, temperature, patchdata, mask, bounds)

        @test temperature.t_veg_day_patch[1] == 0.0
        @test temperature.t_veg_night_patch[1] == 0.0
        @test solarabs.par24d_z_patch[1, 1] == 0.0
        @test solarabs.par24x_z_patch[1, 1] == 0.0
        @test photosyns.fpsn24_patch[1] == 0.0
        @test temperature.ndaysteps_patch[1] == 0
        @test temperature.nnightsteps_patch[1] == 0
    end

    # =====================================================================
    # Acc24 climate
    # =====================================================================
    @testset "acc24_climate_luna!" begin
        CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
        np = 2
        nc = 5
        nl = 2
        ng = 1

        canopystate = CLM.CanopyStateData()
        CLM.canopystate_init!(canopystate, np)
        photosyns = CLM.PhotosynthesisData()
        CLM.photosynthesis_data_init!(photosyns, np; use_luna=true)
        surfalb = CLM.SurfaceAlbedoData()
        CLM.surfalb_init!(surfalb, np, nc, ng)
        solarabs = CLM.SolarAbsorbedData()
        CLM.solarabs_init!(solarabs, np, nl; use_luna=true)
        temperature = CLM.TemperatureData()
        CLM.temperature_init!(temperature, np, nc, nl, ng)
        patchdata = CLM.PatchData()
        CLM.patch_init!(patchdata, np)

        # Set valid state (not first day)
        temperature.t_veg_day_patch .= 0.0
        temperature.t_veg_night_patch .= 0.0
        temperature.t_veg_patch .= 300.0
        temperature.ndaysteps_patch .= 0
        temperature.nnightsteps_patch .= 0
        solarabs.sabv_patch .= 100.0  # daytime (sabv > 0)
        surfalb.nrad_patch .= 1
        canopystate.laisun_z_patch .= 1.0
        canopystate.laisha_z_patch .= 0.5
        solarabs.parsun_z_patch .= 200.0
        solarabs.parsha_z_patch .= 50.0
        solarabs.par24d_z_patch .= 0.0
        solarabs.par24x_z_patch .= 0.0
        photosyns.fpsn_patch .= 10.0
        photosyns.fpsn24_patch .= 0.0

        mask = trues(np)
        bounds = 1:np
        dtime = 1800.0

        CLM.acc24_climate_luna!(canopystate, photosyns, surfalb, solarabs,
            temperature, patchdata, mask, bounds, dtime)

        # Daytime: t_veg_day accumulated, ndaysteps increased
        @test temperature.t_veg_day_patch[1] ≈ 300.0
        @test temperature.ndaysteps_patch[1] == 1
        # PAR accumulated
        @test solarabs.par24d_z_patch[1, 1] > 0.0
        # GPP accumulated
        @test photosyns.fpsn24_patch[1] > 0.0
    end

    # =====================================================================
    # Nitrogen investments
    # =====================================================================
    @testset "nitrogen_investments!" begin
        lp = make_test_luna_params()

        FNCa = 2.0
        Nlc = 0.3
        Fc = CLM.vcmx_t_kattge(25.0, 28.0) * CLM.LUNA_Fc25
        Fj = CLM.jmx_t_kattge(25.0, 28.0) * CLM.LUNA_Fj25
        NUEr = CLM.LUNA_Cv * CLM.LUNA_NUEr25 * (CLM.resp_t_bernacchi(28.0) * 12.0 + CLM.resp_t_bernacchi(18.0) * 12.0)
        NUEj, NUEc, Kj2Kc = CLM.nue_calc(21232.0, 26.6, 25.0, 28.0, lp)
        NUEjref, NUEcref, _ = CLM.nue_ref(lp)
        JmaxCoef = 0.22 * 1.0 * (1.0 - exp(-lp.relhExp * max(0.6 - lp.minrelh, 0.0) / (1.0 - lp.minrelh)))

        Vcmax, Jmax, JmeanL, JmaxL, Net, Ncb, Nresp, PSN, RESP, Kc_out, Kj_out, ci_out =
            CLM.nitrogen_investments!(0, FNCa, Nlc, 101325.0, 0.6, 38.0, 21232.0,
                500.0, 800.0, 25.0, 12.0, 25.0, 28.0, 18.0,
                Kj2Kc, JmaxCoef, Fc, Fj, NUEc, NUEj, NUEcref, NUEjref, NUEr, 1.0,
                0.03, 0.8, 0.0, 0.0, 26.6, lp)

        @test Vcmax > 0.0
        @test Jmax > 0.0
        @test JmeanL > 0.0
        @test JmaxL > 0.0
        @test Net > 0.0
        @test Ncb > 0.0
        @test Nresp > 0.0
        @test PSN > 0.0
        @test RESP > 0.0
        @test isfinite(Vcmax)
        @test isfinite(PSN)
    end

end
