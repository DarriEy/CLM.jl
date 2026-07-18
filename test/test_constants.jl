@testset "Constants & Parameters" begin

    @testset "Landunit constants" begin
        @test CLM.ISTSOIL == 1
        @test CLM.ISTCROP == 2
        @test CLM.ISTICE == 4
        @test CLM.ISTDLAK == 5
        @test CLM.ISTURB_TBD == 7
        @test CLM.ISTURB_HD == 8
        @test CLM.ISTURB_MD == 9
        @test CLM.MAX_LUNIT == 9
        @test CLM.NUMURBL == 3
        @test length(CLM.LANDUNIT_NAMES) == CLM.MAX_LUNIT
    end

    @testset "landunit_is_special" begin
        @test CLM.landunit_is_special(CLM.ISTSOIL) == false
        @test CLM.landunit_is_special(CLM.ISTCROP) == false
        @test CLM.landunit_is_special(CLM.ISTOCN) == true
        @test CLM.landunit_is_special(CLM.ISTICE) == true
        @test CLM.landunit_is_special(CLM.ISTDLAK) == true
        @test CLM.landunit_is_special(CLM.ISTWET) == true
        @test CLM.landunit_is_special(CLM.ISTURB_TBD) == true
        @test CLM.landunit_is_special(CLM.ISTURB_HD) == true
        @test CLM.landunit_is_special(CLM.ISTURB_MD) == true
    end

    @testset "Column constants" begin
        @test CLM.ICOL_ROOF == 71
        @test CLM.ICOL_SUNWALL == 72
        @test CLM.ICOL_SHADEWALL == 73
        @test CLM.ICOL_ROAD_IMPERV == 74
        @test CLM.ICOL_ROAD_PERV == 75
    end

    @testset "is_hydrologically_active" begin
        # Soil and crop are always active
        @test CLM.is_hydrologically_active(1, CLM.ISTSOIL) == true
        @test CLM.is_hydrologically_active(1, CLM.ISTCROP) == true
        # Ice is active
        @test CLM.is_hydrologically_active(1, CLM.ISTICE) == true
        # Lake is active
        @test CLM.is_hydrologically_active(1, CLM.ISTDLAK) == true
        # Wetland is active
        @test CLM.is_hydrologically_active(1, CLM.ISTWET) == true
        # Urban: only pervious road is active
        @test CLM.is_hydrologically_active(CLM.ICOL_ROAD_PERV, CLM.ISTURB_TBD) == true
        @test CLM.is_hydrologically_active(CLM.ICOL_ROOF, CLM.ISTURB_TBD) == false
        @test CLM.is_hydrologically_active(CLM.ICOL_SUNWALL, CLM.ISTURB_HD) == false
        @test CLM.is_hydrologically_active(CLM.ICOL_ROAD_IMPERV, CLM.ISTURB_MD) == false
    end

    @testset "Physical constants" begin
        @test CLM.TFRZ ≈ 273.15
        @test CLM.GRAV ≈ 9.80616
        @test CLM.DENH2O ≈ 1000.0
        @test CLM.DENICE ≈ 917.0
        @test CLM.CPAIR ≈ 1004.64
        @test CLM.HVAP ≈ 2.501e6
        @test CLM.HFUS ≈ 3.337e5
        @test CLM.RGAS_KMOL == 6.02214e26 * 1.38065e-23
        @test CLM.RAIR == CLM.RGAS_KMOL / 28.966
        @test CLM.RWAT == CLM.RGAS_KMOL / 18.016
        @test CLM.RGAS == CLM.RGAS_KMOL * 1.0e-3
        @test CLM.RPI ≈ π atol=1e-10
        @test CLM.SECSPDAY ≈ 86400.0
        @test CLM.VKC ≈ 0.4
    end

    @testset "Varctl defaults" begin
        vc = CLM.VarCtl()
        @test vc.nsrest == CLM.IUNDEF
        @test vc.use_crop == false
        @test vc.use_fates == false
        @test vc.use_cn == false
        @test vc.use_lch4 == true
        @test vc.use_nitrif_denitrif == true
        @test vc.co2_ppmv ≈ 355.0
        @test vc.spinup_state == 0
        @test vc.irrigate == false
    end

    @testset "VarPar defaults" begin
        vp = CLM.VarPar()
        @test vp.nlevsno == -1
        @test vp.maxpatch_urb == 5
    end

    @testset "Radiation band indices" begin
        @test CLM.NUMRAD == 2
        @test CLM.IVIS == 1
        @test CLM.INIR == 2
    end

    @testset "Gas constants arrays" begin
        @test size(CLM.S_CON) == (3, 4)
        @test size(CLM.D_CON_W) == (3, 3)
        @test size(CLM.D_CON_G) == (3, 2)
        @test length(CLM.C_H_INV) == CLM.NGASES
        @test length(CLM.KH_THETA) == CLM.NGASES

        # VALUE checks (not just shape): these are the raw polynomial coefficients
        # from clm_varcon.F90 d_con_w / d_con_g. The consumers (ch4Mod / nitrif)
        # apply the 1e-9 (water) and 1e-4 (gas) scale in-formula and use ALL
        # coefficients — so a mis-ported table (e.g. a 1.0e3 placeholder for the
        # higher-order term) silently explodes the saturated-column CH4 diffusivity
        # by ~1e5x. Guard the exact Fortran values.
        @test CLM.D_CON_W ≈ [0.9798 0.02986 0.0004381;
                             1.172   0.03443 0.0005048;
                             0.939   0.02671 0.0004095]
        @test CLM.D_CON_G ≈ [0.1875 0.0013;
                             0.1759 0.00117;
                             0.1325 0.0009]

        # Henry's-law constants (CH4, O2, CO2) from clm_varcon.F90:
        #   data c_h_inv(1:3)  /1600., 1500., 2400./    (K) temp-dependence const
        #   data kh_theta(1:3) /714.29, 769.23, 29.4/   (L.atm/mol) Henry const @298K
        # A prior port had these mangled (KH_THETA held the c_h_inv values, C_H_INV
        # held garbage), corrupting every CH4/O2 gas-liquid partition in ch4Mod.
        @test CLM.C_H_INV  ≈ [1600.0, 1500.0, 2400.0]
        @test CLM.KH_THETA ≈ [714.29, 769.23, 29.4]
    end

end
