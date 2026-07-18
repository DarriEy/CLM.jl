@testset "Methane (ch4Mod)" begin

    # Helper to create minimal test data
    function make_ch4_test_data(; nc=3, np=4, ng=2, nlevsoi=5)
        params = CLM.CH4Params()
        ch4vc = CLM.CH4VarCon()

        ch4 = CLM.CH4Data(
            ch4_prod_depth_sat_col      = zeros(nc, nlevsoi),
            ch4_prod_depth_unsat_col    = zeros(nc, nlevsoi),
            ch4_prod_depth_lake_col     = zeros(nc, nlevsoi),
            ch4_oxid_depth_sat_col      = zeros(nc, nlevsoi),
            ch4_oxid_depth_unsat_col    = zeros(nc, nlevsoi),
            ch4_oxid_depth_lake_col     = zeros(nc, nlevsoi),
            ch4_aere_depth_sat_col      = zeros(nc, nlevsoi),
            ch4_aere_depth_unsat_col    = zeros(nc, nlevsoi),
            ch4_tran_depth_sat_col      = zeros(nc, nlevsoi),
            ch4_tran_depth_unsat_col    = zeros(nc, nlevsoi),
            ch4_ebul_depth_sat_col      = zeros(nc, nlevsoi),
            ch4_ebul_depth_unsat_col    = zeros(nc, nlevsoi),
            o2_oxid_depth_sat_col       = zeros(nc, nlevsoi),
            o2_oxid_depth_unsat_col     = zeros(nc, nlevsoi),
            o2_aere_depth_sat_col       = zeros(nc, nlevsoi),
            o2_aere_depth_unsat_col     = zeros(nc, nlevsoi),
            co2_decomp_depth_sat_col    = zeros(nc, nlevsoi),
            co2_decomp_depth_unsat_col  = zeros(nc, nlevsoi),
            co2_oxid_depth_sat_col      = zeros(nc, nlevsoi),
            co2_oxid_depth_unsat_col    = zeros(nc, nlevsoi),
            co2_aere_depth_sat_col      = zeros(nc, nlevsoi),
            co2_aere_depth_unsat_col    = zeros(nc, nlevsoi),

            ch4_ebul_total_sat_col      = zeros(nc),
            ch4_ebul_total_unsat_col    = zeros(nc),
            ch4_surf_aere_sat_col       = zeros(nc),
            ch4_surf_aere_unsat_col     = zeros(nc),
            ch4_surf_ebul_sat_col       = zeros(nc),
            ch4_surf_ebul_unsat_col     = zeros(nc),
            ch4_surf_ebul_lake_col      = zeros(nc),
            ch4_surf_diff_sat_col       = zeros(nc),
            ch4_surf_diff_unsat_col     = zeros(nc),
            ch4_surf_diff_lake_col      = zeros(nc),
            ch4_dfsat_flux_col          = zeros(nc),
            ch4_surf_flux_tot_col       = zeros(nc),

            conc_ch4_sat_col            = fill(1.0e-4, nc, nlevsoi),
            conc_ch4_unsat_col          = fill(1.0e-5, nc, nlevsoi),
            conc_ch4_lake_col           = zeros(nc, nlevsoi),
            conc_o2_sat_col             = fill(0.01, nc, nlevsoi),
            conc_o2_unsat_col           = fill(0.05, nc, nlevsoi),
            conc_o2_lake_col            = zeros(nc, nlevsoi),
            o2_decomp_depth_sat_col     = zeros(nc, nlevsoi),
            o2_decomp_depth_unsat_col   = zeros(nc, nlevsoi),

            o2stress_sat_col            = ones(nc, nlevsoi),
            o2stress_unsat_col          = ones(nc, nlevsoi),
            ch4stress_sat_col           = ones(nc, nlevsoi),
            ch4stress_unsat_col         = ones(nc, nlevsoi),

            zwt_ch4_unsat_col           = zeros(nc),
            lake_soilc_col              = fill(100.0, nc, nlevsoi),
            totcolch4_col               = zeros(nc),
            totcolch4_bef_col           = zeros(nc),
            annsum_counter_col          = zeros(nc),
            tempavg_somhr_col           = zeros(nc),
            annavg_somhr_col            = fill(1.0e-6, nc),
            tempavg_finrw_col           = zeros(nc),
            annavg_finrw_col            = fill(0.1, nc),
            sif_col                     = ones(nc),
            qflx_surf_lag_col           = zeros(nc),
            finundated_col              = fill(0.1, nc),
            finundated_pre_snow_col     = fill(0.1, nc),
            finundated_lag_col          = fill(0.1, nc),
            layer_sat_lag_col           = fill(0.5, nc, nlevsoi),
            pH_col                      = fill(6.5, nc),

            c_atm_grc                   = fill(0.03, ng, 3),
            ch4co2f_grc                 = zeros(ng),
            ch4prodg_grc                = zeros(ng),
            totcolch4_grc               = zeros(ng),
            totcolch4_bef_grc           = zeros(ng),

            annavg_agnpp_patch          = fill(1.0e-5, np),
            annavg_bgnpp_patch          = fill(1.0e-5, np),
            tempavg_agnpp_patch         = fill(1.0e-6, np),
            tempavg_bgnpp_patch         = fill(1.0e-6, np),

            grnd_ch4_cond_patch         = fill(0.01, np),
            grnd_ch4_cond_col           = fill(0.01, nc),

            ch4_first_time_grc          = fill(true, ng),
        )

        mask_soil = trues(nc)
        mask_soilp = trues(np)
        mask_nolake = trues(nc)
        mask_lake = falses(nc)

        # Soil properties
        watsat = fill(0.45, nc, nlevsoi)
        h2osoi_vol = fill(0.3, nc, nlevsoi)
        t_soisno = fill(CLM.TFRZ + 15.0, nc, nlevsoi)
        smp_l = fill(-1000.0, nc, nlevsoi)
        dz = fill(0.1, nc, nlevsoi)
        z = zeros(nc, nlevsoi)
        zi = zeros(nc, nlevsoi)
        for j in 1:nlevsoi
            for c in 1:nc
                z[c, j] = (j - 0.5) * 0.1
                zi[c, j] = j * 0.1
            end
        end

        return (;
            params, ch4vc, ch4, mask_soil, mask_soilp, mask_nolake, mask_lake,
            watsat, h2osoi_vol, t_soisno, smp_l, dz, z, zi,
            nc, np, ng, nlevsoi
        )
    end

    @testset "CH4Params construction" begin
        p = CLM.CH4Params()
        @test p.q10ch4 == 1.5
        @test p.f_ch4 == 0.2
        @test p.vmax_ch4_oxid == 0.0125
        @test p.vgc_max == 0.15
        @test p.satpow == 2.0
        @test p.f_sat == 0.95
    end

    @testset "CH4VarCon construction" begin
        # These MUST be the CTSM `ch4varcon.F90` module defaults. Three of them used to
        # be the OPPOSITE of Fortran's, which silently made the Julia methane model a
        # different model (production above the water table on, the CN/moisture HR
        # limitation removed, aerenchyma oxidation prescribed rather than prognostic).
        # The Fortran reference independently confirms anoxicmicrosites=.false.:
        # CH4_PROD_UNSAT is identically zero in the CTSM dump. See docs/CH4_FIRE_PARITY.md.
        vc = CLM.CH4VarCon()
        @test vc.allowlakeprod == false        # ch4varcon.F90:27
        @test vc.ch4offline == true            # ch4varcon.F90:54
        @test vc.transpirationloss == true     # ch4varcon.F90:21
        @test vc.anoxicmicrosites == false     # ch4varcon.F90:63  (was wrongly true)
        @test vc.ch4rmcnlim == false           # ch4varcon.F90:58  (was wrongly true)
        @test vc.use_aereoxid_prog == true     # ch4varcon.F90:17  (was wrongly false)
        @test vc.replenishlakec == true        # ch4varcon.F90:35
        @test vc.usephfact == false            # ch4varcon.F90:33
        @test vc.ch4frzout == false            # ch4varcon.F90:66
        @test vc.usefrootc == false            # ch4varcon.F90:47
        # finundated is now COMPUTED (CalcFinundated), not frozen at its cold-start 0.1.
        @test vc.finundation_mtd == CLM.FINUNDATION_MTD_H2OSFC
    end

    @testset "CH4Data construction" begin
        d = make_ch4_test_data()
        @test size(d.ch4.conc_ch4_sat_col) == (d.nc, d.nlevsoi)
        @test size(d.ch4.c_atm_grc) == (d.ng, 3)
        @test length(d.ch4.finundated_col) == d.nc
        @test length(d.ch4.ch4_first_time_grc) == d.ng
    end

    @testset "get_jwt!" begin
        d = make_ch4_test_data()
        jwt = zeros(Int, d.nc)
        params = CLM.CH4Params(f_sat=0.95)

        # With h2osoi_vol < f_sat * watsat everywhere, jwt should be at top (all unsaturated)
        h2osoi_vol_low = fill(0.3, d.nc, d.nlevsoi)
        CLM.get_jwt!(jwt, d.mask_soil, d.watsat, h2osoi_vol_low, d.t_soisno, d.nlevsoi, params)
        # All columns should have jwt = nlevsoi (fully unsaturated)
        for c in 1:d.nc
            @test jwt[c] == d.nlevsoi
        end

        # Now make bottom layers saturated
        h2osoi_vol_wet = fill(0.3, d.nc, d.nlevsoi)
        for c in 1:d.nc
            h2osoi_vol_wet[c, 4] = 0.44  # near saturation
            h2osoi_vol_wet[c, 5] = 0.44
        end
        CLM.get_jwt!(jwt, d.mask_soil, d.watsat, h2osoi_vol_wet, d.t_soisno, d.nlevsoi, params)
        for c in 1:d.nc
            @test jwt[c] == 3  # layer above water table
        end
    end

    @testset "ch4_totcolch4!" begin
        d = make_ch4_test_data()
        totcolch4 = zeros(d.nc)
        CLM.ch4_totcolch4!(totcolch4, d.ch4, d.mask_nolake, d.mask_lake,
                           d.dz, d.nlevsoi, false)
        # Each column: sum over j of (fin * conc_sat + (1-fin) * conc_unsat) * dz * CATOMW
        for c in 1:d.nc
            expected = 0.0
            for j in 1:d.nlevsoi
                expected += (0.1 * 1.0e-4 + 0.9 * 1.0e-5) * 0.1 * CLM.CATOMW
            end
            @test isapprox(totcolch4[c], expected, atol=1e-12)
        end
    end

    @testset "ch4_init_column_balance_check!" begin
        d = make_ch4_test_data()
        CLM.ch4_init_column_balance_check!(d.ch4, d.mask_nolake, d.mask_lake,
                                            d.dz, d.nlevsoi, false)
        for c in 1:d.nc
            @test d.ch4.totcolch4_bef_col[c] > 0.0  # Should have some CH4
        end
    end

    @testset "ch4_annualupdate!" begin
        d = make_ch4_test_data()
        patch_column = [1, 1, 2, 3]
        is_fates = falses(d.nc)
        somhr = fill(1.0e-6, d.nc)
        agnpp = fill(1.0e-5, d.np)
        bgnpp = fill(1.0e-5, d.np)
        dt = 1800.0
        secsperyear = 365.0 * 86400.0

        # Set counter below threshold
        d.ch4.annsum_counter_col .= 0.0

        CLM.ch4_annualupdate!(d.ch4, d.mask_soil, d.mask_soilp, patch_column,
                              is_fates, somhr, agnpp, bgnpp, dt, secsperyear)

        for c in 1:d.nc
            @test d.ch4.annsum_counter_col[c] == dt
            @test d.ch4.tempavg_somhr_col[c] > 0.0
        end
    end

    @testset "ch4_oxid! basic" begin
        d = make_ch4_test_data()
        jwt = fill(2, d.nc)  # water table at layer 2

        CLM.ch4_oxid!(d.ch4, d.params, d.mask_soil, d.watsat, d.h2osoi_vol,
                      d.smp_l, d.t_soisno, jwt, 0, false, d.nlevsoi, 1800.0)

        # Above WT (layers 1-2): should have oxidation since temp > TFRZ and CH4 present
        for c in 1:d.nc
            for j in 1:2
                @test d.ch4.ch4_oxid_depth_unsat_col[c, j] >= 0.0
            end
            # Below WT (layers 3-5): also should have oxidation
            for j in 3:d.nlevsoi
                @test d.ch4.ch4_oxid_depth_unsat_col[c, j] >= 0.0
            end
            # O2 oxidation should be 2x CH4 oxidation
            for j in 1:d.nlevsoi
                @test isapprox(d.ch4.o2_oxid_depth_unsat_col[c, j],
                               d.ch4.ch4_oxid_depth_unsat_col[c, j] * 2.0, atol=1e-15)
            end
        end
    end

    @testset "ch4_oxid! frozen soil" begin
        d = make_ch4_test_data()
        jwt = fill(2, d.nc)
        t_frozen = fill(CLM.TFRZ - 5.0, d.nc, d.nlevsoi)

        CLM.ch4_oxid!(d.ch4, d.params, d.mask_soil, d.watsat, d.h2osoi_vol,
                      d.smp_l, t_frozen, jwt, 0, false, d.nlevsoi, 1800.0)

        # All oxidation should be zero when frozen
        for c in 1:d.nc
            for j in 1:d.nlevsoi
                @test d.ch4.ch4_oxid_depth_unsat_col[c, j] == 0.0
            end
        end
    end

    @testset "ch4_ebul! basic" begin
        d = make_ch4_test_data()
        jwt = fill(0, d.nc)  # full saturation
        forc_pbot = fill(101325.0, d.nc)
        h2osfc = zeros(d.nc)
        frac_h2osfc = zeros(d.nc)
        lake_icefrac = zeros(d.nc, d.nlevsoi)
        lakedepth = fill(5.0, d.nc)

        # Set high CH4 concentration to trigger ebullition
        d.ch4.conc_ch4_sat_col .= 100.0  # Very high concentration

        CLM.ch4_ebul!(d.ch4, d.params, d.mask_soil, d.watsat, d.h2osoi_vol,
                      d.t_soisno, forc_pbot, h2osfc, frac_h2osfc, lake_icefrac,
                      lakedepth, d.z, d.dz, d.zi, jwt, 1, false, d.nlevsoi, 1800.0)

        # With very high concentrations, ebullition should occur
        has_ebul = false
        for c in 1:d.nc
            for j in 1:d.nlevsoi
                if d.ch4.ch4_ebul_depth_sat_col[c, j] > 0.0
                    has_ebul = true
                end
            end
        end
        @test has_ebul
    end

    @testset "site_ox_aere! basic" begin
        nlevsoi = 5
        tranloss = zeros(nlevsoi)
        aere = zeros(nlevsoi)
        oxaere = zeros(nlevsoi)

        watsat = fill(0.45, nlevsoi)
        h2osoi_vol = fill(0.3, nlevsoi)
        t_soisno = fill(CLM.TFRZ + 15.0, nlevsoi)
        conc_ch4 = fill(1.0e-4, nlevsoi)
        rootr = fill(0.2, nlevsoi)
        rootfr = fill(0.2, nlevsoi)
        conc_o2 = fill(0.05, nlevsoi)
        c_atm = [0.03, 8.0]
        z = [(j - 0.5) * 0.1 for j in 1:nlevsoi]
        dz = fill(0.1, nlevsoi)
        jwt = 0  # fully saturated
        params = CLM.CH4Params()
        ch4vc = CLM.CH4VarCon()

        CLM.site_ox_aere!(tranloss, aere, oxaere,
                          true, watsat, h2osoi_vol, t_soisno,
                          conc_ch4, rootr, 1.0e-3, jwt,
                          500.0, 1.0e-5, 1.0e-5,
                          2.0, 0.3, rootfr, 0.01,
                          conc_o2, c_atm, z, dz,
                          1, nlevsoi, params, ch4vc)

        # All layers below jwt=0 should have some aerenchyma flux
        for j in 1:nlevsoi
            @test aere[j] >= 0.0
        end
        # Transpiration loss should be non-negative
        for j in 1:nlevsoi
            @test tranloss[j] >= 0.0
        end
    end

    @testset "ch4_prod! basic" begin
        d = make_ch4_test_data()
        jwt = fill(2, d.nc)
        patch_column = [1, 1, 2, 3]
        patch_itype = [1, 2, 1, 1]
        patch_wtcol = [0.5, 0.5, 1.0, 1.0]
        is_fates = falses(d.nc)
        crootfr = fill(0.2, d.np, d.nlevsoi)
        rootfr_col = fill(0.2, d.nc, d.nlevsoi)
        somhr = fill(1.0e-6, d.nc)
        lithr = fill(1.0e-7, d.nc)
        hr_vr = fill(1.0e-8, d.nc, d.nlevsoi)
        o_scalar = ones(d.nc, d.nlevsoi)
        fphr = ones(d.nc, d.nlevsoi)
        pot_f_nit_vr = zeros(d.nc, d.nlevsoi)
        rr = fill(1.0e-6, d.np)

        CLM.ch4_prod!(d.ch4, d.params, d.ch4vc, d.mask_soil, d.mask_soilp,
                      patch_column, patch_itype, patch_wtcol, is_fates,
                      crootfr, rootfr_col, d.watsat, d.h2osoi_vol, d.t_soisno,
                      somhr, lithr, hr_vr, o_scalar, fphr, pot_f_nit_vr,
                      rr, jwt, 0, false,
                      d.dz, d.z, d.zi, d.nlevsoi, 1, d.nlevsoi, 5, 0.01,
                      1800.0, 0, true, false, false)

        # Below WT (j > 2): should have production
        for c in 1:d.nc
            for j in 3:d.nlevsoi
                @test d.ch4.ch4_prod_depth_unsat_col[c, j] >= 0.0
            end
            # Above WT with anoxic microsites: may have small production
            for j in 1:2
                @test d.ch4.ch4_prod_depth_unsat_col[c, j] >= 0.0
            end
        end
    end

    @testset "ch4_init_gridcell_balance_check!" begin
        d = make_ch4_test_data()
        col_gridcell = [1, 1, 2]
        col_wtgcell = [0.5, 0.5, 1.0]

        CLM.ch4_init_gridcell_balance_check!(d.ch4, d.mask_nolake, d.mask_lake,
                                              col_gridcell, col_wtgcell,
                                              d.dz, d.nlevsoi, d.ng, false)

        # Gridcell 1 should have positive CH4 (weighted avg of columns 1 and 2)
        @test d.ch4.totcolch4_bef_grc[1] > 0.0
        # Gridcell 2 should have positive CH4 (column 3 only)
        @test d.ch4.totcolch4_bef_grc[2] > 0.0
    end

    @testset "ch4! main driver" begin
        d = make_ch4_test_data()
        patch_column = [1, 1, 2, 3]
        patch_itype = [1, 2, 1, 1]
        patch_wtcol = [0.5, 0.5, 1.0, 1.0]
        col_gridcell = [1, 1, 2]
        col_wtgcell = [0.5, 0.5, 1.0]
        is_fates = falses(d.nc)
        latdeg = [45.0, 30.0]

        forc_pbot = fill(101325.0, d.ng)
        forc_t_grc = fill(CLM.TFRZ + 15.0, d.ng)
        forc_po2 = fill(101325.0 * 0.209, d.ng)
        forc_pco2 = fill(40.0, d.ng)
        forc_pch4 = fill(101325.0 * 1.7e-6, d.ng)

        h2osoi_liq = fill(0.3 * CLM.DENH2O * 0.1, d.nc, d.nlevsoi)
        h2osoi_ice = zeros(d.nc, d.nlevsoi)
        h2osfc_v = zeros(d.nc)
        bsw = fill(5.0, d.nc, d.nlevsoi)
        cellorg = fill(10.0, d.nc, d.nlevsoi)
        t_grnd = fill(CLM.TFRZ + 15.0, d.nc)
        t_h2osfc = fill(CLM.TFRZ + 15.0, d.nc)
        frac_h2osfc = zeros(d.nc)
        snow_depth_v = zeros(d.nc)
        snl_v = zeros(Int, d.nc)
        qflx_surf = zeros(d.nc)

        rootfr_p = fill(0.2, d.np, d.nlevsoi)
        rootfr_col = fill(0.2, d.nc, d.nlevsoi)
        crootfr = fill(0.2, d.np, d.nlevsoi)
        rootr_p = fill(0.2, d.np, d.nlevsoi)
        elai_v = fill(2.0, d.np)
        qflx_tran_veg = fill(1.0e-3, d.np)
        annsum_npp_v = fill(500.0, d.np)
        rr_v = fill(1.0e-6, d.np)

        somhr = fill(1.0e-6, d.nc)
        lithr = fill(1.0e-7, d.nc)
        hr_vr = fill(1.0e-8, d.nc, d.nlevsoi)
        o_scalar = ones(d.nc, d.nlevsoi)
        fphr = ones(d.nc, d.nlevsoi)
        pot_f_nit_vr = zeros(d.nc, d.nlevsoi)
        lake_icefrac = zeros(d.nc, d.nlevsoi)
        lakedepth = fill(5.0, d.nc)

        agnpp = fill(1.0e-5, d.np)
        bgnpp = fill(1.0e-5, d.np)

        # Initialize balance check
        d.ch4.ch4_first_time_grc .= true
        d.ch4.totcolch4_bef_col .= 0.0

        CLM.ch4!(d.ch4, d.params, d.ch4vc,
                 d.mask_soil, d.mask_soilp, d.mask_lake, d.mask_nolake,
                 col_gridcell, col_wtgcell, patch_column, patch_itype, patch_wtcol,
                 is_fates, latdeg,
                 forc_pbot, forc_t_grc, forc_po2, forc_pco2, forc_pch4,
                 d.watsat, d.h2osoi_vol, h2osoi_liq, h2osoi_ice, h2osfc_v,
                 bsw, cellorg, d.smp_l, d.t_soisno, t_grnd, t_h2osfc,
                 frac_h2osfc, snow_depth_v, snl_v, qflx_surf,
                 rootfr_p, rootfr_col, crootfr, rootr_p, elai_v,
                 qflx_tran_veg, annsum_npp_v, rr_v,
                 somhr, lithr, hr_vr, o_scalar, fphr, pot_f_nit_vr,
                 lake_icefrac, lakedepth,
                 d.z, d.dz, d.zi,
                 d.nlevsoi, 5, 1, d.nlevsoi,
                 5, 0.01, 130.0,
                 d.ng, 0, 1800.0,
                 true, false, false,
                 agnpp, bgnpp, 365.0 * 86400.0)

        # After first step, first_time should be false
        for g in 1:d.ng
            @test d.ch4.ch4_first_time_grc[g] == false
        end

        # Atmospheric concentrations should be set
        @test d.ch4.c_atm_grc[1, 1] > 0.0  # CH4
        @test d.ch4.c_atm_grc[1, 2] > 0.0  # O2

        # Surface flux total should be finite
        for c in 1:d.nc
            @test isfinite(d.ch4.ch4_surf_flux_tot_col[c])
        end

        # Total column CH4 should be computed
        for c in 1:d.nc
            @test isfinite(d.ch4.totcolch4_col[c])
        end

        # Gridcell-level totals should be non-negative
        for g in 1:d.ng
            @test d.ch4.totcolch4_grc[g] >= 0.0
        end
    end

    # ------------------------------------------------------------------
    # Satellite-inversion finundation (ch4FInundatedStreamType.F90
    # CalcFinundated, TWS_inversion + ZWT_inversion). Runs the full ch4!
    # driver so the whole dispatch (kernel gating + snow hold + redox lag)
    # is exercised, and asserts finundated_col vs the hand-computed
    # regression -- the DECOUPLED wetland-inundation path that produces
    # finundated > 0 where the h2osfc method drains to 0.
    # ------------------------------------------------------------------
    # Build the full ch4! call for a given finundation config and return
    # finundated_col. col_gridcell = [1,1,2], nc=3, ng=2.
    function run_ch4_finund(d; method, stream, tws, zwt, zwt_perched,
                            snow_depth=zeros(d.nc), qflx_surf=zeros(d.nc))
        d.ch4vc.finundation_mtd = method
        patch_column = [1, 1, 2, 3]; patch_itype = [1, 2, 1, 1]
        patch_wtcol = [0.5, 0.5, 1.0, 1.0]
        col_gridcell = [1, 1, 2]; col_wtgcell = [0.5, 0.5, 1.0]
        is_fates = falses(d.nc); latdeg = [45.0, 30.0]
        forc_pbot = fill(101325.0, d.ng); forc_t_grc = fill(CLM.TFRZ + 15.0, d.ng)
        forc_po2 = fill(101325.0 * 0.209, d.ng); forc_pco2 = fill(40.0, d.ng)
        forc_pch4 = fill(101325.0 * 1.7e-6, d.ng)
        h2osoi_liq = fill(0.3 * CLM.DENH2O * 0.1, d.nc, d.nlevsoi)
        h2osoi_ice = zeros(d.nc, d.nlevsoi); h2osfc_v = zeros(d.nc)
        bsw = fill(5.0, d.nc, d.nlevsoi); cellorg = fill(10.0, d.nc, d.nlevsoi)
        t_grnd = fill(CLM.TFRZ + 15.0, d.nc); t_h2osfc = fill(CLM.TFRZ + 15.0, d.nc)
        frac_h2osfc = zeros(d.nc); snl_v = zeros(Int, d.nc)
        rootfr_p = fill(0.2, d.np, d.nlevsoi); rootfr_col = fill(0.2, d.nc, d.nlevsoi)
        crootfr = fill(0.2, d.np, d.nlevsoi); rootr_p = fill(0.2, d.np, d.nlevsoi)
        elai_v = fill(2.0, d.np); qflx_tran_veg = fill(1.0e-3, d.np)
        annsum_npp_v = fill(500.0, d.np); rr_v = fill(1.0e-6, d.np)
        somhr = fill(1.0e-6, d.nc); lithr = fill(1.0e-7, d.nc)
        hr_vr = fill(1.0e-8, d.nc, d.nlevsoi); o_scalar = ones(d.nc, d.nlevsoi)
        fphr = ones(d.nc, d.nlevsoi); pot_f_nit_vr = zeros(d.nc, d.nlevsoi)
        lake_icefrac = zeros(d.nc, d.nlevsoi); lakedepth = fill(5.0, d.nc)
        agnpp = fill(1.0e-5, d.np); bgnpp = fill(1.0e-5, d.np)
        d.ch4.ch4_first_time_grc .= true; d.ch4.totcolch4_bef_col .= 0.0
        # qflx_surf_lag starts at 0; with redoxlag on, finundated_col is the
        # pre-lag clamped regression (snow_depth<=0 branch) -> deterministic.
        d.ch4.qflx_surf_lag_col .= 0.0
        d.ch4.finundated_lag_col .= 0.0
        CLM.ch4!(d.ch4, d.params, d.ch4vc,
                 d.mask_soil, d.mask_soilp, d.mask_lake, d.mask_nolake,
                 col_gridcell, col_wtgcell, patch_column, patch_itype, patch_wtcol,
                 is_fates, latdeg,
                 forc_pbot, forc_t_grc, forc_po2, forc_pco2, forc_pch4,
                 d.watsat, d.h2osoi_vol, h2osoi_liq, h2osoi_ice, h2osfc_v,
                 bsw, cellorg, d.smp_l, d.t_soisno, t_grnd, t_h2osfc,
                 frac_h2osfc, snow_depth, snl_v, qflx_surf,
                 rootfr_p, rootfr_col, crootfr, rootr_p, elai_v,
                 qflx_tran_veg, annsum_npp_v, rr_v,
                 somhr, lithr, hr_vr, o_scalar, fphr, pot_f_nit_vr,
                 lake_icefrac, lakedepth,
                 d.z, d.dz, d.zi,
                 d.nlevsoi, 5, 1, d.nlevsoi,
                 5, 0.01, 130.0,
                 d.ng, 0, 1800.0,
                 true, false, false,
                 agnpp, bgnpp, 365.0 * 86400.0;
                 finundated_stream=stream, tws=tws, zwt=zwt, zwt_perched=zwt_perched)
        return copy(d.ch4.finundated_col)
    end

    @testset "CalcFinundated TWS_inversion" begin
        d = make_ch4_test_data()
        # finundated = FWS_TWS_A * TWS + FWS_TWS_B (per gridcell), clamped [0,1].
        stream = CLM.CH4FInundatedStream(active=true,
                                         fws_slope_gdc=[0.01, 0.02],
                                         fws_intercept_gdc=[0.05, 0.10])
        tws = [20.0, 30.0]   # gridcell TWS (kg/m^2)
        fin = run_ch4_finund(d; method=CLM.FINUNDATION_MTD_TWS_INVERSION,
                             stream=stream, tws=tws, zwt=nothing, zwt_perched=nothing)
        # g1: 0.01*20+0.05 = 0.25 (cols 1,2);  g2: 0.02*30+0.10 = 0.70 (col 3)
        @test fin[1] ≈ 0.25 atol=1e-12
        @test fin[2] ≈ 0.25 atol=1e-12
        @test fin[3] ≈ 0.70 atol=1e-12
        @test all(fin .> 0.0)   # the whole point: finundated > 0 at a wetland col

        # Clamp: huge TWS -> 1, negative regression -> 0.
        d2 = make_ch4_test_data()
        stream2 = CLM.CH4FInundatedStream(active=true,
                                          fws_slope_gdc=[0.05, -0.01],
                                          fws_intercept_gdc=[0.0, 0.05])
        fin2 = run_ch4_finund(d2; method=CLM.FINUNDATION_MTD_TWS_INVERSION,
                              stream=stream2, tws=[100.0, 100.0],
                              zwt=nothing, zwt_perched=nothing)
        @test fin2[1] ≈ 1.0 atol=1e-12   # 0.05*100 = 5 -> clamp to 1
        @test fin2[3] ≈ 0.0 atol=1e-12   # -0.01*100+0.05 = -0.95 -> clamp to 0
    end

    @testset "CalcFinundated ZWT_inversion" begin
        d = make_ch4_test_data()
        # finundated = F0*exp(-zwt_actual/ZWT0) + P3*qflx_surf_lag, clamped.
        stream = CLM.CH4FInundatedStream(active=true,
                                         zwt0_gdc=[0.5, 0.5], f0_gdc=[0.9, 0.9],
                                         p3_gdc=[0.0, 0.0])
        # zwt_perched large -> condition false -> zwt_actual = zwt.
        zwt = [0.1, 0.1, 2.0]; zwt_perched = fill(10.0, d.nc)
        fin = run_ch4_finund(d; method=CLM.FINUNDATION_MTD_ZWT_INVERSION,
                             stream=stream, tws=nothing, zwt=zwt, zwt_perched=zwt_perched)
        @test fin[1] ≈ 0.9 * exp(-0.1 / 0.5) atol=1e-12   # ~0.7368
        @test fin[3] ≈ 0.9 * exp(-2.0 / 0.5) atol=1e-12   # ~0.0165
        @test all(fin .> 0.0)

        # Perched branch: zwt_perched < z[nlevsoi]-1e-5 (=0.44995) and < zwt.
        d2 = make_ch4_test_data()
        zwt2 = [0.4, 0.4, 0.4]; zwtp2 = [0.05, 0.05, 0.05]
        fin2 = run_ch4_finund(d2; method=CLM.FINUNDATION_MTD_ZWT_INVERSION,
                              stream=stream, tws=nothing, zwt=zwt2, zwt_perched=zwtp2)
        @test fin2[1] ≈ 0.9 * exp(-0.05 / 0.5) atol=1e-12   # uses perched, not zwt
    end

    @testset "read_ch4_finundated_stream! (real inversiondata file)" begin
        # The Prigent/Swenson regression stream, if present. Absence is a valid
        # partial outcome (the PORT is the deliverable); the test then skips.
        candidates = [
            joinpath(@__DIR__, "..", "..", "..", "..", "..", "..",
                     "compHydro", "SYMFLUENCE_data", "installs", "cesm-inputdata",
                     "lnd", "clm2", "paramdata",
                     "finundated_inversiondata_0.9x1.25_c170706.nc"),
            "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/installs/cesm-inputdata/lnd/clm2/paramdata/finundated_inversiondata_0.9x1.25_c170706.nc",
        ]
        path = ""
        for cand in candidates
            if isfile(cand); path = cand; break; end
        end
        if isempty(path)
            @info "finundated_inversiondata stream not found; skipping real-file read test"
            @test_skip false
        else
            # MerBleue peatland ~45.41 N, 75.52 W (=284.48 E).
            stream = CLM.CH4FInundatedStream()
            CLM.read_ch4_finundated_stream!(stream, path, [45.41], [284.48];
                                            method=CLM.FINUNDATION_MTD_TWS_INVERSION)
            @test stream.active
            @test length(stream.fws_slope_gdc) == 1
            @test isfinite(stream.fws_slope_gdc[1])
            @test isfinite(stream.fws_intercept_gdc[1])
            # Drive ch4! at MerBleue coeffs with a plausibly-wet TWS and confirm
            # the model can now produce finundated > 0 at a draining single point.
            # Use the default nc=3/ng=2 topology; put the MerBleue coeffs on both
            # gridcells so run_ch4_finund's hardcoded col_gridcell=[1,1,2] applies.
            a = stream.fws_slope_gdc[1]; b = stream.fws_intercept_gdc[1]
            d = make_ch4_test_data()
            stream1 = CLM.CH4FInundatedStream(active=true,
                        fws_slope_gdc=[a, a], fws_intercept_gdc=[b, b])
            # pick a TWS that yields an in-range regression for this cell
            tws_target = a != 0.0 ? clamp((0.5 - b) / a, 0.0, 5000.0) : 0.0
            fin = run_ch4_finund(d; method=CLM.FINUNDATION_MTD_TWS_INVERSION,
                                 stream=stream1, tws=[tws_target, tws_target],
                                 zwt=nothing, zwt_perched=nothing)
            @test fin[1] ≈ clamp(a * tws_target + b, 0.0, 1.0) atol=1e-10
        end
    end

end
