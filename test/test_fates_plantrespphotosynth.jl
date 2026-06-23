# test_fates_plantrespphotosynth.jl
# Tests for FATES Batch 16 (Tier F): FatesPlantRespPhotosynthMod — the FATES
# leaf-scale Farquhar/Collatz photosynthesis (C3/C4), stomatal conductance
# (Ball-Berry / Medlyn), leaf maintenance respiration (Ryan 1991 / Atkin 2017),
# and the scaling of leaf-layer fluxes up to the cohort.
#
# Strategy: exercise the standalone leaf-scale procedures directly with
# hand-computed / invariant-based checks (these do not require standing up a
# full site/patch/cohort linked list). We populate the EDPftvarcon_inst and
# ed_params() globals with a known single-PFT parameter set and restore them.

using Test
using CLM

@testset "FATES PlantRespPhotosynth (Batch 16)" begin

    tfrz = CLM.t_water_freeze_k_1atm

    # ---------------------------------------------------------------------
    # Temperature-response functions vs hand values
    # ---------------------------------------------------------------------
    @testset "ft1_f / fth_f / fth25_f" begin
        rgas_mol = CLM.rgas_J_K_kmol * 1.0e-3   # J/K/mol

        # ft1_f == 1.0 exactly at 25 C (the reference temperature)
        @test isapprox(CLM.ft1_f(tfrz + 25.0, 46390.0), 1.0; atol = 1e-12)

        # Arrhenius increases with temperature (positive activation energy)
        @test CLM.ft1_f(tfrz + 35.0, 46390.0) > CLM.ft1_f(tfrz + 25.0, 46390.0)
        @test CLM.ft1_f(tfrz + 15.0, 46390.0) < CLM.ft1_f(tfrz + 25.0, 46390.0)

        # Hand value at 35 C
        ha = 46390.0
        tl = tfrz + 35.0
        expect_ft1 = exp(ha / (rgas_mol * (tfrz + 25.0)) * (1.0 - (tfrz + 25.0) / tl))
        @test isapprox(CLM.ft1_f(tl, ha), expect_ft1; rtol = 1e-12)

        # fth25_f hand value
        hd = 150650.0; se = 490.0
        expect_fth25 = 1.0 + exp((-hd + se * (tfrz + 25.0)) / (rgas_mol * (tfrz + 25.0)))
        @test isapprox(CLM.fth25_f(hd, se), expect_fth25; rtol = 1e-12)

        # fth_f with scaleFactor = fth25_f(hd,se) gives ~1.0 at 25 C (normalized)
        sf = CLM.fth25_f(hd, se)
        @test isapprox(CLM.fth_f(tfrz + 25.0, hd, se, sf), 1.0; atol = 1e-12)

        # fth_f is a high-temperature inhibition: decreases at very high T
        @test CLM.fth_f(tfrz + 55.0, hd, se, sf) < CLM.fth_f(tfrz + 25.0, hd, se, sf)
    end

    # ---------------------------------------------------------------------
    # Set up a single-PFT EDPftvarcon + ed_params parameter set
    # ---------------------------------------------------------------------
    old_edpft = CLM.EDPftvarcon_inst[]
    old_params = CLM.EDParams[]
    old_numpft = CLM.numpft[]
    try
        CLM.numpft[] = 1

        edpft = CLM.EDPftvarcon_type()
        edpft.c3psn               = [1.0]            # C3
        edpft.bb_slope            = [9.0]            # Ball-Berry slope
        edpft.medlyn_slope        = [4.1]
        edpft.stomatal_intercept  = [10000.0]       # umol/m2/s (C3 default)
        edpft.maintresp_leaf_ryan1991_baserate  = [2.525e-6]
        edpft.maintresp_leaf_atkin2017_baserate = [1.756]
        edpft.maintresp_leaf_vert_scaler_coeff1 = [1.5]
        edpft.maintresp_leaf_vert_scaler_coeff2 = [3.0]
        edpft.maintresp_reduction_curvature     = [0.01]
        edpft.maintresp_reduction_intercept     = [1.0]
        edpft.vcmaxha = [65330.0]; edpft.jmaxha = [43540.0]
        edpft.vcmaxhd = [149250.0]; edpft.jmaxhd = [152040.0]
        edpft.vcmaxse = [485.0];    edpft.jmaxse = [495.0]
        edpft.hydr_k_lwp = [0.0]    # plant hydraulics off
        CLM.EDPftvarcon_inst[] = edpft

        edp = CLM.ed_params_type()
        edp.stomatal_model        = CLM.ballberry_model
        edp.stomatal_assim_model  = CLM.net_assim_model
        edp.photo_tempsens_model  = CLM.photosynth_acclim_model_none
        edp.dayl_switch           = 0
        edp.maintresp_leaf_model  = CLM.lmrmodel_ryan_1991
        edp.maintresp_nonleaf_baserate = 2.525e-6
        edp.q10_mr                = 1.5
        edp.theta_cj_c3           = 0.999
        edp.theta_cj_c4           = 0.999
        edp.radiation_model       = CLM.norman_solver
        CLM.EDParams[] = edp

        # Make sure the hlm flags this module reads are in benign states.
        CLM.hlm_use_planthydro[] = 0
        CLM.hlm_parteh_mode[]    = CLM.prt_carbon_allom_hyp
        CLM.nleafage[]           = 1

        ft = 1

        # Canopy gas environment shared by photosynthesis tests
        can_press = 1.0e5                 # Pa
        oair      = 0.209 * can_press     # O2 partial pressure (Pa)
        veg_tempk = tfrz + 25.0
        air_tempk = veg_tempk
        veg_esat  = 3160.0                # ~saturation vp at 25 C (Pa)
        air_vpress = 0.7 * veg_esat
        rb        = 50.0                  # boundary layer resistance s/m

        mm_kco2, mm_ko2, co2_cpoint, cf, gb_mol, ceair =
            CLM.GetCanopyGasParameters(can_press, oair, veg_tempk, air_tempk,
                                       air_vpress, veg_esat, rb)

        # ---------------------------------------------------------------------
        # GetCanopyGasParameters sanity
        # ---------------------------------------------------------------------
        @testset "GetCanopyGasParameters" begin
            @test mm_kco2 > 0.0 && mm_ko2 > 0.0 && co2_cpoint > 0.0
            @test cf > 0.0
            @test isapprox(gb_mol, cf / rb; rtol = 1e-12)
            # ceair constrained to [0.05*esat, esat]
            @test 0.05 * veg_esat <= ceair <= veg_esat
            # hand cf
            @test isapprox(cf, can_press / (CLM.rgas_J_K_kmol * air_tempk) * CLM.umol_per_kmol; rtol = 1e-12)
        end

        # Biophysical rates at canopy top, day time
        can_co2 = 40.0   # Pa (~400 ppm at 1e5 Pa)
        btran   = 1.0
        sib     = max(cf / 2.0e8, edpft.stomatal_intercept[ft] * btran)

        vcmax25top = 50.0
        jmax25top  = 85.0
        kp25top    = 20.0     # realistic C4 PEP initial-slope value (~vcmax scale)

        # ---------------------------------------------------------------------
        # LeafLayerBiophysicalRates
        # ---------------------------------------------------------------------
        @testset "LeafLayerBiophysicalRates" begin
            # Night (no PAR) -> zero rates, then vcmax*btran is still 0
            vc0, jm0, kp0 = CLM.LeafLayerBiophysicalRates(0.0, ft, vcmax25top,
                jmax25top, kp25top, 1.0, veg_tempk, 1.0, veg_tempk, veg_tempk, btran)
            @test vc0 == 0.0 && jm0 == 0.0 && kp0 == 0.0

            # Day: at 25 C with nscaler=1 and btran=1, vcmax ~ vcmax25top
            vc, jm, kp = CLM.LeafLayerBiophysicalRates(100.0, ft, vcmax25top,
                jmax25top, kp25top, 1.0, veg_tempk, 1.0, veg_tempk, veg_tempk, btran)
            @test isapprox(vc, vcmax25top; rtol = 1e-6)
            @test isapprox(jm, jmax25top; rtol = 1e-6)

            # btran scales vcmax linearly
            vc_half, _, _ = CLM.LeafLayerBiophysicalRates(100.0, ft, vcmax25top,
                jmax25top, kp25top, 1.0, veg_tempk, 1.0, veg_tempk, veg_tempk, 0.5)
            @test isapprox(vc_half, 0.5 * vc; rtol = 1e-10)
        end

        # ---------------------------------------------------------------------
        # LeafLayerPhotosynthesis — C3 light response, co-limitation, gs floor
        # ---------------------------------------------------------------------
        @testset "C3 photosynthesis light response + co-limitation" begin
            lmr = 1.0  # umol/m2/s

            # Build vcmax/jmax at the test temperature
            vcmax, jmax, kp = CLM.LeafLayerBiophysicalRates(100.0, ft, vcmax25top,
                jmax25top, kp25top, 1.0, veg_tempk, 1.0, veg_tempk, veg_tempk, btran)

            # helper to run at a given absorbed PAR (W/m2 per unit ground over laisun)
            function run_psn(par)
                psn, rs, anet, c13 = CLM.LeafLayerPhotosynthesis(
                    1.0,        # f_sun_lsl (all sunlit)
                    par,        # parsun
                    0.0,        # parsha
                    1.0,        # laisun
                    0.0,        # laisha
                    1.0,        # canopy area
                    ft, vcmax, jmax, kp, veg_tempk, veg_esat, can_press, can_co2,
                    oair, btran, sib, cf, gb_mol, ceair, mm_kco2, mm_ko2,
                    co2_cpoint, lmr, CLM.fates_unset_r8, rb)
                return psn, rs, anet, c13
            end

            psn_lo, rs_lo, anet_lo, c13_lo = run_psn(20.0)
            psn_mid, _, _, _               = run_psn(100.0)
            psn_hi, rs_hi, anet_hi, c13_hi = run_psn(800.0)
            psn_sat, _, _, _               = run_psn(2000.0)

            # GPP increases with light
            @test psn_lo < psn_mid < psn_hi
            # GPP saturates (Rubisco/export co-limit): high->very-high gain shrinks
            @test (psn_sat - psn_hi) < (psn_hi - psn_mid)
            @test psn_hi > 0.0

            # An = Agross - Rd (net = gross minus dark respiration)
            @test isapprox(anet_hi, psn_hi - lmr; rtol = 1e-9)
            @test anet_hi > 0.0   # in the light, net assimilation is positive

            # Stomatal resistance positive and bounded below rsmax0
            @test 0.0 < rs_hi < 2.0e8
            # Implied stomatal conductance gs_mol >= cuticular minimum (intercept)
            gs_mol_hi = cf / rs_hi
            @test gs_mol_hi >= edpft.stomatal_intercept[ft] * (1 - 1e-9)

            # c13 discrimination in a physical band
            @test 4.4 <= c13_hi <= 27.0
        end

        # ---------------------------------------------------------------------
        # Night time: psn = 0, anet = -lmr, rstoma = cf/intercept
        # ---------------------------------------------------------------------
        @testset "night time path" begin
            lmr = 1.3
            psn, rs, anet, c13 = CLM.LeafLayerPhotosynthesis(
                1.0, 0.0, 0.0, 1.0, 0.0, 1.0, ft, 30.0, 60.0, 0.0,
                veg_tempk, veg_esat, can_press, can_co2, oair, btran, sib,
                cf, gb_mol, ceair, mm_kco2, mm_ko2, co2_cpoint, lmr,
                CLM.fates_unset_r8, rb)
            @test psn == 0.0
            @test isapprox(anet, -lmr; rtol = 1e-12)
            @test isapprox(rs, cf / sib; rtol = 1e-12)
            @test c13 == 0.0
        end

        # ---------------------------------------------------------------------
        # C4 photosynthesis path (light-limited aj from quantum efficiency)
        # ---------------------------------------------------------------------
        @testset "C4 photosynthesis path" begin
            edpft.c3psn[ft] = 0.0              # switch to C4
            edpft.stomatal_intercept[ft] = 40000.0  # C4 Ball-Berry intercept
            CLM.EDPftvarcon_inst[] = edpft

            # Realistic C4 PEP initial-slope (FATES default ~ 20000) so the
            # export limit doesn't dominate the assimilation.
            kp25top_c4 = 20000.0

            # C4 vcmax temperature response differs; recompute rates
            vcmax, jmax, kp = CLM.LeafLayerBiophysicalRates(100.0, ft, vcmax25top,
                jmax25top, kp25top_c4, 1.0, veg_tempk, 1.0, veg_tempk, veg_tempk, btran)
            @test vcmax > 0.0

            sib_c4 = max(cf / 2.0e8, edpft.stomatal_intercept[ft] * btran)
            lmr = 1.0

            # NOTE: the FATES C4 shaded RuBP term has no div-by-zero guard
            # (preserved verbatim), so a leaf layer must carry BOTH sun and
            # shade leaf area. Use a partially sunlit layer (f_sun = 0.6).
            function run_c4(par_total)
                fsun = 0.6
                CLM.LeafLayerPhotosynthesis(
                    fsun, 0.7 * par_total, 0.3 * par_total, fsun, 1.0 - fsun, 1.0,
                    ft, vcmax, jmax, kp, veg_tempk, veg_esat, can_press, can_co2,
                    oair, btran, sib_c4, cf, gb_mol, ceair, mm_kco2, mm_ko2,
                    co2_cpoint, lmr, CLM.fates_unset_r8, rb)
            end

            psn_lo, _, _, _            = run_c4(80.0)
            psn_hi, rs_hi, anet_hi, _  = run_c4(600.0)
            @test psn_hi > psn_lo > 0.0
            # An = Agross - Rd at high light (positive net)
            @test anet_hi > 0.0
            @test 0.0 < rs_hi < 2.0e8

            edpft.c3psn[ft] = 1.0              # restore C3
            edpft.stomatal_intercept[ft] = 10000.0
            CLM.EDPftvarcon_inst[] = edpft
        end

        # ---------------------------------------------------------------------
        # Leaf maintenance respiration scaling with leaf N (Ryan 1991)
        # ---------------------------------------------------------------------
        @testset "maintenance respiration scaling (Ryan 1991)" begin
            lmr_lo = CLM.LeafLayerMaintenanceRespiration_Ryan_1991(1.0, 1.0, ft, veg_tempk)
            lmr_hi = CLM.LeafLayerMaintenanceRespiration_Ryan_1991(2.0, 1.0, ft, veg_tempk)
            # lmr linear in leaf N content (lnc_top)
            @test isapprox(lmr_hi, 2.0 * lmr_lo; rtol = 1e-10)
            @test lmr_lo > 0.0

            # nscaler attenuates respiration in lower canopy
            lmr_shade = CLM.LeafLayerMaintenanceRespiration_Ryan_1991(2.0, 0.5, ft, veg_tempk)
            @test isapprox(lmr_shade, 0.5 * lmr_hi; rtol = 1e-10)

            # warmer -> larger lmr (within the inhibition-free band)
            lmr_warm = CLM.LeafLayerMaintenanceRespiration_Ryan_1991(2.0, 1.0, ft, tfrz + 30.0)
            @test lmr_warm > lmr_hi
        end

        # ---------------------------------------------------------------------
        # Atkin 2017 maintenance respiration: increases with leaf N
        # ---------------------------------------------------------------------
        @testset "maintenance respiration (Atkin 2017)" begin
            lmr_lo = CLM.LeafLayerMaintenanceRespiration_Atkin_etal_2017(
                1.0, 0.0, vcmax25top, ft, veg_tempk, veg_tempk)
            lmr_hi = CLM.LeafLayerMaintenanceRespiration_Atkin_etal_2017(
                3.0, 0.0, vcmax25top, ft, veg_tempk, veg_tempk)
            @test lmr_hi > lmr_lo > 0.0

            # rdark scaler: deeper in canopy (larger cumulative_lai) -> less resp
            lmr_deep = CLM.LeafLayerMaintenanceRespiration_Atkin_etal_2017(
                3.0, 2.0, vcmax25top, ft, veg_tempk, veg_tempk)
            @test lmr_deep < lmr_hi
        end

        # ---------------------------------------------------------------------
        # lowstorage_maintresp_reduction: in [0,1], == 1 when storage >= target
        # ---------------------------------------------------------------------
        @testset "lowstorage_maintresp_reduction" begin
            @test CLM.lowstorage_maintresp_reduction(1.0, ft) == 1.0
            @test CLM.lowstorage_maintresp_reduction(2.0, ft) == 1.0
            r_half = CLM.lowstorage_maintresp_reduction(0.5, ft)
            r_low  = CLM.lowstorage_maintresp_reduction(0.1, ft)
            @test 0.0 <= r_low <= r_half <= 1.0   # monotone reduction
        end

        # ---------------------------------------------------------------------
        # ScaleLeafLayerFluxToCohort: integrates leaf-layer rates -> cohort
        # ---------------------------------------------------------------------
        @testset "ScaleLeafLayerFluxToCohort" begin
            nv = 2
            psn_llz = [10.0, 6.0]        # umolC/m2leaf/s
            lmr_llz = [1.0, 0.8]
            rs_llz  = [200.0, 300.0]
            elai    = [1.5, 1.0]
            c13     = [20.0, 18.0]
            c_area  = 2.0
            nplant  = 0.5
            rb_l    = 50.0
            mrf     = 0.9

            g_sb, gpp, rdark, c13clm, eleaf =
                CLM.ScaleLeafLayerFluxToCohort(nv, psn_llz, lmr_llz, rs_llz,
                    elai, c13, c_area, nplant, rb_l, mrf)

            # cohort effective leaf area = sum(elai * c_area)
            @test isapprox(eleaf, sum(elai) * c_area; rtol = 1e-12)

            # GPP hand computation [kgC/plant/s]
            gpp_expect = (psn_llz[1] * elai[1] * c_area + psn_llz[2] * elai[2] * c_area) *
                         CLM.umolC_to_kgC / nplant
            @test isapprox(gpp, gpp_expect; rtol = 1e-12)
            @test gpp > 0.0

            # Dark respiration hand computation (with mr reduction factor)
            rdark_expect = (lmr_llz[1] * elai[1] * c_area + lmr_llz[2] * elai[2] * c_area) *
                           CLM.umolC_to_kgC * mrf / nplant
            @test isapprox(rdark, rdark_expect; rtol = 1e-12)

            # g_sb_laweight = sum(1/(rs+rb) * cohort_layer_eleaf_area), NOT normalized
            g_expect = 1.0 / (rs_llz[1] + rb_l) * (elai[1] * c_area) +
                       1.0 / (rs_llz[2] + rb_l) * (elai[2] * c_area)
            @test isapprox(g_sb, g_expect; rtol = 1e-12)
            @test g_sb > 0.0
        end

        # ---------------------------------------------------------------------
        # RootLayerNFixation: positive MR surcharge & N fixed > 0
        # ---------------------------------------------------------------------
        @testset "RootLayerNFixation" begin
            # need nfix_mresp_scfrac on prt_params (PFT 1)
            old_scfrac = copy(CLM.prt_params.nfix_mresp_scfrac)
            try
                CLM.prt_params.nfix_mresp_scfrac = [0.1]
                fnrt_mr_layer = 1.0e-9   # kgC/s
                dtime = 1800.0
                t_soil = tfrz + 25.0
                surcharge, nfix = CLM.RootLayerNFixation(t_soil, ft, dtime, fnrt_mr_layer)
                @test isapprox(surcharge, fnrt_mr_layer * 0.1; rtol = 1e-12)
                # warm soil -> positive carbon cost -> positive N fixed
                @test nfix > 0.0
            finally
                CLM.prt_params.nfix_mresp_scfrac = old_scfrac
            end
        end

    finally
        CLM.EDPftvarcon_inst[] = old_edpft
        CLM.EDParams[]         = old_params
        CLM.numpft[]           = old_numpft
    end
end
