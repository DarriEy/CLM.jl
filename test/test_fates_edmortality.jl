# test_fates_edmortality.jl
# Tests for FATES Batch 13 (Tier F): EDMortalityFunctionsMod — the per-cohort
# mortality-rate functions (carbon starvation, hydraulic failure, background,
# freezing/cold, size/age senescence) plus the `mortality_rates` aggregator and
# `Mortality_Derivative` (which feeds d(n)/dt for the cohort-number ODE).
#
# Strategy (mirrors test_fates_cohort / test_fates_edcohortdynamics setup):
#   * Build a minimal single-PFT (woody, carbon-only) prt_params / FATES PFT /
#     ed_params table with hand-chosen mortality parameters so each rate is
#     analytically predictable.
#   * Seed a carbon PRT object with allometric pools, build a cohort, and a
#     minimal bc_in (soil temperature) + ed_site.
#   * Exercise each rate independently with hand-computed checks, the aggregate
#     `mortality_rates`, `ExemptTreefallDist`, and the sign/scale of
#     `Mortality_Derivative` for canopy-woody / understory / non-woody paths.
#   * Logging is disabled (hlm_use_logging = ifalse) so the logging terms are
#     identically zero and Mortality_Derivative reduces to the analytic rate sum.

using Test
using CLM

# Configure a single woody carbon PFT with explicit mortality parameters.
function _setup_mort_pft!()
    npft = 1
    p = CLM.prt_params
    CLM.allocate_prt_params!(p, npft, CLM.num_organ_types, 1)  # 1 leaf age class

    # Allometry needed for bleaf (callom) + the carbon PRT seed.
    p.c2b          .= 2.0
    p.wood_density .= 0.6
    p.slatop       .= 0.012
    p.slamax       .= 0.020
    p.allom_agb_frac .= 0.6
    p.allom_dbh_maxheight .= 90.0
    p.allom_la_per_sa_int .= 0.8e3
    p.allom_la_per_sa_slp .= 0.0
    p.allom_sai_scaler    .= 0.1
    p.allom_l2fr          .= 1.0
    p.cushion             .= 1.0
    p.allom_blca_expnt_diff      .= 0.0
    p.allom_d2ca_coefficient_min .= 0.3
    p.allom_d2ca_coefficient_max .= 0.6
    p.allom_h2cd1 .= 0.5
    p.allom_h2cd2 .= 1.0

    p.allom_dmode .= 1
    p.woody       .= 1
    p.allom_cmode .= 1
    p.allom_smode .= 1
    p.allom_fmode .= 1
    p.allom_stmode .= 1
    p.allom_hmode .= 1
    p.allom_amode .= 1
    p.allom_lmode .= 1

    p.allom_d2h1 .= 0.64
    p.allom_d2h2 .= 0.37
    p.allom_d2h3 .= -0.034

    p.allom_agb1 .= 0.06896; p.allom_agb2 .= 0.572; p.allom_agb3 .= 1.94; p.allom_agb4 .= 0.931
    p.allom_d2bl1 .= 0.07; p.allom_d2bl2 .= 1.3; p.allom_d2bl3 .= 0.55

    # Deciduous flags: evergreen (not deciduous, not dormant) for the default path.
    p.stress_decid .= 0
    p.season_decid .= 0
    p.evergreen    .= CLM.itrue
    p.leaf_stor_priority .= 0.8
    p.leaf_long          .= 1.5
    p.leaf_long_ustory   .= 1.5

    # param_derived: branch_frac (used by bsap_allom) + canopy-top derived rates.
    pd = CLM.param_derived_type()
    pd.branch_frac = fill(0.25, npft)
    pd.jmax25top = fill(85.0, npft, 1)
    pd.tpu25top  = fill(8.0, npft, 1)
    pd.kp25top   = fill(0.6, npft, 1)
    CLM.ParamDerived[] = pd

    # FATES PFT table: mortality parameters chosen so each rate is predictable.
    ev = CLM.EDPftvarcon_type()
    ev.vcmax25top = fill(50.0, npft, 1)

    ev.bmort                   = fill(0.014, npft)   # background rate [/yr]
    ev.mort_scalar_coldstress  = fill(3.0, npft)     # max cold-stress rate
    ev.mort_scalar_cstarvation = fill(0.6, npft)     # max C-starvation rate
    ev.mort_scalar_hydrfailure = fill(0.5, npft)     # max hydraulic-failure rate
    ev.mort_upthresh_cstarvation = fill(1.0, npft)   # C-starvation upper threshold / e-fold
    ev.freezetol               = fill(-5.0, npft)    # freezing tolerance [degC]
    ev.hf_sm_threshold         = fill(1.0e-6, npft)  # btran hyd-failure threshold
    ev.hf_flc_threshold        = fill(0.5, npft)
    # Senescence OFF by default (ip >= fates_check_param_set => term = 0).
    ev.mort_ip_size_senescence = fill(CLM.fates_check_param_set, npft)
    ev.mort_r_size_senescence  = fill(0.0, npft)
    ev.mort_ip_age_senescence  = fill(CLM.fates_check_param_set, npft)
    ev.mort_r_age_senescence   = fill(0.0, npft)
    ev.prescribed_mortality_canopy     = fill(0.02, npft)
    ev.prescribed_mortality_understory = fill(0.04, npft)
    # Damage parameters (only used if hlm_use_tree_damage == itrue).
    ev.damage_mort_p1 = fill(0.4, npft)
    ev.damage_mort_p2 = fill(5.0, npft)
    CLM.EDPftvarcon_inst[] = ev

    # ed_params: cstarvation model + disturbance fraction + damage bins.
    edp = CLM.ed_params_type()
    edp.mort_cstarvation_model = CLM.cstarvation_model_lin
    edp.fates_mortality_disturbance_fraction = 1.0
    edp.ED_val_history_damage_bin_edges = [0.0, 20.0, 50.0, 80.0]
    edp.logging_event_code = 1.0  # logging off
    CLM.EDParams[] = edp

    return npft
end

# Carbon PRT object seeded with allometric pools; storage scaled by store_mult.
function _seed_prt_mort(ipft::Int, dbh0::Float64; store_mult::Float64 = 1.0)
    prt = CLM.callom_prt_vartypes()
    CLM.InitPRTVartype!(prt)

    canopy_trim = 1.0; crowndamage = 1; ef = 1.0
    l2fr = CLM.prt_params.allom_l2fr[ipft]

    tgt_leaf, _    = CLM.bleaf(dbh0, ipft, crowndamage, canopy_trim, ef)
    tgt_fnrt, _    = CLM.bfineroot(dbh0, ipft, canopy_trim, l2fr, ef)
    _, tgt_sapw, _ = CLM.bsap_allom(dbh0, ipft, crowndamage, canopy_trim, ef)
    tgt_store, _   = CLM.bstore_allom(dbh0, ipft, crowndamage, canopy_trim)
    tgt_agw, _     = CLM.bagw_allom(dbh0, ipft, crowndamage, ef)
    tgt_bgw, _     = CLM.bbgw_allom(dbh0, ipft, ef)
    tgt_struct, _  = CLM.bdead_allom(tgt_agw, tgt_bgw, tgt_sapw, ipft)

    CLM.SetState!(prt, CLM.leaf_organ,   CLM.carbon12_element, tgt_leaf)
    CLM.SetState!(prt, CLM.fnrt_organ,   CLM.carbon12_element, tgt_fnrt)
    CLM.SetState!(prt, CLM.sapw_organ,   CLM.carbon12_element, tgt_sapw)
    CLM.SetState!(prt, CLM.store_organ,  CLM.carbon12_element, tgt_store * store_mult)
    CLM.SetState!(prt, CLM.struct_organ, CLM.carbon12_element, tgt_struct)
    CLM.SetState!(prt, CLM.repro_organ,  CLM.carbon12_element, 0.0)

    return prt
end

# Build a bare cohort with the given structural/phenology state.
function _build_cohort(ipft::Int, dbh0::Float64, prt; n = 1.0, clayer = 1,
                       status = 2, crowndamage = 1)
    coh = CLM.fates_cohort_type()
    coh.prt          = prt
    coh.pft          = ipft
    coh.n            = n
    coh.dbh          = dbh0
    coh.coage        = 5.0
    coh.canopy_trim  = 1.0
    coh.crowndamage  = crowndamage
    coh.canopy_layer = clayer
    coh.status_coh   = status
    coh.efleaf_coh   = 1.0
    return coh
end

@testset "FATES Batch 13: EDMortalityFunctionsMod" begin

    old_global   = CLM.prt_global[]
    old_numel    = CLM.num_elements[]
    old_ellist   = copy(CLM.element_list)
    old_parteh   = CLM.hlm_parteh_mode[]
    old_use_sp   = CLM.hlm_use_sp[]
    old_hydro    = CLM.hlm_use_planthydro[]
    old_presc    = CLM.hlm_use_ed_prescribed_phys[]
    old_damage   = CLM.hlm_use_tree_damage[]
    old_logging  = CLM.hlm_use_logging[]
    old_luh      = CLM.hlm_use_luh[]
    old_nharv    = CLM.hlm_num_lu_harvest_cats[]
    old_freqday  = CLM.hlm_freq_day[]
    old_nleafage = CLM.nleafage[]
    old_pd       = CLM.ParamDerived[]
    old_edp      = CLM.EDParams[]
    old_ev       = CLM.EDPftvarcon_inst[]

    try
        _setup_mort_pft!()
        ipft = 1
        CLM.InitPRTGlobalAllometricCarbon!()

        # Carbon-only element registry (num_elements = 1, [carbon12_element]).
        CLM.num_elements[] = 1
        empty!(CLM.element_list)
        push!(CLM.element_list, CLM.carbon12_element)

        # Carbon-only flags.
        CLM.hlm_parteh_mode[]            = CLM.prt_carbon_allom_hyp
        CLM.hlm_use_sp[]                 = CLM.ifalse
        CLM.hlm_use_planthydro[]         = CLM.ifalse
        CLM.hlm_use_ed_prescribed_phys[] = CLM.ifalse
        CLM.hlm_use_tree_damage[]        = CLM.ifalse
        CLM.hlm_use_logging[]            = CLM.ifalse   # logging terms -> 0
        CLM.hlm_use_luh[]                = CLM.ifalse
        CLM.hlm_num_lu_harvest_cats[]    = 5
        CLM.hlm_freq_day[]               = 1.0 / 365.0
        CLM.nleafage[]                   = 1

        dbh0 = 15.0

        # bc_in: warm, unfrozen soil so the soil-tfrz gate is open.
        bc_in = CLM.bc_in_type()
        bc_in.t_soisno_sl = fill(290.0, 4)
        bc_in.hlm_harvest_rates    = zeros(Float64, CLM.hlm_num_lu_harvest_cats[])
        bc_in.hlm_harvest_catnames = fill("", CLM.hlm_num_lu_harvest_cats[])
        bc_in.hlm_harvest_units    = 1

        # btran: above hf_sm_threshold (no hyd-failure) unless we lower it.
        btran_hi = fill(0.9, CLM.maxpft)
        btran_lo = fill(0.0, CLM.maxpft)

        # mean_temp warm (well above freezetol => no cold stress).
        T_warm = 290.0   # ~16.85 degC

        # -------------------------------------------------------------------
        # 1. Carbon starvation: increases as storage -> 0 (linear model).
        # -------------------------------------------------------------------
        # Full storage (store=target) => frac large => cmort = 0.
        coh_full = _build_cohort(ipft, dbh0, _seed_prt_mort(ipft, dbh0; store_mult = 1.0))
        cm_full = CLM.mortality_rates(coh_full, bc_in, btran_hi, T_warm)[1]
        @test cm_full == 0.0

        # Half storage. With upthresh=1.0 and frac = store/target, frac<1 so
        # cmort = scalar * (upthresh - frac). Decreasing store increases cmort.
        coh_half = _build_cohort(ipft, dbh0, _seed_prt_mort(ipft, dbh0; store_mult = 0.5))
        coh_low  = _build_cohort(ipft, dbh0, _seed_prt_mort(ipft, dbh0; store_mult = 0.1))
        coh_zero = _build_cohort(ipft, dbh0, _seed_prt_mort(ipft, dbh0; store_mult = 0.0))

        cm_half = CLM.mortality_rates(coh_half, bc_in, btran_hi, T_warm)[1]
        cm_low  = CLM.mortality_rates(coh_low,  bc_in, btran_hi, T_warm)[1]
        cm_zero = CLM.mortality_rates(coh_zero, bc_in, btran_hi, T_warm)[1]

        @test cm_half > 0.0
        @test cm_low  > cm_half           # monotone increasing as storage drops
        @test cm_zero > cm_low
        # At store=0, frac=0, linear model => cmort = scalar * upthresh = 0.6.
        @test cm_zero ≈ 0.6 atol = 1e-12

        # Hand check the half case: frac = store_c / target_leaf_c.
        tgt_leaf, _ = CLM.bleaf(dbh0, ipft, 1, 1.0, 1.0)
        store_half  = CLM.GetState(coh_half.prt, CLM.store_organ, CLM.carbon12_element)
        frac_half   = max(0.0, store_half / max(tgt_leaf, CLM.nearzero))
        cm_expected = frac_half >= 1.0 ? 0.0 : 0.6 * (1.0 - frac_half)
        @test cm_half ≈ cm_expected atol = 1e-12

        # -------------------------------------------------------------------
        #  Carbon starvation: exponential model.
        # -------------------------------------------------------------------
        CLM.EDParams[].mort_cstarvation_model = CLM.cstarvation_model_exp
        cm_exp_zero = CLM.mortality_rates(coh_zero, bc_in, btran_hi, T_warm)[1]
        cm_exp_low  = CLM.mortality_rates(coh_low,  bc_in, btran_hi, T_warm)[1]
        # exp model: cmort = scalar * exp(-frac/upthresh); at frac=0 => scalar.
        @test cm_exp_zero ≈ 0.6 atol = 1e-12
        @test cm_exp_zero > cm_exp_low     # still monotone
        CLM.EDParams[].mort_cstarvation_model = CLM.cstarvation_model_lin

        # -------------------------------------------------------------------
        # 2. Cold / freezing mortality vs temperature.
        # -------------------------------------------------------------------
        # freezetol = -5 degC, buffer = 5 degC, scalar = 3.0.
        # temp_dep = clamp(1 - (T_C - freezetol)/buffer, 0, 1).
        coh = _build_cohort(ipft, dbh0, _seed_prt_mort(ipft, dbh0; store_mult = 1.0))

        # Warm: T_C >> freezetol => temp_dep = 0 => frmort = 0.
        fr_warm = CLM.mortality_rates(coh, bc_in, btran_hi, T_warm)[4]
        @test fr_warm == 0.0

        # Exactly at freezetol (-5 degC = 268.15 K): temp_dep = 1 => frmort = 3.0.
        T_at = CLM.t_water_freeze_k_1atm + (-5.0)
        fr_at = CLM.mortality_rates(coh, bc_in, btran_hi, T_at)[4]
        @test fr_at ≈ 3.0 atol = 1e-12

        # Mid-buffer (-2.5 degC): temp_dep = 1 - (2.5)/5 = 0.5 => frmort = 1.5.
        T_mid = CLM.t_water_freeze_k_1atm + (-2.5)
        fr_mid = CLM.mortality_rates(coh, bc_in, btran_hi, T_mid)[4]
        @test fr_mid ≈ 1.5 atol = 1e-12

        # Very cold (below freezetol - buffer): clamped to 1 => frmort = 3.0 (max).
        T_cold = CLM.t_water_freeze_k_1atm + (-20.0)
        fr_cold = CLM.mortality_rates(coh, bc_in, btran_hi, T_cold)[4]
        @test fr_cold ≈ 3.0 atol = 1e-12

        # Monotone decreasing rate as temperature rises across the buffer.
        @test fr_cold >= fr_at >= fr_mid >= fr_warm

        # -------------------------------------------------------------------
        # 3. Background mortality = constant PFT rate.
        # -------------------------------------------------------------------
        bm = CLM.mortality_rates(coh, bc_in, btran_hi, T_warm)[3]
        @test bm ≈ 0.014 atol = 1e-12

        # -------------------------------------------------------------------
        # 4. Hydraulic-failure mortality (btran path, hydro off).
        # -------------------------------------------------------------------
        # btran above threshold => no hyd failure.
        hm_hi = CLM.mortality_rates(coh, bc_in, btran_hi, T_warm)[2]
        @test hm_hi == 0.0
        # btran below threshold (0 <= 1e-6), leaves on, soil unfrozen => full rate.
        hm_lo = CLM.mortality_rates(coh, bc_in, btran_lo, T_warm)[2]
        @test hm_lo ≈ 0.5 atol = 1e-12

        # Deciduous + dormant (leaves off) suppresses hyd failure even at low btran.
        CLM.prt_params.season_decid[ipft] = CLM.itrue
        coh_dorm = _build_cohort(ipft, dbh0, _seed_prt_mort(ipft, dbh0);
                                 status = CLM.leaves_off)
        hm_dorm = CLM.mortality_rates(coh_dorm, bc_in, btran_lo, T_warm)[2]
        @test hm_dorm == 0.0
        CLM.prt_params.season_decid[ipft] = 0  # restore evergreen

        # Frozen soil also suppresses hyd failure (soil_tfrz_thresh = -2 degC).
        bc_frozen = CLM.bc_in_type()
        bc_frozen.t_soisno_sl = fill(269.0, 4)  # ~ -4 degC < -2 threshold
        bc_frozen.hlm_harvest_rates    = zeros(Float64, CLM.hlm_num_lu_harvest_cats[])
        bc_frozen.hlm_harvest_catnames = fill("", CLM.hlm_num_lu_harvest_cats[])
        bc_frozen.hlm_harvest_units    = 1
        hm_frozen = CLM.mortality_rates(coh, bc_frozen, btran_lo, T_warm)[2]
        @test hm_frozen == 0.0

        # -------------------------------------------------------------------
        # 5. Size- and age-dependent senescence (turn on).
        # -------------------------------------------------------------------
        ev = CLM.EDPftvarcon_inst[]
        # Size senescence ON: ip = 10 cm, r = 0.1. dbh0 = 15 => logistic > 0.5.
        ev.mort_ip_size_senescence[ipft] = 10.0
        ev.mort_r_size_senescence[ipft]  = 0.1
        sm = CLM.mortality_rates(coh, bc_in, btran_hi, T_warm)[5]
        sm_expected = 1.0 / (1.0 + exp(-0.1 * (dbh0 - 10.0)))
        @test sm ≈ sm_expected atol = 1e-12
        @test 0.5 < sm < 1.0
        # Larger tree => larger smort.
        coh_big = _build_cohort(ipft, 40.0, _seed_prt_mort(ipft, 40.0))
        sm_big = CLM.mortality_rates(coh_big, bc_in, btran_hi, T_warm)[5]
        @test sm_big > sm

        # Age senescence ON: ip = 3 yr, r = 0.2. coage = 5 => logistic > 0.5.
        ev.mort_ip_age_senescence[ipft] = 3.0
        ev.mort_r_age_senescence[ipft]  = 0.2
        asm = CLM.mortality_rates(coh, bc_in, btran_hi, T_warm)[6]
        asm_expected = 1.0 / (1.0 + exp(-0.2 * (5.0 - 3.0)))
        @test asm ≈ asm_expected atol = 1e-12
        @test asm > 0.5

        # Restore senescence-off.
        ev.mort_ip_size_senescence[ipft] = CLM.fates_check_param_set
        ev.mort_r_size_senescence[ipft]  = 0.0
        ev.mort_ip_age_senescence[ipft]  = CLM.fates_check_param_set
        ev.mort_r_age_senescence[ipft]   = 0.0
        # With senescence off again, both are zero.
        sm_off, asm_off = CLM.mortality_rates(coh, bc_in, btran_hi, T_warm)[[5, 6]]
        @test sm_off == 0.0
        @test asm_off == 0.0

        # -------------------------------------------------------------------
        # 6. Damage mortality gating (hlm_use_tree_damage).
        # -------------------------------------------------------------------
        dg_off = CLM.mortality_rates(coh, bc_in, btran_hi, T_warm)[7]
        @test dg_off == 0.0   # damage flag off

        # -------------------------------------------------------------------
        # 7. Aggregate mortality_rates returns the 7-tuple consistently.
        # -------------------------------------------------------------------
        cm, hm, bmo, frm, smo, asmo, dgm =
            CLM.mortality_rates(coh, bc_in, btran_lo, T_cold)
        @test bmo ≈ 0.014 atol = 1e-12      # background
        @test hm  ≈ 0.5   atol = 1e-12      # hydraulic failure (low btran)
        @test frm ≈ 3.0   atol = 1e-12      # cold (very cold)
        @test cm  == 0.0                     # storage full
        @test smo == 0.0 && asmo == 0.0 && dgm == 0.0
        @test all(isfinite, (cm, hm, bmo, frm, smo, asmo, dgm))

        # -------------------------------------------------------------------
        # 8. ExemptTreefallDist: woody => false, non-woody => true.
        # -------------------------------------------------------------------
        @test CLM.ExemptTreefallDist(coh) == false   # woody PFT
        CLM.prt_params.woody[ipft] = CLM.ifalse
        @test CLM.ExemptTreefallDist(coh) == true     # non-woody
        CLM.prt_params.woody[ipft] = CLM.itrue        # restore

        # -------------------------------------------------------------------
        # 9. Mortality_Derivative sign / scale.
        # -------------------------------------------------------------------
        site = CLM.ed_site_type()
        site.transition_landuse_from_off_to_on = false
        site.min_allowed_landuse_fraction = 0.0
        harvestable = zeros(Float64, CLM.hlm_num_lu_harvest_cats[])
        htag = fill(2, CLM.hlm_num_lu_harvest_cats[])

        # Canopy woody cohort, disturbance_fraction = 1.0 => non-exempt canopy
        # death is ALL disturbance, so dndt = (1 - 1.0)*(...) = 0.
        nn = 100.0
        coh_can = _build_cohort(ipft, dbh0, _seed_prt_mort(ipft, dbh0; store_mult = 0.5);
                                n = nn, clayer = 1)
        CLM.Mortality_Derivative(site, coh_can, bc_in, btran_lo, T_cold, 1, 0.0,
                                 1.0, -1.0, harvestable, htag)
        rate_sum = coh_can.cmort + coh_can.hmort + coh_can.bmort + coh_can.frmort +
                   coh_can.smort + coh_can.asmort + coh_can.dgmort
        @test rate_sum > 0.0
        @test coh_can.dndt ≈ 0.0 atol = 1e-12          # (1 - 1.0) * (...) = 0
        @test coh_can.dndt <= 0.0                        # never positive

        # Now disturbance_fraction = 0.0 => canopy non-exempt death fully here:
        # dndt = -(sum) * n.
        CLM.EDParams[].fates_mortality_disturbance_fraction = 0.0
        coh_can2 = _build_cohort(ipft, dbh0, _seed_prt_mort(ipft, dbh0; store_mult = 0.5);
                                 n = nn, clayer = 1)
        CLM.Mortality_Derivative(site, coh_can2, bc_in, btran_lo, T_cold, 1, 0.0,
                                 1.0, -1.0, harvestable, htag)
        rate_sum2 = coh_can2.cmort + coh_can2.hmort + coh_can2.bmort + coh_can2.frmort +
                    coh_can2.smort + coh_can2.asmort + coh_can2.dgmort
        @test coh_can2.dndt ≈ -rate_sum2 * nn atol = 1e-10
        @test coh_can2.dndt < 0.0

        # Understory woody cohort (canopy_layer = 2): all mortality here, no
        # disturbance-fraction reduction (logging off => dndt_logging = 0).
        coh_us = _build_cohort(ipft, dbh0, _seed_prt_mort(ipft, dbh0; store_mult = 0.5);
                               n = nn, clayer = 2)
        CLM.Mortality_Derivative(site, coh_us, bc_in, btran_lo, T_cold, 1, 0.0,
                                 1.0, -1.0, harvestable, htag)
        rate_sum_us = coh_us.cmort + coh_us.hmort + coh_us.bmort + coh_us.frmort +
                      coh_us.smort + coh_us.asmort + coh_us.dgmort
        @test coh_us.dndt ≈ -rate_sum_us * nn atol = 1e-10
        @test coh_us.dndt < 0.0
        # Logging fields zeroed (logging off).
        @test coh_us.lmort_direct == 0.0
        @test coh_us.lmort_collateral == 0.0
        @test coh_us.lmort_infra == 0.0

        # Non-woody (grass) canopy cohort is exempt from treefall disturbance, so
        # ALL its mortality stays here (no disturbance-fraction reduction), even
        # at disturbance_fraction != 0.
        CLM.EDParams[].fates_mortality_disturbance_fraction = 1.0
        CLM.prt_params.woody[ipft] = CLM.ifalse
        coh_grass = _build_cohort(ipft, dbh0, _seed_prt_mort(ipft, dbh0; store_mult = 0.5);
                                  n = nn, clayer = 1)
        CLM.Mortality_Derivative(site, coh_grass, bc_in, btran_lo, T_cold, 1, 0.0,
                                 1.0, -1.0, harvestable, htag)
        rate_sum_g = coh_grass.cmort + coh_grass.hmort + coh_grass.bmort + coh_grass.frmort +
                     coh_grass.smort + coh_grass.asmort + coh_grass.dgmort
        @test coh_grass.dndt ≈ -rate_sum_g * nn atol = 1e-10
        @test coh_grass.dndt < 0.0
        CLM.prt_params.woody[ipft] = CLM.itrue

        # harvest_tag passed through (logging off => all 2 = not applicable).
        @test all(htag .== 2)

    finally
        CLM.prt_global[]                 = old_global
        CLM.num_elements[]               = old_numel
        empty!(CLM.element_list); append!(CLM.element_list, old_ellist)
        CLM.hlm_parteh_mode[]            = old_parteh
        CLM.hlm_use_sp[]                 = old_use_sp
        CLM.hlm_use_planthydro[]         = old_hydro
        CLM.hlm_use_ed_prescribed_phys[] = old_presc
        CLM.hlm_use_tree_damage[]        = old_damage
        CLM.hlm_use_logging[]            = old_logging
        CLM.hlm_use_luh[]                = old_luh
        CLM.hlm_num_lu_harvest_cats[]    = old_nharv
        CLM.hlm_freq_day[]               = old_freqday
        CLM.nleafage[]                   = old_nleafage
        CLM.ParamDerived[]               = old_pd
        CLM.EDParams[]                   = old_edp
        CLM.EDPftvarcon_inst[]           = old_ev
    end
end
