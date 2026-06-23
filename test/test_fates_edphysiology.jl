# test_fates_edphysiology.jl
# Tests for FATES Batch 13 (Tier F): EDPhysiologyMod — the physiology hub.
#
# Strategy (mirrors test_fates_edcohortdynamics.jl's setup harness):
#   * Build a synthetic 3-PFT carbon-only prt_params table (PFT1 woody evergreen,
#     PFT2 woody cold-deciduous, PFT3 woody drought-hard-deciduous), install
#     param_derived, the FATES PFT table, ed_params (phenology + regeneration +
#     VAI bins + trim params) and sf_params (CWD/decomp fractions).
#   * Register the carbon-only PARTEH hypothesis; flags = carbon, no SP, no plant
#     hydro, nleafage = 1; numpft = 3; carbon-only element registry.
#   * Exercise:
#       - phenology cold-deciduous: GDD exceedance flips cstatus -> notcold
#         (leaf-on) given enough warmth + chilling history; then a cold snap flips
#         it back to iscold (leaf-off), and the cohort sheds leaves.
#       - phenology drought-deciduous: dry rooting-zone soil drives elong_factor
#         to 0 (leaf-off) and dstatus to moistoff after the minimum on-period.
#       - trim_canopy: a cohort whose bottom leaf layers are in negative annual
#         carbon balance has its canopy_trim reduced toward carbon balance.
#       - recruitment: a seed_germ pool produces a cohort of the right juvenile
#         size (height = hgt_min, dbh = h2d_allom(hgt_min)).
#       - litter / seed mass conservation: PreDisturbanceIntegrateLitter conserves
#         (state + frag-out) given a known input flux.

using Test
using CLM

# ---------------------------------------------------------------------------
# Parameter setup: 3 woody carbon PFTs differing only in phenology strategy.
#   PFT1: evergreen          (season_decid=0, stress_decid=0)
#   PFT2: cold-deciduous     (season_decid=1, stress_decid=0)
#   PFT3: drought-hard-decid (season_decid=0, stress_decid=ihard_stress_decid)
# ---------------------------------------------------------------------------
function _setup_edphys_pfts!()
    npft = 3
    p = CLM.prt_params
    CLM.allocate_prt_params!(p, npft, CLM.num_organ_types, 1)  # 1 leaf age class

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

    # Root profile (Jackson beta, type 1) for set_root_fraction.
    p.fnrt_prof_mode .= 1
    p.fnrt_prof_a    .= 0.976
    p.fnrt_prof_b    .= 0.0
    p.root_long      .= 1.0

    # Leaf-N vertical scaler (used by trim_canopy's decay_coeff_vcmax) and grperc.
    p.leafn_vert_scaler_coeff1 .= 0.0
    p.leafn_vert_scaler_coeff2 .= 1.0
    p.grperc                   .= 0.11

    # Phenology strategy flags per PFT.
    p.season_decid .= [0, CLM.itrue, 0]
    p.stress_decid .= [0, 0, CLM.ihard_stress_decid]
    p.evergreen    .= [CLM.itrue, 0, 0]

    p.leaf_long          .= 1.5
    p.leaf_long_ustory   .= 1.5
    p.leaf_stor_priority .= 0.8

    # Drought phenology (negative thresholds => matric-potential semantics, mm).
    p.phen_drought_threshold .= -152900.0  # below this SMP => abscise
    p.phen_moist_threshold   .= -122000.0
    p.phen_doff_time         .= 100.0
    p.phen_fnrt_drop_fraction .= 0.0
    p.phen_stem_drop_fraction .= 0.0

    p.seed_alloc         .= 0.1
    p.seed_alloc_mature  .= 0.0
    p.dbh_repro_threshold .= 1000.0
    p.repro_alloc_a      .= 0.0
    p.repro_alloc_b      .= 0.0

    # param_derived: branch_frac + canopy-top derived rates.
    pd = CLM.param_derived_type()
    pd.branch_frac = fill(0.25, npft)
    pd.jmax25top = fill(85.0, npft, 1)
    pd.tpu25top  = fill(8.0, npft, 1)
    pd.kp25top   = fill(0.6, npft, 1)
    CLM.ParamDerived[] = pd

    # FATES PFT table.
    ev = CLM.EDPftvarcon_type()
    ev.vcmax25top = fill(50.0, npft, 1)
    ev.damage_recovery_scalar = fill(0.0, npft)
    ev.lf_flab = fill(0.25, npft); ev.lf_fcel = fill(0.50, npft); ev.lf_flig = fill(0.25, npft)
    ev.fr_flab = fill(0.25, npft); ev.fr_fcel = fill(0.50, npft); ev.fr_flig = fill(0.25, npft)
    ev.hgt_min                 = fill(1.5, npft)
    ev.trim_limit              = fill(0.30, npft)
    ev.trim_inc                = fill(0.03, npft)
    ev.phen_cold_size_threshold = fill(0.0, npft)
    ev.phenflush_fraction      = fill(0.5, npft)
    ev.allom_frbstor_repro     = fill(0.0, npft)
    ev.seed_decay_rate         = fill(0.51, npft)
    ev.germination_rate        = fill(0.5, npft)
    ev.seed_dispersal_fraction = fill(0.0, npft)
    ev.repro_frac_seed         = fill(1.0, npft)
    ev.seed_suppl              = fill(0.0, npft)
    ev.prescribed_recruitment  = fill(-1.0, npft)
    CLM.EDPftvarcon_inst[] = ev

    # ed_params: phenology constants + regeneration + VAI bins + logging/q10.
    edp = CLM.ed_params_type()
    edp.regeneration_model = CLM.default_regeneration
    edp.ED_val_phen_a        = -68.0
    edp.ED_val_phen_b        = 638.0
    edp.ED_val_phen_c        = -0.01
    edp.ED_val_phen_chiltemp = 5.0
    edp.ED_val_phen_mindayson = 30.0
    edp.ED_val_phen_ncolddayslim = 5.0
    edp.ED_val_phen_coldtemp = 7.5
    edp.logging_export_frac  = 0.8
    edp.q10_mr               = 1.5
    edp.q10_froz             = 1.5
    nlv = CLM.nlevleaf
    edp.dinc_vai   = fill(1.0, nlv)
    edp.dlower_vai = collect(1.0:1.0:Float64(nlv))
    CLM.EDParams[] = edp

    # SPITFIRE CWD partitioning + decomposition fractions.
    sfp = CLM.sf_params_type()
    sfp.SF_val_CWD_frac = [0.045, 0.075, 0.21, 0.67]
    sfp.SF_val_max_decomp = fill(0.5, CLM.num_fuel_classes)
    CLM.SFParams[] = sfp

    return npft
end

# Carbon PRT object seeded with allometrically-consistent pool carbon at dbh0,
# at a chosen leaf elongation factor (so leaf carbon can be set "on" or "off").
function _seed_prt(ipft::Int, dbh0::Float64; efleaf::Float64=1.0)
    prt = CLM.callom_prt_vartypes()
    CLM.InitPRTVartype!(prt)

    canopy_trim = 1.0; crowndamage = 1
    l2fr = CLM.prt_params.allom_l2fr[ipft]
    effnrt = 1.0 - (1.0 - efleaf) * CLM.prt_params.phen_fnrt_drop_fraction[ipft]
    efstem = 1.0 - (1.0 - efleaf) * CLM.prt_params.phen_stem_drop_fraction[ipft]

    tgt_leaf, _    = CLM.bleaf(dbh0, ipft, crowndamage, canopy_trim, efleaf)
    tgt_fnrt, _    = CLM.bfineroot(dbh0, ipft, canopy_trim, l2fr, effnrt)
    _, tgt_sapw, _ = CLM.bsap_allom(dbh0, ipft, crowndamage, canopy_trim, efstem)
    tgt_store, _   = CLM.bstore_allom(dbh0, ipft, crowndamage, canopy_trim)
    tgt_agw, _     = CLM.bagw_allom(dbh0, ipft, crowndamage, efstem)
    tgt_bgw, _     = CLM.bbgw_allom(dbh0, ipft, efstem)
    tgt_struct, _  = CLM.bdead_allom(tgt_agw, tgt_bgw, tgt_sapw, ipft)

    CLM.SetState!(prt, CLM.leaf_organ,   CLM.carbon12_element, tgt_leaf)
    CLM.SetState!(prt, CLM.fnrt_organ,   CLM.carbon12_element, tgt_fnrt)
    CLM.SetState!(prt, CLM.sapw_organ,   CLM.carbon12_element, tgt_sapw)
    CLM.SetState!(prt, CLM.store_organ,  CLM.carbon12_element, tgt_store)
    CLM.SetState!(prt, CLM.struct_organ, CLM.carbon12_element, tgt_struct)
    CLM.SetState!(prt, CLM.repro_organ,  CLM.carbon12_element, 0.0)

    return prt
end

# A tveg24 running mean fixed at a chosen value (degrees K).
function _make_tveg24(value_k::Float64)
    rmd = CLM.rmean_def_type()
    CLM.define!(rmd, CLM.sec_per_day, CLM.sec_per_day, CLM.moving_ema_window)
    rm = CLM.rmean_type(c_mean=value_k, l_mean=value_k, c_index=1, def_type=rmd)
    return rm
end

# Build a minimal ed_site_type with the phenology + flux arrays touched here.
function _build_site(npft::Int, nlevsoil::Int)
    site = CLM.ed_site_type()
    site.spread = 1.0
    site.lat = 45.0  # Northern Hemisphere

    # Soil layering + root scratch.
    site.nlevsoil     = nlevsoil
    site.zi_soil      = collect(0.0:0.1:0.1 * nlevsoil)   # nlevsoil+1 interfaces
    site.rootfrac_scr = zeros(nlevsoil)

    # Phenology state (start cold / leaves-off, populated memories).
    site.cstatus        = CLM.phen_cstat_iscold
    site.dstatus        = fill(CLM.phen_dstat_moiston, CLM.maxpft)
    site.grow_deg_days  = 0.0
    site.nchilldays     = 10
    site.vegtemp_memory = zeros(CLM.num_vegtemp_mem)
    site.cleafondate    = 1
    site.cleafoffdate   = 1
    site.cndaysleafon   = 0
    site.cndaysleafoff  = 0
    site.dleafondate    = fill(1, CLM.maxpft)
    site.dleafoffdate   = fill(1, CLM.maxpft)
    site.dndaysleafon   = fill(0, CLM.maxpft)
    site.dndaysleafoff  = fill(0, CLM.maxpft)
    site.elong_factor   = fill(1.0, CLM.maxpft)
    site.phen_model_date = 400  # past the first year so the GDD guards behave
    site.liqvol_memory  = zeros(CLM.numWaterMem, CLM.maxpft)
    site.smp_memory     = zeros(CLM.numWaterMem, CLM.maxpft)

    site.use_this_pft   = fill(CLM.itrue, npft)
    site.rec_l2fr       = fill(1.0, CLM.maxpft, CLM.nclmax)
    site.recruitment_rate = zeros(CLM.maxpft)
    site.seed_in        = zeros(npft)
    site.seed_out       = zeros(npft)

    # Per-element mass balance + flux diagnostics.
    site.mass_balance = [CLM.site_massbal_type() for _ in 1:CLM.num_elements[]]
    site.flux_diags = CLM.site_fluxdiags_type()
    site.flux_diags.elem = [CLM.elem_diag_type() for _ in 1:CLM.num_elements[]]
    for el in 1:CLM.num_elements[]
        site.flux_diags.elem[el].surf_fine_litter_input = zeros(npft)
        site.flux_diags.elem[el].root_litter_input      = zeros(npft)
    end

    return site
end

# Build a patch with allocated per-element litter + a fixed tveg24.
function _build_patch(area::Float64, npft::Int, nlevsoil::Int; tveg_k::Float64=283.15)
    patch = CLM.fates_patch_type()
    patch.area = area
    patch.canopy_layer_tlai = zeros(Float64, CLM.nclmax)
    patch.nocomp_pft_label = 1          # not bare ground (so frag scaler runs)
    patch.ncl_p = 1
    patch.land_use_label = CLM.nocomp_bareground_land
    patch.fragmentation_scaler = zeros(nlevsoil)
    patch.tveg24 = _make_tveg24(tveg_k)
    patch.nitr_repro_stoich = fill(0.0, CLM.maxpft)
    patch.phos_repro_stoich = fill(0.0, CLM.maxpft)

    patch.litter = [CLM.litter_type() for _ in 1:CLM.num_elements[]]
    for el in 1:CLM.num_elements[]
        CLM.InitAllocate!(patch.litter[el], npft, nlevsoil, CLM.element_list[el])
        CLM.InitConditions!(patch.litter[el], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    end
    return patch
end

# Link a single patch onto a site (it is both the oldest and youngest).
function _attach_patch!(site, patch)
    site.oldest_patch = patch
    site.youngest_patch = patch
    patch.older = nothing
    patch.younger = nothing
    return nothing
end

# bc_in with the soil fields phenology/recruitment touch.
function _make_bc_in(nlevsoil::Int; liqvol::Float64=0.4, smp::Float64=-50000.0,
                     tempk::Float64=290.0)
    bc = CLM.bc_in_type()
    bc.max_rooting_depth_index_col = nlevsoil
    bc.h2o_liqvol_sl = fill(liqvol, nlevsoil)
    bc.smp_sl        = fill(smp, nlevsoil)
    bc.tempk_sl      = fill(tempk, nlevsoil)
    bc.t_scalar_sisl = fill(1.0, nlevsoil)
    bc.w_scalar_sisl = fill(1.0, nlevsoil)
    bc.z_sisl        = collect(0.05:0.1:0.1 * nlevsoil)
    return bc
end

# Add a cohort node directly to a patch's height-ordered list (single cohort,
# tallest == shortest). Avoids create_cohort here so the leaf elongation factor
# (and thus leaf carbon) is controllable for phenology checks.
function _add_single_cohort!(patch, ipft, dbh0, nn; efleaf=1.0, status=CLM.leaves_on)
    coh = CLM.fates_cohort_type()
    coh.pft = ipft
    coh.n = nn
    coh.dbh = dbh0
    h, _ = CLM.h_allom(dbh0, ipft)
    coh.height = h
    coh.canopy_trim = 1.0
    coh.canopy_layer = 1
    coh.crowndamage = 1
    coh.status_coh = status
    coh.efleaf_coh = efleaf
    coh.effnrt_coh = 1.0
    coh.efstem_coh = 1.0
    coh.l2fr = CLM.prt_params.allom_l2fr[ipft]
    coh.vcmax25top = 50.0
    coh.isnew = false
    coh.dndt = 0.0
    coh.lmort_direct = 0.0
    coh.lmort_collateral = 0.0
    coh.lmort_infra = 0.0
    coh.year_net_uptake = fill(999.0, CLM.nlevleaf)
    coh.prt = _seed_prt(ipft, dbh0; efleaf=efleaf)
    coh.taller = nothing
    coh.shorter = nothing
    patch.tallest = coh
    patch.shortest = coh
    return coh
end

@testset "FATES Batch 13: EDPhysiologyMod" begin

    old_global   = CLM.prt_global[]
    old_numel    = CLM.num_elements[]
    old_ellist   = copy(CLM.element_list)
    old_elpos    = CLM.element_pos[CLM.carbon12_element]
    old_parteh   = CLM.hlm_parteh_mode[]
    old_numpft   = CLM.numpft[]
    old_nleafage = CLM.nleafage[]
    old_doy      = CLM.hlm_day_of_year[]
    old_freqday  = CLM.hlm_freq_day[]
    old_disp     = CLM.hlm_seeddisp_cadence[]
    old_nocomp   = CLM.hlm_use_nocomp[]
    old_luh      = CLM.hlm_use_luh[]
    old_predphys = CLM.hlm_use_ed_prescribed_phys[]
    old_treedmg  = CLM.hlm_use_tree_damage[]
    old_hydro    = CLM.hlm_use_planthydro[]
    old_sfp      = CLM.SFParams[]

    try
        npft = _setup_edphys_pfts!()
        CLM.InitPRTGlobalAllometricCarbon!()

        # Carbon-only element registry.
        CLM.num_elements[] = 1
        empty!(CLM.element_list)
        push!(CLM.element_list, CLM.carbon12_element)
        CLM.element_pos[CLM.carbon12_element] = 1   # carbon12 lives at litter slot 1

        CLM.hlm_parteh_mode[]   = CLM.prt_carbon_allom_hyp
        CLM.numpft[]            = npft
        CLM.nleafage[]          = 1
        CLM.hlm_day_of_year[]   = 100
        CLM.hlm_freq_day[]      = 1.0 / 365.0
        CLM.hlm_seeddisp_cadence[] = CLM.fates_dispersal_cadence_none
        CLM.hlm_use_nocomp[]    = CLM.ifalse
        CLM.hlm_use_luh[]       = CLM.ifalse
        CLM.hlm_use_ed_prescribed_phys[] = CLM.ifalse
        CLM.hlm_use_tree_damage[] = CLM.ifalse
        CLM.hlm_use_planthydro[] = CLM.ifalse

        nlevsoil = 4
        area     = 100.0

        # ==================================================================
        # 1. phenology — COLD DECIDUOUS leaf-on via GDD exceedance.
        #    Site starts iscold (leaves off). Feed warm temps, accumulated GDD,
        #    and a chilling history so the GDD threshold is exceeded -> notcold.
        # ==================================================================
        @testset "phenology: cold-deciduous leaf-on" begin
            site  = _build_site(npft, nlevsoil)
            patch = _build_patch(area, npft, nlevsoil; tveg_k=CLM.t_water_freeze_k_1atm + 20.0)
            _attach_patch!(site, patch)
            bc_in = _make_bc_in(nlevsoil)

            # A PFT2 (cold-deciduous) cohort, currently leaves OFF (efleaf=0).
            coh = _add_single_cohort!(patch, 2, 12.0, 0.05; efleaf=0.0, status=CLM.leaves_off)
            # Zero storage carbon so the leaf-flush PARTEH transfer (which needs a
            # registered leaf stoichiometry, not configured in the carbon-only test
            # PARTEH map) is skipped; we still test the site/cohort state transition.
            CLM.SetState!(coh.prt, CLM.store_organ, CLM.carbon12_element, 0.0)

            # Pre-load the GDD above threshold so leaf-on fires this step.
            # gdd_threshold = -68 + 638*exp(-0.01*nchilldays); at nchilldays=10
            # that is about 509; warm temp keeps accumulating. Set GDD high.
            site.cstatus       = CLM.phen_cstat_iscold
            site.nchilldays    = 10
            site.grow_deg_days = 600.0
            site.cleafoffdate  = 50           # leaf-off was 50 days ago-ish
            site.phen_model_date = 100        # day 101 after advance
            CLM.hlm_day_of_year[] = 100       # not a reset date

            CLM.phenology(site, bc_in)

            @test site.cstatus == CLM.phen_cstat_notcold      # leaves can come on
            @test site.elong_factor[2] == 1.0                 # cold PFT flushed
            @test coh.status_coh == CLM.leaves_on             # cohort flushed
        end

        # ==================================================================
        # 2. phenology — COLD DECIDUOUS leaf-off via cold snap.
        #    Site starts notcold (leaves on); a full cold memory window pushes
        #    ncolddays over the limit -> iscold, and the cohort sheds leaves.
        # ==================================================================
        @testset "phenology: cold-deciduous leaf-off" begin
            site  = _build_site(npft, nlevsoil)
            cold_k = CLM.t_water_freeze_k_1atm - 5.0  # -5 C, below coldtemp(7.5C)
            patch = _build_patch(area, npft, nlevsoil; tveg_k=cold_k)
            _attach_patch!(site, patch)
            bc_in = _make_bc_in(nlevsoil)

            coh = _add_single_cohort!(patch, 2, 12.0, 0.05; efleaf=1.0, status=CLM.leaves_on)
            leaf_c_before = CLM.GetState(coh.prt, CLM.leaf_organ, CLM.carbon12_element)
            @test leaf_c_before > 0.0

            site.cstatus       = CLM.phen_cstat_notcold
            site.cndaysleafon  = 60                  # exceeds mindayson(30)
            site.cleafondate   = 40
            site.phen_model_date = 200
            # Fill the temperature memory with cold so ncolddays exceeds the limit.
            site.vegtemp_memory .= -5.0
            CLM.hlm_day_of_year[] = 200

            CLM.phenology(site, bc_in)

            @test site.cstatus == CLM.phen_cstat_iscold      # leaves off
            @test site.elong_factor[2] == 0.0
            @test (coh.status_coh == CLM.leaves_off || coh.status_coh == CLM.leaves_shedding)
            leaf_c_after = CLM.GetState(coh.prt, CLM.leaf_organ, CLM.carbon12_element)
            @test leaf_c_after < leaf_c_before               # leaves were shed
        end

        # ==================================================================
        # 3. phenology — DROUGHT (hard) DECIDUOUS leaf-off via dry soil.
        #    A hard-deciduous PFT (PFT3) on dry soil (SMP below the drought
        #    threshold) after the minimum on-period should abscise (elong=0).
        # ==================================================================
        @testset "phenology: drought-deciduous leaf-off" begin
            site  = _build_site(npft, nlevsoil)
            patch = _build_patch(area, npft, nlevsoil; tveg_k=CLM.t_water_freeze_k_1atm + 20.0)
            _attach_patch!(site, patch)
            # Very dry: SMP well below phen_drought_threshold (-152900 mm).
            bc_in = _make_bc_in(nlevsoil; smp=-500000.0, liqvol=0.05, tempk=295.0)

            coh = _add_single_cohort!(patch, 3, 12.0, 0.05; efleaf=1.0, status=CLM.leaves_on)

            # Leaves have been on a while (exceed_min_on_period), past spinup.
            site.dstatus[3]      = CLM.phen_dstat_moiston
            site.elong_factor[3] = 1.0
            site.dndaysleafon[3] = 150            # > dleafon_drycheck(100)
            site.dleafondate[3]  = 50
            site.phen_model_date = 200            # > numWaterMem
            # Pre-load dry soil memory so the 10-day mean is below threshold.
            site.smp_memory[:, 3] .= -500000.0
            CLM.hlm_day_of_year[] = 200

            CLM.phenology(site, bc_in)

            @test site.dstatus[3] == CLM.phen_dstat_moistoff  # abscission flagged
            @test site.elong_factor[3] == 0.0                 # leaves dropped
            @test (coh.status_coh == CLM.leaves_off || coh.status_coh == CLM.leaves_shedding)
        end

        # ==================================================================
        # 4. trim_canopy — reduce canopy_trim toward carbon balance.
        #    Give a multi-leaf-layer cohort year_net_uptake values that are far
        #    below the per-layer leaf cost; trim_canopy must reduce canopy_trim.
        # ==================================================================
        @testset "trim_canopy: trims negative-balance leaves" begin
            site  = _build_site(npft, nlevsoil)
            patch = _build_patch(area, npft, nlevsoil)
            _attach_patch!(site, patch)

            # Evergreen PFT1, big enough to have several leaf layers.
            coh = _add_single_cohort!(patch, 1, 40.0, 0.02; efleaf=1.0)
            coh.canopy_trim = 1.0
            # All active leaf layers running a (very) negative carbon balance.
            fill!(coh.year_net_uptake, 0.0)         # net uptake ~ 0 << leaf cost

            trim_before = coh.canopy_trim
            CLM.trim_canopy(site)
            @test coh.canopy_trim < trim_before       # trimmed down
            @test coh.canopy_trim >= CLM.EDPftvarcon_inst[].trim_limit[1]  # not below limit
            # year_net_uptake reset for the next year.
            @test all(coh.year_net_uptake .== 999.0)
        end

        # trim_canopy expands a fully-profitable cohort (no negative layers).
        @testset "trim_canopy: expands profitable canopy" begin
            site  = _build_site(npft, nlevsoil)
            patch = _build_patch(area, npft, nlevsoil)
            _attach_patch!(site, patch)

            coh = _add_single_cohort!(patch, 1, 40.0, 0.02; efleaf=1.0)
            coh.canopy_trim = 0.8
            # No leaf layer had activity (all 999) -> not trimmed -> expand.
            fill!(coh.year_net_uptake, 999.0)

            CLM.trim_canopy(site)
            @test coh.canopy_trim ≈ 0.8 + CLM.EDPftvarcon_inst[].trim_inc[1]
        end

        # ==================================================================
        # 5. recruitment — produce a juvenile cohort of the right size.
        # ==================================================================
        @testset "recruitment: juvenile cohort size" begin
            site  = _build_site(npft, nlevsoil)
            patch = _build_patch(area, npft, nlevsoil)
            _attach_patch!(site, patch)
            patch.tallest = nothing
            patch.shortest = nothing
            bc_in = _make_bc_in(nlevsoil)

            # Only allow PFT1 to recruit; load its germinated-seed pool.
            site.use_this_pft .= CLM.ifalse
            site.use_this_pft[1] = CLM.itrue
            # A generous seed_germ pool [kg/m2] so cohort_n > min_n_safemath.
            patch.litter[1].seed_germ[1] = 1.0

            CLM.recruitment(site, patch, bc_in)

            @test patch.tallest !== nothing
            rec = patch.tallest
            @test rec.pft == 1
            @test rec.height ≈ CLM.EDPftvarcon_inst[].hgt_min[1]
            expected_dbh, _ = CLM.h2d_allom(CLM.EDPftvarcon_inst[].hgt_min[1], 1)
            @test rec.dbh ≈ expected_dbh
            @test rec.n > CLM.min_n_safemath
            @test site.recruitment_rate[1] ≈ rec.n
            # Germinated-seed pool was drawn down by the recruits' mass.
            @test patch.litter[1].seed_germ[1] < 1.0
        end

        # ==================================================================
        # 6. litter / seed mass conservation through the integrate step.
        #    Set known input fluxes and frag fluxes, then verify state changes
        #    are exactly (in - frag) for the seed bank and CWD pools.
        # ==================================================================
        @testset "litter integration: mass conservation" begin
            site  = _build_site(npft, nlevsoil)
            patch = _build_patch(area, npft, nlevsoil)
            _attach_patch!(site, patch)

            litt = patch.litter[1]
            # Seed bank: only an in-flux (no decay, no germination) for PFT1.
            litt.seed[1]           = 0.0
            litt.seed_in_local[1]  = 0.3
            litt.seed_in_extern[1] = 0.2
            litt.seed_decay[1]     = 0.0
            litt.seed_germ_in[1]   = 0.1
            litt.seed_germ[1]      = 0.0
            litt.seed_germ_decay[1] = 0.0
            # AG CWD pool 1: in-flux minus frag.
            litt.ag_cwd[1]      = 1.0
            litt.ag_cwd_in[1]   = 0.4
            litt.ag_cwd_frag[1] = 0.1
            # Leaf fines pool 1.
            litt.leaf_fines[1]      = 2.0
            litt.leaf_fines_in[1]   = 0.5
            litt.leaf_fines_frag[1] = 0.2

            CLM.PreDisturbanceIntegrateLitter(patch)

            # seed[pft] += in_local + in_extern - decay - germ_in
            @test litt.seed[1] ≈ 0.3 + 0.2 - 0.0 - 0.1
            # seed_germ[pft] += germ_in - germ_decay
            @test litt.seed_germ[1] ≈ 0.1 - 0.0
            # ag_cwd += in - frag
            @test litt.ag_cwd[1] ≈ 1.0 + 0.4 - 0.1
            # leaf_fines += in - frag
            @test litt.leaf_fines[1] ≈ 2.0 + 0.5 - 0.2
        end

        # ==================================================================
        # 7. PreDisturbanceLitterFluxes — end-to-end flux generation runs and
        #    conserves: site frag_out equals the patch-area-weighted sum of the
        #    fragmentation fluxes it generated.
        # ==================================================================
        @testset "PreDisturbanceLitterFluxes: frag_out bookkeeping" begin
            site  = _build_site(npft, nlevsoil)
            patch = _build_patch(area, npft, nlevsoil)
            _attach_patch!(site, patch)
            bc_in = _make_bc_in(nlevsoil)

            # A cohort so CWDInput has something to do.
            _add_single_cohort!(patch, 1, 20.0, 0.05; efleaf=1.0)
            # Seed bank + CWD so SeedDecay/CWDOut produce non-trivial fluxes.
            patch.litter[1].seed[1]   = 0.5
            patch.litter[1].ag_cwd   .= 1.0
            patch.litter[1].leaf_fines .= 1.0

            CLM.PreDisturbanceLitterFluxes(site, patch, bc_in)

            litt = patch.litter[1]
            expected_frag = patch.area * (
                sum(litt.ag_cwd_frag) + sum(litt.bg_cwd_frag) +
                sum(litt.leaf_fines_frag) + sum(litt.root_fines_frag) +
                sum(litt.seed_decay) + sum(litt.seed_germ_decay))
            @test site.mass_balance[1].frag_out ≈ expected_frag
            @test site.mass_balance[1].frag_out > 0.0
            # Fragmentation scaler was set (HLM scalars = 1.0 here).
            @test all(patch.fragmentation_scaler .≈ 1.0)
        end

    finally
        CLM.prt_global[]        = old_global
        CLM.num_elements[]      = old_numel
        empty!(CLM.element_list); append!(CLM.element_list, old_ellist)
        CLM.element_pos[CLM.carbon12_element] = old_elpos
        CLM.hlm_parteh_mode[]   = old_parteh
        CLM.numpft[]            = old_numpft
        CLM.nleafage[]          = old_nleafage
        CLM.hlm_day_of_year[]   = old_doy
        CLM.hlm_freq_day[]      = old_freqday
        CLM.hlm_seeddisp_cadence[] = old_disp
        CLM.hlm_use_nocomp[]    = old_nocomp
        CLM.hlm_use_luh[]       = old_luh
        CLM.hlm_use_ed_prescribed_phys[] = old_predphys
        CLM.hlm_use_tree_damage[] = old_treedmg
        CLM.hlm_use_planthydro[] = old_hydro
        CLM.SFParams[]          = old_sfp
    end
end
