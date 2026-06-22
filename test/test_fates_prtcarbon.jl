# test_fates_prtcarbon.jl
# Tests for FATES Batch 7 (Tier F): PRTAllometricCarbonMod — the carbon-only,
# allometric-growth PARTEH allocation hypothesis.
#
# Strategy:
#   * Build a synthetic single-PFT prt_params table using the (well-tested) PFT1
#     allometry configuration from the allometry suite (obrien height + saldarriaga
#     AGB/leaf, constant sapwood/fineroot/storage).
#   * Install param_derived (branch_frac) and ed_params (default regeneration).
#   * Call InitPRTGlobalAllometricCarbon!() to register the 6-variable carbon
#     hypothesis and install it as the global descriptor.
#   * Allocate a concrete callom_prt_vartypes via the generic InitAllocate!, seed
#     known leaf/fnrt/sapw/store/struct/repro carbon for an evergreen PFT, wire the
#     boundary conditions (dbh, carbon_balance, pft, ctrim, lstat, cdamage, elong),
#     and run DailyPRT! phases 1->2->3 with a small POSITIVE daily carbon gain.
#
# Assertions:
#   * Carbon mass conservation: total allocated C (sum over pools, incl. repro)
#     equals the carbon-gain input (carbon_balance fully spent to ~calloc_abs_error).
#   * Pools move toward (not past, within tolerance) their allometric targets.
#   * Storage behaves per the target logic (storage receives carbon when below
#     target / when shedding).

using Test
using CLM

# Configure a single woody evergreen PFT with the PFT1 allometry modes.
function _setup_prtcarbon_pft!()
    npft = 1
    p = CLM.prt_params
    # one param-file organ slot is enough here (carbon hypothesis); 1 leaf age class
    CLM.allocate_prt_params!(p, npft, CLM.num_organ_types, 1)

    # Common physically reasonable defaults.
    p.c2b          .= 2.0
    p.wood_density .= 0.6
    p.slatop       .= 0.012
    p.slamax       .= 0.020
    p.allom_agb_frac .= 0.6
    p.allom_dbh_maxheight .= 90.0       # > min_max_dbh_for_trees (tree)
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

    # Allometry modes (PFT1 family from the allometry suite).
    p.allom_dmode .= 1
    p.woody       .= 1
    p.allom_cmode .= 1
    p.allom_smode .= 1
    p.allom_fmode .= 1
    p.allom_stmode .= 1
    p.allom_hmode .= 1
    p.allom_amode .= 1
    p.allom_lmode .= 1

    # Height (obrien): h = 10^(log10(d)*p1 + p2)
    p.allom_d2h1 .= 0.64
    p.allom_d2h2 .= 0.37
    p.allom_d2h3 .= -0.034

    # AGB (saldarriaga) params
    p.allom_agb1 .= 0.06896; p.allom_agb2 .= 0.572; p.allom_agb3 .= 1.94; p.allom_agb4 .= 0.931
    # Leaf (saldarriaga)
    p.allom_d2bl1 .= 0.07; p.allom_d2bl2 .= 1.3; p.allom_d2bl3 .= 0.55

    # Phenology / allocation params used by DailyPRT.
    p.stress_decid .= 0
    p.season_decid .= 0
    p.evergreen    .= CLM.itrue
    p.leaf_stor_priority .= 0.8
    p.leaf_long          .= 1.5
    p.leaf_long_ustory   .= 1.5

    # Reproductive allocation (default regeneration path).
    p.seed_alloc         .= 0.1
    p.seed_alloc_mature  .= 0.0
    p.dbh_repro_threshold .= 1000.0   # below threshold -> seed_alloc only
    p.repro_alloc_a      .= 0.0
    p.repro_alloc_b      .= 0.0

    # param_derived (branch_frac) for bsap.
    pd = CLM.param_derived_type()
    pd.branch_frac = fill(0.25, npft)
    CLM.ParamDerived[] = pd

    # ed_params: default regeneration + damage/vai bins (unused here but required).
    edp = CLM.ed_params_type()
    edp.regeneration_model = CLM.default_regeneration
    edp.ED_val_history_damage_bin_edges = [0.0, 20.0, 50.0, 80.0]
    nlv = CLM.nlevleaf
    edp.dinc_vai   = fill(1.0, nlv)
    edp.dlower_vai = collect(1.0:1.0:Float64(nlv))
    CLM.EDParams[] = edp

    return npft
end

@testset "FATES Batch 7: PRTAllometricCarbonMod" begin

    old_global = CLM.prt_global[]
    _setup_prtcarbon_pft!()
    ipft = 1

    # --- Register the carbon hypothesis global descriptor. ---------------------
    CLM.InitPRTGlobalAllometricCarbon!()
    g = CLM.prt_global[]
    @test g !== nothing
    @test g.hyp_id == CLM.prt_carbon_allom_hyp
    @test g.num_vars == 6
    @test g.num_bc_in == 7
    @test g.num_bc_inout == 2
    @test g.num_bc_out == 0
    # leaf carbon registered with 1 age position (matches leaf_long dim2 == 1)
    @test g.state_descriptor[CLM.ac_leaf_c_id].num_pos == 1
    # prt_global and prt_global_ac point at the same descriptor
    @test CLM.prt_global_ac[] === g

    # --- Build a concrete callom_prt_vartypes and allocate it. -----------------
    prt = CLM.callom_prt_vartypes()
    CLM.InitPRTVartype!(prt)   # allocate + nan-init + flush BCs

    # Determine allometric targets at dbh0 so we can seed pools deliberately
    # BELOW target (so phase 2/3 growth has something to do).
    dbh0        = 10.0
    canopy_trim = 1.0
    crowndamage = 1
    leaf_status = CLM.leaves_on
    ef          = 1.0
    l2fr        = CLM.prt_params.allom_l2fr[ipft]

    tgt_leaf, _   = CLM.bleaf(dbh0, ipft, crowndamage, canopy_trim, ef)
    tgt_fnrt, _   = CLM.bfineroot(dbh0, ipft, canopy_trim, l2fr, ef)
    _, tgt_sapw, _ = CLM.bsap_allom(dbh0, ipft, crowndamage, canopy_trim, ef)
    tgt_store, _  = CLM.bstore_allom(dbh0, ipft, crowndamage, canopy_trim)
    tgt_agw, _    = CLM.bagw_allom(dbh0, ipft, crowndamage, ef)
    tgt_bgw, _    = CLM.bbgw_allom(dbh0, ipft, ef)
    tgt_struct, _ = CLM.bdead_allom(tgt_agw, tgt_bgw, tgt_sapw, ipft)

    @test tgt_leaf  > 0.0
    @test tgt_sapw  > 0.0
    @test tgt_struct > 0.0

    # Seed pools slightly below target so allocation pushes them up.
    seed_frac = 0.90
    CLM.SetState!(prt, CLM.leaf_organ,   CLM.carbon12_element, seed_frac * tgt_leaf)
    CLM.SetState!(prt, CLM.fnrt_organ,   CLM.carbon12_element, seed_frac * tgt_fnrt)
    CLM.SetState!(prt, CLM.sapw_organ,   CLM.carbon12_element, seed_frac * tgt_sapw)
    CLM.SetState!(prt, CLM.store_organ,  CLM.carbon12_element, seed_frac * tgt_store)
    CLM.SetState!(prt, CLM.struct_organ, CLM.carbon12_element, seed_frac * tgt_struct)
    CLM.SetState!(prt, CLM.repro_organ,  CLM.carbon12_element, 0.0)

    # Zero rates and snapshot val0 (sets up mass-conservation bookkeeping).
    CLM.ZeroRates!(prt)

    # Set a small per-day leaf/fineroot turnover so phase-1 demand is exercised.
    prt.variables[CLM.ac_leaf_c_id].turnover[1] = 0.001 * tgt_leaf
    prt.variables[CLM.ac_fnrt_c_id].turnover[1] = 0.001 * tgt_fnrt

    # --- Wire boundary conditions. --------------------------------------------
    dbh_ref = Ref(dbh0)
    # Positive daily carbon gain (net NPP). Make it generous so phase-3 growth runs.
    carbon_gain = 0.05 * (tgt_leaf + tgt_fnrt + tgt_sapw + tgt_struct + tgt_store)
    @test carbon_gain > CLM.calloc_abs_error
    cb_ref  = Ref(carbon_gain)

    CLM.RegisterBCInout!(prt, CLM.ac_bc_inout_id_dbh;   bc_rval = dbh_ref)
    CLM.RegisterBCInout!(prt, CLM.ac_bc_inout_id_netdc; bc_rval = cb_ref)

    CLM.RegisterBCIn!(prt, CLM.ac_bc_in_id_pft;     bc_ival = Ref(ipft))
    CLM.RegisterBCIn!(prt, CLM.ac_bc_in_id_ctrim;   bc_rval = Ref(canopy_trim))
    CLM.RegisterBCIn!(prt, CLM.ac_bc_in_id_lstat;   bc_ival = Ref(leaf_status))
    CLM.RegisterBCIn!(prt, CLM.ac_bc_in_id_cdamage; bc_ival = Ref(crowndamage))
    CLM.RegisterBCIn!(prt, CLM.ac_bc_in_id_efleaf;  bc_rval = Ref(ef))
    CLM.RegisterBCIn!(prt, CLM.ac_bc_in_id_effnrt;  bc_rval = Ref(ef))
    CLM.RegisterBCIn!(prt, CLM.ac_bc_in_id_efstem;  bc_rval = Ref(ef))

    # Snapshot total carbon mass before allocation.
    total_c_before = CLM.GetState(prt, CLM.leaf_organ,   CLM.carbon12_element) +
                     CLM.GetState(prt, CLM.fnrt_organ,   CLM.carbon12_element) +
                     CLM.GetState(prt, CLM.sapw_organ,   CLM.carbon12_element) +
                     CLM.GetState(prt, CLM.store_organ,  CLM.carbon12_element) +
                     CLM.GetState(prt, CLM.struct_organ, CLM.carbon12_element) +
                     CLM.GetState(prt, CLM.repro_organ,  CLM.carbon12_element)
    dbh_before = dbh_ref[]

    # --- Run the daily allocation in its three phases. ------------------------
    CLM.DailyPRT!(prt, 1)
    CLM.DailyPRT!(prt, 2)
    CLM.DailyPRT!(prt, 3)

    # carbon_balance should be (essentially) fully spent.
    @test abs(cb_ref[]) <= CLM.calloc_abs_error * 10

    # Total carbon mass after.
    total_c_after = CLM.GetState(prt, CLM.leaf_organ,   CLM.carbon12_element) +
                    CLM.GetState(prt, CLM.fnrt_organ,   CLM.carbon12_element) +
                    CLM.GetState(prt, CLM.sapw_organ,   CLM.carbon12_element) +
                    CLM.GetState(prt, CLM.store_organ,  CLM.carbon12_element) +
                    CLM.GetState(prt, CLM.struct_organ, CLM.carbon12_element) +
                    CLM.GetState(prt, CLM.repro_organ,  CLM.carbon12_element)

    # --- Mass conservation: all allocated C == the daily carbon gain. ---------
    allocated = total_c_after - total_c_before
    @test isapprox(allocated, carbon_gain; atol = 1e-7, rtol = 1e-7)

    # No pool went NaN/Inf.
    for org in (CLM.leaf_organ, CLM.fnrt_organ, CLM.sapw_organ,
                CLM.store_organ, CLM.struct_organ, CLM.repro_organ)
        @test isfinite(CLM.GetState(prt, org, CLM.carbon12_element))
    end

    # --- Pools moved toward (and not meaningfully past) their allometric targets.
    leaf_after   = CLM.GetState(prt, CLM.leaf_organ,   CLM.carbon12_element)
    fnrt_after   = CLM.GetState(prt, CLM.fnrt_organ,   CLM.carbon12_element)
    sapw_after   = CLM.GetState(prt, CLM.sapw_organ,   CLM.carbon12_element)
    struct_after = CLM.GetState(prt, CLM.struct_organ, CLM.carbon12_element)

    # Each below-target pool should have increased toward its target.
    @test leaf_after   >= seed_frac * tgt_leaf
    @test fnrt_after   >= seed_frac * tgt_fnrt
    @test sapw_after   >= seed_frac * tgt_sapw
    @test struct_after >= seed_frac * tgt_struct

    # Pools should not exceed their (possibly grown) target by more than tolerance
    # at the starting dbh's targets minus growth; check leaf does not blow past
    # its starting target hugely (it can grow with dbh, so allow a margin).
    @test leaf_after <= tgt_leaf * 1.5

    # dbh should have grown (carbon was available for phase-3 growth).
    @test dbh_ref[] >= dbh_before

    # --- net_alloc bookkeeping mirrors the pool changes (mass conservation). ---
    na_total = 0.0
    for org in (CLM.leaf_organ, CLM.fnrt_organ, CLM.sapw_organ,
                CLM.store_organ, CLM.struct_organ, CLM.repro_organ)
        na_total += CLM.GetNetAlloc(prt, org, CLM.carbon12_element)
    end
    @test isapprox(na_total, carbon_gain; atol = 1e-7, rtol = 1e-7)

    # --- Storage-target behavior: a SHEDDING plant stashes all gain to storage.
    prt2 = CLM.callom_prt_vartypes()
    CLM.InitPRTVartype!(prt2)
    CLM.SetState!(prt2, CLM.leaf_organ,   CLM.carbon12_element, tgt_leaf)
    CLM.SetState!(prt2, CLM.fnrt_organ,   CLM.carbon12_element, tgt_fnrt)
    CLM.SetState!(prt2, CLM.sapw_organ,   CLM.carbon12_element, tgt_sapw)
    CLM.SetState!(prt2, CLM.store_organ,  CLM.carbon12_element, tgt_store)
    CLM.SetState!(prt2, CLM.struct_organ, CLM.carbon12_element, tgt_struct)
    CLM.SetState!(prt2, CLM.repro_organ,  CLM.carbon12_element, 0.0)
    CLM.ZeroRates!(prt2)

    store_before = CLM.GetState(prt2, CLM.store_organ, CLM.carbon12_element)
    gain2 = 0.01 * tgt_leaf
    cb2_ref  = Ref(gain2)
    dbh2_ref = Ref(dbh0)

    CLM.RegisterBCInout!(prt2, CLM.ac_bc_inout_id_dbh;   bc_rval = dbh2_ref)
    CLM.RegisterBCInout!(prt2, CLM.ac_bc_inout_id_netdc; bc_rval = cb2_ref)
    CLM.RegisterBCIn!(prt2, CLM.ac_bc_in_id_pft;     bc_ival = Ref(ipft))
    CLM.RegisterBCIn!(prt2, CLM.ac_bc_in_id_ctrim;   bc_rval = Ref(canopy_trim))
    CLM.RegisterBCIn!(prt2, CLM.ac_bc_in_id_lstat;   bc_ival = Ref(CLM.leaves_shedding))
    CLM.RegisterBCIn!(prt2, CLM.ac_bc_in_id_cdamage; bc_ival = Ref(crowndamage))
    CLM.RegisterBCIn!(prt2, CLM.ac_bc_in_id_efleaf;  bc_rval = Ref(ef))
    CLM.RegisterBCIn!(prt2, CLM.ac_bc_in_id_effnrt;  bc_rval = Ref(ef))
    CLM.RegisterBCIn!(prt2, CLM.ac_bc_in_id_efstem;  bc_rval = Ref(ef))

    CLM.DailyPRT!(prt2, 3)
    store_after = CLM.GetState(prt2, CLM.store_organ, CLM.carbon12_element)
    # All of the gain went to storage (shedding stash branch).
    @test isapprox(store_after - store_before, gain2; atol = 1e-9, rtol = 1e-9)
    @test abs(cb2_ref[]) <= CLM.calloc_abs_error

    # Restore the previous global so we don't leak into other suites.
    CLM.prt_global[] = old_global
end
