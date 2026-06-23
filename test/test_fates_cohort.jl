# test_fates_cohort.jl
# Tests for FATES Batch 8 (Tier F): FatesCohortMod — the COHORT type, the
# fundamental FATES demographic unit.
#
# Strategy:
#   * Build a synthetic single-PFT prt_params table (PFT1 allometry config from
#     the allometry/PRT-carbon suites), install param_derived (with non-trivial
#     jmax/tpu/kp), ed_params (default regeneration + vai bins), and the FATES
#     PFT table (EDPftvarcon_inst.vcmax25top).
#   * Register the carbon-only PARTEH hypothesis and set the interface control
#     flags (parteh mode = carbon, no SP, no plant hydro, nleafage = 1).
#   * Create a cohort with a freshly allocated callom_prt_vartypes carbon object,
#     plus an attached ed_cohort_hydr_type, via Create().
#   * Assert field initialization (structure / classes / biophysical rates / new
#     flag / BC wiring).
#   * Copy the cohort into a second cohort and assert scalar-state deep-equality.
#   * Exercise ZeroValues and assert the special initializations.
#   * Build a 3-cohort taller/shorter linked list and assert the ordering links.

using Test
using CLM

# Configure a single woody evergreen PFT with the PFT1 allometry modes, the
# carbon hypothesis global, and the derived/PFT param tables.
function _setup_cohort_pft!()
    npft = 1
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

    p.stress_decid .= 0
    p.season_decid .= 0
    p.evergreen    .= CLM.itrue
    p.leaf_stor_priority .= 0.8
    p.leaf_long          .= 1.5
    p.leaf_long_ustory   .= 1.5

    p.seed_alloc         .= 0.1
    p.seed_alloc_mature  .= 0.0
    p.dbh_repro_threshold .= 1000.0
    p.repro_alloc_a      .= 0.0
    p.repro_alloc_b      .= 0.0

    # param_derived: branch_frac (used by bsap) + canopy-top derived rates.
    pd = CLM.param_derived_type()
    pd.branch_frac = fill(0.25, npft)
    pd.jmax25top = fill(85.0, npft, 1)
    pd.tpu25top  = fill(8.0, npft, 1)
    pd.kp25top   = fill(0.6, npft, 1)
    CLM.ParamDerived[] = pd

    # FATES PFT table: vcmax25top (numpft x nleafage).
    ev = CLM.EDPftvarcon_type()
    ev.vcmax25top = fill(50.0, npft, 1)
    CLM.EDPftvarcon_inst[] = ev

    # ed_params: default regeneration + vai bins (required by tree_sai).
    edp = CLM.ed_params_type()
    edp.regeneration_model = CLM.default_regeneration
    edp.ED_val_history_damage_bin_edges = [0.0, 20.0, 50.0, 80.0]
    # Size- and cohort-age class binning (used by sizetype/coagetype indexing).
    edp.ED_val_history_sizeclass_bin_edges  = [0.0, 5.0, 10.0, 20.0, 50.0]
    edp.ED_val_history_coageclass_bin_edges = [0.0, 3.0, 6.0, 9.0, 12.0]
    nlv = CLM.nlevleaf
    edp.dinc_vai   = fill(1.0, nlv)
    edp.dlower_vai = collect(1.0:1.0:Float64(nlv))
    CLM.EDParams[] = edp

    # Number of size / cohort-age class bins (history dimensions).
    CLM.nlevsclass[] = length(edp.ED_val_history_sizeclass_bin_edges)
    CLM.nlevcoage[]  = length(edp.ED_val_history_coageclass_bin_edges)

    return npft
end

# Build a freshly allocated carbon PRT object seeded with physically reasonable
# pool carbon for a cohort at the given dbh.
function _build_seeded_prt(ipft::Int, dbh0::Float64)
    prt = CLM.callom_prt_vartypes()
    CLM.InitPRTVartype!(prt)

    canopy_trim = 1.0
    crowndamage = 1
    ef          = 1.0
    l2fr        = CLM.prt_params.allom_l2fr[ipft]

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
    CLM.SetState!(prt, CLM.store_organ,  CLM.carbon12_element, tgt_store)
    CLM.SetState!(prt, CLM.struct_organ, CLM.carbon12_element, tgt_struct)
    CLM.SetState!(prt, CLM.repro_organ,  CLM.carbon12_element, 0.0)

    return prt
end

@testset "FATES Batch 8: FatesCohortMod" begin

    # Preserve and restore module-global state mutated by this suite.
    old_global   = CLM.prt_global[]
    old_parteh   = CLM.hlm_parteh_mode[]
    old_use_sp   = CLM.hlm_use_sp[]
    old_hydro    = CLM.hlm_use_planthydro[]
    old_nleafage = CLM.nleafage[]
    old_nlevsc   = CLM.nlevsclass[]
    old_nlevca   = CLM.nlevcoage[]

    try
        _setup_cohort_pft!()
        ipft = 1

        CLM.InitPRTGlobalAllometricCarbon!()

        # Interface control flags: carbon hypothesis, full (non-SP) physiology,
        # no plant hydraulics, one leaf age class.
        CLM.hlm_parteh_mode[]    = CLM.prt_carbon_allom_hyp
        CLM.hlm_use_sp[]         = CLM.ifalse
        CLM.hlm_use_planthydro[] = CLM.ifalse
        CLM.nleafage[]           = 1

        dbh0  = 10.0
        nn    = 0.2
        coage = 0.0
        clayer      = 1
        crowndamage = 1
        status      = CLM.leaves_on
        ctrim       = 1.0
        spread      = 1.0
        ef          = 1.0
        can_tlai    = zeros(Float64, CLM.nclmax)

        # height from PFT allometry so it is internally consistent.
        height, _ = CLM.h_allom(dbh0, ipft)

        # ----------------------------------------------------------------------
        # Create a cohort.
        # ----------------------------------------------------------------------
        prt    = _build_seeded_prt(ipft, dbh0)
        cohort = CLM.fates_cohort_type()
        CLM.Create(cohort, prt, ipft, nn, height, coage, dbh0, status, ctrim,
                   0.0, clayer, crowndamage, spread, can_tlai, ef, ef, ef)

        # --- Field initialization -------------------------------------------
        @test cohort.prt === prt
        @test cohort.isnew == true
        @test cohort.pft == ipft
        @test cohort.n == nn
        @test cohort.dbh == dbh0
        @test cohort.height == height
        @test cohort.coage == coage
        @test cohort.canopy_layer == clayer
        @test cohort.canopy_layer_yesterday == Float64(clayer)
        @test cohort.crowndamage == crowndamage
        @test cohort.status_coh == status
        @test cohort.canopy_trim == ctrim
        @test cohort.efleaf_coh == ef
        @test cohort.effnrt_coh == ef
        @test cohort.efstem_coh == ef
        @test cohort.l2fr == CLM.prt_params.allom_l2fr[ipft]

        # Pointers nulled by Create/Init.
        @test cohort.taller === nothing
        @test cohort.shorter === nothing
        @test cohort.co_hydr === nothing

        # Size/age classification populated.
        @test cohort.size_class >= 1
        @test cohort.size_by_pft_class >= 1
        @test cohort.coage_class >= 1
        @test cohort.coage_by_pft_class >= 1

        # Biophysical rates: single leaf age class -> exactly the PFT/derived value.
        @test cohort.vcmax25top ≈ CLM.edpftvarcon_inst().vcmax25top[ipft, 1]
        @test cohort.jmax25top  ≈ CLM.param_derived().jmax25top[ipft, 1]
        @test cohort.tpu25top   ≈ CLM.param_derived().tpu25top[ipft, 1]
        @test cohort.kp25top    ≈ CLM.param_derived().kp25top[ipft, 1]

        # Crown area / LAI / SAI computed and physically reasonable.
        @test cohort.c_area > 0.0
        @test cohort.treelai > 0.0
        @test cohort.treesai > 0.0

        # ZeroValues-driven fields.
        @test cohort.g_sb_laweight == 0.0
        @test cohort.daily_n_demand == -9.0
        @test cohort.daily_p_demand == -9.0
        @test all(cohort.year_net_uptake .== 999.0)
        @test all(cohort.ts_net_uptake .== 0.0)

        # BCs wired into the PRT object (carbon hypothesis -> 7 in / 2 inout).
        @test cohort.prt.bc_in[CLM.ac_bc_in_id_pft].ival[] == ipft
        @test cohort.prt.bc_inout[CLM.ac_bc_inout_id_dbh].rval[] == dbh0
        @test cohort.prt.bc_in[CLM.ac_bc_in_id_ctrim].rval[] == ctrim

        # CanUpperUnder diagnostic.
        @test CLM.CanUpperUnder(cohort) == CLM.ican_upper
        cohort.canopy_layer = 2
        @test CLM.CanUpperUnder(cohort) == CLM.ican_ustory
        cohort.canopy_layer = clayer

        # ----------------------------------------------------------------------
        # Copy into a second cohort and assert deep-equality of scalar state.
        # ----------------------------------------------------------------------
        # Set some additional scalar state to verify it copies.
        cohort.gpp_acc_hold = 1.23
        cohort.npp_acc_hold = 4.56
        cohort.resp_acc_hold = 7.89
        cohort.cmort = 0.01
        cohort.bmort = 0.02

        prt2   = _build_seeded_prt(ipft, dbh0)
        cohort2 = CLM.fates_cohort_type()
        CLM.Init(cohort2, prt2)   # allocate the copy's PRT carrier
        CLM.Copy(cohort, cohort2)

        @test cohort2.indexnumber == CLM.fates_unset_int
        @test cohort2.taller === nothing
        @test cohort2.shorter === nothing
        # The copy keeps its own PRT object (not the source's).
        @test cohort2.prt === prt2
        @test cohort2.prt !== cohort.prt

        for f in (:pft, :n, :dbh, :coage, :height, :canopy_layer,
                  :canopy_layer_yesterday, :crowndamage, :status_coh, :canopy_trim,
                  :efleaf_coh, :effnrt_coh, :efstem_coh, :l2fr, :c_area, :treelai,
                  :treesai, :isnew, :size_class, :size_by_pft_class, :coage_class,
                  :coage_by_pft_class, :vcmax25top, :jmax25top, :tpu25top, :kp25top,
                  :gpp_acc_hold, :npp_acc_hold, :resp_acc_hold, :cmort, :bmort,
                  :daily_n_demand, :daily_p_demand)
            @test getfield(cohort2, f) == getfield(cohort, f)
        end
        @test cohort2.year_net_uptake == cohort.year_net_uptake
        @test cohort2.ts_net_uptake == cohort.ts_net_uptake
        # PRT state copied element-wise.
        @test CLM.GetState(cohort2.prt, CLM.leaf_organ, CLM.carbon12_element) ≈
              CLM.GetState(cohort.prt, CLM.leaf_organ, CLM.carbon12_element)

        # ----------------------------------------------------------------------
        # ZeroValues! resets the accountable fields.
        # ----------------------------------------------------------------------
        CLM.ZeroValues(cohort2)
        @test cohort2.g_sb_laweight == 0.0
        @test cohort2.nv == 0
        @test cohort2.status_coh == 0
        @test cohort2.treesai == 0.0
        @test cohort2.size_class == 1
        @test cohort2.coage_class == 1
        @test cohort2.gpp_tstep == 0.0
        @test cohort2.daily_n_demand == -9.0
        @test all(cohort2.year_net_uptake .== 999.0)
        @test all(cohort2.ts_net_uptake .== 0.0)
        @test cohort2.fire_mort == 0.0
        # acc_hold deliberately left untouched by ZeroValues (keeps the copied
        # value, NOT zeroed/NaNed), matching Fortran.
        @test cohort2.gpp_acc_hold == 1.23

        # ----------------------------------------------------------------------
        # Build a 3-cohort taller/shorter linked list and assert ordering.
        # cA (tallest) <-> cB (mid) <-> cC (shortest)
        # ----------------------------------------------------------------------
        cA = CLM.fates_cohort_type(height = 30.0)
        cB = CLM.fates_cohort_type(height = 20.0)
        cC = CLM.fates_cohort_type(height = 10.0)

        # taller points to the next-tallest, shorter to the next-shorter.
        cA.shorter = cB
        cB.taller  = cA
        cB.shorter = cC
        cC.taller  = cB

        # Top of list has no taller; bottom has no shorter.
        @test cA.taller === nothing
        @test cC.shorter === nothing

        # Walk down (taller -> shorter) and check descending heights / identity.
        @test cA.shorter === cB
        @test cA.shorter.shorter === cC
        @test cA.height > cA.shorter.height > cA.shorter.shorter.height

        # Walk up (shorter -> taller) and check ascending heights / identity.
        @test cC.taller === cB
        @test cC.taller.taller === cA
        @test cC.height < cC.taller.height < cC.taller.taller.height

        # ----------------------------------------------------------------------
        # FreeMemory drops the PRT/hydraulics references.
        # ----------------------------------------------------------------------
        CLM.FreeMemory(cohort2)
        @test cohort2.prt === nothing

    finally
        CLM.prt_global[]          = old_global
        CLM.hlm_parteh_mode[]     = old_parteh
        CLM.hlm_use_sp[]          = old_use_sp
        CLM.hlm_use_planthydro[]  = old_hydro
        CLM.nleafage[]            = old_nleafage
        CLM.nlevsclass[]          = old_nlevsc
        CLM.nlevcoage[]           = old_nlevca
    end
end
