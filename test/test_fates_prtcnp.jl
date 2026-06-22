# test_fates_prtcnp.jl
# Tests for FATES Batch 7 (Tier F): PRTAllometricCNPMod — the Carbon-Nitrogen-
# Phosphorus prioritized allometric allocation hypothesis.
#
# We build ONE synthetic woody PFT with a complete allometry + stoichiometry
# parameter set, register the CNP hypothesis global, construct a concrete
# `cnp_allom_prt_vartypes` plant with known C/N/P pools, then drive `DailyPRT!`
# with a given C gain and N/P uptake. We assert:
#   (1) C, N and P mass conservation across the allocation (the full PARTEH
#       balance: input == Δstate + efflux (+ respired excess for C) + Δturnover);
#   (2) the plant grows toward its allometric/stoichiometric target (stature
#       growth: dbh increases) when carbon and nutrients are available;
#   (3) nutrient-limited downregulation: when phosphorus is scarce the plant
#       grows less (stature growth is suppressed) and more carbon/nitrogen is
#       disposed of, relative to the ample-nutrient case;
#   (4) excess carbon is respired/effluxed when storage cannot absorb it (the
#       storage-overflow "burn" path).
#
# These are physics-behaviour assertions (conservation + direction), not a
# byte-for-byte Fortran comparison.

using Test
using CLM

# -------------------------------------------------------------------------------------
# Build a synthetic single-PFT prt_params table tuned for a small woody plant.
# -------------------------------------------------------------------------------------
function _setup_cnp_prt_params!()
    npft     = 1
    norgan   = 6
    nleafage = 1

    p = CLM.prt_params
    CLM.allocate_prt_params!(p, npft, norgan, nleafage)

    # Allometry: simple power-law modes (no height needed for amode 2 AGB).
    p.c2b          .= 2.0
    p.wood_density .= 0.6
    p.slatop       .= 0.012
    p.slamax       .= 0.020
    p.allom_agb_frac .= 0.6
    p.allom_dbh_maxheight .= 90.0       # > min_max_dbh_for_trees -> "tree"
    p.allom_sai_scaler    .= 0.1
    p.cushion             .= 1.0
    p.woody       .= 1
    p.evergreen   .= 1                   # allow leaf replacement of turnover

    p.allom_hmode .= 1
    p.allom_d2h1 .= 0.64;  p.allom_d2h2 .= 0.37;  p.allom_d2h3 .= -0.034

    # AGB: 2-parameter power law  bagw = p1*d^p2/c2b  (amode 2, no height dep).
    p.allom_amode .= 2
    p.allom_agb1 .= 0.1;  p.allom_agb2 .= 2.4

    p.allom_cmode .= 1                   # bbgw constant proportionality
    p.allom_smode .= 1                   # sapwood from leaf area
    p.allom_fmode .= 1                   # fineroot proportional to trimmed bleaf
    p.allom_stmode .= 1                  # storage = cushion * trimmed bleaf

    # Leaf: 2-parameter power  blmax = p1*d^p2/c2b (lmode 2).
    p.allom_lmode .= 2
    p.allom_d2bl1 .= 0.07;  p.allom_d2bl2 .= 1.3;  p.allom_d2bl3 .= 0.55

    p.leafn_vert_scaler_coeff1 .= 0.00963
    p.leafn_vert_scaler_coeff2 .= 2.43
    p.fnrt_prof_mode .= 1.0;  p.fnrt_prof_a .= 0.976;  p.fnrt_prof_b .= 2.0
    p.allom_l2fr .= 1.0
    p.allom_la_per_sa_int .= 0.8e3;  p.allom_la_per_sa_slp .= 0.0

    # Turnover.
    p.leaf_long  .= 1.5
    p.leaf_long_ustory .= 1.5
    p.root_long  .= 1.0
    p.branch_long .= 50.0

    # Organ mapping for the allocation routines.
    # organ_id (param-file index -> global organ id); pick the natural order.
    p.organ_id .= [CLM.leaf_organ, CLM.fnrt_organ, CLM.sapw_organ,
                   CLM.store_organ, CLM.repro_organ, CLM.struct_organ]
    # organ_param_id (global organ id -> param-file column). Here the param column
    # equals the position of that organ in organ_id.
    p.organ_param_id .= 0
    for (col, gid) in enumerate(p.organ_id)
        p.organ_param_id[gid] = col
    end

    # Allocation priority (PFT x param-organ). Priority 1 = leaf/fnrt (replace
    # turnover first), priority 3 = sapw/struct, priority 4 = repro; storage uses
    # its hard-coded priority level 2. Indexed in the organ_id column ordering.
    #              leaf fnrt sapw store repro struct
    p.alloc_priority[1, :] .= [1.0, 1.0, 3.0, 4.0, 4.0, 3.0]
    p.leaf_stor_priority .= 1.0

    # Stoichiometry (N:C and P:C per organ, indexed by param column).
    # leaf richest, wood poorest.
    #                          leaf  fnrt  sapw  store repro struct
    p.nitr_stoich_p1[1, :] .= [0.040, 0.030, 0.004, 0.010, 0.020, 0.002]
    p.phos_stoich_p1[1, :] .= [0.004, 0.003, 0.0004,0.001, 0.002, 0.0002]

    p.nitr_store_ratio .= 2.0
    p.phos_store_ratio .= 2.0
    p.store_ovrflw_frac .= 1.0          # default: generous storage overflow

    # Reproduction (default regeneration so it does not require the TRS).
    p.seed_alloc         .= 0.1
    p.seed_alloc_mature  .= 0.0
    p.dbh_repro_threshold .= 1.0e6      # always "immature" -> seed_alloc only
    p.repro_alloc_a      .= 0.0
    p.repro_alloc_b      .= 0.0

    # PID controller for l2fr. Default OFF (inert) so l2fr is held fixed, which
    # keeps the fine-root turnover deterministic for the conservation checks.
    # Individual testsets re-enable it.
    p.pid_kp .= 0.0;  p.pid_ki .= 0.0;  p.pid_kd .= 0.0

    return p
end

# Build the per-organ derived param (branch_frac) and ed_params globals.
function _setup_cnp_globals!()
    npft = 1
    pd = CLM.param_derived_type()
    pd.branch_frac = fill(0.25, npft)
    CLM.ParamDerived[] = pd

    edp = CLM.ed_params_type()
    edp.ED_val_history_damage_bin_edges = [0.0, 20.0, 50.0, 80.0]
    nlv = CLM.nlevleaf
    edp.dinc_vai   = fill(1.0, nlv)
    edp.dlower_vai = collect(1.0:1.0:Float64(nlv))
    # Predictive (not prescribed) N and P uptake so efflux/limitation are exercised.
    edp.n_uptake_mode      = 2
    edp.p_uptake_mode      = 2
    edp.regeneration_model = CLM.default_regeneration
    CLM.EDParams[] = edp
    return edp
end

# Construct a fresh CNP plant, allocate, and set state to allometric targets
# scaled by `cscale` (carbon), with nutrients at the C target's growth-min
# stoichiometry. `turnfrac` is the fraction of leaf+fnrt pools seeded as turnover
# (paid first in the replacement step). Returns the plant plus the Ref BCs.
function _build_cnp_plant(; dbh0=2.0, ipft=1, c_gain=1.0e-3, n_gain=1.0e-4,
                          p_gain=1.0e-5, cscale=1.0, turnfrac=1.0e-3)
    plant = CLM.cnp_allom_prt_vartypes()
    CLM.InitPRTVartype!(plant)

    canopy_trim = 1.0
    elong       = 1.0
    crown_damage = 1
    leaf_status = CLM.leaves_on

    # Allometric C targets at dbh0 (used to seed the pools on-allometry).
    _, sapw_t, _ = CLM.bsap_allom(dbh0, ipft, crown_damage, canopy_trim, elong)
    agw_t, _ = CLM.bagw_allom(dbh0, ipft, crown_damage, elong)
    bgw_t, _ = CLM.bbgw_allom(dbh0, ipft, elong)
    struct_t, _ = CLM.bdead_allom(agw_t, bgw_t, sapw_t, ipft)
    leaf_t, _   = CLM.bleaf(dbh0, ipft, crown_damage, canopy_trim, elong)
    fnrt_t, _   = CLM.bfineroot(dbh0, ipft, canopy_trim, 1.0, elong)
    store_t, _  = CLM.bstore_allom(dbh0, ipft, crown_damage, canopy_trim)

    p = CLM.prt_params
    opid = p.organ_param_id

    cpools = Dict(CLM.leaf_organ => leaf_t, CLM.fnrt_organ => fnrt_t,
                  CLM.sapw_organ => sapw_t, CLM.store_organ => store_t,
                  CLM.struct_organ => struct_t, CLM.repro_organ => 0.0)

    for gid in (CLM.leaf_organ, CLM.fnrt_organ, CLM.sapw_organ,
                CLM.store_organ, CLM.repro_organ, CLM.struct_organ)
        cval = cpools[gid] * cscale
        CLM.SetState!(plant, gid, CLM.carbon12_element, cval)
        ncoef = p.nitr_stoich_p1[ipft, opid[gid]]
        pcoef = p.phos_stoich_p1[ipft, opid[gid]]
        CLM.SetState!(plant, gid, CLM.nitrogen_element, cval * ncoef)
        CLM.SetState!(plant, gid, CLM.phosphorus_element, cval * pcoef)
    end

    # Reset rates (val0 := val), then seed a small turnover on leaf+fnrt so the
    # prioritized-replacement path has demand.
    CLM.ZeroRates!(plant)
    for gid in (CLM.leaf_organ, CLM.fnrt_organ)
        for el in (CLM.carbon12_element, CLM.nitrogen_element, CLM.phosphorus_element)
            iv = CLM.sp_organ_map_get(CLM.get_prt_global(), gid, el)
            plant.variables[iv].turnover[1] = turnfrac * plant.variables[iv].val[1]
        end
    end

    # Boundary conditions (Refs).
    r_dbh   = Ref(dbh0)
    r_resp  = Ref(0.0)
    r_l2fr  = Ref(1.0)
    r_netdn = Ref(n_gain)
    r_netdp = Ref(p_gain)
    r_cxint = Ref(0.0)
    r_cx0   = Ref(0.0)
    r_emad  = Ref(0.0)

    CLM.RegisterBCInout!(plant, CLM.acnp_bc_inout_id_dbh;         bc_rval=r_dbh)
    CLM.RegisterBCInout!(plant, CLM.acnp_bc_inout_id_resp_excess; bc_rval=r_resp)
    CLM.RegisterBCInout!(plant, CLM.acnp_bc_inout_id_l2fr;        bc_rval=r_l2fr)
    CLM.RegisterBCInout!(plant, CLM.acnp_bc_inout_id_netdn;       bc_rval=r_netdn)
    CLM.RegisterBCInout!(plant, CLM.acnp_bc_inout_id_netdp;       bc_rval=r_netdp)
    CLM.RegisterBCInout!(plant, CLM.acnp_bc_inout_id_cx_int;      bc_rval=r_cxint)
    CLM.RegisterBCInout!(plant, CLM.acnp_bc_inout_id_cx0;         bc_rval=r_cx0)
    CLM.RegisterBCInout!(plant, CLM.acnp_bc_inout_id_emadcxdt;    bc_rval=r_emad)

    r_pft    = Ref(ipft)
    r_ctrim  = Ref(canopy_trim)
    r_lstat  = Ref(leaf_status)
    r_netdc  = Ref(c_gain)
    r_ncrep  = Ref(p.nitr_stoich_p1[ipft, opid[CLM.leaf_organ]])
    r_pcrep  = Ref(p.phos_stoich_p1[ipft, opid[CLM.leaf_organ]])
    r_cdam   = Ref(crown_damage)
    r_efl    = Ref(elong)
    r_eff    = Ref(elong)
    r_efs    = Ref(elong)

    CLM.RegisterBCIn!(plant, CLM.acnp_bc_in_id_pft;      bc_ival=r_pft)
    CLM.RegisterBCIn!(plant, CLM.acnp_bc_in_id_ctrim;    bc_rval=r_ctrim)
    CLM.RegisterBCIn!(plant, CLM.acnp_bc_in_id_lstat;    bc_ival=r_lstat)
    CLM.RegisterBCIn!(plant, CLM.acnp_bc_in_id_netdc;    bc_rval=r_netdc)
    CLM.RegisterBCIn!(plant, CLM.acnp_bc_in_id_nc_repro; bc_rval=r_ncrep)
    CLM.RegisterBCIn!(plant, CLM.acnp_bc_in_id_pc_repro; bc_rval=r_pcrep)
    CLM.RegisterBCIn!(plant, CLM.acnp_bc_in_id_cdamage;  bc_ival=r_cdam)
    CLM.RegisterBCIn!(plant, CLM.acnp_bc_in_id_efleaf;   bc_rval=r_efl)
    CLM.RegisterBCIn!(plant, CLM.acnp_bc_in_id_effnrt;   bc_rval=r_eff)
    CLM.RegisterBCIn!(plant, CLM.acnp_bc_in_id_efstem;   bc_rval=r_efs)

    r_cef = Ref(0.0);  r_nef = Ref(0.0);  r_pef = Ref(0.0);  r_lim = Ref(0)
    CLM.RegisterBCOut!(plant, CLM.acnp_bc_out_id_cefflux; bc_rval=r_cef)
    CLM.RegisterBCOut!(plant, CLM.acnp_bc_out_id_nefflux; bc_rval=r_nef)
    CLM.RegisterBCOut!(plant, CLM.acnp_bc_out_id_pefflux; bc_rval=r_pef)
    CLM.RegisterBCOut!(plant, CLM.acnp_bc_out_id_limiter; bc_ival=r_lim)

    bcs = (dbh=r_dbh, resp=r_resp, l2fr=r_l2fr, netdn=r_netdn, netdp=r_netdp,
           cefflux=r_cef, nefflux=r_nef, pefflux=r_pef, limiter=r_lim)
    inputs = (c_gain=c_gain, n_gain=n_gain, p_gain=p_gain)
    return plant, bcs, inputs
end

# Total mass of `element` summed over all organs/positions.
function _total_mass(plant, element)
    tot = 0.0
    for gid in (CLM.leaf_organ, CLM.fnrt_organ, CLM.sapw_organ,
                CLM.store_organ, CLM.repro_organ, CLM.struct_organ)
        tot += CLM.GetState(plant, gid, element)
    end
    return tot
end

# Total turnover (accumulated this control period) of `element` over all organs.
function _total_turnover(plant, element)
    g = CLM.get_prt_global()
    tot = 0.0
    for gid in (CLM.leaf_organ, CLM.fnrt_organ, CLM.sapw_organ,
                CLM.store_organ, CLM.repro_organ, CLM.struct_organ)
        iv = CLM.sp_organ_map_get(g, gid, element)
        tot += sum(plant.variables[iv].turnover)
    end
    return tot
end

@testset "FATES Batch 7: PRTAllometricCNPMod (CNP allocation)" begin

    _setup_cnp_prt_params!()
    _setup_cnp_globals!()
    CLM.InitPRTGlobalAllometricCNP()

    @testset "hypothesis registration" begin
        g = CLM.get_prt_global()
        @test g.hyp_id == CLM.prt_cnp_flex_allom_hyp
        @test g.num_vars == CLM.acnp_num_vars
        @test g.num_bc_in == CLM.acnp_num_bc_in
        @test g.num_bc_out == CLM.acnp_num_bc_out
        @test g.num_bc_inout == CLM.acnp_num_bc_inout
        @test CLM.prt_global_acnp[] === g
        @test CLM.sp_organ_map_get(g, CLM.leaf_organ, CLM.carbon12_element) == CLM.leaf_c_id
        @test CLM.sp_organ_map_get(g, CLM.struct_organ, CLM.phosphorus_element) == CLM.struct_p_id
    end

    # --------------------------------------------------------------------------------
    # (1) Mass conservation, including the dynamic l2fr/TrimFineRoot path.
    #     The full PARTEH balance over the control period is:
    #       gain == Δstate + efflux (+ respired excess for C) + Δturnover
    # --------------------------------------------------------------------------------
    @testset "mass conservation (C, N, P) with PID l2fr active" begin
        # Turn the l2fr PID controller ON so the fine-root contraction / forceful
        # turnover (TrimFineRoot!) path is exercised, then verify conservation
        # still closes once turnover is accounted for.
        CLM.prt_params.pid_kp .= 1.0
        plant, bcs, inp = _build_cnp_plant(; dbh0=2.0, c_gain=2.0e-3,
                                           n_gain=1.0e-2, p_gain=1.0e-3,
                                           cscale=1.0, turnfrac=5.0e-2)

        c0 = _total_mass(plant, CLM.carbon12_element)
        n0 = _total_mass(plant, CLM.nitrogen_element)
        p0 = _total_mass(plant, CLM.phosphorus_element)
        tc0 = _total_turnover(plant, CLM.carbon12_element)
        tn0 = _total_turnover(plant, CLM.nitrogen_element)
        tp0 = _total_turnover(plant, CLM.phosphorus_element)

        CLM.DailyPRT!(plant, 1)

        c1 = _total_mass(plant, CLM.carbon12_element)
        n1 = _total_mass(plant, CLM.nitrogen_element)
        p1 = _total_mass(plant, CLM.phosphorus_element)
        tc1 = _total_turnover(plant, CLM.carbon12_element)
        tn1 = _total_turnover(plant, CLM.nitrogen_element)
        tp1 = _total_turnover(plant, CLM.phosphorus_element)

        @test isapprox(inp.c_gain,
                       (c1 - c0) + bcs.cefflux[] + bcs.resp[] + (tc1 - tc0);
                       atol=10 * CLM.calloc_abs_error)
        @test isapprox(inp.n_gain,
                       (n1 - n0) + bcs.nefflux[] + (tn1 - tn0);
                       atol=10 * CLM.calloc_abs_error)
        @test isapprox(inp.p_gain,
                       (p1 - p0) + bcs.pefflux[] + (tp1 - tp0);
                       atol=10 * CLM.calloc_abs_error)

        # No negative pools.
        for gid in (CLM.leaf_organ, CLM.fnrt_organ, CLM.sapw_organ,
                    CLM.store_organ, CLM.repro_organ, CLM.struct_organ)
            for el in (CLM.carbon12_element, CLM.nitrogen_element, CLM.phosphorus_element)
                @test CLM.GetState(plant, gid, el) >= -CLM.calloc_abs_error
            end
        end

        CLM.prt_params.pid_kp .= 0.0   # restore inert controller
    end

    # --------------------------------------------------------------------------------
    # (2) Growth toward target: ample C and nutrients -> stature growth (dbh up).
    # --------------------------------------------------------------------------------
    @testset "stature growth toward allometric target" begin
        # Pools start ON allometry (cscale=1) with small turnover; surplus carbon
        # plus ample nutrients drive stature growth.
        plant, bcs, inp = _build_cnp_plant(; dbh0=2.0, c_gain=5.0e-2,
                                           n_gain=5.0e-1, p_gain=5.0e-2,
                                           cscale=1.0, turnfrac=1.0e-3)
        dbh_before = bcs.dbh[]

        # Reproductive pool starts empty; after growth it should receive carbon.
        repro_c0 = CLM.GetState(plant, CLM.repro_organ, CLM.carbon12_element)

        CLM.DailyPRT!(plant, 1)

        # dbh strictly increased -> the plant grew its stature.
        @test bcs.dbh[] > dbh_before
        # Carbon was the binding constraint (nutrients ample).
        @test bcs.limiter[] == CLM.c_limited
        # Reproductive allocation occurred during growth.
        @test CLM.GetState(plant, CLM.repro_organ, CLM.carbon12_element) > repro_c0
        # The new dbh is consistent with the allometric leaf target moving up.
        new_leaf_target, _ = CLM.bleaf(bcs.dbh[], 1, 1, 1.0, 1.0)
        old_leaf_target, _ = CLM.bleaf(dbh_before, 1, 1, 1.0, 1.0)
        @test new_leaf_target > old_leaf_target
    end

    # --------------------------------------------------------------------------------
    # (3) Nutrient-limited downregulation: scarce P suppresses growth.
    # --------------------------------------------------------------------------------
    @testset "phosphorus-limited downregulation" begin
        # Same C and N; contrast ample vs. scarce phosphorus.
        pa, ba, ia = _build_cnp_plant(; dbh0=2.0, c_gain=5.0e-2, n_gain=5.0e-1,
                                      p_gain=5.0e-2, cscale=1.0, turnfrac=1.0e-3)
        ps, bs, is_ = _build_cnp_plant(; dbh0=2.0, c_gain=5.0e-2, n_gain=5.0e-1,
                                       p_gain=1.0e-9, cscale=1.0, turnfrac=1.0e-3)

        # Conservation closes in both runs.
        for (pl, b, inp) in ((pa, ba, ia), (ps, bs, is_))
            c0 = _total_mass(pl, CLM.carbon12_element)
            n0 = _total_mass(pl, CLM.nitrogen_element)
            p0 = _total_mass(pl, CLM.phosphorus_element)
            tc0 = _total_turnover(pl, CLM.carbon12_element)
            tn0 = _total_turnover(pl, CLM.nitrogen_element)
            tp0 = _total_turnover(pl, CLM.phosphorus_element)
            CLM.DailyPRT!(pl, 1)
            c1 = _total_mass(pl, CLM.carbon12_element)
            n1 = _total_mass(pl, CLM.nitrogen_element)
            p1 = _total_mass(pl, CLM.phosphorus_element)
            tc1 = _total_turnover(pl, CLM.carbon12_element)
            tn1 = _total_turnover(pl, CLM.nitrogen_element)
            tp1 = _total_turnover(pl, CLM.phosphorus_element)
            @test isapprox(inp.c_gain, (c1 - c0) + b.cefflux[] + b.resp[] + (tc1 - tc0);
                           atol=10 * CLM.calloc_abs_error)
            @test isapprox(inp.n_gain, (n1 - n0) + b.nefflux[] + (tn1 - tn0);
                           atol=10 * CLM.calloc_abs_error)
            @test isapprox(inp.p_gain, (p1 - p0) + b.pefflux[] + (tp1 - tp0);
                           atol=10 * CLM.calloc_abs_error)
        end

        # (3) With ample P the plant grew stature; with scarce P growth is
        # suppressed (P is consumed by replacement before stature growth, which
        # then cannot proceed). The scarce-P plant therefore grows strictly less.
        @test ba.dbh[] > bs.dbh[]
        @test isapprox(bs.dbh[], 2.0; atol=1e-12)        # no stature growth
        @test ba.limiter[] == CLM.c_limited              # ample-P -> carbon limits
        @test bs.limiter[] != CLM.c_limited              # scarce-P -> growth skipped
        # The under-utilized (because un-balanceable) nitrogen is effluxed more
        # heavily in the scarce-P case.
        @test bs.nefflux[] > ba.nefflux[]
    end

    # --------------------------------------------------------------------------------
    # (4) Excess-carbon disposal via the storage-overflow "burn" path.
    # --------------------------------------------------------------------------------
    @testset "excess carbon respired (storage overflow burn)" begin
        # No storage overflow room + scarce nutrients (so no biomass can absorb
        # the carbon) + a large carbon gain -> the excess is respired as
        # resp_excess (store_c_overflow == burn_c_store_overflow).
        CLM.prt_params.store_ovrflw_frac .= 0.0
        plant, bcs, inp = _build_cnp_plant(; dbh0=2.0, c_gain=5.0e-1,
                                           n_gain=1.0e-3, p_gain=1.0e-4,
                                           cscale=1.0, turnfrac=1.0e-3)

        c0 = _total_mass(plant, CLM.carbon12_element)
        tc0 = _total_turnover(plant, CLM.carbon12_element)

        CLM.DailyPRT!(plant, 1)

        c1 = _total_mass(plant, CLM.carbon12_element)
        tc1 = _total_turnover(plant, CLM.carbon12_element)

        # Conservation: the carbon input is split between biomass, storage,
        # turnover, efflux and respired excess.
        @test isapprox(inp.c_gain,
                       (c1 - c0) + bcs.cefflux[] + bcs.resp[] + (tc1 - tc0);
                       atol=20 * CLM.calloc_abs_error)

        # The bulk of the carbon could not be allocated and was respired.
        @test bcs.resp[] > 0.0
        @test bcs.resp[] > 0.5 * inp.c_gain

        CLM.prt_params.store_ovrflw_frac .= 1.0   # restore default
    end
end
