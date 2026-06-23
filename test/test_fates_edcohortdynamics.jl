# test_fates_edcohortdynamics.jl
# Tests for FATES Batch 12 (Tier F): EDCohortDynamicsMod — the demographic
# ENGINE that manages the cohort lifecycle on the patch's height-ordered
# `taller`/`shorter` linked list.
#
# Strategy (mirrors test_fates_cohort.jl's setup):
#   * Build a synthetic single-PFT (woody, carbon-only) prt_params table, install
#     param_derived, the FATES PFT table, and ed_params (regeneration + bins +
#     fusion tolerances + max_cohort_per_patch).
#   * Register the carbon-only PARTEH hypothesis; flags = carbon, no SP, no plant
#     hydro, no cohort-age-tracking, nleafage = 1.
#   * Build a minimal-but-complete ed_site_type (termination/growth tallies,
#     flux_diags, soil layering + root scratch) and a patch with allocated litter.
#   * Exercise:
#       - InitPRTObject!  -> allocates a callom carbon PRT object.
#       - create_cohort   -> inserts cohorts; assert height ordering + count.
#       - sort_cohorts    -> re-sorts a scrambled list.
#       - count_cohorts   -> forward/backward symmetric count.
#       - fuse_cohorts    -> two near-identical cohorts merge; carbon + number
#                            conserved across fusion.
#       - terminate_cohorts -> a sub-threshold cohort is removed and the list
#                            relinked; its mass lands in litter.
#       - copy/zero        -> via the cohort Copy/ZeroValues on a list node.

using Test
using CLM

# Configure a single woody evergreen carbon PFT (PFT1 allometry), like the cohort
# suite, but additionally set the EDCohortDynamics-relevant params.
function _setup_edcohort_pft!()
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

    # Root profile (Jackson beta, type 1) for set_root_fraction (SendCohortToLitter).
    p.fnrt_prof_mode .= 1
    p.fnrt_prof_a    .= 0.976
    p.fnrt_prof_b    .= 0.0

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

    # FATES PFT table: vcmax25top + damage_recovery_scalar.
    ev = CLM.EDPftvarcon_type()
    ev.vcmax25top = fill(50.0, npft, 1)
    ev.damage_recovery_scalar = fill(0.0, npft)
    # Litter decomposability fractions (leaf + fine root) for GetDecompyFrac.
    ev.lf_flab = fill(0.25, npft); ev.lf_fcel = fill(0.50, npft); ev.lf_flig = fill(0.25, npft)
    ev.fr_flab = fill(0.25, npft); ev.fr_fcel = fill(0.50, npft); ev.fr_flig = fill(0.25, npft)
    CLM.EDPftvarcon_inst[] = ev

    # ed_params: regeneration + bins + fusion tolerances + max cohorts.
    edp = CLM.ed_params_type()
    edp.regeneration_model = CLM.default_regeneration
    edp.ED_val_history_damage_bin_edges = [0.0, 20.0, 50.0, 80.0]
    edp.ED_val_history_sizeclass_bin_edges  = [0.0, 5.0, 10.0, 20.0, 50.0]
    edp.ED_val_history_coageclass_bin_edges = [0.0, 3.0, 6.0, 9.0, 12.0]
    edp.ED_val_cohort_size_fusion_tol = 0.06
    edp.ED_val_cohort_age_fusion_tol  = 0.08
    edp.max_cohort_per_patch          = 100
    nlv = CLM.nlevleaf
    edp.dinc_vai   = fill(1.0, nlv)
    edp.dlower_vai = collect(1.0:1.0:Float64(nlv))
    CLM.EDParams[] = edp

    CLM.nlevsclass[] = length(edp.ED_val_history_sizeclass_bin_edges)
    CLM.nlevcoage[]  = length(edp.ED_val_history_coageclass_bin_edges)

    # SPITFIRE CWD partitioning fractions (twig / small branch / large branch /
    # trunk) used by SendCohortToLitter -> adjust_SF_CWD_frac. Must sum to 1.
    sfp = CLM.sf_params_type()
    sfp.SF_val_CWD_frac = [0.045, 0.075, 0.21, 0.67]
    CLM.SFParams[] = sfp

    return npft
end

# Carbon PRT object seeded with allometrically-consistent pool carbon at dbh0.
function _seed_prt(ipft::Int, dbh0::Float64)
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
    CLM.SetState!(prt, CLM.store_organ,  CLM.carbon12_element, tgt_store)
    CLM.SetState!(prt, CLM.struct_organ, CLM.carbon12_element, tgt_struct)
    CLM.SetState!(prt, CLM.repro_organ,  CLM.carbon12_element, 0.0)

    return prt
end

# Build a minimal ed_site_type with the arrays the termination / fusion paths
# touch (sized by [n_term_mort_types, nlevsclass, npft] / [nlevsclass, npft] etc.).
function _build_site(npft::Int, nlevsoil::Int)
    site = CLM.ed_site_type()
    nmt = CLM.n_term_mort_types
    nsc = CLM.nlevsclass[]

    site.term_nindivs_canopy    = zeros(nmt, nsc, npft)
    site.term_nindivs_ustory    = zeros(nmt, nsc, npft)
    site.term_carbonflux_canopy = zeros(nmt, npft)
    site.term_carbonflux_ustory = zeros(nmt, npft)
    site.term_abg_flux          = zeros(nsc, npft)
    site.growthflux_fusion      = zeros(nsc, npft)
    site.spread                 = 1.0

    # Soil layering + root scratch (SendCohortToLitter).
    site.nlevsoil     = nlevsoil
    site.zi_soil      = collect(0.0:0.1:0.1 * nlevsoil)   # nlevsoil+1 interfaces
    site.rootfrac_scr = zeros(nlevsoil)

    # Flux diagnostics: one elem_diag per element, with the allocatable inputs.
    site.flux_diags = CLM.site_fluxdiags_type()
    site.flux_diags.elem = [CLM.elem_diag_type() for _ in 1:CLM.num_elements[]]
    for el in 1:CLM.num_elements[]
        site.flux_diags.elem[el].surf_fine_litter_input = zeros(npft)
        site.flux_diags.elem[el].root_litter_input      = zeros(npft)
    end

    return site
end

# Build a patch with allocated per-element litter (so SendCohortToLitter has pools).
function _build_patch(area::Float64, npft::Int, nlevsoil::Int)
    patch = CLM.fates_patch_type()
    patch.area = area
    patch.canopy_layer_tlai = zeros(Float64, CLM.nclmax)

    patch.litter = [CLM.litter_type() for _ in 1:CLM.num_elements[]]
    for el in 1:CLM.num_elements[]
        CLM.InitAllocate!(patch.litter[el], npft, nlevsoil, CLM.element_list[el])
        CLM.InitConditions!(patch.litter[el], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    end
    return patch
end

# Count cohorts walking shortest -> taller.
function _list_count(patch)
    n = 0
    c = patch.shortest
    while c !== nothing
        n += 1
        c = c.taller
    end
    return n
end

# Sum carbon over the live pools of all cohorts (number-weighted, kgC).
function _patch_live_c(patch)
    tot = 0.0
    c = patch.tallest
    while c !== nothing
        for org in (CLM.leaf_organ, CLM.fnrt_organ, CLM.sapw_organ,
                    CLM.store_organ, CLM.struct_organ, CLM.repro_organ)
            tot += c.n * CLM.GetState(c.prt, org, CLM.carbon12_element)
        end
        c = c.shorter
    end
    return tot
end

@testset "FATES Batch 12: EDCohortDynamicsMod" begin

    old_global   = CLM.prt_global[]
    old_numel    = CLM.num_elements[]
    old_ellist   = copy(CLM.element_list)
    old_parteh   = CLM.hlm_parteh_mode[]
    old_use_sp   = CLM.hlm_use_sp[]
    old_hydro    = CLM.hlm_use_planthydro[]
    old_agetrk   = CLM.hlm_use_cohort_age_tracking[]
    old_nleafage = CLM.nleafage[]
    old_nlevsc   = CLM.nlevsclass[]
    old_nlevca   = CLM.nlevcoage[]
    old_sfp      = CLM.SFParams[]

    try
        _setup_edcohort_pft!()
        ipft = 1
        CLM.InitPRTGlobalAllometricCarbon!()

        # Carbon-only element registry (num_elements = 1, [carbon12_element]).
        CLM.num_elements[] = 1
        empty!(CLM.element_list)
        push!(CLM.element_list, CLM.carbon12_element)

        CLM.hlm_parteh_mode[]            = CLM.prt_carbon_allom_hyp
        CLM.hlm_use_sp[]                 = CLM.ifalse
        CLM.hlm_use_planthydro[]         = CLM.ifalse
        CLM.hlm_use_cohort_age_tracking[] = CLM.ifalse
        CLM.nleafage[]                   = 1

        nlevsoil = 3
        area     = 100.0
        spread   = 1.0
        status   = CLM.leaves_on
        ctrim    = 1.0
        ef       = 1.0
        coage    = 0.0
        clayer   = 1
        crowndmg = 1

        site  = _build_site(ipft, nlevsoil)
        patch = _build_patch(area, ipft, nlevsoil)

        # ------------------------------------------------------------------
        # InitPRTObject! allocates a carbon PRT object.
        # ------------------------------------------------------------------
        obj = CLM.InitPRTObject!()
        @test obj isa CLM.callom_prt_vartypes

        # ------------------------------------------------------------------
        # create_cohort: insert three cohorts of DIFFERENT sizes, in an order
        # that is NOT height-sorted (small, big, mid). The list must come out
        # height-ordered regardless.
        # ------------------------------------------------------------------
        function _add_cohort!(dbh0, nn)
            prt = _seed_prt(ipft, dbh0)
            h, _ = CLM.h_allom(dbh0, ipft)
            CLM.create_cohort(site, patch, ipft, nn, h, coage, dbh0, prt,
                              ef, ef, ef, status, 0, ctrim, 0.0, clayer, crowndmg,
                              spread)
        end

        _add_cohort!(5.0,  0.20)   # small
        _add_cohort!(30.0, 0.05)   # big
        _add_cohort!(15.0, 0.10)   # mid

        # The cohorts are flagged new by Create; mark them non-new so the level-2
        # termination + fusion flux branches (and is-new fusion gate) behave like
        # established cohorts.
        c = patch.tallest
        while c !== nothing
            c.isnew = false
            c = c.shorter
        end

        # --- List invariants: ordering, head/tail, count ---
        @test patch.tallest !== nothing
        @test patch.shortest !== nothing
        @test patch.tallest.taller === nothing
        @test patch.shortest.shorter === nothing

        # Descending height walking tallest -> shorter.
        heights = Float64[]
        c = patch.tallest
        while c !== nothing
            push!(heights, c.height)
            c = c.shorter
        end
        @test length(heights) == 3
        @test issorted(heights; rev=true)

        # count_cohorts sets countcohorts and is forward/backward symmetric.
        CLM.count_cohorts(patch)
        @test patch.countcohorts == 3
        @test _list_count(patch) == 3

        # ------------------------------------------------------------------
        # sort_cohorts: scramble the height field of the middle cohort so the
        # list is briefly out of order, then re-sort.
        # ------------------------------------------------------------------
        # Grab the middle (shorter of tallest) and bump its height above tallest.
        mid = patch.tallest.shorter
        saved_h = mid.height
        mid.height = patch.tallest.height + 5.0
        CLM.sort_cohorts(patch)
        # After sort, list is height-ordered again with `mid` now at the top.
        @test patch.tallest === mid
        heights2 = Float64[]
        c = patch.tallest
        while c !== nothing
            push!(heights2, c.height)
            c = c.shorter
        end
        @test issorted(heights2; rev=true)
        # restore for the remaining tests
        mid.height = saved_h
        CLM.sort_cohorts(patch)
        @test issorted([patch.tallest.height,
                        patch.tallest.shorter.height,
                        patch.tallest.shorter.shorter.height]; rev=true)

        # ------------------------------------------------------------------
        # fuse_cohorts: add a second cohort that is essentially identical to the
        # mid one (dbh within the size fusion tolerance). After fusion the patch
        # should hold one FEWER cohort, with number density and total live carbon
        # conserved.
        # ------------------------------------------------------------------
        n_before  = _list_count(patch)
        c_before  = _patch_live_c(patch)
        # total number density before
        ndens_before = 0.0
        c = patch.tallest
        while c !== nothing; ndens_before += c.n; c = c.shorter; end

        _add_cohort!(15.3, 0.10)   # ~within 0.06 rel of dbh=15 -> fuses with mid
        c = patch.tallest
        while c !== nothing; c.isnew = false; c = c.shorter; end

        @test _list_count(patch) == n_before + 1   # added one
        CLM.fuse_cohorts(site, patch)
        @test _list_count(patch) == n_before        # fused back down

        c_after = _patch_live_c(patch)
        ndens_after = 0.0
        c = patch.tallest
        while c !== nothing; ndens_after += c.n; c = c.shorter; end

        # number density conserved exactly; total live carbon conserved closely
        # (cohort fusion is by design only approximately carbon-conserving for the
        # allometric pools, but the WeightedFuse of the PRT object conserves the
        # per-area pool mass to round-off).
        @test ndens_after ≈ (ndens_before + 0.10) rtol=1e-12
        @test c_after ≈ (c_before + 0.10 *
              sum(CLM.GetState(_seed_prt(ipft, 15.3), org, CLM.carbon12_element)
                  for org in (CLM.leaf_organ, CLM.fnrt_organ, CLM.sapw_organ,
                              CLM.store_organ, CLM.struct_organ, CLM.repro_organ))) rtol=2e-2

        # list still ordered + symmetric after fusion (fuse calls sort)
        CLM.count_cohorts(patch)
        @test _list_count(patch) == patch.countcohorts

        # ------------------------------------------------------------------
        # terminate_cohorts: drive one cohort's number density below min_npm2 and
        # confirm it is removed (level 2), the neighbors relink, and its mass is
        # deposited into the litter pools (mass conservation: live -> litter).
        # ------------------------------------------------------------------
        n_pre_term = _list_count(patch)

        # Pick the shortest cohort and starve its number density.
        victim = patch.shortest
        victim.isnew = false
        # Per-m2 density just under min_npm2 (-> terminates), but with enough
        # absolute number that its biomass produces a visible litter deposit.
        victim.n     = CLM.min_npm2 * area * 0.5

        # Litter AG-CWD mass before (carbon element = 1).
        litt_before = sum(patch.litter[1].ag_cwd) + sum(patch.litter[1].bg_cwd) +
                      sum(patch.litter[1].leaf_fines) + sum(patch.litter[1].root_fines)

        CLM.terminate_cohorts(site, patch, 2, 1)

        @test _list_count(patch) == n_pre_term - 1     # one removed
        # The removed cohort is gone from the list (no node still has its identity).
        c = patch.tallest
        found_victim = false
        while c !== nothing
            c === victim && (found_victim = true)
            c = c.shorter
        end
        @test !found_victim

        # list head/tail still well-formed
        @test patch.tallest.taller === nothing
        @test patch.shortest.shorter === nothing
        CLM.count_cohorts(patch)
        @test _list_count(patch) == patch.countcohorts

        # Litter increased (cohort biomass -> fragmenting pools).
        litt_after = sum(patch.litter[1].ag_cwd) + sum(patch.litter[1].bg_cwd) +
                     sum(patch.litter[1].leaf_fines) + sum(patch.litter[1].root_fines)
        @test litt_after > litt_before

        # Site-level termination tally recorded the removed individuals.
        @test sum(site.term_nindivs_canopy) ≈ victim.n rtol=1e-12

        # ------------------------------------------------------------------
        # insert_cohort directly: build a fresh 2-node list and splice a 3rd in
        # the middle by height.
        # ------------------------------------------------------------------
        p2 = _build_patch(area, ipft, nlevsoil)
        ca = CLM.fates_cohort_type(height = 30.0, n = 0.1, dbh = 30.0, pft = ipft)
        cc = CLM.fates_cohort_type(height = 10.0, n = 0.1, dbh = 10.0, pft = ipft)
        # insert tall then short
        CLM.insert_cohort(p2, ca, p2.tallest, p2.shortest, 1, 1,
                          Ref{Union{CLM.fates_cohort_type,Nothing}}(nothing),
                          Ref{Union{CLM.fates_cohort_type,Nothing}}(nothing))
        sb = Ref{Union{CLM.fates_cohort_type,Nothing}}(p2.tallest)
        ss = Ref{Union{CLM.fates_cohort_type,Nothing}}(p2.shortest)
        CLM.insert_cohort(p2, cc, p2.tallest, p2.shortest, 0, 0, sb, ss)
        p2.tallest = sb[]; p2.shortest = ss[]
        # now insert a middle cohort
        cm = CLM.fates_cohort_type(height = 20.0, n = 0.1, dbh = 20.0, pft = ipft)
        sb = Ref{Union{CLM.fates_cohort_type,Nothing}}(p2.tallest)
        ss = Ref{Union{CLM.fates_cohort_type,Nothing}}(p2.shortest)
        CLM.insert_cohort(p2, cm, p2.tallest, p2.shortest, 0, 0, sb, ss)
        p2.tallest = sb[]; p2.shortest = ss[]

        @test p2.tallest === ca
        @test p2.tallest.shorter === cm
        @test p2.shortest === cc
        @test cm.taller === ca
        @test cm.shorter === cc
        CLM.count_cohorts(p2)
        @test p2.countcohorts == 3

        # ------------------------------------------------------------------
        # Copy / ZeroValues on a list node (the create/copy/zero logic reuse).
        # ------------------------------------------------------------------
        src = patch.tallest
        dst = CLM.fates_cohort_type()
        CLM.Init(dst, _seed_prt(ipft, src.dbh))
        CLM.Copy(src, dst)
        @test dst.dbh == src.dbh
        @test dst.pft == src.pft
        @test dst.taller === nothing && dst.shorter === nothing
        @test dst.indexnumber == CLM.fates_unset_int

        CLM.ZeroValues(dst)
        @test dst.size_class == 1
        @test dst.daily_n_demand == -9.0
        @test all(dst.year_net_uptake .== 999.0)

    finally
        CLM.prt_global[]                  = old_global
        CLM.num_elements[]                = old_numel
        empty!(CLM.element_list)
        append!(CLM.element_list, old_ellist)
        CLM.hlm_parteh_mode[]             = old_parteh
        CLM.hlm_use_sp[]                  = old_use_sp
        CLM.hlm_use_planthydro[]          = old_hydro
        CLM.hlm_use_cohort_age_tracking[] = old_agetrk
        CLM.nleafage[]                    = old_nleafage
        CLM.nlevsclass[]                  = old_nlevsc
        CLM.nlevcoage[]                   = old_nlevca
        CLM.SFParams[]                    = old_sfp
    end
end
