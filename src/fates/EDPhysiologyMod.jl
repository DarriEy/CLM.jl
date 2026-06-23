# EDPhysiologyMod.jl
# Julia port of FATES src/fates/biogeochem/EDPhysiologyMod.F90
#
# The FATES physiology hub. Contains, faithfully ported:
#   * ZeroLitterFluxes / ZeroAllocationRates        — daily flux/rate resets
#   * phenology / phenology_leafonoff               — cold + drought deciduous
#                                                     leaf on/off economics
#   * trim_canopy                                   — LAI trimming to carbon
#                                                     balance (LLSF optimum LAI)
#   * recruitment                                   — spawn juvenile cohorts
#   * SeedUpdate / SeedDecay / SeedGermination      — seed-bank bookkeeping
#   * PreDisturbanceLitterFluxes /
#     PreDisturbanceIntegrateLitter / CWDInput /
#     CWDOut / fragmentation_scaler                 — litter / CWD fluxes
#   * GenerateDamageAndLitterFluxes                 — tree-damage litter (gated)
#   * satellite_phenology / calculate_SP_properties /
#     assign_cohort_SP_properties                   — SP (prescribed-LAI) mode
#   * UpdateRecruitL2FR / UpdateRecruitStoich /
#     SetRecruitL2FR                                — CNP recruit l2fr/stoich
#
# Deps (already ported):
#   FatesConstantsMod  : nearzero, sec_per_day, default_regeneration,
#                        TRS_regeneration, TRS_no_seedling_dyn,
#                        min_max_dbh_for_trees, megajoules_per_joule,
#                        mpa_per_mm_suction, g_per_kg, ndays_per_year,
#                        nocomp_bareground, nocomp_bareground_land, is_crop,
#                        area_error_2, area_error_3, itrue, ifalse, years_per_day,
#                        leaves_on, leaves_off, leaves_shedding,
#                        ihard_stress_decid, isemi_stress_decid,
#                        t_water_freeze_k_1atm, pi_const
#   EDTypesMod         : ed_site_type, area, area_inv, numWaterMem,
#                        num_vegtemp_mem, min_n_safemath, init_recruit_trim,
#                        homogenize_seed_pfts, phen_cstat_*, phen_dstat_*
#   FatesPatchMod      : fates_patch_type   FatesCohortMod : fates_cohort_type
#                        (ZeroValues, Copy, InitPRTBoundaryConditions)
#   FatesLitterMod     : litter_type, ncwd, ndcmpy, ZeroFlux!, adjust_SF_CWD_frac
#   EDCohortDynamicsMod: create_cohort, InitPRTObject!
#   FatesAllometryMod  : tree_lai, tree_sai, leafc_from_treelai, decay_coeff_vcmax,
#                        h2d_allom, bagw_allom, bsap_allom, bleaf, bfineroot,
#                        bdead_allom, bstore_allom, bbgw_allom, carea_allom,
#                        set_root_fraction
#   PRTGenericMod      : num_elements, element_list, element_pos, carbon12_element,
#                        nitrogen_element, phosphorus_element, leaf_organ,
#                        fnrt_organ, sapw_organ, store_organ, repro_organ,
#                        struct_organ, GetState, GetTurnover, ZeroRates!, SetState!,
#                        CheckInitialConditions, StorageNutrientTarget,
#                        prt_carbon_allom_hyp, prt_cnp_flex_allom_hyp
#   PRTLossFluxesMod   : PRTPhenologyFlush!, PRTDeciduousTurnover!, PRTReproRelease!,
#                        PRTDamageLosses!
#   PRTParametersMod   : prt_params
#   EDPftvarcon        : EDPftvarcon_inst (Ref), GetDecompyFrac
#   EDParamsMod        : ed_params() — fields regeneration_model,
#                        sdlng_mort_par_timescale, logging_export_frac, q10_mr,
#                        q10_froz, crop_lu_pft_vector, dinc_vai, dlower_vai;
#                        const nclmax, nlevleaf, maxpft
#   FatesParameterDerivedMod : param_derived() (branch_frac)
#   DamageMainMod      : damage_time, GetCrownReduction, GetDamageFrac
#   SFParamsMod        : sf_params() (SF_val_CWD_frac, SF_val_max_decomp)
#   FatesFuelClassesMod: fuel_classes, dead_leaves
#   FatesInterfaceTypesMod : numpft, nleafage, nlevdamage, hlm_* flags + dates,
#                        hlm_parteh_mode, hlm_freq_day, hlm_seeddisp_cadence,
#                        fates_dispersal_cadence_none
#   PRTParamsFATESMod  : NewRecruitTotalStoichiometry
#   EDBtranMod         : check_layer_water
#   FatesRunningMeanMod: GetMean
#   FatesGlobals       : fates_endrun
#
# NOT-yet-ported dependencies (stubbed, reached only on gated paths off by
# default in the test path):
#   FatesPlantHydraulicsMod : InitHydrCohort  (plant-hydro damage path;
#                             hlm_use_planthydro == ifalse by default)
#
# Standalone — NOT added to CLMInstances or any dual-copied struct.

using LinearAlgebra  # for the LLSF optimum-LAI least-squares fit in trim_canopy

# ---------------------------------------------------------------------------
# Module-level parameters (Fortran module-save parameters)
# ---------------------------------------------------------------------------

const dleafon_drycheck = 100         # Drought deciduous leaves max days on check
const decid_leaf_long_max = 1.0      # Max leaf lifespan for deciduous PFTs [yr]
const min_daysoff_dforcedflush = 30  # Days since drop required to force flush
const dd_offon_toler = 30            # Tolerance (days) around last-year's dates
const elongf_min = 0.05              # Minimum elongation factor (full abscission)
const smp_lwr_bound = -1.0e6         # Imposed SMP lower bound for frozen/dry soil

# ===========================================================================
# Plant-hydraulics cohort init on the tree-damage path (gated by
# hlm_use_planthydro). Allocates the new damaged cohort's hydro object.
# ===========================================================================
function _ed_phys_init_hydr_cohort!(currentSite, ndcohort)
    InitHydrCohort(currentSite, ndcohort)
    return nothing
end

# ===========================================================================
# ZeroLitterFluxes
# ===========================================================================
"""
    ZeroLitterFluxes(currentSite)

Loop through every patch in a site and zero the per-patch litter flux terms,
for every litter element. Typically called at the start of the dynamics
sequence. Mirrors the Fortran `ZeroLitterFluxes`.
"""
function ZeroLitterFluxes(currentSite::ed_site_type)
    currentPatch = currentSite.youngest_patch
    while currentPatch !== nothing
        for el in 1:num_elements[]
            ZeroFlux!(currentPatch.litter[el])
        end
        currentPatch = currentPatch.older
    end
    return nothing
end

# ===========================================================================
# ZeroAllocationRates
# ===========================================================================
"""
    ZeroAllocationRates(currentSite)

Zero the turnover/growth (allocation) rates on every cohort's PARTEH object.
Mirrors the Fortran `ZeroAllocationRates`.
"""
function ZeroAllocationRates(currentSite::ed_site_type)
    currentPatch = currentSite.youngest_patch
    while currentPatch !== nothing
        currentCohort = currentPatch.tallest
        while currentCohort !== nothing
            ZeroRates!(currentCohort.prt)
            currentCohort = currentCohort.shorter
        end
        currentPatch = currentPatch.older
    end
    return nothing
end

# ===========================================================================
# GenerateDamageAndLitterFluxes  (gated by hlm_use_tree_damage + damage_time)
# ===========================================================================
"""
    GenerateDamageAndLitterFluxes(csite, cpatch, bc_in)

Spawn damaged-class daughter cohorts and route the lost crown/branch biomass to
the litter input fluxes. Only active when `hlm_use_tree_damage == itrue` and it
is a damage time step (both off by default). Mirrors the Fortran routine.
"""
function GenerateDamageAndLitterFluxes(csite::ed_site_type, cpatch::fates_patch_type, bc_in)
    if hlm_use_tree_damage[] != itrue
        return nothing
    end
    if !damage_time[]
        return nothing
    end

    SF_val_CWD_frac = sf_params().SF_val_CWD_frac
    SF_val_CWD_frac_adj = zeros(ncwd)

    ccohort = cpatch.tallest
    while ccohort !== nothing
        # Ignore damage to new plants and non-woody plants
        if prt_params.woody[ccohort.pft] == ifalse || ccohort.isnew
            ccohort = ccohort.shorter
            continue
        end

        ipft        = ccohort.pft
        agb_frac    = prt_params.allom_agb_frac[ccohort.pft]
        branch_frac = param_derived().branch_frac[ccohort.pft]

        for cd in (ccohort.crowndamage + 1):nlevdamage[]
            cd_frac = GetDamageFrac(ccohort.crowndamage, cd, ipft)
            num_trees_cd = ccohort.n * cd_frac

            if num_trees_cd > nearzero
                # Create a new damaged cohort
                ndcohort = fates_cohort_type()
                if hlm_use_planthydro[] == itrue
                    _ed_phys_init_hydr_cohort!(csite, ndcohort)
                end

                ndcohort.prt = InitPRTObject!()
                InitPRTBoundaryConditions(ndcohort)
                ZeroValues(ndcohort)

                Copy(ccohort, ndcohort)

                ndcohort.n = num_trees_cd
                ndcohort.crowndamage = cd

                # Remove these trees from the donor cohort
                ccohort.n = ccohort.n - num_trees_cd

                _, ndcohort.c_area = carea_allom(ndcohort.dbh, ndcohort.n, csite.spread,
                                                 ipft, ndcohort.crowndamage)

                crown_loss_frac = GetCrownReduction(cd - ccohort.crowndamage)

                for el in 1:num_elements[]
                    litt = cpatch.litter[el]
                    elflux_diags = csite.flux_diags.elem[el]

                    branch_loss_frac = crown_loss_frac * branch_frac * agb_frac

                    leaf_loss   = GetState(ndcohort.prt, leaf_organ, element_list[el]) * crown_loss_frac
                    repro_loss  = GetState(ndcohort.prt, repro_organ, element_list[el]) * crown_loss_frac
                    sapw_loss   = GetState(ndcohort.prt, sapw_organ, element_list[el]) * branch_loss_frac
                    store_loss  = GetState(ndcohort.prt, store_organ, element_list[el]) * branch_loss_frac
                    struct_loss = GetState(ndcohort.prt, struct_organ, element_list[el]) * branch_loss_frac

                    for dcmpy in 1:ndcmpy
                        dcmpy_frac = GetDecompyFrac(ipft, leaf_organ, dcmpy)
                        litt.leaf_fines_in[dcmpy] += (store_loss + leaf_loss + repro_loss) *
                            ndcohort.n * dcmpy_frac / cpatch.area
                    end

                    elflux_diags.surf_fine_litter_input[ipft] +=
                        (store_loss + leaf_loss + repro_loss) * ndcohort.n

                    adjust_SF_CWD_frac(ndcohort.dbh, ncwd, SF_val_CWD_frac, SF_val_CWD_frac_adj)

                    for c in 1:ncwd
                        litt.ag_cwd_in[c] += (sapw_loss + struct_loss) *
                            SF_val_CWD_frac_adj[c] * ndcohort.n / cpatch.area
                        elflux_diags.cwd_ag_input[c] += (struct_loss + sapw_loss) *
                            SF_val_CWD_frac_adj[c] * ndcohort.n
                    end
                end

                # Apply the damage to the cohort PARTEH state.
                PRTDamageLosses!(ndcohort.prt, leaf_organ, crown_loss_frac)
                PRTDamageLosses!(ndcohort.prt, repro_organ, crown_loss_frac)
                PRTDamageLosses!(ndcohort.prt, sapw_organ, branch_loss_frac)
                PRTDamageLosses!(ndcohort.prt, store_organ, branch_loss_frac)
                PRTDamageLosses!(ndcohort.prt, struct_organ, branch_loss_frac)

                # Insert the new cohort into a taller position so the loop does not
                # revisit it.
                ndcohort.shorter = ccohort
                if ccohort.taller !== nothing
                    ndcohort.taller = ccohort.taller
                    ccohort.taller.shorter = ndcohort
                else
                    cpatch.tallest = ndcohort
                    ndcohort.taller = nothing
                end
                ccohort.taller = ndcohort
            end
        end

        ccohort = ccohort.shorter
    end
    return nothing
end

# ===========================================================================
# PreDisturbanceLitterFluxes
# ===========================================================================
"""
    PreDisturbanceLitterFluxes(currentSite, currentPatch, bc_in)

Compute all non-disturbance litter input/output fluxes for a patch: seed decay,
seed germination, CWD input from turnover/non-disturbance mortality, and
fragmentation. Mirrors the Fortran `PreDisturbanceLitterFluxes`.
"""
function PreDisturbanceLitterFluxes(currentSite::ed_site_type, currentPatch::fates_patch_type, bc_in)
    # Calculate the fragmentation rates
    fragmentation_scaler(currentPatch, bc_in)

    for el in 1:num_elements[]
        litt      = currentPatch.litter[el]
        site_mass = currentSite.mass_balance[el]
        diag      = currentSite.flux_diags.elem[el]

        # Loss rate of viable seeds to litter
        SeedDecay(litt, currentPatch, bc_in)

        # Seed germination rate (status flags suppress germination in drought/cold)
        SeedGermination(litt, currentSite.cstatus,
                        view(currentSite.dstatus, 1:numpft[]), bc_in, currentPatch)

        # Fluxes from newly created litter into the litter pools (turnover + non-
        # disturbance mortality + live-tree litterfall)
        CWDInput(currentSite, currentPatch, litt, bc_in)

        # Fragmentation over the active soil layers
        nlev_eff_decomp = max(bc_in.max_rooting_depth_index_col, 1)
        CWDOut(litt, currentPatch.fragmentation_scaler, nlev_eff_decomp)

        # Fragmentation flux to soil decomposition model [kg/site/day]
        site_mass.frag_out += currentPatch.area *
            (sum(litt.ag_cwd_frag) + sum(litt.bg_cwd_frag) +
             sum(litt.leaf_fines_frag) + sum(litt.root_fines_frag) +
             sum(litt.seed_decay) + sum(litt.seed_germ_decay))

        # Total seed decay diagnostic [kg/m2/day]
        diag.tot_seed_turnover += (sum(litt.seed_decay) + sum(litt.seed_germ_decay)) *
            currentPatch.area * area_inv
    end
    return nothing
end

# ===========================================================================
# PreDisturbanceIntegrateLitter
# ===========================================================================
"""
    PreDisturbanceIntegrateLitter(currentPatch)

Apply the (already-computed) litter fluxes to the prognostic litter state for a
patch: seed bank, germinated seed pool, AG/BG CWD, and fine litter. The
integration step is one day (time implied). Mirrors the Fortran routine.
"""
function PreDisturbanceIntegrateLitter(currentPatch::fates_patch_type)
    for el in 1:num_elements[]
        litt = currentPatch.litter[el]

        # Bank of viable seeds + germinated seed pool
        for pft in 1:numpft[]
            litt.seed[pft] += litt.seed_in_local[pft] + litt.seed_in_extern[pft] -
                              litt.seed_decay[pft] - litt.seed_germ_in[pft]

            litt.seed_germ[pft] += litt.seed_germ_in[pft] - litt.seed_germ_decay[pft]
        end

        # Coarse Woody Debris (above and below)
        nlevsoil = size(litt.bg_cwd, 2)
        for c in 1:ncwd
            litt.ag_cwd[c] += litt.ag_cwd_in[c] - litt.ag_cwd_frag[c]
            for ilyr in 1:nlevsoil
                litt.bg_cwd[c, ilyr] += litt.bg_cwd_in[c, ilyr] - litt.bg_cwd_frag[c, ilyr]
            end
        end

        # Fine litter pools (leaves and fine roots)
        for dcmpy in 1:ndcmpy
            litt.leaf_fines[dcmpy] += litt.leaf_fines_in[dcmpy] - litt.leaf_fines_frag[dcmpy]
            for ilyr in 1:nlevsoil
                litt.root_fines[dcmpy, ilyr] += litt.root_fines_in[dcmpy, ilyr] -
                                                litt.root_fines_frag[dcmpy, ilyr]
            end
        end
    end
    return nothing
end

# ===========================================================================
# trim_canopy
# ===========================================================================
"""
    trim_canopy(currentSite)

Canopy trimming / leaf optimisation. Removes leaves from leaf layers in negative
annual carbon balance, and uses a linear least-squares fit of net-net uptake vs
cumulative LAI to estimate the optimum trim. Mirrors the Fortran `trim_canopy`.

The Fortran solves the 2x2 normal-equations least-squares problem with LAPACK
`dgels`; here we solve the same `A x = b` system with Julia's `\\` (least-squares
via QR), which is numerically equivalent for this well-conditioned 2x2 system.
"""
function trim_canopy(currentSite::ed_site_type)
    nll = 3  # Number of leaf layers to fit a regression to

    currentPatch = currentSite.youngest_patch
    while currentPatch !== nothing
        currentCohort = currentPatch.tallest
        while currentCohort !== nothing
            initial_trim = currentCohort.canopy_trim
            trimmed = false
            ipft = currentCohort.pft

            _, currentCohort.c_area = carea_allom(currentCohort.dbh, currentCohort.n,
                currentSite.spread, currentCohort.pft, currentCohort.crowndamage)

            leaf_c = GetState(currentCohort.prt, leaf_organ, carbon12_element)

            currentCohort.treelai = tree_lai(leaf_c, currentCohort.pft, currentCohort.c_area,
                currentCohort.n, currentCohort.canopy_layer,
                currentPatch.canopy_layer_tlai, currentCohort.vcmax25top)

            currentCohort.treesai = tree_sai(currentCohort.pft, currentCohort.dbh,
                currentCohort.crowndamage, currentCohort.canopy_trim,
                currentCohort.efstem_coh, currentCohort.c_area, currentCohort.n,
                currentCohort.canopy_layer, currentPatch.canopy_layer_tlai,
                currentCohort.treelai, currentCohort.vcmax25top, 0)

            dlower_vai = ed_params().dlower_vai
            dinc_vai   = ed_params().dinc_vai

            currentCohort.nv = count(<(currentCohort.treelai + currentCohort.treesai),
                                     view(dlower_vai, 1:nlevleaf)) + 1

            if currentCohort.nv > nlevleaf
                fates_endrun("nv > nlevleaf in trim_canopy: nv=$(currentCohort.nv)")
            end

            # Target leaf biomass (fully flushed, elongation factor = 1)
            tar_bl, _ = bleaf(currentCohort.dbh, ipft, currentCohort.crowndamage,
                              currentCohort.canopy_trim, 1.0)

            bfr_per_bleaf = 0.0
            if Int(prt_params.allom_fmode[ipft]) == 1
                tar_bfr, _ = bfineroot(currentCohort.dbh, ipft, currentCohort.canopy_trim,
                                       currentCohort.l2fr, 1.0)
                bfr_per_bleaf = tar_bfr / tar_bl
            end

            cl = currentCohort.canopy_layer

            # Leaf lifespan depends on canopy layer
            if cl == 1
                leaf_long = sum(view(prt_params.leaf_long, ipft, :))
            else
                leaf_long = sum(view(prt_params.leaf_long_ustory, ipft, :))
            end

            sla_max = prt_params.slamax[ipft]

            # 2x2 normal-equations system for the LLSF fit
            nnu_clai_a = zeros(2, 2)
            nnu_clai_b = zeros(2, 1)

            cumulative_lai_cohort = 0.0  # carried out of the z loop, as in Fortran

            for z in 1:currentCohort.nv
                leaf_inc = dinc_vai[z] *
                    currentCohort.treelai / (currentCohort.treelai + currentCohort.treesai)

                lai_layers_above = (dlower_vai[z] - dinc_vai[z]) *
                    currentCohort.treelai / (currentCohort.treelai + currentCohort.treesai)
                lai_current = min(leaf_inc, currentCohort.treelai - lai_layers_above)
                cumulative_lai_cohort = lai_layers_above + 0.5 * lai_current

                lai_canopy_above = sum(view(currentPatch.canopy_layer_tlai, 1:(cl - 1)))
                cumulative_lai = lai_canopy_above + cumulative_lai_cohort

                if currentCohort.year_net_uptake[z] != 999.0
                    kn = decay_coeff_vcmax(currentCohort.vcmax25top,
                        prt_params.leafn_vert_scaler_coeff1[ipft],
                        prt_params.leafn_vert_scaler_coeff2[ipft])

                    nscaler_levleaf = exp(-kn * cumulative_lai)
                    sla_levleaf = min(sla_max, prt_params.slatop[ipft] / nscaler_levleaf)

                    # Realised leaf lifespan, depending on phenology
                    if prt_params.season_decid[ipft] == itrue
                        pft_leaf_lifespan = decid_leaf_long_max
                    elseif prt_params.stress_decid[ipft] == ihard_stress_decid ||
                           prt_params.stress_decid[ipft] == isemi_stress_decid
                        pft_leaf_lifespan = min(decid_leaf_long_max, leaf_long)
                    else
                        pft_leaf_lifespan = leaf_long
                    end

                    currentCohort.leaf_cost = 1.0 / (sla_levleaf * pft_leaf_lifespan * g_per_kg)

                    if Int(prt_params.allom_fmode[ipft]) == 1
                        currentCohort.leaf_cost += 1.0 / (sla_levleaf * g_per_kg) *
                            bfr_per_bleaf / prt_params.root_long[ipft]
                    end
                    currentCohort.leaf_cost *= (prt_params.grperc[ipft] + 1.0)

                    # Build the LLSF arrays for the bottom nll leaf layers
                    if currentCohort.nv > nll && currentCohort.nv - z < nll
                        x = currentCohort.year_net_uptake[z] - currentCohort.leaf_cost
                        nnu_clai_a[1, 1] += 1
                        nnu_clai_a[1, 2] += x
                        nnu_clai_a[2, 1] = nnu_clai_a[1, 2]
                        nnu_clai_a[2, 2] += x^2
                        nnu_clai_b[1, 1] += cumulative_lai_cohort
                        nnu_clai_b[2, 1] += cumulative_lai_cohort * x
                    end

                    # Check leaf cost against the yearly net uptake for this layer
                    if currentCohort.year_net_uptake[z] < currentCohort.leaf_cost
                        if currentCohort.canopy_trim > EDPftvarcon_inst[].trim_limit[ipft]
                            if currentCohort.height > EDPftvarcon_inst[].hgt_min[ipft]
                                currentCohort.canopy_trim -= EDPftvarcon_inst[].trim_inc[ipft]
                                trimmed = true
                            end
                        end
                    end
                end
            end

            # Compute the optimal cumulative LAI if at least 2 leaf layers
            if nnu_clai_a[1, 1] > 1
                # Solve A x = b in the least-squares sense (== LAPACK dgels).
                # x = [b_intercept; m_slope]
                x = nnu_clai_a \ nnu_clai_b

                if cumulative_lai_cohort > 0.0
                    optimum_trim = (x[1, 1] / cumulative_lai_cohort) * initial_trim
                    if optimum_trim > 0.0 && optimum_trim < 1.0
                        currentCohort.canopy_trim = optimum_trim
                        trimmed = true
                    end
                end
            end

            # Reset activity for the start of the next year
            fill!(currentCohort.year_net_uptake, 999.0)

            # Add to trim fraction if cohort not trimmed at all
            if !trimmed && currentCohort.canopy_trim < 1.0
                currentCohort.canopy_trim += EDPftvarcon_inst[].trim_inc[ipft]
            end

            currentCohort = currentCohort.shorter
        end
        currentPatch = currentPatch.older
    end
    return nothing
end

# ===========================================================================
# phenology  (cold + drought deciduous leaf on/off decision)
# ===========================================================================
"""
    phenology(currentSite, bc_in)

Daily phenology driver. Advances the integer model day, accumulates growing/
chilling degree days, updates the cold-deciduous site status (`cstatus`), tracks
soil-moisture memory, and assigns per-PFT drought-deciduous status (`dstatus`)
and elongation factors. Finally applies the leaf on/off economics via
`phenology_leafonoff`. Mirrors the Fortran `phenology`.
"""
function phenology(currentSite::ed_site_type, bc_in)
    tfrz = t_water_freeze_k_1atm
    p = ed_params()

    # Advance the monotone integer model day
    currentSite.phen_model_date += 1
    model_day_int = currentSite.phen_model_date

    # Patch-area-weighted mean 24h veg temperature, converted to Celsius
    temp_in_C = 0.0
    cpatch = currentSite.oldest_patch
    while cpatch !== nothing
        temp_in_C += GetMean(cpatch.tveg24) * cpatch.area
        cpatch = cpatch.younger
    end
    temp_in_C = temp_in_C * area_inv - tfrz

    # ----------------- Cold Phenology -----------------
    if currentSite.lat > 0
        ncdstart = 270  # NH November
        gddstart = 1    # NH January
    else
        ncdstart = 120  # SH May
        gddstart = 181  # SH July
    end

    if hlm_day_of_year[] == ncdstart
        currentSite.nchilldays = 0
    end

    if temp_in_C < p.ED_val_phen_chiltemp
        currentSite.nchilldays += 1
    end

    # GDD exceedance threshold (depends on chilling days)
    gdd_threshold = p.ED_val_phen_a +
        p.ED_val_phen_b * exp(p.ED_val_phen_c * Float64(currentSite.nchilldays))

    # Roll the last-10-days temperature memory
    currentSite.vegtemp_memory[2:num_vegtemp_mem] = currentSite.vegtemp_memory[1:(num_vegtemp_mem - 1)]
    currentSite.vegtemp_memory[1] = temp_in_C

    # Count cold days
    ncolddays = 0
    for i_tmem in 1:num_vegtemp_mem
        if currentSite.vegtemp_memory[i_tmem] < p.ED_val_phen_coldtemp
            ncolddays += 1
        end
    end

    # Reset GDD on the set date
    if hlm_day_of_year[] == gddstart
        currentSite.grow_deg_days = 0.0
    end

    # Accumulate GDD with daily means (only while in the cold/leaves-off state)
    if temp_in_C > 0.0 && currentSite.cstatus == phen_cstat_iscold
        currentSite.grow_deg_days += temp_in_C
    end

    # Prevent GDD accumulating after leaf fall and before the accumulation period
    if model_day_int > ndays_per_year
        if currentSite.lat > 0.0  # Northern Hemisphere
            if model_day_int > currentSite.cleafoffdate && hlm_day_of_year[] > 180
                currentSite.grow_deg_days = 0.0
            end
        else  # Southern Hemisphere
            if model_day_int > currentSite.cleafoffdate && hlm_day_of_year[] < gddstart
                currentSite.grow_deg_days = 0.0
            end
        end
    end

    # Days since leaves last came on/off (provision for first year)
    if model_day_int < currentSite.cleafoffdate
        currentSite.cndaysleafoff = model_day_int - (currentSite.cleafoffdate - ndays_per_year)
    else
        currentSite.cndaysleafoff = model_day_int - currentSite.cleafoffdate
    end
    if model_day_int < currentSite.cleafondate
        currentSite.cndaysleafon = model_day_int - (currentSite.cleafondate - ndays_per_year)
    else
        currentSite.cndaysleafon = model_day_int - currentSite.cleafondate
    end

    # LEAF ON: COLD DECIDUOUS
    if (currentSite.cstatus == phen_cstat_iscold || currentSite.cstatus == phen_cstat_nevercold) &&
       (currentSite.grow_deg_days > gdd_threshold) &&
       (currentSite.cndaysleafoff > p.ED_val_phen_mindayson) &&
       (currentSite.nchilldays >= 1)
        currentSite.cstatus = phen_cstat_notcold
        currentSite.cleafondate = model_day_int
        currentSite.cndaysleafon = 0
        currentSite.grow_deg_days = 0.0
    end

    # LEAF OFF: COLD THRESHOLD
    if (currentSite.cstatus == phen_cstat_notcold) &&
       (model_day_int > num_vegtemp_mem) &&
       (ncolddays > p.ED_val_phen_ncolddayslim) &&
       (currentSite.cndaysleafon > p.ED_val_phen_mindayson)
        currentSite.grow_deg_days = 0.0
        currentSite.cstatus = phen_cstat_iscold
        currentSite.cleafoffdate = model_day_int
        currentSite.cndaysleafoff = 0
    end

    # LEAF OFF: COLD LIFESPAN THRESHOLD (no-cold-day regions)
    if (currentSite.cstatus == phen_cstat_notcold) && (currentSite.cndaysleafoff > 400)
        currentSite.grow_deg_days = 0.0
        currentSite.cstatus = phen_cstat_nevercold
        currentSite.cleafoffdate = model_day_int
        currentSite.cndaysleafoff = 0
    end

    # ----------------- Drought Phenology (per-PFT elongation factors) ---------
    for ipft in 1:numpft[]
        phen_drought_threshold = prt_params.phen_drought_threshold[ipft]
        phen_moist_threshold   = prt_params.phen_moist_threshold[ipft]
        phen_doff_time         = prt_params.phen_doff_time[ipft]

        # Shift soil moisture memory to previous day
        for i_wmem in numWaterMem:-1:2
            currentSite.liqvol_memory[i_wmem, ipft] = currentSite.liqvol_memory[i_wmem - 1, ipft]
            currentSite.smp_memory[i_wmem, ipft]    = currentSite.smp_memory[i_wmem - 1, ipft]
        end

        # Rooting depth distribution for PFT
        set_root_fraction(currentSite.rootfrac_scr, ipft, currentSite.zi_soil;
                          max_nlevroot=bc_in.max_rooting_depth_index_col)
        nlevroot = max(2, min(length(currentSite.zi_soil) - 1, bc_in.max_rooting_depth_index_col))

        # Ignore the (thin, fast-drying) top layer when forming the memory
        rootfrac_notop = sum(view(currentSite.rootfrac_scr, 2:nlevroot))
        if rootfrac_notop <= nearzero
            currentSite.rootfrac_scr[2] = 1.0
            rootfrac_notop = 1.0
        end

        currentSite.liqvol_memory[1, ipft] =
            sum(view(bc_in.h2o_liqvol_sl, 2:nlevroot) .* view(currentSite.rootfrac_scr, 2:nlevroot)) /
            rootfrac_notop
        currentSite.smp_memory[1, ipft] = 0.0
        for j in 2:nlevroot
            if check_layer_water(bc_in.h2o_liqvol_sl[j], bc_in.tempk_sl[j])
                currentSite.smp_memory[1, ipft] +=
                    bc_in.smp_sl[j] * currentSite.rootfrac_scr[j] / rootfrac_notop
            else
                currentSite.smp_memory[1, ipft] +=
                    smp_lwr_bound * currentSite.rootfrac_scr[j] / rootfrac_notop
            end
        end

        mean_10day_liqvol = sum(view(currentSite.liqvol_memory, 1:numWaterMem, ipft)) / Float64(numWaterMem)
        mean_10day_smp    = sum(view(currentSite.smp_memory, 1:numWaterMem, ipft)) / Float64(numWaterMem)

        if phen_drought_threshold >= 0.0
            smoist_below_threshold = mean_10day_liqvol < phen_drought_threshold
        else
            smoist_below_threshold = mean_10day_smp < phen_drought_threshold
        end

        # Days since last flushing/shedding (provision for first year)
        if model_day_int < currentSite.dleafoffdate[ipft]
            currentSite.dndaysleafoff[ipft] = model_day_int - (currentSite.dleafoffdate[ipft] - ndays_per_year)
        else
            currentSite.dndaysleafoff[ipft] = model_day_int - currentSite.dleafoffdate[ipft]
        end
        if model_day_int < currentSite.dleafondate[ipft]
            currentSite.dndaysleafon[ipft] = model_day_int - (currentSite.dleafondate[ipft] - ndays_per_year)
        else
            currentSite.dndaysleafon[ipft] = model_day_int - currentSite.dleafondate[ipft]
        end

        elongf_prev = currentSite.elong_factor[ipft]

        ndays_pft_leaf_lifespan = round(Int,
            ndays_per_year * min(decid_leaf_long_max, sum(view(prt_params.leaf_long, ipft, :))))

        sd = prt_params.stress_decid[ipft]
        if sd == ihard_stress_decid
            # "Hard" drought deciduous
            exceed_min_on_period =
                (currentSite.dstatus[ipft] == phen_dstat_timeon || currentSite.dstatus[ipft] == phen_dstat_moiston) &&
                (currentSite.dndaysleafon[ipft] > dleafon_drycheck)
            exceed_min_off_period =
                (currentSite.dstatus[ipft] == phen_dstat_timeoff) &&
                (currentSite.dndaysleafoff[ipft] > min_daysoff_dforcedflush)
            prolonged_on_period =
                (currentSite.dstatus[ipft] == phen_dstat_timeon || currentSite.dstatus[ipft] == phen_dstat_moiston) &&
                (currentSite.dndaysleafon[ipft] > ndays_pft_leaf_lifespan)
            prolonged_off_period =
                (currentSite.dstatus[ipft] == phen_dstat_timeoff || currentSite.dstatus[ipft] == phen_dstat_moistoff) &&
                (currentSite.dndaysleafoff[ipft] > phen_doff_time) &&
                (currentSite.dndaysleafon[ipft] >= ndays_per_year - dd_offon_toler) &&
                (currentSite.dndaysleafon[ipft] <= ndays_per_year + dd_offon_toler)
            last_flush_long_ago =
                (currentSite.dstatus[ipft] == phen_dstat_moistoff) &&
                (currentSite.dndaysleafon[ipft] > ndays_per_year + dd_offon_toler)

            if model_day_int > numWaterMem
                if prolonged_off_period && !smoist_below_threshold
                    currentSite.dstatus[ipft] = phen_dstat_moiston
                    currentSite.dleafondate[ipft] = model_day_int
                    currentSite.dndaysleafon[ipft] = 0
                    currentSite.elong_factor[ipft] = 1.0
                elseif last_flush_long_ago
                    currentSite.dstatus[ipft] = phen_dstat_timeon
                    currentSite.dleafondate[ipft] = model_day_int
                    currentSite.dndaysleafon[ipft] = 0
                    currentSite.elong_factor[ipft] = 1.0
                elseif exceed_min_off_period
                    currentSite.dstatus[ipft] = phen_dstat_timeon
                    currentSite.dleafondate[ipft] = model_day_int
                    currentSite.dndaysleafon[ipft] = 0
                    currentSite.elong_factor[ipft] = 1.0
                elseif prolonged_on_period
                    currentSite.dstatus[ipft] = phen_dstat_timeoff
                    currentSite.dleafoffdate[ipft] = model_day_int
                    currentSite.dndaysleafoff[ipft] = 0
                    currentSite.elong_factor[ipft] = 0.0
                elseif exceed_min_on_period && smoist_below_threshold
                    currentSite.dstatus[ipft] = phen_dstat_moistoff
                    currentSite.dleafoffdate[ipft] = model_day_int
                    currentSite.dndaysleafoff[ipft] = 0
                    currentSite.elong_factor[ipft] = 0.0
                end
            end

        elseif sd == isemi_stress_decid
            # Semi-deciduous PFT (ED2-style gradual elongation factor)
            if phen_drought_threshold >= 0.0
                elongf_1st = elongf_min + (1.0 - elongf_min) *
                    (mean_10day_liqvol - phen_drought_threshold) /
                    (phen_moist_threshold - phen_drought_threshold)
            else
                elongf_1st = elongf_min + (1.0 - elongf_min) *
                    (mean_10day_smp - phen_drought_threshold) /
                    (phen_moist_threshold - phen_drought_threshold)
            end
            elongf_1st = max(0.0, min(1.0, elongf_1st))

            recent_flush = elongf_prev >= elongf_min &&
                (currentSite.dndaysleafon[ipft] <= dleafon_drycheck)
            recent_abscission = elongf_prev < elongf_min &&
                (currentSite.dndaysleafoff[ipft] <= min_daysoff_dforcedflush)
            prolonged_on_period = (elongf_prev >= elongf_min && elongf_1st >= elongf_min) &&
                (currentSite.dndaysleafon[ipft] > ndays_pft_leaf_lifespan)
            last_flush_long_ago = (elongf_prev < elongf_min && elongf_1st < elongf_min) &&
                (currentSite.dndaysleafon[ipft] > ndays_per_year + dd_offon_toler)

            if model_day_int <= numWaterMem
                currentSite.elong_factor[ipft] = elongf_prev
            elseif prolonged_on_period
                currentSite.elong_factor[ipft] = 0.0
                currentSite.dstatus[ipft] = phen_dstat_timeoff
                currentSite.dleafoffdate[ipft] = model_day_int
                currentSite.dndaysleafoff[ipft] = 0
            elseif last_flush_long_ago
                currentSite.elong_factor[ipft] = elongf_min
                currentSite.dstatus[ipft] = phen_dstat_timeon
                currentSite.dleafondate[ipft] = model_day_int
                currentSite.dndaysleafon[ipft] = 0
            elseif recent_flush && elongf_1st < elongf_prev
                currentSite.elong_factor[ipft] = elongf_prev
                currentSite.dstatus[ipft] = phen_dstat_timeon
            elseif recent_abscission && elongf_1st > elongf_min
                currentSite.elong_factor[ipft] = 0.0
                currentSite.dstatus[ipft] = phen_dstat_timeoff
            elseif elongf_1st < elongf_min
                currentSite.elong_factor[ipft] = 0.0
                if elongf_prev >= elongf_min
                    currentSite.dstatus[ipft] = phen_dstat_moistoff
                    currentSite.dleafoffdate[ipft] = model_day_int
                    currentSite.dndaysleafoff[ipft] = 0
                end
            else
                currentSite.elong_factor[ipft] = elongf_1st
                if elongf_prev < elongf_min
                    currentSite.dstatus[ipft] = phen_dstat_moiston
                    currentSite.dleafondate[ipft] = model_day_int
                    currentSite.dndaysleafon[ipft] = 0
                elseif elongf_1st < elongf_prev
                    currentSite.dstatus[ipft] = phen_dstat_pshed
                end
            end

        else
            # Not drought deciduous: treat as moist-on; assign elongation from cold
            currentSite.dstatus[ipft] = phen_dstat_moiston
            if prt_params.season_decid[ipft] == ifalse
                currentSite.elong_factor[ipft] = 1.0   # Evergreen
            elseif prt_params.season_decid[ipft] == itrue
                # Cold deciduous: elongation factor from cold status
                if currentSite.cstatus == phen_cstat_nevercold || currentSite.cstatus == phen_cstat_iscold
                    currentSite.elong_factor[ipft] = 0.0
                elseif currentSite.cstatus == phen_cstat_notcold
                    currentSite.elong_factor[ipft] = 1.0
                end
            end
        end
    end

    phenology_leafonoff(currentSite)
    return nothing
end

# ===========================================================================
# phenology_leafonoff
# ===========================================================================
"""
    phenology_leafonoff(currentSite)

Controls the leaf on/off economics for every cohort: flush carbon from storage
into living tissues when it is time to flush leaves, and shed leaves/fine-roots
(and, for non-woody PFTs, sapwood/structure) when it is time to drop. Mirrors the
Fortran `phenology_leafonoff`.
"""
function phenology_leafonoff(currentSite::ed_site_type)
    carbon_store_buffer = 0.10

    currentPatch = currentSite.oldest_patch
    while currentPatch !== nothing
        currentCohort = currentPatch.tallest
        while currentCohort !== nothing
            ipft = currentCohort.pft

            store_c  = GetState(currentCohort.prt, store_organ, carbon12_element)
            leaf_c   = GetState(currentCohort.prt, leaf_organ, carbon12_element)
            fnrt_c   = GetState(currentCohort.prt, fnrt_organ, carbon12_element)
            sapw_c   = GetState(currentCohort.prt, sapw_organ, carbon12_element)
            struct_c = GetState(currentCohort.prt, struct_organ, carbon12_element)

            fnrt_drop_fraction = prt_params.phen_fnrt_drop_fraction[ipft]
            stem_drop_fraction = prt_params.phen_stem_drop_fraction[ipft]
            l2fr               = prt_params.allom_l2fr[ipft]

            if prt_params.season_decid[ipft] == itrue  # Cold deciduous
                is_flushing_time = (currentSite.cstatus == phen_cstat_notcold &&
                                    currentCohort.status_coh == leaves_off)
                is_shedding_time =
                    (currentSite.cstatus == phen_cstat_nevercold || currentSite.cstatus == phen_cstat_iscold) &&
                    currentCohort.status_coh == leaves_on &&
                    (currentCohort.dbh > EDPftvarcon_inst[].phen_cold_size_threshold[ipft] ||
                     prt_params.woody[ipft] == itrue)
            elseif prt_params.stress_decid[ipft] == ihard_stress_decid ||
                   prt_params.stress_decid[ipft] == isemi_stress_decid  # Drought deciduous
                is_flushing_time =
                    (currentSite.dstatus[ipft] == phen_dstat_moiston || currentSite.dstatus[ipft] == phen_dstat_timeon) &&
                    (currentCohort.status_coh == leaves_off || currentCohort.status_coh == leaves_shedding)
                is_shedding_time =
                    (currentSite.dstatus[ipft] == phen_dstat_moistoff || currentSite.dstatus[ipft] == phen_dstat_timeoff ||
                     currentSite.dstatus[ipft] == phen_dstat_pshed) &&
                    (currentCohort.status_coh == leaves_on || currentCohort.status_coh == leaves_shedding)
            else
                is_flushing_time = false
                is_shedding_time = false
            end

            currentCohort.efleaf_coh = currentSite.elong_factor[ipft]
            currentCohort.effnrt_coh = 1.0 - (1.0 - currentCohort.efleaf_coh) * fnrt_drop_fraction
            currentCohort.efstem_coh = 1.0 - (1.0 - currentCohort.efleaf_coh) * stem_drop_fraction

            # Target biomass for each tissue accounting for elongation factor
            target_leaf_c, _ = bleaf(currentCohort.dbh, currentCohort.pft,
                currentCohort.crowndamage, currentCohort.canopy_trim, currentCohort.efleaf_coh)
            target_fnrt_c, _ = bfineroot(currentCohort.dbh, currentCohort.pft,
                currentCohort.canopy_trim, l2fr, currentCohort.effnrt_coh)
            _, target_sapw_c, _ = bsap_allom(currentCohort.dbh, currentCohort.pft,
                currentCohort.crowndamage, currentCohort.canopy_trim, currentCohort.efstem_coh)
            target_agw_c, _ = bagw_allom(currentCohort.dbh, currentCohort.pft,
                currentCohort.crowndamage, currentCohort.efstem_coh)
            target_bgw_c, _ = bbgw_allom(currentCohort.dbh, currentCohort.pft, currentCohort.efstem_coh)
            target_struct_c, _ = bdead_allom(target_agw_c, target_bgw_c, target_sapw_c, currentCohort.pft)

            # A. Flush
            if is_flushing_time
                currentCohort.status_coh = leaves_on

                if store_c > nearzero
                    leaf_deficit_c   = max(0.0, target_leaf_c - leaf_c)
                    fnrt_deficit_c   = max(0.0, target_fnrt_c - fnrt_c)
                    sapw_deficit_c   = max(0.0, target_sapw_c - sapw_c)
                    struct_deficit_c = max(0.0, target_struct_c - struct_c)
                    total_deficit_c  = leaf_deficit_c + fnrt_deficit_c + sapw_deficit_c + struct_deficit_c

                    store_c_transfer_frac = min(
                        EDPftvarcon_inst[].phenflush_fraction[ipft] * total_deficit_c / store_c,
                        1.0 - carbon_store_buffer)

                    if total_deficit_c > nearzero
                        PRTPhenologyFlush!(currentCohort.prt, ipft, leaf_organ,
                            store_c_transfer_frac * leaf_deficit_c / total_deficit_c)
                        PRTPhenologyFlush!(currentCohort.prt, ipft, fnrt_organ,
                            store_c_transfer_frac * fnrt_deficit_c / total_deficit_c)

                        if prt_params.woody[ipft] == ifalse
                            PRTPhenologyFlush!(currentCohort.prt, ipft, sapw_organ,
                                store_c_transfer_frac * sapw_deficit_c / total_deficit_c)
                            PRTPhenologyFlush!(currentCohort.prt, ipft, struct_organ,
                                store_c_transfer_frac * struct_deficit_c / total_deficit_c)
                        end
                    end
                end
            end

            # B. Shed
            if is_shedding_time
                if currentCohort.efleaf_coh > 0.0
                    currentCohort.status_coh = leaves_shedding
                else
                    currentCohort.status_coh = leaves_off
                end

                eff_leaf_drop_fraction   = max(0.0, min(1.0, 1.0 - target_leaf_c / max(leaf_c, nearzero)))
                eff_fnrt_drop_fraction   = max(0.0, min(1.0, 1.0 - target_fnrt_c / max(fnrt_c, nearzero)))
                eff_sapw_drop_fraction   = max(0.0, min(1.0, 1.0 - target_sapw_c / max(sapw_c, nearzero)))
                eff_struct_drop_fraction = max(0.0, min(1.0, 1.0 - target_struct_c / max(struct_c, nearzero)))

                PRTDeciduousTurnover!(currentCohort.prt, ipft, leaf_organ, eff_leaf_drop_fraction)
                PRTDeciduousTurnover!(currentCohort.prt, ipft, fnrt_organ, eff_fnrt_drop_fraction)

                if prt_params.woody[ipft] == ifalse
                    PRTDeciduousTurnover!(currentCohort.prt, ipft, sapw_organ, eff_sapw_drop_fraction)
                    PRTDeciduousTurnover!(currentCohort.prt, ipft, struct_organ, eff_struct_drop_fraction)
                end
            end

            currentCohort = currentCohort.shorter
        end
        currentPatch = currentPatch.younger
    end
    return nothing
end

# ===========================================================================
# satellite_phenology  (SP / prescribed-LAI mode)
# ===========================================================================
"""
    satellite_phenology(currentSite, bc_in)

Translate the daily HLM-asserted LAI/SAI/canopy-height drivers into FATES cohort
structure (one cohort per PFT/patch), weighting the FATES-PFT targets by the
contributing HLM PFTs. Mirrors the Fortran `satellite_phenology`.
"""
function satellite_phenology(currentSite::ed_site_type, bc_in)
    fill!(currentSite.sp_tlai, 0.0)
    fill!(currentSite.sp_tsai, 0.0)
    fill!(currentSite.sp_htop, 0.0)

    hlm_pft_map = EDPftvarcon_inst[].hlm_pft_map

    currentPatch = currentSite.oldest_patch
    while currentPatch !== nothing
        fates_pft = currentPatch.nocomp_pft_label
        if fates_pft != 0
            for hlm_pft in 1:size(hlm_pft_map, 2)
                if bc_in.pft_areafrac[hlm_pft] * hlm_pft_map[fates_pft, hlm_pft] > 0.0
                    currentSite.sp_tlai[fates_pft] += bc_in.hlm_sp_tlai[hlm_pft] *
                        bc_in.pft_areafrac[hlm_pft] * hlm_pft_map[fates_pft, hlm_pft]
                    currentSite.sp_tsai[fates_pft] += bc_in.hlm_sp_tsai[hlm_pft] *
                        bc_in.pft_areafrac[hlm_pft] * hlm_pft_map[fates_pft, hlm_pft]
                    currentSite.sp_htop[fates_pft] += bc_in.hlm_sp_htop[hlm_pft] *
                        bc_in.pft_areafrac[hlm_pft] * hlm_pft_map[fates_pft, hlm_pft]
                end
            end

            if currentPatch.area > 0.0
                currentSite.sp_tlai[fates_pft] /= (currentPatch.area / area)
                currentSite.sp_tsai[fates_pft] /= (currentPatch.area / area)
                currentSite.sp_htop[fates_pft] /= (currentPatch.area / area)
            end
        end
        currentPatch = currentPatch.younger
    end

    currentPatch = currentSite.oldest_patch
    while currentPatch !== nothing
        currentCohort = currentPatch.tallest
        while currentCohort !== nothing
            fates_pft = currentCohort.pft
            if fates_pft != currentPatch.nocomp_pft_label
                fates_endrun("wrong PFT label in cohort in SP mode: $(fates_pft) vs $(currentPatch.nocomp_pft_label)")
            end
            if fates_pft == 0
                fates_endrun("PFT0 in SP mode")
            end

            assign_cohort_SP_properties(currentCohort, currentSite.sp_htop[fates_pft],
                currentSite.sp_tlai[fates_pft], currentSite.sp_tsai[fates_pft],
                currentPatch.area, ifalse)

            currentCohort = currentCohort.shorter
        end
        currentPatch = currentPatch.younger
    end
    return nothing
end

# ===========================================================================
# calculate_SP_properties
# ===========================================================================
"""
    calculate_SP_properties(htop, tlai, tsai, parea, pft, crown_damage, canopy_layer, vcmax25top)
        -> (leaf_c, dbh, cohort_n, c_area)

Invert the SP-mode LAI/SAI/height drivers into cohort dbh, density (`cohort_n`),
crown area, and leaf carbon. Mirrors the Fortran `calculate_SP_properties`
(returns the four out-args as a tuple).
"""
function calculate_SP_properties(htop::Real, tlai::Real, tsai::Real, parea::Real,
                                 pft::Integer, crown_damage::Integer,
                                 canopy_layer::Integer, vcmax25top::Real)
    # DBH from input height
    dbh, _ = h2d_allom(htop, pft)

    # Crown area assuming n = 1, spread = 1
    _, c_area = carea_allom(dbh, 1.0, 1.0, pft, crown_damage)

    # Canopy N assuming patch area is full
    cohort_n = parea / c_area

    # Correct c_area for the new n
    _, c_area = carea_allom(dbh, cohort_n, 1.0, pft, crown_damage)

    # Leaf carbon from the target tree LAI
    canopylai = zeros(nclmax)
    leaf_c = leafc_from_treelai(tlai, pft, c_area, cohort_n, canopy_layer, vcmax25top)

    # Validate inverse
    check_treelai = tree_lai(leaf_c, pft, c_area, cohort_n, canopy_layer, canopylai, vcmax25top)
    if abs(tlai - check_treelai) > area_error_2
        fates_endrun("error in validate treelai: $(tlai) vs $(check_treelai)")
    end

    # Correct small c_area / parea precision errors
    if abs(c_area - parea) > nearzero
        if abs(c_area - parea) < area_error_3
            oldcarea = c_area
            c_area = c_area - (c_area - parea)
            cohort_n = cohort_n * (c_area / oldcarea)
            if abs(c_area - parea) > nearzero
                fates_endrun("SPassign, c_area still broken: $(c_area - parea)")
            end
        else
            fates_endrun("SPassign, big error in c_area: $(c_area - parea), pft=$(pft)")
        end
    end

    return leaf_c, dbh, cohort_n, c_area
end

# ===========================================================================
# assign_cohort_SP_properties
# ===========================================================================
"""
    assign_cohort_SP_properties(currentCohort, htop, tlai, tsai, parea, init) -> leaf_c

Set the SP-mode cohort allometric properties from the inverted drivers. When
`init == ifalse`, also writes the leaf carbon into the cohort's PARTEH state.
Mirrors the Fortran `assign_cohort_SP_properties` (returns leaf_c).
"""
function assign_cohort_SP_properties(currentCohort::fates_cohort_type, htop::Real,
                                     tlai::Real, tsai::Real, parea::Real, init::Integer)
    if currentCohort.shorter !== nothing
        fates_endrun("SP mode has >1 cohort")
    end

    if init == itrue
        currentCohort.canopy_layer = 1
        currentCohort.vcmax25top = EDPftvarcon_inst[].vcmax25top[currentCohort.pft, 1]
    end

    leaf_c, dbh, cohort_n, c_area = calculate_SP_properties(htop, tlai, tsai, parea,
        currentCohort.pft, currentCohort.crowndamage, currentCohort.canopy_layer,
        currentCohort.vcmax25top)

    currentCohort.height = htop
    currentCohort.dbh = dbh
    currentCohort.n = cohort_n
    currentCohort.c_area = c_area
    currentCohort.treelai = tlai
    currentCohort.treesai = tsai

    if init == ifalse
        SetState!(currentCohort.prt, leaf_organ, carbon12_element, leaf_c, 1)
    end

    return leaf_c
end

# ===========================================================================
# SeedUpdate
# ===========================================================================
"""
    SeedUpdate(currentSite)

Flux of carbon (and nutrients) from plants into the seed pool: clonal
reproduction from storage of dying plants, release of reproductive tissue, plus
externally-supplied / dispersed seed rain. Mirrors the Fortran `SeedUpdate`.
"""
function SeedUpdate(currentSite::ed_site_type)
    site_disp_frac = zeros(numpft[])
    if hlm_seeddisp_cadence[] != fates_dispersal_cadence_none
        site_disp_frac[1:numpft[]] = EDPftvarcon_inst[].seed_dispersal_fraction[1:numpft[]]
    end

    for el in 1:num_elements[]
        site_seed_rain = zeros(numpft[])
        element_id = element_list[el]
        site_mass = currentSite.mass_balance[el]

        # Sum seed input for each PFT over all patches
        currentPatch = currentSite.oldest_patch
        while currentPatch !== nothing
            currentCohort = currentPatch.tallest
            while currentCohort !== nothing
                pft = currentCohort.pft

                store_m_to_repro = -GetState(currentCohort.prt, store_organ, element_id) *
                    EDPftvarcon_inst[].allom_frbstor_repro[pft] * currentCohort.dndt * years_per_day

                seed_prod = PRTReproRelease!(currentCohort.prt, repro_organ, element_id, 1.0)

                if element_id == carbon12_element
                    currentCohort.seed_prod = seed_prod
                end

                site_seed_rain[pft] += seed_prod * currentCohort.n + store_m_to_repro

                currentCohort = currentCohort.shorter
            end
            currentPatch = currentPatch.younger
        end

        if homogenize_seed_pfts
            tot = sum(site_seed_rain)
            site_seed_rain[1:numpft[]] .= tot / Float64(numpft[])
        end

        # Disperse seeds into the input flux arrays
        currentPatch = currentSite.oldest_patch
        while currentPatch !== nothing
            litt = currentPatch.litter[el]
            for pft in 1:numpft[]
                if currentSite.use_this_pft[pft] == itrue
                    litt.seed_in_local[pft] += site_seed_rain[pft] *
                        (1.0 - site_disp_frac[pft]) / area

                    if (ed_params().regeneration_model == TRS_regeneration ||
                        ed_params().regeneration_model == TRS_no_seedling_dyn) &&
                       prt_params.allom_dbh_maxheight[pft] > min_max_dbh_for_trees
                        litt.seed_decay[pft] = litt.seed_in_local[pft] *
                            (1.0 - EDPftvarcon_inst[].repro_frac_seed[pft])
                    end

                    if element_id == carbon12_element
                        seed_stoich = 1.0
                    elseif element_id == nitrogen_element
                        seed_stoich = currentPatch.nitr_repro_stoich[pft]
                    elseif element_id == phosphorus_element
                        seed_stoich = currentPatch.phos_repro_stoich[pft]
                    else
                        fates_endrun("undefined element specified while defining forced external seed mass flux")
                        seed_stoich = 0.0
                    end

                    seed_in_external = seed_stoich *
                        (currentSite.seed_in[pft] / area +
                         EDPftvarcon_inst[].seed_suppl[pft] * years_per_day)
                    litt.seed_in_extern[pft] += seed_in_external

                    site_mass.seed_in += seed_in_external * currentPatch.area
                end
            end
            currentPatch = currentPatch.younger
        end

        for pft in 1:numpft[]
            site_mass.seed_out += site_seed_rain[pft] * site_disp_frac[pft]
            currentSite.seed_out[pft] += site_seed_rain[pft] * site_disp_frac[pft]
        end
    end
    return nothing
end

# ===========================================================================
# SeedDecay
# ===========================================================================
"""
    SeedDecay(litt, currentPatch, bc_in)

Flux from the seed pool into leaf litter (and, with the TRS seedling-dynamics
model, seedling mortality from the germinated pool). Mirrors the Fortran
`SeedDecay`. The TRS seedling-dynamics path uses running-mean seedling
diagnostics held on the patch.
"""
function SeedDecay(litt::litter_type, currentPatch::fates_patch_type, bc_in)
    regen = ed_params().regeneration_model
    for pft in 1:numpft[]
        # Default seed decay (TRS off, or PFT cannot become a tree)
        if regen == default_regeneration ||
           prt_params.allom_dbh_maxheight[pft] < min_max_dbh_for_trees
            litt.seed_decay[pft] = litt.seed[pft] *
                EDPftvarcon_inst[].seed_decay_rate[pft] * years_per_day
        end

        # TRS tree path: add the non-seed reproductive biomass decay
        if (regen == TRS_regeneration || regen == TRS_no_seedling_dyn) &&
           prt_params.allom_dbh_maxheight[pft] > min_max_dbh_for_trees
            litt.seed_decay[pft] += litt.seed[pft] *
                EDPftvarcon_inst[].seed_decay_rate[pft] * years_per_day
        end

        # TRS with seedling dynamics: seedling mortality from light + moisture + bg
        if regen == TRS_regeneration &&
           prt_params.allom_dbh_maxheight[pft] > min_max_dbh_for_trees
            seedling_layer_par = GetMean(currentPatch.sdlng_mort_par) *
                megajoules_per_joule * sec_per_day * ed_params().sdlng_mort_par_timescale

            seedling_light_mort_rate = exp(
                EDPftvarcon_inst[].seedling_light_mort_a[pft] * seedling_layer_par +
                EDPftvarcon_inst[].seedling_light_mort_b[pft])

            seedling_mdds = GetMean(currentPatch.sdlng_mdd[pft].p)

            if seedling_mdds < EDPftvarcon_inst[].seedling_mdd_crit[pft]
                seedling_h2o_mort_rate = 0.0
            else
                seedling_h2o_mort_rate =
                    EDPftvarcon_inst[].seedling_h2o_mort_a[pft] * seedling_mdds^2 +
                    EDPftvarcon_inst[].seedling_h2o_mort_b[pft] * seedling_mdds +
                    EDPftvarcon_inst[].seedling_h2o_mort_c[pft]
            end

            litt.seed_germ_decay[pft] =
                litt.seed_germ[pft] * seedling_light_mort_rate +
                litt.seed_germ[pft] * seedling_h2o_mort_rate +
                litt.seed_germ[pft] * EDPftvarcon_inst[].background_seedling_mort[pft] * years_per_day
        else
            litt.seed_germ_decay[pft] = litt.seed_germ[pft] *
                EDPftvarcon_inst[].seed_decay_rate[pft] * years_per_day
        end
    end
    return nothing
end

# ===========================================================================
# SeedGermination
# ===========================================================================
"""
    SeedGermination(litt, cold_stat, drought_stat, bc_in, currentPatch)

Flux from the seed bank into the seedling/germinated pool. Germination is
suppressed in cold (for cold-deciduous PFTs) and in drought leaf-off states (for
drought-deciduous PFTs). Mirrors the Fortran `SeedGermination`.
"""
function SeedGermination(litt::litter_type, cold_stat::Integer, drought_stat,
                         bc_in, currentPatch::fates_patch_type)
    max_germination = 1.0  # KgC/m2/yr cap, Lischke et al. 2009
    regen = ed_params().regeneration_model

    for pft in 1:numpft[]
        if regen == default_regeneration ||
           regen == TRS_no_seedling_dyn ||
           prt_params.allom_dbh_maxheight[pft] < min_max_dbh_for_trees
            litt.seed_germ_in[pft] = min(litt.seed[pft] * EDPftvarcon_inst[].germination_rate[pft],
                                         max_germination) * years_per_day
        elseif regen == TRS_regeneration &&
               prt_params.allom_dbh_maxheight[pft] > min_max_dbh_for_trees
            seedling_layer_par = GetMean(currentPatch.seedling_layer_par24) *
                sec_per_day * megajoules_per_joule
            photoblastic_germ_modifier = seedling_layer_par /
                (seedling_layer_par + EDPftvarcon_inst[].par_crit_germ[pft])

            seedling_layer_smp = GetMean(currentPatch.sdlng_emerg_smp[pft].p)
            wetness_index = 1.0 / (seedling_layer_smp * (-1.0) * mpa_per_mm_suction)

            if seedling_layer_smp >= EDPftvarcon_inst[].seedling_psi_emerg[pft]
                seedling_emerg_rate = photoblastic_germ_modifier *
                    EDPftvarcon_inst[].a_emerg[pft] *
                    wetness_index^EDPftvarcon_inst[].b_emerg[pft]
            else
                seedling_emerg_rate = 0.0
            end

            litt.seed_germ_in[pft] = litt.seed[pft] * seedling_emerg_rate
        end

        # No germination when cold (cold-deciduous)
        if prt_params.season_decid[pft] == itrue &&
           (cold_stat == phen_cstat_nevercold || cold_stat == phen_cstat_iscold)
            litt.seed_germ_in[pft] = 0.0
        end

        # Halt germination in drought leaf-off states (drought-deciduous)
        if prt_params.stress_decid[pft] == ihard_stress_decid ||
           prt_params.stress_decid[pft] == isemi_stress_decid
            if drought_stat[pft] == phen_dstat_timeoff ||
               drought_stat[pft] == phen_dstat_moistoff ||
               drought_stat[pft] == phen_dstat_pshed
                litt.seed_germ_in[pft] = 0.0
            end
        end
    end
    return nothing
end

# ===========================================================================
# recruitment
# ===========================================================================
"""
    recruitment(currentSite, currentPatch, bc_in)

Spawn new cohorts of juveniles of each PFT. The number of recruits is limited by
the most-limiting available element in the germinated-seed pool (carbon, and
under CNP, nitrogen/phosphorus). Mirrors the Fortran `recruitment`.
"""
function recruitment(currentSite::ed_site_type, currentPatch::fates_patch_type, bc_in)
    recruitstatus = 1

    edp = EDPftvarcon_inst[]
    regen = ed_params().regeneration_model

    for ft in 1:numpft[]
        # PFT permitted on this patch? (prescribed biogeography / nocomp / crop LU)
        use_this_pft = false
        if currentSite.use_this_pft[ft] == itrue &&
           (hlm_use_nocomp[] == ifalse || ft == currentPatch.nocomp_pft_label)
            use_this_pft = true
        end

        if currentPatch.land_use_label != nocomp_bareground_land
            if hlm_use_luh[] == itrue && is_crop[currentPatch.land_use_label]
                if ed_params().crop_lu_pft_vector[currentPatch.land_use_label] == ft
                    use_this_pft = true
                else
                    use_this_pft = false
                end
            end
        end

        if !use_this_pft
            continue
        end

        height = edp.hgt_min[ft]
        stem_drop_fraction = prt_params.phen_stem_drop_fraction[ft]
        fnrt_drop_fraction = prt_params.phen_fnrt_drop_fraction[ft]
        l2fr = currentSite.rec_l2fr[ft, currentPatch.ncl_p]
        crowndamage = 1

        dbh, _ = h2d_allom(height, ft)

        efleaf_coh = 1.0
        effnrt_coh = 1.0
        efstem_coh = 1.0
        leaf_status = leaves_on

        # Cold deciduous: leaves off when site is cold
        if prt_params.season_decid[ft] == itrue &&
           (currentSite.cstatus == phen_cstat_nevercold || currentSite.cstatus == phen_cstat_iscold)
            efleaf_coh = 0.0
            effnrt_coh = 1.0 - fnrt_drop_fraction
            efstem_coh = 1.0 - stem_drop_fraction
            leaf_status = leaves_off
        end

        # Drought deciduous: leaf status consistent with elongation factor
        if prt_params.stress_decid[ft] == ihard_stress_decid ||
           prt_params.stress_decid[ft] == isemi_stress_decid
            efleaf_coh = currentSite.elong_factor[ft]
            effnrt_coh = 1.0 - (1.0 - efleaf_coh) * fnrt_drop_fraction
            efstem_coh = 1.0 - (1.0 - efleaf_coh) * stem_drop_fraction
            leaf_status = efleaf_coh > 0.0 ? leaves_on : leaves_off
        end

        # Live pools (target allometry)
        c_leaf, _ = bleaf(dbh, ft, crowndamage, init_recruit_trim, efleaf_coh)
        c_fnrt, _ = bfineroot(dbh, ft, init_recruit_trim, l2fr, effnrt_coh)
        _, c_sapw, _ = bsap_allom(dbh, ft, crowndamage, init_recruit_trim, efstem_coh)
        c_agw, _ = bagw_allom(dbh, ft, crowndamage, efstem_coh)
        c_bgw, _ = bbgw_allom(dbh, ft, efstem_coh)
        c_struct, _ = bdead_allom(c_agw, c_bgw, c_sapw, ft)
        c_store, _ = bstore_allom(dbh, ft, crowndamage, init_recruit_trim)

        # Limiting-element-based number density
        if hlm_use_ed_prescribed_phys[] == ifalse || edp.prescribed_recruitment[ft] < 0.0
            cohort_n = 1.0e20

            for el in 1:num_elements[]
                element_id = element_list[el]
                if element_id == carbon12_element
                    mass_demand = c_struct + c_leaf + c_fnrt + c_sapw + c_store
                elseif element_id == nitrogen_element
                    mass_demand =
                        c_struct * prt_params.nitr_stoich_p1[ft, prt_params.organ_param_id[struct_organ]] +
                        c_leaf * prt_params.nitr_stoich_p1[ft, prt_params.organ_param_id[leaf_organ]] +
                        c_fnrt * prt_params.nitr_stoich_p1[ft, prt_params.organ_param_id[fnrt_organ]] +
                        c_sapw * prt_params.nitr_stoich_p1[ft, prt_params.organ_param_id[sapw_organ]] +
                        StorageNutrientTarget(ft, element_id,
                            c_leaf * prt_params.nitr_stoich_p1[ft, prt_params.organ_param_id[leaf_organ]],
                            c_fnrt * prt_params.nitr_stoich_p1[ft, prt_params.organ_param_id[fnrt_organ]],
                            c_sapw * prt_params.nitr_stoich_p1[ft, prt_params.organ_param_id[sapw_organ]],
                            c_struct * prt_params.nitr_stoich_p1[ft, prt_params.organ_param_id[struct_organ]])
                elseif element_id == phosphorus_element
                    mass_demand =
                        c_struct * prt_params.phos_stoich_p1[ft, prt_params.organ_param_id[struct_organ]] +
                        c_leaf * prt_params.phos_stoich_p1[ft, prt_params.organ_param_id[leaf_organ]] +
                        c_fnrt * prt_params.phos_stoich_p1[ft, prt_params.organ_param_id[fnrt_organ]] +
                        c_sapw * prt_params.phos_stoich_p1[ft, prt_params.organ_param_id[sapw_organ]] +
                        StorageNutrientTarget(ft, element_id,
                            c_leaf * prt_params.phos_stoich_p1[ft, prt_params.organ_param_id[leaf_organ]],
                            c_fnrt * prt_params.phos_stoich_p1[ft, prt_params.organ_param_id[fnrt_organ]],
                            c_sapw * prt_params.phos_stoich_p1[ft, prt_params.organ_param_id[sapw_organ]],
                            c_struct * prt_params.phos_stoich_p1[ft, prt_params.organ_param_id[struct_organ]])
                else
                    fates_endrun("Undefined element type in recruitment")
                    mass_demand = 0.0
                end

                if regen == default_regeneration ||
                   regen == TRS_no_seedling_dyn ||
                   prt_params.allom_dbh_maxheight[ft] < min_max_dbh_for_trees
                    mass_avail = currentPatch.area * currentPatch.litter[el].seed_germ[ft]
                elseif regen == TRS_regeneration &&
                       prt_params.allom_dbh_maxheight[ft] > min_max_dbh_for_trees
                    sdlng2sap_par = GetMean(currentPatch.sdlng2sap_par) *
                        sec_per_day * megajoules_per_joule
                    mass_avail = currentPatch.area * currentPatch.litter[el].seed_germ[ft] *
                        edp.seedling_light_rec_a[ft] *
                        sdlng2sap_par^edp.seedling_light_rec_b[ft]

                    ilayer_seedling_root = argmin(abs.(bc_in.z_sisl .- edp.seedling_root_depth[ft]))
                    seedling_layer_smp = bc_in.smp_sl[ilayer_seedling_root]
                    if seedling_layer_smp < edp.seedling_psi_crit[ft]
                        mass_avail = 0.0
                    end
                else
                    mass_avail = 0.0
                end

                cohort_n = min(cohort_n, mass_avail / mass_demand)
            end
        else
            cohort_n = currentPatch.area * edp.prescribed_recruitment[ft] * hlm_freq_day[]
        end

        # Only allocate a new cohort if there is a reasonable amount of it
        if cohort_n > min_n_safemath
            prt = InitPRTObject!()

            for el in 1:num_elements[]
                element_id = element_list[el]
                if element_id == carbon12_element
                    m_struct = c_struct; m_leaf = c_leaf; m_fnrt = c_fnrt
                    m_sapw = c_sapw; m_store = c_store; m_repro = 0.0
                elseif element_id == nitrogen_element
                    m_struct = c_struct * prt_params.nitr_stoich_p1[ft, prt_params.organ_param_id[struct_organ]]
                    m_leaf   = c_leaf * prt_params.nitr_stoich_p1[ft, prt_params.organ_param_id[leaf_organ]]
                    m_fnrt   = c_fnrt * prt_params.nitr_stoich_p1[ft, prt_params.organ_param_id[fnrt_organ]]
                    m_sapw   = c_sapw * prt_params.nitr_stoich_p1[ft, prt_params.organ_param_id[sapw_organ]]
                    m_store  = StorageNutrientTarget(ft, element_id, m_leaf, m_fnrt, m_sapw, m_struct)
                    m_repro  = 0.0
                elseif element_id == phosphorus_element
                    m_struct = c_struct * prt_params.phos_stoich_p1[ft, prt_params.organ_param_id[struct_organ]]
                    m_leaf   = c_leaf * prt_params.phos_stoich_p1[ft, prt_params.organ_param_id[leaf_organ]]
                    m_fnrt   = c_fnrt * prt_params.phos_stoich_p1[ft, prt_params.organ_param_id[fnrt_organ]]
                    m_sapw   = c_sapw * prt_params.phos_stoich_p1[ft, prt_params.organ_param_id[sapw_organ]]
                    m_store  = StorageNutrientTarget(ft, element_id, m_leaf, m_fnrt, m_sapw, m_struct)
                    m_repro  = 0.0
                else
                    m_struct = m_leaf = m_fnrt = m_sapw = m_store = m_repro = 0.0
                end

                if hlm_parteh_mode[] == prt_carbon_allom_hyp ||
                   hlm_parteh_mode[] == prt_cnp_flex_allom_hyp
                    SetState!(prt, leaf_organ, element_id, m_leaf, 1)
                    for iage in 2:nleafage[]
                        SetState!(prt, leaf_organ, element_id, 0.0, iage)
                    end
                    SetState!(prt, fnrt_organ, element_id, m_fnrt)
                    SetState!(prt, sapw_organ, element_id, m_sapw)
                    SetState!(prt, store_organ, element_id, m_store)
                    SetState!(prt, struct_organ, element_id, m_struct)
                    SetState!(prt, repro_organ, element_id, m_repro)
                else
                    fates_endrun("Unspecified PARTEH module during create_cohort")
                end

                site_mass = currentSite.mass_balance[el]

                if hlm_use_ed_prescribed_phys[] == itrue && edp.prescribed_recruitment[ft] >= 0.0
                    site_mass.flux_generic_in += cohort_n *
                        (m_struct + m_leaf + m_fnrt + m_sapw + m_store + m_repro)
                    site_mass.flux_generic_out += currentPatch.area * currentPatch.litter[el].seed_germ[ft]
                    currentPatch.litter[el].seed_germ[ft] = 0.0
                else
                    currentPatch.litter[el].seed_germ[ft] -= cohort_n / currentPatch.area *
                        (m_struct + m_leaf + m_fnrt + m_sapw + m_store + m_repro)
                end
            end

            CheckInitialConditions(prt)

            create_cohort(currentSite, currentPatch, ft, cohort_n, height, 0.0, dbh,
                prt, efleaf_coh, effnrt_coh, efstem_coh, leaf_status, recruitstatus,
                init_recruit_trim, 0.0, currentPatch.ncl_p, crowndamage, currentSite.spread, bc_in)

            currentSite.recruitment_rate[ft] += cohort_n
        end
    end
    return nothing
end

# ===========================================================================
# CWDInput
# ===========================================================================
"""
    CWDInput(currentSite, currentPatch, litt, bc_in)

Generate litter input fluxes from living-plant tissue turnover and from
non-disturbance-inducing mortality (including direct/indirect logging in the
understory). Mirrors the Fortran `CWDInput`. Note: number density has not yet
been reduced for mortality when this is called, so live-tree turnover already
captures the dying-tree contribution (avoiding double counting).
"""
function CWDInput(currentSite::ed_site_type, currentPatch::fates_patch_type,
                  litt::litter_type, bc_in)
    numlevsoil = currentSite.nlevsoil
    element_id = litt.element_id

    elflux_diags = currentSite.flux_diags.elem[element_pos[element_id]]
    site_mass    = currentSite.mass_balance[element_pos[element_id]]

    SF_val_CWD_frac = sf_params().SF_val_CWD_frac
    SF_val_CWD_frac_adj = zeros(ncwd)
    logging_export = ed_params().logging_export_frac

    currentCohort = currentPatch.shortest
    while currentCohort !== nothing
        pft = currentCohort.pft
        set_root_fraction(currentSite.rootfrac_scr, pft, currentSite.zi_soil;
                          max_nlevroot=bc_in.max_rooting_depth_index_col)

        store_m_turnover = GetTurnover(currentCohort.prt, store_organ, element_id)
        fnrt_m_turnover  = GetTurnover(currentCohort.prt, fnrt_organ, element_id)
        repro_m_turnover = GetTurnover(currentCohort.prt, repro_organ, element_id)

        store_m = GetState(currentCohort.prt, store_organ, element_id)
        fnrt_m  = GetState(currentCohort.prt, fnrt_organ, element_id)
        repro_m = GetState(currentCohort.prt, repro_organ, element_id)

        if prt_params.woody[currentCohort.pft] == itrue
            leaf_m_turnover   = GetTurnover(currentCohort.prt, leaf_organ, element_id)
            sapw_m_turnover   = GetTurnover(currentCohort.prt, sapw_organ, element_id)
            struct_m_turnover = GetTurnover(currentCohort.prt, struct_organ, element_id)
            leaf_m   = GetState(currentCohort.prt, leaf_organ, element_id)
            sapw_m   = GetState(currentCohort.prt, sapw_organ, element_id)
            struct_m = GetState(currentCohort.prt, struct_organ, element_id)
        else
            leaf_m_turnover = GetTurnover(currentCohort.prt, leaf_organ, element_id) +
                              GetTurnover(currentCohort.prt, sapw_organ, element_id) +
                              GetTurnover(currentCohort.prt, struct_organ, element_id)
            sapw_m_turnover   = 0.0
            struct_m_turnover = 0.0
            leaf_m = GetState(currentCohort.prt, leaf_organ, element_id) +
                     GetState(currentCohort.prt, sapw_organ, element_id) +
                     GetState(currentCohort.prt, struct_organ, element_id)
            sapw_m   = 0.0
            struct_m = 0.0
        end

        plant_dens = currentCohort.n / currentPatch.area

        # PART 1: non-mortal tissue turnover
        elflux_diags.surf_fine_litter_input[pft] +=
            (leaf_m_turnover + repro_m_turnover) * currentCohort.n

        root_fines_tot = (fnrt_m_turnover + store_m_turnover) * plant_dens

        for dcmpy in 1:ndcmpy
            dcmpy_frac = GetDecompyFrac(pft, leaf_organ, dcmpy)
            litt.leaf_fines_in[dcmpy] += (leaf_m_turnover + repro_m_turnover) * plant_dens * dcmpy_frac

            dcmpy_frac = GetDecompyFrac(pft, fnrt_organ, dcmpy)
            for ilyr in 1:numlevsoil
                litt.root_fines_in[dcmpy, ilyr] +=
                    currentSite.rootfrac_scr[ilyr] * root_fines_tot * dcmpy_frac
            end
        end

        elflux_diags.root_litter_input[pft] +=
            (fnrt_m_turnover + store_m_turnover) * currentCohort.n

        adjust_SF_CWD_frac(currentCohort.dbh, ncwd, SF_val_CWD_frac, SF_val_CWD_frac_adj)

        for c in 1:ncwd
            litt.ag_cwd_in[c] += (sapw_m_turnover + struct_m_turnover) *
                SF_val_CWD_frac_adj[c] * plant_dens * prt_params.allom_agb_frac[pft]
            elflux_diags.cwd_ag_input[c] += (struct_m_turnover + sapw_m_turnover) *
                SF_val_CWD_frac_adj[c] * prt_params.allom_agb_frac[pft] * currentCohort.n

            bg_cwd_tot = (sapw_m_turnover + struct_m_turnover) *
                SF_val_CWD_frac_adj[c] * plant_dens * (1.0 - prt_params.allom_agb_frac[pft])
            for ilyr in 1:numlevsoil
                litt.bg_cwd_in[c, ilyr] += bg_cwd_tot * currentSite.rootfrac_scr[ilyr]
            end
            elflux_diags.cwd_bg_input[c] += bg_cwd_tot * currentPatch.area
        end

        # PART 2: non-disturbance-inducing mortality
        dead_n = -1.0 * currentCohort.dndt / currentPatch.area * years_per_day

        if currentCohort.canopy_layer > 1
            dead_n_dlogging = currentCohort.lmort_direct * currentCohort.n / currentPatch.area
            dead_n_ilogging = (currentCohort.lmort_collateral + currentCohort.lmort_infra) *
                currentCohort.n / currentPatch.area
        else
            dead_n_dlogging = 0.0
            dead_n_ilogging = 0.0
        end
        dead_n_natural = dead_n - dead_n_dlogging - dead_n_ilogging

        root_fines_tot = dead_n * (fnrt_m +
            store_m * (1.0 - EDPftvarcon_inst[].allom_frbstor_repro[pft]))

        for dcmpy in 1:ndcmpy
            dcmpy_frac = GetDecompyFrac(pft, leaf_organ, dcmpy)
            litt.leaf_fines_in[dcmpy] += (leaf_m + repro_m) * dead_n * dcmpy_frac

            dcmpy_frac = GetDecompyFrac(pft, fnrt_organ, dcmpy)
            for ilyr in 1:numlevsoil
                litt.root_fines_in[dcmpy, ilyr] +=
                    root_fines_tot * currentSite.rootfrac_scr[ilyr] * dcmpy_frac
            end
        end

        elflux_diags.surf_fine_litter_input[pft] += (leaf_m + repro_m) * dead_n * currentPatch.area
        elflux_diags.root_litter_input[pft] += root_fines_tot * currentPatch.area

        trunk_wood = 0.0
        for c in 1:ncwd
            bg_cwd_tot = (struct_m + sapw_m) * SF_val_CWD_frac_adj[c] * dead_n *
                (1.0 - prt_params.allom_agb_frac[pft])
            for ilyr in 1:numlevsoil
                litt.bg_cwd_in[c, ilyr] += currentSite.rootfrac_scr[ilyr] * bg_cwd_tot
            end
            elflux_diags.cwd_bg_input[c] += bg_cwd_tot * currentPatch.area

            if c == ncwd
                trunk_wood = (struct_m + sapw_m) * SF_val_CWD_frac_adj[c] * dead_n_dlogging *
                    prt_params.allom_agb_frac[pft]

                site_mass.wood_product_harvest[pft] += trunk_wood * currentPatch.area * logging_export

                litt.ag_cwd_in[c] += trunk_wood * (1.0 - logging_export)
                elflux_diags.cwd_ag_input[c] += trunk_wood * (1.0 - logging_export) * currentPatch.area

                litt.ag_cwd_in[c] += (struct_m + sapw_m) * SF_val_CWD_frac_adj[c] *
                    (dead_n_natural + dead_n_ilogging) * prt_params.allom_agb_frac[pft]
                elflux_diags.cwd_ag_input[c] += (struct_m + sapw_m) * SF_val_CWD_frac_adj[c] *
                    (dead_n_natural + dead_n_ilogging) * currentPatch.area * prt_params.allom_agb_frac[pft]
            else
                litt.ag_cwd_in[c] += (struct_m + sapw_m) * SF_val_CWD_frac_adj[c] * dead_n *
                    prt_params.allom_agb_frac[pft]
                elflux_diags.cwd_ag_input[c] += SF_val_CWD_frac_adj[c] * dead_n *
                    (struct_m + sapw_m) * currentPatch.area * prt_params.allom_agb_frac[pft]
            end
        end

        # Resource-management diagnostics (carbon only)
        if element_id == carbon12_element
            rm = currentSite.resources_management
            rm.delta_litter_stock += (leaf_m + fnrt_m + store_m) *
                (dead_n_ilogging + dead_n_dlogging) * currentPatch.area
            rm.delta_biomass_stock += (leaf_m + fnrt_m + store_m) *
                (dead_n_ilogging + dead_n_dlogging) * currentPatch.area
            rm.trunk_product_site += trunk_wood * logging_export * currentPatch.area

            for c in 1:ncwd
                rm.delta_litter_stock += (struct_m + sapw_m) * SF_val_CWD_frac_adj[c] *
                    (dead_n_natural + dead_n_ilogging) * currentPatch.area
                rm.delta_biomass_stock += (struct_m + sapw_m) * SF_val_CWD_frac_adj[c] *
                    dead_n * currentPatch.area
            end

            rm.delta_individual += (dead_n_dlogging + dead_n_ilogging) *
                hlm_freq_day[] * currentPatch.area
        end

        currentCohort = currentCohort.taller
    end
    return nothing
end

# ===========================================================================
# fragmentation_scaler
# ===========================================================================
"""
    fragmentation_scaler(currentPatch, bc_in)

Compute the per-soil-layer CWD/litter fragmentation scaler for a patch. By
default uses the HLM-supplied temperature and moisture decomposition fractions.
Mirrors the Fortran `fragmentation_scaler` (writes `currentPatch.fragmentation_scaler`).
"""
function fragmentation_scaler(currentPatch::fates_patch_type, bc_in)
    use_century_tfunc = false
    use_hlm_soil_scalar = true
    tfrz = t_water_freeze_k_1atm

    catanf(t1) = 11.75 + (29.7 / pi_const) * atan(pi_const * 0.031 * (t1 - 15.4))
    catanf_30 = catanf(30.0)

    if currentPatch.nocomp_pft_label != nocomp_bareground
        if use_hlm_soil_scalar
            currentPatch.fragmentation_scaler .=
                min.(1.0, max.(0.0, bc_in.t_scalar_sisl .* bc_in.w_scalar_sisl))
        else
            if !use_century_tfunc
                tv = GetMean(currentPatch.tveg24)
                if tv >= tfrz
                    t_scalar = ed_params().q10_mr^((tv - (tfrz + 25.0)) / 10.0)
                else
                    t_scalar = (ed_params().q10_mr^(-25.0 / 10.0)) *
                               (ed_params().q10_froz^((tv - tfrz) / 10.0))
                end
            else
                t_scalar = max(catanf(GetMean(currentPatch.tveg24) - tfrz) / catanf_30, 0.01)
            end

            w_scalar = sum(view(currentPatch.btran_ft, 1:numpft[])) / Float64(numpft[])
            currentPatch.fragmentation_scaler .= min(1.0, max(0.0, t_scalar * w_scalar))
        end
    end
    return nothing
end

# ===========================================================================
# CWDOut
# ===========================================================================
"""
    CWDOut(litt, fragmentation_scaler, nlev_eff_decomp)

Compute the fragmentation (decomposition-input) fluxes out of the CWD and fine
litter pools for a patch, over the active soil layers. Mirrors the Fortran
`CWDOut`. Above-ground litter uses the top-soil-layer (index 1) scaler.
"""
function CWDOut(litt::litter_type, frag_scaler::AbstractVector, nlev_eff_decomp::Integer)
    SF_val_max_decomp = sf_params().SF_val_max_decomp
    soil_layer_index = 1
    dead_leaves_idx = dead_leaves(fuel_classes)

    for c in 1:ncwd
        litt.ag_cwd_frag[c] = litt.ag_cwd[c] * SF_val_max_decomp[c] *
            years_per_day * frag_scaler[soil_layer_index]
        for ilyr in 1:nlev_eff_decomp
            litt.bg_cwd_frag[c, ilyr] = litt.bg_cwd[c, ilyr] * SF_val_max_decomp[c] *
                years_per_day * frag_scaler[ilyr]
        end
    end

    for dcmpy in 1:ndcmpy
        litt.leaf_fines_frag[dcmpy] = litt.leaf_fines[dcmpy] * years_per_day *
            SF_val_max_decomp[dead_leaves_idx] * frag_scaler[soil_layer_index]
        for ilyr in 1:nlev_eff_decomp
            litt.root_fines_frag[dcmpy, ilyr] = litt.root_fines[dcmpy, ilyr] * years_per_day *
                SF_val_max_decomp[dead_leaves_idx] * frag_scaler[ilyr]
        end
    end
    return nothing
end

# ===========================================================================
# UpdateRecruitL2FR / UpdateRecruitStoich / SetRecruitL2FR  (CNP only)
# ===========================================================================
"""
    UpdateRecruitL2FR(csite)

CNP-only: update the recruit target leaf:fine-root multiplier (`rec_l2fr`) as an
exponential moving average of the l2fr of recently-recruited cohorts. No-op
unless `hlm_parteh_mode == prt_cnp_flex_allom_hyp`. Mirrors the Fortran routine.
"""
function UpdateRecruitL2FR(csite::ed_site_type)
    max_delta = 5.0
    smth_wgt = 1.0 / 300.0
    max_count = 3

    if hlm_parteh_mode[] != prt_cnp_flex_allom_hyp
        return nothing
    end

    rec_n = zeros(maxpft, nclmax)
    rec_l2fr0 = zeros(maxpft, nclmax)

    cpatch = csite.youngest_patch
    while cpatch !== nothing
        rec_count = zeros(Int, maxpft, nclmax)

        ccohort = cpatch.shortest
        while ccohort !== nothing
            ft = ccohort.pft
            cl = ccohort.canopy_layer
            dbh_min, _ = h2d_allom(EDPftvarcon_inst[].hgt_min[ft], ft)

            if !ccohort.isnew
                if rec_count[ft, cl] <= max_count && ccohort.dbh - dbh_min < max_delta
                    rec_count[ft, cl] += 1
                    rec_n[ft, cl] += ccohort.n
                    rec_l2fr0[ft, cl] += ccohort.n * ccohort.l2fr
                end
            end

            ccohort = ccohort.taller
        end
        cpatch = cpatch.older
    end

    for cl in 1:nclmax
        for ft in 1:numpft[]
            if rec_n[ft, cl] > nearzero
                rec_l2fr0[ft, cl] = rec_l2fr0[ft, cl] / rec_n[ft, cl]
                csite.rec_l2fr[ft, cl] =
                    (1.0 - smth_wgt) * csite.rec_l2fr[ft, cl] + smth_wgt * rec_l2fr0[ft, cl]
            end
        end
    end
    return nothing
end

"""
    UpdateRecruitStoich(csite)

CNP-only: update the new-recruit total N:C and P:C stoichiometry (patch- and
cohort-level) from the updated `rec_l2fr`. No-op unless CNP mode. Mirrors the
Fortran routine.
"""
function UpdateRecruitStoich(csite::ed_site_type)
    if hlm_parteh_mode[] != prt_cnp_flex_allom_hyp
        return nothing
    end

    cpatch = csite.youngest_patch
    while cpatch !== nothing
        cl = cpatch.ncl_p

        for ft in 1:numpft[]
            rec_l2fr_pft = csite.rec_l2fr[ft, cl]
            cpatch.nitr_repro_stoich[ft] = NewRecruitTotalStoichiometry(ft, rec_l2fr_pft, nitrogen_element)
            cpatch.phos_repro_stoich[ft] = NewRecruitTotalStoichiometry(ft, rec_l2fr_pft, phosphorus_element)
        end

        ccohort = cpatch.shortest
        while ccohort !== nothing
            rec_l2fr_pft = csite.rec_l2fr[ccohort.pft, cl]
            ccohort.nc_repro = NewRecruitTotalStoichiometry(ccohort.pft, rec_l2fr_pft, nitrogen_element)
            ccohort.pc_repro = NewRecruitTotalStoichiometry(ccohort.pft, rec_l2fr_pft, phosphorus_element)
            ccohort = ccohort.taller
        end

        cpatch = cpatch.older
    end
    return nothing
end

"""
    SetRecruitL2FR(csite)

CNP-only: set the l2fr of brand-new cohorts to the recruit `rec_l2fr` for their
PFT/canopy-layer. No-op unless CNP mode. Mirrors the Fortran routine.
"""
function SetRecruitL2FR(csite::ed_site_type)
    if hlm_parteh_mode[] != prt_cnp_flex_allom_hyp
        return nothing
    end

    cpatch = csite.youngest_patch
    while cpatch !== nothing
        ccohort = cpatch.shortest
        while ccohort !== nothing
            if ccohort.isnew
                ft = ccohort.pft
                cl = ccohort.canopy_layer
                ccohort.l2fr = csite.rec_l2fr[ft, cl]
            end
            ccohort = ccohort.taller
        end
        cpatch = cpatch.older
    end
    return nothing
end
