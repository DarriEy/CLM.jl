# PRTLossFluxesMod.jl
# Julia port of FATES src/fates/parteh/PRTLossFluxesMod.F90
#
# The PARTEH loss-flux routines. They operate on the generic `prt_vartypes`
# state (from PRTGenericMod) and handle the two broad types of turnover:
#   1) event turnover  -- storms, deciduous leaf drop, herbivory, fire, damage
#   2) maintenance turnover -- constant background loss of evergreens + branchfall
# plus the flush of leaves from storage at bud-burst and the release of
# reproductive tissue.
#
# IMPORTANT POINTS (verbatim from the Fortran header):
#   Retranslocation is handled by a single flag per PFT. A deciduous plant does
#   not have maintenance leaf/fine-root turnover; an evergreen plant does not
#   have seasonal/stress phenology. So `turnover_*_retrans` means different
#   things per PFT: for evergreens it is retranslocation during maintenance
#   turnover; for deciduous it is during leaf drop.
#
#   THESE ROUTINES ONLY DEAL WITH LOSSES OF BIOMASS FROM PLANTS THAT SURVIVE AN
#   EVENT. If a plant dies, these routines do not handle its fluxes.
#
# Translation notes:
# - `fates_r8` -> Float64.
# - Fortran subroutine names preserved exactly.
# - `class(prt_vartypes) :: prt` -> `prt::AbstractPRTVartypes` (the merged
#   generic container subtypes this in PRTGenericMod).
# - The module-level singleton `prt_global` is a `Ref` here; obtained via
#   `get_prt_global()`. The 0-based `sp_organ_map` is read through
#   `sp_organ_map_get(g, organ, element)`.
# - The `associate(organ_map => prt_global%organ_map)` block becomes a local
#   binding `om = g.organ_map`.
# - `endrun(msg=errMsg(...))` -> `fates_endrun(...)`.
# - Mass-conservation flux math preserved EXACTLY (operation order kept,
#   including the Fortran quirks at lines 191/195 where the storage net_alloc/
#   val are read at `i_store_pos` while written at `i_pos` in PRTPhenologyFlush).
#
# Deps: PRTGenericMod (prt_vartypes/AbstractPRTVartypes, organ/element consts,
# get_prt_global, sp_organ_map_get), PRTParametersMod (prt_params),
# FatesConstantsMod (years_per_day, nearzero, itrue, un_initialized,
# check_initialized via PRTGenericMod), FatesGlobals (fates_endrun, fates_log).

# =====================================================================================

"""
    PRTPhenologyFlush!(prt, ipft, organ_id, c_store_transfer_frac)

Flush (leaves) from storage upon bud-burst. Leaves are implied, but the routine
allows other pools (fine-roots) to be flushed from storage as well. Carbon is
transferred from the storage pool to the organ pool by the fraction
`c_store_transfer_frac`; nutrients are then topped up toward their stoichiometric
target from storage. Mass moved out of storage is mirrored as a `net_alloc`
diagnostic in both the receiving organ and the donating storage pool.
"""
function PRTPhenologyFlush!(prt::AbstractPRTVartypes, ipft::Integer, organ_id::Integer,
                            c_store_transfer_frac::Real)
    g = get_prt_global()

    # We currently allow the flushing and drop of leaves and fine roots (always)
    # and sapwood and heartwood (non-woody PFTs only). If other organs should be
    # desired, those parameters and clauses need to be added.
    if (!(organ_id == leaf_organ || organ_id == fnrt_organ)) &&
       prt_params.woody[ipft] == itrue
        fates_endrun("When PFT is woody, deciduous drop and re-flushing are only " *
                     "allowed in leaves and fine roots. PFT=$(ipft) " *
                     "leaf_organ=$(leaf_organ) fnrt_organ=$(fnrt_organ) " *
                     "attempted organ=$(organ_id)")
    end

    if prt_params.organ_param_id[organ_id] < 1
        fates_endrun("Attempting to flush an organ that does not have a " *
                     "stoichiometry defined. global organ id (leaf=1)=$(organ_id)")
    end

    if g.hyp_id <= 2
        i_leaf_pos  = 1   # also used for sapwood and structural for grass
        i_store_pos = 1   # hypothesis 1/2 only have 1 storage pool
    else
        fates_endrun("You picked a hypothesis that has not defined how and where " *
                     "flushing interacts with the storage pool (multiple storage pools).")
    end

    om = g.organ_map

    # Flush carbon variables first, as their transfer rates from storage depend
    # on the fraction passed in by the argument. After the values are updated, we
    # can then identify the stoichiometry targets which govern the nutrient fluxes.
    for i_var_of_organ in 1:om[organ_id].num_vars
        i_var = om[organ_id].var_id[i_var_of_organ]
        element_id = g.state_descriptor[i_var].element_id

        # This will filter IN all carbon related variables
        if element_id == carbon12_element
            # Variable id of the storage pool for this element (carbon12)
            i_store = sp_organ_map_get(g, store_organ, element_id)

            for i_pos in 1:i_leaf_pos
                # Mass transferred out of storage into the pool of interest
                mass_transfer = prt.variables[i_store].val[i_store_pos] *
                                c_store_transfer_frac

                # Increment the c pool of interest's allocation flux
                prt.variables[i_var].net_alloc[i_pos] =
                    prt.variables[i_var].net_alloc[i_pos] + mass_transfer

                # Update the c pool
                prt.variables[i_var].val[i_pos] =
                    prt.variables[i_var].val[i_pos] + mass_transfer

                # Increment the storage pool's allocation flux
                prt.variables[i_store].net_alloc[i_pos] =
                    prt.variables[i_store].net_alloc[i_store_pos] - mass_transfer

                # Update the storage c pool
                prt.variables[i_store].val[i_pos] =
                    prt.variables[i_store].val[i_store_pos] - mass_transfer
            end
        end
    end

    # Variable index for leaf carbon used to calculate the targets for nutrient flushing
    i_cvar = sp_organ_map_get(g, organ_id, carbon12_element)
    if i_cvar < 1
        fates_endrun("Could not determine the carbon var id during flushing")
    end

    # Transfer in other elements (nutrients)
    for i_var_of_organ in 1:om[organ_id].num_vars
        i_var = om[organ_id].var_id[i_var_of_organ]
        element_id = g.state_descriptor[i_var].element_id

        # This will filter OUT all carbon related elements
        if !(element_id == carbon12_element)
            # Variable id of the storage pool for this element
            i_store = sp_organ_map_get(g, store_organ, element_id)

            # Calculate the stoichiometry with C for this element
            if element_id == nitrogen_element
                target_stoich = prt_params.nitr_stoich_p1[ipft, prt_params.organ_param_id[organ_id]]
            elseif element_id == phosphorus_element
                target_stoich = prt_params.phos_stoich_p1[ipft, prt_params.organ_param_id[organ_id]]
            else
                fates_endrun("Trying to calculate nutrient flushing target for " *
                             "element that DNE. organ=$(organ_id) element=$(element_id)")
            end

            for i_pos in 1:i_leaf_pos
                # The target quantity for this element is based on the amount of carbon
                sp_target = prt.variables[i_cvar].val[i_pos] * target_stoich

                sp_demand = max(0.0, sp_target - prt.variables[i_var].val[i_pos])

                # Assume that all of the storage is transferrable
                mass_transfer = min(sp_demand, prt.variables[i_store].val[i_store_pos])

                # Increment the pool of interest
                prt.variables[i_var].net_alloc[i_pos] =
                    prt.variables[i_var].net_alloc[i_pos] + mass_transfer

                # Update the pool
                prt.variables[i_var].val[i_pos] =
                    prt.variables[i_var].val[i_pos] + mass_transfer

                # Increment the store pool allocation diagnostic
                prt.variables[i_store].net_alloc[i_store_pos] =
                    prt.variables[i_store].net_alloc[i_store_pos] - mass_transfer

                # Update the store pool
                prt.variables[i_store].val[i_store_pos] =
                    prt.variables[i_store].val[i_store_pos] - mass_transfer
            end
        end
    end

    return nothing
end

# =====================================================================================

"""
    PRTBurnLosses!(prt, organ_id, mass_fraction)

Remove a `mass_fraction` of every element/position in `organ_id` due to burn.
There is no retranslocation: all burned mass leaves the plant (a non-PARTEH fire
model decides its destiny). Tracked in the `burned` diagnostic.
"""
function PRTBurnLosses!(prt::AbstractPRTVartypes, organ_id::Integer, mass_fraction::Real)
    g = get_prt_global()
    om = g.organ_map

    # All state variables associated with this particular organ
    for i_var_of_organ in 1:om[organ_id].num_vars
        i_var = om[organ_id].var_id[i_var_of_organ]

        for i_pos in 1:g.state_descriptor[i_var].num_pos
            # The mass that is leaving the plant
            burned_mass = mass_fraction * prt.variables[i_var].val[i_pos]

            # Track the amount of mass being burned (+ is amount lost)
            prt.variables[i_var].burned[i_pos] =
                prt.variables[i_var].burned[i_pos] + burned_mass

            # Update the state of the pool to reflect the mass lost
            prt.variables[i_var].val[i_pos] =
                prt.variables[i_var].val[i_pos] - burned_mass
        end
    end

    return nothing
end

# =====================================================================================

"""
    PRTDamageLosses!(prt, organ_id, mass_fraction)

Remove a `mass_fraction` of every element/position in `organ_id` due to damage.
There is no retranslocation: all damaged mass leaves the plant (a non-PARTEH
damage model decides its destiny). Tracked in the `damaged` diagnostic.
"""
function PRTDamageLosses!(prt::AbstractPRTVartypes, organ_id::Integer, mass_fraction::Real)
    g = get_prt_global()
    om = g.organ_map

    for i_var_of_organ in 1:om[organ_id].num_vars
        i_var = om[organ_id].var_id[i_var_of_organ]

        for i_pos in 1:g.state_descriptor[i_var].num_pos
            # The mass that is leaving the plant
            damaged_mass = mass_fraction * prt.variables[i_var].val[i_pos]

            # Track the amount of mass being lost (+ is amount lost)
            prt.variables[i_var].damaged[i_pos] =
                prt.variables[i_var].damaged[i_pos] + damaged_mass

            # Update the state of the pool to reflect the mass lost
            prt.variables[i_var].val[i_pos] =
                prt.variables[i_var].val[i_pos] - damaged_mass
        end
    end

    return nothing
end

# =====================================================================================

"""
    PRTReproRelease!(prt, organ_id, element_id, mass_fraction) -> mass_out

Release a `mass_fraction` of reproductive tissue (`organ_id` must be
`repro_organ`) for the given `element_id`. There is no retranslocation and no
dedicated flux: the released mass is returned as `mass_out`, and to avoid a
mass-imbalance check failure `val0` is reset so that `val - val0 == net_alloc`.
Returns the total released mass [kg] (Fortran `intent(out) :: mass_out`).
"""
function PRTReproRelease!(prt::AbstractPRTVartypes, organ_id::Integer,
                          element_id::Integer, mass_fraction::Real)
    g = get_prt_global()

    if organ_id != repro_organ
        fates_endrun("Reproductive tissue releases were called for a " *
                     "non-reproductive organ.")
    end

    i_var = sp_organ_map_get(g, organ_id, element_id)

    # Reproductive mass leaving the plant
    mass_out = 0.0

    for i_pos in 1:g.state_descriptor[i_var].num_pos
        # The mass that is leaving the plant
        mass_out = mass_out + mass_fraction * prt.variables[i_var].val[i_pos]

        # Update the state of the pool to reflect the mass lost
        prt.variables[i_var].val[i_pos] = prt.variables[i_var].val[i_pos] -
            (mass_fraction * prt.variables[i_var].val[i_pos])

        # Update the val0 (because we don't give this dedicated flux). A hack.
        prt.variables[i_var].val0[i_pos] = prt.variables[i_var].val[i_pos] -
            prt.variables[i_var].net_alloc[i_pos]
    end

    return mass_out
end

# =====================================================================================

"""
    PRTDeciduousTurnover!(prt, ipft, organ_id, mass_fraction)

Generic wrapper for deciduous turnover (leaf/fine-root drop). Validates the organ
choice for woody PFTs, then dispatches to
[`DeciduousTurnoverSimpleRetranslocation!`](@ref).
"""
function PRTDeciduousTurnover!(prt::AbstractPRTVartypes, ipft::Integer,
                               organ_id::Integer, mass_fraction::Real)
    # We currently allow the flushing and drop of leaves and fine roots (always)
    # and sapwood and heartwood (non-woody PFTs only).
    if (!(organ_id == leaf_organ || organ_id == fnrt_organ)) &&
       prt_params.woody[ipft] == itrue
        fates_endrun("When PFT is woody, deciduous drop and re-flushing are only " *
                     "allowed in leaves and fine roots. PFT=$(ipft) " *
                     "leaf_organ=$(leaf_organ) fnrt_organ=$(fnrt_organ) " *
                     "attempted organ=$(organ_id)")
    end

    DeciduousTurnoverSimpleRetranslocation!(prt, ipft, organ_id, mass_fraction)

    return nothing
end

# =====================================================================================

"""
    DeciduousTurnoverSimpleRetranslocation!(prt, ipft, organ_id, mass_fraction)

Calculate losses due to deciduous turnover. A fraction `retrans` of the lost
mass (carbon: 0; N/P: `turnover_{nitr,phos}_retrans`) is retranslocated to the
storage pool; the rest (`turnover_mass`) leaves the plant via the `turnover`
diagnostic. The retranslocated mass is debited from the organ's `net_alloc` and
credited to the storage pool's `val`/`net_alloc`.

NOTE (Fortran ALERT): no code limits the amount retranslocated into storage; the
maximum allowable storage may be over-shot.
"""
function DeciduousTurnoverSimpleRetranslocation!(prt::AbstractPRTVartypes, ipft::Integer,
                                                 organ_id::Integer, mass_fraction::Real)
    g = get_prt_global()
    om = g.organ_map

    if (organ_id == store_organ) || (organ_id == struct_organ) || (organ_id == sapw_organ)
        if prt_params.woody[ipft] == itrue
            fates_endrun("Deciduous turnover (leaf drop, etc) was specified for " *
                         "an unexpected organ. organ=$(organ_id)")
        end
    end

    if g.hyp_id <= 2
        i_store_pos = 1   # hypothesis 1&2 only have 1 storage pool
    else
        fates_endrun("You picked a hypothesis that has not defined how and where " *
                     "flushing interacts with the storage pool (multiple storage pools).")
    end

    for i_var_of_organ in 1:om[organ_id].num_vars
        i_var = om[organ_id].var_id[i_var_of_organ]
        element_id = g.state_descriptor[i_var].element_id

        if prt_params.organ_param_id[organ_id] < 1
            retrans = 0.0
        else
            if element_id == carbon12_element
                retrans = 0.0
            elseif element_id == nitrogen_element
                retrans = prt_params.turnover_nitr_retrans[ipft, prt_params.organ_param_id[organ_id]]
            elseif element_id == phosphorus_element
                retrans = prt_params.turnover_phos_retrans[ipft, prt_params.organ_param_id[organ_id]]
            else
                fates_endrun("Please add a new re-translocation clause to your " *
                             "organ x element combination. organ=$(leaf_organ) " *
                             "element=$(element_id)")
            end
        end

        # Variable id of the storage pool for this element
        store_var_id = sp_organ_map_get(g, store_organ, element_id)

        for i_pos in 1:g.state_descriptor[i_var].num_pos
            # The mass that is leaving the plant
            turnover_mass = (1.0 - retrans) * mass_fraction * prt.variables[i_var].val[i_pos]

            # The mass that is going towards storage
            retranslocated_mass = retrans * mass_fraction * prt.variables[i_var].val[i_pos]

            # Track the amount of mass being turned over (+ is amount lost)
            prt.variables[i_var].turnover[i_pos] =
                prt.variables[i_var].turnover[i_pos] + turnover_mass

            # Track the amount of mass being re-translocated (- is amount lost)
            prt.variables[i_var].net_alloc[i_pos] =
                prt.variables[i_var].net_alloc[i_pos] - retranslocated_mass

            # Update the state of the pool to reflect the mass lost
            prt.variables[i_var].val[i_pos] =
                prt.variables[i_var].val[i_pos] - (turnover_mass + retranslocated_mass)

            # Retranslocation is handled by the storage pool: add it there
            prt.variables[store_var_id].net_alloc[i_store_pos] =
                prt.variables[store_var_id].net_alloc[i_store_pos] + retranslocated_mass

            prt.variables[store_var_id].val[i_store_pos] =
                prt.variables[store_var_id].val[i_store_pos] + retranslocated_mass
        end
    end

    return nothing
end

# =====================================================================================

"""
    PRTMaintTurnover!(prt, ipft, icanlayer, is_drought)

Generic wrapper for maintenance turnover; dispatches to
[`MaintTurnoverSimpleRetranslocation!`](@ref).
"""
function PRTMaintTurnover!(prt::AbstractPRTVartypes, ipft::Integer, icanlayer::Integer,
                           is_drought::Bool)
    MaintTurnoverSimpleRetranslocation!(prt, ipft, icanlayer, is_drought)
    return nothing
end

# =====================================================================================

"""
    MaintTurnoverSimpleRetranslocation!(prt, ipft, icanlayer, is_drought)

Remove biomass from all applicable pools due to maintenance (background)
turnover. Assumed called daily. Branchfall reduces sapwood/structure/storage;
fine roots turn over with `root_long`; leaves only for evergreens (with the last
leaf-age "senescing" position generating the loss). A fraction is retranslocated
to storage; the rest enters the `turnover` diagnostic.
"""
function MaintTurnoverSimpleRetranslocation!(prt::AbstractPRTVartypes, ipft::Integer,
                                             icanlayer::Integer, is_drought::Bool)
    g = get_prt_global()

    if g.hyp_id <= 2
        i_store_pos = 1   # hypothesis 1&2 only have 1 storage pool
    else
        fates_endrun("You picked a hypothesis that has not defined how and where " *
                     "turnover re-absorption interacts with the storage pool.")
    end

    # Calculate the turnover rates (a temp for the actual turnover rate per organ)
    base_turnover = fill(un_initialized, num_organ_types)

    # All plants can have branch turnover, if branchfall is non-zero, which will
    # reduce sapwood, structure and storage.
    if prt_params.branch_long[ipft] > nearzero
        base_turnover[sapw_organ]   = years_per_day / prt_params.branch_long[ipft]
        base_turnover[struct_organ] = years_per_day / prt_params.branch_long[ipft]
        base_turnover[store_organ]  = years_per_day / prt_params.branch_long[ipft]
    else
        base_turnover[sapw_organ]   = 0.0
        base_turnover[struct_organ] = 0.0
        base_turnover[store_organ]  = 0.0
    end

    # All plants are allowed to have fine-root turnover if a non-zero life-span is selected
    if prt_params.root_long[ipft] > nearzero
        base_turnover[fnrt_organ] = years_per_day / prt_params.root_long[ipft]
    else
        base_turnover[fnrt_organ] = 0.0
    end

    if icanlayer == 1
        # The last index of the leaf longevity array contains the turnover
        # timescale for the senescent pool.
        aclass_sen_id = size(prt_params.leaf_long, 2)
        leaf_long = prt_params.leaf_long[ipft, aclass_sen_id]
    else
        aclass_sen_id = size(prt_params.leaf_long_ustory, 2)
        leaf_long = prt_params.leaf_long_ustory[ipft, aclass_sen_id]
    end

    # Only evergreens have maintenance turnover (must also change trimming logic
    # if we want to change this)
    if leaf_long > nearzero && prt_params.evergreen[ipft] == itrue
        if is_drought
            base_turnover[leaf_organ] = years_per_day /
                (leaf_long * prt_params.senleaf_long_fdrought[ipft])
        else
            base_turnover[leaf_organ] = years_per_day / leaf_long
        end
    else
        base_turnover[leaf_organ] = 0.0
    end

    base_turnover[repro_organ] = 0.0

    for i_var in 1:g.num_vars
        organ_id   = g.state_descriptor[i_var].organ_id
        element_id = g.state_descriptor[i_var].element_id

        # If this organ does not have a retranslocation rate then it is not valid
        # for turnover
        if prt_params.organ_param_id[organ_id] < 1
            retrans_frac = 0.0
        else
            if element_id == carbon12_element
                retrans_frac = 0.0
            elseif element_id == nitrogen_element
                retrans_frac = prt_params.turnover_nitr_retrans[ipft, prt_params.organ_param_id[organ_id]]
            elseif element_id == phosphorus_element
                retrans_frac = prt_params.turnover_phos_retrans[ipft, prt_params.organ_param_id[organ_id]]
            else
                fates_endrun("Please add a new re-translocation clause to your " *
                             "organ x element combination. organ=$(organ_id) " *
                             "element=$(element_id)")
            end
        end

        if base_turnover[organ_id] < check_initialized
            fates_endrun("A maintenance turnover rate for the organ was not " *
                         "specified. organ=$(organ_id) element=$(element_id) " *
                         "base turnover rate=$(base_turnover[organ_id])")
        end

        if retrans_frac < 0.0 || retrans_frac > 1.0
            fates_endrun("Unacceptable retranslocation calculated. organ=$(organ_id) " *
                         "element=$(element_id) retranslocation fraction=$(retrans_frac)")
        end

        # Hypotheses 1 & 2 assume that the leaf pools are stratified by age.
        # We only generate turnover from the last (senescing) position.
        if organ_id == leaf_organ
            if g.hyp_id <= 2
                ipos_1 = g.state_descriptor[i_var].num_pos
            else
                fates_endrun("Unhandled Leaf maintenance turnover condition for " *
                             "PARTEH hypothesis id: $(g.hyp_id)")
            end
        else
            ipos_1 = 1
        end

        store_var_id = sp_organ_map_get(g, store_organ, element_id)

        for i_pos in ipos_1:g.state_descriptor[i_var].num_pos
            turnover_mass = (1.0 - retrans_frac) * base_turnover[organ_id] *
                            prt.variables[i_var].val[i_pos]

            # Remove turnover mass from the organ of interest
            prt.variables[i_var].turnover[i_pos] =
                prt.variables[i_var].turnover[i_pos] + turnover_mass

            prt.variables[i_var].val[i_pos] =
                prt.variables[i_var].val[i_pos] - turnover_mass

            # If any mass is re-absorbed, send it to storage
            retrans_mass = retrans_frac * base_turnover[organ_id] *
                           prt.variables[i_var].val[i_pos]

            prt.variables[i_var].net_alloc[i_pos] =
                prt.variables[i_var].net_alloc[i_pos] - retrans_mass

            prt.variables[i_var].val[i_pos] =
                prt.variables[i_var].val[i_pos] - retrans_mass

            prt.variables[store_var_id].net_alloc[i_store_pos] =
                prt.variables[store_var_id].net_alloc[i_store_pos] + retrans_mass

            prt.variables[store_var_id].val[i_store_pos] =
                prt.variables[store_var_id].val[i_store_pos] + retrans_mass
        end
    end

    return nothing
end

# =====================================================================================

"""
    PRTDamageRecoveryFluxes!(prt, organ_id, mass_0, mass, cc_mass)

Adjust the `net_alloc` diagnostic of `organ_id` (position 1) for the
damage-recovery copy of state between cohorts: remove the amount copied from the
old cohort (`cc_mass - mass_0`) and add the recovery change (`mass - mass_0`).
Note the Fortran indexes by `organ_id` directly into `prt%variables` (it relies
on the organ ordering coinciding with the variable ordering for this hypothesis).
"""
function PRTDamageRecoveryFluxes!(prt::AbstractPRTVartypes, organ_id::Integer,
                                  mass_0::Real, mass::Real, cc_mass::Real)
    icd = 1

    # Remove the amount that was copied from old cohort
    prt.variables[organ_id].net_alloc[icd] =
        prt.variables[organ_id].net_alloc[icd] - (cc_mass - mass_0)

    # Track the amount of mass being lost (+ is amount lost)
    prt.variables[organ_id].net_alloc[icd] =
        prt.variables[organ_id].net_alloc[icd] + (mass - mass_0)

    return nothing
end
