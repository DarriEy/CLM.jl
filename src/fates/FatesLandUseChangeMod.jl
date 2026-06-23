# FatesLandUseChangeMod.jl
# Julia port of FATES src/fates/biogeochem/FatesLandUseChangeMod.F90
#
# FATES land-use-change transitions. Ingests the LUH2 land-use transition rate
# information that the host model read from a dataset, aggregates the (12) LUH2
# state/transition types to the (5) FATES land use categories, and produces a
# transition matrix [m2/m2/day] that drives patch disturbance rates, plus the
# aggregated land-use state vector. Also defines the clearing-rules matrix
# (whether vegetation is cleared during a given transition) and the special
# spin-up→land-use initialization transitions/harvest.
#
# Deps: ed_site_type (area), EDPftvarcon, EDParamsMod, FatesConstantsMod
# (primaryland/secondaryland/rangeland/pastureland/cropland, n_landuse_cats,
# nearzero, itrue, ifalse, fates_unset_int, years_per_day), FatesInterfaceTypesMod
# (bc_in_type, hlm_use_luh, hlm_num_luh2_states, hlm_num_luh2_transitions,
# hlm_use_potentialveg), FatesUtilsMod (FindIndex).
#
# Standalone — NOT added to CLMInstances or any dual-copied struct.

# Module data: maximum number of luh2 types that map to a single fates lu type
const max_luh2_types_per_fates_lu_type = 5

# ---------------------------------------------------------------------------
# luh2_fates_lutype_map — mapping from luh2 state names to fates lu categories
# ---------------------------------------------------------------------------
# Mirrors the Fortran derived type `luh2_fates_lutype_map`. The 'urban' state
# maps to `fates_unset_int`, signalling that it is handled separately (and that
# the urban fraction is factored out of the remaining states/transitions).
struct luh2_fates_lutype_map
    state_names::Vector{String}
    landuse_categories::Vector{Int}
end

function luh2_fates_lutype_map()
    state_names = ["primf", "primn", "secdf", "secdn",
                   "pastr", "range", "urban",
                   "c3ann", "c4ann", "c3per", "c4per", "c3nfx"]
    landuse_categories = [primaryland, primaryland, secondaryland, secondaryland,
                          pastureland, rangeland, fates_unset_int,
                          cropland, cropland, cropland, cropland, cropland]
    return luh2_fates_lutype_map(state_names, landuse_categories)
end

# ---------------------------------------------------------------------------
# GetLUCategoryFromStateName (the type-bound `GetIndex` procedure)
# ---------------------------------------------------------------------------
"""
    GetLUCategoryFromStateName(this, state_name) -> landuse_category

Map a LUH2 state name (e.g. "primf") to the aggregated FATES land use category
index. Errors if the input name matches none of the known state names.
"""
function GetLUCategoryFromStateName(this::luh2_fates_lutype_map, state_name::AbstractString)
    index_statename = FindIndex(this.state_names, state_name)

    if index_statename == 0
        fates_endrun("The input state name from the HLM does not match the FATES landuse state name options: " *
                     String(state_name))
    end

    return this.landuse_categories[index_statename]
end

# ---------------------------------------------------------------------------
# GetLanduseTransitionRates
# ---------------------------------------------------------------------------
"""
    GetLanduseTransitionRates!(bc_in, min_allowed_landuse_fraction,
                               landuse_transition_matrix, landuse_vector_gt_min)

Ingest LUH2 transition rate data, aggregate to FATES land use types, and fill
`landuse_transition_matrix` [m2/m2/day] (donor x receiver). `landuse_vector_gt_min`
is updated in place. Transitions involving 'urban' or diagonal elements are
skipped; urban fraction is factored out of the remaining transitions.
"""
function GetLanduseTransitionRates!(bc_in::bc_in_type, min_allowed_landuse_fraction::Real,
                                    landuse_transition_matrix::AbstractMatrix{<:Real},
                                    landuse_vector_gt_min::AbstractVector{Bool})
    lumap = luh2_fates_lutype_map()

    # zero the transition matrix and the urban fraction
    fill!(landuse_transition_matrix, 0.0)
    urban_fraction = 0.0

    # if using potential veg only, keep all transitions equal to zero
    if hlm_use_potentialveg[] == ifalse

        # Check incoming LUH transitions for NaN
        temp_vector = copy(bc_in.hlm_luh_transitions)
        modified_flag = CheckLUHData!(temp_vector)
        if !modified_flag
            # identify urban fraction so it can be factored into the lu state output
            urban_fraction = bc_in.hlm_luh_states[FindIndex(bc_in.hlm_luh_state_names, "urban")]
        end

        for i_luh2_transitions in 1:hlm_num_luh2_transitions[]
            # transition names: xxxxx_to_yyyyy (donor 1:5, receiver 10:14)
            transition_name = bc_in.hlm_luh_transition_names[i_luh2_transitions]
            donor_name    = transition_name[1:5]
            receiver_name = transition_name[10:14]

            i_donor    = GetLUCategoryFromStateName(lumap, donor_name)
            i_receiver = GetLUCategoryFromStateName(lumap, receiver_name)

            # Avoid 'urban' transitions and diagonal elements
            if !(i_donor == fates_unset_int || i_receiver == fates_unset_int ||
                 i_donor == i_receiver)
                landuse_transition_matrix[i_donor, i_receiver] =
                    landuse_transition_matrix[i_donor, i_receiver] +
                    temp_vector[i_luh2_transitions] * years_per_day / (1.0 - urban_fraction)
            end
        end

        # zero all transitions where the receiving state is below the minimum allowed,
        # else if this is the first timestep where the minimum was exceeded, apply all
        # transitions from primary to this type and set the flag (skipping secondary).
        state_vector = GetLUHStatedata(bc_in)
        for i_lu in secondaryland:n_landuse_cats
            if state_vector[i_lu] <= min_allowed_landuse_fraction
                landuse_transition_matrix[:, i_lu] .= 0.0
            elseif (!landuse_vector_gt_min[i_lu]) && (i_lu != secondaryland)
                landuse_transition_matrix[:, i_lu] .= 0.0
                landuse_transition_matrix[primaryland, i_lu] = state_vector[i_lu]
                landuse_vector_gt_min[i_lu] = true
            end
        end
    end

    return nothing
end

# ---------------------------------------------------------------------------
# GetLanduseChangeRules
# ---------------------------------------------------------------------------
"""
    GetLanduseChangeRules() -> clearing_matrix

Define the ruleset for when vegetation is cleared during a transition from one
land use type to another (donor row, receiver column). `false` => do not clear;
`true` => clear. Uses ruleset 4 from Table 1 of Ma et al. (2020).
"""
function GetLanduseChangeRules()
    clearing_matrix = fill(false, n_landuse_cats, n_landuse_cats)

    # ruleset to apply from table 1 of Ma et al (2020)
    # https://doi.org/10.5194/gmd-13-3203-2020
    ruleset = 4

    if ruleset == 1
        clearing_matrix[:, cropland] .= true
        clearing_matrix[:, pastureland] .= true
        clearing_matrix[primaryland, rangeland] = true
        clearing_matrix[secondaryland, rangeland] = true
    elseif ruleset == 2
        clearing_matrix[:, cropland] .= true
        clearing_matrix[primaryland, pastureland] = true
        clearing_matrix[secondaryland, pastureland] = true
        clearing_matrix[primaryland, rangeland] = true
        clearing_matrix[secondaryland, rangeland] = true
    elseif ruleset == 3
        clearing_matrix[:, cropland] .= true
        clearing_matrix[:, pastureland] .= true
        clearing_matrix[:, rangeland] .= true
    elseif ruleset == 4
        clearing_matrix[:, cropland] .= true
        clearing_matrix[:, pastureland] .= true
        clearing_matrix[:, rangeland] .= false
    elseif ruleset == 5
        clearing_matrix[:, cropland] .= true
        clearing_matrix[:, pastureland] .= false
        clearing_matrix[:, rangeland] .= true
    elseif ruleset == 6
        clearing_matrix[:, cropland] .= true
        clearing_matrix[:, pastureland] .= false
        clearing_matrix[:, rangeland] .= false
    elseif ruleset == 7
        clearing_matrix[:, cropland] .= false
        clearing_matrix[:, pastureland] .= true
        clearing_matrix[:, rangeland] .= true
    elseif ruleset == 8
        clearing_matrix[:, cropland] .= false
        clearing_matrix[:, pastureland] .= true
        clearing_matrix[:, rangeland] .= false
    elseif ruleset == 9
        clearing_matrix[:, cropland] .= false
        clearing_matrix[:, pastureland] .= false
        clearing_matrix[:, rangeland] .= true
    else
        fates_endrun("unknown clearing ruleset: " * string(ruleset))
    end

    return clearing_matrix
end

# ---------------------------------------------------------------------------
# GetLUHStatedata
# ---------------------------------------------------------------------------
"""
    GetLUHStatedata(bc_in) -> state_vector

Aggregate the LUH2 state vector to the FATES land use categories [m2/m2]. Factors
out the urban fraction, normalizes the total to 1, and defaults to all primary
land when no data is present (or for potential-vegetation runs).
"""
function GetLUHStatedata(bc_in::bc_in_type)
    lumap = luh2_fates_lutype_map()

    state_vector = zeros(Float64, n_landuse_cats)
    urban_fraction = 0.0

    if hlm_use_potentialveg[] == itrue
        state_vector[primaryland] = 1.0
    else
        # Check incoming state vector for NaN
        temp_vector = copy(bc_in.hlm_luh_states)
        modified_flag = CheckLUHData!(temp_vector)
        if !modified_flag
            urban_fraction = bc_in.hlm_luh_states[FindIndex(bc_in.hlm_luh_state_names, "urban")]
        end

        # add up the states that correspond to each fates land use type
        for i_luh2_states in 1:hlm_num_luh2_states[]
            state_name = bc_in.hlm_luh_state_names[i_luh2_states]
            ii = GetLUCategoryFromStateName(lumap, state_name)

            # avoid 'urban' states (unset index)
            if ii != fates_unset_int
                state_vector[ii] = state_vector[ii] +
                    temp_vector[i_luh2_states] / (1.0 - urban_fraction)
            end
        end

        if sum(state_vector) > nearzero
            # ensure total area == 1, correct if not
            if abs(sum(state_vector) - 1.0) > nearzero
                state_vector .= state_vector ./ sum(state_vector)
            end
        else
            state_vector[primaryland] = 1.0
        end
    end

    return state_vector
end

# ---------------------------------------------------------------------------
# CheckLUHData
# ---------------------------------------------------------------------------
"""
    CheckLUHData!(luh_vector) -> modified_flag

Check the incoming LUH2 vector for NaN. If all NaN, zero the vector (and set
primary land to 1 if it is a state vector), returning `modified_flag = true`. A
partially-NaN vector indicates corrupt data and ends the run.
"""
function CheckLUHData!(luh_vector::AbstractVector{<:Real})
    modified_flag = false

    if all(isnan, luh_vector)
        luh_vector .= 0.0
        # if this is a state vector, set primary land to 1
        if length(luh_vector) == hlm_num_luh2_states[]
            luh_vector[primaryland] = 1.0
        end
        modified_flag = true
    elseif any(isnan, luh_vector)
        if any(x -> !isnan(x), luh_vector)
            fates_endrun("ERROR: land use vector has NaN")
        end
    end

    return modified_flag
end

# ---------------------------------------------------------------------------
# GetInitLanduseHarvestRate
# ---------------------------------------------------------------------------
"""
    GetInitLanduseHarvestRate(bc_in, min_allowed_landuse_fraction,
                              landuse_vector_gt_min) -> harvest_rate

Spin-up→land-use initialization: apply the land-use changes needed to reach the
target state vector in a single daily instance, for the harvest rate from
primary lands (primary→secondary). `landuse_vector_gt_min` is updated in place.
Returns the harvest rate [m2/m2/day].
"""
function GetInitLanduseHarvestRate(bc_in::bc_in_type, min_allowed_landuse_fraction::Real,
                                   landuse_vector_gt_min::AbstractVector{Bool})
    harvest_rate = 0.0
    state_vector = GetLUHStatedata(bc_in)

    if state_vector[secondaryland] > min_allowed_landuse_fraction
        harvest_rate = state_vector[secondaryland]
        landuse_vector_gt_min[secondaryland] = true
    end

    return harvest_rate
end

# ---------------------------------------------------------------------------
# GetInitLanduseTransitionRates
# ---------------------------------------------------------------------------
"""
    GetInitLanduseTransitionRates!(bc_in, min_allowed_landuse_fraction,
                                   landuse_transition_matrix, landuse_vector_gt_min)

Spin-up→land-use initialization: apply the land-use changes needed to reach the
target state vector in a single daily instance, for the transitions other than
harvest (primary → all categories except secondary). `landuse_transition_matrix`
and `landuse_vector_gt_min` are updated in place.
"""
function GetInitLanduseTransitionRates!(bc_in::bc_in_type, min_allowed_landuse_fraction::Real,
                                        landuse_transition_matrix::AbstractMatrix{<:Real},
                                        landuse_vector_gt_min::AbstractVector{Bool})
    fill!(landuse_transition_matrix, 0.0)

    state_vector = GetLUHStatedata(bc_in)

    for i in (secondaryland + 1):n_landuse_cats
        if state_vector[i] > min_allowed_landuse_fraction
            landuse_transition_matrix[primaryland, i] = state_vector[i]
            landuse_vector_gt_min[i] = true
        end
    end

    return nothing
end
