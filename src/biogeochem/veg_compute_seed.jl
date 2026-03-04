# ==========================================================================
# Ported from: src/biogeochem/CNVegComputeSeedMod.F90
# Compute seed amounts for new patch areas
#
# When patches increase in area, seed amounts are needed to initialize
# the carbon/nitrogen pools in the new area. This module computes those
# seed amounts for leaf (split into leaf itself, leaf storage, and leaf
# transfer) and deadstem, adjusting for C/N species (C12, C13, C14, N).
#
# Public functions:
#   compute_seed_amounts! -- Compute seed amounts for patches increasing in area
#
# Private helpers:
#   species_type_multiplier -- Convert gC to appropriate species amount
#   leaf_proportions        -- Compute leaf/storage/xfer proportions
# ==========================================================================

# --- Module-level constants ---
# Plant component identifiers (for species_type_multiplier dispatch)
const COMPONENT_LEAF     = 1
const COMPONENT_DEADWOOD = 2

# CN species identifiers (from CNSpeciesMod.F90)
const CN_SPECIES_C12 = 1
const CN_SPECIES_C13 = 2
const CN_SPECIES_C14 = 3
const CN_SPECIES_N   = 4

# ---------------------------------------------------------------------------
# species_type_multiplier -- Convert gC seed to appropriate species amount
# ---------------------------------------------------------------------------

"""
    species_type_multiplier(species, pft_type, component, pftcon_data)

Return a multiplier to convert a seed amount expressed in gC/m2 into the
appropriate value for the given C/N species.

- `species`: one of `CN_SPECIES_C12`, `CN_SPECIES_C13`, `CN_SPECIES_C14`, `CN_SPECIES_N`
- `pft_type`: Fortran 0-based PFT index (used for pftcon lookup via `pft_type + 1`)
- `component`: one of `COMPONENT_LEAF`, `COMPONENT_DEADWOOD`
- `pftcon_data`: a `PftconType` instance

Ported from `SpeciesTypeMultiplier` in `CNVegComputeSeedMod.F90`.
"""
function species_type_multiplier(species::Int, pft_type::Int, component::Int,
                                  pftcon_data::PftconType)
    ji = pft_type + 1  # Julia 1-based index into pftcon arrays

    if species == CN_SPECIES_C12
        return 1.0

    elseif species == CN_SPECIES_C13
        if pftcon_data.c3psn[ji] == 1.0
            return C3_R2
        else
            return C4_R2
        end

    elseif species == CN_SPECIES_C14
        # 14C state initialized assuming initial "modern" 14C of 1.e-12
        return C14RATIO

    elseif species == CN_SPECIES_N
        if component == COMPONENT_LEAF
            return 1.0 / pftcon_data.leafcn[ji]
        elseif component == COMPONENT_DEADWOOD
            return 1.0 / pftcon_data.deadwdcn[ji]
        else
            error("species_type_multiplier: unknown component: $component")
        end

    else
        error("species_type_multiplier: unknown species: $species")
    end
end

# ---------------------------------------------------------------------------
# leaf_proportions -- Compute leaf / storage / xfer proportions
# ---------------------------------------------------------------------------

"""
    leaf_proportions(ignore_current_state, pft_type, leaf, leaf_storage, leaf_xfer,
                     pftcon_data)

Compute proportions of total leaf pool allocated to leaf itself, storage,
and transfer. Returns `(pleaf, pstorage, pxfer)`.

If `ignore_current_state` is `true`, or if total leaf mass is zero, use
default proportions:
  - evergreen PFTs: all in leaf (pleaf = 1)
  - deciduous PFTs: all in storage (pstorage = 1)

Otherwise, proportions are based on current leaf/storage/xfer state.

- `pft_type`: Fortran 0-based PFT index
- `pftcon_data`: a `PftconType` instance

Ported from `LeafProportions` in `CNVegComputeSeedMod.F90`.
"""
function leaf_proportions(ignore_current_state::Bool,
                           pft_type::Int,
                           leaf::Float64,
                           leaf_storage::Float64,
                           leaf_xfer::Float64,
                           pftcon_data::PftconType)
    ji = pft_type + 1  # Julia 1-based index

    tot_leaf = leaf + leaf_storage + leaf_xfer
    pleaf    = 0.0
    pstorage = 0.0
    pxfer    = 0.0

    if tot_leaf == 0.0 || ignore_current_state
        if pftcon_data.evergreen[ji] == 1.0
            pleaf = 1.0
        else
            pstorage = 1.0
        end
    else
        pleaf    = leaf / tot_leaf
        pstorage = leaf_storage / tot_leaf
        pxfer    = leaf_xfer / tot_leaf
    end

    return (pleaf, pstorage, pxfer)
end

# ---------------------------------------------------------------------------
# compute_seed_amounts! -- Main public function
# ---------------------------------------------------------------------------

"""
    compute_seed_amounts!(mask_soilp, bounds, patch, pftcon_data;
        species, leafc_seed, deadstemc_seed,
        leaf_patch, leaf_storage_patch, leaf_xfer_patch,
        compute_here_patch, ignore_current_state_patch,
        seed_leaf_patch, seed_leaf_storage_patch, seed_leaf_xfer_patch,
        seed_deadstem_patch, noveg_val)

Compute seed amounts for patches that increase in area, for various
variables, for the given species (C12, C13, C14, or N).

Output variables (`seed_*_patch`) are only set for patches inside the mask
where `compute_here_patch[p]` is `true`; for other patches they retain
their original values.

Regardless of the species, `leafc_seed` and `deadstemc_seed` are specified
in terms of gC/m2; these amounts are converted to the appropriate species
amount internally.

# Arguments
- `mask_soilp::BitVector`: mask over soil patches (including inactive)
- `bounds::UnitRange{Int}`: patch index range
- `patch::PatchData`: patch data (for itype)
- `pftcon_data::PftconType`: PFT constants

# Keyword arguments
- `species::Int`: CN species identifier (CN_SPECIES_C12, etc.)
- `leafc_seed::Float64`: seed amount for leaf C (gC/m2)
- `deadstemc_seed::Float64`: seed amount for deadstem C (gC/m2)
- `leaf_patch::Vector{Float64}`: current leaf C or N content (g/m2)
- `leaf_storage_patch::Vector{Float64}`: current leaf storage C or N (g/m2)
- `leaf_xfer_patch::Vector{Float64}`: current leaf transfer C or N (g/m2)
- `compute_here_patch::Vector{Bool}`: whether to compute outputs per patch
- `ignore_current_state_patch::Vector{Bool}`: use default proportions if true
- `seed_leaf_patch::Vector{Float64}`: (output) seed for leaf [g/m2]
- `seed_leaf_storage_patch::Vector{Float64}`: (output) seed for leaf storage [g/m2]
- `seed_leaf_xfer_patch::Vector{Float64}`: (output) seed for leaf transfer [g/m2]
- `seed_deadstem_patch::Vector{Float64}`: (output) seed for deadstem [g/m2]
- `noveg_val::Int`: bare ground PFT index (Fortran 0-based, default 0)

Ported from `ComputeSeedAmounts` in `CNVegComputeSeedMod.F90`.
"""
function compute_seed_amounts!(mask_soilp::BitVector,
                                bounds::UnitRange{Int},
                                patch::PatchData,
                                pftcon_data::PftconType;
                                species::Int,
                                leafc_seed::Float64,
                                deadstemc_seed::Float64,
                                leaf_patch::Vector{Float64},
                                leaf_storage_patch::Vector{Float64},
                                leaf_xfer_patch::Vector{Float64},
                                compute_here_patch::Vector{Bool},
                                ignore_current_state_patch::Vector{Bool},
                                seed_leaf_patch::Vector{Float64},
                                seed_leaf_storage_patch::Vector{Float64},
                                seed_leaf_xfer_patch::Vector{Float64},
                                seed_deadstem_patch::Vector{Float64},
                                noveg_val::Int = noveg)
    for p in bounds
        mask_soilp[p] || continue
        compute_here_patch[p] || continue

        my_leaf_seed    = 0.0
        my_deadstem_seed = 0.0

        pft_type = patch.itype[p]

        pleaf, pstor, pxfer = leaf_proportions(
            ignore_current_state_patch[p],
            pft_type,
            leaf_patch[p],
            leaf_storage_patch[p],
            leaf_xfer_patch[p],
            pftcon_data)

        if pft_type != noveg_val
            my_leaf_seed = leafc_seed *
                species_type_multiplier(species, pft_type, COMPONENT_LEAF, pftcon_data)

            ji = pft_type + 1  # Julia 1-based index
            if pftcon_data.woody[ji] == 1.0
                my_deadstem_seed = deadstemc_seed *
                    species_type_multiplier(species, pft_type, COMPONENT_DEADWOOD, pftcon_data)
            end
        end

        seed_leaf_patch[p]         = my_leaf_seed * pleaf
        seed_leaf_storage_patch[p] = my_leaf_seed * pstor
        seed_leaf_xfer_patch[p]    = my_leaf_seed * pxfer
        seed_deadstem_patch[p]     = my_deadstem_seed
    end

    return nothing
end
