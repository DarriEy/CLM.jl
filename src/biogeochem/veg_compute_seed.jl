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
# species_from_string -- Parse a string representation into a CN species ID
# ---------------------------------------------------------------------------

"""
    species_from_string(s::AbstractString)::Int

Convert a string representation of a CN species into one of the
`CN_SPECIES_*` constants. The input should be lowercase, matching the
Fortran convention.

- `"c12"` -> `CN_SPECIES_C12`
- `"c13"` -> `CN_SPECIES_C13`
- `"c14"` -> `CN_SPECIES_C14`
- `"n"`   -> `CN_SPECIES_N`

Errors on an unrecognized string.

Ported from `species_from_string` in `CNSpeciesMod.F90`.
"""
function species_from_string(s::AbstractString)::Int
    if s == "c12"
        return CN_SPECIES_C12
    elseif s == "c13"
        return CN_SPECIES_C13
    elseif s == "c14"
        return CN_SPECIES_C14
    elseif s == "n"
        return CN_SPECIES_N
    else
        error("species_from_string: unknown species string: $s")
    end
end

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
                           leaf::Real,
                           leaf_storage::Real,
                           leaf_xfer::Real,
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
- `leaf_patch::Vector{<:Real}`: current leaf C or N content (g/m2)
- `leaf_storage_patch::Vector{<:Real}`: current leaf storage C or N (g/m2)
- `leaf_xfer_patch::Vector{<:Real}`: current leaf transfer C or N (g/m2)
- `compute_here_patch::Vector{Bool}`: whether to compute outputs per patch
- `ignore_current_state_patch::Vector{Bool}`: use default proportions if true
- `seed_leaf_patch::Vector{<:Real}`: (output) seed for leaf [g/m2]
- `seed_leaf_storage_patch::Vector{<:Real}`: (output) seed for leaf storage [g/m2]
- `seed_leaf_xfer_patch::Vector{<:Real}`: (output) seed for leaf transfer [g/m2]
- `seed_deadstem_patch::Vector{<:Real}`: (output) seed for deadstem [g/m2]
- `noveg_val::Int`: bare ground PFT index (Fortran 0-based, default 0)

Ported from `ComputeSeedAmounts` in `CNVegComputeSeedMod.F90`.
"""
# Per-patch seed kernel. Own-index writes (seed_*[p]); leaf_proportions and
# species_type_multiplier are inlined T-generically (their error() branches are
# unreachable for the valid species (1-4) / component the caller passes, so they
# are dropped — byte-identical on valid input). Only species N (==4) distinguishes
# leaf vs deadwood via leafcn/deadwdcn; C12/C13/C14 multipliers are component-free.
@kernel function _seed_amounts_kernel!(seed_leaf, seed_leaf_storage, seed_leaf_xfer,
        seed_deadstem, @Const(mask_soilp), @Const(compute_here), @Const(ignore_state),
        @Const(itype), @Const(leaf), @Const(leaf_storage), @Const(leaf_xfer),
        @Const(c3psn), @Const(leafcn), @Const(deadwdcn), @Const(evergreen), @Const(woody),
        species::Int, leafc_seed, deadstemc_seed, noveg_val::Int,
        c3_r2, c4_r2, c14ratio, begp::Int, endp::Int)
    p = @index(Global)
    T = eltype(seed_leaf)
    @inbounds if begp <= p <= endp && mask_soilp[p] && compute_here[p]
        pft_type = itype[p]
        ji = pft_type + 1

        # --- leaf_proportions (inlined) ---
        tot_leaf = leaf[p] + leaf_storage[p] + leaf_xfer[p]
        pleaf = zero(T); pstor = zero(T); pxfer = zero(T)
        if tot_leaf == zero(T) || ignore_state[p]
            if evergreen[ji] == one(T)
                pleaf = one(T)
            else
                pstor = one(T)
            end
        else
            pleaf = leaf[p] / tot_leaf
            pstor = leaf_storage[p] / tot_leaf
            pxfer = leaf_xfer[p] / tot_leaf
        end

        my_leaf_seed = zero(T)
        my_deadstem_seed = zero(T)
        if pft_type != noveg_val
            # species_type_multiplier(·, LEAF): C12=1, C13=C3/C4_R2, C14=ratio, N=1/leafcn
            mult_leaf = species == CN_SPECIES_C12 ? one(T) :
                        species == CN_SPECIES_C13 ? (c3psn[ji] == one(T) ? c3_r2 : c4_r2) :
                        species == CN_SPECIES_C14 ? c14ratio :
                                                    one(T) / leafcn[ji]
            my_leaf_seed = leafc_seed * mult_leaf
            if woody[ji] == one(T)
                mult_dead = species == CN_SPECIES_C12 ? one(T) :
                            species == CN_SPECIES_C13 ? (c3psn[ji] == one(T) ? c3_r2 : c4_r2) :
                            species == CN_SPECIES_C14 ? c14ratio :
                                                        one(T) / deadwdcn[ji]
                my_deadstem_seed = deadstemc_seed * mult_dead
            end
        end

        seed_leaf[p]         = my_leaf_seed * pleaf
        seed_leaf_storage[p] = my_leaf_seed * pstor
        seed_leaf_xfer[p]    = my_leaf_seed * pxfer
        seed_deadstem[p]     = my_deadstem_seed
    end
end

function compute_seed_amounts!(mask_soilp::AbstractVector{Bool},
                                bounds::UnitRange{Int},
                                patch::PatchData,
                                pftcon_data::PftconType;
                                species::Int,
                                leafc_seed::Real,
                                deadstemc_seed::Real,
                                leaf_patch::AbstractVector{<:Real},
                                leaf_storage_patch::AbstractVector{<:Real},
                                leaf_xfer_patch::AbstractVector{<:Real},
                                compute_here_patch::AbstractVector{Bool},
                                ignore_current_state_patch::AbstractVector{Bool},
                                seed_leaf_patch::AbstractVector{<:Real},
                                seed_leaf_storage_patch::AbstractVector{<:Real},
                                seed_leaf_xfer_patch::AbstractVector{<:Real},
                                seed_deadstem_patch::AbstractVector{<:Real},
                                noveg_val::Int = noveg)
    isempty(bounds) && return nothing
    T = eltype(seed_leaf_patch)
    # Copy the needed pftcon arrays onto the output's backend/precision (the const
    # global pftcon holds concrete host Vectors). On CPU this is a plain copy.
    onb(a) = seed_leaf_patch isa Array ? a :
             copyto!(similar(seed_leaf_patch, T, length(a)), a)

    _launch!(_seed_amounts_kernel!, seed_leaf_patch, seed_leaf_storage_patch,
             seed_leaf_xfer_patch, seed_deadstem_patch,
             mask_soilp, compute_here_patch, ignore_current_state_patch,
             patch.itype, leaf_patch, leaf_storage_patch, leaf_xfer_patch,
             onb(pftcon_data.c3psn), onb(pftcon_data.leafcn), onb(pftcon_data.deadwdcn),
             onb(pftcon_data.evergreen), onb(pftcon_data.woody),
             species, T(leafc_seed), T(deadstemc_seed), noveg_val,
             T(C3_R2), T(C4_R2), T(C14RATIO), first(bounds), last(bounds);
             ndrange = last(bounds))
    return nothing
end
