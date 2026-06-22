# FatesIOVariableKindMod.jl
# Julia port of FATES src/fates/main/FatesIOVariableKindMod.F90
#
# IO variable-kind descriptors: string constants naming each variable kind,
# history-variable group indices, and a type describing one IO variable kind
# (name, dims, active flag). `Init!`/`set_active!`/`is_active` dispatch on the
# fates_io_variable_kind_type, so they coexist with IODimensions' Init!.
# Deps: FatesConstantsMod (fates_long_string_length, fates_unset_int),
#       FatesGlobals (fates_log, fates_endrun), FatesIODimensionsMod.

# ---------------------------------------------------------------------------
# Variable-kind name constants
# ---------------------------------------------------------------------------
const site_r8 = "SI_R8"
const site_int = "SI_INT"
const site_soil_r8 = "SI_SOIL_R8"
const site_size_pft_r8 = "SI_SCPF_R8"
const site_size_r8 = "SI_SCLS_R8"
const site_coage_pft_r8 = "SI_CAPF_R8"
const site_coage_r8 = "SI_CACLS_R8"
const cohort_r8 = "CO_R8"
const cohort_int = "CO_INT"
const site_pft_r8 = "SI_PFT_R8"
const site_age_r8 = "SI_AGE_R8"
const site_height_r8 = "SI_HEIGHT_R8"
const site_fuel_r8 = "SI_FUEL_R8"
const site_cwdsc_r8 = "SI_CWDSC_R8"
const site_can_r8 = "SI_CAN_R8"
const site_cnlf_r8 = "SI_CNLF_R8"
const site_cdpf_r8 = "SI_CDPF_R8"
const site_cdsc_r8 = "SI_CDSC_R8"
const site_cdam_r8 = "SI_CDAM_R8"
const site_cnlfpft_r8 = "SI_CNLFPFT_R8"
const site_scag_r8 = "SI_SCAG_R8"
const site_scagpft_r8 = "SI_SCAGPFT_R8"
const site_agepft_r8 = "SI_AGEPFT_R8"
const site_agefuel_r8 = "SI_AGEFUEL_R8"
const site_clscpf_r8 = "SI_CLSCPF_R8"
const site_landuse_r8 = "SI_LANDUSE_R8"
const site_lulu_r8 = "SI_LULU_R8"
const site_lupft_r8 = "SI_LUPFT_R8"

# Element and multiplexed element dimensions
const site_elem_r8 = "SI_ELEM_R8"
const site_elpft_r8 = "SI_ELEMPFT_R8"
const site_elcwd_r8 = "SI_ELEMCWD_R8"
const site_elage_r8 = "SI_ELEMAGE_R8"

# ---------------------------------------------------------------------------
# History Variable Group indices (for zero-ing and initializing output vars).
# ---------------------------------------------------------------------------
# active when dimlevel(2)>0
const group_dyna_simple = 1
const group_nflx_simple = 7
# active when dimlevel(2)>1
const group_dyna_complx = 2
const group_nflx_complx = 8
# active when dimlevel(1)>0
const group_hifr_simple = 3
const group_hydr_simple = 5
# active when dimlevel(1)>1
const group_hifr_complx = 4
const group_hydr_complx = 6

# ---------------------------------------------------------------------------
# fates_io_variable_kind_type — one IO variable kind descriptor.
# ---------------------------------------------------------------------------
Base.@kwdef mutable struct fates_io_variable_kind_type
    name::String = ""          # string labelling this IO type
    ndims::Int = 0             # number of dimensions
    dimsize::Vector{Int} = Int[]  # size of each dimension
    active_::Bool = false
    dim1_index::Int = fates_unset_int
    dim2_index::Int = fates_unset_int
end

# =====================================================================================

"""
    Init!(this::fates_io_variable_kind_type, name, num_dims)

Initialize an IO variable-kind descriptor, allocating `dimsize` (filled with
`fates_unset_int`) and marking it inactive.
"""
function Init!(this::fates_io_variable_kind_type, name::AbstractString, num_dims::Integer)
    this.name = strip(name)
    this.ndims = num_dims
    this.dimsize = fill(fates_unset_int, num_dims)
    this.active_ = false
    this.dim1_index = fates_unset_int
    this.dim2_index = fates_unset_int
    return nothing
end

# =====================================================================================

set_active!(this::fates_io_variable_kind_type) = (this.active_ = true; nothing)

is_active(this::fates_io_variable_kind_type) = this.active_

# =====================================================================================

"""
    iotype_index(iotype_name, num_dim_kinds, dim_kinds) -> dk_index

Return the (1-based) index in `dim_kinds` of the kind whose name matches
`iotype_name`. Aborts if no match is found.
"""
function iotype_index(iotype_name::AbstractString, num_dim_kinds::Integer,
                      dim_kinds::AbstractVector{fates_io_variable_kind_type})
    for dk_index in 1:num_dim_kinds
        if strip(iotype_name) == strip(dim_kinds[dk_index].name)
            return dk_index
        end
    end
    @warn "An IOTYPE THAT DOESNT EXIST WAS SPECIFIED"
    fates_endrun("iotype_index: iotype '$iotype_name' does not exist")
    return 0
end
