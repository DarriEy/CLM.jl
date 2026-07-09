# ==========================================================================
# Ported from: src/biogeochem/CNVegCarbonStateType.F90
# Vegetation carbon state data type allocation and initialization
# ==========================================================================

# --- Module-level constants ---
const SPINUP_FACTOR_DEADWOOD_DEFAULT = 1.0   # spinup factor used for this simulation
const SPINUP_FACTOR_AD = 10.0                # spinup factor when in Accelerated Decomposition mode
const INITIAL_VEGC = 20.0                    # initial vegetation carbon for leafc/frootc and storage (gC/m2)

"""
    CNVegCarbonStateData

Vegetation carbon state data structure. Holds carbon pools at patch, column,
and gridcell levels for C12 (and optionally C13/C14 isotopes).

Ported from `cnveg_carbonstate_type` in `CNVegCarbonStateType.F90`.
"""
Base.@kwdef mutable struct CNVegCarbonStateData{FT<:Real,
                                    V<:AbstractVector{FT},
                                    M<:AbstractMatrix{FT}}
    species::Int = 0  # carbon species: 1=C12, 2=C13, 3=C14

    # --- Patch-level 2D (patch x nrepr) ---
    reproductivec_patch               ::M = Matrix{Float64}(undef, 0, 0)  # (gC/m2) reproductive (grain) C
    reproductivec_storage_patch       ::M = Matrix{Float64}(undef, 0, 0)  # (gC/m2) reproductive C storage
    reproductivec_xfer_patch          ::M = Matrix{Float64}(undef, 0, 0)  # (gC/m2) reproductive C transfer

    # --- Patch-level 1D core pools ---
    leafc_patch                       ::V = Float64[]  # (gC/m2) leaf C
    leafc_storage_patch               ::V = Float64[]  # (gC/m2) leaf C storage
    leafc_xfer_patch                  ::V = Float64[]  # (gC/m2) leaf C transfer
    leafc_storage_xfer_acc_patch      ::V = Float64[]  # (gC/m2) accumulated leaf C transfer
    storage_cdemand_patch             ::V = Float64[]  # (gC/m2) C use from storage pool
    frootc_patch                      ::V = Float64[]  # (gC/m2) fine root C
    frootc_storage_patch              ::V = Float64[]  # (gC/m2) fine root C storage
    frootc_xfer_patch                 ::V = Float64[]  # (gC/m2) fine root C transfer
    livestemc_patch                   ::V = Float64[]  # (gC/m2) live stem C
    livestemc_storage_patch           ::V = Float64[]  # (gC/m2) live stem C storage
    livestemc_xfer_patch              ::V = Float64[]  # (gC/m2) live stem C transfer
    deadstemc_patch                   ::V = Float64[]  # (gC/m2) dead stem C
    deadstemc_storage_patch           ::V = Float64[]  # (gC/m2) dead stem C storage
    deadstemc_xfer_patch              ::V = Float64[]  # (gC/m2) dead stem C transfer
    livecrootc_patch                  ::V = Float64[]  # (gC/m2) live coarse root C
    livecrootc_storage_patch          ::V = Float64[]  # (gC/m2) live coarse root C storage
    livecrootc_xfer_patch             ::V = Float64[]  # (gC/m2) live coarse root C transfer
    deadcrootc_patch                  ::V = Float64[]  # (gC/m2) dead coarse root C
    deadcrootc_storage_patch          ::V = Float64[]  # (gC/m2) dead coarse root C storage
    deadcrootc_xfer_patch             ::V = Float64[]  # (gC/m2) dead coarse root C transfer
    gresp_storage_patch               ::V = Float64[]  # (gC/m2) growth respiration storage
    gresp_xfer_patch                  ::V = Float64[]  # (gC/m2) growth respiration transfer
    cpool_patch                       ::V = Float64[]  # (gC/m2) temporary photosynthate C pool
    xsmrpool_patch                    ::V = Float64[]  # (gC/m2) abstract C pool for excess MR demand
    xsmrpool_loss_patch               ::V = Float64[]  # (gC/m2) abstract C pool excess MR demand loss
    ctrunc_patch                      ::V = Float64[]  # (gC/m2) patch-level sink for C truncation
    woodc_patch                       ::V = Float64[]  # (gC/m2) wood C
    leafcmax_patch                    ::V = Float64[]  # (gC/m2) annual max leaf C
    cropseedc_deficit_patch           ::V = Float64[]  # (gC/m2) pool for seeding new crop growth (negative)

    # --- Matrix CN capacity fields (patch-level) ---
    matrix_cap_leafc_patch              ::V = Float64[]  # (gC/m2) capacity of leaf C
    matrix_cap_leafc_storage_patch      ::V = Float64[]  # (gC/m2) capacity of leaf C storage
    matrix_cap_leafc_xfer_patch         ::V = Float64[]  # (gC/m2) capacity of leaf C transfer
    matrix_cap_frootc_patch             ::V = Float64[]  # (gC/m2) capacity of fine root C
    matrix_cap_frootc_storage_patch     ::V = Float64[]  # (gC/m2) capacity of fine root C storage
    matrix_cap_frootc_xfer_patch        ::V = Float64[]  # (gC/m2) capacity of fine root C transfer
    matrix_cap_livestemc_patch          ::V = Float64[]  # (gC/m2) capacity of live stem C
    matrix_cap_livestemc_storage_patch  ::V = Float64[]  # (gC/m2) capacity of live stem C storage
    matrix_cap_livestemc_xfer_patch     ::V = Float64[]  # (gC/m2) capacity of live stem C transfer
    matrix_cap_deadstemc_patch          ::V = Float64[]  # (gC/m2) capacity of dead stem C
    matrix_cap_deadstemc_storage_patch  ::V = Float64[]  # (gC/m2) capacity of dead stem C storage
    matrix_cap_deadstemc_xfer_patch     ::V = Float64[]  # (gC/m2) capacity of dead stem C transfer
    matrix_cap_livecrootc_patch         ::V = Float64[]  # (gC/m2) capacity of live coarse root C
    matrix_cap_livecrootc_storage_patch ::V = Float64[]  # (gC/m2) capacity of live coarse root C storage
    matrix_cap_livecrootc_xfer_patch    ::V = Float64[]  # (gC/m2) capacity of live coarse root C transfer
    matrix_cap_deadcrootc_patch         ::V = Float64[]  # (gC/m2) capacity of dead coarse root C
    matrix_cap_deadcrootc_storage_patch ::V = Float64[]  # (gC/m2) capacity of dead coarse root C storage
    matrix_cap_deadcrootc_xfer_patch    ::V = Float64[]  # (gC/m2) capacity of dead coarse root C transfer
    matrix_cap_reproc_patch             ::V = Float64[]  # (gC/m2) capacity of grain C
    matrix_cap_reproc_storage_patch     ::V = Float64[]  # (gC/m2) capacity of grain storage C
    matrix_cap_reproc_xfer_patch        ::V = Float64[]  # (gC/m2) capacity of grain transfer C

    # --- Initial pool size of year for matrix (SASU) ---
    leafc0_patch                      ::V = Float64[]  # (gC/m2) initial leaf C for SASU
    leafc0_storage_patch              ::V = Float64[]  # (gC/m2) initial leaf C storage for SASU
    leafc0_xfer_patch                 ::V = Float64[]  # (gC/m2) initial leaf C transfer for SASU
    frootc0_patch                     ::V = Float64[]  # (gC/m2) initial fine root C for SASU
    frootc0_storage_patch             ::V = Float64[]  # (gC/m2) initial fine root C storage for SASU
    frootc0_xfer_patch                ::V = Float64[]  # (gC/m2) initial fine root C transfer for SASU
    livestemc0_patch                  ::V = Float64[]  # (gC/m2) initial live stem C for SASU
    livestemc0_storage_patch          ::V = Float64[]  # (gC/m2) initial live stem C storage for SASU
    livestemc0_xfer_patch             ::V = Float64[]  # (gC/m2) initial live stem C transfer for SASU
    deadstemc0_patch                  ::V = Float64[]  # (gC/m2) initial dead stem C for SASU
    deadstemc0_storage_patch          ::V = Float64[]  # (gC/m2) initial dead stem C storage for SASU
    deadstemc0_xfer_patch             ::V = Float64[]  # (gC/m2) initial dead stem C transfer for SASU
    livecrootc0_patch                 ::V = Float64[]  # (gC/m2) initial live coarse root C for SASU
    livecrootc0_storage_patch         ::V = Float64[]  # (gC/m2) initial live coarse root C storage for SASU
    livecrootc0_xfer_patch            ::V = Float64[]  # (gC/m2) initial live coarse root C transfer for SASU
    deadcrootc0_patch                 ::V = Float64[]  # (gC/m2) initial dead coarse root C for SASU
    deadcrootc0_storage_patch         ::V = Float64[]  # (gC/m2) initial dead coarse root C storage for SASU
    deadcrootc0_xfer_patch            ::V = Float64[]  # (gC/m2) initial dead coarse root C transfer for SASU
    reproc0_patch                     ::V = Float64[]  # (gC/m2) initial grain C for SASU
    reproc0_storage_patch             ::V = Float64[]  # (gC/m2) initial grain C storage for SASU
    reproc0_xfer_patch                ::V = Float64[]  # (gC/m2) initial grain C transfer for SASU

    # --- SASU save fields ---
    leafc_SASUsave_patch              ::V = Float64[]
    leafc_storage_SASUsave_patch      ::V = Float64[]
    leafc_xfer_SASUsave_patch         ::V = Float64[]
    frootc_SASUsave_patch             ::V = Float64[]
    frootc_storage_SASUsave_patch     ::V = Float64[]
    frootc_xfer_SASUsave_patch        ::V = Float64[]
    livestemc_SASUsave_patch          ::V = Float64[]
    livestemc_storage_SASUsave_patch  ::V = Float64[]
    livestemc_xfer_SASUsave_patch     ::V = Float64[]
    deadstemc_SASUsave_patch          ::V = Float64[]
    deadstemc_storage_SASUsave_patch  ::V = Float64[]
    deadstemc_xfer_SASUsave_patch     ::V = Float64[]
    livecrootc_SASUsave_patch         ::V = Float64[]
    livecrootc_storage_SASUsave_patch ::V = Float64[]
    livecrootc_xfer_SASUsave_patch    ::V = Float64[]
    deadcrootc_SASUsave_patch         ::V = Float64[]
    deadcrootc_storage_SASUsave_patch ::V = Float64[]
    deadcrootc_xfer_SASUsave_patch    ::V = Float64[]
    grainc_SASUsave_patch             ::V = Float64[]
    grainc_storage_SASUsave_patch     ::V = Float64[]

    # --- Matrix accumulation variables (annual, patch-level) ---
    matrix_calloc_leaf_acc_patch        ::V = Float64[]
    matrix_calloc_leafst_acc_patch      ::V = Float64[]
    matrix_calloc_froot_acc_patch       ::V = Float64[]
    matrix_calloc_frootst_acc_patch     ::V = Float64[]
    matrix_calloc_livestem_acc_patch    ::V = Float64[]
    matrix_calloc_livestemst_acc_patch  ::V = Float64[]
    matrix_calloc_deadstem_acc_patch    ::V = Float64[]
    matrix_calloc_deadstemst_acc_patch  ::V = Float64[]
    matrix_calloc_livecroot_acc_patch   ::V = Float64[]
    matrix_calloc_livecrootst_acc_patch ::V = Float64[]
    matrix_calloc_deadcroot_acc_patch   ::V = Float64[]
    matrix_calloc_deadcrootst_acc_patch ::V = Float64[]
    matrix_calloc_grain_acc_patch       ::V = Float64[]
    matrix_calloc_grainst_acc_patch     ::V = Float64[]

    # --- Matrix C transfer accumulation (annual, patch-level) ---
    matrix_ctransfer_leafst_to_leafxf_acc_patch           ::V = Float64[]
    matrix_ctransfer_leafxf_to_leaf_acc_patch             ::V = Float64[]
    matrix_ctransfer_frootst_to_frootxf_acc_patch         ::V = Float64[]
    matrix_ctransfer_frootxf_to_froot_acc_patch           ::V = Float64[]
    matrix_ctransfer_livestemst_to_livestemxf_acc_patch    ::V = Float64[]
    matrix_ctransfer_livestemxf_to_livestem_acc_patch      ::V = Float64[]
    matrix_ctransfer_deadstemst_to_deadstemxf_acc_patch    ::V = Float64[]
    matrix_ctransfer_deadstemxf_to_deadstem_acc_patch      ::V = Float64[]
    matrix_ctransfer_livecrootst_to_livecrootxf_acc_patch ::V = Float64[]
    matrix_ctransfer_livecrootxf_to_livecroot_acc_patch   ::V = Float64[]
    matrix_ctransfer_deadcrootst_to_deadcrootxf_acc_patch ::V = Float64[]
    matrix_ctransfer_deadcrootxf_to_deadcroot_acc_patch   ::V = Float64[]
    matrix_ctransfer_grainst_to_grainxf_acc_patch         ::V = Float64[]
    matrix_ctransfer_grainxf_to_grain_acc_patch           ::V = Float64[]
    matrix_ctransfer_livestem_to_deadstem_acc_patch        ::V = Float64[]
    matrix_ctransfer_livecroot_to_deadcroot_acc_patch      ::V = Float64[]

    # --- Matrix C turnover accumulation (annual, patch-level) ---
    matrix_cturnover_leaf_acc_patch             ::V = Float64[]
    matrix_cturnover_leafst_acc_patch           ::V = Float64[]
    matrix_cturnover_leafxf_acc_patch           ::V = Float64[]
    matrix_cturnover_froot_acc_patch            ::V = Float64[]
    matrix_cturnover_frootst_acc_patch          ::V = Float64[]
    matrix_cturnover_frootxf_acc_patch          ::V = Float64[]
    matrix_cturnover_livestem_acc_patch         ::V = Float64[]
    matrix_cturnover_livestemst_acc_patch       ::V = Float64[]
    matrix_cturnover_livestemxf_acc_patch       ::V = Float64[]
    matrix_cturnover_deadstem_acc_patch         ::V = Float64[]
    matrix_cturnover_deadstemst_acc_patch       ::V = Float64[]
    matrix_cturnover_deadstemxf_acc_patch       ::V = Float64[]
    matrix_cturnover_livecroot_acc_patch        ::V = Float64[]
    matrix_cturnover_livecrootst_acc_patch      ::V = Float64[]
    matrix_cturnover_livecrootxf_acc_patch      ::V = Float64[]
    matrix_cturnover_deadcroot_acc_patch        ::V = Float64[]
    matrix_cturnover_deadcrootst_acc_patch      ::V = Float64[]
    matrix_cturnover_deadcrootxf_acc_patch      ::V = Float64[]
    matrix_cturnover_grain_acc_patch            ::V = Float64[]
    matrix_cturnover_grainst_acc_patch          ::V = Float64[]
    matrix_cturnover_grainxf_acc_patch          ::V = Float64[]

    # --- Column-level fields ---
    rootc_col                         ::V = Float64[]  # (gC/m2) root carbon at column level (fire)
    leafc_col                         ::V = Float64[]  # (gC/m2) column-level leafc (fire)
    deadstemc_col                     ::V = Float64[]  # (gC/m2) column-level deadstemc (fire)
    fuelc_col                         ::V = Float64[]  # fuel load outside cropland
    fuelc_crop_col                    ::V = Float64[]  # fuel load for cropland

    # --- Gridcell-level ---
    seedc_grc                         ::V = Float64[]  # (gC/m2) pool for seeding new PFTs via dynamic landcover

    # --- Summary (diagnostic) state variables ---
    dispvegc_patch                    ::V = Float64[]  # (gC/m2) displayed veg carbon, excluding storage and cpool
    storvegc_patch                    ::V = Float64[]  # (gC/m2) stored vegetation carbon, excluding cpool
    totc_patch                        ::V = Float64[]  # (gC/m2) total patch-level carbon, including cpool
    totvegc_patch                     ::V = Float64[]  # (gC/m2) total vegetation carbon, excluding cpool
    totvegc_col                       ::V = Float64[]  # (gC/m2) total vegetation carbon at column level
    totc_p2c_col                      ::V = Float64[]  # (gC/m2) totc_patch averaged to column

    # --- Private flag ---
    dribble_crophrv_xsmrpool_2atm     ::Bool = false
end

CNVegCarbonStateData{FT}(; kwargs...) where {FT<:Real} =
    CNVegCarbonStateData{FT, Vector{FT}, Matrix{FT}}(; kwargs...)
Adapt.@adapt_structure CNVegCarbonStateData


"""
    cnveg_carbon_state_init!(cs, np, nc, ng; use_matrixcn=false, nrepr=NREPR)

Allocate and initialize all fields of a `CNVegCarbonStateData` instance for
`np` patches, `nc` columns, and `ng` gridcells.
Matches the Fortran `InitAllocate`.
"""
function cnveg_carbon_state_init!(cs::CNVegCarbonStateData, np::Int, nc::Int, ng::Int;
                                   use_matrixcn::Bool=false, nrepr::Int=NREPR)
    # --- Patch-level 2D (patch x nrepr) ---
    cs.reproductivec_patch          = fill(NaN, np, nrepr)
    cs.reproductivec_storage_patch  = fill(NaN, np, nrepr)
    cs.reproductivec_xfer_patch     = fill(NaN, np, nrepr)

    # --- Patch-level 1D core pools ---
    cs.leafc_patch                  = fill(NaN, np)
    cs.leafc_storage_patch          = fill(NaN, np)
    cs.leafc_xfer_patch             = fill(NaN, np)
    cs.leafc_storage_xfer_acc_patch = fill(NaN, np)
    cs.storage_cdemand_patch        = fill(NaN, np)
    cs.frootc_patch                 = fill(NaN, np)
    cs.frootc_storage_patch         = fill(NaN, np)
    cs.frootc_xfer_patch            = fill(NaN, np)
    cs.livestemc_patch              = fill(NaN, np)
    cs.livestemc_storage_patch      = fill(NaN, np)
    cs.livestemc_xfer_patch         = fill(NaN, np)
    cs.deadstemc_patch              = fill(NaN, np)
    cs.deadstemc_storage_patch      = fill(NaN, np)
    cs.deadstemc_xfer_patch         = fill(NaN, np)
    cs.livecrootc_patch             = fill(NaN, np)
    cs.livecrootc_storage_patch     = fill(NaN, np)
    cs.livecrootc_xfer_patch        = fill(NaN, np)
    cs.deadcrootc_patch             = fill(NaN, np)
    cs.deadcrootc_storage_patch     = fill(NaN, np)
    cs.deadcrootc_xfer_patch        = fill(NaN, np)
    cs.gresp_storage_patch          = fill(NaN, np)
    cs.gresp_xfer_patch             = fill(NaN, np)
    cs.cpool_patch                  = fill(NaN, np)
    cs.xsmrpool_patch               = fill(NaN, np)
    cs.xsmrpool_loss_patch          = fill(NaN, np)
    cs.ctrunc_patch                 = fill(NaN, np)
    cs.woodc_patch                  = fill(NaN, np)
    cs.leafcmax_patch               = fill(NaN, np)
    cs.cropseedc_deficit_patch      = fill(NaN, np)

    # --- Summary fields ---
    cs.dispvegc_patch               = fill(NaN, np)
    cs.storvegc_patch               = fill(NaN, np)
    cs.totc_patch                   = fill(NaN, np)
    cs.totvegc_patch                = fill(NaN, np)

    # --- Column-level ---
    cs.rootc_col                    = fill(NaN, nc)
    cs.leafc_col                    = fill(NaN, nc)
    cs.deadstemc_col                = fill(NaN, nc)
    cs.fuelc_col                    = fill(NaN, nc)
    cs.fuelc_crop_col               = fill(NaN, nc)
    cs.totvegc_col                  = fill(NaN, nc)
    cs.totc_p2c_col                 = fill(NaN, nc)

    # --- Gridcell-level ---
    cs.seedc_grc                    = fill(NaN, ng)

    # --- Matrix CN fields (conditional) ---
    if use_matrixcn
        cs.matrix_cap_leafc_patch              = fill(NaN, np)
        cs.matrix_cap_leafc_storage_patch      = fill(NaN, np)
        cs.matrix_cap_leafc_xfer_patch         = fill(NaN, np)
        cs.matrix_cap_frootc_patch             = fill(NaN, np)
        cs.matrix_cap_frootc_storage_patch     = fill(NaN, np)
        cs.matrix_cap_frootc_xfer_patch        = fill(NaN, np)
        cs.matrix_cap_livestemc_patch          = fill(NaN, np)
        cs.matrix_cap_livestemc_storage_patch  = fill(NaN, np)
        cs.matrix_cap_livestemc_xfer_patch     = fill(NaN, np)
        cs.matrix_cap_deadstemc_patch          = fill(NaN, np)
        cs.matrix_cap_deadstemc_storage_patch  = fill(NaN, np)
        cs.matrix_cap_deadstemc_xfer_patch     = fill(NaN, np)
        cs.matrix_cap_livecrootc_patch         = fill(NaN, np)
        cs.matrix_cap_livecrootc_storage_patch = fill(NaN, np)
        cs.matrix_cap_livecrootc_xfer_patch    = fill(NaN, np)
        cs.matrix_cap_deadcrootc_patch         = fill(NaN, np)
        cs.matrix_cap_deadcrootc_storage_patch = fill(NaN, np)
        cs.matrix_cap_deadcrootc_xfer_patch    = fill(NaN, np)
        cs.matrix_cap_reproc_patch             = fill(NaN, np)
        cs.matrix_cap_reproc_storage_patch     = fill(NaN, np)
        cs.matrix_cap_reproc_xfer_patch        = fill(NaN, np)

        # Initial pool size for SASU
        cs.leafc0_patch                = fill(NaN, np)
        cs.leafc0_storage_patch        = fill(NaN, np)
        cs.leafc0_xfer_patch           = fill(NaN, np)
        cs.frootc0_patch               = fill(NaN, np)
        cs.frootc0_storage_patch       = fill(NaN, np)
        cs.frootc0_xfer_patch          = fill(NaN, np)
        cs.livestemc0_patch            = fill(NaN, np)
        cs.livestemc0_storage_patch    = fill(NaN, np)
        cs.livestemc0_xfer_patch       = fill(NaN, np)
        cs.deadstemc0_patch            = fill(NaN, np)
        cs.deadstemc0_storage_patch    = fill(NaN, np)
        cs.deadstemc0_xfer_patch       = fill(NaN, np)
        cs.livecrootc0_patch           = fill(NaN, np)
        cs.livecrootc0_storage_patch   = fill(NaN, np)
        cs.livecrootc0_xfer_patch      = fill(NaN, np)
        cs.deadcrootc0_patch           = fill(NaN, np)
        cs.deadcrootc0_storage_patch   = fill(NaN, np)
        cs.deadcrootc0_xfer_patch      = fill(NaN, np)
        cs.reproc0_patch               = fill(NaN, np)
        cs.reproc0_storage_patch       = fill(NaN, np)
        cs.reproc0_xfer_patch          = fill(NaN, np)

        # SASU save fields
        cs.leafc_SASUsave_patch              = fill(NaN, np)
        cs.leafc_storage_SASUsave_patch      = fill(NaN, np)
        cs.leafc_xfer_SASUsave_patch         = fill(NaN, np)
        cs.frootc_SASUsave_patch             = fill(NaN, np)
        cs.frootc_storage_SASUsave_patch     = fill(NaN, np)
        cs.frootc_xfer_SASUsave_patch        = fill(NaN, np)
        cs.livestemc_SASUsave_patch          = fill(NaN, np)
        cs.livestemc_storage_SASUsave_patch  = fill(NaN, np)
        cs.livestemc_xfer_SASUsave_patch     = fill(NaN, np)
        cs.deadstemc_SASUsave_patch          = fill(NaN, np)
        cs.deadstemc_storage_SASUsave_patch  = fill(NaN, np)
        cs.deadstemc_xfer_SASUsave_patch     = fill(NaN, np)
        cs.livecrootc_SASUsave_patch         = fill(NaN, np)
        cs.livecrootc_storage_SASUsave_patch = fill(NaN, np)
        cs.livecrootc_xfer_SASUsave_patch    = fill(NaN, np)
        cs.deadcrootc_SASUsave_patch         = fill(NaN, np)
        cs.deadcrootc_storage_SASUsave_patch = fill(NaN, np)
        cs.deadcrootc_xfer_SASUsave_patch    = fill(NaN, np)
        cs.grainc_SASUsave_patch             = fill(NaN, np)
        cs.grainc_storage_SASUsave_patch     = fill(NaN, np)

        # Matrix alloc accumulation
        cs.matrix_calloc_leaf_acc_patch        = fill(NaN, np)
        cs.matrix_calloc_leafst_acc_patch      = fill(NaN, np)
        cs.matrix_calloc_froot_acc_patch       = fill(NaN, np)
        cs.matrix_calloc_frootst_acc_patch     = fill(NaN, np)
        cs.matrix_calloc_livestem_acc_patch    = fill(NaN, np)
        cs.matrix_calloc_livestemst_acc_patch  = fill(NaN, np)
        cs.matrix_calloc_deadstem_acc_patch    = fill(NaN, np)
        cs.matrix_calloc_deadstemst_acc_patch  = fill(NaN, np)
        cs.matrix_calloc_livecroot_acc_patch   = fill(NaN, np)
        cs.matrix_calloc_livecrootst_acc_patch = fill(NaN, np)
        cs.matrix_calloc_deadcroot_acc_patch   = fill(NaN, np)
        cs.matrix_calloc_deadcrootst_acc_patch = fill(NaN, np)
        cs.matrix_calloc_grain_acc_patch       = fill(NaN, np)
        cs.matrix_calloc_grainst_acc_patch     = fill(NaN, np)

        # Matrix C transfer accumulation
        cs.matrix_ctransfer_leafst_to_leafxf_acc_patch           = fill(NaN, np)
        cs.matrix_ctransfer_leafxf_to_leaf_acc_patch             = fill(NaN, np)
        cs.matrix_ctransfer_frootst_to_frootxf_acc_patch         = fill(NaN, np)
        cs.matrix_ctransfer_frootxf_to_froot_acc_patch           = fill(NaN, np)
        cs.matrix_ctransfer_livestemst_to_livestemxf_acc_patch   = fill(NaN, np)
        cs.matrix_ctransfer_livestemxf_to_livestem_acc_patch     = fill(NaN, np)
        cs.matrix_ctransfer_deadstemst_to_deadstemxf_acc_patch   = fill(NaN, np)
        cs.matrix_ctransfer_deadstemxf_to_deadstem_acc_patch     = fill(NaN, np)
        cs.matrix_ctransfer_livecrootst_to_livecrootxf_acc_patch = fill(NaN, np)
        cs.matrix_ctransfer_livecrootxf_to_livecroot_acc_patch   = fill(NaN, np)
        cs.matrix_ctransfer_deadcrootst_to_deadcrootxf_acc_patch = fill(NaN, np)
        cs.matrix_ctransfer_deadcrootxf_to_deadcroot_acc_patch   = fill(NaN, np)
        cs.matrix_ctransfer_grainst_to_grainxf_acc_patch         = fill(NaN, np)
        cs.matrix_ctransfer_grainxf_to_grain_acc_patch           = fill(NaN, np)
        cs.matrix_ctransfer_livestem_to_deadstem_acc_patch        = fill(NaN, np)
        cs.matrix_ctransfer_livecroot_to_deadcroot_acc_patch      = fill(NaN, np)

        # Matrix C turnover accumulation
        cs.matrix_cturnover_leaf_acc_patch             = fill(NaN, np)
        cs.matrix_cturnover_leafst_acc_patch           = fill(NaN, np)
        cs.matrix_cturnover_leafxf_acc_patch           = fill(NaN, np)
        cs.matrix_cturnover_froot_acc_patch            = fill(NaN, np)
        cs.matrix_cturnover_frootst_acc_patch          = fill(NaN, np)
        cs.matrix_cturnover_frootxf_acc_patch          = fill(NaN, np)
        cs.matrix_cturnover_livestem_acc_patch         = fill(NaN, np)
        cs.matrix_cturnover_livestemst_acc_patch       = fill(NaN, np)
        cs.matrix_cturnover_livestemxf_acc_patch       = fill(NaN, np)
        cs.matrix_cturnover_deadstem_acc_patch         = fill(NaN, np)
        cs.matrix_cturnover_deadstemst_acc_patch       = fill(NaN, np)
        cs.matrix_cturnover_deadstemxf_acc_patch       = fill(NaN, np)
        cs.matrix_cturnover_livecroot_acc_patch        = fill(NaN, np)
        cs.matrix_cturnover_livecrootst_acc_patch      = fill(NaN, np)
        cs.matrix_cturnover_livecrootxf_acc_patch      = fill(NaN, np)
        cs.matrix_cturnover_deadcroot_acc_patch        = fill(NaN, np)
        cs.matrix_cturnover_deadcrootst_acc_patch      = fill(NaN, np)
        cs.matrix_cturnover_deadcrootxf_acc_patch      = fill(NaN, np)
        cs.matrix_cturnover_grain_acc_patch            = fill(NaN, np)
        cs.matrix_cturnover_grainst_acc_patch          = fill(NaN, np)
        cs.matrix_cturnover_grainxf_acc_patch          = fill(NaN, np)
    end

    return nothing
end

"""
    cnveg_carbon_state_set_values!(cs, mask_patch, value_patch, mask_col, value_col;
                                   use_matrixcn=false, use_crop=false)

Set carbon state variables for filtered patches and columns.
Ported from `cnveg_carbonstate_type%SetValues`.
"""
function cnveg_carbon_state_set_values!(cs::CNVegCarbonStateData,
                                        mask_patch::BitVector, value_patch::Real,
                                        mask_col::BitVector, value_col::Real;
                                        use_matrixcn::Bool=false, use_crop::Bool=false,
                                        nrepr::Int=NREPR)
    # Patch-level
    for i in eachindex(mask_patch)
        mask_patch[i] || continue
        cs.leafc_patch[i]              = value_patch
        cs.leafc_storage_patch[i]      = value_patch
        cs.leafc_xfer_patch[i]         = value_patch
        cs.leafc_storage_xfer_acc_patch[i] = value_patch
        cs.storage_cdemand_patch[i]    = value_patch
        cs.frootc_patch[i]             = value_patch
        cs.frootc_storage_patch[i]     = value_patch
        cs.frootc_xfer_patch[i]        = value_patch
        cs.livestemc_patch[i]          = value_patch
        cs.livestemc_storage_patch[i]  = value_patch
        cs.livestemc_xfer_patch[i]     = value_patch
        cs.deadstemc_patch[i]          = value_patch
        cs.deadstemc_storage_patch[i]  = value_patch
        cs.deadstemc_xfer_patch[i]     = value_patch
        cs.livecrootc_patch[i]         = value_patch
        cs.livecrootc_storage_patch[i] = value_patch
        cs.livecrootc_xfer_patch[i]    = value_patch
        cs.deadcrootc_patch[i]         = value_patch
        cs.deadcrootc_storage_patch[i] = value_patch
        cs.deadcrootc_xfer_patch[i]    = value_patch

        if use_matrixcn && length(cs.matrix_cap_leafc_patch) > 0
            cs.matrix_cap_leafc_patch[i]              = value_patch
            cs.matrix_cap_leafc_storage_patch[i]      = value_patch
            cs.matrix_cap_leafc_xfer_patch[i]         = value_patch
            cs.matrix_cap_frootc_patch[i]             = value_patch
            cs.matrix_cap_frootc_storage_patch[i]     = value_patch
            cs.matrix_cap_frootc_xfer_patch[i]        = value_patch
            cs.matrix_cap_livestemc_patch[i]          = value_patch
            cs.matrix_cap_livestemc_storage_patch[i]  = value_patch
            cs.matrix_cap_livestemc_xfer_patch[i]     = value_patch
            cs.matrix_cap_deadstemc_patch[i]          = value_patch
            cs.matrix_cap_deadstemc_storage_patch[i]  = value_patch
            cs.matrix_cap_deadstemc_xfer_patch[i]     = value_patch
            cs.matrix_cap_livecrootc_patch[i]         = value_patch
            cs.matrix_cap_livecrootc_storage_patch[i] = value_patch
            cs.matrix_cap_livecrootc_xfer_patch[i]    = value_patch
            cs.matrix_cap_deadcrootc_patch[i]         = value_patch
            cs.matrix_cap_deadcrootc_storage_patch[i] = value_patch
            cs.matrix_cap_deadcrootc_xfer_patch[i]    = value_patch

            cs.leafc0_patch[i]              = value_patch
            cs.leafc0_storage_patch[i]      = value_patch
            cs.leafc0_xfer_patch[i]         = value_patch
            cs.frootc0_patch[i]             = value_patch
            cs.frootc0_storage_patch[i]     = value_patch
            cs.frootc0_xfer_patch[i]        = value_patch
            cs.livestemc0_patch[i]          = value_patch
            cs.livestemc0_storage_patch[i]  = value_patch
            cs.livestemc0_xfer_patch[i]     = value_patch
            cs.deadstemc0_patch[i]          = value_patch
            cs.deadstemc0_storage_patch[i]  = value_patch
            cs.deadstemc0_xfer_patch[i]     = value_patch
            cs.livecrootc0_patch[i]         = value_patch
            cs.livecrootc0_storage_patch[i] = value_patch
            cs.livecrootc0_xfer_patch[i]    = value_patch
            cs.deadcrootc0_patch[i]         = value_patch
            cs.deadcrootc0_storage_patch[i] = value_patch
            cs.deadcrootc0_xfer_patch[i]    = value_patch
            cs.reproc0_patch[i]             = value_patch
            cs.reproc0_storage_patch[i]     = value_patch
            cs.reproc0_xfer_patch[i]        = value_patch

            cs.matrix_calloc_leaf_acc_patch[i]        = value_patch
            cs.matrix_calloc_leafst_acc_patch[i]      = value_patch
            cs.matrix_calloc_froot_acc_patch[i]       = value_patch
            cs.matrix_calloc_frootst_acc_patch[i]     = value_patch
            cs.matrix_calloc_livestem_acc_patch[i]    = value_patch
            cs.matrix_calloc_livestemst_acc_patch[i]  = value_patch
            cs.matrix_calloc_deadstem_acc_patch[i]    = value_patch
            cs.matrix_calloc_deadstemst_acc_patch[i]  = value_patch
            cs.matrix_calloc_livecroot_acc_patch[i]   = value_patch
            cs.matrix_calloc_livecrootst_acc_patch[i] = value_patch
            cs.matrix_calloc_deadcroot_acc_patch[i]   = value_patch
            cs.matrix_calloc_deadcrootst_acc_patch[i] = value_patch

            cs.matrix_ctransfer_leafst_to_leafxf_acc_patch[i]           = value_patch
            cs.matrix_ctransfer_leafxf_to_leaf_acc_patch[i]             = value_patch
            cs.matrix_ctransfer_frootst_to_frootxf_acc_patch[i]         = value_patch
            cs.matrix_ctransfer_frootxf_to_froot_acc_patch[i]           = value_patch
            cs.matrix_ctransfer_livestemst_to_livestemxf_acc_patch[i]   = value_patch
            cs.matrix_ctransfer_livestemxf_to_livestem_acc_patch[i]     = value_patch
            cs.matrix_ctransfer_deadstemst_to_deadstemxf_acc_patch[i]   = value_patch
            cs.matrix_ctransfer_deadstemxf_to_deadstem_acc_patch[i]     = value_patch
            cs.matrix_ctransfer_livecrootst_to_livecrootxf_acc_patch[i] = value_patch
            cs.matrix_ctransfer_livecrootxf_to_livecroot_acc_patch[i]   = value_patch
            cs.matrix_ctransfer_deadcrootst_to_deadcrootxf_acc_patch[i] = value_patch
            cs.matrix_ctransfer_deadcrootxf_to_deadcroot_acc_patch[i]   = value_patch
            cs.matrix_ctransfer_livestem_to_deadstem_acc_patch[i]       = value_patch
            cs.matrix_ctransfer_livecroot_to_deadcroot_acc_patch[i]     = value_patch

            cs.matrix_cturnover_leaf_acc_patch[i]             = value_patch
            cs.matrix_cturnover_leafst_acc_patch[i]           = value_patch
            cs.matrix_cturnover_leafxf_acc_patch[i]           = value_patch
            cs.matrix_cturnover_froot_acc_patch[i]            = value_patch
            cs.matrix_cturnover_frootst_acc_patch[i]          = value_patch
            cs.matrix_cturnover_frootxf_acc_patch[i]          = value_patch
            cs.matrix_cturnover_livestem_acc_patch[i]         = value_patch
            cs.matrix_cturnover_livestemst_acc_patch[i]       = value_patch
            cs.matrix_cturnover_livestemxf_acc_patch[i]       = value_patch
            cs.matrix_cturnover_deadstem_acc_patch[i]         = value_patch
            cs.matrix_cturnover_deadstemst_acc_patch[i]       = value_patch
            cs.matrix_cturnover_deadstemxf_acc_patch[i]       = value_patch
            cs.matrix_cturnover_livecroot_acc_patch[i]        = value_patch
            cs.matrix_cturnover_livecrootst_acc_patch[i]      = value_patch
            cs.matrix_cturnover_livecrootxf_acc_patch[i]      = value_patch
            cs.matrix_cturnover_deadcroot_acc_patch[i]        = value_patch
            cs.matrix_cturnover_deadcrootst_acc_patch[i]      = value_patch
            cs.matrix_cturnover_deadcrootxf_acc_patch[i]      = value_patch
        end

        cs.gresp_storage_patch[i]      = value_patch
        cs.gresp_xfer_patch[i]         = value_patch
        cs.cpool_patch[i]              = value_patch
        cs.xsmrpool_patch[i]           = value_patch
        cs.ctrunc_patch[i]             = value_patch
        cs.dispvegc_patch[i]           = value_patch
        cs.storvegc_patch[i]           = value_patch
        cs.woodc_patch[i]              = value_patch
        cs.totvegc_patch[i]            = value_patch
        cs.totc_patch[i]               = value_patch

        if use_crop
            if use_matrixcn && length(cs.matrix_calloc_grain_acc_patch) > 0
                cs.matrix_calloc_grain_acc_patch[i]                  = value_patch
                cs.matrix_calloc_grainst_acc_patch[i]                = value_patch
                cs.matrix_ctransfer_grainst_to_grainxf_acc_patch[i]  = value_patch
                cs.matrix_ctransfer_grainxf_to_grain_acc_patch[i]    = value_patch
                cs.matrix_cturnover_grain_acc_patch[i]               = value_patch
                cs.matrix_cturnover_grainst_acc_patch[i]             = value_patch
                cs.matrix_cturnover_grainxf_acc_patch[i]             = value_patch
            end
            cs.cropseedc_deficit_patch[i] = value_patch
            cs.xsmrpool_loss_patch[i]     = value_patch
        end
    end

    # Reproductive C (2D) for crop
    if use_crop
        for k in 1:nrepr
            for i in eachindex(mask_patch)
                mask_patch[i] || continue
                cs.reproductivec_patch[i, k]         = value_patch
                cs.reproductivec_storage_patch[i, k]  = value_patch
                cs.reproductivec_xfer_patch[i, k]     = value_patch
            end
        end
        if use_matrixcn && length(cs.matrix_cap_reproc_patch) > 0
            for i in eachindex(mask_col)
                mask_col[i] || continue
                cs.matrix_cap_reproc_patch[i]         = value_patch
                cs.matrix_cap_reproc_storage_patch[i]  = value_patch
                cs.matrix_cap_reproc_xfer_patch[i]     = value_patch
            end
        end
    end

    # Column-level
    for i in eachindex(mask_col)
        mask_col[i] || continue
        cs.rootc_col[i]      = value_col
        cs.leafc_col[i]      = value_col
        cs.deadstemc_col[i]  = value_col
        cs.fuelc_col[i]      = value_col
        cs.fuelc_crop_col[i] = value_col
        cs.totvegc_col[i]    = value_col
        cs.totc_p2c_col[i]   = value_col
    end

    return nothing
end

"""
    cnveg_carbon_state_zero_dwt!(cs, bounds_patch)

Initialize variables needed for dynamic land use.
Ported from `cnveg_carbonstate_type%ZeroDwt`.
"""
function cnveg_carbon_state_zero_dwt!(cs::CNVegCarbonStateData,
                                      bounds_patch::UnitRange{Int})
    for p in bounds_patch
        cs.dispvegc_patch[p] = 0.0
        cs.storvegc_patch[p] = 0.0
        cs.totc_patch[p]     = 0.0
    end
    return nothing
end

# --- GPU kernelization (one thread per patch; every reduction into a LOCAL scalar,
#     each output written ONCE — the KA-CPU-under-check-bounds `arr[p] += …`-in-loop
#     miscompile is avoided). Touched fields grouped into one @adapt_structure bundle
#     (one Metal buffer arg); the KA CPU backend runs it byte-identical on Arrays.
struct _CVCSumView{V,M}
    leafc::V; frootc::V; livestemc::V; deadstemc::V; livecrootc::V; deadcrootc::V
    cpool::V; leafc_storage::V; frootc_storage::V; livestemc_storage::V
    deadstemc_storage::V; livecrootc_storage::V; deadcrootc_storage::V
    leafc_xfer::V; frootc_xfer::V; livestemc_xfer::V; deadstemc_xfer::V
    livecrootc_xfer::V; deadcrootc_xfer::V; gresp_storage::V; gresp_xfer::V
    xsmrpool::V; ctrunc::V; cropseedc_deficit::V; xsmrpool_loss::V
    reproductivec_storage::M; reproductivec_xfer::M; reproductivec::M
    dispvegc::V; storvegc::V; totvegc::V; totc::V; woodc::V
end
Adapt.@adapt_structure _CVCSumView

@kernel function _cvc_summary_kernel!(@Const(mask), b, @Const(patch_itype),
        lo::Int, hi::Int, do_crop::Bool, use_crop::Bool, npcropmin::Int, nrepr::Int)
    p = @index(Global)
    @inbounds if lo <= p <= hi && mask[p]
        # displayed vegetation carbon, excluding storage and cpool (DISPVEGC)
        dispvegc = b.leafc[p] + b.frootc[p] + b.livestemc[p] +
                   b.deadstemc[p] + b.livecrootc[p] + b.deadcrootc[p]

        # stored vegetation carbon, excluding cpool (STORVEGC)
        storvegc = b.cpool[p] + b.leafc_storage[p] + b.frootc_storage[p] +
                   b.livestemc_storage[p] + b.deadstemc_storage[p] +
                   b.livecrootc_storage[p] + b.deadcrootc_storage[p] +
                   b.leafc_xfer[p] + b.frootc_xfer[p] + b.livestemc_xfer[p] +
                   b.deadstemc_xfer[p] + b.livecrootc_xfer[p] + b.deadcrootc_xfer[p] +
                   b.gresp_storage[p] + b.gresp_xfer[p]

        # Accumulate directly into the storvegc/dispvegc LOCAL scalars (safe — the KA
        # trap is only for array-element `+=` in a loop) so the term order is exactly
        # the original per-k sequence → byte-identical.
        if do_crop && patch_itype[p] >= npcropmin
            for k in 1:nrepr
                storvegc += b.reproductivec_storage[p, k] + b.reproductivec_xfer[p, k]
                dispvegc += b.reproductivec[p, k]
            end
        end

        b.dispvegc[p] = dispvegc
        b.storvegc[p] = storvegc

        totvegc = dispvegc + storvegc              # TOTVEGC (excludes cpool)
        b.totvegc[p] = totvegc

        totc = totvegc + b.xsmrpool[p] + b.ctrunc[p]
        if use_crop
            totc += b.cropseedc_deficit[p] + b.xsmrpool_loss[p]
        end
        b.totc[p] = totc

        b.woodc[p] = b.deadstemc[p] + b.livestemc[p] + b.deadcrootc[p] + b.livecrootc[p]
    end
end

"""
    cnveg_carbon_state_summary!(cs, mask_patch, bounds_patch;
                                 use_crop=false, patch_itype=nothing,
                                 npcropmin=0, nrepr=NREPR)

Perform patch-level carbon summary calculations.
Ported from `cnveg_carbonstate_type%Summary_carbonstate`.

Note: column-level p2c averaging is stubbed out pending subgridAveMod port.
"""
function cnveg_carbon_state_summary!(cs::CNVegCarbonStateData,
                                     mask_patch::AbstractVector{Bool},
                                     bounds_patch::UnitRange{Int};
                                     use_crop::Bool=false,
                                     patch_itype::Union{AbstractVector{<:Integer},Nothing}=nothing,
                                     npcropmin::Int=0,
                                     nrepr::Int=NREPR)
    isempty(bounds_patch) && return nothing
    b = _CVCSumView(cs.leafc_patch, cs.frootc_patch, cs.livestemc_patch,
        cs.deadstemc_patch, cs.livecrootc_patch, cs.deadcrootc_patch, cs.cpool_patch,
        cs.leafc_storage_patch, cs.frootc_storage_patch, cs.livestemc_storage_patch,
        cs.deadstemc_storage_patch, cs.livecrootc_storage_patch, cs.deadcrootc_storage_patch,
        cs.leafc_xfer_patch, cs.frootc_xfer_patch, cs.livestemc_xfer_patch,
        cs.deadstemc_xfer_patch, cs.livecrootc_xfer_patch, cs.deadcrootc_xfer_patch,
        cs.gresp_storage_patch, cs.gresp_xfer_patch, cs.xsmrpool_patch, cs.ctrunc_patch,
        cs.cropseedc_deficit_patch, cs.xsmrpool_loss_patch,
        cs.reproductivec_storage_patch, cs.reproductivec_xfer_patch, cs.reproductivec_patch,
        cs.dispvegc_patch, cs.storvegc_patch, cs.totvegc_patch, cs.totc_patch, cs.woodc_patch)
    # Crop-reproductive branch only active with a patch_itype; move it to the state
    # backend (preserving Int) when active, else pass an empty device Int vector.
    T = eltype(cs.dispvegc_patch)
    do_crop = use_crop && patch_itype !== nothing
    pit = do_crop ? _to_backend_like(cs.dispvegc_patch, T, patch_itype) :
                    _to_backend_like(cs.dispvegc_patch, T, Int[])
    _launch!(_cvc_summary_kernel!, mask_patch, b, pit,
        first(bounds_patch), last(bounds_patch), do_crop, use_crop, npcropmin, nrepr)

    # Column-level p2c averaging is stubbed pending subgridAveMod port
    return nothing
end

"""
    cnveg_carbon_state_init_cold!(cs, bounds_patch; ratio=1.0,
                                   use_matrixcn=false, use_crop=false)

Initialize cold-start conditions for CN vegetation carbon state variables.
Simplified version that initializes all patches as soil/crop type.

Ported from `cnveg_carbonstate_type%InitCold`.
"""
function cnveg_carbon_state_init_cold!(cs::CNVegCarbonStateData,
                                       bounds_patch::UnitRange{Int};
                                       ratio::Real=1.0,
                                       use_matrixcn::Bool=false,
                                       use_crop::Bool=false)
    for p in bounds_patch
        cs.leafcmax_patch[p] = 0.0

        # Default: deciduous non-crop (leafc=0, storage=initial_vegC)
        cs.leafc_patch[p]         = 0.0
        cs.leafc_storage_patch[p] = INITIAL_VEGC * ratio
        cs.frootc_patch[p]        = 0.0
        cs.frootc_storage_patch[p] = INITIAL_VEGC * ratio

        cs.leafc_xfer_patch[p]                = 0.0
        cs.leafc_storage_xfer_acc_patch[p]    = 0.0
        cs.storage_cdemand_patch[p]           = 0.0
        cs.frootc_xfer_patch[p]               = 0.0

        cs.livestemc_patch[p]                 = 0.0
        cs.livestemc_storage_patch[p]         = 0.0
        cs.livestemc_xfer_patch[p]            = 0.0

        cs.deadstemc_patch[p]                 = 0.0
        cs.deadstemc_storage_patch[p]         = 0.0
        cs.deadstemc_xfer_patch[p]            = 0.0

        cs.livecrootc_patch[p]                = 0.0
        cs.livecrootc_storage_patch[p]        = 0.0
        cs.livecrootc_xfer_patch[p]           = 0.0

        cs.deadcrootc_patch[p]                = 0.0
        cs.deadcrootc_storage_patch[p]        = 0.0
        cs.deadcrootc_xfer_patch[p]           = 0.0

        cs.gresp_storage_patch[p]             = 0.0
        cs.gresp_xfer_patch[p]                = 0.0
        cs.cpool_patch[p]                     = 0.0
        cs.xsmrpool_patch[p]                  = 0.0
        cs.ctrunc_patch[p]                    = 0.0
        cs.dispvegc_patch[p]                  = 0.0
        cs.storvegc_patch[p]                  = 0.0
        cs.woodc_patch[p]                     = 0.0
        cs.totc_patch[p]                      = 0.0

        if use_matrixcn && length(cs.matrix_cap_leafc_patch) > 0
            cs.matrix_cap_leafc_patch[p]          = 0.0
            cs.matrix_cap_leafc_storage_patch[p]  = INITIAL_VEGC * ratio
            cs.matrix_cap_leafc_xfer_patch[p]     = 0.0
            cs.matrix_cap_frootc_patch[p]         = 0.0
            cs.matrix_cap_frootc_storage_patch[p] = INITIAL_VEGC * ratio
            cs.matrix_cap_frootc_xfer_patch[p]    = 0.0

            cs.matrix_cap_livestemc_patch[p]          = 0.0
            cs.matrix_cap_livestemc_storage_patch[p]  = 0.0
            cs.matrix_cap_livestemc_xfer_patch[p]     = 0.0
            cs.matrix_cap_deadstemc_patch[p]          = 0.0
            cs.matrix_cap_deadstemc_storage_patch[p]  = 0.0
            cs.matrix_cap_deadstemc_xfer_patch[p]     = 0.0
            cs.matrix_cap_livecrootc_patch[p]         = 0.0
            cs.matrix_cap_livecrootc_storage_patch[p] = 0.0
            cs.matrix_cap_livecrootc_xfer_patch[p]    = 0.0
            cs.matrix_cap_deadcrootc_patch[p]         = 0.0
            cs.matrix_cap_deadcrootc_storage_patch[p] = 0.0
            cs.matrix_cap_deadcrootc_xfer_patch[p]    = 0.0

            # Initial pool sizes for SASU
            cs.leafc0_patch[p]              = 1.0e-30
            cs.leafc0_storage_patch[p]      = 1.0e-30
            cs.leafc0_xfer_patch[p]         = 1.0e-30
            cs.frootc0_patch[p]             = 1.0e-30
            cs.frootc0_storage_patch[p]     = 1.0e-30
            cs.frootc0_xfer_patch[p]        = 1.0e-30
            cs.livestemc0_patch[p]          = 1.0e-30
            cs.livestemc0_storage_patch[p]  = 1.0e-30
            cs.livestemc0_xfer_patch[p]     = 1.0e-30
            cs.deadstemc0_patch[p]          = 1.0e-30
            cs.deadstemc0_storage_patch[p]  = 1.0e-30
            cs.deadstemc0_xfer_patch[p]     = 1.0e-30
            cs.livecrootc0_patch[p]         = 1.0e-30
            cs.livecrootc0_storage_patch[p] = 1.0e-30
            cs.livecrootc0_xfer_patch[p]    = 1.0e-30
            cs.deadcrootc0_patch[p]         = 1.0e-30
            cs.deadcrootc0_storage_patch[p] = 1.0e-30
            cs.deadcrootc0_xfer_patch[p]    = 1.0e-30
            cs.reproc0_patch[p]             = 1.0e-30
            cs.reproc0_storage_patch[p]     = 1.0e-30
            cs.reproc0_xfer_patch[p]        = 1.0e-30

            # SASU save fields
            cs.leafc_SASUsave_patch[p]              = 0.0
            cs.leafc_storage_SASUsave_patch[p]      = 0.0
            cs.leafc_xfer_SASUsave_patch[p]         = 0.0
            cs.frootc_SASUsave_patch[p]             = 0.0
            cs.frootc_storage_SASUsave_patch[p]     = 0.0
            cs.frootc_xfer_SASUsave_patch[p]        = 0.0
            cs.livestemc_SASUsave_patch[p]          = 0.0
            cs.livestemc_storage_SASUsave_patch[p]  = 0.0
            cs.livestemc_xfer_SASUsave_patch[p]     = 0.0
            cs.deadstemc_SASUsave_patch[p]          = 0.0
            cs.deadstemc_storage_SASUsave_patch[p]  = 0.0
            cs.deadstemc_xfer_SASUsave_patch[p]     = 0.0
            cs.livecrootc_SASUsave_patch[p]         = 0.0
            cs.livecrootc_storage_SASUsave_patch[p] = 0.0
            cs.livecrootc_xfer_SASUsave_patch[p]    = 0.0
            cs.deadcrootc_SASUsave_patch[p]         = 0.0
            cs.deadcrootc_storage_SASUsave_patch[p] = 0.0
            cs.deadcrootc_xfer_SASUsave_patch[p]    = 0.0
            cs.grainc_SASUsave_patch[p]             = 0.0
            cs.grainc_storage_SASUsave_patch[p]     = 0.0

            # Matrix accumulation fields
            cs.matrix_calloc_leaf_acc_patch[p]        = 0.0
            cs.matrix_calloc_leafst_acc_patch[p]      = 0.0
            cs.matrix_calloc_froot_acc_patch[p]       = 0.0
            cs.matrix_calloc_frootst_acc_patch[p]     = 0.0
            cs.matrix_calloc_livestem_acc_patch[p]    = 0.0
            cs.matrix_calloc_livestemst_acc_patch[p]  = 0.0
            cs.matrix_calloc_deadstem_acc_patch[p]    = 0.0
            cs.matrix_calloc_deadstemst_acc_patch[p]  = 0.0
            cs.matrix_calloc_livecroot_acc_patch[p]   = 0.0
            cs.matrix_calloc_livecrootst_acc_patch[p] = 0.0
            cs.matrix_calloc_deadcroot_acc_patch[p]   = 0.0
            cs.matrix_calloc_deadcrootst_acc_patch[p] = 0.0

            cs.matrix_ctransfer_leafst_to_leafxf_acc_patch[p]           = 0.0
            cs.matrix_ctransfer_leafxf_to_leaf_acc_patch[p]             = 0.0
            cs.matrix_ctransfer_frootst_to_frootxf_acc_patch[p]         = 0.0
            cs.matrix_ctransfer_frootxf_to_froot_acc_patch[p]           = 0.0
            cs.matrix_ctransfer_livestemst_to_livestemxf_acc_patch[p]   = 0.0
            cs.matrix_ctransfer_livestemxf_to_livestem_acc_patch[p]     = 0.0
            cs.matrix_ctransfer_deadstemst_to_deadstemxf_acc_patch[p]   = 0.0
            cs.matrix_ctransfer_deadstemxf_to_deadstem_acc_patch[p]     = 0.0
            cs.matrix_ctransfer_livecrootst_to_livecrootxf_acc_patch[p] = 0.0
            cs.matrix_ctransfer_livecrootxf_to_livecroot_acc_patch[p]   = 0.0
            cs.matrix_ctransfer_deadcrootst_to_deadcrootxf_acc_patch[p] = 0.0
            cs.matrix_ctransfer_deadcrootxf_to_deadcroot_acc_patch[p]   = 0.0
            cs.matrix_ctransfer_livestem_to_deadstem_acc_patch[p]       = 0.0
            cs.matrix_ctransfer_livecroot_to_deadcroot_acc_patch[p]     = 0.0

            cs.matrix_cturnover_leaf_acc_patch[p]             = 0.0
            cs.matrix_cturnover_leafst_acc_patch[p]           = 0.0
            cs.matrix_cturnover_leafxf_acc_patch[p]           = 0.0
            cs.matrix_cturnover_froot_acc_patch[p]            = 0.0
            cs.matrix_cturnover_frootst_acc_patch[p]          = 0.0
            cs.matrix_cturnover_frootxf_acc_patch[p]          = 0.0
            cs.matrix_cturnover_livestem_acc_patch[p]         = 0.0
            cs.matrix_cturnover_livestemst_acc_patch[p]       = 0.0
            cs.matrix_cturnover_livestemxf_acc_patch[p]       = 0.0
            cs.matrix_cturnover_deadstem_acc_patch[p]         = 0.0
            cs.matrix_cturnover_deadstemst_acc_patch[p]       = 0.0
            cs.matrix_cturnover_deadstemxf_acc_patch[p]       = 0.0
            cs.matrix_cturnover_livecroot_acc_patch[p]        = 0.0
            cs.matrix_cturnover_livecrootst_acc_patch[p]      = 0.0
            cs.matrix_cturnover_livecrootxf_acc_patch[p]      = 0.0
            cs.matrix_cturnover_deadcroot_acc_patch[p]        = 0.0
            cs.matrix_cturnover_deadcrootst_acc_patch[p]      = 0.0
            cs.matrix_cturnover_deadcrootxf_acc_patch[p]      = 0.0
        end

        if use_crop
            for k in 1:size(cs.reproductivec_patch, 2)
                cs.reproductivec_patch[p, k]         = 0.0
                cs.reproductivec_storage_patch[p, k] = 0.0
                cs.reproductivec_xfer_patch[p, k]    = 0.0
            end
            cs.cropseedc_deficit_patch[p] = 0.0
            cs.xsmrpool_loss_patch[p]     = 0.0
            if use_matrixcn && length(cs.matrix_cap_reproc_patch) > 0
                cs.matrix_cap_reproc_patch[p]         = 0.0
                cs.matrix_cap_reproc_storage_patch[p] = 0.0
                cs.matrix_cap_reproc_xfer_patch[p]    = 0.0
                cs.matrix_calloc_grain_acc_patch[p]                  = 0.0
                cs.matrix_calloc_grainst_acc_patch[p]                = 0.0
                cs.matrix_ctransfer_grainst_to_grainxf_acc_patch[p]  = 0.0
                cs.matrix_ctransfer_grainxf_to_grain_acc_patch[p]    = 0.0
                cs.matrix_cturnover_grain_acc_patch[p]               = 0.0
                cs.matrix_cturnover_grainst_acc_patch[p]             = 0.0
                cs.matrix_cturnover_grainxf_acc_patch[p]             = 0.0
            end
        end
    end

    return nothing
end

# ==========================================================================
# Stubs for infrastructure-dependent subroutines
# ==========================================================================

"""
    cnveg_carbon_state_init_history!(cs, bounds_patch, bounds_col)

Register CN vegetation carbon state fields for history file output.
Stub: requires histFileMod infrastructure.
"""
function cnveg_carbon_state_init_history!(cs::CNVegCarbonStateData,
                                          bounds_patch::UnitRange{Int},
                                          bounds_col::UnitRange{Int})
    return nothing
end

"""
    cnveg_carbon_state_restart!(cs, bounds_patch, bounds_col; flag="read")

Read/write CN vegetation carbon state from/to restart file.
Stub: requires NetCDF/restart infrastructure.
"""
function cnveg_carbon_state_restart!(cs::CNVegCarbonStateData,
                                     bounds_patch::UnitRange{Int},
                                     bounds_col::UnitRange{Int};
                                     flag::String="read")
    return nothing
end

"""
    cnveg_carbon_state_dynamic_patch_adjustments!(...)

Adjust state variables when patch areas change due to dynamic landuse.
Stub: requires dynPatchStateUpdaterMod and CNVegComputeSeedMod infrastructure.
"""
function cnveg_carbon_state_dynamic_patch_adjustments!(cs::CNVegCarbonStateData,
                                                       bounds_patch::UnitRange{Int})
    return nothing
end
