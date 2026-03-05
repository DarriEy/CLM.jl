# ==========================================================================
# Ported from: src/biogeochem/CNVegNitrogenStateType.F90
# Vegetation nitrogen state data type allocation and initialization
# ==========================================================================

"""
    CNVegNitrogenStateData

Vegetation nitrogen state data structure. Holds nitrogen pools at patch, column,
and gridcell levels.

Ported from `cnveg_nitrogenstate_type` in `CNVegNitrogenStateType.F90`.
"""
Base.@kwdef mutable struct CNVegNitrogenStateData{FT<:AbstractFloat}
    # --- Patch-level 2D (patch x nrepr) ---
    reproductiven_patch               ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # (gN/m2) reproductive (grain) N
    reproductiven_storage_patch       ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # (gN/m2) reproductive N storage
    reproductiven_xfer_patch          ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # (gN/m2) reproductive N transfer

    # --- Patch-level 1D core pools ---
    leafn_patch                       ::Vector{FT} = Float64[]  # (gN/m2) leaf N
    leafn_storage_patch               ::Vector{FT} = Float64[]  # (gN/m2) leaf N storage
    leafn_xfer_patch                  ::Vector{FT} = Float64[]  # (gN/m2) leaf N transfer
    leafn_storage_xfer_acc_patch      ::Vector{FT} = Float64[]  # (gN/m2) accumulated leaf N transfer
    storage_ndemand_patch             ::Vector{FT} = Float64[]  # (gN/m2) N demand during offset period
    frootn_patch                      ::Vector{FT} = Float64[]  # (gN/m2) fine root N
    frootn_storage_patch              ::Vector{FT} = Float64[]  # (gN/m2) fine root N storage
    frootn_xfer_patch                 ::Vector{FT} = Float64[]  # (gN/m2) fine root N transfer
    livestemn_patch                   ::Vector{FT} = Float64[]  # (gN/m2) live stem N
    livestemn_storage_patch           ::Vector{FT} = Float64[]  # (gN/m2) live stem N storage
    livestemn_xfer_patch              ::Vector{FT} = Float64[]  # (gN/m2) live stem N transfer
    deadstemn_patch                   ::Vector{FT} = Float64[]  # (gN/m2) dead stem N
    deadstemn_storage_patch           ::Vector{FT} = Float64[]  # (gN/m2) dead stem N storage
    deadstemn_xfer_patch              ::Vector{FT} = Float64[]  # (gN/m2) dead stem N transfer
    livecrootn_patch                  ::Vector{FT} = Float64[]  # (gN/m2) live coarse root N
    livecrootn_storage_patch          ::Vector{FT} = Float64[]  # (gN/m2) live coarse root N storage
    livecrootn_xfer_patch             ::Vector{FT} = Float64[]  # (gN/m2) live coarse root N transfer
    deadcrootn_patch                  ::Vector{FT} = Float64[]  # (gN/m2) dead coarse root N
    deadcrootn_storage_patch          ::Vector{FT} = Float64[]  # (gN/m2) dead coarse root N storage
    deadcrootn_xfer_patch             ::Vector{FT} = Float64[]  # (gN/m2) dead coarse root N transfer
    retransn_patch                    ::Vector{FT} = Float64[]  # (gN/m2) plant pool of retranslocated N
    npool_patch                       ::Vector{FT} = Float64[]  # (gN/m2) temporary plant N pool
    ntrunc_patch                      ::Vector{FT} = Float64[]  # (gN/m2) patch-level sink for N truncation
    cropseedn_deficit_patch           ::Vector{FT} = Float64[]  # (gN/m2) pool for seeding new crop growth (negative)

    # --- Matrix CN capacity fields (patch-level) ---
    matrix_cap_repron_patch             ::Vector{FT} = Float64[]  # (gN/m2) capacity of grain N
    matrix_cap_repron_storage_patch     ::Vector{FT} = Float64[]  # (gN/m2) capacity of grain N storage
    matrix_cap_repron_xfer_patch        ::Vector{FT} = Float64[]  # (gN/m2) capacity of grain N transfer
    matrix_cap_leafn_patch              ::Vector{FT} = Float64[]  # (gN/m2) capacity of leaf N
    matrix_cap_leafn_storage_patch      ::Vector{FT} = Float64[]  # (gN/m2) capacity of leaf N storage
    matrix_cap_leafn_xfer_patch         ::Vector{FT} = Float64[]  # (gN/m2) capacity of leaf N transfer
    matrix_cap_frootn_patch             ::Vector{FT} = Float64[]  # (gN/m2) capacity of fine root N
    matrix_cap_frootn_storage_patch     ::Vector{FT} = Float64[]  # (gN/m2) capacity of fine root N storage
    matrix_cap_frootn_xfer_patch        ::Vector{FT} = Float64[]  # (gN/m2) capacity of fine root N transfer
    matrix_cap_livestemn_patch          ::Vector{FT} = Float64[]  # (gN/m2) capacity of live stem N
    matrix_cap_livestemn_storage_patch  ::Vector{FT} = Float64[]  # (gN/m2) capacity of live stem N storage
    matrix_cap_livestemn_xfer_patch     ::Vector{FT} = Float64[]  # (gN/m2) capacity of live stem N transfer
    matrix_cap_deadstemn_patch          ::Vector{FT} = Float64[]  # (gN/m2) capacity of dead stem N
    matrix_cap_deadstemn_storage_patch  ::Vector{FT} = Float64[]  # (gN/m2) capacity of dead stem N storage
    matrix_cap_deadstemn_xfer_patch     ::Vector{FT} = Float64[]  # (gN/m2) capacity of dead stem N transfer
    matrix_cap_livecrootn_patch         ::Vector{FT} = Float64[]  # (gN/m2) capacity of live coarse root N
    matrix_cap_livecrootn_storage_patch ::Vector{FT} = Float64[]  # (gN/m2) capacity of live coarse root N storage
    matrix_cap_livecrootn_xfer_patch    ::Vector{FT} = Float64[]  # (gN/m2) capacity of live coarse root N transfer
    matrix_cap_deadcrootn_patch         ::Vector{FT} = Float64[]  # (gN/m2) capacity of dead coarse root N
    matrix_cap_deadcrootn_storage_patch ::Vector{FT} = Float64[]  # (gN/m2) capacity of dead coarse root N storage
    matrix_cap_deadcrootn_xfer_patch    ::Vector{FT} = Float64[]  # (gN/m2) capacity of dead coarse root N transfer

    # --- Initial pool size of year for matrix (SASU) ---
    leafn0_patch                      ::Vector{FT} = Float64[]  # (gN/m2) initial leaf N for SASU
    leafn0_storage_patch              ::Vector{FT} = Float64[]  # (gN/m2) initial leaf N storage for SASU
    leafn0_xfer_patch                 ::Vector{FT} = Float64[]  # (gN/m2) initial leaf N transfer for SASU
    frootn0_patch                     ::Vector{FT} = Float64[]  # (gN/m2) initial fine root N for SASU
    frootn0_storage_patch             ::Vector{FT} = Float64[]  # (gN/m2) initial fine root N storage for SASU
    frootn0_xfer_patch                ::Vector{FT} = Float64[]  # (gN/m2) initial fine root N transfer for SASU
    livestemn0_patch                  ::Vector{FT} = Float64[]  # (gN/m2) initial live stem N for SASU
    livestemn0_storage_patch          ::Vector{FT} = Float64[]  # (gN/m2) initial live stem N storage for SASU
    livestemn0_xfer_patch             ::Vector{FT} = Float64[]  # (gN/m2) initial live stem N transfer for SASU
    deadstemn0_patch                  ::Vector{FT} = Float64[]  # (gN/m2) initial dead stem N for SASU
    deadstemn0_storage_patch          ::Vector{FT} = Float64[]  # (gN/m2) initial dead stem N storage for SASU
    deadstemn0_xfer_patch             ::Vector{FT} = Float64[]  # (gN/m2) initial dead stem N transfer for SASU
    livecrootn0_patch                 ::Vector{FT} = Float64[]  # (gN/m2) initial live coarse root N for SASU
    livecrootn0_storage_patch         ::Vector{FT} = Float64[]  # (gN/m2) initial live coarse root N storage for SASU
    livecrootn0_xfer_patch            ::Vector{FT} = Float64[]  # (gN/m2) initial live coarse root N transfer for SASU
    deadcrootn0_patch                 ::Vector{FT} = Float64[]  # (gN/m2) initial dead coarse root N for SASU
    deadcrootn0_storage_patch         ::Vector{FT} = Float64[]  # (gN/m2) initial dead coarse root N storage for SASU
    deadcrootn0_xfer_patch            ::Vector{FT} = Float64[]  # (gN/m2) initial dead coarse root N transfer for SASU
    retransn0_patch                   ::Vector{FT} = Float64[]  # (gN/m2) initial retranslocated N for SASU
    repron0_patch                     ::Vector{FT} = Float64[]  # (gN/m2) initial grain N for SASU
    repron0_storage_patch             ::Vector{FT} = Float64[]  # (gN/m2) initial grain N storage for SASU
    repron0_xfer_patch                ::Vector{FT} = Float64[]  # (gN/m2) initial grain N transfer for SASU

    # --- SASU save fields ---
    leafn_SASUsave_patch              ::Vector{FT} = Float64[]
    leafn_storage_SASUsave_patch      ::Vector{FT} = Float64[]
    leafn_xfer_SASUsave_patch         ::Vector{FT} = Float64[]
    frootn_SASUsave_patch             ::Vector{FT} = Float64[]
    frootn_storage_SASUsave_patch     ::Vector{FT} = Float64[]
    frootn_xfer_SASUsave_patch        ::Vector{FT} = Float64[]
    livestemn_SASUsave_patch          ::Vector{FT} = Float64[]
    livestemn_storage_SASUsave_patch  ::Vector{FT} = Float64[]
    livestemn_xfer_SASUsave_patch     ::Vector{FT} = Float64[]
    deadstemn_SASUsave_patch          ::Vector{FT} = Float64[]
    deadstemn_storage_SASUsave_patch  ::Vector{FT} = Float64[]
    deadstemn_xfer_SASUsave_patch     ::Vector{FT} = Float64[]
    livecrootn_SASUsave_patch         ::Vector{FT} = Float64[]
    livecrootn_storage_SASUsave_patch ::Vector{FT} = Float64[]
    livecrootn_xfer_SASUsave_patch    ::Vector{FT} = Float64[]
    deadcrootn_SASUsave_patch         ::Vector{FT} = Float64[]
    deadcrootn_storage_SASUsave_patch ::Vector{FT} = Float64[]
    deadcrootn_xfer_SASUsave_patch    ::Vector{FT} = Float64[]
    grainn_SASUsave_patch             ::Vector{FT} = Float64[]
    grainn_storage_SASUsave_patch     ::Vector{FT} = Float64[]

    # --- Matrix N allocation accumulation variables (annual, patch-level) ---
    matrix_nalloc_leaf_acc_patch        ::Vector{FT} = Float64[]
    matrix_nalloc_leafst_acc_patch      ::Vector{FT} = Float64[]
    matrix_nalloc_froot_acc_patch       ::Vector{FT} = Float64[]
    matrix_nalloc_frootst_acc_patch     ::Vector{FT} = Float64[]
    matrix_nalloc_livestem_acc_patch    ::Vector{FT} = Float64[]
    matrix_nalloc_livestemst_acc_patch  ::Vector{FT} = Float64[]
    matrix_nalloc_deadstem_acc_patch    ::Vector{FT} = Float64[]
    matrix_nalloc_deadstemst_acc_patch  ::Vector{FT} = Float64[]
    matrix_nalloc_livecroot_acc_patch   ::Vector{FT} = Float64[]
    matrix_nalloc_livecrootst_acc_patch ::Vector{FT} = Float64[]
    matrix_nalloc_deadcroot_acc_patch   ::Vector{FT} = Float64[]
    matrix_nalloc_deadcrootst_acc_patch ::Vector{FT} = Float64[]
    matrix_nalloc_grain_acc_patch       ::Vector{FT} = Float64[]
    matrix_nalloc_grainst_acc_patch     ::Vector{FT} = Float64[]

    # --- Matrix N transfer accumulation (annual, patch-level) ---
    matrix_ntransfer_leafst_to_leafxf_acc_patch           ::Vector{FT} = Float64[]
    matrix_ntransfer_leafxf_to_leaf_acc_patch             ::Vector{FT} = Float64[]
    matrix_ntransfer_frootst_to_frootxf_acc_patch         ::Vector{FT} = Float64[]
    matrix_ntransfer_frootxf_to_froot_acc_patch           ::Vector{FT} = Float64[]
    matrix_ntransfer_livestemst_to_livestemxf_acc_patch    ::Vector{FT} = Float64[]
    matrix_ntransfer_livestemxf_to_livestem_acc_patch      ::Vector{FT} = Float64[]
    matrix_ntransfer_deadstemst_to_deadstemxf_acc_patch    ::Vector{FT} = Float64[]
    matrix_ntransfer_deadstemxf_to_deadstem_acc_patch      ::Vector{FT} = Float64[]
    matrix_ntransfer_livecrootst_to_livecrootxf_acc_patch ::Vector{FT} = Float64[]
    matrix_ntransfer_livecrootxf_to_livecroot_acc_patch   ::Vector{FT} = Float64[]
    matrix_ntransfer_deadcrootst_to_deadcrootxf_acc_patch ::Vector{FT} = Float64[]
    matrix_ntransfer_deadcrootxf_to_deadcroot_acc_patch   ::Vector{FT} = Float64[]
    matrix_ntransfer_grainst_to_grainxf_acc_patch         ::Vector{FT} = Float64[]
    matrix_ntransfer_grainxf_to_grain_acc_patch           ::Vector{FT} = Float64[]
    matrix_ntransfer_livestem_to_deadstem_acc_patch        ::Vector{FT} = Float64[]
    matrix_ntransfer_livecroot_to_deadcroot_acc_patch      ::Vector{FT} = Float64[]

    # --- Matrix N retranslocation transfer accumulation (annual, patch-level) ---
    matrix_ntransfer_retransn_to_leaf_acc_patch           ::Vector{FT} = Float64[]
    matrix_ntransfer_retransn_to_leafst_acc_patch         ::Vector{FT} = Float64[]
    matrix_ntransfer_retransn_to_froot_acc_patch          ::Vector{FT} = Float64[]
    matrix_ntransfer_retransn_to_frootst_acc_patch        ::Vector{FT} = Float64[]
    matrix_ntransfer_retransn_to_livestem_acc_patch       ::Vector{FT} = Float64[]
    matrix_ntransfer_retransn_to_livestemst_acc_patch     ::Vector{FT} = Float64[]
    matrix_ntransfer_retransn_to_deadstem_acc_patch       ::Vector{FT} = Float64[]
    matrix_ntransfer_retransn_to_deadstemst_acc_patch     ::Vector{FT} = Float64[]
    matrix_ntransfer_retransn_to_livecroot_acc_patch      ::Vector{FT} = Float64[]
    matrix_ntransfer_retransn_to_livecrootst_acc_patch    ::Vector{FT} = Float64[]
    matrix_ntransfer_retransn_to_deadcroot_acc_patch      ::Vector{FT} = Float64[]
    matrix_ntransfer_retransn_to_deadcrootst_acc_patch    ::Vector{FT} = Float64[]
    matrix_ntransfer_retransn_to_grain_acc_patch          ::Vector{FT} = Float64[]
    matrix_ntransfer_retransn_to_grainst_acc_patch        ::Vector{FT} = Float64[]

    # --- Matrix N from-pool to retranslocation transfer accumulation ---
    matrix_ntransfer_leaf_to_retransn_acc_patch           ::Vector{FT} = Float64[]
    matrix_ntransfer_froot_to_retransn_acc_patch          ::Vector{FT} = Float64[]
    matrix_ntransfer_livestem_to_retransn_acc_patch       ::Vector{FT} = Float64[]
    matrix_ntransfer_livecroot_to_retransn_acc_patch      ::Vector{FT} = Float64[]

    # --- Matrix N turnover accumulation (annual, patch-level) ---
    matrix_nturnover_leaf_acc_patch             ::Vector{FT} = Float64[]
    matrix_nturnover_leafst_acc_patch           ::Vector{FT} = Float64[]
    matrix_nturnover_leafxf_acc_patch           ::Vector{FT} = Float64[]
    matrix_nturnover_froot_acc_patch            ::Vector{FT} = Float64[]
    matrix_nturnover_frootst_acc_patch          ::Vector{FT} = Float64[]
    matrix_nturnover_frootxf_acc_patch          ::Vector{FT} = Float64[]
    matrix_nturnover_livestem_acc_patch         ::Vector{FT} = Float64[]
    matrix_nturnover_livestemst_acc_patch       ::Vector{FT} = Float64[]
    matrix_nturnover_livestemxf_acc_patch       ::Vector{FT} = Float64[]
    matrix_nturnover_deadstem_acc_patch         ::Vector{FT} = Float64[]
    matrix_nturnover_deadstemst_acc_patch       ::Vector{FT} = Float64[]
    matrix_nturnover_deadstemxf_acc_patch       ::Vector{FT} = Float64[]
    matrix_nturnover_livecroot_acc_patch        ::Vector{FT} = Float64[]
    matrix_nturnover_livecrootst_acc_patch      ::Vector{FT} = Float64[]
    matrix_nturnover_livecrootxf_acc_patch      ::Vector{FT} = Float64[]
    matrix_nturnover_deadcroot_acc_patch        ::Vector{FT} = Float64[]
    matrix_nturnover_deadcrootst_acc_patch      ::Vector{FT} = Float64[]
    matrix_nturnover_deadcrootxf_acc_patch      ::Vector{FT} = Float64[]
    matrix_nturnover_grain_acc_patch            ::Vector{FT} = Float64[]
    matrix_nturnover_grainst_acc_patch          ::Vector{FT} = Float64[]
    matrix_nturnover_grainxf_acc_patch          ::Vector{FT} = Float64[]
    matrix_nturnover_retransn_acc_patch         ::Vector{FT} = Float64[]

    # --- Summary (diagnostic) state variables ---
    dispvegn_patch                    ::Vector{FT} = Float64[]  # (gN/m2) displayed veg nitrogen, excluding storage
    storvegn_patch                    ::Vector{FT} = Float64[]  # (gN/m2) stored vegetation nitrogen
    totvegn_patch                     ::Vector{FT} = Float64[]  # (gN/m2) total vegetation nitrogen
    totn_patch                        ::Vector{FT} = Float64[]  # (gN/m2) total patch-level nitrogen

    # --- Column-level fields ---
    totvegn_col                       ::Vector{FT} = Float64[]  # (gN/m2) total vegetation nitrogen (p2c)
    totn_p2c_col                      ::Vector{FT} = Float64[]  # (gN/m2) totn_patch averaged to column

    # --- Gridcell-level ---
    seedn_grc                         ::Vector{FT} = Float64[]  # (gN/m2) pool for seeding new PFTs via dynamic landcover
end

"""
    cnveg_nitrogen_state_init!(ns, np, nc, ng; use_matrixcn=false, nrepr=NREPR)

Allocate and initialize all fields of a `CNVegNitrogenStateData` instance for
`np` patches, `nc` columns, and `ng` gridcells.
Matches the Fortran `InitAllocate`.
"""
function cnveg_nitrogen_state_init!(ns::CNVegNitrogenStateData, np::Int, nc::Int, ng::Int;
                                     use_matrixcn::Bool=false, nrepr::Int=NREPR)
    # --- Patch-level 2D (patch x nrepr) ---
    ns.reproductiven_patch          = fill(NaN, np, nrepr)
    ns.reproductiven_storage_patch  = fill(NaN, np, nrepr)
    ns.reproductiven_xfer_patch     = fill(NaN, np, nrepr)

    # --- Patch-level 1D core pools ---
    ns.leafn_patch                  = fill(NaN, np)
    ns.leafn_storage_patch          = fill(NaN, np)
    ns.leafn_xfer_patch             = fill(NaN, np)
    ns.leafn_storage_xfer_acc_patch = fill(NaN, np)
    ns.storage_ndemand_patch        = fill(NaN, np)
    ns.frootn_patch                 = fill(NaN, np)
    ns.frootn_storage_patch         = fill(NaN, np)
    ns.frootn_xfer_patch            = fill(NaN, np)
    ns.livestemn_patch              = fill(NaN, np)
    ns.livestemn_storage_patch      = fill(NaN, np)
    ns.livestemn_xfer_patch         = fill(NaN, np)
    ns.deadstemn_patch              = fill(NaN, np)
    ns.deadstemn_storage_patch      = fill(NaN, np)
    ns.deadstemn_xfer_patch         = fill(NaN, np)
    ns.livecrootn_patch             = fill(NaN, np)
    ns.livecrootn_storage_patch     = fill(NaN, np)
    ns.livecrootn_xfer_patch        = fill(NaN, np)
    ns.deadcrootn_patch             = fill(NaN, np)
    ns.deadcrootn_storage_patch     = fill(NaN, np)
    ns.deadcrootn_xfer_patch        = fill(NaN, np)
    ns.retransn_patch               = fill(NaN, np)
    ns.npool_patch                  = fill(NaN, np)
    ns.ntrunc_patch                 = fill(NaN, np)
    ns.cropseedn_deficit_patch      = fill(NaN, np)

    # --- Summary fields ---
    ns.dispvegn_patch               = fill(NaN, np)
    ns.storvegn_patch               = fill(NaN, np)
    ns.totvegn_patch                = fill(NaN, np)
    ns.totn_patch                   = fill(NaN, np)

    # --- Column-level ---
    ns.totvegn_col                  = fill(NaN, nc)
    ns.totn_p2c_col                 = fill(NaN, nc)

    # --- Gridcell-level ---
    ns.seedn_grc                    = fill(NaN, ng)

    # --- Matrix CN fields (conditional) ---
    if use_matrixcn
        ns.matrix_cap_repron_patch             = fill(NaN, np)
        ns.matrix_cap_repron_storage_patch     = fill(NaN, np)
        ns.matrix_cap_repron_xfer_patch        = fill(NaN, np)
        ns.matrix_cap_leafn_patch              = fill(NaN, np)
        ns.matrix_cap_leafn_storage_patch      = fill(NaN, np)
        ns.matrix_cap_leafn_xfer_patch         = fill(NaN, np)
        ns.matrix_cap_frootn_patch             = fill(NaN, np)
        ns.matrix_cap_frootn_storage_patch     = fill(NaN, np)
        ns.matrix_cap_frootn_xfer_patch        = fill(NaN, np)
        ns.matrix_cap_livestemn_patch          = fill(NaN, np)
        ns.matrix_cap_livestemn_storage_patch  = fill(NaN, np)
        ns.matrix_cap_livestemn_xfer_patch     = fill(NaN, np)
        ns.matrix_cap_deadstemn_patch          = fill(NaN, np)
        ns.matrix_cap_deadstemn_storage_patch  = fill(NaN, np)
        ns.matrix_cap_deadstemn_xfer_patch     = fill(NaN, np)
        ns.matrix_cap_livecrootn_patch         = fill(NaN, np)
        ns.matrix_cap_livecrootn_storage_patch = fill(NaN, np)
        ns.matrix_cap_livecrootn_xfer_patch    = fill(NaN, np)
        ns.matrix_cap_deadcrootn_patch         = fill(NaN, np)
        ns.matrix_cap_deadcrootn_storage_patch = fill(NaN, np)
        ns.matrix_cap_deadcrootn_xfer_patch    = fill(NaN, np)

        # Initial pool size for SASU
        ns.leafn0_patch                = fill(NaN, np)
        ns.leafn0_storage_patch        = fill(NaN, np)
        ns.leafn0_xfer_patch           = fill(NaN, np)
        ns.frootn0_patch               = fill(NaN, np)
        ns.frootn0_storage_patch       = fill(NaN, np)
        ns.frootn0_xfer_patch          = fill(NaN, np)
        ns.livestemn0_patch            = fill(NaN, np)
        ns.livestemn0_storage_patch    = fill(NaN, np)
        ns.livestemn0_xfer_patch       = fill(NaN, np)
        ns.deadstemn0_patch            = fill(NaN, np)
        ns.deadstemn0_storage_patch    = fill(NaN, np)
        ns.deadstemn0_xfer_patch       = fill(NaN, np)
        ns.livecrootn0_patch           = fill(NaN, np)
        ns.livecrootn0_storage_patch   = fill(NaN, np)
        ns.livecrootn0_xfer_patch      = fill(NaN, np)
        ns.deadcrootn0_patch           = fill(NaN, np)
        ns.deadcrootn0_storage_patch   = fill(NaN, np)
        ns.deadcrootn0_xfer_patch      = fill(NaN, np)
        ns.retransn0_patch             = fill(NaN, np)
        ns.repron0_patch               = fill(NaN, np)
        ns.repron0_storage_patch       = fill(NaN, np)
        ns.repron0_xfer_patch          = fill(NaN, np)

        # SASU save fields
        ns.leafn_SASUsave_patch              = fill(NaN, np)
        ns.leafn_storage_SASUsave_patch      = fill(NaN, np)
        ns.leafn_xfer_SASUsave_patch         = fill(NaN, np)
        ns.frootn_SASUsave_patch             = fill(NaN, np)
        ns.frootn_storage_SASUsave_patch     = fill(NaN, np)
        ns.frootn_xfer_SASUsave_patch        = fill(NaN, np)
        ns.livestemn_SASUsave_patch          = fill(NaN, np)
        ns.livestemn_storage_SASUsave_patch  = fill(NaN, np)
        ns.livestemn_xfer_SASUsave_patch     = fill(NaN, np)
        ns.deadstemn_SASUsave_patch          = fill(NaN, np)
        ns.deadstemn_storage_SASUsave_patch  = fill(NaN, np)
        ns.deadstemn_xfer_SASUsave_patch     = fill(NaN, np)
        ns.livecrootn_SASUsave_patch         = fill(NaN, np)
        ns.livecrootn_storage_SASUsave_patch = fill(NaN, np)
        ns.livecrootn_xfer_SASUsave_patch    = fill(NaN, np)
        ns.deadcrootn_SASUsave_patch         = fill(NaN, np)
        ns.deadcrootn_storage_SASUsave_patch = fill(NaN, np)
        ns.deadcrootn_xfer_SASUsave_patch    = fill(NaN, np)
        ns.grainn_SASUsave_patch             = fill(NaN, np)
        ns.grainn_storage_SASUsave_patch     = fill(NaN, np)

        # Matrix N alloc accumulation
        ns.matrix_nalloc_leaf_acc_patch        = fill(NaN, np)
        ns.matrix_nalloc_leafst_acc_patch      = fill(NaN, np)
        ns.matrix_nalloc_froot_acc_patch       = fill(NaN, np)
        ns.matrix_nalloc_frootst_acc_patch     = fill(NaN, np)
        ns.matrix_nalloc_livestem_acc_patch    = fill(NaN, np)
        ns.matrix_nalloc_livestemst_acc_patch  = fill(NaN, np)
        ns.matrix_nalloc_deadstem_acc_patch    = fill(NaN, np)
        ns.matrix_nalloc_deadstemst_acc_patch  = fill(NaN, np)
        ns.matrix_nalloc_livecroot_acc_patch   = fill(NaN, np)
        ns.matrix_nalloc_livecrootst_acc_patch = fill(NaN, np)
        ns.matrix_nalloc_deadcroot_acc_patch   = fill(NaN, np)
        ns.matrix_nalloc_deadcrootst_acc_patch = fill(NaN, np)
        ns.matrix_nalloc_grain_acc_patch       = fill(NaN, np)
        ns.matrix_nalloc_grainst_acc_patch     = fill(NaN, np)

        # Matrix N transfer accumulation
        ns.matrix_ntransfer_leafst_to_leafxf_acc_patch           = fill(NaN, np)
        ns.matrix_ntransfer_leafxf_to_leaf_acc_patch             = fill(NaN, np)
        ns.matrix_ntransfer_frootst_to_frootxf_acc_patch         = fill(NaN, np)
        ns.matrix_ntransfer_frootxf_to_froot_acc_patch           = fill(NaN, np)
        ns.matrix_ntransfer_livestemst_to_livestemxf_acc_patch   = fill(NaN, np)
        ns.matrix_ntransfer_livestemxf_to_livestem_acc_patch     = fill(NaN, np)
        ns.matrix_ntransfer_deadstemst_to_deadstemxf_acc_patch   = fill(NaN, np)
        ns.matrix_ntransfer_deadstemxf_to_deadstem_acc_patch     = fill(NaN, np)
        ns.matrix_ntransfer_livecrootst_to_livecrootxf_acc_patch = fill(NaN, np)
        ns.matrix_ntransfer_livecrootxf_to_livecroot_acc_patch   = fill(NaN, np)
        ns.matrix_ntransfer_deadcrootst_to_deadcrootxf_acc_patch = fill(NaN, np)
        ns.matrix_ntransfer_deadcrootxf_to_deadcroot_acc_patch   = fill(NaN, np)
        ns.matrix_ntransfer_grainst_to_grainxf_acc_patch         = fill(NaN, np)
        ns.matrix_ntransfer_grainxf_to_grain_acc_patch           = fill(NaN, np)
        ns.matrix_ntransfer_livestem_to_deadstem_acc_patch        = fill(NaN, np)
        ns.matrix_ntransfer_livecroot_to_deadcroot_acc_patch      = fill(NaN, np)

        # Matrix N retranslocation transfer accumulation
        ns.matrix_ntransfer_retransn_to_leaf_acc_patch           = fill(NaN, np)
        ns.matrix_ntransfer_retransn_to_leafst_acc_patch         = fill(NaN, np)
        ns.matrix_ntransfer_retransn_to_froot_acc_patch          = fill(NaN, np)
        ns.matrix_ntransfer_retransn_to_frootst_acc_patch        = fill(NaN, np)
        ns.matrix_ntransfer_retransn_to_livestem_acc_patch       = fill(NaN, np)
        ns.matrix_ntransfer_retransn_to_livestemst_acc_patch     = fill(NaN, np)
        ns.matrix_ntransfer_retransn_to_deadstem_acc_patch       = fill(NaN, np)
        ns.matrix_ntransfer_retransn_to_deadstemst_acc_patch     = fill(NaN, np)
        ns.matrix_ntransfer_retransn_to_livecroot_acc_patch      = fill(NaN, np)
        ns.matrix_ntransfer_retransn_to_livecrootst_acc_patch    = fill(NaN, np)
        ns.matrix_ntransfer_retransn_to_deadcroot_acc_patch      = fill(NaN, np)
        ns.matrix_ntransfer_retransn_to_deadcrootst_acc_patch    = fill(NaN, np)
        ns.matrix_ntransfer_retransn_to_grain_acc_patch          = fill(NaN, np)
        ns.matrix_ntransfer_retransn_to_grainst_acc_patch        = fill(NaN, np)

        # Matrix N from-pool to retranslocation
        ns.matrix_ntransfer_leaf_to_retransn_acc_patch           = fill(NaN, np)
        ns.matrix_ntransfer_froot_to_retransn_acc_patch          = fill(NaN, np)
        ns.matrix_ntransfer_livestem_to_retransn_acc_patch       = fill(NaN, np)
        ns.matrix_ntransfer_livecroot_to_retransn_acc_patch      = fill(NaN, np)

        # Matrix N turnover accumulation
        ns.matrix_nturnover_leaf_acc_patch             = fill(NaN, np)
        ns.matrix_nturnover_leafst_acc_patch           = fill(NaN, np)
        ns.matrix_nturnover_leafxf_acc_patch           = fill(NaN, np)
        ns.matrix_nturnover_froot_acc_patch            = fill(NaN, np)
        ns.matrix_nturnover_frootst_acc_patch          = fill(NaN, np)
        ns.matrix_nturnover_frootxf_acc_patch          = fill(NaN, np)
        ns.matrix_nturnover_livestem_acc_patch         = fill(NaN, np)
        ns.matrix_nturnover_livestemst_acc_patch       = fill(NaN, np)
        ns.matrix_nturnover_livestemxf_acc_patch       = fill(NaN, np)
        ns.matrix_nturnover_deadstem_acc_patch         = fill(NaN, np)
        ns.matrix_nturnover_deadstemst_acc_patch       = fill(NaN, np)
        ns.matrix_nturnover_deadstemxf_acc_patch       = fill(NaN, np)
        ns.matrix_nturnover_livecroot_acc_patch        = fill(NaN, np)
        ns.matrix_nturnover_livecrootst_acc_patch      = fill(NaN, np)
        ns.matrix_nturnover_livecrootxf_acc_patch      = fill(NaN, np)
        ns.matrix_nturnover_deadcroot_acc_patch        = fill(NaN, np)
        ns.matrix_nturnover_deadcrootst_acc_patch      = fill(NaN, np)
        ns.matrix_nturnover_deadcrootxf_acc_patch      = fill(NaN, np)
        ns.matrix_nturnover_grain_acc_patch            = fill(NaN, np)
        ns.matrix_nturnover_grainst_acc_patch          = fill(NaN, np)
        ns.matrix_nturnover_grainxf_acc_patch          = fill(NaN, np)
        ns.matrix_nturnover_retransn_acc_patch         = fill(NaN, np)
    end

    return nothing
end

"""
    cnveg_nitrogen_state_set_values!(ns, mask_patch, value_patch, mask_col, value_col;
                                     use_matrixcn=false, use_crop=false)

Set nitrogen state variables for filtered patches and columns.
Ported from `cnveg_nitrogenstate_type%SetValues`.
"""
function cnveg_nitrogen_state_set_values!(ns::CNVegNitrogenStateData,
                                          mask_patch::BitVector, value_patch::Float64,
                                          mask_col::BitVector, value_col::Float64;
                                          use_matrixcn::Bool=false, use_crop::Bool=false,
                                          nrepr::Int=NREPR)
    # Patch-level
    for i in eachindex(mask_patch)
        mask_patch[i] || continue
        ns.leafn_patch[i]              = value_patch
        ns.leafn_storage_patch[i]      = value_patch
        ns.leafn_xfer_patch[i]         = value_patch
        ns.leafn_storage_xfer_acc_patch[i] = value_patch
        ns.frootn_patch[i]             = value_patch
        ns.frootn_storage_patch[i]     = value_patch
        ns.frootn_xfer_patch[i]        = value_patch
        ns.livestemn_patch[i]          = value_patch
        ns.livestemn_storage_patch[i]  = value_patch
        ns.livestemn_xfer_patch[i]     = value_patch
        ns.deadstemn_patch[i]          = value_patch
        ns.deadstemn_storage_patch[i]  = value_patch
        ns.deadstemn_xfer_patch[i]     = value_patch
        ns.livecrootn_patch[i]         = value_patch
        ns.livecrootn_storage_patch[i] = value_patch
        ns.livecrootn_xfer_patch[i]    = value_patch
        ns.deadcrootn_patch[i]         = value_patch
        ns.deadcrootn_storage_patch[i] = value_patch
        ns.deadcrootn_xfer_patch[i]    = value_patch

        if use_matrixcn && length(ns.matrix_cap_leafn_patch) > 0
            ns.matrix_cap_leafn_patch[i]              = value_patch
            ns.matrix_cap_leafn_storage_patch[i]      = value_patch
            ns.matrix_cap_leafn_xfer_patch[i]         = value_patch
            ns.matrix_cap_frootn_patch[i]             = value_patch
            ns.matrix_cap_frootn_storage_patch[i]     = value_patch
            ns.matrix_cap_frootn_xfer_patch[i]        = value_patch
            ns.matrix_cap_livestemn_patch[i]          = value_patch
            ns.matrix_cap_livestemn_storage_patch[i]  = value_patch
            ns.matrix_cap_livestemn_xfer_patch[i]     = value_patch
            ns.matrix_cap_deadstemn_patch[i]          = value_patch
            ns.matrix_cap_deadstemn_storage_patch[i]  = value_patch
            ns.matrix_cap_deadstemn_xfer_patch[i]     = value_patch
            ns.matrix_cap_livecrootn_patch[i]         = value_patch
            ns.matrix_cap_livecrootn_storage_patch[i] = value_patch
            ns.matrix_cap_livecrootn_xfer_patch[i]    = value_patch
            ns.matrix_cap_deadcrootn_patch[i]         = value_patch
            ns.matrix_cap_deadcrootn_storage_patch[i] = value_patch
            ns.matrix_cap_deadcrootn_xfer_patch[i]    = value_patch

            ns.leafn0_patch[i]              = value_patch
            ns.leafn0_storage_patch[i]      = value_patch
            ns.leafn0_xfer_patch[i]         = value_patch
            ns.frootn0_patch[i]             = value_patch
            ns.frootn0_storage_patch[i]     = value_patch
            ns.frootn0_xfer_patch[i]        = value_patch
            ns.livestemn0_patch[i]          = value_patch
            ns.livestemn0_storage_patch[i]  = value_patch
            ns.livestemn0_xfer_patch[i]     = value_patch
            ns.deadstemn0_patch[i]          = value_patch
            ns.deadstemn0_storage_patch[i]  = value_patch
            ns.deadstemn0_xfer_patch[i]     = value_patch
            ns.livecrootn0_patch[i]         = value_patch
            ns.livecrootn0_storage_patch[i] = value_patch
            ns.livecrootn0_xfer_patch[i]    = value_patch
            ns.deadcrootn0_patch[i]         = value_patch
            ns.deadcrootn0_storage_patch[i] = value_patch
            ns.deadcrootn0_xfer_patch[i]    = value_patch
            ns.retransn0_patch[i]           = value_patch
            if use_crop
                ns.repron0_patch[i]            = value_patch
                ns.repron0_storage_patch[i]    = value_patch
                ns.repron0_xfer_patch[i]       = value_patch
            end

            ns.matrix_nalloc_leaf_acc_patch[i]        = value_patch
            ns.matrix_nalloc_leafst_acc_patch[i]      = value_patch
            ns.matrix_nalloc_froot_acc_patch[i]       = value_patch
            ns.matrix_nalloc_frootst_acc_patch[i]     = value_patch
            ns.matrix_nalloc_livestem_acc_patch[i]    = value_patch
            ns.matrix_nalloc_livestemst_acc_patch[i]  = value_patch
            ns.matrix_nalloc_deadstem_acc_patch[i]    = value_patch
            ns.matrix_nalloc_deadstemst_acc_patch[i]  = value_patch
            ns.matrix_nalloc_livecroot_acc_patch[i]   = value_patch
            ns.matrix_nalloc_livecrootst_acc_patch[i] = value_patch
            ns.matrix_nalloc_deadcroot_acc_patch[i]   = value_patch
            ns.matrix_nalloc_deadcrootst_acc_patch[i] = value_patch
            ns.matrix_nalloc_grain_acc_patch[i]       = value_patch
            ns.matrix_nalloc_grainst_acc_patch[i]     = value_patch

            ns.matrix_ntransfer_leafst_to_leafxf_acc_patch[i]           = value_patch
            ns.matrix_ntransfer_leafxf_to_leaf_acc_patch[i]             = value_patch
            ns.matrix_ntransfer_frootst_to_frootxf_acc_patch[i]         = value_patch
            ns.matrix_ntransfer_frootxf_to_froot_acc_patch[i]           = value_patch
            ns.matrix_ntransfer_livestemst_to_livestemxf_acc_patch[i]   = value_patch
            ns.matrix_ntransfer_livestemxf_to_livestem_acc_patch[i]     = value_patch
            ns.matrix_ntransfer_deadstemst_to_deadstemxf_acc_patch[i]   = value_patch
            ns.matrix_ntransfer_deadstemxf_to_deadstem_acc_patch[i]     = value_patch
            ns.matrix_ntransfer_livecrootst_to_livecrootxf_acc_patch[i] = value_patch
            ns.matrix_ntransfer_livecrootxf_to_livecroot_acc_patch[i]   = value_patch
            ns.matrix_ntransfer_deadcrootst_to_deadcrootxf_acc_patch[i] = value_patch
            ns.matrix_ntransfer_deadcrootxf_to_deadcroot_acc_patch[i]   = value_patch
            if use_crop
                ns.matrix_ntransfer_grainst_to_grainxf_acc_patch[i]    = value_patch
                ns.matrix_ntransfer_grainxf_to_grain_acc_patch[i]      = value_patch
            end
            ns.matrix_ntransfer_livestem_to_deadstem_acc_patch[i]       = value_patch
            ns.matrix_ntransfer_livecroot_to_deadcroot_acc_patch[i]     = value_patch

            ns.matrix_ntransfer_retransn_to_leaf_acc_patch[i]           = value_patch
            ns.matrix_ntransfer_retransn_to_leafst_acc_patch[i]         = value_patch
            ns.matrix_ntransfer_retransn_to_froot_acc_patch[i]          = value_patch
            ns.matrix_ntransfer_retransn_to_frootst_acc_patch[i]        = value_patch
            ns.matrix_ntransfer_retransn_to_livestem_acc_patch[i]       = value_patch
            ns.matrix_ntransfer_retransn_to_livestemst_acc_patch[i]     = value_patch
            ns.matrix_ntransfer_retransn_to_deadstem_acc_patch[i]       = value_patch
            ns.matrix_ntransfer_retransn_to_deadstemst_acc_patch[i]     = value_patch
            ns.matrix_ntransfer_retransn_to_livecroot_acc_patch[i]      = value_patch
            ns.matrix_ntransfer_retransn_to_livecrootst_acc_patch[i]    = value_patch
            ns.matrix_ntransfer_retransn_to_deadcroot_acc_patch[i]      = value_patch
            ns.matrix_ntransfer_retransn_to_deadcrootst_acc_patch[i]    = value_patch
            ns.matrix_ntransfer_retransn_to_grain_acc_patch[i]          = value_patch
            ns.matrix_ntransfer_retransn_to_grainst_acc_patch[i]        = value_patch

            ns.matrix_ntransfer_leaf_to_retransn_acc_patch[i]           = value_patch
            ns.matrix_ntransfer_froot_to_retransn_acc_patch[i]          = value_patch
            ns.matrix_ntransfer_livestem_to_retransn_acc_patch[i]       = value_patch
            ns.matrix_ntransfer_livecroot_to_retransn_acc_patch[i]      = value_patch

            ns.matrix_nturnover_leaf_acc_patch[i]             = value_patch
            ns.matrix_nturnover_leafst_acc_patch[i]           = value_patch
            ns.matrix_nturnover_leafxf_acc_patch[i]           = value_patch
            ns.matrix_nturnover_froot_acc_patch[i]            = value_patch
            ns.matrix_nturnover_frootst_acc_patch[i]          = value_patch
            ns.matrix_nturnover_frootxf_acc_patch[i]          = value_patch
            ns.matrix_nturnover_livestem_acc_patch[i]         = value_patch
            ns.matrix_nturnover_livestemst_acc_patch[i]       = value_patch
            ns.matrix_nturnover_livestemxf_acc_patch[i]       = value_patch
            ns.matrix_nturnover_deadstem_acc_patch[i]         = value_patch
            ns.matrix_nturnover_deadstemst_acc_patch[i]       = value_patch
            ns.matrix_nturnover_deadstemxf_acc_patch[i]       = value_patch
            ns.matrix_nturnover_livecroot_acc_patch[i]        = value_patch
            ns.matrix_nturnover_livecrootst_acc_patch[i]      = value_patch
            ns.matrix_nturnover_livecrootxf_acc_patch[i]      = value_patch
            ns.matrix_nturnover_deadcroot_acc_patch[i]        = value_patch
            ns.matrix_nturnover_deadcrootst_acc_patch[i]      = value_patch
            ns.matrix_nturnover_deadcrootxf_acc_patch[i]      = value_patch
            ns.matrix_nturnover_retransn_acc_patch[i]         = value_patch
            if use_crop
                ns.matrix_nturnover_grain_acc_patch[i]        = value_patch
                ns.matrix_nturnover_grainst_acc_patch[i]      = value_patch
                ns.matrix_nturnover_grainxf_acc_patch[i]      = value_patch
            end
        end

        ns.retransn_patch[i]           = value_patch
        ns.npool_patch[i]              = value_patch
        ns.ntrunc_patch[i]             = value_patch
        ns.dispvegn_patch[i]           = value_patch
        ns.storvegn_patch[i]           = value_patch
        ns.totvegn_patch[i]            = value_patch
        ns.totn_patch[i]               = value_patch
    end

    # Reproductive N (2D) for crop
    if use_crop
        for k in 1:nrepr
            for i in eachindex(mask_patch)
                mask_patch[i] || continue
                ns.reproductiven_patch[i, k]         = value_patch
                ns.reproductiven_storage_patch[i, k]  = value_patch
                ns.reproductiven_xfer_patch[i, k]     = value_patch
            end
        end
        for i in eachindex(mask_patch)
            mask_patch[i] || continue
            ns.cropseedn_deficit_patch[i] = value_patch
        end
        if use_matrixcn && length(ns.matrix_cap_repron_patch) > 0
            for i in eachindex(mask_patch)
                mask_patch[i] || continue
                ns.matrix_cap_repron_patch[i]         = value_patch
                ns.matrix_cap_repron_storage_patch[i]  = value_patch
                ns.matrix_cap_repron_xfer_patch[i]     = value_patch
            end
        end
    end

    # Column-level
    for i in eachindex(mask_col)
        mask_col[i] || continue
        ns.totvegn_col[i]    = value_col
        ns.totn_p2c_col[i]   = value_col
    end

    return nothing
end

"""
    cnveg_nitrogen_state_zero_dwt!(ns, bounds_patch)

Initialize variables needed for dynamic land use.
Ported from `cnveg_nitrogenstate_type%ZeroDwt`.
"""
function cnveg_nitrogen_state_zero_dwt!(ns::CNVegNitrogenStateData,
                                         bounds_patch::UnitRange{Int})
    for p in bounds_patch
        ns.dispvegn_patch[p] = 0.0
        ns.storvegn_patch[p] = 0.0
        ns.totvegn_patch[p]  = 0.0
        ns.totn_patch[p]     = 0.0
    end
    return nothing
end

"""
    cnveg_nitrogen_state_summary!(ns, mask_patch, bounds_patch;
                                   use_crop=false, patch_itype=nothing,
                                   npcropmin=0, nrepr=NREPR)

Perform patch-level nitrogen summary calculations.
Ported from `cnveg_nitrogenstate_type%Summary_nitrogenstate`.

Note: column-level p2c averaging is stubbed out pending subgridAveMod port.
"""
function cnveg_nitrogen_state_summary!(ns::CNVegNitrogenStateData,
                                        mask_patch::BitVector,
                                        bounds_patch::UnitRange{Int};
                                        use_crop::Bool=false,
                                        patch_itype::Union{Vector{Int},Nothing}=nothing,
                                        npcropmin::Int=0,
                                        nrepr::Int=NREPR)
    for p in bounds_patch
        mask_patch[p] || continue

        # displayed vegetation nitrogen, excluding storage (DISPVEGN)
        ns.dispvegn_patch[p] =
            ns.leafn_patch[p]      +
            ns.frootn_patch[p]     +
            ns.livestemn_patch[p]  +
            ns.deadstemn_patch[p]  +
            ns.livecrootn_patch[p] +
            ns.deadcrootn_patch[p]

        # stored vegetation nitrogen, including retranslocated N pool (STORVEGN)
        ns.storvegn_patch[p] =
            ns.leafn_storage_patch[p]      +
            ns.frootn_storage_patch[p]     +
            ns.livestemn_storage_patch[p]  +
            ns.deadstemn_storage_patch[p]  +
            ns.livecrootn_storage_patch[p] +
            ns.deadcrootn_storage_patch[p] +
            ns.leafn_xfer_patch[p]         +
            ns.frootn_xfer_patch[p]        +
            ns.livestemn_xfer_patch[p]     +
            ns.deadstemn_xfer_patch[p]     +
            ns.livecrootn_xfer_patch[p]    +
            ns.deadcrootn_xfer_patch[p]    +
            ns.npool_patch[p]              +
            ns.retransn_patch[p]

        if use_crop && patch_itype !== nothing && patch_itype[p] >= npcropmin
            for k in 1:nrepr
                ns.dispvegn_patch[p] +=
                    ns.reproductiven_patch[p, k]

                ns.storvegn_patch[p] +=
                    ns.reproductiven_storage_patch[p, k] +
                    ns.reproductiven_xfer_patch[p, k]
            end

            ns.storvegn_patch[p] +=
                ns.cropseedn_deficit_patch[p]
        end

        # total vegetation nitrogen (TOTVEGN)
        ns.totvegn_patch[p] =
            ns.dispvegn_patch[p] +
            ns.storvegn_patch[p]

        # total patch-level nitrogen (add ntrunc)
        ns.totn_patch[p] =
            ns.totvegn_patch[p] +
            ns.ntrunc_patch[p]
    end

    # Column-level p2c averaging is stubbed pending subgridAveMod port
    return nothing
end

"""
    cnveg_nitrogen_state_init_cold!(ns, bounds_patch; use_matrixcn=false, use_crop=false)

Initialize cold-start conditions for CN vegetation nitrogen state variables.
Simplified version that initializes all patches as soil/crop type.

Ported from `cnveg_nitrogenstate_type%InitCold`.
"""
function cnveg_nitrogen_state_init_cold!(ns::CNVegNitrogenStateData,
                                          bounds_patch::UnitRange{Int};
                                          use_matrixcn::Bool=false,
                                          use_crop::Bool=false)
    for p in bounds_patch
        ns.leafn_patch[p]         = 0.0
        ns.leafn_storage_patch[p] = 0.0
        ns.leafn_xfer_patch[p]    = 0.0
        ns.leafn_storage_xfer_acc_patch[p] = 0.0
        ns.storage_ndemand_patch[p]        = 0.0

        ns.frootn_patch[p]            = 0.0
        ns.frootn_storage_patch[p]    = 0.0
        ns.frootn_xfer_patch[p]       = 0.0

        ns.livestemn_patch[p]         = 0.0
        ns.livestemn_storage_patch[p] = 0.0
        ns.livestemn_xfer_patch[p]    = 0.0

        ns.deadstemn_patch[p]         = 0.0
        ns.deadstemn_storage_patch[p] = 0.0
        ns.deadstemn_xfer_patch[p]    = 0.0

        ns.livecrootn_patch[p]         = 0.0
        ns.livecrootn_storage_patch[p] = 0.0
        ns.livecrootn_xfer_patch[p]    = 0.0

        ns.deadcrootn_patch[p]         = 0.0
        ns.deadcrootn_storage_patch[p] = 0.0
        ns.deadcrootn_xfer_patch[p]    = 0.0

        ns.retransn_patch[p]           = 0.0
        ns.npool_patch[p]              = 0.0
        ns.ntrunc_patch[p]             = 0.0
        ns.dispvegn_patch[p]           = 0.0
        ns.storvegn_patch[p]           = 0.0
        ns.totvegn_patch[p]            = 0.0
        ns.totn_patch[p]               = 0.0

        if use_crop
            for k in 1:size(ns.reproductiven_patch, 2)
                ns.reproductiven_patch[p, k]         = 0.0
                ns.reproductiven_storage_patch[p, k] = 0.0
                ns.reproductiven_xfer_patch[p, k]    = 0.0
            end
            ns.cropseedn_deficit_patch[p] = 0.0
            if use_matrixcn && length(ns.matrix_cap_repron_patch) > 0
                ns.matrix_cap_repron_patch[p]         = 0.0
                ns.matrix_cap_repron_storage_patch[p] = 0.0
                ns.matrix_cap_repron_xfer_patch[p]    = 0.0
            end
        end

        if use_matrixcn && length(ns.matrix_cap_leafn_patch) > 0
            ns.matrix_cap_leafn_patch[p]          = 0.0
            ns.matrix_cap_leafn_storage_patch[p]  = 0.0
            ns.matrix_cap_leafn_xfer_patch[p]     = 0.0
            ns.matrix_cap_frootn_patch[p]         = 0.0
            ns.matrix_cap_frootn_storage_patch[p] = 0.0
            ns.matrix_cap_frootn_xfer_patch[p]    = 0.0
            ns.matrix_cap_livestemn_patch[p]          = 0.0
            ns.matrix_cap_livestemn_storage_patch[p]  = 0.0
            ns.matrix_cap_livestemn_xfer_patch[p]     = 0.0
            ns.matrix_cap_deadstemn_patch[p]          = 0.0
            ns.matrix_cap_deadstemn_storage_patch[p]  = 0.0
            ns.matrix_cap_deadstemn_xfer_patch[p]     = 0.0
            ns.matrix_cap_livecrootn_patch[p]         = 0.0
            ns.matrix_cap_livecrootn_storage_patch[p] = 0.0
            ns.matrix_cap_livecrootn_xfer_patch[p]    = 0.0
            ns.matrix_cap_deadcrootn_patch[p]         = 0.0
            ns.matrix_cap_deadcrootn_storage_patch[p] = 0.0
            ns.matrix_cap_deadcrootn_xfer_patch[p]    = 0.0

            # Initial pool sizes for SASU
            ns.leafn0_patch[p]              = 1.0e-30
            ns.leafn0_storage_patch[p]      = 1.0e-30
            ns.leafn0_xfer_patch[p]         = 1.0e-30
            ns.frootn0_patch[p]             = 1.0e-30
            ns.frootn0_storage_patch[p]     = 1.0e-30
            ns.frootn0_xfer_patch[p]        = 1.0e-30
            ns.livestemn0_patch[p]          = 1.0e-30
            ns.livestemn0_storage_patch[p]  = 1.0e-30
            ns.livestemn0_xfer_patch[p]     = 1.0e-30
            ns.deadstemn0_patch[p]          = 1.0e-30
            ns.deadstemn0_storage_patch[p]  = 1.0e-30
            ns.deadstemn0_xfer_patch[p]     = 1.0e-30
            ns.livecrootn0_patch[p]         = 1.0e-30
            ns.livecrootn0_storage_patch[p] = 1.0e-30
            ns.livecrootn0_xfer_patch[p]    = 1.0e-30
            ns.deadcrootn0_patch[p]         = 1.0e-30
            ns.deadcrootn0_storage_patch[p] = 1.0e-30
            ns.deadcrootn0_xfer_patch[p]    = 1.0e-30
            ns.retransn0_patch[p]           = 1.0e-30
            ns.repron0_patch[p]             = 1.0e-30
            ns.repron0_storage_patch[p]     = 1.0e-30
            ns.repron0_xfer_patch[p]        = 1.0e-30

            # SASU save fields
            ns.leafn_SASUsave_patch[p]              = 0.0
            ns.leafn_storage_SASUsave_patch[p]      = 0.0
            ns.leafn_xfer_SASUsave_patch[p]         = 0.0
            ns.frootn_SASUsave_patch[p]             = 0.0
            ns.frootn_storage_SASUsave_patch[p]     = 0.0
            ns.frootn_xfer_SASUsave_patch[p]        = 0.0
            ns.livestemn_SASUsave_patch[p]          = 0.0
            ns.livestemn_storage_SASUsave_patch[p]  = 0.0
            ns.livestemn_xfer_SASUsave_patch[p]     = 0.0
            ns.deadstemn_SASUsave_patch[p]          = 0.0
            ns.deadstemn_storage_SASUsave_patch[p]  = 0.0
            ns.deadstemn_xfer_SASUsave_patch[p]     = 0.0
            ns.livecrootn_SASUsave_patch[p]         = 0.0
            ns.livecrootn_storage_SASUsave_patch[p] = 0.0
            ns.livecrootn_xfer_SASUsave_patch[p]    = 0.0
            ns.deadcrootn_SASUsave_patch[p]         = 0.0
            ns.deadcrootn_storage_SASUsave_patch[p] = 0.0
            ns.deadcrootn_xfer_SASUsave_patch[p]    = 0.0
            ns.grainn_SASUsave_patch[p]             = 0.0
            ns.grainn_storage_SASUsave_patch[p]     = 0.0

            # Matrix accumulation fields
            ns.matrix_nalloc_leaf_acc_patch[p]        = 0.0
            ns.matrix_nalloc_leafst_acc_patch[p]      = 0.0
            ns.matrix_nalloc_froot_acc_patch[p]       = 0.0
            ns.matrix_nalloc_frootst_acc_patch[p]     = 0.0
            ns.matrix_nalloc_livestem_acc_patch[p]    = 0.0
            ns.matrix_nalloc_livestemst_acc_patch[p]  = 0.0
            ns.matrix_nalloc_deadstem_acc_patch[p]    = 0.0
            ns.matrix_nalloc_deadstemst_acc_patch[p]  = 0.0
            ns.matrix_nalloc_livecroot_acc_patch[p]   = 0.0
            ns.matrix_nalloc_livecrootst_acc_patch[p] = 0.0
            ns.matrix_nalloc_deadcroot_acc_patch[p]   = 0.0
            ns.matrix_nalloc_deadcrootst_acc_patch[p] = 0.0
            ns.matrix_nalloc_grain_acc_patch[p]       = 0.0
            ns.matrix_nalloc_grainst_acc_patch[p]     = 0.0

            ns.matrix_ntransfer_leafst_to_leafxf_acc_patch[p]           = 0.0
            ns.matrix_ntransfer_leafxf_to_leaf_acc_patch[p]             = 0.0
            ns.matrix_ntransfer_frootst_to_frootxf_acc_patch[p]         = 0.0
            ns.matrix_ntransfer_frootxf_to_froot_acc_patch[p]           = 0.0
            ns.matrix_ntransfer_livestemst_to_livestemxf_acc_patch[p]   = 0.0
            ns.matrix_ntransfer_livestemxf_to_livestem_acc_patch[p]     = 0.0
            ns.matrix_ntransfer_deadstemst_to_deadstemxf_acc_patch[p]   = 0.0
            ns.matrix_ntransfer_deadstemxf_to_deadstem_acc_patch[p]     = 0.0
            ns.matrix_ntransfer_livecrootst_to_livecrootxf_acc_patch[p] = 0.0
            ns.matrix_ntransfer_livecrootxf_to_livecroot_acc_patch[p]   = 0.0
            ns.matrix_ntransfer_deadcrootst_to_deadcrootxf_acc_patch[p] = 0.0
            ns.matrix_ntransfer_deadcrootxf_to_deadcroot_acc_patch[p]   = 0.0
            ns.matrix_ntransfer_grainst_to_grainxf_acc_patch[p]         = 0.0
            ns.matrix_ntransfer_grainxf_to_grain_acc_patch[p]           = 0.0
            ns.matrix_ntransfer_livestem_to_deadstem_acc_patch[p]       = 0.0
            ns.matrix_ntransfer_livecroot_to_deadcroot_acc_patch[p]     = 0.0

            ns.matrix_ntransfer_retransn_to_leaf_acc_patch[p]           = 0.0
            ns.matrix_ntransfer_retransn_to_leafst_acc_patch[p]         = 0.0
            ns.matrix_ntransfer_retransn_to_froot_acc_patch[p]          = 0.0
            ns.matrix_ntransfer_retransn_to_frootst_acc_patch[p]        = 0.0
            ns.matrix_ntransfer_retransn_to_livestem_acc_patch[p]       = 0.0
            ns.matrix_ntransfer_retransn_to_livestemst_acc_patch[p]     = 0.0
            ns.matrix_ntransfer_retransn_to_deadstem_acc_patch[p]       = 0.0
            ns.matrix_ntransfer_retransn_to_deadstemst_acc_patch[p]     = 0.0
            ns.matrix_ntransfer_retransn_to_livecroot_acc_patch[p]      = 0.0
            ns.matrix_ntransfer_retransn_to_livecrootst_acc_patch[p]    = 0.0
            ns.matrix_ntransfer_retransn_to_deadcroot_acc_patch[p]      = 0.0
            ns.matrix_ntransfer_retransn_to_deadcrootst_acc_patch[p]    = 0.0
            ns.matrix_ntransfer_retransn_to_grain_acc_patch[p]          = 0.0
            ns.matrix_ntransfer_retransn_to_grainst_acc_patch[p]        = 0.0

            ns.matrix_ntransfer_leaf_to_retransn_acc_patch[p]           = 0.0
            ns.matrix_ntransfer_froot_to_retransn_acc_patch[p]          = 0.0
            ns.matrix_ntransfer_livestem_to_retransn_acc_patch[p]       = 0.0
            ns.matrix_ntransfer_livecroot_to_retransn_acc_patch[p]      = 0.0

            ns.matrix_nturnover_leaf_acc_patch[p]             = 0.0
            ns.matrix_nturnover_leafst_acc_patch[p]           = 0.0
            ns.matrix_nturnover_leafxf_acc_patch[p]           = 0.0
            ns.matrix_nturnover_froot_acc_patch[p]            = 0.0
            ns.matrix_nturnover_frootst_acc_patch[p]          = 0.0
            ns.matrix_nturnover_frootxf_acc_patch[p]          = 0.0
            ns.matrix_nturnover_livestem_acc_patch[p]         = 0.0
            ns.matrix_nturnover_livestemst_acc_patch[p]       = 0.0
            ns.matrix_nturnover_livestemxf_acc_patch[p]       = 0.0
            ns.matrix_nturnover_deadstem_acc_patch[p]         = 0.0
            ns.matrix_nturnover_deadstemst_acc_patch[p]       = 0.0
            ns.matrix_nturnover_deadstemxf_acc_patch[p]       = 0.0
            ns.matrix_nturnover_livecroot_acc_patch[p]        = 0.0
            ns.matrix_nturnover_livecrootst_acc_patch[p]      = 0.0
            ns.matrix_nturnover_livecrootxf_acc_patch[p]      = 0.0
            ns.matrix_nturnover_deadcroot_acc_patch[p]        = 0.0
            ns.matrix_nturnover_deadcrootst_acc_patch[p]      = 0.0
            ns.matrix_nturnover_deadcrootxf_acc_patch[p]      = 0.0
            ns.matrix_nturnover_grain_acc_patch[p]            = 0.0
            ns.matrix_nturnover_grainst_acc_patch[p]          = 0.0
            ns.matrix_nturnover_grainxf_acc_patch[p]          = 0.0
            ns.matrix_nturnover_retransn_acc_patch[p]         = 0.0
        end
    end

    return nothing
end

# ==========================================================================
# Stubs for infrastructure-dependent subroutines
# ==========================================================================

"""
    cnveg_nitrogen_state_init_history!(ns, bounds_patch, bounds_col)

Register CN vegetation nitrogen state fields for history file output.
Stub: requires histFileMod infrastructure.
"""
function cnveg_nitrogen_state_init_history!(ns::CNVegNitrogenStateData,
                                             bounds_patch::UnitRange{Int},
                                             bounds_col::UnitRange{Int})
    return nothing
end

"""
    cnveg_nitrogen_state_restart!(ns, bounds_patch, bounds_col; flag="read")

Read/write CN vegetation nitrogen state from/to restart file.
Stub: requires NetCDF/restart infrastructure.
"""
function cnveg_nitrogen_state_restart!(ns::CNVegNitrogenStateData,
                                        bounds_patch::UnitRange{Int},
                                        bounds_col::UnitRange{Int};
                                        flag::String="read")
    return nothing
end

"""
    cnveg_nitrogen_state_dynamic_patch_adjustments!(...)

Adjust state variables when patch areas change due to dynamic landuse.
Stub: requires dynPatchStateUpdaterMod and CNVegComputeSeedMod infrastructure.
"""
function cnveg_nitrogen_state_dynamic_patch_adjustments!(ns::CNVegNitrogenStateData,
                                                          bounds_patch::UnitRange{Int})
    return nothing
end
