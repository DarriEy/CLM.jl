# FatesLitterMod.jl
# Julia port of FATES src/fates/biogeochem/FatesLitterMod.F90
#
# The litter_type holds all "un-fragmented", "un-decomposed" organic material
# no longer associated with a live plant: coarse woody debris (CWD), fine
# materials (leaves/roots), and reproductive materials (seeds). One litter_type
# is allocated per element (C, N, P, ...). Type-bound procedures -> bang
# functions dispatching on litter_type (InitAllocate!, ZeroFlux!, ...).
# Deps: FatesConstantsMod (fates_r8 -> Float64, fates_unset_r8).

# ---------------------------------------------------------------------------
# Litter-class index parameters
# ---------------------------------------------------------------------------
const ncwd = 4   # number of coarse woody debris pools (twig, s branch, l branch, trunk)
const ndcmpy = 3 # number of "decomposability" pools in fines (lignin, cellulose, labile)

const ilabile = 1     # array index for labile portion
const icellulose = 2  # array index for cellulose portion
const ilignin = 3     # array index for the lignin portion

"""
    litter_type

CWD + fine-litter (+ seed) pools by element. Prognostic state, fluxes-in
(dying trees / seed rain) and fluxes-out (fragmentation, seed decay). Allocated
per element via [`InitAllocate!`](@ref). All pools area-normalized [kg/m2] (rates
[kg/m2/day]). Allocatable arrays -> SoA Julia arrays sized at allocate time.
"""
Base.@kwdef mutable struct litter_type
    element_id::Int = fates_unset_int   # PARTEH-compliant element index

    # -- Prognostic variables (litter and coarse woody debris) ----------------
    ag_cwd::Vector{Float64} = Vector{Float64}(undef, ncwd)  # above ground cwd            [kg/m2]
    bg_cwd::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # below ground cwd (cwd x soil)
    leaf_fines::Vector{Float64} = Float64[]                 # above ground leaf litter (dcmpy)
    root_fines::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # below ground root litter (dcmpy x soil)

    seed::Vector{Float64} = Float64[]        # the seed pool (viable)   (pft) [kg/m2]
    seed_germ::Vector{Float64} = Float64[]   # the germinated seed pool (pft) [kg/m2]

    # -- Fluxes in - dying trees / seed rain (no disturbance fluxes) ----------
    ag_cwd_in::Vector{Float64} = Vector{Float64}(undef, ncwd)  # (cwd)         [kg/m2/day]
    bg_cwd_in::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # (cwd x soil)  [kg/m2/day]
    leaf_fines_in::Vector{Float64} = Float64[]                 # (dcmpy)       [kg/m2/day]
    root_fines_in::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # (dcmpy x soil) [kg/m2/day]

    seed_in_local::Vector{Float64} = Float64[]   # (pft) [kg/m2/day] (from local sources)
    seed_in_extern::Vector{Float64} = Float64[]  # (pft) [kg/m2/day] (from outside cell)

    # -- Fluxes out - fragmentation, seed decay (no disturbance) --------------
    ag_cwd_frag::Vector{Float64} = Vector{Float64}(undef, ncwd)  # ag cwd fragmentation flux  [kg/m2/day]
    bg_cwd_frag::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # bg cwd fragmentation flux  [kg/m2/day]
    leaf_fines_frag::Vector{Float64} = Float64[]                 # ag fines fragmentation flux [kg/m2/day]
    root_fines_frag::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # [kg/m2/day]

    seed_decay::Vector{Float64} = Float64[]       # decay of viable seeds to litter     [kg/m2/day]
    seed_germ_decay::Vector{Float64} = Float64[]  # decay of germinated seeds to litter [kg/m2/day]
    seed_germ_in::Vector{Float64} = Float64[]     # flux from viable to germinated seed [kg/m2/day]
end

# =====================================================================================

"""
    InitAllocate!(this::litter_type, numpft, numlevsoil, element_id)

Allocate all soil-/pft-/dcmpy-dimensioned pool arrays and set everything to the
unset sentinel (`fates_unset_r8`).
"""
function InitAllocate!(this::litter_type, numpft::Integer, numlevsoil::Integer,
                       element_id::Integer)
    this.element_id = Int(element_id)

    this.bg_cwd_in   = Matrix{Float64}(undef, ncwd, numlevsoil)
    this.bg_cwd      = Matrix{Float64}(undef, ncwd, numlevsoil)
    this.bg_cwd_frag = Matrix{Float64}(undef, ncwd, numlevsoil)

    this.leaf_fines      = Vector{Float64}(undef, ndcmpy)
    this.root_fines      = Matrix{Float64}(undef, ndcmpy, numlevsoil)
    this.leaf_fines_in   = Vector{Float64}(undef, ndcmpy)
    this.root_fines_in   = Matrix{Float64}(undef, ndcmpy, numlevsoil)
    this.leaf_fines_frag = Vector{Float64}(undef, ndcmpy)
    this.root_fines_frag = Matrix{Float64}(undef, ndcmpy, numlevsoil)

    this.seed_in_local   = Vector{Float64}(undef, numpft)
    this.seed_in_extern  = Vector{Float64}(undef, numpft)
    this.seed            = Vector{Float64}(undef, numpft)
    this.seed_germ       = Vector{Float64}(undef, numpft)
    this.seed_germ_in    = Vector{Float64}(undef, numpft)
    this.seed_germ_decay = Vector{Float64}(undef, numpft)
    this.seed_decay      = Vector{Float64}(undef, numpft)

    # Initialize everything to a nonsense flag
    this.ag_cwd            .= fates_unset_r8
    this.bg_cwd            .= fates_unset_r8
    this.leaf_fines        .= fates_unset_r8
    this.root_fines        .= fates_unset_r8
    this.seed              .= fates_unset_r8
    this.seed_germ         .= fates_unset_r8

    this.ag_cwd_in         .= fates_unset_r8
    this.bg_cwd_in         .= fates_unset_r8
    this.leaf_fines_in     .= fates_unset_r8
    this.root_fines_in     .= fates_unset_r8
    this.seed_in_local     .= fates_unset_r8
    this.seed_in_extern    .= fates_unset_r8

    this.ag_cwd_frag       .= fates_unset_r8
    this.bg_cwd_frag       .= fates_unset_r8
    this.leaf_fines_frag   .= fates_unset_r8
    this.root_fines_frag   .= fates_unset_r8

    this.seed_decay        .= fates_unset_r8
    this.seed_germ_decay   .= fates_unset_r8
    this.seed_germ_in      .= fates_unset_r8

    return nothing
end

# =====================================================================================

"""
    InitConditions!(this::litter_type, init_leaf_fines, init_root_fines,
                    init_ag_cwd, init_bg_cwd, init_seed, init_seed_germ)

Uniform-value initialization of the prognostic pools (cold-starts / pre-fusion
zeroing). Does not set per-layer or per-dcmpy values individually.
"""
function InitConditions!(this::litter_type, init_leaf_fines::Real,
                         init_root_fines::Real, init_ag_cwd::Real,
                         init_bg_cwd::Real, init_seed::Real, init_seed_germ::Real)
    this.ag_cwd      .= init_ag_cwd
    this.bg_cwd      .= init_bg_cwd
    this.leaf_fines  .= init_leaf_fines
    this.root_fines  .= init_root_fines
    this.seed        .= init_seed
    this.seed_germ   .= init_seed_germ
    return nothing
end

# =====================================================================================

"""
    DeallocateLitt!(this::litter_type)

Release all pool arrays (Fortran `deallocate`). In Julia we simply reset to
empty arrays so the storage is reclaimable.
"""
function DeallocateLitt!(this::litter_type)
    this.bg_cwd          = Matrix{Float64}(undef, 0, 0)
    this.leaf_fines      = Float64[]
    this.root_fines      = Matrix{Float64}(undef, 0, 0)
    this.seed            = Float64[]
    this.seed_germ       = Float64[]

    this.bg_cwd_in       = Matrix{Float64}(undef, 0, 0)
    this.leaf_fines_in   = Float64[]
    this.root_fines_in   = Matrix{Float64}(undef, 0, 0)
    this.seed_in_local   = Float64[]
    this.seed_in_extern  = Float64[]

    this.bg_cwd_frag     = Matrix{Float64}(undef, 0, 0)
    this.leaf_fines_frag = Float64[]
    this.root_fines_frag = Matrix{Float64}(undef, 0, 0)

    this.seed_decay      = Float64[]
    this.seed_germ_decay = Float64[]
    this.seed_germ_in    = Float64[]
    return nothing
end

# =====================================================================================

"""
    ZeroFlux!(this::litter_type)

Zero all flux-in and flux-out (fragmentation / seed) pools.
"""
function ZeroFlux!(this::litter_type)
    this.ag_cwd_in         .= 0.0
    this.bg_cwd_in         .= 0.0
    this.leaf_fines_in     .= 0.0
    this.root_fines_in     .= 0.0
    this.seed_in_local     .= 0.0
    this.seed_in_extern    .= 0.0

    this.ag_cwd_frag       .= 0.0
    this.bg_cwd_frag       .= 0.0
    this.leaf_fines_frag   .= 0.0
    this.root_fines_frag   .= 0.0

    this.seed_germ_in      .= 0.0
    this.seed_decay        .= 0.0
    this.seed_germ_decay   .= 0.0
    return nothing
end

# =====================================================================================

"""
    FuseLitter!(this::litter_type, self_area, donor_area, donor_litt)

Area-weighted fusion of `donor_litt` into `this`. All pools are area-normalized,
so the result is the area-weighted density of `this` and the donor.
"""
function FuseLitter!(this::litter_type, self_area::Real, donor_area::Real,
                     donor_litt::litter_type)
    nlevsoil = size(this.bg_cwd, 2)
    npft     = size(this.seed, 1)

    self_weight  = self_area / (donor_area + self_area)
    donor_weight = 1.0 - self_weight

    for c in 1:ncwd
        this.ag_cwd[c]      = this.ag_cwd[c]      * self_weight + donor_litt.ag_cwd[c]      * donor_weight
        this.ag_cwd_in[c]   = this.ag_cwd_in[c]   * self_weight + donor_litt.ag_cwd_in[c]   * donor_weight
        this.ag_cwd_frag[c] = this.ag_cwd_frag[c] * self_weight + donor_litt.ag_cwd_frag[c] * donor_weight
        for ilyr in 1:nlevsoil
            this.bg_cwd[c, ilyr]      = this.bg_cwd[c, ilyr]      * self_weight + donor_litt.bg_cwd[c, ilyr]      * donor_weight
            this.bg_cwd_in[c, ilyr]   = this.bg_cwd_in[c, ilyr]   * self_weight + donor_litt.bg_cwd_in[c, ilyr]   * donor_weight
            this.bg_cwd_frag[c, ilyr] = this.bg_cwd_frag[c, ilyr] * self_weight + donor_litt.bg_cwd_frag[c, ilyr] * donor_weight
        end
    end

    for pft in 1:npft
        this.seed[pft]            = this.seed[pft]            * self_weight + donor_litt.seed[pft]            * donor_weight
        this.seed_germ[pft]       = this.seed_germ[pft]       * self_weight + donor_litt.seed_germ[pft]       * donor_weight
        this.seed_in_local[pft]   = this.seed_in_local[pft]   * self_weight + donor_litt.seed_in_local[pft]   * donor_weight
        this.seed_in_extern[pft]  = this.seed_in_extern[pft]  * self_weight + donor_litt.seed_in_extern[pft]  * donor_weight
        this.seed_decay[pft]      = this.seed_decay[pft]      * self_weight + donor_litt.seed_decay[pft]      * donor_weight
        this.seed_germ_decay[pft] = this.seed_germ_decay[pft] * self_weight + donor_litt.seed_germ_decay[pft] * donor_weight
        this.seed_germ_in[pft]    = this.seed_germ_in[pft]    * self_weight + donor_litt.seed_germ_in[pft]    * donor_weight
    end

    for dcmpy in 1:ndcmpy
        this.leaf_fines[dcmpy]      = this.leaf_fines[dcmpy]      * self_weight + donor_litt.leaf_fines[dcmpy]      * donor_weight
        this.leaf_fines_in[dcmpy]   = this.leaf_fines_in[dcmpy]   * self_weight + donor_litt.leaf_fines_in[dcmpy]   * donor_weight
        this.leaf_fines_frag[dcmpy] = this.leaf_fines_frag[dcmpy] * self_weight + donor_litt.leaf_fines_frag[dcmpy] * donor_weight
        for ilyr in 1:nlevsoil
            this.root_fines[dcmpy, ilyr]      = this.root_fines[dcmpy, ilyr]      * self_weight + donor_litt.root_fines[dcmpy, ilyr]      * donor_weight
            this.root_fines_in[dcmpy, ilyr]   = this.root_fines_in[dcmpy, ilyr]   * self_weight + donor_litt.root_fines_in[dcmpy, ilyr]   * donor_weight
            this.root_fines_frag[dcmpy, ilyr] = this.root_fines_frag[dcmpy, ilyr] * self_weight + donor_litt.root_fines_frag[dcmpy, ilyr] * donor_weight
        end
    end

    return nothing
end

# =====================================================================================

"""
    CopyLitter!(this::litter_type, donor_litt::litter_type)

Copy all pools from `donor_litt` into `this`.
"""
function CopyLitter!(this::litter_type, donor_litt::litter_type)
    this.ag_cwd      .= donor_litt.ag_cwd
    this.ag_cwd_in   .= donor_litt.ag_cwd_in
    this.ag_cwd_frag .= donor_litt.ag_cwd_frag

    this.bg_cwd      .= donor_litt.bg_cwd
    this.bg_cwd_in   .= donor_litt.bg_cwd_in
    this.bg_cwd_frag .= donor_litt.bg_cwd_frag

    this.leaf_fines    .= donor_litt.leaf_fines
    this.seed          .= donor_litt.seed
    this.seed_germ     .= donor_litt.seed_germ
    this.leaf_fines_in .= donor_litt.leaf_fines_in
    this.seed_in_local .= donor_litt.seed_in_local

    this.seed_in_extern  .= donor_litt.seed_in_extern
    this.leaf_fines_frag .= donor_litt.leaf_fines_frag

    this.seed_decay      .= donor_litt.seed_decay
    this.seed_germ_decay .= donor_litt.seed_germ_decay
    this.seed_germ_in    .= donor_litt.seed_germ_in
    this.root_fines      .= donor_litt.root_fines
    this.root_fines_in   .= donor_litt.root_fines_in
    this.root_fines_frag .= donor_litt.root_fines_frag

    return nothing
end

# =====================================================================================

"""
    GetTotalLitterMass(this::litter_type) -> Float64

Total litter mass: sum of ag/bg cwd, root/leaf fines, and seed (+ germ) pools.
"""
function GetTotalLitterMass(this::litter_type)
    return sum(this.ag_cwd) +
           sum(this.bg_cwd) +
           sum(this.root_fines) +
           sum(this.leaf_fines) +
           sum(this.seed) +
           sum(this.seed_germ)
end

# =====================================================================================

"""
    adjust_SF_CWD_frac(dbh, ncwd_in, SF_val_CWD_frac, SF_val_CWD_frac_adj)

Adjust the partitioning of struct + sapw into CWD pools based on cohort `dbh`,
so small cohorts do not send biomass to fuel classes larger than their dbh
(Fosberg et al., 1971; Rothermel, 1983). Writes results into
`SF_val_CWD_frac_adj` (kept as an out-argument to match the Fortran signature).
"""
function adjust_SF_CWD_frac(dbh::Real, ncwd_in::Integer,
                            SF_val_CWD_frac::AbstractVector{<:Real},
                            SF_val_CWD_frac_adj::AbstractVector{<:Real})
    # Diameter ranges (cm) from Fosberg et al., 1971
    lb_max_diam   = 7.6  # max diameter for large branch (100 hr fuel)
    sb_max_diam   = 2.5  # max diameter for small branch (10 hr fuel)
    twig_max_diam = 0.6  # max diameter for twig (1 hr fuel)

    SF_val_CWD_frac_adj .= SF_val_CWD_frac

    if dbh > lb_max_diam
        # main stem is 1,000 hr fuel; leave partitioning unchanged
        return nothing
    elseif dbh > sb_max_diam && dbh <= lb_max_diam
        SF_val_CWD_frac_adj[ncwd_in]     = 0.0
        SF_val_CWD_frac_adj[ncwd_in-1]   = SF_val_CWD_frac[ncwd_in-1] / (1.0 - SF_val_CWD_frac[ncwd_in])
        SF_val_CWD_frac_adj[ncwd_in-2]   = SF_val_CWD_frac[ncwd_in-2] / (1.0 - SF_val_CWD_frac[ncwd_in])
        SF_val_CWD_frac_adj[ncwd_in-3]   = SF_val_CWD_frac[ncwd_in-3] / (1.0 - SF_val_CWD_frac[ncwd_in])
    elseif dbh > twig_max_diam && dbh <= sb_max_diam
        SF_val_CWD_frac_adj[ncwd_in]   = 0.0
        SF_val_CWD_frac_adj[ncwd_in-1] = 0.0
        SF_val_CWD_frac_adj[ncwd_in-2] = SF_val_CWD_frac[ncwd_in-2] /
            (1.0 - (SF_val_CWD_frac[ncwd_in] + SF_val_CWD_frac[ncwd_in-1]))
        SF_val_CWD_frac_adj[ncwd_in-3] = SF_val_CWD_frac[ncwd_in-3] /
            (1.0 - (SF_val_CWD_frac[ncwd_in] + SF_val_CWD_frac[ncwd_in-1]))
    elseif dbh <= twig_max_diam
        SF_val_CWD_frac_adj[ncwd_in]   = 0.0
        SF_val_CWD_frac_adj[ncwd_in-1] = 0.0
        SF_val_CWD_frac_adj[ncwd_in-2] = 0.0
        SF_val_CWD_frac_adj[ncwd_in-3] = sum(SF_val_CWD_frac)
    end
    return nothing
end
