# FatesParameterDerivedMod.jl
# Julia port of FATES src/fates/main/FatesParameterDerivedMod.F90
#
# This module holds quantities that are STATICALLY DERIVED from the raw,
# user-settable model parameters (EDPftvarcon_inst / EDParams / SFParams).
# They are unchanging once the parameters are read in and are computed once,
# early in initialization, immediately after the FATES parameters are loaded.
#
# Derived quantities:
#   * jmax25top / tpu25top / kp25top  — per-(PFT, leafage) photosynthesis
#     constants derived from vcmax25top.
#   * branch_frac                     — per-PFT fraction of aboveground woody
#     biomass in branches (sum of the first 3 CWD fractions).
#   * damage_transitions              — per-PFT damage-class transition matrix.
#
# Translation notes:
#   * fates_r8 -> Float64.
#   * The Fortran `type(param_derived_type),public :: param_derived` module
#     singleton -> a const Ref{param_derived_type} (`ParamDerived`), mirroring
#     the EDParamsMod `EDParams` Ref pattern. `param_derived()` accesses it.
#   * allocatable real(:,:) -> Matrix{Float64}; real(:) -> Vector{Float64};
#     allocatable real(:,:,:) -> Array{Float64,3}. PFT is the LAST dimension of
#     damage_transitions (i, j, ft) exactly as Fortran (nlevdamage, nlevdamage,
#     numpft); jmax25top etc. are (numpft, nleafage).
#   * Raw params read via the merged accessors edpftvarcon_inst(), sf_params(),
#     ed_params() and the module globals nleafage[]/nlevdamage[].
#
# Standalone: nothing here is added to CLMInstances or any dual-copied struct.

# =====================================================================================
# param_derived_type — container for the derived parameters.
# =====================================================================================

"""
    param_derived_type

Container for FATES parameters statically derived from the raw parameter set.
Mirrors the Fortran `param_derived_type`. Field names preserved.

Fields:
- `jmax25top::Matrix{Float64}`  (numpft × nleafage) canopy-top maximum electron
  transport rate at 25C [umol electrons/m^2/s]
- `tpu25top::Matrix{Float64}`   (numpft × nleafage) canopy-top triose phosphate
  utilization rate at 25C [umol CO2/m^2/s]
- `kp25top::Matrix{Float64}`    (numpft × nleafage) canopy-top initial slope of
  the CO2 response curve (C4 plants) at 25C
- `branch_frac::Vector{Float64}` (numpft) fraction of aboveground woody biomass
  in branches (vs. stems) — for use in damage allometries
- `damage_transitions::Array{Float64,3}` (nlevdamage × nlevdamage × numpft)
  matrix of transition probabilities between damage classes, one per PFT
"""
Base.@kwdef mutable struct param_derived_type
    jmax25top::Matrix{Float64}          = Matrix{Float64}(undef, 0, 0)
    tpu25top::Matrix{Float64}           = Matrix{Float64}(undef, 0, 0)
    kp25top::Matrix{Float64}            = Matrix{Float64}(undef, 0, 0)
    branch_frac::Vector{Float64}        = Float64[]
    damage_transitions::Array{Float64,3} = Array{Float64,3}(undef, 0, 0, 0)
end

"""
    ParamDerived

Live module-global instance of [`param_derived_type`](@ref) (Fortran
`param_derived`), held in a `Ref` so the held object can be swapped/reset (e.g.
across tests / re-init). Access the current values with [`param_derived`](@ref).
"""
const ParamDerived = Ref{param_derived_type}(param_derived_type())

"""
    param_derived() -> param_derived_type

Convenience accessor for the live derived-parameter object.
"""
param_derived() = ParamDerived[]

# Module-level debug flag (Fortran `logical :: debug = .false.`).
const fates_param_derived_debug = Ref{Bool}(false)

# =====================================================================================
# InitAllocate — allocate the per-(PFT, leafage) and per-PFT arrays.
# =====================================================================================

"""
    InitAllocate!(this::param_derived_type, numpft::Integer)

Allocate `jmax25top`, `tpu25top`, `kp25top` (numpft × nleafage) and
`branch_frac` (numpft). Mirrors Fortran `InitAllocate`. `nleafage` is read from
the FatesInterfaceTypesMod module global.
"""
function InitAllocate!(this::param_derived_type, numpft::Integer)
    nla = nleafage[]
    this.jmax25top   = zeros(Float64, numpft, nla)
    this.tpu25top    = zeros(Float64, numpft, nla)
    this.kp25top     = zeros(Float64, numpft, nla)
    this.branch_frac = zeros(Float64, numpft)
    return nothing
end

# =====================================================================================
# InitAllocateDamageTransitions — allocate the per-PFT transition matrices.
# =====================================================================================

"""
    InitAllocateDamageTransitions!(this::param_derived_type, numpft::Integer)

Allocate `damage_transitions` (nlevdamage × nlevdamage × numpft). Mirrors
Fortran `InitAllocateDamageTransitions`. `nlevdamage` is read from the
FatesInterfaceTypesMod module global.
"""
function InitAllocateDamageTransitions!(this::param_derived_type, numpft::Integer)
    nld = nlevdamage[]
    this.damage_transitions = zeros(Float64, nld, nld, numpft)
    return nothing
end

# =====================================================================================
# Init — compute and store the derived parameters.
# =====================================================================================

"""
    Init!(this::param_derived_type, numpft::Integer)

Compute and store the statically-derived parameters. Mirrors Fortran `Init`.

`jmax25top = 1.67 * vcmax25top`, `tpu25top = 0.167 * vcmax25top`,
`kp25top = 20000 * vcmax25top` per (PFT, leafage). `branch_frac` is the sum of
the first three CWD fractions (`SF_val_CWD_frac[1:3]`). Also fills the
damage-transition matrices via [`InitDamageTransitions!`](@ref).

Reads `vcmax25top` from `edpftvarcon_inst()` and `SF_val_CWD_frac` from
`sf_params()`.
"""
function Init!(this::param_derived_type, numpft::Integer)
    vcmax25top = edpftvarcon_inst().vcmax25top
    SF_val_CWD_frac = sf_params().SF_val_CWD_frac

    InitAllocate!(this, numpft)
    InitDamageTransitions!(this, numpft)

    nla = nleafage[]
    for ft in 1:numpft
        for iage in 1:nla
            # Parameters derived from vcmax25top.
            # Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 used
            # jmax25 = 1.97 vcmax25 (Wullschleger 1993); here a factor 1.67
            # from Medlyn et al (2002) Plant, Cell and Environment 25:1167-1179.
            this.jmax25top[ft, iage] = 1.67   * vcmax25top[ft, iage]
            this.tpu25top[ft, iage]  = 0.167  * vcmax25top[ft, iage]
            this.kp25top[ft, iage]   = 20000.0 * vcmax25top[ft, iage]
        end

        # Fraction of aboveground woody biomass in branches.
        this.branch_frac[ft] = sum(@view SF_val_CWD_frac[1:3])
    end

    return nothing
end

# =====================================================================================
# InitDamageTransitions — build the per-PFT damage-class transition matrices.
# =====================================================================================

"""
    InitDamageTransitions!(this::param_derived_type, numpft::Integer)

Build the per-PFT annual damage-class transition matrices. Mirrors Fortran
`InitDamageTransitions`.

For each PFT, for each source damage class `i`:
- the diagonal (stay) probability is `1 - damage_frac`,
- the remaining `damage_frac` is split among the more-damaged classes
  (`i+1 .. nlevdamage`) in proportion to their class widths,
- the row is renormalized to sum to one.

Class widths come from `ED_val_history_damage_bin_edges` (from `ed_params()`),
extended with an upper bound of 100. `damage_frac` is read per-PFT from
`edpftvarcon_inst()`.
"""
function InitDamageTransitions!(this::param_derived_type, numpft::Integer)
    InitAllocateDamageTransitions!(this, numpft)

    nld = nlevdamage[]
    bin_edges = ed_params().ED_val_history_damage_bin_edges

    # class widths: append 100 to ED_val_history_damage_bin_edges
    damage_bin_edges_ex = Vector{Float64}(undef, nld + 1)
    for j in 1:nld
        damage_bin_edges_ex[j] = bin_edges[j]
    end
    damage_bin_edges_ex[nld + 1] = 100.0

    # widths of each damage class
    class_widths = damage_bin_edges_ex[2:(nld + 1)] .- damage_bin_edges_ex[1:nld]

    damage_frac_arr = edpftvarcon_inst().damage_frac

    for ft in 1:numpft
        damage_frac = damage_frac_arr[ft]

        for i in 1:nld
            # zero the row
            this.damage_transitions[i, :, ft] .= 0.0
            # damage rate stays the same (diagonal)
            this.damage_transitions[i, i, ft] = 1.0 - damage_frac

            if i < nld
                # fraction damaged gets split according to class width
                denom = sum(@view class_widths[(i + 1):nld])
                this.damage_transitions[i, (i + 1):nld, ft] .=
                    damage_frac .* class_widths[(i + 1):nld] ./ denom
            end

            # Make sure it sums to one — they have to go somewhere.
            rowsum = sum(@view this.damage_transitions[i, :, ft])
            this.damage_transitions[i, :, ft] ./= rowsum
        end
    end

    return nothing
end
