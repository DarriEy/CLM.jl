# FatesDispersalMod.jl
# Julia port of FATES src/fates/main/FatesDispersalMod.F90
#
# Seed dispersal across grid cells. Holds the neighbor / neighborhood linked-list
# types, the dispersal buffer type, the probability-density dispersal kernels
# (exponential / exponential-power / log-sech), and the calendar-driven
# `IsItDispersalTime` / `GetCadenceDate` logic.
#
# Translation notes:
#   * fates_r8 -> Float64; pi_const from FatesConstantsMod.
#   * The Fortran linked-list `neighbor_type`/`neighborhood_type` use `pointer`
#     fields. In Julia these become mutable structs with nullable (Union{...,
#     Nothing}) `next_neighbor` / first/last pointers — the same singly-linked
#     traversal works. The module-global `lneighbors(:)` pointer array becomes a
#     const Ref to a Vector{neighborhood_type} (initially empty).
#   * allocatable arrays -> Vector/Matrix{Float64} or {Int}. Buffers are
#     (numpft × numgc) exactly as Fortran (pft is dim 1).
#   * The module `save` state `dispersal_date`/`dispersal_flag` -> const Refs.
#   * Raw params read via edpftvarcon_inst(); kernel mode / cadence / current
#     date read via the FatesInterfaceTypesMod module-global Refs.
#   * GetCadenceDate uses the hlm_current_* Refs; init() sets dispersal_date.
#
# Standalone: nothing here is added to CLMInstances or any dual-copied struct.

# =====================================================================================
# neighbor_type — a single grid-cell neighbor node in the dispersal linked list.
# =====================================================================================

"""
    neighbor_type

A grid-cell neighbor node (Fortran `neighbor_type`). `next_neighbor` is the
nullable link to the next neighbor in the source cell's neighborhood list.

Fields:
- `next_neighbor::Union{neighbor_type,Nothing}` — next neighbor (or `nothing`)
- `gindex::Int`                  — grid cell index
- `gc_dist::Float64`             — distance between source and this neighbor
- `density_prob::Vector{Float64}`— probability density from source, per PFT
"""
Base.@kwdef mutable struct neighbor_type
    next_neighbor::Union{neighbor_type,Nothing} = nothing
    gindex::Int                                  = 0
    gc_dist::Float64                             = 0.0
    density_prob::Vector{Float64}                = Float64[]
end

# =====================================================================================
# neighborhood_type — the linked list of neighbors for one source grid cell.
# =====================================================================================

"""
    neighborhood_type

Linked list of neighbors for a given source grid cell (Fortran
`neighborhood_type`).

Fields:
- `first_neighbor::Union{neighbor_type,Nothing}` — head of the neighbor list
- `last_neighbor::Union{neighbor_type,Nothing}`  — tail of the neighbor list
- `neighbor_count::Int`            — total neighbors near source
- `neighbor_indices::Vector{Int}`  — list of gridcell indices
"""
Base.@kwdef mutable struct neighborhood_type
    first_neighbor::Union{neighbor_type,Nothing} = nothing
    last_neighbor::Union{neighbor_type,Nothing}  = nothing
    neighbor_count::Int                           = 0
    neighbor_indices::Vector{Int}                 = Int[]
end

# =====================================================================================
# dispersal_type — the per-process seed-dispersal buffers.
# =====================================================================================

"""
    dispersal_type

Seed-dispersal buffers (Fortran `dispersal_type`). All buffers are PFT-major
(pft is dimension 1), matching Fortran.

Fields:
- `outgoing_local::Matrix{Float64}`  (numpft × numgc_local)  local outgoing seeds
- `outgoing_global::Matrix{Float64}` (numpft × numgc_global) global accumulation
- `incoming_global::Matrix{Float64}` (numpft × numgc_global) incoming-seed buffer
- `ncells_array::Vector{Int}`        (numprocs) gridcells per process per rank
- `begg_array::Vector{Int}`          (numprocs) starting gridcell index per rank
"""
Base.@kwdef mutable struct dispersal_type
    outgoing_local::Matrix{Float64}  = Matrix{Float64}(undef, 0, 0)
    outgoing_global::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    incoming_global::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    ncells_array::Vector{Int}        = Int[]
    begg_array::Vector{Int}          = Int[]
end

"""
    lneighbors

Module-global array of per-source-gridcell neighborhoods (Fortran
`type(neighborhood_type),pointer :: lneighbors(:)`). Held in a `Ref` to a
`Vector{neighborhood_type}` (initially empty).
"""
const lneighbors = Ref{Vector{neighborhood_type}}(neighborhood_type[])

# Module `save` state (Fortran module variables).
# Last cadence date in which there was a dispersal.
const dispersal_date = Ref{Int}(0)
# Have seeds been dispersed globally?
const dispersal_flag = Ref{Bool}(false)

# =====================================================================================
# init — allocate / zero the dispersal buffers and set the initial dispersal date.
# =====================================================================================

"""
    init!(this::dispersal_type, numprocs, numgc_global, numgc_local, numpft)

Allocate and zero the dispersal buffers and set the initial `dispersal_date` to
the current cadence date. Mirrors Fortran `init`. Integer arrays are filled with
`fates_unset_int`.
"""
function init!(this::dispersal_type, numprocs::Integer, numgc_global::Integer,
               numgc_local::Integer, numpft::Integer)
    this.outgoing_local  = zeros(Float64, numpft, numgc_local)
    this.outgoing_global = zeros(Float64, numpft, numgc_global)
    this.incoming_global = zeros(Float64, numpft, numgc_global)
    this.begg_array      = fill(fates_unset_int, numprocs)
    this.ncells_array    = fill(fates_unset_int, numprocs)

    # Dispersal will start at the end of the current initial date.
    dispersal_date[] = GetCadenceDate()
    return nothing
end

# =====================================================================================
# ProbabilityDensity — dispatch to the configured dispersal kernel.
# =====================================================================================

"""
    ProbabilityDensity(ipft::Integer, dist::Real) -> Float64

Return the dispersal probability density at distance `dist` for PFT `ipft`,
using the kernel selected by `fates_dispersal_kernel_mode[]`. Mirrors Fortran
`ProbabilityDensity` (returns `pd` rather than using an `out` argument).
"""
function ProbabilityDensity(ipft::Integer, dist::Real)
    mode = fates_dispersal_kernel_mode[]
    if mode == fates_dispersal_kernel_exponential
        return PD_exponential(dist, ipft)
    elseif mode == fates_dispersal_kernel_exppower
        return PD_exppower(dist, ipft)
    elseif mode == fates_dispersal_kernel_logsech
        return PD_logsech(dist, ipft)
    else
        println(fates_log(),
                " ERROR: An undefined dispersal kernel was specified: ", mode)
        fates_endrun("FatesDispersalMod: undefined dispersal kernel")
    end
end

# =====================================================================================
# Gamma function (Lanczos approximation) — avoids a SpecialFunctions.jl
# dependency, matching the codebase's local `erf` pattern in constants/varcon.jl.
# Fortran's intrinsic `gamma` is the (complete) Euler gamma function; the
# exppower kernel only ever evaluates it at positive arguments (2/b, b > 0).
# =====================================================================================

"""
    fates_gamma(x::Real) -> Float64

Complete Euler gamma function Γ(x) via the Lanczos approximation (g=7, n=9),
accurate to ~1e-15 for positive real `x`. Local replacement for Fortran's
intrinsic `gamma` so FATES dispersal carries no SpecialFunctions.jl dependency.
"""
function fates_gamma(x::Real)
    # Lanczos coefficients (g = 7).
    g = 7.0
    c = (0.99999999999980993,
         676.5203681218851,
         -1259.1392167224028,
         771.32342877765313,
         -176.61502916214059,
         12.507343278686905,
         -0.13857109526572012,
         9.9843695780195716e-6,
         1.5056327351493116e-7)
    if x < 0.5
        # Reflection formula.
        return pi_const / (sin(pi_const * x) * fates_gamma(1.0 - x))
    else
        xx = x - 1.0
        a = c[1]
        t = xx + g + 0.5
        for i in 2:9
            a += c[i] / (xx + (i - 1))
        end
        return sqrt(2 * pi_const) * t^(xx + 0.5) * exp(-t) * a
    end
end

# =====================================================================================
# Dispersal kernels.
# =====================================================================================

"""
    PD_exponential(dist::Real, ipft::Integer) -> Float64

Simple exponential decay kernel: `exp(-scale * dist)` with
`scale = seed_dispersal_pdf_scale[ipft]`. Mirrors Fortran `PD_exponential`.
"""
function PD_exponential(dist::Real, ipft::Integer)
    return exp(-edpftvarcon_inst().seed_dispersal_pdf_scale[ipft] * dist)
end

"""
    PD_exppower(dist::Real, ipft::Integer) -> Float64

Exponential-power (ExP) kernel. Mirrors Fortran `PD_exppower`:

    (b / (2π Γ(2/b))) * exp(-(dist^b)/(a^b))

with `a = seed_dispersal_pdf_scale[ipft]`, `b = seed_dispersal_pdf_shape[ipft]`.
"""
function PD_exppower(dist::Real, ipft::Integer)
    param_a = edpftvarcon_inst().seed_dispersal_pdf_scale[ipft]
    param_b = edpftvarcon_inst().seed_dispersal_pdf_shape[ipft]
    return (param_b / (2 * pi_const * fates_gamma(2 / param_b))) *
           exp(-(dist^param_b) / (param_a^param_b))
end

"""
    PD_logsech(dist::Real, ipft::Integer) -> Float64

Logistic-sech (LogS) kernel. Mirrors Fortran `PD_logsech`:

    (1 / (π^2 b dist^2)) / ((dist/a)^(1/b) + (dist/a)^(-1/b))

with `a = seed_dispersal_pdf_scale[ipft]`, `b = seed_dispersal_pdf_shape[ipft]`.
"""
function PD_logsech(dist::Real, ipft::Integer)
    param_a = edpftvarcon_inst().seed_dispersal_pdf_scale[ipft]
    param_b = edpftvarcon_inst().seed_dispersal_pdf_shape[ipft]
    return (1 / (pi_const^2 * param_b * dist^2)) /
           ((dist / param_a)^(1 / param_b) + (dist / param_a)^(-1 / param_b))
end

# =====================================================================================
# IsItDispersalTime — calendar-driven decision of whether to disperse seeds.
# =====================================================================================

"""
    IsItDispersalTime(; setdispersedflag::Union{Bool,Nothing}=nothing) -> Bool

Determine if seeds should be dispersed across gridcells. Mirrors Fortran
`IsItDispersalTime` (calendar-based for now).

Logic:
- If the global dispersal flag is set, return `true` and clear the flag
  (regardless of date — pass dispersed seeds to FATES).
- Otherwise, if the current cadence date differs from `dispersal_date`, return
  `true`; if `setdispersedflag` is `true`, also set the global flag and update
  `dispersal_date` to the current cadence date.

`setdispersedflag` is the optional Fortran argument (default treated as false).
"""
function IsItDispersalTime(; setdispersedflag::Union{Bool,Nothing}=nothing)
    is_time = false

    # Default for the optional flag is false.
    setflag = setdispersedflag === nothing ? false : setdispersedflag

    if dispersal_flag[]
        is_time = true
        dispersal_flag[] = false
    else
        if GetCadenceDate() != dispersal_date[]
            is_time = true
            if setflag
                dispersal_flag[] = true
                dispersal_date[] = GetCadenceDate()
            end
        end
    end

    return is_time
end

# =====================================================================================
# GetCadenceDate — current date relevant to the configured dispersal cadence.
# =====================================================================================

"""
    GetCadenceDate() -> Int

Return the current date value matching the configured seed-dispersal cadence
(`hlm_seeddisp_cadence[]`): day-of-month (daily), month-of-year (monthly), or
year (yearly). Mirrors Fortran `GetCadenceDate`.
"""
function GetCadenceDate()
    cadence = hlm_seeddisp_cadence[]
    if cadence == fates_dispersal_cadence_daily
        return hlm_current_day[]
    elseif cadence == fates_dispersal_cadence_monthly
        return hlm_current_month[]
    elseif cadence == fates_dispersal_cadence_yearly
        return hlm_current_year[]
    else
        println(fates_log(),
                " ERROR: An undefined dispersal cadence was specified: ", cadence)
        fates_endrun("FatesDispersalMod: undefined dispersal cadence")
    end
end
