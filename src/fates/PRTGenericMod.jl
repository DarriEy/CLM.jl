# PRTGenericMod.jl
# Julia port of FATES src/fates/parteh/PRTGenericMod.F90
#
# Plant Allocation and Reactive Transport (PART) + Extensible Hypotheses (EH)
# = PARTEH. This is the Non-Specific (Generic) layer: the base classes for the
# state variables (`prt_vartype`/`prt_vartypes`) and the global descriptor class
# (`prt_global_type`), plus the element/organ indexing system and the generic
# getters/setters/initializers that operate on ANY allocation hypothesis.
#
# General idea: PARTEH treats its state variables as objects. Each object maps to
#   1) an organ, 2) a spatial position associated with that organ, 3) a chemical
#   element (carbon isotope or nutrient), aka chemical species.
#
# THIS MODULE SHOULD NOT HAVE TO BE MODIFIED TO ACCOMODATE NEW HYPOTHESES.
# Concrete hypothesis strategies (CarbonMod / CNPMod) extend it in a later batch.
#
# Translation notes:
# - `fates_r8` -> Float64; `fates_int` -> Int.
# - The element/organ index `const`s -> Julia `const`.
# - The abstract `prt_vartypes` + deferred type-bound procedures -> the abstract
#   type `AbstractPRTVartypes` + generic functions whose base methods `error()`
#   (the Fortran "must be extended by a child class" stubs).
# - The generic `prt_vartype` container -> `@kwdef mutable struct` with SoA
#   Vector fields (one position dimension per variable).
# - Boundary-condition pointers (Fortran `pointer`) -> `Ref` wrappers that may be
#   `nothing` when flushed.
# - The module-level `prt_global` singleton -> a `Ref{Union{Nothing,prt_global_type}}`.
#
# Deps: FatesConstantsMod (fates_r8, fates_int, nearzero, calloc_abs_error,
# years_per_day, days_per_sec), FatesGlobals (fates_endrun, fates_log),
# PRTParametersMod (prt_params).

# ---------------------------------------------------------------------------
# String length parameters
# ---------------------------------------------------------------------------
const maxlen_varname   = 128
const maxlen_varsymbol = 32
const maxlen_varunits  = 32
const len_baseunit     = 6

# We use this parameter as the value for which we set un-initialized values
const un_initialized = -9.9e32

# We use this parameter as the value for which we check un-initialized values
const check_initialized = -8.8e32

# ---------------------------------------------------------------------------
# IMPORTANT! All elements in all organs are expressed in KILOGRAMS, all rates
# of change in kilograms/day. This assumption cannot be broken!
# ---------------------------------------------------------------------------
const mass_unit      = "kg"
const mass_rate_unit = "kg/day"

# ---------------------------------------------------------------------------
# Allocation Hypothesis Types (each gets its own module)
# ---------------------------------------------------------------------------
const prt_carbon_allom_hyp   = 1
const prt_cnp_flex_allom_hyp = 2

# ---------------------------------------------------------------------------
# Organ types — public indices mapping the organs in each hypothesis to organs
# acknowledged in the calling model.
# ---------------------------------------------------------------------------
const num_organ_types = 6
const all_organs   = 0   # index for all organs
const leaf_organ   = 1   # index for leaf organs
const fnrt_organ   = 2   # index for fine-root organs
const sapw_organ   = 3   # index for sapwood organs
const store_organ  = 4   # index for storage organs
const repro_organ  = 5   # index for reproductive organs
const struct_organ = 6   # index for structure (dead) organs

# ---------------------------------------------------------------------------
# Element types — public indices mapping the elements (chem species) in each
# hypothesis to the elements acknowledged in the calling model.
# ---------------------------------------------------------------------------
const num_element_types = 6   # Total number of unique elements currently
                              # recognized by PARTEH (max index in list below)

const carbon12_element   = 1
const carbon13_element   = 2
const carbon14_element   = 3
const nitrogen_element   = 4
const phosphorus_element = 5
const potassium_element  = 6

# We may query lists of elements/organs; the carbon elements are the biggest
# list right now. Set a global max for scratch arrays.
const max_spec_per_group = 3

# List of all carbon elements; the special index "all_carbon_elements" implies
# this list of carbon organs.
const carbon_elements_list = [carbon12_element, carbon13_element, carbon14_element]

# Maximum number of leaf age pools allowed on each plant (scratch space).
const max_nleafage = 4

# Minimum allowable L2FR — keeps understory plants from shrinking their roots so
# far they disappear and cause numerical issues.
const l2fr_min = 0.01

# ---------------------------------------------------------------------------
# prt_vartype — holds the state (carbon, nutrients, etc.) for ONE pool of ONE
# plant: a specific organ x element combination. These are vectors (not scalars)
# because a variable may have more than one discrete spatial position (e.g.
# multiple leaf age bins). Nested inside `prt_vartypes`.
#
# Over the control period (probably 1 day), the change in current state (`val`)
# relative to the value at period start (`val0`) equals the time-integrated flux
# terms (net_alloc, turnover, burned, damaged).
# All in [kg].
# ---------------------------------------------------------------------------
Base.@kwdef mutable struct prt_vartype
    val::Vector{Float64}       = Float64[]   # Instantaneous state variable        [kg]
    val0::Vector{Float64}      = Float64[]   # State at the beginning of control   [kg]
    net_alloc::Vector{Float64} = Float64[]   # Net change due to alloc/transport   [kg]
    turnover::Vector{Float64}  = Float64[]   # Losses due to turnover (-> litter)  [kg]
    burned::Vector{Float64}    = Float64[]   # Losses due to burn                  [kg]
    damaged::Vector{Float64}   = Float64[]   # Losses due to damage                [kg]
end

# ---------------------------------------------------------------------------
# prt_bctype — input boundary conditions. Allocated as an array for each plant.
# In Fortran these are pointers to a real and an integer in the calling model;
# here we hold `Ref`s that are `nothing` when flushed/unset.
# ---------------------------------------------------------------------------
Base.@kwdef mutable struct prt_bctype
    rval::Union{Nothing,Base.RefValue{Float64}} = nothing
    ival::Union{Nothing,Base.RefValue{Int}}     = nothing
end

# ---------------------------------------------------------------------------
# prt_vartypes — the object directly attached to each plant (the parent object).
# It contains the state variable objects (`variables`) and the boundary-condition
# arrays (`bc_inout`, `bc_in`, `bc_out`).
#
# Fortran defines this as an extendable abstract base class with deferred,
# hypothesis-specific procedures (DailyPRT, FastPRT, DamageRecovery,
# GetNutrientTarget) plus non-overridable generic procedures. In Julia, the
# abstract base type is `AbstractPRTVartypes`; the generic container below is the
# concrete carrier of the state. Concrete hypothesis subtypes (ported later) will
# subtype `AbstractPRTVartypes` and specialize the deferred generic functions.
# ---------------------------------------------------------------------------
abstract type AbstractPRTVartypes end

Base.@kwdef mutable struct prt_vartypes <: AbstractPRTVartypes
    variables::Vector{prt_vartype} = prt_vartype[]   # state variables and fluxes
    bc_inout::Vector{prt_bctype}   = prt_bctype[]    # boundaries that may be changed
    bc_in::Vector{prt_bctype}      = prt_bctype[]    # protected
    bc_out::Vector{prt_bctype}     = prt_bctype[]    # overwritten
    ode_opt_step::Float64          = 0.0
end

# ---------------------------------------------------------------------------
# state_descriptor_type — packs the names/symbols associated with each variable
# for any given hypothesis.
# ---------------------------------------------------------------------------
Base.@kwdef mutable struct state_descriptor_type
    longname::String  = ""   # verbose name
    symbol::String    = ""   # short symbol
    organ_id::Int     = 0    # global id for organ
    element_id::Int   = 0    # global id for element
    num_pos::Int      = 0    # number of discrete spatial positions
end

# ---------------------------------------------------------------------------
# organ_map_type — helps loop through all variables associated with a specific
# organ. The number of unique variables per organ is capped at num_element_types.
# `var_id` is indexed 1:num_element_types.
# ---------------------------------------------------------------------------
Base.@kwdef mutable struct organ_map_type
    var_id::Vector{Int} = zeros(Int, num_element_types)
    num_vars::Int       = 0
end

# ---------------------------------------------------------------------------
# prt_global_type — describes the mapping for each specific hypothesis. Globally
# true (not per-plant). Instanced once per memory space, then used read-only.
#
#   sp_organ_map : variable id for each (organ, element) combination.
#                  Fortran dims 0:num_organ_types x 0:num_element_types; here a
#                  (num_organ_types+1) x (num_element_types+1) Matrix accessed
#                  through a +1 offset (index 0 -> slot 1).
#   organ_map    : list of variable ids associated with each organ (1:num_organ_types)
# ---------------------------------------------------------------------------
Base.@kwdef mutable struct prt_global_type
    hyp_name::String = ""               # hypothesis name (index 0 reserved for "all")
    hyp_id::Int      = 0                 # hypothesis index (internal dispatch)

    # variable id for each (organ, element); 0-based in Fortran -> +1 offset here.
    sp_organ_map::Matrix{Int} = zeros(Int, num_organ_types + 1, num_element_types + 1)

    # verbose variable descriptors (length num_vars)
    state_descriptor::Vector{state_descriptor_type} = state_descriptor_type[]

    # variable ids associated with each organ (1:num_organ_types)
    organ_map::Vector{organ_map_type} =
        [organ_map_type() for _ in 1:num_organ_types]

    num_bc_in::Int    = 0   # number of input boundary conditions
    num_bc_out::Int   = 0   # number of output boundary conditions
    num_bc_inout::Int = 0   # number of combo input-output boundary conditions
    num_vars::Int     = 0   # number of variables set by the hypothesis
end

# ---------------------------------------------------------------------------
# Offset accessors for the 0-based sp_organ_map (Fortran indices 0:N).
# ---------------------------------------------------------------------------
@inline sp_organ_map_get(g::prt_global_type, organ_id::Integer, element_id::Integer) =
    g.sp_organ_map[organ_id + 1, element_id + 1]
@inline function sp_organ_map_set!(g::prt_global_type, organ_id::Integer,
                                   element_id::Integer, var_id::Integer)
    g.sp_organ_map[organ_id + 1, element_id + 1] = var_id
    return var_id
end

# ---------------------------------------------------------------------------
# Module-level globals (Fortran module variables -> Refs so they are mutable).
#   num_elements : number of elements in this simulation (C, N, P, K, etc.)
#   element_list : global element identifiers (e.g. carbon12_element, ...)
#   element_pos  : reverse lookup (element global index -> position in element_list)
#   prt_global   : the singleton descriptor object.
# ---------------------------------------------------------------------------
const num_elements = Ref{Int}(0)
const element_list = Int[]
const element_pos  = zeros(Int, num_element_types)
const prt_global   = Ref{Union{Nothing,prt_global_type}}(nothing)

# Convenience accessor that errors if the global has not been initialized.
@inline function get_prt_global()
    g = prt_global[]
    g === nothing && fates_endrun("PRTGenericMod: prt_global has not been initialized")
    return g::prt_global_type
end

# =====================================================================================
# prt_global_type methods
# =====================================================================================

"""
    ZeroGlobal!(this::prt_global_type)

Zero out the map between variable indexes and the elements/organs they are
associated with. Also sets the counts of variables and boundary conditions to a
nonsense value (-9) that will trigger a fail if not specified later.
"""
function ZeroGlobal!(this::prt_global_type)
    # First zero out the maps
    for io in 1:num_organ_types
        for is in 1:num_element_types
            sp_organ_map_set!(this, io, is, 0)
            this.organ_map[io].var_id[is] = 0
        end
        this.organ_map[io].num_vars = 0
    end

    # Set the number of boundary conditions as a bogus value
    this.num_bc_in    = -9
    this.num_bc_out   = -9
    this.num_bc_inout = -9

    # Set the number of variables to a bogus value (must be overwritten by caller)
    this.num_vars = -9

    return nothing
end

"""
    RegisterVarInGlobal!(this, var_id, long_name, symbol, organ_id, element_id, num_pos)

Called once for each variable defined in a specific hypothesis (e.g. six times in
the carbon-only hypothesis), providing names, symbols, and the associated
organ/element/position count for each pool. Populates the descriptor and the
mapping tables.
"""
function RegisterVarInGlobal!(this::prt_global_type, var_id::Integer,
                              long_name::AbstractString, symbol::AbstractString,
                              organ_id::Integer, element_id::Integer, num_pos::Integer)
    sd = this.state_descriptor[var_id]
    sd.longname   = String(long_name)
    sd.symbol     = String(symbol)
    sd.organ_id   = organ_id
    sd.element_id = element_id
    sd.num_pos    = num_pos

    # Set the mapping table for the external model
    sp_organ_map_set!(this, organ_id, element_id, var_id)

    # Set the map that locates all relevant pools associated with an organ
    om = this.organ_map[organ_id]
    om.num_vars += 1
    om.var_id[om.num_vars] = var_id

    return nothing
end

# =====================================================================================
# prt_vartypes generic methods
# =====================================================================================

"""
    InitPRTVartype!(this::AbstractPRTVartypes)

First call whenever a `prt_vartypes` object is instantiated (typically when a new
plant/cohort is created). Allocates memory, initializes states to a nan-like
starter value, and flushes all boundary conditions.
"""
function InitPRTVartype!(this::AbstractPRTVartypes)
    InitAllocate!(this)                  # Allocate memory spaces
    InitializeInitialConditions!(this)   # Set states to a nan-like starter value
    FlushBCs!(this)                      # Set all BC pointers to null
    return nothing
end

"""
    InitAllocate!(this::AbstractPRTVartypes)

Called every time a plant/cohort is newly recruited; allocates space for the
boundary-condition arrays and the per-variable SoA state vectors, sized by
`prt_global`.
"""
function InitAllocate!(this::AbstractPRTVartypes)
    g = get_prt_global()

    # Allocate the boundary condition arrays
    if g.num_bc_in > 0
        this.bc_in = [prt_bctype() for _ in 1:g.num_bc_in]
    end
    if g.num_bc_inout > 0
        this.bc_inout = [prt_bctype() for _ in 1:g.num_bc_inout]
    end
    if g.num_bc_out > 0
        this.bc_out = [prt_bctype() for _ in 1:g.num_bc_out]
    end

    # Allocate the state variables
    this.variables = [prt_vartype() for _ in 1:g.num_vars]

    for i_var in 1:g.num_vars
        num_pos = g.state_descriptor[i_var].num_pos
        v = this.variables[i_var]
        v.val       = Vector{Float64}(undef, num_pos)
        v.val0      = Vector{Float64}(undef, num_pos)
        v.turnover  = Vector{Float64}(undef, num_pos)
        v.net_alloc = Vector{Float64}(undef, num_pos)
        v.burned    = Vector{Float64}(undef, num_pos)
        v.damaged   = Vector{Float64}(undef, num_pos)
    end

    return nothing
end

"""
    InitializeInitialConditions!(this::AbstractPRTVartypes)

Set all PARTEH variables to a nonsense value (`un_initialized`) so that a fail is
triggered if a value is not initialized correctly. Also resets `ode_opt_step`.
"""
function InitializeInitialConditions!(this::AbstractPRTVartypes)
    g = get_prt_global()
    for i_var in 1:g.num_vars
        v = this.variables[i_var]
        fill!(v.val,       un_initialized)
        fill!(v.val0,      un_initialized)
        fill!(v.turnover,  un_initialized)
        fill!(v.burned,    un_initialized)
        fill!(v.damaged,   un_initialized)
        fill!(v.net_alloc, un_initialized)
    end

    # Initialize the optimum step size as very large.
    this.ode_opt_step = 1e6

    return nothing
end

"""
    CheckInitialConditions(this::AbstractPRTVartypes)

Make sure every variable defined in the hypothesis has been given an initial
value (i.e. is no longer below `check_initialized`). Errors otherwise.
"""
function CheckInitialConditions(this::AbstractPRTVartypes)
    g = get_prt_global()
    for i_var in 1:g.num_vars
        n_cor_ids = length(this.variables[i_var].val)
        for i_cor in 1:n_cor_ids
            if this.variables[i_var].val[i_cor] < check_initialized
                i_organ   = g.state_descriptor[i_var].organ_id
                i_element = g.state_descriptor[i_var].element_id
                fates_endrun(
                    "Not all initial conditions for state variables in PRT " *
                    "hypothesis $(g.hyp_name) were written out. " *
                    "i_var=$(i_var) i_cor=$(i_cor) organ_id=$(i_organ) " *
                    "element_id=$(i_element)")
            end
        end
    end
    return nothing
end

"""
    FlushBCs!(this::AbstractPRTVartypes)

Boundary conditions are pointers to the calling model's reals/integers. To flush
them, set all pointers to null (`nothing`).
"""
function FlushBCs!(this::AbstractPRTVartypes)
    for i in eachindex(this.bc_in)
        this.bc_in[i].rval = nothing
        this.bc_in[i].ival = nothing
    end
    for i in eachindex(this.bc_out)
        this.bc_out[i].rval = nothing
        this.bc_out[i].ival = nothing
    end
    for i in eachindex(this.bc_inout)
        this.bc_inout[i].rval = nothing
        this.bc_inout[i].ival = nothing
    end
    return nothing
end

"""
    RegisterBCIn!(this, bc_id; bc_rval=nothing, bc_ival=nothing)

Register one "input only" boundary condition. Called once per input-only BC of
each hypothesis, after `InitPRTVartype!`. `bc_rval`/`bc_ival` are `Ref`s into the
calling model's storage.
"""
function RegisterBCIn!(this::AbstractPRTVartypes, bc_id::Integer;
                       bc_rval::Union{Nothing,Base.RefValue{Float64}}=nothing,
                       bc_ival::Union{Nothing,Base.RefValue{Int}}=nothing)
    bc_ival === nothing || (this.bc_in[bc_id].ival = bc_ival)
    bc_rval === nothing || (this.bc_in[bc_id].rval = bc_rval)
    return nothing
end

"""
    RegisterBCOut!(this, bc_id; bc_rval=nothing, bc_ival=nothing)

Register one "output only" boundary condition (cf. [`RegisterBCIn!`](@ref)).
"""
function RegisterBCOut!(this::AbstractPRTVartypes, bc_id::Integer;
                        bc_rval::Union{Nothing,Base.RefValue{Float64}}=nothing,
                        bc_ival::Union{Nothing,Base.RefValue{Int}}=nothing)
    bc_ival === nothing || (this.bc_out[bc_id].ival = bc_ival)
    bc_rval === nothing || (this.bc_out[bc_id].rval = bc_rval)
    return nothing
end

"""
    RegisterBCInout!(this, bc_id; bc_rval=nothing, bc_ival=nothing)

Register one "input-output" boundary condition: passed into PARTEH, expected to
be updated (or not), and passed back to the host.
"""
function RegisterBCInout!(this::AbstractPRTVartypes, bc_id::Integer;
                          bc_rval::Union{Nothing,Base.RefValue{Float64}}=nothing,
                          bc_ival::Union{Nothing,Base.RefValue{Int}}=nothing)
    bc_ival === nothing || (this.bc_inout[bc_id].ival = bc_ival)
    bc_rval === nothing || (this.bc_inout[bc_id].rval = bc_rval)
    return nothing
end

"""
    CopyPRTVartypes!(this, donor_prt_obj)

Copy all state information from a donor PRT object into `this` (assumed already
initialized via [`InitAllocate!`](@ref)). `val0` is copied too (it is ephemeral
but cheap), matching the Fortran which copies every flux/state field.
"""
function CopyPRTVartypes!(this::AbstractPRTVartypes, donor_prt_obj::AbstractPRTVartypes)
    g = get_prt_global()
    for i_var in 1:g.num_vars
        d = donor_prt_obj.variables[i_var]
        t = this.variables[i_var]
        t.val       .= d.val
        t.val0      .= d.val0
        t.net_alloc .= d.net_alloc
        t.turnover  .= d.turnover
        t.burned    .= d.burned
        t.damaged   .= d.damaged
    end
    this.ode_opt_step = donor_prt_obj.ode_opt_step
    return nothing
end

"""
    WeightedFusePRTVartypes!(this, donor_prt_obj, recipient_fuse_weight)

Fuse two PRT objects based on a fusion weighting assigned to the recipient
(`this`). Each position of each variable becomes the weighted average of the
recipient and donor values.
"""
function WeightedFusePRTVartypes!(this::AbstractPRTVartypes,
                                  donor_prt_obj::AbstractPRTVartypes,
                                  recipient_fuse_weight::Real)
    g = get_prt_global()
    w  = recipient_fuse_weight
    dw = 1.0 - recipient_fuse_weight
    for i_var in 1:g.num_vars
        t = this.variables[i_var]
        d = donor_prt_obj.variables[i_var]
        for pos_id in 1:g.state_descriptor[i_var].num_pos
            t.val[pos_id]       = w * t.val[pos_id]       + dw * d.val[pos_id]
            t.val0[pos_id]      = w * t.val0[pos_id]      + dw * d.val0[pos_id]
            t.net_alloc[pos_id] = w * t.net_alloc[pos_id] + dw * d.net_alloc[pos_id]
            t.turnover[pos_id]  = w * t.turnover[pos_id]  + dw * d.turnover[pos_id]
            t.burned[pos_id]    = w * t.burned[pos_id]    + dw * d.burned[pos_id]
            t.damaged[pos_id]   = w * t.damaged[pos_id]   + dw * d.damaged[pos_id]
        end
    end
    this.ode_opt_step = w * this.ode_opt_step + dw * donor_prt_obj.ode_opt_step
    return nothing
end

"""
    DeallocatePRTVartypes!(this::AbstractPRTVartypes)

Release the per-plant PRT storage. Julia is garbage-collected, so this simply
drops the arrays (sets them to empty), mirroring the Fortran deallocate.
"""
function DeallocatePRTVartypes!(this::AbstractPRTVartypes)
    this.variables = prt_vartype[]
    this.bc_in     = prt_bctype[]
    this.bc_out    = prt_bctype[]
    this.bc_inout  = prt_bctype[]
    return nothing
end

"""
    ZeroRates!(this::AbstractPRTVartypes)

Zero all rates of change and set the initial value (`val0`) to the current state
(`val`). This enables mass-conservation checks, where over the control period
`val - val0 == net_alloc - turnover - burned - damaged`. Called each day in FATES.
"""
function ZeroRates!(this::AbstractPRTVartypes)
    g = get_prt_global()
    for i_var in 1:g.num_vars
        v = this.variables[i_var]
        v.val0 .= v.val
        fill!(v.net_alloc, 0.0)
        fill!(v.turnover,  0.0)
        fill!(v.burned,    0.0)
        fill!(v.damaged,   0.0)
    end
    return nothing
end

"""
    CheckMassConservation(this, ipft, position_id)

At any time, the sum of fluxes should equal the difference between `val` and
`val0`. Loop over all variables/positions and verify
`|(val-val0) - (net_alloc - turnover - burned - damaged)| <= calloc_abs_error`,
erroring otherwise. `ipft`/`position_id` identify the failing call site.
"""
function CheckMassConservation(this::AbstractPRTVartypes, ipft::Integer, position_id::Integer)
    g = get_prt_global()
    for i_var in 1:g.num_vars
        v = this.variables[i_var]
        for i_pos in 1:g.state_descriptor[i_var].num_pos
            err = abs((v.val[i_pos] - v.val0[i_pos]) -
                      (v.net_alloc[i_pos] - v.turnover[i_pos] -
                       v.burned[i_pos] - v.damaged[i_pos]))

            rel_err = v.val[i_pos] > nearzero ? err / v.val[i_pos] : 0.0  # noqa: kept for parity

            if abs(err) > calloc_abs_error
                sd = g.state_descriptor[i_var]
                fates_endrun(
                    "PARTEH mass conservation check failed. " *
                    "Change in mass over control period should equal the " *
                    "integrated fluxes. pft id=$(ipft) position id=$(position_id) " *
                    "organ id=$(sd.organ_id) element_id=$(sd.element_id) " *
                    "i_pos=$(i_pos) symbol=$(sd.symbol) longname=$(sd.longname) " *
                    "err=$(err) max error=$(calloc_abs_error) terms=" *
                    "[$(v.val[i_pos]), $(v.val0[i_pos]), $(v.net_alloc[i_pos]), " *
                    "$(v.turnover[i_pos]), $(v.burned[i_pos]), $(v.damaged[i_pos])]")
            end
        end
    end
    return nothing
end

"""
    GetState(this, organ_id, element_id [, position_id]) -> Float64

Return the current mass [kg] for an (organ, element) combination. If
`position_id` is given, return just that position; otherwise sum over all
positions.
"""
function GetState(this::AbstractPRTVartypes, organ_id::Integer, element_id::Integer,
                  position_id::Union{Nothing,Integer}=nothing)
    g = get_prt_global()
    i_var = sp_organ_map_get(g, organ_id, element_id)
    if position_id !== nothing
        return this.variables[i_var].val[position_id]
    else
        state_val = 0.0
        for i_pos in 1:g.state_descriptor[i_var].num_pos
            state_val += this.variables[i_var].val[i_pos]
        end
        return state_val
    end
end

"""
    GetTurnover(this, organ_id, element_id [, position_id]) -> Float64

Return the turnover mass [kg] so far during the control period. Query only — does
NOT specify turnover.
"""
function GetTurnover(this::AbstractPRTVartypes, organ_id::Integer, element_id::Integer,
                     position_id::Union{Nothing,Integer}=nothing)
    g = get_prt_global()
    i_var = sp_organ_map_get(g, organ_id, element_id)
    if position_id !== nothing
        return this.variables[i_var].turnover[position_id]
    else
        turnover_val = 0.0
        for i_pos in 1:g.state_descriptor[i_var].num_pos
            turnover_val += this.variables[i_var].turnover[i_pos]
        end
        return turnover_val
    end
end

"""
    GetBurned(this, organ_id, element_id [, position_id]) -> Float64

Return the burned mass [kg] so far during the control period. Query only.
"""
function GetBurned(this::AbstractPRTVartypes, organ_id::Integer, element_id::Integer,
                   position_id::Union{Nothing,Integer}=nothing)
    g = get_prt_global()
    i_var = sp_organ_map_get(g, organ_id, element_id)
    if position_id !== nothing
        return this.variables[i_var].burned[position_id]
    else
        burned_val = 0.0
        for i_pos in 1:g.state_descriptor[i_var].num_pos
            burned_val += this.variables[i_var].burned[i_pos]
        end
        return burned_val
    end
end

"""
    GetNetAlloc(this, organ_id, element_id [, position_id]) -> Float64

Return the net change due to allocation/reactions/transport [kg] in that pool.
Query only.
"""
function GetNetAlloc(this::AbstractPRTVartypes, organ_id::Integer, element_id::Integer,
                     position_id::Union{Nothing,Integer}=nothing)
    g = get_prt_global()
    i_var = sp_organ_map_get(g, organ_id, element_id)
    if position_id !== nothing
        return this.variables[i_var].net_alloc[position_id]
    else
        val_netalloc = 0.0
        for i_pos in 1:g.state_descriptor[i_var].num_pos
            val_netalloc += this.variables[i_var].net_alloc[i_pos]
        end
        return val_netalloc
    end
end

"""
    SetState!(prt, organ_id, element_id, state_val [, position_id])

Initialize the state value of a plant's pool for an (organ, element) couplet
(and optional position; defaults to 1). Errors if the position exceeds the
allocated space or if the organ/element combination does not exist.
"""
function SetState!(prt::AbstractPRTVartypes, organ_id::Integer, element_id::Integer,
                   state_val::Real, position_id::Union{Nothing,Integer}=nothing)
    g = get_prt_global()
    i_pos = position_id === nothing ? 1 : position_id

    i_var = sp_organ_map_get(g, organ_id, element_id)

    if i_var > 0 && i_pos > g.state_descriptor[i_var].num_pos
        fates_endrun(
            "A position index was specified that is greater than the allocated " *
            "position space. i_pos=$(i_pos) num_pos=$(g.state_descriptor[i_var].num_pos)")
    end

    if i_var > 0
        prt.variables[i_var].val[i_pos] = state_val
    else
        fates_endrun(
            "A mass was sent to PARTEH to over-write a pool with a specie x organ " *
            "combination that does not exist. organ_id=$(organ_id) element_id=$(element_id)")
    end

    return nothing
end

# =====================================================================================
# Deferred (hypothesis-specific) generic functions — base methods error, to be
# specialized by concrete subtypes of AbstractPRTVartypes (ported in a later batch).
# =====================================================================================

"""
    DailyPRT!(this::AbstractPRTVartypes, phase)

Daily PARTEH allocation. Base implementation errors — must be extended by a
concrete hypothesis (child class).
"""
function DailyPRT!(::AbstractPRTVartypes, phase::Integer)
    fates_endrun("Daily PRT Allocation must be extended")
end

"""
    FastPRT!(this::AbstractPRTVartypes)

Fast (sub-daily) reactive transport. Base implementation errors — must be
extended by a concrete hypothesis.
"""
function FastPRT!(::AbstractPRTVartypes)
    fates_endrun("FastReactiveTransport must be extended by a child class")
end

"""
    DamageRecovery!(this::AbstractPRTVartypes)

Damage recovery. Base implementation errors — must be extended.
"""
function DamageRecovery!(::AbstractPRTVartypes)
    fates_endrun("DamageRecovery must be extended by a child class")
end

"""
    GetNutrientTarget(this, element_id, organ_id [, stoich_mode]) -> Float64

Target nutrient mass for an organ. Base implementation errors — must be extended.
"""
function GetNutrientTarget(::AbstractPRTVartypes, element_id::Integer, organ_id::Integer,
                           stoich_mode::Union{Nothing,Integer}=nothing)
    fates_endrun("GetNutrientTargetBase must be extended by a child class.")
end

"""
    GetCoordVal(this, organ_id, element_id) -> Float64

Support code for variables with multiple discrete spatial positions. Base
implementation errors — must be extended by a child class.
"""
function GetCoordVal(::AbstractPRTVartypes, organ_id::Integer, element_id::Integer)
    fates_endrun("Init must be extended by a child class.")
end

# =====================================================================================
# AgeLeaves — generic, but over-ridable.
# =====================================================================================

"""
    AgeLeaves!(this, ipft, icanlayer, period_sec)

If there is more than one leaf age class, allow some leaf biomass to transition
from younger to older bins (called once per day). No mass leaves the plant; mass
is moved from a young bin to the next older bin (none out of the oldest). Updates
the daily `net_alloc` diagnostic. `period_sec` is the call period [s] (daily=86400).
"""
function AgeLeaves!(this::AbstractPRTVartypes, ipft::Integer, icanlayer::Integer,
                    period_sec::Real)
    g = get_prt_global()
    leaf_m0 = zeros(Float64, max_nleafage)

    for el in 1:num_elements[]
        element_id = element_list[el]

        # Global position of leaf variable
        i_var = sp_organ_map_get(g, leaf_organ, element_id)

        # Number of leaf age class bins
        nleafage = g.state_descriptor[i_var].num_pos

        leaf_m = this.variables[i_var].val

        @views leaf_m0[1:nleafage] .= leaf_m[1:nleafage]

        if nleafage > 1
            for i_age in 1:(nleafage - 1)
                if prt_params.leaf_long[ipft, i_age] > nearzero
                    leaf_long = icanlayer == 1 ? prt_params.leaf_long[ipft, i_age] :
                                                 prt_params.leaf_long_ustory[ipft, i_age]

                    # Units: [-] = [sec] * [day/sec] * [years/day] * [1/years]
                    leaf_age_flux_frac = period_sec * days_per_sec * years_per_day / leaf_long

                    leaf_m[i_age]     -= leaf_m0[i_age] * leaf_age_flux_frac
                    leaf_m[i_age + 1] += leaf_m0[i_age] * leaf_age_flux_frac
                end
            end
        end

        # Update the diagnostic on daily rate of change
        for i_age in 1:nleafage
            this.variables[i_var].net_alloc[i_age] += (leaf_m[i_age] - leaf_m0[i_age])
        end
    end

    return nothing
end

# =====================================================================================
# StorageNutrientTarget — module function (not type-bound).
# =====================================================================================

# Choice of how nutrient storage target is proportioned. Each choice makes the
# nutrient storage proportional to the "in-tissue" total nitrogen content of one
# or more sets of organs.
const lfs_store_prop    = 1   # leaf-sapwood proportional storage
const lfss_store_prop   = 2   # leaf-fnrt-sapw-struct proportional storage
const fnrt_store_prop   = 3   # fineroot proportional storage
const cstore_store_prop = 4   # proportional to carbon storage times mean CN
const lf_store_prop     = 5   # leaf proportional storage
const store_prop        = lf_store_prop

"""
    StorageNutrientTarget(pft, element_id, leaf_target, fnrt_target, sapw_target, struct_target) -> Float64

Target storage nutrient mass [kg] for a plant, proportioned according to the
`store_prop` strategy (currently `lf_store_prop`, leaf-proportional). Errors for
carbon, and for the unimplemented `cstore_store_prop` method.
"""
function StorageNutrientTarget(pft::Integer, element_id::Integer,
                               leaf_target::Real, fnrt_target::Real,
                               sapw_target::Real, struct_target::Real)
    store_target = 0.0

    if element_id == carbon12_element
        fates_endrun("Cannot call StorageNutrientTarget() for carbon")

    elseif element_id == nitrogen_element
        if store_prop == lfs_store_prop
            store_target = prt_params.nitr_store_ratio[pft] * (leaf_target + sapw_target)
        elseif store_prop == lf_store_prop
            store_target = prt_params.nitr_store_ratio[pft] * leaf_target
        elseif store_prop == lfss_store_prop
            store_target = prt_params.nitr_store_ratio[pft] *
                           (leaf_target + fnrt_target + sapw_target + struct_target)
        elseif store_prop == fnrt_store_prop
            store_target = prt_params.nitr_store_ratio[pft] * fnrt_target
        elseif store_prop == cstore_store_prop
            fates_endrun("cstore_store_prop method of calculating target nutrient stores not available")
        end

    elseif element_id == phosphorus_element
        if store_prop == lfs_store_prop
            store_target = prt_params.phos_store_ratio[pft] *
                           (leaf_target + fnrt_target + sapw_target)
        elseif store_prop == lf_store_prop
            store_target = prt_params.phos_store_ratio[pft] * leaf_target
        elseif store_prop == lfss_store_prop
            store_target = prt_params.phos_store_ratio[pft] *
                           (leaf_target + fnrt_target + sapw_target + struct_target)
        elseif store_prop == fnrt_store_prop
            store_target = prt_params.phos_store_ratio[pft] * fnrt_target
        end
    end

    return store_target
end
