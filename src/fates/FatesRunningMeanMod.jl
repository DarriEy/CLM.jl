# FatesRunningMeanMod.jl
# Julia port of FATES src/fates/main/FatesRunningMeanMod.F90
#
# Running-mean accumulator types and update methods. Supports exponential moving
# average (EMA, indefinite) and fixed-window (zero'd each interval) means.
# Fortran type-bound procedures -> bang-functions taking the struct as 1st arg.
# Pointers (`def_type`) -> a direct reference to an rmean_def_type instance.
# Deps: FatesConstantsMod (nearzero), FatesGlobals (fates_log, fates_endrun).

# Averaging-window method flags.
const moving_ema_window = 0   # exponential moving average
const fixed_window = 1

const fates_running_mean_debug = true

"""
    rmean_def_type

Defines a *type* of mean: how often it is updated, how long its memory period
is, and whether it is zero'd (fixed window) or indefinite (EMA). These are
globally defined on the proc.
"""
Base.@kwdef mutable struct rmean_def_type
    mem_period::Float64 = 0.0  # total integration period (s)
    up_period::Float64 = 0.0   # period between updates (s)
    n_mem::Int = 0             # how many updates per integration period
    method::Int = moving_ema_window  # fixed or moving window
end

"""
    rmean_type

Holds the time-varying state for a running mean (instantiated on sites, patches,
cohorts). `def_type` references the `rmean_def_type` defining this mean's nature
(`nothing` until initialized).
"""
Base.@kwdef mutable struct rmean_type
    c_mean::Float64 = NaN   # current mean (full for EMA, partial for fixed window)
    l_mean::Float64 = NaN   # latest reportable mean value
    c_index::Int = 0        # number of values added so far
    def_type::Union{rmean_def_type,Nothing} = nothing
end

"""
    rmean_arr_type

A wrapper holding a reference to an `rmean_type` (for per-PFT etc. arrays).
"""
Base.@kwdef mutable struct rmean_arr_type
    p::Union{rmean_type,Nothing} = nothing
end

# =====================================================================================

"""
    define!(this::rmean_def_type, mem_period, up_period, method)

Configure a running-mean definition. Errors (in debug) if the update period is
not an exact fraction of the memory period.
"""
function define!(this::rmean_def_type, mem_period::Real, up_period::Real, method::Integer)
    if fates_running_mean_debug
        if abs(round(Int, mem_period / up_period) - mem_period / up_period) > nearzero
            @warn "While defining a running mean definition, the update period " *
                  "is not an exact fraction of the memory period " *
                  "(mem_period=$mem_period, up_period=$up_period)"
            fates_endrun("running mean define: up_period not an exact fraction of mem_period")
        end
    end

    this.mem_period = mem_period
    this.up_period = up_period
    this.method = method
    this.n_mem = round(Int, mem_period / up_period)
    return nothing
end

# =====================================================================================

"""
    GetMean(this::rmean_type) -> Float64

Return the latest reportable mean. Errors (in debug) for an EMA window that has
never been given a value.
"""
function GetMean(this::rmean_type)
    if this.def_type.method == moving_ema_window
        if this.c_index == 0 && fates_running_mean_debug
            @warn "attempting to get a running mean from a variable that has not " *
                  "been given a value yet"
            fates_endrun("GetMean on uninitialized EMA")
        end
    end
    return this.l_mean
end

# =====================================================================================

"""
    InitRMean!(this::rmean_type, rmean_def; init_value=nothing, init_offset=nothing)

Initialize a running mean, pointing it at its definition `rmean_def`. For a
fixed window both `init_value` and `init_offset` are required; for an EMA both
are optional.
"""
function InitRMean!(this::rmean_type, rmean_def::rmean_def_type;
                    init_value::Union{Real,Nothing}=nothing,
                    init_offset::Union{Real,Nothing}=nothing)

    # Point to the averaging type.
    this.def_type = rmean_def

    if this.def_type.method == fixed_window

        if fates_running_mean_debug
            if !(init_offset !== nothing && init_value !== nothing)
                @warn "when initializing a temporal mean on a fixed window there " *
                      "must be an initial value and a time offset specified."
                fates_endrun("InitRMean fixed window missing init_value/init_offset")
            end

            # Check the offset is an even increment of the update frequency.
            if abs(Float64(round(Int, init_offset / rmean_def.up_period)) -
                   (init_offset / rmean_def.up_period)) > nearzero
                @warn "when initializing a temporal mean on a fixed window the time " *
                      "offset must be an increment of the update frequency " *
                      "(offset=$init_offset, up freq=$(rmean_def.up_period))"
                fates_endrun("InitRMean fixed window offset not an increment")
            end

            if init_offset < -nearzero
                @warn "offset must be positive: $init_offset"
                fates_endrun("InitRMean fixed window negative offset")
            end
        end

        this.c_index = mod(round(Int, init_offset / rmean_def.up_period), rmean_def.n_mem)
        this.c_mean = Float64(this.c_index) / Float64(rmean_def.n_mem) * init_value
        this.l_mean = init_value

    elseif this.def_type.method == moving_ema_window

        if init_value !== nothing
            this.c_mean = init_value
            this.l_mean = init_value
            if init_offset !== nothing
                this.c_index = min(round(Int, init_offset / rmean_def.up_period), rmean_def.n_mem)
            else
                this.c_index = 1
            end
        else
            this.c_mean = NaN
            this.l_mean = NaN
            this.c_index = 0
        end
    end

    return nothing
end

# =====================================================================================

"""
    CopyFromDonor!(this::rmean_type, donor::rmean_type)

Copy the running-mean state from `donor` into `this`. Errors if `this` has no
`def_type` associated.
"""
function CopyFromDonor!(this::rmean_type, donor::rmean_type)
    if this.def_type === nothing
        @warn "Attempted to copy over running mean info from a donor into a new " *
              "structure but the new structure did not have its def_type associated"
        fates_endrun("CopyFromDonor with unassociated def_type")
    end

    this.c_mean = donor.c_mean
    this.l_mean = donor.l_mean
    this.c_index = donor.c_index
    return nothing
end

# =====================================================================================

"""
    UpdateRMean!(this::rmean_type, new_value)

Add `new_value` to the running mean. For EMA the weight is 1/c_index (capped at
n_mem); for a fixed window values accumulate with equal weight and the window is
reported and zero'd once `n_mem` updates have occurred.
"""
function UpdateRMean!(this::rmean_type, new_value::Real)
    if this.def_type.method == moving_ema_window

        this.c_index = min(this.def_type.n_mem, this.c_index + 1)

        if this.c_index == 1
            this.c_mean = new_value
        else
            wgt = 1.0 / Float64(this.c_index)
            this.c_mean = this.c_mean * (1.0 - wgt) + wgt * new_value
        end

        this.l_mean = this.c_mean

    else
        # Fixed window.
        this.c_index = this.c_index + 1
        wgt = this.def_type.up_period / this.def_type.mem_period
        this.c_mean = this.c_mean + new_value * wgt

        if this.c_index == this.def_type.n_mem
            this.l_mean = this.c_mean
            this.c_mean = 0.0
            this.c_index = 0
        end
    end

    return nothing
end

# =====================================================================================

"""
    FuseRMean!(this::rmean_type, donor::rmean_type, recip_wgt)

Fuse `donor` into `this` with recipient weight `recip_wgt` (0-1). For a fixed
window the two indices must match. Uninitialized donors are left as no-ops.
"""
function FuseRMean!(this::rmean_type, donor::rmean_type, recip_wgt::Real)
    if this.def_type.method == fixed_window
        if this.c_index != donor.c_index
            @warn "trying to fuse two fixed-window averages that are at different " *
                  "points in the window? c_mean=($(this.c_mean),$(donor.c_mean)) " *
                  "l_mean=($(this.l_mean),$(donor.l_mean)) " *
                  "c_index=($(this.c_index),$(donor.c_index))"
            fates_endrun("FuseRMean fixed-window index mismatch")
        end
    end

    # Avoid doing math on uninitialized values.
    if !(donor.c_index == 0)
        if this.c_index == 0
            this.c_mean = donor.c_mean
            this.l_mean = donor.l_mean
            this.c_index = donor.c_index
        else
            this.c_mean = this.c_mean * recip_wgt + donor.c_mean * (1.0 - recip_wgt)
            this.l_mean = this.l_mean * recip_wgt + donor.l_mean * (1.0 - recip_wgt)
            this.c_index = max(this.c_index, donor.c_index)
        end
    end

    return nothing
end
