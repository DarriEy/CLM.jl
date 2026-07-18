# dyn_file_io.jl — Time-stepping framework for transient (dynamic) land-use files.
#
# Julia port of three tightly-coupled Fortran modules from CLM/CTSM's
# src/dyn_subgrid/:
#   - dynTimeInfoMod.F90          -> DynTimeInfo  (maps model year -> file YEAR axis)
#   - dynFileMod.F90              -> DynFile      (NetCDF handle + embedded time_info)
#   - dynVarTimeUninterpMod.F90.in-> DynVarTimeUninterp (no-interp variable reader)
#     (with the metadata/read-data structure from the abstract base dynVarMod.F90.in)
#
# Only the NO-interpolation variable variant is ported (the interpolation variant
# from dynVarTimeInterpMod.F90.in is a later optional enhancement). For an
# uninterpolated variable, the data snap to their new value at the beginning of
# each year and then stay fixed; before the start of the time series the data are
# fixed at the first year's value, and after the end they are fixed at the last
# year's value.
#
# Fortran variable/subroutine names are preserved for traceability. The model has
# no Fortran clm_time_manager coupling here, so `set_current_year!` takes the model
# year explicitly (the caller supplies the year that get_prev_date / get_curr_date
# would have produced, per the stored year_position).
#
# NetCDF I/O uses NCDatasets (already a CLM.jl dependency; `using NCDatasets` lives
# in src/CLM.jl), mirroring src/infrastructure/history_io.jl and restart_io.jl.

# ======================================================================
# Year-position enum (dynTimeInfoMod year_position_type)
# ======================================================================
# The Fortran defines a private year_position_type holding an integer flag, with
# two public parameter instances. We mirror that as a small @enum.
@enum YearPosition begin
    YEAR_POSITION_START_OF_TIMESTEP = 1
    YEAR_POSITION_END_OF_TIMESTEP   = 2
end

# ======================================================================
# DynTimeInfo  (dynTimeInfoMod :: time_info_type)
# ======================================================================
# Stores time information for a single dynamic landuse file, assuming a single
# time sample per year. Maps the current model year to a position in the file's
# YEAR axis via (time_index_lower, time_index_upper), with year-bounding logic for
# years before / on / after the file's range.
Base.@kwdef mutable struct DynTimeInfo
    # Static information about the file:
    nyears::Int = 0                  # number of years in the file
    years::Vector{Int} = Int[]       # all years in this file

    # Other static information:
    year_position::YearPosition = YEAR_POSITION_START_OF_TIMESTEP

    # Information that potentially changes each time step:
    time_index_lower::Int = 1        # lower bound index of the current interval
    time_index_upper::Int = 1        # upper bound index of the current interval
end

"""
    dyn_time_info(years::Vector{Int}, year_position::YearPosition; current_year=nothing)

Initialize a `DynTimeInfo` object (dynTimeInfoMod constructor).

`years` are all years in the file; `year_position` says how to obtain the model
year relative to the current timestep. `time_index_lower`/`time_index_upper` are
set arbitrarily to 1 first (mirroring Fortran) and then set correctly via
`set_current_year!`. The Fortran constructor calls `set_current_year()` with no
argument (which pulls the year from the time manager). Here, if `current_year` is
supplied it is used; otherwise the indices are left at the constructor defaults
(both = 1, i.e. "before the time series") until the caller invokes
`set_current_year!`.
"""
function dyn_time_info(years::Vector{Int}, year_position::YearPosition;
                       current_year::Union{Int,Nothing} = nothing)
    nyears = length(years)
    ti = DynTimeInfo(
        nyears = nyears,
        years = copy(years),
        year_position = year_position,
        # Set arbitrarily; corrected by set_current_year! below.
        time_index_lower = 1,
        time_index_upper = 1,
    )
    if current_year !== nothing
        set_current_year!(ti, current_year)
    end
    return ti
end

# ----------------------------------------------------------------------
# Public method: set_current_year! (dynTimeInfoMod :: set_current_year)
# ----------------------------------------------------------------------
"""
    set_current_year!(ti::DynTimeInfo, year::Int)

Update time information (`time_index_lower` and `time_index_upper`) based on the
current model `year`. Should be called every time step.

In Fortran, `set_current_year` obtains `year` from the time manager
(`get_prev_date` for START_OF_TIMESTEP, `get_curr_date` for END_OF_TIMESTEP) and
then calls `set_info_from_year`. Here the year is passed in directly — the caller
is responsible for choosing prev vs. curr year according to `ti.year_position`.
"""
function set_current_year!(ti::DynTimeInfo, year::Int)
    _set_info_from_year!(ti, year)
    return nothing
end

# ----------------------------------------------------------------------
# Getter routines
# ----------------------------------------------------------------------

# get_time_index_lower / get_time_index_upper
get_time_index_lower(ti::DynTimeInfo) = ti.time_index_lower
get_time_index_upper(ti::DynTimeInfo) = ti.time_index_upper

"""
    get_year(ti::DynTimeInfo, nt::Int)

Get the year associated with time index `nt` (dynTimeInfoMod :: get_year).
"""
function get_year(ti::DynTimeInfo, nt::Int)
    @assert 1 <= nt <= ti.nyears "get_year: nt out of bounds"
    return ti.years[nt]
end

"""
    is_before_time_series(ti::DynTimeInfo)

Returns true if we are currently prior to the bounds of this file
(dynTimeInfoMod :: is_before_time_series). True iff `time_index_upper == 1`.
"""
is_before_time_series(ti::DynTimeInfo) = (ti.time_index_upper == 1)

"""
    is_after_time_series(ti::DynTimeInfo)

Returns true if we are currently after the bounds of this file
(dynTimeInfoMod :: is_after_time_series). If the last year of the file is (e.g.)
2005, then this is TRUE if the current year is 2005. True iff
`time_index_lower == nyears`.
"""
is_after_time_series(ti::DynTimeInfo) = (ti.time_index_lower == ti.nyears)

"""
    is_within_bounds(ti::DynTimeInfo)

Returns true if we are currently within the bounds of this file
(dynTimeInfoMod :: is_within_bounds).
"""
is_within_bounds(ti::DynTimeInfo) =
    (!is_before_time_series(ti)) && (!is_after_time_series(ti))

# The Fortran type does not store the "current year" explicitly — it stores only
# the resulting (lower, upper) indices. We expose `get_current_year` as the year
# at the current lower index, which is the model's notion of the active year on
# the file's axis (clamped to [years[1], years[end]] by the bounding logic).
"""
    get_current_year(ti::DynTimeInfo)

Get the year at the current lower time index (the active year on the file axis,
clamped to the file's range by the year-bounding logic).
"""
function get_current_year(ti::DynTimeInfo)
    return ti.years[ti.time_index_lower]
end

# ----------------------------------------------------------------------
# Private methods
# ----------------------------------------------------------------------

# year_in_current_interval (dynTimeInfoMod :: year_in_current_interval)
function _year_in_current_interval(ti::DynTimeInfo, cur_year::Int)
    if ti.years[ti.time_index_lower] == cur_year &&
       ti.years[ti.time_index_upper] == (cur_year + 1)
        # Normal case: within the time series, in the same interval as before
        return true
    elseif is_before_time_series(ti) && cur_year < ti.years[1]
        # We were and still are before the time series
        return true
    elseif is_after_time_series(ti) && cur_year >= ti.years[ti.nyears]
        # We were and still are after the time series
        return true
    else
        return false
    end
end

# set_info_from_year (dynTimeInfoMod :: set_info_from_year)
#
# Given the current model year, set time_index_lower and time_index_upper.
#   - year < years(1)        -> lower=upper=1            (constant pre-period)
#   - year >= years(nyears)  -> lower=upper=nyears       (constant post-period)
#   - otherwise              -> find n with years(n)==year; lower=n, upper=n+1
function _set_info_from_year!(ti::DynTimeInfo, cur_year::Int)
    nyears = ti.nyears
    years = ti.years

    if _year_in_current_interval(ti, cur_year)
        # DO NOTHING — NT1 AND NT2 ARE ALREADY CORRECT
    else
        if cur_year < years[1]
            # prior to the first interval
            ti.time_index_lower = 1
            ti.time_index_upper = 1
        elseif cur_year >= years[nyears]
            # past the last interval
            ti.time_index_lower = nyears
            ti.time_index_upper = nyears
        else
            # within the time bounds of the file
            found = false
            for n in 1:(nyears - 1)
                if cur_year == years[n]
                    ti.time_index_lower = n
                    ti.time_index_upper = n + 1
                    found = true
                    break
                end
            end
            if !found
                error("set_info_from_year ERROR: model year not found in pftdyn " *
                      "timeseries; model year = $cur_year")
            end
        end
    end

    @assert ti.time_index_upper <= nyears "set_info_from_year: time_index_upper should not be greater than nyears"
    return nothing
end

# ======================================================================
# DynFile  (dynFileMod :: dyn_file_type)
# ======================================================================
# Wraps a NetCDF file (path + open dataset handle) and embeds a DynTimeInfo built
# from the file's YEAR variable (assumed to lie along the 'time' dimension).
Base.@kwdef mutable struct DynFile
    filename::String = ""
    ds::Union{NCDataset,Nothing} = nothing   # open NetCDF dataset handle
    time_info::DynTimeInfo = DynTimeInfo()    # time information for this file
end

"""
    dyn_file_open(filename::String, year_position::YearPosition; current_year=nothing)

Initialize a `DynFile` object (dynFileMod constructor).

Opens the file for reading, reads the `YEAR` variable (assumed to have dimension
`time`), and builds the embedded `DynTimeInfo` from this YEAR variable and (if
given) the current model year.
"""
function dyn_file_open(filename::String, year_position::YearPosition;
                       current_year::Union{Int,Nothing} = nothing)
    isfile(filename) || error("dyn_file_open: file not found: $filename")

    ds = NCDataset(filename, "r")
    # Obtain years: read the YEAR variable along the 'time' dimension.
    haskey(ds, "YEAR") || error("dyn_file_open: 'YEAR' variable not on file $filename")
    years = Int.(Array(ds["YEAR"]))

    time_info = dyn_time_info(years, year_position; current_year = current_year)

    return DynFile(filename = filename, ds = ds, time_info = time_info)
end

"""
    dyn_file_close!(df::DynFile)

Close the underlying NetCDF dataset (no Fortran analogue beyond file_desc_t close,
but needed for clean resource handling in Julia).
"""
function dyn_file_close!(df::DynFile)
    if df.ds !== nothing
        close(df.ds)
        df.ds = nothing
    end
    return nothing
end

# Convenience: forward set_current_year! to the embedded time_info.
set_current_year!(df::DynFile, year::Int) = set_current_year!(df.time_info, year)

# ======================================================================
# DynVarTimeUninterp  (dynVarTimeUninterpMod :: dyn_var_time_uninterp_type)
# ======================================================================
# A single dynamic subgrid variable that is NOT interpolated in time. Holds the
# metadata common to all dyn_var_type objects (from the abstract dynVarMod base)
# plus the uninterpolated data buffer `data_at_tlower`.
#
# Data are stored as a flat 1-d vector (`data_at_tlower`) regardless of the true
# dimensionality, then reshaped to `data_shape` by the get_current_data accessors,
# mirroring the Fortran implementation.
Base.@kwdef mutable struct DynVarTimeUninterp
    # --- metadata (dynVarMod :: set_metadata) ---
    dyn_file::DynFile                       # file containing this variable
    varname::String = ""                    # variable name on file
    dim1name::String = ""                    # dim1name on file (spatial dim)
    conversion_factor::Float64 = 1.0         # data are DIVIDED by this right after reading
    do_check_sums_equal_1::Bool = false      # only relevant for 2-d vars
    data_shape::Vector{Int} = Int[]          # shape of data; first dim is spatial
    allow_nodata::Bool = false               # allow field to be absent on file
    data_on_file::Bool = false               # whether the data was actually on the file

    # --- uninterpolated data (dyn_var_time_uninterp_type) ---
    data_at_tlower::Vector{Float64} = Float64[]  # data at time time_index_lower
    time_index_lower::Int = 0                    # current lower index of the stored data
end

# dyn_var_max_dims (dynVarMod): maximum number of real dimensions allowed.
const DYN_VAR_MAX_DIMS = 2

"""
    dyn_var_time_uninterp(dyn_file, varname, dim1name, conversion_factor,
                          do_check_sums_equal_1, data_shape; allow_nodata=false)

Create a `DynVarTimeUninterp` object (dynVarTimeUninterpMod constructor). This also
reads the first set of data. Assumes `dyn_file` has already been initialized
(its `time_info` reflects the current model year).

`data_shape` gives the shape of the variable (first dim = spatial); its product is
the number of elements. `conversion_factor`: data are divided by it immediately
after reading. `allow_nodata`: if true, a missing field on the file is set to zero
rather than causing an error.
"""
function dyn_var_time_uninterp(dyn_file::DynFile,
                               varname::String,
                               dim1name::String,
                               conversion_factor::Float64,
                               do_check_sums_equal_1::Bool,
                               data_shape::Vector{Int};
                               allow_nodata::Bool = false)
    ndims = length(data_shape)
    # set_metadata error checking (dynVarMod :: set_metadata)
    @assert ndims <= DYN_VAR_MAX_DIMS "set_metadata ERROR: ndims must be <= dyn_var_max_dims"
    if do_check_sums_equal_1
        @assert ndims == 2 "set_metadata ERROR: do_check_sums_equal_1 only valid for ndims==2"
    end

    var = DynVarTimeUninterp(
        dyn_file = dyn_file,
        varname = varname,
        dim1name = dim1name,
        conversion_factor = conversion_factor,
        do_check_sums_equal_1 = do_check_sums_equal_1,
        data_shape = copy(data_shape),
        allow_nodata = allow_nodata,
        data_on_file = false,
        # Allocate space for data (product of data_shape).
        data_at_tlower = zeros(Float64, prod(data_shape)),
        time_index_lower = 0,
    )

    # Read first set of data.
    var.time_index_lower = get_time_index_lower(dyn_file.time_info)
    _read_variable!(var, var.time_index_lower)
    return var
end

# ----------------------------------------------------------------------
# read_variable (dynVarMod :: read_variable / read_variable_{1,2}d)
# ----------------------------------------------------------------------
# Read a single time slice `nt` of the variable from the file into
# `var.data_at_tlower` (a flat 1-d buffer). Applies the conversion factor and,
# for 2-d variables, optionally checks that sums equal 1. Honors allow_nodata.
function _read_variable!(var::DynVarTimeUninterp, nt::Int)
    ndims = length(var.data_shape)
    @assert 1 <= ndims <= DYN_VAR_MAX_DIMS "read_variable can only handle 1 or 2 dimensions"

    ds = var.dyn_file.ds
    die_on_error = !var.allow_nodata

    if !haskey(ds, var.varname)
        # Variable not on file.
        if die_on_error
            error("ERROR: $(var.varname) NOT on file")
        else
            @warn "WARNING: $(var.varname) NOT on file set to zero" maxlog = 1
            fill!(var.data_at_tlower, 0.0)
            var.data_on_file = false
            return nothing
        end
    end

    var.data_on_file = true

    # Read the nt-th time slice. The file variable is laid out with the spatial
    # dimension(s) first and 'time' last (CLM/NetCDF convention), so we slice the
    # trailing time index. The result is reshaped/flattened to data_at_tlower.
    ncvar = ds[var.varname]
    rawslice = _read_time_slice(ncvar, nt, ndims, var.data_shape)
    arrayl = Float64.(replace(rawslice, missing => NaN))

    # Apply conversion: data are DIVIDED by conversion_factor.
    arrayl ./= var.conversion_factor

    if ndims == 2 && var.do_check_sums_equal_1
        _check_sums_equal_1(arrayl, var.varname)
    end

    # Flatten into the 1-d buffer (column-major; matches reshape semantics).
    var.data_at_tlower .= vec(arrayl)
    return nothing
end

# Read the nt-th time slice from a NetCDF variable whose trailing dimension is
# 'time'. `ndims` is the number of dimensions of the in-memory array
# (`data_shape`); `data_shape` is that shape.
#
# The gridcell axis may be stored on file EITHER as a single 'lndgrid' dimension
# OR — as in every real CTSM surface / landuse_timeseries dataset — as a pair of
# ('lsmlon', 'lsmlat') dimensions. Fortran reads both through the same `grlnd`
# decomposition (ncdio_pio collapses the horizontal dims into the gridcell axis),
# so we do the same: any extra non-time dimensions are folded into the gridcell
# axis by reshaping the slice to `data_shape`.
#
# NB: this requires the gridcell dims to be the SLOWEST-varying non-time dims,
# i.e. `PCT_NAT_PFT(time, lsmlat, lsmlon, natpft)` in ncdump order — the same
# natpft-fastest convention the surfdata reader already assumes (and that Bow's
# real surfdata uses). Rejecting these files outright meant CLM.jl could not read
# ANY real CTSM flanduse_timeseries.
function _read_time_slice(ncvar, nt::Int, ndims::Int, data_shape)
    nd_total = Base.ndims(ncvar)
    expected = prod(data_shape)

    slice = if nd_total == ndims && length(ncvar) == expected
        # No explicit time dimension on the variable (single time sample).
        Array(ncvar)
    elseif nd_total >= ndims
        # Slice the trailing time index; keep every other dimension whole.
        idx = Any[Colon() for _ in 1:(nd_total - 1)]
        push!(idx, nt)
        Array(view(ncvar, idx...))
    else
        error("_read_time_slice: variable has $nd_total dims; expected at least $ndims")
    end

    length(slice) == expected || error(
        "_read_time_slice: time slice has $(length(slice)) elements " *
        "($(size(slice))); expected $expected for data_shape $(Tuple(data_shape)). " *
        "The gridcell dimension(s) must be the slowest-varying non-time dims.")

    return length(size(slice)) == ndims && size(slice) == Tuple(data_shape) ?
           slice : reshape(slice, Tuple(data_shape))
end

# check_sums_equal_1 (surfrdUtilsMod): for a 2-d array, verify that the sum over
# dimension 1 equals 1 for each column. arrayl here is shaped (data_shape...).
function _check_sums_equal_1(arrayl::AbstractArray, varname::String)
    sums = sum(arrayl, dims = 1)
    if any(s -> abs(s - 1.0) > 1e-13, sums)
        error("check_sums_equal_1 ERROR: sums not equal to 1 for variable $varname")
    end
    return nothing
end

# ----------------------------------------------------------------------
# read_data_if_needed! (dynVarTimeUninterpMod :: read_data_if_needed)
# ----------------------------------------------------------------------
# Determine if new data need to be read from the file; if so, read them. We need
# to read new data if the current lower time index on dyn_file disagrees with the
# lower index for which we currently have stored data.
function read_data_if_needed!(var::DynVarTimeUninterp)
    time_index_lower_cur = get_time_index_lower(var.dyn_file.time_info)
    if time_index_lower_cur != var.time_index_lower
        _read_variable!(var, time_index_lower_cur)
        var.time_index_lower = time_index_lower_cur
    end
    return nothing
end

# ----------------------------------------------------------------------
# get_current_data_{1,2}d (dynVarTimeUninterpMod :: get_current_data_{1,2}d)
# ----------------------------------------------------------------------
"""
    get_current_data_1d(var::DynVarTimeUninterp) -> Vector{Float64}

Get the current value of the data for a 1-d variable. If necessary, new data are
read from the file. Should be called once per time step AFTER calling
`set_current_year!` on the underlying `DynFile`.
"""
function get_current_data_1d(var::DynVarTimeUninterp)
    ndims = length(var.data_shape)
    @assert ndims == 1 "get_current_data_1d ERROR: # dims of output must match ndims"
    read_data_if_needed!(var)
    return copy(var.data_at_tlower)
end

"""
    get_current_data_2d(var::DynVarTimeUninterp) -> Matrix{Float64}

Get the current value of the data for a 2-d variable, reshaped to `data_shape`. If
necessary, new data are read from the file. Should be called once per time step
AFTER calling `set_current_year!` on the underlying `DynFile`.
"""
function get_current_data_2d(var::DynVarTimeUninterp)
    ndims = length(var.data_shape)
    @assert ndims == 2 "get_current_data_2d ERROR: # dims of output must match ndims"
    read_data_if_needed!(var)
    return reshape(copy(var.data_at_tlower), var.data_shape[1], var.data_shape[2])
end
