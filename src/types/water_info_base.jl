# ==========================================================================
# Ported from: src/biogeophys/WaterInfoBaseType.F90
# Base type for working with information describing a given water instance
# (bulk or tracer), such as building history and restart field names.
# ==========================================================================

"""
    WaterInfoBaseType

Abstract base type for water information (bulk or tracer).

Subtypes must have a `ratio::Float64` field and implement:
- `get_name(info)::String`
- `fname(info, basename::String)::String`
- `lname(info, basename::String)::String`
- `is_communicated_with_coupler(info)::Bool`
- `is_included_in_consistency_check(info)::Bool`

Ported from `water_info_base_type` in `WaterInfoBaseType.F90`.
"""
abstract type WaterInfoBaseType end

# ---------------------------------------------------------------------------
# Deferred interface — subtypes must implement these
# ---------------------------------------------------------------------------

"""
    get_name(info::WaterInfoBaseType) -> String

Get the name of this tracer (or bulk).
"""
function get_name(::WaterInfoBaseType)
    error("get_name must be implemented by subtypes of WaterInfoBaseType")
end

"""
    fname(info::WaterInfoBaseType, basename::String) -> String

Get a history/restart field name for this tracer (or bulk).
`basename` gives the base name of the history/restart field.
"""
function fname(::WaterInfoBaseType, ::String)
    error("fname must be implemented by subtypes of WaterInfoBaseType")
end

"""
    lname(info::WaterInfoBaseType, basename::String) -> String

Get a history/restart long name for this tracer (or bulk).
`basename` gives the base name of the history/restart long name.
"""
function lname(::WaterInfoBaseType, ::String)
    error("lname must be implemented by subtypes of WaterInfoBaseType")
end

"""
    is_communicated_with_coupler(info::WaterInfoBaseType) -> Bool

Returns true if this tracer is received from and sent to the coupler.
Returns false if this tracer is just used internally in CTSM, and is
set to some fixed ratio times the bulk water.
"""
function is_communicated_with_coupler(::WaterInfoBaseType)
    error("is_communicated_with_coupler must be implemented by subtypes of WaterInfoBaseType")
end

"""
    is_included_in_consistency_check(info::WaterInfoBaseType) -> Bool

Returns true if this tracer is included in water consistency checks.
"""
function is_included_in_consistency_check(::WaterInfoBaseType)
    error("is_included_in_consistency_check must be implemented by subtypes of WaterInfoBaseType")
end

# ---------------------------------------------------------------------------
# Concrete methods
# ---------------------------------------------------------------------------

"""
    set_metadata!(info::WaterInfoBaseType, ratio::Float64)

Set the ratio metadata on a water info instance.
"""
function set_metadata!(info::WaterInfoBaseType, ratio::Real)
    info.ratio = ratio
    return nothing
end

"""
    get_ratio(info::WaterInfoBaseType) -> Float64

Get the ratio from a water info instance.
"""
function get_ratio(info::WaterInfoBaseType)
    return info.ratio
end
