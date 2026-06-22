# ==========================================================================
# init_interp — start a run from a finidat restart that is NOT on the run's
# exact grid.
#
# Ported (a faithful, single-column-capable subset) from:
#   src/init_interp/initInterpMindist.F90        (subgrid source→dest matching)
#   src/init_interp/initInterp1dData.F90          (1-d field copy via sgridindex)
#   src/init_interp/initInterpMultilevelInterp.F90 (vertical level remap)
#
# The full Fortran init_interp handles arbitrary spatial regridding across many
# gridcells (nearest-neighbour in lat/lon among same-type subgrid elements) plus
# vertical interpolation of column profiles when the soil/snow layer structure
# differs. This port implements the pieces that matter for the single-point /
# few-column case that CLM.jl currently exercises:
#
#   * SUBGRID MATCHING (set_mindist / set_single_match): for each destination
#     column/patch, pick the source element of the SAME TYPE (column.itype /
#     patch.itype) that is nearest in lat/lon. With a single matching type this
#     reduces to "same type, copy".
#   * VERTICAL INTERP (interp_multilevel / interp_onelevel): when the source
#     restart has a different number of soil/snow levels than the destination,
#     linearly interpolate each column profile from the source coordinate grid
#     (the source `z` levels) onto the destination `z` levels, with copy/clamp
#     at the ends — mirroring initInterpMultilevelInterp.
#
# GATE: if the source and destination grids match exactly (same column/patch
# counts and same number of levels), init_interp! reduces to read_restart!.
#
# DEFERRED (noted, not needed for single-point): true multi-gridcell spatial
# regridding with the glacier/topo tie-breaking, fill_missing_with_natveg, and
# the level_class soil-vs-bedrock partitioning. The vertical remap here treats
# all levels as a single class (initInterpMultilevelInterp's no-levclasses
# constructor), which is what the soil/snow profiles use.
# ==========================================================================

"""
    _interp_onelevel(coord_dest, coords_source, data_source) -> Float64

Linear interpolation of a single destination level from a source profile,
mirroring `initInterpMultilevelInterp::interp_onelevel` with
`scale_by_thickness = .false.` (the common case for intensive quantities such
as temperature and volumetric water).

`coords_source` must be monotonically increasing. Below the first / above the
last source coordinate, the nearest source value is copied (clamp). An exact
coordinate match (within roundoff) copies that level. Otherwise the value is a
linear blend of the bracketing source levels.
"""
function _interp_onelevel(coord_dest::Float64,
                          coords_source::AbstractVector{Float64},
                          data_source::AbstractVector{Float64})
    eps = 1.0e-13
    n = length(coords_source)
    n == 0 && return NaN  # no source info → caller keeps original value

    if coord_dest < coords_source[1]
        return data_source[1]                      # copy first (clamp low)
    elseif coord_dest > coords_source[n]
        return data_source[n]                      # copy last (clamp high)
    end

    # exact match within roundoff
    for lev in 1:n
        if abs(coord_dest - coords_source[lev]) < eps
            return data_source[lev]
        end
    end

    # bracketing interval
    for lev in 1:(n - 1)
        if coord_dest > coords_source[lev] && coord_dest < coords_source[lev + 1]
            c0 = coords_source[lev]
            c1 = coords_source[lev + 1]
            return data_source[lev + 1] * (coord_dest - c0) / (c1 - c0) +
                   data_source[lev]     * (c1 - coord_dest) / (c1 - c0)
        end
    end

    return NaN  # unreachable for monotone input bracketed above
end

_is_missing(v) = isnan(v) || v == SPVAL

"""
    _interp_profile(prof_source, z_source, z_dest) -> Vector{Float64}

Vertically remap one column profile `prof_source` (defined on the source
coordinate levels `z_source`) onto the destination levels `z_dest`. Missing
source levels are dropped before interpolation, mirroring the "static spval"
rule. Returns a vector of length `length(z_dest)`.
"""
function _interp_profile(prof_source::AbstractVector{Float64},
                         z_source::AbstractVector{Float64},
                         z_dest::AbstractVector{Float64})
    # pack out missing source levels (data or coordinate)
    keep = [(!_is_missing(prof_source[k]) && !_is_missing(z_source[k]))
            for k in eachindex(prof_source)]
    cs = Float64[z_source[k] for k in eachindex(z_source) if keep[k]]
    ds = Float64[prof_source[k] for k in eachindex(prof_source) if keep[k]]

    out = Vector{Float64}(undef, length(z_dest))
    for j in eachindex(z_dest)
        out[j] = _interp_onelevel(z_dest[j], cs, ds)
    end
    return out
end

"""
    _build_sgridindex(src_itype, dst_itype, src_lat, src_lon, dst_lat, dst_lon) -> Vector{Int}

For each destination element, find the source element of the same `itype` that
is nearest in lat/lon (the `set_mindist` metric). Returns a 1-based index into
the source arrays, or 0 if no same-type source exists (the field is then left
untouched, as in `interp_1d_data` with `sgridindex <= 0`).

This is the single-/few-point reduction of `initInterpMindist::set_mindist`:
same-type-then-nearest. lat/lon default to 0 when not available, in which case
all same-type sources are equidistant and the first is chosen (the natural
"single match" behaviour).
"""
function _build_sgridindex(src_itype::AbstractVector{<:Integer},
                           dst_itype::AbstractVector{<:Integer},
                           src_lat::AbstractVector{Float64},
                           src_lon::AbstractVector{Float64},
                           dst_lat::AbstractVector{Float64},
                           dst_lon::AbstractVector{Float64})
    ni = length(src_itype)
    no = length(dst_itype)
    idx = zeros(Int, no)
    for o in 1:no
        nmin = 0
        distmin = Inf
        for i in 1:ni
            src_itype[i] == dst_itype[o] || continue
            dy = abs(dst_lat[o] - src_lat[i]) * RE
            coslat = 0.5 * (cos(dst_lat[o]) + cos(src_lat[i]))
            dx = abs(dst_lon[o] - src_lon[i]) * RE * coslat
            dist = dx * dx + dy * dy
            if dist < distmin
                distmin = dist
                nmin = i
            end
        end
        idx[o] = nmin
    end
    return idx
end

# Subgrid metadata (types + coords) extracted from a live inst, for matching.
struct _SubgridMeta
    col_itype::Vector{Int}
    col_lat::Vector{Float64}
    col_lon::Vector{Float64}
    pft_itype::Vector{Int}
    pft_lat::Vector{Float64}
    pft_lon::Vector{Float64}
end

# Pull lat/lon for each column/patch from the gridcell they map to (radians).
# Falls back to zeros when gridcell coords are unavailable (single-point runs).
function _subgrid_meta(inst::CLMInstances, bounds::BoundsType)
    nc = bounds.endc
    np = bounds.endp
    grc = inst.gridcell
    have_grc = grc !== nothing && length(grc.lat) > 0

    col_lat = zeros(Float64, nc); col_lon = zeros(Float64, nc)
    col_it  = zeros(Int, nc)
    for c in 1:nc
        c <= length(inst.column.itype) && (col_it[c] = inst.column.itype[c])
        if have_grc && c <= length(inst.column.gridcell)
            g = inst.column.gridcell[c]
            if 1 <= g <= length(grc.lat)
                col_lat[c] = grc.lat[g]; col_lon[c] = grc.lon[g]
            end
        end
    end

    pft_lat = zeros(Float64, np); pft_lon = zeros(Float64, np)
    pft_it  = zeros(Int, np)
    for p in 1:np
        p <= length(inst.patch.itype) && (pft_it[p] = inst.patch.itype[p])
        if have_grc && p <= length(inst.patch.gridcell)
            g = inst.patch.gridcell[p]
            if 1 <= g <= length(grc.lat)
                pft_lat[p] = grc.lat[g]; pft_lon[p] = grc.lon[g]
            end
        end
    end
    return _SubgridMeta(col_it, col_lat, col_lon, pft_it, pft_lat, pft_lon)
end

# Read every restart variable from `filename` into a Dict{name => Array}, with
# fill values restored to NaN (the same cleanup read_restart! applies).
function _read_restart_raw(filename::String, registry::Vector{RestartVarDef})
    out = Dict{String,Array{Float64}}()
    NCDataset(filename, "r") do ds
        for rv in registry
            haskey(ds, rv.name) || continue
            raw = Array(ds[rv.name])
            data = Float64.(replace(raw, missing => NaN))
            _restore_after_read!(data)
            out[rv.name] = data
        end
    end
    return out
end

"""
    init_interp!(inst, finidat; bounds, use_cn=false, src_meta=nothing)

Initialise `inst` from the restart file `finidat`, interpolating across grids
when the restart's subgrid layout or vertical level structure differs from
`inst`'s. Mirrors CTSM's `initInterpMod` for the single-point / few-column case.

Behaviour:
  * If the restart's subgrid counts AND vertical-level counts match `inst`
    exactly (the same-grid case), delegates to [`read_restart!`] (identity).
  * Otherwise: for each destination column/patch, matches a source element of
    the same type (nearest in lat/lon) and copies its restart values; rank-2/3
    column fields whose level count differs are vertically interpolated onto the
    destination `z` levels (read from the restart's `ZSNO` variable).

`src_meta` (optional) supplies the SOURCE subgrid types/coords for matching when
they differ from the destination. When omitted, the source is assumed to share
the destination's subgrid layout (1:1 by index) — the common single-point case.
"""
function init_interp!(inst::CLMInstances, finidat::String;
                      bounds::BoundsType,
                      use_cn::Bool = false,
                      src_meta::Union{_SubgridMeta,Nothing} = nothing)
    isfile(finidat) || error("init_interp: finidat not found: $finidat")

    registry = _restart_registry_biogeophys()
    use_cn && append!(registry, _restart_registry_cn())

    raw = _read_restart_raw(finidat, registry)

    # Destination subgrid sizes / level structure.
    nc_d = bounds.endc
    np_d = bounds.endp
    dst_meta = _subgrid_meta(inst, bounds)

    # Probe ZSNO (the `z` coordinate) for the source level count; T_SOISNO etc.
    # share that level count. ZSNO is the vertical coordinate for interpolation.
    src_z = get(raw, "ZSNO", nothing)
    nc_s = src_z === nothing ? nc_d : size(src_z, 1)
    nlev_s = src_z === nothing ? 0 : size(src_z, 2)
    nlev_d = size(inst.column.z, 2)

    same_subgrid = (src_meta === nothing) && (nc_s == nc_d) &&
                   (nc_s == length(dst_meta.col_itype))
    same_levels = (nlev_s == 0) || (nlev_s == nlev_d)

    if same_subgrid && same_levels
        # Exact-grid case → identity restart read.
        read_restart!(inst, finidat; bounds = bounds, use_cn = use_cn)
        return nothing
    end

    # Build source→dest index maps (column- and patch-level).
    #
    # When no source metadata is supplied AND the element counts match, the
    # subgrid layouts are identical (1:1 by index) — this is the common case of
    # the SAME single point/columns at a different vertical resolution. Mapping
    # by index here mirrors set_single_match's "same location, same type" result
    # and avoids collapsing all columns onto the first via degenerate (all-zero)
    # type/coordinate metadata. Otherwise fall back to same-type-then-nearest.
    np_s = _patch_count(raw, np_d)
    if src_meta === nothing && nc_s == nc_d
        col_sidx = collect(1:nc_d)
    else
        sm = src_meta === nothing ? dst_meta : src_meta
        col_sidx = _build_sgridindex(sm.col_itype[1:min(nc_s, length(sm.col_itype))],
                                     dst_meta.col_itype,
                                     sm.col_lat[1:min(nc_s, length(sm.col_lat))],
                                     sm.col_lon[1:min(nc_s, length(sm.col_lon))],
                                     dst_meta.col_lat, dst_meta.col_lon)
    end
    if src_meta === nothing && np_s == np_d
        pft_sidx = collect(1:np_d)
    else
        sm = src_meta === nothing ? dst_meta : src_meta
        pft_sidx = _build_sgridindex(sm.pft_itype[1:min(np_s, length(sm.pft_itype))],
                                     dst_meta.pft_itype,
                                     sm.pft_lat[1:min(np_s, length(sm.pft_lat))],
                                     sm.pft_lon[1:min(np_s, length(sm.pft_lon))],
                                     dst_meta.pft_lat, dst_meta.pft_lon)
    end

    # Destination per-column z levels (interp target coordinate).
    z_dest = inst.column.z

    for rv in registry
        haskey(raw, rv.name) || continue
        data = raw[rv.name]

        if rv.dims == 1
            sidx = rv.level == "column" ? col_sidx : pft_sidx
            _apply_1d!(inst, rv, data, sidx)
        elseif rv.dims == 2
            _apply_2d!(inst, rv, data, col_sidx, pft_sidx, src_z, z_dest, nlev_d)
        else # dims == 3 — pool axis copied as-is, level axis interpolated
            _apply_3d!(inst, rv, data, col_sidx, src_z, z_dest, nlev_d)
        end
    end
    return nothing
end

# Patch count from a representative rank-1 patch variable, else dest count.
function _patch_count(raw, np_d::Int)
    for nm in ("T_VEG", "T_REF2M", "TLAI", "ELAI", "HTOP", "leafc")
        haskey(raw, nm) && return length(raw[nm])
    end
    return np_d
end

# Copy a rank-1 field through the sgridindex (interp_1d_data_double).
function _apply_1d!(inst, rv::RestartVarDef, data::Array{Float64}, sidx::Vector{Int})
    no = length(sidx)
    out = fill(NaN, no)
    for o in 1:no
        ni = sidx[o]
        if ni > 0 && ni <= length(data)
            out[o] = _is_missing(data[ni]) ? NaN : data[ni]
        end
    end
    try
        rv.setter!(inst, out)
    catch e
        @warn "init_interp: skipping $(rv.name): $e" maxlog = 1
    end
    return nothing
end

# Destination level count for a registered field, read from the live field's
# shape (lake/decomp axes differ from the soil `z` axis, so we cannot assume a
# single global value).
function _dest_nlev(inst, rv::RestartVarDef, fallback::Int)
    try
        d = rv.getter(inst)
        return ndims(d) >= 2 ? size(d, 2) : fallback
    catch
        return fallback
    end
end

# Copy a rank-2 (column,level) or (patch,level) field, vertically interpolating
# soil/snow column profiles (those whose level axis is the `z` coordinate) when
# the level count differs. Non-soil axes (lake, decomp) are level-copied onto
# the destination field's own level count.
function _apply_2d!(inst, rv::RestartVarDef, data::Array{Float64},
                    col_sidx::Vector{Int}, pft_sidx::Vector{Int},
                    src_z, z_dest, nlev_z::Int)
    sidx = rv.level == "column" ? col_sidx : pft_sidx
    no = length(sidx)
    nlev_s = size(data, 2)
    nlev_out = _dest_nlev(inst, rv, nlev_s)
    out = fill(NaN, no, nlev_out)

    # Vertical interp only for soil/snow column fields: the source level count
    # must match the `z` coordinate axis, and the destination field must too.
    do_vinterp = rv.level == "column" && src_z !== nothing &&
                 nlev_s == size(src_z, 2) && nlev_out == nlev_z &&
                 nlev_s != nlev_out

    for o in 1:no
        ni = sidx[o]
        (ni > 0 && ni <= size(data, 1)) || continue
        if do_vinterp
            prof_s = Float64[data[ni, k] for k in 1:nlev_s]
            zs     = Float64[src_z[ni, k] for k in 1:nlev_s]
            zd     = Float64[z_dest[o, j] for j in 1:nlev_out]
            out[o, :] .= _interp_profile(prof_s, zs, zd)
        else
            ncopy = min(nlev_s, nlev_out)
            for k in 1:ncopy
                out[o, k] = _is_missing(data[ni, k]) ? NaN : data[ni, k]
            end
        end
    end
    try
        rv.setter!(inst, out)
    catch e
        @warn "init_interp: skipping $(rv.name): $e" maxlog = 1
    end
    return nothing
end

# Rank-3 (column,level,pool): interpolate each pool's profile independently
# when the level axis is the `z` coordinate; otherwise level-copy.
function _apply_3d!(inst, rv::RestartVarDef, data::Array{Float64},
                    col_sidx::Vector{Int}, src_z, z_dest, nlev_z::Int)
    no = length(col_sidx)
    nlev_s = size(data, 2)
    npool = size(data, 3)
    nlev_out = _dest_nlev(inst, rv, nlev_s)
    out = fill(NaN, no, nlev_out, npool)
    do_vinterp = src_z !== nothing && nlev_s == size(src_z, 2) &&
                 nlev_out == nlev_z && nlev_s != nlev_out

    for o in 1:no
        ni = col_sidx[o]
        (ni > 0 && ni <= size(data, 1)) || continue
        for q in 1:npool
            if do_vinterp
                prof_s = Float64[data[ni, k, q] for k in 1:nlev_s]
                zs     = Float64[src_z[ni, k] for k in 1:nlev_s]
                zd     = Float64[z_dest[o, j] for j in 1:nlev_out]
                out[o, :, q] .= _interp_profile(prof_s, zs, zd)
            else
                ncopy = min(nlev_s, nlev_out)
                for k in 1:ncopy
                    out[o, k, q] = _is_missing(data[ni, k, q]) ? NaN : data[ni, k, q]
                end
            end
        end
    end
    try
        rv.setter!(inst, out)
    catch e
        @warn "init_interp: skipping $(rv.name): $e" maxlog = 1
    end
    return nothing
end
