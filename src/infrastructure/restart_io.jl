# ==========================================================================
# Restart I/O — write/read the model's prognostic state to a CLM-style
# restart NetCDF file and round-trip it back.
#
# Ported from: src/main/restFileMod.F90 (file-level driver only).
#
# restFileMod's job in Fortran is to open the NetCDF, define the subgrid
# dimensions (gridcell/landunit/column/pft + the level dims levgrnd/levsno/
# levtot/levlak/levdcmp), then call each type's `*_Restart` subroutine which
# in turn issues `restartvar` calls (define on write, read on read). We mirror
# that structure here: `write_restart` / `read_restart!` are the file-level
# driver, and a per-variable registry plays the role of the scattered
# `restartvar` calls — each entry names the CLM variable, its subgrid level and
# vertical dimension, and a getter/setter onto the live `CLMInstances` tree.
#
# Variable names and dimension names follow CLM (T_SOISNO/T_GRND/T_VEG; dims
# column/pft/levtot/levlak/levgrnd/levdcmp/levsoi) where practical, for
# traceability back to the Fortran restart files.
#
# Scope: the prognostic biogeophysical state (temperature, soil/snow water,
# snow-layer geometry, canopy water, lake state, canopy state) plus, when CN
# is active, the soil-BGC and CN-vegetation carbon/nitrogen pools. This is a
# faithful subset — history/diagnostic-only fields are intentionally omitted.
#
# Public functions:
#   write_restart   — create the NetCDF and write the prognostic state
#   read_restart!   — read it back into a (fresh) CLMInstances
# ==========================================================================

const _RESTART_FILL = -9999.0

"""
    RestartVarDef

Definition of a single restart variable. `dims` is the rank of the array
(1=Vector, 2=Matrix, 3=3-D); `level` is the subgrid dimension the first axis
lives on ("column" or "patch"); `levdim` is the NetCDF name for the (optional)
vertical/pool axes — empty for rank-1 variables. `getter`/`setter!` read/write
the corresponding field on the `CLMInstances` tree.
"""
struct RestartVarDef
    name::String
    dims::Int
    level::String          # "column" or "patch"
    levdim::String         # NetCDF level-dim name (or "")
    leveldim2::String      # second level dim for rank-3 vars (or "")
    getter::Function       # (inst) -> AbstractArray
    setter!::Function      # (inst, data) -> nothing
end

# convenience constructors -----------------------------------------------------
_rv1(name, level, get, set) = RestartVarDef(name, 1, level, "", "", get, set)
_rv2(name, level, levdim, get, set) = RestartVarDef(name, 2, level, levdim, "", get, set)
_rv3(name, level, levdim, levdim2, get, set) =
    RestartVarDef(name, 3, level, levdim, levdim2, get, set)

"""
    _restart_registry_biogeophys() -> Vector{RestartVarDef}

Prognostic biogeophysical state: temperature, soil/snow water, snow-layer
geometry, canopy water, lake and canopy state. Always written.
"""
function _restart_registry_biogeophys()
    return RestartVarDef[
        # --- Temperature state ---
        _rv2("T_SOISNO", "column", "levtot",
            i -> i.temperature.t_soisno_col,
            (i, d) -> i.temperature.t_soisno_col .= d),
        _rv1("T_GRND", "column",
            i -> i.temperature.t_grnd_col,
            (i, d) -> i.temperature.t_grnd_col .= d),
        _rv2("T_LAKE", "column", "levlak",
            i -> i.temperature.t_lake_col,
            (i, d) -> i.temperature.t_lake_col .= d),
        _rv1("T_H2OSFC", "column",
            i -> i.temperature.t_h2osfc_col,
            (i, d) -> i.temperature.t_h2osfc_col .= d),
        _rv1("T_VEG", "patch",
            i -> i.temperature.t_veg_patch,
            (i, d) -> i.temperature.t_veg_patch .= d),
        _rv1("T_REF2M", "patch",
            i -> i.temperature.t_ref2m_patch,
            (i, d) -> i.temperature.t_ref2m_patch .= d),

        # --- Water state ---
        _rv2("H2OSOI_LIQ", "column", "levtot",
            i -> i.water.waterstatebulk_inst.ws.h2osoi_liq_col,
            (i, d) -> i.water.waterstatebulk_inst.ws.h2osoi_liq_col .= d),
        _rv2("H2OSOI_ICE", "column", "levtot",
            i -> i.water.waterstatebulk_inst.ws.h2osoi_ice_col,
            (i, d) -> i.water.waterstatebulk_inst.ws.h2osoi_ice_col .= d),
        _rv1("H2OSNO", "column",
            i -> i.water.waterstatebulk_inst.ws.h2osno_no_layers_col,
            (i, d) -> i.water.waterstatebulk_inst.ws.h2osno_no_layers_col .= d),
        _rv1("H2OSFC", "column",
            i -> i.water.waterstatebulk_inst.ws.h2osfc_col,
            (i, d) -> i.water.waterstatebulk_inst.ws.h2osfc_col .= d),
        _rv1("WA", "column",
            i -> i.water.waterstatebulk_inst.ws.wa_col,
            (i, d) -> i.water.waterstatebulk_inst.ws.wa_col .= d),
        _rv1("INT_SNOW", "column",
            i -> i.water.waterstatebulk_inst.int_snow_col,
            (i, d) -> i.water.waterstatebulk_inst.int_snow_col .= d),

        # --- Soil hydrology ---
        _rv1("ZWT", "column",
            i -> i.soilhydrology.zwt_col,
            (i, d) -> i.soilhydrology.zwt_col .= d),
        _rv1("ZWT_PERCH", "column",
            i -> i.soilhydrology.zwt_perched_col,
            (i, d) -> i.soilhydrology.zwt_perched_col .= d),

        # --- Snow state (counts + diagnostics + geometry) ---
        _rv1("SNLSNO", "column",
            i -> Float64.(i.column.snl),
            # snl is an integer count; its natural missing marker is ISPVAL, which
            # equals the float fill (-9999) and is restored to NaN on read — map
            # NaN back to ISPVAL rather than letting round(Int, NaN) throw.
            (i, d) -> i.column.snl .= (x -> isnan(x) ? ISPVAL : round(Int, x)).(d)),
        _rv1("SNOW_DEPTH", "column",
            i -> i.water.waterdiagnosticbulk_inst.snow_depth_col,
            (i, d) -> i.water.waterdiagnosticbulk_inst.snow_depth_col .= d),
        _rv1("frac_sno", "column",
            i -> i.water.waterdiagnosticbulk_inst.frac_sno_col,
            (i, d) -> i.water.waterdiagnosticbulk_inst.frac_sno_col .= d),
        _rv1("frac_sno_eff", "column",
            i -> i.water.waterdiagnosticbulk_inst.frac_sno_eff_col,
            (i, d) -> i.water.waterdiagnosticbulk_inst.frac_sno_eff_col .= d),
        _rv2("DZSNO", "column", "levtot",
            i -> i.column.dz,
            (i, d) -> i.column.dz .= d),
        _rv2("ZSNO", "column", "levtot",
            i -> i.column.z,
            (i, d) -> i.column.z .= d),
        _rv2("ZISNO", "column", "leviface",
            i -> i.column.zi,
            (i, d) -> i.column.zi .= d),

        # --- Canopy water ---
        _rv1("LIQCAN", "patch",
            i -> i.water.waterstatebulk_inst.ws.liqcan_patch,
            (i, d) -> i.water.waterstatebulk_inst.ws.liqcan_patch .= d),
        _rv1("SNOCAN", "patch",
            i -> i.water.waterstatebulk_inst.ws.snocan_patch,
            (i, d) -> i.water.waterstatebulk_inst.ws.snocan_patch .= d),

        # --- Lake state ---
        _rv2("LAKE_ICEFRAC", "column", "levlak",
            i -> i.lakestate.lake_icefrac_col,
            (i, d) -> i.lakestate.lake_icefrac_col .= d),
        _rv1("SAVEDTKE1", "column",
            i -> i.lakestate.savedtke1_col,
            (i, d) -> i.lakestate.savedtke1_col .= d),

        # --- Canopy state ---
        _rv1("TLAI", "patch",
            i -> i.canopystate.tlai_patch,
            (i, d) -> i.canopystate.tlai_patch .= d),
        _rv1("ELAI", "patch",
            i -> i.canopystate.elai_patch,
            (i, d) -> i.canopystate.elai_patch .= d),
        _rv1("HTOP", "patch",
            i -> i.canopystate.htop_patch,
            (i, d) -> i.canopystate.htop_patch .= d),
    ]
end

"""
    _restart_registry_cn() -> Vector{RestartVarDef}

CN prognostic state: soil-BGC vertically-resolved carbon/nitrogen pools and
mineral N, plus CN-vegetation carbon/nitrogen pools. Written only when CN is
active. Names follow CLM (`decomp_cpools_vr`, `sminn_vr`, `leafc`, …).
"""
function _restart_registry_cn()
    return RestartVarDef[
        # --- Soil BGC carbon (vertically + pool resolved) ---
        _rv3("decomp_cpools_vr", "column", "levdcmp", "ndecomp_pools",
            i -> i.soilbiogeochem_carbonstate.decomp_cpools_vr_col,
            (i, d) -> i.soilbiogeochem_carbonstate.decomp_cpools_vr_col .= d),
        _rv2("ctrunc_vr", "column", "levdcmp",
            i -> i.soilbiogeochem_carbonstate.ctrunc_vr_col,
            (i, d) -> i.soilbiogeochem_carbonstate.ctrunc_vr_col .= d),

        # --- Soil BGC nitrogen ---
        _rv3("decomp_npools_vr", "column", "levdcmp", "ndecomp_pools",
            i -> i.soilbiogeochem_nitrogenstate.decomp_npools_vr_col,
            (i, d) -> i.soilbiogeochem_nitrogenstate.decomp_npools_vr_col .= d),
        _rv2("sminn_vr", "column", "levdcmp",
            i -> i.soilbiogeochem_nitrogenstate.sminn_vr_col,
            (i, d) -> i.soilbiogeochem_nitrogenstate.sminn_vr_col .= d),
        _rv2("smin_no3_vr", "column", "levdcmp",
            i -> i.soilbiogeochem_nitrogenstate.smin_no3_vr_col,
            (i, d) -> i.soilbiogeochem_nitrogenstate.smin_no3_vr_col .= d),
        _rv2("smin_nh4_vr", "column", "levdcmp",
            i -> i.soilbiogeochem_nitrogenstate.smin_nh4_vr_col,
            (i, d) -> i.soilbiogeochem_nitrogenstate.smin_nh4_vr_col .= d),

        # --- CN vegetation carbon pools ---
        _rv1("leafc", "patch",
            i -> i.bgc_vegetation.cnveg_carbonstate_inst.leafc_patch,
            (i, d) -> i.bgc_vegetation.cnveg_carbonstate_inst.leafc_patch .= d),
        _rv1("frootc", "patch",
            i -> i.bgc_vegetation.cnveg_carbonstate_inst.frootc_patch,
            (i, d) -> i.bgc_vegetation.cnveg_carbonstate_inst.frootc_patch .= d),
        _rv1("livestemc", "patch",
            i -> i.bgc_vegetation.cnveg_carbonstate_inst.livestemc_patch,
            (i, d) -> i.bgc_vegetation.cnveg_carbonstate_inst.livestemc_patch .= d),
        _rv1("deadstemc", "patch",
            i -> i.bgc_vegetation.cnveg_carbonstate_inst.deadstemc_patch,
            (i, d) -> i.bgc_vegetation.cnveg_carbonstate_inst.deadstemc_patch .= d),
        _rv1("cpool", "patch",
            i -> i.bgc_vegetation.cnveg_carbonstate_inst.cpool_patch,
            (i, d) -> i.bgc_vegetation.cnveg_carbonstate_inst.cpool_patch .= d),

        # --- CN vegetation nitrogen pools ---
        _rv1("leafn", "patch",
            i -> i.bgc_vegetation.cnveg_nitrogenstate_inst.leafn_patch,
            (i, d) -> i.bgc_vegetation.cnveg_nitrogenstate_inst.leafn_patch .= d),
        _rv1("frootn", "patch",
            i -> i.bgc_vegetation.cnveg_nitrogenstate_inst.frootn_patch,
            (i, d) -> i.bgc_vegetation.cnveg_nitrogenstate_inst.frootn_patch .= d),
        _rv1("npool", "patch",
            i -> i.bgc_vegetation.cnveg_nitrogenstate_inst.npool_patch,
            (i, d) -> i.bgc_vegetation.cnveg_nitrogenstate_inst.npool_patch .= d),
    ]
end

# Replace non-finite entries with the NetCDF fill value (in place, on a copy).
function _clean_for_write!(a::AbstractArray)
    @inbounds for i in eachindex(a)
        isfinite(a[i]) || (a[i] = _RESTART_FILL)
    end
    return a
end

# Restore fill values back to NaN after read (in place).
function _restore_after_read!(a::AbstractArray)
    @inbounds for i in eachindex(a)
        a[i] == _RESTART_FILL && (a[i] = NaN)
    end
    return a
end

# Ensure a NetCDF dim of name `nm` exists with size `n`; if it exists with a
# different size, fall back to a unique per-variable name and define that.
function _ensure_dim!(ds, nm::AbstractString, n::Int, fallback::AbstractString)
    if haskey(ds.dim, nm)
        ds.dim[nm] == n && return nm
        # size clash → use a variable-unique dim name
        if !haskey(ds.dim, fallback)
            defDim(ds, fallback, n)
        end
        return fallback
    end
    defDim(ds, nm, n)
    return nm
end

"""
    write_restart(inst, filename; bounds, use_cn=false, time=nothing)

Create a CLM-style restart NetCDF at `filename` and write the prognostic state
held in `inst`. Mirrors `restFileMod::restFile_write`: defines the subgrid and
level dimensions, then writes each registered variable (the role of the
per-type `restartvar` calls).

- `bounds`  — `BoundsType` giving `endc` (#columns) and `endp` (#patches).
- `use_cn`  — also write the CN carbon/nitrogen pools.
- `time`    — optional valid time, stored as a scalar `timemgr_rst_curr_date`.
"""
function write_restart(inst::CLMInstances, filename::String;
                       bounds::BoundsType,
                       use_cn::Bool = false,
                       time::Union{DateTime,Nothing} = nothing)
    nc = bounds.endc
    np = bounds.endp

    registry = _restart_registry_biogeophys()
    use_cn && append!(registry, _restart_registry_cn())

    ds = NCDataset(filename, "c")
    try
        # Subgrid dimensions (restFileMod defines column/pft; we add what we use).
        defDim(ds, "column", nc)
        defDim(ds, "pft", np)

        if time !== nothing
            tvar = defVar(ds, "timemgr_rst_curr_date", Float64, ())
            tvar[] = Dates.value(time - DateTime(2000, 1, 1)) / (1000 * 86400)
        end

        for rv in registry
            data = try
                rv.getter(inst)
            catch
                continue
            end
            (data === nothing || length(data) == 0) && continue

            dim1 = rv.level == "column" ? "column" : "pft"

            if rv.dims == 1
                v = defVar(ds, rv.name, Float64, (dim1,))
                clean = _clean_for_write!(collect(Float64, data))
                v[:] = clean
            elseif rv.dims == 2
                n2 = size(data, 2)
                d2 = _ensure_dim!(ds, rv.levdim, n2, "lev_$(rv.name)")
                v = defVar(ds, rv.name, Float64, (dim1, d2))
                clean = _clean_for_write!(collect(Float64, data))
                v[:, :] = clean
            else # rv.dims == 3
                n2 = size(data, 2); n3 = size(data, 3)
                d2 = _ensure_dim!(ds, rv.levdim,  n2, "lev_$(rv.name)")
                d3 = _ensure_dim!(ds, rv.leveldim2, n3, "pool_$(rv.name)")
                v = defVar(ds, rv.name, Float64, (dim1, d2, d3))
                clean = _clean_for_write!(collect(Float64, data))
                v[:, :, :] = clean
            end
        end
    finally
        close(ds)
    end
    return nothing
end

"""
    read_restart!(inst, filename; bounds, use_cn=false)

Read the prognostic state from a restart NetCDF written by [`write_restart`]
back into `inst`, overwriting each registered field. Mirrors
`restFileMod::restFile_read`. Variables absent from the file are left untouched;
fill values are restored to `NaN`.
"""
function read_restart!(inst::CLMInstances, filename::String;
                       bounds::BoundsType,
                       use_cn::Bool = false)
    isfile(filename) || error("Restart file not found: $filename")

    registry = _restart_registry_biogeophys()
    use_cn && append!(registry, _restart_registry_cn())

    NCDataset(filename, "r") do ds
        for rv in registry
            haskey(ds, rv.name) || continue
            raw = Array(ds[rv.name])
            data = Float64.(replace(raw, missing => NaN))
            _restore_after_read!(data)
            try
                rv.setter!(inst, data)
            catch e
                @warn "Skipping restart variable $(rv.name): $e" maxlog = 1
            end
        end
    end
    return nothing
end

# --- Backwards-compatible aliases (old arg order: filepath, inst, bounds) -----
# The earlier API was `write_restart!(filepath, inst, bounds)`. Keep thin
# shims so existing callers/scripts continue to work.
function write_restart!(filepath::String, inst::CLMInstances, bounds::BoundsType;
                        use_cn::Bool = false, time::Union{DateTime,Nothing} = nothing)
    return write_restart(inst, filepath; bounds = bounds, use_cn = use_cn, time = time)
end

function read_restart!(filepath::String, inst::CLMInstances, bounds::BoundsType;
                       use_cn::Bool = false)
    return read_restart!(inst, filepath; bounds = bounds, use_cn = use_cn)
end
