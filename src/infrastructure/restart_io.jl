# ==========================================================================
# Restart I/O — Save and restore CLM state for spinup
#
# Public functions:
#   write_restart!  — Write all state arrays to a NetCDF file
#   read_restart!   — Read state arrays from a restart file
# ==========================================================================

"""
    RestartVarDef

Definition of a single restart variable.
"""
struct RestartVarDef
    name::String
    dims::Int       # 1 = Vector, 2 = Matrix
    level::String   # "column" or "patch"
    getter::Function  # (inst) -> AbstractArray
    setter!::Function # (inst, data) -> nothing
end

"""
    _restart_registry() -> Vector{RestartVarDef}

Return the list of state variables saved/restored in restart files.
"""
function _restart_registry()
    return RestartVarDef[
        # --- Temperature state ---
        RestartVarDef("T_SOISNO", 2, "column",
            inst -> inst.temperature.t_soisno_col,
            (inst, d) -> inst.temperature.t_soisno_col .= d),
        RestartVarDef("T_GRND", 1, "column",
            inst -> inst.temperature.t_grnd_col,
            (inst, d) -> inst.temperature.t_grnd_col .= d),
        RestartVarDef("T_LAKE", 2, "column",
            inst -> inst.temperature.t_lake_col,
            (inst, d) -> inst.temperature.t_lake_col .= d),
        RestartVarDef("T_H2OSFC", 1, "column",
            inst -> inst.temperature.t_h2osfc_col,
            (inst, d) -> inst.temperature.t_h2osfc_col .= d),
        RestartVarDef("T_VEG", 1, "patch",
            inst -> inst.temperature.t_veg_patch,
            (inst, d) -> inst.temperature.t_veg_patch .= d),
        RestartVarDef("T_REF2M", 1, "patch",
            inst -> inst.temperature.t_ref2m_patch,
            (inst, d) -> inst.temperature.t_ref2m_patch .= d),
        RestartVarDef("T_VEG24", 1, "patch",
            inst -> inst.temperature.t_veg24_patch,
            (inst, d) -> inst.temperature.t_veg24_patch .= d),
        RestartVarDef("T_VEG240", 1, "patch",
            inst -> inst.temperature.t_veg240_patch,
            (inst, d) -> inst.temperature.t_veg240_patch .= d),
        RestartVarDef("T_A10", 1, "patch",
            inst -> inst.temperature.t_a10_patch,
            (inst, d) -> inst.temperature.t_a10_patch .= d),

        # --- Water state ---
        RestartVarDef("H2OSOI_LIQ", 2, "column",
            inst -> inst.water.waterstatebulk_inst.ws.h2osoi_liq_col,
            (inst, d) -> inst.water.waterstatebulk_inst.ws.h2osoi_liq_col .= d),
        RestartVarDef("H2OSOI_ICE", 2, "column",
            inst -> inst.water.waterstatebulk_inst.ws.h2osoi_ice_col,
            (inst, d) -> inst.water.waterstatebulk_inst.ws.h2osoi_ice_col .= d),
        RestartVarDef("H2OSOI_VOL", 2, "column",
            inst -> inst.water.waterstatebulk_inst.ws.h2osoi_vol_col,
            (inst, d) -> inst.water.waterstatebulk_inst.ws.h2osoi_vol_col .= d),
        RestartVarDef("H2OSNO", 1, "column",
            inst -> inst.water.waterstatebulk_inst.ws.h2osno_no_layers_col,
            (inst, d) -> inst.water.waterstatebulk_inst.ws.h2osno_no_layers_col .= d),
        RestartVarDef("H2OSFC", 1, "column",
            inst -> inst.water.waterstatebulk_inst.ws.h2osfc_col,
            (inst, d) -> inst.water.waterstatebulk_inst.ws.h2osfc_col .= d),
        RestartVarDef("WA", 1, "column",
            inst -> inst.water.waterstatebulk_inst.ws.wa_col,
            (inst, d) -> inst.water.waterstatebulk_inst.ws.wa_col .= d),

        # --- Soil hydrology ---
        RestartVarDef("ZWT", 1, "column",
            inst -> inst.soilhydrology.zwt_col,
            (inst, d) -> inst.soilhydrology.zwt_col .= d),
        RestartVarDef("ZWT_PERCHED", 1, "column",
            inst -> inst.soilhydrology.zwt_perched_col,
            (inst, d) -> inst.soilhydrology.zwt_perched_col .= d),

        # --- Snow state ---
        RestartVarDef("SNL", 1, "column",
            inst -> Float64.(inst.column.snl),
            (inst, d) -> inst.column.snl .= Int.(d)),
        RestartVarDef("SNOW_DEPTH", 1, "column",
            inst -> inst.water.waterdiagnosticbulk_inst.snow_depth_col,
            (inst, d) -> inst.water.waterdiagnosticbulk_inst.snow_depth_col .= d),
        RestartVarDef("FRAC_SNO", 1, "column",
            inst -> inst.water.waterdiagnosticbulk_inst.frac_sno_col,
            (inst, d) -> inst.water.waterdiagnosticbulk_inst.frac_sno_col .= d),
        RestartVarDef("FRAC_SNO_EFF", 1, "column",
            inst -> inst.water.waterdiagnosticbulk_inst.frac_sno_eff_col,
            (inst, d) -> inst.water.waterdiagnosticbulk_inst.frac_sno_eff_col .= d),
        RestartVarDef("INT_SNOW", 1, "column",
            inst -> inst.water.waterstatebulk_inst.int_snow_col,
            (inst, d) -> inst.water.waterstatebulk_inst.int_snow_col .= d),

        # Snow layer geometry (full snow+soil arrays)
        RestartVarDef("COL_DZ", 2, "column",
            inst -> inst.column.dz,
            (inst, d) -> inst.column.dz .= d),
        RestartVarDef("COL_Z", 2, "column",
            inst -> inst.column.z,
            (inst, d) -> inst.column.z .= d),
        RestartVarDef("COL_ZI", 2, "column",
            inst -> inst.column.zi,
            (inst, d) -> inst.column.zi .= d),

        # --- Canopy water ---
        RestartVarDef("LIQCAN", 1, "patch",
            inst -> inst.water.waterstatebulk_inst.ws.liqcan_patch,
            (inst, d) -> inst.water.waterstatebulk_inst.ws.liqcan_patch .= d),
        RestartVarDef("SNOCAN", 1, "patch",
            inst -> inst.water.waterstatebulk_inst.ws.snocan_patch,
            (inst, d) -> inst.water.waterstatebulk_inst.ws.snocan_patch .= d),

        # --- Lake state ---
        RestartVarDef("LAKE_ICEFRAC", 2, "column",
            inst -> inst.lakestate.lake_icefrac_col,
            (inst, d) -> inst.lakestate.lake_icefrac_col .= d),
        RestartVarDef("SAVEDTKE1", 1, "column",
            inst -> inst.lakestate.savedtke1_col,
            (inst, d) -> inst.lakestate.savedtke1_col .= d),

        # --- Canopy state ---
        RestartVarDef("TLAI", 1, "patch",
            inst -> inst.canopystate.tlai_patch,
            (inst, d) -> inst.canopystate.tlai_patch .= d),
        RestartVarDef("ELAI", 1, "patch",
            inst -> inst.canopystate.elai_patch,
            (inst, d) -> inst.canopystate.elai_patch .= d),
        RestartVarDef("HTOP", 1, "patch",
            inst -> inst.canopystate.htop_patch,
            (inst, d) -> inst.canopystate.htop_patch .= d),
        RestartVarDef("FSUN24", 1, "patch",
            inst -> inst.canopystate.fsun24_patch,
            (inst, d) -> inst.canopystate.fsun24_patch .= d),
        RestartVarDef("FSUN240", 1, "patch",
            inst -> inst.canopystate.fsun240_patch,
            (inst, d) -> inst.canopystate.fsun240_patch .= d),
        RestartVarDef("ELAI240", 1, "patch",
            inst -> inst.canopystate.elai240_patch,
            (inst, d) -> inst.canopystate.elai240_patch .= d),
    ]
end

"""
    write_restart!(filepath, inst, bounds; time=nothing)

Write all state variables to a restart NetCDF file.
"""
function write_restart!(filepath::String, inst::CLMInstances, bounds::BoundsType;
                         time::Union{DateTime,Nothing} = nothing)
    nc = bounds.endc
    np = bounds.endp

    ds = NCDataset(filepath, "c")
    try
        # Define dimensions
        defDim(ds, "column", nc)
        defDim(ds, "patch", np)

        # Add 2D dimensions based on array sizes
        nlevsno_nlevgrnd = size(inst.temperature.t_soisno_col, 2)
        nlevlak = size(inst.temperature.t_lake_col, 2)
        nlev_h2o = size(inst.water.waterstatebulk_inst.ws.h2osoi_liq_col, 2)
        nlev_vol = size(inst.water.waterstatebulk_inst.ws.h2osoi_vol_col, 2)

        defDim(ds, "levsno_levgrnd", nlevsno_nlevgrnd)
        defDim(ds, "levlak", nlevlak)
        defDim(ds, "lev_h2o", nlev_h2o)
        defDim(ds, "lev_vol", nlev_vol)

        # Write time if provided
        if time !== nothing
            defDim(ds, "time", 1)
            tvar = defVar(ds, "time", Float64, ("time",);
                          attrib = Dict("units" => "days since 2000-01-01"))
            tvar[1] = Dates.value(time - DateTime(2000, 1, 1)) / (1000 * 86400)
        end

        # Write each restart variable
        for rv in _restart_registry()
            data = try
                rv.getter(inst)
            catch
                continue
            end
            data === nothing && continue
            length(data) == 0 && continue

            dim_name = rv.level == "column" ? "column" : "patch"

            if rv.dims == 1
                v = defVar(ds, rv.name, Float64, (dim_name,))
                clean = collect(Float64, data)
                for i in eachindex(clean)
                    isfinite(clean[i]) || (clean[i] = -9999.0)
                end
                v[:] = clean
            else
                # Determine appropriate level dimension
                nlev = size(data, 2)
                lev_dim = if rv.name in ("T_SOISNO",)
                    "levsno_levgrnd"
                elseif rv.name in ("T_LAKE", "LAKE_ICEFRAC")
                    "levlak"
                elseif rv.name in ("H2OSOI_VOL",)
                    "lev_vol"
                else
                    "lev_h2o"
                end
                # Ensure level dim exists with right size
                if !haskey(ds.dim, lev_dim) || ds.dim[lev_dim] != nlev
                    lev_dim = "lev_$(rv.name)"
                    defDim(ds, lev_dim, nlev)
                end
                v = defVar(ds, rv.name, Float64, (dim_name, lev_dim))
                clean = collect(Float64, data)
                for i in eachindex(clean)
                    isfinite(clean[i]) || (clean[i] = -9999.0)
                end
                v[:, :] = clean
            end
        end
    finally
        close(ds)
    end

    return nothing
end

"""
    read_restart!(filepath, inst, bounds)

Read state variables from a restart NetCDF file and overwrite
the corresponding fields in `inst`.
"""
function read_restart!(filepath::String, inst::CLMInstances, bounds::BoundsType)
    isfile(filepath) || error("Restart file not found: $filepath")

    NCDataset(filepath, "r") do ds
        for rv in _restart_registry()
            haskey(ds, rv.name) || continue

            data = Array(ds[rv.name])
            data = replace(data, missing => NaN)

            # Replace fill values with NaN
            for i in eachindex(data)
                if data[i] == -9999.0
                    data[i] = NaN
                end
            end

            try
                rv.setter!(inst, Float64.(data))
            catch e
                @warn "Skipping restart variable $(rv.name): $e" maxlog=1
            end
        end
    end

    return nothing
end
