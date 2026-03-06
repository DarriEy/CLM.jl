# ==========================================================================
# Minimal history output writer — NetCDF time series
# No direct Fortran equivalent (Fortran histFileMod is very complex)
#
# Aggregates patch- and column-level fields to gridcell level using subgrid
# area weights (patch%wtgcell, col%wtgcell), matching Fortran CLM h0 output.
#
# Public functions:
#   history_writer_init!   — Create output file with dimensions
#   history_write_step!    — Write one timestep of output fields
#   history_writer_close!  — Close output file
#   default_hist_fields    — Return list of default output fields
# ==========================================================================

"""
    HistFieldDef

Definition of a single history output field.
"""
struct HistFieldDef
    name::String
    long_name::String
    units::String
    level::String  # "gridcell", "column", or "patch"
    getter::Function  # (inst::CLMInstances) -> AbstractArray
end

"""
    HistoryWriter

Stateful NetCDF history output writer. All fields are aggregated to gridcell
level using subgrid area weights before writing.
"""
Base.@kwdef mutable struct HistoryWriter
    ds::Union{NCDataset, Nothing} = nothing
    fields::Vector{HistFieldDef} = HistFieldDef[]
    time_index::Int = 0
    filepath::String = ""
end

"""
    default_hist_fields() -> Vector{HistFieldDef}

Return default set of history output fields.
"""
function default_hist_fields()
    return HistFieldDef[
        HistFieldDef("T_GRND", "ground temperature", "K", "column",
            inst -> inst.temperature.t_grnd_col),
        HistFieldDef("QFLX_EVAP_TOT", "total evaporation", "mm/s", "column",
            inst -> inst.water.waterfluxbulk_inst.wf.qflx_evap_tot_col),
        HistFieldDef("EFLX_LH_TOT", "total latent heat flux", "W/m2", "patch",
            inst -> inst.energyflux.eflx_lh_tot_patch),
        HistFieldDef("EFLX_SH_TOT", "total sensible heat flux", "W/m2", "patch",
            inst -> inst.energyflux.eflx_sh_tot_patch),
        HistFieldDef("FSA", "absorbed solar radiation", "W/m2", "patch",
            inst -> inst.solarabs.fsa_patch),
        HistFieldDef("H2OSNO", "snow water equivalent", "mm", "column",
            history_h2osno_total_col),
        HistFieldDef("TSA", "2m air temperature", "K", "patch",
            inst -> inst.temperature.t_ref2m_patch),
        HistFieldDef("RAIN", "atmospheric rain", "mm/s", "column",
            inst -> length(inst.atm2lnd.forc_rain_downscaled_col) > 0 ?
                inst.atm2lnd.forc_rain_downscaled_col :
                zeros(length(inst.column.gridcell))),
        HistFieldDef("SNOW", "atmospheric snow", "mm/s", "column",
            inst -> length(inst.atm2lnd.forc_snow_downscaled_col) > 0 ?
                inst.atm2lnd.forc_snow_downscaled_col :
                zeros(length(inst.column.gridcell))),
        HistFieldDef("QRUNOFF", "total runoff", "mm/s", "column",
            inst -> inst.water.waterfluxbulk_inst.wf.qflx_runoff_col),
        HistFieldDef("ZWT", "water table depth", "m", "column",
            history_zwt_veg_col),
        HistFieldDef("TLAI", "total projected leaf area index", "m2/m2", "patch",
            inst -> inst.canopystate.tlai_patch),
        HistFieldDef("ELAI", "exposed one-sided leaf area index", "m2/m2", "patch",
            inst -> inst.canopystate.elai_patch),
        HistFieldDef("QFLX_DRAIN", "sub-surface drainage", "mm/s", "column",
            inst -> inst.water.waterfluxbulk_inst.wf.qflx_drain_col),
        HistFieldDef("QFLX_INFL", "infiltration", "mm/s", "column",
            inst -> inst.water.waterfluxbulk_inst.wf.qflx_infl_col),
        HistFieldDef("QFLX_TRAN_VEG", "vegetation transpiration", "mm/s", "column",
            inst -> inst.water.waterfluxbulk_inst.wf.qflx_tran_veg_col),
        HistFieldDef("SNOW_DEPTH", "snow depth", "m", "column",
            inst -> inst.water.waterdiagnosticbulk_inst.snow_depth_col),
        HistFieldDef("SNOWDP", "area-averaged snow depth", "m", "column",
            inst -> inst.water.waterdiagnosticbulk_inst.snowdp_col),
        HistFieldDef("FRAC_SNO", "fraction of ground covered by snow", "unitless", "column",
            inst -> inst.water.waterdiagnosticbulk_inst.frac_sno_col),
        HistFieldDef("BTRAN", "transpiration beta factor", "unitless", "patch",
            history_btran_daily_min_patch),
    ]
end

"""
    history_h2osno_total_col(inst) -> Vector{Float64}

Return total snow water equivalent (mm) per column, including both unresolved
snow (`h2osno_no_layers_col`) and explicit snow layers.
"""
function history_h2osno_total_col(inst::CLMInstances)
    col = inst.column
    ws = inst.water.waterstatebulk_inst.ws
    nc = length(col.snl)
    h2osno_total = zeros(nc)
    waterstate_calculate_total_h2osno!(ws, trues(nc), 1:nc, col.snl, h2osno_total)
    return h2osno_total
end

"""
    history_btran_daily_min_patch(inst) -> Vector{Float64}

Return daily minimum BTRAN when available to align with Fortran `BTRANMN`
diagnostics.

Preference order:
1. `btran_min_inst_patch` when finite (running minimum for current day)
2. `btran_min_patch` when finite (completed-day minimum)
3. fallback to instantaneous `btran_patch`
"""
function history_btran_daily_min_patch(inst::CLMInstances)
    ef = inst.energyflux
    pch = inst.patch
    np = length(ef.btran_patch)

    out = fill(NaN, np)
    for p in 1:np
        p <= length(pch.itype) || continue
        # Fortran BTRANMN uses a vegetation mask: exclude bare-ground patches.
        if pch.itype[p] != noveg
            bt = ef.btran_patch[p]
            if length(ef.btran_min_inst_patch) >= p && isfinite(ef.btran_min_inst_patch[p]) &&
               ef.btran_min_inst_patch[p] != SPVAL
                bt = ef.btran_min_inst_patch[p]
            elseif length(ef.btran_min_patch) >= p && isfinite(ef.btran_min_patch[p]) &&
                   ef.btran_min_patch[p] != SPVAL
                bt = ef.btran_min_patch[p]
            end
            out[p] = bt
        end
    end
    return out
end

"""
    history_zwt_veg_col(inst) -> Vector{Float64}

Return water table depth for vegetated/crop landunits only, matching Fortran
`ZWT` landunit mask semantics.
"""
function history_zwt_veg_col(inst::CLMInstances)
    sh = inst.soilhydrology
    col = inst.column
    nc = length(sh.zwt_col)
    out = fill(NaN, nc)

    # Fortran ZWT uses a vegetation mask at column level.
    # Mark columns with at least one vegetated patch (itype != noveg).
    has_veg_col = falses(nc)
    pch = inst.patch
    for p in eachindex(pch.column)
        p <= length(pch.itype) || continue
        c = pch.column[p]
        c >= 1 && c <= nc || continue
        if pch.itype[p] != noveg
            has_veg_col[c] = true
        end
    end

    for c in 1:nc
        c <= length(col.landunit) || continue
        if has_veg_col[c]
            out[c] = sh.zwt_col[c]
        end
    end
    return out
end

"""
    history_writer_init!(hw, filepath, fields, ng, nc, np)

Create the NetCDF output file with dimensions and variable definitions.
All fields are written on the gridcell ("lndgrid") dimension after area-weighted
aggregation, matching Fortran CLM h0 output format.
"""
function history_writer_init!(hw::HistoryWriter, filepath::String,
                              fields::Vector{HistFieldDef},
                              ng::Int, nc::Int, np::Int)
    hw.filepath = filepath
    hw.fields = fields
    hw.time_index = 0

    hw.ds = NCDataset(filepath, "c")
    ds = hw.ds

    # Define dimensions — all output is on lndgrid (gridcell) level
    defDim(ds, "lndgrid", ng)
    defDim(ds, "time", Inf)  # unlimited

    # Time variable
    defVar(ds, "time", Float64, ("time",);
           attrib = Dict("units" => "days since 2000-01-01",
                         "calendar" => "noleap"))

    # Define each field on the gridcell dimension
    for f in fields
        defVar(ds, f.name, Float64, ("lndgrid", "time");
               attrib = Dict("long_name" => f.long_name,
                             "units" => f.units,
                             "subgrid_level" => f.level))
    end

    return nothing
end

"""
    history_write_step!(hw, inst, current_time)

Write one timestep of all registered fields to the output file.
Patch- and column-level fields are aggregated to gridcell level using
area weights (patch%wtgcell, col%wtgcell) before writing.
"""
function history_write_step!(hw::HistoryWriter, inst::CLMInstances,
                             current_time::DateTime)
    ds = hw.ds
    ds === nothing && return nothing

    hw.time_index += 1
    ti = hw.time_index

    # Write time (days since 2000-01-01)
    ref_date = DateTime(2000, 1, 1)
    days = Dates.value(current_time - ref_date) / (1000 * 86400)  # ms → days
    ds["time"][ti] = days

    # Subgrid metadata for aggregation
    pch = inst.patch
    col = inst.column
    ng = length(col.gridcell) > 0 ? maximum(col.gridcell) : 0

    # Write each field, aggregated to gridcell level
    for f in hw.fields
        data = try
            f.getter(inst)
        catch
            nothing
        end
        data === nothing && continue
        length(data) == 0 && continue

        # Aggregate to gridcell level with renormalization
        # If some subgrid units have non-finite values (e.g. lake NaN),
        # renormalize by the sum of contributing weights.
        gcell_val = zeros(ng)
        gcell_wt  = zeros(ng)  # sum of contributing weights

        if f.level == "patch"
            # Area-weighted aggregation: sum(patch_val * patch%wtgcell)
            np = length(data)
            for p in 1:np
                val = Float64(data[p])
                if isfinite(val) && p <= length(pch.gridcell) && p <= length(pch.wtgcell)
                    g = pch.gridcell[p]
                    if g >= 1 && g <= ng
                        gcell_val[g] += val * pch.wtgcell[p]
                        gcell_wt[g]  += pch.wtgcell[p]
                    end
                end
            end
        elseif f.level == "column"
            # Area-weighted aggregation: sum(col_val * col%wtgcell)
            nc = length(data)
            for c in 1:nc
                val = Float64(data[c])
                if isfinite(val) && c <= length(col.gridcell) && c <= length(col.wtgcell)
                    g = col.gridcell[c]
                    if g >= 1 && g <= ng
                        gcell_val[g] += val * col.wtgcell[c]
                        gcell_wt[g]  += col.wtgcell[c]
                    end
                end
            end
        else
            # Already gridcell level
            for g in 1:min(ng, length(data))
                val = Float64(data[g])
                gcell_val[g] = isfinite(val) ? val : -9999.0
                gcell_wt[g]  = isfinite(val) ? 1.0 : 0.0
            end
        end

        # Renormalize by contributing weight sum (handles missing subgrid units)
        for g in 1:ng
            if gcell_wt[g] > 0.0
                gcell_val[g] /= gcell_wt[g]
            else
                gcell_val[g] = -9999.0
            end
        end

        ds[f.name][:, ti] = gcell_val
    end

    return nothing
end

"""
    history_writer_close!(hw)

Close the history NetCDF file.
"""
function history_writer_close!(hw::HistoryWriter)
    if hw.ds !== nothing
        close(hw.ds)
        hw.ds = nothing
    end
    return nothing
end
