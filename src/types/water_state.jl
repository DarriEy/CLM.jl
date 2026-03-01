# ==========================================================================
# Ported from: src/biogeophys/WaterStateType.F90
# Water state variables that apply to both bulk water and water tracers.
# ==========================================================================

"""
    WaterStateData

Water state data structure. Holds water state variables at the column, patch,
landunit, and gridcell levels, including soil liquid/ice water content,
snow water, canopy water, aquifer water, and excess ice.

Ported from `waterstate_type` in `WaterStateType.F90`.
"""
Base.@kwdef mutable struct WaterStateData
    # --- Column-level 1D fields ---
    h2osno_no_layers_col   ::Vector{Float64} = Float64[]   # col snow not resolved into layers (mm H2O)
    h2osfc_col             ::Vector{Float64} = Float64[]   # col surface water (mm H2O)
    wa_col                 ::Vector{Float64} = Float64[]   # col water in unconfined aquifer (mm)
    dynbal_baseline_liq_col::Vector{Float64} = Float64[]   # col baseline liquid water for dynbal (mm H2O)
    dynbal_baseline_ice_col::Vector{Float64} = Float64[]   # col baseline ice for dynbal (mm H2O)
    exice_bulk_init        ::Vector{Float64} = Float64[]   # col initial excess ice concentration (unitless)

    # --- Column-level 2D fields ---
    h2osoi_liq_col         ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col liquid water (kg/m2) (-nlevsno+1:nlevmaxurbgrnd)
    h2osoi_ice_col         ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col ice lens (kg/m2) (-nlevsno+1:nlevmaxurbgrnd)
    h2osoi_vol_col         ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col volumetric soil water [m3/m3] (1:nlevmaxurbgrnd)
    excess_ice_col         ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col excess ice (kg/m2) (-nlevsno+1:nlevmaxurbgrnd)

    # --- Gridcell-level 2D fields ---
    h2osoi_vol_prs_grc     ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # grc prescribed volumetric soil water [m3/m3] (1:nlevgrnd)

    # --- Patch-level 1D fields ---
    snocan_patch           ::Vector{Float64} = Float64[]   # patch canopy snow water (mm H2O)
    liqcan_patch           ::Vector{Float64} = Float64[]   # patch canopy liquid water (mm H2O)

    # --- Landunit-level 1D fields ---
    stream_water_volume_lun::Vector{Float64} = Float64[]   # landunit volume of water in streams (m3)

    # --- Scalar ---
    aquifer_water_baseline ::Float64 = 0.0                 # baseline aquifer water for this bulk/tracer (mm)
end

"""
    waterstate_init!(ws::WaterStateData, nc::Int, np::Int, nl::Int, ng::Int)

Allocate and initialize all fields of a `WaterStateData` instance for
`nc` columns, `np` patches, `nl` landunits, and `ng` gridcells.
Real fields are initialized to `NaN`.

Ported from `waterstate_type%InitAllocate` in `WaterStateType.F90`.
"""
function waterstate_init!(ws::WaterStateData, nc::Int, np::Int, nl::Int, ng::Int)
    nlevgrnd       = varpar.nlevgrnd
    nlevsno        = varpar.nlevsno
    nlevmaxurbgrnd = varpar.nlevmaxurbgrnd
    nlevtot        = nlevsno + nlevmaxurbgrnd  # total snow+soil levels

    # --- Column 1D ---
    ws.h2osno_no_layers_col    = fill(NaN, nc)
    ws.h2osfc_col              = fill(NaN, nc)
    ws.wa_col                  = fill(NaN, nc)
    ws.dynbal_baseline_liq_col = fill(NaN, nc)
    ws.dynbal_baseline_ice_col = fill(NaN, nc)
    ws.exice_bulk_init         = fill(NaN, nc)

    # --- Column 2D ---
    ws.h2osoi_liq_col  = fill(NaN, nc, nlevtot)   # (-nlevsno+1:nlevmaxurbgrnd)
    ws.h2osoi_ice_col  = fill(NaN, nc, nlevtot)   # (-nlevsno+1:nlevmaxurbgrnd)
    ws.h2osoi_vol_col  = fill(NaN, nc, nlevmaxurbgrnd)  # (1:nlevmaxurbgrnd)
    ws.excess_ice_col  = fill(NaN, nc, nlevtot)    # (-nlevsno+1:nlevmaxurbgrnd)

    # --- Gridcell 2D ---
    ws.h2osoi_vol_prs_grc = fill(NaN, ng, nlevgrnd)

    # --- Patch 1D ---
    ws.snocan_patch = fill(NaN, np)
    ws.liqcan_patch = fill(NaN, np)

    # --- Landunit 1D ---
    ws.stream_water_volume_lun = fill(NaN, nl)

    return nothing
end

"""
    waterstate_clean!(ws::WaterStateData)

Deallocate (reset to empty) all fields of a `WaterStateData` instance.
"""
function waterstate_clean!(ws::WaterStateData)
    ws.h2osno_no_layers_col    = Float64[]
    ws.h2osfc_col              = Float64[]
    ws.wa_col                  = Float64[]
    ws.dynbal_baseline_liq_col = Float64[]
    ws.dynbal_baseline_ice_col = Float64[]
    ws.exice_bulk_init         = Float64[]

    ws.h2osoi_liq_col     = Matrix{Float64}(undef, 0, 0)
    ws.h2osoi_ice_col     = Matrix{Float64}(undef, 0, 0)
    ws.h2osoi_vol_col     = Matrix{Float64}(undef, 0, 0)
    ws.excess_ice_col     = Matrix{Float64}(undef, 0, 0)
    ws.h2osoi_vol_prs_grc = Matrix{Float64}(undef, 0, 0)

    ws.snocan_patch            = Float64[]
    ws.liqcan_patch            = Float64[]
    ws.stream_water_volume_lun = Float64[]

    return nothing
end

# ---------------------------------------------------------------------------
# Helper: convert Fortran snow-layer index (-nlevsno+1:nlevmaxurbgrnd) to
# Julia 1-based index (1:nlevtot).
# Fortran j ∈ [-nlevsno+1, nlevmaxurbgrnd] → Julia j + nlevsno
# ---------------------------------------------------------------------------
@inline _snow_idx(j::Int) = j + varpar.nlevsno

"""
    waterstate_init_cold!(ws::WaterStateData, bounds_col::UnitRange{Int},
                          bounds_patch::UnitRange{Int}, bounds_lun::UnitRange{Int},
                          bounds_grc::UnitRange{Int};
                          h2osno_input_col, watsat_col, t_soisno_col,
                          use_aquifer_layer, ratio,
                          snl_col, dz_col, landunit_col, lakpoi, urbpoi,
                          lun_itype, col_itype, nbedrock_col, gridcell_col,
                          exice_coldstart_depth, exice_init_conc_col)

Initialize cold-start conditions for water state variables.

Ported from `waterstate_type%InitCold` in `WaterStateType.F90`.
"""
function waterstate_init_cold!(ws::WaterStateData,
                                bounds_col::UnitRange{Int},
                                bounds_patch::UnitRange{Int},
                                bounds_lun::UnitRange{Int},
                                bounds_grc::UnitRange{Int};
                                h2osno_input_col::Vector{Float64},
                                watsat_col::Matrix{Float64},
                                t_soisno_col::Matrix{Float64},
                                use_aquifer_layer::Bool,
                                ratio::Float64 = 1.0,
                                snl_col::Vector{Int},
                                dz_col::Matrix{Float64},
                                landunit_col::Vector{Int},
                                lakpoi::BitVector,
                                urbpoi::BitVector,
                                lun_itype::Vector{Int},
                                col_itype::Vector{Int},
                                nbedrock_col::Vector{Int},
                                gridcell_col::Vector{Int},
                                exice_coldstart_depth::Float64 = 0.5,
                                exice_init_conc_col::Vector{Float64})
    nlevgrnd       = varpar.nlevgrnd
    nlevsoi        = varpar.nlevsoi
    nlevurb        = varpar.nlevurb
    nlevmaxurbgrnd = varpar.nlevmaxurbgrnd
    nlevsno        = varpar.nlevsno

    # Zero-init surface / canopy / stream water
    for c in bounds_col
        ws.h2osfc_col[c] = 0.0
    end
    for p in bounds_patch
        ws.snocan_patch[p] = 0.0
        ws.liqcan_patch[p] = 0.0
    end
    for l in bounds_lun
        ws.stream_water_volume_lun[l] = 0.0
    end

    # Initialize vol/liq/ice to spval
    for c in bounds_col
        for j in 1:nlevmaxurbgrnd
            ws.h2osoi_vol_col[c, j] = SPVAL
        end
    end
    for g in bounds_grc
        for j in 1:nlevgrnd
            ws.h2osoi_vol_prs_grc[g, j] = SPVAL
        end
    end
    nlevtot = nlevsno + nlevmaxurbgrnd
    for c in bounds_col
        for jj in 1:nlevtot
            ws.h2osoi_liq_col[c, jj] = SPVAL
            ws.h2osoi_ice_col[c, jj] = SPVAL
        end
    end

    # -------------------------------------------
    # Set soil water (non-lake columns)
    # -------------------------------------------
    for c in bounds_col
        l = landunit_col[c]
        if !lakpoi[l]  # not lake
            # volumetric water
            nlevs = 0
            if lun_itype[l] == ISTSOIL || lun_itype[l] == ISTCROP
                nlevs = nlevgrnd
                for j in 1:nlevs
                    if varctl.use_bedrock && nbedrock_col[c] <= nlevsoi
                        nbr = nbedrock_col[c]
                    else
                        nbr = nlevsoi
                    end
                    if j > nbr
                        ws.h2osoi_vol_col[c, j] = 0.0
                    else
                        if varctl.use_fates
                            ws.h2osoi_vol_col[c, j] = 0.75 * watsat_col[c, j] * ratio
                        else
                            ws.h2osoi_vol_col[c, j] = 0.15 * ratio
                        end
                    end
                end
            elseif urbpoi[l]
                if col_itype[c] == ICOL_ROAD_PERV
                    nlevs = nlevgrnd
                    for j in 1:nlevs
                        if j <= nlevsoi
                            ws.h2osoi_vol_col[c, j] = 0.3 * ratio
                        else
                            ws.h2osoi_vol_col[c, j] = 0.0
                        end
                    end
                elseif col_itype[c] == ICOL_ROAD_IMPERV
                    nlevs = nlevgrnd
                    for j in 1:nlevs
                        ws.h2osoi_vol_col[c, j] = 0.0
                    end
                else
                    nlevs = nlevurb
                    for j in 1:nlevs
                        ws.h2osoi_vol_col[c, j] = 0.0
                    end
                end
            elseif lun_itype[l] == ISTWET
                nlevs = nlevgrnd
                for j in 1:nlevs
                    if j > nlevsoi
                        ws.h2osoi_vol_col[c, j] = 0.0
                    else
                        ws.h2osoi_vol_col[c, j] = 1.0 * ratio
                    end
                end
            elseif lun_itype[l] == ISTICE
                nlevs = nlevgrnd
                for j in 1:nlevs
                    ws.h2osoi_vol_col[c, j] = 1.0 * ratio
                end
            else
                error("waterstate_init_cold!: unhandled landunit type $(lun_itype[l])")
            end

            # Clamp vol, compute liq/ice from vol
            for j in 1:nlevs
                ws.h2osoi_vol_col[c, j] = min(ws.h2osoi_vol_col[c, j], watsat_col[c, j] * ratio)
                jj = _snow_idx(j)  # convert soil index to combined array index
                if t_soisno_col[c, jj] <= TFRZ
                    ws.h2osoi_ice_col[c, jj] = dz_col[c, jj] * DENICE * ws.h2osoi_vol_col[c, j]
                    ws.h2osoi_liq_col[c, jj] = 0.0
                else
                    ws.h2osoi_ice_col[c, jj] = 0.0
                    ws.h2osoi_liq_col[c, jj] = dz_col[c, jj] * DENH2O * ws.h2osoi_vol_col[c, j]
                end
            end

            # Snow (non-layer)
            if snl_col[c] == 0
                ws.h2osno_no_layers_col[c] = h2osno_input_col[c] * ratio
            else
                ws.h2osno_no_layers_col[c] = 0.0
            end

            # Snow layers
            for j in (-nlevsno + 1):0
                jj = _snow_idx(j)
                if j > snl_col[c]
                    ws.h2osoi_ice_col[c, jj] = dz_col[c, jj] * 250.0 * ratio
                    ws.h2osoi_liq_col[c, jj] = 0.0
                end
            end
        end
    end

    # -------------------------------------------
    # Set lake water
    # -------------------------------------------
    for c in bounds_col
        l = landunit_col[c]
        if lakpoi[l]
            if snl_col[c] == 0
                ws.h2osno_no_layers_col[c] = h2osno_input_col[c] * ratio
            else
                ws.h2osno_no_layers_col[c] = 0.0
            end
            for j in (-nlevsno + 1):0
                jj = _snow_idx(j)
                if j > snl_col[c]
                    ws.h2osoi_ice_col[c, jj] = dz_col[c, jj] * BDSNO * ratio
                    ws.h2osoi_liq_col[c, jj] = 0.0
                end
            end
            for j in 1:nlevgrnd
                jj = _snow_idx(j)
                if j <= nlevsoi
                    ws.h2osoi_vol_col[c, j] = watsat_col[c, j] * ratio
                    ws.h2osoi_liq_col[c, jj] = SPVAL
                    ws.h2osoi_ice_col[c, jj] = SPVAL
                else
                    ws.h2osoi_vol_col[c, j] = 0.0
                end
            end
        end
    end

    # -------------------------------------------
    # For frozen layers (overwrites liq/ice based on temperature)
    # -------------------------------------------
    for c in bounds_col
        for j in 1:nlevmaxurbgrnd
            if ws.h2osoi_vol_col[c, j] != SPVAL
                jj = _snow_idx(j)
                if t_soisno_col[c, jj] <= TFRZ
                    ws.h2osoi_ice_col[c, jj] = dz_col[c, jj] * DENICE * ws.h2osoi_vol_col[c, j]
                    ws.h2osoi_liq_col[c, jj] = 0.0
                else
                    ws.h2osoi_ice_col[c, jj] = 0.0
                    ws.h2osoi_liq_col[c, jj] = dz_col[c, jj] * DENH2O * ws.h2osoi_vol_col[c, j]
                end
            end
        end
    end

    # -------------------------------------------
    # Aquifer water
    # -------------------------------------------
    ws.aquifer_water_baseline = AQUIFER_WATER_BASELINE * ratio
    for c in bounds_col
        ws.wa_col[c] = ws.aquifer_water_baseline
    end
    if use_aquifer_layer
        for c in bounds_col
            l = landunit_col[c]
            if !lakpoi[l]
                if urbpoi[l]
                    if col_itype[c] == ICOL_ROAD_PERV
                        ws.wa_col[c] = 4800.0 * ratio
                    else
                        ws.wa_col[c] = SPVAL
                    end
                else
                    ws.wa_col[c] = 4000.0 * ratio
                end
            end
        end
    end

    # -------------------------------------------
    # Dynamic balance baselines
    # -------------------------------------------
    for c in bounds_col
        ws.dynbal_baseline_liq_col[c] = 0.0
        ws.dynbal_baseline_ice_col[c] = 0.0
    end

    # -------------------------------------------
    # Excess ice initialization
    # -------------------------------------------
    for c in bounds_col
        ws.exice_bulk_init[c] = exice_init_conc_col[c]
    end
    nlevtot = nlevsno + nlevmaxurbgrnd
    for c in bounds_col
        for jj in 1:nlevtot
            ws.excess_ice_col[c, jj] = 0.0
        end
    end
    if varctl.use_excess_ice
        for c in bounds_col
            l = landunit_col[c]
            if !lakpoi[l]
                if lun_itype[l] == ISTSOIL || lun_itype[l] == ISTCROP
                    # Find layer containing exice_coldstart_depth
                    if zisoi[][nlevsoi] >= exice_coldstart_depth
                        nexice = _find_soil_layer_containing_depth(exice_coldstart_depth, nlevsoi)
                    else
                        nexice = nlevsoi - 1
                    end
                    if varctl.use_bedrock && nbedrock_col[c] <= nlevsoi
                        nbr = nbedrock_col[c]
                    else
                        nbr = nlevsoi
                    end
                    for j in 2:nlevmaxurbgrnd  # ignore first layer
                        jj = _snow_idx(j)
                        if nexice < nbr  # bedrock below coldstart depth
                            if j >= nexice && j < nbr && t_soisno_col[c, jj] <= TFRZ
                                ws.excess_ice_col[c, jj] = dz_col[c, jj] * DENICE * ws.exice_bulk_init[c]
                            else
                                ws.excess_ice_col[c, jj] = 0.0
                            end
                        else
                            ws.excess_ice_col[c, jj] = 0.0
                        end
                    end
                end
            else
                # Lakes and other columns: zero excess ice
                for jj in 1:nlevtot
                    ws.excess_ice_col[c, jj] = 0.0
                end
            end
        end
    end

    return nothing
end

"""
    _find_soil_layer_containing_depth(depth, nlevsoi) -> Int

Find the soil layer index containing the given depth.
Equivalent to Fortran `find_soil_layer_containing_depth`.
"""
function _find_soil_layer_containing_depth(depth::Float64, nlevsoi::Int)
    zi = zisoi[]
    for j in 1:nlevsoi
        if zi[j+1] > depth
            return j
        end
    end
    return nlevsoi
end

"""
    waterstate_calculate_total_h2osno!(ws::WaterStateData,
                                       mask::BitVector,
                                       bounds_col::UnitRange{Int},
                                       snl_col::Vector{Int},
                                       h2osno_total::Vector{Float64})

Calculate total snow water (h2osno_total) over the given column mask.
Sums h2osno_no_layers_col plus liquid and ice in snow layers.

Ported from `waterstate_type%CalculateTotalH2osno` in `WaterStateType.F90`.
"""
function waterstate_calculate_total_h2osno!(ws::WaterStateData,
                                             mask::BitVector,
                                             bounds_col::UnitRange{Int},
                                             snl_col::Vector{Int},
                                             h2osno_total::Vector{Float64})
    nlevsno = varpar.nlevsno
    for c in bounds_col
        mask[c] || continue
        h2osno_total[c] = ws.h2osno_no_layers_col[c]
        for j in (snl_col[c] + 1):0
            jj = _snow_idx(j)
            h2osno_total[c] += ws.h2osoi_ice_col[c, jj] + ws.h2osoi_liq_col[c, jj]
        end
    end
    return nothing
end

"""
    waterstate_check_snow_consistency!(ws::WaterStateData,
                                       mask::BitVector,
                                       bounds_col::UnitRange{Int},
                                       snl_col::Vector{Int},
                                       caller::String)

Check that unresolved snow exists only where snl==0 and that resolved snow
layers have no water outside the active snow pack.

Ported from `waterstate_type%CheckSnowConsistency` in `WaterStateType.F90`.
"""
function waterstate_check_snow_consistency!(ws::WaterStateData,
                                             mask::BitVector,
                                             bounds_col::UnitRange{Int},
                                             snl_col::Vector{Int},
                                             caller::String)
    nlevsno = varpar.nlevsno
    for c in bounds_col
        mask[c] || continue
        # Columns with snow layers should have zero h2osno_no_layers
        if snl_col[c] < 0
            if ws.h2osno_no_layers_col[c] != 0.0
                error("CheckSnowConsistency (called from $caller): " *
                      "col $c has snow layers (snl=$(snl_col[c])) but non-zero " *
                      "h2osno_no_layers=$(ws.h2osno_no_layers_col[c])")
            end
        end
        # No water outside resolved snow layers
        for j in (-nlevsno + 1):snl_col[c]
            jj = _snow_idx(j)
            ice_bad = (ws.h2osoi_ice_col[c, jj] != 0.0 && ws.h2osoi_ice_col[c, jj] != SPVAL)
            liq_bad = (ws.h2osoi_liq_col[c, jj] != 0.0 && ws.h2osoi_liq_col[c, jj] != SPVAL)
            if ice_bad || liq_bad
                error("CheckSnowConsistency (called from $caller): " *
                      "col $c has non-zero h2osoi_ice/liq outside resolved snow layers " *
                      "(j=$j, snl=$(snl_col[c]), ice=$(ws.h2osoi_ice_col[c, jj]), " *
                      "liq=$(ws.h2osoi_liq_col[c, jj]))")
            end
        end
    end
    return nothing
end

# ==========================================================================
# Stubs for infrastructure-dependent subroutines
# ==========================================================================

"""
    waterstate_init_history!(ws, bounds_col; use_aquifer_layer=true)

Register water state fields for history file output.

Ported from `waterstate_type%InitHistory` in `WaterStateType.F90`.
Requires history infrastructure — stub until that module is ported.
"""
function waterstate_init_history!(ws::WaterStateData,
                                   bounds_col::UnitRange{Int};
                                   use_aquifer_layer::Bool = true)
    return nothing
end

"""
    waterstate_restart!(ws, bounds_col; flag="read")

Read/write water state from/to restart file.

Ported from `waterstate_type%Restart` in `WaterStateType.F90`.
Requires NetCDF/restart infrastructure — stub until that module is ported.
"""
function waterstate_restart!(ws::WaterStateData,
                               bounds_col::UnitRange{Int};
                               flag::String = "read")
    return nothing
end
