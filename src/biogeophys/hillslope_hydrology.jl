# ==========================================================================
# Ported from: src/biogeophys/HillslopeHydrologyMod.F90 (~1149 lines)
# Hillslope hydrology: geomorphological parameters, stream outflow,
# stream water update, soil thickness profiles, PFT distributions.
#
# Public functions:
#   hillslope_properties_init!           — Initialize hillslope config (stub)
#   check_aquifer_layer!                 — Validate aquifer/hillslope compatibility
#   init_hillslope!                      — Initialize hillslope column properties (stub)
#   set_hillslope_soil_thickness!        — Set soil thickness for hillslope columns
#   hillslope_soil_thickness_profile!    — Modify nbedrock by profile method
#   hillslope_soil_thickness_profile_linear! — Linear soil thickness profile
#   hillslope_set_lowland_upland_pfts!   — Reassign PFTs by lowland/upland
#   hillslope_dominant_lowland_pft!      — Assign dominant PFTs (stub)
#   hillslope_pft_from_file!             — Assign PFTs from index array
#   hillslope_stream_outflow!            — Stream discharge via Manning equation
#   hillslope_update_stream_water!       — Update stream water volume
# ==========================================================================

# -------------------------------------------------------------------------
# Module-level constants
# -------------------------------------------------------------------------

# Streamflow methods
const STREAMFLOW_MANNING = 0

# PFT distribution methods
const PFT_STANDARD              = 0
const PFT_FROM_FILE             = 1
const PFT_UNIFORM_DOMINANT_PFT  = 2
const PFT_LOWLAND_DOMINANT_PFT  = 3
const PFT_LOWLAND_UPLAND        = 4

# Soil profile methods (private in Fortran)
const SOIL_PROFILE_UNIFORM              = 0
const SOIL_PROFILE_FROM_FILE            = 1
const SOIL_PROFILE_SET_LOWLAND_UPLAND   = 2
const SOIL_PROFILE_LINEAR               = 3

# Mutable module-level state (mirrors Fortran module variables)
const hillslope_config = Ref((
    pft_distribution_method = PFT_STANDARD,
    soil_profile_method     = SOIL_PROFILE_UNIFORM,
))

# =========================================================================
# hillslope_properties_init!
# =========================================================================

"""
    hillslope_properties_init!(;
        pft_distribution_method::String="Standard",
        soil_profile_method::String="Uniform")

Initialize hillslope hydrology configuration by setting PFT distribution
and soil profile methods from string arguments.

Stub: In Fortran this reads a namelist and broadcasts via MPI. Here we
accept the method strings directly.

Ported from `hillslope_properties_init` in `HillslopeHydrologyMod.F90`.
"""
function hillslope_properties_init!(;
    pft_distribution_method::String = "Standard",
    soil_profile_method::String = "Uniform"
)
    pft_method = if pft_distribution_method == "Standard"
        PFT_STANDARD
    elseif pft_distribution_method == "FromFile"
        PFT_FROM_FILE
    elseif pft_distribution_method == "DominantPftUniform"
        PFT_UNIFORM_DOMINANT_PFT
    elseif pft_distribution_method == "DominantPftLowland"
        PFT_LOWLAND_DOMINANT_PFT
    elseif pft_distribution_method == "PftLowlandUpland"
        PFT_LOWLAND_UPLAND
    else
        error("ERROR: bad value for pft_distribution_method: $pft_distribution_method")
    end

    soil_method = if soil_profile_method == "Uniform"
        SOIL_PROFILE_UNIFORM
    elseif soil_profile_method == "FromFile"
        SOIL_PROFILE_FROM_FILE
    elseif soil_profile_method == "SetLowlandUpland"
        SOIL_PROFILE_SET_LOWLAND_UPLAND
    elseif soil_profile_method == "Linear"
        SOIL_PROFILE_LINEAR
    else
        error("ERROR: bad value for soil_profile_method: $soil_profile_method")
    end

    hillslope_config[] = (
        pft_distribution_method = pft_method,
        soil_profile_method     = soil_method,
    )

    return nothing
end

# =========================================================================
# check_aquifer_layer!
# =========================================================================

"""
    check_aquifer_layer!(use_hillslope::Bool, use_aquifer_layer::Bool)

Validate that `use_hillslope` and `use_aquifer_layer` are not both enabled.

Ported from `check_aquifer_layer` in `HillslopeHydrologyMod.F90`.
"""
function check_aquifer_layer!(use_hillslope::Bool, use_aquifer_layer::Bool)
    if use_hillslope && use_aquifer_layer
        error("ERROR: use_hillslope and use_aquifer_layer cannot both be set to true")
    end
    return nothing
end

# =========================================================================
# init_hillslope!
# =========================================================================

"""
    init_hillslope!(col, lun, grc, bounds_l, bounds_c;
        nhillslope, max_columns_hillslope,
        ncolumns_hillslope, pct_hillslope,
        hill_ndx, col_ndx, col_dndx,
        hill_slope, hill_aspect, hill_area,
        hill_dist, hill_width, hill_elev,
        hill_pftndx=nothing,
        stream_channel_depth=nothing,
        stream_channel_width=nothing,
        stream_channel_slope=nothing,
        use_hillslope_routing=false)

Initialize hillslope geomorphology from pre-loaded arrays rather than
reading from a NetCDF file. Sets column-level hillslope properties
(distances, widths, elevations, slopes, areas, neighbor indices) and
optionally stream channel properties on landunits.

Stub: The Fortran version reads from NetCDF. Here the arrays are provided
directly. Weight recalculation and subgrid computations are omitted.

Ported from `InitHillslope` in `HillslopeHydrologyMod.F90`.
"""
function init_hillslope!(
    col::ColumnData,
    lun::LandunitData,
    grc::GridcellData,
    bounds_l::UnitRange{Int},
    bounds_c::UnitRange{Int};
    nhillslope::Int,
    max_columns_hillslope::Int,
    ncolumns_hillslope::Vector{Int},
    pct_hillslope::Matrix{Float64},
    hill_ndx::Matrix{Int},
    col_ndx_in::Matrix{Int},
    col_dndx::Matrix{Int},
    hill_slope_in::Matrix{Float64},
    hill_aspect_in::Matrix{Float64},
    hill_area_in::Matrix{Float64},
    hill_dist_in::Matrix{Float64},
    hill_width_in::Matrix{Float64},
    hill_elev_in::Matrix{Float64},
    hill_pftndx::Union{Matrix{Int}, Nothing} = nothing,
    stream_channel_depth::Union{Vector{Float64}, Nothing} = nothing,
    stream_channel_width::Union{Vector{Float64}, Nothing} = nothing,
    stream_channel_slope::Union{Vector{Float64}, Nothing} = nothing,
    use_hillslope_routing::Bool = false
)
    # Mark hillslope columns
    for l in bounds_l
        if lun.itype[l] == ISTSOIL && ncolumns_hillslope[l] > 0 && lun.wtgcell[l] > 0.0
            for c in lun.coli[l]:lun.colf[l]
                col.is_hillslope_column[c] = true
            end
        end
    end

    # Set column-level hillslope properties
    for l in bounds_l
        if lun.itype[l] == ISTSOIL

            # Map external to internal column index and set downhill neighbor
            for c in lun.coli[l]:lun.colf[l]
                ci = c - lun.coli[l] + 1
                if col_dndx[l, ci] <= -999
                    col.cold[c] = ISPVAL
                else
                    col.cold[c] = c + (col_dndx[l, ci] - col_ndx_in[l, ci])
                end
            end

            # Set remaining column properties
            for c in lun.coli[l]:lun.colf[l]
                ci = c - lun.coli[l] + 1

                col.hillslope_ndx[c] = hill_ndx[l, ci]

                # Find uphill neighbor
                col.colu[c] = ISPVAL
                for i in lun.coli[l]:lun.colf[l]
                    if c == col.cold[i]
                        col.colu[c] = i
                    end
                end

                col.hill_distance[c] = hill_dist_in[l, ci]
                col.hill_width[c]    = hill_width_in[l, ci]
                col.hill_elev[c]     = hill_elev_in[l, ci]
                col.hill_slope[c]    = hill_slope_in[l, ci]
                col.hill_area[c]     = hill_area_in[l, ci]
                col.hill_aspect[c]   = hill_aspect_in[l, ci]
            end

            # Calculate total hillslope area and columns per hillslope
            ncol_per_hillslope = zeros(nhillslope)
            hillslope_area     = zeros(nhillslope)
            for c in lun.coli[l]:lun.colf[l]
                nh = col.hillslope_ndx[c]
                if nh > 0
                    ncol_per_hillslope[nh] += 1.0
                    hillslope_area[nh] += col.hill_area[c]
                end
            end

            if use_hillslope_routing && stream_channel_depth !== nothing
                g = lun.gridcell[l]
                lun.stream_channel_depth[l]  = stream_channel_depth[l]
                lun.stream_channel_width[l]  = stream_channel_width[l]
                lun.stream_channel_slope[l]  = stream_channel_slope[l]

                # Number of channels
                nhill_per_landunit = zeros(nhillslope)
                lun.stream_channel_number[l] = 0.0
                for nh in 1:nhillslope
                    if hillslope_area[nh] > 0.0
                        nhill_per_landunit[nh] = grc.area[g] * 1.0e6 * lun.wtgcell[l] *
                            pct_hillslope[l, nh] * 0.01 / hillslope_area[nh]
                        lun.stream_channel_number[l] += 0.5 * nhill_per_landunit[nh]
                    end
                end

                # Channel length
                lun.stream_channel_length[l] = 0.0
                for c in lun.coli[l]:lun.colf[l]
                    if col.cold[c] == ISPVAL
                        lun.stream_channel_length[l] += col.hill_width[c] *
                            0.5 * nhill_per_landunit[col.hillslope_ndx[c]]
                    end
                end
            end

            # Recalculate column weights
            for c in lun.coli[l]:lun.colf[l]
                nh = col.hillslope_ndx[c]
                if col.is_hillslope_column[c] && hillslope_area[nh] > 0.0
                    col.wtlunit[c] = (col.hill_area[c] / hillslope_area[nh]) *
                        (pct_hillslope[l, nh] * 0.01)
                end
            end
        end
    end

    return nothing
end

# =========================================================================
# set_hillslope_soil_thickness!
# =========================================================================

"""
    set_hillslope_soil_thickness!(col, lun, bounds_c;
        soil_profile_method, nlevsoi,
        soil_depth_lowland=8.0, soil_depth_upland=8.0,
        bedrock_from_file=nothing)

Set hillslope column nbedrock values based on the selected soil profile method.

For `SOIL_PROFILE_FROM_FILE`, uses bedrock depths from the `bedrock_from_file` vector.
For `SOIL_PROFILE_SET_LOWLAND_UPLAND` or `SOIL_PROFILE_LINEAR`, calls
`hillslope_soil_thickness_profile!`.

Ported from `SetHillslopeSoilThickness` in `HillslopeHydrologyMod.F90`.
"""
function set_hillslope_soil_thickness!(
    col::ColumnData,
    lun::LandunitData,
    bounds_c::UnitRange{Int};
    soil_profile_method::Int,
    nlevsoi::Int,
    soil_depth_lowland::Float64 = 8.0,
    soil_depth_upland::Float64 = 8.0,
    bedrock_from_file::Union{Vector{Float64}, Nothing} = nothing
)
    zisoi_val = zisoi[]

    if soil_profile_method == SOIL_PROFILE_FROM_FILE
        if bedrock_from_file === nothing
            error("ERROR: soil_profile_method = FromFile, but bedrock_from_file not provided")
        end
        for c in bounds_c
            if col.is_hillslope_column[c] && col.active[c]
                depth = bedrock_from_file[c]
                for j in 1:nlevsoi
                    if zisoi_val[j] > ZMIN_BEDROCK
                        if zisoi_val[j] < depth && zisoi_val[j+1] >= depth
                            col.nbedrock[c] = j
                        end
                    end
                end
            end
        end

    elseif soil_profile_method == SOIL_PROFILE_SET_LOWLAND_UPLAND ||
           soil_profile_method == SOIL_PROFILE_LINEAR
        hillslope_soil_thickness_profile!(
            col, bounds_c,
            soil_profile_method = soil_profile_method,
            nlevsoi = nlevsoi,
            soil_depth_lowland = soil_depth_lowland,
            soil_depth_upland = soil_depth_upland)

    elseif soil_profile_method != SOIL_PROFILE_UNIFORM
        error("ERROR: unrecognized hillslope_soil_profile_method: $soil_profile_method")
    end

    return nothing
end

# =========================================================================
# hillslope_soil_thickness_profile!
# =========================================================================

"""
    hillslope_soil_thickness_profile!(col, bounds_c;
        soil_profile_method, nlevsoi,
        soil_depth_lowland=8.0, soil_depth_upland=8.0)

Modify soil thickness across hillslope columns by changing `col.nbedrock`.

For `SOIL_PROFILE_SET_LOWLAND_UPLAND`: lowland columns (cold == ISPVAL) get
`soil_depth_lowland`, upland columns get `soil_depth_upland`.

For `SOIL_PROFILE_LINEAR`: delegates to `hillslope_soil_thickness_profile_linear!`.

Ported from `HillslopeSoilThicknessProfile` in `HillslopeHydrologyMod.F90`.
"""
function hillslope_soil_thickness_profile!(
    col::ColumnData,
    bounds_c::UnitRange{Int};
    soil_profile_method::Int,
    nlevsoi::Int,
    soil_depth_lowland::Float64 = 8.0,
    soil_depth_upland::Float64 = 8.0
)
    zisoi_val = zisoi[]

    if soil_profile_method == SOIL_PROFILE_SET_LOWLAND_UPLAND
        for c in bounds_c
            if col.is_hillslope_column[c] && col.active[c]
                if col.cold[c] != ISPVAL
                    # Upland
                    for j in 1:nlevsoi
                        if zisoi_val[j] > ZMIN_BEDROCK
                            if zisoi_val[j] < soil_depth_upland && zisoi_val[j+1] >= soil_depth_upland
                                col.nbedrock[c] = j
                            end
                        end
                    end
                else
                    # Lowland
                    for j in 1:nlevsoi
                        if zisoi_val[j] > ZMIN_BEDROCK
                            if zisoi_val[j] < soil_depth_lowland && zisoi_val[j+1] >= soil_depth_lowland
                                col.nbedrock[c] = j
                            end
                        end
                    end
                end
            end
        end

    elseif soil_profile_method == SOIL_PROFILE_LINEAR
        hillslope_soil_thickness_profile_linear!(
            col, bounds_c,
            nlevsoi = nlevsoi,
            soil_depth_lowland = soil_depth_lowland,
            soil_depth_upland = soil_depth_upland)

    else
        error("ERROR: invalid soil_profile_method: $soil_profile_method")
    end

    return nothing
end

# =========================================================================
# hillslope_soil_thickness_profile_linear!
# =========================================================================

"""
    hillslope_soil_thickness_profile_linear!(col, bounds_c;
        nlevsoi, soil_depth_lowland=8.0, soil_depth_upland=8.0)

Apply a linear soil depth profile across hillslope columns within each
landunit. The depth varies linearly from `soil_depth_lowland` at the
lowermost column (largest hill_distance) to `soil_depth_upland` at the
uppermost column (smallest hill_distance).

Ported from `HillslopeSoilThicknessProfile_linear` in
`HillslopeHydrologyUtilsMod.F90`.
"""
function hillslope_soil_thickness_profile_linear!(
    col::ColumnData,
    bounds_c::UnitRange{Int};
    nlevsoi::Int,
    soil_depth_lowland::Float64 = 8.0,
    soil_depth_upland::Float64 = 8.0
)
    zisoi_val = zisoi[]

    # Find min/max hill_distance across active hillslope columns
    min_dist =  Inf
    max_dist = -Inf
    for c in bounds_c
        if col.is_hillslope_column[c] && col.active[c]
            d = col.hill_distance[c]
            if d < min_dist
                min_dist = d
            end
            if d > max_dist
                max_dist = d
            end
        end
    end

    if max_dist > min_dist
        # Linear interpolation slope and intercept
        m = (soil_depth_upland - soil_depth_lowland) / (max_dist - min_dist)
        b = soil_depth_lowland - m * min_dist
    else
        m = 0.0
        b = soil_depth_lowland
    end

    for c in bounds_c
        if col.is_hillslope_column[c] && col.active[c]
            soil_depth_col = m * col.hill_distance[c] + b
            for j in 1:nlevsoi
                if zisoi_val[j] > ZMIN_BEDROCK
                    if zisoi_val[j] < soil_depth_col && zisoi_val[j+1] >= soil_depth_col
                        col.nbedrock[c] = j
                    end
                end
            end
        end
    end

    return nothing
end

# =========================================================================
# hillslope_set_lowland_upland_pfts!
# =========================================================================

"""
    hillslope_set_lowland_upland_pfts!(col, pch, bounds_c;
        lowland_ivt, upland_ivt, natpft_lb=0)

Reassign patch type of each hillslope column based on whether it is lowland
(cold == ISPVAL) or upland. Assumes each column has a single PFT.

Ported from `HillslopeSetLowlandUplandPfts` in `HillslopeHydrologyMod.F90`.
"""
function hillslope_set_lowland_upland_pfts!(
    col::ColumnData,
    pch::PatchData,
    bounds_c::UnitRange{Int};
    lowland_ivt::Int,
    upland_ivt::Int,
    natpft_lb::Int = 0
)
    for c in bounds_c
        if col.is_hillslope_column[c]
            npatches_per_column = 0
            for p in col.patchi[c]:col.patchf[c]
                if col.cold[c] == ISPVAL
                    # Lowland
                    pch.itype[p] = lowland_ivt
                else
                    # Upland
                    pch.itype[p] = upland_ivt
                end
                pch.mxy[p] = pch.itype[p] + (1 - natpft_lb)
                npatches_per_column += 1
            end
            if npatches_per_column != 1
                error("ERROR: number of patches per hillslope column not equal to 1 (c=$c, n=$npatches_per_column)")
            end
        end
    end

    return nothing
end

# =========================================================================
# hillslope_dominant_lowland_pft!
# =========================================================================

"""
    hillslope_dominant_lowland_pft!(col, pch, bounds_c)

Reassign patch weights of each hillslope column based on each gridcell's
two most dominant PFTs. Places the highest stature vegetation on the lowland
column and the second on upland columns.

Stub: This function depends on `find_k_max_indices` and `pftcon` (PFT
properties like is_tree, is_grass, is_shrub) which are not yet ported.
When those modules are available, this should be implemented.

Ported from `HillslopeDominantLowlandPft` in `HillslopeHydrologyMod.F90`.
"""
function hillslope_dominant_lowland_pft!(
    col::ColumnData,
    pch::PatchData,
    bounds_c::UnitRange{Int}
)
    # Stub: depends on find_k_max_indices and pftcon (PFT properties)
    # which are not yet ported. No-op until those modules are available.
    return nothing
end

# =========================================================================
# hillslope_pft_from_file!
# =========================================================================

"""
    hillslope_pft_from_file!(col, pch, bounds_c, col_pftndx;
        natpft_lb=0)

Reassign patch type using indices from a data array. Assumes one patch per
hillslope column.

Ported from `HillslopePftFromFile` in `HillslopeHydrologyMod.F90`.
"""
function hillslope_pft_from_file!(
    col::ColumnData,
    pch::PatchData,
    bounds_c::UnitRange{Int},
    col_pftndx::Vector{Int};
    natpft_lb::Int = 0
)
    for c in bounds_c
        if col.is_hillslope_column[c]
            npatches_per_column = 0
            for p in col.patchi[c]:col.patchf[c]
                pch.itype[p] = col_pftndx[c]
                pch.mxy[p] = pch.itype[p] + (1 - natpft_lb)
                npatches_per_column += 1
            end
            if npatches_per_column != 1
                error("ERROR: number of patches per hillslope column not equal to 1 (c=$c, n=$npatches_per_column)")
            end
        end
    end

    return nothing
end

# =========================================================================
# hillslope_stream_outflow!
# =========================================================================

"""
    hillslope_stream_outflow!(
        stream_water_volume, volumetric_streamflow,
        lun, bounds_l, dtime;
        streamflow_method=STREAMFLOW_MANNING)

Calculate discharge from stream channel using Manning's equation.
Handles overbank flow when stream depth exceeds channel depth.

Ported from `HillslopeStreamOutflow` in `HillslopeHydrologyMod.F90`.
"""
function hillslope_stream_outflow!(
    stream_water_volume::Vector{Float64},
    volumetric_streamflow::Vector{Float64},
    lun::LandunitData,
    bounds_l::UnitRange{Int},
    dtime::Float64;
    streamflow_method::Int = STREAMFLOW_MANNING
)
    manning_roughness = 0.03
    manning_exponent  = 0.667
    overbank_method   = 1   # method 1: increase dynamic slope

    for l in bounds_l
        volumetric_streamflow[l] = 0.0

        # Check for vegetated landunits with initialized stream channel properties
        active_stream = false
        if lun.itype[l] == ISTSOIL &&
           lun.stream_channel_length[l] > 0.0 &&
           lun.stream_channel_width[l] > 0.0
            active_stream = true
        end

        if lun.active[l] && active_stream
            if streamflow_method == STREAMFLOW_MANNING
                cross_sectional_area = stream_water_volume[l] /
                    lun.stream_channel_length[l]
                stream_depth = cross_sectional_area /
                    lun.stream_channel_width[l]
                hydraulic_radius = cross_sectional_area /
                    (lun.stream_channel_width[l] + 2.0 * stream_depth)

                if hydraulic_radius <= 0.0
                    volumetric_streamflow[l] = 0.0
                else
                    flow_velocity = (hydraulic_radius)^manning_exponent *
                        sqrt(lun.stream_channel_slope[l]) /
                        manning_roughness

                    # Overbank flow
                    if stream_depth > lun.stream_channel_depth[l]
                        if overbank_method == 1
                            # Increase dynamic slope
                            volumetric_streamflow[l] = cross_sectional_area * flow_velocity *
                                (stream_depth / lun.stream_channel_depth[l])
                        elseif overbank_method == 2
                            # Increase flow area cross section
                            overbank_area = (stream_depth - lun.stream_channel_depth[l]) *
                                30.0 * lun.stream_channel_width[l]
                            volumetric_streamflow[l] = (cross_sectional_area + overbank_area) * flow_velocity
                        elseif overbank_method == 3
                            # Remove all overbank flow instantly
                            volumetric_streamflow[l] = cross_sectional_area * flow_velocity +
                                (stream_depth - lun.stream_channel_depth[l]) *
                                lun.stream_channel_width[l] * lun.stream_channel_length[l] / dtime
                        else
                            error("ERROR: invalid overbank_method: $overbank_method")
                        end
                    else
                        volumetric_streamflow[l] = cross_sectional_area * flow_velocity
                    end

                    # Scale by number of channel reaches
                    volumetric_streamflow[l] *= lun.stream_channel_number[l]

                    # Clamp to [0, available_water/dtime]
                    volumetric_streamflow[l] = max(0.0,
                        min(volumetric_streamflow[l], stream_water_volume[l] / dtime))
                end
            else
                error("ERROR: invalid streamflow_method: $streamflow_method")
            end
        end
    end

    return nothing
end

# =========================================================================
# hillslope_update_stream_water!
# =========================================================================

"""
    hillslope_update_stream_water!(
        stream_water_volume, volumetric_streamflow, stream_water_depth,
        qflx_drain, qflx_drain_perched, qflx_surf,
        col, lun, grc,
        bounds_l, dtime)

Update stream water volume by accumulating column-level fluxes (drainage,
perched drainage, surface runoff) and removing stream discharge. Adjusts
for negative drainage. Computes stream water depth.

Ported from `HillslopeUpdateStreamWater` in `HillslopeHydrologyMod.F90`.
"""
function hillslope_update_stream_water!(
    stream_water_volume::Vector{Float64},
    volumetric_streamflow::Vector{Float64},
    stream_water_depth::Vector{Float64},
    qflx_drain::Vector{Float64},
    qflx_drain_perched::Vector{Float64},
    qflx_surf::Vector{Float64},
    col::ColumnData,
    lun::LandunitData,
    grc::GridcellData,
    bounds_l::UnitRange{Int},
    dtime::Float64
)
    for l in bounds_l
        # Check for vegetated landunits with initialized stream channel properties
        active_stream = false
        if lun.itype[l] == ISTSOIL &&
           lun.stream_channel_length[l] > 0.0 &&
           lun.stream_channel_width[l] > 0.0
            active_stream = true
        end

        if lun.active[l] && active_stream
            g = lun.gridcell[l]

            # Accumulate volumetric fluxes from all hillslope columns
            for c in lun.coli[l]:lun.colf[l]
                if col.is_hillslope_column[c] && col.active[c]
                    # Convert mm/s to m3/s using column area (m2)
                    col_area_m2 = grc.area[g] * 1.0e6 * col.wtgcell[c]

                    qflx_surf_vol          = qflx_surf[c] * 1.0e-3 * col_area_m2
                    qflx_drain_perched_vol = qflx_drain_perched[c] * 1.0e-3 * col_area_m2
                    qflx_drain_vol         = qflx_drain[c] * 1.0e-3 * col_area_m2

                    stream_water_volume[l] += (qflx_drain_perched_vol +
                        qflx_drain_vol + qflx_surf_vol) * dtime
                end
            end

            # Remove streamflow discharge
            stream_water_volume[l] -= volumetric_streamflow[l] * dtime

            # Account for negative drainage (via searchforwater in soilhydrology)
            if stream_water_volume[l] < 0.0
                volumetric_streamflow[l] += stream_water_volume[l] / dtime
                stream_water_volume[l] = 0.0
            end

            # Compute stream water depth
            stream_water_depth[l] = stream_water_volume[l] /
                lun.stream_channel_length[l] /
                lun.stream_channel_width[l]
        end
    end

    return nothing
end
