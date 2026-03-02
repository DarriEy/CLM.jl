# ==========================================================================
# Ported from: src/biogeophys/SnowHydrologyMod.F90 (~4070 lines)
# Calculate snow hydrology.
#
# Public functions:
#   new_snow_bulk_density!       — Compute bulk density of newly-fallen snow
#   bulkdiag_new_snow_diagnostics! — Update snow diagnostics for new snow
#   update_state_add_new_snow!   — Update h2osno_no_layers or h2osoi_ice for new snow
#   build_filter_thawed_wetland_thin_snowpack! — Build filter of thawed wetlands with thin snow
#   update_state_remove_snow_thawed_wetlands! — Remove snow from thawed wetlands (state)
#   bulk_remove_snow_thawed_wetlands! — Remove snow from thawed wetlands (bulk)
#   build_filter_snowpack_initialized! — Build filter for snowpack initialization
#   update_state_initialize_snow_pack! — Initialize water state for new snow packs
#   bulk_initialize_snow_pack!   — Initialize bulk quantities for new snow packs
#   update_state_top_layer_fluxes! — Update top snow layer with fluxes
#   bulk_flux_snow_percolation!  — Calculate liquid percolation through snow (bulk)
#   update_state_snow_percolation! — Update h2osoi_liq for snow percolation
#   calc_and_apply_aerosol_fluxes! — Aerosol fluxes through snow pack
#   post_percolation_adjust_layer_thicknesses! — Adjust dz after percolation
#   bulkdiag_snow_water_accumulated_snow! — Update int_snow, reset when no snow
#   sum_flux_add_snow_percolation! — Calculate summed fluxes for percolation
#   snow_compaction!             — Snow layer thickness change due to compaction
#   combine_snow_layers!         — Combine snow layers less than min thickness
#   divide_snow_layers!          — Subdivide snow layers exceeding max thickness
#   zero_empty_snow_layers!      — Set empty snow layers to zero
#   init_snow_layers!            — Initialize cold-start snow layer thickness
#   build_snow_filter!           — Construct snow/no-snow masks
#   snow_capping_excess!         — Determine excess snow for capping
#   init_flux_snow_capping!      — Initialize snow capping fluxes to 0
#   bulk_flux_snow_capping_fluxes! — Calculate snow capping fluxes (bulk)
#   update_state_remove_snow_capping_fluxes! — Remove capping fluxes from state
#   snow_capping_update_dz_and_aerosols! — Adjust dz and aerosols after capping
#   overburden_compaction_anderson1976 — Anderson 1976 overburden compaction
#   overburden_compaction_vionnet2012  — Vionnet 2012 overburden compaction
#   wind_drift_compaction!       — Wind drift compaction
#   combo!                       — Combine two snow elements
#   mass_weighted_snow_radius    — Mass-weighted snow grain radius
#   snow_hydrology_set_control_for_testing! — Set controls for unit testing
# ==========================================================================

# ---- Module-level parameters (from params_type in Fortran) ----

Base.@kwdef mutable struct SnowHydrologyParams
    wimp::Float64 = 0.05                      # water impermeable if porosity < wimp (unitless)
    ssi::Float64 = 0.033                      # irreducible water saturation of snow (unitless)
    drift_gs::Float64 = 0.35                  # wind drift compaction / grain size (unitless)
    eta0_anderson::Float64 = 9.0e5            # viscosity coefficient Anderson1976 (kg*s/m2)
    eta0_vionnet::Float64 = 7.62237e6         # viscosity coefficient Vionnet2012 (kg*s/m2)
    wind_snowcompact_fact::Float64 = 5.0      # reference wind for snow density increase (m/s)
    rho_max::Float64 = 350.0                  # wind drift compaction / max density (kg/m3)
    tau_ref::Float64 = 172800.0               # wind drift compaction / reference time (s)
    scvng_fct_mlt_sf::Float64 = 1.0           # scaling factor for scavenging (-)
    scvng_fct_mlt_bcphi::Float64 = 0.20       # scavenging factor for hydrophilic BC [frc]
    scvng_fct_mlt_bcpho::Float64 = 0.03       # scavenging factor for hydrophobic BC [frc]
    scvng_fct_mlt_dst1::Float64 = 0.02        # scavenging factor dust 1 [frc]
    scvng_fct_mlt_dst2::Float64 = 0.02        # scavenging factor dust 2 [frc]
    scvng_fct_mlt_dst3::Float64 = 0.01        # scavenging factor dust 3 [frc]
    scvng_fct_mlt_dst4::Float64 = 0.01        # scavenging factor dust 4 [frc]
    ceta::Float64 = 250.0                     # overburden compaction constant (kg/m3)
    snw_rds_min::Float64 = 54.526             # min allowed snow effective radius [microns]
    upplim_destruct_metamorph::Float64 = 100.0 # upper limit on destructive metamorphism compaction (kg/m3)
end

const snowhydrology_params = SnowHydrologyParams()

# ---- Module-level constants ----

# Aerosol scavenging factors for OC (public constants in Fortran)
const SCVNG_FCT_MLT_OCPHI = 0.20  # scavenging factor for hydrophilic OC [frc]
const SCVNG_FCT_MLT_OCPHO = 0.03  # scavenging factor for hydrophobic OC [frc]

# Low-temperature snow density method enumerations
const LO_TMP_DNS_SLATER2017 = 2
const LO_TMP_DNS_TRUNCATED_ANDERSON1976 = 1

# Overburden compaction method enumerations
const OVERBURDEN_COMPACTION_ANDERSON1976 = 1
const OVERBURDEN_COMPACTION_VIONNET2012  = 2

# Module-level control variables (Ref for mutability)
const WIND_DEPENDENT_SNOW_DENSITY = Ref(false)
const NEW_SNOW_DENSITY = Ref(LO_TMP_DNS_SLATER2017)
const OVERBURDEN_COMPACTION_METHOD = Ref(OVERBURDEN_COMPACTION_ANDERSON1976)
const OVERBURDEN_COMPRESS_TFACTOR = Ref(0.08)

# Snow layer thickness parameters
const SNOW_DZMIN_1 = Ref(0.010)
const SNOW_DZMIN_2 = Ref(0.015)
const SNOW_DZMAX_L_1 = Ref(0.03)
const SNOW_DZMAX_L_2 = Ref(0.07)
const SNOW_DZMAX_U_1 = Ref(0.02)
const SNOW_DZMAX_U_2 = Ref(0.05)

# Snow layer thickness arrays (initialized by init_snow_layers!)
const SNOW_DZMIN   = Ref(Float64[])
const SNOW_DZMAX_L = Ref(Float64[])
const SNOW_DZMAX_U = Ref(Float64[])

# Snow resetting parameters
const RESET_SNOW = Ref(false)
const RESET_SNOW_GLC = Ref(false)
const RESET_SNOW_GLC_ELA = Ref(1.0e9)
const RESET_SNOW_H2OSNO = 35.0   # mm SWE to reset the snow pack to
const RESET_SNOW_TIMESTEPS_PER_LAYER = 4

# Lake snow additional thickness (from LakeCon)
const LSADZ = 0.03  # additional minimum thickness for lake snow layers (m)

# =========================================================================
# NewSnowBulkDensity
# =========================================================================

"""
    new_snow_bulk_density!(bifall, forc_t, forc_wind, col_gridcell,
                           mask, bounds)

Compute the bulk density of newly-fallen snow (kg/m3).
Uses Alta relationship from Anderson(1976).
"""
function new_snow_bulk_density!(
    bifall::Vector{Float64},        # output: bulk density [kg/m3]
    forc_t::Vector{Float64},        # input: atmospheric temperature [K]
    forc_wind::Vector{Float64},     # input: atmospheric wind speed [m/s] (gridcell)
    col_gridcell::Vector{Int},      # input: column-to-gridcell mapping
    mask::BitVector,                # input: column mask
    bounds::UnitRange{Int}          # input: column bounds
)
    params = snowhydrology_params

    for c in bounds
        mask[c] || continue
        g = col_gridcell[c]

        if forc_t[c] > TFRZ + 2.0
            bifall[c] = 50.0 + 1.7 * (17.0)^1.5
        elseif forc_t[c] > TFRZ - 15.0
            bifall[c] = 50.0 + 1.7 * (forc_t[c] - TFRZ + 15.0)^1.5
        elseif NEW_SNOW_DENSITY[] == LO_TMP_DNS_TRUNCATED_ANDERSON1976
            bifall[c] = 50.0
        elseif NEW_SNOW_DENSITY[] == LO_TMP_DNS_SLATER2017
            if forc_t[c] > TFRZ - 57.55
                t_for_bifall_degC = forc_t[c] - TFRZ
            else
                t_for_bifall_degC = -57.55
            end
            bifall[c] = -(50.0/15.0 + 0.0333*15.0)*t_for_bifall_degC - 0.0333*t_for_bifall_degC^2
        end

        if WIND_DEPENDENT_SNOW_DENSITY[] && forc_wind[g] > 0.1
            bifall[c] = bifall[c] + 266.861 * ((1.0 +
                        tanh(forc_wind[g]/params.wind_snowcompact_fact))/2.0)^8.8
        end
    end
end

# =========================================================================
# BulkDiag_NewSnowDiagnostics
# =========================================================================

"""
    bulkdiag_new_snow_diagnostics!(dz, int_snow, swe_old, frac_sno,
        frac_sno_eff, snow_depth, snomelt_accum,
        dtime, lun_itype_col, urbpoi, snl, bifall, h2osno_total,
        h2osoi_ice, h2osoi_liq, qflx_snow_grnd, qflx_snow_drain,
        mask, bounds, nlevsno)

Update snow depth, fractional area, integrated snow and SWE diagnostics
for newly-fallen snow.
"""
function bulkdiag_new_snow_diagnostics!(
    # Outputs (modified in place)
    dz::Matrix{Float64},            # layer depth [m]
    int_snow::Vector{Float64},      # integrated snowfall [mm H2O]
    swe_old::Matrix{Float64},       # SWE before update [mm H2O]
    frac_sno::Vector{Float64},      # fraction of ground covered by snow
    frac_sno_eff::Vector{Float64},  # effective fraction
    snow_depth::Vector{Float64},    # snow height [m]
    snomelt_accum::Vector{Float64}, # accumulated snow melt for z0m [m H2O]
    # Inputs
    dtime::Float64,                 # time step [s]
    lun_itype_col::Vector{Int},     # landunit type per column
    urbpoi::Vector{Bool},           # urban point flags (landunit-level, indexed by col)
    snl::Vector{Int},               # number of snow layers (negative)
    bifall::Vector{Float64},        # bulk density of new snow [kg/m3]
    h2osno_total::Vector{Float64},  # total snow water [mm H2O]
    h2osoi_ice::Matrix{Float64},    # ice lens [kg/m2]
    h2osoi_liq::Matrix{Float64},    # liquid water [kg/m2]
    qflx_snow_grnd::Vector{Float64}, # snow on ground [mm H2O/s]
    qflx_snow_drain::Vector{Float64}, # drainage from snow pack [mm H2O/s]
    mask::BitVector,                # column mask
    bounds::UnitRange{Int},         # column bounds
    nlevsno::Int                    # max snow layers
)
    newsnow = zeros(Float64, length(snow_depth))
    snowmelt = zeros(Float64, length(snow_depth))
    temp_snow_depth = zeros(Float64, length(snow_depth))

    for c in bounds
        mask[c] || continue

        # Save initial snow depth
        temp_snow_depth[c] = snow_depth[c]

        # Save initial snow water equivalent
        for j in (-nlevsno+1):snl[c]
            jj = j + nlevsno  # Julia index
            swe_old[c, jj] = 0.0
        end
        for j in (snl[c]+1):0
            jj = j + nlevsno
            swe_old[c, jj] = h2osoi_liq[c, jj] + h2osoi_ice[c, jj]
        end

        # All snow falls on ground
        newsnow[c] = qflx_snow_grnd[c] * dtime
        snomelt_accum[c] = max(0.0, snomelt_accum[c] - newsnow[c] * 1.0e-3)

        # Update int_snow
        int_snow[c] = max(int_snow[c], h2osno_total[c])

        # Snowmelt from previous time step
        snowmelt[c] = qflx_snow_drain[c] * dtime
    end

    # Update snow depth and fractional snow cover
    # (Simplified: direct density-based depth update without scf_method polymorphism)
    for c in bounds
        mask[c] || continue

        if lun_itype_col[c] == ISTICE || urbpoi[c]
            # For glaciers and urban, snow depth computed from bifall
            if newsnow[c] > 0.0
                snow_depth[c] = snow_depth[c] + newsnow[c] / bifall[c]
            end
            frac_sno[c] = 1.0
            frac_sno_eff[c] = 1.0
        else
            # For other landunits: snow depth update from bifall
            if newsnow[c] > 0.0
                snow_depth[c] = snow_depth[c] + newsnow[c] / bifall[c]
            end
            # Simple snow cover fraction update
            if snow_depth[c] > 0.0
                frac_sno[c] = min(1.0, snow_depth[c] / 0.05)
                frac_sno_eff[c] = frac_sno[c]
            else
                frac_sno[c] = 0.0
                frac_sno_eff[c] = 0.0
            end
        end
    end

    # Reset int_snow if snow pack started at 0
    for c in bounds
        mask[c] || continue
        if h2osno_total[c] == 0.0 && newsnow[c] > 0.0
            int_snow[c] = 0.0
        end
    end

    # Add newsnow to int_snow
    for c in bounds
        mask[c] || continue
        int_snow[c] = int_snow[c] + newsnow[c]
    end

    # Update top snow layer thickness
    for c in bounds
        mask[c] || continue
        if snl[c] < 0
            dz_snowf = (snow_depth[c] - temp_snow_depth[c]) / dtime
            jj = (snl[c] + 1) + nlevsno  # Julia index for top snow layer
            dz[c, jj] = dz[c, jj] + dz_snowf * dtime
        end
    end
end

# =========================================================================
# UpdateState_AddNewSnow
# =========================================================================

"""
    update_state_add_new_snow!(h2osno_no_layers, h2osoi_ice,
        dtime, snl, qflx_snow_grnd, mask, bounds, nlevsno)

Update h2osno_no_layers or h2osoi_ice based on new snow.
"""
function update_state_add_new_snow!(
    h2osno_no_layers::Vector{Float64},
    h2osoi_ice::Matrix{Float64},
    dtime::Float64,
    snl::Vector{Int},
    qflx_snow_grnd::Vector{Float64},
    mask::BitVector,
    bounds::UnitRange{Int},
    nlevsno::Int
)
    for c in bounds
        mask[c] || continue
        if snl[c] == 0
            h2osno_no_layers[c] = h2osno_no_layers[c] + qflx_snow_grnd[c] * dtime
        else
            jj = (snl[c] + 1) + nlevsno  # top snow layer Julia index
            h2osoi_ice[c, jj] = h2osoi_ice[c, jj] + qflx_snow_grnd[c] * dtime
        end
    end
end

# =========================================================================
# BuildFilter_ThawedWetlandThinSnowpack
# =========================================================================

"""
    build_filter_thawed_wetland_thin_snowpack!(mask_out,
        t_grnd, lun_itype_col, snl, mask_nolake, bounds)

Build a column-level mask of thawed wetland columns with thin (no-layer) snow pack.
"""
function build_filter_thawed_wetland_thin_snowpack!(
    mask_out::BitVector,
    t_grnd::Vector{Float64},
    lun_itype_col::Vector{Int},
    snl::Vector{Int},
    mask_nolake::BitVector,
    bounds::UnitRange{Int}
)
    fill!(mask_out, false)
    for c in bounds
        mask_nolake[c] || continue
        if lun_itype_col[c] == ISTWET && t_grnd[c] > TFRZ && snl[c] == 0
            mask_out[c] = true
        end
    end
end

# =========================================================================
# UpdateState_RemoveSnowFromThawedWetlands
# =========================================================================

"""
    update_state_remove_snow_thawed_wetlands!(h2osno_no_layers, mask, bounds)

Remove snow from thawed wetlands — set h2osno_no_layers to zero.
"""
function update_state_remove_snow_thawed_wetlands!(
    h2osno_no_layers::Vector{Float64},
    mask::BitVector,
    bounds::UnitRange{Int}
)
    for c in bounds
        mask[c] || continue
        h2osno_no_layers[c] = 0.0
    end
end

# =========================================================================
# Bulk_RemoveSnowFromThawedWetlands
# =========================================================================

"""
    bulk_remove_snow_thawed_wetlands!(snow_depth, mask, bounds)

Remove snow from thawed wetlands — set snow_depth to zero.
"""
function bulk_remove_snow_thawed_wetlands!(
    snow_depth::Vector{Float64},
    mask::BitVector,
    bounds::UnitRange{Int}
)
    for c in bounds
        mask[c] || continue
        snow_depth[c] = 0.0
    end
end

# =========================================================================
# BuildFilter_SnowpackInitialized
# =========================================================================

"""
    build_filter_snowpack_initialized!(mask_out,
        snl, lun_itype_col, frac_sno_eff, snow_depth, qflx_snow_grnd,
        mask, bounds)

Build a column-level mask of columns where an explicit snow pack needs initialization.
"""
function build_filter_snowpack_initialized!(
    mask_out::BitVector,
    snl::Vector{Int},
    lun_itype_col::Vector{Int},
    frac_sno_eff::Vector{Float64},
    snow_depth::Vector{Float64},
    qflx_snow_grnd::Vector{Float64},
    mask::BitVector,
    bounds::UnitRange{Int}
)
    dzmin = SNOW_DZMIN[]

    fill!(mask_out, false)
    for c in bounds
        mask[c] || continue

        if lun_itype_col[c] == ISTDLAK
            mask_out[c] = (snl[c] == 0 &&
                           frac_sno_eff[c] * snow_depth[c] >= (dzmin[1] + LSADZ) &&
                           qflx_snow_grnd[c] > 0.0)
        else
            mask_out[c] = (snl[c] == 0 &&
                           frac_sno_eff[c] * snow_depth[c] >= dzmin[1])
        end
    end
end

# =========================================================================
# UpdateState_InitializeSnowPack
# =========================================================================

"""
    update_state_initialize_snow_pack!(h2osno_no_layers, h2osoi_ice, h2osoi_liq,
        mask, bounds, nlevsno)

Initialize water state variables for columns with newly-initialized snow pack.
"""
function update_state_initialize_snow_pack!(
    h2osno_no_layers::Vector{Float64},
    h2osoi_ice::Matrix{Float64},
    h2osoi_liq::Matrix{Float64},
    mask::BitVector,
    bounds::UnitRange{Int},
    nlevsno::Int
)
    jj_zero = 0 + nlevsno  # Julia index for layer 0
    for c in bounds
        mask[c] || continue
        h2osoi_ice[c, jj_zero] = h2osno_no_layers[c]
        h2osoi_liq[c, jj_zero] = 0.0
        h2osno_no_layers[c] = 0.0
    end
end

# =========================================================================
# Bulk_InitializeSnowPack
# =========================================================================

"""
    bulk_initialize_snow_pack!(snl, zi, dz, z, t_soisno, frac_iceold,
        snomelt_accum, forc_t, snow_depth, mask, bounds, nlevsno)

Initialize an explicit snow pack in columns warranted by snow depth.
"""
function bulk_initialize_snow_pack!(
    snl::Vector{Int},
    zi::Matrix{Float64},
    dz::Matrix{Float64},
    z::Matrix{Float64},
    t_soisno::Matrix{Float64},
    frac_iceold::Matrix{Float64},
    snomelt_accum::Vector{Float64},
    forc_t::Vector{Float64},
    snow_depth::Vector{Float64},
    mask::BitVector,
    bounds::UnitRange{Int},
    nlevsno::Int
)
    jj_zero = 0 + nlevsno  # Julia index for layer 0
    jj_m1   = -1 + nlevsno # Julia index for layer -1

    for c in bounds
        mask[c] || continue

        snl[c] = -1
        dz[c, jj_zero] = snow_depth[c]
        z[c, jj_zero] = -0.5 * dz[c, jj_zero]
        zi[c, jj_m1] = -dz[c, jj_zero]

        t_soisno[c, jj_zero] = min(TFRZ, forc_t[c])
        frac_iceold[c, jj_zero] = 1.0
        snomelt_accum[c] = 0.0
    end
end

# =========================================================================
# UpdateState_TopLayerFluxes
# =========================================================================

"""
    update_state_top_layer_fluxes!(h2osoi_ice, h2osoi_liq,
        dtime, snl, frac_sno_eff,
        qflx_soliddew_to_top_layer, qflx_solidevap_from_top_layer,
        qflx_liq_grnd, qflx_liqdew_to_top_layer, qflx_liqevap_from_top_layer,
        mask_snow, bounds, nlevsno)

Update top layer of snow pack with various fluxes into and out of the top layer.
"""
function update_state_top_layer_fluxes!(
    h2osoi_ice::Matrix{Float64},
    h2osoi_liq::Matrix{Float64},
    dtime::Float64,
    snl::Vector{Int},
    frac_sno_eff::Vector{Float64},
    qflx_soliddew_to_top_layer::Vector{Float64},
    qflx_solidevap_from_top_layer::Vector{Float64},
    qflx_liq_grnd::Vector{Float64},
    qflx_liqdew_to_top_layer::Vector{Float64},
    qflx_liqevap_from_top_layer::Vector{Float64},
    mask_snow::BitVector,
    bounds::UnitRange{Int},
    nlevsno::Int
)
    for c in bounds
        mask_snow[c] || continue

        lev_top = snl[c] + 1
        jj = lev_top + nlevsno  # Julia index for top snow layer

        h2osoi_ice[c, jj] = h2osoi_ice[c, jj] +
            frac_sno_eff[c] * (qflx_soliddew_to_top_layer[c] -
            qflx_solidevap_from_top_layer[c]) * dtime

        h2osoi_liq[c, jj] = h2osoi_liq[c, jj] +
            frac_sno_eff[c] * (qflx_liq_grnd[c] + qflx_liqdew_to_top_layer[c] -
            qflx_liqevap_from_top_layer[c]) * dtime

        # Truncate near-zero negatives
        if abs(h2osoi_ice[c, jj]) < 1.0e-10
            h2osoi_ice[c, jj] = max(0.0, h2osoi_ice[c, jj])
        end
        if abs(h2osoi_liq[c, jj]) < 1.0e-10
            h2osoi_liq[c, jj] = max(0.0, h2osoi_liq[c, jj])
        end

        # Error check
        if h2osoi_ice[c, jj] < 0.0
            error("update_state_top_layer_fluxes!: h2osoi_ice negative at c=$c, lev_top=$lev_top, val=$(h2osoi_ice[c, jj])")
        end
        if h2osoi_liq[c, jj] < 0.0
            error("update_state_top_layer_fluxes!: h2osoi_liq negative at c=$c, lev_top=$lev_top, val=$(h2osoi_liq[c, jj])")
        end
    end
end

# =========================================================================
# BulkFlux_SnowPercolation
# =========================================================================

"""
    bulk_flux_snow_percolation!(qflx_snow_percolation,
        dtime, snl, dz, frac_sno_eff, h2osoi_ice, h2osoi_liq,
        mask_snow, bounds, nlevsno)

Calculate liquid percolation through the snow pack for bulk water.
qflx_snow_percolation[c,j] gives percolation out of bottom of layer j.
"""
function bulk_flux_snow_percolation!(
    qflx_snow_percolation::Matrix{Float64},
    dtime::Float64,
    snl::Vector{Int},
    dz::Matrix{Float64},
    frac_sno_eff::Vector{Float64},
    h2osoi_ice::Matrix{Float64},
    h2osoi_liq::Matrix{Float64},
    mask_snow::BitVector,
    bounds::UnitRange{Int},
    nlevsno::Int
)
    params = snowhydrology_params

    # Temporary arrays for porosity/volume calculations
    nc = length(snl)
    vol_liq = zeros(Float64, nc, nlevsno)
    vol_ice = zeros(Float64, nc, nlevsno)
    eff_porosity = zeros(Float64, nc, nlevsno)

    # Porosity and partial volume
    for j in (-nlevsno+1):0
        jj = j + nlevsno
        for c in bounds
            mask_snow[c] || continue
            if j >= snl[c] + 1
                vol_ice[c, jj] = min(1.0, h2osoi_ice[c, jj] / (dz[c, jj] * frac_sno_eff[c] * DENICE))
                eff_porosity[c, jj] = 1.0 - vol_ice[c, jj]
                vol_liq[c, jj] = min(eff_porosity[c, jj],
                    h2osoi_liq[c, jj] / (dz[c, jj] * frac_sno_eff[c] * DENH2O))
            end
        end
    end

    # Calculate percolation
    for j in (-nlevsno+1):0
        jj = j + nlevsno
        for c in bounds
            mask_snow[c] || continue
            if j >= snl[c] + 1
                if j <= -1
                    jj_next = (j + 1) + nlevsno
                    if eff_porosity[c, jj] < params.wimp || eff_porosity[c, jj_next] < params.wimp
                        qflx_snow_percolation[c, jj] = 0.0
                    else
                        qflx_snow_percolation[c, jj] = max(0.0,
                            (vol_liq[c, jj] - params.ssi * eff_porosity[c, jj]) *
                            dz[c, jj] * frac_sno_eff[c])
                        qflx_snow_percolation[c, jj] = min(qflx_snow_percolation[c, jj],
                            (1.0 - vol_ice[c, jj_next] - vol_liq[c, jj_next]) *
                            dz[c, jj_next] * frac_sno_eff[c])
                    end
                else  # j == 0, bottom layer
                    qflx_snow_percolation[c, jj] = max(0.0,
                        (vol_liq[c, jj] - params.ssi * eff_porosity[c, jj]) *
                        dz[c, jj] * frac_sno_eff[c])
                end
                # Convert from m to mm H2O/s
                qflx_snow_percolation[c, jj] = (qflx_snow_percolation[c, jj] * 1000.0) / dtime
            end
        end
    end
end

# =========================================================================
# UpdateState_SnowPercolation
# =========================================================================

"""
    update_state_snow_percolation!(h2osoi_liq,
        dtime, snl, qflx_snow_percolation, mask_snow, bounds, nlevsno)

Update h2osoi_liq for snow percolation (bulk or one tracer).
"""
function update_state_snow_percolation!(
    h2osoi_liq::Matrix{Float64},
    dtime::Float64,
    snl::Vector{Int},
    qflx_snow_percolation::Matrix{Float64},
    mask_snow::BitVector,
    bounds::UnitRange{Int},
    nlevsno::Int
)
    for j in (-nlevsno+1):0
        jj = j + nlevsno
        for c in bounds
            mask_snow[c] || continue
            if j >= snl[c] + 1
                # Add percolation from layer above
                if j >= snl[c] + 2
                    jj_above = (j - 1) + nlevsno
                    h2osoi_liq[c, jj] = h2osoi_liq[c, jj] + qflx_snow_percolation[c, jj_above] * dtime
                end
                # Subtract percolation out of this layer
                h2osoi_liq[c, jj] = h2osoi_liq[c, jj] - qflx_snow_percolation[c, jj] * dtime
            end
        end
    end
end

# =========================================================================
# CalcAndApplyAerosolFluxes
# =========================================================================

"""
    calc_and_apply_aerosol_fluxes!(aerosol, dtime, snl, h2osoi_ice, h2osoi_liq,
        qflx_snow_percolation, mask_snow, bounds, nlevsno)

Calculate and apply fluxes of aerosols through the snow pack.
"""
function calc_and_apply_aerosol_fluxes!(
    aerosol::AerosolData,
    dtime::Float64,
    snl::Vector{Int},
    h2osoi_ice::Matrix{Float64},
    h2osoi_liq::Matrix{Float64},
    qflx_snow_percolation::Matrix{Float64},
    mask_snow::BitVector,
    bounds::UnitRange{Int},
    nlevsno::Int
)
    params = snowhydrology_params
    nc = length(snl)

    # Initialize incoming fluxes
    qin_bc_phi = zeros(Float64, nc)
    qin_bc_pho = zeros(Float64, nc)
    qin_oc_phi = zeros(Float64, nc)
    qin_oc_pho = zeros(Float64, nc)
    qin_dst1   = zeros(Float64, nc)
    qin_dst2   = zeros(Float64, nc)
    qin_dst3   = zeros(Float64, nc)
    qin_dst4   = zeros(Float64, nc)

    for j in (-nlevsno+1):0
        jj = j + nlevsno
        for c in bounds
            mask_snow[c] || continue
            if j >= snl[c] + 1

                # Add incoming flux
                aerosol.mss_bcphi_col[c, jj] += qin_bc_phi[c] * dtime
                aerosol.mss_bcpho_col[c, jj] += qin_bc_pho[c] * dtime
                aerosol.mss_ocphi_col[c, jj] += qin_oc_phi[c] * dtime
                aerosol.mss_ocpho_col[c, jj] += qin_oc_pho[c] * dtime
                aerosol.mss_dst1_col[c, jj]  += qin_dst1[c] * dtime
                aerosol.mss_dst2_col[c, jj]  += qin_dst2[c] * dtime
                aerosol.mss_dst3_col[c, jj]  += qin_dst3[c] * dtime
                aerosol.mss_dst4_col[c, jj]  += qin_dst4[c] * dtime

                # Mass of ice+water (avoid division by zero)
                mss_liqice = h2osoi_liq[c, jj] + h2osoi_ice[c, jj]
                if mss_liqice < 1.0e-30
                    mss_liqice = 1.0e-30
                end

                # Compute outgoing fluxes for each aerosol species
                # BCPHI
                qout = qflx_snow_percolation[c, jj] * params.scvng_fct_mlt_sf *
                       params.scvng_fct_mlt_bcphi * (aerosol.mss_bcphi_col[c, jj] / mss_liqice)
                if qout * dtime > aerosol.mss_bcphi_col[c, jj]
                    qout = aerosol.mss_bcphi_col[c, jj] / dtime
                    aerosol.mss_bcphi_col[c, jj] = 0.0
                else
                    aerosol.mss_bcphi_col[c, jj] -= qout * dtime
                end
                qin_bc_phi[c] = qout

                # BCPHO
                qout = qflx_snow_percolation[c, jj] * params.scvng_fct_mlt_sf *
                       params.scvng_fct_mlt_bcpho * (aerosol.mss_bcpho_col[c, jj] / mss_liqice)
                if qout * dtime > aerosol.mss_bcpho_col[c, jj]
                    qout = aerosol.mss_bcpho_col[c, jj] / dtime
                    aerosol.mss_bcpho_col[c, jj] = 0.0
                else
                    aerosol.mss_bcpho_col[c, jj] -= qout * dtime
                end
                qin_bc_pho[c] = qout

                # OCPHI
                qout = qflx_snow_percolation[c, jj] * params.scvng_fct_mlt_sf *
                       SCVNG_FCT_MLT_OCPHI * (aerosol.mss_ocphi_col[c, jj] / mss_liqice)
                if qout * dtime > aerosol.mss_ocphi_col[c, jj]
                    qout = aerosol.mss_ocphi_col[c, jj] / dtime
                    aerosol.mss_ocphi_col[c, jj] = 0.0
                else
                    aerosol.mss_ocphi_col[c, jj] -= qout * dtime
                end
                qin_oc_phi[c] = qout

                # OCPHO
                qout = qflx_snow_percolation[c, jj] * params.scvng_fct_mlt_sf *
                       SCVNG_FCT_MLT_OCPHO * (aerosol.mss_ocpho_col[c, jj] / mss_liqice)
                if qout * dtime > aerosol.mss_ocpho_col[c, jj]
                    qout = aerosol.mss_ocpho_col[c, jj] / dtime
                    aerosol.mss_ocpho_col[c, jj] = 0.0
                else
                    aerosol.mss_ocpho_col[c, jj] -= qout * dtime
                end
                qin_oc_pho[c] = qout

                # DST1
                qout = qflx_snow_percolation[c, jj] * params.scvng_fct_mlt_sf *
                       params.scvng_fct_mlt_dst1 * (aerosol.mss_dst1_col[c, jj] / mss_liqice)
                if qout * dtime > aerosol.mss_dst1_col[c, jj]
                    qout = aerosol.mss_dst1_col[c, jj] / dtime
                    aerosol.mss_dst1_col[c, jj] = 0.0
                else
                    aerosol.mss_dst1_col[c, jj] -= qout * dtime
                end
                qin_dst1[c] = qout

                # DST2
                qout = qflx_snow_percolation[c, jj] * params.scvng_fct_mlt_sf *
                       params.scvng_fct_mlt_dst2 * (aerosol.mss_dst2_col[c, jj] / mss_liqice)
                if qout * dtime > aerosol.mss_dst2_col[c, jj]
                    qout = aerosol.mss_dst2_col[c, jj] / dtime
                    aerosol.mss_dst2_col[c, jj] = 0.0
                else
                    aerosol.mss_dst2_col[c, jj] -= qout * dtime
                end
                qin_dst2[c] = qout

                # DST3
                qout = qflx_snow_percolation[c, jj] * params.scvng_fct_mlt_sf *
                       params.scvng_fct_mlt_dst3 * (aerosol.mss_dst3_col[c, jj] / mss_liqice)
                if qout * dtime > aerosol.mss_dst3_col[c, jj]
                    qout = aerosol.mss_dst3_col[c, jj] / dtime
                    aerosol.mss_dst3_col[c, jj] = 0.0
                else
                    aerosol.mss_dst3_col[c, jj] -= qout * dtime
                end
                qin_dst3[c] = qout

                # DST4
                qout = qflx_snow_percolation[c, jj] * params.scvng_fct_mlt_sf *
                       params.scvng_fct_mlt_dst4 * (aerosol.mss_dst4_col[c, jj] / mss_liqice)
                if qout * dtime > aerosol.mss_dst4_col[c, jj]
                    qout = aerosol.mss_dst4_col[c, jj] / dtime
                    aerosol.mss_dst4_col[c, jj] = 0.0
                else
                    aerosol.mss_dst4_col[c, jj] -= qout * dtime
                end
                qin_dst4[c] = qout
            end
        end
    end
end

# =========================================================================
# PostPercolation_AdjustLayerThicknesses
# =========================================================================

"""
    post_percolation_adjust_layer_thicknesses!(dz,
        snl, h2osoi_ice, h2osoi_liq, mask_snow, bounds, nlevsno)

Adjust layer thickness for any water+ice content changes after percolation.
"""
function post_percolation_adjust_layer_thicknesses!(
    dz::Matrix{Float64},
    snl::Vector{Int},
    h2osoi_ice::Matrix{Float64},
    h2osoi_liq::Matrix{Float64},
    mask_snow::BitVector,
    bounds::UnitRange{Int},
    nlevsno::Int
)
    for j in (-nlevsno+1):0
        jj = j + nlevsno
        for c in bounds
            mask_snow[c] || continue
            if j >= snl[c] + 1
                dz[c, jj] = max(dz[c, jj],
                    h2osoi_liq[c, jj] / DENH2O + h2osoi_ice[c, jj] / DENICE)
            end
        end
    end
end

# =========================================================================
# BulkDiag_SnowWaterAccumulatedSnow
# =========================================================================

"""
    bulkdiag_snow_water_accumulated_snow!(int_snow, frac_sno, snow_depth,
        dtime, frac_sno_eff, qflx_soliddew_to_top_layer, qflx_liqdew_to_top_layer,
        qflx_liq_grnd, h2osno_no_layers,
        mask_snow, mask_nosnow, bounds)

Update int_snow, and reset accumulated snow when no snow present.
"""
function bulkdiag_snow_water_accumulated_snow!(
    int_snow::Vector{Float64},
    frac_sno::Vector{Float64},
    snow_depth::Vector{Float64},
    dtime::Float64,
    frac_sno_eff::Vector{Float64},
    qflx_soliddew_to_top_layer::Vector{Float64},
    qflx_liqdew_to_top_layer::Vector{Float64},
    qflx_liq_grnd::Vector{Float64},
    h2osno_no_layers::Vector{Float64},
    mask_snow::BitVector,
    mask_nosnow::BitVector,
    bounds::UnitRange{Int}
)
    for c in bounds
        if mask_snow[c]
            int_snow[c] = int_snow[c] + frac_sno_eff[c] *
                (qflx_soliddew_to_top_layer[c] + qflx_liqdew_to_top_layer[c] +
                 qflx_liq_grnd[c]) * dtime
        end
    end

    for c in bounds
        if mask_nosnow[c]
            if h2osno_no_layers[c] <= 0.0
                int_snow[c] = 0.0
                frac_sno[c] = 0.0
                snow_depth[c] = 0.0
            end
        end
    end
end

# =========================================================================
# SumFlux_AddSnowPercolation
# =========================================================================

"""
    sum_flux_add_snow_percolation!(qflx_snow_drain, qflx_rain_plus_snomelt,
        frac_sno_eff, qflx_snow_percolation_bottom, qflx_liq_grnd, qflx_snomelt,
        mask_snow, mask_nosnow, bounds)

Calculate summed fluxes accounting for snow percolation.
"""
function sum_flux_add_snow_percolation!(
    qflx_snow_drain::Vector{Float64},
    qflx_rain_plus_snomelt::Vector{Float64},
    frac_sno_eff::Vector{Float64},
    qflx_snow_percolation_bottom::Vector{Float64},
    qflx_liq_grnd::Vector{Float64},
    qflx_snomelt::Vector{Float64},
    mask_snow::BitVector,
    mask_nosnow::BitVector,
    bounds::UnitRange{Int}
)
    for c in bounds
        if mask_snow[c]
            qflx_snow_drain[c] = qflx_snow_drain[c] + qflx_snow_percolation_bottom[c]
            qflx_rain_plus_snomelt[c] = qflx_snow_percolation_bottom[c] +
                (1.0 - frac_sno_eff[c]) * qflx_liq_grnd[c]
        end
    end

    for c in bounds
        if mask_nosnow[c]
            qflx_snow_drain[c] = qflx_snomelt[c]
            qflx_rain_plus_snomelt[c] = qflx_liq_grnd[c] + qflx_snomelt[c]
        end
    end
end

# =========================================================================
# SnowCompaction
# =========================================================================

"""
    snow_compaction!(dz, dtime, snl, t_soisno, h2osoi_ice, h2osoi_liq,
        imelt, frac_sno, frac_h2osfc, swe_old, int_snow, frac_iceold,
        forc_wind, col_gridcell, col_landunit,
        lakpoi, urbpoi,
        mask_snow, bounds, nlevsno)

Determine change in snow layer thickness due to compaction and settling.
"""
function snow_compaction!(
    dz::Matrix{Float64},
    dtime::Float64,
    snl::Vector{Int},
    t_soisno::Matrix{Float64},
    h2osoi_ice::Matrix{Float64},
    h2osoi_liq::Matrix{Float64},
    imelt::Matrix{Int},
    frac_sno::Vector{Float64},
    frac_h2osfc::Vector{Float64},
    swe_old::Matrix{Float64},
    int_snow::Vector{Float64},
    frac_iceold::Matrix{Float64},
    forc_wind::Vector{Float64},
    col_gridcell::Vector{Int},
    col_landunit::Vector{Int},
    lakpoi::Vector{Bool},
    urbpoi::Vector{Bool},
    mask_snow::BitVector,
    bounds::UnitRange{Int},
    nlevsno::Int
)
    params = snowhydrology_params

    # Compaction constants
    c3 = 2.777e-6   # [1/s]
    c4 = 0.04        # [1/K]
    c5 = 2.0

    nc = length(snl)
    burden = zeros(Float64, nc)
    zpseudo = zeros(Float64, nc)
    mobile = trues(nc)

    for c in bounds
        mask_snow[c] || continue
        burden[c] = 0.0
        zpseudo[c] = 0.0
        mobile[c] = true
    end

    for j in (-nlevsno+1):0
        jj = j + nlevsno
        for c in bounds
            mask_snow[c] || continue
            if j >= snl[c] + 1
                g = col_gridcell[c]

                wx = h2osoi_ice[c, jj] + h2osoi_liq[c, jj]
                void = 1.0 - (h2osoi_ice[c, jj] / DENICE + h2osoi_liq[c, jj] / DENH2O) /
                       (frac_sno[c] * dz[c, jj])

                if void > 0.001 && h2osoi_ice[c, jj] > 0.1
                    bi = h2osoi_ice[c, jj] / (frac_sno[c] * dz[c, jj])
                    fi = h2osoi_ice[c, jj] / wx
                    td = TFRZ - t_soisno[c, jj]
                    dexpf = exp(-c4 * td)

                    # Destructive metamorphism
                    ddz1 = -c3 * dexpf
                    if bi > params.upplim_destruct_metamorph
                        ddz1 = ddz1 * exp(-46.0e-3 * (bi - params.upplim_destruct_metamorph))
                    end

                    # Liquid water term
                    if h2osoi_liq[c, jj] > 0.01 * dz[c, jj] * frac_sno[c]
                        ddz1 = ddz1 * c5
                    end

                    # Overburden compaction
                    if OVERBURDEN_COMPACTION_METHOD[] == OVERBURDEN_COMPACTION_ANDERSON1976
                        ddz2 = overburden_compaction_anderson1976(burden[c], wx, td, bi)
                    elseif OVERBURDEN_COMPACTION_METHOD[] == OVERBURDEN_COMPACTION_VIONNET2012
                        ddz2 = overburden_compaction_vionnet2012(
                            h2osoi_liq[c, jj], dz[c, jj], burden[c], wx, td, bi)
                    else
                        error("Unknown overburden_compaction_method")
                    end

                    # Melt compaction
                    if imelt[c, jj] == 1
                        l = col_landunit[c]
                        if !lakpoi[l] && !urbpoi[l]
                            ddz3 = max(0.0, min(1.0, (swe_old[c, jj] - wx) / wx))
                            if (swe_old[c, jj] - wx) > 0.0
                                # Compute total snow water
                                wsum = 0.0
                                for jj2 in (snl[c]+1):0
                                    jj2_idx = jj2 + nlevsno
                                    wsum += h2osoi_liq[c, jj2_idx] + h2osoi_ice[c, jj2_idx]
                                end
                                # Simple fractional snow during melt
                                if int_snow[c] > 0.0
                                    fsno_melt = min(1.0, wsum / int_snow[c])
                                else
                                    fsno_melt = 1.0
                                end
                                if (fsno_melt + frac_h2osfc[c]) > 1.0
                                    fsno_melt = 1.0 - frac_h2osfc[c]
                                end
                                ddz3 = ddz3 - max(0.0, (fsno_melt - frac_sno[c]) / frac_sno[c])
                            end
                            ddz3 = -1.0 / dtime * ddz3
                        else
                            ddz3 = -1.0 / dtime * max(0.0,
                                (frac_iceold[c, jj] - fi) / frac_iceold[c, jj])
                        end
                    else
                        ddz3 = 0.0
                    end

                    # Wind drift compaction
                    if WIND_DEPENDENT_SNOW_DENSITY[]
                        ddz4 = Ref(0.0)
                        wind_drift_compaction!(bi, forc_wind[g], dz[c, jj],
                            zpseudo, c, mobile, c, ddz4)
                        ddz4_val = ddz4[]
                    else
                        ddz4_val = 0.0
                    end

                    # Total compaction rate
                    pdzdtc = ddz1 + ddz2 + ddz3 + ddz4_val

                    # Apply compaction
                    dz[c, jj] = max(
                        dz[c, jj] * (1.0 + pdzdtc * dtime),
                        (h2osoi_ice[c, jj] / DENICE + h2osoi_liq[c, jj] / DENH2O) / frac_sno[c])
                else
                    mobile[c] = false
                end

                # Pressure of overlying snow
                burden[c] = burden[c] + wx
            end
        end
    end
end

# =========================================================================
# OverburdenCompactionAnderson1976
# =========================================================================

"""
    overburden_compaction_anderson1976(burden, wx, td, bi) -> Float64

Compute snow overburden compaction using Anderson 1976 formula.
"""
function overburden_compaction_anderson1976(burden::Float64, wx::Float64,
                                            td::Float64, bi::Float64)
    c2 = 23.0e-3  # [m3/kg]
    return -(burden + wx / 2.0) * exp(-OVERBURDEN_COMPRESS_TFACTOR[] * td - c2 * bi) /
           snowhydrology_params.eta0_anderson
end

# =========================================================================
# OverburdenCompactionVionnet2012
# =========================================================================

"""
    overburden_compaction_vionnet2012(h2osoi_liq, dz, burden, wx, td, bi) -> Float64

Compute snow overburden compaction using Vionnet et al. 2012 formula.
"""
function overburden_compaction_vionnet2012(h2osoi_liq::Float64, dz::Float64,
                                           burden::Float64, wx::Float64,
                                           td::Float64, bi::Float64)
    params = snowhydrology_params
    aeta = 0.1
    beta = 0.023

    f1 = 1.0 / (1.0 + 60.0 * h2osoi_liq / (DENH2O * dz))
    f2 = 4.0  # fixed to maximum value
    eta = f1 * f2 * (bi / params.ceta) * exp(aeta * td + beta * bi) * params.eta0_vionnet
    return -(burden + wx / 2.0) / eta
end

# =========================================================================
# WindDriftCompaction
# =========================================================================

"""
    wind_drift_compaction!(bi, forc_wind, dz,
        zpseudo, zpseudo_idx, mobile, mobile_idx, compaction_rate)

Compute wind drift compaction for a single column and level.
Updates zpseudo and mobile. Must be called top-to-bottom.
"""
function wind_drift_compaction!(
    bi::Float64,
    forc_wind::Float64,
    dz::Float64,
    zpseudo::Vector{Float64}, zpseudo_idx::Int,
    mobile::Vector{Bool}, mobile_idx::Int,
    compaction_rate::Ref{Float64}
)
    params = snowhydrology_params
    rho_min = 50.0
    drift_sph = 1.0

    if mobile[mobile_idx]
        Frho = 1.25 - 0.0042 * (max(rho_min, bi) - rho_min)
        MO = 0.34 * (-0.583 * params.drift_gs - 0.833 * drift_sph + 0.833) + 0.66 * Frho
        SI = -2.868 * exp(-0.085 * forc_wind) + 1.0 + MO

        if SI > 0.0
            SI = min(SI, 3.25)
            zpseudo[zpseudo_idx] = zpseudo[zpseudo_idx] + 0.5 * dz * (3.25 - SI)
            gamma_drift = SI * exp(-zpseudo[zpseudo_idx] / 0.1)
            tau_inverse = gamma_drift / params.tau_ref
            compaction_rate[] = -max(0.0, params.rho_max - bi) * tau_inverse
            zpseudo[zpseudo_idx] = zpseudo[zpseudo_idx] + 0.5 * dz * (3.25 - SI)
        else
            mobile[mobile_idx] = false
            compaction_rate[] = 0.0
        end
    else
        compaction_rate[] = 0.0
    end
end

# =========================================================================
# CombineSnowLayers
# =========================================================================

"""
    combine_snow_layers!(snl, dz, zi, z, t_soisno,
        h2osoi_ice, h2osoi_liq, h2osno_no_layers,
        snow_depth, frac_sno, frac_sno_eff, int_snow,
        snw_rds, aerosol,
        lun_itype, urbpoi,
        mask_snow, bounds, nlevsno)

Combine snow layers that are less than minimum thickness or mass.
"""
function combine_snow_layers!(
    snl::Vector{Int},
    dz::Matrix{Float64},
    zi::Matrix{Float64},
    z::Matrix{Float64},
    t_soisno::Matrix{Float64},
    h2osoi_ice::Matrix{Float64},
    h2osoi_liq::Matrix{Float64},
    h2osno_no_layers::Vector{Float64},
    snow_depth::Vector{Float64},
    frac_sno::Vector{Float64},
    frac_sno_eff::Vector{Float64},
    int_snow::Vector{Float64},
    snw_rds::Matrix{Float64},
    aerosol::AerosolData,
    lun_itype::Vector{Int},     # landunit type indexed by column
    urbpoi::Vector{Bool},       # urban point indexed by landunit
    col_landunit::Vector{Int},  # column-to-landunit mapping
    mask_snow::BitVector,
    bounds::UnitRange{Int},
    nlevsno::Int
)
    dzmin = SNOW_DZMIN[]
    nc = length(snl)

    # Determine dzmin adjustments for lake
    dzminloc = copy(dzmin)

    for c in bounds
        mask_snow[c] || continue

        msn_old = snl[c]

        # Remove thin layers with very little ice
        j = msn_old + 1
        while j <= 0
            jj = j + nlevsno
            if h2osoi_ice[c, jj] <= 0.01
                l = col_landunit[c]
                if j < 0 || (lun_itype[c] == ISTSOIL || urbpoi[l] || lun_itype[c] == ISTCROP)
                    # Transfer water to layer below
                    jj_next = (j + 1) + nlevsno
                    if jj_next <= size(h2osoi_liq, 2)
                        h2osoi_liq[c, jj_next] += h2osoi_liq[c, jj]
                        h2osoi_ice[c, jj_next] += h2osoi_ice[c, jj]
                    end
                end

                if j < 0
                    jj_next = (j + 1) + nlevsno
                    dz[c, jj_next] += dz[c, jj]
                    aerosol.mss_bcphi_col[c, jj_next] += aerosol.mss_bcphi_col[c, jj]
                    aerosol.mss_bcpho_col[c, jj_next] += aerosol.mss_bcpho_col[c, jj]
                    aerosol.mss_ocphi_col[c, jj_next] += aerosol.mss_ocphi_col[c, jj]
                    aerosol.mss_ocpho_col[c, jj_next] += aerosol.mss_ocpho_col[c, jj]
                    aerosol.mss_dst1_col[c, jj_next]  += aerosol.mss_dst1_col[c, jj]
                    aerosol.mss_dst2_col[c, jj_next]  += aerosol.mss_dst2_col[c, jj]
                    aerosol.mss_dst3_col[c, jj_next]  += aerosol.mss_dst3_col[c, jj]
                    aerosol.mss_dst4_col[c, jj_next]  += aerosol.mss_dst4_col[c, jj]
                end

                # Shift all elements above this down one
                if j > snl[c] + 1 && snl[c] < -1
                    for i in j:-1:(snl[c]+2)
                        ii = i + nlevsno
                        ii_above = (i - 1) + nlevsno
                        h2osoi_liq[c, ii] = h2osoi_liq[c, ii_above]
                        h2osoi_ice[c, ii] = h2osoi_ice[c, ii_above]
                        t_soisno[c, ii] = t_soisno[c, ii_above]
                        dz[c, ii] = dz[c, ii_above]
                        aerosol.mss_bcphi_col[c, ii] = aerosol.mss_bcphi_col[c, ii_above]
                        aerosol.mss_bcpho_col[c, ii] = aerosol.mss_bcpho_col[c, ii_above]
                        aerosol.mss_ocphi_col[c, ii] = aerosol.mss_ocphi_col[c, ii_above]
                        aerosol.mss_ocpho_col[c, ii] = aerosol.mss_ocpho_col[c, ii_above]
                        aerosol.mss_dst1_col[c, ii] = aerosol.mss_dst1_col[c, ii_above]
                        aerosol.mss_dst2_col[c, ii] = aerosol.mss_dst2_col[c, ii_above]
                        aerosol.mss_dst3_col[c, ii] = aerosol.mss_dst3_col[c, ii_above]
                        aerosol.mss_dst4_col[c, ii] = aerosol.mss_dst4_col[c, ii_above]
                        snw_rds[c, ii] = snw_rds[c, ii_above]
                    end
                end
                snl[c] = snl[c] + 1
            end
            j += 1
        end

        # Recompute snow depth and total water
        snow_depth[c] = 0.0
        h2osno_total = 0.0
        for j in (snl[c]+1):0
            jj = j + nlevsno
            snow_depth[c] += dz[c, jj]
            h2osno_total += h2osoi_ice[c, jj] + h2osoi_liq[c, jj]
        end

        # Check if all snow gone
        if snow_depth[c] > 0.0
            l = col_landunit[c]
            if (lun_itype[c] == ISTDLAK && snow_depth[c] < dzmin[1] + LSADZ) ||
               (lun_itype[c] != ISTDLAK &&
                (frac_sno_eff[c] * snow_depth[c] < dzmin[1] ||
                 h2osno_total / (frac_sno_eff[c] * snow_depth[c]) < 50.0))

                # Transfer ice to h2osno_no_layers, liquid to top soil layer
                zwice = 0.0
                zwliq = 0.0
                for j in (snl[c]+1):0
                    jj = j + nlevsno
                    zwice += h2osoi_ice[c, jj]
                    zwliq += h2osoi_liq[c, jj]
                end

                h2osno_no_layers[c] = zwice
                if lun_itype[c] == ISTSOIL || urbpoi[l] || lun_itype[c] == ISTCROP
                    # Transfer liquid to soil layer 1
                    jj_soil1 = 0 + nlevsno + 1  # layer index 1
                    if jj_soil1 <= size(h2osoi_liq, 2)
                        h2osoi_liq[c, jj_soil1] += zwliq
                    end
                end

                snl[c] = 0
                h2osno_total = h2osno_no_layers[c]

                # Zero out aerosol masses
                for jj in 1:nlevsno
                    aerosol.mss_bcphi_col[c, jj] = 0.0
                    aerosol.mss_bcpho_col[c, jj] = 0.0
                    aerosol.mss_ocphi_col[c, jj] = 0.0
                    aerosol.mss_ocpho_col[c, jj] = 0.0
                    aerosol.mss_dst1_col[c, jj] = 0.0
                    aerosol.mss_dst2_col[c, jj] = 0.0
                    aerosol.mss_dst3_col[c, jj] = 0.0
                    aerosol.mss_dst4_col[c, jj] = 0.0
                end

                if h2osno_no_layers[c] <= 0.0
                    snow_depth[c] = 0.0
                end
            end
        end

        if h2osno_total <= 0.0
            snow_depth[c] = 0.0
            frac_sno[c] = 0.0
            frac_sno_eff[c] = 0.0
            int_snow[c] = 0.0
        end

        # Combine thin layers (two or more layers)
        if snl[c] < -1
            mssi = 1
            i = snl[c] + 1
            while i <= 0
                jj_i = i + nlevsno
                if (frac_sno_eff[c] * dz[c, jj_i] < dzminloc[mssi]) ||
                   ((h2osoi_ice[c, jj_i] + h2osoi_liq[c, jj_i]) /
                    (frac_sno_eff[c] * dz[c, jj_i]) < 50.0)

                    if i == snl[c] + 1
                        neibor = i + 1
                    elseif i == 0
                        neibor = i - 1
                    else
                        neibor = i + 1
                        jj_im1 = (i - 1) + nlevsno
                        jj_ip1 = (i + 1) + nlevsno
                        if (dz[c, jj_im1] + dz[c, jj_i]) < (dz[c, jj_ip1] + dz[c, jj_i])
                            neibor = i - 1
                        end
                    end

                    # Node l and j are combined, stored as node j
                    if neibor > i
                        j_keep = neibor
                        l_remove = i
                    else
                        j_keep = i
                        l_remove = neibor
                    end

                    jj_keep = j_keep + nlevsno
                    jj_remove = l_remove + nlevsno

                    # Combine aerosol masses
                    aerosol.mss_bcphi_col[c, jj_keep] += aerosol.mss_bcphi_col[c, jj_remove]
                    aerosol.mss_bcpho_col[c, jj_keep] += aerosol.mss_bcpho_col[c, jj_remove]
                    aerosol.mss_ocphi_col[c, jj_keep] += aerosol.mss_ocphi_col[c, jj_remove]
                    aerosol.mss_ocpho_col[c, jj_keep] += aerosol.mss_ocpho_col[c, jj_remove]
                    aerosol.mss_dst1_col[c, jj_keep]  += aerosol.mss_dst1_col[c, jj_remove]
                    aerosol.mss_dst2_col[c, jj_keep]  += aerosol.mss_dst2_col[c, jj_remove]
                    aerosol.mss_dst3_col[c, jj_keep]  += aerosol.mss_dst3_col[c, jj_remove]
                    aerosol.mss_dst4_col[c, jj_keep]  += aerosol.mss_dst4_col[c, jj_remove]

                    # Mass-weighted snow grain size
                    total_mass = h2osoi_liq[c, jj_keep] + h2osoi_ice[c, jj_keep] +
                                 h2osoi_liq[c, jj_remove] + h2osoi_ice[c, jj_remove]
                    if total_mass > 0.0
                        snw_rds[c, jj_keep] = (snw_rds[c, jj_keep] *
                            (h2osoi_liq[c, jj_keep] + h2osoi_ice[c, jj_keep]) +
                            snw_rds[c, jj_remove] *
                            (h2osoi_liq[c, jj_remove] + h2osoi_ice[c, jj_remove])) / total_mass
                    end

                    # Combine elements using combo
                    combo!(dz, c, jj_keep, h2osoi_liq, h2osoi_ice, t_soisno,
                           dz[c, jj_remove], h2osoi_liq[c, jj_remove],
                           h2osoi_ice[c, jj_remove], t_soisno[c, jj_remove])

                    # Shift elements above
                    if j_keep - 1 > snl[c] + 1
                        for k in (j_keep-1):-1:(snl[c]+2)
                            kk = k + nlevsno
                            kk_above = (k - 1) + nlevsno
                            h2osoi_ice[c, kk] = h2osoi_ice[c, kk_above]
                            h2osoi_liq[c, kk] = h2osoi_liq[c, kk_above]
                            t_soisno[c, kk] = t_soisno[c, kk_above]
                            dz[c, kk] = dz[c, kk_above]
                            aerosol.mss_bcphi_col[c, kk] = aerosol.mss_bcphi_col[c, kk_above]
                            aerosol.mss_bcpho_col[c, kk] = aerosol.mss_bcpho_col[c, kk_above]
                            aerosol.mss_ocphi_col[c, kk] = aerosol.mss_ocphi_col[c, kk_above]
                            aerosol.mss_ocpho_col[c, kk] = aerosol.mss_ocpho_col[c, kk_above]
                            aerosol.mss_dst1_col[c, kk] = aerosol.mss_dst1_col[c, kk_above]
                            aerosol.mss_dst2_col[c, kk] = aerosol.mss_dst2_col[c, kk_above]
                            aerosol.mss_dst3_col[c, kk] = aerosol.mss_dst3_col[c, kk_above]
                            aerosol.mss_dst4_col[c, kk] = aerosol.mss_dst4_col[c, kk_above]
                            snw_rds[c, kk] = snw_rds[c, kk_above]
                        end
                    end

                    snl[c] = snl[c] + 1
                    if snl[c] >= -1
                        break
                    end
                else
                    mssi = mssi + 1
                end
                i += 1
            end
        end

        # Reset node depth and interface depth
        for j in 0:-1:(-nlevsno+1)
            jj = j + nlevsno
            if j >= snl[c] + 1
                z[c, jj] = zi[c, jj] - 0.5 * dz[c, jj]
                jj_m1 = (j - 1) + nlevsno
                zi[c, jj_m1] = zi[c, jj] - dz[c, jj]
            end
        end
    end
end

# =========================================================================
# Combo (helper for combining snow layers)
# =========================================================================

"""
    combo!(dz_mat, c, jj, wliq_mat, wice_mat, t_mat,
           dz2, wliq2, wice2, t2)

Combine two snow elements. Updates element at [c, jj] in-place.
"""
function combo!(dz_mat::Matrix{Float64}, c::Int, jj::Int,
                wliq_mat::Matrix{Float64}, wice_mat::Matrix{Float64},
                t_mat::Matrix{Float64},
                dz2::Float64, wliq2::Float64, wice2::Float64, t2::Float64)
    dz1 = dz_mat[c, jj]
    wliq1 = wliq_mat[c, jj]
    wice1 = wice_mat[c, jj]
    t1 = t_mat[c, jj]

    dzc = dz1 + dz2
    wicec = wice1 + wice2
    wliqc = wliq1 + wliq2

    h1 = (CPICE * wice1 + CPLIQ * wliq1) * (t1 - TFRZ) + HFUS * wliq1
    h2 = (CPICE * wice2 + CPLIQ * wliq2) * (t2 - TFRZ) + HFUS * wliq2
    hc = h1 + h2

    denom = CPICE * wicec + CPLIQ * wliqc
    if denom > 0.0
        tc = TFRZ + (hc - HFUS * wliqc) / denom
    else
        tc = TFRZ
    end

    dz_mat[c, jj] = dzc
    wice_mat[c, jj] = wicec
    wliq_mat[c, jj] = wliqc
    t_mat[c, jj] = tc
end

"""
    combo_scalar(dz, wliq, wice, t, dz2, wliq2, wice2, t2) -> (dz, wliq, wice, t)

Scalar version of combo for use in divide_snow_layers.
"""
function combo_scalar(dz::Float64, wliq::Float64, wice::Float64, t::Float64,
                      dz2::Float64, wliq2::Float64, wice2::Float64, t2::Float64)
    dzc = dz + dz2
    wicec = wice + wice2
    wliqc = wliq + wliq2

    h = (CPICE * wice + CPLIQ * wliq) * (t - TFRZ) + HFUS * wliq
    h2val = (CPICE * wice2 + CPLIQ * wliq2) * (t2 - TFRZ) + HFUS * wliq2
    hc = h + h2val

    denom = CPICE * wicec + CPLIQ * wliqc
    if denom > 0.0
        tc = TFRZ + (hc - HFUS * wliqc) / denom
    else
        tc = TFRZ
    end

    return (dzc, wliqc, wicec, tc)
end

# =========================================================================
# MassWeightedSnowRadius
# =========================================================================

"""
    mass_weighted_snow_radius(rds1, rds2, swtot, zwtot) -> Float64

Calculate mass-weighted snow radius when two layers are combined.
"""
function mass_weighted_snow_radius(rds1::Float64, rds2::Float64,
                                    swtot::Float64, zwtot::Float64)
    params = snowhydrology_params
    result = (rds2 * swtot + rds1 * zwtot) / (swtot + zwtot)
    if result > SNW_RDS_MAX
        result = SNW_RDS_MAX
    elseif result < params.snw_rds_min
        result = params.snw_rds_min
    end
    return result
end

# =========================================================================
# DivideSnowLayers
# =========================================================================

"""
    divide_snow_layers!(snl, dz, zi, z, t_soisno,
        h2osoi_ice, h2osoi_liq, frac_sno, snw_rds, aerosol,
        is_lake, mask_snow, bounds, nlevsno)

Subdivide snow layers that exceed their prescribed maximum thickness.
"""
function divide_snow_layers!(
    snl::Vector{Int},
    dz::Matrix{Float64},
    zi::Matrix{Float64},
    z::Matrix{Float64},
    t_soisno::Matrix{Float64},
    h2osoi_ice::Matrix{Float64},
    h2osoi_liq::Matrix{Float64},
    frac_sno::Vector{Float64},
    snw_rds::Matrix{Float64},
    aerosol::AerosolData,
    is_lake::Bool,
    mask_snow::BitVector,
    bounds::UnitRange{Int},
    nlevsno::Int
)
    dzmax_l = SNOW_DZMAX_L[]
    dzmax_u = SNOW_DZMAX_U[]

    for c in bounds
        mask_snow[c] || continue

        msno = abs(snl[c])

        # Copy to local arrays (1-indexed, top=1)
        dzsno  = zeros(Float64, nlevsno)
        swice  = zeros(Float64, nlevsno)
        swliq  = zeros(Float64, nlevsno)
        tsno   = zeros(Float64, nlevsno)
        mbc_phi = zeros(Float64, nlevsno)
        mbc_pho = zeros(Float64, nlevsno)
        moc_phi = zeros(Float64, nlevsno)
        moc_pho = zeros(Float64, nlevsno)
        mdst1  = zeros(Float64, nlevsno)
        mdst2  = zeros(Float64, nlevsno)
        mdst3  = zeros(Float64, nlevsno)
        mdst4  = zeros(Float64, nlevsno)
        rds    = zeros(Float64, nlevsno)

        for k in 1:msno
            jj = (k + snl[c]) + nlevsno  # Fortran j+snl(c) -> Julia index
            if is_lake
                dzsno[k] = dz[c, jj]
            else
                dzsno[k] = frac_sno[c] * dz[c, jj]
            end
            swice[k] = h2osoi_ice[c, jj]
            swliq[k] = h2osoi_liq[c, jj]
            tsno[k]  = t_soisno[c, jj]
            mbc_phi[k] = aerosol.mss_bcphi_col[c, jj]
            mbc_pho[k] = aerosol.mss_bcpho_col[c, jj]
            moc_phi[k] = aerosol.mss_ocphi_col[c, jj]
            moc_pho[k] = aerosol.mss_ocpho_col[c, jj]
            mdst1[k] = aerosol.mss_dst1_col[c, jj]
            mdst2[k] = aerosol.mss_dst2_col[c, jj]
            mdst3[k] = aerosol.mss_dst3_col[c, jj]
            mdst4[k] = aerosol.mss_dst4_col[c, jj]
            rds[k]   = snw_rds[c, jj]
        end

        # Traverse layers top to bottom
        k = 1
        while k <= msno && k < nlevsno
            # Bottom layer case
            if k == msno
                offset = is_lake ? 2.0 * LSADZ : 0.0
                if dzsno[k] > dzmax_l[k] + offset
                    msno += 1
                    dzsno[k] /= 2.0
                    dzsno[k+1] = dzsno[k]
                    swice[k] /= 2.0; swice[k+1] = swice[k]
                    swliq[k] /= 2.0; swliq[k+1] = swliq[k]

                    if k == 1
                        tsno[k+1] = tsno[k]
                    else
                        dtdz = (tsno[k-1] - tsno[k]) / ((dzsno[k-1] + 2*dzsno[k]) / 2.0)
                        tsno[k+1] = tsno[k] - dtdz * dzsno[k] / 2.0
                        if tsno[k+1] >= TFRZ
                            tsno[k+1] = tsno[k]
                        else
                            tsno[k] = tsno[k] + dtdz * dzsno[k] / 2.0
                        end
                    end

                    mbc_phi[k] /= 2.0; mbc_phi[k+1] = mbc_phi[k]
                    mbc_pho[k] /= 2.0; mbc_pho[k+1] = mbc_pho[k]
                    moc_phi[k] /= 2.0; moc_phi[k+1] = moc_phi[k]
                    moc_pho[k] /= 2.0; moc_pho[k+1] = moc_pho[k]
                    mdst1[k] /= 2.0; mdst1[k+1] = mdst1[k]
                    mdst2[k] /= 2.0; mdst2[k+1] = mdst2[k]
                    mdst3[k] /= 2.0; mdst3[k+1] = mdst3[k]
                    mdst4[k] /= 2.0; mdst4[k+1] = mdst4[k]
                    rds[k+1] = rds[k]
                end
            end

            # Layers below exist
            if k < msno
                offset = is_lake ? LSADZ : 0.0
                if dzsno[k] > dzmax_u[k] + offset
                    drr = dzsno[k] - dzmax_u[k] - offset
                    propor = drr / dzsno[k]

                    zwice = propor * swice[k]
                    zwliq = propor * swliq[k]
                    zmbc_phi = propor * mbc_phi[k]
                    zmbc_pho = propor * mbc_pho[k]
                    zmoc_phi = propor * moc_phi[k]
                    zmoc_pho = propor * moc_pho[k]
                    zmdst1 = propor * mdst1[k]
                    zmdst2 = propor * mdst2[k]
                    zmdst3 = propor * mdst3[k]
                    zmdst4 = propor * mdst4[k]

                    propor_keep = (dzmax_u[k] + offset) / dzsno[k]
                    swice[k] *= propor_keep
                    swliq[k] *= propor_keep
                    mbc_phi[k] *= propor_keep
                    mbc_pho[k] *= propor_keep
                    moc_phi[k] *= propor_keep
                    moc_pho[k] *= propor_keep
                    mdst1[k] *= propor_keep
                    mdst2[k] *= propor_keep
                    mdst3[k] *= propor_keep
                    mdst4[k] *= propor_keep

                    dzsno[k] = dzmax_u[k] + offset

                    mbc_phi[k+1] += zmbc_phi
                    mbc_pho[k+1] += zmbc_pho
                    moc_phi[k+1] += zmoc_phi
                    moc_pho[k+1] += zmoc_pho
                    mdst1[k+1] += zmdst1
                    mdst2[k+1] += zmdst2
                    mdst3[k+1] += zmdst3
                    mdst4[k+1] += zmdst4

                    rds[k+1] = mass_weighted_snow_radius(rds[k], rds[k+1],
                        swliq[k+1] + swice[k+1], zwliq + zwice)

                    (dzsno[k+1], swliq[k+1], swice[k+1], tsno[k+1]) =
                        combo_scalar(dzsno[k+1], swliq[k+1], swice[k+1], tsno[k+1],
                                     drr, zwliq, zwice, tsno[k])
                end
            end
            k += 1
        end

        snl[c] = -msno

        # Copy back from local arrays to column arrays
        for j in (-nlevsno+1):0
            jj = j + nlevsno
            if j >= snl[c] + 1
                k_idx = j - snl[c]
                if is_lake
                    dz[c, jj] = dzsno[k_idx]
                else
                    dz[c, jj] = dzsno[k_idx] / frac_sno[c]
                end
                h2osoi_ice[c, jj] = swice[k_idx]
                h2osoi_liq[c, jj] = swliq[k_idx]
                t_soisno[c, jj] = tsno[k_idx]
                aerosol.mss_bcphi_col[c, jj] = mbc_phi[k_idx]
                aerosol.mss_bcpho_col[c, jj] = mbc_pho[k_idx]
                aerosol.mss_ocphi_col[c, jj] = moc_phi[k_idx]
                aerosol.mss_ocpho_col[c, jj] = moc_pho[k_idx]
                aerosol.mss_dst1_col[c, jj] = mdst1[k_idx]
                aerosol.mss_dst2_col[c, jj] = mdst2[k_idx]
                aerosol.mss_dst3_col[c, jj] = mdst3[k_idx]
                aerosol.mss_dst4_col[c, jj] = mdst4[k_idx]
                snw_rds[c, jj] = rds[k_idx]
            end
        end

        # Reset node depth and interface depth
        for j in 0:-1:(-nlevsno+1)
            jj = j + nlevsno
            if j >= snl[c] + 1
                z[c, jj] = zi[c, jj] - 0.5 * dz[c, jj]
                jj_m1 = (j - 1) + nlevsno
                zi[c, jj_m1] = zi[c, jj] - dz[c, jj]
            end
        end
    end
end

# =========================================================================
# ZeroEmptySnowLayers
# =========================================================================

"""
    zero_empty_snow_layers!(snl, dz, z, zi, t_soisno,
        h2osoi_ice, h2osoi_liq, mask_snow, bounds, nlevsno)

Set empty snow layers to zero.
"""
function zero_empty_snow_layers!(
    snl::Vector{Int},
    dz::Matrix{Float64},
    z::Matrix{Float64},
    zi::Matrix{Float64},
    t_soisno::Matrix{Float64},
    h2osoi_ice::Matrix{Float64},
    h2osoi_liq::Matrix{Float64},
    mask_snow::BitVector,
    bounds::UnitRange{Int},
    nlevsno::Int
)
    for j in (-nlevsno+1):0
        jj = j + nlevsno
        for c in bounds
            mask_snow[c] || continue
            if j <= snl[c] && snl[c] > -nlevsno
                h2osoi_ice[c, jj] = 0.0
                h2osoi_liq[c, jj] = 0.0
                t_soisno[c, jj] = 0.0
                dz[c, jj] = 0.0
                z[c, jj] = 0.0
                jj_m1 = (j - 1) + nlevsno
                zi[c, jj_m1] = 0.0
            end
        end
    end
end

# =========================================================================
# InitSnowLayers
# =========================================================================

"""
    init_snow_layers!(snl, dz, z, zi, snow_depth, col_landunit, lakpoi, bounds, nlevsno)

Initialize snow layer depth from specified total depth.
Also initializes the module-level dzmin, dzmax_l, dzmax_u arrays.
"""
function init_snow_layers!(
    snl::Vector{Int},
    dz::Matrix{Float64},
    z::Matrix{Float64},
    zi::Matrix{Float64},
    snow_depth::Vector{Float64},
    col_landunit::Vector{Int},
    lakpoi::Vector{Bool},
    bounds::UnitRange{Int},
    nlevsno::Int
)
    # Initialize snow layer thickness arrays
    dzmin = zeros(Float64, nlevsno)
    dzmax_l_arr = zeros(Float64, nlevsno)
    dzmax_u_arr = zeros(Float64, nlevsno)

    dzmin[1] = SNOW_DZMIN_1[]
    dzmax_l_arr[1] = SNOW_DZMAX_L_1[]
    dzmax_u_arr[1] = SNOW_DZMAX_U_1[]
    dzmin[2] = SNOW_DZMIN_2[]
    dzmax_l_arr[2] = SNOW_DZMAX_L_2[]
    dzmax_u_arr[2] = SNOW_DZMAX_U_2[]

    for j in 3:nlevsno
        dzmin[j] = dzmax_u_arr[j-1] * 0.5
        dzmax_u_arr[j] = 2.0 * dzmax_u_arr[j-1] + 0.01
        dzmax_l_arr[j] = dzmax_u_arr[j] + dzmax_l_arr[j-1]
        if j == nlevsno
            dzmax_u_arr[j] = floatmax(Float64)
            dzmax_l_arr[j] = floatmax(Float64)
        end
    end

    # Store module-level arrays
    SNOW_DZMIN[] = dzmin
    SNOW_DZMAX_L[] = dzmax_l_arr
    SNOW_DZMAX_U[] = dzmax_u_arr

    for c in bounds
        l = col_landunit[c]

        # Initialize snow layers to spval/zero
        for jj in 1:nlevsno
            dz[c, jj] = SPVAL
            z[c, jj] = SPVAL
        end
        for jj in 1:nlevsno
            if jj <= size(zi, 2)
                zi[c, jj] = SPVAL
            end
        end

        # Lake: no snow layers
        if lakpoi[l]
            snl[c] = 0
            for jj in 1:nlevsno
                dz[c, jj] = 0.0
                z[c, jj] = 0.0
            end
            for jj in 1:nlevsno
                if jj <= size(zi, 2)
                    zi[c, jj] = 0.0
                end
            end
            continue
        end

        # Too little snow
        if snow_depth[c] < dzmin[1]
            snl[c] = 0
            for jj in 1:nlevsno
                dz[c, jj] = 0.0
                z[c, jj] = 0.0
            end
            for jj in 1:nlevsno
                if jj <= size(zi, 2)
                    zi[c, jj] = 0.0
                end
            end
            continue
        end

        # At least one layer
        snl[c] = -1
        minbound = dzmin[1]
        maxbound = dzmax_l_arr[1]

        jj_zero = 0 + nlevsno  # Julia index for layer 0

        if snow_depth[c] >= minbound && snow_depth[c] <= maxbound
            dz[c, jj_zero] = snow_depth[c]
        else
            snl[c] = snl[c] - 1
            minbound = maxbound
            maxbound = sum(dzmax_u_arr[1:(-snl[c])])

            while snow_depth[c] > maxbound && (-snl[c]) < nlevsno
                snl[c] = snl[c] - 1
                minbound = maxbound
                maxbound = sum(dzmax_u_arr[1:(-snl[c])])
            end

            # Set thickness of all layers except bottom two
            for j in 1:(-snl[c]-2)
                jj = (j + snl[c]) + nlevsno
                dz[c, jj] = dzmax_u_arr[j]
            end

            # Determine bottom two layers
            if snow_depth[c] <= sum(dzmax_u_arr[1:(-snl[c]-2)]) + 2 * dzmax_u_arr[-snl[c]-1]
                jj_m1 = (-1) + nlevsno
                dz[c, jj_m1] = (snow_depth[c] - sum(dzmax_u_arr[1:(-snl[c]-2)])) / 2.0
                dz[c, jj_zero] = dz[c, jj_m1]
            else
                jj_m1 = (-1) + nlevsno
                dz[c, jj_m1] = dzmax_u_arr[-snl[c]-1]
                dz[c, jj_zero] = snow_depth[c] - sum(dzmax_u_arr[1:(-snl[c]-1)])
            end
        end

        # Initialize node depth and interface depth
        for j in 0:-1:(snl[c]+1)
            jj = j + nlevsno
            z[c, jj] = zi[c, jj] - 0.5 * dz[c, jj]
            jj_m1 = (j - 1) + nlevsno
            zi[c, jj_m1] = zi[c, jj] - dz[c, jj]
        end
    end
end

# =========================================================================
# BuildSnowFilter
# =========================================================================

"""
    build_snow_filter!(mask_snow, mask_nosnow, snl, mask_nolake, bounds)

Construct snow/no-snow masks for snow hydrology.
"""
function build_snow_filter!(
    mask_snow::BitVector,
    mask_nosnow::BitVector,
    snl::Vector{Int},
    mask_nolake::BitVector,
    bounds::UnitRange{Int}
)
    fill!(mask_snow, false)
    fill!(mask_nosnow, false)
    for c in bounds
        mask_nolake[c] || continue
        if snl[c] < 0
            mask_snow[c] = true
        else
            mask_nosnow[c] = true
        end
    end
end

# =========================================================================
# SnowCappingExcess
# =========================================================================

"""
    snow_capping_excess!(h2osno_excess, apply_runoff,
        h2osno, topo, snl, col_landunit, lun_itype,
        mask_snow, bounds, nstep)

Determine the excess snow that needs to be capped.
"""
function snow_capping_excess!(
    h2osno_excess::Vector{Float64},
    apply_runoff::Vector{Bool},
    h2osno::Vector{Float64},
    topo::Vector{Float64},
    snl::Vector{Int},
    col_landunit::Vector{Int},
    lun_itype::Vector{Int},   # landunit type
    mask_snow::BitVector,
    bounds::UnitRange{Int},
    nstep::Int,
    nlevsno::Int
)
    h2osno_max_val = 10000.0  # Maximum allowed SWE [mm H2O]

    for c in bounds
        mask_snow[c] || continue
        h2osno_excess[c] = 0.0
        if h2osno[c] > h2osno_max_val
            h2osno_excess[c] = h2osno[c] - h2osno_max_val
            apply_runoff[c] = true
        end
    end

    # Snow resetting
    is_reset_snow_active = false
    if RESET_SNOW[] || RESET_SNOW_GLC[]
        reset_snow_timesteps = RESET_SNOW_TIMESTEPS_PER_LAYER * nlevsno
        if nstep <= reset_snow_timesteps
            is_reset_snow_active = true
        end
    end

    if is_reset_snow_active
        for c in bounds
            mask_snow[c] || continue
            l = col_landunit[c]
            if lun_itype[l] != ISTICE && RESET_SNOW[] && h2osno[c] > RESET_SNOW_H2OSNO
                h2osno_excess[c] = h2osno[c] - RESET_SNOW_H2OSNO
                apply_runoff[c] = false
            elseif lun_itype[l] == ISTICE && RESET_SNOW_GLC[] &&
                   h2osno[c] > RESET_SNOW_H2OSNO && topo[c] <= RESET_SNOW_GLC_ELA[]
                h2osno_excess[c] = h2osno[c] - RESET_SNOW_H2OSNO
                apply_runoff[c] = false
            end
        end
    end
end

# =========================================================================
# InitFlux_SnowCapping
# =========================================================================

"""
    init_flux_snow_capping!(qflx_snwcp_ice, qflx_snwcp_liq,
        qflx_snwcp_discarded_ice, qflx_snwcp_discarded_liq,
        mask, bounds)

Initialize snow capping fluxes to zero.
"""
function init_flux_snow_capping!(
    qflx_snwcp_ice::Vector{Float64},
    qflx_snwcp_liq::Vector{Float64},
    qflx_snwcp_discarded_ice::Vector{Float64},
    qflx_snwcp_discarded_liq::Vector{Float64},
    mask::BitVector,
    bounds::UnitRange{Int}
)
    for c in bounds
        mask[c] || continue
        qflx_snwcp_ice[c] = 0.0
        qflx_snwcp_liq[c] = 0.0
        qflx_snwcp_discarded_ice[c] = 0.0
        qflx_snwcp_discarded_liq[c] = 0.0
    end
end

# =========================================================================
# BulkFlux_SnowCappingFluxes
# =========================================================================

"""
    bulk_flux_snow_capping_fluxes!(mask_capping, rho_orig_bottom, frac_adjust,
        qflx_snwcp_ice, qflx_snwcp_liq,
        qflx_snwcp_discarded_ice, qflx_snwcp_discarded_liq,
        dtime, dz_bottom_jj, topo, h2osno_total,
        h2osoi_ice_bottom, h2osoi_liq_bottom,
        col_landunit, lun_itype,
        mask_snow, bounds, nlevsno, nstep)

Calculate snow capping fluxes and related terms for bulk water.
"""
function bulk_flux_snow_capping_fluxes!(
    mask_capping::BitVector,
    rho_orig_bottom::Vector{Float64},
    frac_adjust::Vector{Float64},
    qflx_snwcp_ice::Vector{Float64},
    qflx_snwcp_liq::Vector{Float64},
    qflx_snwcp_discarded_ice::Vector{Float64},
    qflx_snwcp_discarded_liq::Vector{Float64},
    dtime::Float64,
    dz_bottom::Vector{Float64},
    topo::Vector{Float64},
    h2osno_total::Vector{Float64},
    h2osoi_ice_bottom::Vector{Float64},
    h2osoi_liq_bottom::Vector{Float64},
    col_landunit::Vector{Int},
    lun_itype::Vector{Int},
    mask_snow::BitVector,
    bounds::UnitRange{Int},
    nlevsno::Int,
    nstep::Int
)
    h2osno_max_val = 10000.0
    min_snow_to_keep = 1.0e-3

    nc = length(h2osno_total)
    h2osno_excess = zeros(Float64, nc)
    apply_runoff = falses(nc)

    snow_capping_excess!(h2osno_excess, apply_runoff,
        h2osno_total, topo, zeros(Int, nc), col_landunit, lun_itype,
        mask_snow, bounds, nstep, nlevsno)

    fill!(mask_capping, false)
    for c in bounds
        mask_snow[c] || continue
        if h2osno_excess[c] > 0.0
            mask_capping[c] = true
        end
    end

    for c in bounds
        mask_capping[c] || continue

        rho_orig_bottom[c] = h2osoi_ice_bottom[c] / dz_bottom[c]

        mss_snow_bottom_lyr = h2osoi_ice_bottom[c] + h2osoi_liq_bottom[c]
        mss_snwcp_tot = min(h2osno_excess[c],
            mss_snow_bottom_lyr * (1.0 - min_snow_to_keep))

        icefrac = h2osoi_ice_bottom[c] / mss_snow_bottom_lyr
        snwcp_flux_ice = mss_snwcp_tot / dtime * icefrac
        snwcp_flux_liq = mss_snwcp_tot / dtime * (1.0 - icefrac)

        if apply_runoff[c]
            qflx_snwcp_ice[c] = snwcp_flux_ice
            qflx_snwcp_liq[c] = snwcp_flux_liq
        else
            qflx_snwcp_discarded_ice[c] = snwcp_flux_ice
            qflx_snwcp_discarded_liq[c] = snwcp_flux_liq
        end

        frac_adjust[c] = (mss_snow_bottom_lyr - mss_snwcp_tot) / mss_snow_bottom_lyr
    end
end

# =========================================================================
# UpdateState_RemoveSnowCappingFluxes
# =========================================================================

"""
    update_state_remove_snow_capping_fluxes!(h2osoi_ice_bottom, h2osoi_liq_bottom,
        dtime, qflx_snwcp_ice, qflx_snwcp_liq,
        qflx_snwcp_discarded_ice, qflx_snwcp_discarded_liq,
        mask_capping, bounds)

Remove snow capping fluxes from h2osoi_ice and h2osoi_liq.
"""
function update_state_remove_snow_capping_fluxes!(
    h2osoi_ice_bottom::Vector{Float64},
    h2osoi_liq_bottom::Vector{Float64},
    dtime::Float64,
    qflx_snwcp_ice::Vector{Float64},
    qflx_snwcp_liq::Vector{Float64},
    qflx_snwcp_discarded_ice::Vector{Float64},
    qflx_snwcp_discarded_liq::Vector{Float64},
    mask_capping::BitVector,
    bounds::UnitRange{Int}
)
    for c in bounds
        mask_capping[c] || continue
        h2osoi_ice_bottom[c] -= (qflx_snwcp_ice[c] + qflx_snwcp_discarded_ice[c]) * dtime
        h2osoi_liq_bottom[c] -= (qflx_snwcp_liq[c] + qflx_snwcp_discarded_liq[c]) * dtime

        if h2osoi_ice_bottom[c] < 0.0 || h2osoi_liq_bottom[c] < 0.0
            error("snow capping failed: negative mass at c=$c, ice=$(h2osoi_ice_bottom[c]), liq=$(h2osoi_liq_bottom[c])")
        end
    end
end

# =========================================================================
# SnowCappingUpdateDzAndAerosols
# =========================================================================

"""
    snow_capping_update_dz_and_aerosols!(dz_bottom, aerosol, jj_bottom,
        rho_orig_bottom, h2osoi_ice_bottom, frac_adjust,
        mask_capping, bounds)

After snow capping, adjust dz and aerosol masses in bottom snow layer.
"""
function snow_capping_update_dz_and_aerosols!(
    dz_bottom::Vector{Float64},
    aerosol::AerosolData,
    jj_bottom::Int,
    rho_orig_bottom::Vector{Float64},
    h2osoi_ice_bottom::Vector{Float64},
    frac_adjust::Vector{Float64},
    mask_capping::BitVector,
    bounds::UnitRange{Int}
)
    for c in bounds
        mask_capping[c] || continue

        if rho_orig_bottom[c] > 1.0
            dz_bottom[c] = h2osoi_ice_bottom[c] / rho_orig_bottom[c]
        end

        aerosol.mss_bcphi_col[c, jj_bottom] *= frac_adjust[c]
        aerosol.mss_bcpho_col[c, jj_bottom] *= frac_adjust[c]
        aerosol.mss_ocphi_col[c, jj_bottom] *= frac_adjust[c]
        aerosol.mss_ocpho_col[c, jj_bottom] *= frac_adjust[c]
        aerosol.mss_dst1_col[c, jj_bottom]  *= frac_adjust[c]
        aerosol.mss_dst2_col[c, jj_bottom]  *= frac_adjust[c]
        aerosol.mss_dst3_col[c, jj_bottom]  *= frac_adjust[c]
        aerosol.mss_dst4_col[c, jj_bottom]  *= frac_adjust[c]
    end
end

# =========================================================================
# SnowHydrologySetControlForTesting
# =========================================================================

"""
    snow_hydrology_set_control_for_testing!(;
        wind_dep_snow_density=nothing, new_snow_density_method=nothing,
        reset_snow_flag=nothing, reset_snow_glc_flag=nothing,
        reset_snow_glc_ela_val=nothing)

Set some of the control settings for SnowHydrology.
NOTE: This is just for unit testing.
"""
function snow_hydrology_set_control_for_testing!(;
    wind_dep_snow_density::Union{Nothing,Bool} = nothing,
    new_snow_density_method::Union{Nothing,Int} = nothing,
    reset_snow_flag::Union{Nothing,Bool} = nothing,
    reset_snow_glc_flag::Union{Nothing,Bool} = nothing,
    reset_snow_glc_ela_val::Union{Nothing,Float64} = nothing
)
    if !isnothing(wind_dep_snow_density)
        WIND_DEPENDENT_SNOW_DENSITY[] = wind_dep_snow_density
    end
    if !isnothing(new_snow_density_method)
        NEW_SNOW_DENSITY[] = new_snow_density_method
    end
    if !isnothing(reset_snow_flag)
        RESET_SNOW[] = reset_snow_flag
    end
    if !isnothing(reset_snow_glc_flag)
        RESET_SNOW_GLC[] = reset_snow_glc_flag
    end
    if !isnothing(reset_snow_glc_ela_val)
        RESET_SNOW_GLC_ELA[] = reset_snow_glc_ela_val
    end

    # Set default namelist values
    SNOW_DZMIN_1[] = 0.010
    SNOW_DZMIN_2[] = 0.015
    SNOW_DZMAX_L_1[] = 0.03
    SNOW_DZMAX_L_2[] = 0.07
    SNOW_DZMAX_U_1[] = 0.02
    SNOW_DZMAX_U_2[] = 0.05

    snowhydrology_params.wind_snowcompact_fact = 5.0
end
