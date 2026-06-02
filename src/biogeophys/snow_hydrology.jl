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
    rho_min::Float64 = 50.0                   # minimum snow density for compaction (kg/m3)
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
    bifall::Vector{<:Real},        # output: bulk density [kg/m3]
    forc_t::Vector{<:Real},        # input: atmospheric temperature [K]
    forc_wind::Vector{<:Real},     # input: atmospheric wind speed [m/s] (gridcell)
    col_gridcell::Vector{Int},      # input: column-to-gridcell mapping
    mask::BitVector,                # input: column mask
    bounds::UnitRange{Int}          # input: column bounds
)
    params = snowhydrology_params

    for c in bounds
        mask[c] || continue
        g = col_gridcell[c]

        # Smooth TFRZ branches for AD compatibility
        # Weight for T > TFRZ + 2.0
        w_high = smooth_heaviside(forc_t[c] - (TFRZ + 2.0))
        # Weight for TFRZ - 15.0 < T <= TFRZ + 2.0
        w_mid = smooth_heaviside(forc_t[c] - (TFRZ - 15.0)) * (1.0 - w_high)
        # Weight for T <= TFRZ - 15.0
        w_low = 1.0 - w_high - w_mid

        bifall_high = 50.0 + 1.7 * (17.0)^1.5
        bifall_mid = 50.0 + 1.7 * smooth_max(forc_t[c] - TFRZ + 15.0, 0.0)^1.5

        if NEW_SNOW_DENSITY[] == LO_TMP_DNS_TRUNCATED_ANDERSON1976
            bifall_low = 50.0
        else  # LO_TMP_DNS_SLATER2017
            t_for_bifall_degC = smooth_max(forc_t[c] - TFRZ, -57.55)
            bifall_low = -(50.0/15.0 + 0.0333*15.0)*t_for_bifall_degC - 0.0333*t_for_bifall_degC^2
        end

        bifall[c] = w_high * bifall_high + w_mid * bifall_mid + w_low * bifall_low

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
    dz::Matrix{<:Real},            # layer depth [m]
    int_snow::Vector{<:Real},      # integrated snowfall [mm H2O]
    swe_old::Matrix{<:Real},       # SWE before update [mm H2O]
    frac_sno::Vector{<:Real},      # fraction of ground covered by snow
    frac_sno_eff::Vector{<:Real},  # effective fraction
    snow_depth::Vector{<:Real},    # snow height [m]
    snomelt_accum::Vector{<:Real}, # accumulated snow melt for z0m [m H2O]
    # Inputs
    dtime::Real,                 # time step [s]
    lun_itype_col::Vector{Int},     # landunit type per column
    urbpoi::Vector{Bool},           # urban point flags (landunit-level, indexed by col)
    snl::Vector{Int},               # number of snow layers (negative)
    bifall::Vector{<:Real},        # bulk density of new snow [kg/m3]
    h2osno_total::Vector{<:Real},  # total snow water [mm H2O]
    h2osoi_ice::Matrix{<:Real},    # ice lens [kg/m2]
    h2osoi_liq::Matrix{<:Real},    # liquid water [kg/m2]
    qflx_snow_grnd::Vector{<:Real}, # snow on ground [mm H2O/s]
    qflx_snow_drain::Vector{<:Real}, # drainage from snow pack [mm H2O/s]
    mask::BitVector,                # column mask
    bounds::UnitRange{Int},         # column bounds
    nlevsno::Int;                   # max snow layers
    scf_method::SnowCoverFractionBase = SnowCoverFractionSwensonLawrence2012()
)
    FT = eltype(snow_depth)
    newsnow = zeros(FT, length(snow_depth))
    snowmelt = zeros(FT, length(snow_depth))
    temp_snow_depth = zeros(FT, length(snow_depth))

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
        snomelt_accum[c] = smooth_max(0.0, snomelt_accum[c] - newsnow[c] * 1.0e-3)

        # Update int_snow
        int_snow[c] = smooth_max(int_snow[c], h2osno_total[c])

        # Snowmelt from previous time step
        snowmelt[c] = qflx_snow_drain[c] * dtime
    end

    # Update snow depth and fractional snow cover via scf_method
    update_snow_depth_and_frac!(scf_method,
        frac_sno, frac_sno_eff, snow_depth,
        lun_itype_col, urbpoi, h2osno_total, snowmelt, int_snow,
        newsnow, bifall, mask, bounds)

    # Reset int_snow if snow pack started at 0
    for c in bounds
        mask[c] || continue
        if h2osno_total[c] == 0.0 && newsnow[c] > 0.0
            int_snow[c] = 0.0
        end
    end

    # Add newsnow to int_snow via scf_method
    add_newsnow_to_intsnow!(scf_method,
        int_snow, newsnow, h2osno_total, frac_sno,
        mask, bounds)

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
@kernel function _snowhyd_add_new_snow_kernel!(h2osno_no_layers, @Const(mask),
        @Const(snl), @Const(qflx_snow_grnd), h2osoi_ice, dtime, nlevsno::Int,
        cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && mask[c]
        if snl[c] == 0
            h2osno_no_layers[c] = h2osno_no_layers[c] + qflx_snow_grnd[c] * dtime
        else
            jj = (snl[c] + 1) + nlevsno  # top snow layer Julia index
            h2osoi_ice[c, jj] = h2osoi_ice[c, jj] + qflx_snow_grnd[c] * dtime
        end
    end
end

snowhyd_add_new_snow!(h2osno_no_layers, mask, snl, qflx_snow_grnd, h2osoi_ice, dtime, nlevsno, bounds) =
    _launch!(_snowhyd_add_new_snow_kernel!, h2osno_no_layers, mask, snl, qflx_snow_grnd,
             h2osoi_ice, dtime, nlevsno, first(bounds), last(bounds))

function update_state_add_new_snow!(
    h2osno_no_layers::Vector{<:Real},
    h2osoi_ice::Matrix{<:Real},
    dtime::Real,
    snl::Vector{Int},
    qflx_snow_grnd::Vector{<:Real},
    mask::BitVector,
    bounds::UnitRange{Int},
    nlevsno::Int
)
    snowhyd_add_new_snow!(h2osno_no_layers, mask, snl, qflx_snow_grnd,
                          h2osoi_ice, dtime, nlevsno, bounds)
end

# =========================================================================
# BuildFilter_ThawedWetlandThinSnowpack
# =========================================================================

"""
    build_filter_thawed_wetland_thin_snowpack!(mask_out,
        t_grnd, lun_itype_col, snl, mask_nolake, bounds)

Build a column-level mask of thawed wetland columns with thin (no-layer) snow pack.
"""
@kernel function _snowhyd_thawed_wetland_filter_kernel!(mask_out, @Const(t_grnd),
        @Const(lun_itype_col), @Const(snl), @Const(mask_nolake), cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && mask_nolake[c]
        mask_out[c] = (lun_itype_col[c] == ISTWET && t_grnd[c] > TFRZ && snl[c] == 0)
    end
end

snowhyd_thawed_wetland_filter!(mask_out, t_grnd, lun_itype_col, snl, mask_nolake, bounds) =
    _launch!(_snowhyd_thawed_wetland_filter_kernel!, mask_out, t_grnd, lun_itype_col,
             snl, mask_nolake, first(bounds), last(bounds))

function build_filter_thawed_wetland_thin_snowpack!(
    mask_out::BitVector,
    t_grnd::Vector{<:Real},
    lun_itype_col::Vector{Int},
    snl::Vector{Int},
    mask_nolake::BitVector,
    bounds::UnitRange{Int}
)
    fill!(mask_out, false)
    snowhyd_thawed_wetland_filter!(mask_out, t_grnd, lun_itype_col, snl, mask_nolake, bounds)
end

# =========================================================================
# UpdateState_RemoveSnowFromThawedWetlands
# =========================================================================

"""
    update_state_remove_snow_thawed_wetlands!(h2osno_no_layers, mask, bounds)

Remove snow from thawed wetlands — set h2osno_no_layers to zero.
"""
@kernel function _snowhyd_remove_snow_wetlands_kernel!(h2osno_no_layers, @Const(mask),
        cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && mask[c]
        h2osno_no_layers[c] = 0.0
    end
end

snowhyd_remove_snow_wetlands!(h2osno_no_layers, mask, bounds) =
    _launch!(_snowhyd_remove_snow_wetlands_kernel!, h2osno_no_layers, mask,
             first(bounds), last(bounds))

function update_state_remove_snow_thawed_wetlands!(
    h2osno_no_layers::Vector{<:Real},
    mask::BitVector,
    bounds::UnitRange{Int}
)
    snowhyd_remove_snow_wetlands!(h2osno_no_layers, mask, bounds)
end

# =========================================================================
# Bulk_RemoveSnowFromThawedWetlands
# =========================================================================

"""
    bulk_remove_snow_thawed_wetlands!(snow_depth, mask, bounds)

Remove snow from thawed wetlands — set snow_depth to zero.
"""
@kernel function _snowhyd_bulk_remove_snow_wetlands_kernel!(snow_depth, @Const(mask),
        cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && mask[c]
        snow_depth[c] = 0.0
    end
end

snowhyd_bulk_remove_snow_wetlands!(snow_depth, mask, bounds) =
    _launch!(_snowhyd_bulk_remove_snow_wetlands_kernel!, snow_depth, mask,
             first(bounds), last(bounds))

function bulk_remove_snow_thawed_wetlands!(
    snow_depth::Vector{<:Real},
    mask::BitVector,
    bounds::UnitRange{Int}
)
    snowhyd_bulk_remove_snow_wetlands!(snow_depth, mask, bounds)
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
@kernel function _snowhyd_snowpack_init_filter_kernel!(mask_out, @Const(snl),
        @Const(lun_itype_col), @Const(frac_sno_eff), @Const(snow_depth),
        @Const(mask), dzmin1, cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && mask[c]
        if lun_itype_col[c] == ISTDLAK
            mask_out[c] = (snl[c] == 0 &&
                           frac_sno_eff[c] * snow_depth[c] >= (dzmin1 + LSADZ))
        else
            mask_out[c] = (snl[c] == 0 &&
                           frac_sno_eff[c] * snow_depth[c] >= dzmin1)
        end
    end
end

snowhyd_snowpack_init_filter!(mask_out, snl, lun_itype_col, frac_sno_eff, snow_depth, mask, dzmin1, bounds) =
    _launch!(_snowhyd_snowpack_init_filter_kernel!, mask_out, snl, lun_itype_col,
             frac_sno_eff, snow_depth, mask, dzmin1, first(bounds), last(bounds))

function build_filter_snowpack_initialized!(
    mask_out::BitVector,
    snl::Vector{Int},
    lun_itype_col::Vector{Int},
    frac_sno_eff::Vector{<:Real},
    snow_depth::Vector{<:Real},
    qflx_snow_grnd::Vector{<:Real},
    mask::BitVector,
    bounds::UnitRange{Int}
)
    dzmin = SNOW_DZMIN[]

    fill!(mask_out, false)
    snowhyd_snowpack_init_filter!(mask_out, snl, lun_itype_col, frac_sno_eff,
                                  snow_depth, mask, dzmin[1], bounds)
end

# =========================================================================
# UpdateState_InitializeSnowPack
# =========================================================================

"""
    update_state_initialize_snow_pack!(h2osno_no_layers, h2osoi_ice, h2osoi_liq,
        mask, bounds, nlevsno)

Initialize water state variables for columns with newly-initialized snow pack.
"""
@kernel function _snowhyd_initialize_snow_pack_kernel!(h2osno_no_layers, @Const(mask),
        h2osoi_ice, h2osoi_liq, jj_zero::Int, cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && mask[c]
        h2osoi_ice[c, jj_zero] = h2osno_no_layers[c]
        h2osoi_liq[c, jj_zero] = 0.0
        h2osno_no_layers[c] = 0.0
    end
end

snowhyd_initialize_snow_pack!(h2osno_no_layers, mask, h2osoi_ice, h2osoi_liq, jj_zero, bounds) =
    _launch!(_snowhyd_initialize_snow_pack_kernel!, h2osno_no_layers, mask,
             h2osoi_ice, h2osoi_liq, jj_zero, first(bounds), last(bounds))

function update_state_initialize_snow_pack!(
    h2osno_no_layers::Vector{<:Real},
    h2osoi_ice::Matrix{<:Real},
    h2osoi_liq::Matrix{<:Real},
    mask::BitVector,
    bounds::UnitRange{Int},
    nlevsno::Int
)
    jj_zero = 0 + nlevsno  # Julia index for layer 0
    snowhyd_initialize_snow_pack!(h2osno_no_layers, mask, h2osoi_ice, h2osoi_liq,
                                  jj_zero, bounds)
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
    zi::Matrix{<:Real},
    dz::Matrix{<:Real},
    z::Matrix{<:Real},
    t_soisno::Matrix{<:Real},
    frac_iceold::Matrix{<:Real},
    snomelt_accum::Vector{<:Real},
    forc_t::Vector{<:Real},
    snow_depth::Vector{<:Real},
    mask::BitVector,
    bounds::UnitRange{Int},
    nlevsno::Int
)
    jj_zero = 0 + nlevsno  # Julia index for layer 0 (dz/z arrays)
    jj_surface_zi = nlevsno + 1  # Julia index for interface j=0 (surface)

    snowhyd_initialize_snow_pack_bulk!(snl, mask, zi, dz, z, t_soisno, frac_iceold,
        snomelt_accum, forc_t, snow_depth, jj_zero, jj_surface_zi, bounds)
end

@kernel function _snowhyd_initialize_snow_pack_bulk_kernel!(snl, @Const(mask),
        zi, dz, z, t_soisno, frac_iceold, snomelt_accum, @Const(forc_t),
        @Const(snow_depth), jj_zero::Int, jj_surface_zi::Int, cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && mask[c]
        snl[c] = -1
        dz[c, jj_zero] = snow_depth[c]
        # Keep snow-layer interfaces anchored to the fixed surface interface.
        zi[c, jj_surface_zi] = 0.0
        # zi has one extra interface level: Fortran j -> Julia j + nlevsno + 1
        # j=-1 (below layer 0) maps to jj_zero; j=0 maps to jj_surface_zi.
        zi[c, jj_zero] = zi[c, jj_surface_zi] - dz[c, jj_zero]
        z[c, jj_zero] = zi[c, jj_surface_zi] - 0.5 * dz[c, jj_zero]

        t_soisno[c, jj_zero] = smooth_min(TFRZ, forc_t[c])
        frac_iceold[c, jj_zero] = 1.0
        snomelt_accum[c] = 0.0
    end
end

snowhyd_initialize_snow_pack_bulk!(snl, mask, zi, dz, z, t_soisno, frac_iceold,
        snomelt_accum, forc_t, snow_depth, jj_zero, jj_surface_zi, bounds) =
    _launch!(_snowhyd_initialize_snow_pack_bulk_kernel!, snl, mask, zi, dz, z,
             t_soisno, frac_iceold, snomelt_accum, forc_t, snow_depth,
             jj_zero, jj_surface_zi, first(bounds), last(bounds))

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
    h2osoi_ice::Matrix{<:Real},
    h2osoi_liq::Matrix{<:Real},
    dtime::Real,
    snl::Vector{Int},
    frac_sno_eff::Vector{<:Real},
    qflx_soliddew_to_top_layer::Vector{<:Real},
    qflx_solidevap_from_top_layer::Vector{<:Real},
    qflx_liq_grnd::Vector{<:Real},
    qflx_liqdew_to_top_layer::Vector{<:Real},
    qflx_liqevap_from_top_layer::Vector{<:Real},
    mask_snow::BitVector,
    bounds::UnitRange{Int},
    nlevsno::Int
)
    snowhyd_top_layer_fluxes!(h2osoi_ice, h2osoi_liq, mask_snow, snl, frac_sno_eff,
        qflx_soliddew_to_top_layer, qflx_solidevap_from_top_layer, qflx_liq_grnd,
        qflx_liqdew_to_top_layer, qflx_liqevap_from_top_layer, dtime, nlevsno, bounds)
end

@kernel function _snowhyd_top_layer_fluxes_kernel!(h2osoi_ice, @Const(mask_snow),
        @Const(snl), @Const(frac_sno_eff), @Const(qflx_soliddew_to_top_layer),
        @Const(qflx_solidevap_from_top_layer), @Const(qflx_liq_grnd),
        @Const(qflx_liqdew_to_top_layer), @Const(qflx_liqevap_from_top_layer),
        h2osoi_liq, dtime, nlevsno::Int, cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && mask_snow[c]
        lev_top = snl[c] + 1
        jj = lev_top + nlevsno  # Julia index for top snow layer

        h2osoi_ice[c, jj] = h2osoi_ice[c, jj] +
            frac_sno_eff[c] * (qflx_soliddew_to_top_layer[c] -
            qflx_solidevap_from_top_layer[c]) * dtime

        h2osoi_liq[c, jj] = h2osoi_liq[c, jj] +
            frac_sno_eff[c] * (qflx_liq_grnd[c] + qflx_liqdew_to_top_layer[c] -
            qflx_liqevap_from_top_layer[c]) * dtime

        # Clamp to non-negative (physical constraint, not a smooth approximation)
        h2osoi_ice[c, jj] = max(0.0, h2osoi_ice[c, jj])
        h2osoi_liq[c, jj] = max(0.0, h2osoi_liq[c, jj])
    end
end

snowhyd_top_layer_fluxes!(h2osoi_ice, h2osoi_liq, mask_snow, snl, frac_sno_eff,
        qflx_soliddew_to_top_layer, qflx_solidevap_from_top_layer, qflx_liq_grnd,
        qflx_liqdew_to_top_layer, qflx_liqevap_from_top_layer, dtime, nlevsno, bounds) =
    _launch!(_snowhyd_top_layer_fluxes_kernel!, h2osoi_ice, mask_snow, snl, frac_sno_eff,
             qflx_soliddew_to_top_layer, qflx_solidevap_from_top_layer, qflx_liq_grnd,
             qflx_liqdew_to_top_layer, qflx_liqevap_from_top_layer, h2osoi_liq,
             dtime, nlevsno, first(bounds), last(bounds))

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
    qflx_snow_percolation::Matrix{<:Real},
    dtime::Real,
    snl::Vector{Int},
    dz::Matrix{<:Real},
    frac_sno_eff::Vector{<:Real},
    h2osoi_ice::Matrix{<:Real},
    h2osoi_liq::Matrix{<:Real},
    mask_snow::BitVector,
    bounds::UnitRange{Int},
    nlevsno::Int
)
    params = snowhydrology_params

    # Temporary arrays for porosity/volume calculations
    nc = length(snl)
    FT = eltype(h2osoi_liq)
    vol_liq = zeros(FT, nc, nlevsno)
    vol_ice = zeros(FT, nc, nlevsno)
    eff_porosity = zeros(FT, nc, nlevsno)

    # Porosity and partial volume
    for j in (-nlevsno+1):0
        jj = j + nlevsno
        for c in bounds
            mask_snow[c] || continue
            if j >= snl[c] + 1
                vol_ice[c, jj] = smooth_min(1.0, h2osoi_ice[c, jj] / (dz[c, jj] * frac_sno_eff[c] * DENICE))
                eff_porosity[c, jj] = 1.0 - vol_ice[c, jj]
                vol_liq[c, jj] = smooth_min(eff_porosity[c, jj],
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
                        qflx_snow_percolation[c, jj] = smooth_max(0.0,
                            (vol_liq[c, jj] - params.ssi * eff_porosity[c, jj]) *
                            dz[c, jj] * frac_sno_eff[c])
                        qflx_snow_percolation[c, jj] = smooth_min(qflx_snow_percolation[c, jj],
                            (1.0 - vol_ice[c, jj_next] - vol_liq[c, jj_next]) *
                            dz[c, jj_next] * frac_sno_eff[c])
                    end
                else  # j == 0, bottom layer
                    qflx_snow_percolation[c, jj] = smooth_max(0.0,
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
    h2osoi_liq::Matrix{<:Real},
    dtime::Real,
    snl::Vector{Int},
    qflx_snow_percolation::Matrix{<:Real},
    mask_snow::BitVector,
    bounds::UnitRange{Int},
    nlevsno::Int
)
    snowhyd_snow_percolation!(h2osoi_liq, mask_snow, snl, qflx_snow_percolation,
        dtime, nlevsno, bounds)
end

@kernel function _snowhyd_snow_percolation_kernel!(h2osoi_liq, @Const(mask_snow),
        @Const(snl), @Const(qflx_snow_percolation), dtime, nlevsno::Int,
        cmin::Int, cmax::Int)
    c, jj = @index(Global, NTuple)
    @inbounds if cmin <= c <= cmax && mask_snow[c]
        j = jj - nlevsno  # Fortran layer index
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

snowhyd_snow_percolation!(h2osoi_liq, mask_snow, snl, qflx_snow_percolation, dtime, nlevsno, bounds) =
    _launch!(_snowhyd_snow_percolation_kernel!, h2osoi_liq, mask_snow, snl,
             qflx_snow_percolation, dtime, nlevsno, first(bounds), last(bounds);
             ndrange = (length(snl), nlevsno))

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
    dtime::Real,
    snl::Vector{Int},
    h2osoi_ice::Matrix{<:Real},
    h2osoi_liq::Matrix{<:Real},
    qflx_snow_percolation::Matrix{<:Real},
    mask_snow::BitVector,
    bounds::UnitRange{Int},
    nlevsno::Int
)
    params = snowhydrology_params
    nc = length(snl)

    # Initialize incoming fluxes
    FT = eltype(h2osoi_liq)
    qin_bc_phi = zeros(FT, nc)
    qin_bc_pho = zeros(FT, nc)
    qin_oc_phi = zeros(FT, nc)
    qin_oc_pho = zeros(FT, nc)
    qin_dst1   = zeros(FT, nc)
    qin_dst2   = zeros(FT, nc)
    qin_dst3   = zeros(FT, nc)
    qin_dst4   = zeros(FT, nc)

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
    dz::Matrix{<:Real},
    snl::Vector{Int},
    h2osoi_ice::Matrix{<:Real},
    h2osoi_liq::Matrix{<:Real},
    mask_snow::BitVector,
    bounds::UnitRange{Int},
    nlevsno::Int
)
    snowhyd_adjust_layer_thicknesses!(dz, mask_snow, snl, h2osoi_ice, h2osoi_liq,
        nlevsno, bounds)
end

@kernel function _snowhyd_adjust_layer_thicknesses_kernel!(dz, @Const(mask_snow),
        @Const(snl), @Const(h2osoi_ice), @Const(h2osoi_liq), nlevsno::Int,
        cmin::Int, cmax::Int)
    c, jj = @index(Global, NTuple)
    @inbounds if cmin <= c <= cmax && mask_snow[c]
        j = jj - nlevsno  # Fortran layer index
        if j >= snl[c] + 1
            dz[c, jj] = smooth_max(dz[c, jj],
                h2osoi_liq[c, jj] / DENH2O + h2osoi_ice[c, jj] / DENICE)
        end
    end
end

snowhyd_adjust_layer_thicknesses!(dz, mask_snow, snl, h2osoi_ice, h2osoi_liq, nlevsno, bounds) =
    _launch!(_snowhyd_adjust_layer_thicknesses_kernel!, dz, mask_snow, snl,
             h2osoi_ice, h2osoi_liq, nlevsno, first(bounds), last(bounds);
             ndrange = (length(snl), nlevsno))

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
    int_snow::Vector{<:Real},
    frac_sno::Vector{<:Real},
    snow_depth::Vector{<:Real},
    dtime::Real,
    frac_sno_eff::Vector{<:Real},
    qflx_soliddew_to_top_layer::Vector{<:Real},
    qflx_liqdew_to_top_layer::Vector{<:Real},
    qflx_liq_grnd::Vector{<:Real},
    h2osno_no_layers::Vector{<:Real},
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
    qflx_snow_drain::Vector{<:Real},
    qflx_rain_plus_snomelt::Vector{<:Real},
    frac_sno_eff::Vector{<:Real},
    qflx_snow_percolation_bottom::Vector{<:Real},
    qflx_liq_grnd::Vector{<:Real},
    qflx_snomelt::Vector{<:Real},
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
    dz::Matrix{<:Real},
    dtime::Real,
    snl::Vector{Int},
    t_soisno::Matrix{<:Real},
    h2osoi_ice::Matrix{<:Real},
    h2osoi_liq::Matrix{<:Real},
    imelt::Matrix{Int},
    frac_sno::Vector{<:Real},
    frac_h2osfc::Vector{<:Real},
    swe_old::Matrix{<:Real},
    int_snow::Vector{<:Real},
    frac_iceold::Matrix{<:Real},
    forc_wind::Vector{<:Real},
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
    FT = eltype(t_soisno)
    burden = zeros(FT, nc)
    zpseudo = zeros(FT, nc)
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
                            ddz3 = smooth_max(0.0, smooth_min(1.0, (swe_old[c, jj] - wx) / wx))
                            if (swe_old[c, jj] - wx) > 0.0
                                # Compute total snow water
                                wsum = 0.0
                                for jj2 in (snl[c]+1):0
                                    jj2_idx = jj2 + nlevsno
                                    wsum += h2osoi_liq[c, jj2_idx] + h2osoi_ice[c, jj2_idx]
                                end
                                # Simple fractional snow during melt
                                if int_snow[c] > 0.0
                                    fsno_melt = smooth_min(1.0, wsum / int_snow[c])
                                else
                                    fsno_melt = 1.0
                                end
                                if (fsno_melt + frac_h2osfc[c]) > 1.0
                                    fsno_melt = 1.0 - frac_h2osfc[c]
                                end
                                ddz3 = ddz3 - smooth_max(0.0, (fsno_melt - frac_sno[c]) / frac_sno[c])
                            end
                            ddz3 = -1.0 / dtime * ddz3
                        else
                            ddz3 = -1.0 / dtime * smooth_max(0.0,
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
                    dz[c, jj] = smooth_max(
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
function overburden_compaction_anderson1976(burden::Real, wx::Real,
                                            td::Real, bi::Real)
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
function overburden_compaction_vionnet2012(h2osoi_liq::Real, dz::Real,
                                           burden::Real, wx::Real,
                                           td::Real, bi::Real)
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
    bi::Real,
    forc_wind::Real,
    dz::Real,
    zpseudo::Vector{<:Real}, zpseudo_idx::Int,
    mobile::Vector{Bool}, mobile_idx::Int,
    compaction_rate::Ref{Float64}
)
    params = snowhydrology_params
    rho_min = snowhydrology_params.rho_min
    drift_sph = 1.0

    if mobile[mobile_idx]
        Frho = 1.25 - 0.0042 * (smooth_max(rho_min, bi) - rho_min)
        MO = 0.34 * (-0.583 * params.drift_gs - 0.833 * drift_sph + 0.833) + 0.66 * Frho
        SI = -2.868 * exp(-0.085 * forc_wind) + 1.0 + MO

        if SI > 0.0
            SI = smooth_min(SI, 3.25)
            zpseudo[zpseudo_idx] = zpseudo[zpseudo_idx] + 0.5 * dz * (3.25 - SI)
            gamma_drift = SI * exp(-zpseudo[zpseudo_idx] / 0.1)
            tau_inverse = gamma_drift / params.tau_ref
            compaction_rate[] = -smooth_max(0.0, params.rho_max - bi) * tau_inverse
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

# GPU helper: move a small host float vector (the per-layer dzmin/dzmax thresholds)
# onto the backend of a device prototype array so a kernel indexing device arrays can
# read it. On the plain CPU backend this is the identity (byte-identical). Private to
# the combine/divide snow-layer kernels below.
@inline _snow_thresh_to_device(proto, v::AbstractVector) =
    _snow_thresh_to_device(KA.get_backend(proto), proto, v)
@inline _snow_thresh_to_device(::KA.CPU, proto, v::AbstractVector) = v
@inline function _snow_thresh_to_device(backend, proto, v::AbstractVector)
    d = similar(proto, eltype(v), length(v))
    copyto!(d, v)
    return d
end

# Device-view bundles for combine_snow_layers! (placed directly above the function
# so insertions stay in a disjoint line-region for clean auto-merge). Each groups
# same-typed (float) arrays so the per-column kernel takes a few struct args instead
# of ~21 loose ones (Metal caps total kernel args ~31). Adapt-registered.
Base.@kwdef struct CombineSnowMat{M}      # per-(col,layer) snow-state float matrices
    dz::M; zi::M; z::M; t_soisno::M
    h2osoi_ice::M; h2osoi_liq::M; snw_rds::M
end
Base.@kwdef struct CombineSnowAero{M}     # per-(col,layer) aerosol mass matrices
    mss_bcphi::M; mss_bcpho::M; mss_ocphi::M; mss_ocpho::M
    mss_dst1::M; mss_dst2::M; mss_dst3::M; mss_dst4::M
end
Base.@kwdef struct CombineSnowCol{V}      # per-column float vectors
    h2osno_no_layers::V; snow_depth::V; frac_sno::V; frac_sno_eff::V; int_snow::V
end
Adapt.@adapt_structure CombineSnowMat
Adapt.@adapt_structure CombineSnowAero
Adapt.@adapt_structure CombineSnowCol

# One thread per column runs the full sequential layer search + merge in-thread on its
# own padded slice (writes only column c). dzmin/dzminloc passed as a device array (small,
# indexed by mssi). Literals converted to the working eltype T (oftype) so no Float64
# reaches a Float32-only backend (Metal) while staying byte-identical on Float64 CPU.
@kernel function _snowhyd_combine_kernel!(snl, m::CombineSnowMat, a::CombineSnowAero,
        cv::CombineSnowCol, @Const(mask_snow), @Const(lun_itype), @Const(urbpoi),
        @Const(col_landunit), @Const(dzmin), nlevsno::Int, cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && mask_snow[c]
        T = eltype(m.dz)
        zr = zero(T); thresh50 = T(50.0); ice_thresh = T(0.01); half = T(0.5)
        dzmin1 = T(dzmin[1]); lsadz = T(LSADZ)
        ncol_lyr = size(m.h2osoi_liq, 2)

        msn_old = snl[c]

        # Remove thin layers with very little ice
        j = msn_old + 1
        while j <= 0
            jj = j + nlevsno
            if m.h2osoi_ice[c, jj] <= ice_thresh
                l = col_landunit[c]
                ltype = lun_itype[l]
                if j < 0 || (ltype == ISTSOIL || urbpoi[l] || ltype == ISTCROP)
                    # Transfer water to layer below
                    jj_next = (j + 1) + nlevsno
                    if jj_next <= ncol_lyr
                        m.h2osoi_liq[c, jj_next] += m.h2osoi_liq[c, jj]
                        m.h2osoi_ice[c, jj_next] += m.h2osoi_ice[c, jj]
                    end
                end

                if j < 0
                    jj_next = (j + 1) + nlevsno
                    m.dz[c, jj_next] += m.dz[c, jj]
                    a.mss_bcphi[c, jj_next] += a.mss_bcphi[c, jj]
                    a.mss_bcpho[c, jj_next] += a.mss_bcpho[c, jj]
                    a.mss_ocphi[c, jj_next] += a.mss_ocphi[c, jj]
                    a.mss_ocpho[c, jj_next] += a.mss_ocpho[c, jj]
                    a.mss_dst1[c, jj_next]  += a.mss_dst1[c, jj]
                    a.mss_dst2[c, jj_next]  += a.mss_dst2[c, jj]
                    a.mss_dst3[c, jj_next]  += a.mss_dst3[c, jj]
                    a.mss_dst4[c, jj_next]  += a.mss_dst4[c, jj]
                end

                # Shift all elements above this down one
                if j > snl[c] + 1 && snl[c] < -1
                    for i in j:-1:(snl[c]+2)
                        ii = i + nlevsno
                        ii_above = (i - 1) + nlevsno
                        m.h2osoi_liq[c, ii] = m.h2osoi_liq[c, ii_above]
                        m.h2osoi_ice[c, ii] = m.h2osoi_ice[c, ii_above]
                        m.t_soisno[c, ii] = m.t_soisno[c, ii_above]
                        m.dz[c, ii] = m.dz[c, ii_above]
                        a.mss_bcphi[c, ii] = a.mss_bcphi[c, ii_above]
                        a.mss_bcpho[c, ii] = a.mss_bcpho[c, ii_above]
                        a.mss_ocphi[c, ii] = a.mss_ocphi[c, ii_above]
                        a.mss_ocpho[c, ii] = a.mss_ocpho[c, ii_above]
                        a.mss_dst1[c, ii] = a.mss_dst1[c, ii_above]
                        a.mss_dst2[c, ii] = a.mss_dst2[c, ii_above]
                        a.mss_dst3[c, ii] = a.mss_dst3[c, ii_above]
                        a.mss_dst4[c, ii] = a.mss_dst4[c, ii_above]
                        m.snw_rds[c, ii] = m.snw_rds[c, ii_above]
                    end
                end
                snl[c] = snl[c] + 1
            end
            j += 1
        end

        # Recompute snow depth and total water
        cv.snow_depth[c] = zr
        h2osno_total = zr
        for jl in (snl[c]+1):0
            jj = jl + nlevsno
            cv.snow_depth[c] += m.dz[c, jj]
            h2osno_total += m.h2osoi_ice[c, jj] + m.h2osoi_liq[c, jj]
        end

        # Check if all snow gone
        if cv.snow_depth[c] > zr
            l = col_landunit[c]
            ltype = lun_itype[l]
            if (ltype == ISTDLAK && cv.snow_depth[c] < dzmin1 + lsadz) ||
               (ltype != ISTDLAK &&
                (cv.frac_sno_eff[c] * cv.snow_depth[c] < dzmin1 ||
                 h2osno_total / (cv.frac_sno_eff[c] * cv.snow_depth[c]) < thresh50))

                # Transfer ice to h2osno_no_layers, liquid to top soil layer
                zwice = zr
                zwliq = zr
                for jl in (snl[c]+1):0
                    jj = jl + nlevsno
                    zwice += m.h2osoi_ice[c, jj]
                    zwliq += m.h2osoi_liq[c, jj]
                end

                cv.h2osno_no_layers[c] = zwice
                if ltype == ISTSOIL || urbpoi[l] || ltype == ISTCROP
                    # Transfer liquid to soil layer 1
                    jj_soil1 = 0 + nlevsno + 1  # layer index 1
                    if jj_soil1 <= ncol_lyr
                        m.h2osoi_liq[c, jj_soil1] += zwliq
                    end
                end

                snl[c] = 0
                h2osno_total = cv.h2osno_no_layers[c]

                # Zero out aerosol masses
                for jj in 1:nlevsno
                    a.mss_bcphi[c, jj] = zr
                    a.mss_bcpho[c, jj] = zr
                    a.mss_ocphi[c, jj] = zr
                    a.mss_ocpho[c, jj] = zr
                    a.mss_dst1[c, jj] = zr
                    a.mss_dst2[c, jj] = zr
                    a.mss_dst3[c, jj] = zr
                    a.mss_dst4[c, jj] = zr
                end

                if cv.h2osno_no_layers[c] <= zr
                    cv.snow_depth[c] = zr
                end
            end
        end

        if h2osno_total <= zr
            cv.snow_depth[c] = zr
            cv.frac_sno[c] = zr
            cv.frac_sno_eff[c] = zr
            cv.int_snow[c] = zr
        end

        # Combine thin layers (two or more layers)
        if snl[c] < -1
            mssi = 1
            i = snl[c] + 1
            while i <= 0
                jj_i = i + nlevsno
                if (cv.frac_sno_eff[c] * m.dz[c, jj_i] < T(dzmin[mssi])) ||
                   ((m.h2osoi_ice[c, jj_i] + m.h2osoi_liq[c, jj_i]) /
                    (cv.frac_sno_eff[c] * m.dz[c, jj_i]) < thresh50)

                    if i == snl[c] + 1
                        neibor = i + 1
                    elseif i == 0
                        neibor = i - 1
                    else
                        neibor = i + 1
                        jj_im1 = (i - 1) + nlevsno
                        jj_ip1 = (i + 1) + nlevsno
                        if (m.dz[c, jj_im1] + m.dz[c, jj_i]) < (m.dz[c, jj_ip1] + m.dz[c, jj_i])
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
                    a.mss_bcphi[c, jj_keep] += a.mss_bcphi[c, jj_remove]
                    a.mss_bcpho[c, jj_keep] += a.mss_bcpho[c, jj_remove]
                    a.mss_ocphi[c, jj_keep] += a.mss_ocphi[c, jj_remove]
                    a.mss_ocpho[c, jj_keep] += a.mss_ocpho[c, jj_remove]
                    a.mss_dst1[c, jj_keep]  += a.mss_dst1[c, jj_remove]
                    a.mss_dst2[c, jj_keep]  += a.mss_dst2[c, jj_remove]
                    a.mss_dst3[c, jj_keep]  += a.mss_dst3[c, jj_remove]
                    a.mss_dst4[c, jj_keep]  += a.mss_dst4[c, jj_remove]

                    # Mass-weighted snow grain size
                    total_mass = m.h2osoi_liq[c, jj_keep] + m.h2osoi_ice[c, jj_keep] +
                                 m.h2osoi_liq[c, jj_remove] + m.h2osoi_ice[c, jj_remove]
                    if total_mass > zr
                        m.snw_rds[c, jj_keep] = (m.snw_rds[c, jj_keep] *
                            (m.h2osoi_liq[c, jj_keep] + m.h2osoi_ice[c, jj_keep]) +
                            m.snw_rds[c, jj_remove] *
                            (m.h2osoi_liq[c, jj_remove] + m.h2osoi_ice[c, jj_remove])) / total_mass
                    end

                    # Combine elements using combo
                    combo!(m.dz, c, jj_keep, m.h2osoi_liq, m.h2osoi_ice, m.t_soisno,
                           m.dz[c, jj_remove], m.h2osoi_liq[c, jj_remove],
                           m.h2osoi_ice[c, jj_remove], m.t_soisno[c, jj_remove])

                    # Shift elements above
                    if j_keep - 1 > snl[c] + 1
                        for k in (j_keep-1):-1:(snl[c]+2)
                            kk = k + nlevsno
                            kk_above = (k - 1) + nlevsno
                            m.h2osoi_ice[c, kk] = m.h2osoi_ice[c, kk_above]
                            m.h2osoi_liq[c, kk] = m.h2osoi_liq[c, kk_above]
                            m.t_soisno[c, kk] = m.t_soisno[c, kk_above]
                            m.dz[c, kk] = m.dz[c, kk_above]
                            a.mss_bcphi[c, kk] = a.mss_bcphi[c, kk_above]
                            a.mss_bcpho[c, kk] = a.mss_bcpho[c, kk_above]
                            a.mss_ocphi[c, kk] = a.mss_ocphi[c, kk_above]
                            a.mss_ocpho[c, kk] = a.mss_ocpho[c, kk_above]
                            a.mss_dst1[c, kk] = a.mss_dst1[c, kk_above]
                            a.mss_dst2[c, kk] = a.mss_dst2[c, kk_above]
                            a.mss_dst3[c, kk] = a.mss_dst3[c, kk_above]
                            a.mss_dst4[c, kk] = a.mss_dst4[c, kk_above]
                            m.snw_rds[c, kk] = m.snw_rds[c, kk_above]
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
        for jl in 0:-1:(-nlevsno+1)
            jj = jl + nlevsno
            if jl >= snl[c] + 1
                m.z[c, jj] = m.zi[c, jj + 1] - half * m.dz[c, jj]
                m.zi[c, jj] = m.zi[c, jj + 1] - m.dz[c, jj]
            end
        end
    end
end

"""
    combine_snow_layers!(snl, dz, zi, z, t_soisno,
        h2osoi_ice, h2osoi_liq, h2osno_no_layers,
        snow_depth, frac_sno, frac_sno_eff, int_snow,
        snw_rds, aerosol,
        lun_itype, urbpoi,
        mask_snow, bounds, nlevsno)

Combine snow layers that are less than minimum thickness or mass. One per-column
kernel; backend-agnostic (CPU loop or whole-function on GPU).
"""
function combine_snow_layers!(
    snl::AbstractVector{Int},
    dz::AbstractMatrix{<:Real},
    zi::AbstractMatrix{<:Real},
    z::AbstractMatrix{<:Real},
    t_soisno::AbstractMatrix{<:Real},
    h2osoi_ice::AbstractMatrix{<:Real},
    h2osoi_liq::AbstractMatrix{<:Real},
    h2osno_no_layers::AbstractVector{<:Real},
    snow_depth::AbstractVector{<:Real},
    frac_sno::AbstractVector{<:Real},
    frac_sno_eff::AbstractVector{<:Real},
    int_snow::AbstractVector{<:Real},
    snw_rds::AbstractMatrix{<:Real},
    aerosol::AerosolData,
    lun_itype::AbstractVector{Int},     # landunit type indexed by landunit
    urbpoi::AbstractVector{Bool},       # urban point indexed by landunit
    col_landunit::AbstractVector{Int},  # column-to-landunit mapping
    mask_snow::AbstractVector{Bool},
    bounds::UnitRange{Int},
    nlevsno::Int
)
    # dzmin/dzminloc as a device-resident array (small, indexed by layer); built to the
    # working precision so no Float64 reaches a Float32-only backend.
    FT = eltype(dz)
    dzmin = _snow_thresh_to_device(dz, FT.(SNOW_DZMIN[]))

    m = CombineSnowMat(; dz, zi, z, t_soisno, h2osoi_ice, h2osoi_liq, snw_rds)
    a = CombineSnowAero(;
        mss_bcphi = aerosol.mss_bcphi_col, mss_bcpho = aerosol.mss_bcpho_col,
        mss_ocphi = aerosol.mss_ocphi_col, mss_ocpho = aerosol.mss_ocpho_col,
        mss_dst1 = aerosol.mss_dst1_col, mss_dst2 = aerosol.mss_dst2_col,
        mss_dst3 = aerosol.mss_dst3_col, mss_dst4 = aerosol.mss_dst4_col)
    cv = CombineSnowCol(; h2osno_no_layers, snow_depth, frac_sno, frac_sno_eff, int_snow)

    _launch!(_snowhyd_combine_kernel!, snl, m, a, cv, mask_snow, lun_itype, urbpoi,
             col_landunit, dzmin, nlevsno, first(bounds), last(bounds);
             ndrange = length(snl))
    return nothing
end

# =========================================================================
# Combo (helper for combining snow layers)
# =========================================================================

"""
    combo!(dz_mat, c, jj, wliq_mat, wice_mat, t_mat,
           dz2, wliq2, wice2, t2)

Combine two snow elements. Updates element at [c, jj] in-place.
"""
function combo!(dz_mat::AbstractMatrix{<:Real}, c::Int, jj::Int,
                wliq_mat::AbstractMatrix{<:Real}, wice_mat::AbstractMatrix{<:Real},
                t_mat::AbstractMatrix{<:Real},
                dz2::Real, wliq2::Real, wice2::Real, t2::Real)
    dz1 = dz_mat[c, jj]
    wliq1 = wliq_mat[c, jj]
    wice1 = wice_mat[c, jj]
    t1 = t_mat[c, jj]

    (dzc, wliqc, wicec, tc) = combo_scalar(dz1, wliq1, wice1, t1, dz2, wliq2, wice2, t2)

    dz_mat[c, jj] = dzc
    wice_mat[c, jj] = wicec
    wliq_mat[c, jj] = wliqc
    t_mat[c, jj] = tc
end

"""
    combo_scalar(dz, wliq, wice, t, dz2, wliq2, wice2, t2) -> (dz, wliq, wice, t)

Scalar version of combo for use in divide_snow_layers.
"""
@inline function combo_scalar(dz::Real, wliq::Real, wice::Real, t::Real,
                      dz2::Real, wliq2::Real, wice2::Real, t2::Real)
    # eltype-generic so the same code is device-safe (Metal/Float32) and byte-identical
    # on Float64 CPU: each Float64 constant is converted to the working precision via
    # oftype (oftype(::Float64, x) === x).
    T = promote_type(typeof(dz), typeof(wliq), typeof(wice), typeof(t),
                     typeof(dz2), typeof(wliq2), typeof(wice2), typeof(t2))
    cpice = T(CPICE); cpliq = T(CPLIQ); hfus = T(HFUS); tfrz = T(TFRZ)

    dzc = dz + dz2
    wicec = wice + wice2
    wliqc = wliq + wliq2

    h = (cpice * wice + cpliq * wliq) * (t - tfrz) + hfus * wliq
    h2val = (cpice * wice2 + cpliq * wliq2) * (t2 - tfrz) + hfus * wliq2
    hc = h + h2val

    denom = cpice * wicec + cpliq * wliqc
    if denom > zero(T)
        tc = tfrz + (hc - hfus * wliqc) / denom
    else
        tc = tfrz
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
# Host-facing convenience method: reads the module-level params/const (host only).
function mass_weighted_snow_radius(rds1::Real, rds2::Real,
                                    swtot::Real, zwtot::Real)
    params = snowhydrology_params
    return mass_weighted_snow_radius(rds1, rds2, swtot, zwtot,
                                     params.snw_rds_min, SNW_RDS_MAX)
end

# Device-safe eltype-generic core: snw_rds_min / snw_rds_max passed in (no host-global
# deref inside a GPU kernel). Byte-identical on Float64 (oftype(::Float64, x) === x).
@inline function mass_weighted_snow_radius(rds1::Real, rds2::Real,
                                    swtot::Real, zwtot::Real,
                                    snw_rds_min::Real, snw_rds_max::Real)
    T = promote_type(typeof(rds1), typeof(rds2), typeof(swtot), typeof(zwtot))
    total_wt = swtot + zwtot
    result = total_wt > zero(T) ? (rds2 * swtot + rds1 * zwtot) / total_wt :
                                  T(0.5) * (rds1 + rds2)
    result = smooth_clamp(result, T(snw_rds_min), T(snw_rds_max))
    return result
end

# =========================================================================
# DivideSnowLayers
# =========================================================================

# Per-column device-resident scratch for divide_snow_layers! (the Fortran "local arrays"
# dzsno/swice/… indexed 1..nlevsno, top=1). One row per column so each thread owns a
# disjoint slice; bundled into two structs to stay under the Metal ~31-arg kernel cap.
# Placed directly above the function for clean auto-merge. Adapt-registered.
Base.@kwdef struct DivideSnowScratch1{M}
    dzsno::M; swice::M; swliq::M; tsno::M; rds::M
end
Base.@kwdef struct DivideSnowScratch2{M}
    mbc_phi::M; mbc_pho::M; moc_phi::M; moc_pho::M
    mdst1::M; mdst2::M; mdst3::M; mdst4::M
end
Adapt.@adapt_structure DivideSnowScratch1
Adapt.@adapt_structure DivideSnowScratch2

# One thread per column runs the full sequential top-to-bottom split in-thread, using
# its own scratch row (s1/s2) and writing only column c. dzmax_l/dzmax_u as device arrays;
# is_lake/snw_rds_min/snw_rds_max resolved on host and passed as scalars. Literals converted
# to the working eltype T so no Float64 reaches a Float32-only backend (Metal) while staying
# byte-identical on Float64 CPU.
@kernel function _snowhyd_divide_kernel!(snl, m::CombineSnowMat, a::CombineSnowAero,
        s1::DivideSnowScratch1, s2::DivideSnowScratch2, @Const(frac_sno),
        @Const(mask_snow), @Const(dzmax_l), @Const(dzmax_u), is_lake::Bool,
        snw_rds_min, snw_rds_max, nlevsno::Int, cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && mask_snow[c]
        T = eltype(m.dz)
        zr = zero(T); half = T(0.5); two = T(2.0)
        lsadz = T(LSADZ); tfrz = T(TFRZ)
        rmin = T(snw_rds_min); rmax = T(snw_rds_max)

        msno = abs(snl[c])

        # Copy to scratch (1-indexed, top=1). Zero-fill the whole scratch row first
        # (similar() is uninitialized), mirroring the original zeros(FT, nlevsno).
        for k in 1:nlevsno
            s1.dzsno[c, k] = zr; s1.swice[c, k] = zr; s1.swliq[c, k] = zr
            s1.tsno[c, k] = zr; s1.rds[c, k] = zr
            s2.mbc_phi[c, k] = zr; s2.mbc_pho[c, k] = zr
            s2.moc_phi[c, k] = zr; s2.moc_pho[c, k] = zr
            s2.mdst1[c, k] = zr; s2.mdst2[c, k] = zr; s2.mdst3[c, k] = zr; s2.mdst4[c, k] = zr
        end

        for k in 1:msno
            jj = (k + snl[c]) + nlevsno  # Fortran j+snl(c) -> Julia index
            if is_lake
                s1.dzsno[c, k] = m.dz[c, jj]
            else
                s1.dzsno[c, k] = frac_sno[c] * m.dz[c, jj]
            end
            s1.swice[c, k] = m.h2osoi_ice[c, jj]
            s1.swliq[c, k] = m.h2osoi_liq[c, jj]
            s1.tsno[c, k]  = m.t_soisno[c, jj]
            s2.mbc_phi[c, k] = a.mss_bcphi[c, jj]
            s2.mbc_pho[c, k] = a.mss_bcpho[c, jj]
            s2.moc_phi[c, k] = a.mss_ocphi[c, jj]
            s2.moc_pho[c, k] = a.mss_ocpho[c, jj]
            s2.mdst1[c, k] = a.mss_dst1[c, jj]
            s2.mdst2[c, k] = a.mss_dst2[c, jj]
            s2.mdst3[c, k] = a.mss_dst3[c, jj]
            s2.mdst4[c, k] = a.mss_dst4[c, jj]
            s1.rds[c, k]   = m.snw_rds[c, jj]
        end

        # Traverse layers top to bottom
        k = 1
        while k <= msno && k < nlevsno
            # Bottom layer case
            if k == msno
                offset = is_lake ? two * lsadz : zr
                if s1.dzsno[c, k] > dzmax_l[k] + offset
                    msno += 1
                    s1.dzsno[c, k] /= two
                    s1.dzsno[c, k+1] = s1.dzsno[c, k]
                    s1.swice[c, k] /= two; s1.swice[c, k+1] = s1.swice[c, k]
                    s1.swliq[c, k] /= two; s1.swliq[c, k+1] = s1.swliq[c, k]

                    if k == 1
                        s1.tsno[c, k+1] = s1.tsno[c, k]
                    else
                        dtdz = (s1.tsno[c, k-1] - s1.tsno[c, k]) /
                               ((s1.dzsno[c, k-1] + two*s1.dzsno[c, k]) / two)
                        tsno_kp1_candidate = s1.tsno[c, k] - dtdz * s1.dzsno[c, k] / two
                        tsno_k_candidate = s1.tsno[c, k] + dtdz * s1.dzsno[c, k] / two
                        # Smooth blend: if candidate >= TFRZ, use tsno[k]; else use candidates
                        w_frozen = smooth_heaviside(tfrz - tsno_kp1_candidate)
                        s1.tsno[c, k+1] = w_frozen * tsno_kp1_candidate + (one(T) - w_frozen) * s1.tsno[c, k]
                        s1.tsno[c, k] = w_frozen * tsno_k_candidate + (one(T) - w_frozen) * s1.tsno[c, k]
                    end

                    s2.mbc_phi[c, k] /= two; s2.mbc_phi[c, k+1] = s2.mbc_phi[c, k]
                    s2.mbc_pho[c, k] /= two; s2.mbc_pho[c, k+1] = s2.mbc_pho[c, k]
                    s2.moc_phi[c, k] /= two; s2.moc_phi[c, k+1] = s2.moc_phi[c, k]
                    s2.moc_pho[c, k] /= two; s2.moc_pho[c, k+1] = s2.moc_pho[c, k]
                    s2.mdst1[c, k] /= two; s2.mdst1[c, k+1] = s2.mdst1[c, k]
                    s2.mdst2[c, k] /= two; s2.mdst2[c, k+1] = s2.mdst2[c, k]
                    s2.mdst3[c, k] /= two; s2.mdst3[c, k+1] = s2.mdst3[c, k]
                    s2.mdst4[c, k] /= two; s2.mdst4[c, k+1] = s2.mdst4[c, k]
                    s1.rds[c, k+1] = s1.rds[c, k]
                end
            end

            # Layers below exist
            if k < msno
                offset = is_lake ? lsadz : zr
                if s1.dzsno[c, k] > dzmax_u[k] + offset
                    drr = s1.dzsno[c, k] - dzmax_u[k] - offset
                    propor = drr / s1.dzsno[c, k]

                    zwice = propor * s1.swice[c, k]
                    zwliq = propor * s1.swliq[c, k]
                    zmbc_phi = propor * s2.mbc_phi[c, k]
                    zmbc_pho = propor * s2.mbc_pho[c, k]
                    zmoc_phi = propor * s2.moc_phi[c, k]
                    zmoc_pho = propor * s2.moc_pho[c, k]
                    zmdst1 = propor * s2.mdst1[c, k]
                    zmdst2 = propor * s2.mdst2[c, k]
                    zmdst3 = propor * s2.mdst3[c, k]
                    zmdst4 = propor * s2.mdst4[c, k]

                    propor_keep = (dzmax_u[k] + offset) / s1.dzsno[c, k]
                    s1.swice[c, k] *= propor_keep
                    s1.swliq[c, k] *= propor_keep
                    s2.mbc_phi[c, k] *= propor_keep
                    s2.mbc_pho[c, k] *= propor_keep
                    s2.moc_phi[c, k] *= propor_keep
                    s2.moc_pho[c, k] *= propor_keep
                    s2.mdst1[c, k] *= propor_keep
                    s2.mdst2[c, k] *= propor_keep
                    s2.mdst3[c, k] *= propor_keep
                    s2.mdst4[c, k] *= propor_keep

                    s1.dzsno[c, k] = dzmax_u[k] + offset

                    s2.mbc_phi[c, k+1] += zmbc_phi
                    s2.mbc_pho[c, k+1] += zmbc_pho
                    s2.moc_phi[c, k+1] += zmoc_phi
                    s2.moc_pho[c, k+1] += zmoc_pho
                    s2.mdst1[c, k+1] += zmdst1
                    s2.mdst2[c, k+1] += zmdst2
                    s2.mdst3[c, k+1] += zmdst3
                    s2.mdst4[c, k+1] += zmdst4

                    s1.rds[c, k+1] = mass_weighted_snow_radius(s1.rds[c, k], s1.rds[c, k+1],
                        s1.swliq[c, k+1] + s1.swice[c, k+1], zwliq + zwice, rmin, rmax)

                    (s1.dzsno[c, k+1], s1.swliq[c, k+1], s1.swice[c, k+1], s1.tsno[c, k+1]) =
                        combo_scalar(s1.dzsno[c, k+1], s1.swliq[c, k+1], s1.swice[c, k+1], s1.tsno[c, k+1],
                                     drr, zwliq, zwice, s1.tsno[c, k])
                end
            end
            k += 1
        end

        snl[c] = -msno

        # Copy back from scratch to column arrays
        for jl in (-nlevsno+1):0
            jj = jl + nlevsno
            if jl >= snl[c] + 1
                k_idx = jl - snl[c]
                if is_lake
                    m.dz[c, jj] = s1.dzsno[c, k_idx]
                else
                    m.dz[c, jj] = s1.dzsno[c, k_idx] / frac_sno[c]
                end
                m.h2osoi_ice[c, jj] = s1.swice[c, k_idx]
                m.h2osoi_liq[c, jj] = s1.swliq[c, k_idx]
                m.t_soisno[c, jj] = s1.tsno[c, k_idx]
                a.mss_bcphi[c, jj] = s2.mbc_phi[c, k_idx]
                a.mss_bcpho[c, jj] = s2.mbc_pho[c, k_idx]
                a.mss_ocphi[c, jj] = s2.moc_phi[c, k_idx]
                a.mss_ocpho[c, jj] = s2.moc_pho[c, k_idx]
                a.mss_dst1[c, jj] = s2.mdst1[c, k_idx]
                a.mss_dst2[c, jj] = s2.mdst2[c, k_idx]
                a.mss_dst3[c, jj] = s2.mdst3[c, k_idx]
                a.mss_dst4[c, jj] = s2.mdst4[c, k_idx]
                m.snw_rds[c, jj] = s1.rds[c, k_idx]
            end
        end

        # Reset node depth and interface depth
        for jl in 0:-1:(-nlevsno+1)
            jj = jl + nlevsno
            if jl >= snl[c] + 1
                m.z[c, jj] = m.zi[c, jj + 1] - half * m.dz[c, jj]
                m.zi[c, jj] = m.zi[c, jj + 1] - m.dz[c, jj]
            end
        end
    end
end

"""
    divide_snow_layers!(snl, dz, zi, z, t_soisno,
        h2osoi_ice, h2osoi_liq, frac_sno, snw_rds, aerosol,
        is_lake, mask_snow, bounds, nlevsno)

Subdivide snow layers that exceed their prescribed maximum thickness. One per-column
kernel; backend-agnostic (CPU loop or whole-function on GPU).
"""
function divide_snow_layers!(
    snl::AbstractVector{Int},
    dz::AbstractMatrix{<:Real},
    zi::AbstractMatrix{<:Real},
    z::AbstractMatrix{<:Real},
    t_soisno::AbstractMatrix{<:Real},
    h2osoi_ice::AbstractMatrix{<:Real},
    h2osoi_liq::AbstractMatrix{<:Real},
    frac_sno::AbstractVector{<:Real},
    snw_rds::AbstractMatrix{<:Real},
    aerosol::AerosolData,
    is_lake::Bool,
    mask_snow::AbstractVector{Bool},
    bounds::UnitRange{Int},
    nlevsno::Int
)
    FT = eltype(dz)
    dzmax_l = _snow_thresh_to_device(dz, FT.(SNOW_DZMAX_L[]))
    dzmax_u = _snow_thresh_to_device(dz, FT.(SNOW_DZMAX_U[]))

    nc = length(snl)
    # Per-column device-resident scratch (one row per column; columns [1..nlevsno]).
    # similar() is device-resident (NOT zeros()); the kernel zero-fills its own row.
    sc1 = DivideSnowScratch1(;
        dzsno = similar(dz, FT, nc, nlevsno), swice = similar(dz, FT, nc, nlevsno),
        swliq = similar(dz, FT, nc, nlevsno), tsno = similar(dz, FT, nc, nlevsno),
        rds = similar(dz, FT, nc, nlevsno))
    sc2 = DivideSnowScratch2(;
        mbc_phi = similar(dz, FT, nc, nlevsno), mbc_pho = similar(dz, FT, nc, nlevsno),
        moc_phi = similar(dz, FT, nc, nlevsno), moc_pho = similar(dz, FT, nc, nlevsno),
        mdst1 = similar(dz, FT, nc, nlevsno), mdst2 = similar(dz, FT, nc, nlevsno),
        mdst3 = similar(dz, FT, nc, nlevsno), mdst4 = similar(dz, FT, nc, nlevsno))

    m = CombineSnowMat(; dz, zi, z, t_soisno, h2osoi_ice, h2osoi_liq, snw_rds)
    a = CombineSnowAero(;
        mss_bcphi = aerosol.mss_bcphi_col, mss_bcpho = aerosol.mss_bcpho_col,
        mss_ocphi = aerosol.mss_ocphi_col, mss_ocpho = aerosol.mss_ocpho_col,
        mss_dst1 = aerosol.mss_dst1_col, mss_dst2 = aerosol.mss_dst2_col,
        mss_dst3 = aerosol.mss_dst3_col, mss_dst4 = aerosol.mss_dst4_col)

    snw_rds_min = FT(snowhydrology_params.snw_rds_min)
    snw_rds_max = FT(SNW_RDS_MAX)

    _launch!(_snowhyd_divide_kernel!, snl, m, a, sc1, sc2, frac_sno, mask_snow,
             dzmax_l, dzmax_u, is_lake, snw_rds_min, snw_rds_max, nlevsno,
             first(bounds), last(bounds); ndrange = length(snl))
    return nothing
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
    dz::Matrix{<:Real},
    z::Matrix{<:Real},
    zi::Matrix{<:Real},
    t_soisno::Matrix{<:Real},
    h2osoi_ice::Matrix{<:Real},
    h2osoi_liq::Matrix{<:Real},
    mask_snow::BitVector,
    bounds::UnitRange{Int},
    nlevsno::Int
)
    snowhyd_zero_empty_snow_layers!(snl, dz, z, zi, t_soisno, h2osoi_ice, h2osoi_liq,
        mask_snow, nlevsno, bounds)
end

@kernel function _snowhyd_zero_empty_layers_kernel!(h2osoi_ice, @Const(mask_snow),
        @Const(snl), h2osoi_liq, t_soisno, dz, z, zi, nlevsno::Int,
        cmin::Int, cmax::Int)
    c, jj = @index(Global, NTuple)
    @inbounds if cmin <= c <= cmax && mask_snow[c]
        j = jj - nlevsno  # Fortran layer index
        if j <= snl[c] && snl[c] > -nlevsno
            h2osoi_ice[c, jj] = 0.0
            h2osoi_liq[c, jj] = 0.0
            t_soisno[c, jj] = 0.0
            dz[c, jj] = 0.0
            z[c, jj] = 0.0
            zi[c, jj] = 0.0
        end
    end
end

snowhyd_zero_empty_snow_layers!(snl, dz, z, zi, t_soisno, h2osoi_ice, h2osoi_liq, mask_snow, nlevsno, bounds) =
    _launch!(_snowhyd_zero_empty_layers_kernel!, h2osoi_ice, mask_snow, snl,
             h2osoi_liq, t_soisno, dz, z, zi, nlevsno, first(bounds), last(bounds);
             ndrange = (length(snl), nlevsno))

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
    dz::Matrix{<:Real},
    z::Matrix{<:Real},
    zi::Matrix{<:Real},
    snow_depth::Vector{<:Real},
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
        jj_surface_zi = nlevsno + 1

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
            zi[c, jj_surface_zi] = 0.0
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
            zi[c, jj_surface_zi] = 0.0
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
        zi[c, jj_surface_zi] = 0.0
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
            z[c, jj] = zi[c, jj + 1] - 0.5 * dz[c, jj]
            zi[c, jj] = zi[c, jj + 1] - dz[c, jj]
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
    snowhyd_build_snow_filter!(mask_snow, mask_nosnow, snl, mask_nolake, bounds)
end

@kernel function _snowhyd_build_snow_filter_kernel!(mask_snow, @Const(snl),
        @Const(mask_nolake), mask_nosnow, cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && mask_nolake[c]
        if snl[c] < 0
            mask_snow[c] = true
        else
            mask_nosnow[c] = true
        end
    end
end

snowhyd_build_snow_filter!(mask_snow, mask_nosnow, snl, mask_nolake, bounds) =
    _launch!(_snowhyd_build_snow_filter_kernel!, mask_snow, snl, mask_nolake,
             mask_nosnow, first(bounds), last(bounds))

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
    h2osno_excess::Vector{<:Real},
    apply_runoff::AbstractVector{Bool},
    h2osno::Vector{<:Real},
    topo::Vector{<:Real},
    snl::Vector{Int},
    col_landunit::Vector{Int},
    lun_itype::Vector{Int},   # landunit type
    mask_snow::BitVector,
    bounds::UnitRange{Int},
    nstep::Int,
    nlevsno::Int
)
    for c in bounds
        mask_snow[c] || continue
        h2osno_excess[c] = 0.0
        if h2osno[c] > H2OSNO_MAX
            h2osno_excess[c] = h2osno[c] - H2OSNO_MAX
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
    qflx_snwcp_ice::Vector{<:Real},
    qflx_snwcp_liq::Vector{<:Real},
    qflx_snwcp_discarded_ice::Vector{<:Real},
    qflx_snwcp_discarded_liq::Vector{<:Real},
    mask::BitVector,
    bounds::UnitRange{Int}
)
    snowhyd_init_flux_snow_capping!(qflx_snwcp_ice, qflx_snwcp_liq,
        qflx_snwcp_discarded_ice, qflx_snwcp_discarded_liq, mask, bounds)
end

@kernel function _snowhyd_init_flux_capping_kernel!(qflx_snwcp_ice, @Const(mask),
        qflx_snwcp_liq, qflx_snwcp_discarded_ice, qflx_snwcp_discarded_liq,
        cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && mask[c]
        qflx_snwcp_ice[c] = 0.0
        qflx_snwcp_liq[c] = 0.0
        qflx_snwcp_discarded_ice[c] = 0.0
        qflx_snwcp_discarded_liq[c] = 0.0
    end
end

snowhyd_init_flux_snow_capping!(qflx_snwcp_ice, qflx_snwcp_liq, qflx_snwcp_discarded_ice, qflx_snwcp_discarded_liq, mask, bounds) =
    _launch!(_snowhyd_init_flux_capping_kernel!, qflx_snwcp_ice, mask, qflx_snwcp_liq,
             qflx_snwcp_discarded_ice, qflx_snwcp_discarded_liq,
             first(bounds), last(bounds))

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
    rho_orig_bottom::Vector{<:Real},
    frac_adjust::Vector{<:Real},
    qflx_snwcp_ice::Vector{<:Real},
    qflx_snwcp_liq::Vector{<:Real},
    qflx_snwcp_discarded_ice::Vector{<:Real},
    qflx_snwcp_discarded_liq::Vector{<:Real},
    dtime::Real,
    dz_bottom::Vector{<:Real},
    topo::Vector{<:Real},
    h2osno_total::Vector{<:Real},
    h2osoi_ice_bottom::Vector{<:Real},
    h2osoi_liq_bottom::Vector{<:Real},
    col_landunit::Vector{Int},
    lun_itype::Vector{Int},
    mask_snow::BitVector,
    bounds::UnitRange{Int},
    nlevsno::Int,
    nstep::Int
)

    nc = length(h2osno_total)
    FT = eltype(h2osno_total)
    h2osno_excess = zeros(FT, nc)
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

        rho_orig_bottom[c] = dz_bottom[c] > eps(Float64) ?
            h2osoi_ice_bottom[c] / dz_bottom[c] : 0.0

        mss_snow_bottom_lyr = h2osoi_ice_bottom[c] + h2osoi_liq_bottom[c]
        mss_snwcp_tot = smooth_min(h2osno_excess[c],
            mss_snow_bottom_lyr * (1.0 - MIN_SNOW_TO_KEEP))

        if mss_snow_bottom_lyr > eps(Float64)
            icefrac = h2osoi_ice_bottom[c] / mss_snow_bottom_lyr
        else
            icefrac = 0.0
        end
        snwcp_flux_ice = mss_snwcp_tot / dtime * icefrac
        snwcp_flux_liq = mss_snwcp_tot / dtime * (1.0 - icefrac)

        if apply_runoff[c]
            qflx_snwcp_ice[c] = snwcp_flux_ice
            qflx_snwcp_liq[c] = snwcp_flux_liq
        else
            qflx_snwcp_discarded_ice[c] = snwcp_flux_ice
            qflx_snwcp_discarded_liq[c] = snwcp_flux_liq
        end

        frac_adjust[c] = mss_snow_bottom_lyr > eps(Float64) ?
            (mss_snow_bottom_lyr - mss_snwcp_tot) / mss_snow_bottom_lyr : 1.0
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
    h2osoi_ice_bottom::Vector{<:Real},
    h2osoi_liq_bottom::Vector{<:Real},
    dtime::Real,
    qflx_snwcp_ice::Vector{<:Real},
    qflx_snwcp_liq::Vector{<:Real},
    qflx_snwcp_discarded_ice::Vector{<:Real},
    qflx_snwcp_discarded_liq::Vector{<:Real},
    mask_capping::BitVector,
    bounds::UnitRange{Int}
)
    for c in bounds
        mask_capping[c] || continue
        h2osoi_ice_bottom[c] -= (qflx_snwcp_ice[c] + qflx_snwcp_discarded_ice[c]) * dtime
        h2osoi_liq_bottom[c] -= (qflx_snwcp_liq[c] + qflx_snwcp_discarded_liq[c]) * dtime

        if h2osoi_ice_bottom[c] < 0.0 || h2osoi_liq_bottom[c] < 0.0
            if _is_ad_type(eltype(h2osoi_ice_bottom))
                @warn "snow capping: negative mass (AD mode, clamping)" maxlog=1
                h2osoi_ice_bottom[c] = smooth_max(0.0, h2osoi_ice_bottom[c])
                h2osoi_liq_bottom[c] = smooth_max(0.0, h2osoi_liq_bottom[c])
            else
                error("snow capping failed: negative mass at c=$c, ice=$(h2osoi_ice_bottom[c]), liq=$(h2osoi_liq_bottom[c])")
            end
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
    dz_bottom::Vector{<:Real},
    aerosol::AerosolData,
    jj_bottom::Int,
    rho_orig_bottom::Vector{<:Real},
    h2osoi_ice_bottom::Vector{<:Real},
    frac_adjust::Vector{<:Real},
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
