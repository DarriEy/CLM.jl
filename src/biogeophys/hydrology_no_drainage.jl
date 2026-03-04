# ==========================================================================
# Ported from: src/biogeophys/HydrologyNoDrainageMod.F90 (~727 lines)
# Calculate snow and soil hydrology without drainage.
#
# Public functions:
#   update_snow_persistence!               — Update snow persistence counters
#   accumulate_snow_ice_liq!               — Sum snowice/snowliq over snow layers
#   update_snowdp!                         — Column-average snow depth
#   compute_snow_internal_temperature!     — Snow internal temperature (t_sno_mul_mss)
#   update_ground_and_soil_temperatures!   — Ground temp, tsl, t_soi_10cm, tsoi17
#   update_h2osoi_vol!                     — Volumetric soil water content
#   update_soilpsi!                        — Soil water potential (MPa)
#   update_smp_l!                          — Matric potential for history/ch4
#   compute_wf!                            — Soil water fraction of WHC (top 0.05m)
#   compute_wf2!                           — Soil water fraction of WHC (top 0.17m)
#   update_snow_top_layer_diagnostics!     — Top-layer snow diagnostics
#   hydrology_no_drainage!                 — Main orchestrator (calls sub-functions)
#   calc_and_withdraw_irrigation_fluxes!   — Irrigation withdrawal (stub)
#   handle_new_snow!                       — Handle new snow falling on ground
# ==========================================================================

# =========================================================================
# update_snow_persistence!
# =========================================================================

"""
    update_snow_persistence!(snow_persistence, dtime,
        mask_snow, mask_nosnow, bounds)

For columns with snow, accumulate time-covered-by-snow counter.
For columns without snow, reset counter to zero.

Ported from inline code in `HydrologyNoDrainage` in `HydrologyNoDrainageMod.F90`.
"""
function update_snow_persistence!(
    snow_persistence::Vector{Float64},
    dtime::Float64,
    mask_snow::BitVector,
    mask_nosnow::BitVector,
    bounds::UnitRange{Int}
)
    for c in bounds
        if mask_snow[c]
            snow_persistence[c] = snow_persistence[c] + dtime
        elseif mask_nosnow[c]
            snow_persistence[c] = 0.0
        end
    end
    return nothing
end

# =========================================================================
# accumulate_snow_ice_liq!
# =========================================================================

"""
    accumulate_snow_ice_liq!(snowice, snowliq, h2osoi_ice, h2osoi_liq,
        snl, mask_nolake, mask_snow, bounds, nlevsno)

Vertically sum h2osoi_ice and h2osoi_liq over all snow layers for history output.

Ported from inline code in `HydrologyNoDrainage` in `HydrologyNoDrainageMod.F90`.
"""
function accumulate_snow_ice_liq!(
    snowice::Vector{Float64},
    snowliq::Vector{Float64},
    h2osoi_ice::Matrix{Float64},
    h2osoi_liq::Matrix{Float64},
    snl::Vector{Int},
    mask_nolake::BitVector,
    mask_snow::BitVector,
    bounds::UnitRange{Int},
    nlevsno::Int
)
    # Zero for all non-lake columns
    for c in bounds
        mask_nolake[c] || continue
        snowice[c] = 0.0
        snowliq[c] = 0.0
    end

    # Accumulate over snow layers
    # Fortran: j from -nlevsno+1 to 0 → Julia offset: j_julia = j_fortran + nlevsno
    for j_fortran in (-nlevsno + 1):0
        j = j_fortran + nlevsno  # Julia 1-based index
        for c in bounds
            mask_snow[c] || continue
            if j_fortran >= snl[c] + 1
                snowice[c] = snowice[c] + h2osoi_ice[c, j]
                snowliq[c] = snowliq[c] + h2osoi_liq[c, j]
            end
        end
    end

    return nothing
end

# =========================================================================
# update_snowdp!
# =========================================================================

"""
    update_snowdp!(snowdp, snow_depth, frac_sno_eff, bounds)

Calculate column average snow depth = snow_depth * frac_sno_eff.

Ported from inline code in `HydrologyNoDrainage` in `HydrologyNoDrainageMod.F90`.
"""
function update_snowdp!(
    snowdp::Vector{Float64},
    snow_depth::Vector{Float64},
    frac_sno_eff::Vector{Float64},
    bounds::UnitRange{Int}
)
    for c in bounds
        snowdp[c] = snow_depth[c] * frac_sno_eff[c]
    end
    return nothing
end

# =========================================================================
# compute_snow_internal_temperature!
# =========================================================================

"""
    compute_snow_internal_temperature!(t_sno_mul_mss, h2osoi_ice, h2osoi_liq,
        t_soisno, snl, mask_nolake, mask_snow, bounds, nlevsno, tfrz)

Compute the snow internal temperature diagnostic: sum of layer mass × temperature.

The snow internal temperature (SIT) is the weighted average temperature of the
snowpack, weighted by layer mass. This function computes the numerator:
  Sum_i [ (ice_mass(i) * T(i)) + (liq_mass(i) * tfrz) ]

Ported from inline code in `HydrologyNoDrainage` in `HydrologyNoDrainageMod.F90`.
"""
function compute_snow_internal_temperature!(
    t_sno_mul_mss::Vector{Float64},
    h2osoi_ice::Matrix{Float64},
    h2osoi_liq::Matrix{Float64},
    t_soisno::Matrix{Float64},
    snl::Vector{Int},
    mask_nolake::BitVector,
    mask_snow::BitVector,
    bounds::UnitRange{Int},
    nlevsno::Int,
    tfrz::Float64
)
    # Zero for all non-lake columns
    for c in bounds
        mask_nolake[c] || continue
        t_sno_mul_mss[c] = 0.0
    end

    # Accumulate over snow layers
    for j_fortran in (-nlevsno + 1):0
        j = j_fortran + nlevsno  # Julia 1-based index
        for c in bounds
            mask_snow[c] || continue
            if j_fortran >= snl[c] + 1
                t_sno_mul_mss[c] = t_sno_mul_mss[c] + h2osoi_ice[c, j] * t_soisno[c, j]
                t_sno_mul_mss[c] = t_sno_mul_mss[c] + h2osoi_liq[c, j] * tfrz
            end
        end
    end

    return nothing
end

# =========================================================================
# update_ground_and_soil_temperatures!
# =========================================================================

"""
    update_ground_and_soil_temperatures!(t_grnd, t_grnd_u, t_grnd_r,
        tsl, t_soi_10cm, tsoi17,
        t_soisno, t_h2osfc, frac_sno_eff, frac_h2osfc,
        snl, dz, zi, col_landunit, col_itype,
        lun_urbpoi, lun_itype,
        mask_nolake, bounds, nlevsoi, nlevsno)

Compute ground temperature, near-surface soil temperature, and
depth-weighted soil temperatures for top 10cm and 17cm.

- t_grnd: weighted average of exposed soil, snow surface, and surface water
- t_grnd_u: urban ground temperature (= top snow/soil layer)
- t_grnd_r: rural ground temperature (for soil/crop landunits)
- tsl: temperature of near-surface soil layer (layer 1)
- t_soi_10cm: depth-weighted soil temperature in top 10cm
- tsoi17: depth-weighted soil temperature in top 17cm

Ported from inline code in `HydrologyNoDrainage` in `HydrologyNoDrainageMod.F90`.
"""
function update_ground_and_soil_temperatures!(
    t_grnd::Vector{Float64},
    t_grnd_u::Vector{Float64},
    t_grnd_r::Vector{Float64},
    tsl::Vector{Float64},
    t_soi_10cm::Vector{Float64},
    tsoi17::Vector{Float64},
    t_soisno::Matrix{Float64},
    t_h2osfc::Vector{Float64},
    frac_sno_eff::Vector{Float64},
    frac_h2osfc::Vector{Float64},
    snl::Vector{Int},
    dz::Matrix{Float64},
    zi::Matrix{Float64},
    col_landunit::Vector{Int},
    col_itype::Vector{Int},
    lun_urbpoi::AbstractVector{Bool},
    lun_itype::Vector{Int},
    mask_nolake::AbstractVector{Bool},
    bounds::UnitRange{Int},
    nlevsoi::Int,
    nlevsno::Int
)
    # Julia offset for z/dz/t_soisno: Fortran j → Julia j + nlevsno
    # (arrays dimensioned -nlevsno+1:nlevmaxurbgrnd in Fortran)
    joff = nlevsno
    # Julia offset for zi: Fortran j → Julia j + nlevsno + 1
    # (zi dimensioned -nlevsno:nlevmaxurbgrnd in Fortran, one extra entry)
    joff_zi = nlevsno + 1

    # Initialize 10cm and 17cm accumulators for non-urban
    for c in bounds
        mask_nolake[c] || continue
        l = col_landunit[c]
        if !lun_urbpoi[l]
            t_soi_10cm[c] = 0.0
            tsoi17[c] = 0.0
        end
    end

    # Accumulate soil temperature weighted by depth
    for j in 1:nlevsoi
        for c in bounds
            mask_nolake[c] || continue
            l = col_landunit[c]
            if !lun_urbpoi[l]
                # Near-surface soil layer temperature
                if j == 1
                    tsl[c] = t_soisno[c, j + joff]
                end

                # Soil T at top 17 cm (F. Li and S. Levis)
                # Fortran: zi(c,j) → Julia: zi[c, j + joff_zi]
                # Fortran: zi(c,j-1) → Julia: zi[c, (j-1) + joff_zi] = zi[c, j + joff]
                if zi[c, j + joff_zi] <= 0.17
                    fracl = 1.0
                    tsoi17[c] = tsoi17[c] + t_soisno[c, j + joff] * dz[c, j + joff] * fracl
                else
                    if zi[c, j + joff_zi] > 0.17 && zi[c, j + joff] < 0.17
                        fracl = (0.17 - zi[c, j + joff]) / dz[c, j + joff]
                        tsoi17[c] = tsoi17[c] + t_soisno[c, j + joff] * dz[c, j + joff] * fracl
                    end
                end

                # Soil T at top 10 cm
                if zi[c, j + joff_zi] <= 0.1
                    fracl = 1.0
                    t_soi_10cm[c] = t_soi_10cm[c] + t_soisno[c, j + joff] * dz[c, j + joff] * fracl
                else
                    if zi[c, j + joff_zi] > 0.1 && zi[c, j + joff] < 0.1
                        fracl = (0.1 - zi[c, j + joff]) / dz[c, j + joff]
                        t_soi_10cm[c] = t_soi_10cm[c] + t_soisno[c, j + joff] * dz[c, j + joff] * fracl
                    end
                end
            end
        end
    end

    # Ground temperature and finalize 10cm/17cm temperatures
    for c in bounds
        mask_nolake[c] || continue
        l = col_landunit[c]

        # t_grnd is weighted average of exposed soil and snow
        if snl[c] < 0
            t_grnd[c] = frac_sno_eff[c] * t_soisno[c, snl[c] + 1 + joff] +
                         (1.0 - frac_sno_eff[c] - frac_h2osfc[c]) * t_soisno[c, 1 + joff] +
                         frac_h2osfc[c] * t_h2osfc[c]
        else
            t_grnd[c] = (1.0 - frac_h2osfc[c]) * t_soisno[c, 1 + joff] +
                         frac_h2osfc[c] * t_h2osfc[c]
        end

        if lun_urbpoi[l]
            t_grnd_u[c] = t_soisno[c, snl[c] + 1 + joff]
        else
            t_soi_10cm[c] = t_soi_10cm[c] / 0.1
            tsoi17[c] = tsoi17[c] / 0.17
        end

        if lun_itype[l] == ISTSOIL || lun_itype[l] == ISTCROP
            t_grnd_r[c] = t_soisno[c, snl[c] + 1 + joff]
        end
    end

    return nothing
end

# =========================================================================
# update_h2osoi_vol!
# =========================================================================

"""
    update_h2osoi_vol!(h2osoi_vol, h2osoi_liq, h2osoi_ice,
        dz, col_itype, mask_nolake, mask_urban, bounds,
        nlevgrnd, nlevurb, nlevsno)

Update volumetric soil water content from liquid and ice masses.

For non-urban non-wall/roof columns: over nlevgrnd layers.
For urban wall/roof columns: over nlevurb layers.

Ported from inline code in `HydrologyNoDrainage` in `HydrologyNoDrainageMod.F90`.
"""
function update_h2osoi_vol!(
    h2osoi_vol::Matrix{Float64},
    h2osoi_liq::Matrix{Float64},
    h2osoi_ice::Matrix{Float64},
    dz::Matrix{Float64},
    col_itype::Vector{Int},
    mask_nolake::BitVector,
    mask_urban::BitVector,
    bounds::UnitRange{Int},
    nlevgrnd::Int,
    nlevurb::Int,
    nlevsno::Int
)
    joff = nlevsno

    # Non-urban columns (excluding sunwall, shadewall, roof)
    # h2osoi_vol is soil-only (1:nlevmaxurbgrnd), no offset needed
    # h2osoi_liq/ice and dz are snow+soil arrays, need joff offset for soil layer j
    for j in 1:nlevgrnd
        for c in bounds
            mask_nolake[c] || continue
            if col_itype[c] != ICOL_SUNWALL && col_itype[c] != ICOL_SHADEWALL &&
               col_itype[c] != ICOL_ROOF
                h2osoi_vol[c, j] = h2osoi_liq[c, j + joff] / (dz[c, j + joff] * DENH2O) +
                                   h2osoi_ice[c, j + joff] / (dz[c, j + joff] * DENICE)
            end
        end
    end

    # Urban columns (sunwall, shadewall, roof) — only nlevurb layers
    for j in 1:nlevurb
        for c in bounds
            mask_urban[c] || continue
            if col_itype[c] == ICOL_SUNWALL || col_itype[c] == ICOL_SHADEWALL ||
               col_itype[c] == ICOL_ROOF
                h2osoi_vol[c, j] = h2osoi_liq[c, j + joff] / (dz[c, j + joff] * DENH2O) +
                                   h2osoi_ice[c, j + joff] / (dz[c, j + joff] * DENICE)
            end
        end
    end

    return nothing
end

# =========================================================================
# update_soilpsi!
# =========================================================================

"""
    update_soilpsi!(soilpsi, h2osoi_liq, dz, watsat, sucsat, bsw,
        mask_hydrology, bounds, nlevgrnd, nlevsno)

Update soil water potential (soilpsi) in each soil layer [MPa].

Ported from inline code in `HydrologyNoDrainage` in `HydrologyNoDrainageMod.F90`.
"""
function update_soilpsi!(
    soilpsi::Matrix{Float64},
    h2osoi_liq::Matrix{Float64},
    dz::Matrix{Float64},
    watsat::Matrix{Float64},
    sucsat::Matrix{Float64},
    bsw::Matrix{Float64},
    mask_hydrology::BitVector,
    bounds::UnitRange{Int},
    nlevgrnd::Int,
    nlevsno::Int
)
    joff = nlevsno

    # h2osoi_liq and dz are snow+soil arrays, need joff offset for soil layer j
    # soilpsi, watsat, sucsat, bsw are soil-only arrays (1-based)
    for j in 1:nlevgrnd
        for c in bounds
            mask_hydrology[c] || continue

            if h2osoi_liq[c, j + joff] > 0.0
                vwc = h2osoi_liq[c, j + joff] / (dz[c, j + joff] * DENH2O)

                # Limit fractional saturation to avoid numerical crash
                fsattmp = max(vwc / watsat[c, j], 0.001)
                psi = sucsat[c, j] * (-9.8e-6) * (fsattmp)^(-bsw[c, j])  # MPa
                soilpsi[c, j] = min(max(psi, -15.0), 0.0)
            else
                soilpsi[c, j] = -15.0
            end
        end
    end

    return nothing
end

# =========================================================================
# update_smp_l!
# =========================================================================

"""
    update_smp_l!(smp_l, h2osoi_vol, watsat, sucsat, bsw, smpmin,
        mask_hydrology, bounds, nlevgrnd)

Update soil matric potential (smp_l) for history output and ch4Mod [mm].

Ported from inline code in `HydrologyNoDrainage` in `HydrologyNoDrainageMod.F90`.
"""
function update_smp_l!(
    smp_l::Matrix{Float64},
    h2osoi_vol::Matrix{Float64},
    watsat::Matrix{Float64},
    sucsat::Matrix{Float64},
    bsw::Matrix{Float64},
    smpmin::Vector{Float64},
    mask_hydrology::BitVector,
    bounds::UnitRange{Int},
    nlevgrnd::Int
)
    for j in 1:nlevgrnd
        for c in bounds
            mask_hydrology[c] || continue

            s_node = max(h2osoi_vol[c, j] / watsat[c, j], 0.01)
            s_node = min(1.0, s_node)

            smp_l[c, j] = -sucsat[c, j] * s_node^(-bsw[c, j])
            smp_l[c, j] = max(smpmin[c], smp_l[c, j])
        end
    end

    return nothing
end

# =========================================================================
# compute_wf!
# =========================================================================

"""
    compute_wf!(wf, h2osoi_vol, watsat, sucsat, bsw,
        z, dz, mask_hydrology, bounds, nlevgrnd, nlevsno,
        depth_limit)

Compute soil water as fraction of water holding capacity (WHC)
integrated to a given depth_limit [m].

Available soil water = h2osoi_vol - watdry (air-dry water content)
Potentially available soil water (WHC) = watsat - watdry

Returns wf = available / potentially_available for each column.

Ported from inline code in `HydrologyNoDrainage` in `HydrologyNoDrainageMod.F90`.
"""
function compute_wf!(
    wf_out::Vector{Float64},
    h2osoi_vol::Matrix{Float64},
    watsat::Matrix{Float64},
    sucsat::Matrix{Float64},
    bsw::Matrix{Float64},
    z::Matrix{Float64},
    dz::Matrix{Float64},
    mask_hydrology::BitVector,
    bounds::UnitRange{Int},
    nlevgrnd::Int,
    nlevsno::Int,
    depth_limit::Float64
)
    joff = nlevsno

    # Temporary accumulators
    nc = length(wf_out)
    rwat = zeros(Float64, nc)
    swat = zeros(Float64, nc)
    rz   = zeros(Float64, nc)

    for c in bounds
        mask_hydrology[c] || continue
        rwat[c] = 0.0
        swat[c] = 0.0
        rz[c]   = 0.0
    end

    for j in 1:nlevgrnd
        for c in bounds
            mask_hydrology[c] || continue
            if z[c, j + joff] + 0.5 * dz[c, j + joff] <= depth_limit
                watdry = watsat[c, j] * (316230.0 / sucsat[c, j])^(-1.0 / bsw[c, j])
                rwat[c] = rwat[c] + (h2osoi_vol[c, j] - watdry) * dz[c, j + joff]
                swat[c] = swat[c] + (watsat[c, j]      - watdry) * dz[c, j + joff]
                rz[c]   = rz[c] + dz[c, j + joff]
            end
        end
    end

    for c in bounds
        mask_hydrology[c] || continue
        if rz[c] != 0.0
            tsw  = rwat[c] / rz[c]
            stsw = swat[c] / rz[c]
        else
            watdry = watsat[c, 1] * (316230.0 / sucsat[c, 1])^(-1.0 / bsw[c, 1])
            tsw  = h2osoi_vol[c, 1] - watdry
            stsw = watsat[c, 1] - watdry
        end
        wf_out[c] = tsw / stsw
    end

    return nothing
end

# =========================================================================
# update_snow_top_layer_diagnostics!
# =========================================================================

"""
    update_snow_top_layer_diagnostics!(h2osno_top, snw_rds, snot_top,
        dTdz_top, snw_rds_top, sno_liq_top,
        h2osoi_ice, h2osoi_liq, snl,
        mask_snow, mask_nosnow, bounds, nlevsno, spval)

Update top-layer snow diagnostics:
- h2osno_top: mass of snow in top layer
- Zero snow variables and set spval diagnostics for columns without snow.

Ported from inline code in `HydrologyNoDrainage` in `HydrologyNoDrainageMod.F90`.
"""
function update_snow_top_layer_diagnostics!(
    h2osno_top::Vector{Float64},
    snw_rds::Matrix{Float64},
    snot_top::Vector{Float64},
    dTdz_top::Vector{Float64},
    snw_rds_top::Vector{Float64},
    sno_liq_top::Vector{Float64},
    h2osoi_ice::Matrix{Float64},
    h2osoi_liq::Matrix{Float64},
    snl::Vector{Int},
    mask_snow::BitVector,
    mask_nosnow::BitVector,
    bounds::UnitRange{Int},
    nlevsno::Int,
    spval::Float64
)
    # Top-layer diagnostics for columns with snow
    for c in bounds
        mask_snow[c] || continue
        j_top = snl[c] + 1 + nlevsno  # Julia index of top snow layer
        h2osno_top[c] = h2osoi_ice[c, j_top] + h2osoi_liq[c, j_top]
    end

    # Zero variables for columns without snow
    for c in bounds
        mask_nosnow[c] || continue

        h2osno_top[c] = 0.0

        # Zero all snow grain radii layers
        for j in 1:nlevsno
            snw_rds[c, j] = 0.0
        end

        # Top-layer diagnostics set to spval (not averaged in history fields)
        snot_top[c]    = spval
        dTdz_top[c]    = spval
        snw_rds_top[c] = spval
        sno_liq_top[c] = spval
    end

    return nothing
end

# =========================================================================
# calc_and_withdraw_irrigation_fluxes!
# =========================================================================

"""
    calc_and_withdraw_irrigation_fluxes!(soilhydrology, soilstate,
        waterfluxbulk, waterstatebulk, water,
        mask_soil, bounds, nlevsno, nlevgrnd;
        use_groundwater_irrigation=false)

Calculates irrigation withdrawal fluxes and withdraws from groundwater.

In Fortran this calls `irrigation_inst%CalcIrrigationFluxes` and
`WithdrawGroundwaterIrrigation`. Since IrrigationMod is not yet ported,
this is a stub that only performs the groundwater withdrawal step when
`use_groundwater_irrigation` is true.

Ported from `CalcAndWithdrawIrrigationFluxes` in `HydrologyNoDrainageMod.F90`.
"""
function calc_and_withdraw_irrigation_fluxes!(
    waterflux::WaterFluxData,
    waterstate::WaterStateData,
    mask_soil::BitVector,
    bounds::UnitRange{Int},
    nlevsoi::Int,
    dtime::Float64;
    use_groundwater_irrigation::Bool = false
)
    # Stub: CalcIrrigationFluxes requires IrrigationMod (not yet ported)
    # When IrrigationMod is ported, add the irrigation flux calculation here.

    # Groundwater irrigation withdrawal
    if use_groundwater_irrigation
        withdraw_groundwater_irrigation!(
            waterflux, waterstate,
            mask_soil, bounds, nlevsoi, dtime)
    end

    return nothing
end

# =========================================================================
# handle_new_snow!
# =========================================================================

"""
    handle_new_snow!(temperature, waterstatebulk, waterdiagbulk,
        col, lun, mask_nolake, bounds, dtime, nlevsno;
        forc_t, qflx_snow_grnd)

Handle new snow falling on the ground.

Coordinates:
1. Update quantities for new snow (bulk density, diagnostics, state)
2. Remove snow from thawed wetlands
3. Initialize explicit snow pack where warranted

Ported from `HandleNewSnow` in `HydrologyNoDrainageMod.F90`.
"""
function handle_new_snow!(
    temperature::TemperatureData,
    waterstatebulk::WaterStateBulkData,
    waterdiagbulk::WaterDiagnosticBulkData,
    col::ColumnData,
    lun::LandunitData,
    mask_nolake::BitVector,
    bounds::UnitRange{Int},
    dtime::Float64,
    nlevsno::Int;
    forc_t::Vector{Float64},
    forc_wind::Vector{Float64} = Float64[],
    qflx_snow_grnd::Vector{Float64},
    qflx_snow_drain::Vector{Float64} = zeros(length(mask_nolake)),
    int_snow::Vector{Float64} = zeros(length(mask_nolake))
)
    nc = length(mask_nolake)

    # --- 1. Update quantities for new snow ---
    # Step 1a: Compute bulk density of new snow
    bifall = fill(50.0, nc)
    if !isempty(forc_wind)
        new_snow_bulk_density!(bifall, forc_t, forc_wind, col.gridcell,
                               mask_nolake, bounds)
    else
        # Fallback: use temperature-only bulk density
        for c in bounds
            mask_nolake[c] || continue
            if forc_t[c] > TFRZ + 2.0
                bifall[c] = 50.0 + 1.7 * (17.0)^1.5
            elseif forc_t[c] > TFRZ - 15.0
                bifall[c] = 50.0 + 1.7 * (forc_t[c] - TFRZ + 15.0)^1.5
            else
                bifall[c] = 50.0
            end
        end
    end

    # Step 1b: Calculate total H2O in snow (h2osno_total)
    h2osno_total = zeros(nc)
    for c in bounds
        mask_nolake[c] || continue
        h2osno_total[c] = waterstatebulk.ws.h2osno_no_layers_col[c]
        for j in (col.snl[c]+1):0
            jj = j + nlevsno
            h2osno_total[c] += waterstatebulk.ws.h2osoi_ice_col[c, jj] +
                               waterstatebulk.ws.h2osoi_liq_col[c, jj]
        end
    end

    # Step 1c: Update snow diagnostics (snow_depth, frac_sno, int_snow, swe_old)
    lun_itype_col = Vector{Int}(undef, nc)
    urbpoi_col = Vector{Bool}(undef, nc)
    for c in bounds
        mask_nolake[c] || continue
        l = col.landunit[c]
        lun_itype_col[c] = lun.itype[l]
        urbpoi_col[c] = lun.urbpoi[l]
    end

    bulkdiag_new_snow_diagnostics!(
        col.dz, int_snow, waterdiagbulk.swe_old_col,
        waterdiagbulk.frac_sno_col, waterdiagbulk.frac_sno_eff_col,
        waterdiagbulk.snow_depth_col, waterdiagbulk.snomelt_accum_col,
        dtime, lun_itype_col, urbpoi_col, col.snl, bifall,
        h2osno_total,
        waterstatebulk.ws.h2osoi_ice_col, waterstatebulk.ws.h2osoi_liq_col,
        qflx_snow_grnd, qflx_snow_drain,
        mask_nolake, bounds, nlevsno)

    # Step 1d: Add new snow mass to appropriate state variable
    update_state_add_new_snow!(
        waterstatebulk.ws.h2osno_no_layers_col,
        waterstatebulk.ws.h2osoi_ice_col,
        dtime, col.snl, qflx_snow_grnd,
        mask_nolake, bounds, nlevsno)

    # --- 2. Remove snow from thawed wetlands ---
    # Build landunit itype lookup for each column
    lun_itype_col = Vector{Int}(undef, nc)
    for c in bounds
        mask_nolake[c] || continue
        lun_itype_col[c] = lun.itype[col.landunit[c]]
    end

    mask_thawed_wetland = falses(nc)
    build_filter_thawed_wetland_thin_snowpack!(
        mask_thawed_wetland,
        temperature.t_grnd_col, lun_itype_col, col.snl,
        mask_nolake, bounds)

    update_state_remove_snow_thawed_wetlands!(
        waterstatebulk.ws.h2osno_no_layers_col,
        mask_thawed_wetland, bounds)

    bulk_remove_snow_thawed_wetlands!(
        waterdiagbulk.snow_depth_col,
        mask_thawed_wetland, bounds)

    # --- 3. Initialize explicit snow pack ---
    # TODO: Re-enable explicit snow layer creation once all downstream snow
    # physics code paths (compaction, water movement, SNICAR grain aging)
    # are properly initialized from cold start. Currently, creating explicit
    # layers (snl=-1) triggers NaN cascade in uninitialized code paths.
    # For now, snow mass accumulates in h2osno_no_layers_col and diagnostics
    # (snow_depth, frac_sno) are updated by bulkdiag_new_snow_diagnostics!.
    #=
    mask_init_snowpack = falses(nc)
    build_filter_snowpack_initialized!(
        mask_init_snowpack,
        col.snl, lun_itype_col,
        waterdiagbulk.frac_sno_eff_col,
        waterdiagbulk.snow_depth_col,
        qflx_snow_grnd,
        mask_nolake, bounds)

    update_state_initialize_snow_pack!(
        waterstatebulk.ws.h2osno_no_layers_col,
        waterstatebulk.ws.h2osoi_ice_col,
        waterstatebulk.ws.h2osoi_liq_col,
        mask_init_snowpack, bounds, nlevsno)

    bulk_initialize_snow_pack!(
        col.snl, col.zi, col.dz, col.z,
        temperature.t_soisno_col,
        waterdiagbulk.frac_iceold_col,
        waterdiagbulk.snomelt_accum_col,
        forc_t, waterdiagbulk.snow_depth_col,
        mask_init_snowpack, bounds, nlevsno)
    =#

    return nothing
end

# =========================================================================
# hydrology_no_drainage! — Main orchestrator
# =========================================================================

"""
    hydrology_no_drainage!(temperature, soilstate, waterstatebulk, waterfluxbulk,
        waterdiagbulk, soilhydrology, col, lun, water,
        mask_nolake, mask_hydrology, mask_urban, mask_snow, mask_nosnow,
        bounds, dtime, nlevsno, nlevsoi, nlevgrnd, nlevurb)

Main subroutine to execute the calculation of soil/snow hydrology
without drainage.

This is the top-level orchestrator that calls sub-functions for:
  - Snow persistence tracking
  - Snow ice/liquid accumulation
  - Ground temperature calculation
  - Volumetric soil water update
  - Soil water potential update
  - Soil water fraction diagnostics
  - Top-layer snow diagnostics

Note: Many sub-calls from the Fortran original (SnowWater, SoilWater,
WaterTable, etc.) are called separately in the driver. This function
focuses on the inline diagnostic/state updates that follow.

Ported from `HydrologyNoDrainage` in `HydrologyNoDrainageMod.F90`.
"""
function hydrology_no_drainage!(
    temperature::TemperatureData,
    soilstate::SoilStateData,
    waterstatebulk::WaterStateBulkData,
    waterdiagbulk::WaterDiagnosticBulkData,
    col::ColumnData,
    lun::LandunitData,
    mask_nolake::BitVector,
    mask_hydrology::BitVector,
    mask_urban::BitVector,
    mask_snow::BitVector,
    mask_nosnow::BitVector,
    bounds::UnitRange{Int},
    dtime::Float64,
    nlevsno::Int,
    nlevsoi::Int,
    nlevgrnd::Int,
    nlevurb::Int;
    soilhydrology::Union{SoilHydrologyData, Nothing} = nothing
)
    # Note: set_soil_water_fractions! is now called in the driver before physics

    # --- Snow persistence ---
    update_snow_persistence!(
        waterstatebulk.snow_persistence_col,
        dtime, mask_snow, mask_nosnow, bounds)

    # --- Accumulate snowice and snowliq over snow layers ---
    accumulate_snow_ice_liq!(
        waterdiagbulk.snowice_col,
        waterdiagbulk.snowliq_col,
        waterstatebulk.ws.h2osoi_ice_col,
        waterstatebulk.ws.h2osoi_liq_col,
        col.snl, mask_nolake, mask_snow, bounds, nlevsno)

    # --- Column-average snow depth ---
    update_snowdp!(
        waterdiagbulk.snowdp_col,
        waterdiagbulk.snow_depth_col,
        waterdiagbulk.frac_sno_eff_col,
        bounds)

    # --- Snow internal temperature ---
    compute_snow_internal_temperature!(
        temperature.t_sno_mul_mss_col,
        waterstatebulk.ws.h2osoi_ice_col,
        waterstatebulk.ws.h2osoi_liq_col,
        temperature.t_soisno_col,
        col.snl, mask_nolake, mask_snow, bounds, nlevsno, TFRZ)

    # --- Ground and soil temperatures ---
    update_ground_and_soil_temperatures!(
        temperature.t_grnd_col,
        temperature.t_grnd_u_col,
        temperature.t_grnd_r_col,
        temperature.tsl_col,
        temperature.t_soi10cm_col,
        temperature.t_soi17cm_col,
        temperature.t_soisno_col,
        temperature.t_h2osfc_col,
        waterdiagbulk.frac_sno_eff_col,
        waterdiagbulk.frac_h2osfc_col,
        col.snl, col.dz, col.zi,
        col.landunit, col.itype,
        lun.urbpoi, lun.itype,
        mask_nolake, bounds, nlevsoi, nlevsno)

    # --- Volumetric soil water ---
    update_h2osoi_vol!(
        waterstatebulk.ws.h2osoi_vol_col,
        waterstatebulk.ws.h2osoi_liq_col,
        waterstatebulk.ws.h2osoi_ice_col,
        col.dz, col.itype,
        mask_nolake, mask_urban, bounds,
        nlevgrnd, nlevurb, nlevsno)

    # --- Soil water potential ---
    update_soilpsi!(
        soilstate.soilpsi_col,
        waterstatebulk.ws.h2osoi_liq_col,
        col.dz,
        soilstate.watsat_col,
        soilstate.sucsat_col,
        soilstate.bsw_col,
        mask_hydrology, bounds, nlevgrnd, nlevsno)

    # --- Matric potential for history/ch4 ---
    update_smp_l!(
        soilstate.smp_l_col,
        waterstatebulk.ws.h2osoi_vol_col,
        soilstate.watsat_col,
        soilstate.sucsat_col,
        soilstate.bsw_col,
        soilstate.smpmin_col,
        mask_hydrology, bounds, nlevgrnd)

    # --- Soil water as fraction of WHC (top 0.05m) ---
    compute_wf!(
        waterdiagbulk.wf_col,
        waterstatebulk.ws.h2osoi_vol_col,
        soilstate.watsat_col,
        soilstate.sucsat_col,
        soilstate.bsw_col,
        col.z, col.dz,
        mask_hydrology, bounds, nlevgrnd, nlevsno, 0.05)

    # --- Soil water as fraction of WHC (top 0.17m) ---
    compute_wf!(
        waterdiagbulk.wf2_col,
        waterstatebulk.ws.h2osoi_vol_col,
        soilstate.watsat_col,
        soilstate.sucsat_col,
        soilstate.bsw_col,
        col.z, col.dz,
        mask_hydrology, bounds, nlevgrnd, nlevsno, 0.17)

    # --- Top-layer snow diagnostics ---
    update_snow_top_layer_diagnostics!(
        waterdiagbulk.h2osno_top_col,
        waterdiagbulk.snw_rds_col,
        temperature.snot_top_col,
        temperature.dTdz_top_col,
        waterdiagbulk.snw_rds_top_col,
        waterdiagbulk.sno_liq_top_col,
        waterstatebulk.ws.h2osoi_ice_col,
        waterstatebulk.ws.h2osoi_liq_col,
        col.snl,
        mask_snow, mask_nosnow, bounds, nlevsno, SPVAL)

    return nothing
end
