# ==========================================================================
# Ported from: src/biogeophys/SoilTemperatureMod.F90 (2975 lines)
# Calculates snow and soil temperatures including phase change.
#
# Public subroutines:
#   soil_temperature!                   — main driver
#   compute_ground_heat_flux_and_deriv! — ground heat flux and dG/dT
#   compute_heat_diff_flux_and_factor!  — heat diffusion at interfaces
#   set_rhs_vec!                        — RHS vector for temperature solve
#   set_rhs_vec_snow!                   — RHS for snow layers
#   set_rhs_vec_soil!                   — RHS for soil layers
#   set_rhs_vec_standing_surface_water! — RHS for standing water
#   set_matrix!                         — band diagonal matrix setup
#   assemble_matrix_from_submatrices!   — assemble full matrix
#   set_matrix_snow!                    — matrix entries for snow
#   set_matrix_soil!                    — matrix entries for soil
#   set_matrix_standing_surface_water!  — matrix entries for standing water
#
# Private subroutines:
#   soil_therm_prop!     — thermal conductivities and heat capacities
#   phase_change_h2osfc! — phase change of surface water
#   phase_change_beta!   — phase change within snow/soil layers
#   building_hac!        — building heating/cooling (simple method)
# ==========================================================================

const THIN_SFCLAYER = 1.0e-6  # Threshold for thin surface layer

# =========================================================================
# soil_temperature! — Main driver
# =========================================================================
"""
    soil_temperature!(col, lun, patch_data, temperature, energyflux, soilstate,
                      waterstatebulk, waterdiagbulk, waterfluxbulk,
                      solarabs, canopystate, urbanparams, urbantv, atm2lnd_forc_lwrad,
                      mask_nolakec, mask_nolakep, mask_urbanl, mask_urbanc,
                      bounds_col, bounds_lun, bounds_patch, dtime)

Snow and soil temperatures including phase change.
Ported from `SoilTemperature` in `SoilTemperatureMod.F90`.
"""
function soil_temperature!(col::ColumnData, lun::LandunitData, patch_data::PatchData,
                           temperature::TemperatureData, energyflux::EnergyFluxData,
                           soilstate::SoilStateData,
                           waterstatebulk::WaterStateBulkData,
                           waterdiagbulk::WaterDiagnosticBulkData,
                           waterfluxbulk::WaterFluxBulkData,
                           solarabs::SolarAbsorbedData,
                           canopystate::CanopyStateData,
                           urbanparams::UrbanParamsData,
                           urbantv_t_building_max::Vector{Float64},
                           atm2lnd_forc_lwrad::Vector{Float64},
                           mask_nolakec::BitVector, mask_nolakep::BitVector,
                           mask_urbanl::BitVector, mask_urbanc::BitVector,
                           bounds_col::UnitRange{Int}, bounds_lun::UnitRange{Int},
                           bounds_patch::UnitRange{Int},
                           dtime::Float64)

    nlevsno = varpar.nlevsno
    nlevgrnd = varpar.nlevgrnd
    nlevurb = varpar.nlevurb
    nlevmaxurbgrnd = varpar.nlevmaxurbgrnd
    nlevsoi = varpar.nlevsoi
    joff = nlevsno  # offset: Fortran j → Julia j+joff

    nc = length(bounds_col)
    nband = 5

    # Local arrays
    jtop = fill(-9999, nc)
    jbot = fill(0, nc)
    cv = zeros(nc, nlevsno + nlevmaxurbgrnd)
    tk = zeros(nc, nlevsno + nlevmaxurbgrnd)
    fn = zeros(nc, nlevsno + nlevmaxurbgrnd)
    fn1 = zeros(nc, nlevsno + nlevmaxurbgrnd)
    tk_h2osfc = fill(NaN, nc)
    hs_top = zeros(nc)
    dhsdT = zeros(nc)
    hs_soil = zeros(nc)
    hs_top_snow = zeros(nc)
    hs_h2osfc = zeros(nc)
    sabg_lyr_col = zeros(nc, nlevsno + 1)  # -nlevsno+1:1
    fn_h2osfc = zeros(nc)
    dz_h2osfc = zeros(nc)
    c_h2osfc = zeros(nc)
    cool_on = falses(length(bounds_lun))
    heat_on = falses(length(bounds_lun))

    # Band matrix and vectors (Fortran: -nlevsno:nlevmaxurbgrnd → nlevsno+1+nlevmaxurbgrnd levels)
    nlev_total = nlevsno + 1 + nlevmaxurbgrnd  # -nlevsno to nlevmaxurbgrnd
    bmatrix = zeros(nc, nband, nlev_total)
    tvector = fill(NaN, nc, nlev_total)
    rvector = fill(NaN, nc, nlev_total)

    # Aliases for temperature fields
    snl = col.snl
    t_soisno = temperature.t_soisno_col
    t_h2osfc = temperature.t_h2osfc_col
    t_grnd = temperature.t_grnd_col
    t_building = temperature.t_building_lun
    t_roof_inner = temperature.t_roof_inner_lun
    t_sunw_inner = temperature.t_sunw_inner_lun
    t_shdw_inner = temperature.t_shdw_inner_lun
    tssbef = temperature.t_ssbef_col
    xmf = temperature.xmf_col
    xmf_h2osfc_arr = temperature.xmf_h2osfc_col
    fact = temperature.fact_col
    c_h2osfc_out = temperature.c_h2osfc_col
    imelt = temperature.imelt_col
    emg = temperature.emg_col

    frac_sno_eff = waterdiagbulk.frac_sno_eff_col
    frac_h2osfc = waterdiagbulk.frac_h2osfc_col
    snow_depth = waterdiagbulk.snow_depth_col
    h2osfc = waterstatebulk.ws.h2osfc_col
    excess_ice = waterstatebulk.ws.excess_ice_col

    htvp = energyflux.htvp_col
    eflx_bot = energyflux.eflx_bot_col
    eflx_fgr12 = energyflux.eflx_fgr12_col
    eflx_fgr = energyflux.eflx_fgr_col
    eflx_building_heat_errsoi = energyflux.eflx_building_heat_errsoi_col
    eflx_urban_ac_col = energyflux.eflx_urban_ac_col
    eflx_urban_heat_col = energyflux.eflx_urban_heat_col

    # Simple building temperature
    if is_simple_build_temp()
        building_hac!(lun, temperature, urbanparams, urbantv_t_building_max,
                      mask_urbanl, bounds_lun, cool_on, heat_on)
    end

    # Set jtop and jbot
    for c in bounds_col
        mask_nolakec[c] || continue
        jtop[c] = snl[c]
        if col.itype[c] == ICOL_SUNWALL || col.itype[c] == ICOL_SHADEWALL || col.itype[c] == ICOL_ROOF
            jbot[c] = nlevurb
        else
            jbot[c] = nlevgrnd
        end
    end

    # Excess ice vertical coordinate adjustment (save originals)
    dz_0 = nothing
    zi_0 = nothing
    z_0 = nothing
    if varctl.use_excess_ice
        dz_0 = copy(col.dz)
        zi_0 = copy(col.zi)
        z_0 = copy(col.z)
        for c in bounds_col
            mask_nolakec[c] || continue
            l = col.landunit[c]
            if lun.itype[l] == ISTSOIL || lun.itype[l] == ISTCROP
                for j in 1:nlevmaxurbgrnd
                    col.dz[c, j + joff] += excess_ice[c, j] / DENICE
                end
                for j in 1:nlevmaxurbgrnd
                    col.zi[c, j + joff + 1] = col.zi[c, joff + 1] # zi(c,0) offset
                    # Recompute: zi(c,j) = zi(c,j) + sum(excess_ice(c,1:j))/denice
                    ei_sum = 0.0
                    for jj in 1:j
                        ei_sum += excess_ice[c, jj]
                    end
                    col.zi[c, j + joff + 1] = zi_0[c, j + joff + 1] + ei_sum / DENICE
                    col.z[c, j + joff] = (col.zi[c, j - 1 + joff + 1] + col.zi[c, j + joff + 1]) * 0.5
                end
            end
        end
    end

    # Thermal conductivity and heat capacity
    soil_therm_prop!(col, lun, urbanparams, temperature, waterstatebulk, waterdiagbulk, soilstate,
                     mask_nolakec, mask_urbanc, bounds_col,
                     tk, cv, tk_h2osfc)

    # Ground heat flux
    compute_ground_heat_flux_and_deriv!(col, lun, patch_data, temperature, energyflux,
                                        solarabs, canopystate, waterdiagbulk, waterfluxbulk,
                                        urbanparams, atm2lnd_forc_lwrad,
                                        mask_nolakec, mask_nolakep,
                                        bounds_col, bounds_patch,
                                        hs_h2osfc, hs_top_snow, hs_soil, hs_top, dhsdT,
                                        sabg_lyr_col)

    # Heat diffusion and factor
    compute_heat_diff_flux_and_factor!(col, lun, temperature, energyflux, urbanparams,
                                       mask_nolakec, bounds_col, dtime,
                                       tk, cv, fn, fact)

    # Thermal properties of h2osfc
    for c in bounds_col
        mask_nolakec[c] || continue
        if h2osfc[c] > THIN_SFCLAYER && frac_h2osfc[c] > THIN_SFCLAYER
            c_h2osfc_out[c] = max(THIN_SFCLAYER, CPLIQ * h2osfc[c] / frac_h2osfc[c])
            dz_h2osfc[c] = max(THIN_SFCLAYER, 1.0e-3 * h2osfc[c] / frac_h2osfc[c])
        else
            c_h2osfc_out[c] = THIN_SFCLAYER
            dz_h2osfc[c] = THIN_SFCLAYER
        end
    end

    # Set up RHS vector
    set_rhs_vec!(col, lun, temperature, waterdiagbulk,
                 mask_nolakec, bounds_col, dtime,
                 hs_h2osfc, hs_top_snow, hs_soil, hs_top, dhsdT,
                 sabg_lyr_col, tk, tk_h2osfc, fact, fn,
                 c_h2osfc_out, dz_h2osfc,
                 rvector)

    # Set up band diagonal matrix
    set_matrix!(col, lun, waterdiagbulk, mask_nolakec, bounds_col, dtime, nband,
                dhsdT, tk, tk_h2osfc, fact, c_h2osfc_out, dz_h2osfc,
                bmatrix)

    # Initialize temperature vector
    for c in bounds_col
        mask_nolakec[c] || continue
        # Snow layers: tvector index = j-1 in Fortran → j-1 + nlevsno + 1 in Julia
        for j in (snl[c] + 1):0
            # Fortran tvector(c, j-1) = t_soisno(c, j)
            # j-1 maps to Julia index: (j-1) + nlevsno + 1
            tvector[c, j - 1 + nlevsno + 1] = t_soisno[c, j + joff]
        end
        # Surface water: tvector(c, 0) → index nlevsno+1
        tvector[c, nlevsno + 1] = t_h2osfc[c]
        # Soil layers: tvector(c, 1:nlevmaxurbgrnd)
        for j in 1:nlevmaxurbgrnd
            tvector[c, j + nlevsno + 1] = t_soisno[c, j + joff]
        end
    end

    # Solve the band diagonal system
    # Convert jtop/jbot to indices in the tvector space
    # Fortran: jtop = snl(c), so system runs from snl(c) to jbot(c) in the -nlevsno:nlevmaxurbgrnd space
    # Julia tvector index: Fortran_idx + nlevsno + 1
    kl = div(nband - 1, 2)
    ku = kl
    m_ab = 2 * kl + ku + 1

    for c in bounds_col
        mask_nolakec[c] || continue
        jt = jtop[c] + nlevsno + 1  # map to Julia 1-based
        jb = jbot[c] + nlevsno + 1
        n = jb - jt + 1
        if n <= 0
            continue
        end

        ab = zeros(m_ab, n)

        for jj in 1:n
            for band_idx in 1:nband
                row_offset = band_idx - (kl + 1)
                i = jj + row_offset
                if 1 <= i <= n
                    ab_row = kl + ku + 1 + i - jj
                    # Fortran BandDiagonal shifts source index by row_offset for off-diagonals:
                    # diagonal reads from jtop+j-1, sub/super-diagonals shift ±1, ±2
                    src_idx = jt + jj - 1 + row_offset
                    ab[ab_row, jj] = bmatrix[c, band_idx, src_idx]
                end
            end
        end

        rhs = rvector[c, jt:jb]
        ipiv = zeros(Int64, n)
        rhs_mat = reshape(rhs, n, 1)

        try
            (ab, ipiv) = LinearAlgebra.LAPACK.gbtrf!(kl, ku, n, ab)
            LinearAlgebra.LAPACK.gbtrs!('N', kl, ku, n, ab, ipiv, rhs_mat)
            tvector[c, jt:jb] .= rhs
        catch e
            e isa LinearAlgebra.LAPACKException || rethrow()
            # Singular matrix: leave tvector unchanged (preserves previous temperatures)
        end
    end

    # Return temperatures from tvector
    for c in bounds_col
        mask_nolakec[c] || continue
        for j in (snl[c] + 1):0
            t_soisno[c, j + joff] = tvector[c, j - 1 + nlevsno + 1]
        end
        for j in 1:nlevmaxurbgrnd
            t_soisno[c, j + joff] = tvector[c, j + nlevsno + 1]
        end
        if frac_h2osfc[c] == 0.0
            t_h2osfc[c] = t_soisno[c, 1 + joff]
        else
            t_h2osfc[c] = tvector[c, nlevsno + 1]
        end
    end

    # Compute fn1 for melting/freezing
    for j in (-nlevsno + 1):nlevmaxurbgrnd
        for c in bounds_col
            mask_nolakec[c] || continue
            l = col.landunit[c]
            jj = j + joff  # Julia index for t_soisno, z, etc.

            if (col.itype[c] == ICOL_SUNWALL || col.itype[c] == ICOL_SHADEWALL ||
                col.itype[c] == ICOL_ROOF) && j <= nlevurb
                if j >= snl[c] + 1
                    if j <= nlevurb - 1
                        fn1[c, jj] = tk[c, jj] * (t_soisno[c, jj + 1] - t_soisno[c, jj]) /
                                     (col.z[c, jj + 1] - col.z[c, jj])
                    elseif j == nlevurb
                        if is_simple_build_temp()
                            fn1[c, jj] = tk[c, jj] * (t_building[l] - t_soisno[c, jj]) /
                                         (col.zi[c, jj + 1] - col.z[c, jj])
                            fn[c, jj] = tk[c, jj] * (t_building[l] - tssbef[c, jj]) /
                                        (col.zi[c, jj + 1] - col.z[c, jj])
                        else
                            if col.itype[c] == ICOL_SUNWALL
                                fn1[c, jj] = tk[c, jj] * (t_sunw_inner[l] - t_soisno[c, jj]) /
                                             (col.zi[c, jj + 1] - col.z[c, jj])
                                fn[c, jj] = tk[c, jj] * (t_sunw_inner[l] - tssbef[c, jj]) /
                                            (col.zi[c, jj + 1] - col.z[c, jj])
                            elseif col.itype[c] == ICOL_SHADEWALL
                                fn1[c, jj] = tk[c, jj] * (t_shdw_inner[l] - t_soisno[c, jj]) /
                                             (col.zi[c, jj + 1] - col.z[c, jj])
                                fn[c, jj] = tk[c, jj] * (t_shdw_inner[l] - tssbef[c, jj]) /
                                            (col.zi[c, jj + 1] - col.z[c, jj])
                            elseif col.itype[c] == ICOL_ROOF
                                fn1[c, jj] = tk[c, jj] * (t_roof_inner[l] - t_soisno[c, jj]) /
                                             (col.zi[c, jj + 1] - col.z[c, jj])
                                fn[c, jj] = tk[c, jj] * (t_roof_inner[l] - tssbef[c, jj]) /
                                            (col.zi[c, jj + 1] - col.z[c, jj])
                            end
                        end
                    end
                end
            elseif col.itype[c] != ICOL_SUNWALL && col.itype[c] != ICOL_SHADEWALL &&
                   col.itype[c] != ICOL_ROOF
                if j >= snl[c] + 1
                    if j <= nlevgrnd - 1
                        fn1[c, jj] = tk[c, jj] * (t_soisno[c, jj + 1] - t_soisno[c, jj]) /
                                     (col.z[c, jj + 1] - col.z[c, jj])
                    elseif j == nlevgrnd
                        fn1[c, jj] = 0.0
                    end
                end
            end
        end
    end

    # Urban building heat flux
    for c in bounds_col
        mask_nolakec[c] || continue
        l = col.landunit[c]
        if lun.urbpoi[l]
            if col.itype[c] == ICOL_SUNWALL || col.itype[c] == ICOL_SHADEWALL || col.itype[c] == ICOL_ROOF
                eflx_building_heat_errsoi[c] = CNFAC * fn[c, nlevurb + joff] +
                    (1.0 - CNFAC) * fn1[c, nlevurb + joff]
            else
                eflx_building_heat_errsoi[c] = 0.0
            end
            if is_simple_build_temp()
                if cool_on[l]
                    eflx_urban_ac_col[c] = abs(eflx_building_heat_errsoi[c])
                    eflx_urban_heat_col[c] = 0.0
                elseif heat_on[l]
                    eflx_urban_ac_col[c] = 0.0
                    eflx_urban_heat_col[c] = abs(eflx_building_heat_errsoi[c])
                else
                    eflx_urban_ac_col[c] = 0.0
                    eflx_urban_heat_col[c] = 0.0
                end
            end
        end
    end

    # Phase change of h2osfc
    for c in bounds_col
        mask_nolakec[c] || continue
        xmf_h2osfc_arr[c] = 0.0
    end

    phase_change_h2osfc!(col, temperature, energyflux, waterstatebulk, waterdiagbulk, waterfluxbulk,
                         mask_nolakec, bounds_col, dtime, dhsdT)

    # Phase change within snow/soil
    phase_change_beta!(col, lun, temperature, energyflux, soilstate,
                       waterstatebulk, waterdiagbulk, waterfluxbulk,
                       mask_nolakec, bounds_col, dtime, dhsdT)

    # Restore excess ice vertical coordinates
    if varctl.use_excess_ice
        for c in bounds_col
            mask_nolakec[c] || continue
            l = col.landunit[c]
            if lun.itype[l] == ISTSOIL || lun.itype[l] == ISTCROP
                for j in 1:nlevmaxurbgrnd
                    col.dz[c, j + joff] = dz_0[c, j + joff]
                    col.zi[c, j + joff + 1] = zi_0[c, j + joff + 1]
                    col.z[c, j + joff] = z_0[c, j + joff]
                end
            end
        end
    end

    # Ground temperature
    for c in bounds_col
        mask_nolakec[c] || continue
        if snl[c] < 0
            if frac_h2osfc[c] != 0.0
                t_grnd[c] = frac_sno_eff[c] * t_soisno[c, snl[c] + 1 + joff] +
                    (1.0 - frac_sno_eff[c] - frac_h2osfc[c]) * t_soisno[c, 1 + joff] +
                    frac_h2osfc[c] * t_h2osfc[c]
            else
                t_grnd[c] = frac_sno_eff[c] * t_soisno[c, snl[c] + 1 + joff] +
                    (1.0 - frac_sno_eff[c]) * t_soisno[c, 1 + joff]
            end
        else
            if frac_h2osfc[c] != 0.0
                t_grnd[c] = (1.0 - frac_h2osfc[c]) * t_soisno[c, 1 + joff] +
                    frac_h2osfc[c] * t_h2osfc[c]
            else
                t_grnd[c] = t_soisno[c, 1 + joff]
            end
        end
    end

    # Soil heat content
    for c in bounds_col
        mask_nolakec[c] || continue
        eflx_fgr12[c] = 0.0
    end

    for j in (-nlevsno + 1):nlevgrnd
        for c in bounds_col
            mask_nolakec[c] || continue
            l = col.landunit[c]
            jj = j + joff

            if j == 1
                eflx_fgr12[c] = -CNFAC * fn[c, 1 + joff] - (1.0 - CNFAC) * fn1[c, 1 + joff]
            end
            if j > 0 && j < nlevgrnd && (lun.itype[l] == ISTSOIL || lun.itype[l] == ISTCROP)
                eflx_fgr[c, j] = -CNFAC * fn[c, jj] - (1.0 - CNFAC) * fn1[c, jj]
            elseif j == nlevgrnd && (lun.itype[l] == ISTSOIL || lun.itype[l] == ISTCROP)
                eflx_fgr[c, j] = 0.0
            end
        end
    end

    return nothing
end

# =========================================================================
# soil_therm_prop! — Thermal conductivities and heat capacities
# =========================================================================
function soil_therm_prop!(col::ColumnData, lun::LandunitData,
                          urbanparams::UrbanParamsData, temperature::TemperatureData,
                          waterstatebulk::WaterStateBulkData,
                          waterdiagbulk::WaterDiagnosticBulkData,
                          soilstate::SoilStateData,
                          mask_nolakec::BitVector, mask_urbanc::BitVector,
                          bounds_col::UnitRange{Int},
                          tk_out::Matrix{Float64}, cv_out::Matrix{Float64},
                          tk_h2osfc_out::Vector{Float64})
    nlevsno = varpar.nlevsno
    nlevgrnd = varpar.nlevgrnd
    nlevurb = varpar.nlevurb
    nlevmaxurbgrnd = varpar.nlevmaxurbgrnd
    nlevsoi = varpar.nlevsoi
    joff = nlevsno

    t_soisno = temperature.t_soisno_col
    h2osoi_liq = waterstatebulk.ws.h2osoi_liq_col
    h2osoi_ice = waterstatebulk.ws.h2osoi_ice_col
    excess_ice = waterstatebulk.ws.excess_ice_col
    h2osno_no_layers = waterstatebulk.ws.h2osno_no_layers_col
    h2osfc = waterstatebulk.ws.h2osfc_col
    frac_sno = waterdiagbulk.frac_sno_eff_col
    bw = waterdiagbulk.bw_col
    thk = soilstate.thk_col
    tkmg = soilstate.tkmg_col
    tkdry = soilstate.tkdry_col
    csol = soilstate.csol_col
    watsat = soilstate.watsat_col
    tksatu = soilstate.tksatu_col

    # Thermal conductivity of soil (Farouki 1981)
    for j in (-nlevsno + 1):nlevgrnd
        for c in bounds_col
            mask_nolakec[c] || continue
            jj = j + joff

            if j >= 1
                l = col.landunit[c]
                nlev_improad_l = urbanparams.nlev_improad[l]

                if (lun.itype[l] != ISTWET && lun.itype[l] != ISTICE &&
                    col.itype[c] != ICOL_SUNWALL && col.itype[c] != ICOL_SHADEWALL &&
                    col.itype[c] != ICOL_ROOF && col.itype[c] != ICOL_ROAD_IMPERV) ||
                   (col.itype[c] == ICOL_ROAD_IMPERV && j > nlev_improad_l)

                    satw = (h2osoi_liq[c, jj] / DENH2O + h2osoi_ice[c, jj] / DENICE +
                            excess_ice[c, j] / DENICE) / (col.dz[c, jj] * watsat[c, j])
                    satw = min(1.0, satw)
                    if satw > 1.0e-7
                        if t_soisno[c, jj] >= TFRZ
                            dke = max(0.0, log10(satw) + 1.0)
                        else
                            dke = satw
                        end
                        fl_denom = h2osoi_liq[c, jj] / (DENH2O * col.dz[c, jj]) +
                                   h2osoi_ice[c, jj] / (DENICE * col.dz[c, jj]) +
                                   excess_ice[c, j] / (DENICE * col.dz[c, jj])
                        fl = (h2osoi_liq[c, jj] / (DENH2O * col.dz[c, jj])) / fl_denom
                        dksat = tkmg[c, j] * TKWAT^(fl * watsat[c, j]) * TKICE^((1.0 - fl) * watsat[c, j])
                        thk[c, jj] = dke * dksat + (1.0 - dke) * tkdry[c, j]
                    else
                        thk[c, jj] = tkdry[c, j]
                    end
                    if j > col.nbedrock[c]
                        thk[c, jj] = THK_BEDROCK
                    end
                elseif lun.itype[l] == ISTICE
                    thk[c, jj] = TKWAT
                    if t_soisno[c, jj] < TFRZ
                        thk[c, jj] = TKICE
                    end
                elseif lun.itype[l] == ISTWET
                    if j > nlevsoi
                        thk[c, jj] = THK_BEDROCK
                    else
                        thk[c, jj] = TKWAT
                        if t_soisno[c, jj] < TFRZ
                            thk[c, jj] = TKICE
                        end
                    end
                end
            end

            # Thermal conductivity of snow
            if col.snl[c] + 1 < 1 && j >= col.snl[c] + 1 && j <= 0
                denom = max(frac_sno[c], 1.0e-6) * max(col.dz[c, jj], 1.0e-6)
                bw[c, jj] = (h2osoi_ice[c, jj] + h2osoi_liq[c, jj]) / denom
                if varctl.snow_thermal_cond_method == "Jordan1991"
                    thk[c, jj] = TKAIR + (7.75e-5 * bw[c, jj] + 1.105e-6 * bw[c, jj]^2) * (TKICE - TKAIR)
                elseif varctl.snow_thermal_cond_method == "Sturm1997"
                    if bw[c, jj] <= 156
                        thk[c, jj] = 0.023 + 0.234 * (bw[c, jj] / 1000)
                    else
                        thk[c, jj] = 0.138 - 1.01 * (bw[c, jj] / 1000) +
                                     3.233 * (bw[c, jj] / 1000)^2
                    end
                else
                    error("Unknown snow_thermal_cond_method: $(varctl.snow_thermal_cond_method)")
                end
            end
        end
    end

    # Urban columns
    for j in 1:nlevurb
        for c in bounds_col
            mask_urbanc[c] || continue
            l = col.landunit[c]
            jj = j + joff
            if col.itype[c] == ICOL_SUNWALL || col.itype[c] == ICOL_SHADEWALL
                thk[c, jj] = urbanparams.tk_wall[l, j]
            elseif col.itype[c] == ICOL_ROOF
                thk[c, jj] = urbanparams.tk_roof[l, j]
            elseif col.itype[c] == ICOL_ROAD_IMPERV && j <= urbanparams.nlev_improad[l]
                thk[c, jj] = urbanparams.tk_improad[l, j]
            end
        end
    end

    # Thermal conductivity at layer interface
    for j in (-nlevsno + 1):nlevmaxurbgrnd
        for c in bounds_col
            mask_nolakec[c] || continue
            jj = j + joff
            if (col.itype[c] == ICOL_SUNWALL || col.itype[c] == ICOL_SHADEWALL ||
                col.itype[c] == ICOL_ROOF) && j <= nlevurb
                if j >= col.snl[c] + 1 && j <= nlevurb - 1
                    tk_out[c, jj] = thk[c, jj] * thk[c, jj + 1] * (col.z[c, jj + 1] - col.z[c, jj]) /
                        (thk[c, jj] * (col.z[c, jj + 1] - col.zi[c, jj + 1]) +
                         thk[c, jj + 1] * (col.zi[c, jj + 1] - col.z[c, jj]))
                elseif j == nlevurb
                    tk_out[c, jj] = thk[c, jj]
                end
            elseif col.itype[c] != ICOL_SUNWALL && col.itype[c] != ICOL_SHADEWALL &&
                   col.itype[c] != ICOL_ROOF
                if j >= col.snl[c] + 1 && j <= nlevgrnd - 1
                    tk_out[c, jj] = thk[c, jj] * thk[c, jj + 1] * (col.z[c, jj + 1] - col.z[c, jj]) /
                        (thk[c, jj] * (col.z[c, jj + 1] - col.zi[c, jj + 1]) +
                         thk[c, jj + 1] * (col.zi[c, jj + 1] - col.z[c, jj]))
                elseif j == nlevgrnd
                    tk_out[c, jj] = 0.0
                end
            end
        end
    end

    # Thermal conductivity of h2osfc
    for c in bounds_col
        mask_nolakec[c] || continue
        zh2osfc = 1.0e-3 * (0.5 * h2osfc[c])
        tk_h2osfc_out[c] = TKWAT * thk[c, 1 + joff] * (col.z[c, 1 + joff] + zh2osfc) /
            (TKWAT * col.z[c, 1 + joff] + thk[c, 1 + joff] * zh2osfc)
    end

    # Soil heat capacity (de Vries 1963)
    for j in 1:nlevgrnd
        for c in bounds_col
            mask_nolakec[c] || continue
            l = col.landunit[c]
            jj = j + joff
            nlev_improad_l = urbanparams.nlev_improad[l]

            if (lun.itype[l] != ISTWET && lun.itype[l] != ISTICE &&
                col.itype[c] != ICOL_SUNWALL && col.itype[c] != ICOL_SHADEWALL &&
                col.itype[c] != ICOL_ROOF && col.itype[c] != ICOL_ROAD_IMPERV) ||
               (col.itype[c] == ICOL_ROAD_IMPERV && j > nlev_improad_l)
                cv_out[c, jj] = csol[c, j] * (1.0 - watsat[c, j]) * col.dz[c, jj] +
                    (h2osoi_ice[c, jj] * CPICE + h2osoi_liq[c, jj] * CPLIQ) +
                    excess_ice[c, j] * CPICE
                if j > col.nbedrock[c]
                    cv_out[c, jj] = CSOL_BEDROCK * col.dz[c, jj]
                end
            elseif lun.itype[l] == ISTWET
                cv_out[c, jj] = h2osoi_ice[c, jj] * CPICE + h2osoi_liq[c, jj] * CPLIQ
                if j > col.nbedrock[c]
                    cv_out[c, jj] = CSOL_BEDROCK * col.dz[c, jj]
                end
            elseif lun.itype[l] == ISTICE
                cv_out[c, jj] = h2osoi_ice[c, jj] * CPICE + h2osoi_liq[c, jj] * CPLIQ
            end
        end
    end

    # Urban heat capacity
    for j in 1:nlevurb
        for c in bounds_col
            mask_urbanc[c] || continue
            l = col.landunit[c]
            jj = j + joff
            if col.itype[c] == ICOL_SUNWALL || col.itype[c] == ICOL_SHADEWALL
                cv_out[c, jj] = urbanparams.cv_wall[l, j] * col.dz[c, jj]
            elseif col.itype[c] == ICOL_ROOF
                cv_out[c, jj] = urbanparams.cv_roof[l, j] * col.dz[c, jj]
            elseif col.itype[c] == ICOL_ROAD_IMPERV && j <= urbanparams.nlev_improad[l]
                cv_out[c, jj] = urbanparams.cv_improad[l, j] * col.dz[c, jj]
            end
        end
    end

    # Unresolved snow on top layer
    for c in bounds_col
        mask_nolakec[c] || continue
        if h2osno_no_layers[c] > 0.0
            cv_out[c, 1 + joff] += CPICE * h2osno_no_layers[c]
        end
    end

    # Snow heat capacity
    for j in (-nlevsno + 1):0
        for c in bounds_col
            mask_nolakec[c] || continue
            jj = j + joff
            if col.snl[c] + 1 < 1 && j >= col.snl[c] + 1
                if frac_sno[c] > 0.0
                    cv_out[c, jj] = max(THIN_SFCLAYER,
                        (CPLIQ * h2osoi_liq[c, jj] + CPICE * h2osoi_ice[c, jj]) / frac_sno[c])
                else
                    cv_out[c, jj] = THIN_SFCLAYER
                end
            end
        end
    end

    return nothing
end

# =========================================================================
# compute_ground_heat_flux_and_deriv!
# =========================================================================
function compute_ground_heat_flux_and_deriv!(
        col::ColumnData, lun::LandunitData, patch_data::PatchData,
        temperature::TemperatureData, energyflux::EnergyFluxData,
        solarabs::SolarAbsorbedData, canopystate::CanopyStateData,
        waterdiagbulk::WaterDiagnosticBulkData, waterfluxbulk::WaterFluxBulkData,
        urbanparams::UrbanParamsData, forc_lwrad::Vector{Float64},
        mask_nolakec::BitVector, mask_nolakep::BitVector,
        bounds_col::UnitRange{Int}, bounds_patch::UnitRange{Int},
        hs_h2osfc::Vector{Float64}, hs_top_snow::Vector{Float64},
        hs_soil::Vector{Float64}, hs_top::Vector{Float64},
        dhsdT::Vector{Float64}, sabg_lyr_col::Matrix{Float64})

    nlevsno = varpar.nlevsno
    joff = nlevsno

    snl = col.snl
    t_grnd = temperature.t_grnd_col
    t_soisno = temperature.t_soisno_col
    t_h2osfc = temperature.t_h2osfc_col
    emg = temperature.emg_col
    frac_sno_eff = waterdiagbulk.frac_sno_eff_col
    htvp = energyflux.htvp_col
    eflx_gnet = energyflux.eflx_gnet_patch
    dgnetdT = energyflux.dgnetdT_patch

    nc = length(bounds_col)
    lwrad_emit = zeros(nc)
    dlwrad_emit = zeros(nc)
    lwrad_emit_snow = zeros(nc)
    lwrad_emit_soil = zeros(nc)
    lwrad_emit_h2osfc_arr = zeros(nc)
    hs = zeros(nc)

    for c in bounds_col
        mask_nolakec[c] || continue
        lwrad_emit[c] = emg[c] * SB * t_grnd[c]^4
        dlwrad_emit[c] = 4.0 * emg[c] * SB * t_grnd[c]^3
        lwrad_emit_snow[c] = emg[c] * SB * t_soisno[c, snl[c] + 1 + joff]^4
        lwrad_emit_soil[c] = emg[c] * SB * t_soisno[c, 1 + joff]^4
        lwrad_emit_h2osfc_arr[c] = emg[c] * SB * t_h2osfc[c]^4
    end

    hs_soil .= 0.0
    hs_h2osfc .= 0.0
    hs .= 0.0
    dhsdT .= 0.0

    for p in bounds_patch
        mask_nolakep[p] || continue
        c = patch_data.column[p]
        l = patch_data.landunit[p]

        frac_veg_nosno = canopystate.frac_veg_nosno_patch[p]
        sabg = solarabs.sabg_patch[p]
        sabg_soil = solarabs.sabg_soil_patch[p]
        sabg_snow = solarabs.sabg_snow_patch[p]
        dlrad = energyflux.dlrad_patch[p]
        cgrnd = energyflux.cgrnd_patch[p]
        eflx_sh_grnd = energyflux.eflx_sh_grnd_patch[p]
        eflx_sh_snow = energyflux.eflx_sh_snow_patch[p]
        eflx_sh_soil = energyflux.eflx_sh_soil_patch[p]
        eflx_sh_h2osfc = energyflux.eflx_sh_h2osfc_patch[p]
        qflx_evap_soi = waterfluxbulk.wf.qflx_evap_soi_patch[p]
        qflx_ev_snow = waterfluxbulk.qflx_ev_snow_patch[p]
        qflx_ev_soil = waterfluxbulk.qflx_ev_soil_patch[p]
        qflx_ev_h2osfc = waterfluxbulk.qflx_ev_h2osfc_patch[p]

        if !lun.urbpoi[l]
            eflx_gnet[p] = sabg + dlrad +
                (1.0 - frac_veg_nosno) * emg[c] * forc_lwrad[c] - lwrad_emit[c] -
                (eflx_sh_grnd + qflx_evap_soi * htvp[c])

            solarabs.sabg_chk_patch[p] = frac_sno_eff[c] * sabg_snow +
                (1.0 - frac_sno_eff[c]) * sabg_soil

            eflx_gnet_snow = sabg_snow + dlrad +
                (1.0 - frac_veg_nosno) * emg[c] * forc_lwrad[c] - lwrad_emit_snow[c] -
                (eflx_sh_snow + qflx_ev_snow * htvp[c])

            eflx_gnet_soil = sabg_soil + dlrad +
                (1.0 - frac_veg_nosno) * emg[c] * forc_lwrad[c] - lwrad_emit_soil[c] -
                (eflx_sh_soil + qflx_ev_soil * htvp[c])

            eflx_gnet_h2osfc = sabg_soil + dlrad +
                (1.0 - frac_veg_nosno) * emg[c] * forc_lwrad[c] - lwrad_emit_h2osfc_arr[c] -
                (eflx_sh_h2osfc + qflx_ev_h2osfc * htvp[c])
        else
            eflx_lwrad_net = energyflux.eflx_lwrad_net_patch[p]
            qflx_tran_veg = waterfluxbulk.wf.qflx_tran_veg_patch[p]

            if col.itype[c] == ICOL_ROAD_PERV || col.itype[c] == ICOL_ROAD_IMPERV
                energyflux.eflx_wasteheat_patch[p] = energyflux.eflx_wasteheat_lun[l] /
                    (1.0 - lun.wtlunit_roof[l])
                if is_simple_build_temp()
                    energyflux.eflx_ventilation_patch[p] = 0.0
                elseif is_prog_build_temp()
                    energyflux.eflx_ventilation_patch[p] = energyflux.eflx_ventilation_lun[l] /
                        (1.0 - lun.wtlunit_roof[l])
                end
                energyflux.eflx_heat_from_ac_patch[p] = energyflux.eflx_heat_from_ac_lun[l] /
                    (1.0 - lun.wtlunit_roof[l])
                energyflux.eflx_traffic_patch[p] = energyflux.eflx_traffic_lun[l] /
                    (1.0 - lun.wtlunit_roof[l])
            else
                energyflux.eflx_wasteheat_patch[p] = 0.0
                energyflux.eflx_ventilation_patch[p] = 0.0
                energyflux.eflx_heat_from_ac_patch[p] = 0.0
                energyflux.eflx_traffic_patch[p] = 0.0
            end

            eflx_gnet[p] = sabg + dlrad - eflx_lwrad_net -
                (eflx_sh_grnd + qflx_evap_soi * htvp[c] + qflx_tran_veg * HVAP) +
                energyflux.eflx_wasteheat_patch[p] + energyflux.eflx_heat_from_ac_patch[p] +
                energyflux.eflx_traffic_patch[p] + energyflux.eflx_ventilation_patch[p]

            if is_simple_build_temp()
                energyflux.eflx_anthro_patch[p] = energyflux.eflx_wasteheat_patch[p] +
                    energyflux.eflx_traffic_patch[p]
            end
            eflx_gnet_snow = eflx_gnet[p]
            eflx_gnet_soil = eflx_gnet[p]
            eflx_gnet_h2osfc = eflx_gnet[p]
        end

        dgnetdT[p] = -cgrnd - dlwrad_emit[c]
        hs[c] += eflx_gnet[p] * patch_data.wtcol[p]
        dhsdT[c] += dgnetdT[p] * patch_data.wtcol[p]
        hs_soil[c] += eflx_gnet_soil * patch_data.wtcol[p]
        hs_h2osfc[c] += eflx_gnet_h2osfc * patch_data.wtcol[p]
    end

    # SNICAR: sabg_lyr_col and hs_top
    sabg_lyr_col .= 0.0
    hs_top .= 0.0
    hs_top_snow .= 0.0

    for p in bounds_patch
        mask_nolakep[p] || continue
        c = patch_data.column[p]
        l = patch_data.landunit[p]
        lyr_top = snl[c] + 1

        frac_veg_nosno = canopystate.frac_veg_nosno_patch[p]
        dlrad = energyflux.dlrad_patch[p]
        eflx_sh_grnd = energyflux.eflx_sh_grnd_patch[p]
        qflx_evap_soi = waterfluxbulk.wf.qflx_evap_soi_patch[p]
        sabg_lyr = solarabs.sabg_lyr_patch

        if !lun.urbpoi[l]
            eflx_gnet_top = sabg_lyr[p, lyr_top + joff] + dlrad +
                (1.0 - frac_veg_nosno) * emg[c] * forc_lwrad[c] - lwrad_emit[c] -
                (eflx_sh_grnd + qflx_evap_soi * htvp[c])

            hs_top[c] += eflx_gnet_top * patch_data.wtcol[p]

            eflx_sh_snow = energyflux.eflx_sh_snow_patch[p]
            qflx_ev_snow = waterfluxbulk.qflx_ev_snow_patch[p]
            eflx_gnet_snow = sabg_lyr[p, lyr_top + joff] + dlrad +
                (1.0 - frac_veg_nosno) * emg[c] * forc_lwrad[c] - lwrad_emit_snow[c] -
                (eflx_sh_snow + qflx_ev_snow * htvp[c])

            hs_top_snow[c] += eflx_gnet_snow * patch_data.wtcol[p]

            for j in lyr_top:1
                jj_slyr = j + joff  # index into sabg_lyr_col
                # sabg_lyr_col is indexed [c, j_in_col_space] where j ranges -nlevsno+1:1
                # Map to matrix column: j - (-nlevsno+1) + 1 = j + nlevsno
                sabg_lyr_col[c, j + nlevsno] += sabg_lyr[p, j + joff] * patch_data.wtcol[p]
            end
        else
            hs_top[c] += eflx_gnet[p] * patch_data.wtcol[p]
            hs_top_snow[c] += eflx_gnet[p] * patch_data.wtcol[p]
            sabg_lyr_col[c, lyr_top + nlevsno] += solarabs.sabg_patch[p] * patch_data.wtcol[p]
        end
    end

    return nothing
end

# =========================================================================
# compute_heat_diff_flux_and_factor!
# =========================================================================
function compute_heat_diff_flux_and_factor!(
        col::ColumnData, lun::LandunitData,
        temperature::TemperatureData, energyflux::EnergyFluxData,
        urbanparams::UrbanParamsData,
        mask_nolakec::BitVector, bounds_col::UnitRange{Int},
        dtime::Float64,
        tk::Matrix{Float64}, cv::Matrix{Float64},
        fn::Matrix{Float64}, fact::Matrix{Float64})

    nlevsno = varpar.nlevsno
    nlevgrnd = varpar.nlevgrnd
    nlevurb = varpar.nlevurb
    nlevmaxurbgrnd = varpar.nlevmaxurbgrnd
    joff = nlevsno

    t_soisno = temperature.t_soisno_col
    t_building = temperature.t_building_lun
    t_roof_inner = temperature.t_roof_inner_lun
    t_sunw_inner = temperature.t_sunw_inner_lun
    t_shdw_inner = temperature.t_shdw_inner_lun
    eflx_bot = energyflux.eflx_bot_col

    for j in (-nlevsno + 1):nlevmaxurbgrnd
        for c in bounds_col
            mask_nolakec[c] || continue
            l = col.landunit[c]
            jj = j + joff

            if (col.itype[c] == ICOL_SUNWALL || col.itype[c] == ICOL_SHADEWALL ||
                col.itype[c] == ICOL_ROOF) && j <= nlevurb
                if j >= col.snl[c] + 1
                    if j == col.snl[c] + 1
                        fact[c, jj] = dtime / cv[c, jj]
                        fn[c, jj] = tk[c, jj] * (t_soisno[c, jj + 1] - t_soisno[c, jj]) /
                                    (col.z[c, jj + 1] - col.z[c, jj])
                    elseif j <= nlevurb - 1
                        fact[c, jj] = dtime / cv[c, jj]
                        fn[c, jj] = tk[c, jj] * (t_soisno[c, jj + 1] - t_soisno[c, jj]) /
                                    (col.z[c, jj + 1] - col.z[c, jj])
                    elseif j == nlevurb
                        fact[c, jj] = dtime / cv[c, jj]
                        if is_simple_build_temp()
                            fn[c, jj] = tk[c, jj] * (t_building[l] - CNFAC * t_soisno[c, jj]) /
                                        (col.zi[c, jj + 1] - col.z[c, jj])
                        else
                            if col.itype[c] == ICOL_SUNWALL
                                fn[c, jj] = tk[c, jj] * (t_sunw_inner[l] - CNFAC * t_soisno[c, jj]) /
                                            (col.zi[c, jj + 1] - col.z[c, jj])
                            elseif col.itype[c] == ICOL_SHADEWALL
                                fn[c, jj] = tk[c, jj] * (t_shdw_inner[l] - CNFAC * t_soisno[c, jj]) /
                                            (col.zi[c, jj + 1] - col.z[c, jj])
                            elseif col.itype[c] == ICOL_ROOF
                                fn[c, jj] = tk[c, jj] * (t_roof_inner[l] - CNFAC * t_soisno[c, jj]) /
                                            (col.zi[c, jj + 1] - col.z[c, jj])
                            end
                        end
                    end
                end
            elseif col.itype[c] != ICOL_SUNWALL && col.itype[c] != ICOL_SHADEWALL &&
                   col.itype[c] != ICOL_ROOF && j <= nlevgrnd
                if j >= col.snl[c] + 1
                    if j == col.snl[c] + 1
                        fact[c, jj] = dtime / cv[c, jj] * col.dz[c, jj] /
                            (0.5 * (col.z[c, jj] - col.zi[c, jj] + CAPR * (col.z[c, jj + 1] - col.zi[c, jj])))
                        fn[c, jj] = tk[c, jj] * (t_soisno[c, jj + 1] - t_soisno[c, jj]) /
                                    (col.z[c, jj + 1] - col.z[c, jj])
                    elseif j <= nlevgrnd - 1
                        fact[c, jj] = dtime / cv[c, jj]
                        fn[c, jj] = tk[c, jj] * (t_soisno[c, jj + 1] - t_soisno[c, jj]) /
                                    (col.z[c, jj + 1] - col.z[c, jj])
                    elseif j == nlevgrnd
                        fact[c, jj] = dtime / cv[c, jj]
                        fn[c, jj] = eflx_bot[c]
                    end
                end
            end
        end
    end

    return nothing
end

# =========================================================================
# set_rhs_vec! — Combines snow, SSW, and soil RHS vectors
# =========================================================================
function set_rhs_vec!(col::ColumnData, lun::LandunitData,
                      temperature::TemperatureData,
                      waterdiagbulk::WaterDiagnosticBulkData,
                      mask_nolakec::BitVector, bounds_col::UnitRange{Int},
                      dtime::Float64,
                      hs_h2osfc::Vector{Float64}, hs_top_snow::Vector{Float64},
                      hs_soil::Vector{Float64}, hs_top::Vector{Float64},
                      dhsdT::Vector{Float64}, sabg_lyr_col::Matrix{Float64},
                      tk::Matrix{Float64}, tk_h2osfc::Vector{Float64},
                      fact::Matrix{Float64}, fn::Matrix{Float64},
                      c_h2osfc::Vector{Float64}, dz_h2osfc::Vector{Float64},
                      rvector::Matrix{Float64})

    nlevsno = varpar.nlevsno
    nlevgrnd = varpar.nlevgrnd
    nlevurb = varpar.nlevurb
    nlevmaxurbgrnd = varpar.nlevmaxurbgrnd
    joff = nlevsno

    t_soisno = temperature.t_soisno_col
    t_h2osfc = temperature.t_h2osfc_col
    frac_h2osfc = waterdiagbulk.frac_h2osfc_col
    frac_sno_eff = waterdiagbulk.frac_sno_eff_col

    rvector .= NaN

    fn_h2osfc = zeros(length(bounds_col))

    # Snow layers RHS
    for j in (-nlevsno + 1):0
        for c in bounds_col
            mask_nolakec[c] || continue
            jj = j + joff
            # rv_idx: Fortran rvector index j-1 → Julia j-1+nlevsno+1
            rv_idx = j - 1 + nlevsno + 1

            hs_top_lev = hs_top_snow[c]
            if col.itype[c] == ICOL_SUNWALL || col.itype[c] == ICOL_SHADEWALL || col.itype[c] == ICOL_ROOF
                hs_top_lev = hs_top[c]
            end

            if j == col.snl[c] + 1
                rvector[c, rv_idx] = t_soisno[c, jj] + fact[c, jj] * (hs_top_lev -
                    dhsdT[c] * t_soisno[c, jj] + CNFAC * fn[c, jj])
            elseif j > col.snl[c] + 1
                rvector[c, rv_idx] = t_soisno[c, jj] + CNFAC * fact[c, jj] * (fn[c, jj] - fn[c, jj - 1])
                rvector[c, rv_idx] += fact[c, jj] * sabg_lyr_col[c, j + nlevsno]
            end
        end
    end

    # Standing surface water RHS
    for c in bounds_col
        mask_nolakec[c] || continue
        dzm = 0.5 * dz_h2osfc[c] + col.z[c, 1 + joff]
        fn_h2osfc[c] = tk_h2osfc[c] * (t_soisno[c, 1 + joff] - t_h2osfc[c]) / dzm
        # SSW RHS at index nlevsno+1 (Fortran index 0)
        rvector[c, nlevsno + 1] = t_h2osfc[c] + (dtime / c_h2osfc[c]) *
            (hs_h2osfc[c] - dhsdT[c] * t_h2osfc[c] + CNFAC * fn_h2osfc[c])
    end

    # Soil layers RHS — urban non-road
    for j in 1:nlevurb
        for c in bounds_col
            mask_nolakec[c] || continue
            if col.itype[c] == ICOL_SUNWALL || col.itype[c] == ICOL_SHADEWALL || col.itype[c] == ICOL_ROOF
                jj = j + joff
                rv_idx = j + nlevsno + 1
                if j == col.snl[c] + 1
                    rvector[c, rv_idx] = t_soisno[c, jj] + fact[c, jj] * (hs_top[c] -
                        dhsdT[c] * t_soisno[c, jj] + CNFAC * fn[c, jj])
                elseif j <= nlevurb - 1
                    rvector[c, rv_idx] = t_soisno[c, jj] + CNFAC * fact[c, jj] * (fn[c, jj] - fn[c, jj - 1])
                    if j == 1
                        rvector[c, rv_idx] += fact[c, jj] * sabg_lyr_col[c, j + nlevsno]
                    end
                elseif j == nlevurb
                    rvector[c, rv_idx] = t_soisno[c, jj] + fact[c, jj] * (fn[c, jj] - CNFAC * fn[c, jj - 1])
                end
            end
        end
    end

    # Soil layers RHS — non-urban and urban road
    for j in 1:nlevgrnd
        for c in bounds_col
            mask_nolakec[c] || continue
            l = col.landunit[c]
            if !lun.urbpoi[l] || col.itype[c] == ICOL_ROAD_IMPERV || col.itype[c] == ICOL_ROAD_PERV
                jj = j + joff
                rv_idx = j + nlevsno + 1
                if j == col.snl[c] + 1
                    rvector[c, rv_idx] = t_soisno[c, jj] + fact[c, jj] * (hs_top_snow[c] -
                        dhsdT[c] * t_soisno[c, jj] + CNFAC * fn[c, jj])
                elseif j == 1
                    rvector[c, rv_idx] = t_soisno[c, jj] + fact[c, jj] *
                        ((1.0 - frac_sno_eff[c]) * (hs_soil[c] - dhsdT[c] * t_soisno[c, jj]) +
                         CNFAC * (fn[c, jj] - frac_sno_eff[c] * fn[c, jj - 1]))
                    rvector[c, rv_idx] += frac_sno_eff[c] * fact[c, jj] * sabg_lyr_col[c, j + nlevsno]
                elseif j <= nlevgrnd - 1
                    rvector[c, rv_idx] = t_soisno[c, jj] + CNFAC * fact[c, jj] * (fn[c, jj] - fn[c, jj - 1])
                elseif j == nlevgrnd
                    rvector[c, rv_idx] = t_soisno[c, jj] - CNFAC * fact[c, jj] * fn[c, jj - 1] +
                        fact[c, jj] * fn[c, jj]
                end
            end
        end
    end

    # Surface water correction to soil layer 1
    for c in bounds_col
        mask_nolakec[c] || continue
        if frac_h2osfc[c] != 0.0
            rv_idx = 1 + nlevsno + 1
            jj = 1 + joff
            rvector[c, rv_idx] -= frac_h2osfc[c] * fact[c, jj] *
                ((hs_soil[c] - dhsdT[c] * t_soisno[c, jj]) + CNFAC * fn_h2osfc[c])
        end
    end

    return nothing
end

# =========================================================================
# set_matrix! — Assembles band diagonal matrix from submatrices
# =========================================================================
function set_matrix!(col::ColumnData, lun::LandunitData,
                     waterdiagbulk::WaterDiagnosticBulkData,
                     mask_nolakec::BitVector, bounds_col::UnitRange{Int},
                     dtime::Float64, nband::Int,
                     dhsdT::Vector{Float64}, tk::Matrix{Float64},
                     tk_h2osfc::Vector{Float64}, fact::Matrix{Float64},
                     c_h2osfc::Vector{Float64}, dz_h2osfc::Vector{Float64},
                     bmatrix::Array{Float64,3})

    nlevsno = varpar.nlevsno
    nlevgrnd = varpar.nlevgrnd
    nlevurb = varpar.nlevurb
    nlevmaxurbgrnd = varpar.nlevmaxurbgrnd
    joff = nlevsno

    frac_h2osfc = waterdiagbulk.frac_h2osfc_col
    frac_sno_eff = waterdiagbulk.frac_sno_eff_col

    bmatrix .= 0.0

    # ---- Snow submatrix ----
    for j in (-nlevsno + 1):0
        for c in bounds_col
            mask_nolakec[c] || continue
            jj = j + joff
            # bmatrix index: j-1+nlevsno+1 = j+nlevsno
            bm_idx = j - 1 + nlevsno + 1  # maps to tvector index space

            # Determine nlev_thresh
            nlev_thresh = nlevgrnd
            if col.itype[c] == ICOL_SUNWALL || col.itype[c] == ICOL_SHADEWALL || col.itype[c] == ICOL_ROOF
                nlev_thresh = nlevurb
            end

            if j >= col.snl[c] + 1
                dzp = col.z[c, jj + 1] - col.z[c, jj]
                if j == col.snl[c] + 1
                    bmatrix[c, 3, bm_idx] = 1.0 + (1.0 - CNFAC) * fact[c, jj] * tk[c, jj] / dzp -
                        fact[c, jj] * dhsdT[c]
                else
                    dzm = col.z[c, jj] - col.z[c, jj - 1]
                    bmatrix[c, 4, bm_idx] = -(1.0 - CNFAC) * fact[c, jj] * tk[c, jj - 1] / dzm
                    bmatrix[c, 3, bm_idx] = 1.0 + (1.0 - CNFAC) * fact[c, jj] *
                        (tk[c, jj] / dzp + tk[c, jj - 1] / dzm)
                end
                if j != 0
                    bmatrix[c, 2, bm_idx] = -(1.0 - CNFAC) * fact[c, jj] * tk[c, jj] / dzp
                else
                    # snow-soil coupling: goes to soil layer 1 (2 rows below in band)
                    # bmatrix_snow_soil → band 1
                    bmatrix[c, 1, bm_idx] = -(1.0 - CNFAC) * fact[c, jj] * tk[c, jj] / dzp
                end
            end
        end
    end

    # ---- Standing surface water submatrix ----
    ssw_idx = nlevsno + 1  # index for SSW row in bmatrix
    for c in bounds_col
        mask_nolakec[c] || continue
        dzm = 0.5 * dz_h2osfc[c] + col.z[c, 1 + joff]
        bmatrix[c, 3, ssw_idx] = 1.0 + (1.0 - CNFAC) * (dtime / c_h2osfc[c]) *
            tk_h2osfc[c] / dzm - (dtime / c_h2osfc[c]) * dhsdT[c]
        # SSW-Soil coupling (band 2)
        bmatrix[c, 2, ssw_idx] = -(1.0 - CNFAC) * (dtime / c_h2osfc[c]) * tk_h2osfc[c] / dzm
    end

    # ---- Soil submatrix — urban non-road ----
    for j in 1:nlevurb
        for c in bounds_col
            mask_nolakec[c] || continue
            if col.itype[c] == ICOL_SUNWALL || col.itype[c] == ICOL_SHADEWALL || col.itype[c] == ICOL_ROOF
                jj = j + joff
                bm_idx = j + nlevsno + 1

                if j == col.snl[c] + 1
                    dzp = col.z[c, jj + 1] - col.z[c, jj]
                    bmatrix[c, 3, bm_idx] = 1.0 + (1.0 - CNFAC) * fact[c, jj] * tk[c, jj] / dzp -
                        fact[c, jj] * dhsdT[c]
                    bmatrix[c, 2, bm_idx] = -(1.0 - CNFAC) * fact[c, jj] * tk[c, jj] / dzp
                elseif j <= nlevurb - 1
                    dzm = col.z[c, jj] - col.z[c, jj - 1]
                    dzp = col.z[c, jj + 1] - col.z[c, jj]
                    if j != 1
                        bmatrix[c, 4, bm_idx] = -(1.0 - CNFAC) * fact[c, jj] * tk[c, jj - 1] / dzm
                    else
                        bmatrix[c, 5, bm_idx] = -(1.0 - CNFAC) * fact[c, jj] * tk[c, jj - 1] / dzm
                    end
                    bmatrix[c, 3, bm_idx] = 1.0 + (1.0 - CNFAC) * fact[c, jj] *
                        (tk[c, jj] / dzp + tk[c, jj - 1] / dzm)
                    bmatrix[c, 2, bm_idx] = -(1.0 - CNFAC) * fact[c, jj] * tk[c, jj] / dzp
                elseif j == nlevurb
                    dzm = col.z[c, jj] - col.z[c, jj - 1]
                    dzp = col.zi[c, jj + 1] - col.z[c, jj]
                    bmatrix[c, 4, bm_idx] = -(1.0 - CNFAC) * fact[c, jj] * tk[c, jj - 1] / dzm
                    bmatrix[c, 3, bm_idx] = 1.0 + (1.0 - CNFAC) * fact[c, jj] *
                        (tk[c, jj - 1] / dzm + tk[c, jj] / dzp)
                end
            end
        end
    end

    # ---- Soil submatrix — non-urban and urban road ----
    for j in 1:nlevgrnd
        for c in bounds_col
            mask_nolakec[c] || continue
            l = col.landunit[c]
            if col.itype[c] == ICOL_ROAD_IMPERV || col.itype[c] == ICOL_ROAD_PERV || !lun.urbpoi[l]
                jj = j + joff
                bm_idx = j + nlevsno + 1

                if j == col.snl[c] + 1
                    dzp = col.z[c, jj + 1] - col.z[c, jj]
                    bmatrix[c, 3, bm_idx] = 1.0 + (1.0 - CNFAC) * fact[c, jj] * tk[c, jj] / dzp -
                        fact[c, jj] * dhsdT[c]
                    bmatrix[c, 2, bm_idx] = -(1.0 - CNFAC) * fact[c, jj] * tk[c, jj] / dzp
                elseif j == 1
                    dzm = col.z[c, jj] - col.z[c, jj - 1]
                    dzp = col.z[c, jj + 1] - col.z[c, jj]
                    bmatrix[c, 2, bm_idx] = -(1.0 - CNFAC) * fact[c, jj] * tk[c, jj] / dzp
                    bmatrix[c, 3, bm_idx] = 1.0 + (1.0 - CNFAC) * fact[c, jj] *
                        (tk[c, jj] / dzp + frac_sno_eff[c] * tk[c, jj - 1] / dzm) -
                        (1.0 - frac_sno_eff[c]) * fact[c, jj] * dhsdT[c]
                    bmatrix[c, 5, bm_idx] = -frac_sno_eff[c] * (1.0 - CNFAC) * fact[c, jj] *
                        tk[c, jj - 1] / dzm
                elseif j <= nlevgrnd - 1
                    dzm = col.z[c, jj] - col.z[c, jj - 1]
                    dzp = col.z[c, jj + 1] - col.z[c, jj]
                    bmatrix[c, 2, bm_idx] = -(1.0 - CNFAC) * fact[c, jj] * tk[c, jj] / dzp
                    bmatrix[c, 3, bm_idx] = 1.0 + (1.0 - CNFAC) * fact[c, jj] *
                        (tk[c, jj] / dzp + tk[c, jj - 1] / dzm)
                    bmatrix[c, 4, bm_idx] = -(1.0 - CNFAC) * fact[c, jj] * tk[c, jj - 1] / dzm
                elseif j == nlevgrnd
                    dzm = col.z[c, jj] - col.z[c, jj - 1]
                    bmatrix[c, 3, bm_idx] = 1.0 + (1.0 - CNFAC) * fact[c, jj] * tk[c, jj - 1] / dzm
                    bmatrix[c, 4, bm_idx] = -(1.0 - CNFAC) * fact[c, jj] * tk[c, jj - 1] / dzm
                end
            end
        end
    end

    # h2osfc correction to soil layer 1 diagonal
    for c in bounds_col
        mask_nolakec[c] || continue
        if frac_h2osfc[c] != 0.0
            dzm = 0.5 * dz_h2osfc[c] + col.z[c, 1 + joff]
            bm_idx = 1 + nlevsno + 1
            bmatrix[c, 3, bm_idx] += frac_h2osfc[c] *
                ((1.0 - CNFAC) * fact[c, 1 + joff] * tk_h2osfc[c] / dzm +
                 fact[c, 1 + joff] * dhsdT[c])
            # Soil-SSW coupling (band 4)
            bmatrix[c, 4, bm_idx] = -frac_h2osfc[c] * (1.0 - CNFAC) * fact[c, 1 + joff] *
                tk_h2osfc[c] / dzm
        end
    end

    return nothing
end

# =========================================================================
# phase_change_h2osfc! — Phase change of surface water
# =========================================================================
function phase_change_h2osfc!(col::ColumnData, temperature::TemperatureData,
                              energyflux::EnergyFluxData,
                              waterstatebulk::WaterStateBulkData,
                              waterdiagbulk::WaterDiagnosticBulkData,
                              waterfluxbulk::WaterFluxBulkData,
                              mask_nolakec::BitVector, bounds_col::UnitRange{Int},
                              dtime::Float64, dhsdT::Vector{Float64})

    nlevsno = varpar.nlevsno
    joff = nlevsno

    snl = col.snl
    t_soisno = temperature.t_soisno_col
    t_h2osfc = temperature.t_h2osfc_col
    fact = temperature.fact_col
    c_h2osfc = temperature.c_h2osfc_col
    xmf_h2osfc = temperature.xmf_h2osfc_col

    frac_sno = waterdiagbulk.frac_sno_eff_col
    frac_h2osfc = waterdiagbulk.frac_h2osfc_col
    snow_depth = waterdiagbulk.snow_depth_col
    h2osfc = waterstatebulk.ws.h2osfc_col
    h2osno_no_layers = waterstatebulk.ws.h2osno_no_layers_col
    h2osoi_ice = waterstatebulk.ws.h2osoi_ice_col
    int_snow = waterstatebulk.int_snow_col
    qflx_h2osfc_to_ice = waterfluxbulk.wf.qflx_h2osfc_to_ice_col
    eflx_h2osfc_to_snow = energyflux.eflx_h2osfc_to_snow_col

    # Initialize
    for c in bounds_col
        mask_nolakec[c] || continue
        xmf_h2osfc[c] = 0.0
        qflx_h2osfc_to_ice[c] = 0.0
        eflx_h2osfc_to_snow[c] = 0.0
    end

    # Freezing identification
    for c in bounds_col
        mask_nolakec[c] || continue
        if frac_h2osfc[c] > 0.0 && t_h2osfc[c] <= TFRZ
            tinc = TFRZ - t_h2osfc[c]
            t_h2osfc[c] = TFRZ

            hm = frac_h2osfc[c] * (dhsdT[c] * tinc - tinc * c_h2osfc[c] / dtime)
            xm = hm * dtime / HFUS
            temp1 = h2osfc[c] + xm

            # Compute h2osno_total
            h2osno_total = h2osno_no_layers[c]
            if snl[c] < 0
                for j in (snl[c] + 1):0
                    h2osno_total += h2osoi_ice[c, j + joff] + waterstatebulk.ws.h2osoi_liq_col[c, j + joff]
                end
            end

            z_avg = frac_sno[c] * snow_depth[c]
            rho_avg = z_avg > 0.0 ? min(800.0, h2osno_total / z_avg) : 200.0

            if temp1 >= 0.0
                int_snow[c] -= xm
                if snl[c] == 0
                    h2osno_no_layers[c] -= xm
                else
                    h2osoi_ice[c, joff] -= xm  # layer 0
                end
                h2osno_total -= xm
                h2osfc[c] += xm
                xmf_h2osfc[c] = hm
                qflx_h2osfc_to_ice[c] = -xm / dtime

                if frac_sno[c] > 0 && snl[c] < 0
                    snow_depth[c] = h2osno_total / (rho_avg * frac_sno[c])
                else
                    snow_depth[c] = h2osno_total / DENICE
                end

                if snl[c] == 0
                    t_soisno[c, joff] = t_h2osfc[c]
                    eflx_h2osfc_to_snow[c] = 0.0
                else
                    if snl[c] == -1
                        c1 = frac_sno[c] * (dtime / fact[c, joff] - dhsdT[c] * dtime)
                    else
                        c1 = frac_sno[c] / fact[c, joff] * dtime
                    end
                    c2 = frac_h2osfc[c] != 0.0 ? (-CPLIQ * xm - frac_h2osfc[c] * dhsdT[c] * dtime) : 0.0
                    t_soisno[c, joff] = (c1 * t_soisno[c, joff] + c2 * t_h2osfc[c]) / (c1 + c2)
                    eflx_h2osfc_to_snow[c] = (t_h2osfc[c] - t_soisno[c, joff]) * c2 / dtime
                end
            else
                rho_avg = (h2osno_total * rho_avg + h2osfc[c] * DENICE) / (h2osno_total + h2osfc[c])
                int_snow[c] += h2osfc[c]
                if snl[c] == 0
                    h2osno_no_layers[c] += h2osfc[c]
                else
                    h2osoi_ice[c, joff] += h2osfc[c]
                end
                h2osno_total += h2osfc[c]
                qflx_h2osfc_to_ice[c] = h2osfc[c] / dtime
                t_h2osfc[c] -= temp1 * HFUS / (dtime * dhsdT[c] - c_h2osfc[c])
                xmf_h2osfc[c] = hm - frac_h2osfc[c] * temp1 * HFUS / dtime

                if snl[c] == 0
                    t_soisno[c, joff] = t_h2osfc[c]
                elseif snl[c] == -1
                    c1 = frac_sno[c] * (dtime / fact[c, joff] - dhsdT[c] * dtime)
                    c2 = frac_h2osfc[c] != 0.0 ? frac_h2osfc[c] * (c_h2osfc[c] - dtime * dhsdT[c]) : 0.0
                    t_soisno[c, joff] = (c1 * t_soisno[c, joff] + c2 * t_h2osfc[c]) / (c1 + c2)
                    t_h2osfc[c] = t_soisno[c, joff]
                else
                    c1 = frac_sno[c] / fact[c, joff] * dtime
                    c2 = frac_h2osfc[c] != 0.0 ? frac_h2osfc[c] * (c_h2osfc[c] - dtime * dhsdT[c]) : 0.0
                    t_soisno[c, joff] = (c1 * t_soisno[c, joff] + c2 * t_h2osfc[c]) / (c1 + c2)
                    t_h2osfc[c] = t_soisno[c, joff]
                end

                h2osfc[c] = 0.0
                if frac_sno[c] > 0 && snl[c] < 0
                    snow_depth[c] = h2osno_total / (rho_avg * frac_sno[c])
                else
                    snow_depth[c] = h2osno_total / DENICE
                end
            end
        end
    end

    return nothing
end

# =========================================================================
# phase_change_beta! — Phase change within snow and soil layers
# =========================================================================
function phase_change_beta!(col::ColumnData, lun::LandunitData,
                            temperature::TemperatureData, energyflux::EnergyFluxData,
                            soilstate::SoilStateData,
                            waterstatebulk::WaterStateBulkData,
                            waterdiagbulk::WaterDiagnosticBulkData,
                            waterfluxbulk::WaterFluxBulkData,
                            mask_nolakec::BitVector, bounds_col::UnitRange{Int},
                            dtime::Float64, dhsdT::Vector{Float64})

    nlevsno = varpar.nlevsno
    nlevgrnd = varpar.nlevgrnd
    nlevurb = varpar.nlevurb
    nlevmaxurbgrnd = varpar.nlevmaxurbgrnd
    joff = nlevsno

    snl = col.snl
    t_soisno = temperature.t_soisno_col
    xmf = temperature.xmf_col
    fact = temperature.fact_col
    imelt = temperature.imelt_col

    frac_sno_eff = waterdiagbulk.frac_sno_eff_col
    frac_h2osfc = waterdiagbulk.frac_h2osfc_col
    snow_depth = waterdiagbulk.snow_depth_col
    h2osno_no_layers = waterstatebulk.ws.h2osno_no_layers_col
    h2osoi_liq = waterstatebulk.ws.h2osoi_liq_col
    h2osoi_ice = waterstatebulk.ws.h2osoi_ice_col
    excess_ice = waterstatebulk.ws.excess_ice_col

    bsw = soilstate.bsw_col
    sucsat = soilstate.sucsat_col
    watsat = soilstate.watsat_col

    qflx_snomelt = waterfluxbulk.wf.qflx_snomelt_col
    qflx_snofrz = waterfluxbulk.wf.qflx_snofrz_col
    qflx_snow_drain = waterfluxbulk.wf.qflx_snow_drain_col
    qflx_snofrz_lyr = waterfluxbulk.wf.qflx_snofrz_lyr_col
    qflx_snomelt_lyr = waterfluxbulk.wf.qflx_snomelt_lyr_col

    eflx_snomelt = energyflux.eflx_snomelt_col
    eflx_snomelt_r = energyflux.eflx_snomelt_r_col
    eflx_snomelt_u = energyflux.eflx_snomelt_u_col
    exice_subs = waterdiagbulk.exice_subs_col

    nc = length(bounds_col)

    # Initialization
    for c in bounds_col
        mask_nolakec[c] || continue
        xmf[c] = 0.0
        qflx_snomelt[c] = 0.0
        qflx_snofrz[c] = 0.0
        qflx_snow_drain[c] = 0.0
    end

    # Local arrays
    hm = zeros(nc, nlevsno + nlevmaxurbgrnd)
    xm = zeros(nc, nlevsno + nlevmaxurbgrnd)
    xm2 = zeros(nc, nlevsno + nlevmaxurbgrnd)
    wice0 = zeros(nc, nlevsno + nlevmaxurbgrnd)
    wliq0 = zeros(nc, nlevsno + nlevmaxurbgrnd)
    wexice0 = zeros(nc, nlevsno + nlevmaxurbgrnd)
    wmass0 = zeros(nc, nlevsno + nlevmaxurbgrnd)
    supercool = zeros(nc, nlevmaxurbgrnd)
    tinc = zeros(nc, nlevsno + nlevmaxurbgrnd)

    # Initialize layer variables
    for j in (-nlevsno + 1):nlevmaxurbgrnd
        for c in bounds_col
            mask_nolakec[c] || continue
            jj = j + joff
            if j >= snl[c] + 1
                imelt[c, jj] = 0
                hm[c, jj] = 0.0
                xm[c, jj] = 0.0
                xm2[c, jj] = 0.0
                wice0[c, jj] = h2osoi_ice[c, jj]
                wliq0[c, jj] = h2osoi_liq[c, jj]
                wexice0[c, jj] = j >= 1 ? excess_ice[c, j] : 0.0
                wmass0[c, jj] = h2osoi_ice[c, jj] + h2osoi_liq[c, jj] + wexice0[c, jj]
                if j >= 1
                    exice_subs[c, j] = 0.0
                end
            end
            if j <= 0
                qflx_snomelt_lyr[c, jj] = 0.0
                qflx_snofrz_lyr[c, jj] = 0.0
            end
        end
    end

    # Snow layers: melting/freezing identification
    for j in (-nlevsno + 1):0
        for c in bounds_col
            mask_nolakec[c] || continue
            jj = j + joff
            if j >= snl[c] + 1
                if h2osoi_ice[c, jj] > 0.0 && t_soisno[c, jj] > TFRZ
                    imelt[c, jj] = 1
                    tinc[c, jj] = TFRZ - t_soisno[c, jj]
                    t_soisno[c, jj] = TFRZ
                end
                if h2osoi_liq[c, jj] > 0.0 && t_soisno[c, jj] < TFRZ
                    imelt[c, jj] = 2
                    tinc[c, jj] = TFRZ - t_soisno[c, jj]
                    t_soisno[c, jj] = TFRZ
                end
            end
        end
    end

    # Soil layers: melting/freezing identification
    for j in 1:nlevmaxurbgrnd
        for c in bounds_col
            mask_nolakec[c] || continue
            l = col.landunit[c]
            jj = j + joff
            supercool[c, j] = 0.0

            active = (col.itype[c] != ICOL_SUNWALL && col.itype[c] != ICOL_SHADEWALL &&
                      col.itype[c] != ICOL_ROOF && j <= nlevgrnd) ||
                     ((col.itype[c] == ICOL_SUNWALL || col.itype[c] == ICOL_SHADEWALL ||
                       col.itype[c] == ICOL_ROOF) && j <= nlevurb)

            if active
                if h2osoi_ice[c, jj] > 0.0 && t_soisno[c, jj] > TFRZ
                    imelt[c, jj] = 1
                    tinc[c, jj] = TFRZ - t_soisno[c, jj]
                    t_soisno[c, jj] = TFRZ
                end

                if excess_ice[c, j] > 0.0 && t_soisno[c, jj] > TFRZ
                    imelt[c, jj] = 1
                    tinc[c, jj] = TFRZ - t_soisno[c, jj]
                    t_soisno[c, jj] = TFRZ
                end

                supercool[c, j] = 0.0
                if lun.itype[l] == ISTSOIL || lun.itype[l] == ISTCROP || col.itype[c] == ICOL_ROAD_PERV
                    if t_soisno[c, jj] < TFRZ
                        smp = HFUS * (TFRZ - t_soisno[c, jj]) / (GRAV * t_soisno[c, jj]) * 1000.0
                        supercool[c, j] = watsat[c, j] * (smp / sucsat[c, j])^(-1.0 / bsw[c, j])
                        supercool[c, j] *= col.dz[c, jj] * 1000.0
                    end
                end

                if h2osoi_liq[c, jj] > supercool[c, j] && t_soisno[c, jj] < TFRZ
                    imelt[c, jj] = 2
                    tinc[c, jj] = TFRZ - t_soisno[c, jj]
                    t_soisno[c, jj] = TFRZ
                end

                if h2osno_no_layers[c] > 0.0 && j == 1
                    if t_soisno[c, jj] > TFRZ
                        imelt[c, jj] = 1
                        tinc[c, jj] = TFRZ - t_soisno[c, jj]
                        t_soisno[c, jj] = TFRZ
                    end
                end
            end
        end
    end

    # Energy surplus/loss and rate of phase change
    for j in (-nlevsno + 1):nlevmaxurbgrnd
        for c in bounds_col
            mask_nolakec[c] || continue
            jj = j + joff

            active = (col.itype[c] != ICOL_SUNWALL && col.itype[c] != ICOL_SHADEWALL &&
                      col.itype[c] != ICOL_ROOF && j <= nlevgrnd) ||
                     ((col.itype[c] == ICOL_SUNWALL || col.itype[c] == ICOL_SHADEWALL ||
                       col.itype[c] == ICOL_ROOF) && j <= nlevurb)

            if active && j >= snl[c] + 1
                if imelt[c, jj] > 0
                    if j == snl[c] + 1
                        if j > 0
                            hm[c, jj] = dhsdT[c] * tinc[c, jj] - tinc[c, jj] / fact[c, jj]
                        else
                            hm[c, jj] = frac_sno_eff[c] * (dhsdT[c] * tinc[c, jj] - tinc[c, jj] / fact[c, jj])
                        end
                        if j == 1 && frac_h2osfc[c] != 0.0
                            hm[c, jj] -= frac_h2osfc[c] * dhsdT[c] * tinc[c, jj]
                        end
                    elseif j == 1
                        hm[c, jj] = (1.0 - frac_sno_eff[c] - frac_h2osfc[c]) *
                            dhsdT[c] * tinc[c, jj] - tinc[c, jj] / fact[c, jj]
                    else
                        if j < 1
                            hm[c, jj] = -frac_sno_eff[c] * tinc[c, jj] / fact[c, jj]
                        else
                            hm[c, jj] = -tinc[c, jj] / fact[c, jj]
                        end
                    end
                end

                if imelt[c, jj] == 1 && hm[c, jj] < 0.0
                    hm[c, jj] = 0.0
                    imelt[c, jj] = 0
                end
                if imelt[c, jj] == 2 && hm[c, jj] > 0.0
                    hm[c, jj] = 0.0
                    imelt[c, jj] = 0
                end

                if imelt[c, jj] > 0 && abs(hm[c, jj]) > 0.0
                    xm[c, jj] = hm[c, jj] * dtime / HFUS

                    # Unresolved snow melting
                    if j == 1
                        if h2osno_no_layers[c] > 0.0 && xm[c, jj] > 0.0
                            temp1 = h2osno_no_layers[c]
                            h2osno_no_layers[c] = max(0.0, temp1 - xm[c, jj])
                            propor = h2osno_no_layers[c] / temp1
                            snow_depth[c] = propor * snow_depth[c]
                            heatr = hm[c, jj] - HFUS * (temp1 - h2osno_no_layers[c]) / dtime
                            if heatr > 0.0
                                xm[c, jj] = heatr * dtime / HFUS
                                hm[c, jj] = heatr
                            else
                                xm[c, jj] = 0.0
                                hm[c, jj] = 0.0
                            end
                            qflx_snomelt[c] = max(0.0, temp1 - h2osno_no_layers[c]) / dtime
                            xmf[c] = HFUS * qflx_snomelt[c]
                            qflx_snow_drain[c] = qflx_snomelt[c]
                        end
                    end

                    heatr = 0.0
                    if xm[c, jj] > 0.0
                        h2osoi_ice[c, jj] = max(0.0, wice0[c, jj] - xm[c, jj])
                        heatr = hm[c, jj] - HFUS * (wice0[c, jj] - h2osoi_ice[c, jj]) / dtime
                        xm2[c, jj] = xm[c, jj] - h2osoi_ice[c, jj]
                        if h2osoi_ice[c, jj] == 0.0
                            if wexice0[c, jj] >= 0.0 && xm2[c, jj] > 0.0 && j >= 2
                                excess_ice[c, j] = max(0.0, wexice0[c, jj] - xm2[c, jj])
                                heatr = hm[c, jj] - HFUS * (wexice0[c, jj] - excess_ice[c, j] +
                                        wice0[c, jj] - h2osoi_ice[c, jj]) / dtime
                            end
                        end
                    elseif xm[c, jj] < 0.0
                        if j <= 0
                            h2osoi_ice[c, jj] = min(wmass0[c, jj], wice0[c, jj] - xm[c, jj])
                        else
                            if wmass0[c, jj] - wexice0[c, jj] < supercool[c, j]
                                h2osoi_ice[c, jj] = 0.0
                            else
                                h2osoi_ice[c, jj] = min(wmass0[c, jj] - wexice0[c, jj] - supercool[c, j],
                                    wice0[c, jj] - xm[c, jj])
                            end
                        end
                        heatr = hm[c, jj] - HFUS * (wice0[c, jj] - h2osoi_ice[c, jj]) / dtime
                    end

                    ei_val = j >= 1 ? excess_ice[c, j] : 0.0
                    h2osoi_liq[c, jj] = max(0.0, wmass0[c, jj] - h2osoi_ice[c, jj] - ei_val)

                    if abs(heatr) > 0.0
                        if j == snl[c] + 1
                            if j == 1
                                t_soisno[c, jj] += fact[c, jj] * heatr /
                                    (1.0 - (1.0 - frac_h2osfc[c]) * fact[c, jj] * dhsdT[c])
                            else
                                if frac_sno_eff[c] > 0.0
                                    t_soisno[c, jj] += (fact[c, jj] / frac_sno_eff[c]) * heatr /
                                        (1.0 - fact[c, jj] * dhsdT[c])
                                end
                            end
                        elseif j == 1
                            t_soisno[c, jj] += fact[c, jj] * heatr /
                                (1.0 - (1.0 - frac_sno_eff[c] - frac_h2osfc[c]) * fact[c, jj] * dhsdT[c])
                        else
                            if j > 0
                                t_soisno[c, jj] += fact[c, jj] * heatr
                            else
                                if frac_sno_eff[c] > 0.0
                                    t_soisno[c, jj] += (fact[c, jj] / frac_sno_eff[c]) * heatr
                                end
                            end
                        end

                        if j <= 0
                            if h2osoi_liq[c, jj] * h2osoi_ice[c, jj] > 0.0
                                t_soisno[c, jj] = TFRZ
                            end
                        end
                    end

                    if j >= 1
                        xmf[c] += HFUS * (wice0[c, jj] - h2osoi_ice[c, jj]) / dtime +
                            HFUS * (wexice0[c, jj] - (j >= 1 ? excess_ice[c, j] : 0.0)) / dtime
                        exice_subs[c, j] = max(0.0, (wexice0[c, jj] - excess_ice[c, j]) / DENICE)
                    else
                        xmf[c] += HFUS * (wice0[c, jj] - h2osoi_ice[c, jj]) / dtime
                    end

                    if imelt[c, jj] == 1 && j < 1
                        qflx_snomelt_lyr[c, jj] = max(0.0, wice0[c, jj] - h2osoi_ice[c, jj]) / dtime
                        qflx_snomelt[c] += qflx_snomelt_lyr[c, jj]
                        waterdiagbulk.snomelt_accum_col[c] += qflx_snomelt_lyr[c, jj] * dtime * 1.0e-3
                    end
                    if imelt[c, jj] == 2 && j < 1
                        qflx_snofrz_lyr[c, jj] = max(0.0, h2osoi_ice[c, jj] - wice0[c, jj]) / dtime
                        qflx_snofrz[c] += qflx_snofrz_lyr[c, jj]
                    end
                end
            end
        end
    end

    # History output
    for c in bounds_col
        mask_nolakec[c] || continue
        eflx_snomelt[c] = qflx_snomelt[c] * HFUS
        l = col.landunit[c]
        if lun.urbpoi[l]
            eflx_snomelt_u[c] = eflx_snomelt[c]
        elseif lun.itype[l] == ISTSOIL || lun.itype[l] == ISTCROP
            eflx_snomelt_r[c] = eflx_snomelt[c]
        end
    end

    return nothing
end

# =========================================================================
# building_hac! — Building Heating and Cooling (simple method)
# =========================================================================
function building_hac!(lun::LandunitData, temperature::TemperatureData,
                       urbanparams::UrbanParamsData,
                       t_building_max::Vector{Float64},
                       mask_urbanl::BitVector, bounds_lun::UnitRange{Int},
                       cool_on::BitVector, heat_on::BitVector)

    t_building = temperature.t_building_lun
    t_building_min = urbanparams.t_building_min

    for l in bounds_lun
        mask_urbanl[l] || continue
        if lun.urbpoi[l]
            cool_on[l] = false
            heat_on[l] = false
            if t_building[l] > t_building_max[l]
                t_building[l] = t_building_max[l]
                cool_on[l] = true
                heat_on[l] = false
            elseif t_building[l] < t_building_min[l]
                t_building[l] = t_building_min[l]
                cool_on[l] = false
                heat_on[l] = true
            end
        end
    end

    return nothing
end
