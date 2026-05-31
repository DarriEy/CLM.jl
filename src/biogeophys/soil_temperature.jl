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
                           urbantv_t_building_max::Vector{<:Real},
                           atm2lnd_forc_lwrad::Vector{<:Real},
                           mask_nolakec::BitVector, mask_nolakep::BitVector,
                           mask_urbanl::BitVector, mask_urbanc::BitVector,
                           bounds_col::UnitRange{Int}, bounds_lun::UnitRange{Int},
                           bounds_patch::UnitRange{Int},
                           dtime::Real)

    nlevsno = varpar.nlevsno
    nlevgrnd = varpar.nlevgrnd
    nlevurb = varpar.nlevurb
    nlevmaxurbgrnd = varpar.nlevmaxurbgrnd
    nlevsoi = varpar.nlevsoi
    joff = nlevsno  # offset: Fortran j → Julia j+joff

    nc = length(bounds_col)
    nband = 5

    # Local arrays — infer FT from temperature to support Dual
    FT = eltype(temperature.t_soisno_col)
    jtop = fill(-9999, nc)
    jbot = fill(0, nc)
    cv = zeros(FT, nc, nlevsno + nlevmaxurbgrnd)
    tk = zeros(FT, nc, nlevsno + nlevmaxurbgrnd)
    fn = zeros(FT, nc, nlevsno + nlevmaxurbgrnd)
    fn1 = zeros(FT, nc, nlevsno + nlevmaxurbgrnd)
    tk_h2osfc = fill(FT(NaN), nc)
    hs_top = zeros(FT, nc)
    dhsdT = zeros(FT, nc)
    hs_soil = zeros(FT, nc)
    hs_top_snow = zeros(FT, nc)
    hs_h2osfc = zeros(FT, nc)
    sabg_lyr_col = zeros(FT, nc, nlevsno + 1)  # -nlevsno+1:1
    fn_h2osfc = zeros(FT, nc)
    dz_h2osfc = zeros(FT, nc)
    c_h2osfc = zeros(FT, nc)
    cool_on = falses(length(bounds_lun))
    heat_on = falses(length(bounds_lun))

    # Band matrix and vectors (Fortran: -nlevsno:nlevmaxurbgrnd → nlevsno+1+nlevmaxurbgrnd levels)
    nlev_total = nlevsno + 1 + nlevmaxurbgrnd  # -nlevsno to nlevmaxurbgrnd
    bmatrix = zeros(FT, nc, nband, nlev_total)
    tvector = fill(FT(NaN), nc, nlev_total)
    rvector = fill(FT(NaN), nc, nlev_total)

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

    # Excess ice vertical coordinate adjustment (save originals).
    # Use typed empty arrays rather than `nothing`: a Union{Nothing,Matrix}
    # local breaks Enzyme reverse-mode codegen (SSA "does not dominate all uses").
    # These are only indexed when use_excess_ice is true (both save and restore
    # blocks are guarded), where they hold concrete copies.
    dz_0 = similar(col.dz, 0, 0)
    zi_0 = similar(col.zi, 0, 0)
    z_0 = similar(col.z, 0, 0)
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
            c_h2osfc_out[c] = smooth_max(THIN_SFCLAYER, CPLIQ * h2osfc[c] / frac_h2osfc[c])
            dz_h2osfc[c] = smooth_max(THIN_SFCLAYER, 1.0e-3 * h2osfc[c] / frac_h2osfc[c])
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

    # Batched pentadiagonal solve: one thread per column. Builds each column's dense
    # system directly from bmatrix (A[i,j] = bmatrix[c,(i-j)+kl+1,jt+i-1]) and solves
    # by partial-pivoted Gaussian elimination on per-column scratch — backend-agnostic
    # (CPU loop on Arrays, GPU kernel on device arrays), matching the prior per-column
    # LAPACK/pure-Julia band_solve! to round-off. ForwardDiff flows through the GE; the
    # host band_solve! + its Enzyme reverse rule remain for the non-kernel/Enzyme path.
    # (The prior NaN guard here was a workaround for an overflow bug in the smooth-AD
    # sigmoid, since fixed in _stable_sigmoid.)
    batched_band_solve!(tvector, bmatrix, rvector, jtop, jbot, mask_nolakec,
                        kl, ku, nlevsno)

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
# Harmonic-mean thermal conductivity at the interface below layer jj (Patankar form).
@inline _tk_iface(tkj, tkj1, zj, zj1, zij1) =
    tkj * tkj1 * (zj1 - zj) / (tkj * (zj1 - zij1) + tkj1 * (zij1 - zj))

# Per-(column, level) interface thermal conductivity tk_out from the per-layer thk.
# Backend-agnostic 2D kernel (one thread per (c, jj)); jj = j + nlevsno so the global
# level index maps directly. Replaces the host `for j … for c` loop in soil_therm_prop!.
@kernel function _soil_tk_interface_kernel!(tk_out, @Const(mask), @Const(itype),
        @Const(snl), @Const(z), @Const(zi), @Const(thk),
        nlevsno::Int, nlevurb::Int, nlevgrnd::Int)
    c, jj = @index(Global, NTuple)
    @inbounds if mask[c]
        j = jj - nlevsno
        it = itype[c]
        wallroof = (it == ICOL_SUNWALL || it == ICOL_SHADEWALL || it == ICOL_ROOF)
        if wallroof
            if j <= nlevurb
                if j >= snl[c] + 1 && j <= nlevurb - 1
                    tk_out[c, jj] = _tk_iface(thk[c,jj], thk[c,jj+1], z[c,jj], z[c,jj+1], zi[c,jj+1])
                elseif j == nlevurb
                    tk_out[c, jj] = thk[c, jj]
                end
            end
        else
            if j >= snl[c] + 1 && j <= nlevgrnd - 1
                tk_out[c, jj] = _tk_iface(thk[c,jj], thk[c,jj+1], z[c,jj], z[c,jj+1], zi[c,jj+1])
            elseif j == nlevgrnd
                tk_out[c, jj] = zero(eltype(tk_out))
            end
        end
    end
end

"""
    compute_soil_tk_interface!(tk_out, mask, itype, snl, z, zi, thk,
                               nlevsno, nlevurb, nlevgrnd, njlev)

Interface thermal conductivity for all active (column, level) pairs. Backend-agnostic.
"""
function compute_soil_tk_interface!(tk_out, mask, itype, snl, z, zi, thk,
                                    nlevsno::Int, nlevurb::Int, nlevgrnd::Int, njlev::Int)
    _launch!(_soil_tk_interface_kernel!, tk_out, mask, itype, snl, z, zi, thk,
             nlevsno, nlevurb, nlevgrnd; ndrange = (size(tk_out, 1), njlev))
end

# Per-layer soil + snow thermal conductivity thk (Farouki 1981 soil; Jordan1991 or
# Sturm1997 snow). One thread per (c, jj); j = jj - nlevsno. Physical constants are
# converted to the working element type `T` so the kernel carries no Float64 on a
# Float32-only backend. `snow_code`: 1=Jordan1991, 2=Sturm1997 (resolved on the host).
@kernel function _soil_snow_tk_kernel!(thk, bw, @Const(mask), @Const(landunit),
        @Const(itype), @Const(lun_itype), @Const(nlev_improad), @Const(dz),
        @Const(nbedrock), @Const(snl), @Const(t_soisno),
        @Const(h2osoi_liq), @Const(h2osoi_ice), @Const(excess_ice),
        @Const(watsat), @Const(tkmg), @Const(tkdry), @Const(frac_sno),
        nlevsno::Int, nlevgrnd::Int, nlevsoi::Int, snow_code::Int)
    c, jj = @index(Global, NTuple)
    @inbounds if mask[c]
        T = eltype(thk)
        j = jj - nlevsno
        if j >= 1
            l = landunit[c]
            it = itype[c]
            lit = lun_itype[l]
            is_special = (lit == ISTWET || lit == ISTICE ||
                          it == ICOL_SUNWALL || it == ICOL_SHADEWALL ||
                          it == ICOL_ROOF || it == ICOL_ROAD_IMPERV)
            if (!is_special) || (it == ICOL_ROAD_IMPERV && j > nlev_improad[l])
                satw = (h2osoi_liq[c, jj] / T(DENH2O) + h2osoi_ice[c, jj] / T(DENICE) +
                        excess_ice[c, j] / T(DENICE)) / (dz[c, jj] * watsat[c, j])
                satw = smooth_min(one(T), satw)
                if satw > T(1.0e-7)
                    if t_soisno[c, jj] >= T(TFRZ)
                        dke = smooth_max(zero(T), log10(satw) + one(T))
                    else
                        dke = satw
                    end
                    fl_denom = h2osoi_liq[c, jj] / (T(DENH2O) * dz[c, jj]) +
                               h2osoi_ice[c, jj] / (T(DENICE) * dz[c, jj]) +
                               excess_ice[c, j] / (T(DENICE) * dz[c, jj])
                    fl = (h2osoi_liq[c, jj] / (T(DENH2O) * dz[c, jj])) / fl_denom
                    dksat = tkmg[c, j] * T(TKWAT)^(fl * watsat[c, j]) *
                            T(TKICE)^((one(T) - fl) * watsat[c, j])
                    thk[c, jj] = dke * dksat + (one(T) - dke) * tkdry[c, j]
                else
                    thk[c, jj] = tkdry[c, j]
                end
                if j > nbedrock[c]
                    thk[c, jj] = T(THK_BEDROCK)
                end
            elseif lit == ISTICE
                thk[c, jj] = t_soisno[c, jj] < T(TFRZ) ? T(TKICE) : T(TKWAT)
            elseif lit == ISTWET
                if j > nlevsoi
                    thk[c, jj] = T(THK_BEDROCK)
                else
                    thk[c, jj] = t_soisno[c, jj] < T(TFRZ) ? T(TKICE) : T(TKWAT)
                end
            end
        end

        # Thermal conductivity of snow
        if snl[c] + 1 < 1 && j >= snl[c] + 1 && j <= 0
            denom = smooth_max(frac_sno[c], T(1.0e-6)) * smooth_max(dz[c, jj], T(1.0e-6))
            bwv = (h2osoi_ice[c, jj] + h2osoi_liq[c, jj]) / denom
            bw[c, jj] = bwv
            if snow_code == 1            # Jordan1991
                thk[c, jj] = T(TKAIR) + (T(7.75e-5) * bwv + T(1.105e-6) * bwv^2) *
                                        (T(TKICE) - T(TKAIR))
            elseif snow_code == 2        # Sturm1997
                if bwv <= T(156.0)
                    thk[c, jj] = T(0.023) + T(0.234) * (bwv / T(1000.0))
                else
                    thk[c, jj] = T(0.138) - T(1.01) * (bwv / T(1000.0)) +
                                 T(3.233) * (bwv / T(1000.0))^2
                end
            end
        end
    end
end

"""
    compute_soil_snow_tk!(thk, bw, mask, landunit, itype, lun_itype, nlev_improad, dz,
                          nbedrock, snl, t_soisno, h2osoi_liq, h2osoi_ice, excess_ice,
                          watsat, tkmg, tkdry, frac_sno, nlevsno, nlevgrnd, nlevsoi,
                          snow_thermal_cond_method)

Soil + snow per-layer thermal conductivity over all active (column, level) pairs.
Resolves the `snow_thermal_cond_method` String to an integer code on the host (kernels
cannot compare Strings), then launches the backend-agnostic kernel.
"""
function compute_soil_snow_tk!(thk, bw, mask, landunit, itype, lun_itype, nlev_improad,
        dz, nbedrock, snl, t_soisno, h2osoi_liq, h2osoi_ice, excess_ice,
        watsat, tkmg, tkdry, frac_sno, nlevsno::Int, nlevgrnd::Int, nlevsoi::Int,
        snow_thermal_cond_method::AbstractString)
    snow_code = snow_thermal_cond_method == "Jordan1991" ? 1 :
                snow_thermal_cond_method == "Sturm1997"  ? 2 :
                error("Unknown snow_thermal_cond_method: $(snow_thermal_cond_method)")
    _launch!(_soil_snow_tk_kernel!, thk, bw, mask, landunit, itype, lun_itype, nlev_improad,
             dz, nbedrock, snl, t_soisno, h2osoi_liq, h2osoi_ice, excess_ice,
             watsat, tkmg, tkdry, frac_sno, nlevsno, nlevgrnd, nlevsoi, snow_code;
             ndrange = (size(thk, 1), nlevsno + nlevgrnd))
end

# Urban per-layer thermal conductivity from urban params (thread per (c, j), j in 1:nlevurb).
@kernel function _urban_tk_kernel!(thk, @Const(mask), @Const(itype), @Const(landunit),
        @Const(tk_wall), @Const(tk_roof), @Const(tk_improad), @Const(nlev_improad), joff::Int)
    c, j = @index(Global, NTuple)
    @inbounds if mask[c]
        l = landunit[c]; jj = j + joff; it = itype[c]
        if it == ICOL_SUNWALL || it == ICOL_SHADEWALL
            thk[c, jj] = tk_wall[l, j]
        elseif it == ICOL_ROOF
            thk[c, jj] = tk_roof[l, j]
        elseif it == ICOL_ROAD_IMPERV && j <= nlev_improad[l]
            thk[c, jj] = tk_improad[l, j]
        end
    end
end
compute_urban_tk!(thk, mask, itype, landunit, tk_wall, tk_roof, tk_improad, nlev_improad,
                  joff::Int, nlevurb::Int) =
    _launch!(_urban_tk_kernel!, thk, mask, itype, landunit, tk_wall, tk_roof, tk_improad,
             nlev_improad, joff; ndrange = (size(thk, 1), nlevurb))

# Surface-water (h2osfc) thermal conductivity (thread per column).
@kernel function _tk_h2osfc_kernel!(tk_h2osfc, @Const(mask), @Const(h2osfc), @Const(z),
                                    @Const(thk), joff::Int)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(tk_h2osfc)
        zh2osfc = T(1.0e-3) * (T(0.5) * h2osfc[c])
        z1 = z[c, 1 + joff]; tk1 = thk[c, 1 + joff]
        tk_h2osfc[c] = T(TKWAT) * tk1 * (z1 + zh2osfc) / (T(TKWAT) * z1 + tk1 * zh2osfc)
    end
end
compute_tk_h2osfc!(tk_h2osfc, mask, h2osfc, z, thk, joff::Int) =
    _launch!(_tk_h2osfc_kernel!, tk_h2osfc, mask, h2osfc, z, thk, joff; ndrange = length(tk_h2osfc))

# Soil heat capacity, de Vries 1963 (thread per (c, j), j in 1:nlevgrnd).
@kernel function _soil_cv_kernel!(cv, @Const(mask), @Const(landunit), @Const(itype),
        @Const(lun_itype), @Const(nlev_improad), @Const(csol), @Const(watsat),
        @Const(dz), @Const(h2osoi_ice), @Const(h2osoi_liq), @Const(excess_ice),
        @Const(nbedrock), nlevsno::Int)
    c, j = @index(Global, NTuple)
    @inbounds if mask[c]
        T = eltype(cv); jj = j + nlevsno; l = landunit[c]; it = itype[c]; lit = lun_itype[l]
        is_special = (lit == ISTWET || lit == ISTICE || it == ICOL_SUNWALL ||
                      it == ICOL_SHADEWALL || it == ICOL_ROOF || it == ICOL_ROAD_IMPERV)
        if (!is_special) || (it == ICOL_ROAD_IMPERV && j > nlev_improad[l])
            cv[c, jj] = csol[c, j] * (one(T) - watsat[c, j]) * dz[c, jj] +
                (h2osoi_ice[c, jj] * T(CPICE) + h2osoi_liq[c, jj] * T(CPLIQ)) +
                excess_ice[c, j] * T(CPICE)
            if j > nbedrock[c]
                cv[c, jj] = T(CSOL_BEDROCK) * dz[c, jj]
            end
        elseif lit == ISTWET
            cv[c, jj] = h2osoi_ice[c, jj] * T(CPICE) + h2osoi_liq[c, jj] * T(CPLIQ)
            if j > nbedrock[c]
                cv[c, jj] = T(CSOL_BEDROCK) * dz[c, jj]
            end
        elseif lit == ISTICE
            cv[c, jj] = h2osoi_ice[c, jj] * T(CPICE) + h2osoi_liq[c, jj] * T(CPLIQ)
        end
    end
end
compute_soil_cv!(cv, mask, landunit, itype, lun_itype, nlev_improad, csol, watsat, dz,
                 h2osoi_ice, h2osoi_liq, excess_ice, nbedrock, nlevsno::Int, nlevgrnd::Int) =
    _launch!(_soil_cv_kernel!, cv, mask, landunit, itype, lun_itype, nlev_improad, csol,
             watsat, dz, h2osoi_ice, h2osoi_liq, excess_ice, nbedrock, nlevsno;
             ndrange = (size(cv, 1), nlevgrnd))

# Urban heat capacity (thread per (c, j), j in 1:nlevurb).
@kernel function _urban_cv_kernel!(cv, @Const(mask), @Const(itype), @Const(landunit),
        @Const(cv_wall), @Const(cv_roof), @Const(cv_improad), @Const(nlev_improad),
        @Const(dz), joff::Int)
    c, j = @index(Global, NTuple)
    @inbounds if mask[c]
        l = landunit[c]; jj = j + joff; it = itype[c]
        if it == ICOL_SUNWALL || it == ICOL_SHADEWALL
            cv[c, jj] = cv_wall[l, j] * dz[c, jj]
        elseif it == ICOL_ROOF
            cv[c, jj] = cv_roof[l, j] * dz[c, jj]
        elseif it == ICOL_ROAD_IMPERV && j <= nlev_improad[l]
            cv[c, jj] = cv_improad[l, j] * dz[c, jj]
        end
    end
end
compute_urban_cv!(cv, mask, itype, landunit, cv_wall, cv_roof, cv_improad, nlev_improad,
                  dz, joff::Int, nlevurb::Int) =
    _launch!(_urban_cv_kernel!, cv, mask, itype, landunit, cv_wall, cv_roof, cv_improad,
             nlev_improad, dz, joff; ndrange = (size(cv, 1), nlevurb))

# Add unresolved (no-layer) snow heat capacity to the top soil layer (thread per column).
@kernel function _unresolved_snow_cv_kernel!(cv, @Const(mask), @Const(h2osno_no_layers), joff::Int)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(cv)
        if h2osno_no_layers[c] > zero(T)
            cv[c, 1 + joff] += T(CPICE) * h2osno_no_layers[c]
        end
    end
end
add_unresolved_snow_cv!(cv, mask, h2osno_no_layers, joff::Int) =
    _launch!(_unresolved_snow_cv_kernel!, cv, mask, h2osno_no_layers, joff; ndrange = size(cv, 1))

# Snow-layer heat capacity (thread per (c, jj), jj in 1:nlevsno → j = jj - nlevsno).
@kernel function _snow_cv_kernel!(cv, @Const(mask), @Const(snl), @Const(frac_sno),
        @Const(h2osoi_liq), @Const(h2osoi_ice), nlevsno::Int)
    c, jj = @index(Global, NTuple)
    @inbounds if mask[c]
        T = eltype(cv); j = jj - nlevsno
        if snl[c] + 1 < 1 && j >= snl[c] + 1
            if frac_sno[c] > zero(T)
                cv[c, jj] = smooth_max(T(THIN_SFCLAYER),
                    (T(CPLIQ) * h2osoi_liq[c, jj] + T(CPICE) * h2osoi_ice[c, jj]) / frac_sno[c])
            else
                cv[c, jj] = T(THIN_SFCLAYER)
            end
        end
    end
end
compute_snow_cv!(cv, mask, snl, frac_sno, h2osoi_liq, h2osoi_ice, nlevsno::Int) =
    _launch!(_snow_cv_kernel!, cv, mask, snl, frac_sno, h2osoi_liq, h2osoi_ice, nlevsno;
             ndrange = (size(cv, 1), nlevsno))

function soil_therm_prop!(col::ColumnData, lun::LandunitData,
                          urbanparams::UrbanParamsData, temperature::TemperatureData,
                          waterstatebulk::WaterStateBulkData,
                          waterdiagbulk::WaterDiagnosticBulkData,
                          soilstate::SoilStateData,
                          mask_nolakec::BitVector, mask_urbanc::BitVector,
                          bounds_col::UnitRange{Int},
                          tk_out::Matrix{<:Real}, cv_out::Matrix{<:Real},
                          tk_h2osfc_out::Vector{<:Real})
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

    # Soil + snow thermal conductivity, Farouki 1981 (kernelized; thread per (c, jj)).
    compute_soil_snow_tk!(thk, bw, mask_nolakec, col.landunit, col.itype, lun.itype,
                          urbanparams.nlev_improad, col.dz, col.nbedrock, col.snl,
                          t_soisno, h2osoi_liq, h2osoi_ice, excess_ice,
                          watsat, tkmg, tkdry, frac_sno, nlevsno, nlevgrnd, nlevsoi,
                          varctl.snow_thermal_cond_method)

    # Urban per-layer thermal conductivity (kernelized).
    compute_urban_tk!(thk, mask_urbanc, col.itype, col.landunit, urbanparams.tk_wall,
                      urbanparams.tk_roof, urbanparams.tk_improad, urbanparams.nlev_improad,
                      joff, nlevurb)

    # Thermal conductivity at layer interface (kernelized; one thread per (c, jj)).
    compute_soil_tk_interface!(tk_out, mask_nolakec, col.itype, col.snl, col.z, col.zi,
                               thk, nlevsno, nlevurb, nlevgrnd, nlevsno + nlevmaxurbgrnd)

    # Surface-water thermal conductivity (kernelized).
    compute_tk_h2osfc!(tk_h2osfc_out, mask_nolakec, h2osfc, col.z, thk, joff)

    # Soil heat capacity, de Vries 1963 (kernelized).
    compute_soil_cv!(cv_out, mask_nolakec, col.landunit, col.itype, lun.itype,
                     urbanparams.nlev_improad, csol, watsat, col.dz, h2osoi_ice,
                     h2osoi_liq, excess_ice, col.nbedrock, nlevsno, nlevgrnd)

    # Urban heat capacity (kernelized).
    compute_urban_cv!(cv_out, mask_urbanc, col.itype, col.landunit, urbanparams.cv_wall,
                      urbanparams.cv_roof, urbanparams.cv_improad, urbanparams.nlev_improad,
                      col.dz, joff, nlevurb)

    # Unresolved (no-layer) snow heat capacity on the top soil layer (kernelized;
    # runs after the soil heat capacity, which it adds to).
    add_unresolved_snow_cv!(cv_out, mask_nolakec, h2osno_no_layers, joff)

    # Snow-layer heat capacity (kernelized).
    compute_snow_cv!(cv_out, mask_nolakec, col.snl, frac_sno, h2osoi_liq, h2osoi_ice, nlevsno)

    return nothing
end

# Longwave radiation emitted by ground / snow / soil / surface-water (thread per column).
# Constants converted to the working element type so no Float64 reaches a Float32 backend.
@kernel function _ground_lwrad_emit_kernel!(lwrad_emit, dlwrad_emit, lwrad_emit_snow,
        lwrad_emit_soil, lwrad_emit_h2osfc, @Const(mask), @Const(emg), @Const(t_grnd),
        @Const(t_soisno), @Const(t_h2osfc), @Const(snl), joff::Int)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(lwrad_emit)
        sb = T(SB)
        eg = emg[c]
        lwrad_emit[c]        = eg * sb * t_grnd[c]^4
        dlwrad_emit[c]       = T(4) * eg * sb * t_grnd[c]^3
        lwrad_emit_snow[c]   = eg * sb * t_soisno[c, snl[c] + 1 + joff]^4
        lwrad_emit_soil[c]   = eg * sb * t_soisno[c, 1 + joff]^4
        lwrad_emit_h2osfc[c] = eg * sb * t_h2osfc[c]^4
    end
end
function compute_ground_lwrad_emit!(lwrad_emit, dlwrad_emit, lwrad_emit_snow,
        lwrad_emit_soil, lwrad_emit_h2osfc, mask, emg, t_grnd, t_soisno, t_h2osfc, snl, joff::Int)
    _launch!(_ground_lwrad_emit_kernel!, lwrad_emit, dlwrad_emit, lwrad_emit_snow,
             lwrad_emit_soil, lwrad_emit_h2osfc, mask, emg, t_grnd, t_soisno, t_h2osfc, snl, joff;
             ndrange = length(lwrad_emit))
end

# --------------------------------------------------------------------------
# Device-view structs for the ground-heat-flux patch kernel. Each groups a set of
# same-typed (float) arrays so a kernel touching ~49 fields takes a few struct args
# instead of ~49 loose ones. Adapt-registered, so adapt(backend, s) moves the arrays
# to the device and the struct passes into a KernelAbstractions kernel (its fields are
# indexed inside). @kwdef gives name-based (mis-order-proof) construction.
# --------------------------------------------------------------------------
Base.@kwdef struct GnetPatchIn{V}      # per-patch float inputs (frac_veg_nosno is Int → loose arg)
    sabg::V; sabg_soil::V; sabg_snow::V
    dlrad::V; cgrnd::V
    eflx_sh_grnd::V; eflx_sh_snow::V; eflx_sh_soil::V; eflx_sh_h2osfc::V
    qflx_evap_soi::V; qflx_ev_snow::V; qflx_ev_soil::V; qflx_ev_h2osfc::V
    eflx_lwrad_net::V; qflx_tran_veg::V; wtcol::V
end
Base.@kwdef struct GnetColIn{V}        # per-column float inputs
    lwrad_emit::V; lwrad_emit_snow::V; lwrad_emit_soil::V; lwrad_emit_h2osfc::V
    dlwrad_emit::V; emg::V; forc_lwrad::V; htvp::V; frac_sno_eff::V
end
Base.@kwdef struct GnetLunIn{V}        # per-landunit float inputs
    eflx_wasteheat_lun::V; eflx_ventilation_lun::V; eflx_heat_from_ac_lun::V
    eflx_traffic_lun::V; wtlunit_roof::V
end
Base.@kwdef struct GnetOut{V}          # per-patch float outputs
    eflx_gnet::V; dgnetdT::V; sabg_chk::V
    eflx_wasteheat::V; eflx_ventilation::V; eflx_heat_from_ac::V; eflx_traffic::V; eflx_anthro::V
end
Adapt.@adapt_structure GnetPatchIn
Adapt.@adapt_structure GnetColIn
Adapt.@adapt_structure GnetLunIn
Adapt.@adapt_structure GnetOut

# Per-patch ground net heat flux + atomic scatter of hs/dhsdT/hs_soil/hs_h2osfc into
# columns. One thread per patch; c = column[p], l = landunit[p]. HVAP converted to the
# working type. `simple_build`/`prog_build`: building-temperature config resolved on host.
@kernel function _gnet_patch_kernel!(out::GnetOut, pin::GnetPatchIn, cin::GnetColIn,
        lin::GnetLunIn, @Const(mask), @Const(column), @Const(landunit), @Const(itype),
        @Const(urbpoi), @Const(frac_veg_nosno), hs, dhsdT, hs_soil, hs_h2osfc,
        simple_build::Bool, prog_build::Bool)
    p = @index(Global)
    @inbounds if mask[p]
        T = eltype(hs)
        c = column[p]; l = landunit[p]; wt = pin.wtcol[p]
        local eflx_gnet_soil, eflx_gnet_h2osfc
        if !urbpoi[l]
            base = pin.dlrad[p] + (one(T) - frac_veg_nosno[p]) * cin.emg[c] * cin.forc_lwrad[c]
            out.eflx_gnet[p] = pin.sabg[p] + base - cin.lwrad_emit[c] -
                (pin.eflx_sh_grnd[p] + pin.qflx_evap_soi[p] * cin.htvp[c])
            out.sabg_chk[p] = cin.frac_sno_eff[c] * pin.sabg_snow[p] +
                (one(T) - cin.frac_sno_eff[c]) * pin.sabg_soil[p]
            eflx_gnet_soil = pin.sabg_soil[p] + base - cin.lwrad_emit_soil[c] -
                (pin.eflx_sh_soil[p] + pin.qflx_ev_soil[p] * cin.htvp[c])
            eflx_gnet_h2osfc = pin.sabg_soil[p] + base - cin.lwrad_emit_h2osfc[c] -
                (pin.eflx_sh_h2osfc[p] + pin.qflx_ev_h2osfc[p] * cin.htvp[c])
        else
            it = itype[c]
            denom = one(T) - lin.wtlunit_roof[l]
            if it == ICOL_ROAD_PERV || it == ICOL_ROAD_IMPERV
                out.eflx_wasteheat[p] = lin.eflx_wasteheat_lun[l] / denom
                if simple_build
                    out.eflx_ventilation[p] = zero(T)
                elseif prog_build
                    out.eflx_ventilation[p] = lin.eflx_ventilation_lun[l] / denom
                end
                out.eflx_heat_from_ac[p] = lin.eflx_heat_from_ac_lun[l] / denom
                out.eflx_traffic[p] = lin.eflx_traffic_lun[l] / denom
            else
                out.eflx_wasteheat[p] = zero(T)
                out.eflx_ventilation[p] = zero(T)
                out.eflx_heat_from_ac[p] = zero(T)
                out.eflx_traffic[p] = zero(T)
            end
            out.eflx_gnet[p] = pin.sabg[p] + pin.dlrad[p] - pin.eflx_lwrad_net[p] -
                (pin.eflx_sh_grnd[p] + pin.qflx_evap_soi[p] * cin.htvp[c] +
                 pin.qflx_tran_veg[p] * T(HVAP)) +
                out.eflx_wasteheat[p] + out.eflx_heat_from_ac[p] +
                out.eflx_traffic[p] + out.eflx_ventilation[p]
            if simple_build
                out.eflx_anthro[p] = out.eflx_wasteheat[p] + out.eflx_traffic[p]
            end
            eflx_gnet_soil = out.eflx_gnet[p]
            eflx_gnet_h2osfc = out.eflx_gnet[p]
        end
        out.dgnetdT[p] = -pin.cgrnd[p] - cin.dlwrad_emit[c]
        _scatter_add!(hs,        c, out.eflx_gnet[p] * wt)
        _scatter_add!(dhsdT,     c, out.dgnetdT[p]   * wt)
        _scatter_add!(hs_soil,   c, eflx_gnet_soil   * wt)
        _scatter_add!(hs_h2osfc, c, eflx_gnet_h2osfc * wt)
    end
end

function compute_gnet_patch!(energyflux, solarabs, canopystate, waterfluxbulk, col, lun,
        patch_data, mask_nolakep, lwrad_emit, dlwrad_emit, lwrad_emit_snow, lwrad_emit_soil,
        lwrad_emit_h2osfc, emg, forc_lwrad, htvp, frac_sno_eff, hs, dhsdT, hs_soil, hs_h2osfc,
        simple_build::Bool, prog_build::Bool)
    pin = GnetPatchIn(;
        sabg = solarabs.sabg_patch, sabg_soil = solarabs.sabg_soil_patch,
        sabg_snow = solarabs.sabg_snow_patch, dlrad = energyflux.dlrad_patch,
        cgrnd = energyflux.cgrnd_patch, eflx_sh_grnd = energyflux.eflx_sh_grnd_patch,
        eflx_sh_snow = energyflux.eflx_sh_snow_patch, eflx_sh_soil = energyflux.eflx_sh_soil_patch,
        eflx_sh_h2osfc = energyflux.eflx_sh_h2osfc_patch,
        qflx_evap_soi = waterfluxbulk.wf.qflx_evap_soi_patch,
        qflx_ev_snow = waterfluxbulk.qflx_ev_snow_patch, qflx_ev_soil = waterfluxbulk.qflx_ev_soil_patch,
        qflx_ev_h2osfc = waterfluxbulk.qflx_ev_h2osfc_patch,
        eflx_lwrad_net = energyflux.eflx_lwrad_net_patch,
        qflx_tran_veg = waterfluxbulk.wf.qflx_tran_veg_patch, wtcol = patch_data.wtcol)
    cin = GnetColIn(; lwrad_emit, lwrad_emit_snow, lwrad_emit_soil, lwrad_emit_h2osfc,
        dlwrad_emit, emg, forc_lwrad, htvp, frac_sno_eff)
    lin = GnetLunIn(; eflx_wasteheat_lun = energyflux.eflx_wasteheat_lun,
        eflx_ventilation_lun = energyflux.eflx_ventilation_lun,
        eflx_heat_from_ac_lun = energyflux.eflx_heat_from_ac_lun,
        eflx_traffic_lun = energyflux.eflx_traffic_lun, wtlunit_roof = lun.wtlunit_roof)
    out = GnetOut(; eflx_gnet = energyflux.eflx_gnet_patch, dgnetdT = energyflux.dgnetdT_patch,
        sabg_chk = solarabs.sabg_chk_patch, eflx_wasteheat = energyflux.eflx_wasteheat_patch,
        eflx_ventilation = energyflux.eflx_ventilation_patch,
        eflx_heat_from_ac = energyflux.eflx_heat_from_ac_patch,
        eflx_traffic = energyflux.eflx_traffic_patch, eflx_anthro = energyflux.eflx_anthro_patch)
    backend = _kernel_backend(hs)
    np = length(mask_nolakep)
    np == 0 && return nothing
    _gnet_patch_kernel!(backend)(out, pin, cin, lin, mask_nolakep, patch_data.column,
        patch_data.landunit, col.itype, lun.urbpoi, canopystate.frac_veg_nosno_patch,
        hs, dhsdT, hs_soil, hs_h2osfc, simple_build, prog_build; ndrange = np)
    KernelAbstractions.synchronize(backend)
    return nothing
end

# SNICAR per-patch top-layer net heat flux + per-snow-layer solar; atomically scatters
# hs_top / hs_top_snow (per column) and sabg_lyr_col (per column, layer). Reuses the
# GnetPatchIn/GnetColIn device-view structs (its fields are a subset). One thread per patch.
@kernel function _gnet_snicar_kernel!(pin::GnetPatchIn, cin::GnetColIn, @Const(mask),
        @Const(column), @Const(landunit), @Const(urbpoi), @Const(frac_veg_nosno),
        @Const(snl), @Const(sabg_lyr), @Const(eflx_gnet),
        hs_top, hs_top_snow, sabg_lyr_col, nlevsno::Int)
    p = @index(Global)
    @inbounds if mask[p]
        T = eltype(hs_top)
        c = column[p]; l = landunit[p]; wt = pin.wtcol[p]; joff = nlevsno
        lyr_top = snl[c] + 1
        if !urbpoi[l]
            base = pin.dlrad[p] + (one(T) - frac_veg_nosno[p]) * cin.emg[c] * cin.forc_lwrad[c]
            sab = sabg_lyr[p, lyr_top + joff]
            eflx_gnet_top = sab + base - cin.lwrad_emit[c] -
                (pin.eflx_sh_grnd[p] + pin.qflx_evap_soi[p] * cin.htvp[c])
            _scatter_add!(hs_top, c, eflx_gnet_top * wt)
            eflx_gnet_snow = sab + base - cin.lwrad_emit_snow[c] -
                (pin.eflx_sh_snow[p] + pin.qflx_ev_snow[p] * cin.htvp[c])
            _scatter_add!(hs_top_snow, c, eflx_gnet_snow * wt)
            for j in lyr_top:1
                _scatter_add!(sabg_lyr_col, c, j + nlevsno, sabg_lyr[p, j + joff] * wt)
            end
        else
            g = eflx_gnet[p]
            _scatter_add!(hs_top, c, g * wt)
            _scatter_add!(hs_top_snow, c, g * wt)
            _scatter_add!(sabg_lyr_col, c, lyr_top + nlevsno, pin.sabg[p] * wt)
        end
    end
end

function compute_gnet_snicar!(energyflux, solarabs, canopystate, waterfluxbulk, col, lun,
        patch_data, mask_nolakep, lwrad_emit, dlwrad_emit, lwrad_emit_snow, lwrad_emit_soil,
        lwrad_emit_h2osfc, emg, forc_lwrad, htvp, frac_sno_eff, hs_top, hs_top_snow,
        sabg_lyr_col, nlevsno::Int)
    pin = GnetPatchIn(;
        sabg = solarabs.sabg_patch, sabg_soil = solarabs.sabg_soil_patch,
        sabg_snow = solarabs.sabg_snow_patch, dlrad = energyflux.dlrad_patch,
        cgrnd = energyflux.cgrnd_patch, eflx_sh_grnd = energyflux.eflx_sh_grnd_patch,
        eflx_sh_snow = energyflux.eflx_sh_snow_patch, eflx_sh_soil = energyflux.eflx_sh_soil_patch,
        eflx_sh_h2osfc = energyflux.eflx_sh_h2osfc_patch,
        qflx_evap_soi = waterfluxbulk.wf.qflx_evap_soi_patch,
        qflx_ev_snow = waterfluxbulk.qflx_ev_snow_patch, qflx_ev_soil = waterfluxbulk.qflx_ev_soil_patch,
        qflx_ev_h2osfc = waterfluxbulk.qflx_ev_h2osfc_patch,
        eflx_lwrad_net = energyflux.eflx_lwrad_net_patch,
        qflx_tran_veg = waterfluxbulk.wf.qflx_tran_veg_patch, wtcol = patch_data.wtcol)
    cin = GnetColIn(; lwrad_emit, lwrad_emit_snow, lwrad_emit_soil, lwrad_emit_h2osfc,
        dlwrad_emit, emg, forc_lwrad, htvp, frac_sno_eff)
    backend = _kernel_backend(hs_top)
    np = length(mask_nolakep)
    np == 0 && return nothing
    _gnet_snicar_kernel!(backend)(pin, cin, mask_nolakep, patch_data.column,
        patch_data.landunit, lun.urbpoi, canopystate.frac_veg_nosno_patch, col.snl,
        solarabs.sabg_lyr_patch, energyflux.eflx_gnet_patch,
        hs_top, hs_top_snow, sabg_lyr_col, nlevsno; ndrange = np)
    KernelAbstractions.synchronize(backend)
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
        urbanparams::UrbanParamsData, forc_lwrad::Vector{<:Real},
        mask_nolakec::BitVector, mask_nolakep::BitVector,
        bounds_col::UnitRange{Int}, bounds_patch::UnitRange{Int},
        hs_h2osfc::Vector{<:Real}, hs_top_snow::Vector{<:Real},
        hs_soil::Vector{<:Real}, hs_top::Vector{<:Real},
        dhsdT::Vector{<:Real}, sabg_lyr_col::Matrix{<:Real})

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
    FT = eltype(temperature.t_soisno_col)
    lwrad_emit = zeros(FT, nc)
    dlwrad_emit = zeros(FT, nc)
    lwrad_emit_snow = zeros(FT, nc)
    lwrad_emit_soil = zeros(FT, nc)
    lwrad_emit_h2osfc_arr = zeros(FT, nc)
    hs = zeros(FT, nc)

    compute_ground_lwrad_emit!(lwrad_emit, dlwrad_emit, lwrad_emit_snow, lwrad_emit_soil,
                               lwrad_emit_h2osfc_arr, mask_nolakec, emg, t_grnd, t_soisno,
                               t_h2osfc, snl, joff)

    hs_soil .= 0.0
    hs_h2osfc .= 0.0
    hs .= 0.0
    dhsdT .= 0.0

    # Per-patch ground net heat flux + atomic patch→column scatter (kernelized).
    compute_gnet_patch!(energyflux, solarabs, canopystate, waterfluxbulk, col, lun,
        patch_data, mask_nolakep, lwrad_emit, dlwrad_emit, lwrad_emit_snow, lwrad_emit_soil,
        lwrad_emit_h2osfc_arr, emg, forc_lwrad, htvp, frac_sno_eff, hs, dhsdT, hs_soil,
        hs_h2osfc, is_simple_build_temp(), is_prog_build_temp())

    # SNICAR: sabg_lyr_col and hs_top
    sabg_lyr_col .= 0.0
    hs_top .= 0.0
    hs_top_snow .= 0.0

    # SNICAR top-layer net heat flux + per-layer solar (kernelized; atomic scatter).
    compute_gnet_snicar!(energyflux, solarabs, canopystate, waterfluxbulk, col, lun,
        patch_data, mask_nolakep, lwrad_emit, dlwrad_emit, lwrad_emit_snow, lwrad_emit_soil,
        lwrad_emit_h2osfc_arr, emg, forc_lwrad, htvp, frac_sno_eff, hs_top, hs_top_snow,
        sabg_lyr_col, nlevsno)

    return nothing
end

# =========================================================================
# compute_heat_diff_flux_and_factor!
# =========================================================================
# Per-(column, level) heat-diffusion flux fn and time-step factor fact. Pure elementwise
# (no scatter); one thread per (c, jj). Constants (CNFAC, CAPR) and the dtime scalar are at
# the working element type. simple_build: building-temp config resolved on host.
@kernel function _heat_diff_kernel!(fact, fn, @Const(mask), @Const(itype), @Const(snl),
        @Const(landunit), @Const(z), @Const(zi), @Const(dz), @Const(tk), @Const(cv),
        @Const(t_soisno), @Const(eflx_bot), @Const(t_building), @Const(t_roof_inner),
        @Const(t_sunw_inner), @Const(t_shdw_inner), dtime, nlevsno::Int, nlevurb::Int,
        nlevgrnd::Int, simple_build::Bool)
    c, jj = @index(Global, NTuple)
    @inbounds if mask[c]
        T = eltype(fact)
        j = jj - nlevsno
        l = landunit[c]
        it = itype[c]
        cnfac = T(CNFAC)
        wallroof = (it == ICOL_SUNWALL || it == ICOL_SHADEWALL || it == ICOL_ROOF)
        if wallroof && j <= nlevurb
            if j >= snl[c] + 1
                if j == snl[c] + 1
                    fact[c, jj] = dtime / cv[c, jj]
                    fn[c, jj] = tk[c, jj] * (t_soisno[c, jj+1] - t_soisno[c, jj]) / (z[c, jj+1] - z[c, jj])
                elseif j <= nlevurb - 1
                    fact[c, jj] = dtime / cv[c, jj]
                    fn[c, jj] = tk[c, jj] * (t_soisno[c, jj+1] - t_soisno[c, jj]) / (z[c, jj+1] - z[c, jj])
                elseif j == nlevurb
                    fact[c, jj] = dtime / cv[c, jj]
                    if simple_build
                        fn[c, jj] = tk[c, jj] * (t_building[l] - cnfac * t_soisno[c, jj]) / (zi[c, jj+1] - z[c, jj])
                    else
                        tinner = it == ICOL_SUNWALL ? t_sunw_inner[l] :
                                 it == ICOL_SHADEWALL ? t_shdw_inner[l] : t_roof_inner[l]
                        fn[c, jj] = tk[c, jj] * (tinner - cnfac * t_soisno[c, jj]) / (zi[c, jj+1] - z[c, jj])
                    end
                end
            end
        elseif !wallroof && j <= nlevgrnd
            if j >= snl[c] + 1
                if j == snl[c] + 1
                    fact[c, jj] = dtime / cv[c, jj] * dz[c, jj] /
                        (T(0.5) * (z[c, jj] - zi[c, jj] + T(CAPR) * (z[c, jj+1] - zi[c, jj])))
                    fn[c, jj] = tk[c, jj] * (t_soisno[c, jj+1] - t_soisno[c, jj]) / (z[c, jj+1] - z[c, jj])
                elseif j <= nlevgrnd - 1
                    fact[c, jj] = dtime / cv[c, jj]
                    fn[c, jj] = tk[c, jj] * (t_soisno[c, jj+1] - t_soisno[c, jj]) / (z[c, jj+1] - z[c, jj])
                elseif j == nlevgrnd
                    fact[c, jj] = dtime / cv[c, jj]
                    fn[c, jj] = eflx_bot[c]
                end
            end
        end
    end
end
function compute_heat_diff!(fact, fn, mask, itype, snl, landunit, z, zi, dz, tk, cv,
        t_soisno, eflx_bot, t_building, t_roof_inner, t_sunw_inner, t_shdw_inner,
        dtime, nlevsno::Int, nlevurb::Int, nlevgrnd::Int, nlevmaxurbgrnd::Int, simple_build::Bool)
    dt = convert(eltype(fact), dtime)
    _launch!(_heat_diff_kernel!, fact, fn, mask, itype, snl, landunit, z, zi, dz, tk, cv,
        t_soisno, eflx_bot, t_building, t_roof_inner, t_sunw_inner, t_shdw_inner,
        dt, nlevsno, nlevurb, nlevgrnd, simple_build;
        ndrange = (size(fact, 1), nlevsno + nlevmaxurbgrnd))
end

function compute_heat_diff_flux_and_factor!(
        col::ColumnData, lun::LandunitData,
        temperature::TemperatureData, energyflux::EnergyFluxData,
        urbanparams::UrbanParamsData,
        mask_nolakec::BitVector, bounds_col::UnitRange{Int},
        dtime::Real,
        tk::Matrix{<:Real}, cv::Matrix{<:Real},
        fn::Matrix{<:Real}, fact::Matrix{<:Real})

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

    # Heat-diffusion flux + factor (kernelized; one thread per (c, jj)).
    compute_heat_diff!(fact, fn, mask_nolakec, col.itype, col.snl, col.landunit, col.z,
        col.zi, col.dz, tk, cv, t_soisno, eflx_bot, t_building, t_roof_inner, t_sunw_inner,
        t_shdw_inner, dtime, nlevsno, nlevurb, nlevgrnd, nlevmaxurbgrnd, is_simple_build_temp())

    return nothing
end

# =========================================================================
# set_rhs_vec! — Combines snow, SSW, and soil RHS vectors
# =========================================================================
function set_rhs_vec!(col::ColumnData, lun::LandunitData,
                      temperature::TemperatureData,
                      waterdiagbulk::WaterDiagnosticBulkData,
                      mask_nolakec::BitVector, bounds_col::UnitRange{Int},
                      dtime::Real,
                      hs_h2osfc::Vector{<:Real}, hs_top_snow::Vector{<:Real},
                      hs_soil::Vector{<:Real}, hs_top::Vector{<:Real},
                      dhsdT::Vector{<:Real}, sabg_lyr_col::Matrix{<:Real},
                      tk::Matrix{<:Real}, tk_h2osfc::Vector{<:Real},
                      fact::Matrix{<:Real}, fn::Matrix{<:Real},
                      c_h2osfc::Vector{<:Real}, dz_h2osfc::Vector{<:Real},
                      rvector::Matrix{<:Real})

    nlevsno = varpar.nlevsno
    nlevgrnd = varpar.nlevgrnd
    nlevurb = varpar.nlevurb
    nlevmaxurbgrnd = varpar.nlevmaxurbgrnd
    joff = nlevsno

    t_soisno = temperature.t_soisno_col
    t_h2osfc = temperature.t_h2osfc_col
    frac_h2osfc = waterdiagbulk.frac_h2osfc_col
    frac_sno_eff = waterdiagbulk.frac_sno_eff_col

    FT = eltype(t_soisno)
    rvector .= FT(NaN)

    fn_h2osfc = zeros(FT, length(bounds_col))

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
                     dtime::Real, nband::Int,
                     dhsdT::Vector{<:Real}, tk::Matrix{<:Real},
                     tk_h2osfc::Vector{<:Real}, fact::Matrix{<:Real},
                     c_h2osfc::Vector{<:Real}, dz_h2osfc::Vector{<:Real},
                     bmatrix::Array{<:Real,3})

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
                              dtime::Real, dhsdT::Vector{<:Real})

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
            rho_avg = z_avg > 0.0 ? smooth_min(800.0, h2osno_total / z_avg) : 200.0

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
                            dtime::Real, dhsdT::Vector{<:Real})

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
    FT = eltype(temperature.t_soisno_col)
    hm = zeros(FT, nc, nlevsno + nlevmaxurbgrnd)
    xm = zeros(FT, nc, nlevsno + nlevmaxurbgrnd)
    xm2 = zeros(FT, nc, nlevsno + nlevmaxurbgrnd)
    wice0 = zeros(FT, nc, nlevsno + nlevmaxurbgrnd)
    wliq0 = zeros(FT, nc, nlevsno + nlevmaxurbgrnd)
    wexice0 = zeros(FT, nc, nlevsno + nlevmaxurbgrnd)
    wmass0 = zeros(FT, nc, nlevsno + nlevmaxurbgrnd)
    supercool = zeros(FT, nc, nlevmaxurbgrnd)
    tinc = zeros(FT, nc, nlevsno + nlevmaxurbgrnd)

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
                            h2osno_no_layers[c] = smooth_max(0.0, temp1 - xm[c, jj])
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
                            qflx_snomelt[c] = smooth_max(0.0, temp1 - h2osno_no_layers[c]) / dtime
                            xmf[c] = HFUS * qflx_snomelt[c]
                            qflx_snow_drain[c] = qflx_snomelt[c]
                        end
                    end

                    heatr = 0.0
                    if xm[c, jj] > 0.0
                        h2osoi_ice[c, jj] = smooth_max(0.0, wice0[c, jj] - xm[c, jj])
                        heatr = hm[c, jj] - HFUS * (wice0[c, jj] - h2osoi_ice[c, jj]) / dtime
                        xm2[c, jj] = xm[c, jj] - h2osoi_ice[c, jj]
                        if h2osoi_ice[c, jj] == 0.0
                            if wexice0[c, jj] >= 0.0 && xm2[c, jj] > 0.0 && j >= 2
                                excess_ice[c, j] = smooth_max(0.0, wexice0[c, jj] - xm2[c, jj])
                                heatr = hm[c, jj] - HFUS * (wexice0[c, jj] - excess_ice[c, j] +
                                        wice0[c, jj] - h2osoi_ice[c, jj]) / dtime
                            end
                        end
                    elseif xm[c, jj] < 0.0
                        if j <= 0
                            h2osoi_ice[c, jj] = smooth_min(wmass0[c, jj], wice0[c, jj] - xm[c, jj])
                        else
                            if wmass0[c, jj] - wexice0[c, jj] < supercool[c, j]
                                h2osoi_ice[c, jj] = 0.0
                            else
                                h2osoi_ice[c, jj] = smooth_min(wmass0[c, jj] - wexice0[c, jj] - supercool[c, j],
                                    wice0[c, jj] - xm[c, jj])
                            end
                        end
                        heatr = hm[c, jj] - HFUS * (wice0[c, jj] - h2osoi_ice[c, jj]) / dtime
                    end

                    ei_val = j >= 1 ? excess_ice[c, j] : 0.0
                    h2osoi_liq[c, jj] = smooth_max(0.0, wmass0[c, jj] - h2osoi_ice[c, jj] - ei_val)

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
                        exice_subs[c, j] = smooth_max(0.0, (wexice0[c, jj] - excess_ice[c, j]) / DENICE)
                    else
                        xmf[c] += HFUS * (wice0[c, jj] - h2osoi_ice[c, jj]) / dtime
                    end

                    if imelt[c, jj] == 1 && j < 1
                        qflx_snomelt_lyr[c, jj] = smooth_max(0.0, wice0[c, jj] - h2osoi_ice[c, jj]) / dtime
                        qflx_snomelt[c] += qflx_snomelt_lyr[c, jj]
                        waterdiagbulk.snomelt_accum_col[c] += qflx_snomelt_lyr[c, jj] * dtime * 1.0e-3
                    end
                    if imelt[c, jj] == 2 && j < 1
                        qflx_snofrz_lyr[c, jj] = smooth_max(0.0, h2osoi_ice[c, jj] - wice0[c, jj]) / dtime
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
                       t_building_max::Vector{<:Real},
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
