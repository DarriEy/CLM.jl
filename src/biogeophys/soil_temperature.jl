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
# Snow-layer RHS. One thread per (column, snow level s = 1..nlevsno); j = s - nlevsno
# (Fortran -nlevsno+1..0) and rv_idx = jj = s (j + nlevsno). CNFAC at the working type.
@kernel function _rhs_snow_kernel!(rvector, @Const(mask), @Const(itype), @Const(snl),
        @Const(t_soisno), @Const(fact), @Const(fn), @Const(dhsdT), @Const(hs_top_snow),
        @Const(hs_top), @Const(sabg_lyr_col), nlevsno::Int)
    c, s = @index(Global, NTuple)
    @inbounds if mask[c]
        T = eltype(rvector)
        cnfac = T(CNFAC)
        j = s - nlevsno
        jj = s
        it = itype[c]
        hs_top_lev = hs_top_snow[c]
        if it == ICOL_SUNWALL || it == ICOL_SHADEWALL || it == ICOL_ROOF
            hs_top_lev = hs_top[c]
        end
        if j == snl[c] + 1
            rvector[c, s] = t_soisno[c, jj] + fact[c, jj] * (hs_top_lev -
                dhsdT[c] * t_soisno[c, jj] + cnfac * fn[c, jj])
        elseif j > snl[c] + 1
            rvector[c, s] = t_soisno[c, jj] + cnfac * fact[c, jj] * (fn[c, jj] - fn[c, jj - 1])
            rvector[c, s] += fact[c, jj] * sabg_lyr_col[c, s]
        end
    end
end

# Standing surface water RHS (one thread per column). Also fills fn_h2osfc[c] for the
# surface-water correction kernel below. dtime arrives converted to the working type.
@kernel function _rhs_ssw_kernel!(rvector, @Const(mask), @Const(t_soisno), @Const(t_h2osfc),
        @Const(z), @Const(tk_h2osfc), @Const(dz_h2osfc), @Const(c_h2osfc),
        @Const(hs_h2osfc), @Const(dhsdT), fn_h2osfc, dtime, nlevsno::Int)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(rvector)
        cnfac = T(CNFAC)
        joff = nlevsno
        dzm = T(0.5) * dz_h2osfc[c] + z[c, 1 + joff]
        fnh = tk_h2osfc[c] * (t_soisno[c, 1 + joff] - t_h2osfc[c]) / dzm
        fn_h2osfc[c] = fnh
        rvector[c, nlevsno + 1] = t_h2osfc[c] + (dtime / c_h2osfc[c]) *
            (hs_h2osfc[c] - dhsdT[c] * t_h2osfc[c] + cnfac * fnh)
    end
end

# Soil-layer RHS — urban non-road (wall/roof). One thread per (column, soil level j=1..nlevurb).
@kernel function _rhs_soil_urban_kernel!(rvector, @Const(mask), @Const(itype), @Const(snl),
        @Const(t_soisno), @Const(fact), @Const(fn), @Const(dhsdT), @Const(hs_top),
        @Const(sabg_lyr_col), nlevsno::Int, nlevurb::Int)
    c, j = @index(Global, NTuple)
    @inbounds if mask[c]
        it = itype[c]
        if it == ICOL_SUNWALL || it == ICOL_SHADEWALL || it == ICOL_ROOF
            T = eltype(rvector)
            cnfac = T(CNFAC)
            jj = j + nlevsno
            rv_idx = j + nlevsno + 1
            if j == snl[c] + 1
                rvector[c, rv_idx] = t_soisno[c, jj] + fact[c, jj] * (hs_top[c] -
                    dhsdT[c] * t_soisno[c, jj] + cnfac * fn[c, jj])
            elseif j <= nlevurb - 1
                rvector[c, rv_idx] = t_soisno[c, jj] + cnfac * fact[c, jj] * (fn[c, jj] - fn[c, jj - 1])
                if j == 1
                    rvector[c, rv_idx] += fact[c, jj] * sabg_lyr_col[c, j + nlevsno]
                end
            elseif j == nlevurb
                rvector[c, rv_idx] = t_soisno[c, jj] + fact[c, jj] * (fn[c, jj] - cnfac * fn[c, jj - 1])
            end
        end
    end
end

# Soil-layer RHS — non-urban and urban road. One thread per (column, soil level j=1..nlevgrnd).
@kernel function _rhs_soil_kernel!(rvector, @Const(mask), @Const(itype), @Const(snl),
        @Const(landunit), @Const(urbpoi), @Const(t_soisno), @Const(fact), @Const(fn),
        @Const(dhsdT), @Const(frac_sno_eff), @Const(hs_top_snow), @Const(hs_soil),
        @Const(sabg_lyr_col), nlevsno::Int, nlevgrnd::Int)
    c, j = @index(Global, NTuple)
    @inbounds if mask[c]
        l = landunit[c]
        it = itype[c]
        if !urbpoi[l] || it == ICOL_ROAD_IMPERV || it == ICOL_ROAD_PERV
            T = eltype(rvector)
            cnfac = T(CNFAC)
            jj = j + nlevsno
            rv_idx = j + nlevsno + 1
            if j == snl[c] + 1
                rvector[c, rv_idx] = t_soisno[c, jj] + fact[c, jj] * (hs_top_snow[c] -
                    dhsdT[c] * t_soisno[c, jj] + cnfac * fn[c, jj])
            elseif j == 1
                rvector[c, rv_idx] = t_soisno[c, jj] + fact[c, jj] *
                    ((one(T) - frac_sno_eff[c]) * (hs_soil[c] - dhsdT[c] * t_soisno[c, jj]) +
                     cnfac * (fn[c, jj] - frac_sno_eff[c] * fn[c, jj - 1]))
                rvector[c, rv_idx] += frac_sno_eff[c] * fact[c, jj] * sabg_lyr_col[c, j + nlevsno]
            elseif j <= nlevgrnd - 1
                rvector[c, rv_idx] = t_soisno[c, jj] + cnfac * fact[c, jj] * (fn[c, jj] - fn[c, jj - 1])
            elseif j == nlevgrnd
                rvector[c, rv_idx] = t_soisno[c, jj] - cnfac * fact[c, jj] * fn[c, jj - 1] +
                    fact[c, jj] * fn[c, jj]
            end
        end
    end
end

# Surface-water correction to soil layer 1 (one thread per column). Reads fn_h2osfc filled
# by the SSW kernel and the rvector[soil-1] value written by the soil kernel; must run last.
@kernel function _rhs_h2osfc_corr_kernel!(rvector, @Const(mask), @Const(frac_h2osfc),
        @Const(fact), @Const(dhsdT), @Const(t_soisno), @Const(hs_soil),
        @Const(fn_h2osfc), nlevsno::Int)
    c = @index(Global)
    @inbounds if mask[c]
        if frac_h2osfc[c] != 0
            T = eltype(rvector)
            cnfac = T(CNFAC)
            jj = 1 + nlevsno
            rv_idx = 1 + nlevsno + 1
            rvector[c, rv_idx] -= frac_h2osfc[c] * fact[c, jj] *
                ((hs_soil[c] - dhsdT[c] * t_soisno[c, jj]) + cnfac * fn_h2osfc[c])
        end
    end
end

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

    t_soisno = temperature.t_soisno_col
    t_h2osfc = temperature.t_h2osfc_col
    frac_h2osfc = waterdiagbulk.frac_h2osfc_col
    frac_sno_eff = waterdiagbulk.frac_sno_eff_col

    FT = eltype(t_soisno)
    rvector .= FT(NaN)
    nc = size(rvector, 1)
    fn_h2osfc = zeros(FT, nc)
    dt = convert(FT, dtime)

    # Snow-layer RHS (one thread per (column, snow level)).
    _launch!(_rhs_snow_kernel!, rvector, mask_nolakec, col.itype, col.snl, t_soisno,
        fact, fn, dhsdT, hs_top_snow, hs_top, sabg_lyr_col, nlevsno; ndrange = (nc, nlevsno))
    # Standing surface water RHS (also fills fn_h2osfc for the correction below).
    _launch!(_rhs_ssw_kernel!, rvector, mask_nolakec, t_soisno, t_h2osfc, col.z, tk_h2osfc,
        dz_h2osfc, c_h2osfc, hs_h2osfc, dhsdT, fn_h2osfc, dt, nlevsno; ndrange = nc)
    # Soil RHS — urban non-road (wall/roof).
    _launch!(_rhs_soil_urban_kernel!, rvector, mask_nolakec, col.itype, col.snl, t_soisno,
        fact, fn, dhsdT, hs_top, sabg_lyr_col, nlevsno, nlevurb; ndrange = (nc, nlevurb))
    # Soil RHS — non-urban and urban road.
    _launch!(_rhs_soil_kernel!, rvector, mask_nolakec, col.itype, col.snl, col.landunit,
        lun.urbpoi, t_soisno, fact, fn, dhsdT, frac_sno_eff, hs_top_snow, hs_soil,
        sabg_lyr_col, nlevsno, nlevgrnd; ndrange = (nc, nlevgrnd))
    # Surface-water correction to soil layer 1 (reads fn_h2osfc; must run last).
    _launch!(_rhs_h2osfc_corr_kernel!, rvector, mask_nolakec, frac_h2osfc, fact, dhsdT,
        t_soisno, hs_soil, fn_h2osfc, nlevsno; ndrange = nc)

    return nothing
end

# =========================================================================
# set_matrix! — Assembles band diagonal matrix from submatrices
# =========================================================================
# Snow submatrix. One thread per (column, snow level s=1..nlevsno); j=s-nlevsno, bm_idx=jj=s.
# omc = 1-CNFAC folded once at the working type (byte-identical to inline (1.0-CNFAC) on Float64).
@kernel function _mat_snow_kernel!(bmatrix, @Const(mask), @Const(snl), @Const(z),
        @Const(tk), @Const(fact), @Const(dhsdT), nlevsno::Int)
    c, s = @index(Global, NTuple)
    @inbounds if mask[c]
        T = eltype(bmatrix)
        omc = one(T) - T(CNFAC)
        j = s - nlevsno
        jj = s
        bm_idx = s
        if j >= snl[c] + 1
            dzp = z[c, jj + 1] - z[c, jj]
            if j == snl[c] + 1
                bmatrix[c, 3, bm_idx] = one(T) + omc * fact[c, jj] * tk[c, jj] / dzp -
                    fact[c, jj] * dhsdT[c]
            else
                dzm = z[c, jj] - z[c, jj - 1]
                bmatrix[c, 4, bm_idx] = -omc * fact[c, jj] * tk[c, jj - 1] / dzm
                bmatrix[c, 3, bm_idx] = one(T) + omc * fact[c, jj] *
                    (tk[c, jj] / dzp + tk[c, jj - 1] / dzm)
            end
            if j != 0
                bmatrix[c, 2, bm_idx] = -omc * fact[c, jj] * tk[c, jj] / dzp
            else
                # snow-soil coupling (band 1)
                bmatrix[c, 1, bm_idx] = -omc * fact[c, jj] * tk[c, jj] / dzp
            end
        end
    end
end

# Standing surface water submatrix (one thread per column). dtime arrives converted.
@kernel function _mat_ssw_kernel!(bmatrix, @Const(mask), @Const(z), @Const(tk_h2osfc),
        @Const(dz_h2osfc), @Const(c_h2osfc), @Const(dhsdT), dtime, nlevsno::Int)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(bmatrix)
        omc = one(T) - T(CNFAC)
        joff = nlevsno
        ssw_idx = nlevsno + 1
        dzm = T(0.5) * dz_h2osfc[c] + z[c, 1 + joff]
        bmatrix[c, 3, ssw_idx] = one(T) + omc * (dtime / c_h2osfc[c]) *
            tk_h2osfc[c] / dzm - (dtime / c_h2osfc[c]) * dhsdT[c]
        bmatrix[c, 2, ssw_idx] = -omc * (dtime / c_h2osfc[c]) * tk_h2osfc[c] / dzm
    end
end

# Soil submatrix — urban non-road (wall/roof). One thread per (column, soil level j=1..nlevurb).
@kernel function _mat_soil_urban_kernel!(bmatrix, @Const(mask), @Const(itype), @Const(snl),
        @Const(z), @Const(zi), @Const(tk), @Const(fact), @Const(dhsdT),
        nlevsno::Int, nlevurb::Int)
    c, j = @index(Global, NTuple)
    @inbounds if mask[c]
        it = itype[c]
        if it == ICOL_SUNWALL || it == ICOL_SHADEWALL || it == ICOL_ROOF
            T = eltype(bmatrix)
            omc = one(T) - T(CNFAC)
            jj = j + nlevsno
            bm_idx = j + nlevsno + 1
            if j == snl[c] + 1
                dzp = z[c, jj + 1] - z[c, jj]
                bmatrix[c, 3, bm_idx] = one(T) + omc * fact[c, jj] * tk[c, jj] / dzp -
                    fact[c, jj] * dhsdT[c]
                bmatrix[c, 2, bm_idx] = -omc * fact[c, jj] * tk[c, jj] / dzp
            elseif j <= nlevurb - 1
                dzm = z[c, jj] - z[c, jj - 1]
                dzp = z[c, jj + 1] - z[c, jj]
                if j != 1
                    bmatrix[c, 4, bm_idx] = -omc * fact[c, jj] * tk[c, jj - 1] / dzm
                else
                    bmatrix[c, 5, bm_idx] = -omc * fact[c, jj] * tk[c, jj - 1] / dzm
                end
                bmatrix[c, 3, bm_idx] = one(T) + omc * fact[c, jj] *
                    (tk[c, jj] / dzp + tk[c, jj - 1] / dzm)
                bmatrix[c, 2, bm_idx] = -omc * fact[c, jj] * tk[c, jj] / dzp
            elseif j == nlevurb
                dzm = z[c, jj] - z[c, jj - 1]
                dzp = zi[c, jj + 1] - z[c, jj]
                bmatrix[c, 4, bm_idx] = -omc * fact[c, jj] * tk[c, jj - 1] / dzm
                bmatrix[c, 3, bm_idx] = one(T) + omc * fact[c, jj] *
                    (tk[c, jj - 1] / dzm + tk[c, jj] / dzp)
            end
        end
    end
end

# Soil submatrix — non-urban and urban road. One thread per (column, soil level j=1..nlevgrnd).
@kernel function _mat_soil_kernel!(bmatrix, @Const(mask), @Const(itype), @Const(snl),
        @Const(landunit), @Const(urbpoi), @Const(z), @Const(tk), @Const(fact),
        @Const(dhsdT), @Const(frac_sno_eff), nlevsno::Int, nlevgrnd::Int)
    c, j = @index(Global, NTuple)
    @inbounds if mask[c]
        l = landunit[c]
        it = itype[c]
        if it == ICOL_ROAD_IMPERV || it == ICOL_ROAD_PERV || !urbpoi[l]
            T = eltype(bmatrix)
            omc = one(T) - T(CNFAC)
            jj = j + nlevsno
            bm_idx = j + nlevsno + 1
            if j == snl[c] + 1
                dzp = z[c, jj + 1] - z[c, jj]
                bmatrix[c, 3, bm_idx] = one(T) + omc * fact[c, jj] * tk[c, jj] / dzp -
                    fact[c, jj] * dhsdT[c]
                bmatrix[c, 2, bm_idx] = -omc * fact[c, jj] * tk[c, jj] / dzp
            elseif j == 1
                dzm = z[c, jj] - z[c, jj - 1]
                dzp = z[c, jj + 1] - z[c, jj]
                bmatrix[c, 2, bm_idx] = -omc * fact[c, jj] * tk[c, jj] / dzp
                bmatrix[c, 3, bm_idx] = one(T) + omc * fact[c, jj] *
                    (tk[c, jj] / dzp + frac_sno_eff[c] * tk[c, jj - 1] / dzm) -
                    (one(T) - frac_sno_eff[c]) * fact[c, jj] * dhsdT[c]
                bmatrix[c, 5, bm_idx] = -frac_sno_eff[c] * omc * fact[c, jj] *
                    tk[c, jj - 1] / dzm
            elseif j <= nlevgrnd - 1
                dzm = z[c, jj] - z[c, jj - 1]
                dzp = z[c, jj + 1] - z[c, jj]
                bmatrix[c, 2, bm_idx] = -omc * fact[c, jj] * tk[c, jj] / dzp
                bmatrix[c, 3, bm_idx] = one(T) + omc * fact[c, jj] *
                    (tk[c, jj] / dzp + tk[c, jj - 1] / dzm)
                bmatrix[c, 4, bm_idx] = -omc * fact[c, jj] * tk[c, jj - 1] / dzm
            elseif j == nlevgrnd
                dzm = z[c, jj] - z[c, jj - 1]
                bmatrix[c, 3, bm_idx] = one(T) + omc * fact[c, jj] * tk[c, jj - 1] / dzm
                bmatrix[c, 4, bm_idx] = -omc * fact[c, jj] * tk[c, jj - 1] / dzm
            end
        end
    end
end

# h2osfc correction to the soil-layer-1 diagonal (one thread per column). Reads the soil-1
# cell written by the soil kernels above; must run last.
@kernel function _mat_h2osfc_corr_kernel!(bmatrix, @Const(mask), @Const(frac_h2osfc),
        @Const(z), @Const(tk_h2osfc), @Const(dz_h2osfc), @Const(fact), @Const(dhsdT),
        nlevsno::Int)
    c = @index(Global)
    @inbounds if mask[c]
        if frac_h2osfc[c] != 0
            T = eltype(bmatrix)
            omc = one(T) - T(CNFAC)
            joff = nlevsno
            dzm = T(0.5) * dz_h2osfc[c] + z[c, 1 + joff]
            bm_idx = 1 + nlevsno + 1
            bmatrix[c, 3, bm_idx] += frac_h2osfc[c] *
                (omc * fact[c, 1 + joff] * tk_h2osfc[c] / dzm +
                 fact[c, 1 + joff] * dhsdT[c])
            bmatrix[c, 4, bm_idx] = -frac_h2osfc[c] * omc * fact[c, 1 + joff] *
                tk_h2osfc[c] / dzm
        end
    end
end

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

    frac_h2osfc = waterdiagbulk.frac_h2osfc_col
    frac_sno_eff = waterdiagbulk.frac_sno_eff_col

    FT = eltype(bmatrix)
    bmatrix .= zero(FT)
    nc = size(bmatrix, 1)
    dt = convert(FT, dtime)

    # Snow submatrix (one thread per (column, snow level)).
    _launch!(_mat_snow_kernel!, bmatrix, mask_nolakec, col.snl, col.z, tk, fact, dhsdT,
        nlevsno; ndrange = (nc, nlevsno))
    # Standing surface water submatrix.
    _launch!(_mat_ssw_kernel!, bmatrix, mask_nolakec, col.z, tk_h2osfc, dz_h2osfc,
        c_h2osfc, dhsdT, dt, nlevsno; ndrange = nc)
    # Soil submatrix — urban non-road (wall/roof).
    _launch!(_mat_soil_urban_kernel!, bmatrix, mask_nolakec, col.itype, col.snl, col.z,
        col.zi, tk, fact, dhsdT, nlevsno, nlevurb; ndrange = (nc, nlevurb))
    # Soil submatrix — non-urban and urban road.
    _launch!(_mat_soil_kernel!, bmatrix, mask_nolakec, col.itype, col.snl, col.landunit,
        lun.urbpoi, col.z, tk, fact, dhsdT, frac_sno_eff, nlevsno, nlevgrnd;
        ndrange = (nc, nlevgrnd))
    # h2osfc correction to soil layer 1 (reads soil-1 diagonal; must run last).
    _launch!(_mat_h2osfc_corr_kernel!, bmatrix, mask_nolakec, frac_h2osfc, col.z,
        tk_h2osfc, dz_h2osfc, fact, dhsdT, nlevsno; ndrange = nc)

    return nothing
end

# =========================================================================
# phase_change_h2osfc! — Phase change of surface water
# =========================================================================
# Surface-water phase change, one thread per column (no cross-column coupling: every
# write is to column c's own cells, layer 0 = index joff). The init zeroing is folded in
# (each thread zeros its outputs first). Constants (TFRZ/HFUS/CPLIQ/DENICE) and the dtime
# scalar are at the working element type; smooth_min is the same GPU-safe form as the loop.
@kernel function _phase_change_h2osfc_kernel!(t_h2osfc, t_soisno, h2osfc, h2osno_no_layers,
        h2osoi_ice, snow_depth, int_snow, xmf_h2osfc, qflx_h2osfc_to_ice, eflx_h2osfc_to_snow,
        @Const(mask), @Const(snl), @Const(fact), @Const(c_h2osfc), @Const(frac_sno),
        @Const(frac_h2osfc), @Const(h2osoi_liq), @Const(dhsdT), dtime, nlevsno::Int)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(t_soisno)
        joff = nlevsno
        tfrz = T(TFRZ); hfus = T(HFUS); cpliq = T(CPLIQ); denice = T(DENICE)
        xmf_h2osfc[c] = zero(T)
        qflx_h2osfc_to_ice[c] = zero(T)
        eflx_h2osfc_to_snow[c] = zero(T)
        if frac_h2osfc[c] > zero(T) && t_h2osfc[c] <= tfrz
            tinc = tfrz - t_h2osfc[c]
            t_h2osfc[c] = tfrz

            hm = frac_h2osfc[c] * (dhsdT[c] * tinc - tinc * c_h2osfc[c] / dtime)
            xm = hm * dtime / hfus
            temp1 = h2osfc[c] + xm

            h2osno_total = h2osno_no_layers[c]
            if snl[c] < 0
                for j in (snl[c] + 1):0
                    h2osno_total += h2osoi_ice[c, j + joff] + h2osoi_liq[c, j + joff]
                end
            end

            z_avg = frac_sno[c] * snow_depth[c]
            rho_avg = z_avg > zero(T) ? smooth_min(T(800.0), h2osno_total / z_avg) : T(200.0)

            if temp1 >= zero(T)
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

                if frac_sno[c] > zero(T) && snl[c] < 0
                    snow_depth[c] = h2osno_total / (rho_avg * frac_sno[c])
                else
                    snow_depth[c] = h2osno_total / denice
                end

                if snl[c] == 0
                    t_soisno[c, joff] = t_h2osfc[c]
                    eflx_h2osfc_to_snow[c] = zero(T)
                else
                    if snl[c] == -1
                        c1 = frac_sno[c] * (dtime / fact[c, joff] - dhsdT[c] * dtime)
                    else
                        c1 = frac_sno[c] / fact[c, joff] * dtime
                    end
                    c2 = frac_h2osfc[c] != zero(T) ? (-cpliq * xm - frac_h2osfc[c] * dhsdT[c] * dtime) : zero(T)
                    t_soisno[c, joff] = (c1 * t_soisno[c, joff] + c2 * t_h2osfc[c]) / (c1 + c2)
                    eflx_h2osfc_to_snow[c] = (t_h2osfc[c] - t_soisno[c, joff]) * c2 / dtime
                end
            else
                rho_avg = (h2osno_total * rho_avg + h2osfc[c] * denice) / (h2osno_total + h2osfc[c])
                int_snow[c] += h2osfc[c]
                if snl[c] == 0
                    h2osno_no_layers[c] += h2osfc[c]
                else
                    h2osoi_ice[c, joff] += h2osfc[c]
                end
                h2osno_total += h2osfc[c]
                qflx_h2osfc_to_ice[c] = h2osfc[c] / dtime
                t_h2osfc[c] -= temp1 * hfus / (dtime * dhsdT[c] - c_h2osfc[c])
                xmf_h2osfc[c] = hm - frac_h2osfc[c] * temp1 * hfus / dtime

                if snl[c] == 0
                    t_soisno[c, joff] = t_h2osfc[c]
                elseif snl[c] == -1
                    c1 = frac_sno[c] * (dtime / fact[c, joff] - dhsdT[c] * dtime)
                    c2 = frac_h2osfc[c] != zero(T) ? frac_h2osfc[c] * (c_h2osfc[c] - dtime * dhsdT[c]) : zero(T)
                    t_soisno[c, joff] = (c1 * t_soisno[c, joff] + c2 * t_h2osfc[c]) / (c1 + c2)
                    t_h2osfc[c] = t_soisno[c, joff]
                else
                    c1 = frac_sno[c] / fact[c, joff] * dtime
                    c2 = frac_h2osfc[c] != zero(T) ? frac_h2osfc[c] * (c_h2osfc[c] - dtime * dhsdT[c]) : zero(T)
                    t_soisno[c, joff] = (c1 * t_soisno[c, joff] + c2 * t_h2osfc[c]) / (c1 + c2)
                    t_h2osfc[c] = t_soisno[c, joff]
                end

                h2osfc[c] = zero(T)
                if frac_sno[c] > zero(T) && snl[c] < 0
                    snow_depth[c] = h2osno_total / (rho_avg * frac_sno[c])
                else
                    snow_depth[c] = h2osno_total / denice
                end
            end
        end
    end
end

function phase_change_h2osfc!(col::ColumnData, temperature::TemperatureData,
                              energyflux::EnergyFluxData,
                              waterstatebulk::WaterStateBulkData,
                              waterdiagbulk::WaterDiagnosticBulkData,
                              waterfluxbulk::WaterFluxBulkData,
                              mask_nolakec::BitVector, bounds_col::UnitRange{Int},
                              dtime::Real, dhsdT::Vector{<:Real})

    nlevsno = varpar.nlevsno

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
    h2osoi_liq = waterstatebulk.ws.h2osoi_liq_col
    int_snow = waterstatebulk.int_snow_col
    qflx_h2osfc_to_ice = waterfluxbulk.wf.qflx_h2osfc_to_ice_col
    eflx_h2osfc_to_snow = energyflux.eflx_h2osfc_to_snow_col

    dt = convert(eltype(t_soisno), dtime)
    # One thread per column (no cross-column coupling).
    _launch!(_phase_change_h2osfc_kernel!, t_h2osfc, t_soisno, h2osfc, h2osno_no_layers,
        h2osoi_ice, snow_depth, int_snow, xmf_h2osfc, qflx_h2osfc_to_ice, eflx_h2osfc_to_snow,
        mask_nolakec, col.snl, fact, c_h2osfc, frac_sno, frac_h2osfc, h2osoi_liq, dhsdT,
        dt, nlevsno)

    return nothing
end

# =========================================================================
# phase_change_beta! — Phase change within snow and soil layers
# =========================================================================
# Device-view structs grouping the phase-change arrays by role so the kernel takes a few
# struct args instead of ~30 loose ones (same Adapt pattern as the gnet kernel). Int/Bool
# arrays (imelt, masks, snl, landunit/itype indices) stay loose. PcbTmp holds the per-call
# scratch the scalar version allocated locally (hm/xm/.../tinc), passed in so each column
# thread owns its row; allocated with `similar` so it lands on the same backend as the state.
Base.@kwdef struct PcbIn{M,V}    # read-only float inputs (M = per-layer matrix, V = per-column vector)
    dz::M; fact::M; bsw::M; sucsat::M; watsat::M
    dhsdT::V; frac_sno_eff::V; frac_h2osfc::V
end
Base.@kwdef struct PcbLyr{V}     # read-write per-(column, layer) float arrays
    t_soisno::V; h2osoi_ice::V; h2osoi_liq::V; excess_ice::V; exice_subs::V
    qflx_snomelt_lyr::V; qflx_snofrz_lyr::V
end
Base.@kwdef struct PcbCol{V}     # read-write per-column float arrays
    xmf::V; qflx_snomelt::V; qflx_snofrz::V; qflx_snow_drain::V
    h2osno_no_layers::V; snow_depth::V; snomelt_accum::V
    eflx_snomelt::V; eflx_snomelt_r::V; eflx_snomelt_u::V
end
Base.@kwdef struct PcbTmp{V}     # per-call scratch (was local arrays in the scalar version)
    hm::V; xm::V; xm2::V; wice0::V; wliq0::V; wexice0::V; wmass0::V; supercool::V; tinc::V
end
Adapt.@adapt_structure PcbIn
Adapt.@adapt_structure PcbLyr
Adapt.@adapt_structure PcbCol
Adapt.@adapt_structure PcbTmp

# One thread per column: the four original (j-outer, c-inner) loops run as four sequential
# j-loops inside the per-column thread — a loop interchange of independent column iterations,
# byte-identical to the scalar version (no column reads another column). The per-column
# accumulators (xmf, qflx_snomelt/snofrz, snomelt_accum) sum over j in the same ascending
# order, preserving the floating-point result. Constants and the dtime scalar are at the
# working element type. supercool is indexed [c, j]; everything else [c, jj], jj = j + nlevsno.
@kernel function _phase_change_beta_kernel!(lyr::PcbLyr, colv::PcbCol, pin::PcbIn,
        tmp::PcbTmp, imelt, @Const(mask), @Const(urbpoi), @Const(snl), @Const(landunit),
        @Const(itype), @Const(lun_itype), dtime, nlevsno::Int, nlevgrnd::Int,
        nlevurb::Int, nlevmaxurbgrnd::Int)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(lyr.t_soisno)
        joff = nlevsno
        tfrz = T(TFRZ); hfus = T(HFUS); grav = T(GRAV); denice = T(DENICE); thou = T(1000.0)
        l = landunit[c]
        it = itype[c]
        wallroof = (it == ICOL_SUNWALL || it == ICOL_SHADEWALL || it == ICOL_ROOF)

        # Initialization (per-column scalars)
        colv.xmf[c] = zero(T)
        colv.qflx_snomelt[c] = zero(T)
        colv.qflx_snofrz[c] = zero(T)
        colv.qflx_snow_drain[c] = zero(T)

        # Initialize layer variables
        for j in (-nlevsno + 1):nlevmaxurbgrnd
            jj = j + joff
            if j >= snl[c] + 1
                imelt[c, jj] = 0
                tmp.hm[c, jj] = zero(T)
                tmp.xm[c, jj] = zero(T)
                tmp.xm2[c, jj] = zero(T)
                tmp.wice0[c, jj] = lyr.h2osoi_ice[c, jj]
                tmp.wliq0[c, jj] = lyr.h2osoi_liq[c, jj]
                tmp.wexice0[c, jj] = j >= 1 ? lyr.excess_ice[c, j] : zero(T)
                tmp.wmass0[c, jj] = lyr.h2osoi_ice[c, jj] + lyr.h2osoi_liq[c, jj] + tmp.wexice0[c, jj]
                if j >= 1
                    lyr.exice_subs[c, j] = zero(T)
                end
            end
            if j <= 0
                lyr.qflx_snomelt_lyr[c, jj] = zero(T)
                lyr.qflx_snofrz_lyr[c, jj] = zero(T)
            end
        end

        # Snow layers: melting/freezing identification
        for j in (-nlevsno + 1):0
            jj = j + joff
            if j >= snl[c] + 1
                if lyr.h2osoi_ice[c, jj] > zero(T) && lyr.t_soisno[c, jj] > tfrz
                    imelt[c, jj] = 1
                    tmp.tinc[c, jj] = tfrz - lyr.t_soisno[c, jj]
                    lyr.t_soisno[c, jj] = tfrz
                end
                if lyr.h2osoi_liq[c, jj] > zero(T) && lyr.t_soisno[c, jj] < tfrz
                    imelt[c, jj] = 2
                    tmp.tinc[c, jj] = tfrz - lyr.t_soisno[c, jj]
                    lyr.t_soisno[c, jj] = tfrz
                end
            end
        end

        # Soil layers: melting/freezing identification
        for j in 1:nlevmaxurbgrnd
            jj = j + joff
            tmp.supercool[c, j] = zero(T)
            active = (!wallroof && j <= nlevgrnd) || (wallroof && j <= nlevurb)
            if active
                if lyr.h2osoi_ice[c, jj] > zero(T) && lyr.t_soisno[c, jj] > tfrz
                    imelt[c, jj] = 1
                    tmp.tinc[c, jj] = tfrz - lyr.t_soisno[c, jj]
                    lyr.t_soisno[c, jj] = tfrz
                end

                if lyr.excess_ice[c, j] > zero(T) && lyr.t_soisno[c, jj] > tfrz
                    imelt[c, jj] = 1
                    tmp.tinc[c, jj] = tfrz - lyr.t_soisno[c, jj]
                    lyr.t_soisno[c, jj] = tfrz
                end

                tmp.supercool[c, j] = zero(T)
                if lun_itype[l] == ISTSOIL || lun_itype[l] == ISTCROP || it == ICOL_ROAD_PERV
                    if lyr.t_soisno[c, jj] < tfrz
                        smp = hfus * (tfrz - lyr.t_soisno[c, jj]) / (grav * lyr.t_soisno[c, jj]) * thou
                        tmp.supercool[c, j] = pin.watsat[c, j] * (smp / pin.sucsat[c, j])^(-one(T) / pin.bsw[c, j])
                        tmp.supercool[c, j] *= pin.dz[c, jj] * thou
                    end
                end

                if lyr.h2osoi_liq[c, jj] > tmp.supercool[c, j] && lyr.t_soisno[c, jj] < tfrz
                    imelt[c, jj] = 2
                    tmp.tinc[c, jj] = tfrz - lyr.t_soisno[c, jj]
                    lyr.t_soisno[c, jj] = tfrz
                end

                if colv.h2osno_no_layers[c] > zero(T) && j == 1
                    if lyr.t_soisno[c, jj] > tfrz
                        imelt[c, jj] = 1
                        tmp.tinc[c, jj] = tfrz - lyr.t_soisno[c, jj]
                        lyr.t_soisno[c, jj] = tfrz
                    end
                end
            end
        end

        # Energy surplus/loss and rate of phase change
        for j in (-nlevsno + 1):nlevmaxurbgrnd
            jj = j + joff
            active = (!wallroof && j <= nlevgrnd) || (wallroof && j <= nlevurb)
            if active && j >= snl[c] + 1
                if imelt[c, jj] > 0
                    if j == snl[c] + 1
                        if j > 0
                            tmp.hm[c, jj] = pin.dhsdT[c] * tmp.tinc[c, jj] - tmp.tinc[c, jj] / pin.fact[c, jj]
                        else
                            tmp.hm[c, jj] = pin.frac_sno_eff[c] * (pin.dhsdT[c] * tmp.tinc[c, jj] - tmp.tinc[c, jj] / pin.fact[c, jj])
                        end
                        if j == 1 && pin.frac_h2osfc[c] != zero(T)
                            tmp.hm[c, jj] -= pin.frac_h2osfc[c] * pin.dhsdT[c] * tmp.tinc[c, jj]
                        end
                    elseif j == 1
                        tmp.hm[c, jj] = (one(T) - pin.frac_sno_eff[c] - pin.frac_h2osfc[c]) *
                            pin.dhsdT[c] * tmp.tinc[c, jj] - tmp.tinc[c, jj] / pin.fact[c, jj]
                    else
                        if j < 1
                            tmp.hm[c, jj] = -pin.frac_sno_eff[c] * tmp.tinc[c, jj] / pin.fact[c, jj]
                        else
                            tmp.hm[c, jj] = -tmp.tinc[c, jj] / pin.fact[c, jj]
                        end
                    end
                end

                if imelt[c, jj] == 1 && tmp.hm[c, jj] < zero(T)
                    tmp.hm[c, jj] = zero(T)
                    imelt[c, jj] = 0
                end
                if imelt[c, jj] == 2 && tmp.hm[c, jj] > zero(T)
                    tmp.hm[c, jj] = zero(T)
                    imelt[c, jj] = 0
                end

                if imelt[c, jj] > 0 && abs(tmp.hm[c, jj]) > zero(T)
                    tmp.xm[c, jj] = tmp.hm[c, jj] * dtime / hfus

                    # Unresolved snow melting
                    if j == 1
                        if colv.h2osno_no_layers[c] > zero(T) && tmp.xm[c, jj] > zero(T)
                            temp1 = colv.h2osno_no_layers[c]
                            colv.h2osno_no_layers[c] = smooth_max(zero(T), temp1 - tmp.xm[c, jj])
                            propor = colv.h2osno_no_layers[c] / temp1
                            colv.snow_depth[c] = propor * colv.snow_depth[c]
                            heatr = tmp.hm[c, jj] - hfus * (temp1 - colv.h2osno_no_layers[c]) / dtime
                            if heatr > zero(T)
                                tmp.xm[c, jj] = heatr * dtime / hfus
                                tmp.hm[c, jj] = heatr
                            else
                                tmp.xm[c, jj] = zero(T)
                                tmp.hm[c, jj] = zero(T)
                            end
                            colv.qflx_snomelt[c] = smooth_max(zero(T), temp1 - colv.h2osno_no_layers[c]) / dtime
                            colv.xmf[c] = hfus * colv.qflx_snomelt[c]
                            colv.qflx_snow_drain[c] = colv.qflx_snomelt[c]
                        end
                    end

                    heatr = zero(T)
                    if tmp.xm[c, jj] > zero(T)
                        lyr.h2osoi_ice[c, jj] = smooth_max(zero(T), tmp.wice0[c, jj] - tmp.xm[c, jj])
                        heatr = tmp.hm[c, jj] - hfus * (tmp.wice0[c, jj] - lyr.h2osoi_ice[c, jj]) / dtime
                        tmp.xm2[c, jj] = tmp.xm[c, jj] - lyr.h2osoi_ice[c, jj]
                        if lyr.h2osoi_ice[c, jj] == zero(T)
                            if tmp.wexice0[c, jj] >= zero(T) && tmp.xm2[c, jj] > zero(T) && j >= 2
                                lyr.excess_ice[c, j] = smooth_max(zero(T), tmp.wexice0[c, jj] - tmp.xm2[c, jj])
                                heatr = tmp.hm[c, jj] - hfus * (tmp.wexice0[c, jj] - lyr.excess_ice[c, j] +
                                        tmp.wice0[c, jj] - lyr.h2osoi_ice[c, jj]) / dtime
                            end
                        end
                    elseif tmp.xm[c, jj] < zero(T)
                        if j <= 0
                            lyr.h2osoi_ice[c, jj] = smooth_min(tmp.wmass0[c, jj], tmp.wice0[c, jj] - tmp.xm[c, jj])
                        else
                            if tmp.wmass0[c, jj] - tmp.wexice0[c, jj] < tmp.supercool[c, j]
                                lyr.h2osoi_ice[c, jj] = zero(T)
                            else
                                lyr.h2osoi_ice[c, jj] = smooth_min(tmp.wmass0[c, jj] - tmp.wexice0[c, jj] - tmp.supercool[c, j],
                                    tmp.wice0[c, jj] - tmp.xm[c, jj])
                            end
                        end
                        heatr = tmp.hm[c, jj] - hfus * (tmp.wice0[c, jj] - lyr.h2osoi_ice[c, jj]) / dtime
                    end

                    ei_val = j >= 1 ? lyr.excess_ice[c, j] : zero(T)
                    lyr.h2osoi_liq[c, jj] = smooth_max(zero(T), tmp.wmass0[c, jj] - lyr.h2osoi_ice[c, jj] - ei_val)

                    if abs(heatr) > zero(T)
                        if j == snl[c] + 1
                            if j == 1
                                lyr.t_soisno[c, jj] += pin.fact[c, jj] * heatr /
                                    (one(T) - (one(T) - pin.frac_h2osfc[c]) * pin.fact[c, jj] * pin.dhsdT[c])
                            else
                                if pin.frac_sno_eff[c] > zero(T)
                                    lyr.t_soisno[c, jj] += (pin.fact[c, jj] / pin.frac_sno_eff[c]) * heatr /
                                        (one(T) - pin.fact[c, jj] * pin.dhsdT[c])
                                end
                            end
                        elseif j == 1
                            lyr.t_soisno[c, jj] += pin.fact[c, jj] * heatr /
                                (one(T) - (one(T) - pin.frac_sno_eff[c] - pin.frac_h2osfc[c]) * pin.fact[c, jj] * pin.dhsdT[c])
                        else
                            if j > 0
                                lyr.t_soisno[c, jj] += pin.fact[c, jj] * heatr
                            else
                                if pin.frac_sno_eff[c] > zero(T)
                                    lyr.t_soisno[c, jj] += (pin.fact[c, jj] / pin.frac_sno_eff[c]) * heatr
                                end
                            end
                        end

                        if j <= 0
                            if lyr.h2osoi_liq[c, jj] * lyr.h2osoi_ice[c, jj] > zero(T)
                                lyr.t_soisno[c, jj] = tfrz
                            end
                        end
                    end

                    if j >= 1
                        colv.xmf[c] += hfus * (tmp.wice0[c, jj] - lyr.h2osoi_ice[c, jj]) / dtime +
                            hfus * (tmp.wexice0[c, jj] - (j >= 1 ? lyr.excess_ice[c, j] : zero(T))) / dtime
                        lyr.exice_subs[c, j] = smooth_max(zero(T), (tmp.wexice0[c, jj] - lyr.excess_ice[c, j]) / denice)
                    else
                        colv.xmf[c] += hfus * (tmp.wice0[c, jj] - lyr.h2osoi_ice[c, jj]) / dtime
                    end

                    if imelt[c, jj] == 1 && j < 1
                        lyr.qflx_snomelt_lyr[c, jj] = smooth_max(zero(T), tmp.wice0[c, jj] - lyr.h2osoi_ice[c, jj]) / dtime
                        colv.qflx_snomelt[c] += lyr.qflx_snomelt_lyr[c, jj]
                        colv.snomelt_accum[c] += lyr.qflx_snomelt_lyr[c, jj] * dtime * T(1.0e-3)
                    end
                    if imelt[c, jj] == 2 && j < 1
                        lyr.qflx_snofrz_lyr[c, jj] = smooth_max(zero(T), lyr.h2osoi_ice[c, jj] - tmp.wice0[c, jj]) / dtime
                        colv.qflx_snofrz[c] += lyr.qflx_snofrz_lyr[c, jj]
                    end
                end
            end
        end

        # History output
        colv.eflx_snomelt[c] = colv.qflx_snomelt[c] * hfus
        if urbpoi[l]
            colv.eflx_snomelt_u[c] = colv.eflx_snomelt[c]
        elseif lun_itype[l] == ISTSOIL || lun_itype[l] == ISTCROP
            colv.eflx_snomelt_r[c] = colv.eflx_snomelt[c]
        end
    end
end

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

    t_soisno = temperature.t_soisno_col
    nc = length(bounds_col)
    FT = eltype(t_soisno)
    nc == 0 && return nothing

    lyr = PcbLyr(; t_soisno = t_soisno,
        h2osoi_ice = waterstatebulk.ws.h2osoi_ice_col,
        h2osoi_liq = waterstatebulk.ws.h2osoi_liq_col,
        excess_ice = waterstatebulk.ws.excess_ice_col,
        exice_subs = waterdiagbulk.exice_subs_col,
        qflx_snomelt_lyr = waterfluxbulk.wf.qflx_snomelt_lyr_col,
        qflx_snofrz_lyr = waterfluxbulk.wf.qflx_snofrz_lyr_col)
    colv = PcbCol(; xmf = temperature.xmf_col,
        qflx_snomelt = waterfluxbulk.wf.qflx_snomelt_col,
        qflx_snofrz = waterfluxbulk.wf.qflx_snofrz_col,
        qflx_snow_drain = waterfluxbulk.wf.qflx_snow_drain_col,
        h2osno_no_layers = waterstatebulk.ws.h2osno_no_layers_col,
        snow_depth = waterdiagbulk.snow_depth_col,
        snomelt_accum = waterdiagbulk.snomelt_accum_col,
        eflx_snomelt = energyflux.eflx_snomelt_col,
        eflx_snomelt_r = energyflux.eflx_snomelt_r_col,
        eflx_snomelt_u = energyflux.eflx_snomelt_u_col)
    pin = PcbIn(; dz = col.dz, fact = temperature.fact_col, dhsdT = dhsdT,
        frac_sno_eff = waterdiagbulk.frac_sno_eff_col,
        frac_h2osfc = waterdiagbulk.frac_h2osfc_col,
        bsw = soilstate.bsw_col, sucsat = soilstate.sucsat_col, watsat = soilstate.watsat_col)

    # Per-call scratch on the same backend as the state (was local zeros() in the loop).
    nlev = nlevsno + nlevmaxurbgrnd
    mk(d2) = fill!(similar(t_soisno, FT, nc, d2), zero(FT))
    tmp = PcbTmp(; hm = mk(nlev), xm = mk(nlev), xm2 = mk(nlev), wice0 = mk(nlev),
        wliq0 = mk(nlev), wexice0 = mk(nlev), wmass0 = mk(nlev),
        supercool = mk(nlevmaxurbgrnd), tinc = mk(nlev))

    dt = convert(FT, dtime)
    backend = _kernel_backend(t_soisno)
    _phase_change_beta_kernel!(backend)(lyr, colv, pin, tmp, temperature.imelt_col,
        mask_nolakec, lun.urbpoi, col.snl, col.landunit, col.itype, lun.itype,
        dt, nlevsno, nlevgrnd, nlevurb, nlevmaxurbgrnd; ndrange = nc)
    KernelAbstractions.synchronize(backend)

    return nothing
end

# =========================================================================
# building_hac! — Building Heating and Cooling (simple method)
# =========================================================================
# Per-landunit building HAC. cool_on/heat_on are BitVector scratch, so this kernel runs on
# the CPU backend (_kernel_backend(::BitArray) = KA.CPU()) — matching the documented
# BitVector-writing-kernel convention; moving the build-temp landunit scratch to device is
# the same follow-up as the other masks. One thread per landunit.
@kernel function _building_hac_kernel!(cool_on, heat_on, t_building, @Const(mask),
        @Const(urbpoi), @Const(t_building_max), @Const(t_building_min))
    l = @index(Global)
    @inbounds if mask[l]
        if urbpoi[l]
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
end

function building_hac!(lun::LandunitData, temperature::TemperatureData,
                       urbanparams::UrbanParamsData,
                       t_building_max::Vector{<:Real},
                       mask_urbanl::BitVector, bounds_lun::UnitRange{Int},
                       cool_on::BitVector, heat_on::BitVector)

    t_building = temperature.t_building_lun
    t_building_min = urbanparams.t_building_min

    # One thread per landunit; keyed on cool_on (BitVector) → CPU backend.
    _launch!(_building_hac_kernel!, cool_on, heat_on, t_building, mask_urbanl,
        lun.urbpoi, t_building_max, t_building_min)

    return nothing
end
