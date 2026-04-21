# ==========================================================================
# Ported from: src/biogeophys/SurfaceAlbedoMod.F90
# Performs surface albedo calculations
#
# Public functions:
#   surface_albedo_init_time_const!  — initialize time-constant albedo data
#   soil_albedo!                     — ground surface albedo
#   two_stream!                      — two-stream canopy radiative transfer
#   surface_albedo!                  — top-level surface albedo driver
# ==========================================================================

# ---- Module-level constants and data ----

# albedo land ice by waveband (1=vis, 2=nir)
const ALBICE = [0.80, 0.55]

# albedo frozen lakes by waveband (1=vis, 2=nir)
const ALBLAK = [0.60, 0.40]

# Coefficient for calculating ice "fraction" for lake surface albedo
# From D. Mironov (2010) Boreal Env. Research
const CALB = 95.6

# prevents overflow for division by zero
const MPE_ALBEDO = 1.0e-06

"""
    SurfaceAlbedoConstants

Module-level data that is initialized once and remains constant during
the simulation. Corresponds to Fortran module variables `albsat`, `albdry`,
`isoicol`, `alblakwi`, and `snowveg_affects_radiation`.
"""
Base.@kwdef mutable struct SurfaceAlbedoConstants{FT<:Real}
    albsat::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # wet soil albedo by color class and waveband
    albdry::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # dry soil albedo by color class and waveband
    isoicol::Vector{Int} = Int[]                              # column soil color class
    alblakwi::Vector{Float64} = [0.10, 0.10]                 # albedo of melting lakes (namelist-settable)
    lake_melt_icealb::Vector{Float64} = [0.10, 0.10]         # namelist default for alblakwi
    snowveg_affects_radiation::Bool = true                    # whether canopy snow affects radiation
end

const surfalb_con = SurfaceAlbedoConstants()

# --------------------------------------------------------------------------
# surface_albedo_init_time_const!
# --------------------------------------------------------------------------

"""
    surface_albedo_init_time_const!(con, mxsoil_color, soic2d,
                                    col_gridcell, bounds_col, bounds_grc)

Initialize time-constant albedo data: soil color indices, saturated/dry soil
albedos, and melting lake albedos.

Ported from `SurfaceAlbedoInitTimeConst` in `SurfaceAlbedoMod.F90`.

# Arguments
- `con::SurfaceAlbedoConstants` : module constants struct to populate
- `mxsoil_color::Int`          : maximum number of soil color classes (8 or 20)
- `soic2d::Vector{Int}`        : soil color class per gridcell
- `col_gridcell::Vector{Int}`  : gridcell index for each column
- `bounds_col::UnitRange{Int}` : column bounds
- `bounds_grc::UnitRange{Int}` : gridcell bounds (unused but kept for API symmetry)
"""
function surface_albedo_init_time_const!(con::SurfaceAlbedoConstants,
                                          mxsoil_color::Int,
                                          soic2d::Vector{Int},
                                          col_gridcell::Vector{Int},
                                          bounds_col::UnitRange{Int},
                                          bounds_grc::UnitRange{Int})
    # Map gridcell soil color to column
    con.isoicol = zeros(Int, length(bounds_col))
    for c in bounds_col
        g = col_gridcell[c]
        con.isoicol[c] = soic2d[g]
    end

    # Allocate and fill saturated/dry soil albedos
    con.albsat = zeros(mxsoil_color, NUMRAD)
    con.albdry = zeros(mxsoil_color, NUMRAD)

    if mxsoil_color == 8
        con.albsat[1:8, 1] = [0.12, 0.11, 0.10, 0.09, 0.08, 0.07, 0.06, 0.05]
        con.albsat[1:8, 2] = [0.24, 0.22, 0.20, 0.18, 0.16, 0.14, 0.12, 0.10]
        con.albdry[1:8, 1] = [0.24, 0.22, 0.20, 0.18, 0.16, 0.14, 0.12, 0.10]
        con.albdry[1:8, 2] = [0.48, 0.44, 0.40, 0.36, 0.32, 0.28, 0.24, 0.20]
    elseif mxsoil_color == 20
        con.albsat[1:20, 1] = [0.25, 0.23, 0.21, 0.20, 0.19, 0.18, 0.17, 0.16,
                               0.15, 0.14, 0.13, 0.12, 0.11, 0.10, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04]
        con.albsat[1:20, 2] = [0.50, 0.46, 0.42, 0.40, 0.38, 0.36, 0.34, 0.32,
                               0.30, 0.28, 0.26, 0.24, 0.22, 0.20, 0.18, 0.16, 0.14, 0.12, 0.10, 0.08]
        con.albdry[1:20, 1] = [0.36, 0.34, 0.32, 0.31, 0.30, 0.29, 0.28, 0.27,
                               0.26, 0.25, 0.24, 0.23, 0.22, 0.20, 0.18, 0.16, 0.14, 0.12, 0.10, 0.08]
        con.albdry[1:20, 2] = [0.61, 0.57, 0.53, 0.51, 0.49, 0.48, 0.45, 0.43,
                               0.41, 0.39, 0.37, 0.35, 0.33, 0.31, 0.29, 0.27, 0.25, 0.23, 0.21, 0.16]
    else
        error("maximum color class = $mxsoil_color is not supported")
    end

    # Set alblakwi from namelist value
    con.alblakwi .= con.lake_melt_icealb

    return nothing
end

# --------------------------------------------------------------------------
# soil_albedo!
# --------------------------------------------------------------------------

"""
    soil_albedo!(surfalb, con, col, lun,
                 coszen_col, temperature, waterstatebulk, lakestate,
                 mask_nourbanc, bounds_col;
                 lakepuddling=false)

Determine ground surface albedo, accounting for snow/lake/glacier/wetland.

Ported from subroutine `SoilAlbedo` in `SurfaceAlbedoMod.F90`.
"""
function soil_albedo!(surfalb::SurfaceAlbedoData,
                      con::SurfaceAlbedoConstants,
                      col::ColumnData,
                      lun::LandunitData,
                      coszen_col::Vector{<:Real},
                      temperature::TemperatureData,
                      waterstatebulk::WaterStateBulkData,
                      lakestate::LakeStateData,
                      mask_nourbanc::BitVector,
                      bounds_col::UnitRange{Int};
                      lakepuddling::Bool = false)

    nband = NUMRAD

    for ib in 1:nband
        for c in bounds_col
            mask_nourbanc[c] || continue
            if coszen_col[c] > 0.0
                l = col.landunit[c]

                if lun.itype[l] == ISTSOIL || lun.itype[l] == ISTCROP
                    # Soil
                    inc = smooth_max(0.11 - 0.40 * waterstatebulk.ws.h2osoi_vol_col[c, 1], 0.0)
                    soilcol = con.isoicol[c]
                    surfalb.albsod_col[c, ib] = min(con.albsat[soilcol, ib] + inc,
                                                    con.albdry[soilcol, ib])
                    surfalb.albsoi_col[c, ib] = surfalb.albsod_col[c, ib]

                elseif lun.itype[l] == ISTICE
                    # Land ice
                    surfalb.albsod_col[c, ib] = ALBICE[ib]
                    surfalb.albsoi_col[c, ib] = surfalb.albsod_col[c, ib]

                elseif temperature.t_grnd_col[c] > TFRZ ||
                       (lakepuddling && lun.itype[l] == ISTDLAK &&
                        temperature.t_grnd_col[c] == TFRZ &&
                        lakestate.lake_icefrac_col[c, 1] < 1.0 &&
                        lakestate.lake_icefrac_col[c, 2] > 0.0)
                    # Unfrozen lake/wetland
                    albsod_unfrozen = 0.05 / (smooth_max(0.001, coszen_col[c]) + 0.15)
                    if lun.itype[l] == ISTDLAK
                        albsoi_unfrozen = 0.10
                    else
                        albsoi_unfrozen = albsod_unfrozen
                    end

                    # Blend with frozen albedo via smooth_heaviside for AD
                    if lun.itype[l] == ISTDLAK && !lakepuddling && col.snl[c] == 0
                        sicefr = 1.0 - exp(-CALB * (TFRZ - temperature.t_grnd_col[c]) / TFRZ)
                        albsod_frozen = sicefr * ALBLAK[ib] +
                            (1.0 - sicefr) * smooth_max(con.alblakwi[ib],
                                                  0.05 / (smooth_max(0.001, coszen_col[c]) + 0.15))
                        albsoi_frozen = sicefr * ALBLAK[ib] +
                            (1.0 - sicefr) * smooth_max(con.alblakwi[ib], 0.10)
                    else
                        albsod_frozen = convert(eltype(temperature.t_grnd_col), ALBLAK[ib])
                        albsoi_frozen = albsod_frozen
                    end

                    w_unfrozen = smooth_heaviside(temperature.t_grnd_col[c] - TFRZ)
                    surfalb.albsod_col[c, ib] = w_unfrozen * albsod_unfrozen + (1.0 - w_unfrozen) * albsod_frozen
                    surfalb.albsoi_col[c, ib] = w_unfrozen * albsoi_unfrozen + (1.0 - w_unfrozen) * albsoi_frozen

                else
                    # Frozen lake/wetland
                    if lun.itype[l] == ISTDLAK && !lakepuddling && col.snl[c] == 0
                        sicefr = 1.0 - exp(-CALB * (TFRZ - temperature.t_grnd_col[c]) / TFRZ)
                        surfalb.albsod_col[c, ib] = sicefr * ALBLAK[ib] +
                            (1.0 - sicefr) * smooth_max(con.alblakwi[ib],
                                                  0.05 / (smooth_max(0.001, coszen_col[c]) + 0.15))
                        surfalb.albsoi_col[c, ib] = sicefr * ALBLAK[ib] +
                            (1.0 - sicefr) * smooth_max(con.alblakwi[ib], 0.10)
                    else
                        surfalb.albsod_col[c, ib] = ALBLAK[ib]
                        surfalb.albsoi_col[c, ib] = surfalb.albsod_col[c, ib]
                    end
                end
            end
        end
    end

    return nothing
end

# --------------------------------------------------------------------------
# two_stream!
# --------------------------------------------------------------------------

"""
    two_stream!(surfalb, patchdata, col, canopystate, temperature,
                waterdiagbulk, coszen_patch, rho, tau,
                mask_vegsol, bounds_patch;
                SFonly=false)

Two-stream fluxes for canopy radiative transfer.

Uses two-stream approximation of Dickinson (1983) and Sellers (1985) to
calculate fluxes absorbed, reflected, and transmitted by vegetation.
Calculates sunlit and shaded fluxes as described by Bonan et al (2011).

Ported from subroutine `TwoStream` in `SurfaceAlbedoMod.F90`.

# Arguments
- `mask_vegsol::BitVector` : mask for vegetated patches with coszen > 0
- `SFonly::Bool`           : if true, only calculate snow-free albedos
"""
function two_stream!(surfalb::SurfaceAlbedoData,
                     patchdata::PatchData,
                     col::ColumnData,
                     canopystate::CanopyStateData,
                     temperature::TemperatureData,
                     waterdiagbulk::WaterDiagnosticBulkData,
                     coszen_patch::Vector{<:Real},
                     rho_in::Matrix{<:Real},
                     tau_in::Matrix{<:Real},
                     pftcon_xl::Vector{<:Real},
                     mask_vegsol::BitVector,
                     bounds_patch::UnitRange{Int};
                     SFonly::Bool = false)

    numrad = NUMRAD
    nlevcan = NLEVCAN
    lSFonly = SFonly

    # Aliases
    elai = canopystate.elai_patch
    esai = canopystate.esai_patch
    t_veg = temperature.t_veg_patch
    fwet = waterdiagbulk.fwet_patch
    fcansno = waterdiagbulk.fcansno_patch
    tlai_z = surfalb.tlai_z_patch
    tsai_z = surfalb.tsai_z_patch
    nrad = surfalb.nrad_patch
    albgrd = surfalb.albgrd_col
    albgri = surfalb.albgri_col
    albsod = surfalb.albsod_col
    albsoi = surfalb.albsoi_col

    # Pre-allocate per-patch arrays
    np = length(bounds_patch)
    FT = eltype(coszen_patch)
    chil = zeros(FT, np)
    gdir_arr = zeros(FT, np)
    twostext = zeros(FT, np)
    avmu = zeros(FT, np)
    temp0_arr = zeros(FT, np)
    temp2_arr = zeros(FT, np)
    omega = zeros(FT, np, numrad)

    # Calculate two-stream parameters independent of waveband
    for p in bounds_patch
        mask_vegsol[p] || continue

        cosz = smooth_max(0.001, coszen_patch[p])

        chil[p] = smooth_clamp(pftcon_xl[patchdata.itype[p] + 1], -0.4, 0.6)
        if abs(chil[p]) <= 0.01
            chil[p] = 0.01
        end
        phi1 = 0.5 - 0.633 * chil[p] - 0.330 * chil[p] * chil[p]
        phi2 = 0.877 * (1.0 - 2.0 * phi1)
        gdir_arr[p] = phi1 + phi2 * cosz
        twostext[p] = gdir_arr[p] / cosz
        avmu[p] = (1.0 - phi1 / phi2 * log((phi1 + phi2) / phi1)) / phi2
        temp0_arr[p] = smooth_max(gdir_arr[p] + phi2 * cosz, 1.0e-6)
        temp1_val = phi1 * cosz
        temp2_arr[p] = (1.0 - temp1_val / temp0_arr[p] * log((temp1_val + temp0_arr[p]) / temp1_val))
    end

    # Loop over wavebands
    for ib in 1:numrad
        for p in bounds_patch
            mask_vegsol[p] || continue
            c = patchdata.column[p]

            cosz = smooth_max(0.001, coszen_patch[p])

            # Two-stream parameters omega, betad, betai
            omegal = rho_in[p, ib] + tau_in[p, ib]
            asu = 0.5 * omegal * gdir_arr[p] / temp0_arr[p] * temp2_arr[p]
            betadl = (1.0 + avmu[p] * twostext[p]) / (omegal * avmu[p] * twostext[p]) * asu
            betail = 0.5 * ((rho_in[p, ib] + tau_in[p, ib]) + (rho_in[p, ib] - tau_in[p, ib]) *
                     ((1.0 + chil[p]) / 2.0)^2) / omegal

            # Adjust omega, betad, betai for intercepted snow
            if lSFonly || (!surfalb_con.snowveg_affects_radiation && t_veg[p] > TFRZ)
                tmp0 = omegal
                tmp1 = betadl
                tmp2 = betail
            else
                if surfalb_con.snowveg_affects_radiation
                    tmp0 = (1.0 - fcansno[p]) * omegal + fcansno[p] * OMEGAS[ib]
                    tmp1 = ((1.0 - fcansno[p]) * omegal * betadl + fcansno[p] * OMEGAS[ib] * BETADS) / tmp0
                    tmp2 = ((1.0 - fcansno[p]) * omegal * betail + fcansno[p] * OMEGAS[ib] * BETAIS) / tmp0
                else
                    tmp0 = (1.0 - fwet[p]) * omegal + fwet[p] * OMEGAS[ib]
                    tmp1 = ((1.0 - fwet[p]) * omegal * betadl + fwet[p] * OMEGAS[ib] * BETADS) / tmp0
                    tmp2 = ((1.0 - fwet[p]) * omegal * betail + fwet[p] * OMEGAS[ib] * BETAIS) / tmp0
                end
            end

            omega[p, ib] = tmp0
            betad = tmp1
            betai = tmp2

            # Common terms
            b = 1.0 - omega[p, ib] + omega[p, ib] * betai
            c1 = omega[p, ib] * betai
            tmp0_val = avmu[p] * twostext[p]
            d = tmp0_val * omega[p, ib] * betad
            f = tmp0_val * omega[p, ib] * (1.0 - betad)
            tmp1_val = b * b - c1 * c1
            h = sqrt(tmp1_val) / avmu[p]
            sigma = tmp0_val * tmp0_val - tmp1_val

            p1 = b + avmu[p] * h
            p2 = b - avmu[p] * h
            p3 = b + tmp0_val
            p4 = b - tmp0_val

            # Absorbed, reflected, transmitted fluxes for full canopy
            t1 = smooth_min(h * (elai[p] + esai[p]), 40.0)
            s1 = exp(-t1)
            t1 = smooth_min(twostext[p] * (elai[p] + esai[p]), 40.0)
            s2 = exp(-t1)

            # ---- Direct beam ----
            if !lSFonly
                u1 = b - c1 / albgrd[c, ib]
                u2 = b - c1 * albgrd[c, ib]
                u3 = f + c1 * albgrd[c, ib]
            else
                u1 = b - c1 / albsod[c, ib]
                u2 = b - c1 * albsod[c, ib]
                u3 = f + c1 * albsod[c, ib]
            end
            tmp2_val = u1 - avmu[p] * h
            tmp3 = u1 + avmu[p] * h
            d1 = p1 * tmp2_val / s1 - p2 * tmp3 * s1
            tmp4 = u2 + avmu[p] * h
            tmp5 = u2 - avmu[p] * h
            d2 = tmp4 / s1 - tmp5 * s1
            h1 = -d * p4 - c1 * f
            tmp6 = d - h1 * p3 / sigma
            tmp7 = (d - c1 - h1 / sigma * (u1 + tmp0_val)) * s2
            h2 = (tmp6 * tmp2_val / s1 - p2 * tmp7) / d1
            h3 = -(tmp6 * tmp3 * s1 - p1 * tmp7) / d1
            h4 = -f * p3 - c1 * d
            tmp8 = h4 / sigma
            tmp9 = (u3 - tmp8 * (u2 - tmp0_val)) * s2
            h5 = -(tmp8 * tmp4 / s1 + tmp9) / d2
            h6 = (tmp8 * tmp5 * s1 + tmp9) / d2

            if !lSFonly
                surfalb.albd_patch[p, ib] = h1 / sigma + h2 + h3
                surfalb.ftid_patch[p, ib] = h4 * s2 / sigma + h5 * s1 + h6 / s1
                surfalb.ftdd_patch[p, ib] = s2
                surfalb.fabd_patch[p, ib] = 1.0 - surfalb.albd_patch[p, ib] -
                    (1.0 - albgrd[c, ib]) * surfalb.ftdd_patch[p, ib] -
                    (1.0 - albgri[c, ib]) * surfalb.ftid_patch[p, ib]
            else
                surfalb.albdSF_patch[p, ib] = h1 / sigma + h2 + h3
            end

            a1 = h1 / sigma * (1.0 - s2 * s2) / (2.0 * twostext[p]) +
                 h2 * (1.0 - s2 * s1) / (twostext[p] + h) +
                 h3 * (1.0 - s2 / s1) / (twostext[p] - h)

            a2 = h4 / sigma * (1.0 - s2 * s2) / (2.0 * twostext[p]) +
                 h5 * (1.0 - s2 * s1) / (twostext[p] + h) +
                 h6 * (1.0 - s2 / s1) / (twostext[p] - h)

            if !lSFonly
                surfalb.fabd_sun_patch[p, ib] = (1.0 - omega[p, ib]) *
                    (1.0 - s2 + 1.0 / avmu[p] * (a1 + a2))
                surfalb.fabd_sha_patch[p, ib] = surfalb.fabd_patch[p, ib] -
                    surfalb.fabd_sun_patch[p, ib]
            end

            # ---- Diffuse ----
            if !lSFonly
                u1 = b - c1 / albgri[c, ib]
                u2 = b - c1 * albgri[c, ib]
            else
                u1 = b - c1 / albsoi[c, ib]
                u2 = b - c1 * albsoi[c, ib]
            end
            tmp2_val = u1 - avmu[p] * h
            tmp3 = u1 + avmu[p] * h
            d1 = p1 * tmp2_val / s1 - p2 * tmp3 * s1
            tmp4 = u2 + avmu[p] * h
            tmp5 = u2 - avmu[p] * h
            d2 = tmp4 / s1 - tmp5 * s1
            h7 = (c1 * tmp2_val) / (d1 * s1)
            h8 = (-c1 * tmp3 * s1) / d1
            h9 = tmp4 / (d2 * s1)
            h10 = (-tmp5 * s1) / d2

            if lSFonly
                surfalb.albiSF_patch[p, ib] = h7 + h8
            else
                surfalb.albi_patch[p, ib] = h7 + h8
                surfalb.ftii_patch[p, ib] = h9 * s1 + h10 / s1
                surfalb.fabi_patch[p, ib] = 1.0 - surfalb.albi_patch[p, ib] -
                    (1.0 - albgri[c, ib]) * surfalb.ftii_patch[p, ib]

                a1 = h7 * (1.0 - s2 * s1) / (twostext[p] + h) +
                     h8 * (1.0 - s2 / s1) / (twostext[p] - h)
                a2 = h9 * (1.0 - s2 * s1) / (twostext[p] + h) +
                     h10 * (1.0 - s2 / s1) / (twostext[p] - h)

                surfalb.fabi_sun_patch[p, ib] = (1.0 - omega[p, ib]) / avmu[p] * (a1 + a2)
                surfalb.fabi_sha_patch[p, ib] = surfalb.fabi_patch[p, ib] -
                    surfalb.fabi_sun_patch[p, ib]

                # Per-layer derivatives for PAR (ib == 1)
                if ib == 1
                    if nlevcan == 1
                        # Sun/shade big leaf: single layer
                        fsun_z = surfalb.fsun_z_patch
                        fsun_z[p, 1] = (1.0 - s2) / t1

                        laisum = elai[p] + esai[p]
                        surfalb.fabd_sun_z_patch[p, 1] = surfalb.fabd_sun_patch[p, ib] /
                            (fsun_z[p, 1] * laisum)
                        surfalb.fabi_sun_z_patch[p, 1] = surfalb.fabi_sun_patch[p, ib] /
                            (fsun_z[p, 1] * laisum)
                        surfalb.fabd_sha_z_patch[p, 1] = surfalb.fabd_sha_patch[p, ib] /
                            ((1.0 - fsun_z[p, 1]) * laisum)
                        surfalb.fabi_sha_z_patch[p, 1] = surfalb.fabi_sha_patch[p, ib] /
                            ((1.0 - fsun_z[p, 1]) * laisum)

                        # Leaf to canopy scaling coefficients
                        extkn = 0.30
                        extkb = twostext[p]
                        surfalb.vcmaxcintsun_patch[p] = (1.0 - exp(-(extkn + extkb) * elai[p])) /
                            (extkn + extkb)
                        surfalb.vcmaxcintsha_patch[p] = (1.0 - exp(-extkn * elai[p])) / extkn -
                            surfalb.vcmaxcintsun_patch[p]
                        if elai[p] > 0.0
                            surfalb.vcmaxcintsun_patch[p] /= (fsun_z[p, 1] * elai[p])
                            surfalb.vcmaxcintsha_patch[p] /= ((1.0 - fsun_z[p, 1]) * elai[p])
                        else
                            surfalb.vcmaxcintsun_patch[p] = 0.0
                            surfalb.vcmaxcintsha_patch[p] = 0.0
                        end

                    elseif nlevcan > 1
                        # Multi-layer canopy
                        for iv in 1:nrad[p]
                            # Cumulative lai+sai at center of layer
                            if iv == 1
                                laisum_l = 0.5 * (tlai_z[p, iv] + tsai_z[p, iv])
                            else
                                laisum_l = laisum_l + 0.5 * ((tlai_z[p, iv-1] + tsai_z[p, iv-1]) +
                                    (tlai_z[p, iv] + tsai_z[p, iv]))
                            end

                            t1_l = smooth_min(h * laisum_l, 40.0)
                            s1_l = exp(-t1_l)
                            t1_l = smooth_min(twostext[p] * laisum_l, 40.0)
                            s2_l = exp(-t1_l)
                            surfalb.fsun_z_patch[p, iv] = s2_l

                            # ---- Direct beam derivatives ----
                            u1_l = b - c1 / albgrd[c, ib]
                            u2_l = b - c1 * albgrd[c, ib]
                            u3_l = f + c1 * albgrd[c, ib]

                            tmp2_l = u1_l - avmu[p] * h
                            tmp3_l = u1_l + avmu[p] * h
                            d1_l = p1 * tmp2_l / s1_l - p2 * tmp3_l * s1_l
                            tmp4_l = u2_l + avmu[p] * h
                            tmp5_l = u2_l - avmu[p] * h
                            d2_l = tmp4_l / s1_l - tmp5_l * s1_l

                            h1_l = -d * p4 - c1 * f
                            tmp6_l = d - h1_l * p3 / sigma
                            tmp7_l = (d - c1 - h1_l / sigma * (u1_l + tmp0_val)) * s2_l

                            h2_l = (tmp6_l * tmp2_l / s1_l - p2 * tmp7_l) / d1_l
                            h3_l = -(tmp6_l * tmp3_l * s1_l - p1 * tmp7_l) / d1_l

                            h4_l = -f * p3 - c1 * d
                            tmp8_l = h4_l / sigma
                            tmp9_l = (u3_l - tmp8_l * (u2_l - tmp0_val)) * s2_l
                            h5_l = -(tmp8_l * tmp4_l / s1_l + tmp9_l) / d2_l
                            h6_l = (tmp8_l * tmp5_l * s1_l + tmp9_l) / d2_l

                            # Derivatives
                            v_val = d1_l
                            dv = h * p1 * tmp2_l / s1_l + h * p2 * tmp3_l * s1_l

                            u_val = tmp6_l * tmp2_l / s1_l - p2 * tmp7_l
                            du = h * tmp6_l * tmp2_l / s1_l + twostext[p] * p2 * tmp7_l
                            dh2 = (v_val * du - u_val * dv) / (v_val * v_val)

                            u_val = -tmp6_l * tmp3_l * s1_l + p1 * tmp7_l
                            du = h * tmp6_l * tmp3_l * s1_l - twostext[p] * p1 * tmp7_l
                            dh3 = (v_val * du - u_val * dv) / (v_val * v_val)

                            v_val = d2_l
                            dv = h * tmp4_l / s1_l + h * tmp5_l * s1_l

                            u_val = -h4_l / sigma * tmp4_l / s1_l - tmp9_l
                            du = -h * h4_l / sigma * tmp4_l / s1_l + twostext[p] * tmp9_l
                            dh5 = (v_val * du - u_val * dv) / (v_val * v_val)

                            u_val = h4_l / sigma * tmp5_l * s1_l + tmp9_l
                            du = -h * h4_l / sigma * tmp5_l * s1_l - twostext[p] * tmp9_l
                            dh6 = (v_val * du - u_val * dv) / (v_val * v_val)

                            da1 = h1_l / sigma * s2_l * s2_l + h2_l * s2_l * s1_l + h3_l * s2_l / s1_l +
                                  (1.0 - s2_l * s1_l) / (twostext[p] + h) * dh2 +
                                  (1.0 - s2_l / s1_l) / (twostext[p] - h) * dh3
                            da2 = h4_l / sigma * s2_l * s2_l + h5_l * s2_l * s1_l + h6_l * s2_l / s1_l +
                                  (1.0 - s2_l * s1_l) / (twostext[p] + h) * dh5 +
                                  (1.0 - s2_l / s1_l) / (twostext[p] - h) * dh6

                            d_ftid = -twostext[p] * h4_l / sigma * s2_l - h * h5_l * s1_l + h * h6_l / s1_l +
                                     dh5 * s1_l + dh6 / s1_l
                            d_fabd = -(dh2 + dh3) + (1.0 - albgrd[c, ib]) * twostext[p] * s2_l -
                                     (1.0 - albgri[c, ib]) * d_ftid
                            d_fabd_sun = (1.0 - omega[p, ib]) *
                                (twostext[p] * s2_l + 1.0 / avmu[p] * (da1 + da2))
                            d_fabd_sha = d_fabd - d_fabd_sun

                            surfalb.fabd_sun_z_patch[p, iv] = smooth_max(d_fabd_sun, 0.0) /
                                surfalb.fsun_z_patch[p, iv]
                            surfalb.fabd_sha_z_patch[p, iv] = smooth_max(d_fabd_sha, 0.0) /
                                (1.0 - surfalb.fsun_z_patch[p, iv])

                            # ---- Diffuse derivatives ----
                            u1_l = b - c1 / albgri[c, ib]
                            u2_l = b - c1 * albgri[c, ib]

                            tmp2_l = u1_l - avmu[p] * h
                            tmp3_l = u1_l + avmu[p] * h
                            d1_l = p1 * tmp2_l / s1_l - p2 * tmp3_l * s1_l
                            tmp4_l = u2_l + avmu[p] * h
                            tmp5_l = u2_l - avmu[p] * h
                            d2_l = tmp4_l / s1_l - tmp5_l * s1_l

                            h7_l = (c1 * tmp2_l) / (d1_l * s1_l)
                            h8_l = (-c1 * tmp3_l * s1_l) / d1_l
                            h9_l = tmp4_l / (d2_l * s1_l)
                            h10_l = (-tmp5_l * s1_l) / d2_l

                            a1_l = h7_l * (1.0 - s2_l * s1_l) / (twostext[p] + h) +
                                   h8_l * (1.0 - s2_l / s1_l) / (twostext[p] - h)
                            a2_l = h9_l * (1.0 - s2_l * s1_l) / (twostext[p] + h) +
                                   h10_l * (1.0 - s2_l / s1_l) / (twostext[p] - h)

                            v_val = d1_l
                            dv = h * p1 * tmp2_l / s1_l + h * p2 * tmp3_l * s1_l

                            u_val = c1 * tmp2_l / s1_l
                            du = h * c1 * tmp2_l / s1_l
                            dh7 = (v_val * du - u_val * dv) / (v_val * v_val)

                            u_val = -c1 * tmp3_l * s1_l
                            du = h * c1 * tmp3_l * s1_l
                            dh8 = (v_val * du - u_val * dv) / (v_val * v_val)

                            v_val = d2_l
                            dv = h * tmp4_l / s1_l + h * tmp5_l * s1_l

                            u_val = tmp4_l / s1_l
                            du = h * tmp4_l / s1_l
                            dh9 = (v_val * du - u_val * dv) / (v_val * v_val)

                            u_val = -tmp5_l * s1_l
                            du = h * tmp5_l * s1_l
                            dh10 = (v_val * du - u_val * dv) / (v_val * v_val)

                            da1 = h7_l * s2_l * s1_l + h8_l * s2_l / s1_l +
                                  (1.0 - s2_l * s1_l) / (twostext[p] + h) * dh7 +
                                  (1.0 - s2_l / s1_l) / (twostext[p] - h) * dh8
                            da2 = h9_l * s2_l * s1_l + h10_l * s2_l / s1_l +
                                  (1.0 - s2_l * s1_l) / (twostext[p] + h) * dh9 +
                                  (1.0 - s2_l / s1_l) / (twostext[p] - h) * dh10

                            d_ftii = -h * h9_l * s1_l + h * h10_l / s1_l + dh9 * s1_l + dh10 / s1_l
                            d_fabi = -(dh7 + dh8) - (1.0 - albgri[c, ib]) * d_ftii
                            d_fabi_sun = (1.0 - omega[p, ib]) / avmu[p] * (da1 + da2)
                            d_fabi_sha = d_fabi - d_fabi_sun

                            surfalb.fabi_sun_z_patch[p, iv] = smooth_max(d_fabi_sun, 0.0) /
                                surfalb.fsun_z_patch[p, iv]
                            surfalb.fabi_sha_z_patch[p, iv] = smooth_max(d_fabi_sha, 0.0) /
                                (1.0 - surfalb.fsun_z_patch[p, iv])
                        end  # iv loop
                    end  # nlevcan
                end  # ib == 1
            end  # !lSFonly
        end  # patch loop
    end  # waveband loop

    return nothing
end

# --------------------------------------------------------------------------
# surface_albedo!  (top-level driver)
# --------------------------------------------------------------------------

"""
    surface_albedo!(surfalb, con, grc, col, lun, patchdata,
                    canopystate, temperature, waterstatebulk,
                    waterdiagbulk, lakestate, aerosol,
                    mask_nourbanc, mask_nourbanp,
                    nextsw_cday, declinp1,
                    bounds_grc, bounds_col, bounds_patch,
                    pftcon_rhol, pftcon_rhos, pftcon_taul, pftcon_taus,
                    pftcon_xl;
                    use_SSRE=false, use_snicar_frc=false, use_fates=false,
                    use_subgrid_fluxes=true, do_sno_oc=false,
                    lakepuddling=false,
                    downscale_hillslope_meteorology=false)

Surface albedo and two-stream fluxes (top-level driver).

Ported from subroutine `SurfaceAlbedo` in `SurfaceAlbedoMod.F90`.
Uses SNICAR_RT for snow albedo when optics tables are loaded.
FATES calls are represented as stubs.
"""
function surface_albedo!(surfalb::SurfaceAlbedoData,
                         con::SurfaceAlbedoConstants,
                         grc::GridcellData,
                         col::ColumnData,
                         lun::LandunitData,
                         patchdata::PatchData,
                         canopystate::CanopyStateData,
                         temperature::TemperatureData,
                         waterstatebulk::WaterStateBulkData,
                         waterdiagbulk::WaterDiagnosticBulkData,
                         lakestate::LakeStateData,
                         aerosol::AerosolData,
                         mask_nourbanc::BitVector,
                         mask_nourbanp::BitVector,
                         nextsw_cday::Real,
                         declinp1::Real,
                         bounds_grc::UnitRange{Int},
                         bounds_col::UnitRange{Int},
                         bounds_patch::UnitRange{Int},
                         pftcon_rhol::Matrix{<:Real},
                         pftcon_rhos::Matrix{<:Real},
                         pftcon_taul::Matrix{<:Real},
                         pftcon_taus::Matrix{<:Real},
                         pftcon_xl::Vector{<:Real},
                         coszen_func::Function;
                         use_SSRE::Bool = false,
                         use_snicar_frc::Bool = false,
                         use_fates::Bool = false,
                         use_subgrid_fluxes::Bool = true,
                         do_sno_oc::Bool = false,
                         lakepuddling::Bool = false,
                         downscale_hillslope_meteorology::Bool = false)

    numrad = NUMRAD
    nlevcan = NLEVCAN
    nlevsno = varpar.nlevsno
    mpe = MPE_ALBEDO

    # Aliases
    coszen_grc = surfalb.coszen_grc
    coszen_col_arr = surfalb.coszen_col
    albgrd = surfalb.albgrd_col
    albgri = surfalb.albgri_col
    albsod = surfalb.albsod_col
    albsoi = surfalb.albsoi_col
    albd = surfalb.albd_patch
    albi = surfalb.albi_patch
    frac_sno = waterdiagbulk.frac_sno_col
    tlai = canopystate.tlai_patch
    tsai = canopystate.tsai_patch
    elai = canopystate.elai_patch
    esai = canopystate.esai_patch
    nrad = surfalb.nrad_patch
    ncan = surfalb.ncan_patch
    tlai_z = surfalb.tlai_z_patch
    tsai_z = surfalb.tsai_z_patch
    fsun_z = surfalb.fsun_z_patch

    # --- Cosine solar zenith angle ---
    for g in bounds_grc
        coszen_grc[g] = coszen_func(nextsw_cday, grc.lat[g], grc.lon[g], declinp1)
    end

    FT = eltype(coszen_col_arr)
    coszen_patch_arr = zeros(FT, length(bounds_patch))

    for c in bounds_col
        g = col.gridcell[c]
        if col.is_hillslope_column[c] && downscale_hillslope_meteorology
            # Hillslope meteorology downscaling (stub — use gridcell coszen)
            coszen_col_arr[c] = coszen_grc[g]
        else
            coszen_col_arr[c] = coszen_grc[g]
        end
    end

    for p in bounds_patch
        mask_nourbanp[p] || continue
        c = patchdata.column[p]
        coszen_patch_arr[p] = coszen_col_arr[c]
    end

    # --- Initialize output ---
    for ib in 1:numrad
        for c in bounds_col
            mask_nourbanc[c] || continue
            albsod[c, ib] = 0.0
            albsoi[c, ib] = 0.0
            albgrd[c, ib] = 0.0
            albgri[c, ib] = 0.0
            surfalb.albgrd_pur_col[c, ib] = 0.0
            surfalb.albgri_pur_col[c, ib] = 0.0
            surfalb.albgrd_bc_col[c, ib] = 0.0
            surfalb.albgri_bc_col[c, ib] = 0.0
            surfalb.albgrd_oc_col[c, ib] = 0.0
            surfalb.albgri_oc_col[c, ib] = 0.0
            surfalb.albgrd_dst_col[c, ib] = 0.0
            surfalb.albgri_dst_col[c, ib] = 0.0
            surfalb.albgrd_hst_col[c, ib] = SPVAL
            surfalb.albgri_hst_col[c, ib] = SPVAL
            surfalb.albgrd_pur_hst_col[c, ib] = SPVAL
            surfalb.albgri_pur_hst_col[c, ib] = SPVAL
            surfalb.albgrd_bc_hst_col[c, ib] = SPVAL
            surfalb.albgri_bc_hst_col[c, ib] = SPVAL
            surfalb.albgrd_oc_hst_col[c, ib] = SPVAL
            surfalb.albgri_oc_hst_col[c, ib] = SPVAL
            surfalb.albgrd_dst_hst_col[c, ib] = SPVAL
            surfalb.albgri_dst_hst_col[c, ib] = SPVAL
            surfalb.albsnd_hst2_col[c, ib] = SPVAL
            surfalb.albsni_hst2_col[c, ib] = SPVAL
            surfalb.flx_absdv_col[c, :] .= 0.0
            surfalb.flx_absdn_col[c, :] .= 0.0
            surfalb.flx_absiv_col[c, :] .= 0.0
            surfalb.flx_absin_col[c, :] .= 0.0
        end

        for p in bounds_patch
            mask_nourbanp[p] || continue
            albd[p, ib] = 1.0
            albi[p, ib] = 1.0
            surfalb.albd_hst_patch[p, ib] = SPVAL
            surfalb.albi_hst_patch[p, ib] = SPVAL
            if use_SSRE
                surfalb.albdSF_patch[p, ib] = 1.0
                surfalb.albiSF_patch[p, ib] = 1.0
            end
            surfalb.fabd_patch[p, ib] = 0.0
            surfalb.fabd_sun_patch[p, ib] = 0.0
            surfalb.fabd_sha_patch[p, ib] = 0.0
            surfalb.fabi_patch[p, ib] = 0.0
            surfalb.fabi_sun_patch[p, ib] = 0.0
            surfalb.fabi_sha_patch[p, ib] = 0.0
            surfalb.ftdd_patch[p, ib] = 0.0
            surfalb.ftid_patch[p, ib] = 0.0
            surfalb.ftii_patch[p, ib] = 0.0
        end
    end

    # --- Soil albedo ---
    soil_albedo!(surfalb, con, col, lun,
                 coszen_col_arr, temperature, waterstatebulk, lakestate,
                 mask_nourbanc, bounds_col; lakepuddling=lakepuddling)

    # --- Snow albedos via SNICAR radiative transfer ---
    nc = length(bounds_col)
    nlevsno_val = nlevsno

    # Check if SNICAR optics tables are loaded
    use_snicar = length(snicar_optics.ext_cff_mss_snw_drc) > 0 &&
                 any(x -> x != 0.0, @view snicar_optics.ext_cff_mss_snw_drc[1:min(10, end), 1:min(1, end)])

    # Prepare SNICAR input arrays
    h2osno_liq_snw = zeros(FT, nc, nlevsno_val)
    h2osno_ice_snw = zeros(FT, nc, nlevsno_val)
    h2osno_total_arr = zeros(FT, nc)
    snw_rds_int = fill(round(Int, snicar_params.snw_rds_min), nc, nlevsno_val)
    mss_cnc_aer = zeros(FT, nc, nlevsno_val, SNO_NBR_AER)
    albsfc_snw = zeros(FT, nc, numrad)

    # Compute total h2osno using the proper method (includes h2osno_no_layers_col)
    waterstate_calculate_total_h2osno!(waterstatebulk.ws, mask_nourbanc,
                                        bounds_col, col.snl, h2osno_total_arr)

    for c in bounds_col
        mask_nourbanc[c] || continue
        # Snow layer water content (first nlevsno slots are snow layers)
        for j in 1:nlevsno_val
            h2osno_liq_snw[c, j] = smooth_max(0.0, waterstatebulk.ws.h2osoi_liq_col[c, j])
            h2osno_ice_snw[c, j] = smooth_max(0.0, waterstatebulk.ws.h2osoi_ice_col[c, j])
        end
        # Snow grain radius (Float64 → Int, microns)
        for j in 1:nlevsno_val
            rds = waterdiagbulk.snw_rds_col[c, j]
            if !isnan(rds) && rds > 0.0
                snw_rds_int[c, j] = round(Int, clamp(rds, SNW_RDS_MIN_TBL, SNW_RDS_MAX_TBL))
            end
        end
        # Aerosol mass concentrations (8 species: bcphi, bcpho, ocphi, ocpho, dst1-4)
        if size(aerosol.mss_cnc_bcphi_col, 1) >= c
            for j in 1:nlevsno_val
                mss_cnc_aer[c, j, 1] = aerosol.mss_cnc_bcphi_col[c, j]
                mss_cnc_aer[c, j, 2] = aerosol.mss_cnc_bcpho_col[c, j]
                mss_cnc_aer[c, j, 3] = aerosol.mss_cnc_ocphi_col[c, j]
                mss_cnc_aer[c, j, 4] = aerosol.mss_cnc_ocpho_col[c, j]
                mss_cnc_aer[c, j, 5] = aerosol.mss_cnc_dst1_col[c, j]
                mss_cnc_aer[c, j, 6] = aerosol.mss_cnc_dst2_col[c, j]
                mss_cnc_aer[c, j, 7] = aerosol.mss_cnc_dst3_col[c, j]
                mss_cnc_aer[c, j, 8] = aerosol.mss_cnc_dst4_col[c, j]
            end
        end
        # Underlying surface albedo (diffuse soil albedo)
        for ib in 1:numrad
            albsfc_snw[c, ib] = albsoi[c, ib]
        end
    end

    # Output arrays
    albsnd = zeros(FT, nc, numrad)
    albsni = zeros(FT, nc, numrad)
    flx_absd_snw = zeros(FT, nc, nlevsno_val + 1, numrad)
    flx_absi_snw = zeros(FT, nc, nlevsno_val + 1, numrad)

    if use_snicar
        # Call SNICAR_RT for direct beam (flg_slr=1)
        snicar_rt!(coszen_col_arr, 1,
                   h2osno_liq_snw, h2osno_ice_snw, h2osno_total_arr,
                   snw_rds_int, mss_cnc_aer, albsfc_snw,
                   col.snl, frac_sno,
                   albsnd, flx_absd_snw, nlevsno_val;
                   mask_nourbanc=mask_nourbanc)

        # Call SNICAR_RT for diffuse (flg_slr=2)
        snicar_rt!(coszen_col_arr, 2,
                   h2osno_liq_snw, h2osno_ice_snw, h2osno_total_arr,
                   snw_rds_int, mss_cnc_aer, albsfc_snw,
                   col.snl, frac_sno,
                   albsni, flx_absi_snw, nlevsno_val;
                   mask_nourbanc=mask_nourbanc)
    else
        # Fallback: age-based snow albedo (used when SNICAR optics not loaded)
        snow_persist = waterstatebulk.snow_persistence_col
        for c in bounds_col
            mask_nourbanc[c] || continue
            if coszen_col_arr[c] > 0.0 && frac_sno[c] > 0.0
                age_days = c <= length(snow_persist) ? snow_persist[c] / 86400.0 : 0.0
                fage = 1.0 - exp(-age_days / 5.0)
                albsnd[c, 1] = clamp(0.95 * (1.0 - 0.15 * fage), 0.1, 0.99)
                albsni[c, 1] = clamp(0.95 * (1.0 - 0.15 * fage), 0.1, 0.99)
                albsnd[c, 2] = clamp(0.65 * (1.0 - 0.50 * fage), 0.1, 0.99)
                albsni[c, 2] = clamp(0.65 * (1.0 - 0.50 * fage), 0.1, 0.99)
            end
        end
    end

    # --- Ground albedos (snow-fraction weighting) ---
    for ib in 1:numrad
        for c in bounds_col
            mask_nourbanc[c] || continue
            if coszen_col_arr[c] > 0.0
                albgrd[c, ib] = albsod[c, ib] * (1.0 - frac_sno[c]) + albsnd[c, ib] * frac_sno[c]
                albgri[c, ib] = albsoi[c, ib] * (1.0 - frac_sno[c]) + albsni[c, ib] * frac_sno[c]

                if use_snicar_frc
                    surfalb.albgrd_bc_col[c, ib] = albsod[c, ib] * (1.0 - frac_sno[c])
                    surfalb.albgri_bc_col[c, ib] = albsoi[c, ib] * (1.0 - frac_sno[c])
                    surfalb.albgrd_dst_col[c, ib] = albsod[c, ib] * (1.0 - frac_sno[c])
                    surfalb.albgri_dst_col[c, ib] = albsoi[c, ib] * (1.0 - frac_sno[c])
                    surfalb.albgrd_pur_col[c, ib] = albsod[c, ib] * (1.0 - frac_sno[c])
                    surfalb.albgri_pur_col[c, ib] = albsoi[c, ib] * (1.0 - frac_sno[c])
                    if do_sno_oc
                        surfalb.albgrd_oc_col[c, ib] = albsod[c, ib] * (1.0 - frac_sno[c])
                        surfalb.albgri_oc_col[c, ib] = albsoi[c, ib] * (1.0 - frac_sno[c])
                    end
                end
            end
        end
    end

    # --- Map SNICAR flux absorption to surfalb fields ---
    if use_snicar
        for c in bounds_col
            mask_nourbanc[c] || continue
            l = col.landunit[c]
            is_lake = lun.lakpoi[l]
            for i in 1:(nlevsno_val + 1)
                if use_subgrid_fluxes && !is_lake
                    # Subgrid fluxes: SNICAR output used directly for snow-covered area
                    surfalb.flx_absdv_col[c, i] = flx_absd_snw[c, i, IVIS]
                    surfalb.flx_absdn_col[c, i] = flx_absd_snw[c, i, INIR]
                    surfalb.flx_absiv_col[c, i] = flx_absi_snw[c, i, IVIS]
                    surfalb.flx_absin_col[c, i] = flx_absi_snw[c, i, INIR]
                else
                    # No subgrid fluxes or lake: weight by snow fraction
                    fsnow = frac_sno[c]
                    for (flx_field, flx_snw, alb_s, alb_d, ib) in [
                        (:flx_absdv_col, flx_absd_snw, albsnd, albsod, IVIS),
                        (:flx_absdn_col, flx_absd_snw, albsnd, albsod, INIR),
                        (:flx_absiv_col, flx_absi_snw, albsni, albsoi, IVIS),
                        (:flx_absin_col, flx_absi_snw, albsni, albsoi, INIR),
                    ]
                        f_snw = flx_snw[c, i, ib]
                        a_snw = alb_s[c, ib]
                        a_soil = alb_d[c, ib]
                        flx_out = getfield(surfalb, flx_field)
                        if a_snw < 1.0
                            flx_out[c, i] = f_snw * fsnow +
                                (1.0 - fsnow) * (1.0 - a_soil) * (f_snw / max(1.0 - a_snw, 1.0e-6))
                        else
                            flx_out[c, i] = 0.0
                        end
                    end
                end
            end
        end
    end

    # --- Snow albedo history diagnostics ---
    for ib in 1:numrad
        for c in bounds_col
            mask_nourbanc[c] || continue
            if coszen_col_arr[c] > 0.0 && h2osno_total_arr[c] > 0.0
                surfalb.albsnd_hst_col[c, ib] = albsnd[c, ib]
                surfalb.albsni_hst_col[c, ib] = albsni[c, ib]
            else
                surfalb.albsnd_hst_col[c, ib] = 0.0
                surfalb.albsni_hst_col[c, ib] = 0.0
            end
        end
    end

    # --- Canopy layering ---
    dincmax = 0.25
    for p in bounds_patch
        mask_nourbanp[p] || continue

        if nlevcan == 1
            nrad[p] = 1
            ncan[p] = 1
            tlai_z[p, 1] = elai[p]
            tsai_z[p, 1] = esai[p]
        elseif nlevcan > 1
            if elai[p] + esai[p] == 0.0
                nrad[p] = 0
            else
                dincmax_sum = 0.0
                for iv in 1:nlevcan
                    dincmax_sum += dincmax
                    if ((elai[p] + esai[p]) - dincmax_sum) > 1.0e-06
                        nrad[p] = iv
                        dinc = dincmax
                        tlai_z[p, iv] = dinc * elai[p] / smooth_max(elai[p] + esai[p], mpe)
                        tsai_z[p, iv] = dinc * esai[p] / smooth_max(elai[p] + esai[p], mpe)
                    else
                        nrad[p] = iv
                        dinc = dincmax - (dincmax_sum - (elai[p] + esai[p]))
                        tlai_z[p, iv] = dinc * elai[p] / smooth_max(elai[p] + esai[p], mpe)
                        tsai_z[p, iv] = dinc * esai[p] / smooth_max(elai[p] + esai[p], mpe)
                        break
                    end
                end

                # Minimum of 4 canopy layers
                if nrad[p] < 4
                    nrad[p] = 4
                    for iv in 1:nrad[p]
                        tlai_z[p, iv] = elai[p] / nrad[p]
                        tsai_z[p, iv] = esai[p] / nrad[p]
                    end
                end
            end

            # Repeat for buried canopy layers
            blai = tlai[p] - elai[p]
            bsai = tsai[p] - esai[p]
            if blai + bsai == 0.0
                ncan[p] = nrad[p]
            else
                dincmax_sum = 0.0
                for iv in (nrad[p] + 1):nlevcan
                    dincmax_sum += dincmax
                    if ((blai + bsai) - dincmax_sum) > 1.0e-06
                        ncan[p] = iv
                        dinc = dincmax
                        tlai_z[p, iv] = dinc * blai / smooth_max(blai + bsai, mpe)
                        tsai_z[p, iv] = dinc * bsai / smooth_max(blai + bsai, mpe)
                    else
                        ncan[p] = iv
                        dinc = dincmax - (dincmax_sum - (blai + bsai))
                        tlai_z[p, iv] = dinc * blai / smooth_max(blai + bsai, mpe)
                        tsai_z[p, iv] = dinc * bsai / smooth_max(blai + bsai, mpe)
                        break
                    end
                end
            end
        end
    end

    # --- Zero fluxes for active canopy layers ---
    for p in bounds_patch
        mask_nourbanp[p] || continue
        for iv in 1:nrad[p]
            surfalb.fabd_sun_z_patch[p, iv] = 0.0
            surfalb.fabd_sha_z_patch[p, iv] = 0.0
            surfalb.fabi_sun_z_patch[p, iv] = 0.0
            surfalb.fabi_sha_z_patch[p, iv] = 0.0
            fsun_z[p, iv] = 0.0
        end
    end

    # --- Default vcmax scaling (coszen <= 0 case) ---
    extkn = 0.30
    for p in bounds_patch
        mask_nourbanp[p] || continue
        if nlevcan == 1
            surfalb.vcmaxcintsun_patch[p] = 0.0
            surfalb.vcmaxcintsha_patch[p] = (1.0 - exp(-extkn * elai[p])) / extkn
            if elai[p] > 0.0
                surfalb.vcmaxcintsha_patch[p] /= elai[p]
            else
                surfalb.vcmaxcintsha_patch[p] = 0.0
            end
        elseif nlevcan > 1
            surfalb.vcmaxcintsun_patch[p] = 0.0
            surfalb.vcmaxcintsha_patch[p] = 0.0
        end
    end

    # --- Create solar-vegetated filter ---
    mask_vegsol = falses(length(bounds_patch))
    mask_novegsol = falses(length(bounds_patch))

    for p in bounds_patch
        mask_nourbanp[p] || continue
        if coszen_patch_arr[p] > 0.0
            l_p = patchdata.landunit[p]
            if (lun.itype[l_p] == ISTSOIL || lun.itype[l_p] == ISTCROP) &&
               (elai[p] + esai[p]) > 0.0
                mask_vegsol[p] = true
            else
                mask_novegsol[p] = true
            end
        end
    end

    # --- Weight reflectance/transmittance by LAI and SAI ---
    np = length(bounds_patch)
    wl = zeros(FT, np)
    ws_arr = zeros(FT, np)
    rho = zeros(FT, np, numrad)
    tau = zeros(FT, np, numrad)

    for p in bounds_patch
        mask_vegsol[p] || continue
        wl[p] = elai[p] / smooth_max(elai[p] + esai[p], mpe)
        ws_arr[p] = esai[p] / smooth_max(elai[p] + esai[p], mpe)
    end

    for ib in 1:numrad
        for p in bounds_patch
            mask_vegsol[p] || continue
            itype = patchdata.itype[p] + 1  # +1: Fortran 0-based PFT → Julia 1-based array
            rho[p, ib] = smooth_max(pftcon_rhol[itype, ib] * wl[p] + pftcon_rhos[itype, ib] * ws_arr[p], mpe)
            tau[p, ib] = smooth_max(pftcon_taul[itype, ib] * wl[p] + pftcon_taus[itype, ib] * ws_arr[p], mpe)
        end
    end

    # --- Two-stream calculation ---
    if !use_fates
        two_stream!(surfalb, patchdata, col, canopystate, temperature,
                    waterdiagbulk, coszen_patch_arr, rho, tau,
                    pftcon_xl, mask_vegsol, bounds_patch)

        if use_SSRE
            if nlevcan > 1
                error("use_ssre option was NOT developed with allowance for multi-layer canopy")
            end
            two_stream!(surfalb, patchdata, col, canopystate, temperature,
                        waterdiagbulk, coszen_patch_arr, rho, tau,
                        pftcon_xl, mask_vegsol, bounds_patch; SFonly=true)
        end
    end
    # (FATES canopy radiation stub — not yet ported)

    # --- Non-vegetated patches where coszen > 0 ---
    for ib in 1:numrad
        for p in bounds_patch
            mask_novegsol[p] || continue
            c = patchdata.column[p]
            surfalb.fabd_patch[p, ib] = 0.0
            surfalb.fabd_sun_patch[p, ib] = 0.0
            surfalb.fabd_sha_patch[p, ib] = 0.0
            surfalb.fabi_patch[p, ib] = 0.0
            surfalb.fabi_sun_patch[p, ib] = 0.0
            surfalb.fabi_sha_patch[p, ib] = 0.0
            surfalb.ftdd_patch[p, ib] = 1.0
            surfalb.ftid_patch[p, ib] = 0.0
            surfalb.ftii_patch[p, ib] = 1.0
            albd[p, ib] = albgrd[c, ib]
            albi[p, ib] = albgri[c, ib]
            if use_SSRE
                surfalb.albdSF_patch[p, ib] = albsod[c, ib]
                surfalb.albiSF_patch[p, ib] = albsoi[c, ib]
            end
        end
    end

    # --- History output variables ---
    for ib in 1:numrad
        for c in bounds_col
            mask_nourbanc[c] || continue
            if coszen_col_arr[c] > 0.0
                surfalb.albgrd_hst_col[c, ib] = albgrd[c, ib]
                surfalb.albgri_hst_col[c, ib] = albgri[c, ib]
                surfalb.albgrd_pur_hst_col[c, ib] = surfalb.albgrd_pur_col[c, ib]
                surfalb.albgri_pur_hst_col[c, ib] = surfalb.albgri_pur_col[c, ib]
                surfalb.albgrd_bc_hst_col[c, ib] = surfalb.albgrd_bc_col[c, ib]
                surfalb.albgri_bc_hst_col[c, ib] = surfalb.albgri_bc_col[c, ib]
                surfalb.albgrd_oc_hst_col[c, ib] = surfalb.albgrd_oc_col[c, ib]
                surfalb.albgri_oc_hst_col[c, ib] = surfalb.albgri_oc_col[c, ib]
                surfalb.albgrd_dst_hst_col[c, ib] = surfalb.albgrd_dst_col[c, ib]
                surfalb.albgri_dst_hst_col[c, ib] = surfalb.albgri_dst_col[c, ib]
                if h2osno_total_arr[c] > 0.0
                    surfalb.albsnd_hst2_col[c, ib] = surfalb.albsnd_hst_col[c, ib]
                    surfalb.albsni_hst2_col[c, ib] = surfalb.albsni_hst_col[c, ib]
                end
            end
        end

        for p in bounds_patch
            mask_nourbanp[p] || continue
            if coszen_patch_arr[p] > 0.0
                surfalb.albd_hst_patch[p, ib] = albd[p, ib]
                surfalb.albi_hst_patch[p, ib] = albi[p, ib]
            end
        end
    end

    return nothing
end
