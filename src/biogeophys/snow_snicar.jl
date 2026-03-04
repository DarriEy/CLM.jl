# ==========================================================================
# Ported from: src/biogeophys/SnowSnicarMod.F90
# Calculate albedo of snow containing impurities and the evolution of snow
# effective radius.
#
# Public functions:
#   snicar_rt!          — Snow albedo and vertically-resolved solar absorption
#   snowage_grain!      — Snow effective grain size evolution
#   fresh_snow_radius   — Temperature-dependent fresh snow grain radius
#   snow_optics_init!   — Initialize snow optics lookup tables (stub)
#   snowage_init!       — Initialize snow aging lookup tables (stub)
#   piecewise_linear_interp1d — 1D piecewise linear interpolation
# ==========================================================================

# ---- Module-level constants ----

const SNO_NBR_AER = 8       # number of aerosol species in snowpack
const DO_SNO_AER = true     # include aerosols in snowpack radiative calculations

const DEFAULT_NUMBER_BANDS = 5
const HIGHEST_DEFAULT_BAND = 5
const SEC_HIGHEST_DEFAULT_BAND = 4
const HIGH_NUMBER_BANDS = 480

const IDX_MIE_SNW_MX = 1471  # number of effective radius indices in Mie lookup table
const IDX_T_MAX       = 11
const IDX_T_MIN       = 1
const IDX_TGRD_MAX    = 31
const IDX_TGRD_MIN    = 1
const IDX_RHOS_MAX    = 8
const IDX_RHOS_MIN    = 1

const SNW_RDS_MAX_TBL = 1500   # maximum effective radius defined in Mie lookup table [microns]
const SNW_RDS_MIN_TBL = 30     # minimum effective radius defined in Mie lookup table [microns]
const SNW_RDS_MAX     = 1500.0 # maximum allowed snow effective radius [microns]

const MIN_SNW = 1.0e-30     # minimum snow mass required for SNICAR RT calculation [kg/m2]

const C1_LIQ_BRUN89 = 0.0   # constant for liquid water grain growth [m3/s], zeroed for dry snow aging

const TIM_CNS_BC_RMV  = 2.2e-8  # time constant for removal of BC in snow on sea-ice [s-1]
const TIM_CNS_OC_RMV  = 2.2e-8  # time constant for removal of OC in snow on sea-ice [s-1]
const TIM_CNS_DST_RMV = 2.2e-8  # time constant for removal of dust in snow on sea-ice [s-1]

# ---- SNICAR Parameters struct ----

Base.@kwdef mutable struct SnicarParams
    xdrdt::Float64 = 1.0              # arbitrary factor applied to snow aging rate [-]
    snw_rds_refrz::Float64 = 1000.0   # effective radius of re-frozen snow [microns]
    C2_liq_Brun89::Float64 = 4.22e-13 # constant for liquid water grain growth [m3/s]
    fresh_snw_rds_max::Float64 = 204.526 # maximum warm fresh snow effective radius [microns]
    snw_rds_min::Float64 = 54.526     # minimum allowed snow effective radius [microns]
end

const snicar_params = SnicarParams()

# ---- Snow Optics Data ----

Base.@kwdef mutable struct SnicarOpticsData
    # direct-beam weighted ice optical properties (idx_Mie_snw_mx, numrad_snw)
    ss_alb_snw_drc::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    asm_prm_snw_drc::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    ext_cff_mss_snw_drc::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)

    # diffuse radiation weighted ice optical properties
    ss_alb_snw_dfs::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    asm_prm_snw_dfs::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    ext_cff_mss_snw_dfs::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)

    # hydrophilic BC (numrad_snw)
    ss_alb_bc_hphil::Vector{Float64} = Float64[]
    asm_prm_bc_hphil::Vector{Float64} = Float64[]
    ext_cff_mss_bc_hphil::Vector{Float64} = Float64[]

    # hydrophobic BC
    ss_alb_bc_hphob::Vector{Float64} = Float64[]
    asm_prm_bc_hphob::Vector{Float64} = Float64[]
    ext_cff_mss_bc_hphob::Vector{Float64} = Float64[]

    # hydrophilic OC
    ss_alb_oc_hphil::Vector{Float64} = Float64[]
    asm_prm_oc_hphil::Vector{Float64} = Float64[]
    ext_cff_mss_oc_hphil::Vector{Float64} = Float64[]

    # hydrophobic OC
    ss_alb_oc_hphob::Vector{Float64} = Float64[]
    asm_prm_oc_hphob::Vector{Float64} = Float64[]
    ext_cff_mss_oc_hphob::Vector{Float64} = Float64[]

    # dust species 1-4
    ss_alb_dst1::Vector{Float64} = Float64[]
    asm_prm_dst1::Vector{Float64} = Float64[]
    ext_cff_mss_dst1::Vector{Float64} = Float64[]

    ss_alb_dst2::Vector{Float64} = Float64[]
    asm_prm_dst2::Vector{Float64} = Float64[]
    ext_cff_mss_dst2::Vector{Float64} = Float64[]

    ss_alb_dst3::Vector{Float64} = Float64[]
    asm_prm_dst3::Vector{Float64} = Float64[]
    ext_cff_mss_dst3::Vector{Float64} = Float64[]

    ss_alb_dst4::Vector{Float64} = Float64[]
    asm_prm_dst4::Vector{Float64} = Float64[]
    ext_cff_mss_dst4::Vector{Float64} = Float64[]

    # downward solar radiation spectral weights
    flx_wgt_dir::Vector{Float64} = Float64[]  # direct
    flx_wgt_dif::Vector{Float64} = Float64[]  # diffuse
end

const snicar_optics = SnicarOpticsData()

# ---- Snow Aging Data ----

Base.@kwdef mutable struct SnicarAgingData
    snowage_tau::Array{Float64,3} = Array{Float64}(undef, 0, 0, 0)    # (idx_rhos_max, idx_Tgrd_max, idx_T_max)
    snowage_kappa::Array{Float64,3} = Array{Float64}(undef, 0, 0, 0)
    snowage_drdt0::Array{Float64,3} = Array{Float64}(undef, 0, 0, 0)
end

const snicar_aging = SnicarAgingData()

# ---- Gaussian quadrature constants for diffuse integration ----
const SNICAR_NGMAX = 8
const SNICAR_DIFGAUSPT = [0.9894009, 0.9445750, 0.8656312, 0.7554044,
                          0.6178762, 0.4580168, 0.2816036, 0.0950125]
const SNICAR_DIFGAUSWT = [0.0271525, 0.0622535, 0.0951585, 0.1246290,
                          0.1495960, 0.1691565, 0.1826034, 0.1894506]

# ---- Aspherical ice particle constants (He et al. 2017) ----
const SNICAR_SEVEN_BANDS = 7
const SNICAR_G_WVL = [0.25, 0.70, 1.41, 1.90, 2.50, 3.50, 4.00, 5.00]
const SNICAR_G_B0 = [9.76029e-1, 9.67798e-1, 1.00111, 1.00224,
                     9.64295e-1, 9.97475e-1, 9.97475e-1]
const SNICAR_G_B1 = [5.21042e-1, 4.96181e-1, 1.83711e-1, 1.37082e-1,
                     5.50598e-2, 8.48743e-2, 8.48743e-2]
const SNICAR_G_B2 = [-2.66792e-4, 1.14088e-3, 2.37011e-4, -2.35905e-4,
                      8.40449e-4, -4.71484e-4, -4.71484e-4]
const SNICAR_G_F07_C2 = [1.349959e-1, 1.115697e-1, 9.853958e-2, 5.557793e-2,
                         -1.233493e-1, 0.0, 0.0]
const SNICAR_G_F07_C1 = [-3.987320e-1, -3.723287e-1, -3.924784e-1, -3.259404e-1,
                          4.429054e-2, -1.726586e-1, -1.726586e-1]
const SNICAR_G_F07_C0 = [7.938904e-1, 8.030084e-1, 8.513932e-1, 8.692241e-1,
                         7.085850e-1, 6.412701e-1, 6.412701e-1]
const SNICAR_G_F07_P2 = [3.165543e-3, 2.014810e-3, 1.780838e-3, 6.987734e-4,
                         -1.882932e-2, -2.277872e-2, -2.277872e-2]
const SNICAR_G_F07_P1 = [1.140557e-1, 1.143152e-1, 1.143814e-1, 1.071238e-1,
                         1.353873e-1, 1.914431e-1, 1.914431e-1]
const SNICAR_G_F07_P0 = [5.292852e-1, 5.425909e-1, 5.601598e-1, 6.023407e-1,
                         6.473899e-1, 4.634944e-1, 4.634944e-1]

# Shape factor defaults
const FS_SPHD_DEFAULT = 0.929
const FS_HEX_REF = 0.788
const FS_KOCH_DEFAULT = 0.712
const AR_TMP_DEFAULT_1 = 0.5
const AR_TMP_DEFAULT_2 = 2.5

# ---- BC-snow internal mixing constants (He et al. 2017) ----
const SNICAR_SIXTEEN_BANDS = 16
const BCINT_WVL = [0.20, 0.25, 0.30, 0.33, 0.36, 0.40, 0.44, 0.48,
                   0.52, 0.57, 0.64, 0.69, 0.75, 0.78, 0.87, 1.0, 1.2]
const BCINT_D0 = [2.48045, 4.70305, 4.68619, 4.67369, 4.65040,
                  2.40364, 7.95408e-1, 2.92745e-1, 8.63396e-2, 2.76299e-2,
                  1.40864e-2, 8.65705e-3, 6.12971e-3, 4.45697e-3, 3.06648e-2,
                  7.96544e-1]
const BCINT_D1 = [9.77209e-1, 9.73317e-1, 9.79650e-1, 9.84579e-1, 9.93537e-1,
                  9.95955e-1, 9.95218e-1, 9.74284e-1, 9.81193e-1, 9.81239e-1,
                  9.55515e-1, 9.10491e-1, 8.74196e-1, 8.27238e-1, 4.82870e-1,
                  4.36649e-2]
const BCINT_D2 = [3.95960e-1, 2.04820e-1, 2.07410e-1, 2.09390e-1, 2.13030e-1,
                  4.18570e-1, 1.29682, 3.75514, 1.27372e+1, 3.93293e+1,
                  8.78918e+1, 1.86969e+2, 3.45600e+2, 7.08637e+2, 1.41067e+3,
                  2.57288e+2]
const DEN_BC = 1.7
const DEN_BC_TARGET = 1.49
const RE_BC = 0.045
const RADIUS_1_BC = 0.1
const RADIUS_2_BC = 0.05
const SNICAR_THREE_BANDS = 3
const BCINT_M = [-0.8724, -0.1866, -0.0046]
const BCINT_N = [-0.0072, -0.1918, -0.5177]

# ---- Dust-snow internal mixing constants (He et al. 2019) ----
const SNICAR_SIZE_BINS = 6
const DSTINT_WVL = [0.2, 0.2632, 0.3448, 0.4415, 0.625, 0.7782, 1.2422]
const DSTINT_A1 = [-2.1307e+1, -1.5815e+1, -9.2880, 1.1115, 1.0307, 1.0185]
const DSTINT_A2 = [1.1746e+2, 9.3241e+1, 4.0605e+1, 3.7389e-1, 1.4800e-2, 2.8921e-4]
const DSTINT_A3 = [9.9701e-1, 9.9781e-1, 9.9848e-1, 1.0035, 1.0024, 1.0356]

const ENH_OMG_MAX = 1.0e5
const KG_KG_TO_PPM = 1.0e6
const KG_TO_UG = 1.0e9

# SZA parameterization constants
const SZA_A0 =  0.085730
const SZA_A1 = -0.630883
const SZA_A2 =  1.303723
const SZA_B0 =  1.467291
const SZA_B1 = -3.338043
const SZA_B2 =  6.807489
const SNICAR_PUNY = 1.0e-11
const MU_75 = 0.2588  # cosine of 75 degrees

# --------------------------------------------------------------------------
# piecewise_linear_interp1d
# --------------------------------------------------------------------------

"""
    piecewise_linear_interp1d(xd, yd, xi)

Piecewise linear interpolation for 1D data.
Returns interpolated value at `xi` given data points `(xd, yd)`.

Ported from `piecewise_linear_interp1d` in `SnowSnicarMod.F90`.
"""
function piecewise_linear_interp1d(xd::AbstractVector{Float64}, yd::AbstractVector{Float64}, xi::Float64)
    nd = length(xd)
    yi = 0.0

    if nd == 1
        return yd[1]
    end

    if xi < xd[1]
        t = (xi - xd[1]) / (xd[2] - xd[1])
        yi = (1.0 - t) * yd[1] + t * yd[2]
    elseif xi > xd[nd]
        t = (xi - xd[nd-1]) / (xd[nd] - xd[nd-1])
        yi = (1.0 - t) * yd[nd-1] + t * yd[nd]
    else
        for k in 2:nd
            if xd[k-1] <= xi && xi <= xd[k]
                t = (xi - xd[k-1]) / (xd[k] - xd[k-1])
                yi = (1.0 - t) * yd[k-1] + t * yd[k]
                break
            end
        end
    end

    return yi
end

# --------------------------------------------------------------------------
# fresh_snow_radius
# --------------------------------------------------------------------------

"""
    fresh_snow_radius(forc_t; params=snicar_params)

Returns fresh snow grain radius [microns], linearly dependent on temperature.
`forc_t` is the atmospheric forcing temperature [K].

Ported from `FreshSnowRadius` in `SnowSnicarMod.F90`.
"""
function fresh_snow_radius(forc_t::Float64; params::SnicarParams=snicar_params)
    tmin = TFRZ - 30.0
    tmax = TFRZ

    if params.fresh_snw_rds_max <= params.snw_rds_min
        return params.snw_rds_min
    end

    gs_max = params.fresh_snw_rds_max
    gs_min = params.snw_rds_min

    if forc_t < tmin
        return gs_min
    elseif forc_t > tmax
        return gs_max
    else
        return (tmax - forc_t) / (tmax - tmin) * gs_min +
               (forc_t - tmin) / (tmax - tmin) * gs_max
    end
end

# --------------------------------------------------------------------------
# snicar_rt!
# --------------------------------------------------------------------------

"""
    snicar_rt!(coszen, flg_slr_in, h2osno_liq, h2osno_ice, h2osno_total,
               snw_rds, mss_cnc_aer_in, albsfc, snl, frac_sno,
               albout, flx_abs, nlevsno;
               optics=snicar_optics, params=snicar_params,
               mask_nourbanc=nothing)

Determine reflectance of, and vertically-resolved solar absorption in,
snow with impurities using the Adding-Doubling radiative transfer solver
(SNICAR-ADv3).

Ported from `SNICAR_RT` in `SnowSnicarMod.F90`.

# Arguments
- `coszen::Vector{Float64}`: cosine of solar zenith angle (col)
- `flg_slr_in::Int`: 1=direct-beam, 2=diffuse incident flux
- `h2osno_liq::Matrix{Float64}`: liquid water content (col,lyr) [kg/m2], 1-based snow layers
- `h2osno_ice::Matrix{Float64}`: ice content (col,lyr) [kg/m2], 1-based snow layers
- `h2osno_total::Vector{Float64}`: total snow content (col) [kg/m2]
- `snw_rds::Matrix{Int}`: snow effective radius (col,lyr) [microns], 1-based snow layers
- `mss_cnc_aer_in::Array{Float64,3}`: aerosol mass concentration (col,lyr,aer) [kg/kg]
- `albsfc::Matrix{Float64}`: albedo of underlying surface (col,bnd) [frc]
- `snl::Vector{Int}`: negative number of snow layers (col)
- `frac_sno::Vector{Float64}`: fraction of ground covered by snow (col)
- `albout::Matrix{Float64}`: output snow albedo (col,bnd) [frc]
- `flx_abs::Array{Float64,3}`: output absorbed flux (col,lyr,bnd) per unit incident flux
- `nlevsno::Int`: maximum number of snow layers
- `optics::SnicarOpticsData`: snow optics lookup tables
- `params::SnicarParams`: SNICAR parameters
- `mask_nourbanc::Union{BitVector,Nothing}`: column mask (if nothing, process all)
"""
function snicar_rt!(coszen::Vector{Float64},
                    flg_slr_in::Int,
                    h2osno_liq::Matrix{Float64},
                    h2osno_ice::Matrix{Float64},
                    h2osno_total::Vector{Float64},
                    snw_rds::Matrix{Int},
                    mss_cnc_aer_in::Array{Float64,3},
                    albsfc::Matrix{Float64},
                    snl::Vector{Int},
                    frac_sno::Vector{Float64},
                    albout::Matrix{Float64},
                    flx_abs::Array{Float64,3},
                    nlevsno::Int;
                    optics::SnicarOpticsData=snicar_optics,
                    params::SnicarParams=snicar_params,
                    mask_nourbanc::Union{BitVector,Nothing}=nothing)

    numrad_snw = varctl.snicar_numrad_snw
    snicar_snw_shape = varctl.snicar_snw_shape
    snicar_snobc_intmix = varctl.snicar_snobc_intmix
    snicar_snodst_intmix = varctl.snicar_snodst_intmix
    pi_val = Float64(π)

    # Snow layer offset: Fortran index i maps to Julia index i + nlevsno
    joff = nlevsno

    # Determine NIR band beginning
    if numrad_snw == DEFAULT_NUMBER_BANDS
        nir_bnd_bgn = 2
    elseif numrad_snw == HIGH_NUMBER_BANDS
        nir_bnd_bgn = 51
    else
        error("SNICAR ERROR: unknown snicar_numrad_snw value: $numrad_snw")
    end
    nir_bnd_end = numrad_snw

    # Band center wavelengths
    wvl_ct = zeros(numrad_snw)
    if numrad_snw == DEFAULT_NUMBER_BANDS
        wvl_ct .= [0.5, 0.85, 1.1, 1.35, 3.25]
    elseif numrad_snw == HIGH_NUMBER_BANDS
        for igb in 1:numrad_snw
            wvl_ct[igb] = 0.205 + 0.01 * (igb - 1.0)
        end
    end

    # Nonspherical snow grain constants
    g_wvl_ct = zeros(SNICAR_SEVEN_BANDS)
    for igb in 1:SNICAR_SEVEN_BANDS
        g_wvl_ct[igb] = SNICAR_G_WVL[igb+1] * 0.5 + SNICAR_G_WVL[igb] * 0.5
    end
    dstint_wvl_ct = zeros(SNICAR_SIZE_BINS)
    for idb in 1:SNICAR_SIZE_BINS
        dstint_wvl_ct[idb] = DSTINT_WVL[idb+1] * 0.5 + DSTINT_WVL[idb] * 0.5
    end
    bcint_wvl_ct = zeros(SNICAR_SIXTEEN_BANDS)
    for ibb in 1:SNICAR_SIXTEEN_BANDS
        bcint_wvl_ct[ibb] = BCINT_WVL[ibb+1] * 0.5 + BCINT_WVL[ibb] * 0.5
    end

    # Always use Delta approximation
    DELTA = 1

    # Constants used in Adding-Doubling solver
    c0 = 0.0; c1 = 1.0; c3 = 3.0; c4 = 4.0; c6 = 6.0
    cp01 = 0.01; cp5 = 0.5; cp75 = 0.75; c1p5 = 1.5
    trmin = 0.001; argmax_val = 10.0

    ncols = length(coszen)

    # Loop over columns
    for c_idx in 1:ncols
        if mask_nourbanc !== nothing && !mask_nourbanc[c_idx]
            continue
        end

        # Zero absorbed radiative fluxes
        for i in 1:(nlevsno+1)
            for bnd in 1:NUMRAD
                flx_abs[c_idx, i, bnd] = 0.0
            end
        end

        h2osno_lcl = h2osno_total[c_idx]

        if coszen[c_idx] > 0.0 && h2osno_lcl > MIN_SNW

            # Local snow layer setup
            if snl[c_idx] > -1
                flg_nosnl = 1
                snl_lcl = -1
                # Allocate local arrays for snow layers (1:nlevsno), index offset joff
                h2osno_ice_lcl = zeros(nlevsno)
                h2osno_liq_lcl = zeros(nlevsno)
                snw_rds_lcl = zeros(Int, nlevsno)
                h2osno_ice_lcl[joff] = h2osno_lcl  # layer 0 -> index joff
                h2osno_liq_lcl[joff] = 0.0
                snw_rds_lcl[joff] = round(Int, params.snw_rds_min)
            else
                flg_nosnl = 0
                snl_lcl = snl[c_idx]
                h2osno_liq_lcl = h2osno_liq[c_idx, :]
                h2osno_ice_lcl = h2osno_ice[c_idx, :]
                snw_rds_lcl = snw_rds[c_idx, :]
            end

            snl_btm = 0         # Fortran index of bottom snow layer
            snl_top = snl_lcl + 1  # Fortran index of top snow layer

            # Julia indices for snow layer range
            snl_btm_j = snl_btm + joff
            snl_top_j = snl_top + joff

            # Local aerosol array
            mss_cnc_aer_lcl = zeros(nlevsno, SNO_NBR_AER)
            for j in 1:SNO_NBR_AER
                mss_cnc_aer_lcl[:, j] = mss_cnc_aer_in[c_idx, :, j]
            end

            # Set spectral underlying surface albedos
            albsfc_lcl = zeros(numrad_snw)
            albsfc_lcl[1:(nir_bnd_bgn-1)] .= albsfc[c_idx, IVIS]
            albsfc_lcl[nir_bnd_bgn:nir_bnd_end] .= albsfc[c_idx, INIR]

            # Error check snow grain size
            for i in snl_top_j:snl_btm_j
                if snw_rds_lcl[i] < SNW_RDS_MIN_TBL || snw_rds_lcl[i] > SNW_RDS_MAX_TBL
                    error("SNICAR ERROR: snow grain radius of $(snw_rds_lcl[i]) out of bounds at column $c_idx layer $i")
                end
            end

            # Flux weights
            if flg_slr_in == 1
                flx_wgt = copy(optics.flx_wgt_dir)
            else
                flx_wgt = copy(optics.flx_wgt_dif)
            end

            exp_min = exp(-argmax_val)

            # Local snow optical properties
            ss_alb_snw_lcl = zeros(nlevsno)
            asm_prm_snw_lcl = zeros(nlevsno)
            ext_cff_mss_snw_lcl = zeros(nlevsno)
            ss_alb_aer_lcl = zeros(SNO_NBR_AER)
            asm_prm_aer_lcl = zeros(SNO_NBR_AER)
            ext_cff_mss_aer_lcl = zeros(SNO_NBR_AER)

            # Snow shape arrays
            sno_shp = fill(snicar_snw_shape, nlevsno)
            sno_fs = zeros(nlevsno)
            sno_AR = zeros(nlevsno)

            # Local output arrays for all spectral bands
            albout_lcl = zeros(numrad_snw)
            flx_abs_lcl = zeros(nlevsno + 1, numrad_snw)  # layers + ground

            # Adding-doubling variables (snow layers + 1 interface)
            nlyr_itf = nlevsno + 1  # number of interfaces

            mu_not = 0.0

            # Loop over spectral bands
            for bnd_idx in 1:numrad_snw
                mu_not = max(coszen[c_idx], cp01)

                # Restore aerosol concentrations that may have been zeroed
                if flg_nosnl == 0
                    for j in 1:SNO_NBR_AER
                        mss_cnc_aer_lcl[:, j] = mss_cnc_aer_in[c_idx, :, j]
                    end
                else
                    fill!(mss_cnc_aer_lcl, 0.0)
                end

                # Set direct or diffuse incident irradiance to 1
                if flg_slr_in == 1
                    flx_slrd_lcl = 1.0 / (mu_not * pi_val)
                    flx_slri_lcl = 0.0
                else
                    flx_slrd_lcl = 0.0
                    flx_slri_lcl = 1.0
                end

                # Zero aerosols for highly absorptive bands
                if numrad_snw == DEFAULT_NUMBER_BANDS && (bnd_idx == HIGHEST_DEFAULT_BAND || bnd_idx == SEC_HIGHEST_DEFAULT_BAND)
                    fill!(mss_cnc_aer_lcl, 0.0)
                end
                if numrad_snw == HIGH_NUMBER_BANDS && bnd_idx > 100
                    fill!(mss_cnc_aer_lcl, 0.0)
                end

                # ---- Snow & aerosol optics ----

                # Snow optical properties from lookup table
                if flg_slr_in == 1
                    for i in snl_top_j:snl_btm_j
                        rds_idx = snw_rds_lcl[i] - SNW_RDS_MIN_TBL + 1
                        ss_alb_snw_lcl[i] = optics.ss_alb_snw_drc[rds_idx, bnd_idx]
                        ext_cff_mss_snw_lcl[i] = optics.ext_cff_mss_snw_drc[rds_idx, bnd_idx]
                        if sno_shp[i] == "sphere"
                            asm_prm_snw_lcl[i] = optics.asm_prm_snw_drc[rds_idx, bnd_idx]
                        end
                    end
                else
                    for i in snl_top_j:snl_btm_j
                        rds_idx = snw_rds_lcl[i] - SNW_RDS_MIN_TBL + 1
                        ss_alb_snw_lcl[i] = optics.ss_alb_snw_dfs[rds_idx, bnd_idx]
                        ext_cff_mss_snw_lcl[i] = optics.ext_cff_mss_snw_dfs[rds_idx, bnd_idx]
                        if sno_shp[i] == "sphere"
                            asm_prm_snw_lcl[i] = optics.asm_prm_snw_dfs[rds_idx, bnd_idx]
                        end
                    end
                end

                # Nonspherical snow: shape-dependent asymmetry factors
                g_ice_Cg_tmp = zeros(SNICAR_SEVEN_BANDS)
                gg_ice_F07_tmp = zeros(SNICAR_SEVEN_BANDS)

                for i in snl_top_j:snl_btm_j
                    if sno_shp[i] == "spheroid"
                        diam_ice = 2.0 * snw_rds_lcl[i]
                        fs_sphd = sno_fs[i] == 0.0 ? FS_SPHD_DEFAULT : sno_fs[i]
                        fs_hex = FS_HEX_REF
                        AR_tmp = sno_AR[i] == 0.0 ? AR_TMP_DEFAULT_1 : sno_AR[i]
                        for igb in 1:SNICAR_SEVEN_BANDS
                            g_ice_Cg_tmp[igb] = SNICAR_G_B0[igb] * ((fs_sphd/fs_hex)^SNICAR_G_B1[igb]) * (diam_ice^SNICAR_G_B2[igb])
                            gg_ice_F07_tmp[igb] = SNICAR_G_F07_C0[igb] + SNICAR_G_F07_C1[igb] * AR_tmp + SNICAR_G_F07_C2[igb] * (AR_tmp * AR_tmp)
                        end

                    elseif sno_shp[i] == "hexagonal_plate"
                        diam_ice = 2.0 * snw_rds_lcl[i]
                        fs_hex0 = sno_fs[i] == 0.0 ? FS_HEX_REF : sno_fs[i]
                        fs_hex = FS_HEX_REF
                        AR_tmp = sno_AR[i] == 0.0 ? AR_TMP_DEFAULT_2 : sno_AR[i]
                        for igb in 1:SNICAR_SEVEN_BANDS
                            g_ice_Cg_tmp[igb] = SNICAR_G_B0[igb] * ((fs_hex0/fs_hex)^SNICAR_G_B1[igb]) * (diam_ice^SNICAR_G_B2[igb])
                            gg_ice_F07_tmp[igb] = SNICAR_G_F07_P0[igb] + SNICAR_G_F07_P1[igb] * log(AR_tmp) + SNICAR_G_F07_P2[igb] * (log(AR_tmp) * log(AR_tmp))
                        end

                    elseif sno_shp[i] == "koch_snowflake"
                        diam_ice = 2.0 * snw_rds_lcl[i] / 0.544
                        fs_koch = sno_fs[i] == 0.0 ? FS_KOCH_DEFAULT : sno_fs[i]
                        fs_hex = FS_HEX_REF
                        AR_tmp = sno_AR[i] == 0.0 ? AR_TMP_DEFAULT_2 : sno_AR[i]
                        for igb in 1:SNICAR_SEVEN_BANDS
                            g_ice_Cg_tmp[igb] = SNICAR_G_B0[igb] * ((fs_koch/fs_hex)^SNICAR_G_B1[igb]) * (diam_ice^SNICAR_G_B2[igb])
                            gg_ice_F07_tmp[igb] = SNICAR_G_F07_P0[igb] + SNICAR_G_F07_P1[igb] * log(AR_tmp) + SNICAR_G_F07_P2[igb] * (log(AR_tmp) * log(AR_tmp))
                        end

                    elseif sno_shp[i] == "sphere"
                        # DO NOTHING
                    else
                        error("SNICAR ERROR: unknown sno_shp: $(sno_shp[i])")
                    end

                    # Compute nonspherical snow asymmetry factor
                    if sno_shp[i] != "sphere"
                        g_Cg_intp = piecewise_linear_interp1d(g_wvl_ct, g_ice_Cg_tmp, wvl_ct[bnd_idx])
                        gg_F07_intp = piecewise_linear_interp1d(g_wvl_ct, gg_ice_F07_tmp, wvl_ct[bnd_idx])
                        g_ice_F07 = gg_F07_intp + 0.5 * (1.0 - gg_F07_intp) / ss_alb_snw_lcl[i]
                        asm_prm_snw_lcl[i] = g_ice_F07 * g_Cg_intp
                    end

                    asm_prm_snw_lcl[i] = min(0.99, asm_prm_snw_lcl[i])
                end

                # Aerosol optical properties from lookup table
                ss_alb_aer_lcl[1] = optics.ss_alb_bc_hphil[bnd_idx]
                asm_prm_aer_lcl[1] = optics.asm_prm_bc_hphil[bnd_idx]
                ext_cff_mss_aer_lcl[1] = optics.ext_cff_mss_bc_hphil[bnd_idx]
                ss_alb_aer_lcl[2] = optics.ss_alb_bc_hphob[bnd_idx]
                asm_prm_aer_lcl[2] = optics.asm_prm_bc_hphob[bnd_idx]
                ext_cff_mss_aer_lcl[2] = optics.ext_cff_mss_bc_hphob[bnd_idx]
                ss_alb_aer_lcl[3] = optics.ss_alb_oc_hphil[bnd_idx]
                asm_prm_aer_lcl[3] = optics.asm_prm_oc_hphil[bnd_idx]
                ext_cff_mss_aer_lcl[3] = optics.ext_cff_mss_oc_hphil[bnd_idx]
                ss_alb_aer_lcl[4] = optics.ss_alb_oc_hphob[bnd_idx]
                asm_prm_aer_lcl[4] = optics.asm_prm_oc_hphob[bnd_idx]
                ext_cff_mss_aer_lcl[4] = optics.ext_cff_mss_oc_hphob[bnd_idx]
                ss_alb_aer_lcl[5] = optics.ss_alb_dst1[bnd_idx]
                asm_prm_aer_lcl[5] = optics.asm_prm_dst1[bnd_idx]
                ext_cff_mss_aer_lcl[5] = optics.ext_cff_mss_dst1[bnd_idx]
                ss_alb_aer_lcl[6] = optics.ss_alb_dst2[bnd_idx]
                asm_prm_aer_lcl[6] = optics.asm_prm_dst2[bnd_idx]
                ext_cff_mss_aer_lcl[6] = optics.ext_cff_mss_dst2[bnd_idx]
                ss_alb_aer_lcl[7] = optics.ss_alb_dst3[bnd_idx]
                asm_prm_aer_lcl[7] = optics.asm_prm_dst3[bnd_idx]
                ext_cff_mss_aer_lcl[7] = optics.ext_cff_mss_dst3[bnd_idx]
                ss_alb_aer_lcl[8] = optics.ss_alb_dst4[bnd_idx]
                asm_prm_aer_lcl[8] = optics.asm_prm_dst4[bnd_idx]
                ext_cff_mss_aer_lcl[8] = optics.ext_cff_mss_dst4[bnd_idx]

                # Layer mass, optical depths, weighted Mie properties
                wvl_doint = wvl_ct[bnd_idx]
                L_snw = zeros(nlevsno)
                tau_snw = zeros(nlevsno)
                L_aer = zeros(nlevsno, SNO_NBR_AER)
                tau_aer = zeros(nlevsno, SNO_NBR_AER)
                tau = zeros(nlevsno)
                omega = zeros(nlevsno)
                g_arr = zeros(nlevsno)

                enh_omg_bcint_tmp = zeros(SNICAR_SIXTEEN_BANDS)
                enh_omg_bcint_tmp2 = zeros(SNICAR_SIXTEEN_BANDS)
                enh_omg_dstint_tmp = zeros(SNICAR_SIZE_BINS)
                enh_omg_dstint_tmp2 = zeros(SNICAR_SIZE_BINS)

                for i in snl_top_j:snl_btm_j
                    # BC/dust-snow internal mixing for wavelength <= 1.2um
                    if wvl_doint <= 1.2

                        # BC-snow internal mixing
                        if snicar_snobc_intmix && mss_cnc_aer_lcl[i, 1] > 0.0
                            for ibb in 1:SNICAR_SIXTEEN_BANDS
                                enh_omg_bcint_tmp[ibb] = BCINT_D0[ibb] *
                                    ((mss_cnc_aer_lcl[i, 1] * KG_TO_UG * DEN_BC / DEN_BC_TARGET + BCINT_D2[ibb])^BCINT_D1[ibb])
                                if ibb < 3
                                    bcint_m_tmp = BCINT_M[1]; bcint_n_tmp = BCINT_N[1]
                                elseif ibb >= 3 && ibb <= 11
                                    bcint_m_tmp = BCINT_M[2]; bcint_n_tmp = BCINT_N[2]
                                else
                                    bcint_m_tmp = BCINT_M[3]; bcint_n_tmp = BCINT_N[3]
                                end
                                bcint_dd = (RE_BC / RADIUS_2_BC)^bcint_m_tmp
                                bcint_dd2 = (RADIUS_1_BC / RADIUS_2_BC)^bcint_m_tmp
                                bcint_f = (RE_BC / RADIUS_1_BC)^bcint_n_tmp
                                enh_omg_bcint_tmp2[ibb] = log10(max(1.0, bcint_dd * ((enh_omg_bcint_tmp[ibb] / bcint_dd2)^bcint_f)))
                            end
                            enh_omg_bcint_intp = piecewise_linear_interp1d(bcint_wvl_ct, enh_omg_bcint_tmp2, wvl_doint)
                            enh_omg_bcint_intp2 = 10.0^enh_omg_bcint_intp
                            enh_omg_bcint_intp2 = min(ENH_OMG_MAX, max(enh_omg_bcint_intp2, 1.0))
                            ss_alb_snw_lcl[i] = 1.0 - (1.0 - ss_alb_snw_lcl[i]) * enh_omg_bcint_intp2
                            ss_alb_snw_lcl[i] = max(0.5, min(ss_alb_snw_lcl[i], 1.0))
                            ss_alb_aer_lcl[1] = 0.0
                            asm_prm_aer_lcl[1] = 0.0
                            ext_cff_mss_aer_lcl[1] = 0.0
                        end

                        # Dust-snow internal mixing
                        tot_dst_snw_conc = (mss_cnc_aer_lcl[i, 5] + mss_cnc_aer_lcl[i, 6] +
                                            mss_cnc_aer_lcl[i, 7] + mss_cnc_aer_lcl[i, 8]) * KG_KG_TO_PPM
                        if snicar_snodst_intmix && tot_dst_snw_conc > 0.0
                            for idb in 1:SNICAR_SIZE_BINS
                                enh_omg_dstint_tmp[idb] = DSTINT_A1[idb] + DSTINT_A2[idb] * (tot_dst_snw_conc^DSTINT_A3[idb])
                                enh_omg_dstint_tmp2[idb] = log10(max(enh_omg_dstint_tmp[idb], 1.0))
                            end
                            enh_omg_dstint_intp = piecewise_linear_interp1d(dstint_wvl_ct, enh_omg_dstint_tmp2, wvl_doint)
                            enh_omg_dstint_intp2 = 10.0^enh_omg_dstint_intp
                            enh_omg_dstint_intp2 = min(ENH_OMG_MAX, max(enh_omg_dstint_intp2, 1.0))
                            ss_alb_snw_lcl[i] = 1.0 - (1.0 - ss_alb_snw_lcl[i]) * enh_omg_dstint_intp2
                            ss_alb_snw_lcl[i] = max(0.5, min(ss_alb_snw_lcl[i], 1.0))
                            ss_alb_aer_lcl[5:8] .= 0.0
                            asm_prm_aer_lcl[5:8] .= 0.0
                            ext_cff_mss_aer_lcl[5:8] .= 0.0
                        end
                    end # BC/dust internal mixing

                    L_snw[i] = h2osno_ice_lcl[i] + h2osno_liq_lcl[i]
                    tau_snw[i] = L_snw[i] * ext_cff_mss_snw_lcl[i]

                    for j in 1:SNO_NBR_AER
                        L_aer[i, j] = L_snw[i] * mss_cnc_aer_lcl[i, j]
                        tau_aer[i, j] = L_aer[i, j] * ext_cff_mss_aer_lcl[j]
                    end

                    tau_sum = 0.0
                    omega_sum = 0.0
                    g_sum_val = 0.0
                    for j in 1:SNO_NBR_AER
                        tau_sum += tau_aer[i, j]
                        omega_sum += tau_aer[i, j] * ss_alb_aer_lcl[j]
                        g_sum_val += tau_aer[i, j] * ss_alb_aer_lcl[j] * asm_prm_aer_lcl[j]
                    end

                    tau[i] = tau_sum + tau_snw[i]
                    omega[i] = (1.0 / tau[i]) * (omega_sum + ss_alb_snw_lcl[i] * tau_snw[i])
                    g_arr[i] = (1.0 / (tau[i] * omega[i])) * (g_sum_val + asm_prm_snw_lcl[i] * ss_alb_snw_lcl[i] * tau_snw[i])
                end

                # Delta transformations
                tau_star = zeros(nlevsno)
                omega_star = zeros(nlevsno)
                g_star = zeros(nlevsno)

                if DELTA == 1
                    for i in snl_top_j:snl_btm_j
                        g_star[i] = g_arr[i] / (1.0 + g_arr[i])
                        omega_star[i] = (1.0 - g_arr[i]^2) * omega[i] / (1.0 - omega[i] * g_arr[i]^2)
                        tau_star[i] = (1.0 - omega[i] * g_arr[i]^2) * tau[i]
                    end
                else
                    for i in snl_top_j:snl_btm_j
                        g_star[i] = g_arr[i]
                        omega_star[i] = omega[i]
                        tau_star[i] = tau[i]
                    end
                end

                # ---- Adding-Doubling RT solver ----
                snl_btm_itf = snl_btm_j + 1  # bottom interface (1-based)

                # Interface arrays: indices snl_top_j to snl_btm_itf
                trndir = zeros(nlyr_itf)
                trntdr = zeros(nlyr_itf)
                trndif = zeros(nlyr_itf)
                rupdir = zeros(nlyr_itf)
                rupdif = zeros(nlyr_itf)
                rdndif = zeros(nlyr_itf)
                dfdir = zeros(nlyr_itf)
                dfdif = zeros(nlyr_itf)
                dftmp = zeros(nlyr_itf)

                # Layer arrays
                rdir = zeros(nlevsno)
                rdif_a = zeros(nlevsno)
                rdif_b = zeros(nlevsno)
                tdir_arr = zeros(nlevsno)
                tdif_a = zeros(nlevsno)
                tdif_b = zeros(nlevsno)
                trnlay = zeros(nlevsno)

                # Initialize top interface
                trndir[snl_top_j] = c1
                trntdr[snl_top_j] = c1
                trndif[snl_top_j] = c1

                # Main level loop
                for i in snl_top_j:snl_btm_j
                    if trntdr[i] > trmin
                        ts = tau_star[i]
                        ws = omega_star[i]
                        gs = g_star[i]

                        lm = sqrt(c3 * (c1 - ws) * (c1 - ws * gs))
                        ue = c1p5 * (c1 - ws * gs) / lm
                        extins = max(exp_min, exp(-lm * ts))
                        ne_val = ((ue + c1)^2 / extins) - ((ue - c1)^2 * extins)

                        rdif_a[i] = (ue^2 - c1) * (1.0/extins - extins) / ne_val
                        tdif_a[i] = c4 * ue / ne_val

                        trnlay[i] = max(exp_min, exp(-ts / mu_not))

                        alp = cp75 * ws * mu_not * ((c1 + gs * (c1 - ws)) / (c1 - lm^2 * mu_not^2))
                        gam = cp5 * ws * ((c1 + c3 * gs * (c1 - ws) * mu_not^2) / (c1 - lm^2 * mu_not^2))
                        apg = alp + gam
                        amg = alp - gam
                        rdir[i] = apg * rdif_a[i] + amg * (tdif_a[i] * trnlay[i] - c1)
                        tdir_arr[i] = apg * tdif_a[i] + (amg * rdif_a[i] - apg + c1) * trnlay[i]

                        # Recalculate rdif, tdif using Gaussian integration
                        R1 = rdif_a[i]
                        T1 = tdif_a[i]
                        swt = c0; smr = c0; smt = c0
                        for ng in 1:SNICAR_NGMAX
                            mu = SNICAR_DIFGAUSPT[ng]
                            gwt = SNICAR_DIFGAUSWT[ng]
                            swt += mu * gwt
                            trn = max(exp_min, exp(-ts / mu))
                            alp = cp75 * ws * mu * ((c1 + gs * (c1 - ws)) / (c1 - lm^2 * mu^2))
                            gam = cp5 * ws * ((c1 + c3 * gs * (c1 - ws) * mu^2) / (c1 - lm^2 * mu^2))
                            apg = alp + gam
                            amg = alp - gam
                            rdr = apg * R1 + amg * T1 * trn - amg
                            tdr = apg * T1 + amg * R1 * trn - apg * trn + trn
                            smr += mu * rdr * gwt
                            smt += mu * tdr * gwt
                        end
                        rdif_a[i] = smr / swt
                        tdif_a[i] = smt / swt

                        rdif_b[i] = rdif_a[i]
                        tdif_b[i] = tdif_a[i]
                    end

                    # Interface calculations
                    trndir[i+1] = trndir[i] * trnlay[i]
                    refkm1 = c1 / (c1 - rdndif[i] * rdif_a[i])
                    tdrrdir = trndir[i] * rdir[i]
                    tdndif = trntdr[i] - trndir[i]
                    trntdr[i+1] = trndir[i] * tdir_arr[i] +
                                  (tdndif + tdrrdir * rdndif[i]) * refkm1 * tdif_a[i]
                    rdndif[i+1] = rdif_b[i] + tdif_b[i] * rdndif[i] * refkm1 * tdif_a[i]
                    trndif[i+1] = trndif[i] * refkm1 * tdif_a[i]
                end

                # Bottom boundary: set underlying ground albedo
                rupdir[snl_btm_itf] = albsfc[c_idx, INIR]
                rupdif[snl_btm_itf] = albsfc[c_idx, INIR]
                if bnd_idx < nir_bnd_bgn
                    rupdir[snl_btm_itf] = albsfc[c_idx, IVIS]
                    rupdif[snl_btm_itf] = albsfc[c_idx, IVIS]
                end

                # Upward sweep
                for i in snl_btm_j:-1:snl_top_j
                    refkp1 = c1 / (c1 - rdif_b[i] * rupdif[i+1])
                    rupdir[i] = rdir[i] +
                                (trnlay[i] * rupdir[i+1] +
                                 (tdir_arr[i] - trnlay[i]) * rupdif[i+1]) * refkp1 * tdif_b[i]
                    rupdif[i] = rdif_a[i] + tdif_a[i] * rupdif[i+1] * refkp1 * tdif_b[i]
                end

                # Net flux at each interface
                for i in snl_top_j:snl_btm_itf
                    refk = c1 / (c1 - rdndif[i] * rupdif[i])
                    dfdir[i] = trndir[i] +
                               (trntdr[i] - trndir[i]) * (c1 - rupdif[i]) * refk -
                               trndir[i] * rupdir[i] * (c1 - rdndif[i]) * refk
                    if dfdir[i] < SNICAR_PUNY
                        dfdir[i] = c0
                    end
                    dfdif[i] = trndif[i] * (c1 - rupdif[i]) * refk
                    if dfdif[i] < SNICAR_PUNY
                        dfdif[i] = c0
                    end
                end

                # Select direct/diffuse results
                if flg_slr_in == 1
                    albedo = rupdir[snl_top_j]
                    dftmp .= dfdir
                    refk = c1 / (c1 - rdndif[snl_top_j] * rupdif[snl_top_j])
                    F_sfc_pls = (trndir[snl_top_j] * rupdir[snl_top_j] +
                                (trntdr[snl_top_j] - trndir[snl_top_j]) * rupdif[snl_top_j]) * refk
                else
                    albedo = rupdif[snl_top_j]
                    dftmp .= dfdif
                    refk = c1 / (c1 - rdndif[snl_top_j] * rupdif[snl_top_j])
                    F_sfc_pls = trndif[snl_top_j] * rupdif[snl_top_j] * refk
                end

                # Absorbed flux in each layer
                F_abs = zeros(nlevsno)
                F_abs_sum = 0.0
                for i in snl_top_j:snl_btm_j
                    F_abs[i] = dftmp[i] - dftmp[i+1]
                    flx_abs_lcl[i, bnd_idx] = F_abs[i]
                    F_abs_sum += F_abs[i]
                end

                # Absorbed flux by underlying ground
                F_btm_net = dftmp[snl_btm_itf]
                flx_abs_lcl[snl_btm_itf, bnd_idx] = F_btm_net

                if flg_nosnl == 1
                    flx_abs_lcl[snl_btm_j, bnd_idx] = F_abs[snl_btm_j]
                    flx_abs_lcl[snl_btm_itf, bnd_idx] = F_btm_net
                end

                # Underflow check
                for i in snl_top_j:snl_btm_itf
                    flx_abs_lcl[i, bnd_idx] = max(0.0, flx_abs_lcl[i, bnd_idx])
                end

                albout_lcl[bnd_idx] = albedo
            end  # loop over spectral bands

            # Weight output albedo into VIS/NIR bands
            if numrad_snw == DEFAULT_NUMBER_BANDS
                albout[c_idx, IVIS] = albout_lcl[IVIS]
            elseif numrad_snw == HIGH_NUMBER_BANDS
                flx_sum = 0.0
                for bnd_idx in 1:(nir_bnd_bgn-1)
                    flx_sum += flx_wgt[bnd_idx] * albout_lcl[bnd_idx]
                end
                albout[c_idx, IVIS] = flx_sum / sum(flx_wgt[1:(nir_bnd_bgn-1)])
            end

            # NIR band average
            flx_sum = 0.0
            for bnd_idx in nir_bnd_bgn:nir_bnd_end
                flx_sum += flx_wgt[bnd_idx] * albout_lcl[bnd_idx]
            end
            albout[c_idx, INIR] = flx_sum / sum(flx_wgt[nir_bnd_bgn:nir_bnd_end])

            # Weight output flx_abs into VIS/NIR bands
            if numrad_snw == DEFAULT_NUMBER_BANDS
                for i in snl_top_j:(nlevsno+1)
                    flx_abs[c_idx, i, IVIS] = flx_abs_lcl[i, 1]
                end
            elseif numrad_snw == HIGH_NUMBER_BANDS
                for i in snl_top_j:(nlevsno+1)
                    flx_sum = 0.0
                    for bnd_idx in 1:(nir_bnd_bgn-1)
                        flx_sum += flx_wgt[bnd_idx] * flx_abs_lcl[i, bnd_idx]
                    end
                    flx_abs[c_idx, i, IVIS] = flx_sum / sum(flx_wgt[1:(nir_bnd_bgn-1)])
                end
            end

            # NIR flux absorption
            for i in snl_top_j:(nlevsno+1)
                flx_sum = 0.0
                for bnd_idx in nir_bnd_bgn:nir_bnd_end
                    flx_sum += flx_wgt[bnd_idx] * flx_abs_lcl[i, bnd_idx]
                end
                flx_abs[c_idx, i, INIR] = flx_sum / sum(flx_wgt[nir_bnd_bgn:nir_bnd_end])
            end

            # High solar zenith angle adjustment for NIR direct albedo
            if mu_not < MU_75 && flg_slr_in == 1
                sza_c1_val = SZA_A0 + SZA_A1 * mu_not + SZA_A2 * mu_not^2
                sza_c0_val = SZA_B0 + SZA_B1 * mu_not + SZA_B2 * mu_not^2
                sza_factor = sza_c1_val * (log10(snw_rds_lcl[snl_top_j] * 1.0) - 6.0) + sza_c0_val
                flx_sza_adjust = albout[c_idx, INIR] * (sza_factor - 1.0) * sum(flx_wgt[nir_bnd_bgn:nir_bnd_end])
                albout[c_idx, INIR] *= sza_factor
                flx_abs[c_idx, snl_top_j, INIR] -= flx_sza_adjust
            end

        elseif coszen[c_idx] > 0.0 && h2osno_lcl < MIN_SNW && h2osno_lcl > 0.0
            # Snow < minimum but > 0 and there is sun: set to surface albedo
            albout[c_idx, IVIS] = albsfc[c_idx, IVIS]
            albout[c_idx, INIR] = albsfc[c_idx, INIR]

        else
            # Zero snow or no sun
            albout[c_idx, IVIS] = 0.0
            albout[c_idx, INIR] = 0.0
        end
    end  # loop over columns

    return nothing
end

# --------------------------------------------------------------------------
# snowage_grain!
# --------------------------------------------------------------------------

"""
    snowage_grain!(snl, dz, frac_sno, h2osoi_liq, h2osoi_ice,
                   t_soisno, t_grnd, qflx_snow_grnd_col, qflx_snofrz_lyr,
                   h2osno_no_layers, forc_t,
                   snw_rds, snw_rds_top, sno_liq_top,
                   snot_top, dTdz_top,
                   nlevsno, dtime;
                   mask_snowc=nothing, mask_nosnowc=nothing,
                   params=snicar_params, aging=snicar_aging)

Updates the snow effective grain size (radius).
Contributions to grain size evolution are from:
1. vapor redistribution (dry snow)
2. liquid water redistribution (wet snow)
3. re-freezing of liquid water

Ported from `SnowAge_grain` in `SnowSnicarMod.F90`.

# Arguments
- `snl::Vector{Int}`: negative number of snow layers (col)
- `dz::Matrix{Float64}`: layer thickness (col,lyr) [m]
- `frac_sno::Vector{Float64}`: fraction of ground covered by snow (col)
- `h2osoi_liq::Matrix{Float64}`: liquid water content (col,lyr) [kg/m2]
- `h2osoi_ice::Matrix{Float64}`: ice content (col,lyr) [kg/m2]
- `t_soisno::Matrix{Float64}`: soil and snow temperature (col,lyr) [K]
- `t_grnd::Vector{Float64}`: ground temperature (col) [K]
- `qflx_snow_grnd_col::Vector{Float64}`: snow on ground after interception [kg/m2/s]
- `qflx_snofrz_lyr::Matrix{Float64}`: snow freezing rate (col,lyr) [kg/m2/s]
- `h2osno_no_layers::Vector{Float64}`: snow not resolved into layers [mm H2O]
- `forc_t::Vector{Float64}`: atmospheric temperature (col) [K]
- `snw_rds::Matrix{Float64}`: (in/out) effective grain radius (col,lyr) [microns]
- `snw_rds_top::Vector{Float64}`: (out) effective grain radius top layer [microns]
- `sno_liq_top::Vector{Float64}`: (out) liquid water fraction in top layer [frc]
- `snot_top::Vector{Float64}`: (out) temperature in top snow layer [K]
- `dTdz_top::Vector{Float64}`: (out) temperature gradient in top layer [K/m]
- `nlevsno::Int`: maximum number of snow layers
- `dtime::Float64`: time step [sec]
"""
function snowage_grain!(snl::Vector{Int},
                        dz::Matrix{Float64},
                        frac_sno::Vector{Float64},
                        h2osoi_liq::Matrix{Float64},
                        h2osoi_ice::Matrix{Float64},
                        t_soisno::Matrix{Float64},
                        t_grnd::Vector{Float64},
                        qflx_snow_grnd_col::Vector{Float64},
                        qflx_snofrz_lyr::Matrix{Float64},
                        h2osno_no_layers::Vector{Float64},
                        forc_t::Vector{Float64},
                        snw_rds::Matrix{Float64},
                        snw_rds_top::Vector{Float64},
                        sno_liq_top::Vector{Float64},
                        snot_top::Vector{Float64},
                        dTdz_top::Vector{Float64},
                        nlevsno::Int,
                        dtime::Float64;
                        mask_snowc::Union{BitVector,Nothing}=nothing,
                        mask_nosnowc::Union{BitVector,Nothing}=nothing,
                        params::SnicarParams=snicar_params,
                        aging::SnicarAgingData=snicar_aging)

    joff = nlevsno  # offset: Fortran index i maps to Julia index i + joff
    ncols = length(snl)

    # Loop over columns with snow layers
    for c_idx in 1:ncols
        if mask_snowc !== nothing && !mask_snowc[c_idx]
            continue
        end
        if snl[c_idx] >= 0
            continue
        end

        snl_btm = 0
        snl_top = snl[c_idx] + 1

        snl_btm_j = snl_btm + joff
        snl_top_j = snl_top + joff

        # Column average layer thickness
        cdz = zeros(nlevsno)
        for i in snl_top_j:snl_btm_j
            cdz[i] = frac_sno[c_idx] * dz[c_idx, i]
        end

        for i in snl_top_j:snl_btm_j
            # 1. DRY SNOW AGING
            h2osno_lyr = h2osoi_liq[c_idx, i] + h2osoi_ice[c_idx, i]

            # Temperature gradient
            if i == snl_top_j
                t_snotop = t_soisno[c_idx, snl_top_j]
                t_snobtm = (t_soisno[c_idx, i+1] * dz[c_idx, i] +
                            t_soisno[c_idx, i] * dz[c_idx, i+1]) /
                           (dz[c_idx, i] + dz[c_idx, i+1])
            else
                t_snotop = (t_soisno[c_idx, i-1] * dz[c_idx, i] +
                            t_soisno[c_idx, i] * dz[c_idx, i-1]) /
                           (dz[c_idx, i] + dz[c_idx, i-1])
                t_snobtm = (t_soisno[c_idx, i+1] * dz[c_idx, i] +
                            t_soisno[c_idx, i] * dz[c_idx, i+1]) /
                           (dz[c_idx, i] + dz[c_idx, i+1])
            end

            dTdz_val = abs((t_snotop - t_snobtm) / cdz[i])

            # Snow density
            rhos = (h2osoi_liq[c_idx, i] + h2osoi_ice[c_idx, i]) / cdz[i]
            rhos = max(50.0, rhos)

            # Best-fit table indices
            T_idx = round(Int, (t_soisno[c_idx, i] - 223) / 5) + 1
            Tgrd_idx = round(Int, dTdz_val / 10) + 1
            rhos_idx = round(Int, (rhos - 50) / 50) + 1

            T_idx = clamp(T_idx, IDX_T_MIN, IDX_T_MAX)
            Tgrd_idx = clamp(Tgrd_idx, IDX_TGRD_MIN, IDX_TGRD_MAX)
            rhos_idx = clamp(rhos_idx, IDX_RHOS_MIN, IDX_RHOS_MAX)

            bst_tau = aging.snowage_tau[rhos_idx, Tgrd_idx, T_idx]
            bst_kappa = aging.snowage_kappa[rhos_idx, Tgrd_idx, T_idx]
            bst_drdt0 = aging.snowage_drdt0[rhos_idx, Tgrd_idx, T_idx]

            # Boundary check
            snw_rds[c_idx, i] = max(snw_rds[c_idx, i], params.snw_rds_min)

            # Change in effective radius
            dr_fresh = snw_rds[c_idx, i] - params.snw_rds_min
            dr = (bst_drdt0 * (bst_tau / (dr_fresh + bst_tau))^(1.0 / bst_kappa)) * (dtime / SECSPHR)

            # 2. WET SNOW AGING
            frc_liq = min(0.1, h2osoi_liq[c_idx, i] / (h2osoi_liq[c_idx, i] + h2osoi_ice[c_idx, i]))
            dr_wet = 1.0e18 * (dtime * (params.C2_liq_Brun89 * frc_liq^3) /
                     (4.0 * π * snw_rds[c_idx, i] * snw_rds[c_idx, i]))

            dr += dr_wet

            # 3. SNOWAGE SCALING
            dr *= params.xdrdt

            # 4. INCREMENT EFFECTIVE RADIUS
            newsnow = max(0.0, qflx_snow_grnd_col[c_idx] * dtime)
            refrzsnow = max(0.0, qflx_snofrz_lyr[c_idx, i] * dtime)

            frc_refrz = refrzsnow / h2osno_lyr

            if i == snl_top_j
                frc_newsnow = newsnow / h2osno_lyr
            else
                frc_newsnow = 0.0
            end

            if (frc_refrz + frc_newsnow) > 1.0
                frc_refrz = frc_refrz / (frc_refrz + frc_newsnow)
                frc_newsnow = 1.0 - frc_refrz
                frc_oldsnow = 0.0
            else
                frc_oldsnow = 1.0 - frc_refrz - frc_newsnow
            end

            snw_rds_fresh = fresh_snow_radius(forc_t[c_idx]; params=params)

            snw_rds[c_idx, i] = (snw_rds[c_idx, i] + dr) * frc_oldsnow +
                                 snw_rds_fresh * frc_newsnow +
                                 params.snw_rds_refrz * frc_refrz

            # 5. BOUNDARY CHECK
            snw_rds[c_idx, i] = max(snw_rds[c_idx, i], params.snw_rds_min)
            snw_rds[c_idx, i] = min(snw_rds[c_idx, i], SNW_RDS_MAX)

            # Top layer variables for history
            if i == snl_top_j
                snot_top[c_idx] = t_soisno[c_idx, i]
                dTdz_top[c_idx] = dTdz_val
                snw_rds_top[c_idx] = snw_rds[c_idx, i]
                sno_liq_top[c_idx] = h2osoi_liq[c_idx, i] / (h2osoi_liq[c_idx, i] + h2osoi_ice[c_idx, i])
            end
        end
    end

    # Columns with no snow layers but some snow mass
    for c_idx in 1:ncols
        if mask_nosnowc !== nothing && !mask_nosnowc[c_idx]
            continue
        end
        if snl[c_idx] < 0
            continue
        end
        if h2osno_no_layers[c_idx] > 0.0
            snw_rds[c_idx, nlevsno] = params.snw_rds_min
        end
    end

    return nothing
end

# --------------------------------------------------------------------------
# snow_optics_init!
# --------------------------------------------------------------------------

"""
    snow_optics_init!(optics; numrad_snw=5, fsnowoptics="")

Initialize snow optics lookup tables. If `fsnowoptics` is provided, reads
Mie lookup tables from the NetCDF file. Otherwise allocates zero arrays.

Ported from `SnowOptics_init` in `SnowSnicarMod.F90`.
"""
function snow_optics_init!(optics::SnicarOpticsData;
                            numrad_snw::Int=varctl.snicar_numrad_snw,
                            fsnowoptics::String="")

    # Allocate arrays
    optics.ss_alb_snw_drc = zeros(IDX_MIE_SNW_MX, numrad_snw)
    optics.asm_prm_snw_drc = zeros(IDX_MIE_SNW_MX, numrad_snw)
    optics.ext_cff_mss_snw_drc = zeros(IDX_MIE_SNW_MX, numrad_snw)
    optics.ss_alb_snw_dfs = zeros(IDX_MIE_SNW_MX, numrad_snw)
    optics.asm_prm_snw_dfs = zeros(IDX_MIE_SNW_MX, numrad_snw)
    optics.ext_cff_mss_snw_dfs = zeros(IDX_MIE_SNW_MX, numrad_snw)

    optics.ss_alb_bc_hphil = zeros(numrad_snw)
    optics.asm_prm_bc_hphil = zeros(numrad_snw)
    optics.ext_cff_mss_bc_hphil = zeros(numrad_snw)
    optics.ss_alb_bc_hphob = zeros(numrad_snw)
    optics.asm_prm_bc_hphob = zeros(numrad_snw)
    optics.ext_cff_mss_bc_hphob = zeros(numrad_snw)
    optics.ss_alb_oc_hphil = zeros(numrad_snw)
    optics.asm_prm_oc_hphil = zeros(numrad_snw)
    optics.ext_cff_mss_oc_hphil = zeros(numrad_snw)
    optics.ss_alb_oc_hphob = zeros(numrad_snw)
    optics.asm_prm_oc_hphob = zeros(numrad_snw)
    optics.ext_cff_mss_oc_hphob = zeros(numrad_snw)

    optics.ss_alb_dst1 = zeros(numrad_snw)
    optics.asm_prm_dst1 = zeros(numrad_snw)
    optics.ext_cff_mss_dst1 = zeros(numrad_snw)
    optics.ss_alb_dst2 = zeros(numrad_snw)
    optics.asm_prm_dst2 = zeros(numrad_snw)
    optics.ext_cff_mss_dst2 = zeros(numrad_snw)
    optics.ss_alb_dst3 = zeros(numrad_snw)
    optics.asm_prm_dst3 = zeros(numrad_snw)
    optics.ext_cff_mss_dst3 = zeros(numrad_snw)
    optics.ss_alb_dst4 = zeros(numrad_snw)
    optics.asm_prm_dst4 = zeros(numrad_snw)
    optics.ext_cff_mss_dst4 = zeros(numrad_snw)

    optics.flx_wgt_dir = zeros(numrad_snw)
    optics.flx_wgt_dif = zeros(numrad_snw)

    # Read from NetCDF if file provided
    if !isempty(fsnowoptics) && isfile(fsnowoptics)
        _snow_optics_read_nc!(optics, fsnowoptics, numrad_snw)
    end

    return nothing
end

"""
    _snow_optics_read_nc!(optics, fsnowoptics, numrad_snw)

Read SNICAR optics from NetCDF file. Uses 5-band mid-latitude-winter
Saharan dust defaults matching CLM5.
"""
function _snow_optics_read_nc!(optics::SnicarOpticsData, fsnowoptics::String, numrad_snw::Int)
    # Solar spectrum and dust optics short codes (CLM5 defaults)
    ss = "mlw"   # mid_latitude_winter
    ds = "sah"   # sahara

    NCDataset(fsnowoptics, "r") do nc
        # Helper to read a variable, handling missing values
        function _read_var(varname)
            if haskey(nc, varname)
                data = Array(nc[varname])
                return replace(data, missing => 0.0)
            else
                return nothing
            end
        end

        # Ice optics (2D: idx_Mie_snw_mx × numrad_snw)
        for (field, prop, beam) in [
            (:ss_alb_snw_drc,      "ss_alb",      "dir"),
            (:asm_prm_snw_drc,     "asm_prm",     "dir"),
            (:ext_cff_mss_snw_drc, "ext_cff_mss", "dir"),
            (:ss_alb_snw_dfs,      "ss_alb",      "dif"),
            (:asm_prm_snw_dfs,     "asm_prm",     "dif"),
            (:ext_cff_mss_snw_dfs, "ext_cff_mss", "dif"),
        ]
            varname = "$(prop)_ice_pic16_$(beam)_$(ss)"
            data = _read_var(varname)
            if data !== nothing
                arr = getfield(optics, field)
                n1 = min(size(data, 1), size(arr, 1))
                n2 = min(size(data, 2), size(arr, 2))
                arr[1:n1, 1:n2] .= Float64.(data[1:n1, 1:n2])
            end
        end

        # BC optics (1D: numrad_snw) — both hphil and hphob read from bcphob
        for (field, species) in [
            (:ss_alb_bc_hphil, "bcphob"), (:ss_alb_bc_hphob, "bcphob"),
            (:asm_prm_bc_hphil, "bcphob"), (:asm_prm_bc_hphob, "bcphob"),
            (:ext_cff_mss_bc_hphil, "bcphob"), (:ext_cff_mss_bc_hphob, "bcphob"),
        ]
            prop = String(field)
            # Extract property type from field name
            if startswith(prop, "ss_alb")
                pname = "ss_alb"
            elseif startswith(prop, "asm_prm")
                pname = "asm_prm"
            else
                pname = "ext_cff_mss"
            end
            varname = "$(pname)_$(species)_dif_$(ss)"
            data = _read_var(varname)
            if data !== nothing
                arr = getfield(optics, field)
                n = min(length(data), length(arr))
                arr[1:n] .= Float64.(data[1:n])
            end
        end

        # OC optics — both hphil and hphob read from ocphob
        for (field, species) in [
            (:ss_alb_oc_hphil, "ocphob"), (:ss_alb_oc_hphob, "ocphob"),
            (:asm_prm_oc_hphil, "ocphob"), (:asm_prm_oc_hphob, "ocphob"),
            (:ext_cff_mss_oc_hphil, "ocphob"), (:ext_cff_mss_oc_hphob, "ocphob"),
        ]
            prop = String(field)
            if startswith(prop, "ss_alb")
                pname = "ss_alb"
            elseif startswith(prop, "asm_prm")
                pname = "asm_prm"
            else
                pname = "ext_cff_mss"
            end
            varname = "$(pname)_$(species)_dif_$(ss)"
            data = _read_var(varname)
            if data !== nothing
                arr = getfield(optics, field)
                n = min(length(data), length(arr))
                arr[1:n] .= Float64.(data[1:n])
            end
        end

        # Dust optics (4 species)
        for (didx, fsuffix) in [(1, :dst1), (2, :dst2), (3, :dst3), (4, :dst4)]
            dust_num = lpad(didx, 2, '0')
            for pname in ["ss_alb", "asm_prm", "ext_cff_mss"]
                varname = "$(pname)_dust$(dust_num)_$(ds)_dif_$(ss)"
                field = Symbol("$(pname)_$(fsuffix)")
                data = _read_var(varname)
                if data !== nothing
                    arr = getfield(optics, field)
                    n = min(length(data), length(arr))
                    arr[1:n] .= Float64.(data[1:n])
                end
            end
        end

        # Flux weights
        bnd = numrad_snw == 5 ? "5" : ""
        for (field, beam) in [(:flx_wgt_dir, "dir"), (:flx_wgt_dif, "dif")]
            varname = "flx_wgt_$(beam)$(bnd)_$(ss)"
            data = _read_var(varname)
            if data !== nothing
                arr = getfield(optics, field)
                n = min(length(data), length(arr))
                arr[1:n] .= Float64.(data[1:n])
            end
        end
    end

    return nothing
end

# --------------------------------------------------------------------------
# snowage_init!
# --------------------------------------------------------------------------

"""
    snowage_init!(aging; fsnowaging="")

Initialize snow aging lookup tables. If `fsnowaging` is provided, reads
tau/kappa/drdt0 3D arrays from the NetCDF file. Otherwise allocates zero arrays.

Ported from `SnowAge_init` in `SnowSnicarMod.F90`.
"""
function snowage_init!(aging::SnicarAgingData; fsnowaging::String="")
    aging.snowage_tau = zeros(IDX_RHOS_MAX, IDX_TGRD_MAX, IDX_T_MAX)
    aging.snowage_kappa = zeros(IDX_RHOS_MAX, IDX_TGRD_MAX, IDX_T_MAX)
    aging.snowage_drdt0 = zeros(IDX_RHOS_MAX, IDX_TGRD_MAX, IDX_T_MAX)

    if !isempty(fsnowaging) && isfile(fsnowaging)
        NCDataset(fsnowaging, "r") do nc
            for (field, varname) in [
                (:snowage_tau,   "tau"),
                (:snowage_kappa, "kappa"),
                (:snowage_drdt0, "drdsdt0"),  # note: NetCDF var has extra 's'
            ]
                if haskey(nc, varname)
                    data = Array(nc[varname])
                    data = replace(data, missing => 0.0)
                    arr = getfield(aging, field)
                    n1 = min(size(data, 1), size(arr, 1))
                    n2 = min(size(data, 2), size(arr, 2))
                    n3 = min(size(data, 3), size(arr, 3))
                    arr[1:n1, 1:n2, 1:n3] .= Float64.(data[1:n1, 1:n2, 1:n3])
                end
            end
        end
    end

    return nothing
end
