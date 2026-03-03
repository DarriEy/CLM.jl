# ==========================================================================
# Ported from: src/biogeochem/DustEmisBase.F90
# Dust emission base module
#
# Routines for dust mobilization and dry deposition.
# Simulates dust mobilization due to wind from the surface into the
# lowest atmospheric layer. On output flx_mss_vrt_dst(ndst) is the
# surface dust emission (kg/m**2/s) [+ = to atm].
# Calculates the turbulent component of dust dry deposition velocity
# through the lowest atmospheric layer. CAM will calculate the settling
# velocity through the whole atmospheric column.
#
# Public types:
#   DustEmisBaseData  — Dust emission state/flux data
#
# Public functions:
#   dust_emis_init!          — Initialize dust emission data
#   dust_emis_clean!         — Deallocate dust emission data
#   init_dust_vars!          — Initialize dust size distributions and factors
#   dust_dry_dep!            — Turbulent dry deposition for dust
#   check_dust_emis_is_valid — Validate dust emission consistency
#   get_dust_patch_vars      — Get patch-level dust variables
#   get_dust_const_vars      — Get constant dust variables
# ==========================================================================

# ---------------------------------------------------------------------------
# Module-level constants
# ---------------------------------------------------------------------------

const DNS_AER = 2.5e3           # [kg m-3] Aerosol density
const BOLTZ = 1.38065e-23       # [J/K/molecule] Boltzmann constant (SHR_CONST_BOLTZ)

# ---------------------------------------------------------------------------
# Local erf approximation (Abramowitz & Stegun 7.1.26, ~1.5e-7 accuracy)
# Used in init_dust_vars! for overlap factor computation.
# Avoids adding SpecialFunctions.jl as a dependency.
# ---------------------------------------------------------------------------

function _dust_erf(x::Float64)
    a1 =  0.254829592
    a2 = -0.284496736
    a3 =  1.421413741
    a4 = -1.453152027
    a5 =  1.061405429
    p  =  0.3275911
    sgn = x >= 0.0 ? 1.0 : -1.0
    ax = abs(x)
    t = 1.0 / (1.0 + p * ax)
    y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-ax * ax)
    return sgn * y
end

# ---------------------------------------------------------------------------
# DustEmisBaseData — Dust emission state and flux data
# Ported from: dust_emis_base_type in DustEmisBase.F90
# ---------------------------------------------------------------------------

"""
    DustEmisBaseData

Dust emission state and flux data structure.
Contains overlap factors, size distribution parameters, saltation factor,
emission fluxes, and turbulent deposition velocities.

Ported from `dust_emis_base_type` in `DustEmisBase.F90`.
"""
Base.@kwdef mutable struct DustEmisBaseData
    # Overlap factors between source and sink distributions [dst_src_nbr × ndst]
    ovr_src_snk_mss::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    # [m] Mass-weighted mean diameter resolved [ndst]
    dmt_vwr::Vector{Float64} = Float64[]
    # [frc] Correction to Stokes settling velocity [ndst]
    stk_crc::Vector{Float64} = Float64[]
    # Factor in saltation computation
    saltation_factor::Float64 = 0.0
    # Surface dust emission (kg/m**2/s) [np × ndst] [+ = to atm]
    flx_mss_vrt_dst_patch::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    # Total dust flux into atmosphere [np]
    flx_mss_vrt_dst_tot_patch::Vector{Float64} = Float64[]
    # Turbulent deposition velocity (m/s) [np × ndst]
    vlc_trb_patch::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    # Turbulent deposition velocity 1 (m/s) [np]
    vlc_trb_1_patch::Vector{Float64} = Float64[]
    # Turbulent deposition velocity 2 (m/s) [np]
    vlc_trb_2_patch::Vector{Float64} = Float64[]
    # Turbulent deposition velocity 3 (m/s) [np]
    vlc_trb_3_patch::Vector{Float64} = Float64[]
    # Turbulent deposition velocity 4 (m/s) [np]
    vlc_trb_4_patch::Vector{Float64} = Float64[]
end

# ---------------------------------------------------------------------------
# dust_emis_init! — Allocate and initialize
# Ported from: InitBase / InitAllocateBase in DustEmisBase.F90
# ---------------------------------------------------------------------------

"""
    dust_emis_init!(dust, np)

Allocate and initialize all dust emission data arrays for `np` patches.
Also computes size distributions, overlap factors, saltation factor, and
Stokes correction via `init_dust_vars!`.

Ported from `InitBase` / `InitAllocateBase` in `DustEmisBase.F90`.
"""
function dust_emis_init!(dust::DustEmisBaseData, np::Int)
    dust.flx_mss_vrt_dst_patch     = fill(NaN, np, NDST)
    dust.flx_mss_vrt_dst_tot_patch = fill(NaN, np)
    dust.vlc_trb_patch             = fill(NaN, np, NDST)
    dust.vlc_trb_1_patch           = fill(NaN, np)
    dust.vlc_trb_2_patch           = fill(NaN, np)
    dust.vlc_trb_3_patch           = fill(NaN, np)
    dust.vlc_trb_4_patch           = fill(NaN, np)

    dust.ovr_src_snk_mss = fill(NaN, DST_SRC_NBR, NDST)
    dust.dmt_vwr         = fill(NaN, NDST)
    dust.stk_crc         = fill(NaN, NDST)

    # Compute size distributions, overlap factors, saltation factor, stk_crc
    init_dust_vars!(dust)

    return nothing
end

# ---------------------------------------------------------------------------
# dust_emis_clean! — Deallocate
# Ported from: CleanBase in DustEmisBase.F90
# ---------------------------------------------------------------------------

"""
    dust_emis_clean!(dust)

Deallocate (reset to empty) all dust emission data arrays.

Ported from `CleanBase` in `DustEmisBase.F90`.
"""
function dust_emis_clean!(dust::DustEmisBaseData)
    dust.flx_mss_vrt_dst_patch     = Matrix{Float64}(undef, 0, 0)
    dust.flx_mss_vrt_dst_tot_patch = Float64[]
    dust.vlc_trb_patch             = Matrix{Float64}(undef, 0, 0)
    dust.vlc_trb_1_patch           = Float64[]
    dust.vlc_trb_2_patch           = Float64[]
    dust.vlc_trb_3_patch           = Float64[]
    dust.vlc_trb_4_patch           = Float64[]

    dust.ovr_src_snk_mss = Matrix{Float64}(undef, 0, 0)
    dust.dmt_vwr         = Float64[]
    dust.stk_crc         = Float64[]

    return nothing
end

# ---------------------------------------------------------------------------
# init_dust_vars! — Compute size distributions, overlap factors,
#                   saltation factor, and Stokes correction
# Ported from: subroutine InitDustVars in DustEmisBase.F90
# ---------------------------------------------------------------------------

"""
    init_dust_vars!(dust)

Compute source efficiency factor from topography, overlap factors
(`ovr_src_snk_mss`), `saltation_factor`, particle diameters (`dmt_vwr`),
and Stokes correction factors (`stk_crc`).

Source: Paul Ginoux (source efficiency factor), C. Zender (dust model),
S. Levis (modifications).

Ported from `InitDustVars` in `DustEmisBase.F90`.
"""
function init_dust_vars!(dust::DustEmisBaseData)
    ndst = NDST
    dst_src_nbr = DST_SRC_NBR
    sz_nbr = SZ_NBR

    # Source distribution parameters (BSM96 p.73 Table 2)
    dmt_vma_src = [0.832e-6, 4.82e-6, 19.38e-6]
    gsd_anl_src = [2.10, 1.90, 1.60]
    mss_frc_src = [0.036, 0.957, 0.007]

    # Particle diameter grid [m]
    dmt_grd = [0.1e-6, 1.0e-6, 2.5e-6, 5.0e-6, 10.0e-6]

    # Saltation parameters
    dmt_slt_opt = 75.0e-6   # [m] Optimal diameter for saltation
    dns_slt = 2650.0         # [kg m-3] Density of optimal saltation particles

    # ---- Compute overlap factors between source and sink distributions ----
    # From szdstlgn.F subroutine ovr_src_snk_frc_get
    for m in 1:dst_src_nbr
        sqrt2lngsdi = sqrt(2.0) * log(gsd_anl_src[m])
        for n in 1:ndst
            lndmaxjovrdmdni = log(dmt_grd[n+1] / dmt_vma_src[m])
            lndminjovrdmdni = log(dmt_grd[n]   / dmt_vma_src[m])
            ovr_src_snk_frc = 0.5 * (_dust_erf(lndmaxjovrdmdni / sqrt2lngsdi) -
                                      _dust_erf(lndminjovrdmdni / sqrt2lngsdi))
            dust.ovr_src_snk_mss[m, n] = ovr_src_snk_frc * mss_frc_src[m]
        end
    end

    # ---- Compute saltation factor ----
    # From subroutine wnd_frc_thr_slt_get
    ryn_nbr_frc_thr_prx_opt = 0.38 + 1331.0 * (100.0 * dmt_slt_opt)^1.56

    if ryn_nbr_frc_thr_prx_opt < 0.03
        error("dstmbl: ryn_nbr_frc_thr_prx_opt < 0.03")
    elseif ryn_nbr_frc_thr_prx_opt < 10.0
        ryn_nbr_frc_thr_opt_fnc = -1.0 + 1.928 * (ryn_nbr_frc_thr_prx_opt^0.0922)
        ryn_nbr_frc_thr_opt_fnc = 0.1291 * 0.1291 / ryn_nbr_frc_thr_opt_fnc
    else
        ryn_nbr_frc_thr_opt_fnc = 1.0 - 0.0858 * exp(-0.0617 * (ryn_nbr_frc_thr_prx_opt - 10.0))
        ryn_nbr_frc_thr_opt_fnc = 0.120 * 0.120 * ryn_nbr_frc_thr_opt_fnc * ryn_nbr_frc_thr_opt_fnc
    end

    icf_fct = 1.0 + 6.0e-07 / (dns_slt * GRAV * (dmt_slt_opt^2.5))
    dns_fct = dns_slt * GRAV * dmt_slt_opt
    dust.saltation_factor = sqrt(icf_fct * dns_fct * ryn_nbr_frc_thr_opt_fnc)

    # ---- Particle diameter and settling velocity ----
    # From Charlie Zender's subroutines dst_psd_ini, dst_sz_rsl, grd_mk
    if ndst != 4
        error("Dustini error: ndst must equal 4 with current code")
    end

    dmt_min = zeros(ndst)
    dmt_max = zeros(ndst)
    dmt_ctr = zeros(ndst)
    dmt_dlt = zeros(ndst)

    for n in 1:ndst
        dmt_min[n] = dmt_grd[n]         # [m] Min diameter in bin
        dmt_max[n] = dmt_grd[n+1]       # [m] Max diameter in bin
        dmt_ctr[n] = 0.5 * (dmt_min[n] + dmt_max[n])  # [m] Diameter at bin center
        dmt_dlt[n] = dmt_max[n] - dmt_min[n]           # [m] Width of size bin
    end

    # Bin physical properties
    gsd_anl = 2.0       # [frc] Geometric std dev PaG77 p. 2080 Table1
    ln_gsd = log(gsd_anl)

    # Mass median diameter analytic She84 p.75 Table1
    dmt_vma = 3.5000e-6  # [m]

    # Convert mass median diameter to number median diameter
    dmt_nma = dmt_vma * exp(-3.0 * ln_gsd * ln_gsd)  # [m]

    # Compute resolved size statistics for each size distribution
    vlm_rsl = zeros(ndst)
    sz_min = zeros(sz_nbr)
    sz_max = zeros(sz_nbr)
    sz_ctr = zeros(sz_nbr)
    sz_dlt = zeros(sz_nbr)

    for n in 1:ndst
        series_ratio = (dmt_max[n] / dmt_min[n])^(1.0 / sz_nbr)
        sz_min[1] = dmt_min[n]
        for m_idx in 2:sz_nbr
            sz_min[m_idx] = sz_min[m_idx-1] * series_ratio
        end

        # Derived grid values
        for m_idx in 1:(sz_nbr-1)
            sz_max[m_idx] = sz_min[m_idx+1]
        end
        sz_max[sz_nbr] = dmt_max[n]

        # Final derived grid values
        for m_idx in 1:sz_nbr
            sz_ctr[m_idx] = 0.5 * (sz_min[m_idx] + sz_max[m_idx])
            sz_dlt[m_idx] = sz_max[m_idx] - sz_min[m_idx]
        end

        lngsdsqrttwopi_rcp = 1.0 / (ln_gsd * sqrt(2.0 * π))
        dust.dmt_vwr[n] = 0.0  # [m] Mass weighted diameter resolved
        vlm_rsl[n] = 0.0       # [m3 m-3] Volume concentration resolved

        for m_idx in 1:sz_nbr
            # Evaluate lognormal distribution for these sizes
            tmp = log(sz_ctr[m_idx] / dmt_nma) / ln_gsd
            lgn_dst = lngsdsqrttwopi_rcp * exp(-0.5 * tmp * tmp) / sz_ctr[m_idx]

            # Integrate moments of size distribution
            dust.dmt_vwr[n] += sz_ctr[m_idx] *
                π / 6.0 * (sz_ctr[m_idx]^3.0) *  # [m3] Volume
                lgn_dst * sz_dlt[m_idx]            # [# m-3] Number concentration

            vlm_rsl[n] +=
                π / 6.0 * (sz_ctr[m_idx]^3.0) *  # [m3] Volume
                lgn_dst * sz_dlt[m_idx]            # [# m-3] Number concentration
        end

        dust.dmt_vwr[n] = dust.dmt_vwr[n] / vlm_rsl[n]  # [m] Mass weighted diameter
    end

    # ---- Calculate correction to Stokes' settling velocity ----
    # From subroutine stk_crc_get
    eps_max = 1.0e-4
    dns_mdp = 100000.0 / (295.0 * RAIR)  # [kg m-3] const prs_mdp & tpt_vrt

    # Size-independent thermokinetic properties
    vsc_dyn_atm = 1.72e-5 * ((295.0 / 273.0)^1.5) * 393.0 /
        (295.0 + 120.0)  # [kg m-1 s-1] RoY94 p.102 tpt_mdp=295.0
    mfp_atm = 2.0 * vsc_dyn_atm /
        (100000.0 * sqrt(8.0 / (π * RAIR * 295.0)))  # [m] SeP97 p.455
    vsc_knm_atm = vsc_dyn_atm / dns_mdp  # [m2 s-1] Kinematic viscosity of air

    slp_crc = zeros(ndst)
    vlc_stk = zeros(ndst)
    vlc_grv = zeros(ndst)
    ryn_nbr_grv = zeros(ndst)
    cff_drg_grv = zeros(ndst)

    for m in 1:ndst
        slp_crc[m] = 1.0 + 2.0 * mfp_atm *
            (1.257 + 0.4 * exp(-1.1 * dust.dmt_vwr[m] / (2.0 * mfp_atm))) /
            dust.dmt_vwr[m]  # [frc] Slip correction factor SeP97 p.464
        vlc_stk[m] = (1.0 / 18.0) * dust.dmt_vwr[m] * dust.dmt_vwr[m] * DNS_AER *
            GRAV * slp_crc[m] / vsc_dyn_atm  # [m s-1] SeP97 p.466
    end

    # Iterative solution for drag coefficient, Reynolds number, terminal velocity
    for m in 1:ndst
        eps_crr = eps_max + 1.0  # [frc] Current relative accuracy
        itr_idx = 0              # [idx] Counting index

        # Initial guess for vlc_grv is exact for Re < 0.1
        vlc_grv[m] = vlc_stk[m]  # [m s-1]

        while eps_crr > eps_max
            # Save terminal velocity for convergence test
            vlc_grv_old = vlc_grv[m]  # [m s-1]
            ryn_nbr_grv[m] = vlc_grv[m] * dust.dmt_vwr[m] / vsc_knm_atm  # SeP97 p.460

            # Update drag coefficient based on new Reynolds number
            if ryn_nbr_grv[m] < 0.1
                cff_drg_grv[m] = 24.0 / ryn_nbr_grv[m]  # Stokes' law Sep97 p.463 (8.32)
            elseif ryn_nbr_grv[m] < 2.0
                cff_drg_grv[m] = (24.0 / ryn_nbr_grv[m]) *
                    (1.0 + 3.0 * ryn_nbr_grv[m] / 16.0 +
                     9.0 * ryn_nbr_grv[m] * ryn_nbr_grv[m] *
                     log(2.0 * ryn_nbr_grv[m]) / 160.0)  # Sep97 p.463 (8.32)
            elseif ryn_nbr_grv[m] < 500.0
                cff_drg_grv[m] = (24.0 / ryn_nbr_grv[m]) *
                    (1.0 + 0.15 * ryn_nbr_grv[m]^0.687)  # Sep97 p.463 (8.32)
            elseif ryn_nbr_grv[m] < 2.0e5
                cff_drg_grv[m] = 0.44  # Sep97 p.463 (8.32)
            else
                error("Dustini error: Reynolds number too large in stk_crc_get(): $(ryn_nbr_grv[m])")
            end

            # Update terminal velocity based on new Reynolds number and drag coeff
            # [m s-1] Terminal velocity SeP97 p.467 (8.44)
            vlc_grv[m] = sqrt(4.0 * GRAV * dust.dmt_vwr[m] * slp_crc[m] * DNS_AER /
                (3.0 * cff_drg_grv[m] * dns_mdp))
            eps_crr = abs((vlc_grv[m] - vlc_grv_old) / vlc_grv[m])  # Relative convergence

            if itr_idx == 12
                # Numerical pingpong may occur when Re = 0.1, 2.0, or 500.0
                vlc_grv[m] = 0.5 * (vlc_grv[m] + vlc_grv_old)  # [m s-1]
            end
            if itr_idx > 20
                @warn "Dustini: Terminal velocity not converging in stk_crc_get(), breaking loop..."
                break
            end
            itr_idx += 1
        end  # while
    end  # m loop

    # Compute factors to convert Stokes' settling velocities to actual
    for m in 1:ndst
        dust.stk_crc[m] = vlc_grv[m] / vlc_stk[m]
    end

    return nothing
end

# ---------------------------------------------------------------------------
# dust_dry_dep! — Turbulent dry deposition for dust
# Ported from: subroutine DustDryDep in DustEmisBase.F90
#
# Determine turbulent dry deposition for dust. Calculate the turbulent
# component of dust dry deposition (the turbulent deposition velocity
# through the lowest atmospheric layer). CAM will calculate the settling
# velocity through the whole atmospheric column.
# Source: C. Zender's dry deposition code
# ---------------------------------------------------------------------------

"""
    dust_dry_dep!(dust, patch_active, patch_column, bounds_p,
                  forc_pbot, forc_rho, forc_t, ram1, fv)

Calculate turbulent dry deposition velocity for dust.

Arguments:
- `dust`: DustEmisBaseData instance (modified in-place)
- `patch_active`: BitVector of active patches
- `patch_column`: Vector mapping patch → column
- `bounds_p`: patch index range
- `forc_pbot`: atmospheric pressure per column [Pa]
- `forc_rho`: atmospheric density per column [kg/m3]
- `forc_t`: atmospheric temperature per column [K]
- `ram1`: aerodynamical resistance per patch [s/m]
- `fv`: friction velocity per patch [m/s]

Ported from `DustDryDep` in `DustEmisBase.F90`.
"""
function dust_dry_dep!(
    dust::DustEmisBaseData,
    patch_active::BitVector,
    patch_column::Vector{Int},
    bounds_p::UnitRange{Int},
    # Atmospheric forcing (column-level)
    forc_pbot::Vector{Float64},      # [Pa] atm pressure
    forc_rho::Vector{Float64},       # [kg/m3] atm density
    forc_t::Vector{Float64},         # [K] atm temperature
    # Friction velocity (patch-level)
    ram1::Vector{Float64},           # [s/m] aerodynamical resistance
    fv::Vector{Float64}              # [m/s] friction velocity
)
    ndst = NDST
    shm_nbr_xpn_lnd = -2.0 / 3.0  # [frc] shm_nbr_xpn over land

    np = length(patch_active)
    vsc_dyn_atm = zeros(np)
    vsc_knm_atm = zeros(np)
    slp_crc = zeros(np, ndst)
    vlc_grv = zeros(np, ndst)
    rss_lmn = zeros(np, ndst)

    # First pass: size-independent thermokinetic properties and settling velocity
    for p in bounds_p
        patch_active[p] || continue
        c = patch_column[p]

        # Quasi-laminar layer resistance: size-independent thermokinetic properties
        vsc_dyn_atm[p] = 1.72e-5 * ((forc_t[c] / 273.0)^1.5) * 393.0 /
            (forc_t[c] + 120.0)  # [kg m-1 s-1] RoY94 p. 102
        mfp_atm = 2.0 * vsc_dyn_atm[p] /
            (forc_pbot[c] * sqrt(8.0 / (π * RAIR * forc_t[c])))  # [m] SeP97 p. 455
        vsc_knm_atm[p] = vsc_dyn_atm[p] / forc_rho[c]  # [m2 s-1] Kinematic viscosity

        for m in 1:ndst
            slp_crc[p, m] = 1.0 + 2.0 * mfp_atm *
                (1.257 + 0.4 * exp(-1.1 * dust.dmt_vwr[m] / (2.0 * mfp_atm))) /
                dust.dmt_vwr[m]  # [frc] Slip correction factor SeP97 p. 464
            vlc_grv[p, m] = (1.0 / 18.0) * dust.dmt_vwr[m] * dust.dmt_vwr[m] * DNS_AER *
                GRAV * slp_crc[p, m] / vsc_dyn_atm[p]  # [m s-1] Stokes' settling velocity
            vlc_grv[p, m] = vlc_grv[p, m] * dust.stk_crc[m]  # [m s-1] Corrected settling velocity
        end
    end

    # Second pass: quasi-laminar layer resistance
    for m in 1:ndst
        for p in bounds_p
            patch_active[p] || continue
            c = patch_column[p]

            stk_nbr = vlc_grv[p, m] * fv[p] * fv[p] / (GRAV * vsc_knm_atm[p])  # [frc] SeP97 p.965
            dff_aer = BOLTZ * forc_t[c] * slp_crc[p, m] /
                (3.0 * π * vsc_dyn_atm[p] * dust.dmt_vwr[m])  # [m2 s-1] SeP97 p.474
            shm_nbr = vsc_knm_atm[p] / dff_aer  # [frc] SeP97 p.972
            shm_nbr_xpn = shm_nbr_xpn_lnd  # [frc]

            tmp = shm_nbr^shm_nbr_xpn + 10.0^(-3.0 / stk_nbr)
            rss_lmn[p, m] = 1.0 / (tmp * fv[p])  # [s m-1] SeP97 p.972,965
        end
    end

    # Third pass: lowest layer turbulent deposition
    for m in 1:ndst
        for p in bounds_p
            patch_active[p] || continue
            rss_trb = ram1[p] + rss_lmn[p, m] + ram1[p] * rss_lmn[p, m] * vlc_grv[p, m]  # [s m-1]
            dust.vlc_trb_patch[p, m] = 1.0 / rss_trb  # [m s-1]
        end
    end

    # Copy to individual bin arrays
    for p in bounds_p
        patch_active[p] || continue
        dust.vlc_trb_1_patch[p] = dust.vlc_trb_patch[p, 1]
        dust.vlc_trb_2_patch[p] = dust.vlc_trb_patch[p, 2]
        dust.vlc_trb_3_patch[p] = dust.vlc_trb_patch[p, 3]
        dust.vlc_trb_4_patch[p] = dust.vlc_trb_patch[p, 4]
    end

    return nothing
end

# ---------------------------------------------------------------------------
# check_dust_emis_is_valid — Validate dust emission consistency
# Ported from: subroutine CheckDustEmisIsValid in DustEmisBase.F90
# ---------------------------------------------------------------------------

"""
    check_dust_emis_is_valid(dust, p)

Check that dust emission state for patch `p` is valid.
Ensures that total dust is the sum of all dust bins, and that dry
deposition for each bin agrees with the array of all bins.

Returns `true` if valid, throws an error otherwise.

Ported from `CheckDustEmisIsValid` in `DustEmisBase.F90`.
"""
function check_dust_emis_is_valid(dust::DustEmisBaseData, p::Int)
    if abs(sum(dust.flx_mss_vrt_dst_patch[p, :]) - dust.flx_mss_vrt_dst_tot_patch[p]) > 0.0
        error("Sum over dust bins does NOT equal total dust for patch $p")
    end
    if dust.vlc_trb_patch[p, 1] != dust.vlc_trb_1_patch[p]
        error("Dry deposition for dust bin 1 not equal to the array bin for patch $p")
    end
    if dust.vlc_trb_patch[p, 2] != dust.vlc_trb_2_patch[p]
        error("Dry deposition for dust bin 2 not equal to the array bin for patch $p")
    end
    if dust.vlc_trb_patch[p, 3] != dust.vlc_trb_3_patch[p]
        error("Dry deposition for dust bin 3 not equal to the array bin for patch $p")
    end
    if dust.vlc_trb_patch[p, 4] != dust.vlc_trb_4_patch[p]
        error("Dry deposition for dust bin 4 not equal to the array bin for patch $p")
    end

    return true
end

# ---------------------------------------------------------------------------
# get_dust_patch_vars — Get patch-level dust variables
# Ported from: subroutine GetPatchVars in DustEmisBase.F90
# ---------------------------------------------------------------------------

"""
    get_dust_patch_vars(dust, p)

Get important dust variables on the given patch. Returns a NamedTuple
with all flux and deposition velocity fields.

Ported from `GetPatchVars` in `DustEmisBase.F90`.
"""
function get_dust_patch_vars(dust::DustEmisBaseData, p::Int)
    return (
        flx_mss_vrt_dst     = dust.flx_mss_vrt_dst_patch[p, :],
        flx_mss_vrt_dst_tot = dust.flx_mss_vrt_dst_tot_patch[p],
        vlc_trb             = dust.vlc_trb_patch[p, :],
        vlc_trb_1           = dust.vlc_trb_1_patch[p],
        vlc_trb_2           = dust.vlc_trb_2_patch[p],
        vlc_trb_3           = dust.vlc_trb_3_patch[p],
        vlc_trb_4           = dust.vlc_trb_4_patch[p],
    )
end

# ---------------------------------------------------------------------------
# get_dust_const_vars — Get constant dust variables
# Ported from: subroutine GetConstVars in DustEmisBase.F90
# ---------------------------------------------------------------------------

"""
    get_dust_const_vars(dust)

Get the saltation factor constant.

Ported from `GetConstVars` in `DustEmisBase.F90`.
"""
function get_dust_const_vars(dust::DustEmisBaseData)
    return dust.saltation_factor
end

# ---------------------------------------------------------------------------
# write_dust_patch_to_log — Write patch info to log
# Ported from: subroutine WritePatchToLog in DustEmisBase.F90
# ---------------------------------------------------------------------------

"""
    write_dust_patch_to_log(dust, p)

Write out information on dust emissions for the given patch.

Ported from `WritePatchToLog` in `DustEmisBase.F90`.
"""
function write_dust_patch_to_log(dust::DustEmisBaseData, p::Int)
    @info "flx_mss_vrt_dst_tot" dust.flx_mss_vrt_dst_tot_patch[p]
    @info "vlc_trb_1" dust.vlc_trb_1_patch[p]
    @info "vlc_trb_2" dust.vlc_trb_2_patch[p]
    @info "vlc_trb_3" dust.vlc_trb_3_patch[p]
    @info "vlc_trb_4" dust.vlc_trb_4_patch[p]
    return nothing
end
