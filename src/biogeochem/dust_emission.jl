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
# _dust_erf: alias for shared erf() from varcon.jl
const _dust_erf = erf

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
Base.@kwdef mutable struct DustEmisBaseData{FT<:Real,
                                V<:AbstractVector{FT},
                                M<:AbstractMatrix{FT}}
    # Overlap factors between source and sink distributions [dst_src_nbr × ndst]
    ovr_src_snk_mss::M = Matrix{Float64}(undef, 0, 0)
    # [m] Mass-weighted mean diameter resolved [ndst]
    dmt_vwr::V = Float64[]
    # [frc] Correction to Stokes settling velocity [ndst]
    stk_crc::V = Float64[]
    # Factor in saltation computation
    saltation_factor::FT = 0.0
    # Surface dust emission (kg/m**2/s) [np × ndst] [+ = to atm]
    flx_mss_vrt_dst_patch::M = Matrix{Float64}(undef, 0, 0)
    # Total dust flux into atmosphere [np]
    flx_mss_vrt_dst_tot_patch::V = Float64[]
    # Turbulent deposition velocity (m/s) [np × ndst]
    vlc_trb_patch::M = Matrix{Float64}(undef, 0, 0)
    # Turbulent deposition velocity 1 (m/s) [np]
    vlc_trb_1_patch::V = Float64[]
    # Turbulent deposition velocity 2 (m/s) [np]
    vlc_trb_2_patch::V = Float64[]
    # Turbulent deposition velocity 3 (m/s) [np]
    vlc_trb_3_patch::V = Float64[]
    # Turbulent deposition velocity 4 (m/s) [np]
    vlc_trb_4_patch::V = Float64[]
end

DustEmisBaseData{FT}(; kwargs...) where {FT<:Real} =
    DustEmisBaseData{FT, Vector{FT}, Matrix{FT}}(; kwargs...)
Adapt.@adapt_structure DustEmisBaseData


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
function dust_emis_init!(dust::DustEmisBaseData{FT}, np::Int) where {FT}
    dust.flx_mss_vrt_dst_patch     = fill(FT(NaN), np, NDST)
    dust.flx_mss_vrt_dst_tot_patch = fill(FT(NaN), np)
    dust.vlc_trb_patch             = fill(FT(NaN), np, NDST)
    dust.vlc_trb_1_patch           = fill(FT(NaN), np)
    dust.vlc_trb_2_patch           = fill(FT(NaN), np)
    dust.vlc_trb_3_patch           = fill(FT(NaN), np)
    dust.vlc_trb_4_patch           = fill(FT(NaN), np)

    dust.ovr_src_snk_mss = fill(FT(NaN), DST_SRC_NBR, NDST)
    dust.dmt_vwr         = fill(FT(NaN), NDST)
    dust.stk_crc         = fill(FT(NaN), NDST)

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
function dust_emis_clean!(dust::DustEmisBaseData{FT}) where {FT}
    dust.flx_mss_vrt_dst_patch     = Matrix{FT}(undef, 0, 0)
    dust.flx_mss_vrt_dst_tot_patch = FT[]
    dust.vlc_trb_patch             = Matrix{FT}(undef, 0, 0)
    dust.vlc_trb_1_patch           = FT[]
    dust.vlc_trb_2_patch           = FT[]
    dust.vlc_trb_3_patch           = FT[]
    dust.vlc_trb_4_patch           = FT[]

    dust.ovr_src_snk_mss = Matrix{FT}(undef, 0, 0)
    dust.dmt_vwr         = FT[]
    dust.stk_crc         = FT[]

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
function init_dust_vars!(dust::DustEmisBaseData{FT}) where {FT}
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

# --------------------------------------------------------------------------
# Kernels for the per-patch passes of dust_dry_dep!. Float64 literals are
# converted to the working element type T (T(...)/one/zero) so no Float64
# reaches a Float32-only backend (Metal); on Float64 CPU this is byte-identical.
# Each thread owns one patch index p = lo + i - 1; the ndst (=4) bin loop runs
# in-thread. `BOLTZ`, `RAIR`, `GRAV`, `DNS_AER` are module constants (Float64);
# they are converted at the point of use.
# --------------------------------------------------------------------------
@kernel function _dust_pass1_kernel!(vsc_dyn_atm, vsc_knm_atm, slp_crc, vlc_grv,
                                     @Const(patch_active), @Const(patch_column),
                                     @Const(forc_pbot), @Const(forc_rho), @Const(forc_t),
                                     @Const(dmt_vwr), @Const(stk_crc), ndst::Int, lo::Int)
    i = @index(Global)
    @inbounds begin
        p = lo + i - 1
        if patch_active[p]
            T = eltype(vsc_dyn_atm)
            c = patch_column[p]
            tc = forc_t[c]
            vsc_dyn_atm[p] = T(1.72e-5) * ((tc / T(273.0))^T(1.5)) * T(393.0) /
                (tc + T(120.0))
            mfp_atm = T(2.0) * vsc_dyn_atm[p] /
                (forc_pbot[c] * sqrt(T(8.0) / (T(π) * T(RAIR) * tc)))
            vsc_knm_atm[p] = vsc_dyn_atm[p] / forc_rho[c]

            for m in 1:ndst
                dm = dmt_vwr[m]
                slp_crc[p, m] = one(T) + T(2.0) * mfp_atm *
                    (T(1.257) + T(0.4) * exp(-T(1.1) * dm / (T(2.0) * mfp_atm))) / dm
                vg = (one(T) / T(18.0)) * dm * dm * T(DNS_AER) *
                    T(GRAV) * slp_crc[p, m] / vsc_dyn_atm[p]
                vlc_grv[p, m] = vg * stk_crc[m]
            end
        end
    end
end

@kernel function _dust_pass2_kernel!(rss_lmn, @Const(patch_active), @Const(patch_column),
                                     @Const(forc_t), @Const(fv), @Const(vlc_grv),
                                     @Const(vsc_knm_atm), @Const(vsc_dyn_atm),
                                     @Const(slp_crc), @Const(dmt_vwr), ndst::Int, lo::Int)
    i = @index(Global)
    @inbounds begin
        p = lo + i - 1
        if patch_active[p]
            T = eltype(rss_lmn)
            c = patch_column[p]
            shm_nbr_xpn = -T(2.0) / T(3.0)  # shm_nbr_xpn over land
            for m in 1:ndst
                stk_nbr = vlc_grv[p, m] * fv[p] * fv[p] / (T(GRAV) * vsc_knm_atm[p])
                dff_aer = T(BOLTZ) * forc_t[c] * slp_crc[p, m] /
                    (T(3.0) * T(π) * vsc_dyn_atm[p] * dmt_vwr[m])
                shm_nbr = vsc_knm_atm[p] / dff_aer
                tmp = shm_nbr^shm_nbr_xpn + T(10.0)^(-T(3.0) / stk_nbr)
                rss_lmn[p, m] = one(T) / (tmp * fv[p])
            end
        end
    end
end

@kernel function _dust_pass3_kernel!(vlc_trb_patch, vlc_trb_1_patch, vlc_trb_2_patch,
                                     vlc_trb_3_patch, vlc_trb_4_patch,
                                     @Const(patch_active), @Const(ram1), @Const(rss_lmn),
                                     @Const(vlc_grv), ndst::Int, lo::Int)
    i = @index(Global)
    @inbounds begin
        p = lo + i - 1
        if patch_active[p]
            T = eltype(vlc_trb_patch)
            for m in 1:ndst
                rss_trb = ram1[p] + rss_lmn[p, m] + ram1[p] * rss_lmn[p, m] * vlc_grv[p, m]
                vlc_trb_patch[p, m] = one(T) / rss_trb
            end
            vlc_trb_1_patch[p] = vlc_trb_patch[p, 1]
            vlc_trb_2_patch[p] = vlc_trb_patch[p, 2]
            vlc_trb_3_patch[p] = vlc_trb_patch[p, 3]
            vlc_trb_4_patch[p] = vlc_trb_patch[p, 4]
        end
    end
end

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
    patch_active::AbstractVector{Bool},
    patch_column::AbstractVector{<:Integer},
    bounds_p::UnitRange{Int},
    # Atmospheric forcing (column-level)
    forc_pbot::AbstractVector{<:Real},      # [Pa] atm pressure
    forc_rho::AbstractVector{<:Real},       # [kg/m3] atm density
    forc_t::AbstractVector{<:Real},         # [K] atm temperature
    # Friction velocity (patch-level)
    ram1::AbstractVector{<:Real},           # [s/m] aerodynamical resistance
    fv::AbstractVector{<:Real}              # [m/s] friction velocity
)
    ndst = NDST

    np = length(patch_active)
    FT = eltype(dust.vlc_trb_patch)

    # Scratch is device-resident (similar() off a device struct field) and zeroed
    # via fill! so no host `zeros` array reaches a GPU kernel.
    vsc_dyn_atm = fill!(similar(dust.vlc_trb_1_patch, np), zero(FT))
    vsc_knm_atm = fill!(similar(dust.vlc_trb_1_patch, np), zero(FT))
    slp_crc     = fill!(similar(dust.vlc_trb_patch, np, ndst), zero(FT))
    vlc_grv     = fill!(similar(dust.vlc_trb_patch, np, ndst), zero(FT))
    rss_lmn     = fill!(similar(dust.vlc_trb_patch, np, ndst), zero(FT))

    lo = first(bounds_p)
    n  = length(bounds_p)

    # Pass 1: thermokinetic properties + Stokes settling velocity (per patch).
    _launch!(_dust_pass1_kernel!, vsc_dyn_atm, vsc_knm_atm, slp_crc, vlc_grv,
             patch_active, patch_column, forc_pbot, forc_rho, forc_t,
             dust.dmt_vwr, dust.stk_crc, ndst, lo; ndrange = n)

    # Pass 2: quasi-laminar layer resistance (per patch).
    _launch!(_dust_pass2_kernel!, rss_lmn, patch_active, patch_column, forc_t,
             fv, vlc_grv, vsc_knm_atm, vsc_dyn_atm, slp_crc, dust.dmt_vwr, ndst, lo;
             ndrange = n)

    # Pass 3: lowest-layer turbulent deposition + copy to bin arrays (per patch).
    _launch!(_dust_pass3_kernel!, dust.vlc_trb_patch, dust.vlc_trb_1_patch,
             dust.vlc_trb_2_patch, dust.vlc_trb_3_patch, dust.vlc_trb_4_patch,
             patch_active, ram1, rss_lmn, vlc_grv, ndst, lo; ndrange = n)

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
function check_dust_emis_is_valid(dust::DustEmisBaseData{FT}, p::Int) where {FT}
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
function get_dust_patch_vars(dust::DustEmisBaseData{FT}, p::Int) where {FT}
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
function write_dust_patch_to_log(dust::DustEmisBaseData{FT}, p::Int) where {FT}
    @info "flx_mss_vrt_dst_tot" dust.flx_mss_vrt_dst_tot_patch[p]
    @info "vlc_trb_1" dust.vlc_trb_1_patch[p]
    @info "vlc_trb_2" dust.vlc_trb_2_patch[p]
    @info "vlc_trb_3" dust.vlc_trb_3_patch[p]
    @info "vlc_trb_4" dust.vlc_trb_4_patch[p]
    return nothing
end
