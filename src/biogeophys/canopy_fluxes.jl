# ==========================================================================
# Ported from: src/biogeophys/CanopyFluxesMod.F90 (1756 lines)
# Calculates leaf temperature and surface fluxes using Newton-Raphson
# iteration to solve the surface energy budget:
#   f(t_veg) = Net radiation - Sensible - Latent = 0
# ==========================================================================

# --- Local constants (from CanopyFluxes subroutine) ---
const BTRAN0             = 0.0     # initial btran value
const ZII_CANOPY         = 1000.0  # convective boundary layer height [m]
const BETA_CANOPY        = 1.0     # coefficient of convective velocity [-]
const DELMAX_CANOPY      = 1.0     # max change in leaf temperature [K]
const DLEMIN_CANOPY      = 0.1     # max limit for energy flux convergence [W/m2]
const DTMIN_CANOPY       = 0.01    # max limit for temperature convergence [K]
const ITMIN_CANOPY       = 2       # minimum number of iterations [-]
const RIA_CANOPY         = 0.5     # free parameter for stable formulation
const ABOVE_CANOPY       = 1       # index for above-canopy resistances
const BELOW_CANOPY       = 2       # index for below-canopy resistances

# Biomass heat storage tuning parameters
const K_VERT_CANOPY      = 0.1     # vertical distribution of stem
const K_CYL_VOL_CANOPY   = 1.0     # departure from cylindrical volume
const K_CYL_AREA_CANOPY  = 1.0     # departure from cylindrical area
const K_INTERNAL_CANOPY  = 0.0     # self-absorption of leaf/stem longwave
const MIN_STEM_DIAMETER  = 0.05    # minimum stem diameter [m]

# --- Parameters type (from params_type in Fortran) ---

"""
    CanopyFluxesParamsData

Parameters for canopy fluxes, read from parameter file in Fortran.
Ported from `params_type` in `CanopyFluxesMod.F90`.
"""
Base.@kwdef mutable struct CanopyFluxesParamsData
    lai_dl   ::Float64 = 0.5    # plant litter area index [m2/m2]
    z_dl     ::Float64 = 0.05   # litter layer thickness [m]
    a_coef   ::Float64 = 0.5    # drag coefficient under less dense canopy [-]
    a_exp    ::Float64 = 1.0    # drag exponent under less dense canopy [-]
    csoilc   ::Float64 = 0.004  # soil drag coefficient under dense canopy [-]
    cv       ::Float64 = 0.01   # turbulent transfer coeff canopy surface-to-air [m/s^(1/2)]
    wind_min ::Float64 = 1.0    # minimum wind speed at forcing height [m/s]
end

const canopy_fluxes_params = CanopyFluxesParamsData()

# --- Module control state ---

"""
    CanopyFluxesControl

Module-level control flags for canopy fluxes.
Ported from module variables in `CanopyFluxesMod.F90`.
"""
Base.@kwdef mutable struct CanopyFluxesControl
    perchroot                ::Bool = false   # btran based only on unfrozen soil levels
    perchroot_alt            ::Bool = false   # btran based on active layer
    use_undercanopy_stability::Bool = false   # use undercanopy stability term
    itmax_canopy_fluxes      ::Int  = 40      # max iterations in CanopyFluxes
    use_biomass_heat_storage ::Bool = false   # use biomass heat storage
end

const canopy_fluxes_ctrl = CanopyFluxesControl()

# --- Init / clean / read functions ---

"""
    canopy_fluxes_read_nml!(; use_undercanopy_stability, use_biomass_heat_storage, itmax_canopy_fluxes)

Set namelist parameters for canopy fluxes.
Ported from `CanopyFluxesReadNML` in `CanopyFluxesMod.F90`.
"""
function canopy_fluxes_read_nml!(; use_undercanopy_stability::Bool = false,
                                   use_biomass_heat_storage::Bool = false,
                                   itmax_canopy_fluxes::Int = 40)
    if itmax_canopy_fluxes < 1
        error("canopy_fluxes_read_nml!: itmax_canopy_fluxes must be > 0")
    end
    canopy_fluxes_ctrl.use_undercanopy_stability = use_undercanopy_stability
    canopy_fluxes_ctrl.use_biomass_heat_storage = use_biomass_heat_storage
    canopy_fluxes_ctrl.itmax_canopy_fluxes = itmax_canopy_fluxes
    return nothing
end

"""
    canopy_fluxes_read_params!(; lai_dl, z_dl, a_coef, a_exp, csoilc, cv, wind_min)

Set parameter values for canopy fluxes.
Ported from `readParams` in `CanopyFluxesMod.F90`.
"""
function canopy_fluxes_read_params!(; lai_dl::Real = 0.5,
                                      z_dl::Real = 0.05,
                                      a_coef::Real = 0.5,
                                      a_exp::Real = 1.0,
                                      csoilc::Real = 0.004,
                                      cv::Real = 0.01,
                                      wind_min::Real = 1.0)
    canopy_fluxes_params.lai_dl = lai_dl
    canopy_fluxes_params.z_dl = z_dl
    canopy_fluxes_params.a_coef = a_coef
    canopy_fluxes_params.a_exp = a_exp
    canopy_fluxes_params.csoilc = csoilc
    canopy_fluxes_params.cv = cv
    canopy_fluxes_params.wind_min = wind_min
    return nothing
end

# =====================================================================
# Main canopy_fluxes! function
# =====================================================================

"""
    canopy_fluxes!(canopystate, energyflux, frictionvel, temperature,
                   solarabs, soilstate, waterfluxbulk, waterstatebulk,
                   waterdiagbulk, photosyns, patch_data, col_data, gridcell_data,
                   mask_exposedvegp, bounds_patch, bounds_col,
                   forc_lwrad_col, forc_q_col, forc_pbot_col, forc_th_col,
                   forc_rho_col, forc_t_col, forc_u_grc, forc_v_grc,
                   forc_pco2_grc, forc_po2_grc, forc_hgt_t_grc, forc_hgt_u_grc,
                   forc_hgt_q_grc, dayl_grc, max_dayl_grc,
                   downreg_patch, leafn_patch, dtime;
                   kwargs...)

Main canopy fluxes calculation. Solves the leaf energy balance iteratively
using Newton-Raphson to find leaf temperature (t_veg).

Ported from `subroutine CanopyFluxes` in `CanopyFluxesMod.F90`.
"""
function canopy_fluxes!(
        canopystate     ::CanopyStateData,
        energyflux      ::EnergyFluxData,
        frictionvel     ::FrictionVelocityData,
        temperature     ::TemperatureData,
        solarabs        ::SolarAbsorbedData,
        soilstate       ::SoilStateData,
        waterfluxbulk   ::WaterFluxBulkData,
        waterstatebulk  ::WaterStateBulkData,
        waterdiagbulk   ::WaterDiagnosticBulkData,
        photosyns       ::PhotosynthesisData,
        patch_data      ::PatchData,
        col_data        ::ColumnData,
        gridcell_data   ::GridcellData,
        mask_exposedvegp::BitVector,
        bounds_patch    ::UnitRange{Int},
        bounds_col      ::UnitRange{Int},
        # Atmospheric forcing (column-level)
        forc_lwrad_col  ::Vector{<:Real},
        forc_q_col      ::Vector{<:Real},
        forc_pbot_col   ::Vector{<:Real},
        forc_th_col     ::Vector{<:Real},
        forc_rho_col    ::Vector{<:Real},
        forc_t_col      ::Vector{<:Real},
        # Atmospheric forcing (gridcell-level)
        forc_u_grc      ::Vector{<:Real},
        forc_v_grc      ::Vector{<:Real},
        forc_pco2_grc   ::Vector{<:Real},
        forc_po2_grc    ::Vector{<:Real},
        forc_hgt_t_grc  ::Vector{<:Real},
        forc_hgt_u_grc  ::Vector{<:Real},
        forc_hgt_q_grc  ::Vector{<:Real},
        # Daylength (gridcell-level)
        dayl_grc        ::Vector{<:Real},
        max_dayl_grc    ::Vector{<:Real},
        # CN inputs (patch-level)
        downreg_patch   ::Vector{<:Real},
        leafn_patch     ::Vector{<:Real},
        # Time step
        dtime           ::Real;
        # PFT parameters (patch-indexed via ivt)
        dleaf_pft       ::Vector{<:Real} = fill(0.04, MXPFT+1),
        dbh_pft         ::Vector{<:Real} = fill(0.1, MXPFT+1),
        slatop_pft      ::Vector{<:Real} = fill(0.01, MXPFT+1),
        fbw_pft         ::Vector{<:Real} = fill(0.1, MXPFT+1),
        nstem_pft       ::Vector{<:Real} = fill(1.0, MXPFT+1),
        woody_pft       ::Vector{<:Real} = fill(0.0, MXPFT+1),
        rstem_per_dbh_pft::Vector{<:Real} = fill(0.0, MXPFT+1),
        wood_density_pft::Vector{<:Real} = fill(500.0, MXPFT+1),
        is_tree_pft     ::Vector{Bool} = fill(false, MXPFT+1),
        is_shrub_pft    ::Vector{Bool} = fill(false, MXPFT+1),
        c3psn_pft       ::Vector{<:Real} = fill(1.0, MXPFT+1),
        leafcn_pft      ::Vector{<:Real} = fill(25.0, MXPFT+1),
        flnr_pft        ::Vector{<:Real} = fill(0.08, MXPFT+1),
        fnitr_pft       ::Vector{<:Real} = fill(0.1, MXPFT+1),
        mbbopt_pft      ::Vector{<:Real} = fill(0.0, MXPFT+1),
        medlynintercept_pft::Vector{<:Real} = fill(100.0, MXPFT+1),
        medlynslope_pft ::Vector{<:Real} = fill(6.0, MXPFT+1),
        crop_pft        ::Vector{<:Real} = Float64[],
        z0v_Cr_pft      ::Vector{<:Real} = fill(0.35, MXPFT+1),
        z0v_Cs_pft      ::Vector{<:Real} = fill(0.003, MXPFT+1),
        z0v_c_pft       ::Vector{<:Real} = fill(0.25, MXPFT+1),
        z0v_cw_pft      ::Vector{<:Real} = fill(2.0, MXPFT+1),
        z0v_LAImax_pft  ::Vector{<:Real} = fill(8.0, MXPFT+1),
        # Feature flags
        use_cn          ::Bool = false,
        use_lch4        ::Bool = false,
        use_c13         ::Bool = false,
        use_hydrstress  ::Bool = false,
        use_fates       ::Bool = false,
        use_luna        ::Bool = false,
        z0param_method  ::String = "ZengWang2007",
        # Ozone stress arrays (optional)
        o3coefv_patch   ::Vector{<:Real} = Float64[],
        o3coefg_patch   ::Vector{<:Real} = Float64[],
        # CH4 conductance output (optional)
        grnd_ch4_cond_patch::Vector{<:Real} = Float64[],
        # Photosynthesis extras
        forc_pc13o2_grc ::Vector{<:Real} = Float64[],
        t10_patch       ::Vector{<:Real} = Float64[],
        nrad_patch      ::Vector{Int} = Int[],
        tlai_z_patch    ::Matrix{<:Real} = Matrix{Float64}(undef,0,0),
        vcmaxcint_sun_patch ::Vector{<:Real} = Float64[],
        vcmaxcint_sha_patch ::Vector{<:Real} = Float64[],
        parsun_z_patch  ::Matrix{<:Real} = Matrix{Float64}(undef,0,0),
        parsha_z_patch  ::Matrix{<:Real} = Matrix{Float64}(undef,0,0),
        laisun_z_patch  ::Matrix{<:Real} = Matrix{Float64}(undef,0,0),
        laisha_z_patch  ::Matrix{<:Real} = Matrix{Float64}(undef,0,0),
        leaf_mr_vcm     ::Real = 0.015,
        # Calibration overrides (NaN = use defaults)
        overrides       ::CalibrationOverrides = CalibrationOverrides())

    np = length(bounds_patch)
    begp = first(bounds_patch)
    endp = last(bounds_patch)
    nlevgrnd = varpar.nlevgrnd
    nlevsno  = varpar.nlevsno

    # --- Convenience aliases ---
    params = canopy_fluxes_params
    ctrl   = canopy_fluxes_ctrl

    # Local work arrays (patch-indexed)
    FT = eltype(forc_t_col)
    zldis     = zeros(FT, endp)
    dth       = zeros(FT, endp)
    dthv      = zeros(FT, endp)
    dqh       = zeros(FT, endp)
    ur        = zeros(FT, endp)
    temp1     = zeros(FT, endp)
    temp12m   = zeros(FT, endp)
    temp2     = zeros(FT, endp)
    temp22m   = zeros(FT, endp)
    rb        = zeros(FT, endp)
    rah       = zeros(FT, endp, 2)
    raw       = zeros(FT, endp, 2)
    wtg       = zeros(FT, endp)
    wta0      = zeros(FT, endp)
    wtl0      = zeros(FT, endp)
    wtstem0   = zeros(FT, endp)
    wtal      = zeros(FT, endp)
    wtga      = zeros(FT, endp)
    wtgq      = zeros(FT, endp)
    wtaq0     = zeros(FT, endp)
    wtlq0     = zeros(FT, endp)
    wtalq     = zeros(FT, endp)
    el        = zeros(FT, endp)
    qsatl     = zeros(FT, endp)
    qsatldT   = zeros(FT, endp)
    air       = zeros(FT, endp)
    bir       = zeros(FT, endp)
    cir       = zeros(FT, endp)
    del_arr   = zeros(FT, endp)  # "del" in Fortran (reserved word in Julia context for clarity)
    del2      = zeros(FT, endp)
    dele      = zeros(FT, endp)
    delq      = zeros(FT, endp)
    det_arr   = zeros(FT, endp)
    efeb      = zeros(FT, endp)
    efe       = zeros(FT, endp)
    err_arr   = zeros(FT, endp)
    obuold    = zeros(FT, endp)
    tlbef     = zeros(FT, endp)
    tl_ini    = zeros(FT, endp)
    ts_ini    = zeros(FT, endp)
    co2_arr   = zeros(FT, endp)
    o2_arr    = zeros(FT, endp)
    svpts     = zeros(FT, endp)
    eah       = zeros(FT, endp)
    dt_veg    = zeros(FT, endp)
    fm        = zeros(FT, endp)
    nmozsgn   = zeros(Int, endp)
    dayl_factor = zeros(FT, endp)
    rootsum   = zeros(FT, endp)
    dbh       = zeros(FT, endp)
    cp_leaf   = zeros(FT, endp)
    cp_stem   = zeros(FT, endp)
    rstem     = zeros(FT, endp)
    dt_stem   = zeros(FT, endp)
    frac_rad_abs_by_stem = zeros(FT, endp)
    lw_stem   = zeros(FT, endp)
    lw_leaf   = zeros(FT, endp)
    sa_stem   = zeros(FT, endp)
    sa_leaf   = zeros(FT, endp)
    sa_internal = zeros(FT, endp)
    uuc       = zeros(FT, endp)
    snocan_baseline = zeros(FT, endp)

    # Integer filter arrays (Fortran-style, for convergence testing)
    filterp = zeros(Int, np)
    fporig  = zeros(Int, np)

    # Build initial filter from mask
    fn = 0
    for p in bounds_patch
        if mask_exposedvegp[p]
            fn += 1
            filterp[fn] = p
        end
    end

    # Time step initialization of photosynthesis variables
    photosynthesis_timestep_init!(photosyns, mask_exposedvegp, bounds_patch)

    # =========================================================================
    # Phase 1: Initialization
    # =========================================================================

    for fi in 1:fn
        p = filterp[fi]
        c = patch_data.column[p]
        g = patch_data.gridcell[p]

        del_arr[p] = 0.0
        efeb[p]    = 0.0
        wtlq0[p]   = 0.0
        wtalq[p]   = 0.0
        wtgq[p]    = 0.0
        wtaq0[p]   = 0.0
        obuold[p]  = 0.0
        # btran_patch is computed by calc_root_moist_stress! — do NOT reset here
        energyflux.dhsdt_canopy_patch[p] = 0.0
        energyflux.eflx_sh_stem_patch[p] = 0.0
    end

    # --- Calculate biomass heat capacities ---
    if ctrl.use_biomass_heat_storage
        for fi in 1:fn
            p = filterp[fi]
            c = patch_data.column[p]
            ft = patch_data.itype[p] + 1  # 0-based Fortran PFT → 1-based Julia

            elai_p = canopystate.elai_patch[p]
            esai_p = canopystate.esai_patch[p]

            # fraction of stem receiving incoming radiation
            if (elai_p + esai_p) > 0.0
                frac_rad_abs_by_stem[p] = esai_p / (elai_p + esai_p)
            else
                frac_rad_abs_by_stem[p] = 0.0
            end
            if elai_p > 0.0
                frac_rad_abs_by_stem[p] = K_VERT_CANOPY * frac_rad_abs_by_stem[p]
            end

            # calculate dbh from stem biomass or parameter file
            if use_cn
                if canopystate.stem_biomass_patch[p] > 0.0
                    dbh[p] = 2.0 * sqrt(canopystate.stem_biomass_patch[p] * (1.0 - fbw_pft[ft]) /
                        (RPI * canopystate.htop_patch[p] * K_CYL_VOL_CANOPY *
                         nstem_pft[ft] * wood_density_pft[ft]))
                else
                    dbh[p] = 0.0
                end
            else
                dbh[p] = dbh_pft[ft]
            end

            # leaf and stem surface area
            sa_leaf[p] = 2.0 * elai_p  # double for full surface area (sensible heat)

            # stem surface area
            sa_stem[p] = nstem_pft[ft] * canopystate.htop_patch[p] * RPI * dbh[p]
            sa_stem[p] = K_CYL_AREA_CANOPY * sa_stem[p]

            # only separate leaf/stem heat capacity for trees/shrubs with sufficient dbh
            if !(is_tree_pft[ft] || is_shrub_pft[ft]) || dbh[p] < MIN_STEM_DIAMETER
                frac_rad_abs_by_stem[p] = 0.0
                sa_stem[p] = 0.0
                sa_leaf[p] = sa_leaf[p] + esai_p
            end

            # SP mode: calculate leaf and stem biomass
            if !use_cn
                canopystate.leaf_biomass_patch[p] = (1.0e-3 * C_TO_B / slatop_pft[ft]) *
                    smooth_max(0.01, 0.5 * sa_leaf[p]) / (1.0 - fbw_pft[ft])
                carea_stem = RPI * (dbh[p] * 0.5)^2
                canopystate.stem_biomass_patch[p] = carea_stem * canopystate.htop_patch[p] *
                    K_CYL_VOL_CANOPY * nstem_pft[ft] * wood_density_pft[ft] /
                    (1.0 - fbw_pft[ft])
            end

            # internal longwave fluxes between leaf and stem
            sa_internal[p] = min(sa_leaf[p], sa_stem[p])
            sa_internal[p] = K_INTERNAL_CANOPY * sa_internal[p]

            # heat capacity of vegetation
            cp_leaf[p] = canopystate.leaf_biomass_patch[p] *
                (C_DRY_BIOMASS * (1.0 - fbw_pft[ft]) + fbw_pft[ft] * C_WATER)
            cp_stem[p] = canopystate.stem_biomass_patch[p] *
                (C_DRY_BIOMASS * (1.0 - fbw_pft[ft]) + fbw_pft[ft] * C_WATER)
            cp_stem[p] = K_CYL_VOL_CANOPY * cp_stem[p]

            # resistance between internal stem temperature and canopy air
            rstem[p] = rstem_per_dbh_pft[ft] * dbh[p]
        end
    else
        # No biomass heat storage: set terms to zero
        for fi in 1:fn
            p = filterp[fi]
            sa_leaf[p] = canopystate.elai_patch[p] + canopystate.esai_patch[p]
            frac_rad_abs_by_stem[p] = 0.0
            sa_stem[p] = 0.0
            sa_internal[p] = 0.0
            cp_leaf[p] = 0.0
            cp_stem[p] = 0.0
            rstem[p] = 0.0
        end
    end

    # --- Daylength control for Vcmax ---
    for fi in 1:fn
        p = filterp[fi]
        g = patch_data.gridcell[p]
        dayl_factor[p] = smooth_clamp((dayl_grc[g] * dayl_grc[g]) / (max_dayl_grc[g] * max_dayl_grc[g]), 0.01, 1.0)
    end

    # Zero boundary layer resistance
    for p in bounds_patch
        frictionvel.rb1_patch[p] = 0.0
    end

    # --- Compute effective soil porosity and volumetric liquid water ---
    # (These are done inline rather than calling separate functions,
    #  since the separate module functions have different signatures)

    # --- Set perchroot options ---
    set_perchroot_opt!(ctrl.perchroot, ctrl.perchroot_alt)

    # --- Root moisture stress ---
    # (Simplified: set btran=1 and rresis=1 if detailed calc not available.
    #  The full calc_root_moist_stress call would go here.)
    # In a full implementation, this would call calc_root_moist_stress.

    # --- Modify aerodynamic parameters for sparse/dense canopy (X. Zeng) ---
    for fi in 1:fn
        p = filterp[fi]
        c = patch_data.column[p]
        g = patch_data.gridcell[p]
        ft = patch_data.itype[p] + 1  # 0-based Fortran PFT → 1-based Julia

        if z0param_method == "ZengWang2007"
            lt = smooth_min(canopystate.elai_patch[p] + canopystate.esai_patch[p], TLSAI_CRIT)
            egvf = (1.0 - ALPHA_AERO * exp(-lt)) / (1.0 - ALPHA_AERO * exp(-TLSAI_CRIT))
            canopystate.displa_patch[p] = egvf * canopystate.displa_patch[p]
            frictionvel.z0mv_patch[p] = exp(egvf * log(frictionvel.z0mv_patch[p]) +
                (1.0 - egvf) * log(frictionvel.z0mg_col[c]))

        elseif z0param_method == "Meier2022"
            lt = smooth_max(1.0e-5, canopystate.elai_patch[p] + canopystate.esai_patch[p])
            canopystate.displa_patch[p] = canopystate.htop_patch[p] *
                (1.0 - (1.0 - exp(-(CD1_PARAM * lt)^0.5)) / (CD1_PARAM * lt)^0.5)

            lt = smooth_min(lt, z0v_LAImax_pft[ft])
            delt_iter = 2.0
            U_ustar_ini = (z0v_Cs_pft[ft] + z0v_Cr_pft[ft] * lt * 0.5)^(-0.5) *
                z0v_c_pft[ft] * lt * 0.25
            U_ustar = U_ustar_ini
            while delt_iter > 1.0e-4
                U_ustar_prev = U_ustar
                U_ustar = U_ustar_ini * exp(U_ustar_prev)
                delt_iter = abs(U_ustar - U_ustar_prev)
            end
            U_ustar = 4.0 * U_ustar / lt / z0v_c_pft[ft]

            frictionvel.z0mv_patch[p] = canopystate.htop_patch[p] *
                (1.0 - canopystate.displa_patch[p] / canopystate.htop_patch[p]) *
                exp(-VKC * U_ustar + log(z0v_cw_pft[ft]) - 1.0 + z0v_cw_pft[ft]^(-1.0))
        else
            error("canopy_fluxes!: unknown z0param_method: $z0param_method")
        end

        frictionvel.z0hv_patch[p] = frictionvel.z0mv_patch[p]
        frictionvel.z0qv_patch[p] = frictionvel.z0mv_patch[p]

        # Update forcing heights
        frictionvel.forc_hgt_u_patch[p] = forc_hgt_u_grc[g] + frictionvel.z0mv_patch[p] + canopystate.displa_patch[p]
        frictionvel.forc_hgt_t_patch[p] = forc_hgt_t_grc[g] + frictionvel.z0hv_patch[p] + canopystate.displa_patch[p]
        frictionvel.forc_hgt_q_patch[p] = forc_hgt_q_grc[g] + frictionvel.z0qv_patch[p] + canopystate.displa_patch[p]
    end

    # --- Net absorbed longwave radiation, QSat, CO2/O2, flux initialization ---
    found = false
    found_index = 0
    for fi in 1:fn
        p = filterp[fi]
        c = patch_data.column[p]
        g = patch_data.gridcell[p]

        emv_p = temperature.emv_patch[p]
        emg_c = temperature.emg_col[c]

        # Net absorbed longwave = air + bir*t_veg^4 + cir*t_grnd^4
        air[p] =  emv_p * (1.0 + (1.0 - emv_p) * (1.0 - emg_c)) * forc_lwrad_col[c]
        bir[p] = -(2.0 - emv_p * (1.0 - emg_c)) * emv_p * SB
        cir[p] =  emv_p * emg_c * SB

        # Saturated vapor pressure at leaf surface
        (qs_tmp, es_tmp, dqsdT_tmp, _) = qsat(temperature.t_veg_patch[p], forc_pbot_col[c])
        qsatl[p] = qs_tmp
        el[p] = es_tmp
        qsatldT[p] = dqsdT_tmp

        # Atmospheric CO2 and O2
        co2_arr[p] = forc_pco2_grc[g]
        o2_arr[p]  = forc_po2_grc[g]

        # Initialize flux profile
        nmozsgn[p] = 0
        frictionvel.taf_patch[p] = (temperature.t_grnd_col[c] + temperature.thm_patch[p]) / 2.0
        frictionvel.qaf_patch[p] = (forc_q_col[c] + waterdiagbulk.qg_col[c]) / 2.0

        ur[p] = smooth_max(params.wind_min, sqrt(forc_u_grc[g]^2 + forc_v_grc[g]^2))
        dth[p] = temperature.thm_patch[p] - frictionvel.taf_patch[p]
        dqh[p] = forc_q_col[c] - frictionvel.qaf_patch[p]
        delq[p] = waterdiagbulk.qg_col[c] - frictionvel.qaf_patch[p]
        dthv[p] = dth[p] * (1.0 + 0.61 * forc_q_col[c]) + 0.61 * forc_th_col[c] * dqh[p]
        zldis[p] = frictionvel.forc_hgt_u_patch[p] - canopystate.displa_patch[p]

        if zldis[p] < 0.0
            found = true
            found_index = p
        end
    end

    if found && !use_fates
        @warn "canopy_fluxes!: forcing height below canopy height at patch $found_index"
    end

    # --- Initialize Monin-Obukhov length ---
    for fi in 1:fn
        p = filterp[fi]
        c = patch_data.column[p]

        (um_val, obu_val) = monin_obuk_ini(frictionvel.zetamaxstable,
            ur[p], temperature.thv_col[c], dthv[p], zldis[p], frictionvel.z0mv_patch[p])
        frictionvel.um_patch[p] = um_val
        frictionvel.obu_patch[p] = obu_val
        frictionvel.num_iter_patch[p] = 0.0

        tl_ini[p] = temperature.t_veg_patch[p]
        ts_ini[p] = temperature.t_stem_patch[p]
    end

    # --- Save original filter for post-iteration ---
    itlef = 0
    fnorig = fn
    fporig[1:fn] .= filterp[1:fn]

    # =========================================================================
    # Phase 2: Stability iteration (Newton-Raphson)
    # =========================================================================

    while itlef <= ctrl.itmax_canopy_fluxes && fn > 0

        # --- Friction velocity and boundary layer profiles ---
        # Build filter arrays for friction_velocity! call
        filt_arr = filterp[1:fn]
        disp_vec  = canopystate.displa_patch
        z0m_vec   = frictionvel.z0mv_patch
        z0h_vec   = frictionvel.z0hv_patch
        z0q_vec   = frictionvel.z0qv_patch
        obu_vec   = frictionvel.obu_patch
        um_vec    = frictionvel.um_patch
        ustar_vec = frictionvel.ustar_patch

        friction_velocity!(frictionvel, fn, filt_arr,
            disp_vec, z0m_vec, z0h_vec, z0q_vec,
            obu_vec, itlef + 1, ur, um_vec, ustar_vec,
            temp1, temp2, temp12m, temp22m, fm)

        for fi in 1:fn
            p = filterp[fi]
            c = patch_data.column[p]
            g = patch_data.gridcell[p]
            ft = patch_data.itype[p] + 1  # 0-based Fortran PFT → 1-based Julia

            tlbef[p] = temperature.t_veg_patch[p]
            del2[p] = del_arr[p]

            # Aerodynamic resistances
            frictionvel.ram1_patch[p] = 1.0 / (frictionvel.ustar_patch[p]^2 / frictionvel.um_patch[p])
            rah[p, ABOVE_CANOPY] = 1.0 / (temp1[p] * frictionvel.ustar_patch[p])
            raw[p, ABOVE_CANOPY] = 1.0 / (temp2[p] * frictionvel.ustar_patch[p])

            # Bulk boundary layer resistance
            frictionvel.uaf_patch[p] = frictionvel.um_patch[p] *
                sqrt(1.0 / (frictionvel.ram1_patch[p] * frictionvel.um_patch[p]))

            # Empirical undercanopy wind speed
            uuc[p] = smooth_min(0.4, 0.03 * frictionvel.um_patch[p] / frictionvel.ustar_patch[p])

            # Leaf characteristic width
            canopystate.dleaf_patch[p] = dleaf_pft[ft]

            cf = params.cv / (sqrt(frictionvel.uaf_patch[p]) * sqrt(canopystate.dleaf_patch[p]))
            rb[p] = 1.0 / (cf * frictionvel.uaf_patch[p])
            frictionvel.rb1_patch[p] = rb[p]

            # Soil drag coefficient (X. Zeng parameterization)
            w = exp(-(canopystate.elai_patch[p] + canopystate.esai_patch[p]))

            csoilb = VKC / (params.a_coef * (frictionvel.z0mg_col[c] * frictionvel.uaf_patch[p] / NU_PARAM)^params.a_exp)

            # Stability parameter for ricsoilc
            # Use calibration override if set, else default from params
            csoilc_val = isnan(overrides.csoilc) ? params.csoilc : overrides.csoilc

            ri = (GRAV * canopystate.htop_patch[p] *
                (frictionvel.taf_patch[p] - temperature.t_grnd_col[c])) /
                (frictionvel.taf_patch[p] * frictionvel.uaf_patch[p]^2)

            if ctrl.use_undercanopy_stability && (frictionvel.taf_patch[p] - temperature.t_grnd_col[c]) > 0.0
                ricsoilc = csoilc_val / (1.0 + RIA_CANOPY * smooth_min(ri, 10.0))
                csoilcn = csoilb * w + ricsoilc * (1.0 - w)
            else
                csoilcn = csoilb * w + csoilc_val * (1.0 - w)
            end

            if ctrl.use_biomass_heat_storage
                rah[p, BELOW_CANOPY] = 1.0 / (csoilcn * uuc[p])
            else
                rah[p, BELOW_CANOPY] = 1.0 / (csoilcn * frictionvel.uaf_patch[p])
            end

            raw[p, BELOW_CANOPY] = rah[p, BELOW_CANOPY]

            if use_lch4 && length(grnd_ch4_cond_patch) >= p
                grnd_ch4_cond_patch[p] = 1.0 / (raw[p, ABOVE_CANOPY] + raw[p, BELOW_CANOPY])
            end

            # Stomatal resistance intermediates
            svpts[p] = el[p]
            eah[p] = forc_pbot_col[c] * frictionvel.qaf_patch[p] / 0.622
            waterdiagbulk.rh_af_patch[p] = eah[p] / svpts[p]

            # History outputs
            frictionvel.rah1_patch[p] = rah[p, ABOVE_CANOPY]
            frictionvel.raw1_patch[p] = raw[p, ABOVE_CANOPY]
            frictionvel.rah2_patch[p] = rah[p, BELOW_CANOPY]
            frictionvel.raw2_patch[p] = raw[p, BELOW_CANOPY]
            frictionvel.vpd_patch[p]  = smooth_max((svpts[p] - eah[p]), 50.0) * 0.001
        end

        # --- Photosynthesis ---
        # Call photosynthesis for sunlit and shaded leaves to update rssun/rssha
        if !isempty(nrad_patch) && !isempty(parsun_z_patch)
            # Build patch-indexed forc_pbot from column-level data
            # Use full mask (not reduced filterp) since photosynthesis uses mask_exposedvegp
            forc_pbot_patch = zeros(FT, endp)
            for p in bounds_patch
                mask_exposedvegp[p] || continue
                c = patch_data.column[p]
                forc_pbot_patch[p] = forc_pbot_col[c]
            end

            # PFT index vector (+1 for 0-based Fortran → 1-based Julia)
            ivt_vec = patch_data.itype .+ 1

            # Sunlit leaves
            photosynthesis!(photosyns,
                svpts, eah, o2_arr, co2_arr, rb,
                energyflux.btran_patch, dayl_factor, leafn_patch,
                forc_pbot_patch, temperature.t_veg_patch, t10_patch,
                temperature.thm_patch, nrad_patch,
                tlai_z_patch, canopystate.tlai_patch,
                parsun_z_patch, laisun_z_patch,
                vcmaxcint_sun_patch,
                o3coefv_patch, o3coefg_patch,
                c3psn_pft, leafcn_pft, flnr_pft, fnitr_pft, slatop_pft,
                mbbopt_pft, medlynintercept_pft, medlynslope_pft,
                ivt_vec, patch_data.column,
                mask_exposedvegp, bounds_patch, "sun";
                use_cn=use_cn, use_luna=use_luna, use_c13=use_c13,
                leaf_mr_vcm=leaf_mr_vcm, crop_pft=crop_pft,
                overrides=overrides)

            # Shaded leaves
            photosynthesis!(photosyns,
                svpts, eah, o2_arr, co2_arr, rb,
                energyflux.btran_patch, dayl_factor, leafn_patch,
                forc_pbot_patch, temperature.t_veg_patch, t10_patch,
                temperature.thm_patch, nrad_patch,
                tlai_z_patch, canopystate.tlai_patch,
                parsha_z_patch, laisha_z_patch,
                vcmaxcint_sha_patch,
                o3coefv_patch, o3coefg_patch,
                c3psn_pft, leafcn_pft, flnr_pft, fnitr_pft, slatop_pft,
                mbbopt_pft, medlynintercept_pft, medlynslope_pft,
                ivt_vec, patch_data.column,
                mask_exposedvegp, bounds_patch, "sha";
                use_cn=use_cn, use_luna=use_luna, use_c13=use_c13,
                leaf_mr_vcm=leaf_mr_vcm, crop_pft=crop_pft,
                overrides=overrides)

        end

        # --- Heat transfer conductances ---
        for fi in 1:fn
            p = filterp[fi]
            c = patch_data.column[p]
            g = patch_data.gridcell[p]

            # Sensible heat conductance for air, leaf, ground, stem
            wta  = 1.0 / rah[p, ABOVE_CANOPY]       # air
            wtl  = sa_leaf[p] / rb[p]                # leaf
            wtg[p] = 1.0 / rah[p, BELOW_CANOPY]     # ground
            wtstem = sa_stem[p] / (rstem[p] + rb[p]) # stem

            wtshi = 1.0 / (wta + wtl + wtstem + wtg[p])

            wtl0[p]    = wtl * wtshi
            wtg0       = wtg[p] * wtshi
            wta0[p]    = wta * wtshi
            wtstem0[p] = wtstem * wtshi
            wtga[p]    = wta0[p] + wtg0 + wtstem0[p]
            wtal[p]    = wta0[p] + wtl0[p] + wtstem0[p]

            # Internal longwave between leaf and stem
            lw_stem[p] = sa_internal[p] * temperature.emv_patch[p] * SB * temperature.t_stem_patch[p]^4
            lw_leaf[p] = sa_internal[p] * temperature.emv_patch[p] * SB * temperature.t_veg_patch[p]^4

            # Fraction of potential evaporation from leaf
            fdry_p = waterdiagbulk.fdry_patch[p]
            elai_p = canopystate.elai_patch[p]
            esai_p = canopystate.esai_patch[p]
            laisun_p = canopystate.laisun_patch[p]
            laisha_p = canopystate.laisha_patch[p]
            rssun_p = photosyns.rssun_patch[p]
            rssha_p = photosyns.rssha_patch[p]

            if fdry_p > 0.0
                rppdry = fdry_p * rb[p] *
                    (laisun_p / (rb[p] + rssun_p) + laisha_p / (rb[p] + rssha_p)) / elai_p
            else
                rppdry = 0.0
            end

            # Canopy conductance for methane
            if use_lch4 && length(energyflux.canopy_cond_patch) >= p
                energyflux.canopy_cond_patch[p] = (laisun_p / (rb[p] + rssun_p) +
                    laisha_p / (rb[p] + rssha_p)) / smooth_max(elai_p, 0.01)
            end

            efpot = forc_rho_col[c] * ((elai_p + esai_p) / rb[p]) *
                (qsatl[p] - frictionvel.qaf_patch[p])
            h2ocan = waterstatebulk.ws.liqcan_patch[p] + waterstatebulk.ws.snocan_patch[p]

            fwet_p = waterdiagbulk.fwet_patch[p]
            btran_p = energyflux.btran_patch[p]
            qflx_tran_veg_p = waterfluxbulk.wf.qflx_tran_veg_patch[p]

            if use_hydrstress
                if efpot > 0.0
                    if btran_p > BTRAN0
                        rpp = rppdry + fwet_p
                    else
                        rpp = fwet_p
                    end
                    rpp = smooth_min(rpp, (qflx_tran_veg_p + h2ocan / dtime) / efpot)
                else
                    rpp = 1.0
                end
            else
                if efpot > 0.0
                    if btran_p > BTRAN0
                        waterfluxbulk.wf.qflx_tran_veg_patch[p] = efpot * rppdry
                        rpp = rppdry + fwet_p
                    else
                        rpp = fwet_p
                        waterfluxbulk.wf.qflx_tran_veg_patch[p] = 0.0
                    end
                    rpp = min(rpp, (waterfluxbulk.wf.qflx_tran_veg_patch[p] + h2ocan / dtime) / efpot)
                else
                    rpp = 1.0
                    waterfluxbulk.wf.qflx_tran_veg_patch[p] = 0.0
                end
            end

            # Latent heat conductances
            fvn = canopystate.frac_veg_nosno_patch[p]
            wtaq  = fvn / raw[p, ABOVE_CANOPY]
            wtlq  = fvn * (elai_p + esai_p) / rb[p] * rpp

            # Litter layer resistance (Sakaguchi)
            snow_depth_c = params.z_dl
            fsno_dl = waterdiagbulk.snow_depth_col[c] / snow_depth_c
            elai_dl = params.lai_dl * (1.0 - smooth_min(fsno_dl, 1.0))
            rdl = (1.0 - exp(-elai_dl)) / (0.004 * frictionvel.uaf_patch[p])

            if delq[p] < 0.0
                wtgq[p] = fvn / (raw[p, BELOW_CANOPY] + rdl)
            else
                if do_soilevap_beta()
                    wtgq[p] = soilstate.soilbeta_col[c] * fvn / (raw[p, BELOW_CANOPY] + rdl)
                end
                if do_soil_resistance_sl14()
                    wtgq[p] = fvn / (raw[p, BELOW_CANOPY] + soilstate.soilresis_col[c])
                end
            end

            wtsqi  = 1.0 / (wtaq + wtlq + wtgq[p])
            wtgq0  = wtgq[p] * wtsqi
            wtlq0[p] = wtlq * wtsqi
            wtaq0[p] = wtaq * wtsqi
            wtgaq  = wtaq0[p] + wtgq0
            wtalq[p] = wtaq0[p] + wtlq0[p]

            dc1 = forc_rho_col[c] * CPAIR * wtl
            dc2 = HVAP * forc_rho_col[c] * wtlq

            efsh = dc1 * (wtga[p] * temperature.t_veg_patch[p] - wtg0 * temperature.t_grnd_col[c] -
                wta0[p] * temperature.thm_patch[p] - wtstem0[p] * temperature.t_stem_patch[p])
            energyflux.eflx_sh_stem_patch[p] = forc_rho_col[c] * CPAIR * wtstem *
                ((wta0[p] + wtg0 + wtl0[p]) * temperature.t_stem_patch[p] -
                 wtg0 * temperature.t_grnd_col[c] - wta0[p] * temperature.thm_patch[p] -
                 wtl0[p] * temperature.t_veg_patch[p])
            efe[p] = dc2 * (wtgaq * qsatl[p] - wtgq0 * waterdiagbulk.qg_col[c] -
                wtaq0[p] * forc_q_col[c])

            # Evaporation flux sign change limiter
            erre = 0.0
            if efe[p] * efeb[p] < 0.0
                efeold = efe[p]
                efe[p] = 0.1 * efeold
                erre = efe[p] - efeold
            end

            # Fractionate ground emitted longwave
            snl_c = col_data.snl[c]
            lw_grnd = (waterdiagbulk.frac_sno_eff_col[c] * temperature.t_soisno_col[c, snl_c + 1 + nlevsno]^4 +
                (1.0 - waterdiagbulk.frac_sno_eff_col[c] - waterdiagbulk.frac_h2osfc_col[c]) *
                    temperature.t_soisno_col[c, 1 + nlevsno]^4 +
                waterdiagbulk.frac_h2osfc_col[c] * temperature.t_h2osfc_col[c]^4)

            # Newton-Raphson: dt_veg
            _numer = ((1.0 - frac_rad_abs_by_stem[p]) * (solarabs.sabv_patch[p] + air[p] +
                bir[p] * temperature.t_veg_patch[p]^4 + cir[p] * lw_grnd) -
                efsh - efe[p] - lw_leaf[p] + lw_stem[p] -
                (cp_leaf[p] / dtime) * (temperature.t_veg_patch[p] - tl_ini[p]))
            _denom = ((1.0 - frac_rad_abs_by_stem[p]) * (-4.0 * bir[p] * temperature.t_veg_patch[p]^3) +
                 4.0 * sa_internal[p] * temperature.emv_patch[p] * SB * temperature.t_veg_patch[p]^3 +
                 dc1 * wtga[p] + dc2 * wtgaq * qsatldT[p] + cp_leaf[p] / dtime)
            dt_veg[p] = _numer / _denom

            temperature.t_veg_patch[p] = tlbef[p] + dt_veg[p]

            dels = dt_veg[p]
            del_arr[p] = abs(dels)
            err_arr[p] = 0.0
            if del_arr[p] > DELMAX_CANOPY
                dt_veg[p] = DELMAX_CANOPY * dels / del_arr[p]
                temperature.t_veg_patch[p] = tlbef[p] + dt_veg[p]
                err_arr[p] = (1.0 - frac_rad_abs_by_stem[p]) * (solarabs.sabv_patch[p] + air[p] +
                    bir[p] * tlbef[p]^3 * (tlbef[p] + 4.0 * dt_veg[p]) + cir[p] * lw_grnd) -
                    sa_internal[p] * temperature.emv_patch[p] * SB * tlbef[p]^3 * (tlbef[p] + 4.0 * dt_veg[p]) +
                    lw_stem[p] -
                    (efsh + dc1 * wtga[p] * dt_veg[p]) -
                    (efe[p] + dc2 * wtgaq * qsatldT[p] * dt_veg[p]) -
                    (cp_leaf[p] / dtime) * (temperature.t_veg_patch[p] - tl_ini[p])
            end

            # Updated fluxes
            efpot = forc_rho_col[c] * ((elai_p + esai_p) / rb[p]) *
                (wtgaq * (qsatl[p] + qsatldT[p] * dt_veg[p]) -
                 wtgq0 * waterdiagbulk.qg_col[c] - wtaq0[p] * forc_q_col[c])
            waterfluxbulk.wf.qflx_evap_veg_patch[p] = rpp * efpot

            # Interception losses / ecidif
            if use_hydrstress
                ecidif = max(0.0, waterfluxbulk.wf.qflx_evap_veg_patch[p] -
                    waterfluxbulk.wf.qflx_tran_veg_patch[p] - h2ocan / dtime)
                waterfluxbulk.wf.qflx_evap_veg_patch[p] = min(waterfluxbulk.wf.qflx_evap_veg_patch[p],
                    waterfluxbulk.wf.qflx_tran_veg_patch[p] + h2ocan / dtime)
            else
                ecidif = 0.0
                if efpot > 0.0 && energyflux.btran_patch[p] > BTRAN0
                    waterfluxbulk.wf.qflx_tran_veg_patch[p] = efpot * rppdry
                else
                    waterfluxbulk.wf.qflx_tran_veg_patch[p] = 0.0
                end
                ecidif = max(0.0, waterfluxbulk.wf.qflx_evap_veg_patch[p] -
                    waterfluxbulk.wf.qflx_tran_veg_patch[p] - h2ocan / dtime)
                waterfluxbulk.wf.qflx_evap_veg_patch[p] = min(waterfluxbulk.wf.qflx_evap_veg_patch[p],
                    waterfluxbulk.wf.qflx_tran_veg_patch[p] + h2ocan / dtime)
            end

            # Sensible heat from leaves
            energyflux.eflx_sh_veg_patch[p] = efsh + dc1 * wtga[p] * dt_veg[p] +
                err_arr[p] + erre + HVAP * ecidif

            # Update SH and lw_leaf for changes in t_veg
            energyflux.eflx_sh_stem_patch[p] = energyflux.eflx_sh_stem_patch[p] +
                forc_rho_col[c] * CPAIR * wtstem * (-wtl0[p] * dt_veg[p])
            lw_leaf[p] = sa_internal[p] * temperature.emv_patch[p] * SB *
                tlbef[p]^3 * (tlbef[p] + 4.0 * dt_veg[p])

            # Re-calculate QSat at updated leaf temperature
            (qs_tmp, es_tmp, dqsdT_tmp, _) = qsat(temperature.t_veg_patch[p], forc_pbot_col[c])
            qsatl[p] = qs_tmp
            el[p] = es_tmp
            qsatldT[p] = dqsdT_tmp

            # Update canopy air temperature and humidity
            frictionvel.taf_patch[p] = wtg0 * temperature.t_grnd_col[c] +
                wta0[p] * temperature.thm_patch[p] +
                wtl0[p] * temperature.t_veg_patch[p] +
                wtstem0[p] * temperature.t_stem_patch[p]
            frictionvel.qaf_patch[p] = wtlq0[p] * qsatl[p] +
                wtgq0 * waterdiagbulk.qg_col[c] +
                forc_q_col[c] * wtaq0[p]

            # Update Monin-Obukhov length and wind speed
            dth[p] = temperature.thm_patch[p] - frictionvel.taf_patch[p]
            dqh[p] = forc_q_col[c] - frictionvel.qaf_patch[p]
            delq[p] = wtalq[p] * waterdiagbulk.qg_col[c] - wtlq0[p] * qsatl[p] - wtaq0[p] * forc_q_col[c]

            tstar = temp1[p] * dth[p]
            qstar = temp2[p] * dqh[p]
            thvstar = tstar * (1.0 + 0.61 * forc_q_col[c]) + 0.61 * forc_th_col[c] * qstar
            frictionvel.zeta_patch[p] = zldis[p] * VKC * GRAV * thvstar /
                (frictionvel.ustar_patch[p]^2 * temperature.thv_col[c])

            if frictionvel.zeta_patch[p] >= 0.0  # stable
                frictionvel.zeta_patch[p] = smooth_clamp(frictionvel.zeta_patch[p], 0.01, frictionvel.zetamaxstable)
                frictionvel.um_patch[p] = smooth_max(ur[p], 0.1)
            else  # unstable
                frictionvel.zeta_patch[p] = smooth_clamp(frictionvel.zeta_patch[p], -100.0, -0.01)
                if frictionvel.ustar_patch[p] * thvstar > 0.0
                    wc = 0.0
                else
                    wc_arg = smooth_max(-GRAV * frictionvel.ustar_patch[p] * thvstar *
                        ZII_CANOPY / temperature.thv_col[c], 0.0)
                    wc = BETA_CANOPY * wc_arg^0.333
                end
                frictionvel.um_patch[p] = sqrt(ur[p]^2 + wc^2)
            end
            frictionvel.obu_patch[p] = zldis[p] / frictionvel.zeta_patch[p]

            if obuold[p] * frictionvel.obu_patch[p] < 0.0
                nmozsgn[p] += 1
            end
            if nmozsgn[p] >= 4
                frictionvel.obu_patch[p] = zldis[p] / (-0.01)
            end
            obuold[p] = frictionvel.obu_patch[p]
        end  # end filtered patch loop

        # --- Test for convergence ---
        itlef += 1
        if itlef > ITMIN_CANOPY
            for fi in 1:fn
                p = filterp[fi]
                dele[p] = abs(efe[p] - efeb[p])
                efeb[p] = efe[p]
                det_arr[p] = max(del_arr[p], del2[p])
                frictionvel.num_iter_patch[p] = Float64(itlef)
            end
            fnold = fn
            fn = 0
            for fi in 1:fnold
                p = filterp[fi]
                if !(det_arr[p] < DTMIN_CANOPY && dele[p] < DLEMIN_CANOPY)
                    fn += 1
                    filterp[fn] = p
                end
            end
        end
    end  # End stability iteration

    # =========================================================================
    # Phase 3: Post-iteration diagnostics
    # =========================================================================

    # Restore original filter
    fn = fnorig
    filterp[1:fn] .= fporig[1:fn]

    for fi in 1:fn
        p = filterp[fi]
        c = patch_data.column[p]
        g = patch_data.gridcell[p]

        # Energy balance check
        snl_c = col_data.snl[c]
        lw_grnd = (waterdiagbulk.frac_sno_eff_col[c] * temperature.t_soisno_col[c, snl_c + 1 + nlevsno]^4 +
            (1.0 - waterdiagbulk.frac_sno_eff_col[c] - waterdiagbulk.frac_h2osfc_col[c]) *
                temperature.t_soisno_col[c, 1 + nlevsno]^4 +
            waterdiagbulk.frac_h2osfc_col[c] * temperature.t_h2osfc_col[c]^4)

        err_arr[p] = (1.0 - frac_rad_abs_by_stem[p]) * (solarabs.sabv_patch[p] + air[p] +
            bir[p] * tlbef[p]^3 * (tlbef[p] + 4.0 * dt_veg[p]) + cir[p] * lw_grnd) -
            lw_leaf[p] + lw_stem[p] - energyflux.eflx_sh_veg_patch[p] -
            HVAP * waterfluxbulk.wf.qflx_evap_veg_patch[p] -
            ((temperature.t_veg_patch[p] - tl_ini[p]) * cp_leaf[p] / dtime)

        # Update stem temperature
        if ctrl.use_biomass_heat_storage
            if canopystate.stem_biomass_patch[p] > 0.0
                dt_stem[p] = (frac_rad_abs_by_stem[p] * (solarabs.sabv_patch[p] + air[p] +
                    bir[p] * ts_ini[p]^4 + cir[p] * lw_grnd) -
                    energyflux.eflx_sh_stem_patch[p] +
                    lw_leaf[p] - lw_stem[p]) /
                    (cp_stem[p] / dtime - frac_rad_abs_by_stem[p] * bir[p] * 4.0 * ts_ini[p]^3)
            else
                dt_stem[p] = 0.0
            end

            energyflux.dhsdt_canopy_patch[p] = dt_stem[p] * cp_stem[p] / dtime +
                (temperature.t_veg_patch[p] - tl_ini[p]) * cp_leaf[p] / dtime
            temperature.t_stem_patch[p] = temperature.t_stem_patch[p] + dt_stem[p]
        else
            dt_stem[p] = 0.0
        end

        delt_val = wtal[p] * temperature.t_grnd_col[c] - wtl0[p] * temperature.t_veg_patch[p] -
            wta0[p] * temperature.thm_patch[p] - wtstem0[p] * temperature.t_stem_patch[p]

        # Ground fluxes
        energyflux.taux_patch[p] = -forc_rho_col[c] * forc_u_grc[g] / frictionvel.ram1_patch[p]
        energyflux.tauy_patch[p] = -forc_rho_col[c] * forc_v_grc[g] / frictionvel.ram1_patch[p]
        energyflux.eflx_sh_grnd_patch[p] = CPAIR * forc_rho_col[c] * wtg[p] * delt_val

        # Individual sensible heat fluxes
        delt_snow = wtal[p] * temperature.t_soisno_col[c, snl_c + 1 + nlevsno] -
            wtl0[p] * temperature.t_veg_patch[p] - wta0[p] * temperature.thm_patch[p] -
            wtstem0[p] * temperature.t_stem_patch[p]
        delt_soil = wtal[p] * temperature.t_soisno_col[c, 1 + nlevsno] -
            wtl0[p] * temperature.t_veg_patch[p] - wta0[p] * temperature.thm_patch[p] -
            wtstem0[p] * temperature.t_stem_patch[p]
        delt_h2osfc = wtal[p] * temperature.t_h2osfc_col[c] -
            wtl0[p] * temperature.t_veg_patch[p] - wta0[p] * temperature.thm_patch[p] -
            wtstem0[p] * temperature.t_stem_patch[p]

        energyflux.eflx_sh_snow_patch[p] = CPAIR * forc_rho_col[c] * wtg[p] * delt_snow
        energyflux.eflx_sh_soil_patch[p] = CPAIR * forc_rho_col[c] * wtg[p] * delt_soil
        energyflux.eflx_sh_h2osfc_patch[p] = CPAIR * forc_rho_col[c] * wtg[p] * delt_h2osfc
        waterfluxbulk.wf.qflx_evap_soi_patch[p] = forc_rho_col[c] * wtgq[p] * delq[p]

        # Individual latent heat fluxes
        delq_snow = wtalq[p] * waterdiagbulk.qg_snow_col[c] - wtlq0[p] * qsatl[p] - wtaq0[p] * forc_q_col[c]
        waterfluxbulk.qflx_ev_snow_patch[p] = forc_rho_col[c] * wtgq[p] * delq_snow

        delq_soil = wtalq[p] * waterdiagbulk.qg_soil_col[c] - wtlq0[p] * qsatl[p] - wtaq0[p] * forc_q_col[c]
        waterfluxbulk.qflx_ev_soil_patch[p] = forc_rho_col[c] * wtgq[p] * delq_soil

        delq_h2osfc = wtalq[p] * waterdiagbulk.qg_h2osfc_col[c] - wtlq0[p] * qsatl[p] - wtaq0[p] * forc_q_col[c]
        waterfluxbulk.qflx_ev_h2osfc_patch[p] = forc_rho_col[c] * wtgq[p] * delq_h2osfc

        # 2m reference height temperature
        temperature.t_ref2m_patch[p] = temperature.thm_patch[p] +
            temp1[p] * dth[p] * (1.0 / temp12m[p] - 1.0 / temp1[p])
        temperature.t_ref2m_r_patch[p] = temperature.t_ref2m_patch[p]

        # 2m specific humidity
        waterdiagbulk.q_ref2m_patch[p] = forc_q_col[c] +
            temp2[p] * dqh[p] * (1.0 / temp22m[p] - 1.0 / temp2[p])

        # 2m relative humidity
        (qsat_ref2m, e_ref2m, _, _) = qsat(temperature.t_ref2m_patch[p], forc_pbot_col[c])
        waterdiagbulk.rh_ref2m_patch[p] = smooth_min(100.0,
            waterdiagbulk.q_ref2m_patch[p] / qsat_ref2m * 100.0)
        waterdiagbulk.rh_ref2m_r_patch[p] = waterdiagbulk.rh_ref2m_patch[p]

        # 2m vapor pressure deficit
        waterdiagbulk.vpd_ref2m_patch[p] = e_ref2m * (1.0 - waterdiagbulk.rh_ref2m_patch[p] / 100.0)

        # Human heat stress indices (stub — HumanIndexMod not yet ported)
        # Would be called here in full implementation

        # Downward longwave below canopy
        emv_p = temperature.emv_patch[p]
        emg_c = temperature.emg_col[c]
        energyflux.dlrad_patch[p] = (1.0 - emv_p) * emg_c * forc_lwrad_col[c] +
            emv_p * emg_c * SB * tlbef[p]^3 * (tlbef[p] + 4.0 * dt_veg[p]) *
            (1.0 - frac_rad_abs_by_stem[p]) +
            emv_p * emg_c * SB * ts_ini[p]^3 * (ts_ini[p] + 4.0 * dt_stem[p]) *
            frac_rad_abs_by_stem[p]

        # Upward longwave above canopy
        energyflux.ulrad_patch[p] = ((1.0 - emg_c) * (1.0 - emv_p) * (1.0 - emv_p) * forc_lwrad_col[c] +
            emv_p * (1.0 + (1.0 - emg_c) * (1.0 - emv_p)) * SB *
            tlbef[p]^3 * (tlbef[p] + 4.0 * dt_veg[p]) * (1.0 - frac_rad_abs_by_stem[p]) +
            emv_p * (1.0 + (1.0 - emg_c) * (1.0 - emv_p)) * SB *
            ts_ini[p]^3 * (ts_ini[p] + 4.0 * dt_stem[p]) * frac_rad_abs_by_stem[p] +
            emg_c * (1.0 - emv_p) * SB * lw_grnd)

        # Skin temperature
        temperature.t_skin_patch[p] = emv_p * temperature.t_veg_patch[p] +
            (1.0 - emv_p) * sqrt(sqrt(lw_grnd))

        # Derivative of soil energy flux
        energyflux.cgrnds_patch[p] = energyflux.cgrnds_patch[p] +
            CPAIR * forc_rho_col[c] * wtg[p] * wtal[p]
        energyflux.cgrndl_patch[p] = energyflux.cgrndl_patch[p] +
            forc_rho_col[c] * wtgq[p] * wtalq[p] * waterdiagbulk.dqgdT_col[c]
        energyflux.cgrnd_patch[p] = energyflux.cgrnds_patch[p] +
            energyflux.cgrndl_patch[p] * energyflux.htvp_col[c]

        # Save baseline snocan
        snocan_baseline[p] = waterstatebulk.ws.snocan_patch[p]

        # Update dew accumulation
        t_veg_p = temperature.t_veg_patch[p]
        qflx_evap_veg_p = waterfluxbulk.wf.qflx_evap_veg_patch[p]
        qflx_tran_veg_p = waterfluxbulk.wf.qflx_tran_veg_patch[p]

        # Smooth freeze/thaw partitioning for AD
        _frac_liq = smooth_heaviside(t_veg_p - TFRZ; k=200.0)
        _frac_ice = one(t_veg_p) - _frac_liq
        _net_evap = (qflx_tran_veg_p - qflx_evap_veg_p) * dtime
        waterstatebulk.ws.liqcan_patch[p] = smooth_max(zero(t_veg_p),
            waterstatebulk.ws.liqcan_patch[p] + _frac_liq * _net_evap)
        waterstatebulk.ws.snocan_patch[p] = smooth_max(zero(t_veg_p),
            waterstatebulk.ws.snocan_patch[p] + _frac_ice * _net_evap)
    end

    # --- Post-photosynthesis: PhotosynthesisTotal ---
    photosynthesis_total!(photosyns,
        canopystate.laisun_patch, canopystate.laisha_patch,
        mask_exposedvegp, bounds_patch)

    # Water use efficiency (iwue) — stub for local noon check
    # LUNA and ozone updates would go here in full implementation

    # Filter out patches with small energy balance errors, report large ones
    fnold = fn
    fn = 0
    for fi in 1:fnold
        p = filterp[fi]
        if abs(err_arr[p]) > 0.1
            fn += 1
            filterp[fn] = p
        end
    end
    for fi in 1:fn
        p = filterp[fi]
        @warn "canopy_fluxes!: energy balance error at patch $p, err=$(err_arr[p])"
    end

    return nothing
end
