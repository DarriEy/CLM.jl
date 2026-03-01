# ==========================================================================
# Ported from: src/biogeophys/FrictionVelocityMod.F90
# Friction velocity data type, stability functions, Monin-Obukhov
# initialization, roughness length setting, and friction velocity calculation.
# ==========================================================================

"""
    FrictionVelocityData

Friction velocity data structure. Holds scalar parameters (roughness lengths,
stability limit) and per-patch / per-column arrays for friction velocity
computations, roughness lengths, forcing heights, and diagnostic variables.

Ported from `frictionvel_type` in `FrictionVelocityMod.F90`.
"""
Base.@kwdef mutable struct FrictionVelocityData
    # --- Scalar parameters ---
    zetamaxstable::Float64 = -999.0   # Max zeta under stable conditions
    zsno::Float64          = -999.0   # Momentum roughness length for snow (m)
    zlnd::Float64          = -999.0   # Momentum roughness length for soil, glacier, wetland (m)
    zglc::Float64          = -999.0   # Momentum roughness length for glacier (Meier2022 only) (m)

    # --- Patch-level fields (1D vectors) ---
    forc_hgt_u_patch ::Vector{Float64} = Float64[]   # patch wind forcing height (10m+z0m+d) (m)
    forc_hgt_t_patch ::Vector{Float64} = Float64[]   # patch temperature forcing height (10m+z0m+d) (m)
    forc_hgt_q_patch ::Vector{Float64} = Float64[]   # patch specific humidity forcing height (10m+z0m+d) (m)
    u10_patch        ::Vector{Float64} = Float64[]   # patch 10-m wind (m/s) (for dust model)
    u10_clm_patch    ::Vector{Float64} = Float64[]   # patch 10-m wind (m/s) (for clm_map2gcell)
    va_patch         ::Vector{Float64} = Float64[]   # patch atmospheric wind speed plus convective velocity (m/s)
    vds_patch        ::Vector{Float64} = Float64[]   # patch deposition velocity term (m/s)
    fv_patch         ::Vector{Float64} = Float64[]   # patch friction velocity (m/s) (for dust model)
    rb1_patch        ::Vector{Float64} = Float64[]   # patch aerodynamical resistance (s/m)
    rb10_patch       ::Vector{Float64} = Float64[]   # 10-day mean patch aerodynamical resistance (s/m)
    ram1_patch       ::Vector{Float64} = Float64[]   # patch aerodynamical resistance (s/m)
    z0mv_patch       ::Vector{Float64} = Float64[]   # patch roughness length over vegetation, momentum [m]
    z0hv_patch       ::Vector{Float64} = Float64[]   # patch roughness length over vegetation, sensible heat [m]
    z0qv_patch       ::Vector{Float64} = Float64[]   # patch roughness length over vegetation, latent heat [m]
    z0mg_patch       ::Vector{Float64} = Float64[]   # patch roughness length over ground, momentum [m]
    z0hg_patch       ::Vector{Float64} = Float64[]   # patch roughness length over ground, sensible heat [m]
    z0qg_patch       ::Vector{Float64} = Float64[]   # patch roughness length over ground, latent heat [m]
    kbm1_patch       ::Vector{Float64} = Float64[]   # natural logarithm of z0mg_p/z0hg_p [-]
    rah1_patch       ::Vector{Float64} = Float64[]   # patch sensible heat flux resistance [s/m]
    rah2_patch       ::Vector{Float64} = Float64[]   # patch below-canopy sensible heat flux resistance [s/m]
    raw1_patch       ::Vector{Float64} = Float64[]   # patch moisture flux resistance [s/m]
    raw2_patch       ::Vector{Float64} = Float64[]   # patch below-canopy moisture flux resistance [s/m]
    ustar_patch      ::Vector{Float64} = Float64[]   # patch friction velocity [m/s]
    um_patch         ::Vector{Float64} = Float64[]   # patch wind speed including the stability effect [m/s]
    uaf_patch        ::Vector{Float64} = Float64[]   # patch canopy air speed [m/s]
    taf_patch        ::Vector{Float64} = Float64[]   # patch canopy air temperature [K]
    qaf_patch        ::Vector{Float64} = Float64[]   # patch canopy humidity [kg/kg]
    obu_patch        ::Vector{Float64} = Float64[]   # patch Monin-Obukhov length [m]
    zeta_patch       ::Vector{Float64} = Float64[]   # patch dimensionless stability parameter
    vpd_patch        ::Vector{Float64} = Float64[]   # patch vapor pressure deficit [Pa]
    num_iter_patch   ::Vector{Float64} = Float64[]   # patch number of iterations
    z0m_actual_patch ::Vector{Float64} = Float64[]   # patch roughness length actually used, momentum [m]

    # --- Column-level fields (1D vectors) ---
    z0mg_col    ::Vector{Float64} = Float64[]   # col roughness length over ground, momentum [m]
    z0mg_2D_col ::Vector{Float64} = Float64[]   # 2-D field of input col roughness length over ground, momentum [m]
    z0hg_col    ::Vector{Float64} = Float64[]   # col roughness length over ground, sensible heat [m]
    z0qg_col    ::Vector{Float64} = Float64[]   # col roughness length over ground, latent heat [m]
end

# ==========================================================================
# Allocation / initialization
# ==========================================================================

"""
    frictionvel_init!(fv::FrictionVelocityData, np::Int, nc::Int)

Allocate and initialize all fields of a `FrictionVelocityData` instance for
`np` patches and `nc` columns. Float fields initialized to NaN, except
rb10_patch which is initialized to SPVAL (matching Fortran).

Ported from `frictionvel_type%InitAllocate` in `FrictionVelocityMod.F90`.
"""
function frictionvel_init!(fv::FrictionVelocityData, np::Int, nc::Int)
    # Patch-level fields
    fv.forc_hgt_u_patch  = fill(NaN, np)
    fv.forc_hgt_t_patch  = fill(NaN, np)
    fv.forc_hgt_q_patch  = fill(NaN, np)
    fv.u10_patch         = fill(NaN, np)
    fv.u10_clm_patch     = fill(NaN, np)
    fv.va_patch          = fill(NaN, np)
    fv.vds_patch         = fill(NaN, np)
    fv.fv_patch          = fill(NaN, np)
    fv.rb1_patch         = fill(NaN, np)
    fv.rb10_patch        = fill(SPVAL, np)
    fv.ram1_patch        = fill(NaN, np)
    fv.z0mv_patch        = fill(NaN, np)
    fv.z0hv_patch        = fill(NaN, np)
    fv.z0qv_patch        = fill(NaN, np)
    fv.z0mg_patch        = fill(NaN, np)
    fv.z0hg_patch        = fill(NaN, np)
    fv.z0qg_patch        = fill(NaN, np)
    fv.kbm1_patch        = fill(NaN, np)
    fv.rah1_patch        = fill(NaN, np)
    fv.rah2_patch        = fill(NaN, np)
    fv.raw1_patch        = fill(NaN, np)
    fv.raw2_patch        = fill(NaN, np)
    fv.ustar_patch       = fill(NaN, np)
    fv.um_patch          = fill(NaN, np)
    fv.uaf_patch         = fill(NaN, np)
    fv.taf_patch         = fill(NaN, np)
    fv.qaf_patch         = fill(NaN, np)
    fv.obu_patch         = fill(NaN, np)
    fv.zeta_patch        = fill(NaN, np)
    fv.vpd_patch         = fill(NaN, np)
    fv.num_iter_patch    = fill(NaN, np)
    fv.z0m_actual_patch  = fill(NaN, np)

    # Column-level fields
    fv.z0mg_col    = fill(NaN, nc)
    fv.z0mg_2D_col = fill(NaN, nc)
    fv.z0hg_col    = fill(NaN, nc)
    fv.z0qg_col    = fill(NaN, nc)

    return nothing
end

"""
    frictionvel_clean!(fv::FrictionVelocityData)

Deallocate (reset to empty) all fields of a `FrictionVelocityData` instance.
"""
function frictionvel_clean!(fv::FrictionVelocityData)
    fv.forc_hgt_u_patch  = Float64[]
    fv.forc_hgt_t_patch  = Float64[]
    fv.forc_hgt_q_patch  = Float64[]
    fv.u10_patch         = Float64[]
    fv.u10_clm_patch     = Float64[]
    fv.va_patch          = Float64[]
    fv.vds_patch         = Float64[]
    fv.fv_patch          = Float64[]
    fv.rb1_patch         = Float64[]
    fv.rb10_patch        = Float64[]
    fv.ram1_patch        = Float64[]
    fv.z0mv_patch        = Float64[]
    fv.z0hv_patch        = Float64[]
    fv.z0qv_patch        = Float64[]
    fv.z0mg_patch        = Float64[]
    fv.z0hg_patch        = Float64[]
    fv.z0qg_patch        = Float64[]
    fv.kbm1_patch        = Float64[]
    fv.rah1_patch        = Float64[]
    fv.rah2_patch        = Float64[]
    fv.raw1_patch        = Float64[]
    fv.raw2_patch        = Float64[]
    fv.ustar_patch       = Float64[]
    fv.um_patch          = Float64[]
    fv.uaf_patch         = Float64[]
    fv.taf_patch         = Float64[]
    fv.qaf_patch         = Float64[]
    fv.obu_patch         = Float64[]
    fv.zeta_patch        = Float64[]
    fv.vpd_patch         = Float64[]
    fv.num_iter_patch    = Float64[]
    fv.z0m_actual_patch  = Float64[]
    fv.z0mg_col          = Float64[]
    fv.z0mg_2D_col       = Float64[]
    fv.z0hg_col          = Float64[]
    fv.z0qg_col          = Float64[]
    return nothing
end

"""
    frictionvel_init_for_testing!(fv::FrictionVelocityData, np::Int, nc::Int)

Initialize for unit testing with hardcoded namelist and parameter values.

Ported from `frictionvel_type%InitForTesting` in `FrictionVelocityMod.F90`.
"""
function frictionvel_init_for_testing!(fv::FrictionVelocityData, np::Int, nc::Int)
    frictionvel_init!(fv, np, nc)
    fv.zetamaxstable = 0.5
    fv.zsno = 0.00085
    fv.zlnd = 0.000775
    fv.zglc = 0.00230000005
    return nothing
end

"""
    frictionvel_init_cold!(fv::FrictionVelocityData, bounds_col::UnitRange{Int},
                           bounds_patch::UnitRange{Int};
                           col_landunit::Union{Vector{Int},Nothing}=nothing,
                           lun_lakpoi::Union{Vector{Bool},Nothing}=nothing)

Cold-start initialization. Sets forc_hgt_u_patch to 30m for CN mode,
and z0mg_col to 0.0004 for lake columns.

Ported from `frictionvel_type%InitCold` in `FrictionVelocityMod.F90`.
"""
function frictionvel_init_cold!(fv::FrictionVelocityData,
                                 bounds_col::UnitRange{Int},
                                 bounds_patch::UnitRange{Int};
                                 col_landunit::Union{Vector{Int},Nothing} = nothing,
                                 lun_lakpoi::Union{Vector{Bool},Nothing} = nothing,
                                 use_cn::Bool = false)
    if use_cn
        for p in bounds_patch
            fv.forc_hgt_u_patch[p] = 30.0
        end
    end

    if col_landunit !== nothing && lun_lakpoi !== nothing
        for c in bounds_col
            l = col_landunit[c]
            if lun_lakpoi[l]
                fv.z0mg_col[c] = 0.0004
            end
        end
    end

    return nothing
end

# ==========================================================================
# History / Restart stubs
# ==========================================================================

"""
    frictionvel_init_history!(fv::FrictionVelocityData, bounds_patch::UnitRange{Int},
                               bounds_col::UnitRange{Int})

Register friction velocity fields for history output.
Stub until histFileMod is ported.

Ported from `frictionvel_type%InitHistory` in `FrictionVelocityMod.F90`.
"""
function frictionvel_init_history!(fv::FrictionVelocityData,
                                    bounds_patch::UnitRange{Int},
                                    bounds_col::UnitRange{Int})
    return nothing
end

"""
    frictionvel_restart!(fv::FrictionVelocityData, bounds::UnitRange{Int};
                          flag::String="read")

Read/write friction velocity restart variables.
Stub until restart infrastructure is ported.
Variables: Z0MG, OBU, rb10.

Ported from `frictionvel_type%Restart` in `FrictionVelocityMod.F90`.
"""
function frictionvel_restart!(fv::FrictionVelocityData,
                               bounds::UnitRange{Int};
                               flag::String = "read")
    return nothing
end

# ==========================================================================
# ReadNamelist — Set namelist parameters
# ==========================================================================

"""
    frictionvel_read_nml!(fv::FrictionVelocityData; zetamaxstable::Float64=0.5)

Set namelist parameters for friction velocity.
In Julia, namelist values are passed directly instead of reading from a file.

Ported from `frictionvel_type%ReadNamelist` in `FrictionVelocityMod.F90`.
"""
function frictionvel_read_nml!(fv::FrictionVelocityData;
                                zetamaxstable::Float64 = 0.5)
    fv.zetamaxstable = zetamaxstable
    return nothing
end

# ==========================================================================
# ReadParams — Read parameters from file (stub)
# ==========================================================================

"""
    frictionvel_read_params!(fv::FrictionVelocityData;
                              zsno::Float64=0.00085,
                              zlnd::Float64=0.000775,
                              zglc::Union{Float64,Nothing}=nothing,
                              z0param_method::String="")

Set roughness length parameters. In Julia, values are passed directly
instead of reading from a NetCDF parameter file.

Parameters:
- `zsno`: Momentum roughness length for snow (m)
- `zlnd`: Momentum roughness length for soil, glacier, wetland (m)
- `zglc`: Momentum roughness length for glacier (m), only used with Meier2022

Ported from `frictionvel_type%ReadParams` in `FrictionVelocityMod.F90`.
"""
function frictionvel_read_params!(fv::FrictionVelocityData;
                                   zsno::Float64 = 0.00085,
                                   zlnd::Float64 = 0.000775,
                                   zglc::Union{Float64,Nothing} = nothing,
                                   z0param_method::String = "")
    fv.zsno = zsno
    fv.zlnd = zlnd
    if z0param_method == "Meier2022" && zglc !== nothing
        fv.zglc = zglc
    end
    return nothing
end

# ==========================================================================
# Stability functions (pure math, no global state)
# ==========================================================================

"""
    stability_func1(zeta::Float64) -> Float64

Stability function for rib < 0 (wind profile).
Ported from `StabilityFunc1` in `FrictionVelocityMod.F90`.
"""
function stability_func1(zeta::Float64)
    chik2 = sqrt(1.0 - 16.0 * zeta)
    chik = sqrt(chik2)
    return 2.0 * log((1.0 + chik) * 0.5) +
           log((1.0 + chik2) * 0.5) -
           2.0 * atan(chik) + π * 0.5
end

"""
    stability_func2(zeta::Float64) -> Float64

Stability function for rib < 0 (temperature/humidity profile).
Ported from `StabilityFunc2` in `FrictionVelocityMod.F90`.
"""
function stability_func2(zeta::Float64)
    chik2 = sqrt(1.0 - 16.0 * zeta)
    return 2.0 * log((1.0 + chik2) * 0.5)
end

# ==========================================================================
# MoninObukIni — Monin-Obukhov length initialization
# ==========================================================================

"""
    monin_obuk_ini(zetamaxstable::Float64, ur::Float64, thv::Float64,
                   dthv::Float64, zldis::Float64, z0m::Float64)
                   -> (um::Float64, obu::Float64)

Initialize the Monin-Obukhov length. Based on Zeng et al. (1998).

Returns `(um, obu)` — wind speed including stability effect and
the Monin-Obukhov length.

Ported from `frictionvel_type%MoninObukIni` in `FrictionVelocityMod.F90`.
"""
function monin_obuk_ini(zetamaxstable::Float64, ur::Float64, thv::Float64,
                        dthv::Float64, zldis::Float64, z0m::Float64)
    ustar = 0.06
    wc = 0.5
    if dthv >= 0.0
        um = max(ur, 0.1)
    else
        um = sqrt(ur * ur + wc * wc)
    end

    rib = GRAV * zldis * dthv / (thv * um * um)

    if rib >= 0.0  # neutral or stable
        zeta = rib * log(zldis / z0m) / (1.0 - 5.0 * min(rib, 0.19))
        zeta = min(zetamaxstable, max(zeta, 0.01))
    else           # unstable
        zeta = rib * log(zldis / z0m)
        zeta = max(-100.0, min(zeta, -0.01))
    end

    obu = zldis / zeta

    return (um, obu)
end

# ==========================================================================
# SetRoughnessLengthsAndForcHeightsNonLake
# ==========================================================================

"""
    set_roughness_and_forc_heights_nonlake!(
        fv, mask_nolakec, mask_nolakep, bounds_col, bounds_patch,
        frac_sno, snomelt_accum, frac_veg_nosno, z0m, displa,
        forc_hgt_u, forc_hgt_t, forc_hgt_q,
        col_landunit, patch_gridcell, patch_landunit, patch_column,
        lun_itype, lun_urbpoi, lun_z_0_town, lun_z_d_town;
        z0param_method, use_z0m_snowmelt)

Set roughness lengths and forcing heights for non-lake points.

Ported from `frictionvel_type%SetRoughnessLengthsAndForcHeightsNonLake`
in `FrictionVelocityMod.F90`.
"""
function set_roughness_and_forc_heights_nonlake!(
        fv::FrictionVelocityData,
        mask_nolakec::BitVector,
        mask_nolakep::BitVector,
        bounds_col::UnitRange{Int},
        bounds_patch::UnitRange{Int},
        frac_sno::Vector{Float64},
        snomelt_accum::Vector{Float64},
        frac_veg_nosno::Vector{Int},
        z0m::Vector{Float64},
        displa::Vector{Float64},
        forc_hgt_u::Vector{Float64},
        forc_hgt_t::Vector{Float64},
        forc_hgt_q::Vector{Float64},
        col_landunit::Vector{Int},
        patch_gridcell::Vector{Int},
        patch_landunit::Vector{Int},
        patch_column::Vector{Int},
        lun_itype::Vector{Int},
        lun_urbpoi::Vector{Bool},
        lun_z_0_town::Vector{Float64},
        lun_z_d_town::Vector{Float64};
        z0param_method::String = "",
        use_z0m_snowmelt::Bool = false)

    # Column loop: set ground roughness lengths
    for c in bounds_col
        mask_nolakec[c] || continue

        if z0param_method == "ZengWang2007"
            if frac_sno[c] > 0.0
                fv.z0mg_col[c] = fv.zsno
            else
                fv.z0mg_col[c] = fv.zlnd
            end
        elseif z0param_method == "Meier2022"
            l = col_landunit[c]
            if frac_sno[c] > 0.0
                if use_z0m_snowmelt
                    if snomelt_accum[c] < 1.0e-5
                        fv.z0mg_col[c] = exp(-B1_PARAM * RPI * 0.5 + B4_PARAM) * 1.0e-3
                    else
                        fv.z0mg_col[c] = exp(B1_PARAM * (atan((log10(snomelt_accum[c]) + MEIER_PARAM1) / MEIER_PARAM2)) + B4_PARAM) * 1.0e-3
                    end
                else
                    fv.z0mg_col[c] = fv.zsno
                end
            elseif lun_itype[l] == ISTICE
                fv.z0mg_col[c] = fv.zglc
            else
                fv.z0mg_col[c] = fv.zlnd
            end
        end

        fv.z0hg_col[c] = fv.z0mg_col[c]  # initial set only
        fv.z0qg_col[c] = fv.z0mg_col[c]  # initial set only
    end

    # Patch loop: set vegetation roughness and initial ground roughness at patch level
    for p in bounds_patch
        mask_nolakep[p] || continue

        # Roughness lengths over vegetation
        fv.z0mv_patch[p] = z0m[p]
        fv.z0hv_patch[p] = fv.z0mv_patch[p]
        fv.z0qv_patch[p] = fv.z0mv_patch[p]

        # Set to arbitrary value (will be overwritten by respective modules)
        fv.z0mg_patch[p] = SPVAL
        fv.z0hg_patch[p] = SPVAL
        fv.z0qg_patch[p] = SPVAL
        fv.kbm1_patch[p] = SPVAL
    end

    # Patch loop: forcing height = atmospheric forcing height + z0m + displa
    for p in bounds_patch
        mask_nolakep[p] || continue
        g = patch_gridcell[p]
        l = patch_landunit[p]
        c = patch_column[p]

        if lun_itype[l] == ISTSOIL || lun_itype[l] == ISTCROP
            if frac_veg_nosno[p] == 0
                fv.forc_hgt_u_patch[p] = forc_hgt_u[g] + fv.z0mg_col[c] + displa[p]
                fv.forc_hgt_t_patch[p] = forc_hgt_t[g] + fv.z0hg_col[c] + displa[p]
                fv.forc_hgt_q_patch[p] = forc_hgt_q[g] + fv.z0qg_col[c] + displa[p]
            else
                fv.forc_hgt_u_patch[p] = forc_hgt_u[g] + fv.z0mv_patch[p] + displa[p]
                fv.forc_hgt_t_patch[p] = forc_hgt_t[g] + fv.z0hv_patch[p] + displa[p]
                fv.forc_hgt_q_patch[p] = forc_hgt_q[g] + fv.z0qv_patch[p] + displa[p]
            end
        elseif lun_itype[l] == ISTWET || lun_itype[l] == ISTICE
            fv.forc_hgt_u_patch[p] = forc_hgt_u[g] + fv.z0mg_col[c] + displa[p]
            fv.forc_hgt_t_patch[p] = forc_hgt_t[g] + fv.z0hg_col[c] + displa[p]
            fv.forc_hgt_q_patch[p] = forc_hgt_q[g] + fv.z0qg_col[c] + displa[p]
        elseif lun_urbpoi[l]
            fv.forc_hgt_u_patch[p] = forc_hgt_u[g] + lun_z_0_town[l] + lun_z_d_town[l]
            fv.forc_hgt_t_patch[p] = forc_hgt_t[g] + lun_z_0_town[l] + lun_z_d_town[l]
            fv.forc_hgt_q_patch[p] = forc_hgt_q[g] + lun_z_0_town[l] + lun_z_d_town[l]
        end
    end

    return nothing
end

# ==========================================================================
# SetActualRoughnessLengths
# ==========================================================================

"""
    set_actual_roughness_lengths!(
        fv, mask_exposedvegp, mask_noexposedvegp, mask_urbanp, mask_lakep,
        bounds_patch, patch_column, patch_landunit, lun_z_0_town)

Set the roughness lengths actually used in flux calculations.

Ported from `frictionvel_type%SetActualRoughnessLengths`
in `FrictionVelocityMod.F90`.
"""
function set_actual_roughness_lengths!(
        fv::FrictionVelocityData,
        mask_exposedvegp::BitVector,
        mask_noexposedvegp::BitVector,
        mask_urbanp::BitVector,
        mask_lakep::BitVector,
        bounds_patch::UnitRange{Int},
        patch_column::Vector{Int},
        patch_landunit::Vector{Int},
        lun_z_0_town::Vector{Float64})

    for p in bounds_patch
        if mask_exposedvegp[p]
            fv.z0m_actual_patch[p] = fv.z0mv_patch[p]
        elseif mask_noexposedvegp[p]
            c = patch_column[p]
            fv.z0m_actual_patch[p] = fv.z0mg_col[c]
        elseif mask_urbanp[p]
            l = patch_landunit[p]
            fv.z0m_actual_patch[p] = lun_z_0_town[l]
        elseif mask_lakep[p]
            c = patch_column[p]
            fv.z0m_actual_patch[p] = fv.z0mg_col[c]
        end
    end

    return nothing
end

# ==========================================================================
# FrictionVelocity — main friction velocity calculation
# ==========================================================================

"""
    friction_velocity!(fv, fn, filtern, displa, z0m, z0h, z0q, obu,
                       iter, ur, um, ustar, temp1, temp2, temp12m,
                       temp22m, fm;
                       landunit_index=false,
                       lun_gridcell=nothing, lun_patchi=nothing,
                       lun_patchf=nothing,
                       patch_gridcell=nothing)

Calculate friction velocity, and profiles for potential temperature and
humidity in the surface boundary layer. Based on Zeng et al. (1998).

When `landunit_index=true`, indices in filtern refer to landunits (and
results are broadcast to constituent patches). Otherwise they refer
to patches directly.

All array arguments are indexed from 1.

Ported from `frictionvel_type%FrictionVelocity` in `FrictionVelocityMod.F90`.
"""
function friction_velocity!(
        fv::FrictionVelocityData,
        fn::Int,
        filtern::Vector{Int},
        displa::Vector{Float64},
        z0m::Vector{Float64},
        z0h::Vector{Float64},
        z0q::Vector{Float64},
        obu::Vector{Float64},
        iter::Int,
        ur::Vector{Float64},
        um::Vector{Float64},
        ustar::Vector{Float64},
        temp1::Vector{Float64},
        temp2::Vector{Float64},
        temp12m::Vector{Float64},
        temp22m::Vector{Float64},
        fm::Vector{Float64};
        landunit_index::Bool = false,
        lun_gridcell::Union{Vector{Int},Nothing} = nothing,
        lun_patchi::Union{Vector{Int},Nothing} = nothing,
        lun_patchf::Union{Vector{Int},Nothing} = nothing,
        patch_gridcell::Union{Vector{Int},Nothing} = nothing)

    zetam = 1.574  # transition point of flux-gradient relation (wind profile)
    zetat = 0.465  # transition point of flux-gradient relation (temp. profile)

    for f in 1:fn
        n = filtern[f]

        # ------- Wind profile -------
        if landunit_index
            zldis = fv.forc_hgt_u_patch[lun_patchi[n]] - displa[n]
        else
            zldis = fv.forc_hgt_u_patch[n] - displa[n]
        end
        zeta = zldis / obu[n]

        if zeta < -zetam
            ustar[n] = VKC * um[n] / (log(-zetam * obu[n] / z0m[n]) -
                stability_func1(-zetam) +
                stability_func1(z0m[n] / obu[n]) +
                1.14 * ((-zeta)^0.333 - (zetam)^0.333))
        elseif zeta < 0.0
            ustar[n] = VKC * um[n] / (log(zldis / z0m[n]) -
                stability_func1(zeta) +
                stability_func1(z0m[n] / obu[n]))
        elseif zeta <= 1.0
            ustar[n] = VKC * um[n] / (log(zldis / z0m[n]) + 5.0 * zeta - 5.0 * z0m[n] / obu[n])
        else
            ustar[n] = VKC * um[n] / (log(obu[n] / z0m[n]) + 5.0 - 5.0 * z0m[n] / obu[n] +
                (5.0 * log(zeta) + zeta - 1.0))
        end

        # Deposition velocity
        if zeta < 0.0
            vds_tmp = 2.0e-3 * ustar[n] * (1.0 + (300.0 / (-obu[n]))^0.666)
        else
            vds_tmp = 2.0e-3 * ustar[n]
        end

        if landunit_index
            for pp in lun_patchi[n]:lun_patchf[n]
                fv.vds_patch[pp] = vds_tmp
            end
        else
            fv.vds_patch[n] = vds_tmp
        end

        # ------- 10-m wind (CLM) -------
        if landunit_index
            for pp in lun_patchi[n]:lun_patchf[n]
                if zldis - z0m[n] <= 10.0
                    fv.u10_clm_patch[pp] = um[n]
                else
                    if zeta < -zetam
                        fv.u10_clm_patch[pp] = um[n] - (ustar[n] / VKC * (log(-zetam * obu[n] / (10.0 + z0m[n])) -
                            stability_func1(-zetam) +
                            stability_func1((10.0 + z0m[n]) / obu[n]) +
                            1.14 * ((-zeta)^0.333 - (zetam)^0.333)))
                    elseif zeta < 0.0
                        fv.u10_clm_patch[pp] = um[n] - (ustar[n] / VKC * (log(zldis / (10.0 + z0m[n])) -
                            stability_func1(zeta) +
                            stability_func1((10.0 + z0m[n]) / obu[n])))
                    elseif zeta <= 1.0
                        fv.u10_clm_patch[pp] = um[n] - (ustar[n] / VKC * (log(zldis / (10.0 + z0m[n])) +
                            5.0 * zeta - 5.0 * (10.0 + z0m[n]) / obu[n]))
                    else
                        fv.u10_clm_patch[pp] = um[n] - (ustar[n] / VKC * (log(obu[n] / (10.0 + z0m[n])) +
                            5.0 - 5.0 * (10.0 + z0m[n]) / obu[n] +
                            (5.0 * log(zeta) + zeta - 1.0)))
                    end
                end
                fv.va_patch[pp] = um[n]
            end
        else
            if zldis - z0m[n] <= 10.0
                fv.u10_clm_patch[n] = um[n]
            else
                if zeta < -zetam
                    fv.u10_clm_patch[n] = um[n] - (ustar[n] / VKC * (log(-zetam * obu[n] / (10.0 + z0m[n])) -
                        stability_func1(-zetam) +
                        stability_func1((10.0 + z0m[n]) / obu[n]) +
                        1.14 * ((-zeta)^0.333 - (zetam)^0.333)))
                elseif zeta < 0.0
                    fv.u10_clm_patch[n] = um[n] - (ustar[n] / VKC * (log(zldis / (10.0 + z0m[n])) -
                        stability_func1(zeta) +
                        stability_func1((10.0 + z0m[n]) / obu[n])))
                elseif zeta <= 1.0
                    fv.u10_clm_patch[n] = um[n] - (ustar[n] / VKC * (log(zldis / (10.0 + z0m[n])) +
                        5.0 * zeta - 5.0 * (10.0 + z0m[n]) / obu[n]))
                else
                    fv.u10_clm_patch[n] = um[n] - (ustar[n] / VKC * (log(obu[n] / (10.0 + z0m[n])) +
                        5.0 - 5.0 * (10.0 + z0m[n]) / obu[n] +
                        (5.0 * log(zeta) + zeta - 1.0)))
                end
            end
            fv.va_patch[n] = um[n]
        end

        # ------- Temperature profile -------
        if landunit_index
            zldis = fv.forc_hgt_t_patch[lun_patchi[n]] - displa[n]
        else
            zldis = fv.forc_hgt_t_patch[n] - displa[n]
        end
        zeta = zldis / obu[n]

        if zeta < -zetat
            temp1[n] = VKC / (log(-zetat * obu[n] / z0h[n]) -
                stability_func2(-zetat) +
                stability_func2(z0h[n] / obu[n]) +
                0.8 * ((zetat)^(-0.333) - (-zeta)^(-0.333)))
        elseif zeta < 0.0
            temp1[n] = VKC / (log(zldis / z0h[n]) -
                stability_func2(zeta) +
                stability_func2(z0h[n] / obu[n]))
        elseif zeta <= 1.0
            temp1[n] = VKC / (log(zldis / z0h[n]) + 5.0 * zeta - 5.0 * z0h[n] / obu[n])
        else
            temp1[n] = VKC / (log(obu[n] / z0h[n]) + 5.0 - 5.0 * z0h[n] / obu[n] +
                (5.0 * log(zeta) + zeta - 1.0))
        end

        # ------- Humidity profile -------
        if landunit_index
            hgt_q = fv.forc_hgt_q_patch[lun_patchi[n]]
            hgt_t = fv.forc_hgt_t_patch[lun_patchi[n]]
        else
            hgt_q = fv.forc_hgt_q_patch[n]
            hgt_t = fv.forc_hgt_t_patch[n]
        end

        if hgt_q == hgt_t && z0q[n] == z0h[n]
            temp2[n] = temp1[n]
        else
            zldis = hgt_q - displa[n]
            zeta = zldis / obu[n]
            if zeta < -zetat
                temp2[n] = VKC / (log(-zetat * obu[n] / z0q[n]) -
                    stability_func2(-zetat) +
                    stability_func2(z0q[n] / obu[n]) +
                    0.8 * ((zetat)^(-0.333) - (-zeta)^(-0.333)))
            elseif zeta < 0.0
                temp2[n] = VKC / (log(zldis / z0q[n]) -
                    stability_func2(zeta) +
                    stability_func2(z0q[n] / obu[n]))
            elseif zeta <= 1.0
                temp2[n] = VKC / (log(zldis / z0q[n]) + 5.0 * zeta - 5.0 * z0q[n] / obu[n])
            else
                temp2[n] = VKC / (log(obu[n] / z0q[n]) + 5.0 - 5.0 * z0q[n] / obu[n] +
                    (5.0 * log(zeta) + zeta - 1.0))
            end
        end

        # ------- Temperature profile at 2m -------
        zldis = 2.0 + z0h[n]
        zeta = zldis / obu[n]
        if zeta < -zetat
            temp12m[n] = VKC / (log(-zetat * obu[n] / z0h[n]) -
                stability_func2(-zetat) +
                stability_func2(z0h[n] / obu[n]) +
                0.8 * ((zetat)^(-0.333) - (-zeta)^(-0.333)))
        elseif zeta < 0.0
            temp12m[n] = VKC / (log(zldis / z0h[n]) -
                stability_func2(zeta) +
                stability_func2(z0h[n] / obu[n]))
        elseif zeta <= 1.0
            temp12m[n] = VKC / (log(zldis / z0h[n]) + 5.0 * zeta - 5.0 * z0h[n] / obu[n])
        else
            temp12m[n] = VKC / (log(obu[n] / z0h[n]) + 5.0 - 5.0 * z0h[n] / obu[n] +
                (5.0 * log(zeta) + zeta - 1.0))
        end

        # ------- Humidity profile at 2m -------
        if z0q[n] == z0h[n]
            temp22m[n] = temp12m[n]
        else
            zldis = 2.0 + z0q[n]
            zeta = zldis / obu[n]
            if zeta < -zetat
                temp22m[n] = VKC / (log(-zetat * obu[n] / z0q[n]) -
                    stability_func2(-zetat) +
                    stability_func2(z0q[n] / obu[n]) +
                    0.8 * ((zetat)^(-0.333) - (-zeta)^(-0.333)))
            elseif zeta < 0.0
                temp22m[n] = VKC / (log(zldis / z0q[n]) -
                    stability_func2(zeta) +
                    stability_func2(z0q[n] / obu[n]))
            elseif zeta <= 1.0
                temp22m[n] = VKC / (log(zldis / z0q[n]) + 5.0 * zeta - 5.0 * z0q[n] / obu[n])
            else
                temp22m[n] = VKC / (log(obu[n] / z0q[n]) + 5.0 - 5.0 * z0q[n] / obu[n] +
                    (5.0 * log(zeta) + zeta - 1.0))
            end
        end

        # ------- 10-m wind for dust model -------
        if landunit_index
            zldis = fv.forc_hgt_u_patch[lun_patchi[n]] - displa[n]
        else
            zldis = fv.forc_hgt_u_patch[n] - displa[n]
        end
        zeta = zldis / obu[n]

        if min(zeta, 1.0) < 0.0
            tmp1 = (1.0 - 16.0 * min(zeta, 1.0))^0.25
            tmp2 = log((1.0 + tmp1 * tmp1) / 2.0)
            tmp3 = log((1.0 + tmp1) / 2.0)
            fmnew = 2.0 * tmp3 + tmp2 - 2.0 * atan(tmp1) + 1.5707963
        else
            fmnew = -5.0 * min(zeta, 1.0)
        end
        if iter == 1
            fm[n] = fmnew
        else
            fm[n] = 0.5 * (fm[n] + fmnew)
        end

        zeta10 = min(10.0 / obu[n], 1.0)
        if zeta == 0.0
            zeta10 = 0.0
        end
        if zeta10 < 0.0
            tmp1 = (1.0 - 16.0 * zeta10)^0.25
            tmp2 = log((1.0 + tmp1 * tmp1) / 2.0)
            tmp3 = log((1.0 + tmp1) / 2.0)
            fm10 = 2.0 * tmp3 + tmp2 - 2.0 * atan(tmp1) + 1.5707963
        else
            fm10 = -5.0 * zeta10
        end

        if landunit_index
            tmp4 = log(max(1.0, fv.forc_hgt_u_patch[lun_patchi[n]] / 10.0))
        else
            tmp4 = log(max(1.0, fv.forc_hgt_u_patch[n] / 10.0))
        end

        if landunit_index
            for pp in lun_patchi[n]:lun_patchf[n]
                fv.u10_patch[pp] = ur[n] - ustar[n] / VKC * (tmp4 - fm[n] + fm10)
                fv.fv_patch[pp] = ustar[n]
            end
        else
            fv.u10_patch[n] = ur[n] - ustar[n] / VKC * (tmp4 - fm[n] + fm10)
            fv.fv_patch[n] = ustar[n]
        end

    end  # end filter loop

    return nothing
end
