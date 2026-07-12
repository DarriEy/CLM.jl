# ==========================================================================
# Ported from: src/biogeophys/FrictionVelocityMod.F90
# Friction velocity data type, stability functions, Monin-Obukhov
# initialization, roughness length setting, and friction velocity calculation.
# ==========================================================================

# Snow-cover fraction below which a column is treated as snow-free when picking
# the ground momentum roughness (zsno vs zlnd). Fortran's frac_sno is exactly 0
# when snow is absent; Julia's snow-water bookkeeping can leave a subnormal
# (~1e-18) that must not read as snow. Far below physical snow-cover onset (~1e-2).
const FRAC_SNO_TRACE_FLOOR = 1.0e-11

"""
    FrictionVelocityData

Friction velocity data structure. Holds scalar parameters (roughness lengths,
stability limit) and per-patch / per-column arrays for friction velocity
computations, roughness lengths, forcing heights, and diagnostic variables.

Ported from `frictionvel_type` in `FrictionVelocityMod.F90`.
"""
Base.@kwdef mutable struct FrictionVelocityData{FT<:Real,
                                    V<:AbstractVector{FT}}
    # --- Scalar parameters ---
    zetamaxstable::FT = -999.0   # Max zeta under stable conditions
    zsno::FT          = -999.0   # Momentum roughness length for snow (m)
    zlnd::FT          = -999.0   # Momentum roughness length for soil, glacier, wetland (m)
    zglc::FT          = -999.0   # Momentum roughness length for glacier (Meier2022 only) (m)

    # --- Patch-level fields (1D vectors) ---
    forc_hgt_u_patch ::V = Float64[]   # patch wind forcing height (10m+z0m+d) (m)
    forc_hgt_t_patch ::V = Float64[]   # patch temperature forcing height (10m+z0m+d) (m)
    forc_hgt_q_patch ::V = Float64[]   # patch specific humidity forcing height (10m+z0m+d) (m)
    u10_patch        ::V = Float64[]   # patch 10-m wind (m/s) (for dust model)
    u10_clm_patch    ::V = Float64[]   # patch 10-m wind (m/s) (for clm_map2gcell)
    va_patch         ::V = Float64[]   # patch atmospheric wind speed plus convective velocity (m/s)
    vds_patch        ::V = Float64[]   # patch deposition velocity term (m/s)
    fv_patch         ::V = Float64[]   # patch friction velocity (m/s) (for dust model)
    rb1_patch        ::V = Float64[]   # patch aerodynamical resistance (s/m)
    rb10_patch       ::V = Float64[]   # 10-day mean patch aerodynamical resistance (s/m)
    ram1_patch       ::V = Float64[]   # patch aerodynamical resistance (s/m)
    z0mv_patch       ::V = Float64[]   # patch roughness length over vegetation, momentum [m]
    z0hv_patch       ::V = Float64[]   # patch roughness length over vegetation, sensible heat [m]
    z0qv_patch       ::V = Float64[]   # patch roughness length over vegetation, latent heat [m]
    z0mg_patch       ::V = Float64[]   # patch roughness length over ground, momentum [m]
    z0hg_patch       ::V = Float64[]   # patch roughness length over ground, sensible heat [m]
    z0qg_patch       ::V = Float64[]   # patch roughness length over ground, latent heat [m]
    kbm1_patch       ::V = Float64[]   # natural logarithm of z0mg_p/z0hg_p [-]
    rah1_patch       ::V = Float64[]   # patch sensible heat flux resistance [s/m]
    rah2_patch       ::V = Float64[]   # patch below-canopy sensible heat flux resistance [s/m]
    raw1_patch       ::V = Float64[]   # patch moisture flux resistance [s/m]
    raw2_patch       ::V = Float64[]   # patch below-canopy moisture flux resistance [s/m]
    ustar_patch      ::V = Float64[]   # patch friction velocity [m/s]
    um_patch         ::V = Float64[]   # patch wind speed including the stability effect [m/s]
    uaf_patch        ::V = Float64[]   # patch canopy air speed [m/s]
    taf_patch        ::V = Float64[]   # patch canopy air temperature [K]
    qaf_patch        ::V = Float64[]   # patch canopy humidity [kg/kg]
    obu_patch        ::V = Float64[]   # patch Monin-Obukhov length [m]
    zeta_patch       ::V = Float64[]   # patch dimensionless stability parameter
    vpd_patch        ::V = Float64[]   # patch vapor pressure deficit [Pa]
    num_iter_patch   ::V = Float64[]   # patch number of iterations
    z0m_actual_patch ::V = Float64[]   # patch roughness length actually used, momentum [m]

    # --- Column-level fields (1D vectors) ---
    z0mg_col    ::V = Float64[]   # col roughness length over ground, momentum [m]
    z0mg_2D_col ::V = Float64[]   # 2-D field of input col roughness length over ground, momentum [m]
    z0hg_col    ::V = Float64[]   # col roughness length over ground, sensible heat [m]
    z0qg_col    ::V = Float64[]   # col roughness length over ground, latent heat [m]
end

FrictionVelocityData{FT}(; kwargs...) where {FT<:Real} =
    FrictionVelocityData{FT, Vector{FT}}(; kwargs...)
Adapt.@adapt_structure FrictionVelocityData


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
function frictionvel_init!(fv::FrictionVelocityData{FT}, np::Int, nc::Int) where {FT}
    # Patch-level fields
    fv.forc_hgt_u_patch  = fill(FT(NaN), np)
    fv.forc_hgt_t_patch  = fill(FT(NaN), np)
    fv.forc_hgt_q_patch  = fill(FT(NaN), np)
    fv.u10_patch         = fill(FT(NaN), np)
    fv.u10_clm_patch     = fill(FT(NaN), np)
    fv.va_patch          = fill(FT(NaN), np)
    fv.vds_patch         = fill(FT(NaN), np)
    fv.fv_patch          = fill(FT(NaN), np)
    fv.rb1_patch         = fill(FT(NaN), np)
    fv.rb10_patch        = fill(FT(SPVAL), np)
    fv.ram1_patch        = fill(FT(NaN), np)
    fv.z0mv_patch        = fill(FT(NaN), np)
    fv.z0hv_patch        = fill(FT(NaN), np)
    fv.z0qv_patch        = fill(FT(NaN), np)
    fv.z0mg_patch        = fill(FT(NaN), np)
    fv.z0hg_patch        = fill(FT(NaN), np)
    fv.z0qg_patch        = fill(FT(NaN), np)
    fv.kbm1_patch        = fill(FT(NaN), np)
    fv.rah1_patch        = fill(FT(NaN), np)
    fv.rah2_patch        = fill(FT(NaN), np)
    fv.raw1_patch        = fill(FT(NaN), np)
    fv.raw2_patch        = fill(FT(NaN), np)
    fv.ustar_patch       = fill(FT(NaN), np)
    fv.um_patch          = fill(FT(NaN), np)
    fv.uaf_patch         = fill(FT(NaN), np)
    fv.taf_patch         = fill(FT(NaN), np)
    fv.qaf_patch         = fill(FT(NaN), np)
    fv.obu_patch         = fill(FT(NaN), np)
    fv.zeta_patch        = fill(FT(NaN), np)
    fv.vpd_patch         = fill(FT(NaN), np)
    fv.num_iter_patch    = fill(FT(NaN), np)
    fv.z0m_actual_patch  = fill(FT(NaN), np)

    # Column-level fields
    fv.z0mg_col    = fill(FT(NaN), nc)
    fv.z0mg_2D_col = fill(FT(NaN), nc)
    fv.z0hg_col    = fill(FT(NaN), nc)
    fv.z0qg_col    = fill(FT(NaN), nc)

    return nothing
end

"""
    frictionvel_clean!(fv::FrictionVelocityData)

Deallocate (reset to empty) all fields of a `FrictionVelocityData` instance.
"""
function frictionvel_clean!(fv::FrictionVelocityData{FT}) where {FT}
    fv.forc_hgt_u_patch  = FT[]
    fv.forc_hgt_t_patch  = FT[]
    fv.forc_hgt_q_patch  = FT[]
    fv.u10_patch         = FT[]
    fv.u10_clm_patch     = FT[]
    fv.va_patch          = FT[]
    fv.vds_patch         = FT[]
    fv.fv_patch          = FT[]
    fv.rb1_patch         = FT[]
    fv.rb10_patch        = FT[]
    fv.ram1_patch        = FT[]
    fv.z0mv_patch        = FT[]
    fv.z0hv_patch        = FT[]
    fv.z0qv_patch        = FT[]
    fv.z0mg_patch        = FT[]
    fv.z0hg_patch        = FT[]
    fv.z0qg_patch        = FT[]
    fv.kbm1_patch        = FT[]
    fv.rah1_patch        = FT[]
    fv.rah2_patch        = FT[]
    fv.raw1_patch        = FT[]
    fv.raw2_patch        = FT[]
    fv.ustar_patch       = FT[]
    fv.um_patch          = FT[]
    fv.uaf_patch         = FT[]
    fv.taf_patch         = FT[]
    fv.qaf_patch         = FT[]
    fv.obu_patch         = FT[]
    fv.zeta_patch        = FT[]
    fv.vpd_patch         = FT[]
    fv.num_iter_patch    = FT[]
    fv.z0m_actual_patch  = FT[]
    fv.z0mg_col          = FT[]
    fv.z0mg_2D_col       = FT[]
    fv.z0hg_col          = FT[]
    fv.z0qg_col          = FT[]
    return nothing
end

"""
    frictionvel_init_for_testing!(fv::FrictionVelocityData, np::Int, nc::Int)

Initialize for unit testing with hardcoded namelist and parameter values.

Ported from `frictionvel_type%InitForTesting` in `FrictionVelocityMod.F90`.
"""
function frictionvel_init_for_testing!(fv::FrictionVelocityData{FT}, np::Int, nc::Int) where {FT}
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
                           col_landunit::Union{AbstractVector{<:Integer},Nothing}=nothing,
                           lun_lakpoi::Union{Vector{Bool},Nothing}=nothing)

Cold-start initialization. Sets forc_hgt_u_patch to 30m for CN mode,
and z0mg_col to 0.0004 for lake columns.

Ported from `frictionvel_type%InitCold` in `FrictionVelocityMod.F90`.
"""
function frictionvel_init_cold!(fv::FrictionVelocityData,
                                 bounds_col::UnitRange{Int},
                                 bounds_patch::UnitRange{Int};
                                 col_landunit::Union{AbstractVector{<:Integer},Nothing} = nothing,
                                 lun_lakpoi::Union{Vector{Bool},Nothing} = nothing,
                                 use_cn::Bool = false)
    for p in bounds_patch
        if isnan(fv.forc_hgt_u_patch[p])
            fv.forc_hgt_u_patch[p] = 30.0
        end
        if isnan(fv.forc_hgt_t_patch[p])
            fv.forc_hgt_t_patch[p] = 30.0
        end
        if isnan(fv.forc_hgt_q_patch[p])
            fv.forc_hgt_q_patch[p] = 30.0
        end
    end

    if col_landunit !== nothing && lun_lakpoi !== nothing
        for c in bounds_col
            l = col_landunit[c]
            if lun_lakpoi[l]
                fv.z0mg_col[c] = 0.0004
            else
                fv.z0mg_col[c] = fv.zlnd
            end
            fv.z0hg_col[c] = fv.z0mg_col[c]
            fv.z0qg_col[c] = fv.z0mg_col[c]
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
Not implemented (no-op stub). History I/O IS ported
(`src/infrastructure/history_io.jl`); fields are registered centrally.

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
Not implemented (no-op stub). Restart I/O IS ported
(`src/infrastructure/restart_io.jl`); variables are declared centrally.
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
                                zetamaxstable::Real = 0.5)
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
                                   zsno::Real = 0.00085,
                                   zlnd::Real = 0.000775,
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
function stability_func1(zeta::Real)
    # Guard: this function is only valid for zeta <= 0 (unstable conditions).
    # Type-preserving (eltype-generic) so it lowers to valid Metal IR under
    # Float32; byte-identical to the Float64 literals on the CPU path.
    T = typeof(zeta)
    z = min(zeta, zero(T))
    chik2 = sqrt(one(T) - T(16.0) * z)
    chik = sqrt(chik2)
    return T(2.0) * log((one(T) + chik) * T(0.5)) +
           log((one(T) + chik2) * T(0.5)) -
           T(2.0) * atan(chik) + T(π) * T(0.5)
end

"""
    stability_func2(zeta::Float64) -> Float64

Stability function for rib < 0 (temperature/humidity profile).
Ported from `StabilityFunc2` in `FrictionVelocityMod.F90`.
"""
function stability_func2(zeta::Real)
    # Guard: this function is only valid for zeta <= 0 (unstable conditions).
    # Type-preserving (see stability_func1).
    T = typeof(zeta)
    z = min(zeta, zero(T))
    chik2 = sqrt(one(T) - T(16.0) * z)
    return T(2.0) * log((one(T) + chik2) * T(0.5))
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
function monin_obuk_ini(zetamaxstable::Real, ur::Real, thv::Real,
                        dthv::Real, zldis::Real, z0m::Real)
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
# --------------------------------------------------------------------------
# Kernel: ground (column-level) roughness lengths for non-lake columns.
# Each column writes only z0mg_col[c]/z0hg_col[c]/z0qg_col[c]; independent.
# `z0method` resolves z0param_method on the host (1=ZengWang2007, 2=Meier2022)
# and `use_snowmelt` flags use_z0m_snowmelt; zsno/zlnd/zglc and the Meier
# Float64 globals (B1/B4/RPI/MEIER1/2) arrive eltype-converted so no Float64
# reaches a Float32-only backend. `lo` = first(bounds_col).
# --------------------------------------------------------------------------
@kernel function _set_ground_roughness_kernel!(z0mg_col, z0hg_col, z0qg_col,
                                               @Const(mask), @Const(col_landunit),
                                               @Const(frac_sno), @Const(snomelt_accum),
                                               @Const(lun_itype),
                                               z0method::Int, use_snowmelt::Bool,
                                               zsno, zlnd, zglc,
                                               b1_param, b4_param, rpi,
                                               meier1, meier2, lo::Int)
    i = @index(Global)
    @inbounds begin
        c = lo + i - 1
        if mask[c]
            T = eltype(z0mg_col)
            # Fortran tests frac_sno > 0 on a value that is exactly 0 when snow is
            # absent. Julia's snow-water bookkeeping can leave a physically-
            # negligible subnormal frac_sno (~1e-18) that would flip this discrete
            # branch to the snow roughness (zsno ≈ 0.0024) instead of the soil
            # roughness (zlnd ≈ 0.01), ~halving the below-canopy resistance rah2 at
            # low-LAI/dry canopies and warm-biasing the coupled canopy-air solve.
            # A tiny threshold reproduces Fortran's exact-zero behaviour on every
            # real case (physical snow cover onset is O(1e-2)).
            has_snow = frac_sno[c] > T(FRAC_SNO_TRACE_FLOOR)
            if z0method == 1  # ZengWang2007
                if has_snow
                    z0mg_col[c] = zsno
                else
                    z0mg_col[c] = zlnd
                end
            elseif z0method == 2  # Meier2022
                l = col_landunit[c]
                if has_snow
                    if use_snowmelt
                        if snomelt_accum[c] < T(1.0e-5)
                            z0mg_col[c] = exp(-b1_param * rpi * T(0.5) + b4_param) * T(1.0e-3)
                        else
                            z0mg_col[c] = exp(b1_param * (atan((log10(snomelt_accum[c]) + meier1) / meier2)) + b4_param) * T(1.0e-3)
                        end
                    else
                        z0mg_col[c] = zsno
                    end
                elseif lun_itype[l] == ISTICE
                    z0mg_col[c] = zglc
                else
                    z0mg_col[c] = zlnd
                end
            end

            z0hg_col[c] = z0mg_col[c]  # initial set only
            z0qg_col[c] = z0mg_col[c]  # initial set only
        end
    end
end

# --------------------------------------------------------------------------
# Kernel: vegetation roughness lengths + forcing heights for non-lake patches.
# Each patch writes only its own indices; reads column ground roughness
# (must be set by _set_ground_roughness_kernel! first). SPVAL stores are bare
# (constant-folded). `lo` = first(bounds_patch).
# --------------------------------------------------------------------------
@kernel function _set_veg_roughness_forcheights_kernel!(
        z0mv_patch, z0hv_patch, z0qv_patch, z0mg_patch, z0hg_patch, z0qg_patch,
        kbm1_patch, forc_hgt_u_patch, forc_hgt_t_patch, forc_hgt_q_patch,
        @Const(mask), @Const(z0m), @Const(displa),
        @Const(z0mg_col), @Const(z0hg_col), @Const(z0qg_col),
        @Const(frac_veg_nosno), @Const(forc_hgt_u), @Const(forc_hgt_t), @Const(forc_hgt_q),
        @Const(patch_gridcell), @Const(patch_landunit), @Const(patch_column),
        @Const(lun_itype), @Const(lun_urbpoi), @Const(lun_z_0_town), @Const(lun_z_d_town),
        lo::Int)
    i = @index(Global)
    @inbounds begin
        p = lo + i - 1
        if mask[p]
            # Roughness lengths over vegetation
            z0mv_patch[p] = z0m[p]
            z0hv_patch[p] = z0mv_patch[p]
            z0qv_patch[p] = z0mv_patch[p]

            # Set to arbitrary value (will be overwritten by respective modules)
            z0mg_patch[p] = SPVAL
            z0hg_patch[p] = SPVAL
            z0qg_patch[p] = SPVAL
            kbm1_patch[p] = SPVAL

            # Forcing height = atmospheric forcing height + z0m + displa
            g = patch_gridcell[p]
            l = patch_landunit[p]
            c = patch_column[p]

            if lun_itype[l] == ISTSOIL || lun_itype[l] == ISTCROP
                if frac_veg_nosno[p] == 0
                    forc_hgt_u_patch[p] = forc_hgt_u[g] + z0mg_col[c] + displa[p]
                    forc_hgt_t_patch[p] = forc_hgt_t[g] + z0hg_col[c] + displa[p]
                    forc_hgt_q_patch[p] = forc_hgt_q[g] + z0qg_col[c] + displa[p]
                else
                    forc_hgt_u_patch[p] = forc_hgt_u[g] + z0mv_patch[p] + displa[p]
                    forc_hgt_t_patch[p] = forc_hgt_t[g] + z0hv_patch[p] + displa[p]
                    forc_hgt_q_patch[p] = forc_hgt_q[g] + z0qv_patch[p] + displa[p]
                end
            elseif lun_itype[l] == ISTWET || lun_itype[l] == ISTICE
                forc_hgt_u_patch[p] = forc_hgt_u[g] + z0mg_col[c] + displa[p]
                forc_hgt_t_patch[p] = forc_hgt_t[g] + z0hg_col[c] + displa[p]
                forc_hgt_q_patch[p] = forc_hgt_q[g] + z0qg_col[c] + displa[p]
            elseif lun_urbpoi[l]
                forc_hgt_u_patch[p] = forc_hgt_u[g] + lun_z_0_town[l] + lun_z_d_town[l]
                forc_hgt_t_patch[p] = forc_hgt_t[g] + lun_z_0_town[l] + lun_z_d_town[l]
                forc_hgt_q_patch[p] = forc_hgt_q[g] + lun_z_0_town[l] + lun_z_d_town[l]
            end
        end
    end
end

function set_roughness_and_forc_heights_nonlake!(
        fv::FrictionVelocityData,
        mask_nolakec::AbstractVector{Bool},
        mask_nolakep::AbstractVector{Bool},
        bounds_col::UnitRange{Int},
        bounds_patch::UnitRange{Int},
        frac_sno::AbstractVector{<:Real},
        snomelt_accum::AbstractVector{<:Real},
        frac_veg_nosno::AbstractVector{<:Integer},
        z0m::AbstractVector{<:Real},
        displa::AbstractVector{<:Real},
        forc_hgt_u::AbstractVector{<:Real},
        forc_hgt_t::AbstractVector{<:Real},
        forc_hgt_q::AbstractVector{<:Real},
        col_landunit::AbstractVector{<:Integer},
        patch_gridcell::AbstractVector{<:Integer},
        patch_landunit::AbstractVector{<:Integer},
        patch_column::AbstractVector{<:Integer},
        lun_itype::AbstractVector{<:Integer},
        lun_urbpoi::AbstractVector{Bool},
        lun_z_0_town::AbstractVector{<:Real},
        lun_z_d_town::AbstractVector{<:Real};
        z0param_method::String = "",
        use_z0m_snowmelt::Bool = false)

    # Resolve the String config to an Int flag on the host (no String compares
    # inside the kernel): 1 = ZengWang2007, 2 = Meier2022, 0 = neither.
    z0method = z0param_method == "ZengWang2007" ? 1 :
               (z0param_method == "Meier2022" ? 2 : 0)

    T = eltype(fv.z0mg_col)

    # Column loop → kernel: set ground roughness lengths
    if !isempty(bounds_col)
        _launch!(_set_ground_roughness_kernel!, fv.z0mg_col, fv.z0hg_col, fv.z0qg_col,
                 mask_nolakec, col_landunit, frac_sno, snomelt_accum, lun_itype,
                 z0method, use_z0m_snowmelt,
                 T(fv.zsno), T(fv.zlnd), T(fv.zglc),
                 T(B1_PARAM), T(B4_PARAM), T(RPI), T(MEIER_PARAM1), T(MEIER_PARAM2),
                 first(bounds_col); ndrange = length(bounds_col))
    end

    # Patch loop → kernel: vegetation roughness + forcing heights (reads the
    # ground roughness set above, so it launches after the column kernel).
    if !isempty(bounds_patch)
        _launch!(_set_veg_roughness_forcheights_kernel!,
                 fv.z0mv_patch, fv.z0hv_patch, fv.z0qv_patch,
                 fv.z0mg_patch, fv.z0hg_patch, fv.z0qg_patch, fv.kbm1_patch,
                 fv.forc_hgt_u_patch, fv.forc_hgt_t_patch, fv.forc_hgt_q_patch,
                 mask_nolakep, z0m, displa,
                 fv.z0mg_col, fv.z0hg_col, fv.z0qg_col,
                 frac_veg_nosno, forc_hgt_u, forc_hgt_t, forc_hgt_q,
                 patch_gridcell, patch_landunit, patch_column,
                 lun_itype, lun_urbpoi, lun_z_0_town, lun_z_d_town,
                 first(bounds_patch); ndrange = length(bounds_patch))
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
# --------------------------------------------------------------------------
# Kernel: select the actual momentum roughness length per patch by priority
# (exposed-veg > no-exposed-veg > urban > lake). Each patch writes only its own
# index; fully independent. `lo` = first(bounds_patch).
# --------------------------------------------------------------------------
@kernel function _set_actual_z0m_kernel!(z0m_actual_patch, @Const(z0mv_patch),
                                         @Const(z0mg_col), @Const(lun_z_0_town),
                                         @Const(mask_exposedvegp), @Const(mask_noexposedvegp),
                                         @Const(mask_urbanp), @Const(mask_lakep),
                                         @Const(patch_column), @Const(patch_landunit),
                                         lo::Int)
    i = @index(Global)
    @inbounds begin
        p = lo + i - 1
        if mask_exposedvegp[p]
            z0m_actual_patch[p] = z0mv_patch[p]
        elseif mask_noexposedvegp[p]
            z0m_actual_patch[p] = z0mg_col[patch_column[p]]
        elseif mask_urbanp[p]
            z0m_actual_patch[p] = lun_z_0_town[patch_landunit[p]]
        elseif mask_lakep[p]
            z0m_actual_patch[p] = z0mg_col[patch_column[p]]
        end
    end
end

function set_actual_roughness_lengths!(
        fv::FrictionVelocityData,
        mask_exposedvegp::AbstractVector{Bool},
        mask_noexposedvegp::AbstractVector{Bool},
        mask_urbanp::AbstractVector{Bool},
        mask_lakep::AbstractVector{Bool},
        bounds_patch::UnitRange{Int},
        patch_column::AbstractVector{<:Integer},
        patch_landunit::AbstractVector{<:Integer},
        lun_z_0_town::AbstractVector{<:Real})

    isempty(bounds_patch) && return nothing
    _launch!(_set_actual_z0m_kernel!, fv.z0m_actual_patch, fv.z0mv_patch,
             fv.z0mg_col, lun_z_0_town, mask_exposedvegp, mask_noexposedvegp,
             mask_urbanp, mask_lakep, patch_column, patch_landunit,
             first(bounds_patch); ndrange = length(bounds_patch))

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
        filtern::AbstractVector{<:Integer},
        displa::AbstractVector{<:Real},
        z0m::AbstractVector{<:Real},
        z0h::AbstractVector{<:Real},
        z0q::AbstractVector{<:Real},
        obu::AbstractVector{<:Real},
        iter::Int,
        ur::AbstractVector{<:Real},
        um::AbstractVector{<:Real},
        ustar::AbstractVector{<:Real},
        temp1::AbstractVector{<:Real},
        temp2::AbstractVector{<:Real},
        temp12m::AbstractVector{<:Real},
        temp22m::AbstractVector{<:Real},
        fm::AbstractVector{<:Real};
        landunit_index::Bool = false,
        active = nothing,
        lun_gridcell::Union{AbstractVector{<:Integer},Nothing} = nothing,
        lun_patchi::Union{AbstractVector{<:Integer},Nothing} = nothing,
        lun_patchf::Union{AbstractVector{<:Integer},Nothing} = nothing,
        patch_gridcell::Union{AbstractVector{<:Integer},Nothing} = nothing)

    zetam = 1.574  # transition point of flux-gradient relation (wind profile)
    zetat = 0.465  # transition point of flux-gradient relation (temp. profile)

    if landunit_index
    for f in 1:fn
        n = filtern[f]

        # ------- Wind profile -------
        if landunit_index
            zldis = fv.forc_hgt_u_patch[lun_patchi[n]] - displa[n]
        else
            zldis = fv.forc_hgt_u_patch[n] - displa[n]
        end
        # Defense-in-depth: keep the wind reference height above the roughness length so
        # log(zldis/z0m) can't collapse to 0 (→ ustar = VKC·um/0 blow-up). See the
        # non-landunit branch for the full rationale; valid inputs are unaffected.
        zldis = max(zldis, z0m[n] + 1.0)
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
            vds_tmp = 2.0e-3 * ustar[n] * (1.0 + (300.0 / max(-obu[n], 1.0e-10))^0.666)
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
    else
        fv_profile_update!(fv, fn, filtern, displa, z0m, z0h, z0q, obu, iter, ur, um,
                           ustar, temp1, temp2, temp12m, temp22m, fm; active = active)
    end

    return nothing
end

# ---------------------------------------------------------------------------
# Kernelized default (non-landunit) Monin-Obukhov profile update
# ---------------------------------------------------------------------------

@kernel function _fv_profile_kernel!(
        # written arrays (ustar is written then read; fm is read+written when iter>1)
        ustar, temp1, temp2, temp12m, temp22m, fm,
        vds_patch, u10_clm_patch, va_patch, u10_patch, fv_patch,
        # read-only arrays
        @Const(filtern), @Const(active), @Const(displa), @Const(z0m), @Const(z0h), @Const(z0q),
        @Const(obu), @Const(ur), @Const(um),
        @Const(forc_hgt_u_patch), @Const(forc_hgt_t_patch), @Const(forc_hgt_q_patch),
        # scalars
        iter::Int, VKC, zetam, zetat)

    f = @index(Global)
    @inbounds if active[filtern[f]]
        n = filtern[f]
        T = eltype(ustar)

        # ------- Wind profile -------
        zldis = forc_hgt_u_patch[n] - displa[n]
        # Defense-in-depth: the wind reference height must sit above the roughness
        # length. A degenerate forc_hgt_u (e.g. an obs height of 0 → zldis collapses to
        # z0m) makes log(zldis/z0m)→0 and ustar = VKC·um/0 → blows up to ~1e15. Clamp
        # zldis to ≥ z0m + 1 m; valid inputs (zldis ≫ z0m) are unaffected.
        zldis = max(zldis, z0m[n] + one(T))
        zeta = zldis / obu[n]

        if zeta < -zetam
            ustar[n] = VKC * um[n] / (log(-zetam * obu[n] / z0m[n]) -
                stability_func1(-zetam) +
                stability_func1(z0m[n] / obu[n]) +
                T(1.14) * ((-zeta)^T(0.333) - (zetam)^T(0.333)))
        elseif zeta < zero(T)
            ustar[n] = VKC * um[n] / (log(zldis / z0m[n]) -
                stability_func1(zeta) +
                stability_func1(z0m[n] / obu[n]))
        elseif zeta <= one(T)
            ustar[n] = VKC * um[n] / (log(zldis / z0m[n]) + T(5.0) * zeta - T(5.0) * z0m[n] / obu[n])
        else
            ustar[n] = VKC * um[n] / (log(obu[n] / z0m[n]) + T(5.0) - T(5.0) * z0m[n] / obu[n] +
                (T(5.0) * log(zeta) + zeta - one(T)))
        end

        # Deposition velocity
        if zeta < zero(T)
            vds_tmp = T(2.0e-3) * ustar[n] * (one(T) + (T(300.0) / max(-obu[n], T(1.0e-10)))^T(0.666))
        else
            vds_tmp = T(2.0e-3) * ustar[n]
        end
        vds_patch[n] = vds_tmp

        # ------- 10-m wind (CLM) -------
        if zldis - z0m[n] <= T(10.0)
            u10_clm_patch[n] = um[n]
        else
            if zeta < -zetam
                u10_clm_patch[n] = um[n] - (ustar[n] / VKC * (log(-zetam * obu[n] / (T(10.0) + z0m[n])) -
                    stability_func1(-zetam) +
                    stability_func1((T(10.0) + z0m[n]) / obu[n]) +
                    T(1.14) * ((-zeta)^T(0.333) - (zetam)^T(0.333))))
            elseif zeta < zero(T)
                u10_clm_patch[n] = um[n] - (ustar[n] / VKC * (log(zldis / (T(10.0) + z0m[n])) -
                    stability_func1(zeta) +
                    stability_func1((T(10.0) + z0m[n]) / obu[n])))
            elseif zeta <= one(T)
                u10_clm_patch[n] = um[n] - (ustar[n] / VKC * (log(zldis / (T(10.0) + z0m[n])) +
                    T(5.0) * zeta - T(5.0) * (T(10.0) + z0m[n]) / obu[n]))
            else
                u10_clm_patch[n] = um[n] - (ustar[n] / VKC * (log(obu[n] / (T(10.0) + z0m[n])) +
                    T(5.0) - T(5.0) * (T(10.0) + z0m[n]) / obu[n] +
                    (T(5.0) * log(zeta) + zeta - one(T))))
            end
        end
        va_patch[n] = um[n]

        # ------- Temperature profile -------
        zldis = forc_hgt_t_patch[n] - displa[n]
        zeta = zldis / obu[n]

        if zeta < -zetat
            temp1[n] = VKC / (log(-zetat * obu[n] / z0h[n]) -
                stability_func2(-zetat) +
                stability_func2(z0h[n] / obu[n]) +
                T(0.8) * ((zetat)^(-T(0.333)) - (-zeta)^(-T(0.333))))
        elseif zeta < zero(T)
            temp1[n] = VKC / (log(zldis / z0h[n]) -
                stability_func2(zeta) +
                stability_func2(z0h[n] / obu[n]))
        elseif zeta <= one(T)
            temp1[n] = VKC / (log(zldis / z0h[n]) + T(5.0) * zeta - T(5.0) * z0h[n] / obu[n])
        else
            temp1[n] = VKC / (log(obu[n] / z0h[n]) + T(5.0) - T(5.0) * z0h[n] / obu[n] +
                (T(5.0) * log(zeta) + zeta - one(T)))
        end

        # ------- Humidity profile -------
        hgt_q = forc_hgt_q_patch[n]
        hgt_t = forc_hgt_t_patch[n]

        if hgt_q == hgt_t && z0q[n] == z0h[n]
            temp2[n] = temp1[n]
        else
            zldis = hgt_q - displa[n]
            zeta = zldis / obu[n]
            if zeta < -zetat
                temp2[n] = VKC / (log(-zetat * obu[n] / z0q[n]) -
                    stability_func2(-zetat) +
                    stability_func2(z0q[n] / obu[n]) +
                    T(0.8) * ((zetat)^(-T(0.333)) - (-zeta)^(-T(0.333))))
            elseif zeta < zero(T)
                temp2[n] = VKC / (log(zldis / z0q[n]) -
                    stability_func2(zeta) +
                    stability_func2(z0q[n] / obu[n]))
            elseif zeta <= one(T)
                temp2[n] = VKC / (log(zldis / z0q[n]) + T(5.0) * zeta - T(5.0) * z0q[n] / obu[n])
            else
                temp2[n] = VKC / (log(obu[n] / z0q[n]) + T(5.0) - T(5.0) * z0q[n] / obu[n] +
                    (T(5.0) * log(zeta) + zeta - one(T)))
            end
        end

        # ------- Temperature profile at 2m -------
        zldis = T(2.0) + z0h[n]
        zeta = zldis / obu[n]
        if zeta < -zetat
            temp12m[n] = VKC / (log(-zetat * obu[n] / z0h[n]) -
                stability_func2(-zetat) +
                stability_func2(z0h[n] / obu[n]) +
                T(0.8) * ((zetat)^(-T(0.333)) - (-zeta)^(-T(0.333))))
        elseif zeta < zero(T)
            temp12m[n] = VKC / (log(zldis / z0h[n]) -
                stability_func2(zeta) +
                stability_func2(z0h[n] / obu[n]))
        elseif zeta <= one(T)
            temp12m[n] = VKC / (log(zldis / z0h[n]) + T(5.0) * zeta - T(5.0) * z0h[n] / obu[n])
        else
            temp12m[n] = VKC / (log(obu[n] / z0h[n]) + T(5.0) - T(5.0) * z0h[n] / obu[n] +
                (T(5.0) * log(zeta) + zeta - one(T)))
        end

        # ------- Humidity profile at 2m -------
        if z0q[n] == z0h[n]
            temp22m[n] = temp12m[n]
        else
            zldis = T(2.0) + z0q[n]
            zeta = zldis / obu[n]
            if zeta < -zetat
                temp22m[n] = VKC / (log(-zetat * obu[n] / z0q[n]) -
                    stability_func2(-zetat) +
                    stability_func2(z0q[n] / obu[n]) +
                    T(0.8) * ((zetat)^(-T(0.333)) - (-zeta)^(-T(0.333))))
            elseif zeta < zero(T)
                temp22m[n] = VKC / (log(zldis / z0q[n]) -
                    stability_func2(zeta) +
                    stability_func2(z0q[n] / obu[n]))
            elseif zeta <= one(T)
                temp22m[n] = VKC / (log(zldis / z0q[n]) + T(5.0) * zeta - T(5.0) * z0q[n] / obu[n])
            else
                temp22m[n] = VKC / (log(obu[n] / z0q[n]) + T(5.0) - T(5.0) * z0q[n] / obu[n] +
                    (T(5.0) * log(zeta) + zeta - one(T)))
            end
        end

        # ------- 10-m wind for dust model -------
        zldis = forc_hgt_u_patch[n] - displa[n]
        zeta = zldis / obu[n]

        if min(zeta, one(T)) < zero(T)
            tmp1 = (one(T) - T(16.0) * min(zeta, one(T)))^T(0.25)
            tmp2 = log((one(T) + tmp1 * tmp1) / T(2.0))
            tmp3 = log((one(T) + tmp1) / T(2.0))
            fmnew = T(2.0) * tmp3 + tmp2 - T(2.0) * atan(tmp1) + T(1.5707963)
        else
            fmnew = -T(5.0) * min(zeta, one(T))
        end
        if iter == 1
            fm[n] = fmnew
        else
            fm[n] = T(0.5) * (fm[n] + fmnew)
        end

        zeta10 = min(T(10.0) / obu[n], one(T))
        if zeta == zero(T)
            zeta10 = zero(T)
        end
        if zeta10 < zero(T)
            tmp1 = (one(T) - T(16.0) * zeta10)^T(0.25)
            tmp2 = log((one(T) + tmp1 * tmp1) / T(2.0))
            tmp3 = log((one(T) + tmp1) / T(2.0))
            fm10 = T(2.0) * tmp3 + tmp2 - T(2.0) * atan(tmp1) + T(1.5707963)
        else
            fm10 = -T(5.0) * zeta10
        end

        tmp4 = log(max(one(T), forc_hgt_u_patch[n] / T(10.0)))

        u10_patch[n] = ur[n] - ustar[n] / VKC * (tmp4 - fm[n] + fm10)
        fv_patch[n] = ustar[n]
    end
end

"""
    fv_profile_update!(fv, fn, filtern, displa, z0m, z0h, z0q, obu, iter, ur, um,
                       ustar, temp1, temp2, temp12m, temp22m, fm)

Launch the kernelized default (non-landunit) Monin-Obukhov wind/temperature/
humidity profile update over the `fn` filtered patches. One thread per filtered
patch; backend-agnostic (CPU loop or GPU). Implements the `landunit_index=false`
path of `friction_velocity!`.
"""
function fv_profile_update!(fv::FrictionVelocityData, fn::Int,
        filtern::AbstractVector{<:Integer}, displa, z0m, z0h, z0q, obu, iter::Int, ur, um,
        ustar, temp1, temp2, temp12m, temp22m, fm; active = nothing)
    # Transition points of the flux-gradient relation, eltype-converted (with VKC)
    # so the kernel lowers to valid Metal IR under Float32; byte-identical on CPU.
    T = eltype(ustar)
    zetam = T(1.574)  # wind profile
    zetat = T(0.465)  # temperature profile
    # Active mask: when unset (single-pass callers), all filtered patches are
    # active. The canopy Newton loop passes a persistent mask so converged
    # patches are skipped — fm[n] is non-idempotent (fm = 0.5*(fm+fmnew)).
    act = active === nothing ? fill(true, length(ustar)) : active
    _launch!(_fv_profile_kernel!, ustar,
        temp1, temp2, temp12m, temp22m, fm,
        fv.vds_patch, fv.u10_clm_patch, fv.va_patch, fv.u10_patch, fv.fv_patch,
        filtern, act, displa, z0m, z0h, z0q, obu, ur, um,
        fv.forc_hgt_u_patch, fv.forc_hgt_t_patch, fv.forc_hgt_q_patch,
        iter, T(VKC), zetam, zetat;
        ndrange = fn)
    return nothing
end
