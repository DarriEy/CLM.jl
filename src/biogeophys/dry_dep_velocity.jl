# ==========================================================================
# Ported from: src/biogeophys/DryDepVelocity.F90
# Dry deposition velocity calculation using a resistance-in-series approach.
#
# This module calculates dry deposition velocities for trace gas species
# using the Wesely (1989) resistance model:
#   V_d = 1 / (R_a + R_b + R_c)
# where:
#   R_a = aerodynamic resistance (from the friction velocity module)
#   R_b = quasi-laminar boundary layer resistance
#   R_c = surface (canopy) resistance composed of:
#         - stomatal resistance (R_s)
#         - mesophyll resistance (R_m)
#         - cuticular resistance (R_lu)
#         - in-canopy aerodynamic resistance (R_ac)
#         - ground surface resistance (R_gs)
#
# References:
#   Wesely, M.L. (1989). Parameterization of surface resistances to
#     gaseous dry deposition in regional-scale numerical models.
#     Atmos. Environ., 23, 1293-1304.
#   Emmons et al. (2010) Description and evaluation of the Model for
#     Ozone and Related chemical Tracers, version 4 (MOZART-4).
#     Geosci. Model Dev., 3, 43-67.
#
# Public functions:
#   drydep_init!        -- Allocate and initialize DryDepVelocityData
#   drydep_clean!       -- Deallocate all fields
#   depvel_compute!     -- Calculate dry deposition velocity for all species
#
# ==========================================================================

# =========================================================================
# Dry deposition velocity constants
# =========================================================================

# Number of Wesely landuse types used in the deposition parameterization
const N_LAND_TYPE = 11

# Season indices
const MIDSUMMER_LUSH_VEG = 1
const AUTUMN_UNHARVESTED  = 2
const LATE_AUTUMN_LEAFLESS = 3
const WINTER_SNOW_ICE      = 4
const TRANSITIONAL_SPRING  = 5

# Reactivity factors (Wesely Table 2)
# f0: reactivity factor for oxidation of biological substances
# Used in computing mesophyll and cuticular resistances
const REACT_FACT_SO2  = 1.0   # SO2 reactivity factor (f0)
const REACT_FACT_O3   = 1.0   # O3 reactivity factor (f0)

# Diffusivity of H2O in air at 25 deg C [cm^2/s]
const DH2O = 0.25

# =========================================================================
# Seasonal resistance table data (Wesely 1989, Table 3)
# Dimensions: (n_landtype, n_seasons) - values in s/m
# =========================================================================

# Minimum stomatal resistance (R_i) [s/m]
const RI_TABLE = [
    9999.0  9999.0  9999.0  9999.0  9999.0   # 1: Urban land
    60.0    9999.0  9999.0  9999.0  120.0     # 2: Agricultural land
    120.0   200.0   9999.0  9999.0  240.0     # 3: Range land
    70.0    150.0   9999.0  9999.0  130.0     # 4: Deciduous forest
    130.0   200.0   9999.0  9999.0  250.0     # 5: Coniferous forest
    100.0   150.0   9999.0  9999.0  190.0     # 6: Mixed forest incl wetland
    9999.0  9999.0  9999.0  9999.0  9999.0    # 7: Water (inland + ocean)
    9999.0  9999.0  9999.0  9999.0  9999.0    # 8: Barren land (desert)
    9999.0  9999.0  9999.0  9999.0  9999.0    # 9: Non-forested wetland
    80.0    100.0   9999.0  9999.0  160.0     # 10: Mixed agricultural + range
    60.0    150.0   9999.0  9999.0  120.0     # 11: Rocky open with low growth
]

# Ground surface resistance for SO2 [s/m]
const RGS_SO2_TABLE = [
    400.0  150.0  350.0  500.0  200.0     # 1
    200.0  150.0  200.0  100.0  200.0     # 2
    350.0  200.0  200.0  100.0  300.0     # 3
    500.0  200.0  200.0  100.0  300.0     # 4
    500.0  200.0  200.0  100.0  300.0     # 5
    100.0  100.0  200.0  100.0  200.0     # 6
    0.0    1000.0 0.0    0.0    0.0       # 7
    1000.0 400.0  400.0  400.0  600.0     # 8
    0.0    1000.0 0.0    0.0    0.0       # 9
    300.0  150.0  200.0  100.0  200.0     # 10
    400.0  200.0  200.0  100.0  300.0     # 11
]

# Ground surface resistance for O3 [s/m]
const RGS_O3_TABLE = [
    300.0  150.0  200.0  200.0  200.0     # 1
    200.0  150.0  200.0  200.0  200.0     # 2
    200.0  200.0  200.0  200.0  200.0     # 3
    200.0  200.0  200.0  400.0  200.0     # 4
    200.0  200.0  200.0  400.0  200.0     # 5
    200.0  200.0  200.0  400.0  200.0     # 6
    2000.0 1000.0 2000.0 2000.0 2000.0    # 7
    400.0  200.0  200.0  200.0  300.0     # 8
    2000.0 1000.0 2000.0 2000.0 2000.0    # 9
    200.0  150.0  200.0  200.0  200.0     # 10
    200.0  200.0  200.0  200.0  200.0     # 11
]

# Cuticular resistance (R_lu) [s/m]
const RLU_TABLE = [
    9999.0  9999.0  9999.0  9999.0  9999.0    # 1
    2000.0  9000.0  9999.0  9999.0  4000.0    # 2
    2000.0  9000.0  9999.0  9999.0  4000.0    # 3
    2000.0  4000.0  9999.0  9999.0  8000.0    # 4
    2000.0  4000.0  9999.0  9999.0  8000.0    # 5
    2000.0  4000.0  9999.0  9999.0  8000.0    # 6
    9999.0  9999.0  9999.0  9999.0  9999.0    # 7
    9999.0  9999.0  9999.0  9999.0  9999.0    # 8
    9999.0  9999.0  9999.0  9999.0  9999.0    # 9
    2000.0  9000.0  9999.0  9999.0  4000.0    # 10
    4000.0  9000.0  9999.0  9999.0  8000.0    # 11
]

# In-canopy aerodynamic resistance (R_ac) [s/m]
const RAC_TABLE = [
    100.0   100.0   100.0   100.0   100.0     # 1
    200.0   150.0   100.0   100.0   150.0     # 2
    100.0   100.0   100.0   100.0   100.0     # 3
    2000.0  1500.0  1000.0  1000.0  1500.0    # 4
    2000.0  1500.0  1000.0  1000.0  1500.0    # 5
    2000.0  1500.0  1000.0  1000.0  1500.0    # 6
    0.0     0.0     0.0     0.0     0.0       # 7
    0.0     0.0     0.0     0.0     0.0       # 8
    0.0     0.0     0.0     0.0     0.0       # 9
    200.0   150.0   100.0   100.0   150.0     # 10
    100.0   100.0   100.0   100.0   100.0     # 11
]

# Canopy height [m] for each land type
const Z0_TABLE = [
    1.0, 0.25, 0.05, 1.0, 1.0, 1.0, 0.0006, 0.002, 0.15, 0.3, 0.1
]

# =========================================================================
# Data type
# =========================================================================

"""
    DryDepVelocityData

Dry deposition velocity data structure. Holds per-patch deposition velocities
and per-species parameters needed for the Wesely resistance model.

Ported from module-level arrays in `DryDepVelocity.F90`.
"""
Base.@kwdef mutable struct DryDepVelocityData{FT<:Real,
                                  V<:AbstractVector{FT},
                                  M<:AbstractMatrix{FT},
                                  VI<:AbstractVector{<:Integer}}
    # --- Configuration ---
    n_drydep::Int = 0                     # number of dry deposition species

    # --- Per-species properties (length n_drydep) ---
    foxd::V = Float64[]     # reactivity factor for oxidation [0-1]
    dv::V = Float64[]     # diffusivity in air [cm^2/s]
    mapping::VI  = Int[]         # mapping from species to Wesely land type category

    # --- Per-patch output (np x n_drydep) ---
    velocity_patch::M = Matrix{Float64}(undef, 0, 0)  # dry deposition velocity [cm/s]
end

DryDepVelocityData{FT}(; kwargs...) where {FT<:Real} =
    DryDepVelocityData{FT, Vector{FT}, Matrix{FT}, Vector{Int}}(; kwargs...)
Adapt.@adapt_structure DryDepVelocityData


# =========================================================================
# Allocation / initialization
# =========================================================================

"""
    drydep_init!(dd::DryDepVelocityData, np::Int, n_drydep::Int;
                 foxd::Vector{<:Real}=Float64[],
                 dv::Vector{<:Real}=Float64[],
                 mapping::Vector{Int}=Int[])

Allocate and initialize all fields of a `DryDepVelocityData` instance for
`np` patches and `n_drydep` dry deposition species.

Arguments:
- `foxd`: reactivity factor for oxidation for each species (0-1)
- `dv`: molecular diffusivity in air for each species [cm^2/s]
- `mapping`: Wesely land type index for each species

Ported from `DryDepVelocity%Init` in `DryDepVelocity.F90`.
"""
function drydep_init!(dd::DryDepVelocityData, np::Int, n_drydep::Int;
                      foxd::Vector{<:Real} = Float64[],
                      dv::Vector{<:Real} = Float64[],
                      mapping::Vector{Int} = Int[])
    dd.n_drydep = n_drydep

    if n_drydep > 0
        FT = eltype(dd.foxd)
        if isempty(foxd)
            dd.foxd = zeros(FT, n_drydep)
        else
            dd.foxd = copy(foxd)
        end

        if isempty(dv)
            dd.dv = fill(FT(0.2), n_drydep)  # default ~0.2 cm^2/s
        else
            dd.dv = copy(dv)
        end

        if isempty(mapping)
            dd.mapping = fill(1, n_drydep)  # default: urban
        else
            dd.mapping = copy(mapping)
        end

        dd.velocity_patch = fill(zero(FT), np, n_drydep)
    end

    return nothing
end

"""
    drydep_clean!(dd::DryDepVelocityData)

Deallocate (reset to empty) all fields of a `DryDepVelocityData` instance.
"""
function drydep_clean!(dd::DryDepVelocityData{FT}) where {FT}
    dd.n_drydep = 0
    dd.foxd = FT[]
    dd.dv = FT[]
    dd.mapping = Int[]
    dd.velocity_patch = Matrix{FT}(undef, 0, 0)
    return nothing
end

# =========================================================================
# Season determination
# =========================================================================

"""
    _wesely_season(lat::Float64, month::Int) -> Int

Determine the Wesely season index (1-5) based on latitude and month.
Follows the seasonal mapping used in MOZART / CLM.

Returns one of:
  1 = Midsummer with lush vegetation
  2 = Autumn with unharvested cropland
  3 = Late autumn with leafless deciduous
  4 = Winter with snow/ice
  5 = Transitional spring
"""
function _wesely_season(lat::Real, month::Int)
    # Northern hemisphere months -> seasons
    #   Dec/Jan/Feb -> 4 (winter)
    #   Mar/Apr     -> 5 (transitional spring)
    #   May/Jun/Jul/Aug -> 1 (midsummer)
    #   Sep         -> 2 (autumn)
    #   Oct/Nov     -> 3 (late autumn)
    nh_season = [4, 4, 5, 5, 1, 1, 1, 1, 2, 3, 3, 4]

    if lat >= 0.0
        return nh_season[month]
    else
        # Southern hemisphere: shift by 6 months
        sh_month = mod(month + 5, 12) + 1
        return nh_season[sh_month]
    end
end

"""
    _wesely_index_season(lun_type, is_urban, has_snow, elai, minlai, maxlai, mlaidiff)
        -> (wesveg_override::Int, index_season::Int)

Select the Wesely seasonal category (1-5) exactly as CTSM
`DryDepVelocity.F90:373-405` does — from the patch's annual min/max LAI
(`annlai`), the current exposed LAI, the month-to-month LAI difference
(`mlaidiff`), snow presence (`snow_depth(c) > 0`), and the landunit type — rather
than from latitude+month.

Returns `(wesveg_override, index_season)`. `wesveg_override == 0` means keep the
PFT-derived Wesely land type; a positive value is the landunit-type override CTSM
sets jointly with the season in the same block (ice→8, deep-lake→7, wetland→9,
urban→1). `index_season < 1` signals an indeterminate season (only reachable with
NaN inputs, e.g. `annlai`/`mlaidiff` never populated); the caller then falls back
to the latitude+month heuristic instead of aborting (Fortran `endrun`s).

Faithful transcription of the Fortran (soil landunit is the common case):

    minlai = minval(annlai(:,p)); maxlai = maxval(annlai(:,p)); index_season = -1
    if lun%itype /= istsoil:
        istice   -> wesveg=8, index_season=4
        istdlak  -> wesveg=7, index_season=4
        istwet   -> wesveg=9, index_season=2
        urbpoi   -> wesveg=1, index_season=2
        (else e.g. crop: leave index_season = -1)
    else if snow_depth(c) > 0:            index_season = 4
    else if elai > 0.5*maxlai:            index_season = 1
    if index_season<0 and elai < minlai+0.05*(maxlai-minlai): index_season = 3
    if index_season<0: mlaidiff>0 -> 2, <0 -> 5, ==0 -> 3
"""
@inline function _wesely_index_season(lun_type::Integer, is_urban::Bool, has_snow::Bool,
                                      elai::T, minlai::T, maxlai::T,
                                      mlaidiff::T) where {T<:Real}
    wesveg_ov = 0
    season = -1

    # Fortran: `if ( lun%itype(l) /= istsoil )` — NOTE only soil (not crop) skips
    # this branch, so a crop landunit falls through with index_season still -1.
    if lun_type != ISTSOIL
        if lun_type == ISTICE
            wesveg_ov = 8
            season = WINTER_SNOW_ICE
        elseif lun_type == ISTDLAK
            wesveg_ov = 7
            season = WINTER_SNOW_ICE
        elseif lun_type == ISTWET
            wesveg_ov = 9
            season = AUTUMN_UNHARVESTED
        elseif is_urban
            wesveg_ov = 1
            season = AUTUMN_UNHARVESTED
        end
        # any other non-soil landunit (e.g. crop): leave season = -1
    elseif has_snow
        season = WINTER_SNOW_ICE
    elseif elai > T(0.5) * maxlai
        season = MIDSUMMER_LUSH_VEG
    end

    if season < 0
        if elai < (minlai + T(0.05) * (maxlai - minlai))
            season = LATE_AUTUMN_LEAFLESS
        end
    end

    if season < 0
        if mlaidiff > zero(T)
            season = AUTUMN_UNHARVESTED
        elseif mlaidiff < zero(T)
            season = TRANSITIONAL_SPRING
        elseif mlaidiff == zero(T)
            season = LATE_AUTUMN_LEAFLESS
        end
    end

    return wesveg_ov, season
end

# =========================================================================
# Wesely land type mapping from CLM PFT
# =========================================================================

"""
    _pft_to_wesely(itype::Int) -> Int

Map a CLM PFT type index (0-based Fortran convention) to a Wesely
land use type index (1-11).

Mapping follows CLM conventions:
  0 (not vegetated / bare soil) -> 8 (barren land)
  1-3 (needleleaf trees)        -> 5 (coniferous forest)
  4-8 (broadleaf trees)         -> 4 (deciduous forest)
  9-11 (shrubs)                 -> 3 (range land)
  12-14 (grasses)               -> 3 (range land)
  15+ (crops)                   -> 2 (agricultural land)
"""
@inline function _pft_to_wesely(itype::Integer)
    if itype == 0
        return 8   # barren land
    elseif itype <= 3
        return 5   # coniferous forest
    elseif itype <= 8
        return 4   # deciduous forest
    elseif itype <= 14
        return 3   # range land
    else
        return 2   # agricultural land
    end
end

# =========================================================================
# Quasi-laminar boundary layer resistance
# =========================================================================

"""
    _calc_rb(ustar::Float64, dv_species::Float64) -> Float64

Calculate the quasi-laminar boundary layer resistance R_b [s/m].
Uses the formula from Wesely (1989):
    R_b = 2 / (kappa * u*) * (Sc/Pr)^(2/3)
where:
    Sc = nu / D_x    (Schmidt number)
    Pr = nu / D_H2O  (Prandtl number, ~0.72 for air)
    kappa = von Karman constant

`dv_species` is the molecular diffusivity [cm^2/s].
`ustar` is the friction velocity [m/s].

Returns R_b in s/m.
"""
@inline function _calc_rb(ustar::Real, dv_species::Real)
    # Working type carries the caller's precision (Float64 on host, Float32 on a
    # device kernel) so no Float64 literal reaches a Metal kernel; byte-identical
    # to the original on the Float64 host path (T(x)==x).
    T = promote_type(typeof(ustar), typeof(dv_species))
    # Avoid division by zero
    ustar_safe = max(T(ustar), T(0.001))

    # Schmidt number / Prandtl number ratio
    # D_H2O ~ 0.25 cm^2/s for water vapor at 25 C
    sc_over_pr = T(DH2O) / max(T(dv_species), T(1.0e-10))

    # R_b = 2/(kappa * u*) * (Sc/Pr)^(2/3)
    rb = (T(2.0) / (T(VKC) * ustar_safe)) * sc_over_pr^(T(2.0) / T(3.0))

    return rb
end

# =========================================================================
# Surface (canopy) resistance
# =========================================================================

"""
    _calc_rc(season::Int, lt::Int, lai::Float64,
             foxd_i::Float64, heff_i::Float64,
             t_sfc::Float64, solar_flux::Float64,
             has_snow::Bool) -> Float64

Calculate the surface (canopy) resistance R_c [s/m] for a single
species at a single point, following the Wesely (1989) parameterization.

Arguments:
- `season`: Wesely season index (1-5)
- `lt`: Wesely land type (1-11)
- `lai`: leaf area index
- `foxd_i`: reactivity factor for species
- `heff_i`: effective Henry's law constant [M/atm] (used scaled as foxd)
- `t_sfc`: surface temperature [K]
- `solar_flux`: incident PAR/solar flux [W/m^2]
- `has_snow`: true if snow is present

Returns R_c in s/m (minimum 1 s/m).
"""
function _calc_rc(season::Int, lt::Int, lai::Real,
                  foxd_i::Real, heff_i::Real,
                  t_sfc::Real, solar_flux::Real,
                  has_snow::Bool)
    # Host wrapper: look the season/land-type resistances up in the module tables,
    # then delegate to the device-callable arithmetic core.
    return _calc_rc_core(RI_TABLE[lt, season], RLU_TABLE[lt, season],
                         RAC_TABLE[lt, season], RGS_SO2_TABLE[lt, season],
                         RGS_O3_TABLE[lt, season],
                         lai, foxd_i, heff_i, t_sfc, solar_flux, has_snow)
end

"""
    _calc_rc_core(ri, rlu, rac, rgs_so2, rgs_o3, lai, foxd_i, heff_i,
                  t_sfc, solar_flux, has_snow) -> R_c

Device-callable surface-resistance arithmetic (the body of `_calc_rc`), taking the
five already-looked-up Wesely resistances as scalars. Generic in the working type
`T` with every literal carried at `T`, so it compiles on a Metal-Float32 kernel and
is byte-identical to the original on the Float64 host path (`T(x) == x`). The LAI
scale factor of the original was dead (computed, never applied) and is omitted.
"""
@inline function _calc_rc_core(ri::T, rlu::T, rac::T, rgs_so2::T, rgs_o3::T,
                               lai::T, foxd_i::T, heff_i::T,
                               t_sfc::T, solar_flux::T, has_snow::Bool) where {T<:Real}
    tfrz = T(TFRZ)

    # --- Stomatal resistance (light-dependent opening) ---
    if ri < T(9999.0) && solar_flux > one(T) && t_sfc > tfrz - T(5.0)
        rs = ri * (one(T) + (T(200.0) / (solar_flux + T(0.1)))^2) * (T(400.0) / (t_sfc - (tfrz - T(40.0))))
        rs = max(rs, ri)
    else
        rs = T(1.0e10)  # effectively infinite (stomata closed)
    end

    # --- Temperature correction (stomata close below 0 C) ---
    if t_sfc < tfrz
        rs = T(1.0e10)
    end

    # --- Mesophyll resistance ---
    rm_denom = heff_i / T(3000.0) + T(100.0) * foxd_i
    if rm_denom > T(1.0e-10)
        rm = one(T) / rm_denom
    else
        rm = T(1.0e10)
    end
    rm = max(rm, zero(T))

    # --- Cuticular resistance ---
    if rlu < T(9999.0)
        rlu_denom = T(1.0e-5) * heff_i + foxd_i
        if rlu_denom > T(1.0e-10)
            rlu_eff = rlu / rlu_denom
        else
            rlu_eff = T(1.0e10)
        end
        rlu_eff = max(rlu_eff, one(T))
    else
        rlu_eff = T(1.0e10)
    end

    # --- Ground surface resistance ---
    if rgs_so2 > zero(T) && rgs_o3 > zero(T)
        rgs_denom = heff_i / (T(1.0e5) * rgs_so2) + foxd_i / rgs_o3
        if rgs_denom > T(1.0e-10)
            rgs = one(T) / rgs_denom
        else
            rgs = T(1.0e10)
        end
        rgs = max(rgs, one(T))
    else
        # Water/wetland surface (a ground component is zero): use the non-zero part.
        rgs_denom = zero(T)
        if rgs_so2 > zero(T)
            rgs_denom += heff_i / (T(1.0e5) * rgs_so2)
        elseif heff_i > T(0.01)
            rgs_denom += heff_i / T(1.0e5)
        end
        if rgs_o3 > zero(T)
            rgs_denom += foxd_i / rgs_o3
        end
        if rgs_denom > T(1.0e-10)
            rgs = one(T) / rgs_denom
        else
            rgs = T(1000.0)
        end
        rgs = max(rgs, one(T))
    end

    # --- In-canopy aerodynamic resistance ---
    rac_eff = max(rac, zero(T))

    # --- Snow correction ---
    if has_snow
        rgs = max(rgs, T(500.0))
    end

    # --- Combine: stomatal + mesophyll in series ---
    if rs < T(1.0e9)
        r_stom = rs + rm
    else
        r_stom = T(1.0e10)
    end

    # Upper canopy: parallel of stomatal and cuticular pathways
    if r_stom < T(1.0e9) && rlu_eff < T(1.0e9)
        r_upper = one(T) / (one(T) / r_stom + one(T) / rlu_eff)
    elseif r_stom < T(1.0e9)
        r_upper = r_stom
    elseif rlu_eff < T(1.0e9)
        r_upper = rlu_eff
    else
        r_upper = T(1.0e10)
    end

    # Lower canopy (in-canopy transport + ground)
    r_lower = rac_eff + rgs

    # Total canopy resistance: parallel of upper and lower
    if r_upper < T(1.0e9) && r_lower > zero(T)
        rc = one(T) / (one(T) / r_upper + one(T) / r_lower)
    elseif r_upper < T(1.0e9)
        rc = r_upper
    elseif r_lower > zero(T)
        rc = r_lower
    else
        rc = one(T)
    end

    return max(rc, one(T))
end

# =========================================================================
# Main dry deposition velocity calculation
# =========================================================================

"""
    depvel_compute!(dd::DryDepVelocityData,
                    mask_patch::AbstractVector{Bool},
                    bounds_patch::UnitRange{Int},
                    patch_gridcell::Vector{Int},
                    patch_column::Vector{Int},
                    patch_landunit::Vector{Int},
                    patch_itype::Vector{Int},
                    ram1_patch::Vector{<:Real},
                    rb1_patch::Vector{<:Real},
                    fv_patch::Vector{<:Real},
                    elai_patch::Vector{<:Real},
                    forc_t_downscaled_col::Vector{<:Real},
                    forc_solar_col::Vector{<:Real},
                    frac_sno::Vector{<:Real},
                    lat_grc::Vector{<:Real},
                    month::Int;
                    heff::Vector{<:Real}=Float64[])

Calculate dry deposition velocities for all species at all active patches.
Results are stored in `dd.velocity_patch` in cm/s.

The deposition velocity is:
    V_d = 1 / (R_a + R_b + R_c)

where:
- R_a: aerodynamic resistance = ram1_patch [s/m]
- R_b: quasi-laminar boundary layer resistance [s/m]
- R_c: surface (canopy) resistance [s/m]

Arguments:
- `dd`: DryDepVelocityData instance
- `mask_patch`: BitVector mask for active non-lake patches
- `bounds_patch`: range of patch indices
- `patch_gridcell`: gridcell index for each patch
- `patch_column`: column index for each patch
- `patch_landunit`: landunit index for each patch
- `patch_itype`: PFT type for each patch (0-based Fortran convention)
- `ram1_patch`: aerodynamic resistance [s/m]
- `rb1_patch`: boundary layer resistance [s/m] (used as backup)
- `fv_patch`: friction velocity [m/s]
- `elai_patch`: exposed LAI [-]
- `forc_t_downscaled_col`: surface temperature [K] (column-level)
- `forc_solar_col`: incident solar radiation [W/m^2] (column-level)
- `frac_sno`: fraction of ground covered by snow (column-level)
- `lat_grc`: latitude [degrees] (gridcell-level)
- `month`: current month (1-12) — only used for the degenerate season fallback
- `heff`: effective Henry's law constant [M/atm] for each species (optional)

Keyword inputs for the Fortran-faithful Wesely `index_season` selection
(`DryDepVelocity.F90:373-405`). When `annlai_patch` is non-empty the season is
picked per-patch from the annual min/max LAI rather than latitude+month; when
omitted the kernel keeps the latitude+month heuristic (legacy behavior):
- `annlai_patch`: 12 months of monthly LAI per patch (np x 12) — min/max drive the season
- `mlaidiff_patch`: month-to-month LAI difference per patch (transitional-season sign)
- `lun_itype`: landunit-level type (ISTSOIL/ISTICE/ISTDLAK/ISTWET/urban) per landunit
- `snow_depth_col`: snow depth [m] (column-level); `> 0` forces the winter season

Ported from `DryDepVelocity` subroutine in `DryDepVelocity.F90`.
"""
# --- Kernel: one thread per patch; inner loop over species. Every write is to the
# patch's own (p, i) slot (distinct per species -> race-free, no fixed-index +=).
# When `use_annlai` is true the Wesely season is picked per-patch from the annual
# min/max LAI (`_wesely_index_season`, matching DryDepVelocity.F90:373-405); the
# lat+month heuristic (precomputed for both hemispheres as two Int scalars) is kept
# only as the degenerate fallback (NaN inputs) and for legacy callers that supply
# no annual-LAI phenology.
@kernel function _depvel_kernel!(velocity_patch, @Const(mask_patch),
        @Const(patch_gridcell), @Const(patch_column), @Const(patch_itype),
        @Const(ram1_patch), @Const(fv_patch), @Const(elai_patch),
        @Const(forc_t_col), @Const(forc_solar_col), @Const(frac_sno), @Const(lat_grc),
        @Const(dv), @Const(foxd), @Const(heff_vals),
        @Const(ri_tab), @Const(rlu_tab), @Const(rac_tab), @Const(rgs_so2_tab), @Const(rgs_o3_tab),
        @Const(annlai), @Const(mlaidiff), @Const(snow_depth), @Const(patch_lun_type),
        lo::Int, hi::Int, n_drydep::Int, season_nh::Int, season_sh::Int, use_annlai::Bool)
    p = @index(Global)
    @inbounds if lo <= p <= hi && mask_patch[p]
        T = eltype(velocity_patch)
        g = patch_gridcell[p]
        c = patch_column[p]

        t_sfc = forc_t_col[c]
        solar = forc_solar_col[c]
        has_snow = frac_sno[c] > T(0.5)
        lat = lat_grc[g]
        lai = max(elai_patch[p], zero(T))

        if use_annlai
            # Fortran-faithful index_season from the patch's annual min/max LAI.
            lun_type = patch_lun_type[p]
            is_urb = (ISTURB_MIN <= lun_type) & (lun_type <= ISTURB_MAX)
            snow_s = snow_depth[c] > zero(T)
            minlai = annlai[p, 1]
            maxlai = annlai[p, 1]
            for k in 2:12
                v = annlai[p, k]
                minlai = min(minlai, v)
                maxlai = max(maxlai, v)
            end
            wesveg_ov, season0 = _wesely_index_season(lun_type, is_urb, snow_s,
                                                      lai, minlai, maxlai, mlaidiff[p])
            if season0 < 1
                # indeterminate (NaN inputs) -> lat+month fallback
                season = lat >= zero(T) ? season_nh : season_sh
                lt = _pft_to_wesely(patch_itype[p])
            else
                season = season0
                lt = wesveg_ov > 0 ? wesveg_ov : _pft_to_wesely(patch_itype[p])
            end
        else
            season = lat >= zero(T) ? season_nh : season_sh   # == _wesely_season(lat, month)
            lt = _pft_to_wesely(patch_itype[p])
        end

        ra = max(ram1_patch[p], one(T))
        ustar = max(fv_patch[p], T(0.001))

        for i in 1:n_drydep
            rb = _calc_rb(ustar, dv[i])
            rc = _calc_rc_core(ri_tab[lt, season], rlu_tab[lt, season],
                               rac_tab[lt, season], rgs_so2_tab[lt, season],
                               rgs_o3_tab[lt, season],
                               lai, foxd[i], heff_vals[i], t_sfc, solar, has_snow)
            vd = one(T) / (ra + rb + rc)
            velocity_patch[p, i] = vd * T(100.0)
        end
    end
end

function depvel_compute!(dd::DryDepVelocityData,
                         mask_patch::AbstractVector{Bool},
                         bounds_patch::UnitRange{Int},
                         patch_gridcell::AbstractVector{<:Integer},
                         patch_column::AbstractVector{<:Integer},
                         patch_landunit::AbstractVector{<:Integer},
                         patch_itype::AbstractVector{<:Integer},
                         ram1_patch::AbstractVector{<:Real},
                         rb1_patch::AbstractVector{<:Real},
                         fv_patch::AbstractVector{<:Real},
                         elai_patch::AbstractVector{<:Real},
                         forc_t_downscaled_col::AbstractVector{<:Real},
                         forc_solar_col::AbstractVector{<:Real},
                         frac_sno::AbstractVector{<:Real},
                         lat_grc::AbstractVector{<:Real},
                         month::Int;
                         heff::AbstractVector{<:Real} = Float64[],
                         annlai_patch::AbstractMatrix{<:Real} = Matrix{Float64}(undef, 0, 0),
                         mlaidiff_patch::AbstractVector{<:Real} = Float64[],
                         lun_itype::AbstractVector{<:Integer} = Int[],
                         snow_depth_col::AbstractVector{<:Real} = Float64[])
    n_drydep = dd.n_drydep
    if n_drydep == 0
        return nothing
    end
    isempty(bounds_patch) && return nothing

    FT = eltype(dd.velocity_patch)
    np = size(dd.velocity_patch, 1)
    # Default effective Henry's law constants (small default), on host.
    heff_vals_h = isempty(heff) ? fill(FT(1.0e-2), n_drydep) : heff

    # Wesely season fallback: a function of lat SIGN and the global `month`;
    # precompute both hemispheres on the host. Used for the degenerate
    # (indeterminate-season) case and for legacy callers that pass no annual LAI.
    season_nh = _wesely_season(1.0, month)   # lat >= 0
    season_sh = _wesely_season(-1.0, month)  # lat < 0

    # Fortran-faithful annlai-based index_season when the annual-LAI phenology is
    # supplied (the live driver passes annlai/mlaidiff/lun_itype/snow_depth);
    # otherwise the kernel keeps the latitude+month heuristic.
    use_annlai = !isempty(annlai_patch)
    if use_annlai
        # Gather the landunit TYPE per patch: lun_itype is landunit-level, and
        # patch_landunit is the 1-based landunit index of each patch.
        patch_lun_type_h = Vector{Int}(undef, np)
        @inbounds for pp in 1:np
            li = Int(patch_landunit[pp])
            patch_lun_type_h[pp] = (1 <= li <= length(lun_itype)) ? Int(lun_itype[li]) : ISTSOIL
        end
        annlai_h   = annlai_patch
        mlaidiff_h = mlaidiff_patch
        snowdep_h  = snow_depth_col
    else
        # Unused when use_annlai=false, but must be valid device arrays.
        patch_lun_type_h = fill(ISTSOIL, np)
        annlai_h   = zeros(FT, np, 12)
        mlaidiff_h = zeros(FT, np)
        snowdep_h  = frac_sno
    end

    # Move host constants (resistance tables + species params) and any host-resident
    # forcing/index arrays onto the state's backend + precision. All no-ops when
    # `dd.velocity_patch` is a host Array, so the host path stays byte-identical.
    p = dd.velocity_patch
    ri_tab  = _to_backend_like(p, FT, RI_TABLE)
    rlu_tab = _to_backend_like(p, FT, RLU_TABLE)
    rac_tab = _to_backend_like(p, FT, RAC_TABLE)
    rgs_so2 = _to_backend_like(p, FT, RGS_SO2_TABLE)
    rgs_o3  = _to_backend_like(p, FT, RGS_O3_TABLE)
    heff_v  = _to_backend_like(p, FT, heff_vals_h)
    pgrid   = _to_backend_like(p, FT, patch_gridcell)
    pcol    = _to_backend_like(p, FT, patch_column)
    pitype  = _to_backend_like(p, FT, patch_itype)
    ram1    = _to_backend_like(p, FT, ram1_patch)
    fv      = _to_backend_like(p, FT, fv_patch)
    elai    = _to_backend_like(p, FT, elai_patch)
    ft_col  = _to_backend_like(p, FT, forc_t_downscaled_col)
    fsol    = _to_backend_like(p, FT, forc_solar_col)
    fsno    = _to_backend_like(p, FT, frac_sno)
    lat     = _to_backend_like(p, FT, lat_grc)
    maskd   = _to_backend_like(p, FT, mask_patch)  # Bool -> integer overload keeps Bool
    annlai_b   = _to_backend_like(p, FT, annlai_h)
    mlaidiff_b = _to_backend_like(p, FT, mlaidiff_h)
    snowdep_b  = _to_backend_like(p, FT, snowdep_h)
    plun_b     = _to_backend_like(p, FT, patch_lun_type_h)  # Int -> integer overload keeps Int

    _launch!(_depvel_kernel!, dd.velocity_patch, maskd,
        pgrid, pcol, pitype, ram1, fv, elai, ft_col, fsol, fsno, lat,
        dd.dv, dd.foxd, heff_v,
        ri_tab, rlu_tab, rac_tab, rgs_so2, rgs_o3,
        annlai_b, mlaidiff_b, snowdep_b, plun_b,
        first(bounds_patch), last(bounds_patch), n_drydep, season_nh, season_sh, use_annlai;
        ndrange = size(dd.velocity_patch, 1))

    return nothing
end

# =========================================================================
# Weighted average to gridcell (for coupling to atmosphere)
# =========================================================================

"""
    drydep_p2g!(dd::DryDepVelocityData,
                ddvel_grc::Matrix{<:Real},
                bounds_patch::UnitRange{Int},
                mask_patch::AbstractVector{Bool},
                patch_gridcell::Vector{Int},
                wtgcell::Vector{<:Real},
                ng::Int)

Average patch-level dry deposition velocities to gridcell level using
the gridcell weights. Results stored in `ddvel_grc` (ng x n_drydep) [cm/s].

Ported from `DryDepVelocity%p2g` in `DryDepVelocity.F90`.
"""
function drydep_p2g!(dd::DryDepVelocityData,
                     ddvel_grc::Matrix{<:Real},
                     bounds_patch::UnitRange{Int},
                     mask_patch::AbstractVector{Bool},
                     patch_gridcell::Vector{Int},
                     wtgcell::Vector{<:Real},
                     ng::Int)
    n_drydep = dd.n_drydep
    if n_drydep == 0
        return nothing
    end

    # Zero out gridcell accumulation
    for g in 1:ng
        for i in 1:n_drydep
            ddvel_grc[g, i] = 0.0
        end
    end

    # Weighted average
    for p in bounds_patch
        mask_patch[p] || continue
        g = patch_gridcell[p]
        wt = wtgcell[p]

        for i in 1:n_drydep
            ddvel_grc[g, i] += dd.velocity_patch[p, i] * wt
        end
    end

    return nothing
end
