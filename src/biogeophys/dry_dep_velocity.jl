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
Base.@kwdef mutable struct DryDepVelocityData
    # --- Configuration ---
    n_drydep::Int = 0                     # number of dry deposition species

    # --- Per-species properties (length n_drydep) ---
    foxd::Vector{Float64} = Float64[]     # reactivity factor for oxidation [0-1]
    dv::Vector{Float64}   = Float64[]     # diffusivity in air [cm^2/s]
    mapping::Vector{Int}  = Int[]         # mapping from species to Wesely land type category

    # --- Per-patch output (np x n_drydep) ---
    velocity_patch::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # dry deposition velocity [cm/s]
end

# =========================================================================
# Allocation / initialization
# =========================================================================

"""
    drydep_init!(dd::DryDepVelocityData, np::Int, n_drydep::Int;
                 foxd::Vector{Float64}=Float64[],
                 dv::Vector{Float64}=Float64[],
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
                      foxd::Vector{Float64} = Float64[],
                      dv::Vector{Float64} = Float64[],
                      mapping::Vector{Int} = Int[])
    dd.n_drydep = n_drydep

    if n_drydep > 0
        if isempty(foxd)
            dd.foxd = zeros(n_drydep)
        else
            dd.foxd = copy(foxd)
        end

        if isempty(dv)
            dd.dv = fill(0.2, n_drydep)  # default ~0.2 cm^2/s
        else
            dd.dv = copy(dv)
        end

        if isempty(mapping)
            dd.mapping = fill(1, n_drydep)  # default: urban
        else
            dd.mapping = copy(mapping)
        end

        dd.velocity_patch = fill(0.0, np, n_drydep)
    end

    return nothing
end

"""
    drydep_clean!(dd::DryDepVelocityData)

Deallocate (reset to empty) all fields of a `DryDepVelocityData` instance.
"""
function drydep_clean!(dd::DryDepVelocityData)
    dd.n_drydep = 0
    dd.foxd = Float64[]
    dd.dv = Float64[]
    dd.mapping = Int[]
    dd.velocity_patch = Matrix{Float64}(undef, 0, 0)
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
function _wesely_season(lat::Float64, month::Int)
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
function _pft_to_wesely(itype::Int)
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
function _calc_rb(ustar::Float64, dv_species::Float64)
    # Avoid division by zero
    ustar_safe = max(ustar, 0.001)

    # Schmidt number / Prandtl number ratio
    # D_H2O ~ 0.25 cm^2/s for water vapor at 25 C
    sc_over_pr = DH2O / max(dv_species, 1.0e-10)

    # R_b = 2/(kappa * u*) * (Sc/Pr)^(2/3)
    rb = (2.0 / (VKC * ustar_safe)) * sc_over_pr^(2.0 / 3.0)

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
function _calc_rc(season::Int, lt::Int, lai::Float64,
                  foxd_i::Float64, heff_i::Float64,
                  t_sfc::Float64, solar_flux::Float64,
                  has_snow::Bool)

    # Look up season-dependent resistances from tables
    ri  = RI_TABLE[lt, season]
    rlu = RLU_TABLE[lt, season]
    rac = RAC_TABLE[lt, season]
    rgs_so2 = RGS_SO2_TABLE[lt, season]
    rgs_o3  = RGS_O3_TABLE[lt, season]

    # --- Stomatal resistance ---
    # Adjust for solar radiation (stomata open with light)
    if ri < 9999.0 && solar_flux > 1.0 && t_sfc > TFRZ - 5.0
        # Light-dependent stomatal opening
        rs = ri * (1.0 + (200.0 / (solar_flux + 0.1))^2) * (400.0 / (t_sfc - (TFRZ - 40.0)))
        rs = max(rs, ri)
    else
        rs = 1.0e10  # effectively infinite (stomata closed)
    end

    # --- Temperature correction for stomatal resistance ---
    # Stomata close below 0 C and above ~40 C
    if t_sfc < TFRZ
        rs = 1.0e10
    end

    # --- Mesophyll resistance ---
    # R_m = 1 / (H*/3000 + 100*f0)
    rm_denom = heff_i / 3000.0 + 100.0 * foxd_i
    if rm_denom > 1.0e-10
        rm = 1.0 / rm_denom
    else
        rm = 1.0e10  # effectively infinite for unreactive species
    end
    rm = max(rm, 0.0)

    # --- Cuticular resistance ---
    # R_lu for O3/SO2 reactivity scaling
    if rlu < 9999.0
        rlu_denom = 1.0e-5 * heff_i + foxd_i
        if rlu_denom > 1.0e-10
            rlu_eff = rlu / rlu_denom
        else
            rlu_eff = 1.0e10
        end
        rlu_eff = max(rlu_eff, 1.0)
    else
        rlu_eff = 1.0e10
    end

    # --- Ground surface resistance ---
    # Combine SO2 and O3 ground resistances based on species properties
    if rgs_so2 > 0.0 && rgs_o3 > 0.0
        rgs_denom = heff_i / (1.0e5 * rgs_so2) + foxd_i / rgs_o3
        if rgs_denom > 1.0e-10
            rgs = 1.0 / rgs_denom
        else
            rgs = 1.0e10  # unreactive species
        end
        rgs = max(rgs, 1.0)
    else
        # Water or wetland surface (rgs_so2 or rgs_o3 is zero)
        # Use the non-zero component if available
        rgs_denom = 0.0
        if rgs_so2 > 0.0
            rgs_denom += heff_i / (1.0e5 * rgs_so2)
        elseif heff_i > 0.01
            rgs_denom += heff_i / 1.0e5  # water surface dissolves soluble gases
        end
        if rgs_o3 > 0.0
            rgs_denom += foxd_i / rgs_o3
        end
        if rgs_denom > 1.0e-10
            rgs = 1.0 / rgs_denom
        else
            rgs = 1000.0  # default for unreactive species on water
        end
        rgs = max(rgs, 1.0)
    end

    # --- In-canopy aerodynamic resistance ---
    rac_eff = max(rac, 0.0)

    # --- Snow correction ---
    # Under snow conditions, increase ground resistance
    if has_snow
        rgs = max(rgs, 500.0)
    end

    # --- Combine resistances ---
    # Canopy pathway: stomatal + mesophyll in series
    if rs < 1.0e9
        r_stom = rs + rm
    else
        r_stom = 1.0e10
    end

    # Upper canopy: parallel of stomatal and cuticular pathways
    if r_stom < 1.0e9 && rlu_eff < 1.0e9
        r_upper = 1.0 / (1.0 / r_stom + 1.0 / rlu_eff)
    elseif r_stom < 1.0e9
        r_upper = r_stom
    elseif rlu_eff < 1.0e9
        r_upper = rlu_eff
    else
        r_upper = 1.0e10
    end

    # Lower canopy (in-canopy transport + ground)
    r_lower = rac_eff + rgs

    # Total canopy resistance: parallel of upper canopy and lower canopy
    if r_upper < 1.0e9 && r_lower > 0.0
        rc = 1.0 / (1.0 / r_upper + 1.0 / r_lower)
    elseif r_upper < 1.0e9
        rc = r_upper
    elseif r_lower > 0.0
        rc = r_lower
    else
        rc = 1.0
    end

    # Scale by LAI effect: reduce rc when LAI is large
    # (more leaf area = more deposition surface)
    if lai > 0.0 && lai < 20.0
        lai_factor = max(0.5, 1.0 - 0.02 * lai)
        # don't scale rc by LAI factor for now -- the Wesely tables
        # already account for typical LAI by season
    end

    return max(rc, 1.0)
end

# =========================================================================
# Main dry deposition velocity calculation
# =========================================================================

"""
    depvel_compute!(dd::DryDepVelocityData,
                    mask_patch::BitVector,
                    bounds_patch::UnitRange{Int},
                    patch_gridcell::Vector{Int},
                    patch_column::Vector{Int},
                    patch_landunit::Vector{Int},
                    patch_itype::Vector{Int},
                    ram1_patch::Vector{Float64},
                    rb1_patch::Vector{Float64},
                    fv_patch::Vector{Float64},
                    elai_patch::Vector{Float64},
                    forc_t_downscaled_col::Vector{Float64},
                    forc_solar_col::Vector{Float64},
                    frac_sno::Vector{Float64},
                    lat_grc::Vector{Float64},
                    month::Int;
                    heff::Vector{Float64}=Float64[])

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
- `month`: current month (1-12)
- `heff`: effective Henry's law constant [M/atm] for each species (optional)

Ported from `DryDepVelocity` subroutine in `DryDepVelocity.F90`.
"""
function depvel_compute!(dd::DryDepVelocityData,
                         mask_patch::BitVector,
                         bounds_patch::UnitRange{Int},
                         patch_gridcell::Vector{Int},
                         patch_column::Vector{Int},
                         patch_landunit::Vector{Int},
                         patch_itype::Vector{Int},
                         ram1_patch::Vector{Float64},
                         rb1_patch::Vector{Float64},
                         fv_patch::Vector{Float64},
                         elai_patch::Vector{Float64},
                         forc_t_downscaled_col::Vector{Float64},
                         forc_solar_col::Vector{Float64},
                         frac_sno::Vector{Float64},
                         lat_grc::Vector{Float64},
                         month::Int;
                         heff::Vector{Float64} = Float64[])
    n_drydep = dd.n_drydep
    if n_drydep == 0
        return nothing
    end

    # Default effective Henry's law constants
    heff_vals = if isempty(heff)
        fill(1.0e-2, n_drydep)  # small default
    else
        heff
    end

    for p in bounds_patch
        mask_patch[p] || continue

        g = patch_gridcell[p]
        c = patch_column[p]

        # Get environmental conditions
        t_sfc = forc_t_downscaled_col[c]
        solar = forc_solar_col[c]
        has_snow = frac_sno[c] > 0.5
        lat = lat_grc[g]
        lai = max(elai_patch[p], 0.0)

        # Determine Wesely season
        season = _wesely_season(lat, month)

        # Map PFT to Wesely land type
        lt = _pft_to_wesely(patch_itype[p])

        # Aerodynamic resistance
        ra = max(ram1_patch[p], 1.0)

        # Friction velocity for boundary layer resistance
        ustar = max(fv_patch[p], 0.001)

        for i in 1:n_drydep
            # Quasi-laminar boundary layer resistance
            rb = _calc_rb(ustar, dd.dv[i])

            # Surface resistance
            rc = _calc_rc(season, lt, lai,
                          dd.foxd[i], heff_vals[i],
                          t_sfc, solar, has_snow)

            # Total deposition velocity [m/s -> cm/s]
            # V_d = 1 / (Ra + Rb + Rc)
            vd = 1.0 / (ra + rb + rc)

            # Convert m/s to cm/s
            dd.velocity_patch[p, i] = vd * 100.0
        end
    end

    return nothing
end

# =========================================================================
# Weighted average to gridcell (for coupling to atmosphere)
# =========================================================================

"""
    drydep_p2g!(dd::DryDepVelocityData,
                ddvel_grc::Matrix{Float64},
                bounds_patch::UnitRange{Int},
                mask_patch::BitVector,
                patch_gridcell::Vector{Int},
                wtgcell::Vector{Float64},
                ng::Int)

Average patch-level dry deposition velocities to gridcell level using
the gridcell weights. Results stored in `ddvel_grc` (ng x n_drydep) [cm/s].

Ported from `DryDepVelocity%p2g` in `DryDepVelocity.F90`.
"""
function drydep_p2g!(dd::DryDepVelocityData,
                     ddvel_grc::Matrix{Float64},
                     bounds_patch::UnitRange{Int},
                     mask_patch::BitVector,
                     patch_gridcell::Vector{Int},
                     wtgcell::Vector{Float64},
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
