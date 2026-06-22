# FatesConstantsMod.jl
# Julia port of FATES src/fates/main/FatesConstantsMod.F90
#
# This module is used to define global _immutable_ data. Everything in this
# module is a `const` (Fortran `parameter`). TRUE LEAF — no FATES dependencies.
# `fates_r8` maps to Float64.

# ---------------------------------------------------------------------------
# kinds
# ---------------------------------------------------------------------------
const fates_r8 = Float64    # 8 byte real (selected_real_kind(12))
const fates_int = Int32     # 4 byte int  (selected_int_kind(8))

# ---------------------------------------------------------------------------
# string lengths
# ---------------------------------------------------------------------------
const fates_avg_flag_length = 3
const fates_short_string_length = 32
const fates_long_string_length = 199

# ---------------------------------------------------------------------------
# Unset sentinels
# ---------------------------------------------------------------------------
# Used to initialize and test unset integers
const fates_unset_int = -9999

# Used to initialize and test unset r8s
const fates_unset_r8 = -1.0e36

# Used to check if a parameter was specified in the parameter file (or left as _)
const fates_check_param_set = 9.9e32

# Integer equivalent of true / false (in case some compilers dont auto convert)
const itrue = 1
const ifalse = 0

# ---------------------------------------------------------------------------
# Enumerations
# ---------------------------------------------------------------------------
# dbh bins used when comparing patches
const N_DBH_BINS = 6
const patchfusion_dbhbin_loweredges = [0.0, 5.0, 20.0, 50.0, 100.0, 150.0]

# Disturbance modes
const N_DIST_TYPES = 4          # 1) tree-fall, 2) fire, 3) logging, 4) land-use change
const dtype_ifall = 1           # naturally occurring tree-fall generated event
const dtype_ifire = 2           # fire generated disturbance event
const dtype_ilog = 3            # logging generated disturbance event
const dtype_ilandusechange = 4  # land use change disturbance (not including logging)

# Labels for patch land use type information
const n_landuse_cats = 5
const primaryland = 1
const secondaryland = 2
const rangeland = 3
const pastureland = 4
const cropland = 5
const is_crop = [false, false, false, false, true]
const n_crop_lu_types = 1

# Bareground nocomp land use label (not a real land use type, only for labeling)
const nocomp_bareground_land = 0

# Bareground nocomp PFT label for no competition mode
const nocomp_bareground = 0

# Deciduous leaf status flags
const leaves_on = 2        # deciduous plant has leaves, allocating to them
const leaves_off = 1       # deciduous plant has dropped its leaves
const leaves_shedding = 3  # deciduous plant has leaves but is shedding (partial)

# Stress (drought) deciduous flavors
const ihard_stress_decid = 1  # "hard" deciduous (two statuses)
const isemi_stress_decid = 2  # semi-deciduous (can downregulate)

# Canopy layer indices
const ican_upper = 1   # nominal index for the upper canopy
const ican_ustory = 2  # diagnostics referring to understory layers

# Phosphorus uptake interaction flags
const prescribed_p_uptake = 1
const coupled_p_uptake = 2

# Nitrogen uptake interaction flags
const prescribed_n_uptake = 1
const coupled_n_uptake = 2
const coupled_np_comp_scaling = 1  # at least 1 chemical element (N or P)

# Tree regeneration scheme flags
const TRS_no_seedling_dyn = 3   # size-based reproductive allocation, no seedling dynamics
const TRS_regeneration = 2      # full TRS
const default_regeneration = 1  # FATES's default regeneration scheme
const min_max_dbh_for_trees = 15.0  # below this max dbh, use default regen (shrubs/grasses)

const trivial_np_comp_scaling = 2  # nutrients off or plants not coupled with belowground chem

# Scaling of nutrient competitors presented to the HLM's soil BGC model.
# This one is mutable in Fortran (no `parameter`); use a Ref so it can be set.
const fates_np_comp_scaling = Ref{Int}(fates_unset_int)

const secondary_age_threshold = 94.0  # young secondary land below this avg age (Hurtt-2011)

# Photosynthesis acclimation models
const photosynth_acclim_model_none = 1
const photosynth_acclim_model_kumarathunge_etal_2019 = 2

# Harvest unit labels
const hlm_harvest_area_fraction = 1  # harvesting by area
const hlm_harvest_carbon = 2         # harvesting based on carbon extracted

# Leaf maintenance respiration models
const lmrmodel_ryan_1991 = 1
const lmrmodel_atkin_etal_2017 = 2

# Carbon starvation models
const cstarvation_model_lin = 1  # Linear scaling
const cstarvation_model_exp = 2  # Exponential scaling

# ---------------------------------------------------------------------------
# Error Tolerances
# ---------------------------------------------------------------------------
const calloc_abs_error = 1.0e-9   # carbon allocation error tol (kgC/plant ~ 1 microgram)
const area_error_1 = 1.0e-16      # area checks (canopy, patch)
const area_error_2 = 1.0e-12      # tree lai checks
const area_error_3 = 10.0e-9      # area checks (canopy, patch)
const area_error_4 = 1.0e-10      # area checks
const rsnbl_math_prec = 1.0e-12   # multiplication/division rounding-error expectation
const min_nocomp_pftfrac_perlanduse = 0.01  # min PFT fraction for any land use type (nocomp)

const tinyr8 = floatmin(Float64)  # precision of 8byte reals (~1e-308)
const nearzero = 1.0e-30          # near-zero comparison threshold

# ---------------------------------------------------------------------------
# Unit conversion constants
# ---------------------------------------------------------------------------
const umolC_to_kgC = 12.0e-9   # umols C -> kg C (1 mol = 12g)
const mg_per_kg = 1.0e6
const g_per_kg = 1000.0
const kg_per_g = 0.001
const mg_per_g = 1000.0
const kg_per_Megag = 1000.0
const umol_per_mmol = 1000.0
const mmol_per_mol = 1000.0
const umol_per_mol = 1.0e6
const mol_per_umol = 1.0e-6
const umol_per_kmol = 1.0e9
const m_per_mm = 1.0e-3
const mm_per_m = 1.0e3
const mm_per_cm = 10.0
const m_per_cm = 1.0e-2
const m2_per_ha = 1.0e4
const m2_per_km2 = 1.0e6
const cm2_per_m2 = 10000.0
const m3_per_mm3 = 1.0e-9
const m3_per_cm3 = 1.0e-6
const cm3_per_m3 = 1.0e6
const ha_per_m2 = 1.0e-4
const sec_per_min = 60.0
const sec_per_day = 86400.0
const megajoules_per_joule = 1.0e-6
const days_per_sec = 1.0 / 86400.0
const days_per_year = 365.00  # assume HLM uses 365 day calendar
const ndays_per_year = round(Int, days_per_year)
const years_per_day = 1.0 / 365.00
const months_per_year = 12.0
const J_per_kJ = 1000.0

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
const dewpoint_a = 17.62
const dewpoint_b = 243.12  # [degrees C]

const rgas_J_K_kmol = 8314.4598     # universal gas constant [J/K/kmol]
const rgas_J_K_mol = 8.3144598      # universal gas constant [J/K/mol]
const t_water_freeze_k_1atm = 273.15    # freezing point of water at 1 atm (K)
const t_water_freeze_k_triple = 273.16  # freezing point at triple point (K)
const dens_fresh_liquid_water = 1.0e3   # density of fresh liquid water (kg/m3)
const molar_mass_water = 18.0           # molar mass of water (g/mol)
const molar_mass_ratio_vapdry = 0.622   # molar mass ratio water vapor to dry air (-)
const grav_earth = 9.8                  # gravity on earth [m/s2]
const pa_per_mpa = 1.0e6                # megapascals to pascals
const mpa_per_pa = 1.0e-6               # pascals to megapascals
const mpa_per_mm_suction = dens_fresh_liquid_water * grav_earth * 1.0e-9

const fates_huge = floatmax(Float64)  # huge(g_per_kg)
const fates_tiny = floatmin(Float64)  # tiny(g_per_kg)

# Geodesy constants (WGS 84)
const earth_radius_eq = 6378137.0           # equatorial radius, earth [m]
const earth_flattening = 1.0 / 298.257223563  # flattening [non-dimensional]

# ---------------------------------------------------------------------------
# Geometric Constants
# ---------------------------------------------------------------------------
const pi_const = 3.14159265359
const rad_per_deg = pi_const / 180.0

# Rdark constants from Atkin et al., 2017 and Heskel et al., 2016
const lmr_b = 0.1012      # (degrees C^-1)
const lmr_c = -0.0005     # (degrees C^-2)
const lmr_TrefC = 25.0    # (degrees C)
const lmr_r_1 = 0.2061    # (umol CO2/m^2/s / (gN/(m2 leaf)))
const lmr_r_2 = -0.0402   # (umol CO2/m^2/s/degree C)

# Termination mortality types
const n_term_mort_types = 3
const i_term_mort_type_cstarv = 1
const i_term_mort_type_canlev = 2
const i_term_mort_type_numdens = 3
