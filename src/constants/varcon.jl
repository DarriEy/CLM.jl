# ==========================================================================
# Ported from: src/main/clm_varcon.F90 (332 lines)
# Physical constants and model parameters
# ==========================================================================

# --- Fundamental physical constants ---
# Source: CESM shared constants (shr_const_mod)

const GRAV   = 9.80616         # gravitational acceleration [m/s^2]
const SB     = 5.67e-8         # Stefan-Boltzmann constant [W/m^2/K^4]
const VKC    = 0.4             # von Karman constant [-]
const RWAT   = 461.505         # gas constant for water vapor [J/K/kg]
const RAIR   = 287.04          # gas constant for dry air [J/K/kg]
const ROVERG = RWAT / GRAV * 1000.0  # Rv/g [K/m] * 1000 for [K/km]
const CPLIQ  = 4188.0          # specific heat of liquid water [J/kg/K]
const CPICE  = 2117.27         # specific heat of ice [J/kg/K]
const CPAIR  = 1004.64         # specific heat of dry air [J/kg/K]
const HVAP   = 2.501e6         # latent heat of vaporization [J/kg]
const HSUB   = 2.501e6 + 3.337e5  # latent heat of sublimation [J/kg]
const HFUS   = 3.337e5         # latent heat of fusion [J/kg]
const DENH2O = 1000.0          # density of fresh water [kg/m^3]
const DENICE = 917.0           # density of ice [kg/m^3]
const RGAS   = 8.31446         # universal gas constant [J/K/mol]
const PSTD   = 101325.0        # standard pressure [Pa]

# --- Temperature constants ---
const TFRZ     = 273.15        # freezing point of water [K]
const TCRIT    = 2.5           # critical temperature offset for snow/rain [K]

# --- Mathematical constants ---
const RPI      = 3.14159265358979323846  # pi

# --- Time constants ---
const SECSPHR   = 3600.0       # seconds per hour
const ISECSPHR  = 3600         # integer seconds per hour
const ISECSPMIN = 60           # integer seconds per minute
const SECSPDAY  = 86400.0      # seconds per day
const ISECSPDAY = 86400        # integer seconds per day
const DEGPSEC   = 15.0 / 3600.0  # degrees per second (Earth rotation)

# --- Molecular weights ---
const MWDAIR = 28.966          # molecular weight of dry air [g/mol]
const MWWV   = 18.016          # molecular weight of water vapor [g/mol]
const WV_TO_DAIR_WEIGHT_RATIO = MWWV / MWDAIR

# --- Thermal conductivities ---
const TKAIR = 0.023            # thermal conductivity of air [W/m/K]
const TKICE = 2.290            # thermal conductivity of ice [W/m/K]
const TKWAT = 0.57             # thermal conductivity of water [W/m/K]

# --- Other physical constants ---
const O2_MOLAR_CONST = 0.209   # O2 mole fraction [-]
const ONEATM = 1.01325e5       # one atmosphere [Pa]
const BDSNO  = 250.0           # bulk density of new snow [kg/m^3]
const ALPHA_AERO = 1.0         # aerodynamic factor [-]
const TLSAI_CRIT = 2.0         # critical LAI+SAI for roughness [-]
const WATMIN = 0.01            # minimum soil water [mm]
const RE     = 6.371e3         # Earth radius [km]

# --- Soil/snow parameters ---
const CAPR  = 0.34             # soil thermal capacity tuning parameter
const CNFAC = 0.5              # Crank-Nicholson factor
const PONDMX = 0.0            # ponding depth for soil [mm]
const PONDMX_URBAN = 1.0      # ponding depth for urban [mm]
const THK_BEDROCK  = 3.0      # thermal conductivity of bedrock [W/m/K]
const CSOL_BEDROCK = 2.0e6    # heat capacity of bedrock [J/m^3/K]

# --- Biomass parameters ---
const C_WATER = CPLIQ          # heat capacity of water (alias) [J/kg/K]
const C_DRY_BIOMASS = 1400.0   # heat capacity of dry biomass [J/kg/K]

# --- Snow parameters ---
const H2OSNO_MAX_DEFAULT = -999.0  # default (unlimited) max snow water
const BETADS = 0.5             # snow cover depletion curve parameter (dry)
const BETAIS = 0.5             # snow cover depletion curve parameter (icy)

# --- Sentinel / tolerance values ---
const SMALLVALUE    = 1.0e-12  # small value for floating point comparisons
const SUM_TO_1_TOL  = 1.0e-13 # tolerance for weights summing to 1
const SPVAL         = 1.0e36  # special (missing) value for reals
const ISPVAL        = -9999   # special (missing) value for integers
const FUN_PERIOD    = 1       # FUN update period

# --- Bedrock ---
const ZMIN_BEDROCK = 0.4       # minimum depth to bedrock [m]
const AQUIFER_WATER_BASELINE = 5000.0  # baseline aquifer water [mm]

# --- Error function approximation ---
# Abramowitz & Stegun 7.1.26, ~1.5e-7 accuracy
# Avoids SpecialFunctions.jl dependency
function erf(x::Real)
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

# --- Carbon ---
const C_TO_B = 2.0             # carbon to biomass ratio
const CATOMW = 12.011          # atomic weight of carbon [g/mol]
const G_TO_MG = 1.0e3
const CM3_TO_M3 = 1.0e-6
const PCT_TO_FRAC = 1.0e-2

# --- Carbon isotope constants ---
const PDB = 0.0112372          # Pee Dee Belemnite standard ratio
const PREIND_ATM_DEL13C = -6.0 # pre-industrial atmospheric delta-13C [permil]
const PREIND_ATM_RATIO = PDB + (PREIND_ATM_DEL13C * PDB) / 1000.0
const C3_DEL13C = -28.0        # C3 plant delta-13C [permil]
const C4_DEL13C = -13.0        # C4 plant delta-13C [permil]
const C3_R1 = PDB * (1.0 + C3_DEL13C / 1000.0) / (1.0 + PDB * (1.0 + C3_DEL13C / 1000.0))
const C3_R2 = 1.0 - C3_R1
const C4_R1 = PDB * (1.0 + C4_DEL13C / 1000.0) / (1.0 + PDB * (1.0 + C4_DEL13C / 1000.0))
const C4_R2 = 1.0 - C4_R1
const C14RATIO = 1.0e-12       # baseline C14/C ratio

# --- Nitrogen constants ---
const NITRIF_N2O_LOSS_FRAC = 6.0e-4      # fraction of nitrification lost as N2O
const FRAC_MINRLZTN_TO_NO3 = 0.2         # fraction of mineralization to NO3

# --- Turbulence parameters ---
const BETA_PARAM  = 7.2
const NU_PARAM    = 1.5e-5
const B1_PARAM    = 1.4
const B4_PARAM    = -0.31
const CD1_PARAM   = 7.5
const MEIER_PARAM1 = 0.23
const MEIER_PARAM2 = 0.08
const MEIER_PARAM3 = 70.0

# --- Urban building parameters ---
const HT_WASTEHEAT_FACTOR  = 0.2
const AC_WASTEHEAT_FACTOR  = 0.6
const EM_ROOF_INT  = 0.9
const EM_SUNW_INT  = 0.9
const EM_SHDW_INT  = 0.9
const EM_FLOOR_INT = 0.9
const HCV_ROOF     = 0.948
const HCV_ROOF_ENHANCED = 4.040
const HCV_FLOOR    = 0.948
const HCV_FLOOR_ENHANCED = 4.040
const HCV_SUNW     = 3.076
const HCV_SHDW     = 3.076
const DZ_FLOOR     = 0.1
const VENT_ACH     = 0.3
const WASTEHEAT_LIMIT = 100.0
const DENS_FLOOR   = 2.35e3   # floor density [kg/m^3]
const SH_FLOOR     = 880.0    # floor specific heat [J/kg/K]
const CP_FLOOR     = DENS_FLOOR * SH_FLOOR  # floor volumetric heat capacity

# --- Temperature reference ---
const KH_TBASE = 298.0         # base temperature for Henry's law [K]

# --- Omega (single-scattering albedo) for snow ---
# omegas[IVIS], omegas[INIR] — set at init
const OMEGAS = [0.8, 0.4]

# --- Grid/subgrid name constants ---
const GRLND      = "lndgrid"
const NAMEG      = "gridcell"
const NAMEL      = "landunit"
const NAMEC      = "column"
const NAMEP      = "pft"
const NAMECOHORT = "cohort"

# --- Gas constants for methane/O2/CO2 ---
# s_con: solubility function coefficients [ngases × 4]
# Gases: 1=CH4, 2=O2, 3=CO2
const S_CON = [
    1.3e-3    1700.0   5.7     0.0     # CH4
    2.1e-3    1500.0   6.7     0.0     # O2
    8.6e-6    2400.0   6.5     0.0     # CO2
]

# d_con_w: diffusivity in water [ngases × 3]
const D_CON_W = [
    1.5e-9   1.7e-9   1.0e3   # CH4
    2.4e-9   2.0e-9   1.0e3   # O2
    1.9e-9   1.8e-9   1.0e3   # CO2
]

# d_con_g: diffusivity in air [ngases × 2]
const D_CON_G = [
    2.0e-5   1.0e3    # CH4
    2.0e-5   1.0e3    # O2
    1.5e-5   1.0e3    # CO2
]

# Henry's law constants [ngases]
const C_H_INV = [600.0, 1.3, 36.0]
const KH_THETA = [1600.0, 1500.0, 2400.0]

# Snow capping constants
const H2OSNO_MAX = 10000.0       # Maximum allowed SWE [mm H2O]
const MIN_SNOW_TO_KEEP = 1.0e-3  # Minimum fraction of snow mass to retain during capping

# --- Runtime-allocated vertical coordinate arrays ---
# These get populated during initialization
const zlak  = Ref(Float64[])    # lake level depths [m]
const dzlak = Ref(Float64[])    # lake level thicknesses [m]
const zsoi  = Ref(Float64[])    # soil level depths [m]
const dzsoi = Ref(Float64[])    # soil level thicknesses [m]
const zisoi = Ref(Float64[])    # soil interface depths [m]
const dzsoi_decomp = Ref(Float64[])  # decomposition level thicknesses [m]

"""
    varcon_init!()

Initialize runtime-dependent constants (vertical coordinate arrays, etc.).
Called after varpar_init!() has set level counts.
"""
function varcon_init!()
    nlevsoi = varpar.nlevsoi
    nlevgrnd = varpar.nlevgrnd
    nlevlak = varpar.nlevlak

    # CLM5 "20SL_8.5m" soil layer structure
    # Piecewise linear layer thicknesses from initVerticalMod.F90
    dzsoi_val = zeros(nlevgrnd)
    zsoi_val = zeros(nlevgrnd)
    zisoi_val = zeros(nlevgrnd + 1)

    # Layer thicknesses: linearly increasing within 3 tiers
    for j in 1:min(4, nlevsoi)
        dzsoi_val[j] = j * 0.02                              # 0.02, 0.04, 0.06, 0.08 m
    end
    for j in 5:min(13, nlevsoi)
        dzsoi_val[j] = dzsoi_val[4] + (j - 4) * 0.04        # 0.12 → 0.44 m
    end
    for j in 14:nlevsoi
        dzsoi_val[j] = dzsoi_val[min(13, nlevsoi)] + (j - 13) * 0.10  # 0.54 → 1.14 m
    end
    # Bedrock layers: exponential thickening
    for j in (nlevsoi+1):nlevgrnd
        dzsoi_val[j] = dzsoi_val[nlevsoi] + (((j - nlevsoi) * 25.0)^1.5) / 100.0
    end

    # Interface depths (cumulative sum of layer thicknesses)
    zisoi_val[1] = 0.0
    for j in 1:nlevgrnd
        zisoi_val[j+1] = zisoi_val[j] + dzsoi_val[j]
    end

    # Node depths (layer midpoints)
    for j in 1:nlevgrnd
        zsoi_val[j] = 0.5 * (zisoi_val[j] + zisoi_val[j+1])
    end

    zsoi[]  = zsoi_val
    dzsoi[] = dzsoi_val
    zisoi[] = zisoi_val

    # Lake levels
    if varctl.use_extralakelayers
        dzlak_val = fill(1.0, nlevlak)  # 25 × 1m layers
    else
        dzlak_val = [1.0, 2.0, 3.0, 4.0, 5.0, 7.0, 7.0, 7.0, 7.0, 7.0]
    end

    zlak_val = zeros(nlevlak)
    zlak_val[1] = 0.5 * dzlak_val[1]
    for j in 2:nlevlak
        zlak_val[j] = zlak_val[j-1] + 0.5 * (dzlak_val[j-1] + dzlak_val[j])
    end

    zlak[]  = zlak_val
    dzlak[] = dzlak_val

    # Decomposition level thicknesses (use soil level thicknesses)
    dzsoi_decomp[] = dzsoi_val[1:min(varpar.nlevdecomp_full, nlevsoi)]

    nothing
end
