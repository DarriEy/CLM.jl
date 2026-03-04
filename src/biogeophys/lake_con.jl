# ==========================================================================
# Ported from: src/biogeophys/LakeCon.F90
# Lake constants and parameters for CLM lake code (CLM4-LISSS)
# Documented in Subin et al. 2011, JAMES
#
# Contains non-tuneable and tuneable constants, namelist-controlled
# parameters, and the initialization routine for lake physics.
#
# Public constants:
#   tdmax, emg_lake, z0frzlake, za_lake, cur0, cus, curm, n2min
#
# Public mutable parameters (LakeConParams):
#   betavis, lsadz, fcrit, minz0lake, pudz, depthcrit, mixfact,
#   lake_use_old_fcrit_minz0, deepmixing_depthcrit, deepmixing_mixfact,
#   lake_no_ed, lakepuddling, lake_puddle_thick
#
# Public functions:
#   lake_con_init!  — Initialize time-invariant lake variables
# ==========================================================================

# ------------------------------------------------------------------
# Lake Model non-tuneable constants
# ------------------------------------------------------------------

# Temperature of maximum water density (K)
# From Hostetler and Bartlein (1990); more updated sources suggest 277.13 K.
const tdmax = 277.0

# ------------------------------------------------------------------
# Lake Model tuneable constants (Fortran parameters)
# ------------------------------------------------------------------

# Lake emissivity (used for both frozen and unfrozen lakes)
const emg_lake = 0.97

# Momentum roughness length over frozen lakes without snow (m)
# Typical value from literature, consistent with Mironov expressions.
# See e.g. Morris 1989, Andreas 1987, Guest & Davidson 1991 (as cited in Vavrus 1996)
const z0frzlake = 0.001

# Base of surface light absorption layer for lakes (m)
const za_lake = 0.6

# Charnock parameter constants for prognostic roughness length
const cur0 = 0.01    # minimum Charnock parameter
const cus  = 0.1     # empirical constant for roughness under smooth flow
const curm = 0.1     # maximum Charnock parameter

# Minimum Brunt-Vaisala frequency squared for enhanced diffusivity (s^-2)
# Yields diffusivity about 6 times molecular; from Fang & Stefan 1996
const n2min = 7.5e-5

# ------------------------------------------------------------------
# Mutable lake parameters (namelist-controlled and derived)
# ------------------------------------------------------------------

"""
    LakeConParams

Mutable container for lake model parameters that are set at runtime
via namelists or derived in `lake_con_init!`.

Fields correspond to Fortran module-level variables in `LakeCon.F90`.
"""
Base.@kwdef mutable struct LakeConParams
    # --- Namelist-controlled parameters (with defaults) ---

    # Fraction of visible sunlight absorbed in ~1 m of water (the surface layer za_lake).
    # As long as NIR = 700 nm and up, this can be zero.
    betavis::Float64 = 0.0

    # Additional snow layer thickness for lake numerics (m).
    # Adjusted in lake_con_init! if timestep is not 1800 s.
    # See LakeCon.F90 comments for CFL analysis of minimum snow layer thickness.
    lsadz::Float64 = 0.03

    # true => use old fcrit & minz0 as per Subin et al 2011 form.
    # See lake_con_init! for details. Currently hardwired off.
    lake_use_old_fcrit_minz0::Bool = false

    # (m) Minimum lake depth to invoke deep mixing
    deepmixing_depthcrit::Float64 = 25.0

    # Factor to increase mixing by for deep lakes
    deepmixing_mixfact::Float64 = 10.0

    # true => suppress enhanced diffusion. Small differences.
    # Currently hardwired false. See Subin et al 2011 for details.
    lake_no_ed::Bool = false

    # true => suppress convection when greater than minimum amount of ice present.
    # Also effectively sets lake_no_melt_icealb. Not extensively tested.
    lakepuddling::Bool = false

    # (m) Minimum total ice nominal thickness before convection is suppressed
    lake_puddle_thick::Float64 = 0.2

    # --- Derived parameters (set by lake_con_init!) ---

    # Critical dimensionless fetch for Charnock parameter
    fcrit::Float64 = NaN

    # (m) Minimum allowed roughness length for unfrozen lakes
    minz0lake::Float64 = NaN

    # (m) Minimum total ice thickness required to allow lake puddling
    pudz::Float64 = NaN

    # (m) Depth beneath which to increase mixing. See Subin et al. 2011
    depthcrit::Float64 = NaN

    # Mixing increase factor for deep lakes
    mixfact::Float64 = NaN
end

# Global instance of lake parameters
const lake_con = LakeConParams()

"""
    lake_con_init!(params::LakeConParams)

Initialize time-invariant variables for the lake model.

Sets `fcrit`, `minz0lake`, `pudz`, `depthcrit`, and `mixfact` based on
namelist-controlled fields in `params`.

Ported from `LakeConInit` in `LakeCon.F90`.
"""
function lake_con_init!(params::LakeConParams = lake_con)

    # Set fcrit and minz0lake based on lake_use_old_fcrit_minz0
    if params.lake_use_old_fcrit_minz0
        # Critical dimensionless fetch for Charnock parameter.
        # From Vickers & Mahrt 1997 but converted to use u instead of u*
        # (Form used in Subin et al. 2011)
        params.fcrit = 22.0

        # (m) Minimum allowed roughness length for unfrozen lakes.
        # (Used in Subin et al. 2011)
        params.minz0lake = 1.0e-5
    else
        # Vickers & Mahrt 1997
        params.fcrit = 100.0

        # (m) Minimum allowed roughness length for unfrozen lakes.
        # Now set low so it is only to avoid floating point exceptions.
        params.minz0lake = 1.0e-10
    end

    # Set pudz if lake puddling is enabled
    if params.lakepuddling
        # (m) Minimum total ice thickness required to allow lake puddling.
        # Default is 0.2m. This option has not been extensively tested.
        # This option turns on lake_no_melt_icealb, as the decrease in
        # albedo will be based on whether there is water over ice, not
        # purely a function of ice top temperature.
        params.pudz = params.lake_puddle_thick
    else
        params.pudz = 0.0
    end

    # (m) Depth beneath which to increase mixing. See Subin et al. 2011
    params.depthcrit = params.deepmixing_depthcrit

    # Mixing increase factor. Defaults are 25 m depth, increase by 10.
    params.mixfact = params.deepmixing_mixfact

    return nothing
end
