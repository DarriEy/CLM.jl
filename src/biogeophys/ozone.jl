# ==========================================================================
# Ported from: src/biogeophys/OzoneMod.F90
# Calculates ozone-induced stress.
#
# Note that the ozone calculations need to happen AFTER rssun and rssha are
# computed by the Photosynthesis routine. However, Photosynthesis also uses
# the ozone stress computed here. Thus, the ozone stress computed in timestep
# i is applied in timestep (i+1), requiring these stresses to be saved on
# the restart file.
#
# Developed by Danica Lombardozzi: Lombardozzi, D., S. Levis, G. Bonan,
# P. G. Hess, and J. P. Sparks (2015), The Influence of Chronic Ozone
# Exposure on Global Carbon and Water Cycles, J Climate, 28(1), 292–305.
#
# Public functions:
#   calc_ozone_uptake!              — Calculate ozone uptake for all patches
#   calc_ozone_stress!              — Calculate ozone stress (dispatch)
#
# Private functions:
#   calc_ozone_uptake_one_point!    — Ozone uptake for a single point
#   calc_ozone_stress_lombardozzi2015!          — Lombardozzi2015 stress
#   calc_ozone_stress_lombardozzi2015_one_point — Lombardozzi2015 one point
#   calc_ozone_stress_falk!                     — Falk stress
#   calc_ozone_stress_falk_one_point            — Falk one point
# ==========================================================================

# --- Stress method constants ---
const STRESS_METHOD_LOMBARDOZZI2015 = 1
const STRESS_METHOD_FALK = 2

# --- Ozone parameters (module-level constants from Fortran) ---

# o3:h2o resistance ratio defined by Sitch et al. 2007
const O3_KO3 = 1.67

# LAI threshold for LAIs that asymptote and don't reach 0
const O3_LAI_THRESH = 0.5

# threshold below which o3flux is set to 0 (nmol m^-2 s^-1)
const O3_FLUX_THRESHOLD = 0.8

# o3 intercepts and slopes for photosynthesis
const O3_NEEDLELEAF_PHOTO_INT   = 0.8390   # unitless
const O3_NEEDLELEAF_PHOTO_SLOPE = 0.0      # per mmol m^-2
const O3_BROADLEAF_PHOTO_INT    = 0.8752   # unitless
const O3_BROADLEAF_PHOTO_SLOPE  = 0.0      # per mmol m^-2
const O3_NONWOODY_PHOTO_INT     = 0.8021   # unitless
const O3_NONWOODY_PHOTO_SLOPE   = -0.0009  # per mmol m^-2

# o3 intercepts and slopes for conductance
const O3_NEEDLELEAF_COND_INT    = 0.7823   # unitless
const O3_NEEDLELEAF_COND_SLOPE  = 0.0048   # per mmol m^-2
const O3_BROADLEAF_COND_INT     = 0.9125   # unitless
const O3_BROADLEAF_COND_SLOPE   = 0.0      # per mmol m^-2
const O3_NONWOODY_COND_INT      = 0.7511   # unitless
const O3_NONWOODY_COND_SLOPE    = 0.0      # per mmol m^-2

# o3 intercepts and slopes for JmaxO3/Jmax0
const O3_NEEDLELEAF_JMAX_INT    = 1.0      # unitless
const O3_NEEDLELEAF_JMAX_SLOPE  = 0.0      # per mmol m^-2
const O3_BROADLEAF_JMAX_INT     = 1.0      # unitless
const O3_BROADLEAF_JMAX_SLOPE   = -0.0037  # per mmol m^-2
const O3_NONWOODY_JMAX_INT      = 1.0      # unitless
const O3_NONWOODY_JMAX_SLOPE    = 0.0      # per mmol m^-2

# ==========================================================================
# Ozone data type
# ==========================================================================

"""
    OzoneData

Ozone state variables. Combines both the base type fields (ozone_base_type)
and the derived type fields (ozone_type) from Fortran.

Fields from ozone_base_type:
  - o3coefvsha_patch, o3coefvsun_patch: ozone coefs for photosynthesis (0-1)
  - o3coefgsha_patch, o3coefgsun_patch: ozone coefs for conductance (0-1)
  - o3coefjmaxsha_patch, o3coefjmaxsun_patch: ozone coefs for Jmax (0-1)

Fields from ozone_type:
  - stress_method: which stress parameterization to use
  - o3uptakesha_patch, o3uptakesun_patch: ozone dose (mmol O3/m^2)
  - tlai_old_patch: tlai from last time step
"""
Base.@kwdef mutable struct OzoneData{FT<:Real}
    # Stress method selector
    stress_method::Int = STRESS_METHOD_LOMBARDOZZI2015

    # --- Base type fields (from ozone_base_type) ---
    o3coefvsha_patch    ::Vector{FT} = Float64[]   # ozone coef for photosynthesis, shaded (0-1)
    o3coefvsun_patch    ::Vector{FT} = Float64[]   # ozone coef for photosynthesis, sunlit (0-1)
    o3coefgsha_patch    ::Vector{FT} = Float64[]   # ozone coef for conductance, shaded (0-1)
    o3coefgsun_patch    ::Vector{FT} = Float64[]   # ozone coef for conductance, sunlit (0-1)
    o3coefjmaxsha_patch ::Vector{FT} = Float64[]   # ozone coef for Jmax, shaded (0-1)
    o3coefjmaxsun_patch ::Vector{FT} = Float64[]   # ozone coef for Jmax, sunlit (0-1)

    # --- Derived type fields (from ozone_type) ---
    o3uptakesha_patch   ::Vector{FT} = Float64[]   # ozone dose, shaded leaves (mmol O3/m^2)
    o3uptakesun_patch   ::Vector{FT} = Float64[]   # ozone dose, sunlit leaves (mmol O3/m^2)
    tlai_old_patch      ::Vector{FT} = Float64[]   # tlai from last time step
end

# ==========================================================================
# Init / clean functions
# ==========================================================================

"""
    ozone_init!(oz, npatches; stress_method, use_luna)

Allocate and initialize ozone data for `npatches` patches.

Ported from `Init`, `InitAllocate`, `InitHistory`, and `InitCold` in OzoneMod.F90.
"""
function ozone_init!(oz::OzoneData{FT}, npatches::Int;
                     stress_method::String = "stress_lombardozzi2015",
                     use_luna::Bool = false) where {FT}
    # Set stress method
    if stress_method == "stress_lombardozzi2015"
        oz.stress_method = STRESS_METHOD_LOMBARDOZZI2015
    elseif stress_method == "stress_falk"
        oz.stress_method = STRESS_METHOD_FALK
        if !use_luna
            error("use_luna=true is required when stress_method = stress_falk")
        end
    else
        error("ozone_init!: unknown ozone stress method: $stress_method")
    end

    # InitAllocateBase: base type fields initialized to NaN, then cold-start to 1.0
    oz.o3coefvsha_patch    = fill(one(FT), npatches)
    oz.o3coefvsun_patch    = fill(one(FT), npatches)
    oz.o3coefgsha_patch    = fill(one(FT), npatches)
    oz.o3coefgsun_patch    = fill(one(FT), npatches)
    oz.o3coefjmaxsha_patch = fill(one(FT), npatches)
    oz.o3coefjmaxsun_patch = fill(one(FT), npatches)

    # InitAllocate + InitCold: derived type fields
    oz.o3uptakesha_patch   = fill(zero(FT), npatches)
    oz.o3uptakesun_patch   = fill(zero(FT), npatches)
    oz.tlai_old_patch      = fill(zero(FT), npatches)

    return nothing
end

"""
    ozone_clean!(oz)

Deallocate (reset to empty) all fields of an `OzoneData` instance.
"""
function ozone_clean!(oz::OzoneData{FT}) where {FT}
    oz.o3coefvsha_patch    = FT[]
    oz.o3coefvsun_patch    = FT[]
    oz.o3coefgsha_patch    = FT[]
    oz.o3coefgsun_patch    = FT[]
    oz.o3coefjmaxsha_patch = FT[]
    oz.o3coefjmaxsun_patch = FT[]
    oz.o3uptakesha_patch   = FT[]
    oz.o3uptakesun_patch   = FT[]
    oz.tlai_old_patch      = FT[]
    return nothing
end

# ==========================================================================
# Ozone uptake
# ==========================================================================

"""
    calc_ozone_uptake!(oz, patchdata, mask_exposedvegp, bounds,
                       forc_pbot, forc_th, rssun, rssha, rb, ram, tlai, forc_o3,
                       evergreen_pft, leaf_long_pft, dtime)

Calculate ozone uptake for all exposed vegetation patches.

Ported from `CalcOzoneUptake` in OzoneMod.F90.

# Arguments
- `oz::OzoneData`: ozone state (modified in place)
- `patchdata::PatchData`: patch hierarchy data
- `mask_exposedvegp::BitVector`: mask for non-snow-covered vegetation patches
- `bounds::UnitRange{Int}`: patch index range
- `forc_pbot::Vector{<:Real}`: atmospheric pressure (Pa), column-indexed
- `forc_th::Vector{<:Real}`: atmospheric potential temperature (K), column-indexed
- `rssun::Vector{<:Real}`: leaf stomatal resistance, sunlit leaves (s/m)
- `rssha::Vector{<:Real}`: leaf stomatal resistance, shaded leaves (s/m)
- `rb::Vector{<:Real}`: boundary layer resistance (s/m)
- `ram::Vector{<:Real}`: aerodynamical resistance (s/m)
- `tlai::Vector{<:Real}`: one-sided leaf area index, no burying by snow
- `forc_o3::Vector{<:Real}`: ozone partial pressure (mol/mol), gridcell-indexed
- `evergreen_pft::Vector{<:Real}`: evergreen flag per PFT (1=evergreen)
- `leaf_long_pft::Vector{<:Real}`: leaf longevity per PFT (years)
- `dtime::Float64`: time step size (seconds)
"""
function calc_ozone_uptake!(oz::OzoneData,
                            patchdata::PatchData,
                            mask_exposedvegp::BitVector,
                            bounds::UnitRange{Int},
                            forc_pbot::Vector{<:Real},
                            forc_th::Vector{<:Real},
                            rssun::Vector{<:Real},
                            rssha::Vector{<:Real},
                            rb::Vector{<:Real},
                            ram::Vector{<:Real},
                            tlai::Vector{<:Real},
                            forc_o3::Vector{<:Real},
                            evergreen_pft::Vector{<:Real},
                            leaf_long_pft::Vector{<:Real},
                            dtime::Real)
    for p in bounds
        mask_exposedvegp[p] || continue

        c = patchdata.column[p]
        g = patchdata.gridcell[p]
        pft_type = patchdata.itype[p]

        # Ozone uptake for shaded leaves
        oz.o3uptakesha_patch[p] = calc_ozone_uptake_one_point(
            forc_ozone = forc_o3[g],
            forc_pbot  = forc_pbot[c],
            forc_th    = forc_th[c],
            rs         = rssha[p],
            rb         = rb[p],
            ram        = ram[p],
            tlai_val   = tlai[p],
            tlai_old   = oz.tlai_old_patch[p],
            pft_type   = pft_type,
            o3uptake   = oz.o3uptakesha_patch[p],
            evergreen_pft = evergreen_pft,
            leaf_long_pft = leaf_long_pft,
            dtime      = dtime)

        # Ozone uptake for sunlit leaves
        oz.o3uptakesun_patch[p] = calc_ozone_uptake_one_point(
            forc_ozone = forc_o3[g],
            forc_pbot  = forc_pbot[c],
            forc_th    = forc_th[c],
            rs         = rssun[p],
            rb         = rb[p],
            ram        = ram[p],
            tlai_val   = tlai[p],
            tlai_old   = oz.tlai_old_patch[p],
            pft_type   = pft_type,
            o3uptake   = oz.o3uptakesun_patch[p],
            evergreen_pft = evergreen_pft,
            leaf_long_pft = leaf_long_pft,
            dtime      = dtime)

        oz.tlai_old_patch[p] = tlai[p]
    end
    return nothing
end

"""
    calc_ozone_uptake_one_point(; forc_ozone, forc_pbot, forc_th,
        rs, rb, ram, tlai_val, tlai_old, pft_type, o3uptake,
        evergreen_pft, leaf_long_pft, dtime) -> Float64

Calculate ozone uptake for a single point, for just sunlit or shaded leaves.
Returns the updated o3uptake value.

Ported from `CalcOzoneUptakeOnePoint` in OzoneMod.F90.
"""
function calc_ozone_uptake_one_point(;
        forc_ozone::Real,
        forc_pbot::Real,
        forc_th::Real,
        rs::Real,
        rb::Real,
        ram::Real,
        tlai_val::Real,
        tlai_old::Real,
        pft_type::Int,
        o3uptake::Real,
        evergreen_pft::Vector{<:Real},
        leaf_long_pft::Vector{<:Real},
        dtime::Real)

    # convert o3 from mol/mol to nmol m^-3
    # SHR_CONST_RGAS = 8.31446 (same as RGAS in varcon.jl)
    o3concnmolm3 = forc_ozone * 1.0e9 * (forc_pbot / (forc_th * RGAS))

    # calculate instantaneous flux
    o3flux = o3concnmolm3 / (O3_KO3 * rs + rb + ram)

    # apply o3 flux threshold
    if o3flux < O3_FLUX_THRESHOLD
        o3fluxcrit = 0.0
    else
        o3fluxcrit = o3flux - O3_FLUX_THRESHOLD
    end

    dtimeh = dtime / 3600.0

    # calculate o3 flux per timestep
    o3fluxperdt = o3fluxcrit * dtime * 0.000001

    if tlai_val > O3_LAI_THRESH
        # checking if new leaf area was added
        if tlai_val - tlai_old > 0.0
            # minimizing o3 damage to new leaves
            heal = max(0.0, ((tlai_val - tlai_old) / tlai_val) * o3fluxperdt)
        else
            heal = 0.0
        end

        if evergreen_pft[pft_type] == 1.0
            leafturn = 1.0 / (leaf_long_pft[pft_type] * 365.0 * 24.0)
        else
            leafturn = 0.0
        end

        # o3 uptake decay based on leaf lifetime for evergreen plants
        decay = o3uptake * leafturn * dtimeh
        # cumulative uptake (mmol m^-2)
        o3uptake_new = max(0.0, o3uptake + o3fluxperdt - decay - heal)
    else
        o3uptake_new = 0.0
    end

    return o3uptake_new
end

# ==========================================================================
# Ozone stress (dispatcher)
# ==========================================================================

"""
    calc_ozone_stress!(oz, mask_exposedvegp, mask_noexposedvegp, bounds,
                       patchdata, woody_pft; is_time_to_run_luna)

Calculate ozone stress, dispatching to the appropriate method.

Ported from `CalcOzoneStress` in OzoneMod.F90.
"""
function calc_ozone_stress!(oz::OzoneData,
                            mask_exposedvegp::BitVector,
                            mask_noexposedvegp::BitVector,
                            bounds::UnitRange{Int},
                            patchdata::PatchData,
                            woody_pft::Vector{<:Real};
                            is_time_to_run_luna::Bool = false)
    if oz.stress_method == STRESS_METHOD_LOMBARDOZZI2015
        calc_ozone_stress_lombardozzi2015!(oz, mask_exposedvegp, mask_noexposedvegp,
                                           bounds, patchdata, woody_pft)
    elseif oz.stress_method == STRESS_METHOD_FALK
        calc_ozone_stress_falk!(oz, mask_exposedvegp, mask_noexposedvegp,
                                bounds, patchdata, woody_pft,
                                is_time_to_run_luna = is_time_to_run_luna)
    else
        error("calc_ozone_stress!: unknown ozone stress method: $(oz.stress_method)")
    end
    return nothing
end

# ==========================================================================
# Lombardozzi 2015 stress
# ==========================================================================

"""
    calc_ozone_stress_lombardozzi2015!(oz, mask_exposedvegp, mask_noexposedvegp,
                                       bounds, patchdata, woody_pft)

Calculate ozone stress using the Lombardozzi 2015 formulation.

Ported from `CalcOzoneStressLombardozzi2015` in OzoneMod.F90.
"""
function calc_ozone_stress_lombardozzi2015!(oz::OzoneData,
                                            mask_exposedvegp::BitVector,
                                            mask_noexposedvegp::BitVector,
                                            bounds::UnitRange{Int},
                                            patchdata::PatchData,
                                            woody_pft::Vector{<:Real})
    for p in bounds
        if mask_exposedvegp[p]
            pft_type = patchdata.itype[p]

            # Ozone stress for shaded leaves
            o3coefv_sha, o3coefg_sha = calc_ozone_stress_lombardozzi2015_one_point(
                pft_type = pft_type,
                o3uptake = oz.o3uptakesha_patch[p],
                woody_pft = woody_pft)
            oz.o3coefvsha_patch[p] = o3coefv_sha
            oz.o3coefgsha_patch[p] = o3coefg_sha

            # Ozone stress for sunlit leaves
            o3coefv_sun, o3coefg_sun = calc_ozone_stress_lombardozzi2015_one_point(
                pft_type = pft_type,
                o3uptake = oz.o3uptakesun_patch[p],
                woody_pft = woody_pft)
            oz.o3coefvsun_patch[p] = o3coefv_sun
            oz.o3coefgsun_patch[p] = o3coefg_sun

        elseif mask_noexposedvegp[p]
            # Reset to 1 over non-exposed veg points
            oz.o3coefvsha_patch[p] = 1.0
            oz.o3coefgsha_patch[p] = 1.0
            oz.o3coefvsun_patch[p] = 1.0
            oz.o3coefgsun_patch[p] = 1.0
        end
    end
    return nothing
end

"""
    calc_ozone_stress_lombardozzi2015_one_point(; pft_type, o3uptake, woody_pft)
        -> (o3coefv, o3coefg)

Calculates ozone stress for a single point using Lombardozzi 2015 formulation.
Returns (o3coefv, o3coefg) — coefficients for photosynthesis and conductance.

Ported from `CalcOzoneStressLombardozzi2015OnePoint` in OzoneMod.F90.
"""
function calc_ozone_stress_lombardozzi2015_one_point(;
        pft_type::Int,
        o3uptake::Real,
        woody_pft::Vector{<:Real})

    if o3uptake == 0.0
        # No o3 damage if no o3 uptake
        o3coefv = 1.0
        o3coefg = 1.0
    else
        # Determine parameter values for this pft
        if pft_type > 3
            if woody_pft[pft_type] == 0.0
                photoInt   = O3_NONWOODY_PHOTO_INT
                photoSlope = O3_NONWOODY_PHOTO_SLOPE
                condInt    = O3_NONWOODY_COND_INT
                condSlope  = O3_NONWOODY_COND_SLOPE
            else
                photoInt   = O3_BROADLEAF_PHOTO_INT
                photoSlope = O3_BROADLEAF_PHOTO_SLOPE
                condInt    = O3_BROADLEAF_COND_INT
                condSlope  = O3_BROADLEAF_COND_SLOPE
            end
        else
            photoInt   = O3_NEEDLELEAF_PHOTO_INT
            photoSlope = O3_NEEDLELEAF_PHOTO_SLOPE
            condInt    = O3_NEEDLELEAF_COND_INT
            condSlope  = O3_NEEDLELEAF_COND_SLOPE
        end

        # Apply parameter values to compute o3 coefficients
        o3coefv = max(0.0, min(1.0, photoInt + photoSlope * o3uptake))
        o3coefg = max(0.0, min(1.0, condInt  + condSlope  * o3uptake))
    end

    return (o3coefv, o3coefg)
end

# ==========================================================================
# Falk stress
# ==========================================================================

"""
    calc_ozone_stress_falk!(oz, mask_exposedvegp, mask_noexposedvegp,
                            bounds, patchdata, woody_pft; is_time_to_run_luna)

Calculate ozone stress using the Falk formulation.
Only runs when it is time to run LUNA.

Ported from `CalcOzoneStressFalk` in OzoneMod.F90.
"""
function calc_ozone_stress_falk!(oz::OzoneData,
                                 mask_exposedvegp::BitVector,
                                 mask_noexposedvegp::BitVector,
                                 bounds::UnitRange{Int},
                                 patchdata::PatchData,
                                 woody_pft::Vector{<:Real};
                                 is_time_to_run_luna::Bool = false)
    if !is_time_to_run_luna
        return nothing
    end

    for p in bounds
        if mask_exposedvegp[p]
            pft_type = patchdata.itype[p]

            # Ozone stress for shaded leaves
            oz.o3coefjmaxsha_patch[p] = calc_ozone_stress_falk_one_point(
                pft_type = pft_type,
                o3uptake = oz.o3uptakesha_patch[p],
                o3coefjmax = oz.o3coefjmaxsha_patch[p],
                woody_pft = woody_pft)

            # Ozone stress for sunlit leaves
            oz.o3coefjmaxsun_patch[p] = calc_ozone_stress_falk_one_point(
                pft_type = pft_type,
                o3uptake = oz.o3uptakesun_patch[p],
                o3coefjmax = oz.o3coefjmaxsun_patch[p],
                woody_pft = woody_pft)

        elseif mask_noexposedvegp[p]
            # Reset to 1 over non-exposed veg points
            oz.o3coefjmaxsha_patch[p] = 1.0
            oz.o3coefjmaxsun_patch[p] = 1.0
        end
    end
    return nothing
end

"""
    calc_ozone_stress_falk_one_point(; pft_type, o3uptake, o3coefjmax, woody_pft)
        -> Float64

Calculates ozone stress for a single point using the Falk formulation.
Returns the updated o3coefjmax coefficient.

Ported from `CalcOzoneStressFalkOnePoint` in OzoneMod.F90.
"""
function calc_ozone_stress_falk_one_point(;
        pft_type::Int,
        o3uptake::Real,
        o3coefjmax::Real,
        woody_pft::Vector{<:Real})

    if o3uptake == 0.0
        return 1.0
    else
        # Determine parameter values for this pft
        if pft_type > 3
            if woody_pft[pft_type] == 0.0
                jmaxInt   = O3_NONWOODY_JMAX_INT
                jmaxSlope = O3_NONWOODY_JMAX_SLOPE
            else
                jmaxInt   = O3_BROADLEAF_JMAX_INT
                jmaxSlope = O3_BROADLEAF_JMAX_SLOPE
            end
        else
            jmaxInt   = O3_NEEDLELEAF_JMAX_INT
            jmaxSlope = O3_NEEDLELEAF_JMAX_SLOPE
        end

        # Apply parameter values to compute o3 coefficient
        return max(0.0, min(1.0, jmaxInt + jmaxSlope * o3uptake))
    end
end
