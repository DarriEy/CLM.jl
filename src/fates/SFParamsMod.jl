# SFParamsMod.jl
# Julia port of FATES src/fates/fire/SFParamsMod.F90
#
# The SPITFIRE fire-model parameters (fuel-energy content, mineral fractions,
# moisture-of-extinction / drying ratios, fire-spread + fire-intensity
# coefficients, etc.), plus their registration with / retrieval from the FATES
# parameters interface (fates_parameters_type) and init/consistency checks.
#
# In Fortran these are module-level `protected, save` scalars and per-fuel-class
# / per-CWD arrays. Here we collect every settable value into a single mutable
# holder (`SFParams`, a const Ref to an `sf_params_type` instance), mirroring the
# Ref-holder pattern used by EDParamsMod / FatesGlobals. SoA: per-fuel-class
# arrays -> `Vector{Float64}` of length `num_fuel_classes`; the per-CWD fraction
# -> `Vector{Float64}` of length `ncwd`.
#
# The parameter-file READ stays stubbed (no NetCDF), exactly like the foundation
# FatesParametersInterface: SpitFireRegisterParams! registers the names + shapes,
# SpitFireReceiveParams! retrieves values from an already-populated
# fates_parameters_type, and the parameter names are preserved verbatim
# ("fates_fire_*" / "fates_frag_*").
#
# Deps: FatesConstantsMod (fates_r8 -> Float64, fates_check_param_set,
#       fates_unset_r8), FatesFuelClassesMod (num_fuel_classes), FatesLitterMod
#       (ncwd), FatesParametersInterface (fates_parameters_type,
#       param_string_length, dimension_*, RegisterParameter!/RetrieveParameter*),
#       FatesGlobals (fates_log, fates_endrun).

# ---------------------------------------------------------------------------
# Module-private parameters (Fortran `parameter, private`) -> const.
# ---------------------------------------------------------------------------

# The minimum reasonable fire intensity threshold [kW/m]
const min_fire_threshold = 0.0001

# ---------------------------------------------------------------------------
# Parameter NAMES (Fortran character `parameter`s) -> const Strings.
# Preserve the FATES parameter-file names verbatim.
# ---------------------------------------------------------------------------
const SF_name_fdi_alpha          = "fates_fire_fdi_alpha"
const SF_name_miner_total        = "fates_fire_miner_total"
const SF_name_fuel_energy        = "fates_fire_fuel_energy"
const SF_name_part_dens          = "fates_fire_part_dens"
const SF_name_miner_damp         = "fates_fire_miner_damp"
const SF_name_max_durat          = "fates_fire_max_durat"
const SF_name_durat_slope        = "fates_fire_durat_slope"
const SF_name_drying_ratio       = "fates_fire_drying_ratio"
const SF_name_fire_threshold     = "fates_fire_threshold"
const SF_name_CWD_frac           = "fates_frag_cwd_frac"
const SF_name_max_decomp         = "fates_frag_maxdecomp"
const SF_name_SAV                = "fates_fire_SAV"
const SF_name_FBD                = "fates_fire_FBD"
const SF_name_min_moisture       = "fates_fire_min_moisture"
const SF_name_mid_moisture       = "fates_fire_mid_moisture"
const SF_name_low_moisture_Coeff = "fates_fire_low_moisture_Coeff"
const SF_name_low_moisture_Slope = "fates_fire_low_moisture_Slope"
const SF_name_mid_moisture_Coeff = "fates_fire_mid_moisture_Coeff"
const SF_name_mid_moisture_Slope = "fates_fire_mid_moisture_Slope"

# ---------------------------------------------------------------------------
# sf_params_type — the settable SPITFIRE parameters.
#
# Holds every Fortran module variable that is set inside SpitFireParamsInit /
# SpitFireReceiveParams. Scalars default to the unset sentinel; per-fuel-class
# (num_fuel_classes) and per-CWD (ncwd) arrays are sized + filled with the
# sentinel. The const `SFParams` Ref below is the live module-global instance
# the API mutates (mirrors EDParams).
# ---------------------------------------------------------------------------
Base.@kwdef mutable struct sf_params_type
    # --- scalars ---
    SF_val_fdi_alpha::Float64       = fates_unset_r8
    SF_val_miner_total::Float64     = fates_unset_r8
    SF_val_fuel_energy::Float64     = fates_unset_r8
    SF_val_part_dens::Float64       = fates_unset_r8
    SF_val_miner_damp::Float64      = fates_unset_r8
    SF_val_max_durat::Float64       = fates_unset_r8
    SF_val_durat_slope::Float64     = fates_unset_r8
    SF_val_drying_ratio::Float64    = fates_unset_r8
    # threshold for fires that spread or go out. kW/m (Pyne 1996)
    SF_val_fire_threshold::Float64  = fates_unset_r8

    # --- per-CWD (ncwd) ---
    SF_val_CWD_frac::Vector{Float64} = fill(fates_unset_r8, ncwd)

    # --- per-fuel-class (num_fuel_classes) ---
    SF_val_max_decomp::Vector{Float64}        = fill(fates_unset_r8, num_fuel_classes)
    SF_val_SAV::Vector{Float64}               = fill(fates_unset_r8, num_fuel_classes)
    SF_val_FBD::Vector{Float64}               = fill(fates_unset_r8, num_fuel_classes)
    SF_val_min_moisture::Vector{Float64}      = fill(fates_unset_r8, num_fuel_classes)
    SF_val_mid_moisture::Vector{Float64}      = fill(fates_unset_r8, num_fuel_classes)
    SF_val_low_moisture_Coeff::Vector{Float64} = fill(fates_unset_r8, num_fuel_classes)
    SF_val_low_moisture_Slope::Vector{Float64} = fill(fates_unset_r8, num_fuel_classes)
    SF_val_mid_moisture_Coeff::Vector{Float64} = fill(fates_unset_r8, num_fuel_classes)
    SF_val_mid_moisture_Slope::Vector{Float64} = fill(fates_unset_r8, num_fuel_classes)
end

"""
    SFParams

Live module-global instance of [`sf_params_type`](@ref), held in a `Ref` so the
SPITFIRE parameters can be re-initialized at runtime (mirrors the FatesGlobals /
EDParams Ref-holder pattern). Access the current values with [`sf_params`](@ref).
"""
const SFParams = Ref{sf_params_type}(sf_params_type())

"""
    sf_params() -> sf_params_type

Return the live SPITFIRE parameter instance.
"""
sf_params() = SFParams[]

# =====================================================================================

"""
    SpitFireParamsInit!()

Initialize all SPITFIRE parameters to the unset sentinel (Fortran initializes to
NaN to guarantee valid values come back from the host). Resets the module-global
[`SFParams`](@ref) instance.
"""
function SpitFireParamsInit!()
    SFParams[] = sf_params_type()
    return nothing
end

# =====================================================================================

"""
    SpitFireRegisterParams!(fates_params::fates_parameters_type)

Register all SPITFIRE parameters (scalars + per-CWD + per-fuel-class) with the
FATES parameter reader. Calls [`SpitFireParamsInit!`](@ref) first, then the
scalar / NCWD / NFSC registration helpers.
"""
function SpitFireRegisterParams!(fates_params::fates_parameters_type)
    SpitFireParamsInit!()
    SpitFireRegisterScalars!(fates_params)
    SpitFireRegisterNCWD!(fates_params)
    SpitFireRegisterNFSC!(fates_params)
    return nothing
end

# =====================================================================================

"""
    SpitFireReceiveParams!(fates_params::fates_parameters_type)

Retrieve all SPITFIRE parameter values from an already-populated
`fates_parameters_type` into the module-global [`SFParams`](@ref) instance.
"""
function SpitFireReceiveParams!(fates_params::fates_parameters_type)
    SpitFireReceiveScalars!(fates_params)
    SpitFireReceiveNCWD!(fates_params)
    SpitFireReceiveNFSC!(fates_params)
    return nothing
end

# =====================================================================================

"""
    SpitFireRegisterScalars!(fates_params::fates_parameters_type)

Register the scalar-shaped SPITFIRE parameters.
"""
function SpitFireRegisterScalars!(fates_params::fates_parameters_type)
    dim_names_scalar = [dimension_name_scalar]

    scalar_names = (
        SF_name_fdi_alpha,
        SF_name_miner_total,
        SF_name_fuel_energy,
        SF_name_part_dens,
        SF_name_miner_damp,
        SF_name_max_durat,
        SF_name_durat_slope,
        SF_name_drying_ratio,
        SF_name_fire_threshold,
    )
    for nm in scalar_names
        RegisterParameter!(fates_params, nm, dimension_shape_scalar, dim_names_scalar)
    end
    return nothing
end

# =====================================================================================

"""
    SpitFireReceiveScalars!(fates_params::fates_parameters_type)

Retrieve the scalar-shaped SPITFIRE parameters.
"""
function SpitFireReceiveScalars!(fates_params::fates_parameters_type)
    p = SFParams[]
    p.SF_val_fdi_alpha      = RetrieveParameterScalar(fates_params, SF_name_fdi_alpha)
    p.SF_val_miner_total    = RetrieveParameterScalar(fates_params, SF_name_miner_total)
    p.SF_val_fuel_energy    = RetrieveParameterScalar(fates_params, SF_name_fuel_energy)
    p.SF_val_part_dens      = RetrieveParameterScalar(fates_params, SF_name_part_dens)
    p.SF_val_miner_damp     = RetrieveParameterScalar(fates_params, SF_name_miner_damp)
    p.SF_val_max_durat      = RetrieveParameterScalar(fates_params, SF_name_max_durat)
    p.SF_val_durat_slope    = RetrieveParameterScalar(fates_params, SF_name_durat_slope)
    p.SF_val_drying_ratio   = RetrieveParameterScalar(fates_params, SF_name_drying_ratio)
    p.SF_val_fire_threshold = RetrieveParameterScalar(fates_params, SF_name_fire_threshold)
    return nothing
end

# =====================================================================================

"""
    SpitFireRegisterNCWD!(fates_params::fates_parameters_type)

Register the per-CWD (1-D, `ncwd`) SPITFIRE parameter (`fates_frag_cwd_frac`).
"""
function SpitFireRegisterNCWD!(fates_params::fates_parameters_type)
    dim_names_cwd = [dimension_name_cwd]
    RegisterParameter!(fates_params, SF_name_CWD_frac, dimension_shape_1d, dim_names_cwd)
    return nothing
end

# =====================================================================================

"""
    SpitFireReceiveNCWD!(fates_params::fates_parameters_type)

Retrieve the per-CWD (1-D, `ncwd`) SPITFIRE parameter.
"""
function SpitFireReceiveNCWD!(fates_params::fates_parameters_type)
    p = SFParams[]
    RetrieveParameter1D(fates_params, SF_name_CWD_frac, p.SF_val_CWD_frac)
    return nothing
end

# =====================================================================================

"""
    SpitFireRegisterNFSC!(fates_params::fates_parameters_type)

Register the per-fuel-class (1-D, `num_fuel_classes`) SPITFIRE parameters. Note
that the Fortran uses the `fates_litterclass` dimension (`dimension_name_fsc`)
for these.
"""
function SpitFireRegisterNFSC!(fates_params::fates_parameters_type)
    dim_names = [dimension_name_fsc]

    nfsc_names = (
        SF_name_SAV,
        SF_name_FBD,
        SF_name_min_moisture,
        SF_name_mid_moisture,
        SF_name_low_moisture_Coeff,
        SF_name_low_moisture_Slope,
        SF_name_mid_moisture_Coeff,
        SF_name_mid_moisture_Slope,
        SF_name_max_decomp,
    )
    for nm in nfsc_names
        RegisterParameter!(fates_params, nm, dimension_shape_1d, dim_names)
    end
    return nothing
end

# =====================================================================================

"""
    SpitFireReceiveNFSC!(fates_params::fates_parameters_type)

Retrieve the per-fuel-class (1-D, `num_fuel_classes`) SPITFIRE parameters.
"""
function SpitFireReceiveNFSC!(fates_params::fates_parameters_type)
    p = SFParams[]
    RetrieveParameter1D(fates_params, SF_name_SAV,                p.SF_val_SAV)
    RetrieveParameter1D(fates_params, SF_name_FBD,                p.SF_val_FBD)
    RetrieveParameter1D(fates_params, SF_name_min_moisture,       p.SF_val_min_moisture)
    RetrieveParameter1D(fates_params, SF_name_mid_moisture,       p.SF_val_mid_moisture)
    RetrieveParameter1D(fates_params, SF_name_low_moisture_Coeff, p.SF_val_low_moisture_Coeff)
    RetrieveParameter1D(fates_params, SF_name_low_moisture_Slope, p.SF_val_low_moisture_Slope)
    RetrieveParameter1D(fates_params, SF_name_mid_moisture_Coeff, p.SF_val_mid_moisture_Coeff)
    RetrieveParameter1D(fates_params, SF_name_mid_moisture_Slope, p.SF_val_mid_moisture_Slope)
    RetrieveParameter1D(fates_params, SF_name_max_decomp,         p.SF_val_max_decomp)
    return nothing
end

# =====================================================================================

"""
    SpitFireCheckParams(is_master::Bool)

Logical consistency checks on the user-supplied SPITFIRE parameters (only logs /
acts if `is_master`). Cross-compares parameters and aborts if they do not make
sense:
  * max-decomposition rates must be ≥ 0;
  * the per-CWD fractions must sum to unity (a tiny residual is corrected into
    the largest pool for tight mass conservation; a large mismatch is a fatal
    error);
  * the fire-intensity threshold must be set and above `min_fire_threshold`.

Mutates the module-global [`SFParams`](@ref) when applying the CWD-fraction
correction (matching the Fortran in-place adjustment of `SF_val_CWD_frac`).
"""
function SpitFireCheckParams(is_master::Bool)
    is_master || return nothing

    p = SFParams[]

    # Decomposition rates should not be less than zero.
    for c in 1:num_fuel_classes
        if p.SF_val_max_decomp[c] < 0.0
            println("Decomposition rates should be >0")
            println("c = ", c, " SF_val_max_decomp(c) = ", p.SF_val_max_decomp[c])
            fates_endrun("SFParamsMod: SF_val_max_decomp < 0")
        end
    end

    # Check if the CWD fraction sums to unity; if it is not wayyy off, add a
    # small correction to the largest pool (important for tight mass
    # conservation checks).
    cwd_sum = sum(@view p.SF_val_CWD_frac[1:ncwd])
    if abs(1.0 - cwd_sum) > 1.0e-5
        println("The CWD fractions from index 1:4 must sum to unity")
        println("SF_val_CWD_frac(1:ncwd) = ", p.SF_val_CWD_frac[1:ncwd])
        println("error = ", 1.0 - cwd_sum)
        fates_endrun("SFParamsMod: CWD fractions do not sum to unity")
    else
        correction = 1.0 - cwd_sum
        corr_id = argmax(@view p.SF_val_CWD_frac[1:ncwd])
        p.SF_val_CWD_frac[corr_id] = p.SF_val_CWD_frac[corr_id] + correction
    end

    # Check to see if the fire threshold is above the minimum and set at all.
    if p.SF_val_fire_threshold < min_fire_threshold ||
       p.SF_val_fire_threshold > fates_check_param_set
        println("The fates_fire_threshold parameter must be set, and > ", min_fire_threshold)
        println("The value is set at :", p.SF_val_fire_threshold)
        println("Please provide a reasonable value, aborting.")
        fates_endrun("SFParamsMod: fates_fire_threshold not set or out of range")
    end

    return nothing
end
