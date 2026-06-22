# FatesSynchronizedParamsMod.jl
# Julia port of FATES src/fates/main/FatesSynchronizedParamsMod.F90
#
# Helper for parameters synchronized between FATES and its host. FATES currently
# shares NO synchronized parameters (previously q10). The type + register/receive
# scaffolding is kept as a template for future shared parameters.
# Deps: FatesConstantsMod, FatesParametersInterface.

"""
    FatesSynchronizedParamsType

Holder for host-synchronized parameters. Currently empty (no shared params);
fields like Q10 / froz_q10 are kept commented in the Fortran as a template.
"""
Base.@kwdef mutable struct FatesSynchronizedParamsType
    # Q10::Float64 = NaN      # temperature dependence (currently unused)
    # froz_q10::Float64 = NaN # separate q10 for frozen soil respiration (unused)
end

# Module-level instance (Fortran public protected instance).
const FatesSynchronizedParamsInst = FatesSynchronizedParamsType()

# =====================================================================================

"""
    InitSynchronizedParams!(this::FatesSynchronizedParamsType)

Initialize all synchronized parameters (to NaN) to ensure valid values come back
from the host. No-op while no parameters are shared.
"""
function InitSynchronizedParams!(this::FatesSynchronizedParamsType)
    # this.Q10 = NaN
    # this.froz_q10 = NaN
    return nothing
end

# =====================================================================================

"""
    RegisterParams!(this::FatesSynchronizedParamsType, fates_params::fates_parameters_type)

Register the parameters we want the host to provide. No-op while no parameters
are shared.
"""
function RegisterParams!(this::FatesSynchronizedParamsType, fates_params::fates_parameters_type)
    InitSynchronizedParams!(this)
    RegisterParamsScalar!(this, fates_params)
    return nothing
end

# =====================================================================================

"""
    ReceiveParams!(this::FatesSynchronizedParamsType, fates_params::fates_parameters_type)

Retrieve synchronized parameter values from the host. No-op while no parameters
are shared.
"""
function ReceiveParams!(this::FatesSynchronizedParamsType, fates_params::fates_parameters_type)
    ReceiveParamsScalar!(this, fates_params)
    return nothing
end

# =====================================================================================

function RegisterParamsScalar!(this::FatesSynchronizedParamsType, fates_params::fates_parameters_type)
    InitSynchronizedParams!(this)
    # No scalar params currently registered.
    return nothing
end

# =====================================================================================

function ReceiveParamsScalar!(this::FatesSynchronizedParamsType, fates_params::fates_parameters_type)
    # name = "q10_mr"; RetrieveParameter(..., this.Q10)
    # name = "froz_q10"; RetrieveParameter(..., this.froz_q10)
    return nothing
end
