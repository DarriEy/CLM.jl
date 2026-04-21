# ==========================================================================
# LUNA parameters type — extracted from luna.jl for early include ordering
# ==========================================================================

"""
    LunaParamsData{FT<:Real}

Parameters for the LUNA model (from `params_type` in Fortran).
Parameterized on float type `FT` for AD compatibility.
"""
mutable struct LunaParamsData{FT<:Real}
    cp25_yr2000::FT               # CO2 compensation point at 25°C at present day O2 (mol/mol)
    kc25_coef::FT                 # Michaelis-Menten const. at 25°C for CO2 (unitless)
    ko25_coef::FT                 # Michaelis-Menten const. at 25°C for O2 (unitless)
    luna_theta_cj::FT             # LUNA empirical curvature parameter for ac, aj photosynthesis co-limitation
    enzyme_turnover_daily::FT     # Daily turnover rate for photosynthetic enzyme at 25oC
    relhExp::FT                   # Impact of relative humidity on electron transport rate
    minrelh::FT                   # Minimum relative humidity for nitrogen optimization (fraction)
    jmaxb0::Vector{FT}            # Baseline proportion of nitrogen allocated for electron transport
    jmaxb1::Vector{FT}            # Coefficient for electron transport rate response to light
    wc2wjb0::Vector{FT}           # Baseline ratio of rubisco limited vs light limited rate (Wc:Wj)
end

# Default constructor: LunaParamsData() creates Float64 version with NaN scalars and empty vectors
function LunaParamsData{FT}() where {FT<:Real}
    LunaParamsData{FT}(
        FT(NaN), FT(NaN), FT(NaN), FT(NaN), FT(NaN), FT(NaN), FT(NaN),
        FT[], FT[], FT[]
    )
end
LunaParamsData() = LunaParamsData{Float64}()

"""
    luna_params_init!(lp, mxpft)

Allocate and initialize LUNA parameters with NaN.
"""
function luna_params_init!(lp::LunaParamsData{FT}, mxpft::Int) where {FT}
    lp.jmaxb0  = fill(FT(NaN), mxpft + 1)
    lp.jmaxb1  = fill(FT(NaN), mxpft + 1)
    lp.wc2wjb0 = fill(FT(NaN), mxpft + 1)
    return nothing
end

"""
    luna_params_clean!(lp)

Deallocate LUNA parameters.
"""
function luna_params_clean!(lp::LunaParamsData{FT}) where {FT}
    lp.jmaxb0  = FT[]
    lp.jmaxb1  = FT[]
    lp.wc2wjb0 = FT[]
    return nothing
end
