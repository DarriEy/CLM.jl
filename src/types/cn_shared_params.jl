# ==========================================================================
# Ported from: src/biogeochem/CNSharedParamsMod.F90
# Parameters shared by the Carbon Nitrogen Biogeochemistry modules
# ==========================================================================

"""
    CNSharedParamsData

Shared CN/BGC parameters. Holds biogeochemical rate constants and control
flags used across multiple CN modules.

Ported from `CNParamsShareType` and module-level variables in
`CNSharedParamsMod.F90`.
"""
Base.@kwdef mutable struct CNSharedParamsData{FT<:Real}
    # --- Parameters from CNParamsShareType ---
    Q10                              ::FT = 1.5     # temperature dependence
    minpsi                           ::FT = -10.0   # minimum soil water potential for heterotrophic resp (MPa)
    maxpsi                           ::FT = -0.1    # maximum soil water potential for heterotrophic resp (MPa)
    rf_cwdl2                         ::FT = 0.0     # respiration fraction in CWD to litter2 transition (frac)
    tau_cwd                          ::FT = 0.0     # fragmentation rate constant CWD (1/yr)
    cwd_flig                         ::FT = 0.0     # lignin fraction of coarse woody debris
    froz_q10                         ::FT = 1.5     # separate q10 for frozen soil respiration rates
    decomp_depth_efolding            ::FT = 0.5     # e-folding depth for reduction in decomposition (m)
    mino2lim                         ::FT = 0.0     # minimum anaerobic decomposition rate as fraction of potential aerobic rate
    organic_max                      ::FT = 0.0     # organic matter content (kg/m3) where soil acts like peat
    constrain_stress_deciduous_onset ::Bool    = false   # if true, use additional constraint on stress deciduous onset trigger

    # --- Module-level variables ---
    use_matrixcn                     ::Bool    = false   # true => use cn matrix solution
    use_fun                          ::Bool    = false   # Use the FUN2.0 model
    nlev_soildecomp_standard         ::Int     = 5       # standard number of soil decomposition levels
    upper_soil_layer                 ::Int     = -1      # upper soil layer for 10-day average in CNPhenology
end

"""
    cn_shared_params_read_netcdf!(params; kwargs...)

Read CN shared parameters from a Dict-like source (replaces NetCDF file reading).
Corresponds to `CNParamsReadShared_netcdf` in the Fortran source.

Each keyword argument maps to a NetCDF variable name from the Fortran source:
- `q10_mr` → Q10
- `minpsi_hr` → minpsi
- `maxpsi_hr` → maxpsi
- `rf_cwdl2` → rf_cwdl2
- `tau_cwd` → tau_cwd
- `cwd_flig` → cwd_flig
- `decomp_depth_efolding` → decomp_depth_efolding
- `froz_q10` → froz_q10
- `mino2lim` → mino2lim
- `organic_max` → organic_max
"""
function cn_shared_params_read_netcdf!(params::CNSharedParamsData;
                                       q10_mr::Real,
                                       minpsi_hr::Real,
                                       maxpsi_hr::Real,
                                       rf_cwdl2::Real,
                                       tau_cwd::Real,
                                       cwd_flig::Real,
                                       decomp_depth_efolding::Real,
                                       froz_q10::Real,
                                       mino2lim::Real,
                                       organic_max::Real)
    params.Q10                   = q10_mr
    params.minpsi                = minpsi_hr
    params.maxpsi                = maxpsi_hr
    params.rf_cwdl2              = rf_cwdl2
    params.tau_cwd               = tau_cwd
    params.cwd_flig              = cwd_flig
    params.decomp_depth_efolding = decomp_depth_efolding
    params.froz_q10              = froz_q10
    params.mino2lim              = mino2lim
    params.organic_max           = organic_max
    return nothing
end

"""
    cn_shared_params_read_namelist!(params; constrain_stress_deciduous_onset=false)

Read CN shared namelist parameters.
Corresponds to `CNParamsReadShared_namelist` in the Fortran source.
"""
function cn_shared_params_read_namelist!(params::CNSharedParamsData;
                                         constrain_stress_deciduous_onset::Bool=false)
    params.constrain_stress_deciduous_onset = constrain_stress_deciduous_onset
    return nothing
end

"""
    cn_shared_params_read!(params; netcdf_params, constrain_stress_deciduous_onset=false)

Read all CN shared parameters (NetCDF + namelist).
Corresponds to `CNParamsReadShared` in the Fortran source.
"""
function cn_shared_params_read!(params::CNSharedParamsData;
                                 q10_mr::Real,
                                 minpsi_hr::Real,
                                 maxpsi_hr::Real,
                                 rf_cwdl2::Real,
                                 tau_cwd::Real,
                                 cwd_flig::Real,
                                 decomp_depth_efolding::Real,
                                 froz_q10::Real,
                                 mino2lim::Real,
                                 organic_max::Real,
                                 constrain_stress_deciduous_onset::Bool=false)
    cn_shared_params_read_netcdf!(params;
        q10_mr=q10_mr, minpsi_hr=minpsi_hr, maxpsi_hr=maxpsi_hr,
        rf_cwdl2=rf_cwdl2, tau_cwd=tau_cwd, cwd_flig=cwd_flig,
        decomp_depth_efolding=decomp_depth_efolding,
        froz_q10=froz_q10, mino2lim=mino2lim, organic_max=organic_max)
    cn_shared_params_read_namelist!(params;
        constrain_stress_deciduous_onset=constrain_stress_deciduous_onset)
    return nothing
end

"""
    cn_shared_params_set_soil_depth!(params, depth_to_layer_fn)

Set the upper soil layer for CNPhenology by finding the soil layer
containing depth 0.12 m.
Corresponds to `CNParamsSetSoilDepth` in the Fortran source.

`depth_to_layer_fn` should be a callable `f(depth) -> layer_index`,
corresponding to `find_soil_layer_containing_depth` in the Fortran source.
"""
function cn_shared_params_set_soil_depth!(params::CNSharedParamsData,
                                           depth_to_layer_fn)
    params.upper_soil_layer = depth_to_layer_fn(0.12)
    return nothing
end
