# ==========================================================================
# Ported from: src/main/lnd2atmType.F90
# Land-to-atmosphere data types and initialization
# ==========================================================================

# --------------------------------------------------------------------------
# Parameters type
# --------------------------------------------------------------------------

"""
    Lnd2AtmParamsData

Parameters controlling land-to-atmosphere coupling:
ice runoff conversion behavior.

Ported from `lnd2atm_params_type` in `lnd2atmType.F90`.
"""
Base.@kwdef mutable struct Lnd2AtmParamsData
    # true => ice runoff generated from non-glacier columns and glacier columns
    # outside icesheet regions is converted to liquid, with an appropriate
    # sensible heat flux
    melt_non_icesheet_ice_runoff::Bool = false
end

# --------------------------------------------------------------------------
# Main lnd2atm data type
# --------------------------------------------------------------------------

"""
    Lnd2AtmData

Land-to-atmosphere data. Contains gridcell-level flux fields sent from the
land model to the atmosphere, plus column-level sensible heat from ice-to-liquid
conversion.

Ported from `lnd2atm_type` in `lnd2atmType.F90`.
"""
Base.@kwdef mutable struct Lnd2AtmData{FT<:AbstractFloat}
    params::Lnd2AtmParamsData = Lnd2AtmParamsData()

    # lnd->atm (gridcell-level)
    t_rad_grc                      ::Vector{FT} = Float64[]   # radiative temperature (K)
    t_ref2m_grc                    ::Vector{FT} = Float64[]   # 2m surface air temperature (K)
    u_ref10m_grc                   ::Vector{FT} = Float64[]   # 10m surface wind speed (m/s)
    albd_grc                       ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # surface albedo (direct) (numrad)
    albi_grc                       ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # surface albedo (diffuse) (numrad)
    taux_grc                       ::Vector{FT} = Float64[]   # wind stress: e-w (kg/m/s**2)
    tauy_grc                       ::Vector{FT} = Float64[]   # wind stress: n-s (kg/m/s**2)
    eflx_lh_tot_grc                ::Vector{FT} = Float64[]   # total latent HF (W/m**2) [+ to atm]
    eflx_sh_tot_grc                ::Vector{FT} = Float64[]   # total sensible HF (W/m**2) [+ to atm]
    eflx_sh_precip_conversion_grc  ::Vector{FT} = Float64[]   # sensible HF from precip conversion (W/m**2) [+ to atm]
    eflx_lwrad_out_grc             ::Vector{FT} = Float64[]   # IR (longwave) radiation (W/m**2)
    fsa_grc                        ::Vector{FT} = Float64[]   # solar rad absorbed (total) (W/m**2)
    z0m_grc                        ::Vector{FT} = Float64[]   # roughness length, momentum (m)
    net_carbon_exchange_grc        ::Vector{FT} = Float64[]   # net CO2 flux (kg CO2/m**2/s) [+ to atm]
    nem_grc                        ::Vector{FT} = Float64[]   # gridcell average net methane correction to CO2 flux (g C/m^2/s)
    ram1_grc                       ::Vector{FT} = Float64[]   # aerodynamical resistance (s/m)
    fv_grc                         ::Vector{FT} = Float64[]   # friction velocity (m/s) (for dust model)
    flxdst_grc                     ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # dust flux (size bins)
    ddvel_grc                      ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # dry deposition velocities
    flxvoc_grc                     ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # VOC flux (size bins)
    fireflx_grc                    ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # Wild Fire Emissions
    fireztop_grc                   ::Vector{FT} = Float64[]   # Wild Fire Emissions vertical distribution top
    ch4_surf_flux_tot_grc          ::Vector{FT} = Float64[]   # net CH4 flux (kg C/m**2/s) [+ to atm]

    # lnd->atm (column-level)
    eflx_sh_ice_to_liq_col         ::Vector{FT} = Float64[]   # sensible HF from ice runoff to liquid conversion (W/m**2) [+ to atm]
end

# --------------------------------------------------------------------------
# Parameters constructor (mirrors lnd2atm_params_constructor)
# --------------------------------------------------------------------------

"""
    lnd2atm_params_init!(params; melt_non_icesheet_ice_runoff)

Initialize lnd2atm parameters.

Ported from `lnd2atm_params_constructor` in `lnd2atmType.F90`.
"""
function lnd2atm_params_init!(params::Lnd2AtmParamsData;
        melt_non_icesheet_ice_runoff::Bool)
    params.melt_non_icesheet_ice_runoff = melt_non_icesheet_ice_runoff
    nothing
end

# --------------------------------------------------------------------------
# Allocation (mirrors InitAllocate)
# --------------------------------------------------------------------------

"""
    lnd2atm_init!(l2a, ng, nc;
        n_megan_comps=0, n_fire_emis_comps=0, n_drydep=0)

Allocate and zero-initialize all arrays in `Lnd2AtmData`.
- `ng`: number of gridcells
- `nc`: number of columns
- `n_megan_comps`: number of MEGAN VOC compounds (shr_megan_mechcomps_n)
- `n_fire_emis_comps`: number of fire emission compounds (shr_fire_emis_mechcomps_n)
- `n_drydep`: number of dry deposition species

Ported from `InitAllocate` in `lnd2atmType.F90`.
"""
function lnd2atm_init!(l2a::Lnd2AtmData, ng::Int, nc::Int;
        n_megan_comps::Int = 0,
        n_fire_emis_comps::Int = 0,
        n_drydep::Int = 0)
    ival = 0.0

    # gridcell-level 1D fields
    l2a.t_rad_grc                     = fill(ival, ng)
    l2a.t_ref2m_grc                   = fill(ival, ng)
    l2a.u_ref10m_grc                  = fill(ival, ng)
    l2a.taux_grc                      = fill(ival, ng)
    l2a.tauy_grc                      = fill(ival, ng)
    l2a.eflx_lwrad_out_grc            = fill(ival, ng)
    l2a.eflx_sh_tot_grc               = fill(ival, ng)
    l2a.eflx_sh_precip_conversion_grc = fill(ival, ng)
    l2a.eflx_lh_tot_grc               = fill(ival, ng)
    l2a.fsa_grc                       = fill(ival, ng)
    l2a.z0m_grc                       = fill(ival, ng)
    l2a.net_carbon_exchange_grc       = fill(ival, ng)
    l2a.nem_grc                       = fill(ival, ng)
    l2a.ram1_grc                      = fill(ival, ng)
    l2a.fv_grc                        = fill(ival, ng)
    l2a.ch4_surf_flux_tot_grc         = fill(ival, ng)

    # gridcell-level 2D fields
    l2a.albd_grc                      = fill(ival, ng, NUMRAD)
    l2a.albi_grc                      = fill(ival, ng, NUMRAD)
    l2a.flxdst_grc                    = fill(ival, ng, NDST)

    # column-level
    l2a.eflx_sh_ice_to_liq_col        = fill(ival, nc)

    # conditional allocations
    if n_megan_comps > 0
        l2a.flxvoc_grc = fill(ival, ng, n_megan_comps)
    end
    if n_fire_emis_comps > 0
        l2a.fireflx_grc  = fill(ival, ng, n_fire_emis_comps)
        l2a.fireztop_grc = fill(ival, ng)
    end
    if n_drydep > 0
        l2a.ddvel_grc = fill(ival, ng, n_drydep)
    end

    nothing
end

# --------------------------------------------------------------------------
# ReadNamelist stub (Fortran I/O not applicable)
# --------------------------------------------------------------------------

"""
    lnd2atm_read_namelist!(l2a; melt_non_icesheet_ice_runoff=false)

Set parameters from keyword arguments (replaces Fortran namelist read).

Ported from `ReadNamelist` in `lnd2atmType.F90`.
"""
function lnd2atm_read_namelist!(l2a::Lnd2AtmData;
        melt_non_icesheet_ice_runoff::Bool = false)
    lnd2atm_params_init!(l2a.params;
        melt_non_icesheet_ice_runoff = melt_non_icesheet_ice_runoff)
    nothing
end

# --------------------------------------------------------------------------
# InitHistory stub (Fortran history output not applicable in Julia)
# --------------------------------------------------------------------------

"""
    lnd2atm_init_history!(l2a, ng, nc)

Stub for history field registration. No-op in Julia port.

Ported from `InitHistory` in `lnd2atmType.F90`.
"""
function lnd2atm_init_history!(l2a::Lnd2AtmData{FT}, ng::Int, nc::Int) where {FT}
    nothing
end

# --------------------------------------------------------------------------
# Clean (mirrors deallocate logic)
# --------------------------------------------------------------------------

"""
    lnd2atm_clean!(l2a)

Reset all array fields to empty. Mirrors Fortran deallocation.
"""
function lnd2atm_clean!(l2a::Lnd2AtmData{FT}) where {FT}
    l2a.t_rad_grc                     = FT[]
    l2a.t_ref2m_grc                   = FT[]
    l2a.u_ref10m_grc                  = FT[]
    l2a.albd_grc                      = Matrix{FT}(undef, 0, 0)
    l2a.albi_grc                      = Matrix{FT}(undef, 0, 0)
    l2a.taux_grc                      = FT[]
    l2a.tauy_grc                      = FT[]
    l2a.eflx_lh_tot_grc               = FT[]
    l2a.eflx_sh_tot_grc               = FT[]
    l2a.eflx_sh_precip_conversion_grc = FT[]
    l2a.eflx_lwrad_out_grc            = FT[]
    l2a.fsa_grc                       = FT[]
    l2a.z0m_grc                       = FT[]
    l2a.net_carbon_exchange_grc       = FT[]
    l2a.nem_grc                       = FT[]
    l2a.ram1_grc                      = FT[]
    l2a.fv_grc                        = FT[]
    l2a.flxdst_grc                    = Matrix{FT}(undef, 0, 0)
    l2a.ddvel_grc                     = Matrix{FT}(undef, 0, 0)
    l2a.flxvoc_grc                    = Matrix{FT}(undef, 0, 0)
    l2a.fireflx_grc                   = Matrix{FT}(undef, 0, 0)
    l2a.fireztop_grc                  = FT[]
    l2a.ch4_surf_flux_tot_grc         = FT[]
    l2a.eflx_sh_ice_to_liq_col        = FT[]
    nothing
end

# --------------------------------------------------------------------------
# InitForTesting (mirrors atm2lnd pattern)
# --------------------------------------------------------------------------

"""
    lnd2atm_init_for_testing!(l2a, ng, nc; params=nothing, kwargs...)

Initialize for unit testing. Allocates arrays and optionally sets params.
"""
function lnd2atm_init_for_testing!(l2a::Lnd2AtmData, ng::Int, nc::Int;
        params::Union{Lnd2AtmParamsData, Nothing} = nothing,
        n_megan_comps::Int = 0,
        n_fire_emis_comps::Int = 0,
        n_drydep::Int = 0)
    lnd2atm_init!(l2a, ng, nc;
        n_megan_comps = n_megan_comps,
        n_fire_emis_comps = n_fire_emis_comps,
        n_drydep = n_drydep)
    if params !== nothing
        l2a.params = params
    else
        l2a.params = Lnd2AtmParamsData(
            melt_non_icesheet_ice_runoff = false)
    end
    nothing
end
