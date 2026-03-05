# ==========================================================================
# Ported from: src/biogeophys/UrbanParamsType.F90
# Urban parameters data type, input data, and initialization
# ==========================================================================

# ---------------------------------------------------------------------------
# Urban control variables (module-level in Fortran)
# ---------------------------------------------------------------------------

const URBAN_HAC_OFF        = "OFF"
const URBAN_HAC_ON         = "ON"
const URBAN_WASTEHEAT_ON   = "ON_WASTEHEAT"

Base.@kwdef mutable struct UrbanControl
    urban_hac::String          = URBAN_HAC_OFF
    urban_explicit_ac::Bool    = true
    urban_traffic::Bool        = false
    building_temp_method::Int  = 1   # 0=simple (CLM4.5), 1=prognostic (CLM5.0)
    read_namelist::Bool        = false
end

const urban_ctrl = UrbanControl()

const BUILDING_TEMP_METHOD_SIMPLE = 0
const BUILDING_TEMP_METHOD_PROG   = 1

"""
    is_simple_build_temp() -> Bool

Returns true if the simple building temperature method is being used.

Ported from `IsSimpleBuildTemp` in `UrbanParamsType.F90`.
"""
function is_simple_build_temp()
    if !urban_ctrl.read_namelist
        error("Testing on building_temp_method before urban namelist was read in")
    end
    return urban_ctrl.building_temp_method == BUILDING_TEMP_METHOD_SIMPLE
end

"""
    is_prog_build_temp() -> Bool

Returns true if the prognostic building temperature method is being used.

Ported from `IsProgBuildTemp` in `UrbanParamsType.F90`.
"""
function is_prog_build_temp()
    if !urban_ctrl.read_namelist
        error("Testing on building_temp_method before urban namelist was read in")
    end
    return urban_ctrl.building_temp_method == BUILDING_TEMP_METHOD_PROG
end

# ---------------------------------------------------------------------------
# Urban input data type (urbinp_type in Fortran)
# Holds raw input data read from the surface dataset, indexed by
# (gridcell, density_type) or (gridcell, density_type, band/level).
# ---------------------------------------------------------------------------

"""
    UrbanInputData

Raw urban input data read from the surface dataset.
Indexed by (gridcell, density_type) for 2D fields and
(gridcell, density_type, band_or_level) for 3D fields.

Ported from `urbinp_type` in `UrbanParamsType.F90`.
"""
Base.@kwdef mutable struct UrbanInputData{FT<:AbstractFloat}
    canyon_hwr      ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    wtlunit_roof    ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    wtroad_perv     ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    em_roof         ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    em_improad      ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    em_perroad      ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    em_wall         ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    alb_roof_dir    ::Array{FT,3} = Array{Float64}(undef, 0, 0, 0)
    alb_roof_dif    ::Array{FT,3} = Array{Float64}(undef, 0, 0, 0)
    alb_improad_dir ::Array{FT,3} = Array{Float64}(undef, 0, 0, 0)
    alb_improad_dif ::Array{FT,3} = Array{Float64}(undef, 0, 0, 0)
    alb_perroad_dir ::Array{FT,3} = Array{Float64}(undef, 0, 0, 0)
    alb_perroad_dif ::Array{FT,3} = Array{Float64}(undef, 0, 0, 0)
    alb_wall_dir    ::Array{FT,3} = Array{Float64}(undef, 0, 0, 0)
    alb_wall_dif    ::Array{FT,3} = Array{Float64}(undef, 0, 0, 0)
    ht_roof         ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    wind_hgt_canyon ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    tk_wall         ::Array{FT,3} = Array{Float64}(undef, 0, 0, 0)
    tk_roof         ::Array{FT,3} = Array{Float64}(undef, 0, 0, 0)
    tk_improad      ::Array{FT,3} = Array{Float64}(undef, 0, 0, 0)
    cv_wall         ::Array{FT,3} = Array{Float64}(undef, 0, 0, 0)
    cv_roof         ::Array{FT,3} = Array{Float64}(undef, 0, 0, 0)
    cv_improad      ::Array{FT,3} = Array{Float64}(undef, 0, 0, 0)
    thick_wall      ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    thick_roof      ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    nlev_improad    ::Matrix{Int}      = Matrix{Int}(undef, 0, 0)
    t_building_min  ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
end

"""
    urbinp_init!(ui::UrbanInputData, ng::Int, numurbl::Int, numrad::Int, nlevurb::Int)

Allocate all fields of an `UrbanInputData` instance.

Ported from the 'initialize' branch of `UrbanInput` in `UrbanParamsType.F90`.
"""
function urbinp_init!(ui::UrbanInputData{FT}, ng::Int, numurbl::Int, numrad::Int, nlevurb::Int) where {FT}
    ui.canyon_hwr      = fill(FT(NaN), ng, numurbl)
    ui.wtlunit_roof    = fill(FT(NaN), ng, numurbl)
    ui.wtroad_perv     = fill(FT(NaN), ng, numurbl)
    ui.em_roof         = fill(FT(NaN), ng, numurbl)
    ui.em_improad      = fill(FT(NaN), ng, numurbl)
    ui.em_perroad      = fill(FT(NaN), ng, numurbl)
    ui.em_wall         = fill(FT(NaN), ng, numurbl)
    ui.alb_roof_dir    = fill(FT(NaN), ng, numurbl, numrad)
    ui.alb_roof_dif    = fill(FT(NaN), ng, numurbl, numrad)
    ui.alb_improad_dir = fill(FT(NaN), ng, numurbl, numrad)
    ui.alb_improad_dif = fill(FT(NaN), ng, numurbl, numrad)
    ui.alb_perroad_dir = fill(FT(NaN), ng, numurbl, numrad)
    ui.alb_perroad_dif = fill(FT(NaN), ng, numurbl, numrad)
    ui.alb_wall_dir    = fill(FT(NaN), ng, numurbl, numrad)
    ui.alb_wall_dif    = fill(FT(NaN), ng, numurbl, numrad)
    ui.ht_roof         = fill(FT(NaN), ng, numurbl)
    ui.wind_hgt_canyon = fill(FT(NaN), ng, numurbl)
    ui.tk_wall         = fill(FT(NaN), ng, numurbl, nlevurb)
    ui.tk_roof         = fill(FT(NaN), ng, numurbl, nlevurb)
    ui.tk_improad      = fill(FT(NaN), ng, numurbl, nlevurb)
    ui.cv_wall         = fill(FT(NaN), ng, numurbl, nlevurb)
    ui.cv_roof         = fill(FT(NaN), ng, numurbl, nlevurb)
    ui.cv_improad      = fill(FT(NaN), ng, numurbl, nlevurb)
    ui.thick_wall      = fill(FT(NaN), ng, numurbl)
    ui.thick_roof      = fill(FT(NaN), ng, numurbl)
    ui.nlev_improad    = fill(0,   ng, numurbl)
    ui.t_building_min  = fill(FT(NaN), ng, numurbl)

    return nothing
end

"""
    urbinp_clean!(ui::UrbanInputData)

Deallocate (reset to empty) all fields of an `UrbanInputData` instance.

Ported from the 'finalize' branch of `UrbanInput` in `UrbanParamsType.F90`.
"""
function urbinp_clean!(ui::UrbanInputData{FT}) where {FT}
    ui.canyon_hwr      = Matrix{FT}(undef, 0, 0)
    ui.wtlunit_roof    = Matrix{FT}(undef, 0, 0)
    ui.wtroad_perv     = Matrix{FT}(undef, 0, 0)
    ui.em_roof         = Matrix{FT}(undef, 0, 0)
    ui.em_improad      = Matrix{FT}(undef, 0, 0)
    ui.em_perroad      = Matrix{FT}(undef, 0, 0)
    ui.em_wall         = Matrix{FT}(undef, 0, 0)
    ui.alb_roof_dir    = Array{Float64}(undef, 0, 0, 0)
    ui.alb_roof_dif    = Array{Float64}(undef, 0, 0, 0)
    ui.alb_improad_dir = Array{Float64}(undef, 0, 0, 0)
    ui.alb_improad_dif = Array{Float64}(undef, 0, 0, 0)
    ui.alb_perroad_dir = Array{Float64}(undef, 0, 0, 0)
    ui.alb_perroad_dif = Array{Float64}(undef, 0, 0, 0)
    ui.alb_wall_dir    = Array{Float64}(undef, 0, 0, 0)
    ui.alb_wall_dif    = Array{Float64}(undef, 0, 0, 0)
    ui.ht_roof         = Matrix{FT}(undef, 0, 0)
    ui.wind_hgt_canyon = Matrix{FT}(undef, 0, 0)
    ui.tk_wall         = Array{Float64}(undef, 0, 0, 0)
    ui.tk_roof         = Array{Float64}(undef, 0, 0, 0)
    ui.tk_improad      = Array{Float64}(undef, 0, 0, 0)
    ui.cv_wall         = Array{Float64}(undef, 0, 0, 0)
    ui.cv_roof         = Array{Float64}(undef, 0, 0, 0)
    ui.cv_improad      = Array{Float64}(undef, 0, 0, 0)
    ui.thick_wall      = Matrix{FT}(undef, 0, 0)
    ui.thick_roof      = Matrix{FT}(undef, 0, 0)
    ui.nlev_improad    = Matrix{Int}(undef, 0, 0)
    ui.t_building_min  = Matrix{FT}(undef, 0, 0)

    return nothing
end

# ---------------------------------------------------------------------------
# Urban parameters type (urbanparams_type in Fortran)
# Landunit-level parameters derived from the input data.
# ---------------------------------------------------------------------------

"""
    UrbanParamsData

Landunit-level urban parameters. Derived from `UrbanInputData` during
initialization for each urban landunit.

Ported from `urbanparams_type` in `UrbanParamsType.F90`.
"""
Base.@kwdef mutable struct UrbanParamsData{FT<:AbstractFloat}
    # --- Emissivities (landunit-level 1D) ---
    wind_hgt_canyon     ::Vector{FT} = Float64[]   # lun height above road at which wind in canyon is computed (m)
    em_roof             ::Vector{FT} = Float64[]   # lun roof emissivity
    em_improad          ::Vector{FT} = Float64[]   # lun impervious road emissivity
    em_perroad          ::Vector{FT} = Float64[]   # lun pervious road emissivity
    em_wall             ::Vector{FT} = Float64[]   # lun wall emissivity

    # --- Albedos (landunit-level 2D: nlun × numrad) ---
    alb_roof_dir        ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # lun direct roof albedo
    alb_roof_dif        ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # lun diffuse roof albedo
    alb_improad_dir     ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # lun direct impervious road albedo
    alb_improad_dif     ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # lun diffuse impervious road albedo
    alb_perroad_dir     ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # lun direct pervious road albedo
    alb_perroad_dif     ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # lun diffuse pervious road albedo
    alb_wall_dir        ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # lun direct wall albedo
    alb_wall_dif        ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # lun diffuse wall albedo

    # --- Thermal properties (landunit-level 2D: nlun × nlevurb) ---
    nlev_improad        ::Vector{Int}     = Int[]       # lun number of impervious road layers (-)
    tk_wall             ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # lun thermal conductivity of wall (W/m/K)
    tk_roof             ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # lun thermal conductivity of roof (W/m/K)
    tk_improad          ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # lun thermal conductivity of impervious road (W/m/K)
    cv_wall             ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # lun heat capacity of wall (J/m^3/K)
    cv_roof             ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # lun heat capacity of roof (J/m^3/K)
    cv_improad          ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # lun heat capacity of impervious road (J/m^3/K)
    thick_wall          ::Vector{FT} = Float64[]   # lun total thickness of wall (m)
    thick_roof          ::Vector{FT} = Float64[]   # lun total thickness of roof (m)

    # --- View factors (landunit-level 1D) ---
    vf_sr               ::Vector{FT} = Float64[]   # lun view factor of sky for road
    vf_wr               ::Vector{FT} = Float64[]   # lun view factor of one wall for road
    vf_sw               ::Vector{FT} = Float64[]   # lun view factor of sky for one wall
    vf_rw               ::Vector{FT} = Float64[]   # lun view factor of road for one wall
    vf_ww               ::Vector{FT} = Float64[]   # lun view factor of opposing wall for one wall

    # --- Building / traffic parameters (landunit-level 1D) ---
    t_building_min      ::Vector{FT} = Float64[]   # lun minimum internal building air temperature (K)
    eflx_traffic_factor ::Vector{FT} = Float64[]   # lun multiplicative traffic factor for sensible heat flux (-)
end

"""
    urbanparams_init!(up::UrbanParamsData, nl::Int; nlevurb::Int=0, numrad::Int=NUMRAD)

Allocate all fields of an `UrbanParamsData` instance for `nl` landunits.
Real fields are initialized to `NaN`, integer fields to a large sentinel.

Ported from the allocation portion of `urbanparams_type%Init` in `UrbanParamsType.F90`.
"""
function urbanparams_init!(up::UrbanParamsData, nl::Int;
                            nlevurb::Int = 0, numrad::Int = NUMRAD)
    # 1D landunit fields
    up.wind_hgt_canyon     = fill(NaN, nl)
    up.em_roof             = fill(NaN, nl)
    up.em_improad          = fill(NaN, nl)
    up.em_perroad          = fill(NaN, nl)
    up.em_wall             = fill(NaN, nl)
    up.t_building_min      = fill(NaN, nl)
    up.thick_wall          = fill(NaN, nl)
    up.thick_roof          = fill(NaN, nl)
    up.nlev_improad        = fill(typemax(Int), nl)
    up.vf_sr               = fill(NaN, nl)
    up.vf_wr               = fill(NaN, nl)
    up.vf_sw               = fill(NaN, nl)
    up.vf_rw               = fill(NaN, nl)
    up.vf_ww               = fill(NaN, nl)
    up.eflx_traffic_factor = fill(NaN, nl)

    # 2D landunit × numrad (albedos)
    up.alb_roof_dir    = fill(NaN, nl, numrad)
    up.alb_roof_dif    = fill(NaN, nl, numrad)
    up.alb_improad_dir = fill(NaN, nl, numrad)
    up.alb_improad_dif = fill(NaN, nl, numrad)
    up.alb_perroad_dir = fill(NaN, nl, numrad)
    up.alb_perroad_dif = fill(NaN, nl, numrad)
    up.alb_wall_dir    = fill(NaN, nl, numrad)
    up.alb_wall_dif    = fill(NaN, nl, numrad)

    # 2D landunit × nlevurb (thermal properties)
    if nlevurb > 0
        up.tk_wall     = fill(NaN, nl, nlevurb)
        up.tk_roof     = fill(NaN, nl, nlevurb)
        up.cv_wall     = fill(NaN, nl, nlevurb)
        up.cv_roof     = fill(NaN, nl, nlevurb)
    else
        up.tk_wall     = Matrix{Float64}(undef, 0, 0)
        up.tk_roof     = Matrix{Float64}(undef, 0, 0)
        up.cv_wall     = Matrix{Float64}(undef, 0, 0)
        up.cv_roof     = Matrix{Float64}(undef, 0, 0)
    end
    up.tk_improad  = fill(NaN, nl, nlevurb)
    up.cv_improad  = fill(NaN, nl, nlevurb)

    return nothing
end

"""
    urbanparams_clean!(up::UrbanParamsData)

Deallocate (reset to empty) all fields of an `UrbanParamsData` instance.
"""
function urbanparams_clean!(up::UrbanParamsData{FT}) where {FT}
    up.wind_hgt_canyon     = FT[]
    up.em_roof             = FT[]
    up.em_improad          = FT[]
    up.em_perroad          = FT[]
    up.em_wall             = FT[]
    up.t_building_min      = FT[]
    up.thick_wall          = FT[]
    up.thick_roof          = FT[]
    up.nlev_improad        = Int[]
    up.vf_sr               = FT[]
    up.vf_wr               = FT[]
    up.vf_sw               = FT[]
    up.vf_rw               = FT[]
    up.vf_ww               = FT[]
    up.eflx_traffic_factor = FT[]

    up.alb_roof_dir    = Matrix{FT}(undef, 0, 0)
    up.alb_roof_dif    = Matrix{FT}(undef, 0, 0)
    up.alb_improad_dir = Matrix{FT}(undef, 0, 0)
    up.alb_improad_dif = Matrix{FT}(undef, 0, 0)
    up.alb_perroad_dir = Matrix{FT}(undef, 0, 0)
    up.alb_perroad_dif = Matrix{FT}(undef, 0, 0)
    up.alb_wall_dir    = Matrix{FT}(undef, 0, 0)
    up.alb_wall_dif    = Matrix{FT}(undef, 0, 0)

    up.tk_wall     = Matrix{FT}(undef, 0, 0)
    up.tk_roof     = Matrix{FT}(undef, 0, 0)
    up.tk_improad  = Matrix{FT}(undef, 0, 0)
    up.cv_wall     = Matrix{FT}(undef, 0, 0)
    up.cv_roof     = Matrix{FT}(undef, 0, 0)
    up.cv_improad  = Matrix{FT}(undef, 0, 0)

    return nothing
end

"""
    urbanparams_populate!(up::UrbanParamsData, lun::LandunitData,
                          urbinp::UrbanInputData, bounds_lun::UnitRange{Int};
                          urban_traffic::Bool=false,
                          use_vancouver::Bool=false,
                          use_mexicocity::Bool=false,
                          urban_hac::String=URBAN_HAC_OFF)

Populate urban parameters from input data for each urban landunit.
Computes view factors, aerodynamic constants, traffic factors,
and sets landunit-level urban properties (canyon_hwr, ht_roof, etc.).

Ported from the initialization loop in `urbanparams_type%Init` in `UrbanParamsType.F90`.
"""
function urbanparams_populate!(up::UrbanParamsData, lun::LandunitData,
                                urbinp::UrbanInputData, bounds_lun::UnitRange{Int};
                                urban_traffic::Bool = false,
                                use_vancouver::Bool = false,
                                use_mexicocity::Bool = false,
                                urban_hac::String = URBAN_HAC_OFF)
    numrad = NUMRAD
    nlevurb = size(up.tk_improad, 2)  # infer from allocated size

    # Constants for Macdonald (1998) aerodynamic calculation
    alpha = 4.43
    beta  = 1.0
    C_d   = 1.2

    for l in bounds_lun
        if lun.urbpoi[l]
            g = lun.gridcell[l]
            dindx = lun.itype[l] - ISTURB_MIN + 1

            # Copy emissivities and wind height from input data
            up.wind_hgt_canyon[l] = urbinp.wind_hgt_canyon[g, dindx]
            for ib in 1:numrad
                up.alb_roof_dir[l, ib]    = urbinp.alb_roof_dir[g, dindx, ib]
                up.alb_roof_dif[l, ib]    = urbinp.alb_roof_dif[g, dindx, ib]
                up.alb_improad_dir[l, ib] = urbinp.alb_improad_dir[g, dindx, ib]
                up.alb_perroad_dir[l, ib] = urbinp.alb_perroad_dir[g, dindx, ib]
                up.alb_improad_dif[l, ib] = urbinp.alb_improad_dif[g, dindx, ib]
                up.alb_perroad_dif[l, ib] = urbinp.alb_perroad_dif[g, dindx, ib]
                up.alb_wall_dir[l, ib]    = urbinp.alb_wall_dir[g, dindx, ib]
                up.alb_wall_dif[l, ib]    = urbinp.alb_wall_dif[g, dindx, ib]
            end
            up.em_roof[l]    = urbinp.em_roof[g, dindx]
            up.em_improad[l] = urbinp.em_improad[g, dindx]
            up.em_perroad[l] = urbinp.em_perroad[g, dindx]
            up.em_wall[l]    = urbinp.em_wall[g, dindx]

            # Set landunit-level urban properties
            lun.canyon_hwr[l]   = urbinp.canyon_hwr[g, dindx]
            lun.wtroad_perv[l]  = urbinp.wtroad_perv[g, dindx]
            lun.ht_roof[l]      = urbinp.ht_roof[g, dindx]
            lun.wtlunit_roof[l] = urbinp.wtlunit_roof[g, dindx]

            # Thermal properties
            if nlevurb > 0
                for k in 1:nlevurb
                    up.tk_wall[l, k]    = urbinp.tk_wall[g, dindx, k]
                    up.tk_roof[l, k]    = urbinp.tk_roof[g, dindx, k]
                    up.tk_improad[l, k] = urbinp.tk_improad[g, dindx, k]
                    up.cv_wall[l, k]    = urbinp.cv_wall[g, dindx, k]
                    up.cv_roof[l, k]    = urbinp.cv_roof[g, dindx, k]
                    up.cv_improad[l, k] = urbinp.cv_improad[g, dindx, k]
                end
            end
            up.thick_wall[l]     = urbinp.thick_wall[g, dindx]
            up.thick_roof[l]     = urbinp.thick_roof[g, dindx]
            up.nlev_improad[l]   = urbinp.nlev_improad[g, dindx]
            up.t_building_min[l] = urbinp.t_building_min[g, dindx]

            # Traffic factor — inferred from Sailor and Lu 2004
            if urban_traffic
                up.eflx_traffic_factor[l] = 3.6 * (lun.canyon_hwr[l] - 0.5) + 1.0
            else
                up.eflx_traffic_factor[l] = 0.0
            end

            # Building temperature minimum override
            if use_vancouver || use_mexicocity
                up.t_building_min[l] = 200.0
            else
                if urban_hac == URBAN_HAC_OFF
                    up.t_building_min[l] = 200.0
                end
            end

            # -----------------------------------------------------------------
            # View factors for road and one wall in urban canyon
            # Source: Masson (2000) Boundary-Layer Meteorology 94:357-397
            # -----------------------------------------------------------------

            # Road: sky view factor -> 1 as H/W -> 0, -> 0 as H/W -> inf
            up.vf_sr[l] = sqrt(lun.canyon_hwr[l]^2 + 1.0) - lun.canyon_hwr[l]
            up.vf_wr[l] = 0.5 * (1.0 - up.vf_sr[l])

            # Wall: sky view factor -> 0.5 as H/W -> 0, -> 0 as H/W -> inf
            up.vf_sw[l] = 0.5 * (lun.canyon_hwr[l] + 1.0 - sqrt(lun.canyon_hwr[l]^2 + 1.0)) / lun.canyon_hwr[l]
            up.vf_rw[l] = up.vf_sw[l]
            up.vf_ww[l] = 1.0 - up.vf_sw[l] - up.vf_rw[l]

            # Error check: view factors must sum to one
            sumvf = up.vf_sr[l] + 2.0 * up.vf_wr[l]
            if abs(sumvf - 1.0) > 1.0e-06
                error("urban road view factor error: sumvf=$sumvf for landunit $l")
            end
            sumvf = up.vf_sw[l] + up.vf_rw[l] + up.vf_ww[l]
            if abs(sumvf - 1.0) > 1.0e-06
                error("urban wall view factor error: sumvf=$sumvf for landunit $l")
            end

            # -----------------------------------------------------------------
            # Aerodynamic constants using Macdonald (1998) / Grimmond & Oke (1999)
            # -----------------------------------------------------------------

            # Plan area index
            plan_ai = lun.canyon_hwr[l] / (lun.canyon_hwr[l] + 1.0)

            # Building shape: shortside/longside ratio
            build_lw_ratio = plan_ai

            # Frontal area index
            frontal_ai = (1.0 - plan_ai) * lun.canyon_hwr[l]
            frontal_ai = frontal_ai * sqrt(1.0 / build_lw_ratio) * sqrt(plan_ai)

            # Displacement height
            if use_vancouver
                lun.z_d_town[l] = 3.5
            elseif use_mexicocity
                lun.z_d_town[l] = 10.9
            else
                lun.z_d_town[l] = (1.0 + alpha^(-plan_ai) * (plan_ai - 1.0)) * lun.ht_roof[l]
            end

            # Roughness length
            if use_vancouver
                lun.z_0_town[l] = 0.35
            elseif use_mexicocity
                lun.z_0_town[l] = 2.2
            else
                lun.z_0_town[l] = lun.ht_roof[l] * (1.0 - lun.z_d_town[l] / lun.ht_roof[l]) *
                    exp(-1.0 * (0.5 * beta * C_d / VKC^2 *
                    (1.0 - lun.z_d_town[l] / lun.ht_roof[l]) * frontal_ai)^(-0.5))
            end

        else
            # Non-urban landunit: set to SPVAL
            up.eflx_traffic_factor[l] = SPVAL
            up.t_building_min[l]      = SPVAL
            up.vf_sr[l] = SPVAL
            up.vf_wr[l] = SPVAL
            up.vf_sw[l] = SPVAL
            up.vf_rw[l] = SPVAL
            up.vf_ww[l] = SPVAL
        end
    end

    return nothing
end

"""
    check_urban(urbinp::UrbanInputData, pcturb::Matrix{Float64},
                bounds_grid::UnitRange{Int}; urban_valid::Union{BitVector,Nothing}=nothing,
                caller::String="")

Validate that all grid cells with pct urban > 0 have valid urban data.
Throws an error with diagnostic information if invalid data is found.

Ported from `CheckUrban` in `UrbanParamsType.F90`.
"""
function check_urban(urbinp::UrbanInputData, pcturb::Matrix{Float64},
                     bounds_grid::UnitRange{Int};
                     urban_valid::Union{BitVector,Nothing} = nothing,
                     caller::String = "")
    numurbl_check = size(pcturb, 2)

    for nl in bounds_grid
        for n in 1:numurbl_check
            if pcturb[nl, n] > 0.0
                valid = urban_valid === nothing ? true : urban_valid[nl]
                if !valid ||
                   urbinp.canyon_hwr[nl, n]      <= 0.0 ||
                   urbinp.em_improad[nl, n]      <= 0.0 ||
                   urbinp.em_perroad[nl, n]      <= 0.0 ||
                   urbinp.em_roof[nl, n]         <= 0.0 ||
                   urbinp.em_wall[nl, n]         <= 0.0 ||
                   urbinp.ht_roof[nl, n]         <= 0.0 ||
                   urbinp.thick_roof[nl, n]      <= 0.0 ||
                   urbinp.thick_wall[nl, n]      <= 0.0 ||
                   urbinp.t_building_min[nl, n]  <= 0.0 ||
                   urbinp.wind_hgt_canyon[nl, n]  <= 0.0 ||
                   urbinp.wtlunit_roof[nl, n]    <= 0.0 ||
                   urbinp.wtroad_perv[nl, n]     <= 0.0 ||
                   any(x -> x <= 0.0, urbinp.alb_improad_dir[nl, n, :]) ||
                   any(x -> x <= 0.0, urbinp.alb_improad_dif[nl, n, :]) ||
                   any(x -> x <= 0.0, urbinp.alb_perroad_dir[nl, n, :]) ||
                   any(x -> x <= 0.0, urbinp.alb_perroad_dif[nl, n, :]) ||
                   any(x -> x <= 0.0, urbinp.alb_roof_dir[nl, n, :]) ||
                   any(x -> x <= 0.0, urbinp.alb_roof_dif[nl, n, :]) ||
                   any(x -> x <= 0.0, urbinp.alb_wall_dir[nl, n, :]) ||
                   any(x -> x <= 0.0, urbinp.alb_wall_dif[nl, n, :]) ||
                   any(x -> x <= 0.0, urbinp.tk_roof[nl, n, :]) ||
                   any(x -> x <= 0.0, urbinp.tk_wall[nl, n, :]) ||
                   any(x -> x <= 0.0, urbinp.cv_roof[nl, n, :]) ||
                   any(x -> x <= 0.0, urbinp.cv_wall[nl, n, :])
                    error("$caller ERROR: no valid urban data for nl=$nl, density type=$n")
                else
                    if urbinp.nlev_improad[nl, n] > 0
                        nlev = urbinp.nlev_improad[nl, n]
                        if any(x -> x <= 0.0, urbinp.tk_improad[nl, n, 1:nlev]) ||
                           any(x -> x <= 0.0, urbinp.cv_improad[nl, n, 1:nlev])
                            error("$caller ERROR: no valid urban data for nl=$nl, density type=$n")
                        end
                    end
                end
            end
        end
    end

    return nothing
end

"""
    urban_read_nml!(ctrl::UrbanControl;
                    urban_hac::String=URBAN_HAC_OFF,
                    urban_explicit_ac::Bool=true,
                    urban_traffic::Bool=false,
                    building_temp_method::Int=BUILDING_TEMP_METHOD_PROG)

Set urban namelist parameters. In the Julia port, values are passed directly
rather than read from a namelist file.

Ported from `UrbanReadNML` in `UrbanParamsType.F90`.
"""
function urban_read_nml!(ctrl::UrbanControl;
                          urban_hac::String = URBAN_HAC_OFF,
                          urban_explicit_ac::Bool = true,
                          urban_traffic::Bool = false,
                          building_temp_method::Int = BUILDING_TEMP_METHOD_PROG)
    if urban_traffic
        error("Urban traffic fluxes are not implemented currently")
    end

    ctrl.urban_hac          = urban_hac
    ctrl.urban_explicit_ac  = urban_explicit_ac
    ctrl.urban_traffic      = urban_traffic
    ctrl.building_temp_method = building_temp_method
    ctrl.read_namelist      = true

    return nothing
end
