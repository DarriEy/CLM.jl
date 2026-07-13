# ==========================================================================
# Ported from: src/biogeophys/TemperatureType.F90
# Temperature state data type allocation and initialization
# ==========================================================================

"""
    TemperatureData

Temperature state data structure. Holds all temperature variables at patch,
column, landunit, and gridcell levels, plus melt flags, emissivities,
heat content, and accumulated quantities.

Ported from `temperature_type` in `TemperatureType.F90`.
"""
Base.@kwdef mutable struct TemperatureData{FT<:Real,
                                           V<:AbstractVector{FT},
                                           M<:AbstractMatrix{FT},
                                           VI<:AbstractVector{<:Integer},
                                           MI<:AbstractMatrix{<:Integer}}
    # --- Temperatures (patch-level) ---
    t_stem_patch             ::V = Float64[]   # patch stem temperature (K)
    t_veg_patch              ::V = Float64[]   # patch vegetation temperature (K)
    t_skin_patch             ::V = Float64[]   # patch skin temperature (K)
    t_veg_day_patch          ::V = Float64[]   # patch daytime accumulative vegetation temperature (K*nsteps), LUNA
    t_veg_night_patch        ::V = Float64[]   # patch nighttime accumulative vegetation temperature (K*nsteps), LUNA
    t_veg10_day_patch        ::V = Float64[]   # 10-day running mean of patch daytime vegetation temperature (K), LUNA
    t_veg10_night_patch      ::V = Float64[]   # 10-day running mean of patch nighttime vegetation temperature (K), LUNA
    ndaysteps_patch          ::VI     = Int[]       # number of daytime steps from midnight, LUNA
    nnightsteps_patch        ::VI     = Int[]       # number of nighttime steps from midnight, LUNA
    dt_veg_patch             ::V = Float64[]   # patch change in t_veg, last iteration (K)
    thm_patch                ::V = Float64[]   # patch intermediate variable (forc_t+0.0098*forc_hgt_t_patch)

    # --- Temperatures (column-level, 1D) ---
    t_h2osfc_col             ::V = Float64[]   # col surface water temperature (K)
    t_h2osfc_bef_col         ::V = Float64[]   # col surface water temperature from time-step before (K)
    tsl_col                  ::V = Float64[]   # col temperature of near-surface soil layer (K)
    t_soi10cm_col            ::V = Float64[]   # col soil temperature in top 10cm (K)
    t_soi17cm_col            ::V = Float64[]   # col soil temperature in top 17cm (K)
    t_sno_mul_mss_col        ::V = Float64[]   # col snow temp * layer mass, layer sum (K*kg/m2)
    t_grnd_col               ::V = Float64[]   # col ground temperature (K)
    t_grnd_r_col             ::V = Float64[]   # col rural ground temperature (K)
    t_grnd_u_col             ::V = Float64[]   # col urban ground temperature (K)
    snot_top_col             ::V = Float64[]   # col temperature of top snow layer (K)
    dTdz_top_col             ::V = Float64[]   # col temperature gradient in top layer (K/m)
    dt_grnd_col              ::V = Float64[]   # col change in t_grnd, last iteration (K)
    thv_col                  ::V = Float64[]   # col virtual potential temperature (K)
    soila10_col              ::V = Float64[]   # col 10-day running mean of 12cm soil layer temperature (K)

    # --- Temperatures (column-level, 2D) ---
    t_ssbef_col              ::M = Matrix{Float64}(undef, 0, 0)  # col soil/snow temp before update (ncols, nlevsno+nlevmaxurbgrnd)
    t_soisno_col             ::M = Matrix{Float64}(undef, 0, 0)  # col soil/snow temperature (K) (ncols, nlevsno+nlevmaxurbgrnd)
    t_lake_col               ::M = Matrix{Float64}(undef, 0, 0)  # col lake temperature (K) (ncols, nlevlak)

    # --- Temperatures (landunit-level) ---
    t_building_lun           ::V = Float64[]   # lun internal building air temperature (K)
    t_roof_inner_lun         ::V = Float64[]   # lun roof inside surface temperature (K)
    t_sunw_inner_lun         ::V = Float64[]   # lun sunwall inside surface temperature (K)
    t_shdw_inner_lun         ::V = Float64[]   # lun shadewall inside surface temperature (K)
    t_floor_lun              ::V = Float64[]   # lun floor temperature (K)
    taf_lun                  ::V = Float64[]   # lun urban canopy air temperature (K)

    # --- 2m reference temperatures (patch-level) ---
    t_ref2m_patch            ::V = Float64[]   # patch 2m height surface air temperature (K)
    t_ref2m_r_patch          ::V = Float64[]   # patch rural 2m height surface air temperature (K)
    t_ref2m_u_patch          ::V = Float64[]   # patch urban 2m height surface air temperature (K)
    t_ref2m_min_patch        ::V = Float64[]   # patch daily minimum of average 2m height surface air temperature (K)
    t_ref2m_min_r_patch      ::V = Float64[]   # patch daily minimum of average 2m height surface air temperature - rural (K)
    t_ref2m_min_u_patch      ::V = Float64[]   # patch daily minimum of average 2m height surface air temperature - urban (K)
    t_ref2m_max_patch        ::V = Float64[]   # patch daily maximum of average 2m height surface air temperature (K)
    t_ref2m_max_r_patch      ::V = Float64[]   # patch daily maximum of average 2m height surface air temperature - rural (K)
    t_ref2m_max_u_patch      ::V = Float64[]   # patch daily maximum of average 2m height surface air temperature - urban (K)
    t_ref2m_min_inst_patch   ::V = Float64[]   # patch instantaneous daily min of average 2m height surface air temp (K)
    t_ref2m_min_inst_r_patch ::V = Float64[]   # patch instantaneous daily min - rural (K)
    t_ref2m_min_inst_u_patch ::V = Float64[]   # patch instantaneous daily min - urban (K)
    t_ref2m_max_inst_patch   ::V = Float64[]   # patch instantaneous daily max of average 2m height surface air temp (K)
    t_ref2m_max_inst_r_patch ::V = Float64[]   # patch instantaneous daily max - rural (K)
    t_ref2m_max_inst_u_patch ::V = Float64[]   # patch instantaneous daily max - urban (K)

    # TREFAV / TREFAV_U / TREFAV_R are hourly `timeavg` accumulators
    # (accum_period = nint(3600/dtime)); the daily min/max above track the min/max
    # of the HOURLY AVERAGES, not of the raw per-step t_ref2m. These scratch pairs
    # are the in-array equivalent of accumulMod's per-field val/nsteps buffers
    # (same pattern as BTRANAV in EnergyFluxType).
    trefav_accum_patch       ::V = Float64[]   # running sum of t_ref2m over the current hour
    trefav_naccum_patch      ::V = Float64[]   # steps accumulated in the current hour
    trefav_u_accum_patch     ::V = Float64[]   # ... urban
    trefav_u_naccum_patch    ::V = Float64[]
    trefav_r_accum_patch     ::V = Float64[]   # ... rural
    trefav_r_naccum_patch    ::V = Float64[]

    # --- Running mean temperatures (patch-level) ---
    t_a10_patch              ::V = Float64[]   # patch 10-day running mean of 2m temperature (K)
    t_a10min_patch           ::V = Float64[]   # patch 10-day running mean of min 2m temperature (K)
    t_a5min_patch            ::V = Float64[]   # patch 5-day running mean of min 2m temperature (K)

    # --- Accumulated quantities (patch-level) ---
    t_veg24_patch            ::V = Float64[]   # patch 24hr average vegetation temperature (K)
    t_veg240_patch           ::V = Float64[]   # patch 240hr average vegetation temperature (K)
    gdd0_patch               ::V = Float64[]   # patch growing degree-days base  0C from planting (ddays)
    gdd8_patch               ::V = Float64[]   # patch growing degree-days base  8C from planting (ddays)
    gdd10_patch              ::V = Float64[]   # patch growing degree-days base 10C from planting (ddays)
    gdd020_patch             ::V = Float64[]   # patch 20-year average of gdd0 (ddays)
    gdd820_patch             ::V = Float64[]   # patch 20-year average of gdd8 (ddays)
    gdd1020_patch            ::V = Float64[]   # patch 20-year average of gdd10 (ddays)

    # --- Heat content ---
    beta_col                 ::V = Float64[]   # col coefficient of convective velocity (-)
    dynbal_baseline_heat_col ::V = Float64[]   # col baseline heat content subtracted from total (J/m^2)
    heat1_grc                ::V = Float64[]   # grc initial gridcell total heat content (J/m^2)
    heat2_grc                ::V = Float64[]   # grc post land cover change total heat content (J/m^2)
    liquid_water_temp1_grc   ::V = Float64[]   # grc initial weighted average liquid water temperature (K)
    liquid_water_temp2_grc   ::V = Float64[]   # grc post land cover change weighted average liquid water temperature (K)

    # --- Flags ---
    imelt_col                ::MI     = Matrix{Int}(undef, 0, 0)  # flag for melting (=1), freezing (=2), Not=0 (ncols, nlevsno+nlevmaxurbgrnd)

    # --- Emissivities ---
    emv_patch                ::V = Float64[]   # patch vegetation emissivity
    emg_col                  ::V = Float64[]   # col ground emissivity

    # --- Misc ---
    xmf_col                  ::V = Float64[]   # col total latent heat of phase change of ground water
    xmf_h2osfc_col           ::V = Float64[]   # col latent heat of phase change of surface water
    fact_col                 ::M = Matrix{Float64}(undef, 0, 0)  # col used in computing tridiagonal matrix (ncols, nlevsno+nlevmaxurbgrnd)
    c_h2osfc_col             ::V = Float64[]   # col heat capacity of surface water

    # --- Namelist parameters ---
    excess_ice_coldstart_depth ::FT = 0.5   # depth below which excess ice will be present (m)
    excess_ice_coldstart_temp  ::FT = -5.0  # coldstart temperature of layers with excess ice (deg C)
end

# Convenience constructor preserving the TemperatureData{FT}() call sites (cold
# start, make_dual_copy): defaults the array type params to CPU Arrays. Adapt
# (below) rebuilds the struct with device array types (e.g. CuArray) on demand.
TemperatureData{FT}(; kwargs...) where {FT<:Real} =
    TemperatureData{FT, Vector{FT}, Matrix{FT}, Vector{Int}, Matrix{Int}}(; kwargs...)

Adapt.@adapt_structure TemperatureData

"""
    temperature_init!(temp::TemperatureData, np::Int, nc::Int, nl::Int, ng::Int)

Allocate and initialize all fields of a `TemperatureData` instance for
`np` patches, `nc` columns, `nl` landunits, and `ng` gridcells.
Real fields are initialized to `NaN` (or `SPVAL` where Fortran uses `spval`),
integer fields to `ISPVAL` or `typemax(Int)`.

Ported from `temperature_type%InitAllocate` in `TemperatureType.F90`.
"""
function temperature_init!(temp::TemperatureData{FT}, np::Int, nc::Int, nl::Int, ng::Int) where {FT}
    nlevsno        = varpar.nlevsno
    nlevmaxurbgrnd = varpar.nlevmaxurbgrnd
    nlevlak        = varpar.nlevlak

    # Vertical dimension size (Fortran offset -nlevsno+1:nlevmaxurbgrnd → nlevsno + nlevmaxurbgrnd levels)
    nlev_soisno = nlevsno + nlevmaxurbgrnd

    # --- Temperatures (patch-level) ---
    temp.t_stem_patch             = fill(FT(NaN), np)
    temp.t_veg_patch              = fill(FT(NaN), np)
    temp.t_skin_patch             = fill(FT(NaN), np)
    temp.t_veg_day_patch          = fill(FT(SPVAL), np)
    temp.t_veg_night_patch        = fill(FT(SPVAL), np)
    temp.t_veg10_day_patch        = fill(FT(SPVAL), np)
    temp.t_veg10_night_patch      = fill(FT(SPVAL), np)
    temp.ndaysteps_patch          = fill(ISPVAL, np)
    temp.nnightsteps_patch        = fill(ISPVAL, np)
    temp.dt_veg_patch             = fill(FT(NaN), np)
    temp.thm_patch                = fill(FT(NaN), np)

    # --- Temperatures (column-level, 1D) ---
    temp.t_h2osfc_col             = fill(FT(NaN), nc)
    temp.t_h2osfc_bef_col         = fill(FT(NaN), nc)
    temp.tsl_col                  = fill(FT(NaN), nc)
    temp.t_soi10cm_col            = fill(FT(NaN), nc)
    temp.t_soi17cm_col            = fill(FT(SPVAL), nc)
    temp.t_sno_mul_mss_col        = fill(FT(NaN), nc)
    temp.t_grnd_col               = fill(FT(NaN), nc)
    temp.t_grnd_r_col             = fill(FT(NaN), nc)
    temp.t_grnd_u_col             = fill(FT(NaN), nc)
    temp.snot_top_col             = fill(FT(NaN), nc)
    temp.dTdz_top_col             = fill(FT(NaN), nc)
    temp.dt_grnd_col              = fill(FT(NaN), nc)
    temp.thv_col                  = fill(FT(NaN), nc)
    temp.soila10_col              = fill(FT(NaN), nc)

    # --- Temperatures (column-level, 2D) ---
    temp.t_ssbef_col              = fill(FT(NaN), nc, nlev_soisno)
    temp.t_soisno_col             = fill(FT(NaN), nc, nlev_soisno)
    temp.t_lake_col               = fill(FT(NaN), nc, nlevlak)

    # --- Temperatures (landunit-level) ---
    temp.t_building_lun           = fill(FT(NaN), nl)
    temp.t_roof_inner_lun         = fill(FT(NaN), nl)
    temp.t_sunw_inner_lun         = fill(FT(NaN), nl)
    temp.t_shdw_inner_lun         = fill(FT(NaN), nl)
    temp.t_floor_lun              = fill(FT(NaN), nl)
    temp.taf_lun                  = fill(FT(NaN), nl)

    # --- 2m reference temperatures (patch-level) ---
    temp.t_ref2m_patch            = fill(FT(NaN), np)
    temp.t_ref2m_r_patch          = fill(FT(NaN), np)
    temp.t_ref2m_u_patch          = fill(FT(NaN), np)
    temp.t_ref2m_min_patch        = fill(FT(NaN), np)
    temp.t_ref2m_min_r_patch      = fill(FT(NaN), np)
    temp.t_ref2m_min_u_patch      = fill(FT(NaN), np)
    temp.t_ref2m_max_patch        = fill(FT(NaN), np)
    temp.t_ref2m_max_r_patch      = fill(FT(NaN), np)
    temp.t_ref2m_max_u_patch      = fill(FT(NaN), np)
    temp.t_ref2m_min_inst_patch   = fill(FT(NaN), np)
    temp.t_ref2m_min_inst_r_patch = fill(FT(NaN), np)
    temp.t_ref2m_min_inst_u_patch = fill(FT(NaN), np)
    temp.t_ref2m_max_inst_patch   = fill(FT(NaN), np)
    temp.t_ref2m_max_inst_r_patch = fill(FT(NaN), np)
    temp.t_ref2m_max_inst_u_patch = fill(FT(NaN), np)

    # TREFAV/_U/_R hourly-timeavg scratch (see the field declarations).
    temp.trefav_accum_patch       = fill(FT(0), np)
    temp.trefav_naccum_patch      = fill(FT(0), np)
    temp.trefav_u_accum_patch     = fill(FT(0), np)
    temp.trefav_u_naccum_patch    = fill(FT(0), np)
    temp.trefav_r_accum_patch     = fill(FT(0), np)
    temp.trefav_r_naccum_patch    = fill(FT(0), np)

    # --- Running mean temperatures (patch-level) ---
    temp.t_a10_patch              = fill(FT(NaN), np)
    temp.t_a10min_patch           = fill(FT(NaN), np)
    temp.t_a5min_patch            = fill(FT(NaN), np)

    # --- Accumulated quantities (patch-level) ---
    temp.t_veg24_patch            = fill(FT(NaN), np)
    temp.t_veg240_patch           = fill(FT(NaN), np)
    temp.gdd0_patch               = fill(FT(SPVAL), np)
    temp.gdd8_patch               = fill(FT(SPVAL), np)
    temp.gdd10_patch              = fill(FT(SPVAL), np)
    temp.gdd020_patch             = fill(FT(SPVAL), np)
    temp.gdd820_patch             = fill(FT(SPVAL), np)
    temp.gdd1020_patch            = fill(FT(SPVAL), np)

    # --- Heat content ---
    temp.beta_col                 = fill(FT(NaN), nc)
    temp.dynbal_baseline_heat_col = fill(FT(NaN), nc)
    temp.heat1_grc                = fill(FT(NaN), ng)
    temp.heat2_grc                = fill(FT(NaN), ng)
    temp.liquid_water_temp1_grc   = fill(FT(NaN), ng)
    temp.liquid_water_temp2_grc   = fill(FT(NaN), ng)

    # --- Flags ---
    temp.imelt_col                = fill(typemax(Int), nc, nlev_soisno)

    # --- Emissivities ---
    temp.emv_patch                = fill(FT(NaN), np)
    temp.emg_col                  = fill(FT(NaN), nc)

    # --- Misc ---
    temp.xmf_col                  = fill(FT(NaN), nc)
    temp.xmf_h2osfc_col           = fill(FT(NaN), nc)
    temp.fact_col                 = fill(FT(NaN), nc, nlev_soisno)
    temp.c_h2osfc_col             = fill(FT(NaN), nc)

    return nothing
end

"""
    temperature_clean!(temp::TemperatureData)

Deallocate (reset to empty) all fields of a `TemperatureData` instance.

Ported from `temperature_type%Clean` in `TemperatureType.F90`.
"""
function temperature_clean!(temp::TemperatureData{FT}) where {FT}
    # Patch-level
    temp.t_stem_patch             = FT[]
    temp.t_veg_patch              = FT[]
    temp.t_skin_patch             = FT[]
    temp.t_veg_day_patch          = FT[]
    temp.t_veg_night_patch        = FT[]
    temp.t_veg10_day_patch        = FT[]
    temp.t_veg10_night_patch      = FT[]
    temp.ndaysteps_patch          = Int[]
    temp.nnightsteps_patch        = Int[]
    temp.dt_veg_patch             = FT[]
    temp.thm_patch                = FT[]
    temp.t_ref2m_patch            = FT[]
    temp.t_ref2m_r_patch          = FT[]
    temp.t_ref2m_u_patch          = FT[]
    temp.t_ref2m_min_patch        = FT[]
    temp.t_ref2m_min_r_patch      = FT[]
    temp.t_ref2m_min_u_patch      = FT[]
    temp.t_ref2m_max_patch        = FT[]
    temp.t_ref2m_max_r_patch      = FT[]
    temp.t_ref2m_max_u_patch      = FT[]
    temp.t_ref2m_min_inst_patch   = FT[]
    temp.t_ref2m_min_inst_r_patch = FT[]
    temp.t_ref2m_min_inst_u_patch = FT[]
    temp.t_ref2m_max_inst_patch   = FT[]
    temp.t_ref2m_max_inst_r_patch = FT[]
    temp.t_ref2m_max_inst_u_patch = FT[]
    temp.trefav_accum_patch       = FT[]
    temp.trefav_naccum_patch      = FT[]
    temp.trefav_u_accum_patch     = FT[]
    temp.trefav_u_naccum_patch    = FT[]
    temp.trefav_r_accum_patch     = FT[]
    temp.trefav_r_naccum_patch    = FT[]
    temp.t_a10_patch              = FT[]
    temp.t_a10min_patch           = FT[]
    temp.t_a5min_patch            = FT[]
    temp.t_veg24_patch            = FT[]
    temp.t_veg240_patch           = FT[]
    temp.gdd0_patch               = FT[]
    temp.gdd8_patch               = FT[]
    temp.gdd10_patch              = FT[]
    temp.gdd020_patch             = FT[]
    temp.gdd820_patch             = FT[]
    temp.gdd1020_patch            = FT[]
    temp.emv_patch                = FT[]

    # Column-level
    temp.t_h2osfc_col             = FT[]
    temp.t_h2osfc_bef_col         = FT[]
    temp.tsl_col                  = FT[]
    temp.t_soi10cm_col            = FT[]
    temp.t_soi17cm_col            = FT[]
    temp.t_sno_mul_mss_col        = FT[]
    temp.t_grnd_col               = FT[]
    temp.t_grnd_r_col             = FT[]
    temp.t_grnd_u_col             = FT[]
    temp.snot_top_col             = FT[]
    temp.dTdz_top_col             = FT[]
    temp.dt_grnd_col              = FT[]
    temp.thv_col                  = FT[]
    temp.soila10_col              = FT[]
    temp.t_ssbef_col              = Matrix{FT}(undef, 0, 0)
    temp.t_soisno_col             = Matrix{FT}(undef, 0, 0)
    temp.t_lake_col               = Matrix{FT}(undef, 0, 0)
    temp.beta_col                 = FT[]
    temp.dynbal_baseline_heat_col = FT[]
    temp.emg_col                  = FT[]
    temp.xmf_col                  = FT[]
    temp.xmf_h2osfc_col           = FT[]
    temp.fact_col                 = Matrix{FT}(undef, 0, 0)
    temp.c_h2osfc_col             = FT[]
    temp.imelt_col                = Matrix{Int}(undef, 0, 0)

    # Landunit-level
    temp.t_building_lun           = FT[]
    temp.t_roof_inner_lun         = FT[]
    temp.t_sunw_inner_lun         = FT[]
    temp.t_shdw_inner_lun         = FT[]
    temp.t_floor_lun              = FT[]
    temp.taf_lun                  = FT[]

    # Gridcell-level
    temp.heat1_grc                = FT[]
    temp.heat2_grc                = FT[]
    temp.liquid_water_temp1_grc   = FT[]
    temp.liquid_water_temp2_grc   = FT[]

    return nothing
end

"""
    temperature_init_cold!(temp, col, lun, patch_data, bounds_col, bounds_patch, bounds_lun;
                           em_roof_lun, em_wall_lun, em_improad_lun, em_perroad_lun,
                           is_prog_buildtemp=false)

Initialize cold-start conditions for temperature state variables.
Sets soil/snow temperatures, ground temperatures, vegetation temperatures,
reference temperatures, urban building temperatures, and ground emissivities.

Arguments:
- `temp`: TemperatureData instance
- `col`: ColumnData instance (provides snl, landunit, itype)
- `lun`: LandunitData instance (provides lakpoi, itype, urbpoi, ifspecial, coli, colf)
- `patch_data`: PatchData instance (provides column, landunit)
- `bounds_col`: UnitRange{Int} for columns (1:ncols)
- `bounds_patch`: UnitRange{Int} for patches (1:npatches)
- `bounds_lun`: UnitRange{Int} for landunits (1:nlandunits)
- `em_roof_lun`: emissivity of roof per landunit
- `em_wall_lun`: emissivity of wall per landunit
- `em_improad_lun`: emissivity of impervious road per landunit
- `em_perroad_lun`: emissivity of pervious road per landunit
- `is_prog_buildtemp`: whether prognostic building temperature is used

Ported from `temperature_type%InitCold` in `TemperatureType.F90`.
"""
function temperature_init_cold!(temp::TemperatureData,
                                col::ColumnData,
                                lun::LandunitData,
                                patch_data::PatchData,
                                bounds_col::UnitRange{Int},
                                bounds_patch::UnitRange{Int},
                                bounds_lun::UnitRange{Int};
                                em_roof_lun::Vector{<:Real},
                                em_wall_lun::Vector{<:Real},
                                em_improad_lun::Vector{<:Real},
                                em_perroad_lun::Vector{<:Real},
                                is_prog_buildtemp::Bool = false)

    nlevgrnd       = varpar.nlevgrnd
    nlevsno        = varpar.nlevsno
    nlevmaxurbgrnd = varpar.nlevmaxurbgrnd
    nlevurb        = varpar.nlevurb
    nlevlak        = varpar.nlevlak

    # Helper to convert Fortran layer index (-nlevsno+1:nlevmaxurbgrnd) to Julia 1-based index
    # Fortran j maps to Julia j + nlevsno
    joff = nlevsno

    # Set snow/soil temperature
    for c in bounds_col
        l = col.landunit[c]

        # Initialize all levels to SPVAL
        temp.t_soisno_col[c, :] .= SPVAL

        # Snow level temperatures — all land points
        if col.snl[c] < 0
            for j in (col.snl[c] + 1):0
                temp.t_soisno_col[c, j + joff] = 250.0
            end
        end

        # Below-snow temperatures — nonlake points (lake points set below)
        if !lun.lakpoi[l]
            if lun.itype[l] == ISTICE
                for j in 1:nlevgrnd
                    temp.t_soisno_col[c, j + joff] = 250.0
                end

            elseif lun.itype[l] == ISTWET
                for j in 1:nlevgrnd
                    temp.t_soisno_col[c, j + joff] = 277.0
                end

            elseif lun.urbpoi[l]
                if col.itype[c] == ICOL_ROAD_PERV || col.itype[c] == ICOL_ROAD_IMPERV
                    for j in 1:nlevgrnd
                        temp.t_soisno_col[c, j + joff] = 274.0
                    end
                elseif col.itype[c] == ICOL_SUNWALL || col.itype[c] == ICOL_SHADEWALL || col.itype[c] == ICOL_ROOF
                    for j in 1:nlevurb
                        temp.t_soisno_col[c, j + joff] = 292.0
                    end
                end

            else
                for j in 1:nlevgrnd
                    temp.t_soisno_col[c, j + joff] = 272.0
                end
            end
        end
    end

    # Initialize internal building temperature from roof temperature
    if is_prog_buildtemp
        for l in bounds_lun
            for c in lun.coli[l]:lun.colf[l]
                if col.itype[c] == ICOL_ROOF
                    temp.t_roof_inner_lun[l] = temp.t_soisno_col[c, nlevurb + joff]
                    temp.t_building_lun[l]   = temp.t_soisno_col[c, nlevurb + joff]
                    temp.t_floor_lun[l]      = temp.t_soisno_col[c, nlevurb + joff]
                elseif col.itype[c] == ICOL_SUNWALL
                    temp.t_sunw_inner_lun[l] = temp.t_soisno_col[c, nlevurb + joff]
                elseif col.itype[c] == ICOL_SHADEWALL
                    temp.t_shdw_inner_lun[l] = temp.t_soisno_col[c, nlevurb + joff]
                end
            end
        end
    end

    # Set ground temperatures
    for c in bounds_col
        l = col.landunit[c]
        if lun.lakpoi[l]
            temp.t_grnd_col[c] = 277.0
        else
            temp.t_grnd_col[c] = temp.t_soisno_col[c, col.snl[c] + 1 + joff]
        end
        temp.t_soi17cm_col[c] = temp.t_grnd_col[c]
    end

    # Lake columns: set lake temp and soil temp to ground temp
    for c in bounds_col
        l = col.landunit[c]
        if lun.lakpoi[l]
            for j in 1:nlevlak
                temp.t_lake_col[c, j] = temp.t_grnd_col[c]
            end
            for j in 1:nlevgrnd
                temp.t_soisno_col[c, j + joff] = temp.t_grnd_col[c]
            end
        end
    end

    # Set t_h2osfc_col
    for c in bounds_col
        temp.t_h2osfc_col[c] = 274.0
    end

    # Set t_veg, t_stem, t_ref2m, t_ref2m_u, t_ref2m_r
    for p in bounds_patch
        c = patch_data.column[p]
        l = patch_data.landunit[p]

        temp.t_veg_patch[p]  = 283.0
        temp.t_stem_patch[p] = temp.t_veg_patch[p]
        temp.t_ref2m_patch[p] = 283.0

        if lun.urbpoi[l]
            temp.t_ref2m_u_patch[p] = 283.0
        else
            if !lun.ifspecial[l]
                temp.t_ref2m_r_patch[p] = 283.0
            else
                temp.t_ref2m_r_patch[p] = SPVAL
            end
        end
    end

    # Set urban canopy air temperature
    for l in bounds_lun
        if lun.urbpoi[l]
            temp.taf_lun[l] = 283.0
        end
    end

    # Set ground emissivity from urban parameters
    for c in bounds_col
        l = col.landunit[c]
        if col.itype[c] == ICOL_ROOF
            temp.emg_col[c] = em_roof_lun[l]
        elseif col.itype[c] == ICOL_SUNWALL
            temp.emg_col[c] = em_wall_lun[l]
        elseif col.itype[c] == ICOL_SHADEWALL
            temp.emg_col[c] = em_wall_lun[l]
        elseif col.itype[c] == ICOL_ROAD_IMPERV
            temp.emg_col[c] = em_improad_lun[l]
        elseif col.itype[c] == ICOL_ROAD_PERV
            temp.emg_col[c] = em_perroad_lun[l]
        end
    end

    # Initialize dynbal_baseline_heat_col
    for c in bounds_col
        temp.dynbal_baseline_heat_col[c] = 0.0
    end

    return nothing
end

# ==========================================================================
# The following per-type `*_InitHistory` / `*_Restart` / `*_InitAccBuffer` /
# `*_UpdateAccVars` methods are NOT implemented — they are no-op stubs that
# document the Fortran interface.
#
# This is NOT a missing-infrastructure gap: history I/O
# (`src/infrastructure/history_io.jl`, `history_writer.jl`), restart I/O
# (`src/infrastructure/restart_io.jl`, `fortran_restart.jl`) and the
# accumulator (`src/infrastructure/accumul.jl`) are all ported and live. CLM.jl
# does not route them through per-type methods the way Fortran does: history
# fields and restart variables are declared in a CENTRAL registry that reads and
# writes the `CLMInstances` tree directly. These stubs are therefore a
# structural artifact of the port, not an unported capability.
# ==========================================================================

"""
    temperature_init_history!(temp, bounds_col, bounds_patch, bounds_lun, bounds_grc;
                              is_simple_buildtemp=false, is_prog_buildtemp=false)

Register temperature fields for history file output.

Ported from `temperature_type%InitHistory` in `TemperatureType.F90`.
Not implemented (no-op stub). History I/O IS ported
(`src/infrastructure/history_io.jl`); fields are registered in a central
registry rather than per-type methods.
"""
function temperature_init_history!(temp::TemperatureData,
                                   bounds_col::UnitRange{Int},
                                   bounds_patch::UnitRange{Int},
                                   bounds_lun::UnitRange{Int},
                                   bounds_grc::UnitRange{Int};
                                   is_simple_buildtemp::Bool = false,
                                   is_prog_buildtemp::Bool = false)
    # No-op: history fields are registered centrally (infrastructure/history_io.jl).
    # All fields that would be registered:
    #   TH2OSFC, TG_U, TLAKE, SNO_T, TSA, TSA_R, TREFMNAV, TREFMXAV,
    #   TSA_U, TREFMNAV_U, TREFMXAV_U, TSTEM, TV, TSKIN, TG, TG_R,
    #   TSOI, TSOI_ICE, TSOI_10CM, TSL, SNOTXMASS, T10, SOIL10, A5TMIN,
    #   A10TMIN, TBUILD, TROOF_INNER, TSUNW_INNER, TSHDW_INNER, TFLOOR,
    #   HEAT_CONTENT1, HEAT_CONTENT2, LIQUID_WATER_TEMP1, SNOTTOPL,
    #   SNOdTdzL, DT_VEG, EMV, EMG, BETA, TV24, TV240,
    #   GDD0, GDD8, GDD10, GDD020, GDD820, GDD1020, TVEGD10, TVEGN10
    return nothing
end

"""
    temperature_restart!(temp, bounds_col, bounds_patch, bounds_lun;
                         flag="read", is_simple_buildtemp=false, is_prog_buildtemp=false)

Read/write temperature state from/to restart file.

Ported from `temperature_type%Restart` in `TemperatureType.F90`.
Not implemented (no-op stub). Restart I/O IS ported
(`src/infrastructure/restart_io.jl`, `fortran_restart.jl`); restart variables
are declared in a central registry rather than per-type methods.
"""
function temperature_restart!(temp::TemperatureData,
                              bounds_col::UnitRange{Int},
                              bounds_patch::UnitRange{Int},
                              bounds_lun::UnitRange{Int};
                              flag::String = "read",
                              is_simple_buildtemp::Bool = false,
                              is_prog_buildtemp::Bool = false)
    # No-op: restart variables are declared centrally (infrastructure/restart_io.jl).
    # Variables that would be read/written:
    #   T_SOISNO, T_VEG, T_STEM, TH2OSFC, T_LAKE, T_GRND, T_GRND_R, T_GRND_U,
    #   T_REF2M, T_REF2M_R, T_REF2M_U, T_REF2M_MIN/MAX (incl. _R, _U),
    #   T_REF2M_MIN/MAX_INST (incl. _R, _U), taf, DYNBAL_BASELINE_HEAT,
    #   gdd1020, gdd820, gdd020, tvegd10, tvegd, tvegn10, tvegn, tair10,
    #   ndaysteps, nnightsteps, t_building, t_roof_inner, t_sunw_inner,
    #   t_shdw_inner, t_floor
    return nothing
end

"""
    temperature_init_acc_buffer!(mgr, npts; active=fill(true, npts),
                                 use_crop=false, step_size=1800)

Register the crop growing-degree-day accumulation fields on `mgr`
(`AccumManager`), following `temperature_type%InitAccBuffer`.

Only the crop-GDD fields are registered here, because the other temperature
accumulators (T_VEG24/T_VEG240/T10/SOIL10/TDM5/…) are ported as direct
in-array running means in [`temperature_update_acc_vars!`](@ref) (the
established pattern for the always-on temperature accumulators in this port)
rather than going through the generic `AccumManager`. The crop GDDs go
through `AccumManager` because they use the `runaccum` (GDD0/8/10) and 20-year
`runmean` (GDD020/820/1020) machinery with Jan-1 / end-of-year resets.

When `use_crop` is true this registers (all `subgrid_type="pft"`, `numlev=1`):
- `GDD0`, `GDD8`, `GDD10` — `runaccum`, `init_value=0`
- `GDD020`, `GDD820`, `GDD1020` — `runmean`, 20-year window (`accum_period=-20*365`)

`active` is the per-patch active mask, `step_size` the model timestep in
seconds (used to convert the negative day-based runmean period to timesteps).

Ported from `temperature_type%InitAccBuffer` in `TemperatureType.F90`.
"""
function temperature_init_acc_buffer!(mgr::AccumManager,
                                      npts::Int;
                                      active::Vector{Bool} = fill(true, npts),
                                      use_crop::Bool = false,
                                      step_size::Int = 1800)
    use_crop || return nothing

    # All GDD summations are relative to the planting date (Kucharik & Brye 2003).
    # runaccum: never auto-reset by period (Fortran accum_period = not_used = huge).
    for (name, base) in (("GDD0", 0), ("GDD8", 8), ("GDD10", 10))
        init_accum_field!(mgr;
            name = name, units = "K",
            desc = "growing degree-days base $(base)C from planting",
            accum_type = "runaccum", accum_period = typemax(Int),
            subgrid_type = "pft", numlev = 1, init_value = 0.0,
            active = active, npts = npts, step_size = step_size)
    end

    # 20-year running means (20*365 days).
    for (name, base) in (("GDD020", 0), ("GDD820", 8), ("GDD1020", 10))
        init_accum_field!(mgr;
            name = name, units = "K",
            desc = "20-year running mean of growing degree days base $(base)C from planting",
            accum_type = "runmean", accum_period = -20 * 365,
            subgrid_type = "pft", numlev = 1, init_value = 0.0,
            active = active, npts = npts, step_size = step_size)
    end

    return nothing
end

"""
    temperature_init_acc_vars!(temp, bounds_col, bounds_patch; is_startup=false)

Initialize accumulated variables from the accumulation buffer.
Called for both initial and restart runs.

Ported from `temperature_type%InitAccVars` in `TemperatureType.F90`.
Not implemented (no-op stub). The accumulator IS ported
(`src/infrastructure/accumul.jl`, `AccumManager`) and is used by the live
driver (crop GDD, `t_mo_min`); these particular fields are simply not
registered with it.
"""
function temperature_init_acc_vars!(temp::TemperatureData,
                                   bounds_col::UnitRange{Int},
                                   bounds_patch::UnitRange{Int};
                                   is_startup::Bool = false)
    # When fully ported, this will:
    # 1. Extract T_VEG24, T_VEG240, T10, SOIL10, TDM5, TDM10 from accum buffer
    # 2. On startup, initialize t_ref2m_max/min and inst variants
    # 3. Extract GDD0, GDD8, GDD10, GDD020, GDD820, GDD1020 (crop only)

    if is_startup
        for p in bounds_patch
            temp.t_ref2m_max_patch[p]        =  SPVAL
            temp.t_ref2m_max_r_patch[p]      =  SPVAL
            temp.t_ref2m_max_u_patch[p]      =  SPVAL
            temp.t_ref2m_min_patch[p]        =  SPVAL
            temp.t_ref2m_min_r_patch[p]      =  SPVAL
            temp.t_ref2m_min_u_patch[p]      =  SPVAL
            temp.t_ref2m_max_inst_patch[p]   = -SPVAL
            temp.t_ref2m_max_inst_r_patch[p] = -SPVAL
            temp.t_ref2m_max_inst_u_patch[p] = -SPVAL
            temp.t_ref2m_min_inst_patch[p]   =  SPVAL
            temp.t_ref2m_min_inst_r_patch[p] =  SPVAL
            temp.t_ref2m_min_inst_u_patch[p] =  SPVAL
        end
    end

    return nothing
end

"""
    temperature_update_acc_vars!(temp, bounds_col, bounds_patch, lun, patch_data;
                                 end_cd=false, secs=0, dtime=0, nstep=0,
                                 mgr=nothing, use_crop=false,
                                 month=1, day=1, jday=1, is_end_curr_year=false,
                                 latdeg=Float64[], gridcell=Int[], active=Bool[],
                                 itype=Int[], npcropmin=typemax(Int),
                                 gdd20_season_start=Float64[],
                                 gdd20_season_end=Float64[])

Update accumulated temperature variables each timestep.
Handles 24/240-step running means of vegetation temperature, the 10-day
running mean of 2m temperature, and — when `mgr` (an `AccumManager`) is
supplied and `use_crop` is true — the crop growing-degree-day accumulation
(GDD0/GDD8/GDD10) plus the end-of-year 20-year running means
(GDD020/GDD820/GDD1020).

The always-on running means (T_VEG24/T_VEG240/T10) are done with the in-array
kernel `_temp_acc_kernel!` (the established temperature-accumulator pattern in
this port). The crop GDDs are routed through the generic `AccumManager` via
[`temperature_update_acc_vars_crop_gdds!`](@ref); the GDD fields must already be
registered with [`temperature_init_acc_buffer!`](@ref).

End-of-year GDD20 runmeans (mirroring Fortran's `is_end_curr_year` block):
when `is_end_curr_year`, if `varctl.flush_gdd20` is set the GDD020/820/1020
accumulators are flushed (reset) once, then each year's GDD0/8/10 is folded
into the 20-year running means.

Ported from `temperature_type%UpdateAccVars` in `TemperatureType.F90`.
"""
function temperature_update_acc_vars!(temp::TemperatureData,
                                     bounds_col::UnitRange{Int},
                                     bounds_patch::UnitRange{Int},
                                     lun::LandunitData,
                                     patch_data::PatchData;
                                     end_cd::Bool = false,
                                     secs::Int = 0,
                                     dtime::Int = 0,
                                     nstep::Int = 0,
                                     mgr::Union{AccumManager, Nothing} = nothing,
                                     use_crop::Bool = false,
                                     month::Int = 1,
                                     day::Int = 1,
                                     jday::Int = 1,
                                     is_end_curr_year::Bool = false,
                                     latdeg::AbstractVector{<:Real} = Float64[],
                                     gridcell::AbstractVector{Int} = Int[],
                                     active::AbstractVector{Bool} = Bool[],
                                     itype::AbstractVector{Int} = Int[],
                                     npcropmin::Int = typemax(Int),
                                     upper_soil_layer::Int = 0,
                                     gdd20_season_start::AbstractVector{<:Real} = Float64[],
                                     gdd20_season_end::AbstractVector{<:Real} = Float64[])
    # Running-mean window lengths, in TIMESTEPS, derived from the accumulation
    # PERIOD IN DAYS ÷ dtime — exactly as CTSM accumulMod converts a negative
    # accum_period (in days) via `-period * SHR_CONST_CDAY / get_step_size()`
    # (accumulMod.F90:246-247). CTSM periods: T_VEG24 = -1 day, T_VEG240 = -10
    # days, T10 = -10 days. These MUST scale with dtime — hardcoding them to the
    # dtime=1800s values (48/480/480) mis-sizes every window at dtime≠1800. A too-
    # long T10 window (e.g. 480 steps = 20 days at dtime=3600) lags the seasonal
    # 2-m temperature low, over-cooling the vcmax high-T acclimation term
    # (vcmaxse = 668.39 - 1.07·t10) → excessive high-T inhibition of vcmax at hot
    # midday → low assimilation → low gs → low transpiration.
    _dt = dtime > 0 ? dtime : 1800
    _win_1day  = accum_window_steps(1,  _dt)
    _win_10day = accum_window_steps(10, _dt)
    _win_5day  = accum_window_steps(5,  _dt)

    if !isempty(bounds_patch)
        _launch!(_temp_acc_kernel!, temp.t_veg24_patch, temp.t_veg240_patch,
            temp.t_a10_patch, temp.t_veg_patch, temp.t_ref2m_patch,
            nstep, _win_1day, _win_10day, first(bounds_patch), last(bounds_patch);
            ndrange = length(temp.t_veg24_patch))

        # --- TREFAV / TREFAV_U / TREFAV_R -> daily min/max of 2-m temperature ---
        # Fortran averages t_ref2m over each HOUR (the TREFAV `timeavg` accumulator,
        # accum_period = nint(3600/dtime)) and tracks the daily min/max of those
        # HOURLY AVERAGES. These four fields (and their _r/_u splits) were allocated
        # NaN and NEVER written — the accumulator was simply absent — so every
        # consumer read NaN: CNPhenology's `vernalization!` (t_ref2m_max/min drive
        # `cumvd` -> the winter-cereal vernalization factor `vf`), TDM5/TDM10 below,
        # and the TREFMXAV/TREFMNAV history fields.
        _hour_period = max(1, round(Int, 3600 / _dt, RoundNearestTiesAway))
        _hour_end    = (nstep % _hour_period == 0)
        _first_step_of_day = (secs == dtime)
        if length(temp.trefav_accum_patch) >= last(bounds_patch)
            _launch!(_trefav_kernel!,
                temp.t_ref2m_max_patch, temp.t_ref2m_min_patch,
                temp.t_ref2m_max_inst_patch, temp.t_ref2m_min_inst_patch,
                temp.trefav_accum_patch, temp.trefav_naccum_patch,
                temp.t_ref2m_patch, _hour_end, end_cd, _first_step_of_day,
                first(bounds_patch), last(bounds_patch);
                ndrange = length(temp.t_ref2m_patch))

            # Urban split: the end-of-day copy happens only on urban patches
            # (Fortran `if (lun%urbpoi(l))`); elsewhere max/min_u stay at their
            # missing-value flag. Same for rural with `.not. lun%ifspecial(l)`.
            _launch!(_trefav_masked_kernel!,
                temp.t_ref2m_max_u_patch, temp.t_ref2m_min_u_patch,
                temp.t_ref2m_max_inst_u_patch, temp.t_ref2m_min_inst_u_patch,
                temp.trefav_u_accum_patch, temp.trefav_u_naccum_patch,
                temp.t_ref2m_u_patch, patch_data.landunit, lun.urbpoi, true,
                _hour_end, end_cd, _first_step_of_day,
                first(bounds_patch), last(bounds_patch);
                ndrange = length(temp.t_ref2m_patch))

            _launch!(_trefav_masked_kernel!,
                temp.t_ref2m_max_r_patch, temp.t_ref2m_min_r_patch,
                temp.t_ref2m_max_inst_r_patch, temp.t_ref2m_min_inst_r_patch,
                temp.trefav_r_accum_patch, temp.trefav_r_naccum_patch,
                temp.t_ref2m_r_patch, patch_data.landunit, lun.ifspecial, false,
                _hour_end, end_cd, _first_step_of_day,
                first(bounds_patch), last(bounds_patch);
                ndrange = length(temp.t_ref2m_patch))

            # --- TDM5 (always) / TDM10 (use_crop) ---
            # Fortran: rbufslp = min(t_ref2m_min, t_ref2m_min_inst), and if that is
            # > 1e30 (i.e. both are still the spval flag, early in the run) it falls
            # back to TFRZ. Then a 5-day (TDM5) / 10-day (TDM10) running mean.
            # t_a5min feeds CNPhenology's seasonal-deciduous onset test; it was NaN.
            _launch!(_tdm_kernel!, temp.t_a5min_patch, temp.t_a10min_patch,
                temp.t_ref2m_min_patch, temp.t_ref2m_min_inst_patch,
                nstep, _win_5day, _win_10day, use_crop,
                first(bounds_patch), last(bounds_patch);
                ndrange = length(temp.t_a5min_patch))
        end
    end

    # --- SOIL10: 10-day running mean of the `upper_soil_layer` soil temperature ---
    # Fortran reads t_soisno_col(c, upper_soil_layer) where upper_soil_layer comes
    # from CNSharedParams (the layer containing 0.12 m). t_soisno is snow-padded in
    # this port, so the soil layer lives at varpar.nlevsno + upper_soil_layer. soila10_col
    # feeds CNPhenology's onset criteria and was NaN for the whole run.
    if upper_soil_layer > 0 && !isempty(bounds_col) &&
            length(temp.soila10_col) >= last(bounds_col)
        _launch!(_soil10_kernel!, temp.soila10_col, temp.t_soisno_col,
            varpar.nlevsno + upper_soil_layer, nstep, _win_10day,
            first(bounds_col), last(bounds_col);
            ndrange = length(temp.soila10_col))
    end

    # --- Crop growing-degree-day accumulation (use_crop only) ---
    # Mirrors the `if ( use_crop )then` block of Fortran UpdateAccVars: GDD0,
    # GDD8, GDD10 each timestep, then the end-of-year 20-year running means.
    if use_crop && mgr !== nothing && !isempty(bounds_patch)
        begp = first(bounds_patch)
        endp = last(bounds_patch)

        # GDD0 / GDD8 / GDD10 (runaccum, reset on Jan 1).
        for (base, gddx) in ((0,  temp.gdd0_patch),
                             (8,  temp.gdd8_patch),
                             (10, temp.gdd10_patch))
            temperature_update_acc_vars_crop_gdds!(temp, mgr, gddx, base;
                begp = begp, endp = endp, month = month, day = day,
                secs = secs, dtime = dtime, nstep = nstep, jday = jday,
                latdeg = latdeg, gridcell = gridcell, active = active,
                itype = itype, npcropmin = npcropmin,
                gdd20_season_start = gdd20_season_start,
                gdd20_season_end = gdd20_season_end)
        end

        # 20-year running means, updated once at end of year.
        if is_end_curr_year
            if varctl.flush_gdd20
                @info "Flushing GDD20 variables"
                markreset_accum_field!(mgr, "GDD020")
                markreset_accum_field!(mgr, "GDD820")
                markreset_accum_field!(mgr, "GDD1020")
                varctl.flush_gdd20 = false
            end
            update_accum_field!(mgr, "GDD020", temp.gdd0_patch, nstep)
            extract_accum_field!(mgr, "GDD020", temp.gdd020_patch, nstep)
            update_accum_field!(mgr, "GDD820", temp.gdd8_patch, nstep)
            extract_accum_field!(mgr, "GDD820", temp.gdd820_patch, nstep)
            update_accum_field!(mgr, "GDD1020", temp.gdd10_patch, nstep)
            extract_accum_field!(mgr, "GDD1020", temp.gdd1020_patch, nstep)
        end
    end

    return nothing
end

# Per-patch running means of veg temperature (24/240 ts) and 2m T (480 ts).
# Own-index read-modify-write; one thread per patch.
@kernel function _temp_acc_kernel!(t_veg24_patch, t_veg240_patch, t_a10_patch,
        @Const(t_veg_patch), @Const(t_ref2m_patch),
        nstep::Int, win_1day::Int, win_10day::Int, pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax
        tv = t_veg_patch[p]
        if isfinite(tv)
            # T_VEG24: 1-day (win_1day timesteps) running mean of veg temperature
            if isnan(t_veg24_patch[p]) || nstep <= 1
                t_veg24_patch[p] = tv
            else
                n24 = min(nstep, win_1day)
                t_veg24_patch[p] += (tv - t_veg24_patch[p]) / n24
            end

            # T_VEG240: 10-day (win_10day timesteps) running mean of veg temperature
            if isnan(t_veg240_patch[p]) || nstep <= 1
                t_veg240_patch[p] = tv
            else
                n240 = min(nstep, win_10day)
                t_veg240_patch[p] += (tv - t_veg240_patch[p]) / n240
            end

            # T10: 10-day (win_10day timesteps) running mean of 2m reference temperature
            t2m = t_ref2m_patch[p]
            if isfinite(t2m)
                if isnan(t_a10_patch[p]) || nstep <= 1
                    t_a10_patch[p] = t2m
                else
                    n480 = min(nstep, win_10day)
                    t_a10_patch[p] += (t2m - t_a10_patch[p]) / n480
                end
            end
        end
    end
end

# TREFAV: hourly `timeavg` of t_ref2m, feeding the daily min/max. Mirrors the
# BTRANAV pattern in EnergyFluxType. `_first_step_of_day` reproduces Fortran's
# `else if (secs == dtime)` branch, which flags max/min as missing at the start of
# a day so that a partial day never reports a bogus extreme.
@kernel function _trefav_kernel!(t_ref2m_max_patch, t_ref2m_min_patch,
        t_ref2m_max_inst_patch, t_ref2m_min_inst_patch,
        trefav_accum_patch, trefav_naccum_patch, @Const(t_ref2m_patch),
        hour_end::Bool, end_cd::Bool, first_step_of_day::Bool,
        pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax
        T = eltype(t_ref2m_max_patch)
        t2m = t_ref2m_patch[p]
        if isfinite(t2m) && t2m != T(SPVAL)
            trefav_accum_patch[p]  += t2m
            trefav_naccum_patch[p] += one(T)
        end

        # Hour boundary: extract the hourly average and fold into the daily extremes.
        if hour_end
            n = trefav_naccum_patch[p]
            if n > zero(T)
                avg = trefav_accum_patch[p] / n
                cmax = t_ref2m_max_inst_patch[p]
                cmin = t_ref2m_min_inst_patch[p]
                t_ref2m_max_inst_patch[p] =
                    (!isfinite(cmax) || cmax == -T(SPVAL)) ? avg : max(avg, cmax)
                t_ref2m_min_inst_patch[p] =
                    (!isfinite(cmin) || cmin ==  T(SPVAL)) ? avg : min(avg, cmin)
            end
            trefav_accum_patch[p]  = zero(T)
            trefav_naccum_patch[p] = zero(T)
        end

        if end_cd
            t_ref2m_max_patch[p]      =  t_ref2m_max_inst_patch[p]
            t_ref2m_min_patch[p]      =  t_ref2m_min_inst_patch[p]
            t_ref2m_max_inst_patch[p] = -T(SPVAL)
            t_ref2m_min_inst_patch[p] =  T(SPVAL)
        elseif first_step_of_day
            t_ref2m_max_patch[p] = T(SPVAL)
            t_ref2m_min_patch[p] = T(SPVAL)
        end
    end
end

# TREFAV_U / TREFAV_R: same as above, but the end-of-day copy is gated on a
# landunit predicate — urban (`lun%urbpoi`) for _U, non-special (`.not.
# lun%ifspecial`) for _R. `want` selects which sense of the flag passes.
@kernel function _trefav_masked_kernel!(t_ref2m_max_x, t_ref2m_min_x,
        t_ref2m_max_inst_x, t_ref2m_min_inst_x,
        accum_x, naccum_x, @Const(t_ref2m_x),
        @Const(patch_landunit), @Const(lun_flag), want::Bool,
        hour_end::Bool, end_cd::Bool, first_step_of_day::Bool,
        pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax
        T = eltype(t_ref2m_max_x)
        t2m = t_ref2m_x[p]
        if isfinite(t2m) && t2m != T(SPVAL)
            accum_x[p]  += t2m
            naccum_x[p] += one(T)
        end

        if hour_end
            n = naccum_x[p]
            if n > zero(T)
                avg = accum_x[p] / n
                cmax = t_ref2m_max_inst_x[p]
                cmin = t_ref2m_min_inst_x[p]
                t_ref2m_max_inst_x[p] =
                    (!isfinite(cmax) || cmax == -T(SPVAL)) ? avg : max(avg, cmax)
                t_ref2m_min_inst_x[p] =
                    (!isfinite(cmin) || cmin ==  T(SPVAL)) ? avg : min(avg, cmin)
            end
            accum_x[p]  = zero(T)
            naccum_x[p] = zero(T)
        end

        l = patch_landunit[p]
        passes = (lun_flag[l] == want)
        if end_cd
            if passes
                t_ref2m_max_x[p]      =  t_ref2m_max_inst_x[p]
                t_ref2m_min_x[p]      =  t_ref2m_min_inst_x[p]
                t_ref2m_max_inst_x[p] = -T(SPVAL)
                t_ref2m_min_inst_x[p] =  T(SPVAL)
            end
        elseif first_step_of_day
            t_ref2m_max_x[p] = T(SPVAL)
            t_ref2m_min_x[p] = T(SPVAL)
        end
    end
end

# TDM5 (5-day) / TDM10 (10-day, use_crop): running means of the daily minimum
# 2-m temperature. Fortran takes min(t_ref2m_min, t_ref2m_min_inst) and falls
# back to TFRZ while both are still the (huge) spval flag.
@kernel function _tdm_kernel!(t_a5min_patch, t_a10min_patch,
        @Const(t_ref2m_min_patch), @Const(t_ref2m_min_inst_patch),
        nstep::Int, win_5day::Int, win_10day::Int, use_crop::Bool,
        pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax
        T = eltype(t_a5min_patch)
        a = t_ref2m_min_patch[p]
        b = t_ref2m_min_inst_patch[p]
        va = isfinite(a) ? a : T(SPVAL)
        vb = isfinite(b) ? b : T(SPVAL)
        v  = min(va, vb)
        if v > T(1.0e30)
            v = T(TFRZ)
        end
        t_a5min_patch[p] = accum_runmean(t_a5min_patch[p], v, nstep, win_5day)
        if use_crop
            t_a10min_patch[p] = accum_runmean(t_a10min_patch[p], v, nstep, win_10day)
        end
    end
end

# SOIL10: 10-day running mean of the soil temperature at `lev` (= NLEVSNO +
# upper_soil_layer, the snow-padded index of the layer containing 0.12 m).
@kernel function _soil10_kernel!(soila10_col, @Const(t_soisno_col),
        lev::Int, nstep::Int, win_10day::Int, cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax
        ts = t_soisno_col[c, lev]
        if isfinite(ts)
            soila10_col[c] = accum_runmean(soila10_col[c], ts, nstep, win_10day)
        end
    end
end

"""
    temperature_update_acc_vars_crop_gdds!(temp, mgr, gddx_patch, basetemp_int;
        begp, endp, month, day, secs, dtime, nstep, jday,
        latdeg, gridcell, active, itype, npcropmin=typemax(Int),
        gdd20_season_start=Float64[], gdd20_season_end=Float64[])

Accumulate and extract one crop growing-degree-day field (e.g. GDD0/GDD8/GDD10).

For each active patch this computes the daily GDD contribution

    max(0, min(max_accum, t_ref2m - (TFRZ + basetemp))) * dtime/SECSPDAY

added into a `runaccum` accumulation field named `"GDD<basetemp>"`, resetting the
accumulator on Jan 1 (first step of the day) and contributing zero outside the
GDD20 accumulation season. The accumulated value is then extracted into
`gddx_patch`. `max_accum` is 26 for `basetemp_int == 0` and 30 otherwise.

The GDD20 accumulation season is, by default, derived from latitude
(Apr–Sep in the NH, Oct–Mar in the SH). If valid per-patch read-in season
boundaries are supplied via `gdd20_season_start`/`gdd20_season_end` (values in
[1, 366]), those override the latitude fallback for crop patches
(`itype >= npcropmin`).

The accumulation field must already be registered on `mgr` as a `runaccum`
field (see `init_accum_field!`).

Ported from `temperature_type%UpdateAccVars_CropGDDs` in `TemperatureType.F90`.
"""
function temperature_update_acc_vars_crop_gdds!(temp::TemperatureData,
                                                mgr::AccumManager,
                                                gddx_patch::AbstractVector{Float64},
                                                basetemp_int::Int;
                                                begp::Int,
                                                endp::Int,
                                                month::Int,
                                                day::Int,
                                                secs::Int,
                                                dtime::Int,
                                                nstep::Int,
                                                jday::Int,
                                                latdeg::AbstractVector{<:Real},
                                                gridcell::AbstractVector{Int},
                                                active::AbstractVector{Bool},
                                                itype::AbstractVector{Int},
                                                npcropmin::Int = typemax(Int),
                                                gdd20_season_start::AbstractVector{<:Real} = Float64[],
                                                gdd20_season_end::AbstractVector{<:Real} = Float64[])

    basetemp_r8 = Float64(basetemp_int)

    # Get maximum daily accumulation
    # (SSR 2024-05-31: unsure why base 0 differs, but preserved as in Fortran)
    max_accum = basetemp_int == 0 ? 26.0 : 30.0

    # Field name, e.g. "GDD0", "GDD8", "GDD10"
    field_name = "GDD" * string(basetemp_int)

    # Are valid read-in GDD20 seasons available?
    # Fortran: any(starts(begp:endp) > 0.5) .and. any(starts(begp:endp) < 366.5)
    have_seasons = !isempty(gdd20_season_start) && !isempty(gdd20_season_end)
    stream_gdd20_seasons_tt = false
    if have_seasons
        any_gt = false
        any_lt = false
        for p in begp:endp
            s = gdd20_season_start[p]
            any_gt |= s > 0.5
            any_lt |= s < 366.5
        end
        stream_gdd20_seasons_tt = any_gt && any_lt
    end

    # Per-patch single-level buffer for the daily increment (Fortran rbufslp).
    rbufslp = zeros(Float64, endp - begp + 1)

    for p in begp:endp
        # Avoid unnecessary calculations over inactive points
        active[p] || continue

        # Is this patch in its gdd20 accumulation season?
        # First, latitude-based fallback.
        lat = Float64(latdeg[gridcell[p]])
        in_accumulation_season =
            ((month > 3 && month < 10) && lat >= 0.0) ||
            ((month > 9 || month < 4) && lat < 0.0)

        # Replace with read-in gdd20 accumulation season, if needed and valid.
        if stream_gdd20_seasons_tt && itype[p] >= npcropmin
            gdd20_start = Int(gdd20_season_start[p])
            gdd20_end   = Int(gdd20_season_end[p])
            if gdd20_start >= 1 && gdd20_end >= 1
                if gdd20_start > 366 || gdd20_end > 366
                    error("invalid gdd20 season! start: $gdd20_start  end: $gdd20_end")
                end
                in_accumulation_season =
                    _is_doy_in_interval(gdd20_start, gdd20_end, jday)
            end
        end

        kf = p - begp + 1
        if month == 1 && day == 1 && secs == dtime
            # Jan 1, first timestep of the day: reset the accumulator.
            markreset_accum_field!(mgr, field_name; kf = kf)
            rbufslp[kf] = 0.0
        elseif in_accumulation_season
            rbufslp[kf] = max(0.0, min(max_accum,
                temp.t_ref2m_patch[p] - (TFRZ + basetemp_r8))) * dtime / SECSPDAY
        else
            rbufslp[kf] = 0.0      # keeps gdd unchanged outside accumulation season
        end
    end

    # Save: accumulate the increment, then extract the running total.
    # Note: runaccum cannot reset AND accumulate in the same call, so on the
    # reset step rbufslp is zero and the increment is applied the next step.
    update_accum_field!(mgr, field_name, rbufslp, nstep)
    extract_accum_field!(mgr, field_name, gddx_patch, nstep)

    return nothing
end

"""
    temperature_read_namelist!(temp; excess_ice_coldstart_depth=0.5, excess_ice_coldstart_temp=-5.0)

Set namelist parameters for temperature initialization.

Ported from `temperature_type%ReadNL` in `TemperatureType.F90`.
In Julia, namelist values are passed directly instead of reading from a file.
"""
function temperature_read_namelist!(temp::TemperatureData;
                                   excess_ice_coldstart_depth::Real = 0.5,
                                   excess_ice_coldstart_temp::Real = -5.0)
    if excess_ice_coldstart_depth <= 0.0
        error("excess_ice_coldstart_depth must be positive, got $excess_ice_coldstart_depth")
    end
    if excess_ice_coldstart_temp >= 0.0
        error("excess_ice_coldstart_temp must be below freezing, got $excess_ice_coldstart_temp")
    end
    temp.excess_ice_coldstart_depth = excess_ice_coldstart_depth
    temp.excess_ice_coldstart_temp  = excess_ice_coldstart_temp
    return nothing
end
