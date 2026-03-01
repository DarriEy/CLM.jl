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
Base.@kwdef mutable struct TemperatureData
    # --- Temperatures (patch-level) ---
    t_stem_patch             ::Vector{Float64} = Float64[]   # patch stem temperature (K)
    t_veg_patch              ::Vector{Float64} = Float64[]   # patch vegetation temperature (K)
    t_skin_patch             ::Vector{Float64} = Float64[]   # patch skin temperature (K)
    t_veg_day_patch          ::Vector{Float64} = Float64[]   # patch daytime accumulative vegetation temperature (K*nsteps), LUNA
    t_veg_night_patch        ::Vector{Float64} = Float64[]   # patch nighttime accumulative vegetation temperature (K*nsteps), LUNA
    t_veg10_day_patch        ::Vector{Float64} = Float64[]   # 10-day running mean of patch daytime vegetation temperature (K), LUNA
    t_veg10_night_patch      ::Vector{Float64} = Float64[]   # 10-day running mean of patch nighttime vegetation temperature (K), LUNA
    ndaysteps_patch          ::Vector{Int}     = Int[]       # number of daytime steps from midnight, LUNA
    nnightsteps_patch        ::Vector{Int}     = Int[]       # number of nighttime steps from midnight, LUNA
    dt_veg_patch             ::Vector{Float64} = Float64[]   # patch change in t_veg, last iteration (K)
    thm_patch                ::Vector{Float64} = Float64[]   # patch intermediate variable (forc_t+0.0098*forc_hgt_t_patch)

    # --- Temperatures (column-level, 1D) ---
    t_h2osfc_col             ::Vector{Float64} = Float64[]   # col surface water temperature (K)
    t_h2osfc_bef_col         ::Vector{Float64} = Float64[]   # col surface water temperature from time-step before (K)
    tsl_col                  ::Vector{Float64} = Float64[]   # col temperature of near-surface soil layer (K)
    t_soi10cm_col            ::Vector{Float64} = Float64[]   # col soil temperature in top 10cm (K)
    t_soi17cm_col            ::Vector{Float64} = Float64[]   # col soil temperature in top 17cm (K)
    t_sno_mul_mss_col        ::Vector{Float64} = Float64[]   # col snow temp * layer mass, layer sum (K*kg/m2)
    t_grnd_col               ::Vector{Float64} = Float64[]   # col ground temperature (K)
    t_grnd_r_col             ::Vector{Float64} = Float64[]   # col rural ground temperature (K)
    t_grnd_u_col             ::Vector{Float64} = Float64[]   # col urban ground temperature (K)
    snot_top_col             ::Vector{Float64} = Float64[]   # col temperature of top snow layer (K)
    dTdz_top_col             ::Vector{Float64} = Float64[]   # col temperature gradient in top layer (K/m)
    dt_grnd_col              ::Vector{Float64} = Float64[]   # col change in t_grnd, last iteration (K)
    thv_col                  ::Vector{Float64} = Float64[]   # col virtual potential temperature (K)
    soila10_col              ::Vector{Float64} = Float64[]   # col 10-day running mean of 12cm soil layer temperature (K)

    # --- Temperatures (column-level, 2D) ---
    t_ssbef_col              ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col soil/snow temp before update (ncols, nlevsno+nlevmaxurbgrnd)
    t_soisno_col             ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col soil/snow temperature (K) (ncols, nlevsno+nlevmaxurbgrnd)
    t_lake_col               ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col lake temperature (K) (ncols, nlevlak)

    # --- Temperatures (landunit-level) ---
    t_building_lun           ::Vector{Float64} = Float64[]   # lun internal building air temperature (K)
    t_roof_inner_lun         ::Vector{Float64} = Float64[]   # lun roof inside surface temperature (K)
    t_sunw_inner_lun         ::Vector{Float64} = Float64[]   # lun sunwall inside surface temperature (K)
    t_shdw_inner_lun         ::Vector{Float64} = Float64[]   # lun shadewall inside surface temperature (K)
    t_floor_lun              ::Vector{Float64} = Float64[]   # lun floor temperature (K)
    taf_lun                  ::Vector{Float64} = Float64[]   # lun urban canopy air temperature (K)

    # --- 2m reference temperatures (patch-level) ---
    t_ref2m_patch            ::Vector{Float64} = Float64[]   # patch 2m height surface air temperature (K)
    t_ref2m_r_patch          ::Vector{Float64} = Float64[]   # patch rural 2m height surface air temperature (K)
    t_ref2m_u_patch          ::Vector{Float64} = Float64[]   # patch urban 2m height surface air temperature (K)
    t_ref2m_min_patch        ::Vector{Float64} = Float64[]   # patch daily minimum of average 2m height surface air temperature (K)
    t_ref2m_min_r_patch      ::Vector{Float64} = Float64[]   # patch daily minimum of average 2m height surface air temperature - rural (K)
    t_ref2m_min_u_patch      ::Vector{Float64} = Float64[]   # patch daily minimum of average 2m height surface air temperature - urban (K)
    t_ref2m_max_patch        ::Vector{Float64} = Float64[]   # patch daily maximum of average 2m height surface air temperature (K)
    t_ref2m_max_r_patch      ::Vector{Float64} = Float64[]   # patch daily maximum of average 2m height surface air temperature - rural (K)
    t_ref2m_max_u_patch      ::Vector{Float64} = Float64[]   # patch daily maximum of average 2m height surface air temperature - urban (K)
    t_ref2m_min_inst_patch   ::Vector{Float64} = Float64[]   # patch instantaneous daily min of average 2m height surface air temp (K)
    t_ref2m_min_inst_r_patch ::Vector{Float64} = Float64[]   # patch instantaneous daily min - rural (K)
    t_ref2m_min_inst_u_patch ::Vector{Float64} = Float64[]   # patch instantaneous daily min - urban (K)
    t_ref2m_max_inst_patch   ::Vector{Float64} = Float64[]   # patch instantaneous daily max of average 2m height surface air temp (K)
    t_ref2m_max_inst_r_patch ::Vector{Float64} = Float64[]   # patch instantaneous daily max - rural (K)
    t_ref2m_max_inst_u_patch ::Vector{Float64} = Float64[]   # patch instantaneous daily max - urban (K)

    # --- Running mean temperatures (patch-level) ---
    t_a10_patch              ::Vector{Float64} = Float64[]   # patch 10-day running mean of 2m temperature (K)
    t_a10min_patch           ::Vector{Float64} = Float64[]   # patch 10-day running mean of min 2m temperature (K)
    t_a5min_patch            ::Vector{Float64} = Float64[]   # patch 5-day running mean of min 2m temperature (K)

    # --- Accumulated quantities (patch-level) ---
    t_veg24_patch            ::Vector{Float64} = Float64[]   # patch 24hr average vegetation temperature (K)
    t_veg240_patch           ::Vector{Float64} = Float64[]   # patch 240hr average vegetation temperature (K)
    gdd0_patch               ::Vector{Float64} = Float64[]   # patch growing degree-days base  0C from planting (ddays)
    gdd8_patch               ::Vector{Float64} = Float64[]   # patch growing degree-days base  8C from planting (ddays)
    gdd10_patch              ::Vector{Float64} = Float64[]   # patch growing degree-days base 10C from planting (ddays)
    gdd020_patch             ::Vector{Float64} = Float64[]   # patch 20-year average of gdd0 (ddays)
    gdd820_patch             ::Vector{Float64} = Float64[]   # patch 20-year average of gdd8 (ddays)
    gdd1020_patch            ::Vector{Float64} = Float64[]   # patch 20-year average of gdd10 (ddays)

    # --- Heat content ---
    beta_col                 ::Vector{Float64} = Float64[]   # col coefficient of convective velocity (-)
    dynbal_baseline_heat_col ::Vector{Float64} = Float64[]   # col baseline heat content subtracted from total (J/m^2)
    heat1_grc                ::Vector{Float64} = Float64[]   # grc initial gridcell total heat content (J/m^2)
    heat2_grc                ::Vector{Float64} = Float64[]   # grc post land cover change total heat content (J/m^2)
    liquid_water_temp1_grc   ::Vector{Float64} = Float64[]   # grc initial weighted average liquid water temperature (K)
    liquid_water_temp2_grc   ::Vector{Float64} = Float64[]   # grc post land cover change weighted average liquid water temperature (K)

    # --- Flags ---
    imelt_col                ::Matrix{Int}     = Matrix{Int}(undef, 0, 0)  # flag for melting (=1), freezing (=2), Not=0 (ncols, nlevsno+nlevmaxurbgrnd)

    # --- Emissivities ---
    emv_patch                ::Vector{Float64} = Float64[]   # patch vegetation emissivity
    emg_col                  ::Vector{Float64} = Float64[]   # col ground emissivity

    # --- Misc ---
    xmf_col                  ::Vector{Float64} = Float64[]   # col total latent heat of phase change of ground water
    xmf_h2osfc_col           ::Vector{Float64} = Float64[]   # col latent heat of phase change of surface water
    fact_col                 ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col used in computing tridiagonal matrix (ncols, nlevsno+nlevmaxurbgrnd)
    c_h2osfc_col             ::Vector{Float64} = Float64[]   # col heat capacity of surface water

    # --- Namelist parameters ---
    excess_ice_coldstart_depth ::Float64 = 0.5   # depth below which excess ice will be present (m)
    excess_ice_coldstart_temp  ::Float64 = -5.0  # coldstart temperature of layers with excess ice (deg C)
end

"""
    temperature_init!(temp::TemperatureData, np::Int, nc::Int, nl::Int, ng::Int)

Allocate and initialize all fields of a `TemperatureData` instance for
`np` patches, `nc` columns, `nl` landunits, and `ng` gridcells.
Real fields are initialized to `NaN` (or `SPVAL` where Fortran uses `spval`),
integer fields to `ISPVAL` or `typemax(Int)`.

Ported from `temperature_type%InitAllocate` in `TemperatureType.F90`.
"""
function temperature_init!(temp::TemperatureData, np::Int, nc::Int, nl::Int, ng::Int)
    nlevsno        = varpar.nlevsno
    nlevmaxurbgrnd = varpar.nlevmaxurbgrnd
    nlevlak        = varpar.nlevlak

    # Vertical dimension size (Fortran offset -nlevsno+1:nlevmaxurbgrnd → nlevsno + nlevmaxurbgrnd levels)
    nlev_soisno = nlevsno + nlevmaxurbgrnd

    # --- Temperatures (patch-level) ---
    temp.t_stem_patch             = fill(NaN, np)
    temp.t_veg_patch              = fill(NaN, np)
    temp.t_skin_patch             = fill(NaN, np)
    temp.t_veg_day_patch          = fill(SPVAL, np)
    temp.t_veg_night_patch        = fill(SPVAL, np)
    temp.t_veg10_day_patch        = fill(SPVAL, np)
    temp.t_veg10_night_patch      = fill(SPVAL, np)
    temp.ndaysteps_patch          = fill(ISPVAL, np)
    temp.nnightsteps_patch        = fill(ISPVAL, np)
    temp.dt_veg_patch             = fill(NaN, np)
    temp.thm_patch                = fill(NaN, np)

    # --- Temperatures (column-level, 1D) ---
    temp.t_h2osfc_col             = fill(NaN, nc)
    temp.t_h2osfc_bef_col         = fill(NaN, nc)
    temp.tsl_col                  = fill(NaN, nc)
    temp.t_soi10cm_col            = fill(NaN, nc)
    temp.t_soi17cm_col            = fill(SPVAL, nc)
    temp.t_sno_mul_mss_col        = fill(NaN, nc)
    temp.t_grnd_col               = fill(NaN, nc)
    temp.t_grnd_r_col             = fill(NaN, nc)
    temp.t_grnd_u_col             = fill(NaN, nc)
    temp.snot_top_col             = fill(NaN, nc)
    temp.dTdz_top_col             = fill(NaN, nc)
    temp.dt_grnd_col              = fill(NaN, nc)
    temp.thv_col                  = fill(NaN, nc)
    temp.soila10_col              = fill(NaN, nc)

    # --- Temperatures (column-level, 2D) ---
    temp.t_ssbef_col              = fill(NaN, nc, nlev_soisno)
    temp.t_soisno_col             = fill(NaN, nc, nlev_soisno)
    temp.t_lake_col               = fill(NaN, nc, nlevlak)

    # --- Temperatures (landunit-level) ---
    temp.t_building_lun           = fill(NaN, nl)
    temp.t_roof_inner_lun         = fill(NaN, nl)
    temp.t_sunw_inner_lun         = fill(NaN, nl)
    temp.t_shdw_inner_lun         = fill(NaN, nl)
    temp.t_floor_lun              = fill(NaN, nl)
    temp.taf_lun                  = fill(NaN, nl)

    # --- 2m reference temperatures (patch-level) ---
    temp.t_ref2m_patch            = fill(NaN, np)
    temp.t_ref2m_r_patch          = fill(NaN, np)
    temp.t_ref2m_u_patch          = fill(NaN, np)
    temp.t_ref2m_min_patch        = fill(NaN, np)
    temp.t_ref2m_min_r_patch      = fill(NaN, np)
    temp.t_ref2m_min_u_patch      = fill(NaN, np)
    temp.t_ref2m_max_patch        = fill(NaN, np)
    temp.t_ref2m_max_r_patch      = fill(NaN, np)
    temp.t_ref2m_max_u_patch      = fill(NaN, np)
    temp.t_ref2m_min_inst_patch   = fill(NaN, np)
    temp.t_ref2m_min_inst_r_patch = fill(NaN, np)
    temp.t_ref2m_min_inst_u_patch = fill(NaN, np)
    temp.t_ref2m_max_inst_patch   = fill(NaN, np)
    temp.t_ref2m_max_inst_r_patch = fill(NaN, np)
    temp.t_ref2m_max_inst_u_patch = fill(NaN, np)

    # --- Running mean temperatures (patch-level) ---
    temp.t_a10_patch              = fill(NaN, np)
    temp.t_a10min_patch           = fill(NaN, np)
    temp.t_a5min_patch            = fill(NaN, np)

    # --- Accumulated quantities (patch-level) ---
    temp.t_veg24_patch            = fill(NaN, np)
    temp.t_veg240_patch           = fill(NaN, np)
    temp.gdd0_patch               = fill(SPVAL, np)
    temp.gdd8_patch               = fill(SPVAL, np)
    temp.gdd10_patch              = fill(SPVAL, np)
    temp.gdd020_patch             = fill(SPVAL, np)
    temp.gdd820_patch             = fill(SPVAL, np)
    temp.gdd1020_patch            = fill(SPVAL, np)

    # --- Heat content ---
    temp.beta_col                 = fill(NaN, nc)
    temp.dynbal_baseline_heat_col = fill(NaN, nc)
    temp.heat1_grc                = fill(NaN, ng)
    temp.heat2_grc                = fill(NaN, ng)
    temp.liquid_water_temp1_grc   = fill(NaN, ng)
    temp.liquid_water_temp2_grc   = fill(NaN, ng)

    # --- Flags ---
    temp.imelt_col                = fill(typemax(Int), nc, nlev_soisno)

    # --- Emissivities ---
    temp.emv_patch                = fill(NaN, np)
    temp.emg_col                  = fill(NaN, nc)

    # --- Misc ---
    temp.xmf_col                  = fill(NaN, nc)
    temp.xmf_h2osfc_col           = fill(NaN, nc)
    temp.fact_col                 = fill(NaN, nc, nlev_soisno)
    temp.c_h2osfc_col             = fill(NaN, nc)

    return nothing
end

"""
    temperature_clean!(temp::TemperatureData)

Deallocate (reset to empty) all fields of a `TemperatureData` instance.

Ported from `temperature_type%Clean` in `TemperatureType.F90`.
"""
function temperature_clean!(temp::TemperatureData)
    # Patch-level
    temp.t_stem_patch             = Float64[]
    temp.t_veg_patch              = Float64[]
    temp.t_skin_patch             = Float64[]
    temp.t_veg_day_patch          = Float64[]
    temp.t_veg_night_patch        = Float64[]
    temp.t_veg10_day_patch        = Float64[]
    temp.t_veg10_night_patch      = Float64[]
    temp.ndaysteps_patch          = Int[]
    temp.nnightsteps_patch        = Int[]
    temp.dt_veg_patch             = Float64[]
    temp.thm_patch                = Float64[]
    temp.t_ref2m_patch            = Float64[]
    temp.t_ref2m_r_patch          = Float64[]
    temp.t_ref2m_u_patch          = Float64[]
    temp.t_ref2m_min_patch        = Float64[]
    temp.t_ref2m_min_r_patch      = Float64[]
    temp.t_ref2m_min_u_patch      = Float64[]
    temp.t_ref2m_max_patch        = Float64[]
    temp.t_ref2m_max_r_patch      = Float64[]
    temp.t_ref2m_max_u_patch      = Float64[]
    temp.t_ref2m_min_inst_patch   = Float64[]
    temp.t_ref2m_min_inst_r_patch = Float64[]
    temp.t_ref2m_min_inst_u_patch = Float64[]
    temp.t_ref2m_max_inst_patch   = Float64[]
    temp.t_ref2m_max_inst_r_patch = Float64[]
    temp.t_ref2m_max_inst_u_patch = Float64[]
    temp.t_a10_patch              = Float64[]
    temp.t_a10min_patch           = Float64[]
    temp.t_a5min_patch            = Float64[]
    temp.t_veg24_patch            = Float64[]
    temp.t_veg240_patch           = Float64[]
    temp.gdd0_patch               = Float64[]
    temp.gdd8_patch               = Float64[]
    temp.gdd10_patch              = Float64[]
    temp.gdd020_patch             = Float64[]
    temp.gdd820_patch             = Float64[]
    temp.gdd1020_patch            = Float64[]
    temp.emv_patch                = Float64[]

    # Column-level
    temp.t_h2osfc_col             = Float64[]
    temp.t_h2osfc_bef_col         = Float64[]
    temp.tsl_col                  = Float64[]
    temp.t_soi10cm_col            = Float64[]
    temp.t_soi17cm_col            = Float64[]
    temp.t_sno_mul_mss_col        = Float64[]
    temp.t_grnd_col               = Float64[]
    temp.t_grnd_r_col             = Float64[]
    temp.t_grnd_u_col             = Float64[]
    temp.snot_top_col             = Float64[]
    temp.dTdz_top_col             = Float64[]
    temp.dt_grnd_col              = Float64[]
    temp.thv_col                  = Float64[]
    temp.soila10_col              = Float64[]
    temp.t_ssbef_col              = Matrix{Float64}(undef, 0, 0)
    temp.t_soisno_col             = Matrix{Float64}(undef, 0, 0)
    temp.t_lake_col               = Matrix{Float64}(undef, 0, 0)
    temp.beta_col                 = Float64[]
    temp.dynbal_baseline_heat_col = Float64[]
    temp.emg_col                  = Float64[]
    temp.xmf_col                  = Float64[]
    temp.xmf_h2osfc_col           = Float64[]
    temp.fact_col                 = Matrix{Float64}(undef, 0, 0)
    temp.c_h2osfc_col             = Float64[]
    temp.imelt_col                = Matrix{Int}(undef, 0, 0)

    # Landunit-level
    temp.t_building_lun           = Float64[]
    temp.t_roof_inner_lun         = Float64[]
    temp.t_sunw_inner_lun         = Float64[]
    temp.t_shdw_inner_lun         = Float64[]
    temp.t_floor_lun              = Float64[]
    temp.taf_lun                  = Float64[]

    # Gridcell-level
    temp.heat1_grc                = Float64[]
    temp.heat2_grc                = Float64[]
    temp.liquid_water_temp1_grc   = Float64[]
    temp.liquid_water_temp2_grc   = Float64[]

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
                                em_roof_lun::Vector{Float64},
                                em_wall_lun::Vector{Float64},
                                em_improad_lun::Vector{Float64},
                                em_perroad_lun::Vector{Float64},
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
# The following subroutines depend on infrastructure modules that are not yet
# ported (history, restart/IO, accumulator). They are provided as stubs that
# document the Fortran interface and can be filled in when those modules
# become available.
# ==========================================================================

"""
    temperature_init_history!(temp, bounds_col, bounds_patch, bounds_lun, bounds_grc;
                              is_simple_buildtemp=false, is_prog_buildtemp=false)

Register temperature fields for history file output.

Ported from `temperature_type%InitHistory` in `TemperatureType.F90`.
Requires history infrastructure (histFileMod) — stub until that module is ported.
"""
function temperature_init_history!(temp::TemperatureData,
                                   bounds_col::UnitRange{Int},
                                   bounds_patch::UnitRange{Int},
                                   bounds_lun::UnitRange{Int},
                                   bounds_grc::UnitRange{Int};
                                   is_simple_buildtemp::Bool = false,
                                   is_prog_buildtemp::Bool = false)
    # Stub: history field registration will be added when histFileMod is ported.
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
Requires NetCDF/restart infrastructure — stub until that module is ported.
"""
function temperature_restart!(temp::TemperatureData,
                              bounds_col::UnitRange{Int},
                              bounds_patch::UnitRange{Int},
                              bounds_lun::UnitRange{Int};
                              flag::String = "read",
                              is_simple_buildtemp::Bool = false,
                              is_prog_buildtemp::Bool = false)
    # Stub: restart variable I/O will be added when restUtilMod/ncdio_pio is ported.
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
    temperature_init_acc_buffer!(temp, bounds_patch)

Initialize accumulation buffer for temperature-related accumulated fields.

Ported from `temperature_type%InitAccBuffer` in `TemperatureType.F90`.
Requires accumulation infrastructure (accumulMod) — stub until that module is ported.
"""
function temperature_init_acc_buffer!(temp::TemperatureData,
                                     bounds_patch::UnitRange{Int})
    # Stub: accumulation field definitions will be added when accumulMod is ported.
    # Fields that would be initialized:
    #   T_VEG24 (runmean, -1 day), T_VEG240 (runmean, -10 days),
    #   TREFAV (timeavg, 1 hour), TREFAV_U, TREFAV_R,
    #   T10 (runmean, -10 days), SOIL10 (runmean, -10 days), TDM5 (runmean, -5 days),
    #   TDM10 (runmean, -10 days, crop only),
    #   GDD0, GDD8, GDD10 (runaccum, crop only),
    #   GDD020, GDD820, GDD1020 (runmean, -20*365 days, crop only),
    #   TDA (timeavg, -30 days, cndv only)
    return nothing
end

"""
    temperature_init_acc_vars!(temp, bounds_col, bounds_patch; is_startup=false)

Initialize accumulated variables from the accumulation buffer.
Called for both initial and restart runs.

Ported from `temperature_type%InitAccVars` in `TemperatureType.F90`.
Requires accumulation infrastructure (accumulMod) — stub until that module is ported.
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
                                t_ref2m_patch, t_ref2m_u_patch, t_ref2m_r_patch,
                                end_cd, secs, dtime)

Update accumulated temperature variables each timestep.
Handles hourly averages of 2m temperature, daily min/max tracking,
10-day running means, and growing degree-day accumulation.

Ported from `temperature_type%UpdateAccVars` in `TemperatureType.F90`.
Core logic is ported; accumulator calls are stubs until accumulMod is ported.
"""
function temperature_update_acc_vars!(temp::TemperatureData,
                                     bounds_col::UnitRange{Int},
                                     bounds_patch::UnitRange{Int},
                                     lun::LandunitData,
                                     patch_data::PatchData;
                                     end_cd::Bool = false,
                                     secs::Int = 0,
                                     dtime::Int = 0)
    # --- T_VEG24 / T_VEG240 ---
    # Stub: update_accum_field/extract_accum_field calls go here

    # --- TREFAV: hourly 2m air temperature average ---
    # When accumulator is available, extract hourly average into rbufslp,
    # then update daily min/max tracking:
    #
    # for p in bounds_patch
    #     if rbufslp[p] != SPVAL
    #         temp.t_ref2m_max_inst_patch[p] = max(rbufslp[p], temp.t_ref2m_max_inst_patch[p])
    #         temp.t_ref2m_min_inst_patch[p] = min(rbufslp[p], temp.t_ref2m_min_inst_patch[p])
    #     end
    #     if end_cd
    #         temp.t_ref2m_max_patch[p] = temp.t_ref2m_max_inst_patch[p]
    #         temp.t_ref2m_min_patch[p] = temp.t_ref2m_min_inst_patch[p]
    #         temp.t_ref2m_max_inst_patch[p] = -SPVAL
    #         temp.t_ref2m_min_inst_patch[p] =  SPVAL
    #     elseif secs == dtime
    #         temp.t_ref2m_max_patch[p] = SPVAL
    #         temp.t_ref2m_min_patch[p] = SPVAL
    #     end
    # end

    # Similar logic for TREFAV_U (urban) and TREFAV_R (rural)
    # and for T10, SOIL10, TDM5, TDM10, GDD0/8/10, GDD020/820/1020

    return nothing
end

"""
    temperature_read_namelist!(temp; excess_ice_coldstart_depth=0.5, excess_ice_coldstart_temp=-5.0)

Set namelist parameters for temperature initialization.

Ported from `temperature_type%ReadNL` in `TemperatureType.F90`.
In Julia, namelist values are passed directly instead of reading from a file.
"""
function temperature_read_namelist!(temp::TemperatureData;
                                   excess_ice_coldstart_depth::Float64 = 0.5,
                                   excess_ice_coldstart_temp::Float64 = -5.0)
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
