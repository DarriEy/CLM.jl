# ==========================================================================
# Ported from: src/main/clm_initializeMod.F90
# Master initialization orchestrator
# ==========================================================================

"""
    clm_initialize!(; fsurdat, paramfile, kwargs...)
    -> (inst, bounds, filt, tm)

Master initialization function. Reads surface data and parameter files,
builds the subgrid hierarchy, initializes vertical structure, sets initial
state, and returns all data structures ready for `clm_drv!()`.

# Required arguments
- `fsurdat::String`   — path to surface data NetCDF file (surfdata_clm.nc)
- `paramfile::String` — path to parameter NetCDF file (clm5_params.nc)

# Optional keyword arguments
- `start_date::DateTime`  — simulation start date (default 2000-01-01)
- `dtime::Int`            — timestep in seconds (default 1800)
- `use_cn::Bool`          — use CN biogeochemistry (default false)
- `use_crop::Bool`        — use crop model (default false)
- `use_bedrock::Bool`     — use bedrock (default true)
- `use_aquifer_layer::Bool` — use aquifer lower boundary (default true)
- `all_active::Bool`      — all points active (default false)
- `lat::Float64`          — gridcell latitude  (default: read from surfdata)
- `lon::Float64`          — gridcell longitude (default: read from surfdata)
- `soil_layerstruct::String` — soil layer structure (default "20SL_8.5m")
- `h2osfcflag::Int`       — surface water flag for soil hydrology (default 0)

# Returns
- `inst::CLMInstances`  — all physics/BGC data instances
- `bounds::BoundsType`  — subgrid bounds
- `filt::ClumpFilter`   — active/inactive filters
- `tm::TimeManager`     — time manager
"""
function clm_initialize!(;
    fsurdat::String,
    paramfile::String,
    start_date::DateTime = DateTime(2000, 1, 1),
    dtime::Int = 1800,
    use_cn::Bool = false,
    use_crop::Bool = false,
    use_bedrock::Bool = true,
    use_aquifer_layer::Bool = true,
    all_active::Bool = false,
    lat::Real = NaN,
    lon::Real = NaN,
    soil_layerstruct::String = "20SL_8.5m",
    h2osfcflag::Int = 0,
    fsnowoptics::String = "",
    fsnowaging::String = "",
    int_snow_max::Real = 2000.0)

    # ---- Step 1: Set control flags ----
    varctl.use_cn = use_cn
    varctl.use_crop = use_crop
    varctl.create_crop_landunit = use_crop
    varctl.use_bedrock = use_bedrock
    varctl.all_active = all_active
    varctl.soil_layerstruct_predefined = soil_layerstruct

    # ---- Step 2: Read surface file metadata ----
    (numpft, numcft) = surfrd_get_num_patches(fsurdat)
    nlevurb = surfrd_get_nlevurb(fsurdat)

    # ---- Step 3: Initialize varpar dimensions ----
    varpar_init!(varpar, 1, numpft, numcft, nlevurb)

    # ---- Step 4: Initialize varcon (vertical coordinate arrays) ----
    varcon_init!()

    # ---- Step 5: Read surface data ----
    ng = 1  # single gridcell
    surf = SurfaceInputData()
    surfrd_get_data!(surf, 1, ng, fsurdat)

    # ---- Step 6: Count subgrid elements ----
    (nl, nc, np) = count_subgrid_elements(surf, ng)

    # ---- Step 7: Build bounds ----
    bounds = BoundsType(
        begg=1, endg=ng,
        begl=1, endl=nl,
        begc=1, endc=nc,
        begp=1, endp=np,
        level=BOUNDS_LEVEL_CLUMP,
        clump_index=1
    )

    # ---- Step 8: Allocate all instances ----
    inst = CLMInstances()
    inst.surfdata = surf
    # Propagate CN flag to vegetation facade config before init
    inst.bgc_vegetation.config.use_cn = use_cn
    nlevdecomp_full = varpar.nlevdecomp_full
    ndecomp_cascade_transitions = use_cn ? 10 : 5
    clm_instInit!(inst; ng=ng, nl=nl, nc=nc, np=np,
                  nlevdecomp_full=nlevdecomp_full,
                  ndecomp_cascade_transitions=ndecomp_cascade_transitions)

    # ---- Step 9: Build subgrid hierarchy ----
    # Get lat/lon from surface file if not provided
    actual_lat = isnan(lat) ? _read_latlon(fsurdat, "LATIXY", 0.0) : lat
    actual_lon = isnan(lon) ? _read_latlon(fsurdat, "LONGXY", 0.0) : lon

    initGridCells!(bounds, surf,
                   inst.gridcell, inst.landunit, inst.column, inst.patch;
                   lat=actual_lat, lon=actual_lon)

    # ---- Step 10: Set active flags and check weights ----
    set_active!(bounds, inst.landunit, inst.column, inst.patch)
    check_weights!(bounds, inst.gridcell, inst.landunit, inst.column, inst.patch;
                   active_only=false)

    # ---- Step 11: Build filters ----
    alloc_all_filters!(nc, np, nl)
    set_filters!(bounds, inst.column, inst.landunit, inst.patch, inst.gridcell)
    filt = clump_filter

    # ---- Step 12: Read parameters ----
    readParameters!(paramfile)

    # ---- Step 13: Initialize vertical structure ----
    initVertical!(bounds, inst.gridcell, inst.landunit, inst.column, surf)

    # ---- Step 14: Cold-start initialization ----
    cold_start_initialize!(inst, bounds, filt, surf; use_aquifer_layer=use_aquifer_layer)

    # ---- Step 14a: Soil hydrology runtime controls ----
    soilhydrology_read_nl!(inst.soilhydrology; h2osfcflag=h2osfcflag)

    # ---- Step 15: Initialize snow layer constants ----
    _init_snow_layer_constants!()

    # ---- Step 15a: Initialize lake constants ----
    lake_con_init!()

    # ---- Step 15b: Initialize SNICAR optics and aging tables ----
    snow_optics_init!(snicar_optics; fsnowoptics=fsnowoptics)
    snowage_init!(snicar_aging; fsnowaging=fsnowaging)

    # ---- Step 15c: Initialize surface albedo constants ----
    mxsoil_color = 20  # CLM5 default
    surface_albedo_init_time_const!(inst.surfalb_con, mxsoil_color, surf.soil_color,
                                     inst.column.gridcell,
                                     bounds.begc:bounds.endc,
                                     bounds.begg:bounds.endg)

    # ---- Step 15d: Set control defaults for snow/soil physics ----
    if isempty(varctl.snow_thermal_cond_method)
        varctl.snow_thermal_cond_method = "Jordan1991"
    end
    if isempty(varctl.snow_cover_fraction_method)
        varctl.snow_cover_fraction_method = "SwensonLawrence2012"
    end

    # ---- Step 15e: Initialize snow cover fraction method ----
    scf = SnowCoverFractionSwensonLawrence2012()
    snow_cover_fraction_init!(scf, nc;
        col_lun_itype = inst.column.lun_itype,
        col_gridcell = inst.column.gridcell,
        col_topo_std = inst.column.topo_std,
        int_snow_max = int_snow_max)
    inst.scf_method = scf

    # ---- Step 15f: Initialize decomposition cascade for CN mode ----
    if use_cn
        # Read BGC decomposition parameters from param file
        NCDatasets.Dataset(paramfile) do ds
            ndecomp_pools_max = 7
            cstocks_raw = haskey(ds, "bgc_initial_Cstocks") ?
                Float64.(ds["bgc_initial_Cstocks"][:]) :
                [0.0, 0.0, 0.0, 200.0, 200.0, 200.0, 0.0]
            # Clamp fill values to 0
            cstocks = [v > 1e30 ? 0.0 : v for v in cstocks_raw[1:min(end, ndecomp_pools_max)]]
            while length(cstocks) < ndecomp_pools_max
                push!(cstocks, 0.0)
            end
            cstocks_depth = haskey(ds, "bgc_initial_Cstocks_depth") ?
                Float64(ds["bgc_initial_Cstocks_depth"][1]) : 0.3

            decomp_bgc_read_params!(inst.decomp_bgc_params;
                tau_l1=1.0/18.5, tau_l2_l3=1.0/4.9, tau_s1=1.0/7.3,
                tau_s2=1.0/0.2, tau_s3=1.0/0.0045,
                cn_s1=12.0, cn_s2=12.0, cn_s3=10.0,
                rf_l1s1=0.39, rf_l2s1=0.55, rf_l3s2=0.29,
                rf_s2s1=0.55, rf_s2s3=0.55, rf_s3s1=0.55,
                rf_cwdl3=0.0, cwd_fcel=0.0,
                bgc_initial_Cstocks=cstocks,
                bgc_initial_Cstocks_depth=cstocks_depth)
        end

        # Build cellsand matrix from soil state (sand fraction)
        nlevdecomp_val = varpar.nlevdecomp
        cellsand = zeros(nc, nlevdecomp_val)
        for c in 1:nc
            for j in 1:nlevdecomp_val
                cellsand[c, j] = inst.soilstate.cellsand_col[c, j]
            end
        end
        init_decomp_cascade_bgc!(inst.decomp_bgc_state,
                                  inst.decomp_cascade,
                                  inst.decomp_bgc_params,
                                  inst.cn_shared_params;
                                  cellsand=cellsand,
                                  bounds=1:nc,
                                  nlevdecomp=nlevdecomp_val,
                                  ndecomp_pools_max=7,
                                  ndecomp_cascade_transitions_max=10,
                                  spinup_state=0,
                                  use_fates=false)

        # Initialize competition state
        soil_bgc_competition_init!(inst.competition_state,
                                    inst.competition_params;
                                    dt=Float64(dtime))
    end

    # Initialize urban namelist if available
    if isdefined(CLM, :urban_read_nml!) && isdefined(CLM, :urban_ctrl)
        urban_read_nml!(urban_ctrl)
    end

    # ---- Step 16: Create time manager ----
    tm = TimeManager(
        start_date = start_date,
        current_date = start_date,
        dtime = dtime,
        nstep = 0,
        calendar = "NO_LEAP"
    )

    return (inst, bounds, filt, tm)
end

"""
    _init_snow_layer_constants!()

Initialize the snow layer thickness constants (SNOW_DZMIN, SNOW_DZMAX_U, SNOW_DZMAX_L).
"""
function _init_snow_layer_constants!()
    nlevsno = varpar.nlevsno
    dzmin = zeros(nlevsno)
    dzmax_u = zeros(nlevsno)
    dzmax_l = zeros(nlevsno)
    dzmin[1] = 0.010; dzmax_u[1] = 0.02; dzmax_l[1] = 0.03
    if nlevsno >= 2
        dzmin[2] = 0.015; dzmax_u[2] = 0.05; dzmax_l[2] = 0.07
    end
    for j in 3:nlevsno
        dzmin[j] = dzmax_u[j-1] * 0.5
        dzmax_u[j] = 2.0 * dzmax_u[j-1] + 0.01
        dzmax_l[j] = dzmax_u[j] + dzmax_l[j-1]
        if j == nlevsno
            dzmax_u[j] = floatmax(Float64)
            dzmax_l[j] = floatmax(Float64)
        end
    end
    SNOW_DZMIN[] = dzmin
    SNOW_DZMAX_U[] = dzmax_u
    SNOW_DZMAX_L[] = dzmax_l
    nothing
end

"""
    _read_latlon(fsurdat, varname, default) -> Float64

Read lat/lon from surface file, return scalar for single-gridcell case.
"""
function _read_latlon(fsurdat::String, varname::String, default::Real)
    ds = NCDataset(fsurdat, "r")
    try
        if haskey(ds, varname)
            data = Array(ds[varname])
            if data isa AbstractArray && length(data) > 0
                return Float64(data[1])
            end
        end
        return default
    finally
        close(ds)
    end
end
