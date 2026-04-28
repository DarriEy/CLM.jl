# ==========================================================================
# Fortran CLM5 Restart Reader
#
# Reads a Fortran CLM5 restart file and injects state into CLM.jl instances.
# Handles dimension mismatches (Fortran has 1 col / 3 patches;
# CLM.jl may have 2 cols / 4 patches with different PFT mapping).
#
# Public function:
#   read_fortran_restart!(filepath, inst, bounds)
# ==========================================================================

using NCDatasets

"""
    read_fortran_restart!(filepath, inst, bounds)

Read a Fortran CLM5 restart file and inject all compatible state variables
into CLM.jl instances. Handles dimension mismatches gracefully.
"""
function read_fortran_restart!(filepath::String, inst::CLMInstances, bounds::BoundsType)
    isfile(filepath) || error("Restart file not found: $filepath")

    ds = NCDataset(filepath, "r")
    nlevsno = varpar.nlevsno
    nlevgrnd = varpar.nlevgrnd
    nlevsoi = varpar.nlevsoi

    nc_jl = bounds.endc  # CLM.jl column count
    np_jl = bounds.endp  # CLM.jl patch count

    # Fortran grid info
    nc_f90 = haskey(ds.dim, "column") ? ds.dim["column"] : 1
    np_f90 = haskey(ds.dim, "pft") ? ds.dim["pft"] : 3
    nlev_f90 = haskey(ds.dim, "levtot") ? ds.dim["levtot"] : 37

    n_ok = 0; n_skip = 0; n_fail = 0

    # Helper: inject 1D column variable
    function set_col_1d!(name, setter!)
        haskey(ds, name) || return
        raw = ds[name][:]; data = [ismissing(v) ? NaN : Float64(v) for v in raw]
        try
            # Map Fortran col 1 → CLM.jl col 1
            setter!(data[1])
            n_ok += 1
        catch e
            n_fail += 1
            @debug "Skip $name: $e"
        end
    end

    # Helper: inject 2D column variable (levtot × column)
    function set_col_2d!(name, target_array)
        haskey(ds, name) || return
        raw = ds[name][:]; data = [ismissing(v) ? NaN : Float64(v) for v in raw]
        try
            nlev = min(size(data, 1), size(target_array, 2))
            for j in 1:nlev
                target_array[1, j] = ndims(data) == 1 ? data[j] : data[j, 1]
            end
            n_ok += 1
        catch e
            n_fail += 1
            @debug "Skip $name: $e"
        end
    end

    # Helper: inject 1D patch variable (handling 3→4 patch mapping)
    function set_patch_1d!(name, target_array)
        haskey(ds, name) || return
        raw = ds[name][:]
        data = Float64[ismissing(v) ? NaN : Float64(v) for v in raw]
        try
            # Map Fortran patches to CLM.jl patches by PFT type
            f90_pfts = haskey(ds, "pfts1d_itypveg") ? Int.(ds["pfts1d_itypveg"][:]) : Int[]
            jl_pfts = inst.patch.itype

            for pf in 1:min(np_f90, length(data))
                # Find matching CLM.jl patch by PFT type
                for pj in 1:np_jl
                    if !isempty(f90_pfts) && f90_pfts[pf] == jl_pfts[pj]
                        target_array[pj] = data[pf]
                        break
                    elseif pf <= np_jl
                        target_array[pf] = data[pf]
                        break
                    end
                end
            end
            n_ok += 1
        catch e
            n_fail += 1
            @debug "Skip $name: $e"
        end
    end

    # =====================================================================
    # COLUMN-LEVEL STATE (most critical for parity)
    # =====================================================================

    # Temperature (37 levels: 12 snow + 25 soil)
    set_col_2d!("T_SOISNO", inst.temperature.t_soisno_col)
    set_col_1d!("T_GRND", v -> (inst.temperature.t_grnd_col[1] = v))
    set_col_1d!("T_GRND_R", v -> (inst.temperature.t_grnd_r_col[1] = v))

    # Water state (37 levels)
    ws = inst.water.waterstatebulk_inst.ws
    set_col_2d!("H2OSOI_LIQ", ws.h2osoi_liq_col)
    set_col_2d!("H2OSOI_ICE", ws.h2osoi_ice_col)
    set_col_1d!("WA", v -> (ws.wa_col[1] = v))
    set_col_1d!("H2OSFC", v -> (ws.h2osfc_col[1] = v))

    # Snow
    set_col_1d!("SNLSNO", v -> (inst.column.snl[1] = Int(v)))
    wd = inst.water.waterdiagnosticbulk_inst
    set_col_1d!("SNOW_DEPTH", v -> (wd.snow_depth_col[1] = v))
    set_col_1d!("frac_sno", v -> (wd.frac_sno_col[1] = v))
    set_col_1d!("frac_sno_eff", v -> (wd.frac_sno_eff_col[1] = v))

    # Snow layer geometry
    if haskey(ds, "DZSNO")
        dz = Float64.(ds["DZSNO"][:, 1])
        for j in 1:min(nlevsno, length(dz))
            inst.column.dz[1, j] = dz[j]
        end
        n_ok += 1
    end
    if haskey(ds, "ZSNO")
        z = Float64.(ds["ZSNO"][:, 1])
        for j in 1:min(nlevsno, length(z))
            inst.column.z[1, j] = z[j]
        end
        n_ok += 1
    end
    if haskey(ds, "ZISNO")
        zi = Float64.(ds["ZISNO"][:, 1])
        for j in 1:min(nlevsno, length(zi))
            inst.column.zi[1, j] = zi[j]
        end
        n_ok += 1
    end

    # Snow grain radii
    if haskey(ds, "snw_rds")
        sr = Float64.(ds["snw_rds"][:, 1])
        for j in 1:min(nlevsno, length(sr))
            wd.snw_rds_col[1, j] = sr[j]
        end
        n_ok += 1
    end

    # Snow aerosols
    for (vname, field) in [("mss_bcphi", :mss_bcphi_col), ("mss_bcpho", :mss_bcpho_col),
                           ("mss_ocphi", :mss_ocphi_col), ("mss_ocpho", :mss_ocpho_col),
                           ("mss_dst1", :mss_dst1_col), ("mss_dst2", :mss_dst2_col),
                           ("mss_dst3", :mss_dst3_col), ("mss_dst4", :mss_dst4_col)]
        if haskey(ds, vname) && hasproperty(inst.aerosol, field)
            data = Float64.(ds[vname][:, 1])
            arr = getfield(inst.aerosol, field)
            for j in 1:min(nlevsno, length(data), size(arr, 2))
                arr[1, j] = data[j]
            end
            n_ok += 1
        end
    end

    # Integrated snowfall
    wsb = inst.water.waterstatebulk_inst
    set_col_1d!("INT_SNOW", v -> (wsb.int_snow_col[1] = v))

    # Water table
    set_col_1d!("ZWT", v -> (inst.soilhydrology.zwt_col[1] = v))
    set_col_1d!("ZWT_PERCH", v -> (inst.soilhydrology.zwt_perched_col[1] = v))
    set_col_1d!("FROST_TABLE", v -> begin
        if hasproperty(inst.soilhydrology, :frost_table_col)
            inst.soilhydrology.frost_table_col[1] = v
        end
    end)

    # Albedo state
    if haskey(ds, "albgrd")
        data = Float64.(ds["albgrd"][:, 1])
        for b in 1:min(2, length(data))
            inst.surfalb.albgrd_col[1, b] = data[b]
        end
        n_ok += 1
    end
    if haskey(ds, "albgri")
        data = Float64.(ds["albgri"][:, 1])
        for b in 1:min(2, length(data))
            inst.surfalb.albgri_col[1, b] = data[b]
        end
        n_ok += 1
    end

    # Soil resistance
    set_col_1d!("SOILRESIS", v -> begin
        if hasproperty(inst.soilstate, :soilresis_col)
            inst.soilstate.soilresis_col[1] = v
        end
    end)

    # Ground roughness
    set_col_1d!("Z0MG", v -> (inst.frictionvel.z0mg_col[1] = v))

    # =====================================================================
    # PATCH-LEVEL STATE
    # =====================================================================

    # Vegetation temperature
    set_patch_1d!("T_VEG", inst.temperature.t_veg_patch)
    set_patch_1d!("T_REF2M", inst.temperature.t_ref2m_patch)

    # Canopy state
    set_patch_1d!("elai", inst.canopystate.elai_patch)
    set_patch_1d!("esai", inst.canopystate.esai_patch)
    set_patch_1d!("tlai", inst.canopystate.tlai_patch)
    set_patch_1d!("tsai", inst.canopystate.tsai_patch)
    set_patch_1d!("htop", inst.canopystate.htop_patch)
    set_patch_1d!("hbot", inst.canopystate.hbot_patch)
    set_patch_1d!("fsun", inst.canopystate.fsun_patch)

    # Canopy water
    set_patch_1d!("LIQCAN", ws.liqcan_patch)
    set_patch_1d!("SNOCAN", ws.snocan_patch)
    set_patch_1d!("FWET", wd.fwet_patch)
    set_patch_1d!("FCANSNO", wd.fcansno_patch)

    # Friction velocity
    set_patch_1d!("OBU", inst.frictionvel.obu_patch)

    # Compute H2OSNO from snow layer totals
    h2osno_total = 0.0
    for j in 1:nlevsno
        h2osno_total += ws.h2osoi_liq_col[1, j] + ws.h2osoi_ice_col[1, j]
    end
    ws.h2osno_no_layers_col[1] = h2osno_total

    close(ds)

    @info "Fortran restart loaded: $n_ok vars set, $n_fail failed, $n_skip skipped"
    return n_ok
end
