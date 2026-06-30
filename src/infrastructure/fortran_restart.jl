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
function read_fortran_restart!(filepath::String, inst::CLMInstances, bounds::BoundsType;
                               inject_veg_struct::Bool = true)
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

    # Helper: inject a named per-pool vertical var. The Fortran restart writes the
    # decomposition pools as SEPARATE 2D vars (column, levgrnd) — litr1c_vr,
    # litr2c_vr, ... cwdc_vr — whereas CLM.jl stores one 3D array
    # decomp_*pools_vr_col[col, lev, pool]. This routes one named var into pool
    # slot `pidx`. NCDatasets returns (lev, column), so data[j,1] is level j.
    function set_col_pool_vr!(name, target3d, pidx)
        haskey(ds, name) || return
        (ndims(target3d) == 3 && size(target3d, 1) >= 1 &&
         size(target3d, 2) >= 1 && size(target3d, 3) >= pidx) || return
        raw = ds[name][:]; data = [ismissing(v) ? NaN : Float64(v) for v in raw]
        try
            nlev = min(size(data, 1), size(target3d, 2))
            for j in 1:nlev
                target3d[1, j, pidx] = ndims(data) == 1 ? data[j] : data[j, 1]
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

    # Helper: inject 2D patch variable. Fortran layout is (sec × patch) where
    # sec is a radiation band (numrad) or a canopy layer (nlevcan); CLM.jl stores
    # the transpose [patch, sec]. Reuses the PFT-type patch remap of set_patch_1d!.
    function set_patch_2d!(name, target_array)
        haskey(ds, name) || return
        try
            raw = ds[name][:, :]                       # (nsec, npatch_f90)
            f90_pfts = haskey(ds, "pfts1d_itypveg") ? Int.(ds["pfts1d_itypveg"][:]) : Int[]
            jl_pfts = inst.patch.itype
            nsec = min(size(raw, 1), size(target_array, 2))
            for pf in 1:min(np_f90, size(raw, 2))
                pj_found = 0
                for pj in 1:np_jl
                    if !isempty(f90_pfts) && f90_pfts[pf] == jl_pfts[pj]
                        pj_found = pj; break
                    elseif pf <= np_jl
                        pj_found = pf; break
                    end
                end
                pj_found == 0 && continue
                for s in 1:nsec
                    v = raw[s, pf]
                    target_array[pj_found, s] = ismissing(v) ? NaN : Float64(v)
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
    # Recompute h2osoi_vol_col from the injected liq/ice. h2osoi_vol is normally a
    # mid-step diagnostic refreshed AFTER SurfaceAlbedo runs, so on a restart-read
    # it otherwise keeps the cold-start default (0.75*watsat, wet) — which biases
    # the soil-albedo moisture blend (inc = 0.11 - 0.40*h2osoi_vol[c,1]) by ~6.5%.
    # Formula matches lake_hydrology.jl / CLM (jj = j + nlevsno snow padding).
    if !isempty(ws.h2osoi_vol_col) && !isempty(inst.column.dz)
        dz = inst.column.dz
        for c in 1:nc_jl, j in 1:nlevgrnd
            jj = j + nlevsno
            (jj <= size(ws.h2osoi_liq_col, 2) && jj <= size(dz, 2)) || continue
            ws.h2osoi_vol_col[c, j] = ws.h2osoi_liq_col[c, jj] / (dz[c, jj] * DENH2O) +
                                      ws.h2osoi_ice_col[c, jj] / (dz[c, jj] * DENICE)
        end
    end
    # Init eff_porosity from the injected ice. eff_porosity is normally computed mid-
    # step by HydrologyNoDrainage, but compute_h2osoi_liqvol! AND the use_hydrstress
    # PHS plant-water solve (canopy_fluxes) READ it at step 1 BEFORE that. On a bare
    # restart it would otherwise be NaN → liqvol = min(lv, NaN) = NaN → smp_l →
    # k_soil_root → vegwp Newton all NaN (a use_hydrstress restart starting on frozen
    # soil blew up here). Formula matches soil_hydrology.jl:104-108.
    let ss = inst.soilstate
        if !isempty(ss.eff_porosity_col) && !isempty(ss.watsat_col) && !isempty(inst.column.dz)
            dz = inst.column.dz
            has_exice = isdefined(ws, :excess_ice_col) && !isempty(ws.excess_ice_col)
            for c in 1:nc_jl, j in 1:nlevgrnd
                jj = j + nlevsno
                (jj <= size(ws.h2osoi_ice_col, 2) && jj <= size(dz, 2) &&
                 j <= size(ss.eff_porosity_col, 2) && j <= size(ss.watsat_col, 2)) || continue
                ei = has_exice ? ws.excess_ice_col[c, j] : 0.0
                dz_ext = dz[c, jj] + ei / DENICE
                vol_ice = min(ss.watsat_col[c, j], (ws.h2osoi_ice_col[c, jj] + ei) / (dz_ext * DENICE))
                ss.eff_porosity_col[c, j] = max(0.01, ss.watsat_col[c, j] - vol_ice)
            end
        end
    end
    set_col_1d!("WA", v -> (ws.wa_col[1] = v))
    set_col_1d!("H2OSFC", v -> (ws.h2osfc_col[1] = v))

    # Active-layer annual-max trackers (carried prognostic state). altmax/altmax_indx
    # are re-tracked each step by alt_calc!, but altmax_lastyear is only refreshed at
    # the year boundary, so it must be injected to match the Fortran restart.
    al = inst.active_layer
    set_col_1d!("altmax",               v -> (al.altmax_col[1]               = v))
    set_col_1d!("altmax_lastyear",      v -> (al.altmax_lastyear_col[1]      = v))
    set_col_1d!("altmax_indx",          v -> (al.altmax_indx_col[1]          = Int(round(v))))
    set_col_1d!("altmax_lastyear_indx", v -> (al.altmax_lastyear_indx_col[1] = Int(round(v))))

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
    if hasproperty(inst.temperature, :t_stem_patch)
        set_patch_1d!("T_STEM", inst.temperature.t_stem_patch)
    end
    # 10-day mean air temperature (photosynthesis jmax25 input; a restart-persisted
    # accumulator mean in CTSM). Without it t_a10=NaN → vcmax/jmax NaN → t_veg NaN.
    if hasproperty(inst.temperature, :t_a10_patch)
        set_patch_1d!("tair10", inst.temperature.t_a10_patch)
    end

    # Canopy state. The restart stores the PRE-SatellitePhenology veg structure
    # (the value persisted at the END of the previous step). At a cold-start FIRST
    # step that is all-zero, whereas the flux calcs need the POST-SP structure that
    # SatellitePhenology sets in-step from the surfdata monthly LAI/SAI/height —
    # which cold_start! already computed identically. Injecting the restart zeros
    # would clobber that into a leafless "canopy" (elai=0 but frac_veg_nosno_alb=1,
    # htop=0 → log(0) in the z0 blend → t_veg=NaN). `inject_veg_struct=false` keeps
    # the cold_start surfdata structure so the single-step flux matches Fortran's
    # post-SP state (needed off-Bow, where dumps carry pre-SP zeros).
    if inject_veg_struct
        set_patch_1d!("elai", inst.canopystate.elai_patch)
        set_patch_1d!("esai", inst.canopystate.esai_patch)
        set_patch_1d!("tlai", inst.canopystate.tlai_patch)
        set_patch_1d!("tsai", inst.canopystate.tsai_patch)
        set_patch_1d!("htop", inst.canopystate.htop_patch)
        set_patch_1d!("hbot", inst.canopystate.hbot_patch)
    end
    set_patch_1d!("fsun", inst.canopystate.fsun_patch)

    # PHS plant water potential per veg segment (sun/sha/xyl/root). Restart-persisted
    # in CTSM; without it the hydraulic-stress solver starts from the cold default
    # (vegwp≈0 → plc≈1 → unstressed) and BTRANMN spikes on day 1 to ~0.98 vs Fortran's
    # ~0.13 — a pure IC artifact (day 2 onward already matches once vegwp draws down).
    # Restart `vegwp` is (pft, vegwcs) in CDL → NCDatasets (seg, pft) = the (nsec,
    # npatch) layout set_patch_2d! expects. No-op when absent / non-PHS (empty target).
    if hasproperty(inst.canopystate, :vegwp_patch)
        set_patch_2d!("vegwp", inst.canopystate.vegwp_patch)
    end

    # BTRANMN accumulator state (restart-persisted in CTSM). BTRANMN history points at
    # the COMPLETED daily-min (btran_min); on day 1 there is no completed previous day,
    # so Fortran reports its restart-persisted BTRAN_MIN (≈0.13 here). Julia left these
    # at NaN → the history fell back to the instantaneous btran whose day-1 average is
    # ~0.98 (a pure IC spike; day 2 onward already matches). Inject BTRAN_MIN + the
    # hourly-average accumulator (BTRANAV_VALUE/NSTEPS). btran_min_inst starts the day
    # fresh, so seed it to SPVAL (the restart stores NaN = "no running min yet"); a NaN
    # there would poison the day-1 min(avg, NaN) reduction.
    ef = inst.energyflux
    if hasproperty(ef, :btran_min_patch) && !isempty(ef.btran_min_patch)
        set_patch_1d!("BTRAN_MIN", ef.btran_min_patch)
        if hasproperty(ef, :btranav_accum_patch) && !isempty(ef.btranav_accum_patch)
            set_patch_1d!("BTRANAV_VALUE", ef.btranav_accum_patch)
            set_patch_1d!("BTRANAV_NSTEPS", ef.btranav_naccum_patch)
        end
        if hasproperty(ef, :btran_min_inst_patch) && !isempty(ef.btran_min_inst_patch)
            for p in eachindex(ef.btran_min_inst_patch)
                ef.btran_min_inst_patch[p] = SPVAL
            end
        end
    end

    # Exposed-veg flag from the PREVIOUS step's SurfaceAlbedo (persisted to the
    # restart). drv_init copies it into frac_veg_nosno → the exposed/noexposed-veg
    # filters, deciding canopy-solve vs bareground for THIS step. CLM ramps it on
    # over the first ~2–3 cold-start steps (elai/SAI establishing), so injecting it
    # makes Julia's canopy/bareground split — and the resulting surface fluxes —
    # match Fortran step-for-step instead of activating the canopy immediately from
    # cold_start's elai. Integer field → round (SPVAL/missing → leave cold_start).
    if haskey(ds, "FRAC_VEG_NOSNO_ALB")
        fva = ds["FRAC_VEG_NOSNO_ALB"][:]
        for p in 1:min(np_jl, length(fva))
            ismissing(fva[p]) && continue
            inst.canopystate.frac_veg_nosno_alb_patch[p] = Int(round(Float64(fva[p])))
        end
    end

    # SurfaceAlbedo radiative-transfer outputs from the PREVIOUS step. CLM calls
    # SurfaceAlbedo at the end of clm_driver and saves these to the restart; the
    # NEXT step's SurfaceRadiation consumes them to split absorbed solar between
    # canopy (sabv) and ground (sabg). Without them fabd/fabi=0 → the canopy is
    # radiatively invisible, all solar hits the ground, and the ground over-heats.
    sa = inst.surfalb
    set_patch_2d!("fabd",     sa.fabd_patch)
    set_patch_2d!("fabi",     sa.fabi_patch)
    set_patch_2d!("fabd_sun", sa.fabd_sun_patch)
    set_patch_2d!("fabd_sha", sa.fabd_sha_patch)
    set_patch_2d!("fabi_sun", sa.fabi_sun_patch)
    set_patch_2d!("fabi_sha", sa.fabi_sha_patch)
    set_patch_2d!("ftdd",     sa.ftdd_patch)
    set_patch_2d!("ftid",     sa.ftid_patch)
    set_patch_2d!("ftii",     sa.ftii_patch)
    set_patch_2d!("albd",     sa.albd_patch)
    set_patch_2d!("albi",     sa.albi_patch)
    set_patch_2d!("fabd_sun_z", sa.fabd_sun_z_patch)
    set_patch_2d!("fabd_sha_z", sa.fabd_sha_z_patch)
    set_patch_2d!("fabi_sun_z", sa.fabi_sun_z_patch)
    set_patch_2d!("fabi_sha_z", sa.fabi_sha_z_patch)
    set_patch_2d!("tlai_z",   sa.tlai_z_patch)
    set_patch_2d!("tsai_z",   sa.tsai_z_patch)
    set_patch_2d!("fsun_z",   sa.fsun_z_patch)
    set_patch_1d!("vcmaxcintsun", sa.vcmaxcintsun_patch)
    set_patch_1d!("vcmaxcintsha", sa.vcmaxcintsha_patch)
    if haskey(ds, "nrad") && hasproperty(sa, :nrad_patch)
        set_patch_1d!("nrad", sa.nrad_patch)
    end
    # LUNA-acclimated vcmax25/jmax25 (dump dims (levcan, pft) → Julia [patch, levcan]).
    # Only injected when the photosyns LUNA fields are allocated (use_luna harness).
    if haskey(ds, "vcmx25_z") && hasproperty(inst.photosyns, :vcmx25_z_patch) &&
       !isempty(inst.photosyns.vcmx25_z_patch)
        set_patch_2d!("vcmx25_z", inst.photosyns.vcmx25_z_patch)
        set_patch_2d!("jmx25_z",  inst.photosyns.jmx25_z_patch)
    end

    # Canopy water
    set_patch_1d!("LIQCAN", ws.liqcan_patch)
    set_patch_1d!("SNOCAN", ws.snocan_patch)
    set_patch_1d!("FWET", wd.fwet_patch)
    set_patch_1d!("FCANSNO", wd.fcansno_patch)

    # Friction velocity
    set_patch_1d!("OBU", inst.frictionvel.obu_patch)

    # H2OSNO_NO_LAYERS holds snow water ONLY for columns with no explicit snow
    # layers (snl == 0). For a layered column (snl < 0) the snow water lives in
    # the layers, so this must be 0 — otherwise downstream h2osno_total
    # (= h2osno_no_layers + Σ layers) double-counts the pack and, via the
    # snow-cover-fraction overwrite, doubles int_snow.
    if inst.column.snl[1] < 0
        ws.h2osno_no_layers_col[1] = 0.0
    elseif haskey(ds, "H2OSNO_NO_LAYERS") && !ismissing(ds["H2OSNO_NO_LAYERS"][1])
        # No explicit snow layers (snl == 0): the trace snow water IS the dump's
        # H2OSNO_NO_LAYERS field — use it directly. Do NOT recompute it from the
        # snow-layer slots h2osoi_{liq,ice}[1, 1:nlevsno]: those slots are not
        # injected for a no-layer column, so they hold the cold-start NaN, which
        # would NaN h2osno_total → route the snow-depth kernel into its reset
        # branch (the pack never accumulates) and NaN int_snow.
        ws.h2osno_no_layers_col[1] = Float64(ds["H2OSNO_NO_LAYERS"][1])
    else
        # Dump lacks the field: fall back to the snow-layer sum, guarding NaN slots.
        h2osno_total = 0.0
        for j in 1:nlevsno
            liq = ws.h2osoi_liq_col[1, j]; ice = ws.h2osoi_ice_col[1, j]
            h2osno_total += (isnan(liq) ? 0.0 : liq) + (isnan(ice) ? 0.0 : ice)
        end
        ws.h2osno_no_layers_col[1] = h2osno_total
    end

    # =====================================================================
    # CN / BGC POOLS (use_cn). Only populated when the inst was built with
    # use_cn=true; on an SP inst the target arrays are empty (size 0) and the
    # guards below skip the whole block (no spurious n_fail).
    # =====================================================================
    if hasproperty(inst, :bgc_vegetation) &&
       hasproperty(inst.bgc_vegetation, :cnveg_carbonstate_inst)
        ccs = inst.bgc_vegetation.cnveg_carbonstate_inst
        cns = inst.bgc_vegetation.cnveg_nitrogenstate_inst
        if length(getfield(ccs, :leafc_patch)) > 0
            # --- vegetation carbon pools (gC/m2, per pft). Dump var name == Julia
            #     field stem; field = <stem>_patch.
            for stem in ("leafc","leafc_storage","leafc_xfer",
                         "frootc","frootc_storage","frootc_xfer",
                         "livestemc","livestemc_storage","livestemc_xfer",
                         "deadstemc","deadstemc_storage","deadstemc_xfer",
                         "livecrootc","livecrootc_storage","livecrootc_xfer",
                         "deadcrootc","deadcrootc_storage","deadcrootc_xfer",
                         "gresp_storage","gresp_xfer","cpool","xsmrpool")
                fld = Symbol(stem * "_patch")
                hasproperty(ccs, fld) && set_patch_1d!(stem, getfield(ccs, fld))
            end
            hasproperty(ccs, :ctrunc_patch) && set_patch_1d!("pft_ctrunc", ccs.ctrunc_patch)
            # --- vegetation nitrogen pools (gN/m2, per pft)
            for stem in ("leafn","leafn_storage","leafn_xfer",
                         "frootn","frootn_storage","frootn_xfer",
                         "livestemn","livestemn_storage","livestemn_xfer",
                         "deadstemn","deadstemn_storage","deadstemn_xfer",
                         "livecrootn","livecrootn_storage","livecrootn_xfer",
                         "deadcrootn","deadcrootn_storage","deadcrootn_xfer",
                         "retransn","npool")
                fld = Symbol(stem * "_patch")
                hasproperty(cns, fld) && set_patch_1d!(stem, getfield(cns, fld))
            end
            hasproperty(cns, :ntrunc_patch) && set_patch_1d!("pft_ntrunc", cns.ntrunc_patch)
        end
        # CNVeg phenology accumulators (restart-persisted; needed for cn_phenology!
        # parity). Dump var name == Julia field stem; field = <stem>_patch.
        if hasproperty(inst.bgc_vegetation, :cnveg_state_inst)
            cvs = inst.bgc_vegetation.cnveg_state_inst
            if hasproperty(cvs, :dormant_flag_patch) && length(cvs.dormant_flag_patch) > 0
                for stem in ("dormant_flag","days_active","onset_flag","onset_counter",
                             "onset_gddflag","onset_fdd","onset_gdd","onset_swi",
                             "offset_flag","offset_counter","offset_fdd","offset_swi",
                             "lgsf","bglfr","bgtr","annavg_t2m","tempavg_t2m",
                             "leafcn_offset",   # FUN leaf C:N target (prognostic, restart-persisted)
                             # GPP/retrans accumulators (were NaN → drive the seasonal
                             # allocation + N-demand scaling; restart-persisted)
                             "tempsum_potential_gpp","annsum_potential_gpp",
                             "tempmax_retransn","annmax_retransn")
                    fld = Symbol(stem * "_patch")
                    hasproperty(cvs, fld) && set_patch_1d!(stem, getfield(cvs, fld))
                end
            end
            # annsum_npp lives on the carbon-flux struct (24h/annual NPP accumulator)
            if hasproperty(inst.bgc_vegetation, :cnveg_carbonflux_inst)
                cvf = inst.bgc_vegetation.cnveg_carbonflux_inst
                hasproperty(cvf, :annsum_npp_patch) && set_patch_1d!("annsum_npp", cvf.annsum_npp_patch)
                # leaf litterfall used by FUN for retranslocation accounting. A
                # restart var: FUN (which runs before phenology phase-2) reads the
                # PREVIOUS step's value, so it must be injected, not left at 0.
                hasproperty(cvf, :leafc_to_litter_fun_patch) &&
                    set_patch_1d!("leafc_to_litter_fun", cvf.leafc_to_litter_fun_patch)
            end
        end
    end

    # Soil decomposition pools (vertically resolved). Pool index order (CNMode):
    # 1 litr1, 2 litr2, 3 litr3, 4 soil1, 5 soil2, 6 soil3, 7 cwd.
    if hasproperty(inst, :soilbiogeochem_carbonstate)
        scs = inst.soilbiogeochem_carbonstate
        sns = inst.soilbiogeochem_nitrogenstate
        cpools = ("litr1c_vr"=>1,"litr2c_vr"=>2,"litr3c_vr"=>3,
                  "soil1c_vr"=>4,"soil2c_vr"=>5,"soil3c_vr"=>6,"cwdc_vr"=>7)
        for (vn, pidx) in cpools
            set_col_pool_vr!(vn, scs.decomp_cpools_vr_col, pidx)
        end
        npools = ("litr1n_vr"=>1,"litr2n_vr"=>2,"litr3n_vr"=>3,
                  "soil1n_vr"=>4,"soil2n_vr"=>5,"soil3n_vr"=>6,"cwdn_vr"=>7)
        for (vn, pidx) in npools
            set_col_pool_vr!(vn, sns.decomp_npools_vr_col, pidx)
        end
        # mineral N + truncation sinks (column, levgrnd) → matrices [nc, nlevdecomp]
        size(sns.sminn_vr_col, 1) > 0    && set_col_2d!("sminn_vr",    sns.sminn_vr_col)
        size(sns.smin_no3_vr_col, 1) > 0 && set_col_2d!("smin_no3_vr", sns.smin_no3_vr_col)
        size(sns.smin_nh4_vr_col, 1) > 0 && set_col_2d!("smin_nh4_vr", sns.smin_nh4_vr_col)
        size(sns.ntrunc_vr_col, 1) > 0   && set_col_2d!("ntrunc_vr",   sns.ntrunc_vr_col)
        size(scs.ctrunc_vr_col, 1) > 0   && set_col_2d!("ctrunc_vr",   scs.ctrunc_vr_col)
    end

    close(ds)

    @info "Fortran restart loaded: $n_ok vars set, $n_fail failed, $n_skip skipped"
    return n_ok
end
