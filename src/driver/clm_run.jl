# ==========================================================================
# Top-level simulation loop
# Orchestrates initialization, forcing, driver, and output
#
# Public functions:
#   clm_run!  — Run a complete CLM simulation from start to end
# ==========================================================================

"""
    clm_run!(; fsurdat, paramfile, fforcing, kwargs...)
    -> CLMInstances

Run a complete CLM simulation from initialization through time integration
to output. This is the top-level entry point for offline CLM simulations.

# Required Arguments
- `fsurdat::String`   — Path to surface data NetCDF file
- `paramfile::String` — Path to CLM parameter NetCDF file
- `fforcing::String`  — Path to atmospheric forcing NetCDF file

# Optional Arguments
- `fhistory::String`     — Output history file (default "clm_history.nc")
- `start_date::DateTime` — Simulation start date (default 2000-01-01)
- `end_date::DateTime`   — Simulation end date (default 2000-02-01)
- `dtime::Int`           — Timestep in seconds (default 1800)
- `use_cn::Bool`         — Use CN biogeochemistry (default false)
- `use_bedrock::Bool`    — Use bedrock-limited soil column (default true)
- `use_aquifer_layer::Bool` — Use aquifer lower boundary (default true)
- `hist_fields`          — Custom history fields (default: default_hist_fields())
- `verbose::Bool`        — Print progress messages (default true)
- `h2osfcflag::Int`      — Surface water flag for soil hydrology (default 0)

# Returns
- `inst::CLMInstances` — Final state of all CLM data instances
"""
function clm_run!(;
    fsurdat::String,
    paramfile::String,
    fforcing::String,
    fhistory::String = "clm_history.nc",
    start_date::DateTime = DateTime(2000, 1, 1),
    end_date::DateTime = DateTime(2000, 2, 1),
    dtime::Int = 1800,
    use_cn::Bool = false,
    use_bedrock::Bool = true,
    use_aquifer_layer::Bool = true,
    h2osfcflag::Int = 0,
    hist_fields::Union{Vector{HistFieldDef}, Nothing} = nothing,
    verbose::Bool = true,
    fsnowoptics::String = "",
    fsnowaging::String = "",
    frestart::String = "",
    baseflow_scalar::Real = 1.0e-2,
    int_snow_max::Real = 2000.0,
    overrides::Union{CalibrationOverrides, Nothing} = nothing)

    # ========================================================================
    # Phase 1: Initialize
    # ========================================================================
    verbose && println("CLM.jl: Initializing...")
    (inst, bounds, filt, tm) = clm_initialize!(;
        fsurdat=fsurdat, paramfile=paramfile,
        start_date=start_date, dtime=dtime, use_cn=use_cn,
        use_bedrock=use_bedrock,
        use_aquifer_layer=use_aquifer_layer,
        h2osfcflag=h2osfcflag,
        fsnowoptics=fsnowoptics, fsnowaging=fsnowaging,
        int_snow_max=int_snow_max)

    if overrides !== nothing
        inst.overrides = overrides
    end

    config = CLMDriverConfig(use_cn=use_cn, use_aquifer_layer=use_aquifer_layer)

    # Set baseflow scalar (namelist-adjustable parameter)
    bf = (overrides !== nothing && !isnan(overrides.baseflow_scalar)) ?
        overrides.baseflow_scalar : baseflow_scalar
    init_soil_hydrology_config(baseflow_scalar=bf)

    # Wire fff override into sat_excess_runoff_params
    if overrides !== nothing && !isnan(overrides.fff)
        sat_excess_runoff_params.fff = overrides.fff
    end

    # Wire ksat_scale override: multiply hksat in soil state
    if overrides !== nothing && !isnan(overrides.ksat_scale)
        inst.soilstate.hksat_col .*= overrides.ksat_scale
    end

    filt_ia = clump_filter_inactive_and_active

    ng = bounds.endg
    nc = bounds.endc
    np = bounds.endp

    # Apply soil property multipliers to soilstate (post-init, since CLM.jl
    # computes bsw/watsat/sucsat from pedotransfer, ignoring surfdata values)
    if overrides !== nothing
        nlevsoi_l = varpar.nlevsoi
        for (field, mult_val, lo, hi) in [
            (:bsw_col, overrides.bsw_mult, 1.0, 30.0),
            (:watsat_col, overrides.watsat_mult, 0.01, 0.95),
            (:sucsat_col, overrides.sucsat_mult, 1.0, 1e6),
        ]
            if !isnan(mult_val) && mult_val != 1.0
                arr = getfield(inst.soilstate, field)
                for c in 1:nc, j in 1:nlevsoi_l
                    arr[c, j] = clamp(arr[c, j] * mult_val, lo, hi)
                end
            end
        end
    end

    # ---- Wire params from params.nc into runtime structs ----
    # Must be after nc is defined. readParameters! handles PFT params and some
    # scalars; snow/hydrology params need explicit wiring here.
    try
        ds_p = NCDataset(paramfile, "r")
        scf = inst.scf_method
        if haskey(ds_p, "n_melt_coef")
            nmc = Float64(ds_p["n_melt_coef"][1])
            for c in 1:nc
                topo_std = max(10.0, inst.column.topo_std[c])
                scf.n_melt[c] = nmc / topo_std
            end
        end
        haskey(ds_p, "accum_factor") && (scf.accum_factor = Float64(ds_p["accum_factor"][1]))
        haskey(ds_p, "SNOW_DENSITY_MAX") && (snowhydrology_params.rho_max = Float64(ds_p["SNOW_DENSITY_MAX"][1]))
        haskey(ds_p, "SNOW_DENSITY_MIN") && (snowhydrology_params.rho_min = Float64(ds_p["SNOW_DENSITY_MIN"][1]))
        haskey(ds_p, "fresh_snw_rds_max") && (snowhydrology_params.snw_rds_min = Float64(ds_p["fresh_snw_rds_max"][1]))
        haskey(ds_p, "SNO_Z0MV") && (inst.frictionvel.zsno = Float64(ds_p["SNO_Z0MV"][1]))
        if haskey(ds_p, "snw_aging_bst")
            snicar_params.xdrdt = Float64(ds_p["snw_aging_bst"][1])
        elseif haskey(ds_p, "xdrdt")
            snicar_params.xdrdt = Float64(ds_p["xdrdt"][1])
        end
        if haskey(ds_p, "pc")
            pc_val = Float64(ds_p["pc"][1])
            pc_val > 0 && (for c in 1:nc; inst.soilhydrology.hkdepth_col[c] = 1.0 / pc_val; end)
        end
        close(ds_p)
    catch e
        @warn "Runtime param wiring failed: $e" maxlog=1
    end

    # ---- Read restart if provided ----
    if !isempty(frestart) && isfile(frestart)
        verbose && println("CLM.jl: Reading restart from: ", frestart)
        read_restart!(frestart, inst, bounds)
    end

    # ---- Initialize water table near bedrock for cold-start runs ----
    # Without restart, wa_col defaults to 4000mm which puts zwt ~13m deep
    # (well below bedrock). This disables baseflow for years. Set wa to put
    # zwt at bedrock depth so hydrology params are immediately active.
    if isempty(frestart) || !isfile(frestart)
        nlevsno_l = varpar.nlevsno
        nlevsoi_l = varpar.nlevsoi
        for c in 1:nc
            nbr = inst.column.nbedrock[c]
            zi_bed = inst.column.zi[c, nbr + nlevsno_l + 1]
            zi_bot = inst.column.zi[c, nlevsoi_l + nlevsno_l + 1]
            target_zwt = zi_bed  # match Fortran equilibrium: zwt at bedrock
            wa_target = (25.0 + zi_bot - target_zwt) * 0.2 * 1000.0
            inst.water.waterstatebulk_inst.ws.wa_col[c] = wa_target
            inst.soilhydrology.zwt_col[c] = target_zwt
        end
    end

    # ========================================================================
    # Phase 2: Open forcing file
    # ========================================================================
    verbose && println("CLM.jl: Opening forcing file: ", fforcing)
    fr = ForcingReader()
    forcing_reader_init!(fr, fforcing)

    # Read topo forcing to set atmospheric topography for lapse-rate correction.
    # Without this, forc_topo_grc=0 and downscale_forcings! doesn't correct T.
    topo_file = replace(fforcing, r"clmforc\.[^/]*\.nc$" => "topo_forcing.nc")
    if isfile(topo_file)
        try
            ds_topo = NCDataset(topo_file, "r")
            if haskey(ds_topo, "TOPO")
                forc_topo = Float64(ds_topo["TOPO"][1])
                for g in 1:ng
                    inst.atm2lnd.forc_topo_grc[g] = forc_topo
                end
                # Set topo_col to surfdata TOPO (the physical column elevation).
                # The forcing is at forc_topo (basin-average elevation from EASYMORE);
                # downscale_forcings! corrects: T_col = T_forc - lapse*(topo_col - forc_topo)
                try
                    ds_sf2 = NCDataset(fsurdat, "r")
                    surf_topo = haskey(ds_sf2, "TOPO") ? Float64(ds_sf2["TOPO"][1]) : forc_topo
                    for c2 in 1:nc
                        inst.topo.topo_col[c2] = surf_topo
                    end
                    close(ds_sf2)
                    verbose && println("CLM.jl: Lapse: forcing=$(round(forc_topo,digits=0))m → surface=$(round(surf_topo,digits=0))m (ΔT=$(round(0.006*(surf_topo-forc_topo),digits=1))K)")
                catch; end
                verbose && println("CLM.jl: Forcing elevation: $(round(forc_topo, digits=0))m, surface: $(round(inst.topo.topo_col[1], digits=0))m")
            end
            close(ds_topo)
        catch e
            @debug "Topo forcing read: $e"
        end
    end

    # ========================================================================
    # Phase 3: Open history writer
    # ========================================================================
    fields = hist_fields !== nothing ? hist_fields : default_hist_fields()
    hw = HistoryWriter()
    history_writer_init!(hw, fhistory, fields, ng, nc, np)
    verbose && println("CLM.jl: Writing output to: ", fhistory)

    # ========================================================================
    # Phase 4: Time integration loop
    # ========================================================================
    total_steps = Int(ceil(Dates.value(end_date - start_date) / (1000 * dtime)))
    step_count = 0

    verbose && println("CLM.jl: Starting time loop ($total_steps steps, dtime=$(dtime)s)")

    while tm.current_date < end_date
        advance_timestep!(tm)
        step_count += 1

        # --- Read and downscale forcings ---
        read_forcing_step!(fr, inst.atm2lnd, tm.current_date, ng, nc)
        downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)

        # --- Compute orbital parameters ---
        calday = get_curr_calday(tm)
        (declin, eccf) = compute_orbital(calday)
        declinp1 = declin  # simplified: same as current step
        obliqr = ORB_OBLIQR_DEFAULT

        # --- Time flags ---
        doalb = true  # compute albedo every step for simplicity
        nextsw_cday = calday + Float64(dtime) / SECSPDAY

        (yr, mon, d, tod) = get_curr_date(tm)
        first_step = (tm.nstep == 1)
        beg_day = is_beg_curr_day(tm)
        end_day = is_end_curr_day(tm)
        beg_year = is_beg_curr_year(tm)

        # --- Run driver ---
        clm_drv!(config, inst, filt, filt_ia, bounds,
                 doalb, nextsw_cday, declinp1, declin, obliqr,
                 false, false, "", false;
                 nstep=tm.nstep,
                 is_first_step=first_step,
                 is_beg_curr_day=beg_day,
                 is_end_curr_day=end_day,
                 is_beg_curr_year=beg_year,
                 dtime=Float64(dtime),
                 mon=mon,
                 day=d,
                 photosyns=inst.photosyns)

        # --- Land to atmosphere aggregation ---
        lnd2atm!(bounds, inst)

        # --- Write history output ---
        history_write_step!(hw, inst, tm.current_date; is_end_curr_day=end_day)

        # --- Progress ---
        if verbose && (step_count % 48 == 0 || step_count == total_steps)
            println("CLM.jl: Step $step_count/$total_steps — $(tm.current_date)")
        end
    end

    # ========================================================================
    # Phase 5: Write restart and cleanup
    # ========================================================================

    # Write restart file for continuation runs
    if !isempty(frestart)
        restart_out = replace(frestart, r"\.nc$" => "") * "_out.nc"
        verbose && println("CLM.jl: Writing restart to: ", restart_out)
        write_restart!(restart_out, inst, bounds; time=tm.current_date)
    end

    forcing_reader_close!(fr)
    history_writer_close!(hw)

    verbose && println("CLM.jl: Simulation complete. $step_count timesteps written to $fhistory")

    return inst
end
