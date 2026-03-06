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
    frestart::String = "")

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
        fsnowoptics=fsnowoptics, fsnowaging=fsnowaging)

    config = CLMDriverConfig(use_cn=use_cn, use_aquifer_layer=use_aquifer_layer)
    filt_ia = clump_filter_inactive_and_active

    ng = bounds.endg
    nc = bounds.endc
    np = bounds.endp

    # ---- Read restart if provided ----
    if !isempty(frestart) && isfile(frestart)
        verbose && println("CLM.jl: Reading restart from: ", frestart)
        read_restart!(frestart, inst, bounds)
    end

    # ========================================================================
    # Phase 2: Open forcing file
    # ========================================================================
    verbose && println("CLM.jl: Opening forcing file: ", fforcing)
    fr = ForcingReader()
    forcing_reader_init!(fr, fforcing)

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
        history_write_step!(hw, inst, tm.current_date)

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
