# ==========================================================================
# run_clm.jl — standalone run harness.
#
# A thin orchestration of the now-existing pieces so a run produces history
# output + a restart WITHOUT the Fortran parity harness:
#
#   clm_initialize!   (caller-built `inst`/`bounds`/`filt`/`tm`)
#   forcing_reader_*  (read one forcing timestep)
#   downscale_forcings!
#   clm_drv!          (one biogeophysical timestep)
#   history_io.jl     (HistoryTape: hist_accumulate! / hist_write!)
#   restart_io.jl     (write_restart / read_restart!)
#
# This file does NOT define new physics and does NOT modify clm_drv! or
# CLMInstances — there is no AD risk. It mirrors the driver-call pattern in
# scripts/run_clm_streamflow.jl + scripts/longhorizon_conservation.jl exactly.
#
# Public function:
#   run_clm!  — start-from-cold-or-restart, loop the driver, write h0 + restart
# ==========================================================================

"""
    run_clm!(config, inst, bounds, filt, tm, fr;
             nsteps, dtime=3600.0, start_date=nothing,
             hist_tape=nothing, hist_interval=0, hist_path="clm.h0.nc",
             restart_path=nothing, finidat=nothing, use_cn=false,
             orb_obliqr=ORB_OBLIQR_DEFAULT, verbose=false)

Run the CLM.jl biogeophysical driver for `nsteps` timesteps and (optionally)
produce a CLM-style h0 history file and a restart file — a standalone run that
needs no parity/forcing-dump machinery beyond a `ForcingReader` already opened
on a forcing NetCDF.

The caller supplies the live state built by [`clm_initialize!`]:
`config::CLMDriverConfig`, `inst::CLMInstances`, `bounds::BoundsType`, the
filter `filt`, the time manager `tm::TimeManager`, and a `ForcingReader` `fr`
already initialized (`forcing_reader_init!`) on the run's forcing file.

Keyword arguments
- `nsteps`        — number of timesteps to advance.
- `dtime`         — timestep length in seconds (default 3600).
- `start_date`    — the first timestep's date; defaults to `tm`'s current date.
- `finidat`       — if a path, `read_restart!` it into `inst` BEFORE the loop
                    (start-from-restart); the state is loaded and the run
                    continues from it.
- `hist_tape`     — a `HistoryTape`; if given, `hist_accumulate!` every step and
                    `hist_write!` to `hist_path` every `hist_interval` steps.
- `hist_interval` — write the h0 record every N steps (0 ⇒ write once at the end
                    if a tape is present and any samples remain).
- `hist_path`     — h0 NetCDF filename (records appended across writes).
- `restart_path`  — if a path, `write_restart` the final state there at the end.
- `use_cn`        — pass CN flag through to restart read/write (CN pools).
- `orb_obliqr`    — Earth obliquity for the solar-geometry / daylength update.
- `verbose`       — print a per-step / per-write progress line.

Returns a NamedTuple `(; nsteps, last_date, hist_path, restart_path)` recording
what was run and written (paths are `nothing` when not produced).
"""
function run_clm!(config::CLMDriverConfig, inst::CLMInstances, bounds::BoundsType,
                  filt, tm::TimeManager, fr::ForcingReader;
                  nsteps::Int,
                  dtime::Real = 3600.0,
                  start_date::Union{DateTime,Nothing} = nothing,
                  hist_tape = nothing,
                  hist_interval::Int = 0,
                  hist_path::AbstractString = "clm.h0.nc",
                  restart_path::Union{AbstractString,Nothing} = nothing,
                  finidat::Union{AbstractString,Nothing} = nothing,
                  use_cn::Bool = false,
                  orb_obliqr::Real = ORB_OBLIQR_DEFAULT,
                  verbose::Bool = false)

    nsteps >= 1 || error("run_clm!: nsteps must be >= 1, got $nsteps")
    ng, nc, np = bounds.endg, bounds.endc, bounds.endp

    # --- start-from-restart: load prognostic state before the loop ----------
    if finidat !== nothing
        read_restart!(inst, String(finidat); bounds = bounds, use_cn = use_cn)
        verbose && @info "run_clm!: loaded restart" finidat
    end

    # First timestep date (the loop uses `tm`'s notion of "now" otherwise).
    sdate = start_date === nothing ? DateTime(get_curr_date(tm)[1:3]...) : start_date

    filt_ia = clump_filter_inactive_and_active
    secs_per_day = SECSPDAY
    dt = Float64(dtime)
    hist_created = false   # whether the h0 file has been created yet (append after)
    last_written_step = 0
    last_date = sdate

    for step in 1:nsteps
        cur = sdate + Second(round(Int, (step - 1) * dt))
        last_date = cur

        # Solar geometry + daylength (mirror run_clm_streamflow.jl).
        calday = get_curr_calday(tm)
        (declin, _)   = compute_orbital(calday)
        (declinm1, _) = compute_orbital(calday - dt / secs_per_day)
        nextsw = calday + dt / secs_per_day
        init_daylength!(inst.gridcell, declin, declinm1, orb_obliqr, 1:ng)

        # Advance time, read + downscale forcing for this step.
        advance_timestep!(tm)
        read_forcing_step!(fr, inst.atm2lnd, cur, ng, nc;
                           gridcell_latdeg = inst.gridcell.latdeg,
                           gridcell_londeg = inst.gridcell.londeg)
        downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)

        (yr, mon, dy, tod) = get_curr_date(tm)
        clm_drv!(config, inst, filt, filt_ia, bounds, true, nextsw, declin, declin,
                 orb_obliqr, false, false, "", false;
                 nstep = tm.nstep, is_first_step = (step == 1),
                 is_beg_curr_day = is_beg_curr_day(tm),
                 is_end_curr_day = is_end_curr_day(tm),
                 is_beg_curr_year = is_beg_curr_year(tm),
                 dtime = dt, mon = mon, day = dy,
                 photosyns = inst.photosyns)

        # --- history accumulate + periodic write -----------------------------
        if hist_tape !== nothing
            hist_accumulate!(hist_tape, inst)
            if hist_interval > 0 && step % hist_interval == 0
                hist_write!(hist_tape, String(hist_path);
                            time = Dates.value(cur - DateTime(2000, 1, 1)) / 86400_000,
                            append = hist_created)
                hist_created = true
                last_written_step = step
                verbose && @info "run_clm!: wrote h0 record" step hist_path
            end
        end

        verbose && step % max(1, nsteps ÷ 10) == 0 &&
            @info "run_clm!: step" step nsteps date=cur
    end

    # Flush any remaining accumulated samples (final partial interval / no interval).
    if hist_tape !== nothing && hist_tape.nsamples > 0 && last_written_step < nsteps
        hist_write!(hist_tape, String(hist_path);
                    time = Dates.value(last_date - DateTime(2000, 1, 1)) / 86400_000,
                    append = hist_created)
        hist_created = true
        verbose && @info "run_clm!: wrote final h0 record" hist_path
    end

    # --- write restart at the end -------------------------------------------
    if restart_path !== nothing
        write_restart(inst, String(restart_path);
                      bounds = bounds, use_cn = use_cn, time = last_date)
        verbose && @info "run_clm!: wrote restart" restart_path
    end

    return (; nsteps = nsteps,
            last_date = last_date,
            hist_path = (hist_tape !== nothing && hist_created) ? String(hist_path) : nothing,
            restart_path = restart_path === nothing ? nothing : String(restart_path))
end
