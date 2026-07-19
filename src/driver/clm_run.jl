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
- `use_bedrock::Union{Bool,Nothing}` — Use bedrock-limited soil column. `nothing` (default)
  defers to clm_initialize!'s CTSM-conditional default; non-FATES resolves to `true`.
- `use_aquifer_layer::Bool` — Use aquifer lower boundary. Default `false` (CTSM clm5_0
  derivation); `true` selects the ZD09 + BC_AQUIFER solver.
- `hist_fields`          — Custom history fields (default: default_hist_fields())
- `verbose::Bool`        — Print progress messages (default true)
- `h2osfcflag::Int`      — Surface water flag for soil hydrology (default 1 = CTSM default)

# Returns
- `inst::CLMInstances` — Final state of all CLM data instances
"""
# Push every (host) array field of `src` into the matching (device) field of `dst`,
# converting to the destination element type (Float64 host forcing -> Float32 device).
# Used by the optional GPU path of clm_run! to hand the freshly-downscaled forcing to
# the device mirror each step.
function _sync_arrays_convert!(dst, src)
    T = typeof(src)
    @inbounds for i in 1:fieldcount(T)
        s = getfield(src, i); d = getfield(dst, i)
        if s isa AbstractArray && d isa AbstractArray && !isempty(s) && length(s) == length(d)
            copyto!(d, eltype(d).(s))
        end
    end
    return nothing
end

# Refresh a PERSISTENT host view `dst` (built once via `Adapt.adapt(Array, inst_d)`) from
# the evolving DEVICE tree `src`, IN-PLACE: recurse the struct tree and `copyto!` every
# matching leaf array (device -> host, same eltype -> no allocation). This replaces a
# per-step full-tree `Adapt.adapt(Array, ·)` SNAPSHOT whose per-step allocation OOM'd a
# full-year GPU run. `dst`/`src` share structure by construction (dst = adapt(Array, src)),
# so we walk fields positionally; arrays are copied, mutable sub-structs are recursed.
function _refresh_host_tree!(dst, src)
    @inbounds for i in 1:fieldcount(typeof(src))
        s = getfield(src, i); d = getfield(dst, i)
        if s isa AbstractArray && d isa AbstractArray
            (!isempty(s) && length(s) == length(d)) && copyto!(d, s)
        elseif ismutable(s) && !(s isa AbstractArray) && fieldcount(typeof(s)) > 0 &&
               fieldcount(typeof(s)) == fieldcount(typeof(d))
            _refresh_host_tree!(d, s)
        end
    end
    return nothing
end

function clm_run!(;
    fsurdat::String,
    paramfile::String,
    fforcing::String,
    fhistory::String = "clm_history.nc",
    start_date::DateTime = DateTime(2000, 1, 1),
    end_date::DateTime = DateTime(2000, 2, 1),
    dtime::Int = 1800,
    use_cn::Bool = false,
    # `nothing` defers to clm_initialize!'s CTSM-conditional default (.false. under FATES,
    # .true. otherwise). This driver is non-FATES, so it resolves to .true. as before.
    use_bedrock::Union{Bool,Nothing} = nothing,
    # CTSM DERIVES this from lower_boundary_condition rather than reading a
    # namelist flag: for phys=clm5_0 the chain is soilwater_movement_method=1 ->
    # lower_boundary_condition=2 (bc_zero_flux) -> use_aquifer_layer() = .false.
    # (namelist_defaults_ctsm.xml:419-428; SoilWaterMovementMod.F90:221-236).
    # The port defaulted `true`, which is CTSM's *code fallback* (bc_aquifer) --
    # the value the namelist overrides -- and `true` together with use_bedrock
    # is a combination CTSM itself endrun's on (SoilWaterMovementMod.F90:181).
    # It selects the whole soil-water solver: ZD09+BC_AQUIFER vs the CLM5
    # moisture-form solver + BC_ZERO_FLUX (see clm_driver.jl's swm_cfg switch).
    # 71 of the 75 explicit call sites in scripts/+test/ already pass `false`.
    use_aquifer_layer::Bool = false,
    use_luna::Bool = false,
    use_hydrstress::Bool = false,
    # CTSM namelist default (SoilHydrologyType.F90 clm_soilhydrology_inparm) is 1.
    # See the note on clm_initialize!'s h2osfcflag: this entry point defaulted to 0
    # and silently disabled the surface-water store.
    h2osfcflag::Int = 1,
    hist_fields::Union{Vector{HistFieldDef}, Nothing} = nothing,
    verbose::Bool = true,
    fsnowoptics::String = "",
    fsnowaging::String = "",
    frestart::String = "",
    ffortran_restart::String = "",
    baseflow_scalar::Real = 1.0e-2,
    int_snow_max::Real = 2000.0,
    interp_forcing::Bool = false,
    forcing_phase_shift_s::Int = 0,
    overrides::Union{CalibrationOverrides, Nothing} = nothing,
    step_probe::Union{Function, Nothing} = nothing,
    device_adapt::Union{Function, Nothing} = nothing)

    # ========================================================================
    # Phase 1: Initialize
    # ========================================================================
    verbose && println("CLM.jl: Initializing...")
    (inst, bounds, filt, tm) = clm_initialize!(;
        fsurdat=fsurdat, paramfile=paramfile,
        start_date=start_date, dtime=dtime, use_cn=use_cn,
        use_bedrock=use_bedrock,
        use_aquifer_layer=use_aquifer_layer,
        use_luna=use_luna,
        use_hydrstress=use_hydrstress,
        h2osfcflag=h2osfcflag,
        fsnowoptics=fsnowoptics, fsnowaging=fsnowaging,
        int_snow_max=int_snow_max)

    active_fields = Symbol[]
    if overrides !== nothing
        validate_overrides!(overrides)
        active_fields = active_override_fields(overrides)
        inst.overrides = overrides
    end

    # Light inhibition of leaf respiration is on by default in CLM5 (namelist
    # light_inhibit=.true. for the clm5 physics; the reference lnd_in sets it too).
    # The Julia PhotosynsData field defaults to false, so enable it here to match.
    inst.photosyns.light_inhibit = true

    config = CLMDriverConfig(use_cn=use_cn, use_aquifer_layer=use_aquifer_layer,
                             use_luna=use_luna, use_hydrstress=use_hydrstress)

    # Configure atm2lnd downscaling to match Fortran lnd_in defaults
    atm2lnd_read_namelist!(inst.atm2lnd;
        repartition_rain_snow=true,
        lapse_rate=0.006,
        lapse_rate_longwave=0.032,
        precip_repartition_nonglc_all_snow_t=0.0,
        precip_repartition_nonglc_all_rain_t=2.0,
        precip_repartition_glc_all_snow_t=-2.0,
        precip_repartition_glc_all_rain_t=0.0)

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
        # SNO_Z0MV is per-instance (frictionvel), so it must be wired here rather than
        # in readParameters!.
        haskey(ds_p, "SNO_Z0MV") && (inst.frictionvel.zsno = Float64(ds_p["SNO_Z0MV"][1]))
        # NOTE: the snow/SNICAR scalars that used to be re-read here
        # (SNOW_DENSITY_MAX/MIN, fresh_snw_rds_max, snw_aging_bst/xdrdt,
        # interception_fraction) are now wired once, in readParameters!
        # (infrastructure/read_params.jl), which clm_initialize! already calls.
        # Keeping two readers of the same file was how params got silently dropped.
        # NOTE: `pc` is the surface-water/frost connectivity threshold (Fortran reads it
        # into params_inst%pc for the frac_h2osfc fd calc + the SurfaceWater infiltration
        # cluster), NOT the baseflow decay depth. Fortran fixes hkdepth = 1/2.5 (0.4 m);
        # the old `hkdepth = 1/pc` override gave fff=1/2.5=0.4 instead of 2.5, exploding
        # the topographic baseflow exp(-fff*zwt) at wet sites (Aripuana drainage 14745
        # mm/d). Bow/Stillwater hid it (zero drainage). Leave hkdepth at its 1/2.5 default.
        # Canopy interception: CLM5 (clm5_0 physics, matching the Fortran reference
        # lnd_in use_clm5_fpi=.true.) uses fpiliq = interception_fraction*tanh(elai+esai),
        # NOT the CLM4 default 0.25*(1-exp(-0.5*lai)). The CLM4 form intercepts ~20% at
        # lai~3 vs ~64% for CLM5 with the calibrated interception_fraction=0.6455 → Julia
        # under-evaporated the canopy all summer (FCEV -17%). Enable CLM5 fpi + read the
        # calibrated interception_fraction (the latter now read in readParameters!).
        canopy_hydrology_read_nml!(use_clm5_fpi=true)
        close(ds_p)
    catch e
        @warn "Runtime param wiring failed: $e" maxlog=1
    end

    # ---- Read restart if provided ----
    if !isempty(frestart) && isfile(frestart)
        verbose && println("CLM.jl: Reading restart from: ", frestart)
        read_restart!(frestart, inst, bounds)
    end

    # ---- Inject a FORTRAN restart (spun-up IC) if provided ----
    # Lets the free-run start from the Fortran model's equilibrium prognostic
    # state (soil water/temps/snow/...), so an annual comparison against a
    # calibrated spun-up Fortran reference is spun-up-vs-spun-up rather than
    # cold-start-vs-spun-up (the latter leaves the cold-start soil too wet for
    # a year with deep drainage off → spurious summer BTRAN inversion).
    have_ic_restart = (!isempty(frestart) && isfile(frestart)) ||
                      (!isempty(ffortran_restart) && isfile(ffortran_restart))
    if !isempty(ffortran_restart) && isfile(ffortran_restart)
        verbose && println("CLM.jl: Injecting Fortran restart IC from: ", ffortran_restart)
        read_fortran_restart!(ffortran_restart, inst, bounds)
    end

    # ---- Initialize water table near bedrock for cold-start runs ----
    # Without restart, wa_col defaults to 4000mm which puts zwt ~13m deep
    # (well below bedrock). This disables baseflow for years. Set wa to put
    # zwt at bedrock depth so hydrology params are immediately active.
    if !have_ic_restart
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
    fr.interp_time = interp_forcing

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
                # Match Fortran TopoMod.UpdateTopo: for columns without glacier/hillslope
                # downscaling, topo_col = atm_topo (line 282 of TopoMod.F90).
                # This means no lapse correction (hsurf_c = hsurf_g → ΔT = 0).
                for c2 in 1:nc
                    inst.topo.topo_col[c2] = forc_topo
                end
                verbose && println("CLM.jl: topo_col=forc_topo=$(round(forc_topo,digits=0))m (no lapse correction)")
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

    # ---- Optional device (GPU) mirror ----
    # When `device_adapt` is supplied (e.g. moves the tree to a Metal MtlArray/Float32
    # backend), keep `inst` as the host forcing-in / history-out shuttle and run clm_drv!
    # on a device mirror. downscale_forcings! + history_write_step! are host-only, so per
    # step: read+downscale forcing on host -> push atm2lnd to device -> clm_drv! on device
    # -> device state -> host view -> lnd2atm! + history on the host view.
    _use_dev = device_adapt !== nothing
    inst_d = inst; filt_d = filt; filt_ia_d = filt_ia
    if _use_dev
        verbose && println("CLM.jl: building device mirror (device_adapt) …")
        inst_d = device_adapt(inst)
        filt_d = device_adapt(filt)
        filt_ia_d = device_adapt(filt_ia)
    end
    # Persistent host view of the device tree, allocated ONCE and refreshed in-place each
    # step (see _refresh_host_tree!) — the host shuttle for the host-only lnd2atm! +
    # history_write_step!. Same Float32 leaves as inst_d so the per-step copy is a plain
    # device->host memcpy (no allocation).
    hview = _use_dev ? Adapt.adapt(Array, inst_d) : inst

    # ========================================================================
    # Phase 4: Time integration loop
    # ========================================================================
    total_steps = Int(ceil(Dates.value(end_date - start_date) / (1000 * dtime)))
    step_count = 0

    verbose && println("CLM.jl: Starting time loop ($total_steps steps, dtime=$(dtime)s)")

    while tm.current_date < end_date
        # The forcing for a step is the value at the step's START time (verified
        # against the Fortran datm: step 8761 from 2003-01-01 00:00 uses
        # PRECTmms[0], not [1]). Capture it before advancing the clock, which
        # otherwise reads one interval ahead.
        step_start = tm.current_date
        advance_timestep!(tm)
        step_count += 1

        # --- Read and downscale forcings (at step-start time) ---
        # forcing_phase_shift_s shifts the forcing-evaluation time. The Fortran
        # datm evaluates the LINEAR-interp state fields (esp. temperature) one
        # coupling interval behind CLM.jl's step-start read (verified at Kherlen
        # n9500: model 20:00 uses the 19:00 temperature record) — but its coszen
        # SOLAR tracks the true astronomical time, which CLM.jl already matches.
        # A uniform shift=-dtime therefore fixes the temperature phase but breaks
        # the solar/albedo diurnal alignment, and A/B on Kherlen made it worse
        # (67/69 -> 42/69: reflected-SW/albedo blow up). So the default is 0
        # (step-start): faithful for solar/precip, and the residual temperature-
        # phase offset is a night-time single-step floor that doesn't gate the
        # annual scorecard. Knob retained for datm-exact single-step oracle work.
        read_forcing_step!(fr, inst.atm2lnd, step_start + Second(forcing_phase_shift_s), ng, nc;
                           gridcell_latdeg = inst.gridcell.latdeg,
                           gridcell_londeg = inst.gridcell.londeg,
                           dtime = dtime)
        downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)

        # --- Compute orbital parameters ---
        calday = get_curr_calday(tm)
        (declin, eccf) = compute_orbital(calday)
        declinp1 = declin  # simplified: same as current step
        obliqr = ORB_OBLIQR_DEFAULT

        # --- Time flags ---
        doalb = true  # compute albedo every step for simplicity
        # `calday` is already post-advance_timestep! (= step_start + dtime), i.e. the
        # start of the NEXT radiation step — exactly the forcing time this step's
        # albedo is paired with in the next surface_radiation!. Adding another dtime
        # put the albedo coszen one full step ahead of its forcing, so at each
        # sunset the last sunlit forcing step got the next step's coszen<=0 -> night
        # albedo=1.0 fallback, spiking daily-mean reflected SW (dense-canopy FSR
        # +14-20%). Pairing nextsw_cday = calday removes the mismatch.
        nextsw_cday = calday

        (yr, mon, d, tod) = get_curr_date(tm)
        first_step = (tm.nstep == 1)
        beg_day = is_beg_curr_day(tm)
        end_day = is_end_curr_day(tm)
        beg_year = is_beg_curr_year(tm)

        # --- Run driver (host, or on the device mirror) ---
        if _use_dev
            # Hand the freshly-downscaled forcing to the device, run physics on-device,
            # then pull a host view for the host-only lnd2atm! + history_write_step!.
            _sync_arrays_convert!(inst_d.atm2lnd, inst.atm2lnd)
            clm_drv!(config, inst_d, filt_d, filt_ia_d, bounds,
                     doalb, nextsw_cday, declinp1, declin, obliqr,
                     false, false, "", false;
                     nstep=tm.nstep, is_first_step=first_step,
                     is_beg_curr_day=beg_day, is_end_curr_day=end_day,
                     is_beg_curr_year=beg_year, dtime=Float64(dtime),
                     mon=mon, day=d, jday=dayofyear(tm.current_date), secs=tod,
                     year=yr, photosyns=inst_d.photosyns)
            _refresh_host_tree!(hview, inst_d)   # in-place device -> host sync (no per-step alloc)
            lnd2atm!(bounds, hview)
            history_write_step!(hw, hview, tm.current_date; is_end_curr_day=end_day)
            step_probe === nothing || step_probe(hview, tm)
            # Device kernel scratch (KA temporaries) still churns Apple's unified memory;
            # collect periodically so a full-year run reclaims it. The host side no longer
            # allocates per step (in-place refresh above), so this is the only pressure.
            (step_count % 48 == 0) && GC.gc()
        else
            clm_drv!(config, inst, filt, filt_ia, bounds,
                     doalb, nextsw_cday, declinp1, declin, obliqr,
                     false, false, "", false;
                     nstep=tm.nstep, is_first_step=first_step,
                     is_beg_curr_day=beg_day, is_end_curr_day=end_day,
                     is_beg_curr_year=beg_year, dtime=Float64(dtime),
                     mon=mon, day=d, jday=dayofyear(tm.current_date), secs=tod,
                     year=yr, photosyns=inst.photosyns)
            lnd2atm!(bounds, inst)
            history_write_step!(hw, inst, tm.current_date; is_end_curr_day=end_day)
            step_probe === nothing || step_probe(inst, tm)
        end

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

    if overrides !== nothing && !isempty(active_fields)
        validate_overrides_post!(overrides, active_fields)
    end

    return inst
end
