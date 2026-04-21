# ==========================================================================
# Calibration framework for CLM.jl
#
# Provides differentiable calibration using ForwardDiff through clm_drv!.
# Parameters are injected as scaling factors: calibrated_value = default * theta[i].
# This avoids Dual-ifying global parameter structs.
#
# Public API:
#   CalibrationParameter  — parameter descriptor
#   CalibrationTarget     — observation/output to match
#   CalibrationProblem    — full calibration setup
#   calibration_objective — compute weighted SSE
#   calibration_gradient  — ForwardDiff.gradient of objective
# ==========================================================================

import ForwardDiff

# =====================================================================
# CalibrationParameter — describes one tunable parameter
# =====================================================================

"""
    CalibrationParameter

Descriptor for a single calibration parameter.

- `name`: human-readable identifier
- `default_value`: baseline value (Float64)
- `apply!`: function `(inst, val) -> nothing` that injects the parameter value into CLM state
- `bounds`: (lo, hi) bounds in transformed space
- `transform`: :identity, :log, or :logit — maps unconstrained θ to physical space
"""
struct CalibrationParameter
    name::String
    default_value::Float64
    apply!::Function        # (inst, val) -> nothing
    bounds::Tuple{Float64, Float64}
    transform::Symbol       # :identity, :log, :logit
end

"""
    transform_param(cp::CalibrationParameter, theta_raw)

Transform unconstrained θ to physical parameter space.
- `:identity` — val = default * theta_raw  (theta_raw ≈ 1.0 at default)
- `:log`      — val = default * exp(theta_raw)  (theta_raw ≈ 0.0 at default)
- `:logit`    — val = lo + (hi - lo) * sigmoid(theta_raw)
"""
function transform_param(cp::CalibrationParameter, theta_raw)
    if cp.transform === :identity
        return cp.default_value * theta_raw
    elseif cp.transform === :log
        return cp.default_value * exp(theta_raw)
    elseif cp.transform === :logit
        lo, hi = cp.bounds
        s = one(theta_raw) / (one(theta_raw) + exp(-theta_raw))
        return lo + (hi - lo) * s
    else
        error("Unknown transform: $(cp.transform)")
    end
end

"""
    default_theta(cp::CalibrationParameter)

Return the θ value that reproduces the default parameter value.
"""
function default_theta(cp::CalibrationParameter)
    if cp.transform === :identity
        return 1.0
    elseif cp.transform === :log
        return 0.0
    elseif cp.transform === :logit
        lo, hi = cp.bounds
        # Solve: lo + (hi-lo)*σ(θ) = default → σ(θ) = (default-lo)/(hi-lo) → θ = logit(...)
        p = (cp.default_value - lo) / (hi - lo)
        p = clamp(p, 1e-6, 1.0 - 1e-6)
        return log(p / (1.0 - p))
    else
        error("Unknown transform: $(cp.transform)")
    end
end

# =====================================================================
# CalibrationTarget — an observed quantity to match
# =====================================================================

"""
    CalibrationTarget

Descriptor for a calibration target (observation).

- `name`: human-readable identifier
- `extract`: function `(inst) -> value` that extracts the model output
- `observed`: target value
- `weight`: weight in the SSE objective
"""
struct CalibrationTarget
    name::String
    extract::Function       # (inst) -> scalar
    observed::Float64
    weight::Float64
end

# =====================================================================
# CalibrationProblem — full calibration setup
# =====================================================================

"""
    CalibrationProblem

Bundles parameters, targets, and simulation configuration for calibration.

- `params`: vector of CalibrationParameter
- `targets`: vector of CalibrationTarget
- `fsurdat`, `paramfile`: input file paths
- `n_warmup`: number of warmup timesteps before the evaluation step
- `T_forc`: forcing temperature [K]
- `dtime`: timestep [s]
- `calday`: calendar day for orbital calculations
"""
Base.@kwdef struct CalibrationProblem
    params::Vector{CalibrationParameter}
    targets::Vector{CalibrationTarget}
    fsurdat::String
    paramfile::String
    n_warmup::Int = 3
    n_eval::Int = 1                    # number of evaluation timesteps
    eval_aggregation::Symbol = :last   # :last or :mean
    T_forc::Float64 = 285.0
    dtime::Float64 = 1800.0
    calday::Float64 = 172.5             # noon of summer solstice (max sunlight)
    month::Int = 6                     # month for interp_monthly_veg! (June)
    day::Int = 21                      # day for interp_monthly_veg!
    use_smooth_fd::Bool = true         # Float64 FD evals use smooth path for FD/AD consistency
end

# =====================================================================
# Core calibration functions
# =====================================================================

"""
    run_clm_with_params(prob::CalibrationProblem, theta::AbstractVector{T}) where T<:Real

Run CLM with parameter vector θ and return the CLMInstances after simulation.
When T is a Dual type, derivatives flow through the physics.
"""
function run_clm_with_params(prob::CalibrationProblem, theta::AbstractVector{T}) where T<:Real

    # When use_smooth_fd is true, force smooth evaluation for ALL types
    # so that Float64 (FD) and Dual (AD) evaluate the same mathematical function.
    # This must be set here (not just in calibration_objective) so that observation
    # generation and objective evaluation use the same code path.
    smooth_override = prob.use_smooth_fd && T === Float64
    if smooth_override
        old_mode = SMOOTH_MODE[]
        SMOOTH_MODE[] = :always
    end

    try  # finally block restores SMOOTH_MODE

    # Step 1: Initialize Float64 state
    (inst_f64, bounds, filt, tm) = clm_initialize!(;
        fsurdat=prob.fsurdat, paramfile=prob.paramfile)

    ng = bounds.endg
    nc = bounds.endc
    config = CLMDriverConfig()
    filt_ia = clump_filter_inactive_and_active
    (declin, eccf) = compute_orbital(prob.calday)
    nextsw_cday = prob.calday + prob.dtime / SECSPDAY

    # Initialize soil moisture (cold-start leaves h2osoi_liq = 0 → NaN in physics)
    _init_calib_soil_moisture!(inst_f64, bounds)

    # Set up forcing
    _setup_calib_forcing!(inst_f64.atm2lnd, prob.T_forc, ng)
    downscale_forcings!(bounds, inst_f64.atm2lnd, inst_f64.column,
                        inst_f64.landunit, inst_f64.topo)

    # Initialize timwt so satellite_phenology! computes non-zero LAI
    interp_monthly_veg!(inst_f64.satellite_phenology; kmo=prob.month, kda=prob.day)

    # Run satellite_phenology! to set LAI and frac_veg_nosno_alb
    cs = inst_f64.canopystate
    wdb = inst_f64.water.waterdiagnosticbulk_inst
    pch = inst_f64.patch
    satellite_phenology!(inst_f64.satellite_phenology, cs, wdb, pch,
                         filt.nolakep, bounds.begp:bounds.endp)

    # Copy frac_veg_nosno_alb to frac_veg_nosno and rebuild exposedvegp filter
    for p in bounds.begp:bounds.endp
        cs.frac_veg_nosno_patch[p] = cs.frac_veg_nosno_alb_patch[p]
    end
    set_exposedvegp_filter!(filt, bounds, cs.frac_veg_nosno_patch)

    # Warmup steps (Float64 only)
    for n in 1:prob.n_warmup
        clm_drv!(config, inst_f64, filt, filt_ia, bounds,
            true, nextsw_cday, declin, declin, ORB_OBLIQR_DEFAULT,
            false, false, "", false;
            nstep=n, is_first_step=(n==1), is_beg_curr_day=(n==1),
            dtime=prob.dtime, mon=1, day=1,
            photosyns=inst_f64.photosyns)
    end

    n_eval = prob.n_eval

    if T === Float64
        # Pure Float64 path — inject parameters and run
        for (i, cp) in enumerate(prob.params)
            val = transform_param(cp, theta[i])
            cp.apply!(inst_f64, val)
        end

        for k in 1:n_eval
            clm_drv!(config, inst_f64, filt, filt_ia, bounds,
                true, nextsw_cday, declin, declin, ORB_OBLIQR_DEFAULT,
                false, false, "", false;
                nstep=prob.n_warmup + k, is_first_step=false,
                dtime=prob.dtime, mon=1, day=1,
                photosyns=inst_f64.photosyns)
        end

        return inst_f64
    else
        # Dual path — create Dual-typed copy for AD
        inst_d = _make_dual_instances(inst_f64, T)

        # Inject Dual-typed parameters
        for (i, cp) in enumerate(prob.params)
            val = transform_param(cp, theta[i])
            cp.apply!(inst_d, val)
        end

        # Set up forcing for Dual instances
        _setup_calib_forcing!(inst_d.atm2lnd, prob.T_forc, ng)
        downscale_forcings!(bounds, inst_d.atm2lnd, inst_d.column,
                            inst_d.landunit, inst_d.topo)

        for k in 1:n_eval
            clm_drv!(config, inst_d, filt, filt_ia, bounds,
                true, nextsw_cday, declin, declin, ORB_OBLIQR_DEFAULT,
                false, false, "", false;
                nstep=prob.n_warmup + k, is_first_step=false,
                dtime=prob.dtime, mon=1, day=1,
                photosyns=inst_d.photosyns)
        end

        return inst_d
    end

    finally
        if smooth_override
            SMOOTH_MODE[] = old_mode
        end
    end  # try/finally
end

"""
    calibration_objective(prob::CalibrationProblem, theta::AbstractVector{T}) where T<:Real

Run CLM with θ and compute weighted sum of squared errors vs targets.
"""
function calibration_objective(prob::CalibrationProblem, theta::AbstractVector{T}) where T<:Real
    # Smooth mode override is handled by run_clm_with_params (so observation
    # generation and objective evaluation use the same code path).
    inst = run_clm_with_params(prob, theta)

    sse = zero(T)
    for target in prob.targets
        model_val = target.extract(inst)
        residual = model_val - T(target.observed)
        sse = sse + T(target.weight) * residual * residual
    end
    return sse
end

"""
    calibration_gradient(prob::CalibrationProblem, theta::Vector{Float64})

Compute ∇_θ objective using ForwardDiff. Returns gradient vector.
"""
function calibration_gradient(prob::CalibrationProblem, theta::Vector{Float64})
    return ForwardDiff.gradient(th -> calibration_objective(prob, th), theta)
end

# =====================================================================
# Internal helpers
# =====================================================================

"""Set up standard forcing fields for calibration runs."""
function _setup_calib_forcing!(a2l, T, ng)
    for g in 1:ng
        a2l.forc_t_not_downscaled_grc[g] = T
        a2l.forc_pbot_not_downscaled_grc[g] = 85000.0
        a2l.forc_th_not_downscaled_grc[g] = T * (100000.0 / 85000.0)^(RAIR / CPAIR)
        a2l.forc_rho_not_downscaled_grc[g] = 85000.0 / (RAIR * T)
        a2l.forc_lwrad_not_downscaled_grc[g] = 250.0
        a2l.forc_vp_grc[g] = 800.0
        a2l.forc_hgt_grc[g] = 30.0
        a2l.forc_hgt_u_grc[g] = 30.0
        a2l.forc_hgt_t_grc[g] = 30.0
        a2l.forc_hgt_q_grc[g] = 30.0
        a2l.forc_topo_grc[g] = 0.0
        a2l.forc_wind_grc[g] = 3.0
        a2l.forc_u_grc[g] = 3.0
        a2l.forc_v_grc[g] = 0.0
        for b in 1:NUMRAD
            a2l.forc_solad_not_downscaled_grc[g, b] = 150.0
            a2l.forc_solai_grc[g, b] = 75.0
        end
        a2l.forc_solar_not_downscaled_grc[g] = 450.0
        a2l.forc_rain_not_downscaled_grc[g] = 1e-5
        a2l.forc_snow_not_downscaled_grc[g] = 0.0001
        # CO2 and O2 partial pressures (required for photosynthesis)
        a2l.forc_pco2_grc[g] = 40.0     # ~400 ppm CO2 at 100kPa
        a2l.forc_po2_grc[g] = 20900.0   # ~209000 ppm O2 at 100kPa
    end
end

"""
Create a Dual-typed copy of CLMInstances from a Float64 instance.
Reuses the pattern from test_ad_e2e.jl Level 8.
"""
function _make_dual_instances(inst_f64::CLMInstances, ::Type{D}) where D
    inst_d = CLMInstances()

    # Copy all parameterized structs as Dual-typed copies
    for name in fieldnames(CLMInstances)
        name === :water && continue
        name === :surfdata && continue
        src = getfield(inst_f64, name)
        setfield!(inst_d, name, _calib_dual_copy(src, D))
    end
    inst_d.surfdata = inst_f64.surfdata

    # Copy overrides as Dual-typed
    inst_d.overrides = CalibrationOverrides{D}()

    # Copy integer/bool fields that aren't auto-converted
    if hasproperty(inst_d.photosyns, :stomatalcond_mtd)
        inst_d.photosyns.stomatalcond_mtd = inst_f64.photosyns.stomatalcond_mtd
    end

    # Handle WaterData facade
    water_d = WaterData()
    for name in fieldnames(WaterData)
        sv = getfield(inst_f64.water, name)
        try
            setfield!(water_d, name, sv)
        catch e
            @debug "Skipping WaterData field $name during dual copy: $e"
        end
    end

    wsb_d = _calib_dual_copy(inst_f64.water.waterstatebulk_inst, D)
    ws_d  = _calib_dual_copy(inst_f64.water.waterstatebulk_inst.ws, D)
    wsb_d.ws = ws_d
    water_d.waterstatebulk_inst = wsb_d

    wfb_d = _calib_dual_copy(inst_f64.water.waterfluxbulk_inst, D)
    wf_d  = _calib_dual_copy(inst_f64.water.waterfluxbulk_inst.wf, D)
    wfb_d.wf = wf_d
    water_d.waterfluxbulk_inst = wfb_d

    water_d.waterdiagnosticbulk_inst = _calib_dual_copy(
        inst_f64.water.waterdiagnosticbulk_inst, D)
    water_d.waterbalancebulk_inst = _calib_dual_copy(
        inst_f64.water.waterbalancebulk_inst, D)

    if !isempty(water_d.bulk_and_tracers)
        bt = water_d.bulk_and_tracers[water_d.i_bulk]
        bt.waterflux = wfb_d
        bt.waterstate = wsb_d
        bt.waterdiagnostic = water_d.waterdiagnosticbulk_inst
        bt.waterbalance = water_d.waterbalancebulk_inst
    end

    inst_d.water = water_d
    return inst_d
end

"""
Make a Dual-typed copy of a parameterized struct.
Non-parameterized structs returned by reference.
"""
function _calib_dual_copy(src, ::Type{D}) where D
    T = typeof(src)
    wrapper = T.name.wrapper
    if wrapper === T
        return src
    end
    dst = wrapper{D}()
    for name in fieldnames(T)
        sv = getfield(src, name)
        try
            if sv isa Array{Float64}
                setfield!(dst, name, D.(sv))
            elseif sv isa Float64
                setfield!(dst, name, D(sv))
            else
                setfield!(dst, name, sv)
            end
        catch e
            e isa TypeError || rethrow()
        end
    end
    return dst
end

"""
Initialize soil/snow/temperature state for calibration runs.
Cold-start leaves many fields at default values that cause NaN in physics.
"""
function _init_calib_soil_moisture!(inst, bounds)
    nlevsno = varpar.nlevsno
    nlevsoi = varpar.nlevsoi
    nlevgrnd = varpar.nlevgrnd
    ws = inst.water.waterstatebulk_inst.ws
    watsat = inst.soilstate.watsat_col
    col = inst.column
    temp = inst.temperature

    for c in bounds.begc:bounds.endc
        # Initialize soil moisture to 50% of porosity
        for j in 1:nlevsoi
            jj = j + nlevsno
            dz_j = col.dz[c, jj]
            vol_frac = 0.5 * watsat[c, j]
            ws.h2osoi_liq_col[c, jj] = vol_frac * dz_j * 1000.0
            ws.h2osoi_ice_col[c, jj] = 0.0
        end
        # Bedrock layers below nlevsoi
        for j in (nlevsoi+1):nlevgrnd
            jj = j + nlevsno
            ws.h2osoi_liq_col[c, jj] = 0.0
            ws.h2osoi_ice_col[c, jj] = 0.0
        end
        # Reset surface water and snow
        ws.h2osfc_col[c] = 0.0
        ws.h2osno_no_layers_col[c] = 0.0

        # Set ground temperature to something reasonable (will be overwritten by forcing)
        temp.t_grnd_col[c] = 280.0
        for j in 1:nlevgrnd
            jj = j + nlevsno
            temp.t_soisno_col[c, jj] = 280.0
        end
    end

    # Reset snow diagnostic
    wdb = inst.water.waterdiagnosticbulk_inst
    for c in bounds.begc:bounds.endc
        wdb.snow_depth_col[c] = 0.0
        wdb.frac_sno_col[c] = 0.0
    end

    return nothing
end
