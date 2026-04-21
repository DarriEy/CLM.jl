# ==========================================================================
# Multi-site calibration for CLM.jl
#
# Extends the single-site CalibrationProblem to support:
# - Multiple FLUXNET sites sharing a parameter vector
# - Site-specific forcing and observations
# - Weighted multi-site objective aggregation
#
# Public API:
#   SiteCalibrationProblem   — single site calibration setup
#   MultiSiteCalibration     — multi-site aggregation
#   calibration_objective_multisite — multi-site SSE
#   calibration_gradient_multisite  — multi-site gradient
# ==========================================================================

"""
    SiteCalibrationProblem

Calibration problem for a single FLUXNET site.
Wraps CalibrationProblem with site-specific forcing and observations.
"""
Base.@kwdef struct SiteCalibrationProblem
    site::FluxnetSiteInfo = FluxnetSiteInfo(site_id="synthetic", lat=51.0, lon=-116.0)
    problem::CalibrationProblem
    weight::Float64 = 1.0           # weight in multi-site objective
    forcing_series::Vector{FluxnetForcing} = FluxnetForcing[]  # time series forcing
    obs_series::Vector{FluxnetTarget} = FluxnetTarget[]        # time series observations
end

"""
    MultiSiteCalibration

Aggregates multiple site calibration problems sharing a common parameter vector.

- `sites`: vector of SiteCalibrationProblem
- `params`: shared CalibrationParameter vector
- `aggregation`: :sum (default) or :mean of per-site objectives
"""
Base.@kwdef struct MultiSiteCalibration
    sites::Vector{SiteCalibrationProblem}
    params::Vector{CalibrationParameter}
    aggregation::Symbol = :sum
end

"""
    calibration_objective_multisite(msc::MultiSiteCalibration, theta::AbstractVector{T}) where T<:Real

Compute weighted sum of per-site objectives.
Each site is initialized independently but shares the parameter vector theta.
"""
function calibration_objective_multisite(msc::MultiSiteCalibration, theta::AbstractVector{T}) where T<:Real
    total = zero(T)
    for sp in msc.sites
        site_obj = calibration_objective(sp.problem, theta)
        total = total + T(sp.weight) * site_obj
    end

    if msc.aggregation === :mean && !isempty(msc.sites)
        total = total / T(length(msc.sites))
    end

    return total
end

"""
    calibration_gradient_multisite(msc::MultiSiteCalibration, theta::Vector{Float64})

Compute gradient of multi-site objective using ForwardDiff.
"""
function calibration_gradient_multisite(msc::MultiSiteCalibration, theta::Vector{Float64})
    return ForwardDiff.gradient(th -> calibration_objective_multisite(msc, th), theta)
end

"""
    calibrate_multisite(msc::MultiSiteCalibration;
                        theta0::Vector{Float64} = default_theta(msc.params),
                        maxiter::Int = 50,
                        verbose::Bool = true)

Run gradient descent on the multi-site objective.
"""
function calibrate_multisite(msc::MultiSiteCalibration;
                             theta0::Vector{Float64} = default_theta(msc.params),
                             maxiter::Int = 50,
                             gtol::Float64 = 1e-6,
                             ftol::Float64 = 1e-8,
                             verbose::Bool = true,
                             max_grad::Float64 = 1e3)
    # Create a thin CalibrationProblem wrapper that delegates to multisite
    # We reuse the gradient descent machinery from optimize.jl
    theta = copy(theta0)
    n = length(theta)
    trajectory = Tuple{Int, Float64, Vector{Float64}}[]

    f_prev = calibration_objective_multisite(msc, theta)
    push!(trajectory, (0, f_prev, copy(theta)))

    if verbose
        println("MultiSite: iter=0, objective=$(round(f_prev, sigdigits=6))")
    end

    converged = false
    iter = 0

    for k in 1:maxiter
        iter = k
        grad = calibration_gradient_multisite(msc, theta)
        grad = clamp.(grad, -max_grad, max_grad)
        gnorm = sqrt(sum(g^2 for g in grad))

        if gnorm < gtol
            converged = true
            break
        end

        alpha = min(1.0, 1.0 / max(1.0, gnorm))
        direction = -grad
        c1 = 1e-4; rho = 0.5

        f_new = Inf
        for _ in 1:20
            theta_trial = theta + alpha * direction
            f_new = calibration_objective_multisite(msc, theta_trial)
            if f_new <= f_prev + c1 * alpha * dot(grad, direction)
                theta = theta_trial
                break
            end
            alpha *= rho
        end

        push!(trajectory, (k, f_new, copy(theta)))

        if verbose
            println("MultiSite: iter=$k, obj=$(round(f_new, sigdigits=6)), |grad|=$(round(gnorm, sigdigits=4))")
        end

        if abs(f_prev - f_new) < ftol
            converged = true
            break
        end
        f_prev = f_new
    end

    f_final = calibration_objective_multisite(msc, theta)
    return CalibrationResult(theta, f_final, trajectory, converged, iter)
end
