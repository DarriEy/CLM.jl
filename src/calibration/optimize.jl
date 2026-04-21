# ==========================================================================
# Optimization wrapper for CLM.jl calibration
#
# Provides a built-in gradient descent with backtracking line search
# (no external dependencies) and optional Optim.jl integration.
# ==========================================================================

"""
    CalibrationResult

Result of a calibration optimization run.

- `theta`: optimized parameter vector
- `objective`: final objective value
- `trajectory`: vector of (iteration, objective, theta) tuples
- `converged`: whether optimization converged
- `iterations`: number of iterations taken
"""
struct CalibrationResult
    theta::Vector{Float64}
    objective::Float64
    trajectory::Vector{Tuple{Int, Float64, Vector{Float64}}}
    converged::Bool
    iterations::Int
end

"""
    calibrate(prob::CalibrationProblem;
              theta0::Vector{Float64} = default_theta(prob.params),
              method::Symbol = :gradient_descent,
              maxiter::Int = 50,
              gtol::Float64 = 1e-6,
              ftol::Float64 = 1e-8,
              verbose::Bool = true)

Run calibration optimization.

Methods:
- `:gradient_descent` — built-in gradient descent with backtracking line search
- `:lbfgsb` — L-BFGS-B via Optim.jl (requires `using Optim`)
"""
function calibrate(prob::CalibrationProblem;
                   theta0::Vector{Float64} = default_theta(prob.params),
                   method::Symbol = :gradient_descent,
                   maxiter::Int = 50,
                   gtol::Float64 = 1e-6,
                   ftol::Float64 = 1e-8,
                   verbose::Bool = true,
                   max_grad::Float64 = 1e3)

    if method === :gradient_descent
        return _gradient_descent(prob, theta0; maxiter, gtol, ftol, verbose, max_grad)
    elseif method === :lbfgsb
        return _optim_lbfgsb(prob, theta0; maxiter, gtol, ftol, verbose)
    else
        error("Unknown optimization method: $method. Use :gradient_descent or :lbfgsb")
    end
end

"""
Built-in gradient descent with backtracking line search (Armijo condition).
"""
function _gradient_descent(prob::CalibrationProblem, theta0::Vector{Float64};
                           maxiter::Int = 50, gtol::Float64 = 1e-6,
                           ftol::Float64 = 1e-8, verbose::Bool = true,
                           max_grad::Float64 = 1e3)
    theta = copy(theta0)
    n = length(theta)
    trajectory = Tuple{Int, Float64, Vector{Float64}}[]

    f_prev = calibration_objective(prob, theta)
    push!(trajectory, (0, f_prev, copy(theta)))

    if verbose
        @info "Calibration: iter=0, objective=$(round(f_prev, sigdigits=6))"
    end

    converged = false
    iter = 0

    for k in 1:maxiter
        iter = k
        grad = calibration_gradient(prob, theta)
        gnorm = sqrt(sum(g^2 for g in grad))

        if gnorm < gtol
            converged = true
            if verbose
                @info "Calibration: converged (gradient norm $(round(gnorm, sigdigits=4)) < $gtol)"
            end
            break
        end

        # Gradient clipping
        grad = clamp.(grad, -max_grad, max_grad)
        gnorm = sqrt(sum(g^2 for g in grad))

        # Backtracking line search (Armijo) with normalized initial step
        alpha = min(1.0, 1.0 / max(1.0, gnorm))
        c1 = 1e-4
        rho = 0.5
        direction = -grad

        f_new = Inf
        for _ in 1:20  # max line search iterations
            theta_trial = theta + alpha * direction
            f_new = calibration_objective(prob, theta_trial)
            if f_new <= f_prev + c1 * alpha * dot(grad, direction)
                theta = theta_trial
                break
            end
            alpha *= rho
        end

        push!(trajectory, (k, f_new, copy(theta)))

        if verbose
            @info "Calibration: iter=$k, objective=$(round(f_new, sigdigits=6)), " *
                  "|grad|=$(round(gnorm, sigdigits=4)), step=$(round(alpha, sigdigits=4))"
        end

        if abs(f_prev - f_new) < ftol
            converged = true
            if verbose
                @info "Calibration: converged (objective change < $ftol)"
            end
            break
        end

        f_prev = f_new
    end

    f_final = calibration_objective(prob, theta)
    return CalibrationResult(theta, f_final, trajectory, converged, iter)
end

"""
L-BFGS-B optimization via Optim.jl (requires Optim to be loaded).
"""
function _optim_lbfgsb(prob::CalibrationProblem, theta0::Vector{Float64};
                       maxiter::Int = 50, gtol::Float64 = 1e-6,
                       ftol::Float64 = 1e-8, verbose::Bool = true)
    # Check if Optim is available
    if !isdefined(Main, :Optim)
        error("Optim.jl not loaded. Run `using Optim` first, or use method=:gradient_descent")
    end

    Optim = Main.Optim

    # Set up bounds
    n = length(prob.params)
    lower = fill(-Inf, n)
    upper = fill(Inf, n)

    f = th -> calibration_objective(prob, th)
    g! = (G, th) -> begin
        grad = calibration_gradient(prob, th)
        G .= grad
    end

    result = Optim.optimize(f, g!, lower, upper, theta0,
                            Optim.Fminbox(Optim.LBFGS()),
                            Optim.Options(iterations=maxiter, g_tol=gtol,
                                          f_tol=ftol, show_trace=verbose))

    trajectory = [(0, Optim.minimum(result), copy(Optim.minimizer(result)))]
    return CalibrationResult(
        Optim.minimizer(result),
        Optim.minimum(result),
        trajectory,
        Optim.converged(result),
        Optim.iterations(result)
    )
end
