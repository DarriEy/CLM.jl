# ==========================================================================
# Enzyme.jl reverse-mode AD wrapper for CLM.jl calibration
#
# Provides calibration_gradient_enzyme() as a drop-in alternative to
# calibration_gradient() which uses ForwardDiff forward-mode.
# Enzyme reverse-mode is O(1) in number of parameters vs O(n) for ForwardDiff.
#
# Requires: `using Enzyme` before calling these functions.
# ==========================================================================

"""
    calibration_gradient_enzyme(prob::CalibrationProblem, theta::Vector{Float64})

Compute ∇_θ objective using Enzyme reverse-mode AD.
Falls back to ForwardDiff if Enzyme is not loaded.

This is more efficient than ForwardDiff when nparams > 1, since
reverse-mode computes the full gradient in a single backward pass.
"""
function calibration_gradient_enzyme(prob::CalibrationProblem, theta::Vector{Float64})
    if !isdefined(Main, :Enzyme)
        @warn "Enzyme.jl not loaded, falling back to ForwardDiff" maxlog=1
        return calibration_gradient(prob, theta)
    end

    Enz = Main.Enzyme

    n = length(theta)
    grad = zeros(n)

    # Use Enzyme's gradient API
    # We differentiate calibration_objective w.r.t. theta
    function obj_wrapper(th::Vector{Float64})
        return calibration_objective(prob, th)
    end

    dtheta = zeros(n)
    Enz.autodiff(Enz.Reverse, obj_wrapper, Enz.Active, Enz.Duplicated(copy(theta), dtheta))
    return dtheta
end
