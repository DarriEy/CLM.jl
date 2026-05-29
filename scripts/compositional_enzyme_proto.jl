# Prototype: compositional (checkpointed) reverse-mode over a sequence of
# in-place "phase" functions, with a SEPARATE Enzyme call per phase.
#
# This is the orchestration pattern for differentiating clm_drv! without
# compiling the whole monolith at once (which segfaults). For state s and
# phases f1;f2;f3 producing scalar L = g(s):
#   forward:  checkpoint the INPUT state to each phase
#   reverse:  for k = N..1, run Enzyme.Reverse on phase k alone with the
#             checkpointed input and the running adjoint shadow s̄. Enzyme
#             transforms s̄ from the phase's output-adjoint to its input-adjoint.
#
# Validates: chained-per-phase gradient == one-shot Enzyme == ForwardDiff.

using Enzyme, ForwardDiff, LinearAlgebra

# --- three nonlinear in-place "phases" on a state vector (CLM-like mutation) ---
function phase1!(s)
    n = length(s)
    @inbounds for i in 1:n
        s[i] = tanh(s[i]) + 0.1 * s[i]^2
    end
    return nothing
end
function phase2!(s)
    n = length(s)
    tmp = copy(s)
    @inbounds for i in 1:n
        j = mod1(i + 1, n)
        s[i] = tmp[i] * 0.7 + 0.3 * tmp[j] + sin(tmp[i])
    end
    return nothing
end
function phase3!(s)
    n = length(s)
    @inbounds for i in 1:n
        s[i] = exp(-abs(s[i]) * 0.5) + s[i]
    end
    return nothing
end
const PHASES = (phase1!, phase2!, phase3!)

# loss = sum of squares of final state
loss_of(s) = sum(abs2, s)

# full composition (for reference gradients)
function full!(s)
    phase1!(s); phase2!(s); phase3!(s)
    return nothing
end
function full_loss(s0)
    s = copy(s0); full!(s); return loss_of(s)
end

s0 = collect(range(-1.0, 1.0, length=6))

# --- reference 1: ForwardDiff ---
g_fd = ForwardDiff.gradient(full_loss, s0)

# --- reference 2: one-shot Enzyme (small enough here) ---
g_enz_full = zeros(length(s0))
let s = copy(s0), ds = zeros(length(s0))
    Enzyme.autodiff(Enzyme.Reverse, (s) -> (full!(s); nothing), Enzyme.Const,
                    Enzyme.Duplicated(s, ds))
    # ds is gradient of... nothing (full! returns nothing). Need loss seed.
end
# proper one-shot: differentiate full_loss
g_enz_full = Enzyme.gradient(Enzyme.Reverse, full_loss, s0)[1]

# --- compositional chained reverse-mode ---
function compositional_grad(s0, phases)
    # forward with input checkpoints
    s = copy(s0)
    checkpoints = Vector{Vector{Float64}}(undef, length(phases))
    for (k, f) in enumerate(phases)
        checkpoints[k] = copy(s)   # input to phase k
        f(s)
    end
    # seed: ∂loss/∂s_final  (loss = sum s^2  ->  2s)
    s̄ = 2.0 .* s
    # reverse, phase by phase, separate Enzyme call each
    for k in length(phases):-1:1
        s_work = copy(checkpoints[k])
        Enzyme.autodiff(Enzyme.Reverse, (x) -> (phases[k](x); nothing),
                        Enzyme.Const, Enzyme.Duplicated(s_work, s̄))
        # s̄ now holds the input-adjoint of phase k (= output-adjoint of phase k-1)
    end
    return s̄
end

g_comp = compositional_grad(s0, PHASES)

println("ForwardDiff grad   = ", round.(g_fd, digits=6))
println("Enzyme one-shot    = ", round.(g_enz_full, digits=6))
println("Compositional grad = ", round.(g_comp, digits=6))
err_fd   = maximum(abs.(g_comp .- g_fd))
err_enz  = maximum(abs.(g_comp .- g_enz_full))
println("max|comp - FD|      = ", err_fd)
println("max|comp - 1shot|   = ", err_enz)
println((err_fd < 1e-8 && err_enz < 1e-8) ? "\nCOMPOSITIONAL PATTERN PASSED ✓" : "\nFAILED ✗")
