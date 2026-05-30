# PoC for the convergence-mask kernel pattern — the GPU-friendly shape for CLM's
# iterative per-element solvers (canopy_fluxes Newton, photosynthesis Brent, ...).
#
# A `while not converged` loop diverges across GPU threads (each element needs a
# different iteration count). The pattern: run a FIXED number of iterations and
# carry a per-thread `done` flag — once an element converges it stops updating but
# the uniform loop keeps running, so all threads share the same control flow.
#
# Demonstrated on a real physics root-find: invert CLM.qsat per column (find T
# such that qsat(T, p) == q_target) via Newton, using the analytic dqs/dT that
# qsat already returns. Validated against a scalar while-loop reference.
#
#   julia --project=. scripts/ka_convergence_kernel.jl [ncolumns]

using CLM
using KernelAbstractions
using Printf

# --- scalar reference: per-column Newton with an early-exit while-style loop ---
function invert_qsat_scalar(T0, p, qtarget; max_iter=100, tol=1e-13)
    T = copy(T0)
    for c in eachindex(T)
        Tc = T[c]
        for _ in 1:max_iter
            (qs, _, dqsdT, _) = CLM.qsat(Tc, p[c])
            r = qs - qtarget[c]
            abs(r) < tol && break
            Tc = Tc - r / dqsdT
        end
        T[c] = Tc
    end
    return T
end

# --- convergence-mask kernel: fixed iterations + per-thread `done` flag ---
@kernel function _invert_qsat_kernel!(T, @Const(p), @Const(qtarget), max_iter::Int, tol::Float64)
    c = @index(Global)
    @inbounds begin
        Tc = T[c]
        done = false
        for _ in 1:max_iter
            if !done
                (qs, _, dqsdT, _) = CLM.qsat(Tc, p[c])
                r = qs - qtarget[c]
                if abs(r) < tol
                    done = true            # converged: stop updating, loop continues uniformly
                else
                    Tc = Tc - r / dqsdT
                end
            end
        end
        T[c] = Tc
    end
end

function invert_qsat_kernel!(T, p, qtarget; max_iter=100, tol=1e-13)
    backend = KernelAbstractions.get_backend(T)
    _invert_qsat_kernel!(backend)(T, p, qtarget, max_iter, tol; ndrange=length(T))
    KernelAbstractions.synchronize(backend)
    return T
end

# ---------------- validation ----------------
n = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1000
p = fill(85000.0, n)
# targets spanning a range of saturation humidities (so elements need different
# iteration counts — the whole point of the convergence mask)
Ttrue = collect(range(260.0, 300.0, length=n))
qtarget = [CLM.qsat(Ttrue[c], p[c])[1] for c in 1:n]

T0 = fill(280.0, n)                          # common initial guess
Tref = invert_qsat_scalar(copy(T0), p, qtarget)
Tker = invert_qsat_kernel!(copy(T0), p, qtarget)

err_vs_scalar = maximum(abs.(Tker .- Tref))
err_vs_true   = maximum(abs.(Tker .- Ttrue))
@printf("columns = %d\n", n)
@printf("max |kernel - scalar Newton| = %.3e\n", err_vs_scalar)
@printf("max |kernel - true T|        = %.3e\n", err_vs_true)
println((err_vs_scalar < 1e-10 && err_vs_true < 1e-6) ?
    "\nCONVERGENCE-MASK KERNEL PASSED ✓ (matches scalar; GPU-friendly fixed-iteration shape)" :
    "\nFAILED ✗")
