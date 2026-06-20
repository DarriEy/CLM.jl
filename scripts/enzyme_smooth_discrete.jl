# =============================================================================
# Phase 3 — DISCRETE phases reverse-differentiate UNDER ENZYME when smoothed. The discrete CLM
# constructs (phenology onset/offset triggers, fire thresholds, mortality steps) are written
# with the smooth primitives (smooth_heaviside / smooth_max / smooth_min), but those are EXACT
# for Float64 under the default SMOOTH_MODE=:auto — so a HARD step/max has a ZERO (or one-sided)
# reverse-mode gradient and the discrete phase is INVISIBLE to Enzyme reverse-AD (and thus to
# gradient calibration). With SMOOTH_MODE[]=:always the Float64 methods evaluate the SMOOTH
# function (the _smooth_f64() host override added to smooth_ad.jl), and Enzyme gets a finite,
# FD-consistent gradient. This contrasts a hard vs a smoothed mini phenology ONSET trigger,
# built from the EXACT primitives the real phenology FSM uses (phenology.jl:1814/1823).
#
# GOTCHA: SMOOTH_MODE must be set to :always BEFORE the first Enzyme compile of a given
# function — Enzyme caches the adjoint keyed by function type, not by the runtime global, so
# compiling a function once under :auto poisons it for :always. (Hence hard/smooth use SEPARATE
# functions here.) For real reverse-AD, set :always once at process start.
#   julia +1.10 --project=/tmp/clm_jl10_<id> scripts/enzyme_smooth_discrete.jl
# =============================================================================
using CLM, Enzyme, Printf
const C = CLM
const TFRZ = C.TFRZ
C.SMOOTH_MODE[] = :always           # set up front, before any Enzyme compilation

# HARD onset trigger: literal step + max (what the model is WITHOUT smoothing, and what Enzyme
# reverse-AD sees for Float64 under the default :auto). onset_fdd = heaviside(TFRZ−soilt)·fracday,
# onset_gdd = max(soilt−TFRZ, 0)·fracday ; L = fdd² + gdd².
function onset_hard(soilt::Vector{Float64})
    fracday = 1.0/48.0
    fdd = (TFRZ - soilt[1] >= 0.0 ? 1.0 : 0.0) * fracday
    gdd = max(soilt[1] - TFRZ, 0.0) * fracday
    return fdd^2 + gdd^2
end
# SMOOTH onset trigger: the SAME math via the CLM smooth primitives (sigmoid / LogSumExp).
function onset_smooth(soilt::Vector{Float64})
    fracday = 1.0/48.0
    fdd = C.smooth_heaviside(TFRZ - soilt[1]) * fracday
    gdd = C.smooth_max(soilt[1] - TFRZ, 0.0) * fracday
    return fdd^2 + gdd^2
end

rev(f, x0) = (x=copy(x0); dx=zero(x); Enzyme.autodiff(set_runtime_activity(Reverse), f, Active, Duplicated(x,dx)); dx[1])
fd(f, x0; h=1e-5) = (a=copy(x0); a[1]+=h; b=copy(x0); b[1]-=h; (f(a)-f(b))/(2h))

# Right at the freeze threshold — where the hard step JUMPS (so its smooth gradient is real
# information the hard model throws away).
soilt0 = [TFRZ - 0.01]
println("="^72)
for (name, f, note) in (("onset_HARD  ", onset_hard,   "literal step/max → reverse-AD BLIND (discrete phase invisible)"),
                        ("onset_SMOOTH", onset_smooth, "CLM smooth primitives → reverse-AD SEES the trigger"))
    gr = rev(f, soilt0); gf = fd(f, soilt0)
    rel = abs(gf) > 1e-12 ? abs(gr-gf)/abs(gf) : abs(gr-gf)
    verdict = name == "onset_HARD  " ? (abs(gr) < 1e-12 ? "grad ZERO ✓ (uncalibratable)" : "") :
                                       (rel < 1e-4 ? "PASS ✓ (Enzyme reverse matches FD)" : "FAIL ✗")
    @printf("[%s] %s\n             rev d/d(soilt)=% .6e  FD=% .6e  rel=%.2e  %s\n", name, note, gr, gf, rel, verdict)
end
println("="^72)
println("→ Smoothing (SMOOTH_MODE=:always + the _smooth_f64 Float64 override) makes the discrete")
println("  phenology/fire/mortality triggers reverse-differentiable under Enzyme: the hard model")
println("  hands calibration a zero gradient; the smoothed model exposes a finite, FD-exact one.")
