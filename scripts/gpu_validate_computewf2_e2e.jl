# ==========================================================================
# gpu_validate_computewf2_e2e.jl — Metal parity for compute_wf2!
#
# `compute_wf2!` (src/biogeophys/hydrology_no_drainage.jl) was added by #248 to
# reproduce the CTSM CARRYOVER quirk: HydrologyNoDrainageMod.F90:640-699 computes
# wf (top 0.05 m) then wf2 (top 0.17 m) in ONE routine and does NOT reset
# rwat/swat/rz between the two depth loops, so the top-0.05 m layers are
# accumulated TWICE into wf2. Omitting that double-count biases wf2 by O(0.01) —
# enough to flip it across the `max(wf2,0)` kink in the peatland-fire baf_peatf
# branch (it made BAF_PEATF run 7.7% high until #248).
#
# It is a real KernelAbstractions kernel (`_hydrond_compute_wf2_kernel!`) and had
# NO GPU coverage at all: the last GPU sweep (#242) predates it. A kernel whose
# whole purpose is an accumulation-order subtlety is exactly where a device port
# can silently diverge, so this pins device == host.
#
#   julia --project=scripts scripts/gpu_validate_computewf2_e2e.jl
#
# Metal has no Float64, so the device runs Float32 and the comparison is at the
# Float32 floor (~1e-7), matching the other gpu_validate_* harnesses.
# ==========================================================================

using CLM
using Printf
import Metal
include(joinpath(@__DIR__, "gpu_backends.jl"))

function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(A[i] - B[i]) / (1.0 + max(abs(A[i]), abs(B[i]))))
    end
    return m
end

struct MetalF32 end
CLM.Adapt.adapt_storage(::MetalF32, x::AbstractArray{<:AbstractFloat}) = Metal.MtlArray(Float32.(x))
CLM.Adapt.adapt_storage(::MetalF32, x::AbstractArray{<:Integer})       = Metal.MtlArray(x)
CLM.Adapt.adapt_storage(::MetalF32, x::AbstractArray{Bool})            = Metal.MtlArray(x)

const NC       = 6      # columns (mix of masked-in / masked-out)
const NLEVGRND = 5
const NLEVSNO  = 5
const INNER    = 0.05   # CTSM wf  depth
const OUTER    = 0.17   # CTSM wf2 depth

# CLM 20SL_8.5m top-5 node depths / thicknesses. With these, z+0.5dz <= 0.05 holds
# for layers 1-2 and <= 0.17 for layers 1-4 — so the inner sums really are carried
# into the outer loop and 2 layers are double-counted. If a future layer structure
# made the inner set EMPTY this harness would be vacuous, so we assert it isn't.
const Z_NODE = [0.00710, 0.02792, 0.06226, 0.11886, 0.21222]
const DZ_LAY = [0.01751, 0.02757, 0.04547, 0.07496, 0.12360]

function build()
    z  = zeros(Float64, NC, NLEVSNO + NLEVGRND)
    dz = zeros(Float64, NC, NLEVSNO + NLEVGRND)
    for c in 1:NC, j in 1:NLEVGRND
        z[c,  j + NLEVSNO] = Z_NODE[j]
        dz[c, j + NLEVSNO] = DZ_LAY[j]
    end

    watsat = fill(0.0, NC, NLEVGRND)
    sucsat = fill(0.0, NC, NLEVGRND)
    bsw    = fill(0.0, NC, NLEVGRND)
    h2osoi_vol = fill(0.0, NC, NLEVGRND)
    for c in 1:NC, j in 1:NLEVGRND
        # vary by column so every column has a distinct answer (a kernel that
        # wrote a single broadcast value would still pass a uniform fixture)
        watsat[c, j] = 0.40 + 0.02 * c + 0.005 * j
        sucsat[c, j] = 100.0 + 10.0 * c + 5.0 * j
        bsw[c, j]    = 4.0 + 0.3 * c + 0.1 * j
        # span wet -> dry, incl. a column below watdry (drives wf2 negative, the
        # regime that matters for the max(wf2,0) kink baf_peatf sits on)
        h2osoi_vol[c, j] = c == NC ? 0.02 : (0.35 - 0.04 * c + 0.01 * j)
    end

    mask = fill(true, NC)
    mask[2] = false          # masked-out column must be left untouched
    return (; z, dz, watsat, sucsat, bsw, h2osoi_vol, mask)
end

function main()
    s = build()
    bounds = 1:NC

    # ---- sanity: the fixture must actually exercise the carryover ----
    n_inner = count(j -> Z_NODE[j] + 0.5 * DZ_LAY[j] <= INNER, 1:NLEVGRND)
    n_outer = count(j -> Z_NODE[j] + 0.5 * DZ_LAY[j] <= OUTER, 1:NLEVGRND)
    @printf("layers within inner %.2fm: %d | within outer %.2fm: %d\n", INNER, n_inner, OUTER, n_outer)
    if n_inner == 0 || n_outer <= n_inner
        println("[FAIL] fixture does not exercise the double-count carryover — vacuous")
        return 1
    end

    # ---- CPU reference (Float64) ----
    wf2_cpu = fill(-999.0, NC)
    CLM.compute_wf2!(wf2_cpu, s.h2osoi_vol, s.watsat, s.sucsat, s.bsw, s.z, s.dz,
                     s.mask, bounds, NLEVGRND, NLEVSNO, INNER, OUTER)

    # ---- device (Metal, Float32) ----
    A = MetalF32()
    wf2_dev = CLM.Adapt.adapt(A, fill(-999.0, NC))
    CLM.compute_wf2!(wf2_dev,
                     CLM.Adapt.adapt(A, s.h2osoi_vol), CLM.Adapt.adapt(A, s.watsat),
                     CLM.Adapt.adapt(A, s.sucsat),     CLM.Adapt.adapt(A, s.bsw),
                     CLM.Adapt.adapt(A, s.z),          CLM.Adapt.adapt(A, s.dz),
                     CLM.Adapt.adapt(A, s.mask), bounds, NLEVGRND, NLEVSNO, INNER, OUTER)
    wf2_host = Array(wf2_dev)

    ok = true

    # masked-out column must be untouched on BOTH sides
    if wf2_cpu[2] != -999.0 || wf2_host[2] != -999.0f0
        @printf("[FAIL] masked-out column written: cpu=%g dev=%g\n", wf2_cpu[2], wf2_host[2])
        ok = false
    else
        println("[PASS] masked-out column untouched")
    end

    # the active columns must be finite, non-degenerate, and device == host
    active = [c for c in 1:NC if s.mask[c]]
    if !all(isfinite, wf2_cpu[active])
        println("[FAIL] CPU wf2 non-finite on an active column"); ok = false
    end
    if length(unique(round.(wf2_cpu[active], digits=9))) < length(active)
        println("[FAIL] wf2 identical across columns — fixture/kernel not column-varying"); ok = false
    end
    if !any(<(0.0), wf2_cpu[active])
        println("[FAIL] no column landed below the max(wf2,0) kink — regime untested"); ok = false
    else
        println("[PASS] a column exercises the negative-wf2 (kink) regime")
    end

    r = reldiff(wf2_cpu[active], wf2_host[active])
    tol = 1e-6   # Float32 device vs Float64 host
    @printf("%s wf2 device-vs-host  rel = %.3e  (tol %.0e)\n", r <= tol ? "[PASS]" : "[FAIL]", r, tol)
    r <= tol || (ok = false)

    for c in active
        @printf("    col %d: cpu=% .8f  dev=% .8f\n", c, wf2_cpu[c], wf2_host[c])
    end

    println(ok ? "\ncompute_wf2! Metal parity: PASS" : "\ncompute_wf2! Metal parity: FAIL")
    return ok ? 0 : 1
end

exit(main())
