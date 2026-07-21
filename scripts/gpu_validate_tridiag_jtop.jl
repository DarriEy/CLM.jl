# ==========================================================================
# gpu_validate_tridiag_jtop.jl — device parity for `tridiagonal_multi!` on the
# RAGGED jtop > 1 path, against an independent dense oracle.
#
# WHY THIS EXISTS
# ---------------
# #263 (CUDA bring-up) found `tridiagonal_multi!` was WRONG ON THE CPU at -O2 for
# nlevs >= 5: the forward sweep stored cp/dp and re-read them on the next
# iteration, and LLVM reordered that store->load recurrence. The commit message
# records that it "corrupt[ed] the jtop > 1 path outright". The device was correct
# all along. The bug survived for months because the `multi!` unit tests used
# nlevs = 3 and 4 (test/test_tridiagonal.jl:71,101) — too shallow to trigger it.
#
# Two things follow, and this harness exists for both:
#
# 1. The existing device coverage (scripts/gpu_validate.jl:81, test_tridiagonal_multi)
#    runs nlevs = 20 but pins `jtop = ones(Int, ncols)`. The jtop > 1 sub-path —
#    the one #263 names as corrupted — has NO direct device coverage. It is a LIVE
#    production configuration: lake_temperature.jl:1219 passes
#    `jtop_ext = jtop .+ nlevsno`, i.e. strictly > 1. Its only device exercise
#    today is the whole-lake composite gpu_validate_lake_e2e.jl, where a
#    solver-local error is diluted across lt1..lt7 and is hard to attribute.
#
# 2. The reference is a DENSE per-column `LinearAlgebra.Tridiagonal` solve, NOT the
#    KA CPU result. This is the whole methodological point: #263's bug made the
#    CPU wrong and the device right, so a device-vs-host comparison would have
#    reported "parity" while both were measured against a broken yardstick. An
#    independent oracle catches a regression on EITHER side.
#
#   julia --project=scripts scripts/gpu_validate_tridiag_jtop.jl
#
# Metal has no Float64, so the device runs Float32 and is compared to the Float64
# oracle at the Float32 solver floor, matching the other gpu_validate_* harnesses.
# ==========================================================================

using CLM
using Printf
using Random
using LinearAlgebra
include(joinpath(@__DIR__, "gpu_backends.jl"))

# Sentinel prefill for `u`: every entry the kernel must NOT touch (masked-out
# columns, and rows above jtop) has to still hold this afterwards. A plain zero
# fill would make "untouched" indistinguishable from "written zero".
const SENTINEL = -987654.0

"""
Build a batched, diagonally dominant tridiagonal system with a deliberately
RAGGED jtop and a partly-false mask.

`jtop` cycles through 1, 4, and nlevsno+1 so the harness spans the CTSM-realistic
range (soil-only start, mid-column start, and the lake/snow offset start that
lake_temperature.jl actually produces).
"""
function build(::Type{FT}; ncols=48, nlevs=21, nlevsno=5, seed=2) where {FT}
    rng = MersenneTwister(seed)
    a = FT.(-0.5 .* ones(ncols, nlevs))
    c = FT.(-0.5 .* ones(ncols, nlevs))
    b = FT.(2 .+ 0.5 .* rand(rng, ncols, nlevs))
    r = FT.(rand(rng, ncols, nlevs))

    jtop = [ (col % 3 == 0) ? 1 : (col % 3 == 1) ? 4 : nlevsno + 1 for col in 1:ncols ]
    # Every column's first active row has no sub-diagonal, its last no super-diagonal.
    for col in 1:ncols
        a[col, jtop[col]] = 0
        c[col, nlevs] = 0
    end

    mask = trues(ncols)
    mask[2] = false          # a masked column at a jtop > 1 offset
    mask[ncols] = false
    return (; a, b, c, r, jtop, mask, ncols, nlevs)
end

"""Dense per-column oracle: solve rows jtop..nlevs with LinearAlgebra."""
function dense_reference(d)
    u = fill(SENTINEL, d.ncols, d.nlevs)
    for col in 1:d.ncols
        d.mask[col] || continue
        j = d.jtop[col]
        n = d.nlevs - j + 1
        n <= 0 && continue
        dl = Float64[d.a[col, k] for k in (j+1):d.nlevs]   # sub-diagonal
        dg = Float64[d.b[col, k] for k in j:d.nlevs]       # diagonal
        du = Float64[d.c[col, k] for k in j:(d.nlevs-1)]   # super-diagonal
        rhs = Float64[d.r[col, k] for k in j:d.nlevs]
        u[col, j:d.nlevs] = Tridiagonal(dl, dg, du) \ rhs
    end
    return u
end

"""Max relative deviation over only the ACTIVE (mask && row>=jtop) entries."""
function active_reldiff(got, ref, d)
    G = Array(got); R = Array(ref); m = 0.0
    for col in 1:d.ncols
        d.mask[col] || continue
        for k in d.jtop[col]:d.nlevs
            g = Float64(G[col, k]); rr = Float64(R[col, k])
            m = max(m, abs(g - rr) / (1.0 + abs(rr)))
        end
    end
    return m
end

"""Verify the kernel left every entry it must not touch at SENTINEL."""
function untouched_ok(got, d)
    G = Array(got); bad = 0
    for col in 1:d.ncols
        if !d.mask[col]
            for k in 1:d.nlevs
                Float64(G[col, k]) == Float32(SENTINEL) || Float64(G[col, k]) == SENTINEL || (bad += 1)
            end
        else
            for k in 1:(d.jtop[col]-1)
                Float64(G[col, k]) == Float32(SENTINEL) || Float64(G[col, k]) == SENTINEL || (bad += 1)
            end
        end
    end
    return bad
end

function run_case(backend, nlevs)
    name, to, FT = backend
    nfail = 0
    @printf("\n--- nlevs = %d ---------------------------------------------\n", nlevs)

    d64 = build(Float64; nlevs=nlevs)
    ref = dense_reference(d64)

    # ---- anti-vacuity guards: this harness must fail loudly if its own fixture
    # ---- ever stops exercising the thing it was written to exercise.
    nd = length(unique(d64.jtop))
    if nd < 3
        println("  [FAIL] fixture is not ragged — jtop has $nd distinct values (need >= 3)"); nfail += 1
    else
        println("  [PASS] jtop is ragged: distinct values = ", sort(unique(d64.jtop)))
    end
    if !any(d64.jtop .> 1)
        println("  [FAIL] fixture has no jtop > 1 column — the whole point of this harness"); nfail += 1
    else
        @printf("  [PASS] %d of %d columns have jtop > 1\n", count(d64.jtop .> 1), d64.ncols)
    end
    if all(d64.mask)
        println("  [FAIL] fixture has no masked-out column"); nfail += 1
    else
        println("  [PASS] fixture masks out $(count(.!d64.mask)) column(s)")
    end
    # The oracle must be non-trivial and column-dependent, else any solver "passes".
    act = [ref[col, k] for col in 1:d64.ncols if d64.mask[col] for k in d64.jtop[col]:nlevs]
    if !(maximum(abs, act) > 1e-3) || length(unique(round.(act; digits=9))) < 10
        println("  [FAIL] oracle solution is trivial/degenerate — test would be vacuous"); nfail += 1
    else
        @printf("  [PASS] oracle non-trivial: max|u| = %.4f over %d active entries\n",
                maximum(abs, act), length(act))
    end

    # ---- HOST (KA CPU backend) vs dense oracle. This is the leg that would have
    # ---- caught #263 on its own.
    u_host = fill(SENTINEL, d64.ncols, nlevs)
    CLM.tridiagonal_multi!(u_host, d64.a, d64.b, d64.c, d64.r, d64.jtop, d64.mask,
                           d64.ncols, nlevs)
    dh = active_reldiff(u_host, ref, d64)
    okh = dh < 1e-10
    @printf("  [%s] host  vs dense oracle   rel = %.3e  (tol 1e-10)\n", okh ? "PASS" : "FAIL", dh)
    okh || (nfail += 1)
    bh = untouched_ok(u_host, d64)
    @printf("  [%s] host  left %d must-not-touch entries modified\n", bh == 0 ? "PASS" : "FAIL", bh)
    bh == 0 || (nfail += 1)

    if backend === nothing
        return nfail
    end

    # ---- DEVICE vs dense oracle, at the backend's working precision.
    dF = build(FT; nlevs=nlevs)
    u_dev = to(fill(FT(SENTINEL), dF.ncols, nlevs))
    CLM.tridiagonal_multi!(u_dev, _dev(to, dF.a), _dev(to, dF.b), _dev(to, dF.c),
                           _dev(to, dF.r), _dev(to, dF.jtop), _dev(to, dF.mask),
                           dF.ncols, nlevs)
    device_synchronize()

    tol = FT === Float32 ? 5e-3 : 1e-10
    dd = active_reldiff(u_dev, ref, dF)
    okd = dd < tol
    @printf("  [%s] %-5s vs dense oracle   rel = %.3e  (tol %.0e)\n",
            okd ? "PASS" : "FAIL", name, dd, tol)
    okd || (nfail += 1)
    bd = untouched_ok(u_dev, dF)
    @printf("  [%s] %-5s left %d must-not-touch entries modified\n",
            bd == 0 ? "PASS" : "FAIL", name, bd)
    bd == 0 || (nfail += 1)

    # ---- device vs host, for continuity with the other harnesses.
    dhd = active_reldiff(u_dev, u_host, dF)
    okhd = dhd < tol
    @printf("  [%s] %-5s vs host           rel = %.3e  (tol %.0e)\n",
            okhd ? "PASS" : "FAIL", name, dhd, tol)
    okhd || (nfail += 1)

    return nfail
end

function main(backend)
    println("=" ^ 72)
    println("tridiagonal_multi! — RAGGED jtop>1 device parity vs dense oracle")
    println("=" ^ 72)
    if backend === nothing
        println("  No GPU backend detected — running the host-vs-oracle legs only.")
    else
        name, _, FT = backend
        @printf("  Backend: %s   (working precision: %s)\n", name, FT)
    end

    # 21 = production soil depth (nlevsoi=20 -> nlevsoi+1, soil_water_movement.jl:1717)
    # 25 = the padded snow+soil depth used by the lake/band configurations.
    nfail = run_case(backend, 21) + run_case(backend, 25)

    println()
    println(nfail == 0 ? "  tridiagonal_multi! ragged-jtop parity: PASS" :
                         "  DIVERGENCE — investigate ($nfail failed).")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
