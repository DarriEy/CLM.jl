# ==========================================================================
# gpu_ad_validate.jl — real-hardware AUTODIFF validation for CLM.jl's kernels.
#
# Phase 5: verify that automatic differentiation works THROUGH the Phase-4 kernels
# when they run on the GPU. Forward-mode (ForwardDiff `Dual` numbers) flows through
# a kernel as a custom isbits element type: an array of `Dual`s on the device runs
# the same kernel and the derivative rides along in the partial component.
#
# For each kernel we differentiate the output w.r.t. a scalar θ added to one input
# (so every output carries a partial), then check:
#   1. PARITY  — device-AD partials vs CPU-AD partials (the core question: does AD
#                give the SAME answer on the GPU as on the CPU?). PASS/FAIL on this.
#   2. CORRECT — AD partials vs a finite-difference estimate (reported; confirms the
#                derivative itself is right, not just that CPU and GPU agree).
#
# RUN
#   julia --project=scripts scripts/gpu_ad_validate.jl
#
# REVERSE MODE (Enzyme): forward-mode is what is validated here. Reverse-mode
# through these kernels works on the CPU backend, but Enzyme does not differentiate
# the Metal device kernel launch (no Enzyme GPU support for the Metal/Apple backend);
# see the note printed at the end. On CUDA, Enzyme-over-KA is the intended reverse path.
# ==========================================================================

using CLM
using Printf
using Random
using ForwardDiff

include(joinpath(@__DIR__, "gpu_backends.jl"))

# --- forward-mode helpers ---------------------------------------------------
_dualT(::Type{FT}) where {FT} = ForwardDiff.Dual{Nothing,FT,1}
_seed(v, p) = ForwardDiff.Dual(v, p)                 # value v, single partial p
_todual(::Type{FT}, A) where {FT} = _dualT(FT).(A)   # lift to Dual, partial 0
_part(A) = ForwardDiff.partials.(Array(A), 1)        # extract the 1st partial
_val(A)  = ForwardDiff.value.(Array(A))

maxabs(x) = maximum(abs.(x); init = 0.0)
relerr(a, b) = maxabs(a .- b) / max(maxabs(b), eps(Float64))

# PASS thresholds. Parity (CPU-AD vs device-AD) must be tight; the finite-difference
# correctness check is generous (Float32 FD is noisy, esp. through the solver sweeps).
_parity_tol(::Type{Float32}) = 1f-3
_parity_tol(::Type{Float64}) = 1e-9
const FD_REL_TOL = 0.05    # 5% relative — finite-difference sanity, not a tight bound

# Each test returns (name, parity, fd_rel, ok).
function adtest_forc_q(to, ::Type{FT}) where {FT}
    rng = MersenneTwister(1); nc = 32
    gc = collect(1:nc)
    vp  = FT.(500 .+ 600 .* rand(rng, nc))
    pb  = FT.(80000 .+ 10000 .* rand(rng, nc))
    # θ added to every vp → each out carries d(out)/d(vp)
    vpd = _todual(FT, vp); for i in eachindex(vpd); vpd[i] = _seed(vp[i], one(FT)); end
    pbd = _todual(FT, pb)
    oc = zeros(_dualT(FT), nc);    CLM.compute_forc_q!(oc, gc, vpd, pbd)
    od = to(zeros(_dualT(FT), nc)); CLM.compute_forc_q!(od, _dev(to,gc), _dev(to,vpd), _dev(to,pbd))
    parity = maxabs(_part(oc) .- _part(od))
    # finite difference wrt θ: perturb all vp
    h = FT(1f-1)
    o0 = zeros(FT, nc); CLM.compute_forc_q!(o0, gc, vp, pb)
    o1 = zeros(FT, nc); CLM.compute_forc_q!(o1, gc, vp .+ h, pb)
    fd = (o1 .- o0) ./ h
    fdrel = relerr(_part(od), fd)
    return ("compute_forc_q!", parity, fdrel, parity < _parity_tol(FT) && fdrel < FD_REL_TOL)
end

function adtest_tridiagonal(to, ::Type{FT}) where {FT}
    rng = MersenneTwister(2); nc, nl = 24, 16
    a = FT.(-0.5 .* ones(nc, nl)); c = FT.(-0.5 .* ones(nc, nl))
    b = FT.(2 .+ 0.5 .* rand(rng, nc, nl)); r = FT.(rand(rng, nc, nl))
    a[:,1] .= 0; c[:,nl] .= 0; jt = ones(Int, nc); mk = trues(nc)
    rd = _todual(FT, r); for i in eachindex(rd); rd[i] = _seed(r[i], one(FT)); end
    ad = _todual(FT, a); bd = _todual(FT, b); cd = _todual(FT, c)
    uc = zeros(_dualT(FT), nc, nl); CLM.tridiagonal_multi!(uc, ad, bd, cd, rd, jt, mk, nc, nl)
    ud = to(zeros(_dualT(FT), nc, nl))
    CLM.tridiagonal_multi!(ud, _dev(to,ad), _dev(to,bd), _dev(to,cd), _dev(to,rd), _dev(to,jt), _dev(to,mk), nc, nl)
    parity = maxabs(_part(uc) .- _part(ud))
    h = FT(1f-3)
    u0 = zeros(FT, nc, nl); CLM.tridiagonal_multi!(u0, a, b, c, r,        jt, mk, nc, nl)
    u1 = zeros(FT, nc, nl); CLM.tridiagonal_multi!(u1, a, b, c, r .+ h,   jt, mk, nc, nl)
    fd = (u1 .- u0) ./ h
    fdrel = relerr(_part(ud), fd)
    return ("tridiagonal_multi! (batched Thomas)", parity, fdrel, parity < _parity_tol(FT) && fdrel < FD_REL_TOL)
end

function adtest_band(to, ::Type{FT}) where {FT}
    rng = MersenneTwister(3); nc = 24; kl = ku = 2; nband = 5; nlevsno = 5; nlev = 25
    jt = ones(Int, nc); jb = fill(15, nc); mk = trues(nc)
    bm = zeros(FT, nc, nband, nlev)
    for c in 1:nc, lev in 1:nlev
        bm[c, kl+1, lev] = FT(4) + rand(rng, FT)
        bm[c, kl, lev] = FT(-0.7); bm[c, kl+2, lev] = FT(-0.7)
        bm[c, 1, lev] = FT(-0.3);  bm[c, nband, lev] = FT(-0.3)
    end
    rv = FT.(rand(rng, nc, nlev))
    rvd = _todual(FT, rv); for i in eachindex(rvd); rvd[i] = _seed(rv[i], one(FT)); end
    bmd = _todual(FT, bm)
    tc = zeros(_dualT(FT), nc, nlev); CLM.batched_band_solve!(tc, bmd, rvd, jt, jb, mk, kl, ku, nlevsno)
    td = to(zeros(_dualT(FT), nc, nlev))
    CLM.batched_band_solve!(td, _dev(to,bmd), _dev(to,rvd), _dev(to,jt), _dev(to,jb), _dev(to,mk), kl, ku, nlevsno)
    parity = maxabs(_part(tc) .- _part(td))
    h = FT(1f-3)
    t0 = zeros(FT, nc, nlev); CLM.batched_band_solve!(t0, bm, rv,      jt, jb, mk, kl, ku, nlevsno)
    t1 = zeros(FT, nc, nlev); CLM.batched_band_solve!(t1, bm, rv .+ h, jt, jb, mk, kl, ku, nlevsno)
    fd = (t1 .- t0) ./ h
    fdrel = relerr(_part(td), fd)
    return ("batched_band_solve! (pentadiagonal GE)", parity, fdrel, parity < _parity_tol(FT) && fdrel < FD_REL_TOL)
end

const AD_TESTS = [adtest_forc_q, adtest_tridiagonal, adtest_band]

function main(backend)
    println("=" ^ 74)
    println("GPU AUTODIFF VALIDATION for CLM.jl kernels (forward-mode, CPU vs device)")
    println("=" ^ 74)
    if backend === nothing
        println("""
  No GPU backend detected (CUDA / AMDGPU / oneAPI / Metal).
  Forward-mode AD through the kernels is exercised on the CPU by the test suite.
  Install a backend package and re-run on real GPU hardware.
""")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n", name, FT)
    println("  PARITY = max|device-AD partial − CPU-AD partial|   (must be ~0)")
    println("  FD-rel = relative error of AD partial vs finite difference (sanity)\n")

    npass = 0; nfail = 0
    for t in AD_TESTS
        nm, parity, fdrel, ok = try
            t(to, FT)
        catch e
            (string(t) * "  (ERROR: " * first(split(sprint(showerror, e), Char(10))) * ")", NaN, NaN, false)
        end
        @printf("  [%s] %-42s  PARITY=%.2e  FD-rel=%.2e\n", ok ? "PASS" : "FAIL", nm, parity, fdrel)
        ok ? (npass += 1) : (nfail += 1)
    end
    @printf("\n  %d passed, %d failed\n", npass, nfail)
    if nfail == 0
        println("\n  FORWARD-MODE AD MATCHES CPU ON $(name) ($(FT)) ✓")
        println("  Reverse-mode (Enzyme) through these kernels works on CPU; Enzyme does not")
        println("  differentiate the Metal device launch (no Enzyme GPU support for Apple).")
    else
        println("\n  SOME AD CHECKS FAILED — investigate before trusting GPU gradients.")
    end
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()   # top-level: advance world age before main()
exit(main(BACKEND))
