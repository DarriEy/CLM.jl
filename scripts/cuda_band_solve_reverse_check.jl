# Is Enzyme's reverse of the batched band-solve kernel itself right on CUDA?
#
# Differentiates the raw kernel launch (bypassing the batched_band_solve! rule) on
# CPU and CUDA and compares the adjoints of rvector and bmatrix with the exact
# linear-solve adjoint (_batched_band_solve_pullback!, FD-checked in
# test/test_band_solve_adjoint.jl).
#   julia --project=scripts scripts/cuda_band_solve_reverse_check.jl
using CLM, Enzyme, CUDA, Printf, Random

function systems()
    rng = MersenneTwister(7)
    nc, ns, nt, kl = 3, 2, 9, 2
    jtop = [-1, 0, -2]; jbot = [5, 6, 4]; mask = [true, true, true]
    bm = zeros(nc, 2kl + 1, nt); rv = zeros(nc, nt)
    for c in 1:nc
        jt = jtop[c] + ns + 1; n = jbot[c] - jtop[c] + 1
        for i in 1:n, j in max(1, i - kl):min(n, i + kl)
            bm[c, (i - j) + kl + 1, jt + i - 1] = i == j ? 6 + rand(rng) : rand(rng) - 0.5
        end
        for i in 1:n
            rv[c, jt + i - 1] = randn(rng)
        end
    end
    w = zeros(nc, nt)
    for c in 1:nc, i in 1:(jbot[c] - jtop[c] + 1)
        w[c, jtop[c] + ns + i] = randn(rng)
    end
    (; nc, ns, nt, kl, jtop, jbot, mask, bm, rv, w)
end

# The raw launch, as batched_band_solve! does it, but with its scratch passed in.
raw_solve!(t, bm, rv, jtop, jbot, mask, A, x, kl, ns) =
    (CLM._launch!(CLM._batched_band_solve_kernel!, t, bm, rv, jtop, jbot, mask, A, x,
                  kl, kl, ns; ndrange = size(t, 1)); nothing)

function check(AT, S)
    t = AT(zeros(S.nc, S.nt)); dt = AT(copy(S.w))
    bm = AT(copy(S.bm)); dbm = AT(zeros(size(S.bm)))
    rv = AT(copy(S.rv)); drv = AT(zeros(size(S.rv)))
    A = AT(zeros(S.nc, S.nt, S.nt)); dA = AT(zeros(S.nc, S.nt, S.nt))
    x = AT(zeros(S.nc, S.nt)); dx = AT(zeros(S.nc, S.nt))
    jtop = AT(S.jtop); jbot = AT(S.jbot); mask = AT(S.mask)
    Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.Reverse), raw_solve!, Enzyme.Const,
        Enzyme.Duplicated(t, dt), Enzyme.Duplicated(bm, dbm), Enzyme.Duplicated(rv, drv),
        Enzyme.Const(jtop), Enzyme.Const(jbot), Enzyme.Const(mask),
        Enzyme.Duplicated(A, dA), Enzyme.Duplicated(x, dx),
        Enzyme.Const(S.kl), Enzyme.Const(S.ns))
    # exact adjoint on the host
    th = zeros(S.nc, S.nt)
    CLM.batched_band_solve!(th, S.bm, S.rv, S.jtop, S.jbot, S.mask, S.kl, S.kl, S.ns)
    edt = copy(S.w); edrv = zeros(size(S.rv)); edbm = zeros(size(S.bm))
    CLM._batched_band_solve_pullback!(edt, edrv, edbm, th, S.bm, S.jtop, S.jbot, S.mask,
                                      S.kl, S.kl, S.ns)
    er = maximum(abs, Array(drv) .- edrv) / maximum(abs, edrv)
    eb = maximum(abs, Array(dbm) .- edbm) / maximum(abs, edbm)
    @printf("%-5s rel err  d(rvector)=%.3e  d(bmatrix)=%.3e\n", AT === Array ? "CPU" : "CUDA", er, eb)
    max(er, eb)
end

println("Julia ", VERSION, ", Enzyme ", pkgversion(Enzyme), ", CUDA.jl ", pkgversion(CUDA))
S = systems()
ecpu = check(Array, S)
egpu = check(CuArray, S)
ok = ecpu < 1e-10 && egpu < 1e-10
println(ok ? "both match the exact adjoint" : "MISMATCH")
exit(ok ? 0 : 1)
