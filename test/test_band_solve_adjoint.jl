# The batched_band_solve! adjoint used by its GPU reverse rule
# (_batched_band_solve_pullback!, src/infrastructure/band_diagonal.jl), checked on
# host arrays against central finite differences of L = w ⋅ x(bmatrix, rvector).
using Test, Random
import CLM

@testset "batched_band_solve! adjoint vs FD" begin
    rng = MersenneTwister(7)
    nc, ns, nt, kl = 3, 2, 9, 2
    jtop = [-1, 0, -2]            # snl-like tops; system spans jtop..jbot
    jbot = [5, 6, 4]
    mask = [true, true, true]
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
    w = randn(rng, nc, nt)
    function L(bm, rv)
        t = zeros(nc, nt)
        CLM.batched_band_solve!(t, bm, rv, jtop, jbot, mask, kl, kl, ns)
        s = 0.0
        for c in 1:nc
            jt = jtop[c] + ns + 1
            for i in 1:(jbot[c] - jtop[c] + 1)
                s += w[c, jt + i - 1] * t[c, jt + i - 1]
            end
        end
        s
    end
    x = zeros(nc, nt)
    CLM.batched_band_solve!(x, bm, rv, jtop, jbot, mask, kl, kl, ns)
    dtv = copy(w); drv = zeros(nc, nt); dbm = zeros(size(bm))
    CLM._batched_band_solve_pullback!(dtv, drv, dbm, x, bm, jtop, jbot, mask, kl, kl, ns)
    h = 1e-6; errs = Float64[]
    for c in 1:nc
        jt = jtop[c] + ns + 1; n = jbot[c] - jtop[c] + 1
        for i in 1:n
            r1 = copy(rv); r1[c, jt + i - 1] += h; r2 = copy(rv); r2[c, jt + i - 1] -= h
            push!(errs, abs((L(bm, r1) - L(bm, r2)) / 2h - drv[c, jt + i - 1]))
            @test dtv[c, jt + i - 1] == 0
            for j in max(1, i - kl):min(n, i + kl)
                k = (i - j) + kl + 1
                b1 = copy(bm); b1[c, k, jt + i - 1] += h; b2 = copy(bm); b2[c, k, jt + i - 1] -= h
                push!(errs, abs((L(b1, rv) - L(b2, rv)) / 2h - dbm[c, k, jt + i - 1]))
            end
        end
    end
    @test maximum(errs) < 1e-7
end
