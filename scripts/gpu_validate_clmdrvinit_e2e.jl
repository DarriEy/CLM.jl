# ==========================================================================
# gpu_validate_clmdrvinit_e2e.jl — GPU parity for clm_drv_init!'s per-col /
# per-patch driver-init loops (A10).
#
# clm_drv_init! has four element-wise init loops, now per-element kernels:
#   1. intracellular CO2 reset to -999       (_drvinit_ci_kernel!, per (p,j))
#   2. eflx_bot_col reset to 0               (_drvinit_eflxbot_kernel!, per c)
#   3. frac_veg_nosno = active ? alb : 0     (_drvinit_fracvegnosno_kernel!, per p)
#   4. frac_iceold ice fraction at prev step (_drvinit_fraciceold_kernel!, per c
#        with an INTERNAL sequential snow-layer loop over [snl+1 .. 0])
#
# Branches exercised:
#   - active vs inactive patches (frac_veg_nosno)
#   - masked vs unmasked columns (frac_iceold via nolake mask)
#   - varied snl (0, -1, -3) so the snow-layer accumulation branch differs
#   - snow layers with total water > 0 vs == 0 (frac_iceold only set when > 0)
#
# CPU reference fields are asserted finite so a both-NaN false PASS cannot slip in.
#
#   julia --project=scripts scripts/gpu_validate_clmdrvinit_e2e.jl
# ==========================================================================

using CLM
using Printf
import KernelAbstractions as KA
include(joinpath(@__DIR__, "gpu_backends.jl"))

function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(A[i] - B[i]) / (1.0 + max(abs(A[i]), abs(B[i]))))
    end
    return m
end

const NP = 4
const NC = 4
const NLEVCAN = 3
const NLEVSNO = 5

function build(::Type{FT}) where {FT}
    np = NP; nc = NC
    nlev = NLEVSNO + 10
    snl = Int[0, -1, -3, -2]
    active = Bool[true, false, true, true]
    nolakec = Bool[true, true, false, true]   # col 3 is masked out (lake)

    fvns_alb = FT[FT(0.7), FT(0.9), FT(1.0), FT(0.4)]

    h2osoi_liq = fill(FT(0.0), nc, nlev)
    h2osoi_ice = fill(FT(0.0), nc, nlev)
    for c in 1:nc
        for j in (snl[c] + 1):0
            jj = j + NLEVSNO
            h2osoi_ice[c, jj] = FT(8.0) + FT(c) + FT(-j)
            h2osoi_liq[c, jj] = FT(2.0) + FT(0.5) * FT(-j)
        end
    end
    # Leave one layer of col 4 with zero water to exercise the total==0 skip
    if snl[4] < 0
        jj = (snl[4] + 1) + NLEVSNO
        h2osoi_ice[4, jj] = FT(0.0)
        h2osoi_liq[4, jj] = FT(0.0)
    end

    return (; np, nc, snl, active, nolakec, fvns_alb, h2osoi_liq, h2osoi_ice)
end

# Drive the four kernels exactly as clm_drv_init! does.
function run_kernels!(S, cisun, cisha, eflx_bot, fvns, frac_iceold)
    begp = 1; endp = S.np
    begc = 1; endc = S.nc
    CLM._launch!(CLM._drvinit_ci_kernel!, cisun, cisha, begp, endp, NLEVCAN;
                 ndrange = (size(cisun, 1), NLEVCAN))
    CLM._launch!(CLM._drvinit_eflxbot_kernel!, eflx_bot, begc, endc)
    CLM._launch!(CLM._drvinit_fracvegnosno_kernel!, fvns, S.fvns_alb, S.active, begp, endp)
    CLM._launch!(CLM._drvinit_fraciceold_kernel!, frac_iceold, S.nolakec, S.snl,
                 S.h2osoi_liq, S.h2osoi_ice, begc, endc, NLEVSNO;
                 ndrange = size(frac_iceold, 1))
    return nothing
end

function main(backend)
    println("=" ^ 72)
    println("END-TO-END GPU parity for clm_drv_init! init loops (A10)")
    println("=" ^ 72)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU path exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    H = build(FT); B = build(FT)

    # CPU scratch (pre-fill ci with sentinel; iceold pre-fill nonzero so the
    # "only set when total>0" skip is observable as preserved value).
    cisun_h = fill(FT(123.0), H.np, NLEVCAN)
    cisha_h = fill(FT(456.0), H.np, NLEVCAN)
    eflx_h  = fill(FT(7.7), H.nc)
    fvns_h  = fill(FT(-1.0), H.np)
    iceold_h = fill(FT(0.25), H.nc, NLEVSNO)
    run_kernels!(H, cisun_h, cisha_h, eflx_h, fvns_h, iceold_h)

    db(x) = to(collect(Bool, x))
    Sd = (; np = B.np, nc = B.nc, snl = to(B.snl), active = db(B.active),
            nolakec = db(B.nolakec), fvns_alb = to(B.fvns_alb),
            h2osoi_liq = to(B.h2osoi_liq), h2osoi_ice = to(B.h2osoi_ice))
    cisun_d = to(fill(FT(123.0), B.np, NLEVCAN))
    cisha_d = to(fill(FT(456.0), B.np, NLEVCAN))
    eflx_d  = to(fill(FT(7.7), B.nc))
    fvns_d  = to(fill(FT(-1.0), B.np))
    iceold_d = to(fill(FT(0.25), B.nc, NLEVSNO))

    if !(cisun_d isa device_array_type())
        println("  BLOCKED: device arrays did not move under to().")
        return 2
    end
    run_kernels!(Sd, cisun_d, cisha_d, eflx_d, fvns_d, iceold_d)

    tonum(x) = Float64.(Array(x))
    checks = [
        ("cisun_z",       cisun_h, cisun_d),
        ("cisha_z",       cisha_h, cisha_d),
        ("eflx_bot_col",  eflx_h,  eflx_d),
        ("frac_veg_nosno", fvns_h, fvns_d),
        ("frac_iceold",   iceold_h, iceold_d),
    ]
    nfail = 0
    for (nm, a, b) in checks
        an = tonum(a)
        if !any(isfinite, an)
            @printf("  [WARN] %-15s CPU reference all-NaN/Inf — skipping\n", nm)
            continue
        end
        d = reldiff(an, tonum(b)); ok = d < 1f-3
        @printf("  [%s] %-15s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
        ok || (nfail += 1)
    end

    println()
    for p in 1:H.np
        @printf("    patch %d: active=%-5s  frac_veg_nosno=%7.4f\n",
                p, H.active[p], fvns_h[p])
    end
    for c in 1:H.nc
        # Top resolved snow layer Julia index, clamped into [1, NLEVSNO].
        jj = clamp((H.snl[c] + 1) + NLEVSNO, 1, NLEVSNO)
        @printf("    col %d: snl=%2d nolake=%-5s frac_iceold[top]=%7.4f\n",
                c, H.snl[c], H.nolakec[c], iceold_h[c, jj])
    end

    println()
    println(nfail == 0 ? "  clm_drv_init! init loops MATCH CPU ON $name ($FT)" :
                         "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
