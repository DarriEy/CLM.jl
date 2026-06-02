# ==========================================================================
# gpu_validate_daylength_e2e.jl — end-to-end Metal parity for the WHOLE
# update_daylength! driver (per-grid-cell advance of prev_dayl/dayl +
# compute_max_daylength!), plus the underlying elemental daylength().
#
# Builds a multi-grid-cell GridcellData at Float32 with latitudes spanning the
# pole, mid-NH, equator, mid-SH (so the >0 / <=0 hemisphere branch in
# compute_max_daylength! and the pole-clamp branch in daylength() are both
# exercised), runs update_daylength! on the CPU, adapts grc to Metal, runs the
# SAME call on the device, and compares dayl / prev_dayl / max_dayl with reldiff.
#
# CRITICAL: latitudes are real (radians) so the CPU reference is FINITE on the
# non-polar cells; we assert finiteness before trusting parity.
#
#   julia --project=scripts scripts/gpu_validate_daylength_e2e.jl
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

function build(::Type{FT}) where {FT}
    ng = 5
    grc = CLM.GridcellData{FT}()
    CLM.gridcell_init!(grc, ng)
    # latitudes (radians): near N pole, mid-NH, equator, mid-SH, near S pole
    lats = FT[1.5, 0.7, 0.0, -0.7, -1.5]
    for g in 1:ng
        grc.lat[g]       = lats[g]
        grc.dayl[g]      = FT(43200.0)   # 12 h, a finite "previous step" value
        grc.prev_dayl[g] = FT(40000.0)
        grc.max_dayl[g]  = FT(NaN)
    end
    declin    = 0.2     # solar declination (radians), ~ mid-spring
    obliquity = 0.409   # Earth's obliquity (radians)
    return (; ng, grc, declin, obliquity)
end

run_ud!(H) = CLM.update_daylength!(H.grc, H.declin, H.obliquity, false, 1:H.ng)

function main(backend)
    println("=" ^ 70)
    println("END-TO-END Metal parity for update_daylength! (whole driver)")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU driver exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    H = build(FT)
    B = build(FT)

    grc_d = CLM.Adapt.adapt(Metal.MtlArray, B.grc)
    if !(grc_d.dayl isa Metal.MtlArray)
        println("  BLOCKED: GridcellData did not move to the device under adapt.")
        return 2
    end

    run_ud!(H)
    CLM.update_daylength!(grc_d, B.declin, B.obliquity, false, 1:B.ng)

    # Guard: the non-polar cells (2,3,4) must be finite in the CPU reference.
    for g in (2, 3, 4)
        if !(isfinite(H.grc.dayl[g]) && isfinite(H.grc.max_dayl[g]))
            @printf("  BLOCKED: CPU dayl/max_dayl[%d] not finite.\n", g)
            return 2
        end
    end

    checks = [
        ("dayl",      H.grc.dayl,      grc_d.dayl),
        ("prev_dayl", H.grc.prev_dayl, grc_d.prev_dayl),
        ("max_dayl",  H.grc.max_dayl,  grc_d.max_dayl),
    ]
    nfail = 0; worst = 0.0
    for (nm, a, b) in checks
        d = reldiff(a, b); ok = d < 1f-2
        worst = max(worst, d)
        @printf("  [%s] %-15s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
        ok || (nfail += 1)
    end
    @printf("\n  worst reldiff = %.3e\n", worst)
    println(nfail == 0 ? "  WHOLE update_daylength! MATCHES CPU ON $name ($FT)" :
                         "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
