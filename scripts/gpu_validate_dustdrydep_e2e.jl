# ==========================================================================
# gpu_validate_dustdrydep_e2e.jl — end-to-end Metal parity for the WHOLE
# dust_dry_dep! turbulent dry-deposition driver (3 per-patch passes:
# thermokinetic props + Stokes settling, quasi-laminar resistance, lowest-layer
# turbulent deposition + per-bin copy).
#
# Builds a multi-patch / multi-column DustEmisBaseData at Float32, runs
# dust_dry_dep! on the CPU, adapts the dust struct + every forcing array the
# kernels touch to Metal, runs the SAME call on the device, and compares the
# mutated outputs with reldiff.
#
# Exercises active AND inactive patches (the active mask), several columns at
# different pbot/rho/t, and all 4 dust bins.
#
# CRITICAL: every input is a real CLM default so the CPU reference is FINITE —
# reldiff silently PASSES when both sides are NaN, so we assert the CPU
# reference fields are finite before trusting parity.
#
#   julia --project=scripts scripts/gpu_validate_dustdrydep_e2e.jl
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
    np = 4; nc = 3
    dust = CLM.DustEmisBaseData{FT}()
    CLM.dust_emis_init!(dust, np)

    patch_active = trues(np)
    patch_active[3] = false                 # one inactive patch (mask branch)
    patch_column = [1, 2, 2, 3]
    bounds_p = 1:np

    forc_pbot = FT[101325.0, 85000.0, 101325.0]   # sea level + ~1500 m
    forc_rho  = FT[1.225, 1.05, 1.225]
    forc_t    = FT[300.0, 280.0, 293.0]
    ram1 = fill(FT(40.0), np)
    fv   = FT[0.25, 0.30, 0.20, 0.35]

    return (; np, nc, dust, patch_active, patch_column, bounds_p,
              forc_pbot, forc_rho, forc_t, ram1, fv)
end

run_ddd!(H) = CLM.dust_dry_dep!(H.dust, H.patch_active, H.patch_column, H.bounds_p,
                                H.forc_pbot, H.forc_rho, H.forc_t, H.ram1, H.fv)

function main(backend)
    println("=" ^ 70)
    println("END-TO-END Metal parity for dust_dry_dep! (whole driver)")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU driver exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    H = build(FT)   # CPU reference
    B = build(FT)   # device copy

    ad(x) = CLM.Adapt.adapt(Metal.MtlArray, x)
    dust_d = ad(B.dust)
    forc_pbot_d = ad(B.forc_pbot); forc_rho_d = ad(B.forc_rho); forc_t_d = ad(B.forc_t)
    ram1_d = ad(B.ram1); fv_d = ad(B.fv)
    # Device arrays disallow BitArray; move masks/index vectors as plain arrays.
    active_d = to(collect(Bool, B.patch_active))
    col_d    = to(collect(Int32, B.patch_column))

    if !(dust_d.vlc_trb_patch isa Metal.MtlArray)
        println("  BLOCKED: dust struct did not move to the device under adapt.")
        return 2
    end

    run_ddd!(H)
    CLM.dust_dry_dep!(dust_d, active_d, col_d, B.bounds_p,
                      forc_pbot_d, forc_rho_d, forc_t_d, ram1_d, fv_d)

    # Guard against a false PASS: the CPU reference must be FINITE on active patches.
    for p in (1, 2, 4), m in 1:CLM.NDST
        if !isfinite(H.dust.vlc_trb_patch[p, m])
            @printf("  BLOCKED: CPU vlc_trb_patch[%d,%d] not finite.\n", p, m)
            return 2
        end
    end

    checks = [
        ("vlc_trb_patch", H.dust.vlc_trb_patch,   dust_d.vlc_trb_patch),
        ("vlc_trb_1",     H.dust.vlc_trb_1_patch, dust_d.vlc_trb_1_patch),
        ("vlc_trb_2",     H.dust.vlc_trb_2_patch, dust_d.vlc_trb_2_patch),
        ("vlc_trb_3",     H.dust.vlc_trb_3_patch, dust_d.vlc_trb_3_patch),
        ("vlc_trb_4",     H.dust.vlc_trb_4_patch, dust_d.vlc_trb_4_patch),
    ]
    nfail = 0; worst = 0.0
    for (nm, a, b) in checks
        d = reldiff(a, b); ok = d < 1f-2
        worst = max(worst, d)
        @printf("  [%s] %-15s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
        ok || (nfail += 1)
    end
    @printf("\n  worst reldiff = %.3e\n", worst)
    println(nfail == 0 ? "  WHOLE dust_dry_dep! MATCHES CPU ON $name ($FT)" :
                         "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
