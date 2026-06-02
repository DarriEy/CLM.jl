# ==========================================================================
# gpu_validate_roughlengths_e2e.jl — end-to-end Metal parity for the WHOLE
# set_actual_roughness_lengths! driver (per-patch priority selection of the
# momentum roughness length: exposed-veg > no-exposed-veg > urban > lake).
#
# Builds a multi-patch FrictionVelocityData at Float32 with one patch per
# priority class (so every branch of the kernel is exercised), runs
# set_actual_roughness_lengths! on the CPU, adapts fv + the masks/index vectors
# to Metal, runs the SAME call on the device, and compares z0m_actual_patch.
#
# CRITICAL: source roughness lengths are real defaults so the CPU reference is
# FINITE; we assert finiteness before trusting parity.
#
#   julia --project=scripts scripts/gpu_validate_roughlengths_e2e.jl
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
    np = 4; nc = 4; nl = 4
    fv = CLM.FrictionVelocityData{FT}()
    CLM.frictionvel_init!(fv, np, nc)
    for p in 1:np
        fv.z0mv_patch[p]       = FT(0.10 + 0.01 * p)   # veg roughness
        fv.z0m_actual_patch[p] = FT(NaN)
    end
    for c in 1:nc
        fv.z0mg_col[c] = FT(0.001 * c)                 # ground roughness
    end
    lun_z_0_town = FT[0.5, 0.6, 0.7, 0.8]              # urban canyon roughness

    patch_column   = [1, 2, 3, 4]
    patch_landunit = [1, 2, 3, 4]

    # one patch per priority class
    m_exposedveg   = falses(np); m_exposedveg[1]   = true
    m_noexposedveg = falses(np); m_noexposedveg[2] = true
    m_urban        = falses(np); m_urban[3]        = true
    m_lake         = falses(np); m_lake[4]         = true

    return (; np, fv, lun_z_0_town, patch_column, patch_landunit,
              m_exposedveg, m_noexposedveg, m_urban, m_lake)
end

run_sar!(H, fv, mev, mnev, mu, ml, col, lun, ltown) =
    CLM.set_actual_roughness_lengths!(fv, mev, mnev, mu, ml, 1:H.np, col, lun, ltown)

function main(backend)
    println("=" ^ 70)
    println("END-TO-END Metal parity for set_actual_roughness_lengths! (whole driver)")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU driver exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    H = build(FT)
    B = build(FT)

    ad(x) = CLM.Adapt.adapt(Metal.MtlArray, x)
    fv_d   = ad(B.fv)
    ltown_d = ad(B.lun_z_0_town)
    col_d  = to(collect(Int32, B.patch_column))
    lun_d  = to(collect(Int32, B.patch_landunit))
    mev_d  = to(collect(Bool, B.m_exposedveg))
    mnev_d = to(collect(Bool, B.m_noexposedveg))
    mu_d   = to(collect(Bool, B.m_urban))
    ml_d   = to(collect(Bool, B.m_lake))

    if !(fv_d.z0m_actual_patch isa Metal.MtlArray)
        println("  BLOCKED: FrictionVelocityData did not move to the device under adapt.")
        return 2
    end

    run_sar!(H, H.fv, H.m_exposedveg, H.m_noexposedveg, H.m_urban, H.m_lake,
             H.patch_column, H.patch_landunit, H.lun_z_0_town)
    run_sar!(B, fv_d, mev_d, mnev_d, mu_d, ml_d, col_d, lun_d, ltown_d)

    if !all(isfinite, Array(H.fv.z0m_actual_patch))
        println("  BLOCKED: CPU z0m_actual_patch not finite — parity would be a false PASS.")
        return 2
    end

    checks = [("z0m_actual", H.fv.z0m_actual_patch, fv_d.z0m_actual_patch)]
    nfail = 0; worst = 0.0
    for (nm, a, b) in checks
        d = reldiff(a, b); ok = d < 1f-2
        worst = max(worst, d)
        @printf("  [%s] %-15s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
        ok || (nfail += 1)
    end
    @printf("\n  worst reldiff = %.3e\n", worst)
    println(nfail == 0 ? "  WHOLE set_actual_roughness_lengths! MATCHES CPU ON $name ($FT)" :
                         "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
