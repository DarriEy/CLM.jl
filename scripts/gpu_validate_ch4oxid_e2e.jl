# ==========================================================================
# gpu_validate_ch4oxid_e2e.jl — end-to-end GPU parity for the WHOLE ch4_oxid!
# function (the single per-(column, soil layer) _ch4oxid_kernel! computing
# ch4_oxid_depth / o2_oxid_depth via double Michaelis-Menten kinetics).
#
# Builds the test_methane.jl CH4Data + ch4_oxid! arg setup at Float32, runs the
# CPU reference, then drives the SAME kernel on the GPU device via the loose-
# array launcher CLM.meth_oxid! with every input array adapted to device_array_type(), and
# compares the mutated outputs (ch4_oxid_depth, o2_oxid_depth) with reldiff.
#
# Exercises all the branches of ch4_oxid!:
#   * sat == 1   (saturated; all layers use the sat k_m / vmax; below-WT Henry)
#   * sat == 0   (unsaturated; layers above jwt switch to the unsat k_m/vmax and
#                 the >watsat Henry branch, below-WT use the smp_fact path)
#   * a frozen column (t_soisno <= TFRZ) so the oxid_a := 0 branch fires
#
# Note (worktree foundation): on this branch CH4Data is NOT yet @adapt_structure'd
# (Float64 fields), so the device path drives the loose-array kernel launcher
# CLM.meth_oxid! directly with Float32 device_array_type() inputs — exactly the arrays the
# kernelized ch4_oxid! extracts from the struct. This validates the kernel itself.
#
# IMPORTANT (A1 lesson): reldiff silently PASSES when both sides are NaN, so this
# harness first asserts the CPU-F32 reference outputs are FINITE before trusting
# any parity number.
#
#   julia --project=scripts scripts/gpu_validate_ch4oxid_e2e.jl
# ==========================================================================

using CLM
using Printf
include(joinpath(@__DIR__, "gpu_backends.jl"))

function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(A[i] - B[i]) / (1.0 + max(abs(A[i]), abs(B[i]))))
    end
    return m
end

allfinite(a) = all(isfinite, Array(a))

# Build the inputs at precision FT. Returns the CH4Data, the loose input arrays,
# the per-column jwt, and dims. jwt is chosen so both above-WT and below-WT
# layers exist in the unsat run (column water table partway down the profile).
function build(::Type{FT}) where {FT}
    nc = 3; nlevsoi = 5

    ch4 = CLM.CH4Data{FT}(
        ch4_oxid_depth_sat_col      = zeros(FT, nc, nlevsoi),
        ch4_oxid_depth_unsat_col    = zeros(FT, nc, nlevsoi),
        o2_oxid_depth_sat_col       = zeros(FT, nc, nlevsoi),
        o2_oxid_depth_unsat_col     = zeros(FT, nc, nlevsoi),
        conc_ch4_sat_col            = fill(FT(1.0e-4), nc, nlevsoi),
        conc_ch4_unsat_col          = fill(FT(1.0e-5), nc, nlevsoi),
        conc_o2_sat_col             = fill(FT(0.01),   nc, nlevsoi),
        conc_o2_unsat_col           = fill(FT(0.05),   nc, nlevsoi),
    )

    params = CLM.CH4Params()

    mask_soil = trues(nc)
    watsat     = fill(FT(0.45), nc, nlevsoi)
    h2osoi_vol = fill(FT(0.30), nc, nlevsoi)
    smp_l      = fill(FT(-1000.0), nc, nlevsoi)
    # column 3 frozen so the oxid_a := 0 branch fires
    t_soisno = fill(FT(CLM.TFRZ) + FT(15.0), nc, nlevsoi)
    for j in 1:nlevsoi
        t_soisno[3, j] = FT(CLM.TFRZ) - FT(5.0)
    end
    # water table partway down the profile so each unsat column spans both
    # the above-WT (unsat k_m/vmax, >watsat Henry) and below-WT branches.
    jwt = Int[2, 3, 2]

    return (; nc, nlevsoi, ch4, params, mask_soil,
            watsat, h2osoi_vol, smp_l, t_soisno, jwt)
end

function main(backend)
    println("=" ^ 70)
    println("END-TO-END GPU parity for ch4_oxid! (per-(column,layer) kernel)")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU driver exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    nfail = 0
    # Both config branches: sat==0 (unsaturated) and sat==1 (saturated).
    for sat in (0, 1)
        lake = false
        B = build(FT)
        nc, nlevsoi = B.nc, B.nlevsoi
        p = B.params

        # --- CPU reference: the whole kernelized ch4_oxid! on a Float32 CH4Data ---
        CLM.ch4_oxid!(B.ch4, p, B.mask_soil, B.watsat, B.h2osoi_vol,
                      B.smp_l, B.t_soisno, B.jwt, sat, lake, nlevsoi, FT(1800.0))

        # Pick the variant this `sat` writes.
        if sat == 0
            cpu_ch4 = B.ch4.ch4_oxid_depth_unsat_col
            cpu_o2  = B.ch4.o2_oxid_depth_unsat_col
            conc_ch4 = B.ch4.conc_ch4_unsat_col
            conc_o2  = B.ch4.conc_o2_unsat_col
        else
            cpu_ch4 = B.ch4.ch4_oxid_depth_sat_col
            cpu_o2  = B.ch4.o2_oxid_depth_sat_col
            conc_ch4 = B.ch4.conc_ch4_sat_col
            conc_o2  = B.ch4.conc_o2_sat_col
        end

        if !allfinite(cpu_ch4) || !allfinite(cpu_o2)
            @printf("  BLOCKED [sat=%d]: CPU-F32 reference output non-finite — parity untrustworthy.\n", sat)
            return 2
        end

        # --- Device: drive the SAME kernel via the loose-array launcher ---
        d_ch4 = to(zeros(FT, nc, nlevsoi))
        d_o2  = to(zeros(FT, nc, nlevsoi))
        dmask = to(collect(Bool, B.mask_soil))
        d_jwt = to(collect(Int, B.jwt))
        d_watsat = to(B.watsat); d_h2osoi = to(B.h2osoi_vol)
        d_smp = to(B.smp_l); d_tsoi = to(B.t_soisno)
        d_cch4 = to(conc_ch4); d_co2 = to(conc_o2)

        if !(d_ch4 isa device_array_type())
            println("  BLOCKED: output array did not move to the device under to().")
            return 2
        end

        CLM.meth_oxid!(d_ch4, d_o2, dmask, d_jwt, d_watsat, d_h2osoi, d_smp,
                       d_tsoi, d_cch4, d_co2, sat,
                       p.vmax_ch4_oxid, p.k_m, p.q10_ch4oxid, p.smp_crit,
                       p.k_m_o2, p.k_m_unsat, p.vmax_oxid_unsat, nc, nlevsoi)

        for (nm, a, b) in (("ch4_oxid_depth", cpu_ch4, d_ch4),
                           ("o2_oxid_depth",  cpu_o2,  d_o2))
            d = reldiff(a, b); ok = d < 1f-2
            @printf("  [%s] sat=%d %-15s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", sat, nm, d)
            ok || (nfail += 1)
        end
    end

    println()
    println(nfail == 0 ? "  WHOLE ch4_oxid! MATCHES CPU ON $name ($FT) ✓" :
                         "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
