# ==========================================================================
# gpu_validate_infexcessrunoff_e2e.jl — end-to-end Metal parity for the WHOLE
# infiltration_excess_runoff! routine (both per-column kernels: the
# compute_qinmax_hksat! top-3-layer min reduction AND the column-average +
# Hortonian-excess kernel).
#
# Builds a small multi-column Float32 instance (exercising the major branches:
# zero ice, partial ice, fully-frozen soil, partial fsat, partial frac_h2osfc,
# and a masked-out column), runs infiltration_excess_runoff! on the CPU, adapts
# every state struct (+ masks/inputs) to Metal, runs the SAME call on the device,
# and compares the mutated outputs (qinmax_col, qflx_infl_excess_col) field by
# field. Also exercises the QINMAX_METHOD_NONE path in a separate device run.
#
#   julia --project=scripts scripts/gpu_validate_infexcessrunoff_e2e.jl
# ==========================================================================

using CLM
using Printf
import Metal   # Metal-specific; MtlArray is the Adapt adaptor type for the structs
include(joinpath(@__DIR__, "gpu_backends.jl"))

# reldiff: NaN-aware (both-NaN agrees; one-sided NaN flags divergence). Also asserts
# (via cpu_has_finite at the call site) that the CPU reference is FINITE so a
# both-NaN false PASS can't slip through.
function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        # Matching infinities agree exactly (the NONE qinmax_unlimited path yields
        # Inf qinmax on a Float32 backend; abs(Inf-Inf) would be NaN otherwise).
        (isinf(A[i]) && isinf(B[i]) && sign(A[i]) == sign(B[i])) && continue
        m = max(m, abs(A[i] - B[i]) / (1.0 + max(abs(A[i]), abs(B[i]))))
    end
    return m
end
# True if the field carries any usable (non-NaN) parity signal: a finite value OR a
# matching infinity (qinmax is legitimately Inf on the NONE / Float32 path).
cpu_has_signal(a) = any(x -> isfinite(x) || isinf(x), Array(a))

# Build a state covering the major branches. `method` selects HKSAT (default) or
# NONE (qinmax_unlimited path). `nc` columns; column 5 is masked OUT.
function build(::Type{FT}, method::Int) where {FT}
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
    nc = 5
    nlevgrnd = CLM.varpar.nlevgrnd

    ier = CLM.InfiltrationExcessRunoffData(qinmax_col = Vector{FT}(undef, 0))
    CLM.infilt_excess_runoff_init!(ier, nc)
    ier.qinmax_method = method

    sh = CLM.SoilHydrologyData{FT}(); CLM.soilhydrology_init!(sh, nc)
    ss = CLM.SoilStateData{FT}();     CLM.soilstate_init!(ss, nc, nc)

    # Initialize the full ice/hksat arrays to finite values (avoid NaN/SPVAL
    # leaking into the top-3-layer min reduction).
    sh.icefrac_col .= FT(0)
    ss.hksat_col   .= FT(0.01)

    # Column 1: zero ice, uniform hksat (max infiltration branch).
    # Column 2: high conductivity.
    for j in 1:3; ss.hksat_col[2, j] = FT(0.1); end
    # Column 3: partial ice (impedance branch).
    for j in 1:3; sh.icefrac_col[3, j] = FT(0.3); end
    # Column 4: fully frozen top 3 layers (minimum infiltration branch).
    for j in 1:3; sh.icefrac_col[4, j] = FT(1.0); end
    # Column 5: masked out (left at init values; must stay NaN).

    fsat = FT[0.1, 0.2, 0.0, 0.0, 0.5]

    wfb = CLM.WaterFluxBulkData{FT}(); CLM.waterfluxbulk_init!(wfb, nc, nc, 1, 1)
    wfb.qflx_in_soil_col .= FT[0.02, 0.005, 0.005, 0.001, 0.005]
    wfb.qflx_infl_excess_col .= FT(NaN)

    wdb = CLM.WaterDiagnosticBulkData{FT}(); CLM.waterdiagnosticbulk_init!(wdb, nc, nc, 1, 1)
    wdb.frac_h2osfc_col .= FT[0.0, 0.0, 0.0, 0.0, 0.5]

    mask = BitVector([true, true, true, true, false])

    params = CLM.InfiltrationExcessRunoffParams(e_ice = 6.0)
    return (; nc, ier, sh, ss, fsat, wfb, wdb, mask, params)
end

run_ier!(S) = CLM.infiltration_excess_runoff!(
    S.ier, S.sh, S.ss, S.fsat, S.wfb, S.wdb, S.mask, 1:S.nc; params = S.params)

# Device run: structs are adapted by the caller; the device variant takes already-
# device-resident states/masks.
function run_ier_dev!(ier, sh, ss, fsat, wfb, wdb, mask, nc, params)
    CLM.infiltration_excess_runoff!(ier, sh, ss, fsat, wfb, wdb, mask, 1:nc; params = params)
end

function validate(name, to, FT, method)
    methodname = method == CLM.QINMAX_METHOD_NONE ? "NONE" : "HKSAT"
    @printf("\n--- qinmax_method = %s ---\n", methodname)

    H = build(FT, method)   # CPU reference
    B = build(FT, method)   # source for the device snapshot

    ad(x) = CLM.Adapt.adapt(Metal.MtlArray, x)
    dmask(m) = to(collect(Bool, m))

    sh_d  = ad(B.sh)
    ss_d  = ad(B.ss)
    wfb_d = ad(B.wfb)
    wdb_d = ad(B.wdb)
    ier_d = ad(B.ier)
    fsat_d = to(B.fsat)
    mask_d = dmask(B.mask)

    if !(wfb_d.qflx_infl_excess_col isa Metal.MtlArray)
        println("  BLOCKED: a state struct did not move to the device under adapt.")
        return 2
    end

    run_ier!(H)                                                          # CPU
    run_ier_dev!(ier_d, sh_d, ss_d, fsat_d, wfb_d, wdb_d, mask_d, B.nc, B.params)  # device

    checks = [
        ("qinmax_col",       H.ier.qinmax_col,             ier_d.qinmax_col),
        ("qflx_infl_excess", H.wfb.qflx_infl_excess_col,   wfb_d.qflx_infl_excess_col),
    ]
    nfail = 0
    for (nm, a, b) in checks
        if !cpu_has_signal(a)
            @printf("  [WARN] %-18s CPU reference is all-NaN — no parity signal\n", nm)
            nfail += 1
            continue
        end
        d = reldiff(a, b); ok = d < 1f-3
        @printf("  [%s] %-18s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
        ok || (nfail += 1)
    end

    # The masked-out column 5 must remain untouched (NaN on both backends).
    q_cpu = Array(H.ier.qinmax_col)
    q_dev = Array(ier_d.qinmax_col)
    masked_ok = isnan(q_cpu[5]) && isnan(q_dev[5])
    @printf("  [%s] masked-out col 5 untouched (NaN both sides)\n", masked_ok ? "PASS" : "FAIL")
    masked_ok || (nfail += 1)

    return nfail
end

function main(backend)
    println("=" ^ 70)
    println("END-TO-END Metal parity for infiltration_excess_runoff! (whole routine)")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU routine exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n", name, FT)

    nfail = 0
    nfail += validate(name, to, FT, CLM.QINMAX_METHOD_HKSAT)
    nfail += validate(name, to, FT, CLM.QINMAX_METHOD_NONE)

    println()
    println(nfail == 0 ? "  WHOLE infiltration_excess_runoff! MATCHES CPU ON $name ($FT) ✓" :
                         "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
