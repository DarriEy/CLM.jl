# ==========================================================================
# gpu_validate_decprecision_e2e.jl — end-to-end GPU parity for the WHOLE
# soil_bgc_precision_control! BGC routine (decomposition C/N precision control).
#
# Builds a small multi-column / multi-level / multi-pool instance mirroring
# test/test_decomp_precision_control.jl, runs soil_bgc_precision_control! on the
# CPU, converts every state struct to Float32 + adapts to the GPU, runs the SAME
# call on the device, and compares the mutated outputs field-by-field. The
# scenario deliberately exercises the branchy config paths in one shot:
#   * the base C/N pool truncation + own-index ctrunc/ntrunc sinks
#   * the C13 + C14 isotope truncation (use_c13 / use_c14)
#   * the use_nitrif_denitrif mineral-N (NO3/NH4) small-negative reset block
#   * a masked-out column (mask filtering)
#
#   julia --project=scripts scripts/gpu_validate_decprecision_e2e.jl
#
# The CN state structs are concretely Float64-typed, so we build at Float64 and
# down-convert array fields to Float32 for the device snapshot (Metal has no
# Float64).
# ==========================================================================

using CLM
using Printf
include(joinpath(@__DIR__, "gpu_backends.jl"))

# NaN-aware mixed abs/rel diff; asserts the CPU reference is finite so a both-NaN
# false PASS can't slip through.
function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(A[i] - B[i]) / (1.0 + max(abs(A[i]), abs(B[i]))))
    end
    return m
end
cpu_has_finite(a) = any(isfinite, Array(a))

# Float32-down-converting Metal adaptor. The CN state structs are concretely
# Float64-typed and *_init! fills Float64, so we adapt with a custom storage rule
# that down-converts float arrays to Float32 while reconstructing the struct.
struct MetalF32 end
CLM.Adapt.adapt_storage(::MetalF32, x::AbstractArray{<:AbstractFloat}) = device_array_type()(Float32.(x))
CLM.Adapt.adapt_storage(::MetalF32, x::AbstractArray{<:Integer})       = device_array_type()(x)
CLM.Adapt.adapt_storage(::MetalF32, x::AbstractArray{Bool})            = device_array_type()(x)

const NC = 4
const NLEVDECOMP = 3
const NDECOMP_POOLS = 7

# Build the CPU reference state (Float64). Returns a NamedTuple of the state
# structs (cs, ns, c13cs, c14cs) + the column mask.
function build()
    cs = CLM.SoilBiogeochemCarbonStateData()
    ns = CLM.SoilBiogeochemNitrogenStateData()

    cs.decomp_cpools_vr_col = zeros(NC, NLEVDECOMP, NDECOMP_POOLS)
    cs.ctrunc_vr_col        = zeros(NC, NLEVDECOMP)

    ns.decomp_npools_vr_col = zeros(NC, NLEVDECOMP, NDECOMP_POOLS)
    ns.ntrunc_vr_col        = zeros(NC, NLEVDECOMP)
    ns.smin_no3_vr_col      = zeros(NC, NLEVDECOMP)
    ns.smin_nh4_vr_col      = zeros(NC, NLEVDECOMP)

    c13cs = CLM.SoilBiogeochemCarbonStateData()
    c13cs.decomp_cpools_vr_col = zeros(NC, NLEVDECOMP, NDECOMP_POOLS)
    c13cs.ctrunc_vr_col        = zeros(NC, NLEVDECOMP)

    c14cs = CLM.SoilBiogeochemCarbonStateData()
    c14cs.decomp_cpools_vr_col = zeros(NC, NLEVDECOMP, NDECOMP_POOLS)
    c14cs.ctrunc_vr_col        = zeros(NC, NLEVDECOMP)

    # --- Fill pools with a mix of below- and above-threshold values, plus a
    #     pre-existing truncation sink, so every branch moves something. ---
    tiny = 5.0e-9   # below ccrit = 1e-8
    for c in 1:NC, j in 1:NLEVDECOMP, k in 1:NDECOMP_POOLS
        # alternate sign + magnitude; some tiny (truncate), some large (keep)
        base = ((c + j + k) % 3 == 0) ? 1.0 : tiny * (iseven(k) ? -1.0 : 1.0)
        cs.decomp_cpools_vr_col[c, j, k] = base
        ns.decomp_npools_vr_col[c, j, k] = base * 0.1
        c13cs.decomp_cpools_vr_col[c, j, k] = base * 0.011
        c14cs.decomp_cpools_vr_col[c, j, k] = base * 1.2e-12
    end
    # pre-existing truncation accumulators (exercise the += into existing values)
    cs.ctrunc_vr_col    .= 1.0e-5
    ns.ntrunc_vr_col    .= 2.0e-5
    c13cs.ctrunc_vr_col .= 3.0e-7
    c14cs.ctrunc_vr_col .= 4.0e-19

    # --- Mineral N: small negatives (reset), small positives (keep), large (keep) ---
    for c in 1:NC, j in 1:NLEVDECOMP
        ns.smin_no3_vr_col[c, j] = (j == 1) ? -5.0e-13 : (j == 2 ? 5.0e-13 : -5.0e-8)
        ns.smin_nh4_vr_col[c, j] = (j == 1) ? -3.0e-13 : (j == 2 ? 3.0e-13 : -3.0e-8)
    end

    # Column 3 masked out (verify mask filtering on device too)
    mask = trues(NC)
    mask[3] = false

    return (; cs, ns, c13cs, c14cs, mask)
end

run_dpc!(cs, ns, mask; c13cs, c14cs) =
    CLM.soil_bgc_precision_control!(cs, ns;
        mask_bgc_soilc=mask,
        nlevdecomp=NLEVDECOMP,
        ndecomp_pools=NDECOMP_POOLS,
        c13cs=c13cs, c14cs=c14cs,
        use_c13=true, use_c14=true,
        use_nitrif_denitrif=true, use_fun=false)

function main(backend)
    println("=" ^ 70)
    println("END-TO-END GPU parity for soil_bgc_precision_control! (BGC C/N precision)")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU path exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    H = build()    # CPU reference (Float64)
    B = build()    # source for the device snapshot

    # Adapt the device snapshot to Metal, down-converting float arrays to Float32.
    ad(x) = CLM.Adapt.adapt(MetalF32(), x)
    cs_d    = ad(B.cs)
    ns_d    = ad(B.ns)
    c13cs_d = ad(B.c13cs)
    c14cs_d = ad(B.c14cs)

    if !(cs_d.decomp_cpools_vr_col isa device_array_type())
        println("  BLOCKED: a state struct did not move to the device under adapt.")
        return 2
    end

    dmask(m) = to(collect(Bool, m))

    # CPU
    run_dpc!(H.cs, H.ns, H.mask; c13cs=H.c13cs, c14cs=H.c14cs)
    # Device
    run_dpc!(cs_d, ns_d, dmask(B.mask); c13cs=c13cs_d, c14cs=c14cs_d)

    checks = [
        ("decomp_cpools",  H.cs.decomp_cpools_vr_col,    cs_d.decomp_cpools_vr_col),
        ("ctrunc",         H.cs.ctrunc_vr_col,           cs_d.ctrunc_vr_col),
        ("decomp_npools",  H.ns.decomp_npools_vr_col,    ns_d.decomp_npools_vr_col),
        ("ntrunc",         H.ns.ntrunc_vr_col,           ns_d.ntrunc_vr_col),
        ("c13 decomp_cp",  H.c13cs.decomp_cpools_vr_col, c13cs_d.decomp_cpools_vr_col),
        ("c13 ctrunc",     H.c13cs.ctrunc_vr_col,        c13cs_d.ctrunc_vr_col),
        ("c14 decomp_cp",  H.c14cs.decomp_cpools_vr_col, c14cs_d.decomp_cpools_vr_col),
        ("c14 ctrunc",     H.c14cs.ctrunc_vr_col,        c14cs_d.ctrunc_vr_col),
        ("smin_no3",       H.ns.smin_no3_vr_col,         ns_d.smin_no3_vr_col),
        ("smin_nh4",       H.ns.smin_nh4_vr_col,         ns_d.smin_nh4_vr_col),
    ]
    nfail = 0
    for (nm, a, b) in checks
        if !cpu_has_finite(a)
            @printf("  [WARN] %-16s CPU reference is all-NaN/Inf — skipping (no parity signal)\n", nm)
            continue
        end
        d = reldiff(a, b); ok = d < 1f-3
        @printf("  [%s] %-16s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
        ok || (nfail += 1)
    end
    println()
    println(nfail == 0 ? "  WHOLE soil_bgc_precision_control! MATCHES CPU ON $name ($FT) ✓" :
                         "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
