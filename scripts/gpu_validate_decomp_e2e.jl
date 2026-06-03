# ==========================================================================
# gpu_validate_decomp_e2e.jl — end-to-end Metal parity for the WHOLE
# soil_biogeochem_decomp! routine (every per-(c)/internal-l/j/k loop now a
# KernelAbstractions kernel: the C:N-ratio loop, the decomposition-cascade
# application, the methane fphr loop, and the vertical integration).
#
# Builds a small instance mirroring test/test_soil_biogeochem_decomp.jl's BGC
# setup at Float64, runs soil_biogeochem_decomp! on the CPU, down-converts every
# array to Float32 / Metal (MetalF32 adaptor below), runs the SAME call on the
# device, and compares the mutated outputs field-by-field with a NaN-aware
# reldiff. Exercises the branchy config paths (use_nitrif_denitrif / use_lch4 /
# use_mimics-off) plus the masked-column skip.
#
#   julia --project=scripts scripts/gpu_validate_decomp_e2e.jl
# ==========================================================================

using CLM
using Printf
import Metal   # MtlArray is the Adapt adaptor type for the device arrays
const Adapt = CLM.Adapt   # Adapt isn't in scripts/Project.toml; reach it via CLM
import KernelAbstractions as KA
include(joinpath(@__DIR__, "gpu_backends.jl"))

# reldiff: NaN-aware (both-NaN agrees; one-sided NaN flags divergence).
function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(Float64(A[i]) - Float64(B[i])) / (1.0 + max(abs(Float64(A[i])), abs(Float64(B[i])))))
    end
    return m
end
cpu_has_finite(a) = any(isfinite, Array(a))

# --------------------------------------------------------------------------
# MetalF32 down-convert adaptor: the decomp state structs are built at Float64
# (the oracle test uses default constructors). Metal has no hardware Float64, so
# every Float64 array is down-converted to Float32 on the device; integer index
# vectors keep an integer eltype (Int32 — never coerced to float, which would
# break device indexing); BitVector masks/flags become Bool device arrays.
# --------------------------------------------------------------------------
struct MetalF32 end
Adapt.adapt_storage(::MetalF32, x::AbstractArray{<:AbstractFloat}) =
    Metal.MtlArray(Float32.(Array(x)))
Adapt.adapt_storage(::MetalF32, x::AbstractArray{<:Integer}) =
    Metal.MtlArray(Int32.(Array(x)))
Adapt.adapt_storage(::MetalF32, x::BitArray) =
    Metal.MtlArray(collect(Bool, x))
Adapt.adapt_storage(::MetalF32, x::AbstractArray{Bool}) =
    Metal.MtlArray(collect(Bool, x))
adf(x) = Adapt.adapt(MetalF32(), x)

# --------------------------------------------------------------------------
# Build a BGC decomposition instance (clone of test/test_soil_biogeochem_decomp.jl
# make_decomp_data), at Float64.
# --------------------------------------------------------------------------
function build(; nc=4, nlevdecomp=3, ndecomp_pools=7, ndecomp_cascade_transitions=10)
    bounds = 1:nc
    nlevdecomp_full = nlevdecomp

    cs = CLM.SoilBiogeochemCarbonStateData()
    CLM.soil_bgc_carbon_state_init!(cs, nc, 1, nlevdecomp_full, ndecomp_pools)
    for c in 1:nc, j in 1:nlevdecomp, l in 1:ndecomp_pools
        cs.decomp_cpools_vr_col[c, j, l] = 100.0 + 10.0 * l
    end

    ns = CLM.SoilBiogeochemNitrogenStateData()
    CLM.soil_bgc_nitrogen_state_init!(ns, nc, 1, nlevdecomp_full, ndecomp_pools)
    for c in 1:nc, j in 1:nlevdecomp, l in 1:ndecomp_pools
        ns.decomp_npools_vr_col[c, j, l] = cs.decomp_cpools_vr_col[c, j, l] / 12.0
    end

    cf = CLM.SoilBiogeochemCarbonFluxData()
    CLM.soil_bgc_carbon_flux_init!(cf, nc, nlevdecomp_full, ndecomp_pools,
                                    ndecomp_cascade_transitions)
    for c in 1:nc, j in 1:nlevdecomp, k in 1:ndecomp_cascade_transitions
        cf.rf_decomp_cascade_col[c, j, k] = 0.5
        cf.c_overflow_vr[c, j, k] = 0.0
    end
    for c in 1:nc, j in 1:nlevdecomp
        cf.w_scalar_col[c, j] = 0.8
        cf.phr_vr_col[c, j] = 1.0e-5
    end

    nf = CLM.SoilBiogeochemNitrogenFluxData()
    CLM.soil_bgc_nitrogen_flux_init!(nf, nc, nlevdecomp_full, ndecomp_pools,
                                      ndecomp_cascade_transitions)
    for c in 1:nc
        nf.net_nmin_col[c] = 0.0
        nf.gross_nmin_col[c] = 0.0
        for j in 1:nlevdecomp
            nf.net_nmin_vr_col[c, j] = 0.0
            nf.gross_nmin_vr_col[c, j] = 0.0
        end
    end

    st = CLM.SoilBiogeochemStateData()
    CLM.soil_bgc_state_init!(st, nc, nc, nlevdecomp_full, ndecomp_cascade_transitions)
    for c in 1:nc, j in 1:nlevdecomp
        st.fpi_vr_col[c, j] = 0.5
    end

    params_bgc = CLM.DecompBGCParams(bgc_initial_Cstocks = fill(200.0, ndecomp_pools))
    cn_params = CLM.CNSharedParamsData()
    bgc_state = CLM.DecompBGCState()
    cascade_con = CLM.DecompCascadeConData()
    cellsand = fill(50.0, nc, max(nlevdecomp, 5))
    CLM.init_decomp_cascade_bgc!(
        bgc_state, cascade_con, params_bgc, cn_params;
        cellsand=cellsand, bounds=bounds, nlevdecomp=nlevdecomp,
        ndecomp_pools_max=ndecomp_pools,
        ndecomp_cascade_transitions_max=ndecomp_cascade_transitions,
        use_fates=false)

    decomp_params = CLM.DecompParams(dnp=0.01)
    mask_bgc_soilc = trues(nc)

    cn_decomp_pools = zeros(nc, nlevdecomp, ndecomp_pools)
    p_decomp_cpool_loss = fill(1.0e-6, nc, nlevdecomp, ndecomp_cascade_transitions)
    pmnf_decomp_cascade = fill(-1.0e-7, nc, nlevdecomp, ndecomp_cascade_transitions)
    p_decomp_npool_to_din = zeros(nc, nlevdecomp, ndecomp_cascade_transitions)
    dzsoi_decomp = fill(0.1, nlevdecomp)

    # Non-zero gross_nmin_vr so the vertical-integration kernel has a signal.
    for c in 1:nc, j in 1:nlevdecomp
        nf.gross_nmin_vr_col[c, j] = 1.0e-6
    end

    return (; cf, cs, nf, ns, st, cascade_con, decomp_params, mask_bgc_soilc,
              bounds, nlevdecomp, ndecomp_pools, ndecomp_cascade_transitions,
              cn_decomp_pools, p_decomp_cpool_loss, pmnf_decomp_cascade,
              p_decomp_npool_to_din, dzsoi_decomp, nc)
end

# Run on CPU (mask is a BitVector) — H is the in-place state bundle.
function run_cpu!(H; use_nitrif_denitrif, use_lch4, use_mimics)
    CLM.soil_biogeochem_decomp!(
        H.cf, H.cs, H.nf, H.ns, H.st, H.cascade_con, H.decomp_params;
        mask_bgc_soilc=H.mask_bgc_soilc, bounds=H.bounds, nlevdecomp=H.nlevdecomp,
        ndecomp_pools=H.ndecomp_pools,
        ndecomp_cascade_transitions=H.ndecomp_cascade_transitions,
        cn_decomp_pools=H.cn_decomp_pools, p_decomp_cpool_loss=H.p_decomp_cpool_loss,
        pmnf_decomp_cascade=H.pmnf_decomp_cascade,
        p_decomp_npool_to_din=H.p_decomp_npool_to_din, dzsoi_decomp=H.dzsoi_decomp,
        use_nitrif_denitrif=use_nitrif_denitrif, use_lch4=use_lch4,
        use_mimics=use_mimics)
end

# Build a device snapshot: down-convert every state struct + work array + the
# cascade-con index/flag vectors, run the SAME call on the device.
function run_dev(B; use_nitrif_denitrif, use_lch4, use_mimics)
    cf  = adf(B.cf); cs = adf(B.cs); nf = adf(B.nf); ns = adf(B.ns); st = adf(B.st)
    # cascade_con stays host-resident: it holds String metadata that cannot move
    # to a device, and soil_biogeochem_decomp! itself copies its index/flag/cn
    # vectors onto the state backend (preserving Int/Bool eltype). The prototype
    # for that copy is cn_decomp_pools, which IS a device array below.
    cc = B.cascade_con
    cn  = adf(B.cn_decomp_pools)
    pcl = adf(B.p_decomp_cpool_loss)
    pmnf = adf(B.pmnf_decomp_cascade)
    pdin = adf(B.p_decomp_npool_to_din)
    dz   = adf(B.dzsoi_decomp)
    mask = adf(B.mask_bgc_soilc)

    CLM.soil_biogeochem_decomp!(
        cf, cs, nf, ns, st, cc, B.decomp_params;
        mask_bgc_soilc=mask, bounds=B.bounds, nlevdecomp=B.nlevdecomp,
        ndecomp_pools=B.ndecomp_pools,
        ndecomp_cascade_transitions=B.ndecomp_cascade_transitions,
        cn_decomp_pools=cn, p_decomp_cpool_loss=pcl, pmnf_decomp_cascade=pmnf,
        p_decomp_npool_to_din=pdin, dzsoi_decomp=dz,
        use_nitrif_denitrif=use_nitrif_denitrif, use_lch4=use_lch4,
        use_mimics=use_mimics)

    return (; cf, cs, nf, ns, st, cn, pcl, pmnf, pdin)
end

function compare(H, D, label)
    checks = [
        ("cn_decomp_pools",        H.cn_decomp_pools,                       D.cn),
        ("p_decomp_cpool_loss",    H.p_decomp_cpool_loss,                   D.pcl),
        ("pmnf_decomp_cascade",    H.pmnf_decomp_cascade,                   D.pmnf),
        ("decomp_cascade_hr_vr",   H.cf.decomp_cascade_hr_vr_col,           D.cf.decomp_cascade_hr_vr_col),
        ("decomp_cascade_ctrans",  H.cf.decomp_cascade_ctransfer_vr_col,    D.cf.decomp_cascade_ctransfer_vr_col),
        ("decomp_cascade_ntrans",  H.nf.decomp_cascade_ntransfer_vr_col,    D.nf.decomp_cascade_ntransfer_vr_col),
        ("decomp_cascade_sminn",   H.nf.decomp_cascade_sminn_flux_vr_col,   D.nf.decomp_cascade_sminn_flux_vr_col),
        ("sminn_to_denit",         H.nf.sminn_to_denit_decomp_cascade_vr_col, D.nf.sminn_to_denit_decomp_cascade_vr_col),
        ("net_nmin_vr",            H.nf.net_nmin_vr_col,                    D.nf.net_nmin_vr_col),
        ("net_nmin",               H.nf.net_nmin_col,                       D.nf.net_nmin_col),
        ("gross_nmin",             H.nf.gross_nmin_col,                     D.nf.gross_nmin_col),
        ("fphr",                   H.cf.fphr_col,                           D.cf.fphr_col),
    ]
    nfail = 0
    for (nm, a, b) in checks
        if !cpu_has_finite(a)
            @printf("    [WARN] %-22s CPU ref all-NaN/Inf — skipping\n", nm)
            continue
        end
        d = reldiff(a, b); ok = d < 1f-3
        @printf("    [%s] %-22s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
        ok || (nfail += 1)
    end
    println("    ", nfail == 0 ? "$label MATCHES ✓" : "$label DIVERGENCE")
    return nfail
end

function main(backend)
    println("=" ^ 70)
    println("END-TO-END Metal parity for soil_biogeochem_decomp! (whole routine)")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU routine exercised by the suite).")
        return 0
    end
    name, _to, FT = backend
    @printf("  Backend: %s   (device precision: %s)\n\n", name, FT)
    if name != "Metal"
        println("  This harness targets the Metal Float32 down-convert path; backend is $name. Skipping.")
        return 0
    end

    total = 0

    # --- Config A: default (use_nitrif_denitrif=false, no lch4, no mimics) ---
    println("  [A] default config (simple denitrification path, pmnf<0)")
    H = build(); B = build()
    fill!(H.pmnf_decomp_cascade, -1.0e-7); fill!(B.pmnf_decomp_cascade, -1.0e-7)
    run_cpu!(H; use_nitrif_denitrif=false, use_lch4=false, use_mimics=false)
    D = run_dev(B; use_nitrif_denitrif=false, use_lch4=false, use_mimics=false)
    total += compare(H, D, "config-A")

    # --- Config B: immobilization branch (pmnf>0) + use_lch4 fphr ---
    println("  [B] pmnf>0 (fpi_vr-limited immobilization) + use_lch4 fphr")
    H = build(); B = build()
    fill!(H.pmnf_decomp_cascade, 1.0e-7); fill!(B.pmnf_decomp_cascade, 1.0e-7)
    run_cpu!(H; use_nitrif_denitrif=false, use_lch4=true, use_mimics=false)
    D = run_dev(B; use_nitrif_denitrif=false, use_lch4=true, use_mimics=false)
    total += compare(H, D, "config-B")

    # --- Config C: use_nitrif_denitrif=true (sminn_to_denit untouched) ---
    println("  [C] use_nitrif_denitrif=true")
    H = build(); B = build()
    fill!(H.pmnf_decomp_cascade, -1.0e-7); fill!(B.pmnf_decomp_cascade, -1.0e-7)
    run_cpu!(H; use_nitrif_denitrif=true, use_lch4=false, use_mimics=false)
    D = run_dev(B; use_nitrif_denitrif=true, use_lch4=false, use_mimics=false)
    total += compare(H, D, "config-C")

    # --- Config D: masked columns (skip 2 & 3) ---
    println("  [D] masked columns 2,3 skipped")
    H = build(); B = build()
    H.mask_bgc_soilc[2] = false; H.mask_bgc_soilc[3] = false
    B.mask_bgc_soilc[2] = false; B.mask_bgc_soilc[3] = false
    fill!(H.pmnf_decomp_cascade, -1.0e-7); fill!(B.pmnf_decomp_cascade, -1.0e-7)
    run_cpu!(H; use_nitrif_denitrif=false, use_lch4=true, use_mimics=false)
    D = run_dev(B; use_nitrif_denitrif=false, use_lch4=true, use_mimics=false)
    total += compare(H, D, "config-D")

    println()
    println(total == 0 ? "  WHOLE soil_biogeochem_decomp! MATCHES CPU ON $name ($FT) ✓" :
                         "  DIVERGENCE — investigate.")
    return total == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
