# ==========================================================================
# gpu_validate_decomppotential_e2e.jl — end-to-end GPU parity for the WHOLE
# soil_bgc_potential! BGC driver (potential decomposition / immobilization demand).
#
# Builds two small multi-column / multi-level instances mirroring
# test/test_decomp_potential.jl — one standard (use_mimics=false) and one MIMICS
# (use_mimics=true, microbial cop/oli pools) — runs soil_bgc_potential! on the CPU,
# converts every state struct + loose output array to Float32 + adapts to Metal,
# runs the SAME call on the device, and compares the mutated outputs field-by-field.
# The MIMICS scenario exercises the per-column in-thread receiver-pool gain cascade
# (_decpot_mimics_kernel!) + the immob/gross_nmin DIN accumulation; the standard
# scenario exercises the fixed/floating C:N, atmosphere-respiration, and
# immob/gross_nmin/phr accumulation paths.
#
#   julia --project=scripts scripts/gpu_validate_decomppotential_e2e.jl
#
# The CN/soil-BGC structs are concretely Float64-typed (default ctor pins
# ::Array{Float64}); Metal has no Float64, so we build at Float64 and down-convert
# array fields to Float32 for the device snapshot via the MetalF32 adaptor.
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

# Float32-down-converting Metal adaptor: down-converts float arrays to Float32 as it
# reconstructs each struct (integer/bool arrays move as-is). @adapt_structure rebuilds
# each struct positionally, inferring the new {Float32,…} params from adapted fields.
struct MetalF32 end
CLM.Adapt.adapt_storage(::MetalF32, x::AbstractArray{<:AbstractFloat}) = device_array_type()(Float32.(x))
CLM.Adapt.adapt_storage(::MetalF32, x::AbstractArray{<:Integer})       = device_array_type()(x)
CLM.Adapt.adapt_storage(::MetalF32, x::AbstractArray{Bool})            = device_array_type()(x)

# --------------------------------------------------------------------------
# Standard (non-MIMICS) instance — mirrors test_decomp_potential.jl make_test_data
# --------------------------------------------------------------------------
function build_std(; nc=4, nlevdecomp=3, ndecomp_pools=4, ntrans=5)
    cascade_con = CLM.DecompCascadeConData()
    cascade_con.cascade_donor_pool    = [1, 2, 3, 3, 4]
    cascade_con.cascade_receiver_pool = [3, 3, 4, CLM.I_ATM, 3]
    cascade_con.floating_cn_ratio_decomp_pools = BitVector([true, true, false, false])
    cascade_con.initial_cn_ratio = [20.0, 25.0, 12.0, 10.0]

    cs = CLM.SoilBiogeochemCarbonStateData()
    cs.decomp_cpools_vr_col = zeros(nc, nlevdecomp, ndecomp_pools)
    for c in 1:nc, j in 1:nlevdecomp
        cs.decomp_cpools_vr_col[c, j, 1] = 100.0 + 10.0 * c
        cs.decomp_cpools_vr_col[c, j, 2] = 80.0 + 5.0 * c
        cs.decomp_cpools_vr_col[c, j, 3] = 50.0 + 3.0 * c
        cs.decomp_cpools_vr_col[c, j, 4] = 30.0 + 2.0 * c
    end

    ns = CLM.SoilBiogeochemNitrogenStateData()
    ns.decomp_npools_vr_col = zeros(nc, nlevdecomp, ndecomp_pools)
    for c in 1:nc, j in 1:nlevdecomp
        ns.decomp_npools_vr_col[c, j, 1] = cs.decomp_cpools_vr_col[c, j, 1] / 20.0
        ns.decomp_npools_vr_col[c, j, 2] = cs.decomp_cpools_vr_col[c, j, 2] / 16.0
        ns.decomp_npools_vr_col[c, j, 3] = cs.decomp_cpools_vr_col[c, j, 3] / 12.0
        ns.decomp_npools_vr_col[c, j, 4] = cs.decomp_cpools_vr_col[c, j, 4] / 10.0
    end

    cf = CLM.SoilBiogeochemCarbonFluxData()
    cf.rf_decomp_cascade_col       = zeros(nc, nlevdecomp, ntrans)
    cf.pathfrac_decomp_cascade_col = ones(nc, nlevdecomp, ntrans)
    cf.decomp_k_col                = zeros(nc, nlevdecomp, ndecomp_pools)
    cf.cn_col                      = zeros(nc, ndecomp_pools)
    cf.phr_vr_col                  = zeros(nc, nlevdecomp)
    for c in 1:nc, j in 1:nlevdecomp
        cf.rf_decomp_cascade_col[c, j, 1] = 0.39
        cf.rf_decomp_cascade_col[c, j, 2] = 0.55
        cf.rf_decomp_cascade_col[c, j, 3] = 0.28
        cf.rf_decomp_cascade_col[c, j, 4] = 0.55
        cf.rf_decomp_cascade_col[c, j, 5] = 0.55
        cf.decomp_k_col[c, j, 1] = 0.01
        cf.decomp_k_col[c, j, 2] = 0.005
        cf.decomp_k_col[c, j, 3] = 0.002
        cf.decomp_k_col[c, j, 4] = 0.001
    end

    st = CLM.SoilBiogeochemStateData()
    st.nue_decomp_cascade_col = ones(ntrans) * 0.5

    nf = CLM.SoilBiogeochemNitrogenFluxData()
    nf.potential_immob_vr_col = zeros(nc, nlevdecomp)
    nf.gross_nmin_vr_col      = zeros(nc, nlevdecomp)

    return (; cf, cs, nf, ns, st, cascade_con,
              mask = trues(nc), nc, nlevdecomp, ndecomp_pools, ntrans,
              use_mimics=false, i_cop_mic=0, i_oli_mic=0)
end

# --------------------------------------------------------------------------
# MIMICS instance — mirrors test_decomp_potential.jl "MIMICS pathway" testset
# --------------------------------------------------------------------------
function build_mimics(; nc=2, nlevdecomp=2)
    ndecomp_pools = 6
    ntrans = 3
    i_cop = 4; i_oli = 5

    cascade_con = CLM.DecompCascadeConData()
    cascade_con.cascade_donor_pool    = [1, 2, 3]
    cascade_con.cascade_receiver_pool = [i_cop, i_oli, i_cop]
    cascade_con.floating_cn_ratio_decomp_pools = BitVector([true, true, false, true, true, true])
    cascade_con.initial_cn_ratio = [20.0, 25.0, 12.0, 10.0, 10.0, 500.0]

    cs = CLM.SoilBiogeochemCarbonStateData()
    cs.decomp_cpools_vr_col = zeros(nc, nlevdecomp, ndecomp_pools)
    cs.decomp_cpools_vr_col[:, :, 1] .= 100.0
    cs.decomp_cpools_vr_col[:, :, 2] .= 80.0
    cs.decomp_cpools_vr_col[:, :, 3] .= 50.0
    cs.decomp_cpools_vr_col[:, :, 4] .= 20.0
    cs.decomp_cpools_vr_col[:, :, 5] .= 15.0

    ns = CLM.SoilBiogeochemNitrogenStateData()
    ns.decomp_npools_vr_col = zeros(nc, nlevdecomp, ndecomp_pools)
    ns.decomp_npools_vr_col[:, :, 1] .= 5.0
    ns.decomp_npools_vr_col[:, :, 2] .= 3.2
    ns.decomp_npools_vr_col[:, :, 3] .= 4.167
    ns.decomp_npools_vr_col[:, :, 4] .= 2.0
    ns.decomp_npools_vr_col[:, :, 5] .= 1.5

    cf = CLM.SoilBiogeochemCarbonFluxData()
    cf.rf_decomp_cascade_col = zeros(nc, nlevdecomp, ntrans)
    cf.rf_decomp_cascade_col[:, :, 1] .= 0.3
    cf.rf_decomp_cascade_col[:, :, 2] .= 0.4
    cf.rf_decomp_cascade_col[:, :, 3] .= 0.2
    cf.pathfrac_decomp_cascade_col = ones(nc, nlevdecomp, ntrans)
    cf.decomp_k_col = zeros(nc, nlevdecomp, ndecomp_pools)
    cf.decomp_k_col[:, :, 1] .= 0.01
    cf.decomp_k_col[:, :, 2] .= 0.005
    cf.decomp_k_col[:, :, 3] .= 0.002
    cf.cn_col = zeros(nc, ndecomp_pools)
    cf.cn_col[:, i_cop] .= 8.0
    cf.cn_col[:, i_oli] .= 8.0
    cf.phr_vr_col = zeros(nc, nlevdecomp)

    st = CLM.SoilBiogeochemStateData()
    st.nue_decomp_cascade_col = [0.5, 0.6, 0.5]

    nf = CLM.SoilBiogeochemNitrogenFluxData()
    nf.potential_immob_vr_col = zeros(nc, nlevdecomp)
    nf.gross_nmin_vr_col      = zeros(nc, nlevdecomp)

    return (; cf, cs, nf, ns, st, cascade_con,
              mask = trues(nc), nc, nlevdecomp, ndecomp_pools, ntrans,
              use_mimics=true, i_cop_mic=i_cop, i_oli_mic=i_oli)
end

# Output arrays soil_bgc_potential! mutates (passed loose as kwargs).
make_outs(D, FT) = (
    cn_decomp_pools       = zeros(FT, D.nc, D.nlevdecomp, D.ndecomp_pools),
    p_decomp_cpool_loss   = zeros(FT, D.nc, D.nlevdecomp, D.ntrans),
    p_decomp_cn_gain      = zeros(FT, D.nc, D.nlevdecomp, D.ndecomp_pools),
    pmnf_decomp_cascade   = zeros(FT, D.nc, D.nlevdecomp, D.ntrans),
    p_decomp_npool_to_din = zeros(FT, D.nc, D.nlevdecomp, D.ntrans),
)

run_potential!(D, mask, O) = CLM.soil_bgc_potential!(
    D.cf, D.cs, D.nf, D.ns, D.st, D.cascade_con;
    mask_bgc_soilc=mask, bounds=1:D.nc, nlevdecomp=D.nlevdecomp,
    ndecomp_pools=D.ndecomp_pools, ndecomp_cascade_transitions=D.ntrans,
    cn_decomp_pools=O.cn_decomp_pools, p_decomp_cpool_loss=O.p_decomp_cpool_loss,
    p_decomp_cn_gain=O.p_decomp_cn_gain, pmnf_decomp_cascade=O.pmnf_decomp_cascade,
    p_decomp_npool_to_din=O.p_decomp_npool_to_din,
    use_mimics=D.use_mimics, i_cop_mic=D.i_cop_mic, i_oli_mic=D.i_oli_mic)

function run_case(name, builder, backend)
    name_b, to, FT = backend
    println("-" ^ 70)
    @printf("  Scenario: %s   (backend %s, precision %s)\n", name, name_b, FT)

    H  = builder()                       # CPU reference (Float64)
    B  = builder()                       # source for the device snapshot
    Ho = make_outs(H, Float64)
    Bo = make_outs(B, Float64)

    ad(x) = CLM.Adapt.adapt(MetalF32(), x)
    cf_d = ad(B.cf); cs_d = ad(B.cs); nf_d = ad(B.nf); ns_d = ad(B.ns)
    st_d = ad(B.st); cc_d = ad(B.cascade_con)
    Dd = (; cf=cf_d, cs=cs_d, nf=nf_d, ns=ns_d, st=st_d, cascade_con=cc_d,
            nc=B.nc, nlevdecomp=B.nlevdecomp, ndecomp_pools=B.ndecomp_pools,
            ntrans=B.ntrans, use_mimics=B.use_mimics,
            i_cop_mic=B.i_cop_mic, i_oli_mic=B.i_oli_mic)

    if !(cs_d.decomp_cpools_vr_col isa device_array_type())
        println("  BLOCKED: a state struct did not move to the device under adapt.")
        return 2
    end

    # device-resident output arrays (Float32)
    Od = (
        cn_decomp_pools       = to(Float32.(Bo.cn_decomp_pools)),
        p_decomp_cpool_loss   = to(Float32.(Bo.p_decomp_cpool_loss)),
        p_decomp_cn_gain      = to(Float32.(Bo.p_decomp_cn_gain)),
        pmnf_decomp_cascade   = to(Float32.(Bo.pmnf_decomp_cascade)),
        p_decomp_npool_to_din = to(Float32.(Bo.p_decomp_npool_to_din)),
    )
    dmask(m) = to(collect(Bool, m))

    run_potential!(H,  H.mask,        Ho)   # CPU
    run_potential!(Dd, dmask(B.mask), Od)   # device

    checks = [
        ("cn_decomp_pools",     Ho.cn_decomp_pools,       Od.cn_decomp_pools),
        ("p_decomp_cpool_loss", Ho.p_decomp_cpool_loss,   Od.p_decomp_cpool_loss),
        ("p_decomp_cn_gain",    Ho.p_decomp_cn_gain,      Od.p_decomp_cn_gain),
        ("pmnf_decomp_cascade", Ho.pmnf_decomp_cascade,   Od.pmnf_decomp_cascade),
        ("p_decomp_npool_to_din", Ho.p_decomp_npool_to_din, Od.p_decomp_npool_to_din),
        ("potential_immob_vr",  H.nf.potential_immob_vr_col, Dd.nf.potential_immob_vr_col),
        ("gross_nmin_vr",       H.nf.gross_nmin_vr_col,      Dd.nf.gross_nmin_vr_col),
        ("phr_vr",              H.cf.phr_vr_col,             Dd.cf.phr_vr_col),
    ]
    nfail = 0
    for (nm, a, b) in checks
        if !cpu_has_finite(a)
            @printf("    [WARN] %-22s CPU reference all-NaN/Inf — skipping\n", nm)
            continue
        end
        d = reldiff(a, b); ok = d < 1f-3
        @printf("    [%s] %-22s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
        ok || (nfail += 1)
    end
    return nfail
end

function main(backend)
    println("=" ^ 70)
    println("END-TO-END GPU parity for soil_bgc_potential! (potential decomposition)")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU path exercised by the suite).")
        return 0
    end
    nfail = 0
    nfail += run_case("standard (use_mimics=false)", build_std,    backend)
    nfail += run_case("MIMICS  (use_mimics=true)",  build_mimics, backend)
    println()
    println(nfail == 0 ? "  WHOLE soil_bgc_potential! MATCHES CPU ON DEVICE ✓" :
                         "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
