# ==========================================================================
# gpu_validate_littervert_e2e.jl — end-to-end GPU parity for the WHOLE
# litter_vert_transp! routine (the layer-coupled advection-diffusion vertical
# transport of all decomposing C and N pools).
#
# The only layer-coupled step is the per-column tridiagonal (Thomas) solve, which
# a single column thread runs sequentially on its own row — so the whole routine is
# two per-column kernels (a one-shot coefficient kernel + a transport kernel launched
# per (tracer-type, pool) with the in-thread Thomas solve, inlined to match
# tridiagonal_solve! exactly). This harness clones test/test_litter_vert_transp.jl's
# setup at Float64, runs the function on the CPU, adapts every state struct to Float32
# device arrays (MetalF32), runs the SAME call on the device, and compares the
# mutated outputs with a NaN-aware reldiff.
#
#   julia --project=scripts scripts/gpu_validate_littervert_e2e.jl
# ==========================================================================

using CLM
using Printf
include(joinpath(@__DIR__, "gpu_backends.jl"))

function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(Float64(A[i]) - Float64(B[i])) / (1.0 + max(abs(Float64(A[i])), abs(Float64(B[i])))))
    end
    return m
end
cpu_has_finite(a) = any(isfinite, Array(a))

# Float32-down-converting Metal adaptor (CN state structs build Float64).
struct MetalF32 end
CLM.Adapt.adapt_storage(::MetalF32, x::AbstractArray{<:AbstractFloat}) = device_array_type()(Float32.(x))
CLM.Adapt.adapt_storage(::MetalF32, x::AbstractArray{<:Integer})       = device_array_type()(x)
CLM.Adapt.adapt_storage(::MetalF32, x::AbstractArray{Bool})            = device_array_type()(x)

# Build a Float64 instance mirroring test/test_litter_vert_transp.jl make_test_data.
# `altmax_val` selects the mixing regime: <= max_altdepth_cryoturbation -> cryoturbation,
# else bioturbation. spinup_state>=1 exercises the spinup-factor/latitude term.
function build(; nc=3, nlevdecomp=5, ndecomp_pools=7, cwd_pool=5, altmax_val=1.0)
    scalez = 0.025
    nlevgrnd = nlevdecomp + 5
    zsoi_vals = [scalez * (exp(0.5 * (j - 0.5)) - 1.0) for j in 1:nlevgrnd]
    zisoi_vals = zeros(nlevgrnd + 1)
    for j in 1:nlevgrnd
        zisoi_vals[j + 1] = (j < nlevgrnd) ? 0.5 * (zsoi_vals[j] + zsoi_vals[j + 1]) :
                                             zsoi_vals[j] + 0.5 * (zsoi_vals[j] - zisoi_vals[j])
    end
    dzsoi_decomp_vals = [zisoi_vals[j + 1] - zisoi_vals[j] for j in 1:nlevdecomp]

    cs = CLM.SoilBiogeochemCarbonStateData()
    cs.decomp_cpools_vr_col = zeros(nc, nlevdecomp, ndecomp_pools)
    for s in 1:ndecomp_pools, j in 1:nlevdecomp, c in 1:nc
        cs.decomp_cpools_vr_col[c, j, s] = 100.0 * exp(-0.5 * j) * (1.0 + 0.1 * s)
    end

    cf = CLM.SoilBiogeochemCarbonFluxData()
    cf.decomp_cpools_sourcesink_col = fill(1.0e-3, nc, nlevdecomp, ndecomp_pools)
    cf.decomp_cpools_transport_tendency_col = zeros(nc, nlevdecomp, ndecomp_pools)

    ns = CLM.SoilBiogeochemNitrogenStateData()
    ns.decomp_npools_vr_col = cs.decomp_cpools_vr_col ./ 15.0

    nf = CLM.SoilBiogeochemNitrogenFluxData()
    nf.decomp_npools_sourcesink_col = fill(1.0e-4, nc, nlevdecomp, ndecomp_pools)
    nf.decomp_npools_transport_tendency_col = zeros(nc, nlevdecomp, ndecomp_pools)

    st = CLM.SoilBiogeochemStateData()
    st.som_adv_coef_col = zeros(nc, nlevdecomp + 1)
    st.som_diffus_coef_col = zeros(nc, nlevdecomp + 1)

    col = CLM.ColumnData()
    col.nbedrock = fill(nlevdecomp, nc)
    col.gridcell = ones(Int, nc)

    grc = CLM.GridcellData()
    grc.latdeg = fill(45.0, 1)

    active_layer = CLM.ActiveLayerData()
    active_layer.altmax_col = fill(altmax_val, nc)
    active_layer.altmax_lastyear_col = fill(altmax_val, nc)

    cascade_con = CLM.DecompCascadeConData()
    is_cwd = falses(ndecomp_pools); is_cwd[cwd_pool] = true
    cascade_con.is_cwd = is_cwd
    cascade_con.spinup_factor = fill(2.0, ndecomp_pools)  # nontrivial spinup factor

    params = CLM.LitterVertTranspParams()
    params.som_diffus = 1.0e-4
    params.cryoturb_diffusion_k = 5.0e-4
    params.max_altdepth_cryoturbation = 2.0

    return (; cs, cf, ns, nf, st, col, grc, cascade_con, active_layer, params,
              mask_bgc_soilc=trues(nc), bounds=1:nc, dtime=1800.0,
              nlevdecomp, ndecomp_pools, zsoi_vals, dzsoi_decomp_vals, zisoi_vals)
end

run_lvt!(S; spinup_state) = CLM.litter_vert_transp!(
    S.cs, S.cf, S.ns, S.nf, S.st, S.col, S.grc, S.cascade_con, S.active_layer, S.params;
    mask_bgc_soilc=S.mask_bgc_soilc, bounds=S.bounds, dtime=S.dtime,
    nlevdecomp=S.nlevdecomp, ndecomp_pools=S.ndecomp_pools,
    zsoi_vals=S.zsoi_vals, dzsoi_decomp_vals=S.dzsoi_decomp_vals,
    zisoi_vals=S.zisoi_vals, spinup_state=spinup_state)

function check_case(to, FT, label; altmax_val, spinup_state)
    @printf("\n  --- %s (altmax=%.1f, spinup_state=%d) ---\n", label, altmax_val, spinup_state)
    H = build(; altmax_val=altmax_val)   # CPU reference (Float64)
    B = build(; altmax_val=altmax_val)   # device source snapshot

    ad(x) = CLM.Adapt.adapt(MetalF32(), x)
    Sd = (; cs=ad(B.cs), cf=ad(B.cf), ns=ad(B.ns), nf=ad(B.nf), st=ad(B.st),
            col=ad(B.col), grc=ad(B.grc), active_layer=ad(B.active_layer),
            B.cascade_con, B.params, mask_bgc_soilc=B.mask_bgc_soilc, B.bounds,
            dtime=FT(B.dtime), B.nlevdecomp, B.ndecomp_pools,
            zsoi_vals=B.zsoi_vals, dzsoi_decomp_vals=B.dzsoi_decomp_vals, zisoi_vals=B.zisoi_vals)

    if !(Sd.cs.decomp_cpools_vr_col isa device_array_type())
        println("  BLOCKED: carbon state did not move to the device under adapt.")
        return 2
    end

    run_lvt!(H; spinup_state=spinup_state)
    run_lvt!(Sd; spinup_state=spinup_state)

    checks = [
        ("decomp_cpools_vr", H.cs.decomp_cpools_vr_col,                 Sd.cs.decomp_cpools_vr_col),
        ("decomp_npools_vr", H.ns.decomp_npools_vr_col,                 Sd.ns.decomp_npools_vr_col),
        ("c_transp_tend",    H.cf.decomp_cpools_transport_tendency_col, Sd.cf.decomp_cpools_transport_tendency_col),
        ("n_transp_tend",    H.nf.decomp_npools_transport_tendency_col, Sd.nf.decomp_npools_transport_tendency_col),
        ("som_adv_coef",     H.st.som_adv_coef_col,                     Sd.st.som_adv_coef_col),
        ("som_diffus_coef",  H.st.som_diffus_coef_col,                  Sd.st.som_diffus_coef_col),
    ]
    nfail = 0
    for (nm, a, b) in checks
        if !cpu_has_finite(a)
            @printf("  [WARN] %-16s CPU reference all-NaN/Inf — skipping\n", nm)
            continue
        end
        d = reldiff(a, b); ok = d < 1f-3
        @printf("  [%s] %-16s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
        ok || (nfail += 1)
    end
    return nfail
end

function main(backend)
    println("=" ^ 70)
    println("END-TO-END GPU parity for litter_vert_transp! (whole vertical transport)")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU driver exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n", name, FT)

    nfail = 0
    nfail += check_case(to, FT, "cryoturbation";  altmax_val=1.0, spinup_state=0)
    nfail += check_case(to, FT, "bioturbation";   altmax_val=3.0, spinup_state=0)
    nfail += check_case(to, FT, "spinup+latterm"; altmax_val=1.0, spinup_state=1)

    println()
    println(nfail == 0 ? "  WHOLE litter_vert_transp! MATCHES CPU ON $name ($FT)" :
                         "  DIVERGENCE — investigate ($nfail field(s) failed).")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
