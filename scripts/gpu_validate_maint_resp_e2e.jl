# ==========================================================================
# gpu_validate_maint_resp_e2e.jl — end-to-end GPU parity for cn_mresp!
# (maintenance respiration: tcsoi column kernel + leaf/livewood patch kernel
# + froot-by-layer kernel). Validates the device-resident tcsoi scratch fix.
#
#   julia --project=scripts scripts/gpu_validate_maint_resp_e2e.jl
# ==========================================================================
using CLM
using Printf
include(joinpath(@__DIR__, "gpu_backends.jl"))

struct _F32 end
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:AbstractFloat}) = Float32.(x)
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:Integer}) = x
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{Bool}) = x
struct _F32S end   # TemperatureData pins ::FT scalars -> convert scalars too
CLM.Adapt.adapt_storage(::_F32S, x::AbstractArray{<:AbstractFloat}) = Float32.(x)
CLM.Adapt.adapt_storage(::_F32S, x::AbstractArray{<:Integer}) = x
CLM.Adapt.adapt_storage(::_F32S, x::AbstractArray{Bool}) = x
CLM.Adapt.adapt_storage(::_F32S, x::Float64) = Float32(x)

reldiff(a, b) = begin
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(Float64(A[i]) - Float64(B[i])) / (1.0 + max(abs(Float64(A[i])), abs(Float64(B[i])))))
    end
    m
end

function make_mresp(; np=2, nc=2)
    # Initialize varpar so the *_init! allocators get valid layer dims, then
    # derive nlevgrnd/nlevsno from it so every array + fill stays consistent.
    CLM.varpar_init!(CLM.varpar, 17, 17, 0, 5)
    nlevgrnd = CLM.varpar.nlevgrnd; nlevsno = CLM.varpar.nlevsno
    params = CLM.MaintRespParams(br=2.525e-6, br_root=2.525e-6)
    cn_params = CLM.CNSharedParamsData(Q10=1.5)
    pftcon = CLM.PftConMaintResp(woody=[0.0, 1.0])
    patch = CLM.PatchData(); CLM.patch_init!(patch, np)
    for p in 1:np; patch.itype[p] = 1; patch.column[p] = p; end
    cs = CLM.CanopyStateData(); CLM.canopystate_init!(cs, np)
    for p in 1:np
        cs.frac_veg_nosno_patch[p] = 1
        cs.laisun_patch[p] = 2.0 + p; cs.laisha_patch[p] = 1.0 + p
    end
    ss = CLM.SoilStateData(); CLM.soilstate_init!(ss, np, nc)
    for j in 1:nlevgrnd, p in 1:np; ss.crootfr_patch[p, j] = 1.0 / nlevgrnd; end
    temp = CLM.TemperatureData(); CLM.temperature_init!(temp, np, nc, 1, 1)
    for p in 1:np; temp.t_ref2m_patch[p] = 293.15; temp.t_a10_patch[p] = 293.15; end
    for j in 1:nlevgrnd, c in 1:nc; temp.t_soisno_col[c, j + nlevsno] = 290.0 + c; end
    ps = CLM.PhotosynthesisData(); CLM.photosynthesis_data_init!(ps, np)
    for p in 1:np; ps.lmrsun_patch[p] = 1.0 * p; ps.lmrsha_patch[p] = 0.5 * p; end
    ps.rootstem_acc = false
    cf = CLM.CNVegCarbonFluxData(); CLM.cnveg_carbon_flux_init!(cf, np, nc, 1)
    ns = CLM.CNVegNitrogenStateData(); CLM.cnveg_nitrogen_state_init!(ns, np, nc, 1)
    for p in 1:np
        ns.frootn_patch[p] = 0.01 * p; ns.livestemn_patch[p] = 0.05 * p
        ns.livecrootn_patch[p] = 0.03 * p
    end
    return (; params, cn_params, pftcon, patch, cs, ss, temp, ps, cf, ns,
            mask_c=BitVector(fill(true, nc)), mask_p=BitVector(fill(true, np)),
            bounds_c=1:nc, bounds_p=1:np)
end

run_mresp!(d) = CLM.cn_mresp!(d.mask_c, d.mask_p, d.bounds_c, d.bounds_p,
    d.params, d.cn_params, d.pftcon, d.patch, d.cs, d.ss, d.temp, d.ps, d.cf, d.ns)

function main(backend)
    println("="^66); println("END-TO-END GPU parity for cn_mresp! (maint_resp)"); println("="^66)
    if backend === nothing; println("  No GPU backend."); return 0; end
    name, _, FT = backend
    @printf("  Backend: %s   (precision: %s)\n", name, FT)
    H = make_mresp(); B = make_mresp()
    run_mresp!(H)
    mf(x)  = CLM.Adapt.adapt(device_array_type(), CLM.Adapt.adapt(_F32(), x))
    mfS(x) = CLM.Adapt.adapt(device_array_type(), CLM.Adapt.adapt(_F32S(), x))
    mm(b)  = device_array_type()(collect(b))
    D = (; params=B.params, cn_params=B.cn_params, pftcon=mf(B.pftcon), patch=mf(B.patch),
         cs=mf(B.cs), ss=mf(B.ss), temp=mfS(B.temp), ps=mf(B.ps), cf=mf(B.cf), ns=mf(B.ns),
         mask_c=mm(B.mask_c), mask_p=mm(B.mask_p), bounds_c=B.bounds_c, bounds_p=B.bounds_p)
    if !(D.cf.leaf_mr_patch isa device_array_type()); println("  BLOCKED."); return 2; end
    run_mresp!(D)
    checks = [("leaf_mr", H.cf.leaf_mr_patch, D.cf.leaf_mr_patch),
              ("froot_mr", H.cf.froot_mr_patch, D.cf.froot_mr_patch),
              ("livestem_mr", H.cf.livestem_mr_patch, D.cf.livestem_mr_patch),
              ("livecroot_mr", H.cf.livecroot_mr_patch, D.cf.livecroot_mr_patch),
              ("reproductive_mr", H.cf.reproductive_mr_patch, D.cf.reproductive_mr_patch)]
    nfail = 0
    for (nm, a, b) in checks
        dd = reldiff(a, b); ok = dd < 1f-3
        @printf("  [%s] %-18s rel = %.3e\n", ok ? "PASS" : "FAIL", nm, dd); ok || (nfail += 1)
    end
    println(); println(nfail == 0 ? "  cn_mresp! MATCHES CPU ON $name ($FT)" : "  DIVERGENCE ($nfail).")
    return nfail == 0 ? 0 : 1
end
const BACKEND = detect_backend()
exit(main(BACKEND))
