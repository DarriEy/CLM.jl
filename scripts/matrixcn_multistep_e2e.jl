# ==========================================================================
# matrixcn_multistep_e2e.jl — sustained multi-step numerical stability of the
# matrix-CN solve on Metal (Float32) vs the host (Float64). Single-step parity is
# ~1e-8 (gpu_validate_matrixcn_e2e); here we run N steps of the SAME solve with
# steady forcing and track the host-vs-device pool drift + a mass-conservation
# check, confirming the device path stays bounded/stable (Float32 error does not
# blow up) over a long run.
#
#   julia --project=scripts scripts/matrixcn_multistep_e2e.jl
# ==========================================================================
using CLM
using Printf
import Metal
using Metal: MtlArray
dev = MtlArray

struct DevF32{F}; f::F; end
CLM.Adapt.adapt_storage(d::DevF32, x::AbstractArray{<:AbstractFloat}) = d.f(Float32.(x))
CLM.Adapt.adapt_storage(d::DevF32, x::AbstractArray{<:Integer}) = d.f(collect(x))
CLM.Adapt.adapt_storage(d::DevF32, x::AbstractArray{Bool}) = d.f(collect(x))
mvf(x) = CLM.Adapt.adapt(DevF32(dev), x)
poolsum(cs) = sum(Float64.(Array(cs.leafc_patch))) + sum(Float64.(Array(cs.frootc_patch))) +
    sum(Float64.(Array(cs.livestemc_patch))) + sum(Float64.(Array(cs.deadstemc_patch))) +
    sum(Float64.(Array(cs.livecrootc_patch))) + sum(Float64.(Array(cs.deadcrootc_patch)))
leaf(cs) = Float64.(Array(cs.leafc_patch))

const NVEG = CLM.NVEGPOOL_NATVEG
const CNT = CLM.veg_matrix_transfer_counts(false)
const VEGF = (:leafc_patch, :leafc_storage_patch, :leafc_xfer_patch, :frootc_patch, :frootc_storage_patch,
    :frootc_xfer_patch, :livestemc_patch, :livestemc_storage_patch, :livestemc_xfer_patch, :deadstemc_patch,
    :deadstemc_storage_patch, :deadstemc_xfer_patch, :livecrootc_patch, :livecrootc_storage_patch,
    :livecrootc_xfer_patch, :deadcrootc_patch, :deadcrootc_storage_patch, :deadcrootc_xfer_patch)

function mk(np)
    cs = CLM.CNVegCarbonStateData(); CLM.cnveg_carbon_state_init!(cs, np, 1, 1; use_matrixcn=false, nrepr=1)
    for (k, f) in enumerate(VEGF); getfield(cs, f) .= Float64[20.0k + (p % 5) for p in 1:np]; end
    cf = CLM.CNVegCarbonFluxData(); CLM.cnveg_carbon_flux_init!(cf, np, 1, 1; use_matrixcn=true)
    CLM.cn_veg_matrix_c_topology!(cf; use_crop=false, nvegcpool=NVEG)
    # steady moderate forcing: some transfers + turnover + allocation input each step.
    fill!(cf.matrix_phtransfer_patch, 1.0e-6); fill!(cf.matrix_phturnover_patch, 0.02)
    fill!(cf.matrix_gmtransfer_patch, 5.0e-7); fill!(cf.matrix_gmturnover_patch, 0.01)
    fill!(cf.matrix_alloc_patch, 0.02); fill!(cf.matrix_Cinput_patch, 5.0e-4)
    return cs, cf
end

function main()
    np = 2000; N = 2000
    bk = Metal.functional() ? "Metal (Float32)" : "no Metal"
    println("=" ^ 66); println("matrix-CN multi-step stability: HOST (Float64) vs $bk"); @printf("  np=%d, N=%d steps\n", np, N); println("=" ^ 66)
    args = (; mask_soilp=trues(np), bounds_patch=1:np, ivt=ones(Int, np), woody=ones(np),
        npcropmin=15, nvegcpool=NVEG, counts=CNT, dt=1800.0, num_actfirep=0)

    # HOST trajectory
    st = CLM.CNVegMatrixSolveState(); csh, cfh = mk(np)
    CLM.cn_veg_matrix_solve_c!(csh, cfh; args..., state=st)      # step 1 (build structure)
    for _ in 2:N; CLM.cn_veg_matrix_solve_c!(csh, cfh; args..., state=st); end
    sumh = poolsum(csh); leafh = leaf(csh)

    # DEVICE trajectory (same initial state, memoized structure from the host warm-up)
    st_d = mvf(st); csd, _ = mk(np); csd = mvf(csd); cfd = mvf(cfh); ivt_d = dev(ones(Int, np))
    argsd = (; args..., ivt=ivt_d, ref=dev(zeros(Float32, 1, 1)), FT=Float32, state=st_d)
    st_d.init_ph = true; st_d.init_fi = true; st_d.ab_ready = true   # reuse the built structure
    for _ in 1:N; CLM.cn_veg_matrix_solve_c!(csd, cfd; argsd...); end
    Metal.synchronize()
    sumd = poolsum(csd); leafd = leaf(csd)

    reld = maximum(abs.(leafd .- leafh) ./ (1.0 .+ abs.(leafh)))
    @printf("  total veg-C after %d steps:  host %.6e   device %.6e   rel %.3e\n", N, sumh, sumd, abs(sumd - sumh) / sumh)
    @printf("  per-patch leafc drift (max rel host-vs-device): %.3e\n", reld)
    @printf("  all finite (device): %s\n", all(isfinite, leafd))
    ok = reld < 1.0e-3 && isfinite(sumd) && all(isfinite, leafd)
    println(); println(ok ? "  STABLE — device Float32 tracks host Float64 over $N steps (no blow-up)" : "  UNSTABLE — drift/NaN exceeded tolerance")
    return ok ? 0 : 1
end
exit(main())
