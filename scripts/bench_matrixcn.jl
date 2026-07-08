# ==========================================================================
# bench_matrixcn.jl — wall-time of the matrix-CN solve on the HOST (KA CPU
# backend, Float64) vs METAL (Float32), at realistic column/patch counts. The
# solve is per-patch/column parallel, so the GPU payoff grows with the unit
# count. Times the STEADY-STATE per-step cost: structure is prebuilt once (host
# warm-up), then the memoized solve is timed (that is what runs every timestep).
#
#   julia --project=scripts scripts/bench_matrixcn.jl
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

# min elapsed over K reps (after a warm-up); Metal solve is synchronous (each
# _launch! calls KA.synchronize), and we add a final device sync to be safe.
function timeit(f, K; sync=false)
    f(); sync && Metal.synchronize()
    t = Inf
    for _ in 1:K
        dt = @elapsed (f(); sync && Metal.synchronize())
        t = min(t, dt)
    end
    return t
end

const NVEG = CLM.NVEGPOOL_NATVEG
const CNT = CLM.veg_matrix_transfer_counts(false)
const VEGF = (:leafc_patch, :leafc_storage_patch, :leafc_xfer_patch, :frootc_patch, :frootc_storage_patch,
    :frootc_xfer_patch, :livestemc_patch, :livestemc_storage_patch, :livestemc_xfer_patch, :deadstemc_patch,
    :deadstemc_storage_patch, :deadstemc_xfer_patch, :livecrootc_patch, :livecrootc_storage_patch,
    :livecrootc_xfer_patch, :deadcrootc_patch, :deadcrootc_storage_patch, :deadcrootc_xfer_patch)

function mk_vegc(np)
    cs = CLM.CNVegCarbonStateData(); CLM.cnveg_carbon_state_init!(cs, np, 1, 1; use_matrixcn=false, nrepr=1)
    for (k, f) in enumerate(VEGF); getfield(cs, f) .= Float64[10.0k + (p % 7) for p in 1:np]; end
    cf = CLM.CNVegCarbonFluxData(); CLM.cnveg_carbon_flux_init!(cf, np, 1, 1; use_matrixcn=true)
    CLM.cn_veg_matrix_c_topology!(cf; use_crop=false, nvegcpool=NVEG)
    fill!(cf.matrix_phtransfer_patch, 1.0e-6); fill!(cf.matrix_phturnover_patch, 0.05)
    fill!(cf.matrix_gmtransfer_patch, 5.0e-7); fill!(cf.matrix_gmturnover_patch, 0.02)
    fill!(cf.matrix_alloc_patch, 0.01); fill!(cf.matrix_Cinput_patch, 1.0e-3)
    return cs, cf
end

function bench_veg(np, K)
    ivt = ones(Int, np); woody = ones(np); args = (; mask_soilp=trues(np), bounds_patch=1:np,
        ivt=ivt, woody=woody, npcropmin=15, nvegcpool=NVEG, counts=CNT, dt=1800.0, num_actfirep=0)
    st = CLM.CNVegMatrixSolveState()
    cs, cf = mk_vegc(np); CLM.cn_veg_matrix_solve_c!(cs, cf; args..., state=st)   # host warm-up (build)
    th = timeit(K) do
        CLM.cn_veg_matrix_solve_c!(cs, cf; args..., state=st)
    end
    cs_d = mvf(mk_vegc(np)[1]); cf_d = mvf(cf); st_d = mvf(st); ivt_d = dev(ivt); ref = dev(zeros(Float32, 1, 1))
    argsd = (; args..., ivt=ivt_d, ref=ref, FT=Float32, state=st_d)
    td = timeit(K; sync=true) do
        CLM.cn_veg_matrix_solve_c!(cs_d, cf_d; argsd...)
    end
    return th, td
end

# ---- soil solve (per-column) ----
function mk_soil(nc)
    nlev = 2; npool = 3; ndct = 3; nout = 1; cn = [20.0, 15.0, 90.0]
    cc = CLM.DecompCascadeConData(cascade_donor_pool=[1, 2, 3], cascade_receiver_pool=[2, 0, 1],
        floating_cn_ratio_decomp_pools=BitVector([false, false, false]),
        is_cwd=BitVector([false, false, true]), initial_cn_ratio=cn)
    CLM.init_soil_transfer!(cc; ndecomp_pools=npool, nlevdecomp=nlev, ndecomp_cascade_transitions=ndct, ndecomp_cascade_outtransitions=nout)
    ndpvr = npool * nlev
    rf = zeros(nc, nlev, ndct); rf[:, :, 1] .= 0.3; rf[:, :, 2] .= 1.0
    pf = ones(nc, nlev, ndct)
    Ks = zeros(nc, ndpvr); for i in 1:npool, j in 1:nlev; Ks[:, j + (i - 1) * nlev] .= 0.02i; end
    Cin = zeros(nc, ndpvr); Nin = zeros(nc, ndpvr)
    for i in 1:npool, j in 1:nlev; Cin[:, j + (i - 1) * nlev] .= 2.0i; Nin[:, j + (i - 1) * nlev] .= 2.0i / cn[i]; end
    mkC() = (C = zeros(nc, nlev, npool); for i in 1:npool, j in 1:nlev; C[:, j, i] .= 100.0i + 10j; end; C)
    mkN() = (N = zeros(nc, nlev, npool); for i in 1:npool, j in 1:nlev; N[:, j, i] .= (100.0i + 10j) / cn[i]; end; N)
    return (; cc, nlev, npool, ndct, nout, ndpvr, rf, pf, Ks, Cin, Nin, mkC, mkN)
end

function bench_soil(nc, K)
    d = mk_soil(nc)
    hargs() = (; Ksoil=copy(d.Ks), tri_ma_vr=zeros(nc, d.cc.Ntri_setup), matrix_Cinput=copy(d.Cin),
        matrix_Ninput=copy(d.Nin), rf_decomp_cascade=d.rf, pathfrac_decomp_cascade=d.pf, mask_soilc=trues(nc),
        begc=1, endc=nc, nlevdecomp=d.nlev, ndecomp_pools=d.npool, ndecomp_cascade_transitions=d.ndct, ndecomp_cascade_outtransitions=d.nout)
    CLM.cn_soil_matrix!(CLM.CNSoilMatrixState(), d.cc; decomp_cpools_vr=d.mkC(), decomp_npools_vr=d.mkN(), hargs()...)  # warm-up (build cc)
    msh = CLM.CNSoilMatrixState(); Ch = d.mkC(); Nh = d.mkN()
    CLM.cn_soil_matrix!(msh, d.cc; decomp_cpools_vr=Ch, decomp_npools_vr=Nh, hargs()...)  # host warm-up (ms flags)
    th = timeit(K) do; CLM.cn_soil_matrix!(msh, d.cc; decomp_cpools_vr=Ch, decomp_npools_vr=Nh, hargs()...); end
    msd = CLM.CNSoilMatrixState(); msd.init_readyAsoilc = true; msd.init_readyAsoiln = true
    msd.list_ready1_nofire = true; msd.list_ready2_nofire = true
    Cd = mvf(d.mkC()); Nd = mvf(d.mkN()); ref = dev(zeros(Float32, 1, 1))
    dargs() = (; Ksoil=mvf(d.Ks), tri_ma_vr=dev(zeros(Float32, nc, d.cc.Ntri_setup)), matrix_Cinput=mvf(d.Cin),
        matrix_Ninput=mvf(d.Nin), rf_decomp_cascade=mvf(d.rf), pathfrac_decomp_cascade=mvf(d.pf), mask_soilc=dev(trues(nc)),
        begc=1, endc=nc, nlevdecomp=d.nlev, ndecomp_pools=d.npool, ndecomp_cascade_transitions=d.ndct, ndecomp_cascade_outtransitions=d.nout, ref=ref, FT=Float32)
    td = timeit(K; sync=true) do; CLM.cn_soil_matrix!(msd, d.cc; decomp_cpools_vr=Cd, decomp_npools_vr=Nd, dargs()...); end
    return th, td
end

function main()
    bk = Metal.functional() ? "Metal (Float32)" : "no Metal"
    println("=" ^ 62); println("matrix-CN solve: HOST (KA CPU, Float64) vs $bk"); println("=" ^ 62)
    println("  VEG-C solve (per patch):")
    @printf("  %-10s %12s %12s %9s %14s\n", "np", "host (ms)", "metal (ms)", "speedup", "M patch/s")
    for np in (1_000, 10_000, 100_000, 500_000)
        th, td = bench_veg(np, np <= 10_000 ? 50 : 10)
        @printf("  %-10d %12.3f %12.3f %8.2fx %14.1f\n", np, th * 1e3, td * 1e3, th / td, np / td / 1e6)
    end
    println("\n  SOIL-C/N solve (per column):")
    @printf("  %-10s %12s %12s %9s %14s\n", "nc", "host (ms)", "metal (ms)", "speedup", "M col/s")
    for nc in (1_000, 10_000, 100_000, 500_000)
        th, td = bench_soil(nc, nc <= 10_000 ? 50 : 10)
        @printf("  %-10d %12.3f %12.3f %8.2fx %14.1f\n", nc, th * 1e3, td * 1e3, th / td, nc / td / 1e6)
    end
    println("\n  host = KernelAbstractions CPU backend, JULIA_NUM_THREADS=$(Threads.nthreads()).")
    println("  Payoff grows with the unit count (real global runs = 10^5–10^6 patches/cols).")
    println("  Optimizations applied: sync-fusion (1 GPU sync/solve, not ~20) + a cached")
    println("  device filter on the solve state. Residual small-size overhead is now the")
    println("  per-op kernel enqueue/command-buffer latency (~20 launches) — would need")
    println("  kernel FUSION (merging the glue kernels) to cut further.")
end
main()
