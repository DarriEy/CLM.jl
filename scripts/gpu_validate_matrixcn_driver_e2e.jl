# ==========================================================================
# gpu_validate_matrixcn_driver_e2e.jl — WHOLE-DRIVER GPU parity for the
# matrix-CN path. Runs cn_driver_no_leaching! with use_matrixcn=true (the veg-C
# + veg-N matrix solves, fed by the kernelized wiring accumulate helpers) on the
# CPU and on Metal over the SAME synthetic CN state, and compares the veg pools.
# The topology-static sparse structure is built once on the host (a warm-up call
# fills the CNVegMatrixSolveStates); the device run reuses it via the kernelized
# memoized fills. This is the capstone: the driver reaches the matrix branch and
# every piece (wiring accumulate -> solve) executes on the GPU.
#
#   julia --project=scripts scripts/gpu_validate_matrixcn_driver_e2e.jl
# ==========================================================================
using CLM
using Printf
include(joinpath(@__DIR__, "gpu_backends.jl"))

struct _F32 end
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:AbstractFloat}) = Float32.(x)
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:Integer}) = x
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{Bool}) = x

function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(Float64(A[i]) - Float64(B[i])) / (1.0 + max(abs(Float64(A[i])), abs(Float64(B[i])))))
    end
    return m
end

function make_data(; nc=4, np=6, ng=2, nlevdecomp=1, ndecomp_pools=7, ndecomp_cascade_transitions=5, nrepr=1)
    i_litr_min=1; i_litr_max=3; i_cwd=4; dt=1800.0
    mask_soilc = trues(nc); mask_soilp = trues(np)
    cs_veg = CLM.CNVegCarbonStateData(); CLM.cnveg_carbon_state_init!(cs_veg, np, nc, ng; nrepr=nrepr)
    CLM.cnveg_carbon_state_set_values!(cs_veg, mask_soilp, 10.0, mask_soilc, 10.0; nrepr=nrepr)
    cs_veg.cpool_patch .= 100.0; cs_veg.leafc_patch .= 50.0; cs_veg.frootc_patch .= 40.0
    cs_veg.livestemc_patch .= 30.0; cs_veg.deadstemc_patch .= 200.0; cs_veg.xsmrpool_patch .= 15.0
    cf_veg = CLM.CNVegCarbonFluxData()
    CLM.cnveg_carbon_flux_init!(cf_veg, np, nc, ng; nrepr=nrepr, nlevdecomp_full=nlevdecomp, ndecomp_pools=ndecomp_pools)
    CLM.cnveg_carbon_flux_set_values!(cf_veg, mask_soilp, 0.0, mask_soilc, 0.0; nrepr=nrepr, nlevdecomp_full=nlevdecomp, ndecomp_pools=ndecomp_pools)
    cf_veg.cpool_to_leafc_patch .= 1.0e-3; cf_veg.cpool_to_frootc_patch .= 5.0e-4; cf_veg.cpool_to_livestemc_patch .= 3.0e-4
    cf_veg.leafc_to_litter_patch .= 2.0e-4; cf_veg.frootc_to_litter_patch .= 1.5e-4
    cf_veg.leafc_storage_to_xfer_patch .= 1.0e-4; cf_veg.livestemc_to_deadstemc_patch .= 5.0e-5
    ns_veg = CLM.CNVegNitrogenStateData(); CLM.cnveg_nitrogen_state_init!(ns_veg, np, nc, ng; nrepr=nrepr)
    CLM.cnveg_nitrogen_state_set_values!(ns_veg, mask_soilp, 5.0, mask_soilc, 5.0; nrepr=nrepr)
    ns_veg.leafn_patch .= 5.0; ns_veg.frootn_patch .= 4.0; ns_veg.npool_patch .= 50.0; ns_veg.retransn_patch .= 2.0
    nf_veg = CLM.CNVegNitrogenFluxData()
    CLM.cnveg_nitrogen_flux_init!(nf_veg, np, nc, ng; nrepr=nrepr, nlevdecomp_full=nlevdecomp, ndecomp_pools=ndecomp_pools, i_litr_max=i_litr_max)
    CLM.cnveg_nitrogen_flux_set_values!(nf_veg, mask_soilp, 0.0, mask_soilc, 0.0; nrepr=nrepr, nlevdecomp_full=nlevdecomp, ndecomp_pools=ndecomp_pools, i_litr_max=i_litr_max)
    nf_veg.npool_to_leafn_patch .= 1.0e-4; nf_veg.npool_to_frootn_patch .= 5.0e-5; nf_veg.leafn_to_litter_patch .= 2.0e-5
    cs_soil = CLM.SoilBiogeochemCarbonStateData(); CLM.soil_bgc_carbon_state_init!(cs_soil, nc, ng, nlevdecomp, ndecomp_pools)
    CLM.soil_bgc_carbon_state_set_values!(cs_soil, mask_soilc, 100.0)
    cf_soil = CLM.SoilBiogeochemCarbonFluxData(); CLM.soil_bgc_carbon_flux_init!(cf_soil, nc, nlevdecomp, ndecomp_pools, ndecomp_cascade_transitions)
    CLM.soil_bgc_carbon_flux_set_values!(cf_soil, mask_soilc, 0.0)
    ns_soil = CLM.SoilBiogeochemNitrogenStateData(); CLM.soil_bgc_nitrogen_state_init!(ns_soil, nc, ng, nlevdecomp, ndecomp_pools)
    CLM.soil_bgc_nitrogen_state_set_values!(ns_soil, mask_soilc, 5.0)
    nf_soil = CLM.SoilBiogeochemNitrogenFluxData(); CLM.soil_bgc_nitrogen_flux_init!(nf_soil, nc, nlevdecomp, ndecomp_pools, ndecomp_cascade_transitions)
    CLM.soil_bgc_nitrogen_flux_set_values!(nf_soil, mask_soilc, 0.0)
    soilbgc_st = CLM.SoilBiogeochemStateData(); CLM.soil_bgc_state_init!(soilbgc_st, nc, np, nlevdecomp, ndecomp_cascade_transitions)
    cfg = CLM.CNDriverConfig(); cfg.use_matrixcn = true
    return (; config=cfg, cs_veg, cf_veg, ns_veg, nf_veg, cs_soil, cf_soil, ns_soil, nf_soil, soilbgc_st,
            cascade_donor_pool=[1,2,3,1,2], cascade_receiver_pool=[2,3,0,4,4],
            mask_soilc, mask_soilp, patch_column=[1,1,2,2,3,4], ivt=fill(1,np),
            woody=(w=zeros(80); w[2]=1.0; w), harvdate=fill(999,np), col_is_fates=fill(false,nc),
            nc, np, ng, dt, nlevdecomp, ndecomp_pools, ndecomp_cascade_transitions, nrepr, i_litr_min, i_litr_max, i_cwd)
end

function run_driver!(d, stc, stn)
    CLM.cn_driver_no_leaching!(d.config;
        mask_bgc_soilc=d.mask_soilc, mask_bgc_vegp=d.mask_soilp, bounds_col=1:d.nc, bounds_patch=1:d.np,
        nlevdecomp=d.nlevdecomp, ndecomp_pools=d.ndecomp_pools, ndecomp_cascade_transitions=d.ndecomp_cascade_transitions,
        i_litr_min=d.i_litr_min, i_litr_max=d.i_litr_max, i_cwd=d.i_cwd,
        patch_column=d.patch_column, ivt=d.ivt, woody=d.woody, harvdate=d.harvdate,
        col_is_fates=d.col_is_fates, cascade_donor_pool=d.cascade_donor_pool, cascade_receiver_pool=d.cascade_receiver_pool, dt=d.dt,
        cnveg_cs=d.cs_veg, cnveg_cf=d.cf_veg, cnveg_ns=d.ns_veg, cnveg_nf=d.nf_veg,
        soilbgc_cs=d.cs_soil, soilbgc_cf=d.cf_soil, soilbgc_ns=d.ns_soil, soilbgc_nf=d.nf_soil, soilbgc_state=d.soilbgc_st,
        veg_c_solve_state=stc, veg_n_solve_state=stn)
end

function main(backend)
    println("=" ^ 70); println("WHOLE-DRIVER GPU parity for cn_driver_no_leaching! (use_matrixcn)"); println("=" ^ 70)
    if backend === nothing; println("  No GPU backend."); return 0; end
    name, _, FT = backend
    @printf("  Backend: %s   (precision: %s)\n", name, FT)

    # Host warm-up: build the veg-C/N matrix structure into the solve states.
    W = make_data(); run_driver!(W, CLM.CNVegMatrixSolveState(), CLM.CNVegMatrixSolveState())
    stc = CLM.CNVegMatrixSolveState(); stn = CLM.CNVegMatrixSolveState()
    Wb = make_data(); run_driver!(Wb, stc, stn)            # this warm-up FILLS stc/stn

    # Golden host solve (fresh states → rebuild = same result).
    H = make_data(); run_driver!(H, CLM.CNVegMatrixSolveState(), CLM.CNVegMatrixSolveState())

    # Device solve: adapt state + the prebuilt solve states to Metal.
    mf(x) = CLM.Adapt.adapt(device_array_type(), CLM.Adapt.adapt(_F32(), x))
    mm(b) = device_array_type()(collect(b))
    B = make_data()
    D = (; config=B.config, cs_veg=mf(B.cs_veg), cf_veg=mf(B.cf_veg), ns_veg=mf(B.ns_veg), nf_veg=mf(B.nf_veg),
         cs_soil=mf(B.cs_soil), cf_soil=mf(B.cf_soil), ns_soil=mf(B.ns_soil), nf_soil=mf(B.nf_soil), soilbgc_st=mf(B.soilbgc_st),
         mask_soilc=mm(B.mask_soilc), mask_soilp=mm(B.mask_soilp), patch_column=mf(B.patch_column), ivt=mf(B.ivt),
         woody=mf(B.woody), harvdate=mf(B.harvdate), col_is_fates=mf(B.col_is_fates),
         cascade_donor_pool=mf(B.cascade_donor_pool), cascade_receiver_pool=mf(B.cascade_receiver_pool),
         nc=B.nc, np=B.np, ng=B.ng, dt=B.dt, nlevdecomp=B.nlevdecomp, ndecomp_pools=B.ndecomp_pools,
         ndecomp_cascade_transitions=B.ndecomp_cascade_transitions, nrepr=B.nrepr, i_litr_min=B.i_litr_min, i_litr_max=B.i_litr_max, i_cwd=B.i_cwd)
    if !(D.cs_veg.leafc_patch isa device_array_type()); println("  BLOCKED (adapt)."); return 2; end
    run_driver!(D, CLM.Adapt.adapt(_F32(), stc) |> x -> CLM.Adapt.adapt(device_array_type(), x),
                   CLM.Adapt.adapt(_F32(), stn) |> x -> CLM.Adapt.adapt(device_array_type(), x))

    checks = [("cs_veg.leafc", H.cs_veg.leafc_patch, D.cs_veg.leafc_patch),
              ("cs_veg.frootc", H.cs_veg.frootc_patch, D.cs_veg.frootc_patch),
              ("cs_veg.livestemc", H.cs_veg.livestemc_patch, D.cs_veg.livestemc_patch),
              ("cs_veg.deadstemc", H.cs_veg.deadstemc_patch, D.cs_veg.deadstemc_patch),
              ("cs_veg.leafc_xfer", H.cs_veg.leafc_xfer_patch, D.cs_veg.leafc_xfer_patch),
              ("ns_veg.leafn", H.ns_veg.leafn_patch, D.ns_veg.leafn_patch),
              ("ns_veg.retransn", H.ns_veg.retransn_patch, D.ns_veg.retransn_patch)]
    nfail = 0
    for (nm, a, b) in checks
        dd = reldiff(a, b); ok = dd < 1f-3
        @printf("  [%s] %-22s rel = %.3e\n", ok ? "PASS" : "FAIL", nm, dd); ok || (nfail += 1)
    end
    println(); println(nfail == 0 ? "  WHOLE veg matrix driver path MATCHES CPU ON $name ($FT)" : "  DIVERGENCE ($nfail).")
    return nfail == 0 ? 0 : 1
end
exit(main(detect_backend()))
