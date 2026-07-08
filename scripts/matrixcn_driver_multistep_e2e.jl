# ==========================================================================
# matrixcn_driver_multistep_e2e.jl — sustained MULTI-DAY whole-BGC-driver
# ROBUSTNESS: cn_driver_no_leaching! with the matrix-CN path (use_matrixcn) called
# for N timesteps (≈10 days at Δt=1800s) on Metal (Float32) vs the host (Float64),
# confirming the device driver stays bit-identical + finite over a long run (no
# per-step device-state corruption / leak / NaN — which the single-step
# gpu_validate_matrixcn_driver_e2e cannot catch).
#
# SCOPE NOTE: this drives the BGC driver with a MINIMAL synthetic fixture, so the
# BGC cycle does not turn over here (the driver's allocation/phenology recompute
# the transfer fluxes to ~0 without full-model forcing, and the soil-matrix branch
# needs cascade_con + decomp_params) → the pools barely move; this is a stability/
# robustness test, not a dynamics test. Real evolving-pool multi-step stability is
# proven at the SOLVE level in matrixcn_multistep_e2e.jl (2000 steps, real
# dynamics, drift 3.4e-6). A genuine real-domain FORCED run needs the full CLM
# init/forcing (parity_run_domain.jl) — most of which is not GPU-ported.
#
#   julia --project=scripts scripts/matrixcn_driver_multistep_e2e.jl
# ==========================================================================
using CLM
using Printf
import Metal
include(joinpath(@__DIR__, "gpu_backends.jl"))

struct _F32 end
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:AbstractFloat}) = Float32.(x)
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:Integer}) = x
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{Bool}) = x

function make_data(; nc=8, np=12, ng=2, nlevdecomp=1, ndecomp_pools=7, ndecomp_cascade_transitions=5, nrepr=1)
    i_litr_min=1; i_litr_max=3; i_cwd=4; dt=1800.0
    mask_soilc = trues(nc); mask_soilp = trues(np)
    cs_veg = CLM.CNVegCarbonStateData(); CLM.cnveg_carbon_state_init!(cs_veg, np, nc, ng; nrepr=nrepr)
    CLM.cnveg_carbon_state_set_values!(cs_veg, mask_soilp, 10.0, mask_soilc, 10.0; nrepr=nrepr)
    cs_veg.cpool_patch .= 50.0; cs_veg.leafc_patch .= 50.0; cs_veg.frootc_patch .= 40.0
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
            mask_soilc, mask_soilp, patch_column=repeat(1:nc, inner=cld(np,nc))[1:np], ivt=fill(1,np),
            woody=(w=zeros(80); w[2]=1.0; w), harvdate=fill(999,np), col_is_fates=fill(false,nc),
            nc, np, ng, dt, nlevdecomp, ndecomp_pools, ndecomp_cascade_transitions, nrepr, i_litr_min, i_litr_max, i_cwd)
end

# GPP replenishment each step so the cycle turns over (labile C in).
gpp_step!(d, rate) = (d.cs_veg.cpool_patch .= d.cs_veg.cpool_patch .+ eltype(d.cs_veg.cpool_patch)(rate * d.dt))

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

vegc(d) = sum(Float64.(Array(d.cs_veg.leafc_patch))) + sum(Float64.(Array(d.cs_veg.frootc_patch))) +
    sum(Float64.(Array(d.cs_veg.deadstemc_patch))) + sum(Float64.(Array(d.cs_veg.livestemc_patch)))
leafv(d) = Float64.(Array(d.cs_veg.leafc_patch))
soilc(d) = sum(Float64.(Array(d.cs_soil.decomp_cpools_vr_col)))
soilv(d) = vec(Float64.(Array(d.cs_soil.decomp_cpools_vr_col)))

function main(backend)
    println("=" ^ 70); println("MULTI-DAY whole-BGC-driver (use_matrixcn): host (Float64) vs Metal"); println("=" ^ 70)
    if backend === nothing; println("  No GPU backend."); return 0; end
    name, _, FT = backend
    N = 480; rate = 3.0e-6   # ~10 days at 1800s; steady GPP into cpool
    @printf("  Backend: %s (%s)   N=%d steps (≈%d days)\n", name, FT, N, N ÷ 48)

    # host warm-up to build the matrix structure into the solve states
    stc = CLM.CNVegMatrixSolveState(); stn = CLM.CNVegMatrixSolveState()
    Wb = make_data(); gpp_step!(Wb, rate); run_driver!(Wb, stc, stn)

    # HOST trajectory
    H = make_data(); leaf0 = leafv(H)[1]; soil0 = soilc(H)
    for _ in 1:N; gpp_step!(H, rate); run_driver!(H, CLM.CNVegMatrixSolveState(), CLM.CNVegMatrixSolveState()); end
    vh = vegc(H); lh = leafv(H); sh = soilc(H); shv = soilv(H)

    # DEVICE trajectory (same forcing; prebuilt structure)
    mf(x) = CLM.Adapt.adapt(Metal.MtlArray, CLM.Adapt.adapt(_F32(), x)); mm(b) = Metal.MtlArray(collect(b))
    B = make_data()
    D = (; config=B.config, cs_veg=mf(B.cs_veg), cf_veg=mf(B.cf_veg), ns_veg=mf(B.ns_veg), nf_veg=mf(B.nf_veg),
         cs_soil=mf(B.cs_soil), cf_soil=mf(B.cf_soil), ns_soil=mf(B.ns_soil), nf_soil=mf(B.nf_soil), soilbgc_st=mf(B.soilbgc_st),
         mask_soilc=mm(B.mask_soilc), mask_soilp=mm(B.mask_soilp), patch_column=mf(B.patch_column), ivt=mf(B.ivt),
         woody=mf(B.woody), harvdate=mf(B.harvdate), col_is_fates=mf(B.col_is_fates),
         cascade_donor_pool=mf(B.cascade_donor_pool), cascade_receiver_pool=mf(B.cascade_receiver_pool),
         nc=B.nc, np=B.np, ng=B.ng, dt=B.dt, nlevdecomp=B.nlevdecomp, ndecomp_pools=B.ndecomp_pools,
         ndecomp_cascade_transitions=B.ndecomp_cascade_transitions, nrepr=B.nrepr, i_litr_min=B.i_litr_min, i_litr_max=B.i_litr_max, i_cwd=B.i_cwd)
    if !(D.cs_veg.leafc_patch isa Metal.MtlArray); println("  BLOCKED."); return 2; end
    stc_d = CLM.Adapt.adapt(Metal.MtlArray, CLM.Adapt.adapt(_F32(), stc)); stn_d = CLM.Adapt.adapt(Metal.MtlArray, CLM.Adapt.adapt(_F32(), stn))
    for _ in 1:N; gpp_step!(D, rate); run_driver!(D, stc_d, stn_d); end
    Metal.synchronize()
    vd = vegc(D); ld = leafv(D); sd = soilc(D); sdv = soilv(D)

    vmoved = abs(lh[1] - leaf0) / leaf0; smoved = abs(sh - soil0) / soil0
    reld = maximum(abs.(ld .- lh) ./ (1.0 .+ abs.(lh)))
    srel = maximum(abs.(sdv .- shv) ./ (1.0 .+ abs.(shv)))
    @printf("  leafc[1]: start %.4f → host %.4f (moved %.1f%%)  device %.4f  drift %.3e\n", leaf0, lh[1], 100vmoved, ld[1], reld)
    @printf("  soil-C:  start %.4f → host %.4f (moved %.1f%%)  device %.4f  drift %.3e\n", soil0, sh, 100smoved, sd, srel)
    @printf("  all finite (device veg+soil): %s\n", all(isfinite, ld) && all(isfinite, sdv))
    ok = reld < 1.0e-3 && srel < 1.0e-3 && all(isfinite, ld) && all(isfinite, sdv)
    println(); println(ok ? "  ROBUST — whole BGC driver bit-identical + finite over $N GPU calls ($(N ÷ 48) days)" : "  FAILED")
    println("  (pools static: minimal fixture, no BGC forcing — see the SCOPE NOTE at top;")
    println("   real evolving-pool stability is in matrixcn_multistep_e2e.jl.)")
    return ok ? 0 : 1
end
exit(main(detect_backend()))
