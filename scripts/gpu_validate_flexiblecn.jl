# ==========================================================================
# gpu_validate_flexiblecn.jl — GPU parity for the FlexibleCN nutrient path:
# calc_plant_nutrient_demand_flexiblecn! (5 per-patch kernels) +
# calc_plant_nutrient_competition_flexiblecn! -> calc_plant_cn_alloc_flexiblecn!
# (one per-patch allocation kernel with _FlexAllocOut/_FlexAllocIn device-view
# bundles). Exercises carbon_resp_opt==1 (the flexible-C:N turnover branch).
#
#   julia --project=scripts scripts/gpu_validate_flexiblecn.jl
# ==========================================================================
using CLM
using Printf
include(joinpath(@__DIR__, "gpu_backends.jl"))
include(joinpath(@__DIR__, "gpu_adapt.jl"))
mf(x)  = mf(device_array_type(), x)
mfs(x) = mfs(device_array_type(), x)

function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(Float64(A[i]) - Float64(B[i])) / (1.0 + max(abs(Float64(A[i])), abs(Float64(B[i])))))
    end
    return m
end
allfinite(a) = all(isfinite, Array(a))

# Single non-crop woody patch with the extra fields FlexibleCN reads (mirrors
# test/test_flexiblecn.jl make_flexcn_data). leafc_storage/leafn_storage are set
# so actual_storage_leafcn > leafcn+15 → carbon_resp_opt=1 turns over extra C.
function build()
    np = 1; nc = 1; nlevdecomp = 1; nrepr = CLM.NREPR
    pftcon = CLM.PftConNutrientCompetition(
        woody=[1.0], froot_leaf=[1.0], croot_stem=[0.3], stem_leaf=[-1.0],
        flivewd=[0.1], leafcn=[25.0], frootcn=[42.0], livewdcn=[50.0], deadwdcn=[500.0],
        fcur=[0.5], graincn=[50.0], grperc=[0.3], grpnow=[0.5],
        fleafcn=[65.0], ffrootcn=[100.0], fstemcn=[130.0],
        astemf=[0.0], season_decid=[0.0], stress_decid=[0.0], evergreen=[0.0])
    cn_shared_params = CLM.CNSharedParamsData(use_matrixcn=false, use_fun=false, use_flexiblecn=true)
    patch = CLM.PatchData(); patch.column=[1]; patch.itype=[0]
    crop = CLM.CropData(); crop.croplive_patch=[false]
    vs = CLM.CNVegStateData()
    vs.c_allometry_patch=[1.5]; vs.n_allometry_patch=[0.06]; vs.downreg_patch=[0.0]
    vs.aleaf_patch=[NaN]; vs.astem_patch=[NaN]; vs.aroot_patch=[NaN]; vs.arepr_patch=fill(NaN,np,nrepr)
    vs.peaklai_patch=[0]; vs.tempsum_potential_gpp_patch=[0.0]; vs.annsum_potential_gpp_patch=[5000.0]
    vs.tempmax_retransn_patch=[0.0]; vs.annmax_retransn_patch=[1.0]; vs.grain_flag_patch=[0.0]
    cs = CLM.CNVegCarbonStateData()
    cs.leafc_patch=[10.0]; cs.frootc_patch=[5.0]; cs.livestemc_patch=[20.0]; cs.livecrootc_patch=[15.0]
    cs.leafc_storage_patch=[80.0]; cs.frootc_storage_patch=[40.0]  # high C:N → resp turnover
    cs.livestemc_storage_patch=[60.0]; cs.livecrootc_storage_patch=[50.0]
    cf = CLM.CNVegCarbonFluxData()
    cf.gpp_before_downreg_patch=[10.0]; cf.availc_patch=[8.0]; cf.npp_growth_patch=[8.0]
    cf.excess_cflux_patch=[0.0]; cf.plant_calloc_patch=[0.0]
    cf.psnsun_to_cpool_patch=[6.0]; cf.psnshade_to_cpool_patch=[4.0]; cf.annsum_npp_patch=[400.0]
    for fld in (:cpool_to_leafc_patch,:cpool_to_leafc_storage_patch,:cpool_to_frootc_patch,
                :cpool_to_frootc_storage_patch,:cpool_to_livestemc_patch,:cpool_to_livestemc_storage_patch,
                :cpool_to_deadstemc_patch,:cpool_to_deadstemc_storage_patch,:cpool_to_livecrootc_patch,
                :cpool_to_livecrootc_storage_patch,:cpool_to_deadcrootc_patch,:cpool_to_deadcrootc_storage_patch,
                :cpool_to_gresp_storage_patch,:cpool_to_resp_patch,:cpool_to_leafc_resp_patch,
                :cpool_to_leafc_storage_resp_patch,:cpool_to_frootc_resp_patch,:cpool_to_frootc_storage_resp_patch,
                :cpool_to_livecrootc_resp_patch,:cpool_to_livecrootc_storage_resp_patch,
                :cpool_to_livestemc_resp_patch,:cpool_to_livestemc_storage_resp_patch)
        setfield!(cf, fld, [0.0])
    end
    cf.cpool_to_reproductivec_patch=fill(0.0,np,nrepr); cf.cpool_to_reproductivec_storage_patch=fill(0.0,np,nrepr)
    ns = CLM.CNVegNitrogenStateData()
    ns.retransn_patch=[0.5]; ns.npool_patch=[0.4]; ns.leafn_patch=[0.4]; ns.frootn_patch=[0.12]; ns.livestemn_patch=[0.4]
    ns.leafn_storage_patch=[0.32]; ns.frootn_storage_patch=[0.1]; ns.livestemn_storage_patch=[0.2]; ns.livecrootn_storage_patch=[0.14]
    nf = CLM.CNVegNitrogenFluxData()
    nf.plant_ndemand_patch=[0.2]; nf.avail_retransn_patch=[0.0]; nf.retransn_to_npool_patch=[0.05]
    nf.sminn_to_npool_patch=[0.0]; nf.plant_nalloc_patch=[0.0]
    for fld in (:npool_to_leafn_patch,:npool_to_leafn_storage_patch,:npool_to_frootn_patch,
                :npool_to_frootn_storage_patch,:npool_to_livestemn_patch,:npool_to_livestemn_storage_patch,
                :npool_to_deadstemn_patch,:npool_to_deadstemn_storage_patch,:npool_to_livecrootn_patch,
                :npool_to_livecrootn_storage_patch,:npool_to_deadcrootn_patch,:npool_to_deadcrootn_storage_patch,
                :sminn_to_plant_fun_patch,:leafn_to_retransn_patch,:frootn_to_retransn_patch,:livestemn_to_retransn_patch)
        setfield!(nf, fld, [0.0])
    end
    nf.npool_to_reproductiven_patch=fill(0.0,np,nrepr); nf.npool_to_reproductiven_storage_patch=fill(0.0,np,nrepr)
    canopystate = CLM.CanopyStateData(); canopystate.laisun_patch=[1.0]; canopystate.laisha_patch=[0.5]
    soilbgc_ns = CLM.SoilBiogeochemNitrogenStateData(); soilbgc_ns.sminn_vr_col=fill(2.0,nc,nlevdecomp)
    soilbgc_cf = CLM.SoilBiogeochemCarbonFluxData(); soilbgc_cf.t_scalar_col=fill(0.7,nc,nlevdecomp)
    return (; pftcon, cn_shared_params, patch, crop, canopystate, cnveg_state=vs, cnveg_cs=cs,
        cnveg_cf=cf, cnveg_ns=ns, cnveg_nf=nf, soilbgc_ns, soilbgc_cf,
        mask_soilp=BitVector([true]), fpg_col=[0.8], bounds=1:np, dzsoi_decomp=[1.0], nlevdecomp)
end

function run_flex!(d, dt)
    CLM.calc_plant_nutrient_demand_flexiblecn!(d.mask_soilp, d.bounds, false,
        d.pftcon, d.cn_shared_params, d.patch, d.crop, d.canopystate,
        d.cnveg_state, d.cnveg_cs, d.cnveg_cf, d.cnveg_ns, d.cnveg_nf,
        d.soilbgc_ns, d.soilbgc_cf;
        dt=dt, dzsoi_decomp=d.dzsoi_decomp, nlevdecomp=d.nlevdecomp)
    CLM.calc_plant_nutrient_competition_flexiblecn!(d.mask_soilp, d.bounds,
        d.pftcon, d.cn_shared_params, d.patch, d.crop, d.canopystate,
        d.cnveg_state, d.cnveg_cs, d.cnveg_cf, d.cnveg_ns, d.cnveg_nf;
        fpg_col=d.fpg_col, dt=dt, carbon_resp_opt=1)
end

function main(backend)
    println("="^70); println("GPU parity — FlexibleCN nutrient demand + competition"); println("="^70)
    if backend === nothing; println("  No GPU backend."); return 0; end
    name, _, FT = backend
    @printf("  Backend: %s (%s)\n", name, FT)
    dt = 1800.0
    H = build(); B = build()
    run_flex!(H, dt)

    D = (; pftcon=mf(B.pftcon), cn_shared_params=B.cn_shared_params, patch=mf(B.patch),
        crop=mfs(B.crop), canopystate=mf(B.canopystate), cnveg_state=mfs(B.cnveg_state),
        cnveg_cs=mf(B.cnveg_cs), cnveg_cf=mf(B.cnveg_cf), cnveg_ns=mf(B.cnveg_ns), cnveg_nf=mf(B.cnveg_nf),
        soilbgc_ns=mf(B.soilbgc_ns), soilbgc_cf=mf(B.soilbgc_cf),
        mask_soilp=device_array_type()(collect(Bool, B.mask_soilp)), fpg_col=mf(B.fpg_col),
        bounds=B.bounds, dzsoi_decomp=B.dzsoi_decomp, nlevdecomp=B.nlevdecomp)
    if !(D.cnveg_cf.plant_calloc_patch isa device_array_type()); println("  BLOCKED: adapt."); return 2; end
    run_flex!(D, dt)

    checks = [
        ("plant_ndemand",    H.cnveg_nf.plant_ndemand_patch,    D.cnveg_nf.plant_ndemand_patch),
        ("avail_retransn",   H.cnveg_nf.avail_retransn_patch,   D.cnveg_nf.avail_retransn_patch),
        ("retransn_to_npool",H.cnveg_nf.retransn_to_npool_patch,D.cnveg_nf.retransn_to_npool_patch),
        ("plant_calloc",     H.cnveg_cf.plant_calloc_patch,     D.cnveg_cf.plant_calloc_patch),
        ("plant_nalloc",     H.cnveg_nf.plant_nalloc_patch,     D.cnveg_nf.plant_nalloc_patch),
        ("sminn_to_npool",   H.cnveg_nf.sminn_to_npool_patch,   D.cnveg_nf.sminn_to_npool_patch),
        ("cpool_to_leafc",   H.cnveg_cf.cpool_to_leafc_patch,   D.cnveg_cf.cpool_to_leafc_patch),
        ("cpool_to_livestemc",H.cnveg_cf.cpool_to_livestemc_patch,D.cnveg_cf.cpool_to_livestemc_patch),
        ("cpool_to_deadcrootc",H.cnveg_cf.cpool_to_deadcrootc_patch,D.cnveg_cf.cpool_to_deadcrootc_patch),
        ("cpool_to_gresp_storage",H.cnveg_cf.cpool_to_gresp_storage_patch,D.cnveg_cf.cpool_to_gresp_storage_patch),
        ("npool_to_leafn",   H.cnveg_nf.npool_to_leafn_patch,   D.cnveg_nf.npool_to_leafn_patch),
        ("npool_to_livestemn",H.cnveg_nf.npool_to_livestemn_patch,D.cnveg_nf.npool_to_livestemn_patch),
        ("cpool_to_resp",    H.cnveg_cf.cpool_to_resp_patch,    D.cnveg_cf.cpool_to_resp_patch),
        ("cpool_to_leafc_resp",H.cnveg_cf.cpool_to_leafc_resp_patch,D.cnveg_cf.cpool_to_leafc_resp_patch),
        ("cpool_to_frootc_resp",H.cnveg_cf.cpool_to_frootc_resp_patch,D.cnveg_cf.cpool_to_frootc_resp_patch),
    ]
    nfail = 0; ncmp = 0
    for (nm, a, b) in checks
        allfinite(a) || (@printf("  [WARN] %-24s non-finite host\n", nm); continue)
        dd = reldiff(a, b); ok = dd < 1f-3; ncmp += 1
        @printf("  [%s] %-24s rel = %.3e\n", ok ? "PASS" : "FAIL", nm, dd); ok || (nfail += 1)
    end
    println()
    println(nfail == 0 ? "  FlexibleCN kernels MATCH host ON $name over $ncmp fields" : "  DIVERGENCE ($nfail).")
    return nfail == 0 ? 0 : 1
end
exit(main(detect_backend()))
