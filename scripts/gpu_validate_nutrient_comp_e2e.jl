# ==========================================================================
# gpu_validate_nutrient_comp_e2e.jl — end-to-end Metal parity for the CN
# allocation pipeline calc_plant_nutrient_competition! -> calc_plant_cn_alloc!
# (per-patch C/N allocation; _cnalloc_main_kernel! ~69 args grouped into
# _CnAllocOut{V,M}/_CnAllocIn{V,M,VB} device-view bundles).
#
#   julia --project=scripts scripts/gpu_validate_nutrient_comp_e2e.jl
# ==========================================================================
using CLM
using Printf
import Metal
include(joinpath(@__DIR__, "gpu_backends.jl"))

struct _F32 end
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:AbstractFloat}) = Float32.(x)
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:Integer}) = x
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{Bool}) = x
struct _F32S end   # CropData pins ::FT scalars
CLM.Adapt.adapt_storage(::_F32S, x::AbstractArray{<:AbstractFloat}) = Float32.(x)
CLM.Adapt.adapt_storage(::_F32S, x::AbstractArray{<:Integer}) = x
CLM.Adapt.adapt_storage(::_F32S, x::AbstractArray{Bool}) = x
CLM.Adapt.adapt_storage(::_F32S, x::Float64) = Float32(x)

function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(Float64(A[i]) - Float64(B[i])) / (1.0 + max(abs(Float64(A[i])), abs(Float64(B[i])))))
    end
    return m
end
allfinite(a) = all(isfinite, Array(a))

function build(; np=2)
    nrepr = CLM.NREPR
    pftcon = CLM.PftConNutrientCompetition(
        woody=[1.0,1.0], froot_leaf=[1.0,1.0], croot_stem=[0.3,0.3], stem_leaf=[-1.0,-1.0],
        flivewd=[0.1,0.1], leafcn=[25.0,25.0], frootcn=[42.0,42.0], livewdcn=[50.0,50.0],
        deadwdcn=[500.0,500.0], fcur=[0.5,0.5], graincn=[50.0,50.0], grperc=[0.3,0.3],
        grpnow=[0.5,0.5], fleafcn=[65.0,65.0], ffrootcn=[100.0,100.0], fstemcn=[130.0,130.0],
        astemf=[0.0,0.0], season_decid=[0.0,0.0], stress_decid=[0.0,0.0])
    cn_shared_params = CLM.CNSharedParamsData(use_matrixcn=false, use_fun=false)
    patch = CLM.PatchData(); patch.column=[1,1]; patch.itype=[0,0]
    crop = CLM.CropData(); crop.croplive_patch=[false,false]
    cnveg_state = CLM.CNVegStateData()
    cnveg_state.c_allometry_patch=[1.5,1.5]; cnveg_state.n_allometry_patch=[0.06,0.06]
    cnveg_state.downreg_patch=zeros(np)
    cnveg_state.aleaf_patch=fill(NaN,np); cnveg_state.astem_patch=fill(NaN,np); cnveg_state.aroot_patch=fill(NaN,np)
    cnveg_state.arepr_patch=fill(NaN,np,nrepr)
    cnveg_cs = CLM.CNVegCarbonStateData()
    cnveg_cs.leafc_patch=[10.,10.]; cnveg_cs.frootc_patch=[5.,5.]; cnveg_cs.livestemc_patch=[20.,20.]
    cnveg_cf = CLM.CNVegCarbonFluxData()
    cnveg_cf.gpp_before_downreg_patch=[10.,10.]; cnveg_cf.availc_patch=[8.,8.]
    cnveg_cf.psnsun_to_cpool_patch=[6.,6.]; cnveg_cf.psnshade_to_cpool_patch=[4.,4.]
    cnveg_cf.annsum_npp_patch=[400.,400.]
    for f in (:excess_cflux_patch,:plant_calloc_patch,:cpool_to_leafc_patch,:cpool_to_leafc_storage_patch,
              :cpool_to_frootc_patch,:cpool_to_frootc_storage_patch,:cpool_to_livestemc_patch,
              :cpool_to_livestemc_storage_patch,:cpool_to_deadstemc_patch,:cpool_to_deadstemc_storage_patch,
              :cpool_to_livecrootc_patch,:cpool_to_livecrootc_storage_patch,:cpool_to_deadcrootc_patch,
              :cpool_to_deadcrootc_storage_patch,:cpool_to_gresp_storage_patch)
        setfield!(cnveg_cf, f, zeros(np))
    end
    cnveg_cf.cpool_to_reproductivec_patch=zeros(np,nrepr); cnveg_cf.cpool_to_reproductivec_storage_patch=zeros(np,nrepr)
    cnveg_nf = CLM.CNVegNitrogenFluxData()
    cnveg_nf.plant_ndemand_patch=[0.2,0.2]; cnveg_nf.retransn_to_npool_patch=[0.05,0.05]
    cnveg_nf.sminn_to_plant_fun_patch=zeros(np)
    for f in (:sminn_to_npool_patch,:plant_nalloc_patch,:npool_to_leafn_patch,:npool_to_leafn_storage_patch,
              :npool_to_frootn_patch,:npool_to_frootn_storage_patch,:npool_to_livestemn_patch,
              :npool_to_livestemn_storage_patch,:npool_to_deadstemn_patch,:npool_to_deadstemn_storage_patch,
              :npool_to_livecrootn_patch,:npool_to_livecrootn_storage_patch,:npool_to_deadcrootn_patch,
              :npool_to_deadcrootn_storage_patch)
        setfield!(cnveg_nf, f, zeros(np))
    end
    cnveg_nf.npool_to_reproductiven_patch=zeros(np,nrepr); cnveg_nf.npool_to_reproductiven_storage_patch=zeros(np,nrepr)
    return (; pftcon, cn_shared_params, patch, crop, cnveg_state, cnveg_cs, cnveg_cf, cnveg_nf,
        mask_soilp=BitVector(fill(true,np)), bounds=1:np, fpg_col=[0.8], nc=1)
end

run_c!(d) = CLM.calc_plant_nutrient_competition!(d.mask_soilp, d.bounds, d.pftcon, d.cn_shared_params,
    d.patch, d.crop, d.cnveg_state, d.cnveg_cs, d.cnveg_cf, d.cnveg_nf; fpg_col=d.fpg_col)

function main(backend)
    println("=" ^ 66); println("END-TO-END Metal parity for calc_plant_nutrient_competition!"); println("=" ^ 66)
    if backend === nothing; println("  No GPU backend."); return 0; end
    name, _, FT = backend
    @printf("  Backend: %s   (precision: %s)\n", name, FT)
    H = build(); B = build()
    run_c!(H)
    mf(x) = CLM.Adapt.adapt(Metal.MtlArray, CLM.Adapt.adapt(_F32(), x))
    mfs(x) = CLM.Adapt.adapt(Metal.MtlArray, CLM.Adapt.adapt(_F32S(), x))
    D = (; pftcon=mf(B.pftcon), cn_shared_params=B.cn_shared_params, patch=mf(B.patch), crop=mfs(B.crop),
         cnveg_state=mf(B.cnveg_state), cnveg_cs=mf(B.cnveg_cs), cnveg_cf=mf(B.cnveg_cf), cnveg_nf=mf(B.cnveg_nf),
         mask_soilp=Metal.MtlArray(collect(B.mask_soilp)), bounds=B.bounds, fpg_col=mf(B.fpg_col), nc=B.nc)
    if !(D.cnveg_cf.plant_calloc_patch isa Metal.MtlArray); println("  BLOCKED."); return 2; end
    run_c!(D)
    checks = [("plant_calloc", H.cnveg_cf.plant_calloc_patch, D.cnveg_cf.plant_calloc_patch),
              ("plant_nalloc", H.cnveg_nf.plant_nalloc_patch, D.cnveg_nf.plant_nalloc_patch),
              ("cpool_to_leafc", H.cnveg_cf.cpool_to_leafc_patch, D.cnveg_cf.cpool_to_leafc_patch),
              ("cpool_to_livestemc", H.cnveg_cf.cpool_to_livestemc_patch, D.cnveg_cf.cpool_to_livestemc_patch),
              ("cpool_to_deadcrootc", H.cnveg_cf.cpool_to_deadcrootc_patch, D.cnveg_cf.cpool_to_deadcrootc_patch),
              ("npool_to_leafn", H.cnveg_nf.npool_to_leafn_patch, D.cnveg_nf.npool_to_leafn_patch),
              ("npool_to_deadstemn", H.cnveg_nf.npool_to_deadstemn_patch, D.cnveg_nf.npool_to_deadstemn_patch),
              ("sminn_to_npool", H.cnveg_nf.sminn_to_npool_patch, D.cnveg_nf.sminn_to_npool_patch),
              ("excess_cflux", H.cnveg_cf.excess_cflux_patch, D.cnveg_cf.excess_cflux_patch),
              ("psnsun_to_cpool", H.cnveg_cf.psnsun_to_cpool_patch, D.cnveg_cf.psnsun_to_cpool_patch)]
    nfail = 0
    for (nm,a,b) in checks
        allfinite(a) || (@printf("  [WARN] %-22s non-finite\n", nm); continue)
        dd = reldiff(a,b); ok = dd < 1f-3
        @printf("  [%s] %-22s rel = %.3e\n", ok ? "PASS" : "FAIL", nm, dd); ok || (nfail += 1)
    end
    println(); println(nfail == 0 ? "  calc_plant_cn_alloc! MATCHES CPU ON $name ($FT)" : "  DIVERGENCE ($nfail).")
    return nfail == 0 ? 0 : 1
end
const BACKEND = detect_backend()
exit(main(BACKEND))
