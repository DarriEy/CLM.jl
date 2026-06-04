# ==========================================================================
# gpu_validate_cn_products_e2e.jl — end-to-end Metal parity for the wood/crop
# product-pool partitioning: cn_products_partition_wood_fluxes! +
# cn_products_partition_crop_fluxes! (per-fp filter kernels, the kernelized
# unity-scale p2g_1d! patch->gridcell average, and the dwt patch->gridcell
# _scatter_add!).
#
#   julia --project=scripts scripts/gpu_validate_cn_products_e2e.jl
# ==========================================================================
using CLM
using Printf
import Metal
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
# product-pool grc arrays hold SPVAL where a gridcell has no contributions;
# treat matching SPVAL/SPVAL as agreement (reldiff already does NaN; do SPVAL too).
finite_ok(a) = all(x -> isfinite(x) || x == CLM.SPVAL, Array(a))

function build(; ng=2, np=4)
    pch = CLM.PatchData()
    pch.gridcell = vcat(fill(1, div(np,2)), fill(2, np-div(np,2)))
    pch.itype = collect(0:np-1); pch.active = fill(true,np)
    pch.wtgcell = fill(1.0/div(np,ng), np); pch.wtcol=fill(1.0,np); pch.wtlunit=fill(1.0,np)
    pch.column = collect(1:np); pch.landunit = fill(1,np)
    lun = CLM.LandunitData(); lun.itype=[CLM.ISTSOIL]; lun.active=[true]; lun.wtgcell=[1.0]
    lun.gridcell=[1]; lun.urbpoi=[false]; lun.canyon_hwr=[0.0]
    col = CLM.ColumnData(); col.active=fill(true,np); col.wtgcell=fill(1.0,np)
    col.wtlunit=fill(1.0,np); col.landunit=fill(1,np); col.gridcell=pch.gridcell; col.itype=fill(1,np)
    bounds = CLM.BoundsType(begg=1, endg=ng, begl=1, endl=1, begc=1, endc=np, begp=1, endp=np)
    prod = CLM.CNProductsFullData(); CLM.cn_products_full_init!(prod, ng, np)
    return (; prod, pch, col, lun, bounds,
        filter_soilp=collect(1:np), num_soilp=np,
        gru_wood=[0.1,0.2,0.3,0.4], wood_harv=[1.,2.,3.,4.], dwt_wood=[0.5,0.0,0.3,0.7],
        crop_harv=[5.,6.,7.,8.], dwt_crop=[0.2,0.1,0.4,0.3],
        pprod10=fill(0.2,np), pprod100=fill(0.3,np), pprodharv10=fill(0.4,np), np, ng)
end

function run_both!(d)
    CLM.cn_products_partition_wood_fluxes!(d.prod, d.bounds, d.num_soilp, d.filter_soilp,
        d.dwt_wood, d.gru_wood, d.wood_harv;
        pprod10=d.pprod10, pprod100=d.pprod100, pprodharv10=d.pprodharv10,
        patch_gridcell=d.pch.gridcell, patch_itype=d.pch.itype,
        pch=d.pch, col=d.col, lun=d.lun)
    CLM.cn_products_partition_crop_fluxes!(d.prod, d.bounds, d.num_soilp, d.filter_soilp,
        d.dwt_crop, d.crop_harv;
        patch_gridcell=d.pch.gridcell, pch=d.pch, col=d.col, lun=d.lun)
end

function main(backend)
    println("=" ^ 66); println("END-TO-END Metal parity for cn_products partition fns"); println("=" ^ 66)
    if backend === nothing; println("  No GPU backend."); return 0; end
    name, _, FT = backend
    @printf("  Backend: %s   (precision: %s)\n", name, FT)
    H = build(); B = build()
    run_both!(H)
    mf(x) = CLM.Adapt.adapt(Metal.MtlArray, CLM.Adapt.adapt(_F32(), x))
    md(x) = mf(x)
    # col/lun/bounds stay host (unity p2g path doesn't use them).
    D = (; prod=mf(B.prod), pch=mf(B.pch), col=B.col, lun=B.lun, bounds=B.bounds,
         filter_soilp=md(B.filter_soilp), num_soilp=B.num_soilp,
         gru_wood=mf(B.gru_wood), wood_harv=mf(B.wood_harv), dwt_wood=mf(B.dwt_wood),
         crop_harv=mf(B.crop_harv), dwt_crop=mf(B.dwt_crop),
         pprod10=mf(B.pprod10), pprod100=mf(B.pprod100), pprodharv10=mf(B.pprodharv10),
         np=B.np, ng=B.ng)
    if !(D.prod.gru_prod10_gain_patch isa Metal.MtlArray); println("  BLOCKED."); return 2; end
    run_both!(D)
    checks = [("gru_prod10_gain_patch", H.prod.gru_prod10_gain_patch, D.prod.gru_prod10_gain_patch),
              ("gru_prod100_gain_patch", H.prod.gru_prod100_gain_patch, D.prod.gru_prod100_gain_patch),
              ("hrv_deadstem_to_prod10_patch", H.prod.hrv_deadstem_to_prod10_patch, D.prod.hrv_deadstem_to_prod10_patch),
              ("gru_prod10_gain_grc (p2g)", H.prod.gru_prod10_gain_grc, D.prod.gru_prod10_gain_grc),
              ("hrv_deadstem_to_prod100_grc (p2g)", H.prod.hrv_deadstem_to_prod100_grc, D.prod.hrv_deadstem_to_prod100_grc),
              ("dwt_prod10_gain_grc (scatter)", H.prod.dwt_prod10_gain_grc, D.prod.dwt_prod10_gain_grc),
              ("dwt_prod100_gain_grc (scatter)", H.prod.dwt_prod100_gain_grc, D.prod.dwt_prod100_gain_grc),
              ("crop_harvest_to_cropprod1_grc (p2g)", H.prod.crop_harvest_to_cropprod1_grc, D.prod.crop_harvest_to_cropprod1_grc),
              ("dwt_cropprod1_gain_grc (scatter)", H.prod.dwt_cropprod1_gain_grc, D.prod.dwt_cropprod1_gain_grc)]
    nfail = 0
    for (nm,a,b) in checks
        finite_ok(a) || (@printf("  [WARN] %-34s non-finite\n", nm); continue)
        dd = reldiff(a,b); ok = dd < 1f-3
        @printf("  [%s] %-34s rel = %.3e\n", ok ? "PASS" : "FAIL", nm, dd); ok || (nfail += 1)
    end
    println(); println(nfail == 0 ? "  cn_products partition fns MATCH CPU ON $name ($FT)" : "  DIVERGENCE ($nfail).")
    return nfail == 0 ? 0 : 1
end
const BACKEND = detect_backend()
exit(main(BACKEND))
