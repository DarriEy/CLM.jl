# ==========================================================================
# gpu_validate_ndemand_e2e.jl — end-to-end Metal parity for
# calc_plant_nitrogen_demand! (plant N demand + crop grain-fill retransn +
# avail-retransn kernels). Run with call_is_for_pcrop=true to exercise the
# device-resident crop_phase_vals scratch fix + crop_phase! on Metal.
#
#   julia --project=scripts scripts/gpu_validate_ndemand_e2e.jl
# ==========================================================================
using CLM
using Printf
import Metal
include(joinpath(@__DIR__, "gpu_backends.jl"))

struct _F32 end
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:AbstractFloat}) = Float32.(x)
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:Integer}) = x
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{Bool}) = x
struct _F32S end
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

function make_data()
    np = 1; nrepr = CLM.NREPR
    pftcon = CLM.PftConNutrientCompetition(
        woody=[1.0], froot_leaf=[1.0], croot_stem=[0.3], stem_leaf=[-1.0],
        flivewd=[0.1], leafcn=[25.0], frootcn=[42.0], livewdcn=[50.0],
        deadwdcn=[500.0], fcur=[0.5], graincn=[50.0], grperc=[0.3], grpnow=[0.5],
        fleafcn=[65.0], ffrootcn=[100.0], fstemcn=[130.0], astemf=[0.0],
        season_decid=[0.0], stress_decid=[0.0])
    cn_shared = CLM.CNSharedParamsData(use_matrixcn=false, use_fun=false)
    patch = CLM.PatchData(); patch.column=[1]; patch.itype=[0]
    crop = CLM.CropData(); crop.croplive_patch=[false]; crop.hui_patch=[0.0]
    crop.gddtsoi_patch=[0.0]
    vs = CLM.CNVegStateData()
    vs.c_allometry_patch=[1.5]; vs.n_allometry_patch=[0.06]; vs.downreg_patch=[0.0]
    vs.aleaf_patch=[NaN]; vs.astem_patch=[NaN]; vs.aroot_patch=[NaN]
    vs.arepr_patch=fill(NaN,np,nrepr); vs.peaklai_patch=[0]
    vs.tempsum_potential_gpp_patch=[0.0]; vs.annsum_potential_gpp_patch=[5000.0]
    vs.tempmax_retransn_patch=[0.0]; vs.annmax_retransn_patch=[1.0]; vs.grain_flag_patch=[0.0]
    vs.huileaf_patch=[0.0]; vs.huigrain_patch=[0.0]; vs.gddmaturity_patch=[0.0]
    cs = CLM.CNVegCarbonStateData()
    cs.leafc_patch=[10.0]; cs.frootc_patch=[5.0]; cs.livestemc_patch=[20.0]
    cf = CLM.CNVegCarbonFluxData()
    cf.gpp_before_downreg_patch=[10.0]; cf.availc_patch=[8.0]; cf.annsum_npp_patch=[400.0]
    ns = CLM.CNVegNitrogenStateData(); ns.retransn_patch=[0.5]
    nf = CLM.CNVegNitrogenFluxData()
    nf.plant_ndemand_patch=[0.2]; nf.avail_retransn_patch=[0.0]
    nf.retransn_to_npool_patch=[0.05]
    nf.leafn_to_retransn_patch=[0.0]; nf.frootn_to_retransn_patch=[0.0]
    nf.livestemn_to_retransn_patch=[0.0]
    return (; pftcon, cn_shared, patch, crop, vs, cs, cf, ns, nf,
            mask=trues(np), bounds=1:np, dt=1800.0)
end

run_nd!(d) = CLM.calc_plant_nitrogen_demand!(d.mask, d.bounds, true,
    d.pftcon, d.cn_shared, d.patch, d.crop, d.vs, d.cs, d.cf, d.ns, d.nf; dt=d.dt)

function main(backend)
    println("="^66); println("END-TO-END Metal parity for calc_plant_nitrogen_demand!"); println("="^66)
    if backend === nothing; println("  No GPU backend."); return 0; end
    name, _, FT = backend
    @printf("  Backend: %s   (precision: %s)\n", name, FT)
    H = make_data(); B = make_data()
    run_nd!(H)
    mf(x)  = CLM.Adapt.adapt(Metal.MtlArray, CLM.Adapt.adapt(_F32(), x))
    mfS(x) = CLM.Adapt.adapt(Metal.MtlArray, CLM.Adapt.adapt(_F32S(), x))
    D = (; pftcon=mf(B.pftcon), cn_shared=B.cn_shared, patch=mf(B.patch), crop=mfS(B.crop),
         vs=mf(B.vs), cs=mf(B.cs), cf=mf(B.cf), ns=mf(B.ns), nf=mf(B.nf),
         mask=Metal.MtlArray(collect(B.mask)), bounds=B.bounds, dt=B.dt)
    if !(D.nf.plant_ndemand_patch isa Metal.MtlArray); println("  BLOCKED."); return 2; end
    run_nd!(D)
    checks = [("plant_ndemand", H.nf.plant_ndemand_patch, D.nf.plant_ndemand_patch),
              ("tempsum_potential_gpp", H.vs.tempsum_potential_gpp_patch, D.vs.tempsum_potential_gpp_patch),
              ("tempmax_retransn", H.vs.tempmax_retransn_patch, D.vs.tempmax_retransn_patch),
              ("grain_flag", H.vs.grain_flag_patch, D.vs.grain_flag_patch),
              ("avail_retransn", H.nf.avail_retransn_patch, D.nf.avail_retransn_patch),
              ("retransn_to_npool", H.nf.retransn_to_npool_patch, D.nf.retransn_to_npool_patch),
              ("leafn_to_retransn", H.nf.leafn_to_retransn_patch, D.nf.leafn_to_retransn_patch)]
    nfail = 0
    for (nm, a, b) in checks
        dd = reldiff(a, b); ok = dd < 1f-3
        @printf("  [%s] %-22s rel = %.3e\n", ok ? "PASS" : "FAIL", nm, dd); ok || (nfail += 1)
    end
    println(); println(nfail == 0 ? "  calc_plant_nitrogen_demand! MATCHES CPU ON $name ($FT)" : "  DIVERGENCE ($nfail).")
    return nfail == 0 ? 0 : 1
end
const BACKEND = detect_backend()
exit(main(BACKEND))
