# ==========================================================================
# gpu_validate_alloc_pipeline_e2e.jl — end-to-end GPU parity for the three
# already-kernelized allocation.jl functions:
#   calc_gpp_mr_availc!            (GPP / maintenance-resp / available-C)
#   calc_allometry!                (C/N allometry coefficients)
#   calc_crop_allocation_fractions!(crop aleaf/astem/aroot/arepr fractions)
#
#   julia --project=scripts scripts/gpu_validate_alloc_pipeline_e2e.jl
# ==========================================================================
using CLM
using Printf
include(joinpath(@__DIR__, "gpu_backends.jl"))

# _F32: convert float arrays only (for structs with concrete ::Float64 scalar
# fields such as CanopyStateData.leaf_mr_vcm — those must stay Float64).
struct _F32 end
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:AbstractFloat}) = Float32.(x)
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:Integer}) = x
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{Bool}) = x
# _F32S: also convert Float64 scalar fields (CropData pins ::FT scalars that must
# move in lockstep with their arrays).
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

function make_woody(; psnsun, psnsha, laisun, laisha, leaf_mr, froot_mr,
                    livestem_mr, livecroot_mr, xsmrpool)
    np = 1; nrepr = CLM.NREPR
    alloc_params = CLM.AllocationParams(dayscrecover=30.0)
    pftcon = CLM.PftConAllocation(
        woody=[1.0], froot_leaf=[1.0], croot_stem=[0.3], stem_leaf=[-1.0],
        flivewd=[0.1], leafcn=[25.0], frootcn=[42.0], livewdcn=[50.0],
        deadwdcn=[500.0], graincn=[50.0], grperc=[0.3], arooti=[0.0], arootf=[0.0],
        bfact=[0.0], fleafi=[0.0], aleaff=[0.0], astemf=[0.0], allconss=[0.0],
        allconsl=[0.0], declfact=[0.0])
    cn_shared = CLM.CNSharedParamsData(use_matrixcn=false, use_fun=false)
    patch = CLM.PatchData(); CLM.patch_init!(patch, np); patch.itype[1] = 0
    crop = CLM.CropData(); CLM.crop_init!(crop, np); crop.croplive_patch[1] = false
    photo = CLM.PhotosynthesisData(); CLM.photosynthesis_data_init!(photo, np)
    photo.psnsun_patch[1] = psnsun; photo.psnsha_patch[1] = psnsha
    photo.c13_psnsun_patch[1] = 0.0; photo.c13_psnsha_patch[1] = 0.0
    photo.c14_psnsun_patch[1] = 0.0; photo.c14_psnsha_patch[1] = 0.0
    cstate = CLM.CanopyStateData(); CLM.canopystate_init!(cstate, np)
    cstate.laisun_patch[1] = laisun; cstate.laisha_patch[1] = laisha
    cs = CLM.CNVegCarbonStateData(); CLM.cnveg_carbon_state_init!(cs, np, 1, 1)
    cs.xsmrpool_patch[1] = xsmrpool
    cf = CLM.CNVegCarbonFluxData(); CLM.cnveg_carbon_flux_init!(cf, np, 1, 1)
    cf.leaf_mr_patch[1] = leaf_mr; cf.froot_mr_patch[1] = froot_mr
    cf.livestem_mr_patch[1] = livestem_mr; cf.livecroot_mr_patch[1] = livecroot_mr
    for k in 1:nrepr; cf.reproductive_mr_patch[1, k] = 0.0; end
    cf.annsum_npp_patch[1] = 400.0; cf.xsmrpool_recover_patch[1] = 0.0
    cf.cpool_to_xsmrpool_patch[1] = 0.0
    vs = CLM.CNVegStateData(); CLM.cnveg_state_init!(vs, np, 1)
    return (; alloc_params, pftcon, cn_shared, patch, crop, photo, cstate, cs, cf, vs,
            mask=BitVector([true]), bounds=1:np)
end

run_gpp!(d) = CLM.calc_gpp_mr_availc!(d.mask, d.bounds, d.alloc_params, d.pftcon,
    d.cn_shared, d.patch, d.crop, d.photo, d.cstate, d.cs, d.cf)
run_allom!(d) = CLM.calc_allometry!(d.mask, d.bounds, d.pftcon, d.cn_shared, d.patch, d.cf, d.vs)
run_crop!(d) = CLM.calc_crop_allocation_fractions!(d.mask, d.bounds, d.pftcon, d.patch, d.crop, d.vs)

function main(backend)
    println("="^66); println("END-TO-END GPU parity for allocation.jl (gpp/availc)"); println("="^66)
    if backend === nothing; println("  No GPU backend."); return 0; end
    name, _, FT = backend
    @printf("  Backend: %s   (precision: %s)\n", name, FT)
    cfg = (; psnsun=5.0, psnsha=2.0, laisun=2.0, laisha=3.0, leaf_mr=0.001,
           froot_mr=0.0005, livestem_mr=0.0003, livecroot_mr=0.0002, xsmrpool=-100.0)
    H = make_woody(; cfg...); B = make_woody(; cfg...)
    run_gpp!(H); run_allom!(H); run_crop!(H)
    mf(x)  = CLM.Adapt.adapt(device_array_type(), CLM.Adapt.adapt(_F32(), x))
    mfS(x) = CLM.Adapt.adapt(device_array_type(), CLM.Adapt.adapt(_F32S(), x))
    D = (; alloc_params=B.alloc_params, pftcon=mf(B.pftcon), cn_shared=B.cn_shared,
         patch=mf(B.patch), crop=mfS(B.crop), photo=mf(B.photo), cstate=mf(B.cstate),
         cs=mf(B.cs), cf=mf(B.cf), vs=mf(B.vs), mask=device_array_type()(collect(B.mask)), bounds=B.bounds)
    if !(D.cf.availc_patch isa device_array_type()); println("  BLOCKED."); return 2; end
    run_gpp!(D); run_allom!(D); run_crop!(D)
    checks = [("psnsun_to_cpool", H.cf.psnsun_to_cpool_patch, D.cf.psnsun_to_cpool_patch),
              ("psnshade_to_cpool", H.cf.psnshade_to_cpool_patch, D.cf.psnshade_to_cpool_patch),
              ("gpp_before_downreg", H.cf.gpp_before_downreg_patch, D.cf.gpp_before_downreg_patch),
              ("availc", H.cf.availc_patch, D.cf.availc_patch),
              ("leaf_curmr", H.cf.leaf_curmr_patch, D.cf.leaf_curmr_patch),
              ("leaf_xsmr", H.cf.leaf_xsmr_patch, D.cf.leaf_xsmr_patch),
              ("froot_curmr", H.cf.froot_curmr_patch, D.cf.froot_curmr_patch),
              ("xsmrpool_recover", H.cf.xsmrpool_recover_patch, D.cf.xsmrpool_recover_patch),
              ("allom: c_allometry", H.vs.c_allometry_patch, D.vs.c_allometry_patch),
              ("allom: n_allometry", H.vs.n_allometry_patch, D.vs.n_allometry_patch),
              ("crop: aleaf", H.vs.aleaf_patch, D.vs.aleaf_patch),
              ("crop: astem", H.vs.astem_patch, D.vs.astem_patch),
              ("crop: aroot", H.vs.aroot_patch, D.vs.aroot_patch),
              ("crop: arepr", H.vs.arepr_patch, D.vs.arepr_patch)]
    nfail = 0
    for (nm, a, b) in checks
        dd = reldiff(a, b); ok = dd < 1f-3
        @printf("  [%s] %-22s rel = %.3e\n", ok ? "PASS" : "FAIL", nm, dd); ok || (nfail += 1)
    end
    println(); println(nfail == 0 ? "  gpp/availc + allometry + crop_fractions MATCH CPU ON $name ($FT)" : "  DIVERGENCE ($nfail).")
    return nfail == 0 ? 0 : 1
end
const BACKEND = detect_backend()
exit(main(BACKEND))
