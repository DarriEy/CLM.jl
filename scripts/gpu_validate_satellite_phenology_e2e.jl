# ==========================================================================
# gpu_validate_satellite_phenology_e2e.jl — end-to-end Metal parity for the
# prescribed-LAI (satellite phenology) path: satellite_phenology! (LAI/SAI/
# height interpolation + snow burial), read_annual_vegetation!, and
# read_monthly_vegetation! — all per-patch kernels.
#
#   julia --project=scripts scripts/gpu_validate_satellite_phenology_e2e.jl
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
allfinite(a) = all(isfinite, Array(a))

function build(; np=4, nc=2, ng=1, maxveg=78)
    sp = CLM.SatellitePhenologyData(); CLM.satellite_phenology_init!(sp, np)
    sp.timwt[1]=0.6; sp.timwt[2]=0.4
    sp.mlai2t[1,:].=[0.,0.]; sp.mlai2t[2,:].=[3.,4.]; sp.mlai2t[3,:].=[1.5,2.]; sp.mlai2t[4,:].=[2.,2.5]
    sp.msai2t[1,:].=[0.,0.]; sp.msai2t[2,:].=[.5,.6]; sp.msai2t[3,:].=[.3,.4]; sp.msai2t[4,:].=[.4,.5]
    sp.mhvt2t[1,:].=[0.,0.]; sp.mhvt2t[2,:].=[17.,18.]; sp.mhvt2t[3,:].=[.5,.6]; sp.mhvt2t[4,:].=[.3,.4]
    sp.mhvb2t[1,:].=[0.,0.]; sp.mhvb2t[2,:].=[8.,9.]; sp.mhvb2t[3,:].=[.1,.15]; sp.mhvb2t[4,:].=[.01,.02]
    cs = CLM.CanopyStateData(); CLM.canopystate_init!(cs, np)
    wd = CLM.WaterDiagnosticBulkData(); wd.frac_sno_col=[0.0,0.3]; wd.snow_depth_col=[0.0,0.15]
    patch = CLM.PatchData(); CLM.patch_init!(patch, np)
    patch.itype .= [0,1,11,12]; patch.column .= [1,1,2,2]; patch.gridcell .= [1,1,1,1]
    # monthly_* are [g, l+1, month]
    mlai = fill(2.0, ng, maxveg+1, 12); msai = fill(0.5, ng, maxveg+1, 12)
    mht = fill(15.0, ng, maxveg+1, 12); mhb = fill(7.0, ng, maxveg+1, 12)
    return (; sp, cs, wd, patch, mlai, msai, mht, mhb,
        mask=BitVector(fill(true,np)), bounds=1:np, np, nc, ng)
end

function run_all!(d)
    CLM.read_monthly_vegetation!(d.sp, d.cs, d.patch, d.bounds;
        monthly_lai=d.mlai, monthly_sai=d.msai, monthly_height_top=d.mht,
        monthly_height_bot=d.mhb, months=(1,2), noveg=0, maxveg=78)
    CLM.read_annual_vegetation!(d.cs, d.patch, d.bounds; monthly_lai=d.mlai, noveg=0, maxveg=78)
    CLM.satellite_phenology!(d.sp, d.cs, d.wd, d.patch, d.mask, d.bounds;
        noveg=0, nbrdlf_dcd_brl_shrub=11, use_lai_streams=false, use_fates_sp=false)
end

function main(backend)
    println("=" ^ 66); println("END-TO-END Metal parity for satellite_phenology! (prescribed LAI)"); println("=" ^ 66)
    if backend === nothing; println("  No GPU backend."); return 0; end
    name, _, FT = backend
    @printf("  Backend: %s   (precision: %s)\n", name, FT)
    H = build(); B = build()
    run_all!(H)
    mf(x) = CLM.Adapt.adapt(Metal.MtlArray, CLM.Adapt.adapt(_F32(), x))
    D = (; sp=mf(B.sp), cs=mf(B.cs), wd=mf(B.wd), patch=mf(B.patch),
         mlai=mf(B.mlai), msai=mf(B.msai), mht=mf(B.mht), mhb=mf(B.mhb),
         mask=Metal.MtlArray(collect(B.mask)), bounds=B.bounds, np=B.np, nc=B.nc, ng=B.ng)
    if !(D.cs.elai_patch isa Metal.MtlArray); println("  BLOCKED."); return 2; end
    run_all!(D)
    checks = [("tlai", H.cs.tlai_patch, D.cs.tlai_patch),
              ("tsai", H.cs.tsai_patch, D.cs.tsai_patch),
              ("htop", H.cs.htop_patch, D.cs.htop_patch),
              ("elai (snow burial)", H.cs.elai_patch, D.cs.elai_patch),
              ("esai (snow burial)", H.cs.esai_patch, D.cs.esai_patch),
              ("frac_veg_nosno_alb", H.cs.frac_veg_nosno_alb_patch, D.cs.frac_veg_nosno_alb_patch),
              ("annlai (read_annual)", H.cs.annlai_patch, D.cs.annlai_patch),
              ("mlai2t (read_monthly)", H.sp.mlai2t, D.sp.mlai2t),
              ("mlaidiff", H.cs.mlaidiff_patch, D.cs.mlaidiff_patch)]
    nfail = 0
    for (nm,a,b) in checks
        allfinite(a) || (@printf("  [WARN] %-22s non-finite\n", nm); continue)
        dd = reldiff(a,b); ok = dd < 1f-3
        @printf("  [%s] %-22s rel = %.3e\n", ok ? "PASS" : "FAIL", nm, dd); ok || (nfail += 1)
    end
    println(); println(nfail == 0 ? "  satellite_phenology! pipeline MATCHES CPU ON $name ($FT)" : "  DIVERGENCE ($nfail).")
    return nfail == 0 ? 0 : 1
end
const BACKEND = detect_backend()
exit(main(BACKEND))
