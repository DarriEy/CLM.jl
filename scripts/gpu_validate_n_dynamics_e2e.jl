# ==========================================================================
# gpu_validate_n_dynamics_e2e.jl — end-to-end GPU parity for the N-dynamics
# column-input fluxes: n_deposition!, n_free_living_fixation!, n_fixation!,
# n_fert! (patch->col scatter), n_soyfix! (patch->col scatter). Exercises the
# per-column gather/fixation kernels + the two new _scatter_add! reductions.
#
#   julia --project=scripts scripts/gpu_validate_n_dynamics_e2e.jl
# ==========================================================================

using CLM
using Printf
include(joinpath(@__DIR__, "gpu_backends.jl"))

struct _F32 end   # arrays only — leaves scalar Float64 fields (N structs pin ::Float64)
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:AbstractFloat}) = Float32.(x)
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:Integer}) = x
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{Bool}) = x
struct _F32S end  # arrays + scalars — CropData pins its scalar fields to ::FT
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

const NTMP_SOY = 23

function build(; nc=4, np=8, ng=2)
    ndecomp_pools=7; ndecomp_cascade_transitions=10; nlevsoi=3
    nf = CLM.SoilBiogeochemNitrogenFluxData()
    CLM.soil_bgc_nitrogen_flux_init!(nf, nc, nlevsoi, ndecomp_pools, ndecomp_cascade_transitions)
    ns = CLM.SoilBiogeochemNitrogenStateData()
    CLM.soil_bgc_nitrogen_state_init!(ns, nc, 1, nlevsoi, ndecomp_pools); ns.sminn_col .= 20.0
    soilbgc_st = CLM.SoilBiogeochemStateData()
    CLM.soil_bgc_state_init!(soilbgc_st, nc, np, nlevsoi, ndecomp_cascade_transitions); soilbgc_st.fpg_col .= 0.5
    cf = CLM.CNVegCarbonFluxData(); cf.annsum_npp_col=fill(200.0,nc); cf.lag_npp_col=fill(5.0e-6,nc)
    cnveg_nf = CLM.CNVegNitrogenFluxData()
    cnveg_nf.fert_patch=fill(1.0e-7,np); cnveg_nf.soyfixn_patch=zeros(np); cnveg_nf.plant_ndemand_patch=fill(2.0e-7,np)
    cnveg_state = CLM.CNVegStateData(); cnveg_state.gddmaturity_patch=fill(1500.0,np)
    crop = CLM.CropData(); crop.hui_patch=fill(600.0,np); crop.croplive_patch=fill(true,np)
    wdiag = CLM.WaterDiagnosticBulkData(); wdiag.wf_col=fill(0.6,nc)
    wfb = CLM.WaterFluxBulkData(); wfb.AnnET=fill(2.0e-5,nc)
    col = CLM.ColumnData(); col.gridcell=[1,1,2,2]; col.is_fates=fill(false,nc)
    patch = CLM.PatchData(); patch.column=[1,1,2,2,3,3,4,4]; patch.wtcol=fill(0.5,np)
    patch.itype=[1,1,1,1,1,1, NTMP_SOY, NTMP_SOY]   # last 2 patches soybean -> exercise soyfix
    return (; nf, ns, soilbgc_st, cf, cnveg_nf, cnveg_state, crop, wdiag, wfb, col, patch,
        params=CLM.NDynamicsParams(), forc_ndep=[1.0e-8, 2.0e-8],
        mask_soilc=trues(nc), mask_soilp=trues(np), np, nc, ng)
end

function run_all!(d)
    CLM.n_deposition!(d.nf; forc_ndep=d.forc_ndep, col_gridcell=d.col.gridcell, bounds=1:d.nc)
    CLM.n_free_living_fixation!(d.nf, d.params; mask_soilc=d.mask_soilc, bounds=1:d.nc,
        AnnET=d.wfb.AnnET, dayspyr=365.0)
    CLM.n_fixation!(d.nf, d.cf; mask_soilc=d.mask_soilc, bounds=1:d.nc,
        col_is_fates=d.col.is_fates, dayspyr=365.0, nfix_timeconst=10.0, use_fun=false)
    CLM.n_fert!(d.nf, d.cnveg_nf; mask_soilc=d.mask_soilc, bounds=1:d.nc,
        patch=d.patch, mask_soilp=d.mask_soilp, bounds_p=1:d.np)
    CLM.n_soyfix!(d.nf, d.cnveg_nf, d.soilbgc_st, d.ns, d.cnveg_state, d.crop, d.wdiag;
        mask_soilc=d.mask_soilc, bounds=1:d.nc, mask_soilp=d.mask_soilp, bounds_p=1:d.np,
        patch=d.patch, ntmp_soybean=23, nirrig_tmp_soybean=58, ntrp_soybean=24, nirrig_trp_soybean=59)
end

const OUT = (:ndep_to_sminn_col, :ffix_to_sminn_col, :nfix_to_sminn_col,
             :fert_to_sminn_col, :soyfixn_to_sminn_col)

function main(backend)
    println("=" ^ 72); println("END-TO-END GPU parity for N-dynamics (deposition/fixation/fert/soyfix)"); println("=" ^ 72)
    if backend === nothing; println("  No GPU backend."); return 0; end
    name, _, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n", name, FT)
    H = build(); B = build()
    run_all!(H)
    mf(x)  = CLM.Adapt.adapt(device_array_type(), CLM.Adapt.adapt(_F32(), x))
    mfs(x) = CLM.Adapt.adapt(device_array_type(), CLM.Adapt.adapt(_F32S(), x))  # crop: ::FT scalars
    D = (; nf=mf(B.nf), ns=mf(B.ns), soilbgc_st=mf(B.soilbgc_st), cf=mf(B.cf),
        cnveg_nf=mf(B.cnveg_nf), cnveg_state=mfs(B.cnveg_state), crop=mfs(B.crop),
        wdiag=mf(B.wdiag), wfb=mf(B.wfb), col=mf(B.col), patch=mf(B.patch),
        params=B.params, forc_ndep=mf(B.forc_ndep),
        mask_soilc=device_array_type()(collect(B.mask_soilc)), mask_soilp=device_array_type()(collect(B.mask_soilp)),
        np=B.np, nc=B.nc, ng=B.ng)
    if !(D.nf.fert_to_sminn_col isa device_array_type()); println("  BLOCKED."); return 2; end
    run_all!(D)
    nfail = 0
    for f in OUT
        a = getfield(H.nf, f); b = getfield(D.nf, f)
        allfinite(a) || (@printf("  [WARN] %-24s non-finite\n", f); continue)
        dd = reldiff(a, b); ok = dd < 1f-3
        @printf("  [%s] %-24s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", f, dd)
        ok || (nfail += 1)
    end
    println()
    println(nfail == 0 ? "  N-dynamics col fluxes MATCH CPU ON $name ($FT)" :
                         "  DIVERGENCE — investigate ($nfail failed).")
    return nfail == 0 ? 0 : 1
end
const BACKEND = detect_backend()
exit(main(BACKEND))
