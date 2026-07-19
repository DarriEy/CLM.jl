# ==========================================================================
# gpu_validate_gap_mortality_e2e.jl — end-to-end GPU parity for cn_gap_mortality!
# (per-patch mortality fluxes, device-view-grouped) + cn_gap_patch_to_column!
# (3-level patch->column litter/CWD scatter via _scatter_add!).
#
#   julia --project=scripts scripts/gpu_validate_gap_mortality_e2e.jl
# ==========================================================================
using CLM
using Printf
include(joinpath(@__DIR__, "gpu_backends.jl"))

struct _F32 end     # arrays only
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:AbstractFloat}) = Float32.(x)
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:Integer}) = x
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{Bool}) = x
struct _F32S end    # arrays + scalars (GapMortalityParams pins k_mort::FT)
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

function build(; np=3, nc=2, nlevdecomp=1, n_litr=3)
    params = CLM.GapMortalityParams(k_mort=0.3, r_mort=fill(0.02, 21))
    pftcon = CLM.PftConGapMort(woody=vcat(fill(1.0,9),fill(0.0,12)), leafcn=fill(25.0,21),
        livewdcn=fill(50.0,21), lf_f=fill(1.0/n_litr,21,n_litr), fr_f=fill(1.0/n_litr,21,n_litr))
    dgvs = CLM.DgvsGapMortData(greffic_patch=fill(0.5,np), heatstress_patch=fill(0.0,np), nind_patch=fill(100.0,np))
    patch = CLM.PatchData(); patch.itype=[2,10,5]; patch.column=[1,1,2]; patch.wtcol=[0.5,0.3,0.2]
    canopystate = CLM.CanopyStateData(); canopystate.laisun_patch=fill(2.0,np); canopystate.laisha_patch=fill(1.0,np)
    cnveg_cs = CLM.CNVegCarbonStateData()
    for (f,v) in ((:leafc_patch,[10.,5.,8.]),(:frootc_patch,[4.,2.,3.]),(:livestemc_patch,[20.,0.,15.]),
                  (:deadstemc_patch,[50.,0.,40.]),(:livecrootc_patch,[10.,0.,8.]),(:deadcrootc_patch,[25.,0.,20.]),
                  (:leafc_storage_patch,[1.,.5,.8]),(:frootc_storage_patch,[.5,.2,.3]),(:livestemc_storage_patch,[2.,0.,1.5]),
                  (:deadstemc_storage_patch,[3.,0.,2.5]),(:livecrootc_storage_patch,[1.,0.,.8]),(:deadcrootc_storage_patch,[1.5,0.,1.2]),
                  (:gresp_storage_patch,[.2,.1,.15]),(:leafc_xfer_patch,[.5,.25,.4]),(:frootc_xfer_patch,[.2,.1,.15]),
                  (:livestemc_xfer_patch,[1.,0.,.8]),(:deadstemc_xfer_patch,[1.5,0.,1.2]),(:livecrootc_xfer_patch,[.5,0.,.4]),
                  (:deadcrootc_xfer_patch,[.8,0.,.6]),(:gresp_xfer_patch,[.1,.05,.08]))
        setfield!(cnveg_cs, f, Float64.(v))
    end
    cnveg_cf = CLM.CNVegCarbonFluxData()
    for f in (:m_leafc_to_litter_patch,:m_frootc_to_litter_patch,:m_livestemc_to_litter_patch,:m_deadstemc_to_litter_patch,
              :m_livecrootc_to_litter_patch,:m_deadcrootc_to_litter_patch,:m_leafc_storage_to_litter_patch,
              :m_frootc_storage_to_litter_patch,:m_livestemc_storage_to_litter_patch,:m_deadstemc_storage_to_litter_patch,
              :m_livecrootc_storage_to_litter_patch,:m_deadcrootc_storage_to_litter_patch,:m_gresp_storage_to_litter_patch,
              :m_leafc_xfer_to_litter_patch,:m_frootc_xfer_to_litter_patch,:m_livestemc_xfer_to_litter_patch,
              :m_deadstemc_xfer_to_litter_patch,:m_livecrootc_xfer_to_litter_patch,:m_deadcrootc_xfer_to_litter_patch,
              :m_gresp_xfer_to_litter_patch)
        setfield!(cnveg_cf, f, zeros(np))
    end
    cnveg_cf.gap_mortality_c_to_litr_c_col = zeros(nc,nlevdecomp,n_litr)
    cnveg_cf.gap_mortality_c_to_cwdc_col = zeros(nc,nlevdecomp)
    cnveg_ns = CLM.CNVegNitrogenStateData()
    for (f,v) in ((:leafn_patch,[.4,.2,.32]),(:frootn_patch,[.16,.08,.12]),(:livestemn_patch,[.4,0.,.3]),
                  (:deadstemn_patch,[.5,0.,.4]),(:livecrootn_patch,[.2,0.,.16]),(:deadcrootn_patch,[.25,0.,.2]),
                  (:retransn_patch,[.1,.05,.08]),(:leafn_storage_patch,[.04,.02,.032]),(:frootn_storage_patch,[.02,.008,.012]),
                  (:livestemn_storage_patch,[.04,0.,.03]),(:deadstemn_storage_patch,[.03,0.,.025]),(:livecrootn_storage_patch,[.02,0.,.016]),
                  (:deadcrootn_storage_patch,[.015,0.,.012]),(:leafn_xfer_patch,[.02,.01,.016]),(:frootn_xfer_patch,[.008,.004,.006]),
                  (:livestemn_xfer_patch,[.02,0.,.016]),(:deadstemn_xfer_patch,[.015,0.,.012]),(:livecrootn_xfer_patch,[.01,0.,.008]),
                  (:deadcrootn_xfer_patch,[.008,0.,.006]))
        setfield!(cnveg_ns, f, Float64.(v))
    end
    cnveg_nf = CLM.CNVegNitrogenFluxData()
    for f in (:m_leafn_to_litter_patch,:m_frootn_to_litter_patch,:m_livestemn_to_litter_patch,:m_deadstemn_to_litter_patch,
              :m_livecrootn_to_litter_patch,:m_deadcrootn_to_litter_patch,:m_retransn_to_litter_patch,
              :m_leafn_storage_to_litter_patch,:m_frootn_storage_to_litter_patch,:m_livestemn_storage_to_litter_patch,
              :m_deadstemn_storage_to_litter_patch,:m_livecrootn_storage_to_litter_patch,:m_deadcrootn_storage_to_litter_patch,
              :m_leafn_xfer_to_litter_patch,:m_frootn_xfer_to_litter_patch,:m_livestemn_xfer_to_litter_patch,
              :m_deadstemn_xfer_to_litter_patch,:m_livecrootn_xfer_to_litter_patch,:m_deadcrootn_xfer_to_litter_patch)
        setfield!(cnveg_nf, f, zeros(np))
    end
    cnveg_nf.gap_mortality_n_to_litr_n_col = zeros(nc,nlevdecomp,n_litr)
    cnveg_nf.gap_mortality_n_to_cwdn_col = zeros(nc,nlevdecomp)
    return (; params, pftcon, dgvs, patch, canopystate, cnveg_cs, cnveg_cf, cnveg_ns, cnveg_nf,
        mask=BitVector(fill(true,np)), bounds=1:np,
        leaf_prof=ones(np,nlevdecomp), froot_prof=ones(np,nlevdecomp),
        croot_prof=ones(np,nlevdecomp), stem_prof=ones(np,nlevdecomp), nc, nlevdecomp, n_litr)
end

function run_both!(d)
    CLM.cn_gap_mortality!(d.mask, d.bounds, d.params, d.pftcon, d.dgvs, d.patch, d.canopystate,
        d.cnveg_cs, d.cnveg_cf, d.cnveg_ns, d.cnveg_nf;
        dt=1800.0, days_per_year=365.0, use_cndv=false, use_matrixcn=false,
        spinup_state=0, npcropmin=17, spinup_factor_deadwood=1.0)
    CLM.cn_gap_patch_to_column!(d.mask, d.bounds, d.pftcon, d.patch, d.cnveg_cf, d.cnveg_nf,
        d.leaf_prof, d.froot_prof, d.croot_prof, d.stem_prof;
        nlevdecomp=d.nlevdecomp, i_litr_min=1, i_litr_max=d.n_litr, i_met_lit=1)
end

function main(backend)
    println("=" ^ 64); println("END-TO-END GPU parity for gap_mortality (+ patch->col scatter)"); println("=" ^ 64)
    if backend === nothing; println("  No GPU backend."); return 0; end
    name, _, FT = backend
    @printf("  Backend: %s   (precision: %s)\n", name, FT)
    H = build(); B = build()
    run_both!(H)
    mf(x) = CLM.Adapt.adapt(device_array_type(), CLM.Adapt.adapt(_F32(), x))
    mfs(x) = CLM.Adapt.adapt(device_array_type(), CLM.Adapt.adapt(_F32S(), x))
    D = (; params=mfs(B.params), pftcon=mf(B.pftcon), dgvs=mf(B.dgvs), patch=mf(B.patch),
         canopystate=mf(B.canopystate), cnveg_cs=mf(B.cnveg_cs), cnveg_cf=mf(B.cnveg_cf),
         cnveg_ns=mf(B.cnveg_ns), cnveg_nf=mf(B.cnveg_nf),
         mask=device_array_type()(collect(B.mask)), bounds=B.bounds,
         leaf_prof=mf(B.leaf_prof), froot_prof=mf(B.froot_prof), croot_prof=mf(B.croot_prof),
         stem_prof=mf(B.stem_prof), nc=B.nc, nlevdecomp=B.nlevdecomp, n_litr=B.n_litr)
    if !(D.cnveg_cf.m_leafc_to_litter_patch isa device_array_type()); println("  BLOCKED."); return 2; end
    run_both!(D)
    checks = [("m_leafc_to_litter", H.cnveg_cf.m_leafc_to_litter_patch, D.cnveg_cf.m_leafc_to_litter_patch),
              ("m_deadstemc_to_litter", H.cnveg_cf.m_deadstemc_to_litter_patch, D.cnveg_cf.m_deadstemc_to_litter_patch),
              ("m_leafn_to_litter", H.cnveg_nf.m_leafn_to_litter_patch, D.cnveg_nf.m_leafn_to_litter_patch),
              ("c_to_litr_c_col (scatter)", H.cnveg_cf.gap_mortality_c_to_litr_c_col, D.cnveg_cf.gap_mortality_c_to_litr_c_col),
              ("c_to_cwdc_col (scatter)", H.cnveg_cf.gap_mortality_c_to_cwdc_col, D.cnveg_cf.gap_mortality_c_to_cwdc_col),
              ("n_to_litr_n_col (scatter)", H.cnveg_nf.gap_mortality_n_to_litr_n_col, D.cnveg_nf.gap_mortality_n_to_litr_n_col),
              ("n_to_cwdn_col (scatter)", H.cnveg_nf.gap_mortality_n_to_cwdn_col, D.cnveg_nf.gap_mortality_n_to_cwdn_col)]
    nfail = 0
    for (nm,a,b) in checks
        allfinite(a) || (@printf("  [WARN] %-26s non-finite\n", nm); continue)
        dd = reldiff(a,b); ok = dd < 1f-3
        @printf("  [%s] %-26s rel = %.3e\n", ok ? "PASS" : "FAIL", nm, dd); ok || (nfail += 1)
    end
    println(); println(nfail == 0 ? "  gap_mortality + patch->col MATCH CPU ON $name ($FT)" : "  DIVERGENCE ($nfail).")
    return nfail == 0 ? 0 : 1
end
const BACKEND = detect_backend()
exit(main(BACKEND))
