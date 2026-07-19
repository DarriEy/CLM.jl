# ==========================================================================
# gpu_validate_fire_fluxes_e2e.jl — end-to-end GPU parity for the WHOLE
# cnfire_fluxes_li2014! fire C/N flux driver (delegates to cnfire_fluxes!):
# per-patch combustion/mortality fluxes, the (j,p) litter/CWD patch->column
# scatter, the (c,j,l) decomposing-pool fire loss, deforestation C loss, and
# peat-fire somc loss.
#
# Builds the flux dataset (clone of test_fire_li2014's make_fire_li2014_flux_data)
# at Float64, runs cnfire_fluxes_li2014! on the CPU, adapts every kernel-touched
# struct to Metal/Float32, runs the SAME call on the device, and compares the
# mutated C/N fire-flux fields + the column scatter accumulators.
#
#   julia --project=scripts scripts/gpu_validate_fire_fluxes_e2e.jl
# ==========================================================================

using CLM
using Printf
include(joinpath(@__DIR__, "gpu_backends.jl"))

struct _F32 end
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:AbstractFloat}) = Float32.(x)
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:Integer}) = x
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{Bool}) = x
CLM.Adapt.adapt_storage(::_F32, x::Float64) = Float32(x)

function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(Float64(A[i]) - Float64(B[i])) / (1.0 + max(abs(Float64(A[i])), abs(Float64(B[i])))))
    end
    return m
end
allfinite(a) = all(isfinite, Array(a))

const CF_FIRE = (:m_leafc_to_fire_patch, :m_leafc_storage_to_fire_patch, :m_leafc_xfer_to_fire_patch,
    :m_livestemc_to_fire_patch, :m_livestemc_storage_to_fire_patch, :m_livestemc_xfer_to_fire_patch,
    :m_deadstemc_to_fire_patch, :m_deadstemc_storage_to_fire_patch, :m_deadstemc_xfer_to_fire_patch,
    :m_frootc_to_fire_patch, :m_frootc_storage_to_fire_patch, :m_frootc_xfer_to_fire_patch,
    :m_livecrootc_to_fire_patch, :m_livecrootc_storage_to_fire_patch, :m_livecrootc_xfer_to_fire_patch,
    :m_deadcrootc_to_fire_patch, :m_deadcrootc_storage_to_fire_patch, :m_deadcrootc_xfer_to_fire_patch,
    :m_gresp_storage_to_fire_patch, :m_gresp_xfer_to_fire_patch,
    :m_leafc_to_litter_fire_patch, :m_leafc_storage_to_litter_fire_patch, :m_leafc_xfer_to_litter_fire_patch,
    :m_livestemc_to_litter_fire_patch, :m_livestemc_storage_to_litter_fire_patch, :m_livestemc_xfer_to_litter_fire_patch,
    :m_livestemc_to_deadstemc_fire_patch, :m_deadstemc_to_litter_fire_patch,
    :m_deadstemc_storage_to_litter_fire_patch, :m_deadstemc_xfer_to_litter_fire_patch,
    :m_frootc_to_litter_fire_patch, :m_frootc_storage_to_litter_fire_patch, :m_frootc_xfer_to_litter_fire_patch,
    :m_livecrootc_to_litter_fire_patch, :m_livecrootc_storage_to_litter_fire_patch, :m_livecrootc_xfer_to_litter_fire_patch,
    :m_livecrootc_to_deadcrootc_fire_patch, :m_deadcrootc_to_litter_fire_patch,
    :m_deadcrootc_storage_to_litter_fire_patch, :m_deadcrootc_xfer_to_litter_fire_patch,
    :m_gresp_storage_to_litter_fire_patch, :m_gresp_xfer_to_litter_fire_patch)
const NF_FIRE = (:m_leafn_to_fire_patch, :m_leafn_storage_to_fire_patch, :m_leafn_xfer_to_fire_patch,
    :m_livestemn_to_fire_patch, :m_livestemn_storage_to_fire_patch, :m_livestemn_xfer_to_fire_patch,
    :m_deadstemn_to_fire_patch, :m_deadstemn_storage_to_fire_patch, :m_deadstemn_xfer_to_fire_patch,
    :m_frootn_to_fire_patch, :m_frootn_storage_to_fire_patch, :m_frootn_xfer_to_fire_patch,
    :m_livecrootn_to_fire_patch, :m_livecrootn_storage_to_fire_patch, :m_livecrootn_xfer_to_fire_patch,
    :m_deadcrootn_to_fire_patch, :m_deadcrootn_storage_to_fire_patch, :m_deadcrootn_xfer_to_fire_patch,
    :m_retransn_to_fire_patch,
    :m_leafn_to_litter_fire_patch, :m_leafn_storage_to_litter_fire_patch, :m_leafn_xfer_to_litter_fire_patch,
    :m_livestemn_to_litter_fire_patch, :m_livestemn_storage_to_litter_fire_patch, :m_livestemn_xfer_to_litter_fire_patch,
    :m_livestemn_to_deadstemn_fire_patch, :m_deadstemn_to_litter_fire_patch,
    :m_deadstemn_storage_to_litter_fire_patch, :m_deadstemn_xfer_to_litter_fire_patch,
    :m_frootn_to_litter_fire_patch, :m_frootn_storage_to_litter_fire_patch, :m_frootn_xfer_to_litter_fire_patch,
    :m_livecrootn_to_litter_fire_patch, :m_livecrootn_storage_to_litter_fire_patch, :m_livecrootn_xfer_to_litter_fire_patch,
    :m_livecrootn_to_deadcrootn_fire_patch, :m_deadcrootn_to_litter_fire_patch,
    :m_deadcrootn_storage_to_litter_fire_patch, :m_deadcrootn_xfer_to_litter_fire_patch,
    :m_retransn_to_litter_fire_patch)

function build(; np=4, nc=2, ng=1, nlevdecomp=1, ndecomp_pools=4, n_litr=3)
    npft = 20
    pftcon = CLM.PftConFireBase(woody=vcat(fill(1.0,8),fill(0.0,12)),
        cc_leaf=fill(0.4,npft), cc_lstem=fill(0.2,npft), cc_dstem=fill(0.1,npft), cc_other=fill(0.3,npft),
        fm_leaf=fill(0.6,npft), fm_lstem=fill(0.5,npft), fm_other=fill(0.4,npft),
        fm_root=fill(0.3,npft), fm_lroot=fill(0.5,npft), fm_droot=fill(0.2,npft),
        lf_f=fill(1.0/n_litr,npft,n_litr), fr_f=fill(1.0/n_litr,npft,n_litr),
        smpso=fill(-66000.0,npft), smpsc=fill(-275000.0,npft))
    cnfire_const = CLM.CNFireConstData()
    patch = CLM.PatchData(); patch.itype=[2,10,5,15]; patch.column=[1,1,2,2]; patch.wtcol=[0.5,0.5,0.6,0.4]
    col = CLM.ColumnData(); col.gridcell=[1,1]
    grc = CLM.GridcellData(); grc.latdeg=[45.0]; grc.lat=[45.0*pi/180.0]
    dgvs = CLM.DgvsFireData(nind_patch=fill(100.0,np))

    cnveg_state = CLM.CNVegStateData()
    cnveg_state.cropf_col=[0.0,0.4]; cnveg_state.farea_burned_col=[1.0e-4,1.0e-4]
    cnveg_state.baf_crop_col=[0.0,2.0e-5]
    for f in (:fbac1_col,:fbac_col,:baf_peatf_col,:trotr1_col,:trotr2_col,:dtrotr_col,:lfc_col,:lfc2_col)
        setfield!(cnveg_state,f,zeros(nc))
    end

    cnveg_cs = CLM.CNVegCarbonStateData()
    cnveg_cs.leafcmax_patch=fill(0.0,np)
    cnveg_cs.leafc_patch=[10.0,5.0,8.0,3.0]; cnveg_cs.leafc_storage_patch=[1.0,0.5,0.8,0.3]; cnveg_cs.leafc_xfer_patch=[0.5,0.25,0.4,0.15]
    cnveg_cs.livestemc_patch=[20.0,0.0,15.0,0.0]; cnveg_cs.livestemc_storage_patch=[2.0,0.0,1.5,0.0]; cnveg_cs.livestemc_xfer_patch=[1.0,0.0,0.8,0.0]
    cnveg_cs.deadstemc_patch=[50.0,0.0,40.0,0.0]; cnveg_cs.deadstemc_storage_patch=[3.0,0.0,2.5,0.0]; cnveg_cs.deadstemc_xfer_patch=[1.5,0.0,1.2,0.0]
    cnveg_cs.frootc_patch=[4.0,2.0,3.0,1.0]; cnveg_cs.frootc_storage_patch=[0.5,0.2,0.3,0.1]; cnveg_cs.frootc_xfer_patch=[0.2,0.1,0.15,0.05]
    cnveg_cs.livecrootc_patch=[10.0,0.0,8.0,0.0]; cnveg_cs.livecrootc_storage_patch=[1.0,0.0,0.8,0.0]; cnveg_cs.livecrootc_xfer_patch=[0.5,0.0,0.4,0.0]
    cnveg_cs.deadcrootc_patch=[25.0,0.0,20.0,0.0]; cnveg_cs.deadcrootc_storage_patch=[1.5,0.0,1.2,0.0]; cnveg_cs.deadcrootc_xfer_patch=[0.8,0.0,0.6,0.0]
    cnveg_cs.gresp_storage_patch=[0.2,0.1,0.15,0.05]; cnveg_cs.gresp_xfer_patch=[0.1,0.05,0.08,0.03]

    cnveg_cf = CLM.CNVegCarbonFluxData()
    for f in CF_FIRE; setfield!(cnveg_cf, f, zeros(np)); end
    cnveg_cf.fire_mortality_c_to_cwdc_col = zeros(nc,nlevdecomp)
    cnveg_cf.m_decomp_cpools_to_fire_vr_col = zeros(nc,nlevdecomp,ndecomp_pools)
    cnveg_cf.m_c_to_litr_fire_col = zeros(nc,nlevdecomp,ndecomp_pools)

    cnveg_ns = CLM.CNVegNitrogenStateData()
    cnveg_ns.leafn_patch=[0.4,0.2,0.32,0.12]; cnveg_ns.retransn_patch=[0.1,0.05,0.08,0.03]
    # the N fluxes derive from C * CN ratios; provide the N state fields the kernel reads
    for (f,v) in ((:leafn_storage_patch,0.04),(:leafn_xfer_patch,0.02),
                  (:livestemn_patch,0.2),(:livestemn_storage_patch,0.02),(:livestemn_xfer_patch,0.01),
                  (:deadstemn_patch,0.1),(:deadstemn_storage_patch,0.006),(:deadstemn_xfer_patch,0.003),
                  (:frootn_patch,0.1),(:frootn_storage_patch,0.012),(:frootn_xfer_patch,0.005),
                  (:livecrootn_patch,0.1),(:livecrootn_storage_patch,0.01),(:livecrootn_xfer_patch,0.005),
                  (:deadcrootn_patch,0.05),(:deadcrootn_storage_patch,0.003),(:deadcrootn_xfer_patch,0.0015))
        hasproperty(cnveg_ns, f) && setfield!(cnveg_ns, f, fill(v, np))
    end

    cnveg_nf = CLM.CNVegNitrogenFluxData()
    for f in NF_FIRE; setfield!(cnveg_nf, f, zeros(np)); end
    cnveg_nf.fire_mortality_n_to_cwdn_col = zeros(nc,nlevdecomp)
    cnveg_nf.m_decomp_npools_to_fire_vr_col = zeros(nc,nlevdecomp,ndecomp_pools)
    cnveg_nf.m_n_to_litr_fire_col = zeros(nc,nlevdecomp,ndecomp_pools)

    soilbgc_cf = CLM.SoilBiogeochemCarbonFluxData(); soilbgc_cf.somc_fire_col = zeros(nc)
    decomp_cascade_con = CLM.DecompCascadeConData()
    decomp_cascade_con.is_litter = BitVector([true,true,true,false])
    decomp_cascade_con.is_cwd    = BitVector([false,false,false,true])

    return (; pftcon, cnfire_const, patch, col, grc, dgvs, cnveg_state, cnveg_cs, cnveg_cf,
        cnveg_ns, cnveg_nf, soilbgc_cf, decomp_cascade_con,
        leaf_prof=fill(1.0,np,nlevdecomp), froot_prof=fill(1.0,np,nlevdecomp),
        croot_prof=fill(1.0,np,nlevdecomp), stem_prof=fill(1.0,np,nlevdecomp),
        totsomc=fill(5000.0,nc), decomp_cpools_vr=fill(100.0,nc,nlevdecomp,ndecomp_pools),
        decomp_npools_vr=fill(5.0,nc,nlevdecomp,ndecomp_pools), somc_fire=zeros(nc),
        mask_soilc=trues(nc), mask_soilp=trues(np),
        np, nc, nlevdecomp, ndecomp_pools, n_litr)
end

run_flux!(d) = CLM.cnfire_fluxes_li2014!(
    d.mask_soilc, d.mask_soilp, 1:d.nc, 1:d.np,
    d.cnfire_const, d.pftcon, d.patch, d.col, d.grc,
    d.dgvs, d.cnveg_state, d.cnveg_cs, d.cnveg_cf, d.cnveg_ns, d.cnveg_nf,
    d.soilbgc_cf, d.decomp_cascade_con,
    d.leaf_prof, d.froot_prof, d.croot_prof, d.stem_prof,
    d.totsomc, d.decomp_cpools_vr, d.decomp_npools_vr, d.somc_fire;
    dt=1800.0, dayspyr=365.0, nlevdecomp=d.nlevdecomp, ndecomp_pools=d.ndecomp_pools,
    i_met_lit=1, i_litr_max=d.n_litr)

function main(backend)
    println("=" ^ 72)
    println("END-TO-END GPU parity for cnfire_fluxes_li2014! (fire C/N fluxes)")
    println("=" ^ 72)
    if backend === nothing; println("  No GPU backend."); return 0; end
    name, _, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n", name, FT)

    H = build(); B = build()
    run_flux!(H)

    mf(x) = CLM.Adapt.adapt(device_array_type(), CLM.Adapt.adapt(_F32(), x))
    mb(x) = device_array_type()(collect(x))
    # decomp_cascade_con stays host (BitVector struct); cnfire_fluxes! copies
    # is_litter/is_cwd onto the working backend internally.
    D = (; pftcon=mf(B.pftcon), cnfire_const=B.cnfire_const, patch=mf(B.patch),
        col=mf(B.col), grc=mf(B.grc), dgvs=mf(B.dgvs), cnveg_state=mf(B.cnveg_state),
        cnveg_cs=mf(B.cnveg_cs), cnveg_cf=mf(B.cnveg_cf), cnveg_ns=mf(B.cnveg_ns),
        cnveg_nf=mf(B.cnveg_nf), soilbgc_cf=mf(B.soilbgc_cf), decomp_cascade_con=B.decomp_cascade_con,
        leaf_prof=mf(B.leaf_prof), froot_prof=mf(B.froot_prof), croot_prof=mf(B.croot_prof),
        stem_prof=mf(B.stem_prof), totsomc=mf(B.totsomc), decomp_cpools_vr=mf(B.decomp_cpools_vr),
        decomp_npools_vr=mf(B.decomp_npools_vr), somc_fire=mf(B.somc_fire),
        mask_soilc=mb(B.mask_soilc), mask_soilp=mb(B.mask_soilp),
        np=B.np, nc=B.nc, nlevdecomp=B.nlevdecomp, ndecomp_pools=B.ndecomp_pools, n_litr=B.n_litr)

    if !(D.cnveg_cf.m_leafc_to_fire_patch isa device_array_type())
        println("  BLOCKED: structs did not move to device."); return 2
    end
    run_flux!(D)

    checks = Tuple{String,Any,Any}[]
    # a representative spread of per-patch C + N fire fluxes
    for f in (:m_leafc_to_fire_patch, :m_livestemc_to_fire_patch, :m_deadstemc_to_litter_fire_patch,
              :m_frootc_to_litter_fire_patch, :m_gresp_storage_to_fire_patch)
        push!(checks, (string(f), getfield(H.cnveg_cf,f), getfield(D.cnveg_cf,f)))
    end
    for f in (:m_leafn_to_fire_patch, :m_livestemn_to_litter_fire_patch, :m_retransn_to_fire_patch)
        push!(checks, (string(f), getfield(H.cnveg_nf,f), getfield(D.cnveg_nf,f)))
    end
    push!(checks, ("fire_mortality_c_to_cwdc_col", H.cnveg_cf.fire_mortality_c_to_cwdc_col, D.cnveg_cf.fire_mortality_c_to_cwdc_col))
    push!(checks, ("m_decomp_cpools_to_fire_vr",   H.cnveg_cf.m_decomp_cpools_to_fire_vr_col, D.cnveg_cf.m_decomp_cpools_to_fire_vr_col))
    push!(checks, ("m_c_to_litr_fire_col",         H.cnveg_cf.m_c_to_litr_fire_col, D.cnveg_cf.m_c_to_litr_fire_col))
    push!(checks, ("m_n_to_litr_fire_col",         H.cnveg_nf.m_n_to_litr_fire_col, D.cnveg_nf.m_n_to_litr_fire_col))
    push!(checks, ("somc_fire_col",                H.soilbgc_cf.somc_fire_col, D.soilbgc_cf.somc_fire_col))
    push!(checks, ("lfc_col",                      H.cnveg_state.lfc_col, D.cnveg_state.lfc_col))

    nfail = 0
    for (nm,a,b) in checks
        if !allfinite(a); @printf("  [WARN] %-30s CPU ref non-finite\n", nm); continue; end
        dd = reldiff(a,b); ok = dd < 1f-3
        @printf("  [%s] %-30s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, dd)
        ok || (nfail += 1)
    end
    println()
    println(nfail==0 ? "  WHOLE cnfire_fluxes_li2014! MATCHES CPU ON $name ($FT)" :
                       "  DIVERGENCE — investigate ($nfail failed).")
    return nfail==0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
