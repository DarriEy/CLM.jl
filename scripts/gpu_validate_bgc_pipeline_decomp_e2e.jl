# ==========================================================================
# gpu_validate_bgc_pipeline_decomp_e2e.jl — WHOLE-DRIVER GPU parity for the
# FULL-DECOMP branch of cn_driver_no_leaching! (_has_decomp = true). Drives the
# real BGC decomposition chain through the orchestrator on the CPU and on Metal:
#   decomp_rate_constants_bgc! -> soil_bgc_potential! -> soil_biogeochem_decomp!
#   -> C/N state-update cascade -> precision control -> summarize
# over one synthetic CN+soil state, and compares the full resulting state.
# cascade_con + param structs stay host (the modules copy needed fields onto the
# backend); the state structs + decomp_bgc_state + forcing move to the device.
#
#   julia --project=scripts scripts/gpu_validate_bgc_pipeline_decomp_e2e.jl
# ==========================================================================
using CLM
using Printf
include(joinpath(@__DIR__, "gpu_backends.jl"))

# _F32: arrays only (the soil/cnveg state structs pin concrete ::Float64 scalars
# that must stay Float64). _F32S also converts scalar floats (DecompBGCState mixes
# ::FT scalars with its arrays under one {FT,M} param).
struct _F32 end
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:AbstractFloat}) = Float32.(x)
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:Integer}) = x
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{Bool}) = x
struct _F32S end
CLM.Adapt.adapt_storage(::_F32S, x::AbstractArray{<:AbstractFloat}) = Float32.(x)
CLM.Adapt.adapt_storage(::_F32S, x::AbstractArray{<:Integer}) = x
CLM.Adapt.adapt_storage(::_F32S, x::AbstractArray{Bool}) = x
CLM.Adapt.adapt_storage(::_F32S, x::AbstractFloat) = Float32(x)

function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(Float64(A[i]) - Float64(B[i])) / (1.0 + max(abs(Float64(A[i])), abs(Float64(B[i])))))
    end
    return m
end
finite_ok(a) = any(isfinite, Array(a))

make_bgc_params() = CLM.DecompBGCParams(
    cn_s1_bgc=12.0, cn_s2_bgc=12.0, cn_s3_bgc=10.0,
    rf_l1s1_bgc=0.39, rf_l2s1_bgc=0.55, rf_l3s2_bgc=0.29,
    rf_s2s1_bgc=0.55, rf_s2s3_bgc=0.55, rf_s3s1_bgc=0.55, rf_cwdl3_bgc=0.0,
    tau_l1_bgc=1.0/18.5, tau_l2_l3_bgc=1.0/4.9, tau_s1_bgc=1.0/7.3,
    tau_s2_bgc=1.0/0.2, tau_s3_bgc=1.0/0.0045, cwd_fcel_bgc=0.45,
    bgc_initial_Cstocks=fill(200.0, 7), bgc_initial_Cstocks_depth=0.3)
make_cn_params() = CLM.CNSharedParamsData(Q10=1.5, minpsi=-10.0, maxpsi=-0.1,
    rf_cwdl2=0.0, tau_cwd=10.0, cwd_flig=0.24, froz_q10=1.5,
    decomp_depth_efolding=0.5, mino2lim=0.0)

function make_data(; nlevdecomp=5)
    # nlevdecomp=5 → multi-level decomp (reads per-layer t_soisno/soilpsi directly).
    # nlevdecomp=1 → single-level decomp, which needs col_dz — the driver builds it
    # from a ColumnData (col.dz slice), so we supply a col below.
    nc=4; np=6; ng=2; ndecomp_pools=7; ncascade_max=10; nrepr=1
    nlevfull=5
    i_litr_min=1; i_litr_max=3; i_cwd=4; dt=1800.0
    mask_soilc = trues(nc); mask_soilp = trues(np)

    # --- decomposition infrastructure (real BGC cascade) ---
    bgc_params = make_bgc_params(); cn_params = make_cn_params()
    bgc_state = CLM.DecompBGCState()
    bgc_state.use_century_tfunc = false; bgc_state.normalize_q10_to_century_tfunc = true
    cascade_con = CLM.DecompCascadeConData()
    CLM.init_decomp_cascade_bgc!(bgc_state, cascade_con, bgc_params, cn_params;
        cellsand=fill(50.0, nc, nlevfull), bounds=1:nc, nlevdecomp=nlevdecomp,
        ndecomp_pools_max=ndecomp_pools, ndecomp_cascade_transitions_max=ncascade_max, use_fates=false)
    ndct = length(cascade_con.cascade_donor_pool)          # real # transitions

    decomp_params = CLM.DecompParams(dnp=0.01)

    # --- CN veg state/flux ---
    cs_veg = CLM.CNVegCarbonStateData(); CLM.cnveg_carbon_state_init!(cs_veg, np, nc, ng; nrepr=nrepr)
    CLM.cnveg_carbon_state_set_values!(cs_veg, mask_soilp, 10.0, mask_soilc, 10.0; nrepr=nrepr)
    cs_veg.cpool_patch .= 100.0; cs_veg.leafc_patch .= 50.0; cs_veg.frootc_patch .= 40.0
    cf_veg = CLM.CNVegCarbonFluxData()
    CLM.cnveg_carbon_flux_init!(cf_veg, np, nc, ng; nrepr=nrepr, nlevdecomp_full=nlevdecomp, ndecomp_pools=ndecomp_pools)
    CLM.cnveg_carbon_flux_set_values!(cf_veg, mask_soilp, 0.0, mask_soilc, 0.0; nrepr=nrepr, nlevdecomp_full=nlevdecomp, ndecomp_pools=ndecomp_pools)
    cf_veg.cpool_to_leafc_patch .= 1.0e-3; cf_veg.leafc_to_litter_patch .= 2.0e-4
    ns_veg = CLM.CNVegNitrogenStateData(); CLM.cnveg_nitrogen_state_init!(ns_veg, np, nc, ng; nrepr=nrepr)
    CLM.cnveg_nitrogen_state_set_values!(ns_veg, mask_soilp, 5.0, mask_soilc, 5.0; nrepr=nrepr)
    nf_veg = CLM.CNVegNitrogenFluxData()
    CLM.cnveg_nitrogen_flux_init!(nf_veg, np, nc, ng; nrepr=nrepr, nlevdecomp_full=nlevdecomp, ndecomp_pools=ndecomp_pools, i_litr_max=i_litr_max)
    CLM.cnveg_nitrogen_flux_set_values!(nf_veg, mask_soilp, 0.0, mask_soilc, 0.0; nrepr=nrepr, nlevdecomp_full=nlevdecomp, ndecomp_pools=ndecomp_pools, i_litr_max=i_litr_max)

    # --- soil BGC state/flux (cascade dims = ndct) ---
    cs_soil = CLM.SoilBiogeochemCarbonStateData(); CLM.soil_bgc_carbon_state_init!(cs_soil, nc, ng, nlevdecomp, ndecomp_pools)
    CLM.soil_bgc_carbon_state_set_values!(cs_soil, mask_soilc, 100.0)
    cf_soil = CLM.SoilBiogeochemCarbonFluxData(); CLM.soil_bgc_carbon_flux_init!(cf_soil, nc, nlevdecomp, ndecomp_pools, ndct)
    CLM.soil_bgc_carbon_flux_set_values!(cf_soil, mask_soilc, 0.0)
    ns_soil = CLM.SoilBiogeochemNitrogenStateData(); CLM.soil_bgc_nitrogen_state_init!(ns_soil, nc, ng, nlevdecomp, ndecomp_pools)
    CLM.soil_bgc_nitrogen_state_set_values!(ns_soil, mask_soilc, 5.0)
    nf_soil = CLM.SoilBiogeochemNitrogenFluxData(); CLM.soil_bgc_nitrogen_flux_init!(nf_soil, nc, nlevdecomp, ndecomp_pools, ndct)
    CLM.soil_bgc_nitrogen_flux_set_values!(nf_soil, mask_soilc, 0.0)
    soilbgc_st = CLM.SoilBiogeochemStateData(); CLM.soil_bgc_state_init!(soilbgc_st, nc, np, nlevdecomp, ndct)

    # single-level decomp needs a ColumnData so the driver can build col_dz from
    # col.dz[:, nlevsno+1:end]; only col.dz is read on the decomp path.
    col = nothing
    if nlevdecomp == 1
        CLM.varpar_init!(CLM.varpar, 17, 17, 0, 5)
        nls = CLM.varpar.nlevsno
        col = CLM.ColumnData()
        col.dz = fill(0.1, nc, nls + 10)   # slice [:, nls+1:end] = (nc,10), first 5 used
    end

    return (; config=CLM.CNDriverConfig(), bgc_params, cn_params, bgc_state, cascade_con, decomp_params, col,
            cs_veg, cf_veg, ns_veg, nf_veg, cs_soil, cf_soil, ns_soil, nf_soil, soilbgc_st,
            cascade_donor_pool=copy(cascade_con.cascade_donor_pool),
            cascade_receiver_pool=copy(cascade_con.cascade_receiver_pool),
            mask_soilc, mask_soilp, patch_column=[1,1,2,2,3,4], ivt=fill(1,np),
            woody=(w=zeros(80); w[2]=1.0; w), harvdate=fill(999,np), col_is_fates=fill(false,nc),
            t_soisno=fill(CLM.TFRZ+15.0, nc, nlevfull), soilpsi=fill(-1.0, nc, nlevfull),
            dzsoi_decomp=fill(0.1, max(nlevdecomp,10)), zsoi_vals=[0.01,0.04,0.09,0.16,0.26,0.40,0.58,0.80,1.06,1.36],
            nc, np, ng, dt, nlevdecomp, ndecomp_pools, ndct, nrepr, i_litr_min, i_litr_max, i_cwd)
end

function run_driver!(d)
    CLM.cn_driver_no_leaching!(d.config;
        mask_bgc_soilc=d.mask_soilc, mask_bgc_vegp=d.mask_soilp,
        bounds_col=1:d.nc, bounds_patch=1:d.np,
        nlevdecomp=d.nlevdecomp, ndecomp_pools=d.ndecomp_pools,
        ndecomp_cascade_transitions=d.ndct, i_litr_min=d.i_litr_min, i_litr_max=d.i_litr_max, i_cwd=d.i_cwd,
        patch_column=d.patch_column, ivt=d.ivt, woody=d.woody, harvdate=d.harvdate,
        col_is_fates=d.col_is_fates, cascade_donor_pool=d.cascade_donor_pool,
        cascade_receiver_pool=d.cascade_receiver_pool, dt=d.dt,
        cnveg_cs=d.cs_veg, cnveg_cf=d.cf_veg, cnveg_ns=d.ns_veg, cnveg_nf=d.nf_veg,
        soilbgc_cs=d.cs_soil, soilbgc_cf=d.cf_soil, soilbgc_ns=d.ns_soil, soilbgc_nf=d.nf_soil,
        soilbgc_state=d.soilbgc_st,
        # --- full-decomp infrastructure (triggers _has_decomp) ---
        cascade_con=d.cascade_con, decomp_bgc_state=d.bgc_state, decomp_bgc_params=d.bgc_params,
        cn_shared_params=d.cn_params, decomp_params=d.decomp_params, col=d.col,
        t_soisno=d.t_soisno, soilpsi=d.soilpsi, dzsoi_decomp=d.dzsoi_decomp, zsoi_vals=d.zsoi_vals)
end

function scenario(backend, nlevdecomp)
    name, _, FT = backend
    @printf("\n  -- scenario: nlevdecomp=%d (%s) --\n", nlevdecomp,
            nlevdecomp == 1 ? "single-level, col_dz path" : "multi-level")
    H = make_data(; nlevdecomp); B = make_data(; nlevdecomp)
    run_driver!(H)
    mf(x) = CLM.Adapt.adapt(device_array_type(), CLM.Adapt.adapt(_F32(), x))
    mfS(x) = CLM.Adapt.adapt(device_array_type(), CLM.Adapt.adapt(_F32S(), x))
    mm(b) = device_array_type()(collect(b))
    # state structs + decomp_bgc_state + forcing + col -> device; cascade_con + params stay host
    D = (; config=B.config, bgc_params=B.bgc_params, cn_params=B.cn_params,
         bgc_state=mfS(B.bgc_state), cascade_con=B.cascade_con, decomp_params=B.decomp_params,
         col=(B.col === nothing ? nothing : mf(B.col)),
         cs_veg=mf(B.cs_veg), cf_veg=mf(B.cf_veg), ns_veg=mf(B.ns_veg), nf_veg=mf(B.nf_veg),
         cs_soil=mf(B.cs_soil), cf_soil=mf(B.cf_soil), ns_soil=mf(B.ns_soil), nf_soil=mf(B.nf_soil),
         soilbgc_st=mf(B.soilbgc_st),
         cascade_donor_pool=mf(B.cascade_donor_pool), cascade_receiver_pool=mf(B.cascade_receiver_pool),
         mask_soilc=mm(B.mask_soilc), mask_soilp=mm(B.mask_soilp), patch_column=mf(B.patch_column),
         ivt=mf(B.ivt), woody=mf(B.woody), harvdate=mf(B.harvdate), col_is_fates=mf(B.col_is_fates),
         t_soisno=mf(B.t_soisno), soilpsi=mf(B.soilpsi), dzsoi_decomp=mf(B.dzsoi_decomp), zsoi_vals=mf(B.zsoi_vals),
         nc=B.nc, np=B.np, ng=B.ng, dt=B.dt, nlevdecomp=B.nlevdecomp, ndecomp_pools=B.ndecomp_pools,
         ndct=B.ndct, nrepr=B.nrepr, i_litr_min=B.i_litr_min, i_litr_max=B.i_litr_max, i_cwd=B.i_cwd)
    if !(D.cf_soil.decomp_k_col isa device_array_type()); println("  BLOCKED."); return 1; end
    run_driver!(D)
    checks = [
        ("cf_soil.decomp_k (rates)", H.cf_soil.decomp_k_col, D.cf_soil.decomp_k_col),
        ("cf_soil.t_scalar", H.cf_soil.t_scalar_col, D.cf_soil.t_scalar_col),
        ("cf_soil.w_scalar", H.cf_soil.w_scalar_col, D.cf_soil.w_scalar_col),
        ("cf_soil.hr_vr (potential)", H.cf_soil.decomp_cascade_hr_vr_col, D.cf_soil.decomp_cascade_hr_vr_col),
        ("cf_soil.ctransfer_vr", H.cf_soil.decomp_cascade_ctransfer_vr_col, D.cf_soil.decomp_cascade_ctransfer_vr_col),
        ("nf_soil.sminn_to_denit_vr", H.nf_soil.sminn_to_denit_decomp_cascade_vr_col, D.nf_soil.sminn_to_denit_decomp_cascade_vr_col),
        ("nf_soil.gross_nmin_vr", H.nf_soil.gross_nmin_vr_col, D.nf_soil.gross_nmin_vr_col),
        ("cs_soil.decomp_cpools_vr", H.cs_soil.decomp_cpools_vr_col, D.cs_soil.decomp_cpools_vr_col),
        ("cs_veg.cpool", H.cs_veg.cpool_patch, D.cs_veg.cpool_patch),
        ("cs_veg.leafc", H.cs_veg.leafc_patch, D.cs_veg.leafc_patch),
        ("ns_veg.leafn", H.ns_veg.leafn_patch, D.ns_veg.leafn_patch)]
    nfail = 0
    for (nm, a, b) in checks
        if !finite_ok(a); @printf("    [WARN] %-28s CPU all-NaN, skip\n", nm); continue; end
        dd = reldiff(a, b); ok = dd < 1f-3
        @printf("    [%s] %-28s rel = %.3e\n", ok ? "PASS" : "FAIL", nm, dd); ok || (nfail += 1)
    end
    return nfail
end

function main(backend)
    println("="^66); println("WHOLE-DRIVER GPU parity — FULL-DECOMP branch (_has_decomp)"); println("="^66)
    if backend === nothing; println("  No GPU backend."); return 0; end
    name, _, FT = backend
    @printf("  Backend: %s   (precision: %s)\n", name, FT)
    nfail = scenario(backend, 5) + scenario(backend, 1)
    println(); println(nfail == 0 ? "  FULL-DECOMP BGC DRIVER (single+multi level) MATCHES CPU ON $name ($FT)" :
                       "  DIVERGENCE ($nfail).")
    return nfail == 0 ? 0 : 1
end
const BACKEND = detect_backend()
exit(main(BACKEND))
