# ==========================================================================
# gpu_validate_fire_area_li2016_e2e.jl — end-to-end GPU parity for the WHOLE
# cnfire_area_li2016! burned-area driver (revised fire-occurrence formulation:
# afuel/arh/arh30 + bt btran term, revised GDP/pop/ag-fire terms, degrees
# lightning normalization, and the AD-spinup fuel branch).
#
# Builds the fire dataset (clone of test_fire_methods' make_fire_methods_area_data)
# at Float64, runs cnfire_area_li2016! on the CPU, adapts every kernel-touched
# struct to Metal/Float32, runs the SAME call on the device, and compares the
# mutated column burned-area diagnostics.
#
#   julia --project=scripts scripts/gpu_validate_fire_area_li2016_e2e.jl
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

function build(; np=4, nc=2, ng=1, nlevgrnd=3, nlevdecomp=1, ndecomp_pools=4, n_litr=3)
    npft = 20
    pftcon = CLM.PftConFireBase(
        woody = vcat(fill(1.0,8), fill(0.0,12)),
        cc_leaf=fill(0.4,npft), cc_lstem=fill(0.2,npft), cc_dstem=fill(0.1,npft), cc_other=fill(0.3,npft),
        fm_leaf=fill(0.6,npft), fm_lstem=fill(0.5,npft), fm_other=fill(0.4,npft),
        fm_root=fill(0.3,npft), fm_lroot=fill(0.5,npft), fm_droot=fill(0.2,npft),
        lf_f=fill(1.0/n_litr,npft,n_litr), fr_f=fill(1.0/n_litr,npft,n_litr),
        smpso=fill(-66000.0,npft), smpsc=fill(-275000.0,npft),
        rswf_min=fill(0.1,npft), rswf_max=fill(0.9,npft))
    pftcon_li2014 = CLM.PftConFireLi2014(fsr_pft=fill(0.2,npft), fd_pft=fill(1.0,npft))
    cnfire_const = CLM.CNFireConstData()
    cnfire_params = CLM.CNFireParams(prh30=0.05, ignition_efficiency=0.02)
    fire_data = CLM.CNFireBaseData(btran2_patch=zeros(np))
    fire_li2014 = CLM.CNFireLi2014Data(forc_hdm=[50.0], forc_lnfm=[0.05],
        gdp_lf_col=[10.0,10.0], peatf_lf_col=[0.0,0.0], abm_lf_col=[6,6])

    patch = CLM.PatchData(); patch.itype=[2,10,5,15]; patch.column=[1,1,2,2]
    patch.wtcol=[0.5,0.5,0.6,0.4]; patch.wtgcell=[0.5,0.5,0.6,0.4]
    col = CLM.ColumnData(); col.gridcell=[1,1]; col.wtgcell=[1.0,1.0]
    grc = CLM.GridcellData(); grc.latdeg=[45.0]; grc.lat=[45.0*pi/180.0]

    soilstate = CLM.SoilStateData()
    soilstate.watsat_col=fill(0.45,nc,nlevgrnd); soilstate.rootfr_patch=fill(1.0/nlevgrnd,np,nlevgrnd)
    soilstate.sucsat_col=fill(200.0,nc,nlevgrnd); soilstate.bsw_col=fill(5.0,nc,nlevgrnd)
    h2osoi_vol_col = fill(0.1,nc,nlevgrnd)

    cnveg_state = CLM.CNVegStateData()
    cnveg_state.dwt_smoothed_patch=zeros(np)
    for f in (:cropf_col,:baf_crop_col,:baf_peatf_col,:fbac_col,:fbac1_col,:farea_burned_col,
              :nfire_col,:fsr_col,:fd_col,:lgdp_col,:lgdp1_col,:lpop_col,:lfwt_col,
              :trotr1_col,:trotr2_col,:dtrotr_col,:lfc_col,:wtlf_col)
        setfield!(cnveg_state, f, zeros(nc))
    end
    cnveg_state.burndate_patch = fill(10000, np)

    cnveg_cs = CLM.CNVegCarbonStateData()
    cnveg_cs.totvegc_col=[500.0,400.0]
    for f in (:rootc_col,:leafc_col,:deadstemc_col,:fuelc_col,:fuelc_crop_col)
        setfield!(cnveg_cs,f,zeros(nc))
    end
    cnveg_cs.leafc_patch=[10.0,5.0,8.0,3.0]; cnveg_cs.leafc_storage_patch=[1.0,0.5,0.8,0.3]
    cnveg_cs.leafc_xfer_patch=[0.5,0.25,0.4,0.15]; cnveg_cs.frootc_patch=[4.0,2.0,3.0,1.0]
    cnveg_cs.frootc_storage_patch=[0.5,0.2,0.3,0.1]; cnveg_cs.frootc_xfer_patch=[0.2,0.1,0.15,0.05]
    cnveg_cs.deadcrootc_patch=[25.0,0.0,20.0,0.0]; cnveg_cs.deadcrootc_storage_patch=[1.5,0.0,1.2,0.0]
    cnveg_cs.deadcrootc_xfer_patch=[0.8,0.0,0.6,0.0]; cnveg_cs.livecrootc_patch=[10.0,0.0,8.0,0.0]
    cnveg_cs.livecrootc_storage_patch=[1.0,0.0,0.8,0.0]; cnveg_cs.livecrootc_xfer_patch=[0.5,0.0,0.4,0.0]
    cnveg_cs.deadstemc_patch=[50.0,0.0,40.0,0.0]

    decomp_cascade_con = CLM.DecompCascadeConData()
    decomp_cascade_con.is_litter = BitVector([true,true,true,false])
    decomp_cascade_con.is_cwd    = BitVector([false,false,false,true])
    decomp_cascade_con.spinup_factor = ones(ndecomp_pools)

    return (; pftcon, pftcon_li2014, cnfire_const, cnfire_params, fire_data, fire_li2014,
        patch, col, grc, soilstate, h2osoi_vol_col, cnveg_state, cnveg_cs, decomp_cascade_con,
        totlitc_col=[200.0,150.0], decomp_cpools_vr_col=fill(100.0,nc,nlevdecomp,ndecomp_pools),
        t_soi17cm_col=fill(280.0,nc),
        forc_rh_grc=[50.0], forc_wind_grc=[3.0], forc_t_col=fill(290.0,nc),
        forc_rain_col=fill(0.0,nc), forc_snow_col=fill(0.0,nc),
        prec60_patch=fill(2.0e-5,np), prec10_patch=fill(3.0e-5,np),
        rh30_patch=fill(60.0,np),
        fsat_col=fill(0.1,nc), wf2_col=fill(0.25,nc),
        mask_soilc=trues(nc), mask_soilp=trues(np),
        mask_exposedveg=trues(np), mask_noexposedveg=falses(np),
        np, nc, ng, nlevgrnd, nlevdecomp, ndecomp_pools)
end

run_area!(d) = CLM.cnfire_area_li2016!(
    d.fire_li2014, d.pftcon_li2014, d.fire_data, d.cnfire_const, d.cnfire_params, d.pftcon,
    d.mask_soilc, d.mask_soilp, d.mask_exposedveg, d.mask_noexposedveg,
    1:d.nc, 1:d.np, d.patch, d.col, d.grc, d.soilstate, d.h2osoi_vol_col,
    d.cnveg_state, d.cnveg_cs, d.decomp_cascade_con,
    d.totlitc_col, d.decomp_cpools_vr_col, d.t_soi17cm_col;
    forc_rh_grc=d.forc_rh_grc, forc_wind_grc=d.forc_wind_grc, forc_t_col=d.forc_t_col,
    forc_rain_col=d.forc_rain_col, forc_snow_col=d.forc_snow_col,
    prec60_patch=d.prec60_patch, prec10_patch=d.prec10_patch, rh30_patch=d.rh30_patch,
    fsat_col=d.fsat_col, wf2_col=d.wf2_col,
    dt=1800.0, dayspyr=365.0, kmo=6, kda=15, mcsec=3600, nstep=10,
    nlevgrnd=d.nlevgrnd, nlevdecomp=d.nlevdecomp, ndecomp_pools=d.ndecomp_pools)

const OUT_FIELDS = (:farea_burned_col, :nfire_col, :fbac_col, :fbac1_col, :baf_crop_col,
                    :baf_peatf_col, :fsr_col, :fd_col, :lgdp_col, :lgdp1_col, :lpop_col,
                    :trotr1_col, :lfwt_col, :wtlf_col, :cropf_col)

function main(backend)
    println("=" ^ 72)
    println("END-TO-END GPU parity for cnfire_area_li2016! (Li2016 burned area)")
    println("=" ^ 72)
    if backend === nothing
        println("  No GPU backend — nothing to validate."); return 0
    end
    name, _, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n", name, FT)

    CLM.dzsoi_decomp[] = [0.1]
    H = build()
    B = build()

    run_area!(H)   # CPU reference

    mf(x) = CLM.Adapt.adapt(device_array_type(), CLM.Adapt.adapt(_F32(), x))
    mb(x) = device_array_type()(collect(x))

    D = (; fire_li2014=mf(B.fire_li2014), pftcon_li2014=mf(B.pftcon_li2014),
        fire_data=mf(B.fire_data), cnfire_const=B.cnfire_const, cnfire_params=B.cnfire_params,
        pftcon=mf(B.pftcon), patch=mf(B.patch), col=mf(B.col), grc=mf(B.grc),
        soilstate=mf(B.soilstate), h2osoi_vol_col=mf(B.h2osoi_vol_col),
        cnveg_state=mf(B.cnveg_state), cnveg_cs=mf(B.cnveg_cs), decomp_cascade_con=B.decomp_cascade_con,
        totlitc_col=mf(B.totlitc_col), decomp_cpools_vr_col=mf(B.decomp_cpools_vr_col),
        t_soi17cm_col=mf(B.t_soi17cm_col),
        forc_rh_grc=mf(B.forc_rh_grc), forc_wind_grc=mf(B.forc_wind_grc), forc_t_col=mf(B.forc_t_col),
        forc_rain_col=mf(B.forc_rain_col), forc_snow_col=mf(B.forc_snow_col),
        prec60_patch=mf(B.prec60_patch), prec10_patch=mf(B.prec10_patch), rh30_patch=mf(B.rh30_patch),
        fsat_col=mf(B.fsat_col), wf2_col=mf(B.wf2_col),
        mask_soilc=mb(B.mask_soilc), mask_soilp=mb(B.mask_soilp),
        mask_exposedveg=mb(B.mask_exposedveg), mask_noexposedveg=mb(B.mask_noexposedveg),
        np=B.np, nc=B.nc, ng=B.ng, nlevgrnd=B.nlevgrnd, nlevdecomp=B.nlevdecomp,
        ndecomp_pools=B.ndecomp_pools)

    if !(D.cnveg_state.farea_burned_col isa device_array_type())
        println("  BLOCKED: cnveg_state did not move to the device."); return 2
    end

    run_area!(D)

    nfail = 0
    for f in OUT_FIELDS
        a = getfield(H.cnveg_state, f); b = getfield(D.cnveg_state, f)
        if !allfinite(a)
            @printf("  [WARN] %-20s CPU ref non-finite — skipping\n", f); continue
        end
        dd = reldiff(a, b); ok = dd < 1f-3
        @printf("  [%s] %-20s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", f, dd)
        ok || (nfail += 1)
    end
    println()
    println(nfail == 0 ? "  WHOLE cnfire_area_li2016! MATCHES CPU ON $name ($FT)" :
                         "  DIVERGENCE — investigate ($nfail failed).")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
