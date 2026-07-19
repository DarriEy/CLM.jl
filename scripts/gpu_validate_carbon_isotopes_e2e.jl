# ==========================================================================
# gpu_validate_carbon_isotopes_e2e.jl — GPU parity for the C13/C14 tracer
# routines: c14_decay! (radioactive decay of all C14 pools — gridcell seedc,
# soil decomp vr pools incl. the spinup-accelerated branch, and the veg patch
# pools) and c13_c14_photosynthesis! (C13 fractionation + C14 bomb-factor flux).
#
# Mirrors test_carbon_isotopes' fixtures. Runs each on CPU, adapts the C14 state
# / photosynthesis struct to Metal/Float32, runs the SAME call on-device, and
# compares. (The C14 flux structs are inert on the default non-matrix path, so
# they stay host.)
#
#   julia --project=scripts scripts/gpu_validate_carbon_isotopes_e2e.jl
# ==========================================================================

using CLM
using Printf
include(joinpath(@__DIR__, "gpu_backends.jl"))

struct _F32 end
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:AbstractFloat}) = Float32.(x)
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:Integer}) = x
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{Bool}) = x
# leave scalar Float64 fields as-is (avoids breaking Adapt reconstruction of
# structs with concrete ::Float64 fields, e.g. SoilBiogeochemCarbonStateData).

reldiff(a, b) = begin
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(Float64(A[i]) - Float64(B[i])) / (1.0 + max(abs(Float64(A[i])), abs(Float64(B[i])))))
    end
    m
end

mf(x) = CLM.Adapt.adapt(device_array_type(), CLM.Adapt.adapt(_F32(), x))
mb(x) = device_array_type()(collect(x))
mi(x) = device_array_type()(collect(Int32.(x)))

function make_c14(nc, np, ng, nlevdecomp, ndecomp_pools)
    cs = CLM.CNVegCarbonStateData(); CLM.cnveg_carbon_state_init!(cs, np, nc, ng)
    cs.leafc_patch .= [100.0, 200.0, 300.0]; cs.cpool_patch .= [50.0, 60.0, 70.0]
    cs.xsmrpool_patch .= [10.0, 20.0, 30.0]; cs.deadstemc_patch .= [1000.0, 2000.0, 3000.0]
    cs.gresp_storage_patch .= [5.0, 10.0, 15.0]; cs.gresp_xfer_patch .= [1.0, 2.0, 3.0]
    cs.ctrunc_patch .= [0.1, 0.2, 0.3]; cs.seedc_grc .= [500.0]
    scs = CLM.SoilBiogeochemCarbonStateData()
    CLM.soil_bgc_carbon_state_init!(scs, nc, ng, nlevdecomp, ndecomp_pools)
    scs.decomp_cpools_vr_col .= 100.0
    return cs, scs
end

function make_ps(np)
    ps = CLM.PhotosynthesisData(); CLM.photosynthesis_data_init!(ps, np; nlevcan=1)
    ps.psnsun_patch .= [5.0, 10.0, 0.0]; ps.psnsha_patch .= [3.0, 6.0, 0.0]
    ps.alphapsnsun_patch .= [1.02, 1.03, 1.0]; ps.alphapsnsha_patch .= [1.01, 1.02, 1.0]
    return ps
end

function main(backend)
    println("=" ^ 72)
    println("GPU parity for the C13/C14 tracer routines (c14_decay! + photosynthesis)")
    println("=" ^ 72)
    if backend === nothing
        println("  No GPU backend — nothing to validate."); return 0
    end
    name, _, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n", name, FT)
    nc, np, ng, nlev, ndp = 2, 3, 1, 2, 3
    nfail = 0
    report(nm, a, b) = begin
        dd = reldiff(a, b); ok = dd < 1f-3
        @printf("  [%s] %-34s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, dd)
        ok || (global nfail += 1)
    end

    # --- c14_decay! basic (spinup_state=0) ---
    csH, scsH = make_c14(nc, np, ng, nlev, ndp)
    cfH = CLM.CNVegCarbonFluxData(); CLM.cnveg_carbon_flux_init!(cfH, np, nc, ng)
    scfH = CLM.SoilBiogeochemCarbonFluxData(); CLM.soil_bgc_carbon_flux_init!(scfH, nc, nlev, ndp, 3)
    CLM.c14_decay!(csH, cfH, scsH, scfH; mask_soilc=trues(nc), mask_soilp=trues(np),
        bounds_col=1:nc, bounds_patch=1:np, bounds_gridcell=1:ng, dt=1800.0,
        nlevdecomp=nlev, ndecomp_pools=ndp)

    csD, scsD = make_c14(nc, np, ng, nlev, ndp); csD = mf(csD); scsD = mf(scsD)
    if !(csD.leafc_patch isa device_array_type())
        println("  BLOCKED: C14 state did not move to the device."); return 2
    end
    CLM.c14_decay!(csD, cfH, scsD, scfH; mask_soilc=mb(trues(nc)), mask_soilp=mb(trues(np)),
        bounds_col=1:nc, bounds_patch=1:np, bounds_gridcell=1:ng, dt=1800.0,
        nlevdecomp=nlev, ndecomp_pools=ndp)
    report("decay basic: seedc_grc", csH.seedc_grc, csD.seedc_grc)
    report("decay basic: leafc_patch", csH.leafc_patch, csD.leafc_patch)
    report("decay basic: deadstemc_patch", csH.deadstemc_patch, csD.deadstemc_patch)
    report("decay basic: decomp_cpools_vr", scsH.decomp_cpools_vr_col, scsD.decomp_cpools_vr_col)

    # --- c14_decay! spinup (spinup_state=1, spinup_factor != 1 → latitude term) ---
    spf = [1.0, 8.0, 20.0]; colg = [1, 1]; latd = [55.0]
    csH2, scsH2 = make_c14(nc, np, ng, nlev, ndp)
    CLM.c14_decay!(csH2, cfH, scsH2, scfH; mask_soilc=trues(nc), mask_soilp=trues(np),
        bounds_col=1:nc, bounds_patch=1:np, bounds_gridcell=1:ng, dt=1800.0,
        nlevdecomp=nlev, ndecomp_pools=ndp, spinup_state=1, spinup_factor=spf,
        col_gridcell=colg, latdeg_grc=latd)
    csD2, scsD2 = make_c14(nc, np, ng, nlev, ndp); csD2 = mf(csD2); scsD2 = mf(scsD2)
    CLM.c14_decay!(csD2, cfH, scsD2, scfH; mask_soilc=mb(trues(nc)), mask_soilp=mb(trues(np)),
        bounds_col=1:nc, bounds_patch=1:np, bounds_gridcell=1:ng, dt=1800.0,
        nlevdecomp=nlev, ndecomp_pools=ndp, spinup_state=1, spinup_factor=mf(spf),
        col_gridcell=mi(colg), latdeg_grc=mf(latd))
    report("decay spinup: decomp_cpools_vr", scsH2.decomp_cpools_vr_col, scsD2.decomp_cpools_vr_col)

    # --- c13_c14_photosynthesis! (use_c13 + use_c14) ---
    psH = make_ps(np)
    CLM.c13_c14_photosynthesis!(psH, [400.0e-6], [4.4e-6], [1,1,1], [51.0], trues(np), 1:np;
        use_c13=true, use_c14=true)
    psD = mf(make_ps(np))
    CLM.c13_c14_photosynthesis!(psD, mf([400.0e-6]), mf([4.4e-6]), mi([1,1,1]), mf([51.0]),
        mb(trues(np)), 1:np; use_c13=true, use_c14=true)
    report("photo: rc13_canair", psH.rc13_canair_patch, psD.rc13_canair_patch)
    report("photo: c13_psnsun", psH.c13_psnsun_patch, psD.c13_psnsun_patch)
    report("photo: c13_psnsha", psH.c13_psnsha_patch, psD.c13_psnsha_patch)
    report("photo: c14_psnsun", psH.c14_psnsun_patch, psD.c14_psnsun_patch)

    println()
    println(nfail == 0 ? "  C13/C14 tracer routines MATCH CPU ON $name ($FT)" :
                         "  DIVERGENCE — investigate ($nfail failed).")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
