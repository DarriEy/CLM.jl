# ==========================================================================
# gpu_validate_luna.jl — host-vs-Metal parity for the LUNA 24hr/240hr climate
# accumulation kernels (clear24_climate_luna! / acc24_climate_luna! /
# acc240_climate_luna!). These run every timestep / end-of-day under use_luna.
#
# (The LUNA optimization core update_photosynthesis_capacity! is NOT covered here:
# its deep scalar chain — nitrogen_allocation!→nitrogen_investments!→
# photosynthesis_luna! — carries ~250 bare Float64 module constants across ~320
# lines of iterative optimization, which would need a full T-generic rewrite to
# compile on Metal-Float32. Left as a scoped follow-up.)
#
#   julia --project=scripts scripts/gpu_validate_luna.jl
# ==========================================================================
using CLM
using Printf
import Metal
include(joinpath(@__DIR__, "gpu_backends.jl"))
include(joinpath(@__DIR__, "gpu_adapt.jl"))
mf(x) = mf(Metal.MtlArray, x)

const NP, NC, NL, NG = 2, 5, 2, 1

function build_fixture()
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
    cs = CLM.CanopyStateData();  CLM.canopystate_init!(cs, NP)
    ps = CLM.PhotosynthesisData(); CLM.photosynthesis_data_init!(ps, NP; use_luna=true)
    alb = CLM.SurfaceAlbedoData(); CLM.surfalb_init!(alb, NP, NC, NG)
    sa = CLM.SolarAbsorbedData(); CLM.solarabs_init!(sa, NP, NL; use_luna=true)
    temp = CLM.TemperatureData(); CLM.temperature_init!(temp, NP, NC, NL, NG)
    pch = CLM.PatchData(); CLM.patch_init!(pch, NP)
    wdb = CLM.WaterDiagnosticBulkData(); CLM.waterdiagnosticbulk_init!(wdb, NC, NP, NL, NG)
    fv = CLM.FrictionVelocityData(); CLM.frictionvel_init!(fv, NP, NC)
    # populate a valid (non-first-day) daytime climate state
    temp.t_veg_day_patch .= 100.0;  temp.t_veg_night_patch .= 50.0
    temp.t_veg_patch .= 300.0;      temp.t_a10_patch .= 295.0
    temp.ndaysteps_patch .= 3;      temp.nnightsteps_patch .= 2
    temp.t_veg10_day_patch .= 298.0; temp.t_veg10_night_patch .= 290.0
    sa.sabv_patch .= 100.0;         alb.nrad_patch .= size(sa.par24d_z_patch, 2)
    cs.laisun_z_patch .= 1.0;       cs.laisha_z_patch .= 0.5
    sa.parsun_z_patch .= 200.0;     sa.parsha_z_patch .= 50.0
    sa.par24d_z_patch .= 30.0;      sa.par24x_z_patch .= 60.0
    sa.par240d_z_patch .= 25.0;     sa.par240x_z_patch .= 55.0
    ps.fpsn_patch .= 10.0;          ps.fpsn24_patch .= 5.0
    wdb.rh10_af_patch .= 0.6;       fv.rb10_patch .= 40.0
    return cs, ps, alb, sa, temp, pch, wdb, fv
end

reldiff(H, D) = begin
    A = Array(H); B = Array(D); m = 0.0; n = 0
    for i in eachindex(A, B)
        (isfinite(A[i]) && isfinite(B[i])) || continue
        m = max(m, abs(Float64(A[i]) - Float64(B[i])) / (1.0 + abs(Float64(A[i])))); n += 1
    end
    (m, n)
end

function main(backend)
    println("="^64); println("LUNA climate-accumulation kernels — host vs Metal"); println("="^64)
    backend === nothing && (println("  no GPU backend"); return 0)
    name, _, FT = backend
    @printf("  Backend: %s (%s)\n\n", name, FT)
    dtime = 1800.0; mask = trues(NP); bounds = 1:NP
    nfail = 0

    # ---- acc24 ----
    cs,ps,alb,sa,temp,pch,wdb,fv = build_fixture()
    csD,psD,albD,saD,tempD,pchD = mf(cs),mf(ps),mf(alb),mf(sa),mf(temp),mf(pch)
    maskD = mf(mask)
    CLM.acc24_climate_luna!(cs, ps, alb, sa, temp, pch, mask, bounds, dtime)
    CLM.acc24_climate_luna!(csD, psD, albD, saD, tempD, pchD, maskD, bounds, dtime); Metal.synchronize()
    for (nm,h,d) in (("t_veg_day",temp.t_veg_day_patch,tempD.t_veg_day_patch),
                     ("par24d_z",sa.par24d_z_patch,saD.par24d_z_patch),
                     ("par24x_z",sa.par24x_z_patch,saD.par24x_z_patch),
                     ("fpsn24",ps.fpsn24_patch,psD.fpsn24_patch))
        m,n = reldiff(h,d); ok = m < 1f-3
        @printf("  [%s] acc24.%-12s rel=%.2e over %d\n", ok ? "PASS" : "FAIL", nm, m, n); ok || (nfail+=1)
    end

    # ---- acc240 ----
    cs,ps,alb,sa,temp,pch,wdb,fv = build_fixture()
    saD,tempD,wdbD,fvD = mf(sa),mf(temp),mf(wdb),mf(fv)
    maskD = mf(mask)
    rb = fill(40.0, NP); rh = fill(0.6, NP); oair = fill(21000.0, NP); cair = fill(38.0, NP)
    CLM.acc240_climate_luna!(temp, ps, alb, sa, wdb, fv, pch, mask, bounds, oair, cair, rb, rh, dtime)
    CLM.acc240_climate_luna!(tempD, mf(ps), mf(alb), saD, wdbD, fvD, mf(pch), maskD,
        bounds, oair, cair, rb, rh, dtime); Metal.synchronize()
    for (nm,h,d) in (("par240d_z",sa.par240d_z_patch,saD.par240d_z_patch),
                     ("par240x_z",sa.par240x_z_patch,saD.par240x_z_patch),
                     ("t_veg10_day",temp.t_veg10_day_patch,tempD.t_veg10_day_patch),
                     ("t_veg10_night",temp.t_veg10_night_patch,tempD.t_veg10_night_patch),
                     ("rh10_af",wdb.rh10_af_patch,wdbD.rh10_af_patch),
                     ("rb10",fv.rb10_patch,fvD.rb10_patch))
        m,n = reldiff(h,d); ok = m < 1f-3
        @printf("  [%s] acc240.%-12s rel=%.2e over %d\n", ok ? "PASS" : "FAIL", nm, m, n); ok || (nfail+=1)
    end

    # ---- clear24 ----
    cs,ps,alb,sa,temp,pch,wdb,fv = build_fixture()
    psD,saD,tempD = mf(ps),mf(sa),mf(temp)
    CLM.clear24_climate_luna!(sa, ps, temp, pch, mask, bounds)
    CLM.clear24_climate_luna!(saD, psD, tempD, mf(pch), mf(mask), bounds); Metal.synchronize()
    for (nm,h,d) in (("par24d_z",sa.par24d_z_patch,saD.par24d_z_patch),
                     ("fpsn24",ps.fpsn24_patch,psD.fpsn24_patch))
        m,n = reldiff(h,d); ok = m < 1f-3
        @printf("  [%s] clear24.%-11s rel=%.2e over %d\n", ok ? "PASS" : "FAIL", nm, m, n); ok || (nfail+=1)
    end

    println()
    println(nfail == 0 ? "  LUNA climate kernels MATCH host on $name" : "  DIVERGENCE ($nfail)")
    return nfail == 0 ? 0 : 1
end
exit(main(detect_backend()))
