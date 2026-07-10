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

# Optimization-triggering fixture for update_photosynthesis_capacity! (C3, fpsn24>0,
# tlai>0, valid t_veg_day) — self-contained (manual pft params + luna_params).
function build_luna_opt_fixture()
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
    npft = 20
    cs = CLM.CanopyStateData();  CLM.canopystate_init!(cs, NP)
    ps = CLM.PhotosynthesisData(); CLM.photosynthesis_data_init!(ps, NP; use_luna=true)
    alb = CLM.SurfaceAlbedoData(); CLM.surfalb_init!(alb, NP, NC, NG)
    sa = CLM.SolarAbsorbedData(); CLM.solarabs_init!(sa, NP, NL; use_luna=true)
    temp = CLM.TemperatureData(); CLM.temperature_init!(temp, NP, NC, NL, NG)
    pch = CLM.PatchData(); CLM.patch_init!(pch, NP)
    wdb = CLM.WaterDiagnosticBulkData(); CLM.waterdiagnosticbulk_init!(wdb, NC, NP, NL, NG)
    fv = CLM.FrictionVelocityData(); CLM.frictionvel_init!(fv, NP, NC)
    grc = CLM.GridcellData(); CLM.gridcell_init!(grc, NG)
    @inbounds for p in 1:NP
        pch.itype[p] = 0; pch.gridcell[p] = 1; pch.column[p] = 1
        temp.t_veg_day_patch[p] = 290.0
        temp.t_veg10_day_patch[p] = 296.8; temp.t_veg10_night_patch[p] = 280.0; temp.t_a10_patch[p] = 285.5
        ps.lnca_patch[p] = 2.98; ps.fpsn24_patch[p] = 10424.0
        ps.pnlc_z_patch[p,1] = 0.01; ps.enzs_z_patch[p,1] = 1.0
        ps.vcmx25_z_patch[p,1] = 49.56; ps.jmx25_z_patch[p,1] = 92.98
        ps.vcmx25_z_last_valid_patch[p,1] = 49.56; ps.jmx25_z_last_valid_patch[p,1] = 92.98
        cs.tlai_patch[p] = 0.0476
        sa.par240d_z_patch[p,1] = 164.4; sa.par240x_z_patch[p,1] = 249.0
        alb.nrad_patch[p] = 1; alb.tlai_z_patch[p,1] = 0.0476
        wdb.rh10_af_patch[p] = 0.199; fv.rb10_patch[p] = 42.6
    end
    grc.dayl[1] = 57263.8; grc.max_dayl[1] = 57263.8/0.9716
    c3psn = fill(1.0, npft); slatop = fill(0.012, npft); leafcn = fill(25.0, npft)
    rhol = fill(0.1, npft, 2); taul = fill(0.05, npft, 2)
    o3 = fill(1.0, NP); daylf = fill(0.944, NP)
    pbot = fill(78915.0, NP); co2 = fill(28.96, NP); o2 = fill(16493.0, NP)
    lp = CLM.LunaParamsData()
    lp.cp25_yr2000 = 42.75e-6*1e5; lp.kc25_coef = 404.9e-6*1e5; lp.ko25_coef = 278.4e-3*1e5
    lp.luna_theta_cj = 0.98; lp.enzyme_turnover_daily = 0.0114; lp.relhExp = 6.0686; lp.minrelh = 0.65
    lp.jmaxb0 = fill(0.0311, npft); lp.jmaxb1 = fill(0.1745, npft); lp.wc2wjb0 = fill(0.8054, npft)
    return (cs,ps,alb,sa,temp,pch,wdb,fv,grc, c3psn,slatop,leafcn,rhol,taul,o3,daylf,pbot,co2,o2, lp)
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

    # ---- update_photosynthesis_capacity! (the optimization core) ----
    F = build_luna_opt_fixture()
    (cs,ps,alb,sa,temp,pch,wdb,fv,grc, c3psn,slatop,leafcn,rhol,taul,o3,daylf,pbot,co2,o2, lp) = F
    maskU = trues(NP)
    CLM.update_photosynthesis_capacity!(ps, temp, cs, alb, sa, wdb, fv, pch, grc, maskU, 1:NP,
        daylf, pbot, co2, o2, c3psn, slatop, leafcn, rhol, taul, o3, lp, 3600.0, CLM.NLEVCAN)
    G = build_luna_opt_fixture()
    (cs2,ps2,alb2,sa2,temp2,pch2,wdb2,fv2,grc2, c3psn2,slatop2,leafcn2,rhol2,taul2,o32,daylf2,pbot2,co22,o22, lp2) = G
    psD = mf(ps2)
    CLM.update_photosynthesis_capacity!(psD, mf(temp2), mf(cs2), mf(alb2), mf(sa2), mf(wdb2),
        mf(fv2), mf(pch2), mf(grc2), mf(maskU), 1:NP,
        daylf2, pbot2, co22, o22, c3psn2, slatop2, leafcn2, rhol2, taul2, o32, lp2, 3600.0, CLM.NLEVCAN)
    Metal.synchronize()
    for (nm,h,d) in (("vcmx25_z",ps.vcmx25_z_patch,psD.vcmx25_z_patch),
                     ("jmx25_z",ps.jmx25_z_patch,psD.jmx25_z_patch),
                     ("pnlc_z",ps.pnlc_z_patch,psD.pnlc_z_patch),
                     ("enzs_z",ps.enzs_z_patch,psD.enzs_z_patch))
        m,n = reldiff(h,d); ok = m < 1f-4
        @printf("  [%s] update.%-10s rel=%.2e over %d  (host[1]=%.4f dev[1]=%.4f)\n",
                ok ? "PASS" : "FAIL", nm, m, n, Float64(Array(h)[1]), Float64(Array(d)[1])); ok || (nfail+=1)
    end

    println()
    println(nfail == 0 ? "  LUNA climate + capacity-update kernels MATCH host on $name" : "  DIVERGENCE ($nfail)")
    return nfail == 0 ? 0 : 1
end
exit(main(detect_backend()))
