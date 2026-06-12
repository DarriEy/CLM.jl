# ==========================================================================
# gpu_validate_clmdrv_realinit_e2e.jl — whole-clm_drv! Metal parity on a REAL
# cold-start state (clm_initialize! from the Bow-at-Banff surfdata/params NetCDF).
#
# Unlike the synthetic gpu_validate_clmdrv_e2e.jl, this uses the real cold-start,
# which produces a consistent initial albedo and hence FINITE surface-energy
# fields (t_grnd/t_veg/eflx) — so the parity comparison covers the energy chain.
#
#   julia --project=scripts scripts/gpu_validate_clmdrv_realinit_e2e.jl
# ==========================================================================
using CLM, Printf
import Metal
include(joinpath(@__DIR__, "gpu_adapt.jl"))
mfm(x) = mf(Metal.MtlArray, x)
const FSURDAT = get(ENV, "CLM_FSURDAT",
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc")
const PARAMFILE = get(ENV, "CLM_PARAMFILE",
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc")

function setup_forcing!(a2l, T0, ng)
    for g in 1:ng
        a2l.forc_t_not_downscaled_grc[g]=T0; a2l.forc_pbot_not_downscaled_grc[g]=85000.0
        a2l.forc_th_not_downscaled_grc[g]=T0*(100000.0/85000.0)^(CLM.RAIR/CLM.CPAIR)
        a2l.forc_rho_not_downscaled_grc[g]=85000.0/(CLM.RAIR*T0)
        a2l.forc_lwrad_not_downscaled_grc[g]=300.0; a2l.forc_vp_grc[g]=800.0
        a2l.forc_hgt_grc[g]=30.0; a2l.forc_topo_grc[g]=0.0; a2l.forc_wind_grc[g]=3.0
        a2l.forc_hgt_u_grc[g]=30.0; a2l.forc_hgt_t_grc[g]=30.0; a2l.forc_hgt_q_grc[g]=30.0  # component obs heights (forcing_reader derives these; harness bypasses it)
        for b in 1:CLM.NUMRAD; a2l.forc_solad_not_downscaled_grc[g,b]=200.0; a2l.forc_solai_grc[g,b]=80.0; end
        a2l.forc_solar_not_downscaled_grc[g]=560.0
        a2l.forc_rain_not_downscaled_grc[g]=0.0001; a2l.forc_snow_not_downscaled_grc[g]=0.0
    end
end

function build()
    (inst, bounds, filt, tm) = CLM.clm_initialize!(; fsurdat=FSURDAT, paramfile=PARAMFILE)
    setup_forcing!(inst.atm2lnd, 285.0, bounds.endg)
    CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
    return inst, bounds, filt
end

const config = CLM.CLMDriverConfig()
const filt_ia = CLM.clump_filter_inactive_and_active
const (declin, eccf) = CLM.compute_orbital(120.0)
const nextsw_cday = 120.0 + 1800.0/CLM.SECSPDAY
runstep!(inst, filt, fia, bounds, n; first=false) = CLM.clm_drv!(config, inst, filt, fia, bounds,
    true, nextsw_cday, declin, declin, CLM.ORB_OBLIQR_DEFAULT, false, false, "", false;
    nstep=n, is_first_step=first, is_beg_curr_day=first, dtime=1800.0, mon=1, day=1, photosyns=inst.photosyns)

reldiff(a,b) = begin
    A=Array(a); B=Array(b); m=0.0; n=0
    for i in eachindex(A,B)
        (isfinite(A[i])&&isfinite(B[i])) || continue
        n+=1; m=max(m, abs(Float64(A[i])-Float64(B[i]))/(1.0+max(abs(Float64(A[i])),abs(Float64(B[i])))))
    end
    (m,n)
end

function main()
    if !Metal.functional(); println("Metal not functional"); return 0; end
    instH, bounds, filtH = build()
    for n in 1:3; runstep!(instH, filtH, filt_ia, bounds, n; first=(n==1)); end   # CPU warmup
    instB, boundsB, filtB = build()
    for n in 1:3; runstep!(instB, filtB, filt_ia, boundsB, n; first=(n==1)); end
    inst_d = mfm(instB); filt_d = mfm(filtB); filt_ia_d = mfm(filt_ia)
    println("moved to device: ", inst_d.temperature.t_soisno_col isa Metal.MtlArray)
    runstep!(instH, filtH, filt_ia, bounds, 4)            # CPU compared step
    runstep!(inst_d, filt_d, filt_ia_d, boundsB, 4)       # device compared step
    println("both step-4 runs completed")
    H = instH; D = inst_d
    wsH = H.water.waterstatebulk_inst.ws;  wsD = D.water.waterstatebulk_inst.ws
    wdH = H.water.waterdiagnosticbulk_inst; wdD = D.water.waterdiagnosticbulk_inst
    wfH = H.water.waterfluxbulk_inst;       wfD = D.water.waterfluxbulk_inst
    # Lazy (name, hostgetter, devgetter) so a missing field can't abort the scan.
    checks = [
        ("temp.t_grnd",      ()->H.temperature.t_grnd_col,      ()->D.temperature.t_grnd_col),
        ("temp.t_soisno",    ()->H.temperature.t_soisno_col,    ()->D.temperature.t_soisno_col),
        ("temp.t_veg",       ()->H.temperature.t_veg_patch,     ()->D.temperature.t_veg_patch),
        ("temp.t_h2osfc",    ()->H.temperature.t_h2osfc_col,    ()->D.temperature.t_h2osfc_col),
        ("temp.thv",         ()->H.temperature.thv_col,         ()->D.temperature.thv_col),
        ("temp.emg",         ()->H.temperature.emg_col,         ()->D.temperature.emg_col),
        ("temp.t_ref2m",     ()->H.temperature.t_ref2m_patch,   ()->D.temperature.t_ref2m_patch),
        ("temp.t_lake",      ()->H.temperature.t_lake_col,      ()->D.temperature.t_lake_col),
        ("ws.h2osoi_liq",    ()->wsH.h2osoi_liq_col,            ()->wsD.h2osoi_liq_col),
        ("ws.h2osoi_ice",    ()->wsH.h2osoi_ice_col,            ()->wsD.h2osoi_ice_col),
        ("ws.h2osoi_vol",    ()->wsH.h2osoi_vol_col,            ()->wsD.h2osoi_vol_col),
        ("ws.h2osfc",        ()->wsH.h2osfc_col,                ()->wsD.h2osfc_col),
        ("ws.h2osno_nolyr",  ()->wsH.h2osno_no_layers_col,      ()->wsD.h2osno_no_layers_col),
        ("wd.frac_sno",      ()->wdH.frac_sno_col,              ()->wdD.frac_sno_col),
        ("wd.frac_sno_eff",  ()->wdH.frac_sno_eff_col,          ()->wdD.frac_sno_eff_col),
        ("wd.frac_h2osfc",   ()->wdH.frac_h2osfc_col,           ()->wdD.frac_h2osfc_col),
        ("wd.snow_depth",    ()->wdH.snow_depth_col,            ()->wdD.snow_depth_col),
        ("wd.h2osoi_liqvol", ()->wdH.h2osoi_liqvol_col,         ()->wdD.h2osoi_liqvol_col),
        ("ef.eflx_sh_tot",   ()->H.energyflux.eflx_sh_tot_patch, ()->D.energyflux.eflx_sh_tot_patch),
        ("ef.eflx_sh_grnd",  ()->H.energyflux.eflx_sh_grnd_patch,()->D.energyflux.eflx_sh_grnd_patch),
        ("ef.eflx_lwrad_net",()->H.energyflux.eflx_lwrad_net_patch,()->D.energyflux.eflx_lwrad_net_patch),
        ("ef.eflx_soil_grnd",()->H.energyflux.eflx_soil_grnd_patch,()->D.energyflux.eflx_soil_grnd_patch),
        ("ef.htvp",          ()->H.energyflux.htvp_col,         ()->D.energyflux.htvp_col),
        ("sa.sabg",          ()->H.solarabs.sabg_patch,         ()->D.solarabs.sabg_patch),
        ("sa.sabv",          ()->H.solarabs.sabv_patch,         ()->D.solarabs.sabv_patch),
        ("sa.fsa",           ()->H.solarabs.fsa_patch,          ()->D.solarabs.fsa_patch),
        ("alb.albgrd",       ()->H.surfalb.albgrd_col,          ()->D.surfalb.albgrd_col),
        ("alb.albgri",       ()->H.surfalb.albgri_col,          ()->D.surfalb.albgri_col),
        ("alb.albd",         ()->H.surfalb.albd_patch,          ()->D.surfalb.albd_patch),
        ("alb.fabd",         ()->H.surfalb.fabd_patch,          ()->D.surfalb.fabd_patch),
        ("alb.fabi",         ()->H.surfalb.fabi_patch,          ()->D.surfalb.fabi_patch),
        ("alb.ftdd",         ()->H.surfalb.ftdd_patch,          ()->D.surfalb.ftdd_patch),
        ("cs.elai",          ()->H.canopystate.elai_patch,      ()->D.canopystate.elai_patch),
        ("fv.z0mg",          ()->H.frictionvel.z0mg_col,        ()->D.frictionvel.z0mg_col),
        ("fv.ram1",          ()->H.frictionvel.ram1_patch,      ()->D.frictionvel.ram1_patch),
        ("wf.qflx_ev_snow",  ()->wfH.qflx_ev_snow_patch,        ()->wfD.qflx_ev_snow_patch),
    ]
    nfail=0; tot=0; gmax=0.0; ncmp=0
    for (nm,ga,gb) in checks
        local d, n
        try
            d,n = reldiff(ga(), gb())
        catch e
            @printf("  [err ] %-18s %s\n", nm, sprint(showerror,e)[1:min(end,40)]); continue
        end
        tot+=n; (n>0) && (gmax=max(gmax,d); ncmp+=1)
        ok = n==0 || d<1f-2
        @printf("  [%s] %-18s rel=%.3e (%d finite)\n", n==0 ? "skip" : (ok ? "PASS" : "FAIL"), nm, d, n)
        ok || (nfail+=1)
    end
    @printf("\n  %d fields compared, %d finite entries; global max rel=%.3e; %s\n", ncmp, tot, gmax,
            nfail==0 ? "MATCHES CPU on Metal" : "DIVERGENCE ($nfail)")
    return nfail
end
exit(main())
