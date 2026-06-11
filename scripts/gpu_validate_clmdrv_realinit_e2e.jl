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
runstep!(inst, filt, bounds, n; first=false) = CLM.clm_drv!(config, inst, filt, filt_ia, bounds,
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
    for n in 1:3; runstep!(instH, filtH, bounds, n; first=(n==1)); end   # CPU warmup
    instB, boundsB, filtB = build()
    for n in 1:3; runstep!(instB, filtB, boundsB, n; first=(n==1)); end
    inst_d = mfm(instB); filt_d = mfm(filtB)
    println("moved to device: ", inst_d.temperature.t_soisno_col isa Metal.MtlArray)
    runstep!(instH, filtH, bounds, 4)                 # CPU compared step
    runstep!(inst_d, filt_d, boundsB, 4)              # device compared step
    println("both step-4 runs completed")
    checks = [
        ("t_grnd", instH.temperature.t_grnd_col, inst_d.temperature.t_grnd_col),
        ("t_soisno", instH.temperature.t_soisno_col, inst_d.temperature.t_soisno_col),
        ("t_veg", instH.temperature.t_veg_patch, inst_d.temperature.t_veg_patch),
        ("eflx_sh_tot", instH.energyflux.eflx_sh_tot_patch, inst_d.energyflux.eflx_sh_tot_patch),
        ("sabg", instH.solarabs.sabg_patch, inst_d.solarabs.sabg_patch),
        ("fabd", instH.surfalb.fabd_patch, inst_d.surfalb.fabd_patch),
    ]
    nfail=0; tot=0; gmax=0.0
    for (nm,a,b) in checks
        d,n = reldiff(a,b); tot+=n; gmax=max(gmax,d)
        ok = n==0 || d<1f-2
        @printf("  [%s] %-12s rel=%.3e (%d finite)\n", n==0 ? "skip" : (ok ? "PASS" : "FAIL"), nm, d, n)
        ok || (nfail+=1)
    end
    @printf("\n  %d finite entries; global max rel=%.3e; %s\n", tot, gmax,
            nfail==0 ? "MATCHES CPU on Metal" : "DIVERGENCE ($nfail)")
    return nfail
end
exit(main())
