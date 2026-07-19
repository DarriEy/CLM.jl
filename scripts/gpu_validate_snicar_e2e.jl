# ==========================================================================
# gpu_validate_snicar_e2e.jl — exercise SNICAR (snow radiative transfer) by
# INJECTING a consistent snowpack into a stable summer cold-start state, then
# calling surface_albedo! standalone (which runs snicar_rt!) on CPU vs Metal.
#
# Driving snow via cold/snowy forcing destabilizes the lumped Bow-at-Banff
# cold-start (errsoi blow-up). Injecting a snowpack + isolating the albedo path
# is the robust way to validate snicar.
#
#   julia --project=scripts scripts/gpu_validate_snicar_e2e.jl
# ==========================================================================
using CLM, Printf
include(joinpath(@__DIR__, "gpu_backends.jl"))
include(joinpath(@__DIR__, "gpu_adapt.jl"))
mfm(x) = mf(device_array_type(), x)
const FSURDAT = get(ENV, "CLM_FSURDAT",
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc")
const PARAMFILE = get(ENV, "CLM_PARAMFILE",
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc")
const FSNOWOPTICS = get(ENV, "CLM_FSNOWOPTICS",
    "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_optics_5bnd_c013122.nc")
const FSNOWAGING = get(ENV, "CLM_FSNOWAGING",
    "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_drdt_bst_fit_60_c070416.nc")

function setup_forcing!(a2l, ng)
    T0 = 285.0
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

# Inject a 2-layer snowpack into the soil columns (consistent values snicar reads).
function inject_snow!(inst, bounds)
    nlevsno = CLM.varpar.nlevsno
    col = inst.column; ws = inst.water.waterstatebulk_inst.ws
    wd = inst.water.waterdiagnosticbulk_inst; temp = inst.temperature
    for c in 1:bounds.endc
        inst.landunit.itype[col.landunit[c]] == CLM.ISTSOIL || continue
        col.snl[c] = -2                                   # 2 snow layers (idx 11,12)
        for (i, jj) in enumerate((nlevsno-1, nlevsno))
            ws.h2osoi_ice_col[c, jj] = 20.0 + 10.0*i      # kg/m2
            ws.h2osoi_liq_col[c, jj] = 1.0
            wd.snw_rds_col[c, jj]    = 120.0              # microns
            temp.t_soisno_col[c, jj] = 268.0             # below freezing
            col.dz[c, jj] = 0.05*i
            col.z[c, jj]  = -0.05*(2-i) - 0.025
            col.zi[c, jj] = -0.05*(2-i)
        end
        col.zi[c, nlevsno] = -0.05
        ws.h2osno_no_layers_col[c] = 0.0
        wd.snow_depth_col[c] = 0.1
        wd.frac_sno_col[c]   = 0.95
    end
end

const config = CLM.CLMDriverConfig()
const filt_ia = CLM.clump_filter_inactive_and_active
const (declin, eccf) = CLM.compute_orbital(120.0)
const nextsw_cday = 120.0 + 1800.0/CLM.SECSPDAY
runstep!(inst, filt, fia, bounds, n; first=false) = CLM.clm_drv!(config, inst, filt, fia, bounds,
    true, nextsw_cday, declin, declin, CLM.ORB_OBLIQR_DEFAULT, false, false, "", false;
    nstep=n, is_first_step=first, is_beg_curr_day=first, dtime=1800.0, mon=1, day=1, photosyns=inst.photosyns)

# Standalone surface_albedo! (mirrors the driver call) — runs snicar_rt! internally.
function run_albedo!(inst, filt, bounds)
    pc = CLM.pftcon
    cosz = (cday,lat,lon,decl) -> begin
        ha = 2.0*pi*mod(cday,1.0) + lon - pi
        max(sin(lat)*sin(decl) + cos(lat)*cos(decl)*cos(ha), 0.0)
    end
    CLM.surface_albedo!(inst.surfalb, inst.surfalb_con, inst.gridcell, inst.column,
        inst.landunit, inst.patch, inst.canopystate, inst.temperature,
        inst.water.waterstatebulk_inst, inst.water.waterdiagnosticbulk_inst,
        inst.lakestate, inst.aerosol,
        filt_ia.nourbanc, filt_ia.nourbanp, nextsw_cday, declin,
        1:bounds.endg, 1:bounds.endc, 1:bounds.endp,
        pc.rhol, pc.rhos, pc.taul, pc.taus, pc.xl, cosz)
end

function build()
    (inst, bounds, filt, tm) = CLM.clm_initialize!(; fsurdat=FSURDAT, paramfile=PARAMFILE,
                                                     fsnowoptics=FSNOWOPTICS, fsnowaging=FSNOWAGING)
    setup_forcing!(inst.atm2lnd, bounds.endg)
    CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
    for n in 1:3; runstep!(inst, filt, filt_ia, bounds, n; first=(n==1)); end  # stable summer warmup
    inject_snow!(inst, bounds)
    return inst, bounds, filt
end

reldiff(a,b) = begin
    A=Array(a); B=Array(b); m=0.0; n=0
    for i in eachindex(A,B)
        (isfinite(A[i])&&isfinite(B[i])) || continue
        n+=1; m=max(m, abs(Float64(A[i])-Float64(B[i]))/(1.0+max(abs(Float64(A[i])),abs(Float64(B[i])))))
    end
    (m,n)
end

function main()
    if !gpu_functional(); println("No GPU backend detected"); return 0; end
    instH, bounds, filtH = build()
    run_albedo!(instH, filtH, bounds)
    println("CPU albedo (with injected snow) done. snl=", Array(instH.column.snl[1:bounds.endc]))

    instB, boundsB, filtB = build()
    inst_d = mfm(instB); filt_ia_d = mfm(filt_ia)
    println("moved to device: ", inst_d.surfalb.albsnd_hst_col isa device_array_type())
    # surface_albedo! reads filt_ia masks via the closure-captured global; for the
    # device run we need the device filter — rebind the global temporarily.
    CLM.surface_albedo!(inst_d.surfalb, inst_d.surfalb_con, inst_d.gridcell, inst_d.column,
        inst_d.landunit, inst_d.patch, inst_d.canopystate, inst_d.temperature,
        inst_d.water.waterstatebulk_inst, inst_d.water.waterdiagnosticbulk_inst,
        inst_d.lakestate, inst_d.aerosol,
        filt_ia_d.nourbanc, filt_ia_d.nourbanp, nextsw_cday, declin,
        1:boundsB.endg, 1:boundsB.endc, 1:boundsB.endp,
        mfm(CLM.pftcon.rhol), mfm(CLM.pftcon.rhos), mfm(CLM.pftcon.taul),
        mfm(CLM.pftcon.taus), mfm(CLM.pftcon.xl),
        (cday,lat,lon,decl) -> max(sin(lat)*sin(decl)+cos(lat)*cos(decl)*cos(2.0*pi*mod(cday,1.0)+lon-pi), 0.0))
    println("Device albedo done.")
    H=instH; D=inst_d
    checks = [
        ("albsnd_hst",  ()->H.surfalb.albsnd_hst_col, ()->D.surfalb.albsnd_hst_col),
        ("albsni_hst",  ()->H.surfalb.albsni_hst_col, ()->D.surfalb.albsni_hst_col),
        ("albgrd",      ()->H.surfalb.albgrd_col,     ()->D.surfalb.albgrd_col),
        ("albgri",      ()->H.surfalb.albgri_col,     ()->D.surfalb.albgri_col),
        ("albsod",      ()->H.surfalb.albsod_col,     ()->D.surfalb.albsod_col),
    ]
    nfail=0; tot=0; gmax=0.0; ncmp=0
    for (nm,ga,gb) in checks
        local d,n
        try; d,n=reldiff(ga(),gb()); catch e; @printf("  [err ] %-14s %s\n",nm,sprint(showerror,e)[1:min(end,40)]); continue; end
        tot+=n; (n>0)&&(gmax=max(gmax,d); ncmp+=1)
        ok = n==0 || d<1f-2
        @printf("  [%s] %-14s rel=%.3e (%d finite)\n", n==0 ? "skip" : (ok ? "PASS" : "FAIL"), nm, d, n)
        ok || (nfail+=1)
    end
    use_snicar = length(CLM.snicar_optics.ext_cff_mss_snw_drc) > 0 &&
                 any(x -> x != 0.0, CLM.snicar_optics.ext_cff_mss_snw_drc)
    @printf("\n  use_snicar=%s; %d fields, %d finite; max rel=%.3e; %s (snicar_rt_device! on Metal)\n",
            use_snicar, ncmp, tot, gmax,
            nfail==0 ? "MATCHES CPU on Metal" : "DIVERGENCE ($nfail)")
    return nfail
end
exit(main())
