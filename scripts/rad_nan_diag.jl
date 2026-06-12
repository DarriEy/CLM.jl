# Diagnose the canopy radiation NaN: print coszen + albedo + solar-absorbed cascade
# after real-init + warmup, and test a DAYTIME nextsw_cday (coszen>0).
using CLM, Printf
const FSURDAT = get(ENV, "CLM_FSURDAT",
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc")
const PARAMFILE = get(ENV, "CLM_PARAMFILE",
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc")

function setup_forcing!(a2l, T, ng)
    for g in 1:ng
        a2l.forc_t_not_downscaled_grc[g]=T; a2l.forc_pbot_not_downscaled_grc[g]=85000.0
        a2l.forc_th_not_downscaled_grc[g]=T*(100000.0/85000.0)^(CLM.RAIR/CLM.CPAIR)
        a2l.forc_rho_not_downscaled_grc[g]=85000.0/(CLM.RAIR*T); a2l.forc_lwrad_not_downscaled_grc[g]=300.0
        a2l.forc_vp_grc[g]=800.0; a2l.forc_hgt_grc[g]=30.0; a2l.forc_topo_grc[g]=0.0; a2l.forc_wind_grc[g]=3.0
        # The component obs heights (set from forc_hgt_grc by forcing_reader, which this
        # harness bypasses) — without these, forc_hgt_u_grc=0 → zldis=z0m → ustar=0.4/0.
        a2l.forc_hgt_u_grc[g]=30.0; a2l.forc_hgt_t_grc[g]=30.0; a2l.forc_hgt_q_grc[g]=30.0
        for b in 1:CLM.NUMRAD; a2l.forc_solad_not_downscaled_grc[g,b]=200.0; a2l.forc_solai_grc[g,b]=80.0; end
        a2l.forc_solar_not_downscaled_grc[g]=560.0; a2l.forc_rain_not_downscaled_grc[g]=0.0001; a2l.forc_snow_not_downscaled_grc[g]=0.0
    end
end

function run(cday_for_orbital, nextsw_off_days, label)
    (inst, bounds, filt, tm) = CLM.clm_initialize!(; fsurdat=FSURDAT, paramfile=PARAMFILE)
    ng = bounds.endg
    config = CLM.CLMDriverConfig(); filt_ia = CLM.clump_filter_inactive_and_active
    (declin, eccf) = CLM.compute_orbital(cday_for_orbital)
    nextsw_cday = cday_for_orbital + nextsw_off_days
    setup_forcing!(inst.atm2lnd, 285.0, ng)
    CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
    for n in 1:3
        CLM.clm_drv!(config, inst, filt, filt_ia, bounds, true, nextsw_cday, declin, declin,
            CLM.ORB_OBLIQR_DEFAULT, false, false, "", false; nstep=n, is_first_step=(n==1),
            is_beg_curr_day=(n==1), dtime=1800.0, mon=1, day=1, photosyns=inst.photosyns)
    end
    sa = inst.surfalb; ev = filt.exposedvegp
    ps = [p for p in bounds.begp:bounds.endp if ev[p]]
    @printf("\n[%s]  cday_orbital=%.4f  nextsw_cday=%.4f  declin=%.4f  exposed-veg patches=%s\n",
        label, cday_for_orbital, nextsw_cday, declin, string(ps))
    @printf("  coszen_grc = %s\n", string(round.(sa.coszen_grc, sigdigits=4)))
    f(name, a) = @printf("  %-16s = %s\n", name, string([round(a[p], sigdigits=5) for p in ps]))
    fcol(name, a) = @printf("  %-16s = %s\n", name, string([round(a[inst.patch.column[p]], sigdigits=5) for p in ps]))
    pcol(name, M, ib) = @printf("  %-18s = %s\n", name, string([round(M[inst.patch.column[p], ib], sigdigits=5) for p in ps]))
    ppatch(name, M, ib) = @printf("  %-18s = %s\n", name, string([round(M[p, ib], sigdigits=5) for p in ps]))
    f("t_veg", inst.temperature.t_veg_patch)
    f("sabv", inst.solarabs.sabv_patch)
    ppatch("albd_patch[,1]", sa.albd_patch, 1)
    ppatch("albi_patch[,1]", sa.albi_patch, 1)
    pcol("albgrd_col[,1]", sa.albgrd_col, 1)   # ground (direct)
    pcol("albgri_col[,1]", sa.albgri_col, 1)   # ground (diffuse)
    pcol("albsod_col[,1]", sa.albsod_col, 1)   # soil (direct)
    pcol("albsoi_col[,1]", sa.albsoi_col, 1)   # soil (diffuse)
    pcol("albsnd_hst[,1]", sa.albsnd_hst_col, 1)  # snow (direct, history)
    ppatch("parsun_z[,1]", inst.solarabs.parsun_z_patch, 1)
    @printf("  nrad               = %s\n", string([sa.nrad_patch[p] for p in ps]))
    f("btran", inst.energyflux.btran_patch)
    @printf("  frac_sno_eff(col)  = %s\n", string([round(inst.water.waterdiagnosticbulk_inst.frac_sno_eff_col[inst.patch.column[p]], sigdigits=5) for p in ps]))
    @printf("  snl(col)           = %s\n", string([inst.column.snl[inst.patch.column[p]] for p in ps]))
    return nothing
end

# original (the canopy-probe setup): orbital day 120, nextsw +0.5 timestep
run(120.0, 1800.0/CLM.SECSPDAY, "ORIG day120 +1step")
# daytime: put nextsw_cday near local solar noon (cday + 0.5 day)
run(120.0, 0.5, "DAYTIME day120 +0.5day")
