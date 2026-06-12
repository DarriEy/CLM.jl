# Make the real cold-start inst FINITE for end-to-end reverse-AD: the Bow-at-Banff
# lumped cold-start leaves NaN in soil water (h2osoi_*/liqvol top layers) and deep-
# layer hydraulic params (bsw/watsat below nlevsoi), which cascade to smp_l→soilpsi→
# btran and to the soil albedo → canopy t_veg/sabv NaN. We scrub those NaN to physical
# values once after init, warm up, and check canopy is finite. (This is a harness
# unblock for AD validation, not a cold-start physics fix.)
using CLM, Printf
const FSURDAT = get(ENV, "CLM_FSURDAT",
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc")
const PARAMFILE = get(ENV, "CLM_PARAMFILE",
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc")

# Replace NaN entries of array `A` with `v` (returns count scrubbed).
function scrub!(A, v)
    n = 0
    @inbounds for i in eachindex(A)
        if A[i] isa AbstractFloat && isnan(A[i]); A[i] = v; n += 1; end
    end
    return n
end

function finite_inst()
    (inst, bounds, filt, tm) = CLM.clm_initialize!(; fsurdat=FSURDAT, paramfile=PARAMFILE)
    ng = bounds.endg
    a2l = inst.atm2lnd; T0 = 285.0
    for g in 1:ng
        a2l.forc_t_not_downscaled_grc[g]=T0; a2l.forc_pbot_not_downscaled_grc[g]=85000.0
        a2l.forc_th_not_downscaled_grc[g]=T0*(100000.0/85000.0)^(CLM.RAIR/CLM.CPAIR)
        a2l.forc_rho_not_downscaled_grc[g]=85000.0/(CLM.RAIR*T0); a2l.forc_lwrad_not_downscaled_grc[g]=300.0
        a2l.forc_vp_grc[g]=800.0; a2l.forc_hgt_grc[g]=30.0; a2l.forc_topo_grc[g]=0.0; a2l.forc_wind_grc[g]=3.0
        for b in 1:CLM.NUMRAD; a2l.forc_solad_not_downscaled_grc[g,b]=200.0; a2l.forc_solai_grc[g,b]=80.0; end
        a2l.forc_solar_not_downscaled_grc[g]=560.0; a2l.forc_rain_not_downscaled_grc[g]=0.0001; a2l.forc_snow_not_downscaled_grc[g]=0.0
    end
    CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)

    # --- scrub soil hydraulic params + soil water to physical values ---
    ss = inst.soilstate; ws = inst.water.waterstatebulk_inst.ws; wd = inst.water.waterdiagnosticbulk_inst
    tot = 0
    tot += scrub!(ss.watsat_col, 0.45);  tot += scrub!(ss.bsw_col, 6.0)
    tot += scrub!(ss.sucsat_col, 150.0); tot += scrub!(ss.hksat_col, 5.0e-3)
    tot += scrub!(ss.smp_l_col, -50.0);  tot += scrub!(ss.soilpsi_col, -0.05)
    # mineral thermal-conductivity / heat-capacity params (init to NaN; cold-start
    # apparently leaves some layers unset → soil_temperature! tk/cv NaN → t_soisno NaN)
    for (f, v) in ((:tkmg_col, 2.0), (:tkdry_col, 0.3), (:tksatu_col, 2.5),
                   (:csol_col, 2.0e6), (:bd_col, 1300.0), (:thk_col, 1.5))
        hasproperty(ss, f) && (n = scrub!(getproperty(ss, f), v); tot += n;
            n > 0 && @printf("  scrubbed %d NaN in ss.%s\n", n, f))
    end
    for f in (:hk_l_col, :watdry_col, :watopt_col, :watfc_col, :eff_porosity_col, :rootr_patch, :rootr_col)
        hasproperty(ss, f) && (tot += scrub!(getproperty(ss, f), 0.3))
    end
    tot += scrub!(ws.h2osoi_liq_col, 30.0); tot += scrub!(ws.h2osoi_ice_col, 0.0)
    tot += scrub!(ws.h2osoi_vol_col, 0.3);  tot += scrub!(wd.h2osoi_liqvol_col, 0.3)
    # excess_ice + surface water (both init NaN; excess_ice feeds the tk satw term →
    # NaN thk → singular tridiagonal → whole soil column NaN)
    for (obj, f, v) in ((ws, :excess_ice_col, 0.0), (ws, :h2osfc_col, 0.0))
        hasproperty(obj, f) && (n = scrub!(getproperty(obj, f), v); tot += n;
            n > 0 && @printf("  scrubbed %d NaN in ws.%s\n", n, f))
    end
    # snow / surface water diagnostics sometimes NaN at cold start too
    for f in (:frac_sno_col, :frac_sno_eff_col, :frac_h2osfc_col, :snow_depth_col)
        hasproperty(wd, f) && (tot += scrub!(getproperty(wd, f), 0.0))
    end
    # albedo constant tables (soil wet/dry albedo by color) — check + scrub
    con = inst.surfalb_con
    @printf("albsat finite=%s  albdry finite=%s  isoicol=%s\n",
        isempty(con.albsat) ? "empty" : string(all(isfinite, con.albsat)),
        isempty(con.albdry) ? "empty" : string(all(isfinite, con.albdry)), string(con.isoicol))
    tot += scrub!(con.albsat, 0.2); tot += scrub!(con.albdry, 0.3)
    # temperatures (soil/ground/veg) — scrub NaN to a sane cold-start value
    tp = inst.temperature
    for (f, v) in ((:t_soisno_col, 283.0), (:t_grnd_col, 285.0), (:t_veg_patch, 285.0),
                   (:t_stem_patch, 285.0), (:thv_col, 285.0), (:t_h2osfc_col, 285.0),
                   (:t_skin_patch, 285.0), (:thm_patch, 285.0), (:t_a10_patch, 285.0),
                   (:emv_patch, 0.97), (:emg_col, 0.96))
        hasproperty(tp, f) && (tot += scrub!(getproperty(tp, f), v))
    end
    # column geometry z/dz/zi (init to NaN; snow slots NaN when snl=0). Scan first,
    # then scrub: z = cumulative depth, dz = thickness, zi = interface depth. Inactive
    # snow slots get small finite placeholders so the solve sees no NaN.
    col = inst.column
    fnl(M, c) = (for j in 1:size(M,2); isnan(M[c,j]) && return j; end; 0)
    @printf("col.z first-NaN layer per col = %s\n", string([fnl(col.z, c) for c in 1:bounds.endc]))
    @printf("col.dz first-NaN layer per col= %s\n", string([fnl(col.dz, c) for c in 1:bounds.endc]))
    @printf("col.zi first-NaN layer per col= %s\n", string([fnl(col.zi, c) for c in 1:bounds.endc]))
    tot += scrub!(col.dz, 0.1)
    # z/zi: fill NaN with a monotonic placeholder so layer spacings stay positive
    for c in 1:bounds.endc
        for j in 1:size(col.z, 2); isnan(col.z[c,j]) && (col.z[c,j] = 0.05 + 0.1*(j-1); tot += 1); end
        for j in 1:size(col.zi, 2); isnan(col.zi[c,j]) && (col.zi[c,j] = 0.1*(j-1); tot += 1); end
    end
    @printf("scrubbed %d NaN entries total\n", tot)
    return inst, bounds, filt
end

inst, bounds, filt = finite_inst()
config = CLM.CLMDriverConfig(); filt_ia = CLM.clump_filter_inactive_and_active
(declin, eccf) = CLM.compute_orbital(120.0); nextsw_cday = 120.0 + 1800.0/CLM.SECSPDAY
for n in 1:3
    CLM.clm_drv!(config, inst, filt, filt_ia, bounds, true, nextsw_cday, declin, declin,
        CLM.ORB_OBLIQR_DEFAULT, false, false, "", false; nstep=n, is_first_step=(n==1),
        is_beg_curr_day=(n==1), dtime=1800.0, mon=1, day=1, photosyns=inst.photosyns)
end

ev = filt.exposedvegp
ps = [p for p in bounds.begp:bounds.endp if ev[p]]
@printf("\nexposed-veg patches = %s\n", string(ps))
chk(name, a) = @printf("  %-14s = %s   allfinite=%s\n", name, string([round(a[p],sigdigits=6) for p in ps]), all(isfinite, a[p] for p in ps))
chk("t_veg", inst.temperature.t_veg_patch)
chk("sabv", inst.solarabs.sabv_patch)
chk("btran", inst.energyflux.btran_patch)
chk("t_grnd(col)", [inst.temperature.t_grnd_col[inst.patch.column[p]] for p in 1:bounds.endp])
@printf("  t_soisno allfinite (col1) = %s\n", all(isfinite, inst.temperature.t_soisno_col))
# which combined layers of t_soisno are NaN, per exposed column (nlevsno=12 snow slots)
let ts = inst.temperature.t_soisno_col, nls = CLM.varpar.nlevsno
    for p in ps
        c = inst.patch.column[p]
        nanlayers = [j for j in 1:size(ts,2) if isnan(ts[c,j])]
        @printf("  col %d snl=%d  t_soisno NaN combined-layers = %s  (snow=1:%d, soil=%d:%d)\n",
            c, inst.column.snl[c], string(nanlayers), nls, nls+1, size(ts,2))
        @printf("    t_h2osfc=%.4g  frac_h2osfc=%.4g  frac_sno_eff=%.4g\n",
            inst.temperature.t_h2osfc_col[c],
            inst.water.waterdiagnosticbulk_inst.frac_h2osfc_col[c],
            inst.water.waterdiagnosticbulk_inst.frac_sno_eff_col[c])
    end
end
allfin = all(isfinite, inst.temperature.t_veg_patch[p] for p in ps) &&
         all(isfinite, inst.solarabs.sabv_patch[p] for p in ps) &&
         all(isfinite, inst.energyflux.btran_patch[p] for p in ps)
println(allfin ? "\nCANOPY IS FINITE after scrub ✓" : "\nstill NaN — more scrubbing needed ✗")
