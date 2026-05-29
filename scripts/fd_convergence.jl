using ForwardDiff
using ForwardDiff: Dual, value, partials
using CLM

function setup_forcing_conv!(a2l, T_forc, ng)
    for g in 1:ng
        a2l.forc_t_not_downscaled_grc[g] = T_forc
        a2l.forc_pbot_not_downscaled_grc[g] = 85000.0
        a2l.forc_th_not_downscaled_grc[g] = T_forc * (100000.0 / 85000.0)^(CLM.RAIR / CLM.CPAIR)
        a2l.forc_rho_not_downscaled_grc[g] = 85000.0 / (CLM.RAIR * T_forc)
        a2l.forc_lwrad_not_downscaled_grc[g] = 250.0
        a2l.forc_vp_grc[g] = 300.0
        a2l.forc_hgt_grc[g] = 30.0
        a2l.forc_topo_grc[g] = 0.0
        a2l.forc_wind_grc[g] = 3.0
        for b in 1:CLM.NUMRAD
            a2l.forc_solad_not_downscaled_grc[g, b] = 100.0
            a2l.forc_solai_grc[g, b] = 50.0
        end
        a2l.forc_solar_not_downscaled_grc[g] = 300.0
        a2l.forc_rain_not_downscaled_grc[g] = 0.0
        a2l.forc_snow_not_downscaled_grc[g] = 0.0001
    end
end

fsurdat = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc"
paramfile = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc"

T_base = 270.0
config = CLM.CLMDriverConfig()

function run_one(T_step4, fsurdat, paramfile, T_warmup)
    (inst, bounds, filt, tm) = CLM.clm_initialize!(; fsurdat=fsurdat, paramfile=paramfile)
    ng = bounds.endg
    filt_ia = CLM.clump_filter_inactive_and_active
    calday = 1.0
    (declin, eccf) = CLM.compute_orbital(calday)
    nextsw_cday = calday + 1800.0 / CLM.SECSPDAY

    setup_forcing_conv!(inst.atm2lnd, T_warmup, ng)
    CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
    for n in 1:3
        CLM.clm_drv!(config, inst, filt, filt_ia, bounds,
            true, nextsw_cday, declin, declin, CLM.ORB_OBLIQR_DEFAULT,
            false, false, "", false;
            nstep=n, is_first_step=(n==1), is_beg_curr_day=(n==1),
            dtime=1800.0, mon=1, day=1, photosyns=inst.photosyns)
    end

    setup_forcing_conv!(inst.atm2lnd, T_step4, ng)
    CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
    CLM.clm_drv!(config, inst, filt, filt_ia, bounds,
        true, nextsw_cday, declin, declin, CLM.ORB_OBLIQR_DEFAULT,
        false, false, "", false;
        nstep=4, is_first_step=false, dtime=1800.0, mon=1, day=1, photosyns=inst.photosyns)

    np = bounds.endp
    test_p = 0
    for p in 1:np
        lh = inst.energyflux.eflx_lh_tot_patch[p]
        if isfinite(lh) && abs(lh) > 0.0
            test_p = p; break
        end
    end
    test_p == 0 && (test_p = 1)

    return (inst.energyflux.eflx_lh_tot_patch[test_p],
            inst.energyflux.eflx_sh_tot_patch[test_p],
            inst.temperature.t_grnd_col[1])
end

println("Computing reference run at T=$T_base...")
lh_ref, sh_ref, tg_ref = run_one(T_base, fsurdat, paramfile, T_base)
println("  LH_ref=$(round(lh_ref, digits=4)), SH_ref=$(round(sh_ref, digits=4)), Tg_ref=$(round(tg_ref, digits=6))")

println("\n=== FD Convergence Study ===")
println("eps        | dlh_fd         | dsh_fd          | dtg_fd")
println("-"^70)

for eps_val in [1.0, 0.1, 0.01, 0.001, 0.0001]
    lh_p, sh_p, tg_p = run_one(T_base + eps_val, fsurdat, paramfile, T_base)
    dlh = (lh_p - lh_ref) / eps_val
    dsh = (sh_p - sh_ref) / eps_val
    dtg = (tg_p - tg_ref) / eps_val
    println("$(lpad(eps_val, 10)) | $(lpad(round(dlh, digits=4), 14)) | $(lpad(round(dsh, digits=4), 14)) | $(round(dtg, digits=6))")
end

println("\nExpected AD values (from ForwardDiff):")
println("  dlh_ad ≈ 65.88, dsh_ad ≈ -534.84, dtg_ad ≈ 0.0701")
println("\nFirst-order accuracy: FD values should converge to AD values as eps → 0")
