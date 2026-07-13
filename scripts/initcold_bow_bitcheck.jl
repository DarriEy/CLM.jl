# =============================================================================
# initcold_bow_bitcheck.jl — full-precision multi-step Bow dump, for proving the
# InitCold wiring left the validated default domain BIT-IDENTICAL.
#
# Prints every diagnostic at %.17g (round-trippable Float64) each step, so a
# plain `diff` of two runs is an exact bit comparison. Run once on main, once on
# the branch, diff.
#
#   julia +1.12 --project=. scripts/initcold_bow_bitcheck.jl > /tmp/bow_X.txt
# =============================================================================
using CLM, NCDatasets, Dates, Printf

const BOW_CAL = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/" *
                "domain_Bow_at_Banff_lumped/settings/CLM/parameters"
const BOW_FS  = joinpath(BOW_CAL, "surfdata_clm.nc")
const BOW_FP  = joinpath(BOW_CAL, "clm5_params.nc")
const SNOWOPT = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_optics_5bnd_c013122.nc"
const SNOWAGE = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_drdt_bst_fit_60_c070416.nc"

function set_forcing!(inst, ng, nc, step)
    a = inst.atm2lnd
    # Deterministic, mildly varying forcing so the run exercises melt/freeze,
    # canopy wet/dry and a day/night radiation cycle rather than a fixed state.
    T    = 278.0 + 8.0 * sin(2π * step / 24)
    sw   = max(0.0, 500.0 * sin(2π * (step - 6) / 24))
    rain = (step % 7 == 0) ? 2.0e-5 : 0.0
    snow = (step % 11 == 0) ? 1.0e-5 : 0.0
    pbot = 88000.0; q = 0.004
    th  = T * (1.0e5 / pbot)^0.286
    rho = pbot / (287.058 * T * (1.0 + 0.61 * q))
    vp  = q * pbot / (0.622 + 0.378 * q)
    for g in 1:ng
        a.forc_t_not_downscaled_grc[g]     = T
        a.forc_th_not_downscaled_grc[g]    = th
        a.forc_pbot_not_downscaled_grc[g]  = pbot
        a.forc_q_not_downscaled_grc[g]     = q
        a.forc_rho_not_downscaled_grc[g]   = rho
        a.forc_lwrad_not_downscaled_grc[g] = 280.0
        a.forc_rain_not_downscaled_grc[g]  = rain
        a.forc_snow_not_downscaled_grc[g]  = snow
        a.forc_u_grc[g] = 3.0; a.forc_v_grc[g] = 1.0
        isempty(a.forc_wind_grc) || (a.forc_wind_grc[g] = sqrt(10.0))
        a.forc_hgt_grc[g] = 30.0
        a.forc_hgt_u_grc[g] = 30.0; a.forc_hgt_t_grc[g] = 30.0; a.forc_hgt_q_grc[g] = 30.0
        a.forc_vp_grc[g] = vp
        a.forc_pco2_grc[g] = 367.0e-6 * pbot
        a.forc_po2_grc[g]  = 0.209 * pbot
        a.forc_topo_grc[g] = 1400.0
        for b in 1:size(a.forc_solad_not_downscaled_grc, 2)
            a.forc_solad_not_downscaled_grc[g, b] = b == 1 ? sw : 0.8 * sw
            a.forc_solai_grc[g, b]                = b == 1 ? 0.3 * sw : 0.25 * sw
        end
    end
    for c in 1:nc; inst.topo.topo_col[c] = 1400.0; end
    return nothing
end

# Every field a Bow parity/validation run actually looks at, at full precision.
function dump(inst, bounds, step)
    t   = inst.temperature
    ef  = inst.energyflux
    ws  = inst.water.waterstatebulk_inst.ws
    wf  = inst.water.waterfluxbulk_inst.wf
    wfb = inst.water.waterfluxbulk_inst
    wdb = inst.water.waterdiagnosticbulk_inst
    sa  = inst.solarabs
    sb  = inst.surfalb
    cs  = inst.canopystate
    ss  = inst.soilstate
    l2a = inst.lnd2atm
    emit(tag, v) = for (i, x) in enumerate(v)
        @printf("%04d %-26s %3d %.17g\n", step, tag, i, Float64(x))
    end
    emit("t_grnd",         t.t_grnd_col)
    emit("t_veg",          t.t_veg_patch)
    emit("t_ref2m",        t.t_ref2m_patch)
    emit("t_soisno",       vec(t.t_soisno_col))
    emit("t_lake",         vec(t.t_lake_col))
    emit("eflx_sh_tot",    ef.eflx_sh_tot_patch)
    emit("eflx_lh_tot",    ef.eflx_lh_tot_patch)
    emit("eflx_lwrad_out", ef.eflx_lwrad_out_patch)
    emit("eflx_soil_grnd", ef.eflx_soil_grnd_patch)
    emit("btran",          ef.btran_patch)
    emit("h2osoi_liq",     vec(ws.h2osoi_liq_col))
    emit("h2osoi_ice",     vec(ws.h2osoi_ice_col))
    emit("h2osoi_vol",     vec(ws.h2osoi_vol_col))
    emit("h2osno",         ws.h2osno_no_layers_col)
    emit("h2osfc",         ws.h2osfc_col)
    emit("wa",             ws.wa_col)
    emit("snow_depth",     wdb.snow_depth_col)
    emit("frac_sno",       wdb.frac_sno_col)
    emit("qflx_evap_tot",  wf.qflx_evap_tot_col)
    emit("qflx_tran_veg",  wf.qflx_tran_veg_patch)
    emit("qflx_surf",      wf.qflx_surf_col)
    emit("qflx_drain",     wf.qflx_drain_col)
    emit("qflx_snomelt",   wf.qflx_snomelt_col)
    emit("qflx_infl",      wf.qflx_infl_col)
    emit("sabv",           sa.sabv_patch)
    emit("sabg",           sa.sabg_patch)
    emit("fsa",            sa.fsa_patch)
    emit("albgrd",         vec(sb.albgrd_col))
    emit("albgri",         vec(sb.albgri_col))
    emit("albd",           vec(sb.albd_patch))
    emit("albi",           vec(sb.albi_patch))
    emit("elai",           cs.elai_patch)
    emit("esai",           cs.esai_patch)
    emit("laisun",         cs.laisun_patch)
    emit("laisha",         cs.laisha_patch)
    emit("smp_l",          vec(ss.smp_l_col))
    emit("eff_porosity",   vec(ss.eff_porosity_col))
    emit("psnsun",         inst.photosyns.psnsun_patch)
    emit("psnsha",         inst.photosyns.psnsha_patch)
    for f in (:eflx_sh_tot_grc, :eflx_lh_tot_grc, :eflx_lwrad_out_grc, :t_rad_grc,
              :t_ref2m_grc, :fsa_grc, :albd_grc, :albi_grc, :taux_grc, :tauy_grc,
              :qflx_ice_runoff_col, :qflx_liq_from_ice_col)
        v = getfield(l2a, f)
        isempty(v) || emit("l2a_" * String(f), vec(v))
    end
    return nothing
end

function main(; nsteps::Int = 48, use_cn::Bool = false)
    (isfile(BOW_FS) && isfile(BOW_FP)) || (println("SKIP: Bow inputs absent"); return)
    CLM.rooting_profile_config.rooting_profile_method_water  = CLM.JACKSON_1996_ROOT
    CLM.rooting_profile_config.rooting_profile_method_carbon = CLM.JACKSON_1996_ROOT
    dtime = 3600
    (inst, bounds, filt, tm) = CLM.clm_initialize!(; fsurdat=BOW_FS, paramfile=BOW_FP,
        start_date=DateTime(2006, 6, 15, 0), dtime=dtime, use_cn=use_cn, use_luna=false,
        use_bedrock=true, use_aquifer_layer=false, h2osfcflag=0,
        fsnowoptics=SNOWOPT, fsnowaging=SNOWAGE, int_snow_max=2000.0)
    ng, nc, np = bounds.endg, bounds.endc, bounds.endp
    @printf("# subgrid ng=%d nl=%d nc=%d np=%d use_cn=%s\n", ng, bounds.endl, nc, np, use_cn)
    dump(inst, bounds, 0)   # post-init state

    config  = CLM.CLMDriverConfig(use_cn=use_cn, use_aquifer_layer=false,
                                  use_hydrstress=false, use_luna=false)
    filt_ia = CLM.clump_filter_inactive_and_active
    for i in 1:nsteps
        set_forcing!(inst, ng, nc, i)
        calday = CLM.get_curr_calday(tm); (declin, _) = CLM.compute_orbital(calday)
        CLM.init_daylength!(inst.gridcell, declin, declin, CLM.ORB_OBLIQR_DEFAULT, 1:ng)
        CLM.advance_timestep!(tm)
        CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
        (yr, mon, d, tod) = CLM.get_curr_date(tm)
        nextsw = calday + dtime / CLM.SECSPDAY
        CLM.clm_drv!(config, inst, filt, filt_ia, bounds, true, nextsw, declin, declin,
            CLM.ORB_OBLIQR_DEFAULT, false, false, "", false;
            nstep=tm.nstep, is_first_step=(i == 1),
            is_beg_curr_day=CLM.is_beg_curr_day(tm), is_end_curr_day=CLM.is_end_curr_day(tm),
            is_beg_curr_year=CLM.is_beg_curr_year(tm), dtime=Float64(dtime), mon=mon, day=d,
            photosyns=inst.photosyns)
        dump(inst, bounds, i)
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(; nsteps = 48)
end
