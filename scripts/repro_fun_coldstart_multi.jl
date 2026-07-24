# Multi-step use_fun=true COLD START: does a NaN ever appear as LAI develops?
include(joinpath(@__DIR__, "fortran_parity_common.jl"))
using Printf

function main()
    NSTEPS = get(ENV, "NSTEPS", "72")
    nsteps = parse(Int, NSTEPS)
    reintroduce_hgt_bug = get(ENV, "HGT_BUG", "0") == "1"

    CLM.cnvegcstate_const.initial_vegC = 100.0
    CLM.rooting_profile_config.rooting_profile_method_water  = CLM.JACKSON_1996_ROOT
    CLM.rooting_profile_config.rooting_profile_method_carbon = CLM.JACKSON_1996_ROOT

    start_date = DateTime(2003, 7, 1)
    (inst, bounds, filt, tm) = CLM.clm_initialize!(;
        fsurdat=FSURDAT, paramfile=FPARAM,
        start_date=start_date, dtime=3600,
        use_cn=true, use_luna=true, use_bedrock=true, use_aquifer_layer=false,
        h2osfcflag=0, fsnowoptics=FSNOWOPT, fsnowaging=FSNOWAGE,
        int_snow_max=INT_SNOW_MAX)

    inst.bgc_vegetation.config.use_fun = true
    inst.bgc_vegetation.config.use_flexiblecn = true
    inst.bgc_vegetation.driver_config.use_fun = true
    inst.bgc_vegetation.driver_config.use_flexiblecn = true
    inst.cn_shared_params.use_fun = true
    inst.cn_shared_params.use_flexiblecn = true
    inst.cn_shared_params.br_root = 0.83e-6

    ccs = inst.bgc_vegetation.cnveg_carbonstate_inst
    cvf = inst.bgc_vegetation.cnveg_carbonflux_inst
    ps  = inst.photosyns

    filt_ia = CLM.clump_filter_inactive_and_active
    ng, nc, np = bounds.endg, bounds.endc, bounds.endp
    fr = CLM.ForcingReader(); CLM.forcing_reader_init!(fr, FFORCING)
    config = CLM.CLMDriverConfig(use_cn=true, use_aquifer_layer=false,
                                 use_hydrstress=true, use_luna=true)
    obliqr = CLM.ORB_OBLIQR_DEFAULT

    println("reintroduce_hgt_bug=", reintroduce_hgt_bug, "  nsteps=", nsteps)
    first_nan_step = 0
    for k in 1:nsteps
        step_start = tm.current_date
        CLM.advance_timestep!(tm)
        calday = CLM.get_curr_calday(tm)
        (declin, _) = CLM.compute_orbital(calday)
        (declinm1, _) = CLM.compute_orbital(calday - 3600.0 / CLM.SECSPDAY)
        CLM.init_daylength!(inst.gridcell, declin, declinm1, obliqr, 1:bounds.endg)
        CLM.read_forcing_step!(fr, inst.atm2lnd, step_start, ng, nc)
        CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
        if reintroduce_hgt_bug
            fill!(inst.atm2lnd.forc_hgt_u_grc, 0.0)
            fill!(inst.atm2lnd.forc_hgt_t_grc, 0.0)
            fill!(inst.atm2lnd.forc_hgt_q_grc, 0.0)
        end
        (yr, mon, d, tod) = CLM.get_curr_date(tm)
        CLM.clm_drv!(config, inst, filt, filt_ia, bounds,
                     true, calday, declin, declin, obliqr,
                     false, false, "", false;
                     nstep=tm.nstep, is_first_step=(tm.nstep==1),
                     is_beg_curr_day=CLM.is_beg_curr_day(tm),
                     is_end_curr_day=CLM.is_end_curr_day(tm),
                     is_beg_curr_year=CLM.is_beg_curr_year(tm),
                     dtime=3600.0, mon=mon, day=d,
                     jday=dayofyear(tm.current_date), secs=tod, year=yr,
                     photosyns=inst.photosyns)

        av = Array(cvf.availc_patch); ngw = Array(cvf.npp_growth_patch)
        lc = Array(ccs.leafc_patch); psn = Array(ps.psnsun_patch)
        lai = Array(inst.canopystate.laisun_patch)
        anynan = !all(isfinite, av) || !all(isfinite, ngw) || !all(isfinite, lc)
        if anynan && first_nan_step == 0
            first_nan_step = k
        end
        if k <= 5 || k % 12 == 0 || anynan
            @printf("step %3d: laisun[2]=%.4f psnsun[2]=%.4e availc[2]=%.4e npp_growth[2]=%.4e leafc[2]=%.4f  NaN=%s\n",
                    k, lai[2], psn[2], av[2], ngw[2], lc[2], anynan)
        end
        anynan && break
    end
    CLM.forcing_reader_close!(fr)
    println(first_nan_step == 0 ? "\nNO NaN over $nsteps steps." : "\nFIRST NaN at step $first_nan_step.")
    return 0
end
main()
