using Dates, Printf, CLM

basedir = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped"
init_fn = getfield(CLM, Symbol("clm_initialize!"))
(inst, bounds, filt, tm) = init_fn(;
    fsurdat=joinpath(basedir, "settings/CLM/parameters/surfdata_clm.nc"),
    paramfile=joinpath(basedir, "settings/CLM/parameters/clm5_params.nc"),
    start_date=DateTime(2002,7,1), dtime=1800)

nc = bounds.endc; np = bounds.endp
temp = inst.temperature
joff = CLM.varpar.nlevsno

println("=== INITIAL STATE ===")
@printf("T_GRND col 1 = %.4f K\n", temp.t_grnd_col[1])
@printf("t_soisno col 1, layer 1 = %.4f K\n", temp.t_soisno_col[1, 1 + joff])

# Run 3 steps
fr = CLM.ForcingReader()
CLM.forcing_reader_init!(fr, joinpath(basedir, "data/forcing/CLM_input/clmforc.2002.nc"))
config = CLM.CLMDriverConfig()
filt_ia = CLM.clump_filter_inactive_and_active

for step in 1:10
    CLM.advance_timestep!(tm)
    CLM.read_forcing_step!(fr, inst.atm2lnd, tm.current_date, bounds.endg, nc)
    CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
    calday = CLM.get_curr_calday(tm)
    (declin, eccf) = CLM.compute_orbital(calday)
    (yr, mon, d, tod) = CLM.get_curr_date(tm)
    CLM.clm_drv!(config, inst, filt, filt_ia, bounds,
                 true, calday + 1800.0/86400, declin, declin, CLM.ORB_OBLIQR_DEFAULT,
                 false, false, "", false; nstep=tm.nstep,
                 is_first_step=(tm.nstep==1), is_beg_curr_day=CLM.is_beg_curr_day(tm),
                 is_beg_curr_year=CLM.is_beg_curr_year(tm), dtime=1800.0,
                 mon=mon, day=d, photosyns=inst.photosyns)

    ef = inst.energyflux
    p = 2  # PFT 1 patch
    @printf("Step %d: T_GRND=%.4f, btran=%.6f\n", step, temp.t_grnd_col[1], ef.btran_patch[p])
    @printf("  eflx_gnet=%.4g, sabg=%.4g\n", ef.eflx_gnet_patch[p], inst.solarabs.sabg_patch[p])
    @printf("  eflx_sh_grnd=%.4g, eflx_lh_grnd=%.4g\n", ef.eflx_sh_grnd_patch[p], ef.eflx_lh_grnd_patch[p])
    @printf("  eflx_sh_veg=%.4g, eflx_sh_tot=%.4g\n", ef.eflx_sh_veg_patch[p], ef.eflx_sh_tot_patch[p])
    @printf("  eflx_lh_vegt=%.4g, eflx_lh_vege=%.4g, eflx_lh_tot=%.4g\n",
            ef.eflx_lh_vegt_patch[p], ef.eflx_lh_vege_patch[p], ef.eflx_lh_tot_patch[p])
    @printf("  lwrad_net=%.4g, lwrad_out=%.4g\n", ef.eflx_lwrad_net_patch[p], ef.eflx_lwrad_out_patch[p])
    @printf("  cgrnds=%.4g, cgrndl=%.4g, cgrnd=%.4g\n", ef.cgrnds_patch[p], ef.cgrndl_patch[p], ef.cgrnd_patch[p])
    @printf("  dlrad=%.4g, ulrad=%.4g\n", ef.dlrad_patch[p], ef.ulrad_patch[p])
    @printf("  t_veg=%.4f\n", temp.t_veg_patch[p])
end

println("\nDone.")
