using Dates, Printf, CLM

basedir = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped"
init_fn = getfield(CLM, Symbol("clm_initialize!"))
(inst, bounds, filt, tm) = init_fn(;
    fsurdat=joinpath(basedir, "settings/CLM/parameters/surfdata_clm.nc"),
    paramfile=joinpath(basedir, "settings/CLM/parameters/clm5_params.nc"),
    start_date=DateTime(2002,7,1), dtime=1800)

nc = bounds.endc; ng = bounds.endg; np = bounds.endp

# Check rootfr after initialization
ss = inst.soilstate
println("=== ROOT FRACTIONS ===")
for p in 1:np
    ptype = inst.patch.itype[p]
    total_rootfr = sum(ss.rootfr_patch[p, j] for j in 1:CLM.varpar.nlevgrnd)
    @printf("Patch %d (PFT %d): rootfr sum = %.6f, first 5 layers: %s\n",
            p, ptype, total_rootfr,
            join([@sprintf("%.4f", ss.rootfr_patch[p,j]) for j in 1:5], ", "))
end

# Run 3 steps
fr = CLM.ForcingReader()
CLM.forcing_reader_init!(fr, joinpath(basedir, "data/forcing/CLM_input/clmforc.2002.nc"))
config = CLM.CLMDriverConfig()
filt_ia = CLM.clump_filter_inactive_and_active

for step in 1:3
    CLM.advance_timestep!(tm)
    CLM.read_forcing_step!(fr, inst.atm2lnd, tm.current_date, ng, nc)
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
end

println("\n=== RUNOFF COMPONENTS (after 3 steps) ===")
wf = inst.water.waterfluxbulk_inst.wf
for c in 1:nc
    @printf("Col %d: qflx_runoff=%.6g, drain=%.6g, surf=%.6g, qrgwl=%.6g, drain_p=%.6g, snwcp_liq=%.6g\n",
            c, wf.qflx_runoff_col[c], wf.qflx_drain_col[c], wf.qflx_surf_col[c],
            wf.qflx_qrgwl_col[c], wf.qflx_drain_perched_col[c], wf.qflx_snwcp_liq_col[c])
end

println("\n=== BTRAN (after 3 steps) ===")
for p in 1:np
    @printf("Patch %d (PFT %d): btran=%.6f\n", p, inst.patch.itype[p], inst.energyflux.btran_patch[p])
end

println("\nDone.")
