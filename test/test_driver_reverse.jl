# Guards that the productionized clm_drv! reverse PHASES (src/driver/driver_reverse.jl)
# — soiltemp_rev_phase!, soilwater_rev_phase!, watertable_rev_phase! and the
# driver_rev_phases assembler — stay in sync with the production physics signatures.
# NO Enzyme here: a fast, suite-safe forward guard (the wrappers just call the physics
# functions, so a signature/aux-builder drift shows up as an error or a parity break).
# The reverse-AD gradients are validated on Julia 1.10 + Enzyme in
# scripts/enzyme_driver_reverse_full.jl (canopy+soil_temp+soil_water) and
# scripts/enzyme_driver_reverse_hydro.jl (+ water_table!).

@testset "driver reverse-phase forward parity" begin
    fsurdat = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc"
    paramfile = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc"
    if !isfile(fsurdat) || !isfile(paramfile)
        @warn "Skipping driver-reverse forward parity: input files not found"
        @test true
        return
    end

    (inst, bounds, filt, tm) = CLM.clm_initialize!(; fsurdat=fsurdat, paramfile=paramfile)
    for g in 1:bounds.endg
        a = inst.atm2lnd; a.forc_t_not_downscaled_grc[g] = 285.0
        a.forc_pbot_not_downscaled_grc[g] = 85000.0
        a.forc_th_not_downscaled_grc[g] = 285.0 * (100000.0/85000.0)^(CLM.RAIR/CLM.CPAIR)
        a.forc_rho_not_downscaled_grc[g] = 85000.0/(CLM.RAIR*285.0)
        a.forc_lwrad_not_downscaled_grc[g] = 300.0; a.forc_vp_grc[g] = 800.0
        a.forc_hgt_grc[g] = 30.0; a.forc_wind_grc[g] = 3.0
        a.forc_hgt_u_grc[g] = 30.0; a.forc_hgt_t_grc[g] = 30.0; a.forc_hgt_q_grc[g] = 30.0
        a.forc_topo_grc[g] = 0.0
        for b in 1:CLM.NUMRAD; a.forc_solad_not_downscaled_grc[g,b] = 200.0; a.forc_solai_grc[g,b] = 80.0; end
        a.forc_solar_not_downscaled_grc[g] = 560.0
        a.forc_rain_not_downscaled_grc[g] = 0.0001; a.forc_snow_not_downscaled_grc[g] = 0.0
    end
    CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
    config = CLM.CLMDriverConfig(); filt_ia = CLM.clump_filter_inactive_and_active
    (declin, _) = CLM.compute_orbital(120.0); nextsw = 120.0 + 1800.0/CLM.SECSPDAY
    for n in 1:3
        CLM.clm_drv!(config, inst, filt, filt_ia, bounds, true, nextsw, declin, declin,
            CLM.ORB_OBLIQR_DEFAULT, false, false, "", false;
            nstep=n, is_first_step=(n==1), is_beg_curr_day=(n==1), dtime=1800.0,
            mon=1, day=1, photosyns=inst.photosyns)
    end

    hc = filt.hydrologyc
    c0 = first(c for c in bounds.begc:bounds.endc if hc[c])
    j1 = CLM.varpar.nlevsno + 1

    # assembler (forward order): soil_temp + 10-phase surface-hydrology block +
    # soil_water + water_table + hydrology_no_drainage = 14 (15 with a canopy block);
    # include_surface=false reverts to the soil_temp→soil_water direct jump (4).
    ph = CLM.driver_rev_phases(bounds, filt, config)
    @test length(ph) == 14
    @test length(CLM.driver_rev_phases(bounds, filt, config; canopy_aux=(;), n_canopy=1)) == 15
    @test length(CLM.driver_rev_phases(bounds, filt, config; include_surface=false)) == 4
    @test length(CLM.surface_hydrology_rev_phases(bounds, filt)) == 10

    # run the full assembled chain forward on the shared-inst bundle.
    b = CLM.driver_rev_bundle(inst)
    for (f, cargs) in ph
        f(b, cargs...)
    end
    @test b.inst === inst                                   # bundle aliases the inst
    @test isfinite(b.inst.temperature.t_soisno_col[c0, j1])
    @test all(isfinite, b.inst.water.waterstatebulk_inst.ws.h2osoi_liq_col)
    @test isfinite(b.inst.soilhydrology.zwt_col[c0])
    @test isfinite(b.inst.water.waterstatebulk_inst.ws.h2osoi_vol_col[c0, CLM.varpar.nlevsno+4])
    @test isfinite(b.inst.water.waterfluxbulk_inst.wf.qflx_infl_col[c0])   # surface block ran

    # parity: soiltemp_rev_phase! reproduces a raw soil_temperature! call (i.e. the
    # soiltemp aux builder passes the same args the driver does).
    ref = deepcopy(inst)
    let i = ref
        CLM.soil_temperature!(i.column, i.landunit, i.patch, i.temperature, i.energyflux,
            i.soilstate, i.water.waterstatebulk_inst, i.water.waterdiagnosticbulk_inst,
            i.water.waterfluxbulk_inst, i.solarabs, i.canopystate, i.urbanparams,
            fill(323.15, bounds.endl), i.atm2lnd.forc_lwrad_downscaled_col,
            filt.nolakec, filt.nolakep, filt.urbanl, filt.urbanc,
            bounds.begc:bounds.endc, bounds.begl:bounds.endl, bounds.begp:bounds.endp, 1800.0)
    end
    b2 = CLM.driver_rev_bundle(deepcopy(inst))
    CLM.soiltemp_rev_phase!(b2, CLM.soiltemp_rev_aux(bounds, filt))
    @test b2.inst.temperature.t_soisno_col[c0, j1] ≈ ref.temperature.t_soisno_col[c0, j1] atol = 1e-12

    # pre-soil_water surface-hydrology phase wrappers run + produce finite outputs.
    bs = CLM.driver_rev_bundle(deepcopy(inst))
    sh = CLM.surfhydro_rev_aux(bounds, filt)
    CLM.satexcess_rev_phase!(bs, sh)      # fsat
    CLM.inflexcess_rev_phase!(bs, sh)     # qinmax (reads fsat)
    CLM.infil_rev_phase!(bs, sh)          # qflx_infl
    @test isfinite(bs.inst.sat_excess_runoff.fsat_col[c0])
    @test isfinite(bs.inst.infilt_excess_runoff.qinmax_col[c0])
    @test isfinite(bs.inst.water.waterfluxbulk_inst.wf.qflx_infl_col[c0])
end
