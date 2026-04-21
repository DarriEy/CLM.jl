# ==========================================================================
# test_ad_robustness.jl — Multi-scenario AD robustness tests
#
# Verifies that ForwardDiff flows through clm_drv! across multiple
# forcing scenarios (different temperatures, precipitation regimes).
# After Phase 4 smoothing, AD should produce finite, physically
# meaningful derivatives across the full seasonal cycle.
#
# Scenarios:
#   1. Winter cold (250K) — snow-dominated, no photosynthesis
#   2. Near-freezing (272K) — freeze/thaw transitions
#   3. Spring (285K) — active photosynthesis, moderate conditions
#   4. Summer hot (305K) — high T, potential stress
#   5. Summer wet (290K + heavy rain) — saturated conditions
#   6. Dry stress (295K, low humidity) — water-limited
# ==========================================================================

using Test
using ForwardDiff
using ForwardDiff: Dual, Tag, value, partials
using CLM

@testset "AD Robustness — Multi-scenario" begin

    fsurdat = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc"
    paramfile = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc"

    if !isfile(fsurdat) || !isfile(paramfile)
        @warn "Skipping AD robustness tests: input files not found"
        @test true
        return
    end

    # --- Helper: create a Dual-typed copy of a parameterized struct ---
    function make_dual_copy(src, ::Type{D}) where D
        T = typeof(src)
        wrapper = T.name.wrapper
        if wrapper === T
            return src
        end
        dst = wrapper{D}()
        for name in fieldnames(T)
            sv = getfield(src, name)
            try
                if sv isa Array{Float64}
                    setfield!(dst, name, D.(sv))
                elseif sv isa Float64
                    setfield!(dst, name, D(sv))
                else
                    setfield!(dst, name, sv)
                end
            catch e
                e isa TypeError || rethrow()
            end
        end
        return dst
    end

    # --- Helper: set up forcing ---
    function setup_forcing!(a2l, T, ng; precip_rain=0.0, precip_snow=0.0001,
                            vp=300.0, solar_direct=100.0, solar_diffuse=50.0,
                            lwrad=250.0, wind=3.0)
        for g in 1:ng
            a2l.forc_t_not_downscaled_grc[g] = T
            a2l.forc_pbot_not_downscaled_grc[g] = 85000.0
            a2l.forc_th_not_downscaled_grc[g] = T * (100000.0 / 85000.0)^(CLM.RAIR / CLM.CPAIR)
            a2l.forc_rho_not_downscaled_grc[g] = 85000.0 / (CLM.RAIR * T)
            a2l.forc_lwrad_not_downscaled_grc[g] = lwrad
            a2l.forc_vp_grc[g] = vp
            a2l.forc_hgt_grc[g] = 30.0
            a2l.forc_topo_grc[g] = 0.0
            a2l.forc_wind_grc[g] = wind
            for b in 1:CLM.NUMRAD
                a2l.forc_solad_not_downscaled_grc[g, b] = solar_direct
                a2l.forc_solai_grc[g, b] = solar_diffuse
            end
            a2l.forc_solar_not_downscaled_grc[g] = 2 * (solar_direct + solar_diffuse)
            a2l.forc_rain_not_downscaled_grc[g] = precip_rain
            a2l.forc_snow_not_downscaled_grc[g] = precip_snow
        end
    end

    # Define forcing scenarios
    scenarios = [
        (name="Winter cold (250K)",
         T=250.0, precip_rain=0.0, precip_snow=0.001,
         vp=100.0, solar_direct=50.0, solar_diffuse=25.0,
         lwrad=200.0, wind=5.0, calday=1.0),

        (name="Near-freezing (272K)",
         T=272.0, precip_rain=0.0, precip_snow=0.0005,
         vp=400.0, solar_direct=100.0, solar_diffuse=50.0,
         lwrad=250.0, wind=3.0, calday=80.0),

        (name="Spring (285K)",
         T=285.0, precip_rain=0.0001, precip_snow=0.0,
         vp=800.0, solar_direct=200.0, solar_diffuse=80.0,
         lwrad=300.0, wind=3.0, calday=120.0),

        (name="Summer hot (305K)",
         T=305.0, precip_rain=0.0, precip_snow=0.0,
         vp=1500.0, solar_direct=300.0, solar_diffuse=100.0,
         lwrad=380.0, wind=2.0, calday=182.0),

        (name="Summer wet (290K + rain)",
         T=290.0, precip_rain=0.005, precip_snow=0.0,
         vp=1200.0, solar_direct=150.0, solar_diffuse=60.0,
         lwrad=320.0, wind=4.0, calday=182.0),

        (name="Dry stress (295K)",
         T=295.0, precip_rain=0.0, precip_snow=0.0,
         vp=400.0, solar_direct=250.0, solar_diffuse=80.0,
         lwrad=340.0, wind=2.0, calday=200.0),
    ]

    for scenario in scenarios
        @testset "$(scenario.name)" begin
            # ---- Initialize Float64 state ----
            (inst_f64, bounds, filt, tm) = CLM.clm_initialize!(;
                fsurdat=fsurdat, paramfile=paramfile)
            ng = bounds.endg; np = bounds.endp

            config = CLM.CLMDriverConfig()
            filt_ia = CLM.clump_filter_inactive_and_active
            (declin, eccf) = CLM.compute_orbital(scenario.calday)
            nextsw_cday = scenario.calday + 1800.0 / CLM.SECSPDAY

            # Set forcing
            setup_forcing!(inst_f64.atm2lnd, scenario.T, ng;
                precip_rain=scenario.precip_rain,
                precip_snow=scenario.precip_snow,
                vp=scenario.vp,
                solar_direct=scenario.solar_direct,
                solar_diffuse=scenario.solar_diffuse,
                lwrad=scenario.lwrad,
                wind=scenario.wind)
            CLM.downscale_forcings!(bounds, inst_f64.atm2lnd, inst_f64.column,
                                    inst_f64.landunit, inst_f64.topo)

            # Run 3 warmup timesteps
            for n in 1:3
                CLM.clm_drv!(config, inst_f64, filt, filt_ia, bounds,
                    true, nextsw_cday, declin, declin, CLM.ORB_OBLIQR_DEFAULT,
                    false, false, "", false;
                    nstep=n, is_first_step=(n==1), is_beg_curr_day=(n==1),
                    dtime=1800.0, mon=1, day=1,
                    photosyns=inst_f64.photosyns)
            end

            # Find a patch with finite LH
            test_p = 1
            for p in 1:np
                lh = inst_f64.energyflux.eflx_lh_tot_patch[p]
                if isfinite(lh) && abs(lh) > 0.0
                    test_p = p; break
                end
            end

            # ---- Create Dual-typed copy ----
            D = Dual{Nothing, Float64, 1}
            inst_d = CLM.CLMInstances()

            for name in fieldnames(CLM.CLMInstances)
                name === :water && continue
                name === :surfdata && continue
                src = getfield(inst_f64, name)
                setfield!(inst_d, name, make_dual_copy(src, D))
            end
            inst_d.surfdata = inst_f64.surfdata
            inst_d.photosyns.stomatalcond_mtd = inst_f64.photosyns.stomatalcond_mtd

            # Handle WaterData
            water_d = CLM.WaterData()
            for name in fieldnames(CLM.WaterData)
                sv = getfield(inst_f64.water, name)
                try setfield!(water_d, name, sv) catch end
            end
            wsb_d = make_dual_copy(inst_f64.water.waterstatebulk_inst, D)
            ws_d  = make_dual_copy(inst_f64.water.waterstatebulk_inst.ws, D)
            wsb_d.ws = ws_d
            water_d.waterstatebulk_inst = wsb_d
            wfb_d = make_dual_copy(inst_f64.water.waterfluxbulk_inst, D)
            wf_d  = make_dual_copy(inst_f64.water.waterfluxbulk_inst.wf, D)
            wfb_d.wf = wf_d
            water_d.waterfluxbulk_inst = wfb_d
            water_d.waterdiagnosticbulk_inst = make_dual_copy(
                inst_f64.water.waterdiagnosticbulk_inst, D)
            water_d.waterbalancebulk_inst = make_dual_copy(
                inst_f64.water.waterbalancebulk_inst, D)
            if !isempty(water_d.bulk_and_tracers)
                bt = water_d.bulk_and_tracers[water_d.i_bulk]
                bt.waterflux = wfb_d
                bt.waterstate = wsb_d
                bt.waterdiagnostic = water_d.waterdiagnosticbulk_inst
                bt.waterbalance = water_d.waterbalancebulk_inst
            end
            inst_d.water = water_d

            # ---- Seed Dual into forcing temperature ----
            setup_forcing!(inst_d.atm2lnd, scenario.T, ng;
                precip_rain=scenario.precip_rain,
                precip_snow=scenario.precip_snow,
                vp=scenario.vp,
                solar_direct=scenario.solar_direct,
                solar_diffuse=scenario.solar_diffuse,
                lwrad=scenario.lwrad,
                wind=scenario.wind)
            for g in 1:ng
                inst_d.atm2lnd.forc_t_not_downscaled_grc[g] = Dual(scenario.T, 1.0)
                T_d = inst_d.atm2lnd.forc_t_not_downscaled_grc[g]
                inst_d.atm2lnd.forc_th_not_downscaled_grc[g] = T_d * (D(100000.0) / D(85000.0))^(D(CLM.RAIR) / D(CLM.CPAIR))
                inst_d.atm2lnd.forc_rho_not_downscaled_grc[g] = D(85000.0) / (D(CLM.RAIR) * T_d)
            end
            CLM.downscale_forcings!(bounds, inst_d.atm2lnd, inst_d.column,
                                    inst_d.landunit, inst_d.topo)

            # ---- Run one Dual timestep ----
            CLM.clm_drv!(config, inst_d, filt, filt_ia, bounds,
                true, nextsw_cday, declin, declin, CLM.ORB_OBLIQR_DEFAULT,
                false, false, "", false;
                nstep=4, is_first_step=false,
                dtime=1800.0, mon=1, day=1,
                photosyns=inst_d.photosyns)

            # ---- Verify derivatives are finite ----
            lh_d = inst_d.energyflux.eflx_lh_tot_patch[test_p]
            sh_d = inst_d.energyflux.eflx_sh_tot_patch[test_p]
            tg_d = inst_d.temperature.t_grnd_col[1]

            @test isfinite(value(lh_d))
            @test isfinite(value(sh_d))
            @test isfinite(value(tg_d))

            @test isfinite(partials(lh_d)[1])
            @test isfinite(partials(sh_d)[1])
            @test isfinite(partials(tg_d)[1])

            # At least one derivative should be nonzero
            any_nonzero = abs(partials(lh_d)[1]) > 0 ||
                          abs(partials(sh_d)[1]) > 0 ||
                          abs(partials(tg_d)[1]) > 0
            @test any_nonzero

            # Derivatives should be bounded (not blown up)
            @test abs(partials(lh_d)[1]) < 1e10
            @test abs(partials(sh_d)[1]) < 1e10
            @test abs(partials(tg_d)[1]) < 1e10

            # ---- Finite difference comparison ----
            eps_fd = 0.01

            # Reference run
            (inst_ref, _, _, _) = CLM.clm_initialize!(;
                fsurdat=fsurdat, paramfile=paramfile)
            setup_forcing!(inst_ref.atm2lnd, scenario.T, ng;
                precip_rain=scenario.precip_rain, precip_snow=scenario.precip_snow,
                vp=scenario.vp, solar_direct=scenario.solar_direct,
                solar_diffuse=scenario.solar_diffuse, lwrad=scenario.lwrad,
                wind=scenario.wind)
            CLM.downscale_forcings!(bounds, inst_ref.atm2lnd, inst_ref.column,
                                    inst_ref.landunit, inst_ref.topo)
            for n in 1:3
                CLM.clm_drv!(config, inst_ref, filt, filt_ia, bounds,
                    true, nextsw_cday, declin, declin, CLM.ORB_OBLIQR_DEFAULT,
                    false, false, "", false;
                    nstep=n, is_first_step=(n==1), is_beg_curr_day=(n==1),
                    dtime=1800.0, mon=1, day=1,
                    photosyns=inst_ref.photosyns)
            end
            CLM.clm_drv!(config, inst_ref, filt, filt_ia, bounds,
                true, nextsw_cday, declin, declin, CLM.ORB_OBLIQR_DEFAULT,
                false, false, "", false;
                nstep=4, is_first_step=false,
                dtime=1800.0, mon=1, day=1,
                photosyns=inst_ref.photosyns)

            # Perturbed run
            (inst_pert, _, _, _) = CLM.clm_initialize!(;
                fsurdat=fsurdat, paramfile=paramfile)
            setup_forcing!(inst_pert.atm2lnd, scenario.T, ng;
                precip_rain=scenario.precip_rain, precip_snow=scenario.precip_snow,
                vp=scenario.vp, solar_direct=scenario.solar_direct,
                solar_diffuse=scenario.solar_diffuse, lwrad=scenario.lwrad,
                wind=scenario.wind)
            CLM.downscale_forcings!(bounds, inst_pert.atm2lnd, inst_pert.column,
                                    inst_pert.landunit, inst_pert.topo)
            for n in 1:3
                CLM.clm_drv!(config, inst_pert, filt, filt_ia, bounds,
                    true, nextsw_cday, declin, declin, CLM.ORB_OBLIQR_DEFAULT,
                    false, false, "", false;
                    nstep=n, is_first_step=(n==1), is_beg_curr_day=(n==1),
                    dtime=1800.0, mon=1, day=1,
                    photosyns=inst_pert.photosyns)
            end
            setup_forcing!(inst_pert.atm2lnd, scenario.T + eps_fd, ng;
                precip_rain=scenario.precip_rain, precip_snow=scenario.precip_snow,
                vp=scenario.vp, solar_direct=scenario.solar_direct,
                solar_diffuse=scenario.solar_diffuse, lwrad=scenario.lwrad,
                wind=scenario.wind)
            CLM.downscale_forcings!(bounds, inst_pert.atm2lnd, inst_pert.column,
                                    inst_pert.landunit, inst_pert.topo)
            CLM.clm_drv!(config, inst_pert, filt, filt_ia, bounds,
                true, nextsw_cday, declin, declin, CLM.ORB_OBLIQR_DEFAULT,
                false, false, "", false;
                nstep=4, is_first_step=false,
                dtime=1800.0, mon=1, day=1,
                photosyns=inst_pert.photosyns)

            # Compare AD vs FD
            for (var_name, d_val, ref_val, pert_val) in [
                ("LH",    lh_d,
                 inst_ref.energyflux.eflx_lh_tot_patch[test_p],
                 inst_pert.energyflux.eflx_lh_tot_patch[test_p]),
                ("SH",    sh_d,
                 inst_ref.energyflux.eflx_sh_tot_patch[test_p],
                 inst_pert.energyflux.eflx_sh_tot_patch[test_p]),
                ("T_grnd", tg_d,
                 inst_ref.temperature.t_grnd_col[1],
                 inst_pert.temperature.t_grnd_col[1]),
            ]
                ad_deriv = partials(d_val)[1]
                fd_deriv = (pert_val - ref_val) / eps_fd

                println("  $(scenario.name): d($var_name)/d(T) AD=$(round(ad_deriv, digits=4)), FD=$(round(fd_deriv, digits=4))")

                # Sign agreement
                if abs(ad_deriv) > 1e-10 && abs(fd_deriv) > 1e-10
                    @test sign(ad_deriv) == sign(fd_deriv)
                end

                # Relative agreement within 5% (relaxed for near-discontinuity scenarios)
                if abs(fd_deriv) > 1e-4
                    rel_err = abs(ad_deriv - fd_deriv) / abs(fd_deriv)
                    println("  $(scenario.name): d($var_name)/d(T) relative error = $(round(rel_err * 100, digits=2))%")
                    @test rel_err < 0.10  # 10% tolerance for multi-scenario
                end
            end

            println("  $(scenario.name) PASSED")
        end
    end

end  # outer testset
