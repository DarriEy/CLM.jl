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

include(joinpath(@__DIR__, "testdata.jl"))

@testset "AD Robustness — Multi-scenario" begin

    # Data root is env-overridable (see test/testdata.jl) so this test still RUNS
    # off the Drive migration copy once the original Calgary-laptop path is gone.
    # Without this, the skip below turned the whole testset into a permanent
    # `@test true` — vacuously green forever.
    fsurdat, paramfile = bow_params()

    if !isfile(fsurdat) || !isfile(paramfile)
        testdata_missing("AD robustness tests", fsurdat, paramfile) && return
    end

    # Runs on the global default cold start (Fortran-matching). No pin is needed: the FD
    # is evaluated on the SAME smoothed physics as the AD (Dual values, see build_run_dual),
    # so the freeze/thaw phase change no longer makes a hard-FD reference disagree with the
    # smooth-AD derivative — they match to FD truncation error in every scenario.

    # --- Helper: create a Dual-typed copy of a parameterized struct ---
    # Delegate to the PRODUCTION routine rather than keeping a private copy of it.
    # This file used to carry a duplicate that drifted: when #238 added
    # CH4FInundatedStream{V<:AbstractVector} to the instance tree, `wrapper{D}()`
    # started throwing `TypeError: in V, expected V<:AbstractVector, got Dual`, and
    # only calibration.jl's copy got the by-ref guard — so this whole testset errored
    # while the production path was fine. One implementation, one guard.
    make_dual_copy(src, ::Type{D}) where {D} = CLM._calib_dual_copy(src, D)

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
            a2l.forc_hgt_u_grc[g] = 30.0
            a2l.forc_hgt_t_grc[g] = 30.0
            a2l.forc_hgt_q_grc[g] = 30.0
            a2l.forc_topo_grc[g] = 0.0
            a2l.forc_wind_grc[g] = wind
            a2l.forc_u_grc[g] = wind
            a2l.forc_v_grc[g] = 0.0
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

            # ---- Build a Dual-typed inst from the warmed-up Float64 state, seed the forcing
            #      temperature with `seedp`, and run one step. AD seeds partial 1; the FD
            #      below reuses this with partial 0. This is deliberate: the Dual path always
            #      evaluates the SMOOTHED physics (the smooth_* primitives dispatch on AD
            #      element types — _use_smooth is purely type-based, so a plain Float64 re-run
            #      evaluates the HARD physics and would diverge from AD by the smoothing offset
            #      near the phase-change front, e.g. 284.3 K hard vs 283.0 K smooth at 295 K).
            #      Evaluating the FD on Dual VALUES makes FD and AD the SAME smoothed function,
            #      so the check validates the AD derivative, not the smooth-vs-hard gap. ----
            D = Dual{Nothing, Float64, 1}
            function build_run_dual(Tseed, seedp)
                inst_d = CLM.CLMInstances()
                for name in fieldnames(CLM.CLMInstances)
                    # dgv_ecophyscon is parametric on its *vector* type (not a Real),
                    # so make_dual_copy's wrapper{D}() rebuild does not apply — share
                    # it by reference like surfdata.
                    (name === :water || name === :surfdata ||
                     name === :dgv_ecophyscon) && continue
                    setfield!(inst_d, name, make_dual_copy(getfield(inst_f64, name), D))
                end
                inst_d.surfdata = inst_f64.surfdata
                inst_d.dgv_ecophyscon = inst_f64.dgv_ecophyscon
                inst_d.photosyns.stomatalcond_mtd = inst_f64.photosyns.stomatalcond_mtd
                water_d = CLM.WaterData()
                for name in fieldnames(CLM.WaterData)
                    try setfield!(water_d, name, getfield(inst_f64.water, name)) catch end
                end
                wsb_d = make_dual_copy(inst_f64.water.waterstatebulk_inst, D)
                wsb_d.ws = make_dual_copy(inst_f64.water.waterstatebulk_inst.ws, D)
                water_d.waterstatebulk_inst = wsb_d
                wfb_d = make_dual_copy(inst_f64.water.waterfluxbulk_inst, D)
                wfb_d.wf = make_dual_copy(inst_f64.water.waterfluxbulk_inst.wf, D)
                water_d.waterfluxbulk_inst = wfb_d
                water_d.waterdiagnosticbulk_inst = make_dual_copy(inst_f64.water.waterdiagnosticbulk_inst, D)
                water_d.waterbalancebulk_inst = make_dual_copy(inst_f64.water.waterbalancebulk_inst, D)
                if !isempty(water_d.bulk_and_tracers)
                    bt = water_d.bulk_and_tracers[water_d.i_bulk]
                    bt.waterflux = wfb_d; bt.waterstate = wsb_d
                    bt.waterdiagnostic = water_d.waterdiagnosticbulk_inst
                    bt.waterbalance = water_d.waterbalancebulk_inst
                end
                inst_d.water = water_d
                setup_forcing!(inst_d.atm2lnd, Tseed, ng;
                    precip_rain=scenario.precip_rain, precip_snow=scenario.precip_snow,
                    vp=scenario.vp, solar_direct=scenario.solar_direct,
                    solar_diffuse=scenario.solar_diffuse, lwrad=scenario.lwrad, wind=scenario.wind)
                for g in 1:ng
                    T_d = Dual(Tseed, seedp)
                    inst_d.atm2lnd.forc_t_not_downscaled_grc[g] = T_d
                    inst_d.atm2lnd.forc_th_not_downscaled_grc[g] = T_d * (D(100000.0) / D(85000.0))^(D(CLM.RAIR) / D(CLM.CPAIR))
                    inst_d.atm2lnd.forc_rho_not_downscaled_grc[g] = D(85000.0) / (D(CLM.RAIR) * T_d)
                end
                CLM.downscale_forcings!(bounds, inst_d.atm2lnd, inst_d.column, inst_d.landunit, inst_d.topo)
                CLM.clm_drv!(config, inst_d, filt, filt_ia, bounds,
                    true, nextsw_cday, declin, declin, CLM.ORB_OBLIQR_DEFAULT,
                    false, false, "", false;
                    nstep=4, is_first_step=false, dtime=1800.0, mon=1, day=1,
                    photosyns=inst_d.photosyns)
                return inst_d
            end

            inst_d = build_run_dual(scenario.T, 1.0)

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

            # ---- Finite difference on the SAME smoothed physics (Dual VALUE, partial 0) ----
            # Central difference of the Dual-evaluated function at T ± eps. Because both the
            # AD derivative (partials at T) and this FD secant come from the identical smoothed
            # physics, they agree to the FD truncation error — the check now validates the AD
            # derivative itself rather than the smooth-vs-hard offset. (Re-running in Float64
            # would evaluate the HARD physics and disagree near the phase-change front; see the
            # build_run_dual note above.)
            eps_fd = 0.01
            inst_dp = build_run_dual(scenario.T + eps_fd, 0.0)
            inst_dm = build_run_dual(scenario.T - eps_fd, 0.0)

            for (var_name, d_val, plus_val, minus_val) in [
                ("LH",    lh_d,
                 inst_dp.energyflux.eflx_lh_tot_patch[test_p],
                 inst_dm.energyflux.eflx_lh_tot_patch[test_p]),
                ("SH",    sh_d,
                 inst_dp.energyflux.eflx_sh_tot_patch[test_p],
                 inst_dm.energyflux.eflx_sh_tot_patch[test_p]),
                ("T_grnd", tg_d,
                 inst_dp.temperature.t_grnd_col[1],
                 inst_dm.temperature.t_grnd_col[1]),
            ]
                ad_deriv = partials(d_val)[1]
                fd_deriv = (value(plus_val) - value(minus_val)) / (2 * eps_fd)

                println("  $(scenario.name): d($var_name)/d(T) AD=$(round(ad_deriv, digits=4)), FD=$(round(fd_deriv, digits=4))")

                # Sign agreement
                if abs(ad_deriv) > 1e-10 && abs(fd_deriv) > 1e-10
                    @test sign(ad_deriv) == sign(fd_deriv)
                end

                # Relative agreement: AD and the smoothed-physics FD evaluate the same
                # function, so they agree tightly; a small allowance covers FD truncation
                # (O(eps^2)) and any residual near-freezing curvature.
                if abs(fd_deriv) > 1e-4
                    rel_err = abs(ad_deriv - fd_deriv) / abs(fd_deriv)
                    println("  $(scenario.name): d($var_name)/d(T) relative error = $(round(rel_err * 100, digits=2))%")
                    near_freeze = abs(scenario.T - CLM.TFRZ) < 5.0
                    tol = near_freeze ? 0.35 : 0.10
                    @test rel_err < tol
                end
            end

            println("  $(scenario.name) PASSED")
        end
    end
end  # outer testset
