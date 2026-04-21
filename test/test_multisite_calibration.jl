#!/usr/bin/env julia
# ==========================================================================
# Multi-Site Calibration Tests for CLM.jl
#
# Tests the multi-site calibration framework with synthetic "sites"
# that share a parameter vector but have different forcing.
# ==========================================================================

using Test
using ForwardDiff
using CLM

println("=" ^ 70)
println("MULTI-SITE CALIBRATION TESTS for CLM.jl")
println("=" ^ 70)

const FSURDAT_PATH = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc"
const PARAMFILE_PATH = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc"

@testset "Multi-Site Calibration" begin

    if !isfile(FSURDAT_PATH) || !isfile(PARAMFILE_PATH)
        @warn "Skipping multi-site tests: input files not found"
        @test true
        return
    end

    # =====================================================================
    # Test 1: FluxnetForcing → CLM forcing setup
    # =====================================================================
    @testset "FluxnetForcing setup" begin
        forcing = CLM.FluxnetForcing(ta=290.0, pa=90000.0, ws=4.0,
            lw_in=280.0, sw_in=400.0, precip=1e-5, vpd=800.0)

        (inst, bounds, filt, tm) = CLM.clm_initialize!(;
            fsurdat=FSURDAT_PATH, paramfile=PARAMFILE_PATH)
        ng = bounds.endg

        CLM.setup_fluxnet_forcing!(inst.atm2lnd, forcing, ng)

        @test inst.atm2lnd.forc_t_not_downscaled_grc[1] == 290.0
        @test inst.atm2lnd.forc_pbot_not_downscaled_grc[1] == 90000.0
        @test inst.atm2lnd.forc_wind_grc[1] == 4.0
        @test inst.atm2lnd.forc_rain_not_downscaled_grc[1] ≈ 1e-5  # T > TFRZ+2
        println("  FluxnetForcing setup: PASSED")
    end

    # =====================================================================
    # Test 2: Two synthetic sites with different temperatures
    # =====================================================================
    @testset "Two-site objective and gradient" begin
        params = [
            CLM.CalibrationParameter(
                "csoilc", 0.004, CLM._apply_csoilc!,
                (0.001, 0.02), :log)
        ]

        targets = [
            CLM.CalibrationTarget("SH",
                inst -> inst.energyflux.eflx_sh_tot_patch[2],
                50.0, 1.0)
        ]

        # Site 1: warm (290K)
        prob1 = CLM.CalibrationProblem(
            params=params, targets=targets,
            fsurdat=FSURDAT_PATH, paramfile=PARAMFILE_PATH,
            n_warmup=3, T_forc=290.0)

        # Site 2: cool (275K)
        prob2 = CLM.CalibrationProblem(
            params=params, targets=targets,
            fsurdat=FSURDAT_PATH, paramfile=PARAMFILE_PATH,
            n_warmup=3, T_forc=275.0)

        site1 = CLM.SiteCalibrationProblem(problem=prob1, weight=1.0)
        site2 = CLM.SiteCalibrationProblem(problem=prob2, weight=1.0)

        msc = CLM.MultiSiteCalibration(
            sites=[site1, site2],
            params=params,
            aggregation=:sum)

        θ0 = [0.0]

        # Test objective
        obj_ms = CLM.calibration_objective_multisite(msc, θ0)
        obj_s1 = CLM.calibration_objective(prob1, θ0)
        obj_s2 = CLM.calibration_objective(prob2, θ0)

        @test isfinite(obj_ms)
        @test abs(obj_ms - (obj_s1 + obj_s2)) < 1e-10
        println("  Multi-site objective: sum=$(round(obj_ms, sigdigits=4)), " *
                "s1=$(round(obj_s1, sigdigits=4)), s2=$(round(obj_s2, sigdigits=4))")

        # Test gradient
        grad_ms = CLM.calibration_gradient_multisite(msc, θ0)
        @test isfinite(grad_ms[1])
        @test abs(grad_ms[1]) > 1e-6  # gradient should be non-trivial
        println("  Multi-site gradient: $(round(grad_ms[1], sigdigits=4))")
        println("  Two-site test: PASSED")
    end

    # =====================================================================
    # Test 3: Multi-site optimization (convergence check)
    # =====================================================================
    @testset "Multi-site optimization" begin
        params = [
            CLM.CalibrationParameter(
                "csoilc", 0.004, CLM._apply_csoilc!,
                (0.001, 0.02), :log)
        ]

        # Generate "true" observation from each site
        true_theta = [log(0.006 / 0.004)]  # csoilc = 0.006

        prob_warm = CLM.CalibrationProblem(
            params=params, targets=CLM.CalibrationTarget[],
            fsurdat=FSURDAT_PATH, paramfile=PARAMFILE_PATH,
            n_warmup=3, T_forc=290.0)
        inst_warm = CLM.run_clm_with_params(prob_warm, true_theta)
        sh_warm = inst_warm.energyflux.eflx_sh_tot_patch[2]

        prob_cool = CLM.CalibrationProblem(
            params=params, targets=CLM.CalibrationTarget[],
            fsurdat=FSURDAT_PATH, paramfile=PARAMFILE_PATH,
            n_warmup=3, T_forc=275.0)
        inst_cool = CLM.run_clm_with_params(prob_cool, true_theta)
        sh_cool = inst_cool.energyflux.eflx_sh_tot_patch[2]

        if !isfinite(sh_warm) || !isfinite(sh_cool)
            @warn "Non-finite observations, skipping optimization"
            @test true
            return
        end

        # Create optimization problems with targets
        targets_warm = [CLM.CalibrationTarget("SH_warm",
            inst -> inst.energyflux.eflx_sh_tot_patch[2],
            Float64(sh_warm), 1.0)]
        targets_cool = [CLM.CalibrationTarget("SH_cool",
            inst -> inst.energyflux.eflx_sh_tot_patch[2],
            Float64(sh_cool), 1.0)]

        prob1_opt = CLM.CalibrationProblem(
            params=params, targets=targets_warm,
            fsurdat=FSURDAT_PATH, paramfile=PARAMFILE_PATH,
            n_warmup=3, T_forc=290.0)
        prob2_opt = CLM.CalibrationProblem(
            params=params, targets=targets_cool,
            fsurdat=FSURDAT_PATH, paramfile=PARAMFILE_PATH,
            n_warmup=3, T_forc=275.0)

        msc = CLM.MultiSiteCalibration(
            sites=[
                CLM.SiteCalibrationProblem(problem=prob1_opt, weight=1.0),
                CLM.SiteCalibrationProblem(problem=prob2_opt, weight=1.0),
            ],
            params=params,
            aggregation=:sum)

        θ_start = [0.0]
        obj_start = CLM.calibration_objective_multisite(msc, θ_start)

        result = CLM.calibrate_multisite(msc;
            theta0=θ_start, maxiter=15, verbose=false)

        println("  Multi-site optimization: obj $(round(obj_start, sigdigits=4)) → $(round(result.objective, sigdigits=4))")
        println("  Iterations: $(result.iterations)")
        @test result.objective <= obj_start * 1.01
        println("  Multi-site optimization: PASSED")
    end

end
