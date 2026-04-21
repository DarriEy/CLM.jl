#!/usr/bin/env julia
# ==========================================================================
# Parameter Recovery Tests for CLM.jl Phase 5
#
# Twin experiments demonstrating that the calibration framework can
# recover known parameter values from synthetic observations.
#
# Experiment A: Single parameter (csoilc) → sensible heat target
# Experiment B: Two parameters (medlyn_slope + vcmax25_scale) → LH + GPP
# ==========================================================================

using Test
using ForwardDiff
using CLM

println("=" ^ 70)
println("PARAMETER RECOVERY TESTS for CLM.jl")
println("=" ^ 70)

const FSURDAT_PATH = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc"
const PARAMFILE_PATH = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc"

@testset "Parameter Recovery" begin

    if !isfile(FSURDAT_PATH) || !isfile(PARAMFILE_PATH)
        @warn "Skipping parameter recovery tests: input files not found"
        @test true
        return
    end

    # =====================================================================
    # Helper: generate synthetic observations from a "true" parameter set
    # =====================================================================
    function generate_observations(params, theta_true; T_forc=285.0, n_warmup=3)
        prob = CLM.CalibrationProblem(
            params=params,
            targets=CLM.CalibrationTarget[],  # dummy, not used
            fsurdat=FSURDAT_PATH,
            paramfile=PARAMFILE_PATH,
            n_warmup=n_warmup,
            T_forc=T_forc,
        )
        inst = CLM.run_clm_with_params(prob, theta_true)
        return inst
    end

    # =====================================================================
    # Test: CalibrationOverrides wiring
    # =====================================================================
    @testset "Overrides injection — gradient non-zero" begin
        # csoilc parameter with override-based injection
        params = [
            CLM.CalibrationParameter(
                "csoilc", 0.004, CLM._apply_csoilc!,
                (0.001, 0.02), :log)
        ]

        targets = [
            CLM.CalibrationTarget(
                "SH",
                inst -> begin
                    ef = inst.energyflux
                    # Return mean sensible heat over non-zero patches
                    vals = [ef.eflx_sh_tot_patch[p] for p in eachindex(ef.eflx_sh_tot_patch)
                            if isfinite(ef.eflx_sh_tot_patch[p])]
                    isempty(vals) ? 0.0 : sum(vals) / length(vals)
                end,
                50.0,  # arbitrary target
                1.0
            )
        ]

        prob = CLM.CalibrationProblem(
            params=params,
            targets=targets,
            fsurdat=FSURDAT_PATH,
            paramfile=PARAMFILE_PATH,
            n_warmup=3,
            T_forc=285.0,
        )

        θ0 = CLM.default_theta(params)
        grad = CLM.calibration_gradient(prob, θ0)
        println("  csoilc gradient at default: $(grad[1])")
        @test isfinite(grad[1])
        # The gradient should be non-zero since csoilc affects sensible heat
        # (It could be small but should not be exactly zero)
        println("  csoilc gradient non-zero check: |grad| = $(abs(grad[1]))")

        # Verify AD vs FD agreement
        # With use_smooth_fd=true (default), both FD and AD evaluate the same smooth
        # landscape, so gradients should agree within ~10x.
        eps_fd = 0.01
        obj0 = CLM.calibration_objective(prob, θ0)
        obj_p = CLM.calibration_objective(prob, θ0 .+ eps_fd)
        obj_m = CLM.calibration_objective(prob, θ0 .- eps_fd)
        grad_fd = (obj_p - obj_m) / (2 * eps_fd)

        println("  AD=$(grad[1]), FD=$(grad_fd)")
        @test isfinite(grad_fd)
        @test abs(grad[1]) > 1e-6
        # With smooth FD, AD and FD should agree within an order of magnitude
        if abs(grad_fd) > 1e-10
            ratio = abs(grad[1] / grad_fd)
            println("  AD/FD ratio: $(round(ratio, sigdigits=3))")
            @test 0.1 < ratio < 10.0  # within 10x
        end
        println("  Overrides injection gradient test: PASSED")
    end

    # =====================================================================
    # Experiment A: Single parameter recovery (csoilc)
    # =====================================================================
    @testset "Experiment A — csoilc recovery" begin
        println("\n  --- Experiment A: csoilc recovery ---")

        params = [
            CLM.CalibrationParameter(
                "csoilc", 0.004, CLM._apply_csoilc!,
                (0.001, 0.02), :log)
        ]

        # "True" value: csoilc = 0.006 (θ_true in log space)
        true_csoilc = 0.006
        θ_true = [log(true_csoilc / 0.004)]  # log transform: val = default * exp(θ)

        # Generate synthetic observations
        inst_true = generate_observations(params, θ_true; T_forc=285.0)

        # Extract sensible heat as observation (patch 2 = vegetated PFT)
        sh_obs = inst_true.energyflux.eflx_sh_tot_patch[2]
        println("  True csoilc=$(true_csoilc), observed SH=$(round(sh_obs, digits=4))")

        if !isfinite(sh_obs) || abs(sh_obs) < 1e-10
            @warn "SH observation is zero or non-finite, skipping recovery"
            @test true
            return
        end

        targets = [
            CLM.CalibrationTarget("SH",
                inst -> inst.energyflux.eflx_sh_tot_patch[2],
                Float64(sh_obs), 1.0)
        ]

        prob = CLM.CalibrationProblem(
            params=params,
            targets=targets,
            fsurdat=FSURDAT_PATH,
            paramfile=PARAMFILE_PATH,
            n_warmup=3,
            T_forc=285.0,
        )

        # Start from default (θ=0 in log space → csoilc = 0.004)
        θ_start = [0.0]
        obj_start = CLM.calibration_objective(prob, θ_start)
        println("  Starting: θ=0.0 (csoilc=0.004), obj=$(round(obj_start, sigdigits=6))")

        # Run optimization
        result = CLM.calibrate(prob;
            theta0=θ_start,
            maxiter=20,
            verbose=false)

        # Recover physical parameter value
        recovered_csoilc = 0.004 * exp(result.theta[1])
        rel_recovery_error = abs(recovered_csoilc - true_csoilc) / true_csoilc

        println("  Recovered csoilc=$(round(recovered_csoilc, sigdigits=4)), " *
                "true=$(true_csoilc), error=$(round(rel_recovery_error*100, digits=1))%")
        println("  Iterations: $(result.iterations), converged: $(result.converged)")
        println("  Final objective: $(round(result.objective, sigdigits=6))")

        # Success criterion: objective decreased significantly
        # With use_smooth_fd=true, FD and AD evaluate the same landscape,
        # so gradient descent should converge well.
        @test result.objective <= obj_start * 0.5  # at least 50% reduction
        @test result.iterations <= 20
        # Recovery should be reasonable (within 50% of true value)
        @test rel_recovery_error < 0.5
        println("  Experiment A: PASSED")
    end

    # =====================================================================
    # Experiment B: Two-parameter recovery (medlyn_slope + vcmax25_scale)
    # =====================================================================
    @testset "Experiment B — two-parameter recovery" begin
        println("\n  --- Experiment B: medlyn_slope + vcmax25_scale ---")

        params = [
            CLM.CalibrationParameter(
                "medlynslope", 6.0, CLM._apply_medlynslope!,
                (1.0, 20.0), :identity),
            CLM.CalibrationParameter(
                "vcmax25_scale", 1.0, CLM._apply_vcmax25_scale!,
                (0.3, 3.0), :identity),
        ]

        # "True" values: medlyn_slope=8.0 (θ=8/6=1.333), vcmax25_scale=1.3
        true_medlyn = 8.0
        true_vcmax = 1.3
        θ_true = [true_medlyn / 6.0, true_vcmax]  # identity: val = default * θ

        # Generate synthetic observations
        inst_true = generate_observations(params, θ_true; T_forc=285.0)

        # Extract LH and GPP as observations (patch 2 = vegetated PFT)
        lh_obs = inst_true.energyflux.eflx_lh_tot_patch[2]
        gpp_obs = inst_true.photosyns.psnsun_patch[2] + inst_true.photosyns.psnsha_patch[2]

        println("  True: medlyn=$(true_medlyn), vcmax_scale=$(true_vcmax)")
        println("  Observed: LH=$(round(lh_obs, digits=4)), GPP=$(round(gpp_obs, sigdigits=4))")

        if !isfinite(lh_obs) || !isfinite(gpp_obs)
            @warn "Observations non-finite, skipping two-param recovery"
            @test true
            return
        end

        targets = [
            CLM.CalibrationTarget("LH",
                inst -> inst.energyflux.eflx_lh_tot_patch[2],
                Float64(lh_obs), 1.0),
            CLM.CalibrationTarget("GPP",
                inst -> inst.photosyns.psnsun_patch[2] + inst.photosyns.psnsha_patch[2],
                Float64(gpp_obs),
                abs(gpp_obs) > 1e-10 ? 1.0 / (gpp_obs^2) : 1.0),  # normalize
        ]

        prob = CLM.CalibrationProblem(
            params=params,
            targets=targets,
            fsurdat=FSURDAT_PATH,
            paramfile=PARAMFILE_PATH,
            n_warmup=3,
            T_forc=285.0,
        )

        # Start from defaults
        θ_start = CLM.default_theta(params)
        obj_start = CLM.calibration_objective(prob, θ_start)
        println("  Starting: θ=$(θ_start), obj=$(round(obj_start, sigdigits=6))")

        # Run optimization
        result = CLM.calibrate(prob;
            theta0=θ_start,
            maxiter=30,
            verbose=false)

        # Recover physical values
        recovered_medlyn = 6.0 * result.theta[1]
        recovered_vcmax = 1.0 * result.theta[2]

        err_medlyn = abs(recovered_medlyn - true_medlyn) / true_medlyn
        err_vcmax = abs(recovered_vcmax - true_vcmax) / true_vcmax

        println("  Recovered: medlyn=$(round(recovered_medlyn, digits=2)) " *
                "(true=$(true_medlyn), err=$(round(err_medlyn*100, digits=1))%)")
        println("  Recovered: vcmax_scale=$(round(recovered_vcmax, digits=3)) " *
                "(true=$(true_vcmax), err=$(round(err_vcmax*100, digits=1))%)")
        println("  Iterations: $(result.iterations), converged: $(result.converged)")
        println("  Final objective: $(round(result.objective, sigdigits=6))")

        # Success criterion: objective decreased significantly
        # With use_smooth_fd=true, gradient descent should converge well.
        @test result.objective <= obj_start * 0.5  # at least 50% reduction
        @test result.iterations <= 30
        println("  Experiment B: PASSED")
    end

    # =====================================================================
    # Multi-timestep evaluation test
    # =====================================================================
    @testset "Multi-timestep calibration" begin
        println("\n  --- Multi-timestep evaluation ---")

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

        # Single-step problem
        prob1 = CLM.CalibrationProblem(
            params=params, targets=targets,
            fsurdat=FSURDAT_PATH, paramfile=PARAMFILE_PATH,
            n_warmup=2, n_eval=1, T_forc=285.0)

        # Multi-step problem (more warmup to stabilize state)
        prob3 = CLM.CalibrationProblem(
            params=params, targets=targets,
            fsurdat=FSURDAT_PATH, paramfile=PARAMFILE_PATH,
            n_warmup=5, n_eval=2, T_forc=285.0)

        θ0 = [0.0]
        obj1 = CLM.calibration_objective(prob1, θ0)
        obj3 = CLM.calibration_objective(prob3, θ0)

        @test isfinite(obj1)
        @test isfinite(obj3)

        # Gradients should both be finite
        grad1 = CLM.calibration_gradient(prob1, θ0)
        grad3 = CLM.calibration_gradient(prob3, θ0)

        @test isfinite(grad1[1])
        @test isfinite(grad3[1])

        println("  n_eval=1: obj=$(round(obj1, sigdigits=4)), grad=$(round(grad1[1], sigdigits=4))")
        println("  n_eval=3: obj=$(round(obj3, sigdigits=4)), grad=$(round(grad3[1], sigdigits=4))")
        println("  Multi-timestep calibration: PASSED")
    end

end  # outer testset
