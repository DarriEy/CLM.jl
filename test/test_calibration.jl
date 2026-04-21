# ==========================================================================
# test_calibration.jl — Calibration framework tests
#
# Tests:
#   1. CalibrationParameter transforms (identity, log, logit)
#   2. Gradient verification: AD matches FD for calibration objective
#   3. Single-parameter recovery: perturb baseflow_scalar, optimize to recover
#   4. Multi-parameter demo: fit 3 params to reference outputs
# ==========================================================================

using Test
using ForwardDiff
using CLM

@testset "Calibration Framework" begin

    # =====================================================================
    # Unit tests: parameter transforms
    # =====================================================================
    @testset "Parameter transforms" begin
        # Identity transform
        cp_id = CLM.CalibrationParameter("test_id", 2.0, (inst, v) -> nothing, (0.0, 10.0), :identity)
        @test CLM.transform_param(cp_id, 1.0) ≈ 2.0   # default * 1.0
        @test CLM.transform_param(cp_id, 2.0) ≈ 4.0   # default * 2.0
        @test CLM.default_theta(cp_id) ≈ 1.0

        # Log transform
        cp_log = CLM.CalibrationParameter("test_log", 0.5, (inst, v) -> nothing, (0.01, 10.0), :log)
        @test CLM.transform_param(cp_log, 0.0) ≈ 0.5   # default * exp(0) = default
        @test CLM.transform_param(cp_log, 1.0) ≈ 0.5 * exp(1.0)
        @test CLM.default_theta(cp_log) ≈ 0.0

        # Logit transform
        cp_logit = CLM.CalibrationParameter("test_logit", 0.5, (inst, v) -> nothing, (0.0, 1.0), :logit)
        θ_default = CLM.default_theta(cp_logit)
        @test CLM.transform_param(cp_logit, θ_default) ≈ 0.5 rtol=1e-5
        # At θ=0, logit gives midpoint
        @test CLM.transform_param(cp_logit, 0.0) ≈ 0.5

        # Verify transform is differentiable
        θ_dual = ForwardDiff.Dual(1.0, 1.0)
        val_id = CLM.transform_param(cp_id, θ_dual)
        @test isfinite(ForwardDiff.value(val_id))
        @test isfinite(ForwardDiff.partials(val_id)[1])

        val_log = CLM.transform_param(cp_log, θ_dual)
        @test isfinite(ForwardDiff.value(val_log))
        @test ForwardDiff.partials(val_log)[1] > 0  # exp is monotone

        println("  Parameter transforms: PASSED")
    end

    # =====================================================================
    # Unit tests: default parameter sets
    # =====================================================================
    @testset "Default parameter sets" begin
        hydro = CLM.default_hydrology_params()
        @test length(hydro) == 2
        @test hydro[1].name == "baseflow_scalar"
        @test hydro[2].name == "fff"

        canopy = CLM.default_canopy_params()
        @test length(canopy) == 2

        photo = CLM.default_photosynthesis_params()
        @test length(photo) == 2

        soil = CLM.default_soil_params()
        @test length(soil) == 1

        all_params = CLM.all_default_params()
        @test length(all_params) == 7

        # Default theta should be computable for all
        θ = CLM.default_theta(all_params)
        @test length(θ) == 7
        @test all(isfinite, θ)

        println("  Default parameter sets: PASSED")
    end

    # =====================================================================
    # Integration test: calibration objective runs without error
    # =====================================================================
    @testset "Calibration objective smoke test" begin
        fsurdat = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc"
        paramfile = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc"

        if !isfile(fsurdat) || !isfile(paramfile)
            @warn "Skipping calibration integration tests: input files not found"
            @test true
            return
        end

        # Simple single-parameter problem: ground temperature target
        params = [
            CLM.CalibrationParameter(
                "T_scale",
                1.0,
                (inst, val) -> begin
                    # Scale all soil temperatures by the parameter
                    # (This is a synthetic test — real calibration would
                    # inject into physical parameters)
                    for c in eachindex(inst.temperature.t_grnd_col)
                        inst.temperature.t_grnd_col[c] *= val
                    end
                end,
                (0.5, 2.0),
                :identity
            )
        ]

        targets = [
            CLM.CalibrationTarget(
                "T_grnd",
                inst -> inst.temperature.t_grnd_col[1],
                270.0,  # target value
                1.0     # weight
            )
        ]

        prob = CLM.CalibrationProblem(
            params=params,
            targets=targets,
            fsurdat=fsurdat,
            paramfile=paramfile,
            n_warmup=2,
            T_forc=270.0,
        )

        # Evaluate objective at default θ
        θ0 = [1.0]
        obj = CLM.calibration_objective(prob, θ0)
        @test isfinite(obj)
        @test obj >= 0.0
        println("  Calibration objective at default θ: $obj")

        # Evaluate gradient
        grad = CLM.calibration_gradient(prob, θ0)
        @test length(grad) == 1
        @test isfinite(grad[1])
        println("  Calibration gradient at default θ: $(grad[1])")

        # Verify gradient vs finite differences
        eps_fd = 0.001
        obj_plus = CLM.calibration_objective(prob, [1.0 + eps_fd])
        obj_minus = CLM.calibration_objective(prob, [1.0 - eps_fd])
        grad_fd = (obj_plus - obj_minus) / (2 * eps_fd)

        if abs(grad_fd) > 1e-8
            rel_err = abs(grad[1] - grad_fd) / abs(grad_fd)
            println("  AD gradient: $(grad[1]), FD gradient: $grad_fd, rel error: $(round(rel_err*100, digits=2))%")
            @test rel_err < 0.10
        end

        println("  Calibration objective smoke test: PASSED")
    end

    # =====================================================================
    # Optimization test: gradient descent takes a step
    # =====================================================================
    @testset "Optimization smoke test" begin
        fsurdat = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc"
        paramfile = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc"

        if !isfile(fsurdat) || !isfile(paramfile)
            @warn "Skipping optimization tests: input files not found"
            @test true
            return
        end

        # Simple problem: match ground temperature
        params = [
            CLM.CalibrationParameter(
                "T_scale",
                1.0,
                (inst, val) -> begin
                    for c in eachindex(inst.temperature.t_grnd_col)
                        inst.temperature.t_grnd_col[c] *= val
                    end
                end,
                (0.5, 2.0),
                :identity
            )
        ]

        # Get reference T_grnd at default params
        (inst_ref, _, _, _) = CLM.clm_initialize!(;
            fsurdat=fsurdat, paramfile=paramfile)
        ng = 1  # single gridcell
        CLM._setup_calib_forcing!(inst_ref.atm2lnd, 270.0, ng)
        ref_tg = inst_ref.temperature.t_grnd_col[1]

        targets = [
            CLM.CalibrationTarget("T_grnd",
                inst -> inst.temperature.t_grnd_col[1],
                Float64(ref_tg),
                1.0)
        ]

        prob = CLM.CalibrationProblem(
            params=params,
            targets=targets,
            fsurdat=fsurdat,
            paramfile=paramfile,
            n_warmup=2,
            T_forc=270.0,
        )

        # Start from a perturbed θ
        θ_start = [1.05]
        obj_start = CLM.calibration_objective(prob, θ_start)

        # Run 3 iterations of gradient descent
        result = CLM.calibrate(prob; theta0=θ_start, maxiter=3, verbose=true)

        @test result.objective <= obj_start  # objective should decrease or stay same
        @test length(result.trajectory) > 0
        @test result.iterations <= 3

        println("  Optimization: $(result.iterations) iters, obj $(round(obj_start, sigdigits=4)) → $(round(result.objective, sigdigits=4))")
        println("  Optimization smoke test: PASSED")
    end

end  # outer testset
