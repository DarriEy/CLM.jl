using Test
using CLM

# ============================================================================
# FUN-path crop N parity -- soybean N fixation in a COUPLED crop column with
# use_fun = .true. (ctsm5.3.012, use_cn + use_crop + use_fertilizer +
# use_flexiblecn + use_fun), against per-step dumps of the instrumented
# crop_fun_macos case (scripts/validation/crop_parity/build_crop_fun_case.sh;
# recipe: docs/CROP_PARITY.md § "FUN crop N").
#
# With FUN, CNSoyfix is never called (CNDriverMod.F90:324): soybean fixation is
# CNFUN's symbiotic Nfix (FUN_fracfixers = 1 for soybean). Every check injects the
# Fortran inputs of a step into the Julia routine and compares with the Fortran
# outputs of the same step (scripts/fortran_parity_crop_fun.jl).
#
# DATA-GATED: skipped when $SYMFLUENCE_DATA/runs/crop_fun_ref (or CROP_FUN_REF)
# holds no dumps.
# ============================================================================

ENV["CROP_FUN_NO_MAIN"] = "1"
ENV["CROP_TAILS_NO_MAIN"] = "1"
include(joinpath(@__DIR__, "..", "scripts", "fortran_parity_crop_fun.jl"))

const _CF_HAVE = isdir(CF_REF) && !isempty(cf_steps(CF_REF)) && isfile(CT_PARAMFILE)

@testset "FUN crop N vs ctsm5.3.012 (soybean fixation, use_fun=.true.)" begin
    if !_CF_HAVE
        @test_skip "crop FUN reference dumps not found at $CF_REF"
    else
        pc = ct_pftcon()

        @testset "cnfun! (incl. soybean Nfix) per step" begin
            r = cf_check_cnfun(CF_REF; pc)
            @test r.nsteps >= 300
            # soybean genuinely fixes N in the reference (non-vacuous)
            @test r.soy.nfix_pos >= 300
            @test r.soy.nfix_max > 1e-7
            for (k, s) in r.stats
                @test s.n > 0
                # ULP-level only, relative to the field's scale (the harness's
                # dzsoi_decomp is recomputed, not dumped). The COST_* diagnostics
                # are ratios N/npp of near-zero fluxes on some patches, which
                # amplify ULP noise (worst 5e-7 of scale).
                @test cf_normerr(s) < (startswith(k, "COST_") ? 1e-5 : 1e-12)
            end
            for k in ("NFIX", "NUPTAKE", "SMINN_TO_PLANT_FUN", "NNONMYC", "NNONMYC_NH4", "NPP_GROWTH",
                      "NPP_BURNEDOFF", "NFIX_TO_SMINN", "RETRANSN_TO_NPOOL", "FREE_RETRANSN_TO_NPOOL")
                @test cf_flips(r.stats[k]) == 0
            end
        end

        @testset "FUN reads SOIL-layer temperature (not the snow layers)" begin
            # The driver used to hand cnfun! the full t_soisno_col, whose first
            # nlevsno columns are SNOW layers; fun_cost_fix then saw SPVAL and
            # soybean fixation collapsed. Reproduce both layouts on a live step.
            n = first(n for n in cf_steps(CF_REF) if begin
                a = cf_load(CF_REF, "after_cnfun", n)
                any(a["NFIX"][p] > 1e-8 for p in eachindex(a["ITYPE"]) if a["ITYPE"][p] in CF_SOY)
            end)
            b = cf_load(CF_REF, "before_cnfun", n); a = cf_load(CF_REF, "after_cnfun", n)
            psoy = first(p for p in eachindex(b["ITYPE"]) if b["ITYPE"][p] in CF_SOY && a["NFIX"][p] > 1e-8)
            good = cf_run_cnfun(b; pc)
            @test good.nf.Nfix_patch[psoy] ≈ a["NFIX"][psoy] rtol = 1e-12
            nsno = 12
            padded = hcat(fill(CLM.SPVAL, size(b["T_SOISNO"], 1), nsno), b["T_SOISNO"])
            bad = copy(b); bad["T_SOISNO"] = padded[:, 1:size(b["T_SOISNO"], 2)]
            wrong = cf_run_cnfun(bad; pc)
            @test wrong.nf.Nfix_patch[psoy] < 1e-6 * a["NFIX"][psoy]   # the old wiring's signature
        end

        @testset "FlexibleCN N demand, prognostic-crop call (FUN)" begin
            r = cf_check_ndemand(CF_REF; pc)
            if r.available
                @test r.ngrainfill_onset >= 1        # grain-fill retranslocation exercised
                for (k, s) in r.stats
                    @test s.n > 0
                    @test cf_normerr(s) < 1e-12
                end
            else
                @test_skip "dumps predate GRAIN_FLAG"
            end
        end

        @testset "FlexibleCN allocation under FUN (sminn_to_npool, npp_growth)" begin
            r = cf_check_alloc(CF_REF; pc)
            if r.available
                @test r.nsoy > 100
                for (k, s) in r.stats
                    @test s.n > 0
                    @test cf_normerr(s) < 1e-12
                end
            else
                @test_skip "dumps predate the FUN allocation fields"
            end
        end

        @testset "npool / retransn update" begin
            r = cf_check_npool(CF_REF)
            if r.available
                for (k, s) in r.stats
                    @test s.n > 0
                    @test cf_normerr(s) < 1e-12
                end
            else
                @test_skip "dumps predate the FUN allocation fields"
            end
        end

        @testset "SOYFIXN and FPG under FUN" begin
            r = cf_check_soyfixn_fpg(CF_REF)
            @test r.nsoy > 1000
            @test r.soyfixn_nonzero == 0        # CNSoyfix never runs with FUN
        end

        @testset "harvest step of CNOffsetLitterfall with FUN" begin
            r = cf_check_harvest(CF_REF; pc)
            @test length(r.events) >= 1
            for (k, s) in r.stats
                @test s.n > 0
                @test cf_normerr(s) < 1e-12
            end
        end
    end
end
