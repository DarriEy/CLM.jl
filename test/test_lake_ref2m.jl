# ==========================================================================
# test_lake_ref2m.jl — the two PORT defects behind the lake surface-flux
# residual (docs/LAKE_FLUX_RESIDUAL.md findings B and C).
#
# B. `t_ref2m_patch` was NEVER written on a lake patch. Only
#    bareground_fluxes.jl and urban_fluxes.jl touched it, and neither runs on a
#    lake, so TSA on a 100%-lake gridcell reported the untouched cold-start
#    283.000 K at every one of 48 steps while Fortran ranged over 255-270 K.
#    A DEAD WRITE reported as a "7-11% physics residual".
#
# C. The lake Monin-Obukhov kernel implemented ONE of CTSM's FOUR `zeta`
#    stability regimes (the plain unstable one). A 273 K lake under 255 K air is
#    strongly convective, so the reference run sat in the VERY-unstable regime,
#    where CTSM caps the log at the transition point and adds a free-convection
#    correction that the port simply did not have.
#
# WHAT THESE TESTS ASSERT — and deliberately do NOT assert.
# The defect class this repo keeps hitting is a green test that is structurally
# blind (see MEMORY: conservation-is-not-accuracy, vacuous-checks-bug-class). A
# test that asserted `isfinite(t_ref2m)`, or `t_ref2m != 283.0`, would have
# passed on the broken code the moment ANY value was written, and would have
# proved nothing about the physics. So:
#   * the 2-m diagnostic is asserted against the FORTRAN REFERENCE VALUE, per
#     step, in Kelvin;
#   * each of the four stability regimes is asserted to (a) be REACHED by a
#     stated `zeta`, and (b) produce a value that DIFFERS materially from the
#     single regime the port used to apply everywhere — so a regression to the
#     old one-branch code fails these, rather than sliding through.
# ==========================================================================
using Test, CLM

@testset "Lake 2-m diagnostics + Monin-Obukhov stability regimes" begin

    # ------------------------------------------------------------------
    # C1. All four regimes are reachable, distinct, and agree with the
    #     independently-ported ARRAY form of CTSM's FrictionVelocity.
    #
    # `mo_profile_denom_m` / `mo_profile_denom_h` are the new scalar (GPU-safe)
    # form used by the lake kernel. `friction_velocity!` is the pre-existing
    # array/filter port of the SAME Fortran subroutine, written separately and
    # covered by test_friction_velocity.jl. If the scalar helper reproduces the
    # array port to round-off across all four branches, both are the same
    # Fortran — that is a real cross-check, not a restatement of the new code.
    # ------------------------------------------------------------------
    zetam = CLM.ZETAM_PROFILE   # 1.574, wind-profile transition
    zetat = CLM.ZETAT_PROFILE   # 0.465, temperature-profile transition

    # zldis is fixed at 30 m (the lake reference's forcing height); obu is chosen
    # to place zeta = zldis/obu squarely inside each regime.
    zldis = 30.0
    z0m   = 1.0e-3
    z0h   = 2.0e-4

    # (label, obu, expected regime index for the HEAT profile)
    heat_cases = [
        ("very unstable", -10.0,   1),   # zeta = -3.0   < -zetat
        ("unstable",      -100.0,  2),   # zeta = -0.3   in [-zetat, 0)
        ("stable",         60.0,   3),   # zeta = +0.5   in [0, 1]
        ("very stable",    10.0,   4),   # zeta = +3.0   > 1
    ]
    # Same obu values give the wind profile: zeta = -3.0 (< -zetam=1.574, regime 1),
    # -0.3 (regime 2), +0.5 (regime 3), +3.0 (regime 4) — all four also covered.
    wind_regime_of(zeta) = zeta < -zetam ? 1 : zeta < 0.0 ? 2 : zeta <= 1.0 ? 3 : 4
    heat_regime_of(zeta) = zeta < -zetat ? 1 : zeta < 0.0 ? 2 : zeta <= 1.0 ? 3 : 4

    for (label, obu, want_h) in heat_cases
        zeta = zldis / obu
        @testset "regime $want_h ($label), zeta = $(round(zeta, digits=3))" begin
            # (a) the branch is genuinely REACHED — not a test that only ever
            #     exercises the one regime the port already had.
            @test heat_regime_of(zeta) == want_h
            @test wind_regime_of(zeta) == want_h

            dh = CLM.mo_profile_denom_h(zeta, zldis, z0h, obu)
            dm = CLM.mo_profile_denom_m(zeta, zldis, z0m, obu)
            @test isfinite(dh) && dh > 0.0
            @test isfinite(dm) && dm > 0.0

            # (b) cross-check against the array port of the same Fortran.
            #     friction_velocity! computes temp1 (heat, at forcing height),
            #     temp12m (heat, at 2 m) and ustar (momentum) from obu/z0*.
            fv = CLM.FrictionVelocityData()
            CLM.frictionvel_init!(fv, 1, 1)
            fv.zetamaxstable = 0.5
            fv.forc_hgt_u_patch[1] = zldis
            fv.forc_hgt_t_patch[1] = zldis
            fv.forc_hgt_q_patch[1] = zldis
            displa  = [0.0]
            z0m_v   = [z0m]; z0h_v = [z0h]; z0q_v = [z0h]
            obu_v   = [obu]
            ur      = [3.0];  um = [3.0]
            ustar   = [0.0]
            temp1   = [0.0];  temp2 = [0.0]; temp12m = [0.0]; temp22m = [0.0]
            fmv     = [0.0]
            CLM.friction_velocity!(fv, 1, [1], displa, z0m_v, z0h_v, z0q_v,
                                   obu_v, 1, ur, um, ustar, temp1, temp2,
                                   temp12m, temp22m, fmv)

            # temp1 = VKC / denominator(heat, at zldis)
            @test CLM.VKC / dh ≈ temp1[1] rtol=1e-12
            # ustar = VKC * um / denominator(momentum, at zldis)
            @test CLM.VKC * um[1] / dm ≈ ustar[1] rtol=1e-12

            # the 2-m form is the same branch structure at zldis = 2 + z0h
            zldis2 = 2.0 + z0h
            dh2 = CLM.mo_profile_denom_h(zldis2 / obu, zldis2, z0h, obu)
            @test CLM.VKC / dh2 ≈ temp12m[1] rtol=1e-12

            # (c) the branch MATTERS: what the lake kernel used to apply
            #     everywhere is the regime-2 expression. Outside regime 2 it must
            #     give a materially different answer, or "porting the other three
            #     regimes" would be a no-op and this whole fix cosmetic.
            regime2_h = log(zldis / z0h) - CLM.stability_func2(min(zeta, 0.0)) +
                        CLM.stability_func2(min(z0h / obu, 0.0))
            if want_h == 2
                @test dh ≈ regime2_h rtol=1e-12       # regime 2 IS the old form
            else
                @test abs(dh - regime2_h) / abs(dh) > 0.05
            end
        end
    end

    # C2. The transition is continuous where CTSM makes it continuous: the
    # regime 1/2 boundary at zeta = -zetat is a matched join by construction
    # (regime 1 evaluates regime 2 at -zetat plus a correction that vanishes
    # there). A sign error in the correction term shows up here.
    @testset "regime 1/2 join at zeta = -zetat is continuous" begin
        obu = -zldis / zetat            # puts zeta exactly at -zetat
        eps = 1e-9
        below = CLM.mo_profile_denom_h(-zetat - eps, zldis, z0h, obu)
        above = CLM.mo_profile_denom_h(-zetat + eps, zldis, z0h, obu)
        @test isapprox(below, above; rtol=1e-6)
    end
    @testset "regime 3/4 join at zeta = 1 is continuous" begin
        obu = zldis                     # zeta = 1 exactly
        eps = 1e-9
        below = CLM.mo_profile_denom_h(1.0 - eps, zldis, z0h, obu)
        above = CLM.mo_profile_denom_h(1.0 + eps, zldis, z0h, obu)
        @test isapprox(below, above; rtol=1e-6)
        belowm = CLM.mo_profile_denom_m(1.0 - eps, zldis, z0m, obu)
        abovem = CLM.mo_profile_denom_m(1.0 + eps, zldis, z0m, obu)
        @test isapprox(belowm, abovem; rtol=1e-6)
    end

    # ------------------------------------------------------------------
    # B. The lake 2-m air temperature against the FORTRAN REFERENCE, per step.
    #
    # GATED: the reference h0 and the forcing/param files are machine-local.
    # Runs in its own module (the harness mutates module globals).
    # ------------------------------------------------------------------
    @testset "lake TSA matches Fortran per step (gated)" begin
        script = joinpath(@__DIR__, "..", "scripts", "fortran_parity_lake.jl")
        if !isfile(script)
            @info "lake parity harness absent, skipping"
            @test_skip isfile(script)
        else
            mod = Module(:LAKE_REF2M)
            Core.eval(mod, :(using Test, CLM, NCDatasets, Dates, Printf))
            Base.include(mod, script)   # PROGRAM_FILE-guarded: no auto-run
            main_fn = Base.invokelatest(getfield, mod, :main)
            rep = Base.invokelatest(main_fn; nsteps=24)
            if rep === missing
                @info "lake reference inputs absent, skipping"
                @test_skip rep === true
            else
                jt = rep.jv["TSA"]; ft = rep.fv["TSA"]
                @test length(jt) == 24

                # The assertion that matters: the VALUE, in Kelvin, against the
                # Fortran reference at every step. Before the fix this was a flat
                # 283.000 against a 255-270 K reference — errors of 13-28 K.
                # The bound is the MEASURED post-fix residual (max ~1.9 K, on the
                # cold-start transient at step 1-2; ~0.1 K by the end of the run),
                # not a target that was tuned until it passed.
                @test all(isfinite, jt)
                @test maximum(abs.(jt .- ft)) < 2.5

                # ...and it is a genuinely evolving field, not a new constant:
                # the reference TSA moves ~10 K over the first day.
                @test (maximum(jt) - minimum(jt)) > 3.0
                @test all(!=(283.0), jt)   # the exact cold-start value never recurs

                # The reference exercises regime 1 (very unstable) at EVERY step —
                # zeta stays in [-100, -15.7], never above the -zetat = -0.465
                # transition. So the branch this fix adds is the branch this
                # reference actually runs in, and a regression to the old
                # regime-2-only kernel changes every one of these steps.
                @test all(z -> z < -0.465, rep.zetas)
                @test all(==(1), rep.regimes)

                # Surface fluxes: this is the CHARACTERISED residual, not a pass
                # target. See docs/LAKE_FLUX_RESIDUAL.md finding D — the port's
                # lake surface never drops below freezing (TG 273-276 K vs the
                # reference's 267-269 K), so it stays on the UNFROZEN roughness
                # branch and gets z0hg ~25x smaller than the frozen z0frzlake
                # branch Fortran is on, leaving rah ~300 s/m against Fortran's
                # ~150. Tightening these bounds requires fixing that, not this.
                jh = rep.jv["FSH"]; fh = rep.fv["FSH"]
                @test maximum(abs.(jh .- fh) ./ (1 .+ abs.(fh))) < 0.30
                jl = rep.jv["EFLX_LH"]; fl = rep.fv["EFLX_LH"]
                @test maximum(abs.(jl .- fl) ./ (1 .+ abs.(fl))) < 0.30
            end
        end
    end
end
