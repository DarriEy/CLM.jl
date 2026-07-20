# ==========================================================================
# test_lake_eddy_ks.jl — docs/LAKE_FLUX_RESIDUAL.md finding D.
#
# THE DEFECT. `LakeFluxesMod.F90:755` sets the depth-decay rate of the lake's
# wind-driven eddy diffusivity:
#
#     ks = 6.6 * sqrt(|sin(lat)|) * u2m^-1.84        [1/m], lat in RADIANS
#
# and `LakeTemperatureMod.F90:373` uses it as
#
#     ke = vkc*ws*z/p0 * exp(-ks*z) / (1 + 37*ri^2)
#
# The port hardwired `ks_col = 0` because the gridcell latitude had never been
# wired into the lake kernel. `exp(-0*z) == 1` at EVERY depth, so the eddy
# diffusivity never decayed: `ke` grew linearly all the way to the lake bed and
# the whole 50 m column stayed coupled to the skin. Every W/m2 the surface lost
# was replenished from deep 277 K water instead of from the top ~0.3 m, so the
# surface could not reach freezing — and because it never froze it stayed on the
# UNFROZEN roughness branch (`z0mg ~3e-5`) instead of `z0frzlake = 1e-3`, leaving
# `rah` ~305 s/m against the reference's ~145.
#
# WHAT THESE TESTS ASSERT — and deliberately do NOT assert.
# `ks != 0` would pass on any nonzero write and prove nothing (see MEMORY:
# vacuous-checks-bug-class). So every assertion below is a VALUE:
#   * `ks` in 1/m against the Fortran formula evaluated at the reference's own
#     latitude and 2-m wind — a number, not a sign test;
#   * the physical CONSEQUENCE, exp(-ks*z) at 5 m depth, which the old code got
#     exactly wrong by 4+ orders of magnitude (it was exactly 1.0);
#   * `t_grnd`, `z0mg` and `rah` against the FORTRAN REFERENCE, per step, in
#     their own units — the three numbers that decide whether D is closed.
# ==========================================================================
using Test, CLM

@testset "Lake eddy-diffusivity depth decay (ks) and the frozen roughness branch" begin

    # ------------------------------------------------------------------
    # D1. The Fortran formula, evaluated independently of the port, at the
    #     reference's conditions. This fixes the EXPECTED MAGNITUDE that the
    #     driver test below is checked against, so a regression to ks = 0 (or
    #     to degrees-instead-of-radians) fails on a number rather than sliding
    #     through a sign test.
    # ------------------------------------------------------------------
    ks_fortran(lat_rad, u2m) = 6.6 * sqrt(abs(sin(lat_rad))) * u2m^(-1.84)

    lat_bow = 51.17 * pi / 180        # the lake reference gridcell
    u2m_ref = 1.6                     # ustar/vkc*log(2/z0mg) at the reference

    ks_ref = ks_fortran(lat_bow, u2m_ref)
    @test 2.0 < ks_ref < 3.0          # ~2.4 /m

    # The consequence, which is what the lake column actually feels. At 5 m the
    # broken code applied the FULL surface eddy diffusivity; the correct code
    # suppresses it by more than four orders of magnitude.
    @test exp(-ks_ref * 5.0) < 1.0e-4
    @test exp(-0.0     * 5.0) == 1.0     # what the port did before this fix

    # Degrees-for-radians is a silent unit error in the argument of sin. Pin the
    # unit by ANALYTIC VALUE at a latitude whose sine is exact:
    #   lat = pi/6 rad = 30 deg -> |sin| = 0.5 exactly.
    # Reading pi/6 as DEGREES gives sin(0.5236 deg) = 0.00914, i.e. sqrt() is 7.4x
    # smaller — so this single equality separates the two readings unambiguously.
    @test ks_fortran(pi/6, 1.0) ≈ 6.6 * sqrt(0.5) rtol=1e-12
    @test ks_fortran(pi/6, 1.0) / ks_fortran(deg2rad(pi/6), 1.0) > 7.0
    #
    # NOTE ON THE PREMISE, because the first version of this pin FAILED and the
    # failure was correct: it compared ks(51.17) against ks(51.17*pi/180) at the
    # reference's own latitude and asserted they differ by >10%. They differ by
    # 0.46% — 51.17 rad wraps to 0.905 rad, whose sine (0.786) is coincidentally
    # almost exactly sin(51.17 deg) = 0.779. Testing a unit convention at a value
    # where the two conventions happen to agree proves nothing; the fix is to move
    # the test to where the distinction actually lives, not to lower the bound.
    # Bow is still checked below, but for MAGNITUDE, not for units.

    # The equator is the one place the sqrt has an infinite derivative, which is
    # why the port had written a literal zero. The branch must still return an
    # exact, finite zero there (and never a NaN) so AD stays clean.
    @test ks_fortran(0.0, u2m_ref) == 0.0

    # ------------------------------------------------------------------
    # D2. The port, against the FORTRAN REFERENCE, per step.
    #
    # GATED: the reference h0 and the forcing/param files are machine-local.
    # Runs in its own module (the harness mutates module globals).
    # ------------------------------------------------------------------
    @testset "lake freezes and picks up z0frzlake, vs Fortran (gated)" begin
        script = joinpath(@__DIR__, "..", "scripts", "fortran_parity_lake.jl")
        if !isfile(script)
            @info "lake parity harness absent, skipping"
            @test_skip isfile(script)
        else
            mod = Module(:LAKE_EDDY_KS)
            Core.eval(mod, :(using Test, CLM, NCDatasets, Dates, Printf))
            Base.include(mod, script)   # PROGRAM_FILE-guarded: no auto-run
            main_fn = Base.invokelatest(getfield, mod, :main)
            rep = Base.invokelatest(main_fn; nsteps=16)
            if rep === missing
                @info "lake reference inputs absent, skipping"
                @test_skip rep === true
            else
                # --- the fix itself, as a value in 1/m -------------------
                # ks was IDENTICALLY ZERO at every step before this fix.
                @test length(rep.kss) == 16
                @test all(k -> 1.5 < k < 4.0, rep.kss)

                # --- the surface reaches freezing ------------------------
                # Before: min(TG_J) = 273.246 K, i.e. it NEVER went below tfrz,
                # against a reference that is at 267-269 K from its first record.
                tgJ = rep.jv["TG"]; tgF = rep.fv["TG"]
                @test minimum(tgJ) < CLM.TFRZ

                # --- and once frozen it tracks the reference in KELVIN ---
                # Steps 7-16 are the frozen window (the port cold-starts from a
                # uniform 277 K lake while the reference's first record already
                # shows a partly-cooled 274.81 K top layer, so the port takes ~6
                # steps to catch up; that lag is the remaining residual and is
                # NOT tuned away here).
                @test maximum(abs.(tgJ[7:16] .- tgF[7:16])) < 0.15

                # --- the frozen ROUGHNESS branch is the one being taken --
                # z0frzlake = 1e-3 (LakeCon.F90:50) vs the unfrozen Charnock
                # z0mg ~3.2e-5 the port was stuck on. Assert the value, not
                # "changed": these differ by 31x.
                @test all(z -> isapprox(z, CLM.z0frzlake; rtol=1e-12), rep.z0mgs[7:16])
                @test all(z -> z < 1e-4, rep.z0mgs[1:5])   # unfrozen early steps

                # --- which puts rah onto the reference -------------------
                # rah is inverted from each model's OWN flux definition
                # (LakeFluxesMod.F90:526), so this compares like with like.
                # Before the fix: rah_J 296-317 against rah_F 145-160, i.e. 2x.
                @test all(r -> 100.0 < r < 200.0, rep.rahs[7:16])

                # --- and lake ice forms ----------------------------------
                # Before: LAKEICEFRAC_SURF was 0.0 at every one of 47 steps
                # against a reference reaching 0.36.
                iceJ = rep.jv["LAKEICE_SRF"]; iceF = rep.fv["LAKEICE_SRF"]
                @test maximum(iceJ) > 0.2
                @test maximum(abs.(iceJ .- iceF)) < 0.10
            end
        end
    end
end
