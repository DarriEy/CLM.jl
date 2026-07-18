# =============================================================================
# test_smooth_axis_guard.jl
#
# GUARD AGAINST THE "k IS DIMENSIONAL" BUG CLASS.
#
# The AD-smoothing primitives in src/infrastructure/smooth_ad.jl are LogSumExp /
# sigmoid surrogates with a sharpness `k`. The critical, easily-forgotten fact:
#
#     smooth_max(0, x)  OVERSHOOTS  max(0, x)  by up to  log(2)/k
#     ...AND THAT OVERSHOOT IS IN THE PHYSICAL UNITS OF THE AXIS BEING SMOOTHED.
#
# The generic default k = 50 is calibrated for an O(1) axis, where log(2)/50 = 0.0139
# rounds a corner. On an axis that is NOT O(1) it does not round a corner — it
# fabricates physics:
#
#   * metres of water        -> 13.9 mm of water conjured from nothing
#   * a stomatal conductance -> a 0.0139 mol/m2/s = 13 900 umol/m2/s FLOOR on gs
#   * a layer transmittance  -> thick snow becomes semi-transparent
#   * an albedo              -> a flat +0.0139 additive albedo offset
#   * a clamp box NARROWER than 0.0139 -> the clamp is DESTROYED, not softened
#
# And a `smooth_max(0, ·)` / `smooth_min(·, cap)` ReLU GUARD sits AT its kink in the
# COMMON case, so the error fires every step on every column — it is not a rare-branch
# rounding.
#
# This file is the regression guard. It enforces three things:
#   (1) the error law itself, so the hazard is stated in executable form;
#   (2) every NAMED axis-scaled sharpness constant satisfies the nondimensional
#       criterion  log(2)/k <= REL_TOL * (characteristic scale of its axis);
#   (3) the specific guards that were found biasing the physics stay HARD, and no NEW
#       default-k smooth_* call appears without a human looking at its axis (a census pin).
# =============================================================================
using Test
using CLM

const _SRCDIR = normpath(joinpath(@__DIR__, "..", "src"))

@testset "smooth axis guard" begin

    # -------------------------------------------------------------------------
    # (1) The error law, in executable form.
    # -------------------------------------------------------------------------
    @testset "error law: smooth_max(0,x) overshoots by log(2)/k at the kink" begin
        old = CLM.SMOOTH_MODE[]
        CLM.SMOOTH_MODE[] = :always
        try
            for k in (50.0, 500.0, 1.0e6)
                # AT the kink (x == 0) the overshoot is exactly log(2)/k, in axis units.
                @test CLM.smooth_max(0.0, 0.0; k = k) ≈ log(2) / k rtol = 1e-10
                # The overshoot is scale-covariant: it does NOT shrink just because the
                # axis is small. This is the whole bug.
                @test CLM.smooth_max(0.0, 1.0e-9; k = k) > 1.0e-9
            end
            # The headline number: the generic default puts the width at 1.39e-2 AXIS UNITS.
            @test CLM.smooth_max(0.0, 0.0) ≈ 0.013862943611198907 rtol = 1e-12

            # A clamp box narrower than the smoothing width is DESTROYED, not softened.
            # (This is what happened to qcharge: box +/-10/1800 = +/-5.56e-3 mm/s, k=50.)
            bound = 10.0 / 1800.0
            @test log(2) / 50.0 > bound                       # width exceeds the HALF-box
            destroyed = CLM.smooth_min(bound, CLM.smooth_max(-bound, 0.0))
            @test abs(destroyed - 0.0) > 0.5 * bound          # a true clamp would return ~0
            # ...and scaling k to the axis restores it.
            ok = CLM.smooth_min(bound, CLM.smooth_max(-bound, 0.0; k = 1e9); k = 1e9)
            @test isapprox(ok, 0.0; atol = 1e-8)
        finally
            CLM.SMOOTH_MODE[] = old
        end
    end

    # -------------------------------------------------------------------------
    # (2) Registry: every named axis-scaled sharpness must be small on ITS OWN axis.
    #
    # Each entry declares the axis, its units, and the characteristic scale that the
    # smoothing width must be small against. The criterion is NONDIMENSIONAL:
    #
    #       log(2)/k  <=  REL_TOL * characteristic_scale
    #
    # This is what "k must be scaled to the axis" means, and it is why a bare k = 50 is
    # only ever defensible on an axis whose characteristic scale is O(1).
    # -------------------------------------------------------------------------
    @testset "named sharpness constants are scaled to their axis" begin
        REL_TOL = 1.0e-3   # smoothing width must be <= 0.1% of the axis's characteristic scale

        # (name, k, axis description, units, characteristic scale in those units)
        registry = [
            ("SOIL_HYDRAULIC_K", CLM.SOIL_HYDRAULIC_K,
             "relative saturation / volumetric water / layer liquid mass / aquifer recharge",
             "[-], m3/m3, mm, mm/s", 1.0e-2),
            ("SNOW_PERC_WATER_K", CLM.SNOW_PERC_WATER_K[],
             "snow-layer drainage water depth", "m", 1.0e-4),
            ("PHASE_CHANGE_MASS_K", CLM.PHASE_CHANGE_MASS_K[],
             "phase-change water mass", "kg/m2", 1.0e-2),
        ]

        for (name, k, axis, units, scale) in registry
            width = log(2) / k
            @test width <= REL_TOL * scale
            # Guard the guard: a scale of 0 or a k of 50 would silently pass nothing.
            @test k > 50.0
            @info "axis-scaled sharpness OK" name k axis units scale width
        end

        # BTRAN_SMOOTH_K is DELIBERATELY left at the generic k=50: after this PR its only
        # remaining uses are on a relative-saturation [0,1] axis and a matric-potential axis
        # whose magnitude is ~1e5 mm, i.e. axes where 0.0139 is genuinely negligible or is
        # the intended physical wilting-point transition width. The rresis cap at the CONSTANT
        # 1.0 — which cost a flat 1.4% of btran on every well-watered layer — is now a hard min.
        @test CLM.BTRAN_SMOOTH_K[] == 50.0
    end

    # -------------------------------------------------------------------------
    # (3a) Regression pins: the guards that were fabricating physics must stay HARD.
    #
    # Each of these clamps a quantity against a CONSTANT. The clamped branch therefore
    # carries no derivative information at all, so a smooth surrogate buys NOTHING for AD
    # and only injects a log(2)/k bias into the primal. Re-smoothing any of them would
    # silently reintroduce the bug, with no balance check to catch it.
    # -------------------------------------------------------------------------
    @testset "hardened guards stay hard" begin
        # file => list of (banned substring, why)
        banned = Dict(
            "biogeophys/photosynthesis.jl" => [
                ("smooth_max(smooth_max(r1, r2)", "gs_mol: quadratic_solve already SORTS its roots (r1>=r2), so this is a no-op in exact physics — but on the Medlyn conductance axis (mol/m2/s) it imposed a 0.0139 mol/m2/s = 13900 umol/m2/s floor on gs"),
                ("smooth_max(zero(r1), smooth_min(r1, r2))", "ag: fabricated ~0.008-0.024 umol/m2/s of photosynthesis on patches whose true assimilation is exactly 0"),
                ("smooth_max(soilflux", "PHS transpiration: a ReLU on a kg/m2/s axis, sitting at its kink -> 0.0139 mm/s = 1200 mm/day"),
            ],
            "biogeophys/canopy_fluxes.jl" => [
                ("smooth_clamp(zeta_patch", "Monin-Obukhov zeta: the +/-0.01 box is NARROWER than the k=50 width 0.0139, so smoothing DESTROYED the limiter"),
            ],
            "biogeophys/bareground_fluxes.jl" => [
                ("smooth_clamp(zeta_patch", "as canopy_fluxes: the zeta box is narrower than the smoothing width"),
            ],
            "biogeophys/lake_fluxes.jl" => [
                ("smooth_max(z0mg, T(1.0e-10))", "lake roughness: the 1e-10 m floor became 0.0139 m, giving a lake the roughness of a forest canopy"),
            ],
            "biogeophys/lake_temperature.jl" => [
                ("smooth_max(n2, n2min_p)", "Brunt-Vaisala floor on an s^-2 axis"),
                ("smooth_max(ws[c]^T(2.0)", "squared friction velocity (m2/s2, ~1e-6): the 1e-10 floor became 0.0139, collapsing the Richardson number and inflating lake eddy diffusivity ~56x"),
            ],
            "biogeophys/snow_snicar.jl" => [
                ("smooth_max(exp_min", "SNICAR layer transmittance: an optically thick layer transmits 1e-9..1e-2, far below the 0.0139 width -> thick snow became semi-transparent"),
                ("smooth_max(0.0, flx_abs_lcl", "absorbed-flux ReLU: +0.0139 of incident SW spuriously absorbed in every snow layer"),
            ],
            "biogeophys/surface_albedo.jl" => [
                ("smooth_max(T(0.11) - T(0.40)", "soil albedo increment: a ReLU that sits at its kink for ALL WET SOIL, adding a flat +0.0139 of ALBEDO"),
                ("smooth_max(cn.alblakwi[ib], T(0.10))", "lake albedo: alblakwi == 0.10 in the params, so this ReLU sat exactly on its kink and added exactly log(2)/50 of albedo"),
            ],
            "biogeophys/soil_temperature.jl" => [
                ("smooth_max(frac_sno[c], T(1.0e-6))", "snow bulk density denominator: raised a 1e-6 floor to 0.0139 in TWO unit systems at once (a [0,1] fraction and metres)"),
                ("smooth_max(zero(T), log10(satw)", "Kersten number: constant-0 branch, and log10 has a singular derivative in the clamped tail"),
            ],
            "biogeophys/soil_moist_stress.jl" => [
                ("rresis_j = smooth_min", "rresis cap at the CONSTANT 1.0 cost a flat 1.4% of btran on every well-watered rooted layer -> -1.4% GPP/transpiration"),
            ],
            "biogeochem/methane.jl" => [
                # The ch4_aere! KERNEL (the live runtime path) had smoothed the aerenchyma/
                # transpiration CH4 & O2 fluxes, which live on a mol/m3/s axis ~1e-9. The
                # log(2)/50 = 0.0139 floor is ~7 ORDERS OF MAGNITUDE larger than the flux, so
                # the smoothed (AD / SMOOTH_MODE=:always) methane transport was fabricated. The
                # site-level reference site_ox_aere! already used hard max on the SAME axis
                # (with a "0.0139 is 7 orders too big" comment); the kernel now matches it.
                ("aere = smooth_max(aere", "CH4 aerenchyma flux (mol/m3/s ~1e-9): 0-floor became 0.0139, 7 orders > flux; site_ox_aere! l.1271 is hard"),
                ("oxaere = smooth_max(oxaere", "O2 aerenchyma flux (mol/m3/s ~1e-9): same 7-orders floor; site_ox_aere! l.1280 is hard"),
                ("tranloss = smooth_max(tranloss", "CH4 transpiration flux (mol/m3/s ~1e-9): same 7-orders floor; site_ox_aere! l.1238 is hard"),
                ("smooth_min(aere + tranloss", "scattered CH4 aerenchyma flux capped by availability — smoothing biases the TRANSPORTED CH4 by 0.0139 mol/m3/s"),
                ("smooth_min(tranloss, aeretran)", "scattered CH4 transpiration flux — same mol/m3/s flux-axis bias"),
            ],
        )

        # Match CODE only — these files now document the bug class in prose, and the comments
        # deliberately quote the offending expressions.
        strip_comments(s) = join((replace(ln, r"#.*$" => "") for ln in split(s, '\n')), "\n")

        for (relpath, entries) in banned
            path = joinpath(_SRCDIR, relpath)
            @test isfile(path)
            src = strip_comments(read(path, String))
            for (needle, why) in entries
                if occursin(needle, src)
                    @warn "REINTRODUCED a dimensional-k smoothing bug" file=relpath pattern=needle reason=why
                end
                @test !occursin(needle, src)
            end
        end
    end

    # -------------------------------------------------------------------------
    # (3b) Census pin.
    #
    # Any smooth_* call that does NOT pass an explicit `k` inherits the generic k = 50,
    # which is ONLY defensible on an O(1) axis. We cannot infer an axis's units statically,
    # so instead we PIN THE COUNT: a new default-k smoothed guard makes this test fail, and
    # whoever adds it must either (a) confirm the axis really is O(1) and bump the baseline,
    # or (b) pass an axis-scaled `k` / use a hard guard. That forces the units question to be
    # asked exactly once per call site, which is the whole lesson of this bug class.
    # -------------------------------------------------------------------------
    @testset "no NEW default-k smooth_* call sites" begin
        baseline = Dict(
            # Fully hardened by this PR — every smooth_* guard in these files clamped against a
            # CONSTANT, so they are pinned at ZERO: any reintroduction is a regression.
            "biogeophys/lake_fluxes.jl" => 0,
            "biogeophys/soil_moist_stress.jl" => 0,
            # Remaining default-k sites. Each still needs its axis checked before it is trusted;
            # the ones left are either on genuinely O(1) axes, or are latent (biogeochem is off the
            # default SP path). See the PR that added this file for the full ranked audit.
            "biogeochem/decomp_competition.jl" => 1,
            "biogeochem/decomp_mimics.jl" => 20,
            "biogeochem/fire_li2014.jl" => 28,
            "biogeochem/fire_li2016.jl" => 24,
            "biogeochem/fire_li2021.jl" => 21,
            "biogeochem/fire_li2024.jl" => 21,
            "biogeochem/methane.jl" => 18,   # 5 aerenchyma/transpiration mol/m3/s flux guards hardened (see "hardened guards stay hard")
            "biogeochem/phenology.jl" => 34,
            "biogeochem/veg_struct_update.jl" => 23,
            "biogeophys/bareground_fluxes.jl" => 4,
            "biogeophys/canopy_fluxes.jl" => 10,
            "biogeophys/lake_temperature.jl" => 18,
            "biogeophys/luna.jl" => 39,
            "biogeophys/photosynthesis.jl" => 63,
            "biogeophys/qsat.jl" => 3,
            "biogeophys/snow_cover_fraction.jl" => 1,
            "biogeophys/snow_hydrology.jl" => 26,
            "biogeophys/snow_snicar.jl" => 18,
            "biogeophys/snow_snicar_device.jl" => 16,
            "biogeophys/soil_hydrology.jl" => 1,
            "biogeophys/soil_temperature.jl" => 2,
            "biogeophys/soil_water_movement.jl" => 24,
            "biogeophys/surface_albedo.jl" => 23,
            "biogeophys/urban_fluxes.jl" => 4,
            "driver/clm_driver.jl" => 1,
        )

        call_re = r"\bsmooth_(?:max|min|clamp|heaviside|abs|ifelse)\s*\("
        # a call is "axis-scaled" if it passes an explicit k on the same line
        kwarg_re = r"\bsmooth_(?:max|min|clamp|heaviside|abs|ifelse)\s*\([^\n]*;\s*k\s*="

        # Count CODE only. These files document the bug class in prose, and a `smooth_max(...)`
        # written inside a comment must not be counted as a call site.
        strip_comments(s) = join((replace(ln, r"#.*$" => "") for ln in split(s, '\n')), "\n")

        for (relpath, expected) in sort(collect(baseline))
            path = joinpath(_SRCDIR, relpath)
            @test isfile(path)
            src = strip_comments(read(path, String))
            n_default = length(collect(eachmatch(call_re, src))) -
                        length(collect(eachmatch(kwarg_re, src)))
            if n_default > expected
                @warn """
                      NEW default-k smooth_* call site(s) added.
                      `k` IS DIMENSIONAL: smooth_max(0,x) overshoots by log(2)/k = 0.0139 IN THE
                      UNITS OF THE AXIS. Before bumping this baseline, state the axis and its
                      characteristic scale and confirm 0.0139 is negligible on it. If it is not,
                      pass an axis-scaled `k` (see SOIL_HYDRAULIC_K) or use a hard guard — a clamp
                      against a CONSTANT has zero derivative on the clamped branch anyway, so
                      smoothing it buys nothing.
                      """ file=relpath expected n_default
            end
            @test n_default <= expected
        end
    end
end
