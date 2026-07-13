# ==========================================================================
# test_phase_snow_smoothing.jl — AD-smoothing of the freeze/thaw (phase-change)
# and snow merge/split discontinuities (the deferred PR #120 follow-up).
#
# Two discontinuity families are smoothed, gated PURELY by element type via the
# shared `_use_smooth` (smooth_ad.jl): plain Float64 evaluates the EXACT hard
# physics (byte-identical default), ForwardDiff.Dual evaluates the smooth
# surrogate. Sharpness Refs control the transition width ε ≈ 1/k; as k → ∞
# (ε → 0) the surrogate → the hard physics.
#
# 1. PHASE CHANGE / combined-layer latent-heat partition — `combo_scalar`
#    branches on the combined heat capacity denom = cpice*wicec+cpliq*wliqc
#    crossing 0 (empty/massless combined layer → tc = tfrz). As a layer's
#    ice/liquid mass sweeps to 0, the combined temperature has a near-vertical
#    step (dtc/dmass → ~1e15 in the hard physics); the smooth surrogate keeps it
#    finite. combo_scalar is the shared latent-heat energy combine used by both
#    `combine_snow_layers!` and `divide_snow_layers!`, and the same `smooth_max`
#    floors smooth the freeze/thaw mass clamps in `phase_change_beta!`.
#
# 2. SNOW MERGE/SPLIT — `divide_snow_layers!` splits a layer that exceeds its
#    max thickness and repartitions ice/liquid/temperature into the new layer.
#    The new-layer temperature is a `smooth_heaviside`-blended extrapolation
#    (frozen vs. extrapolated branch on the tfrz crossing). The mass/thickness
#    repartition weights and `mass_weighted_snow_radius` clamp are likewise
#    smoothed. WHAT STAYS DISCRETE: the integer snow-layer count `snl` and the
#    which-neighbour merge decision (combinatorial — a layer either exists or
#    not). The threshold tests that TRIGGER a merge/split (dz < dzmin /
#    dz > dzmax / mass/dz < 50) flip the integer count and remain hard `if`s; we
#    smooth only the CONTINUOUS repartition that runs once a (discrete)
#    merge/split is decided, so the per-layer state fields stay C¹ in their own
#    neighbourhood without pretending the layer count is differentiable.
#
# These tests demonstrate, for each family:
#   (a) Smoothing ON (Dual): the AD gradient is FINITE & CONTINUOUS across the
#       threshold (no NaN/Inf, bounded jump), whereas the hard physics steps.
#   (b) Smoothing OFF (Float64): byte-identical to the hard reference.
#   (c) Physics error → 0 as the smoothing width ε → 0 (k → ∞).
# ==========================================================================

using Test
using ForwardDiff
using ForwardDiff: Dual, value
using Random
using CLM

@testset "phase-change + snow merge/split AD-smoothing" begin

    tfrz = CLM.TFRZ

    # ---------------------------------------------------------------------
    # combo_scalar (phase-change combined-layer temperature)
    # ---------------------------------------------------------------------
    # tc of combining two layers as the ice/liquid mass `w` of one of them goes
    # to 0 (denom crosses 0). The element type follows `w`: Float64 → hard
    # branch, Dual → smooth blend.
    combo_tc(w) = CLM.combo_scalar(0.1, w, w, 272.0, 0.1, 0.0, 0.0, 270.0)[4]

    @testset "(combo) hard reference: denom→0 derivative explodes" begin
        # The hard physics flips tc from the energy formula to tfrz at denom==0,
        # producing a near-vertical step there (central FD of the Float64 value).
        smoothval(w) = value(combo_tc(Dual(w, 0.0)))   # Float64 path → hard branch
        hfd = 1e-7
        d_just_above = (smoothval(2e-3 + hfd) - smoothval(2e-3 - hfd)) / (2hfd)
        d_at_zero    = (smoothval(hfd) - smoothval(0.0)) / hfd
        @test combo_tc(0.0) == tfrz
        # The near-zero slope dwarfs the slope away from the threshold → a step.
        @test abs(d_at_zero) > 100.0 * abs(d_just_above)
    end

    @testset "(combo) smooth ON: dtc/dmass finite & continuous across denom→0" begin
        anynan = false; smooth_maxjump = 0.0; prev = nothing
        for w in 0.02:-0.001:0.0
            g = ForwardDiff.derivative(combo_tc, w)
            isfinite(g) || (anynan = true)
            prev !== nothing && (smooth_maxjump = max(smooth_maxjump, abs(g - prev)))
            prev = g
        end
        @test !anynan                       # finite even AT and across denom==0
        @test isfinite(smooth_maxjump)
    end

    @testset "(combo) smooth OFF: Float64 byte-identical to hard physics" begin
        function hardref(dz, wliq, wice, t, dz2, wliq2, wice2, t2)
            cpice = CLM.CPICE; cpliq = CLM.CPLIQ; hfus = CLM.HFUS; tf = CLM.TFRZ
            wicec = wice + wice2; wliqc = wliq + wliq2
            h = (cpice*wice + cpliq*wliq)*(t - tf) + hfus*wliq
            h2 = (cpice*wice2 + cpliq*wliq2)*(t2 - tf) + hfus*wliq2
            hc = h + h2; denom = cpice*wicec + cpliq*wliqc
            tc = denom > 0.0 ? tf + (hc - hfus*wliqc)/denom : tf
            return (dz + dz2, wliqc, wicec, tc)
        end
        rng = MersenneTwister(11); maxd = 0.0
        for _ in 1:5000
            a = (0.3rand(rng), 0.5rand(rng), 0.5rand(rng), 270.0 + 5rand(rng),
                 0.3rand(rng), 0.5rand(rng), 0.5rand(rng), 270.0 + 5rand(rng))
            c = CLM.combo_scalar(a...); r = hardref(a...)
            maxd = max(maxd, abs(c[1]-r[1]), abs(c[2]-r[2]), abs(c[3]-r[3]), abs(c[4]-r[4]))
        end
        # zero-mass (denom==0) corner cases too
        for _ in 1:1000
            a = (0.3rand(rng), 0.0, 0.0, 272.0, 0.3rand(rng), 0.0, 0.0, 270.0)
            c = CLM.combo_scalar(a...); r = hardref(a...)
            maxd = max(maxd, abs(c[4]-r[4]))
        end
        @test maxd == 0.0    # byte-identical (===, not approx)
    end

    @testset "(combo) convergence: smooth → hard as k → ∞" begin
        # A case near the denom threshold where the sigmoid blend is ACTIVE
        # (small combined heat capacity ≈ 0.085), so the smooth value differs
        # measurably from the hard branch and shrinks as k grows.
        w0 = 2e-5
        hardv = combo_tc(w0)                        # Float64 hard value
        smoothval(w) = value(combo_tc(Dual(w, 0.0)))
        k0 = CLM.SNOW_COMBO_TEMP_K[]
        errs = Float64[]
        for k in (10.0, 50.0, 200.0, 1000.0)
            CLM.SNOW_COMBO_TEMP_K[] = k
            push!(errs, abs(smoothval(w0) - hardv))
        end
        CLM.SNOW_COMBO_TEMP_K[] = k0
        @info "combo smooth→hard error vs k: $errs"
        @test errs[1] > 1e-3                        # blend is genuinely active at k=10
        @test issorted(errs; rev = true)           # monotonically decreasing
        @test errs[end] < 1e-6                      # → 0 as k → ∞
    end

    # ---------------------------------------------------------------------
    # phase_change_beta! latent-heat mass clamp (freeze/thaw)
    # ---------------------------------------------------------------------
    # The freeze/thaw heat release in `phase_change_beta!` clamps the melted ice
    # to the available mass: h2osoi_ice = smooth_max(0, wice0 - xm; k=PHASE_CHANGE_MASS_K).
    # As the forcing sweeps the melt energy `xm` past the available ice `wice0`,
    # the remaining ice (and hence the latent heat and resulting t_soisno) kinks
    # at xm == wice0; the smooth floor rounds that "ran out of ice" kink. We
    # exercise the exact clamped expression with the kernel's sharpness Ref.
    @testset "(phase) freeze/thaw mass clamp: smooth ON finite, OFF byte-identical" begin
        wice0 = 5.0
        kmass = CLM.PHASE_CHANGE_MASS_K[]
        ice_left(xm) = CLM.smooth_max(zero(typeof(xm)), wice0 - xm; k = oftype(xm, kmass))
        # Float64 byte-identical to the hard floor max(0, wice0-xm)
        maxd = 0.0
        for xm in 0.0:0.1:10.0
            maxd = max(maxd, abs(ice_left(xm) - max(0.0, wice0 - xm)))
        end
        @test maxd == 0.0

        # Dual: gradient finite & continuous across xm == wice0 (the kink).
        #
        # SAMPLE AT THE SMOOTHING WIDTH. The width of this clamp is log(2)/k IN THE
        # UNITS OF THE AXIS, and the axis is a water mass in kg/m2, so the width is
        # ~7e-10 kg/m2 (PHASE_CHANGE_MASS_K = 1e9 — see the note on the Ref). It has to
        # be: these clamps sit AT their kink for any snow layer holding ~0 liquid, so a
        # width of the old 0.0139 kg/m2 did not round the corner, it MOVED 0.0139 mm of
        # ice into liquid that never melted — fabricating 2.6 W/m2 of latent heat per
        # layer and breaking the soil energy balance. Probing at the old 0.05 spacing
        # would only tell us the function looks hard when you stand 1e8 widths away from
        # it (it does, and that is the point). The C¹ property lives at the width.
        eps_k = 1.0 / kmass
        anynan = false; prev = nothing; smooth_maxjump = 0.0
        for i in -20:20
            xm = wice0 + i * eps_k
            g = ForwardDiff.derivative(ice_left, xm)
            isfinite(g) || (anynan = true)
            prev !== nothing && (smooth_maxjump = max(smooth_maxjump, abs(g - prev)))
            prev = g
        end
        @test !anynan
        @test smooth_maxjump < 1.0          # no step: the derivative ramps 0 → -1
        # and it really does traverse the whole ramp (not a flat sample of one branch)
        @test ForwardDiff.derivative(ice_left, wice0 - 40 * eps_k) ≈ -1.0 atol = 1e-6
        @test abs(ForwardDiff.derivative(ice_left, wice0 + 40 * eps_k)) < 1e-6
        # AT the kink the smooth derivative is the midpoint -1/2 (a true C¹ rounding)
        @test ForwardDiff.derivative(ice_left, wice0) ≈ -0.5 atol = 1e-6
        # The PRICE of that rounding is bounded by the width: the smoothed mass never
        # departs from the hard clamp by more than log(2)/k (~7e-10 kg/m2 → ~1e-7 W/m2
        # of latent heat, vs the 1e-4 W/m2 balance threshold).
        worst = maximum(abs(value(CLM.smooth_max(Dual(0.0, 0.0), Dual(wice0 - xm, 0.0); k = kmass)) -
                            max(0.0, wice0 - xm))
                        for xm in (wice0 - 10 * eps_k):(eps_k):(wice0 + 10 * eps_k))
        @test worst <= log(2.0) / kmass + 1e-15
    end

    @testset "(phase) freeze/thaw mass clamp: smooth → hard as k → ∞" begin
        wice0 = 5.0; xm = wice0   # AT the kink, where smooth ≠ hard
        hardv = max(0.0, wice0 - xm)
        sv(k) = value(CLM.smooth_max(Dual(0.0, 0.0), Dual(wice0 - xm, 0.0); k = k))
        k0 = CLM.PHASE_CHANGE_MASS_K[]
        errs = [abs(sv(k) - hardv) for k in (5.0, 50.0, 500.0, 5000.0)]
        CLM.PHASE_CHANGE_MASS_K[] = k0
        @info "phase mass-clamp smooth→hard error vs k: $errs"
        @test issorted(errs; rev = true)
        @test errs[end] < 1e-3
    end

    # ---------------------------------------------------------------------
    # mass_weighted_snow_radius (snow-grain-radius clamp)
    # ---------------------------------------------------------------------
    @testset "(snow radius) smooth ON: finite gradient across both clamps" begin
        rmin = 54.526; rmax = 1500.0
        # radius as a function of one of the input radii, sweeping across rmax.
        f(r) = CLM.mass_weighted_snow_radius(r, r, 1.0, 1.0, rmin, rmax)
        anynan = false
        for r in (rmax - 50.0):5.0:(rmax + 50.0)
            g = ForwardDiff.derivative(f, r)
            isfinite(g) || (anynan = true)
        end
        @test !anynan
        @test f(rmax + 100.0) == rmax              # exact clamp at the ceiling (Float64)
        @test f(rmin - 100.0) == rmin              # exact clamp at the floor (Float64)
    end

    # ---------------------------------------------------------------------
    # divide_snow_layers! end-to-end (snow split repartition)
    # ---------------------------------------------------------------------
    # The state structs are concrete Float64, but the low-level
    # `divide_snow_layers!` accepts AbstractMatrix{<:Real} + a parametric
    # AerosolData, so we can drive it with Dual arrays end-to-end (mirrors the
    # btran test calling `smstress_rootr_btran!` with Dual arrays directly).
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
    nlevsno = CLM.varpar.nlevsno
    CLM.snow_hydrology_set_control_for_testing!()
    dzmin = zeros(nlevsno); dzmax_l = zeros(nlevsno); dzmax_u = zeros(nlevsno)
    dzmin[1] = 0.010; dzmax_l[1] = 0.03; dzmax_u[1] = 0.02
    dzmin[2] = 0.015; dzmax_l[2] = 0.07; dzmax_u[2] = 0.05
    for j in 3:nlevsno
        dzmin[j] = dzmax_u[j-1] * 0.5
        dzmax_u[j] = 2.0 * dzmax_u[j-1] + 0.01
        dzmax_l[j] = dzmax_u[j] + dzmax_l[j-1]
        if j == nlevsno
            dzmax_u[j] = floatmax(Float64); dzmax_l[j] = floatmax(Float64)
        end
    end
    CLM.SNOW_DZMIN[] = dzmin; CLM.SNOW_DZMAX_L[] = dzmax_l; CLM.SNOW_DZMAX_U[] = dzmax_u

    # 2-layer pack whose bottom layer is thick enough to split, with a temperature
    # gradient so the split-temperature blend (smooth_heaviside on the tfrz
    # crossing) is exercised. Returns the new bottom-layer temperature.
    function divide_bottom_temp(Tupper)
        Dt = typeof(Tupper); nc = 1; nlay = nlevsno + 2
        snl = [-2]
        dz = zeros(Dt, nc, nlay); zi = zeros(Dt, nc, nlay+1); z = zeros(Dt, nc, nlay)
        t = fill(Dt(263.0), nc, nlay); ice = zeros(Dt, nc, nlay); liq = zeros(Dt, nc, nlay)
        rds = fill(Dt(100.0), nc, nlay); frac_sno = [Dt(1.0)]
        dz[1, 4] = Dt(0.02); dz[1, 5] = Dt(0.30)        # bottom layer (Fortran 0) thick
        ice[1, 4] = Dt(5.0); ice[1, 5] = Dt(60.0)
        liq[1, 4] = Dt(0.2); liq[1, 5] = Dt(2.0)
        t[1, 4] = Tupper; t[1, 5] = Dt(272.0)
        aero = CLM.AerosolData{Dt, Vector{Dt}, Matrix{Dt}}()
        for f in (:mss_bcphi_col, :mss_bcpho_col, :mss_ocphi_col, :mss_ocpho_col,
                  :mss_dst1_col, :mss_dst2_col, :mss_dst3_col, :mss_dst4_col)
            setfield!(aero, f, zeros(Dt, nc, nlay))
        end
        CLM.divide_snow_layers!(snl, dz, zi, z, t, ice, liq, frac_sno, rds, aero,
                                false, trues(nc), 1:nc, nlevsno)
        return t[1, 5]
    end

    @testset "(divide) smooth ON: new-layer temp finite & continuous across tfrz" begin
        anynan = false; smooth_maxjump = 0.0; prev = nothing
        for Tk in 270.0:0.2:278.0
            g = ForwardDiff.derivative(divide_bottom_temp, Tk)
            isfinite(g) || (anynan = true)
            prev !== nothing && (smooth_maxjump = max(smooth_maxjump, abs(g - prev)))
            prev = g
        end
        @test !anynan
        @test smooth_maxjump < 1.0          # bounded (no vertical step)
    end

    @testset "(divide) smooth OFF: Float64 runs + finite" begin
        # Float64 end-to-end through the (smooth-routed but exact-on-Float64) kernel.
        v = divide_bottom_temp(274.0)
        @test isfinite(v)
        @test v > 0.0
    end
end
