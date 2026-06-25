# ==========================================================================
# test_btran_smoothing.jl — AD-smoothing of the btran (root soil-moisture
# stress) discontinuities in SoilMoistStressMod (`_smstress_rootr_btran_kernel!`).
#
# The CLM4.5-default btran kernel has a HARD per-layer gate that zeroes a
# layer's root-water contribution the instant the layer goes dry
# (h2osoi_liqvol <= 0) or near-frozen (t_soisno <= tfrz-2), plus hard clamps on
# s_node / smp_node / rresis / the rootr floor. As temperature or soil water
# sweep across those thresholds, btran[p] — the factor that scales stomatal
# conductance / photosynthesis — has KINKS, so its forward-AD derivative jumps
# discontinuously.
#
# These tests demonstrate:
#   (a) Smoothing ON (ForwardDiff.Dual): d(btran)/dT and d(btran)/d(h2o) are
#       FINITE and CONTINUOUS across the gate (bounded jump between adjacent
#       points), whereas the hard physics has a near-vertical step there.
#   (b) Smoothing OFF (plain Float64): byte-identical to the hard reference.
#   (c) Physics error → 0 as the smoothing width ε → 0 (k → ∞).
# ==========================================================================

using Test
using ForwardDiff
using ForwardDiff: Dual, value
using Random
using CLM

@testset "btran AD-smoothing (soil moisture stress)" begin

    nlevsno = 5; nlevgrnd = 4; joff = nlevsno
    np = 1; nc = 1; nlevtot = nlevsno + nlevgrnd
    tfrz = CLM.TFRZ

    smpso = [-35000.0]; smpsc = [-275000.0]
    watsat = fill(0.45, nc, nlevgrnd)
    sucsat = fill(100.0, nc, nlevgrnd)
    bsw = fill(5.0, nc, nlevgrnd)
    eff_porosity = fill(0.40, nc, nlevgrnd)
    patch_column = [1]; patch_itype = [0]
    mask_patch = trues(np)
    rootfr = fill(1.0 / nlevgrnd, np, nlevgrnd)
    rootfr_unf = zeros(np, nlevgrnd)

    CLM.set_perchroot_opt!(false, false)
    CLM.init_root_moist_stress!()

    # btran as a function of layer-1 temperature (sweeps across the tfrz-2 gate).
    # The element type follows Tval: Float64 → hard physics, Dual → smooth.
    function btran_of_T(Tval)
        Dt = typeof(Tval)
        rootr = zeros(Dt, np, nlevgrnd); rresis = zeros(Dt, np, nlevgrnd); btran = zeros(Dt, np)
        t = fill(Dt(290.0), nc, nlevtot); t[1, 1 + joff] = Tval
        h = zeros(Dt, nc, nlevtot)
        for j in 1:nlevgrnd; h[1, j + joff] = Dt(0.25); end
        CLM.smstress_rootr_btran!(rootr, rresis, btran, mask_patch, patch_column,
            patch_itype, Dt.(rootfr), Dt.(rootfr_unf), Dt.(smpso), Dt.(smpsc), t,
            Dt.(watsat), Dt.(sucsat), Dt.(bsw), Dt.(eff_porosity), h,
            nlevgrnd, joff, tfrz, false)
        return btran[1]
    end

    # btran as a function of layer-1 liquid water (sweeps across the >0 gate).
    function btran_of_H(Hval)
        Dt = typeof(Hval)
        rootr = zeros(Dt, np, nlevgrnd); rresis = zeros(Dt, np, nlevgrnd); btran = zeros(Dt, np)
        t = fill(Dt(290.0), nc, nlevtot)
        h = zeros(Dt, nc, nlevtot)
        for j in 1:nlevgrnd; h[1, j + joff] = Dt(0.25); end
        h[1, 1 + joff] = Hval
        CLM.smstress_rootr_btran!(rootr, rresis, btran, mask_patch, patch_column,
            patch_itype, Dt.(rootfr), Dt.(rootfr_unf), Dt.(smpso), Dt.(smpsc), t,
            Dt.(watsat), Dt.(sucsat), Dt.(bsw), Dt.(eff_porosity), h,
            nlevgrnd, joff, tfrz, false)
        return btran[1]
    end

    # ---------------------------------------------------------------------
    # (a) Smoothing ON: AD gradient is finite & continuous across the gate
    # ---------------------------------------------------------------------
    @testset "smooth ON: d(btran)/dT finite & continuous across freeze gate" begin
        gate = tfrz - 2.0

        # HARD reference: FD of the Float64 physics shows a near-vertical step.
        hard_maxjump = 0.0; prev = nothing
        for Tk in (gate - 0.3):0.05:(gate + 0.3)
            hfd = 1e-6
            d = (btran_of_T(Tk + hfd) - btran_of_T(Tk - hfd)) / (2hfd)
            prev !== nothing && (hard_maxjump = max(hard_maxjump, abs(d - prev)))
            prev = d
        end

        # SMOOTH: ForwardDiff derivative is finite everywhere and only mildly varying.
        smooth_maxjump = 0.0; prev2 = nothing
        for Tk in (gate - 0.3):0.05:(gate + 0.3)
            g = ForwardDiff.derivative(btran_of_T, Tk)
            @test isfinite(g)                  # no NaN/Inf across the discontinuity
            prev2 !== nothing && (smooth_maxjump = max(smooth_maxjump, abs(g - prev2)))
            prev2 = g
        end

        @info "btran/T gate: hard max |Δ(dbtran/dT)| = $hard_maxjump, smooth = $smooth_maxjump"
        # The smooth derivative jump is bounded and orders of magnitude below the hard step.
        @test smooth_maxjump < 1.0
        @test hard_maxjump > 100.0 * max(smooth_maxjump, 1e-12)
    end

    @testset "smooth ON: d(btran)/d(h2o) finite & continuous across dry gate" begin
        smooth_maxjump = 0.0; prev = nothing; anynan = false
        for Hk in 0.0:0.005:0.06
            g = ForwardDiff.derivative(btran_of_H, Hk)
            isfinite(g) || (anynan = true)
            prev !== nothing && (smooth_maxjump = max(smooth_maxjump, abs(g - prev)))
            prev = g
        end
        @test !anynan                            # finite even AT and across h2o=0
        @test smooth_maxjump < 1e4               # bounded (no vertical step)
    end

    # ---------------------------------------------------------------------
    # (b) Smoothing OFF (Float64): byte-identical to the hard reference
    # ---------------------------------------------------------------------
    @testset "smooth OFF: Float64 byte-identical to hard physics" begin
        nlevg = 6; nct = 3; npt = 3; jo = nlevsno; nt = nlevsno + nlevg
        smpso_b = [-35000.0, -83000.0]; smpsc_b = [-275000.0, -428000.0]
        watsat_b = fill(0.45, nct, nlevg); sucsat_b = fill(100.0, nct, nlevg)
        bsw_b = fill(5.0, nct, nlevg); eff_b = fill(0.40, nct, nlevg)
        pcol = [1, 2, 3]; pity = [0, 1, 0]; mask_b = trues(npt)
        rf = fill(1.0 / nlevg, npt, nlevg); rfu = fill(0.5 / nlevg, npt, nlevg)

        rng = MersenneTwister(7)
        t_b = fill(290.0, nct, nt); h_b = zeros(nct, nt)
        for c in 1:nct, j in 1:nlevg
            t_b[c, j + jo] = 290.0 - 25.0 * rand(rng)
            h_b[c, j + jo] = j == 2 ? 0.0 : 0.05 + 0.4 * rand(rng)
        end
        # exact-boundary cells: t == tfrz-2 and h2o == 0
        t_b[1, 1 + jo] = tfrz - 2.0

        function hard_ref()
            rootr = zeros(npt, nlevg); rresis = fill(-9.9, npt, nlevg); btran = zeros(npt)
            for p in 1:npt
                c = pcol[p]; it = pity[p] + 1; acc = 0.0
                for j in 1:nlevg
                    if h_b[c, j + jo] <= 0.0 || t_b[c, j + jo] <= tfrz - 2.0
                        rootr[p, j] = 0.0
                    else
                        s = max(h_b[c, j + jo] / eff_b[c, j], 0.01); s = min(s, 1.0)
                        smp = -sucsat_b[c, j] * s^(-bsw_b[c, j]); smp = max(smpsc_b[it], smp)
                        rresis[p, j] = min((eff_b[c, j] / watsat_b[c, j]) *
                            (smp - smpsc_b[it]) / (smpso_b[it] - smpsc_b[it]), 1.0)
                        rootr[p, j] = rf[p, j] * rresis[p, j]
                        acc += max(rootr[p, j], 0.0)
                    end
                end
                btran[p] += acc
            end
            for p in 1:npt, j in 1:nlevg
                rootr[p, j] = btran[p] > 0.0 ? rootr[p, j] / btran[p] : 0.0
            end
            return rootr, rresis, btran
        end

        function via_code()
            rootr = zeros(npt, nlevg); rresis = fill(-9.9, npt, nlevg); btran = zeros(npt)
            CLM.calc_root_moist_stress_clm45default!(rfu, rf, rootr, btran, rresis,
                smpso_b, smpsc_b, t_b, watsat_b, sucsat_b, bsw_b, eff_b, h_b,
                pcol, pity, mask_b, 1:npt, nlevg, nlevsno)
            return rootr, rresis, btran
        end

        CLM.set_perchroot_opt!(false, false); CLM.init_root_moist_stress!()
        r1, re1, b1 = hard_ref(); r2, re2, b2 = via_code()
        @test r1 == r2     # byte-identical (===, not approx)
        @test re1 == re2
        @test b1 == b2
    end

    # ---------------------------------------------------------------------
    # (c) Physics error → 0 as the smoothing width ε → 0 (k → ∞)
    # ---------------------------------------------------------------------
    @testset "convergence: smooth → hard as k → ∞" begin
        Tk = 280.0                       # gate-on region; both gates ≈ 1
        hardv = btran_of_T(Tk)           # Float64 hard value
        smoothval(t) = value(btran_of_T(Dual(t, 0.0)))  # smooth value at zero partial

        # Save defaults, sweep all smoothing widths together, restore.
        k0c = CLM.BTRAN_SMOOTH_K[]; k0w = CLM.BTRAN_GATE_K_WATER[]; k0t = CLM.BTRAN_GATE_K_TEMP[]
        errs = Float64[]
        for k in (50.0, 200.0, 1000.0, 5000.0)
            CLM.BTRAN_SMOOTH_K[] = k
            CLM.BTRAN_GATE_K_WATER[] = k * 4
            CLM.BTRAN_GATE_K_TEMP[] = k / 10
            push!(errs, abs(smoothval(Tk) - hardv))
        end
        CLM.BTRAN_SMOOTH_K[] = k0c; CLM.BTRAN_GATE_K_WATER[] = k0w; CLM.BTRAN_GATE_K_TEMP[] = k0t

        @info "btran smooth→hard error vs k: $errs"
        @test issorted(errs; rev = true)   # monotonically decreasing
        @test errs[end] < 1e-6             # → 0 as k → ∞
    end
end
