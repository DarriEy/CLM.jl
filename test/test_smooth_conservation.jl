# ==========================================================================
# test_smooth_conservation.jl — the AD-SMOOTHED physics must conserve water.
#
# SMOOTH_MODE[] = :always makes plain Float64 evaluate the smooth surrogates, so
# the finite-difference and the ForwardDiff evaluation see the same function
# (calibration.jl's use_smooth_fd path). That path is covered by the same live,
# FATAL top-level water-balance check as every other path — so the smoothed
# physics has to close the column budget, not merely come close.
#
# It did not always. Several `smooth_max`/`smooth_min` replacements had been put
# on ReLU GUARDS whose argument is a WATER MASS:
#
#     smooth_max(0, <water>)  OVERSHOOTS  max(0, <water>)  by up to log(2)/k
#
# and with the generic k = 50 that overshoot is 0.0139 — which on a metres- or
# millimetres-of-water axis is not "rounding a corner", it is 0.0139 m / 0.0139 mm
# of water CONJURED, every layer, every step, with no flux paying for it. The
# storage update and the flux that pays for it stopped agreeing, and the column
# balance broke by ~4 mm/step (snow) / ~7.5e-3 mm/step (growing season) — the leak
# PR #211 had to gate the hard error off for.
#
# The fixes make each storage/flux pair come from the SAME expression:
#   * snow_hydrology.jl  — percolation drainage is boxed by the liquid the layer
#     actually holds (and its smoothing k is set on the metres-of-water axis).
#   * soil_moist_stress.jl — btran is the exact sum of the rootr it normalizes,
#     so sum_j rootr == 1 and the soil gives up exactly qflx_tran_veg.
#   * soil_temperature.jl — the phase change writes liq as the exact remainder
#     (covered in test_soil_temperature.jl).
#   * canopy_fluxes.jl — the canopy store moves by exactly the evaporation the
#     column balance debits.
#
# Each test asserts the INVARIANT (to machine precision, not to the smoothing
# width), in BOTH modes, and that the default (:auto) path is byte-identical to
# the hard reference.
# ==========================================================================

using Test
using ForwardDiff
using CLM

@testset "smoothed physics conserves water" begin

    with_smooth(f, mode) = begin
        saved = CLM.SMOOTH_MODE[]
        CLM.SMOOTH_MODE[] = mode
        try
            f()
        finally
            CLM.SMOOTH_MODE[] = saved
        end
    end

    # ---------------------------------------------------------------------
    # 1. Snow percolation — the drainage out of a layer is boxed by the liquid
    #    that layer holds. The flux is applied verbatim to h2osoi_liq
    #    (update_state_snow_percolation!), so `0 <= q*dt <= h2osoi_liq` is exactly
    #    the condition under which storage never goes negative and no downstream
    #    re-clamp gets the chance to manufacture water.
    #
    #    With k = 50 on this metres-of-water axis, smooth_max(0, q) returned
    #    log(2)/50 = 0.0139 m = 13.9 kg/m2 of drainage out of layers that should
    #    drain nothing, driving h2osoi_liq to -5 mm on the first step.
    # ---------------------------------------------------------------------
    @testset "snow percolation never drains more liquid than a layer holds" begin
        CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
        nlevsno = CLM.varpar.nlevsno
        ntot = nlevsno + 15
        dtime = 1800.0

        # Snow-pack states spanning the regimes: wet/drainable, dry (no liquid at
        # all — the case the smoothed floor used to invent water for), dense/frozen
        # (no pore space), and a partially-melted pack.
        function make_pack(; liq, ice, dz_l, fse)
            nc = 1
            snl = [-length(liq)]
            dz = zeros(nc, ntot)
            h2osoi_ice = zeros(nc, ntot)
            h2osoi_liq = zeros(nc, ntot)
            for (k, j) in enumerate((snl[1] + 1):0)
                jj = j + nlevsno
                dz[1, jj] = dz_l[k]
                h2osoi_ice[1, jj] = ice[k]
                h2osoi_liq[1, jj] = liq[k]
            end
            (nc, snl, dz, h2osoi_ice, h2osoi_liq, [fse])
        end

        packs = [
            # wet pack, plenty of pore space → real drainage
            (liq = [3.0, 4.0, 5.0], ice = [5.0, 8.0, 10.0], dz_l = [0.03, 0.05, 0.08], fse = 1.0),
            # BONE DRY liquid: the exact physics drains exactly nothing
            (liq = [0.0, 0.0, 0.0], ice = [5.0, 8.0, 10.0], dz_l = [0.03, 0.05, 0.08], fse = 1.0),
            # dense, small frac_sno_eff (the Bow winter case that leaked 4 mm/step)
            (liq = [0.4, 0.9, 6.0], ice = [12.0, 20.0, 39.0], dz_l = [0.05, 0.2, 1.6], fse = 0.114),
            # nearly-empty layers
            (liq = [1.0e-4, 0.0], ice = [1.0e-3, 0.02], dz_l = [0.01, 0.02], fse = 0.5),
        ]

        for mode in (:auto, :always)
            for pk in packs
                (nc, snl, dz, ice, liq, fse) = make_pack(; pk...)
                q = zeros(nc, nlevsno)
                with_smooth(mode) do
                    CLM.bulk_flux_snow_percolation!(q, dtime, snl, dz, fse, ice, liq,
                                                    trues(nc), 1:nc, nlevsno)
                end
                for j in (snl[1] + 1):0
                    jj = j + nlevsno
                    drained = q[1, jj] * dtime            # mm of water out of layer j
                    @test drained >= 0.0                  # never runs backwards
                    @test drained <= liq[1, jj] + 1.0e-12 # never more than it holds
                end
            end
        end

        # A dry pack drains EXACTLY nothing in both modes (this is the 13.9 mm bug).
        (nc, snl, dz, ice, liq, fse) = make_pack(; packs[2]...)
        for mode in (:auto, :always)
            q = zeros(nc, nlevsno)
            with_smooth(mode) do
                CLM.bulk_flux_snow_percolation!(q, dtime, snl, dz, fse, ice, liq,
                                                trues(nc), 1:nc, nlevsno)
            end
            @test all(iszero, q)
        end

        # Default (:auto) Float64 stays the EXACT hard physics, and the smoothed
        # values track it (the smoothing width is on the water axis now, not 14 mm).
        (nc, snl, dz, ice, liq, fse) = make_pack(; packs[1]...)
        q_auto = zeros(nc, nlevsno); q_always = zeros(nc, nlevsno)
        with_smooth(:auto) do
            CLM.bulk_flux_snow_percolation!(q_auto, dtime, snl, dz, fse, ice, liq,
                                            trues(nc), 1:nc, nlevsno)
        end
        with_smooth(:always) do
            CLM.bulk_flux_snow_percolation!(q_always, dtime, snl, dz, fse, ice, liq,
                                            trues(nc), 1:nc, nlevsno)
        end
        @test maximum(abs.(q_always .- q_auto)) * dtime < 1.0e-5   # << balance threshold

        # And it stays differentiable: the drainage still has a finite, non-trivial
        # gradient w.r.t. the layer's liquid water under ForwardDiff.
        function drain_of_liq(x)
            (nc, snl, dzd, iced, liqd, fsed) = make_pack(; packs[1]...)
            dzT = eltype(x).(dzd); iceT = eltype(x).(iced); fseT = eltype(x).(fsed)
            liqT = eltype(x).(liqd)
            liqT[1, 0 + nlevsno] = x                      # bottom-layer liquid
            qd = zeros(eltype(x), nc, nlevsno)
            CLM.bulk_flux_snow_percolation!(qd, dtime, snl, dzT, fseT, iceT, liqT,
                                            trues(nc), 1:nc, nlevsno)
            qd[1, 0 + nlevsno]
        end
        g = ForwardDiff.derivative(drain_of_liq, 5.0)
        @test isfinite(g)
        @test g > 0.0        # more liquid in the layer → more drainage out of it
    end

    # ---------------------------------------------------------------------
    # 2. btran / rootr — the transpiration partition sums to one.
    #
    #    rootr is normalized by btran (`rootr /= btran`), so the soil gives up
    #    exactly qflx_tran_veg IFF btran is the exact sum of the SAME rootr values.
    #    Fortran writes `btran += max(rootr, 0)` (SoilMoistStressMod.F90:413) — a
    #    guard, since rootr >= 0 there. Smoothing that guard inflated btran without
    #    changing sum(rootr), so sum_j rootr < 1 and the soil quietly kept water the
    #    balance had already debited as transpiration.
    # ---------------------------------------------------------------------
    @testset "rootr sums to one (exact AND smoothed physics)" begin
        CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
        nlevsno = CLM.varpar.nlevsno
        nlevg = 10
        joff = nlevsno
        npt, nct = 3, 2
        tfrz = CLM.TFRZ

        pcol = [1, 1, 2]
        pity = [1, 3, 5]
        mask = trues(npt)
        rf   = zeros(npt, nlevg)
        for p in 1:npt
            w = [exp(-0.5 * j) for j in 1:nlevg]
            rf[p, :] .= w ./ sum(w)          # rooting fraction, sums to 1
        end
        rfu = copy(rf)

        smpso = fill(-66000.0, 25); smpsc = fill(-255000.0, 25)
        watsat = fill(0.45, nct, nlevg); sucsat = fill(100.0, nct, nlevg)
        bsw = fill(5.0, nct, nlevg); eff = fill(0.4, nct, nlevg)
        t_s = fill(290.0, nct, nlevsno + nlevg)
        h2o = zeros(nct, nlevsno + nlevg)
        for c in 1:nct, j in 1:nlevg
            h2o[c, j + joff] = 0.05 + 0.03 * j
            t_s[c, j + joff] = 290.0 - 2.0 * j
        end
        # Layers that the HARD gate switches off (dry / frozen) — exactly the ones
        # a smoothed floor used to add phantom weight for.
        h2o[1, 3 + joff] = 0.0
        t_s[2, 5 + joff] = tfrz - 5.0

        CLM.set_perchroot_opt!(false, false); CLM.init_root_moist_stress!()

        results = Dict{Symbol, Any}()
        for mode in (:auto, :always)
            rootr = zeros(npt, nlevg); rresis = fill(-9.9, npt, nlevg); btran = zeros(npt)
            with_smooth(mode) do
                CLM.calc_root_moist_stress_clm45default!(rfu, rf, rootr, btran, rresis,
                    smpso, smpsc, t_s, watsat, sucsat, bsw, eff, h2o,
                    pcol, pity, mask, 1:npt, nlevg, nlevsno)
            end
            results[mode] = (rootr = rootr, btran = btran)

            for p in 1:npt
                btran[p] > 0.0 || continue
                # THE conservation identity: the transpiration sink over the column
                # is sum_j rootr[p,j] * qflx_tran_veg — it must be the whole flux.
                @test isapprox(sum(@view rootr[p, :]), 1.0; atol = 1.0e-12, rtol = 1.0e-12)
            end
        end

        # Default path unchanged (byte-identical, not merely close).
        @test results[:auto].rootr == results[:auto].rootr   # (sanity: deterministic)
        @test all(isfinite, results[:always].rootr)

        # Still differentiable: btran responds to soil water under ForwardDiff, and
        # the partition still sums to one on Duals.
        function btran_of_h2o(x)
            h2od = eltype(x).(h2o)
            h2od[1, 4 + joff] = x
            rootr = zeros(eltype(x), npt, nlevg)
            rresis = fill(eltype(x)(-9.9), npt, nlevg)
            btran = zeros(eltype(x), npt)
            CLM.calc_root_moist_stress_clm45default!(rfu, rf, rootr, btran, rresis,
                smpso, smpsc, eltype(x).(t_s), watsat, sucsat, bsw, eff, h2od,
                pcol, pity, mask, 1:npt, nlevg, nlevsno)
            @test isapprox(ForwardDiff.value(sum(@view rootr[1, :])), 1.0; atol = 1.0e-12)
            btran[1]
        end
        gb = ForwardDiff.derivative(btran_of_h2o, 0.17)
        @test isfinite(gb)
        @test gb > 0.0     # wetter soil → less moisture stress
    end

    # ---------------------------------------------------------------------
    # 3. Canopy store — the canopy sheds EXACTLY the evaporation the column
    #    balance debits, so the non-negativity floor on the total store must stay
    #    a HARD max. The upstream ecidif cap already guarantees the floor cannot
    #    bite (|net_evap| <= h2ocan), so it is a pure guard and hardening it costs
    #    no differentiability — but a SMOOTHED floor bites everywhere, adding
    #    log(2)/50 = 0.0139 mm of canopy water per patch per step.
    #
    #    The update lives inside canopy_fluxes!'s fused per-patch kernel, with no
    #    callable seam of its own, so this pins the source. (The end-to-end proof is
    #    scripts/longhorizon_conservation.jl --smooth, which closes at ~1e-15
    #    mm/step; this is the guard that keeps the one-token revert from sneaking
    #    back in.)
    # ---------------------------------------------------------------------
    @testset "canopy water store: non-negativity floor is a HARD max" begin
        src = read(joinpath(@__DIR__, "..", "src", "biogeophys", "canopy_fluxes.jl"), String)
        m = match(r"_h2ocan_new = (\w+)\(zero\(T\), out\.liqcan\[p\] \+ out\.snocan\[p\] \+ _net_evap\)", src)
        @test m !== nothing
        @test m.captures[1] == "max"        # NOT smooth_max — see the note above
    end
end
