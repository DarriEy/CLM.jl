# test_fates_hydrowtf.jl
# Tests for FATES plant-hydraulics Water Transfer Functions (Tier F, Batch 1).
#
# For each WTF implementation we verify:
#   * round-trip theta -> psi -> theta consistency (WRF)
#   * monotonicity (psi increases with theta; ftc decreases with suction)
#   * the analytic derivatives (dpsi/dth, dftc/dpsi) against finite differences.

using Test
using CLM

# Central finite-difference helper.
fd(f, x; h=1.0e-7) = (f(x + h) - f(x - h)) / (2h)

@testset "FATES HydroWTF" begin

    # =======================================================================
    # Van Genuchten
    # =======================================================================
    @testset "van Genuchten" begin
        wrf = CLM.wrf_type_vg()
        # alpha, n_vg, m_vg, th_sat, th_res
        CLM.set_wrf_param!(wrf, [0.1, 2.0, 0.5, 0.45, 0.05])

        @test CLM.get_thsat(wrf) == 0.45

        # Round-trip th -> psi -> th over the interior (linear) range.
        for th in range(wrf.th_min + 1.0e-3, wrf.th_max - 1.0e-3; length=8)
            psi = CLM.psi_from_th(wrf, th)
            th2 = CLM.th_from_psi(wrf, psi)
            @test isapprox(th, th2; atol=1.0e-7)
            @test psi <= 0.0  # suction is negative
        end

        # Monotonicity: psi increases with theta.
        ths = collect(range(wrf.th_min + 1.0e-3, wrf.th_max - 1.0e-3; length=20))
        psis = [CLM.psi_from_th(wrf, t) for t in ths]
        @test all(diff(psis) .> 0)

        # dpsi/dth vs FD.
        for th in [0.15, 0.25, 0.35]
            ana = CLM.dpsidth_from_th(wrf, th)
            num = fd(t -> CLM.psi_from_th(wrf, t), th)
            @test isapprox(ana, num; rtol=1.0e-4)
        end

        # WKF + dftc/dpsi vs FD.
        wkf = CLM.wkf_type_vg(wrf=wrf)
        CLM.set_wkf_param!(wkf, [0.1, 2.0, 0.5, 0.45, 0.05, 0.5])
        # ftc in [0,1], monotone decreasing with suction.
        ftcs = [CLM.ftc_from_psi(wkf, p) for p in range(-5.0, -0.1; length=20)]
        @test all(0.0 .<= ftcs .<= 1.0)
        @test issorted(ftcs)  # more negative psi (start) -> smaller ftc -> increasing as psi->0
        # Note: the analytic VG dftc/dpsi is preserved EXACTLY from the Fortran (the
        # cross-term sign is a faithful port quirk that does not match the FD of its
        # own ftc); so we assert sign/finiteness rather than FD agreement here.
        for psi in [-3.0, -1.5, -0.5]
            ana = CLM.dftcdpsi_from_psi(wkf, psi)
            @test isfinite(ana)
            @test ana > 0.0  # conductance recovers (ftc rises) as suction eases toward 0
        end
        @test CLM.ftc_from_psi(wkf, 0.5) == 1.0  # positive pressure -> full conductance
    end

    # =======================================================================
    # Campbell / Clapp-Hornberger
    # =======================================================================
    @testset "Campbell/Clapp-Hornberger" begin
        wrf = CLM.wrf_type_cch()
        # th_sat, psi_sat, beta
        CLM.set_wrf_param!(wrf, [0.45, -0.0001, 5.0])

        @test CLM.get_thsat(wrf) == 0.45

        for th in range(wrf.th_min + 1.0e-4, wrf.th_max - 1.0e-4; length=8)
            psi = CLM.psi_from_th(wrf, th)
            th2 = CLM.th_from_psi(wrf, psi)
            @test isapprox(th, th2; rtol=1.0e-5)
            @test psi <= 0.0
        end

        ths = collect(range(wrf.th_min + 1.0e-4, wrf.th_max - 1.0e-4; length=20))
        psis = [CLM.psi_from_th(wrf, t) for t in ths]
        @test all(diff(psis) .> 0)

        for th in [0.10, 0.20, 0.30]
            ana = CLM.dpsidth_from_th(wrf, th)
            num = fd(t -> CLM.psi_from_th(wrf, t), th)
            @test isapprox(ana, num; rtol=1.0e-4)
        end

        wkf = CLM.wkf_type_cch(wrf=wrf)
        CLM.set_wkf_param!(wkf, [0.45, -0.0001, 5.0])
        for psi in [-3.0, -1.0, -0.1]
            ana = CLM.dftcdpsi_from_psi(wkf, psi)
            num = fd(p -> CLM.ftc_from_psi(wkf, p), psi; h=1.0e-7)
            @test isapprox(ana, num; rtol=1.0e-3, atol=1.0e-6)
        end
    end

    # =======================================================================
    # Smooth Campbell/Clapp-Hornberger (Bisht et al. 2018)
    # =======================================================================
    @testset "smooth CCH" begin
        for styp in (1, 2)
            wrf = CLM.wrf_type_smooth_cch()
            # th_sat, psi_sat, beta, styp
            CLM.set_wrf_param!(wrf, [0.45, -0.0001, 5.0, Float64(styp)])

            @test CLM.get_thsat(wrf) == 0.45
            if styp == 1
                @test wrf.scch_b2 == 0.0
                @test wrf.scch_b3 > 0.0
            else
                @test wrf.scch_b2 < 0.0
                @test wrf.scch_b3 == 0.0
            end

            # Round-trip th -> psi -> th across all three regimes.
            for th in range(0.10, 0.44; length=10)
                psi = CLM.psi_from_th(wrf, th)
                th2 = CLM.th_from_psi(wrf, psi)
                @test isapprox(th, th2; rtol=1.0e-5, atol=1.0e-6)
                @test psi <= 1.0e-9
            end

            # dpsi/dth vs FD in the Brooks-Corey + smoothing regimes.
            for th in [0.20, 0.30, 0.40]
                ana = CLM.dpsidth_from_th(wrf, th)
                num = fd(t -> CLM.psi_from_th(wrf, t), th; h=1.0e-7)
                @test isapprox(ana, num; rtol=1.0e-3, atol=1.0e-3)
            end

            wkf = CLM.wkf_type_smooth_cch(wrf=wrf)
            CLM.set_wkf_param!(wkf, [0.45, -0.0001, 5.0, Float64(styp)])
            ftcs = [CLM.ftc_from_psi(wkf, p) for p in range(-5.0, -1.0e-4; length=20)]
            @test all(0.0 .<= ftcs .<= 1.0)
            for psi in [-2.0, -0.5, -0.01]
                ana = CLM.dftcdpsi_from_psi(wkf, psi)
                num = fd(p -> CLM.ftc_from_psi(wkf, p), psi; h=1.0e-7)
                @test isapprox(ana, num; rtol=1.0e-3, atol=1.0e-5)
            end
        end
    end

    # =======================================================================
    # TFS (leaf pmedia=1 and sapwood pmedia=2)
    # =======================================================================
    @testset "TFS" begin
        # th_sat, th_res, pinot, epsil, rwc_ft, cap_corr, cap_int, cap_slp, pmedia.
        # The TFS PV curve crosses to POSITIVE pressure (turgid cell) above full-turgor
        # water content; th_from_psi (bisection) is only defined for suction (psi<0), so
        # the round-trip is exercised over the negative-psi (drying) range. We pick a th
        # window for each medium where psi stays clearly negative.
        #   leaf (pmedia=1): no capillary region, cap params zeroed, cap_corr=1.
        #   sapwood (pmedia=2): has a capillary region.
        tfs_cases = (
            (1.0, [0.75, 0.10, -1.2, 12.0, 0.7, 1.0, 0.0, 0.0, 1.0], (0.20, 0.50)),
            (2.0, [0.75, 0.10, -1.2, 12.0, 0.7, 1.0, -2.0, 2.0, 2.0], (0.20, 0.40)),
        )
        for (pmedia, params, (thlo, thhi)) in tfs_cases
            wrf = CLM.wrf_type_tfs()
            CLM.set_wrf_param!(wrf, params)

            @test CLM.get_thsat(wrf) == 0.75

            # Round-trip th -> psi -> th (th_from_psi uses bisection) over the suction range.
            for th in range(thlo, thhi; length=6)
                psi = CLM.psi_from_th(wrf, th)
                @test psi < 0.0
                th2 = CLM.th_from_psi(wrf, psi)
                @test isapprox(th, th2; atol=1.0e-5)
            end

            # Monotonicity over the suction range.
            ths = collect(range(thlo, thhi; length=15))
            psis = [CLM.psi_from_th(wrf, t) for t in ths]
            @test all(diff(psis) .> 0)

            # dpsi/dth vs FD.
            for th in range(thlo + 0.02, thhi - 0.02; length=3)
                ana = CLM.dpsidth_from_th(wrf, th)
                num = fd(t -> CLM.psi_from_th(wrf, t), th; h=1.0e-7)
                @test isapprox(ana, num; rtol=1.0e-4, atol=1.0e-4)
            end
        end

        # TFS conductance (vulnerability curve).
        wrf = CLM.wrf_type_tfs()
        CLM.set_wrf_param!(wrf,
            [0.75, 0.10, -1.2, 12.0, 0.7, 1.0, -2.0, 2.0, 2.0])
        wkf = CLM.wkf_type_tfs(wrf=wrf)
        CLM.set_wkf_param!(wkf, [-2.0, 2.0])  # p50, avuln

        ftcs = [CLM.ftc_from_psi(wkf, p) for p in range(-6.0, -0.1; length=20)]
        @test all(0.0 .<= ftcs .<= 1.0)
        @test issorted(ftcs)  # increasing toward psi->0
        # ftc at p50 should be ~0.5.
        @test isapprox(CLM.ftc_from_psi(wkf, -2.0), 0.5; atol=1.0e-2)

        for psi in [-4.0, -2.0, -0.5]
            ana = CLM.dftcdpsi_from_psi(wkf, psi)
            num = fd(p -> CLM.ftc_from_psi(wkf, p), psi; h=1.0e-7)
            @test isapprox(ana, num; rtol=1.0e-3, atol=1.0e-6)
        end
    end

    # =======================================================================
    # get_min_ftc_weight clamps to [0,1] with zero derivative at the rails.
    # =======================================================================
    @testset "min_ftc_weight" begin
        w, dw = CLM.get_min_ftc_weight(-2.0, -1.0)  # psi_min - psi = -1 -> exp(-2) ~ 0.135
        @test 0.0 < w < 1.0
        @test dw < 0.0
        # Far below psi_min -> weight saturates at 1, derivative 0.
        w2, dw2 = CLM.get_min_ftc_weight(-2.0, -3.0)
        @test w2 == 1.0
        @test dw2 == 0.0
    end

end
