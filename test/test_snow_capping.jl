# Tests for snow capping (SnowCapping in SnowHydrologyMod.F90).
#
# Snow capping removes excess snow mass from the bottom snow layer for columns
# whose total snow water exceeds h2osno_max. The removed mass leaves the snow
# pack as the qflx_snwcp_* fluxes (runoff) or qflx_snwcp_discarded_* (reset).
# These tests verify the end-to-end mass conservation of the capping sequence:
#   bottom-layer mass removed == capping fluxes * dtime.

@testset "Snow capping" begin

    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
    nlevsno = CLM.varpar.nlevsno

    dtime = 1800.0
    nc = 3
    bounds = 1:nc

    # Three columns:
    #   c=1: snow below h2osno_max  -> no capping
    #   c=2: snow above h2osno_max  -> capping, runoff
    #   c=3: snow well above max     -> capping, runoff
    h2osno_max = CLM.H2OSNO_MAX

    @testset "mass conservation (capping -> runoff)" begin
        # Bottom-layer ice/liquid that hold the (large) snow mass.
        ice_bottom = [10.0, h2osno_max + 500.0, h2osno_max + 5000.0]
        liq_bottom = [2.0, 50.0, 200.0]
        # Total snow water = bottom layer here (single-layer columns).
        h2osno_total = ice_bottom .+ liq_bottom

        # Save originals for conservation check.
        ice0 = copy(ice_bottom)
        liq0 = copy(liq_bottom)

        dz_bottom = [0.5, 1.5, 5.0]
        topo = zeros(nc)
        col_landunit = collect(1:nc)
        lun_itype = fill(CLM.ISTSOIL, nc)   # not glacier -> standard capping
        mask_snow = trues(nc)

        # Flux arrays initialised to zero.
        qflx_snwcp_ice = zeros(nc)
        qflx_snwcp_liq = zeros(nc)
        qflx_snwcp_discarded_ice = zeros(nc)
        qflx_snwcp_discarded_liq = zeros(nc)
        CLM.init_flux_snow_capping!(qflx_snwcp_ice, qflx_snwcp_liq,
            qflx_snwcp_discarded_ice, qflx_snwcp_discarded_liq, mask_snow, bounds)

        mask_capping = falses(nc)
        rho_orig_bottom = zeros(nc)
        frac_adjust = zeros(nc)

        CLM.bulk_flux_snow_capping_fluxes!(
            mask_capping, rho_orig_bottom, frac_adjust,
            qflx_snwcp_ice, qflx_snwcp_liq,
            qflx_snwcp_discarded_ice, qflx_snwcp_discarded_liq,
            dtime, dz_bottom, topo, h2osno_total,
            ice_bottom, liq_bottom,
            col_landunit, lun_itype,
            mask_snow, bounds, nlevsno, 100)

        # c=1 has no excess -> not capping; c=2,c=3 capping.
        @test mask_capping == [false, true, true]
        @test all(rho_orig_bottom[2:3] .> 0.0)
        # Fraction remaining is in (0,1) where capping happened.
        @test all(0.0 .< frac_adjust[2:3] .< 1.0)

        # All capping went to runoff (apply_runoff = true), none discarded.
        @test all(qflx_snwcp_discarded_ice .== 0.0)
        @test all(qflx_snwcp_discarded_liq .== 0.0)

        # Now remove the capped mass from the bottom layer.
        CLM.update_state_remove_snow_capping_fluxes!(
            ice_bottom, liq_bottom, dtime,
            qflx_snwcp_ice, qflx_snwcp_liq,
            qflx_snwcp_discarded_ice, qflx_snwcp_discarded_liq,
            mask_capping, bounds)

        # ---- Mass conservation: removed mass == flux * dtime ----
        for c in 1:nc
            removed_ice = ice0[c] - ice_bottom[c]
            removed_liq = liq0[c] - liq_bottom[c]
            @test removed_ice ≈ (qflx_snwcp_ice[c] + qflx_snwcp_discarded_ice[c]) * dtime atol=1e-9
            @test removed_liq ≈ (qflx_snwcp_liq[c] + qflx_snwcp_discarded_liq[c]) * dtime atol=1e-9
        end

        # Untouched column conserved exactly.
        @test ice_bottom[1] == ice0[1]
        @test liq_bottom[1] == liq0[1]

        # Capped columns: mass actually decreased and stays non-negative.
        @test ice_bottom[2] < ice0[2]
        @test ice_bottom[3] < ice0[3]
        @test all(ice_bottom .>= 0.0)
        @test all(liq_bottom .>= 0.0)

        # The total snow water removed should not exceed the original excess
        # (we cap at most (1 - min_snow_to_keep) of the bottom layer).
        for c in 2:3
            removed_tot = (ice0[c] - ice_bottom[c]) + (liq0[c] - liq_bottom[c])
            excess = h2osno_total[c] - h2osno_max
            @test removed_tot <= excess + 1e-6
            @test removed_tot > 0.0
        end
    end

    @testset "dz and aerosol scaling after capping" begin
        # Single capping column to check dz + aerosol scaling.
        nc2 = 1
        bounds2 = 1:nc2
        h2osno_total = [h2osno_max + 1000.0]
        ice_bottom = [h2osno_max + 1000.0]
        liq_bottom = [0.0]
        dz_bottom = [3.0]
        topo = [0.0]
        col_landunit = [1]
        lun_itype = [CLM.ISTSOIL]
        mask_snow = trues(nc2)

        qflx_snwcp_ice = zeros(nc2)
        qflx_snwcp_liq = zeros(nc2)
        qflx_snwcp_discarded_ice = zeros(nc2)
        qflx_snwcp_discarded_liq = zeros(nc2)
        mask_capping = falses(nc2)
        rho_orig_bottom = zeros(nc2)
        frac_adjust = zeros(nc2)

        CLM.bulk_flux_snow_capping_fluxes!(
            mask_capping, rho_orig_bottom, frac_adjust,
            qflx_snwcp_ice, qflx_snwcp_liq,
            qflx_snwcp_discarded_ice, qflx_snwcp_discarded_liq,
            dtime, dz_bottom, topo, h2osno_total,
            ice_bottom, liq_bottom,
            col_landunit, lun_itype,
            mask_snow, bounds2, nlevsno, 100)

        @test mask_capping[1]
        ice0 = copy(ice_bottom)
        CLM.update_state_remove_snow_capping_fluxes!(
            ice_bottom, liq_bottom, dtime,
            qflx_snwcp_ice, qflx_snwcp_liq,
            qflx_snwcp_discarded_ice, qflx_snwcp_discarded_liq,
            mask_capping, bounds2)

        # Build an aerosol struct with a known bottom-layer mass.
        jj_bottom = nlevsno
        aer = CLM.AerosolData()
        ncol = nc2
        nlay = jj_bottom
        aer.mss_bcphi_col = zeros(ncol, nlay); aer.mss_bcphi_col[1, jj_bottom] = 5.0
        aer.mss_bcpho_col = zeros(ncol, nlay); aer.mss_bcpho_col[1, jj_bottom] = 6.0
        aer.mss_ocphi_col = zeros(ncol, nlay); aer.mss_ocphi_col[1, jj_bottom] = 7.0
        aer.mss_ocpho_col = zeros(ncol, nlay); aer.mss_ocpho_col[1, jj_bottom] = 8.0
        aer.mss_dst1_col  = zeros(ncol, nlay); aer.mss_dst1_col[1, jj_bottom]  = 1.0
        aer.mss_dst2_col  = zeros(ncol, nlay); aer.mss_dst2_col[1, jj_bottom]  = 2.0
        aer.mss_dst3_col  = zeros(ncol, nlay); aer.mss_dst3_col[1, jj_bottom]  = 3.0
        aer.mss_dst4_col  = zeros(ncol, nlay); aer.mss_dst4_col[1, jj_bottom]  = 4.0

        dz0 = dz_bottom[1]
        CLM.snow_capping_update_dz_and_aerosols!(
            dz_bottom, aer, jj_bottom,
            rho_orig_bottom, ice_bottom, frac_adjust,
            mask_capping, bounds2)

        # dz should be scaled to conserve ice density: dz = ice/rho_orig.
        @test dz_bottom[1] ≈ ice_bottom[1] / rho_orig_bottom[1] atol=1e-9
        @test dz_bottom[1] < dz0   # thinner after removing mass

        # Aerosols scaled by frac_adjust (mass-concentration conserving).
        fa = frac_adjust[1]
        @test aer.mss_bcphi_col[1, jj_bottom] ≈ 5.0 * fa atol=1e-12
        @test aer.mss_dst4_col[1, jj_bottom]  ≈ 4.0 * fa atol=1e-12
    end
end
