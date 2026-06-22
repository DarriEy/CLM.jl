# test_fates_twostream.jl
# Tests for the FATES TwoStreamMLPEMod port (Tier F, Batch 1, radiation).
#
# Set up a simple single-layer, single-PFT canopy with known optical properties,
# run the two-stream solve, and assert energy conservation (absorbed + reflected
# + transmitted-to-ground ≈ incident) and sane bounds (all fractions in [0,1]).

using Test
using CLM

# Helper: build a fresh rad_params for a single PFT (pft index 1) over `nbands`
# broad bands with the given leaf/stem optical properties, then derive terms.
function setup_rad_params!(; nbands::Int=2,
                           rhol=0.10, taul=0.05, rhos=0.16, taus=0.001,
                           xl=0.01, clumping=1.0)
    n_pft = 1
    CLM.AllocateRadParams(n_pft, nbands)
    for ib in 1:nbands
        CLM.rad_params.rhol[ib, 1] = rhol
        CLM.rad_params.rhos[ib, 1] = rhos
        CLM.rad_params.taul[ib, 1] = taul
        CLM.rad_params.taus[ib, 1] = taus
    end
    CLM.rad_params.xl[1] = xl
    CLM.rad_params.clumping_index[1] = clumping
    CLM.RadParamPrep()
    return nothing
end

@testset "FATES TwoStreamMLPEMod" begin

    @testset "constants / API surface" begin
        @test CLM.normalized_upper_boundary == 1
        @test CLM.absolute_upper_boundary == 2
        @test CLM.air_ft == 0
        @test CLM.twostr_vis == 1
        @test CLM.twostr_nir == 2
        @test CLM.rel_err_thresh == 1.0e-6
        @test length(CLM.om_snow) == 2
    end

    @testset "RadParamPrep derives sane parameters" begin
        setup_rad_params!()
        @test CLM.rad_params.kd_leaf[1] > 0.0
        @test CLM.rad_params.kd_stem[1] == 1.0
        @test CLM.rad_params.avmu[1] > 0.0
        # om = rho + tau for leaves / stems
        @test CLM.rad_params.om_leaf[1, 1] ≈ 0.10 + 0.05
        @test CLM.rad_params.om_stem[1, 1] ≈ 0.16 + 0.001
        # phi2 = 0.877*(1 - 2*phi1)
        @test CLM.rad_params.phi2[1] ≈ 0.877 * (1.0 - 2.0 * CLM.rad_params.phi1[1])
    end

    @testset "single-layer single-PFT solve: energy conservation + bounds" begin
        setup_rad_params!()

        # One canopy layer, two columns (one vegetated PFT-1 element + an air
        # ghost element) so the layer area sums to exactly 1.0.
        ncan = 1
        ncol = 2
        ts = CLM.twostream_type()
        CLM.AllocInitTwoStream!(ts, [1, 2], ncan, ncol)
        ts.n_col[1] = 2
        ts.force_prep = true

        # Vegetated element
        veg = ts.scelg[1, 1]
        veg.pft = 1
        veg.area = 0.7
        veg.lai = 2.0
        veg.sai = 0.5

        # Air ghost element fills the rest of the footprint
        air = ts.scelg[1, 2]
        air.pft = CLM.air_ft
        air.area = 0.3
        air.lai = 0.0
        air.sai = 0.0

        CLM.GetNSCel!(ts)
        @test ts.n_scel == 2

        # Ground albedos
        for b in ts.band
            b.albedo_grnd_diff = 0.2
            b.albedo_grnd_beam = 0.2
        end

        CLM.CanopyPrep!(ts, 0.0)        # no snow
        CLM.ZenithPrep!(ts, 0.7)        # cos(zenith) = 0.7

        n_eq = 2 * ts.n_scel
        taulamb = zeros(Float64, n_eq)
        omega = zeros(Float64, n_eq, n_eq)
        ipiv = zeros(Int, n_eq)

        # The solver works on a normalized boundary condition: unit beam and unit
        # diffuse incident. The returned beam/diff fractions are each relative to
        # a unit incident in that radiation stream.
        Rbeam_atm = 1.0
        Rdiff_atm = 1.0
        ib = 1

        (albedo_beam, albedo_diff, consv_err,
         frac_abs_can_beam, frac_abs_can_diff,
         frac_beam_grnd_beam, frac_diff_grnd_beam, frac_diff_grnd_diff) =
            CLM.Solve!(ts, ib, CLM.absolute_upper_boundary,
                       Rbeam_atm, Rdiff_atm, taulamb, omega, ipiv)

        # Conservation error (already absolute, normalized by total incident)
        @test consv_err < CLM.rel_err_thresh

        # All component fluxes are finite
        for v in (albedo_beam, albedo_diff, frac_abs_can_beam, frac_abs_can_diff,
                  frac_beam_grnd_beam, frac_diff_grnd_beam, frac_diff_grnd_diff)
            @test isfinite(v)
        end

        # Albedos (reflected fractions) are non-negative and in [0,1] per stream
        @test 0.0 - CLM.rel_err_thresh <= albedo_beam <= 1.0 + CLM.rel_err_thresh
        @test 0.0 - CLM.rel_err_thresh <= albedo_diff <= 1.0 + CLM.rel_err_thresh

        # Canopy-absorbed fractions are non-negative and in [0,1]
        @test 0.0 - CLM.rel_err_thresh <= frac_abs_can_beam <= 1.0 + CLM.rel_err_thresh
        @test 0.0 - CLM.rel_err_thresh <= frac_abs_can_diff <= 1.0 + CLM.rel_err_thresh

        # Ground-incident fractions are non-negative and in [0,1]
        @test 0.0 - CLM.rel_err_thresh <= frac_beam_grnd_beam <= 1.0 + CLM.rel_err_thresh
        @test 0.0 - CLM.rel_err_thresh <= frac_diff_grnd_beam <= 1.0 + CLM.rel_err_thresh
        @test 0.0 - CLM.rel_err_thresh <= frac_diff_grnd_diff <= 1.0 + CLM.rel_err_thresh

        # Explicit per-stream energy balance: unit incident = reflected +
        # canopy-absorbed + ground-absorbed (= incident-at-ground * (1 - albedo)).
        beam_bal = albedo_beam + frac_abs_can_beam +
                   frac_diff_grnd_beam * (1.0 - ts.band[ib].albedo_grnd_diff) +
                   frac_beam_grnd_beam * (1.0 - ts.band[ib].albedo_grnd_beam)
        diff_bal = albedo_diff + frac_abs_can_diff +
                   frac_diff_grnd_diff * (1.0 - ts.band[ib].albedo_grnd_diff)
        @test beam_bal ≈ 1.0 atol = 1.0e-5
        @test diff_bal ≈ 1.0 atol = 1.0e-5

        # Absolute boundary type leaves the actual BCs on the band record
        @test ts.band[ib].Rbeam_atm ≈ Rbeam_atm
        @test ts.band[ib].Rdiff_atm ≈ Rdiff_atm
    end

    @testset "two-layer canopy solve closes balance" begin
        setup_rad_params!()

        ncan = 2
        ncol = 2
        ts = CLM.twostream_type()
        CLM.AllocInitTwoStream!(ts, [1, 2], ncan, ncol)
        ts.n_col[1] = 2
        ts.n_col[2] = 2
        ts.force_prep = true

        # Upper layer: vegetated + air
        u_veg = ts.scelg[1, 1]; u_veg.pft = 1; u_veg.area = 0.6; u_veg.lai = 1.5; u_veg.sai = 0.3
        u_air = ts.scelg[1, 2]; u_air.pft = CLM.air_ft; u_air.area = 0.4; u_air.lai = 0.0; u_air.sai = 0.0
        # Lower layer: vegetated + air
        l_veg = ts.scelg[2, 1]; l_veg.pft = 1; l_veg.area = 0.5; l_veg.lai = 1.0; l_veg.sai = 0.2
        l_air = ts.scelg[2, 2]; l_air.pft = CLM.air_ft; l_air.area = 0.5; l_air.lai = 0.0; l_air.sai = 0.0

        CLM.GetNSCel!(ts)
        @test ts.n_scel == 4

        for b in ts.band
            b.albedo_grnd_diff = 0.15
            b.albedo_grnd_beam = 0.15
        end

        CLM.CanopyPrep!(ts, 0.0)
        CLM.ZenithPrep!(ts, 0.5)

        n_eq = 2 * ts.n_scel
        taulamb = zeros(Float64, n_eq)
        omega = zeros(Float64, n_eq, n_eq)
        ipiv = zeros(Int, n_eq)

        ib = 1

        (_, _, consv_err, _, _, _, _, _) =
            CLM.Solve!(ts, ib, CLM.normalized_upper_boundary,
                       1.0, 1.0, taulamb, omega, ipiv)

        @test consv_err < CLM.rel_err_thresh

        # Normalized boundary type resets the band BCs to the unset sentinel
        @test ts.band[ib].Rbeam_atm == CLM.twostr_unset_r8
        @test ts.band[ib].Rdiff_atm == CLM.twostr_unset_r8
    end

    @testset "snow-cover solve still conserves" begin
        setup_rad_params!()

        ts = CLM.twostream_type()
        CLM.AllocInitTwoStream!(ts, [1, 2], 1, 1)
        ts.n_col[1] = 1
        ts.force_prep = true

        veg = ts.scelg[1, 1]
        veg.pft = 1; veg.area = 1.0; veg.lai = 3.0; veg.sai = 0.5

        CLM.GetNSCel!(ts)
        for b in ts.band
            b.albedo_grnd_diff = 0.3
            b.albedo_grnd_beam = 0.3
        end

        CLM.CanopyPrep!(ts, 0.5)   # 50% snow cover
        CLM.ZenithPrep!(ts, 0.8)

        n_eq = 2 * ts.n_scel
        taulamb = zeros(Float64, n_eq)
        omega = zeros(Float64, n_eq, n_eq)
        ipiv = zeros(Int, n_eq)

        (_, _, consv_err, _, _, _, _, _) =
            CLM.Solve!(ts, 1, CLM.absolute_upper_boundary, 1.0, 1.0, taulamb, omega, ipiv)

        @test consv_err < CLM.rel_err_thresh
    end
end
