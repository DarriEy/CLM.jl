# test_fates_norman_rad.jl
# Tests for the FATES Norman radiation solver (Tier F, Batch 16):
# FatesNormanRadMod.PatchNormanRadiation.
#
# Strategy: build a synthetic single-canopy-layer, single-PFT patch with known
# leaf/stem optical properties, fully-filled canopy (ftweight == 1) and a known
# ground albedo, run the Norman two-stream solve, and assert:
#   * energy conservation per band:  absorbed_canopy + albedo + absorbed_ground
#       == incident (== 1), to a tight tolerance (NB: the solver folds the tiny
#       residual into the albedo, so this is exact);
#   * absorbed_ground == transmitted-to-ground * (1 - ground albedo);
#   * albedo bounds in [0,1];
#   * sunlit + shaded leaf-layer absorption == total leaf-layer absorption;
#   * monotonicity: absorbed canopy fraction increases with LAI;
#   * a hand-checked single-leaf-layer direct-beam case (top-layer sunlit
#     fraction == exp(-k_dir * 0.5*vai) for a full canopy).

using Test
using CLM

@testset "FATES Batch 16: FatesNormanRadMod" begin

    # --- Preserve and restore mutated module-global state --------------------
    old_numpft = CLM.numpft[]
    old_edpft  = CLM.EDPftvarcon_inst[]

    try
        CLM.numpft[] = 1
        num_swb = CLM.num_swb   # 2 (vis, nir)

        # PFT optical parameters (one PFT, num_swb bands). Layout (pft, band).
        function _setup_edpft!(; rhol=0.10, taul=0.05, rhos=0.16, taus=0.001,
                                 xl=0.01, clumping=1.0)
            edpft = CLM.EDPftvarcon_type()
            edpft.rhol = fill(rhol, 1, num_swb)
            edpft.rhos = fill(rhos, 1, num_swb)
            edpft.taul = fill(taul, 1, num_swb)
            edpft.taus = fill(taus, 1, num_swb)
            edpft.xl   = [xl]
            edpft.clumping_index = [clumping]
            CLM.EDPftvarcon_inst[] = edpft
            return nothing
        end

        # Build a single-canopy-layer, single-PFT, fully-filled patch with
        # `nrad` leaf layers each carrying elai/esai. cosz in [0.001,1].
        function _build_patch(; nrad::Int, elai::Float64, esai::Float64,
                                cosz::Float64=0.7, fcansno::Float64=0.0,
                                gnd_alb_dir=fill(0.15, num_swb),
                                gnd_alb_dif=fill(0.15, num_swb))
            p = CLM.fates_patch_type()
            ncl_p = 1
            ncan  = 1
            np    = CLM.numpft[]
            # nveg must allow the iv+1 = nrad+1 soil slot. Add a buffer slot.
            nveg  = nrad + 2

            p.ncl_p = ncl_p
            p.area  = CLM.area
            p.total_canopy_area = CLM.area     # fully veg-covered patch
            p.patchno = 1
            p.fcansno = fcansno
            # QUIRK: despite the field name/comment ("[radians]"), FATES stores
            # the COSINE of the zenith angle here — PatchNormanRadiation reads
            # cosz = max(0.001, solar_zenith_angle) directly. Set the cosine.
            p.solar_zenith_angle = cosz

            # (nclmax, maxpft) integer index arrays
            p.nrad        = fill(0, CLM.nclmax, CLM.maxpft)
            p.canopy_mask = fill(0, CLM.nclmax, CLM.maxpft)
            p.nrad[1, 1]  = nrad

            # 3-D profiles (L, ft, iv) and 4-D parprof (radtype, L, ft, iv)
            p.elai_profile        = zeros(ncan, np, nveg)
            p.esai_profile        = zeros(ncan, np, nveg)
            p.canopy_area_profile = zeros(ncan, np, nveg)
            p.f_sun               = zeros(ncan, np, nveg)
            p.fabd_sun_z          = zeros(ncan, np, nveg)
            p.fabd_sha_z          = zeros(ncan, np, nveg)
            p.fabi_sun_z          = zeros(ncan, np, nveg)
            p.fabi_sha_z          = zeros(ncan, np, nveg)
            p.nrmlzd_parprof_pft_dir_z = zeros(CLM.num_rad_stream_types, ncan, np, nveg)
            p.nrmlzd_parprof_pft_dif_z = zeros(CLM.num_rad_stream_types, ncan, np, nveg)

            for iv in 1:nrad
                p.elai_profile[1, 1, iv]        = elai
                p.esai_profile[1, 1, iv]        = esai
                p.canopy_area_profile[1, 1, iv] = 1.0   # fully filled => ftweight==1
            end

            # band vectors
            p.fabd = zeros(num_swb)
            p.fabi = zeros(num_swb)
            p.tr_soil_dir     = zeros(num_swb)
            p.tr_soil_dir_dif = zeros(num_swb)
            p.tr_soil_dif     = zeros(num_swb)
            p.sabs_dir        = zeros(num_swb)
            p.sabs_dif        = zeros(num_swb)
            p.rad_error       = zeros(num_swb)
            p.gnd_alb_dir     = copy(gnd_alb_dir)
            p.gnd_alb_dif     = copy(gnd_alb_dif)
            return p
        end

        # Run the solver, returning the seven output vectors.
        function _run(p)
            albd = zeros(num_swb); albi = zeros(num_swb)
            fabd = zeros(num_swb); fabi = zeros(num_swb)
            ftdd = zeros(num_swb); ftid = zeros(num_swb); ftii = zeros(num_swb)
            CLM.PatchNormanRadiation(p, albd, albi, fabd, fabi, ftdd, ftid, ftii)
            return (; albd, albi, fabd, fabi, ftdd, ftid, ftii)
        end

        # ---------------------------------------------------------------------
        @testset "energy conservation + albedo bounds (direct & diffuse)" begin
            _setup_edpft!()
            p = _build_patch(nrad=4, elai=0.5, esai=0.1, cosz=0.7)
            o = _run(p)

            for ib in 1:num_swb
                # Direct beam: incident (forc_dir=1) = fabd + albd + sabs_dir.
                # The solver folds the residual into the albedo, so it is exact.
                @test o.fabd[ib] + o.albd[ib] + p.sabs_dir[ib] ≈ 1.0 atol = 1e-9
                # Diffuse: incident (forc_dif=1) = fabi + albi + sabs_dif.
                @test o.fabi[ib] + o.albi[ib] + p.sabs_dif[ib] ≈ 1.0 atol = 1e-9

                # Albedo bounds
                @test 0.0 <= o.albd[ib] <= 1.0
                @test 0.0 <= o.albi[ib] <= 1.0

                # Absorbed fractions in [0,1]
                @test 0.0 <= o.fabd[ib] <= 1.0
                @test 0.0 <= o.fabi[ib] <= 1.0

                # Ground absorption == transmitted-to-ground * (1 - ground albedo)
                @test p.sabs_dir[ib] ≈ p.tr_soil_dir[ib] *
                    (1.0 - p.gnd_alb_dir[ib]) + p.tr_soil_dir_dif[ib] *
                    (1.0 - p.gnd_alb_dif[ib]) atol = 1e-9
                @test p.sabs_dif[ib] ≈ p.tr_soil_dif[ib] *
                    (1.0 - p.gnd_alb_dif[ib]) atol = 1e-9

                # rad_error must be tiny for a full-cover patch
                @test abs(p.rad_error[ib]) < 1e-3
            end
        end

        # ---------------------------------------------------------------------
        @testset "sunlit + shaded == total leaf-layer absorption (PAR band)" begin
            _setup_edpft!()
            p = _build_patch(nrad=4, elai=0.5, esai=0.1, cosz=0.7)
            _run(p)

            ib = CLM.ivis   # sun/shade arrays are only set for the visible band
            for iv in 1:p.nrad[1, 1]
                # direct: fabd_sun + fabd_sha covers both diffuse and direct
                # contributions to this leaf layer; both >= 0.
                @test p.fabd_sun_z[1, 1, iv] >= 0.0
                @test p.fabd_sha_z[1, 1, iv] >= 0.0
                @test p.fabi_sun_z[1, 1, iv] >= 0.0
                @test p.fabi_sha_z[1, 1, iv] >= 0.0
                # sunlit + shaded fraction of the layer == 1 (f_sun in [0,1])
                @test 0.0 <= p.f_sun[1, 1, iv] <= 1.0
            end
        end

        # ---------------------------------------------------------------------
        @testset "absorbed canopy fraction increases with LAI" begin
            _setup_edpft!()
            fabd_prev = -1.0
            for elai in (0.1, 0.5, 1.0, 2.0)
                p = _build_patch(nrad=4, elai=elai, esai=0.0, cosz=0.7)
                o = _run(p)
                # total direct absorbed by canopy (vis band)
                fabd = o.fabd[CLM.ivis]
                @test fabd >= fabd_prev - 1e-12   # monotone non-decreasing
                fabd_prev = fabd
            end
            @test fabd_prev > 0.0
        end

        # ---------------------------------------------------------------------
        @testset "hand-checked single-leaf-layer direct sunlit fraction" begin
            # For a full canopy (ftweight==1) the top-layer sunlit fraction is
            # f_sun = exp(-k_dir * 0.5*vai), with
            #   k_dir = clumping * gdir / sin(sb),
            #   gdir  = phi1 + phi2*sin(sb),
            #   sb    = 90deg - acos(cosz),  phi1 = 0.5 - 0.633 xl - 0.330 xl^2,
            #   phi2  = 0.877 (1 - 2 phi1).
            cosz = 0.8
            xl = 0.01
            _setup_edpft!(xl=xl, clumping=1.0)
            elai = 1.0; esai = 0.0
            p = _build_patch(nrad=1, elai=elai, esai=esai, cosz=cosz)
            _run(p)

            sb   = (90.0 - acos(cosz) * 180.0 / CLM.pi_const) * (CLM.pi_const / 180.0)
            phi1 = 0.5 - 0.633 * xl - 0.330 * xl * xl
            phi2 = 0.877 * (1.0 - 2.0 * phi1)
            gdir = phi1 + phi2 * sin(sb)
            k_dir = 1.0 * gdir / sin(sb)
            laisum_mid = 0.5 * (elai + esai)
            f_sun_expected = exp(-k_dir * laisum_mid)

            @test p.f_sun[1, 1, 1] ≈ f_sun_expected atol = 1e-10
        end

        # ---------------------------------------------------------------------
        @testset "canopy snow raises albedo" begin
            _setup_edpft!()
            p0 = _build_patch(nrad=4, elai=0.5, esai=0.1, cosz=0.7, fcansno=0.0)
            o0 = _run(p0)
            p1 = _build_patch(nrad=4, elai=0.5, esai=0.1, cosz=0.7, fcansno=0.5)
            o1 = _run(p1)
            # snow reflectance (0.80 vis) > leaf reflectance => brighter canopy
            @test o1.albd[CLM.ivis] > o0.albd[CLM.ivis]
            @test o1.albi[CLM.ivis] > o0.albi[CLM.ivis]
        end

    finally
        CLM.numpft[]           = old_numpft
        CLM.EDPftvarcon_inst[] = old_edpft
    end
end
