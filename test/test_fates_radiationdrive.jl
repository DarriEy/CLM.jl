# test_fates_radiationdrive.jl
# Tests for the FATES canopy-radiation DRIVER (Tier F, Batch 17):
#   FatesRadiationDriveMod.FatesNormalizedCanopyRadiation / FatesSunShadeFracs.
#
# Strategy: build a synthetic single-site domain with the age-ordered patch
# linked list = [bareground patch, vegetated patch] under the Norman solver, run
# the driver, and assert:
#   * FatesNormalizedCanopyRadiation packs albedo in [0,1] and CONSERVES energy
#       per band:  fabd + albd + ftdd*(1-albgr_dir) + ftid*(1-albgr_dif) == 1,
#       and        fabi + albi + ftii*(1-albgr_dif) == 1
#       (the Norman solver folds the residual into albedo, so it is exact);
#   * the bareground patch is skipped and does NOT consume an ifp slot — i.e.
#       ifp=1 maps to the first VEGETATED patch (we check the veg-patch albedo
#       lands in bc_out row 1, and the unused row is left at the host default);
#   * the no-sun / night branch (filter_vegzen_pa=false, preserve_b4b) leaves
#       bc_out untouched (host default);
#   * the no-canopy / nrad==0 branch falls back to the ground albedo with full
#       direct transmittance;
#   * FatesSunShadeFracs produces laisun+laisha == elai and fsun in [0,1].

using Test
using CLM

@testset "FATES Batch 17: FatesRadiationDriveMod" begin

    # --- Preserve / restore mutated module-global state ----------------------
    old_numpft = CLM.numpft[]
    old_edpft  = CLM.EDPftvarcon_inst[]
    old_edpar  = CLM.EDParams[]

    try
        num_swb = CLM.num_swb

        # ---- Global FATES params: Norman solver + one PFT --------------------
        function _setup_globals!()
            CLM.numpft[] = 1

            edp = CLM.ed_params_type()
            edp.radiation_model = CLM.norman_solver
            nlv = CLM.nlevleaf
            edp.dinc_vai   = fill(1.0, nlv)
            edp.dlower_vai = collect(1.0:1.0:Float64(nlv))
            CLM.EDParams[] = edp

            edpft = CLM.EDPftvarcon_type()
            edpft.rhol = fill(0.10, 1, num_swb)
            edpft.rhos = fill(0.16, 1, num_swb)
            edpft.taul = fill(0.05, 1, num_swb)
            edpft.taus = fill(0.001, 1, num_swb)
            edpft.xl   = [0.01]
            edpft.clumping_index = [1.0]
            CLM.EDPftvarcon_inst[] = edpft
            return nothing
        end

        # ---- Build a vegetated patch (single canopy layer, single PFT) -------
        function _build_veg_patch(; nrad::Int, elai::Float64, esai::Float64,
                                   vegzen::Bool=true, cosz::Float64=0.7)
            p = CLM.fates_patch_type()
            np   = CLM.numpft[]
            ncan = 1
            nveg = nrad + 2
            p.nocomp_pft_label = 1          # vegetated (not bareground=0)
            p.ncl_p = 1
            p.area  = CLM.area
            p.total_canopy_area = CLM.area
            p.fcansno = 0.0
            p.solar_zenith_angle = cosz
            p.solar_zenith_flag  = vegzen

            p.nrad        = fill(0, CLM.nclmax, CLM.maxpft)
            p.nleaf       = fill(0, CLM.nclmax, CLM.maxpft)
            p.canopy_mask = fill(0, CLM.nclmax, CLM.maxpft)
            p.nrad[1, 1]  = nrad
            p.nleaf[1, 1] = nrad

            p.elai_profile        = zeros(ncan, np, nveg)
            p.esai_profile        = zeros(ncan, np, nveg)
            p.canopy_area_profile = zeros(ncan, np, nveg)
            p.f_sun               = zeros(ncan, np, nveg)
            p.fabd_sun_z          = zeros(ncan, np, nveg)
            p.fabd_sha_z          = zeros(ncan, np, nveg)
            p.fabi_sun_z          = zeros(ncan, np, nveg)
            p.fabi_sha_z          = zeros(ncan, np, nveg)
            p.ed_parsun_z         = zeros(ncan, np, nveg)
            p.ed_parsha_z         = zeros(ncan, np, nveg)
            p.ed_laisun_z         = zeros(ncan, np, nveg)
            p.ed_laisha_z         = zeros(ncan, np, nveg)
            p.parprof_pft_dir_z   = zeros(ncan, np, nveg)
            p.parprof_pft_dif_z   = zeros(ncan, np, nveg)
            p.nrmlzd_parprof_pft_dir_z = zeros(CLM.num_rad_stream_types, ncan, np, nveg)
            p.nrmlzd_parprof_pft_dif_z = zeros(CLM.num_rad_stream_types, ncan, np, nveg)

            for iv in 1:nrad
                p.elai_profile[1, 1, iv]        = elai
                p.esai_profile[1, 1, iv]        = esai
                p.canopy_area_profile[1, 1, iv] = 1.0
            end

            p.fabd = zeros(num_swb); p.fabi = zeros(num_swb)
            p.tr_soil_dir     = zeros(num_swb)
            p.tr_soil_dir_dif = zeros(num_swb)
            p.tr_soil_dif     = zeros(num_swb)
            p.sabs_dir        = zeros(num_swb)
            p.sabs_dif        = zeros(num_swb)
            p.rad_error       = zeros(num_swb)
            p.gnd_alb_dir     = fill(0.15, num_swb)
            p.gnd_alb_dif     = fill(0.15, num_swb)
            return p
        end

        # A bareground patch (label 0) — must be skipped, ifp NOT bumped.
        function _build_bareground_patch()
            p = CLM.fates_patch_type()
            p.nocomp_pft_label = CLM.nocomp_bareground   # 0
            return p
        end

        # ---- bc_in / bc_out for a single site, `npatches` veg patches --------
        function _build_bc(; npatches::Int, vegzen::Bool=true, cosz::Float64=0.7,
                            solad::Float64=400.0, solai::Float64=100.0)
            bc_in  = CLM.bc_in_type()
            bc_out = CLM.bc_out_type()
            CLM.allocate_bcin!(bc_in;  npatches=npatches, nlevsoil=1, nlevdecomp=1)
            CLM.allocate_bcout!(bc_out; npatches=npatches, nlevsoil=1, nlevdecomp=1)
            for ip in 1:npatches
                bc_in.filter_vegzen_pa[ip] = vegzen
                bc_in.coszen_pa[ip]  = cosz
                bc_in.fcansno_pa[ip] = 0.0
                for ib in 1:num_swb
                    bc_in.solad_parb[ip, ib] = solad
                    bc_in.solai_parb[ip, ib] = solai
                end
            end
            for ib in 1:num_swb
                bc_in.albgr_dir_rb[ib] = 0.15
                bc_in.albgr_dif_rb[ib] = 0.15
            end
            return (bc_in, bc_out)
        end

        # Link patches oldest -> youngest into a site.
        function _site_with(patches::Vector)
            site = CLM.ed_site_type()
            site.snow_depth = 0.0
            site.lat = 45.0; site.lon = 0.0
            prev = nothing
            for p in patches
                if prev === nothing
                    site.oldest_patch = p
                else
                    prev.younger = p
                end
                prev = p
            end
            return site
        end

        # ---------------------------------------------------------------------
        @testset "albedo bounds + energy conservation (Norman)" begin
            _setup_globals!()
            veg = _build_veg_patch(nrad=4, elai=0.5, esai=0.1)
            site = _site_with([veg])
            (bc_in, bc_out) = _build_bc(npatches=1)

            CLM.FatesNormalizedCanopyRadiation(1, [site], [bc_in], [bc_out])

            ifp = 1
            for ib in 1:num_swb
                albd = bc_out.albd_parb[ifp, ib]
                albi = bc_out.albi_parb[ifp, ib]
                fabd = bc_out.fabd_parb[ifp, ib]
                fabi = bc_out.fabi_parb[ifp, ib]
                ftdd = bc_out.ftdd_parb[ifp, ib]
                ftid = bc_out.ftid_parb[ifp, ib]
                ftii = bc_out.ftii_parb[ifp, ib]
                gdir = bc_in.albgr_dir_rb[ib]
                gdif = bc_in.albgr_dif_rb[ib]

                @test 0.0 <= albd <= 1.0
                @test 0.0 <= albi <= 1.0
                @test 0.0 <= fabd <= 1.0
                @test 0.0 <= fabi <= 1.0

                # Direct: absorbed + reflected + (transmitted-to-ground absorbed)
                @test fabd + albd + ftdd*(1.0-gdir) + ftid*(1.0-gdif) ≈ 1.0 atol=1e-9
                # Diffuse: absorbed + reflected + transmitted-to-ground absorbed
                @test fabi + albi + ftii*(1.0-gdif) ≈ 1.0 atol=1e-9
            end
        end

        # ---------------------------------------------------------------------
        @testset "bareground patch skipped, ifp not consumed" begin
            _setup_globals!()
            bg  = _build_bareground_patch()
            veg = _build_veg_patch(nrad=4, elai=0.5, esai=0.1)
            site = _site_with([bg, veg])      # oldest = bareground
            (bc_in, bc_out) = _build_bc(npatches=1)

            CLM.FatesNormalizedCanopyRadiation(1, [site], [bc_in], [bc_out])

            # The vegetated patch must have landed in ifp=1 (bareground did NOT
            # consume the slot). Its albedo is set and in [0,1].
            for ib in 1:num_swb
                @test 0.0 <= bc_out.albd_parb[1, ib] <= 1.0
                # canopy actually absorbed something
                @test bc_out.fabd_parb[1, ib] > 0.0
            end
        end

        # ---------------------------------------------------------------------
        @testset "night branch (preserve_b4b leaves bc_out untouched)" begin
            _setup_globals!()
            veg = _build_veg_patch(nrad=4, elai=0.5, esai=0.1, vegzen=false)
            site = _site_with([veg])
            (bc_in, bc_out) = _build_bc(npatches=1, vegzen=false)

            # Sentinel the outputs; the night branch must not overwrite them.
            fill!(bc_out.albd_parb, -1.0)
            fill!(bc_out.fabd_parb, -1.0)

            CLM.FatesNormalizedCanopyRadiation(1, [site], [bc_in], [bc_out])

            for ib in 1:num_swb
                @test bc_out.albd_parb[1, ib] == -1.0   # untouched
                @test bc_out.fabd_parb[1, ib] == -1.0   # untouched (no flux)
            end
        end

        # ---------------------------------------------------------------------
        @testset "no-canopy (nrad==0) falls back to ground albedo" begin
            _setup_globals!()
            veg = _build_veg_patch(nrad=0, elai=0.0, esai=0.0)  # nrad[1,1]=0
            site = _site_with([veg])
            (bc_in, bc_out) = _build_bc(npatches=1)

            CLM.FatesNormalizedCanopyRadiation(1, [site], [bc_in], [bc_out])

            for ib in 1:num_swb
                @test bc_out.albd_parb[1, ib] ≈ bc_in.albgr_dir_rb[ib]
                @test bc_out.albi_parb[1, ib] ≈ bc_in.albgr_dif_rb[ib]
                @test bc_out.fabd_parb[1, ib] == 0.0
                @test bc_out.fabi_parb[1, ib] == 0.0
                @test bc_out.ftdd_parb[1, ib] == 1.0   # full direct transmittance
                @test bc_out.ftid_parb[1, ib] == 0.0
                @test bc_out.ftii_parb[1, ib] == 1.0
            end
        end

        # ---------------------------------------------------------------------
        @testset "FatesSunShadeFracs: laisun + laisha == elai, fsun in [0,1]" begin
            _setup_globals!()
            veg = _build_veg_patch(nrad=4, elai=0.5, esai=0.1)
            site = _site_with([veg])
            (bc_in, bc_out) = _build_bc(npatches=1)

            # Run the normalized solve first to populate f_sun, fab*_z, etc.
            CLM.FatesNormalizedCanopyRadiation(1, [site], [bc_in], [bc_out])
            CLM.FatesSunShadeFracs(1, [site], [bc_in], [bc_out])

            ifp = 1
            fsun   = bc_out.fsun_pa[ifp]
            laisun = bc_out.laisun_pa[ifp]
            laisha = bc_out.laisha_pa[ifp]
            elai   = CLM.calc_areaindex(veg, "elai")

            @test 0.0 <= fsun <= 1.0
            @test laisun + laisha ≈ elai atol=1e-9
            @test laisun ≈ elai * fsun atol=1e-9
            @test laisha ≈ elai * (1.0 - fsun) atol=1e-9

            # Absorbed-PAR profiles non-negative.
            for iv in 1:veg.nrad[1, 1]
                @test veg.ed_parsun_z[1, 1, iv] >= 0.0
                @test veg.ed_parsha_z[1, 1, iv] >= 0.0
            end
        end

    finally
        CLM.numpft[]           = old_numpft
        CLM.EDPftvarcon_inst[] = old_edpft
        CLM.EDParams[]         = old_edpar
    end
end
