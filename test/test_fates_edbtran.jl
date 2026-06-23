# test_fates_edbtran.jl
# Tests for FATES Batch 12 (Tier F): EDBtranMod — the FATES soil-moisture
# limited transpiration wetness factor (btran) and root-soil water uptake
# distribution.
#
# Strategy:
#   * check_layer_water: exercise the two thresholds (positive liquid water AND
#     temperature above soil_tfrz_thresh + freezing point) and the boundary cases.
#   * get_active_suction_layers!: assert per-layer activation follows
#     check_layer_water, and that filter_btran == false zeroes the whole column.
#   * btran_ed!: build a synthetic single-site / single-veg-patch / single-cohort
#     system with a known root profile (1-parameter exponential), known soil
#     suction, porosity and PFT smpsc/smpso, and:
#       - assert btran_ft matches a hand-computed Fortran-equivalent value,
#       - assert rootr_pasl sums to 1 (conservation / re-normalization),
#       - assert btran_pa equals the (single-PFT) btran_ft when hydraulics off,
#       - assert frozen / dry layers contribute nothing,
#       - assert a wetter (less-negative smp) column yields a larger btran
#         (monotonicity), and a bareground patch is skipped (ifp not advanced).

using Test
using CLM

@testset "FATES Batch 12: EDBtranMod" begin

    # --- Preserve and restore mutated module-global state --------------------
    old_numpft   = CLM.numpft[]
    old_hydro    = CLM.hlm_use_planthydro[]
    old_edpft    = CLM.EDPftvarcon_inst[]

    try
        sftz = CLM.soil_tfrz_thresh + CLM.t_water_freeze_k_1atm  # active-water temp threshold [K]

        # Helper: allocate prt_params for a single-PFT exponential root profile.
        # (set_root_fraction reads prt_params.fnrt_prof_mode/a/b per PFT.)
        function _setup_prt_params!(npft)
            CLM.allocate_prt_params!(CLM.prt_params, npft, CLM.num_organ_types, 1)
            CLM.prt_params.woody          .= CLM.itrue
            CLM.prt_params.fnrt_prof_mode .= 2.0   # 1-parameter exponential
            CLM.prt_params.fnrt_prof_a    .= 4.0
            CLM.prt_params.fnrt_prof_b    .= 0.0
            return nothing
        end

        # ---------------------------------------------------------------------
        # check_layer_water
        # ---------------------------------------------------------------------
        @testset "check_layer_water" begin
            # warm + wet -> available
            @test CLM.check_layer_water(0.3, sftz + 1.0) == true
            # warm but no liquid water -> not available
            @test CLM.check_layer_water(0.0, sftz + 1.0) == false
            @test CLM.check_layer_water(-0.1, sftz + 1.0) == false
            # wet but at/below the freeze threshold -> not available (strict >)
            @test CLM.check_layer_water(0.3, sftz) == false
            @test CLM.check_layer_water(0.3, sftz - 1.0) == false
            # just above threshold -> available
            @test CLM.check_layer_water(1.0e-9, sftz + 1.0e-6) == true
        end

        # ---------------------------------------------------------------------
        # get_active_suction_layers!
        # ---------------------------------------------------------------------
        @testset "get_active_suction_layers!" begin
            nlevsoil = 4
            bc_in  = CLM.bc_in_type()
            bc_out = CLM.bc_out_type()
            bc_in.nlevsoil      = nlevsoil
            bc_in.h2o_liqvol_sl = [0.3, 0.0, 0.3, 0.3]            # layer 2 dry
            bc_in.tempk_sl      = [sftz + 5, sftz + 5, sftz - 5, sftz + 5]  # layer 3 frozen
            bc_out.active_suction_sl = falses(nlevsoil)

            bc_in.filter_btran = true
            CLM.get_active_suction_layers!([nothing], [bc_in], [bc_out])
            @test bc_out.active_suction_sl == [true, false, false, true]

            # filter off -> all inactive regardless of water state
            bc_in.filter_btran = false
            CLM.get_active_suction_layers!([nothing], [bc_in], [bc_out])
            @test all(.!bc_out.active_suction_sl)
        end

        # ---------------------------------------------------------------------
        # btran_ed! — synthetic single site / single veg patch / single cohort
        # ---------------------------------------------------------------------
        @testset "btran_ed! single pft" begin
            CLM.numpft[] = 1
            CLM.hlm_use_planthydro[] = CLM.ifalse

            nlevsoil = 3

            # PFT params: smpso (open), smpsc (close), both negative [mm].
            edpft = CLM.EDPftvarcon_type()
            edpft.smpso = [-66000.0]
            edpft.smpsc = [-255000.0]
            CLM.EDPftvarcon_inst[] = edpft

            # Root profile: 1-parameter exponential (mode 2), scale a.
            _setup_prt_params!(1)

            # Site soil layering.
            site = CLM.ed_site_type()
            site.nlevsoil     = nlevsoil
            site.zi_soil      = [0.0, 0.1, 0.3, 0.6]  # zero index for surface
            site.rootfrac_scr = zeros(Float64, nlevsoil)

            # Cohort with a known leaf-area-weighted conductance.
            coh = CLM.fates_cohort_type()
            coh.pft           = 1
            coh.g_sb_laweight = 0.5
            coh.taller        = nothing
            coh.shorter       = nothing

            # One veg patch holding that cohort.
            patch = CLM.fates_patch_type()
            patch.nocomp_pft_label = 1
            patch.area    = CLM.area
            patch.btran_ft = zeros(Float64, CLM.maxpft)
            patch.tallest  = coh
            patch.shortest = coh
            patch.older    = nothing
            patch.younger  = nothing

            site.oldest_patch   = patch
            site.youngest_patch = patch

            # bc_in: all layers warm + wet, known suction / porosity.
            bc_in = CLM.bc_in_type()
            bc_in.nlevsoil                    = nlevsoil
            bc_in.max_rooting_depth_index_col = nlevsoil
            bc_in.h2o_liqvol_sl = fill(0.3, nlevsoil)
            bc_in.tempk_sl      = fill(sftz + 5.0, nlevsoil)
            bc_in.smp_sl        = [-100000.0, -150000.0, -200000.0]
            bc_in.eff_porosity_sl = fill(0.4, nlevsoil)
            bc_in.watsat_sl       = fill(0.5, nlevsoil)

            bc_out = CLM.bc_out_type()
            bc_out.rootr_pasl = zeros(Float64, 1, nlevsoil)
            bc_out.btran_pa   = zeros(Float64, 1)

            CLM.btran_ed!([site], [bc_in], [bc_out])

            # --- Hand-computed Fortran-equivalent btran_ft -------------------
            # Recompute the root profile the same way the module does.
            rootfr = zeros(Float64, nlevsoil)
            CLM.set_root_fraction(rootfr, 1, site.zi_soil; max_nlevroot=nlevsoil)

            smpsc = edpft.smpsc[1]
            smpso = edpft.smpso[1]
            poro_ratio = 0.4 / 0.5
            expect_btran = 0.0
            for j in 1:nlevsoil
                smp_node = max(smpsc, bc_in.smp_sl[j])
                rresis = min(poro_ratio * (smp_node - smpsc) / (smpso - smpsc), 1.0)
                expect_btran += rootfr[j] * rresis
            end

            @test isapprox(patch.btran_ft[1], expect_btran; rtol=1e-12)

            # Single PFT -> btran_pa == btran_ft.
            @test isapprox(bc_out.btran_pa[1], patch.btran_ft[1]; rtol=1e-12)
            @test 0.0 < bc_out.btran_pa[1] <= 1.0

            # rootr_pasl is a distribution -> sums to 1 (renormalized).
            @test isapprox(sum(bc_out.rootr_pasl[1, :]), 1.0; atol=1e-12)
            @test all(bc_out.rootr_pasl[1, :] .>= 0.0)
        end

        # ---------------------------------------------------------------------
        # frozen / dry layers contribute nothing; bareground patch is skipped
        # ---------------------------------------------------------------------
        @testset "btran_ed! frozen layers + bareground skip" begin
            CLM.numpft[] = 1
            CLM.hlm_use_planthydro[] = CLM.ifalse
            nlevsoil = 3

            edpft = CLM.EDPftvarcon_type()
            edpft.smpso = [-66000.0]
            edpft.smpsc = [-255000.0]
            CLM.EDPftvarcon_inst[] = edpft

            _setup_prt_params!(1)

            site = CLM.ed_site_type()
            site.nlevsoil     = nlevsoil
            site.zi_soil      = [0.0, 0.1, 0.3, 0.6]
            site.rootfrac_scr = zeros(Float64, nlevsoil)

            coh = CLM.fates_cohort_type()
            coh.pft = 1; coh.g_sb_laweight = 0.5
            coh.taller = nothing; coh.shorter = nothing

            # A leading bareground patch (must be skipped) then the veg patch.
            bare = CLM.fates_patch_type()
            bare.nocomp_pft_label = CLM.nocomp_bareground
            bare.area = CLM.area
            bare.btran_ft = zeros(Float64, CLM.maxpft)
            bare.tallest = nothing; bare.shortest = nothing

            patch = CLM.fates_patch_type()
            patch.nocomp_pft_label = 1
            patch.area = CLM.area
            patch.btran_ft = zeros(Float64, CLM.maxpft)
            patch.tallest = coh; patch.shortest = coh

            bare.older = nothing; bare.younger = patch
            patch.older = bare;   patch.younger = nothing
            site.oldest_patch = bare; site.youngest_patch = patch

            # All warm/wet EXCEPT layer 2 frozen and layer 3 dry.
            bc_in = CLM.bc_in_type()
            bc_in.nlevsoil = nlevsoil
            bc_in.max_rooting_depth_index_col = nlevsoil
            bc_in.h2o_liqvol_sl = [0.3, 0.3, 0.0]                       # layer 3 dry
            bc_in.tempk_sl      = [sftz + 5, sftz - 5, sftz + 5]        # layer 2 frozen
            bc_in.smp_sl        = [-100000.0, -150000.0, -200000.0]
            bc_in.eff_porosity_sl = fill(0.4, nlevsoil)
            bc_in.watsat_sl       = fill(0.5, nlevsoil)

            # Only one veg patch -> bc_out sized for one patch index.
            bc_out = CLM.bc_out_type()
            bc_out.rootr_pasl = zeros(Float64, 1, nlevsoil)
            bc_out.btran_pa   = zeros(Float64, 1)

            CLM.btran_ed!([site], [bc_in], [bc_out])

            # Only layer 1 is active -> all uptake in layer 1, sum is 1.
            @test isapprox(bc_out.rootr_pasl[1, 1], 1.0; atol=1e-12)
            @test isapprox(bc_out.rootr_pasl[1, 2], 0.0; atol=1e-12)
            @test isapprox(bc_out.rootr_pasl[1, 3], 0.0; atol=1e-12)

            # btran_ft only accumulated the single active layer (rootfr[1]*rresis).
            rootfr = zeros(Float64, nlevsoil)
            CLM.set_root_fraction(rootfr, 1, site.zi_soil; max_nlevroot=nlevsoil)
            smpsc = edpft.smpsc[1]; smpso = edpft.smpso[1]
            poro_ratio = 0.4 / 0.5
            smp_node = max(smpsc, bc_in.smp_sl[1])
            rresis = min(poro_ratio * (smp_node - smpsc) / (smpso - smpsc), 1.0)
            @test isapprox(patch.btran_ft[1], rootfr[1] * rresis; rtol=1e-12)
        end

        # ---------------------------------------------------------------------
        # monotonicity: wetter soil (less-negative smp) -> larger btran
        # ---------------------------------------------------------------------
        @testset "btran_ed! wetness monotonicity" begin
            CLM.numpft[] = 1
            CLM.hlm_use_planthydro[] = CLM.ifalse
            nlevsoil = 3

            edpft = CLM.EDPftvarcon_type()
            edpft.smpso = [-66000.0]
            edpft.smpsc = [-255000.0]
            CLM.EDPftvarcon_inst[] = edpft

            _setup_prt_params!(1)

            function _run_btran(smp_vec)
                site = CLM.ed_site_type()
                site.nlevsoil = nlevsoil
                site.zi_soil  = [0.0, 0.1, 0.3, 0.6]
                site.rootfrac_scr = zeros(Float64, nlevsoil)

                coh = CLM.fates_cohort_type()
                coh.pft = 1; coh.g_sb_laweight = 0.5
                coh.taller = nothing; coh.shorter = nothing

                patch = CLM.fates_patch_type()
                patch.nocomp_pft_label = 1
                patch.area = CLM.area
                patch.btran_ft = zeros(Float64, CLM.maxpft)
                patch.tallest = coh; patch.shortest = coh
                patch.older = nothing; patch.younger = nothing
                site.oldest_patch = patch; site.youngest_patch = patch

                bc_in = CLM.bc_in_type()
                bc_in.nlevsoil = nlevsoil
                bc_in.max_rooting_depth_index_col = nlevsoil
                bc_in.h2o_liqvol_sl = fill(0.3, nlevsoil)
                bc_in.tempk_sl      = fill(sftz + 5.0, nlevsoil)
                bc_in.smp_sl        = copy(smp_vec)
                bc_in.eff_porosity_sl = fill(0.4, nlevsoil)
                bc_in.watsat_sl       = fill(0.5, nlevsoil)

                bc_out = CLM.bc_out_type()
                bc_out.rootr_pasl = zeros(Float64, 1, nlevsoil)
                bc_out.btran_pa   = zeros(Float64, 1)

                CLM.btran_ed!([site], [bc_in], [bc_out])
                return bc_out.btran_pa[1]
            end

            dry_btran = _run_btran([-200000.0, -200000.0, -200000.0])
            wet_btran = _run_btran([-100000.0, -100000.0, -100000.0])
            @test wet_btran > dry_btran
            @test 0.0 < dry_btran <= 1.0
            @test 0.0 < wet_btran <= 1.0
        end

    finally
        # --- restore module globals ------------------------------------------
        CLM.numpft[]             = old_numpft
        CLM.hlm_use_planthydro[] = old_hydro
        CLM.EDPftvarcon_inst[]   = old_edpft
    end
end
