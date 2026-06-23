# test_fates_bstress.jl
# Tests for FATES Batch 13 (Tier F): FatesBstressMod — the FATES salinity
# transpiration-stress factor (bstress_sal_ft), the close sibling of the
# soil-moisture btran in EDBtranMod.
#
# Strategy (mirrors test_fates_edbtran.jl):
#   * btran_sal_stress_fates!: build a synthetic single-site / single-patch /
#     single-PFT system with a known root profile (1-parameter exponential) and
#     known per-layer salinity, then:
#       - assert bstress_sal_ft matches a hand-computed Fortran-equivalent value
#         (root-fraction-weighted sigmoid),
#       - assert frozen / dry layers contribute nothing,
#       - assert salinity monotonicity (lower salinity -> larger bstress_sal_ft,
#         because the sigmoid rresis decreases with salinity),
#       - assert the per-PFT reset to 0 happens, and that a multi-patch site
#         processes every patch (no bareground skip — Fortran has no guard).

using Test
using CLM

# Reproduce the Fortran salinity-stress sigmoid for hand-computed expectations.
_sal_rresis(sal) = min(1.244 / (1 + exp((0.186 - sal) / (-0.132))), 1.0)

@testset "FATES Batch 13: FatesBstressMod" begin

    # --- Preserve and restore mutated module-global state --------------------
    old_numpft = CLM.numpft[]
    old_edpft  = CLM.EDPftvarcon_inst[]

    try
        sftz = CLM.soil_tfrz_thresh + CLM.t_water_freeze_k_1atm  # active-water temp threshold [K]

        # Helper: allocate prt_params for a single-PFT exponential root profile.
        function _setup_prt_params!(npft)
            CLM.allocate_prt_params!(CLM.prt_params, npft, CLM.num_organ_types, 1)
            CLM.prt_params.woody          .= CLM.itrue
            CLM.prt_params.fnrt_prof_mode .= 2.0   # 1-parameter exponential
            CLM.prt_params.fnrt_prof_a    .= 4.0
            CLM.prt_params.fnrt_prof_b    .= 0.0
            return nothing
        end

        # Build a single-site / single-patch / single-PFT system with a given
        # per-layer salinity / liquid-water / temperature column.
        function _build(nlevsoil; salinity, h2o, tempk, npatch=1)
            site = CLM.ed_site_type()
            site.nlevsoil     = nlevsoil
            site.zi_soil      = [0.0, 0.1, 0.3, 0.6][1:nlevsoil+1]
            site.rootfrac_scr = zeros(Float64, nlevsoil)

            # Link npatch identical veg patches oldest->youngest.
            patches = CLM.fates_patch_type[]
            for _ in 1:npatch
                p = CLM.fates_patch_type()
                p.nocomp_pft_label = 1
                p.area     = CLM.area
                p.btran_ft = zeros(Float64, CLM.maxpft)
                # bstress_sal_ft starts as NaN (default); the routine must zero it.
                p.bstress_sal_ft = fill(NaN, CLM.maxpft)
                push!(patches, p)
            end
            for i in 1:npatch
                patches[i].older   = i > 1       ? patches[i-1] : nothing
                patches[i].younger = i < npatch  ? patches[i+1] : nothing
            end
            site.oldest_patch   = patches[1]
            site.youngest_patch = patches[end]

            bc_in = CLM.bc_in_type()
            bc_in.nlevsoil                    = nlevsoil
            bc_in.max_rooting_depth_index_col = nlevsoil
            bc_in.salinity_sl   = copy(salinity)
            bc_in.h2o_liqvol_sl = copy(h2o)
            bc_in.tempk_sl      = copy(tempk)

            return site, bc_in, patches
        end

        # ---------------------------------------------------------------------
        # core: hand-computed root-weighted salinity stress
        # ---------------------------------------------------------------------
        @testset "btran_sal_stress_fates! single pft" begin
            CLM.numpft[] = 1
            nlevsoil = 3

            # PFT params not actually read by this routine, but EDPftvarcon must
            # be populated to keep set_root_fraction / module globals consistent.
            edpft = CLM.EDPftvarcon_type()
            edpft.smpso = [-66000.0]
            edpft.smpsc = [-255000.0]
            CLM.EDPftvarcon_inst[] = edpft
            _setup_prt_params!(1)

            salinity = [0.1, 0.2, 0.5]
            site, bc_in, patches = _build(nlevsoil;
                salinity = salinity,
                h2o      = fill(0.3, nlevsoil),
                tempk    = fill(sftz + 5.0, nlevsoil))

            CLM.btran_sal_stress_fates!([site], [bc_in])

            # Hand-compute the root-fraction-weighted sigmoid.
            rootfr = zeros(Float64, nlevsoil)
            CLM.set_root_fraction(rootfr, 1, site.zi_soil; max_nlevroot=nlevsoil)
            expect = sum(rootfr[j] * _sal_rresis(salinity[j]) for j in 1:nlevsoil)

            @test isapprox(patches[1].bstress_sal_ft[1], expect; rtol=1e-12)
            @test patches[1].bstress_sal_ft[1] > 0.0
            # not NaN -> the per-PFT reset to 0 fired before accumulation.
            @test !isnan(patches[1].bstress_sal_ft[1])
        end

        # ---------------------------------------------------------------------
        # frozen / dry layers contribute nothing
        # ---------------------------------------------------------------------
        @testset "btran_sal_stress_fates! frozen/dry layers" begin
            CLM.numpft[] = 1
            nlevsoil = 3

            edpft = CLM.EDPftvarcon_type()
            edpft.smpso = [-66000.0]; edpft.smpsc = [-255000.0]
            CLM.EDPftvarcon_inst[] = edpft
            _setup_prt_params!(1)

            salinity = [0.1, 0.2, 0.5]
            # layer 2 frozen, layer 3 dry -> only layer 1 contributes.
            site, bc_in, patches = _build(nlevsoil;
                salinity = salinity,
                h2o      = [0.3, 0.3, 0.0],
                tempk    = [sftz + 5, sftz - 5, sftz + 5])

            CLM.btran_sal_stress_fates!([site], [bc_in])

            rootfr = zeros(Float64, nlevsoil)
            CLM.set_root_fraction(rootfr, 1, site.zi_soil; max_nlevroot=nlevsoil)
            expect = rootfr[1] * _sal_rresis(salinity[1])  # only layer 1 active

            @test isapprox(patches[1].bstress_sal_ft[1], expect; rtol=1e-12)
        end

        # ---------------------------------------------------------------------
        # salinity monotonicity: lower salinity -> larger bstress_sal_ft
        # ---------------------------------------------------------------------
        @testset "btran_sal_stress_fates! salinity monotonicity" begin
            CLM.numpft[] = 1
            nlevsoil = 3

            edpft = CLM.EDPftvarcon_type()
            edpft.smpso = [-66000.0]; edpft.smpsc = [-255000.0]
            CLM.EDPftvarcon_inst[] = edpft
            _setup_prt_params!(1)

            function _run(sal)
                site, bc_in, patches = _build(nlevsoil;
                    salinity = fill(sal, nlevsoil),
                    h2o      = fill(0.3, nlevsoil),
                    tempk    = fill(sftz + 5.0, nlevsoil))
                CLM.btran_sal_stress_fates!([site], [bc_in])
                return patches[1].bstress_sal_ft[1]
            end

            fresh = _run(0.05)   # low salinity  (less stress -> larger factor)
            salty = _run(1.0)    # high salinity (more stress -> smaller factor)
            @test fresh > salty
            @test fresh > 0.0
            @test salty > 0.0
        end

        # ---------------------------------------------------------------------
        # multi-patch: every patch is processed (Fortran has no bareground skip)
        # ---------------------------------------------------------------------
        @testset "btran_sal_stress_fates! all patches processed" begin
            CLM.numpft[] = 1
            nlevsoil = 3

            edpft = CLM.EDPftvarcon_type()
            edpft.smpso = [-66000.0]; edpft.smpsc = [-255000.0]
            CLM.EDPftvarcon_inst[] = edpft
            _setup_prt_params!(1)

            salinity = [0.1, 0.2, 0.5]
            site, bc_in, patches = _build(nlevsoil;
                salinity = salinity,
                h2o      = fill(0.3, nlevsoil),
                tempk    = fill(sftz + 5.0, nlevsoil),
                npatch   = 3)

            CLM.btran_sal_stress_fates!([site], [bc_in])

            rootfr = zeros(Float64, nlevsoil)
            CLM.set_root_fraction(rootfr, 1, site.zi_soil; max_nlevroot=nlevsoil)
            expect = sum(rootfr[j] * _sal_rresis(salinity[j]) for j in 1:nlevsoil)

            # All three patches computed -> none left at the NaN default.
            for p in patches
                @test !isnan(p.bstress_sal_ft[1])
                @test isapprox(p.bstress_sal_ft[1], expect; rtol=1e-12)
            end
        end

    finally
        CLM.numpft[]           = old_numpft
        CLM.EDPftvarcon_inst[] = old_edpft
    end
end
