using Test
using CLM

# ==========================================================================
# Regression guards for the two CN bugs that ONLY the autumn / leaf-offset
# window could expose (a summer window has offset_flag == 0 everywhere, so it
# cannot execute a single line of the senescence path).
#
# See docs/BGC_PARITY_SCORECARD.md.
#
#  1. The offset litterfall RAMP MEMORY (prev_{leafc,frootc}_to_litter) must
#     survive the per-step CN flux reset. Fortran's
#     CNVegCarbonFluxType::SetValues does not touch them (they are set only in
#     InitCold and are restart vars); they are reset in exactly two places, both
#     in CNPhenologyMod: the offset trigger and the ramp-end cleanup.
#
#     Zeroing them each step collapses
#         leafc_to_litter = prev + t1*(leafc - prev*offset_counter)
#     to t1*leafc -- only ~2-7% of the correct flux mid-ramp.
#
#  2. The use_cn COLD START must branch on PFT type. Seeding every patch as
#     deciduous gave evergreen PFTs ZERO displayed leaf/root carbon.
# ==========================================================================

@testset "CN offset ramp memory + PFT-aware cold start" begin

    # ---- 1. prev_*_to_litter survives the per-step flux reset ----------------
    @testset "offset ramp memory survives the CN flux reset" begin
        np, nc, ng = 3, 1, 1
        cf = CLM.CNVegCarbonFluxData()
        CLM.cnveg_carbon_flux_init!(cf, np, nc, ng)
        CLM.cnveg_carbon_flux_init_cold!(cf, 1:np, 1:nc)

        # Mid-ramp memory, as CNOffsetLitterfall would have left it.
        cf.prev_leafc_to_litter_patch[2]  = 1.6549e-5
        cf.prev_frootc_to_litter_patch[2] = 3.1e-6
        # A transient flux that SHOULD be cleared each step.
        cf.gpp_patch[2] = 12.34

        CLM._zero_cnveg_flux_arrays!(cf)

        @test cf.prev_leafc_to_litter_patch[2]  == 1.6549e-5   # memory preserved
        @test cf.prev_frootc_to_litter_patch[2] == 3.1e-6
        @test cf.gpp_patch[2] == 0.0                            # ordinary flux zeroed
    end

    # ---- 2. the ramp formula, given the memory, reproduces Fortran ----------
    # Ground truth from the Fortran autumn reference (Bow, grass pft 12,
    # nstep 1760041 = 2202-10-15 00:00):
    #   leafc = 31.8047, offset_counter (post-decrement) = 900000 s,
    #   prev_leafc_to_litter = 1.654925296836598e-5
    #   => Fortran leafc_to_litter = 1.6700e-05
    @testset "offset ramp formula matches the Fortran reference" begin
        dt   = 3600.0
        leafc = 31.804723
        oc    = 900000.0
        prev  = 1.654925296836598e-5

        t1 = dt * 2.0 / (oc * oc)
        leafc_to_litter = prev + t1 * (leafc - prev * oc)
        @test isapprox(leafc_to_litter, 1.6700e-5; rtol = 1e-4)

        # The bug: with prev dropped, the flux is a tiny fraction of the truth.
        broken_flux = t1 * leafc
        @test broken_flux / leafc_to_litter < 0.02      # ~1.7% -- 50x too small
    end

    # ---- 3. cold start branches on PFT type --------------------------------
    @testset "cold start: evergreen vs deciduous vs bare" begin
        np = 3
        itype = [0, 1, 12]          # bare, needleleaf EVERGREEN tree, C3 arctic grass
        nveg  = 25
        evergreen = zeros(nveg); evergreen[1 + 1] = 1.0   # pft 1 is evergreen
        woody     = zeros(nveg); woody[1 + 1]     = 1.0
        leafcn    = fill(30.0, nveg); leafcn[1 + 1] = 58.0;  leafcn[12 + 1] = 20.7
        frootcn   = fill(42.0, nveg)
        deadwdcn  = fill(500.0, nveg)

        cs = CLM.CNVegCarbonStateData()
        CLM.cnveg_carbon_state_init!(cs, np, 1, 1)
        CLM.cnveg_carbon_state_init_cold!(cs, 1:np;
            patch_itype = itype, evergreen = evergreen, woody = woody,
            initial_vegC = 100.0)

        # bare: nothing at all
        @test cs.leafc_patch[1] == 0.0 && cs.leafc_storage_patch[1] == 0.0
        @test cs.frootc_patch[1] == 0.0 && cs.frootc_storage_patch[1] == 0.0
        @test cs.deadstemc_patch[1] == 0.0

        # EVERGREEN: displayed leaf/root C, no storage (this is what was broken --
        # an evergreen PFT used to cold-start with leafc == 0, i.e. no canopy)
        @test cs.leafc_patch[2]          == 100.0
        @test cs.leafc_storage_patch[2]  == 0.0
        @test cs.frootc_patch[2]         == 100.0
        @test cs.frootc_storage_patch[2] == 0.0
        @test cs.deadstemc_patch[2]      == 0.1     # woody seed

        # DECIDUOUS: all in storage
        @test cs.leafc_patch[3]          == 0.0
        @test cs.leafc_storage_patch[3]  == 100.0
        @test cs.frootc_patch[3]         == 0.0
        @test cs.frootc_storage_patch[3] == 100.0
        @test cs.deadstemc_patch[3]      == 0.0     # not woody

        # N is DERIVED from the cold C pools via the PFT C:N ratios.
        ns = CLM.CNVegNitrogenStateData()
        CLM.cnveg_nitrogen_state_init!(ns, np, 1, 1)
        CLM.cnveg_nitrogen_state_init_cold!(ns, 1:np;
            carbonstate = cs, patch_itype = itype,
            leafcn = leafcn, frootcn = frootcn, deadwdcn = deadwdcn, woody = woody)

        @test ns.leafn_patch[1] == 0.0
        @test isapprox(ns.leafn_patch[2],          100.0 / 58.0;  rtol = 1e-12)  # 1.7241
        @test ns.leafn_storage_patch[2] == 0.0
        @test isapprox(ns.frootn_patch[2],         100.0 / 42.0;  rtol = 1e-12)  # 2.3810
        @test isapprox(ns.deadstemn_patch[2],      0.1   / 500.0; rtol = 1e-12)
        @test ns.leafn_patch[3] == 0.0
        @test isapprox(ns.leafn_storage_patch[3],  100.0 / 20.7;  rtol = 1e-12)  # 4.8309
        @test isapprox(ns.frootn_storage_patch[3], 100.0 / 42.0;  rtol = 1e-12)
        @test ns.deadstemn_patch[3] == 0.0
    end

    # ---- 4. initial_vegC is a namelist knob, defaulting to the CTSM value ----
    @testset "initial_vegC default is the CTSM default (20)" begin
        @test CLM.cnvegcstate_const.initial_vegC == CLM.INITIAL_VEGC == 20.0
    end
end
