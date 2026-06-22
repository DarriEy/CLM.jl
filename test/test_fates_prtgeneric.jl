# test_fates_prtgeneric.jl
# Tests for FATES Batch 2: PRTGenericMod — the generic PARTEH framework.
# Exercises the element/organ index constants + counts, builds a small
# prt_global_type descriptor, constructs the generic prt_vartypes container,
# registers a few organ/element pools, and verifies the generic
# GetState/SetState round-trips and the mass accessors.

using Test
using CLM

@testset "FATES PRTGenericMod" begin

    @testset "organ / element index constants + counts" begin
        # Organ indices
        @test CLM.all_organs    == 0
        @test CLM.leaf_organ    == 1
        @test CLM.fnrt_organ    == 2
        @test CLM.sapw_organ    == 3
        @test CLM.store_organ   == 4
        @test CLM.repro_organ   == 5
        @test CLM.struct_organ  == 6
        @test CLM.num_organ_types == 6

        # Element indices
        @test CLM.carbon12_element   == 1
        @test CLM.carbon13_element   == 2
        @test CLM.carbon14_element   == 3
        @test CLM.nitrogen_element   == 4
        @test CLM.phosphorus_element == 5
        @test CLM.potassium_element  == 6
        @test CLM.num_element_types  == 6

        # Carbon element list + max group size
        @test CLM.carbon_elements_list ==
              [CLM.carbon12_element, CLM.carbon13_element, CLM.carbon14_element]
        @test CLM.max_spec_per_group == 3
        @test length(CLM.carbon_elements_list) <= CLM.max_spec_per_group

        # Hypothesis ids
        @test CLM.prt_carbon_allom_hyp   == 1
        @test CLM.prt_cnp_flex_allom_hyp == 2

        # Misc parameters
        @test CLM.max_nleafage == 4
        @test CLM.l2fr_min ≈ 0.01
        @test CLM.mass_unit == "kg"
        @test CLM.mass_rate_unit == "kg/day"
        @test CLM.un_initialized < CLM.check_initialized   # -9.9e32 < -8.8e32
    end

    @testset "prt_global_type: ZeroGlobal! + RegisterVarInGlobal!" begin
        g = CLM.prt_global_type()
        CLM.ZeroGlobal!(g)

        # After zeroing, the BC + var counts are bogus (-9)
        @test g.num_bc_in    == -9
        @test g.num_bc_out   == -9
        @test g.num_bc_inout == -9
        @test g.num_vars     == -9

        # Maps zeroed
        @test all(CLM.sp_organ_map_get(g, io, is) == 0
                  for io in 1:CLM.num_organ_types, is in 1:CLM.num_element_types)

        # Build a small carbon-only hypothesis: leaf, fineroot, sapwood, structure.
        # Use multiple leaf age positions to exercise the sum-over-positions paths.
        nleaf_pos = 2
        g.hyp_name = "test_carbon"
        g.hyp_id   = CLM.prt_carbon_allom_hyp
        g.num_vars     = 4
        g.num_bc_in    = 0
        g.num_bc_out   = 0
        g.num_bc_inout = 0
        g.state_descriptor = [CLM.state_descriptor_type() for _ in 1:g.num_vars]

        el = CLM.carbon12_element
        CLM.RegisterVarInGlobal!(g, 1, "leaf carbon",   "leaf_c",   CLM.leaf_organ,   el, nleaf_pos)
        CLM.RegisterVarInGlobal!(g, 2, "fineroot carbon", "fnrt_c", CLM.fnrt_organ,   el, 1)
        CLM.RegisterVarInGlobal!(g, 3, "sapwood carbon", "sapw_c",  CLM.sapw_organ,   el, 1)
        CLM.RegisterVarInGlobal!(g, 4, "structure carbon", "str_c", CLM.struct_organ, el, 1)

        # sp_organ_map now maps each (organ, element) -> var id
        @test CLM.sp_organ_map_get(g, CLM.leaf_organ,   el) == 1
        @test CLM.sp_organ_map_get(g, CLM.fnrt_organ,   el) == 2
        @test CLM.sp_organ_map_get(g, CLM.sapw_organ,   el) == 3
        @test CLM.sp_organ_map_get(g, CLM.struct_organ, el) == 4

        # Each organ map collected exactly one variable
        @test g.organ_map[CLM.leaf_organ].num_vars == 1
        @test g.organ_map[CLM.leaf_organ].var_id[1] == 1
        @test g.organ_map[CLM.struct_organ].var_id[1] == 4

        # Descriptors recorded
        @test g.state_descriptor[1].longname == "leaf carbon"
        @test g.state_descriptor[1].symbol   == "leaf_c"
        @test g.state_descriptor[1].num_pos  == nleaf_pos
        @test g.state_descriptor[1].organ_id == CLM.leaf_organ
        @test g.state_descriptor[1].element_id == el

        # ----- Install as the module global, then build a plant -----
        old_global = CLM.prt_global[]
        CLM.prt_global[] = g
        try
            prt = CLM.prt_vartypes()
            CLM.InitPRTVartype!(prt)

            # Allocation sized each variable's SoA vectors by num_pos
            @test length(prt.variables) == g.num_vars
            @test length(prt.variables[1].val) == nleaf_pos
            @test length(prt.variables[2].val) == 1
            @test prt.ode_opt_step ≈ 1e6

            # Initial conditions are the "un-initialized" sentinel
            @test all(prt.variables[1].val .== CLM.un_initialized)

            # Before all states are set, CheckInitialConditions must fail
            @test_throws Exception CLM.CheckInitialConditions(prt)

            # ----- SetState! / GetState round-trips -----
            CLM.SetState!(prt, CLM.leaf_organ, el, 1.5, 1)   # leaf bin 1
            CLM.SetState!(prt, CLM.leaf_organ, el, 0.5, 2)   # leaf bin 2
            CLM.SetState!(prt, CLM.fnrt_organ, el, 2.0)      # default position 1
            CLM.SetState!(prt, CLM.sapw_organ, el, 3.0)
            CLM.SetState!(prt, CLM.struct_organ, el, 4.0)

            # Position-specific read
            @test CLM.GetState(prt, CLM.leaf_organ, el, 1) ≈ 1.5
            @test CLM.GetState(prt, CLM.leaf_organ, el, 2) ≈ 0.5
            # Sum over all positions
            @test CLM.GetState(prt, CLM.leaf_organ, el) ≈ 2.0
            @test CLM.GetState(prt, CLM.fnrt_organ, el) ≈ 2.0
            @test CLM.GetState(prt, CLM.struct_organ, el) ≈ 4.0

            # All states now set -> CheckInitialConditions passes
            CLM.SetState!(prt, CLM.fnrt_organ, el, 2.0)  # already set; keep
            @test CLM.CheckInitialConditions(prt) === nothing

            # SetState! out-of-range position errors
            @test_throws Exception CLM.SetState!(prt, CLM.fnrt_organ, el, 1.0, 5)
            # SetState! on a non-existent organ/element combo errors
            @test_throws Exception CLM.SetState!(prt, CLM.repro_organ, el, 1.0)

            # ----- ZeroRates! then exercise the flux accessors + conservation -----
            CLM.ZeroRates!(prt)
            # val0 snapshot equals val; rates zeroed
            @test prt.variables[1].val0[1] ≈ 1.5
            @test all(prt.variables[1].net_alloc .== 0.0)
            @test CLM.GetTurnover(prt, CLM.leaf_organ, el) ≈ 0.0
            @test CLM.GetBurned(prt, CLM.leaf_organ, el) ≈ 0.0
            @test CLM.GetNetAlloc(prt, CLM.leaf_organ, el) ≈ 0.0

            # A consistent flux/state update conserves mass:
            #   val - val0 == net_alloc - turnover - burned - damaged
            prt.variables[2].net_alloc[1] = 0.7
            prt.variables[2].turnover[1]  = 0.2
            prt.variables[2].val[1] = prt.variables[2].val0[1] + 0.7 - 0.2
            @test CLM.GetNetAlloc(prt, CLM.fnrt_organ, el) ≈ 0.7
            @test CLM.GetTurnover(prt, CLM.fnrt_organ, el) ≈ 0.2
            @test CLM.CheckMassConservation(prt, 1, 1) === nothing

            # Break conservation -> error
            prt.variables[2].val[1] += 1.0
            @test_throws Exception CLM.CheckMassConservation(prt, 1, 2)
            prt.variables[2].val[1] -= 1.0  # restore

            # ----- Copy + weighted fuse -----
            donor = CLM.prt_vartypes()
            CLM.InitPRTVartype!(donor)
            CLM.CopyPRTVartypes!(donor, prt)
            @test CLM.GetState(donor, CLM.leaf_organ, el) ≈ CLM.GetState(prt, CLM.leaf_organ, el)
            @test donor.ode_opt_step ≈ prt.ode_opt_step

            # Set the donor leaf bin-1 to a different value, fuse 50/50
            CLM.SetState!(donor, CLM.leaf_organ, el, 3.5, 1)
            recip = CLM.prt_vartypes()
            CLM.InitPRTVartype!(recip)
            CLM.CopyPRTVartypes!(recip, prt)   # recip leaf bin1 = 1.5
            CLM.WeightedFusePRTVartypes!(recip, donor, 0.5)
            @test recip.variables[1].val[1] ≈ 0.5 * 1.5 + 0.5 * 3.5

            # ----- Deallocate drops the arrays -----
            CLM.DeallocatePRTVartypes!(prt)
            @test isempty(prt.variables)
        finally
            CLM.prt_global[] = old_global
        end
    end

    @testset "deferred (must-extend) base methods error" begin
        prt = CLM.prt_vartypes()
        @test_throws Exception CLM.DailyPRT!(prt, 1)
        @test_throws Exception CLM.FastPRT!(prt)
        @test_throws Exception CLM.DamageRecovery!(prt)
        @test_throws Exception CLM.GetNutrientTarget(prt, CLM.nitrogen_element, CLM.leaf_organ)
        @test_throws Exception CLM.GetCoordVal(prt, CLM.leaf_organ, CLM.carbon12_element)
    end

    @testset "StorageNutrientTarget (leaf-proportional default)" begin
        # store_prop defaults to lf_store_prop -> nitr_store_ratio[pft] * leaf_target
        npft = 2
        CLM.allocate_prt_params!(CLM.prt_params, npft, CLM.num_organ_types, CLM.max_nleafage)
        CLM.prt_params.nitr_store_ratio[1] = 2.0
        CLM.prt_params.phos_store_ratio[1] = 0.5

        @test CLM.StorageNutrientTarget(1, CLM.nitrogen_element, 1.0, 9.0, 9.0, 9.0) ≈ 2.0
        @test CLM.StorageNutrientTarget(1, CLM.phosphorus_element, 4.0, 9.0, 9.0, 9.0) ≈ 2.0

        # Carbon is not allowed
        @test_throws Exception CLM.StorageNutrientTarget(1, CLM.carbon12_element,
                                                         1.0, 1.0, 1.0, 1.0)
    end
end
