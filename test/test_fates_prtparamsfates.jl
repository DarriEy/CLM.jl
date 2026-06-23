# test_fates_prtparamsfates.jl
# Tests for FATES PRTParamsFATESMod (Tier F, Batch 11) — the PARTEH parameter
# initialization / registration layer. Drives PRTRegisterParams! against a
# synthetic param reader and asserts the expected fates_* names + counts;
# populates the reader, runs PRTReceiveParams! into the merged prt_params for a
# small synthetic dimension set and asserts array shapes; drives PRTDerivedParams!
# and a couple of PRTCheckParams consistency throws.

using Test
using CLM

@testset "FATES PRTParamsFATESMod" begin

    npft     = 2   # small synthetic PFT count
    norgan   = 2   # parameter-file organ count (leaf, fnrt)
    nleafage = 1   # leaf-age classes

    # -----------------------------------------------------------------------
    @testset "lower-bound constants" begin
        @test CLM.lower_bound_pft == 1
        @test CLM.lower_bound_general == 1
    end

    # -----------------------------------------------------------------------
    @testset "ArrayNint rounds to nearest int" begin
        @test CLM.ArrayNint([0.0, 0.4, 0.6, 1.5, 2.4]) == [0, 0, 1, 2, 2]
        @test eltype(CLM.ArrayNint([1.0, 2.0])) == Int
    end

    # -----------------------------------------------------------------------
    @testset "PRTRegisterParams! lists expected fates_* names + count" begin
        fp = CLM.fates_parameters_type()
        CLM.PRTRegisterParams!(fp)

        registered = Set(fp.parameters[i].name for i in 1:fp.num_parameters)

        # representative names across every registration group, verbatim
        expected = [
            "fates_phen_stress_decid",          # PFT int
            "fates_phen_evergreen",
            "fates_woody",
            "fates_allom_hmode",
            "fates_leaf_slatop",                # PFT real
            "fates_grperc",
            "fates_allom_agb1",
            "fates_cnp_nfix1",
            "fates_turnover_branch",
            "fates_turnover_leaf_canopy",       # PFT x leaf-age (2-D)
            "fates_turnover_leaf_ustory",
            "fates_stoich_nitr",                # PFT x organ (2-D)
            "fates_alloc_organ_priority",
            "fates_cnp_turnover_phos_retrans",
            "fates_alloc_organ_id",             # organ (1-D)
        ]
        for nm in expected
            @test nm in registered
        end

        # count = ints + reals + 2 leaf-age + 5 PFT-organ + 1 organ-id
        @test fp.num_parameters ==
              length(CLM._prt_pft_int_params) +
              length(CLM._prt_pft_real_params) +
              2 + length(CLM._prt_pft_organ_params) + 1

        # the dimension shapes are set correctly for representative names
        idx1 = CLM.FindIndex(fp, "fates_woody")
        @test fp.parameters[idx1].dimension_shape == CLM.dimension_shape_1d
        @test fp.parameters[idx1].dimension_names[1] == CLM.dimension_name_pft

        idx2 = CLM.FindIndex(fp, "fates_stoich_nitr")
        @test fp.parameters[idx2].dimension_shape == CLM.dimension_shape_2d
        @test fp.parameters[idx2].dimension_names[1] == CLM.dimension_name_pft
        @test fp.parameters[idx2].dimension_names[2] == CLM.dimension_name_prt_organs

        idxo = CLM.FindIndex(fp, "fates_alloc_organ_id")
        @test fp.parameters[idxo].dimension_shape == CLM.dimension_shape_1d
        @test fp.parameters[idxo].dimension_names[1] == CLM.dimension_name_prt_organs
    end

    # -----------------------------------------------------------------------
    # Build a fully-populated reader for the synthetic dims, run PRTReceiveParams!,
    # and check the prt_params shapes. We populate every registered parameter with
    # finite, physically-valid data (so PRTCheckParams will pass afterwards).
    #
    # `setup_reader` returns a reader whose data are valid for carbon-only mode.
    function setup_reader()
        fp = CLM.fates_parameters_type()
        CLM.PRTRegisterParams!(fp)
        for i in 1:fp.num_parameters
            par = fp.parameters[i]
            nm = par.name
            if par.dimension_shape == CLM.dimension_shape_1d
                if nm == "fates_alloc_organ_id"
                    # organ-dimensioned: leaf_organ=1, fnrt_organ=2
                    par.dimension_sizes[1] = norgan
                    CLM.SetData1D!(fp, i, Float64[CLM.leaf_organ, CLM.fnrt_organ])
                else
                    par.dimension_sizes[1] = npft
                    # choose a generally-valid default for PFT reals/ints
                    val = if nm == "fates_woody"
                        1.0                       # woody
                    elseif nm == "fates_phen_evergreen"
                        1.0                       # evergreen (mutually exclusive)
                    elseif nm in ("fates_phen_stress_decid", "fates_phen_season_decid")
                        0.0
                    elseif nm == "fates_allom_hmode"
                        1.0                       # O'Brien (so d2h2 sign unconstrained)
                    elseif nm in ("fates_allom_lmode", "fates_allom_fmode",
                                  "fates_allom_amode", "fates_allom_cmode",
                                  "fates_allom_smode", "fates_allom_stmode",
                                  "fates_allom_dmode", "fates_allom_fnrt_prof_mode")
                        1.0
                    elseif nm == "fates_allom_agb1"
                        0.05                      # non-zero AGB intercept (woody)
                    elseif nm == "fates_allom_d2h2"
                        0.5                       # positive ok for hmode=1
                    elseif nm == "fates_grperc"
                        0.11                      # in [0,1]
                    elseif nm == "fates_alloc_store_priority_frac"
                        0.5                       # leaf_stor_priority in [0,1]
                    elseif nm == "fates_turnover_senleaf_fdrought"
                        1.0                       # in (0,1]
                    elseif nm in ("fates_turnover_fnrt", "fates_turnover_branch")
                        1.0                       # 1 yr turnover (< 1/day rate)
                    elseif nm == "fates_recruit_seed_alloc"
                        0.1
                    elseif nm == "fates_recruit_seed_alloc_mature"
                        0.1                       # sum < 1
                    elseif nm == "fates_phen_stem_drop_fraction"
                        0.0                       # woody must be ~0
                    elseif nm == "fates_phen_fnrt_drop_fraction"
                        0.5
                    elseif nm == "fates_cnp_nfix1"
                        0.0                       # in [0,1]
                    elseif nm in ("fates_cnp_nitr_store_ratio", "fates_cnp_phos_store_ratio")
                        1.5
                    else
                        0.1
                    end
                    CLM.SetData1D!(fp, i, fill(val, npft))
                end
            elseif par.dimension_shape == CLM.dimension_shape_2d
                ncol = (par.name in ("fates_turnover_leaf_canopy",
                                     "fates_turnover_leaf_ustory")) ? nleafage : norgan
                # stoich/retrans/priority valid: in [0,1]; leaf_long = 1 yr
                fillval = if par.name in ("fates_turnover_leaf_canopy",
                                          "fates_turnover_leaf_ustory")
                    1.0                            # 1 yr (turnover < 1/day)
                elseif par.name == "fates_alloc_organ_priority"
                    1.0
                else
                    0.05                           # stoich + retrans in [0,1]
                end
                CLM.SetData2D!(fp, i, fill(fillval, npft, ncol))
            end
        end
        return fp
    end

    @testset "PRTReceiveParams! shapes" begin
        fp = setup_reader()
        CLM.PRTReceiveParams!(fp)
        p = CLM.prt_params

        # 1-D PFT arrays
        @test length(p.evergreen) == npft
        @test eltype(p.evergreen) == Int
        @test length(p.woody) == npft
        @test eltype(p.woody) == Int
        @test length(p.grperc) == npft
        @test eltype(p.grperc) == Float64
        @test length(p.allom_hmode) == npft
        @test eltype(p.allom_hmode) == Int

        # 1-D organ array
        @test length(p.organ_id) == norgan
        @test p.organ_id == [CLM.leaf_organ, CLM.fnrt_organ]

        # 2-D PFT x organ
        @test size(p.nitr_stoich_p1) == (npft, norgan)
        @test size(p.turnover_nitr_retrans) == (npft, norgan)
        @test size(p.alloc_priority) == (npft, norgan)

        # 2-D PFT x leaf-age
        @test size(p.leaf_long) == (npft, nleafage)
        @test size(p.leaf_long_ustory) == (npft, nleafage)
    end

    # -----------------------------------------------------------------------
    @testset "PRTDerivedParams! reverse organ lookup" begin
        fp = setup_reader()
        CLM.PRTReceiveParams!(fp)
        CLM.PRTDerivedParams!()
        p = CLM.prt_params

        @test length(p.organ_param_id) == CLM.num_organ_types
        # leaf_organ -> param index 1, fnrt_organ -> param index 2
        @test p.organ_param_id[CLM.leaf_organ] == 1
        @test p.organ_param_id[CLM.fnrt_organ] == 2
        # an organ not in the parameter file stays invalid (-1)
        @test p.organ_param_id[CLM.struct_organ] == -1
    end

    # -----------------------------------------------------------------------
    @testset "PRTCheckParams passes on valid carbon-only params" begin
        CLM.hlm_parteh_mode[] = CLM.prt_carbon_allom_hyp
        fp = setup_reader()
        CLM.PRTReceiveParams!(fp)
        CLM.PRTDerivedParams!()
        @test CLM.PRTCheckParams(true) === nothing
        # is_master=false is a no-op
        @test CLM.PRTCheckParams(false) === nothing
    end

    # -----------------------------------------------------------------------
    @testset "PRTCheckParams consistency throws" begin
        CLM.hlm_parteh_mode[] = CLM.prt_carbon_allom_hyp

        # (1) grperc out of [0,1] -> throw
        fp = setup_reader()
        CLM.PRTReceiveParams!(fp)
        CLM.PRTDerivedParams!()
        CLM.prt_params.grperc[1] = 1.5
        @test_throws ErrorException CLM.PRTCheckParams(true)

        # (2) seed allocation sum > 1 -> throw
        fp = setup_reader()
        CLM.PRTReceiveParams!(fp)
        CLM.PRTDerivedParams!()
        CLM.prt_params.seed_alloc[1] = 0.7
        CLM.prt_params.seed_alloc_mature[1] = 0.7
        @test_throws ErrorException CLM.PRTCheckParams(true)

        # (3) woody plant with zero AGB intercept -> throw
        fp = setup_reader()
        CLM.PRTReceiveParams!(fp)
        CLM.PRTDerivedParams!()
        CLM.prt_params.allom_agb1[1] = 0.0
        @test_throws ErrorException CLM.PRTCheckParams(true)

        # (4) mutually-exclusive phenology habit violated -> throw
        fp = setup_reader()
        CLM.PRTReceiveParams!(fp)
        CLM.PRTDerivedParams!()
        CLM.prt_params.evergreen[1] = CLM.itrue
        CLM.prt_params.season_decid[1] = CLM.itrue
        @test_throws ErrorException CLM.PRTCheckParams(true)

        # (5) invalid organ id (> num_organ_types) -> throw
        fp = setup_reader()
        CLM.PRTReceiveParams!(fp)
        CLM.PRTDerivedParams!()
        CLM.prt_params.organ_id[1] = CLM.num_organ_types + 1
        @test_throws ErrorException CLM.PRTCheckParams(true)
    end
end
