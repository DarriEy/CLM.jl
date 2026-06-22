# test_fates_edparams.jl
# Tests for FATES EDParamsMod (Tier F, Batch 1).

using Test
using CLM

@testset "FATES EDParamsMod" begin

    @testset "compile-time parameters" begin
        @test CLM.soil_tfrz_thresh == -2.0
        @test CLM.nclmax == 2
        @test CLM.nlevleaf == 30
        @test CLM.maxpft == 16
    end

    @testset "parameter names preserved" begin
        @test CLM.fates_name_active_crown_fire == "fates_fire_active_crown_fire"
        @test CLM.fates_name_cg_strikes == "fates_fire_cg_strikes"
        @test CLM.name_dev_arbitrary == "fates_dev_arbitrary"
        @test CLM.ED_name_photo_temp_acclim_timescale == "fates_leaf_photo_temp_acclim_timescale"
        @test CLM.name_maintresp_model == "fates_maintresp_leaf_model"
        @test CLM.name_radiation_model == "fates_rad_model"
        @test CLM.hydr_name_solver == "fates_hydro_solver"
        @test CLM.maxcohort_name == "fates_maxcohort"
        @test CLM.eca_name_plant_escalar == "fates_cnp_eca_plant_escalar"
        @test CLM.ED_name_history_sizeclass_bin_edges == "fates_history_sizeclass_bin_edges"
        @test CLM.ED_name_maxpatches_by_landuse == "fates_maxpatches_by_landuse"
    end

    @testset "FatesParamsInit! defaults" begin
        CLM.FatesParamsInit!()
        p = CLM.ed_params()
        # reals reset to the unset sentinel
        @test p.photo_temp_acclim_timescale == CLM.fates_unset_r8
        @test p.ED_val_phen_a == CLM.fates_unset_r8
        @test p.hydr_kmax_rsurf1 == CLM.fates_unset_r8
        @test p.logging_dbhmin == CLM.fates_unset_r8
        # integer switches reset to -9
        @test p.photo_tempsens_model == -9
        @test p.radiation_model == -9
        @test p.stomatal_model == -9
        @test p.hydr_solver == -9
        @test p.damage_event_code == -9
        # bool default
        @test p.active_crown_fire == false
        # VAI bin arrays sized nlevleaf, filled with the unset sentinel
        @test length(p.dinc_vai) == CLM.nlevleaf
        @test length(p.dlower_vai) == CLM.nlevleaf
        @test all(p.dinc_vai .== CLM.fates_unset_r8)
        # landuse-sized vectors
        @test length(p.maxpatches_by_landuse) == CLM.n_landuse_cats
        @test length(p.crop_lu_pft_vector) == CLM.n_landuse_cats
    end

    @testset "FatesRegisterParams! registers expected names" begin
        fp = CLM.fates_parameters_type()
        CLM.FatesRegisterParams!(fp)

        # 60 scalar + 9 non-scalar (1-D) registrations
        @test fp.num_parameters == 69

        registered = [fp.parameters[i].name for i in 1:fp.num_parameters]

        # a sampling across categories must all be present
        for nm in (CLM.ED_name_photo_temp_acclim_timescale,
                   CLM.name_radiation_model,
                   CLM.name_maintresp_model,
                   CLM.ED_name_phen_a,
                   CLM.ED_name_stomatal_model,
                   CLM.maxcohort_name,
                   CLM.hydr_name_solver,
                   CLM.bgc_name_soil_salinity,
                   CLM.logging_name_dbhmin,
                   CLM.eca_name_plant_escalar,
                   CLM.fates_name_active_crown_fire,
                   CLM.fates_name_cg_strikes,
                   CLM.damage_name_event_code,
                   # non-scalar
                   CLM.ED_name_hydr_htftype_node,
                   CLM.ED_name_history_sizeclass_bin_edges,
                   CLM.ED_name_crop_lu_pft_vector,
                   CLM.ED_name_maxpatches_by_landuse,
                   CLM.ED_name_max_nocomp_pfts_by_landuse)
            @test nm in registered
        end

        # the 1-D params carry the right dimension shape
        i_hist = CLM.FindIndex(fp, CLM.ED_name_history_sizeclass_bin_edges)
        @test fp.parameters[i_hist].dimension_shape == CLM.dimension_shape_1d
        i_scalar = CLM.FindIndex(fp, CLM.ED_name_phen_a)
        @test fp.parameters[i_scalar].dimension_shape == CLM.dimension_shape_scalar

        # FatesRegisterParams! resets the params first (calls FatesParamsInit!)
        @test CLM.ed_params().stomatal_model == -9
    end

    @testset "FatesReceiveParams! round-trip + derived values" begin
        fp = CLM.fates_parameters_type()
        CLM.FatesRegisterParams!(fp)

        # populate stored data for every registered parameter
        function setscalar!(name, val)
            idx = CLM.FindIndex(fp, name)
            CLM.SetDataScalar!(fp, idx, val)
        end
        function set1d!(name, vec)
            idx = CLM.FindIndex(fp, name)
            fp.parameters[idx].dimension_sizes[1] = length(vec)
            CLM.SetData1D!(fp, idx, vec)
        end

        # default every scalar to something finite, then override the ones we check
        for i in 1:fp.num_parameters
            if fp.parameters[i].dimension_shape == CLM.dimension_shape_scalar
                CLM.SetDataScalar!(fp, i, 0.0)
            end
        end

        setscalar!(CLM.ED_name_photo_temp_acclim_timescale, 30.0)
        setscalar!(CLM.ED_name_phen_a, -68.0)
        setscalar!(CLM.ED_name_cwd_fcel, 0.76)
        setscalar!(CLM.hydr_name_kmax_rsurf1, 20.0)
        setscalar!(CLM.eca_name_plant_escalar, 1.25e-5)
        # integer-coded scalars: stored as real, nint-rounded on receive
        setscalar!(CLM.name_radiation_model, 2.0)
        setscalar!(CLM.ED_name_stomatal_model, 1.0)
        setscalar!(CLM.maxcohort_name, 100.0)
        setscalar!(CLM.hydr_name_solver, 1.0)
        setscalar!(CLM.damage_name_event_code, 1.0)
        setscalar!(CLM.name_maintresp_model, 2.0)
        # boolean: abs(val-1) < nearzero
        setscalar!(CLM.fates_name_active_crown_fire, 1.0)

        # 1-D parameters
        set1d!(CLM.ED_name_hydr_htftype_node, Float64[1.0, 2.0, 1.0])
        set1d!(CLM.ED_name_history_sizeclass_bin_edges, Float64[0.0, 10.0, 50.0, 100.0])
        set1d!(CLM.ED_name_history_ageclass_bin_edges, Float64[0.0, 1.0, 5.0])
        set1d!(CLM.ED_name_history_height_bin_edges, Float64[0.0, 5.0])
        set1d!(CLM.ED_name_history_coageclass_bin_edges, Float64[0.0, 2.0])
        set1d!(CLM.ED_name_history_damage_bin_edges, Float64[0.0, 0.5])
        set1d!(CLM.ED_name_crop_lu_pft_vector, fill(0.0, CLM.n_landuse_cats))
        set1d!(CLM.ED_name_maxpatches_by_landuse,
               Float64[10.0, 4.0, 1.0, 1.0, 1.0])
        set1d!(CLM.ED_name_max_nocomp_pfts_by_landuse,
               fill(2.0, CLM.n_landuse_cats))

        CLM.FatesReceiveParams!(fp)
        p = CLM.ed_params()

        # scalars round-tripped
        @test p.photo_temp_acclim_timescale == 30.0
        @test p.ED_val_phen_a == -68.0
        @test p.ED_val_cwd_fcel == 0.76
        @test p.hydr_kmax_rsurf1 == 20.0
        @test p.eca_plant_escalar == 1.25e-5

        # integer-coded scalars nint-rounded to Int
        @test p.radiation_model === 2
        @test p.stomatal_model === 1
        @test p.max_cohort_per_patch === 100
        @test p.hydr_solver === 1
        @test p.damage_event_code === 1
        @test p.maintresp_leaf_model === 2

        # boolean parameter
        @test p.active_crown_fire === true

        # 1-D round-trips, integer vectors rounded
        @test p.hydr_htftype_node == [1, 2, 1]
        @test p.ED_val_history_sizeclass_bin_edges == [0.0, 10.0, 50.0, 100.0]
        @test p.maxpatches_by_landuse == [10, 4, 1, 1, 1]
        @test length(p.max_nocomp_pfts_by_landuse) == CLM.n_landuse_cats

        # derived: maxpatch_total = sum(maxpatches_by_landuse)
        @test p.maxpatch_total == sum(p.maxpatches_by_landuse)
        @test p.maxpatch_total == 17
    end

    @testset "FatesReportParams (no-op unless debug)" begin
        @test CLM.FatesReportParams(true) === nothing
        @test CLM.FatesReportParams(false) === nothing
    end

end
