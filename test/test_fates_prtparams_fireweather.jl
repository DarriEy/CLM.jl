# test_fates_prtparams_fireweather.jl
# Tests for FATES Batch 2: PRTParametersMod (PARTEH/allometry params) and
# SFFireWeatherMod (SPITFIRE fire-weather base type).

using Test
using CLM

@testset "FATES PRTParameters + SFFireWeather" begin

    @testset "prt_param_type allocate + shapes" begin
        # Empty construction: every array starts zero-length.
        p0 = CLM.prt_param_type()
        @test isempty(p0.stress_decid)
        @test isempty(p0.allom_agb1)
        @test size(p0.leaf_long) == (0, 0)
        @test size(p0.nitr_stoich_p1) == (0, 0)

        # Small synthetic dimensions.
        npft     = 4
        norgan   = 3
        nleafage = 2
        nall     = 6   # all PRT organs (for organ_param_id reverse lookup)

        p = CLM.prt_param_type()
        CLM.allocate_prt_params!(p, npft, norgan, nleafage; nall_organs = nall)

        # PFT-indexed vectors.
        for fld in (:stress_decid, :season_decid, :evergreen,
                    :phen_fnrt_drop_fraction, :phen_stem_drop_fraction,
                    :phen_drought_threshold, :phen_moist_threshold, :phen_doff_time,
                    :senleaf_long_fdrought, :root_long, :branch_long,
                    :leafn_vert_scaler_coeff1, :leafn_vert_scaler_coeff2, :grperc,
                    :nitr_store_ratio, :phos_store_ratio, :cushion,
                    :leaf_stor_priority, :dbh_repro_threshold, :seed_alloc_mature,
                    :seed_alloc, :repro_alloc_a, :repro_alloc_b,
                    :fnrt_prof_mode, :fnrt_prof_a, :fnrt_prof_b, :c2b,
                    :wood_density, :woody, :slamax, :slatop, :allom_sai_scaler,
                    :allom_dbh_maxheight, :allom_hmode, :allom_lmode, :allom_fmode,
                    :allom_amode, :allom_cmode, :allom_smode, :allom_stmode,
                    :allom_dmode, :allom_la_per_sa_int, :allom_la_per_sa_slp,
                    :allom_l2fr, :allom_agb_frac, :allom_d2h1, :allom_d2h2,
                    :allom_d2h3, :allom_d2bl1, :allom_d2bl2, :allom_d2bl3,
                    :allom_blca_expnt_diff, :allom_d2ca_coefficient_max,
                    :allom_d2ca_coefficient_min, :allom_agb1, :allom_agb2,
                    :allom_agb3, :allom_agb4, :allom_h2cd1, :allom_h2cd2,
                    :allom_zroot_max_dbh, :allom_zroot_max_z, :allom_zroot_min_dbh,
                    :allom_zroot_min_z, :allom_zroot_k, :pid_kp, :pid_ki, :pid_kd,
                    :store_ovrflw_frac, :nfix_mresp_scfrac)
            v = getfield(p, fld)
            @test length(v) == npft
            @test all(iszero, v)
        end

        # (PFT x organ) matrices.
        for fld in (:turnover_nitr_retrans, :turnover_phos_retrans,
                    :nitr_stoich_p1, :phos_stoich_p1, :alloc_priority)
            @test size(getfield(p, fld)) == (npft, norgan)
        end

        # (PFT x age-class) matrices.
        @test size(p.leaf_long) == (npft, nleafage)
        @test size(p.leaf_long_ustory) == (npft, nleafage)

        # organ-indexed and all-PRT-organs vectors.
        @test length(p.organ_id) == norgan
        @test length(p.organ_param_id) == nall

        # Integer arrays really are integers.
        @test eltype(p.stress_decid) == Int
        @test eltype(p.woody) == Int
        @test eltype(p.allom_hmode) == Int
        @test eltype(p.grperc) == Float64

        # Default nall_organs == norgan.
        p2 = CLM.prt_param_type()
        CLM.allocate_prt_params!(p2, 2, 5, 1)
        @test length(p2.organ_param_id) == 5

        # Module-level singleton exists.
        @test CLM.prt_params isa CLM.prt_param_type
    end

    @testset "fire_weather base type + accessors" begin
        @test isabstracttype(CLM.fire_weather)
        @test CLM.base_fire_weather <: CLM.fire_weather

        # Wind attenuation constants.
        @test CLM.wind_atten_treed == 0.4
        @test CLM.wind_atten_grass == 0.6

        fw = CLM.base_fire_weather()
        @test fw.fire_weather_index == 0.0
        @test fw.effective_windspeed == 0.0

        # Init resets state.
        fw.fire_weather_index = 5.0
        fw.effective_windspeed = 9.0
        CLM.init_fire_weather!(fw)
        @test fw.fire_weather_index == 0.0
        @test fw.effective_windspeed == 0.0

        # UpdateEffectiveWindSpeed: all-tree -> wind * 0.4
        CLM.update_effective_windspeed!(fw, 10.0, 1.0, 0.0, 0.0)
        @test fw.effective_windspeed ≈ 10.0 * 0.4

        # all-grass -> wind * 0.6
        CLM.update_effective_windspeed!(fw, 10.0, 0.0, 1.0, 0.0)
        @test fw.effective_windspeed ≈ 10.0 * 0.6

        # all-bare -> wind * 0.6 (bare grouped with grass)
        CLM.update_effective_windspeed!(fw, 10.0, 0.0, 0.0, 1.0)
        @test fw.effective_windspeed ≈ 10.0 * 0.6

        # mixed cover.
        CLM.update_effective_windspeed!(fw, 8.0, 0.5, 0.3, 0.2)
        @test fw.effective_windspeed ≈ 8.0 * (0.5 * 0.4 + (0.3 + 0.2) * 0.6)

        # Deferred update_index! errors on the bare base subtype.
        @test_throws ErrorException CLM.update_index!(fw, 20.0, 0.0, 50.0, 3.0)
    end

end
