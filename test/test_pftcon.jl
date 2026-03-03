@testset "PFT Constants (pftconMod)" begin

    @testset "Module-level constants" begin
        # PFT index constants (Fortran 0-based values)
        @test CLM.noveg == 0
        @test CLM.ndllf_evr_tmp_tree == 1
        @test CLM.ndllf_evr_brl_tree == 2
        @test CLM.nbrdlf_evr_trp_tree == 4
        @test CLM.nc3_arctic_grass == 12
        @test CLM.nc4_grass == 14
        @test CLM.nc3crop == 15
        @test CLM.nc3irrig == 16
        @test CLM.ntmp_corn == 17
        @test CLM.ntrp_soybean == 77
        @test CLM.nirrig_trp_soybean == 78

        # Real parameters
        @test CLM.REINICKERP ≈ 1.6
        @test CLM.DWOOD_PARAM ≈ 2.5e5
        @test CLM.ALLOM1 ≈ 100.0
        @test CLM.ALLOM2 ≈ 40.0
        @test CLM.ALLOM3 ≈ 0.5
        @test CLM.ALLOM1S ≈ 250.0
        @test CLM.ALLOM2S ≈ 8.0
        @test CLM.ROOT_DENSITY_PARAM ≈ 0.31e6
        @test CLM.ROOT_RADIUS_PARAM ≈ 0.29e-3

        # PFT names
        @test length(CLM.PFTNAMES) == CLM.MXPFT + 1
        @test CLM.PFTNAMES[1] == "not_vegetated"          # Fortran index 0
        @test CLM.PFTNAMES[2] == "needleleaf_evergreen_temperate_tree"  # Fortran index 1
        @test CLM.PFTNAMES[15] == "c4_grass"              # Fortran index 14
        @test CLM.PFTNAMES[16] == "c3_crop"               # Fortran index 15
        @test CLM.PFTNAMES[79] == "irrigated_tropical_soybean"  # Fortran index 78
    end

    @testset "PftconType construction" begin
        p = CLM.PftconType()
        @test isempty(p.dleaf)
        @test isempty(p.is_tree)
        @test isempty(p.noveg)
        @test size(p.rhol) == (0, 0)
    end

    @testset "pftcon_allocate!" begin
        p = CLM.PftconType()
        CLM.pftcon_allocate!(p)
        n = CLM.MXPFT + 1

        # Check 1D vector sizes
        @test length(p.dleaf) == n
        @test length(p.c3psn) == n
        @test length(p.xl) == n
        @test length(p.is_tree) == n
        @test length(p.is_shrub) == n
        @test length(p.is_grass) == n
        @test length(p.noveg) == n
        @test length(p.crop) == n
        @test length(p.irrigated) == n
        @test length(p.mergetoclmpft) == n
        @test length(p.is_pft_known_to_model) == n
        @test length(p.mxmat) == n
        @test length(p.mnNHplantdate) == n

        # Check 2D matrix sizes
        @test size(p.rhol) == (n, CLM.NUMRAD)
        @test size(p.rhos) == (n, CLM.NUMRAD)
        @test size(p.taul) == (n, CLM.NUMRAD)
        @test size(p.taus) == (n, CLM.NUMRAD)
        @test size(p.rootprof_beta) == (n, CLM.NVARIANTS)
        @test size(p.lf_f) == (n, 3)
        @test size(p.fr_f) == (n, 3)
        @test size(p.repr_structure_harvfrac) == (n, CLM.NREPR)

        # Check initial values
        @test all(p.dleaf .== 0.0)
        @test all(p.is_tree .== false)
        @test all(p.is_pft_known_to_model .== false)
        @test all(p.noveg .== typemax(Int))

        # FUN parameters
        @test length(p.a_fix) == n
        @test length(p.FUN_fracfixers) == n

        # Tree structure
        @test length(p.dbh) == n
        @test length(p.wood_density) == n
        @test length(p.crit_onset_gdd_sf) == n
        @test length(p.ndays_on) == n
    end

    @testset "pftcon_init_for_testing!" begin
        p = CLM.PftconType()
        CLM.pftcon_init_for_testing!(p)
        n = CLM.MXPFT + 1

        @test length(p.dleaf) == n
        @test size(p.rhol) == (n, CLM.NUMRAD)

        # Verify we can set values after init
        p.dleaf[1] = 0.04
        @test p.dleaf[1] ≈ 0.04
        p.is_tree[2] = true
        @test p.is_tree[2] == true
    end

    @testset "pftcon_clean!" begin
        p = CLM.PftconType()
        CLM.pftcon_allocate!(p)

        # Verify arrays are allocated
        @test length(p.dleaf) == CLM.MXPFT + 1

        # Clean
        CLM.pftcon_clean!(p)

        # Verify arrays are emptied
        @test isempty(p.dleaf)
        @test isempty(p.is_tree)
        @test isempty(p.mxmat)
        @test size(p.rhol) == (0, 0)
        @test size(p.lf_f) == (0, 0)
        @test isempty(p.FUN_fracfixers)
        @test isempty(p.dbh)
    end

    @testset "set_is_pft_known_to_model!" begin
        p = CLM.PftconType()
        CLM.pftcon_allocate!(p)

        # Set mergetoclmpft: each PFT maps to itself (identity mapping)
        for i in 1:length(p.mergetoclmpft)
            p.mergetoclmpft[i] = i - 1  # Fortran 0-based value
        end

        CLM.set_is_pft_known_to_model!(p)

        # Type 0 (Julia index 1) should always be known
        @test p.is_pft_known_to_model[1] == true

        # All types should be known with identity mapping
        @test all(p.is_pft_known_to_model)

        # Test merge: types 17-78 (Fortran) all merge to type 15 (nc3crop)
        CLM.pftcon_clean!(p)
        CLM.pftcon_allocate!(p)
        for i in 1:length(p.mergetoclmpft)
            p.mergetoclmpft[i] = i - 1  # identity first
        end
        # Now merge crops to nc3crop (Fortran index 15)
        for i in 18:79  # Julia indices for Fortran 17:78
            p.mergetoclmpft[i] = 15
        end

        CLM.set_is_pft_known_to_model!(p)
        @test p.is_pft_known_to_model[1] == true   # noveg always known
        @test p.is_pft_known_to_model[16] == true   # nc3crop (Fortran 15) is known
    end

    @testset "set_num_cfts_known_to_model!" begin
        p = CLM.PftconType()
        CLM.pftcon_allocate!(p)

        # Set all CFTs as known
        p.is_pft_known_to_model .= true

        # With cft_lb=16, cft_ub=78 (Fortran 0-based), that's 63 CFTs
        CLM.set_num_cfts_known_to_model!(p; cft_lb=16, cft_ub=78)
        @test CLM.num_cfts_known_to_model == 63

        # Set none as known
        p.is_pft_known_to_model .= false
        CLM.set_num_cfts_known_to_model!(p; cft_lb=16, cft_ub=78)
        @test CLM.num_cfts_known_to_model == 0

        # Set just a few as known
        p.is_pft_known_to_model[17] = true  # Fortran index 16
        p.is_pft_known_to_model[18] = true  # Fortran index 17
        CLM.set_num_cfts_known_to_model!(p; cft_lb=16, cft_ub=78)
        @test CLM.num_cfts_known_to_model == 2
    end

    @testset "Global pftcon instance" begin
        @test CLM.pftcon isa CLM.PftconType
    end

    @testset "Custom ndecomp_pools allocation" begin
        p = CLM.PftconType()
        CLM.pftcon_allocate!(p; ndecomp_pools=7)
        @test size(p.lf_f, 2) == 7
        @test size(p.fr_f, 2) == 7
        CLM.pftcon_clean!(p)
    end
end
