@testset "LakeStateData" begin

    # Ensure varpar is initialized for dimension parameters
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)

    nlevlak = CLM.varpar.nlevlak

    @testset "default construction" begin
        ls = CLM.LakeStateData()
        @test length(ls.lakefetch_col) == 0
        @test length(ls.etal_col) == 0
        @test length(ls.lake_raw_col) == 0
        @test length(ls.ks_col) == 0
        @test length(ls.ws_col) == 0
        @test length(ls.ust_lake_col) == 0
        @test length(ls.betaprime_col) == 0
        @test length(ls.savedtke1_col) == 0
        @test size(ls.lake_icefrac_col) == (0, 0)
        @test length(ls.lake_icefracsurf_col) == 0
        @test length(ls.lake_icethick_col) == 0
        @test length(ls.lakeresist_col) == 0
        @test length(ls.ram1_lake_patch) == 0
    end

    @testset "lakestate_init!" begin
        nc = 8
        np = 12
        ls = CLM.LakeStateData()
        CLM.lakestate_init!(ls, nc, np)

        # --- Check 1D column sizes ---
        @test length(ls.lakefetch_col) == nc
        @test length(ls.etal_col) == nc
        @test length(ls.lake_raw_col) == nc
        @test length(ls.ks_col) == nc
        @test length(ls.ws_col) == nc
        @test length(ls.ust_lake_col) == nc
        @test length(ls.betaprime_col) == nc
        @test length(ls.savedtke1_col) == nc
        @test length(ls.lake_icefracsurf_col) == nc
        @test length(ls.lake_icethick_col) == nc
        @test length(ls.lakeresist_col) == nc

        # --- Check 2D size ---
        @test size(ls.lake_icefrac_col) == (nc, nlevlak)

        # --- Check patch-level size ---
        @test length(ls.ram1_lake_patch) == np

        # --- Check NaN initialization ---
        @test all(isnan, ls.lakefetch_col)
        @test all(isnan, ls.etal_col)
        @test all(isnan, ls.lake_raw_col)
        @test all(isnan, ls.ks_col)
        @test all(isnan, ls.ws_col)
        @test all(isnan, ls.betaprime_col)
        @test all(isnan, ls.lake_icefrac_col)
        @test all(isnan, ls.lake_icefracsurf_col)
        @test all(isnan, ls.lake_icethick_col)
        @test all(isnan, ls.lakeresist_col)
        @test all(isnan, ls.ram1_lake_patch)

        # --- Check SPVAL initialization ---
        @test all(x -> x == CLM.SPVAL, ls.savedtke1_col)
        @test all(x -> x == CLM.SPVAL, ls.ust_lake_col)
    end

    @testset "lakestate_clean!" begin
        ls = CLM.LakeStateData()
        CLM.lakestate_init!(ls, 5, 8)
        CLM.lakestate_clean!(ls)

        @test length(ls.lakefetch_col) == 0
        @test length(ls.etal_col) == 0
        @test length(ls.lake_raw_col) == 0
        @test length(ls.ks_col) == 0
        @test length(ls.ws_col) == 0
        @test length(ls.ust_lake_col) == 0
        @test length(ls.betaprime_col) == 0
        @test length(ls.savedtke1_col) == 0
        @test size(ls.lake_icefrac_col) == (0, 0)
        @test length(ls.lake_icefracsurf_col) == 0
        @test length(ls.lake_icethick_col) == 0
        @test length(ls.lakeresist_col) == 0
        @test length(ls.ram1_lake_patch) == 0
    end

    @testset "lakestate_init_cold! (all columns)" begin
        nc = 5
        np = 5
        ls = CLM.LakeStateData()
        CLM.lakestate_init!(ls, nc, np)

        CLM.lakestate_init_cold!(ls, 1:nc)

        for c in 1:nc
            # Ice fraction should be zero for all lake levels
            for j in 1:nlevlak
                @test ls.lake_icefrac_col[c, j] ≈ 0.0
            end
            # savedtke1 should be TKWAT
            @test ls.savedtke1_col[c] ≈ CLM.TKWAT
            # ust_lake should be 0.1
            @test ls.ust_lake_col[c] ≈ 0.1
        end
    end

    @testset "lakestate_init_cold! (with mask)" begin
        nc = 6
        np = 6
        ls = CLM.LakeStateData()
        CLM.lakestate_init!(ls, nc, np)

        # Only columns 2, 4, 6 are lake
        mask_lake = BitVector([false, true, false, true, false, true])
        CLM.lakestate_init_cold!(ls, 1:nc; mask_lake=mask_lake)

        # Lake columns should be initialized
        for c in [2, 4, 6]
            for j in 1:nlevlak
                @test ls.lake_icefrac_col[c, j] ≈ 0.0
            end
            @test ls.savedtke1_col[c] ≈ CLM.TKWAT
            @test ls.ust_lake_col[c] ≈ 0.1
        end

        # Non-lake columns should retain SPVAL/NaN from init
        for c in [1, 3, 5]
            @test all(isnan, ls.lake_icefrac_col[c, :])
            @test ls.savedtke1_col[c] == CLM.SPVAL
            @test ls.ust_lake_col[c] == CLM.SPVAL
        end
    end

    @testset "field mutability" begin
        ls = CLM.LakeStateData()
        CLM.lakestate_init!(ls, 3, 4)

        # Write and verify 1D column fields
        ls.lakefetch_col[1] = 500.0
        @test ls.lakefetch_col[1] == 500.0

        ls.etal_col[2] = 0.5
        @test ls.etal_col[2] == 0.5

        ls.savedtke1_col[3] = 1.23
        @test ls.savedtke1_col[3] == 1.23

        ls.lakeresist_col[1] = 42.0
        @test ls.lakeresist_col[1] == 42.0

        # Write and verify 2D field
        ls.lake_icefrac_col[2, 3] = 0.75
        @test ls.lake_icefrac_col[2, 3] == 0.75

        # Write and verify patch-level field
        ls.ram1_lake_patch[4] = 99.0
        @test ls.ram1_lake_patch[4] == 99.0
    end

    @testset "re-init overwrites previous state" begin
        ls = CLM.LakeStateData()
        CLM.lakestate_init!(ls, 3, 4)
        ls.lakefetch_col[1] = 999.0

        CLM.lakestate_init!(ls, 7, 10)
        @test length(ls.lakefetch_col) == 7
        @test all(isnan, ls.lakefetch_col)
        @test size(ls.lake_icefrac_col) == (7, nlevlak)
        @test length(ls.ram1_lake_patch) == 10
    end

    @testset "stub functions run without error" begin
        ls = CLM.LakeStateData()
        CLM.lakestate_init!(ls, 5, 8)

        @test CLM.lakestate_init_history!(ls, 1:5, 1:8) === nothing
        @test CLM.lakestate_restart!(ls, 1:5) === nothing
    end

    @testset "lakestate_init_history! sets SPVAL" begin
        ls = CLM.LakeStateData()
        CLM.lakestate_init!(ls, 4, 6)

        CLM.lakestate_init_history!(ls, 1:4, 1:6)

        for c in 1:4
            @test all(x -> x == CLM.SPVAL, ls.lake_icefrac_col[c, :])
            @test ls.lake_icefracsurf_col[c] == CLM.SPVAL
            @test ls.lake_icethick_col[c] == CLM.SPVAL
            @test ls.savedtke1_col[c] == CLM.SPVAL
            @test ls.ust_lake_col[c] == CLM.SPVAL
        end
        for p in 1:6
            @test ls.ram1_lake_patch[p] == CLM.SPVAL
        end
    end

end
