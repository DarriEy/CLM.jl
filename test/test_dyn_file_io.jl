# Tests for src/infrastructure/dyn_file_io.jl — the transient (dynamic) land-use
# file time-stepping framework: DynTimeInfo, DynFile, DynVarTimeUninterp.
#
# DynTimeInfo tests are file-free: a synthetic `years` vector exercises the
# year-bounding / index / interp-weight behavior at: before the axis, on the first
# year, mid-axis, on a boundary year, and after the axis.
#
# DynFile + DynVarTimeUninterp tests write a small temporary NetCDF file (mirroring
# test/test_history_io.jl) with a YEAR axis and a data variable, then assert that
# stepping the current year reads the correct data_at_tlower slice.

using Test
using NCDatasets

@testset "dyn_file_io" begin

    # ==================================================================
    # DynTimeInfo: year-bounding, indices, current-year accessors
    # ==================================================================
    @testset "DynTimeInfo year-bounding" begin
        years = [2000, 2001, 2002, 2003]   # nyears = 4
        n = length(years)

        # --- before the axis: year < years[1] -> lower=upper=1 ---
        ti = CLM.dyn_time_info(years, CLM.YEAR_POSITION_START_OF_TIMESTEP;
                               current_year = 1995)
        @test CLM.get_time_index_lower(ti) == 1
        @test CLM.get_time_index_upper(ti) == 1
        @test CLM.is_before_time_series(ti)
        @test !CLM.is_after_time_series(ti)
        @test !CLM.is_within_bounds(ti)
        @test CLM.get_current_year(ti) == 2000   # clamped to first year

        # --- on the first year: years[1] -> lower=1, upper=2 (within bounds) ---
        CLM.set_current_year!(ti, 2000)
        @test CLM.get_time_index_lower(ti) == 1
        @test CLM.get_time_index_upper(ti) == 2
        @test !CLM.is_before_time_series(ti)
        @test !CLM.is_after_time_series(ti)
        @test CLM.is_within_bounds(ti)
        @test CLM.get_current_year(ti) == 2000

        # --- mid-axis: years[2] = 2001 -> lower=2, upper=3 ---
        CLM.set_current_year!(ti, 2001)
        @test CLM.get_time_index_lower(ti) == 2
        @test CLM.get_time_index_upper(ti) == 3
        @test CLM.is_within_bounds(ti)
        @test CLM.get_current_year(ti) == 2001
        @test CLM.get_year(ti, 2) == 2001
        @test CLM.get_year(ti, 3) == 2002

        # --- another mid-axis: years[3] = 2002 -> lower=3, upper=4 ---
        CLM.set_current_year!(ti, 2002)
        @test CLM.get_time_index_lower(ti) == 3
        @test CLM.get_time_index_upper(ti) == 4
        @test CLM.is_within_bounds(ti)

        # --- on the last (boundary) year: year >= years[nyears] -> lower=upper=nyears
        #     "after the time series" is TRUE when current year == last file year.
        CLM.set_current_year!(ti, 2003)
        @test CLM.get_time_index_lower(ti) == n
        @test CLM.get_time_index_upper(ti) == n
        @test !CLM.is_before_time_series(ti)
        @test CLM.is_after_time_series(ti)
        @test !CLM.is_within_bounds(ti)
        @test CLM.get_current_year(ti) == 2003

        # --- after the axis: year > years[nyears] -> still lower=upper=nyears ---
        CLM.set_current_year!(ti, 2010)
        @test CLM.get_time_index_lower(ti) == n
        @test CLM.get_time_index_upper(ti) == n
        @test CLM.is_after_time_series(ti)
        @test CLM.get_current_year(ti) == 2003

        # --- step back into the middle works (indices recompute) ---
        CLM.set_current_year!(ti, 2001)
        @test CLM.get_time_index_lower(ti) == 2
        @test CLM.get_time_index_upper(ti) == 3
        @test CLM.is_within_bounds(ti)

        # get_year bounds check
        @test_throws Exception CLM.get_year(ti, 0)
        @test_throws Exception CLM.get_year(ti, n + 1)
    end

    @testset "DynTimeInfo END_OF_TIMESTEP position stored" begin
        ti = CLM.dyn_time_info([2000, 2001], CLM.YEAR_POSITION_END_OF_TIMESTEP;
                               current_year = 2000)
        @test ti.year_position == CLM.YEAR_POSITION_END_OF_TIMESTEP
        @test CLM.get_time_index_lower(ti) == 1
        @test CLM.get_time_index_upper(ti) == 2
    end

    @testset "DynTimeInfo no current_year leaves constructor default" begin
        # Without current_year, indices stay at the arbitrary constructor default.
        ti = CLM.dyn_time_info([2000, 2001, 2002], CLM.YEAR_POSITION_START_OF_TIMESTEP)
        @test CLM.get_time_index_lower(ti) == 1
        @test CLM.get_time_index_upper(ti) == 1
        @test CLM.is_before_time_series(ti)
    end

    # ==================================================================
    # DynFile + DynVarTimeUninterp: write a temp NetCDF, step the year,
    # and assert the correct data slice is read (step-change at year bounds).
    # ==================================================================
    @testset "DynFile + DynVarTimeUninterp 1d" begin
        file_years = [2000, 2001, 2002]
        nspace = 4
        # data[s, t] = 100*s + year   (distinct, checkable per (space, year))
        data1d = Float64[100 * s + file_years[t] for s in 1:nspace, t in 1:length(file_years)]

        mktempdir() do dir
            fn = joinpath(dir, "dyn_test_1d.nc")
            NCDataset(fn, "c") do ds
                defDim(ds, "space", nspace)
                defDim(ds, "time", length(file_years))
                yv = defVar(ds, "YEAR", Int, ("time",))
                yv[:] = file_years
                dv = defVar(ds, "MYVAR", Float64, ("space", "time"))
                dv[:, :] = data1d
            end

            # Open file positioned before the axis (1999) -> lower index 1.
            df = CLM.dyn_file_open(fn, CLM.YEAR_POSITION_START_OF_TIMESTEP;
                                   current_year = 1999)
            @test df.time_info.years == file_years
            @test CLM.get_time_index_lower(df.time_info) == 1
            @test CLM.is_before_time_series(df.time_info)

            var = CLM.dyn_var_time_uninterp(df, "MYVAR", "space", 1.0, false,
                                            [nspace])
            @test var.data_on_file
            # Before the axis -> data fixed at first year (2000).
            @test CLM.get_current_data_1d(var) == data1d[:, 1]
            @test var.time_index_lower == 1

            # Step to 2000 (first year) -> lower still 1 -> same data, no reread change.
            CLM.set_current_year!(df, 2000)
            @test CLM.get_current_data_1d(var) == data1d[:, 1]
            @test var.time_index_lower == 1

            # Step to 2001 -> lower index 2 -> data snaps to year 2001 slice.
            CLM.set_current_year!(df, 2001)
            @test CLM.get_current_data_1d(var) == data1d[:, 2]
            @test var.time_index_lower == 2

            # Step to 2002 (last year) -> after series -> lower = nyears = 3.
            CLM.set_current_year!(df, 2002)
            @test CLM.is_after_time_series(df.time_info)
            @test CLM.get_current_data_1d(var) == data1d[:, 3]
            @test var.time_index_lower == 3

            # Step past the axis (2050) -> still last-year slice.
            CLM.set_current_year!(df, 2050)
            @test CLM.get_current_data_1d(var) == data1d[:, 3]
            @test var.time_index_lower == 3

            # Step back to 2001 -> reread the 2001 slice.
            CLM.set_current_year!(df, 2001)
            @test CLM.get_current_data_1d(var) == data1d[:, 2]
            @test var.time_index_lower == 2

            CLM.dyn_file_close!(df)
            @test df.ds === nothing
        end
    end

    @testset "DynVarTimeUninterp conversion_factor" begin
        file_years = [2000, 2001]
        nspace = 3
        data1d = Float64[10 * s + t for s in 1:nspace, t in 1:length(file_years)]

        mktempdir() do dir
            fn = joinpath(dir, "dyn_test_conv.nc")
            NCDataset(fn, "c") do ds
                defDim(ds, "space", nspace)
                defDim(ds, "time", length(file_years))
                defVar(ds, "YEAR", Int, ("time",))[:] = file_years
                defVar(ds, "MYVAR", Float64, ("space", "time"))[:, :] = data1d
            end

            df = CLM.dyn_file_open(fn, CLM.YEAR_POSITION_START_OF_TIMESTEP;
                                   current_year = 2000)
            # conversion_factor = 2.0 -> data DIVIDED by 2.
            var = CLM.dyn_var_time_uninterp(df, "MYVAR", "space", 2.0, false, [nspace])
            @test CLM.get_current_data_1d(var) ≈ data1d[:, 1] ./ 2.0
            CLM.dyn_file_close!(df)
        end
    end

    @testset "DynVarTimeUninterp 2d + check_sums_equal_1" begin
        file_years = [2000, 2001]
        d1, d2 = 2, 3   # data_shape = [d1, d2]; sums over dim 1 must equal 1
        # Build per-year (d1 x d2) arrays whose columns sum to 1, flattened with time last.
        # File layout: (d1, d2, time).
        arr = Array{Float64}(undef, d1, d2, length(file_years))
        for t in 1:length(file_years)
            for j in 1:d2
                # column (over d1) sums to 1
                arr[1, j, t] = 0.25 + 0.05 * t
                arr[2, j, t] = 1.0 - arr[1, j, t]
            end
        end

        mktempdir() do dir
            fn = joinpath(dir, "dyn_test_2d.nc")
            NCDataset(fn, "c") do ds
                defDim(ds, "d1", d1)
                defDim(ds, "d2", d2)
                defDim(ds, "time", length(file_years))
                defVar(ds, "YEAR", Int, ("time",))[:] = file_years
                defVar(ds, "PCT", Float64, ("d1", "d2", "time"))[:, :, :] = arr
            end

            df = CLM.dyn_file_open(fn, CLM.YEAR_POSITION_START_OF_TIMESTEP;
                                   current_year = 2000)
            var = CLM.dyn_var_time_uninterp(df, "PCT", "d1", 1.0, true, [d1, d2])
            got = CLM.get_current_data_2d(var)
            @test size(got) == (d1, d2)
            @test got ≈ arr[:, :, 1]
            # Columns sum to 1 (check_sums_equal_1 passed during read).
            @test all(abs.(sum(got, dims = 1) .- 1.0) .< 1e-12)

            # Step into year 2001 -> snaps to 2001 slice.
            CLM.set_current_year!(df, 2001)
            @test CLM.get_current_data_2d(var) ≈ arr[:, :, 2]
            CLM.dyn_file_close!(df)
        end
    end

    @testset "DynVarTimeUninterp allow_nodata + missing var error" begin
        mktempdir() do dir
            fn = joinpath(dir, "dyn_test_nodata.nc")
            NCDataset(fn, "c") do ds
                defDim(ds, "space", 2)
                defDim(ds, "time", 2)
                defVar(ds, "YEAR", Int, ("time",))[:] = [2000, 2001]
            end

            df = CLM.dyn_file_open(fn, CLM.YEAR_POSITION_START_OF_TIMESTEP;
                                   current_year = 2000)
            # Missing var with allow_nodata=false -> error.
            @test_throws Exception CLM.dyn_var_time_uninterp(df, "ABSENT", "space",
                                                             1.0, false, [2])
            # Missing var with allow_nodata=true -> data set to zero.
            var = CLM.dyn_var_time_uninterp(df, "ABSENT", "space", 1.0, false, [2];
                                            allow_nodata = true)
            @test !var.data_on_file
            @test CLM.get_current_data_1d(var) == zeros(2)
            CLM.dyn_file_close!(df)
        end
    end

end
