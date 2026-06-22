# Tests for src/infrastructure/dyn_lake_urban_file.jl — the transient lake &
# urban land-use readers (dynlake_init/interp, dynurban_init/interp).
#
# Each test writes a small temporary NetCDF file (mirroring test_dyn_file_io.jl)
# with a YEAR axis plus a PCT_LAKE (space x time) and/or PCT_URBAN
# (space x numurbl x time) variable, builds a minimal gridcell/landunit subgrid
# carrying the lake (ISTDLAK) and urban (ISTURB_TBD/HD/MD) landunits, then runs
# init + interp over successive years and asserts the landunit weights snap to
# the correct per-year values (percent -> fraction via conversion_factor=100).

using Test
using NCDatasets

@testset "dyn_lake_urban_file" begin

    # Build a minimal subgrid with `ngrid` gridcells, each carrying one lake
    # landunit (ISTDLAK) and the three urban landunits (TBD/HD/MD). Returns
    # (grc, lun) with landunit_indices and wtgcell populated.
    function build_lake_urban_subgrid(ngrid::Int)
        # Per gridcell: 4 landunits (lake, tbd, hd, md).
        ltypes = [CLM.ISTDLAK, CLM.ISTURB_TBD, CLM.ISTURB_HD, CLM.ISTURB_MD]
        nlun = ngrid * length(ltypes)

        grc = CLM.GridcellData{Float64}()
        grc.landunit_indices = fill(CLM.ISPVAL, CLM.MAX_LUNIT, ngrid)

        lun = CLM.LandunitData{Float64}()
        lun.gridcell = zeros(Int, nlun)
        lun.itype    = zeros(Int, nlun)
        lun.wtgcell  = fill(NaN, nlun)

        l = 0
        for g in 1:ngrid
            for lt in ltypes
                l += 1
                lun.gridcell[l] = g
                lun.itype[l] = lt
                lun.wtgcell[l] = 0.0
                grc.landunit_indices[lt, g] = l
            end
        end
        return grc, lun
    end

    # ==================================================================
    # set_landunit_weight! — basic set + non-existent-landunit guard
    # ==================================================================
    @testset "set_landunit_weight!" begin
        grc, lun = build_lake_urban_subgrid(1)
        l_lake = grc.landunit_indices[CLM.ISTDLAK, 1]

        CLM.set_landunit_weight!(grc, lun, 1, CLM.ISTDLAK, 0.3)
        @test lun.wtgcell[l_lake] == 0.3

        # No glacier (ISTICE) landunit on this gridcell -> weight 0 is a no-op.
        CLM.set_landunit_weight!(grc, lun, 1, CLM.ISTICE, 0.0)
        # ... but a non-zero weight to a non-existent landunit is an error.
        @test_throws Exception CLM.set_landunit_weight!(grc, lun, 1, CLM.ISTICE, 0.1)
    end

    # ==================================================================
    # dynlake_init / dynlake_interp (PCT_LAKE)
    # ==================================================================
    @testset "dynlake" begin
        file_years = [2000, 2001, 2002]
        ngrid = 2
        # PCT_LAKE[space, time] in PERCENT. Distinct per (g, year).
        pct_lake = Float64[10 * g + (file_years[t] - 2000) * 5
                           for g in 1:ngrid, t in 1:length(file_years)]

        mktempdir() do dir
            fn = joinpath(dir, "dynlake.nc")
            NCDataset(fn, "c") do ds
                defDim(ds, "lndgrid", ngrid)
                defDim(ds, "time", length(file_years))
                defVar(ds, "YEAR", Int, ("time",))[:] = file_years
                defVar(ds, "PCT_LAKE", Float64, ("lndgrid", "time"))[:, :] = pct_lake
            end

            grc, lun = build_lake_urban_subgrid(ngrid)
            l_lake = [grc.landunit_indices[CLM.ISTDLAK, g] for g in 1:ngrid]

            # Init positioned before the axis (1999) -> first-year slice.
            dl = CLM.dynlake_init(fn, 1, ngrid; current_year = 1999)
            @test dl.dynlake_file.time_info.years == file_years

            # interp at 2000 -> fraction = percent/100 of first-year slice.
            CLM.dynlake_interp(dl, grc, lun, 1, ngrid, 2000)
            for g in 1:ngrid
                @test lun.wtgcell[l_lake[g]] ≈ pct_lake[g, 1] / 100.0
            end

            # interp at 2001 -> second-year slice.
            CLM.dynlake_interp(dl, grc, lun, 1, ngrid, 2001)
            for g in 1:ngrid
                @test lun.wtgcell[l_lake[g]] ≈ pct_lake[g, 2] / 100.0
            end

            # interp at 2002 (last year) -> third-year slice (after-series).
            CLM.dynlake_interp(dl, grc, lun, 1, ngrid, 2002)
            for g in 1:ngrid
                @test lun.wtgcell[l_lake[g]] ≈ pct_lake[g, 3] / 100.0
            end

            # interp past the axis (2050) -> still last-year slice.
            CLM.dynlake_interp(dl, grc, lun, 1, ngrid, 2050)
            for g in 1:ngrid
                @test lun.wtgcell[l_lake[g]] ≈ pct_lake[g, 3] / 100.0
            end

            # urban landunits untouched by the lake interp.
            for g in 1:ngrid
                @test lun.wtgcell[grc.landunit_indices[CLM.ISTURB_TBD, g]] == 0.0
            end

            CLM.dyn_file_close!(dl.dynlake_file)
        end
    end

    # ==================================================================
    # dynurban_init / dynurban_interp (PCT_URBAN[numurbl])
    # ==================================================================
    @testset "dynurban" begin
        file_years = [2000, 2001]
        ngrid = 2
        numurbl = CLM.NUMURBL   # == 3 (TBD/HD/MD)
        # PCT_URBAN[space, numurbl, time] in PERCENT, distinct per (g, class, year).
        pct_urban = Array{Float64}(undef, ngrid, numurbl, length(file_years))
        for t in 1:length(file_years), n in 1:numurbl, g in 1:ngrid
            pct_urban[g, n, t] = g + 10 * n + 100 * (t - 1)   # 0..100ish percents
        end

        mktempdir() do dir
            fn = joinpath(dir, "dynurban.nc")
            NCDataset(fn, "c") do ds
                defDim(ds, "lndgrid", ngrid)
                defDim(ds, "numurbl", numurbl)
                defDim(ds, "time", length(file_years))
                defVar(ds, "YEAR", Int, ("time",))[:] = file_years
                defVar(ds, "PCT_URBAN", Float64,
                       ("lndgrid", "numurbl", "time"))[:, :, :] = pct_urban
            end

            grc, lun = build_lake_urban_subgrid(ngrid)
            l_tbd = [grc.landunit_indices[CLM.ISTURB_TBD, g] for g in 1:ngrid]
            l_hd  = [grc.landunit_indices[CLM.ISTURB_HD,  g] for g in 1:ngrid]
            l_md  = [grc.landunit_indices[CLM.ISTURB_MD,  g] for g in 1:ngrid]

            du = CLM.dynurban_init(fn, 1, ngrid; current_year = 2000)
            @test du.dynurban_file.time_info.years == file_years
            @test du.wturban.data_shape == [ngrid, numurbl]

            # interp at 2000 -> first-year slice, per density class.
            CLM.dynurban_interp(du, grc, lun, 1, ngrid, 2000)
            for g in 1:ngrid
                @test lun.wtgcell[l_tbd[g]] ≈ pct_urban[g, 1, 1] / 100.0
                @test lun.wtgcell[l_hd[g]]  ≈ pct_urban[g, 2, 1] / 100.0
                @test lun.wtgcell[l_md[g]]  ≈ pct_urban[g, 3, 1] / 100.0
            end

            # interp at 2001 -> second-year slice.
            CLM.dynurban_interp(du, grc, lun, 1, ngrid, 2001)
            for g in 1:ngrid
                @test lun.wtgcell[l_tbd[g]] ≈ pct_urban[g, 1, 2] / 100.0
                @test lun.wtgcell[l_hd[g]]  ≈ pct_urban[g, 2, 2] / 100.0
                @test lun.wtgcell[l_md[g]]  ≈ pct_urban[g, 3, 2] / 100.0
            end

            # lake landunit untouched by the urban interp.
            for g in 1:ngrid
                @test lun.wtgcell[grc.landunit_indices[CLM.ISTDLAK, g]] == 0.0
            end

            CLM.dyn_file_close!(du.dynurban_file)
        end
    end

    # ==================================================================
    # dynurban_init numurbl dimension-size check
    # ==================================================================
    @testset "dynurban numurbl mismatch errors" begin
        file_years = [2000]
        ngrid = 1
        bad_numurbl = CLM.NUMURBL + 1

        mktempdir() do dir
            fn = joinpath(dir, "dynurban_bad.nc")
            NCDataset(fn, "c") do ds
                defDim(ds, "lndgrid", ngrid)
                defDim(ds, "numurbl", bad_numurbl)
                defDim(ds, "time", length(file_years))
                defVar(ds, "YEAR", Int, ("time",))[:] = file_years
                defVar(ds, "PCT_URBAN", Float64,
                       ("lndgrid", "numurbl", "time"))[:, :, :] =
                    fill(0.0, ngrid, bad_numurbl, length(file_years))
            end

            @test_throws Exception CLM.dynurban_init(fn, 1, ngrid;
                                                     current_year = 2000)
        end
    end

end
