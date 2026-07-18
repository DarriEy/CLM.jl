# Tests for src/infrastructure/dyn_pft_crop_file.jl — the transient natural-PFT
# (dynpft) and crop (dyncrop) land-use data readers.
#
# Each test writes a small temporary NetCDF flanduse_timeseries file with a YEAR
# axis and the relevant PCT_* variable(s) across a couple of years, builds a
# minimal subgrid (gridcell / landunit / column / patch) by hand, then runs
# init + interp for successive years and asserts the patch/column/landunit weights
# snap to the expected per-year values.

using Test
using NCDatasets

@testset "dyn_pft_crop_file" begin

    # ==================================================================
    # dynpft: PCT_NAT_PFT -> pch.wtcol for soil (istsoil) patches
    # ==================================================================
    @testset "dynpft_init + dynpft_interp!" begin
        file_years = [2000, 2001, 2002]
        natpft_size = 3                       # 3 natural PFTs (itype 0,1,2)
        ng = 1                                # single grid cell

        # PCT for each (pft, year): columns (over pft) must sum to 100.
        #   year 2000: [50, 30, 20]
        #   year 2001: [20, 50, 30]
        #   year 2002: [10, 10, 80]
        pct = Array{Float64}(undef, natpft_size, ng, length(file_years))
        pct[:, 1, 1] = [50.0, 30.0, 20.0]
        pct[:, 1, 2] = [20.0, 50.0, 30.0]
        pct[:, 1, 3] = [10.0, 10.0, 80.0]
        # expected fractional weights (PCT / 100)
        wexp = pct ./ 100.0

        mktempdir() do dir
            fn = joinpath(dir, "flanduse_pft.nc")
            NCDataset(fn, "c") do ds
                defDim(ds, "natpft", natpft_size)
                defDim(ds, "lndgrid", ng)
                defDim(ds, "time", length(file_years))
                defVar(ds, "YEAR", Int, ("time",))[:] = file_years
                # File layout: (natpft, lndgrid, time). The 2-d slice read per
                # year is (natpft, lndgrid), which flattens to data_shape
                # [natpft_size, ngridcells] (natpft FIRST, so each grid cell's
                # PFT weights are the dim-1 sum the check_sums machinery verifies).
                pv = defVar(ds, "PCT_NAT_PFT", Float64, ("natpft", "lndgrid", "time"))
                for t in 1:length(file_years)
                    pv[:, :, t] = pct[:, :, t]      # (natpft, ng)
                end
            end

            # --- consistency check: time-1 weights vs surface dataset ---
            wt_nat_patch = wexp[:, :, 1]                # (natpft, ng) fractions

            # Mismatched surface data -> error.
            bad_surf = copy(wt_nat_patch); bad_surf[1, 1] += 0.1
            @test_throws Exception CLM.dynpft_init(fn; ngridcells = ng,
                natpft_size = natpft_size, current_year = 2000,
                wt_nat_patch = bad_surf, check_dynpft_consistency = true)

            # Wrong natpft dim size -> error.
            @test_throws Exception CLM.dynpft_init(fn; ngridcells = ng,
                natpft_size = natpft_size + 1, current_year = 2000,
                check_dynpft_consistency = false)

            # --- good init (consistency passes) ---
            state = CLM.dynpft_init(fn; ngridcells = ng, natpft_size = natpft_size,
                current_year = 1999, wt_nat_patch = wt_nat_patch,
                check_dynpft_consistency = true)

            # Build a minimal subgrid: 1 gridcell, 1 soil landunit, 1 column,
            # 3 patches (one per natural PFT, itype 0,1,2 on the natpft axis).
            np = natpft_size
            bounds = CLM.BoundsType(begg = 1, endg = 1, begl = 1, endl = 1,
                                    begc = 1, endc = 1, begp = 1, endp = np,
                                    level = 1)
            lun = CLM.LandunitData{Float64}()
            lun.itype = [CLM.ISTSOIL]
            pch = CLM.PatchData{Float64}()
            pch.gridcell = fill(1, np)
            pch.landunit = fill(1, np)
            pch.column   = fill(1, np)
            pch.itype    = collect(0:(np - 1))     # natpft itype 0..np-1
            pch.wtcol    = fill(NaN, np)

            # Before the axis -> data fixed at first year (2000) -> [.5 .3 .2].
            CLM.dynpft_interp!(state, bounds, pch, lun; year = 1999)
            @test pch.wtcol ≈ [0.50, 0.30, 0.20]
            @test sum(pch.wtcol) ≈ 1.0

            # Step to 2001 -> snaps to year-2001 weights [.2 .5 .3].
            CLM.dynpft_interp!(state, bounds, pch, lun; year = 2001)
            @test pch.wtcol ≈ [0.20, 0.50, 0.30]
            @test sum(pch.wtcol) ≈ 1.0

            # Step to 2002 (last/after series) -> [.1 .1 .8].
            CLM.dynpft_interp!(state, bounds, pch, lun; year = 2002)
            @test pch.wtcol ≈ [0.10, 0.10, 0.80]
            @test sum(pch.wtcol) ≈ 1.0

            # Step back to 2000 -> re-reads first slice.
            CLM.dynpft_interp!(state, bounds, pch, lun; year = 2000)
            @test pch.wtcol ≈ [0.50, 0.30, 0.20]

            CLM.dynpft_close!(state)
        end
    end

    # ==================================================================
    # dyncrop: PCT_CROP -> lun.wtgcell ; PCT_CFT -> col.wtlunit ;
    #          FERTNITRO_CFT -> crop_inst.fertnitro_patch
    # ==================================================================
    @testset "dyncrop_init + dyncrop_interp!" begin
        file_years = [2000, 2001]
        cft_size = 2                          # 2 crop functional types
        ng = 1

        # PCT_CROP (crop landunit weight on gridcell), per year [%].
        pct_crop = Float64[ 30.0  45.0 ]      # year 2000 -> 0.30, 2001 -> 0.45
        # PCT_CFT (per-CFT weight within crop landunit), columns sum to 100.
        cft = Array{Float64}(undef, cft_size, ng, length(file_years))
        cft[:, 1, 1] = [60.0, 40.0]           # year 2000 -> [.6 .4]
        cft[:, 1, 2] = [25.0, 75.0]           # year 2001 -> [.25 .75]
        # FERTNITRO_CFT (no conversion; absolute values).
        fert = Array{Float64}(undef, cft_size, ng, length(file_years))
        fert[:, 1, 1] = [10.0, 20.0]
        fert[:, 1, 2] = [11.0, 22.0]

        mktempdir() do dir
            fn = joinpath(dir, "flanduse_crop.nc")
            NCDataset(fn, "c") do ds
                defDim(ds, "cft", cft_size)
                defDim(ds, "lndgrid", ng)
                defDim(ds, "time", length(file_years))
                defVar(ds, "YEAR", Int, ("time",))[:] = file_years
                defVar(ds, "PCT_CROP", Float64, ("lndgrid", "time"))[:, :] =
                    reshape(pct_crop, ng, length(file_years))
                cv = defVar(ds, "PCT_CFT", Float64, ("cft", "lndgrid", "time"))
                fv = defVar(ds, "FERTNITRO_CFT", Float64, ("cft", "lndgrid", "time"))
                for t in 1:length(file_years)
                    cv[:, :, t] = cft[:, :, t]      # (cft, ng)
                    fv[:, :, t] = fert[:, :, t]
                end
            end

            # Wrong cft dim size -> error.
            @test_throws Exception CLM.dyncrop_init(fn; ngridcells = ng,
                cft_size = cft_size + 1, current_year = 2000)

            state = CLM.dyncrop_init(fn; ngridcells = ng, cft_size = cft_size,
                                     current_year = 2000)

            # Set cft_lb so the patch itype -> CFT-column mapping is well defined.
            # We give the two crop patches itype = cft_lb, cft_lb+1.
            cft_lb = 15
            old_cft_lb = CLM.varpar.cft_lb
            CLM.varpar.cft_lb = cft_lb

            try
                # Minimal subgrid: 1 gridcell; landunit 1 = soil (no crop),
                # landunit 2 = crop with 2 columns, 1 patch per crop column.
                # patches: p1 -> soil (ignored), p2,p3 -> crop columns c2,c3.
                np = 3
                bounds = CLM.BoundsType(begg = 1, endg = 1, begl = 1, endl = 2,
                                        begc = 1, endc = 3, begp = 1, endp = np,
                                        level = 1)

                grc = CLM.GridcellData{Float64}()
                # landunit_indices[ltype, g]: crop (ISTCROP=2) landunit is index 2.
                grc.landunit_indices = fill(CLM.ISPVAL, CLM.MAX_LUNIT, ng)
                grc.landunit_indices[CLM.ISTSOIL, 1] = 1
                grc.landunit_indices[CLM.ISTCROP, 1] = 2

                lun = CLM.LandunitData{Float64}()
                lun.itype   = [CLM.ISTSOIL, CLM.ISTCROP]
                lun.wtgcell = fill(NaN, 2)

                col = CLM.ColumnData{Float64}()
                col.wtlunit = fill(NaN, 3)

                pch = CLM.PatchData{Float64}()
                #            p1(soil)  p2(crop c2)  p3(crop c3)
                pch.gridcell = [1, 1, 1]
                pch.landunit = [1, 2, 2]
                pch.column   = [1, 2, 3]
                pch.itype    = [0, cft_lb, cft_lb + 1]
                pch.wtcol    = fill(1.0, np)

                # --- year 2000 ---
                CLM.dyncrop_interp!(state, bounds, grc, lun, col, pch; year = 2000)
                @test lun.wtgcell[2] ≈ 0.30          # PCT_CROP / 100
                @test col.wtlunit[2] ≈ 0.60          # CFT 1 weight
                @test col.wtlunit[3] ≈ 0.40          # CFT 2 weight
                @test col.wtlunit[2] + col.wtlunit[3] ≈ 1.0

                # --- year 2001 (snaps) ---
                CLM.dyncrop_interp!(state, bounds, grc, lun, col, pch; year = 2001)
                @test lun.wtgcell[2] ≈ 0.45
                @test col.wtlunit[2] ≈ 0.25
                @test col.wtlunit[3] ≈ 0.75
                @test col.wtlunit[2] + col.wtlunit[3] ≈ 1.0

                # --- fertilizer wiring with use_crop = true ---
                crop_inst = CLM.CropData{Float64}()
                crop_inst.fertnitro_patch = fill(NaN, np)
                CLM.dyncrop_interp!(state, bounds, grc, lun, col, pch; year = 2000,
                                    use_crop = true, crop_inst = crop_inst)
                @test crop_inst.fertnitro_patch[2] ≈ 10.0
                @test crop_inst.fertnitro_patch[3] ≈ 20.0

                # --- duplicate-CFT-on-one-column guard ---
                # Two crop patches mapped to the SAME column -> error.
                pch_bad = CLM.PatchData{Float64}()
                pch_bad.gridcell = [1, 1, 1]
                pch_bad.landunit = [1, 2, 2]
                pch_bad.column   = [1, 2, 2]       # p2,p3 both on column 2
                pch_bad.itype    = [0, cft_lb, cft_lb + 1]
                pch_bad.wtcol    = fill(1.0, np)
                @test_throws Exception CLM.dyncrop_interp!(state, bounds, grc, lun,
                    col, pch_bad; year = 2000)
            finally
                CLM.varpar.cft_lb = old_cft_lb
                CLM.dyncrop_close!(state)
            end
        end
    end

    # ==================================================================
    # REAL CTSM layout: the gridcell axis stored as an (lsmlon, lsmlat) PAIR.
    #
    # Every real CTSM surface / landuse_timeseries dataset carries the
    # horizontal axis as two dimensions, not a single 'lndgrid'. Fortran reads
    # both through the same `grlnd` decomposition. CLM.jl's reader used to
    # hard-reject anything with more than 2 spatial dims ("variable has 4 dims;
    # expected 2 or 3"), which meant it could not read ANY real CTSM
    # flanduse_timeseries — the unit tests above only ever fabricated the
    # 1-D 'lndgrid' flavour, so nothing caught it. Found by generating the
    # dtrotr (tropical deforestation-fire) parity reference, which needs a real
    # flanduse file that BOTH codes read.
    # ==================================================================
    @testset "dynpft reads the (lsmlon, lsmlat) gridcell layout" begin
        file_years = [2000, 2001]
        natpft_size = 3
        nlon, nlat = 2, 3
        ng = nlon * nlat

        # (natpft, lsmlon, lsmlat, time); each gridcell's PFT weights sum to 100.
        pct = Array{Float64}(undef, natpft_size, nlon, nlat, length(file_years))
        for t in 1:length(file_years), j in 1:nlat, i in 1:nlon
            a = 10.0 * i + 5.0 * j + 20.0 * t          # varies per cell AND year
            b = 30.0
            pct[:, i, j, t] = [a, b, 100.0 - a - b]
        end

        mktempdir() do dir
            fn = joinpath(dir, "flanduse_lonlat.nc")
            NCDataset(fn, "c") do ds
                defDim(ds, "natpft", natpft_size)
                defDim(ds, "lsmlon", nlon)
                defDim(ds, "lsmlat", nlat)
                defDim(ds, "time", length(file_years))
                defVar(ds, "YEAR", Int, ("time",))[:] = file_years
                pv = defVar(ds, "PCT_NAT_PFT", Float64,
                            ("natpft", "lsmlon", "lsmlat", "time"))
                for t in 1:length(file_years)
                    pv[:, :, :, t] = pct[:, :, :, t]
                end
            end

            # The (lsmlon, lsmlat) pair must fold into the gridcell axis, giving
            # exactly the same data_shape [natpft, ng] as the 'lndgrid' flavour.
            # Assert the VALUES, not merely that init did not throw.
            expect(t) = reshape(pct[:, :, :, t], natpft_size, ng) ./ 100.0

            state = CLM.dynpft_init(fn; ngridcells = ng, natpft_size = natpft_size,
                current_year = 2000, wt_nat_patch = expect(1),
                check_dynpft_consistency = true)
            got = reshape(state.wtpatch.data_at_tlower, natpft_size, ng)
            @test got ≈ expect(1)
            # ... and it must be a real read, not a coincidence: the year-2001
            # slice differs from year-2000 everywhere.
            @test !isapprox(expect(2), expect(1))

            state2 = CLM.dynpft_init(fn; ngridcells = ng, natpft_size = natpft_size,
                current_year = 2001, check_dynpft_consistency = false)
            @test reshape(state2.wtpatch.data_at_tlower, natpft_size, ng) ≈ expect(2)

            # A genuinely wrong element count must still be rejected, not
            # silently reshaped.
            @test_throws Exception CLM.dynpft_init(fn; ngridcells = ng + 1,
                natpft_size = natpft_size, current_year = 2000,
                check_dynpft_consistency = false)
        end
    end

end
