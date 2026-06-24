@testset "Multi-gridcell surfdata (2D read + init_gridcells)" begin
    using NCDatasets

    # ----------------------------------------------------------------------
    # Build a small synthetic surface dataset on a 2 (lon) x 3 (lat) grid, so
    # ng = 6 gridcells, each with DISTINCT lat/lon + landunit/PFT/soil weights.
    # Spatial dims follow the real CLM surface file: (lsmlon, lsmlat) for
    # scalars, (natpft|nlevsoi, lsmlon, lsmlat) for layered fields.
    # ----------------------------------------------------------------------
    # Ensure varpar/varcon dimensions are initialized (self-contained: does not
    # depend on test ordering). Matches the clm_initialize! setup.
    CLM.varctl.soil_layerstruct_predefined = "20SL_8.5m"
    CLM.varpar_init!(CLM.varpar, 1, 15, 2, 5)
    CLM.varcon_init!()

    nlon = 2
    nlat = 3
    ng   = nlon * nlat
    natpft_size = CLM.varpar.natpft_ub - CLM.varpar.natpft_lb + 1
    nlevsoi = CLM.varpar.nlevsoi

    # Distinct per-cell value: flat index g = i + (j-1)*nlon (lon-fastest).
    gflat(i, j) = i + (j - 1) * nlon

    tmpdir = mktempdir()
    fsurdat = joinpath(tmpdir, "surfdata_synthetic_2x3.nc")

    NCDataset(fsurdat, "c") do ds
        defDim(ds, "lsmlon", nlon)
        defDim(ds, "lsmlat", nlat)
        defDim(ds, "natpft", natpft_size)
        defDim(ds, "nlevsoi", nlevsoi)

        # Per-cell lat/lon: lat varies with j, lon varies with i (distinct cells).
        latixy = [10.0 * j + i for i in 1:nlon, j in 1:nlat]   # (lon,lat)
        longxy = [100.0 * i + j for i in 1:nlon, j in 1:nlat]
        defVar(ds, "LATIXY", latixy, ("lsmlon", "lsmlat"))
        defVar(ds, "LONGXY", longxy, ("lsmlon", "lsmlat"))

        # Special landunits: give each cell a distinct lake fraction.
        pct_lake = [Float64(gflat(i, j)) for i in 1:nlon, j in 1:nlat]  # 1..6 %
        defVar(ds, "PCT_LAKE", pct_lake, ("lsmlon", "lsmlat"))
        defVar(ds, "PCT_WETLAND", zeros(nlon, nlat), ("lsmlon", "lsmlat"))
        defVar(ds, "PCT_GLACIER", zeros(nlon, nlat), ("lsmlon", "lsmlat"))

        # Natural-veg landunit fraction (the remainder goes to soil). Make it
        # distinct: 100 - pct_lake so wt_lunit[ISTSOIL] differs per cell.
        pct_natveg = [100.0 - Float64(gflat(i, j)) for i in 1:nlon, j in 1:nlat]
        defVar(ds, "PCT_NATVEG", pct_natveg, ("lsmlon", "lsmlat"))
        defVar(ds, "PCT_CROP", zeros(nlon, nlat), ("lsmlon", "lsmlat"))

        # PCT_NAT_PFT (natpft, lsmlon, lsmlat): put all weight on a DISTINCT PFT
        # per gridcell so the normalized wt_nat_patch differs cell-to-cell.
        pct_nat_pft = zeros(natpft_size, nlon, nlat)
        for i in 1:nlon, j in 1:nlat
            m = ((gflat(i, j) - 1) % natpft_size) + 1
            pct_nat_pft[m, i, j] = 100.0
        end
        defVar(ds, "PCT_NAT_PFT", pct_nat_pft, ("natpft", "lsmlon", "lsmlat"))

        # Soil: PCT_SAND distinct per cell (constant across depth, varies by cell).
        pct_sand = zeros(nlevsoi, nlon, nlat)
        pct_clay = zeros(nlevsoi, nlon, nlat)
        for i in 1:nlon, j in 1:nlat
            pct_sand[:, i, j] .= 10.0 * gflat(i, j)
            pct_clay[:, i, j] .= 5.0 * gflat(i, j)
        end
        defVar(ds, "PCT_SAND", pct_sand, ("nlevsoi", "lsmlon", "lsmlat"))
        defVar(ds, "PCT_CLAY", pct_clay, ("nlevsoi", "lsmlon", "lsmlat"))

        # Topo: distinct lake depth + slope per cell.
        defVar(ds, "LAKEDEPTH",
               [50.0 + gflat(i, j) for i in 1:nlon, j in 1:nlat],
               ("lsmlon", "lsmlat"))
        defVar(ds, "SLOPE",
               [0.5 + 0.1 * gflat(i, j) for i in 1:nlon, j in 1:nlat],
               ("lsmlon", "lsmlat"))
    end

    # ----------------------------------------------------------------------
    # Helper: allocate a fresh g/l/c/p set big enough for (nl, nc, np).
    # ----------------------------------------------------------------------
    function build_subgrid(surf, ng)
        (nl, nc, np) = CLM.count_subgrid_elements(surf, ng)
        bounds = CLM.BoundsType(begg=1, endg=ng, begl=1, endl=nl,
                                begc=1, endc=nc, begp=1, endp=np,
                                level=CLM.BOUNDS_LEVEL_CLUMP, clump_index=1)
        grc = CLM.GridcellData();  CLM.gridcell_init!(grc, ng)
        lun = CLM.LandunitData();  CLM.landunit_init!(lun, nl)
        col = CLM.ColumnData();    CLM.column_init!(col, nc)
        pch = CLM.PatchData();     CLM.patch_init!(pch, np)
        CLM.initGridCells!(bounds, surf, grc, lun, col, pch)
        return (bounds, grc, lun, col, pch, nl, nc, np)
    end

    @testset "surfrd_get_grid_dims" begin
        (rnlon, rnlat) = CLM.surfrd_get_grid_dims(fsurdat)
        @test rnlon == nlon
        @test rnlat == nlat
    end

    @testset "2D read: g→(i,j) map + distinct per-cell data" begin
        surf = CLM.SurfaceInputData()
        CLM.surfrd_get_data!(surf, 1, ng, fsurdat; full_grid=true)

        @test surf.nlon == nlon
        @test surf.nlat == nlat
        @test length(surf.ixy) == ng
        @test length(surf.jxy) == ng
        @test length(surf.grid_latdeg) == ng
        @test length(surf.grid_londeg) == ng

        # g → (i,j) map: lon-fastest flatten. ixy[g]=lon idx, jxy[g]=lat idx.
        for i in 1:nlon, j in 1:nlat
            g = gflat(i, j)
            @test surf.ixy[g] == i
            @test surf.jxy[g] == j
            # Per-gridcell lat/lon match the file values at (i,j).
            @test surf.grid_latdeg[g] ≈ 10.0 * j + i
            @test surf.grid_londeg[g] ≈ 100.0 * i + j
        end

        # Distinct surfdata per cell: lake weight, soil weight, sand, PFT.
        lake = [surf.wt_lunit[g, CLM.ISTDLAK] for g in 1:ng]
        @test length(unique(round.(lake; digits=6))) == ng   # all distinct
        for i in 1:nlon, j in 1:nlat
            g = gflat(i, j)
            @test surf.wt_lunit[g, CLM.ISTDLAK] ≈ Float64(g) / 100.0
            @test surf.pct_sand[g, 1] ≈ 10.0 * g
            @test surf.pct_clay[g, 1] ≈ 5.0 * g
            @test surf.lakedepth[g] ≈ 50.0 + g
            # The dominant natural PFT (highest weight) is the per-cell choice.
            m = ((g - 1) % natpft_size) + 1
            @test argmax(@view surf.wt_nat_patch[g, :]) == m
        end
    end

    @testset "initGridCells!: per-cell lat/lon + count add-up" begin
        surf = CLM.SurfaceInputData()
        CLM.surfrd_get_data!(surf, 1, ng, fsurdat; full_grid=true)
        (bounds, grc, lun, col, pch, nl, nc, np) = build_subgrid(surf, ng)

        # Each gridcell carries its OWN lat/lon (from surf.grid_latdeg/londeg).
        for i in 1:nlon, j in 1:nlat
            g = gflat(i, j)
            @test grc.latdeg[g] ≈ 10.0 * j + i
            @test grc.londeg[g] ≈ 100.0 * i + j
            @test grc.lat[g] ≈ (10.0 * j + i) * CLM.RPI / 180.0
        end
        @test length(unique(round.(grc.latdeg; digits=6))) == ng
        @test length(unique(round.(grc.londeg; digits=6))) == ng

        # Counts add up across gridcells: the per-gridcell landunit indices map
        # back to exactly nl landunits, and each landunit belongs to a gridcell.
        @test bounds.endl == nl
        @test bounds.endc == nc
        @test bounds.endp == np
        # Every landunit's owning gridcell is in range, and the union covers all g.
        owners = unique(lun.gridcell[1:nl])
        @test sort(owners) == collect(1:ng)
        # Sum of per-gridcell landunit counts equals nl.
        per_g = [count(==(g), lun.gridcell[1:nl]) for g in 1:ng]
        @test sum(per_g) == nl
        @test all(per_g .>= 1)
        # Same add-up for columns and patches.
        @test sum(count(==(g), col.gridcell[1:nc]) for g in 1:ng) == nc
        @test sum(count(==(g), pch.gridcell[1:np]) for g in 1:ng) == np
    end

    @testset "ng=1 single-cell reproduces the old tiled result" begin
        # Old path: read single gridcell (full_grid=false) — byte-identical with
        # the historical reader — then tile with ncopies=1 (a no-op).
        surf_old = CLM.SurfaceInputData()
        CLM.surfrd_get_data!(surf_old, 1, 1, fsurdat)   # legacy, reads cell g=1
        CLM._tile_surface_data!(surf_old, 1)

        # New full-grid read of the SAME file with ng = 1 would read all 6 cells,
        # so to compare the single-cell slice we read just g=1 via full_grid too.
        surf_new = CLM.SurfaceInputData()
        CLM.surfrd_get_data!(surf_new, 1, 1, fsurdat; full_grid=true)

        # The numeric surfdata fields must be byte-identical between the legacy
        # reader and the full_grid reader for the same single gridcell.
        for fld in (:wt_lunit, :wt_nat_patch, :pct_sand, :pct_clay,
                    :lakedepth, :slope)
            @test getfield(surf_old, fld) == getfield(surf_new, fld)
        end

        # initGridCells! with the legacy surf (no grid_latdeg) falls back to the
        # scalar lat/lon, exactly as before.
        (nl, nc, np) = CLM.count_subgrid_elements(surf_old, 1)
        bounds = CLM.BoundsType(begg=1, endg=1, begl=1, endl=nl,
                                begc=1, endc=nc, begp=1, endp=np,
                                level=CLM.BOUNDS_LEVEL_CLUMP, clump_index=1)
        grc = CLM.GridcellData();  CLM.gridcell_init!(grc, 1)
        lun = CLM.LandunitData();  CLM.landunit_init!(lun, nl)
        col = CLM.ColumnData();    CLM.column_init!(col, nc)
        pch = CLM.PatchData();     CLM.patch_init!(pch, np)
        CLM.initGridCells!(bounds, surf_old, grc, lun, col, pch;
                           lat=42.0, lon=-99.0)
        @test grc.latdeg[1] ≈ 42.0       # scalar fallback honored
        @test grc.londeg[1] ≈ -99.0
    end

    rm(tmpdir; recursive=true, force=true)
end
