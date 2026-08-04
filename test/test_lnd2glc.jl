#!/usr/bin/env julia
# ==========================================================================
# lnd2glcMod.F90 port: the land -> ice-sheet (CISM) coupling fields.
#
# CLM.jl carries a genuine glacier elevation-class subgrid (set_landunit_ice! in
# init_gridcells.jl builds one ISTICE column per non-zero wt_glc_mec[g,m], with
# col.itype = the ice class and its own topo_col), so update_lnd2glc! reads the
# per-class surface mass balance straight out of CLM's own state — no lapse-rate
# stand-in is needed by a coupled driver.
#
# What is pinned here:
#   * the (maxpatch_glc+1, ng) shape and the class-0-is-bare-land index map
#   * the defaults sent OUTSIDE the do_smb filter (qice = 0 for conservation,
#     tsrf = TFRZ, topo = 0)
#   * per-elevation-class routing: an ISTICE column lands in the row for its own
#     ice class, an ISTSOIL column in the bare-land row
#   * bareland_normalization CONSERVES mass: the bare-land flux the ice sheet
#     receives, spread over the whole non-glacier fraction, equals the flux CLM
#     computed on the natural-veg landunit times that landunit's area
#   * `init = true` leaves qice at 0 (qflx_glcice is not valid during init)
#   * the double-assignment guard
# ==========================================================================

using Test
using CLM

@testset "lnd2glc (land -> ice sheet)" begin

    # varpar is only populated by varpar_init!; initialize it here so this file
    # does not depend on whether an earlier test happened to call it.
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
    nlevsno  = CLM.varpar.nlevsno
    nlevgrnd = CLM.varpar.nlevgrnd
    nlev     = nlevsno + nlevgrnd

    # ---- domain: 1 gridcell = 60% glacier (2 elevation classes) + 30% natural
    #      veg + 10% lake. This is exactly the worked example in the Fortran
    #      `bareland_normalization` docstring.
    function build_glc_domain()
        col = CLM.ColumnData()
        #        c1: glc class 1   c2: glc class 2   c3: natveg   c4: lake
        col.gridcell = [1, 1, 1, 1]
        col.landunit = [1, 1, 2, 3]
        col.itype    = [CLM.ice_class_to_col_itype(1), CLM.ice_class_to_col_itype(2), 1, 1]
        col.wtgcell  = [0.4, 0.2, 0.3, 0.1]
        col.wtlunit  = [0.4 / 0.6, 0.2 / 0.6, 1.0, 1.0]
        col.active   = [true, true, true, true]
        col.is_hillslope_column = [false, false, false, false]

        lun = CLM.LandunitData()
        lun.itype      = [CLM.ISTICE, CLM.ISTSOIL, CLM.ISTDLAK]
        lun.gridcell   = [1, 1, 1]
        lun.wtgcell    = [0.6, 0.3, 0.1]
        lun.active     = [true, true, true]
        lun.urbpoi     = [false, false, false]
        lun.canyon_hwr = [0.0, 0.0, 0.0]

        grc = CLM.GridcellData()
        grc.lat = [0.7]; grc.lon = [0.1]
        grc.latdeg = [40.0]; grc.londeg = [-105.0]
        grc.area = [100.0]
        # landunit_indices[ltype, g] — what get_landunit_weight looks up.
        li = fill(CLM.ISPVAL, CLM.MAX_LUNIT, 1)
        li[CLM.ISTICE,  1] = 1
        li[CLM.ISTSOIL, 1] = 2
        li[CLM.ISTDLAK, 1] = 3
        grc.landunit_indices = li

        temp = CLM.TemperatureData()
        temp.t_soisno_col = zeros(4, nlev)
        # Fortran t_soisno(c,1) = the TOP GROUND layer = Julia column nlevsno+1.
        temp.t_soisno_col[:, nlevsno + 1] = [265.0, 260.0, 275.0, 280.0]

        topo = CLM.TopoData()
        topo.topo_col = [1500.0, 2500.0, 800.0, 700.0]

        wfb = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(wfb, 4, 4, 3, 1)
        wfb.wf.qflx_glcice_col = zeros(4)

        return (; col, lun, grc, temp, topo, wfb)
    end

    @testset "allocation shape and the class-0 index map" begin
        l2g = CLM.Lnd2GlcData()
        @test size(l2g.qice_grc) == (0, 0)

        CLM.lnd2glc_init!(l2g, 3; maxpatch_glc = 10)
        @test size(l2g.qice_grc) == (11, 3)
        @test size(l2g.tsrf_grc) == (11, 3)
        @test size(l2g.topo_grc) == (11, 3)
        @test all(iszero, l2g.qice_grc)

        # Fortran class n -> Julia row n+1; class 0 (bare land) is row 1.
        @test CLM.lnd2glc_class_index(0) == 1
        @test CLM.lnd2glc_class_index(1) == 2
        @test CLM.lnd2glc_class_index(10) == 11

        CLM.lnd2glc_clean!(l2g)
        @test size(l2g.qice_grc) == (0, 0)
    end

    @testset "defaults outside the do_smb filter" begin
        d = build_glc_domain()
        l2g = CLM.Lnd2GlcData()
        CLM.lnd2glc_init!(l2g, 1; maxpatch_glc = 3)

        # Empty filter => nothing is in do_smb => everything keeps the defaults.
        CLM.update_lnd2glc!(l2g, d.temp, d.wfb, d.topo, d.col, d.lun, d.grc,
                            falses(4), 1:4, 1:1; nlevsno = nlevsno)

        # qice MUST be 0 outside the filter (conservation: the ice sheet must not
        # be handed mass CLM never computed).
        @test all(iszero, l2g.qice_grc)
        @test all(≈(CLM.TFRZ), l2g.tsrf_grc)
        @test all(iszero, l2g.topo_grc)
    end

    @testset "per-elevation-class routing of tsrf / topo / qice" begin
        d = build_glc_domain()
        # SMB: class 1 gains ice, class 2 gains more, bare land loses a little.
        d.wfb.wf.qflx_glcice_col .= [1.0e-5, 3.0e-5, -2.0e-6, 0.0]

        l2g = CLM.Lnd2GlcData()
        CLM.lnd2glc_init!(l2g, 1; maxpatch_glc = 3)

        # do_smb: the two glacier columns + the natural-veg column (the lake
        # column is never in the SMB filter, matching filters.jl).
        do_smb = [true, true, true, false]
        CLM.update_lnd2glc!(l2g, d.temp, d.wfb, d.topo, d.col, d.lun, d.grc,
                            do_smb, 1:4, 1:1; nlevsno = nlevsno)

        i0 = CLM.lnd2glc_class_index(0)   # bare land
        i1 = CLM.lnd2glc_class_index(1)
        i2 = CLM.lnd2glc_class_index(2)

        # tsrf = t_soisno(c,1), topo = topo_col(c), each in its OWN class row.
        @test l2g.tsrf_grc[i1, 1] ≈ 265.0
        @test l2g.tsrf_grc[i2, 1] ≈ 260.0
        @test l2g.tsrf_grc[i0, 1] ≈ 275.0
        @test l2g.topo_grc[i1, 1] ≈ 1500.0
        @test l2g.topo_grc[i2, 1] ≈ 2500.0
        @test l2g.topo_grc[i0, 1] ≈ 800.0

        # Glacier columns: flux normalization is exactly 1.
        @test l2g.qice_grc[i1, 1] ≈ 1.0e-5
        @test l2g.qice_grc[i2, 1] ≈ 3.0e-5

        # Bare land: normalized by wtgcell / (1 - area_glacier) = 0.3/0.4.
        @test l2g.qice_grc[i0, 1] ≈ -2.0e-6 * (0.3 / 0.4)

        # The lake column contributed nothing anywhere; class 3 keeps its default.
        i3 = CLM.lnd2glc_class_index(3)
        @test l2g.qice_grc[i3, 1] == 0.0
        @test l2g.tsrf_grc[i3, 1] ≈ CLM.TFRZ
    end

    @testset "bareland_normalization conserves the bare-land flux" begin
        d = build_glc_domain()
        # natveg column c3 has wtgcell 0.3; the non-glacier fraction is 0.4.
        @test CLM.bareland_normalization(3, d.col, d.lun, d.grc) ≈ 0.3 / 0.4

        # Conservation statement: (flux CLM computed) * (its landunit area) must
        # equal (flux the ice sheet applies over bare land) * (bare-land area).
        flux_clm = 1.25e-5
        f_sent   = flux_clm * CLM.bareland_normalization(3, d.col, d.lun, d.grc)
        area_bare = 1.0 - CLM.get_landunit_weight(d.grc, d.lun, 1, CLM.ISTICE)
        @test f_sent * area_bare ≈ flux_clm * d.col.wtgcell[3]

        # All-glacier gridcell: factor is arbitrary, Fortran returns 1.
        d2 = build_glc_domain()
        d2.lun.wtgcell[1] = 1.0; d2.lun.wtgcell[2] = 0.0; d2.lun.wtgcell[3] = 0.0
        @test CLM.bareland_normalization(3, d2.col, d2.lun, d2.grc) == 1.0
    end

    @testset "init = true leaves qice at its default" begin
        d = build_glc_domain()
        d.wfb.wf.qflx_glcice_col .= [1.0e-5, 3.0e-5, -2.0e-6, 0.0]
        l2g = CLM.Lnd2GlcData()
        CLM.lnd2glc_init!(l2g, 1; maxpatch_glc = 3)

        CLM.update_lnd2glc!(l2g, d.temp, d.wfb, d.topo, d.col, d.lun, d.grc,
                            [true, true, true, false], 1:4, 1:1; init = true, nlevsno = nlevsno)

        # t_soisno / topo_col ARE valid during init, qflx_glcice is NOT.
        @test all(iszero, l2g.qice_grc)
        @test l2g.tsrf_grc[CLM.lnd2glc_class_index(1), 1] ≈ 265.0
        @test l2g.topo_grc[CLM.lnd2glc_class_index(2), 1] ≈ 2500.0
    end

    @testset "double assignment of one (gridcell, class) slot is an error" begin
        d = build_glc_domain()
        # Give the natural-veg landunit a SECOND column: both map to class 0.
        d.col.gridcell = [1, 1, 1, 1, 1]
        d.col.landunit = [1, 1, 2, 3, 2]
        d.col.itype    = [CLM.ice_class_to_col_itype(1), CLM.ice_class_to_col_itype(2), 1, 1, 1]
        d.col.wtgcell  = [0.4, 0.2, 0.15, 0.1, 0.15]
        d.col.active   = [true, true, true, true, true]
        d.col.is_hillslope_column = falses(5)
        temp = CLM.TemperatureData()
        temp.t_soisno_col = zeros(5, nlev)
        temp.t_soisno_col[:, nlevsno + 1] .= 275.0
        topo = CLM.TopoData(); topo.topo_col = [1500.0, 2500.0, 800.0, 700.0, 800.0]
        wfb = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(wfb, 5, 5, 3, 1)
        wfb.wf.qflx_glcice_col = zeros(5)

        l2g = CLM.Lnd2GlcData()
        CLM.lnd2glc_init!(l2g, 1; maxpatch_glc = 3)
        do_smb = [true, true, true, false, true]

        @test_throws ErrorException CLM.update_lnd2glc!(l2g, temp, wfb, topo,
            d.col, d.lun, d.grc, do_smb, 1:5, 1:1; nlevsno = nlevsno)

        # ...unless use_hillslope bypasses it (ESCOMP/CTSM#204).
        @test nothing === CLM.update_lnd2glc!(l2g, temp, wfb, topo,
            d.col, d.lun, d.grc, do_smb, 1:5, 1:1; use_hillslope = true, nlevsno = nlevsno)
    end
end
