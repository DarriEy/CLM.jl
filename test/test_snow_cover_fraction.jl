@testset "Snow Cover Fraction Module" begin
    # ------------------------------------------------------------------
    # Tests for SnowCoverFractionBaseMod, SnowCoverFractionSwensonLawrence2012Mod,
    # SnowCoverFractionNiuYang2007Mod, and SnowCoverFractionFactoryMod port.
    #
    # Verifies:
    #   1. Type construction and defaults
    #   2. Factory function (create_snow_cover_fraction)
    #   3. Initialization (snow_cover_fraction_init!)
    #   4. calc_frac_sno_eff!
    #   5. frac_snow_during_melt (SwensonLawrence2012)
    #   6. frac_snow_during_melt (NiuYang2007 returns NaN)
    #   7. update_snow_depth_and_frac! (SwensonLawrence2012)
    #   8. update_snow_depth_and_frac! (NiuYang2007)
    #   9. add_newsnow_to_intsnow! (SwensonLawrence2012)
    #  10. add_newsnow_to_intsnow! (NiuYang2007)
    #  11. Physics: accumulation then melt cycle (SwensonLawrence2012)
    #  12. Physics: NiuYang2007 frac_sno depends on snow depth and density
    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    # 1. Type construction and defaults
    # ------------------------------------------------------------------
    @testset "SwensonLawrence2012 defaults" begin
        scf = CLM.SnowCoverFractionSwensonLawrence2012()
        @test scf.int_snow_max == 2000.0
        @test scf.accum_factor == 0.1
        @test isempty(scf.n_melt)
        @test scf isa CLM.SnowCoverFractionBase
    end

    @testset "NiuYang2007 defaults" begin
        scf = CLM.SnowCoverFractionNiuYang2007()
        @test scf.zlnd == 0.01
        @test scf isa CLM.SnowCoverFractionBase
    end

    # ------------------------------------------------------------------
    # 2. Factory function
    # ------------------------------------------------------------------
    @testset "Factory: create_snow_cover_fraction" begin
        scf_sl = CLM.create_snow_cover_fraction("SwensonLawrence2012")
        @test scf_sl isa CLM.SnowCoverFractionSwensonLawrence2012

        scf_ny = CLM.create_snow_cover_fraction("NiuYang2007")
        @test scf_ny isa CLM.SnowCoverFractionNiuYang2007

        # Default is SwensonLawrence2012
        scf_def = CLM.create_snow_cover_fraction()
        @test scf_def isa CLM.SnowCoverFractionSwensonLawrence2012

        # Invalid method
        @test_throws ErrorException CLM.create_snow_cover_fraction("InvalidMethod")
    end

    # ------------------------------------------------------------------
    # 3. Initialization
    # ------------------------------------------------------------------
    @testset "SwensonLawrence2012 init" begin
        nc = 4
        scf = CLM.SnowCoverFractionSwensonLawrence2012()

        col_lun_itype = [CLM.ISTSOIL, CLM.ISTSOIL, CLM.ISTICE, CLM.ISTSOIL]
        col_gridcell = [1, 1, 1, 1]
        col_topo_std = [100.0, 50.0, 200.0, 5.0]  # note: 5.0 < 10.0 -> clamp to 10.0

        CLM.snow_cover_fraction_init!(scf, nc;
            col_lun_itype=col_lun_itype,
            col_gridcell=col_gridcell,
            col_topo_std=col_topo_std,
            n_melt_coef=200.0,
            n_melt_glcmec=10.0)

        @test length(scf.n_melt) == nc
        @test scf.n_melt[1] ≈ 200.0 / 100.0  # 2.0
        @test scf.n_melt[2] ≈ 200.0 / 50.0   # 4.0
        # Column 3 is ISTICE but no allow_multiple_columns_grc provided (empty)
        # so it falls through to the topo_std path
        @test scf.n_melt[3] ≈ 200.0 / 200.0  # 1.0
        @test scf.n_melt[4] ≈ 200.0 / 10.0   # 20.0 (clamped from 5.0 to 10.0)
    end

    @testset "SwensonLawrence2012 init with glc_mec" begin
        nc = 2
        scf = CLM.SnowCoverFractionSwensonLawrence2012()

        col_lun_itype = [CLM.ISTICE, CLM.ISTSOIL]
        col_gridcell = [1, 1]
        col_topo_std = [100.0, 50.0]
        allow_mc = [true]  # gridcell 1 allows multiple columns

        CLM.snow_cover_fraction_init!(scf, nc;
            col_lun_itype=col_lun_itype,
            col_gridcell=col_gridcell,
            col_topo_std=col_topo_std,
            allow_multiple_columns_grc=allow_mc,
            n_melt_coef=200.0,
            n_melt_glcmec=10.0)

        @test scf.n_melt[1] == 10.0   # ISTICE with allow_multiple_columns -> n_melt_glcmec
        @test scf.n_melt[2] ≈ 200.0 / 50.0  # regular soil
    end

    @testset "NiuYang2007 init" begin
        scf = CLM.SnowCoverFractionNiuYang2007()
        CLM.snow_cover_fraction_init!(scf; zlnd=0.02, use_subgrid_fluxes=false)
        @test scf.zlnd == 0.02

        # Should error with use_subgrid_fluxes=true
        scf2 = CLM.SnowCoverFractionNiuYang2007()
        @test_throws ErrorException CLM.snow_cover_fraction_init!(scf2;
            use_subgrid_fluxes=true)
    end

    # ------------------------------------------------------------------
    # 4. calc_frac_sno_eff!
    # ------------------------------------------------------------------
    @testset "calc_frac_sno_eff!" begin
        nc = 5
        frac_sno = [0.5, 0.8, 0.3, 0.0, 0.6]
        frac_sno_eff = zeros(nc)
        lun_itype = [CLM.ISTSOIL, CLM.ISTDLAK, CLM.ISTSOIL, CLM.ISTSOIL, CLM.ISTSOIL]
        urbpoi = [false, false, true, false, false]
        mask = trues(nc)
        bounds = 1:nc

        # With use_subgrid_fluxes=true
        CLM.calc_frac_sno_eff!(frac_sno_eff, frac_sno, lun_itype, urbpoi,
            mask, bounds; use_subgrid_fluxes=true)

        @test frac_sno_eff[1] ≈ 0.5   # soil, subgrid -> fractional
        @test frac_sno_eff[2] ≈ 1.0   # lake -> forced to 1 (frac_sno > 0)
        @test frac_sno_eff[3] ≈ 1.0   # urban -> forced to 1 (frac_sno > 0)
        @test frac_sno_eff[4] ≈ 0.0   # soil, frac_sno=0 -> 0
        @test frac_sno_eff[5] ≈ 0.6   # soil, subgrid -> fractional

        # With use_subgrid_fluxes=false -> all forced to 0/1
        frac_sno_eff2 = zeros(nc)
        CLM.calc_frac_sno_eff!(frac_sno_eff2, frac_sno, lun_itype, urbpoi,
            mask, bounds; use_subgrid_fluxes=false)

        @test frac_sno_eff2[1] ≈ 1.0   # forced to 1
        @test frac_sno_eff2[4] ≈ 0.0   # forced to 0
    end

    # ------------------------------------------------------------------
    # 5. frac_snow_during_melt (SwensonLawrence2012)
    # ------------------------------------------------------------------
    @testset "frac_snow_during_melt SwensonLawrence2012" begin
        scf = CLM.SnowCoverFractionSwensonLawrence2012()
        scf.n_melt = [2.0]

        # When h2osno_total == int_snow -> smr = 1.0 -> frac_sno should be ~ 1.0
        fs = CLM.frac_snow_during_melt(scf, 1, 100.0, 100.0)
        @test fs ≈ 1.0 atol=1e-10

        # When h2osno_total is half of int_snow -> smr = 0.5
        fs2 = CLM.frac_snow_during_melt(scf, 1, 50.0, 100.0)
        @test 0.0 < fs2 < 1.0

        # When h2osno_total -> 0 -> smr -> 0 -> frac_sno -> 0
        fs3 = CLM.frac_snow_during_melt(scf, 1, 0.001, 100.0)
        @test fs3 ≈ 0.0 atol=0.01

        # int_snow_max limiting
        scf.int_snow_max = 50.0
        # int_snow=100 -> limited to 50, h2osno_total=50 -> smr=1.0
        fs4 = CLM.frac_snow_during_melt(scf, 1, 50.0, 100.0)
        @test fs4 ≈ 1.0 atol=1e-10
    end

    # ------------------------------------------------------------------
    # 6. frac_snow_during_melt (NiuYang2007 returns NaN)
    # ------------------------------------------------------------------
    @testset "frac_snow_during_melt NiuYang2007" begin
        scf = CLM.SnowCoverFractionNiuYang2007()
        fs = CLM.frac_snow_during_melt(scf, 1, 100.0, 100.0)
        @test isnan(fs)
    end

    # ------------------------------------------------------------------
    # 7. update_snow_depth_and_frac! (SwensonLawrence2012)
    # ------------------------------------------------------------------
    @testset "update_snow_depth_and_frac! SwensonLawrence2012" begin
        nc = 3
        scf = CLM.SnowCoverFractionSwensonLawrence2012(
            int_snow_max=2000.0, accum_factor=0.1,
            n_melt=[2.0, 2.0, 2.0])

        # Case 1: column 1 has no existing snow, new snow arrives
        # Case 2: column 2 has existing snow, melt occurs
        # Case 3: column 3 has existing snow, new snow arrives
        frac_sno = [0.0, 0.5, 0.3]
        frac_sno_eff = [0.0, 0.5, 0.3]
        snow_depth = [0.0, 0.1, 0.05]
        lun_itype = [CLM.ISTSOIL, CLM.ISTSOIL, CLM.ISTSOIL]
        urbpoi = [false, false, false]
        h2osno_total = [0.0, 50.0, 30.0]
        snowmelt = [0.0, 5.0, 0.0]
        int_snow = [0.0, 100.0, 60.0]
        newsnow = [10.0, 0.0, 5.0]
        bifall = [100.0, 100.0, 100.0]
        mask = trues(nc)
        bounds = 1:nc

        CLM.update_snow_depth_and_frac!(scf,
            frac_sno, frac_sno_eff, snow_depth,
            lun_itype, urbpoi, h2osno_total, snowmelt, int_snow,
            newsnow, bifall, mask, bounds;
            use_subgrid_fluxes=true)

        # Column 1: h2osno=0, newsnow=10 -> frac_sno = tanh(0.1 * 10) = tanh(1.0)
        @test frac_sno[1] ≈ tanh(1.0)
        @test snow_depth[1] > 0.0  # new snow added

        # Column 2: h2osno=50, melt=5, no newsnow -> frac_sno via melt curve
        @test 0.0 < frac_sno[2] <= 1.0
        # Snow depth should remain positive (no new snow added, existing depth unchanged)
        @test snow_depth[2] ≈ 0.1  # no new snow, depth unchanged

        # Column 3: h2osno=30, no melt, newsnow=5 -> frac_sno updated
        @test frac_sno[3] > 0.3  # should have increased from accumulation
        @test snow_depth[3] > 0.05  # increased from new snow
    end

    # ------------------------------------------------------------------
    # 8. update_snow_depth_and_frac! (NiuYang2007)
    # ------------------------------------------------------------------
    @testset "update_snow_depth_and_frac! NiuYang2007" begin
        nc = 3
        scf = CLM.SnowCoverFractionNiuYang2007(zlnd=0.01)

        # Case 1: new snow on bare ground
        # Case 2: existing snow with more new snow
        # Case 3: zero snow
        frac_sno = [0.0, 0.5, 0.0]
        frac_sno_eff = [0.0, 0.5, 0.0]
        snow_depth = [0.0, 0.1, 0.0]
        lun_itype = [CLM.ISTSOIL, CLM.ISTSOIL, CLM.ISTSOIL]
        urbpoi = [false, false, false]
        h2osno_total = [0.0, 50.0, 0.0]
        snowmelt = [0.0, 0.0, 0.0]
        int_snow = [0.0, 100.0, 0.0]
        newsnow = [20.0, 5.0, 0.0]
        bifall = [100.0, 100.0, 100.0]
        mask = trues(nc)
        bounds = 1:nc

        CLM.update_snow_depth_and_frac!(scf,
            frac_sno, frac_sno_eff, snow_depth,
            lun_itype, urbpoi, h2osno_total, snowmelt, int_snow,
            newsnow, bifall, mask, bounds;
            use_subgrid_fluxes=false)

        # Column 1: h2osno=0 -> depth reset to 0, then +20/100 = 0.2
        @test snow_depth[1] ≈ 0.2
        @test frac_sno[1] > 0.0

        # Column 2: existing snow + new snow
        @test snow_depth[2] ≈ 0.1 + 5.0/100.0  # 0.15
        @test frac_sno[2] > 0.0

        # Column 3: no snow, no new snow -> snow_depth stays 0
        @test snow_depth[3] ≈ 0.0
        @test frac_sno[3] ≈ 0.0
    end

    # ------------------------------------------------------------------
    # 9. add_newsnow_to_intsnow! (SwensonLawrence2012)
    # ------------------------------------------------------------------
    @testset "add_newsnow_to_intsnow! SwensonLawrence2012" begin
        nc = 2
        scf = CLM.SnowCoverFractionSwensonLawrence2012(
            int_snow_max=2000.0, accum_factor=0.1,
            n_melt=[2.0, 2.0])

        int_snow = [100.0, 50.0]
        newsnow = [10.0, 0.0]
        h2osno_total = [80.0, 50.0]
        frac_sno = [0.7, 0.5]
        mask = trues(nc)
        bounds = 1:nc

        CLM.add_newsnow_to_intsnow!(scf, int_snow, newsnow, h2osno_total,
            frac_sno, mask, bounds)

        # Column 1: newsnow > 0 -> int_snow is reset and then newsnow added
        @test int_snow[1] > 0.0

        # Column 2: newsnow = 0 -> int_snow just gets 0 added
        @test int_snow[2] ≈ 50.0
    end

    # ------------------------------------------------------------------
    # 10. add_newsnow_to_intsnow! (NiuYang2007)
    # ------------------------------------------------------------------
    @testset "add_newsnow_to_intsnow! NiuYang2007" begin
        nc = 2
        scf = CLM.SnowCoverFractionNiuYang2007()

        int_snow = [100.0, 50.0]
        newsnow = [10.0, 5.0]
        h2osno_total = [80.0, 50.0]
        frac_sno = [0.7, 0.5]
        mask = trues(nc)
        bounds = 1:nc

        CLM.add_newsnow_to_intsnow!(scf, int_snow, newsnow, h2osno_total,
            frac_sno, mask, bounds)

        @test int_snow[1] ≈ 110.0
        @test int_snow[2] ≈ 55.0
    end

    # ------------------------------------------------------------------
    # 11. Physics: accumulation then melt cycle (SwensonLawrence2012)
    # ------------------------------------------------------------------
    @testset "SwensonLawrence2012 accumulation-melt cycle" begin
        nc = 1
        scf = CLM.SnowCoverFractionSwensonLawrence2012(
            int_snow_max=2000.0, accum_factor=0.1,
            n_melt=[2.0])

        frac_sno = [0.0]
        frac_sno_eff = [0.0]
        snow_depth = [0.0]
        lun_itype = [CLM.ISTSOIL]
        urbpoi = [false]
        h2osno_total = [0.0]
        snowmelt = [0.0]
        int_snow = [0.0]
        bifall = [100.0]
        mask = trues(nc)
        bounds = 1:nc

        # Phase 1: Accumulate snow
        for _ in 1:10
            newsnow = [5.0]
            CLM.update_snow_depth_and_frac!(scf,
                frac_sno, frac_sno_eff, snow_depth,
                lun_itype, urbpoi, h2osno_total, snowmelt, int_snow,
                newsnow, bifall, mask, bounds; use_subgrid_fluxes=true)
            h2osno_total[1] += 5.0
            CLM.add_newsnow_to_intsnow!(scf, int_snow, newsnow, h2osno_total,
                frac_sno, mask, bounds)
        end

        # After accumulation, frac_sno should be close to 1.0
        @test frac_sno[1] > 0.9
        frac_sno_peak = frac_sno[1]

        # Phase 2: Melt snow
        for i in 1:8
            newsnow = [0.0]
            snowmelt[1] = 5.0
            h2osno_total[1] = max(0.0, h2osno_total[1] - 5.0)
            CLM.update_snow_depth_and_frac!(scf,
                frac_sno, frac_sno_eff, snow_depth,
                lun_itype, urbpoi, h2osno_total, snowmelt, int_snow,
                newsnow, bifall, mask, bounds; use_subgrid_fluxes=true)
        end

        # After melt, frac_sno should have decreased from peak
        @test frac_sno[1] < frac_sno_peak

        # When h2osno_total reaches 0, frac_sno should eventually be 0
        h2osno_total[1] = 0.0
        snowmelt[1] = 0.0
        CLM.update_snow_depth_and_frac!(scf,
            frac_sno, frac_sno_eff, snow_depth,
            lun_itype, urbpoi, [0.0], [0.0], int_snow,
            [0.0], bifall, mask, bounds; use_subgrid_fluxes=true)
        @test frac_sno[1] ≈ 0.0
    end

    # ------------------------------------------------------------------
    # 12. Physics: NiuYang2007 frac_sno depends on snow depth and density
    # ------------------------------------------------------------------
    @testset "NiuYang2007 depth-density relationship" begin
        scf = CLM.SnowCoverFractionNiuYang2007(zlnd=0.01)

        nc = 1
        frac_sno = [0.0]
        frac_sno_eff = [0.0]
        snow_depth = [0.0]
        lun_itype = [CLM.ISTSOIL]
        urbpoi = [false]
        mask = trues(nc)
        bounds = 1:nc

        # With very deep snow, frac_sno should approach 1.0
        CLM.update_snow_depth_and_frac!(scf,
            frac_sno, frac_sno_eff, snow_depth,
            lun_itype, urbpoi,
            [0.0],    # h2osno_total (start at 0 to trigger reset)
            [0.0],    # snowmelt
            [0.0],    # int_snow
            [100.0],  # newsnow: 100 mm
            [100.0],  # bifall: 100 kg/m3
            mask, bounds; use_subgrid_fluxes=false)

        @test snow_depth[1] ≈ 1.0  # 100/100 = 1.0 m
        @test frac_sno[1] > 0.9     # deep snow -> near full coverage

        # With very shallow snow, frac_sno should be small
        frac_sno2 = [0.0]
        frac_sno_eff2 = [0.0]
        snow_depth2 = [0.0]

        CLM.update_snow_depth_and_frac!(scf,
            frac_sno2, frac_sno_eff2, snow_depth2,
            lun_itype, urbpoi,
            [0.0],   # h2osno_total
            [0.0],   # snowmelt
            [0.0],   # int_snow
            [0.1],   # newsnow: 0.1 mm
            [100.0], # bifall
            mask, bounds; use_subgrid_fluxes=false)

        @test snow_depth2[1] ≈ 0.001  # 0.1/100
        @test frac_sno2[1] < frac_sno[1]  # shallower snow -> less coverage
    end

    # ------------------------------------------------------------------
    # 13. Mask filtering
    # ------------------------------------------------------------------
    @testset "Mask filtering" begin
        nc = 3
        scf = CLM.SnowCoverFractionNiuYang2007(zlnd=0.01)

        frac_sno = [0.0, 0.0, 0.0]
        frac_sno_eff = [0.0, 0.0, 0.0]
        snow_depth = [0.0, 0.0, 0.0]
        lun_itype = [CLM.ISTSOIL, CLM.ISTSOIL, CLM.ISTSOIL]
        urbpoi = [false, false, false]
        newsnow = [50.0, 50.0, 50.0]
        bifall = [100.0, 100.0, 100.0]

        # Only process column 2
        mask = BitVector([false, true, false])
        bounds = 1:nc

        CLM.update_snow_depth_and_frac!(scf,
            frac_sno, frac_sno_eff, snow_depth,
            lun_itype, urbpoi,
            [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0],
            newsnow, bifall, mask, bounds; use_subgrid_fluxes=false)

        # Only column 2 should be modified
        @test snow_depth[1] ≈ 0.0
        @test snow_depth[2] ≈ 0.5  # 50/100
        @test snow_depth[3] ≈ 0.0
        @test frac_sno[1] ≈ 0.0
        @test frac_sno[2] > 0.0
        @test frac_sno[3] ≈ 0.0
    end
end
