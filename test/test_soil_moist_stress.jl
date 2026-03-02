@testset "SoilMoistStress" begin
    # Save and restore varpar state
    vp = CLM.varpar
    saved_nlevsno = vp.nlevsno
    saved_nlevgrnd = vp.nlevgrnd
    saved_nlevurb = vp.nlevurb
    saved_nlevmaxurbgrnd = vp.nlevmaxurbgrnd
    saved_nlevsoi = vp.nlevsoi

    vp.nlevsno = 5
    vp.nlevgrnd = 10
    vp.nlevurb = 5
    vp.nlevmaxurbgrnd = 10
    vp.nlevsoi = 8

    nlevsno = vp.nlevsno
    nlevgrnd = vp.nlevgrnd
    nlevmaxurbgrnd = vp.nlevmaxurbgrnd
    joff = nlevsno

    nc = 4   # number of columns
    np = 4   # number of patches

    # Save module control state
    saved_method = CLM.soil_moist_stress_ctrl.root_moist_stress_method
    saved_perchroot = CLM.soil_moist_stress_ctrl.perchroot
    saved_perchroot_alt = CLM.soil_moist_stress_ctrl.perchroot_alt

    @testset "init_root_moist_stress!" begin
        CLM.soil_moist_stress_ctrl.root_moist_stress_method = 99
        CLM.init_root_moist_stress!()
        @test CLM.soil_moist_stress_ctrl.root_moist_stress_method == CLM.MOIST_STRESS_CLM_DEFAULT
    end

    @testset "set_perchroot_opt!" begin
        CLM.set_perchroot_opt!(true, false)
        @test CLM.soil_moist_stress_ctrl.perchroot == true
        @test CLM.soil_moist_stress_ctrl.perchroot_alt == false

        CLM.set_perchroot_opt!(false, true)
        @test CLM.soil_moist_stress_ctrl.perchroot == false
        @test CLM.soil_moist_stress_ctrl.perchroot_alt == true

        CLM.set_perchroot_opt!(false, false)
    end

    @testset "calc_effective_soilporosity!" begin
        watsat = fill(0.45, nc, nlevgrnd)
        h2osoi_ice = zeros(nc, nlevsno + nlevgrnd)
        col_dz = fill(0.1, nc, nlevsno + nlevgrnd)
        eff_por = zeros(nc, nlevgrnd)
        mask = trues(nc)
        bounds = 1:nc

        # No ice → eff_por == watsat
        CLM.calc_effective_soilporosity!(watsat, h2osoi_ice, col_dz, eff_por,
                                         mask, bounds, nlevgrnd, nlevsno)
        for c in 1:nc, j in 1:nlevgrnd
            @test eff_por[c, j] ≈ 0.45
        end

        # With some ice
        h2osoi_ice[1, 1 + joff] = 0.1 * CLM.DENICE * 0.2  # vol_ice = 0.2
        CLM.calc_effective_soilporosity!(watsat, h2osoi_ice, col_dz, eff_por,
                                         mask, bounds, nlevgrnd, nlevsno)
        @test eff_por[1, 1] ≈ 0.45 - 0.2 atol=1e-10

        # With ice exceeding porosity → vol_ice clamped to watsat
        h2osoi_ice[2, 1 + joff] = 0.1 * CLM.DENICE * 0.5  # vol_ice = 0.5 > watsat=0.45
        CLM.calc_effective_soilporosity!(watsat, h2osoi_ice, col_dz, eff_por,
                                         mask, bounds, nlevgrnd, nlevsno)
        @test eff_por[2, 1] ≈ 0.0 atol=1e-10

        # Masked columns are skipped
        eff_por_saved = copy(eff_por)
        mask[3] = false
        h2osoi_ice[3, 1 + joff] = 999.0
        CLM.calc_effective_soilporosity!(watsat, h2osoi_ice, col_dz, eff_por,
                                         mask, bounds, nlevgrnd, nlevsno)
        @test eff_por[3, 1] == eff_por_saved[3, 1]
    end

    @testset "calc_effective_snowporosity!" begin
        h2osoi_ice = zeros(nc, nlevsno + nlevgrnd)
        col_dz = fill(0.1, nc, nlevsno + nlevgrnd)
        jtop = fill(-nlevsno + 1, nc)  # all snow layers active
        eff_por = zeros(nc, nlevsno + nlevgrnd)
        mask = trues(nc)
        bounds = 1:nc
        lbj = -nlevsno + 1

        # No ice → eff_por = 1.0
        CLM.calc_effective_snowporosity!(h2osoi_ice, col_dz, jtop, eff_por,
                                          mask, bounds, lbj, nlevsno)
        for j in lbj:0
            jj = j + nlevsno
            @test eff_por[1, jj] ≈ 1.0
        end

        # With some ice in top snow layer
        j_test = -nlevsno + 1
        jj_test = j_test + nlevsno
        h2osoi_ice[1, jj_test] = 0.1 * CLM.DENICE * 0.3  # vol_ice = 0.3
        CLM.calc_effective_snowporosity!(h2osoi_ice, col_dz, jtop, eff_por,
                                          mask, bounds, lbj, nlevsno)
        @test eff_por[1, jj_test] ≈ 0.7 atol=1e-10
    end

    @testset "calc_volumetric_h2oliq!" begin
        eff_porosity = fill(0.4, nc, nlevsno + nlevgrnd)
        h2osoi_liq = zeros(nc, nlevsno + nlevgrnd)
        col_dz = fill(0.1, nc, nlevsno + nlevgrnd)
        jtop = fill(1, nc)
        vol_liq = zeros(nc, nlevsno + nlevgrnd)
        mask = trues(nc)
        bounds = 1:nc

        # Set some liquid water in soil layer 1
        h2osoi_liq[1, 1 + joff] = 0.1 * CLM.DENH2O * 0.3  # vol = 0.3
        CLM.calc_volumetric_h2oliq!(eff_porosity, h2osoi_liq, col_dz, jtop,
                                     vol_liq, mask, bounds, 1, nlevgrnd, nlevsno)
        @test vol_liq[1, 1 + joff] ≈ 0.3 atol=1e-10

        # vol_liq capped by eff_porosity
        h2osoi_liq[2, 1 + joff] = 0.1 * CLM.DENH2O * 0.6  # vol = 0.6 > eff_por = 0.4
        CLM.calc_volumetric_h2oliq!(eff_porosity, h2osoi_liq, col_dz, jtop,
                                     vol_liq, mask, bounds, 1, nlevgrnd, nlevsno)
        @test vol_liq[2, 1 + joff] ≈ 0.4 atol=1e-10
    end

    @testset "array_normalization!" begin
        arr = zeros(np, nlevgrnd)
        mask = trues(np)
        bounds = 1:np

        arr[1, 1] = 0.3
        arr[1, 2] = 0.6
        arr[1, 3] = 0.1
        CLM.array_normalization!(arr, mask, bounds, nlevgrnd)
        @test sum(arr[1, :]) ≈ 1.0
        @test arr[1, 1] ≈ 0.3 atol=1e-10
        @test arr[1, 2] ≈ 0.6 atol=1e-10

        # All zeros → stays zero
        arr[2, :] .= 0.0
        CLM.array_normalization!(arr, mask, bounds, nlevgrnd)
        @test sum(arr[2, :]) ≈ 0.0

        # Masked patches skipped
        arr[3, 1] = 5.0
        mask[3] = false
        CLM.array_normalization!(arr, mask, bounds, nlevgrnd)
        @test arr[3, 1] ≈ 5.0  # unchanged
    end

    @testset "soil_suction_clapp_hornberger" begin
        sucsat = 100.0  # mm
        bsw = 5.0
        s_node = 0.5
        smp = CLM.soil_suction_clapp_hornberger(sucsat, s_node, bsw)
        @test smp ≈ -100.0 * 0.5^(-5.0)
        @test smp < 0.0  # matric potential is negative

        # Saturated: s_node = 1.0 → smp = -sucsat
        smp_sat = CLM.soil_suction_clapp_hornberger(sucsat, 1.0, bsw)
        @test smp_sat ≈ -100.0
    end

    @testset "calc_root_moist_stress_clm45default! — basic stress" begin
        CLM.set_perchroot_opt!(false, false)
        CLM.init_root_moist_stress!()

        nlevtot = nlevsno + nlevgrnd

        rootfr_unf = zeros(np, nlevgrnd)
        rootfr = zeros(np, nlevgrnd)
        rootr = zeros(np, nlevgrnd)
        btran = zeros(np)
        rresis = zeros(np, nlevgrnd)

        # PFT parameters
        npft = 2
        smpso = fill(-35000.0, npft)  # soil water potential at full stomatal opening (mm)
        smpsc = fill(-275000.0, npft) # soil water potential at full stomatal closure (mm)

        t_soisno = fill(290.0, nc, nlevtot)  # warm soil
        watsat = fill(0.45, nc, nlevgrnd)
        sucsat = fill(100.0, nc, nlevgrnd)
        bsw = fill(5.0, nc, nlevgrnd)
        eff_porosity = fill(0.40, nc, nlevgrnd)
        h2osoi_liqvol = zeros(nc, nlevtot)

        patch_column = [1, 2, 3, 4]
        patch_itype = [1, 1, 2, 1]
        mask_patch = trues(np)
        bounds_patch = 1:np

        # Set moderate soil moisture in all layers
        for c in 1:nc, j in 1:nlevgrnd
            h2osoi_liqvol[c, j + joff] = 0.25
        end

        # Set uniform root fraction
        for p in 1:np, j in 1:nlevgrnd
            rootfr[p, j] = 1.0 / nlevgrnd
        end

        CLM.calc_root_moist_stress_clm45default!(rootfr_unf, rootfr, rootr, btran, rresis,
            smpso, smpsc, t_soisno, watsat, sucsat, bsw, eff_porosity,
            h2osoi_liqvol, patch_column, patch_itype,
            mask_patch, bounds_patch, nlevgrnd, nlevsno)

        # btran should be positive (soil is moist and warm)
        for p in 1:np
            @test btran[p] > 0.0
            @test btran[p] <= 1.0
        end

        # rootr should sum to ~1.0 (normalized)
        for p in 1:np
            @test sum(rootr[p, :]) ≈ 1.0 atol=1e-10
        end

        # rresis should be between 0 and 1
        for p in 1:np, j in 1:nlevgrnd
            @test rresis[p, j] >= 0.0
            @test rresis[p, j] <= 1.0
        end
    end

    @testset "calc_root_moist_stress_clm45default! — dry soil gives zero stress" begin
        CLM.set_perchroot_opt!(false, false)

        nlevtot = nlevsno + nlevgrnd

        rootfr_unf = zeros(np, nlevgrnd)
        rootfr = zeros(np, nlevgrnd)
        rootr = zeros(np, nlevgrnd)
        btran = zeros(np)
        rresis = zeros(np, nlevgrnd)

        npft = 2
        smpso = fill(-35000.0, npft)
        smpsc = fill(-275000.0, npft)

        t_soisno = fill(290.0, nc, nlevtot)
        watsat = fill(0.45, nc, nlevgrnd)
        sucsat = fill(100.0, nc, nlevgrnd)
        bsw = fill(5.0, nc, nlevgrnd)
        eff_porosity = fill(0.40, nc, nlevgrnd)
        h2osoi_liqvol = zeros(nc, nlevtot)  # completely dry

        patch_column = [1, 2, 3, 4]
        patch_itype = [1, 1, 1, 1]
        mask_patch = trues(np)
        bounds_patch = 1:np

        rootfr .= 1.0 / nlevgrnd

        CLM.calc_root_moist_stress_clm45default!(rootfr_unf, rootfr, rootr, btran, rresis,
            smpso, smpsc, t_soisno, watsat, sucsat, bsw, eff_porosity,
            h2osoi_liqvol, patch_column, patch_itype,
            mask_patch, bounds_patch, nlevgrnd, nlevsno)

        # Dry soil → btran = 0, rootr = 0
        for p in 1:np
            @test btran[p] ≈ 0.0
            @test sum(rootr[p, :]) ≈ 0.0
        end
    end

    @testset "calc_root_moist_stress_clm45default! — frozen soil gives zero stress" begin
        CLM.set_perchroot_opt!(false, false)

        nlevtot = nlevsno + nlevgrnd

        rootfr_unf = zeros(np, nlevgrnd)
        rootfr = zeros(np, nlevgrnd)
        rootr = zeros(np, nlevgrnd)
        btran = zeros(np)
        rresis = zeros(np, nlevgrnd)

        npft = 2
        smpso = fill(-35000.0, npft)
        smpsc = fill(-275000.0, npft)

        t_soisno = fill(CLM.TFRZ - 5.0, nc, nlevtot)  # frozen
        watsat = fill(0.45, nc, nlevgrnd)
        sucsat = fill(100.0, nc, nlevgrnd)
        bsw = fill(5.0, nc, nlevgrnd)
        eff_porosity = fill(0.40, nc, nlevgrnd)
        h2osoi_liqvol = fill(0.25, nc, nlevtot)

        patch_column = [1, 2, 3, 4]
        patch_itype = [1, 1, 1, 1]
        mask_patch = trues(np)
        bounds_patch = 1:np

        rootfr .= 1.0 / nlevgrnd

        CLM.calc_root_moist_stress_clm45default!(rootfr_unf, rootfr, rootr, btran, rresis,
            smpso, smpsc, t_soisno, watsat, sucsat, bsw, eff_porosity,
            h2osoi_liqvol, patch_column, patch_itype,
            mask_patch, bounds_patch, nlevgrnd, nlevsno)

        # Frozen → btran = 0
        for p in 1:np
            @test btran[p] ≈ 0.0
        end
    end

    @testset "calc_root_moist_stress! — top-level dispatcher" begin
        CLM.set_perchroot_opt!(false, false)
        CLM.init_root_moist_stress!()

        nlevtot = nlevsno + nlevgrnd

        col = CLM.ColumnData()
        CLM.column_init!(col, nc)
        col.dz .= 0.1

        patchdata = CLM.PatchData()
        CLM.patch_init!(patchdata, np)
        for p in 1:np
            patchdata.column[p] = p
            patchdata.itype[p] = 1
        end

        temp = CLM.TemperatureData()
        CLM.temperature_init!(temp, np, nc, nc, 1)
        temp.t_soisno_col .= 290.0

        soilstate = CLM.SoilStateData()
        CLM.soilstate_init!(soilstate, np, nc)
        soilstate.watsat_col .= 0.45
        soilstate.sucsat_col .= 100.0
        soilstate.bsw_col .= 5.0
        soilstate.eff_porosity_col .= 0.40
        for p in 1:np, j in 1:nlevgrnd
            soilstate.rootfr_patch[p, j] = 1.0 / nlevgrnd
        end
        soilstate.rootr_patch .= 0.0

        energyflux = CLM.EnergyFluxData()
        CLM.energyflux_init!(energyflux, np, nc, nc, 1)
        energyflux.btran_patch .= 0.0
        energyflux.rresis_patch .= 0.0

        waterstatebulk = CLM.WaterStateBulkData()
        CLM.waterstatebulk_init!(waterstatebulk, nc, np, nc, 1)

        waterdiagbulk = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(waterdiagbulk, nc, np, nc, 1)
        # Set moderate liqvol in soil layers
        for c in 1:nc, j in 1:nlevgrnd
            waterdiagbulk.h2osoi_liqvol_col[c, j + joff] = 0.25
        end

        npft = 2
        smpso = fill(-35000.0, npft)
        smpsc = fill(-275000.0, npft)
        altmax_lastyear_indx = fill(10.0, nc)
        altmax_indx = fill(10.0, nc)

        mask_patch = trues(np)
        bounds_patch = 1:np

        CLM.calc_root_moist_stress!(soilstate, energyflux, temp,
            waterstatebulk, waterdiagbulk, col, patchdata,
            smpso, smpsc, altmax_lastyear_indx, altmax_indx,
            mask_patch, bounds_patch, nlevgrnd, nlevsno)

        for p in 1:np
            @test energyflux.btran_patch[p] > 0.0
            @test energyflux.btran_patch[p] <= 1.0
            @test sum(soilstate.rootr_patch[p, 1:nlevgrnd]) ≈ 1.0 atol=1e-10
        end
    end

    @testset "calc_root_moist_stress! — invalid method errors" begin
        CLM.soil_moist_stress_ctrl.root_moist_stress_method = 99

        nlevtot = nlevsno + nlevgrnd

        col = CLM.ColumnData()
        CLM.column_init!(col, nc)
        col.dz .= 0.1

        patchdata = CLM.PatchData()
        CLM.patch_init!(patchdata, np)
        for p in 1:np
            patchdata.column[p] = p
            patchdata.itype[p] = 1
        end

        temp = CLM.TemperatureData()
        CLM.temperature_init!(temp, np, nc, nc, 1)
        temp.t_soisno_col .= 290.0

        soilstate = CLM.SoilStateData()
        CLM.soilstate_init!(soilstate, np, nc)
        soilstate.watsat_col .= 0.45
        soilstate.sucsat_col .= 100.0
        soilstate.bsw_col .= 5.0
        soilstate.eff_porosity_col .= 0.40
        soilstate.rootfr_patch .= 1.0 / nlevgrnd

        energyflux = CLM.EnergyFluxData()
        CLM.energyflux_init!(energyflux, np, nc, nc, 1)

        waterstatebulk = CLM.WaterStateBulkData()
        CLM.waterstatebulk_init!(waterstatebulk, nc, np, nc, 1)

        waterdiagbulk = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(waterdiagbulk, nc, np, nc, 1)

        smpso = fill(-35000.0, 2)
        smpsc = fill(-275000.0, 2)
        altmax_lastyear_indx = fill(10.0, nc)
        altmax_indx = fill(10.0, nc)

        mask_patch = trues(np)
        bounds_patch = 1:np

        @test_throws ErrorException CLM.calc_root_moist_stress!(soilstate, energyflux, temp,
            waterstatebulk, waterdiagbulk, col, patchdata,
            smpso, smpsc, altmax_lastyear_indx, altmax_indx,
            mask_patch, bounds_patch, nlevgrnd, nlevsno)

        # Restore
        CLM.init_root_moist_stress!()
    end

    # Restore module state
    CLM.soil_moist_stress_ctrl.root_moist_stress_method = saved_method
    CLM.soil_moist_stress_ctrl.perchroot = saved_perchroot
    CLM.soil_moist_stress_ctrl.perchroot_alt = saved_perchroot_alt

    # Restore varpar
    vp.nlevsno = saved_nlevsno
    vp.nlevgrnd = saved_nlevgrnd
    vp.nlevurb = saved_nlevurb
    vp.nlevmaxurbgrnd = saved_nlevmaxurbgrnd
    vp.nlevsoi = saved_nlevsoi
end
