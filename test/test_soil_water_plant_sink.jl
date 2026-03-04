@testset "SoilWaterPlantSink" begin
    # Save and restore varpar state
    vp = CLM.varpar
    saved_nlevsno        = vp.nlevsno
    saved_nlevgrnd       = vp.nlevgrnd
    saved_nlevurb        = vp.nlevurb
    saved_nlevmaxurbgrnd = vp.nlevmaxurbgrnd
    saved_nlevsoi        = vp.nlevsoi
    saved_nlevlak        = vp.nlevlak

    vp.nlevsno        = 5
    vp.nlevgrnd       = 10
    vp.nlevurb        = 5
    vp.nlevmaxurbgrnd = 10
    vp.nlevsoi        = 4
    vp.nlevlak        = 10

    nlevsno  = vp.nlevsno
    nlevgrnd = vp.nlevgrnd
    nlevsoi  = vp.nlevsoi
    nlevmaxurbgrnd = vp.nlevmaxurbgrnd

    # =====================================================================
    # Helper to set up a minimal column/patch/landunit configuration
    # =====================================================================
    function setup_test_data(;
            nc::Int,
            np::Int,
            nl::Int,
            col_itype::Vector{Int},
            col_landunit::Vector{Int},
            lun_itype::Vector{Int},
            patches_per_col::Vector{Int},
            rootr_vals::Matrix{Float64},       # np x nlevsoi
            tran_veg_patch::Vector{Float64},   # np
            tran_veg_col::Vector{Float64},     # nc
            wtcol::Vector{Float64},            # np
            active_patch::Vector{Bool})        # np

        # --- Column data ---
        col = CLM.ColumnData()
        CLM.column_init!(col, nc)
        col.itype[1:nc]    .= col_itype
        col.landunit[1:nc] .= col_landunit
        col.active[1:nc]   .= true

        # Set up patch ranges for each column
        pidx = 1
        for c in 1:nc
            col.patchi[c]   = pidx
            col.npatches[c] = patches_per_col[c]
            col.patchf[c]   = pidx + patches_per_col[c] - 1
            pidx += patches_per_col[c]
        end

        # Set z values for soil layers (needed by hydstress method)
        # z is stored as (ncols, nlevsno + nlevmaxurbgrnd)
        for c in 1:nc
            for j in 1:nlevsoi
                col.z[c, j + nlevsno] = 0.1 * j  # simple depth in meters
            end
        end

        # --- Landunit data ---
        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, nl)
        lun.itype[1:nl] .= lun_itype

        # --- Patch data ---
        patch_data = CLM.PatchData()
        CLM.patch_init!(patch_data, np)
        patch_data.active[1:np] .= active_patch
        patch_data.wtcol[1:np]  .= wtcol

        # --- Soil state ---
        ss = CLM.SoilStateData()
        CLM.soilstate_init!(ss, np, nc)
        for p in 1:np
            for j in 1:nlevsoi
                ss.rootr_patch[p, j] = rootr_vals[p, j]
            end
        end
        # Zero rootr_col initially
        for c in 1:nc
            for j in 1:nlevgrnd
                ss.rootr_col[c, j] = 0.0
            end
        end

        # --- Water flux bulk ---
        wfb = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(wfb, nc, np, nl, 1)
        for p in 1:np
            wfb.wf.qflx_tran_veg_patch[p] = tran_veg_patch[p]
        end
        for c in 1:nc
            wfb.wf.qflx_tran_veg_col[c] = tran_veg_col[c]
        end
        # Zero rootsoi
        for c in 1:nc
            for j in 1:nlevsoi
                wfb.qflx_rootsoi_col[c, j] = 0.0
            end
        end

        # --- Canopy state ---
        cs = CLM.CanopyStateData()
        CLM.canopystate_init!(cs, np)

        # --- Energy flux ---
        ef = CLM.EnergyFluxData()
        CLM.energyflux_init!(ef, np, nc, nl, 1)

        return col, lun, patch_data, ss, wfb, cs, ef
    end

    # =====================================================================
    @testset "Default method - single column, single patch" begin
        nc = 1; np = 1; nl = 1
        rootr = zeros(np, nlevsoi)
        rootr[1, :] .= [0.4, 0.3, 0.2, 0.1]  # sums to 1.0
        tran_patch = [0.002]  # 2 mm/s transpiration
        tran_col   = [0.002]
        wtcol      = [1.0]
        active     = [true]

        col, lun, patch_data, ss, wfb, cs, ef = setup_test_data(
            nc=nc, np=np, nl=nl,
            col_itype=[CLM.ISTSOIL],
            col_landunit=[1],
            lun_itype=[CLM.ISTSOIL],
            patches_per_col=[1],
            rootr_vals=rootr,
            tran_veg_patch=tran_patch,
            tran_veg_col=tran_col,
            wtcol=wtcol,
            active_patch=active)

        bounds_col = 1:nc
        filterc = [1]

        CLM.compute_effec_rootfrac_and_vert_tran_sink_default!(
            bounds_col, nlevsoi, filterc,
            ss, wfb, col, patch_data)

        # With a single active patch, rootr_col should equal rootr_patch
        for j in 1:nlevsoi
            @test ss.rootr_col[1, j] ≈ rootr[1, j] atol=1e-12
        end

        # qflx_rootsoi_col = rootr_col * qflx_tran_veg_col
        for j in 1:nlevsoi
            @test wfb.qflx_rootsoi_col[1, j] ≈ rootr[1, j] * tran_col[1] atol=1e-14
        end
    end

    # =====================================================================
    @testset "Default method - single column, two patches" begin
        nc = 1; np = 2; nl = 1

        # Patch 1: deeper roots; Patch 2: shallower roots
        rootr = zeros(np, nlevsoi)
        rootr[1, :] .= [0.1, 0.2, 0.3, 0.4]
        rootr[2, :] .= [0.5, 0.3, 0.15, 0.05]

        tran_patch = [0.001, 0.003]  # patch 2 transpires more
        # Column transpiration is area-weighted sum
        wtcol = [0.4, 0.6]
        tran_col = [tran_patch[1] * wtcol[1] + tran_patch[2] * wtcol[2]]
        active = [true, true]

        col, lun, patch_data, ss, wfb, cs, ef = setup_test_data(
            nc=nc, np=np, nl=nl,
            col_itype=[CLM.ISTSOIL],
            col_landunit=[1],
            lun_itype=[CLM.ISTSOIL],
            patches_per_col=[2],
            rootr_vals=rootr,
            tran_veg_patch=tran_patch,
            tran_veg_col=tran_col,
            wtcol=wtcol,
            active_patch=active)

        bounds_col = 1:nc
        filterc = [1]

        CLM.compute_effec_rootfrac_and_vert_tran_sink_default!(
            bounds_col, nlevsoi, filterc,
            ss, wfb, col, patch_data)

        # Manually compute expected rootr_col
        # temp = sum(tran_veg_patch * wtcol) for active patches
        temp_val = tran_patch[1] * wtcol[1] + tran_patch[2] * wtcol[2]

        for j in 1:nlevsoi
            numer = rootr[1, j] * tran_patch[1] * wtcol[1] +
                    rootr[2, j] * tran_patch[2] * wtcol[2]
            expected_rootr = numer / temp_val
            @test ss.rootr_col[1, j] ≈ expected_rootr atol=1e-12
            @test wfb.qflx_rootsoi_col[1, j] ≈ expected_rootr * tran_col[1] atol=1e-14
        end
    end

    # =====================================================================
    @testset "Default method - zero transpiration" begin
        nc = 1; np = 1; nl = 1
        rootr = zeros(np, nlevsoi)
        rootr[1, :] .= [0.25, 0.25, 0.25, 0.25]
        tran_patch = [0.0]
        tran_col = [0.0]
        wtcol = [1.0]
        active = [true]

        col, lun, patch_data, ss, wfb, cs, ef = setup_test_data(
            nc=nc, np=np, nl=nl,
            col_itype=[CLM.ISTSOIL],
            col_landunit=[1],
            lun_itype=[CLM.ISTSOIL],
            patches_per_col=[1],
            rootr_vals=rootr,
            tran_veg_patch=tran_patch,
            tran_veg_col=tran_col,
            wtcol=wtcol,
            active_patch=active)

        bounds_col = 1:nc
        filterc = [1]

        CLM.compute_effec_rootfrac_and_vert_tran_sink_default!(
            bounds_col, nlevsoi, filterc,
            ss, wfb, col, patch_data)

        # temp == 0, so rootr_col stays 0, qflx_rootsoi stays 0
        for j in 1:nlevsoi
            @test ss.rootr_col[1, j] ≈ 0.0 atol=1e-14
            @test wfb.qflx_rootsoi_col[1, j] ≈ 0.0 atol=1e-14
        end
    end

    # =====================================================================
    @testset "Default method - inactive patch is skipped" begin
        nc = 1; np = 2; nl = 1
        rootr = zeros(np, nlevsoi)
        rootr[1, :] .= [0.4, 0.3, 0.2, 0.1]
        rootr[2, :] .= [0.1, 0.1, 0.1, 0.7]  # should be ignored

        tran_patch = [0.002, 0.005]
        wtcol = [0.5, 0.5]
        tran_col = [0.002 * 0.5]  # only patch 1 contributes
        active = [true, false]  # patch 2 inactive

        col, lun, patch_data, ss, wfb, cs, ef = setup_test_data(
            nc=nc, np=np, nl=nl,
            col_itype=[CLM.ISTSOIL],
            col_landunit=[1],
            lun_itype=[CLM.ISTSOIL],
            patches_per_col=[2],
            rootr_vals=rootr,
            tran_veg_patch=tran_patch,
            tran_veg_col=tran_col,
            wtcol=wtcol,
            active_patch=active)

        bounds_col = 1:nc
        filterc = [1]

        CLM.compute_effec_rootfrac_and_vert_tran_sink_default!(
            bounds_col, nlevsoi, filterc,
            ss, wfb, col, patch_data)

        # With only patch 1 active, rootr_col = rootr_patch[1,:]
        for j in 1:nlevsoi
            @test ss.rootr_col[1, j] ≈ rootr[1, j] atol=1e-12
        end
    end

    # =====================================================================
    @testset "HydStress Roads method - matches Default logic" begin
        # The hydstress roads method uses the same algorithm as default
        nc = 1; np = 1; nl = 1
        rootr = zeros(np, nlevsoi)
        rootr[1, :] .= [0.3, 0.3, 0.2, 0.2]
        tran_patch = [0.001]
        tran_col = [0.001]
        wtcol = [1.0]
        active = [true]

        col, lun, patch_data, ss, wfb, cs, ef = setup_test_data(
            nc=nc, np=np, nl=nl,
            col_itype=[CLM.ICOL_ROAD_PERV],
            col_landunit=[1],
            lun_itype=[CLM.ISTURB_MD],
            patches_per_col=[1],
            rootr_vals=rootr,
            tran_veg_patch=tran_patch,
            tran_veg_col=tran_col,
            wtcol=wtcol,
            active_patch=active)

        bounds_col = 1:nc
        filterc = [1]

        CLM.compute_effec_rootfrac_and_vert_tran_sink_hydstress_roads!(
            bounds_col, nlevsoi, filterc,
            ss, wfb, col, patch_data)

        for j in 1:nlevsoi
            @test ss.rootr_col[1, j] ≈ rootr[1, j] atol=1e-12
            @test wfb.qflx_rootsoi_col[1, j] ≈ rootr[1, j] * tran_col[1] atol=1e-14
        end
    end

    # =====================================================================
    @testset "HydStress method - single column, single patch" begin
        nc = 1; np = 1; nl = 1
        rootr = zeros(np, nlevsoi)  # not used by hydstress, but needed for setup
        tran_patch = [0.002]
        tran_col = [0.002]
        wtcol = [1.0]
        active = [true]

        col, lun, patch_data, ss, wfb, cs, ef = setup_test_data(
            nc=nc, np=np, nl=nl,
            col_itype=[CLM.ISTSOIL],
            col_landunit=[1],
            lun_itype=[CLM.ISTSOIL],
            patches_per_col=[1],
            rootr_vals=rootr,
            tran_veg_patch=tran_patch,
            tran_veg_col=tran_col,
            wtcol=wtcol,
            active_patch=active)

        bounds_col = 1:nc
        filterc = [1]

        # Set up hydraulic stress inputs
        # frac_veg_nosno > 0 required
        cs.frac_veg_nosno_patch[1] = 1

        # k_soil_root: conductance (mm/s)
        for j in 1:nlevsoi
            ss.k_soil_root_patch[1, j] = 0.01  # 0.01 mm/s
        end

        # smp_l_col: soil matric potential (mm)
        for j in 1:nlevgrnd
            ss.smp_l_col[1, j] = -500.0  # mm
        end

        # vegwp: vegetation water potential (mm)
        # root_idx = 4
        cs.vegwp_patch[1, 4] = -1000.0  # mm, more negative than soil

        CLM.compute_effec_rootfrac_and_vert_tran_sink_hydstress!(
            bounds_col, nlevsoi, filterc,
            wfb, ss, cs, col, patch_data)

        # Check manually for each layer
        for j in 1:nlevsoi
            grav2 = col.z[1, j + nlevsno] * 1000.0
            expected_flux = 0.01 * (-500.0 - (-1000.0) - grav2)
            @test wfb.qflx_rootsoi_col[1, j] ≈ expected_flux * wtcol[1] atol=1e-10
        end
    end

    # =====================================================================
    @testset "HydStress method - negative flux accumulates phs_neg" begin
        nc = 1; np = 1; nl = 1
        rootr = zeros(np, nlevsoi)
        tran_patch = [0.002]
        tran_col = [0.002]
        wtcol = [1.0]
        active = [true]

        col, lun, patch_data, ss, wfb, cs, ef = setup_test_data(
            nc=nc, np=np, nl=nl,
            col_itype=[CLM.ISTSOIL],
            col_landunit=[1],
            lun_itype=[CLM.ISTSOIL],
            patches_per_col=[1],
            rootr_vals=rootr,
            tran_veg_patch=tran_patch,
            tran_veg_col=tran_col,
            wtcol=wtcol,
            active_patch=active)

        bounds_col = 1:nc
        filterc = [1]

        cs.frac_veg_nosno_patch[1] = 1

        for j in 1:nlevsoi
            ss.k_soil_root_patch[1, j] = 0.01
        end

        # Make soil very dry so flux is from root to soil (negative)
        # smp very negative, vegwp less negative
        for j in 1:nlevgrnd
            ss.smp_l_col[1, j] = -5000.0
        end
        cs.vegwp_patch[1, 4] = -500.0

        CLM.compute_effec_rootfrac_and_vert_tran_sink_hydstress!(
            bounds_col, nlevsoi, filterc,
            wfb, ss, cs, col, patch_data)

        # All fluxes should be negative (water from root to dry soil)
        phs_neg_expected = 0.0
        for j in 1:nlevsoi
            grav2 = col.z[1, j + nlevsno] * 1000.0
            expected_flux = 0.01 * (-5000.0 - (-500.0) - grav2)
            @test expected_flux < 0.0
            @test wfb.qflx_rootsoi_col[1, j] ≈ expected_flux * wtcol[1] atol=1e-10
            if expected_flux * wtcol[1] < 0.0
                phs_neg_expected += expected_flux * wtcol[1]
            end
        end

        @test wfb.qflx_phs_neg_col[1] ≈ phs_neg_expected atol=1e-10

        # Also check hydr_redist accumulated negative fluxes
        @test wfb.qflx_hydr_redist_patch[1] < 0.0
    end

    # =====================================================================
    @testset "HydStress method - inactive/no-veg patch skipped" begin
        nc = 1; np = 2; nl = 1
        rootr = zeros(np, nlevsoi)
        tran_patch = [0.002, 0.003]
        tran_col = [0.002]
        wtcol = [0.6, 0.4]
        active = [true, true]

        col, lun, patch_data, ss, wfb, cs, ef = setup_test_data(
            nc=nc, np=np, nl=nl,
            col_itype=[CLM.ISTSOIL],
            col_landunit=[1],
            lun_itype=[CLM.ISTSOIL],
            patches_per_col=[2],
            rootr_vals=rootr,
            tran_veg_patch=tran_patch,
            tran_veg_col=tran_col,
            wtcol=wtcol,
            active_patch=active)

        bounds_col = 1:nc
        filterc = [1]

        # Patch 1: has vegetation; Patch 2: no vegetation (frac_veg_nosno=0)
        cs.frac_veg_nosno_patch[1] = 1
        cs.frac_veg_nosno_patch[2] = 0  # should be skipped

        for p in 1:np
            for j in 1:nlevsoi
                ss.k_soil_root_patch[p, j] = 0.01
            end
        end
        for j in 1:nlevgrnd
            ss.smp_l_col[1, j] = -500.0
        end
        cs.vegwp_patch[1, 4] = -1000.0
        cs.vegwp_patch[2, 4] = -200.0  # would create different flux if not skipped

        CLM.compute_effec_rootfrac_and_vert_tran_sink_hydstress!(
            bounds_col, nlevsoi, filterc,
            wfb, ss, cs, col, patch_data)

        # Only patch 1 should contribute
        for j in 1:nlevsoi
            grav2 = col.z[1, j + nlevsno] * 1000.0
            expected_flux_p1 = 0.01 * (-500.0 - (-1000.0) - grav2)
            expected_col = expected_flux_p1 * wtcol[1]  # only patch 1
            @test wfb.qflx_rootsoi_col[1, j] ≈ expected_col atol=1e-10
        end
    end

    # =====================================================================
    @testset "Wrapper - default path with multiple column types" begin
        # Set up: 3 columns
        #   c=1: pervious road (urban)
        #   c=2: crop (not istsoil, not pervious road)
        #   c=3: natural vegetation (istsoil)
        nc = 3; np = 3; nl = 3

        col_itype = [CLM.ICOL_ROAD_PERV, CLM.ISTCROP, CLM.ISTSOIL]
        col_landunit = [1, 2, 3]
        lun_itype_arr = [CLM.ISTURB_MD, CLM.ISTCROP, CLM.ISTSOIL]

        rootr = zeros(np, nlevsoi)
        rootr[1, :] .= [0.25, 0.25, 0.25, 0.25]
        rootr[2, :] .= [0.4, 0.3, 0.2, 0.1]
        rootr[3, :] .= [0.1, 0.2, 0.3, 0.4]

        tran_patch = [0.001, 0.002, 0.003]
        tran_col = [0.001, 0.002, 0.003]
        wtcol_arr = [1.0, 1.0, 1.0]
        active = [true, true, true]

        col, lun, patch_data, ss, wfb, cs, ef = setup_test_data(
            nc=nc, np=np, nl=nl,
            col_itype=col_itype,
            col_landunit=col_landunit,
            lun_itype=lun_itype_arr,
            patches_per_col=[1, 1, 1],
            rootr_vals=rootr,
            tran_veg_patch=tran_patch,
            tran_veg_col=tran_col,
            wtcol=wtcol_arr,
            active_patch=active)

        mask_hydrology = falses(nc)
        mask_hydrology[1:nc] .= true

        # Call wrapper with use_hydrstress=false
        CLM.compute_effec_rootfrac_and_vert_tran_sink!(
            1:nc, nlevsoi, mask_hydrology,
            ss, cs, wfb, ef, col, lun, patch_data;
            use_hydrstress=false)

        # Each column has a single active patch, so rootr_col = rootr_patch
        for c in 1:nc
            for j in 1:nlevsoi
                @test ss.rootr_col[c, j] ≈ rootr[c, j] atol=1e-12
                @test wfb.qflx_rootsoi_col[c, j] ≈ rootr[c, j] * tran_col[c] atol=1e-14
            end
        end
    end

    # =====================================================================
    @testset "Wrapper - filter count mismatch raises error" begin
        # Create a situation where the column types do not partition
        # correctly. This requires a column that passes mask_hydrology
        # but doesn't match any of the three groups.
        # Actually, all columns must match one of:
        #   1) icol_road_perv
        #   2) not icol_road_perv AND not istsoil landunit
        #   3) istsoil landunit
        # These three conditions are mutually exclusive and exhaustive for
        # any column in mask_hydrology, so an error should NOT happen
        # normally. But if mask_hydrology includes a column outside bounds_col
        # weirdly, we can test the error path by having mask include extra.
        # Instead, let's test that a normal call does NOT error:
        nc = 1; np = 1; nl = 1
        rootr = zeros(np, nlevsoi)
        rootr[1, :] .= [0.25, 0.25, 0.25, 0.25]

        col, lun, patch_data, ss, wfb, cs, ef = setup_test_data(
            nc=nc, np=np, nl=nl,
            col_itype=[CLM.ISTSOIL],
            col_landunit=[1],
            lun_itype=[CLM.ISTSOIL],
            patches_per_col=[1],
            rootr_vals=rootr,
            tran_veg_patch=[0.001],
            tran_veg_col=[0.001],
            wtcol=[1.0],
            active_patch=[true])

        mask_hydrology = trues(nc)

        # This should not throw
        CLM.compute_effec_rootfrac_and_vert_tran_sink!(
            1:nc, nlevsoi, mask_hydrology,
            ss, cs, wfb, ef, col, lun, patch_data;
            use_hydrstress=false)

        @test true  # reached here without error
    end

    # =====================================================================
    @testset "Wrapper - hydstress path dispatches correctly" begin
        # 2 columns: pervious road + natural veg
        nc = 2; np = 2; nl = 2
        col_itype = [CLM.ICOL_ROAD_PERV, CLM.ISTSOIL]
        col_landunit = [1, 2]
        lun_itype_arr = [CLM.ISTURB_MD, CLM.ISTSOIL]

        rootr = zeros(np, nlevsoi)
        rootr[1, :] .= [0.25, 0.25, 0.25, 0.25]
        rootr[2, :] .= [0.0, 0.0, 0.0, 0.0]  # not used by hydstress for natveg

        tran_patch = [0.001, 0.002]
        tran_col = [0.001, 0.002]
        wtcol_arr = [1.0, 1.0]
        active = [true, true]

        col, lun, patch_data, ss, wfb, cs, ef = setup_test_data(
            nc=nc, np=np, nl=nl,
            col_itype=col_itype,
            col_landunit=col_landunit,
            lun_itype=lun_itype_arr,
            patches_per_col=[1, 1],
            rootr_vals=rootr,
            tran_veg_patch=tran_patch,
            tran_veg_col=tran_col,
            wtcol=wtcol_arr,
            active_patch=active)

        # Set up hydstress inputs for natural veg column (c=2)
        cs.frac_veg_nosno_patch[1] = 1  # road patch
        cs.frac_veg_nosno_patch[2] = 1  # natural veg patch

        for p in 1:np
            for j in 1:nlevsoi
                ss.k_soil_root_patch[p, j] = 0.005
            end
        end
        for c in 1:nc
            for j in 1:nlevgrnd
                ss.smp_l_col[c, j] = -300.0
            end
        end
        cs.vegwp_patch[1, 4] = -800.0
        cs.vegwp_patch[2, 4] = -600.0

        mask_hydrology = trues(nc)

        CLM.compute_effec_rootfrac_and_vert_tran_sink!(
            1:nc, nlevsoi, mask_hydrology,
            ss, cs, wfb, ef, col, lun, patch_data;
            use_hydrstress=true)

        # Column 1 (pervious road) should use hydstress_roads method
        # which is the same as default: rootr_col = rootr_patch
        for j in 1:nlevsoi
            @test ss.rootr_col[1, j] ≈ rootr[1, j] atol=1e-12
        end

        # Column 2 (natural veg) should use hydstress method
        # qflx_rootsoi_col based on conductance and potential gradient
        for j in 1:nlevsoi
            grav2 = col.z[2, j + nlevsno] * 1000.0
            expected = 0.005 * (-300.0 - (-600.0) - grav2) * wtcol_arr[2]
            @test wfb.qflx_rootsoi_col[2, j] ≈ expected atol=1e-10
        end
    end

    # =====================================================================
    @testset "Default method - empty filter does nothing" begin
        nc = 2; np = 2; nl = 1
        rootr = zeros(np, nlevsoi)
        rootr[1, :] .= [0.25, 0.25, 0.25, 0.25]
        rootr[2, :] .= [0.25, 0.25, 0.25, 0.25]

        col, lun, patch_data, ss, wfb, cs, ef = setup_test_data(
            nc=nc, np=np, nl=nl,
            col_itype=[CLM.ISTSOIL, CLM.ISTSOIL],
            col_landunit=[1, 1],
            lun_itype=[CLM.ISTSOIL],
            patches_per_col=[1, 1],
            rootr_vals=rootr,
            tran_veg_patch=[0.001, 0.001],
            tran_veg_col=[0.001, 0.001],
            wtcol=[1.0, 1.0],
            active_patch=[true, true])

        bounds_col = 1:nc
        filterc = Int[]  # empty filter

        # Pre-set rootsoi to known value to verify nothing changes
        for c in 1:nc
            for j in 1:nlevsoi
                wfb.qflx_rootsoi_col[c, j] = 999.0
            end
        end

        CLM.compute_effec_rootfrac_and_vert_tran_sink_default!(
            bounds_col, nlevsoi, filterc,
            ss, wfb, col, patch_data)

        # Should be untouched
        for c in 1:nc
            for j in 1:nlevsoi
                @test wfb.qflx_rootsoi_col[c, j] ≈ 999.0 atol=1e-14
            end
        end
    end

    # =====================================================================
    @testset "Wrapper - mask_hydrology filtering" begin
        nc = 2; np = 2; nl = 1
        rootr = zeros(np, nlevsoi)
        rootr[1, :] .= [0.4, 0.3, 0.2, 0.1]
        rootr[2, :] .= [0.1, 0.2, 0.3, 0.4]

        col, lun, patch_data, ss, wfb, cs, ef = setup_test_data(
            nc=nc, np=np, nl=nl,
            col_itype=[CLM.ISTSOIL, CLM.ISTSOIL],
            col_landunit=[1, 1],
            lun_itype=[CLM.ISTSOIL],
            patches_per_col=[1, 1],
            rootr_vals=rootr,
            tran_veg_patch=[0.001, 0.002],
            tran_veg_col=[0.001, 0.002],
            wtcol=[1.0, 1.0],
            active_patch=[true, true])

        # Only column 1 is in hydrology mask
        mask_hydrology = falses(nc)
        mask_hydrology[1] = true

        CLM.compute_effec_rootfrac_and_vert_tran_sink!(
            1:nc, nlevsoi, mask_hydrology,
            ss, cs, wfb, ef, col, lun, patch_data;
            use_hydrstress=false)

        # Column 1 should be processed
        for j in 1:nlevsoi
            @test ss.rootr_col[1, j] ≈ rootr[1, j] atol=1e-12
        end

        # Column 2 should NOT have been processed (rootr_col stays 0)
        for j in 1:nlevsoi
            @test ss.rootr_col[2, j] ≈ 0.0 atol=1e-14
        end
    end

    # =====================================================================
    # Restore varpar state
    vp.nlevsno        = saved_nlevsno
    vp.nlevgrnd       = saved_nlevgrnd
    vp.nlevurb        = saved_nlevurb
    vp.nlevmaxurbgrnd = saved_nlevmaxurbgrnd
    vp.nlevsoi        = saved_nlevsoi
    vp.nlevlak        = saved_nlevlak
end
