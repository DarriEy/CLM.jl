@testset "Dynamic Init Columns" begin

    # Build a small synthetic gridcell with:
    #   - landunit 1 (ISTSOIL): 1 natural-veg column (c=1), already active,
    #     with populated soil/water state.
    #   - landunit 2 (ISTCROP): 1 crop column (c=2), newly activating
    #     (was inactive in the prior step), with zeroed state.
    # Run template selection + initialize_new_columns!; assert the crop
    # column's below-ground state was copied from the natural-veg template,
    # the snow layers were left untouched, and the already-active column is
    # unchanged.

    function build_subgrid()
        CLM.varpar_init!(CLM.varpar, 1, 15, 2, 5)
        CLM.varcon_init!()

        ng, nl, nc, np = 1, 2, 2, 2

        grc = CLM.GridcellData()
        lun = CLM.LandunitData()
        col = CLM.ColumnData()
        pch = CLM.PatchData()
        CLM.gridcell_init!(grc, ng)
        CLM.landunit_init!(lun, nl)
        CLM.column_init!(col, nc)
        CLM.patch_init!(pch, np)

        li = Ref(0)
        ci = Ref(0)
        pi = Ref(0)

        # Landunit 1 (soil): 1 column, 1 patch (active, weight > 0)
        CLM.add_landunit!(lun, li, 1, CLM.ISTSOIL, 0.6)
        CLM.add_column!(col, lun, ci, li[], CLM.ISTSOIL, 1.0)
        CLM.add_patch!(pch, col, lun, pi, ci[], 1, 1.0)

        # Landunit 2 (crop): 1 column, 1 patch (newly activating)
        CLM.add_landunit!(lun, li, 1, CLM.ISTCROP, 0.4)
        CLM.add_column!(col, lun, ci, li[], CLM.ISTCROP, 1.0)
        CLM.add_patch!(pch, col, lun, pi, ci[], 1, 1.0)

        bounds = CLM.BoundsType(begg=1, endg=ng, begl=1, endl=nl,
                                begc=1, endc=nc, begp=1, endp=np)

        CLM.clm_ptrs_compdown!(bounds, grc, lun, col, pch)

        return bounds, grc, lun, col, pch
    end

    @testset "template_col_from_landunit" begin
        bounds, grc, lun, col, pch = build_subgrid()
        nlevsno = CLM.varpar.nlevsno

        # Prior active flags: only the natural-veg column (c=1) was active.
        cactive_prior = [true, false]

        # Crop column (c=2): natural-veg template should be column 1.
        c_template = CLM.template_col_from_landunit(bounds, 2, CLM.ISTSOIL,
                                                    cactive_prior, grc, lun, col)
        @test c_template == 1

        # No active crop column in the prior step -> TEMPLATE_NONE_FOUND.
        c_none = CLM.template_col_from_landunit(bounds, 2, CLM.ISTCROP,
                                                cactive_prior, grc, lun, col)
        @test c_none == CLM.TEMPLATE_NONE_FOUND

        # A landunit type not present on the gridcell -> TEMPLATE_NONE_FOUND.
        c_absent = CLM.template_col_from_landunit(bounds, 1, CLM.ISTDLAK,
                                                  cactive_prior, grc, lun, col)
        @test c_absent == CLM.TEMPLATE_NONE_FOUND
    end

    @testset "template_col_from_natveg_array!" begin
        bounds, grc, lun, col, pch = build_subgrid()
        cactive_prior = [true, false]
        c_templates = fill(CLM.ISPVAL, bounds.endc)
        CLM.template_col_from_natveg_array!(bounds, cactive_prior, c_templates,
                                            grc, lun, col)
        # Both columns share the same gridcell -> both template off natveg c=1.
        @test c_templates[1] == 1
        @test c_templates[2] == 1
    end

    @testset "initial_template_col_crop / dispatcher" begin
        bounds, grc, lun, col, pch = build_subgrid()
        cactive_prior = [true, false]

        # Crop column gets the natural-veg column as template.
        c_template = CLM.initial_template_col_crop(bounds, 2, cactive_prior,
                                                   grc, lun, col)
        @test c_template == 1

        # Dispatcher for the crop column routes to the same answer.
        c_disp = CLM.initial_template_col_dispatcher(bounds, 2, cactive_prior,
                                                     grc, lun, col)
        @test c_disp == 1

        # If no natveg is active, fall through to the active crop column:
        # mark the crop column itself prior-active so the ISTCROP lookup finds it.
        cactive_none = [false, false]
        cactive_cropactive = [false, true]
        c_cropfallback = CLM.initial_template_col_crop(bounds, 2,
                                                       cactive_cropactive,
                                                       grc, lun, col)
        @test c_cropfallback == 2

        c_nonefound = CLM.initial_template_col_crop(bounds, 2, cactive_none,
                                                    grc, lun, col)
        @test c_nonefound == CLM.TEMPLATE_NONE_FOUND
    end

    @testset "initial_template_col_soil (0-weight virtual col)" begin
        bounds, grc, lun, col, pch = build_subgrid()
        # A natveg column with wtgcell > 0 violates the assumption -> error.
        col.wtgcell[1] = 0.6
        @test_throws ErrorException CLM.initial_template_col_soil(1, col)

        # A 0-weight virtual natveg column returns TEMPLATE_NONE_FOUND.
        col.wtgcell[1] = 0.0
        @test CLM.initial_template_col_soil(1, col) == CLM.TEMPLATE_NONE_FOUND
    end

    @testset "initialize_new_columns! copies state to new crop column" begin
        bounds, grc, lun, col, pch = build_subgrid()
        nlevsno = CLM.varpar.nlevsno
        nc = bounds.endc

        temp = CLM.TemperatureData()
        ws = CLM.WaterStateData()
        CLM.temperature_init!(temp, bounds.endp, nc, bounds.endl, bounds.endg)
        CLM.waterstate_init!(ws, nc, bounds.endp, bounds.endl, bounds.endg)

        nlev_soisno = size(temp.t_soisno_col, 2)
        nlev_vol    = size(ws.h2osoi_vol_col, 2)

        # --- Populate the natural-veg template column (c=1) ---
        for j in 1:nlev_soisno
            temp.t_soisno_col[1, j] = 280.0 + j      # distinct per level
            ws.h2osoi_liq_col[1, j] = 10.0 + j
            ws.h2osoi_ice_col[1, j] = 5.0 + j
        end
        for j in 1:nlev_vol
            ws.h2osoi_vol_col[1, j] = 0.3 + 0.01 * j
        end
        ws.wa_col[1] = 4321.0

        # --- Zero / sentinel the newly-activating crop column (c=2) ---
        for j in 1:nlev_soisno
            temp.t_soisno_col[2, j] = 0.0
            ws.h2osoi_liq_col[2, j] = 0.0
            ws.h2osoi_ice_col[2, j] = 0.0
        end
        for j in 1:nlev_vol
            ws.h2osoi_vol_col[2, j] = 0.0
        end
        ws.wa_col[2] = 0.0

        # Snapshot the template column to confirm it is untouched afterwards.
        t_soisno_tmpl = copy(temp.t_soisno_col[1, :])
        liq_tmpl      = copy(ws.h2osoi_liq_col[1, :])
        wa_tmpl       = ws.wa_col[1]

        # Snapshot the crop column's snow (above-ground) layers — must NOT be
        # overwritten (cold-start snow is preserved).
        snow_t   = copy(temp.t_soisno_col[2, 1:nlevsno])
        snow_liq = copy(ws.h2osoi_liq_col[2, 1:nlevsno])

        # c=1 active prior + now; c=2 newly active now (was inactive prior).
        col.active[1] = true
        col.active[2] = true
        cactive_prior = [true, false]

        CLM.initialize_new_columns!(bounds, cactive_prior, temp, ws,
                                    grc, lun, col)

        soi_lb = nlevsno + 1

        # Below-ground levels copied from template (c=1) into new col (c=2).
        @test temp.t_soisno_col[2, soi_lb:nlev_soisno] ==
              temp.t_soisno_col[1, soi_lb:nlev_soisno]
        @test ws.h2osoi_liq_col[2, soi_lb:nlev_soisno] ==
              ws.h2osoi_liq_col[1, soi_lb:nlev_soisno]
        @test ws.h2osoi_ice_col[2, soi_lb:nlev_soisno] ==
              ws.h2osoi_ice_col[1, soi_lb:nlev_soisno]
        @test ws.h2osoi_vol_col[2, :] == ws.h2osoi_vol_col[1, :]
        @test ws.wa_col[2] == ws.wa_col[1]
        @test ws.wa_col[2] == 4321.0

        # Snow (above-ground) layers of the new column NOT overwritten.
        @test temp.t_soisno_col[2, 1:nlevsno] == snow_t
        @test ws.h2osoi_liq_col[2, 1:nlevsno] == snow_liq
        @test all(temp.t_soisno_col[2, 1:nlevsno] .== 0.0)

        # Template column (already active) is untouched.
        @test temp.t_soisno_col[1, :] == t_soisno_tmpl
        @test ws.h2osoi_liq_col[1, :] == liq_tmpl
        @test ws.wa_col[1] == wa_tmpl
    end

    @testset "already-active column is untouched" begin
        # If a column was active in the prior step, it is not re-initialized.
        bounds, grc, lun, col, pch = build_subgrid()
        nlevsno = CLM.varpar.nlevsno
        nc = bounds.endc

        temp = CLM.TemperatureData()
        ws = CLM.WaterStateData()
        CLM.temperature_init!(temp, bounds.endp, nc, bounds.endl, bounds.endg)
        CLM.waterstate_init!(ws, nc, bounds.endp, bounds.endl, bounds.endg)

        nlev_soisno = size(temp.t_soisno_col, 2)
        for j in 1:nlev_soisno
            temp.t_soisno_col[1, j] = 280.0 + j
            temp.t_soisno_col[2, j] = 111.0 + j   # distinct, would change if copied
        end
        ws.wa_col[1] = 1.0
        ws.wa_col[2] = 2.0

        crop_before = copy(temp.t_soisno_col[2, :])
        wa2_before  = ws.wa_col[2]

        # Both columns were active in the prior step -> nothing newly active.
        col.active[1] = true
        col.active[2] = true
        cactive_prior = [true, true]

        CLM.initialize_new_columns!(bounds, cactive_prior, temp, ws,
                                    grc, lun, col)

        @test temp.t_soisno_col[2, :] == crop_before
        @test ws.wa_col[2] == wa2_before
    end

end
