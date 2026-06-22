@testset "Dynamic Column State Updater" begin

    # ---------------------------------------------------------------------
    # Build a small synthetic single-gridcell subgrid. We override
    # col.wtgcell / col.active directly per test, since these are the
    # quantities the updater snapshots.
    #
    #   landunit 1 (ISTSOIL, natveg)  -> column 1
    #   landunit 2 (ISTSOIL, natveg)  -> column 2   (a second natveg col)
    #   landunit 3 (ISTDLAK, special) -> column 3
    # All on the same gridcell (g=1).
    # ---------------------------------------------------------------------
    function build_csu_subgrid()
        CLM.varpar_init!(CLM.varpar, 1, 15, 2, 5)
        CLM.varcon_init!()

        ng, nl, nc, np = 1, 2, 3, 3

        grc = CLM.GridcellData()
        lun = CLM.LandunitData()
        col = CLM.ColumnData()
        pch = CLM.PatchData()
        CLM.gridcell_init!(grc, ng)
        CLM.landunit_init!(lun, nl)
        CLM.column_init!(col, nc)
        CLM.patch_init!(pch, np)

        li = Ref(0); ci = Ref(0); pi = Ref(0)

        # NOTE: clm_ptrs_compdown! requires each landunit type to be unique
        # per gridcell, so the two natveg columns live on ONE soil landunit.
        # Landunit 1 (soil): 2 columns
        CLM.add_landunit!(lun, li, 1, CLM.ISTSOIL, 0.6)
        CLM.add_column!(col, lun, ci, li[], CLM.ISTSOIL, 0.5)
        CLM.add_patch!(pch, col, lun, pi, ci[], 1, 1.0)
        CLM.add_column!(col, lun, ci, li[], CLM.ISTSOIL, 0.5)
        CLM.add_patch!(pch, col, lun, pi, ci[], 1, 1.0)

        # Landunit 2 (deep lake, special): 1 column
        CLM.add_landunit!(lun, li, 1, CLM.ISTDLAK, 0.4)
        CLM.add_column!(col, lun, ci, li[], CLM.ISTDLAK, 1.0)
        CLM.add_patch!(pch, col, lun, pi, ci[], 1, 1.0)

        bounds = CLM.BoundsType(begg=1, endg=ng, begl=1, endl=nl,
                                begc=1, endc=nc, begp=1, endp=np,
                                level=CLM.BOUNDS_LEVEL_PROC)

        CLM.clm_ptrs_compdown!(bounds, grc, lun, col, pch)

        return bounds, grc, lun, col, pch
    end

    # ---------------------------------------------------------------------
    # Constructor + set_old_weights! / set_new_weights!
    # ---------------------------------------------------------------------
    @testset "constructor + set_old/new_weights!" begin
        bounds, grc, lun, col, pch = build_csu_subgrid()
        nclumps = 1

        csu = CLM.ColumnStateUpdater(bounds, nclumps)
        @test length(csu.cwtgcell_old) == bounds.endc
        @test length(csu.any_changes) == nclumps
        @test all(csu.natveg_template_col .== CLM.TEMPLATE_NONE_FOUND)
        @test csu.any_changes[1] == false

        # Old weights: c1=0.3, c2=0.3, c3=0.4 (cols active in prior step).
        col.wtgcell .= [0.3, 0.3, 0.4]
        col.active  .= [true, true, true]
        CLM.set_old_weights!(csu, bounds, grc, lun, col)
        @test csu.cwtgcell_old == [0.3, 0.3, 0.4]
        # Both columns share gridcell 1; first active natveg col is c=1.
        @test csu.natveg_template_col == [1, 1, 1]

        # No change yet: new == old.
        col.wtgcell .= [0.3, 0.3, 0.4]
        CLM.set_new_weights!(csu, bounds, 1, col)
        @test csu.any_changes[1] == false
        @test all(csu.area_gained_col .== 0.0)

        # Now change: c1 shrinks, c2 grows.
        col.wtgcell .= [0.2, 0.4, 0.4]
        CLM.set_new_weights!(csu, bounds, 1, col)
        @test csu.any_changes[1] == true
        @test csu.area_gained_col ≈ [-0.1, 0.1, 0.0]
    end

    # ---------------------------------------------------------------------
    # Variant 1: no special handling — mass conservation
    # ---------------------------------------------------------------------
    @testset "variant 1: no_special_handling mass conservation" begin
        bounds, grc, lun, col, pch = build_csu_subgrid()
        csu = CLM.ColumnStateUpdater(bounds, 1)

        # Old weights, all active.
        col.active .= [true, true, true]
        col.wtgcell .= [0.3, 0.3, 0.4]
        CLM.set_old_weights!(csu, bounds, grc, lun, col)

        # c1 shrinks 0.3->0.1, c2 grows 0.3->0.5, c3 unchanged.
        col.wtgcell .= [0.1, 0.5, 0.4]
        CLM.set_new_weights!(csu, bounds, 1, col)

        # State variable (mass per unit area). Track total mass = sum(wt*var).
        var = [10.0, 4.0, 7.0]
        mass_old = sum(csu.cwtgcell_old .* var)

        adjustment = zeros(Float64, bounds.endc)
        CLM.update_column_state_no_special_handling!(csu, bounds, 1, col, var;
            adjustment=adjustment)

        mass_new = sum(csu.cwtgcell_new .* var)
        @test mass_new ≈ mass_old atol=1e-12

        # Shrinking & unchanged columns keep their per-area value; their
        # adjustment is 0. Only the growing column gets a nonzero adjustment.
        @test adjustment[1] == 0.0
        @test adjustment[3] == 0.0
        @test var[1] == 10.0          # shrinking col value unchanged
        @test var[3] == 7.0           # unchanged col value unchanged
        @test adjustment[2] != 0.0

        # Growing-column value: pooled loss redistributed.
        #   total loss = 0.2*10 = 2.0; total area lost = 0.2; gain/area = 10.0
        #   mass_gained(c2) = 0.2*10.0 = 2.0
        #   var[2] = (0.3*4.0 + 2.0)/0.5 = 3.2/0.5 = 6.4
        @test var[2] ≈ 6.4
    end

    # ---------------------------------------------------------------------
    # Variant 1: a column growing into newly-available gridcell area
    # (here, the special lake shrinks and a natveg column absorbs that area).
    # ---------------------------------------------------------------------
    @testset "variant 1: column grows into vacated area" begin
        bounds, grc, lun, col, pch = build_csu_subgrid()
        csu = CLM.ColumnStateUpdater(bounds, 1)

        col.active .= [true, true, true]
        col.wtgcell .= [0.3, 0.3, 0.4]
        CLM.set_old_weights!(csu, bounds, grc, lun, col)

        # Lake (c3) shrinks 0.4->0.2; natveg c1 grows 0.3->0.5.
        col.wtgcell .= [0.5, 0.3, 0.2]
        CLM.set_new_weights!(csu, bounds, 1, col)

        var = [10.0, 4.0, 7.0]
        mass_old = sum(csu.cwtgcell_old .* var)
        CLM.update_column_state_no_special_handling!(csu, bounds, 1, col, var)
        mass_new = sum(csu.cwtgcell_new .* var)
        @test mass_new ≈ mass_old atol=1e-12
        # c1 grew, absorbing lake mass:
        #   loss = 0.2*7 = 1.4; area lost = 0.2; gain/area = 7.0
        #   mass_gained = 0.2*7.0 = 1.4
        #   var[1] = (0.3*10 + 1.4)/0.5 = 4.4/0.5 = 8.8
        @test var[1] ≈ 8.8
    end

    # ---------------------------------------------------------------------
    # Variant 1: fractional areas
    # ---------------------------------------------------------------------
    @testset "variant 1: fractional areas" begin
        bounds, grc, lun, col, pch = build_csu_subgrid()
        csu = CLM.ColumnStateUpdater(bounds, 1)

        col.active .= [true, true, true]
        col.wtgcell .= [0.3, 0.3, 0.4]
        CLM.set_old_weights!(csu, bounds, grc, lun, col)
        col.wtgcell .= [0.1, 0.5, 0.4]
        CLM.set_new_weights!(csu, bounds, 1, col)

        var = [10.0, 4.0, 7.0]
        frac_old = [0.5, 0.5, 1.0]
        frac_new = [0.5, 0.5, 1.0]
        # mass = sum(wt * var * frac)
        mass_old = sum(csu.cwtgcell_old .* var .* frac_old)
        CLM.update_column_state_no_special_handling!(csu, bounds, 1, col, var;
            fractional_area_old=frac_old, fractional_area_new=frac_new)
        mass_new = sum(csu.cwtgcell_new .* var .* frac_new)
        @test mass_new ≈ mass_old atol=1e-12

        # providing only one of the two fractional areas is an error
        var2 = [10.0, 4.0, 7.0]
        @test_throws ErrorException CLM.update_column_state_no_special_handling!(
            csu, bounds, 1, col, var2; fractional_area_old=frac_old)
    end

    # ---------------------------------------------------------------------
    # Variant 2: fill_special_using_natveg — seed special col from natveg
    # ---------------------------------------------------------------------
    @testset "variant 2: fill_special_using_natveg" begin
        bounds, grc, lun, col, pch = build_csu_subgrid()
        csu = CLM.ColumnStateUpdater(bounds, 1)

        col.active .= [true, true, true]
        col.wtgcell .= [0.3, 0.3, 0.4]
        CLM.set_old_weights!(csu, bounds, grc, lun, col)
        @test csu.natveg_template_col == [1, 1, 1]

        # Special lake (c3) shrinks 0.4->0.2; natveg c2 grows 0.3->0.5.
        col.wtgcell .= [0.3, 0.5, 0.2]
        CLM.set_new_weights!(csu, bounds, 1, col)

        # var: natveg cols carry 10 & 4; special lake var is irrelevant (not
        # prognostic) — its seed comes from natveg template col (c1 = 10).
        var = [10.0, 4.0, 999.0]
        non_conserved = zeros(Float64, bounds.endg)
        adjustment = zeros(Float64, bounds.endc)
        CLM.update_column_state_fill_special_using_natveg!(csu, bounds, 1,
            grc, lun, col, var, non_conserved; adjustment=adjustment)

        # Shrinking special col seeds value = var[template=c1] = 10.0 per area.
        #   loss = 0.2 * 10.0 = 2.0; area lost = 0.2; gain/area = 10.0
        #   mass_gained(c2) = 0.2 * 10.0 = 2.0
        #   var[2] = (0.3*4 + 2.0)/0.5 = 3.2/0.5 = 6.4
        @test var[2] ≈ 6.4
        @test var[1] == 10.0     # untouched natveg col
        # The special col's seed (10*0.2 = 2.0) was conjured out of thin air:
        # the grid cell GAINED that virtual mass, so non_conserved is negative
        # (positive = mass lost from cell, negative = mass gained by cell).
        @test non_conserved[1] ≈ -2.0
    end

    # ---------------------------------------------------------------------
    # Variant 3: fill_using_fixed_values
    # ---------------------------------------------------------------------
    @testset "variant 3: fill_using_fixed_values" begin
        bounds, grc, lun, col, pch = build_csu_subgrid()
        csu = CLM.ColumnStateUpdater(bounds, 1)

        col.active .= [true, true, true]
        col.wtgcell .= [0.3, 0.3, 0.4]
        CLM.set_old_weights!(csu, bounds, grc, lun, col)

        # Lake (c3) shrinks 0.4->0.2; natveg c1 grows 0.3->0.5.
        col.wtgcell .= [0.5, 0.3, 0.2]
        CLM.set_new_weights!(csu, bounds, 1, col)

        # landunit_values: soil -> use existing; lake -> fixed 5.0; rest existing.
        landunit_values = fill(CLM.FILLVAL_USE_EXISTING_VALUE, CLM.MAX_LUNIT)
        landunit_values[CLM.ISTDLAK] = 5.0

        var = [10.0, 4.0, 999.0]   # lake var irrelevant; fixed seed is 5.0
        non_conserved = zeros(Float64, bounds.endg)
        CLM.update_column_state_fill_using_fixed_values!(csu, bounds, 1,
            lun, col, var, landunit_values, non_conserved)

        # Shrinking lake seeds 5.0 per area:
        #   loss = 0.2 * 5.0 = 1.0; gain/area = 5.0; mass_gained(c1) = 0.2*5.0 = 1.0
        #   var[1] = (0.3*10 + 1.0)/0.5 = 4.0/0.5 = 8.0
        @test var[1] ≈ 8.0
        # 1.0 of fixed mass conjured into the system -> grid cell gained mass,
        # so non_conserved is negative.
        @test non_conserved[1] ≈ -1.0

        # length check
        @test_throws AssertionError CLM.update_column_state_fill_using_fixed_values!(
            csu, bounds, 1, lun, col, [10.0, 4.0, 999.0],
            [1.0, 2.0], zeros(Float64, bounds.endg))
    end

    # ---------------------------------------------------------------------
    # Variant 4: fill_special_using_fixed_value (wrapper over variant 3)
    # ---------------------------------------------------------------------
    @testset "variant 4: fill_special_using_fixed_value" begin
        bounds, grc, lun, col, pch = build_csu_subgrid()
        csu = CLM.ColumnStateUpdater(bounds, 1)

        col.active .= [true, true, true]
        col.wtgcell .= [0.3, 0.3, 0.4]
        CLM.set_old_weights!(csu, bounds, grc, lun, col)
        col.wtgcell .= [0.5, 0.3, 0.2]
        CLM.set_new_weights!(csu, bounds, 1, col)

        # All special landunits seed value 0.0; natveg keeps existing.
        var = [10.0, 4.0, 999.0]
        non_conserved = zeros(Float64, bounds.endg)
        CLM.update_column_state_fill_special_using_fixed_value!(csu, bounds, 1,
            lun, col, var, 0.0, non_conserved)

        # Lake seeds 0.0: loss = 0; gain/area = 0; var[1] unchanged in value...
        #   var[1] = (0.3*10 + 0)/0.5 = 6.0  (per-area dilution as col grew)
        @test var[1] ≈ 6.0
        @test non_conserved[1] ≈ 0.0   # special seed value 0 -> nothing conjured
    end

    # ---------------------------------------------------------------------
    # No changes -> updates are a no-op (and zero out adjustment)
    # ---------------------------------------------------------------------
    @testset "no changes -> no-op" begin
        bounds, grc, lun, col, pch = build_csu_subgrid()
        csu = CLM.ColumnStateUpdater(bounds, 1)

        col.active .= [true, true, true]
        col.wtgcell .= [0.3, 0.3, 0.4]
        CLM.set_old_weights!(csu, bounds, grc, lun, col)
        col.wtgcell .= [0.3, 0.3, 0.4]   # unchanged
        CLM.set_new_weights!(csu, bounds, 1, col)
        @test csu.any_changes[1] == false

        var = [10.0, 4.0, 7.0]
        adjustment = fill(99.0, bounds.endc)
        CLM.update_column_state_no_special_handling!(csu, bounds, 1, col, var;
            adjustment=adjustment)
        @test var == [10.0, 4.0, 7.0]      # untouched
        @test all(adjustment .== 0.0)      # zeroed even with no work
    end

end
