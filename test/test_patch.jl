@testset "PatchData" begin

    @testset "default construction" begin
        pch = CLM.PatchData()
        @test length(pch.column) == 0
        @test length(pch.wtcol) == 0
        @test length(pch.landunit) == 0
        @test length(pch.wtlunit) == 0
        @test length(pch.gridcell) == 0
        @test length(pch.wtgcell) == 0
        @test length(pch.itype) == 0
        @test length(pch.mxy) == 0
        @test length(pch.active) == 0
        @test length(pch.is_veg) == 0
        @test length(pch.is_bareground) == 0
        @test length(pch.wt_ed) == 0
        @test length(pch.sp_pftorder_index) == 0
        @test length(pch.is_fates) == 0
    end

    @testset "patch_init! (use_fates=false)" begin
        # Ensure FATES is off
        CLM.varctl.use_fates = false

        np = 10
        pch = CLM.PatchData()
        CLM.patch_init!(pch, np)

        # --- Check sizes: g/l/c/p hierarchy ---
        @test length(pch.column) == np
        @test length(pch.wtcol) == np
        @test length(pch.landunit) == np
        @test length(pch.wtlunit) == np
        @test length(pch.gridcell) == np
        @test length(pch.wtgcell) == np

        # --- Check sizes: non-ED fields ---
        @test length(pch.itype) == np
        @test length(pch.mxy) == np
        @test length(pch.active) == np
        @test length(pch.is_fates) == np

        # --- FATES fields should remain empty when use_fates=false ---
        @test length(pch.is_veg) == 0
        @test length(pch.is_bareground) == 0
        @test length(pch.wt_ed) == 0
        @test length(pch.sp_pftorder_index) == 0

        # --- Check integer sentinel values ---
        @test all(==(CLM.ISPVAL), pch.column)
        @test all(==(CLM.ISPVAL), pch.landunit)
        @test all(==(CLM.ISPVAL), pch.gridcell)
        @test all(==(CLM.ISPVAL), pch.itype)
        @test all(==(CLM.ISPVAL), pch.mxy)

        # --- Check NaN initialization for real fields ---
        @test all(isnan, pch.wtcol)
        @test all(isnan, pch.wtlunit)
        @test all(isnan, pch.wtgcell)

        # --- Check boolean defaults ---
        @test all(.!pch.active)
        @test all(.!pch.is_fates)
    end

    @testset "patch_init! (use_fates=true)" begin
        # Enable FATES
        CLM.varctl.use_fates = true

        np = 8
        pch = CLM.PatchData()
        CLM.patch_init!(pch, np)

        # --- FATES fields should be allocated ---
        @test length(pch.is_veg) == np
        @test length(pch.is_bareground) == np
        @test length(pch.wt_ed) == np
        @test length(pch.sp_pftorder_index) == np

        # --- Check FATES defaults ---
        @test all(.!pch.is_veg)
        @test all(.!pch.is_bareground)
        @test all(isnan, pch.wt_ed)
        @test all(isnan, pch.sp_pftorder_index)

        # --- Core fields still valid ---
        @test length(pch.column) == np
        @test all(==(CLM.ISPVAL), pch.column)
        @test all(isnan, pch.wtcol)
        @test all(.!pch.active)

        # Reset FATES flag
        CLM.varctl.use_fates = false
    end

    @testset "patch_clean! (use_fates=false)" begin
        CLM.varctl.use_fates = false
        pch = CLM.PatchData()
        CLM.patch_init!(pch, 10)
        CLM.patch_clean!(pch)

        # All vectors should be empty
        @test length(pch.column) == 0
        @test length(pch.wtcol) == 0
        @test length(pch.landunit) == 0
        @test length(pch.wtlunit) == 0
        @test length(pch.gridcell) == 0
        @test length(pch.wtgcell) == 0
        @test length(pch.itype) == 0
        @test length(pch.mxy) == 0
        @test length(pch.active) == 0
        @test length(pch.is_fates) == 0
    end

    @testset "patch_clean! (use_fates=true)" begin
        CLM.varctl.use_fates = true
        pch = CLM.PatchData()
        CLM.patch_init!(pch, 10)
        CLM.patch_clean!(pch)

        # All vectors should be empty, including FATES fields
        @test length(pch.column) == 0
        @test length(pch.is_veg) == 0
        @test length(pch.is_bareground) == 0
        @test length(pch.wt_ed) == 0
        @test length(pch.sp_pftorder_index) == 0

        # Reset FATES flag
        CLM.varctl.use_fates = false
    end

    @testset "field mutability" begin
        CLM.varctl.use_fates = true
        pch = CLM.PatchData()
        CLM.patch_init!(pch, 5)

        # Write to hierarchy fields and verify
        pch.column[1] = 3
        pch.wtcol[1] = 0.5
        pch.landunit[2] = 7
        pch.wtlunit[2] = 0.25
        pch.gridcell[3] = 1
        pch.wtgcell[3] = 1.0

        # Write to non-ED fields
        pch.itype[1] = 4
        pch.mxy[1] = 2
        pch.active[1] = true

        # Write to FATES fields
        pch.is_veg[1] = true
        pch.is_bareground[2] = true
        pch.wt_ed[1] = 0.8
        pch.sp_pftorder_index[1] = 3.0
        pch.is_fates[1] = true

        @test pch.column[1] == 3
        @test pch.wtcol[1] == 0.5
        @test pch.landunit[2] == 7
        @test pch.wtlunit[2] == 0.25
        @test pch.gridcell[3] == 1
        @test pch.wtgcell[3] == 1.0
        @test pch.itype[1] == 4
        @test pch.mxy[1] == 2
        @test pch.active[1] == true
        @test pch.is_veg[1] == true
        @test pch.is_bareground[2] == true
        @test pch.wt_ed[1] == 0.8
        @test pch.sp_pftorder_index[1] == 3.0
        @test pch.is_fates[1] == true

        # Reset FATES flag
        CLM.varctl.use_fates = false
    end

    @testset "re-init overwrites previous state" begin
        CLM.varctl.use_fates = false
        pch = CLM.PatchData()
        CLM.patch_init!(pch, 5)
        pch.wtgcell[1] = 999.0

        # Re-init with different size
        CLM.patch_init!(pch, 12)
        @test length(pch.wtgcell) == 12
        @test all(isnan, pch.wtgcell)  # old value gone
    end

end
