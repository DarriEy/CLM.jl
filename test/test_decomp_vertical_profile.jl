@testset "SoilBiogeochemVerticalProfile" begin

    # =====================================================================
    # Helper: build minimal data structures for testing
    # =====================================================================

    """
    Build a minimal test setup with `nc` columns and `np` patches.
    Each patch maps to a single column; all patches are vegetated (itype != noveg).
    Active layer index defaults to `nlevdecomp` (fully thawed).

    Returns a NamedTuple with all inputs for soil_bgc_vertical_profile!.
    """
    function make_test_setup(; nc::Int=2, np::Int=2,
                               nlevdecomp::Int=3,
                               nlevdecomp_full::Int=3,
                               alt_indx::Union{Vector{Int}, Nothing}=nothing,
                               patch_itypes::Union{Vector{Int}, Nothing}=nothing,
                               is_fates::Union{Vector{Bool}, Nothing}=nothing,
                               rootfr_vals::Union{Matrix{Float64}, Nothing}=nothing,
                               wtcol_vals::Union{Vector{Float64}, Nothing}=nothing,
                               patch_col_map::Union{Vector{Int}, Nothing}=nothing)

        # Vertical coordinate arrays (simple uniform layers)
        dz = fill(0.1, nlevdecomp)
        zs = [(j - 0.5) * 0.1 for j in 1:nlevdecomp]  # midpoints: 0.05, 0.15, ...

        # --- SoilBiogeochemStateData ---
        bgc_state = CLM.SoilBiogeochemStateData()
        bgc_state.leaf_prof_patch     = zeros(np, nlevdecomp_full)
        bgc_state.froot_prof_patch    = zeros(np, nlevdecomp_full)
        bgc_state.croot_prof_patch    = zeros(np, nlevdecomp_full)
        bgc_state.stem_prof_patch     = zeros(np, nlevdecomp_full)
        bgc_state.nfixation_prof_col  = zeros(nc, nlevdecomp_full)
        bgc_state.ndep_prof_col       = zeros(nc, nlevdecomp_full)

        # --- ActiveLayerData ---
        active_layer = CLM.ActiveLayerData()
        if alt_indx === nothing
            active_layer.altmax_lastyear_indx_col = fill(nlevdecomp, nc)
        else
            active_layer.altmax_lastyear_indx_col = alt_indx
        end

        # --- SoilStateData ---
        soilstate = CLM.SoilStateData()
        if rootfr_vals === nothing
            # Uniform root fraction across layers
            soilstate.crootfr_patch = fill(1.0 / nlevdecomp, np, nlevdecomp_full)
        else
            soilstate.crootfr_patch = rootfr_vals
        end

        # --- ColumnData ---
        col_data = CLM.ColumnData()
        col_data.nbedrock = fill(nlevdecomp, nc)
        if is_fates === nothing
            col_data.is_fates = fill(false, nc)
        else
            col_data.is_fates = is_fates
        end

        # --- PatchData ---
        patch_data = CLM.PatchData()
        if patch_itypes === nothing
            # Default: all vegetated (itype = 1, which is not noveg=0)
            patch_data.itype = fill(1, np)
        else
            patch_data.itype = patch_itypes
        end
        if patch_col_map === nothing
            # Map patch i -> column i (1:1 for simple cases)
            patch_data.column = collect(1:np)
        else
            patch_data.column = patch_col_map
        end
        if wtcol_vals === nothing
            patch_data.wtcol = fill(1.0, np)
        else
            patch_data.wtcol = wtcol_vals
        end

        # --- Masks ---
        mask_soilc = trues(nc)
        mask_vegp  = trues(np)

        return (bgc_state=bgc_state, active_layer=active_layer,
                soilstate=soilstate, col=col_data, patch=patch_data,
                mask_soilc=mask_soilc, mask_vegp=mask_vegp,
                nlevdecomp=nlevdecomp, nlevdecomp_full=nlevdecomp_full,
                dzsoi_decomp=dz, zsoi=zs)
    end

    # =====================================================================
    # Test: SURFPROF_EXP constant
    # =====================================================================
    @testset "SURFPROF_EXP constant" begin
        @test CLM.SURFPROF_EXP == 10.0
    end

    # =====================================================================
    # Test: Basic execution -- profiles integrate to 1
    # =====================================================================
    @testset "profiles integrate to 1 (basic case)" begin
        s = make_test_setup(nc=2, np=2, nlevdecomp=5, nlevdecomp_full=5)

        CLM.soil_bgc_vertical_profile!(
            s.bgc_state, s.active_layer, s.soilstate, s.col, s.patch;
            mask_bgc_soilc=s.mask_soilc,
            mask_bgc_vegp=s.mask_vegp,
            nlevdecomp=s.nlevdecomp,
            nlevdecomp_full=s.nlevdecomp_full,
            dzsoi_decomp=s.dzsoi_decomp,
            zsoi=s.zsoi)

        dz = s.dzsoi_decomp
        nl = s.nlevdecomp

        # Column profiles
        for c in 1:2
            ndep_sum = sum(s.bgc_state.ndep_prof_col[c, j] * dz[j] for j in 1:nl)
            nfix_sum = sum(s.bgc_state.nfixation_prof_col[c, j] * dz[j] for j in 1:nl)
            @test ndep_sum ≈ 1.0 atol=1e-9
            @test nfix_sum ≈ 1.0 atol=1e-9
        end

        # Patch profiles
        for p in 1:2
            froot_sum = sum(s.bgc_state.froot_prof_patch[p, j] * dz[j] for j in 1:nl)
            croot_sum = sum(s.bgc_state.croot_prof_patch[p, j] * dz[j] for j in 1:nl)
            leaf_sum  = sum(s.bgc_state.leaf_prof_patch[p, j]  * dz[j] for j in 1:nl)
            stem_sum  = sum(s.bgc_state.stem_prof_patch[p, j]  * dz[j] for j in 1:nl)
            @test froot_sum ≈ 1.0 atol=1e-9
            @test croot_sum ≈ 1.0 atol=1e-9
            @test leaf_sum  ≈ 1.0 atol=1e-9
            @test stem_sum  ≈ 1.0 atol=1e-9
        end
    end

    # =====================================================================
    # Test: Fully frozen column (altmax_lastyear_indx == 0)
    # =====================================================================
    @testset "fully frozen column -- top layer receives all input" begin
        s = make_test_setup(nc=1, np=1, nlevdecomp=4, nlevdecomp_full=4,
                            alt_indx=[0])

        CLM.soil_bgc_vertical_profile!(
            s.bgc_state, s.active_layer, s.soilstate, s.col, s.patch;
            mask_bgc_soilc=s.mask_soilc,
            mask_bgc_vegp=s.mask_vegp,
            nlevdecomp=s.nlevdecomp,
            nlevdecomp_full=s.nlevdecomp_full,
            dzsoi_decomp=s.dzsoi_decomp,
            zsoi=s.zsoi)

        dz1 = s.dzsoi_decomp[1]
        inv_dz1 = 1.0 / dz1

        # All profiles should be concentrated in layer 1
        @test s.bgc_state.froot_prof_patch[1, 1] ≈ inv_dz1
        @test s.bgc_state.croot_prof_patch[1, 1] ≈ inv_dz1
        @test s.bgc_state.leaf_prof_patch[1, 1]  ≈ inv_dz1
        @test s.bgc_state.stem_prof_patch[1, 1]  ≈ inv_dz1

        @test s.bgc_state.nfixation_prof_col[1, 1] ≈ inv_dz1
        @test s.bgc_state.ndep_prof_col[1, 1]      ≈ inv_dz1

        # Deeper layers should be zero
        for j in 2:4
            @test s.bgc_state.froot_prof_patch[1, j] == 0.0
            @test s.bgc_state.nfixation_prof_col[1, j] == 0.0
        end
    end

    # =====================================================================
    # Test: Noveg patch (itype == noveg) -- everything in top layer
    # =====================================================================
    @testset "noveg patch -- top layer" begin
        # patch type = noveg (0 in Fortran; CLM.noveg is the global)
        s = make_test_setup(nc=1, np=1, nlevdecomp=3, nlevdecomp_full=3,
                            patch_itypes=[CLM.noveg])

        CLM.soil_bgc_vertical_profile!(
            s.bgc_state, s.active_layer, s.soilstate, s.col, s.patch;
            mask_bgc_soilc=s.mask_soilc,
            mask_bgc_vegp=s.mask_vegp,
            nlevdecomp=s.nlevdecomp,
            nlevdecomp_full=s.nlevdecomp_full,
            dzsoi_decomp=s.dzsoi_decomp,
            zsoi=s.zsoi)

        dz1 = s.dzsoi_decomp[1]
        inv_dz1 = 1.0 / dz1

        # Patch profiles: noveg has cinput_rootfr=0, so rootfr_tot=0, fallback to top
        @test s.bgc_state.froot_prof_patch[1, 1] ≈ inv_dz1
        @test s.bgc_state.leaf_prof_patch[1, 1]  ≈ inv_dz1
    end

    # =====================================================================
    # Test: Surface profile exponential shape
    # =====================================================================
    @testset "surface profile is shallower than root profile" begin
        nlevdecomp = 5
        dz = fill(0.1, nlevdecomp)
        zs = [(j - 0.5) * 0.1 for j in 1:nlevdecomp]

        # Uniform root fraction
        rootfr = fill(1.0 / nlevdecomp, 1, nlevdecomp)

        s = make_test_setup(nc=1, np=1, nlevdecomp=nlevdecomp,
                            nlevdecomp_full=nlevdecomp,
                            rootfr_vals=rootfr)

        CLM.soil_bgc_vertical_profile!(
            s.bgc_state, s.active_layer, s.soilstate, s.col, s.patch;
            mask_bgc_soilc=s.mask_soilc,
            mask_bgc_vegp=s.mask_vegp,
            nlevdecomp=s.nlevdecomp,
            nlevdecomp_full=s.nlevdecomp_full,
            dzsoi_decomp=s.dzsoi_decomp,
            zsoi=s.zsoi)

        # With uniform roots, froot_prof should be constant across layers
        # leaf_prof (surface) should decrease with depth (exponential)
        leaf1 = s.bgc_state.leaf_prof_patch[1, 1]
        leaf_last = s.bgc_state.leaf_prof_patch[1, nlevdecomp]
        @test leaf1 > leaf_last  # surface profile decreases with depth

        # froot should be constant (uniform root distribution)
        froot1 = s.bgc_state.froot_prof_patch[1, 1]
        froot_last = s.bgc_state.froot_prof_patch[1, nlevdecomp]
        @test froot1 ≈ froot_last atol=1e-12
    end

    # =====================================================================
    # Test: Partial active layer (only first 2 of 5 levels thawed)
    # =====================================================================
    @testset "partial active layer" begin
        nlevdecomp = 5
        # Active layer index = 2 (only layers 1,2 are thawed)
        s = make_test_setup(nc=1, np=1, nlevdecomp=nlevdecomp,
                            nlevdecomp_full=nlevdecomp,
                            alt_indx=[2])

        CLM.soil_bgc_vertical_profile!(
            s.bgc_state, s.active_layer, s.soilstate, s.col, s.patch;
            mask_bgc_soilc=s.mask_soilc,
            mask_bgc_vegp=s.mask_vegp,
            nlevdecomp=s.nlevdecomp,
            nlevdecomp_full=s.nlevdecomp_full,
            dzsoi_decomp=s.dzsoi_decomp,
            zsoi=s.zsoi)

        # Layers 3..5 should have zero profile values
        for j in 3:5
            @test s.bgc_state.froot_prof_patch[1, j] == 0.0
            @test s.bgc_state.leaf_prof_patch[1, j]  == 0.0
            @test s.bgc_state.ndep_prof_col[1, j]    == 0.0
            @test s.bgc_state.nfixation_prof_col[1, j] == 0.0
        end

        # But the profiles over layers 1..2 should still integrate to 1
        dz = s.dzsoi_decomp
        froot_sum = sum(s.bgc_state.froot_prof_patch[1, j] * dz[j] for j in 1:nlevdecomp)
        ndep_sum  = sum(s.bgc_state.ndep_prof_col[1, j]    * dz[j] for j in 1:nlevdecomp)
        @test froot_sum ≈ 1.0 atol=1e-9
        @test ndep_sum  ≈ 1.0 atol=1e-9
    end

    # =====================================================================
    # Test: Multiple patches per column (aggregation via wtcol)
    # =====================================================================
    @testset "multiple patches aggregate to column" begin
        nc = 1; np = 2
        nlevdecomp = 3; nlevdecomp_full = 3

        # Patch 1: roots in layer 1 only
        # Patch 2: roots in layer 3 only
        rootfr = zeros(2, 3)
        rootfr[1, 1] = 1.0
        rootfr[2, 3] = 1.0

        # Both patches map to column 1, equal weights
        s = make_test_setup(nc=nc, np=np, nlevdecomp=nlevdecomp,
                            nlevdecomp_full=nlevdecomp_full,
                            rootfr_vals=rootfr,
                            wtcol_vals=[0.5, 0.5],
                            patch_col_map=[1, 1])

        CLM.soil_bgc_vertical_profile!(
            s.bgc_state, s.active_layer, s.soilstate, s.col, s.patch;
            mask_bgc_soilc=s.mask_soilc,
            mask_bgc_vegp=s.mask_vegp,
            nlevdecomp=s.nlevdecomp,
            nlevdecomp_full=s.nlevdecomp_full,
            dzsoi_decomp=s.dzsoi_decomp,
            zsoi=s.zsoi)

        dz = s.dzsoi_decomp

        # Column nfixation_prof should reflect aggregated rootfr
        # col_cinput_rootfr[1,1] = rootfr[1,1]/dz[1] * 0.5 = 0.5/dz
        # col_cinput_rootfr[1,3] = rootfr[2,3]/dz[3] * 0.5 = 0.5/dz
        # After normalisation, nfix layer 1 and 3 should be non-zero,
        # layer 2 should be zero
        @test s.bgc_state.nfixation_prof_col[1, 2] == 0.0

        nfix_sum = sum(s.bgc_state.nfixation_prof_col[1, j] * dz[j] for j in 1:nlevdecomp)
        @test nfix_sum ≈ 1.0 atol=1e-9

        # With equal dz and equal weights:
        # col_cinput_rootfr[1,1]/dz = col_cinput_rootfr[1,3]/dz
        # so nfixation_prof layers 1 and 3 should be equal
        @test s.bgc_state.nfixation_prof_col[1, 1] ≈ s.bgc_state.nfixation_prof_col[1, 3] atol=1e-12
    end

    # =====================================================================
    # Test: FATES column uses surface profile for both nfix and ndep
    # =====================================================================
    @testset "FATES column uses surface profile" begin
        nc = 1; np = 1
        nlevdecomp = 3; nlevdecomp_full = 3

        s = make_test_setup(nc=nc, np=np, nlevdecomp=nlevdecomp,
                            nlevdecomp_full=nlevdecomp_full,
                            is_fates=[true])

        CLM.soil_bgc_vertical_profile!(
            s.bgc_state, s.active_layer, s.soilstate, s.col, s.patch;
            mask_bgc_soilc=s.mask_soilc,
            mask_bgc_vegp=s.mask_vegp,
            nlevdecomp=s.nlevdecomp,
            nlevdecomp_full=s.nlevdecomp_full,
            dzsoi_decomp=s.dzsoi_decomp,
            zsoi=s.zsoi)

        dz = s.dzsoi_decomp

        # For FATES, nfixation_prof and ndep_prof should be identical
        for j in 1:nlevdecomp
            @test s.bgc_state.nfixation_prof_col[1, j] ≈ s.bgc_state.ndep_prof_col[1, j] atol=1e-15
        end

        # Both should integrate to 1
        nfix_sum = sum(s.bgc_state.nfixation_prof_col[1, j] * dz[j] for j in 1:nlevdecomp)
        ndep_sum = sum(s.bgc_state.ndep_prof_col[1, j]      * dz[j] for j in 1:nlevdecomp)
        @test nfix_sum ≈ 1.0 atol=1e-9
        @test ndep_sum ≈ 1.0 atol=1e-9

        # Surface profile should decrease with depth
        @test s.bgc_state.nfixation_prof_col[1, 1] > s.bgc_state.nfixation_prof_col[1, nlevdecomp]
    end

    # =====================================================================
    # Test: FATES frozen column -- top layer fallback
    # =====================================================================
    @testset "FATES frozen column -- top layer fallback" begin
        s = make_test_setup(nc=1, np=1, nlevdecomp=3, nlevdecomp_full=3,
                            alt_indx=[0],
                            is_fates=[true])

        CLM.soil_bgc_vertical_profile!(
            s.bgc_state, s.active_layer, s.soilstate, s.col, s.patch;
            mask_bgc_soilc=s.mask_soilc,
            mask_bgc_vegp=s.mask_vegp,
            nlevdecomp=s.nlevdecomp,
            nlevdecomp_full=s.nlevdecomp_full,
            dzsoi_decomp=s.dzsoi_decomp,
            zsoi=s.zsoi)

        inv_dz1 = 1.0 / s.dzsoi_decomp[1]
        @test s.bgc_state.nfixation_prof_col[1, 1] ≈ inv_dz1
        @test s.bgc_state.ndep_prof_col[1, 1]      ≈ inv_dz1
    end

    # =====================================================================
    # Test: Non-FATES: ndep uses surface_prof, nfix uses rootfr
    # =====================================================================
    @testset "non-FATES: ndep vs nfix profiles differ" begin
        nc = 1; np = 1
        nlevdecomp = 4; nlevdecomp_full = 4

        # Deep roots: all weight in bottom layer
        rootfr = zeros(1, 4)
        rootfr[1, 4] = 1.0

        s = make_test_setup(nc=nc, np=np, nlevdecomp=nlevdecomp,
                            nlevdecomp_full=nlevdecomp_full,
                            rootfr_vals=rootfr)

        CLM.soil_bgc_vertical_profile!(
            s.bgc_state, s.active_layer, s.soilstate, s.col, s.patch;
            mask_bgc_soilc=s.mask_soilc,
            mask_bgc_vegp=s.mask_vegp,
            nlevdecomp=s.nlevdecomp,
            nlevdecomp_full=s.nlevdecomp_full,
            dzsoi_decomp=s.dzsoi_decomp,
            zsoi=s.zsoi)

        dz = s.dzsoi_decomp

        # ndep should follow surface profile (shallow)
        # nfix should follow rootfr (deep -- all in layer 4)
        @test s.bgc_state.ndep_prof_col[1, 1] > s.bgc_state.ndep_prof_col[1, 4]

        # nfixation_prof should be zero in layers 1..3 (no roots there)
        for j in 1:3
            @test s.bgc_state.nfixation_prof_col[1, j] == 0.0
        end
        @test s.bgc_state.nfixation_prof_col[1, 4] > 0.0

        # Both still integrate to 1
        ndep_sum = sum(s.bgc_state.ndep_prof_col[1, j] * dz[j] for j in 1:nlevdecomp)
        nfix_sum = sum(s.bgc_state.nfixation_prof_col[1, j] * dz[j] for j in 1:nlevdecomp)
        @test ndep_sum ≈ 1.0 atol=1e-9
        @test nfix_sum ≈ 1.0 atol=1e-9
    end

    # =====================================================================
    # Test: use_bedrock zeroes surface profile below zmin_bedrock
    # =====================================================================
    @testset "use_bedrock zeroes deep surface profile" begin
        nlevdecomp = 5
        # Layer midpoints: 0.05, 0.15, 0.25, 0.35, 0.45
        # With zmin_bedrock=0.3, layers 4,5 are below bedrock
        dz = fill(0.1, nlevdecomp)
        zs = [(j - 0.5) * 0.1 for j in 1:nlevdecomp]

        s = make_test_setup(nc=1, np=1, nlevdecomp=nlevdecomp,
                            nlevdecomp_full=nlevdecomp)

        CLM.soil_bgc_vertical_profile!(
            s.bgc_state, s.active_layer, s.soilstate, s.col, s.patch;
            mask_bgc_soilc=s.mask_soilc,
            mask_bgc_vegp=s.mask_vegp,
            nlevdecomp=s.nlevdecomp,
            nlevdecomp_full=s.nlevdecomp_full,
            dzsoi_decomp=dz,
            zsoi=zs,
            use_bedrock=true,
            zmin_bedrock=0.3)

        # Leaf/stem profile (surface) layers 4,5 should be zero
        @test s.bgc_state.leaf_prof_patch[1, 4] == 0.0
        @test s.bgc_state.leaf_prof_patch[1, 5] == 0.0
        @test s.bgc_state.stem_prof_patch[1, 4] == 0.0
        @test s.bgc_state.stem_prof_patch[1, 5] == 0.0

        # ndep should also have zeros for those layers
        @test s.bgc_state.ndep_prof_col[1, 4] == 0.0
        @test s.bgc_state.ndep_prof_col[1, 5] == 0.0

        # Profiles still integrate to 1
        leaf_sum = sum(s.bgc_state.leaf_prof_patch[1, j] * dz[j] for j in 1:nlevdecomp)
        @test leaf_sum ≈ 1.0 atol=1e-9
    end

    # =====================================================================
    # Test: Masked-out columns/patches are skipped
    # =====================================================================
    @testset "masked-out columns and patches are skipped" begin
        nc = 2; np = 2; nlevdecomp = 3

        s = make_test_setup(nc=nc, np=np, nlevdecomp=nlevdecomp,
                            nlevdecomp_full=nlevdecomp)

        # Mask out column 2 and patch 2
        s.mask_soilc[2] = false
        s.mask_vegp[2]  = false

        CLM.soil_bgc_vertical_profile!(
            s.bgc_state, s.active_layer, s.soilstate, s.col, s.patch;
            mask_bgc_soilc=s.mask_soilc,
            mask_bgc_vegp=s.mask_vegp,
            nlevdecomp=s.nlevdecomp,
            nlevdecomp_full=s.nlevdecomp_full,
            dzsoi_decomp=s.dzsoi_decomp,
            zsoi=s.zsoi)

        # Column 2 ndep/nfix should remain zero (was initialized to 0 in step 2)
        for j in 1:nlevdecomp
            @test s.bgc_state.ndep_prof_col[2, j] == 0.0
            @test s.bgc_state.nfixation_prof_col[2, j] == 0.0
        end

        # Patch 2 profiles should remain zero
        for j in 1:nlevdecomp
            @test s.bgc_state.froot_prof_patch[2, j] == 0.0
            @test s.bgc_state.leaf_prof_patch[2, j]  == 0.0
        end

        # Column 1 and patch 1 should be correctly computed
        dz = s.dzsoi_decomp
        ndep_sum = sum(s.bgc_state.ndep_prof_col[1, j] * dz[j] for j in 1:nlevdecomp)
        @test ndep_sum ≈ 1.0 atol=1e-9
    end

    # =====================================================================
    # Test: croot_prof equals froot_prof
    # =====================================================================
    @testset "croot_prof equals froot_prof" begin
        s = make_test_setup(nc=1, np=1, nlevdecomp=4, nlevdecomp_full=4)

        CLM.soil_bgc_vertical_profile!(
            s.bgc_state, s.active_layer, s.soilstate, s.col, s.patch;
            mask_bgc_soilc=s.mask_soilc,
            mask_bgc_vegp=s.mask_vegp,
            nlevdecomp=s.nlevdecomp,
            nlevdecomp_full=s.nlevdecomp_full,
            dzsoi_decomp=s.dzsoi_decomp,
            zsoi=s.zsoi)

        for j in 1:4
            @test s.bgc_state.croot_prof_patch[1, j] ≈ s.bgc_state.froot_prof_patch[1, j] atol=1e-15
        end
    end

    # =====================================================================
    # Test: leaf_prof equals stem_prof
    # =====================================================================
    @testset "leaf_prof equals stem_prof" begin
        s = make_test_setup(nc=1, np=1, nlevdecomp=4, nlevdecomp_full=4)

        CLM.soil_bgc_vertical_profile!(
            s.bgc_state, s.active_layer, s.soilstate, s.col, s.patch;
            mask_bgc_soilc=s.mask_soilc,
            mask_bgc_vegp=s.mask_vegp,
            nlevdecomp=s.nlevdecomp,
            nlevdecomp_full=s.nlevdecomp_full,
            dzsoi_decomp=s.dzsoi_decomp,
            zsoi=s.zsoi)

        for j in 1:4
            @test s.bgc_state.leaf_prof_patch[1, j] ≈ s.bgc_state.stem_prof_patch[1, j] atol=1e-15
        end
    end

    # =====================================================================
    # Test: Single decomposition level (nlevdecomp=1, default CLM config)
    # =====================================================================
    @testset "single decomposition level (nlevdecomp=1)" begin
        s = make_test_setup(nc=1, np=1, nlevdecomp=1, nlevdecomp_full=1)

        CLM.soil_bgc_vertical_profile!(
            s.bgc_state, s.active_layer, s.soilstate, s.col, s.patch;
            mask_bgc_soilc=s.mask_soilc,
            mask_bgc_vegp=s.mask_vegp,
            nlevdecomp=1,
            nlevdecomp_full=1,
            dzsoi_decomp=s.dzsoi_decomp,
            zsoi=s.zsoi)

        dz = s.dzsoi_decomp[1]

        # With a single level, all profiles should be 1/dz
        @test s.bgc_state.froot_prof_patch[1, 1] * dz ≈ 1.0 atol=1e-12
        @test s.bgc_state.leaf_prof_patch[1, 1]  * dz ≈ 1.0 atol=1e-12
        @test s.bgc_state.ndep_prof_col[1, 1]    * dz ≈ 1.0 atol=1e-12
        @test s.bgc_state.nfixation_prof_col[1, 1] * dz ≈ 1.0 atol=1e-12
    end

    # =====================================================================
    # Test: Unequal patch weights sum correctly at column level
    # =====================================================================
    @testset "unequal patch weights" begin
        nc = 1; np = 3
        nlevdecomp = 3; nlevdecomp_full = 3

        rootfr = zeros(3, 3)
        rootfr[1, 1] = 1.0  # patch 1: roots in layer 1
        rootfr[2, 2] = 1.0  # patch 2: roots in layer 2
        rootfr[3, 3] = 1.0  # patch 3: roots in layer 3

        s = make_test_setup(nc=nc, np=np, nlevdecomp=nlevdecomp,
                            nlevdecomp_full=nlevdecomp_full,
                            rootfr_vals=rootfr,
                            wtcol_vals=[0.5, 0.3, 0.2],
                            patch_col_map=[1, 1, 1])

        CLM.soil_bgc_vertical_profile!(
            s.bgc_state, s.active_layer, s.soilstate, s.col, s.patch;
            mask_bgc_soilc=s.mask_soilc,
            mask_bgc_vegp=s.mask_vegp,
            nlevdecomp=s.nlevdecomp,
            nlevdecomp_full=s.nlevdecomp_full,
            dzsoi_decomp=s.dzsoi_decomp,
            zsoi=s.zsoi)

        dz = s.dzsoi_decomp

        # nfixation profile should integrate to 1
        nfix_sum = sum(s.bgc_state.nfixation_prof_col[1, j] * dz[j] for j in 1:nlevdecomp)
        @test nfix_sum ≈ 1.0 atol=1e-9

        # The relative weights of col_cinput_rootfr should be 0.5:0.3:0.2
        # (since all rootfr are 1.0 and dz is uniform)
        # nfixation_prof[c,j] = col_cinput_rootfr[c,j] / rootfr_tot
        # So nfixation_prof[1,1] * dz[1] / nfixation_prof[1,2] * dz[2] ~ 0.5/0.3
        ratio12 = (s.bgc_state.nfixation_prof_col[1, 1] * dz[1]) /
                  (s.bgc_state.nfixation_prof_col[1, 2] * dz[2])
        @test ratio12 ≈ 0.5 / 0.3 atol=1e-9
    end

    # =====================================================================
    # Test: Deep roots with shallow active layer -- fallback
    # =====================================================================
    @testset "roots below active layer trigger fallback" begin
        nc = 1; np = 1
        nlevdecomp = 5; nlevdecomp_full = 5

        # All roots in layer 5, but active layer only reaches layer 2
        rootfr = zeros(1, 5)
        rootfr[1, 5] = 1.0

        s = make_test_setup(nc=nc, np=np, nlevdecomp=nlevdecomp,
                            nlevdecomp_full=nlevdecomp_full,
                            rootfr_vals=rootfr,
                            alt_indx=[2])

        CLM.soil_bgc_vertical_profile!(
            s.bgc_state, s.active_layer, s.soilstate, s.col, s.patch;
            mask_bgc_soilc=s.mask_soilc,
            mask_bgc_vegp=s.mask_vegp,
            nlevdecomp=s.nlevdecomp,
            nlevdecomp_full=s.nlevdecomp_full,
            dzsoi_decomp=s.dzsoi_decomp,
            zsoi=s.zsoi)

        inv_dz1 = 1.0 / s.dzsoi_decomp[1]

        # rootfr_tot will be 0 (roots layer 5 is below alt_indx=2)
        # => fallback: everything in top layer
        @test s.bgc_state.froot_prof_patch[1, 1] ≈ inv_dz1
        @test s.bgc_state.leaf_prof_patch[1, 1]  ≈ inv_dz1
    end

end
