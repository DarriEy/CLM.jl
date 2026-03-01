@testset "SurfaceAlbedoData" begin

    # Ensure varpar is initialized for dimension parameters
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)

    numrad  = CLM.NUMRAD
    nlevcan = CLM.NLEVCAN
    nlevsno = CLM.varpar.nlevsno

    @testset "default construction" begin
        sa = CLM.SurfaceAlbedoData()
        @test length(sa.azsun_grc) == 0
        @test length(sa.coszen_grc) == 0
        @test length(sa.coszen_col) == 0
        @test size(sa.albd_patch) == (0, 0)
        @test size(sa.albgrd_col) == (0, 0)
        @test size(sa.flx_absdv_col) == (0, 0)
        @test length(sa.ncan_patch) == 0
        @test length(sa.vcmaxcintsun_patch) == 0
    end

    @testset "surfalb_init!" begin
        np = 12
        nc = 8
        ng = 3
        sa = CLM.SurfaceAlbedoData()
        CLM.surfalb_init!(sa, np, nc, ng)

        # --- Check gridcell-level sizes ---
        @test length(sa.azsun_grc) == ng
        @test length(sa.coszen_grc) == ng

        # --- Check column-level 1D sizes ---
        @test length(sa.coszen_col) == nc

        # --- Check column-level 2D sizes ---
        @test size(sa.albgrd_col) == (nc, numrad)
        @test size(sa.albgri_col) == (nc, numrad)
        @test size(sa.albsod_col) == (nc, numrad)
        @test size(sa.albsoi_col) == (nc, numrad)
        @test size(sa.albsnd_hst_col) == (nc, numrad)
        @test size(sa.albsni_hst_col) == (nc, numrad)
        @test size(sa.albgrd_pur_col) == (nc, numrad)
        @test size(sa.albgri_pur_col) == (nc, numrad)
        @test size(sa.albgrd_bc_col) == (nc, numrad)
        @test size(sa.albgri_bc_col) == (nc, numrad)
        @test size(sa.albgrd_oc_col) == (nc, numrad)
        @test size(sa.albgri_oc_col) == (nc, numrad)
        @test size(sa.albgrd_dst_col) == (nc, numrad)
        @test size(sa.albgri_dst_col) == (nc, numrad)

        # --- Check patch-level 2D sizes ---
        @test size(sa.albd_patch) == (np, numrad)
        @test size(sa.albi_patch) == (np, numrad)
        @test size(sa.albdSF_patch) == (np, numrad)
        @test size(sa.albiSF_patch) == (np, numrad)
        @test size(sa.ftdd_patch) == (np, numrad)
        @test size(sa.fabd_patch) == (np, numrad)
        @test size(sa.fabi_patch) == (np, numrad)
        @test size(sa.fabd_sun_patch) == (np, numrad)
        @test size(sa.fabd_sha_patch) == (np, numrad)
        @test size(sa.fabi_sun_patch) == (np, numrad)
        @test size(sa.fabi_sha_patch) == (np, numrad)

        # --- Check canopy layer sizes ---
        @test size(sa.fabd_sun_z_patch) == (np, nlevcan)
        @test size(sa.fabd_sha_z_patch) == (np, nlevcan)
        @test size(sa.fabi_sun_z_patch) == (np, nlevcan)
        @test size(sa.fabi_sha_z_patch) == (np, nlevcan)
        @test size(sa.fsun_z_patch) == (np, nlevcan)
        @test size(sa.tlai_z_patch) == (np, nlevcan)
        @test size(sa.tsai_z_patch) == (np, nlevcan)

        # --- Check snow layer flux sizes ---
        nlev_sno1 = nlevsno + 1
        @test size(sa.flx_absdv_col) == (nc, nlev_sno1)
        @test size(sa.flx_absdn_col) == (nc, nlev_sno1)
        @test size(sa.flx_absiv_col) == (nc, nlev_sno1)
        @test size(sa.flx_absin_col) == (nc, nlev_sno1)

        # --- Check 1D patch sizes ---
        @test length(sa.ncan_patch) == np
        @test length(sa.nrad_patch) == np
        @test length(sa.vcmaxcintsun_patch) == np
        @test length(sa.vcmaxcintsha_patch) == np

        # --- Check SNICAR history sizes ---
        @test size(sa.albgrd_hst_col) == (nc, numrad)
        @test size(sa.albd_hst_patch) == (np, numrad)

        # --- Check NaN initialization ---
        @test all(isnan, sa.azsun_grc)
        @test all(isnan, sa.coszen_grc)
        @test all(isnan, sa.coszen_col)
        @test all(isnan, sa.albgrd_col)
        @test all(isnan, sa.albd_patch)
        @test all(isnan, sa.vcmaxcintsun_patch)

        # --- Check SPVAL initialization ---
        @test all(x -> x == CLM.SPVAL, sa.albsnd_hst_col)
        @test all(x -> x == CLM.SPVAL, sa.albsoi_col)
        @test all(x -> x == CLM.SPVAL, sa.flx_absdv_col)
        @test all(x -> x == CLM.SPVAL, sa.albgrd_hst_col)
        @test all(x -> x == CLM.SPVAL, sa.albd_hst_patch)

        # --- Check zero initialization ---
        @test all(x -> x == 0.0, sa.fabd_sun_z_patch)
        @test all(x -> x == 0.0, sa.fsun_z_patch)
        @test all(x -> x == 0.0, sa.tlai_z_patch)
        @test all(x -> x == 0, sa.ncan_patch)
        @test all(x -> x == 0, sa.nrad_patch)
    end

    @testset "surfalb_clean!" begin
        sa = CLM.SurfaceAlbedoData()
        CLM.surfalb_init!(sa, 10, 6, 3)
        CLM.surfalb_clean!(sa)

        @test length(sa.azsun_grc) == 0
        @test length(sa.coszen_grc) == 0
        @test length(sa.coszen_col) == 0
        @test size(sa.albgrd_col) == (0, 0)
        @test size(sa.albd_patch) == (0, 0)
        @test size(sa.ftdd_patch) == (0, 0)
        @test size(sa.fabd_sun_z_patch) == (0, 0)
        @test size(sa.flx_absdv_col) == (0, 0)
        @test size(sa.fsun_z_patch) == (0, 0)
        @test length(sa.ncan_patch) == 0
        @test length(sa.vcmaxcintsun_patch) == 0
        @test size(sa.albgrd_hst_col) == (0, 0)
        @test size(sa.albd_hst_patch) == (0, 0)
    end

    @testset "surfalb_init_cold! basic" begin
        np = 8
        nc = 5
        ng = 2
        sa = CLM.SurfaceAlbedoData()
        CLM.surfalb_init!(sa, np, nc, ng)
        CLM.surfalb_init_cold!(sa, 1:nc, 1:np)

        # Column-level albedo defaults
        for c in 1:nc
            for ib in 1:numrad
                @test sa.albgrd_col[c, ib] ≈ 0.2
                @test sa.albgri_col[c, ib] ≈ 0.2
                @test sa.albsod_col[c, ib] ≈ 0.2
                @test sa.albsoi_col[c, ib] ≈ 0.2
                @test sa.albsnd_hst_col[c, ib] ≈ 0.6
                @test sa.albsni_hst_col[c, ib] ≈ 0.6
                @test sa.albgrd_pur_col[c, ib] ≈ 0.2
                @test sa.albgri_pur_col[c, ib] ≈ 0.2
                @test sa.albgrd_bc_col[c, ib] ≈ 0.2
                @test sa.albgri_bc_col[c, ib] ≈ 0.2
                @test sa.albgrd_oc_col[c, ib] ≈ 0.2
                @test sa.albgri_oc_col[c, ib] ≈ 0.2
                @test sa.albgrd_dst_col[c, ib] ≈ 0.2
                @test sa.albgri_dst_col[c, ib] ≈ 0.2
            end
        end

        # Patch-level defaults
        for p in 1:np
            for ib in 1:numrad
                @test sa.albd_patch[p, ib] ≈ 0.2
                @test sa.albi_patch[p, ib] ≈ 0.2
                @test sa.fabd_patch[p, ib] ≈ 0.0
                @test sa.fabi_patch[p, ib] ≈ 0.0
                @test sa.ftdd_patch[p, ib] ≈ 1.0
                @test sa.ftid_patch[p, ib] ≈ 0.0
                @test sa.ftii_patch[p, ib] ≈ 1.0
            end
        end

        # Snow-free albedos should NOT be set (use_SSRE=false by default)
        # They remain NaN from init!
        @test all(isnan, sa.albdSF_patch)
        @test all(isnan, sa.albiSF_patch)
    end

    @testset "surfalb_init_cold! with use_SSRE" begin
        np = 4
        nc = 3
        ng = 1
        sa = CLM.SurfaceAlbedoData()
        CLM.surfalb_init!(sa, np, nc, ng)
        CLM.surfalb_init_cold!(sa, 1:nc, 1:np; use_SSRE=true)

        for p in 1:np
            for ib in 1:numrad
                @test sa.albdSF_patch[p, ib] ≈ 0.2
                @test sa.albiSF_patch[p, ib] ≈ 0.2
            end
        end
    end

    @testset "field mutability" begin
        sa = CLM.SurfaceAlbedoData()
        CLM.surfalb_init!(sa, 5, 4, 2)

        # Write and verify gridcell field
        sa.coszen_grc[1] = 0.5
        @test sa.coszen_grc[1] == 0.5

        # Write and verify column 2D field
        sa.albgrd_col[2, 1] = 0.35
        @test sa.albgrd_col[2, 1] == 0.35

        # Write and verify patch 2D field
        sa.albd_patch[3, 2] = 0.42
        @test sa.albd_patch[3, 2] == 0.42

        # Write and verify integer field
        sa.ncan_patch[1] = 5
        @test sa.ncan_patch[1] == 5

        # Write and verify snow flux field
        sa.flx_absdv_col[1, 1] = 0.99
        @test sa.flx_absdv_col[1, 1] == 0.99
    end

    @testset "re-init overwrites previous state" begin
        sa = CLM.SurfaceAlbedoData()
        CLM.surfalb_init!(sa, 5, 4, 2)
        sa.coszen_grc[1] = 999.0

        CLM.surfalb_init!(sa, 10, 8, 3)
        @test length(sa.coszen_grc) == 3
        @test all(isnan, sa.coszen_grc)
        @test size(sa.albgrd_col) == (8, numrad)
        @test size(sa.albd_patch) == (10, numrad)
    end

    @testset "stub functions run without error" begin
        sa = CLM.SurfaceAlbedoData()
        CLM.surfalb_init!(sa, 6, 4, 2)

        @test CLM.surfalb_init_history!(sa, 1:4, 1:6, 1:2) === nothing
        @test CLM.surfalb_restart!(sa, 1:4, 1:6) === nothing
    end

    @testset "surfalb_init_history! sets SPVAL" begin
        sa = CLM.SurfaceAlbedoData()
        CLM.surfalb_init!(sa, 6, 4, 2)

        CLM.surfalb_init_history!(sa, 1:4, 1:6, 1:2)

        for g in 1:2
            @test sa.azsun_grc[g] == CLM.SPVAL
            @test sa.coszen_grc[g] == CLM.SPVAL
        end
        for c in 1:4
            @test sa.coszen_col[c] == CLM.SPVAL
            for ib in 1:numrad
                @test sa.albgrd_col[c, ib] == CLM.SPVAL
                @test sa.albgri_col[c, ib] == CLM.SPVAL
            end
        end
        for p in 1:6
            for ib in 1:numrad
                @test sa.albd_patch[p, ib] == CLM.SPVAL
                @test sa.albi_patch[p, ib] == CLM.SPVAL
            end
        end
    end

end
