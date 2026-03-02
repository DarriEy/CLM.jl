@testset "N Leaching" begin

    # ----------------------------------------------------------------
    # Helper: create minimal data structures for nc columns
    # ----------------------------------------------------------------
    function make_n_leaching_data(; nc=4, nlevdecomp=1, nlevsoi=3)
        ndecomp_pools = 7
        ndecomp_cascade_transitions = 10

        # --- N Leaching Parameters ---
        params = CLM.NLeachingParams(sf=0.1, sf_no3=1.0)

        # --- Nitrogen flux data ---
        nf = CLM.SoilBiogeochemNitrogenFluxData()
        CLM.soil_bgc_nitrogen_flux_init!(nf, nc, nlevsoi, ndecomp_pools,
                                          ndecomp_cascade_transitions)

        # --- Nitrogen state data ---
        ns = CLM.SoilBiogeochemNitrogenStateData()
        CLM.soil_bgc_nitrogen_state_init!(ns, nc, 1, nlevsoi, ndecomp_pools)

        # Set mineral N pools
        for j in 1:nlevsoi
            for c in 1:nc
                ns.sminn_vr_col[c, j]    = 5.0    # gN/m3
                ns.smin_no3_vr_col[c, j] = 2.0    # gN/m3
            end
        end

        # --- Water state arrays ---
        h2osoi_liq = fill(30.0, nc, nlevsoi)  # kg/m2

        # --- Water flux arrays ---
        qflx_drain = fill(1.0e-5, nc)   # mm H2O/s (subsurface drainage)
        qflx_surf  = fill(5.0e-6, nc)   # mm H2O/s (surface runoff)

        # --- Column geometry ---
        col_dz = fill(0.1, nc, nlevsoi)

        # --- Vertical coordinate (interface depths) ---
        # zisoi has length nlevsoi+1, with zisoi[1]=0 (surface) and
        # zisoi[j+1] = interface depth below layer j.
        # Matches the convention in varcon.jl.
        zisoi = zeros(nlevsoi + 1)
        for j in 1:nlevsoi
            zisoi[j+1] = zisoi[j] + col_dz[1, j]
        end

        dt = 1800.0  # seconds (30 minute timestep)

        mask_bgc_soilc = trues(nc)

        return (params=params, nf=nf, ns=ns,
                h2osoi_liq=h2osoi_liq, qflx_drain=qflx_drain, qflx_surf=qflx_surf,
                col_dz=col_dz, zisoi=zisoi, dt=dt,
                mask_bgc_soilc=mask_bgc_soilc,
                nlevdecomp=nlevdecomp, nlevsoi=nlevsoi, nc=nc)
    end

    # ================================================================
    # Test 1: NLeachingParams construction and read
    # ================================================================
    @testset "NLeachingParams construction" begin
        params = CLM.NLeachingParams()
        @test params.sf == 0.1
        @test params.sf_no3 == 1.0
    end

    @testset "n_leaching_read_params!" begin
        params = CLM.NLeachingParams()
        CLM.n_leaching_read_params!(params; sf=0.05, sf_no3=0.8)
        @test params.sf == 0.05
        @test params.sf_no3 == 0.8
    end

    # ================================================================
    # Test 2: n_leaching! with use_nitrif_denitrif=false
    # ================================================================
    @testset "n_leaching! (nitrif_denitrif OFF)" begin
        d = make_n_leaching_data()

        CLM.n_leaching!(d.nf, d.ns, d.params;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:d.nc,
            nlevdecomp=d.nlevdecomp,
            nlevsoi=d.nlevsoi,
            dt=d.dt,
            h2osoi_liq=d.h2osoi_liq,
            qflx_drain=d.qflx_drain,
            qflx_surf=d.qflx_surf,
            col_dz=d.col_dz,
            zisoi=d.zisoi,
            use_nitrif_denitrif=false)

        for c in 1:d.nc
            # sminn_leached_vr should be non-negative
            @test d.nf.sminn_leached_vr_col[c, 1] >= 0.0
            # With positive sminn, drainage, and water, leaching should be positive
            @test d.nf.sminn_leached_vr_col[c, 1] > 0.0
        end
    end

    # ================================================================
    # Test 3: n_leaching! with use_nitrif_denitrif=true
    # ================================================================
    @testset "n_leaching! (nitrif_denitrif ON)" begin
        d = make_n_leaching_data()

        CLM.n_leaching!(d.nf, d.ns, d.params;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:d.nc,
            nlevdecomp=d.nlevdecomp,
            nlevsoi=d.nlevsoi,
            dt=d.dt,
            h2osoi_liq=d.h2osoi_liq,
            qflx_drain=d.qflx_drain,
            qflx_surf=d.qflx_surf,
            col_dz=d.col_dz,
            zisoi=d.zisoi,
            use_nitrif_denitrif=true)

        for c in 1:d.nc
            # smin_no3_leached_vr should be non-negative
            @test d.nf.smin_no3_leached_vr_col[c, 1] >= 0.0
            # With positive NO3, drainage, and water, leaching should be positive
            @test d.nf.smin_no3_leached_vr_col[c, 1] > 0.0
            # smin_no3_runoff_vr should be non-negative
            @test d.nf.smin_no3_runoff_vr_col[c, 1] >= 0.0
        end
    end

    # ================================================================
    # Test 4: Zero drainage produces zero leaching
    # ================================================================
    @testset "Zero drainage produces zero leaching" begin
        d = make_n_leaching_data()
        d.qflx_drain .= 0.0

        CLM.n_leaching!(d.nf, d.ns, d.params;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:d.nc,
            nlevdecomp=d.nlevdecomp,
            nlevsoi=d.nlevsoi,
            dt=d.dt,
            h2osoi_liq=d.h2osoi_liq,
            qflx_drain=d.qflx_drain,
            qflx_surf=d.qflx_surf,
            col_dz=d.col_dz,
            zisoi=d.zisoi,
            use_nitrif_denitrif=false)

        for c in 1:d.nc
            @test d.nf.sminn_leached_vr_col[c, 1] == 0.0
        end
    end

    # ================================================================
    # Test 5: Zero N produces zero leaching
    # ================================================================
    @testset "Zero N produces zero leaching" begin
        d = make_n_leaching_data()
        d.ns.sminn_vr_col .= 0.0
        d.ns.smin_no3_vr_col .= 0.0

        # Test both modes
        CLM.n_leaching!(d.nf, d.ns, d.params;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:d.nc,
            nlevdecomp=d.nlevdecomp,
            nlevsoi=d.nlevsoi,
            dt=d.dt,
            h2osoi_liq=d.h2osoi_liq,
            qflx_drain=d.qflx_drain,
            qflx_surf=d.qflx_surf,
            col_dz=d.col_dz,
            zisoi=d.zisoi,
            use_nitrif_denitrif=false)

        for c in 1:d.nc
            @test d.nf.sminn_leached_vr_col[c, 1] == 0.0
        end

        CLM.n_leaching!(d.nf, d.ns, d.params;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:d.nc,
            nlevdecomp=d.nlevdecomp,
            nlevsoi=d.nlevsoi,
            dt=d.dt,
            h2osoi_liq=d.h2osoi_liq,
            qflx_drain=d.qflx_drain,
            qflx_surf=d.qflx_surf,
            col_dz=d.col_dz,
            zisoi=d.zisoi,
            use_nitrif_denitrif=true)

        for c in 1:d.nc
            @test d.nf.smin_no3_leached_vr_col[c, 1] == 0.0
            @test d.nf.smin_no3_runoff_vr_col[c, 1] == 0.0
        end
    end

    # ================================================================
    # Test 6: Higher drainage produces higher leaching
    # ================================================================
    @testset "Higher drainage produces higher leaching" begin
        d_high = make_n_leaching_data()
        d_low  = make_n_leaching_data()

        d_high.qflx_drain .= 1.0e-4
        d_low.qflx_drain  .= 1.0e-6

        CLM.n_leaching!(d_high.nf, d_high.ns, d_high.params;
            mask_bgc_soilc=d_high.mask_bgc_soilc,
            bounds=1:d_high.nc,
            nlevdecomp=d_high.nlevdecomp,
            nlevsoi=d_high.nlevsoi,
            dt=d_high.dt,
            h2osoi_liq=d_high.h2osoi_liq,
            qflx_drain=d_high.qflx_drain,
            qflx_surf=d_high.qflx_surf,
            col_dz=d_high.col_dz,
            zisoi=d_high.zisoi,
            use_nitrif_denitrif=false)

        CLM.n_leaching!(d_low.nf, d_low.ns, d_low.params;
            mask_bgc_soilc=d_low.mask_bgc_soilc,
            bounds=1:d_low.nc,
            nlevdecomp=d_low.nlevdecomp,
            nlevsoi=d_low.nlevsoi,
            dt=d_low.dt,
            h2osoi_liq=d_low.h2osoi_liq,
            qflx_drain=d_low.qflx_drain,
            qflx_surf=d_low.qflx_surf,
            col_dz=d_low.col_dz,
            zisoi=d_low.zisoi,
            use_nitrif_denitrif=false)

        @test d_high.nf.sminn_leached_vr_col[1, 1] > d_low.nf.sminn_leached_vr_col[1, 1]
    end

    # ================================================================
    # Test 7: Leaching capped at soluble fraction / dt
    # ================================================================
    @testset "Leaching rate capped" begin
        d = make_n_leaching_data()

        # Set extremely high drainage to force the cap
        d.qflx_drain .= 1.0e3

        CLM.n_leaching!(d.nf, d.ns, d.params;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:d.nc,
            nlevdecomp=d.nlevdecomp,
            nlevsoi=d.nlevsoi,
            dt=d.dt,
            h2osoi_liq=d.h2osoi_liq,
            qflx_drain=d.qflx_drain,
            qflx_surf=d.qflx_surf,
            col_dz=d.col_dz,
            zisoi=d.zisoi,
            use_nitrif_denitrif=false)

        for c in 1:d.nc
            # Should be capped at sf * sminn_vr / dt
            expected_cap = d.params.sf * d.ns.sminn_vr_col[c, 1] / d.dt
            @test d.nf.sminn_leached_vr_col[c, 1] <= expected_cap + 1.0e-15
        end
    end

    # ================================================================
    # Test 8: Mask filtering
    # ================================================================
    @testset "Mask filtering" begin
        d = make_n_leaching_data()

        # Mask out columns 2 and 4
        d.mask_bgc_soilc[2] = false
        d.mask_bgc_soilc[4] = false

        CLM.n_leaching!(d.nf, d.ns, d.params;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:d.nc,
            nlevdecomp=d.nlevdecomp,
            nlevsoi=d.nlevsoi,
            dt=d.dt,
            h2osoi_liq=d.h2osoi_liq,
            qflx_drain=d.qflx_drain,
            qflx_surf=d.qflx_surf,
            col_dz=d.col_dz,
            zisoi=d.zisoi,
            use_nitrif_denitrif=false)

        # Active columns should have computed values
        @test d.nf.sminn_leached_vr_col[1, 1] > 0.0
        @test d.nf.sminn_leached_vr_col[3, 1] > 0.0

        # Masked columns should still have NaN (not updated)
        @test isnan(d.nf.sminn_leached_vr_col[2, 1])
        @test isnan(d.nf.sminn_leached_vr_col[4, 1])
    end

    # ================================================================
    # Test 9: Zero surface runoff produces zero NO3 runoff
    # ================================================================
    @testset "Zero surface runoff produces zero NO3 runoff" begin
        d = make_n_leaching_data()
        d.qflx_surf .= 0.0

        CLM.n_leaching!(d.nf, d.ns, d.params;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:d.nc,
            nlevdecomp=d.nlevdecomp,
            nlevsoi=d.nlevsoi,
            dt=d.dt,
            h2osoi_liq=d.h2osoi_liq,
            qflx_drain=d.qflx_drain,
            qflx_surf=d.qflx_surf,
            col_dz=d.col_dz,
            zisoi=d.zisoi,
            use_nitrif_denitrif=true)

        for c in 1:d.nc
            @test d.nf.smin_no3_runoff_vr_col[c, 1] == 0.0
        end
    end

    # ================================================================
    # Test 10: Multi-level computation
    # ================================================================
    @testset "Multi-level computation" begin
        nlevdecomp = 3
        nlevsoi = 3
        d = make_n_leaching_data(; nlevdecomp=nlevdecomp, nlevsoi=nlevsoi)

        CLM.n_leaching!(d.nf, d.ns, d.params;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:d.nc,
            nlevdecomp=nlevdecomp,
            nlevsoi=nlevsoi,
            dt=d.dt,
            h2osoi_liq=d.h2osoi_liq,
            qflx_drain=d.qflx_drain,
            qflx_surf=d.qflx_surf,
            col_dz=d.col_dz,
            zisoi=d.zisoi,
            use_nitrif_denitrif=false)

        for j in 1:nlevdecomp
            for c in 1:d.nc
                @test d.nf.sminn_leached_vr_col[c, j] >= 0.0
                @test isfinite(d.nf.sminn_leached_vr_col[c, j])
            end
        end
    end

end
