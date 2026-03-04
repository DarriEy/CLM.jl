@testset "LitterVertTranspMod" begin

    # ---------------------------------------------------------------
    # Helper: build minimal test data structures for vertical transport
    # ---------------------------------------------------------------
    function make_test_data(; nc=3, nlevdecomp=5, ndecomp_pools=7, ndecomp_pools_cwd=[5])
        # Build soil grid: exponentially spaced nodes (simplified CLM grid)
        scalez = 0.025
        nlevgrnd = nlevdecomp + 5  # extra bedrock levels
        zsoi_vals = zeros(nlevgrnd)
        for j in 1:nlevgrnd
            zsoi_vals[j] = scalez * (exp(0.5 * (j - 0.5)) - 1.0)
        end

        zisoi_vals = zeros(nlevgrnd + 1)
        zisoi_vals[1] = 0.0
        for j in 1:nlevgrnd
            if j < nlevgrnd
                zisoi_vals[j + 1] = 0.5 * (zsoi_vals[j] + zsoi_vals[j + 1])
            else
                zisoi_vals[j + 1] = zsoi_vals[j] + 0.5 * (zsoi_vals[j] - zisoi_vals[j])
            end
        end

        dzsoi_decomp_vals = zeros(nlevdecomp)
        for j in 1:nlevdecomp
            dzsoi_decomp_vals[j] = zisoi_vals[j + 1] - zisoi_vals[j]
        end

        # Carbon state
        cs = CLM.SoilBiogeochemCarbonStateData()
        cs.decomp_cpools_vr_col = zeros(nc, nlevdecomp, ndecomp_pools)
        # Initialize with some concentration profile (decreasing with depth)
        for s in 1:ndecomp_pools
            for j in 1:nlevdecomp
                for c in 1:nc
                    cs.decomp_cpools_vr_col[c, j, s] = 100.0 * exp(-0.5 * j) * (1.0 + 0.1 * s)
                end
            end
        end

        # Carbon flux
        cf = CLM.SoilBiogeochemCarbonFluxData()
        cf.decomp_cpools_sourcesink_col = zeros(nc, nlevdecomp, ndecomp_pools)
        cf.decomp_cpools_transport_tendency_col = zeros(nc, nlevdecomp, ndecomp_pools)
        cf.tri_ma_vr = zeros(nc, 1)

        # Nitrogen state
        ns = CLM.SoilBiogeochemNitrogenStateData()
        ns.decomp_npools_vr_col = zeros(nc, nlevdecomp, ndecomp_pools)
        for s in 1:ndecomp_pools
            for j in 1:nlevdecomp
                for c in 1:nc
                    ns.decomp_npools_vr_col[c, j, s] = cs.decomp_cpools_vr_col[c, j, s] / 15.0
                end
            end
        end

        # Nitrogen flux
        nf = CLM.SoilBiogeochemNitrogenFluxData()
        nf.decomp_npools_sourcesink_col = zeros(nc, nlevdecomp, ndecomp_pools)
        nf.decomp_npools_transport_tendency_col = zeros(nc, nlevdecomp, ndecomp_pools)

        # SoilBiogeochemState
        st = CLM.SoilBiogeochemStateData()
        st.som_adv_coef_col = zeros(nc, nlevdecomp + 1)
        st.som_diffus_coef_col = zeros(nc, nlevdecomp + 1)

        # Column data
        col = CLM.ColumnData()
        col.nbedrock = fill(nlevdecomp, nc)  # bedrock at bottom of decomp layers
        col.gridcell = ones(Int, nc)

        # Gridcell data
        grc = CLM.GridcellData()
        grc.latdeg = fill(45.0, 1)

        # Active layer data
        active_layer = CLM.ActiveLayerData()
        active_layer.altmax_col = fill(1.0, nc)
        active_layer.altmax_lastyear_col = fill(1.0, nc)

        # Decomp cascade configuration
        cascade_con = CLM.DecompCascadeConData()
        is_cwd = falses(ndecomp_pools)
        for idx in ndecomp_pools_cwd
            if idx <= ndecomp_pools
                is_cwd[idx] = true
            end
        end
        cascade_con.is_cwd = is_cwd
        cascade_con.spinup_factor = ones(ndecomp_pools)

        # Parameters
        params = CLM.LitterVertTranspParams()
        params.som_diffus = 1.0e-4            # m^2/s bioturbation diffusivity
        params.cryoturb_diffusion_k = 5.0e-4  # m^2/s cryoturbation
        params.max_altdepth_cryoturbation = 2.0  # m

        # Mask
        mask_bgc_soilc = trues(nc)

        bounds = 1:nc
        dtime = 1800.0  # 30 min timestep

        return (cs=cs, cf=cf, ns=ns, nf=nf, st=st, col=col, grc=grc,
                cascade_con=cascade_con, active_layer=active_layer,
                params=params, mask_bgc_soilc=mask_bgc_soilc,
                bounds=bounds, dtime=dtime,
                nlevdecomp=nlevdecomp, ndecomp_pools=ndecomp_pools,
                zsoi_vals=zsoi_vals, dzsoi_decomp_vals=dzsoi_decomp_vals,
                zisoi_vals=zisoi_vals)
    end

    # ---------------------------------------------------------------
    @testset "LitterVertTranspParams construction and readparams" begin
        params = CLM.LitterVertTranspParams()
        @test params.som_diffus == 0.0
        @test params.cryoturb_diffusion_k == 0.0
        @test params.max_altdepth_cryoturbation == 0.0

        CLM.litter_vert_transp_readparams!(params;
            som_diffus = 1.0e-4,
            cryoturb_diffusion_k = 5.0e-4,
            max_altdepth_cryoturbation = 2.0)

        @test params.som_diffus == 1.0e-4
        @test params.cryoturb_diffusion_k == 5.0e-4
        @test params.max_altdepth_cryoturbation == 2.0
    end

    # ---------------------------------------------------------------
    @testset "patankar_A function" begin
        # pe = 0 -> A = max(0, 1^5) = 1
        @test CLM.patankar_A(0.0) == 1.0

        # pe = 10 -> A = max(0, (1 - 1)^5) = 0
        @test CLM.patankar_A(10.0) == 0.0

        # pe = -10 -> A = max(0, (1 - 1)^5) = 0
        @test CLM.patankar_A(-10.0) == 0.0

        # pe = 5 -> A = max(0, (1 - 0.5)^5) = 0.5^5 = 0.03125
        @test CLM.patankar_A(5.0) == 0.5^5

        # pe = 20 -> A = max(0, (1 - 2)^5) = max(0, -1) = 0
        @test CLM.patankar_A(20.0) == 0.0

        # pe = 2 -> A = max(0, (1 - 0.2)^5) = 0.8^5
        @test CLM.patankar_A(2.0) ≈ 0.8^5
    end

    # ---------------------------------------------------------------
    @testset "runs without error - bioturbation regime" begin
        td = make_test_data()

        # altmax > max_altdepth_cryoturbation triggers bioturbation regime
        td.active_layer.altmax_col .= 3.0
        td.active_layer.altmax_lastyear_col .= 3.0

        # Save initial C state for comparison
        c_init = copy(td.cs.decomp_cpools_vr_col)
        n_init = copy(td.ns.decomp_npools_vr_col)

        CLM.litter_vert_transp!(
            td.cs, td.cf, td.ns, td.nf, td.st, td.col, td.grc,
            td.cascade_con, td.active_layer, td.params;
            mask_bgc_soilc = td.mask_bgc_soilc,
            bounds = td.bounds,
            dtime = td.dtime,
            nlevdecomp = td.nlevdecomp,
            ndecomp_pools = td.ndecomp_pools,
            zsoi_vals = td.zsoi_vals,
            dzsoi_decomp_vals = td.dzsoi_decomp_vals,
            zisoi_vals = td.zisoi_vals,
            spinup_state = 0,
            use_soil_matrixcn = false,
            som_adv_flux_val = 0.0,
            max_depth_cryoturb_val = 3.0)

        # Verify non-CWD pools were modified by transport
        for s in 1:td.ndecomp_pools
            if !td.cascade_con.is_cwd[s]
                # Concentrations should have been updated (transport mixes them)
                @test td.cs.decomp_cpools_vr_col[:, :, s] != c_init[:, :, s] ||
                      all(td.cs.decomp_cpools_vr_col[:, :, s] .== 0.0)
            end
        end

        # CWD pool (pool 5) should just have source added (source is zero here)
        cwd_idx = 5
        @test td.cs.decomp_cpools_vr_col[:, :, cwd_idx] ≈ c_init[:, :, cwd_idx]

        # Nitrogen should also have been modified
        for s in 1:td.ndecomp_pools
            if !td.cascade_con.is_cwd[s]
                @test td.ns.decomp_npools_vr_col[:, :, s] != n_init[:, :, s] ||
                      all(td.ns.decomp_npools_vr_col[:, :, s] .== 0.0)
            end
        end
    end

    # ---------------------------------------------------------------
    @testset "runs without error - cryoturbation regime" begin
        td = make_test_data()

        # altmax < max_altdepth_cryoturbation triggers cryoturbation regime
        td.active_layer.altmax_col .= 0.5
        td.active_layer.altmax_lastyear_col .= 0.5
        td.params.max_altdepth_cryoturbation = 2.0

        c_init = copy(td.cs.decomp_cpools_vr_col)

        CLM.litter_vert_transp!(
            td.cs, td.cf, td.ns, td.nf, td.st, td.col, td.grc,
            td.cascade_con, td.active_layer, td.params;
            mask_bgc_soilc = td.mask_bgc_soilc,
            bounds = td.bounds,
            dtime = td.dtime,
            nlevdecomp = td.nlevdecomp,
            ndecomp_pools = td.ndecomp_pools,
            zsoi_vals = td.zsoi_vals,
            dzsoi_decomp_vals = td.dzsoi_decomp_vals,
            zisoi_vals = td.zisoi_vals,
            spinup_state = 0,
            use_soil_matrixcn = false,
            som_adv_flux_val = 0.0,
            max_depth_cryoturb_val = 3.0)

        # Diffusivity coefficients should reflect cryoturbation pattern
        for c in 1:3
            # In cryoturbation regime, som_diffus_coef should be set
            @test td.st.som_diffus_coef_col[c, 1] ≈ td.params.cryoturb_diffusion_k
            @test td.st.som_adv_coef_col[c, 1] == 0.0
        end
    end

    # ---------------------------------------------------------------
    @testset "frozen soil - no mixing" begin
        td = make_test_data()

        # altmax = 0 triggers frozen soil regime (no mixing)
        td.active_layer.altmax_col .= 0.0
        td.active_layer.altmax_lastyear_col .= 0.0

        CLM.litter_vert_transp!(
            td.cs, td.cf, td.ns, td.nf, td.st, td.col, td.grc,
            td.cascade_con, td.active_layer, td.params;
            mask_bgc_soilc = td.mask_bgc_soilc,
            bounds = td.bounds,
            dtime = td.dtime,
            nlevdecomp = td.nlevdecomp,
            ndecomp_pools = td.ndecomp_pools,
            zsoi_vals = td.zsoi_vals,
            dzsoi_decomp_vals = td.dzsoi_decomp_vals,
            zisoi_vals = td.zisoi_vals,
            spinup_state = 0,
            use_soil_matrixcn = false,
            som_adv_flux_val = 0.0,
            max_depth_cryoturb_val = 3.0)

        # All diffusivity and advection coefficients should be zero
        for c in 1:3
            for j in 1:(td.nlevdecomp + 1)
                @test td.st.som_diffus_coef_col[c, j] == 0.0
                @test td.st.som_adv_coef_col[c, j] == 0.0
            end
        end
    end

    # ---------------------------------------------------------------
    @testset "mass conservation - no source/sink" begin
        nc = 2
        nlevdecomp = 5
        ndecomp_pools = 4
        td = make_test_data(nc=nc, nlevdecomp=nlevdecomp, ndecomp_pools=ndecomp_pools,
                            ndecomp_pools_cwd=[4])

        # Bioturbation regime
        td.active_layer.altmax_col .= 3.0
        td.active_layer.altmax_lastyear_col .= 3.0

        # Zero source/sink
        td.cf.decomp_cpools_sourcesink_col .= 0.0
        td.nf.decomp_npools_sourcesink_col .= 0.0

        # Compute initial total mass per pool (integrated over depth)
        c_total_init = zeros(nc, ndecomp_pools)
        n_total_init = zeros(nc, ndecomp_pools)
        for s in 1:ndecomp_pools
            for j in 1:nlevdecomp
                for c in 1:nc
                    c_total_init[c, s] += td.cs.decomp_cpools_vr_col[c, j, s] * td.dzsoi_decomp_vals[j]
                    n_total_init[c, s] += td.ns.decomp_npools_vr_col[c, j, s] * td.dzsoi_decomp_vals[j]
                end
            end
        end

        CLM.litter_vert_transp!(
            td.cs, td.cf, td.ns, td.nf, td.st, td.col, td.grc,
            td.cascade_con, td.active_layer, td.params;
            mask_bgc_soilc = td.mask_bgc_soilc,
            bounds = td.bounds,
            dtime = td.dtime,
            nlevdecomp = nlevdecomp,
            ndecomp_pools = ndecomp_pools,
            zsoi_vals = td.zsoi_vals,
            dzsoi_decomp_vals = td.dzsoi_decomp_vals,
            zisoi_vals = td.zisoi_vals,
            spinup_state = 0,
            use_soil_matrixcn = false,
            som_adv_flux_val = 0.0,
            max_depth_cryoturb_val = 3.0)

        # Compute final total mass per pool
        c_total_final = zeros(nc, ndecomp_pools)
        n_total_final = zeros(nc, ndecomp_pools)
        for s in 1:ndecomp_pools
            for j in 1:nlevdecomp
                for c in 1:nc
                    c_total_final[c, s] += td.cs.decomp_cpools_vr_col[c, j, s] * td.dzsoi_decomp_vals[j]
                    n_total_final[c, s] += td.ns.decomp_npools_vr_col[c, j, s] * td.dzsoi_decomp_vals[j]
                end
            end
        end

        # Total mass should be conserved (within tolerance) for non-CWD pools
        for s in 1:ndecomp_pools
            for c in 1:nc
                if c_total_init[c, s] > 0.0
                    @test c_total_final[c, s] ≈ c_total_init[c, s] rtol=1e-6
                    @test n_total_final[c, s] ≈ n_total_init[c, s] rtol=1e-6
                end
            end
        end
    end

    # ---------------------------------------------------------------
    @testset "CWD pools only add source" begin
        td = make_test_data(nc=2, nlevdecomp=5, ndecomp_pools=4,
                            ndecomp_pools_cwd=[3, 4])  # pools 3,4 are CWD

        td.active_layer.altmax_col .= 3.0
        td.active_layer.altmax_lastyear_col .= 3.0

        # Set a nonzero source for CWD pools
        for s in [3, 4]
            for j in 1:5
                for c in 1:2
                    td.cf.decomp_cpools_sourcesink_col[c, j, s] = 0.5
                    td.nf.decomp_npools_sourcesink_col[c, j, s] = 0.05
                end
            end
        end

        c_init_cwd3 = copy(td.cs.decomp_cpools_vr_col[:, :, 3])
        c_init_cwd4 = copy(td.cs.decomp_cpools_vr_col[:, :, 4])

        CLM.litter_vert_transp!(
            td.cs, td.cf, td.ns, td.nf, td.st, td.col, td.grc,
            td.cascade_con, td.active_layer, td.params;
            mask_bgc_soilc = td.mask_bgc_soilc,
            bounds = td.bounds,
            dtime = td.dtime,
            nlevdecomp = 5,
            ndecomp_pools = 4,
            zsoi_vals = td.zsoi_vals,
            dzsoi_decomp_vals = td.dzsoi_decomp_vals,
            zisoi_vals = td.zisoi_vals,
            spinup_state = 0,
            use_soil_matrixcn = false,
            som_adv_flux_val = 0.0,
            max_depth_cryoturb_val = 3.0)

        # CWD pools should have initial + source
        for j in 1:5
            for c in 1:2
                @test td.cs.decomp_cpools_vr_col[c, j, 3] ≈ c_init_cwd3[c, j] + 0.5
                @test td.cs.decomp_cpools_vr_col[c, j, 4] ≈ c_init_cwd4[c, j] + 0.5
            end
        end
    end

    # ---------------------------------------------------------------
    @testset "masked columns are skipped" begin
        td = make_test_data(nc=3)

        td.active_layer.altmax_col .= 3.0
        td.active_layer.altmax_lastyear_col .= 3.0

        # Mask out column 2
        td.mask_bgc_soilc[2] = false

        c_init = copy(td.cs.decomp_cpools_vr_col)
        n_init = copy(td.ns.decomp_npools_vr_col)

        CLM.litter_vert_transp!(
            td.cs, td.cf, td.ns, td.nf, td.st, td.col, td.grc,
            td.cascade_con, td.active_layer, td.params;
            mask_bgc_soilc = td.mask_bgc_soilc,
            bounds = td.bounds,
            dtime = td.dtime,
            nlevdecomp = td.nlevdecomp,
            ndecomp_pools = td.ndecomp_pools,
            zsoi_vals = td.zsoi_vals,
            dzsoi_decomp_vals = td.dzsoi_decomp_vals,
            zisoi_vals = td.zisoi_vals,
            spinup_state = 0,
            use_soil_matrixcn = false,
            som_adv_flux_val = 0.0,
            max_depth_cryoturb_val = 3.0)

        # Column 2 should be unchanged
        @test td.cs.decomp_cpools_vr_col[2, :, :] == c_init[2, :, :]
        @test td.ns.decomp_npools_vr_col[2, :, :] == n_init[2, :, :]

        # Columns 1 and 3 should have been modified (for non-CWD pools)
        for s in 1:td.ndecomp_pools
            if !td.cascade_con.is_cwd[s]
                @test td.cs.decomp_cpools_vr_col[1, :, s] != c_init[1, :, s]
                @test td.cs.decomp_cpools_vr_col[3, :, s] != c_init[3, :, s]
            end
        end
    end

    # ---------------------------------------------------------------
    @testset "transport tendency is consistent" begin
        td = make_test_data(nc=1, nlevdecomp=5, ndecomp_pools=3,
                            ndecomp_pools_cwd=[3])

        td.active_layer.altmax_col .= 3.0
        td.active_layer.altmax_lastyear_col .= 3.0

        # Zero source/sink
        td.cf.decomp_cpools_sourcesink_col .= 0.0
        td.nf.decomp_npools_sourcesink_col .= 0.0

        c_init = copy(td.cs.decomp_cpools_vr_col)

        CLM.litter_vert_transp!(
            td.cs, td.cf, td.ns, td.nf, td.st, td.col, td.grc,
            td.cascade_con, td.active_layer, td.params;
            mask_bgc_soilc = td.mask_bgc_soilc,
            bounds = td.bounds,
            dtime = td.dtime,
            nlevdecomp = 5,
            ndecomp_pools = 3,
            zsoi_vals = td.zsoi_vals,
            dzsoi_decomp_vals = td.dzsoi_decomp_vals,
            zisoi_vals = td.zisoi_vals,
            spinup_state = 0,
            use_soil_matrixcn = false,
            som_adv_flux_val = 0.0,
            max_depth_cryoturb_val = 3.0)

        # For non-CWD pools, tendency = (new_conc - old_conc - source) / dtime
        # Since source = 0, tendency = (new_conc - old_conc) / dtime
        for s in 1:2  # non-CWD
            for j in 1:5
                expected_tendency = (td.cs.decomp_cpools_vr_col[1, j, s] - c_init[1, j, s]) / td.dtime
                @test td.cf.decomp_cpools_transport_tendency_col[1, j, s] ≈ expected_tendency rtol=1e-10
            end
        end
    end

    # ---------------------------------------------------------------
    @testset "diffusivity coefficients - cryoturbation profile" begin
        td = make_test_data(nc=1, nlevdecomp=10, ndecomp_pools=3,
                            ndecomp_pools_cwd=[3])

        # Set active layer depth well below cryoturbation threshold
        td.active_layer.altmax_col .= 0.3
        td.active_layer.altmax_lastyear_col .= 0.3
        td.params.max_altdepth_cryoturbation = 2.0
        td.params.cryoturb_diffusion_k = 5.0e-4

        CLM.litter_vert_transp!(
            td.cs, td.cf, td.ns, td.nf, td.st, td.col, td.grc,
            td.cascade_con, td.active_layer, td.params;
            mask_bgc_soilc = td.mask_bgc_soilc,
            bounds = td.bounds,
            dtime = td.dtime,
            nlevdecomp = 10,
            ndecomp_pools = 3,
            zsoi_vals = td.zsoi_vals,
            dzsoi_decomp_vals = td.dzsoi_decomp_vals,
            zisoi_vals = td.zisoi_vals,
            spinup_state = 0,
            use_soil_matrixcn = false,
            som_adv_flux_val = 0.0,
            max_depth_cryoturb_val = 3.0)

        # Above active layer: diffusivity should equal cryoturb_diffusion_k
        # Below: should decrease linearly to zero
        # Note: zisoi_vals[j+1] = Fortran zisoi(j) = interface depth at position j
        alt_max_val = 0.3
        for j in 1:11
            if td.zisoi_vals[j + 1] < alt_max_val && j <= td.col.nbedrock[1] + 1
                @test td.st.som_diffus_coef_col[1, j] ≈ 5.0e-4
            end
            # Advection should be zero in cryoturbation regime
            @test td.st.som_adv_coef_col[1, j] == 0.0
        end
    end

    # ---------------------------------------------------------------
    @testset "bedrock correction redistributes mass" begin
        nc = 1
        nlevdecomp = 5
        ndecomp_pools = 2
        td = make_test_data(nc=nc, nlevdecomp=nlevdecomp, ndecomp_pools=ndecomp_pools,
                            ndecomp_pools_cwd=Int[])

        # Set bedrock at level 3 (so levels 4,5 are below bedrock)
        td.col.nbedrock[1] = 3

        # Give some initial concentration below bedrock
        td.cs.decomp_cpools_vr_col[1, 4, 1] = 10.0
        td.cs.decomp_cpools_vr_col[1, 5, 1] = 5.0
        td.ns.decomp_npools_vr_col[1, 4, 1] = 1.0
        td.ns.decomp_npools_vr_col[1, 5, 1] = 0.5

        td.active_layer.altmax_col .= 3.0
        td.active_layer.altmax_lastyear_col .= 3.0

        CLM.litter_vert_transp!(
            td.cs, td.cf, td.ns, td.nf, td.st, td.col, td.grc,
            td.cascade_con, td.active_layer, td.params;
            mask_bgc_soilc = td.mask_bgc_soilc,
            bounds = td.bounds,
            dtime = td.dtime,
            nlevdecomp = nlevdecomp,
            ndecomp_pools = ndecomp_pools,
            zsoi_vals = td.zsoi_vals,
            dzsoi_decomp_vals = td.dzsoi_decomp_vals,
            zisoi_vals = td.zisoi_vals,
            spinup_state = 0,
            use_soil_matrixcn = false,
            som_adv_flux_val = 0.0,
            max_depth_cryoturb_val = 3.0)

        # Concentrations below bedrock should be zero after correction
        @test td.cs.decomp_cpools_vr_col[1, 4, 1] == 0.0
        @test td.cs.decomp_cpools_vr_col[1, 5, 1] == 0.0
        @test td.ns.decomp_npools_vr_col[1, 4, 1] == 0.0
        @test td.ns.decomp_npools_vr_col[1, 5, 1] == 0.0
    end

    # ---------------------------------------------------------------
    @testset "source/sink affects final concentrations" begin
        nc = 1
        nlevdecomp = 5
        ndecomp_pools = 2
        td = make_test_data(nc=nc, nlevdecomp=nlevdecomp, ndecomp_pools=ndecomp_pools,
                            ndecomp_pools_cwd=Int[])

        td.active_layer.altmax_col .= 3.0
        td.active_layer.altmax_lastyear_col .= 3.0

        # Add a significant positive source term
        for j in 1:nlevdecomp
            td.cf.decomp_cpools_sourcesink_col[1, j, 1] = 5.0
        end

        c_init = copy(td.cs.decomp_cpools_vr_col)

        CLM.litter_vert_transp!(
            td.cs, td.cf, td.ns, td.nf, td.st, td.col, td.grc,
            td.cascade_con, td.active_layer, td.params;
            mask_bgc_soilc = td.mask_bgc_soilc,
            bounds = td.bounds,
            dtime = td.dtime,
            nlevdecomp = nlevdecomp,
            ndecomp_pools = ndecomp_pools,
            zsoi_vals = td.zsoi_vals,
            dzsoi_decomp_vals = td.dzsoi_decomp_vals,
            zisoi_vals = td.zisoi_vals,
            spinup_state = 0,
            use_soil_matrixcn = false,
            som_adv_flux_val = 0.0,
            max_depth_cryoturb_val = 3.0)

        # With positive source, concentrations should generally increase
        # (exact values depend on diffusion, but total mass should reflect added source)
        c_total_final = 0.0
        c_total_init = 0.0
        source_total = 0.0
        for j in 1:nlevdecomp
            c_total_final += td.cs.decomp_cpools_vr_col[1, j, 1] * td.dzsoi_decomp_vals[j]
            c_total_init += c_init[1, j, 1] * td.dzsoi_decomp_vals[j]
            source_total += 5.0 * td.dzsoi_decomp_vals[j]
        end

        # Final total should approximately equal initial + source
        @test c_total_final ≈ (c_total_init + source_total) rtol=1e-5
    end

    # ---------------------------------------------------------------
    @testset "single-column single-level minimal case" begin
        nc = 1
        nlevdecomp = 1
        ndecomp_pools = 1
        nlevgrnd = 6

        scalez = 0.025
        zsoi_vals = zeros(nlevgrnd)
        for j in 1:nlevgrnd
            zsoi_vals[j] = scalez * (exp(0.5 * (j - 0.5)) - 1.0)
        end
        zisoi_vals = zeros(nlevgrnd + 1)
        zisoi_vals[1] = 0.0
        for j in 1:nlevgrnd
            if j < nlevgrnd
                zisoi_vals[j + 1] = 0.5 * (zsoi_vals[j] + zsoi_vals[j + 1])
            else
                zisoi_vals[j + 1] = zsoi_vals[j] + 0.5 * (zsoi_vals[j] - zisoi_vals[j])
            end
        end
        dzsoi_decomp_vals = [zisoi_vals[2] - zisoi_vals[1]]

        cs = CLM.SoilBiogeochemCarbonStateData()
        cs.decomp_cpools_vr_col = fill(50.0, nc, nlevdecomp, ndecomp_pools)
        cf = CLM.SoilBiogeochemCarbonFluxData()
        cf.decomp_cpools_sourcesink_col = zeros(nc, nlevdecomp, ndecomp_pools)
        cf.decomp_cpools_transport_tendency_col = zeros(nc, nlevdecomp, ndecomp_pools)
        cf.tri_ma_vr = zeros(nc, 1)
        ns = CLM.SoilBiogeochemNitrogenStateData()
        ns.decomp_npools_vr_col = fill(5.0, nc, nlevdecomp, ndecomp_pools)
        nf = CLM.SoilBiogeochemNitrogenFluxData()
        nf.decomp_npools_sourcesink_col = zeros(nc, nlevdecomp, ndecomp_pools)
        nf.decomp_npools_transport_tendency_col = zeros(nc, nlevdecomp, ndecomp_pools)

        st = CLM.SoilBiogeochemStateData()
        st.som_adv_coef_col = zeros(nc, nlevdecomp + 1)
        st.som_diffus_coef_col = zeros(nc, nlevdecomp + 1)

        col = CLM.ColumnData()
        col.nbedrock = [nlevdecomp]
        col.gridcell = [1]

        grc = CLM.GridcellData()
        grc.latdeg = [60.0]

        active_layer = CLM.ActiveLayerData()
        active_layer.altmax_col = [1.0]
        active_layer.altmax_lastyear_col = [1.0]

        cascade_con = CLM.DecompCascadeConData()
        cascade_con.is_cwd = falses(ndecomp_pools)
        cascade_con.spinup_factor = ones(ndecomp_pools)

        params = CLM.LitterVertTranspParams(
            som_diffus=1.0e-4,
            cryoturb_diffusion_k=5.0e-4,
            max_altdepth_cryoturbation=2.0)

        mask_bgc_soilc = trues(nc)

        # Should run without error
        CLM.litter_vert_transp!(
            cs, cf, ns, nf, st, col, grc,
            cascade_con, active_layer, params;
            mask_bgc_soilc = mask_bgc_soilc,
            bounds = 1:nc,
            dtime = 3600.0,
            nlevdecomp = nlevdecomp,
            ndecomp_pools = ndecomp_pools,
            zsoi_vals = zsoi_vals,
            dzsoi_decomp_vals = dzsoi_decomp_vals,
            zisoi_vals = zisoi_vals,
            spinup_state = 0,
            use_soil_matrixcn = false,
            som_adv_flux_val = 0.0,
            max_depth_cryoturb_val = 3.0)

        # With nlevdecomp=1, transport should not change concentration much
        # (only 1 level, boundary conditions dominate)
        @test cs.decomp_cpools_vr_col[1, 1, 1] >= 0.0
        @test ns.decomp_npools_vr_col[1, 1, 1] >= 0.0
    end

    # ---------------------------------------------------------------
    @testset "diffusion smooths profiles" begin
        # Test that diffusion tends to smooth out concentration gradients
        nc = 1
        nlevdecomp = 10
        ndecomp_pools = 1
        td = make_test_data(nc=nc, nlevdecomp=nlevdecomp, ndecomp_pools=ndecomp_pools,
                            ndecomp_pools_cwd=Int[])

        td.active_layer.altmax_col .= 5.0
        td.active_layer.altmax_lastyear_col .= 5.0
        td.col.nbedrock .= nlevdecomp

        # Set a sharp concentration spike at level 3
        td.cs.decomp_cpools_vr_col .= 0.0
        td.cs.decomp_cpools_vr_col[1, 3, 1] = 100.0
        td.ns.decomp_npools_vr_col .= 0.0
        td.ns.decomp_npools_vr_col[1, 3, 1] = 10.0

        # Use a relatively large diffusivity
        td.params.som_diffus = 1.0e-2

        CLM.litter_vert_transp!(
            td.cs, td.cf, td.ns, td.nf, td.st, td.col, td.grc,
            td.cascade_con, td.active_layer, td.params;
            mask_bgc_soilc = td.mask_bgc_soilc,
            bounds = td.bounds,
            dtime = td.dtime,
            nlevdecomp = nlevdecomp,
            ndecomp_pools = ndecomp_pools,
            zsoi_vals = td.zsoi_vals,
            dzsoi_decomp_vals = td.dzsoi_decomp_vals,
            zisoi_vals = td.zisoi_vals,
            spinup_state = 0,
            use_soil_matrixcn = false,
            som_adv_flux_val = 0.0,
            max_depth_cryoturb_val = 3.0)

        # After diffusion, concentration at level 3 should have decreased
        @test td.cs.decomp_cpools_vr_col[1, 3, 1] < 100.0

        # And neighboring levels should have gained concentration
        @test td.cs.decomp_cpools_vr_col[1, 2, 1] > 0.0 || td.cs.decomp_cpools_vr_col[1, 4, 1] > 0.0
    end

    # ---------------------------------------------------------------
    @testset "all concentrations remain non-negative" begin
        td = make_test_data(nc=2, nlevdecomp=5, ndecomp_pools=3,
                            ndecomp_pools_cwd=[3])

        td.active_layer.altmax_col .= 3.0
        td.active_layer.altmax_lastyear_col .= 3.0

        # Use small initial concentrations with some source/sink
        td.cs.decomp_cpools_vr_col .= 1.0
        td.ns.decomp_npools_vr_col .= 0.1
        td.cf.decomp_cpools_sourcesink_col .= 0.0
        td.nf.decomp_npools_sourcesink_col .= 0.0

        CLM.litter_vert_transp!(
            td.cs, td.cf, td.ns, td.nf, td.st, td.col, td.grc,
            td.cascade_con, td.active_layer, td.params;
            mask_bgc_soilc = td.mask_bgc_soilc,
            bounds = td.bounds,
            dtime = td.dtime,
            nlevdecomp = 5,
            ndecomp_pools = 3,
            zsoi_vals = td.zsoi_vals,
            dzsoi_decomp_vals = td.dzsoi_decomp_vals,
            zisoi_vals = td.zisoi_vals,
            spinup_state = 0,
            use_soil_matrixcn = false,
            som_adv_flux_val = 0.0,
            max_depth_cryoturb_val = 3.0)

        # All concentrations should remain non-negative (or very close)
        @test all(td.cs.decomp_cpools_vr_col .>= -1e-10)
        @test all(td.ns.decomp_npools_vr_col .>= -1e-10)
    end

end
