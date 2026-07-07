@testset "dwt → persistent soil-matrix B-input (matrix_Cinput/Ninput_col)" begin
    # In matrix mode the dynamic-landcover dwt→litter/CWD inputs must NOT be added
    # directly to the decomp pools (the pools stay at start-of-step for the soil-matrix
    # solve); instead they accumulate into the persistent soilbgc_cf/nf.matrix_Cinput/
    # Ninput_col, which cn_soil_matrix_advance! folds into the B-input and then zeroes.
    nc = 1; ng = 1; nlevdecomp = 2; ndecomp_pools = 3; dt = 1800.0
    i_litr_min = 1; i_litr_max = 2; i_cwd = 3; ndp_vr = ndecomp_pools * nlevdecomp

    function mk()
        csv = CLM.CNVegCarbonStateData(); csv.seedc_grc = zeros(ng)
        cfv = CLM.CNVegCarbonFluxData()
        cfv.dwt_frootc_to_litr_c_col   = zeros(nc, nlevdecomp, ndecomp_pools)
        cfv.dwt_livecrootc_to_cwdc_col = zeros(nc, nlevdecomp)
        cfv.dwt_deadcrootc_to_cwdc_col = zeros(nc, nlevdecomp)
        cfv.dwt_seedc_to_leaf_grc = zeros(ng); cfv.dwt_seedc_to_deadstem_grc = zeros(ng)
        for j in 1:nlevdecomp, i in i_litr_min:i_litr_max
            cfv.dwt_frootc_to_litr_c_col[1, j, i] = 0.5 * i + 0.1 * j
        end
        for j in 1:nlevdecomp
            cfv.dwt_livecrootc_to_cwdc_col[1, j] = 0.3 + 0.05 * j
            cfv.dwt_deadcrootc_to_cwdc_col[1, j] = 0.2 + 0.02 * j
        end
        css = CLM.SoilBiogeochemCarbonStateData()
        css.decomp_cpools_vr_col = fill(100.0, nc, nlevdecomp, ndecomp_pools)
        cfs = CLM.SoilBiogeochemCarbonFluxData()
        cfs.matrix_Cinput_col = zeros(nc, ndp_vr)
        return (csv=csv, cfv=cfv, css=css, cfs=cfs)
    end

    runc(d; mx) = CLM.c_state_update_dyn_patch!(d.csv, d.cfv, d.css;
        mask_soilc_with_inactive = [true], bounds_col = 1:nc, bounds_grc = 1:ng,
        nlevdecomp = nlevdecomp, i_litr_min = i_litr_min, i_litr_max = i_litr_max,
        i_cwd = i_cwd, use_soil_matrixcn = mx, cf_soil = d.cfs, dt = dt)

    # expected B-input from dwt
    expC = zeros(nc, ndp_vr)
    let d0 = mk()
        for j in 1:nlevdecomp
            for i in i_litr_min:i_litr_max
                expC[1, j + (i - 1) * nlevdecomp] = d0.cfv.dwt_frootc_to_litr_c_col[1, j, i] * dt
            end
            expC[1, j + (i_cwd - 1) * nlevdecomp] =
                (d0.cfv.dwt_livecrootc_to_cwdc_col[1, j] + d0.cfv.dwt_deadcrootc_to_cwdc_col[1, j]) * dt
        end
    end

    # --- matrix mode: dwt → matrix_Cinput_col, pools untouched ---
    dm = mk(); pools0 = copy(dm.css.decomp_cpools_vr_col)
    runc(dm; mx = true)
    @test dm.cfs.matrix_Cinput_col ≈ expC atol = 1e-10
    @test dm.css.decomp_cpools_vr_col == pools0          # pools untouched in matrix mode
    # accumulates (second call doubles) — the field is persistent within a step.
    runc(dm; mx = true)
    @test dm.cfs.matrix_Cinput_col ≈ 2 .* expC atol = 1e-10

    # --- non-matrix mode: dwt → pools directly, matrix_Cinput_col untouched ---
    dn = mk()
    runc(dn; mx = false)
    @test all(==(0.0), dn.cfs.matrix_Cinput_col)         # B-input untouched
    @test dn.css.decomp_cpools_vr_col[1, 1, 1] ≈ 100.0 + dn.cfv.dwt_frootc_to_litr_c_col[1, 1, 1] * dt
    @test dn.css.decomp_cpools_vr_col[1, 1, i_cwd] ≈
        100.0 + (dn.cfv.dwt_livecrootc_to_cwdc_col[1, 1] + dn.cfv.dwt_deadcrootc_to_cwdc_col[1, 1]) * dt

    # --- cn_soil_matrix_advance! folds in the persistent dwt B-input then zeroes it ---
    dm2 = mk(); runc(dm2; mx = true)      # matrix_Cinput_col now holds the dwt B-input
    css2 = dm2.css; nss2 = CLM.SoilBiogeochemNitrogenStateData()
    nss2.decomp_npools_vr_col = fill(100.0 / 15.0, nc, nlevdecomp, ndecomp_pools)
    cfs2 = dm2.cfs
    cfs2.rf_decomp_cascade_col = zeros(nc, nlevdecomp, ndecomp_pools)
    cfs2.pathfrac_decomp_cascade_col = ones(nc, nlevdecomp, ndecomp_pools)
    cfs2.tri_ma_vr = zeros(nc, (nlevdecomp * 3 - 2) * (ndecomp_pools - 1))
    nfs2 = CLM.SoilBiogeochemNitrogenFluxData(); nfs2.matrix_Ninput_col = zeros(nc, ndp_vr)
    vcf = CLM.CNVegCarbonFluxData(); CLM.cnveg_carbon_flux_init!(vcf, 1, nc, 1;
        nlevdecomp_full=nlevdecomp, ndecomp_pools=ndecomp_pools)
    vnf = CLM.CNVegNitrogenFluxData(); CLM.cnveg_nitrogen_flux_init!(vnf, 1, nc, 1;
        nlevdecomp_full=nlevdecomp, ndecomp_pools=ndecomp_pools, i_litr_max=i_litr_max)
    for f in (:phenology_c_to_litr_c_col, :gap_mortality_c_to_litr_c_col, :m_c_to_litr_fire_col,
              :gap_mortality_c_to_cwdc_col, :fire_mortality_c_to_cwdc_col); fill!(getfield(vcf, f), 0.0); end
    for f in (:phenology_n_to_litr_n_col, :gap_mortality_n_to_litr_n_col, :m_n_to_litr_fire_col,
              :gap_mortality_n_to_cwdn_col, :fire_mortality_n_to_cwdn_col); fill!(getfield(vnf, f), 0.0); end
    cc = CLM.DecompCascadeConData(cascade_donor_pool=[1,2,3], cascade_receiver_pool=[2,0,1],
        floating_cn_ratio_decomp_pools=BitVector([false,false,false]),
        is_cwd=BitVector([false,false,true]), initial_cn_ratio=[20.0,15.0,90.0])
    X0 = copy(css2.decomp_cpools_vr_col)
    CLM.cn_soil_matrix_advance!(cc, css2, nss2, cfs2, vcf, vnf; Ksoil=zeros(nc, ndp_vr),
        mask_soilc=[true], bounds_col=1:nc, nlevdecomp=nlevdecomp, ndecomp_pools=ndecomp_pools,
        ndecomp_cascade_transitions=3, i_litr_min=i_litr_min, i_litr_max=i_litr_max,
        i_cwd=i_cwd, dt=dt, soilbgc_nf=nfs2)
    # with no cascade/transport/other inputs, the pools gained exactly the dwt B-input.
    for j in 1:nlevdecomp, i in 1:ndecomp_pools
        @test css2.decomp_cpools_vr_col[1, j, i] ≈ X0[1, j, i] + expC[1, j + (i - 1) * nlevdecomp] atol = 1e-9
    end
    @test all(==(0.0), cfs2.matrix_Cinput_col)   # zeroed after the solve
end
