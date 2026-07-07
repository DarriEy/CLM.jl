@testset "cn_soil_matrix_input_accumulate! (B-input from litter/CWD fluxes)" begin
    # Phase-3 soil-matrix B-input: matrix_Cinput/Ninput[c, j+(i-1)*nlev] = Σ (veg→soil
    # litterfall/CWD input fluxes to pool i, level j) · dt. Litter pools get the
    # *_to_litr_* fluxes (indexed by pool); the CWD pool gets the *_to_cwd* fluxes.
    nc = 1; nlevdecomp = 2; ndecomp_pools = 3; dt = 1800.0
    i_litr_min = 1; i_litr_max = 2; i_cwd = 3
    ndecomp_pools_vr = ndecomp_pools * nlevdecomp

    cf = CLM.CNVegCarbonFluxData(); CLM.cnveg_carbon_flux_init!(cf, 1, nc, 1;
        nlevdecomp_full=nlevdecomp, ndecomp_pools=ndecomp_pools)
    nf = CLM.CNVegNitrogenFluxData(); CLM.cnveg_nitrogen_flux_init!(nf, 1, nc, 1;
        nlevdecomp_full=nlevdecomp, ndecomp_pools=ndecomp_pools, i_litr_max=i_litr_max)

    # zero every col field the accumulator reads, then set distinct known values.
    litr_c = (:phenology_c_to_litr_c_col, :dwt_frootc_to_litr_c_col,
              :gap_mortality_c_to_litr_c_col, :m_c_to_litr_fire_col)
    cwd_c  = (:dwt_livecrootc_to_cwdc_col, :dwt_deadcrootc_to_cwdc_col,
              :gap_mortality_c_to_cwdc_col, :fire_mortality_c_to_cwdc_col)
    litr_n = (:phenology_n_to_litr_n_col, :dwt_frootn_to_litr_n_col,
              :gap_mortality_n_to_litr_n_col, :m_n_to_litr_fire_col)
    cwd_n  = (:dwt_livecrootn_to_cwdn_col, :dwt_deadcrootn_to_cwdn_col,
              :gap_mortality_n_to_cwdn_col, :fire_mortality_n_to_cwdn_col)
    for f in litr_c; fill!(getfield(cf, f), 0.0); end
    for f in cwd_c;  fill!(getfield(cf, f), 0.0); end
    for f in litr_n; fill!(getfield(nf, f), 0.0); end
    for f in cwd_n;  fill!(getfield(nf, f), 0.0); end

    # assign distinct values so a mis-mapping (wrong pool/level/field) is caught.
    for (ki, f) in enumerate(litr_c), j in 1:nlevdecomp, i in i_litr_min:i_litr_max
        getfield(cf, f)[1, j, i] = ki * 1.0 + 0.1 * j + 0.01 * i
    end
    for (ki, f) in enumerate(cwd_c), j in 1:nlevdecomp
        getfield(cf, f)[1, j] = 10.0 * ki + j
    end
    for (ki, f) in enumerate(litr_n), j in 1:nlevdecomp, i in i_litr_min:i_litr_max
        getfield(nf, f)[1, j, i] = 0.1 * ki + 0.02 * j + 0.003 * i
    end
    for (ki, f) in enumerate(cwd_n), j in 1:nlevdecomp
        getfield(nf, f)[1, j] = 1.0 * ki + 0.5 * j
    end

    Cin = zeros(nc, ndecomp_pools_vr); Nin = zeros(nc, ndecomp_pools_vr)
    CLM.cn_soil_matrix_input_accumulate!(Cin, Nin, cf, nf;
        mask_soilc=[true], bounds_col=1:nc, nlevdecomp=nlevdecomp,
        i_litr_min=i_litr_min, i_litr_max=i_litr_max, i_cwd=i_cwd, dt=dt,
        transient_landcover=false)

    # independent reference
    Cref = zeros(nc, ndecomp_pools_vr); Nref = zeros(nc, ndecomp_pools_vr)
    for j in 1:nlevdecomp
        for i in i_litr_min:i_litr_max
            vr = j + (i - 1) * nlevdecomp
            Cref[1, vr] = sum(getfield(cf, f)[1, j, i] for f in litr_c) * dt
            Nref[1, vr] = sum(getfield(nf, f)[1, j, i] for f in litr_n) * dt
        end
        vrc = j + (i_cwd - 1) * nlevdecomp
        Cref[1, vrc] = sum(getfield(cf, f)[1, j] for f in cwd_c) * dt
        Nref[1, vrc] = sum(getfield(nf, f)[1, j] for f in cwd_n) * dt
    end
    @test Cin ≈ Cref atol=1e-12
    @test Nin ≈ Nref atol=1e-12
    # CWD pool got only CWD fluxes (not litter); litter pools got only litter fluxes.
    @test Cin[1, 1 + (i_cwd-1)*nlevdecomp] > 0
    @test all(Cin[1, j + (i_cwd-1)*nlevdecomp] == sum(getfield(cf,f)[1,j] for f in cwd_c)*dt for j in 1:nlevdecomp)

    # transient adds harvest + gru; with those fields zero it's identical.
    for f in (:harvest_c_to_litr_c_col, :gru_c_to_litr_c_col, :harvest_c_to_cwdc_col, :gru_c_to_cwdc_col)
        fill!(getfield(cf, f), 0.0)
    end
    for f in (:harvest_n_to_litr_n_col, :gru_n_to_litr_n_col, :harvest_n_to_cwdn_col, :gru_n_to_cwdn_col)
        fill!(getfield(nf, f), 0.0)
    end
    Cin2 = zeros(nc, ndecomp_pools_vr); Nin2 = zeros(nc, ndecomp_pools_vr)
    CLM.cn_soil_matrix_input_accumulate!(Cin2, Nin2, cf, nf;
        mask_soilc=[true], bounds_col=1:nc, nlevdecomp=nlevdecomp,
        i_litr_min=i_litr_min, i_litr_max=i_litr_max, i_cwd=i_cwd, dt=dt,
        transient_landcover=true)
    @test Cin2 ≈ Cin atol=1e-12
    @test Nin2 ≈ Nin atol=1e-12
end
