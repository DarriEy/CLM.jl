@testset "Carbon Isotopes (C13/C14)" begin

    # -----------------------------------------------------------------------
    # Utility functions: delta13C <-> ratio conversions
    # -----------------------------------------------------------------------
    @testset "delta13C_to_ratio and ratio_to_delta13C" begin
        # PDB standard: delta=0 corresponds to R=0.0112372
        r0 = CLM.delta13C_to_ratio(0.0)
        @test r0 ≈ CLM.C13_PDB_RATIO / (1.0 + CLM.C13_PDB_RATIO) atol=1e-10

        # Round-trip: delta -> ratio -> delta
        for delta in [-25.0, -10.0, 0.0, 5.0, 50.0]
            r = CLM.delta13C_to_ratio(delta)
            delta_back = CLM.ratio_to_delta13C(r)
            @test delta_back ≈ delta atol=1e-8
        end

        # Edge cases
        @test CLM.ratio_to_delta13C(0.0) == 0.0
        @test CLM.ratio_to_delta13C(1.0) == 0.0
    end

    # -----------------------------------------------------------------------
    # c14_bomb_factor: latitude sector dispatch
    # -----------------------------------------------------------------------
    @testset "c14_bomb_factor latitude sectors" begin
        rc14 = [1.5, 1.2, 0.9]
        @test CLM.c14_bomb_factor(45.0; rc14_atm=rc14) == 1.5   # sector 1: >= 30
        @test CLM.c14_bomb_factor(30.0; rc14_atm=rc14) == 1.5   # boundary: >= 30
        @test CLM.c14_bomb_factor(0.0; rc14_atm=rc14)  == 1.2   # sector 2: -30 to 30
        @test CLM.c14_bomb_factor(-30.0; rc14_atm=rc14) == 1.2  # boundary: >= -30
        @test CLM.c14_bomb_factor(-45.0; rc14_atm=rc14) == 0.9  # sector 3: < -30
    end

    # -----------------------------------------------------------------------
    # C14 decay constants
    # -----------------------------------------------------------------------
    @testset "C14 decay constant" begin
        days_per_year = 365.0
        half_life = CLM.C14_HALF_LIFE_YEARS * CLM.SECSPDAY * days_per_year
        decay_const = -log(0.5) / half_life

        # The Fortran uses first-order decay: pool *= (1 - k * dt).
        # For small dt (one timestep), this closely approximates exp(-k*dt).
        # Verify for a typical 30-min timestep:
        dt = 1800.0
        factor = 1.0 - decay_const * dt
        factor_exact = exp(-decay_const * dt)
        @test factor ≈ factor_exact atol=1e-15  # nearly identical for small dt

        # After one year, pool should lose ~0.0121% (ln2/5730)
        dt_year = CLM.SECSPDAY * days_per_year
        factor_year_exact = exp(-log(2.0) / CLM.C14_HALF_LIFE_YEARS)
        # First-order approximation is close for one year too
        factor_year_approx = 1.0 - decay_const * dt_year
        @test factor_year_approx ≈ factor_year_exact atol=1e-8
        @test 0.999 < factor_year_approx < 1.0  # very small decay per year
    end

    # -----------------------------------------------------------------------
    # c14_decay! — radioactive decay of C14 pools
    # -----------------------------------------------------------------------
    @testset "c14_decay! basic" begin
        nc = 2; np = 3; ng = 1
        nlevdecomp = 2; ndecomp_pools = 3

        # Create C14 veg carbon state
        c14_cs = CLM.CNVegCarbonStateData()
        CLM.cnveg_carbon_state_init!(c14_cs, np, nc, ng)
        # Set some initial pool values
        c14_cs.leafc_patch .= [100.0, 200.0, 300.0]
        c14_cs.cpool_patch .= [50.0, 60.0, 70.0]
        c14_cs.xsmrpool_patch .= [10.0, 20.0, 30.0]
        c14_cs.deadstemc_patch .= [1000.0, 2000.0, 3000.0]
        c14_cs.gresp_storage_patch .= [5.0, 10.0, 15.0]
        c14_cs.gresp_xfer_patch .= [1.0, 2.0, 3.0]
        c14_cs.ctrunc_patch .= [0.1, 0.2, 0.3]
        c14_cs.seedc_grc .= [500.0]

        # Create C14 veg carbon flux (needed for matrix mode, but we use non-matrix)
        c14_cf = CLM.CNVegCarbonFluxData()
        CLM.cnveg_carbon_flux_init!(c14_cf, np, nc, ng)

        # Create soil C14 state and flux
        c14_soilcs = CLM.SoilBiogeochemCarbonStateData()
        CLM.soil_bgc_carbon_state_init!(c14_soilcs, nc, ng, nlevdecomp, ndecomp_pools)
        c14_soilcs.decomp_cpools_vr_col .= 100.0

        c14_soilcf = CLM.SoilBiogeochemCarbonFluxData()
        CLM.soil_bgc_carbon_flux_init!(c14_soilcf, nc, nlevdecomp, ndecomp_pools, 3)

        mask_c = trues(nc)
        mask_p = trues(np)
        dt = 1800.0  # 30 min timestep

        # Compute expected decay factor
        half_life = CLM.C14_HALF_LIFE_YEARS * CLM.SECSPDAY * 365.0
        decay_const = -log(0.5) / half_life
        decay_factor = 1.0 - decay_const * dt

        leafc_before = copy(c14_cs.leafc_patch)
        seed_before = copy(c14_cs.seedc_grc)
        soil_before = copy(c14_soilcs.decomp_cpools_vr_col)

        CLM.c14_decay!(c14_cs, c14_cf, c14_soilcs, c14_soilcf;
            mask_soilc=mask_c,
            mask_soilp=mask_p,
            bounds_col=1:nc,
            bounds_patch=1:np,
            bounds_gridcell=1:ng,
            dt=dt,
            nlevdecomp=nlevdecomp,
            ndecomp_pools=ndecomp_pools)

        # Check patch-level pools decayed
        for p in 1:np
            @test c14_cs.leafc_patch[p] ≈ leafc_before[p] * decay_factor atol=1e-10
            @test c14_cs.cpool_patch[p] ≈ (p == 1 ? 50.0 : p == 2 ? 60.0 : 70.0) * decay_factor atol=1e-10
        end

        # Check gridcell-level seedc decayed
        @test c14_cs.seedc_grc[1] ≈ seed_before[1] * decay_factor atol=1e-10

        # Check soil pools decayed
        for c in 1:nc, j in 1:nlevdecomp, l in 1:ndecomp_pools
            @test c14_soilcs.decomp_cpools_vr_col[c, j, l] ≈ soil_before[c, j, l] * decay_factor atol=1e-10
        end

        # Check that multiple timesteps accumulate correctly
        soil_before2 = copy(c14_soilcs.decomp_cpools_vr_col)
        CLM.c14_decay!(c14_cs, c14_cf, c14_soilcs, c14_soilcf;
            mask_soilc=mask_c,
            mask_soilp=mask_p,
            bounds_col=1:nc,
            bounds_patch=1:np,
            bounds_gridcell=1:ng,
            dt=dt,
            nlevdecomp=nlevdecomp,
            ndecomp_pools=ndecomp_pools)
        for c in 1:nc, j in 1:nlevdecomp, l in 1:ndecomp_pools
            @test c14_soilcs.decomp_cpools_vr_col[c, j, l] ≈ soil_before2[c, j, l] * decay_factor atol=1e-10
        end
    end

    # -----------------------------------------------------------------------
    # c13_c14_photosynthesis! — isotope ratio and flux computation
    # -----------------------------------------------------------------------
    @testset "c13_c14_photosynthesis!" begin
        np = 3; ng = 1

        ps = CLM.PhotosynthesisData()
        CLM.photosynthesis_data_init!(ps, np; nlevcan=1)

        # Set up photosynthesis state
        ps.psnsun_patch .= [5.0, 10.0, 0.0]
        ps.psnsha_patch .= [3.0, 6.0, 0.0]
        ps.alphapsnsun_patch .= [1.02, 1.03, 1.0]
        ps.alphapsnsha_patch .= [1.01, 1.02, 1.0]

        forc_pco2 = [400.0e-6]     # 400 ppm in Pa (approx)
        forc_pc13o2 = [4.4e-6]     # ~1.1% of CO2
        gridcell_of_patch = [1, 1, 1]
        latdeg = [51.0]
        mask = trues(np)

        CLM.c13_c14_photosynthesis!(ps,
            forc_pco2, forc_pc13o2,
            gridcell_of_patch, latdeg,
            mask, 1:np;
            use_c13=true, use_c14=true)

        # Check C13 canopy air ratio
        expected_rc13 = forc_pc13o2[1] / (forc_pco2[1] - forc_pc13o2[1])
        @test ps.rc13_canair_patch[1] ≈ expected_rc13 atol=1e-10

        # Check C13 fractionation
        @test ps.rc13_psnsun_patch[1] ≈ expected_rc13 / 1.02 atol=1e-10
        @test ps.rc13_psnsha_patch[1] ≈ expected_rc13 / 1.01 atol=1e-10

        # Check C13 psn flux
        rc_sun = ps.rc13_psnsun_patch[1]
        @test ps.c13_psnsun_patch[1] ≈ 5.0 * rc_sun / (1.0 + rc_sun) atol=1e-10

        # Check C14 psn flux (lat >= 30 -> sector 1, default ratio 1.0)
        @test ps.c14_psnsun_patch[1] ≈ 5.0 * 1.0 atol=1e-10
        @test ps.c14_psnsha_patch[1] ≈ 3.0 * 1.0 atol=1e-10

        # Zero photosynthesis patch should produce zero isotope fluxes
        @test ps.c13_psnsun_patch[3] == 0.0
        @test ps.c14_psnsun_patch[3] == 0.0
    end

    # -----------------------------------------------------------------------
    # C13/C14 photosynthesis with custom bomb spike ratios
    # -----------------------------------------------------------------------
    @testset "c14_bomb_spike_custom" begin
        np = 3; ng = 3

        ps = CLM.PhotosynthesisData()
        CLM.photosynthesis_data_init!(ps, np; nlevcan=1)

        ps.psnsun_patch .= [10.0, 10.0, 10.0]
        ps.psnsha_patch .= [5.0, 5.0, 5.0]
        ps.alphapsnsun_patch .= [1.0, 1.0, 1.0]
        ps.alphapsnsha_patch .= [1.0, 1.0, 1.0]

        forc_pco2 = [400e-6, 400e-6, 400e-6]
        forc_pc13o2 = [4.4e-6, 4.4e-6, 4.4e-6]
        gridcell_of_patch = [1, 2, 3]
        latdeg = [45.0, 0.0, -45.0]  # one per sector
        mask = trues(np)

        rc14_custom = [1.5, 1.2, 0.8]

        CLM.c13_c14_photosynthesis!(ps,
            forc_pco2, forc_pc13o2,
            gridcell_of_patch, latdeg,
            mask, 1:np;
            use_c13=false, use_c14=true,
            rc14_atm=rc14_custom)

        @test ps.c14_psnsun_patch[1] ≈ 10.0 * 1.5 atol=1e-10  # NH
        @test ps.c14_psnsun_patch[2] ≈ 10.0 * 1.2 atol=1e-10  # tropics
        @test ps.c14_psnsun_patch[3] ≈ 10.0 * 0.8 atol=1e-10  # SH
    end

end
