@testset "Nitrif-Denitrif" begin

    # ----------------------------------------------------------------
    # Helper: create minimal data structures for nc columns
    # ----------------------------------------------------------------
    function make_nitrif_denitrif_data(; nc=4, nlevdecomp=1)
        ndecomp_pools = 7
        ndecomp_cascade_transitions = 10

        # --- Nitrif-Denitrif Parameters ---
        params = CLM.NitrifDenitrifParams(
            k_nitr_max_perday                = 0.1,
            surface_tension_water            = 0.073,
            rij_kro_a                        = 1.5e-10,
            rij_kro_alpha                    = 1.26,
            rij_kro_beta                     = 0.6,
            rij_kro_gamma                    = 0.6,
            rij_kro_delta                    = 0.85,
            denitrif_respiration_coefficient = 0.1,
            denitrif_respiration_exponent    = 1.3,
            denitrif_nitrateconc_coefficient = 0.1,
            denitrif_nitrateconc_exponent    = 1.3,
            om_frac_sf                       = 1.0,
        )

        # --- CN shared params ---
        cn_params = CLM.CNSharedParamsData(
            organic_max = 130.0,
        )

        # --- Nitrogen flux data ---
        nf = CLM.SoilBiogeochemNitrogenFluxData()
        CLM.soil_bgc_nitrogen_flux_init!(nf, nc, nlevdecomp, ndecomp_pools,
                                          ndecomp_cascade_transitions)

        # --- Nitrogen state data ---
        ns = CLM.SoilBiogeochemNitrogenStateData()
        CLM.soil_bgc_nitrogen_state_init!(ns, nc, 1, nlevdecomp, ndecomp_pools)
        # Set mineral N pools
        for j in 1:nlevdecomp
            for c in 1:nc
                ns.smin_nh4_vr_col[c, j] = 1.0   # gN/m3
                ns.smin_no3_vr_col[c, j] = 2.0   # gN/m3
            end
        end

        # --- Carbon flux data ---
        cf = CLM.SoilBiogeochemCarbonFluxData()
        CLM.soil_bgc_carbon_flux_init!(cf, nc, nlevdecomp, ndecomp_pools,
                                        ndecomp_cascade_transitions)
        # Set scalars and potential HR
        for j in 1:nlevdecomp
            for c in 1:nc
                cf.t_scalar_col[c, j] = 0.8
                cf.w_scalar_col[c, j] = 0.6
                cf.phr_vr_col[c, j]   = 1.0e-6  # gC/m3/s (moderate respiration)
            end
        end

        # --- Soil state arrays ---
        watsat   = fill(0.45, nc, nlevdecomp)   # porosity
        watfc    = fill(0.30, nc, nlevdecomp)   # field capacity
        bd       = fill(1300.0, nc, nlevdecomp) # bulk density (kg/m3)
        bsw      = fill(5.0, nc, nlevdecomp)    # Clapp-Hornberger b
        cellorg  = fill(20.0, nc, nlevdecomp)   # organic matter (kg/m3)
        sucsat   = fill(200.0, nc, nlevdecomp)  # minimum suction (mm)
        soilpsi  = fill(-1.0, nc, nlevdecomp)   # soil water potential (MPa)

        # --- Water state arrays ---
        h2osoi_vol = fill(0.30, nc, nlevdecomp) # volumetric soil water
        h2osoi_liq = fill(30.0, nc, nlevdecomp) # liquid water (kg/m2)

        # --- Temperature ---
        t_soisno = fill(CLM.TFRZ + 15.0, nc, nlevdecomp)

        # --- Column geometry ---
        col_dz = fill(0.1, nc, nlevdecomp)

        # --- CH4 arrays ---
        o2_decomp_depth_unsat = fill(1.0e-6, nc, nlevdecomp)  # mol/m3/s
        conc_o2_unsat         = fill(0.2, nc, nlevdecomp)      # mol/m3

        # --- Mask ---
        mask_bgc_soilc = trues(nc)

        return (params=params, cn_params=cn_params, nf=nf, ns=ns, cf=cf,
                watsat=watsat, watfc=watfc, bd=bd, bsw=bsw, cellorg=cellorg,
                sucsat=sucsat, soilpsi=soilpsi, h2osoi_vol=h2osoi_vol,
                h2osoi_liq=h2osoi_liq, t_soisno=t_soisno, col_dz=col_dz,
                o2_decomp_depth_unsat=o2_decomp_depth_unsat,
                conc_o2_unsat=conc_o2_unsat,
                mask_bgc_soilc=mask_bgc_soilc,
                nlevdecomp=nlevdecomp, nc=nc)
    end

    # ================================================================
    # Test 1: NitrifDenitrifParams construction
    # ================================================================
    @testset "NitrifDenitrifParams construction" begin
        params = CLM.NitrifDenitrifParams()
        @test params.k_nitr_max_perday == 0.1
        @test params.surface_tension_water == 0.073
        @test params.rij_kro_a == 1.5e-10
        @test params.om_frac_sf == 1.0
    end

    # ================================================================
    # Test 2: nitrif_denitrif_read_params!
    # ================================================================
    @testset "nitrif_denitrif_read_params!" begin
        params = CLM.NitrifDenitrifParams()
        CLM.nitrif_denitrif_read_params!(params;
            k_nitr_max_perday                = 0.2,
            surface_tension_water            = 0.08,
            rij_kro_a                        = 2.0e-10,
            rij_kro_alpha                    = 1.3,
            rij_kro_beta                     = 0.7,
            rij_kro_gamma                    = 0.7,
            rij_kro_delta                    = 0.9,
            denitrif_respiration_coefficient = 0.2,
            denitrif_respiration_exponent    = 1.4,
            denitrif_nitrateconc_coefficient = 0.2,
            denitrif_nitrateconc_exponent    = 1.4,
            om_frac_sf                       = 0.8)
        @test params.k_nitr_max_perday == 0.2
        @test params.surface_tension_water == 0.08
        @test params.rij_kro_a == 2.0e-10
        @test params.denitrif_respiration_coefficient == 0.2
        @test params.om_frac_sf == 0.8
    end

    # ================================================================
    # Test 3: nitrif_denitrif! with use_lch4=false
    # ================================================================
    @testset "nitrif_denitrif! (use_lch4=false)" begin
        d = make_nitrif_denitrif_data()

        CLM.nitrif_denitrif!(d.nf, d.ns, d.cf, d.params, d.cn_params;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:d.nc,
            nlevdecomp=d.nlevdecomp,
            watsat=d.watsat, watfc=d.watfc, bd=d.bd, bsw=d.bsw,
            cellorg=d.cellorg, sucsat=d.sucsat, soilpsi=d.soilpsi,
            h2osoi_vol=d.h2osoi_vol, h2osoi_liq=d.h2osoi_liq,
            t_soisno=d.t_soisno, col_dz=d.col_dz,
            use_lch4=false,
            no_frozen_nitrif_denitrif=false)

        for c in 1:d.nc
            # Anaerobic fraction should be zero without CH4 model
            @test d.nf.anaerobic_frac_col[c, 1] == 0.0
            @test d.nf.diffus_col[c, 1] == 0.0

            # Nitrification should be positive (warm, moist, NH4 available)
            @test d.nf.pot_f_nit_vr_col[c, 1] > 0.0

            # Temperature scalar should be positive and <= 1
            @test d.nf.k_nitr_t_vr_col[c, 1] > 0.0
            @test d.nf.k_nitr_t_vr_col[c, 1] <= 1.0

            # Moisture scalar should be positive
            @test d.nf.k_nitr_h2o_vr_col[c, 1] > 0.0

            # pH scalar should be in reasonable range
            @test d.nf.k_nitr_ph_vr_col[c, 1] > 0.0
            @test d.nf.k_nitr_ph_vr_col[c, 1] < 1.0

            # Denitrification potential should be non-negative
            @test d.nf.f_denit_base_vr_col[c, 1] >= 0.0

            # Potential denitrification should be zero (anaerobic_frac=0)
            @test d.nf.pot_f_denit_vr_col[c, 1] == 0.0

            # N2:N2O ratio should be positive
            @test d.nf.n2_n2o_ratio_denit_vr_col[c, 1] > 0.0

            # WFPS should be reasonable
            @test d.nf.wfps_vr_col[c, 1] > 0.0
            @test d.nf.wfps_vr_col[c, 1] <= 100.0

            # Soil bulk density should be positive
            @test d.nf.soil_bulkdensity_col[c, 1] > 0.0

            # NO3 mass density should be non-negative
            @test d.nf.smin_no3_massdens_vr_col[c, 1] >= 0.0
        end
    end

    # ================================================================
    # Test 4: nitrif_denitrif! with use_lch4=true
    # ================================================================
    @testset "nitrif_denitrif! (use_lch4=true)" begin
        d = make_nitrif_denitrif_data()

        CLM.nitrif_denitrif!(d.nf, d.ns, d.cf, d.params, d.cn_params;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:d.nc,
            nlevdecomp=d.nlevdecomp,
            watsat=d.watsat, watfc=d.watfc, bd=d.bd, bsw=d.bsw,
            cellorg=d.cellorg, sucsat=d.sucsat, soilpsi=d.soilpsi,
            h2osoi_vol=d.h2osoi_vol, h2osoi_liq=d.h2osoi_liq,
            t_soisno=d.t_soisno, col_dz=d.col_dz,
            o2_decomp_depth_unsat=d.o2_decomp_depth_unsat,
            conc_o2_unsat=d.conc_o2_unsat,
            use_lch4=true,
            no_frozen_nitrif_denitrif=false)

        for c in 1:d.nc
            # Anaerobic fraction should be in [0, 1]
            @test d.nf.anaerobic_frac_col[c, 1] >= 0.0
            @test d.nf.anaerobic_frac_col[c, 1] <= 1.0

            # r_psi should be positive
            @test d.nf.r_psi_col[c, 1] > 0.0

            # Diffusivity should be positive (now in m2/s)
            @test d.nf.diffus_col[c, 1] >= 0.0

            # Nitrification should be positive
            @test d.nf.pot_f_nit_vr_col[c, 1] > 0.0

            # Denitrification base rate should be non-negative
            @test d.nf.f_denit_base_vr_col[c, 1] >= 0.0

            # N2:N2O ratio should be positive
            @test d.nf.n2_n2o_ratio_denit_vr_col[c, 1] > 0.0

            # ratio_k1 should be >= 1.7 (minimum from the max expression)
            @test d.nf.ratio_k1_col[c, 1] >= 1.7
        end
    end

    # ================================================================
    # Test 5: Frozen soil behavior
    # ================================================================
    @testset "Frozen soil suppression" begin
        d = make_nitrif_denitrif_data()

        # Set temperature below freezing
        t_frozen = fill(CLM.TFRZ - 5.0, d.nc, d.nlevdecomp)

        CLM.nitrif_denitrif!(d.nf, d.ns, d.cf, d.params, d.cn_params;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:d.nc,
            nlevdecomp=d.nlevdecomp,
            watsat=d.watsat, watfc=d.watfc, bd=d.bd, bsw=d.bsw,
            cellorg=d.cellorg, sucsat=d.sucsat, soilpsi=d.soilpsi,
            h2osoi_vol=d.h2osoi_vol, h2osoi_liq=d.h2osoi_liq,
            t_soisno=t_frozen, col_dz=d.col_dz,
            use_lch4=false,
            no_frozen_nitrif_denitrif=true)

        for c in 1:d.nc
            # With no_frozen_nitrif_denitrif=true and T<TFRZ,
            # both nitrification and denitrification should be zero
            @test d.nf.pot_f_nit_vr_col[c, 1] == 0.0
            @test d.nf.f_denit_base_vr_col[c, 1] == 0.0
        end

        # Now test that frozen soil does NOT suppress when flag is false
        d2 = make_nitrif_denitrif_data()

        CLM.nitrif_denitrif!(d2.nf, d2.ns, d2.cf, d2.params, d2.cn_params;
            mask_bgc_soilc=d2.mask_bgc_soilc,
            bounds=1:d2.nc,
            nlevdecomp=d2.nlevdecomp,
            watsat=d2.watsat, watfc=d2.watfc, bd=d2.bd, bsw=d2.bsw,
            cellorg=d2.cellorg, sucsat=d2.sucsat, soilpsi=d2.soilpsi,
            h2osoi_vol=d2.h2osoi_vol, h2osoi_liq=d2.h2osoi_liq,
            t_soisno=t_frozen, col_dz=d2.col_dz,
            use_lch4=false,
            no_frozen_nitrif_denitrif=false)

        for c in 1:d2.nc
            # With no_frozen_nitrif_denitrif=false, nitrification should still be positive
            # (T<TFRZ doesn't suppress if flag is off)
            @test d2.nf.pot_f_nit_vr_col[c, 1] > 0.0
        end
    end

    # ================================================================
    # Test 6: NH4 dependence
    # ================================================================
    @testset "NH4 dependence of nitrification" begin
        d_high = make_nitrif_denitrif_data()
        d_low  = make_nitrif_denitrif_data()

        # High NH4
        for j in 1:d_high.nlevdecomp
            for c in 1:d_high.nc
                d_high.ns.smin_nh4_vr_col[c, j] = 10.0
            end
        end

        # Low NH4
        for j in 1:d_low.nlevdecomp
            for c in 1:d_low.nc
                d_low.ns.smin_nh4_vr_col[c, j] = 0.1
            end
        end

        # Run both
        CLM.nitrif_denitrif!(d_high.nf, d_high.ns, d_high.cf, d_high.params, d_high.cn_params;
            mask_bgc_soilc=d_high.mask_bgc_soilc, bounds=1:d_high.nc,
            nlevdecomp=d_high.nlevdecomp,
            watsat=d_high.watsat, watfc=d_high.watfc, bd=d_high.bd, bsw=d_high.bsw,
            cellorg=d_high.cellorg, sucsat=d_high.sucsat, soilpsi=d_high.soilpsi,
            h2osoi_vol=d_high.h2osoi_vol, h2osoi_liq=d_high.h2osoi_liq,
            t_soisno=d_high.t_soisno, col_dz=d_high.col_dz,
            use_lch4=false)

        CLM.nitrif_denitrif!(d_low.nf, d_low.ns, d_low.cf, d_low.params, d_low.cn_params;
            mask_bgc_soilc=d_low.mask_bgc_soilc, bounds=1:d_low.nc,
            nlevdecomp=d_low.nlevdecomp,
            watsat=d_low.watsat, watfc=d_low.watfc, bd=d_low.bd, bsw=d_low.bsw,
            cellorg=d_low.cellorg, sucsat=d_low.sucsat, soilpsi=d_low.soilpsi,
            h2osoi_vol=d_low.h2osoi_vol, h2osoi_liq=d_low.h2osoi_liq,
            t_soisno=d_low.t_soisno, col_dz=d_low.col_dz,
            use_lch4=false)

        # Higher NH4 should give higher nitrification
        @test d_high.nf.pot_f_nit_vr_col[1, 1] > d_low.nf.pot_f_nit_vr_col[1, 1]
    end

    # ================================================================
    # Test 7: NO3 dependence of denitrification
    # ================================================================
    @testset "NO3 dependence of denitrification" begin
        d_high = make_nitrif_denitrif_data()
        d_low  = make_nitrif_denitrif_data()

        # High NO3
        for j in 1:d_high.nlevdecomp
            for c in 1:d_high.nc
                d_high.ns.smin_no3_vr_col[c, j] = 20.0
            end
        end

        # Low NO3
        for j in 1:d_low.nlevdecomp
            for c in 1:d_low.nc
                d_low.ns.smin_no3_vr_col[c, j] = 0.1
            end
        end

        # Run both
        CLM.nitrif_denitrif!(d_high.nf, d_high.ns, d_high.cf, d_high.params, d_high.cn_params;
            mask_bgc_soilc=d_high.mask_bgc_soilc, bounds=1:d_high.nc,
            nlevdecomp=d_high.nlevdecomp,
            watsat=d_high.watsat, watfc=d_high.watfc, bd=d_high.bd, bsw=d_high.bsw,
            cellorg=d_high.cellorg, sucsat=d_high.sucsat, soilpsi=d_high.soilpsi,
            h2osoi_vol=d_high.h2osoi_vol, h2osoi_liq=d_high.h2osoi_liq,
            t_soisno=d_high.t_soisno, col_dz=d_high.col_dz,
            use_lch4=false)

        CLM.nitrif_denitrif!(d_low.nf, d_low.ns, d_low.cf, d_low.params, d_low.cn_params;
            mask_bgc_soilc=d_low.mask_bgc_soilc, bounds=1:d_low.nc,
            nlevdecomp=d_low.nlevdecomp,
            watsat=d_low.watsat, watfc=d_low.watfc, bd=d_low.bd, bsw=d_low.bsw,
            cellorg=d_low.cellorg, sucsat=d_low.sucsat, soilpsi=d_low.soilpsi,
            h2osoi_vol=d_low.h2osoi_vol, h2osoi_liq=d_low.h2osoi_liq,
            t_soisno=d_low.t_soisno, col_dz=d_low.col_dz,
            use_lch4=false)

        # Higher NO3 should give higher denitrification base rate
        @test d_high.nf.f_denit_base_vr_col[1, 1] >= d_low.nf.f_denit_base_vr_col[1, 1]
    end

    # ================================================================
    # Test 8: Multi-level
    # ================================================================
    @testset "Multi-level computation" begin
        nlevdecomp = 3
        d = make_nitrif_denitrif_data(; nlevdecomp=nlevdecomp)

        CLM.nitrif_denitrif!(d.nf, d.ns, d.cf, d.params, d.cn_params;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:d.nc,
            nlevdecomp=nlevdecomp,
            watsat=d.watsat, watfc=d.watfc, bd=d.bd, bsw=d.bsw,
            cellorg=d.cellorg, sucsat=d.sucsat, soilpsi=d.soilpsi,
            h2osoi_vol=d.h2osoi_vol, h2osoi_liq=d.h2osoi_liq,
            t_soisno=d.t_soisno, col_dz=d.col_dz,
            use_lch4=false)

        # All levels should have valid outputs
        for j in 1:nlevdecomp
            for c in 1:d.nc
                @test d.nf.pot_f_nit_vr_col[c, j] > 0.0
                @test d.nf.f_denit_base_vr_col[c, j] >= 0.0
                @test d.nf.n2_n2o_ratio_denit_vr_col[c, j] > 0.0
                @test d.nf.soil_bulkdensity_col[c, j] > 0.0
            end
        end
    end

    # ================================================================
    # Test 9: WFPS calculation
    # ================================================================
    @testset "WFPS calculation" begin
        d = make_nitrif_denitrif_data()

        # Set h2osoi_vol = watsat (saturated)
        d.h2osoi_vol .= d.watsat

        CLM.nitrif_denitrif!(d.nf, d.ns, d.cf, d.params, d.cn_params;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:d.nc,
            nlevdecomp=d.nlevdecomp,
            watsat=d.watsat, watfc=d.watfc, bd=d.bd, bsw=d.bsw,
            cellorg=d.cellorg, sucsat=d.sucsat, soilpsi=d.soilpsi,
            h2osoi_vol=d.h2osoi_vol, h2osoi_liq=d.h2osoi_liq,
            t_soisno=d.t_soisno, col_dz=d.col_dz,
            use_lch4=false)

        # At saturation, WFPS should be 100%
        for c in 1:d.nc
            @test d.nf.wfps_vr_col[c, 1] ≈ 100.0
        end

        # fr_WFPS at saturation: max(0.1, 0.015*100 - 0.32) = max(0.1, 1.18) = 1.18
        for c in 1:d.nc
            @test d.nf.fr_WFPS_col[c, 1] ≈ 1.18
        end
    end

    # ================================================================
    # Test 10: Mask behavior
    # ================================================================
    @testset "Mask filtering" begin
        d = make_nitrif_denitrif_data()

        # Mask out columns 2 and 4
        d.mask_bgc_soilc[2] = false
        d.mask_bgc_soilc[4] = false

        # Initialize output fields to NaN (already NaN from init)
        # Run the function
        CLM.nitrif_denitrif!(d.nf, d.ns, d.cf, d.params, d.cn_params;
            mask_bgc_soilc=d.mask_bgc_soilc,
            bounds=1:d.nc,
            nlevdecomp=d.nlevdecomp,
            watsat=d.watsat, watfc=d.watfc, bd=d.bd, bsw=d.bsw,
            cellorg=d.cellorg, sucsat=d.sucsat, soilpsi=d.soilpsi,
            h2osoi_vol=d.h2osoi_vol, h2osoi_liq=d.h2osoi_liq,
            t_soisno=d.t_soisno, col_dz=d.col_dz,
            use_lch4=false)

        # Active columns should have computed values
        @test d.nf.pot_f_nit_vr_col[1, 1] > 0.0
        @test d.nf.pot_f_nit_vr_col[3, 1] > 0.0

        # Masked columns should still have NaN (not updated)
        @test isnan(d.nf.pot_f_nit_vr_col[2, 1])
        @test isnan(d.nf.pot_f_nit_vr_col[4, 1])
    end

end
