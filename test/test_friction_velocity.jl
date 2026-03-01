@testset "FrictionVelocityData" begin

    @testset "default construction" begin
        fvd = CLM.FrictionVelocityData()
        @test fvd.zetamaxstable == -999.0
        @test fvd.zsno == -999.0
        @test fvd.zlnd == -999.0
        @test fvd.zglc == -999.0
        @test length(fvd.forc_hgt_u_patch) == 0
        @test length(fvd.z0mg_col) == 0
        @test length(fvd.z0m_actual_patch) == 0
    end

    @testset "frictionvel_init!" begin
        np = 8
        nc = 4
        fvd = CLM.FrictionVelocityData()
        CLM.frictionvel_init!(fvd, np, nc)

        # Check patch-level sizes
        @test length(fvd.forc_hgt_u_patch) == np
        @test length(fvd.forc_hgt_t_patch) == np
        @test length(fvd.forc_hgt_q_patch) == np
        @test length(fvd.u10_patch) == np
        @test length(fvd.u10_clm_patch) == np
        @test length(fvd.va_patch) == np
        @test length(fvd.vds_patch) == np
        @test length(fvd.fv_patch) == np
        @test length(fvd.rb1_patch) == np
        @test length(fvd.rb10_patch) == np
        @test length(fvd.ram1_patch) == np
        @test length(fvd.z0mv_patch) == np
        @test length(fvd.z0hv_patch) == np
        @test length(fvd.z0qv_patch) == np
        @test length(fvd.z0mg_patch) == np
        @test length(fvd.z0hg_patch) == np
        @test length(fvd.z0qg_patch) == np
        @test length(fvd.kbm1_patch) == np
        @test length(fvd.rah1_patch) == np
        @test length(fvd.rah2_patch) == np
        @test length(fvd.raw1_patch) == np
        @test length(fvd.raw2_patch) == np
        @test length(fvd.ustar_patch) == np
        @test length(fvd.um_patch) == np
        @test length(fvd.uaf_patch) == np
        @test length(fvd.taf_patch) == np
        @test length(fvd.qaf_patch) == np
        @test length(fvd.obu_patch) == np
        @test length(fvd.zeta_patch) == np
        @test length(fvd.vpd_patch) == np
        @test length(fvd.num_iter_patch) == np
        @test length(fvd.z0m_actual_patch) == np

        # Check column-level sizes
        @test length(fvd.z0mg_col) == nc
        @test length(fvd.z0mg_2D_col) == nc
        @test length(fvd.z0hg_col) == nc
        @test length(fvd.z0qg_col) == nc

        # Check NaN initialization
        @test all(isnan, fvd.forc_hgt_u_patch)
        @test all(isnan, fvd.z0mg_col)
        @test all(isnan, fvd.ustar_patch)

        # rb10_patch initialized to SPVAL
        @test all(x -> x == CLM.SPVAL, fvd.rb10_patch)
    end

    @testset "frictionvel_clean!" begin
        fvd = CLM.FrictionVelocityData()
        CLM.frictionvel_init!(fvd, 5, 3)
        CLM.frictionvel_clean!(fvd)

        @test length(fvd.forc_hgt_u_patch) == 0
        @test length(fvd.z0mg_col) == 0
        @test length(fvd.ustar_patch) == 0
        @test length(fvd.z0m_actual_patch) == 0
        @test length(fvd.z0mg_2D_col) == 0
    end

    @testset "frictionvel_init_for_testing!" begin
        fvd = CLM.FrictionVelocityData()
        CLM.frictionvel_init_for_testing!(fvd, 4, 2)

        @test fvd.zetamaxstable ≈ 0.5
        @test fvd.zsno ≈ 0.00085
        @test fvd.zlnd ≈ 0.000775
        @test fvd.zglc ≈ 0.00230000005
        @test length(fvd.forc_hgt_u_patch) == 4
        @test length(fvd.z0mg_col) == 2
    end

    @testset "frictionvel_init_cold!" begin
        fvd = CLM.FrictionVelocityData()
        CLM.frictionvel_init!(fvd, 4, 3)

        # With use_cn=true, forc_hgt_u_patch should be set to 30.0
        CLM.frictionvel_init_cold!(fvd, 1:3, 1:4; use_cn=true)
        for p in 1:4
            @test fvd.forc_hgt_u_patch[p] ≈ 30.0
        end

        # Lake column initialization
        col_landunit = [1, 2, 1]
        lun_lakpoi = [false, true]
        CLM.frictionvel_init_cold!(fvd, 1:3, 1:4;
            col_landunit=col_landunit,
            lun_lakpoi=lun_lakpoi,
            use_cn=false)

        @test isnan(fvd.z0mg_col[1])  # non-lake
        @test fvd.z0mg_col[2] ≈ 0.0004  # lake
        @test isnan(fvd.z0mg_col[3])  # non-lake
    end

    @testset "stability_func1" begin
        # At zeta = 0, should return 0
        @test CLM.stability_func1(0.0) ≈ 0.0 atol=1e-12

        # For negative zeta (unstable), should return positive value
        val = CLM.stability_func1(-0.5)
        @test val > 0.0
        @test isfinite(val)

        # Known reference: at zeta = -1
        val1 = CLM.stability_func1(-1.0)
        @test isfinite(val1)
        @test val1 > 0.0
    end

    @testset "stability_func2" begin
        # At zeta = 0, should return 0
        @test CLM.stability_func2(0.0) ≈ 0.0 atol=1e-12

        # For negative zeta (unstable), should return positive value
        val = CLM.stability_func2(-0.5)
        @test val > 0.0
        @test isfinite(val)
    end

    @testset "monin_obuk_ini" begin
        zetamaxstable = 0.5

        # Stable case (dthv >= 0)
        ur = 5.0; thv = 300.0; dthv = 1.0; zldis = 30.0; z0m = 0.01
        um_out, obu_out = CLM.monin_obuk_ini(zetamaxstable, ur, thv, dthv, zldis, z0m)
        @test um_out ≈ max(ur, 0.1)
        @test isfinite(obu_out)
        @test obu_out > 0.0  # stable → positive obu

        # Unstable case (dthv < 0)
        dthv_neg = -1.0
        um_out2, obu_out2 = CLM.monin_obuk_ini(zetamaxstable, ur, thv, dthv_neg, zldis, z0m)
        @test um_out2 ≈ sqrt(ur^2 + 0.5^2)
        @test isfinite(obu_out2)
        @test obu_out2 < 0.0  # unstable → negative obu

        # Near-zero wind speed
        ur_low = 0.01
        um_out3, obu_out3 = CLM.monin_obuk_ini(zetamaxstable, ur_low, thv, 0.5, zldis, z0m)
        @test um_out3 ≈ 0.1  # clamped to 0.1
        @test isfinite(obu_out3)
    end

    @testset "friction_velocity! (patch-level, unstable)" begin
        np = 2
        nc = 2
        fvd = CLM.FrictionVelocityData()
        CLM.frictionvel_init_for_testing!(fvd, np, nc)

        # Set forcing heights
        fvd.forc_hgt_u_patch .= 30.0
        fvd.forc_hgt_t_patch .= 30.0
        fvd.forc_hgt_q_patch .= 30.0

        fn = 1
        filtern = [1]
        displa_v = [1.0, 1.0]
        z0m_v = [0.1, 0.1]
        z0h_v = [0.1, 0.1]
        z0q_v = [0.1, 0.1]
        obu_v = [-50.0, -50.0]  # unstable
        ur_v = [5.0, 5.0]
        um_v = [5.0, 5.0]
        ustar_v = [0.0, 0.0]
        temp1_v = [0.0, 0.0]
        temp2_v = [0.0, 0.0]
        temp12m_v = [0.0, 0.0]
        temp22m_v = [0.0, 0.0]
        fm_v = [0.0, 0.0]
        iter = 1

        CLM.friction_velocity!(fvd, fn, filtern, displa_v, z0m_v, z0h_v, z0q_v,
            obu_v, iter, ur_v, um_v, ustar_v, temp1_v, temp2_v,
            temp12m_v, temp22m_v, fm_v)

        # ustar should be positive and finite
        @test ustar_v[1] > 0.0
        @test isfinite(ustar_v[1])

        # temp1, temp2 should be positive (unstable)
        @test temp1_v[1] > 0.0
        @test temp2_v[1] > 0.0

        # 10m wind and friction velocity should be set
        @test isfinite(fvd.u10_patch[1])
        @test fvd.fv_patch[1] ≈ ustar_v[1]

        # va should be set to um
        @test fvd.va_patch[1] ≈ um_v[1]

        # vds should be positive
        @test fvd.vds_patch[1] > 0.0
    end

    @testset "friction_velocity! (patch-level, stable)" begin
        np = 2
        nc = 2
        fvd = CLM.FrictionVelocityData()
        CLM.frictionvel_init_for_testing!(fvd, np, nc)

        fvd.forc_hgt_u_patch .= 30.0
        fvd.forc_hgt_t_patch .= 30.0
        fvd.forc_hgt_q_patch .= 30.0

        fn = 1
        filtern = [1]
        displa_v = [1.0, 1.0]
        z0m_v = [0.1, 0.1]
        z0h_v = [0.1, 0.1]
        z0q_v = [0.1, 0.1]
        obu_v = [50.0, 50.0]  # stable
        ur_v = [5.0, 5.0]
        um_v = [5.0, 5.0]
        ustar_v = [0.0, 0.0]
        temp1_v = [0.0, 0.0]
        temp2_v = [0.0, 0.0]
        temp12m_v = [0.0, 0.0]
        temp22m_v = [0.0, 0.0]
        fm_v = [0.0, 0.0]
        iter = 1

        CLM.friction_velocity!(fvd, fn, filtern, displa_v, z0m_v, z0h_v, z0q_v,
            obu_v, iter, ur_v, um_v, ustar_v, temp1_v, temp2_v,
            temp12m_v, temp22m_v, fm_v)

        @test ustar_v[1] > 0.0
        @test isfinite(ustar_v[1])
        @test isfinite(temp1_v[1])
        @test isfinite(temp2_v[1])
        @test isfinite(temp12m_v[1])
        @test isfinite(temp22m_v[1])
        @test isfinite(fvd.u10_patch[1])
    end

    @testset "friction_velocity! iter averaging of fm" begin
        np = 1
        nc = 1
        fvd = CLM.FrictionVelocityData()
        CLM.frictionvel_init_for_testing!(fvd, np, nc)
        fvd.forc_hgt_u_patch .= 30.0
        fvd.forc_hgt_t_patch .= 30.0
        fvd.forc_hgt_q_patch .= 30.0

        filtern = [1]
        displa_v = [1.0]
        z0m_v = [0.1]
        z0h_v = [0.1]
        z0q_v = [0.1]
        obu_v = [-50.0]
        ur_v = [5.0]
        um_v = [5.0]
        ustar_v = [0.0]
        temp1_v = [0.0]
        temp2_v = [0.0]
        temp12m_v = [0.0]
        temp22m_v = [0.0]
        fm_v = [0.0]

        # iter 1: fm = fmnew
        CLM.friction_velocity!(fvd, 1, filtern, displa_v, z0m_v, z0h_v, z0q_v,
            obu_v, 1, ur_v, um_v, ustar_v, temp1_v, temp2_v,
            temp12m_v, temp22m_v, fm_v)
        fm_iter1 = fm_v[1]

        # iter 2: fm = 0.5*(fm_old + fmnew) = 0.5*(fm_iter1 + fmnew)
        # Since conditions haven't changed, fmnew is the same, so fm = fm_iter1
        CLM.friction_velocity!(fvd, 1, filtern, displa_v, z0m_v, z0h_v, z0q_v,
            obu_v, 2, ur_v, um_v, ustar_v, temp1_v, temp2_v,
            temp12m_v, temp22m_v, fm_v)
        @test fm_v[1] ≈ fm_iter1  # same conditions → same fmnew → average is unchanged
    end

    @testset "set_roughness_and_forc_heights_nonlake! ZengWang2007" begin
        np = 2
        nc = 2
        fvd = CLM.FrictionVelocityData()
        CLM.frictionvel_init_for_testing!(fvd, np, nc)

        mask_nolakec = BitVector([true, true])
        mask_nolakep = BitVector([true, true])
        frac_sno = [0.5, 0.0]  # first column has snow
        snomelt_accum = [0.0, 0.0]
        frac_veg_nosno = [1, 0]  # first patch has veg, second bare
        z0m_canopy = [0.5, 0.1]
        displa_v = [2.0, 0.0]
        forc_hgt_u = [10.0]
        forc_hgt_t = [10.0]
        forc_hgt_q = [10.0]
        col_landunit = [1, 1]
        patch_gridcell = [1, 1]
        patch_landunit = [1, 1]
        patch_column = [1, 2]
        lun_itype = [CLM.ISTSOIL]
        lun_urbpoi = [false]
        lun_z_0_town = [0.0]
        lun_z_d_town = [0.0]

        CLM.set_roughness_and_forc_heights_nonlake!(fvd,
            mask_nolakec, mask_nolakep, 1:nc, 1:np,
            frac_sno, snomelt_accum, frac_veg_nosno, z0m_canopy, displa_v,
            forc_hgt_u, forc_hgt_t, forc_hgt_q,
            col_landunit, patch_gridcell, patch_landunit, patch_column,
            lun_itype, lun_urbpoi, lun_z_0_town, lun_z_d_town;
            z0param_method="ZengWang2007", use_z0m_snowmelt=false)

        # Column with snow → zsno
        @test fvd.z0mg_col[1] ≈ fvd.zsno
        # Column without snow → zlnd
        @test fvd.z0mg_col[2] ≈ fvd.zlnd

        # z0hg = z0qg = z0mg
        @test fvd.z0hg_col[1] ≈ fvd.z0mg_col[1]
        @test fvd.z0qg_col[1] ≈ fvd.z0mg_col[1]

        # Vegetation roughness = canopy z0m
        @test fvd.z0mv_patch[1] ≈ 0.5
        @test fvd.z0mv_patch[2] ≈ 0.1

        # Patch 1: has vegetation → forc_hgt uses z0mv + displa
        @test fvd.forc_hgt_u_patch[1] ≈ 10.0 + 0.5 + 2.0

        # Patch 2: bare ground → forc_hgt uses z0mg_col + displa
        @test fvd.forc_hgt_u_patch[2] ≈ 10.0 + fvd.zlnd + 0.0
    end

    @testset "set_actual_roughness_lengths!" begin
        np = 4
        nc = 2
        fvd = CLM.FrictionVelocityData()
        CLM.frictionvel_init_for_testing!(fvd, np, nc)

        fvd.z0mv_patch .= [0.5, 0.3, 0.0, 0.0]
        fvd.z0mg_col .= [0.01, 0.02]

        mask_exposedvegp = BitVector([true, false, false, false])
        mask_noexposedvegp = BitVector([false, true, false, false])
        mask_urbanp = BitVector([false, false, true, false])
        mask_lakep = BitVector([false, false, false, true])
        patch_column = [1, 1, 2, 2]
        patch_landunit = [1, 1, 2, 2]
        lun_z_0_town = [0.0, 1.5]

        CLM.set_actual_roughness_lengths!(fvd,
            mask_exposedvegp, mask_noexposedvegp, mask_urbanp, mask_lakep,
            1:np, patch_column, patch_landunit, lun_z_0_town)

        @test fvd.z0m_actual_patch[1] ≈ 0.5   # exposed veg → z0mv
        @test fvd.z0m_actual_patch[2] ≈ 0.01  # no exposed veg → z0mg_col
        @test fvd.z0m_actual_patch[3] ≈ 1.5   # urban → z_0_town
        @test fvd.z0m_actual_patch[4] ≈ 0.02  # lake → z0mg_col
    end

    @testset "frictionvel_read_nml!" begin
        fvd = CLM.FrictionVelocityData()
        @test fvd.zetamaxstable == -999.0

        # Default value
        CLM.frictionvel_read_nml!(fvd)
        @test fvd.zetamaxstable ≈ 0.5

        # Custom value
        CLM.frictionvel_read_nml!(fvd; zetamaxstable=1.0)
        @test fvd.zetamaxstable ≈ 1.0
    end

    @testset "frictionvel_read_params!" begin
        fvd = CLM.FrictionVelocityData()
        @test fvd.zsno == -999.0
        @test fvd.zlnd == -999.0
        @test fvd.zglc == -999.0

        # Default values
        CLM.frictionvel_read_params!(fvd)
        @test fvd.zsno ≈ 0.00085
        @test fvd.zlnd ≈ 0.000775
        @test fvd.zglc == -999.0  # not set without Meier2022

        # With Meier2022
        CLM.frictionvel_read_params!(fvd; zsno=0.001, zlnd=0.002, zglc=0.003,
            z0param_method="Meier2022")
        @test fvd.zsno ≈ 0.001
        @test fvd.zlnd ≈ 0.002
        @test fvd.zglc ≈ 0.003
    end

    @testset "stub functions" begin
        fvd = CLM.FrictionVelocityData()
        CLM.frictionvel_init!(fvd, 3, 2)
        @test CLM.frictionvel_init_history!(fvd, 1:3, 1:2) === nothing
        @test CLM.frictionvel_restart!(fvd, 1:3) === nothing
    end

end
