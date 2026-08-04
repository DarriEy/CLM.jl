@testset "Atm2LndData" begin

    @testset "Atm2LndParamsData default construction" begin
        p = CLM.Atm2LndParamsData()
        @test p.repartition_rain_snow == false
        @test p.glcmec_downscale_longwave == false
        @test isnan(p.lapse_rate)
        @test isnan(p.lapse_rate_longwave)
        @test isnan(p.longwave_downscaling_limit)
        @test isnan(p.precip_repartition_glc_all_snow_t)
        @test isnan(p.precip_repartition_glc_frac_rain_slope)
        @test isnan(p.precip_repartition_nonglc_all_snow_t)
        @test isnan(p.precip_repartition_nonglc_frac_rain_slope)
    end

    @testset "Atm2LndData default construction" begin
        a = CLM.Atm2LndData()
        @test length(a.forc_u_grc) == 0
        @test length(a.forc_t_downscaled_col) == 0
        @test length(a.fsd24_patch) == 0
        @test size(a.forc_solad_not_downscaled_grc) == (0, 0)
        @test size(a.forc_aer_grc) == (0, 0)
    end

    @testset "compute_ramp_params" begin
        # all_snow_t_c = -5, all_rain_t_c = 5 => slope = 1/10 = 0.1
        (snow_k, slope) = CLM.compute_ramp_params(-5.0, 5.0)
        @test snow_k ≈ -5.0 + CLM.TFRZ
        @test slope ≈ 0.1

        # all_snow_t_c = 0, all_rain_t_c = 2 => slope = 0.5
        (snow_k2, slope2) = CLM.compute_ramp_params(0.0, 2.0)
        @test snow_k2 ≈ CLM.TFRZ
        @test slope2 ≈ 0.5
    end

    @testset "atm2lnd_params_init! basic (no repartition, no downscale)" begin
        p = CLM.Atm2LndParamsData()
        CLM.atm2lnd_params_init!(p;
            repartition_rain_snow = false,
            glcmec_downscale_longwave = false,
            lapse_rate = 0.006)

        @test p.repartition_rain_snow == false
        @test p.glcmec_downscale_longwave == false
        @test p.lapse_rate ≈ 0.006
        @test isnan(p.lapse_rate_longwave)
        @test isnan(p.longwave_downscaling_limit)
        @test isnan(p.precip_repartition_glc_all_snow_t)
        @test isnan(p.precip_repartition_glc_frac_rain_slope)
        @test isnan(p.precip_repartition_nonglc_all_snow_t)
        @test isnan(p.precip_repartition_nonglc_frac_rain_slope)
    end

    @testset "atm2lnd_params_init! with longwave downscaling" begin
        p = CLM.Atm2LndParamsData()
        CLM.atm2lnd_params_init!(p;
            repartition_rain_snow = false,
            glcmec_downscale_longwave = true,
            lapse_rate = 0.006,
            lapse_rate_longwave = 0.032,
            longwave_downscaling_limit = 0.5)

        @test p.glcmec_downscale_longwave == true
        @test p.lapse_rate_longwave ≈ 0.032
        @test p.longwave_downscaling_limit ≈ 0.5
    end

    @testset "atm2lnd_params_init! longwave downscale error checks" begin
        p = CLM.Atm2LndParamsData()
        # Missing lapse_rate_longwave
        @test_throws ErrorException CLM.atm2lnd_params_init!(p;
            repartition_rain_snow = false,
            glcmec_downscale_longwave = true,
            lapse_rate = 0.006)

        # Missing longwave_downscaling_limit
        @test_throws ErrorException CLM.atm2lnd_params_init!(p;
            repartition_rain_snow = false,
            glcmec_downscale_longwave = true,
            lapse_rate = 0.006,
            lapse_rate_longwave = 0.032)

        # longwave_downscaling_limit out of range
        @test_throws ErrorException CLM.atm2lnd_params_init!(p;
            repartition_rain_snow = false,
            glcmec_downscale_longwave = true,
            lapse_rate = 0.006,
            lapse_rate_longwave = 0.032,
            longwave_downscaling_limit = 1.5)
    end

    @testset "atm2lnd_params_init! with repartitioning" begin
        p = CLM.Atm2LndParamsData()
        CLM.atm2lnd_params_init!(p;
            repartition_rain_snow = true,
            glcmec_downscale_longwave = false,
            lapse_rate = 0.006,
            precip_repartition_glc_all_snow_t = -2.0,
            precip_repartition_glc_all_rain_t = 0.0,
            precip_repartition_nonglc_all_snow_t = -5.0,
            precip_repartition_nonglc_all_rain_t = 2.5)

        @test p.repartition_rain_snow == true
        # glc: all_snow_t_k = -2.0 + TFRZ, slope = 1/(0-(-2)) = 0.5
        @test p.precip_repartition_glc_all_snow_t ≈ -2.0 + CLM.TFRZ
        @test p.precip_repartition_glc_frac_rain_slope ≈ 0.5
        # nonglc: all_snow_t_k = -5.0 + TFRZ, slope = 1/(2.5-(-5)) = 1/7.5
        @test p.precip_repartition_nonglc_all_snow_t ≈ -5.0 + CLM.TFRZ
        @test p.precip_repartition_nonglc_frac_rain_slope ≈ 1.0 / 7.5
    end

    @testset "atm2lnd_params_init! repartitioning error checks" begin
        p = CLM.Atm2LndParamsData()
        # Missing required param
        @test_throws ErrorException CLM.atm2lnd_params_init!(p;
            repartition_rain_snow = true,
            glcmec_downscale_longwave = false,
            lapse_rate = 0.006)

        # all_rain <= all_snow (glc)
        @test_throws ErrorException CLM.atm2lnd_params_init!(p;
            repartition_rain_snow = true,
            glcmec_downscale_longwave = false,
            lapse_rate = 0.006,
            precip_repartition_glc_all_snow_t = 0.0,
            precip_repartition_glc_all_rain_t = -1.0,
            precip_repartition_nonglc_all_snow_t = -5.0,
            precip_repartition_nonglc_all_rain_t = 2.5)

        # all_rain <= all_snow (nonglc)
        @test_throws ErrorException CLM.atm2lnd_params_init!(p;
            repartition_rain_snow = true,
            glcmec_downscale_longwave = false,
            lapse_rate = 0.006,
            precip_repartition_glc_all_snow_t = -2.0,
            precip_repartition_glc_all_rain_t = 0.0,
            precip_repartition_nonglc_all_snow_t = 2.5,
            precip_repartition_nonglc_all_rain_t = 2.5)
    end

    @testset "atm2lnd_init! allocation" begin
        ng, nc, np = 4, 8, 12
        a = CLM.Atm2LndData()
        CLM.atm2lnd_init!(a, ng, nc, np)

        # Gridcell-level vectors
        @test length(a.forc_u_grc) == ng
        @test all(x -> x == 0.0, a.forc_u_grc)
        @test length(a.forc_v_grc) == ng
        @test length(a.forc_wind_grc) == ng
        @test length(a.forc_hgt_grc) == ng
        @test length(a.forc_topo_grc) == ng
        @test length(a.forc_hgt_u_grc) == ng
        @test length(a.forc_hgt_t_grc) == ng
        @test length(a.forc_hgt_q_grc) == ng
        @test length(a.forc_vp_grc) == ng
        @test length(a.forc_pco2_grc) == ng
        @test length(a.forc_ndep_grc) == ng
        @test length(a.forc_pc13o2_grc) == ng
        @test length(a.forc_po2_grc) == ng
        @test length(a.forc_pch4_grc) == ng
        @test length(a.forc_o3_grc) == ng
        @test length(a.forc_solar_not_downscaled_grc) == ng

        # Gridcell-level not-downscaled
        @test length(a.forc_t_not_downscaled_grc) == ng
        @test length(a.forc_pbot_not_downscaled_grc) == ng
        @test length(a.forc_th_not_downscaled_grc) == ng
        @test length(a.forc_rho_not_downscaled_grc) == ng
        @test length(a.forc_lwrad_not_downscaled_grc) == ng

        # 2D gridcell arrays
        @test size(a.forc_solad_not_downscaled_grc) == (ng, CLM.NUMRAD)
        @test size(a.forc_solai_grc) == (ng, CLM.NUMRAD)
        @test size(a.forc_aer_grc) == (ng, 14)
        @test all(x -> x == 0.0, a.forc_aer_grc)

        # Column-level downscaled
        @test length(a.forc_t_downscaled_col) == nc
        @test all(x -> x == 0.0, a.forc_t_downscaled_col)
        @test length(a.forc_pbot_downscaled_col) == nc
        @test length(a.forc_th_downscaled_col) == nc
        @test length(a.forc_rho_downscaled_col) == nc
        @test length(a.forc_lwrad_downscaled_col) == nc
        @test size(a.forc_solad_downscaled_col) == (nc, CLM.NUMRAD)
        @test length(a.forc_solar_downscaled_col) == nc

        # Patch-level time-averaged
        @test length(a.fsd24_patch) == np
        @test all(isnan, a.fsd24_patch)
        @test length(a.fsd240_patch) == np
        @test all(isnan, a.fsd240_patch)
        @test length(a.fsi24_patch) == np
        @test all(isnan, a.fsi24_patch)
        @test length(a.fsi240_patch) == np
        @test all(isnan, a.fsi240_patch)
        @test length(a.t_mo_patch) == np
        @test all(isnan, a.t_mo_patch)
        @test length(a.t_mo_min_patch) == np
        @test all(x -> x == CLM.SPVAL, a.t_mo_min_patch)
    end

    @testset "atm2lnd_init_for_testing! default params" begin
        ng, nc, np = 3, 5, 7
        a = CLM.Atm2LndData()
        CLM.atm2lnd_init_for_testing!(a, ng, nc, np)

        @test a.params.repartition_rain_snow == false
        @test a.params.glcmec_downscale_longwave == false
        @test a.params.lapse_rate ≈ 0.01
        @test length(a.forc_u_grc) == ng
    end

    @testset "atm2lnd_init_for_testing! with custom params" begin
        ng, nc, np = 2, 3, 4
        custom_params = CLM.Atm2LndParamsData(
            repartition_rain_snow = true,
            glcmec_downscale_longwave = false,
            lapse_rate = 0.005,
            precip_repartition_glc_all_snow_t = 270.0,
            precip_repartition_glc_frac_rain_slope = 0.5)

        a = CLM.Atm2LndData()
        CLM.atm2lnd_init_for_testing!(a, ng, nc, np; params = custom_params)

        @test a.params.repartition_rain_snow == true
        @test a.params.lapse_rate ≈ 0.005
        @test a.params.precip_repartition_glc_all_snow_t ≈ 270.0
        @test length(a.forc_u_grc) == ng
    end

    @testset "atm2lnd_read_namelist!" begin
        a = CLM.Atm2LndData()
        CLM.atm2lnd_init!(a, 3, 5, 7)
        CLM.atm2lnd_read_namelist!(a;
            repartition_rain_snow = false,
            glcmec_downscale_longwave = false,
            lapse_rate = 0.008)
        @test a.params.lapse_rate ≈ 0.008
        @test a.params.repartition_rain_snow == false
    end

    @testset "atm2lnd_update_acc_vars! basic" begin
        ng, nc, np = 2, 4, 6
        a = CLM.Atm2LndData()
        CLM.atm2lnd_init!(a, ng, nc, np)

        # Set some forcing values
        a.forc_solad_not_downscaled_grc[1, 1] = 100.0
        a.forc_solad_not_downscaled_grc[2, 1] = 200.0
        a.forc_solai_grc[1, 1] = 50.0
        a.forc_solai_grc[2, 1] = 75.0

        # patches 1-3 on gridcell 1, patches 4-6 on gridcell 2
        patch_gridcell = [1, 1, 1, 2, 2, 2]
        patch_column   = [1, 2, 3, 4, 4, 4]

        CLM.atm2lnd_update_acc_vars!(a, 1:np, patch_gridcell, patch_column)

        # Check direct beam
        @test a.fsd24_patch[1] ≈ 100.0
        @test a.fsd240_patch[3] ≈ 100.0
        @test a.fsd24_patch[4] ≈ 200.0
        @test a.fsd240_patch[6] ≈ 200.0

        # Check diffuse
        @test a.fsi24_patch[1] ≈ 50.0
        @test a.fsi240_patch[2] ≈ 50.0
        @test a.fsi24_patch[5] ≈ 75.0
        @test a.fsi240_patch[6] ≈ 75.0
    end

    @testset "atm2lnd_clean!" begin
        ng, nc, np = 3, 5, 7
        a = CLM.Atm2LndData()
        CLM.atm2lnd_init!(a, ng, nc, np)

        @test length(a.forc_u_grc) == ng

        CLM.atm2lnd_clean!(a)

        @test length(a.forc_u_grc) == 0
        @test length(a.forc_t_downscaled_col) == 0
        @test length(a.fsd24_patch) == 0
        @test size(a.forc_solad_not_downscaled_grc) == (0, 0)
        @test size(a.forc_aer_grc) == (0, 0)
        @test length(a.t_mo_patch) == 0
        @test length(a.t_mo_min_patch) == 0
    end

    @testset "stub functions run without error" begin
        a = CLM.Atm2LndData()
        CLM.atm2lnd_init!(a, 3, 5, 7)
        @test CLM.atm2lnd_init_history!(a, 3, 5, 7) === nothing
        @test CLM.atm2lnd_init_acc_buffer!(a) === nothing
        @test CLM.atm2lnd_init_acc_vars!(a, 1:7) === nothing
        @test CLM.atm2lnd_restart!(a, 1:7) === nothing
    end

    @testset "field mutability" begin
        a = CLM.Atm2LndData()
        CLM.atm2lnd_init!(a, 3, 5, 7)

        a.forc_u_grc[1] = 5.0
        @test a.forc_u_grc[1] == 5.0

        a.forc_t_downscaled_col[2] = 300.0
        @test a.forc_t_downscaled_col[2] == 300.0

        a.forc_aer_grc[1, 3] = 0.001
        @test a.forc_aer_grc[1, 3] == 0.001
    end

    @testset "re-init overwrites previous state" begin
        a = CLM.Atm2LndData()
        CLM.atm2lnd_init!(a, 3, 5, 7)
        a.forc_u_grc[1] = 999.0

        CLM.atm2lnd_init!(a, 10, 20, 30)
        @test length(a.forc_u_grc) == 10
        @test all(x -> x == 0.0, a.forc_u_grc)
        @test length(a.forc_t_downscaled_col) == 20
        @test length(a.fsd24_patch) == 30
    end

end

# ==========================================================================
# rof -> lnd feedback state (Wateratm2lndType.F90). The coupler fills these from
# the river model's export bundle (Flrr_volr / Flrr_volrmch / Flrr_flood); they
# are the ONLY river->land feedback CTSM has. Zero in an uncoupled run, which is
# the Fortran `ival = 0` initialization and reproduces standalone CLM exactly.
# ==========================================================================
@testset "atm2lnd rof->lnd feedback fields (volr / forc_flood)" begin
    ng, nc, np = 3, 6, 9
    a = CLM.Atm2LndData()
    @test length(a.volr_grc) == 0
    @test length(a.forc_flood_grc) == 0

    CLM.atm2lnd_init!(a, ng, nc, np)
    for f in (:volr_grc, :volrmch_grc, :forc_flood_grc)
        v = getfield(a, f)
        @test length(v) == ng
        @test all(iszero, v)      # uncoupled default
    end

    CLM.atm2lnd_clean!(a)
    @test length(a.volr_grc) == 0
    @test length(a.volrmch_grc) == 0
    @test length(a.forc_flood_grc) == 0
end

# ==========================================================================
# downscale_forcings! rescales forc_pco2 / forc_po2 / forc_pc13o2 — GRIDCELL
# fields — onto the elevation-corrected pressure. It is a read-modify-write, so
# it must be applied exactly ONCE per gridcell however many columns the gridcell
# has. A per-column version compounded the ratio: with n columns the implied
# pressure came out as pbot*(pbot_ds/pbot_nd)^n (a ~50% CO2 error at the 11
# columns of a glacier elevation-class subgrid), and the pco2/po2 molar-ratio
# invariant did NOT catch it because both fields compound identically.
# ==========================================================================
@testset "downscale_forcings!: the pco2 rescale is once per gridcell" begin
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
    CLM.varcon_init!()
    CLM.varpar.nlevsno = 5; CLM.varpar.nlevsoi = 10
    CLM.varpar.nlevgrnd = 10; CLM.varpar.nlevmaxurbgrnd = 10

    # One gridcell; `ncol` columns, all at the same elevation so every column
    # gets the same downscaled pressure and the correct answer is unambiguous.
    function run_downscale(ncol::Int; set_weights::Bool = true, elev = 500.0)
        inst = CLM.CLMInstances()
        CLM.clm_instInit!(inst; ng = 1, nl = 1, nc = ncol, np = ncol,
                          nlevdecomp_full = CLM.varpar.nlevdecomp_full)
        col = inst.column; lun = inst.landunit; a = inst.atm2lnd
        for k in 1:ncol
            col.gridcell[k] = 1; col.landunit[k] = 1; col.active[k] = true
            col.itype[k] = CLM.ISTSOIL
            if set_weights
                col.wtgcell[k] = 1 / ncol; col.wtlunit[k] = 1 / ncol
            end
        end
        lun.itype[1] = CLM.ISTSOIL; lun.gridcell[1] = 1
        lun.active[1] = true; lun.urbpoi[1] = false; lun.wtgcell[1] = 1.0

        pbot = 9.9e4
        a.forc_pbot_not_downscaled_grc[1] = pbot
        a.forc_t_not_downscaled_grc[1]  = 288.0
        a.forc_th_not_downscaled_grc[1] = 288.0
        a.forc_q_not_downscaled_grc[1]  = 0.008
        a.forc_topo_grc[1] = 0.0
        a.forc_pco2_grc[1]   = 367.0e-6 * pbot
        a.forc_po2_grc[1]    = 0.209 * pbot
        a.forc_pc13o2_grc[1] = 367.0e-6 * pbot * 0.01
        inst.topo.topo_col .= elev

        bounds = CLM.BoundsType(begg = 1, endg = 1, begl = 1, endl = 1,
                                begc = 1, endc = ncol, begp = 1, endp = ncol,
                                begCohort = 0, endCohort = 0,
                                level = CLM.BOUNDS_LEVEL_CLUMP, clump_index = 1)
        CLM.downscale_forcings!(bounds, a, col, lun, inst.topo)
        return (pbot = pbot, pbot_ds = a.forc_pbot_downscaled_col[1],
                pco2 = a.forc_pco2_grc[1], po2 = a.forc_po2_grc[1],
                pc13o2 = a.forc_pc13o2_grc[1])
    end

    r1 = run_downscale(1)
    @test r1.pbot_ds < r1.pbot                       # 500 m up: lower pressure
    @test r1.pco2 / 367.0e-6 ≈ r1.pbot_ds rtol = 1e-12

    # The whole point: 11 columns (a glc_nec = 10 subgrid) must give the SAME
    # gridcell CO2 as one column, not the ratio raised to the 11th power.
    r11 = run_downscale(11)
    @test r11.pco2 ≈ r1.pco2 rtol = 1e-12
    @test r11.pco2 / 367.0e-6 ≈ r11.pbot_ds rtol = 1e-12
    ratio = r11.pbot_ds / r11.pbot
    @test !isapprox(r11.pco2 / 367.0e-6, r11.pbot * ratio^11; rtol = 1e-6)

    # O2 and C13 ride the same single rescale.
    @test r11.po2 / 0.209 ≈ r11.pbot_ds rtol = 1e-12
    @test r11.pc13o2 / (367.0e-6 * 0.01) ≈ r11.pbot_ds rtol = 1e-12
    # ...so the molar-ratio invariant holds too (it held before the fix as well,
    # which is exactly why it cannot be the only assertion here).
    @test r11.pco2 / r11.po2 ≈ 367.0e-6 / 0.209 rtol = 1e-12

    # `column_init!` allocates wtgcell as NaN and CLM.jl's hand-built driver
    # fixtures leave it unset; that must degrade to an unweighted mean, not NaN.
    rnan = run_downscale(4; set_weights = false)
    @test isfinite(rnan.pco2)
    @test rnan.pco2 / 367.0e-6 ≈ rnan.pbot_ds rtol = 1e-12

    # ZERO weight is NOT the same as unset. A zero-area column is CTSM's
    # "virtual" column — a glacier elevation-class column that exists only so an
    # SMB can be computed at that class elevation — and it must contribute
    # NOTHING to the gridcell mean. Treating zero like unset (weight 1) drags
    # the mean up the whole elevation ladder and biases CO2 in gridcells that
    # carry no glacier at all.
    @testset "zero-area virtual columns do not pull the gridcell mean" begin
        ncol = 6
        inst = CLM.CLMInstances()
        CLM.clm_instInit!(inst; ng = 1, nl = 1, nc = ncol, np = ncol,
                          nlevdecomp_full = CLM.varpar.nlevdecomp_full)
        col = inst.column; lun = inst.landunit; a = inst.atm2lnd
        for k in 1:ncol
            col.gridcell[k] = 1; col.landunit[k] = 1; col.active[k] = true
            col.itype[k] = CLM.ISTSOIL
            # column 1 is the whole gridcell; 2..ncol are virtual (zero area)
            col.wtgcell[k] = (k == 1) ? 1.0 : 0.0
            col.wtlunit[k] = col.wtgcell[k]
        end
        lun.itype[1] = CLM.ISTSOIL; lun.gridcell[1] = 1
        lun.active[1] = true; lun.urbpoi[1] = false; lun.wtgcell[1] = 1.0

        pbot = 9.9e4
        a.forc_pbot_not_downscaled_grc[1] = pbot
        a.forc_t_not_downscaled_grc[1]  = 288.0
        a.forc_th_not_downscaled_grc[1] = 288.0
        a.forc_q_not_downscaled_grc[1]  = 0.008
        a.forc_topo_grc[1] = 0.0
        a.forc_pco2_grc[1] = 367.0e-6 * pbot
        a.forc_po2_grc[1]  = 0.209 * pbot
        # the real column sits at the atmosphere's own height (no correction);
        # the virtual ones are strung far up the elevation-class ladder
        inst.topo.topo_col[1] = 0.0
        for k in 2:ncol
            inst.topo.topo_col[k] = 1000.0 * (k - 1)
        end

        bounds = CLM.BoundsType(begg = 1, endg = 1, begl = 1, endl = 1,
                                begc = 1, endc = ncol, begp = 1, endp = ncol,
                                begCohort = 0, endCohort = 0,
                                level = CLM.BOUNDS_LEVEL_CLUMP, clump_index = 1)
        CLM.downscale_forcings!(bounds, a, col, lun, inst.topo)

        # zero elevation difference on the only column with area => no rescale
        @test a.forc_pco2_grc[1] / 367.0e-6 ≈ pbot rtol = 1e-12
        @test a.forc_po2_grc[1] / 0.209 ≈ pbot rtol = 1e-12
    end
end
