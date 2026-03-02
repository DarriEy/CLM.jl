@testset "Hydrology No Drainage Module" begin
    # ------------------------------------------------------------------
    # Tests for HydrologyNoDrainageMod port.
    # Verifies:
    #   1. update_snow_persistence!
    #   2. accumulate_snow_ice_liq!
    #   3. update_snowdp!
    #   4. compute_snow_internal_temperature!
    #   5. update_ground_and_soil_temperatures!
    #   6. update_h2osoi_vol!
    #   7. update_soilpsi!
    #   8. update_smp_l!
    #   9. compute_wf!
    #  10. update_snow_top_layer_diagnostics!
    # ------------------------------------------------------------------

    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
    nlevsno  = CLM.varpar.nlevsno
    nlevsoi  = CLM.varpar.nlevsoi
    nlevgrnd = CLM.varpar.nlevgrnd
    nlevurb  = CLM.varpar.nlevurb
    nlevmaxurbgrnd = CLM.varpar.nlevmaxurbgrnd

    nc = 4  # number of test columns
    bounds = 1:nc

    # ------------------------------------------------------------------
    # 1. update_snow_persistence!
    # ------------------------------------------------------------------
    @testset "update_snow_persistence!" begin
        snow_persistence = [100.0, 200.0, 0.0, 50.0]
        dtime = 1800.0
        mask_snow   = BitVector([true, false, false, true])
        mask_nosnow = BitVector([false, true, true, false])

        CLM.update_snow_persistence!(snow_persistence, dtime, mask_snow, mask_nosnow, bounds)

        @test snow_persistence[1] == 100.0 + 1800.0   # snow: accumulate
        @test snow_persistence[2] == 0.0                # nosnow: reset
        @test snow_persistence[3] == 0.0                # nosnow: reset
        @test snow_persistence[4] == 50.0 + 1800.0     # snow: accumulate
    end

    # ------------------------------------------------------------------
    # 2. accumulate_snow_ice_liq!
    # ------------------------------------------------------------------
    @testset "accumulate_snow_ice_liq!" begin
        nc2 = 2
        nlevtot = nlevsno + nlevmaxurbgrnd
        snowice = zeros(nc2)
        snowliq = zeros(nc2)
        h2osoi_ice = zeros(nc2, nlevtot)
        h2osoi_liq = zeros(nc2, nlevtot)
        snl = [-2, 0]  # col 1 has 2 snow layers, col 2 has none
        mask_nolake = BitVector([true, true])
        mask_snow   = BitVector([true, false])
        bounds2 = 1:nc2

        # Set snow layer values for column 1
        # snl=-2 → top snow layer at Fortran index -1 → Julia index nlevsno-1
        # active layers: Fortran -1 and 0 → Julia nlevsno-1 and nlevsno
        h2osoi_ice[1, nlevsno - 1] = 5.0   # Fortran j=-1
        h2osoi_ice[1, nlevsno]     = 10.0   # Fortran j=0
        h2osoi_liq[1, nlevsno - 1] = 1.0
        h2osoi_liq[1, nlevsno]     = 2.0

        CLM.accumulate_snow_ice_liq!(snowice, snowliq, h2osoi_ice, h2osoi_liq,
            snl, mask_nolake, mask_snow, bounds2, nlevsno)

        @test snowice[1] ≈ 15.0   # 5 + 10
        @test snowliq[1] ≈ 3.0    # 1 + 2
        @test snowice[2] ≈ 0.0    # no snow
        @test snowliq[2] ≈ 0.0
    end

    # ------------------------------------------------------------------
    # 3. update_snowdp!
    # ------------------------------------------------------------------
    @testset "update_snowdp!" begin
        snowdp     = zeros(nc)
        snow_depth = [0.5, 0.0, 1.0, 0.3]
        frac_sno_eff = [0.8, 0.0, 1.0, 0.5]

        CLM.update_snowdp!(snowdp, snow_depth, frac_sno_eff, bounds)

        @test snowdp[1] ≈ 0.4    # 0.5 * 0.8
        @test snowdp[2] ≈ 0.0
        @test snowdp[3] ≈ 1.0    # 1.0 * 1.0
        @test snowdp[4] ≈ 0.15   # 0.3 * 0.5
    end

    # ------------------------------------------------------------------
    # 4. compute_snow_internal_temperature!
    # ------------------------------------------------------------------
    @testset "compute_snow_internal_temperature!" begin
        nc2 = 2
        nlevtot = nlevsno + nlevmaxurbgrnd
        t_sno_mul_mss = zeros(nc2)
        h2osoi_ice = zeros(nc2, nlevtot)
        h2osoi_liq = zeros(nc2, nlevtot)
        t_soisno   = fill(250.0, nc2, nlevtot)
        snl = [-1, 0]
        mask_nolake = BitVector([true, true])
        mask_snow   = BitVector([true, false])
        bounds2 = 1:nc2

        # Column 1: 1 snow layer at Fortran j=0 → Julia nlevsno
        h2osoi_ice[1, nlevsno] = 20.0
        h2osoi_liq[1, nlevsno] = 5.0
        t_soisno[1, nlevsno]   = 260.0

        CLM.compute_snow_internal_temperature!(t_sno_mul_mss,
            h2osoi_ice, h2osoi_liq, t_soisno, snl,
            mask_nolake, mask_snow, bounds2, nlevsno, CLM.TFRZ)

        # Expected: 20.0 * 260.0 + 5.0 * 273.15
        expected = 20.0 * 260.0 + 5.0 * CLM.TFRZ
        @test t_sno_mul_mss[1] ≈ expected
        @test t_sno_mul_mss[2] ≈ 0.0  # no snow layers
    end

    # ------------------------------------------------------------------
    # 5. update_ground_and_soil_temperatures!
    # ------------------------------------------------------------------
    @testset "update_ground_and_soil_temperatures!" begin
        nc2 = 2
        nlevtot = nlevsno + nlevmaxurbgrnd
        nlev_zi = nlevsno + nlevmaxurbgrnd + 1
        joff = nlevsno

        t_grnd     = zeros(nc2)
        t_grnd_u   = zeros(nc2)
        t_grnd_r   = zeros(nc2)
        tsl        = zeros(nc2)
        t_soi_10cm = zeros(nc2)
        tsoi17     = zeros(nc2)

        t_soisno   = fill(280.0, nc2, nlevtot)
        t_h2osfc   = fill(275.0, nc2)
        frac_sno_eff = [0.0, 0.0]
        frac_h2osfc  = [0.0, 0.0]
        snl = [0, 0]  # no snow

        # Set up dz and zi for first soil layer
        dz = fill(0.1, nc2, nlevtot)
        zi = zeros(nc2, nlev_zi)

        # Initialize zi for soil layers: interface depths
        # Fortran zi(c, 0) = 0.0 (soil surface) → Julia zi[c, joff+1] = zi[c, nlevsno+1]
        # Fortran zi(c, j) for j=1..nlevsoi → Julia zi[c, j+nlevsno+1]
        for c in 1:nc2
            zi[c, joff + 1] = 0.0  # Fortran zi(c,0) = soil surface
            for j in 1:nlevsoi
                zi[c, j + joff + 1] = zi[c, j + joff] + dz[c, j + joff]
            end
        end

        col_landunit = [1, 1]
        col_itype = [1, 1]  # soil columns
        lun_urbpoi = BitVector([false])
        lun_itype = [CLM.ISTSOIL]
        mask_nolake = BitVector([true, true])
        bounds2 = 1:nc2

        CLM.update_ground_and_soil_temperatures!(t_grnd, t_grnd_u, t_grnd_r,
            tsl, t_soi_10cm, tsoi17,
            t_soisno, t_h2osfc, frac_sno_eff, frac_h2osfc,
            snl, dz, zi, col_landunit, col_itype,
            lun_urbpoi, lun_itype,
            mask_nolake, bounds2, nlevsoi, nlevsno)

        # No snow, no surface water → t_grnd = t_soisno(c,1) = 280.0
        @test t_grnd[1] ≈ 280.0
        @test t_grnd[2] ≈ 280.0

        # tsl = t_soisno(c,1) = 280.0
        @test tsl[1] ≈ 280.0
        @test tsl[2] ≈ 280.0

        # t_grnd_r = t_soisno(c, snl+1) = t_soisno(c,1) = 280.0 for ISTSOIL
        @test t_grnd_r[1] ≈ 280.0

        # t_soi_10cm should be 280.0 (uniform temp, divided by 0.1)
        # Since all layers are 0.1m thick and zi starts at 0, layers within 0.1m
        # contribute. zi[c, 1+joff+1] = 0.1 → zi <= 0.1 → fracl = 1.0
        # t_soi_10cm = 280.0 * 0.1 * 1.0 / 0.1 = 280.0
        @test t_soi_10cm[1] ≈ 280.0
    end

    # ------------------------------------------------------------------
    # 6. update_h2osoi_vol!
    # ------------------------------------------------------------------
    @testset "update_h2osoi_vol!" begin
        nc2 = 2
        nlevtot = nlevsno + nlevmaxurbgrnd
        joff = nlevsno

        h2osoi_vol = zeros(nc2, nlevmaxurbgrnd)
        h2osoi_liq = zeros(nc2, nlevtot)
        h2osoi_ice = zeros(nc2, nlevtot)
        dz = fill(0.1, nc2, nlevtot)

        # Set soil layer 1 (Fortran j=1 → Julia j+nlevsno in liq/ice/dz arrays)
        # 50 kg/m2 liquid, 10 kg/m2 ice in 0.1m thick layer
        h2osoi_liq[1, 1 + joff] = 50.0
        h2osoi_ice[1, 1 + joff] = 10.0
        h2osoi_liq[2, 1 + joff] = 100.0
        h2osoi_ice[2, 1 + joff] = 0.0

        col_itype = [1, 1]  # soil type (not urban)
        mask_nolake = BitVector([true, true])
        mask_urban  = BitVector([false, false])
        bounds2 = 1:nc2

        CLM.update_h2osoi_vol!(h2osoi_vol, h2osoi_liq, h2osoi_ice,
            dz, col_itype, mask_nolake, mask_urban, bounds2,
            nlevgrnd, nlevurb, nlevsno)

        # h2osoi_vol = liq/(dz*denh2o) + ice/(dz*denice)
        # col 1 layer 1: 50/(0.1*1000) + 10/(0.1*917)
        expected1 = 50.0 / (0.1 * CLM.DENH2O) + 10.0 / (0.1 * CLM.DENICE)
        @test h2osoi_vol[1, 1] ≈ expected1

        # col 2 layer 1: 100/(0.1*1000) + 0
        expected2 = 100.0 / (0.1 * CLM.DENH2O)
        @test h2osoi_vol[2, 1] ≈ expected2
    end

    # ------------------------------------------------------------------
    # 7. update_soilpsi!
    # ------------------------------------------------------------------
    @testset "update_soilpsi!" begin
        nc2 = 2
        nlevtot = nlevsno + nlevmaxurbgrnd
        joff = nlevsno

        soilpsi = fill(-999.0, nc2, nlevgrnd)
        h2osoi_liq = zeros(nc2, nlevtot)
        dz = fill(0.1, nc2, nlevtot)
        watsat = fill(0.4, nc2, nlevgrnd)
        sucsat = fill(200.0, nc2, nlevgrnd)
        bsw = fill(5.0, nc2, nlevgrnd)

        # Column 1: positive liquid water → calculate psi
        h2osoi_liq[1, 1 + joff] = 30.0  # 30 kg/m2 in 0.1m → vwc = 0.3
        # Column 2: zero liquid → soilpsi = -15.0
        h2osoi_liq[2, 1 + joff] = 0.0

        mask_hydrology = BitVector([true, true])
        bounds2 = 1:nc2

        CLM.update_soilpsi!(soilpsi, h2osoi_liq, dz, watsat, sucsat, bsw,
            mask_hydrology, bounds2, nlevgrnd, nlevsno)

        # Column 1: vwc = 30/(0.1*1000) = 0.3
        # fsattmp = max(0.3/0.4, 0.001) = 0.75
        # psi = 200 * (-9.8e-6) * (0.75)^(-5.0) = 200 * (-9.8e-6) * (0.75)^(-5)
        vwc = 0.3
        fsattmp = max(vwc / 0.4, 0.001)
        psi = 200.0 * (-9.8e-6) * (fsattmp)^(-5.0)
        expected = min(max(psi, -15.0), 0.0)
        @test soilpsi[1, 1] ≈ expected

        # Column 2: no liquid
        @test soilpsi[2, 1] ≈ -15.0
    end

    # ------------------------------------------------------------------
    # 8. update_smp_l!
    # ------------------------------------------------------------------
    @testset "update_smp_l!" begin
        nc2 = 2
        smp_l = fill(-999.0, nc2, nlevgrnd)
        h2osoi_vol = fill(0.3, nc2, nlevmaxurbgrnd)
        watsat = fill(0.4, nc2, nlevgrnd)
        sucsat = fill(200.0, nc2, nlevgrnd)
        bsw = fill(5.0, nc2, nlevgrnd)
        smpmin = fill(-1.0e8, nc2)

        mask_hydrology = BitVector([true, true])
        bounds2 = 1:nc2

        CLM.update_smp_l!(smp_l, h2osoi_vol, watsat, sucsat, bsw, smpmin,
            mask_hydrology, bounds2, nlevgrnd)

        # s_node = max(0.3/0.4, 0.01) = 0.75
        # smp_l = -200 * 0.75^(-5.0)
        s_node = max(0.3 / 0.4, 0.01)
        s_node = min(1.0, s_node)
        expected = -200.0 * s_node^(-5.0)
        expected = max(-1.0e8, expected)
        @test smp_l[1, 1] ≈ expected
        @test smp_l[2, 1] ≈ expected
    end

    # ------------------------------------------------------------------
    # 9. compute_wf!
    # ------------------------------------------------------------------
    @testset "compute_wf!" begin
        nc2 = 2
        nlevtot = nlevsno + nlevmaxurbgrnd
        joff = nlevsno

        wf_out = zeros(nc2)
        h2osoi_vol = fill(0.3, nc2, nlevmaxurbgrnd)
        watsat = fill(0.4, nc2, nlevgrnd)
        sucsat = fill(200.0, nc2, nlevgrnd)
        bsw = fill(5.0, nc2, nlevgrnd)

        # z and dz are snow+soil arrays
        z  = zeros(nc2, nlevtot)
        dz = fill(0.02, nc2, nlevtot)  # thin layers so they fit within 0.05m

        # Set soil layer depths (Fortran j=1..nlevgrnd, Julia j+nlevsno)
        for c in 1:nc2
            for j in 1:nlevgrnd
                z[c, j + joff] = 0.01 + (j - 1) * 0.02  # center depth
            end
        end

        mask_hydrology = BitVector([true, true])
        bounds2 = 1:nc2

        CLM.compute_wf!(wf_out, h2osoi_vol, watsat, sucsat, bsw,
            z, dz, mask_hydrology, bounds2, nlevgrnd, nlevsno, 0.05)

        # z[c,1+joff]+0.5*dz[c,1+joff] = 0.01 + 0.01 = 0.02 <= 0.05 ✓
        # z[c,2+joff]+0.5*dz[c,2+joff] = 0.03 + 0.01 = 0.04 <= 0.05 ✓
        # z[c,3+joff]+0.5*dz[c,3+joff] = 0.05 + 0.01 = 0.06 > 0.05 ✗
        # So layers 1 and 2 contribute
        watdry = 0.4 * (316230.0 / 200.0)^(-1.0 / 5.0)
        rwat = 2 * (0.3 - watdry) * 0.02
        swat = 2 * (0.4 - watdry) * 0.02
        rz_val = 2 * 0.02
        tsw = rwat / rz_val
        stsw = swat / rz_val
        expected = tsw / stsw
        @test wf_out[1] ≈ expected
    end

    # ------------------------------------------------------------------
    # 10. update_snow_top_layer_diagnostics!
    # ------------------------------------------------------------------
    @testset "update_snow_top_layer_diagnostics!" begin
        nc2 = 2
        nlevtot = nlevsno + nlevmaxurbgrnd

        h2osno_top = zeros(nc2)
        snw_rds = fill(100.0, nc2, nlevsno)
        snot_top = zeros(nc2)
        dTdz_top = zeros(nc2)
        snw_rds_top = zeros(nc2)
        sno_liq_top = zeros(nc2)
        h2osoi_ice = zeros(nc2, nlevtot)
        h2osoi_liq = zeros(nc2, nlevtot)
        snl = [-1, 0]
        mask_snow   = BitVector([true, false])
        mask_nosnow = BitVector([false, true])
        bounds2 = 1:nc2

        # Column 1: 1 snow layer at Fortran j=0 → Julia nlevsno
        h2osoi_ice[1, nlevsno] = 20.0
        h2osoi_liq[1, nlevsno] = 3.0

        CLM.update_snow_top_layer_diagnostics!(h2osno_top, snw_rds,
            snot_top, dTdz_top, snw_rds_top, sno_liq_top,
            h2osoi_ice, h2osoi_liq, snl,
            mask_snow, mask_nosnow, bounds2, nlevsno, CLM.SPVAL)

        # Column 1: snow → h2osno_top = ice + liq
        @test h2osno_top[1] ≈ 23.0

        # Column 2: no snow → zeroed and spval diagnostics
        @test h2osno_top[2] ≈ 0.0
        @test all(snw_rds[2, :] .== 0.0)
        @test snot_top[2] ≈ CLM.SPVAL
        @test dTdz_top[2] ≈ CLM.SPVAL
        @test snw_rds_top[2] ≈ CLM.SPVAL
        @test sno_liq_top[2] ≈ CLM.SPVAL
    end
end
