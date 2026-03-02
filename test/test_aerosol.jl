@testset "Aerosol" begin
    # Setup: initialize varpar so nlevsno is valid
    vp = CLM.varpar
    saved_nlevsno = vp.nlevsno
    vp.nlevsno = 5
    nlevsno = vp.nlevsno

    nc = 4   # number of columns
    ng = 2   # number of gridcells

    @testset "AerosolData struct construction" begin
        aer = CLM.AerosolData()
        @test size(aer.mss_bcpho_col) == (0, 0)
        @test length(aer.flx_bc_dep_col) == 0
    end

    @testset "aerosol_init! allocation" begin
        aer = CLM.AerosolData()
        CLM.aerosol_init!(aer, nc)

        # 2D arrays: (nc, nlevsno)
        @test size(aer.mss_bcpho_col) == (nc, nlevsno)
        @test size(aer.mss_bcphi_col) == (nc, nlevsno)
        @test size(aer.mss_bctot_col) == (nc, nlevsno)
        @test size(aer.mss_ocpho_col) == (nc, nlevsno)
        @test size(aer.mss_dst1_col)  == (nc, nlevsno)
        @test size(aer.mss_dsttot_col) == (nc, nlevsno)
        @test size(aer.mss_cnc_bcphi_col) == (nc, nlevsno)
        @test size(aer.mss_cnc_dst4_col)  == (nc, nlevsno)

        # 1D arrays: (nc,)
        @test length(aer.mss_bc_col_col) == nc
        @test length(aer.mss_bc_top_col) == nc
        @test length(aer.flx_bc_dep_col) == nc
        @test length(aer.flx_dst_dep_col) == nc
        @test length(aer.flx_oc_dep_col) == nc

        # All initialized to NaN
        @test all(isnan, aer.mss_bcpho_col)
        @test all(isnan, aer.flx_bc_dep_col)
    end

    @testset "aerosol_init_cold! zeros masses" begin
        aer = CLM.AerosolData()
        CLM.aerosol_init!(aer, nc)
        CLM.aerosol_init_cold!(aer, 1:nc)

        @test all(aer.mss_bcpho_col .== 0.0)
        @test all(aer.mss_bcphi_col .== 0.0)
        @test all(aer.mss_bctot_col .== 0.0)
        @test all(aer.mss_ocpho_col .== 0.0)
        @test all(aer.mss_dst1_col  .== 0.0)
        @test all(aer.mss_dsttot_col .== 0.0)
        @test all(aer.mss_cnc_bcphi_col .== 0.0)
        @test all(aer.mss_cnc_dst4_col  .== 0.0)
    end

    @testset "aerosol_reset! zeros single column" begin
        aer = CLM.AerosolData()
        CLM.aerosol_init!(aer, nc)

        # Set some non-zero values
        aer.mss_bcpho_col[2, :] .= 1.0
        aer.mss_bc_col_col[2] = 5.0
        aer.mss_dst_top_col[2] = 3.0

        CLM.aerosol_reset!(aer, 2)

        @test all(aer.mss_bcpho_col[2, :] .== 0.0)
        @test aer.mss_bc_col_col[2] == 0.0
        @test aer.mss_dst_top_col[2] == 0.0
        # Column 1 should be unchanged (still NaN from init)
        @test all(isnan, aer.mss_bcpho_col[1, :])
    end

    @testset "aerosol_reset_filter! with mask" begin
        aer = CLM.AerosolData()
        CLM.aerosol_init!(aer, nc)
        CLM.aerosol_init_cold!(aer, 1:nc)

        # Set non-zero values for columns 1 and 3
        aer.mss_bcpho_col[1, :] .= 2.0
        aer.mss_bcpho_col[3, :] .= 4.0
        aer.mss_bc_col_col[1] = 10.0
        aer.mss_bc_col_col[3] = 20.0

        mask = BitVector([true, false, true, false])
        CLM.aerosol_reset_filter!(aer, mask, 1:nc)

        @test all(aer.mss_bcpho_col[1, :] .== 0.0)
        @test all(aer.mss_bcpho_col[3, :] .== 0.0)
        @test aer.mss_bc_col_col[1] == 0.0
        @test aer.mss_bc_col_col[3] == 0.0
        # Columns 2,4 unchanged (0 from cold init)
        @test all(aer.mss_bcpho_col[2, :] .== 0.0)
    end

    @testset "aerosol_clean! deallocates" begin
        aer = CLM.AerosolData()
        CLM.aerosol_init!(aer, nc)
        CLM.aerosol_clean!(aer)

        @test size(aer.mss_bcpho_col) == (0, 0)
        @test length(aer.flx_bc_dep_col) == 0
        @test length(aer.mss_bc_col_col) == 0
    end

    @testset "aerosol_masses! with snow" begin
        aer = CLM.AerosolData()
        CLM.aerosol_init!(aer, nc)
        CLM.aerosol_init_cold!(aer, 1:nc)

        # Setup: column 1 has 3 snow layers (snl = -3)
        # Active layers: Fortran j in [-2, -1, 0] → Julia j in [3, 4, 5]
        snl = [-3, 0, -1, 0]

        # h2osoi arrays: (nc, nlevsno + nlevmaxurbgrnd), but we only need snow part
        # For simplicity, make them (nc, nlevsno) since we only access snow layers
        nlevmaxurbgrnd = 25
        ntot = nlevsno + nlevmaxurbgrnd
        h2osoi_ice = zeros(nc, ntot)
        h2osoi_liq = zeros(nc, ntot)
        h2osno_top = zeros(nc)
        snw_rds = zeros(nc, nlevsno)

        # Set snow mass in active layers of column 1 (Julia indices 3,4,5)
        h2osoi_ice[1, 3] = 10.0
        h2osoi_ice[1, 4] = 20.0
        h2osoi_ice[1, 5] = 30.0
        h2osoi_liq[1, 3] = 1.0
        h2osoi_liq[1, 4] = 2.0
        h2osoi_liq[1, 5] = 3.0

        # Set some aerosol masses in active layers
        aer.mss_bcpho_col[1, 3] = 0.001
        aer.mss_bcphi_col[1, 3] = 0.002
        aer.mss_bcpho_col[1, 4] = 0.003
        aer.mss_bcphi_col[1, 4] = 0.004
        aer.mss_bcpho_col[1, 5] = 0.005
        aer.mss_bcphi_col[1, 5] = 0.006

        aer.mss_dst1_col[1, 4] = 0.01
        aer.mss_dst2_col[1, 4] = 0.02

        # Column 3 has 1 snow layer (snl = -1), active at Julia index 5
        h2osoi_ice[3, 5] = 50.0
        h2osoi_liq[3, 5] = 5.0
        aer.mss_ocpho_col[3, 5] = 0.1
        aer.mss_ocphi_col[3, 5] = 0.2

        # mask_on = columns with snow, mask_off = columns without
        mask_on  = BitVector([true, false, true, false])
        mask_off = BitVector([false, true, false, true])

        CLM.aerosol_masses!(aer, mask_on, mask_off, 1:nc,
                            snl, h2osoi_ice, h2osoi_liq, h2osno_top, snw_rds)

        # Check column 1
        # Total BC per layer = pho + phi
        @test aer.mss_bctot_col[1, 3] ≈ 0.003  # 0.001 + 0.002
        @test aer.mss_bctot_col[1, 4] ≈ 0.007  # 0.003 + 0.004
        @test aer.mss_bctot_col[1, 5] ≈ 0.011  # 0.005 + 0.006

        # Column-integrated BC
        @test aer.mss_bc_col_col[1] ≈ 0.021  # 0.003 + 0.007 + 0.011

        # Top layer = snl+1+nlevsno = -3+1+5 = 3
        @test aer.mss_bc_top_col[1] ≈ 0.003
        @test h2osno_top[1] ≈ 11.0  # 10.0 + 1.0

        # Mass concentration in layer 4: bcphi / snowmass
        snowmass_4 = 20.0 + 2.0
        @test aer.mss_cnc_bcphi_col[1, 4] ≈ 0.004 / snowmass_4
        @test aer.mss_cnc_bcpho_col[1, 4] ≈ 0.003 / snowmass_4

        # Dust total in layer 4
        @test aer.mss_dsttot_col[1, 4] ≈ 0.03  # 0.01 + 0.02

        # Inactive layers (1,2) should be zeroed
        @test aer.mss_bcpho_col[1, 1] == 0.0
        @test aer.mss_bcpho_col[1, 2] == 0.0
        @test snw_rds[1, 1] == 0.0
        @test snw_rds[1, 2] == 0.0

        # Column 3: OC in top (only) layer
        @test aer.mss_octot_col[3, 5] ≈ 0.3  # 0.1 + 0.2
        @test aer.mss_oc_col_col[3] ≈ 0.3
        @test aer.mss_oc_top_col[3] ≈ 0.3

        # Off-snow columns (2,4) should have zero masses
        @test aer.mss_bc_col_col[2] == 0.0
        @test aer.mss_bc_top_col[2] == 0.0
        @test all(aer.mss_bcpho_col[2, :] .== 0.0)
        @test all(aer.mss_dsttot_col[4, :] .== 0.0)
    end

    @testset "aerosol_fluxes! deposition" begin
        aer = CLM.AerosolData()
        CLM.aerosol_init!(aer, nc)
        CLM.aerosol_init_cold!(aer, 1:nc)

        snl = [-2, 0, -1, 0]
        col_gridcell = [1, 1, 2, 2]
        dtime = 1800.0  # 30-minute timestep

        # forc_aer: (ng, 14) aerosol deposition forcing
        forc_aer = zeros(ng, 14)
        # BC: species 1 (BCPHI dry), 2 (BCPHO dry), 3 (BCPHI wet)
        forc_aer[1, 1] = 1.0e-10  # BCPHI dry, gridcell 1
        forc_aer[1, 2] = 2.0e-10  # BCPHO dry, gridcell 1
        forc_aer[1, 3] = 3.0e-10  # BCPHI wet, gridcell 1
        # Dust species 1: wet=7, dry=8
        forc_aer[2, 7] = 5.0e-10  # DST1 wet, gridcell 2
        forc_aer[2, 8] = 6.0e-10  # DST1 dry, gridcell 2

        mask_snow = BitVector([true, false, true, false])

        CLM.aerosol_fluxes!(aer, mask_snow, 1:nc, snl, col_gridcell, forc_aer, dtime)

        # Check deposition fluxes for column 1 (gridcell 1)
        @test aer.flx_bc_dep_dry_col[1] ≈ 3.0e-10  # aer1 + aer2
        @test aer.flx_bc_dep_wet_col[1] ≈ 3.0e-10  # aer3
        @test aer.flx_bc_dep_phi_col[1] ≈ 4.0e-10  # aer1 + aer3
        @test aer.flx_bc_dep_pho_col[1] ≈ 2.0e-10  # aer2
        @test aer.flx_bc_dep_col[1]     ≈ 6.0e-10  # sum of 1,2,3

        # Check deposition into top snow layer for column 1
        # snl=-2, top_j = -2 + 5 + 1 = 4
        top_j_1 = snl[1] + nlevsno + 1  # = 4
        @test aer.mss_bcphi_col[1, top_j_1] ≈ forc_aer[1, 1] * dtime + forc_aer[1, 3] * dtime  # phi dep
        @test aer.mss_bcpho_col[1, top_j_1] ≈ forc_aer[1, 2] * dtime  # pho dep

        # Check column 3 (gridcell 2, snl=-1, top_j = -1+5+1 = 5)
        top_j_3 = snl[3] + nlevsno + 1  # = 5
        @test aer.mss_dst1_col[3, top_j_3] ≈ (forc_aer[2, 7] + forc_aer[2, 8]) * dtime

        # Non-snow columns should have deposition fluxes set but no mass added
        @test aer.flx_dst_dep_col[4] ≈ forc_aer[2, 7] + forc_aer[2, 8]
        # Mass should remain 0 for non-snow columns
        @test all(aer.mss_bcphi_col[2, :] .== 0.0)
    end

    @testset "aerosol_fluxes! with snicar_use_aerosol=false" begin
        aer = CLM.AerosolData()
        CLM.aerosol_init!(aer, nc)
        CLM.aerosol_init_cold!(aer, 1:nc)

        snl = [-1, 0, 0, 0]
        col_gridcell = [1, 1, 1, 1]
        dtime = 1800.0

        forc_aer = ones(ng, 14) * 1.0e-9
        mask_snow = BitVector([true, false, false, false])

        CLM.aerosol_fluxes!(aer, mask_snow, 1:nc, snl, col_gridcell, forc_aer, dtime;
                            snicar_use_aerosol=false)

        # All fluxes should be zero
        @test aer.flx_bc_dep_col[1] == 0.0
        @test aer.flx_oc_dep_col[1] == 0.0
        @test aer.flx_dst_dep_col[1] == 0.0

        # No mass should be deposited
        @test all(aer.mss_bcphi_col[1, :] .== 0.0)
    end

    @testset "aerosol_restart! concentration computation" begin
        aer = CLM.AerosolData()
        CLM.aerosol_init!(aer, nc)

        # Set aerosol masses
        aer.mss_bcpho_col[1, 3] = 0.01
        aer.mss_bcphi_col[1, 3] = 0.02
        aer.mss_ocpho_col[1, 3] = 0.005
        aer.mss_ocphi_col[1, 3] = 0.015
        aer.mss_dst1_col[1, 3]  = 0.1

        # Set h2osoi with water in layer 3 and zero in layer 1
        h2osoi_ice = zeros(nc, nlevsno)
        h2osoi_liq = zeros(nc, nlevsno)
        h2osoi_ice[1, 3] = 8.0
        h2osoi_liq[1, 3] = 2.0  # total = 10.0

        CLM.aerosol_restart!(aer, 1:nc; flag="read",
                             h2osoi_ice_col=h2osoi_ice,
                             h2osoi_liq_col=h2osoi_liq)

        @test aer.mss_cnc_bcpho_col[1, 3] ≈ 0.01 / 10.0
        @test aer.mss_cnc_bcphi_col[1, 3] ≈ 0.02 / 10.0
        @test aer.mss_cnc_ocpho_col[1, 3] ≈ 0.005 / 10.0
        @test aer.mss_cnc_dst1_col[1, 3]  ≈ 0.1 / 10.0

        # Layer with zero water should have zero concentrations
        @test aer.mss_cnc_bcpho_col[1, 1] == 0.0
    end

    @testset "aerosol_init_history! stub" begin
        aer = CLM.AerosolData()
        CLM.aerosol_init!(aer, nc)
        @test CLM.aerosol_init_history!(aer, 1:nc) === nothing
    end

    # Restore
    vp.nlevsno = saved_nlevsno
end
