@testset "Glacier Surface Mass Balance" begin
    # ------------------------------------------------------------------
    # Tests for GlacierSurfaceMassBalanceMod port:
    #   1. handle_ice_melt!              — liquid→ice conversion + melt flux
    #   2. compute_surface_mass_balance! — frz / net glcice / dyn water flux
    #   3. adjust_runoff_terms!          — runoff adjustments (incl. dyn routing)
    #   4. mass conservation on the istice path
    # ------------------------------------------------------------------

    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
    nlevsno  = CLM.varpar.nlevsno
    nlevgrnd = CLM.varpar.nlevgrnd
    ISTICE   = CLM.ISTICE
    SECSPDAY = CLM.SECSPDAY

    # ------------------------------------------------------------------
    # 1. handle_ice_melt! : meltwater converted back to ice, flux accumulated
    # ------------------------------------------------------------------
    @testset "handle_ice_melt!" begin
        nc = 2
        bounds = 1:nc
        dtime = 1800.0
        ncol_lyr = nlevsno + nlevgrnd

        # Column 1: istice with meltwater; column 2: non-istice (ISTSOIL) with liquid.
        col_landunit = [1, 2]
        lun_itype    = [ISTICE, CLM.ISTSOIL]
        mask_do_smb  = BitVector([true, true])

        h2osoi_liq = zeros(nc, ncol_lyr)
        h2osoi_ice = zeros(nc, ncol_lyr)
        # Ground layers 1 and 2 of column 1 have meltwater.
        j1 = 1 + nlevsno
        j2 = 2 + nlevsno
        h2osoi_liq[1, j1] = 5.0
        h2osoi_ice[1, j1] = 100.0
        h2osoi_liq[1, j2] = 3.0
        h2osoi_ice[1, j2] = 50.0
        # Column 2 (non-istice) liquid — must be left untouched.
        h2osoi_liq[2, j1] = 7.0
        h2osoi_ice[2, j1] = 0.0

        qflx_glcice_melt = fill(NaN, nc)

        # Conserve total mass (liq+ice) per layer before the call.
        tot_before_c1_j1 = h2osoi_liq[1, j1] + h2osoi_ice[1, j1]
        tot_before_c1_j2 = h2osoi_liq[1, j2] + h2osoi_ice[1, j2]

        CLM.handle_ice_melt!(h2osoi_liq, h2osoi_ice, qflx_glcice_melt,
            col_landunit, lun_itype, mask_do_smb, bounds, dtime,
            nlevsno, nlevgrnd)

        # Melt flux = sum of melted liquid / dtime (positive definite).
        @test qflx_glcice_melt[1] ≈ (5.0 + 3.0) / dtime
        # Non-istice column: no melt flux generated.
        @test qflx_glcice_melt[2] == 0.0

        # istice layers converted to pure ice; total layer mass conserved.
        @test h2osoi_liq[1, j1] == 0.0
        @test h2osoi_liq[1, j2] == 0.0
        @test h2osoi_ice[1, j1] ≈ tot_before_c1_j1
        @test h2osoi_ice[1, j2] ≈ tot_before_c1_j2

        # Non-istice column untouched.
        @test h2osoi_liq[2, j1] == 7.0
        @test h2osoi_ice[2, j1] == 0.0
    end

    # ------------------------------------------------------------------
    # 2. compute_surface_mass_balance!
    # ------------------------------------------------------------------
    @testset "compute_surface_mass_balance!" begin
        nc = 3
        bounds = 1:nc

        # col 1: istice (always accumulates frz)
        # col 2: non-istice, snow persisted long enough (glacial inception)
        # col 3: non-istice, snow not persistent enough (no frz)
        col_landunit = [1, 2, 3]
        col_gridcell = [1, 1, 1]
        lun_itype    = [ISTICE, CLM.ISTSOIL, CLM.ISTSOIL]
        mask_allc    = BitVector([true, true, true])
        mask_do_smb  = BitVector([true, true, true])

        max_days = 7300
        threshold = max_days * SECSPDAY
        snow_persistence = [0.0, threshold + 1.0, threshold - 1.0]

        qflx_snwcp_ice   = [0.4, 0.5, 0.6]
        qflx_glcice_melt = [0.1, 0.0, 0.0]
        qflx_glcice      = fill(NaN, nc)
        qflx_glcice_frz  = fill(NaN, nc)
        qflx_dyn         = fill(NaN, nc)
        glc_dyn_routing  = [0.0]

        CLM.compute_surface_mass_balance!(
            qflx_glcice, qflx_glcice_frz, qflx_dyn,
            qflx_snwcp_ice, qflx_glcice_melt, snow_persistence,
            glc_dyn_routing, col_landunit, col_gridcell, lun_itype,
            mask_allc, mask_do_smb, bounds;
            glc_snow_persistence_max_days = max_days)

        # col 1 (istice) and col 2 (persistent snow): frz = snwcp_ice
        @test qflx_glcice_frz[1] ≈ 0.4
        @test qflx_glcice_frz[2] ≈ 0.5
        # col 3: not persistent, not istice → frz = 0
        @test qflx_glcice_frz[3] == 0.0

        # net glcice = frz - melt
        @test qflx_glcice[1] ≈ 0.4 - 0.1
        @test qflx_glcice[2] ≈ 0.5
        @test qflx_glcice[3] ≈ 0.0

        # routing = 0 → dyn water flux = 0
        @test all(qflx_dyn .== 0.0)
    end

    # ------------------------------------------------------------------
    # 3. compute_surface_mass_balance! + adjust_runoff_terms! with dyn routing > 0
    #    Mass-balance consistency check across the SMB/runoff term set.
    # ------------------------------------------------------------------
    @testset "dynamic routing + runoff mass balance" begin
        nc = 1
        bounds = 1:nc
        col_landunit = [1]
        col_gridcell = [1]
        lun_itype    = [ISTICE]
        mask_allc    = BitVector([true])
        mask_do_smb  = BitVector([true])

        snow_persistence = [0.0]
        qflx_snwcp_ice   = [1.0]
        qflx_glcice_melt = [0.3]
        qflx_glcice      = [NaN]
        qflx_glcice_frz  = [NaN]
        qflx_dyn         = [NaN]
        routing          = 0.7
        glc_dyn_routing  = [routing]

        CLM.compute_surface_mass_balance!(
            qflx_glcice, qflx_glcice_frz, qflx_dyn,
            qflx_snwcp_ice, qflx_glcice_melt, snow_persistence,
            glc_dyn_routing, col_landunit, col_gridcell, lun_itype,
            mask_allc, mask_do_smb, bounds)

        @test qflx_glcice_frz[1] ≈ 1.0                       # istice → frz = snwcp
        @test qflx_glcice[1] ≈ 1.0 - 0.3                     # frz - melt
        @test qflx_dyn[1] ≈ routing * (0.3 - 1.0)            # routing*(melt-frz)

        qflx_qrgwl = [0.0]
        qflx_ice_runoff_snwcp = [2.0]   # initial snow-capping ice runoff

        CLM.adjust_runoff_terms!(qflx_qrgwl, qflx_ice_runoff_snwcp,
            qflx_glcice_frz, qflx_glcice_melt, glc_dyn_routing,
            col_gridcell, mask_do_smb, bounds)

        # qrgwl gets all the melt as liquid runoff.
        @test qflx_qrgwl[1] ≈ 0.3
        # ice runoff: -= routing*frz  -= (1-routing)*melt
        @test qflx_ice_runoff_snwcp[1] ≈ 2.0 - routing*1.0 - (1.0 - routing)*0.3

        # All outputs finite.
        @test all(isfinite, qflx_glcice)
        @test all(isfinite, qflx_glcice_frz)
        @test all(isfinite, qflx_dyn)
        @test all(isfinite, qflx_qrgwl)
        @test all(isfinite, qflx_ice_runoff_snwcp)
    end

    # ------------------------------------------------------------------
    # 4. Closed system mass conservation on the istice path
    #    (handle_ice_melt borrows ice from below; the runoff/dyn terms must
    #     reconcile the borrowed mass exactly).
    # ------------------------------------------------------------------
    @testset "istice path mass conservation" begin
        nc = 1
        bounds = 1:nc
        dtime = 1800.0
        ncol_lyr = nlevsno + nlevgrnd

        col_landunit = [1]
        col_gridcell = [1]
        lun_itype    = [ISTICE]
        mask_allc    = BitVector([true])
        mask_do_smb  = BitVector([true])

        h2osoi_liq = zeros(nc, ncol_lyr)
        h2osoi_ice = zeros(nc, ncol_lyr)
        j1 = 1 + nlevsno
        melt_mass = 4.0
        h2osoi_liq[1, j1] = melt_mass
        h2osoi_ice[1, j1] = 200.0

        qflx_glcice_melt = [NaN]
        CLM.handle_ice_melt!(h2osoi_liq, h2osoi_ice, qflx_glcice_melt,
            col_landunit, lun_itype, mask_do_smb, bounds, dtime,
            nlevsno, nlevgrnd)

        # Melt flux (mm/s) integrated over dtime = the meltwater mass removed.
        @test qflx_glcice_melt[1] * dtime ≈ melt_mass

        # No dynamic coupling → the borrow (source) and the snwcp reduction (sink)
        # cancel exactly: net water added to the system via the dyn term is 0.
        snow_persistence = [0.0]
        qflx_snwcp_ice   = [1.0]
        qflx_glcice      = [NaN]
        qflx_glcice_frz  = [NaN]
        qflx_dyn         = [NaN]
        glc_dyn_routing  = [0.0]

        CLM.compute_surface_mass_balance!(
            qflx_glcice, qflx_glcice_frz, qflx_dyn,
            qflx_snwcp_ice, qflx_glcice_melt, snow_persistence,
            glc_dyn_routing, col_landunit, col_gridcell, lun_itype,
            mask_allc, mask_do_smb, bounds)

        @test qflx_dyn[1] == 0.0   # routing=0 → no dyn balance correction
    end
end
