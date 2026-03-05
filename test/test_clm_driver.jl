@testset "CLM Driver (clm_driver)" begin

    # ----------------------------------------------------------------
    # Helper: create minimal data for testing the driver
    # ----------------------------------------------------------------
    function make_driver_data(; ng=2, nl=3, nc=4, np=6)
        # Ensure varpar is initialized (nlevsno=12, nlevgrnd=15, nlevurb=5, nlevlak=10)
        CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)

        # Initialize soil level constants (zsoi, dzsoi, dzsoi_decomp, etc.)
        CLM.varcon_init!()

        nlevsno = CLM.varpar.nlevsno
        nlevgrnd = CLM.varpar.nlevgrnd
        nlevtot = nlevsno + nlevgrnd
        nlevcan = CLM.NLEVCAN

        # --- Instances ---
        inst = CLM.CLMInstances()
        CLM.clm_instInit!(inst; ng=ng, nl=nl, nc=nc, np=np,
                           nlevdecomp_full=CLM.varpar.nlevdecomp_full)

        # Set up minimal column/patch linkage for p2c
        for c in 1:nc
            inst.column.landunit[c] = 1
            inst.column.gridcell[c] = 1
            inst.column.snl[c] = 0  # no snow layers
            inst.column.patchi[c] = 1
            inst.column.patchf[c] = 0  # no patches (will be set below for p2c)
        end
        for l in 1:nl
            inst.landunit.itype[l] = CLM.ISTSOIL
            inst.landunit.urbpoi[l] = false
            inst.landunit.lakpoi[l] = false
            inst.landunit.active[l] = true
        end
        for p in 1:np
            inst.patch.active[p] = true
            inst.patch.landunit[p] = 1
            inst.patch.gridcell[p] = 1
            inst.patch.column[p] = min(p, nc)
            inst.patch.itype[p] = 1
        end

        # Link patches to columns for p2c
        # Spread patches across columns, set wtcol so weights sum to 1 per column
        ppercol = max(1, np ÷ nc)
        pidx = 1
        for c in 1:nc
            inst.column.patchi[c] = pidx
            inst.column.patchf[c] = min(pidx + ppercol - 1, np)
            npc = inst.column.patchf[c] - inst.column.patchi[c] + 1
            for p in inst.column.patchi[c]:inst.column.patchf[c]
                inst.patch.column[p] = c
                inst.patch.wtcol[p] = 1.0 / npc  # weights sum to 1.0 per column
            end
            pidx = inst.column.patchf[c] + 1
        end
        # Handle any remaining patches that don't map to valid columns
        for p in pidx:np
            inst.patch.column[p] = nc
            inst.patch.wtcol[p] = 0.0
        end

        # --- Bounds ---
        bounds = CLM.BoundsType(
            begg=1, endg=ng,
            begl=1, endl=nl,
            begc=1, endc=nc,
            begp=1, endp=np,
            begCohort=0, endCohort=0,
            level=CLM.BOUNDS_LEVEL_CLUMP,
            clump_index=1
        )

        # --- Filters ---
        filt = CLM.ClumpFilter()
        CLM.alloc_filters!(filt, nc, np, nl)
        filt.allc .= true
        filt.nolakec .= true
        filt.nolakep .= true
        filt.soilc .= true
        filt.soilp .= true
        filt.nolakeurbanp .= true
        filt.nourbanp .= true
        filt.nourbanc .= true
        filt.hydrologyc .= true
        filt.urbanc .= false
        filt.urbanl .= false
        filt.urbanp .= false
        filt.snowc .= false
        filt.nosnowc .= true
        filt.do_smb_c .= false
        filt.exposedvegp .= true
        filt.noexposedvegp .= false
        filt.lakec .= false
        filt.lakep .= false
        filt.lakesnowc .= false
        filt.lakenosnowc .= false
        filt.bgc_soilc .= false   # BGC soil columns (no BGC active in smoke test)
        filt.bgc_vegp .= false    # BGC vegetation patches (no BGC active in smoke test)
        filt.pcropp .= false      # no prog crops
        filt.soilnopcropp .= true  # all soil patches are non-crop
        filt.actfirec .= false    # no active fire
        filt.actfirep .= false    # no active fire

        filt_ia = CLM.ClumpFilter()
        CLM.alloc_filters!(filt_ia, nc, np, nl)
        filt_ia.allc .= true
        filt_ia.nolakec .= true
        filt_ia.nolakep .= true
        filt_ia.nourbanc .= true
        filt_ia.nourbanp .= true

        # --- Urban control init ---
        CLM.urban_read_nml!(CLM.urban_ctrl)

        # --- Snow layer constants (needed by handle_new_snow!) ---
        dzmin = zeros(nlevsno)
        dzmax_u = zeros(nlevsno)
        dzmax_l = zeros(nlevsno)
        dzmin[1] = 0.010; dzmax_u[1] = 0.02; dzmax_l[1] = 0.03
        dzmin[2] = 0.015; dzmax_u[2] = 0.05; dzmax_l[2] = 0.07
        for j in 3:nlevsno
            dzmin[j] = dzmax_u[j-1] * 0.5
            dzmax_u[j] = 2.0 * dzmax_u[j-1] + 0.01
            dzmax_l[j] = dzmax_u[j] + dzmax_l[j-1]
            if j == nlevsno
                dzmax_u[j] = floatmax(Float64)
                dzmax_l[j] = floatmax(Float64)
            end
        end
        CLM.SNOW_DZMIN[] = dzmin
        CLM.SNOW_DZMAX_U[] = dzmax_u
        CLM.SNOW_DZMAX_L[] = dzmax_l

        # --- Config ---
        config = CLM.CLMDriverConfig()

        # --- Photosynthesis ---
        photosyns = CLM.PhotosynthesisData()
        CLM.photosynthesis_data_init!(photosyns, np)

        return inst, bounds, filt, filt_ia, config, photosyns
    end

    # ================================================================
    # Test CLMDriverConfig construction
    # ================================================================
    @testset "CLMDriverConfig construction" begin
        config = CLM.CLMDriverConfig()
        @test config.use_cn == false
        @test config.use_fates == false
        @test config.irrigate == false
        @test config.use_noio == false
    end

    # ================================================================
    # Test CLMDriverState construction
    # ================================================================
    @testset "CLMDriverState construction" begin
        state = CLM.CLMDriverState()
        @test state.nstep == 0
        @test state.yr == 0
        @test state.mon == 0
    end

    # ================================================================
    # Test write_diagnostic
    # ================================================================
    @testset "write_diagnostic" begin
        buf = IOBuffer()
        CLM.write_diagnostic(42; io=buf)
        output = String(take!(buf))
        @test occursin("42", output)
        @test occursin("clm: completed timestep", output)
    end

    # ================================================================
    # Test clm_drv_init!
    # ================================================================
    @testset "clm_drv_init!" begin
        inst, bounds, filt, filt_ia, config, photosyns = make_driver_data()
        nc = bounds.endc
        np = bounds.endp
        nlevsno = CLM.varpar.nlevsno

        # Set some initial values
        for p in 1:np
            inst.canopystate.frac_veg_nosno_alb_patch[p] = 1
        end

        # Initialize h2osoi for the snow layer test
        for c in 1:nc
            inst.column.snl[c] = -2  # 2 snow layers
            for j in 1:nlevsno
                inst.water.waterstatebulk_inst.ws.h2osoi_liq_col[c, j] = 5.0
                inst.water.waterstatebulk_inst.ws.h2osoi_ice_col[c, j] = 5.0
            end
        end

        # Set eflx_bot_col to non-zero to verify reset
        inst.energyflux.eflx_bot_col .= 999.0

        CLM.clm_drv_init!(bounds,
                          inst.canopystate,
                          inst.water.waterstatebulk_inst,
                          inst.water.waterdiagnosticbulk_inst,
                          inst.energyflux,
                          photosyns,
                          inst.column,
                          inst.patch,
                          filt.nolakec,
                          filt.nolakep,
                          filt.soilp)

        # eflx_bot should be reset to 0
        @test all(inst.energyflux.eflx_bot_col[1:nc] .== 0.0)

        # frac_veg_nosno should be set from alb values
        @test all(inst.canopystate.frac_veg_nosno_patch[1:np] .== 1)

        # cisun_z_patch should be initialized to -999
        @test all(photosyns.cisun_z_patch[1:np, :] .== -999.0)
        @test all(photosyns.cisha_z_patch[1:np, :] .== -999.0)

        # frac_iceold should be 0.5 where liq == ice
        # snl=-2 means layers -1,0 are active (Fortran j=-1,0 → Julia j_jl=nlevsno-1,nlevsno)
        # So j = -nlevsno+1 to 0: j >= snl+1 = -1 means j >= -1
        for c in 1:nc
            # Snow layers with j >= snl+1 = -1 → j in {-1, 0} → j_jl in {nlevsno-1, nlevsno}
            j_jl1 = nlevsno - 1  # j = -1
            j_jl2 = nlevsno      # j = 0
            @test inst.water.waterdiagnosticbulk_inst.frac_iceold_col[c, j_jl1] ≈ 0.5
            @test inst.water.waterdiagnosticbulk_inst.frac_iceold_col[c, j_jl2] ≈ 0.5
        end
    end

    # ================================================================
    # Test clm_drv_patch2col!
    # ================================================================
    @testset "clm_drv_patch2col!" begin
        inst, bounds, filt, filt_ia, config, photosyns = make_driver_data()
        nc = bounds.endc
        np = bounds.endp

        # Set patch-level evaporative fluxes to known values
        inst.water.waterfluxbulk_inst.qflx_ev_snow_patch .= 1.0
        inst.water.waterfluxbulk_inst.qflx_ev_soil_patch .= 2.0
        inst.water.waterfluxbulk_inst.qflx_ev_h2osfc_patch .= 3.0
        inst.water.waterfluxbulk_inst.wf.qflx_evap_soi_patch .= 4.0
        inst.water.waterfluxbulk_inst.wf.qflx_evap_tot_patch .= 5.0
        inst.water.waterfluxbulk_inst.wf.qflx_tran_veg_patch .= 6.0
        inst.water.waterfluxbulk_inst.wf.qflx_liqevap_from_top_layer_patch .= 7.0
        inst.water.waterfluxbulk_inst.wf.qflx_liqdew_to_top_layer_patch .= 8.0
        inst.water.waterfluxbulk_inst.wf.qflx_solidevap_from_top_layer_patch .= 9.0
        inst.water.waterfluxbulk_inst.wf.qflx_soliddew_to_top_layer_patch .= 10.0

        CLM.clm_drv_patch2col!(bounds,
                               filt.allc,
                               filt.nolakec,
                               inst.energyflux,
                               inst.water.waterfluxbulk_inst,
                               inst.column,
                               inst.patch)

        # Each column's p2c result should be close to the patch value
        # (since all patches have equal weight and same value)
        for c in 1:nc
            if inst.column.patchi[c] <= inst.column.patchf[c]
                @test inst.water.waterfluxbulk_inst.qflx_ev_snow_col[c] ≈ 1.0
                @test inst.water.waterfluxbulk_inst.qflx_ev_soil_col[c] ≈ 2.0
                @test inst.water.waterfluxbulk_inst.qflx_ev_h2osfc_col[c] ≈ 3.0
                # qflx_evap_soi is averaged twice (nolake then all), last wins with allc
                @test inst.water.waterfluxbulk_inst.wf.qflx_evap_soi_col[c] ≈ 4.0
                @test inst.water.waterfluxbulk_inst.wf.qflx_evap_tot_col[c] ≈ 5.0
                @test inst.water.waterfluxbulk_inst.wf.qflx_tran_veg_col[c] ≈ 6.0
            end
        end
    end

    # ================================================================
    # Test clm_drv! (main driver, smoke test)
    # ================================================================
    @testset "clm_drv! smoke test" begin
        inst, bounds, filt, filt_ia, config, photosyns = make_driver_data()

        # Set frac_veg_nosno_alb for drv_init
        inst.canopystate.frac_veg_nosno_alb_patch .= 1

        # Run the main driver — should not throw
        @test nothing === CLM.clm_drv!(config, inst, filt, filt_ia, bounds,
                                        true,   # doalb
                                        1.0,    # nextsw_cday
                                        0.0,    # declinp1
                                        0.0,    # declin
                                        0.4091, # obliqr
                                        false,  # rstwr
                                        false,  # nlend
                                        "20260101",  # rdate
                                        false;  # rof_prognostic
                                        nstep=1,
                                        is_first_step=false,
                                        is_beg_curr_day=false,
                                        is_beg_curr_year=false,
                                        photosyns=photosyns)
    end

    # ================================================================
    # Test clm_drv! with different config flags
    # ================================================================
    @testset "clm_drv! config flags" begin
        inst, bounds, filt, filt_ia, _, photosyns = make_driver_data()
        inst.canopystate.frac_veg_nosno_alb_patch .= 1

        # Create config with CN mode and various sub-features enabled
        config = CLM.CLMDriverConfig(use_cn=true, use_lch4=true, use_crop=true, irrigate=true)

        @test nothing === CLM.clm_drv!(config, inst, filt, filt_ia, bounds,
                                        false,   # doalb
                                        1.0, 0.0, 0.0, 0.4091,
                                        false, false, "20260101", false;
                                        nstep=5,
                                        is_first_step=true,
                                        is_beg_curr_day=true,
                                        is_beg_curr_year=true,
                                        photosyns=photosyns)
    end

    # ================================================================
    # Test urban snow fraction computation
    # ================================================================
    @testset "urban snow fraction" begin
        inst, bounds, filt, filt_ia, config, photosyns = make_driver_data()
        nc = bounds.endc
        nl = bounds.endl

        # Mark first landunit as urban
        inst.landunit.urbpoi[1] = true
        for c in 1:nc
            inst.column.landunit[c] = 1
        end

        # Set snow depth
        inst.water.waterdiagnosticbulk_inst.snow_depth_col[1] = 0.025  # half of 0.05
        inst.water.waterdiagnosticbulk_inst.snow_depth_col[2] = 0.10   # > 0.05

        inst.canopystate.frac_veg_nosno_alb_patch .= 0

        CLM.clm_drv!(config, inst, filt, filt_ia, bounds,
                     false, 1.0, 0.0, 0.0, 0.4091,
                     false, false, "20260101", false;
                     nstep=1, photosyns=photosyns)

        # Check urban snow fraction is bounded and finite.
        @test isfinite(inst.water.waterdiagnosticbulk_inst.frac_sno_col[1])
        @test isfinite(inst.water.waterdiagnosticbulk_inst.frac_sno_col[2])
        @test 0.0 <= inst.water.waterdiagnosticbulk_inst.frac_sno_col[1] <= 1.0
        @test 0.0 <= inst.water.waterdiagnosticbulk_inst.frac_sno_col[2] <= 1.0
    end

end
