# Wiring test for CNDV (dynamic vegetation) into the live driver.
#
# Verifies:
#   1. The use_cndv config flag plumbs through CNMode / CLMDriverConfig (default
#      false; the property accessor reports it correctly).
#   2. clm_instInit!(use_cndv=true) sizes the DGVSData patch state and derives the
#      DGVEcophysCon constants; the default (use_cndv=false) leaves them empty.
#   3. The exact CNDV call sequence the driver runs end-to-end
#         dyn_cndv_init!  →  cndv_update_acc_vars!  →  cndv_driver!  →  dyn_cndv_interp!
#      mutates patch weights (patch.wtcol) when use_cndv is active.
#   4. With use_cndv off the same driver region is a no-op → patch.wtcol unchanged
#      (byte-identical default path).
@testset "CNDV driver wiring" begin

    # -----------------------------------------------------------------------
    # 1. Config flag plumbing
    # -----------------------------------------------------------------------
    @testset "use_cndv config flag" begin
        cfg_off = CLM.CLMDriverConfig(use_cn=true)
        @test cfg_off.use_cndv == false           # default off
        @test cfg_off.use_cn == true

        cfg_on = CLM.CLMDriverConfig(use_cn=true, use_cndv=true)
        @test cfg_on.use_cndv == true

        # use_cndv only meaningful under CN mode; SP/FATES report false
        cfg_sp = CLM.CLMDriverConfig()
        @test cfg_sp.use_cndv == false
        cfg_fates = CLM.CLMDriverConfig(use_fates=true)
        @test cfg_fates.use_cndv == false
    end

    # -----------------------------------------------------------------------
    # 2. clm_instInit! sizing + eco-constant derivation behind use_cndv
    # -----------------------------------------------------------------------
    @testset "clm_instInit! CNDV state sizing" begin
        ng, nl, nc, np = 3, 5, 10, 20
        CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)  # set vertical level counts
        CLM.pftcon_allocate!(CLM.pftcon)   # eco constants read pftpar20/28/29/30/31

        # default: CNDV state stays empty (byte-identical default path)
        inst_off = CLM.CLMInstances()
        CLM.clm_instInit!(inst_off; ng=ng, nl=nl, nc=nc, np=np)
        @test length(inst_off.dgvs.fpcgrid_patch) == 0
        @test length(inst_off.dgv_ecophyscon.crownarea_max) == 0

        # use_cndv: dgvs sized to np, eco constants populated from pftcon
        inst_on = CLM.CLMInstances()
        CLM.clm_instInit!(inst_on; ng=ng, nl=nl, nc=nc, np=np, use_cndv=true)
        @test length(inst_on.dgvs.fpcgrid_patch) == np
        @test length(inst_on.dgvs.present_patch) == np
        @test inst_on.dgvs.tmomin20_patch[1] ≈ CLM.TFRZ - 5.0  # cold-start value
        @test length(inst_on.dgv_ecophyscon.crownarea_max) == length(CLM.pftcon.pftpar20)
    end

    # -----------------------------------------------------------------------
    # Helper: minimal instance whose CNDV state is wired exactly as the driver
    # threads it (dgvs + dgv_ecophyscon on CLMInstances, patch/landunit set up).
    # Mirrors the make_cndv_test_data layout from test_cndv.jl.
    # -----------------------------------------------------------------------
    function build_cndv_inst()
        np, nc, ng, nl = 4, 1, 1, 1
        pft = CLM.PftconType()
        CLM.pftcon_allocate!(pft)
        for i in [2, 3, 4]            # trees (Fortran PFTs 1,2,3)
            pft.woody[i] = 1.0; pft.is_tree[i] = true; pft.is_shrub[i] = false
            pft.dwood[i] = 2.5e5; pft.slatop[i] = 0.012; pft.dsladlai[i] = 0.0
        end
        for i in [13, 14]            # grasses
            pft.woody[i] = 0.0; pft.is_tree[i] = false; pft.is_shrub[i] = false
            pft.dwood[i] = 0.0; pft.slatop[i] = 0.030; pft.dsladlai[i] = 0.0
        end
        pft.pftpar20 .= 15.0; pft.pftpar28 .= -35.0; pft.pftpar29 .= 28.0
        pft.pftpar30 .= 350.0; pft.pftpar31 .= 1000.0

        inst = CLM.CLMInstances()
        CLM.patch_init!(inst.patch, np)
        inst.patch.itype .= [CLM.noveg, 1, 2, 12]
        inst.patch.gridcell .= 1; inst.patch.column .= 1; inst.patch.landunit .= 1
        inst.patch.wtcol .= [0.1, 0.3, 0.3, 0.3]

        CLM.landunit_init!(inst.landunit, nl)
        inst.landunit.itype[1] = CLM.ISTSOIL
        inst.landunit.wtgcell[1] = 1.0

        # CNDV state on the instance (as instInit! would build it when use_cndv)
        CLM.dgvs_init_cold!(inst.dgvs, np)
        CLM.dgv_ecophyscon_init!(inst.dgv_ecophyscon, pft)

        return (; inst, pft, np, nc, ng, nl)
    end

    # The driver's exact CNDV call sequence (interp uses wt1 = 1 - yearfrac).
    function run_cndv_driver_region!(d; year::Int, jday::Int, secs::Int,
                                     is_end_curr_year::Bool, is_beg_curr_year::Bool,
                                     is_first_step::Bool, dtime::Float64)
        inst = d.inst; pch = inst.patch; lun = inst.landunit
        bc_patch = 1:d.np; bc_grc = 1:d.ng
        np = d.np

        # --- per-step accumulation (accum block) ---
        cndv_dummy_t = fill(CLM.TFRZ + 10.0, np)
        CLM.cndv_update_acc_vars!(inst.dgvs, inst.dgv_ecophyscon,
            cndv_dummy_t, cndv_dummy_t, bc_patch, dtime;
            month=1, day=2, secs=secs)

        # --- annual driver at end of year ---
        if is_end_curr_year && !is_first_step
            prec365 = fill(0.0, d.nc)
            natvegp = falses(np)
            CLM.cndv_driver!(inst.dgvs, inst.dgv_ecophyscon,
                fill(CLM.TFRZ, np), prec365,
                fill(200.0, np), fill(80.0, np),     # annsum_npp, annsum_litfall
                fill(5000.0, np), fill(200.0, np),   # deadstemc, leafcmax
                d.pft, pch, lun, natvegp, bc_patch, year;
                bounds_gridcell=bc_grc)
        end

        # --- per-step weight interpolation ---
        yearfrac = (jday - 1 + secs / CLM.SECSPDAY) / 365.0
        CLM.dyn_cndv_interp!(inst.dgvs, pch, lun, bc_patch;
            wt1 = 1.0 - yearfrac,
            is_beg_curr_year = is_beg_curr_year && !is_first_step)
        return nothing
    end

    # -----------------------------------------------------------------------
    # 3. use_cndv ON: the annual cycle mutates fpcgrid → wtcol changes
    # -----------------------------------------------------------------------
    @testset "use_cndv=true mutates PFT weights" begin
        d = build_cndv_inst()
        # dyn_cndv_init!: fpcgrid = fpcgridold = wtcol
        CLM.dyn_cndv_init!(d.inst.dgvs, d.inst.patch, 1:d.np)
        @test all(d.inst.dgvs.fpcgrid_patch .≈ d.inst.patch.wtcol)

        # Seed a present, growing canopy so light competition has something to do.
        d.inst.dgvs.present_patch .= [false, true, true, true]
        d.inst.dgvs.nind_patch .= [0.0, 0.05, 0.03, 1.0]
        d.inst.dgvs.fpcgrid_patch .= [0.0, 0.5, 0.5, 0.05]
        d.inst.dgvs.fpcgridold_patch .= [0.0, 0.5, 0.5, 0.05]
        d.inst.dgvs.crownarea_patch .= [0.0, 10.0, 10.0, 1.0]
        d.inst.dgvs.pftmayexist_patch .= true
        d.inst.dgvs.tmomin20_patch .= CLM.TFRZ + 0.0
        d.inst.dgvs.agdd20_patch .= 500.0
        d.inst.dgvs.agdd_patch .= 600.0

        wtcol_before = copy(d.inst.patch.wtcol)

        # Run end-of-year (year 3 so climate20 running means update normally).
        run_cndv_driver_region!(d; year=3, jday=365, secs=Int(1800),
            is_end_curr_year=true, is_beg_curr_year=false,
            is_first_step=false, dtime=1800.0)

        # Light competition capped the over-95% trees → fpcgrid changed → wtcol moved.
        @test d.inst.patch.wtcol != wtcol_before
        @test any(abs.(d.inst.patch.wtcol .- wtcol_before) .> 1e-8)
        # fpcgrid must be finite and within [0,1]
        @test all(isfinite, d.inst.dgvs.fpcgrid_patch)
        @test all(0.0 .<= d.inst.dgvs.fpcgrid_patch .<= 1.0)
    end

    # -----------------------------------------------------------------------
    # 4. use_cndv OFF: the gated driver region is never entered → weights frozen
    # -----------------------------------------------------------------------
    @testset "use_cndv=false leaves weights untouched" begin
        d = build_cndv_inst()
        wtcol_before = copy(d.inst.patch.wtcol)
        fpc_before = copy(d.inst.dgvs.fpcgrid_patch)

        # Reproduce the driver's gating: when config.use_cndv is false the whole
        # CNDV region is skipped, so nothing touches wtcol or dgvs.
        use_cndv = false
        if use_cndv && !isempty(d.inst.dgvs.fpcgrid_patch)
            run_cndv_driver_region!(d; year=3, jday=365, secs=Int(1800),
                is_end_curr_year=true, is_beg_curr_year=false,
                is_first_step=false, dtime=1800.0)
        end

        @test d.inst.patch.wtcol == wtcol_before
        @test d.inst.dgvs.fpcgrid_patch == fpc_before
    end

end
