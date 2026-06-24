# ==========================================================================
# test_clm_wiring.jl — verify the ported-but-unwired modules that were wired
# into the live driver loop (clm_drv_core!) actually run when their config flag
# is ON, and that the default path (flag OFF) is unchanged.
#
# The full clm_drv! cannot run on the minimal synthetic make_driver_data fixture
# (the forcing-dependent canopy_fluxes solve BoundsErrors there — a pre-existing
# limitation of that fixture, see test_clm_driver.jl), so each item is exercised
# through the exact integration path that clm_drv_core! now invokes:
#   * the wired entry function(s),
#   * the new sizing/setup helpers added for the wiring,
#   * the config-flag plumbing,
# asserting that the module's output state changed from its no-op baseline + is
# finite + the relevant balance/conservation holds.
# ==========================================================================
using Test, CLM
using NCDatasets

@testset "Driver wiring (methane/ozone/irrigation/dyn_subgrid)" begin

    # ======================================================================
    # 1. OZONE UPTAKE → STRESS coupling.
    # Before wiring, calc_ozone_stress! ran but o3uptake* was never updated, so
    # stress sat at coef 1.0 forever. Wiring calc_ozone_uptake! BEFORE the next
    # step's stress makes the dose accumulate and the stress coefficients drop
    # below 1.0. This asserts that exact wired sequence has teeth.
    # ======================================================================
    @testset "ozone uptake feeds stress (was a no-op before wiring)" begin
        np = CLM.MXPFT + 1
        oz = CLM.OzoneData()
        CLM.ozone_init!(oz, 2; stress_method="stress_lombardozzi2015")

        pch = CLM.PatchData()
        pch.column   = [1, 1]
        pch.gridcell = [1, 1]
        pch.itype    = [2, 5]

        evergreen_pft = zeros(np); evergreen_pft[2] = 1.0
        leaf_long_pft = fill(1.0, np); leaf_long_pft[2] = 3.0
        woody_pft     = fill(1.0, np)

        mask  = BitVector([true, true])
        nomask = BitVector([false, false])
        bounds = 1:2

        # Baseline: cold-start stress coefficients are all 1.0 (no stress).
        @test all(oz.o3coefvsun_patch .== 1.0)
        @test all(oz.o3uptakesun_patch .== 0.0)

        # Drive uptake with a large ozone dose over many steps so the accumulated
        # dose is high enough to bite (mirrors the per-timestep wired call).
        forc_pbot = [101325.0, 101325.0]
        forc_th   = [300.0, 300.0]
        rssun     = [200.0, 150.0]
        rssha     = [250.0, 180.0]
        rb_vec    = [50.0, 60.0]
        ram_vec   = [100.0, 110.0]
        tlai_vec  = [3.0, 4.0]
        forc_o3   = [200.0e-9]                 # high ambient ozone
        for _ in 1:500
            CLM.calc_ozone_uptake!(oz, pch, mask, bounds,
                forc_pbot, forc_th, rssun, rssha, rb_vec, ram_vec,
                tlai_vec, forc_o3, evergreen_pft, leaf_long_pft, 1800.0)
        end

        # MODULE RAN: uptake accumulated a strictly positive dose + tlai_old updated.
        @test oz.o3uptakesun_patch[1] > 0.0
        @test oz.o3uptakesha_patch[1] > 0.0
        @test all(isfinite, oz.o3uptakesun_patch)
        @test oz.tlai_old_patch[1] == 3.0

        # Now the stress step (the already-wired call) uses that dose → coef < 1.0.
        CLM.calc_ozone_stress!(oz, mask, nomask, bounds, pch, woody_pft;
                               is_time_to_run_luna=false)
        @test oz.o3coefvsun_patch[1] < 1.0   # photosynthesis now stressed
        @test oz.o3coefgsun_patch[1] <= 1.0
        @test all(isfinite, oz.o3coefvsun_patch)
    end

    # ======================================================================
    # 2. METHANE: ch4_init_allocate! sizes the prognostic state, the gridcell +
    # column balance-check inits seed the conservation accumulators, and ch4!
    # produces a finite surface flux + total-column CH4. This is exactly the
    # chain clm_drv_core! now runs under use_lch4.
    # ======================================================================
    @testset "methane ch4! runs + balance check seeded" begin
        nc, np, ng, nlevsoi = 3, 4, 2, 5
        ch4 = CLM.CH4Data{Float64}()
        CLM.ch4_init_allocate!(ch4, nc, np, ng, nlevsoi)

        # MODULE SIZED: every array now has its proper non-zero length.
        @test size(ch4.ch4_prod_depth_sat_col) == (nc, nlevsoi)
        @test length(ch4.ch4_surf_flux_tot_col) == nc
        @test size(ch4.c_atm_grc) == (ng, 3)
        @test length(ch4.ch4_first_time_grc) == ng
        @test all(ch4.ch4_first_time_grc)

        params = CLM.CH4Params()
        ch4vc  = CLM.CH4VarCon()
        mask_soil  = trues(nc); mask_soilp = trues(np)
        mask_nolake = trues(nc); mask_lake = falses(nc)
        watsat = fill(0.45, nc, nlevsoi)
        dz = fill(0.1, nc, nlevsoi)
        z  = zeros(nc, nlevsoi); zi = zeros(nc, nlevsoi + 1)
        for j in 1:nlevsoi, c in 1:nc
            z[c, j] = (j - 0.5) * 0.1; zi[c, j + 1] = j * 0.1
        end

        col_gridcell = [1, 1, 2]; col_wtgcell = [0.5, 0.5, 1.0]

        # Balance-check inits (the ~line-728/792 wirings) seed totcolch4_bef_*.
        ch4.conc_ch4_unsat_col .= 1.0e-5
        CLM.ch4_init_gridcell_balance_check!(ch4, mask_nolake, mask_lake,
            col_gridcell, col_wtgcell, dz, nlevsoi, ng, false)
        CLM.ch4_init_column_balance_check!(ch4, mask_nolake, mask_lake,
            dz, nlevsoi, false)
        @test all(isfinite, ch4.totcolch4_bef_grc)
        @test all(isfinite, ch4.totcolch4_bef_col)

        # Main ch4! driver (the ~line-1927 wiring).
        patch_column = [1, 1, 2, 3]; patch_itype = [1, 2, 1, 1]
        patch_wtcol = [0.5, 0.5, 1.0, 1.0]
        is_fates = falses(nc); latdeg = [45.0, 30.0]
        forc_pbot = fill(101325.0, ng); forc_t_grc = fill(CLM.TFRZ + 15.0, ng)
        forc_po2 = fill(101325.0 * 0.209, ng); forc_pco2 = fill(40.0, ng)
        forc_pch4 = fill(101325.0 * 1.7e-6, ng)
        h2osoi_vol = fill(0.3, nc, nlevsoi)
        h2osoi_liq = fill(0.3 * CLM.DENH2O * 0.1, nc, nlevsoi)
        h2osoi_ice = zeros(nc, nlevsoi); h2osfc_v = zeros(nc)
        bsw = fill(5.0, nc, nlevsoi); cellorg = fill(10.0, nc, nlevsoi)
        smp_l = fill(-1000.0, nc, nlevsoi)
        t_soisno = fill(CLM.TFRZ + 15.0, nc, nlevsoi)
        t_grnd = fill(CLM.TFRZ + 15.0, nc); t_h2osfc = fill(CLM.TFRZ + 15.0, nc)
        frac_h2osfc = zeros(nc); snow_depth = zeros(nc); snl = zeros(Int, nc)
        qflx_surf = zeros(nc)
        rootfr_p = fill(0.2, np, nlevsoi); rootfr_col = fill(0.2, nc, nlevsoi)
        crootfr = fill(0.2, np, nlevsoi); rootr_p = fill(0.2, np, nlevsoi)
        elai = fill(2.0, np); qflx_tran_veg = fill(1.0e-3, np)
        annsum_npp = fill(500.0, np); rr = fill(1.0e-6, np)
        somhr = fill(1.0e-6, nc); lithr = fill(1.0e-7, nc)
        hr_vr = fill(1.0e-8, nc, nlevsoi); o_scalar = ones(nc, nlevsoi)
        fphr = ones(nc, nlevsoi); pot_f_nit_vr = zeros(nc, nlevsoi)
        lake_icefrac = zeros(nc, nlevsoi); lakedepth = fill(5.0, nc)
        agnpp = fill(1.0e-5, np); bgnpp = fill(1.0e-5, np)

        ch4.ch4_surf_flux_tot_col .= 999.0   # sentinel to prove ch4! overwrites it

        CLM.ch4!(ch4, params, ch4vc,
            mask_soil, mask_soilp, mask_lake, mask_nolake,
            col_gridcell, col_wtgcell, patch_column, patch_itype, patch_wtcol,
            is_fates, latdeg,
            forc_pbot, forc_t_grc, forc_po2, forc_pco2, forc_pch4,
            watsat, h2osoi_vol, h2osoi_liq, h2osoi_ice, h2osfc_v,
            bsw, cellorg, smp_l, t_soisno, t_grnd, t_h2osfc,
            frac_h2osfc, snow_depth, snl, qflx_surf,
            rootfr_p, rootfr_col, crootfr, rootr_p, elai,
            qflx_tran_veg, annsum_npp, rr,
            somhr, lithr, hr_vr, o_scalar, fphr, pot_f_nit_vr,
            lake_icefrac, lakedepth, z, dz, zi,
            nlevsoi, 5, 1, nlevsoi, 5, 0.01, 130.0,
            ng, 0, 1800.0, true, false, false,
            agnpp, bgnpp, 365.0 * 86400.0)

        # MODULE RAN: first-time flag cleared, surface flux overwritten + finite.
        @test all(.!ch4.ch4_first_time_grc)
        @test all(ch4.ch4_surf_flux_tot_col .!= 999.0)
        @test all(isfinite, ch4.ch4_surf_flux_tot_col)
        @test all(isfinite, ch4.totcolch4_col)
        @test all(ch4.totcolch4_grc .>= 0.0)
        @test ch4.c_atm_grc[1, 1] > 0.0
    end

    # ======================================================================
    # 3. IRRIGATION demand → withdrawal → flux pipeline. Replicates the two
    # driver wirings: calc_irrigation_needed! (end-of-step demand) then
    # calc_irrigation_fluxes! (next-step withdrawal + application).
    # ======================================================================
    @testset "irrigation needed → fluxes pipeline" begin
        CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
        CLM.varcon_init!()
        nc, np, nlevsoi = 1, 1, 5
        # Build a soil column whose liquid water is well below the irrigation
        # target so a deficit (and thus a positive irrigation rate) is produced.
        irrig = CLM.IrrigationData()
        CLM.irrigation_init_allocate!(irrig, np, nc, nlevsoi)
        irrig.dtime = 1800
        irrig.irrig_nsteps_per_day = CLM.calc_irrig_nsteps_per_day(irrig.params.irrig_length, 1800)
        irrig.relsat_wilting_point_col .= 0.1
        irrig.relsat_target_col .= 0.9
        irrig.irrig_method_patch .= CLM.IRRIG_METHOD_DRIP

        col = CLM.ColumnData{Float64}()
        col.landunit = [1]; col.gridcell = [1]
        col.z  = fill(0.05, nc, nlevsoi); col.dz = fill(0.1, nc, nlevsoi)
        col.nbedrock = [nlevsoi]

        grc = CLM.GridcellData{Float64}(); grc.londeg = [0.0]

        pch = CLM.PatchData{Float64}()
        pch.column = [1]; pch.gridcell = [1]; pch.landunit = [1]
        pch.itype  = [CLM.nc3irrig]      # an irrigated crop PFT
        pch.wtcol  = [1.0]; pch.active = trues(np)

        npft = CLM.MXPFT + 1
        pftcon_irrigated = zeros(npft); pftcon_irrigated[CLM.nc3irrig] = 1.0

        elai = fill(2.0, np)
        t_soisno = fill(CLM.TFRZ + 10.0, nc, nlevsoi)   # unfrozen
        eff_porosity = fill(0.4, nc, nlevsoi)
        h2osoi_liq = fill(0.05 * CLM.DENH2O * 0.1, nc, nlevsoi)  # dry → deficit
        volr = zeros(1)
        mask_exposedveg = BitVector([true])

        # local solar time just after the irrigation start (so seconds_since the
        # (start - dtime) offset is < dtime → point_needs_check_for_irrig fires)
        lts = [irrig.params.irrig_start_time - irrig.dtime + 1]

        CLM.calc_irrigation_needed!(irrig, elai, t_soisno, eff_porosity,
            h2osoi_liq, volr, false, pftcon_irrigated, lts,
            col, grc, pch, mask_exposedveg, 1:nc, 1:np, 1:1, nlevsoi)

        # MODULE RAN (demand): a positive surface irrigation rate + steps queued.
        @test irrig.sfc_irrig_rate_patch[1] > 0.0
        @test irrig.n_irrig_steps_left_patch[1] > 0

        # Now the withdrawal/flux step (the ~line-810 wiring) consumes that demand.
        ss = CLM.SoilStateData(); CLM.soilstate_init!(ss, np, nc)
        sh = CLM.SoilHydrologyData(); CLM.soilhydrology_init!(sh, nc)
        wfb = CLM.WaterFluxBulkData(); CLM.waterfluxbulk_init!(wfb, nc, np, 1, 1)
        wfb.wf.qflx_sfc_irrig_col .= 0.0

        CLM.calc_irrigation_fluxes!(irrig, sh, ss, wfb, col, pch,
            BitVector([true]), BitVector([true]), 1:nc, 1:np, nlevsoi, 1800.0)

        # MODULE RAN (flux): the per-patch drip application flux is now positive,
        # and the column surface-irrigation flux was filled.
        @test wfb.wf.qflx_irrig_drip_patch[1] > 0.0
        @test wfb.wf.qflx_sfc_irrig_col[1] > 0.0
        @test isfinite(wfb.wf.qflx_irrig_drip_patch[1])
    end

    # ======================================================================
    # 4. DYN_SUBGRID: setup_dyn_subgrid! attaches a DynSubgridState to the config
    # (the init counterpart of the driver hooks), and dynSubgrid_driver! snaps the
    # patch weights to the per-year transient-PFT values. This is the exact state
    # object + driver call that clm_drv_core! now invokes at year boundaries.
    # ======================================================================
    @testset "dyn_subgrid setup + driver changes weights" begin
        natpft_size = 3; ng = 1; ncol = 1; np = natpft_size
        file_years = [2000, 2001]
        pct = Array{Float64}(undef, natpft_size, ng, length(file_years))
        pct[:, 1, 1] = [50.0, 30.0, 20.0]
        pct[:, 1, 2] = [10.0, 10.0, 80.0]

        mktempdir() do dir
            fn = joinpath(dir, "flanduse_pft.nc")
            NCDataset(fn, "c") do ds
                defDim(ds, "natpft", natpft_size)
                defDim(ds, "lndgrid", ng)
                defDim(ds, "time", length(file_years))
                defVar(ds, "YEAR", Int, ("time",))[:] = file_years
                pv = defVar(ds, "PCT_NAT_PFT", Float64,
                    ("natpft", "lndgrid", "time"))
                for t in 1:length(file_years)
                    pv[:, :, t] = pct[:, :, t]
                end
            end

            bounds = CLM.BoundsType(begg=1, endg=1, begl=1, endl=1,
                begc=1, endc=ncol, begp=1, endp=np, level=CLM.BOUNDS_LEVEL_PROC)

            inst = CLM.CLMInstances()
            inst.gridcell = CLM.GridcellData{Float64}()
            inst.gridcell.landunit_indices = fill(CLM.ISPVAL, CLM.MAX_LUNIT, ng)
            inst.gridcell.landunit_indices[CLM.ISTSOIL, 1] = 1
            inst.landunit = CLM.LandunitData{Float64}()
            inst.landunit.itype = [CLM.ISTSOIL]; inst.landunit.wtgcell = [1.0]
            inst.landunit.active = falses(1)
            inst.landunit.coli = [1]; inst.landunit.colf = [ncol]
            inst.column = CLM.ColumnData{Float64}()
            inst.column.landunit = [1]; inst.column.gridcell = [1]
            inst.column.wtlunit = [1.0]; inst.column.wtgcell = [1.0]
            inst.column.active = falses(ncol)
            inst.patch = CLM.PatchData{Float64}()
            inst.patch.gridcell = fill(1, np); inst.patch.landunit = fill(1, np)
            inst.patch.column = fill(1, np); inst.patch.itype = collect(0:(np - 1))
            inst.patch.wtcol = fill(NaN, np); inst.patch.wtlunit = fill(NaN, np)
            inst.patch.wtgcell = fill(NaN, np); inst.patch.active = falses(np)

            ctl = CLM.dyn_subgrid_control_init(flanduse_timeseries=fn,
                                               do_transient_pfts=true)
            config = CLM.CLMDriverConfig()
            @test config.dyn_subgrid === nothing   # default: unwired

            state = CLM.setup_dyn_subgrid!(config, ctl, bounds, inst;
                current_year=2000, natpft_size=natpft_size,
                check_dynpft_consistency=false)

            # WIRING IN PLACE: state is built + attached to the config.
            @test state isa CLM.DynSubgridState
            @test config.dyn_subgrid === state
            @test inst.patch.wtcol ≈ [0.5, 0.3, 0.2]   # year-2000 weights

            # The per-year driver call (the ~line-749 wiring) snaps to year 2001.
            CLM.dynSubgrid_driver!(state, bounds, inst.gridcell, inst.landunit,
                inst.column, inst.patch; year=2001)
            @test inst.patch.wtcol ≈ [0.1, 0.1, 0.8]   # MODULE RAN: weights changed
            @test isapprox(sum(CLM.get_landunit_weight(inst.gridcell, inst.landunit, 1, lt)
                               for lt in 1:CLM.MAX_LUNIT), 1.0; atol=1e-12)
        end
    end

    # ======================================================================
    # 5. CONFIG PLUMBING + default-path invariance.
    # ======================================================================
    @testset "config flags + default invariance" begin
        # New use_ozone flag + dyn_subgrid carrier exist and default OFF.
        c0 = CLM.CLMDriverConfig()
        @test c0.use_ozone == false
        @test c0.dyn_subgrid === nothing

        c1 = CLM.CLMDriverConfig(use_ozone=true)
        @test c1.use_ozone == true
        @test c1.dyn_subgrid === nothing      # still off unless setup_dyn_subgrid!

        # clm_instInit!(use_lch4=true) sizes CH4 state; default leaves it empty.
        CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
        CLM.varcon_init!()
        inst_off = CLM.CLMInstances()
        CLM.clm_instInit!(inst_off; ng=1, nl=1, nc=1, np=1,
                          nlevdecomp_full=CLM.varpar.nlevdecomp_full)
        @test isempty(inst_off.ch4.ch4_surf_flux_tot_col)   # default: unsized

        inst_on = CLM.CLMInstances()
        CLM.clm_instInit!(inst_on; ng=1, nl=1, nc=1, np=1,
                          nlevdecomp_full=CLM.varpar.nlevdecomp_full, use_lch4=true)
        @test length(inst_on.ch4.ch4_surf_flux_tot_col) == 1   # CH4 sized when on
        @test all(inst_on.ch4.ch4_first_time_grc)
    end

end
