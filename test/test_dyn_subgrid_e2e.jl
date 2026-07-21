#!/usr/bin/env julia
# ==========================================================================
# End-to-end transient-land-use test: drives the REAL clm_drv! timestep loop
# with config.dyn_subgrid set, exercising the wired dynSubgrid_driver! hook at
# a year boundary (clm_driver.jl:925 — the `if config.dyn_subgrid !== nothing &&
# is_beg_curr_year && !is_first_step` block, mirroring CTSM clm_driver.F90:367
# `call dynSubgrid_driver`).
#
# All the OTHER dyn_subgrid tests (test_dyn_subgrid_driver.jl,
# test_dyn_subgrid_cnbal_wiring.jl, test_dyn_cons_biogeophys.jl) call
# dynSubgrid_driver! DIRECTLY on a synthetic domain. This is the only test that
# drives it through the full clm_drv! step on the real Bow surface dataset,
# proving:
#   * default (no transient dataset) leaves config.dyn_subgrid === nothing
#     (the gate that keeps a non-transient run byte-identical);
#   * setup_dyn_subgrid! snaps the natural-PFT patch weights to the file's
#     year-2000 values at init;
#   * a year-boundary clm_drv! step (is_beg_curr_year=true, !is_first_step)
#     fires the hook and snaps the weights to the file's year-2001 values —
#     the weight change propagates through the WHOLE driver, not a direct call;
#   * the biogeophysics stays water-conserving through the transient hook: the
#     column water balance error after the transient step is at machine
#     precision (the FATAL errh2o check — balance_check.jl — did not trip).
#
# Biogeophys-only (use_cn=false) on purpose: a transient-PFT change redistributes
# patch weights WITHIN a fixed-area column, so this exercises the hook + weight
# propagation + water-balance integrity. Column-area-change dynbal fluxes (C/N and
# heat/water) are covered by the component tests above.
# ==========================================================================

using Test
using CLM
using NCDatasets

include(joinpath(@__DIR__, "testdata.jl"))

@testset "dyn_subgrid end-to-end through clm_drv!" begin
    fsurdat, paramfile = bow_params()
    if !isfile(fsurdat) || !isfile(paramfile)
        testdata_missing("dyn_subgrid e2e", fsurdat, paramfile)
        return
    end

    (inst, bounds, filt, tm) = CLM.clm_initialize!(;
        fsurdat=fsurdat, paramfile=paramfile, use_cn=false)

    config = CLM.CLMDriverConfig(use_cn=false)
    @test config.dyn_subgrid === nothing   # DEFAULT: gate off → byte-identical path

    filt_ia = CLM.clump_filter_inactive_and_active
    dtime = 1800.0
    calday = 172.5
    (declin, eccf) = CLM.compute_orbital(calday)
    nextsw_cday = calday + dtime / CLM.SECSPDAY
    ng = bounds.endg

    # --- forcing / soil moisture / vegetation setup (as in test_cn_integration.jl) ---
    a2l = inst.atm2lnd
    CLM._setup_calib_forcing!(a2l, 285.0, ng)
    CLM.downscale_forcings!(bounds, a2l, inst.column, inst.landunit, inst.topo)
    CLM._init_calib_soil_moisture!(inst, bounds)
    CLM.interp_monthly_veg!(inst.satellite_phenology; kmo=6, kda=21)
    cs = inst.canopystate
    wdb = inst.water.waterdiagnosticbulk_inst
    pch = inst.patch
    CLM.satellite_phenology!(inst.satellite_phenology, cs, wdb, pch,
                             filt.nolakep, bounds.begp:bounds.endp)
    for p in bounds.begp:bounds.endp
        cs.frac_veg_nosno_patch[p] = cs.frac_veg_nosno_alb_patch[p]
    end
    CLM.set_exposedvegp_filter!(filt, bounds, cs.frac_veg_nosno_patch)

    # --- build a transient PCT_NAT_PFT dataset for the natveg landunit ---------
    # The natural-veg patches on this domain have itype [0, 1, 12]; dynpft_interp!
    # maps patch p → file natpft column m = itype+1, so we size natpft=15 and place
    # the three PFTs at columns 1, 2, 13 (summing to 100% per year).
    natpft_size = 15
    file_years = [2000, 2001]
    pct = zeros(Float64, natpft_size, ng, length(file_years))
    pct[1, 1, 1] = 5.0;  pct[2, 1, 1] = 60.0; pct[13, 1, 1] = 35.0   # 2000
    pct[1, 1, 2] = 10.0; pct[2, 1, 2] = 30.0; pct[13, 1, 2] = 60.0   # 2001

    mktempdir() do dir
        fn = joinpath(dir, "flanduse_pft.nc")
        NCDataset(fn, "c") do ds
            defDim(ds, "natpft", natpft_size)
            defDim(ds, "lndgrid", ng)
            defDim(ds, "time", length(file_years))
            defVar(ds, "YEAR", Int, ("time",))[:] = file_years
            pv = defVar(ds, "PCT_NAT_PFT", Float64, ("natpft", "lndgrid", "time"))
            for t in 1:length(file_years)
                pv[:, :, t] = pct[:, :, t]
            end
        end

        # PROC-level bounds view (the dyn_subgrid driver asserts PROC level; the
        # single-clump bounds ARE the proc bounds — same as clm_drv_core!).
        bp = CLM.BoundsType(begg=bounds.begg, endg=bounds.endg,
                            begl=bounds.begl, endl=bounds.endl,
                            begc=bounds.begc, endc=bounds.endc,
                            begp=bounds.begp, endp=bounds.endp,
                            level=CLM.BOUNDS_LEVEL_PROC)

        ctl = CLM.dyn_subgrid_control_init(flanduse_timeseries=fn,
                                           do_transient_pfts=true)
        CLM.setup_dyn_subgrid!(config, ctl, bp, inst; current_year=2000,
                               natpft_size=natpft_size,
                               check_dynpft_consistency=false)

        # WIRING attached; init snapped the natveg patch weights to year 2000.
        @test config.dyn_subgrid isa CLM.DynSubgridState
        natveg = [p for p in bounds.begp:bounds.endp
                  if inst.landunit.itype[inst.patch.landunit[p]] == CLM.ISTSOIL]
        @test pch.wtcol[natveg] ≈ [0.05, 0.6, 0.35]
        wt_before = copy(pch.wtcol[natveg])

        # Step 1: ordinary first step — hook does NOT fire (is_first_step).
        CLM.clm_drv!(config, inst, filt, filt_ia, bounds,
            true, nextsw_cday, declin, declin, CLM.ORB_OBLIQR_DEFAULT,
            false, false, "", false;
            nstep=1, is_first_step=true, is_beg_curr_day=true,
            dtime=dtime, mon=6, day=21, year=2000, is_beg_curr_year=false,
            photosyns=inst.photosyns)
        @test pch.wtcol[natveg] ≈ wt_before   # no change on a non-boundary step

        # Step 2: YEAR BOUNDARY — the dynSubgrid_driver! hook fires and snaps the
        # weights to the file's year-2001 values, THROUGH the full driver.
        CLM.clm_drv!(config, inst, filt, filt_ia, bounds,
            true, nextsw_cday, declin, declin, CLM.ORB_OBLIQR_DEFAULT,
            false, false, "", false;
            nstep=2, is_first_step=false, is_beg_curr_day=false,
            dtime=dtime, mon=6, day=21, year=2001, is_beg_curr_year=true,
            photosyns=inst.photosyns)

        # MODULE RAN end-to-end: weights changed to the year-2001 file values.
        @test pch.wtcol[natveg] ≈ [0.1, 0.3, 0.6]
        @test !(pch.wtcol[natveg] ≈ wt_before)

        # Biogeophysics stayed water-conserving through the transient hook: the
        # column water balance error is at machine precision (the FATAL errh2o
        # check did not trip on the weight-change step).
        errh2o = inst.water.waterbalancebulk_inst.errh2o_col
        natveg_col = inst.patch.column[natveg[1]]
        @test abs(errh2o[natveg_col]) < 1e-8

        println("  dyn_subgrid e2e: weights ", wt_before, " -> ",
                pch.wtcol[natveg], ", errh2o_col=", errh2o[natveg_col])
    end
end
