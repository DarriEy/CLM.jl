#!/usr/bin/env julia
# ==========================================================================
# End-to-end TRANSIENT-land-use test WITH use_cn=true, driving the REAL clm_drv!
# timestep loop across a year boundary. The sibling test_dyn_subgrid_e2e.jl runs
# the same hook use_cn=false (water-only); this one exercises the CN mass-balance
# path through the transient dwt (dynamic-weight) fluxes.
#
# It guards three bugs that only fire on a use_cn transient year-boundary step —
# none of which the non-transient CN audit or the use_cn=false e2e test can see:
#
#   1. dyn_cnbal_col! looped the vertical index over size(decomp_*_vr_col, 2)
#      (== nlevdecomp_full, 25 on this grid) instead of nlevdecomp (== 20 ==
#      length(dzsoi_decomp)), so `dzsoi_decomp[j]` ran off the end and the
#      boundary step CRASHED with a BoundsError before any check.
#      (Fortran DynamicColumnAdjustments loops `do j = 1, nlevdecomp`.)
#
#   2. The per-step CNVeg flux reset (_zero_cnveg_flux_arrays!, the port's
#      SetValues) blanket-zeroed EVERY flux field — including the dwt_* landuse
#      conversion/seed/product fluxes that dynSubgrid_driver had just written at
#      the year boundary. Fortran's SetValues never touches dwt_ fields (they are
#      managed only by ZeroDwt, at step START, before dynSubgrid). Result: the
#      gridcell C/N balance saw the store drop but read every dwt term as 0.
#
#   3. ZeroDwt (CNVegetation InitEachTimeStep) was a documented no-op, so once (2)
#      stopped wiping dwt in the CN driver, the dwt grc fluxes had to be zeroed at
#      step start instead — otherwise a non-boundary step would carry a stale
#      conversion flux into its balance.
#
# The COLUMN-level C and N balances CLOSE to machine precision through the
# transient step, and the dwt conversion flux survives to the end-of-step check
# (it is a REAL term, not wiped). BOTH steps now run with the FATAL thresholds
# armed — including STEP 1: an earlier ~0.086 gN STEP-1 column-N residual proved to
# be a FIXTURE inconsistency (sminn_vr seeded while its NH4/NO3 source pools were
# left at zero, which the N-state update then discarded), not a port gap; the CN
# pool seeding below now seeds NH4/NO3 consistently so step 1 conserves N exactly.
#
# The GRIDCELL balance ALSO closes now that the wood/crop PRODUCT POOLS are wired
# end-to-end (PR #292: cn_products_update! runs on the live cn_driver_no_leaching!
# path, so tot_woodprod_grc / cropprod1_grc advance from the year-boundary
# product-gain fluxes — the sink the gridcell check subtracts). This test therefore
# runs the REAL, FATAL C/N balance check at Fortran's default thresholds (cerror
# 1e-7, nerror 1e-3): it must pass inside clm_drv!, and afterwards we assert
# errcb_grc / errnb_grc < 1e-8 AND that the wood product pool actually advanced
# (~60 gC into prod10/prod100) — proving the closure is NON-VACUOUS (product mass
# really moved), not a trivial pass on zero flux. Were products still unwired the
# driver would throw on the harvest/conversion residual (the ~37 gC / ~0.085 gN
# class guarded directly by test_cn_products_gridcell_balance.jl).
# ==========================================================================

using Test
using CLM
using NCDatasets

include(joinpath(@__DIR__, "testdata.jl"))

@testset "dyn_subgrid + use_cn end-to-end through clm_drv!" begin
    fsurdat, paramfile = bow_params()
    if !isfile(fsurdat) || !isfile(paramfile)
        testdata_missing("dyn_subgrid CN e2e", fsurdat, paramfile)
        return
    end

    (inst, bounds, filt, tm) = CLM.clm_initialize!(;
        fsurdat=fsurdat, paramfile=paramfile, use_cn=true)

    config = CLM.CLMDriverConfig(use_cn=true)
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

    # --- CN pool init (replace NaN inits, seed realistic veg + soil pools) ------
    bgc_cs = inst.bgc_vegetation.cnveg_carbonstate_inst
    bgc_ns = inst.bgc_vegetation.cnveg_nitrogenstate_inst
    bgc_cf = inst.bgc_vegetation.cnveg_carbonflux_inst
    bgc_nf = inst.bgc_vegetation.cnveg_nitrogenflux_inst
    for s in (bgc_cs, bgc_ns, bgc_cf, bgc_nf,
              inst.soilbiogeochem_carbonstate, inst.soilbiogeochem_nitrogenstate,
              inst.soilbiogeochem_carbonflux, inst.soilbiogeochem_nitrogenflux)
        for fn in fieldnames(typeof(s))
            v = getfield(s, fn)
            v isa AbstractArray{Float64} && replace!(v, NaN => 0.0)
        end
    end
    for p in bounds.begp:bounds.endp
        inst.patch.itype[p] > 0 || continue
        bgc_cs.leafc_patch[p] = 100.0;  bgc_cs.frootc_patch[p] = 50.0
        bgc_cs.livestemc_patch[p] = 200.0; bgc_cs.deadstemc_patch[p] = 500.0
        bgc_ns.leafn_patch[p] = 3.0; bgc_ns.frootn_patch[p] = 1.5
        bgc_ns.livestemn_patch[p] = 2.0
    end
    soilbgc_cs = inst.soilbiogeochem_carbonstate
    soilbgc_ns = inst.soilbiogeochem_nitrogenstate
    nlevdecomp = CLM.varpar.nlevdecomp
    for c in bounds.begc:bounds.endc
        for j in 1:nlevdecomp, p in 1:min(7, size(soilbgc_cs.decomp_cpools_vr_col, 3))
            soilbgc_cs.decomp_cpools_vr_col[c, j, p] = 10.0
            soilbgc_ns.decomp_npools_vr_col[c, j, p] = 0.5
        end
        soilbgc_ns.sminn_vr_col[c, 1:nlevdecomp] .= 0.01
        # sminn_vr is a DERIVED diagnostic under use_nitrif_denitrif (the default):
        # every SoilBiogeochem N-state update OVERWRITES it with smin_nh4_vr +
        # smin_no3_vr (nitrif_denitrif.jl:628 == SoilBiogeochemNStateUpdate1Mod.F90),
        # and totn_col sums sminn (soil_bgc_nitrogen_state.jl:559). Seeding sminn_vr
        # WITHOUT seeding the NH4/NO3 pools it is derived from is an inconsistency the
        # real cold start never produces (SoilBiogeochemNitrogenStateType.F90:InitCold
        # sets sminn_vr = smin_nh4_vr = smin_no3_vr = 0). Left unfixed, the seeded
        # mineral N is silently discarded on the first N-state update — ~0.086 gN/m2 of
        # column mineral N vanishing with no matching output flux — so the STEP-1 column
        # (and gridcell) N balance could not close. Seed NH4/NO3 consistently so the
        # first step conserves N exactly.
        soilbgc_ns.smin_nh4_vr_col[c, 1:nlevdecomp] .= 0.01
        soilbgc_ns.smin_no3_vr_col[c, 1:nlevdecomp] .= 0.0
    end

    # --- transient PCT_NAT_PFT dataset (natveg itypes [0,1,12] -> cols 1,2,13) ---
    natpft_size = 15
    file_years = [2000, 2001]
    pct = zeros(Float64, natpft_size, ng, length(file_years))
    pct[1, 1, 1] = 5.0;  pct[2, 1, 1] = 60.0; pct[13, 1, 1] = 35.0   # 2000
    pct[1, 1, 2] = 10.0; pct[2, 1, 2] = 30.0; pct[13, 1, 2] = 60.0   # 2001

    # Save/restore global balance-check state so this test does not leak into others.
    saved_enabled = CLM.cn_balance_check_enabled()
    bal = inst.bgc_vegetation.cn_balance_inst

    mktempdir() do dir
        fn = joinpath(dir, "flanduse_pft.nc")
        NCDataset(fn, "c") do ds
            defDim(ds, "natpft", natpft_size)
            defDim(ds, "lndgrid", ng)
            defDim(ds, "time", length(file_years))
            defVar(ds, "YEAR", Int, ("time",))[:] = file_years
            pv = defVar(ds, "PCT_NAT_PFT", Float64, ("natpft", "lndgrid", "time"))
            for t in 1:length(file_years); pv[:, :, t] = pct[:, :, t]; end
        end

        bp = CLM.BoundsType(begg=bounds.begg, endg=bounds.endg,
                            begl=bounds.begl, endl=bounds.endl,
                            begc=bounds.begc, endc=bounds.endc,
                            begp=bounds.begp, endp=bounds.endp,
                            level=CLM.BOUNDS_LEVEL_PROC)
        ctl = CLM.dyn_subgrid_control_init(flanduse_timeseries=fn, do_transient_pfts=true)
        CLM.setup_dyn_subgrid!(config, ctl, bp, inst; current_year=2000,
                               natpft_size=natpft_size, check_dynpft_consistency=false)
        @test config.dyn_subgrid isa CLM.DynSubgridState

        natveg = [p for p in bounds.begp:bounds.endp
                  if inst.landunit.itype[inst.patch.landunit[p]] == CLM.ISTSOIL]
        wt_before = copy(pch.wtcol[natveg])
        @test pch.wtcol[natveg] ≈ [0.05, 0.6, 0.35]

        # Turn the balance check ON with Fortran's REAL, FATAL thresholds armed for
        # STEP 1 too. The ~0.086 gN STEP-1 column-N residual this test used to scope
        # away (thresholds = Inf) was NOT a first-step-allocation port gap: it was a
        # FIXTURE inconsistency in this very test. Under use_nitrif_denitrif (the
        # default) sminn_vr is a DERIVED diagnostic — the N-state update overwrites it
        # with smin_nh4_vr + smin_no3_vr, and totn_col sums sminn — so seeding sminn_vr
        # while leaving NH4/NO3 at zero silently discarded ~0.086 gN/m2 of mineral N on
        # the first step (the real InitCold keeps all three = 0). The CN pool seeding
        # above now seeds NH4/NO3 consistently, so the first-step column AND gridcell
        # C/N balances close to machine precision; the fatal check must PASS here.
        CLM.cn_balance_check_enabled!(true)
        bal.cerror = 1.0e-7; bal.nerror = 1.0e-3
        bal.cwarning = 1.0e-8; bal.nwarning = 1.0e-7

        # Step 1: ordinary first step — hook does NOT fire (is_first_step).
        CLM.clm_drv!(config, inst, filt, filt_ia, bounds,
            true, nextsw_cday, declin, declin, CLM.ORB_OBLIQR_DEFAULT,
            false, false, "", false;
            nstep=1, is_first_step=true, is_beg_curr_day=true,
            dtime=dtime, mon=6, day=21, year=2000, is_beg_curr_year=false,
            photosyns=inst.photosyns)
        @test pch.wtcol[natveg] ≈ wt_before

        # STEP-1 column C AND N balances close to machine precision — assert the FATAL
        # column-N threshold (nerror = 1e-3) and the tighter nwarning (1e-7), exactly
        # as done for the year-boundary step below. This is the assertion that proves
        # the first-step column-N accounting is CLOSED (previously ~0.086 gN, hidden by
        # the Inf scoping). Note clm_drv! above already ran this check FATALLY, so a
        # regression aborts before reaching these lines.
        soil_s1 = [c for c in bounds.begc:bounds.endc if filt.bgc_soilc[c]]
        errnb_col_s1 = maximum(abs, bal.errnb_col[soil_s1])
        errcb_col_s1 = maximum(abs, bal.errcb_col[soil_s1])
        @test errnb_col_s1 < 1e-3      # << nerror (fatal); here ~6e-15
        @test errnb_col_s1 < 1e-7      # even clears the (tighter) nwarning
        @test errcb_col_s1 < 1e-7      # << cerror; here ~1e-13
        println("    STEP 1 errcb_col=", errcb_col_s1, "  errnb_col=", errnb_col_s1)

        # Keep Fortran's REAL, FATAL C/N balance thresholds for the year-boundary
        # step. With the wood/crop PRODUCT POOLS wired end-to-end (PR #292:
        # cn_products_update! runs on the live cn_driver_no_leaching! path, advancing
        # tot_woodprod_grc from the year-boundary product-gain fluxes — the sink the
        # GRIDCELL check subtracts), both the column AND gridcell C/N balances close,
        # so the fatal check must PASS inside the step below. Were products still
        # unwired the driver would THROW here on the harvest/conversion residual.
        bal.cerror = 1.0e-7; bal.nerror = 1.0e-3
        bal.cwarning = 1.0e-8; bal.nwarning = 1.0e-7

        # Step 2: YEAR BOUNDARY — dynSubgrid_driver! fires, computing dwt fluxes.
        # Before the fixes this CRASHED (BoundsError in dyn_cnbal_col!); with the
        # dwt fluxes wiped it left the balance unaccounted.
        CLM.clm_drv!(config, inst, filt, filt_ia, bounds,
            true, nextsw_cday, declin, declin, CLM.ORB_OBLIQR_DEFAULT,
            false, false, "", false;
            nstep=2, is_first_step=false, is_beg_curr_day=false,
            dtime=dtime, mon=6, day=21, year=2001, is_beg_curr_year=true,
            photosyns=inst.photosyns)

        # Module ran end-to-end: weights snapped to the year-2001 file values.
        @test pch.wtcol[natveg] ≈ [0.1, 0.3, 0.6]

        soil = [c for c in bounds.begc:bounds.endc if filt.bgc_soilc[c]]

        # (1) No crash + (2)/(3) the COLUMN C and N balances CLOSE to machine
        # precision through the transient step (fatal bounds: cerror 1e-7, nerror 1e-3).
        errnb_col = maximum(abs, bal.errnb_col[soil])
        errcb_col = maximum(abs, bal.errcb_col[soil])
        @test errnb_col < 1e-3      # << nerror; here ~1e-15
        @test errnb_col < 1e-7      # even clears the (tighter) nwarning
        @test errcb_col < 1e-7      # << cerror; here ~1e-13

        # The dwt CONVERSION flux survives to the end-of-step balance (it is a REAL
        # term now, not wiped to 0 by the per-step SetValues reset).
        nf = inst.bgc_vegetation.cnveg_nitrogenflux_inst
        cf = inst.bgc_vegetation.cnveg_carbonflux_inst
        @test nf.dwt_conv_nflux_grc[1] > 0.0
        @test cf.dwt_conv_cflux_grc[1] > 0.0

        # (4) The GRIDCELL C and N balances ALSO close to machine precision through
        # the transient step — the product-pool advance is the sink that absorbs the
        # year-boundary product-gain fluxes (without it the gridcell store would lose
        # that mass and errcb_grc/errnb_grc would carry the harvest/conversion
        # residual). These clear even the tighter cwarning (1e-8) / nwarning (1e-7);
        # the fatal gridcell check (enabled above) already passed inside clm_drv!.
        grc = bounds.begg:bounds.endg
        errcb_grc = maximum(abs, bal.errcb_grc[grc])
        errnb_grc = maximum(abs, bal.errnb_grc[grc])
        @test errcb_grc < 1e-8      # << cerror (1e-7); here ~1e-13
        @test errnb_grc < 1e-8      # << nerror (1e-3); here ~1e-16

        # NON-VACUOUS proof: the wood PRODUCT POOL actually advanced this step — the
        # transient PFT contraction fed deadstem/livestem C & N into prod10/prod100
        # (partitioned by pprod10/pprod100). A gridcell balance that "closed" on a
        # zero product gain would be trivially vacuous; assert real mass moved.
        cprod = inst.bgc_vegetation.c_products_inst
        nprod = inst.bgc_vegetation.n_products_inst
        @test cprod.tot_woodprod_grc[1] > 1.0                          # ~60 gC moved into wood products
        @test nprod.tot_woodprod_grc[1] > 0.0                          # N wood product gain (nonzero)
        @test cprod.tot_woodprod_grc[1] ≈ cprod.prod10_grc[1] + cprod.prod100_grc[1]  # sum == pools

        println("  dyn_subgrid CN e2e: weights ", wt_before, " -> ", pch.wtcol[natveg])
        println("    errcb_col=", errcb_col, "  errnb_col=", errnb_col)
        println("    errcb_grc=", errcb_grc, "  errnb_grc=", errnb_grc)
        println("    dwt_conv_nflux_grc*dt=", nf.dwt_conv_nflux_grc[1] * dtime,
                "  dwt_conv_cflux_grc*dt=", cf.dwt_conv_cflux_grc[1] * dtime)
        println("    tot_woodprod_grc: C=", cprod.tot_woodprod_grc[1],
                " (prod10=", cprod.prod10_grc[1], " prod100=", cprod.prod100_grc[1], ")",
                "  N=", nprod.tot_woodprod_grc[1])
    end

    CLM.cn_balance_check_enabled!(saved_enabled)
end
