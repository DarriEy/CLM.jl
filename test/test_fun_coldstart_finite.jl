# ==========================================================================
# test_fun_coldstart_finite.jl — use_fun=true COLD START stays FINITE (gated).
#
# WHAT THIS LOCKS IN
# ------------------
# CLM5's faithful default is use_fun=.true. + use_flexiblecn=.true. The port
# leaves FUN OFF by default only to keep the default path byte-identical, and a
# stale comment (clm_initialize.jl, since corrected) claimed FUN could not run at
# COLD START because "FUN reads availc/canopy which are NaN before the canopy
# spins up". That reason is DISPROVEN here.
#
# On the LIVE CN driver path calc_gpp_mr_availc! (cn_driver.jl) runs over
# mask_bgc_vegp and populates a FINITE availc BEFORE cnfun! (same mask) reads it,
# so the availc NaN-init is never actually read. A real Bow cold start with
# use_fun=true + REAL forcing photosynthesizes (psnsun finite, LAI develops) and
# keeps availc / npp_growth / leafc finite over several steps.
#
# The earlier "FUN cold-start NaN" was the SAME harness bug as the general canopy
# NaN: an UNSET forc_hgt_{u,t,q}_grc=0 → ustar=0 → the canopy Monin-Obukhov solve
# diverges → t_veg NaN → maintenance respiration NaN → availc = gpp - mr = NaN.
#
# NON-VACUITY. Beyond the finite assertions on the real path, this test also
# re-injects the forc_hgt=0 harness bug and asserts availc goes NaN — i.e. the
# test exercises the real failure mode, and the finiteness is genuinely due to the
# forcing_reader forc_hgt derivation, not a masked/clamped number. It would trip
# if that forcing fix regressed.
#
#   julia +1.12 --project=. --check-bounds=yes test/test_fun_coldstart_finite.jl
# ==========================================================================
using Test, CLM, NCDatasets, Dates, Printf

include(joinpath(@__DIR__, "testdata.jl"))

# One cold-start step. `zero_hgt` re-injects the old harness bug (forc_hgt=0).
# Returns (availc, npp_growth, leafc, psnsun, laisun) as plain Vectors.
function _fun_coldstart_run(; nsteps::Int, zero_hgt::Bool,
                              FSURDAT, FPARAM, FFORCING, FSNOWOPT, FSNOWAGE, INT_SNOW_MAX)
    CLM.cnvegcstate_const.initial_vegC = 100.0
    CLM.rooting_profile_config.rooting_profile_method_water  = CLM.JACKSON_1996_ROOT
    CLM.rooting_profile_config.rooting_profile_method_carbon = CLM.JACKSON_1996_ROOT

    (inst, bounds, filt, tm) = CLM.clm_initialize!(;
        fsurdat=FSURDAT, paramfile=FPARAM,
        start_date=DateTime(2003, 7, 1), dtime=3600,
        use_cn=true, use_luna=true, use_bedrock=true, use_aquifer_layer=false,
        h2osfcflag=0, fsnowoptics=FSNOWOPT, fsnowaging=FSNOWAGE,
        int_snow_max=INT_SNOW_MAX)

    # CLM5's real default: FUN + flexibleCN ON (Bow lnd_in).
    inst.bgc_vegetation.config.use_fun = true
    inst.bgc_vegetation.config.use_flexiblecn = true
    inst.bgc_vegetation.driver_config.use_fun = true
    inst.bgc_vegetation.driver_config.use_flexiblecn = true
    inst.cn_shared_params.use_fun = true
    inst.cn_shared_params.use_flexiblecn = true
    inst.cn_shared_params.br_root = 0.83e-6

    cvf = inst.bgc_vegetation.cnveg_carbonflux_inst
    ccs = inst.bgc_vegetation.cnveg_carbonstate_inst
    ps  = inst.photosyns

    filt_ia = CLM.clump_filter_inactive_and_active
    ng, nc, np = bounds.endg, bounds.endc, bounds.endp
    fr = CLM.ForcingReader(); CLM.forcing_reader_init!(fr, FFORCING)
    config = CLM.CLMDriverConfig(use_cn=true, use_aquifer_layer=false,
                                 use_hydrstress=true, use_luna=true)
    obliqr = CLM.ORB_OBLIQR_DEFAULT

    for _ in 1:nsteps
        step_start = tm.current_date
        CLM.advance_timestep!(tm)
        calday = CLM.get_curr_calday(tm)
        (declin, _) = CLM.compute_orbital(calday)
        (declinm1, _) = CLM.compute_orbital(calday - 3600.0 / CLM.SECSPDAY)
        CLM.init_daylength!(inst.gridcell, declin, declinm1, obliqr, 1:bounds.endg)
        CLM.read_forcing_step!(fr, inst.atm2lnd, step_start, ng, nc)
        CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
        if zero_hgt
            fill!(inst.atm2lnd.forc_hgt_u_grc, 0.0)
            fill!(inst.atm2lnd.forc_hgt_t_grc, 0.0)
            fill!(inst.atm2lnd.forc_hgt_q_grc, 0.0)
        end
        (yr, mon, d, tod) = CLM.get_curr_date(tm)
        CLM.clm_drv!(config, inst, filt, filt_ia, bounds,
                     true, calday, declin, declin, obliqr,
                     false, false, "", false;
                     nstep=tm.nstep, is_first_step=(tm.nstep==1),
                     is_beg_curr_day=CLM.is_beg_curr_day(tm),
                     is_end_curr_day=CLM.is_end_curr_day(tm),
                     is_beg_curr_year=CLM.is_beg_curr_year(tm),
                     dtime=3600.0, mon=mon, day=d,
                     jday=dayofyear(tm.current_date), secs=tod, year=yr,
                     photosyns=inst.photosyns)
    end
    CLM.forcing_reader_close!(fr)
    return (Array(cvf.availc_patch), Array(cvf.npp_growth_patch),
            Array(ccs.leafc_patch), Array(ps.psnsun_patch),
            Array(inst.canopystate.laisun_patch))
end

@testset "FUN cold start stays finite (gated)" begin
    common = joinpath(@__DIR__, "..", "scripts", "fortran_parity_common.jl")
    Base.include(@__MODULE__, common)

    if !isfile(FSURDAT) || !isfile(FPARAM) || !isfile(FFORCING) ||
       !isfile(FSNOWOPT) || !isfile(FSNOWAGE)
        @info "FUN cold-start finiteness: Bow data absent, skipping" FSURDAT FFORCING
        @test_skip isfile(FSURDAT)
    else
        # --- real forcing: FUN cold start must stay finite and photosynthesize ---
        av, npg, lc, psn, lai = Base.invokelatest(_fun_coldstart_run;
            nsteps=4, zero_hgt=false,
            FSURDAT=FSURDAT, FPARAM=FPARAM, FFORCING=FFORCING,
            FSNOWOPT=FSNOWOPT, FSNOWAGE=FSNOWAGE, INT_SNOW_MAX=INT_SNOW_MAX)

        @test all(isfinite, av)          # availc — the field FUN reads
        @test all(isfinite, npg)         # npp_growth — FUN's carbon-cost output
        @test all(isfinite, lc)          # leafc — the veg pool that used to NaN
        # canopy actually spins up: some patch develops LAI and photosynthesizes.
        @test any(>(0.0), lai)
        @test any(x -> isfinite(x) && x > 0.0, psn)
        # leafc seeded from initial_vegC=100 and stays physically sane (Fortran ≈ 99+).
        @test maximum(lc) > 90.0 && maximum(lc) < 101.0

        # --- NON-VACUITY: forc_hgt=0 (the old harness bug) DOES NaN availc ---
        av_bug, _, _, _, _ = Base.invokelatest(_fun_coldstart_run;
            nsteps=1, zero_hgt=true,
            FSURDAT=FSURDAT, FPARAM=FPARAM, FFORCING=FFORCING,
            FSNOWOPT=FSNOWOPT, FSNOWAGE=FSNOWAGE, INT_SNOW_MAX=INT_SNOW_MAX)
        @test any(isnan, av_bug)
    end
end
