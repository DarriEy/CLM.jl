# test_fates_photosynthesis.jl
# ==========================================================================
# FATES DAYTIME PHOTOSYNTHESIS / GPP regression guard.
#
# For a long time EVERY FATES live test ran with the solar-zenith cosine
# (inst.surfalb.coszen_col) = NaN — i.e. perpetual "night". Photosynthesis is
# gated on a valid, daylit coszen, so with coszen=NaN the productivity path was
# silently skipped and stand GPP was identically 0. Four real bugs hid behind
# that "runs in the dark" condition. Every existing FATES test (including the
# 15-day test_fates_spinup.jl carbon-only stability run) still exercises only
# the demography under a hand-rolled cold-soil scaffold and never asserts that
# the canopy actually photosynthesizes.
#
# This test closes that gap. It cold-starts a REAL tropical (Aripuana) stand via
# CLM.clm_initialize!(; use_fates=true) — the only setup that primes a consistent
# surface albedo and hence a FINITE leaf temperature and a real, daylit coszen —
# then loops the actual clm_drv! over one simulated day of REAL daytime forcing
# and ASSERTS the things the dark tests never could:
#   * coszen is FINITE (never NaN) every step, and is clearly daylit (>0) at noon,
#   * total stand GPP accumulates > 0 over the day (photosynthesis path exercised),
#   * btran > 0 (roots take up water — the beta-stress path is live),
#   * the daily FATES carbon balance is finite and holds (TotalBalanceCheck),
#   * live stand carbon changes by a finite amount across the day boundary
#     (the daily growth/allocation path actually ran).
#
# It reuses build() + the census/carbon helpers from the validated harness
# scripts/fates_longhorizon.jl (include()-guarded: including it only defines
# functions). It is deliberately FAST — one day + one step (49 half-hourly
# clm_drv! steps), vs the 720 steps of test_fates_spinup.jl.
#
# The stand needs the Aripuana surfdata/params + forcing that clm_initialize!
# already reads. If that data is absent (e.g. bare CI), the test SKIPS cleanly
# rather than faking a pass — the productivity assertions are only meaningful
# against a real daylit stand.
# ==========================================================================
using Test

# Brings in CLM, `const _C = CLM`, `const DATA`, and build() + the census /
# total_site_carbon / max_dbh helpers. The script guards its own entrypoint
# (`if abspath(PROGRAM_FILE)==@__FILE__`), so including it just defines these.
include(joinpath(@__DIR__, "..", "scripts", "fates_longhorizon.jl"))

@testset "FATES daytime photosynthesis / GPP (real clm_drv! stand)" begin
    fsurdat = get(ENV, "FATES_FSURDAT",
        "$DATA/domain_Aripuana_Amazon/settings/CLM/parameters/surfdata_clm.nc")
    paramfile = get(ENV, "FATES_PARAMFILE",
        "$DATA/domain_Aripuana_Amazon/settings/CLM/parameters/clm5_params.nc")
    forcing = get(ENV, "FATES_FORCING",
        "$DATA/domain_Aripuana_Amazon/data/forcing/CLM_input/clmforc.2004.nc")

    if !(isfile(fsurdat) && isfile(paramfile) && isfile(forcing))
        # Data-dependent test; nothing to fake. Report clearly + skip.
        @info "FATES photosynthesis test SKIPPED — Aripuana stand data not present" fsurdat paramfile forcing
        @test_skip true
    else
        # ---- stand-level probes over the FATES site linked list ----
        stand_gpp(site) = begin           # Σ cohort gpp_tstep × n  [kgC/m2/step]
            g = 0.0; cp = site.oldest_patch
            while cp !== nothing
                cc = cp.tallest
                while cc !== nothing
                    isfinite(cc.gpp_tstep) && (g += cc.gpp_tstep * cc.n)
                    cc = cc.shorter
                end
                cp = cp.younger
            end
            return g
        end
        max_btran(site) = begin           # max over patches of the FT beta-stress
            b = 0.0; cp = site.oldest_patch
            while cp !== nothing
                try
                    b = max(b, maximum(x -> isfinite(x) ? x : 0.0, cp.btran_ft))
                catch
                end
                cp = cp.younger
            end
            return b
        end

        # One simulated day + one step so the day-boundary daily dynamics fires once.
        function run_photo_day(inst, fates, config, bounds, filt, filt_ia, forcing)
            site = fates.sites[1]; photosyns = inst.photosyns
            dtime = 1800.0
            steps_per_day = Int(round(86400 / dtime))         # 48
            nsteps = steps_per_day + 1                          # 49 → cross one day boundary
            fr = _C.ForcingReader(); _C.forcing_reader_init!(fr, forcing); fr.interp_time = true
            start_date = DateTime(2004, 1, 1)

            c0 = total_site_carbon(site)
            acc_gpp = 0.0; max_coszen = -Inf; max_btr = 0.0
            coszen_nan = false; day_balance_ok = false; days_advanced = 0
            try
                for i in 1:nsteps
                    step_start = start_date + Second((i - 1) * Int(dtime))
                    is_beg = (Dates.hour(step_start) == 0 && Dates.minute(step_start) == 0 &&
                              Dates.second(step_start) == 0)
                    _C.read_forcing_step!(fr, inst.atm2lnd, step_start, 1, 1;
                        gridcell_latdeg=inst.gridcell.latdeg,
                        gridcell_londeg=inst.gridcell.londeg, dtime=Int(dtime))
                    inst.atm2lnd.forc_topo_grc[1] = 200.0; inst.topo.topo_col[1] = 200.0
                    _C.downscale_forcings!(bounds, inst.atm2lnd, inst.column,
                                           inst.landunit, inst.topo)
                    sod = Dates.hour(step_start) * 3600 + Dates.minute(step_start) * 60
                    calday = Dates.dayofyear(step_start) + sod / 86400.0
                    (declin, _e) = _C.compute_orbital(calday)
                    nextsw_cday = calday + Int(dtime) / 86400.0

                    _C.clm_drv!(config, inst, filt, filt_ia, bounds, true, nextsw_cday,
                        declin, declin, 0.4091, false, false, "20260101", false;
                        nstep=i, is_first_step=(i == 1), is_beg_curr_day=is_beg,
                        is_end_curr_day=false, is_beg_curr_year=false, dtime=dtime,
                        mon=Dates.month(step_start), day=Dates.day(step_start),
                        secs=sod, jday=Dates.dayofyear(step_start), photosyns=photosyns)

                    cz = inst.surfalb.coszen_col[1]
                    isnan(cz) && (coszen_nan = true)
                    isfinite(cz) && (max_coszen = max(max_coszen, cz))
                    acc_gpp += stand_gpp(site)
                    max_btr = max(max_btr, max_btran(site))

                    # Day boundary (after step 1): the daily demographic + balance step
                    # ran inside clm_drv!. Re-check the balance holds + is finite.
                    if is_beg && i > 1
                        days_advanced += 1
                        day_balance_ok = (_C.TotalBalanceCheck(site, -1) === nothing)
                    end
                end
            finally
                _C.forcing_reader_close!(fr)
            end
            cF = total_site_carbon(site)
            return (; c0, cF, acc_gpp, max_coszen, max_btr, coszen_nan,
                    day_balance_ok, days_advanced)
        end

        inst, fates, config, bounds, filt, filt_ia = build()
        site = fates.sites[1]
        r = run_photo_day(inst, fates, config, bounds, filt, filt_ia, forcing)

        # --- 1. coszen is real, not NaN: the stand is NOT "running in the dark" ---
        @test !r.coszen_nan                       # never NaN over the day (the old-bug guard)
        @test isfinite(r.max_coszen)
        @test r.max_coszen > 0.1                  # clearly daylit at solar noon

        # --- 2. the photosynthesis path actually produced carbon uptake ---
        @test isfinite(r.acc_gpp)
        @test r.acc_gpp > 0.0                     # stand GPP > 0 (was silently 0 in the dark)

        # --- 3. roots took up water (beta-stress path is live) ---
        @test isfinite(r.max_btr)
        @test r.max_btr > 0.0

        # --- 4. the daily FATES carbon balance is finite and holds ---
        @test r.days_advanced >= 1
        @test r.day_balance_ok
        @test _C.TotalBalanceCheck(site, -1) === nothing

        # --- 5. the growth path ran: live stand carbon changed by a finite amount ---
        @test isfinite(r.c0) && isfinite(r.cF) && r.cF > 0.0
        @test isfinite(r.cF - r.c0)
        @test r.cF != r.c0                        # daily allocation/growth updated the pools
    end
end
