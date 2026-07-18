# ==========================================================================
# test_fortran_parity_luna_decomp.jl — gated parity for LUNA cadence + the
# soil-BGC N source-term *integrated* effect (mineral-N pool deltas).
#
# Companion to test_fortran_parity_cn.jl. Two independent checks, each over the
# instrumented Bow summer BGC dumps (machine-local; SKIPPED when absent):
#
#  1. LUNA CADENCE + INVARIANT (dump-only, no model run — fast):
#       update_photosynthesis_capacity! (LUNA acclimation) updates vcmx25_z /
#       jmx25_z ONLY at end-of-day (is_time_to_run_luna == is_end_curr_day,
#       luna.jl:559). In the dump window n1757845..n1757872 the only end-of-day
#       step is n1757856 (runs 23:00->00:00). The dumps prove it directly:
#       max|after - before| of vcmx25_z is exactly 0 on every step EXCEPT
#       n1757856. Also the LUNA invariant vcmx25_z(post) == vcmx25_z_last_valid
#       must hold at the EOD step.
#       (NOTE a real residual the probe surfaced: CLM.jl's clm_drv! never calls
#       update_photosynthesis_capacity! — the harness injects vcmx25_z. And the
#       dump omits the LUNA climate accumulators t_veg10_day/night/fpsn24/dayl,
#       so a no-inject recompute-vs-dump is not possible. See
#       scripts/fortran_parity_luna.jl. Hence only the cadence+invariant are
#       asserted here.)
#
#  2. DECOMP / N-CYCLE SOURCE TERMS (one CN step, nstep 1757852):
#       The _P rate-print dumps (GROSS_NMIN_VR_P, ACT_IMMOB_NH4_VR_P, F_NIT_VR_P,
#       POT_F_NIT_VR_P, F_DENIT_VR_P, SMIN_NH4_TO_PLANT_VR_P) are identically
#       ZERO at every boundary/step (a dead instrumentation hook), so no direct
#       rate parity is possible. Instead we validate the rates' INTEGRATED effect
#       via the per-layer mineral-N pool deltas (smin_no3_vr / smin_nh4_vr /
#       sminn_vr): Julia recompute vs the Fortran after-dump.
#         - smin_no3_vr: per-layer rel err < 1e-4, top-10 ratio ~0.985  (clean).
#         - smin_nh4_vr / sminn_vr: bounded; ratio ~1.19 / ~1.14 is the KNOWN,
#           documented closed-GPP/N-uptake residual (not a decomp-rate bug), so
#           tested only for boundedness with headroom.
#       See scripts/fortran_parity_decomprates.jl.
#
# GATED: Fortran BGC reference dumps live outside the repo. When absent the test
# is skipped (not failed) so CI / other machines stay green.
# ==========================================================================

using Test
using CLM
using NCDatasets
using Dates

include(joinpath(@__DIR__, "testdata.jl"))

@testset "Fortran LUNA-cadence + decomp N-source parity (gated)" begin
    common  = joinpath(@__DIR__, "..", "scripts", "fortran_parity_common.jl")
    dumpdir = symfluence_path("clm_bgc_spinup", "bgc_ref_summer")
    date_base = 1753153

    window    = 1757845:1757872
    eod_nstep = 1757856        # the only end-of-day step in the window
    nstep_dec = 1757852        # the daytime step used for the decomp N-delta check

    have_window = all(isfile(joinpath(dumpdir, "pdump_before_step_n$(n).nc")) &&
                      isfile(joinpath(dumpdir, "pdump_after_hydrologydrainage_n$(n).nc"))
                      for n in window)

    if !(isfile(common) && have_window)
        @info "Fortran LUNA/decomp parity test SKIPPED (BGC reference dumps not present)"
        @test_skip have_window
    else
        # ----------------------------------------------------------------------
        # 1. LUNA cadence + invariant (dump-only)
        # ----------------------------------------------------------------------
        @testset "LUNA cadence (update only at end-of-day n$eod_nstep)" begin
            updated = Int[]
            for n in window
                db = NCDataset(joinpath(dumpdir, "pdump_before_step_n$(n).nc"), "r")
                da = NCDataset(joinpath(dumpdir, "pdump_after_hydrologydrainage_n$(n).nc"), "r")
                dv = maximum(abs.(Float64.(da["vcmx25_z"][:, :]) .- Float64.(db["vcmx25_z"][:, :])))
                dj = maximum(abs.(Float64.(da["jmx25_z"][:, :])  .- Float64.(db["jmx25_z"][:, :])))
                (dv > 1e-8 || dj > 1e-8) && push!(updated, n)
                close(db); close(da)
            end
            @info "LUNA updated on step(s)" updated
            # LUNA acclimation fires on exactly one step in the window, the EOD step.
            @test updated == [eod_nstep]

            # LUNA invariant: at the EOD step the post vcmx25_z equals the
            # *_last_valid copy the update writes.
            da = NCDataset(joinpath(dumpdir, "pdump_after_hydrologydrainage_n$(eod_nstep).nc"), "r")
            if haskey(da, "vcmx25_z_last_valid_patch")
                inv_diff = maximum(abs.(Float64.(da["vcmx25_z"][:, :]) .-
                                        Float64.(da["vcmx25_z_last_valid_patch"][:, :])))
                @info "LUNA invariant vcmx25_z == vcmx25_z_last_valid" maxdiff=inv_diff
                @test inv_diff < 1e-8
            end
            close(da)
        end

        # ----------------------------------------------------------------------
        # 2. Decomp / N-cycle source terms (integrated via mineral-N pool deltas)
        # ----------------------------------------------------------------------
        @testset "decomp N-source integrated delta (n$nstep_dec)" begin
            Base.include(@__MODULE__, common)   # build_bow_inst / run_one_parity_step! / FFORCING
            forcing2002 = replace(FFORCING, "clmforc.2003.nc" => "clmforc.2002.nc")

            # Confirm the _P rate-print dumps really are dead (justifies the
            # integrated-delta approach rather than a direct rate diff).
            da0 = NCDataset(joinpath(dumpdir, "pdump_after_hydrologydrainage_n$(nstep_dec).nc"), "r")
            rate_dead = true
            for v in ("GROSS_NMIN_VR_P","ACT_IMMOB_NH4_VR_P","F_NIT_VR_P",
                      "POT_F_NIT_VR_P","F_DENIT_VR_P","SMIN_NH4_TO_PLANT_VR_P")
                haskey(da0, v) && maximum(abs.(Float64.(da0[v][:, :]))) > 0.0 && (rate_dead = false)
            end
            close(da0)
            @info "Fortran _P rate-term dumps dead (all zero)?" rate_dead
            @test rate_dead   # if this ever flips, a direct rate parity becomes possible

            inst, bounds = run_one_parity_step!(nstep_dec; use_cn=true, dumpdir=dumpdir,
                use_hydrstress=true, use_luna=true,
                step_date=DateTime(2002,1,1) + Hour(nstep_dec - date_base),
                forcing_file=forcing2002)
            sns = inst.soilbiogeochem_nitrogenstate

            db = NCDataset(joinpath(dumpdir, "pdump_before_step_n$(nstep_dec).nc"), "r")
            da = NCDataset(joinpath(dumpdir, "pdump_after_hydrologydrainage_n$(nstep_dec).nc"), "r")
            nl = min(10, CLM.varpar.nlevdecomp)

            # per-layer net-delta rel err + top-NL ratio for one mineral-N pool
            function pool_delta(vn, jarr)
                fb = Float64.(db[vn][:, 1]); fa = Float64.(da[vn][:, 1])
                maxrel = 0.0; sF = 0.0; sJ = 0.0
                for j in 1:nl
                    dF = fa[j] - fb[j]; dJ = Float64(jarr[1, j]) - fb[j]
                    maxrel = max(maxrel, abs(dJ - dF) / (1.0 + max(abs(dF), abs(dJ))))
                    sF += dF; sJ += dJ
                end
                (maxrel, sF == 0.0 ? NaN : sJ / sF)
            end

            no3_rel, no3_ratio   = pool_delta("smin_no3_vr", sns.smin_no3_vr_col)
            nh4_rel, nh4_ratio   = pool_delta("smin_nh4_vr", sns.smin_nh4_vr_col)
            sminn_rel, sminn_ratio = pool_delta("sminn_vr",  sns.sminn_vr_col)
            close(db); close(da)

            @info "mineral-N integrated pool-delta parity" no3=(no3_rel, no3_ratio) nh4=(nh4_rel, nh4_ratio) sminn=(sminn_rel, sminn_ratio)

            # NO3 source terms (denit/leaching/nitrif-in/uptake) match cleanly.
            @test no3_rel   < 5e-4          # verified 1.4e-5 (≈35× headroom)
            @test abs(no3_ratio - 1.0) < 0.05   # verified 0.985

            # NH4 / sminn carry the KNOWN closed-GPP/uptake over-drain residual
            # (~1.19 / ~1.14). Bound it with headroom so a regression that blows
            # the source terms up still trips, but the documented residual passes.
            @test nh4_rel   < 5e-3          # verified 8.8e-4
            @test sminn_rel < 5e-3
            @test 1.0 < nh4_ratio   < 1.35  # verified 1.195
            @test 1.0 < sminn_ratio < 1.30  # verified 1.144
        end
    end
end
