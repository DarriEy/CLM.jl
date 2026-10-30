# ==========================================================================
# test_fortran_parity_ch4_source.jl — gated CH4-as-a-source regression.
#
# Drives scripts/fortran_parity_ch4_source.jl (docs/CH4_FIRE_PARITY.md §16): Bow
# TWS_inversion reference (finundated ~0.045), CN + biogeophysics re-injected from
# the Fortran restart every step (the CN confound removed), CH4 state injected ONCE
# and then free-running in Julia. It guards the fixes that closed the source term:
#   * read_ch4_params! (q10ch4 1.33 from the paramfile, was the 1.5 default) — the
#     whole former ~7.9% CH4_PROD_SAT residual;
#   * the use_lch4 CN coupling (anoxia o_scalar, nitrif anaerobic_frac / O2 demand);
#   * watfc (field capacity) in cold start;
#   * the zi(1:nlevsoi) slice into ch4! (zwt_ch4_unsat / ebullition head).
# Tolerances carry ~3-5x headroom over the verified numbers.
#
# GATED on the bow_ref_ch4tws dumps; skipped (not failed) when absent.
# ==========================================================================

using Test
using CLM
using NCDatasets
using Dates

include(joinpath(@__DIR__, "testdata.jl"))

@testset "CH4-as-a-source parity (gated)" begin
    script  = joinpath(@__DIR__, "..", "scripts", "fortran_parity_ch4_source.jl")
    # bow_ref_ch4o2 (the O2-demand bgcdump rebuild, bit-identical to bow_ref_ch4tws in
    # every shared field) is what the harness reads when present.
    o2dir   = symfluence_path("clm_bgc_spinup", "bow_ref_ch4o2")
    dumpdir = isdir(o2dir) ? o2dir : symfluence_path("clm_bgc_spinup", "bow_ref_ch4tws")
    nsteps  = 12
    # The harness also reads two CESM_INPUTDATA files: the uniform-ndep parity stream
    # and the TWS_inversion finundation coefficients. Skip (not error) without them.
    ndep    = joinpath(cesm_inputdata_root(), "lnd", "clm2", "ndepdata",
                       "fndep_clm_UNIFORM1e-8_0.9x1.25_yr2000_CLMjl-parity.nc")
    invfile = joinpath(cesm_inputdata_root(), "lnd", "clm2", "paramdata",
                       "finundated_inversiondata_0.9x1.25_c170706.nc")
    have_dumps = isfile(script) && all(isfile(joinpath(dumpdir, "bgcdump_after_ch4_n$(4721 + i).nc")) &&
                                       isfile(joinpath(dumpdir, "pdump_before_step_n$(4721 + i).nc"))
                                       for i in 0:(nsteps - 1))
    have = have_dumps && isfile(ndep) && isfile(invfile)
    if !have
        missing_what = !have_dumps ? "bow_ref_ch4o2/ch4tws dumps" :
                       !isfile(ndep) ? "ndep file $(ndep)" : "inversion file $(invfile)"
        @info "CH4-as-a-source parity test SKIPPED ($(missing_what) not present)"
        @test_skip have
    else
        Base.include(@__MODULE__, script)
        run_source = Base.invokelatest(getglobal, @__MODULE__, :run_source)
        r = Base.invokelatest(run_source; nsteps = nsteps, mode = :locked)
        w = r.worst

        # CN confound: Julia's hr_vr vs Fortran's, recovered from CH4_PROD_SAT. Was a
        # constant 6.6e-4 = 365.2425/365 - 1 until the harness ran the reference's
        # GREGORIAN calendar (docs/CH4_FIRE_PARITY.md §18); 1.9e-6 over 24 h now.
        @test r.hrworst < 1e-5
        # The source term and everything it feeds, CH4 free-running for 12 h.
        # Verified 24-h worsts in the comments; gates carry ~3-5x headroom.
        @test w["CH4_PROD_SAT"]      < 1e-5    # 1.9e-6
        @test w["CH4_OXID_SAT"]      < 2e-4    # 3.9e-5
        @test w["CH4_AERE_SAT"]      < 5e-5    # 1.1e-5
        @test w["CONC_CH4_SAT"]      < 1e-5    # 2.6e-6
        @test w["CONC_O2_SAT"]       < 1.5e-3  # 3.2e-4; layer-1 aerenchyma sensitivity (§17)
        @test w["O2_DECOMP_SAT"]     < 3e-4    # 6.5e-5 (was 1.3e-3 before lnc_opt + calendar)
        @test w["O2_DECOMP_UNSAT"]   < 3e-4
        @test w["CH4_SURF_AERE_SAT"] < 2e-5    # 4.6e-6
        @test w["CH4_SURF_FLUX_TOT"] < 5e-5    # 1.1e-5
        @test w["TOTCOLCH4"]         < 5e-6    # 1.0e-6
        @test w["ZWT_CH4_UNSAT"]     < 1e-6

        # O2-demand components (docs/CH4_FIRE_PARITY.md §17-18), only in bow_ref_ch4o2.
        if haskey(w, "O2_DEMAND_RR")
            # Day-time RR growth respiration: 3.6e-3 -> 6.9e-5 once the harness runs
            # the reference's lnc_opt=.true. (leaf MR from prognostic leafn).
            @test w["RR"]             < 3e-4
            @test w["FROOT_MR"]       < 1e-8   # 1.3e-10
            @test w["LIVECROOT_MR"]   < 1e-10  # 7.5e-13
            @test w["O2_DEMAND_RR"]   < 3e-4   # 6.5e-5
            @test w["O2_DEMAND_NIT"]  < 1e-5   # 8.1e-7
            @test w["POT_F_NIT_VR"]   < 1e-5
            @test w["O_SCALAR"]       < 1e-12
            @test w["O2_OXID_SAT"]    < 2e-4   # 3.9e-5
        end
    end
end
