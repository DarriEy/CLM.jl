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
    have = isfile(script) && all(isfile(joinpath(dumpdir, "bgcdump_after_ch4_n$(4721 + i).nc")) &&
                                 isfile(joinpath(dumpdir, "pdump_before_step_n$(4721 + i).nc"))
                                 for i in 0:(nsteps - 1))
    if !have
        @info "CH4-as-a-source parity test SKIPPED (bow_ref_ch4tws dumps not present)"
        @test_skip have
    else
        Base.include(@__MODULE__, script)
        run_source = Base.invokelatest(getglobal, @__MODULE__, :run_source)
        r = Base.invokelatest(run_source; nsteps = nsteps, mode = :locked)
        w = r.worst

        # CN confound: Julia's hr_vr vs Fortran's, recovered from CH4_PROD_SAT.
        @test r.hrworst < 5e-3
        # The source term and everything it feeds, CH4 free-running for 12 h.
        @test w["CH4_PROD_SAT"]      < 5e-3
        @test w["CH4_OXID_SAT"]      < 5e-3
        @test w["CH4_AERE_SAT"]      < 5e-3
        @test w["CONC_CH4_SAT"]      < 5e-3
        @test w["CONC_O2_SAT"]       < 5e-3
        @test w["O2_DECOMP_SAT"]     < 4e-3   # 1.3e-3 after the rr livecroot_mr fix
        @test w["O2_DECOMP_UNSAT"]   < 4e-3
        @test w["CH4_SURF_AERE_SAT"] < 5e-3
        @test w["CH4_SURF_FLUX_TOT"] < 5e-3
        @test w["TOTCOLCH4"]         < 5e-3
        @test w["ZWT_CH4_UNSAT"]     < 1e-6

        # O2-demand components (docs/CH4_FIRE_PARITY.md §17), only in bow_ref_ch4o2.
        if haskey(w, "O2_DEMAND_RR")
            @test w["RR"]             < 2e-3   # root resp. incl. livecroot_mr (was ~1.7e-3 at night)
            @test w["LIVECROOT_MR"]   < 1e-2
            @test w["O2_DEMAND_RR"]   < 2e-3
            @test w["O2_DEMAND_NIT"]  < 1e-4
            @test w["POT_F_NIT_VR"]   < 1e-4
            @test w["O_SCALAR"]       < 1e-12
            @test w["O2_OXID_SAT"]    < 1e-2
        end
    end
end
