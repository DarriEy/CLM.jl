# ==========================================================================
# test_run_clm.jl — standalone run harness (src/driver/run_clm.jl).
#
# Exercises run_clm! end-to-end: a SHORT driver run that
#   (1) writes a CLM-style h0 history file + a restart, then
#   (2) starts a SECOND run with finidat= that restart and confirms the
#       prognostic state was loaded (continuation).
#
# GATED on the Bow domain inputs (machine-local, like test_*_robustness.jl):
# if surfdata/params/forcing are absent the test is skipped so CI / other
# machines stay green. Runs in an isolated subprocess (clm_initialize! mutates
# module globals such as the rooting-profile config).
# ==========================================================================
using Test, CLM

const _BOW_ROOT = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data"
const _BOW_DIR  = joinpath(_BOW_ROOT, "domain_Bow_at_Banff_lumped")
const _BOW_FS   = joinpath(_BOW_DIR, "settings", "CLM", "parameters", "surfdata_clm.nc")
const _BOW_FP   = joinpath(_BOW_DIR, "settings", "CLM", "parameters", "clm5_params.nc")
const _BOW_FORC = joinpath(_BOW_DIR, "data", "forcing", "CLM_input", "clmforc.2002_2004.nc")
const _SNOWOPT  = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_optics_5bnd_c013122.nc"
const _SNOWAGE  = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_drdt_bst_fit_60_c070416.nc"

# Build a fresh (inst, bounds, filt, tm) + opened ForcingReader for Bow.
function _bow_setup(start_date)
    CLM.rooting_profile_config.rooting_profile_method_water  = CLM.JACKSON_1996_ROOT
    CLM.rooting_profile_config.rooting_profile_method_carbon = CLM.JACKSON_1996_ROOT
    (inst, bounds, filt, tm) = CLM.clm_initialize!(; fsurdat=_BOW_FS, paramfile=_BOW_FP,
        start_date=start_date, dtime=3600, use_cn=false, use_luna=false,
        use_bedrock=true, use_aquifer_layer=false, h2osfcflag=0,
        fsnowoptics=_SNOWOPT, fsnowaging=_SNOWAGE, int_snow_max=3113.2227)
    CLM.init_soil_hydrology_config(baseflow_scalar=0.0022119554)
    config = CLM.CLMDriverConfig(use_cn=false, use_aquifer_layer=false,
                                 use_hydrstress=false, use_luna=false)
    fr = CLM.ForcingReader()
    CLM.forcing_reader_init!(fr, _BOW_FORC)
    return (config, inst, bounds, filt, tm, fr)
end

# A column-level history tape. Each field is tagged `level="column"` so it
# writes on the native `column` dim. (HistoryTape now supports per-level dims,
# so mixed patch+column tapes also round-trip even when np≠nc — see
# default_history_tape — but this harness keeps a clean single-level tape.)
function _column_tape()
    tape = CLM.HistoryTape()
    CLM.hist_addfld!(tape, "TG", "K",
        inst -> inst.temperature.t_grnd_col; level="column",
        long_name="ground temperature")
    CLM.hist_addfld!(tape, "QRUNOFF", "mm/s",
        inst -> inst.water.waterfluxbulk_inst.wf.qflx_runoff_col;
        level="column", long_name="total liquid runoff")
    CLM.hist_addfld!(tape, "QFLX_EVAP_TOT", "mm/s",
        inst -> inst.water.waterfluxbulk_inst.wf.qflx_evap_tot_col;
        level="column", long_name="total evaporation")
    return tape
end

@testset "run_clm! standalone harness (gated)" begin
    have_inputs = isfile(_BOW_FS) && isfile(_BOW_FP) && isfile(_BOW_FORC) &&
                  isfile(_SNOWOPT) && isfile(_SNOWAGE)
    if !have_inputs
        @info "run_clm!: Bow inputs absent — skipping" _BOW_FS _BOW_FORC
        @test_skip have_inputs
    else
        tmp = mktempdir()
        hist1 = joinpath(tmp, "run1.clm2.h0.nc")
        rst1  = joinpath(tmp, "run1.clm2.r.nc")
        hist2 = joinpath(tmp, "run2.clm2.h0.nc")
        start = CLM.DateTime(2002, 6, 1)
        nsteps = 6

        # --- Run 1: cold start → h0 + restart ----------------------------------
        (config, inst, bounds, filt, tm, fr) = _bow_setup(start)
        tape = _column_tape()
        res1 = CLM.run_clm!(config, inst, bounds, filt, tm, fr;
                            nsteps=nsteps, dtime=3600.0, start_date=start,
                            hist_tape=tape, hist_interval=3, hist_path=hist1,
                            restart_path=rst1, use_cn=false)
        CLM.forcing_reader_close!(fr)

        @test res1.nsteps == nsteps
        @test res1.hist_path == hist1
        @test res1.restart_path == rst1
        @test isfile(hist1)
        @test isfile(rst1)

        # h0 sanity: time records + a finite field present.
        CLM.NCDatasets.NCDataset(hist1, "r") do ds
            @test haskey(ds, "time")
            @test ds.dim["time"] >= 2          # interval=3 over 6 steps → 2 records
            @test haskey(ds, "TG")
            tg = CLM.NCDatasets.Array(ds["TG"])
            @test any(isfinite, tg)            # at least some valid ground-temp samples
        end

        # Capture the final restarted state to verify it loads back.
        tg_end  = copy(inst.temperature.t_grnd_col)
        liq_end = copy(inst.water.waterstatebulk_inst.ws.h2osoi_liq_col)
        @test any(isfinite, tg_end)

        # --- Run 2: start FROM the restart (finidat=) → confirm state loaded ----
        (config2, inst2, bounds2, filt2, tm2, fr2) = _bow_setup(start)
        # Before loading, a fresh cold-start inst2 differs from run-1's end state.
        tg_cold = copy(inst2.temperature.t_grnd_col)

        res2 = CLM.run_clm!(config2, inst2, bounds2, filt2, tm2, fr2;
                            nsteps=2, dtime=3600.0, start_date=start,
                            hist_tape=_column_tape(), hist_interval=0,
                            hist_path=hist2, finidat=rst1, use_cn=false)
        CLM.forcing_reader_close!(fr2)

        @test res2.nsteps == 2
        @test isfile(hist2)                    # final-flush write (interval=0)

        # The restart load must have taken effect: immediately after read_restart!
        # (verified indirectly) the t_grnd that drove run 2 came from rst1, not the
        # cold IC. We re-read the restart into a third inst and compare to tg_end.
        (_, inst3, bounds3, _, _, fr3) = _bow_setup(start)
        CLM.forcing_reader_close!(fr3)
        CLM.read_restart!(inst3, rst1; bounds=bounds3, use_cn=false)
        tg_loaded  = inst3.temperature.t_grnd_col
        liq_loaded = inst3.water.waterstatebulk_inst.ws.h2osoi_liq_col

        # Loaded state equals run-1's written end state (round-trip), and differs
        # from the cold IC (so the load genuinely changed something).
        @test tg_loaded == tg_end
        @test liq_loaded == liq_end
        @test tg_loaded != tg_cold

        rm(tmp; recursive=true, force=true)
    end
end
