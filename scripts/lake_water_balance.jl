# =============================================================================
# Lake column water-balance harness.
#
# Runs the full CLM.jl driver on a pure-LAKE column (Bow surfdata with
# PCT_LAKE=100; test_inputs/lake/surfdata_lake100.nc) from a cold start with the
# real Bow 2003 hourly forcing, and reports the column water-balance error
# errh2o_col = (endwb - begwb) - (sources - sinks)*dtime for the lake column.
#
# WHY: lake_hydrology.jl used to STUB the ending water balance as endwb = begwb
# (ComputeWaterMassLake was unported). That stub zeroed the storage-change term
# in the lake residual runoff qflx_qrgwl:
#
#   qflx_qrgwl = rain + snow - evap - snwcp_ice - snwcp_discarded_*
#                - (endwb - begwb)/dtime + flood
#
# so qflx_qrgwl was missing the real snow/soil storage change, and the column
# water balance (which uses the REAL endwb recomputed later by
# end_water_column_balance!) was off by exactly that storage change. Now that
# compute_water_mass_lake! is called for real, qflx_qrgwl absorbs the true
# storage change and the lake balance closes to roundoff.
#
# The check here is a genuine closure test, NOT a tautology: it asserts BOTH
#   (a) endwb != begwb on at least one step (the storage really moves), and
#   (b) |errh2o_col| stays below the CLM balance threshold (1e-5 mm).
#
# Usage: julia +1.12 --project=. scripts/lake_water_balance.jl [NSTEPS]
# Programmatic: lake_water_balance_report(; nsteps=24) -> NamedTuple or `missing`
#               (missing = machine-local inputs absent -> caller should skip).
# =============================================================================
using CLM, NCDatasets, Dates, Printf

const LWB_LAKE    = joinpath(@__DIR__, "..", "test_inputs", "lake", "surfdata_lake100.nc")
const LWB_FP      = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc"
const LWB_FFORC   = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/data/forcing/CLM_input/clmforc.2003.nc"
const LWB_SNOWOPT = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_optics_5bnd_c013122.nc"
const LWB_SNOWAGE = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_drdt_bst_fit_60_c070416.nc"

"""
    lake_water_balance_report(; nsteps=24, verbose=false)

Cold-start the driver on a 100%-lake column and track the lake column's water
balance for `nsteps` hourly steps.

Returns `missing` if the machine-local inputs are absent, else a NamedTuple:
  `max_abs_errh2o`   — max |errh2o_col| on the lake column over the run [mm]
  `max_abs_dwb`      — max |endwb - begwb| over the run [mm] (storage really moved)
  `endwb_ne_begwb`   — true if endwb differed from begwb on some step
  `all_finite`       — true if begwb/endwb/errh2o stayed finite
"""
function lake_water_balance_report(; nsteps::Int = 24, verbose::Bool = false)
    for f in (LWB_LAKE, LWB_FP, LWB_FFORC, LWB_SNOWOPT, LWB_SNOWAGE)
        isfile(f) || return missing
    end

    base = DateTime(2003, 1, 1)
    (inst, bounds, filt, tm) = CLM.clm_initialize!(; fsurdat=LWB_LAKE, paramfile=LWB_FP,
        start_date=base, dtime=3600, use_cn=false, use_luna=false, use_bedrock=true,
        use_aquifer_layer=false, h2osfcflag=0, fsnowoptics=LWB_SNOWOPT,
        fsnowaging=LWB_SNOWAGE, int_snow_max=2000.0)
    ng, nc = bounds.endg, bounds.endc

    c_lake = findfirst(filt.lakec)
    c_lake === nothing && error("no lake column found (filt.lakec all false)")

    config  = CLM.CLMDriverConfig(use_cn=false, use_aquifer_layer=false,
                                  use_hydrstress=false, use_luna=false)
    filt_ia = CLM.clump_filter_inactive_and_active
    fr = CLM.ForcingReader(); CLM.forcing_reader_init!(fr, LWB_FFORC)
    tf = replace(LWB_FFORC, r"clmforc\.[^/]*\.nc$" => "topo_forcing.nc")
    if isfile(tf)
        dt = NCDataset(tf, "r")
        if haskey(dt, "TOPO")
            ft = Float64(dt["TOPO"][1])
            for g in 1:ng; inst.atm2lnd.forc_topo_grc[g] = ft; end
            for c in 1:nc; inst.topo.topo_col[c] = ft; end
        end
        close(dt)
    end

    wb = inst.water.waterbalancebulk_inst
    max_abs_errh2o = 0.0
    max_abs_dwb    = 0.0
    endwb_ne_begwb = false
    all_finite     = true

    verbose && @printf("%-5s %14s %14s %16s %13s\n",
                       "step", "begwb", "endwb", "endwb-begwb", "errh2o")
    for s in 1:nsteps
        calday = CLM.get_curr_calday(tm); (declin, _) = CLM.compute_orbital(calday)
        (declinm1, _) = CLM.compute_orbital(calday - 3600.0/CLM.SECSPDAY)
        nextsw = calday + 3600.0/CLM.SECSPDAY
        CLM.init_daylength!(inst.gridcell, declin, declinm1, CLM.ORB_OBLIQR_DEFAULT, 1:ng)
        CLM.advance_timestep!(tm)
        force_date = base + Hour(max(s - 1, 1))
        CLM.read_forcing_step!(fr, inst.atm2lnd, force_date, ng, nc)
        for g in 1:ng
            inst.atm2lnd.forc_hgt_u_grc[g] = 30.0
            inst.atm2lnd.forc_hgt_t_grc[g] = 30.0
            inst.atm2lnd.forc_hgt_q_grc[g] = 30.0
        end
        CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
        (_, mon, d, _) = CLM.get_curr_date(tm)
        CLM.clm_drv!(config, inst, filt, filt_ia, bounds, true, nextsw, declin, declin,
            CLM.ORB_OBLIQR_DEFAULT, false, false, "", false;
            nstep=tm.nstep, is_first_step=(s == 1),
            is_beg_curr_day=CLM.is_beg_curr_day(tm), is_end_curr_day=CLM.is_end_curr_day(tm),
            is_beg_curr_year=CLM.is_beg_curr_year(tm), dtime=3600.0, mon=mon, day=d,
            photosyns=inst.photosyns)

        b = Float64(wb.begwb_col[c_lake])
        e = Float64(wb.endwb_col[c_lake])
        r = Float64(wb.errh2o_col[c_lake])
        (isfinite(b) && isfinite(e) && isfinite(r)) || (all_finite = false)
        if isfinite(e) && isfinite(b)
            max_abs_dwb = max(max_abs_dwb, abs(e - b))
            (e != b) && (endwb_ne_begwb = true)
        end
        isfinite(r) && (max_abs_errh2o = max(max_abs_errh2o, abs(r)))
        verbose && @printf("%-5d %14.6f %14.6f %16.6e %13.3e\n", s, b, e, e - b, r)
    end
    CLM.forcing_reader_close!(fr)

    return (max_abs_errh2o = max_abs_errh2o, max_abs_dwb = max_abs_dwb,
            endwb_ne_begwb = endwb_ne_begwb, all_finite = all_finite)
end

"""
    lake_water_balance_ok(; nsteps=24, tol=1e-5)

Boolean gate for the test suite: `missing` if inputs absent, else true iff the
lake column's water balance closes (|errh2o| < tol), endwb is genuinely computed
(differs from begwb on some step) and everything stayed finite.
"""
function lake_water_balance_ok(; nsteps::Int = 24, tol::Real = 1.0e-5)
    rep = lake_water_balance_report(; nsteps=nsteps)
    rep === missing && return missing
    return rep.all_finite && rep.endwb_ne_begwb && rep.max_abs_errh2o < tol
end

if abspath(PROGRAM_FILE) == @__FILE__
    na = filter(a -> occursin(r"^\d+$", a), ARGS)
    rep = lake_water_balance_report(; nsteps = isempty(na) ? 24 : parse(Int, na[1]),
                                    verbose = true)
    if rep === missing
        println("\nInputs absent — skipped.")
    else
        @printf("\nmax |errh2o_col| (lake) = %.6e mm   [threshold %.1e]\n",
                rep.max_abs_errh2o, 1.0e-5)
        @printf("max |endwb - begwb|     = %.6e mm   (endwb genuinely computed: %s)\n",
                rep.max_abs_dwb, rep.endwb_ne_begwb)
        @printf("all finite              = %s\n", rep.all_finite)
    end
end
