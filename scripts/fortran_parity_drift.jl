# =============================================================================
# Multi-step CN/BGC drift harness (free-running).
#
# Single-step parity (fortran_parity_common.jl) re-injects the Fortran state
# before every step, so it measures per-step *translation* error. This harness
# instead injects the Fortran restart ONCE and free-runs N steps, advancing
# forcing + time each step, diffing the Julia trajectory against the Fortran
# per-step `after_hydrologydrainage` dumps. It shows how/where the single-step
# residuals COMPOUND (drift) and whether anything runs away.
#
# Reference dumps: a contiguous window of before_step_n<N> + after_hydrology-
# drainage_n<N> dumps (the summer window n1757845..n1757872, 28 steps).
#
# Findings (2026-06): drift is bounded — the temporary C/N buffer pools drift
# fastest (cpool ~2e-4->2.9e-2, xsmrpool ~8e-4->3.5e-2 over a day), the mineral
# N compounds to ~1% then plateaus, and the large long-memory pools (leafc,
# decomp C/N) stay ≤1e-4. The cpool/xsmrpool drift is a coupled feedback that
# integrates the known single-step residuals (FUN uptake, smin N-balance), not a
# C leak — the single-step cpool budget conserves C.
#
# Usage: julia +1.12 --project=. scripts/fortran_parity_drift.jl
# =============================================================================
include(joinpath(@__DIR__, "fortran_parity_common.jl"))

const DRIFT_SUM    = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_summer"
const DRIFT_AUT    = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_autumn"
const DRIFT_N0     = 1757845          # first step with a dump
const DRIFT_NSTEPS = 28               # contiguous dumps available
const DRIFT_BASE   = 1753153          # nstep -> date base (DateTime(2002,1,1)+Hour(n-base))

# AUTUMN (leaf-offset) window: 480 contiguous dumps, 2202-10-09 00:00 .. 10-28 23:00.
# The seasonal-deciduous offset fires at n1759931 and the ndays_off=15 ramp ends at
# n1760291, so this window free-runs Julia straight through the whole senescence ramp.
const DRIFT_AUT_N0     = 1759897
const DRIFT_AUT_NSTEPS = 480

_relmax(jl, fa; thr = 1e-9) = begin
    mx = 0.0
    @inbounds for k in eachindex(jl)
        j = jl[k]; f = fa[k]
        (isfinite(j) && isfinite(f)) || continue
        d = abs(f) > thr ? abs(j - f) / abs(f) : abs(j - f)
        mx = max(mx, d)
    end
    mx
end

function run_drift(; n0::Int = DRIFT_N0, nsteps::Int = DRIFT_NSTEPS,
                     sumdir::String = DRIFT_SUM)
    start_date = DateTime(2002, 1, 1) + Hour(n0 - DRIFT_BASE)
    ffile = replace(FFORCING, "clmforc.2003.nc" => "clmforc.2002.nc")
    (inst, bounds, filt, tm) = build_bow_inst(; dtime = 3600, start_date = start_date, use_cn = true, use_luna = true)
    # PHS+LUNA (Bow lnd_in): allocate LUNA vcmax/jmax fields so inject fills them,
    # then seed vegwp from the IC dump so the first PHS Newton solve starts finite.
    if isempty(inst.photosyns.vcmx25_z_patch)
        inst.photosyns.vcmx25_z_patch = fill(30.0, bounds.endp, CLM.NLEVCAN)
        inst.photosyns.jmx25_z_patch  = fill(60.0, bounds.endp, CLM.NLEVCAN)
    end
    icdump = joinpath(sumdir, "pdump_before_step_n$(n0).nc")
    inject_dump!(inst, bounds, icdump)
    let ds = NCDataset(icdump, "r")
        if haskey(ds, "vegwp")
            vw = ds["vegwp"][:, :]
            for pd in 1:size(vw, 2), seg in 1:4
                inst.canopystate.vegwp_patch[pd, seg] = Float64(vw[seg, pd])
            end
        end
        close(ds)
    end

    config = CLM.CLMDriverConfig(use_cn = true, use_aquifer_layer = false,
                                 use_hydrstress = true, use_luna = true)
    filt_ia = CLM.clump_filter_inactive_and_active
    ng, nc, np = bounds.endg, bounds.endc, bounds.endp
    fr = CLM.ForcingReader(); CLM.forcing_reader_init!(fr, ffile)
    tf = replace(ffile, r"clmforc\.[^/]*\.nc$" => "topo_forcing.nc")
    if isfile(tf)
        dt = NCDataset(tf, "r")
        if haskey(dt, "TOPO")
            ft = Float64(dt["TOPO"][1])
            for g in 1:ng; inst.atm2lnd.forc_topo_grc[g] = ft; end
            for c in 1:nc; inst.topo.topo_col[c] = ft; end
        end
        close(dt)
    end

    sbns = inst.soilbiogeochem_nitrogenstate; sbcs = inst.soilbiogeochem_carbonstate
    cs = inst.bgc_vegetation.cnveg_carbonstate_inst
    cvs = inst.bgc_vegetation.cnveg_state_inst
    println("step | sminn_vr | smin_nh4 | leafc | soil1c | cpool(tree) | leafc J/F (grass) | offFlag J/F")
    nl = CLM.varpar.nlevdecomp
    report_every = nsteps > 60 ? 24 : 1   # long (autumn) window: report daily
    # Seed dayl/prev_dayl ONCE, then let clm_drv!'s update_daylength! maintain them —
    # which is exactly what the production driver does (clm_run.jl never calls
    # init_daylength! inside the loop).
    #
    # Calling init_daylength! EVERY step is actively wrong: update_daylength! (inside
    # clm_drv!, when !is_first_step) shifts prev_dayl = dayl and then recomputes dayl
    # from the same declin. If init_daylength! has just set dayl = daylength(declin),
    # the shift makes prev_dayl == dayl, so ws_flag = (dayl >= prev_dayl) is stuck at 1
    # and the seasonal-deciduous offset test (ws_flag == 0 && dayl < crit_dayl) can
    # NEVER fire — the leaf-offset ramp never starts in a free run.
    let calday0 = CLM.get_curr_calday(tm)
        (d0, _)  = CLM.compute_orbital(calday0)
        (dm1, _) = CLM.compute_orbital(calday0 - 3600.0 / CLM.SECSPDAY)
        CLM.init_daylength!(inst.gridcell, d0, dm1, CLM.ORB_OBLIQR_DEFAULT, 1:bounds.endg)
    end

    for i in 0:(nsteps - 1)
        cur = start_date + Hour(i)
        # CTSM's get_curr_calday() is the calday at the END of the timestep, and
        # UpdateDaylength runs on that (clm_run.jl does the same: advance, then
        # get_curr_calday).
        CLM.advance_timestep!(tm)
        calday = CLM.get_curr_calday(tm); (declin, _) = CLM.compute_orbital(calday)
        # post-advance calday IS the next radiation time (end of this step) — do not add dtime
        nextsw = calday
        CLM.read_forcing_step!(fr, inst.atm2lnd, cur, ng, nc)
        CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
        (yr, mon, d, tod) = CLM.get_curr_date(tm)
        CLM.clm_drv!(config, inst, filt, filt_ia, bounds, true, nextsw, declin, declin,
            CLM.ORB_OBLIQR_DEFAULT, false, false, "", false;
            nstep = tm.nstep, is_first_step = false,
            is_beg_curr_day = CLM.is_beg_curr_day(tm), is_end_curr_day = CLM.is_end_curr_day(tm),
            is_beg_curr_year = CLM.is_beg_curr_year(tm), dtime = 3600.0, mon = mon, day = d,
            photosyns = inst.photosyns)
        n = n0 + i
        da = NCDataset(joinpath(sumdir, "pdump_after_hydrologydrainage_n$(n).nc"), "r")
        msm = _relmax(sbns.sminn_vr_col[1, 1:nl], Float64.(da["sminn_vr"][1:nl, 1]))
        mnh = _relmax(sbns.smin_nh4_vr_col[1, 1:nl], Float64.(da["smin_nh4_vr"][1:nl, 1]))
        mlc = _relmax(cs.leafc_patch[2:3], Float64.(da["leafc"][2:3]))
        ms1 = _relmax(sbcs.decomp_cpools_vr_col[1, 1:nl, 4], Float64.(da["soil1c_vr"][1:nl, 1]))
        f_lc = Float64(da["leafc"][3]); f_of = Float64(da["offset_flag"][3])
        close(da)
        j_lc = Float64(cs.leafc_patch[3]); j_of = Float64(cvs.offset_flag_patch[3])
        if (i % report_every == 0) || i == nsteps - 1
            println("n", lpad(i + 1, 3), " | ", rpad(round(msm, sigdigits = 3), 9), " | ",
                    rpad(round(mnh, sigdigits = 3), 9), " | ", rpad(round(mlc, sigdigits = 3), 9),
                    " | ", rpad(round(ms1, sigdigits = 3), 9), " | ",
                    rpad(round(cs.cpool_patch[2], sigdigits = 4), 10), " | ",
                    rpad(string(round(j_lc, digits = 3), "/", round(f_lc, digits = 3)), 17), " | ",
                    Int(j_of), "/", Int(f_of), j_of != f_of ? "  <<< PHENOLOGY MISMATCH" : "")
        end
    end
    CLM.forcing_reader_close!(fr)
    return inst, bounds
end

if !isempty(ARGS) && lowercase(ARGS[1]) == "autumn"
    println("=== AUTUMN (leaf-offset) CN drift: inject once at n$DRIFT_AUT_N0, free-run $DRIFT_AUT_NSTEPS steps ===")
    println("    Fortran offset trigger n1759931 (10-10 10:00), ndays_off=15 ramp ends n1760291 (10-25 10:00).\n")
    run_drift(; n0 = DRIFT_AUT_N0, nsteps = DRIFT_AUT_NSTEPS, sumdir = DRIFT_AUT)
else
    println("=== SUMMER CN drift: inject once at n$DRIFT_N0, free-run $DRIFT_NSTEPS steps ===\n")
    run_drift()
end
