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

const DRIFT_SUM    = joinpath(dirname(DUMPDIR), "bgc_ref_summer")
const DRIFT_N0     = 1757845          # first step with a dump
const DRIFT_NSTEPS = 28               # contiguous dumps available
const DRIFT_BASE   = 1753153          # nstep -> date base (DateTime(2002,1,1)+Hour(n-base))

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
    (inst, bounds, filt, tm) = build_bow_inst(; dtime = 3600, start_date = start_date, use_cn = true)
    inject_dump!(inst, bounds, joinpath(sumdir, "pdump_before_step_n$(n0).nc"))

    config = CLM.CLMDriverConfig(use_cn = true, use_aquifer_layer = false)
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
    println("step | sminn_vr | smin_nh4 | leafc | soil1c | cpool(tree)")
    nl = CLM.varpar.nlevdecomp
    for i in 0:(nsteps - 1)
        cur = start_date + Hour(i)
        calday = CLM.get_curr_calday(tm); (declin, _) = CLM.compute_orbital(calday)
        nextsw = calday + 3600.0 / CLM.SECSPDAY
        (declinm1, _) = CLM.compute_orbital(calday - 3600.0 / CLM.SECSPDAY)
        CLM.init_daylength!(inst.gridcell, declin, declinm1, CLM.ORB_OBLIQR_DEFAULT, 1:bounds.endg)
        CLM.advance_timestep!(tm)
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
        close(da)
        println("n", lpad(i + 1, 2), " | ", rpad(round(msm, sigdigits = 3), 9), " | ",
                rpad(round(mnh, sigdigits = 3), 9), " | ", rpad(round(mlc, sigdigits = 3), 9),
                " | ", rpad(round(ms1, sigdigits = 3), 9), " | ", round(cs.cpool_patch[2], sigdigits = 4))
    end
    CLM.forcing_reader_close!(fr)
    return inst, bounds
end

run_drift()
