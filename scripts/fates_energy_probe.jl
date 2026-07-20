#!/usr/bin/env julia
# ==========================================================================
# fates_energy_probe.jl — bound the multi-year FATES longwave failure (#277 §4).
#
# Reuses `fates_longhorizon.jl`'s `build()` (identical cold start, identical
# flags — including the INHERITED `use_bedrock`, which is the whole point) and
# steps the driver, watching `errlon_patch` every step. At the FIRST step whose
# |errlon| crosses a threshold it dumps, per HLM patch, everything needed to
# attribute the imbalance:
#
#   active / wtcol / in-filter membership / elai+esai / t_veg
#   eflx_lwrad_out, eflx_lwrad_net, forc_lwrad, and errlon itself
#
# The load-bearing column is **in-filter**: the BitVector filters are built ONCE
# by `set_filters!` at initialisation, while `fates_set_filters!` flips
# `pch.active[p] = true` on a patch the day FATES's disturbance population grows.
# A patch that is `active` but NOT in the flux filters is skipped by every
# surface-energy routine and still audited by the balance check.
#
#   FATES_NDAYS=120 julia +1.12 --project=. scripts/fates_energy_probe.jl
#
# Env: all of fates_longhorizon.jl's, plus
#   PROBE_THRESH  (W/m2, default 1.0)   first |errlon| that triggers the dump
#   PROBE_NDUMP   (default 3)           how many triggering steps to dump
# ==========================================================================
Base.include(@__MODULE__, joinpath(@__DIR__, "..", "test", "testdata.jl"))
using CLM, Printf, Dates
const _C = CLM

const LH = Module(:LH)
Base.include(LH, joinpath(@__DIR__, "fates_longhorizon.jl"))

function inflt(mask, p)
    isempty(mask) && return "—"
    p <= length(mask) || return "oob"
    mask[p] ? "Y" : "."
end

function dump_patches(inst, filt, bounds, i, day, tag)
    pch = inst.patch; ef = inst.energyflux; a2l = inst.atm2lnd
    cs = inst.canopystate; tp = inst.temperature
    @printf("\n=== %s  step=%d day=%d ===\n", tag, i, day)
    @printf("%4s %3s %6s %9s  %-3s%-3s%-3s%-3s%-3s  %8s %8s %10s %10s %10s %10s\n",
        "p", "c", "act", "wtcol", "nlk","sol","snc","exv","nxv",
        "elai", "esai", "t_veg", "lw_out", "lw_net", "errlon")
    for p in 1:bounds.endp
        c = pch.column[p]
        @printf("%4d %3d %6s %9.5f  %-3s%-3s%-3s%-3s%-3s  %8.4f %8.4f %10.3f %10.3f %10.3f %10.3f\n",
            p, c, string(pch.active[p]), pch.wtcol[p],
            inflt(filt.nolakep, p), inflt(filt.soilp, p),
            inflt(filt.soilnopcropp, p),
            inflt(filt.exposedvegp, p), inflt(filt.noexposedvegp, p),
            cs.elai_patch[p], cs.esai_patch[p], tp.t_veg_patch[p],
            ef.eflx_lwrad_out_patch[p], ef.eflx_lwrad_net_patch[p], ef.errlon_patch[p])
    end
    @printf("forc_lwrad(col1) = %.3f    (errlon == -forc_lwrad  =>  the patch's LW was never computed)\n",
        a2l.forc_lwrad_downscaled_col[1])
    flush(stdout)
end

function main()
    ndays  = parse(Int, get(ENV, "FATES_NDAYS", "120"))
    thresh = parse(Float64, get(ENV, "PROBE_THRESH", "1.0"))
    ndump  = parse(Int, get(ENV, "PROBE_NDUMP", "3"))
    println("="^70)
    println("  FATES longwave-imbalance probe — $ndays days, thresh=$thresh W/m2")
    println("  use_bedrock env = '", get(ENV, "FATES_USE_BEDROCK", "(inherit)"), "'")
    println("="^70)

    inst, fates, config, bounds, filt, filt_ia = LH.build()
    inst.balcheck.hard_error = false            # probe: observe, never abort
    site = fates.sites[1]; photosyns = inst.photosyns
    dtime = 1800.0; nsteps = Int(round(86400/dtime))*ndays
    forcing = get(ENV, "FATES_FORCING",
        "$(LH.DATA)/domain_Aripuana_Amazon/data/forcing/CLM_input/clmforc.2004.nc")
    fyr = parse(Int, get(ENV, "FATES_YEAR", "2004"))
    fr = _C.ForcingReader(); _C.forcing_reader_init!(fr, forcing); fr.interp_time = true
    start_date = DateTime(fyr, 1, 1)

    cen0 = LH.census(inst)
    @printf("  cold start: ncoh=%d npatch=%d  bounds.endp=%d\n", cen0.ncoh, cen0.npatch, bounds.endp)
    @printf("  %6s %6s %7s %12s %12s\n", "day", "ncoh", "npatch", "max|errlon|", "argmax p")
    println("  " * "-"^52)

    dumped = 0; day = 0; prev_npatch = cen0.npatch
    for i in 1:nsteps
        step_start = start_date + Second((i-1)*Int(dtime))
        is_beg = (Dates.hour(step_start)==0 && Dates.minute(step_start)==0 && Dates.second(step_start)==0)
        yr_off = Dates.year(step_start) - fyr
        forcing_time = yr_off == 0 ? step_start : step_start - Dates.Year(yr_off)
        _C.read_forcing_step!(fr, inst.atm2lnd, forcing_time, 1, 1;
            gridcell_latdeg=inst.gridcell.latdeg, gridcell_londeg=inst.gridcell.londeg, dtime=Int(dtime))
        inst.atm2lnd.forc_topo_grc[1]=200.0; inst.topo.topo_col[1]=200.0
        _C.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
        sod = Dates.hour(step_start)*3600 + Dates.minute(step_start)*60
        calday = Dates.dayofyear(step_start) + sod/86400.0
        (declin, _e) = _C.compute_orbital(calday); nextsw_cday = calday + Int(dtime)/86400.0
        try
            _C.clm_drv!(config, inst, filt, filt_ia, bounds, true, nextsw_cday, declin, declin, 0.4091,
                false, false, "20260101", false; nstep=i, is_first_step=(i==1), is_beg_curr_day=is_beg,
                is_end_curr_day=false, is_beg_curr_year=false, dtime=dtime, mon=Dates.month(step_start),
                day=Dates.day(step_start), secs=sod, jday=Dates.dayofyear(step_start), photosyns=photosyns)
        catch e
            @printf("  ✗ ERROR at step %d (day %d): %s\n", i, day, first(split(sprint(showerror, e), "\n")))
            dump_patches(inst, filt, bounds, i, day, "AT ERROR")
            break
        end

        ef = inst.energyflux
        best = 0.0; bestp = 0
        for p in 1:bounds.endp
            v = ef.errlon_patch[p]
            (isfinite(v) && abs(v) < 1e30) || continue
            abs(v) > best && (best = abs(v); bestp = p)
        end
        if best > thresh && dumped < ndump
            dumped += 1
            dump_patches(inst, filt, bounds, i, day, "FIRST |errlon| > $thresh  (#$dumped)")
        end

        if is_beg && i > 1
            day = Dates.dayofyear(step_start) - 1
            cen = LH.census(inst)
            if cen.npatch != prev_npatch
                @printf("  >>> day %d: FATES patch count %d -> %d\n", day, prev_npatch, cen.npatch)
                prev_npatch = cen.npatch
            end
            @printf("  %6d %6d %7d %12.4g %12d\n", day, cen.ncoh, cen.npatch, best, bestp)
            flush(stdout)
        end
    end
    _C.forcing_reader_close!(fr)
    println("\n  probe done (dumped $dumped step(s))")
end

main()
