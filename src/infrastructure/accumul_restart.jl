# ==========================================================================
# Ported from: src/main/accumulMod.F90 (accumulRest) + the InitAccVars restore
# path in atm2lndType.F90 / Wateratm2lndBulkType.F90.
#
# accumulRest persists, for every accumulator field, BOTH its running value AND
# its per-field step counter (`nsteps`); InitAccVars restores them so the running
# means are warm from the first step of a restart run. In CLM.jl the always-on
# per-timestep runmean accumulators (RH30/PREC10/PREC30/PREC60, RH24/PREC24, …)
# are evaluated IN-ARRAY on the Atm2LndData patch fields (see `accum_runmean` /
# `atm2lnd_update_acc_vars!` in types/atm2lnd.jl). For those fields:
#
#   * the VALUE lives directly in the patch vector (a2l.rh30_patch, …), and
#   * the per-point step COUNTER is NOT a separate array — it is `min(nstep, win)`
#     derived from the global model timestep (valid because the driver updates
#     every accumulator on every step from nstep=1; the design comment in
#     accumul.jl spells this out).
#
# So the faithful restore of accumulRest for these fields is:
#   (1) copy the persisted VALUE into the patch vector, and
#   (2) run the model at the true (warm) global `nstep` — which makes
#       `accum_runmean`'s `min(nstep, win)` equal the converged window, exactly as
#       the persisted `nsteps` records. No separate counter array is needed;
#       `restore_summary.warm_nstep` returns the counter so the caller can set the
#       time manager (a fresh restart that starts at nstep=1 would otherwise let
#       `accum_runmean` RESET the mean to the instantaneous sample on step 1).
#
# `read_fortran_restart!` maps ~170 prognostic restart variables but NONE of the
# `accumul` fields; this closes that gap for the runmean forcing accumulators.
# The Fortran dump writes them as `<NAME>_VALUE` / `<NAME>_NSTEPS`.
# ==========================================================================

# (dump VALUE name, dump NSTEPS name, Atm2LndData field, accumulation window [days])
const _ATM2LND_RUNMEAN_ACCUM = (
    ("RH30_VALUE",   "RH30_NSTEPS",   :rh30_patch,   30),
    ("PREC10_VALUE", "PREC10_NSTEPS", :prec10_patch, 10),
    ("PREC30_VALUE", "PREC30_NSTEPS", :prec30_patch, 30),
    ("PREC60_VALUE", "PREC60_NSTEPS", :prec60_patch, 60),
    ("RH24_VALUE",   "RH24_NSTEPS",   :rh24_patch,    1),
    ("PREC24_VALUE", "PREC24_NSTEPS", :prec24_patch,  1),
)

"""
    restore_atm2lnd_runmean_accum!(a2l, dumpfile, bounds; dtime=1800, verbose=false)
      -> NamedTuple

Restore the `atm2lnd` / `wateratm2lndbulk` runmean forcing accumulators
(`RH30`, `PREC10`, `PREC30`, `PREC60`, and — when present — `RH24`, `PREC24`)
from a Fortran restart / dump file, reproducing `accumulRest`'s restore.

For each accumulator whose `<NAME>_VALUE` is present in `dumpfile` AND whose
target `Atm2LndData` patch vector is allocated long enough, the persisted value
is copied into the patch vector. The persisted `<NAME>_NSTEPS` counter is read
and returned; the largest one is reported as `warm_nstep`.

Returns `(; restored, warm_nstep, nsteps)` where
- `restored::Vector{Symbol}` — the Atm2LndData fields that were filled,
- `warm_nstep::Int`          — `max` persisted counter (0 if none): the global
  `nstep` a restart must run at for `accum_runmean` to CONTINUE these running
  means instead of resetting them to the instantaneous sample (see file header),
- `nsteps::Dict{Symbol,Int}` — the persisted per-field counter.

Only touches the listed forcing accumulators; a cold start (no `<NAME>_VALUE` in
the file) leaves every field untouched.
"""
function restore_atm2lnd_runmean_accum!(a2l::Atm2LndData, dumpfile::AbstractString,
                                        bounds; dtime::Int = 1800, verbose::Bool = false)
    restored = Symbol[]
    nsteps = Dict{Symbol,Int}()
    warm_nstep = 0
    endp = _accum_bounds_endp(bounds)

    NCDataset(dumpfile, "r") do ds
        for (vname, nname, field, windays) in _ATM2LND_RUNMEAN_ACCUM
            haskey(ds, vname) || continue
            dst = getfield(a2l, field)
            (dst isa AbstractVector && length(dst) > 0) || continue

            raw = ds[vname][:]
            np = min(endp, length(dst), length(raw))
            filled = false
            @inbounds for p in 1:np
                x = raw[p]
                ismissing(x) && continue
                dst[p] = Float64(x)
                filled = true
            end
            filled || continue
            push!(restored, field)

            # Persisted step counter (accumulRest writes it alongside the value).
            if haskey(ds, nname)
                nraw = ds[nname][:]
                ns = 0
                for v in nraw
                    ismissing(v) || (ns = max(ns, Int(v)))
                end
                # A converged runmean should carry nsteps == its window in steps.
                win = accum_window_steps(windays, dtime)
                nsteps[field] = ns == 0 ? win : ns
            else
                nsteps[field] = accum_window_steps(windays, dtime)
            end
            warm_nstep = max(warm_nstep, nsteps[field])
        end
    end

    if verbose
        @info "restore_atm2lnd_runmean_accum!" restored warm_nstep
    end
    return (; restored, warm_nstep, nsteps)
end

# bounds may be a NamedTuple/struct with `endp`, or a UnitRange of patches.
_accum_bounds_endp(bounds) = hasproperty(bounds, :endp) ? bounds.endp :
    (bounds isa AbstractUnitRange ? last(bounds) : error(
        "restore_atm2lnd_runmean_accum!: cannot infer patch count from bounds::$(typeof(bounds))"))
