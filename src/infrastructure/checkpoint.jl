# ==========================================================================
# checkpoint.jl — complete, bit-exact-resumable model checkpoints.
#
# The NetCDF restart (restart_io.jl: write_restart/read_restart!) is a CURATED,
# Fortran-compatible snapshot — it persists the prognostic registry but OMITS the
# running-mean accumulators, snow-grain/melt diagnostics, and the many small
# forward-feeding fields a bit-exact continuation needs. (Diagnosed via the
# validation harness: a restart-file-only continuation diverges ~1e-4 in t_grnd,
# and the omitted state is distributed across nearly every sub-instance — no small
# curated set closes it; only the COMPLETE numeric state does.)
#
# This module provides the complementary artifact: a COMPLETE checkpoint that
# serializes the entire instance + time-manager cursor, so a resumed run reproduces
# an uninterrupted one BIT-EXACTLY. It is Julia-native (Serialization), not portable
# to Fortran — use write_restart for interop, write_checkpoint for exact resume.
# ==========================================================================

using Serialization

const CHECKPOINT_VERSION = 1

"""
    write_checkpoint(path, inst, tm) -> path

Write a COMPLETE, bit-exact-resumable checkpoint of the model state to `path`. Captures
the entire instance (every prognostic field + running-mean accumulator + snow-grain/melt
diagnostic + …) plus the time-manager cursor (nstep + current_date), so [`read_checkpoint`]
+ continuation reproduces an uninterrupted run BIT-EXACTLY. The forcing reader is NOT
stored — it self-seeks by date once the restored date seeds the TimeManager.

Contrast write_restart (restart_io.jl): that is the curated, Fortran-interoperable NetCDF
restart; this is the Julia-native exact-resume checkpoint.
"""
function write_checkpoint(path::AbstractString, inst, tm)
    open(path, "w") do io
        serialize(io, (CHECKPOINT_VERSION, inst, tm.nstep, tm.current_date))
    end
    return path
end

"""
    read_checkpoint(path) -> (inst, nstep, current_date)

Restore a complete checkpoint written by [`write_checkpoint`]. Returns the deserialized
instance plus the timestep counter and date to seed a fresh TimeManager. Pair the returned
instance with the run's bounds/filter/config and a fresh forcing reader to continue.
"""
function read_checkpoint(path::AbstractString)
    (ver, inst, nstep, date) = open(deserialize, path)
    ver == CHECKPOINT_VERSION ||
        error("checkpoint version mismatch: file=$ver, expected $CHECKPOINT_VERSION")
    return (inst, nstep, date)
end
