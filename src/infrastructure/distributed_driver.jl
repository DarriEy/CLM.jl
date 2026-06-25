# ==========================================================================
# distributed_driver.jl — MPI distributed-execution driver layer
#
# This is the final invasive piece of "Tier E: everything Fortran can do,
# including fully distributed / MPI". It maps CLM's per-clump physics loop onto
# MPI ranks (clumps -> ranks) instead of (or in addition to) OS threads
# (clumps -> threads), mirroring CTSM's two-level parallelism:
#
#   MPI over PROCESSES   (decompInitMod.F90 : pid = mod(n-1, npes))
#   OpenMP over CLUMPS   (clm_driver.F90    : do nc = 1, nclumps)
#
# Fortran call structure mirrored
# --------------------------------
#   * decompInit_lnd / _clumps / _glcp  (main/decompInitMod.F90):
#       round-robin gridcell->clump->rank assignment + per-rank proc bounds +
#       the gindex_* global-index maps. Already ported in decomp_init.jl; here we
#       wrap it as `decompInit_distributed!` which drives it per-rank from the
#       live `SPMDState` (iam/npes).
#   * the `do nc = 1, nclumps` driver loop (clm_driver.F90 / clm_drv): each rank
#       walks ONLY the clumps it owns (procinfo.cid lists exactly the local
#       rank's clumps), so `clm_run_distributed!` runs the existing
#       `clm_run_clump_physics!` over the rank-local clump bounds. Each rank
#       touches only its subdomain — no halo, no cross-rank reads in the physics.
#   * the gather-to-master serial-I/O path (spmdMod.F90 mpi_gatherv +
#       spmdGathScatMod, history/restart write on masterproc): `gather_to_master`
#       collects a decomposed subgrid field onto rank 0 in GLOBAL index order
#       (using the gindex_* maps), exactly as CLM gathers a begg:endg / begc:endc
#       vector for serial NetCDF output. `scatter_from_master` is the inverse.
#
# Single-rank invariant
# ---------------------
# When `npes == 1` (the default — no `mpiexec`), every routine here is a
# transparent passthrough: `clm_run_distributed!` runs the whole proc domain as
# before, `gather_to_master` returns the local field reordered into global order
# (which on one rank is just the field itself when gindex is the identity), and
# `scatter_from_master` returns the master's own block. The default driver and
# the full serial test suite are byte-identical and unaffected.
#
# Bit-identity invariant (the CTSM PE-layout guarantee)
# -----------------------------------------------------
# Because CLM land physics has NO horizontal coupling (no halo, no cross-column
# stencil — see spmd.jl header), a per-subgrid field computed on N ranks and
# gathered back into global index order is BIT-IDENTICAL to the same field
# computed serially on 1 rank over the same global domain. `gather_to_master`
# reassembles strictly by global index (`out[gindex[k]] = local[k]`), so the
# result is independent of how many ranks the domain was split across or what
# order their blocks arrive in. This is asserted by the 2-rank CI smoke.
# ==========================================================================

"""
    decompInit_distributed!(numg; clump_pproc=1, spmd::SPMDState=SPMD,
                            decomp_data::DecompData=decomp, kwargs...) -> DecompData

Initialize the domain decomposition for the *local* MPI rank, driven by the live
`SPMDState`. Thin wrapper over [`decompInit!`](@ref) that supplies `npes` and
`iam` from `spmd` (so the local rank's `procinfo`/bounds are the ones populated).

On a single rank (`npes == 1`) this is exactly `decompInit!(numg; ...)` with
`iam = 0` — byte-identical to the serial decomposition. Extra `kwargs` (the
per-gridcell subgrid-count vectors `ncols_per_g`, `npatches_per_g`, …) are
forwarded unchanged.

Mirrors the per-rank entry into `decompInit_lnd` in `decompInitMod.F90`, where
every process calls the same decomposition routine but records its own
`procinfo` via `iam`.
"""
function decompInit_distributed!(numg::Integer; clump_pproc::Integer = 1,
                                 spmd::SPMDState = SPMD,
                                 decomp_data::DecompData = decomp, kwargs...)
    npes = spmd.npes
    iam  = spmd.iam
    return decompInit!(numg; clump_pproc = clump_pproc, npes = npes, iam = iam,
                       decomp_data = decomp_data, kwargs...)
end

"""
    rank_clump_bounds(; decomp_data::DecompData=decomp) -> Vector{BoundsType}

The list of PROC-RELATIVE clump bounds owned by the local rank, in local clump
order. `procinfo.cid` already lists exactly the clumps this rank owns (the
round-robin assignment in `decompInit!` only records `iam`'s clumps), so this is
`[get_clump_bounds(n) for n in 1:get_proc_clumps()]`.

This is the rank-local analogue of the clump list the threaded driver builds; it
is what `clm_run_distributed!` walks. With `clump_pproc == 1` it is a single
whole-rank clump (== the rank's proc bounds), so the distributed path collapses
to the historical single-clump path on each rank.
"""
function rank_clump_bounds(; decomp_data::DecompData = decomp)
    nclumps_local = get_proc_clumps(decomp_data = decomp_data)
    return [get_clump_bounds(n; decomp_data = decomp_data) for n in 1:nclumps_local]
end

"""
    clm_run_distributed!(phys!, nsteps; decomp_data=decomp, spmd=SPMD,
                         threaded=false) -> Nothing

Run the rank-local distributed driver for `nsteps` timesteps. On each step the
clump-safe per-column physics `phys!(bounds_clump)` is applied over ONLY the
clumps owned by the local rank (via [`clm_run_clump_physics!`](@ref) on
[`rank_clump_bounds`](@ref)). Each rank therefore advances only its subdomain;
there is no cross-rank communication inside the step (CLM land physics has no
halo).

`phys!` follows the same contract as the threaded clump loop: it iterates
`bounds_clump.begc:bounds_clump.endc` (etc.), gates on absolute-indexed
proc-level filter masks, and writes only into the disjoint proc-level state slice
for those columns — no cross-clump / cross-rank reads or shared mutation.

`threaded` enables `Threads.@threads` over the rank's clumps (two-level
MPI+OpenMP parallelism, mirroring CTSM); the result is independent of clump count
and threading.

On a single rank this runs the whole proc domain `nsteps` times exactly as the
serial driver would. A barrier is taken at the end of each step under a real
multi-rank launch (no-op otherwise) so ranks stay loosely in lockstep for any
subsequent gather.

`phys!` may optionally accept the step index as a second argument
(`phys!(bounds_clump, step)`); detected by method existence.
"""
function clm_run_distributed!(phys!, nsteps::Integer;
                              decomp_data::DecompData = decomp,
                              spmd::SPMDState = SPMD,
                              threaded::Bool = false)
    clump_bounds = rank_clump_bounds(decomp_data = decomp_data)
    wants_step = hasmethod(phys!, Tuple{BoundsType, Int})
    for step in 1:nsteps
        if wants_step
            clm_run_clump_physics!(bc -> phys!(bc, step), clump_bounds; threaded = threaded)
        else
            clm_run_clump_physics!(phys!, clump_bounds; threaded = threaded)
        end
        spmd_barrier()
    end
    return nothing
end

# --------------------------------------------------------------------------
# Gather / scatter a decomposed subgrid field to / from the master rank, in
# GLOBAL index order. These mirror the gather-to-master serial-I/O path in
# CTSM (spmdMod mpi_gatherv + history/restart write on masterproc).
# --------------------------------------------------------------------------

# Per-rank gindex for a subgrid level (gridcell/column/patch/...). The local
# field passed to gather/scatter is indexed 1:length(gindex) in the SAME local
# (task) order the gindex lists.
function _level_gindex(level::Int; decomp_data::DecompData = decomp)
    return get_subgrid_level_gindex(level; decomp_data = decomp_data)
end

# Global size for a subgrid level (numg / numc / nump / ...).
function _level_gsize(level::Int; decomp_data::DecompData = decomp)
    return get_subgrid_level_gsize(level; decomp_data = decomp_data)
end

"""
    gather_to_master(local_field, level; decomp_data=decomp, spmd=SPMD, root=0) -> Vector

Gather a decomposed subgrid field onto the `root` (master) rank in GLOBAL index
order. `local_field` is this rank's values for its local points at subgrid
`level` (one of `SUBGRID_LEVEL_GRIDCELL` / `_COLUMN` / `_PATCH` / `_LANDUNIT` /
`_COHORT`), ordered exactly like that level's `gindex_*` map.

Returns, on `root`, a length-`numg`/`numc`/`nump`/… vector `out` with
`out[g] = local_field_of_the_rank_owning_g[local_slot_of_g]` for every global
index `g` (i.e. reassembled by `gindex`). On non-root ranks returns an empty
vector.

Bit-identity: the reassembly is purely `out[gindex[k]] = local[k]`, so the
gathered global field is independent of the rank/clump layout and is
BIT-IDENTICAL to the same field computed serially on one rank over the global
domain. This is the CTSM PE-layout invariant and the workhorse of the 2-rank
smoke.

Single rank: `gindex` is the identity (`1:numg`), so this returns `local_field`
placed into global order — equal to `local_field` itself when the rank owns the
whole domain in order.

Mirrors CTSM's `mpi_gatherv` of a begg:endg / begc:endc vector to the master for
serial NetCDF output, plus the gdc2glo reorder.
"""
function gather_to_master(local_field::AbstractVector{T}, level::Int;
                          decomp_data::DecompData = decomp,
                          spmd::SPMDState = SPMD, root::Int = 0) where {T}
    gindex = _level_gindex(level; decomp_data = decomp_data)
    gsize  = _level_gsize(level; decomp_data = decomp_data)
    length(local_field) == length(gindex) ||
        error("gather_to_master: local_field length $(length(local_field)) != local gindex length $(length(gindex)) at level $level")

    if !is_mpi_active()
        # Single rank: place the local field into global order via gindex.
        out = Vector{T}(undef, gsize)
        @inbounds for k in eachindex(local_field)
            out[gindex[k]] = local_field[k]
        end
        return out
    end

    # Multi-rank: gather the variable-length local blocks AND their gindex maps
    # to the master, then reassemble strictly by global index.
    npes = spmd.npes
    mylen = length(local_field)
    # Each rank reports its local length so root can build the gatherv counts.
    counts = _allgather_counts(mylen, npes, spmd)

    vals = spmd_gatherv(collect(local_field), counts; root = root)
    idxs = spmd_gatherv(collect(gindex), counts; root = root)

    if spmd.iam == root
        out = Vector{T}(undef, gsize)
        @inbounds for k in eachindex(idxs)
            out[idxs[k]] = vals[k]
        end
        return out
    else
        return T[]
    end
end

"""
    scatter_from_master(global_field, level; decomp_data=decomp, spmd=SPMD, root=0) -> Vector

Inverse of [`gather_to_master`](@ref): distribute a global subgrid field held on
`root` to every rank, returning each rank's local block (length
`length(gindex_level)`, in local task order). `global_field` is read on `root`
(length `numg`/`numc`/…); ignored on other ranks.

Returns `local[k] = global_field[gindex[k]]` for the local rank's points — i.e.
each rank picks out exactly its own subdomain by global index.

Single rank: returns `global_field[gindex]` which, with the identity gindex, is
the whole field — the master keeps its own block.

Mirrors CTSM's scatter-from-master used to distribute an initial decomposed
field (e.g. surfdata / restart read on the master then handed out per rank).
"""
function scatter_from_master(global_field::AbstractVector{T}, level::Int;
                             decomp_data::DecompData = decomp,
                             spmd::SPMDState = SPMD, root::Int = 0) where {T}
    gindex = _level_gindex(level; decomp_data = decomp_data)
    gsize  = _level_gsize(level; decomp_data = decomp_data)

    if !is_mpi_active()
        length(global_field) == gsize ||
            error("scatter_from_master: global_field length $(length(global_field)) != gsize $gsize at level $level")
        out = Vector{T}(undef, length(gindex))
        @inbounds for k in eachindex(gindex)
            out[k] = global_field[gindex[k]]
        end
        return out
    end

    # Multi-rank: the master indexes the global field by EACH rank's gindex to
    # form that rank's block, then scatters the concatenation. We mirror CLM's
    # serial-read-then-scatter: gather every rank's gindex to root, build the
    # per-rank blocks in rank order on root, and scatterv them back.
    npes = spmd.npes
    mylen = length(gindex)
    counts = _allgather_counts(mylen, npes, spmd)
    all_idxs = spmd_gatherv(collect(gindex), counts; root = root)

    sendbuf = if spmd.iam == root
        length(global_field) == gsize ||
            error("scatter_from_master: global_field length $(length(global_field)) != gsize $gsize at level $level")
        buf = Vector{T}(undef, sum(counts))
        @inbounds for k in eachindex(all_idxs)
            buf[k] = global_field[all_idxs[k]]
        end
        buf
    else
        T[]
    end
    return spmd_scatterv(sendbuf, counts; root = root)
end

# Allgather each rank's local length into a length-npes counts vector (Julia
# 1-based, indexed by rank+1). Mirrors the `mpi_allgather(mylength,...)` in
# spmdMod's procname gather. On a single rank this is just [mylen].
function _allgather_counts(mylen::Integer, npes::Integer, spmd::SPMDState)
    if !is_mpi_active() || npes == 1
        return [Int(mylen)]
    end
    comm = _spmd_comm()
    counts = MPI.Allgather(Cint(mylen), comm)
    return Int.(counts)
end
