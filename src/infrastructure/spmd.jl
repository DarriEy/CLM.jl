# spmd.jl — SPMD / MPI communication shim
#
# Julia port of the SPMD layer in CTSM/CLM5 (`src/utils/spmdMod.F90` for the
# init/identity bookkeeping, plus the `mpi_*` communication primitives CLM
# calls inline throughout the code: bcast, allreduce[SUM/MIN/MAX], gatherv,
# scatterv, barrier).
#
# CLM's land physics is embarrassingly parallel — there are no horizontal
# stencils across columns, so there is NO halo exchange. The only genuine
# cross-rank operations are:
#   * one global allreduce of index counts in decompInit,
#   * a parameter broadcast from the master rank,
#   * the gather-to-master of gridcell/column/patch vectors for serial I/O,
#   * the inverse scatter-from-master to distribute a decomposed field.
#
# Design: this layer works single-rank WITHOUT an MPI launch. When MPI is not
# initialized (or `npes == 1`), every primitive reduces to the trivial local
# operation:
#   bcast      = identity
#   allreduce  = identity (return the local value/array)
#   gatherv    = copy (master receives its own local block)
#   scatterv   = copy (master keeps its own local block)
#   barrier    = no-op
# This keeps the default single-process path and CI (which has no `mpiexec`)
# fully working, while the real MPI calls activate under a multi-rank launch.

using MPI

# ----------------------------------------------------------------------------
# Module-level SPMD state (mirrors the public vars of Fortran spmdMod)
# ----------------------------------------------------------------------------

"""
    SPMDState

Mutable bookkeeping for the SPMD layer — the Julia analogue of the public
`save`d variables in Fortran `spmdMod` (`masterproc`, `iam`, `npes`, `mpicom`,
`comp_id`), plus an `active` flag recording whether a real MPI environment is
in play.
"""
Base.@kwdef mutable struct SPMDState
    active::Bool = false     # is a real multi-rank MPI environment in play?
    masterproc::Bool = true  # proc-0 logical (true on single rank)
    iam::Int = 0             # processor (rank) number, 0-based like Fortran `iam`
    npes::Int = 1            # number of processors for CLM
    comp_id::Int = 1         # component id (LNDID)
    comm::Any = nothing      # the MPI communicator (MPI.Comm) when active
end

# The single global SPMD state. Initialized to the single-rank defaults so the
# layer is usable immediately without any spmd_init! call.
const SPMD = SPMDState()

# ----------------------------------------------------------------------------
# Init / finalize / query
# ----------------------------------------------------------------------------

"""
    spmd_init!(; comm = nothing, comp_id = 1) -> SPMDState

Initialize the SPMD layer. Port of Fortran `spmd_init`.

If MPI has been initialized in this process (i.e. the program was launched under
`mpiexec` and `MPI.Init()` was called) and the communicator has more than one
rank, this records the rank/size and activates the real MPI path. Otherwise it
leaves the layer in its single-rank identity mode (the default for CI and the
ordinary single-process driver).

`comm` defaults to `MPI.COMM_WORLD` when MPI is initialized.
"""
function spmd_init!(; comm = nothing, comp_id::Int = 1)
    SPMD.comp_id = comp_id
    if MPI.Initialized()
        c = comm === nothing ? MPI.COMM_WORLD : comm
        npes = MPI.Comm_size(c)
        rank = MPI.Comm_rank(c)
        SPMD.comm = c
        SPMD.npes = npes
        SPMD.iam = rank
        SPMD.masterproc = (rank == 0)
        # Only treat as a genuine distributed run when there is more than one
        # rank — a 1-rank launch behaves exactly like the serial path.
        SPMD.active = (npes > 1)
    else
        SPMD.comm = comm   # may still be nothing
        SPMD.npes = 1
        SPMD.iam = 0
        SPMD.masterproc = true
        SPMD.active = false
    end
    return SPMD
end

"""
    spmd_finalize!()

Tear-down hook. No-op in single-rank mode; does not call `MPI.Finalize()`
(lifetime of the MPI environment is owned by the launching program), it just
resets the layer to single-rank identity defaults.
"""
function spmd_finalize!()
    SPMD.active = false
    SPMD.masterproc = true
    SPMD.iam = 0
    SPMD.npes = 1
    SPMD.comm = nothing
    return nothing
end

"""
    is_mpi_active() -> Bool

True only when a genuine multi-rank MPI environment is in play. False for the
default single-process path (and for a degenerate 1-rank `mpiexec` launch).
"""
is_mpi_active() = SPMD.active

"""
    get_rank() -> Int

0-based processor (rank) number — Fortran `iam`. 0 on a single rank.
"""
get_rank() = SPMD.iam

"""
    get_npes() -> Int

Number of processors participating in the CLM computation — Fortran `npes`.
1 on a single rank.
"""
get_npes() = SPMD.npes

"""
    is_master() -> Bool

True on the master (rank-0) process — Fortran `masterproc`. Always true on a
single rank.
"""
is_master() = SPMD.masterproc

# Communicator accessor (falls back to COMM_WORLD when MPI is initialized).
function _spmd_comm()
    if SPMD.comm !== nothing
        return SPMD.comm
    end
    return MPI.Initialized() ? MPI.COMM_WORLD : nothing
end

# ----------------------------------------------------------------------------
# Broadcast — Fortran `mpi_bcast` from `root`
# ----------------------------------------------------------------------------

"""
    spmd_bcast!(buf; root = 0) -> buf

Broadcast `buf` from `root` to all ranks in place (mirrors CLM's parameter
broadcast from the master). On a single rank this is the identity — `buf` is
returned unchanged. `buf` may be an array (broadcast in place) or a scalar
(returned, since scalars can't be mutated in place).
"""
function spmd_bcast!(buf::AbstractArray; root::Int = 0)
    if is_mpi_active()
        MPI.Bcast!(buf, root, _spmd_comm())
    end
    return buf
end

# Scalar form: nothing to mutate in place, so return the (root's) value. On a
# single rank that is just the input. (Multi-rank scalar broadcast should use
# a 1-element array via the AbstractArray method.)
spmd_bcast!(buf; root::Int = 0) = buf

# ----------------------------------------------------------------------------
# Allreduce — Fortran `mpi_allreduce` with MPI_SUM / MPI_MIN / MPI_MAX
# ----------------------------------------------------------------------------

# Map a Julia reduction op symbol to the corresponding MPI op + local reducer.
const _ALLREDUCE_OPS = Dict{Symbol,Any}(
    :SUM => (MPI.SUM, +),
    :MIN => (MPI.MIN, min),
    :MAX => (MPI.MAX, max),
)

"""
    spmd_allreduce(val, op::Symbol = :SUM)

Global reduction across all ranks (Fortran `mpi_allreduce`). `op` is one of
`:SUM`, `:MIN`, `:MAX`.

* If `val` is a scalar, returns the reduced scalar.
* If `val` is an array, returns an array reduced element-wise across ranks
  (each element reduced independently — matches `mpi_allreduce` of a vector).

On a single rank this is the identity: the scalar (or array) is returned
unchanged, since the global reduction over one rank is the local value itself.
"""
function spmd_allreduce(val::Number, op::Symbol = :SUM)
    haskey(_ALLREDUCE_OPS, op) || throw(ArgumentError("spmd_allreduce: unknown op $op (use :SUM/:MIN/:MAX)"))
    if is_mpi_active()
        mpiop, _ = _ALLREDUCE_OPS[op]
        return MPI.Allreduce(val, mpiop, _spmd_comm())
    end
    return val
end

function spmd_allreduce(val::AbstractArray, op::Symbol = :SUM)
    haskey(_ALLREDUCE_OPS, op) || throw(ArgumentError("spmd_allreduce: unknown op $op (use :SUM/:MIN/:MAX)"))
    if is_mpi_active()
        mpiop, _ = _ALLREDUCE_OPS[op]
        return MPI.Allreduce(val, mpiop, _spmd_comm())
    end
    return copy(val)
end

"""
    spmd_allreduce_local(arr, op::Symbol = :SUM)

Reduce a *local* array down to a single scalar with `op`, then allreduce that
scalar across ranks. This is the common CLM idiom of "reduce my local block,
then combine globally" (e.g. a global sum of a per-gridcell field). On a single
rank it is exactly the local reduction of `arr`.
"""
function spmd_allreduce_local(arr::AbstractArray, op::Symbol = :SUM)
    haskey(_ALLREDUCE_OPS, op) || throw(ArgumentError("spmd_allreduce_local: unknown op $op (use :SUM/:MIN/:MAX)"))
    _, localred = _ALLREDUCE_OPS[op]
    local_val = reduce(localred, arr)
    return spmd_allreduce(local_val, op)
end

# ----------------------------------------------------------------------------
# Gatherv / Scatterv — distribute / collect a begg:endg vector across ranks
# ----------------------------------------------------------------------------
#
# `counts` is the per-proc element count (Fortran's `length`/`numg` per proc),
# a 0-based-indexed-by-rank vector here stored as a 1-based Julia Vector of
# length npes. The displacements are the running offsets, exactly as CLM builds
# `displ(i) = sum(length(0:i-1))` for the gridcell/column/patch gatherv/scatterv.

"""
    spmd_gatherv(sendbuf, counts; root = 0)

Gather variable-length local blocks onto `root` (Fortran `mpi_gatherv`).
`sendbuf` is this rank's local block (length `counts[get_rank()+1]`); `counts`
is the per-proc element count for every rank (Julia 1-based, length `npes`).

Returns, on `root`, the concatenated global vector (length `sum(counts)`); on
non-root ranks returns an empty vector of the same eltype. On a single rank
this is `copy(sendbuf)`.
"""
function spmd_gatherv(sendbuf::AbstractVector{T}, counts::AbstractVector{<:Integer}; root::Int = 0) where {T}
    if !is_mpi_active()
        return copy(sendbuf)
    end
    comm = _spmd_comm()
    cnts = Cint.(collect(counts))
    if get_rank() == root
        recvbuf = Vector{T}(undef, sum(cnts))
        MPI.Gatherv!(sendbuf, MPI.VBuffer(recvbuf, cnts), root, comm)
        return recvbuf
    else
        MPI.Gatherv!(sendbuf, nothing, root, comm)
        return T[]
    end
end

"""
    spmd_scatterv(sendbuf, counts; root = 0)

Scatter a global vector held on `root` into variable-length local blocks
(Fortran `mpi_scatterv`). On `root`, `sendbuf` is the full global vector
(length `sum(counts)`); on non-root ranks `sendbuf` is ignored. `counts` is the
per-proc element count for every rank (Julia 1-based, length `npes`).

Returns this rank's local block (length `counts[get_rank()+1]`). On a single
rank this is `copy(sendbuf)`.
"""
function spmd_scatterv(sendbuf::AbstractVector{T}, counts::AbstractVector{<:Integer}; root::Int = 0) where {T}
    if !is_mpi_active()
        return copy(sendbuf)
    end
    comm = _spmd_comm()
    cnts = Cint.(collect(counts))
    mycount = Int(cnts[get_rank() + 1])
    recvbuf = Vector{T}(undef, mycount)
    if get_rank() == root
        MPI.Scatterv!(MPI.VBuffer(sendbuf, cnts), recvbuf, root, comm)
    else
        MPI.Scatterv!(nothing, recvbuf, root, comm)
    end
    return recvbuf
end

# ----------------------------------------------------------------------------
# Barrier — Fortran `mpi_barrier`
# ----------------------------------------------------------------------------

"""
    spmd_barrier()

Synchronization barrier across all ranks (Fortran `mpi_barrier`). No-op on a
single rank.
"""
function spmd_barrier()
    if is_mpi_active()
        MPI.Barrier(_spmd_comm())
    end
    return nothing
end
