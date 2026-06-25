# ==========================================================================
# multigpu.jl — multi-GPU domain split over MPI (one GPU per rank)
#
# This wires the centralized backend layer (`backend.jl`) onto the MPI
# distributed driver (`distributed_driver.jl`) to give CLM's intended multi-GPU
# parallelism: ONE GPU PER MPI RANK. It mirrors CTSM's two-level decomposition
# (MPI over processes, threads/clumps within a process) but replaces the
# intra-process OpenMP with a GPU: each rank owns a subdomain (its gindex slice),
# binds to its LOCAL GPU, holds its rank-local state on that device, and runs its
# clumps' kernels there. Gather-to-master (already in distributed_driver.jl) moves
# results back to the host on rank 0 for serial NetCDF I/O.
#
#   global domain ──decompInit_distributed!──▶ rank r owns gindex_r (a subgrid slice)
#         │
#         ├─ rank 0  ── binds GPU local_rank 0 ── state_0 on device ── kernels ─┐
#         ├─ rank 1  ── binds GPU local_rank 1 ── state_1 on device ── kernels ─┤
#         │   ...                                                               │
#         └─ rank N  ── binds GPU local_rank N ── state_N on device ── kernels ─┘
#                                                                               │
#                          gather_to_master (host, global index order) ◀────────┘
#
# WHY ONE GPU PER RANK (vs many GPUs in one process)
# --------------------------------------------------
# CLM land physics is embarrassingly parallel over gridcells (no halo), so the
# cleanest multi-GPU model reuses the EXACT MPI decomposition that already gives
# bit-identical results (the CTSM PE-layout invariant): split the domain across
# ranks, bind each rank to one device. No new domain-split logic is needed — the
# round-robin gridcell→clump→rank map in `decompInit!` IS the GPU split, and
# `gather_to_master` already reassembles by global index regardless of layout.
# This also matches the standard HPC pattern (`CUDA_VISIBLE_DEVICES` per rank /
# `AMDGPU.device!(local_rank)`), and keeps each rank's working set within one
# GPU's memory.
#
# LOCAL-RANK → DEVICE BINDING
# ---------------------------
# Each node typically launches `ngpu_per_node` ranks. A rank's *node-local* rank
# (0..ngpu_per_node-1) selects which physical GPU it binds. We derive the local
# rank from the MPI shared-memory split when available (the robust way), falling
# back to `global_rank % ngpu_per_node`. The bind itself is delegated to the
# backend extension (CUDA `device!` / AMDGPU `device!`), so the base package
# names no GPU type.
#
# WHAT IS WIRED vs DESIGN-ONLY (HARDWARE HONESTY)
# -----------------------------------------------
# WIRED & runnable on this (Metal-only) machine:
#   * `clm_local_rank` / device-count plumbing and the rank→device index math
#     (validated on CPU: single-rank ⇒ local_rank 0).
#   * `clm_bind_local_gpu!` dispatch through the backend registry — on CPU it is
#     a no-op returning 0; the GPU bind closure is installed by the extension.
#   * `clm_run_multigpu!` — moves rank-local state to the active backend's device
#     (`clm_device_array`), runs `clm_run_distributed!`, and gathers requested
#     fields to the host master. On `:cpu` this is byte-identical to
#     `clm_run_distributed!` (device move is the identity).
# DESIGN-ONLY (cannot be validated here — needs real CUDA/AMD multi-GPU):
#   * the actual `CUDA.device!` / `AMDGPU.device!` binding effect, NCCL-free
#     correctness across >1 physical GPU, and inter-rank gather over a real
#     interconnect. The closures are in place (ext files) but UNEXERCISED on this
#     hardware. See report.
# ==========================================================================

"""
    clm_node_local_rank(; spmd::SPMDState = SPMD) -> Int

This rank's NODE-LOCAL rank (0-based): its index among the ranks sharing the same
physical node. Used to pick which GPU on the node this rank binds. Derived from an
MPI shared-memory communicator split when a real MPI env is active (the robust
way — handles uneven rank/node packing); on a single rank it is `0`.

Mirrors the common HPC idiom `local_rank = Comm_rank(split_shared(COMM_WORLD))`.
"""
function clm_node_local_rank(; spmd::SPMDState = SPMD)
    is_mpi_active() || return 0
    comm = _spmd_comm()
    comm === nothing && return spmd.iam
    # Split COMM_WORLD into per-node shared-memory communicators; this rank's rank
    # within its node comm is exactly the node-local rank.
    nodecomm = MPI.Comm_split_type(comm, MPI.COMM_TYPE_SHARED, spmd.iam)
    lr = MPI.Comm_rank(nodecomm)
    MPI.free(nodecomm)
    return Int(lr)
end

"""
    clm_gpu_device_count(; backend::Symbol = clm_backend()) -> Int

Number of GPUs visible to this process for `backend`. `0` on `:cpu` or when the
backend isn't registered. For a GPU backend this calls the extension-installed
device-count closure (CUDA `ndevices()` / AMDGPU device list length).
"""
function clm_gpu_device_count(; backend::Symbol = clm_backend())
    backend === :cpu && return 0
    haskey(_BACKEND_REGISTRY, backend) || return 0
    entry = _BACKEND_REGISTRY[backend]
    hasproperty(entry, :device_count) || return 0
    try
        return Int(entry.device_count())
    catch
        return 0
    end
end

"""
    clm_bind_local_gpu!(; backend::Symbol = clm_backend(), spmd::SPMDState = SPMD) -> Int

Bind THIS MPI rank to its local GPU for `backend` and return the bound device
index (0-based). The device index is `clm_node_local_rank() % device_count` (wrap
when ranks ≥ GPUs, the standard oversubscription fallback). On `:cpu`, or when no
GPU device is available, this is a no-op returning `0`.

The actual binding (`CUDA.device!(idx)` / `AMDGPU.device!(idx)`) is performed by
the backend extension's installed closure — the base package names no GPU type.
Call this ONCE per rank, right after `spmd_init!`/`decompInit_distributed!` and
before moving state to the device, so all subsequent `clm_device_array` /
kernel launches target the right physical GPU.
"""
function clm_bind_local_gpu!(; backend::Symbol = clm_backend(), spmd::SPMDState = SPMD)
    backend === :cpu && return 0
    haskey(_BACKEND_REGISTRY, backend) || return 0
    ndev = clm_gpu_device_count(backend = backend)
    ndev <= 0 && return 0
    idx = clm_node_local_rank(spmd = spmd) % ndev
    entry = _BACKEND_REGISTRY[backend]
    if hasproperty(entry, :bind_device)
        try
            entry.bind_device(idx)
        catch err
            @warn "clm_bind_local_gpu!: device bind failed for :$backend idx=$idx" exception = err
        end
    end
    return idx
end

"""
    clm_run_multigpu!(phys!, nsteps; gather_levels = (), gather_fields = (),
                      backend = clm_backend(), spmd = SPMD,
                      decomp_data = decomp, move_state! = nothing,
                      threaded = false) -> NamedTuple

One-GPU-per-rank distributed run. Each rank:

 1. binds its local GPU (`clm_bind_local_gpu!`);
 2. optionally moves its rank-local state onto the device — supply `move_state!()`
    to do the actual `clm_device_array` adaptation of your state tree (the base
    package can't know your container), or leave it `nothing` if `phys!` already
    closes over device-resident state;
 3. runs `nsteps` of `phys!` over ONLY its clumps via `clm_run_distributed!`
    (kernels execute on the bound device because the state arrays live there);
 4. gathers each `(host_field, level)` pair in `gather_fields` to the host master
    in global index order via `gather_to_master` (after moving the field to host
    with `move_field_to_host` if it is a device array).

Returns `(; bound_device, gathered)` where `gathered` is a Vector aligned with
`gather_fields`, each entry the master's global-order host vector (empty on
non-master ranks).

On `:cpu` (the only path this machine can run) every device step is the identity,
so the result is byte-identical to `clm_run_distributed!` followed by
`gather_to_master` — which is how the CPU proxy validates the wiring.
"""
function clm_run_multigpu!(phys!, nsteps::Integer;
                           gather_fields = (),
                           backend::Symbol = clm_backend(),
                           spmd::SPMDState = SPMD,
                           decomp_data::DecompData = decomp,
                           move_state! = nothing,
                           threaded::Bool = false)
    bound = clm_bind_local_gpu!(backend = backend, spmd = spmd)
    move_state! === nothing || move_state!()
    clm_run_distributed!(phys!, nsteps; decomp_data = decomp_data, spmd = spmd,
                         threaded = threaded)
    gathered = Vector{Any}(undef, length(gather_fields))
    for (k, (field, level)) in enumerate(gather_fields)
        host_field = _to_host_vector(field)
        gathered[k] = gather_to_master(host_field, level; decomp_data = decomp_data,
                                       spmd = spmd)
    end
    return (; bound_device = bound, gathered = gathered)
end

# Move a (possibly device) array to a plain host Vector for the host-side gather.
# `Array(x)` materializes a CuArray/ROCArray/MtlArray back to the CPU; a plain
# Vector passes through. Kept here (not in the extension) because `Array(::GPUArray)`
# is defined by the GPU package itself — no GPU type is named.
_to_host_vector(x::AbstractVector) = x isa Array ? x : convert(Vector, Array(x))
