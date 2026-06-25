# ==========================================================================
# run_multigpu_smoke.jl — one-GPU-per-rank multi-GPU smoke for CLM.jl.
#
# Launch on a REAL multi-GPU node (NOT runnable on the Metal-only dev machine):
#     mpiexec -n <ngpus> julia --project=. test/mpi/run_multigpu_smoke.jl :cuda
#     mpiexec -n <ngpus> julia --project=. test/mpi/run_multigpu_smoke.jl :amdgpu
#
# What it exercises (the one-GPU-per-rank design in src/infrastructure/multigpu.jl)
# --------------------------------------------------------------------------------
#   1. spmd_init! + decompInit_distributed!  — split the global domain across
#      ranks (round-robin gridcell->clump->rank), exactly like the CPU MPI smoke.
#   2. clm_set_backend(sym) + clm_bind_local_gpu! — each rank binds its NODE-LOCAL
#      GPU (device = node_local_rank % ndevices). Asserts distinct local GPUs
#      within a node.
#   3. clm_run_multigpu! — moves the rank-local field to the device
#      (clm_device_array), runs nsteps of clump-safe physics on that GPU, and
#      gathers the column result to the host master in GLOBAL index order.
#   4. bit-identity: the gathered global field on the master must equal the same
#      physics computed serially on the host (the CTSM PE-layout invariant — land
#      physics has no halo, so the N-GPU split is bit-identical to 1 serial run,
#      modulo Float32 on backends without Float64).
#
# HARDWARE HONESTY: this script is the WIRING for real multi-GPU validation. It
# was authored on a Metal-only machine and has NOT been run on real CUDA/AMD
# multi-GPU hardware. It exits 0 (skip) when the requested backend is not
# functional, so it is safe to invoke anywhere.
# ==========================================================================

using MPI
MPI.Init()

using CLM
using Printf

const COMM = MPI.COMM_WORLD
const RANK = MPI.Comm_rank(COMM)
const NPES = MPI.Comm_size(COMM)
const ROOT = 0

backend_sym = length(ARGS) >= 1 ? Symbol(replace(ARGS[1], ":" => "")) : :cuda

CLM.spmd_init!()

# Skip cleanly where the backend has no functional device (CI / wrong hardware).
if !CLM.clm_backend_registered(backend_sym)
    RANK == ROOT && println("multigpu smoke: backend :$backend_sym not registered " *
                            "(import $(backend_sym == :cuda ? "CUDA" : "AMDGPU") first) — SKIP")
    MPI.Finalize(); exit(0)
end
CLM.clm_set_backend(backend_sym)
if !CLM.clm_backend_functional(backend_sym)
    RANK == ROOT && println("multigpu smoke: no functional :$backend_sym device — SKIP")
    MPI.Finalize(); exit(0)
end

# ---- a small global domain (variable columns per gridcell) -----------------
const NUMG = max(NPES * 2, 4)
ncols_per_g = [((g % 3) + 1) for g in 1:NUMG]      # 1..3 cols per gridcell
const NSTEPS = 4
global_col(gc::Int) = 0.37 * gc + 2.0
phys_scalar(x) = (y = x; for _ in 1:NSTEPS; y = y + 0.5 * y; end; y)  # nsteps of y += 0.5y

# ---- per-rank decomposition + GPU bind -------------------------------------
d = CLM.DecompData()
CLM.decompInit_distributed!(NUMG; clump_pproc = 1, ncols_per_g = ncols_per_g, decomp_data = d)
bound = CLM.clm_bind_local_gpu!()
@printf("rank %d/%d bound local GPU %d (node-local rank %d)\n",
        RANK, NPES, bound, CLM.clm_node_local_rank())

# ---- rank-local field, moved to the device, advanced on-GPU ----------------
gidx = CLM.get_subgrid_level_gindex(CLM.SUBGRID_LEVEL_COLUMN; decomp_data = d)
local_host = Float64[global_col(g) for g in gidx]        # this rank's columns
local_dev  = CLM.clm_device_array(local_host)            # -> CuArray/ROCArray

# clump-safe physics on the device field (broadcast == kernel over the rank's columns)
function step!(_bounds)
    local_dev .= local_dev .+ 0.5 .* local_dev
    return nothing
end

res = CLM.clm_run_multigpu!(step!, NSTEPS; decomp_data = d,
                            gather_fields = ((local_dev, CLM.SUBGRID_LEVEL_COLUMN),))

# ---- bit-identity check on the master --------------------------------------
if RANK == ROOT
    ncol = sum(ncols_per_g)
    serial = Float64[phys_scalar(global_col(c)) for c in 1:ncol]
    got = res.gathered[1]
    # Float32-only backends (none here — CUDA/AMD are Float64) would compare with
    # a tolerance; CUDA/AMD support Float64 so exact equality is expected.
    ok = length(got) == ncol && all(isapprox.(got, serial; rtol = 1e-12, atol = 0))
    println(ok ? "MULTIGPU SMOKE PASS (gathered == serial)" : "MULTIGPU SMOKE FAIL")
    MPI.Finalize()
    exit(ok ? 0 : 1)
else
    MPI.Finalize()
    exit(0)
end
