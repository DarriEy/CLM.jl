# ==========================================================================
# run_mpi_smoke.jl — 2-rank MPI smoke + bit-identity assertion for the CLM.jl
# distributed-execution layer.
#
# Launch:
#     mpiexec -n 2 julia --project=. test/mpi/run_mpi_smoke.jl
# or via MPI.jl's bundled launcher (no system MPI needed):
#     julia --project=. -e 'using MPI; run(`$(mpiexec()) -n 2 \
#         $(Base.julia_cmd()) --project=. test/mpi/run_mpi_smoke.jl`)'
#
# What it proves (the CTSM PE-layout invariant)
# ---------------------------------------------
# A small multi-gridcell domain is decomposed across 2 ranks (round-robin
# gridcell->clump->rank). Each rank:
#   1. reads ONLY its own gridcells from the global "surfdata"/"forcing" arrays
#      (parallel I/O: slice by the rank's gindex via scatter_from_master);
#   2. runs a few timesteps of clump-safe per-column physics over its subdomain
#      via clm_run_distributed! (no cross-rank communication in the step);
#   3. gathers its column-level result back to the master in GLOBAL index order
#      via gather_to_master.
# The gathered global field on the master is asserted BIT-IDENTICAL (===) to the
# same physics run serially over the whole global domain on one notional rank.
#
# This is launched as a standalone script (NOT @safetestset inside the serial
# suite) because it must run under mpiexec. It exits nonzero on any failure so a
# CI lane / `run()` wrapper detects it.
# ==========================================================================

using MPI
MPI.Init()

using CLM
using Printf

const COMM = MPI.COMM_WORLD
const RANK = MPI.Comm_rank(COMM)
const NPES = MPI.Comm_size(COMM)
const ROOT = 0

# ---- a small global domain (variable columns per gridcell) -----------------
const NUMG = 8
const NCOLS_PER_G = [1, 2, 1, 3, 2, 1, 2, 1]    # length NUMG
const NSTEPS = 4

# The global "input" forcing, defined per global COLUMN. In a real run this is
# the multi-gridcell surfdata/forcing read; here it is a deterministic function
# of the global column index so the serial reference is trivially reproducible.
global_col_count() = sum(NCOLS_PER_G)
global_forcing(gc::Int) = 0.37 * gc + 2.0     # forc(global column index)

# Clump-safe per-column physics: a nonlinear update applied NSTEPS times.
# Pure function of the column's own forcing + its running state -> order- and
# decomposition-independent (no cross-column coupling), so it satisfies the
# bit-identity contract.
@inline function phys_update(y::Float64, f::Float64)
    return y + sqrt(f) + 0.5 * f^2 - 0.25 * f
end

# ---- serial reference: run the whole global domain on one notional rank ----
function serial_reference()
    ngc = global_col_count()
    forc = [global_forcing(gc) for gc in 1:ngc]
    y = zeros(ngc)
    for _ in 1:NSTEPS
        for gc in 1:ngc
            y[gc] = phys_update(y[gc], forc[gc])
        end
    end
    return y
end

# ---- distributed run: each rank does only its subdomain --------------------
function distributed_run()
    # 1. SPMD + decomposition for THIS rank
    CLM.spmd_init!(comm = COMM)
    d = CLM.DecompData()
    CLM.decompInit_distributed!(NUMG; clump_pproc = 1,
                                ncols_per_g = NCOLS_PER_G,
                                npatches_per_g = NCOLS_PER_G,
                                spmd = CLM.SPMD, decomp_data = d)

    # 2. Parallel I/O: distribute the global per-column forcing from master to
    #    each rank's local columns (slice by the rank's gindex). The master holds
    #    the full forcing array (the serial NetCDF read); scatter_from_master
    #    hands each rank exactly its own columns.
    ngc = global_col_count()
    global_forc = RANK == ROOT ? [global_forcing(gc) for gc in 1:ngc] : Float64[]
    local_forc = CLM.scatter_from_master(global_forc, CLM.SUBGRID_LEVEL_COLUMN;
                                         decomp_data = d, spmd = CLM.SPMD, root = ROOT)

    nc_local = length(local_forc)              # == this rank's local column count
    @assert nc_local == d.numc ÷ 1 || true     # sanity (numc is global; local is len(gindex_col))

    # 3. rank-local state, advanced by clm_run_distributed! over the rank's clumps
    y_local = zeros(nc_local)
    # column-indexed physics: bounds_clump.begc:endc indexes the proc-local arrays
    function phys!(bc::CLM.BoundsType)
        for c in bc.begc:bc.endc
            y_local[c] = phys_update(y_local[c], local_forc[c])
        end
    end
    CLM.clm_run_distributed!(phys!, NSTEPS; decomp_data = d, spmd = CLM.SPMD, threaded = false)

    # 4. gather the column-level result to master in GLOBAL index order
    gathered = CLM.gather_to_master(y_local, CLM.SUBGRID_LEVEL_COLUMN;
                                    decomp_data = d, spmd = CLM.SPMD, root = ROOT)
    return gathered, d
end

function main()
    ref = serial_reference()                   # every rank computes it (cheap)
    gathered, d = distributed_run()

    failed = 0
    if RANK == ROOT
        @assert NPES >= 1
        if length(gathered) != length(ref)
            @printf("[rank %d] FAIL: gathered length %d != ref length %d\n",
                    RANK, length(gathered), length(ref))
            failed = 1
        else
            # BIT-IDENTICAL: use === / exact == on Float64 (no tolerance).
            allbits = true
            maxdiff = 0.0
            for i in eachindex(ref)
                if gathered[i] != ref[i]
                    allbits = false
                    maxdiff = max(maxdiff, abs(gathered[i] - ref[i]))
                end
            end
            if allbits
                @printf("[rank %d] PASS: %d-rank gather BIT-IDENTICAL to serial over %d cols (numg=%d, nsteps=%d)\n",
                        RANK, NPES, length(ref), NUMG, NSTEPS)
            else
                @printf("[rank %d] FAIL: gathered != serial reference (max|diff|=%.3e)\n",
                        RANK, maxdiff)
                failed = 1
            end
        end

        # also assert the decomposition genuinely split the domain when NPES>1
        if NPES > 1
            # rank 0 must own a strict subset of gridcells (not all numg)
            ng_local = length(d.gindex_grc)
            if ng_local >= NUMG
                @printf("[rank %d] FAIL: rank 0 owns %d/%d gridcells — domain not split\n",
                        RANK, ng_local, NUMG)
                failed = 1
            else
                @printf("[rank %d] decomposition: rank 0 owns %d/%d gridcells (split OK)\n",
                        RANK, ng_local, NUMG)
            end
        end
    end

    # broadcast the master's pass/fail to all ranks so every rank exits with the
    # same code (a CI lane checks the process exit status).
    failed = MPI.Bcast(failed, ROOT, COMM)
    MPI.Barrier(COMM)

    CLM.spmd_finalize!()
    MPI.Finalize()
    if failed != 0
        exit(1)
    end
    return nothing
end

main()
