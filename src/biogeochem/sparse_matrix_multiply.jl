# =============================================================================
# SparseMatrixMultiplyMod — sparse matrix infrastructure for the matrix-CN solver
#
# Julia port of CTSM `src/utils/SparseMatrixMultiplyMod.F90` (Author: Xingjie Lu).
#
# This is the foundation of the matrix-based C/N solution method
# (`use_soil_matrixcn` / `use_matrixcn`). The Veg/Soil matrix solvers
# (CNSoilMatrixMod / CNVegMatrixMod) advance the carbon/nitrogen pools as
#
#     X(t+1) = X(t) + A · X(t) + B · input
#
# where `A` is the sparse transfer/turnover matrix and `B` the sparse input
# matrix. To do this efficiently for every active land-unit (column / patch /
# gridcell) at once, the matrix VALUES are stored per unit (`M[u, k]`) while the
# SPARSITY STRUCTURE (row index RI, column index CI) is SHARED across all units.
# This batched "Structure-of-Arrays" layout matches the Fortran exactly and is
# GPU-friendly.
#
# Representation chosen: a faithful COO (coordinate) port rather than
# `SparseArrays`. `SparseArrays.SparseMatrixCSC` cannot model the
# shared-indices / per-unit-values batching that the Veg/Soil solvers rely on,
# and the spinup-accelerated SPMP/SPMM operators below depend on the exact
# entry ordering (row index varying faster than column index) that the COO
# layout preserves.
#
# COO ordering contract (from the Fortran header comment):
#   "Both row index and column index should be in ascending order. Row index
#    should change faster than Column index to ensure SPMP_AB work properly."
#   i.e. entries are sorted by the linear index  (CI-1)*SM + RI  (column-major).
# =============================================================================

# Sentinel values matching the Fortran `empty_int` / `empty_real`.
const SMM_EMPTY_INT  = -9999
const SMM_EMPTY_REAL = -9999.0

# -----------------------------------------------------------------------------
# Types
# -----------------------------------------------------------------------------

"""
    SparseMatrixType

General sparse matrix in COO (coordinate) format, batched over land units.

Fields (mirroring the Fortran `sparse_matrix_type`):
- `M::Matrix{Float64}`  — non-zero entries, indexed `M[u, k]` (unit, sparse-entry).
- `RI::Vector{Int}`     — row index of entry `k` (shared across units).
- `CI::Vector{Int}`     — column index of entry `k` (shared across units).
- `NE::Int`             — number of non-zero entries actually set.
- `SM::Int`             — size of the (SM × SM) matrix.
- `num_unit::Int`       — number of active units (informational).
- `begu::Int`, `endu::Int` — unit index bounds in the current process.

Units are indexed `begu:endu`; `M` is allocated `(endu-begu+1) × maxsm` and
addressed through the helper `u_idx` so the Fortran `M(u,k)` with `u ∈ begu:endu`
maps to a 1-based row.
"""
Base.@kwdef mutable struct SparseMatrixType
    M::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    RI::Vector{Int}    = Int[]
    CI::Vector{Int}    = Int[]
    NE::Int            = SMM_EMPTY_INT
    SM::Int            = SMM_EMPTY_INT
    num_unit::Int      = 0
    begu::Int          = SMM_EMPTY_INT
    endu::Int          = SMM_EMPTY_INT
end

"""
    DiagMatrixType

Diagonal matrix, storing only the diagonal entries `DM[u, i]`, batched over units.
Mirrors the Fortran `diag_matrix_type`.
"""
Base.@kwdef mutable struct DiagMatrixType
    DM::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    SM::Int             = SMM_EMPTY_INT
    num_unit::Int       = 0
    begu::Int           = SMM_EMPTY_INT
    endu::Int           = SMM_EMPTY_INT
end

"""
    VectorType

Dense vector `V[u, i]`, batched over units. Mirrors the Fortran `vector_type`.
"""
Base.@kwdef mutable struct VectorType
    V::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    SV::Int            = SMM_EMPTY_INT
    num_unit::Int      = 0
    begu::Int          = SMM_EMPTY_INT
    endu::Int          = SMM_EMPTY_INT
end

# Map a Fortran unit index `u ∈ begu:endu` to the 1-based row of the storage matrix.
@inline u_idx(begu::Int, u::Int) = u - begu + 1

# -----------------------------------------------------------------------------
# SparseMatrixType — allocation / status
# -----------------------------------------------------------------------------

"""
    init_sm!(this::SparseMatrixType, SM_in, begu, endu; maxsm=nothing)

Initialize the sparse matrix: set size `SM_in`, unit bounds `begu:endu`, and
allocate the value/index storage. `maxsm` (optional) caps the number of stored
entries (default `SM_in*SM_in`). Matches Fortran `InitSM`.
"""
function init_sm!(this::SparseMatrixType, SM_in::Int, begu::Int, endu::Int; maxsm::Union{Int,Nothing}=nothing)
    if is_alloc_sm(this)
        error("InitSM ERROR: Sparse Matrix was already allocated")
    end
    this.SM = SM_in
    this.begu = begu
    this.endu = endu
    nunit = endu - begu + 1
    if maxsm !== nothing
        @assert maxsm >= 1
        @assert maxsm <= SM_in * SM_in
        this.M = fill(SMM_EMPTY_REAL, nunit, maxsm)
    else
        this.M = fill(SMM_EMPTY_REAL, nunit, SM_in * SM_in)
    end
    this.RI = fill(SMM_EMPTY_INT, SM_in * SM_in)
    this.CI = fill(SMM_EMPTY_INT, SM_in * SM_in)
    this.NE = SMM_EMPTY_INT
    return nothing
end

"""
    release_sm!(this::SparseMatrixType)

Deallocate the sparse matrix data (reset to empty). Matches Fortran `ReleaseSM`.
"""
function release_sm!(this::SparseMatrixType)
    this.SM = SMM_EMPTY_INT
    this.begu = SMM_EMPTY_INT
    this.endu = SMM_EMPTY_INT
    this.M = Matrix{Float64}(undef, 0, 0)
    this.RI = Int[]
    this.CI = Int[]
    return nothing
end

"""
    is_alloc_sm(this::SparseMatrixType) -> Bool

True if `init_sm!` has been called (storage allocated). Matches `IsAllocSM`.
"""
is_alloc_sm(this::SparseMatrixType) = !isempty(this.M) || !isempty(this.RI) || !isempty(this.CI)

"""
    is_values_set_sm(this::SparseMatrixType) -> Bool

True if matrix values have been set (one of the `set_value_*` routines ran).
Matches `IsValuesSetSM`.
"""
function is_values_set_sm(this::SparseMatrixType)
    if !is_alloc_sm(this)
        return false
    elseif this.NE == SMM_EMPTY_INT
        return false
    else
        return true
    end
end

"""
    is_equiv_idx_sm(this::SparseMatrixType, A::SparseMatrixType) -> Bool

True if the two sparse matrices share the same sparsity structure (same SM, NE,
and identical RI/CI ordering). Matches `IsEquivIdxSM`.
"""
function is_equiv_idx_sm(this::SparseMatrixType, A::SparseMatrixType)
    if this.SM != A.SM
        return false
    end
    if this.NE == A.NE
        if all(this.RI[1:this.NE] .== A.RI[1:this.NE]) &&
           all(this.CI[1:this.NE] .== A.CI[1:this.NE])
            return true
        else
            return false
        end
    else
        return false
    end
end

# -----------------------------------------------------------------------------
# SparseMatrixType — set values
# -----------------------------------------------------------------------------

"""
    set_value_sm!(this, begu, endu, num_unit, filter_u, M, I, J, NE_in)

Set sparse-matrix values from non-zero values `M[u,k]` and their row/column
indices `I[k]`, `J[k]` (`k = 1:NE_in`). `filter_u[1:num_unit]` lists the active
unit indices. Matches Fortran `SetValueSM`.
"""
function set_value_sm!(this::SparseMatrixType, begu::Int, endu::Int, num_unit::Int,
                       filter_u::AbstractVector{Int}, M::AbstractMatrix{Float64},
                       I::AbstractVector{Int}, J::AbstractVector{Int}, NE_in::Int)
    if !is_alloc_sm(this)
        error("SetValueSM ERROR: Sparse Matrix was NOT already allocated")
    end
    @assert length(filter_u) >= num_unit
    @assert size(M, 2) >= NE_in
    @assert length(I) >= NE_in
    @assert length(J) >= NE_in
    @assert begu == this.begu
    @assert endu == this.endu
    for k in 1:NE_in
        for fu in 1:num_unit
            u = filter_u[fu]
            this.M[u_idx(this.begu, u), k] = M[u_idx(begu, u), k]
        end
    end
    this.NE = NE_in
    for k in 1:NE_in
        this.RI[k] = I[k]
        this.CI[k] = J[k]
    end
    return nothing
end

"""
    set_value_a_diag!(this, num_unit, filter_u, scaler)

Set `this` to a diagonal sparse matrix with constant value `scaler` on the
diagonal (RI[i]=CI[i]=i). Matches Fortran `SetValueA_diag`.
"""
function set_value_a_diag!(this::SparseMatrixType, num_unit::Int,
                           filter_u::AbstractVector{Int}, scaler::Float64)
    if !is_alloc_sm(this)
        error("SetValueA_diag ERROR: Sparse Matrix was NOT already allocated")
    end
    @assert length(filter_u) >= num_unit
    for i in 1:this.SM
        for fu in 1:num_unit
            u = filter_u[fu]
            this.M[u_idx(this.begu, u), i] = scaler
        end
    end
    for i in 1:this.SM
        this.RI[i] = i
        this.CI[i] = i
    end
    this.NE = this.SM
    return nothing
end

"""
    set_value_a!(this, begu, endu, num_unit, filter_u, M, AI, AJ, NE_NON, init_ready;
                 list=nothing, RI_A=nothing, CI_A=nothing) -> init_ready

Build the matrix `A = (off-diagonal entries M at AI,AJ) + (-1 on the diagonal)`.

On the first call (`init_ready == false`) it constructs A by adding a diagonal
`-1` matrix to the off-diagonal matrix via `spmp_ab!`, and (if `list`/`RI_A`/`CI_A`
are supplied) memorizes the entry map / index structure, returning
`init_ready = true`. On subsequent calls (`init_ready == true`, with the memorized
`list`/`RI_A`/`CI_A`) it reuses the structure for speed. Matches Fortran `SetValueA`.

Returns the (possibly updated) `init_ready` flag.
"""
function set_value_a!(this::SparseMatrixType, begu::Int, endu::Int, num_unit::Int,
                      filter_u::AbstractVector{Int}, M::AbstractMatrix{Float64},
                      AI::AbstractVector{Int}, AJ::AbstractVector{Int}, NE_NON::Int,
                      init_ready::Bool; list::Union{AbstractVector{Int},Nothing}=nothing,
                      RI_A::Union{AbstractVector{Int},Nothing}=nothing,
                      CI_A::Union{AbstractVector{Int},Nothing}=nothing)
    if !is_alloc_sm(this)
        error("SetValueA ERROR: Sparse Matrix was NOT already allocated")
    end
    if init_ready && !(list !== nothing && RI_A !== nothing && CI_A !== nothing)
        error("SetValueA ERROR: required optional arguments were NOT sent in")
    end
    @assert length(filter_u) >= num_unit
    @assert begu == this.begu
    @assert endu == this.endu
    @assert size(M, 2) >= NE_NON

    if init_ready
        for i in 1:(this.SM + NE_NON)
            for fu in 1:num_unit
                u = filter_u[fu]
                this.M[u_idx(this.begu, u), i] = -1.0
            end
        end
        for i in 1:NE_NON
            for fu in 1:num_unit
                u = filter_u[fu]
                this.M[u_idx(this.begu, u), list[i]] = M[u_idx(begu, u), i]
            end
        end
        this.NE = this.SM + NE_NON
        this.RI[1:this.NE] .= RI_A[1:this.NE]
        this.CI[1:this.NE] .= CI_A[1:this.NE]
    else
        A_diag = SparseMatrixType()
        A_nondiag = SparseMatrixType()
        init_sm!(A_diag, this.SM, begu, endu)
        init_sm!(A_nondiag, this.SM, begu, endu)

        set_value_a_diag!(A_diag, num_unit, filter_u, -1.0)
        set_value_sm!(A_nondiag, begu, endu, num_unit, filter_u, M, AI, AJ, NE_NON)

        list_ready = false
        if list !== nothing
            spmp_ab!(this, num_unit, filter_u, A_nondiag, A_diag, list_ready; list_A=list)
        else
            spmp_ab!(this, num_unit, filter_u, A_nondiag, A_diag, list_ready)
        end
        if RI_A !== nothing
            RI_A[1:this.NE] .= this.RI[1:this.NE]
        end
        if CI_A !== nothing
            CI_A[1:this.NE] .= this.CI[1:this.NE]
        end
        init_ready = true
        release_sm!(A_diag)
        release_sm!(A_nondiag)
    end
    return init_ready
end

"""
    set_value_copy_sm!(this, num_unit, filter_u, matrix)

Copy values + indices from `matrix` into `this`. Matches Fortran `SetValueCopySM`.
"""
function set_value_copy_sm!(this::SparseMatrixType, num_unit::Int,
                            filter_u::AbstractVector{Int}, matrix::SparseMatrixType)
    if !is_alloc_sm(this)
        error("SetValueCopySM ERROR: Sparse Matrix was NOT already allocated")
    end
    if !is_values_set_sm(matrix)
        error("SetValueCopySM ERROR: Sparse Matrix data sent in was NOT already set")
    end
    @assert this.SM == matrix.SM
    @assert this.begu == matrix.begu
    @assert this.endu == matrix.endu
    set_value_sm!(this, matrix.begu, matrix.endu, num_unit, filter_u,
                  matrix.M, matrix.RI, matrix.CI, matrix.NE)
    return nothing
end

"""
    copy_idx_sm!(this, matrix)

Copy only the RI/CI indices from `matrix` into `this` (values untouched), after
checking the implied number of non-empty entries matches. Matches `CopyIdxSM`.
"""
function copy_idx_sm!(this::SparseMatrixType, matrix::SparseMatrixType)
    if !is_alloc_sm(this)
        error("CopyIdxSM ERROR: Sparse Matrix was NOT already allocated")
    end
    if !is_values_set_sm(matrix)
        error("CopyIdxSM ERROR: Sparse Matrix data sent in was NOT already set")
    end
    @assert this.SM == matrix.SM
    @assert this.begu == matrix.begu
    @assert this.endu == matrix.endu
    # Figure out the number of non-empty data values.
    this.NE = size(this.M, 2)
    for i in 1:size(this.M, 2)
        if all(this.M[:, i] .== SMM_EMPTY_INT)
            this.NE = i - 1
            break
        end
    end
    if this.NE != matrix.NE
        error("CopyIdxSM ERROR: Sparse Matrix empty data size is different from input one copying the indices from")
    end
    this.RI[1:this.NE] .= matrix.RI[1:matrix.NE]
    this.CI[1:this.NE] .= matrix.CI[1:matrix.NE]
    return nothing
end

# -----------------------------------------------------------------------------
# DiagMatrixType
# -----------------------------------------------------------------------------

"""
    init_dm!(this::DiagMatrixType, SM_in, begu, endu)

Allocate a diagonal matrix of size `SM_in` over units `begu:endu`. Matches `InitDM`.
"""
function init_dm!(this::DiagMatrixType, SM_in::Int, begu::Int, endu::Int)
    if is_alloc_dm(this)
        error("InitDM ERROR: Diagonal Matrix was already allocated")
    end
    this.SM = SM_in
    this.DM = fill(SMM_EMPTY_REAL, endu - begu + 1, SM_in)
    this.begu = begu
    this.endu = endu
    return nothing
end

"""
    release_dm!(this::DiagMatrixType)

Deallocate the diagonal matrix. Matches `ReleaseDM`.
"""
function release_dm!(this::DiagMatrixType)
    this.SM = SMM_EMPTY_INT
    this.begu = SMM_EMPTY_INT
    this.endu = SMM_EMPTY_INT
    this.DM = Matrix{Float64}(undef, 0, 0)
    return nothing
end

"""
    is_alloc_dm(this::DiagMatrixType) -> Bool

True if `init_dm!` has been called. Matches `IsAllocDM`.
"""
is_alloc_dm(this::DiagMatrixType) = !isempty(this.DM)

"""
    is_values_set_dm(this::DiagMatrixType) -> Bool

True if the diagonal matrix has been initialized/set. Matches `IsValuesSetDM`.
"""
function is_values_set_dm(this::DiagMatrixType)
    if !is_alloc_dm(this)
        return false
    elseif this.SM == SMM_EMPTY_INT
        return false
    else
        return true
    end
end

"""
    set_value_dm!(this, begu, endu, num_unit, filter_u, M)

Set the diagonal entries `DM[u,i] = M[u,i]` for active units. Matches `SetValueDM`.
"""
function set_value_dm!(this::DiagMatrixType, begu::Int, endu::Int, num_unit::Int,
                       filter_u::AbstractVector{Int}, M::AbstractMatrix{Float64})
    if !is_alloc_dm(this)
        error("SetValueDM ERROR: Diagonal matrix was NOT already allocated")
    end
    @assert length(filter_u) >= num_unit
    @assert begu == this.begu
    @assert endu == this.endu
    @assert size(M, 2) >= this.SM
    for i in 1:this.SM
        for fu in 1:num_unit
            u = filter_u[fu]
            this.DM[u_idx(this.begu, u), i] = M[u_idx(begu, u), i]
        end
    end
    return nothing
end

"""
    set_value_copy_dm!(this, num_unit, filter_u, matrix)

Copy the diagonal matrix `matrix` into `this`. Matches `SetValueCopyDM`.
"""
function set_value_copy_dm!(this::DiagMatrixType, num_unit::Int,
                            filter_u::AbstractVector{Int}, matrix::DiagMatrixType)
    if !is_alloc_dm(this)
        error("SetValueCopyDM ERROR: Diagonal Matrix was NOT already allocated")
    end
    if !is_values_set_dm(matrix)
        error("SetValueCopyDM ERROR: Diagonal Matrix data sent in was NOT already set")
    end
    @assert this.SM == matrix.SM
    @assert this.begu == matrix.begu
    @assert this.endu == matrix.endu
    set_value_dm!(this, matrix.begu, matrix.endu, num_unit, filter_u, matrix.DM)
    return nothing
end

# -----------------------------------------------------------------------------
# VectorType
# -----------------------------------------------------------------------------

"""
    init_v!(this::VectorType, SV_in, begu, endu)

Allocate a vector of size `SV_in` over units `begu:endu`. Matches `InitV`.
"""
function init_v!(this::VectorType, SV_in::Int, begu::Int, endu::Int)
    if is_alloc_v(this)
        error("InitV ERROR: Vector was already allocated")
    end
    this.SV = SV_in
    this.V = fill(SMM_EMPTY_REAL, endu - begu + 1, SV_in)
    this.begu = begu
    this.endu = endu
    return nothing
end

"""
    release_v!(this::VectorType)

Deallocate the vector. Matches `ReleaseV`.
"""
function release_v!(this::VectorType)
    this.V = Matrix{Float64}(undef, 0, 0)
    this.begu = SMM_EMPTY_INT
    this.endu = SMM_EMPTY_INT
    this.SV = SMM_EMPTY_INT
    return nothing
end

"""
    is_alloc_v(this::VectorType) -> Bool

True if `init_v!` has been called. Matches `IsAllocV`.
"""
is_alloc_v(this::VectorType) = !isempty(this.V)

"""
    set_value_v_scaler!(this, num_unit, filter_u, scaler)

Set all entries of the vector to the constant `scaler`. Matches `SetValueV_scaler`.
"""
function set_value_v_scaler!(this::VectorType, num_unit::Int,
                             filter_u::AbstractVector{Int}, scaler::Float64)
    if !is_alloc_v(this)
        error("SetValueV_scaler ERROR: Vector was NOT already allocated")
    end
    @assert length(filter_u) >= num_unit
    for i in 1:this.SV
        for fu in 1:num_unit
            u = filter_u[fu]
            this.V[u_idx(this.begu, u), i] = scaler
        end
    end
    return nothing
end

"""
    set_value_v!(this, begu, endu, num_unit, filter_u, M)

Set the vector entries `V[u,i] = M[u,i]` for active units. Matches `SetValueV`.
"""
function set_value_v!(this::VectorType, begu::Int, endu::Int, num_unit::Int,
                      filter_u::AbstractVector{Int}, M::AbstractMatrix{Float64})
    if !is_alloc_v(this)
        error("SetValueV ERROR: Vector was NOT already allocated")
    end
    @assert length(filter_u) >= num_unit
    @assert begu == this.begu
    @assert endu == this.endu
    @assert size(M, 2) >= this.SV
    for i in 1:this.SV
        for fu in 1:num_unit
            u = filter_u[fu]
            this.V[u_idx(this.begu, u), i] = M[u_idx(begu, u), i]
        end
    end
    return nothing
end

# -----------------------------------------------------------------------------
# Multiply / add operators (the matrix-CN spinup-accelerated kernels)
# -----------------------------------------------------------------------------

"""
    spmm_ak!(this::SparseMatrixType, num_unit, filter_u, K::DiagMatrixType)

Sparse matrix × diagonal matrix, in place: `A(this) = A(this) · K`. Scales each
non-zero entry by the diagonal element of its column: `M[u,i] *= K.DM[u, CI[i]]`.
This is the AKX operator used to apply turnover/transfer rates. Matches `SPMM_AK`.
"""
function spmm_ak!(this::SparseMatrixType, num_unit::Int,
                  filter_u::AbstractVector{Int}, K::DiagMatrixType)
    @assert length(filter_u) >= num_unit
    @assert this.SM == K.SM
    for i in 1:this.NE
        for fu in 1:num_unit
            u = filter_u[fu]
            ui = u_idx(this.begu, u)
            this.M[ui, i] = this.M[ui, i] * K.DM[u_idx(K.begu, u), this.CI[i]]
        end
    end
    return nothing
end

"""
    spmm_ax!(this::VectorType, num_unit, filter_u, A::SparseMatrixType)

Sparse matrix × vector, accumulating in place: `X(this) = X(this) + A · X(this)`,
using the ORIGINAL `X` on the right (a temporary snapshot is taken first). This is
the pool-advance operator `X(t+1) = (I + A)·X(t)`. Matches `SPMM_AX`.
"""
function spmm_ax!(this::VectorType, num_unit::Int,
                  filter_u::AbstractVector{Int}, A::SparseMatrixType)
    @assert length(filter_u) >= num_unit
    @assert this.SV <= A.SM
    @assert size(A.M, 2) >= A.NE
    # Snapshot of the input vector (right-hand operand).
    V = copy(this.V)
    for i in 1:A.NE
        for fu in 1:num_unit
            u = filter_u[fu]
            ui = u_idx(this.begu, u)
            this.V[ui, A.RI[i]] = this.V[ui, A.RI[i]] + A.M[u_idx(A.begu, u), i] * V[ui, A.CI[i]]
        end
    end
    return nothing
end

"""
    spmp_b_acc!(this::SparseMatrixType, num_unit, filter_u, A::SparseMatrixType)

Sparse matrix accumulation in place: `B(this) = B(this) + A`. Requires A and B to
share the exact same entry locations (same RI/CI). Matches `SPMP_B_ACC`.
"""
function spmp_b_acc!(this::SparseMatrixType, num_unit::Int,
                     filter_u::AbstractVector{Int}, A::SparseMatrixType)
    @assert length(filter_u) >= num_unit
    @assert this.SM == A.SM
    @assert this.NE == A.NE
    @assert all(this.RI[1:this.NE] .== A.RI[1:A.NE])
    @assert all(this.CI[1:this.NE] .== A.CI[1:A.NE])
    for i in 1:A.NE
        for fu in 1:num_unit
            u = filter_u[fu]
            ui = u_idx(this.begu, u)
            this.M[ui, i] = this.M[ui, i] + A.M[u_idx(A.begu, u), i]
        end
    end
    return nothing
end

"""
    spmp_ab!(this::SparseMatrixType, num_unit, filter_u, A, B, list_ready;
             list_A=nothing, list_B=nothing, NE_AB=nothing, RI_AB=nothing, CI_AB=nothing)
        -> (list_ready, NE_AB)

Sparse matrix addition `AB(this) = A + B` via a merge of the two COO entry lists
(ordered by linear index `(CI-1)*SM + RI`). On the first call (`list_ready=false`)
it builds the merged structure and, if the optional memo arrays are supplied,
records the entry map (`list_A`/`list_B`), result indices (`RI_AB`/`CI_AB`) and
count (returned as the second tuple element) and flips `list_ready` to true. On
later calls (`list_ready=true`) it reuses the memorized structure. Matches `SPMP_AB`.

Returns `(list_ready, NE_AB)`. `NE_AB` is `this.NE` after the operation (or the
passed-in value on the memorized path).
"""
function spmp_ab!(this::SparseMatrixType, num_unit::Int, filter_u::AbstractVector{Int},
                  A::SparseMatrixType, B::SparseMatrixType, list_ready::Bool;
                  list_A::Union{AbstractVector{Int},Nothing}=nothing,
                  list_B::Union{AbstractVector{Int},Nothing}=nothing,
                  NE_AB::Union{Int,Nothing}=nothing,
                  RI_AB::Union{AbstractVector{Int},Nothing}=nothing,
                  CI_AB::Union{AbstractVector{Int},Nothing}=nothing)
    if list_ready && !(list_A !== nothing && list_B !== nothing && NE_AB !== nothing &&
                       RI_AB !== nothing && CI_AB !== nothing)
        error("SPMP_AB ERROR: missing required optional arguments")
    end
    @assert length(filter_u) >= num_unit
    @assert A.NE > 0
    @assert B.NE > 0
    @assert this.SM > 0
    @assert this.SM == A.SM
    @assert this.SM == B.SM

    if !list_ready
        # Linear (column-major) indices of A and B entries, with sentinel at end.
        Aindex = Vector{Int}(undef, A.NE + 1)
        Bindex = Vector{Int}(undef, B.NE + 1)
        @inbounds for k in 1:A.NE
            Aindex[k] = (A.CI[k] - 1) * A.SM + A.RI[k]
        end
        @inbounds for k in 1:B.NE
            Bindex[k] = (B.CI[k] - 1) * B.SM + B.RI[k]
        end
        Aindex[A.NE + 1] = A.SM * A.SM + 1
        Bindex[B.NE + 1] = B.SM * B.SM + 1
        ABindex = Vector{Int}(undef, this.SM * this.SM)

        i_a = 1; i_b = 1; i_ab = 1
        while i_a <= A.NE || i_b <= B.NE
            if Aindex[i_a] < Bindex[i_b]
                for fu in 1:num_unit
                    u = filter_u[fu]
                    this.M[u_idx(this.begu, u), i_ab] = A.M[u_idx(A.begu, u), i_a]
                end
                ABindex[i_ab] = Aindex[i_a]
                list_A !== nothing && (list_A[i_a] = i_ab)
                i_a += 1; i_ab += 1
            elseif Aindex[i_a] > Bindex[i_b]
                for fu in 1:num_unit
                    u = filter_u[fu]
                    this.M[u_idx(this.begu, u), i_ab] = B.M[u_idx(B.begu, u), i_b]
                end
                ABindex[i_ab] = Bindex[i_b]
                list_B !== nothing && (list_B[i_b] = i_ab)
                i_b += 1; i_ab += 1
            else
                for fu in 1:num_unit
                    u = filter_u[fu]
                    this.M[u_idx(this.begu, u), i_ab] =
                        A.M[u_idx(A.begu, u), i_a] + B.M[u_idx(B.begu, u), i_b]
                end
                ABindex[i_ab] = Aindex[i_a]
                list_A !== nothing && (list_A[i_a] = i_ab)
                list_B !== nothing && (list_B[i_b] = i_ab)
                i_a += 1; i_b += 1; i_ab += 1
            end
        end

        this.NE = i_ab - 1
        @inbounds for k in 1:this.NE
            this.CI[k] = div(ABindex[k] - 1, this.SM) + 1
            this.RI[k] = ABindex[k] - this.SM * (this.CI[k] - 1)
        end
        out_NE = this.NE
        NE_AB !== nothing && (out_NE = this.NE)
        CI_AB !== nothing && (CI_AB[1:this.NE] .= this.CI[1:this.NE])
        RI_AB !== nothing && (RI_AB[1:this.NE] .= this.RI[1:this.NE])
        if list_A !== nothing && list_B !== nothing && NE_AB !== nothing &&
           RI_AB !== nothing && CI_AB !== nothing
            list_ready = true
        end
        return (list_ready, out_NE)
    else
        for i in 1:NE_AB
            for fu in 1:num_unit
                u = filter_u[fu]
                this.M[u_idx(this.begu, u), i] = 0.0
            end
        end
        for i_a in 1:A.NE
            for fu in 1:num_unit
                u = filter_u[fu]
                this.M[u_idx(this.begu, u), list_A[i_a]] = A.M[u_idx(A.begu, u), i_a]
            end
        end
        for i_b in 1:B.NE
            for fu in 1:num_unit
                u = filter_u[fu]
                ui = u_idx(this.begu, u)
                this.M[ui, list_B[i_b]] = this.M[ui, list_B[i_b]] + B.M[u_idx(B.begu, u), i_b]
            end
        end
        this.NE = NE_AB
        this.CI[1:this.NE] .= CI_AB[1:NE_AB]
        this.RI[1:this.NE] .= RI_AB[1:NE_AB]
        return (list_ready, NE_AB)
    end
end

"""
    spmp_abc!(this::SparseMatrixType, num_unit, filter_u, A, B, C, list_ready;
              list_A=nothing, list_B=nothing, list_C=nothing,
              NE_ABC=nothing, RI_ABC=nothing, CI_ABC=nothing,
              num_actunit_A=nothing, filter_actunit_A=nothing,
              num_actunit_B=nothing, filter_actunit_B=nothing,
              num_actunit_C=nothing, filter_actunit_C=nothing) -> (list_ready, NE_ABC)

Sparse matrix addition `ABC(this) = A + B + C` via a 3-way merge of the COO entry
lists. As with `spmp_ab!`, the first call builds and (optionally) memorizes the
entry map / structure; later calls reuse it. On the memorized path each summand
may use its own active-unit filter (`num_actunit_*`/`filter_actunit_*`). Matches
`SPMP_ABC`. Returns `(list_ready, NE_ABC)`.
"""
function spmp_abc!(this::SparseMatrixType, num_unit::Int, filter_u::AbstractVector{Int},
                   A::SparseMatrixType, B::SparseMatrixType, C::SparseMatrixType,
                   list_ready::Bool;
                   list_A::Union{AbstractVector{Int},Nothing}=nothing,
                   list_B::Union{AbstractVector{Int},Nothing}=nothing,
                   list_C::Union{AbstractVector{Int},Nothing}=nothing,
                   NE_ABC::Union{Int,Nothing}=nothing,
                   RI_ABC::Union{AbstractVector{Int},Nothing}=nothing,
                   CI_ABC::Union{AbstractVector{Int},Nothing}=nothing,
                   num_actunit_A::Union{Int,Nothing}=nothing,
                   filter_actunit_A::Union{AbstractVector{Int},Nothing}=nothing,
                   num_actunit_B::Union{Int,Nothing}=nothing,
                   filter_actunit_B::Union{AbstractVector{Int},Nothing}=nothing,
                   num_actunit_C::Union{Int,Nothing}=nothing,
                   filter_actunit_C::Union{AbstractVector{Int},Nothing}=nothing)
    @assert this.SM > 0
    @assert A.NE > 0
    @assert B.NE > 0
    @assert C.NE > 0
    @assert this.SM == A.SM
    @assert this.SM == B.SM
    @assert this.SM == C.SM
    if list_ready && !(list_A !== nothing && list_B !== nothing && list_C !== nothing &&
                       NE_ABC !== nothing && RI_ABC !== nothing && CI_ABC !== nothing)
        error("SPMP_ABC ERROR: missing required optional arguments")
    end
    if num_actunit_A !== nothing
        num_actunit_A < 0 && error("SPMP_ABC ERROR: bad value for num_actunit_A")
        filter_actunit_A === nothing && error("SPMP_ABC ERROR: missing filter_actunit_A")
    end
    if num_actunit_B !== nothing
        num_actunit_B < 0 && error("SPMP_ABC ERROR: bad value for num_actunit_B")
        filter_actunit_B === nothing && error("SPMP_ABC ERROR: missing filter_actunit_B")
    end
    if num_actunit_C !== nothing
        num_actunit_C < 0 && error("SPMP_ABC ERROR: bad value for num_actunit_C")
        filter_actunit_C === nothing && error("SPMP_ABC ERROR: missing filter_actunit_C")
    end

    if !list_ready
        Aindex = Vector{Int}(undef, A.NE + 1)
        Bindex = Vector{Int}(undef, B.NE + 1)
        Cindex = Vector{Int}(undef, C.NE + 1)
        @inbounds for k in 1:A.NE
            Aindex[k] = (A.CI[k] - 1) * A.SM + A.RI[k]
        end
        @inbounds for k in 1:B.NE
            Bindex[k] = (B.CI[k] - 1) * B.SM + B.RI[k]
        end
        @inbounds for k in 1:C.NE
            Cindex[k] = (C.CI[k] - 1) * C.SM + C.RI[k]
        end
        Aindex[A.NE + 1] = A.SM * A.SM + 1
        Bindex[B.NE + 1] = B.SM * B.SM + 1
        Cindex[C.NE + 1] = C.SM * C.SM + 1
        ABCindex = Vector{Int}(undef, this.SM * this.SM)

        i_a = 1; i_b = 1; i_c = 1; i_abc = 1
        while i_a <= A.NE || i_b <= B.NE || i_c <= C.NE
            if Aindex[i_a] < Bindex[i_b] && Aindex[i_a] < Cindex[i_c]
                for fu in 1:num_unit
                    u = filter_u[fu]
                    this.M[u_idx(this.begu, u), i_abc] = A.M[u_idx(A.begu, u), i_a]
                end
                ABCindex[i_abc] = Aindex[i_a]
                list_A !== nothing && (list_A[i_a] = i_abc)
                i_a += 1; i_abc += 1
            elseif Bindex[i_b] < Aindex[i_a] && Bindex[i_b] < Cindex[i_c]
                for fu in 1:num_unit
                    u = filter_u[fu]
                    this.M[u_idx(this.begu, u), i_abc] = B.M[u_idx(B.begu, u), i_b]
                end
                ABCindex[i_abc] = Bindex[i_b]
                list_B !== nothing && (list_B[i_b] = i_abc)
                i_b += 1; i_abc += 1
            elseif Cindex[i_c] < Aindex[i_a] && Cindex[i_c] < Bindex[i_b]
                for fu in 1:num_unit
                    u = filter_u[fu]
                    this.M[u_idx(this.begu, u), i_abc] = C.M[u_idx(C.begu, u), i_c]
                end
                ABCindex[i_abc] = Cindex[i_c]
                list_C !== nothing && (list_C[i_c] = i_abc)
                i_c += 1; i_abc += 1
            elseif Aindex[i_a] == Bindex[i_b] && Aindex[i_a] < Cindex[i_c]
                for fu in 1:num_unit
                    u = filter_u[fu]
                    this.M[u_idx(this.begu, u), i_abc] =
                        A.M[u_idx(A.begu, u), i_a] + B.M[u_idx(B.begu, u), i_b]
                end
                ABCindex[i_abc] = Aindex[i_a]
                list_A !== nothing && (list_A[i_a] = i_abc)
                list_B !== nothing && (list_B[i_b] = i_abc)
                i_a += 1; i_b += 1; i_abc += 1
            elseif Aindex[i_a] == Cindex[i_c] && Aindex[i_a] < Bindex[i_b]
                for fu in 1:num_unit
                    u = filter_u[fu]
                    this.M[u_idx(this.begu, u), i_abc] =
                        A.M[u_idx(A.begu, u), i_a] + C.M[u_idx(C.begu, u), i_c]
                end
                ABCindex[i_abc] = Aindex[i_a]
                list_A !== nothing && (list_A[i_a] = i_abc)
                list_C !== nothing && (list_C[i_c] = i_abc)
                i_a += 1; i_c += 1; i_abc += 1
            elseif Bindex[i_b] == Cindex[i_c] && Bindex[i_b] < Aindex[i_a]
                for fu in 1:num_unit
                    u = filter_u[fu]
                    this.M[u_idx(this.begu, u), i_abc] =
                        B.M[u_idx(B.begu, u), i_b] + C.M[u_idx(C.begu, u), i_c]
                end
                ABCindex[i_abc] = Bindex[i_b]
                list_B !== nothing && (list_B[i_b] = i_abc)
                list_C !== nothing && (list_C[i_c] = i_abc)
                i_b += 1; i_c += 1; i_abc += 1
            elseif Aindex[i_a] == Bindex[i_b] && Aindex[i_a] == Cindex[i_c]
                for fu in 1:num_unit
                    u = filter_u[fu]
                    this.M[u_idx(this.begu, u), i_abc] =
                        A.M[u_idx(A.begu, u), i_a] + B.M[u_idx(B.begu, u), i_b] + C.M[u_idx(C.begu, u), i_c]
                end
                ABCindex[i_abc] = Bindex[i_b]
                list_A !== nothing && (list_A[i_a] = i_abc)
                list_B !== nothing && (list_B[i_b] = i_abc)
                list_C !== nothing && (list_C[i_c] = i_abc)
                i_a += 1; i_b += 1; i_c += 1; i_abc += 1
            else
                error("Error in subroutine SPMP_ABC: $(Aindex[i_a]) $(Bindex[i_b]) $(Cindex[i_c])")
            end
        end

        this.NE = i_abc - 1
        @inbounds for k in 1:this.NE
            this.CI[k] = div(ABCindex[k] - 1, this.SM) + 1
            this.RI[k] = ABCindex[k] - this.SM * (this.CI[k] - 1)
        end
        out_NE = this.NE
        CI_ABC !== nothing && (CI_ABC[1:this.NE] .= this.CI[1:this.NE])
        RI_ABC !== nothing && (RI_ABC[1:this.NE] .= this.RI[1:this.NE])
        if list_A !== nothing && list_B !== nothing && list_C !== nothing &&
           NE_ABC !== nothing && RI_ABC !== nothing && CI_ABC !== nothing
            list_ready = true
        end
        return (list_ready, out_NE)
    else
        for i in 1:NE_ABC
            for fu in 1:num_unit
                u = filter_u[fu]
                this.M[u_idx(this.begu, u), i] = 0.0
            end
        end
        if num_actunit_A !== nothing
            for i_a in 1:A.NE
                for fu in 1:num_actunit_A
                    u = filter_actunit_A[fu]
                    this.M[u_idx(this.begu, u), list_A[i_a]] = A.M[u_idx(A.begu, u), i_a]
                end
            end
        else
            for i_a in 1:A.NE
                for fu in 1:num_unit
                    u = filter_u[fu]
                    this.M[u_idx(this.begu, u), list_A[i_a]] = A.M[u_idx(A.begu, u), i_a]
                end
            end
        end
        if num_actunit_B !== nothing
            for i_b in 1:B.NE
                for fu in 1:num_actunit_B
                    u = filter_actunit_B[fu]
                    ui = u_idx(this.begu, u)
                    this.M[ui, list_B[i_b]] = this.M[ui, list_B[i_b]] + B.M[u_idx(B.begu, u), i_b]
                end
            end
        else
            for i_b in 1:B.NE
                for fu in 1:num_unit
                    u = filter_u[fu]
                    ui = u_idx(this.begu, u)
                    this.M[ui, list_B[i_b]] = this.M[ui, list_B[i_b]] + B.M[u_idx(B.begu, u), i_b]
                end
            end
        end
        if num_actunit_C !== nothing
            for i_c in 1:C.NE
                for fu in 1:num_actunit_C
                    u = filter_actunit_C[fu]
                    ui = u_idx(this.begu, u)
                    this.M[ui, list_C[i_c]] = this.M[ui, list_C[i_c]] + C.M[u_idx(C.begu, u), i_c]
                end
            end
        else
            for i_c in 1:C.NE
                for fu in 1:num_unit
                    u = filter_u[fu]
                    ui = u_idx(this.begu, u)
                    this.M[ui, list_C[i_c]] = this.M[ui, list_C[i_c]] + C.M[u_idx(C.begu, u), i_c]
                end
            end
        end
        this.NE = NE_ABC
        this.CI[1:this.NE] .= CI_ABC[1:NE_ABC]
        this.RI[1:this.NE] .= RI_ABC[1:NE_ABC]
        return (list_ready, NE_ABC)
    end
end

# -----------------------------------------------------------------------------
# Convenience: dense conversion (for testing / debugging against a dense reference)
# -----------------------------------------------------------------------------

"""
    sm_to_dense(this::SparseMatrixType, u::Int) -> Matrix{Float64}

Materialize the (SM × SM) dense matrix for unit `u` from the COO entries. Used by
tests to compare against a dense reference; not on the production hot path.
"""
function sm_to_dense(this::SparseMatrixType, u::Int)
    D = zeros(Float64, this.SM, this.SM)
    ui = u_idx(this.begu, u)
    for k in 1:this.NE
        D[this.RI[k], this.CI[k]] = this.M[ui, k]
    end
    return D
end
