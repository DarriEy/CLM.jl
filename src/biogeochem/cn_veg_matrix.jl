# =============================================================================
# CNVegMatrixMod — matrix solution for the vegetation C (and N) cycle
#
# ┌────────────────────────────────────────────────────────────────────────┐
# │ GPU STATUS: OUT-OF-SCOPE for kernelization (2026-07). The matrix-CN      │
# │ solver (this file + cn_soil_matrix.jl + sparse_matrix_multiply.jl) is    │
# │ off by default (use_matrixcn/use_soil_matrixcn=false) AND NOT driver-    │
# │ wired: the CN cycle no-op's the matrix path (see c_state_update1.jl      │
# │ "matrix_update_phc path omitted / replaced by a no-op") and runs the ODE │
# │ path. It is unit-tested (test_cn_veg_matrix/soil_matrix/sparse_matrix)   │
# │ but never executes in a run.                                             │
# │ It is also a SPARSE-MATRIX solver (irregular sparse indexing/matmul) —   │
# │ the most GPU-hostile pattern — so kernelizing code that never runs is    │
# │ maximal effort for zero run-path value.                                  │
# │ TODO (if matrix-CN is ever driver-wired + made default-relevant): only   │
# │ then kernelize, and the sparse ops will need a CSR/atomic-scatter design │
# │ distinct from the elementwise/FUN-style playbook.                        │
# └────────────────────────────────────────────────────────────────────────┘
#
# Julia port of (the core of) CTSM `src/biogeochem/CNVegMatrixMod.F90`
# (matrix model by Yiqi Luo EcoLab: Drs. Xingjie Lu, Yuanyuan Huang, Zhengguang Du).
#
# When `use_matrixcn` is on, the per-flux sequential vegetation pool updates
# (`c_state_update1/2/2h/2g/3`) are NOT applied to the veg pools directly (those
# writes are gated off behind `!use_matrixcn` in the c_state_update*! kernels).
# Instead, every allocation / phenology-transfer / turnover / gap-mortality /
# fire flux is ASSEMBLED into a sparse transfer matrix and the pools advance in a
# single matrix solve:
#
#     X(n+1) = X(n) + (Aph·Kph + Agm·Kgm + Afi·Kfi)·X(n)·dt-fraction  +  B·I·dt
#            = (I + AKall)·X(n) + Binput
#
# where, per Fortran lines 1024–1538:
#   * `K*` (diagonal) holds the per-pool turnover FRACTION over the step
#     (`matrix_*turnover[p, pool]`, dimensionless 0..1).
#   * `A*` (off-diagonal at (receiver, doner)) holds the transfer FRACTION
#     `matrix_*transfer[p,k]·dt / matrix_*turnover[p, doner]`, and `-1` on the
#     diagonal.
#   * `AK = A·K` → entry (receiver,doner) = `matrix_*transfer[p,k]·dt`, and
#     diagonal (doner,doner) = `-matrix_*turnover[p,doner]`.
#   * `B·I` = `matrix_alloc[p,i]·matrix_Cinput[p]·dt` (the NPP allocation input).
#
# This is an EXACT reformulation of the sequential update: a pool gains the
# allocation input + all transfers into it, and loses its turnover (which equals
# the sum of the transfers out of it plus the litterfall/respiration losses).
#
# This port implements the CORE C-pool path (and the analogous N path) on the
# `sparse_matrix_multiply.jl` foundation:
#   - `cn_veg_matrix_assemble_ak!`  — build AKph/AKgm/AKfi and sum to AKall
#   - `cn_veg_matrix_solve_c!`      — advance the C pools and write them back
#   - `matrix_update_phc!/gmc!/fic!`— the turnover-rate accumulator helpers used
#     by phenology / gap-mortality / fire to populate the transfer/turnover arrays
#
# DEFERRED (clear TODOs, not needed for the matrix==sequential equivalence of the
# main C fluxes): spinup-acceleration (SASU) AKX / capacity (`matrix_cap_*`) at
# end-of-year, the C13/C14 isotope matrices, the explicit N-matrix solve wiring
# into n_state_update*, the `matrix_c*_acc` annual diagnostic accumulators, and
# the matrix restart (`CNVegMatrixRest`). The C-pool solve below is complete and
# byte-equivalent to the sequential path for the allocation + phenology-transfer +
# turnover + gap-mortality + fire fluxes.
# =============================================================================

# -----------------------------------------------------------------------------
# Transfer-count constants (Fortran clm_varpar.F90:328–353)
#
# These are the number of non-zero transfer entries per process. The matrix
# off-diagonal is built from the FIRST `n*trans - n*outtrans` entries (the
# in-vegetation transfers); the trailing `*outtrans` entries are pure losses to
# `ioutc`/`ioutn` (captured by the diagonal turnover, not a matrix entry).
# -----------------------------------------------------------------------------

"""
    veg_matrix_transfer_counts(use_crop) -> NamedTuple

Number of C/N transfer entries per process (phenology/gap-mortality/fire) and the
number of those that are pure out-of-vegetation losses. Mirrors the assignments in
`clm_varpar_init` (Fortran clm_varpar.F90:328–353).
"""
function veg_matrix_transfer_counts(use_crop::Bool)
    if use_crop
        ncphtrans = 18; nnphtrans = 37; ncphouttrans = 4; nnphouttrans = 5
    else
        ncphtrans = 17; nnphtrans = 34; ncphouttrans = 3; nnphouttrans = 4
    end
    return (
        ncphtrans = ncphtrans, ncphouttrans = ncphouttrans,
        nnphtrans = nnphtrans, nnphouttrans = nnphouttrans,
        ncgmtrans = 18, ncgmouttrans = 18,
        ncfitrans = 20, ncfiouttrans = 18,
        nngmtrans = 19, nngmouttrans = 19,
        nnfitrans = 21, nnfiouttrans = 19,
    )
end

# -----------------------------------------------------------------------------
# Turnover-rate accumulator helpers (Fortran functions matrix_update_ph/gm/fic)
#
# Phenology, gap-mortality and fire call these to register a transfer of rate
# `rate` (1/s) for transfer index `itransfer` (doner→receiver). They accumulate
# the doner's turnover fraction (capped so total turnover ≤ 1 over the step) and
# the transfer's transfer rate. `acc=true` accumulates; `acc=false` replaces.
# `matrixcheck=true` applies the cap. These bang-functions mutate the matrix
# arrays on the carbon/nitrogen flux struct and return the (capped) applied rate.
# -----------------------------------------------------------------------------

"""
    matrix_update_phc!(cf, p, itransfer, rate, dt; matrixcheck=true, acc=true) -> applied_rate

Register a phenology C transfer (doner→receiver, index `itransfer`) at rate `rate`
(gC/gC/s). Updates `matrix_phturnover_patch[p, doner]` and `matrix_phtransfer_patch[p, itransfer]`.
Port of Fortran `matrix_update_phc` (CNVegMatrixMod.F90:3579).
"""
function matrix_update_phc!(cf, p::Int, itransfer::Int, rate::Float64, dt::Float64;
                            matrixcheck::Bool=true, acc::Bool=true)
    phtransfer = cf.matrix_phtransfer_patch
    phturnover = cf.matrix_phturnover_patch
    doner = cf.matrix_phtransfer_doner_patch
    d = doner[itransfer]
    if matrixcheck
        if acc && phturnover[p, d] + rate * dt >= 1.0
            applied = max(0.0, (1.0 - phturnover[p, d]) / dt)
        else
            applied = rate
        end
    else
        applied = rate
    end
    if acc
        phturnover[p, d] += applied * dt
        phtransfer[p, itransfer] += applied
    else
        phturnover[p, d] += -phtransfer[p, itransfer] * dt + applied * dt
        phtransfer[p, itransfer] = applied
    end
    return applied
end

"""
    matrix_update_gmc!(cf, p, itransfer, rate, dt; matrixcheck=true, acc=true) -> applied_rate

Register a gap-mortality C transfer. The cap accounts for the already-accumulated
phenology turnover (`matrix_phturnover`) plus gap-mortality turnover. Port of
Fortran `matrix_update_gmc` (CNVegMatrixMod.F90:3617).
"""
function matrix_update_gmc!(cf, p::Int, itransfer::Int, rate::Float64, dt::Float64;
                            matrixcheck::Bool=true, acc::Bool=true)
    phturnover = cf.matrix_phturnover_patch
    gmtransfer = cf.matrix_gmtransfer_patch
    gmturnover = cf.matrix_gmturnover_patch
    doner = cf.matrix_gmtransfer_doner_patch
    d = doner[itransfer]
    if matrixcheck
        if acc && phturnover[p, d] + gmturnover[p, d] + rate * dt >= 1.0
            applied = max(0.0, (1.0 - phturnover[p, d] - gmturnover[p, d]) / dt)
        else
            applied = rate
        end
    else
        applied = rate
    end
    if acc
        gmturnover[p, d] += applied * dt
        gmtransfer[p, itransfer] += applied
    else
        gmturnover[p, d] += -gmtransfer[p, itransfer] * dt + applied * dt
        gmtransfer[p, itransfer] = applied
    end
    return applied
end

"""
    matrix_update_fic!(cf, p, itransfer, rate, dt; matrixcheck=true, acc=true) -> applied_rate

Register a fire C transfer. The cap accounts for phenology + gap-mortality + fire
turnover already accumulated. Port of Fortran `matrix_update_fic`
(CNVegMatrixMod.F90:3657).
"""
function matrix_update_fic!(cf, p::Int, itransfer::Int, rate::Float64, dt::Float64;
                            matrixcheck::Bool=true, acc::Bool=true)
    phturnover = cf.matrix_phturnover_patch
    gmturnover = cf.matrix_gmturnover_patch
    fitransfer = cf.matrix_fitransfer_patch
    fiturnover = cf.matrix_fiturnover_patch
    doner = cf.matrix_fitransfer_doner_patch
    d = doner[itransfer]
    if matrixcheck
        if acc && phturnover[p, d] + gmturnover[p, d] + fiturnover[p, d] + rate * dt >= 1.0
            applied = max(0.0, (1.0 - phturnover[p, d] - gmturnover[p, d] - fiturnover[p, d]) / dt)
        else
            applied = rate
        end
    else
        applied = rate
    end
    if acc
        fiturnover[p, d] += applied * dt
        fitransfer[p, itransfer] += applied
    else
        fiturnover[p, d] += -fitransfer[p, itransfer] * dt + applied * dt
        fitransfer[p, itransfer] = applied
    end
    return applied
end

# -----------------------------------------------------------------------------
# Nitrogen accumulators — identical mechanics to the C versions but writing the
# matrix_nph*/matrix_ngm*/matrix_nfi* arrays. Port of matrix_update_phn/gmn/fin.
# -----------------------------------------------------------------------------

"""
    matrix_update_phn!(nf, p, itransfer, rate, dt; matrixcheck=true, acc=true) -> applied_rate

Register a phenology N transfer (doner→receiver, index `itransfer`). N analog of
[`matrix_update_phc!`]; writes `matrix_nphturnover_patch`/`matrix_nphtransfer_patch`.
"""
function matrix_update_phn!(nf, p::Int, itransfer::Int, rate::Float64, dt::Float64;
                            matrixcheck::Bool=true, acc::Bool=true)
    phtransfer = nf.matrix_nphtransfer_patch
    phturnover = nf.matrix_nphturnover_patch
    doner = nf.matrix_nphtransfer_doner_patch
    d = doner[itransfer]
    if matrixcheck
        if acc && phturnover[p, d] + rate * dt >= 1.0
            applied = max(0.0, (1.0 - phturnover[p, d]) / dt)
        else
            applied = rate
        end
    else
        applied = rate
    end
    if acc
        phturnover[p, d] += applied * dt
        phtransfer[p, itransfer] += applied
    else
        phturnover[p, d] += -phtransfer[p, itransfer] * dt + applied * dt
        phtransfer[p, itransfer] = applied
    end
    return applied
end

"""
    matrix_update_gmn!(nf, p, itransfer, rate, dt; matrixcheck=true, acc=true) -> applied_rate

Register a gap-mortality N transfer. N analog of [`matrix_update_gmc!`].
"""
function matrix_update_gmn!(nf, p::Int, itransfer::Int, rate::Float64, dt::Float64;
                            matrixcheck::Bool=true, acc::Bool=true)
    phturnover = nf.matrix_nphturnover_patch
    gmtransfer = nf.matrix_ngmtransfer_patch
    gmturnover = nf.matrix_ngmturnover_patch
    doner = nf.matrix_ngmtransfer_doner_patch
    d = doner[itransfer]
    if matrixcheck
        if acc && phturnover[p, d] + gmturnover[p, d] + rate * dt >= 1.0
            applied = max(0.0, (1.0 - phturnover[p, d] - gmturnover[p, d]) / dt)
        else
            applied = rate
        end
    else
        applied = rate
    end
    if acc
        gmturnover[p, d] += applied * dt
        gmtransfer[p, itransfer] += applied
    else
        gmturnover[p, d] += -gmtransfer[p, itransfer] * dt + applied * dt
        gmtransfer[p, itransfer] = applied
    end
    return applied
end

"""
    matrix_update_fin!(nf, p, itransfer, rate, dt; matrixcheck=true, acc=true) -> applied_rate

Register a fire N transfer. N analog of [`matrix_update_fic!`].
"""
function matrix_update_fin!(nf, p::Int, itransfer::Int, rate::Float64, dt::Float64;
                            matrixcheck::Bool=true, acc::Bool=true)
    phturnover = nf.matrix_nphturnover_patch
    gmturnover = nf.matrix_ngmturnover_patch
    fitransfer = nf.matrix_nfitransfer_patch
    fiturnover = nf.matrix_nfiturnover_patch
    doner = nf.matrix_nfitransfer_doner_patch
    d = doner[itransfer]
    if matrixcheck
        if acc && phturnover[p, d] + gmturnover[p, d] + fiturnover[p, d] + rate * dt >= 1.0
            applied = max(0.0, (1.0 - phturnover[p, d] - gmturnover[p, d] - fiturnover[p, d]) / dt)
        else
            applied = rate
        end
    else
        applied = rate
    end
    if acc
        fiturnover[p, d] += applied * dt
        fitransfer[p, itransfer] += applied
    else
        fiturnover[p, d] += -fitransfer[p, itransfer] * dt + applied * dt
        fitransfer[p, itransfer] = applied
    end
    return applied
end

# -----------------------------------------------------------------------------
# Matrix assembly + solve (Fortran CNVegMatrix core, lines 1024–1538 + 2309–2336)
# -----------------------------------------------------------------------------

"""
    _build_ak_process!(AK, begp, endp, num_soilp, filter_soilp, Aoned, transfer, turnover,
                       doner, receiver, ntrans, nouttrans, nvegpool, init_ready, list, RI, CI)
        -> init_ready

Assemble one process matrix `AK = A·K` for the patches in `filter_soilp`, where
`A` has off-diagonal fraction entries `Aoned[p,k] = transfer[p,k]·dt / turnover[p,doner[k]]`
at (receiver[k], doner[k]) and `-1` on the diagonal, and `K` is the diagonal
`turnover[p, pool]`. The result entry (receiver,doner) = `transfer[p,k]·dt` and the
diagonal (pool,pool) = `-turnover[p,pool]`.

`Aoned` is a pre-sized `(np × (ntrans-nouttrans))` scratch matrix that this fills.
Mirrors Fortran lines 1024–1073 (Aoned) + 1407–1483 (SetValueA / SetValueDM / SPMM_AK).
Returns the (possibly updated) `init_ready` flag for the memoized index structure.

CONTRACT: the first `ntrans-nouttrans` entries of `doner`/`receiver` (the in-veg
transfers) MUST be ordered so the COO linear index `(doner[k]-1)·nvegpool + receiver[k]`
is strictly ascending — i.e. column-major, row varying fastest. This is exactly how
Fortran `InitTransfer` (CNVegCarbonFluxType.F90:568) authors them, and the
`spmp_ab!` merge in `set_value_a!` relies on it (the diagonal `-1` matrix is the
sorted operand). The trailing `nouttrans` entries are pure losses to `ioutc`/`ioutn`
and are NOT placed in the matrix (their effect enters via the diagonal turnover).
"""
function _build_ak_process!(AK::SparseMatrixType, begp::Int, endp::Int, num_soilp::Int,
                            filter_soilp::AbstractVector{Int},
                            Aoned::AbstractMatrix{Float64},
                            transfer::AbstractMatrix{Float64}, turnover::AbstractMatrix{Float64},
                            doner::AbstractVector{Int}, receiver::AbstractVector{Int},
                            ntrans::Int, nouttrans::Int, nvegpool::Int,
                            dt::Float64, init_ready::Bool,
                            list::AbstractVector{Int}, RI::AbstractVector{Int}, CI::AbstractVector{Int})
    nnon = ntrans - nouttrans
    Kveg = DiagMatrixType()
    init_dm!(Kveg, nvegpool, begp, endp)

    if nnon > 0
        # Aoned[p,k] = transfer[p,k]*dt / turnover[p, doner[k]]   (Fortran 1028–1032)
        for k in 1:nnon
            d = doner[k]
            for fp in 1:num_soilp
                p = filter_soilp[fp]
                pr = u_idx(begp, p)
                if turnover[pr, d] != 0.0
                    Aoned[pr, k] = transfer[pr, k] * dt / turnover[pr, d]
                else
                    Aoned[pr, k] = 0.0
                end
            end
        end
        AI = receiver[1:nnon]              # row indices  (receivers)
        AJ = doner[1:nnon]                 # column indices (doners)
        init_ready = set_value_a!(AK, begp, endp, num_soilp, filter_soilp, Aoned,
                                  AI, AJ, nnon, init_ready; list=list, RI_A=RI, CI_A=CI)
    else
        set_value_a_diag!(AK, num_soilp, filter_soilp, -1.0)
    end

    # K diagonal = turnover[p, 1:nvegpool]   (Fortran SetValueDM, 1420)
    set_value_dm!(Kveg, begp, endp, num_soilp, filter_soilp, view(turnover, :, 1:nvegpool))
    # AK = A·K   (Fortran SPMM_AK, 1423)
    spmm_ak!(AK, num_soilp, filter_soilp, Kveg)
    release_dm!(Kveg)
    return init_ready
end

"""
    cn_veg_matrix_solve_c!(cs_veg, cf_veg; mask_soilp, bounds_patch, ivt, woody,
                           npcropmin, nvegcpool, counts, dt,
                           num_actfirep=0, irepr=1)

Advance the vegetation CARBON pools by one step via the matrix solution
`X(n+1) = (I + AKph + AKgm + AKfi)·X(n) + matrix_alloc·matrix_Cinput·dt`, and write
the result back to the leaf/froot/livestem/deadstem/livecroot/deadcroot {+ st/xf}
pools (and grain for crops). The transfer/turnover matrices must already have been
populated (via the `matrix_update_*!` helpers / phenology / gap-mortality / fire).

Port of the C-pool path of Fortran `CNVegMatrix` (lines 1116–1538, 2309–2336). The
sequential c_state_update*! veg pool writes are gated off (`!use_matrixcn`) when
this runs, so this is the sole pool advance.

`counts` is `veg_matrix_transfer_counts(use_crop)`. `num_actfirep` selects the
2-process (no fire) vs 3-process (with fire) sum, matching Fortran 1503–1510.
"""
function cn_veg_matrix_solve_c!(cs_veg::CNVegCarbonStateData, cf_veg::CNVegCarbonFluxData;
                                mask_soilp::AbstractVector{Bool},
                                bounds_patch::UnitRange{Int},
                                ivt::AbstractVector{<:Integer},
                                woody::AbstractVector,
                                npcropmin::Int, nvegcpool::Int,
                                counts, dt::Real,
                                num_actfirep::Int=0, irepr::Int=1)
    dt = Float64(dt)
    begp = first(bounds_patch); endp = last(bounds_patch)

    # Build the active-patch filter from the mask (Fortran filter_soilp).
    filter_soilp = Int[p for p in bounds_patch if mask_soilp[p]]
    num_soilp = length(filter_soilp)
    if num_soilp == 0
        return nothing
    end

    cf = cf_veg; cs = cs_veg

    # --- Xvegc input vector and old-state load (Fortran 992, 1119–1149) ---
    Xvegc = VectorType()
    init_v!(Xvegc, nvegcpool, begp, endp)
    Binput = VectorType()
    init_v!(Binput, nvegcpool, begp, endp)

    for fp in 1:num_soilp
        p = filter_soilp[fp]; pr = u_idx(begp, p)
        Xvegc.V[pr, ILEAF]         = cs.leafc_patch[p]
        Xvegc.V[pr, ILEAF_ST]      = cs.leafc_storage_patch[p]
        Xvegc.V[pr, ILEAF_XF]      = cs.leafc_xfer_patch[p]
        Xvegc.V[pr, IFROOT]        = cs.frootc_patch[p]
        Xvegc.V[pr, IFROOT_ST]     = cs.frootc_storage_patch[p]
        Xvegc.V[pr, IFROOT_XF]     = cs.frootc_xfer_patch[p]
        Xvegc.V[pr, ILIVESTEM]     = cs.livestemc_patch[p]
        Xvegc.V[pr, ILIVESTEM_ST]  = cs.livestemc_storage_patch[p]
        Xvegc.V[pr, ILIVESTEM_XF]  = cs.livestemc_xfer_patch[p]
        Xvegc.V[pr, IDEADSTEM]     = cs.deadstemc_patch[p]
        Xvegc.V[pr, IDEADSTEM_ST]  = cs.deadstemc_storage_patch[p]
        Xvegc.V[pr, IDEADSTEM_XF]  = cs.deadstemc_xfer_patch[p]
        Xvegc.V[pr, ILIVECROOT]    = cs.livecrootc_patch[p]
        Xvegc.V[pr, ILIVECROOT_ST] = cs.livecrootc_storage_patch[p]
        Xvegc.V[pr, ILIVECROOT_XF] = cs.livecrootc_xfer_patch[p]
        Xvegc.V[pr, IDEADCROOT]    = cs.deadcrootc_patch[p]
        Xvegc.V[pr, IDEADCROOT_ST] = cs.deadcrootc_storage_patch[p]
        Xvegc.V[pr, IDEADCROOT_XF] = cs.deadcrootc_xfer_patch[p]
        if ivt[p] >= npcropmin && nvegcpool >= IGRAIN_XF
            Xvegc.V[pr, IGRAIN]    = cs.reproductivec_patch[p, irepr]
            Xvegc.V[pr, IGRAIN_ST] = cs.reproductivec_storage_patch[p, irepr]
            Xvegc.V[pr, IGRAIN_XF] = cs.reproductivec_xfer_patch[p, irepr]
        end
    end

    # --- B·I : allocation input (Fortran 1400–1405) ---
    for i in 1:nvegcpool
        for fp in 1:num_soilp
            p = filter_soilp[fp]; pr = u_idx(begp, p)
            Binput.V[pr, i] = cf.matrix_alloc_patch[p, i] * cf.matrix_Cinput_patch[p] * dt
        end
    end

    # --- Assemble AKph, AKgm, AKfi (Fortran 1407–1497) ---
    AKph = SparseMatrixType(); init_sm!(AKph, nvegcpool, begp, endp)
    AKgm = SparseMatrixType(); init_sm!(AKgm, nvegcpool, begp, endp)
    AKfi = SparseMatrixType(); init_sm!(AKfi, nvegcpool, begp, endp)

    Aph = fill(0.0, endp - begp + 1, max(1, counts.ncphtrans - counts.ncphouttrans))
    Agm = fill(0.0, endp - begp + 1, max(1, counts.ncgmtrans - counts.ncgmouttrans))
    Afi = fill(0.0, endp - begp + 1, max(1, counts.ncfitrans - counts.ncfiouttrans))

    nn = nvegcpool * nvegcpool
    list_ph = fill(0, nn); RI_ph = fill(0, nn); CI_ph = fill(0, nn)
    list_gm = fill(0, nn); RI_gm = fill(0, nn); CI_gm = fill(0, nn)
    list_fi = fill(0, nn); RI_fi = fill(0, nn); CI_fi = fill(0, nn)

    _build_ak_process!(AKph, begp, endp, num_soilp, filter_soilp, Aph,
                       cf.matrix_phtransfer_patch, cf.matrix_phturnover_patch,
                       cf.matrix_phtransfer_doner_patch, cf.matrix_phtransfer_receiver_patch,
                       counts.ncphtrans, counts.ncphouttrans, nvegcpool, dt, false,
                       list_ph, RI_ph, CI_ph)
    _build_ak_process!(AKgm, begp, endp, num_soilp, filter_soilp, Agm,
                       cf.matrix_gmtransfer_patch, cf.matrix_gmturnover_patch,
                       cf.matrix_gmtransfer_doner_patch, cf.matrix_gmtransfer_receiver_patch,
                       counts.ncgmtrans, counts.ncgmouttrans, nvegcpool, dt, false,
                       list_gm, RI_gm, CI_gm)
    _build_ak_process!(AKfi, begp, endp, num_soilp, filter_soilp, Afi,
                       cf.matrix_fitransfer_patch, cf.matrix_fiturnover_patch,
                       cf.matrix_fitransfer_doner_patch, cf.matrix_fitransfer_receiver_patch,
                       counts.ncfitrans, counts.ncfiouttrans, nvegcpool, dt, false,
                       list_fi, RI_fi, CI_fi)

    # --- AKall = AKph + AKgm (+ AKfi if fire active) (Fortran 1503–1510) ---
    AKall = SparseMatrixType(); init_sm!(AKall, nvegcpool, begp, endp)
    if num_actfirep == 0
        la = fill(0, nn); lb = fill(0, nn); RIab = fill(0, nn); CIab = fill(0, nn)
        spmp_ab!(AKall, num_soilp, filter_soilp, AKph, AKgm, false;
                 list_A=la, list_B=lb, NE_AB=0, RI_AB=RIab, CI_AB=CIab)
    else
        la = fill(0, nn); lb = fill(0, nn); lc = fill(0, nn)
        RIabc = fill(0, nn); CIabc = fill(0, nn)
        spmp_abc!(AKall, num_soilp, filter_soilp, AKph, AKgm, AKfi, false;
                  list_A=la, list_B=lb, list_C=lc, NE_ABC=0, RI_ABC=RIabc, CI_ABC=CIabc)
    end

    # --- Xvegc = (I + AKall)·Xvegc  (Fortran SPMM_AX, 1530) ---
    spmm_ax!(Xvegc, num_soilp, filter_soilp, AKall)
    # --- Xvegc += B·I  (Fortran 1533–1538) ---
    for i in 1:nvegcpool
        for fp in 1:num_soilp
            p = filter_soilp[fp]; pr = u_idx(begp, p)
            Xvegc.V[pr, i] = Xvegc.V[pr, i] + Binput.V[pr, i]
        end
    end

    # --- Write the advanced pools back (Fortran 2309–2336) ---
    for fp in 1:num_soilp
        p = filter_soilp[fp]; pr = u_idx(begp, p)
        cs.leafc_patch[p]            = Xvegc.V[pr, ILEAF]
        cs.leafc_storage_patch[p]    = Xvegc.V[pr, ILEAF_ST]
        cs.leafc_xfer_patch[p]       = Xvegc.V[pr, ILEAF_XF]
        cs.frootc_patch[p]           = Xvegc.V[pr, IFROOT]
        cs.frootc_storage_patch[p]   = Xvegc.V[pr, IFROOT_ST]
        cs.frootc_xfer_patch[p]      = Xvegc.V[pr, IFROOT_XF]
        cs.livestemc_patch[p]        = Xvegc.V[pr, ILIVESTEM]
        cs.livestemc_storage_patch[p] = Xvegc.V[pr, ILIVESTEM_ST]
        cs.livestemc_xfer_patch[p]   = Xvegc.V[pr, ILIVESTEM_XF]
        cs.deadstemc_patch[p]        = Xvegc.V[pr, IDEADSTEM]
        cs.deadstemc_storage_patch[p] = Xvegc.V[pr, IDEADSTEM_ST]
        cs.deadstemc_xfer_patch[p]   = Xvegc.V[pr, IDEADSTEM_XF]
        cs.livecrootc_patch[p]       = Xvegc.V[pr, ILIVECROOT]
        cs.livecrootc_storage_patch[p] = Xvegc.V[pr, ILIVECROOT_ST]
        cs.livecrootc_xfer_patch[p]  = Xvegc.V[pr, ILIVECROOT_XF]
        cs.deadcrootc_patch[p]       = Xvegc.V[pr, IDEADCROOT]
        cs.deadcrootc_storage_patch[p] = Xvegc.V[pr, IDEADCROOT_ST]
        cs.deadcrootc_xfer_patch[p]  = Xvegc.V[pr, IDEADCROOT_XF]
        if ivt[p] >= npcropmin && nvegcpool >= IGRAIN_XF
            cs.reproductivec_patch[p, irepr]         = Xvegc.V[pr, IGRAIN]
            cs.reproductivec_storage_patch[p, irepr] = Xvegc.V[pr, IGRAIN_ST]
            cs.reproductivec_xfer_patch[p, irepr]    = Xvegc.V[pr, IGRAIN_XF]
        end
    end

    release_sm!(AKph); release_sm!(AKgm); release_sm!(AKfi); release_sm!(AKall)
    release_v!(Xvegc); release_v!(Binput)
    return nothing
end

# =============================================================================
# DEFERRED MATRIX-CN TAILS — N-pool solve, C13/C14 isotope solves, and the
# veg SASU spinup capacity. These extend the matrix model to the remaining pool
# species and the spinup accelerator. They are gated behind use_matrixcn (and
# the corresponding use_c13/use_c14/spinup flags) exactly like the C solve, and
# share the same `_build_ak_process!` assembly + `spmm_ax!` advance, so each is
# an EXACT reformulation of the corresponding sequential update for its species.
# =============================================================================

"""
    _cn_veg_matrix_advance!(X, B, begp, endp, num_soilp, filter_soilp, nvegpool,
                            phtransfer, phturnover, ph_doner, ph_receiver, nphtrans, nphouttrans,
                            gmtransfer, gmturnover, gm_doner, gm_receiver, ngmtrans, ngmouttrans,
                            fitransfer, fiturnover, fi_doner, fi_receiver, nfitrans, nfiouttrans,
                            dt, num_actfirep)

Generic batched matrix advance `X ← (I + AKph + AKgm [+ AKfi])·X + B` for the veg
pool vector `X` (already loaded) and input vector `B` (already = alloc·input·dt).
This is the species-agnostic core shared by the C, N and isotope solves — it
assembles each process matrix `AK = A·K` from its `transfer`/`turnover`/`doner`/
`receiver` arrays via [`_build_ak_process!`] and sums them, then applies the
in-place `spmm_ax!` advance and adds `B`. The pool LOAD and WRITE-BACK (which
differ per species) stay in the species-specific wrappers. Port of the matrix
assemble+solve of `CNVegMatrix` (Fortran 1407–1538), factored so the N and
isotope variants reuse it verbatim.
"""
function _cn_veg_matrix_advance!(X::VectorType, B::VectorType,
        begp::Int, endp::Int, num_soilp::Int, filter_soilp::AbstractVector{Int}, nvegpool::Int,
        phtransfer::AbstractMatrix{Float64}, phturnover::AbstractMatrix{Float64},
        ph_doner::AbstractVector{Int}, ph_receiver::AbstractVector{Int}, nphtrans::Int, nphouttrans::Int,
        gmtransfer::AbstractMatrix{Float64}, gmturnover::AbstractMatrix{Float64},
        gm_doner::AbstractVector{Int}, gm_receiver::AbstractVector{Int}, ngmtrans::Int, ngmouttrans::Int,
        fitransfer::AbstractMatrix{Float64}, fiturnover::AbstractMatrix{Float64},
        fi_doner::AbstractVector{Int}, fi_receiver::AbstractVector{Int}, nfitrans::Int, nfiouttrans::Int,
        dt::Float64, num_actfirep::Int)

    AKph = SparseMatrixType(); init_sm!(AKph, nvegpool, begp, endp)
    AKgm = SparseMatrixType(); init_sm!(AKgm, nvegpool, begp, endp)
    AKfi = SparseMatrixType(); init_sm!(AKfi, nvegpool, begp, endp)

    Aph = fill(0.0, endp - begp + 1, max(1, nphtrans - nphouttrans))
    Agm = fill(0.0, endp - begp + 1, max(1, ngmtrans - ngmouttrans))
    Afi = fill(0.0, endp - begp + 1, max(1, nfitrans - nfiouttrans))

    nn = nvegpool * nvegpool
    list_ph = fill(0, nn); RI_ph = fill(0, nn); CI_ph = fill(0, nn)
    list_gm = fill(0, nn); RI_gm = fill(0, nn); CI_gm = fill(0, nn)
    list_fi = fill(0, nn); RI_fi = fill(0, nn); CI_fi = fill(0, nn)

    _build_ak_process!(AKph, begp, endp, num_soilp, filter_soilp, Aph,
                       phtransfer, phturnover, ph_doner, ph_receiver,
                       nphtrans, nphouttrans, nvegpool, dt, false, list_ph, RI_ph, CI_ph)
    _build_ak_process!(AKgm, begp, endp, num_soilp, filter_soilp, Agm,
                       gmtransfer, gmturnover, gm_doner, gm_receiver,
                       ngmtrans, ngmouttrans, nvegpool, dt, false, list_gm, RI_gm, CI_gm)
    _build_ak_process!(AKfi, begp, endp, num_soilp, filter_soilp, Afi,
                       fitransfer, fiturnover, fi_doner, fi_receiver,
                       nfitrans, nfiouttrans, nvegpool, dt, false, list_fi, RI_fi, CI_fi)

    AKall = SparseMatrixType(); init_sm!(AKall, nvegpool, begp, endp)
    if num_actfirep == 0
        la = fill(0, nn); lb = fill(0, nn); RIab = fill(0, nn); CIab = fill(0, nn)
        spmp_ab!(AKall, num_soilp, filter_soilp, AKph, AKgm, false;
                 list_A=la, list_B=lb, NE_AB=0, RI_AB=RIab, CI_AB=CIab)
    else
        la = fill(0, nn); lb = fill(0, nn); lc = fill(0, nn)
        RIabc = fill(0, nn); CIabc = fill(0, nn)
        spmp_abc!(AKall, num_soilp, filter_soilp, AKph, AKgm, AKfi, false;
                  list_A=la, list_B=lb, list_C=lc, NE_ABC=0, RI_ABC=RIabc, CI_ABC=CIabc)
    end

    spmm_ax!(X, num_soilp, filter_soilp, AKall)
    for i in 1:nvegpool
        for fp in 1:num_soilp
            p = filter_soilp[fp]; pr = u_idx(begp, p)
            X.V[pr, i] = X.V[pr, i] + B.V[pr, i]
        end
    end

    release_sm!(AKph); release_sm!(AKgm); release_sm!(AKfi); release_sm!(AKall)
    return nothing
end

"""
    cn_veg_matrix_solve_n!(ns_veg, nf_veg; mask_soilp, bounds_patch, ivt,
                           npcropmin, nvegnpool, counts, dt, num_actfirep=0, irepr=1)

Advance the vegetation NITROGEN pools by one step via the matrix solution
`X(n+1) = (I + AKnph + AKngm + AKnfi)·X(n) + matrix_nalloc·matrix_Ninput·dt`, and
write the result back to the leaf/froot/livestem/deadstem/livecroot/deadcroot
{+ st/xf} N pools, the retranslocated-N pool (`iretransn = nvegnpool`, the last
pool, which has no storage/transfer compartment), and grain N for crops.

Port of the N-pool path of Fortran `CNVegMatrix` (Xvegn load 1221–1240, the
`n*transfer/n*turnover` assembly, and the Xvegn write-back). This is the EXACT
reformulation of the sequential `n_state_update*` veg pool writes, which are gated
off (`!use_matrixcn`) when this runs. `counts` is `veg_matrix_transfer_counts`.
"""
function cn_veg_matrix_solve_n!(ns_veg::CNVegNitrogenStateData, nf_veg::CNVegNitrogenFluxData;
                                mask_soilp::AbstractVector{Bool},
                                bounds_patch::UnitRange{Int},
                                ivt::AbstractVector{<:Integer},
                                npcropmin::Int, nvegnpool::Int,
                                counts, dt::Real,
                                num_actfirep::Int=0, irepr::Int=1)
    dt = Float64(dt)
    begp = first(bounds_patch); endp = last(bounds_patch)
    filter_soilp = Int[p for p in bounds_patch if mask_soilp[p]]
    num_soilp = length(filter_soilp)
    num_soilp == 0 && return nothing

    nf = nf_veg; ns = ns_veg
    iretransn = nvegnpool   # last N pool (Fortran iretransn = nvegnpool)

    Xvegn = VectorType(); init_v!(Xvegn, nvegnpool, begp, endp)
    Binput = VectorType(); init_v!(Binput, nvegnpool, begp, endp)

    # --- load N pools (Fortran 1221–1247) ---
    for fp in 1:num_soilp
        p = filter_soilp[fp]; pr = u_idx(begp, p)
        Xvegn.V[pr, ILEAF]         = ns.leafn_patch[p]
        Xvegn.V[pr, ILEAF_ST]      = ns.leafn_storage_patch[p]
        Xvegn.V[pr, ILEAF_XF]      = ns.leafn_xfer_patch[p]
        Xvegn.V[pr, IFROOT]        = ns.frootn_patch[p]
        Xvegn.V[pr, IFROOT_ST]     = ns.frootn_storage_patch[p]
        Xvegn.V[pr, IFROOT_XF]     = ns.frootn_xfer_patch[p]
        Xvegn.V[pr, ILIVESTEM]     = ns.livestemn_patch[p]
        Xvegn.V[pr, ILIVESTEM_ST]  = ns.livestemn_storage_patch[p]
        Xvegn.V[pr, ILIVESTEM_XF]  = ns.livestemn_xfer_patch[p]
        Xvegn.V[pr, IDEADSTEM]     = ns.deadstemn_patch[p]
        Xvegn.V[pr, IDEADSTEM_ST]  = ns.deadstemn_storage_patch[p]
        Xvegn.V[pr, IDEADSTEM_XF]  = ns.deadstemn_xfer_patch[p]
        Xvegn.V[pr, ILIVECROOT]    = ns.livecrootn_patch[p]
        Xvegn.V[pr, ILIVECROOT_ST] = ns.livecrootn_storage_patch[p]
        Xvegn.V[pr, ILIVECROOT_XF] = ns.livecrootn_xfer_patch[p]
        Xvegn.V[pr, IDEADCROOT]    = ns.deadcrootn_patch[p]
        Xvegn.V[pr, IDEADCROOT_ST] = ns.deadcrootn_storage_patch[p]
        Xvegn.V[pr, IDEADCROOT_XF] = ns.deadcrootn_xfer_patch[p]
        Xvegn.V[pr, iretransn]     = ns.retransn_patch[p]
        if ivt[p] >= npcropmin && nvegnpool >= IGRAIN_XF
            Xvegn.V[pr, IGRAIN]    = ns.reproductiven_patch[p, irepr]
            Xvegn.V[pr, IGRAIN_ST] = ns.reproductiven_storage_patch[p, irepr]
            Xvegn.V[pr, IGRAIN_XF] = ns.reproductiven_xfer_patch[p, irepr]
        end
    end

    # --- B·I : N allocation input (Fortran matrix_nalloc·matrix_Ninput·dt) ---
    for i in 1:nvegnpool
        for fp in 1:num_soilp
            p = filter_soilp[fp]; pr = u_idx(begp, p)
            Binput.V[pr, i] = nf.matrix_nalloc_patch[p, i] * nf.matrix_Ninput_patch[p] * dt
        end
    end

    _cn_veg_matrix_advance!(Xvegn, Binput, begp, endp, num_soilp, filter_soilp, nvegnpool,
        nf.matrix_nphtransfer_patch, nf.matrix_nphturnover_patch,
        nf.matrix_nphtransfer_doner_patch, nf.matrix_nphtransfer_receiver_patch,
        counts.nnphtrans, counts.nnphouttrans,
        nf.matrix_ngmtransfer_patch, nf.matrix_ngmturnover_patch,
        nf.matrix_ngmtransfer_doner_patch, nf.matrix_ngmtransfer_receiver_patch,
        counts.nngmtrans, counts.nngmouttrans,
        nf.matrix_nfitransfer_patch, nf.matrix_nfiturnover_patch,
        nf.matrix_nfitransfer_doner_patch, nf.matrix_nfitransfer_receiver_patch,
        counts.nnfitrans, counts.nnfiouttrans,
        dt, num_actfirep)

    # --- write the advanced N pools back ---
    for fp in 1:num_soilp
        p = filter_soilp[fp]; pr = u_idx(begp, p)
        ns.leafn_patch[p]             = Xvegn.V[pr, ILEAF]
        ns.leafn_storage_patch[p]     = Xvegn.V[pr, ILEAF_ST]
        ns.leafn_xfer_patch[p]        = Xvegn.V[pr, ILEAF_XF]
        ns.frootn_patch[p]            = Xvegn.V[pr, IFROOT]
        ns.frootn_storage_patch[p]    = Xvegn.V[pr, IFROOT_ST]
        ns.frootn_xfer_patch[p]       = Xvegn.V[pr, IFROOT_XF]
        ns.livestemn_patch[p]         = Xvegn.V[pr, ILIVESTEM]
        ns.livestemn_storage_patch[p] = Xvegn.V[pr, ILIVESTEM_ST]
        ns.livestemn_xfer_patch[p]    = Xvegn.V[pr, ILIVESTEM_XF]
        ns.deadstemn_patch[p]         = Xvegn.V[pr, IDEADSTEM]
        ns.deadstemn_storage_patch[p] = Xvegn.V[pr, IDEADSTEM_ST]
        ns.deadstemn_xfer_patch[p]    = Xvegn.V[pr, IDEADSTEM_XF]
        ns.livecrootn_patch[p]        = Xvegn.V[pr, ILIVECROOT]
        ns.livecrootn_storage_patch[p] = Xvegn.V[pr, ILIVECROOT_ST]
        ns.livecrootn_xfer_patch[p]   = Xvegn.V[pr, ILIVECROOT_XF]
        ns.deadcrootn_patch[p]        = Xvegn.V[pr, IDEADCROOT]
        ns.deadcrootn_storage_patch[p] = Xvegn.V[pr, IDEADCROOT_ST]
        ns.deadcrootn_xfer_patch[p]   = Xvegn.V[pr, IDEADCROOT_XF]
        ns.retransn_patch[p]          = Xvegn.V[pr, iretransn]
        if ivt[p] >= npcropmin && nvegnpool >= IGRAIN_XF
            ns.reproductiven_patch[p, irepr]         = Xvegn.V[pr, IGRAIN]
            ns.reproductiven_storage_patch[p, irepr] = Xvegn.V[pr, IGRAIN_ST]
            ns.reproductiven_xfer_patch[p, irepr]    = Xvegn.V[pr, IGRAIN_XF]
        end
    end

    release_v!(Xvegn); release_v!(Binput)
    return nothing
end

"""
    cn_veg_matrix_solve_iso!(Xiso, Biso; cf_veg, mask_soilp, bounds_patch,
                             nvegcpool, counts, dt, num_actfirep=0)

Advance an ISOTOPE (C13 or C14) veg-carbon pool vector by the matrix solution
`Xiso(n+1) = (I + AKph + AKgm + AKfi)·Xiso(n) + Biso`. The isotope pools ride the
SAME transfer/turnover operator `A` as the bulk carbon (the cascade topology and
fractional transfer rates are isotope-independent — Fortran reuses the same
`matrix_*transfer`/`matrix_*turnover` arrays from `cnveg_carbonflux_inst` for the
C13/C14 advance, lines 99–137, 2845–2857); only the input vector `Biso`
(`matrix_*alloc·matrix_Cinput13/14·dt`) carries the isotopic signal.

`Xiso` and `Biso` are pre-loaded `(np × nvegcpool)` matrices over the active
patches (1-based patch row). Returns the advanced `Xiso` (mutated in place). This
keeps the isotope pool storage decoupled from the bulk C state structs while
faithfully reusing the bulk-C `A` operator. Port of the C13/C14 branch of
`CNVegMatrix`.
"""
function cn_veg_matrix_solve_iso!(Xiso::AbstractMatrix{Float64}, Biso::AbstractMatrix{Float64},
                                  cf_veg::CNVegCarbonFluxData;
                                  mask_soilp::AbstractVector{Bool},
                                  bounds_patch::UnitRange{Int},
                                  nvegcpool::Int, counts, dt::Real,
                                  num_actfirep::Int=0)
    dt = Float64(dt)
    begp = first(bounds_patch); endp = last(bounds_patch)
    filter_soilp = Int[p for p in bounds_patch if mask_soilp[p]]
    num_soilp = length(filter_soilp)
    num_soilp == 0 && return Xiso
    cf = cf_veg

    X = VectorType(); init_v!(X, nvegcpool, begp, endp)
    B = VectorType(); init_v!(B, nvegcpool, begp, endp)
    for fp in 1:num_soilp
        p = filter_soilp[fp]; pr = u_idx(begp, p)
        for i in 1:nvegcpool
            X.V[pr, i] = Xiso[pr, i]
            B.V[pr, i] = Biso[pr, i]
        end
    end

    _cn_veg_matrix_advance!(X, B, begp, endp, num_soilp, filter_soilp, nvegcpool,
        cf.matrix_phtransfer_patch, cf.matrix_phturnover_patch,
        cf.matrix_phtransfer_doner_patch, cf.matrix_phtransfer_receiver_patch,
        counts.ncphtrans, counts.ncphouttrans,
        cf.matrix_gmtransfer_patch, cf.matrix_gmturnover_patch,
        cf.matrix_gmtransfer_doner_patch, cf.matrix_gmtransfer_receiver_patch,
        counts.ncgmtrans, counts.ncgmouttrans,
        cf.matrix_fitransfer_patch, cf.matrix_fiturnover_patch,
        cf.matrix_fitransfer_doner_patch, cf.matrix_fitransfer_receiver_patch,
        counts.ncfitrans, counts.ncfiouttrans,
        dt, num_actfirep)

    for fp in 1:num_soilp
        p = filter_soilp[fp]; pr = u_idx(begp, p)
        for i in 1:nvegcpool
            Xiso[pr, i] = X.V[pr, i]
        end
    end
    release_v!(X); release_v!(B)
    return Xiso
end

# -----------------------------------------------------------------------------
# Veg SASU spinup-capacity accumulator + steady-state solve
# -----------------------------------------------------------------------------

"""
    CNVegMatrixSASU

Per-spin-period SASU accumulators for the vegetation matrix (the Fortran
`matrix_*alloc_acc` / `matrix_*transfer_acc` annual accumulators on the CN-veg
carbon/nitrogen state structs, plus the begin-of-year pool snapshot `X0`).

Accumulated over a SASU period (`nyr_SASU` forcing years) by
[`cn_veg_matrix_sasu_accumulate!`]; the analytic steady-state capacity is computed
at period end by [`cn_veg_matrix_sasu_capacity!`]. All gated behind
`spinup_matrixcn`. Stored `(np × nvegpool)` for the input/alloc accumulator,
`(np × nvegpool × nvegpool)` for the transfer accumulator, `(np × nvegpool)` for
the begin-of-period pool snapshot.
"""
Base.@kwdef mutable struct CNVegMatrixSASU
    alloc_acc::Matrix{Float64}    = Matrix{Float64}(undef, 0, 0)   # Σ B  (np × nveg)
    transfer_acc::Array{Float64,3} = Array{Float64,3}(undef, 0, 0, 0)  # Σ AK·X0 (np × nveg × nveg)
    X0::Matrix{Float64}           = Matrix{Float64}(undef, 0, 0)   # begin-of-period pools (np × nveg)
    nyears::Int = 0
    allocated::Bool = false
end

function cn_veg_matrix_sasu_alloc!(s::CNVegMatrixSASU, np::Int, nveg::Int)
    s.allocated && return nothing
    s.alloc_acc    = zeros(Float64, np, nveg)
    s.transfer_acc = zeros(Float64, np, nveg, nveg)
    s.X0           = zeros(Float64, np, nveg)
    s.nyears = 0
    s.allocated = true
    return nothing
end

"""
    cn_veg_matrix_sasu_save_x0!(s, X0; mask_soilp, bounds_patch, nveg, epsi=1e-30)

Save the begin-of-(SASU)-period pool snapshot `X0[pr, i]` (Fortran `*c0`/`*n0`
variables, set `max(pool, epsi)` at `is_beg_curr_year` — CNVegMatrixMod.F90:1257).
Used to normalize the accumulated transfer matrix into a per-X0 rate matrix.
"""
function cn_veg_matrix_sasu_save_x0!(s::CNVegMatrixSASU, X0::AbstractMatrix{Float64};
                                     mask_soilp::AbstractVector{Bool},
                                     bounds_patch::UnitRange{Int}, nveg::Int, epsi::Float64=1.0e-30)
    begp = first(bounds_patch)
    cn_veg_matrix_sasu_alloc!(s, size(s.X0, 1) == 0 ? last(bounds_patch) - begp + 1 : size(s.X0, 1), nveg)
    for p in bounds_patch
        mask_soilp[p] || continue
        pr = u_idx(begp, p)
        for i in 1:nveg
            s.X0[pr, i] = max(X0[pr, i], epsi)
        end
    end
    return nothing
end

"""
    cn_veg_matrix_sasu_accumulate!(s, A_dense, B; mask_soilp, bounds_patch, nveg)

Accumulate one step (or one year) into the SASU period: `alloc_acc += B`, and
`transfer_acc += A·diag(X0)` (i.e. the actual transfer FLUX matrix `AK·X`, which
is what Fortran accumulates as `matrix_*transfer_acc`). `A_dense[pr, :, :]` is the
per-patch `(I + AKall)` − I operator (the off-diagonal transfers + `-turnover`
diagonal); `B[pr, :]` is the per-patch input `alloc·input·dt`. Port of the
`matrix_*transfer_acc`/`matrix_*alloc_acc` accumulation (CNVegMatrixMod.F90 ~1734,
2056, 2740–2762 where columns are divided by `X0` at period end).
"""
function cn_veg_matrix_sasu_accumulate!(s::CNVegMatrixSASU,
        A_dense::AbstractArray{Float64,3}, B::AbstractMatrix{Float64};
        mask_soilp::AbstractVector{Bool}, bounds_patch::UnitRange{Int}, nveg::Int)
    begp = first(bounds_patch)
    @assert s.allocated "cn_veg_matrix_sasu_accumulate!: SASU not allocated (call save_x0! first)"
    for p in bounds_patch
        mask_soilp[p] || continue
        pr = u_idx(begp, p)
        for i in 1:nveg
            s.alloc_acc[pr, i] += B[pr, i]
            for j in 1:nveg
                # Accumulate the actual transfer FLUX A[i,j]·X0[j].
                s.transfer_acc[pr, i, j] += A_dense[pr, i, j] * s.X0[pr, j]
            end
        end
    end
    s.nyears += 1
    return nothing
end

"""
    cn_veg_matrix_sasu_capacity!(s; mask_soilp, bounds_patch, nveg, reset=true) -> cap

End-of-SASU-period analytic steady-state capacity. Normalizes the accumulated
transfer flux matrix back to a per-X0 rate matrix (`transfer_acc[:,:,j] / X0[j]`,
Fortran 2740–2762), then solves `cap = −A^{-1}·alloc_acc` via [`sasu_steady_state`]
(Fortran `inverse` + `-matmul`, 2845–2862). Returns `cap` `(np × nveg)`; resets the
accumulators when `reset=true`. Gated behind `spinup_matrixcn`.
"""
function cn_veg_matrix_sasu_capacity!(s::CNVegMatrixSASU;
        mask_soilp::AbstractVector{Bool}, bounds_patch::UnitRange{Int}, nveg::Int,
        reset::Bool=true)
    begp = first(bounds_patch)
    @assert s.allocated "cn_veg_matrix_sasu_capacity!: SASU not allocated"
    np = size(s.alloc_acc, 1)
    cap = zeros(Float64, np, nveg)
    for p in bounds_patch
        mask_soilp[p] || continue
        pr = u_idx(begp, p)
        A = zeros(Float64, nveg, nveg)
        for i in 1:nveg, j in 1:nveg
            A[i, j] = s.transfer_acc[pr, i, j] / s.X0[pr, j]
        end
        Xss = sasu_steady_state(A, s.alloc_acc[pr, 1:nveg])
        cap[pr, :] .= Xss
    end
    if reset
        fill!(s.alloc_acc, 0.0)
        fill!(s.transfer_acc, 0.0)
        s.nyears = 0
    end
    return cap
end
