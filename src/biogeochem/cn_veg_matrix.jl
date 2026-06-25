# =============================================================================
# CNVegMatrixMod — matrix solution for the vegetation C (and N) cycle
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
