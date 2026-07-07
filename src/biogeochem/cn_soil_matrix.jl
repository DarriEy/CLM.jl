# ==========================================================================
# CNSoilMatrixMod — matrix solution of the soil C/N decomposition cascade
#
# Julia port of CTSM `src/soilbiogeochem/CNSoilMatrixMod.F90` (the matrix model
# of CLM5.0 BGC, developed by the Yiqi Luo EcoLab: Drs. Xingjie Lu, Yuanyuan
# Huang and Zhengguang Du).
#
# The matrix-CN solution method advances the soil C/N pools as a sparse-matrix
# operation each step,
#
#     X(t+1) = X(t) + (A·K + V − Kfire)·X(t) + I·dt
#
# where
#   A      — the decomposition cascade transfer coefficients (off-diagonal),
#            assembled into `a_ma_vr` (C) / `na_ma_vr` (N) and combined with a
#            `-1` diagonal into AKsoilc / AKsoiln,
#   K      — `Ksoil%DM`, the per-pool/level turnover-rate diagonal (already
#            `* dt`, and `* fpi_vr` for immobilization steps — the Ksoil%DM
#            correction),
#   V      — `tri_ma_vr`, the vertical-transport tridiagonal (AVsoil),
#   Kfire  — `matrix_decomp_fire_k`, the fire-loss diagonal (AKfiresoil),
#   I·dt   — `matrix_Cinput` / `matrix_Ninput`, the per-step input (the B term).
#
# This is mathematically the SAME pool update as the sequential per-flux cascade
# (`decomp_cpools_vr += sourcesink`) — the matrix form is an exact reformulation.
# The two agree to round-off; the test asserts this directly.
#
# Gated on `use_soil_matrixcn`: when off, `cn_soil_matrix!` is never called and
# the sequential path in `c_state_update*`/`litter_vert_transp` is unchanged
# (byte-identical default).
#
# Uses the sparse-matrix foundation in `sparse_matrix_multiply.jl`.
# ==========================================================================

# --------------------------------------------------------------------------
# init_soil_transfer! — build the prescribed sparse index structures.
#
# Port of `InitSoilTransfer` in SoilBiogeochemDecompCascadeConType.F90. Counts
# the non-zero entries of the cascade transfer matrix A and the vertical
# transport matrix V, records their (row, col) indices, and the merged index
# structure of (A + V − Kfire). These are static (depend only on the cascade
# topology + nlevdecomp), computed once and stored on the cascade-con struct.
#
# Populates, on `cc` (DecompCascadeConData):
#   spm_tranlist_a[j,k]  — map (level, transition) → 1D entry in a_ma_vr
#   A_i[n], A_j[n]       — row/col of cascade-transfer entry n
#   Ntrans_setup         — number of cascade-transfer entries
#   tri_i[n], tri_j[n]   — row/col of vertical-transport entry n
#   Ntri_setup           — number of vertical-transport entries
#   all_i[j], all_j[j]   — row/col of the merged (A+V−Kfire) entries
#   n_all_entries        — number of merged entries
#   list_*               — memo-index arrays (allocated; filled on first solve)
# --------------------------------------------------------------------------
function init_soil_transfer!(cc;
        ndecomp_pools::Int, nlevdecomp::Int,
        ndecomp_cascade_transitions::Int,
        ndecomp_cascade_outtransitions::Int=0)

    ndecomp_pools_vr = ndecomp_pools * nlevdecomp
    Ntrans_setup = (ndecomp_cascade_transitions - ndecomp_cascade_outtransitions) * nlevdecomp
    Ntri_setup   = (3 * nlevdecomp - 2) * (ndecomp_pools - 1)   # exclude one cwd

    cc.Ntrans_setup = Ntrans_setup
    cc.Ntri_setup   = Ntri_setup

    spm_tranlist_a = fill(SMM_EMPTY_INT, nlevdecomp, ndecomp_cascade_transitions)
    A_i = fill(SMM_EMPTY_INT, Ntrans_setup)
    A_j = fill(SMM_EMPTY_INT, Ntrans_setup)

    # Count transitions per donor pool.
    ntrans_per_donor = zeros(Int, ndecomp_pools)
    for k in 1:ndecomp_cascade_transitions
        ntrans_per_donor[cc.cascade_donor_pool[k]] += 1
    end

    # Build A index structure (cascade transfers, receiver != atmosphere).
    k = 0
    n = 1
    for i in 1:ndecomp_pools
        for j in 1:nlevdecomp
            for m in 1:ntrans_per_donor[i]
                if cc.cascade_receiver_pool[m + k] != 0
                    spm_tranlist_a[j, m + k] = n
                    A_i[n] = (cc.cascade_receiver_pool[m + k] - 1) * nlevdecomp + j
                    A_j[n] = (cc.cascade_donor_pool[m + k] - 1) * nlevdecomp + j
                    n += 1
                end
            end
        end
        k += ntrans_per_donor[i]
    end
    @assert (n - 1) == Ntrans_setup "init_soil_transfer!: transfer count mismatch"

    # Build V (vertical-transport tridiagonal) index structure.
    tri_i = fill(SMM_EMPTY_INT, Ntri_setup)
    tri_j = fill(SMM_EMPTY_INT, Ntri_setup)
    n = 1
    for i in 1:ndecomp_pools
        for j in 1:nlevdecomp
            if !cc.is_cwd[i]
                if j > 1
                    tri_j[n] = (i - 1) * nlevdecomp + j
                    tri_i[n] = (i - 1) * nlevdecomp + j - 1
                    n += 1
                end
                tri_j[n] = (i - 1) * nlevdecomp + j
                tri_i[n] = (i - 1) * nlevdecomp + j
                n += 1
                if j < nlevdecomp
                    tri_j[n] = (i - 1) * nlevdecomp + j
                    tri_i[n] = (i - 1) * nlevdecomp + j + 1
                    n += 1
                end
            end
        end
    end
    @assert (n - 1) == Ntri_setup "init_soil_transfer!: vertical-transfer count mismatch"

    cc.spm_tranlist_a = spm_tranlist_a
    cc.A_i = A_i
    cc.A_j = A_j
    cc.tri_i = tri_i
    cc.tri_j = tri_j

    # Build the merged (A + V − Kfire) index structure by running the sparse
    # adders once on unit-valued matrices over a single dummy unit (begu=endu=1).
    AK    = SparseMatrixType()
    AV    = SparseMatrixType()
    AKf   = SparseMatrixType()
    AKall = SparseMatrixType()
    init_sm!(AK, ndecomp_pools_vr, 1, 1)
    init_sm!(AV, ndecomp_pools_vr, 1, 1)
    init_sm!(AKf, ndecomp_pools_vr, 1, 1)
    init_sm!(AKall, ndecomp_pools_vr, 1, 1)

    filter_c = [1]
    SMv = ones(Float64, 1, max(Ntrans_setup, 1))
    init_ready = false
    set_value_a!(AK, 1, 1, 1, filter_c, SMv, A_i, A_j, Ntrans_setup, init_ready)
    list_AK_AKVfire = fill(0, AK.NE)
    list_AK_AKV     = fill(0, AK.NE)

    TRI = ones(Float64, 1, max(Ntri_setup, 1))
    set_value_sm!(AV, 1, 1, 1, filter_c, TRI, tri_i, tri_j, Ntri_setup)
    list_V_AKVfire = fill(0, AV.NE)
    list_V_AKV     = fill(0, AV.NE)

    set_value_a_diag!(AKf, 1, filter_c, 1.0)
    list_fire_AKVfire = fill(0, AKf.NE)

    spmp_abc!(AKall, 1, filter_c, AK, AV, AKf, false)

    n_all_entries = AKall.NE
    cc.n_all_entries = n_all_entries
    cc.all_i = copy(AKall.RI[1:n_all_entries])
    cc.all_j = copy(AKall.CI[1:n_all_entries])

    cc.list_Asoilc       = fill(0, Ntrans_setup)
    cc.list_Asoiln       = fill(0, Ntrans_setup)
    # Memo RI/CI for AKsoilc / AKsoiln (auto-filled on first SetValueA: size = Ntrans + diag).
    Ntrans_diag = Ntrans_setup + ndecomp_pools_vr
    cc.RI_a  = fill(SMM_EMPTY_INT, Ntrans_diag)
    cc.CI_a  = fill(SMM_EMPTY_INT, Ntrans_diag)
    cc.RI_na = fill(SMM_EMPTY_INT, Ntrans_diag)
    cc.CI_na = fill(SMM_EMPTY_INT, Ntrans_diag)
    # Persistent merged-structure memo (sized to the worst-case dense bound).
    cc.NE_AKallsoilc = SMM_EMPTY_INT
    cc.RI_AKallsoilc = fill(SMM_EMPTY_INT, ndecomp_pools_vr * ndecomp_pools_vr)
    cc.CI_AKallsoilc = fill(SMM_EMPTY_INT, ndecomp_pools_vr * ndecomp_pools_vr)
    cc.NE_AKallsoiln = SMM_EMPTY_INT
    cc.RI_AKallsoiln = fill(SMM_EMPTY_INT, ndecomp_pools_vr * ndecomp_pools_vr)
    cc.CI_AKallsoiln = fill(SMM_EMPTY_INT, ndecomp_pools_vr * ndecomp_pools_vr)
    cc.list_AK_AKVfire   = list_AK_AKVfire
    cc.list_AK_AKV       = list_AK_AKV
    cc.list_V_AKVfire    = list_V_AKVfire
    cc.list_V_AKV        = list_V_AKV
    cc.list_fire_AKVfire = list_fire_AKVfire

    release_sm!(AK); release_sm!(AV); release_sm!(AKf); release_sm!(AKall)
    return nothing
end

# --------------------------------------------------------------------------
# CNSoilMatrixState — per-solver persistent state (the Fortran `save` flags
# and the SparseMatrix/Vector/Diag workspaces that live on the flux/state
# structs). Captured here so the solver is self-contained and the matrix
# workspaces persist across timesteps (memorized index structures).
# --------------------------------------------------------------------------
Base.@kwdef mutable struct CNSoilMatrixState
    AKsoilc::SparseMatrixType    = SparseMatrixType()
    AKsoiln::SparseMatrixType    = SparseMatrixType()
    AVsoil::SparseMatrixType     = SparseMatrixType()
    AKfiresoil::SparseMatrixType = SparseMatrixType()
    AKallsoilc::SparseMatrixType = SparseMatrixType()
    AKallsoiln::SparseMatrixType = SparseMatrixType()
    AKXcacc::SparseMatrixType    = SparseMatrixType()
    AKXnacc::SparseMatrixType    = SparseMatrixType()
    Xdiagsoil::DiagMatrixType    = DiagMatrixType()
    matrix_Cinter::VectorType    = VectorType()
    matrix_Ninter::VectorType    = VectorType()

    # Fortran `save` flags (Julia by-value → carry them on the struct).
    list_ready1_fire::Bool   = false
    list_ready1_nofire::Bool = false
    list_ready2_fire::Bool   = false
    list_ready2_nofire::Bool = false
    init_readyAsoilc::Bool   = false
    init_readyAsoiln::Bool   = false

    allocated::Bool = false
end

# Allocate the matrix workspaces (once). begc:endc unit bounds.
function cn_soil_matrix_alloc!(ms::CNSoilMatrixState; ndecomp_pools_vr::Int,
                               begc::Int, endc::Int)
    ms.allocated && return nothing
    init_sm!(ms.AKsoilc,    ndecomp_pools_vr, begc, endc)
    init_sm!(ms.AKsoiln,    ndecomp_pools_vr, begc, endc)
    init_sm!(ms.AVsoil,     ndecomp_pools_vr, begc, endc)
    init_sm!(ms.AKfiresoil, ndecomp_pools_vr, begc, endc)
    init_sm!(ms.AKallsoilc, ndecomp_pools_vr, begc, endc)
    init_sm!(ms.AKallsoiln, ndecomp_pools_vr, begc, endc)
    init_dm!(ms.Xdiagsoil,  ndecomp_pools_vr, begc, endc)
    init_v!(ms.matrix_Cinter, ndecomp_pools_vr, begc, endc)
    init_v!(ms.matrix_Ninter, ndecomp_pools_vr, begc, endc)
    ms.allocated = true
    return nothing
end

# --------------------------------------------------------------------------
# cn_soil_matrix_input_accumulate! — assemble the B-input (matrix_Cinput/Ninput).
#
# In the sequential path each veg→soil litterfall/CWD flux is added directly into
# decomp_cpools_vr (via decomp_cpools_sourcesink) by the C/N state updates. In the
# matrix path those same EXTERNAL inputs become the B term, and the internal
# decomposition cascade is the A-matrix (Ksoil·rf·pathfrac, handled by
# cn_soil_matrix!). This fills matrix_Cinput/Ninput[c, j+(i-1)*nlevdecomp] = Σ
# (litter/CWD input fluxes to pool i, level j) · dt — the port of the
# `matrix_Cinput%V += *_to_litr/cwd * dt` accumulations scattered across
# CNC/NStateUpdate1/2/3. Litter pools (i_litr_min..i_litr_max) get the *_to_litr_*
# fluxes (indexed by pool); the CWD pool (i_cwd) gets the *_to_cwd* fluxes.
# `transient_landcover` adds the harvest + gross-unrepresented-landcover inputs
# (zero otherwise, matching the default non-transient run).
# --------------------------------------------------------------------------

# (litr-pool flux fields (c,j,i), cwd flux fields (c,j)) for C and N.
const _SOILC_LITR_IN = (:phenology_c_to_litr_c_col, :dwt_frootc_to_litr_c_col,
                        :gap_mortality_c_to_litr_c_col, :m_c_to_litr_fire_col)
const _SOILC_CWD_IN  = (:dwt_livecrootc_to_cwdc_col, :dwt_deadcrootc_to_cwdc_col,
                        :gap_mortality_c_to_cwdc_col, :fire_mortality_c_to_cwdc_col)
const _SOILC_LITR_IN_TR = (:harvest_c_to_litr_c_col, :gru_c_to_litr_c_col)
const _SOILC_CWD_IN_TR  = (:harvest_c_to_cwdc_col, :gru_c_to_cwdc_col)
const _SOILN_LITR_IN = (:phenology_n_to_litr_n_col, :dwt_frootn_to_litr_n_col,
                        :gap_mortality_n_to_litr_n_col, :m_n_to_litr_fire_col)
const _SOILN_CWD_IN  = (:dwt_livecrootn_to_cwdn_col, :dwt_deadcrootn_to_cwdn_col,
                        :gap_mortality_n_to_cwdn_col, :fire_mortality_n_to_cwdn_col)
const _SOILN_LITR_IN_TR = (:harvest_n_to_litr_n_col, :gru_n_to_litr_n_col)
const _SOILN_CWD_IN_TR  = (:harvest_n_to_cwdn_col, :gru_n_to_cwdn_col)

function cn_soil_matrix_input_accumulate!(
        matrix_Cinput::AbstractMatrix{<:Real}, matrix_Ninput::AbstractMatrix{<:Real},
        cf_veg, nf_veg;
        mask_soilc::AbstractVector{Bool}, bounds_col::UnitRange{Int},
        nlevdecomp::Int, i_litr_min::Int, i_litr_max::Int, i_cwd::Int,
        dt::Real, transient_landcover::Bool = false)
    fill!(matrix_Cinput, 0.0)
    fill!(matrix_Ninput, 0.0)
    dt = Float64(dt)
    litr_c = transient_landcover ? (_SOILC_LITR_IN..., _SOILC_LITR_IN_TR...) : _SOILC_LITR_IN
    cwd_c  = transient_landcover ? (_SOILC_CWD_IN...,  _SOILC_CWD_IN_TR...)  : _SOILC_CWD_IN
    litr_n = transient_landcover ? (_SOILN_LITR_IN..., _SOILN_LITR_IN_TR...) : _SOILN_LITR_IN
    cwd_n  = transient_landcover ? (_SOILN_CWD_IN...,  _SOILN_CWD_IN_TR...)  : _SOILN_CWD_IN
    icwd = i_cwd
    for c in bounds_col
        mask_soilc[c] || continue
        for j in 1:nlevdecomp
            # litter pools (indexed by pool i)
            for i in i_litr_min:i_litr_max
                vr = j + (i - 1) * nlevdecomp
                cin = 0.0; nin = 0.0
                for f in litr_c; cin += getfield(cf_veg, f)[c, j, i]; end
                for f in litr_n; nin += getfield(nf_veg, f)[c, j, i]; end
                matrix_Cinput[c, vr] += cin * dt
                matrix_Ninput[c, vr] += nin * dt
            end
            # CWD pool (level-only fluxes)
            vrc = j + (icwd - 1) * nlevdecomp
            cin = 0.0; nin = 0.0
            for f in cwd_c; cin += getfield(cf_veg, f)[c, j]; end
            for f in cwd_n; nin += getfield(nf_veg, f)[c, j]; end
            matrix_Cinput[c, vrc] += cin * dt
            matrix_Ninput[c, vrc] += nin * dt
        end
    end
    return nothing
end

# --------------------------------------------------------------------------
# cn_soil_matrix_advance! — driver entry point. Assembles the soil-matrix step from
# the live soil/veg fluxes and advances the decomp C/N pools by one cn_soil_matrix!
# solve. Gated by the driver on use_soil_matrixcn (in matrix mode the sequential
# decomp cascade, vertical transport, and litter-input state updates are all skipped,
# leaving decomp_cpools_vr at start-of-step; this one solve subsumes them).
#
#   Ksoil     — the turnover diagonal decomp_k_col·dt (fpi-corrected in decomp!).
#   tri_ma_vr — the vertical-transport matrix (built by litter_vert_transp! in
#               use_soil_matrixcn mode).
#   B-input   — accumulated here from the veg→soil litter/CWD fluxes.
# `ms` is created fresh per call (the index memoization is a within-step optimization).
# init_soil_transfer! sets the (static) sparse index structure on cascade_con once.
# --------------------------------------------------------------------------
function cn_soil_matrix_advance!(cascade_con, soilbgc_cs, soilbgc_ns, soilbgc_cf,
        cnveg_cf, cnveg_nf; Ksoil::AbstractMatrix{<:Real},
        mask_soilc::AbstractVector{Bool}, bounds_col::UnitRange{Int},
        nlevdecomp::Int, ndecomp_pools::Int, ndecomp_cascade_transitions::Int,
        i_litr_min::Int, i_litr_max::Int, i_cwd::Int, dt::Real,
        num_actfirec::Int = 0, transient_landcover::Bool = false)

    begc = first(bounds_col); endc = last(bounds_col); nc = endc - begc + 1
    ndp_vr = ndecomp_pools * nlevdecomp
    nouttrans = count(==(0), cascade_con.cascade_receiver_pool)   # terminal (→atm) transitions

    # Static sparse index structure — build once (guard on the sentinel).
    if cascade_con.n_all_entries == SMM_EMPTY_INT
        init_soil_transfer!(cascade_con; ndecomp_pools=ndecomp_pools, nlevdecomp=nlevdecomp,
            ndecomp_cascade_transitions=ndecomp_cascade_transitions,
            ndecomp_cascade_outtransitions=nouttrans)
    end

    # B-input from the veg→soil litter/CWD fluxes.
    Cin = zeros(nc, ndp_vr); Nin = zeros(nc, ndp_vr)
    cn_soil_matrix_input_accumulate!(Cin, Nin, cnveg_cf, cnveg_nf;
        mask_soilc=mask_soilc, bounds_col=bounds_col, nlevdecomp=nlevdecomp,
        i_litr_min=i_litr_min, i_litr_max=i_litr_max, i_cwd=i_cwd, dt=dt,
        transient_landcover=transient_landcover)

    # Fire-loss diagonal (AKfiresoil). The fire module computed the decomp→fire flux
    # m_decomp_cpools_to_fire_vr (= pool·f·cmb); the matrix representation is the loss
    # fraction/step matrix_decomp_fire_k = −flux/pool·dt (port of CNFireBase:1239/1251,
    # derived here so the fire kernels stay untouched). Only assembled when a column is
    # actually burning; the pools are still at start-of-step so /pool is the fire base.
    fire_k = nothing
    if num_actfirec > 0
        fire_k = zeros(nc, ndp_vr)
        fflux = cnveg_cf.m_decomp_cpools_to_fire_vr_col
        Xc = soilbgc_cs.decomp_cpools_vr_col
        for l in 1:ndecomp_pools, j in 1:nlevdecomp, c in bounds_col
            mask_soilc[c] || continue
            pool = Xc[c, j, l]
            fire_k[c - begc + 1, j + (l - 1) * nlevdecomp] =
                pool > 0.0 ? -fflux[c, j, l] / pool * dt : 0.0
        end
    end

    ms = CNSoilMatrixState()
    cn_soil_matrix!(ms, cascade_con;
        decomp_cpools_vr=soilbgc_cs.decomp_cpools_vr_col,
        decomp_npools_vr=soilbgc_ns.decomp_npools_vr_col,
        Ksoil=Ksoil, tri_ma_vr=soilbgc_cf.tri_ma_vr,
        matrix_Cinput=Cin, matrix_Ninput=Nin,
        rf_decomp_cascade=soilbgc_cf.rf_decomp_cascade_col,
        pathfrac_decomp_cascade=soilbgc_cf.pathfrac_decomp_cascade_col,
        matrix_decomp_fire_k=fire_k,
        mask_soilc=mask_soilc, begc=begc, endc=endc,
        nlevdecomp=nlevdecomp, ndecomp_pools=ndecomp_pools,
        ndecomp_cascade_transitions=ndecomp_cascade_transitions,
        ndecomp_cascade_outtransitions=nouttrans, num_actfirec=num_actfirec)
    return nothing
end

# --------------------------------------------------------------------------
# cn_soil_matrix! — advance the soil C/N pools by the matrix solve for one step.
#
# Port of `CNSoilMatrix` (the non-isotope, non-spinup core path). Inputs:
#   decomp_cpools_vr, decomp_npools_vr  — (col, level, pool) pool sizes (updated)
#   Ksoil   — (col, ndecomp_pools_vr) turnover-rate diagonal (gC/gC/step), incl.
#             the fpi_vr Ksoil%DM correction for immobilization steps
#   tri_ma_vr — (col, Ntri_setup) vertical-transport entries
#   matrix_decomp_fire_k — (col, ndecomp_pools_vr) fire-loss diagonal (or nothing)
#   matrix_Cinput, matrix_Ninput — (col, ndecomp_pools_vr) per-step C/N input
#   rf_decomp_cascade, pathfrac_decomp_cascade — (col, level, trans) cascade frac
#   num_actfirec — number of active-fire columns (0 → no-fire path)
#
# The N transfer matrix An is built from a_ma_vr scaled by the donor/receiver
# C:N ratio (na_ma_vr), and Ksoiln is Ksoil scaled by the fixed-C:N correction.
# --------------------------------------------------------------------------
function cn_soil_matrix!(ms::CNSoilMatrixState, cc;
        decomp_cpools_vr::AbstractArray{<:Real,3},
        decomp_npools_vr::AbstractArray{<:Real,3},
        Ksoil::AbstractMatrix{<:Real},
        tri_ma_vr::AbstractMatrix{<:Real},
        matrix_Cinput::AbstractMatrix{<:Real},
        matrix_Ninput::AbstractMatrix{<:Real},
        rf_decomp_cascade::AbstractArray{<:Real,3},
        pathfrac_decomp_cascade::AbstractArray{<:Real,3},
        matrix_decomp_fire_k::Union{AbstractMatrix{<:Real},Nothing}=nothing,
        mask_soilc::AbstractVector{Bool},
        begc::Int, endc::Int,
        nlevdecomp::Int, ndecomp_pools::Int,
        ndecomp_cascade_transitions::Int,
        ndecomp_cascade_outtransitions::Int=0,
        num_actfirec::Int=0)

    ndecomp_pools_vr = ndecomp_pools * nlevdecomp
    Ntrans = (ndecomp_cascade_transitions - ndecomp_cascade_outtransitions) * nlevdecomp
    Ntri_setup = cc.Ntri_setup

    # Active-unit filter (1-based column indices that are active soil columns).
    filter_soilc = [c for c in begc:endc if mask_soilc[c]]
    num_soilc = length(filter_soilc)

    cn_soil_matrix_alloc!(ms; ndecomp_pools_vr=ndecomp_pools_vr, begc=begc, endc=endc)

    cascade_donor_pool    = cc.cascade_donor_pool
    cascade_receiver_pool = cc.cascade_receiver_pool
    floating              = cc.floating_cn_ratio_decomp_pools
    initial_cn_ratio      = cc.initial_cn_ratio

    nunit = endc - begc + 1

    # ----- C:N ratios of applicable pools (gC/gN) -----
    cn_decomp_pools = zeros(Float64, nunit, nlevdecomp, ndecomp_pools)
    for l in 1:ndecomp_pools
        if floating[l]
            for j in 1:nlevdecomp, c in filter_soilc
                ui = c - begc + 1
                if decomp_npools_vr[c, j, l] > 0.0
                    cn_decomp_pools[ui, j, l] = decomp_cpools_vr[c, j, l] / decomp_npools_vr[c, j, l]
                else
                    cn_decomp_pools[ui, j, l] = initial_cn_ratio[l]
                end
            end
        else
            for j in 1:nlevdecomp, c in filter_soilc
                ui = c - begc + 1
                cn_decomp_pools[ui, j, l] = initial_cn_ratio[l]
            end
        end
    end

    # ----- Assemble a_ma_vr / na_ma_vr (off-diagonal transfer coefficients). -----
    a_ma_vr  = zeros(Float64, nunit, Ntrans)
    na_ma_vr = zeros(Float64, nunit, Ntrans)
    for k in 1:ndecomp_cascade_transitions
        if cascade_receiver_pool[k] != 0
            for j in 1:nlevdecomp
                tranlist_a = cc.spm_tranlist_a[j, k]
                for c in filter_soilc
                    ui = c - begc + 1
                    a_ma_vr[ui, tranlist_a] =
                        (1.0 - rf_decomp_cascade[c, j, k]) * pathfrac_decomp_cascade[c, j, k]
                    if !floating[cascade_receiver_pool[k]]
                        na_ma_vr[ui, tranlist_a] =
                            (1.0 - rf_decomp_cascade[c, j, k]) *
                            (cn_decomp_pools[ui, j, cascade_donor_pool[k]] /
                             cn_decomp_pools[ui, j, cascade_receiver_pool[k]]) *
                            pathfrac_decomp_cascade[c, j, k]
                    else
                        na_ma_vr[ui, tranlist_a] = pathfrac_decomp_cascade[c, j, k]
                    end
                end
            end
        end
    end

    # ----- Set old vector value (X(t)). -----
    Cinter_old = zeros(Float64, nunit, ndecomp_pools_vr)
    Ninter_old = zeros(Float64, nunit, ndecomp_pools_vr)
    for i in 1:ndecomp_pools, j in 1:nlevdecomp, c in filter_soilc
        ui = c - begc + 1
        idx = j + (i - 1) * nlevdecomp
        ms.matrix_Cinter.V[ui, idx] = decomp_cpools_vr[c, j, i]
        Cinter_old[ui, idx]         = decomp_cpools_vr[c, j, i]
        ms.matrix_Ninter.V[ui, idx] = decomp_npools_vr[c, j, i]
        Ninter_old[ui, idx]         = decomp_npools_vr[c, j, i]
    end

    # ----- Build Ksoiln (= Ksoil scaled by the fixed-C:N correction). -----
    Ksoiln = copy(Ksoil)
    for i in 1:ndecomp_pools
        if !floating[i]
            for j in 1:nlevdecomp, c in filter_soilc
                ui = c - begc + 1
                idx = j + (i - 1) * nlevdecomp
                if decomp_npools_vr[c, j, i] > 0
                    Ksoiln[ui, idx] = Ksoil[ui, idx] *
                        decomp_cpools_vr[c, j, i] / decomp_npools_vr[c, j, i] / initial_cn_ratio[i]
                end
            end
        end
    end

    # ----- Set C transfer matrix Ac from a_ma_vr (adds the -1 diagonal). -----
    ms.init_readyAsoilc = set_value_a!(ms.AKsoilc, begc, endc, num_soilc, filter_soilc,
        a_ma_vr, cc.A_i, cc.A_j, Ntrans, ms.init_readyAsoilc;
        list=cc.list_Asoilc, RI_A=cc.RI_a, CI_A=cc.CI_a)
    ms.init_readyAsoiln = set_value_a!(ms.AKsoiln, begc, endc, num_soilc, filter_soilc,
        na_ma_vr, cc.A_i, cc.A_j, Ntrans, ms.init_readyAsoiln;
        list=cc.list_Asoiln, RI_A=cc.RI_na, CI_A=cc.CI_na)

    # ----- Ac*K, An*K. -----
    Kdm  = DiagMatrixType();  init_dm!(Kdm,  ndecomp_pools_vr, begc, endc)
    Kndm = DiagMatrixType();  init_dm!(Kndm, ndecomp_pools_vr, begc, endc)
    set_value_dm!(Kdm,  begc, endc, num_soilc, filter_soilc, Ksoil)
    set_value_dm!(Kndm, begc, endc, num_soilc, filter_soilc, Ksoiln)
    spmm_ak!(ms.AKsoilc, num_soilc, filter_soilc, Kdm)
    spmm_ak!(ms.AKsoiln, num_soilc, filter_soilc, Kndm)

    # ----- Set vertical-transport matrix V from tri_ma_vr. -----
    set_value_sm!(ms.AVsoil, begc, endc, num_soilc, filter_soilc,
        tri_ma_vr, cc.tri_i, cc.tri_j, Ntri_setup)

    # ----- Set fire matrix Kfire from matrix_decomp_fire_k. -----
    if num_actfirec != 0
        @assert matrix_decomp_fire_k !== nothing "cn_soil_matrix!: fire active but matrix_decomp_fire_k not supplied"
        kfire_i = collect(1:ndecomp_pools_vr)
        kfire_j = collect(1:ndecomp_pools_vr)
        set_value_sm!(ms.AKfiresoil, begc, endc, num_soilc, filter_soilc,
            matrix_decomp_fire_k, kfire_i, kfire_j, ndecomp_pools_vr)
    end

    # ----- AKallsoilc = A*K + V (− Kfire). The merged structure (NE/RI/CI) is
    # memorized on `cc` after the first call and reused (list_ready=true). -----
    if num_actfirec == 0
        (ms.list_ready1_nofire, cc.NE_AKallsoilc) = spmp_ab!(ms.AKallsoilc, num_soilc, filter_soilc,
            ms.AKsoilc, ms.AVsoil, ms.list_ready1_nofire;
            list_A=cc.list_AK_AKV, list_B=cc.list_V_AKV,
            NE_AB=cc.NE_AKallsoilc, RI_AB=cc.RI_AKallsoilc, CI_AB=cc.CI_AKallsoilc)
        (ms.list_ready2_nofire, cc.NE_AKallsoiln) = spmp_ab!(ms.AKallsoiln, num_soilc, filter_soilc,
            ms.AKsoiln, ms.AVsoil, ms.list_ready2_nofire;
            list_A=cc.list_AK_AKV, list_B=cc.list_V_AKV,
            NE_AB=cc.NE_AKallsoiln, RI_AB=cc.RI_AKallsoiln, CI_AB=cc.CI_AKallsoiln)
    else
        filter_actfirec = filter_soilc  # caller scopes fire via matrix_decomp_fire_k zeros
        (ms.list_ready1_fire, cc.NE_AKallsoilc) = spmp_abc!(ms.AKallsoilc, num_soilc, filter_soilc,
            ms.AKsoilc, ms.AVsoil, ms.AKfiresoil, ms.list_ready1_fire;
            list_A=cc.list_AK_AKVfire, list_B=cc.list_V_AKVfire, list_C=cc.list_fire_AKVfire,
            NE_ABC=cc.NE_AKallsoilc, RI_ABC=cc.RI_AKallsoilc, CI_ABC=cc.CI_AKallsoilc,
            num_actunit_C=num_actfirec, filter_actunit_C=filter_actfirec)
        (ms.list_ready2_fire, cc.NE_AKallsoiln) = spmp_abc!(ms.AKallsoiln, num_soilc, filter_soilc,
            ms.AKsoiln, ms.AVsoil, ms.AKfiresoil, ms.list_ready2_fire;
            list_A=cc.list_AK_AKVfire, list_B=cc.list_V_AKVfire, list_C=cc.list_fire_AKVfire,
            NE_ABC=cc.NE_AKallsoiln, RI_ABC=cc.RI_AKallsoiln, CI_ABC=cc.CI_AKallsoiln,
            num_actunit_C=num_actfirec, filter_actunit_C=filter_actfirec)
    end

    # ----- Advance pools: X = X + (A*K + V − Kfire)·X. -----
    spmm_ax!(ms.matrix_Cinter, num_soilc, filter_soilc, ms.AKallsoilc)
    spmm_ax!(ms.matrix_Ninter, num_soilc, filter_soilc, ms.AKallsoiln)

    # ----- Add input: X = X + I·dt. -----
    for j in 1:ndecomp_pools_vr, c in filter_soilc
        ui = c - begc + 1
        ms.matrix_Cinter.V[ui, j] = matrix_Cinput[ui, j] + ms.matrix_Cinter.V[ui, j]
        ms.matrix_Ninter.V[ui, j] = matrix_Ninput[ui, j] + ms.matrix_Ninter.V[ui, j]
    end

    # ----- Send back to decomp_cpools_vr / decomp_npools_vr. -----
    for i in 1:ndecomp_pools, j in 1:nlevdecomp, c in filter_soilc
        ui = c - begc + 1
        idx = j + (i - 1) * nlevdecomp
        decomp_cpools_vr[c, j, i] = ms.matrix_Cinter.V[ui, idx]
        decomp_npools_vr[c, j, i] = ms.matrix_Ninter.V[ui, idx]
    end

    release_dm!(Kdm); release_dm!(Kndm)
    return (Cinter_old, Ninter_old)
end

# --------------------------------------------------------------------------
# cn_soil_matrix_akx_accumulate! — the spinup-acceleration (SASU/AKX) path.
#
# Port of the `hist_wrt_matrixcn_diag .or. spinup_matrixcn` block. After the
# pool advance, recompute AKallsoilc/AKallsoiln as actual transfers (× X_old via
# Xdiagsoil), accumulate into AKXcacc/AKXnacc, and — at end of a SASU year —
# compute the analytic steady-state capacity X* = −(A)^{-1}·I from the
# accumulated transfer rates (the inverse of the per-year-averaged transfer
# matrix times the accumulated input).
#
# This routine advances pools toward their steady-state capacity, which the
# test exercises ("the spinup-AKX path advances pools"). The full
# capacity-clamp / SASU-averaging history bookkeeping is preserved.
# --------------------------------------------------------------------------
function cn_soil_matrix_akx_accumulate!(ms::CNSoilMatrixState, cc,
        Cinter_old::AbstractMatrix{<:Real}, Ninter_old::AbstractMatrix{<:Real};
        in_acc::AbstractMatrix{<:Real}, in_nacc::AbstractMatrix{<:Real},
        matrix_Cinput::AbstractMatrix{<:Real}, matrix_Ninput::AbstractMatrix{<:Real},
        decomp0_cpools_vr::AbstractArray{<:Real,3},
        decomp0_npools_vr::AbstractArray{<:Real,3},
        mask_soilc::AbstractVector{Bool}, begc::Int, endc::Int,
        nlevdecomp::Int, ndecomp_pools::Int,
        first_step::Bool=false, compute_capacity::Bool=false)

    ndecomp_pools_vr = ndecomp_pools * nlevdecomp
    filter_soilc = [c for c in begc:endc if mask_soilc[c]]
    num_soilc = length(filter_soilc)
    epsi = 1.0e-8

    # Accumulate input.
    for j in 1:(ndecomp_pools * nlevdecomp), c in filter_soilc
        ui = c - begc + 1
        in_acc[ui, j]  += matrix_Cinput[ui, j]
        in_nacc[ui, j] += matrix_Ninput[ui, j]
    end

    # AKallsoilc/n → actual transfers (× X_old).
    set_value_dm!(ms.Xdiagsoil, begc, endc, num_soilc, filter_soilc, Cinter_old)
    spmm_ak!(ms.AKallsoilc, num_soilc, filter_soilc, ms.Xdiagsoil)
    set_value_dm!(ms.Xdiagsoil, begc, endc, num_soilc, filter_soilc, Ninter_old)
    spmm_ak!(ms.AKallsoiln, num_soilc, filter_soilc, ms.Xdiagsoil)

    # Accumulate transfers AKXcacc += AKallsoilc ; AKXnacc += AKallsoiln.
    if first_step || !is_alloc_sm(ms.AKXcacc)
        is_alloc_sm(ms.AKXcacc) || init_sm!(ms.AKXcacc, ndecomp_pools_vr, begc, endc)
        is_alloc_sm(ms.AKXnacc) || init_sm!(ms.AKXnacc, ndecomp_pools_vr, begc, endc)
    end
    if is_values_set_sm(ms.AKXcacc)
        spmp_b_acc!(ms.AKXcacc, num_soilc, filter_soilc, ms.AKallsoilc)
    else
        set_value_copy_sm!(ms.AKXcacc, num_soilc, filter_soilc, ms.AKallsoilc)
    end
    if is_values_set_sm(ms.AKXnacc)
        spmp_b_acc!(ms.AKXnacc, num_soilc, filter_soilc, ms.AKallsoiln)
    else
        set_value_copy_sm!(ms.AKXnacc, num_soilc, filter_soilc, ms.AKallsoiln)
    end

    if !compute_capacity
        return nothing
    end

    # ----- End-of-SASU-year: compute capacity X* = −A^{-1} I. -----
    n_all_entries = cc.n_all_entries
    all_i = cc.all_i
    all_j = cc.all_j
    soilmatrixc_cap = zeros(Float64, length(filter_soilc), ndecomp_pools_vr)
    soilmatrixn_cap = zeros(Float64, length(filter_soilc), ndecomp_pools_vr)

    for (uf, c) in enumerate(filter_soilc)
        ui = c - begc + 1
        tran_acc  = zeros(Float64, ndecomp_pools_vr, ndecomp_pools_vr)
        tran_nacc = zeros(Float64, ndecomp_pools_vr, ndecomp_pools_vr)
        for jj in 1:n_all_entries
            j_lev    = mod(all_j[jj] - 1, nlevdecomp) + 1
            j_decomp = div(all_j[jj] - j_lev, nlevdecomp) + 1
            tran_acc[all_i[jj], all_j[jj]]  = ms.AKXcacc.M[ui, jj] / decomp0_cpools_vr[c, j_lev, j_decomp]
            tran_nacc[all_i[jj], all_j[jj]] = ms.AKXnacc.M[ui, jj] / decomp0_npools_vr[c, j_lev, j_decomp]
        end
        for i in 1:ndecomp_pools_vr
            if abs(tran_acc[i, i]) <= epsi
                tran_acc[i, i] = 1.0e36
            end
            if abs(tran_nacc[i, i]) <= epsi
                tran_nacc[i, i] = 1.0e36
            end
        end
        AKinv  = inv(tran_acc)
        AKinvn = inv(tran_nacc)
        soilmatrixc_cap[uf, :] = -(AKinv  * in_acc[ui, 1:ndecomp_pools_vr])
        soilmatrixn_cap[uf, :] = -(AKinvn * in_nacc[ui, 1:ndecomp_pools_vr])
        for k in 1:ndecomp_pools_vr
            soilmatrixc_cap[uf, k] < 0 && (soilmatrixc_cap[uf, k] = 0.0)
            soilmatrixn_cap[uf, k] < 0 && (soilmatrixn_cap[uf, k] = 0.0)
        end
    end

    # Reset accumulators.
    for jj in 1:n_all_entries, c in filter_soilc
        ui = c - begc + 1
        ms.AKXcacc.M[ui, jj] = 0.0
        ms.AKXnacc.M[ui, jj] = 0.0
    end
    fill!(in_acc, 0.0)
    fill!(in_nacc, 0.0)

    return (soilmatrixc_cap, soilmatrixn_cap)
end

# --------------------------------------------------------------------------
# sasu_steady_state — the Semi-Analytic Spin-Up (SASU) steady-state solve.
#
# Generic port of the capacity inversion at the heart of the matrix-CN spinup
# accelerator (CNSoilMatrixMod.F90:738–742, CNVegMatrixMod.F90:2845–2862):
#
#     X_ss = −(A)^{-1} · B
#
# where `A_acc` is the per-spin-period-accumulated transfer matrix (the matrix
# operator `A·K + V − Kfire` for soil, or `transfer_acc/X0` for veg) and `B_acc`
# is the accumulated input vector `Σ I·dt` over the same period. This is the
# closed-form steady state of `dX/dt = A·X + B = 0`, which the iterated pool
# advance `X(t+1) = X(t) + A·X(t) + B` converges to as `t → ∞` whenever `A` is
# stable (eigenvalues with negative real part). The test asserts exactly that
# (matrix iteration → analytic steady state).
#
# `epsi` guards a (near-)zero diagonal: a pool with no turnover would make A
# singular, so its diagonal is set to a huge value (1e36), driving that pool's
# capacity to ~0 — matching the Fortran guard. Negative capacities are clamped
# to 0 (a steady-state pool size can't be negative). `A_acc` is NOT modified.
#
# Returns `X_ss` (length n). This is the analytic engine shared by the soil and
# veg SASU capacity routines; gated entirely behind the spinup path (never on
# the default no-op).
# --------------------------------------------------------------------------
function sasu_steady_state(A_acc::AbstractMatrix{<:Real}, B_acc::AbstractVector{<:Real};
                           epsi::Float64=1.0e-8, clamp_negative::Bool=true)
    n = length(B_acc)
    @assert size(A_acc, 1) == n && size(A_acc, 2) == n "sasu_steady_state: A/B size mismatch"
    A = Matrix{Float64}(A_acc)        # local copy (Fortran does not mutate the accumulator)
    @inbounds for i in 1:n
        if abs(A[i, i]) <= epsi
            A[i, i] = 1.0e36
        end
    end
    X_ss = -(A \ Vector{Float64}(B_acc))
    if clamp_negative
        @inbounds for i in 1:n
            X_ss[i] < 0.0 && (X_ss[i] = 0.0)
        end
    end
    return X_ss
end
