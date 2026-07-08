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
                               begc::Int, endc::Int, ref=nothing, FT::Type=Float64)
    ms.allocated && return nothing
    init_sm!(ms.AKsoilc,    ndecomp_pools_vr, begc, endc; ref=ref, FT=FT)
    init_sm!(ms.AKsoiln,    ndecomp_pools_vr, begc, endc; ref=ref, FT=FT)
    init_sm!(ms.AVsoil,     ndecomp_pools_vr, begc, endc; ref=ref, FT=FT)
    init_sm!(ms.AKfiresoil, ndecomp_pools_vr, begc, endc; ref=ref, FT=FT)
    init_sm!(ms.AKallsoilc, ndecomp_pools_vr, begc, endc; ref=ref, FT=FT)
    init_sm!(ms.AKallsoiln, ndecomp_pools_vr, begc, endc; ref=ref, FT=FT)
    init_dm!(ms.Xdiagsoil,  ndecomp_pools_vr, begc, endc; ref=ref, FT=FT)
    init_v!(ms.matrix_Cinter, ndecomp_pools_vr, begc, endc; ref=ref, FT=FT)
    init_v!(ms.matrix_Ninter, ndecomp_pools_vr, begc, endc; ref=ref, FT=FT)
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
# cn_driver-LOCAL litter/CWD inputs (phenology/gap/fire, computed within cn_driver).
# The dwt inputs are NOT here — they are dyn_subgrid-local and enter via the persistent
# matrix_Cinput_col field (see _csu_dyn_col_kernel!), added in cn_soil_matrix_advance!.
const _SOILC_LITR_IN = (:phenology_c_to_litr_c_col,
                        :gap_mortality_c_to_litr_c_col, :m_c_to_litr_fire_col)
const _SOILC_CWD_IN  = (:gap_mortality_c_to_cwdc_col, :fire_mortality_c_to_cwdc_col)
const _SOILC_LITR_IN_TR = (:harvest_c_to_litr_c_col, :gru_c_to_litr_c_col)
const _SOILC_CWD_IN_TR  = (:harvest_c_to_cwdc_col, :gru_c_to_cwdc_col)
const _SOILN_LITR_IN = (:phenology_n_to_litr_n_col,
                        :gap_mortality_n_to_litr_n_col, :m_n_to_litr_fire_col)
const _SOILN_CWD_IN  = (:gap_mortality_n_to_cwdn_col, :fire_mortality_n_to_cwdn_col)
const _SOILN_LITR_IN_TR = (:harvest_n_to_litr_n_col, :gru_n_to_litr_n_col)
const _SOILN_CWD_IN_TR  = (:harvest_n_to_cwdn_col, :gru_n_to_cwdn_col)

# Device-view bundles for the soil B-input accumulate (litter 3D + CWD 2D veg-flux fields).
Base.@kwdef struct _SoilInOut{M}
    cin::M; nin::M
end
Adapt.@adapt_structure _SoilInOut
Base.@kwdef struct _SoilInIn{M3,M2}
    lc1::M3; lc2::M3; lc3::M3; lc4::M3; lc5::M3     # litter C (phenology/gap/fire + harvest/gru)
    ln1::M3; ln2::M3; ln3::M3; ln4::M3; ln5::M3     # litter N
    wc1::M2; wc2::M2; wc3::M2; wc4::M2              # cwd C (gap/fire + harvest/gru)
    wn1::M2; wn2::M2; wn3::M2; wn4::M2              # cwd N
end
Adapt.@adapt_structure _SoilInIn

# One thread per column; sums the litter + CWD contributions into cin/nin (pre-zeroed).
# Byte-identical to the host loops (each [c,vr] written once). transient adds harvest/gru.
@kernel function _soil_input_kernel!(out::_SoilInOut, in::_SoilInIn, @Const(mask),
        nlev::Int, ilmin::Int, ilmax::Int, icwd::Int, transient::Bool, dt, cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && mask[c]
        cin = out.cin; nin = out.nin
        for j in 1:nlev
            for i in ilmin:ilmax
                vr = j + (i - 1) * nlev
                cc = in.lc1[c, j, i] + in.lc2[c, j, i] + in.lc3[c, j, i]
                nn = in.ln1[c, j, i] + in.ln2[c, j, i] + in.ln3[c, j, i]
                if transient
                    cc += in.lc4[c, j, i] + in.lc5[c, j, i]
                    nn += in.ln4[c, j, i] + in.ln5[c, j, i]
                end
                cin[c, vr] = cin[c, vr] + cc * dt
                nin[c, vr] = nin[c, vr] + nn * dt
            end
            vrc = j + (icwd - 1) * nlev
            cc = in.wc1[c, j] + in.wc2[c, j]; nn = in.wn1[c, j] + in.wn2[c, j]
            if transient
                cc += in.wc3[c, j] + in.wc4[c, j]; nn += in.wn3[c, j] + in.wn4[c, j]
            end
            cin[c, vrc] = cin[c, vrc] + cc * dt
            nin[c, vrc] = nin[c, vrc] + nn * dt
        end
    end
end

function cn_soil_matrix_input_accumulate!(
        matrix_Cinput::AbstractMatrix{<:Real}, matrix_Ninput::AbstractMatrix{<:Real},
        cf_veg, nf_veg;
        mask_soilc::AbstractVector{Bool}, bounds_col::UnitRange{Int},
        nlevdecomp::Int, i_litr_min::Int, i_litr_max::Int, i_cwd::Int,
        dt::Real, transient_landcover::Bool = false)
    fill!(matrix_Cinput, zero(eltype(matrix_Cinput)))
    fill!(matrix_Ninput, zero(eltype(matrix_Ninput)))
    cf(s) = getfield(cf_veg, s); nf(s) = getfield(nf_veg, s)
    out = _SoilInOut(; cin=matrix_Cinput, nin=matrix_Ninput)
    in_ = _SoilInIn(;
        lc1=cf(_SOILC_LITR_IN[1]), lc2=cf(_SOILC_LITR_IN[2]), lc3=cf(_SOILC_LITR_IN[3]),
        lc4=cf(_SOILC_LITR_IN_TR[1]), lc5=cf(_SOILC_LITR_IN_TR[2]),
        ln1=nf(_SOILN_LITR_IN[1]), ln2=nf(_SOILN_LITR_IN[2]), ln3=nf(_SOILN_LITR_IN[3]),
        ln4=nf(_SOILN_LITR_IN_TR[1]), ln5=nf(_SOILN_LITR_IN_TR[2]),
        wc1=cf(_SOILC_CWD_IN[1]), wc2=cf(_SOILC_CWD_IN[2]), wc3=cf(_SOILC_CWD_IN_TR[1]), wc4=cf(_SOILC_CWD_IN_TR[2]),
        wn1=nf(_SOILN_CWD_IN[1]), wn2=nf(_SOILN_CWD_IN[2]), wn3=nf(_SOILN_CWD_IN_TR[1]), wn4=nf(_SOILN_CWD_IN_TR[2]))
    backend = _kernel_backend(matrix_Cinput)
    _soil_input_kernel!(backend)(out, in_, mask_soilc, nlevdecomp, i_litr_min, i_litr_max, i_cwd,
        transient_landcover, eltype(matrix_Cinput)(Float64(dt)), first(bounds_col), last(bounds_col);
        ndrange = last(bounds_col))
    KA.synchronize(backend)
    return nothing
end

# Carbon-only isotope B-input: same litter/CWD structure as the C accumulator, reading
# the isotope veg carbon flux (phenology/gap/fire → litr/cwd). Fills matrix_Ciso_input
# (nc, ndecomp_pools_vr, local row ui). Used for the C13/C14 soil-matrix advance.
Base.@kwdef struct _SoilIsoIn{M3,M2}
    lc1::M3; lc2::M3; lc3::M3; lc4::M3; lc5::M3
    wc1::M2; wc2::M2; wc3::M2; wc4::M2
end
Adapt.@adapt_structure _SoilIsoIn

@kernel function _soil_iso_input_kernel!(ciso, in::_SoilIsoIn, @Const(mask), nlev::Int,
        ilmin::Int, ilmax::Int, icwd::Int, transient::Bool, this_begc::Int, dt, cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && mask[c]
        ui = c - this_begc + 1
        for j in 1:nlev
            for i in ilmin:ilmax
                vr = j + (i - 1) * nlev
                cc = in.lc1[c, j, i] + in.lc2[c, j, i] + in.lc3[c, j, i]
                if transient; cc += in.lc4[c, j, i] + in.lc5[c, j, i]; end
                ciso[ui, vr] = ciso[ui, vr] + cc * dt
            end
            vrc = j + (icwd - 1) * nlev
            cc = in.wc1[c, j] + in.wc2[c, j]
            if transient; cc += in.wc3[c, j] + in.wc4[c, j]; end
            ciso[ui, vrc] = ciso[ui, vrc] + cc * dt
        end
    end
end

function cn_soil_matrix_iso_input!(matrix_Ciso_input::AbstractMatrix{<:Real}, cf_iso;
        mask_soilc::AbstractVector{Bool}, bounds_col::UnitRange{Int},
        nlevdecomp::Int, i_litr_min::Int, i_litr_max::Int, i_cwd::Int,
        dt::Real, transient_landcover::Bool = false, begc::Int = first(bounds_col))
    fill!(matrix_Ciso_input, zero(eltype(matrix_Ciso_input)))
    cf(s) = getfield(cf_iso, s)
    in_ = _SoilIsoIn(;
        lc1=cf(_SOILC_LITR_IN[1]), lc2=cf(_SOILC_LITR_IN[2]), lc3=cf(_SOILC_LITR_IN[3]),
        lc4=cf(_SOILC_LITR_IN_TR[1]), lc5=cf(_SOILC_LITR_IN_TR[2]),
        wc1=cf(_SOILC_CWD_IN[1]), wc2=cf(_SOILC_CWD_IN[2]), wc3=cf(_SOILC_CWD_IN_TR[1]), wc4=cf(_SOILC_CWD_IN_TR[2]))
    backend = _kernel_backend(matrix_Ciso_input)
    _soil_iso_input_kernel!(backend)(matrix_Ciso_input, in_, mask_soilc, nlevdecomp, i_litr_min,
        i_litr_max, i_cwd, transient_landcover, begc, eltype(matrix_Ciso_input)(Float64(dt)),
        first(bounds_col), last(bounds_col); ndrange = last(bounds_col))
    KA.synchronize(backend)
    return nothing
end

# Fire-loss diagonal: matrix_decomp_fire_k = −m_decomp_cpools_to_fire / pool · dt (per column).
@kernel function _soil_firek_kernel!(fire_k, @Const(fflux), @Const(Xc), @Const(mask),
        nlev::Int, npool::Int, this_begc::Int, dt, cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && mask[c]
        ui = c - this_begc + 1
        for l in 1:npool, j in 1:nlev
            pool = Xc[c, j, l]
            fire_k[ui, j + (l - 1) * nlev] = pool > 0 ? -fflux[c, j, l] / pool * dt : zero(eltype(fire_k))
        end
    end
end

# Add a per-column persistent B-input field (indexed by global c) into the local Cin/Nin.
@kernel function _soil_col_add_kernel!(dst, @Const(src), @Const(mask), ndpvr::Int, this_begc::Int, cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && mask[c]
        ui = c - this_begc + 1
        for k in 1:ndpvr
            dst[ui, k] = dst[ui, k] + src[c, k]
        end
    end
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
        num_actfirec::Int = 0, transient_landcover::Bool = false,
        soilbgc_nf = nothing, soil_matrix_state = nothing,
        # Optional carbon isotopes: isotope soil-carbon state (decomp_cpools_vr_col) +
        # isotope veg-carbon flux (litter fields, for the isotope B-input). Advanced with
        # the same operator as C. nothing → isotope not carried.
        c13_soilbgc_cs = nothing, c13_cnveg_cf = nothing,
        c14_soilbgc_cs = nothing, c14_cnveg_cf = nothing,
        ref = nothing, FT::Type = Float64)

    begc = first(bounds_col); endc = last(bounds_col); nc = endc - begc + 1
    ndp_vr = ndecomp_pools * nlevdecomp
    nouttrans = count(==(0), cascade_con.cascade_receiver_pool)   # terminal (→atm) transitions

    # Static sparse index structure — build once (guard on the sentinel).
    if cascade_con.n_all_entries == SMM_EMPTY_INT
        init_soil_transfer!(cascade_con; ndecomp_pools=ndecomp_pools, nlevdecomp=nlevdecomp,
            ndecomp_cascade_transitions=ndecomp_cascade_transitions,
            ndecomp_cascade_outtransitions=nouttrans)
    end

    # B-input: the cn_driver-local litter/CWD fluxes (phenology/gap/fire, +harvest/gru
    # when transient) PLUS the persistent matrix_Cinput_col accumulated by the
    # dyn_subgrid driver (the dwt landcover-change inputs). The persistent field is
    # zeroed after the solve so it starts fresh next step.
    Cin = _smm_alloc(ref, FT, 0.0, nc, ndp_vr); Nin = _smm_alloc(ref, FT, 0.0, nc, ndp_vr)
    cn_soil_matrix_input_accumulate!(Cin, Nin, cnveg_cf, cnveg_nf;
        mask_soilc=mask_soilc, bounds_col=bounds_col, nlevdecomp=nlevdecomp,
        i_litr_min=i_litr_min, i_litr_max=i_litr_max, i_cwd=i_cwd, dt=dt,
        transient_landcover=transient_landcover)
    if size(soilbgc_cf.matrix_Cinput_col, 2) == ndp_vr
        _launch!(_soil_col_add_kernel!, Cin, soilbgc_cf.matrix_Cinput_col, mask_soilc, ndp_vr, begc, begc, endc; ndrange = endc)
    end
    if soilbgc_nf !== nothing && size(soilbgc_nf.matrix_Ninput_col, 2) == ndp_vr
        _launch!(_soil_col_add_kernel!, Nin, soilbgc_nf.matrix_Ninput_col, mask_soilc, ndp_vr, begc, begc, endc; ndrange = endc)
    end

    # Fire-loss diagonal (AKfiresoil). The fire module computed the decomp→fire flux
    # m_decomp_cpools_to_fire_vr (= pool·f·cmb); the matrix representation is the loss
    # fraction/step matrix_decomp_fire_k = −flux/pool·dt (port of CNFireBase:1239/1251,
    # derived here so the fire kernels stay untouched). Only assembled when a column is
    # actually burning; the pools are still at start-of-step so /pool is the fire base.
    fire_k = nothing
    if num_actfirec > 0
        fire_k = _smm_alloc(ref, FT, 0.0, nc, ndp_vr)
        _launch!(_soil_firek_kernel!, fire_k, cnveg_cf.m_decomp_cpools_to_fire_vr_col,
            soilbgc_cs.decomp_cpools_vr_col, mask_soilc, nlevdecomp, ndecomp_pools, begc,
            eltype(fire_k)(Float64(dt)), begc, endc; ndrange = endc)
    end

    # Reuse a caller-owned CNSoilMatrixState across steps to keep the memoized sparse
    # index structure (init_ready*/list_ready*); a fresh one per call is correct but
    # rebuilds the indices every step. Reuse is validated in test_cn_soil_matrix
    # ("memoized index reuse"). Falls back to a fresh state when none is supplied.
    # Isotope pools + B-inputs (carbon only). Built from the isotope veg-carbon flux.
    c13pools = c13_soilbgc_cs === nothing ? nothing : c13_soilbgc_cs.decomp_cpools_vr_col
    c14pools = c14_soilbgc_cs === nothing ? nothing : c14_soilbgc_cs.decomp_cpools_vr_col
    C13in = nothing; C14in = nothing
    if c13pools !== nothing && c13_cnveg_cf !== nothing
        C13in = zeros(nc, ndp_vr)
        cn_soil_matrix_iso_input!(C13in, c13_cnveg_cf; mask_soilc=mask_soilc,
            bounds_col=bounds_col, nlevdecomp=nlevdecomp, i_litr_min=i_litr_min,
            i_litr_max=i_litr_max, i_cwd=i_cwd, dt=dt,
            transient_landcover=transient_landcover, begc=begc)
    end
    if c14pools !== nothing && c14_cnveg_cf !== nothing
        C14in = zeros(nc, ndp_vr)
        cn_soil_matrix_iso_input!(C14in, c14_cnveg_cf; mask_soilc=mask_soilc,
            bounds_col=bounds_col, nlevdecomp=nlevdecomp, i_litr_min=i_litr_min,
            i_litr_max=i_litr_max, i_cwd=i_cwd, dt=dt,
            transient_landcover=transient_landcover, begc=begc)
    end

    ms = soil_matrix_state === nothing ? CNSoilMatrixState() : soil_matrix_state
    cn_soil_matrix!(ms, cascade_con;
        decomp_cpools_vr=soilbgc_cs.decomp_cpools_vr_col,
        decomp_npools_vr=soilbgc_ns.decomp_npools_vr_col,
        Ksoil=Ksoil, tri_ma_vr=soilbgc_cf.tri_ma_vr,
        matrix_Cinput=Cin, matrix_Ninput=Nin,
        rf_decomp_cascade=soilbgc_cf.rf_decomp_cascade_col,
        pathfrac_decomp_cascade=soilbgc_cf.pathfrac_decomp_cascade_col,
        matrix_decomp_fire_k=fire_k,
        decomp_c13pools_vr=c13pools, matrix_C13input=C13in,
        decomp_c14pools_vr=c14pools, matrix_C14input=C14in,
        mask_soilc=mask_soilc, begc=begc, endc=endc,
        nlevdecomp=nlevdecomp, ndecomp_pools=ndecomp_pools,
        ndecomp_cascade_transitions=ndecomp_cascade_transitions,
        ndecomp_cascade_outtransitions=nouttrans, num_actfirec=num_actfirec,
        ref=ref, FT=FT)

    # Reset the persistent B-input for the next step (its dwt contribution was consumed).
    size(soilbgc_cf.matrix_Cinput_col, 2) == ndp_vr && fill!(soilbgc_cf.matrix_Cinput_col, 0.0)
    soilbgc_nf !== nothing && size(soilbgc_nf.matrix_Ninput_col, 2) == ndp_vr &&
        fill!(soilbgc_nf.matrix_Ninput_col, 0.0)
    return nothing
end

# --------------------------------------------------------------------------
# GPU glue-loop kernels for cn_soil_matrix! — one thread per active soil column
# (fu → c = filter_c[fu], local row ui = c - begc + 1). Source pool arrays are
# indexed by the GLOBAL column c; per-unit workspaces by ui. Byte-identical to the
# host `for c in filter_soilc` loops on the CPU backend; device-capable when the
# decomp arrays + cc index bundle + workspaces are on the GPU.
# --------------------------------------------------------------------------

# C:N ratio of each pool: floating → C/N (or initial_cn if N==0); fixed → initial_cn.
@kernel function _soil_cndecomp_kernel!(cnd, @Const(dc), @Const(dn), @Const(floating),
        @Const(icn), @Const(filter_c), nlev::Int, npool::Int, this_begc::Int)
    fu = @index(Global)
    @inbounds begin
        c = filter_c[fu]; ui = c - this_begc + 1
        for l in 1:npool
            if floating[l]
                for j in 1:nlev
                    cnd[ui, j, l] = dn[c, j, l] > 0 ? dc[c, j, l] / dn[c, j, l] : icn[l]
                end
            else
                for j in 1:nlev
                    cnd[ui, j, l] = icn[l]
                end
            end
        end
    end
end

# Off-diagonal transfer coefficients a_ma_vr (C) + na_ma_vr (N, scaled by C:N ratio).
@kernel function _soil_ama_kernel!(a_ma, na_ma, @Const(cnd), @Const(rf), @Const(pathfrac),
        @Const(tranlist), @Const(donor), @Const(receiver), @Const(floating), @Const(filter_c),
        ndct::Int, nlev::Int, this_begc::Int)
    fu = @index(Global)
    @inbounds begin
        c = filter_c[fu]; ui = c - this_begc + 1
        for k in 1:ndct
            rp = receiver[k]
            if rp != 0
                dp = donor[k]
                for j in 1:nlev
                    ta = tranlist[j, k]
                    onem = one(eltype(a_ma)) - rf[c, j, k]
                    pf = pathfrac[c, j, k]
                    a_ma[ui, ta] = onem * pf
                    if !floating[rp]
                        na_ma[ui, ta] = onem * (cnd[ui, j, dp] / cnd[ui, j, rp]) * pf
                    else
                        na_ma[ui, ta] = pf
                    end
                end
            end
        end
    end
end

# Load C/N pools into the interim vectors + save the old-state snapshots.
@kernel function _soil_load_kernel!(Cint, Nint, Cold, Nold, @Const(dc), @Const(dn),
        @Const(filter_c), npool::Int, nlev::Int, this_begc::Int)
    fu = @index(Global)
    @inbounds begin
        c = filter_c[fu]; ui = c - this_begc + 1
        for i in 1:npool, j in 1:nlev
            idx = j + (i - 1) * nlev
            Cint[ui, idx] = dc[c, j, i]; Cold[ui, idx] = dc[c, j, i]
            Nint[ui, idx] = dn[c, j, i]; Nold[ui, idx] = dn[c, j, i]
        end
    end
end

# Ksoiln = Ksoil scaled by the fixed-C:N correction (only for non-floating pools, N>0).
# Ksoiln must be pre-initialized as a copy of Ksoil (the kernel only overwrites those).
@kernel function _soil_ksoiln_kernel!(Ksoiln, @Const(Ksoil), @Const(dc), @Const(dn),
        @Const(floating), @Const(icn), @Const(filter_c), npool::Int, nlev::Int, this_begc::Int)
    fu = @index(Global)
    @inbounds begin
        c = filter_c[fu]; ui = c - this_begc + 1
        for i in 1:npool
            if !floating[i]
                for j in 1:nlev
                    idx = j + (i - 1) * nlev
                    if dn[c, j, i] > 0
                        Ksoiln[ui, idx] = Ksoil[ui, idx] * dc[c, j, i] / dn[c, j, i] / icn[i]
                    end
                end
            end
        end
    end
end

# X = X + I·dt  (add the per-step C/N input vectors).
@kernel function _soil_binput_kernel!(Cint, Nint, @Const(Cin), @Const(Nin), @Const(filter_c),
        ndpvr::Int, this_begc::Int)
    fu = @index(Global)
    @inbounds begin
        c = filter_c[fu]; ui = c - this_begc + 1
        for j in 1:ndpvr
            Cint[ui, j] = Cin[ui, j] + Cint[ui, j]
            Nint[ui, j] = Nin[ui, j] + Nint[ui, j]
        end
    end
end

# Send advanced pools back to decomp_cpools_vr / decomp_npools_vr.
@kernel function _soil_writeback_kernel!(dc, dn, @Const(Cint), @Const(Nint), @Const(filter_c),
        npool::Int, nlev::Int, this_begc::Int)
    fu = @index(Global)
    @inbounds begin
        c = filter_c[fu]; ui = c - this_begc + 1
        for i in 1:npool, j in 1:nlev
            idx = j + (i - 1) * nlev
            dc[c, j, i] = Cint[ui, idx]; dn[c, j, i] = Nint[ui, idx]
        end
    end
end

# Isotope soil advance: load one decomp-pool array → vector, then (after the caller's
# spmm_ax! + B-input) write it back. Two small kernels reused for C13/C14.
@kernel function _soil_iso_load_kernel!(Xv, @Const(diso), @Const(filter_c), npool::Int, nlev::Int, this_begc::Int)
    fu = @index(Global)
    @inbounds begin
        c = filter_c[fu]; ui = c - this_begc + 1
        for i in 1:npool, j in 1:nlev
            Xv[ui, j + (i - 1) * nlev] = diso[c, j, i]
        end
    end
end
@kernel function _soil_iso_binput_kernel!(Xv, @Const(Iin), @Const(filter_c), ndpvr::Int, this_begc::Int)
    fu = @index(Global)
    @inbounds begin
        c = filter_c[fu]; ui = c - this_begc + 1
        for jj in 1:ndpvr
            Xv[ui, jj] = Iin[ui, jj] + Xv[ui, jj]
        end
    end
end
@kernel function _soil_iso_writeback_kernel!(diso, @Const(Xv), @Const(filter_c), npool::Int, nlev::Int, this_begc::Int)
    fu = @index(Global)
    @inbounds begin
        c = filter_c[fu]; ui = c - this_begc + 1
        for i in 1:npool, j in 1:nlev
            diso[c, j, i] = Xv[ui, j + (i - 1) * nlev]
        end
    end
end

# SASU/AKX spinup: accumulate the per-step input + reset the transfer accumulators.
@kernel function _soil_akx_inacc_kernel!(in_acc, in_nacc, @Const(Cin), @Const(Nin),
        @Const(filter_c), ndpvr::Int, this_begc::Int)
    fu = @index(Global)
    @inbounds begin
        c = filter_c[fu]; ui = c - this_begc + 1
        for j in 1:ndpvr
            in_acc[ui, j]  = in_acc[ui, j]  + Cin[ui, j]
            in_nacc[ui, j] = in_nacc[ui, j] + Nin[ui, j]
        end
    end
end
@kernel function _soil_akx_reset_kernel!(Mc, Mn, @Const(filter_c), NE::Int, this_begc::Int)
    fu = @index(Global)
    @inbounds begin
        c = filter_c[fu]; ui = c - this_begc + 1
        for jj in 1:NE
            Mc[ui, jj] = zero(eltype(Mc)); Mn[ui, jj] = zero(eltype(Mn))
        end
    end
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
        num_actfirec::Int=0,
        # Optional carbon-isotope decomp pools (C13/C14). When supplied they are advanced
        # with the SAME transfer operator AKallsoilc as C (transfer/vertical/fire fractions
        # are isotope-independent; C14 radioactive decay is applied separately by
        # c14_decay!), each with its own per-step B-input. Port of the matrix_Cinter13/14
        # SPMM_AX(AKallsoilc) path in CNSoilMatrixMod. Gated: nothing → no-op.
        decomp_c13pools_vr::Union{AbstractArray{<:Real,3},Nothing}=nothing,
        matrix_C13input::Union{AbstractMatrix{<:Real},Nothing}=nothing,
        decomp_c14pools_vr::Union{AbstractArray{<:Real,3},Nothing}=nothing,
        matrix_C14input::Union{AbstractMatrix{<:Real},Nothing}=nothing,
        ref=nothing, FT::Type=Float64)

    ndecomp_pools_vr = ndecomp_pools * nlevdecomp
    Ntrans = (ndecomp_cascade_transitions - ndecomp_cascade_outtransitions) * nlevdecomp
    Ntri_setup = cc.Ntri_setup

    # Active-unit filter; moved onto the device backend for the kernels (host otherwise).
    mh = _mask_to_host(mask_soilc); filter_host = [c for c in begc:endc if mh[c]]
    num_soilc = length(filter_host)
    filter_soilc = _backend_vec(ref, filter_host)

    cn_soil_matrix_alloc!(ms; ndecomp_pools_vr=ndecomp_pools_vr, begc=begc, endc=endc, ref=ref, FT=FT)

    # cc index + memoized-structure arrays used by the glue kernels + the memoized sparse
    # fills. On a device solve (ref given) they must be on the device — static after the
    # host warm-up, so device copies are correct. On the host path `_iv`/`_imat`/`_fvv` alias
    # cc directly (ref===nothing), so the first-call SetValueA/SPMP writes land on cc.
    _iv(v)   = _backend_vec(ref, v)                                   # Int/Bool vector → device
    _imat(m) = ref === nothing ? m : typeof(ref).name.wrapper(m)      # Int matrix → device
    _fvv(v)  = ref === nothing ? v : typeof(ref).name.wrapper(FT.(v)) # Float vector → device FT
    cascade_donor_pool    = _iv(cc.cascade_donor_pool)
    cascade_receiver_pool = _iv(cc.cascade_receiver_pool)
    floating              = _iv(collect(cc.floating_cn_ratio_decomp_pools))
    initial_cn_ratio      = _fvv(cc.initial_cn_ratio)
    spm_tranlist_a = _imat(cc.spm_tranlist_a)
    A_i = _iv(cc.A_i); A_j = _iv(cc.A_j); tri_i = _iv(cc.tri_i); tri_j = _iv(cc.tri_j)
    list_Asoilc = _iv(cc.list_Asoilc); RI_a = _iv(cc.RI_a); CI_a = _iv(cc.CI_a)
    list_Asoiln = _iv(cc.list_Asoiln); RI_na = _iv(cc.RI_na); CI_na = _iv(cc.CI_na)
    list_AK_AKV = _iv(cc.list_AK_AKV); list_V_AKV = _iv(cc.list_V_AKV)
    list_AK_AKVfire = _iv(cc.list_AK_AKVfire); list_V_AKVfire = _iv(cc.list_V_AKVfire)
    list_fire_AKVfire = _iv(cc.list_fire_AKVfire)
    RI_AKallc = _iv(cc.RI_AKallsoilc); CI_AKallc = _iv(cc.CI_AKallsoilc)
    RI_AKalln = _iv(cc.RI_AKallsoiln); CI_AKalln = _iv(cc.CI_AKallsoiln)

    nunit = endc - begc + 1

    # ----- C:N ratios of applicable pools (gC/gN) -----
    cn_decomp_pools = _smm_alloc(ref, FT, 0.0, nunit, nlevdecomp, ndecomp_pools)
    _launch!(_soil_cndecomp_kernel!, cn_decomp_pools, decomp_cpools_vr, decomp_npools_vr,
             floating, initial_cn_ratio, filter_soilc, nlevdecomp, ndecomp_pools, begc; ndrange = num_soilc)

    # ----- Assemble a_ma_vr / na_ma_vr (off-diagonal transfer coefficients). -----
    a_ma_vr  = _smm_alloc(ref, FT, 0.0, nunit, Ntrans)
    na_ma_vr = _smm_alloc(ref, FT, 0.0, nunit, Ntrans)
    _launch!(_soil_ama_kernel!, a_ma_vr, na_ma_vr, cn_decomp_pools, rf_decomp_cascade,
             pathfrac_decomp_cascade, spm_tranlist_a, cascade_donor_pool, cascade_receiver_pool,
             floating, filter_soilc, ndecomp_cascade_transitions, nlevdecomp, begc; ndrange = num_soilc)

    # ----- Set old vector value (X(t)) + save old-state snapshots. -----
    Cinter_old = _smm_alloc(ref, FT, 0.0, nunit, ndecomp_pools_vr)
    Ninter_old = _smm_alloc(ref, FT, 0.0, nunit, ndecomp_pools_vr)
    _launch!(_soil_load_kernel!, ms.matrix_Cinter.V, ms.matrix_Ninter.V, Cinter_old, Ninter_old,
             decomp_cpools_vr, decomp_npools_vr, filter_soilc, ndecomp_pools, nlevdecomp, begc; ndrange = num_soilc)

    # ----- Build Ksoiln (= Ksoil scaled by the fixed-C:N correction). -----
    Ksoiln = copy(Ksoil)
    _launch!(_soil_ksoiln_kernel!, Ksoiln, Ksoil, decomp_cpools_vr, decomp_npools_vr,
             floating, initial_cn_ratio, filter_soilc, ndecomp_pools, nlevdecomp, begc; ndrange = num_soilc)

    # ----- Set C transfer matrix Ac from a_ma_vr (adds the -1 diagonal). -----
    ms.init_readyAsoilc = set_value_a!(ms.AKsoilc, begc, endc, num_soilc, filter_soilc,
        a_ma_vr, A_i, A_j, Ntrans, ms.init_readyAsoilc;
        list=list_Asoilc, RI_A=RI_a, CI_A=CI_a)
    ms.init_readyAsoiln = set_value_a!(ms.AKsoiln, begc, endc, num_soilc, filter_soilc,
        na_ma_vr, A_i, A_j, Ntrans, ms.init_readyAsoiln;
        list=list_Asoiln, RI_A=RI_na, CI_A=CI_na)

    # ----- Ac*K, An*K. -----
    Kdm  = DiagMatrixType();  init_dm!(Kdm,  ndecomp_pools_vr, begc, endc; ref=ref, FT=FT)
    Kndm = DiagMatrixType();  init_dm!(Kndm, ndecomp_pools_vr, begc, endc; ref=ref, FT=FT)
    set_value_dm!(Kdm,  begc, endc, num_soilc, filter_soilc, Ksoil)
    set_value_dm!(Kndm, begc, endc, num_soilc, filter_soilc, Ksoiln)
    spmm_ak!(ms.AKsoilc, num_soilc, filter_soilc, Kdm)
    spmm_ak!(ms.AKsoiln, num_soilc, filter_soilc, Kndm)

    # ----- Set vertical-transport matrix V from tri_ma_vr. -----
    set_value_sm!(ms.AVsoil, begc, endc, num_soilc, filter_soilc,
        tri_ma_vr, tri_i, tri_j, Ntri_setup)

    # ----- Set fire matrix Kfire from matrix_decomp_fire_k. -----
    if num_actfirec != 0
        @assert matrix_decomp_fire_k !== nothing "cn_soil_matrix!: fire active but matrix_decomp_fire_k not supplied"
        kfire_i = _iv(collect(1:ndecomp_pools_vr))
        kfire_j = _iv(collect(1:ndecomp_pools_vr))
        set_value_sm!(ms.AKfiresoil, begc, endc, num_soilc, filter_soilc,
            matrix_decomp_fire_k, kfire_i, kfire_j, ndecomp_pools_vr)
    end

    # ----- AKallsoilc = A*K + V (− Kfire). The merged structure (NE/RI/CI) is
    # memorized on `cc` after the first call and reused (list_ready=true). -----
    if num_actfirec == 0
        (ms.list_ready1_nofire, cc.NE_AKallsoilc) = spmp_ab!(ms.AKallsoilc, num_soilc, filter_soilc,
            ms.AKsoilc, ms.AVsoil, ms.list_ready1_nofire;
            list_A=list_AK_AKV, list_B=list_V_AKV,
            NE_AB=cc.NE_AKallsoilc, RI_AB=RI_AKallc, CI_AB=CI_AKallc)
        (ms.list_ready2_nofire, cc.NE_AKallsoiln) = spmp_ab!(ms.AKallsoiln, num_soilc, filter_soilc,
            ms.AKsoiln, ms.AVsoil, ms.list_ready2_nofire;
            list_A=list_AK_AKV, list_B=list_V_AKV,
            NE_AB=cc.NE_AKallsoiln, RI_AB=RI_AKalln, CI_AB=CI_AKalln)
    else
        filter_actfirec = filter_soilc  # caller scopes fire via matrix_decomp_fire_k zeros
        (ms.list_ready1_fire, cc.NE_AKallsoilc) = spmp_abc!(ms.AKallsoilc, num_soilc, filter_soilc,
            ms.AKsoilc, ms.AVsoil, ms.AKfiresoil, ms.list_ready1_fire;
            list_A=list_AK_AKVfire, list_B=list_V_AKVfire, list_C=list_fire_AKVfire,
            NE_ABC=cc.NE_AKallsoilc, RI_ABC=RI_AKallc, CI_ABC=CI_AKallc,
            num_actunit_C=num_actfirec, filter_actunit_C=filter_actfirec)
        (ms.list_ready2_fire, cc.NE_AKallsoiln) = spmp_abc!(ms.AKallsoiln, num_soilc, filter_soilc,
            ms.AKsoiln, ms.AVsoil, ms.AKfiresoil, ms.list_ready2_fire;
            list_A=list_AK_AKVfire, list_B=list_V_AKVfire, list_C=list_fire_AKVfire,
            NE_ABC=cc.NE_AKallsoiln, RI_ABC=RI_AKalln, CI_ABC=CI_AKalln,
            num_actunit_C=num_actfirec, filter_actunit_C=filter_actfirec)
    end

    # ----- Advance pools: X = X + (A*K + V − Kfire)·X. -----
    spmm_ax!(ms.matrix_Cinter, num_soilc, filter_soilc, ms.AKallsoilc)
    spmm_ax!(ms.matrix_Ninter, num_soilc, filter_soilc, ms.AKallsoiln)

    # ----- Add input: X = X + I·dt. -----
    _launch!(_soil_binput_kernel!, ms.matrix_Cinter.V, ms.matrix_Ninter.V, matrix_Cinput,
             matrix_Ninput, filter_soilc, ndecomp_pools_vr, begc; ndrange = num_soilc)

    # ----- Send back to decomp_cpools_vr / decomp_npools_vr. -----
    _launch!(_soil_writeback_kernel!, decomp_cpools_vr, decomp_npools_vr, ms.matrix_Cinter.V,
             ms.matrix_Ninter.V, filter_soilc, ndecomp_pools, nlevdecomp, begc; ndrange = num_soilc)

    # ----- Carbon isotopes: advance with the SAME AKallsoilc as C, own B-input. -----
    _advance_iso_soil!(decomp_c13pools_vr, matrix_C13input, ms.AKallsoilc,
        ndecomp_pools_vr, ndecomp_pools, nlevdecomp, begc, endc, filter_soilc, num_soilc; ref=ref, FT=FT)
    _advance_iso_soil!(decomp_c14pools_vr, matrix_C14input, ms.AKallsoilc,
        ndecomp_pools_vr, ndecomp_pools, nlevdecomp, begc, endc, filter_soilc, num_soilc; ref=ref, FT=FT)

    release_dm!(Kdm); release_dm!(Kndm)
    return (Cinter_old, Ninter_old)
end

# Advance one carbon-isotope decomp-pool array with the (already-built) C transfer
# operator AKallsoilc: X_iso = X_iso + AKallsoilc·X_iso + I_iso·dt. No-op when the pool
# array is nothing (isotope disabled). A missing B-input is treated as zero.
function _advance_iso_soil!(decomp_iso_vr, matrix_Ciso_input, AKallsoilc,
        ndecomp_pools_vr::Int, ndecomp_pools::Int, nlevdecomp::Int,
        begc::Int, endc::Int, filter_soilc::AbstractVector{<:Integer}, num_soilc::Int;
        ref=nothing, FT::Type=Float64)
    decomp_iso_vr === nothing && return nothing
    Xiso = VectorType(); init_v!(Xiso, ndecomp_pools_vr, begc, endc; ref=ref, FT=FT)
    _launch!(_soil_iso_load_kernel!, Xiso.V, decomp_iso_vr, filter_soilc, ndecomp_pools, nlevdecomp, begc; ndrange = num_soilc)
    spmm_ax!(Xiso, num_soilc, filter_soilc, AKallsoilc)
    if matrix_Ciso_input !== nothing
        _launch!(_soil_iso_binput_kernel!, Xiso.V, matrix_Ciso_input, filter_soilc, ndecomp_pools_vr, begc; ndrange = num_soilc)
    end
    _launch!(_soil_iso_writeback_kernel!, decomp_iso_vr, Xiso.V, filter_soilc, ndecomp_pools, nlevdecomp, begc; ndrange = num_soilc)
    return nothing
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
        first_step::Bool=false, compute_capacity::Bool=false,
        ref=nothing, FT::Type=Float64)

    ndecomp_pools_vr = ndecomp_pools * nlevdecomp
    mh = _mask_to_host(mask_soilc); filter_host = [c for c in begc:endc if mh[c]]
    num_soilc = length(filter_host)
    filter_soilc = _backend_vec(ref, filter_host)
    epsi = 1.0e-8

    # Accumulate input.
    _launch!(_soil_akx_inacc_kernel!, in_acc, in_nacc, matrix_Cinput, matrix_Ninput,
             filter_soilc, ndecomp_pools_vr, begc; ndrange = num_soilc)

    # AKallsoilc/n → actual transfers (× X_old).
    set_value_dm!(ms.Xdiagsoil, begc, endc, num_soilc, filter_soilc, Cinter_old)
    spmm_ak!(ms.AKallsoilc, num_soilc, filter_soilc, ms.Xdiagsoil)
    set_value_dm!(ms.Xdiagsoil, begc, endc, num_soilc, filter_soilc, Ninter_old)
    spmm_ak!(ms.AKallsoiln, num_soilc, filter_soilc, ms.Xdiagsoil)

    # Accumulate transfers AKXcacc += AKallsoilc ; AKXnacc += AKallsoiln.
    if first_step || !is_alloc_sm(ms.AKXcacc)
        is_alloc_sm(ms.AKXcacc) || init_sm!(ms.AKXcacc, ndecomp_pools_vr, begc, endc; ref=ref, FT=FT)
        is_alloc_sm(ms.AKXnacc) || init_sm!(ms.AKXnacc, ndecomp_pools_vr, begc, endc; ref=ref, FT=FT)
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

    # ----- End-of-SASU-year: compute capacity X* = −A^{-1} I. This is a per-column DENSE
    # matrix inverse — GPU-hostile and run once per SASU year — so it is done on the host:
    # the (device or host) accumulators are gathered with `Array` first (a no-op copy on the
    # host path, so byte-identical). -----
    n_all_entries = cc.n_all_entries
    all_i = cc.all_i
    all_j = cc.all_j
    Mc = Array(ms.AKXcacc.M); Mn = Array(ms.AKXnacc.M)
    iac = Array(in_acc); ina = Array(in_nacc)
    d0c = Array(decomp0_cpools_vr); d0n = Array(decomp0_npools_vr)
    soilmatrixc_cap = zeros(Float64, num_soilc, ndecomp_pools_vr)
    soilmatrixn_cap = zeros(Float64, num_soilc, ndecomp_pools_vr)

    for (uf, c) in enumerate(filter_host)
        ui = c - begc + 1
        tran_acc  = zeros(Float64, ndecomp_pools_vr, ndecomp_pools_vr)
        tran_nacc = zeros(Float64, ndecomp_pools_vr, ndecomp_pools_vr)
        for jj in 1:n_all_entries
            j_lev    = mod(all_j[jj] - 1, nlevdecomp) + 1
            j_decomp = div(all_j[jj] - j_lev, nlevdecomp) + 1
            tran_acc[all_i[jj], all_j[jj]]  = Float64(Mc[ui, jj]) / d0c[c, j_lev, j_decomp]
            tran_nacc[all_i[jj], all_j[jj]] = Float64(Mn[ui, jj]) / d0n[c, j_lev, j_decomp]
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
        soilmatrixc_cap[uf, :] = -(AKinv  * Float64.(iac[ui, 1:ndecomp_pools_vr]))
        soilmatrixn_cap[uf, :] = -(AKinvn * Float64.(ina[ui, 1:ndecomp_pools_vr]))
        for k in 1:ndecomp_pools_vr
            soilmatrixc_cap[uf, k] < 0 && (soilmatrixc_cap[uf, k] = 0.0)
            soilmatrixn_cap[uf, k] < 0 && (soilmatrixn_cap[uf, k] = 0.0)
        end
    end

    # Reset accumulators.
    _launch!(_soil_akx_reset_kernel!, ms.AKXcacc.M, ms.AKXnacc.M, filter_soilc, n_all_entries, begc; ndrange = num_soilc)
    fill!(in_acc, zero(eltype(in_acc)))
    fill!(in_nacc, zero(eltype(in_nacc)))

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
