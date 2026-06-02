# ==========================================================================
# Ported from: src/soilbiogeochem/SoilBiogeochemDecompMod.F90
# Module holding routines used in litter and soil decomposition model.
#
# Public functions:
#   decomp_read_params!          — Read decomposition parameters (dnp)
#   soil_biogeochem_decomp!      — Main decomposition routine
#
# GPU: soil_biogeochem_decomp! runs whole-function on Metal. The C:N-ratio
# loop, the decomposition-cascade application, the methane fphr loop, and the
# vertical integration are each per-column KernelAbstractions kernels with
# internal sequential l/j/k loops (the c_state_update1! in-thread cascade
# pattern — every write is own-index [c,j,k|l], no cross-thread scatter, so
# race-free and byte-identical on the Float64 CPU path). The static cascade
# config (pool/transition metadata) stays on the host; only the per-timestep
# compute is kernelized.
# ==========================================================================

# ---------------------------------------------------------------------------
# DecompParams — decomposition parameters
# Ported from params_type in SoilBiogeochemDecompMod.F90
# ---------------------------------------------------------------------------

Base.@kwdef mutable struct DecompParams
    dnp::Float64 = 0.01  # denitrification proportion
end

# Constant: atmosphere pool index (receiver_pool == 0 means "to atmosphere")
# Ported from i_atm in SoilBiogeochemDecompCascadeConType.F90
const I_ATM = 0

# ---------------------------------------------------------------------------
# decomp_read_params! — Read decomposition parameters
# Ported from readParams in SoilBiogeochemDecompMod.F90
# ---------------------------------------------------------------------------

function decomp_read_params!(params::DecompParams; dnp::Real)
    params.dnp = dnp
    return nothing
end

# ---------------------------------------------------------------------------
# Device-view bundles for soil_biogeochem_decomp!
#
# This routine touches ~30 column/cascade arrays — over the Metal ~31-arg
# kernel cap once masks + scalars are added — so the heavy arrays are gathered
# into isbits-on-device bundle structs (Adapt.@adapt_structure'd). Field names
# mirror the live-state variable paths so the kernel bodies read verbatim.
# Bundles are built per-launch from the live arrays, so every write flows
# straight back into the originating struct (shared array refs).
#
# Separate type params: V (1D vectors), M2 (2D col,level arrays), M3 (3D
# col,level,transition|pool arrays). Keeping them distinct lets the adaptor
# narrow each independently (a Float64-default work array can stay Float64 even
# when the state arrays downcast).
# ---------------------------------------------------------------------------

# 3D cascade/pool arrays (col, level, transition|pool).
Base.@kwdef struct _DecompCascadeDV{M3}
    decomp_cpools_vr::M3
    decomp_npools_vr::M3
    cn_decomp_pools::M3
    rf_decomp_cascade::M3
    c_overflow_vr::M3
    p_decomp_cpool_loss::M3
    pmnf_decomp_cascade::M3
    p_decomp_npool_to_din::M3
    decomp_cascade_hr_vr::M3
    decomp_cascade_ctransfer_vr::M3
    decomp_cascade_ntransfer_vr::M3
    decomp_cascade_sminn_flux_vr::M3
    sminn_to_denit_decomp_cascade_vr::M3
end
Adapt.@adapt_structure _DecompCascadeDV

# 2D (col, level) arrays.
Base.@kwdef struct _DecompColDV{M2}
    fpi_vr::M2
    w_scalar::M2
    phr_vr::M2
    fphr::M2
    net_nmin_vr::M2
    gross_nmin_vr::M2
end
Adapt.@adapt_structure _DecompColDV

# ---------------------------------------------------------------------------
# Host-vector -> device-backend helpers. Move a host config vector onto the
# backend of a state-array prototype, preserving its element-type class. On the
# CPU path the prototype is an Array so these return ordinary host vectors
# (byte-identical to the originals). On the GPU path they allocate a device
# vector of the matching eltype and copy in — NEVER float-coercing an index/flag.
# ---------------------------------------------------------------------------
@inline function _to_dev_int(proto, v::AbstractVector{<:Integer})
    d = similar(proto, eltype(v), length(v)); copyto!(d, v); return d
end
@inline function _to_dev_bool(proto, v::AbstractVector)
    d = similar(proto, Bool, length(v)); copyto!(d, collect(Bool, v)); return d
end
@inline function _to_dev_real(proto, v::AbstractVector{<:Real})
    d = similar(proto, eltype(proto), length(v)); copyto!(d, v); return d
end

# ---------------------------------------------------------------------------
# Kernel: C:N ratios of applicable pools (per column; internal l/j loops).
# Branch on the per-pool floating-C:N flag; own-index writes [c,j,l].
# ---------------------------------------------------------------------------
@kernel function _decomp_cnratio_kernel!(@Const(mask), dv,
        @Const(floating_cn_ratio_decomp_pools), @Const(initial_cn_ratio),
        ndecomp_pools::Int, nlevdecomp::Int)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(dv.cn_decomp_pools)
        for l in 1:ndecomp_pools
            if floating_cn_ratio_decomp_pools[l]
                for j in 1:nlevdecomp
                    if dv.decomp_npools_vr[c, j, l] > zero(T)
                        dv.cn_decomp_pools[c, j, l] =
                            dv.decomp_cpools_vr[c, j, l] / dv.decomp_npools_vr[c, j, l]
                    end
                end
            else
                for j in 1:nlevdecomp
                    dv.cn_decomp_pools[c, j, l] = initial_cn_ratio[l]
                end
            end
        end
    end
end

# ---------------------------------------------------------------------------
# Kernel: apply the decomposition cascade (per column; internal k/j loops).
# All writes are own-index [c,j,k] (+ own-index [c,j] reduction into
# net_nmin_vr). The fpi_vr limitation, denitrification, HR/transfer split, the
# MIMICS overflow correction, and the N-transfer / sminn-flux bookkeeping are
# carried verbatim. Config flags branch inside the kernel.
# ---------------------------------------------------------------------------
@kernel function _decomp_cascade_kernel!(@Const(mask), dv, cdv,
        @Const(cascade_donor_pool), @Const(cascade_receiver_pool),
        ndecomp_cascade_transitions::Int, nlevdecomp::Int, dnp,
        use_nitrif_denitrif::Bool, use_mimics::Bool)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(dv.decomp_cpools_vr)
        for k in 1:ndecomp_cascade_transitions
            donor = cascade_donor_pool[k]
            recv  = cascade_receiver_pool[k]
            for j in 1:nlevdecomp
                if dv.decomp_cpools_vr[c, j, donor] > zero(T)
                    if dv.pmnf_decomp_cascade[c, j, k] > zero(T)
                        dv.p_decomp_cpool_loss[c, j, k] *= cdv.fpi_vr[c, j]
                        dv.pmnf_decomp_cascade[c, j, k] *= cdv.fpi_vr[c, j]
                        # use_soil_matrixcn Ksoil%DM correction not yet implemented
                        if !use_nitrif_denitrif
                            dv.sminn_to_denit_decomp_cascade_vr[c, j, k] = zero(T)
                        end
                    else
                        if !use_nitrif_denitrif
                            dv.sminn_to_denit_decomp_cascade_vr[c, j, k] =
                                -dnp * dv.pmnf_decomp_cascade[c, j, k]
                        end
                    end

                    dv.decomp_cascade_hr_vr[c, j, k] =
                        dv.rf_decomp_cascade[c, j, k] * dv.p_decomp_cpool_loss[c, j, k]
                    dv.decomp_cascade_ctransfer_vr[c, j, k] =
                        (one(T) - dv.rf_decomp_cascade[c, j, k]) * dv.p_decomp_cpool_loss[c, j, k]

                    if use_mimics
                        dv.decomp_cascade_hr_vr[c, j, k] = min(
                            dv.p_decomp_cpool_loss[c, j, k],
                            dv.decomp_cascade_hr_vr[c, j, k] + dv.c_overflow_vr[c, j, k])
                        dv.decomp_cascade_ctransfer_vr[c, j, k] = max(zero(T),
                            dv.p_decomp_cpool_loss[c, j, k] - dv.decomp_cascade_hr_vr[c, j, k])
                    end

                    if dv.decomp_npools_vr[c, j, donor] > zero(T) && recv != I_ATM
                        dv.decomp_cascade_ntransfer_vr[c, j, k] =
                            dv.p_decomp_cpool_loss[c, j, k] / dv.cn_decomp_pools[c, j, donor]
                    else
                        dv.decomp_cascade_ntransfer_vr[c, j, k] = zero(T)
                    end

                    if recv != 0
                        dv.decomp_cascade_sminn_flux_vr[c, j, k] = dv.pmnf_decomp_cascade[c, j, k]
                    else  # keep sign convention negative for terminal pools
                        dv.decomp_cascade_sminn_flux_vr[c, j, k] = -dv.pmnf_decomp_cascade[c, j, k]
                    end

                    cdv.net_nmin_vr[c, j] -= dv.pmnf_decomp_cascade[c, j, k]

                    if use_mimics
                        dv.decomp_cascade_sminn_flux_vr[c, j, k] -= dv.p_decomp_npool_to_din[c, j, k]
                        cdv.net_nmin_vr[c, j] += dv.p_decomp_npool_to_din[c, j, k]
                    end
                else
                    dv.decomp_cascade_ntransfer_vr[c, j, k] = zero(T)
                    if !use_nitrif_denitrif
                        dv.sminn_to_denit_decomp_cascade_vr[c, j, k] = zero(T)
                    end
                    dv.decomp_cascade_sminn_flux_vr[c, j, k] = zero(T)
                end
            end
        end
    end
end

# ---------------------------------------------------------------------------
# Kernel: total fraction of potential HR for the methane code (per column).
# Internal k/j loop accumulates hrsum into own-index [c,j] (sequential, ascending
# k — byte-identical), then the per-(c,j) fphr is computed in place.
# ---------------------------------------------------------------------------
@kernel function _decomp_fphr_kernel!(@Const(mask), dv, cdv,
        ndecomp_cascade_transitions::Int, nlevdecomp::Int)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(cdv.fphr)
        for j in 1:nlevdecomp
            hrsum = zero(T)
            for k in 1:ndecomp_cascade_transitions
                hrsum += dv.rf_decomp_cascade[c, j, k] * dv.p_decomp_cpool_loss[c, j, k]
            end
            if cdv.phr_vr[c, j] > zero(T)
                f = hrsum / cdv.phr_vr[c, j] * cdv.w_scalar[c, j]
                cdv.fphr[c, j] = max(f, T(0.01))  # Prevent overflow for 0 respiration
            else
                cdv.fphr[c, j] = one(T)
            end
        end
    end
end

# ---------------------------------------------------------------------------
# Kernel: vertical integration of net/gross mineralization (per column).
# Own-column reduction over levels (ascending j — byte-identical) into [c].
# ---------------------------------------------------------------------------
@kernel function _decomp_vint_kernel!(net_nmin, @Const(mask), gross_nmin, cdv,
        @Const(dzsoi_decomp), nlevdecomp::Int)
    c = @index(Global)
    @inbounds if mask[c]
        for j in 1:nlevdecomp
            net_nmin[c]   += cdv.net_nmin_vr[c, j] * dzsoi_decomp[j]
            gross_nmin[c] += cdv.gross_nmin_vr[c, j] * dzsoi_decomp[j]
        end
    end
end

# ---------------------------------------------------------------------------
# soil_biogeochem_decomp! — Main decomposition subroutine
# Ported from SoilBiogeochemDecomp in SoilBiogeochemDecompMod.F90
#
# Calculates actual immobilization and decomp rates, following resolution
# of plant/heterotroph competition for mineral N.
# ---------------------------------------------------------------------------

function soil_biogeochem_decomp!(
        cf::SoilBiogeochemCarbonFluxData,
        cs::SoilBiogeochemCarbonStateData,
        nf::SoilBiogeochemNitrogenFluxData,
        ns::SoilBiogeochemNitrogenStateData,
        st::SoilBiogeochemStateData,
        cascade_con::DecompCascadeConData,
        params::DecompParams;
        mask_bgc_soilc::AbstractVector{Bool},
        bounds::UnitRange{Int},
        nlevdecomp::Int,
        ndecomp_pools::Int,
        ndecomp_cascade_transitions::Int,
        cn_decomp_pools::AbstractArray{<:Real,3},
        p_decomp_cpool_loss::AbstractArray{<:Real,3},
        pmnf_decomp_cascade::AbstractArray{<:Real,3},
        p_decomp_npool_to_din::AbstractArray{<:Real,3},
        dzsoi_decomp::AbstractVector{<:Real},
        use_nitrif_denitrif::Bool=false,
        use_lch4::Bool=false,
        use_mimics::Bool=false,
        use_soil_matrixcn::Bool=false)

    # Cascade configuration lives on the host (the cascade_con struct holds host
    # Int/Bool/String metadata). The per-transition index/flag/cn vectors the
    # kernels index must be device-resident on the GPU path, so copy them onto the
    # backend of a state array (`cn_decomp_pools`), PRESERVING the integer eltype
    # (a float-coerced index breaks device indexing) and turning the BitVector
    # floating-CN flag into a plain Bool device array. On the CPU path `_to_dev_*`
    # return live host copies — byte-identical.
    _proto = cn_decomp_pools                       # backend/precision prototype
    cascade_donor_pool    = _to_dev_int(_proto, cascade_con.cascade_donor_pool)
    cascade_receiver_pool = _to_dev_int(_proto, cascade_con.cascade_receiver_pool)
    floating_cn_ratio_decomp_pools =
        _to_dev_bool(_proto, cascade_con.floating_cn_ratio_decomp_pools)
    initial_cn_ratio      = _to_dev_real(_proto, cascade_con.initial_cn_ratio)

    # Device-view bundles (shared array refs → writes flow back into the structs).
    dv = _DecompCascadeDV(;
        decomp_cpools_vr                 = cs.decomp_cpools_vr_col,
        decomp_npools_vr                 = ns.decomp_npools_vr_col,
        cn_decomp_pools                  = cn_decomp_pools,
        rf_decomp_cascade                = cf.rf_decomp_cascade_col,
        c_overflow_vr                    = cf.c_overflow_vr,
        p_decomp_cpool_loss              = p_decomp_cpool_loss,
        pmnf_decomp_cascade              = pmnf_decomp_cascade,
        p_decomp_npool_to_din            = p_decomp_npool_to_din,
        decomp_cascade_hr_vr             = cf.decomp_cascade_hr_vr_col,
        decomp_cascade_ctransfer_vr      = cf.decomp_cascade_ctransfer_vr_col,
        decomp_cascade_ntransfer_vr      = nf.decomp_cascade_ntransfer_vr_col,
        decomp_cascade_sminn_flux_vr     = nf.decomp_cascade_sminn_flux_vr_col,
        sminn_to_denit_decomp_cascade_vr = nf.sminn_to_denit_decomp_cascade_vr_col)

    cdv = _DecompColDV(;
        fpi_vr        = st.fpi_vr_col,
        w_scalar      = cf.w_scalar_col,
        phr_vr        = cf.phr_vr_col,
        fphr          = cf.fphr_col,
        net_nmin_vr   = nf.net_nmin_vr_col,
        gross_nmin_vr = nf.gross_nmin_vr_col)

    # Scalar parameter at working precision (no Float64 reaches a Float32 backend).
    _T  = eltype(cf.rf_decomp_cascade_col)
    dnp = _T(params.dnp)

    # -------------------------------------------------------------------
    # Calculate c:n ratios of applicable pools
    # -------------------------------------------------------------------
    _launch!(_decomp_cnratio_kernel!, mask_bgc_soilc, dv,
             floating_cn_ratio_decomp_pools, initial_cn_ratio,
             ndecomp_pools, nlevdecomp)

    # -------------------------------------------------------------------
    # Calculate actual immobilization and decomp rates, following
    # resolution of plant/heterotroph competition for mineral N.
    # Only the immobilization steps are limited by fpi_vr (pmnf > 0).
    # Also calculate denitrification losses as a simple proportion of
    # mineralization flux.  (use_soil_matrixcn Ksoil%DM correction TODO)
    # -------------------------------------------------------------------
    _launch!(_decomp_cascade_kernel!, mask_bgc_soilc, dv, cdv,
             cascade_donor_pool, cascade_receiver_pool,
             ndecomp_cascade_transitions, nlevdecomp, dnp,
             use_nitrif_denitrif, use_mimics)

    # -------------------------------------------------------------------
    # Calculate total fraction of potential HR, for methane code
    # -------------------------------------------------------------------
    if use_lch4
        _launch!(_decomp_fphr_kernel!, mask_bgc_soilc, dv, cdv,
                 ndecomp_cascade_transitions, nlevdecomp)
    end

    # -------------------------------------------------------------------
    # Vertically integrate net and gross mineralization fluxes
    # for diagnostic output
    # -------------------------------------------------------------------
    _launch!(_decomp_vint_kernel!, nf.net_nmin_col, mask_bgc_soilc,
             nf.gross_nmin_col, cdv, dzsoi_decomp, nlevdecomp)

    return nothing
end
