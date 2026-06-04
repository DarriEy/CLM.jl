# ==========================================================================
# Ported from: src/biogeochem/CNMRespMod.F90
# Maintenance respiration fluxes for coupled carbon-nitrogen code
# ==========================================================================

# ---------------------------------------------------------------------------
# Module-level parameters (from params_type in CNMRespMod)
# ---------------------------------------------------------------------------

"""
    MaintRespParams

Parameters for maintenance respiration, read from namelist/file.
Corresponds to `params_type` in `CNMRespMod.F90`.
"""
Base.@kwdef mutable struct MaintRespParams
    br      ::Float64 = SPVAL   # base rate for maintenance respiration (gC/gN/s)
    br_root ::Float64 = SPVAL   # base rate for maintenance respiration for roots (gC/gN/s)
end

"""
    maint_resp_read_params!(params; br=2.525e-6, br_root=SPVAL)

Read/set parameters for maintenance respiration.
Combines `CNMRespReadNML` and `readParams` from `CNMRespMod.F90`.

If `br_root` is not provided (remains SPVAL), it defaults to `br`.
"""
function maint_resp_read_params!(params::MaintRespParams;
                                  br::Real = 2.525e-6,
                                  br_root::Real = SPVAL)
    params.br = br
    if br_root == SPVAL
        params.br_root = br
    else
        params.br_root = br_root
    end
    return nothing
end

# ---------------------------------------------------------------------------
# PFT constants needed by maintenance respiration (from pftcon)
# ---------------------------------------------------------------------------

"""
    PftConMaintResp

PFT-level parameters used by the maintenance respiration routine.
Contains a subset of fields from `pftconMod` referenced in `CNMRespMod.F90`.
"""
Base.@kwdef mutable struct PftConMaintResp{FT<:Real, V<:AbstractVector{FT}}
    woody ::V = Float64[]   # binary woody flag (1=woody, 0=not woody)
end
PftConMaintResp{FT}(; kwargs...) where {FT<:Real} = PftConMaintResp{FT, Vector{FT}}(; kwargs...)
Adapt.@adapt_structure PftConMaintResp

# ---------------------------------------------------------------------------
# cn_mresp! — Maintenance respiration
# ---------------------------------------------------------------------------

"""
    cn_mresp!(mask_soilc, mask_soilp, bounds_c, bounds_p,
              params, cn_params_share, pftcon,
              patch, canopystate, soilstate, temperature, photosyns,
              cnveg_cf, cnveg_ns;
              nlevgrnd, npcropmin, nrepr, nlevsno)

Calculate maintenance respiration fluxes (gC/m2/s) for leaf, fine root,
live stem, live coarse root, and crop reproductive tissues.

Ported from `CNMResp` in `CNMRespMod.F90`.
"""
# --- Kernel 1: per-column soil-layer temperature correction factors ------
# tcsoi[c, j] = Q10^((t_soisno[c, j+joff] - TFRZ - 20) / 10) for j = 1:nlevgrnd
@kernel function _mresp_tcsoi_kernel!(tcsoi, @Const(mask_soilc), @Const(t_soisno),
                                      Q10, TFRZ_, joff, nlevgrnd)
    T = eltype(tcsoi)
    c = @index(Global)
    @inbounds if mask_soilc[c]
        for j in 1:nlevgrnd
            tcsoi[c, j] = Q10^((t_soisno[c, j + joff] - TFRZ_ - T(20.0)) / T(10.0))
        end
    end
end

# --- Kernel 2: per-patch leaf and live-wood maintenance respiration ------
@kernel function _mresp_patch_kernel!(leaf_mr, livestem_mr, livecroot_mr,
                                      reproductive_mr, froot_mr,
                                      @Const(mask_soilp), @Const(ivt), @Const(woody),
                                      @Const(frac_veg_nosno),
                                      @Const(lmrsun), @Const(laisun),
                                      @Const(lmrsha), @Const(laisha),
                                      @Const(t_ref2m), @Const(t10),
                                      @Const(livestemn), @Const(livecrootn),
                                      @Const(reproductiven),
                                      Q10, TFRZ_, br, br_root,
                                      rootstem_acc::Bool, npcropmin, nrepr)
    T = eltype(leaf_mr)
    p = @index(Global)
    @inbounds if mask_soilp[p]
        # Temperature correction from 2m air temperature
        tc = Q10^((t_ref2m[p] - TFRZ_ - T(20.0)) / T(10.0))

        # Local copies of base rates (may be modified by acclimation)
        br_local      = br
        br_root_local = br_root

        # Acclimation of root and stem respiration fluxes
        if rootstem_acc
            acc = T(10.0)^(T(-0.00794) * ((t10[p] - TFRZ_) - T(25.0)))
            br_local      = br_local      * acc
            br_root_local = br_root_local * acc
        end

        # Leaf maintenance respiration
        if frac_veg_nosno[p] == 1
            leaf_mr[p] = lmrsun[p] * laisun[p] * T(12.011e-6) +
                         lmrsha[p] * laisha[p] * T(12.011e-6)
        else
            leaf_mr[p] = zero(T)
        end

        # Live stem and live coarse root MR
        if woody[ivt[p] + 1] == one(T)
            livestem_mr[p]  = livestemn[p]  * br_local * tc
            livecroot_mr[p] = livecrootn[p] * br_root_local * tc
        elseif ivt[p] >= npcropmin
            livestem_mr[p] = livestemn[p] * br_local * tc
            for k in 1:nrepr
                reproductive_mr[p, k] = reproductiven[p, k] * br_local * tc
            end
        end

        # --- Initialize fine root MR accumulation ---
        froot_mr[p] = zero(T)
    end
end

# --- Kernel 3: per-patch fine root MR accumulation over soil layers ------
@kernel function _mresp_froot_kernel!(froot_mr, @Const(mask_soilp),
                                      @Const(column), @Const(t10),
                                      @Const(frootn), @Const(crootfr), @Const(tcsoi),
                                      TFRZ_, br_root, rootstem_acc::Bool, nlevgrnd)
    T = eltype(froot_mr)
    p = @index(Global)
    @inbounds if mask_soilp[p]
        c = column[p]

        # Local root base rate (may be modified by acclimation)
        br_root_local = br_root
        if rootstem_acc
            br_root_local = br_root_local * T(10.0)^(T(-0.00794) * ((t10[p] - TFRZ_) - T(25.0)))
        end

        # Ascending-j accumulation into a local, write froot_mr[p] once
        acc = zero(T)
        for j in 1:nlevgrnd
            acc = acc + frootn[p] * br_root_local * tcsoi[c, j] * crootfr[p, j]
        end
        froot_mr[p] = acc
    end
end

function cn_mresp!(mask_soilc::AbstractVector{Bool}, mask_soilp::AbstractVector{Bool},
                   bounds_c::UnitRange{Int}, bounds_p::UnitRange{Int},
                   params::MaintRespParams,
                   cn_params_share::CNSharedParamsData,
                   pftcon::PftConMaintResp,
                   patch::PatchData,
                   canopystate::CanopyStateData,
                   soilstate::SoilStateData,
                   temperature::TemperatureData,
                   photosyns::PhotosynthesisData,
                   cnveg_cf::CNVegCarbonFluxData,
                   cnveg_ns::CNVegNitrogenStateData;
                   nlevgrnd::Int = varpar.nlevgrnd,
                   npcropmin::Int = NPCROPMIN,
                   nrepr::Int = NREPR,
                   nlevsno::Int = varpar.nlevsno)

    # --- Aliases (matching Fortran associate block) ---
    ivt            = patch.itype
    column         = patch.column

    woody          = pftcon.woody

    frac_veg_nosno = canopystate.frac_veg_nosno_patch
    laisun         = canopystate.laisun_patch
    laisha         = canopystate.laisha_patch

    crootfr        = soilstate.crootfr_patch

    # t_soisno is column-level 2D: (ncols, nlevsno+nlevmaxurbgrnd)
    # Fortran indices: (-nlevsno+1:nlevmaxurbgrnd), so soil layer j in Julia
    # is at offset nlevsno + j
    t_soisno       = temperature.t_soisno_col
    t_ref2m        = temperature.t_ref2m_patch
    t10            = temperature.t_a10_patch

    lmrsun         = photosyns.lmrsun_patch
    lmrsha         = photosyns.lmrsha_patch
    rootstem_acc   = photosyns.rootstem_acc

    frootn         = cnveg_ns.frootn_patch
    livestemn      = cnveg_ns.livestemn_patch
    livecrootn     = cnveg_ns.livecrootn_patch
    reproductiven  = cnveg_ns.reproductiven_patch

    leaf_mr        = cnveg_cf.leaf_mr_patch
    froot_mr       = cnveg_cf.froot_mr_patch
    livestem_mr    = cnveg_cf.livestem_mr_patch
    livecroot_mr   = cnveg_cf.livecroot_mr_patch
    reproductive_mr = cnveg_cf.reproductive_mr_patch

    # --- Constants ---
    br      = params.br
    br_root = params.br_root
    Q10     = cn_params_share.Q10

    # Offset for t_soisno column indexing (Fortran -nlevsno+1:.. → Julia 1-based)
    joff = nlevsno

    # --- Column loop: temperature correction factors for each soil layer ---
    # tcsoi scratch (ncols × nlevgrnd), device-resident via similar() off the
    # device-resident t_soisno so the kernel launches on the same backend
    # (a bare Matrix{FT}(undef,…) is host → would force the CPU backend).
    FT = eltype(t_soisno)
    tcsoi = similar(t_soisno, FT, last(bounds_c), nlevgrnd)

    # Convert Float64 scalar args to the working precision so no double is
    # materialized inside the kernel on a Float32-only backend (Metal).
    FTc = eltype(tcsoi)
    _launch!(_mresp_tcsoi_kernel!, tcsoi, mask_soilc, t_soisno,
             FTc(Q10), FTc(TFRZ), joff, nlevgrnd; ndrange = size(tcsoi, 1))

    # --- Patch loop: leaf and live wood maintenance respiration ---
    FTp = eltype(leaf_mr)
    _launch!(_mresp_patch_kernel!, leaf_mr, livestem_mr, livecroot_mr,
             reproductive_mr, froot_mr,
             mask_soilp, ivt, woody, frac_veg_nosno,
             lmrsun, laisun, lmrsha, laisha,
             t_ref2m, t10, livestemn, livecrootn, reproductiven,
             FTp(Q10), FTp(TFRZ), FTp(br), FTp(br_root),
             rootstem_acc, npcropmin, nrepr)

    # --- Soil and patch loop: fine root maintenance respiration ---
    FTf = eltype(froot_mr)
    _launch!(_mresp_froot_kernel!, froot_mr, mask_soilp,
             column, t10, frootn, crootfr, tcsoi,
             FTf(TFRZ), FTf(br_root), rootstem_acc, nlevgrnd)

    return nothing
end
