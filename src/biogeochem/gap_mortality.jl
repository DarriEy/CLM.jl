# ==========================================================================
# Ported from: src/biogeochem/CNGapMortalityMod.F90
# Gap-phase mortality for coupled carbon-nitrogen code (CN).
#
# Public functions:
#   cn_gap_mortality!       — Patch-level gap mortality C & N fluxes
#   cn_gap_patch_to_column! — Aggregate patch mortality fluxes to column
# ==========================================================================

# ---------------------------------------------------------------------------
# Parameters for gap mortality (from CNGapMortalityMod params_type)
# ---------------------------------------------------------------------------

"""
    GapMortalityParams

Module-level parameters for the gap mortality routine.

Fields:
- `k_mort`: coefficient of growth efficiency in the mortality equation
- `r_mort`: per-PFT annual mortality rate (1/yr), indexed by vegetation type

Ported from `params_type` in `CNGapMortalityMod.F90`.
"""
Base.@kwdef mutable struct GapMortalityParams{FT<:Real, V<:AbstractVector{FT}}
    k_mort ::FT        = 0.3     # coeff. of growth efficiency in mortality equation
    r_mort ::V = Float64[]  # mortality rate (1/yr), one per PFT
end
GapMortalityParams{FT}(; kwargs...) where {FT<:Real} = GapMortalityParams{FT, Vector{FT}}(; kwargs...)
Adapt.@adapt_structure GapMortalityParams

# ---------------------------------------------------------------------------
# PFT constants needed by gap mortality (from pftcon)
# ---------------------------------------------------------------------------

"""
    PftConGapMort

PFT-level parameters referenced by the gap mortality and patch-to-column
routines.  Contains the subset of `pftconMod` fields used in
`CNGapMortalityMod.F90`.
"""
Base.@kwdef mutable struct PftConGapMort{FT<:Real, V<:AbstractVector{FT}, M<:AbstractMatrix{FT}}
    woody    ::V = Float64[]  # binary woody flag (1=woody, 0=not woody)
    leafcn   ::V = Float64[]  # leaf C:N (gC/gN)
    livewdcn ::V = Float64[]  # live wood C:N (gC/gN)
    lf_f     ::M = Matrix{Float64}(undef, 0, 0)  # leaf litter fractions (pft, litr)
    fr_f     ::M = Matrix{Float64}(undef, 0, 0)  # fine root litter fractions (pft, litr)
end
PftConGapMort{FT}(; kwargs...) where {FT<:Real} = PftConGapMort{FT, Vector{FT}, Matrix{FT}}(; kwargs...)
Adapt.@adapt_structure PftConGapMort

# ---------------------------------------------------------------------------
# DGVS (Dynamic Global Vegetation) data needed by gap mortality
# ---------------------------------------------------------------------------

"""
    DgvsGapMortData

Dynamic-vegetation fields referenced by gap mortality when `use_cndv` is
enabled.

Ported from `dgvs_type` fields used in `CNGapMortalityMod.F90`.
"""
Base.@kwdef mutable struct DgvsGapMortData{FT<:Real, V<:AbstractVector{FT}}
    greffic_patch    ::V = Float64[]  # growth efficiency
    heatstress_patch ::V = Float64[]  # heat stress mortality
    nind_patch       ::V = Float64[]  # number of individuals (#/m2)
end
DgvsGapMortData{FT}(; kwargs...) where {FT<:Real} = DgvsGapMortData{FT, Vector{FT}}(; kwargs...)
Adapt.@adapt_structure DgvsGapMortData

# ---------------------------------------------------------------------------
# cn_gap_mortality! — Main gap mortality routine
# ---------------------------------------------------------------------------

"""
    cn_gap_mortality!(mask_soilp, bounds, params, pftcon, dgvs,
                      patch, canopystate, cnveg_cs, cnveg_cf, cnveg_ns, cnveg_nf;
                      dt, days_per_year, use_cndv, use_matrixcn, spinup_state,
                      npcropmin, spinup_factor_deadwood, secspday)

Gap-phase mortality routine for coupled carbon-nitrogen code.
Calculates patch-level mortality carbon and nitrogen fluxes.

Ported from `CNGapMortality` in `CNGapMortalityMod.F90`.
"""
# Per-patch gap-mortality kernel: each thread computes the ~20 mortality C & N
# fluxes for one patch p (all independent per-patch outputs) plus the optional
# use_cndv `nind` individual-count update. Config flags (use_cndv, spinup_state,
# npcropmin) are resolved on the host and passed as Bool/Int scalars; no branch
# on a String and no error()/@warn in the kernel.
# =====================================================================
# Device-view structs for the gap-mortality patch kernel (Metal ~31-arg cap).
# Each bundle is ONE kernel arg. Outputs split C/N (20 + 19 fields), inputs
# split C/N (20 + 19 fields), Float64 scalars bundled into _GapScalars{T} at the
# output's working precision. `nind` stays LOOSE (read+written for woody CNDV).
# All are @adapt_structure'd so adapt(MtlArray, …) moves the fields to device.
# =====================================================================

# --- mortality C output flux arrays (20) ---
Base.@kwdef struct _GapMortOutC{V}
    m_leafc_to_litter              ::V
    m_frootc_to_litter             ::V
    m_livestemc_to_litter          ::V
    m_livecrootc_to_litter         ::V
    m_deadstemc_to_litter          ::V
    m_deadcrootc_to_litter         ::V
    m_leafc_storage_to_litter      ::V
    m_frootc_storage_to_litter     ::V
    m_livestemc_storage_to_litter  ::V
    m_deadstemc_storage_to_litter  ::V
    m_livecrootc_storage_to_litter ::V
    m_deadcrootc_storage_to_litter ::V
    m_gresp_storage_to_litter      ::V
    m_leafc_xfer_to_litter         ::V
    m_frootc_xfer_to_litter        ::V
    m_livestemc_xfer_to_litter     ::V
    m_deadstemc_xfer_to_litter     ::V
    m_livecrootc_xfer_to_litter    ::V
    m_deadcrootc_xfer_to_litter    ::V
    m_gresp_xfer_to_litter         ::V
end
Adapt.@adapt_structure _GapMortOutC

# --- mortality N output flux arrays (19) ---
Base.@kwdef struct _GapMortOutN{V}
    m_leafn_to_litter              ::V
    m_frootn_to_litter             ::V
    m_livestemn_to_litter          ::V
    m_livecrootn_to_litter         ::V
    m_deadstemn_to_litter          ::V
    m_deadcrootn_to_litter         ::V
    m_retransn_to_litter           ::V
    m_leafn_storage_to_litter      ::V
    m_frootn_storage_to_litter     ::V
    m_livestemn_storage_to_litter  ::V
    m_deadstemn_storage_to_litter  ::V
    m_livecrootn_storage_to_litter ::V
    m_deadcrootn_storage_to_litter ::V
    m_leafn_xfer_to_litter         ::V
    m_frootn_xfer_to_litter        ::V
    m_livestemn_xfer_to_litter     ::V
    m_deadstemn_xfer_to_litter     ::V
    m_livecrootn_xfer_to_litter    ::V
    m_deadcrootn_xfer_to_litter    ::V
end
Adapt.@adapt_structure _GapMortOutN

# --- C state input arrays (20) ---
Base.@kwdef struct _GapMortInC{V}
    leafc              ::V
    frootc             ::V
    livestemc          ::V
    livecrootc         ::V
    deadstemc          ::V
    deadcrootc         ::V
    leafc_storage      ::V
    frootc_storage     ::V
    livestemc_storage  ::V
    deadstemc_storage  ::V
    livecrootc_storage ::V
    deadcrootc_storage ::V
    gresp_storage      ::V
    leafc_xfer         ::V
    frootc_xfer        ::V
    livestemc_xfer     ::V
    deadstemc_xfer     ::V
    livecrootc_xfer    ::V
    deadcrootc_xfer    ::V
    gresp_xfer         ::V
end
Adapt.@adapt_structure _GapMortInC

# --- N state input arrays (19) ---
Base.@kwdef struct _GapMortInN{V}
    leafn              ::V
    frootn             ::V
    livestemn          ::V
    livecrootn         ::V
    deadstemn          ::V
    deadcrootn         ::V
    retransn           ::V
    leafn_storage      ::V
    frootn_storage     ::V
    livestemn_storage  ::V
    deadstemn_storage  ::V
    livecrootn_storage ::V
    deadcrootn_storage ::V
    leafn_xfer         ::V
    frootn_xfer        ::V
    livestemn_xfer     ::V
    deadstemn_xfer     ::V
    livecrootn_xfer    ::V
    deadcrootn_xfer    ::V
end
Adapt.@adapt_structure _GapMortInN

# --- Float64 scalar params bundled at working precision T (isbits) ---
Base.@kwdef struct _GapScalars{T}
    k_mort                 ::T
    days_per_year          ::T
    secspday               ::T
    spinup_factor_deadwood ::T
end
Adapt.@adapt_structure _GapScalars

# Per-patch gap-mortality kernel: each thread computes the ~20 mortality C & N
# fluxes for one patch p (all independent per-patch outputs) plus the optional
# use_cndv `nind` individual-count update. Array args are grouped into immutable
# device-view structs (outc/outn/inc/inn) to respect the Metal ~31-arg cap; the
# four Float64 mortality-rate/spinup scalars are bundled into _GapScalars{T} at the
# output's working precision so NO Float64 is materialized on a Float32-only backend.
# Every field is aliased to its Fortran-named local at the kernel top so the physics
# body is verbatim except numeric literal -> T(...). Config flags (use_cndv,
# spinup_state, npcropmin) stay loose Bool/Int scalars; no String branch, no error().
@kernel function _gapmort_patch_kernel!(
        outc::_GapMortOutC, outn::_GapMortOutN, nind,
        @Const(mask_soilp), @Const(ivt), @Const(woody),
        @Const(greffic), @Const(heatstress), @Const(r_mort),
        inc::_GapMortInC, inn::_GapMortInN, s::_GapScalars,
        use_cndv::Bool, spinup_state::Int, npcropmin::Int)

    # Working element type taken from an output array (Metal: must be Float32).
    T = eltype(outc.m_leafc_to_litter)

    # ---- alias output flux arrays (verbatim Fortran names) ----
    m_leafc_to_litter              = outc.m_leafc_to_litter
    m_frootc_to_litter             = outc.m_frootc_to_litter
    m_livestemc_to_litter          = outc.m_livestemc_to_litter
    m_livecrootc_to_litter         = outc.m_livecrootc_to_litter
    m_deadstemc_to_litter          = outc.m_deadstemc_to_litter
    m_deadcrootc_to_litter         = outc.m_deadcrootc_to_litter
    m_leafc_storage_to_litter      = outc.m_leafc_storage_to_litter
    m_frootc_storage_to_litter     = outc.m_frootc_storage_to_litter
    m_livestemc_storage_to_litter  = outc.m_livestemc_storage_to_litter
    m_deadstemc_storage_to_litter  = outc.m_deadstemc_storage_to_litter
    m_livecrootc_storage_to_litter = outc.m_livecrootc_storage_to_litter
    m_deadcrootc_storage_to_litter = outc.m_deadcrootc_storage_to_litter
    m_gresp_storage_to_litter      = outc.m_gresp_storage_to_litter
    m_leafc_xfer_to_litter         = outc.m_leafc_xfer_to_litter
    m_frootc_xfer_to_litter        = outc.m_frootc_xfer_to_litter
    m_livestemc_xfer_to_litter     = outc.m_livestemc_xfer_to_litter
    m_deadstemc_xfer_to_litter     = outc.m_deadstemc_xfer_to_litter
    m_livecrootc_xfer_to_litter    = outc.m_livecrootc_xfer_to_litter
    m_deadcrootc_xfer_to_litter    = outc.m_deadcrootc_xfer_to_litter
    m_gresp_xfer_to_litter         = outc.m_gresp_xfer_to_litter

    m_leafn_to_litter              = outn.m_leafn_to_litter
    m_frootn_to_litter             = outn.m_frootn_to_litter
    m_livestemn_to_litter          = outn.m_livestemn_to_litter
    m_livecrootn_to_litter         = outn.m_livecrootn_to_litter
    m_deadstemn_to_litter          = outn.m_deadstemn_to_litter
    m_deadcrootn_to_litter         = outn.m_deadcrootn_to_litter
    m_retransn_to_litter           = outn.m_retransn_to_litter
    m_leafn_storage_to_litter      = outn.m_leafn_storage_to_litter
    m_frootn_storage_to_litter     = outn.m_frootn_storage_to_litter
    m_livestemn_storage_to_litter  = outn.m_livestemn_storage_to_litter
    m_deadstemn_storage_to_litter  = outn.m_deadstemn_storage_to_litter
    m_livecrootn_storage_to_litter = outn.m_livecrootn_storage_to_litter
    m_deadcrootn_storage_to_litter = outn.m_deadcrootn_storage_to_litter
    m_leafn_xfer_to_litter         = outn.m_leafn_xfer_to_litter
    m_frootn_xfer_to_litter        = outn.m_frootn_xfer_to_litter
    m_livestemn_xfer_to_litter     = outn.m_livestemn_xfer_to_litter
    m_deadstemn_xfer_to_litter     = outn.m_deadstemn_xfer_to_litter
    m_livecrootn_xfer_to_litter    = outn.m_livecrootn_xfer_to_litter
    m_deadcrootn_xfer_to_litter    = outn.m_deadcrootn_xfer_to_litter

    # ---- alias C state input arrays ----
    leafc              = inc.leafc
    frootc             = inc.frootc
    livestemc          = inc.livestemc
    livecrootc         = inc.livecrootc
    deadstemc          = inc.deadstemc
    deadcrootc         = inc.deadcrootc
    leafc_storage      = inc.leafc_storage
    frootc_storage     = inc.frootc_storage
    livestemc_storage  = inc.livestemc_storage
    deadstemc_storage  = inc.deadstemc_storage
    livecrootc_storage = inc.livecrootc_storage
    deadcrootc_storage = inc.deadcrootc_storage
    gresp_storage      = inc.gresp_storage
    leafc_xfer         = inc.leafc_xfer
    frootc_xfer        = inc.frootc_xfer
    livestemc_xfer     = inc.livestemc_xfer
    deadstemc_xfer     = inc.deadstemc_xfer
    livecrootc_xfer    = inc.livecrootc_xfer
    deadcrootc_xfer    = inc.deadcrootc_xfer
    gresp_xfer         = inc.gresp_xfer

    # ---- alias N state input arrays ----
    leafn              = inn.leafn
    frootn             = inn.frootn
    livestemn          = inn.livestemn
    livecrootn         = inn.livecrootn
    deadstemn          = inn.deadstemn
    deadcrootn         = inn.deadcrootn
    retransn           = inn.retransn
    leafn_storage      = inn.leafn_storage
    frootn_storage     = inn.frootn_storage
    livestemn_storage  = inn.livestemn_storage
    deadstemn_storage  = inn.deadstemn_storage
    livecrootn_storage = inn.livecrootn_storage
    deadcrootn_storage = inn.deadcrootn_storage
    leafn_xfer         = inn.leafn_xfer
    frootn_xfer        = inn.frootn_xfer
    livestemn_xfer     = inn.livestemn_xfer
    deadstemn_xfer     = inn.deadstemn_xfer
    livecrootn_xfer    = inn.livecrootn_xfer
    deadcrootn_xfer    = inn.deadcrootn_xfer

    # ---- alias scalar params (already at precision T) ----
    k_mort                 = s.k_mort
    days_per_year          = s.days_per_year
    secspday               = s.secspday
    spinup_factor_deadwood = s.spinup_factor_deadwood

    p = @index(Global)
    @inbounds if mask_soilp[p]
        ivtp = ivt[p] + 1   # 0-based Fortran → 1-based Julia

        if use_cndv
            # Stress mortality from lpj's subr Mortality
            if woody[ivtp] == one(T)
                if ivt[p] == 8
                    mort_max = T(0.03)  # BDT boreal
                else
                    mort_max = T(0.01)  # original value for all patches
                end
                # heatstress and greffic calculated in Establishment once/yr
                # Mortality rate inversely related to growth efficiency
                # (Prentice et al 1993)
                am = mort_max / (one(T) + k_mort * greffic[p])
                am = min(one(T), am + heatstress[p])
            else
                # lpj didn't set this for grasses; cn does
                am = r_mort[ivtp]
            end
        else
            am = r_mort[ivtp]
        end

        m = am / (days_per_year * secspday)

        # --------------------------------------------------
        # patch-level gap mortality carbon fluxes
        # --------------------------------------------------

        # displayed pools
        m_leafc_to_litter[p]      = leafc[p]      * m
        m_frootc_to_litter[p]     = frootc[p]     * m
        m_livestemc_to_litter[p]  = livestemc[p]  * m
        m_livecrootc_to_litter[p] = livecrootc[p] * m

        m_deadstemc_to_litter[p]  = deadstemc[p]  * m * spinup_factor_deadwood
        m_deadcrootc_to_litter[p] = deadcrootc[p] * m * spinup_factor_deadwood

        # storage pools
        m_leafc_storage_to_litter[p]      = leafc_storage[p]      * m
        m_frootc_storage_to_litter[p]     = frootc_storage[p]     * m
        m_livestemc_storage_to_litter[p]  = livestemc_storage[p]  * m
        m_deadstemc_storage_to_litter[p]  = deadstemc_storage[p]  * m
        m_livecrootc_storage_to_litter[p] = livecrootc_storage[p] * m
        m_deadcrootc_storage_to_litter[p] = deadcrootc_storage[p] * m
        m_gresp_storage_to_litter[p]      = gresp_storage[p]      * m

        # transfer pools
        m_leafc_xfer_to_litter[p]      = leafc_xfer[p]      * m
        m_frootc_xfer_to_litter[p]     = frootc_xfer[p]     * m
        m_livestemc_xfer_to_litter[p]  = livestemc_xfer[p]  * m
        m_deadstemc_xfer_to_litter[p]  = deadstemc_xfer[p]  * m
        m_livecrootc_xfer_to_litter[p] = livecrootc_xfer[p] * m
        m_deadcrootc_xfer_to_litter[p] = deadcrootc_xfer[p] * m
        m_gresp_xfer_to_litter[p]      = gresp_xfer[p]      * m

        # --------------------------------------------------
        # patch-level gap mortality nitrogen fluxes
        # --------------------------------------------------

        # displayed pools
        m_leafn_to_litter[p]      = leafn[p]      * m
        m_frootn_to_litter[p]     = frootn[p]     * m
        m_livestemn_to_litter[p]  = livestemn[p]  * m
        m_livecrootn_to_litter[p] = livecrootn[p] * m

        if spinup_state == 2 && !use_cndv
            # accelerate mortality of dead woody pools
            m_deadstemn_to_litter[p]  = deadstemn[p]  * m * spinup_factor_deadwood
            m_deadcrootn_to_litter[p] = deadcrootn[p] * m * spinup_factor_deadwood
        else
            m_deadstemn_to_litter[p]  = deadstemn[p]  * m
            m_deadcrootn_to_litter[p] = deadcrootn[p] * m
        end

        if ivt[p] < npcropmin
            m_retransn_to_litter[p] = retransn[p] * m
        end

        # storage pools
        m_leafn_storage_to_litter[p]      = leafn_storage[p]      * m
        m_frootn_storage_to_litter[p]     = frootn_storage[p]     * m
        m_livestemn_storage_to_litter[p]  = livestemn_storage[p]  * m
        m_deadstemn_storage_to_litter[p]  = deadstemn_storage[p]  * m
        m_livecrootn_storage_to_litter[p] = livecrootn_storage[p] * m
        m_deadcrootn_storage_to_litter[p] = deadcrootn_storage[p] * m

        # transfer pools
        m_leafn_xfer_to_litter[p]      = leafn_xfer[p]      * m
        m_frootn_xfer_to_litter[p]     = frootn_xfer[p]     * m
        m_livestemn_xfer_to_litter[p]  = livestemn_xfer[p]  * m
        m_deadstemn_xfer_to_litter[p]  = deadstemn_xfer[p]  * m
        m_livecrootn_xfer_to_litter[p] = livecrootn_xfer[p] * m
        m_deadcrootn_xfer_to_litter[p] = deadcrootn_xfer[p] * m

        # added by F. Li and S. Levis
        if use_cndv
            if woody[ivtp] == one(T)
                if livestemc[p] + deadstemc[p] > zero(T)
                    nind[p] = nind[p] * (one(T) - m)
                else
                    nind[p] = zero(T)
                end
            end
        end
    end
end

function cn_gap_mortality!(mask_soilp::AbstractVector{Bool},
                           bounds::UnitRange{Int},
                           params::GapMortalityParams,
                           pftcon::PftConGapMort,
                           dgvs::DgvsGapMortData,
                           patch::PatchData,
                           canopystate::CanopyStateData,
                           cnveg_cs::CNVegCarbonStateData,
                           cnveg_cf::CNVegCarbonFluxData,
                           cnveg_ns::CNVegNitrogenStateData,
                           cnveg_nf::CNVegNitrogenFluxData;
                           dt::Real = 1800.0,
                           days_per_year::Real = 365.0,
                           use_cndv::Bool = false,
                           use_matrixcn::Bool = false,
                           spinup_state::Int = 0,
                           npcropmin::Int = 17,
                           spinup_factor_deadwood::Real = SPINUP_FACTOR_DEADWOOD_DEFAULT,
                           secspday::Real = SECSPDAY)

    # The Fortran loop runs over a filtered patch range `bounds`; the mask is
    # active only on that range, so launching one thread per patch over the full
    # mask reproduces it (the kernel skips inactive patches). bounds == 1:np here.

    # Working precision taken from an output array (Float32 on Metal, Float64 on host).
    T = eltype(cnveg_cf.m_leafc_to_litter_patch)

    outc = _GapMortOutC(;
        m_leafc_to_litter              = cnveg_cf.m_leafc_to_litter_patch,
        m_frootc_to_litter             = cnveg_cf.m_frootc_to_litter_patch,
        m_livestemc_to_litter          = cnveg_cf.m_livestemc_to_litter_patch,
        m_livecrootc_to_litter         = cnveg_cf.m_livecrootc_to_litter_patch,
        m_deadstemc_to_litter          = cnveg_cf.m_deadstemc_to_litter_patch,
        m_deadcrootc_to_litter         = cnveg_cf.m_deadcrootc_to_litter_patch,
        m_leafc_storage_to_litter      = cnveg_cf.m_leafc_storage_to_litter_patch,
        m_frootc_storage_to_litter     = cnveg_cf.m_frootc_storage_to_litter_patch,
        m_livestemc_storage_to_litter  = cnveg_cf.m_livestemc_storage_to_litter_patch,
        m_deadstemc_storage_to_litter  = cnveg_cf.m_deadstemc_storage_to_litter_patch,
        m_livecrootc_storage_to_litter = cnveg_cf.m_livecrootc_storage_to_litter_patch,
        m_deadcrootc_storage_to_litter = cnveg_cf.m_deadcrootc_storage_to_litter_patch,
        m_gresp_storage_to_litter      = cnveg_cf.m_gresp_storage_to_litter_patch,
        m_leafc_xfer_to_litter         = cnveg_cf.m_leafc_xfer_to_litter_patch,
        m_frootc_xfer_to_litter        = cnveg_cf.m_frootc_xfer_to_litter_patch,
        m_livestemc_xfer_to_litter     = cnveg_cf.m_livestemc_xfer_to_litter_patch,
        m_deadstemc_xfer_to_litter     = cnveg_cf.m_deadstemc_xfer_to_litter_patch,
        m_livecrootc_xfer_to_litter    = cnveg_cf.m_livecrootc_xfer_to_litter_patch,
        m_deadcrootc_xfer_to_litter    = cnveg_cf.m_deadcrootc_xfer_to_litter_patch,
        m_gresp_xfer_to_litter         = cnveg_cf.m_gresp_xfer_to_litter_patch)

    outn = _GapMortOutN(;
        m_leafn_to_litter              = cnveg_nf.m_leafn_to_litter_patch,
        m_frootn_to_litter             = cnveg_nf.m_frootn_to_litter_patch,
        m_livestemn_to_litter          = cnveg_nf.m_livestemn_to_litter_patch,
        m_livecrootn_to_litter         = cnveg_nf.m_livecrootn_to_litter_patch,
        m_deadstemn_to_litter          = cnveg_nf.m_deadstemn_to_litter_patch,
        m_deadcrootn_to_litter         = cnveg_nf.m_deadcrootn_to_litter_patch,
        m_retransn_to_litter           = cnveg_nf.m_retransn_to_litter_patch,
        m_leafn_storage_to_litter      = cnveg_nf.m_leafn_storage_to_litter_patch,
        m_frootn_storage_to_litter     = cnveg_nf.m_frootn_storage_to_litter_patch,
        m_livestemn_storage_to_litter  = cnveg_nf.m_livestemn_storage_to_litter_patch,
        m_deadstemn_storage_to_litter  = cnveg_nf.m_deadstemn_storage_to_litter_patch,
        m_livecrootn_storage_to_litter = cnveg_nf.m_livecrootn_storage_to_litter_patch,
        m_deadcrootn_storage_to_litter = cnveg_nf.m_deadcrootn_storage_to_litter_patch,
        m_leafn_xfer_to_litter         = cnveg_nf.m_leafn_xfer_to_litter_patch,
        m_frootn_xfer_to_litter        = cnveg_nf.m_frootn_xfer_to_litter_patch,
        m_livestemn_xfer_to_litter     = cnveg_nf.m_livestemn_xfer_to_litter_patch,
        m_deadstemn_xfer_to_litter     = cnveg_nf.m_deadstemn_xfer_to_litter_patch,
        m_livecrootn_xfer_to_litter    = cnveg_nf.m_livecrootn_xfer_to_litter_patch,
        m_deadcrootn_xfer_to_litter    = cnveg_nf.m_deadcrootn_xfer_to_litter_patch)

    inc = _GapMortInC(;
        leafc              = cnveg_cs.leafc_patch,
        frootc             = cnveg_cs.frootc_patch,
        livestemc          = cnveg_cs.livestemc_patch,
        livecrootc         = cnveg_cs.livecrootc_patch,
        deadstemc          = cnveg_cs.deadstemc_patch,
        deadcrootc         = cnveg_cs.deadcrootc_patch,
        leafc_storage      = cnveg_cs.leafc_storage_patch,
        frootc_storage     = cnveg_cs.frootc_storage_patch,
        livestemc_storage  = cnveg_cs.livestemc_storage_patch,
        deadstemc_storage  = cnveg_cs.deadstemc_storage_patch,
        livecrootc_storage = cnveg_cs.livecrootc_storage_patch,
        deadcrootc_storage = cnveg_cs.deadcrootc_storage_patch,
        gresp_storage      = cnveg_cs.gresp_storage_patch,
        leafc_xfer         = cnveg_cs.leafc_xfer_patch,
        frootc_xfer        = cnveg_cs.frootc_xfer_patch,
        livestemc_xfer     = cnveg_cs.livestemc_xfer_patch,
        deadstemc_xfer     = cnveg_cs.deadstemc_xfer_patch,
        livecrootc_xfer    = cnveg_cs.livecrootc_xfer_patch,
        deadcrootc_xfer    = cnveg_cs.deadcrootc_xfer_patch,
        gresp_xfer         = cnveg_cs.gresp_xfer_patch)

    inn = _GapMortInN(;
        leafn              = cnveg_ns.leafn_patch,
        frootn             = cnveg_ns.frootn_patch,
        livestemn          = cnveg_ns.livestemn_patch,
        livecrootn         = cnveg_ns.livecrootn_patch,
        deadstemn          = cnveg_ns.deadstemn_patch,
        deadcrootn         = cnveg_ns.deadcrootn_patch,
        retransn           = cnveg_ns.retransn_patch,
        leafn_storage      = cnveg_ns.leafn_storage_patch,
        frootn_storage     = cnveg_ns.frootn_storage_patch,
        livestemn_storage  = cnveg_ns.livestemn_storage_patch,
        deadstemn_storage  = cnveg_ns.deadstemn_storage_patch,
        livecrootn_storage = cnveg_ns.livecrootn_storage_patch,
        deadcrootn_storage = cnveg_ns.deadcrootn_storage_patch,
        leafn_xfer         = cnveg_ns.leafn_xfer_patch,
        frootn_xfer        = cnveg_ns.frootn_xfer_patch,
        livestemn_xfer     = cnveg_ns.livestemn_xfer_patch,
        deadstemn_xfer     = cnveg_ns.deadstemn_xfer_patch,
        livecrootn_xfer    = cnveg_ns.livecrootn_xfer_patch,
        deadcrootn_xfer    = cnveg_ns.deadcrootn_xfer_patch)

    # Float64 scalar params converted to working precision T (no Float64 in-kernel).
    s = _GapScalars{T}(;
        k_mort                 = T(params.k_mort),
        days_per_year          = T(days_per_year),
        secspday               = T(secspday),
        spinup_factor_deadwood = T(spinup_factor_deadwood))

    # Struct-first kernel: the bundle args carry no backend, so launch manually off
    # an output array's backend and synchronize. ndrange = patch count (1:np).
    np = length(mask_soilp)
    if np != 0
        backend = _kernel_backend(cnveg_cf.m_leafc_to_litter_patch)
        _gapmort_patch_kernel!(backend)(
            outc, outn, dgvs.nind_patch,
            mask_soilp, patch.itype, pftcon.woody,
            dgvs.greffic_patch, dgvs.heatstress_patch, params.r_mort,
            inc, inn, s,
            use_cndv, spinup_state, npcropmin; ndrange = np)
        KA.synchronize(backend)
    end

    return nothing
end

# ---------------------------------------------------------------------------
# cn_gap_patch_to_column! — Gather patch mortality fluxes to column level
# ---------------------------------------------------------------------------

"""
    cn_gap_patch_to_column!(mask_soilp, bounds, pftcon, patch,
                            cnveg_cf, cnveg_nf,
                            leaf_prof, froot_prof, croot_prof, stem_prof;
                            nlevdecomp, i_litr_min, i_litr_max, i_met_lit)

Gathers all patch-level gap mortality fluxes to the column level and assigns
them to the appropriate litter and coarse woody debris (CWD) pools.

Ported from `CNGap_PatchToColumn` in `CNGapMortalityMod.F90`.
"""
# Per-patch patch→column scatter kernel: each thread handles one patch p and
# scatters its gap-mortality C & N fluxes into the column litter / CWD pools over
# all decomp levels j and litter pools i. The host loop was j-outer/p-inner; for a
# fixed (c,j,i) accumulator the contributions are summed in ascending-p order in
# both the host loop and the KA-CPU launch (threads iterate p ascending), so the
# floating-point accumulation order — and the result — is byte-identical on CPU.
# The two C-to-CWD (stem then croot) and the two N-to-CWD adds are kept as two
# separate _scatter_add! calls each, matching the host's two separate `+=`.
## --- Device-view bundles for the gap-mortality patch->column scatter ---
## Outputs: column litter (3D: col,level,litr-pool) + CWD (2D: col,level), C and N.
Base.@kwdef struct _GapP2COut{V,M}
    c_to_litr_c::M      # gap_mortality_c_to_litr_c_col  (col, level, litr)
    c_to_cwdc::V        # gap_mortality_c_to_cwdc_col     (col, level)
    n_to_litr_n::M      # gap_mortality_n_to_litr_n_col   (col, level, litr)
    n_to_cwdn::V        # gap_mortality_n_to_cwdn_col      (col, level)
end
Adapt.@adapt_structure _GapP2COut

## Profile matrices + per-PFT litter fractions (all 2D).
Base.@kwdef struct _GapP2CProf{M}
    lf_f::M; fr_f::M
    leaf_prof::M; froot_prof::M; croot_prof::M; stem_prof::M
end
Adapt.@adapt_structure _GapP2CProf

## Per-patch gap-mortality CARBON fluxes (all vectors).
Base.@kwdef struct _GapP2CFluxC{V}
    leafc::V; frootc::V; livestemc::V; deadstemc::V; livecrootc::V; deadcrootc::V
    leafc_storage::V; frootc_storage::V; livestemc_storage::V; deadstemc_storage::V
    livecrootc_storage::V; deadcrootc_storage::V; gresp_storage::V
    leafc_xfer::V; frootc_xfer::V; livestemc_xfer::V; deadstemc_xfer::V
    livecrootc_xfer::V; deadcrootc_xfer::V; gresp_xfer::V
end
Adapt.@adapt_structure _GapP2CFluxC

## Per-patch gap-mortality NITROGEN fluxes (all vectors).
Base.@kwdef struct _GapP2CFluxN{V}
    leafn::V; frootn::V; livestemn::V; deadstemn::V; livecrootn::V; deadcrootn::V
    retransn::V
    leafn_storage::V; frootn_storage::V; livestemn_storage::V; deadstemn_storage::V
    livecrootn_storage::V; deadcrootn_storage::V
    leafn_xfer::V; frootn_xfer::V; livestemn_xfer::V; deadstemn_xfer::V
    livecrootn_xfer::V; deadcrootn_xfer::V
end
Adapt.@adapt_structure _GapP2CFluxN

@kernel function _gap_p2c_kernel!(
        out::_GapP2COut, prof::_GapP2CProf,
        fc::_GapP2CFluxC, fn::_GapP2CFluxN,
        @Const(mask_soilp), @Const(ivt), @Const(column), @Const(wtcol),
        nlevdecomp::Int, i_litr_min::Int, i_litr_max::Int, i_met_lit::Int)

    p = @index(Global)
    @inbounds if mask_soilp[p]
        # Alias every device-view field to its Fortran-named local so the body
        # below stays verbatim w.r.t. the original loose-arg kernel.
        gap_mortality_c_to_litr_c = out.c_to_litr_c
        gap_mortality_c_to_cwdc   = out.c_to_cwdc
        gap_mortality_n_to_litr_n = out.n_to_litr_n
        gap_mortality_n_to_cwdn   = out.n_to_cwdn

        lf_f       = prof.lf_f
        fr_f       = prof.fr_f
        leaf_prof  = prof.leaf_prof
        froot_prof = prof.froot_prof
        croot_prof = prof.croot_prof
        stem_prof  = prof.stem_prof

        m_leafc_to_litter              = fc.leafc
        m_frootc_to_litter             = fc.frootc
        m_livestemc_to_litter          = fc.livestemc
        m_deadstemc_to_litter          = fc.deadstemc
        m_livecrootc_to_litter         = fc.livecrootc
        m_deadcrootc_to_litter         = fc.deadcrootc
        m_leafc_storage_to_litter      = fc.leafc_storage
        m_frootc_storage_to_litter     = fc.frootc_storage
        m_livestemc_storage_to_litter  = fc.livestemc_storage
        m_deadstemc_storage_to_litter  = fc.deadstemc_storage
        m_livecrootc_storage_to_litter = fc.livecrootc_storage
        m_deadcrootc_storage_to_litter = fc.deadcrootc_storage
        m_gresp_storage_to_litter      = fc.gresp_storage
        m_leafc_xfer_to_litter         = fc.leafc_xfer
        m_frootc_xfer_to_litter        = fc.frootc_xfer
        m_livestemc_xfer_to_litter     = fc.livestemc_xfer
        m_deadstemc_xfer_to_litter     = fc.deadstemc_xfer
        m_livecrootc_xfer_to_litter    = fc.livecrootc_xfer
        m_deadcrootc_xfer_to_litter    = fc.deadcrootc_xfer
        m_gresp_xfer_to_litter         = fc.gresp_xfer

        m_leafn_to_litter              = fn.leafn
        m_frootn_to_litter             = fn.frootn
        m_livestemn_to_litter          = fn.livestemn
        m_deadstemn_to_litter          = fn.deadstemn
        m_livecrootn_to_litter         = fn.livecrootn
        m_deadcrootn_to_litter         = fn.deadcrootn
        m_retransn_to_litter           = fn.retransn
        m_leafn_storage_to_litter      = fn.leafn_storage
        m_frootn_storage_to_litter     = fn.frootn_storage
        m_livestemn_storage_to_litter  = fn.livestemn_storage
        m_deadstemn_storage_to_litter  = fn.deadstemn_storage
        m_livecrootn_storage_to_litter = fn.livecrootn_storage
        m_deadcrootn_storage_to_litter = fn.deadcrootn_storage
        m_leafn_xfer_to_litter         = fn.leafn_xfer
        m_frootn_xfer_to_litter        = fn.frootn_xfer
        m_livestemn_xfer_to_litter     = fn.livestemn_xfer
        m_deadstemn_xfer_to_litter     = fn.deadstemn_xfer
        m_livecrootn_xfer_to_litter    = fn.livecrootn_xfer
        m_deadcrootn_xfer_to_litter    = fn.deadcrootn_xfer

        ivtp = ivt[p] + 1
        c    = column[p]
        wt   = wtcol[p]

        for j in 1:nlevdecomp
            lp = leaf_prof[p, j]
            fp = froot_prof[p, j]
            cp = croot_prof[p, j]
            sp = stem_prof[p, j]

            # -- Leaf and fine root gap mortality C to litter pools --
            for i in i_litr_min:i_litr_max
                _scatter_add!(gap_mortality_c_to_litr_c, c, j, i,
                    m_leafc_to_litter[p]  * lf_f[ivtp, i] * wt * lp +
                    m_frootc_to_litter[p] * fr_f[ivtp, i] * wt * fp)
            end

            # -- Wood gap mortality C to CWD --
            _scatter_add!(gap_mortality_c_to_cwdc, c, j,
                (m_livestemc_to_litter[p] + m_deadstemc_to_litter[p]) * wt * sp)
            _scatter_add!(gap_mortality_c_to_cwdc, c, j,
                (m_livecrootc_to_litter[p] + m_deadcrootc_to_litter[p]) * wt * cp)

            # -- Storage and transfer gap mortality C to metabolic litter --
            _scatter_add!(gap_mortality_c_to_litr_c, c, j, i_met_lit,
                (m_leafc_storage_to_litter[p] + m_gresp_storage_to_litter[p]) * wt * lp +
                m_frootc_storage_to_litter[p] * wt * fp +
                (m_livestemc_storage_to_litter[p] + m_deadstemc_storage_to_litter[p]) * wt * sp +
                (m_livecrootc_storage_to_litter[p] + m_deadcrootc_storage_to_litter[p]) * wt * cp +
                # transfer gap mortality carbon fluxes
                (m_leafc_xfer_to_litter[p] + m_gresp_xfer_to_litter[p]) * wt * lp +
                m_frootc_xfer_to_litter[p] * wt * fp +
                (m_livestemc_xfer_to_litter[p] + m_deadstemc_xfer_to_litter[p]) * wt * sp +
                (m_livecrootc_xfer_to_litter[p] + m_deadcrootc_xfer_to_litter[p]) * wt * cp)

            # -- Leaf and fine root gap mortality N to litter pools --
            for i in i_litr_min:i_litr_max
                _scatter_add!(gap_mortality_n_to_litr_n, c, j, i,
                    m_leafn_to_litter[p]  * lf_f[ivtp, i] * wt * lp +
                    m_frootn_to_litter[p] * fr_f[ivtp, i] * wt * fp)
            end

            # -- Wood gap mortality N to CWD --
            _scatter_add!(gap_mortality_n_to_cwdn, c, j,
                (m_livestemn_to_litter[p] + m_deadstemn_to_litter[p]) * wt * sp)
            _scatter_add!(gap_mortality_n_to_cwdn, c, j,
                (m_livecrootn_to_litter[p] + m_deadcrootn_to_litter[p]) * wt * cp)

            # -- Storage, transfer, and retranslocated N to metabolic litter --
            _scatter_add!(gap_mortality_n_to_litr_n, c, j, i_met_lit,
                # retranslocated N pool gap mortality fluxes
                m_retransn_to_litter[p] * wt * lp +
                # storage gap mortality nitrogen fluxes
                m_leafn_storage_to_litter[p] * wt * lp +
                m_frootn_storage_to_litter[p] * wt * fp +
                (m_livestemn_storage_to_litter[p] + m_deadstemn_storage_to_litter[p]) * wt * sp +
                (m_livecrootn_storage_to_litter[p] + m_deadcrootn_storage_to_litter[p]) * wt * cp +
                # transfer gap mortality nitrogen fluxes
                m_leafn_xfer_to_litter[p] * wt * lp +
                m_frootn_xfer_to_litter[p] * wt * fp +
                (m_livestemn_xfer_to_litter[p] + m_deadstemn_xfer_to_litter[p]) * wt * sp +
                (m_livecrootn_xfer_to_litter[p] + m_deadcrootn_xfer_to_litter[p]) * wt * cp)
        end
    end
end

function cn_gap_patch_to_column!(mask_soilp::AbstractVector{Bool},
                                 bounds::UnitRange{Int},
                                 pftcon::PftConGapMort,
                                 patch::PatchData,
                                 cnveg_cf::CNVegCarbonFluxData,
                                 cnveg_nf::CNVegNitrogenFluxData,
                                 leaf_prof::AbstractMatrix{<:Real},
                                 froot_prof::AbstractMatrix{<:Real},
                                 croot_prof::AbstractMatrix{<:Real},
                                 stem_prof::AbstractMatrix{<:Real};
                                 nlevdecomp::Int = 1,
                                 i_litr_min::Int = 1,
                                 i_litr_max::Int = 3,
                                 i_met_lit::Int = 1)

    # The column accumulators are NOT zeroed here (the caller manages that); the
    # kernel only does `+=` scatters, matching the host loop's accumulation.
    #
    # Metal's ~31-arg limit forbids passing ~52 arrays loose, so they are grouped
    # into immutable Adapt-able device-view bundles (one kernel arg each). Index
    # arrays (mask/itype/column) + the Int loop bounds stay loose. The kernel body
    # aliases each bundle field back to its Fortran-named local, so the arithmetic
    # — and the CPU float accumulation order — is byte-identical to the loose form.
    out = _GapP2COut(;
        c_to_litr_c = cnveg_cf.gap_mortality_c_to_litr_c_col,
        c_to_cwdc   = cnveg_cf.gap_mortality_c_to_cwdc_col,
        n_to_litr_n = cnveg_nf.gap_mortality_n_to_litr_n_col,
        n_to_cwdn   = cnveg_nf.gap_mortality_n_to_cwdn_col)

    prof = _GapP2CProf(;
        lf_f = pftcon.lf_f, fr_f = pftcon.fr_f,
        leaf_prof = leaf_prof, froot_prof = froot_prof,
        croot_prof = croot_prof, stem_prof = stem_prof)

    fc = _GapP2CFluxC(;
        leafc              = cnveg_cf.m_leafc_to_litter_patch,
        frootc             = cnveg_cf.m_frootc_to_litter_patch,
        livestemc          = cnveg_cf.m_livestemc_to_litter_patch,
        deadstemc          = cnveg_cf.m_deadstemc_to_litter_patch,
        livecrootc         = cnveg_cf.m_livecrootc_to_litter_patch,
        deadcrootc         = cnveg_cf.m_deadcrootc_to_litter_patch,
        leafc_storage      = cnveg_cf.m_leafc_storage_to_litter_patch,
        frootc_storage     = cnveg_cf.m_frootc_storage_to_litter_patch,
        livestemc_storage  = cnveg_cf.m_livestemc_storage_to_litter_patch,
        deadstemc_storage  = cnveg_cf.m_deadstemc_storage_to_litter_patch,
        livecrootc_storage = cnveg_cf.m_livecrootc_storage_to_litter_patch,
        deadcrootc_storage = cnveg_cf.m_deadcrootc_storage_to_litter_patch,
        gresp_storage      = cnveg_cf.m_gresp_storage_to_litter_patch,
        leafc_xfer         = cnveg_cf.m_leafc_xfer_to_litter_patch,
        frootc_xfer        = cnveg_cf.m_frootc_xfer_to_litter_patch,
        livestemc_xfer     = cnveg_cf.m_livestemc_xfer_to_litter_patch,
        deadstemc_xfer     = cnveg_cf.m_deadstemc_xfer_to_litter_patch,
        livecrootc_xfer    = cnveg_cf.m_livecrootc_xfer_to_litter_patch,
        deadcrootc_xfer    = cnveg_cf.m_deadcrootc_xfer_to_litter_patch,
        gresp_xfer         = cnveg_cf.m_gresp_xfer_to_litter_patch)

    fn = _GapP2CFluxN(;
        leafn              = cnveg_nf.m_leafn_to_litter_patch,
        frootn             = cnveg_nf.m_frootn_to_litter_patch,
        livestemn          = cnveg_nf.m_livestemn_to_litter_patch,
        deadstemn          = cnveg_nf.m_deadstemn_to_litter_patch,
        livecrootn         = cnveg_nf.m_livecrootn_to_litter_patch,
        deadcrootn         = cnveg_nf.m_deadcrootn_to_litter_patch,
        retransn           = cnveg_nf.m_retransn_to_litter_patch,
        leafn_storage      = cnveg_nf.m_leafn_storage_to_litter_patch,
        frootn_storage     = cnveg_nf.m_frootn_storage_to_litter_patch,
        livestemn_storage  = cnveg_nf.m_livestemn_storage_to_litter_patch,
        deadstemn_storage  = cnveg_nf.m_deadstemn_storage_to_litter_patch,
        livecrootn_storage = cnveg_nf.m_livecrootn_storage_to_litter_patch,
        deadcrootn_storage = cnveg_nf.m_deadcrootn_storage_to_litter_patch,
        leafn_xfer         = cnveg_nf.m_leafn_xfer_to_litter_patch,
        frootn_xfer        = cnveg_nf.m_frootn_xfer_to_litter_patch,
        livestemn_xfer     = cnveg_nf.m_livestemn_xfer_to_litter_patch,
        deadstemn_xfer     = cnveg_nf.m_deadstemn_xfer_to_litter_patch,
        livecrootn_xfer    = cnveg_nf.m_livecrootn_xfer_to_litter_patch,
        deadcrootn_xfer    = cnveg_nf.m_deadcrootn_xfer_to_litter_patch)

    # Struct-first kernel: a device-view bundle carries no KA backend, so launch
    # manually off a known device output array + synchronize (mirrors methane.jl).
    ndrange = length(mask_soilp)
    if ndrange != 0
        backend = _kernel_backend(out.c_to_litr_c)
        _gap_p2c_kernel!(backend)(out, prof, fc, fn,
            mask_soilp, patch.itype, patch.column, patch.wtcol,
            nlevdecomp, i_litr_min, i_litr_max, i_met_lit; ndrange = ndrange)
        KA.synchronize(backend)
    end

    return nothing
end
