# ==========================================================================
# Ported from: src/biogeochem/CNFireBaseMod.F90
# Fire dynamics module for coupled carbon-nitrogen code (CN).
# Based on Li et al. (2012a,b; 2013; 2014).
#
# Public types:
#   CNFireConstData   — Fire constants shared by Li versions
#   CNFireParams      — Fire parameters read from params file
#   PftConFireBase    — PFT-level parameters for fire
#   CNFireBaseData    — Fire base data (btran2_patch)
#   DgvsFireData      — DGVS data needed by fire
#
# Public functions:
#   cnfire_calc_fire_root_wetness_li2014!  — Root wetness (original version)
#   cnfire_calc_fire_root_wetness_li2021!  — Root wetness (2021 version)
#   cnfire_fluxes!                         — Fire C/N flux calculations
# ==========================================================================

# ---------------------------------------------------------------------------
# KernelAbstractions kernels for fire base per-element loops.
# @kernel / @index / @Const / _launch! are module-wide (infrastructure/kernels.jl).
# Each kernel below replaces a simple, fully-independent per-element loop.
# `cmin`/`cmax` reproduce the original `for p in bounds` / `for c in bounds`
# iteration range while ndrange spans the full mask length (one thread per index).
# ---------------------------------------------------------------------------

# Zero btran2[p] for every masked patch in [cmin, cmax] (used for both the
# exposed and non-exposed veg zeroing loops in the Li2014/Li2021 root-wetness
# routines — same operation, called per mask).
@kernel function _fireb_zero_btran2_kernel!(btran2, @Const(mask), cmin::Int, cmax::Int)
    p = @index(Global)
    T = eltype(btran2)
    @inbounds if cmin <= p <= cmax && mask[p]
        btran2[p] = zero(T)
    end
end

fireb_zero_btran2!(btran2, mask, cmin::Int, cmax::Int) =
    _launch!(_fireb_zero_btran2_kernel!, btran2, mask, cmin, cmax)

# Clamp btran2[p] to <= 1.0 for every masked patch in [cmin, cmax].
@kernel function _fireb_clamp_btran2_kernel!(btran2, @Const(mask), cmin::Int, cmax::Int)
    p = @index(Global)
    T = eltype(btran2)
    @inbounds if cmin <= p <= cmax && mask[p]
        if btran2[p] > one(T)
            btran2[p] = one(T)
        end
    end
end

fireb_clamp_btran2!(btran2, mask, cmin::Int, cmax::Int) =
    _launch!(_fireb_clamp_btran2_kernel!, btran2, mask, cmin, cmax)

# Carbon loss to peat fires (per-column, fully independent). Branches on
# gridcell latitude vs. borealat; no accumulation.
@kernel function _fireb_peatfire_somc_kernel!(somc_fire, @Const(mask), @Const(gridcell),
                                              @Const(latdeg), @Const(totsomc),
                                              @Const(baf_peatf), cmin::Int, cmax::Int,
                                              borealat)
    c = @index(Global)
    T = eltype(somc_fire)
    @inbounds if cmin <= c <= cmax && mask[c]
        g = gridcell[c]
        if latdeg[g] < borealat
            somc_fire[c] = totsomc[c] * baf_peatf[c] * T(6) / T(33.9)
        else
            somc_fire[c] = baf_peatf[c] * T(2.2e3)
        end
    end
end

fireb_peatfire_somc!(somc_fire, mask, gridcell, latdeg, totsomc, baf_peatf,
                     cmin::Int, cmax::Int, borealat) =
    _launch!(_fireb_peatfire_somc_kernel!, somc_fire, mask, gridcell, latdeg,
             totsomc, baf_peatf, cmin, cmax, borealat)

# ---------------------------------------------------------------------------
# Constants for fire (from cnfire_const_type)
# ---------------------------------------------------------------------------

"""
    CNFireConstData

Fire constants shared by Li versions of the fire model.
Ported from `cnfire_const_type` in `CNFireBaseMod.F90`.
"""
Base.@kwdef mutable struct CNFireConstData
    borealat                    ::Float64 = 40.0       # latitude for boreal peat fires
    lfuel                       ::Float64 = 75.0       # lower threshold fuel mass (gC/m2)
    ufuel                       ::Float64 = 650.0      # upper threshold fuel mass (gC/m2)
    g0                          ::Float64 = 0.05       # g(W) when W=0 m/s
    rh_low                      ::Float64 = 30.0       # relative humidity low (%)
    rh_hgh                      ::Float64 = 80.0       # relative humidity high (%)
    bt_min                      ::Float64 = 0.3        # btran minimum (fraction)
    bt_max                      ::Float64 = 0.7        # btran maximum (fraction)
    cli_scale                   ::Float64 = 0.035      # global constant for deforestation fires (/d)
    boreal_peatfire_c           ::Float64 = 4.2e-5     # c parameter for boreal peatland fire (/hr)
    pot_hmn_ign_counts_alpha    ::Float64 = 0.0035     # potential human ignition counts (/person/month)
    non_boreal_peatfire_c       ::Float64 = 0.001      # c parameter for non-boreal peatland fire (/hr)
    cropfire_a1                 ::Float64 = 0.3        # a1 parameter for cropland fire (/hr)
    occur_hi_gdp_tree           ::Float64 = 0.39       # fire occurrence for high GDP tree-dominated areas
    cmb_cmplt_fact_litter       ::Float64 = 0.5        # combustion completion factor for litter
    cmb_cmplt_fact_cwd          ::Float64 = 0.25       # combustion completion factor for CWD
    max_rh30_affecting_fuel     ::Float64 = 90.0       # max 30-day RH affecting fuel (%)
    defo_fire_precip_thresh_bet ::Float64 = 4.0        # max running mean precip for BET deforestation fire (mm/d)
    defo_fire_precip_thresh_bdt ::Float64 = 1.8        # max running mean precip for BDT deforestation fire (mm/d)
    borpeat_fire_soilmoist_denom::Float64 = 0.3        # denominator in boreal peat fire soil moisture term
    nonborpeat_fire_precip_denom::Float64 = 1.0        # denominator in non-boreal peat fire precip term
end

# ---------------------------------------------------------------------------
# Parameters read from params file (from params_type)
# ---------------------------------------------------------------------------

"""
    CNFireParams

Fire parameters read from the input parameter file.
Ported from `params_type` in `CNFireBaseMod.F90`.
"""
Base.@kwdef mutable struct CNFireParams
    prh30               ::Float64 = 0.0    # factor for fuel combustibility dependence on 30-day RH
    ignition_efficiency ::Float64 = 0.0    # ignition efficiency of cloud-to-ground lightning
end

# ---------------------------------------------------------------------------
# PFT constants needed by fire base (from pftcon)
# ---------------------------------------------------------------------------

"""
    PftConFireBase

PFT-level parameters referenced by the fire base module.
Contains the subset of `pftconMod` fields used in `CNFireBaseMod.F90`.
"""
Base.@kwdef mutable struct PftConFireBase{FT<:Real,
                                          V<:AbstractVector{FT},
                                          M<:AbstractMatrix{FT}}
    woody    ::V = Float64[]  # binary woody flag (1=woody, 0=not woody)
    cc_leaf  ::V = Float64[]  # combustion completeness factor for leaves
    cc_lstem ::V = Float64[]  # combustion completeness factor for live stems
    cc_dstem ::V = Float64[]  # combustion completeness factor for dead stems
    cc_other ::V = Float64[]  # combustion completeness factor for other pools
    fm_leaf  ::V = Float64[]  # fire mortality factor for leaves
    fm_lstem ::V = Float64[]  # fire mortality factor for live stems
    fm_other ::V = Float64[]  # fire mortality factor for other pools
    fm_root  ::V = Float64[]  # fire mortality factor for fine roots
    fm_lroot ::V = Float64[]  # fire mortality factor for live coarse roots
    fm_droot ::V = Float64[]  # fire mortality factor for dead coarse roots
    lf_f     ::M = Matrix{Float64}(undef, 0, 0)  # leaf litter fractions (pft, litr)
    fr_f     ::M = Matrix{Float64}(undef, 0, 0)  # fine root litter fractions (pft, litr)
    smpso    ::V = Float64[]  # soil water potential at full stomatal opening (mm)
    smpsc    ::V = Float64[]  # soil water potential at full stomatal closure (mm)
end
PftConFireBase{FT}(; kwargs...) where {FT<:Real} =
    PftConFireBase{FT, Vector{FT}, Matrix{FT}}(; kwargs...)
Adapt.@adapt_structure PftConFireBase

# ---------------------------------------------------------------------------
# Fire base data (btran2_patch)
# ---------------------------------------------------------------------------

"""
    CNFireBaseData

Fire base data containing root zone soil wetness factor.
Ported from `cnfire_base_type` member data in `CNFireBaseMod.F90`.
"""
Base.@kwdef mutable struct CNFireBaseData{FT<:Real, V<:AbstractVector{FT}}
    btran2_patch ::V = Float64[]  # patch root zone soil wetness factor (0 to 1)
end
CNFireBaseData{FT}(; kwargs...) where {FT<:Real} = CNFireBaseData{FT, Vector{FT}}(; kwargs...)
Adapt.@adapt_structure CNFireBaseData

# ---------------------------------------------------------------------------
# DGVS data needed by fire
# ---------------------------------------------------------------------------

"""
    DgvsFireData

DGVS fields needed by fire routines (nind_patch, leafcmax_patch).
"""
Base.@kwdef mutable struct DgvsFireData{FT<:Real, V<:AbstractVector{FT}}
    nind_patch ::V = Float64[]  # number of individuals (#/m2)
end
DgvsFireData{FT}(; kwargs...) where {FT<:Real} = DgvsFireData{FT, Vector{FT}}(; kwargs...)
Adapt.@adapt_structure DgvsFireData

# ---------------------------------------------------------------------------
# cnfire_calc_fire_root_wetness_li2014! — Li 2014 root wetness
# ---------------------------------------------------------------------------

"""
    cnfire_calc_fire_root_wetness_li2014!(fire_data, mask_exposedveg, mask_noexposedveg,
        bounds, pftcon, patch, soilstate, h2osoi_vol_col, nlevgrnd;
        soil_suction_fn)

Calculate root zone soil wetness for fire model (Li et al. 2014 version).
Uses soil water retention curve to compute soil suction potential.

Ported from `CNFire_calc_fire_root_wetness_Li2014` in `CNFireBaseMod.F90`.
"""
function cnfire_calc_fire_root_wetness_li2014!(
    fire_data::CNFireBaseData,
    mask_exposedveg::BitVector,
    mask_noexposedveg::BitVector,
    bounds::UnitRange{Int},
    pftcon::PftConFireBase,
    patch::PatchData,
    soilstate::SoilStateData,
    h2osoi_vol_col::Matrix{<:Real},
    nlevgrnd::Int;
    soil_suction_fn::Function = default_soil_suction
)
    btran2   = fire_data.btran2_patch
    smpso    = pftcon.smpso
    smpsc    = pftcon.smpsc
    watsat   = soilstate.watsat_col
    rootfr   = soilstate.rootfr_patch
    sucsat   = soilstate.sucsat_col
    bsw      = soilstate.bsw_col
    column   = patch.column
    itype    = patch.itype

    cmin = first(bounds)
    cmax = last(bounds)

    # Zero out non-exposed vegetation patches
    fireb_zero_btran2!(btran2, mask_noexposedveg, cmin, cmax)

    # Initialize exposed veg patches to zero
    fireb_zero_btran2!(btran2, mask_exposedveg, cmin, cmax)

    # Accumulate root-weighted soil wetness across layers.
    # NOTE: the default soil-suction callback (default_soil_suction,
    # Clapp-Hornberger smp = -sucsat * s^(-bsw)) is INLINED into the kernel —
    # KA kernels can't call host callbacks. The test exercises the default.
    fireb_accum_btran2_li2014!(btran2, mask_exposedveg, cmin, cmax, nlevgrnd,
                               column, itype, h2osoi_vol_col, watsat,
                               rootfr, sucsat, bsw, smpso, smpsc)

    # Clamp to [0, 1]
    fireb_clamp_btran2!(btran2, mask_exposedveg, cmin, cmax)

    return nothing
end

"""
    default_soil_suction(c, j, s_node, soilstate)

Default soil suction function placeholder. Returns suction based on
Clapp-Hornberger parameterization. In practice, this should be replaced
with the appropriate soil water retention curve method.
"""
function default_soil_suction(c::Int, j::Int, s_node::Real, soilstate::SoilStateData)
    # Clapp-Hornberger: smp = sucsat * s^(-bsw)
    sucsat = soilstate.sucsat_col[c, j]
    bsw    = soilstate.bsw_col[c, j]
    return -sucsat * s_node^(-bsw)
end

# Per-patch accumulation of root-weighted soil wetness for Li2014.
# One thread per patch p; the j-loop runs INSIDE the thread in ascending order
# so btran2[p] accumulates loop-carried exactly as the host did.
# default_soil_suction's body is inlined (Clapp-Hornberger).
@kernel function _fireb_accum_btran2_li2014_kernel!(btran2, @Const(mask),
                                                    cmin::Int, cmax::Int, nlevgrnd::Int,
                                                    @Const(column), @Const(itype),
                                                    @Const(h2osoi_vol_col), @Const(watsat),
                                                    @Const(rootfr), @Const(sucsat),
                                                    @Const(bsw), @Const(smpso),
                                                    @Const(smpsc))
    p = @index(Global)
    T = eltype(btran2)
    @inbounds if cmin <= p <= cmax && mask[p]
        c = column[p]
        ivt = itype[p] + 1  # 0-based Fortran → 1-based Julia
        acc = btran2[p]
        for j in 1:nlevgrnd
            s_node = max(h2osoi_vol_col[c, j] / watsat[c, j], T(0.01))
            # default_soil_suction inlined: smp = -sucsat * s^(-bsw)
            smp_node_lf = -sucsat[c, j] * s_node^(-bsw[c, j])
            smp_node_lf = max(smpsc[ivt], smp_node_lf)
            acc = acc + rootfr[p, j] * max(zero(T),
                min((smp_node_lf - smpsc[ivt]) / (smpso[ivt] - smpsc[ivt]), one(T)))
        end
        btran2[p] = acc
    end
end

fireb_accum_btran2_li2014!(btran2, mask, cmin::Int, cmax::Int, nlevgrnd::Int,
                           column, itype, h2osoi_vol_col, watsat,
                           rootfr, sucsat, bsw, smpso, smpsc) =
    _launch!(_fireb_accum_btran2_li2014_kernel!, btran2, mask, cmin, cmax, nlevgrnd,
             column, itype, h2osoi_vol_col, watsat, rootfr, sucsat, bsw, smpso, smpsc)

# ---------------------------------------------------------------------------
# cnfire_calc_fire_root_wetness_li2021! — Li 2021 root wetness
# ---------------------------------------------------------------------------

"""
    cnfire_calc_fire_root_wetness_li2021!(fire_data, mask_exposedveg, mask_noexposedveg,
        bounds, patch, soilstate, h2osoi_vol_col, nlevgrnd)

Calculate root zone soil wetness for fire model (Li et al. 2021 version).
Uses simple volumetric water content ratio instead of soil suction.

Ported from `CNFire_calc_fire_root_wetness_Li2021` in `CNFireBaseMod.F90`.
"""
# Per-patch accumulation of root-weighted s_node for Li2021.
# One thread per patch p; the j-loop runs INSIDE the thread in ascending order
# so btran2[p] accumulates loop-carried exactly as the host did.
@kernel function _fireb_accum_btran2_li2021_kernel!(btran2, @Const(mask),
                                                    cmin::Int, cmax::Int, nlevgrnd::Int,
                                                    @Const(column),
                                                    @Const(h2osoi_vol_col), @Const(watsat),
                                                    @Const(rootfr))
    p = @index(Global)
    T = eltype(btran2)
    @inbounds if cmin <= p <= cmax && mask[p]
        c = column[p]
        acc = btran2[p]
        for j in 1:nlevgrnd
            s_node = max(h2osoi_vol_col[c, j] / watsat[c, j], T(0.01))
            acc = acc + rootfr[p, j] * s_node
        end
        btran2[p] = acc
    end
end

fireb_accum_btran2_li2021!(btran2, mask, cmin::Int, cmax::Int, nlevgrnd::Int,
                           column, h2osoi_vol_col, watsat, rootfr) =
    _launch!(_fireb_accum_btran2_li2021_kernel!, btran2, mask, cmin, cmax, nlevgrnd,
             column, h2osoi_vol_col, watsat, rootfr)

function cnfire_calc_fire_root_wetness_li2021!(
    fire_data::CNFireBaseData,
    mask_exposedveg::BitVector,
    mask_noexposedveg::BitVector,
    bounds::UnitRange{Int},
    patch::PatchData,
    soilstate::SoilStateData,
    h2osoi_vol_col::Matrix{<:Real},
    nlevgrnd::Int
)
    btran2 = fire_data.btran2_patch
    watsat = soilstate.watsat_col
    rootfr = soilstate.rootfr_patch
    column = patch.column

    cmin = first(bounds)
    cmax = last(bounds)

    # Zero out non-exposed vegetation patches
    fireb_zero_btran2!(btran2, mask_noexposedveg, cmin, cmax)

    # Initialize exposed veg patches to zero
    fireb_zero_btran2!(btran2, mask_exposedveg, cmin, cmax)

    # Accumulate root-weighted s_node across layers
    fireb_accum_btran2_li2021!(btran2, mask_exposedveg, cmin, cmax, nlevgrnd,
                               column, h2osoi_vol_col, watsat, rootfr)

    # Clamp to [0, 1]
    fireb_clamp_btran2!(btran2, mask_exposedveg, cmin, cmax)

    return nothing
end

# ---------------------------------------------------------------------------
# cnfire_fluxes! — Main fire flux calculation
# ---------------------------------------------------------------------------

# KernelAbstractions kernels for cnfire_fluxes! loop sections.
# Each kernel below reproduces one host loop section byte-identically.
# ---------------------------------------------------------------------------

# Section (1): per-patch fire emission + mortality fluxes. Each masked patch
# writes only its own index [p] (no accumulation). `mask_actfirep[p]` is set
# per-patch. The `!use_matrixcn` block and the `use_cndv` block are gated by
# host-resolved Bool flags. `f` is recomputed inside the thread exactly as the
# host loop did (transient_landcover selects fbac vs farea_burned).
"""
Device-view bundles for the per-patch fire C/N flux kernel.

All per-patch float arrays share one element type -> parameter V.
mask_actfirep is an output Bool/Bit vector -> parameter B (kept in the outputs
bundle, written per-patch). Int index arrays (column, itype) stay LOOSE in the
kernel signature. Scalars (m_spinup, dt) are bundled into _FirePatchScalars{T}.

Every field name is the Fortran/original local var name; the kernel aliases each
struct field back to that local at the top so the body stays verbatim.
"""

# --- C "*_to_fire" outputs + gresp-to-fire (19 fields) ---
Base.@kwdef struct _FireCFireOut{V}
    m_gresp_storage_to_fire::V; m_gresp_xfer_to_fire::V
    m_leafc_to_fire::V; m_leafc_storage_to_fire::V; m_leafc_xfer_to_fire::V
    m_livestemc_to_fire::V; m_livestemc_storage_to_fire::V; m_livestemc_xfer_to_fire::V
    m_deadstemc_to_fire::V; m_deadstemc_storage_to_fire::V; m_deadstemc_xfer_to_fire::V
    m_frootc_to_fire::V; m_frootc_storage_to_fire::V; m_frootc_xfer_to_fire::V
    m_livecrootc_to_fire::V; m_livecrootc_storage_to_fire::V; m_livecrootc_xfer_to_fire::V
    m_deadcrootc_to_fire::V; m_deadcrootc_storage_to_fire::V; m_deadcrootc_xfer_to_fire::V
end
Adapt.@adapt_structure _FireCFireOut

# --- N "*_to_fire" outputs (16 fields) ---
Base.@kwdef struct _FireNFireOut{V}
    m_leafn_to_fire::V; m_leafn_storage_to_fire::V; m_leafn_xfer_to_fire::V
    m_livestemn_to_fire::V; m_livestemn_storage_to_fire::V; m_livestemn_xfer_to_fire::V
    m_deadstemn_to_fire::V; m_deadstemn_storage_to_fire::V; m_deadstemn_xfer_to_fire::V
    m_frootn_to_fire::V; m_frootn_storage_to_fire::V; m_frootn_xfer_to_fire::V
    m_livecrootn_to_fire::V; m_livecrootn_storage_to_fire::V; m_livecrootn_xfer_to_fire::V
    m_deadcrootn_to_fire::V; m_deadcrootn_storage_to_fire::V; m_deadcrootn_xfer_to_fire::V
    m_retransn_to_fire::V
end
Adapt.@adapt_structure _FireNFireOut

# --- C "*_to_litter_fire" outputs (24 fields) ---
Base.@kwdef struct _FireCLitOut{V}
    m_leafc_to_litter_fire::V; m_leafc_storage_to_litter_fire::V; m_leafc_xfer_to_litter_fire::V
    m_livestemc_to_litter_fire::V; m_livestemc_storage_to_litter_fire::V; m_livestemc_xfer_to_litter_fire::V
    m_livestemc_to_deadstemc_fire::V; m_deadstemc_to_litter_fire::V
    m_deadstemc_storage_to_litter_fire::V; m_deadstemc_xfer_to_litter_fire::V
    m_frootc_to_litter_fire::V; m_frootc_storage_to_litter_fire::V; m_frootc_xfer_to_litter_fire::V
    m_livecrootc_to_litter_fire::V; m_livecrootc_storage_to_litter_fire::V; m_livecrootc_xfer_to_litter_fire::V
    m_livecrootc_to_deadcrootc_fire::V; m_deadcrootc_to_litter_fire::V
    m_deadcrootc_storage_to_litter_fire::V; m_deadcrootc_xfer_to_litter_fire::V
    m_gresp_storage_to_litter_fire::V; m_gresp_xfer_to_litter_fire::V
end
Adapt.@adapt_structure _FireCLitOut

# --- N "*_to_litter_fire" outputs (22 fields) ---
Base.@kwdef struct _FireNLitOut{V}
    m_leafn_to_litter_fire::V; m_leafn_storage_to_litter_fire::V; m_leafn_xfer_to_litter_fire::V
    m_livestemn_to_litter_fire::V; m_livestemn_storage_to_litter_fire::V; m_livestemn_xfer_to_litter_fire::V
    m_livestemn_to_deadstemn_fire::V; m_deadstemn_to_litter_fire::V
    m_deadstemn_storage_to_litter_fire::V; m_deadstemn_xfer_to_litter_fire::V
    m_frootn_to_litter_fire::V; m_frootn_storage_to_litter_fire::V; m_frootn_xfer_to_litter_fire::V
    m_livecrootn_to_litter_fire::V; m_livecrootn_storage_to_litter_fire::V; m_livecrootn_xfer_to_litter_fire::V
    m_livecrootn_to_deadcrootn_fire::V; m_deadcrootn_to_litter_fire::V
    m_deadcrootn_storage_to_litter_fire::V; m_deadcrootn_xfer_to_litter_fire::V
    m_retransn_to_litter_fire::V
end
Adapt.@adapt_structure _FireNLitOut

# --- CNDV outputs (read+write) + the actfirep mask (Bool/Bit output) ---
Base.@kwdef struct _FireCNDVOut{V,B}
    nind::V; leafcmax::V
    mask_actfirep::B
end
Adapt.@adapt_structure _FireCNDVOut

# --- read-only forcing: column-indexed scalars + per-patch C state (20 fields) ---
Base.@kwdef struct _FireStateC{V}
    cropf_col::V; farea_burned::V; fbac::V; baf_crop::V
    gresp_storage::V; gresp_xfer::V
    leafc::V; leafc_storage::V; leafc_xfer::V
    livestemc::V; livestemc_storage::V; livestemc_xfer::V
    deadstemc::V; deadstemc_storage::V; deadstemc_xfer::V
    frootc::V; frootc_storage::V; frootc_xfer::V
    livecrootc::V; livecrootc_storage::V; livecrootc_xfer::V
    deadcrootc::V; deadcrootc_storage::V; deadcrootc_xfer::V
end
Adapt.@adapt_structure _FireStateC

# --- read-only per-patch N state (19 fields) ---
Base.@kwdef struct _FireStateN{V}
    leafn::V; leafn_storage::V; leafn_xfer::V
    livestemn::V; livestemn_storage::V; livestemn_xfer::V
    deadstemn::V; deadstemn_storage::V; deadstemn_xfer::V
    frootn::V; frootn_storage::V; frootn_xfer::V
    livecrootn::V; livecrootn_storage::V; livecrootn_xfer::V
    deadcrootn::V; deadcrootn_storage::V; deadcrootn_xfer::V
    retransn::V
end
Adapt.@adapt_structure _FireStateN

# --- read-only pftcon (ivt-indexed) combustion / mortality params (11 fields) ---
Base.@kwdef struct _FirePft{V}
    woody::V
    cc_leaf::V; cc_lstem::V; cc_dstem::V; cc_other::V
    fm_leaf::V; fm_lstem::V; fm_other::V; fm_root::V; fm_lroot::V; fm_droot::V
end
Adapt.@adapt_structure _FirePft

# --- isbits scalar bundle (working precision) ---
Base.@kwdef struct _FirePatchScalars{T}
    m_spinup::T
    dt::T
end

@kernel function _fireb_cnfire_patchflux_kernel!(
    cfo::_FireCFireOut, nfo::_FireNFireOut,
    clit::_FireCLitOut, nlit::_FireNLitOut, cndv::_FireCNDVOut,
    @Const(mask_soilp), @Const(itype), @Const(column),
    stc::_FireStateC, stn::_FireStateN, pft::_FirePft,
    nc3crop::Int, noveg::Int, prm::_FirePatchScalars,
    transient_landcover::Bool, use_matrixcn::Bool, use_cndv::Bool
)
    p = @index(Global)
    @inbounds if mask_soilp[p]
        T = eltype(cfo.m_gresp_storage_to_fire)

        # --- alias output bundles back to Fortran local names ---
        m_gresp_storage_to_fire = cfo.m_gresp_storage_to_fire
        m_gresp_xfer_to_fire = cfo.m_gresp_xfer_to_fire
        m_leafc_to_fire = cfo.m_leafc_to_fire
        m_leafc_storage_to_fire = cfo.m_leafc_storage_to_fire
        m_leafc_xfer_to_fire = cfo.m_leafc_xfer_to_fire
        m_livestemc_to_fire = cfo.m_livestemc_to_fire
        m_livestemc_storage_to_fire = cfo.m_livestemc_storage_to_fire
        m_livestemc_xfer_to_fire = cfo.m_livestemc_xfer_to_fire
        m_deadstemc_to_fire = cfo.m_deadstemc_to_fire
        m_deadstemc_storage_to_fire = cfo.m_deadstemc_storage_to_fire
        m_deadstemc_xfer_to_fire = cfo.m_deadstemc_xfer_to_fire
        m_frootc_to_fire = cfo.m_frootc_to_fire
        m_frootc_storage_to_fire = cfo.m_frootc_storage_to_fire
        m_frootc_xfer_to_fire = cfo.m_frootc_xfer_to_fire
        m_livecrootc_to_fire = cfo.m_livecrootc_to_fire
        m_livecrootc_storage_to_fire = cfo.m_livecrootc_storage_to_fire
        m_livecrootc_xfer_to_fire = cfo.m_livecrootc_xfer_to_fire
        m_deadcrootc_to_fire = cfo.m_deadcrootc_to_fire
        m_deadcrootc_storage_to_fire = cfo.m_deadcrootc_storage_to_fire
        m_deadcrootc_xfer_to_fire = cfo.m_deadcrootc_xfer_to_fire

        m_leafn_to_fire = nfo.m_leafn_to_fire
        m_leafn_storage_to_fire = nfo.m_leafn_storage_to_fire
        m_leafn_xfer_to_fire = nfo.m_leafn_xfer_to_fire
        m_livestemn_to_fire = nfo.m_livestemn_to_fire
        m_livestemn_storage_to_fire = nfo.m_livestemn_storage_to_fire
        m_livestemn_xfer_to_fire = nfo.m_livestemn_xfer_to_fire
        m_deadstemn_to_fire = nfo.m_deadstemn_to_fire
        m_deadstemn_storage_to_fire = nfo.m_deadstemn_storage_to_fire
        m_deadstemn_xfer_to_fire = nfo.m_deadstemn_xfer_to_fire
        m_frootn_to_fire = nfo.m_frootn_to_fire
        m_frootn_storage_to_fire = nfo.m_frootn_storage_to_fire
        m_frootn_xfer_to_fire = nfo.m_frootn_xfer_to_fire
        m_livecrootn_to_fire = nfo.m_livecrootn_to_fire
        m_livecrootn_storage_to_fire = nfo.m_livecrootn_storage_to_fire
        m_livecrootn_xfer_to_fire = nfo.m_livecrootn_xfer_to_fire
        m_deadcrootn_to_fire = nfo.m_deadcrootn_to_fire
        m_deadcrootn_storage_to_fire = nfo.m_deadcrootn_storage_to_fire
        m_deadcrootn_xfer_to_fire = nfo.m_deadcrootn_xfer_to_fire
        m_retransn_to_fire = nfo.m_retransn_to_fire

        m_leafc_to_litter_fire = clit.m_leafc_to_litter_fire
        m_leafc_storage_to_litter_fire = clit.m_leafc_storage_to_litter_fire
        m_leafc_xfer_to_litter_fire = clit.m_leafc_xfer_to_litter_fire
        m_livestemc_to_litter_fire = clit.m_livestemc_to_litter_fire
        m_livestemc_storage_to_litter_fire = clit.m_livestemc_storage_to_litter_fire
        m_livestemc_xfer_to_litter_fire = clit.m_livestemc_xfer_to_litter_fire
        m_livestemc_to_deadstemc_fire = clit.m_livestemc_to_deadstemc_fire
        m_deadstemc_to_litter_fire = clit.m_deadstemc_to_litter_fire
        m_deadstemc_storage_to_litter_fire = clit.m_deadstemc_storage_to_litter_fire
        m_deadstemc_xfer_to_litter_fire = clit.m_deadstemc_xfer_to_litter_fire
        m_frootc_to_litter_fire = clit.m_frootc_to_litter_fire
        m_frootc_storage_to_litter_fire = clit.m_frootc_storage_to_litter_fire
        m_frootc_xfer_to_litter_fire = clit.m_frootc_xfer_to_litter_fire
        m_livecrootc_to_litter_fire = clit.m_livecrootc_to_litter_fire
        m_livecrootc_storage_to_litter_fire = clit.m_livecrootc_storage_to_litter_fire
        m_livecrootc_xfer_to_litter_fire = clit.m_livecrootc_xfer_to_litter_fire
        m_livecrootc_to_deadcrootc_fire = clit.m_livecrootc_to_deadcrootc_fire
        m_deadcrootc_to_litter_fire = clit.m_deadcrootc_to_litter_fire
        m_deadcrootc_storage_to_litter_fire = clit.m_deadcrootc_storage_to_litter_fire
        m_deadcrootc_xfer_to_litter_fire = clit.m_deadcrootc_xfer_to_litter_fire
        m_gresp_storage_to_litter_fire = clit.m_gresp_storage_to_litter_fire
        m_gresp_xfer_to_litter_fire = clit.m_gresp_xfer_to_litter_fire

        m_leafn_to_litter_fire = nlit.m_leafn_to_litter_fire
        m_leafn_storage_to_litter_fire = nlit.m_leafn_storage_to_litter_fire
        m_leafn_xfer_to_litter_fire = nlit.m_leafn_xfer_to_litter_fire
        m_livestemn_to_litter_fire = nlit.m_livestemn_to_litter_fire
        m_livestemn_storage_to_litter_fire = nlit.m_livestemn_storage_to_litter_fire
        m_livestemn_xfer_to_litter_fire = nlit.m_livestemn_xfer_to_litter_fire
        m_livestemn_to_deadstemn_fire = nlit.m_livestemn_to_deadstemn_fire
        m_deadstemn_to_litter_fire = nlit.m_deadstemn_to_litter_fire
        m_deadstemn_storage_to_litter_fire = nlit.m_deadstemn_storage_to_litter_fire
        m_deadstemn_xfer_to_litter_fire = nlit.m_deadstemn_xfer_to_litter_fire
        m_frootn_to_litter_fire = nlit.m_frootn_to_litter_fire
        m_frootn_storage_to_litter_fire = nlit.m_frootn_storage_to_litter_fire
        m_frootn_xfer_to_litter_fire = nlit.m_frootn_xfer_to_litter_fire
        m_livecrootn_to_litter_fire = nlit.m_livecrootn_to_litter_fire
        m_livecrootn_storage_to_litter_fire = nlit.m_livecrootn_storage_to_litter_fire
        m_livecrootn_xfer_to_litter_fire = nlit.m_livecrootn_xfer_to_litter_fire
        m_livecrootn_to_deadcrootn_fire = nlit.m_livecrootn_to_deadcrootn_fire
        m_deadcrootn_to_litter_fire = nlit.m_deadcrootn_to_litter_fire
        m_deadcrootn_storage_to_litter_fire = nlit.m_deadcrootn_storage_to_litter_fire
        m_deadcrootn_xfer_to_litter_fire = nlit.m_deadcrootn_xfer_to_litter_fire
        m_retransn_to_litter_fire = nlit.m_retransn_to_litter_fire

        nind = cndv.nind
        leafcmax = cndv.leafcmax
        mask_actfirep = cndv.mask_actfirep

        # --- alias read-only state / pftcon ---
        cropf_col = stc.cropf_col; farea_burned = stc.farea_burned
        fbac = stc.fbac; baf_crop = stc.baf_crop
        gresp_storage = stc.gresp_storage; gresp_xfer = stc.gresp_xfer
        leafc = stc.leafc; leafc_storage = stc.leafc_storage; leafc_xfer = stc.leafc_xfer
        livestemc = stc.livestemc; livestemc_storage = stc.livestemc_storage; livestemc_xfer = stc.livestemc_xfer
        deadstemc = stc.deadstemc; deadstemc_storage = stc.deadstemc_storage; deadstemc_xfer = stc.deadstemc_xfer
        frootc = stc.frootc; frootc_storage = stc.frootc_storage; frootc_xfer = stc.frootc_xfer
        livecrootc = stc.livecrootc; livecrootc_storage = stc.livecrootc_storage; livecrootc_xfer = stc.livecrootc_xfer
        deadcrootc = stc.deadcrootc; deadcrootc_storage = stc.deadcrootc_storage; deadcrootc_xfer = stc.deadcrootc_xfer

        leafn = stn.leafn; leafn_storage = stn.leafn_storage; leafn_xfer = stn.leafn_xfer
        livestemn = stn.livestemn; livestemn_storage = stn.livestemn_storage; livestemn_xfer = stn.livestemn_xfer
        deadstemn = stn.deadstemn; deadstemn_storage = stn.deadstemn_storage; deadstemn_xfer = stn.deadstemn_xfer
        frootn = stn.frootn; frootn_storage = stn.frootn_storage; frootn_xfer = stn.frootn_xfer
        livecrootn = stn.livecrootn; livecrootn_storage = stn.livecrootn_storage; livecrootn_xfer = stn.livecrootn_xfer
        deadcrootn = stn.deadcrootn; deadcrootn_storage = stn.deadcrootn_storage; deadcrootn_xfer = stn.deadcrootn_xfer
        retransn = stn.retransn

        woody = pft.woody
        cc_leaf = pft.cc_leaf; cc_lstem = pft.cc_lstem; cc_dstem = pft.cc_dstem; cc_other = pft.cc_other
        fm_leaf = pft.fm_leaf; fm_lstem = pft.fm_lstem; fm_other = pft.fm_other
        fm_root = pft.fm_root; fm_lroot = pft.fm_lroot; fm_droot = pft.fm_droot

        c = column[p]
        m = prm.m_spinup
        dt = prm.dt

        if itype[p] < nc3crop && cropf_col[c] < one(T)
            if transient_landcover
                f = (fbac[c] - baf_crop[c]) / (one(T) - cropf_col[c])
            else
                f = (farea_burned[c] - baf_crop[c]) / (one(T) - cropf_col[c])
            end
        else
            if cropf_col[c] > zero(T)
                f = baf_crop[c] / cropf_col[c]
            else
                f = zero(T)
            end
        end

        if f != zero(T)
            mask_actfirep[p] = true
        end

        ivt = itype[p] + 1

        m_gresp_storage_to_fire[p]       = gresp_storage[p]      * f * cc_other[ivt]
        m_gresp_xfer_to_fire[p]          = gresp_xfer[p]         * f * cc_other[ivt]

        if !use_matrixcn
            m_leafc_to_fire[p]               = leafc[p]              * f * cc_leaf[ivt]
            m_leafc_storage_to_fire[p]       = leafc_storage[p]      * f * cc_other[ivt]
            m_leafc_xfer_to_fire[p]          = leafc_xfer[p]         * f * cc_other[ivt]
            m_livestemc_to_fire[p]           = livestemc[p]          * f * cc_lstem[ivt]
            m_livestemc_storage_to_fire[p]   = livestemc_storage[p]  * f * cc_other[ivt]
            m_livestemc_xfer_to_fire[p]      = livestemc_xfer[p]     * f * cc_other[ivt]
            m_deadstemc_to_fire[p]           = deadstemc[p]          * f * cc_dstem[ivt] * m
            m_deadstemc_storage_to_fire[p]   = deadstemc_storage[p]  * f * cc_other[ivt]
            m_deadstemc_xfer_to_fire[p]      = deadstemc_xfer[p]     * f * cc_other[ivt]
            m_frootc_to_fire[p]              = frootc[p]             * f * zero(T)
            m_frootc_storage_to_fire[p]      = frootc_storage[p]     * f * cc_other[ivt]
            m_frootc_xfer_to_fire[p]         = frootc_xfer[p]        * f * cc_other[ivt]
            m_livecrootc_to_fire[p]          = livecrootc[p]         * f * zero(T)
            m_livecrootc_storage_to_fire[p]  = livecrootc_storage[p] * f * cc_other[ivt]
            m_livecrootc_xfer_to_fire[p]     = livecrootc_xfer[p]    * f * cc_other[ivt]
            m_deadcrootc_to_fire[p]          = deadcrootc[p]         * f * zero(T)
            m_deadcrootc_storage_to_fire[p]  = deadcrootc_storage[p] * f * cc_other[ivt]
            m_deadcrootc_xfer_to_fire[p]     = deadcrootc_xfer[p]    * f * cc_other[ivt]

            m_leafn_to_fire[p]               = leafn[p]              * f * cc_leaf[ivt]
            m_leafn_storage_to_fire[p]       = leafn_storage[p]      * f * cc_other[ivt]
            m_leafn_xfer_to_fire[p]          = leafn_xfer[p]         * f * cc_other[ivt]
            m_livestemn_to_fire[p]           = livestemn[p]          * f * cc_lstem[ivt]
            m_livestemn_storage_to_fire[p]   = livestemn_storage[p]  * f * cc_other[ivt]
            m_livestemn_xfer_to_fire[p]      = livestemn_xfer[p]     * f * cc_other[ivt]
            m_deadstemn_to_fire[p]           = deadstemn[p]          * f * cc_dstem[ivt] * m
            m_deadstemn_storage_to_fire[p]   = deadstemn_storage[p]  * f * cc_other[ivt]
            m_deadstemn_xfer_to_fire[p]      = deadstemn_xfer[p]     * f * cc_other[ivt]
            m_frootn_to_fire[p]              = frootn[p]             * f * zero(T)
            m_frootn_storage_to_fire[p]      = frootn_storage[p]     * f * cc_other[ivt]
            m_frootn_xfer_to_fire[p]         = frootn_xfer[p]        * f * cc_other[ivt]
            m_livecrootn_to_fire[p]          = livecrootn[p]         * f * zero(T)
            m_livecrootn_storage_to_fire[p]  = livecrootn_storage[p] * f * cc_other[ivt]
            m_livecrootn_xfer_to_fire[p]     = livecrootn_xfer[p]    * f * cc_other[ivt]
            m_deadcrootn_to_fire[p]          = deadcrootn[p]         * f * zero(T)
            m_deadcrootn_xfer_to_fire[p]     = deadcrootn_xfer[p]    * f * cc_other[ivt]
            m_deadcrootn_storage_to_fire[p]  = deadcrootn_storage[p] * f * cc_other[ivt]
            m_retransn_to_fire[p]            = retransn[p]           * f * cc_other[ivt]
        end

        if !use_matrixcn
            m_leafc_to_litter_fire[p]                   = leafc[p] * f *
                (one(T) - cc_leaf[ivt]) * fm_leaf[ivt]
            m_leafc_storage_to_litter_fire[p]           = leafc_storage[p] * f *
                (one(T) - cc_other[ivt]) * fm_other[ivt]
            m_leafc_xfer_to_litter_fire[p]              = leafc_xfer[p] * f *
                (one(T) - cc_other[ivt]) * fm_other[ivt]
            m_livestemc_to_litter_fire[p]               = livestemc[p] * f *
                (one(T) - cc_lstem[ivt]) * fm_droot[ivt]
            m_livestemc_storage_to_litter_fire[p]       = livestemc_storage[p] * f *
                (one(T) - cc_other[ivt]) * fm_other[ivt]
            m_livestemc_xfer_to_litter_fire[p]          = livestemc_xfer[p] * f *
                (one(T) - cc_other[ivt]) * fm_other[ivt]
            m_livestemc_to_deadstemc_fire[p]            = livestemc[p] * f *
                (one(T) - cc_lstem[ivt]) * (fm_lstem[ivt] - fm_droot[ivt])
            m_deadstemc_to_litter_fire[p]               = deadstemc[p] * f * m *
                (one(T) - cc_dstem[ivt]) * fm_droot[ivt]
            m_deadstemc_storage_to_litter_fire[p]       = deadstemc_storage[p] * f *
                (one(T) - cc_other[ivt]) * fm_other[ivt]
            m_deadstemc_xfer_to_litter_fire[p]          = deadstemc_xfer[p] * f *
                (one(T) - cc_other[ivt]) * fm_other[ivt]
            m_frootc_to_litter_fire[p]                  = frootc[p] * f *
                fm_root[ivt]
            m_frootc_storage_to_litter_fire[p]          = frootc_storage[p] * f *
                (one(T) - cc_other[ivt]) * fm_other[ivt]
            m_frootc_xfer_to_litter_fire[p]             = frootc_xfer[p] * f *
                (one(T) - cc_other[ivt]) * fm_other[ivt]
            m_livecrootc_to_litter_fire[p]              = livecrootc[p] * f *
                fm_droot[ivt]
            m_livecrootc_storage_to_litter_fire[p]      = livecrootc_storage[p] * f *
                (one(T) - cc_other[ivt]) * fm_other[ivt]
            m_livecrootc_xfer_to_litter_fire[p]         = livecrootc_xfer[p] * f *
                (one(T) - cc_other[ivt]) * fm_other[ivt]
            m_livecrootc_to_deadcrootc_fire[p]          = livecrootc[p] * f *
                (fm_lroot[ivt] - fm_droot[ivt])
            m_deadcrootc_to_litter_fire[p]              = deadcrootc[p] * f * m *
                fm_droot[ivt]
            m_deadcrootc_storage_to_litter_fire[p]      = deadcrootc_storage[p] * f *
                (one(T) - cc_other[ivt]) * fm_other[ivt]
            m_deadcrootc_xfer_to_litter_fire[p]         = deadcrootc_xfer[p] * f *
                (one(T) - cc_other[ivt]) * fm_other[ivt]
            m_gresp_storage_to_litter_fire[p]           = gresp_storage[p] * f *
                (one(T) - cc_other[ivt]) * fm_other[ivt]
            m_gresp_xfer_to_litter_fire[p]              = gresp_xfer[p] * f *
                (one(T) - cc_other[ivt]) * fm_other[ivt]

            m_leafn_to_litter_fire[p]                  = leafn[p] * f *
                (one(T) - cc_leaf[ivt]) * fm_leaf[ivt]
            m_leafn_storage_to_litter_fire[p]          = leafn_storage[p] * f *
                (one(T) - cc_other[ivt]) * fm_other[ivt]
            m_leafn_xfer_to_litter_fire[p]             = leafn_xfer[p] * f *
                (one(T) - cc_other[ivt]) * fm_other[ivt]
            m_livestemn_to_litter_fire[p]              = livestemn[p] * f *
                (one(T) - cc_lstem[ivt]) * fm_droot[ivt]
            m_livestemn_storage_to_litter_fire[p]      = livestemn_storage[p] * f *
                (one(T) - cc_other[ivt]) * fm_other[ivt]
            m_livestemn_xfer_to_litter_fire[p]         = livestemn_xfer[p] * f *
                (one(T) - cc_other[ivt]) * fm_other[ivt]
            m_livestemn_to_deadstemn_fire[p]           = livestemn[p] * f *
                (one(T) - cc_lstem[ivt]) * (fm_lstem[ivt] - fm_droot[ivt])
            m_deadstemn_to_litter_fire[p]              = deadstemn[p] * f * m *
                (one(T) - cc_dstem[ivt]) * fm_droot[ivt]
            m_deadstemn_storage_to_litter_fire[p]      = deadstemn_storage[p] * f *
                (one(T) - cc_other[ivt]) * fm_other[ivt]
            m_deadstemn_xfer_to_litter_fire[p]         = deadstemn_xfer[p] * f *
                (one(T) - cc_other[ivt]) * fm_other[ivt]
            m_frootn_to_litter_fire[p]                 = frootn[p] * f *
                fm_root[ivt]
            m_frootn_storage_to_litter_fire[p]         = frootn_storage[p] * f *
                (one(T) - cc_other[ivt]) * fm_other[ivt]
            m_frootn_xfer_to_litter_fire[p]            = frootn_xfer[p] * f *
                (one(T) - cc_other[ivt]) * fm_other[ivt]
            m_livecrootn_to_litter_fire[p]             = livecrootn[p] * f *
                fm_droot[ivt]
            m_livecrootn_storage_to_litter_fire[p]     = livecrootn_storage[p] * f *
                (one(T) - cc_other[ivt]) * fm_other[ivt]
            m_livecrootn_xfer_to_litter_fire[p]        = livecrootn_xfer[p] * f *
                (one(T) - cc_other[ivt]) * fm_other[ivt]
            m_livecrootn_to_deadcrootn_fire[p]         = livecrootn[p] * f *
                (fm_lroot[ivt] - fm_droot[ivt])
            m_deadcrootn_to_litter_fire[p]             = deadcrootn[p] * f * m *
                fm_droot[ivt]
            m_deadcrootn_storage_to_litter_fire[p]     = deadcrootn_storage[p] * f *
                (one(T) - cc_other[ivt]) * fm_other[ivt]
            m_deadcrootn_xfer_to_litter_fire[p]        = deadcrootn_xfer[p] * f *
                (one(T) - cc_other[ivt]) * fm_other[ivt]
            m_retransn_to_litter_fire[p]               = retransn[p] * f *
                (one(T) - cc_other[ivt]) * fm_other[ivt]
        end

        if use_cndv
            if woody[ivt] == one(T)
                if livestemc[p] + deadstemc[p] > zero(T)
                    nind[p] = nind[p] * (one(T) - one(T) * fm_droot[ivt] * f)
                else
                    nind[p] = zero(T)
                end
            end
            leafcmax[p] = max(leafc[p] - m_leafc_to_fire[p] * dt, leafcmax[p])
            if ivt == noveg
                leafcmax[p] = zero(T)
            end
        end
    end
end

# Wrapper: same outer name + same positional arg order as before (callers in
# cnfire_fluxes! pass exactly these). Builds the device-view bundles + scalar
# bundle, then launches the struct-arg kernel via the manual struct-first form.
function fireb_cnfire_patchflux!(
    mask_actfirep,
    m_gresp_storage_to_fire, m_gresp_xfer_to_fire,
    m_leafc_to_fire, m_leafc_storage_to_fire, m_leafc_xfer_to_fire,
    m_livestemc_to_fire, m_livestemc_storage_to_fire, m_livestemc_xfer_to_fire,
    m_deadstemc_to_fire, m_deadstemc_storage_to_fire, m_deadstemc_xfer_to_fire,
    m_frootc_to_fire, m_frootc_storage_to_fire, m_frootc_xfer_to_fire,
    m_livecrootc_to_fire, m_livecrootc_storage_to_fire, m_livecrootc_xfer_to_fire,
    m_deadcrootc_to_fire, m_deadcrootc_storage_to_fire, m_deadcrootc_xfer_to_fire,
    m_leafn_to_fire, m_leafn_storage_to_fire, m_leafn_xfer_to_fire,
    m_livestemn_to_fire, m_livestemn_storage_to_fire, m_livestemn_xfer_to_fire,
    m_deadstemn_to_fire, m_deadstemn_storage_to_fire, m_deadstemn_xfer_to_fire,
    m_frootn_to_fire, m_frootn_storage_to_fire, m_frootn_xfer_to_fire,
    m_livecrootn_to_fire, m_livecrootn_storage_to_fire, m_livecrootn_xfer_to_fire,
    m_deadcrootn_to_fire, m_deadcrootn_storage_to_fire, m_deadcrootn_xfer_to_fire,
    m_retransn_to_fire,
    m_leafc_to_litter_fire, m_leafc_storage_to_litter_fire, m_leafc_xfer_to_litter_fire,
    m_livestemc_to_litter_fire, m_livestemc_storage_to_litter_fire, m_livestemc_xfer_to_litter_fire,
    m_livestemc_to_deadstemc_fire, m_deadstemc_to_litter_fire,
    m_deadstemc_storage_to_litter_fire, m_deadstemc_xfer_to_litter_fire,
    m_frootc_to_litter_fire, m_frootc_storage_to_litter_fire, m_frootc_xfer_to_litter_fire,
    m_livecrootc_to_litter_fire, m_livecrootc_storage_to_litter_fire, m_livecrootc_xfer_to_litter_fire,
    m_livecrootc_to_deadcrootc_fire, m_deadcrootc_to_litter_fire,
    m_deadcrootc_storage_to_litter_fire, m_deadcrootc_xfer_to_litter_fire,
    m_gresp_storage_to_litter_fire, m_gresp_xfer_to_litter_fire,
    m_leafn_to_litter_fire, m_leafn_storage_to_litter_fire, m_leafn_xfer_to_litter_fire,
    m_livestemn_to_litter_fire, m_livestemn_storage_to_litter_fire, m_livestemn_xfer_to_litter_fire,
    m_livestemn_to_deadstemn_fire, m_deadstemn_to_litter_fire,
    m_deadstemn_storage_to_litter_fire, m_deadstemn_xfer_to_litter_fire,
    m_frootn_to_litter_fire, m_frootn_storage_to_litter_fire, m_frootn_xfer_to_litter_fire,
    m_livecrootn_to_litter_fire, m_livecrootn_storage_to_litter_fire, m_livecrootn_xfer_to_litter_fire,
    m_livecrootn_to_deadcrootn_fire, m_deadcrootn_to_litter_fire,
    m_deadcrootn_storage_to_litter_fire, m_deadcrootn_xfer_to_litter_fire,
    m_retransn_to_litter_fire,
    nind, leafcmax,
    mask_soilp, itype, column,
    cropf_col, farea_burned, fbac, baf_crop,
    gresp_storage, gresp_xfer,
    leafc, leafc_storage, leafc_xfer,
    livestemc, livestemc_storage, livestemc_xfer,
    deadstemc, deadstemc_storage, deadstemc_xfer,
    frootc, frootc_storage, frootc_xfer,
    livecrootc, livecrootc_storage, livecrootc_xfer,
    deadcrootc, deadcrootc_storage, deadcrootc_xfer,
    leafn, leafn_storage, leafn_xfer,
    livestemn, livestemn_storage, livestemn_xfer,
    deadstemn, deadstemn_storage, deadstemn_xfer,
    frootn, frootn_storage, frootn_xfer,
    livecrootn, livecrootn_storage, livecrootn_xfer,
    deadcrootn, deadcrootn_storage, deadcrootn_xfer,
    retransn,
    woody, cc_leaf, cc_lstem, cc_dstem, cc_other,
    fm_leaf, fm_lstem, fm_other, fm_root, fm_lroot, fm_droot,
    nc3crop, noveg, m_spinup, dt,
    transient_landcover, use_matrixcn, use_cndv)

    n = length(mask_soilp)
    n == 0 && return mask_actfirep

    cfo = _FireCFireOut(;
        m_gresp_storage_to_fire, m_gresp_xfer_to_fire,
        m_leafc_to_fire, m_leafc_storage_to_fire, m_leafc_xfer_to_fire,
        m_livestemc_to_fire, m_livestemc_storage_to_fire, m_livestemc_xfer_to_fire,
        m_deadstemc_to_fire, m_deadstemc_storage_to_fire, m_deadstemc_xfer_to_fire,
        m_frootc_to_fire, m_frootc_storage_to_fire, m_frootc_xfer_to_fire,
        m_livecrootc_to_fire, m_livecrootc_storage_to_fire, m_livecrootc_xfer_to_fire,
        m_deadcrootc_to_fire, m_deadcrootc_storage_to_fire, m_deadcrootc_xfer_to_fire)

    nfo = _FireNFireOut(;
        m_leafn_to_fire, m_leafn_storage_to_fire, m_leafn_xfer_to_fire,
        m_livestemn_to_fire, m_livestemn_storage_to_fire, m_livestemn_xfer_to_fire,
        m_deadstemn_to_fire, m_deadstemn_storage_to_fire, m_deadstemn_xfer_to_fire,
        m_frootn_to_fire, m_frootn_storage_to_fire, m_frootn_xfer_to_fire,
        m_livecrootn_to_fire, m_livecrootn_storage_to_fire, m_livecrootn_xfer_to_fire,
        m_deadcrootn_to_fire, m_deadcrootn_storage_to_fire, m_deadcrootn_xfer_to_fire,
        m_retransn_to_fire)

    clit = _FireCLitOut(;
        m_leafc_to_litter_fire, m_leafc_storage_to_litter_fire, m_leafc_xfer_to_litter_fire,
        m_livestemc_to_litter_fire, m_livestemc_storage_to_litter_fire, m_livestemc_xfer_to_litter_fire,
        m_livestemc_to_deadstemc_fire, m_deadstemc_to_litter_fire,
        m_deadstemc_storage_to_litter_fire, m_deadstemc_xfer_to_litter_fire,
        m_frootc_to_litter_fire, m_frootc_storage_to_litter_fire, m_frootc_xfer_to_litter_fire,
        m_livecrootc_to_litter_fire, m_livecrootc_storage_to_litter_fire, m_livecrootc_xfer_to_litter_fire,
        m_livecrootc_to_deadcrootc_fire, m_deadcrootc_to_litter_fire,
        m_deadcrootc_storage_to_litter_fire, m_deadcrootc_xfer_to_litter_fire,
        m_gresp_storage_to_litter_fire, m_gresp_xfer_to_litter_fire)

    nlit = _FireNLitOut(;
        m_leafn_to_litter_fire, m_leafn_storage_to_litter_fire, m_leafn_xfer_to_litter_fire,
        m_livestemn_to_litter_fire, m_livestemn_storage_to_litter_fire, m_livestemn_xfer_to_litter_fire,
        m_livestemn_to_deadstemn_fire, m_deadstemn_to_litter_fire,
        m_deadstemn_storage_to_litter_fire, m_deadstemn_xfer_to_litter_fire,
        m_frootn_to_litter_fire, m_frootn_storage_to_litter_fire, m_frootn_xfer_to_litter_fire,
        m_livecrootn_to_litter_fire, m_livecrootn_storage_to_litter_fire, m_livecrootn_xfer_to_litter_fire,
        m_livecrootn_to_deadcrootn_fire, m_deadcrootn_to_litter_fire,
        m_deadcrootn_storage_to_litter_fire, m_deadcrootn_xfer_to_litter_fire,
        m_retransn_to_litter_fire)

    cndv = _FireCNDVOut(; nind, leafcmax, mask_actfirep)

    stc = _FireStateC(;
        cropf_col, farea_burned, fbac, baf_crop,
        gresp_storage, gresp_xfer,
        leafc, leafc_storage, leafc_xfer,
        livestemc, livestemc_storage, livestemc_xfer,
        deadstemc, deadstemc_storage, deadstemc_xfer,
        frootc, frootc_storage, frootc_xfer,
        livecrootc, livecrootc_storage, livecrootc_xfer,
        deadcrootc, deadcrootc_storage, deadcrootc_xfer)

    stn = _FireStateN(;
        leafn, leafn_storage, leafn_xfer,
        livestemn, livestemn_storage, livestemn_xfer,
        deadstemn, deadstemn_storage, deadstemn_xfer,
        frootn, frootn_storage, frootn_xfer,
        livecrootn, livecrootn_storage, livecrootn_xfer,
        deadcrootn, deadcrootn_storage, deadcrootn_xfer,
        retransn)

    pft = _FirePft(;
        woody, cc_leaf, cc_lstem, cc_dstem, cc_other,
        fm_leaf, fm_lstem, fm_other, fm_root, fm_lroot, fm_droot)

    # Scalar bundle at the working precision of the outputs (Float64 on CPU ->
    # byte-identical; Float32 on Metal). m_spinup and dt are Float64 host scalars.
    T = eltype(m_gresp_storage_to_fire)
    prm = _FirePatchScalars{T}(; m_spinup = T(m_spinup), dt = T(dt))

    # Struct-first launch: backend from a concrete device output float array
    # (mask_actfirep may be a BitVector with no backend, so do NOT key off it).
    backend = _kernel_backend(m_gresp_storage_to_fire)
    _fireb_cnfire_patchflux_kernel!(backend)(
        cfo, nfo, clit, nlit, cndv,
        mask_soilp, itype, column,
        stc, stn, pft,
        nc3crop, noveg, prm,
        transient_landcover, use_matrixcn, use_cndv;
        ndrange = n)
    KA.synchronize(backend)
    return mask_actfirep
end

# Section (2): per-patch litter/CWD PATCH->COLUMN SCATTER. Each masked patch
# accumulates (over decomp layers j and litter pools i) into per-column arrays
# via _scatter_add!. The column accumulators are NOT zeroed here — the host loop
# also accumulates into pre-zeroed arrays (no zero loop). The host's outer
# ordering is (j outer, p inner); since each (c,j[,i]) cell is accumulated once
# per patch, reordering to (p outer, j inner) keeps the same per-cell sum; on the
# KA CPU backend _scatter_add! iterates patches in p-order so the running sum is
# byte-identical to the host.
## Device-view structs for the litter/CWD patch->column scatter kernel.
## Grouped to defeat Metal's ~31-argument limit. All per-patch/per-column
## float arrays share eltype -> single V (vector) / M (matrix) / A3 (3D) params.
## NOTE: no scalar-bundle struct is needed here — every scalar arg of this
## kernel (nlevdecomp, i_met_lit, i_litr_max) is an Int, so they stay loose.
## The kernel body contains NO Float64 numeric literals (only integer +1 loop
## offsets), so no zero(T)/one(T)/T(x) conversions are required for byte-identity.

# Output accumulators (column-indexed; pre-zeroed by caller, accumulated via _scatter_add!)
Base.@kwdef struct _FireLSOut{M,A3}
    fire_mortality_c_to_cwdc::M   # [c, j]
    fire_mortality_n_to_cwdn::M   # [c, j]
    m_c_to_litr_fire::A3          # [c, j, i]
    m_n_to_litr_fire::A3          # [c, j, i]
end
Adapt.@adapt_structure _FireLSOut

# Read-only profiles + partition fractions + patch weights
Base.@kwdef struct _FireLSProf{V,M}
    stem_prof::M       # [p, j]
    croot_prof::M      # [p, j]
    leaf_prof::M       # [p, j]
    froot_prof::M      # [p, j]
    lf_f::M            # [ivt, i]
    fr_f::M            # [ivt, i]
    wtcol::V           # [p]
end
Adapt.@adapt_structure _FireLSProf

# Read-only per-patch CARBON mortality->litter fluxes
Base.@kwdef struct _FireLSCFlux{V}
    m_deadstemc_to_litter_fire::V
    m_deadcrootc_to_litter_fire::V
    m_livestemc_to_litter_fire::V
    m_livecrootc_to_litter_fire::V
    m_leafc_to_litter_fire::V
    m_leafc_storage_to_litter_fire::V
    m_leafc_xfer_to_litter_fire::V
    m_gresp_storage_to_litter_fire::V
    m_gresp_xfer_to_litter_fire::V
    m_frootc_to_litter_fire::V
    m_frootc_storage_to_litter_fire::V
    m_frootc_xfer_to_litter_fire::V
    m_livestemc_storage_to_litter_fire::V
    m_livestemc_xfer_to_litter_fire::V
    m_deadstemc_storage_to_litter_fire::V
    m_deadstemc_xfer_to_litter_fire::V
    m_livecrootc_storage_to_litter_fire::V
    m_livecrootc_xfer_to_litter_fire::V
    m_deadcrootc_storage_to_litter_fire::V
    m_deadcrootc_xfer_to_litter_fire::V
end
Adapt.@adapt_structure _FireLSCFlux

# Read-only per-patch NITROGEN mortality->litter fluxes
Base.@kwdef struct _FireLSNFlux{V}
    m_deadstemn_to_litter_fire::V
    m_deadcrootn_to_litter_fire::V
    m_livestemn_to_litter_fire::V
    m_livecrootn_to_litter_fire::V
    m_leafn_to_litter_fire::V
    m_leafn_storage_to_litter_fire::V
    m_leafn_xfer_to_litter_fire::V
    m_retransn_to_litter_fire::V
    m_frootn_to_litter_fire::V
    m_frootn_storage_to_litter_fire::V
    m_frootn_xfer_to_litter_fire::V
    m_livestemn_storage_to_litter_fire::V
    m_livestemn_xfer_to_litter_fire::V
    m_deadstemn_storage_to_litter_fire::V
    m_deadstemn_xfer_to_litter_fire::V
    m_livecrootn_storage_to_litter_fire::V
    m_livecrootn_xfer_to_litter_fire::V
    m_deadcrootn_storage_to_litter_fire::V
    m_deadcrootn_xfer_to_litter_fire::V
end
Adapt.@adapt_structure _FireLSNFlux

@kernel function _fireb_cnfire_litterscatter_kernel!(
    out::_FireLSOut, prof::_FireLSProf, cf::_FireLSCFlux, nf::_FireLSNFlux,
    @Const(mask_soilp), @Const(itype), @Const(column),
    nlevdecomp::Int, i_met_lit::Int, i_litr_max::Int
)
    p = @index(Global)
    @inbounds if mask_soilp[p]
        c = column[p]
        ivt = itype[p] + 1

        # --- alias struct fields to Fortran-named locals (body verbatim) ---
        fire_mortality_c_to_cwdc = out.fire_mortality_c_to_cwdc
        fire_mortality_n_to_cwdn = out.fire_mortality_n_to_cwdn
        m_c_to_litr_fire = out.m_c_to_litr_fire
        m_n_to_litr_fire = out.m_n_to_litr_fire

        stem_prof  = prof.stem_prof
        croot_prof = prof.croot_prof
        leaf_prof  = prof.leaf_prof
        froot_prof = prof.froot_prof
        lf_f = prof.lf_f
        fr_f = prof.fr_f
        w = prof.wtcol[p]

        m_deadstemc_to_litter_fire          = cf.m_deadstemc_to_litter_fire[p]
        m_deadcrootc_to_litter_fire         = cf.m_deadcrootc_to_litter_fire[p]
        m_livestemc_to_litter_fire          = cf.m_livestemc_to_litter_fire[p]
        m_livecrootc_to_litter_fire         = cf.m_livecrootc_to_litter_fire[p]
        m_leafc_to_litter_fire              = cf.m_leafc_to_litter_fire[p]
        m_leafc_storage_to_litter_fire      = cf.m_leafc_storage_to_litter_fire[p]
        m_leafc_xfer_to_litter_fire         = cf.m_leafc_xfer_to_litter_fire[p]
        m_gresp_storage_to_litter_fire      = cf.m_gresp_storage_to_litter_fire[p]
        m_gresp_xfer_to_litter_fire         = cf.m_gresp_xfer_to_litter_fire[p]
        m_frootc_to_litter_fire             = cf.m_frootc_to_litter_fire[p]
        m_frootc_storage_to_litter_fire     = cf.m_frootc_storage_to_litter_fire[p]
        m_frootc_xfer_to_litter_fire        = cf.m_frootc_xfer_to_litter_fire[p]
        m_livestemc_storage_to_litter_fire  = cf.m_livestemc_storage_to_litter_fire[p]
        m_livestemc_xfer_to_litter_fire     = cf.m_livestemc_xfer_to_litter_fire[p]
        m_deadstemc_storage_to_litter_fire  = cf.m_deadstemc_storage_to_litter_fire[p]
        m_deadstemc_xfer_to_litter_fire     = cf.m_deadstemc_xfer_to_litter_fire[p]
        m_livecrootc_storage_to_litter_fire = cf.m_livecrootc_storage_to_litter_fire[p]
        m_livecrootc_xfer_to_litter_fire    = cf.m_livecrootc_xfer_to_litter_fire[p]
        m_deadcrootc_storage_to_litter_fire = cf.m_deadcrootc_storage_to_litter_fire[p]
        m_deadcrootc_xfer_to_litter_fire    = cf.m_deadcrootc_xfer_to_litter_fire[p]

        m_deadstemn_to_litter_fire          = nf.m_deadstemn_to_litter_fire[p]
        m_deadcrootn_to_litter_fire         = nf.m_deadcrootn_to_litter_fire[p]
        m_livestemn_to_litter_fire          = nf.m_livestemn_to_litter_fire[p]
        m_livecrootn_to_litter_fire         = nf.m_livecrootn_to_litter_fire[p]
        m_leafn_to_litter_fire              = nf.m_leafn_to_litter_fire[p]
        m_leafn_storage_to_litter_fire      = nf.m_leafn_storage_to_litter_fire[p]
        m_leafn_xfer_to_litter_fire         = nf.m_leafn_xfer_to_litter_fire[p]
        m_retransn_to_litter_fire           = nf.m_retransn_to_litter_fire[p]
        m_frootn_to_litter_fire             = nf.m_frootn_to_litter_fire[p]
        m_frootn_storage_to_litter_fire     = nf.m_frootn_storage_to_litter_fire[p]
        m_frootn_xfer_to_litter_fire        = nf.m_frootn_xfer_to_litter_fire[p]
        m_livestemn_storage_to_litter_fire  = nf.m_livestemn_storage_to_litter_fire[p]
        m_livestemn_xfer_to_litter_fire     = nf.m_livestemn_xfer_to_litter_fire[p]
        m_deadstemn_storage_to_litter_fire  = nf.m_deadstemn_storage_to_litter_fire[p]
        m_deadstemn_xfer_to_litter_fire     = nf.m_deadstemn_xfer_to_litter_fire[p]
        m_livecrootn_storage_to_litter_fire = nf.m_livecrootn_storage_to_litter_fire[p]
        m_livecrootn_xfer_to_litter_fire    = nf.m_livecrootn_xfer_to_litter_fire[p]
        m_deadcrootn_storage_to_litter_fire = nf.m_deadcrootn_storage_to_litter_fire[p]
        m_deadcrootn_xfer_to_litter_fire    = nf.m_deadcrootn_xfer_to_litter_fire[p]

        for j in 1:nlevdecomp
            _scatter_add!(fire_mortality_c_to_cwdc, c, j,
                m_deadstemc_to_litter_fire * w * stem_prof[p, j])
            _scatter_add!(fire_mortality_c_to_cwdc, c, j,
                m_deadcrootc_to_litter_fire * w * croot_prof[p, j])
            _scatter_add!(fire_mortality_n_to_cwdn, c, j,
                m_deadstemn_to_litter_fire * w * stem_prof[p, j])
            _scatter_add!(fire_mortality_n_to_cwdn, c, j,
                m_deadcrootn_to_litter_fire * w * croot_prof[p, j])

            _scatter_add!(fire_mortality_c_to_cwdc, c, j,
                m_livestemc_to_litter_fire * w * stem_prof[p, j])
            _scatter_add!(fire_mortality_c_to_cwdc, c, j,
                m_livecrootc_to_litter_fire * w * croot_prof[p, j])
            _scatter_add!(fire_mortality_n_to_cwdn, c, j,
                m_livestemn_to_litter_fire * w * stem_prof[p, j])
            _scatter_add!(fire_mortality_n_to_cwdn, c, j,
                m_livecrootn_to_litter_fire * w * croot_prof[p, j])

            _scatter_add!(m_c_to_litr_fire, c, j, i_met_lit,
                ((m_leafc_to_litter_fire * lf_f[ivt, i_met_lit] +
                  m_leafc_storage_to_litter_fire +
                  m_leafc_xfer_to_litter_fire +
                  m_gresp_storage_to_litter_fire +
                  m_gresp_xfer_to_litter_fire) * leaf_prof[p, j] +
                 (m_frootc_to_litter_fire * fr_f[ivt, i_met_lit] +
                  m_frootc_storage_to_litter_fire +
                  m_frootc_xfer_to_litter_fire) * froot_prof[p, j] +
                 (m_livestemc_storage_to_litter_fire +
                  m_livestemc_xfer_to_litter_fire +
                  m_deadstemc_storage_to_litter_fire +
                  m_deadstemc_xfer_to_litter_fire) * stem_prof[p, j] +
                 (m_livecrootc_storage_to_litter_fire +
                  m_livecrootc_xfer_to_litter_fire +
                  m_deadcrootc_storage_to_litter_fire +
                  m_deadcrootc_xfer_to_litter_fire) * croot_prof[p, j]) * w)

            for i in (i_met_lit+1):i_litr_max
                _scatter_add!(m_c_to_litr_fire, c, j, i,
                    (m_leafc_to_litter_fire * lf_f[ivt, i] * leaf_prof[p, j] +
                     m_frootc_to_litter_fire * fr_f[ivt, i] * froot_prof[p, j]) * w)
            end

            _scatter_add!(m_n_to_litr_fire, c, j, i_met_lit,
                ((m_leafn_to_litter_fire * lf_f[ivt, i_met_lit] +
                  m_leafn_storage_to_litter_fire +
                  m_leafn_xfer_to_litter_fire +
                  m_retransn_to_litter_fire) * leaf_prof[p, j] +
                 (m_frootn_to_litter_fire * fr_f[ivt, i_met_lit] +
                  m_frootn_storage_to_litter_fire +
                  m_frootn_xfer_to_litter_fire) * froot_prof[p, j] +
                 (m_livestemn_storage_to_litter_fire +
                  m_livestemn_xfer_to_litter_fire +
                  m_deadstemn_storage_to_litter_fire +
                  m_deadstemn_xfer_to_litter_fire) * stem_prof[p, j] +
                 (m_livecrootn_storage_to_litter_fire +
                  m_livecrootn_xfer_to_litter_fire +
                  m_deadcrootn_storage_to_litter_fire +
                  m_deadcrootn_xfer_to_litter_fire) * croot_prof[p, j]) * w)

            for i in (i_met_lit+1):i_litr_max
                _scatter_add!(m_n_to_litr_fire, c, j, i,
                    (m_leafn_to_litter_fire * lf_f[ivt, i] * leaf_prof[p, j] +
                     m_frootn_to_litter_fire * fr_f[ivt, i] * froot_prof[p, j]) * w)
            end
        end
    end
end

function fireb_cnfire_litterscatter!(
        fire_mortality_c_to_cwdc, fire_mortality_n_to_cwdn,
        m_c_to_litr_fire, m_n_to_litr_fire,
        mask_soilp, itype, column, wtcol,
        stem_prof, croot_prof, leaf_prof, froot_prof,
        lf_f, fr_f,
        m_deadstemc_to_litter_fire, m_deadcrootc_to_litter_fire,
        m_deadstemn_to_litter_fire, m_deadcrootn_to_litter_fire,
        m_livestemc_to_litter_fire, m_livecrootc_to_litter_fire,
        m_livestemn_to_litter_fire, m_livecrootn_to_litter_fire,
        m_leafc_to_litter_fire, m_leafc_storage_to_litter_fire,
        m_leafc_xfer_to_litter_fire, m_gresp_storage_to_litter_fire,
        m_gresp_xfer_to_litter_fire,
        m_frootc_to_litter_fire, m_frootc_storage_to_litter_fire,
        m_frootc_xfer_to_litter_fire,
        m_livestemc_storage_to_litter_fire, m_livestemc_xfer_to_litter_fire,
        m_deadstemc_storage_to_litter_fire, m_deadstemc_xfer_to_litter_fire,
        m_livecrootc_storage_to_litter_fire, m_livecrootc_xfer_to_litter_fire,
        m_deadcrootc_storage_to_litter_fire, m_deadcrootc_xfer_to_litter_fire,
        m_leafn_to_litter_fire, m_leafn_storage_to_litter_fire,
        m_leafn_xfer_to_litter_fire, m_retransn_to_litter_fire,
        m_frootn_to_litter_fire, m_frootn_storage_to_litter_fire,
        m_frootn_xfer_to_litter_fire,
        m_livestemn_storage_to_litter_fire, m_livestemn_xfer_to_litter_fire,
        m_deadstemn_storage_to_litter_fire, m_deadstemn_xfer_to_litter_fire,
        m_livecrootn_storage_to_litter_fire, m_livecrootn_xfer_to_litter_fire,
        m_deadcrootn_storage_to_litter_fire, m_deadcrootn_xfer_to_litter_fire,
        nlevdecomp, i_met_lit, i_litr_max)

    isempty(mask_soilp) && return nothing

    out = _FireLSOut(;
        fire_mortality_c_to_cwdc, fire_mortality_n_to_cwdn,
        m_c_to_litr_fire, m_n_to_litr_fire)

    prof = _FireLSProf(;
        stem_prof, croot_prof, leaf_prof, froot_prof, lf_f, fr_f, wtcol)

    cf = _FireLSCFlux(;
        m_deadstemc_to_litter_fire, m_deadcrootc_to_litter_fire,
        m_livestemc_to_litter_fire, m_livecrootc_to_litter_fire,
        m_leafc_to_litter_fire, m_leafc_storage_to_litter_fire,
        m_leafc_xfer_to_litter_fire, m_gresp_storage_to_litter_fire,
        m_gresp_xfer_to_litter_fire,
        m_frootc_to_litter_fire, m_frootc_storage_to_litter_fire,
        m_frootc_xfer_to_litter_fire,
        m_livestemc_storage_to_litter_fire, m_livestemc_xfer_to_litter_fire,
        m_deadstemc_storage_to_litter_fire, m_deadstemc_xfer_to_litter_fire,
        m_livecrootc_storage_to_litter_fire, m_livecrootc_xfer_to_litter_fire,
        m_deadcrootc_storage_to_litter_fire, m_deadcrootc_xfer_to_litter_fire)

    nf = _FireLSNFlux(;
        m_deadstemn_to_litter_fire, m_deadcrootn_to_litter_fire,
        m_livestemn_to_litter_fire, m_livecrootn_to_litter_fire,
        m_leafn_to_litter_fire, m_leafn_storage_to_litter_fire,
        m_leafn_xfer_to_litter_fire, m_retransn_to_litter_fire,
        m_frootn_to_litter_fire, m_frootn_storage_to_litter_fire,
        m_frootn_xfer_to_litter_fire,
        m_livestemn_storage_to_litter_fire, m_livestemn_xfer_to_litter_fire,
        m_deadstemn_storage_to_litter_fire, m_deadstemn_xfer_to_litter_fire,
        m_livecrootn_storage_to_litter_fire, m_livecrootn_xfer_to_litter_fire,
        m_deadcrootn_storage_to_litter_fire, m_deadcrootn_xfer_to_litter_fire)

    # Struct-first launch: struct args carry no backend, so resolve it from an
    # output array and synchronize manually (mirrors the methane/phenology pattern).
    backend = _kernel_backend(out.fire_mortality_c_to_cwdc)
    _fireb_cnfire_litterscatter_kernel!(backend)(
        out, prof, cf, nf,
        mask_soilp, itype, column,
        nlevdecomp, i_met_lit, i_litr_max;
        ndrange = length(mask_soilp))
    KA.synchronize(backend)
    nothing
end

# Section (3): per-column vertically-resolved decomposing C/N fire loss.
# Independent per-column; internal j,l loops branch on is_litter[l]/is_cwd[l].
# Sets mask_actfirec[c] per column.
@kernel function _fireb_cnfire_decompfire_kernel!(
    m_decomp_cpools_to_fire_vr, m_decomp_npools_to_fire_vr, mask_actfirec,
    @Const(mask_soilc), @Const(farea_burned), @Const(baf_crop),
    @Const(decomp_cpools_vr), @Const(decomp_npools_vr),
    @Const(is_litter), @Const(is_cwd),
    nlevdecomp::Int, ndecomp_pools::Int,
    cmb_cmplt_fact_litter, cmb_cmplt_fact_cwd
)
    c = @index(Global)
    T = eltype(m_decomp_cpools_to_fire_vr)
    @inbounds if mask_soilc[c]
        f = farea_burned[c]
        if f != zero(T) || f != baf_crop[c]
            mask_actfirec[c] = true
        end
        for j in 1:nlevdecomp
            for l in 1:ndecomp_pools
                if is_litter[l]
                    m_decomp_cpools_to_fire_vr[c, j, l] = decomp_cpools_vr[c, j, l] * f *
                        cmb_cmplt_fact_litter
                end
                if is_cwd[l]
                    m_decomp_cpools_to_fire_vr[c, j, l] = decomp_cpools_vr[c, j, l] *
                        (f - baf_crop[c]) * cmb_cmplt_fact_cwd
                end
            end
            for l in 1:ndecomp_pools
                if is_litter[l]
                    m_decomp_npools_to_fire_vr[c, j, l] = decomp_npools_vr[c, j, l] * f *
                        cmb_cmplt_fact_litter
                end
                if is_cwd[l]
                    m_decomp_npools_to_fire_vr[c, j, l] = decomp_npools_vr[c, j, l] *
                        (f - baf_crop[c]) * cmb_cmplt_fact_cwd
                end
            end
        end
    end
end

fireb_cnfire_decompfire!(m_decomp_cpools_to_fire_vr, m_decomp_npools_to_fire_vr,
                         mask_actfirec, mask_soilc, args...) =
    _launch!(_fireb_cnfire_decompfire_kernel!, m_decomp_cpools_to_fire_vr,
             m_decomp_npools_to_fire_vr, mask_actfirec, mask_soilc, args...;
             ndrange = length(mask_soilc))

# Section (4): per-column deforestation-fire carbon loss. Gated by
# transient_landcover (host-resolved Bool) and the date predicate
# (kmo==1 && kda==1 && mcsec==0) which is uniform across columns and so is
# resolved once on the host to `is_newyear_start::Bool` and passed in.
@kernel function _fireb_cnfire_deforest_kernel!(
    lfc, lfc2,
    @Const(mask_soilc), @Const(trotr1_col), @Const(trotr2_col), @Const(dtrotr_col),
    @Const(fbac1), @Const(farea_burned), @Const(baf_crop), @Const(baf_peatf),
    is_newyear_start::Bool, dt, dayspyr, secspday
)
    c = @index(Global)
    T = eltype(lfc)
    @inbounds if mask_soilc[c]
        lfc2[c] = zero(T)
        if !is_newyear_start
            if trotr1_col[c] + trotr2_col[c] > T(0.6) && dtrotr_col[c] > zero(T) &&
               lfc[c] > zero(T) && fbac1[c] == zero(T)
                lfc2[c] = max(zero(T), min(lfc[c], (farea_burned[c] - baf_crop[c] -
                    baf_peatf[c]) / T(2) * dt)) / (dtrotr_col[c] * dayspyr * secspday / dt) / dt
                lfc[c]  = lfc[c] - max(zero(T), min(lfc[c], (farea_burned[c] - baf_crop[c] -
                    baf_peatf[c]) * dt / T(2)))
            end
        end
    end
end

fireb_cnfire_deforest!(lfc, args...) =
    _launch!(_fireb_cnfire_deforest_kernel!, lfc, args...)

"""
    cnfire_fluxes!(mask_soilc, mask_soilp, bounds_c, bounds_p,
        cnfire_const, pftcon, patch, col, grc,
        dgvs, cnveg_state, cnveg_cs, cnveg_cf, cnveg_ns, cnveg_nf,
        soilbgc_cf, decomp_cascade_con,
        leaf_prof_patch, froot_prof_patch, croot_prof_patch, stem_prof_patch,
        totsomc_col, decomp_cpools_vr_col, decomp_npools_vr_col, somc_fire_col;
        dt, dayspyr, nlevdecomp, ndecomp_pools, i_met_lit, i_litr_max,
        transient_landcover, use_cndv, use_matrixcn,
        spinup_factor_deadwood, secspday, nc3crop, noveg, kmo, kda, mcsec)

Fire effects routine for coupled carbon-nitrogen code (CN).
Relies primarily on estimate of fractional area burned.

Returns `(mask_actfirec, mask_actfirep)` — BitVectors indicating active fire
columns and patches.

Ported from `CNFireFluxes` in `CNFireBaseMod.F90`.
"""
function cnfire_fluxes!(
    mask_soilc::BitVector,
    mask_soilp::BitVector,
    bounds_c::UnitRange{Int},
    bounds_p::UnitRange{Int},
    cnfire_const::CNFireConstData,
    pftcon::PftConFireBase,
    patch::PatchData,
    col::ColumnData,
    grc::GridcellData,
    dgvs::DgvsFireData,
    cnveg_state::CNVegStateData,
    cnveg_cs::CNVegCarbonStateData,
    cnveg_cf::CNVegCarbonFluxData,
    cnveg_ns::CNVegNitrogenStateData,
    cnveg_nf::CNVegNitrogenFluxData,
    soilbgc_cf::SoilBiogeochemCarbonFluxData,
    decomp_cascade_con::DecompCascadeConData,
    leaf_prof_patch::Matrix{<:Real},
    froot_prof_patch::Matrix{<:Real},
    croot_prof_patch::Matrix{<:Real},
    stem_prof_patch::Matrix{<:Real},
    totsomc_col::Vector{<:Real},
    decomp_cpools_vr_col::Array{<:Real,3},
    decomp_npools_vr_col::Array{<:Real,3},
    somc_fire_col::Vector{<:Real};
    dt::Real = 1800.0,
    dayspyr::Real = 365.0,
    nlevdecomp::Int = 1,
    ndecomp_pools::Int = 7,
    i_met_lit::Int = 1,
    i_litr_max::Int = 3,
    transient_landcover::Bool = false,
    use_cndv::Bool = false,
    use_matrixcn::Bool = false,
    spinup_factor_deadwood::Real = SPINUP_FACTOR_DEADWOOD_DEFAULT,
    secspday::Real = SECSPDAY,
    nc3crop::Int = 15,
    noveg::Int = 0,
    kmo::Int = 6,
    kda::Int = 15,
    mcsec::Int = 3600
)
    # --- Aliases (matching Fortran associate block) ---
    # Profile arrays
    croot_prof   = croot_prof_patch
    stem_prof    = stem_prof_patch
    froot_prof   = froot_prof_patch
    leaf_prof    = leaf_prof_patch
    totsomc      = totsomc_col
    decomp_cpools_vr = decomp_cpools_vr_col
    decomp_npools_vr = decomp_npools_vr_col
    somc_fire    = somc_fire_col

    is_cwd     = decomp_cascade_con.is_cwd
    is_litter  = decomp_cascade_con.is_litter

    woody      = pftcon.woody
    cc_leaf    = pftcon.cc_leaf
    cc_lstem   = pftcon.cc_lstem
    cc_dstem   = pftcon.cc_dstem
    cc_other   = pftcon.cc_other
    fm_leaf    = pftcon.fm_leaf
    fm_lstem   = pftcon.fm_lstem
    fm_other   = pftcon.fm_other
    fm_root    = pftcon.fm_root
    fm_lroot   = pftcon.fm_lroot
    fm_droot   = pftcon.fm_droot
    lf_f       = pftcon.lf_f
    fr_f       = pftcon.fr_f

    cmb_cmplt_fact_litter = cnfire_const.cmb_cmplt_fact_litter
    cmb_cmplt_fact_cwd    = cnfire_const.cmb_cmplt_fact_cwd

    nind       = dgvs.nind_patch

    cropf_col      = cnveg_state.cropf_col
    farea_burned   = cnveg_state.farea_burned_col
    fbac1          = cnveg_state.fbac1_col
    fbac           = cnveg_state.fbac_col
    baf_crop       = cnveg_state.baf_crop_col
    baf_peatf      = cnveg_state.baf_peatf_col
    trotr1_col     = cnveg_state.trotr1_col
    trotr2_col     = cnveg_state.trotr2_col
    dtrotr_col     = cnveg_state.dtrotr_col
    lfc            = cnveg_state.lfc_col
    lfc2           = cnveg_state.lfc2_col

    leafcmax                    = cnveg_cs.leafcmax_patch
    leafc                       = cnveg_cs.leafc_patch
    leafc_storage               = cnveg_cs.leafc_storage_patch
    leafc_xfer                  = cnveg_cs.leafc_xfer_patch
    livestemc                   = cnveg_cs.livestemc_patch
    livestemc_storage           = cnveg_cs.livestemc_storage_patch
    livestemc_xfer              = cnveg_cs.livestemc_xfer_patch
    deadstemc                   = cnveg_cs.deadstemc_patch
    deadstemc_storage           = cnveg_cs.deadstemc_storage_patch
    deadstemc_xfer              = cnveg_cs.deadstemc_xfer_patch
    frootc                      = cnveg_cs.frootc_patch
    frootc_storage              = cnveg_cs.frootc_storage_patch
    frootc_xfer                 = cnveg_cs.frootc_xfer_patch
    livecrootc                  = cnveg_cs.livecrootc_patch
    livecrootc_storage          = cnveg_cs.livecrootc_storage_patch
    livecrootc_xfer             = cnveg_cs.livecrootc_xfer_patch
    deadcrootc                  = cnveg_cs.deadcrootc_patch
    deadcrootc_storage          = cnveg_cs.deadcrootc_storage_patch
    deadcrootc_xfer             = cnveg_cs.deadcrootc_xfer_patch
    gresp_storage               = cnveg_cs.gresp_storage_patch
    gresp_xfer                  = cnveg_cs.gresp_xfer_patch

    leafn                       = cnveg_ns.leafn_patch
    leafn_storage               = cnveg_ns.leafn_storage_patch
    leafn_xfer                  = cnveg_ns.leafn_xfer_patch
    livestemn                   = cnveg_ns.livestemn_patch
    livestemn_storage           = cnveg_ns.livestemn_storage_patch
    livestemn_xfer              = cnveg_ns.livestemn_xfer_patch
    deadstemn                   = cnveg_ns.deadstemn_patch
    deadstemn_storage           = cnveg_ns.deadstemn_storage_patch
    deadstemn_xfer              = cnveg_ns.deadstemn_xfer_patch
    frootn                      = cnveg_ns.frootn_patch
    frootn_storage              = cnveg_ns.frootn_storage_patch
    frootn_xfer                 = cnveg_ns.frootn_xfer_patch
    livecrootn                  = cnveg_ns.livecrootn_patch
    livecrootn_storage          = cnveg_ns.livecrootn_storage_patch
    livecrootn_xfer             = cnveg_ns.livecrootn_xfer_patch
    deadcrootn                  = cnveg_ns.deadcrootn_patch
    deadcrootn_storage          = cnveg_ns.deadcrootn_storage_patch
    deadcrootn_xfer             = cnveg_ns.deadcrootn_xfer_patch
    retransn                    = cnveg_ns.retransn_patch

    # Carbon fluxes: fire emissions
    m_leafc_to_fire               = cnveg_cf.m_leafc_to_fire_patch
    m_leafc_storage_to_fire       = cnveg_cf.m_leafc_storage_to_fire_patch
    m_leafc_xfer_to_fire          = cnveg_cf.m_leafc_xfer_to_fire_patch
    m_livestemc_to_fire           = cnveg_cf.m_livestemc_to_fire_patch
    m_livestemc_storage_to_fire   = cnveg_cf.m_livestemc_storage_to_fire_patch
    m_livestemc_xfer_to_fire      = cnveg_cf.m_livestemc_xfer_to_fire_patch
    m_deadstemc_to_fire           = cnveg_cf.m_deadstemc_to_fire_patch
    m_deadstemc_storage_to_fire   = cnveg_cf.m_deadstemc_storage_to_fire_patch
    m_deadstemc_xfer_to_fire      = cnveg_cf.m_deadstemc_xfer_to_fire_patch
    m_frootc_to_fire              = cnveg_cf.m_frootc_to_fire_patch
    m_frootc_storage_to_fire      = cnveg_cf.m_frootc_storage_to_fire_patch
    m_frootc_xfer_to_fire         = cnveg_cf.m_frootc_xfer_to_fire_patch
    m_livecrootc_to_fire          = cnveg_cf.m_livecrootc_to_fire_patch
    m_livecrootc_storage_to_fire  = cnveg_cf.m_livecrootc_storage_to_fire_patch
    m_livecrootc_xfer_to_fire     = cnveg_cf.m_livecrootc_xfer_to_fire_patch
    m_deadcrootc_to_fire          = cnveg_cf.m_deadcrootc_to_fire_patch
    m_deadcrootc_storage_to_fire  = cnveg_cf.m_deadcrootc_storage_to_fire_patch
    m_deadcrootc_xfer_to_fire     = cnveg_cf.m_deadcrootc_xfer_to_fire_patch
    m_gresp_storage_to_fire       = cnveg_cf.m_gresp_storage_to_fire_patch
    m_gresp_xfer_to_fire          = cnveg_cf.m_gresp_xfer_to_fire_patch

    # Carbon fluxes: fire mortality to litter
    m_leafc_to_litter_fire              = cnveg_cf.m_leafc_to_litter_fire_patch
    m_leafc_storage_to_litter_fire      = cnveg_cf.m_leafc_storage_to_litter_fire_patch
    m_leafc_xfer_to_litter_fire         = cnveg_cf.m_leafc_xfer_to_litter_fire_patch
    m_livestemc_to_litter_fire          = cnveg_cf.m_livestemc_to_litter_fire_patch
    m_livestemc_storage_to_litter_fire  = cnveg_cf.m_livestemc_storage_to_litter_fire_patch
    m_livestemc_xfer_to_litter_fire     = cnveg_cf.m_livestemc_xfer_to_litter_fire_patch
    m_livestemc_to_deadstemc_fire       = cnveg_cf.m_livestemc_to_deadstemc_fire_patch
    m_deadstemc_to_litter_fire          = cnveg_cf.m_deadstemc_to_litter_fire_patch
    m_deadstemc_storage_to_litter_fire  = cnveg_cf.m_deadstemc_storage_to_litter_fire_patch
    m_deadstemc_xfer_to_litter_fire     = cnveg_cf.m_deadstemc_xfer_to_litter_fire_patch
    m_frootc_to_litter_fire             = cnveg_cf.m_frootc_to_litter_fire_patch
    m_frootc_storage_to_litter_fire     = cnveg_cf.m_frootc_storage_to_litter_fire_patch
    m_frootc_xfer_to_litter_fire        = cnveg_cf.m_frootc_xfer_to_litter_fire_patch
    m_livecrootc_to_litter_fire         = cnveg_cf.m_livecrootc_to_litter_fire_patch
    m_livecrootc_storage_to_litter_fire = cnveg_cf.m_livecrootc_storage_to_litter_fire_patch
    m_livecrootc_xfer_to_litter_fire    = cnveg_cf.m_livecrootc_xfer_to_litter_fire_patch
    m_livecrootc_to_deadcrootc_fire     = cnveg_cf.m_livecrootc_to_deadcrootc_fire_patch
    m_deadcrootc_to_litter_fire         = cnveg_cf.m_deadcrootc_to_litter_fire_patch
    m_deadcrootc_storage_to_litter_fire = cnveg_cf.m_deadcrootc_storage_to_litter_fire_patch
    m_deadcrootc_xfer_to_litter_fire    = cnveg_cf.m_deadcrootc_xfer_to_litter_fire_patch
    m_gresp_storage_to_litter_fire      = cnveg_cf.m_gresp_storage_to_litter_fire_patch
    m_gresp_xfer_to_litter_fire         = cnveg_cf.m_gresp_xfer_to_litter_fire_patch

    m_decomp_cpools_to_fire_vr  = cnveg_cf.m_decomp_cpools_to_fire_vr_col
    m_c_to_litr_fire            = cnveg_cf.m_c_to_litr_fire_col
    fire_mortality_c_to_cwdc    = cnveg_cf.fire_mortality_c_to_cwdc_col

    # Nitrogen fluxes: fire emissions
    m_leafn_to_fire               = cnveg_nf.m_leafn_to_fire_patch
    m_leafn_storage_to_fire       = cnveg_nf.m_leafn_storage_to_fire_patch
    m_leafn_xfer_to_fire          = cnveg_nf.m_leafn_xfer_to_fire_patch
    m_livestemn_to_fire           = cnveg_nf.m_livestemn_to_fire_patch
    m_livestemn_storage_to_fire   = cnveg_nf.m_livestemn_storage_to_fire_patch
    m_livestemn_xfer_to_fire      = cnveg_nf.m_livestemn_xfer_to_fire_patch
    m_deadstemn_to_fire           = cnveg_nf.m_deadstemn_to_fire_patch
    m_deadstemn_storage_to_fire   = cnveg_nf.m_deadstemn_storage_to_fire_patch
    m_deadstemn_xfer_to_fire      = cnveg_nf.m_deadstemn_xfer_to_fire_patch
    m_frootn_to_fire              = cnveg_nf.m_frootn_to_fire_patch
    m_frootn_storage_to_fire      = cnveg_nf.m_frootn_storage_to_fire_patch
    m_frootn_xfer_to_fire         = cnveg_nf.m_frootn_xfer_to_fire_patch
    m_livecrootn_to_fire          = cnveg_nf.m_livecrootn_to_fire_patch
    m_livecrootn_storage_to_fire  = cnveg_nf.m_livecrootn_storage_to_fire_patch
    m_livecrootn_xfer_to_fire     = cnveg_nf.m_livecrootn_xfer_to_fire_patch
    m_deadcrootn_to_fire          = cnveg_nf.m_deadcrootn_to_fire_patch
    m_deadcrootn_storage_to_fire  = cnveg_nf.m_deadcrootn_storage_to_fire_patch
    m_deadcrootn_xfer_to_fire     = cnveg_nf.m_deadcrootn_xfer_to_fire_patch
    m_retransn_to_fire            = cnveg_nf.m_retransn_to_fire_patch

    # Nitrogen fluxes: fire mortality to litter
    m_leafn_to_litter_fire              = cnveg_nf.m_leafn_to_litter_fire_patch
    m_leafn_storage_to_litter_fire      = cnveg_nf.m_leafn_storage_to_litter_fire_patch
    m_leafn_xfer_to_litter_fire         = cnveg_nf.m_leafn_xfer_to_litter_fire_patch
    m_livestemn_to_litter_fire          = cnveg_nf.m_livestemn_to_litter_fire_patch
    m_livestemn_storage_to_litter_fire  = cnveg_nf.m_livestemn_storage_to_litter_fire_patch
    m_livestemn_xfer_to_litter_fire     = cnveg_nf.m_livestemn_xfer_to_litter_fire_patch
    m_livestemn_to_deadstemn_fire       = cnveg_nf.m_livestemn_to_deadstemn_fire_patch
    m_deadstemn_to_litter_fire          = cnveg_nf.m_deadstemn_to_litter_fire_patch
    m_deadstemn_storage_to_litter_fire  = cnveg_nf.m_deadstemn_storage_to_litter_fire_patch
    m_deadstemn_xfer_to_litter_fire     = cnveg_nf.m_deadstemn_xfer_to_litter_fire_patch
    m_frootn_to_litter_fire             = cnveg_nf.m_frootn_to_litter_fire_patch
    m_frootn_storage_to_litter_fire     = cnveg_nf.m_frootn_storage_to_litter_fire_patch
    m_frootn_xfer_to_litter_fire        = cnveg_nf.m_frootn_xfer_to_litter_fire_patch
    m_livecrootn_to_litter_fire         = cnveg_nf.m_livecrootn_to_litter_fire_patch
    m_livecrootn_storage_to_litter_fire = cnveg_nf.m_livecrootn_storage_to_litter_fire_patch
    m_livecrootn_xfer_to_litter_fire    = cnveg_nf.m_livecrootn_xfer_to_litter_fire_patch
    m_livecrootn_to_deadcrootn_fire     = cnveg_nf.m_livecrootn_to_deadcrootn_fire_patch
    m_deadcrootn_to_litter_fire         = cnveg_nf.m_deadcrootn_to_litter_fire_patch
    m_deadcrootn_storage_to_litter_fire = cnveg_nf.m_deadcrootn_storage_to_litter_fire_patch
    m_deadcrootn_xfer_to_litter_fire    = cnveg_nf.m_deadcrootn_xfer_to_litter_fire_patch
    m_retransn_to_litter_fire           = cnveg_nf.m_retransn_to_litter_fire_patch

    m_decomp_npools_to_fire_vr  = cnveg_nf.m_decomp_npools_to_fire_vr_col
    m_n_to_litr_fire            = cnveg_nf.m_n_to_litr_fire_col
    fire_mortality_n_to_cwdn    = cnveg_nf.fire_mortality_n_to_cwdn_col

    # Build output active-fire masks
    mask_actfirep = falses(length(mask_soilp))
    mask_actfirec = falses(length(mask_soilc))

    # ===================================================================
    # Patch loop — fire emissions and mortality fluxes (kernelized)
    # ===================================================================
    m = spinup_factor_deadwood

    fireb_cnfire_patchflux!(
        mask_actfirep,
        m_gresp_storage_to_fire, m_gresp_xfer_to_fire,
        m_leafc_to_fire, m_leafc_storage_to_fire, m_leafc_xfer_to_fire,
        m_livestemc_to_fire, m_livestemc_storage_to_fire, m_livestemc_xfer_to_fire,
        m_deadstemc_to_fire, m_deadstemc_storage_to_fire, m_deadstemc_xfer_to_fire,
        m_frootc_to_fire, m_frootc_storage_to_fire, m_frootc_xfer_to_fire,
        m_livecrootc_to_fire, m_livecrootc_storage_to_fire, m_livecrootc_xfer_to_fire,
        m_deadcrootc_to_fire, m_deadcrootc_storage_to_fire, m_deadcrootc_xfer_to_fire,
        m_leafn_to_fire, m_leafn_storage_to_fire, m_leafn_xfer_to_fire,
        m_livestemn_to_fire, m_livestemn_storage_to_fire, m_livestemn_xfer_to_fire,
        m_deadstemn_to_fire, m_deadstemn_storage_to_fire, m_deadstemn_xfer_to_fire,
        m_frootn_to_fire, m_frootn_storage_to_fire, m_frootn_xfer_to_fire,
        m_livecrootn_to_fire, m_livecrootn_storage_to_fire, m_livecrootn_xfer_to_fire,
        m_deadcrootn_to_fire, m_deadcrootn_storage_to_fire, m_deadcrootn_xfer_to_fire,
        m_retransn_to_fire,
        m_leafc_to_litter_fire, m_leafc_storage_to_litter_fire, m_leafc_xfer_to_litter_fire,
        m_livestemc_to_litter_fire, m_livestemc_storage_to_litter_fire, m_livestemc_xfer_to_litter_fire,
        m_livestemc_to_deadstemc_fire, m_deadstemc_to_litter_fire,
        m_deadstemc_storage_to_litter_fire, m_deadstemc_xfer_to_litter_fire,
        m_frootc_to_litter_fire, m_frootc_storage_to_litter_fire, m_frootc_xfer_to_litter_fire,
        m_livecrootc_to_litter_fire, m_livecrootc_storage_to_litter_fire, m_livecrootc_xfer_to_litter_fire,
        m_livecrootc_to_deadcrootc_fire, m_deadcrootc_to_litter_fire,
        m_deadcrootc_storage_to_litter_fire, m_deadcrootc_xfer_to_litter_fire,
        m_gresp_storage_to_litter_fire, m_gresp_xfer_to_litter_fire,
        m_leafn_to_litter_fire, m_leafn_storage_to_litter_fire, m_leafn_xfer_to_litter_fire,
        m_livestemn_to_litter_fire, m_livestemn_storage_to_litter_fire, m_livestemn_xfer_to_litter_fire,
        m_livestemn_to_deadstemn_fire, m_deadstemn_to_litter_fire,
        m_deadstemn_storage_to_litter_fire, m_deadstemn_xfer_to_litter_fire,
        m_frootn_to_litter_fire, m_frootn_storage_to_litter_fire, m_frootn_xfer_to_litter_fire,
        m_livecrootn_to_litter_fire, m_livecrootn_storage_to_litter_fire, m_livecrootn_xfer_to_litter_fire,
        m_livecrootn_to_deadcrootn_fire, m_deadcrootn_to_litter_fire,
        m_deadcrootn_storage_to_litter_fire, m_deadcrootn_xfer_to_litter_fire,
        m_retransn_to_litter_fire,
        nind, leafcmax,
        mask_soilp, patch.itype, patch.column,
        cropf_col, farea_burned, fbac, baf_crop,
        gresp_storage, gresp_xfer,
        leafc, leafc_storage, leafc_xfer,
        livestemc, livestemc_storage, livestemc_xfer,
        deadstemc, deadstemc_storage, deadstemc_xfer,
        frootc, frootc_storage, frootc_xfer,
        livecrootc, livecrootc_storage, livecrootc_xfer,
        deadcrootc, deadcrootc_storage, deadcrootc_xfer,
        leafn, leafn_storage, leafn_xfer,
        livestemn, livestemn_storage, livestemn_xfer,
        deadstemn, deadstemn_storage, deadstemn_xfer,
        frootn, frootn_storage, frootn_xfer,
        livecrootn, livecrootn_storage, livecrootn_xfer,
        deadcrootn, deadcrootn_storage, deadcrootn_xfer,
        retransn,
        woody, cc_leaf, cc_lstem, cc_dstem, cc_other,
        fm_leaf, fm_lstem, fm_other, fm_root, fm_lroot, fm_droot,
        nc3crop, noveg, m, dt,
        transient_landcover, use_matrixcn, use_cndv)

    # ===================================================================
    # Fire-induced transfer of C and N pools to litter and CWD (kernelized
    # patch->column scatter; accumulators arrive pre-zeroed, not re-zeroed here)
    # ===================================================================
    fireb_cnfire_litterscatter!(
        fire_mortality_c_to_cwdc, fire_mortality_n_to_cwdn,
        m_c_to_litr_fire, m_n_to_litr_fire,
        mask_soilp, patch.itype, patch.column, patch.wtcol,
        stem_prof, croot_prof, leaf_prof, froot_prof,
        lf_f, fr_f,
        m_deadstemc_to_litter_fire, m_deadcrootc_to_litter_fire,
        m_deadstemn_to_litter_fire, m_deadcrootn_to_litter_fire,
        m_livestemc_to_litter_fire, m_livecrootc_to_litter_fire,
        m_livestemn_to_litter_fire, m_livecrootn_to_litter_fire,
        m_leafc_to_litter_fire, m_leafc_storage_to_litter_fire,
        m_leafc_xfer_to_litter_fire, m_gresp_storage_to_litter_fire,
        m_gresp_xfer_to_litter_fire,
        m_frootc_to_litter_fire, m_frootc_storage_to_litter_fire,
        m_frootc_xfer_to_litter_fire,
        m_livestemc_storage_to_litter_fire, m_livestemc_xfer_to_litter_fire,
        m_deadstemc_storage_to_litter_fire, m_deadstemc_xfer_to_litter_fire,
        m_livecrootc_storage_to_litter_fire, m_livecrootc_xfer_to_litter_fire,
        m_deadcrootc_storage_to_litter_fire, m_deadcrootc_xfer_to_litter_fire,
        m_leafn_to_litter_fire, m_leafn_storage_to_litter_fire,
        m_leafn_xfer_to_litter_fire, m_retransn_to_litter_fire,
        m_frootn_to_litter_fire, m_frootn_storage_to_litter_fire,
        m_frootn_xfer_to_litter_fire,
        m_livestemn_storage_to_litter_fire, m_livestemn_xfer_to_litter_fire,
        m_deadstemn_storage_to_litter_fire, m_deadstemn_xfer_to_litter_fire,
        m_livecrootn_storage_to_litter_fire, m_livecrootn_xfer_to_litter_fire,
        m_deadcrootn_storage_to_litter_fire, m_deadcrootn_xfer_to_litter_fire,
        nlevdecomp, i_met_lit, i_litr_max)

    # ===================================================================
    # Vertically-resolved decomposing C/N fire loss — per-column kernel
    # ===================================================================
    fireb_cnfire_decompfire!(
        m_decomp_cpools_to_fire_vr, m_decomp_npools_to_fire_vr, mask_actfirec,
        mask_soilc, farea_burned, baf_crop,
        decomp_cpools_vr, decomp_npools_vr,
        is_litter, is_cwd,
        nlevdecomp, ndecomp_pools,
        cmb_cmplt_fact_litter, cmb_cmplt_fact_cwd)

    # ===================================================================
    # Carbon loss due to deforestation fires (per-column kernel; the
    # new-year-start date predicate is uniform across columns -> host-resolved)
    # ===================================================================
    if transient_landcover
        is_newyear_start = (kmo == 1 && kda == 1 && mcsec == 0)
        fireb_cnfire_deforest!(
            lfc, lfc2,
            mask_soilc, trotr1_col, trotr2_col, dtrotr_col,
            fbac1, farea_burned, baf_crop, baf_peatf,
            is_newyear_start, dt, dayspyr, secspday)
    end

    # ===================================================================
    # Carbon loss due to peat fires
    # ===================================================================
    fireb_peatfire_somc!(somc_fire, mask_soilc, col.gridcell, grc.latdeg,
                         totsomc, baf_peatf, first(bounds_c), last(bounds_c),
                         cnfire_const.borealat)

    return (mask_actfirec, mask_actfirep)
end
