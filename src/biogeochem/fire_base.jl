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
    @inbounds if cmin <= p <= cmax && mask[p]
        btran2[p] = 0.0
    end
end

fireb_zero_btran2!(btran2, mask, cmin::Int, cmax::Int) =
    _launch!(_fireb_zero_btran2_kernel!, btran2, mask, cmin, cmax)

# Clamp btran2[p] to <= 1.0 for every masked patch in [cmin, cmax].
@kernel function _fireb_clamp_btran2_kernel!(btran2, @Const(mask), cmin::Int, cmax::Int)
    p = @index(Global)
    @inbounds if cmin <= p <= cmax && mask[p]
        if btran2[p] > 1.0
            btran2[p] = 1.0
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
    @inbounds if cmin <= c <= cmax && mask[c]
        g = gridcell[c]
        if latdeg[g] < borealat
            somc_fire[c] = totsomc[c] * baf_peatf[c] * 6.0 / 33.9
        else
            somc_fire[c] = baf_peatf[c] * 2.2e3
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
Base.@kwdef mutable struct PftConFireBase
    woody    ::Vector{Float64} = Float64[]  # binary woody flag (1=woody, 0=not woody)
    cc_leaf  ::Vector{Float64} = Float64[]  # combustion completeness factor for leaves
    cc_lstem ::Vector{Float64} = Float64[]  # combustion completeness factor for live stems
    cc_dstem ::Vector{Float64} = Float64[]  # combustion completeness factor for dead stems
    cc_other ::Vector{Float64} = Float64[]  # combustion completeness factor for other pools
    fm_leaf  ::Vector{Float64} = Float64[]  # fire mortality factor for leaves
    fm_lstem ::Vector{Float64} = Float64[]  # fire mortality factor for live stems
    fm_other ::Vector{Float64} = Float64[]  # fire mortality factor for other pools
    fm_root  ::Vector{Float64} = Float64[]  # fire mortality factor for fine roots
    fm_lroot ::Vector{Float64} = Float64[]  # fire mortality factor for live coarse roots
    fm_droot ::Vector{Float64} = Float64[]  # fire mortality factor for dead coarse roots
    lf_f     ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # leaf litter fractions (pft, litr)
    fr_f     ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # fine root litter fractions (pft, litr)
    smpso    ::Vector{Float64} = Float64[]  # soil water potential at full stomatal opening (mm)
    smpsc    ::Vector{Float64} = Float64[]  # soil water potential at full stomatal closure (mm)
end

# ---------------------------------------------------------------------------
# Fire base data (btran2_patch)
# ---------------------------------------------------------------------------

"""
    CNFireBaseData

Fire base data containing root zone soil wetness factor.
Ported from `cnfire_base_type` member data in `CNFireBaseMod.F90`.
"""
Base.@kwdef mutable struct CNFireBaseData
    btran2_patch ::Vector{Float64} = Float64[]  # patch root zone soil wetness factor (0 to 1)
end

# ---------------------------------------------------------------------------
# DGVS data needed by fire
# ---------------------------------------------------------------------------

"""
    DgvsFireData

DGVS fields needed by fire routines (nind_patch, leafcmax_patch).
"""
Base.@kwdef mutable struct DgvsFireData
    nind_patch ::Vector{Float64} = Float64[]  # number of individuals (#/m2)
end

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
    @inbounds if cmin <= p <= cmax && mask[p]
        c = column[p]
        ivt = itype[p] + 1  # 0-based Fortran → 1-based Julia
        acc = btran2[p]
        for j in 1:nlevgrnd
            s_node = max(h2osoi_vol_col[c, j] / watsat[c, j], 0.01)
            # default_soil_suction inlined: smp = -sucsat * s^(-bsw)
            smp_node_lf = -sucsat[c, j] * s_node^(-bsw[c, j])
            smp_node_lf = max(smpsc[ivt], smp_node_lf)
            acc = acc + rootfr[p, j] * max(0.0,
                min((smp_node_lf - smpsc[ivt]) / (smpso[ivt] - smpsc[ivt]), 1.0))
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
    @inbounds if cmin <= p <= cmax && mask[p]
        c = column[p]
        acc = btran2[p]
        for j in 1:nlevgrnd
            s_node = max(h2osoi_vol_col[c, j] / watsat[c, j], 0.01)
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
@kernel function _fireb_cnfire_patchflux_kernel!(
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
    @Const(mask_soilp), @Const(itype), @Const(column),
    @Const(cropf_col), @Const(farea_burned), @Const(fbac), @Const(baf_crop),
    @Const(gresp_storage), @Const(gresp_xfer),
    @Const(leafc), @Const(leafc_storage), @Const(leafc_xfer),
    @Const(livestemc), @Const(livestemc_storage), @Const(livestemc_xfer),
    @Const(deadstemc), @Const(deadstemc_storage), @Const(deadstemc_xfer),
    @Const(frootc), @Const(frootc_storage), @Const(frootc_xfer),
    @Const(livecrootc), @Const(livecrootc_storage), @Const(livecrootc_xfer),
    @Const(deadcrootc), @Const(deadcrootc_storage), @Const(deadcrootc_xfer),
    @Const(leafn), @Const(leafn_storage), @Const(leafn_xfer),
    @Const(livestemn), @Const(livestemn_storage), @Const(livestemn_xfer),
    @Const(deadstemn), @Const(deadstemn_storage), @Const(deadstemn_xfer),
    @Const(frootn), @Const(frootn_storage), @Const(frootn_xfer),
    @Const(livecrootn), @Const(livecrootn_storage), @Const(livecrootn_xfer),
    @Const(deadcrootn), @Const(deadcrootn_storage), @Const(deadcrootn_xfer),
    @Const(retransn),
    @Const(woody), @Const(cc_leaf), @Const(cc_lstem), @Const(cc_dstem), @Const(cc_other),
    @Const(fm_leaf), @Const(fm_lstem), @Const(fm_other), @Const(fm_root),
    @Const(fm_lroot), @Const(fm_droot),
    nc3crop::Int, noveg::Int, m_spinup, dt,
    transient_landcover::Bool, use_matrixcn::Bool, use_cndv::Bool
)
    p = @index(Global)
    @inbounds if mask_soilp[p]
        c = column[p]
        m = m_spinup

        if itype[p] < nc3crop && cropf_col[c] < 1.0
            if transient_landcover
                f = (fbac[c] - baf_crop[c]) / (1.0 - cropf_col[c])
            else
                f = (farea_burned[c] - baf_crop[c]) / (1.0 - cropf_col[c])
            end
        else
            if cropf_col[c] > 0.0
                f = baf_crop[c] / cropf_col[c]
            else
                f = 0.0
            end
        end

        if f != 0
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
            m_frootc_to_fire[p]              = frootc[p]             * f * 0.0
            m_frootc_storage_to_fire[p]      = frootc_storage[p]     * f * cc_other[ivt]
            m_frootc_xfer_to_fire[p]         = frootc_xfer[p]        * f * cc_other[ivt]
            m_livecrootc_to_fire[p]          = livecrootc[p]         * f * 0.0
            m_livecrootc_storage_to_fire[p]  = livecrootc_storage[p] * f * cc_other[ivt]
            m_livecrootc_xfer_to_fire[p]     = livecrootc_xfer[p]    * f * cc_other[ivt]
            m_deadcrootc_to_fire[p]          = deadcrootc[p]         * f * 0.0
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
            m_frootn_to_fire[p]              = frootn[p]             * f * 0.0
            m_frootn_storage_to_fire[p]      = frootn_storage[p]     * f * cc_other[ivt]
            m_frootn_xfer_to_fire[p]         = frootn_xfer[p]        * f * cc_other[ivt]
            m_livecrootn_to_fire[p]          = livecrootn[p]         * f * 0.0
            m_livecrootn_storage_to_fire[p]  = livecrootn_storage[p] * f * cc_other[ivt]
            m_livecrootn_xfer_to_fire[p]     = livecrootn_xfer[p]    * f * cc_other[ivt]
            m_deadcrootn_to_fire[p]          = deadcrootn[p]         * f * 0.0
            m_deadcrootn_xfer_to_fire[p]     = deadcrootn_xfer[p]    * f * cc_other[ivt]
            m_deadcrootn_storage_to_fire[p]  = deadcrootn_storage[p] * f * cc_other[ivt]
            m_retransn_to_fire[p]            = retransn[p]           * f * cc_other[ivt]
        end

        if !use_matrixcn
            m_leafc_to_litter_fire[p]                   = leafc[p] * f *
                (1.0 - cc_leaf[ivt]) * fm_leaf[ivt]
            m_leafc_storage_to_litter_fire[p]           = leafc_storage[p] * f *
                (1.0 - cc_other[ivt]) * fm_other[ivt]
            m_leafc_xfer_to_litter_fire[p]              = leafc_xfer[p] * f *
                (1.0 - cc_other[ivt]) * fm_other[ivt]
            m_livestemc_to_litter_fire[p]               = livestemc[p] * f *
                (1.0 - cc_lstem[ivt]) * fm_droot[ivt]
            m_livestemc_storage_to_litter_fire[p]       = livestemc_storage[p] * f *
                (1.0 - cc_other[ivt]) * fm_other[ivt]
            m_livestemc_xfer_to_litter_fire[p]          = livestemc_xfer[p] * f *
                (1.0 - cc_other[ivt]) * fm_other[ivt]
            m_livestemc_to_deadstemc_fire[p]            = livestemc[p] * f *
                (1.0 - cc_lstem[ivt]) * (fm_lstem[ivt] - fm_droot[ivt])
            m_deadstemc_to_litter_fire[p]               = deadstemc[p] * f * m *
                (1.0 - cc_dstem[ivt]) * fm_droot[ivt]
            m_deadstemc_storage_to_litter_fire[p]       = deadstemc_storage[p] * f *
                (1.0 - cc_other[ivt]) * fm_other[ivt]
            m_deadstemc_xfer_to_litter_fire[p]          = deadstemc_xfer[p] * f *
                (1.0 - cc_other[ivt]) * fm_other[ivt]
            m_frootc_to_litter_fire[p]                  = frootc[p] * f *
                fm_root[ivt]
            m_frootc_storage_to_litter_fire[p]          = frootc_storage[p] * f *
                (1.0 - cc_other[ivt]) * fm_other[ivt]
            m_frootc_xfer_to_litter_fire[p]             = frootc_xfer[p] * f *
                (1.0 - cc_other[ivt]) * fm_other[ivt]
            m_livecrootc_to_litter_fire[p]              = livecrootc[p] * f *
                fm_droot[ivt]
            m_livecrootc_storage_to_litter_fire[p]      = livecrootc_storage[p] * f *
                (1.0 - cc_other[ivt]) * fm_other[ivt]
            m_livecrootc_xfer_to_litter_fire[p]         = livecrootc_xfer[p] * f *
                (1.0 - cc_other[ivt]) * fm_other[ivt]
            m_livecrootc_to_deadcrootc_fire[p]          = livecrootc[p] * f *
                (fm_lroot[ivt] - fm_droot[ivt])
            m_deadcrootc_to_litter_fire[p]              = deadcrootc[p] * f * m *
                fm_droot[ivt]
            m_deadcrootc_storage_to_litter_fire[p]      = deadcrootc_storage[p] * f *
                (1.0 - cc_other[ivt]) * fm_other[ivt]
            m_deadcrootc_xfer_to_litter_fire[p]         = deadcrootc_xfer[p] * f *
                (1.0 - cc_other[ivt]) * fm_other[ivt]
            m_gresp_storage_to_litter_fire[p]           = gresp_storage[p] * f *
                (1.0 - cc_other[ivt]) * fm_other[ivt]
            m_gresp_xfer_to_litter_fire[p]              = gresp_xfer[p] * f *
                (1.0 - cc_other[ivt]) * fm_other[ivt]

            m_leafn_to_litter_fire[p]                  = leafn[p] * f *
                (1.0 - cc_leaf[ivt]) * fm_leaf[ivt]
            m_leafn_storage_to_litter_fire[p]          = leafn_storage[p] * f *
                (1.0 - cc_other[ivt]) * fm_other[ivt]
            m_leafn_xfer_to_litter_fire[p]             = leafn_xfer[p] * f *
                (1.0 - cc_other[ivt]) * fm_other[ivt]
            m_livestemn_to_litter_fire[p]              = livestemn[p] * f *
                (1.0 - cc_lstem[ivt]) * fm_droot[ivt]
            m_livestemn_storage_to_litter_fire[p]      = livestemn_storage[p] * f *
                (1.0 - cc_other[ivt]) * fm_other[ivt]
            m_livestemn_xfer_to_litter_fire[p]         = livestemn_xfer[p] * f *
                (1.0 - cc_other[ivt]) * fm_other[ivt]
            m_livestemn_to_deadstemn_fire[p]           = livestemn[p] * f *
                (1.0 - cc_lstem[ivt]) * (fm_lstem[ivt] - fm_droot[ivt])
            m_deadstemn_to_litter_fire[p]              = deadstemn[p] * f * m *
                (1.0 - cc_dstem[ivt]) * fm_droot[ivt]
            m_deadstemn_storage_to_litter_fire[p]      = deadstemn_storage[p] * f *
                (1.0 - cc_other[ivt]) * fm_other[ivt]
            m_deadstemn_xfer_to_litter_fire[p]         = deadstemn_xfer[p] * f *
                (1.0 - cc_other[ivt]) * fm_other[ivt]
            m_frootn_to_litter_fire[p]                 = frootn[p] * f *
                fm_root[ivt]
            m_frootn_storage_to_litter_fire[p]         = frootn_storage[p] * f *
                (1.0 - cc_other[ivt]) * fm_other[ivt]
            m_frootn_xfer_to_litter_fire[p]            = frootn_xfer[p] * f *
                (1.0 - cc_other[ivt]) * fm_other[ivt]
            m_livecrootn_to_litter_fire[p]             = livecrootn[p] * f *
                fm_droot[ivt]
            m_livecrootn_storage_to_litter_fire[p]     = livecrootn_storage[p] * f *
                (1.0 - cc_other[ivt]) * fm_other[ivt]
            m_livecrootn_xfer_to_litter_fire[p]        = livecrootn_xfer[p] * f *
                (1.0 - cc_other[ivt]) * fm_other[ivt]
            m_livecrootn_to_deadcrootn_fire[p]         = livecrootn[p] * f *
                (fm_lroot[ivt] - fm_droot[ivt])
            m_deadcrootn_to_litter_fire[p]             = deadcrootn[p] * f * m *
                fm_droot[ivt]
            m_deadcrootn_storage_to_litter_fire[p]     = deadcrootn_storage[p] * f *
                (1.0 - cc_other[ivt]) * fm_other[ivt]
            m_deadcrootn_xfer_to_litter_fire[p]        = deadcrootn_xfer[p] * f *
                (1.0 - cc_other[ivt]) * fm_other[ivt]
            m_retransn_to_litter_fire[p]               = retransn[p] * f *
                (1.0 - cc_other[ivt]) * fm_other[ivt]
        end

        if use_cndv
            if woody[ivt] == 1.0
                if livestemc[p] + deadstemc[p] > 0.0
                    nind[p] = nind[p] * (1.0 - 1.0 * fm_droot[ivt] * f)
                else
                    nind[p] = 0.0
                end
            end
            leafcmax[p] = max(leafc[p] - m_leafc_to_fire[p] * dt, leafcmax[p])
            if ivt == noveg
                leafcmax[p] = 0.0
            end
        end
    end
end

fireb_cnfire_patchflux!(mask_actfirep, args...) =
    _launch!(_fireb_cnfire_patchflux_kernel!, mask_actfirep, args...)

# Section (2): per-patch litter/CWD PATCH->COLUMN SCATTER. Each masked patch
# accumulates (over decomp layers j and litter pools i) into per-column arrays
# via _scatter_add!. The column accumulators are NOT zeroed here — the host loop
# also accumulates into pre-zeroed arrays (no zero loop). The host's outer
# ordering is (j outer, p inner); since each (c,j[,i]) cell is accumulated once
# per patch, reordering to (p outer, j inner) keeps the same per-cell sum; on the
# KA CPU backend _scatter_add! iterates patches in p-order so the running sum is
# byte-identical to the host.
@kernel function _fireb_cnfire_litterscatter_kernel!(
    fire_mortality_c_to_cwdc, fire_mortality_n_to_cwdn,
    m_c_to_litr_fire, m_n_to_litr_fire,
    @Const(mask_soilp), @Const(itype), @Const(column), @Const(wtcol),
    @Const(stem_prof), @Const(croot_prof), @Const(leaf_prof), @Const(froot_prof),
    @Const(lf_f), @Const(fr_f),
    @Const(m_deadstemc_to_litter_fire), @Const(m_deadcrootc_to_litter_fire),
    @Const(m_deadstemn_to_litter_fire), @Const(m_deadcrootn_to_litter_fire),
    @Const(m_livestemc_to_litter_fire), @Const(m_livecrootc_to_litter_fire),
    @Const(m_livestemn_to_litter_fire), @Const(m_livecrootn_to_litter_fire),
    @Const(m_leafc_to_litter_fire), @Const(m_leafc_storage_to_litter_fire),
    @Const(m_leafc_xfer_to_litter_fire), @Const(m_gresp_storage_to_litter_fire),
    @Const(m_gresp_xfer_to_litter_fire),
    @Const(m_frootc_to_litter_fire), @Const(m_frootc_storage_to_litter_fire),
    @Const(m_frootc_xfer_to_litter_fire),
    @Const(m_livestemc_storage_to_litter_fire), @Const(m_livestemc_xfer_to_litter_fire),
    @Const(m_deadstemc_storage_to_litter_fire), @Const(m_deadstemc_xfer_to_litter_fire),
    @Const(m_livecrootc_storage_to_litter_fire), @Const(m_livecrootc_xfer_to_litter_fire),
    @Const(m_deadcrootc_storage_to_litter_fire), @Const(m_deadcrootc_xfer_to_litter_fire),
    @Const(m_leafn_to_litter_fire), @Const(m_leafn_storage_to_litter_fire),
    @Const(m_leafn_xfer_to_litter_fire), @Const(m_retransn_to_litter_fire),
    @Const(m_frootn_to_litter_fire), @Const(m_frootn_storage_to_litter_fire),
    @Const(m_frootn_xfer_to_litter_fire),
    @Const(m_livestemn_storage_to_litter_fire), @Const(m_livestemn_xfer_to_litter_fire),
    @Const(m_deadstemn_storage_to_litter_fire), @Const(m_deadstemn_xfer_to_litter_fire),
    @Const(m_livecrootn_storage_to_litter_fire), @Const(m_livecrootn_xfer_to_litter_fire),
    @Const(m_deadcrootn_storage_to_litter_fire), @Const(m_deadcrootn_xfer_to_litter_fire),
    nlevdecomp::Int, i_met_lit::Int, i_litr_max::Int
)
    p = @index(Global)
    @inbounds if mask_soilp[p]
        c = column[p]
        ivt = itype[p] + 1
        w = wtcol[p]
        for j in 1:nlevdecomp
            _scatter_add!(fire_mortality_c_to_cwdc, c, j,
                m_deadstemc_to_litter_fire[p] * w * stem_prof[p, j])
            _scatter_add!(fire_mortality_c_to_cwdc, c, j,
                m_deadcrootc_to_litter_fire[p] * w * croot_prof[p, j])
            _scatter_add!(fire_mortality_n_to_cwdn, c, j,
                m_deadstemn_to_litter_fire[p] * w * stem_prof[p, j])
            _scatter_add!(fire_mortality_n_to_cwdn, c, j,
                m_deadcrootn_to_litter_fire[p] * w * croot_prof[p, j])

            _scatter_add!(fire_mortality_c_to_cwdc, c, j,
                m_livestemc_to_litter_fire[p] * w * stem_prof[p, j])
            _scatter_add!(fire_mortality_c_to_cwdc, c, j,
                m_livecrootc_to_litter_fire[p] * w * croot_prof[p, j])
            _scatter_add!(fire_mortality_n_to_cwdn, c, j,
                m_livestemn_to_litter_fire[p] * w * stem_prof[p, j])
            _scatter_add!(fire_mortality_n_to_cwdn, c, j,
                m_livecrootn_to_litter_fire[p] * w * croot_prof[p, j])

            _scatter_add!(m_c_to_litr_fire, c, j, i_met_lit,
                ((m_leafc_to_litter_fire[p] * lf_f[ivt, i_met_lit] +
                  m_leafc_storage_to_litter_fire[p] +
                  m_leafc_xfer_to_litter_fire[p] +
                  m_gresp_storage_to_litter_fire[p] +
                  m_gresp_xfer_to_litter_fire[p]) * leaf_prof[p, j] +
                 (m_frootc_to_litter_fire[p] * fr_f[ivt, i_met_lit] +
                  m_frootc_storage_to_litter_fire[p] +
                  m_frootc_xfer_to_litter_fire[p]) * froot_prof[p, j] +
                 (m_livestemc_storage_to_litter_fire[p] +
                  m_livestemc_xfer_to_litter_fire[p] +
                  m_deadstemc_storage_to_litter_fire[p] +
                  m_deadstemc_xfer_to_litter_fire[p]) * stem_prof[p, j] +
                 (m_livecrootc_storage_to_litter_fire[p] +
                  m_livecrootc_xfer_to_litter_fire[p] +
                  m_deadcrootc_storage_to_litter_fire[p] +
                  m_deadcrootc_xfer_to_litter_fire[p]) * croot_prof[p, j]) * w)

            for i in (i_met_lit+1):i_litr_max
                _scatter_add!(m_c_to_litr_fire, c, j, i,
                    (m_leafc_to_litter_fire[p] * lf_f[ivt, i] * leaf_prof[p, j] +
                     m_frootc_to_litter_fire[p] * fr_f[ivt, i] * froot_prof[p, j]) * w)
            end

            _scatter_add!(m_n_to_litr_fire, c, j, i_met_lit,
                ((m_leafn_to_litter_fire[p] * lf_f[ivt, i_met_lit] +
                  m_leafn_storage_to_litter_fire[p] +
                  m_leafn_xfer_to_litter_fire[p] +
                  m_retransn_to_litter_fire[p]) * leaf_prof[p, j] +
                 (m_frootn_to_litter_fire[p] * fr_f[ivt, i_met_lit] +
                  m_frootn_storage_to_litter_fire[p] +
                  m_frootn_xfer_to_litter_fire[p]) * froot_prof[p, j] +
                 (m_livestemn_storage_to_litter_fire[p] +
                  m_livestemn_xfer_to_litter_fire[p] +
                  m_deadstemn_storage_to_litter_fire[p] +
                  m_deadstemn_xfer_to_litter_fire[p]) * stem_prof[p, j] +
                 (m_livecrootn_storage_to_litter_fire[p] +
                  m_livecrootn_xfer_to_litter_fire[p] +
                  m_deadcrootn_storage_to_litter_fire[p] +
                  m_deadcrootn_xfer_to_litter_fire[p]) * croot_prof[p, j]) * w)

            for i in (i_met_lit+1):i_litr_max
                _scatter_add!(m_n_to_litr_fire, c, j, i,
                    (m_leafn_to_litter_fire[p] * lf_f[ivt, i] * leaf_prof[p, j] +
                     m_frootn_to_litter_fire[p] * fr_f[ivt, i] * froot_prof[p, j]) * w)
            end
        end
    end
end

fireb_cnfire_litterscatter!(fire_mortality_c_to_cwdc, fire_mortality_n_to_cwdn,
                            m_c_to_litr_fire, m_n_to_litr_fire, mask_soilp, args...) =
    _launch!(_fireb_cnfire_litterscatter_kernel!, fire_mortality_c_to_cwdc,
             fire_mortality_n_to_cwdn, m_c_to_litr_fire, m_n_to_litr_fire,
             mask_soilp, args...; ndrange = length(mask_soilp))

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
    @inbounds if mask_soilc[c]
        f = farea_burned[c]
        if f != 0 || f != baf_crop[c]
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
    @inbounds if mask_soilc[c]
        lfc2[c] = 0.0
        if !is_newyear_start
            if trotr1_col[c] + trotr2_col[c] > 0.6 && dtrotr_col[c] > 0.0 &&
               lfc[c] > 0.0 && fbac1[c] == 0.0
                lfc2[c] = max(0.0, min(lfc[c], (farea_burned[c] - baf_crop[c] -
                    baf_peatf[c]) / 2.0 * dt)) / (dtrotr_col[c] * dayspyr * secspday / dt) / dt
                lfc[c]  = lfc[c] - max(0.0, min(lfc[c], (farea_burned[c] - baf_crop[c] -
                    baf_peatf[c]) * dt / 2.0))
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
