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
    h2osoi_vol_col::Matrix{Float64},
    nlevgrnd::Int;
    soil_suction_fn::Function = default_soil_suction
)
    btran2   = fire_data.btran2_patch
    smpso    = pftcon.smpso
    smpsc    = pftcon.smpsc
    watsat   = soilstate.watsat_col
    rootfr   = soilstate.rootfr_patch

    # Zero out non-exposed vegetation patches
    for p in bounds
        mask_noexposedveg[p] || continue
        btran2[p] = 0.0
    end

    # Initialize exposed veg patches to zero
    for p in bounds
        mask_exposedveg[p] || continue
        btran2[p] = 0.0
    end

    # Accumulate root-weighted soil wetness across layers
    for j in 1:nlevgrnd
        for p in bounds
            mask_exposedveg[p] || continue
            c = patch.column[p]
            s_node = max(h2osoi_vol_col[c, j] / watsat[c, j], 0.01)

            smp_node_lf = soil_suction_fn(c, j, s_node, soilstate)

            ivt = patch.itype[p]
            smp_node_lf = max(smpsc[ivt], smp_node_lf)
            btran2[p] = btran2[p] + rootfr[p, j] * max(0.0,
                min((smp_node_lf - smpsc[ivt]) / (smpso[ivt] - smpsc[ivt]), 1.0))
        end
    end

    # Clamp to [0, 1]
    for p in bounds
        mask_exposedveg[p] || continue
        if btran2[p] > 1.0
            btran2[p] = 1.0
        end
    end

    return nothing
end

"""
    default_soil_suction(c, j, s_node, soilstate)

Default soil suction function placeholder. Returns suction based on
Clapp-Hornberger parameterization. In practice, this should be replaced
with the appropriate soil water retention curve method.
"""
function default_soil_suction(c::Int, j::Int, s_node::Float64, soilstate::SoilStateData)
    # Clapp-Hornberger: smp = sucsat * s^(-bsw)
    sucsat = soilstate.sucsat_col[c, j]
    bsw    = soilstate.bsw_col[c, j]
    return -sucsat * s_node^(-bsw)
end

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
function cnfire_calc_fire_root_wetness_li2021!(
    fire_data::CNFireBaseData,
    mask_exposedveg::BitVector,
    mask_noexposedveg::BitVector,
    bounds::UnitRange{Int},
    patch::PatchData,
    soilstate::SoilStateData,
    h2osoi_vol_col::Matrix{Float64},
    nlevgrnd::Int
)
    btran2 = fire_data.btran2_patch
    watsat = soilstate.watsat_col
    rootfr = soilstate.rootfr_patch

    # Zero out non-exposed vegetation patches
    for p in bounds
        mask_noexposedveg[p] || continue
        btran2[p] = 0.0
    end

    # Initialize exposed veg patches to zero
    for p in bounds
        mask_exposedveg[p] || continue
        btran2[p] = 0.0
    end

    # Accumulate root-weighted s_node across layers
    for j in 1:nlevgrnd
        for p in bounds
            mask_exposedveg[p] || continue
            c = patch.column[p]
            s_node = max(h2osoi_vol_col[c, j] / watsat[c, j], 0.01)
            btran2[p] = btran2[p] + rootfr[p, j] * s_node
        end
    end

    # Clamp to [0, 1]
    for p in bounds
        mask_exposedveg[p] || continue
        if btran2[p] > 1.0
            btran2[p] = 1.0
        end
    end

    return nothing
end

# ---------------------------------------------------------------------------
# cnfire_fluxes! — Main fire flux calculation
# ---------------------------------------------------------------------------

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
    leaf_prof_patch::Matrix{Float64},
    froot_prof_patch::Matrix{Float64},
    croot_prof_patch::Matrix{Float64},
    stem_prof_patch::Matrix{Float64},
    totsomc_col::Vector{Float64},
    decomp_cpools_vr_col::Array{Float64,3},
    decomp_npools_vr_col::Array{Float64,3},
    somc_fire_col::Vector{Float64};
    dt::Float64 = 1800.0,
    dayspyr::Float64 = 365.0,
    nlevdecomp::Int = 1,
    ndecomp_pools::Int = 7,
    i_met_lit::Int = 1,
    i_litr_max::Int = 3,
    transient_landcover::Bool = false,
    use_cndv::Bool = false,
    use_matrixcn::Bool = false,
    spinup_factor_deadwood::Float64 = SPINUP_FACTOR_DEADWOOD_DEFAULT,
    secspday::Float64 = SECSPDAY,
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
    # Patch loop — fire emissions and mortality fluxes
    # ===================================================================
    m = spinup_factor_deadwood

    for p in bounds_p
        mask_soilp[p] || continue
        c = patch.column[p]

        if patch.itype[p] < nc3crop && cropf_col[c] < 1.0
            # Non-crop (bare-soil and natural vegetation)
            if transient_landcover
                f = (fbac[c] - baf_crop[c]) / (1.0 - cropf_col[c])
            else
                f = (farea_burned[c] - baf_crop[c]) / (1.0 - cropf_col[c])
            end
        else
            # Crops
            if cropf_col[c] > 0.0
                f = baf_crop[c] / cropf_col[c]
            else
                f = 0.0
            end
        end

        # Track active fire patches
        if f != 0
            mask_actfirep[p] = true
        end

        ivt = patch.itype[p]

        # Biomass burning — carbon emission fluxes
        m_gresp_storage_to_fire[p]       = gresp_storage[p]      * f * cc_other[ivt]
        m_gresp_xfer_to_fire[p]          = gresp_xfer[p]         * f * cc_other[ivt]

        if !use_matrixcn
            # Non-matrix version of fire emissions
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

            # Nitrogen emission fluxes
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
        # (Matrix CN branch omitted — only non-matrix version ported for Phase 1)

        # Mortality due to fire — carbon pools
        if !use_matrixcn
            m_leafc_to_litter_fire[p]                   = leafc[p] * f *
                (1.0 - cc_leaf[ivt]) * fm_leaf[ivt]
            m_leafc_storage_to_litter_fire[p]           = leafc_storage[p] * f *
                (1.0 - cc_other[ivt]) * fm_other[ivt]
            m_leafc_xfer_to_litter_fire[p]              = leafc_xfer[p] * f *
                (1.0 - cc_other[ivt]) * fm_other[ivt]
            # NOTE: fm_droot used for livestem transport to litter (see Fortran bug 2516)
            m_livestemc_to_litter_fire[p]               = livestemc[p] * f *
                (1.0 - cc_lstem[ivt]) * fm_droot[ivt]
            m_livestemc_storage_to_litter_fire[p]       = livestemc_storage[p] * f *
                (1.0 - cc_other[ivt]) * fm_other[ivt]
            m_livestemc_xfer_to_litter_fire[p]          = livestemc_xfer[p] * f *
                (1.0 - cc_other[ivt]) * fm_other[ivt]
            # NOTE: fm_droot used for fraction of plant-tissue mortality for deadstem (bug 2516)
            m_livestemc_to_deadstemc_fire[p]            = livestemc[p] * f *
                (1.0 - cc_lstem[ivt]) * (fm_lstem[ivt] - fm_droot[ivt])
            # NOTE: fm_droot used for deadstem transport to litter (bug 2516)
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
            # NOTE: fm_droot used for livecroot transport to litter (bug 2516)
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

            # Nitrogen mortality pools
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

        # CNDV updates
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

    end  # end of patches loop

    # ===================================================================
    # Fire-induced transfer of C and N pools to litter and CWD
    # ===================================================================
    for j in 1:nlevdecomp
        for p in bounds_p
            mask_soilp[p] || continue
            c = patch.column[p]
            ivt = patch.itype[p]

            fire_mortality_c_to_cwdc[c, j] = fire_mortality_c_to_cwdc[c, j] +
                m_deadstemc_to_litter_fire[p] * patch.wtcol[p] * stem_prof[p, j]
            fire_mortality_c_to_cwdc[c, j] = fire_mortality_c_to_cwdc[c, j] +
                m_deadcrootc_to_litter_fire[p] * patch.wtcol[p] * croot_prof[p, j]
            fire_mortality_n_to_cwdn[c, j] = fire_mortality_n_to_cwdn[c, j] +
                m_deadstemn_to_litter_fire[p] * patch.wtcol[p] * stem_prof[p, j]
            fire_mortality_n_to_cwdn[c, j] = fire_mortality_n_to_cwdn[c, j] +
                m_deadcrootn_to_litter_fire[p] * patch.wtcol[p] * croot_prof[p, j]

            fire_mortality_c_to_cwdc[c, j] = fire_mortality_c_to_cwdc[c, j] +
                m_livestemc_to_litter_fire[p] * patch.wtcol[p] * stem_prof[p, j]
            fire_mortality_c_to_cwdc[c, j] = fire_mortality_c_to_cwdc[c, j] +
                m_livecrootc_to_litter_fire[p] * patch.wtcol[p] * croot_prof[p, j]
            fire_mortality_n_to_cwdn[c, j] = fire_mortality_n_to_cwdn[c, j] +
                m_livestemn_to_litter_fire[p] * patch.wtcol[p] * stem_prof[p, j]
            fire_mortality_n_to_cwdn[c, j] = fire_mortality_n_to_cwdn[c, j] +
                m_livecrootn_to_litter_fire[p] * patch.wtcol[p] * croot_prof[p, j]

            # Litter: metabolic litter (i_met_lit) is handled specially
            m_c_to_litr_fire[c, j, i_met_lit] =
                m_c_to_litr_fire[c, j, i_met_lit] +
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
                  m_deadcrootc_xfer_to_litter_fire[p]) * croot_prof[p, j]) * patch.wtcol[p]

            # Other litter types
            for i in (i_met_lit+1):i_litr_max
                m_c_to_litr_fire[c, j, i] = m_c_to_litr_fire[c, j, i] +
                    (m_leafc_to_litter_fire[p] * lf_f[ivt, i] * leaf_prof[p, j] +
                     m_frootc_to_litter_fire[p] * fr_f[ivt, i] * froot_prof[p, j]) * patch.wtcol[p]
            end

            # Nitrogen to litter: metabolic
            m_n_to_litr_fire[c, j, i_met_lit] =
                m_n_to_litr_fire[c, j, i_met_lit] +
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
                  m_deadcrootn_xfer_to_litter_fire[p]) * croot_prof[p, j]) * patch.wtcol[p]

            # Other litter N types
            for i in (i_met_lit+1):i_litr_max
                m_n_to_litr_fire[c, j, i] =
                    m_n_to_litr_fire[c, j, i] +
                    (m_leafn_to_litter_fire[p] * lf_f[ivt, i] * leaf_prof[p, j] +
                     m_frootn_to_litter_fire[p] * fr_f[ivt, i] * froot_prof[p, j]) * patch.wtcol[p]
            end
        end
    end

    # ===================================================================
    # Vertically-resolved decomposing C/N fire loss — column loop
    # ===================================================================
    for c in bounds_c
        mask_soilc[c] || continue

        f = farea_burned[c]

        if f != 0 || f != baf_crop[c]
            mask_actfirec[c] = true
        end

        for j in 1:nlevdecomp
            # Carbon fluxes
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

            # Nitrogen fluxes
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

    # ===================================================================
    # Carbon loss due to deforestation fires
    # ===================================================================
    if transient_landcover
        for c in bounds_c
            mask_soilc[c] || continue
            lfc2[c] = 0.0
            if !(kmo == 1 && kda == 1 && mcsec == 0)
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

    # ===================================================================
    # Carbon loss due to peat fires
    # ===================================================================
    for c in bounds_c
        mask_soilc[c] || continue
        g = col.gridcell[c]
        if grc.latdeg[g] < cnfire_const.borealat
            somc_fire[c] = totsomc[c] * baf_peatf[c] * 6.0 / 33.9
        else
            somc_fire[c] = baf_peatf[c] * 2.2e3
        end
    end

    return (mask_actfirec, mask_actfirep)
end
