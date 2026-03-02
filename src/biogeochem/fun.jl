# =============================================================================
# FUN (Fixation and Uptake of Nitrogen) Model
# =============================================================================
# Ported from CNFUNMod.F90
# The FUN model developed by Fisher et al. 2010 and Brzostek et al. 2014.
# Coded by Mingjie Shi 2015; logic altered by Rosie Fisher Oct 2015.
# Critical outputs: sminn_to_plant_fun and npp_Nuptake (N available for growth
# and C spent obtaining it).

# ---- Module-level constants ------------------------------------------------
const COST_METHOD    = 2       # new N uptake resistance method
const NSTP           = 2       # number of calculation substeps (ECM, AM)
const NCOST6         = 6       # number of N transport pathways

const ICOST_FIX          = 1   # process number for fixation
const ICOST_RETRANS      = 2   # process number for retranslocation
const ICOST_ACTIVE_NO3   = 3   # process number for mycorrhizal NO3 uptake
const ICOST_ACTIVE_NH4   = 4   # process number for mycorrhizal NH4 uptake
const ICOST_NONMYC_NO3   = 5   # process number for non-myc NO3 uptake
const ICOST_NONMYC_NH4   = 6   # process number for non-myc NH4 uptake

const BIG_COST = 1.0e9        # arbitrarily large cost (gC/gN)

const PLANTS_ARE_FIXING = 1
const PLANTS_NOT_FIXING = 2
const ECM_STEP = 1
const AM_STEP  = 2

# ---- FUN parameters struct ------------------------------------------------
Base.@kwdef mutable struct FUNParams
    ndays_off::Float64 = 21.0   # number of days to complete leaf offset
end

# ---- PFT constants needed by FUN ------------------------------------------
Base.@kwdef mutable struct PftConFUN
    leafcn          ::Vector{Float64} = Float64[]   # leaf C:N (gC/gN)
    season_decid    ::Vector{Float64} = Float64[]   # binary flag seasonal-deciduous
    stress_decid    ::Vector{Float64} = Float64[]   # binary flag stress-deciduous
    a_fix           ::Vector{Float64} = Float64[]   # BNF parameter
    b_fix           ::Vector{Float64} = Float64[]   # BNF parameter
    c_fix           ::Vector{Float64} = Float64[]   # BNF parameter
    s_fix           ::Vector{Float64} = Float64[]   # BNF parameter
    akc_active      ::Vector{Float64} = Float64[]   # AM mycorrhizal uptake kc
    akn_active      ::Vector{Float64} = Float64[]   # AM mycorrhizal uptake kn
    ekc_active      ::Vector{Float64} = Float64[]   # ECM mycorrhizal uptake kc
    ekn_active      ::Vector{Float64} = Float64[]   # ECM mycorrhizal uptake kn
    kc_nonmyc       ::Vector{Float64} = Float64[]   # non-mycorrhizal uptake kc
    kn_nonmyc       ::Vector{Float64} = Float64[]   # non-mycorrhizal uptake kn
    perecm          ::Vector{Float64} = Float64[]   # fraction ECM-associated
    grperc          ::Vector{Float64} = Float64[]   # growth respiration percentage
    fun_cn_flex_a   ::Vector{Float64} = Float64[]   # flexCN parameter a
    fun_cn_flex_b   ::Vector{Float64} = Float64[]   # flexCN parameter b
    fun_cn_flex_c   ::Vector{Float64} = Float64[]   # flexCN parameter c
    FUN_fracfixers  ::Vector{Float64} = Float64[]   # fraction of C for fixation
    c3psn           ::Vector{Float64} = Float64[]   # C3 photosynthesis flag
end

# =============================================================================
#  Cost functions (pure functions, GPU-friendly)
# =============================================================================

"""
    fun_cost_fix(fixer, a_fix, b_fix, c_fix, big_cost, crootfr, s_fix, tc_soisno)

Cost of fixing N by nodules (gC/gN). Returns `big_cost` for non-fixers or
layers with negligible root fraction.
"""
function fun_cost_fix(fixer::Int, a_fix::Float64, b_fix::Float64,
                      c_fix::Float64, big_cost::Float64,
                      crootfr::Float64, s_fix::Float64,
                      tc_soisno::Float64)
    if fixer == 1 && crootfr > 1.0e-6
        return (-1.0 * s_fix) * 1.0 /
               (1.25 * exp(a_fix + b_fix * tc_soisno * (1.0 - 0.5 * tc_soisno / c_fix)))
    else
        return big_cost
    end
end

"""
    fun_cost_active(sminn_layer, big_cost, kc_active, kn_active, rootc_dens, crootfr, smallValue)

Cost of mycorrhizal active uptake of N (gC/gN).
"""
function fun_cost_active(sminn_layer::Float64, big_cost::Float64,
                         kc_active::Float64, kn_active::Float64,
                         rootc_dens::Float64, crootfr::Float64,
                         smallValue::Float64)
    if rootc_dens > 1.0e-6 && sminn_layer > smallValue
        return kn_active / sminn_layer + kc_active / rootc_dens
    else
        return big_cost
    end
end

"""
    fun_cost_nonmyc(sminn_layer, big_cost, kc_nonmyc, kn_nonmyc, rootc_dens, crootfr, smallValue)

Cost of non-mycorrhizal uptake of N (gC/gN).
"""
function fun_cost_nonmyc(sminn_layer::Float64, big_cost::Float64,
                         kc_nonmyc::Float64, kn_nonmyc::Float64,
                         rootc_dens::Float64, crootfr::Float64,
                         smallValue::Float64)
    if rootc_dens > 1.0e-6 && sminn_layer > smallValue
        return kn_nonmyc / sminn_layer + kc_nonmyc / rootc_dens
    else
        return big_cost
    end
end

# =============================================================================
#  Retranslocation subroutine
# =============================================================================

"""
    fun_retranslocation!(p, dt, npp_to_spend, total_falling_leaf_c,
        total_falling_leaf_n, total_n_resistance, target_leafcn, grperc, plantCN)

Calculate the amount of N absorbed and C spent during retranslocation.
Returns (total_c_spent_retrans, total_c_accounted_retrans, free_n_retrans, paid_for_n_retrans).
"""
function fun_retranslocation(dt::Float64, npp_to_spend::Float64,
                             total_falling_leaf_c::Float64,
                             total_falling_leaf_n::Float64,
                             total_n_resistance::Float64,
                             target_leafcn::Float64,
                             grperc::Float64,
                             plantCN::Float64)
    # Initialize total fluxes
    total_c_spent_retrans     = 0.0
    total_c_accounted_retrans = 0.0
    c_accounted_retrans       = 0.0
    paid_for_n_retrans        = 0.0
    npp_to_spend_temp         = npp_to_spend

    # Initial C and N pools in falling leaves
    falling_leaf_c = total_falling_leaf_c
    falling_leaf_n = total_falling_leaf_n

    # Parameters
    max_falling_leaf_cn = target_leafcn * 3.0
    min_falling_leaf_cn = target_leafcn * 1.5
    cost_escalation     = 1.3

    # Free uptake
    free_n_retrans = max(falling_leaf_n - (falling_leaf_c / min_falling_leaf_cn), 0.0)
    falling_leaf_n = falling_leaf_n - free_n_retrans

    # Initial CN ratio and costs
    falling_leaf_cn   = falling_leaf_c / falling_leaf_n
    kresorb           = 1.0 / target_leafcn
    cost_retrans_temp = kresorb / ((1.0 / falling_leaf_cn)^1.3)

    # Iteration loop
    iter = 0
    exitloop = false
    while !exitloop && cost_retrans_temp < total_n_resistance &&
          falling_leaf_n >= 0.0 && npp_to_spend > 0.0

        # Spend some C on removing N
        c_spent_retrans = cost_retrans_temp * (falling_leaf_n - falling_leaf_c /
                          (falling_leaf_cn + 1.0))
        # don't spend more C than you have
        c_spent_retrans = min(npp_to_spend_temp, c_spent_retrans)
        # N extracted
        leaf_n_ext = c_spent_retrans / cost_retrans_temp
        # Do not empty N pool
        leaf_n_ext = min(falling_leaf_n, leaf_n_ext)
        # C accounted for by N uptake
        c_accounted_retrans = leaf_n_ext * plantCN * (1.0 + grperc)

        # Update leafCN, recalculate costs
        falling_leaf_n = falling_leaf_n - leaf_n_ext
        if falling_leaf_n > 0.0
            falling_leaf_cn   = falling_leaf_c / falling_leaf_n
            cost_retrans_temp = kresorb / ((1.0 / falling_leaf_cn)^1.3)
        else
            exitloop = true
        end

        # Accumulate total fluxes
        total_c_spent_retrans     = total_c_spent_retrans + c_spent_retrans
        total_c_accounted_retrans = total_c_accounted_retrans + c_accounted_retrans
        paid_for_n_retrans        = paid_for_n_retrans + leaf_n_ext
        npp_to_spend_temp         = npp_to_spend_temp - c_spent_retrans - c_accounted_retrans
        iter += 1

        # ran out of C or N
        if npp_to_spend_temp <= 0.0
            exitloop = true
            total_c_spent_retrans = total_c_spent_retrans + npp_to_spend_temp
        end
        # leaf CN is too high
        if falling_leaf_cn >= max_falling_leaf_cn
            exitloop = true
        end
        # safety check
        if iter >= 150
            exitloop = true
        end
    end

    return (total_c_spent_retrans=total_c_spent_retrans,
            total_c_accounted_retrans=total_c_accounted_retrans,
            free_n_retrans=free_n_retrans,
            paid_for_n_retrans=paid_for_n_retrans)
end

# =============================================================================
#  CNFUNInit — initialization called at start of year / FUN period
# =============================================================================

"""
    cnfun_init!(mask_soilp, bounds, fun_params, pftcon, patch,
                cnveg_state, cnveg_cs, cnveg_ns;
                dt, nstep, dayspyr, npcropmin)

Initialize FUN variables at the start of a new FUN period.
"""
function cnfun_init!(mask_soilp::BitVector, bounds::UnitRange{Int},
                     fun_params::FUNParams,
                     pftcon::PftConFUN,
                     patch::PatchData,
                     cnveg_state::CNVegStateData,
                     cnveg_cs::CNVegCarbonStateData,
                     cnveg_ns::CNVegNitrogenStateData;
                     dt::Float64=1800.0,
                     nstep::Int=1,
                     dayspyr::Float64=365.0,
                     npcropmin::Int=17)

    timestep_fun = SECSPDAY * FUN_PERIOD
    nstep_fun    = round(Int, SECSPDAY * dayspyr / dt)

    if nstep_fun > 0 && mod(nstep, nstep_fun) == 0
        for p in bounds
            mask_soilp[p] || continue
            ivt = patch.itype[p]
            cnveg_state.leafcn_offset_patch[p]         = pftcon.leafcn[ivt]
            cnveg_cs.storage_cdemand_patch[p]           = 0.0
            cnveg_ns.storage_ndemand_patch[p]           = 0.0
            cnveg_ns.leafn_storage_xfer_acc_patch[p]    = 0.0
            cnveg_cs.leafc_storage_xfer_acc_patch[p]    = 0.0
        end
    end
end

# =============================================================================
#  CNFUN — main FUN subroutine
# =============================================================================

"""
    cnfun!(mask_soilp, mask_soilc, bounds_p, bounds_c,
           fun_params, pftcon, patch, waterstate, temperature,
           soilstate, cnveg_state, cnveg_cs, cnveg_cf, cnveg_ns, cnveg_nf,
           soilbgc_nf, soilbgc_cf, canopystate, soilbgc_ns;
           dt, nlevdecomp, dzsoi_decomp_vals, use_flexiblecn, use_matrixcn,
           npcropmin, smallValue, spval)

Main FUN calculation. Computes N uptake via fixation, retranslocation,
mycorrhizal and non-mycorrhizal pathways.
"""
function cnfun!(mask_soilp::BitVector, mask_soilc::BitVector,
                bounds_p::UnitRange{Int}, bounds_c::UnitRange{Int},
                fun_params::FUNParams,
                pftcon::PftConFUN,
                patch::PatchData,
                waterstate::WaterStateData,
                temperature::TemperatureData,
                soilstate::SoilStateData,
                cnveg_state::CNVegStateData,
                cnveg_cs::CNVegCarbonStateData,
                cnveg_cf::CNVegCarbonFluxData,
                cnveg_ns::CNVegNitrogenStateData,
                cnveg_nf::CNVegNitrogenFluxData,
                soilbgc_nf::SoilBiogeochemNitrogenFluxData,
                soilbgc_cf::SoilBiogeochemCarbonFluxData,
                canopystate::CanopyStateData,
                soilbgc_ns::SoilBiogeochemNitrogenStateData;
                dt::Float64=1800.0,
                nlevdecomp::Int=1,
                dzsoi_decomp_vals::Vector{Float64}=Float64[],
                use_flexiblecn::Bool=false,
                use_matrixcn::Bool=false,
                npcropmin::Int=17,
                smallValue::Float64=SMALLVALUE,
                spval::Float64=SPVAL)

    # Aliases for state/flux fields
    ivt                    = patch.itype
    leafcn                 = pftcon.leafcn
    season_decid           = pftcon.season_decid
    stress_decid           = pftcon.stress_decid
    a_fix_pft              = pftcon.a_fix
    b_fix_pft              = pftcon.b_fix
    c_fix_pft              = pftcon.c_fix
    s_fix_pft              = pftcon.s_fix
    akc_active_pft         = pftcon.akc_active
    akn_active_pft         = pftcon.akn_active
    ekc_active_pft         = pftcon.ekc_active
    ekn_active_pft         = pftcon.ekn_active
    kc_nonmyc_pft          = pftcon.kc_nonmyc
    kn_nonmyc_pft          = pftcon.kn_nonmyc
    perecm                 = pftcon.perecm
    grperc_pft             = pftcon.grperc
    fun_cn_flex_a_pft      = pftcon.fun_cn_flex_a
    fun_cn_flex_b_pft      = pftcon.fun_cn_flex_b
    fun_cn_flex_c_pft      = pftcon.fun_cn_flex_c
    FUN_fracfixers_pft     = pftcon.FUN_fracfixers
    c3psn_pft              = pftcon.c3psn

    ndays_off = fun_params.ndays_off
    steppday  = 48.0
    local_use_flexibleCN = use_flexiblecn

    np = length(bounds_p)
    nc = length(bounds_c)

    # ---- Allocate local arrays ----
    rootc_dens       = zeros(Float64, last(bounds_p), nlevdecomp)
    rootC            = zeros(Float64, last(bounds_p))
    permyc           = zeros(Float64, last(bounds_p), NSTP)
    kc_active        = zeros(Float64, last(bounds_p), NSTP)
    kn_active        = zeros(Float64, last(bounds_p), NSTP)
    availc_pool      = zeros(Float64, last(bounds_p))
    plantN           = zeros(Float64, last(bounds_p))
    plant_ndemand_pool       = zeros(Float64, last(bounds_p))
    plant_ndemand_pool_step  = zeros(Float64, last(bounds_p), NSTP)
    leafn_step               = zeros(Float64, last(bounds_p), NSTP)
    leafn_retrans_step       = zeros(Float64, last(bounds_p), NSTP)
    litterfall_n             = zeros(Float64, last(bounds_p))
    litterfall_n_step        = zeros(Float64, last(bounds_p), NSTP)
    litterfall_c_step        = zeros(Float64, last(bounds_p), NSTP)
    tc_soisno_local          = zeros(Float64, last(bounds_c), nlevdecomp)
    npp_remaining            = zeros(Float64, last(bounds_p), NSTP)
    n_passive_step           = zeros(Float64, last(bounds_p), NSTP)
    n_passive_acc            = zeros(Float64, last(bounds_p))

    cost_retran_local  = zeros(Float64, last(bounds_p), nlevdecomp)
    cost_fix_local     = zeros(Float64, last(bounds_p), nlevdecomp)
    cost_active_no3    = zeros(Float64, last(bounds_p), nlevdecomp)
    cost_active_nh4    = zeros(Float64, last(bounds_p), nlevdecomp)
    cost_nonmyc_no3    = zeros(Float64, last(bounds_p), nlevdecomp)
    cost_nonmyc_nh4    = zeros(Float64, last(bounds_p), nlevdecomp)

    n_fix_acc              = zeros(Float64, last(bounds_p), NSTP)
    n_fix_acc_total        = zeros(Float64, last(bounds_p))
    npp_fix_acc            = zeros(Float64, last(bounds_p), NSTP)
    npp_fix_acc_total      = zeros(Float64, last(bounds_p))
    n_retrans_acc          = zeros(Float64, last(bounds_p), NSTP)
    n_retrans_acc_total    = zeros(Float64, last(bounds_p))
    free_nretrans_acc      = zeros(Float64, last(bounds_p), NSTP)
    npp_retrans_acc        = zeros(Float64, last(bounds_p), NSTP)
    npp_retrans_acc_total  = zeros(Float64, last(bounds_p))
    nt_uptake              = zeros(Float64, last(bounds_p), NSTP)
    npp_uptake             = zeros(Float64, last(bounds_p), NSTP)

    # NO3/NH4 arrays
    sminn_no3_conc       = zeros(Float64, last(bounds_c), nlevdecomp)
    sminn_no3_conc_step  = zeros(Float64, last(bounds_p), nlevdecomp, NSTP)
    sminn_no3_layer      = zeros(Float64, last(bounds_c), nlevdecomp)
    sminn_no3_layer_step = zeros(Float64, last(bounds_p), nlevdecomp, NSTP)
    sminn_no3_uptake     = zeros(Float64, last(bounds_p), nlevdecomp, NSTP)
    sminn_nh4_conc       = zeros(Float64, last(bounds_c), nlevdecomp)
    sminn_nh4_conc_step  = zeros(Float64, last(bounds_p), nlevdecomp, NSTP)
    sminn_nh4_layer      = zeros(Float64, last(bounds_c), nlevdecomp)
    sminn_nh4_layer_step = zeros(Float64, last(bounds_p), nlevdecomp, NSTP)
    sminn_nh4_uptake     = zeros(Float64, last(bounds_p), nlevdecomp, NSTP)

    n_active_no3_acc       = zeros(Float64, last(bounds_p), NSTP)
    n_active_nh4_acc       = zeros(Float64, last(bounds_p), NSTP)
    n_nonmyc_no3_acc       = zeros(Float64, last(bounds_p), NSTP)
    n_nonmyc_nh4_acc       = zeros(Float64, last(bounds_p), NSTP)
    n_active_no3_acc_total = zeros(Float64, last(bounds_p))
    n_active_nh4_acc_total = zeros(Float64, last(bounds_p))
    n_nonmyc_no3_acc_total = zeros(Float64, last(bounds_p))
    n_nonmyc_nh4_acc_total = zeros(Float64, last(bounds_p))
    npp_active_no3_acc       = zeros(Float64, last(bounds_p), NSTP)
    npp_active_nh4_acc       = zeros(Float64, last(bounds_p), NSTP)
    npp_nonmyc_no3_acc       = zeros(Float64, last(bounds_p), NSTP)
    npp_nonmyc_nh4_acc       = zeros(Float64, last(bounds_p), NSTP)
    npp_active_no3_acc_total = zeros(Float64, last(bounds_p))
    npp_active_nh4_acc_total = zeros(Float64, last(bounds_p))
    npp_nonmyc_no3_acc_total = zeros(Float64, last(bounds_p))
    npp_nonmyc_nh4_acc_total = zeros(Float64, last(bounds_p))

    n_am_no3_acc  = zeros(Float64, last(bounds_p))
    n_am_nh4_acc  = zeros(Float64, last(bounds_p))
    n_ecm_no3_acc = zeros(Float64, last(bounds_p))
    n_ecm_nh4_acc = zeros(Float64, last(bounds_p))
    n_am_no3_retrans  = zeros(Float64, last(bounds_p))
    n_am_nh4_retrans  = zeros(Float64, last(bounds_p))
    n_ecm_no3_retrans = zeros(Float64, last(bounds_p))
    n_ecm_nh4_retrans = zeros(Float64, last(bounds_p))

    n_active_no3_retrans       = zeros(Float64, last(bounds_p), NSTP)
    n_active_nh4_retrans       = zeros(Float64, last(bounds_p), NSTP)
    n_nonmyc_no3_retrans       = zeros(Float64, last(bounds_p), NSTP)
    n_nonmyc_nh4_retrans       = zeros(Float64, last(bounds_p), NSTP)
    n_active_no3_retrans_total = zeros(Float64, last(bounds_p))
    n_active_nh4_retrans_total = zeros(Float64, last(bounds_p))
    n_nonmyc_no3_retrans_total = zeros(Float64, last(bounds_p))
    n_nonmyc_nh4_retrans_total = zeros(Float64, last(bounds_p))
    npp_active_no3_retrans_total = zeros(Float64, last(bounds_p))
    npp_active_nh4_retrans_total = zeros(Float64, last(bounds_p))
    npp_nonmyc_no3_retrans_total = zeros(Float64, last(bounds_p))
    npp_nonmyc_nh4_retrans_total = zeros(Float64, last(bounds_p))

    n_passive_no3_vr = zeros(Float64, last(bounds_p), nlevdecomp)
    n_passive_nh4_vr = zeros(Float64, last(bounds_p), nlevdecomp)
    n_active_no3_vr  = zeros(Float64, last(bounds_p), nlevdecomp)
    n_nonmyc_no3_vr  = zeros(Float64, last(bounds_p), nlevdecomp)
    n_active_nh4_vr  = zeros(Float64, last(bounds_p), nlevdecomp)
    n_nonmyc_nh4_vr  = zeros(Float64, last(bounds_p), nlevdecomp)

    free_Nretrans_local = zeros(Float64, last(bounds_p))

    costNit = fill(BIG_COST, nlevdecomp, NCOST6)

    # Per-layer flux arrays (reused per-patch)
    npp_to_fixation     = zeros(Float64, nlevdecomp)
    npp_to_retrans      = zeros(Float64, nlevdecomp)
    npp_to_active_nh4   = zeros(Float64, nlevdecomp)
    npp_to_nonmyc_nh4   = zeros(Float64, nlevdecomp)
    npp_to_active_no3   = zeros(Float64, nlevdecomp)
    npp_to_nonmyc_no3   = zeros(Float64, nlevdecomp)

    npp_frac_to_fixation     = zeros(Float64, nlevdecomp)
    npp_frac_to_retrans      = zeros(Float64, nlevdecomp)
    npp_frac_to_active_nh4   = zeros(Float64, nlevdecomp)
    npp_frac_to_nonmyc_nh4   = zeros(Float64, nlevdecomp)
    npp_frac_to_active_no3   = zeros(Float64, nlevdecomp)
    npp_frac_to_nonmyc_no3   = zeros(Float64, nlevdecomp)

    n_exch_fixation     = zeros(Float64, nlevdecomp)
    n_exch_active_nh4   = zeros(Float64, nlevdecomp)
    n_exch_nonmyc_nh4   = zeros(Float64, nlevdecomp)
    n_exch_active_no3   = zeros(Float64, nlevdecomp)
    n_exch_nonmyc_no3   = zeros(Float64, nlevdecomp)

    n_from_fixation     = zeros(Float64, nlevdecomp)
    n_from_active_nh4   = zeros(Float64, nlevdecomp)
    n_from_nonmyc_nh4   = zeros(Float64, nlevdecomp)
    n_from_active_no3   = zeros(Float64, nlevdecomp)
    n_from_nonmyc_no3   = zeros(Float64, nlevdecomp)

    # ======================================================================
    # Phase 1: Pre-computation (across all patches)
    # ======================================================================

    # Compute litterfall_n, rootC, plantN, plantCN
    for p in bounds_p
        mask_soilp[p] || continue

        cnveg_cf.npp_burnedoff_patch[p] = 0.0

        litterfall_n[p] = (cnveg_cf.leafc_to_litter_fun_patch[p] /
                           cnveg_state.leafcn_offset_patch[p]) * dt
        rootC[p] = cnveg_cs.frootc_patch[p]

        plantN[p] = cnveg_ns.leafn_patch[p] + cnveg_ns.frootn_patch[p] +
                    cnveg_ns.livestemn_patch[p] + cnveg_ns.livecrootn_patch[p]

        if cnveg_state.n_allometry_patch[p] > 0.0
            cnveg_state.plantCN_patch[p] = cnveg_state.c_allometry_patch[p] /
                                           cnveg_state.n_allometry_patch[p]
        else
            cnveg_state.plantCN_patch[p] = 0.0
        end
    end

    # Set up permyc, kc_active, kn_active, litterfall steps
    for istp in 1:NSTP
        for p in bounds_p
            mask_soilp[p] || continue
            vt = ivt[p]

            if istp == ECM_STEP
                permyc[p, istp]    = perecm[vt]
                kc_active[p, istp] = ekc_active_pft[vt]
                kn_active[p, istp] = ekn_active_pft[vt]
            else
                permyc[p, istp]    = 1.0 - perecm[vt]
                kc_active[p, istp] = akc_active_pft[vt]
                kn_active[p, istp] = akn_active_pft[vt]
            end

            if cnveg_cs.leafc_patch[p] > 0.0
                litterfall_c_step[p, istp] = dt * permyc[p, istp] *
                                              cnveg_cf.leafc_to_litter_fun_patch[p]
                litterfall_n_step[p, istp] = dt * permyc[p, istp] *
                                              cnveg_ns.leafn_patch[p] *
                                              cnveg_cf.leafc_to_litter_fun_patch[p] /
                                              cnveg_cs.leafc_patch[p]
            end

            if (season_decid[vt] == 1.0 || stress_decid[vt] == 1.0)
                if cnveg_state.offset_flag_patch[p] != 1.0
                    litterfall_n_step[p, istp] = 0.0
                    litterfall_c_step[p, istp] = 0.0
                end
            end
        end
    end

    # Compute soil N layers
    for j in 1:nlevdecomp
        for p in bounds_p
            mask_soilp[p] || continue
            c = patch.column[p]

            sminn_no3_layer[c, j] = soilbgc_nf.smin_no3_to_plant_vr_col[c, j] *
                                    dzsoi_decomp_vals[j] * dt
            sminn_nh4_layer[c, j] = soilbgc_nf.smin_nh4_to_plant_vr_col[c, j] *
                                    dzsoi_decomp_vals[j] * dt

            if waterstate.h2osoi_liq_col[c, j] < smallValue
                sminn_no3_layer[c, j] = 0.0
                sminn_nh4_layer[c, j] = 0.0
            end

            sminn_no3_layer[c, j] = max(sminn_no3_layer[c, j], 0.0)
            sminn_nh4_layer[c, j] = max(sminn_nh4_layer[c, j], 0.0)

            if waterstate.h2osoi_liq_col[c, j] > smallValue
                sminn_no3_conc[c, j] = sminn_no3_layer[c, j] /
                                       (waterstate.h2osoi_liq_col[c, j] * 1000.0)
                sminn_nh4_conc[c, j] = sminn_nh4_layer[c, j] /
                                       (waterstate.h2osoi_liq_col[c, j] * 1000.0)
            else
                sminn_no3_conc[c, j] = 0.0
                sminn_nh4_conc[c, j] = 0.0
            end
        end
    end

    # Split by permyc
    for istp in 1:NSTP
        for j in 1:nlevdecomp
            for p in bounds_p
                mask_soilp[p] || continue
                c = patch.column[p]

                sminn_no3_layer_step[p, j, istp] = sminn_no3_layer[c, j] * permyc[p, istp]
                sminn_nh4_layer_step[p, j, istp] = sminn_nh4_layer[c, j] * permyc[p, istp]
                sminn_no3_conc_step[p, j, istp]  = sminn_no3_conc[c, j]  * permyc[p, istp]
                sminn_nh4_conc_step[p, j, istp]  = sminn_nh4_conc[c, j]  * permyc[p, istp]
            end
        end
    end

    # ======================================================================
    # Phase 2: Main PFT loop
    # ======================================================================
    for p in bounds_p
        mask_soilp[p] || continue
        c = patch.column[p]
        vt = ivt[p]
        excess_carbon_acc = 0.0
        burned_off_carbon_local = 0.0

        cnveg_nf.sminn_to_plant_fun_nh4_vr_patch[p, :] .= 0.0
        cnveg_nf.sminn_to_plant_fun_no3_vr_patch[p, :] .= 0.0

        # Deciduous storage demand
        if season_decid[vt] == 1.0 || stress_decid[vt] == 1.0
            if cnveg_state.onset_flag_patch[p] == 1.0
                cnveg_cs.leafc_storage_xfer_acc_patch[p] += cnveg_cf.leafc_storage_to_xfer_patch[p] * dt
                cnveg_ns.leafn_storage_xfer_acc_patch[p] += cnveg_nf.leafn_storage_to_xfer_patch[p] * dt
            end
            if cnveg_state.offset_flag_patch[p] == 1.0
                cnveg_cs.storage_cdemand_patch[p] = cnveg_cs.leafc_storage_patch[p] /
                                                     (ndays_off * steppday)
                cnveg_ns.storage_ndemand_patch[p] = cnveg_ns.leafn_storage_xfer_acc_patch[p] /
                                                     (ndays_off * steppday)
                cnveg_ns.storage_ndemand_patch[p] = max(cnveg_ns.storage_ndemand_patch[p], 0.0)
            else
                cnveg_cs.storage_cdemand_patch[p] = 0.0
                cnveg_ns.storage_ndemand_patch[p] = 0.0
            end
        else
            cnveg_cs.storage_cdemand_patch[p] = 0.0
            cnveg_ns.storage_ndemand_patch[p] = 0.0
        end

        # Available carbon pool
        availc_pool[p] = cnveg_cf.availc_patch[p] * dt

        if availc_pool[p] > 0.0
            for j in 1:nlevdecomp
                rootc_dens[p, j] = soilstate.crootfr_patch[p, j] * rootC[p]
            end
        end

        plant_ndemand_pool[p] = cnveg_nf.plant_ndemand_patch[p] * dt
        plant_ndemand_pool[p] = max(plant_ndemand_pool[p], 0.0)
        cnveg_nf.plant_ndemand_retrans_patch[p] = cnveg_ns.storage_ndemand_patch[p]

        plantCN_p = cnveg_state.plantCN_patch[p]

        # ==================================================================
        # Substep loop (ECM then AM)
        # ==================================================================
        for istp in ECM_STEP:AM_STEP
            # Zero accumulators for substep
            sminn_no3_diff    = 0.0
            sminn_nh4_diff    = 0.0
            active_no3_limit1 = 0.0
            active_nh4_limit1 = 0.0

            n_from_active_no3 .= 0.0
            n_from_active_nh4 .= 0.0
            n_from_nonmyc_no3 .= 0.0
            n_from_nonmyc_nh4 .= 0.0
            n_from_fixation .= 0.0

            n_active_no3_acc[p, istp] = 0.0
            n_active_nh4_acc[p, istp] = 0.0
            n_nonmyc_no3_acc[p, istp] = 0.0
            n_nonmyc_nh4_acc[p, istp] = 0.0
            n_fix_acc[p, istp] = 0.0
            n_retrans_acc[p, istp] = 0.0
            free_nretrans_acc[p, istp] = 0.0

            npp_active_no3_acc[p, istp] = 0.0
            npp_active_nh4_acc[p, istp] = 0.0
            npp_nonmyc_no3_acc[p, istp] = 0.0
            npp_nonmyc_nh4_acc[p, istp] = 0.0
            npp_fix_acc[p, istp] = 0.0
            npp_retrans_acc[p, istp] = 0.0

            npp_to_active_no3 .= 0.0
            npp_to_active_nh4 .= 0.0
            npp_to_nonmyc_no3 .= 0.0
            npp_to_nonmyc_nh4 .= 0.0
            npp_to_fixation .= 0.0
            npp_to_retrans .= 0.0

            plant_ndemand_pool_step[p, istp] = plant_ndemand_pool[p] * permyc[p, istp]
            npp_remaining[p, istp]           = availc_pool[p] * permyc[p, istp]

            # Compute fixation costs
            for j in 1:nlevdecomp
                tc_soisno_local[c, j] = temperature.t_soisno_col[c, j] - TFRZ
                fixer = round(Int, c3psn_pft[vt]) == 1 ? 1 : 0
                costNit[j, ICOST_FIX] = fun_cost_fix(fixer, a_fix_pft[vt],
                    b_fix_pft[vt], c_fix_pft[vt], BIG_COST,
                    soilstate.crootfr_patch[p, j], s_fix_pft[vt],
                    tc_soisno_local[c, j])
            end
            cost_fix_local[p, 1:nlevdecomp] .= @view costNit[:, ICOST_FIX]

            # Mycorrhizal uptake costs
            for j in 1:nlevdecomp
                rootc_dens_step = rootc_dens[p, j] * permyc[p, istp]
                costNit[j, ICOST_ACTIVE_NO3] = fun_cost_active(
                    sminn_no3_layer_step[p, j, istp], BIG_COST,
                    kc_active[p, istp], kn_active[p, istp],
                    rootc_dens_step, soilstate.crootfr_patch[p, j], smallValue)
                costNit[j, ICOST_ACTIVE_NH4] = fun_cost_active(
                    sminn_nh4_layer_step[p, j, istp], BIG_COST,
                    kc_active[p, istp], kn_active[p, istp],
                    rootc_dens_step, soilstate.crootfr_patch[p, j], smallValue)
            end
            cost_active_no3[p, 1:nlevdecomp] .= @view costNit[:, ICOST_ACTIVE_NO3]
            cost_active_nh4[p, 1:nlevdecomp] .= @view costNit[:, ICOST_ACTIVE_NH4]

            # Non-mycorrhizal uptake costs
            for j in 1:nlevdecomp
                rootc_dens_step = rootc_dens[p, j] * permyc[p, istp]
                costNit[j, ICOST_NONMYC_NO3] = fun_cost_nonmyc(
                    sminn_no3_layer_step[p, j, istp], BIG_COST,
                    kc_nonmyc_pft[vt], kn_nonmyc_pft[vt],
                    rootc_dens_step, soilstate.crootfr_patch[p, j], smallValue)
                costNit[j, ICOST_NONMYC_NH4] = fun_cost_nonmyc(
                    sminn_nh4_layer_step[p, j, istp], BIG_COST,
                    kc_nonmyc_pft[vt], kn_nonmyc_pft[vt],
                    rootc_dens_step, soilstate.crootfr_patch[p, j], smallValue)
            end
            cost_nonmyc_no3[p, 1:nlevdecomp] .= @view costNit[:, ICOST_NONMYC_NO3]
            cost_nonmyc_nh4[p, 1:nlevdecomp] .= @view costNit[:, ICOST_NONMYC_NH4]

            # Remove C required to pair with N from passive uptake
            npp_remaining[p, istp] -= n_passive_step[p, istp] * plantCN_p

            # ==============================================================
            # Fix loop: fixers then non-fixers
            # ==============================================================
            for FIX in PLANTS_ARE_FIXING:PLANTS_NOT_FIXING
                if FIX == PLANTS_ARE_FIXING
                    fixerfrac = FUN_fracfixers_pft[vt]
                else
                    fixerfrac = 1.0 - FUN_fracfixers_pft[vt]
                end
                npp_to_spend = npp_remaining[p, istp] * fixerfrac

                # Reset per-layer N uptake arrays
                n_from_active_no3 .= 0.0
                n_from_active_nh4 .= 0.0
                n_from_nonmyc_no3 .= 0.0
                n_from_nonmyc_nh4 .= 0.0

                # ---- Integrated conductance over soil column ----
                sum_n_acquired      = 0.0
                total_N_conductance = 0.0

                for j in 1:nlevdecomp
                    total_N_conductance += 1.0 / cost_active_no3[p, j] +
                                           1.0 / cost_active_nh4[p, j] +
                                           1.0 / cost_nonmyc_no3[p, j] +
                                           1.0 / cost_nonmyc_nh4[p, j]
                    if FIX == PLANTS_ARE_FIXING
                        total_N_conductance += 1.0 / cost_fix_local[p, j]
                    end
                end

                for j in 1:nlevdecomp
                    # NPP fraction allocation
                    npp_frac_to_active_nh4[j] = (1.0 / cost_active_nh4[p, j]) / total_N_conductance
                    npp_frac_to_nonmyc_nh4[j] = (1.0 / cost_nonmyc_nh4[p, j]) / total_N_conductance
                    npp_frac_to_active_no3[j] = (1.0 / cost_active_no3[p, j]) / total_N_conductance
                    npp_frac_to_nonmyc_no3[j] = (1.0 / cost_nonmyc_no3[p, j]) / total_N_conductance

                    if FIX == PLANTS_ARE_FIXING
                        npp_frac_to_fixation[j] = (1.0 / cost_fix_local[p, j]) / total_N_conductance
                    else
                        npp_frac_to_fixation[j] = 0.0
                    end

                    # Hypothetical N exchange
                    if FIX == PLANTS_ARE_FIXING
                        n_exch_fixation[j] = npp_frac_to_fixation[j] / cost_fix_local[p, j]
                    else
                        n_exch_fixation[j] = 0.0
                    end

                    n_exch_active_nh4[j] = npp_frac_to_active_nh4[j] / cost_active_nh4[p, j]
                    n_exch_nonmyc_nh4[j] = npp_frac_to_nonmyc_nh4[j] / cost_nonmyc_nh4[p, j]
                    n_exch_active_no3[j] = npp_frac_to_active_no3[j] / cost_active_no3[p, j]
                    n_exch_nonmyc_no3[j] = npp_frac_to_nonmyc_no3[j] / cost_nonmyc_no3[p, j]

                    sum_n_acquired += n_exch_active_nh4[j] + n_exch_nonmyc_nh4[j] +
                                      n_exch_active_no3[j] + n_exch_nonmyc_no3[j]
                    if FIX == PLANTS_ARE_FIXING
                        sum_n_acquired += n_exch_fixation[j]
                    end
                end

                total_N_resistance = 1.0 / sum_n_acquired

                # ---- Retranslocation ----
                if cnveg_cs.leafc_patch[p] > 0.0 &&
                   litterfall_n_step[p, istp] * fixerfrac > 0.0 &&
                   vt < npcropmin

                    rt_result = fun_retranslocation(
                        dt, npp_to_spend,
                        litterfall_c_step[p, istp] * fixerfrac,
                        litterfall_n_step[p, istp] * fixerfrac,
                        total_N_resistance,
                        leafcn[vt], grperc_pft[vt], plantCN_p)

                    total_c_spent_retrans     = rt_result.total_c_spent_retrans
                    total_c_accounted_retrans = rt_result.total_c_accounted_retrans
                    free_n_retrans_val        = rt_result.free_n_retrans
                    paid_for_n_retrans        = rt_result.paid_for_n_retrans
                else
                    total_c_accounted_retrans = 0.0
                    total_c_spent_retrans     = 0.0
                    paid_for_n_retrans        = 0.0
                    free_n_retrans_val        = 0.0
                end

                # Add retrans fluxes to budgets
                npp_to_spend -= total_c_spent_retrans + total_c_accounted_retrans
                npp_retrans_acc[p, istp] += total_c_spent_retrans
                n_retrans_acc[p, istp]   += paid_for_n_retrans
                free_nretrans_acc[p, istp] += free_n_retrans_val

                # ---- Spend C on extracting N ----
                if plant_ndemand_pool_step[p, istp] > 0.0
                    # FlexibleCN adjustment
                    if local_use_flexibleCN
                        if cnveg_ns.leafn_patch[p] == 0.0
                            delta_CN = fun_cn_flex_c_pft[vt]
                        else
                            delta_CN = (cnveg_cs.leafc_patch[p] + cnveg_cs.leafc_storage_patch[p]) /
                                       (cnveg_ns.leafn_patch[p] + cnveg_ns.leafn_storage_patch[p]) -
                                       leafcn[vt]
                        end
                        frac_ideal_C_use = max(0.0, 1.0 - (total_N_resistance - fun_cn_flex_a_pft[vt]) /
                                           fun_cn_flex_b_pft[vt])
                        if delta_CN < 0.0
                            frac_ideal_C_use += (1.0 - frac_ideal_C_use) *
                                                min(1.0, delta_CN / fun_cn_flex_c_pft[vt])
                        end
                        if delta_CN > 0.0 && frac_ideal_C_use < 1.0
                            frac_ideal_C_use += 0.5 * (1.0 * delta_CN / fun_cn_flex_c_pft[vt])
                        end
                        frac_ideal_C_use = max(min(1.0, frac_ideal_C_use), 0.5)
                    else
                        frac_ideal_C_use = 1.0
                    end

                    excess_carbon = npp_to_spend * (1.0 - frac_ideal_C_use)
                    if excess_carbon * (1.0 + grperc_pft[vt]) > npp_to_spend
                        excess_carbon = npp_to_spend / (1.0 + grperc_pft[vt])
                    end
                    excess_carbon_acc += excess_carbon

                    npp_to_spend -= excess_carbon * (1.0 + grperc_pft[vt])

                    # Main FUN equation
                    dnpp = npp_to_spend / ((1.0 + grperc_pft[vt]) * (plantCN_p / total_N_resistance) + 1.0)
                    dnpp *= frac_ideal_C_use

                    dn = dnpp / total_N_resistance

                    for j in 1:nlevdecomp
                        npp_to_active_nh4[j] = npp_frac_to_active_nh4[j] * dnpp
                        npp_to_nonmyc_nh4[j] = npp_frac_to_nonmyc_nh4[j] * dnpp
                        npp_to_active_no3[j] = npp_frac_to_active_no3[j] * dnpp
                        npp_to_nonmyc_no3[j] = npp_frac_to_nonmyc_no3[j] * dnpp

                        if FIX == PLANTS_ARE_FIXING
                            npp_to_fixation[j] = npp_frac_to_fixation[j] * dnpp
                        else
                            npp_to_fixation[j] = 0.0
                        end

                        n_from_active_nh4[j] = npp_to_active_nh4[j] / cost_active_nh4[p, j]
                        n_from_nonmyc_nh4[j] = npp_to_nonmyc_nh4[j] / cost_nonmyc_nh4[p, j]
                        n_from_active_no3[j] = npp_to_active_no3[j] / cost_active_no3[p, j]
                        n_from_nonmyc_no3[j] = npp_to_nonmyc_no3[j] / cost_nonmyc_no3[p, j]

                        if FIX == PLANTS_ARE_FIXING
                            n_from_fixation[j] = npp_to_fixation[j] / cost_fix_local[p, j]
                        else
                            n_from_fixation[j] = 0.0
                        end
                    end

                    # Check uptake limits
                    for j in 1:nlevdecomp
                        # NO3 limit
                        active_no3_limit1 = sminn_no3_layer_step[p, j, istp] * fixerfrac
                        if n_from_active_no3[j] + n_from_nonmyc_no3[j] > active_no3_limit1
                            sminn_no3_diff = n_from_active_no3[j] + n_from_nonmyc_no3[j] - active_no3_limit1
                            temp_n_flux = n_from_active_no3[j]
                            n_from_active_no3[j] -= sminn_no3_diff *
                                (n_from_active_no3[j] / (n_from_active_no3[j] + n_from_nonmyc_no3[j]))
                            n_from_nonmyc_no3[j] -= sminn_no3_diff *
                                (n_from_nonmyc_no3[j] / (temp_n_flux + n_from_nonmyc_no3[j]))
                            npp_to_active_no3[j] = n_from_active_no3[j] * cost_active_no3[p, j]
                            npp_to_nonmyc_no3[j] = n_from_nonmyc_no3[j] * cost_nonmyc_no3[p, j]
                        end

                        # NH4 limit
                        active_nh4_limit1 = sminn_nh4_layer_step[p, j, istp] * fixerfrac
                        if n_from_active_nh4[j] + n_from_nonmyc_nh4[j] > active_nh4_limit1
                            sminn_nh4_diff = n_from_active_nh4[j] + n_from_nonmyc_nh4[j] - active_nh4_limit1
                            temp_n_flux = n_from_active_nh4[j]
                            n_from_active_nh4[j] -= sminn_nh4_diff *
                                n_from_active_nh4[j] / (n_from_active_nh4[j] + n_from_nonmyc_nh4[j])
                            n_from_nonmyc_nh4[j] -= sminn_nh4_diff *
                                n_from_nonmyc_nh4[j] / (temp_n_flux + n_from_nonmyc_nh4[j])
                            npp_to_active_nh4[j] = n_from_active_nh4[j] * cost_active_nh4[p, j]
                            npp_to_nonmyc_nh4[j] = n_from_nonmyc_nh4[j] * cost_nonmyc_nh4[p, j]
                        end

                        # Total N acquired and C spent in this layer
                        N_acquired = n_from_active_no3[j] + n_from_nonmyc_no3[j] +
                                     n_from_active_nh4[j] + n_from_nonmyc_nh4[j]
                        C_spent = npp_to_active_no3[j] + npp_to_nonmyc_no3[j] +
                                  npp_to_active_nh4[j] + npp_to_nonmyc_nh4[j]

                        if FIX == PLANTS_ARE_FIXING
                            N_acquired += n_from_fixation[j]
                            C_spent    += npp_to_fixation[j]
                        end

                        npp_to_spend -= C_spent + N_acquired * plantCN_p * (1.0 + grperc_pft[vt])

                        nt_uptake[p, istp]  += N_acquired
                        npp_uptake[p, istp] += C_spent

                        # N flux accumulation
                        n_active_no3_acc[p, istp] += n_from_active_no3[j]
                        n_active_nh4_acc[p, istp] += n_from_active_nh4[j]
                        n_nonmyc_no3_acc[p, istp] += n_from_nonmyc_no3[j]
                        n_nonmyc_nh4_acc[p, istp] += n_from_nonmyc_nh4[j]

                        # C flux accumulation
                        npp_active_no3_acc[p, istp] += npp_to_active_no3[j]
                        npp_active_nh4_acc[p, istp] += npp_to_active_nh4[j]
                        npp_nonmyc_no3_acc[p, istp] += npp_to_nonmyc_no3[j]
                        npp_nonmyc_nh4_acc[p, istp] += npp_to_nonmyc_nh4[j]

                        if FIX == PLANTS_ARE_FIXING
                            n_fix_acc[p, istp]   += n_from_fixation[j]
                            npp_fix_acc[p, istp] += npp_to_fixation[j]
                        end
                    end # j

                    # Burn off excess carbon
                    if npp_to_spend >= 1.0e-13
                        burned_off_carbon_local += npp_to_spend
                    end

                    # Vertical fluxes
                    for j in 1:nlevdecomp
                        n_active_no3_vr[p, j] += n_from_active_no3[j]
                        n_active_nh4_vr[p, j] += n_from_active_nh4[j]
                        n_nonmyc_no3_vr[p, j] += n_from_nonmyc_no3[j]
                        n_nonmyc_nh4_vr[p, j] += n_from_nonmyc_nh4[j]
                    end
                end # unmet demand
            end # FIX loop

            # ECM/AM accumulation
            if istp == ECM_STEP
                n_ecm_no3_acc[p] = n_active_no3_acc[p, istp]
                n_ecm_nh4_acc[p] = n_active_nh4_acc[p, istp]
            else
                n_am_no3_acc[p] = n_active_no3_acc[p, istp]
                n_am_nh4_acc[p] = n_active_nh4_acc[p, istp]
            end

            # Accumulate totals over istp
            n_active_no3_acc_total[p]    += n_active_no3_acc[p, istp]
            n_active_nh4_acc_total[p]    += n_active_nh4_acc[p, istp]
            n_nonmyc_no3_acc_total[p]    += n_nonmyc_no3_acc[p, istp]
            n_nonmyc_nh4_acc_total[p]    += n_nonmyc_nh4_acc[p, istp]
            n_fix_acc_total[p]           += n_fix_acc[p, istp]
            n_retrans_acc_total[p]       += n_retrans_acc[p, istp]
            free_Nretrans_local[p]       += free_nretrans_acc[p, istp]

            npp_active_no3_acc_total[p]  += npp_active_no3_acc[p, istp]
            npp_active_nh4_acc_total[p]  += npp_active_nh4_acc[p, istp]
            npp_nonmyc_no3_acc_total[p]  += npp_nonmyc_no3_acc[p, istp]
            npp_nonmyc_nh4_acc_total[p]  += npp_nonmyc_nh4_acc[p, istp]
            npp_fix_acc_total[p]         += npp_fix_acc[p, istp]
            npp_retrans_acc_total[p]     += npp_retrans_acc[p, istp]
        end # istp

        # ==================================================================
        # Convert step-level quantities back to fluxes per second
        # ==================================================================

        # ---- N fluxes ----
        cnveg_nf.Npassive_patch[p]          = n_passive_acc[p] / dt
        cnveg_nf.Nfix_patch[p]              = n_fix_acc_total[p] / dt
        cnveg_nf.retransn_to_npool_patch[p] = n_retrans_acc_total[p] / dt

        if !use_matrixcn
            cnveg_nf.free_retransn_to_npool_patch[p] = free_Nretrans_local[p] / dt
        else
            if cnveg_ns.retransn_patch[p] > 0.0
                # Simplified matrix update path (skip matrix_update_phn for now)
                cnveg_nf.free_retransn_to_npool_patch[p] = free_Nretrans_local[p] / dt
            else
                cnveg_nf.free_retransn_to_npool_patch[p] = 0.0
            end
        end

        cnveg_nf.Nretrans_patch[p] = cnveg_nf.retransn_to_npool_patch[p] +
                                      cnveg_nf.free_retransn_to_npool_patch[p]

        # Extract active uptake N from soil pools
        for j in 1:nlevdecomp
            cnveg_nf.sminn_to_plant_fun_no3_vr_patch[p, j] =
                (n_passive_no3_vr[p, j] + n_active_no3_vr[p, j] +
                 n_nonmyc_no3_vr[p, j]) / (dzsoi_decomp_vals[j] * dt)
            cnveg_nf.sminn_to_plant_fun_nh4_vr_patch[p, j] =
                (n_passive_nh4_vr[p, j] + n_active_nh4_vr[p, j] +
                 n_nonmyc_nh4_vr[p, j]) / (dzsoi_decomp_vals[j] * dt)
        end

        cnveg_nf.Nactive_no3_patch[p] = n_active_no3_acc_total[p] / dt +
                                         n_active_no3_retrans_total[p] / dt
        cnveg_nf.Nactive_nh4_patch[p] = n_active_nh4_acc_total[p] / dt +
                                         n_active_nh4_retrans_total[p] / dt

        cnveg_nf.Necm_no3_patch[p] = n_ecm_no3_acc[p] / dt + n_ecm_no3_retrans[p] / dt
        cnveg_nf.Necm_nh4_patch[p] = n_ecm_nh4_acc[p] / dt + n_ecm_nh4_retrans[p] / dt
        cnveg_nf.Necm_patch[p]     = cnveg_nf.Necm_no3_patch[p] + cnveg_nf.Necm_nh4_patch[p]

        cnveg_nf.Nam_no3_patch[p] = n_am_no3_acc[p] / dt + n_am_no3_retrans[p] / dt
        cnveg_nf.Nam_nh4_patch[p] = n_am_nh4_acc[p] / dt + n_am_nh4_retrans[p] / dt
        cnveg_nf.Nam_patch[p]     = cnveg_nf.Nam_no3_patch[p] + cnveg_nf.Nam_nh4_patch[p]

        cnveg_nf.Nnonmyc_no3_patch[p] = n_nonmyc_no3_acc_total[p] / dt +
                                          n_nonmyc_no3_retrans_total[p] / dt
        cnveg_nf.Nnonmyc_nh4_patch[p] = n_nonmyc_nh4_acc_total[p] / dt +
                                          n_nonmyc_nh4_retrans_total[p] / dt
        cnveg_nf.Nnonmyc_patch[p] = cnveg_nf.Nnonmyc_no3_patch[p] +
                                     cnveg_nf.Nnonmyc_nh4_patch[p]

        cnveg_nf.plant_ndemand_retrans_patch[p] /= dt

        cnveg_nf.Nuptake_patch[p] = cnveg_nf.Nactive_no3_patch[p] +
            cnveg_nf.Nactive_nh4_patch[p] + cnveg_nf.Nnonmyc_no3_patch[p] +
            cnveg_nf.Nnonmyc_nh4_patch[p] + cnveg_nf.Nfix_patch[p] +
            cnveg_nf.Npassive_patch[p] + cnveg_nf.retransn_to_npool_patch[p] +
            cnveg_nf.free_retransn_to_npool_patch[p]

        cnveg_nf.Nactive_patch[p] = cnveg_nf.Nactive_no3_patch[p] +
            cnveg_nf.Nactive_nh4_patch[p] + cnveg_nf.Nnonmyc_no3_patch[p] +
            cnveg_nf.Nnonmyc_nh4_patch[p]

        cnveg_nf.sminn_to_plant_fun_patch[p] = cnveg_nf.Nactive_no3_patch[p] +
            cnveg_nf.Nactive_nh4_patch[p] + cnveg_nf.Nnonmyc_no3_patch[p] +
            cnveg_nf.Nnonmyc_nh4_patch[p] + cnveg_nf.Nfix_patch[p] +
            cnveg_nf.Npassive_patch[p]

        # ---- C fluxes ----
        cnveg_cf.npp_Nactive_no3_patch[p] = npp_active_no3_acc_total[p] / dt +
                                              npp_active_no3_retrans_total[p] / dt
        cnveg_cf.npp_Nactive_nh4_patch[p] = npp_active_nh4_acc_total[p] / dt +
                                              npp_active_nh4_retrans_total[p] / dt
        cnveg_cf.npp_Nnonmyc_no3_patch[p] = npp_nonmyc_no3_acc_total[p] / dt +
                                               npp_nonmyc_no3_retrans_total[p] / dt
        cnveg_cf.npp_Nnonmyc_nh4_patch[p] = npp_nonmyc_nh4_acc_total[p] / dt +
                                               npp_nonmyc_nh4_retrans_total[p] / dt
        cnveg_cf.npp_Nactive_patch[p] = cnveg_cf.npp_Nactive_no3_patch[p] +
            cnveg_cf.npp_Nactive_nh4_patch[p] + cnveg_cf.npp_Nnonmyc_no3_patch[p] +
            cnveg_cf.npp_Nnonmyc_nh4_patch[p]
        cnveg_cf.npp_Nnonmyc_patch[p] = cnveg_cf.npp_Nnonmyc_no3_patch[p] +
            cnveg_cf.npp_Nnonmyc_nh4_patch[p]
        cnveg_cf.npp_Nfix_patch[p]    = npp_fix_acc_total[p] / dt
        cnveg_cf.npp_Nretrans_patch[p] = npp_retrans_acc_total[p] / dt

        # Extra respiration fluxes
        cnveg_cf.soilc_change_patch[p] = (npp_active_no3_acc_total[p] +
            npp_active_nh4_acc_total[p] + npp_nonmyc_no3_acc_total[p] +
            npp_nonmyc_nh4_acc_total[p] + npp_fix_acc_total[p]) / dt +
            cnveg_cf.npp_Nretrans_patch[p]
        cnveg_cf.soilc_change_patch[p] += burned_off_carbon_local / dt
        cnveg_cf.npp_burnedoff_patch[p] = burned_off_carbon_local / dt
        cnveg_cf.npp_Nuptake_patch[p]   = cnveg_cf.soilc_change_patch[p]

        cnveg_cf.npp_growth_patch[p] = (cnveg_nf.Nuptake_patch[p] -
            cnveg_nf.free_retransn_to_npool_patch[p]) * plantCN_p +
            (excess_carbon_acc / dt)

        # ---- Diagnostic fluxes ----
        if cnveg_cf.availc_patch[p] > 0.0
            cnveg_nf.nuptake_npp_fraction_patch[p] = cnveg_cf.npp_Nuptake_patch[p] /
                cnveg_cf.availc_patch[p]
        else
            cnveg_nf.nuptake_npp_fraction_patch[p] = spval
        end

        if cnveg_cf.npp_Nfix_patch[p] > 0.0
            cnveg_nf.cost_nfix_patch[p] = cnveg_nf.Nfix_patch[p] / cnveg_cf.npp_Nfix_patch[p]
        else
            cnveg_nf.cost_nfix_patch[p] = spval
        end

        if cnveg_cf.npp_Nactive_patch[p] > 0.0
            cnveg_nf.cost_nactive_patch[p] = cnveg_nf.Nactive_patch[p] / cnveg_cf.npp_Nactive_patch[p]
        else
            cnveg_nf.cost_nactive_patch[p] = spval
        end

        if cnveg_cf.npp_Nretrans_patch[p] > 0.0
            cnveg_nf.cost_nretrans_patch[p] = cnveg_nf.Nretrans_patch[p] / cnveg_cf.npp_Nretrans_patch[p]
        else
            cnveg_nf.cost_nretrans_patch[p] = spval
        end
    end # PFT loop

    # ======================================================================
    # Phase 3: Patch-to-column aggregation (simplified p2c)
    # ======================================================================
    for c in bounds_c
        mask_soilc[c] || continue
        soilbgc_cf.soilc_change_col[c] = 0.0
        soilbgc_nf.nfix_to_sminn_col[c] = 0.0
    end

    for p in bounds_p
        mask_soilp[p] || continue
        c = patch.column[p]
        wtcol = patch.wtcol[p]
        soilbgc_cf.soilc_change_col[c]  += cnveg_cf.soilc_change_patch[p] * wtcol
        soilbgc_nf.nfix_to_sminn_col[c] += cnveg_nf.Nfix_patch[p] * wtcol
    end

    return nothing
end
