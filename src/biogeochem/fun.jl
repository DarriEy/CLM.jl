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

"""
    pftcon_fun_from(p) -> PftConFUN

Build the FUN per-PFT constant bundle from the main `pftcon` (PftconType). All
20 FUN parameters are already loaded into the main pftcon by `pftcon_read!`
(same NetCDF names), so this is a straight field map — no separate param read.
"""
function pftcon_fun_from(p)
    PftConFUN(
        leafcn = p.leafcn, season_decid = p.season_decid, stress_decid = p.stress_decid,
        a_fix = p.a_fix, b_fix = p.b_fix, c_fix = p.c_fix, s_fix = p.s_fix,
        akc_active = p.akc_active, akn_active = p.akn_active,
        ekc_active = p.ekc_active, ekn_active = p.ekn_active,
        kc_nonmyc = p.kc_nonmyc, kn_nonmyc = p.kn_nonmyc,
        perecm = p.perecm, grperc = p.grperc,
        fun_cn_flex_a = p.fun_cn_flex_a, fun_cn_flex_b = p.fun_cn_flex_b,
        fun_cn_flex_c = p.fun_cn_flex_c,
        FUN_fracfixers = p.FUN_fracfixers, c3psn = p.c3psn)
end

"""
    _fun_p2c!(col_out, patch_in, patch, mask_c, bounds_c, bounds_p, nlev)

Unity-weighted patch→column average of a per-(patch,level) array, used to push
FUN's cost-based plant N uptake (`sminn_to_plant_fun_*_vr_patch`) up to the
column (`*_vr_col`) for the competition re-sum + soil-N state update. Mirrors
Fortran `p2c(..., 'unity')` after the CNFUN call.
"""
function _fun_p2c!(col_out, patch_in, patch, mask_c::AbstractVector{Bool},
                   bounds_c::UnitRange{Int}, bounds_p::UnitRange{Int}, nlev::Int)
    @inbounds for c in bounds_c
        mask_c[c] || continue
        for j in 1:nlev
            col_out[c, j] = 0.0
        end
    end
    @inbounds for p in bounds_p
        c = patch.column[p]
        (first(bounds_c) <= c <= last(bounds_c)) || continue
        mask_c[c] || continue
        w = patch.wtcol[p]
        isfinite(w) || continue
        for j in 1:nlev
            col_out[c, j] += patch_in[p, j] * w
        end
    end
    return nothing
end

# =============================================================================
#  Cost functions (pure functions, GPU-friendly)
# =============================================================================

"""
    fun_cost_fix(fixer, a_fix, b_fix, c_fix, big_cost, crootfr, s_fix, tc_soisno)

Cost of fixing N by nodules (gC/gN). Returns `big_cost` for non-fixers or
layers with negligible root fraction.
"""
function fun_cost_fix(fixer::Int, a_fix::Real, b_fix::Real,
                      c_fix::Real, big_cost::Real,
                      crootfr::Real, s_fix::Real,
                      tc_soisno::Real)
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
function fun_cost_active(sminn_layer::Real, big_cost::Real,
                         kc_active::Real, kn_active::Real,
                         rootc_dens::Real, crootfr::Real,
                         smallValue::Real)
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
function fun_cost_nonmyc(sminn_layer::Real, big_cost::Real,
                         kc_nonmyc::Real, kn_nonmyc::Real,
                         rootc_dens::Real, crootfr::Real,
                         smallValue::Real)
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
function fun_retranslocation(dt::Real, npp_to_spend::Real,
                             total_falling_leaf_c::Real,
                             total_falling_leaf_n::Real,
                             total_n_resistance::Real,
                             target_leafcn::Real,
                             grperc::Real,
                             plantCN::Real)
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
                     dt::Real=1800.0,
                     nstep::Int=1,
                     dayspyr::Real=365.0,
                     npcropmin::Int=17)

    timestep_fun = SECSPDAY * FUN_PERIOD
    nstep_fun    = round(Int, SECSPDAY * dayspyr / dt)

    if nstep_fun > 0 && mod(nstep, nstep_fun) == 0
        for p in bounds
            mask_soilp[p] || continue
            ivt = patch.itype[p] + 1  # 0-based Fortran → 1-based Julia
            cnveg_state.leafcn_offset_patch[p]         = pftcon.leafcn[ivt]
            cnveg_cs.storage_cdemand_patch[p]           = 0.0
            cnveg_ns.storage_ndemand_patch[p]           = 0.0
            cnveg_ns.leafn_storage_xfer_acc_patch[p]    = 0.0
            cnveg_cs.leafc_storage_xfer_acc_patch[p]    = 0.0
        end
    end
end

# =============================================================================
#  FUN Phase-1 pre-computation kernels (one thread per patch/column; internal
#  substep/level loops run sequentially in-thread, every write to the thread's
#  own patch/column index → race-free, CPU byte-identical). Consumed by the
#  main per-patch cnfun! body.
# =============================================================================

# Loop A: per-patch litterfall_n / rootC / plantN / plantCN
@kernel function _fun_p1a_kernel!(@Const(mask_soilp), npp_burnedoff, litterfall_n,
        rootC, plantN, plantCN, @Const(leafc_to_litter_fun), @Const(leafcn_offset),
        @Const(frootc), @Const(leafn), @Const(frootn), @Const(livestemn),
        @Const(livecrootn), @Const(n_allometry), @Const(c_allometry), dt)
    p = @index(Global)
    @inbounds if mask_soilp[p]
        npp_burnedoff[p] = zero(eltype(npp_burnedoff))
        litterfall_n[p] = (leafc_to_litter_fun[p] / leafcn_offset[p]) * dt
        rootC[p] = frootc[p]
        plantN[p] = leafn[p] + frootn[p] + livestemn[p] + livecrootn[p]
        if n_allometry[p] > zero(eltype(n_allometry))
            plantCN[p] = c_allometry[p] / n_allometry[p]
        else
            plantCN[p] = zero(eltype(plantCN))
        end
    end
end

# Loop B: per-patch, internal substep loop — permyc / kc_active / kn_active /
# litterfall_{c,n}_step. `ivt` is the 0-based Fortran PFT index (add 1).
@kernel function _fun_p1b_kernel!(@Const(mask_soilp), permyc, kc_active, kn_active,
        litterfall_c_step, litterfall_n_step, @Const(ivt), @Const(perecm),
        @Const(ekc_active), @Const(ekn_active), @Const(akc_active), @Const(akn_active),
        @Const(season_decid), @Const(stress_decid), @Const(leafc), @Const(leafn),
        @Const(leafc_to_litter_fun), @Const(offset_flag), nstp::Int, ecm_step::Int, dt)
    p = @index(Global)
    @inbounds if mask_soilp[p]
        T = eltype(permyc)
        vt = ivt[p] + 1
        for istp in 1:nstp
            if istp == ecm_step
                permyc[p, istp]    = perecm[vt]
                kc_active[p, istp] = ekc_active[vt]
                kn_active[p, istp] = ekn_active[vt]
            else
                permyc[p, istp]    = one(T) - perecm[vt]
                kc_active[p, istp] = akc_active[vt]
                kn_active[p, istp] = akn_active[vt]
            end
            if leafc[p] > zero(T)
                litterfall_c_step[p, istp] = dt * permyc[p, istp] * leafc_to_litter_fun[p]
                litterfall_n_step[p, istp] = dt * permyc[p, istp] * leafn[p] *
                                              leafc_to_litter_fun[p] / leafc[p]
            end
            if (season_decid[vt] == one(T) || stress_decid[vt] == one(T))
                if offset_flag[p] != one(T)
                    litterfall_n_step[p, istp] = zero(T)
                    litterfall_c_step[p, istp] = zero(T)
                end
            end
        end
    end
end

# Loop C: per-soil-column, internal level loop — sminn NO3/NH4 layer & conc.
@kernel function _fun_p1c_kernel!(@Const(mask_soilc), sminn_no3_layer, sminn_nh4_layer,
        sminn_no3_conc, sminn_nh4_conc, @Const(smin_no3_to_plant_vr),
        @Const(smin_nh4_to_plant_vr), @Const(h2osoi_liq), @Const(dzsoi),
        nlevdecomp::Int, smallValue, dt)
    c = @index(Global)
    @inbounds if mask_soilc[c]
        T = eltype(sminn_no3_layer)
        sv = T(smallValue)
        for j in 1:nlevdecomp
            sminn_no3_layer[c, j] = smin_no3_to_plant_vr[c, j] * dzsoi[j] * dt
            sminn_nh4_layer[c, j] = smin_nh4_to_plant_vr[c, j] * dzsoi[j] * dt
            if h2osoi_liq[c, j] < sv
                sminn_no3_layer[c, j] = zero(T)
                sminn_nh4_layer[c, j] = zero(T)
            end
            sminn_no3_layer[c, j] = max(sminn_no3_layer[c, j], zero(T))
            sminn_nh4_layer[c, j] = max(sminn_nh4_layer[c, j], zero(T))
            if h2osoi_liq[c, j] > sv
                sminn_no3_conc[c, j] = sminn_no3_layer[c, j] / (h2osoi_liq[c, j] * T(1000.0))
                sminn_nh4_conc[c, j] = sminn_nh4_layer[c, j] / (h2osoi_liq[c, j] * T(1000.0))
            else
                sminn_no3_conc[c, j] = zero(T)
                sminn_nh4_conc[c, j] = zero(T)
            end
        end
    end
end

# Loop D: per-patch, internal substep+level — split the column sminn by permyc.
@kernel function _fun_p1d_kernel!(@Const(mask_soilp), sminn_no3_layer_step,
        sminn_nh4_layer_step, sminn_no3_conc_step, sminn_nh4_conc_step,
        @Const(sminn_no3_layer), @Const(sminn_nh4_layer), @Const(sminn_no3_conc),
        @Const(sminn_nh4_conc), @Const(permyc), @Const(column), nlevdecomp::Int, nstp::Int)
    p = @index(Global)
    @inbounds if mask_soilp[p]
        c = column[p]
        for istp in 1:nstp
            for j in 1:nlevdecomp
                sminn_no3_layer_step[p, j, istp] = sminn_no3_layer[c, j] * permyc[p, istp]
                sminn_nh4_layer_step[p, j, istp] = sminn_nh4_layer[c, j] * permyc[p, istp]
                sminn_no3_conc_step[p, j, istp]  = sminn_no3_conc[c, j]  * permyc[p, istp]
                sminn_nh4_conc_step[p, j, istp]  = sminn_nh4_conc[c, j]  * permyc[p, istp]
            end
        end
    end
end

# Phase-2 setup (per-patch): deciduous storage C/N demand, available-C pool,
# root-C density, plant N demand. Independent of the substep solver within a
# patch, so it runs as its own pass (byte-identical to the fused host loop).
@kernel function _fun_p2setup_kernel!(@Const(mask_soilp), storage_cdemand,
        leafc_storage_xfer_acc, storage_ndemand, leafn_storage_xfer_acc,
        availc_pool, rootc_dens, plant_ndemand_pool, plant_ndemand_retrans,
        @Const(ivt), @Const(season_decid), @Const(stress_decid), @Const(onset_flag),
        @Const(offset_flag), @Const(leafc_storage_to_xfer), @Const(leafn_storage_to_xfer),
        @Const(leafc_storage), @Const(availc), @Const(crootfr), @Const(rootC),
        @Const(plant_ndemand), ndays_off, steppday, nlevdecomp::Int, dt)
    p = @index(Global)
    @inbounds if mask_soilp[p]
        T = eltype(availc_pool)
        vt = ivt[p] + 1
        if season_decid[vt] == one(T) || stress_decid[vt] == one(T)
            if onset_flag[p] == one(T)
                leafc_storage_xfer_acc[p] += leafc_storage_to_xfer[p] * dt
                leafn_storage_xfer_acc[p] += leafn_storage_to_xfer[p] * dt
            end
            if offset_flag[p] == one(T)
                storage_cdemand[p] = leafc_storage[p] / (T(ndays_off) * T(steppday))
                storage_ndemand[p] = leafn_storage_xfer_acc[p] / (T(ndays_off) * T(steppday))
                storage_ndemand[p] = max(storage_ndemand[p], zero(T))
            else
                storage_cdemand[p] = zero(T)
                storage_ndemand[p] = zero(T)
            end
        else
            storage_cdemand[p] = zero(T)
            storage_ndemand[p] = zero(T)
        end

        availc_pool[p] = availc[p] * dt
        if availc_pool[p] > zero(T)
            for j in 1:nlevdecomp
                rootc_dens[p, j] = crootfr[p, j] * rootC[p]
            end
        end

        plant_ndemand_pool[p] = max(plant_ndemand[p] * dt, zero(T))
        plant_ndemand_retrans[p] = storage_ndemand[p]
    end
end

# Phase-2 substep solver (one thread per patch): the loop-carried ECM/AM ×
# fixer/non-fixer × soil-level machinery + retranslocation + the flux-conversion
# tail. Per-patch scratch (accumulators, totals, ecm/am, the always-zero retrans
# terms) are THREAD-LOCAL SCALARS (~46 arrays eliminated); only the genuinely
# per-(p,j)-persistent working arrays + inputs + outputs stay in the bundle `B`.
# This both simplifies the body and cuts the array count toward Metal's ~31-buffer
# kernel limit (remaining bundle arrays get dense-tensor packed next). Every write
# is to the thread's own patch → race-free, CPU byte-identical to the host loop.
@kernel function _fun_p2_kernel!(@Const(mask_soilp), @Const(column), @Const(ivt), B,
        dt, nlevdecomp::Int, smallValue, spval, npcropmin::Int,
        use_flexiblecn::Bool, use_matrixcn::Bool)
    p = @index(Global)
    @inbounds if mask_soilp[p]
        (; leafcn, grperc_pft, FUN_fracfixers_pft, c3psn_pft, a_fix_pft, b_fix_pft,
           c_fix_pft, s_fix_pft, kc_nonmyc_pft, kn_nonmyc_pft, fun_cn_flex_a_pft,
           fun_cn_flex_b_pft, fun_cn_flex_c_pft, season_decid, stress_decid,
           leafc, leafc_storage, leafn, leafn_storage, retransn, plantCN, crootfr,
           t_soisno, dzsoi, availc, permyc, kc_active, kn_active, availc_pool,
           plant_ndemand_pool, litterfall_c_step, litterfall_n_step, rootc_dens,
           sminn_no3_layer_step, sminn_nh4_layer_step, n_passive_step, n_passive_acc,
           n_passive_no3_vr, n_passive_nh4_vr, cost_fix_local, cost_active_no3,
           cost_active_nh4, cost_nonmyc_no3, cost_nonmyc_nh4, npp_to_fixation,
           npp_to_active_nh4, npp_to_nonmyc_nh4, npp_to_active_no3, npp_to_nonmyc_no3,
           npp_frac_to_fixation, npp_frac_to_active_nh4, npp_frac_to_nonmyc_nh4,
           npp_frac_to_active_no3, npp_frac_to_nonmyc_no3, n_from_fixation,
           n_from_active_nh4, n_from_nonmyc_nh4, n_from_active_no3, n_from_nonmyc_no3,
           n_active_no3_vr, n_active_nh4_vr, n_nonmyc_no3_vr, n_nonmyc_nh4_vr,
           Npassive, Nfix, retransn_to_npool, free_retransn_to_npool, Nretrans,
           sminn_fun_no3_vr, sminn_fun_nh4_vr, Nactive_no3, Nactive_nh4, Necm_no3,
           Necm_nh4, Necm, Nam_no3, Nam_nh4, Nam, Nnonmyc_no3, Nnonmyc_nh4, Nnonmyc,
           plant_ndemand_retrans, Nuptake, Nactive, sminn_to_plant_fun, nuptake_npp_fraction,
           cost_nfix, cost_nactive, cost_nretrans, npp_Nactive_no3, npp_Nactive_nh4,
           npp_Nnonmyc_no3, npp_Nnonmyc_nh4, npp_Nactive, npp_Nnonmyc, npp_Nfix,
           npp_Nretrans, soilc_change, npp_burnedoff, npp_Nuptake, npp_growth) = B

        c = column[p]
        vt = ivt[p] + 1
        local_use_flexibleCN = use_flexiblecn
        excess_carbon_acc = 0.0
        burned_off_carbon_local = 0.0

        for j in 1:nlevdecomp
            sminn_fun_nh4_vr[p, j] = 0.0
            sminn_fun_no3_vr[p, j] = 0.0
        end

        plantCN_p = plantCN[p]

        # per-patch totals (were [p] scratch arrays) — thread-local scalars
        n_active_no3_acc_total = 0.0; n_active_nh4_acc_total = 0.0
        n_nonmyc_no3_acc_total = 0.0; n_nonmyc_nh4_acc_total = 0.0
        n_fix_acc_total = 0.0; n_retrans_acc_total = 0.0; free_Nretrans_local = 0.0
        npp_active_no3_acc_total = 0.0; npp_active_nh4_acc_total = 0.0
        npp_nonmyc_no3_acc_total = 0.0; npp_nonmyc_nh4_acc_total = 0.0
        npp_fix_acc_total = 0.0; npp_retrans_acc_total = 0.0
        n_ecm_no3_acc = 0.0; n_ecm_nh4_acc = 0.0; n_am_no3_acc = 0.0; n_am_nh4_acc = 0.0

        for istp in ECM_STEP:AM_STEP
            sminn_no3_diff = 0.0; sminn_nh4_diff = 0.0
            active_no3_limit1 = 0.0; active_nh4_limit1 = 0.0

            for j in 1:nlevdecomp
                n_from_active_no3[p, j] = 0.0; n_from_active_nh4[p, j] = 0.0
                n_from_nonmyc_no3[p, j] = 0.0; n_from_nonmyc_nh4[p, j] = 0.0
                n_from_fixation[p, j] = 0.0
            end

            # per-(p,istp) accumulators (were [p,istp] scratch) — reset each substep
            n_active_no3_acc = 0.0; n_active_nh4_acc = 0.0
            n_nonmyc_no3_acc = 0.0; n_nonmyc_nh4_acc = 0.0
            n_fix_acc = 0.0; n_retrans_acc = 0.0; free_nretrans_acc = 0.0
            npp_active_no3_acc = 0.0; npp_active_nh4_acc = 0.0
            npp_nonmyc_no3_acc = 0.0; npp_nonmyc_nh4_acc = 0.0
            npp_fix_acc = 0.0; npp_retrans_acc = 0.0
            nt_uptake = 0.0; npp_uptake = 0.0

            for j in 1:nlevdecomp
                npp_to_active_no3[p, j] = 0.0; npp_to_active_nh4[p, j] = 0.0
                npp_to_nonmyc_no3[p, j] = 0.0; npp_to_nonmyc_nh4[p, j] = 0.0
                npp_to_fixation[p, j] = 0.0
            end

            plant_ndemand_pool_step = plant_ndemand_pool[p] * permyc[p, istp]
            npp_remaining = availc_pool[p] * permyc[p, istp]

            for j in 1:nlevdecomp
                tc = t_soisno[c, j] - TFRZ
                fixer = round(Int, c3psn_pft[vt]) == 1 ? 1 : 0
                cost_fix_local[p, j] = fun_cost_fix(fixer, a_fix_pft[vt], b_fix_pft[vt],
                    c_fix_pft[vt], BIG_COST, crootfr[p, j], s_fix_pft[vt], tc)
            end
            for j in 1:nlevdecomp
                rootc_dens_step = rootc_dens[p, j] * permyc[p, istp]
                cost_active_no3[p, j] = fun_cost_active(sminn_no3_layer_step[p, j, istp],
                    BIG_COST, kc_active[p, istp], kn_active[p, istp], rootc_dens_step,
                    crootfr[p, j], smallValue)
                cost_active_nh4[p, j] = fun_cost_active(sminn_nh4_layer_step[p, j, istp],
                    BIG_COST, kc_active[p, istp], kn_active[p, istp], rootc_dens_step,
                    crootfr[p, j], smallValue)
            end
            for j in 1:nlevdecomp
                rootc_dens_step = rootc_dens[p, j] * permyc[p, istp]
                cost_nonmyc_no3[p, j] = fun_cost_nonmyc(sminn_no3_layer_step[p, j, istp],
                    BIG_COST, kc_nonmyc_pft[vt], kn_nonmyc_pft[vt], rootc_dens_step,
                    crootfr[p, j], smallValue)
                cost_nonmyc_nh4[p, j] = fun_cost_nonmyc(sminn_nh4_layer_step[p, j, istp],
                    BIG_COST, kc_nonmyc_pft[vt], kn_nonmyc_pft[vt], rootc_dens_step,
                    crootfr[p, j], smallValue)
            end

            npp_remaining -= n_passive_step[p, istp] * plantCN_p

            for FIX in PLANTS_ARE_FIXING:PLANTS_NOT_FIXING
                if FIX == PLANTS_ARE_FIXING
                    fixerfrac = FUN_fracfixers_pft[vt]
                else
                    fixerfrac = 1.0 - FUN_fracfixers_pft[vt]
                end
                npp_to_spend = npp_remaining * fixerfrac

                for j in 1:nlevdecomp
                    n_from_active_no3[p, j] = 0.0; n_from_active_nh4[p, j] = 0.0
                    n_from_nonmyc_no3[p, j] = 0.0; n_from_nonmyc_nh4[p, j] = 0.0
                end

                sum_n_acquired = 0.0
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
                    npp_frac_to_active_nh4[p, j] = (1.0 / cost_active_nh4[p, j]) / total_N_conductance
                    npp_frac_to_nonmyc_nh4[p, j] = (1.0 / cost_nonmyc_nh4[p, j]) / total_N_conductance
                    npp_frac_to_active_no3[p, j] = (1.0 / cost_active_no3[p, j]) / total_N_conductance
                    npp_frac_to_nonmyc_no3[p, j] = (1.0 / cost_nonmyc_no3[p, j]) / total_N_conductance
                    if FIX == PLANTS_ARE_FIXING
                        npp_frac_to_fixation[p, j] = (1.0 / cost_fix_local[p, j]) / total_N_conductance
                    else
                        npp_frac_to_fixation[p, j] = 0.0
                    end
                    ne_fix = FIX == PLANTS_ARE_FIXING ? npp_frac_to_fixation[p, j] / cost_fix_local[p, j] : 0.0
                    ne_anh4 = npp_frac_to_active_nh4[p, j] / cost_active_nh4[p, j]
                    ne_nnh4 = npp_frac_to_nonmyc_nh4[p, j] / cost_nonmyc_nh4[p, j]
                    ne_ano3 = npp_frac_to_active_no3[p, j] / cost_active_no3[p, j]
                    ne_nno3 = npp_frac_to_nonmyc_no3[p, j] / cost_nonmyc_no3[p, j]
                    sum_n_acquired += ne_anh4 + ne_nnh4 + ne_ano3 + ne_nno3
                    if FIX == PLANTS_ARE_FIXING
                        sum_n_acquired += ne_fix
                    end
                end

                total_N_resistance = 1.0 / sum_n_acquired

                if leafc[p] > 0.0 && litterfall_n_step[p, istp] * fixerfrac > 0.0 && vt < npcropmin
                    rt_result = fun_retranslocation(dt, npp_to_spend,
                        litterfall_c_step[p, istp] * fixerfrac,
                        litterfall_n_step[p, istp] * fixerfrac, total_N_resistance,
                        leafcn[vt], grperc_pft[vt], plantCN_p)
                    total_c_spent_retrans     = rt_result.total_c_spent_retrans
                    total_c_accounted_retrans = rt_result.total_c_accounted_retrans
                    free_n_retrans_val        = rt_result.free_n_retrans
                    paid_for_n_retrans        = rt_result.paid_for_n_retrans
                else
                    total_c_accounted_retrans = 0.0; total_c_spent_retrans = 0.0
                    paid_for_n_retrans = 0.0; free_n_retrans_val = 0.0
                end

                npp_to_spend -= total_c_spent_retrans + total_c_accounted_retrans
                npp_retrans_acc += total_c_spent_retrans
                n_retrans_acc   += paid_for_n_retrans
                free_nretrans_acc += free_n_retrans_val

                if plant_ndemand_pool_step > 0.0
                    if local_use_flexibleCN
                        if leafn[p] == 0.0
                            delta_CN = fun_cn_flex_c_pft[vt]
                        else
                            delta_CN = (leafc[p] + leafc_storage[p]) /
                                       (leafn[p] + leafn_storage[p]) - leafcn[vt]
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

                    dnpp = npp_to_spend / ((1.0 + grperc_pft[vt]) * (plantCN_p / total_N_resistance) + 1.0)
                    dnpp *= frac_ideal_C_use
                    dn = dnpp / total_N_resistance

                    for j in 1:nlevdecomp
                        npp_to_active_nh4[p, j] = npp_frac_to_active_nh4[p, j] * dnpp
                        npp_to_nonmyc_nh4[p, j] = npp_frac_to_nonmyc_nh4[p, j] * dnpp
                        npp_to_active_no3[p, j] = npp_frac_to_active_no3[p, j] * dnpp
                        npp_to_nonmyc_no3[p, j] = npp_frac_to_nonmyc_no3[p, j] * dnpp
                        if FIX == PLANTS_ARE_FIXING
                            npp_to_fixation[p, j] = npp_frac_to_fixation[p, j] * dnpp
                        else
                            npp_to_fixation[p, j] = 0.0
                        end
                        n_from_active_nh4[p, j] = npp_to_active_nh4[p, j] / cost_active_nh4[p, j]
                        n_from_nonmyc_nh4[p, j] = npp_to_nonmyc_nh4[p, j] / cost_nonmyc_nh4[p, j]
                        n_from_active_no3[p, j] = npp_to_active_no3[p, j] / cost_active_no3[p, j]
                        n_from_nonmyc_no3[p, j] = npp_to_nonmyc_no3[p, j] / cost_nonmyc_no3[p, j]
                        if FIX == PLANTS_ARE_FIXING
                            n_from_fixation[p, j] = npp_to_fixation[p, j] / cost_fix_local[p, j]
                        else
                            n_from_fixation[p, j] = 0.0
                        end
                    end

                    for j in 1:nlevdecomp
                        active_no3_limit1 = sminn_no3_layer_step[p, j, istp] * fixerfrac
                        if n_from_active_no3[p, j] + n_from_nonmyc_no3[p, j] > active_no3_limit1
                            sminn_no3_diff = n_from_active_no3[p, j] + n_from_nonmyc_no3[p, j] - active_no3_limit1
                            temp_n_flux = n_from_active_no3[p, j]
                            n_from_active_no3[p, j] -= sminn_no3_diff *
                                (n_from_active_no3[p, j] / (n_from_active_no3[p, j] + n_from_nonmyc_no3[p, j]))
                            n_from_nonmyc_no3[p, j] -= sminn_no3_diff *
                                (n_from_nonmyc_no3[p, j] / (temp_n_flux + n_from_nonmyc_no3[p, j]))
                            npp_to_active_no3[p, j] = n_from_active_no3[p, j] * cost_active_no3[p, j]
                            npp_to_nonmyc_no3[p, j] = n_from_nonmyc_no3[p, j] * cost_nonmyc_no3[p, j]
                        end

                        active_nh4_limit1 = sminn_nh4_layer_step[p, j, istp] * fixerfrac
                        if n_from_active_nh4[p, j] + n_from_nonmyc_nh4[p, j] > active_nh4_limit1
                            sminn_nh4_diff = n_from_active_nh4[p, j] + n_from_nonmyc_nh4[p, j] - active_nh4_limit1
                            temp_n_flux = n_from_active_nh4[p, j]
                            n_from_active_nh4[p, j] -= sminn_nh4_diff *
                                n_from_active_nh4[p, j] / (n_from_active_nh4[p, j] + n_from_nonmyc_nh4[p, j])
                            n_from_nonmyc_nh4[p, j] -= sminn_nh4_diff *
                                n_from_nonmyc_nh4[p, j] / (temp_n_flux + n_from_nonmyc_nh4[p, j])
                            npp_to_active_nh4[p, j] = n_from_active_nh4[p, j] * cost_active_nh4[p, j]
                            npp_to_nonmyc_nh4[p, j] = n_from_nonmyc_nh4[p, j] * cost_nonmyc_nh4[p, j]
                        end

                        N_acquired = n_from_active_no3[p, j] + n_from_nonmyc_no3[p, j] +
                                     n_from_active_nh4[p, j] + n_from_nonmyc_nh4[p, j]
                        C_spent = npp_to_active_no3[p, j] + npp_to_nonmyc_no3[p, j] +
                                  npp_to_active_nh4[p, j] + npp_to_nonmyc_nh4[p, j]
                        if FIX == PLANTS_ARE_FIXING
                            N_acquired += n_from_fixation[p, j]
                            C_spent    += npp_to_fixation[p, j]
                        end

                        npp_to_spend -= C_spent + N_acquired * plantCN_p * (1.0 + grperc_pft[vt])
                        nt_uptake  += N_acquired
                        npp_uptake += C_spent

                        n_active_no3_acc += n_from_active_no3[p, j]
                        n_active_nh4_acc += n_from_active_nh4[p, j]
                        n_nonmyc_no3_acc += n_from_nonmyc_no3[p, j]
                        n_nonmyc_nh4_acc += n_from_nonmyc_nh4[p, j]
                        npp_active_no3_acc += npp_to_active_no3[p, j]
                        npp_active_nh4_acc += npp_to_active_nh4[p, j]
                        npp_nonmyc_no3_acc += npp_to_nonmyc_no3[p, j]
                        npp_nonmyc_nh4_acc += npp_to_nonmyc_nh4[p, j]
                        if FIX == PLANTS_ARE_FIXING
                            n_fix_acc   += n_from_fixation[p, j]
                            npp_fix_acc += npp_to_fixation[p, j]
                        end
                    end

                    if npp_to_spend >= 1.0e-13
                        burned_off_carbon_local += npp_to_spend
                    end

                    for j in 1:nlevdecomp
                        n_active_no3_vr[p, j] += n_from_active_no3[p, j]
                        n_active_nh4_vr[p, j] += n_from_active_nh4[p, j]
                        n_nonmyc_no3_vr[p, j] += n_from_nonmyc_no3[p, j]
                        n_nonmyc_nh4_vr[p, j] += n_from_nonmyc_nh4[p, j]
                    end
                end
            end # FIX

            if istp == ECM_STEP
                n_ecm_no3_acc = n_active_no3_acc
                n_ecm_nh4_acc = n_active_nh4_acc
            else
                n_am_no3_acc = n_active_no3_acc
                n_am_nh4_acc = n_active_nh4_acc
            end

            n_active_no3_acc_total    += n_active_no3_acc
            n_active_nh4_acc_total    += n_active_nh4_acc
            n_nonmyc_no3_acc_total    += n_nonmyc_no3_acc
            n_nonmyc_nh4_acc_total    += n_nonmyc_nh4_acc
            n_fix_acc_total           += n_fix_acc
            n_retrans_acc_total       += n_retrans_acc
            free_Nretrans_local       += free_nretrans_acc
            npp_active_no3_acc_total  += npp_active_no3_acc
            npp_active_nh4_acc_total  += npp_active_nh4_acc
            npp_nonmyc_no3_acc_total  += npp_nonmyc_no3_acc
            npp_nonmyc_nh4_acc_total  += npp_nonmyc_nh4_acc
            npp_fix_acc_total         += npp_fix_acc
            npp_retrans_acc_total     += npp_retrans_acc
        end # istp

        # ---- flux conversion (retrans-total terms are always 0 → dropped) ----
        Npassive[p]          = n_passive_acc[p] / dt
        Nfix[p]              = n_fix_acc_total / dt
        retransn_to_npool[p] = n_retrans_acc_total / dt
        if !use_matrixcn
            free_retransn_to_npool[p] = free_Nretrans_local / dt
        else
            if retransn[p] > 0.0
                free_retransn_to_npool[p] = free_Nretrans_local / dt
            else
                free_retransn_to_npool[p] = 0.0
            end
        end
        Nretrans[p] = retransn_to_npool[p] + free_retransn_to_npool[p]

        for j in 1:nlevdecomp
            sminn_fun_no3_vr[p, j] = (n_passive_no3_vr[p, j] + n_active_no3_vr[p, j] +
                 n_nonmyc_no3_vr[p, j]) / (dzsoi[j] * dt)
            sminn_fun_nh4_vr[p, j] = (n_passive_nh4_vr[p, j] + n_active_nh4_vr[p, j] +
                 n_nonmyc_nh4_vr[p, j]) / (dzsoi[j] * dt)
        end

        Nactive_no3[p] = n_active_no3_acc_total / dt
        Nactive_nh4[p] = n_active_nh4_acc_total / dt
        Necm_no3[p] = n_ecm_no3_acc / dt
        Necm_nh4[p] = n_ecm_nh4_acc / dt
        Necm[p]     = Necm_no3[p] + Necm_nh4[p]
        Nam_no3[p] = n_am_no3_acc / dt
        Nam_nh4[p] = n_am_nh4_acc / dt
        Nam[p]     = Nam_no3[p] + Nam_nh4[p]
        Nnonmyc_no3[p] = n_nonmyc_no3_acc_total / dt
        Nnonmyc_nh4[p] = n_nonmyc_nh4_acc_total / dt
        Nnonmyc[p] = Nnonmyc_no3[p] + Nnonmyc_nh4[p]
        plant_ndemand_retrans[p] /= dt
        Nuptake[p] = Nactive_no3[p] + Nactive_nh4[p] + Nnonmyc_no3[p] + Nnonmyc_nh4[p] +
            Nfix[p] + Npassive[p] + retransn_to_npool[p] + free_retransn_to_npool[p]
        Nactive[p] = Nactive_no3[p] + Nactive_nh4[p] + Nnonmyc_no3[p] + Nnonmyc_nh4[p]
        sminn_to_plant_fun[p] = Nactive_no3[p] + Nactive_nh4[p] + Nnonmyc_no3[p] +
            Nnonmyc_nh4[p] + Nfix[p] + Npassive[p]

        npp_Nactive_no3[p] = npp_active_no3_acc_total / dt
        npp_Nactive_nh4[p] = npp_active_nh4_acc_total / dt
        npp_Nnonmyc_no3[p] = npp_nonmyc_no3_acc_total / dt
        npp_Nnonmyc_nh4[p] = npp_nonmyc_nh4_acc_total / dt
        npp_Nactive[p] = npp_Nactive_no3[p] + npp_Nactive_nh4[p] + npp_Nnonmyc_no3[p] + npp_Nnonmyc_nh4[p]
        npp_Nnonmyc[p] = npp_Nnonmyc_no3[p] + npp_Nnonmyc_nh4[p]
        npp_Nfix[p]    = npp_fix_acc_total / dt
        npp_Nretrans[p] = npp_retrans_acc_total / dt

        soilc_change[p] = (npp_active_no3_acc_total + npp_active_nh4_acc_total +
            npp_nonmyc_no3_acc_total + npp_nonmyc_nh4_acc_total + npp_fix_acc_total) / dt +
            npp_Nretrans[p]
        soilc_change[p] += burned_off_carbon_local / dt
        npp_burnedoff[p] = burned_off_carbon_local / dt
        npp_Nuptake[p]   = soilc_change[p]
        npp_growth[p] = (Nuptake[p] - free_retransn_to_npool[p]) * plantCN_p + (excess_carbon_acc / dt)

        if availc[p] > 0.0
            nuptake_npp_fraction[p] = npp_Nuptake[p] / availc[p]
        else
            nuptake_npp_fraction[p] = spval
        end
        if npp_Nfix[p] > 0.0
            cost_nfix[p] = Nfix[p] / npp_Nfix[p]
        else
            cost_nfix[p] = spval
        end
        if npp_Nactive[p] > 0.0
            cost_nactive[p] = Nactive[p] / npp_Nactive[p]
        else
            cost_nactive[p] = spval
        end
        if npp_Nretrans[p] > 0.0
            cost_nretrans[p] = Nretrans[p] / npp_Nretrans[p]
        else
            cost_nretrans[p] = spval
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
                dt::Real=1800.0,
                nlevdecomp::Int=1,
                dzsoi_decomp_vals::Vector{<:Real}=Float64[],
                use_flexiblecn::Bool=false,
                use_matrixcn::Bool=false,
                npcropmin::Int=17,
                smallValue::Real=SMALLVALUE,
                spval::Real=SPVAL)

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

    # Infer floating-point type for AD compatibility
    FT = eltype(dzsoi_decomp_vals)

    # ---- Allocate local arrays ----
    rootc_dens       = zeros(FT, last(bounds_p), nlevdecomp)
    rootC            = zeros(FT, last(bounds_p))
    permyc           = zeros(FT, last(bounds_p), NSTP)
    kc_active        = zeros(FT, last(bounds_p), NSTP)
    kn_active        = zeros(FT, last(bounds_p), NSTP)
    availc_pool      = zeros(FT, last(bounds_p))
    plantN           = zeros(FT, last(bounds_p))
    plant_ndemand_pool       = zeros(FT, last(bounds_p))
    plant_ndemand_pool_step  = zeros(FT, last(bounds_p), NSTP)
    leafn_step               = zeros(FT, last(bounds_p), NSTP)
    leafn_retrans_step       = zeros(FT, last(bounds_p), NSTP)
    litterfall_n             = zeros(FT, last(bounds_p))
    litterfall_n_step        = zeros(FT, last(bounds_p), NSTP)
    litterfall_c_step        = zeros(FT, last(bounds_p), NSTP)
    npp_remaining            = zeros(FT, last(bounds_p), NSTP)
    n_passive_step           = zeros(FT, last(bounds_p), NSTP)
    n_passive_acc            = zeros(FT, last(bounds_p))

    cost_retran_local  = zeros(FT, last(bounds_p), nlevdecomp)
    cost_fix_local     = zeros(FT, last(bounds_p), nlevdecomp)
    cost_active_no3    = zeros(FT, last(bounds_p), nlevdecomp)
    cost_active_nh4    = zeros(FT, last(bounds_p), nlevdecomp)
    cost_nonmyc_no3    = zeros(FT, last(bounds_p), nlevdecomp)
    cost_nonmyc_nh4    = zeros(FT, last(bounds_p), nlevdecomp)

    n_fix_acc              = zeros(FT, last(bounds_p), NSTP)
    n_fix_acc_total        = zeros(FT, last(bounds_p))
    npp_fix_acc            = zeros(FT, last(bounds_p), NSTP)
    npp_fix_acc_total      = zeros(FT, last(bounds_p))
    n_retrans_acc          = zeros(FT, last(bounds_p), NSTP)
    n_retrans_acc_total    = zeros(FT, last(bounds_p))
    free_nretrans_acc      = zeros(FT, last(bounds_p), NSTP)
    npp_retrans_acc        = zeros(FT, last(bounds_p), NSTP)
    npp_retrans_acc_total  = zeros(FT, last(bounds_p))
    nt_uptake              = zeros(FT, last(bounds_p), NSTP)
    npp_uptake             = zeros(FT, last(bounds_p), NSTP)

    # NO3/NH4 arrays
    sminn_no3_conc       = zeros(FT, last(bounds_c), nlevdecomp)
    sminn_no3_conc_step  = zeros(FT, last(bounds_p), nlevdecomp, NSTP)
    sminn_no3_layer      = zeros(FT, last(bounds_c), nlevdecomp)
    sminn_no3_layer_step = zeros(FT, last(bounds_p), nlevdecomp, NSTP)
    sminn_no3_uptake     = zeros(FT, last(bounds_p), nlevdecomp, NSTP)
    sminn_nh4_conc       = zeros(FT, last(bounds_c), nlevdecomp)
    sminn_nh4_conc_step  = zeros(FT, last(bounds_p), nlevdecomp, NSTP)
    sminn_nh4_layer      = zeros(FT, last(bounds_c), nlevdecomp)
    sminn_nh4_layer_step = zeros(FT, last(bounds_p), nlevdecomp, NSTP)
    sminn_nh4_uptake     = zeros(FT, last(bounds_p), nlevdecomp, NSTP)

    n_active_no3_acc       = zeros(FT, last(bounds_p), NSTP)
    n_active_nh4_acc       = zeros(FT, last(bounds_p), NSTP)
    n_nonmyc_no3_acc       = zeros(FT, last(bounds_p), NSTP)
    n_nonmyc_nh4_acc       = zeros(FT, last(bounds_p), NSTP)
    n_active_no3_acc_total = zeros(FT, last(bounds_p))
    n_active_nh4_acc_total = zeros(FT, last(bounds_p))
    n_nonmyc_no3_acc_total = zeros(FT, last(bounds_p))
    n_nonmyc_nh4_acc_total = zeros(FT, last(bounds_p))
    npp_active_no3_acc       = zeros(FT, last(bounds_p), NSTP)
    npp_active_nh4_acc       = zeros(FT, last(bounds_p), NSTP)
    npp_nonmyc_no3_acc       = zeros(FT, last(bounds_p), NSTP)
    npp_nonmyc_nh4_acc       = zeros(FT, last(bounds_p), NSTP)
    npp_active_no3_acc_total = zeros(FT, last(bounds_p))
    npp_active_nh4_acc_total = zeros(FT, last(bounds_p))
    npp_nonmyc_no3_acc_total = zeros(FT, last(bounds_p))
    npp_nonmyc_nh4_acc_total = zeros(FT, last(bounds_p))

    n_am_no3_acc  = zeros(FT, last(bounds_p))
    n_am_nh4_acc  = zeros(FT, last(bounds_p))
    n_ecm_no3_acc = zeros(FT, last(bounds_p))
    n_ecm_nh4_acc = zeros(FT, last(bounds_p))
    n_am_no3_retrans  = zeros(FT, last(bounds_p))
    n_am_nh4_retrans  = zeros(FT, last(bounds_p))
    n_ecm_no3_retrans = zeros(FT, last(bounds_p))
    n_ecm_nh4_retrans = zeros(FT, last(bounds_p))

    n_active_no3_retrans       = zeros(FT, last(bounds_p), NSTP)
    n_active_nh4_retrans       = zeros(FT, last(bounds_p), NSTP)
    n_nonmyc_no3_retrans       = zeros(FT, last(bounds_p), NSTP)
    n_nonmyc_nh4_retrans       = zeros(FT, last(bounds_p), NSTP)
    n_active_no3_retrans_total = zeros(FT, last(bounds_p))
    n_active_nh4_retrans_total = zeros(FT, last(bounds_p))
    n_nonmyc_no3_retrans_total = zeros(FT, last(bounds_p))
    n_nonmyc_nh4_retrans_total = zeros(FT, last(bounds_p))
    npp_active_no3_retrans_total = zeros(FT, last(bounds_p))
    npp_active_nh4_retrans_total = zeros(FT, last(bounds_p))
    npp_nonmyc_no3_retrans_total = zeros(FT, last(bounds_p))
    npp_nonmyc_nh4_retrans_total = zeros(FT, last(bounds_p))

    n_passive_no3_vr = zeros(FT, last(bounds_p), nlevdecomp)
    n_passive_nh4_vr = zeros(FT, last(bounds_p), nlevdecomp)
    n_active_no3_vr  = zeros(FT, last(bounds_p), nlevdecomp)
    n_nonmyc_no3_vr  = zeros(FT, last(bounds_p), nlevdecomp)
    n_active_nh4_vr  = zeros(FT, last(bounds_p), nlevdecomp)
    n_nonmyc_nh4_vr  = zeros(FT, last(bounds_p), nlevdecomp)

    free_Nretrans_local = zeros(FT, last(bounds_p))


    # Per-layer working arrays — promoted to per-patch rows [np, nlevdecomp] so
    # each kernel thread owns its own layer scratch (npp_{to,frac}_to_retrans were
    # dead — reset but never read — and are dropped).
    npp_to_fixation     = zeros(FT, last(bounds_p), nlevdecomp)
    npp_to_active_nh4   = zeros(FT, last(bounds_p), nlevdecomp)
    npp_to_nonmyc_nh4   = zeros(FT, last(bounds_p), nlevdecomp)
    npp_to_active_no3   = zeros(FT, last(bounds_p), nlevdecomp)
    npp_to_nonmyc_no3   = zeros(FT, last(bounds_p), nlevdecomp)

    npp_frac_to_fixation     = zeros(FT, last(bounds_p), nlevdecomp)
    npp_frac_to_active_nh4   = zeros(FT, last(bounds_p), nlevdecomp)
    npp_frac_to_nonmyc_nh4   = zeros(FT, last(bounds_p), nlevdecomp)
    npp_frac_to_active_no3   = zeros(FT, last(bounds_p), nlevdecomp)
    npp_frac_to_nonmyc_no3   = zeros(FT, last(bounds_p), nlevdecomp)

    n_exch_fixation     = zeros(FT, last(bounds_p), nlevdecomp)
    n_exch_active_nh4   = zeros(FT, last(bounds_p), nlevdecomp)
    n_exch_nonmyc_nh4   = zeros(FT, last(bounds_p), nlevdecomp)
    n_exch_active_no3   = zeros(FT, last(bounds_p), nlevdecomp)
    n_exch_nonmyc_no3   = zeros(FT, last(bounds_p), nlevdecomp)

    n_from_fixation     = zeros(FT, last(bounds_p), nlevdecomp)
    n_from_active_nh4   = zeros(FT, last(bounds_p), nlevdecomp)
    n_from_nonmyc_nh4   = zeros(FT, last(bounds_p), nlevdecomp)
    n_from_active_no3   = zeros(FT, last(bounds_p), nlevdecomp)
    n_from_nonmyc_no3   = zeros(FT, last(bounds_p), nlevdecomp)

    # ======================================================================
    # Phase 1: Pre-computation (across all patches)
    # ======================================================================

    # Compute litterfall_n, rootC, plantN, plantCN (Loop A → per-patch kernel)
    _launch!(_fun_p1a_kernel!, mask_soilp, cnveg_cf.npp_burnedoff_patch, litterfall_n,
        rootC, plantN, cnveg_state.plantCN_patch, cnveg_cf.leafc_to_litter_fun_patch,
        cnveg_state.leafcn_offset_patch, cnveg_cs.frootc_patch, cnveg_ns.leafn_patch,
        cnveg_ns.frootn_patch, cnveg_ns.livestemn_patch, cnveg_ns.livecrootn_patch,
        cnveg_state.n_allometry_patch, cnveg_state.c_allometry_patch, dt)

    # Set up permyc, kc_active, kn_active, litterfall steps (Loop B → per-patch kernel)
    _launch!(_fun_p1b_kernel!, mask_soilp, permyc, kc_active, kn_active,
        litterfall_c_step, litterfall_n_step, ivt, perecm, ekc_active_pft, ekn_active_pft,
        akc_active_pft, akn_active_pft, season_decid, stress_decid, cnveg_cs.leafc_patch,
        cnveg_ns.leafn_patch, cnveg_cf.leafc_to_litter_fun_patch,
        cnveg_state.offset_flag_patch, NSTP, ECM_STEP, dt)

    # Compute soil N layers (Loop C → per-soil-column kernel; the original wrote
    # column-indexed quantities redundantly per patch — done once per column here).
    _launch!(_fun_p1c_kernel!, mask_soilc, sminn_no3_layer, sminn_nh4_layer,
        sminn_no3_conc, sminn_nh4_conc, soilbgc_nf.smin_no3_to_plant_vr_col,
        soilbgc_nf.smin_nh4_to_plant_vr_col, waterstate.h2osoi_liq_col,
        dzsoi_decomp_vals, nlevdecomp, smallValue, dt)

    # Split by permyc (Loop D → per-patch kernel, internal substep+level)
    _launch!(_fun_p1d_kernel!, mask_soilp, sminn_no3_layer_step, sminn_nh4_layer_step,
        sminn_no3_conc_step, sminn_nh4_conc_step, sminn_no3_layer, sminn_nh4_layer,
        sminn_no3_conc, sminn_nh4_conc, permyc, patch.column, nlevdecomp, NSTP)

    # ======================================================================
    # Phase 2: Main PFT loop — setup (storage demand / availc_pool / rootc_dens /
    # plant N demand) is a per-patch kernel; the substep solver below is
    # loop-carried and stays a per-patch host loop (kernelized in a follow-up).
    # ======================================================================
    _launch!(_fun_p2setup_kernel!, mask_soilp, cnveg_cs.storage_cdemand_patch,
        cnveg_cs.leafc_storage_xfer_acc_patch, cnveg_ns.storage_ndemand_patch,
        cnveg_ns.leafn_storage_xfer_acc_patch, availc_pool, rootc_dens, plant_ndemand_pool,
        cnveg_nf.plant_ndemand_retrans_patch, ivt, season_decid, stress_decid,
        cnveg_state.onset_flag_patch, cnveg_state.offset_flag_patch,
        cnveg_cf.leafc_storage_to_xfer_patch, cnveg_nf.leafn_storage_to_xfer_patch,
        cnveg_cs.leafc_storage_patch, cnveg_cf.availc_patch, soilstate.crootfr_patch,
        rootC, cnveg_nf.plant_ndemand_patch, ndays_off, steppday, nlevdecomp, dt)

    B_fun = (; leafcn, grperc_pft, FUN_fracfixers_pft, c3psn_pft, a_fix_pft, b_fix_pft,
        c_fix_pft, s_fix_pft, kc_nonmyc_pft, kn_nonmyc_pft, fun_cn_flex_a_pft,
        fun_cn_flex_b_pft, fun_cn_flex_c_pft, season_decid, stress_decid,
        leafc = cnveg_cs.leafc_patch, leafc_storage = cnveg_cs.leafc_storage_patch,
        leafn = cnveg_ns.leafn_patch, leafn_storage = cnveg_ns.leafn_storage_patch,
        retransn = cnveg_ns.retransn_patch, plantCN = cnveg_state.plantCN_patch,
        crootfr = soilstate.crootfr_patch, t_soisno = temperature.t_soisno_col,
        dzsoi = dzsoi_decomp_vals, availc = cnveg_cf.availc_patch,
        permyc, kc_active, kn_active, availc_pool, plant_ndemand_pool, litterfall_c_step,
        litterfall_n_step, rootc_dens, sminn_no3_layer_step, sminn_nh4_layer_step,
        n_passive_step, n_passive_acc, n_passive_no3_vr, n_passive_nh4_vr, cost_fix_local,
        cost_active_no3, cost_active_nh4, cost_nonmyc_no3, cost_nonmyc_nh4, npp_to_fixation,
        npp_to_active_nh4, npp_to_nonmyc_nh4, npp_to_active_no3, npp_to_nonmyc_no3,
        npp_frac_to_fixation, npp_frac_to_active_nh4, npp_frac_to_nonmyc_nh4,
        npp_frac_to_active_no3, npp_frac_to_nonmyc_no3, n_exch_fixation, n_exch_active_nh4,
        n_exch_nonmyc_nh4, n_exch_active_no3, n_exch_nonmyc_no3, n_from_fixation,
        n_from_active_nh4, n_from_nonmyc_nh4, n_from_active_no3, n_from_nonmyc_no3,
        n_active_no3_acc, n_active_nh4_acc, n_nonmyc_no3_acc, n_nonmyc_nh4_acc, n_fix_acc,
        n_retrans_acc, free_nretrans_acc, npp_active_no3_acc, npp_active_nh4_acc,
        npp_nonmyc_no3_acc, npp_nonmyc_nh4_acc, npp_fix_acc, npp_retrans_acc, nt_uptake,
        npp_uptake, plant_ndemand_pool_step, npp_remaining, n_active_no3_acc_total,
        n_active_nh4_acc_total, n_nonmyc_no3_acc_total, n_nonmyc_nh4_acc_total,
        n_fix_acc_total, n_retrans_acc_total, free_Nretrans_local, npp_active_no3_acc_total,
        npp_active_nh4_acc_total, npp_nonmyc_no3_acc_total, npp_nonmyc_nh4_acc_total,
        npp_fix_acc_total, npp_retrans_acc_total, n_ecm_no3_acc, n_ecm_nh4_acc, n_am_no3_acc,
        n_am_nh4_acc, n_active_no3_vr, n_active_nh4_vr, n_nonmyc_no3_vr, n_nonmyc_nh4_vr,
        n_active_no3_retrans_total, n_active_nh4_retrans_total, n_nonmyc_no3_retrans_total,
        n_nonmyc_nh4_retrans_total, npp_active_no3_retrans_total, npp_active_nh4_retrans_total,
        npp_nonmyc_no3_retrans_total, npp_nonmyc_nh4_retrans_total, n_am_no3_retrans,
        n_am_nh4_retrans, n_ecm_no3_retrans, n_ecm_nh4_retrans,
        Npassive = cnveg_nf.Npassive_patch, Nfix = cnveg_nf.Nfix_patch,
        retransn_to_npool = cnveg_nf.retransn_to_npool_patch,
        free_retransn_to_npool = cnveg_nf.free_retransn_to_npool_patch,
        Nretrans = cnveg_nf.Nretrans_patch,
        sminn_fun_no3_vr = cnveg_nf.sminn_to_plant_fun_no3_vr_patch,
        sminn_fun_nh4_vr = cnveg_nf.sminn_to_plant_fun_nh4_vr_patch,
        Nactive_no3 = cnveg_nf.Nactive_no3_patch, Nactive_nh4 = cnveg_nf.Nactive_nh4_patch,
        Necm_no3 = cnveg_nf.Necm_no3_patch, Necm_nh4 = cnveg_nf.Necm_nh4_patch,
        Necm = cnveg_nf.Necm_patch, Nam_no3 = cnveg_nf.Nam_no3_patch,
        Nam_nh4 = cnveg_nf.Nam_nh4_patch, Nam = cnveg_nf.Nam_patch,
        Nnonmyc_no3 = cnveg_nf.Nnonmyc_no3_patch, Nnonmyc_nh4 = cnveg_nf.Nnonmyc_nh4_patch,
        Nnonmyc = cnveg_nf.Nnonmyc_patch,
        plant_ndemand_retrans = cnveg_nf.plant_ndemand_retrans_patch,
        Nuptake = cnveg_nf.Nuptake_patch, Nactive = cnveg_nf.Nactive_patch,
        sminn_to_plant_fun = cnveg_nf.sminn_to_plant_fun_patch,
        nuptake_npp_fraction = cnveg_nf.nuptake_npp_fraction_patch,
        cost_nfix = cnveg_nf.cost_nfix_patch, cost_nactive = cnveg_nf.cost_nactive_patch,
        cost_nretrans = cnveg_nf.cost_nretrans_patch,
        npp_Nactive_no3 = cnveg_cf.npp_Nactive_no3_patch,
        npp_Nactive_nh4 = cnveg_cf.npp_Nactive_nh4_patch,
        npp_Nnonmyc_no3 = cnveg_cf.npp_Nnonmyc_no3_patch,
        npp_Nnonmyc_nh4 = cnveg_cf.npp_Nnonmyc_nh4_patch,
        npp_Nactive = cnveg_cf.npp_Nactive_patch, npp_Nnonmyc = cnveg_cf.npp_Nnonmyc_patch,
        npp_Nfix = cnveg_cf.npp_Nfix_patch, npp_Nretrans = cnveg_cf.npp_Nretrans_patch,
        soilc_change = cnveg_cf.soilc_change_patch, npp_burnedoff = cnveg_cf.npp_burnedoff_patch,
        npp_Nuptake = cnveg_cf.npp_Nuptake_patch, npp_growth = cnveg_cf.npp_growth_patch)

    _launch!(_fun_p2_kernel!, mask_soilp, patch.column, ivt, B_fun, dt, nlevdecomp,
        smallValue, spval, npcropmin, use_flexiblecn, use_matrixcn)


    # ======================================================================
    # Phase 3: Patch-to-column aggregation (p2c) — reuse the generic zero-init +
    # atomic scatter helpers (n_dynamics.jl). Zero the column accumulators, then
    # each active patch scatters patch*wtcol into its column (atomic on GPU, CPU
    # ascending-p == the old host += loop).
    # ======================================================================
    ndyn_col_zero!(soilbgc_cf.soilc_change_col, mask_soilc)
    ndyn_col_zero!(soilbgc_nf.nfix_to_sminn_col, mask_soilc)
    ndyn_p2c_scatter!(soilbgc_cf.soilc_change_col, mask_soilp, patch.column,
                      cnveg_cf.soilc_change_patch, patch.wtcol,
                      first(bounds_p), last(bounds_p))
    ndyn_p2c_scatter!(soilbgc_nf.nfix_to_sminn_col, mask_soilp, patch.column,
                      cnveg_nf.Nfix_patch, patch.wtcol,
                      first(bounds_p), last(bounds_p))

    return nothing
end
