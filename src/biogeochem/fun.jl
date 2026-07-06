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
@inline function fun_cost_fix(fixer::Int, a_fix, b_fix, c_fix, big_cost,
                              crootfr, s_fix, tc_soisno)
    T = typeof(tc_soisno)
    if fixer == 1 && crootfr > T(1.0e-6)
        return (-s_fix) /
               (T(1.25) * exp(a_fix + b_fix * tc_soisno * (one(T) - T(0.5) * tc_soisno / c_fix)))
    else
        return big_cost
    end
end

"""
    fun_cost_active(sminn_layer, big_cost, kc_active, kn_active, rootc_dens, crootfr, smallValue)

Cost of mycorrhizal active uptake of N (gC/gN).
"""
@inline function fun_cost_active(sminn_layer, big_cost, kc_active, kn_active,
                                 rootc_dens, crootfr, smallValue)
    if rootc_dens > oftype(rootc_dens, 1.0e-6) && sminn_layer > smallValue
        return kn_active / sminn_layer + kc_active / rootc_dens
    else
        return big_cost
    end
end

"""
    fun_cost_nonmyc(sminn_layer, big_cost, kc_nonmyc, kn_nonmyc, rootc_dens, crootfr, smallValue)

Cost of non-mycorrhizal uptake of N (gC/gN).
"""
@inline function fun_cost_nonmyc(sminn_layer, big_cost, kc_nonmyc, kn_nonmyc,
                                 rootc_dens, crootfr, smallValue)
    if rootc_dens > oftype(rootc_dens, 1.0e-6) && sminn_layer > smallValue
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
@inline function fun_retranslocation(dt, npp_to_spend,
                             total_falling_leaf_c,
                             total_falling_leaf_n,
                             total_n_resistance,
                             target_leafcn,
                             grperc,
                             plantCN)
    T = typeof(npp_to_spend)
    # Initialize total fluxes
    total_c_spent_retrans     = zero(T)
    total_c_accounted_retrans = zero(T)
    c_accounted_retrans       = zero(T)
    paid_for_n_retrans        = zero(T)
    npp_to_spend_temp         = npp_to_spend

    # Initial C and N pools in falling leaves
    falling_leaf_c = total_falling_leaf_c
    falling_leaf_n = total_falling_leaf_n

    # Parameters
    max_falling_leaf_cn = target_leafcn * T(3.0)
    min_falling_leaf_cn = target_leafcn * T(1.5)

    # Free uptake
    free_n_retrans = max(falling_leaf_n - (falling_leaf_c / min_falling_leaf_cn), zero(T))
    falling_leaf_n = falling_leaf_n - free_n_retrans

    # Initial CN ratio and costs
    falling_leaf_cn   = falling_leaf_c / falling_leaf_n
    kresorb           = one(T) / target_leafcn
    cost_retrans_temp = kresorb / ((one(T) / falling_leaf_cn)^T(1.3))

    # Iteration loop
    iter = 0
    exitloop = false
    while !exitloop && cost_retrans_temp < total_n_resistance &&
          falling_leaf_n >= zero(T) && npp_to_spend > zero(T)

        # Spend some C on removing N
        c_spent_retrans = cost_retrans_temp * (falling_leaf_n - falling_leaf_c /
                          (falling_leaf_cn + one(T)))
        # don't spend more C than you have
        c_spent_retrans = min(npp_to_spend_temp, c_spent_retrans)
        # N extracted
        leaf_n_ext = c_spent_retrans / cost_retrans_temp
        # Do not empty N pool
        leaf_n_ext = min(falling_leaf_n, leaf_n_ext)
        # C accounted for by N uptake
        c_accounted_retrans = leaf_n_ext * plantCN * (one(T) + grperc)

        # Update leafCN, recalculate costs
        falling_leaf_n = falling_leaf_n - leaf_n_ext
        if falling_leaf_n > zero(T)
            falling_leaf_cn   = falling_leaf_c / falling_leaf_n
            cost_retrans_temp = kresorb / ((one(T) / falling_leaf_cn)^T(1.3))
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
        if npp_to_spend_temp <= zero(T)
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

# ---- Dense-tensor layout for the Metal-buffer-fit substep kernel -------------
# Metal caps compute-kernel args at ~31 buffers; the substep touches ~103 arrays.
# Pack them into a handful of dense tensors and alias each logical array to a
# @view slice at the top of the kernel (verified: @view slices work inside a
# Metal kernel), so the solver body stays byte-identical. Index maps:
const _fSL = (cost_fix_local=1, cost_active_no3=2, cost_active_nh4=3, cost_nonmyc_no3=4,
    cost_nonmyc_nh4=5, npp_to_fixation=6, npp_to_active_nh4=7, npp_to_nonmyc_nh4=8,
    npp_to_active_no3=9, npp_to_nonmyc_no3=10, npp_frac_to_fixation=11,
    npp_frac_to_active_nh4=12, npp_frac_to_nonmyc_nh4=13, npp_frac_to_active_no3=14,
    npp_frac_to_nonmyc_no3=15, n_from_fixation=16, n_from_active_nh4=17,
    n_from_nonmyc_nh4=18, n_from_active_no3=19, n_from_nonmyc_no3=20, n_active_no3_vr=21,
    n_active_nh4_vr=22, n_nonmyc_no3_vr=23, n_nonmyc_nh4_vr=24)              # [p,j,K]
const _NFSL = 24
const _fPFT = (leafcn=1, grperc=2, FUN_fracfixers=3, c3psn=4, a_fix=5, b_fix=6, c_fix=7,
    s_fix=8, kc_nonmyc=9, kn_nonmyc=10, fun_cn_flex_a=11, fun_cn_flex_b=12,
    fun_cn_flex_c=13, season_decid=14, stress_decid=15)                     # [vt,K]
const _NFPFT = 15
const _fIV = (leafc=1, leafc_storage=2, leafn=3, leafn_storage=4, retransn=5, plantCN=6,
    availc=7, availc_pool=8, plant_ndemand_pool=9, n_passive_acc=10)        # [p,K]
const _NFIV = 10
const _fIJ = (crootfr=1, rootc_dens=2, n_passive_no3_vr=3, n_passive_nh4_vr=4)  # [p,j,K]
const _NFIJ = 4
const _fIPI = (permyc=1, kc_active=2, kn_active=3, litterfall_c_step=4,
    litterfall_n_step=5, n_passive_step=6)                                  # [p,istp,K]
const _NFIPI = 6
const _fOV = (Npassive=1, Nfix=2, retransn_to_npool=3, free_retransn_to_npool=4, Nretrans=5,
    Nactive_no3=6, Nactive_nh4=7, Necm_no3=8, Necm_nh4=9, Necm=10, Nam_no3=11, Nam_nh4=12,
    Nam=13, Nnonmyc_no3=14, Nnonmyc_nh4=15, Nnonmyc=16, plant_ndemand_retrans=17,
    Nuptake=18, Nactive=19, sminn_to_plant_fun=20, nuptake_npp_fraction=21, cost_nfix=22,
    cost_nactive=23, cost_nretrans=24, npp_Nactive_no3=25, npp_Nactive_nh4=26,
    npp_Nnonmyc_no3=27, npp_Nnonmyc_nh4=28, npp_Nactive=29, npp_Nnonmyc=30, npp_Nfix=31,
    npp_Nretrans=32, soilc_change=33, npp_burnedoff=34, npp_Nuptake=35, npp_growth=36)  # [p,K]
const _NFOV = 36  # OVR[p,j,1:2] = sminn_fun_no3_vr, sminn_fun_nh4_vr

# Pack the per-vt pftcon params into PFT[vt,K].
@kernel function _fun_pack_pft!(PFT, @Const(leafcn), @Const(grperc), @Const(FUN_fracfixers),
        @Const(c3psn), @Const(a_fix), @Const(b_fix), @Const(c_fix), @Const(s_fix),
        @Const(kc_nonmyc), @Const(kn_nonmyc), @Const(fun_cn_flex_a), @Const(fun_cn_flex_b),
        @Const(fun_cn_flex_c), @Const(season_decid), @Const(stress_decid))
    v = @index(Global)
    @inbounds begin
        PFT[v,1]=leafcn[v]; PFT[v,2]=grperc[v]; PFT[v,3]=FUN_fracfixers[v]; PFT[v,4]=c3psn[v]
        PFT[v,5]=a_fix[v]; PFT[v,6]=b_fix[v]; PFT[v,7]=c_fix[v]; PFT[v,8]=s_fix[v]
        PFT[v,9]=kc_nonmyc[v]; PFT[v,10]=kn_nonmyc[v]; PFT[v,11]=fun_cn_flex_a[v]
        PFT[v,12]=fun_cn_flex_b[v]; PFT[v,13]=fun_cn_flex_c[v]; PFT[v,14]=season_decid[v]
        PFT[v,15]=stress_decid[v]
    end
end

# Pack per-patch inputs into IV[p,K] and seed OV plant_ndemand_retrans (col 17).
@kernel function _fun_pack_iv!(IV, OV, @Const(mask), @Const(leafc), @Const(leafc_storage),
        @Const(leafn), @Const(leafn_storage), @Const(retransn), @Const(plantCN),
        @Const(availc), @Const(availc_pool), @Const(plant_ndemand_pool),
        @Const(n_passive_acc), @Const(plant_ndemand_retrans))
    p = @index(Global)
    @inbounds if mask[p]
        IV[p,1]=leafc[p]; IV[p,2]=leafc_storage[p]; IV[p,3]=leafn[p]; IV[p,4]=leafn_storage[p]
        IV[p,5]=retransn[p]; IV[p,6]=plantCN[p]; IV[p,7]=availc[p]; IV[p,8]=availc_pool[p]
        IV[p,9]=plant_ndemand_pool[p]; IV[p,10]=n_passive_acc[p]
        OV[p,17]=plant_ndemand_retrans[p]   # substep does OV[p,17] /= dt
    end
end

# Pack per-(p,j) inputs into IJ[p,j,K].
@kernel function _fun_pack_ij!(IJ, @Const(mask), @Const(crootfr), @Const(rootc_dens),
        @Const(n_passive_no3_vr), @Const(n_passive_nh4_vr), nlevdecomp::Int)
    p = @index(Global)
    @inbounds if mask[p]
        for j in 1:nlevdecomp
            IJ[p,j,1]=crootfr[p,j]; IJ[p,j,2]=rootc_dens[p,j]
            IJ[p,j,3]=n_passive_no3_vr[p,j]; IJ[p,j,4]=n_passive_nh4_vr[p,j]
        end
    end
end

# Pack per-(p,istp) inputs into IPI[p,istp,K].
@kernel function _fun_pack_ipi!(IPI, @Const(mask), @Const(permyc), @Const(kc_active),
        @Const(kn_active), @Const(litterfall_c_step), @Const(litterfall_n_step),
        @Const(n_passive_step), nstp::Int)
    p = @index(Global)
    @inbounds if mask[p]
        for istp in 1:nstp
            IPI[p,istp,1]=permyc[p,istp]; IPI[p,istp,2]=kc_active[p,istp]
            IPI[p,istp,3]=kn_active[p,istp]; IPI[p,istp,4]=litterfall_c_step[p,istp]
            IPI[p,istp,5]=litterfall_n_step[p,istp]; IPI[p,istp,6]=n_passive_step[p,istp]
        end
    end
end

# Unpack OV[p,K] / OVR[p,j,K] back into the cnveg flux fields.
@kernel function _fun_unpack!(@Const(mask), @Const(OV), @Const(OVR),
        Npassive, Nfix, retransn_to_npool, free_retransn_to_npool, Nretrans,
        Nactive_no3, Nactive_nh4, Necm_no3, Necm_nh4, Necm, Nam_no3, Nam_nh4, Nam,
        Nnonmyc_no3, Nnonmyc_nh4, Nnonmyc, plant_ndemand_retrans, Nuptake, Nactive,
        sminn_to_plant_fun, nuptake_npp_fraction, cost_nfix, cost_nactive, cost_nretrans,
        npp_Nactive_no3, npp_Nactive_nh4, npp_Nnonmyc_no3, npp_Nnonmyc_nh4, npp_Nactive,
        npp_Nnonmyc, npp_Nfix, npp_Nretrans, soilc_change, npp_burnedoff, npp_Nuptake,
        npp_growth, sminn_fun_no3_vr, sminn_fun_nh4_vr, nlevdecomp::Int)
    p = @index(Global)
    @inbounds if mask[p]
        Npassive[p]=OV[p,1]; Nfix[p]=OV[p,2]; retransn_to_npool[p]=OV[p,3]
        free_retransn_to_npool[p]=OV[p,4]; Nretrans[p]=OV[p,5]; Nactive_no3[p]=OV[p,6]
        Nactive_nh4[p]=OV[p,7]; Necm_no3[p]=OV[p,8]; Necm_nh4[p]=OV[p,9]; Necm[p]=OV[p,10]
        Nam_no3[p]=OV[p,11]; Nam_nh4[p]=OV[p,12]; Nam[p]=OV[p,13]; Nnonmyc_no3[p]=OV[p,14]
        Nnonmyc_nh4[p]=OV[p,15]; Nnonmyc[p]=OV[p,16]; plant_ndemand_retrans[p]=OV[p,17]
        Nuptake[p]=OV[p,18]; Nactive[p]=OV[p,19]; sminn_to_plant_fun[p]=OV[p,20]
        nuptake_npp_fraction[p]=OV[p,21]; cost_nfix[p]=OV[p,22]; cost_nactive[p]=OV[p,23]
        cost_nretrans[p]=OV[p,24]; npp_Nactive_no3[p]=OV[p,25]; npp_Nactive_nh4[p]=OV[p,26]
        npp_Nnonmyc_no3[p]=OV[p,27]; npp_Nnonmyc_nh4[p]=OV[p,28]; npp_Nactive[p]=OV[p,29]
        npp_Nnonmyc[p]=OV[p,30]; npp_Nfix[p]=OV[p,31]; npp_Nretrans[p]=OV[p,32]
        soilc_change[p]=OV[p,33]; npp_burnedoff[p]=OV[p,34]; npp_Nuptake[p]=OV[p,35]
        npp_growth[p]=OV[p,36]
        for j in 1:nlevdecomp
            sminn_fun_no3_vr[p,j]=OVR[p,j,1]; sminn_fun_nh4_vr[p,j]=OVR[p,j,2]
        end
    end
end

# Phase-2 substep solver (one thread per patch). Reads packed dense tensors, aliases
# each logical array to a @view slice (Metal-safe), runs the loop-carried ECM/AM ×
# fixer × level solver + retrans while-loop + flux tail on the thread's own patch.
@kernel function _fun_p2_kernel!(@Const(mask_soilp), @Const(column), @Const(ivt),
        @Const(PFT), @Const(IV), @Const(IJ), @Const(IPI), @Const(sminn_no3_layer_step),
        @Const(sminn_nh4_layer_step), @Const(t_soisno), @Const(dzsoi), SL, OV, OVR,
        dt, nlevdecomp::Int, smallValue, spval, npcropmin::Int,
        use_flexiblecn::Bool, use_matrixcn::Bool)
    p = @index(Global)
    @inbounds if mask_soilp[p]
        # pftcon (PFT[vt,K])
        leafcn=@view PFT[:,1]; grperc_pft=@view PFT[:,2]; FUN_fracfixers_pft=@view PFT[:,3]
        c3psn_pft=@view PFT[:,4]; a_fix_pft=@view PFT[:,5]; b_fix_pft=@view PFT[:,6]
        c_fix_pft=@view PFT[:,7]; s_fix_pft=@view PFT[:,8]; kc_nonmyc_pft=@view PFT[:,9]
        kn_nonmyc_pft=@view PFT[:,10]; fun_cn_flex_a_pft=@view PFT[:,11]
        fun_cn_flex_b_pft=@view PFT[:,12]; fun_cn_flex_c_pft=@view PFT[:,13]
        season_decid=@view PFT[:,14]; stress_decid=@view PFT[:,15]
        # per-patch inputs (IV[p,K])
        leafc=@view IV[:,1]; leafc_storage=@view IV[:,2]; leafn=@view IV[:,3]
        leafn_storage=@view IV[:,4]; retransn=@view IV[:,5]; plantCN=@view IV[:,6]
        availc=@view IV[:,7]; availc_pool=@view IV[:,8]; plant_ndemand_pool=@view IV[:,9]
        n_passive_acc=@view IV[:,10]
        # per-(p,j) inputs (IJ[p,j,K])
        crootfr=@view IJ[:,:,1]; rootc_dens=@view IJ[:,:,2]; n_passive_no3_vr=@view IJ[:,:,3]
        n_passive_nh4_vr=@view IJ[:,:,4]
        # per-(p,istp) inputs (IPI[p,istp,K])
        permyc=@view IPI[:,:,1]; kc_active=@view IPI[:,:,2]; kn_active=@view IPI[:,:,3]
        litterfall_c_step=@view IPI[:,:,4]; litterfall_n_step=@view IPI[:,:,5]
        n_passive_step=@view IPI[:,:,6]
        # per-(p,j) persistent scratch (SL[p,j,K])
        cost_fix_local=@view SL[:,:,1]; cost_active_no3=@view SL[:,:,2]
        cost_active_nh4=@view SL[:,:,3]; cost_nonmyc_no3=@view SL[:,:,4]
        cost_nonmyc_nh4=@view SL[:,:,5]; npp_to_fixation=@view SL[:,:,6]
        npp_to_active_nh4=@view SL[:,:,7]; npp_to_nonmyc_nh4=@view SL[:,:,8]
        npp_to_active_no3=@view SL[:,:,9]; npp_to_nonmyc_no3=@view SL[:,:,10]
        npp_frac_to_fixation=@view SL[:,:,11]; npp_frac_to_active_nh4=@view SL[:,:,12]
        npp_frac_to_nonmyc_nh4=@view SL[:,:,13]; npp_frac_to_active_no3=@view SL[:,:,14]
        npp_frac_to_nonmyc_no3=@view SL[:,:,15]; n_from_fixation=@view SL[:,:,16]
        n_from_active_nh4=@view SL[:,:,17]; n_from_nonmyc_nh4=@view SL[:,:,18]
        n_from_active_no3=@view SL[:,:,19]; n_from_nonmyc_no3=@view SL[:,:,20]
        n_active_no3_vr=@view SL[:,:,21]; n_active_nh4_vr=@view SL[:,:,22]
        n_nonmyc_no3_vr=@view SL[:,:,23]; n_nonmyc_nh4_vr=@view SL[:,:,24]
        # outputs (OV[p,K], OVR[p,j,K])
        Npassive=@view OV[:,1]; Nfix=@view OV[:,2]; retransn_to_npool=@view OV[:,3]
        free_retransn_to_npool=@view OV[:,4]; Nretrans=@view OV[:,5]; Nactive_no3=@view OV[:,6]
        Nactive_nh4=@view OV[:,7]; Necm_no3=@view OV[:,8]; Necm_nh4=@view OV[:,9]
        Necm=@view OV[:,10]; Nam_no3=@view OV[:,11]; Nam_nh4=@view OV[:,12]; Nam=@view OV[:,13]
        Nnonmyc_no3=@view OV[:,14]; Nnonmyc_nh4=@view OV[:,15]; Nnonmyc=@view OV[:,16]
        plant_ndemand_retrans=@view OV[:,17]; Nuptake=@view OV[:,18]; Nactive=@view OV[:,19]
        sminn_to_plant_fun=@view OV[:,20]; nuptake_npp_fraction=@view OV[:,21]
        cost_nfix=@view OV[:,22]; cost_nactive=@view OV[:,23]; cost_nretrans=@view OV[:,24]
        npp_Nactive_no3=@view OV[:,25]; npp_Nactive_nh4=@view OV[:,26]
        npp_Nnonmyc_no3=@view OV[:,27]; npp_Nnonmyc_nh4=@view OV[:,28]; npp_Nactive=@view OV[:,29]
        npp_Nnonmyc=@view OV[:,30]; npp_Nfix=@view OV[:,31]; npp_Nretrans=@view OV[:,32]
        soilc_change=@view OV[:,33]; npp_burnedoff=@view OV[:,34]; npp_Nuptake=@view OV[:,35]
        npp_growth=@view OV[:,36]
        sminn_fun_no3_vr=@view OVR[:,:,1]; sminn_fun_nh4_vr=@view OVR[:,:,2]

        c = column[p]
        vt = ivt[p] + 1
        T = eltype(OV)   # working precision (Float64 CPU / Float32 Metal)
        local_use_flexibleCN = use_flexiblecn
        excess_carbon_acc = zero(T)
        burned_off_carbon_local = zero(T)

        for j in 1:nlevdecomp
            sminn_fun_nh4_vr[p, j] = zero(T)
            sminn_fun_no3_vr[p, j] = zero(T)
        end

        plantCN_p = plantCN[p]

        # per-patch totals (were [p] scratch arrays) — thread-local scalars
        n_active_no3_acc_total = zero(T); n_active_nh4_acc_total = zero(T)
        n_nonmyc_no3_acc_total = zero(T); n_nonmyc_nh4_acc_total = zero(T)
        n_fix_acc_total = zero(T); n_retrans_acc_total = zero(T); free_Nretrans_local = zero(T)
        npp_active_no3_acc_total = zero(T); npp_active_nh4_acc_total = zero(T)
        npp_nonmyc_no3_acc_total = zero(T); npp_nonmyc_nh4_acc_total = zero(T)
        npp_fix_acc_total = zero(T); npp_retrans_acc_total = zero(T)
        n_ecm_no3_acc = zero(T); n_ecm_nh4_acc = zero(T); n_am_no3_acc = zero(T); n_am_nh4_acc = zero(T)

        for istp in ECM_STEP:AM_STEP
            sminn_no3_diff = zero(T); sminn_nh4_diff = zero(T)
            active_no3_limit1 = zero(T); active_nh4_limit1 = zero(T)

            for j in 1:nlevdecomp
                n_from_active_no3[p, j] = zero(T); n_from_active_nh4[p, j] = zero(T)
                n_from_nonmyc_no3[p, j] = zero(T); n_from_nonmyc_nh4[p, j] = zero(T)
                n_from_fixation[p, j] = zero(T)
            end

            # per-(p,istp) accumulators (were [p,istp] scratch) — reset each substep
            n_active_no3_acc = zero(T); n_active_nh4_acc = zero(T)
            n_nonmyc_no3_acc = zero(T); n_nonmyc_nh4_acc = zero(T)
            n_fix_acc = zero(T); n_retrans_acc = zero(T); free_nretrans_acc = zero(T)
            npp_active_no3_acc = zero(T); npp_active_nh4_acc = zero(T)
            npp_nonmyc_no3_acc = zero(T); npp_nonmyc_nh4_acc = zero(T)
            npp_fix_acc = zero(T); npp_retrans_acc = zero(T)
            nt_uptake = zero(T); npp_uptake = zero(T)

            for j in 1:nlevdecomp
                npp_to_active_no3[p, j] = zero(T); npp_to_active_nh4[p, j] = zero(T)
                npp_to_nonmyc_no3[p, j] = zero(T); npp_to_nonmyc_nh4[p, j] = zero(T)
                npp_to_fixation[p, j] = zero(T)
            end

            plant_ndemand_pool_step = plant_ndemand_pool[p] * permyc[p, istp]
            npp_remaining = availc_pool[p] * permyc[p, istp]

            for j in 1:nlevdecomp
                tc = t_soisno[c, j] - T(TFRZ)
                fixer = c3psn_pft[vt] == one(T) ? 1 : 0   # c3psn is a 0/1 flag; avoids GPU-hostile round(Int,·)
                cost_fix_local[p, j] = fun_cost_fix(fixer, a_fix_pft[vt], b_fix_pft[vt],
                    c_fix_pft[vt], T(BIG_COST), crootfr[p, j], s_fix_pft[vt], tc)
            end
            for j in 1:nlevdecomp
                rootc_dens_step = rootc_dens[p, j] * permyc[p, istp]
                cost_active_no3[p, j] = fun_cost_active(sminn_no3_layer_step[p, j, istp],
                    T(BIG_COST), kc_active[p, istp], kn_active[p, istp], rootc_dens_step,
                    crootfr[p, j], smallValue)
                cost_active_nh4[p, j] = fun_cost_active(sminn_nh4_layer_step[p, j, istp],
                    T(BIG_COST), kc_active[p, istp], kn_active[p, istp], rootc_dens_step,
                    crootfr[p, j], smallValue)
            end
            for j in 1:nlevdecomp
                rootc_dens_step = rootc_dens[p, j] * permyc[p, istp]
                cost_nonmyc_no3[p, j] = fun_cost_nonmyc(sminn_no3_layer_step[p, j, istp],
                    T(BIG_COST), kc_nonmyc_pft[vt], kn_nonmyc_pft[vt], rootc_dens_step,
                    crootfr[p, j], smallValue)
                cost_nonmyc_nh4[p, j] = fun_cost_nonmyc(sminn_nh4_layer_step[p, j, istp],
                    T(BIG_COST), kc_nonmyc_pft[vt], kn_nonmyc_pft[vt], rootc_dens_step,
                    crootfr[p, j], smallValue)
            end

            npp_remaining -= n_passive_step[p, istp] * plantCN_p

            for FIX in PLANTS_ARE_FIXING:PLANTS_NOT_FIXING
                if FIX == PLANTS_ARE_FIXING
                    fixerfrac = FUN_fracfixers_pft[vt]
                else
                    fixerfrac = one(T) - FUN_fracfixers_pft[vt]
                end
                npp_to_spend = npp_remaining * fixerfrac

                for j in 1:nlevdecomp
                    n_from_active_no3[p, j] = zero(T); n_from_active_nh4[p, j] = zero(T)
                    n_from_nonmyc_no3[p, j] = zero(T); n_from_nonmyc_nh4[p, j] = zero(T)
                end

                sum_n_acquired = zero(T)
                total_N_conductance = zero(T)
                for j in 1:nlevdecomp
                    total_N_conductance += one(T) / cost_active_no3[p, j] +
                                           one(T) / cost_active_nh4[p, j] +
                                           one(T) / cost_nonmyc_no3[p, j] +
                                           one(T) / cost_nonmyc_nh4[p, j]
                    if FIX == PLANTS_ARE_FIXING
                        total_N_conductance += one(T) / cost_fix_local[p, j]
                    end
                end

                for j in 1:nlevdecomp
                    npp_frac_to_active_nh4[p, j] = (one(T) / cost_active_nh4[p, j]) / total_N_conductance
                    npp_frac_to_nonmyc_nh4[p, j] = (one(T) / cost_nonmyc_nh4[p, j]) / total_N_conductance
                    npp_frac_to_active_no3[p, j] = (one(T) / cost_active_no3[p, j]) / total_N_conductance
                    npp_frac_to_nonmyc_no3[p, j] = (one(T) / cost_nonmyc_no3[p, j]) / total_N_conductance
                    if FIX == PLANTS_ARE_FIXING
                        npp_frac_to_fixation[p, j] = (one(T) / cost_fix_local[p, j]) / total_N_conductance
                    else
                        npp_frac_to_fixation[p, j] = zero(T)
                    end
                    ne_fix = FIX == PLANTS_ARE_FIXING ? npp_frac_to_fixation[p, j] / cost_fix_local[p, j] : zero(T)
                    ne_anh4 = npp_frac_to_active_nh4[p, j] / cost_active_nh4[p, j]
                    ne_nnh4 = npp_frac_to_nonmyc_nh4[p, j] / cost_nonmyc_nh4[p, j]
                    ne_ano3 = npp_frac_to_active_no3[p, j] / cost_active_no3[p, j]
                    ne_nno3 = npp_frac_to_nonmyc_no3[p, j] / cost_nonmyc_no3[p, j]
                    sum_n_acquired += ne_anh4 + ne_nnh4 + ne_ano3 + ne_nno3
                    if FIX == PLANTS_ARE_FIXING
                        sum_n_acquired += ne_fix
                    end
                end

                total_N_resistance = one(T) / sum_n_acquired

                if leafc[p] > zero(T) && litterfall_n_step[p, istp] * fixerfrac > zero(T) && vt < npcropmin
                    rt_result = fun_retranslocation(dt, npp_to_spend,
                        litterfall_c_step[p, istp] * fixerfrac,
                        litterfall_n_step[p, istp] * fixerfrac, total_N_resistance,
                        leafcn[vt], grperc_pft[vt], plantCN_p)
                    total_c_spent_retrans     = rt_result.total_c_spent_retrans
                    total_c_accounted_retrans = rt_result.total_c_accounted_retrans
                    free_n_retrans_val        = rt_result.free_n_retrans
                    paid_for_n_retrans        = rt_result.paid_for_n_retrans
                else
                    total_c_accounted_retrans = zero(T); total_c_spent_retrans = zero(T)
                    paid_for_n_retrans = zero(T); free_n_retrans_val = zero(T)
                end

                npp_to_spend -= total_c_spent_retrans + total_c_accounted_retrans
                npp_retrans_acc += total_c_spent_retrans
                n_retrans_acc   += paid_for_n_retrans
                free_nretrans_acc += free_n_retrans_val

                if plant_ndemand_pool_step > zero(T)
                    if local_use_flexibleCN
                        if leafn[p] == zero(T)
                            delta_CN = fun_cn_flex_c_pft[vt]
                        else
                            delta_CN = (leafc[p] + leafc_storage[p]) /
                                       (leafn[p] + leafn_storage[p]) - leafcn[vt]
                        end
                        frac_ideal_C_use = max(zero(T), one(T) - (total_N_resistance - fun_cn_flex_a_pft[vt]) /
                                           fun_cn_flex_b_pft[vt])
                        if delta_CN < zero(T)
                            frac_ideal_C_use += (one(T) - frac_ideal_C_use) *
                                                min(one(T), delta_CN / fun_cn_flex_c_pft[vt])
                        end
                        if delta_CN > zero(T) && frac_ideal_C_use < one(T)
                            frac_ideal_C_use += T(0.5) * (one(T) * delta_CN / fun_cn_flex_c_pft[vt])
                        end
                        frac_ideal_C_use = max(min(one(T), frac_ideal_C_use), T(0.5))
                    else
                        frac_ideal_C_use = one(T)
                    end

                    excess_carbon = npp_to_spend * (one(T) - frac_ideal_C_use)
                    if excess_carbon * (one(T) + grperc_pft[vt]) > npp_to_spend
                        excess_carbon = npp_to_spend / (one(T) + grperc_pft[vt])
                    end
                    excess_carbon_acc += excess_carbon
                    npp_to_spend -= excess_carbon * (one(T) + grperc_pft[vt])

                    dnpp = npp_to_spend / ((one(T) + grperc_pft[vt]) * (plantCN_p / total_N_resistance) + one(T))
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
                            npp_to_fixation[p, j] = zero(T)
                        end
                        n_from_active_nh4[p, j] = npp_to_active_nh4[p, j] / cost_active_nh4[p, j]
                        n_from_nonmyc_nh4[p, j] = npp_to_nonmyc_nh4[p, j] / cost_nonmyc_nh4[p, j]
                        n_from_active_no3[p, j] = npp_to_active_no3[p, j] / cost_active_no3[p, j]
                        n_from_nonmyc_no3[p, j] = npp_to_nonmyc_no3[p, j] / cost_nonmyc_no3[p, j]
                        if FIX == PLANTS_ARE_FIXING
                            n_from_fixation[p, j] = npp_to_fixation[p, j] / cost_fix_local[p, j]
                        else
                            n_from_fixation[p, j] = zero(T)
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

                        npp_to_spend -= C_spent + N_acquired * plantCN_p * (one(T) + grperc_pft[vt])
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

                    if npp_to_spend >= T(1.0e-13)
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
            if retransn[p] > zero(T)
                free_retransn_to_npool[p] = free_Nretrans_local / dt
            else
                free_retransn_to_npool[p] = zero(T)
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

        if availc[p] > zero(T)
            nuptake_npp_fraction[p] = npp_Nuptake[p] / availc[p]
        else
            nuptake_npp_fraction[p] = spval
        end
        if npp_Nfix[p] > zero(T)
            cost_nfix[p] = Nfix[p] / npp_Nfix[p]
        else
            cost_nfix[p] = spval
        end
        if npp_Nactive[p] > zero(T)
            cost_nactive[p] = Nactive[p] / npp_Nactive[p]
        else
            cost_nactive[p] = spval
        end
        if npp_Nretrans[p] > zero(T)
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

    # ---- Phase-2 substep: pack ~103 arrays into dense tensors (Metal ≤31 buffers),
    #      run the solver, unpack. @view slices alias each logical array in the kernel
    #      so the solver body is byte-identical to the host loop. ----
    npft_ct = length(leafcn)
    npb = last(bounds_p)
    PFT_t = similar(availc_pool, FT, npft_ct, _NFPFT)
    IV_t  = similar(availc_pool, FT, npb, _NFIV);            fill!(IV_t, zero(FT))
    IJ_t  = similar(availc_pool, FT, npb, nlevdecomp, _NFIJ)
    IPI_t = similar(availc_pool, FT, npb, NSTP, _NFIPI)
    SL_t  = similar(availc_pool, FT, npb, nlevdecomp, _NFSL); fill!(SL_t, zero(FT))
    OV_t  = similar(availc_pool, FT, npb, _NFOV);             fill!(OV_t, zero(FT))
    OVR_t = similar(availc_pool, FT, npb, nlevdecomp, 2);     fill!(OVR_t, zero(FT))

    _launch!(_fun_pack_pft!, PFT_t, leafcn, grperc_pft, FUN_fracfixers_pft, c3psn_pft,
        a_fix_pft, b_fix_pft, c_fix_pft, s_fix_pft, kc_nonmyc_pft, kn_nonmyc_pft,
        fun_cn_flex_a_pft, fun_cn_flex_b_pft, fun_cn_flex_c_pft, season_decid, stress_decid;
        ndrange=npft_ct)
    _launch!(_fun_pack_iv!, IV_t, OV_t, mask_soilp, cnveg_cs.leafc_patch,
        cnveg_cs.leafc_storage_patch, cnveg_ns.leafn_patch, cnveg_ns.leafn_storage_patch,
        cnveg_ns.retransn_patch, cnveg_state.plantCN_patch, cnveg_cf.availc_patch,
        availc_pool, plant_ndemand_pool, n_passive_acc, cnveg_nf.plant_ndemand_retrans_patch;
        ndrange=npb)
    _launch!(_fun_pack_ij!, IJ_t, mask_soilp, soilstate.crootfr_patch, rootc_dens,
        n_passive_no3_vr, n_passive_nh4_vr, nlevdecomp; ndrange=npb)
    _launch!(_fun_pack_ipi!, IPI_t, mask_soilp, permyc, kc_active, kn_active,
        litterfall_c_step, litterfall_n_step, n_passive_step, NSTP; ndrange=npb)

    _launch!(_fun_p2_kernel!, mask_soilp, patch.column, ivt, PFT_t, IV_t, IJ_t, IPI_t,
        sminn_no3_layer_step, sminn_nh4_layer_step, temperature.t_soisno_col,
        dzsoi_decomp_vals, SL_t, OV_t, OVR_t, dt, nlevdecomp, smallValue, spval,
        npcropmin, use_flexiblecn, use_matrixcn)

    _launch!(_fun_unpack!, mask_soilp, OV_t, OVR_t, cnveg_nf.Npassive_patch,
        cnveg_nf.Nfix_patch, cnveg_nf.retransn_to_npool_patch,
        cnveg_nf.free_retransn_to_npool_patch, cnveg_nf.Nretrans_patch,
        cnveg_nf.Nactive_no3_patch, cnveg_nf.Nactive_nh4_patch, cnveg_nf.Necm_no3_patch,
        cnveg_nf.Necm_nh4_patch, cnveg_nf.Necm_patch, cnveg_nf.Nam_no3_patch,
        cnveg_nf.Nam_nh4_patch, cnveg_nf.Nam_patch, cnveg_nf.Nnonmyc_no3_patch,
        cnveg_nf.Nnonmyc_nh4_patch, cnveg_nf.Nnonmyc_patch,
        cnveg_nf.plant_ndemand_retrans_patch, cnveg_nf.Nuptake_patch, cnveg_nf.Nactive_patch,
        cnveg_nf.sminn_to_plant_fun_patch, cnveg_nf.nuptake_npp_fraction_patch,
        cnveg_nf.cost_nfix_patch, cnveg_nf.cost_nactive_patch, cnveg_nf.cost_nretrans_patch,
        cnveg_cf.npp_Nactive_no3_patch, cnveg_cf.npp_Nactive_nh4_patch,
        cnveg_cf.npp_Nnonmyc_no3_patch, cnveg_cf.npp_Nnonmyc_nh4_patch,
        cnveg_cf.npp_Nactive_patch, cnveg_cf.npp_Nnonmyc_patch, cnveg_cf.npp_Nfix_patch,
        cnveg_cf.npp_Nretrans_patch, cnveg_cf.soilc_change_patch, cnveg_cf.npp_burnedoff_patch,
        cnveg_cf.npp_Nuptake_patch, cnveg_cf.npp_growth_patch,
        cnveg_nf.sminn_to_plant_fun_no3_vr_patch, cnveg_nf.sminn_to_plant_fun_nh4_vr_patch,
        nlevdecomp; ndrange=npb)

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
