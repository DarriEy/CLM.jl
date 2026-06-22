# PRTAllometricCNPMod.jl
# Julia port of FATES src/fates/parteh/PRTAllometricCNPMod.F90
#
# The Carbon-Nitrogen-Phosphorus (CNP) Prioritized Allometric Allocation
# hypothesis — the most complex PARTEH allocation strategy. This is a concrete
# subtype of the merged `AbstractPRTVartypes` (PRTGenericMod.jl) implementing
# the deferred generics `DailyPRT!`, `FastPRT!`, `GetNutrientTarget` and
# `GetCoordVal`, plus the CNP-specific extended methods (CNPPrioritizedReplacement,
# CNPStatureGrowth, EstimateGrowthNC, CNPAdjustFRootTargets, CNPAllocateRemainder,
# GetDeficit, TrimFineRoot).
#
# Ryan Knox Aug 2018 (Fortran original).
#
# Translation notes:
# - `fates_r8`     -> Float64; `fates_int` -> Int.
# - The concrete Fortran type `cnp_allom_prt_vartypes <: prt_vartypes` becomes a
#   Julia mutable struct `cnp_allom_prt_vartypes <: AbstractPRTVartypes` that
#   carries the SAME fields as the generic carrier `prt_vartypes` (variables,
#   bc_inout, bc_in, bc_out, ode_opt_step). All the generic PRTGenericMod
#   functions accept `::AbstractPRTVartypes` and so work on this subtype directly.
# - Type-bound procedures (`this%DailyPRT`, `this%CNPStatureGrowth`, ...) become
#   plain Julia functions taking `this::cnp_allom_prt_vartypes` first.
# - Boundary-condition pointers (`this%bc_inout(id)%rval`) are `Ref`s; we read /
#   write through `[]` (e.g. `this.bc_inout[id].rval[]`).
# - The module-level singletons `n_uptake_mode`, `p_uptake_mode`,
#   `regeneration_model` come from the merged EDParamsMod via `ed_params()`.
# - `prt_global` is the merged Ref singleton; the module-level Fortran
#   `prt_global_acnp` is mirrored by a Julia Ref (`prt_global_acnp`).
# - Generic functions `GetState`, `GetTurnover`, `GetNetAlloc`, `SetState!`,
#   `StorageNutrientTarget`, the allometry `(value,deriv)` routines, and the ODE
#   integrators (Euler/RKF45) are reused from the merged modules verbatim.
#
# Deps: PRTGenericMod, FatesAllometryMod, PRTParametersMod (prt_params),
# EDParamsMod (ed_params), FatesConstantsMod, FatesIntegratorsMod, FatesGlobals.

# -------------------------------------------------------------------------------------
# State variable indices (organ x element). These are the *local* variable ids of
# this hypothesis (1:num_vars); they match the registration order below.
# -------------------------------------------------------------------------------------
const leaf_c_id   = 1       # leaf carbon index
const fnrt_c_id   = 2       # fine-root carbon index
const sapw_c_id   = 3       # sapwood carbon index
const store_c_id  = 4       # storage carbon index
const repro_c_id  = 5       # reproductive carbon index
const struct_c_id = 6       # structural carbon index

const leaf_n_id   = 7
const fnrt_n_id   = 8
const sapw_n_id   = 9
const store_n_id  = 10
const repro_n_id  = 11
const struct_n_id = 12

const leaf_p_id   = 13
const fnrt_p_id   = 14
const sapw_p_id   = 15
const store_p_id  = 16
const repro_p_id  = 17
const struct_p_id = 18

# Total number of state variables
const acnp_num_vars = 18

# Global identifiers for the two stoichiometry "modes"
const stoich_growth_min = 1    # stoichiometry associated with the minimum needed for growth
const stoich_max        = 2    # deprecated (no reasonable hypothesis yet)

# The ordered list of organs used in this module
const acnp_num_organs = 6

# Converting from local (1:num_organs) to global organ id (PRTGenericMod)
const l2g_organ_list = [leaf_organ, fnrt_organ, sapw_organ,
                        store_organ, repro_organ, struct_organ]

# Local indices for the integrable state array (pools + dbh)
const leaf_id   = 1
const fnrt_id   = 2
const sapw_id   = 3
const store_id  = 4
const repro_id  = 5
const struct_id = 6
const dbh_id    = 7
const num_intgr_vars = 7

# -------------------------------------------------------------------------------------
# Input/Output Boundary Indices (in-out)
# -------------------------------------------------------------------------------------
const acnp_bc_inout_id_dbh         = 1  # Plant DBH
const acnp_bc_inout_id_resp_excess = 2  # Respiration of excess storage
const acnp_bc_inout_id_l2fr        = 3  # leaf-2-fineroot scalar (dynamic with CNP)
const acnp_bc_inout_id_netdn       = 4  # net daily NH4 input BC
const acnp_bc_inout_id_netdp       = 5  # net daily P input BC
const acnp_bc_inout_id_cx_int      = 6  # EMA log storage ratio max(N,P)/C (integral)
const acnp_bc_inout_id_cx0         = 7  # previous step's log storage ratio max(N,P)/C
const acnp_bc_inout_id_emadcxdt    = 8  # EMA log storage ratio derivative
const acnp_num_bc_inout            = 8

# Input only Boundary Indices
const acnp_bc_in_id_pft      = 1   # PFT input BC
const acnp_bc_in_id_ctrim    = 2   # canopy trimming function
const acnp_bc_in_id_lstat    = 3   # phenology status
const acnp_bc_in_id_netdc    = 4   # net daily C input BC
const acnp_bc_in_id_nc_repro = 5
const acnp_bc_in_id_pc_repro = 6
const acnp_bc_in_id_cdamage  = 7   # crowndamage input BC
const acnp_bc_in_id_efleaf   = 8   # Leaf elongation factor
const acnp_bc_in_id_effnrt   = 9   # Fine-root "elongation factor"
const acnp_bc_in_id_efstem   = 10  # Stem "elongation factor"
const acnp_num_bc_in         = 10

# Output Boundary Indices
const acnp_bc_out_id_cefflux = 1  # Daily exudation of C  [kg]
const acnp_bc_out_id_nefflux = 2  # Daily exudation of N  [kg]
const acnp_bc_out_id_pefflux = 3  # Daily exudation of P  [kg]
const acnp_bc_out_id_limiter = 4  # The limiting element flag
const acnp_num_bc_out        = 4

# Indices for parameters passed to the integrator
const intgr_parm_ctrim   = 1
const intgr_parm_pft     = 2
const intgr_parm_l2fr    = 3
const intgr_parm_cdamage = 4
const intgr_parm_efleaf  = 5
const intgr_parm_effnrt  = 6
const intgr_parm_efstem  = 7
const num_intgr_parm     = 7

# Leaf coordinate (only one growing position by definition: the youngest)
const icd = 1

# Excess-carbon disposition methods.
const exude_c_store_overflow  = 1
const retain_c_store_overflow = 2
const burn_c_store_overflow   = 3
const store_c_overflow = burn_c_store_overflow

# Growth-limitation flags.
const cnp_limited = 0
const c_limited   = 1
const n_limited   = 2
const p_limited   = 3

# Growth limitation strategy.
const grow_lim_conly = 1   # Just use C to decide stature on this step
const grow_lim_estNP = 2   # Estimate equivalent C from N and P
const grow_lim_type  = grow_lim_estNP

# Following growth, optionally prioritize balanced CNP for reproductive tissue.
const prioritize_repro_nutr_growth = true

# Unrestricted contraction of fine-roots when l2fr drops.
const use_unrestricted_contraction = true

const acnp_debug = false

# =====================================================================================
# Concrete CNP hypothesis type
# =====================================================================================

"""
    cnp_allom_prt_vartypes

Concrete PARTEH state carrier for the Carbon-Nitrogen-Phosphorus prioritized
allometric allocation hypothesis. Subtypes `AbstractPRTVartypes` and carries the
same SoA state/boundary fields as the generic `prt_vartypes` carrier (so all the
generic PRTGenericMod functions apply). The CNP-specific allocation logic is
implemented as methods dispatching on this type.
"""
Base.@kwdef mutable struct cnp_allom_prt_vartypes <: AbstractPRTVartypes
    variables::Vector{prt_vartype} = prt_vartype[]
    bc_inout::Vector{prt_bctype}   = prt_bctype[]
    bc_in::Vector{prt_bctype}      = prt_bctype[]
    bc_out::Vector{prt_bctype}     = prt_bctype[]
    ode_opt_step::Float64          = 0.0
end

# The module-level instance of the mapping table / variable definitions for this
# hypothesis (mirrors Fortran `prt_global_acnp`). Allocated once per node.
const prt_global_acnp = Ref{Union{Nothing,prt_global_type}}(nothing)

# =====================================================================================
# InitPRTGlobalAllometricCNP — set up the global descriptor for this hypothesis
# =====================================================================================

"""
    InitPRTGlobalAllometricCNP()

Initialize and populate the general mapping table that organizes the specific
state variables of this hypothesis into the pre-ordained organ/element groups, so
they can be used to inform the rest of the model. Sets the module-level
`prt_global_acnp` AND points the merged `prt_global[]` singleton at it (the rest
of PARTEH dispatches through `prt_global`).
"""
function InitPRTGlobalAllometricCNP()
    g = prt_global_type()
    g.state_descriptor = [state_descriptor_type() for _ in 1:acnp_num_vars]

    g.hyp_name = "Allometric Flexible C+N+P"
    g.hyp_id   = prt_cnp_flex_allom_hyp

    ZeroGlobal!(g)

    # The number of leaf age classes is the size of leaf_long's second dim.
    nleafage = size(prt_params.leaf_long, 2)

    if nleafage > max_nleafage
        fates_endrun("The allometric CNP PARTEH hypothesis sets a maximum number " *
                     "of leaf age classes ($(max_nleafage)) for scratch space; the " *
                     "model wants $(nleafage). Increase max_nleafage in PRTGenericMod.")
    end

    RegisterVarInGlobal!(g, leaf_c_id,   "Leaf Carbon",       "leaf_c",   leaf_organ,   carbon12_element,   nleafage)
    RegisterVarInGlobal!(g, fnrt_c_id,   "Fine Root Carbon",  "fnrt_c",   fnrt_organ,   carbon12_element,   icd)
    RegisterVarInGlobal!(g, sapw_c_id,   "Sapwood Carbon",    "sapw_c",   sapw_organ,   carbon12_element,   icd)
    RegisterVarInGlobal!(g, store_c_id,  "Storage Carbon",    "store_c",  store_organ,  carbon12_element,   icd)
    RegisterVarInGlobal!(g, struct_c_id, "Structural Carbon", "struct_c", struct_organ, carbon12_element,   icd)
    RegisterVarInGlobal!(g, repro_c_id,  "Reproductive Carbon","repro_c", repro_organ,  carbon12_element,   icd)

    RegisterVarInGlobal!(g, leaf_n_id,   "Leaf Nitrogen",       "leaf_n",   leaf_organ,   nitrogen_element, nleafage)
    RegisterVarInGlobal!(g, fnrt_n_id,   "Fine Root Nitrogen",  "fnrt_n",   fnrt_organ,   nitrogen_element, icd)
    RegisterVarInGlobal!(g, sapw_n_id,   "Sapwood Nitrogen",    "sapw_n",   sapw_organ,   nitrogen_element, icd)
    RegisterVarInGlobal!(g, store_n_id,  "Storage Nitrogen",    "store_n",  store_organ,  nitrogen_element, icd)
    RegisterVarInGlobal!(g, struct_n_id, "Structural Nitrogen", "struct_n", struct_organ, nitrogen_element, icd)
    RegisterVarInGlobal!(g, repro_n_id,  "Reproductive Nitrogen","repro_n", repro_organ,  nitrogen_element, icd)

    RegisterVarInGlobal!(g, leaf_p_id,   "Leaf Phosphorus",       "leaf_p",   leaf_organ,   phosphorus_element, nleafage)
    RegisterVarInGlobal!(g, fnrt_p_id,   "Fine Root Phosphorus",  "fnrt_p",   fnrt_organ,   phosphorus_element, icd)
    RegisterVarInGlobal!(g, sapw_p_id,   "Sapwood Phosphorus",    "sapw_p",   sapw_organ,   phosphorus_element, icd)
    RegisterVarInGlobal!(g, store_p_id,  "Storage Phosphorus",    "store_p",  store_organ,  phosphorus_element, icd)
    RegisterVarInGlobal!(g, struct_p_id, "Structural Phosphorus", "struct_p", struct_organ, phosphorus_element, icd)
    RegisterVarInGlobal!(g, repro_p_id,  "Reproductive Phosphorus","repro_p", repro_organ,  phosphorus_element, icd)

    g.num_bc_in    = acnp_num_bc_in
    g.num_bc_out   = acnp_num_bc_out
    g.num_bc_inout = acnp_num_bc_inout
    g.num_vars     = acnp_num_vars

    prt_global_acnp[] = g
    prt_global[]      = g

    return nothing
end

# =====================================================================================
# Small helpers
# =====================================================================================

"""
    SafeLog(val) -> Float64

The log functions used to transform storage ratios need not be large; a ratio of
10 already sends a strong signal to the root adaptation algorithm. We bound the
argument to [1e-3, 1e3] to prevent numerical overflow/underflow.
"""
function SafeLog(val::Real)
    safelog_min = 0.001
    safelog_max = 1000.0
    return log(max(safelog_min, min(safelog_max, val)))
end

# Convenience BC accessors (read/write the Ref payloads).
@inline _rin(this, id)  = this.bc_in[id].rval[]
@inline _iin(this, id)  = this.bc_in[id].ival[]
@inline _rio(this, id)  = this.bc_inout[id].rval[]
@inline _iio(this, id)  = this.bc_inout[id].ival[]
@inline _rout(this, id) = this.bc_out[id].rval[]
@inline _iout(this, id) = this.bc_out[id].ival[]
@inline _set_rio!(this, id, v) = (this.bc_inout[id].rval[] = v)
@inline _set_rout!(this, id, v) = (this.bc_out[id].rval[] = v)
@inline _set_iout!(this, id, v) = (this.bc_out[id].ival[] = v)

# =====================================================================================
# GetDeficit
# =====================================================================================

"""
    GetDeficit(this, element_id, organ_id, target_m) -> deficit_m

The mass needed to bring (organ, element) up to `target_m`. For carbon the deficit
is taken against the sum over all positions; for nutrients, against position 1.
"""
function GetDeficit(this::cnp_allom_prt_vartypes, element_id::Integer,
                    organ_id::Integer, target_m::Real)
    g = get_prt_global()
    i_var = sp_organ_map_get(g, organ_id, element_id)
    if element_id == carbon12_element
        return target_m - sum(this.variables[i_var].val)
    else
        return target_m - this.variables[i_var].val[1]
    end
end

# =====================================================================================
# GetNutrientTarget (the CNP override of the deferred generic)
# =====================================================================================

"""
    GetNutrientTarget(this::cnp_allom_prt_vartypes, element_id, organ_id [, stoich_mode]) -> target_m

Target nutrient mass [kg] for the given organ. Storage targets are a fraction of
the maximum nutrient that can be bound in non-reproductive tissues; reproductive
targets follow the seed N:C / P:C ratios; all other organs follow the organ's
growth-minimum stoichiometry.
"""
function GetNutrientTarget(this::cnp_allom_prt_vartypes, element_id::Integer,
                           organ_id::Integer,
                           stoich_mode::Union{Nothing,Integer}=nothing)
    g = get_prt_global()

    dbh         = _rio(this, acnp_bc_inout_id_dbh)
    canopy_trim = _rin(this, acnp_bc_in_id_ctrim)
    ipft        = _iin(this, acnp_bc_in_id_pft)
    elongf_leaf = _rin(this, acnp_bc_in_id_efleaf)
    elongf_fnrt = _rin(this, acnp_bc_in_id_effnrt)
    elongf_stem = _rin(this, acnp_bc_in_id_efstem)
    i_cvar      = sp_organ_map_get(g, organ_id, carbon12_element)
    l2fr        = _rio(this, acnp_bc_inout_id_l2fr)
    nc_repro    = _rin(this, acnp_bc_in_id_nc_repro)
    pc_repro    = _rin(this, acnp_bc_in_id_pc_repro)
    crown_damage = _iin(this, acnp_bc_in_id_cdamage)

    target_m = 0.0

    if organ_id == store_organ

        leaf_c_target, _ = bleaf(dbh, ipft, crown_damage, canopy_trim, elongf_leaf)
        fnrt_c_target, _ = bfineroot(dbh, ipft, canopy_trim, l2fr, elongf_fnrt)
        _, sapw_c_target, _ = bsap_allom(dbh, ipft, crown_damage, canopy_trim, elongf_stem)
        agw_c_target, _  = bagw_allom(dbh, ipft, crown_damage, elongf_stem)
        bgw_c_target, _  = bbgw_allom(dbh, ipft, elongf_stem)
        struct_c_target, _ = bdead_allom(agw_c_target, bgw_c_target, sapw_c_target, ipft)

        # Target for storage is a fraction of the sum target of all
        # non-reproductive organs.
        if element_id == nitrogen_element
            target_m = StorageNutrientTarget(ipft, element_id,
                leaf_c_target  * prt_params.nitr_stoich_p1[ipft, prt_params.organ_param_id[leaf_organ]],
                fnrt_c_target  * prt_params.nitr_stoich_p1[ipft, prt_params.organ_param_id[fnrt_organ]],
                sapw_c_target  * prt_params.nitr_stoich_p1[ipft, prt_params.organ_param_id[sapw_organ]],
                struct_c_target* prt_params.nitr_stoich_p1[ipft, prt_params.organ_param_id[struct_organ]])
        else
            target_m = StorageNutrientTarget(ipft, element_id,
                leaf_c_target  * prt_params.phos_stoich_p1[ipft, prt_params.organ_param_id[leaf_organ]],
                fnrt_c_target  * prt_params.phos_stoich_p1[ipft, prt_params.organ_param_id[fnrt_organ]],
                sapw_c_target  * prt_params.phos_stoich_p1[ipft, prt_params.organ_param_id[sapw_organ]],
                struct_c_target* prt_params.phos_stoich_p1[ipft, prt_params.organ_param_id[struct_organ]])
        end

        # This is only called during phase 3 (remainder) and allows some overflow.
        if stoich_mode == stoich_max
            target_m = target_m * (1.0 + prt_params.store_ovrflw_frac[ipft])
        end

    elseif organ_id == repro_organ

        target_c = this.variables[i_cvar].val[1]
        if element_id == nitrogen_element
            target_m = target_c * nc_repro
        else
            target_m = target_c * pc_repro
        end

    else

        if stoich_mode === nothing
            fates_endrun("Must specify if nutrient target is growthmin or max " *
                         "for non-reproductive and non-storage organs")
        end

        # We always want the first index: for non-leaves it is the only index, and
        # for leaves it is the newly-growing index.
        target_c = this.variables[i_cvar].val[1]

        if stoich_mode == stoich_growth_min
            if element_id == nitrogen_element
                target_m = target_c * prt_params.nitr_stoich_p1[ipft, prt_params.organ_param_id[organ_id]]
            else
                target_m = target_c * prt_params.phos_stoich_p1[ipft, prt_params.organ_param_id[organ_id]]
            end
        else
            # stoich_max (or any other) is unsupported for these organs.
            fates_endrun("invalid stoichiometry mode specified while getting " *
                         "nutrient targets: stoich_mode=$(stoich_mode)")
        end
    end

    return target_m
end

# =====================================================================================
# GetCoordVal — support for variables with multiple discrete positions
# =====================================================================================

"""
    GetCoordVal(this::cnp_allom_prt_vartypes, organ_id, element_id) -> Float64

Return the value at the (single) growing coordinate for an (organ, element). For
this hypothesis all allocation occurs in position 1 (the youngest leaf bin / the
sole pool for other organs), so this returns `val[1]`.
"""
function GetCoordVal(this::cnp_allom_prt_vartypes, organ_id::Integer, element_id::Integer)
    g = get_prt_global()
    i_var = sp_organ_map_get(g, organ_id, element_id)
    return this.variables[i_var].val[1]
end

# =====================================================================================
# ProportionalNutrAllocation — module function (not type-bound in Fortran)
# =====================================================================================

"""
    ProportionalNutrAllocation!(this, deficit_m, gain_m, element_id, list) -> gain_m

Allocate the nutrient `gain_m` to a set of organs (`list`, global organ ids) in
proportion to their deficits (`deficit_m`, modified in place). Returns the updated
`gain_m` (caller must reassign). Mirrors the Fortran in-out semantics on `gain_m`.
"""
function ProportionalNutrAllocation!(this::cnp_allom_prt_vartypes,
                                     deficit_m::AbstractVector{Float64},
                                     gain_m::Float64, element_id::Integer,
                                     list::AbstractVector{<:Integer})
    g = get_prt_global()
    n_org = length(list)

    sum_deficit = 0.0
    for i in 1:n_org
        sum_deficit += max(0.0, deficit_m[i])
    end

    if sum_deficit > nearzero
        sum_flux = min(gain_m, sum_deficit)
        for i in 1:n_org
            i_org = list[i]
            flux  = sum_flux * max(0.0, deficit_m[i]) / sum_deficit
            i_var = sp_organ_map_get(g, i_org, element_id)
            this.variables[i_var].val[1] += flux
            deficit_m[i] -= flux
            gain_m       -= flux
        end
    end

    if acnp_debug && gain_m < -calloc_abs_error
        fates_endrun("Negative nutrient gain during proportional allocation: " *
                     "gain=$(gain_m) element=$(element_id)")
    end

    return gain_m
end

# =====================================================================================
# AllomCNPGrowthDeriv — the ODE derivative function for stature growth
# =====================================================================================

"""
    AllomCNPGrowthDeriv(l_state_array, l_state_mask, cbalance, intgr_params) -> dCdx

Derivatives of the carbon pools (and dbh) with respect to the carbon balance
spent on stature growth. Allometry-based: dC_pool/dx is the pool's share of the
total per-dbh carbon cost, with a fraction directed to reproduction. The masked
pools are the ones currently on allometry. `cbalance` is unused (matches Fortran
signature for the integrator's deriv-function contract).
"""
function AllomCNPGrowthDeriv(l_state_array::AbstractVector, l_state_mask::AbstractVector{Bool},
                             cbalance::Real, intgr_params::AbstractVector)
    dbh      = l_state_array[dbh_id]
    mask_dbh    = l_state_mask[dbh_id]
    mask_leaf   = l_state_mask[leaf_id]
    mask_fnrt   = l_state_mask[fnrt_id]
    mask_sapw   = l_state_mask[sapw_id]
    mask_store  = l_state_mask[store_id]
    mask_struct = l_state_mask[struct_id]
    mask_repro  = l_state_mask[repro_id]

    canopy_trim  = intgr_params[intgr_parm_ctrim]
    ipft         = Int(intgr_params[intgr_parm_pft])
    l2fr         = intgr_params[intgr_parm_l2fr]
    crown_damage = Int(intgr_params[intgr_parm_cdamage])
    elongf_leaf  = intgr_params[intgr_parm_efleaf]
    elongf_fnrt  = intgr_params[intgr_parm_effnrt]
    elongf_stem  = intgr_params[intgr_parm_efstem]

    leaf_c_target,  leaf_dcdd_target  = bleaf(dbh, ipft, crown_damage, canopy_trim, elongf_leaf)
    fnrt_c_target,  fnrt_dcdd_target  = bfineroot(dbh, ipft, canopy_trim, l2fr, elongf_fnrt)
    _, sapw_c_target, sapw_dcdd_target = bsap_allom(dbh, ipft, crown_damage, canopy_trim, elongf_stem)
    agw_c_target,   agw_dcdd_target   = bagw_allom(dbh, ipft, crown_damage, elongf_stem)
    bgw_c_target,   bgw_dcdd_target   = bbgw_allom(dbh, ipft, elongf_stem)
    struct_c_target, struct_dcdd_target = bdead_allom(agw_c_target, bgw_c_target, sapw_c_target, ipft;
                                                      dbagwdd=agw_dcdd_target, dbbgwdd=bgw_dcdd_target,
                                                      dbsapdd=sapw_dcdd_target)
    store_c_target, store_dcdd_target = bstore_allom(dbh, ipft, crown_damage, canopy_trim)

    repro_fraction = 0.0
    if mask_repro
        regen = ed_params().regeneration_model
        if regen == default_regeneration ||
           prt_params.allom_dbh_maxheight[ipft] < min_max_dbh_for_trees
            if dbh <= prt_params.dbh_repro_threshold[ipft]
                repro_fraction = prt_params.seed_alloc[ipft]
            else
                repro_fraction = prt_params.seed_alloc[ipft] + prt_params.seed_alloc_mature[ipft]
            end
        elseif (regen == TRS_regeneration || regen == TRS_no_seedling_dyn) &&
               prt_params.allom_dbh_maxheight[ipft] > min_max_dbh_for_trees
            ex = exp(prt_params.repro_alloc_b[ipft] + prt_params.repro_alloc_a[ipft] * dbh * mm_per_cm)
            repro_fraction = prt_params.seed_alloc[ipft] * (ex / (1.0 + ex))
        else
            fates_endrun("unknown seed allocation and regeneration model: $(regen)")
        end
    end

    total_dcostdd = 0.0
    mask_struct && (total_dcostdd += struct_dcdd_target)
    mask_leaf   && (total_dcostdd += leaf_dcdd_target)
    mask_fnrt   && (total_dcostdd += fnrt_dcdd_target)
    mask_sapw   && (total_dcostdd += sapw_dcdd_target)
    mask_store  && (total_dcostdd += store_dcdd_target)

    dCdx = zeros(Float64, length(l_state_array))

    # With asymptotic / hard-capped allometries all growth rates can reach zero.
    # In that case, give the carbon to reproduction.
    if total_dcostdd > nearzero
        mask_struct && (dCdx[struct_id] = struct_dcdd_target / total_dcostdd * (1.0 - repro_fraction))
        mask_leaf   && (dCdx[leaf_id]   = leaf_dcdd_target   / total_dcostdd * (1.0 - repro_fraction))
        mask_fnrt   && (dCdx[fnrt_id]   = fnrt_dcdd_target   / total_dcostdd * (1.0 - repro_fraction))
        mask_sapw   && (dCdx[sapw_id]   = sapw_dcdd_target   / total_dcostdd * (1.0 - repro_fraction))
        mask_store  && (dCdx[store_id]  = store_dcdd_target  / total_dcostdd * (1.0 - repro_fraction))
        mask_repro  && (dCdx[repro_id]  = repro_fraction)

        if abs(sum(dCdx) - 1.0) > rsnbl_math_prec
            fates_endrun("dCdx should sum to 1: dCdx=$(dCdx) repro_fraction=$(repro_fraction)")
        end

        dCdx[dbh_id] = (1.0 / total_dcostdd) * (1.0 - repro_fraction)
    else
        if repro_fraction < nearzero
            fates_endrun("A plant has reached a stature where allometry dictates no " *
                         "further growth in any non-reproductive pool, yet the " *
                         "reproductive fraction is also zero. dbh=$(dbh) pft=$(ipft)")
        end
        dCdx[repro_id] = 1.0
    end

    return dCdx
end

# =====================================================================================
# EstimateGrowthNC
# =====================================================================================

"""
    EstimateGrowthNC(this, target_c, target_dcdd, state_mask) -> (avg_nc, avg_pc)

Predict the effective N:C and P:C allocation ratios of the forthcoming stature
growth step (a dC/dd-weighted average of each masked organ's stoichiometry). Used
by the growth step to decide which element will be limiting.
"""
function EstimateGrowthNC(this::cnp_allom_prt_vartypes, target_c::AbstractVector,
                          target_dcdd::AbstractVector, state_mask::AbstractVector{Bool})
    dbh      = _rio(this, acnp_bc_inout_id_dbh)
    ipft     = _iin(this, acnp_bc_in_id_pft)
    nc_repro = _rin(this, acnp_bc_in_id_nc_repro)
    pc_repro = _rin(this, acnp_bc_in_id_pc_repro)

    repro_c_frac = 0.0
    if state_mask[repro_id]
        regen = ed_params().regeneration_model
        if regen == default_regeneration ||
           prt_params.allom_dbh_maxheight[ipft] < min_max_dbh_for_trees
            if dbh <= prt_params.dbh_repro_threshold[ipft]
                repro_c_frac = prt_params.seed_alloc[ipft]
            else
                repro_c_frac = prt_params.seed_alloc[ipft] + prt_params.seed_alloc_mature[ipft]
            end
        elseif (regen == TRS_regeneration || regen == TRS_no_seedling_dyn) &&
               prt_params.allom_dbh_maxheight[ipft] > min_max_dbh_for_trees
            ex = exp(prt_params.repro_alloc_b[ipft] + prt_params.repro_alloc_a[ipft] * dbh * mm_per_cm)
            repro_c_frac = prt_params.seed_alloc[ipft] * (ex / (1.0 + ex))
        else
            fates_endrun("unknown seed allocation and regeneration model: $(regen)")
        end
    end

    total_w = 0.0
    avg_nc  = 0.0
    avg_pc  = 0.0

    if state_mask[leaf_id]
        leaf_w = target_dcdd[leaf_organ] * (1.0 - repro_c_frac)
        total_w += leaf_w
        avg_nc += leaf_w * prt_params.nitr_stoich_p1[ipft, prt_params.organ_param_id[leaf_organ]]
        avg_pc += leaf_w * prt_params.phos_stoich_p1[ipft, prt_params.organ_param_id[leaf_organ]]
    end
    if state_mask[fnrt_id]
        fnrt_w = target_dcdd[fnrt_organ] * (1.0 - repro_c_frac)
        total_w += fnrt_w
        avg_nc += fnrt_w * prt_params.nitr_stoich_p1[ipft, prt_params.organ_param_id[fnrt_organ]]
        avg_pc += fnrt_w * prt_params.phos_stoich_p1[ipft, prt_params.organ_param_id[fnrt_organ]]
    end
    if state_mask[sapw_id]
        sapw_w = target_dcdd[sapw_organ] * (1.0 - repro_c_frac)
        total_w += sapw_w
        avg_nc += sapw_w * prt_params.nitr_stoich_p1[ipft, prt_params.organ_param_id[sapw_organ]]
        avg_pc += sapw_w * prt_params.phos_stoich_p1[ipft, prt_params.organ_param_id[sapw_organ]]
    end
    if state_mask[struct_id]
        struct_w = target_dcdd[struct_organ] * (1.0 - repro_c_frac)
        total_w += struct_w
        avg_nc += struct_w * prt_params.nitr_stoich_p1[ipft, prt_params.organ_param_id[struct_organ]]
        avg_pc += struct_w * prt_params.phos_stoich_p1[ipft, prt_params.organ_param_id[struct_organ]]
    end
    if state_mask[store_id]
        store_w = target_dcdd[store_organ] * (1.0 - repro_c_frac)
        total_w += store_w
        store_nc = GetNutrientTarget(this, nitrogen_element, store_organ, stoich_growth_min) / target_c[store_organ]
        store_pc = GetNutrientTarget(this, phosphorus_element, store_organ, stoich_growth_min) / target_c[store_organ]
        avg_nc += store_w * store_nc
        avg_pc += store_w * store_pc
    end
    if state_mask[repro_id]
        # repro_w = total_w * repro_c_frac/(1 - repro_c_frac)  (see Fortran derivation)
        repro_w = total_w * repro_c_frac / (1.0 - repro_c_frac)
        total_w += repro_w
        avg_nc += repro_w * nc_repro
        avg_pc += repro_w * pc_repro
    end

    avg_nc = avg_nc / total_w
    avg_pc = avg_pc / total_w

    return avg_nc, avg_pc
end

# =====================================================================================
# CNPAdjustFRootTargets — PID controller updating l2fr (and the fnrt target)
# =====================================================================================

"""
    CNPAdjustFRootTargets!(this, target_c, target_dcdd)

Update the leaf-to-fineroot target-biomass scalar (`l2fr`) via a PID controller
on the (log) ratio of relative carbon storage over relative nutrient storage,
then recompute the fine-root carbon target. Modifies `target_c[fnrt_organ]` and
`target_dcdd[fnrt_organ]` (and the l2fr / EMA boundary conditions) in place.
"""
function CNPAdjustFRootTargets!(this::cnp_allom_prt_vartypes,
                                target_c::AbstractVector, target_dcdd::AbstractVector)
    pid_drv_wgt = 1.0 / 20.0   # n-day smoothing of the PID derivative

    leaf_status = _iin(this, acnp_bc_in_id_lstat)
    ipft        = _iin(this, acnp_bc_in_id_pft)
    elongf_fnrt = _rin(this, acnp_bc_in_id_effnrt)
    l2fr        = _rio(this, acnp_bc_inout_id_l2fr)
    dbh         = _rio(this, acnp_bc_inout_id_dbh)
    canopy_trim = _rin(this, acnp_bc_in_id_ctrim)

    # Abort if leaves are off.
    leaf_status == leaves_off && return nothing

    npmode = ed_params().n_uptake_mode
    ppmode = ed_params().p_uptake_mode

    store_c_max = target_c[store_organ]
    store_c_act = max(0.001 * store_c_max,
                      GetState(this, store_organ, carbon12_element) +
                      _rin(this, acnp_bc_in_id_netdc))

    if npmode == prescribed_n_uptake
        cn_ratio = -1.0
    else
        store_nut_max = GetNutrientTarget(this, nitrogen_element, store_organ, stoich_growth_min)
        store_nut_act = max(0.001 * store_nut_max,
                            GetState(this, store_organ, nitrogen_element) +
                            _rio(this, acnp_bc_inout_id_netdn))
        cn_ratio = (store_c_act / store_c_max) / (store_nut_act / store_nut_max)
    end

    if ppmode == prescribed_p_uptake
        cp_ratio = -1.0
    else
        store_nut_max = GetNutrientTarget(this, phosphorus_element, store_organ, stoich_growth_min)
        store_nut_act = max(0.001 * store_nut_max,
                            GetState(this, store_organ, phosphorus_element) +
                            _rio(this, acnp_bc_inout_id_netdp))
        cp_ratio = (store_c_act / store_c_max) / (store_nut_act / store_nut_max)
    end

    if npmode == prescribed_n_uptake && ppmode == prescribed_p_uptake
        _set_rio!(this, acnp_bc_inout_id_cx_int, 0.0)
        _set_rio!(this, acnp_bc_inout_id_emadcxdt, 0.0)
        _set_rio!(this, acnp_bc_inout_id_cx0, 0.0)
        return nothing
    end

    if npmode == prescribed_n_uptake
        cx_logratio = SafeLog(cp_ratio)
    elseif ppmode == prescribed_p_uptake
        cx_logratio = SafeLog(cn_ratio)
    else
        cx_logratio = SafeLog(max(cp_ratio, cn_ratio))
    end

    cx_int = _rio(this, acnp_bc_inout_id_cx_int)
    cx0    = _rio(this, acnp_bc_inout_id_cx0)
    ema_dcxdt = _rio(this, acnp_bc_inout_id_emadcxdt)

    cx_int = cx_int + cx_logratio

    # Reset the integrator if the sign of the ratio changed.
    if abs(cx_logratio) > nearzero && abs(cx0) > nearzero
        if abs(cx_logratio / abs(cx_logratio) - cx0 / abs(cx0)) > nearzero
            cx_int = cx_logratio
        end
    end

    dcxdt_ratio = cx_logratio - cx0
    ema_dcxdt = pid_drv_wgt * dcxdt_ratio + (1.0 - pid_drv_wgt) * ema_dcxdt
    cx0 = cx_logratio

    _set_rio!(this, acnp_bc_inout_id_cx_int, cx_int)
    _set_rio!(this, acnp_bc_inout_id_emadcxdt, ema_dcxdt)
    _set_rio!(this, acnp_bc_inout_id_cx0, cx0)

    l2fr_delta = prt_params.pid_kp[ipft] * cx_logratio +
                 prt_params.pid_ki[ipft] * cx_int +
                 prt_params.pid_kd[ipft] * ema_dcxdt

    # Apply the delta, avoiding incredibly small l2fr's.
    l2fr = max(l2fr_min, l2fr + l2fr_delta)
    _set_rio!(this, acnp_bc_inout_id_l2fr, l2fr)

    # Find the updated target fineroot biomass.
    target_c[fnrt_organ], target_dcdd[fnrt_organ] =
        bfineroot(dbh, ipft, canopy_trim, l2fr, elongf_fnrt)

    return nothing
end

# =====================================================================================
# TrimFineRoot — forceful turnover when a lower l2fr is generated
# =====================================================================================

"""
    TrimFineRoot!(this)

If a new (smaller) l2fr leaves fine-roots above their target, force the excess
biomass into turnover (→ litter) immediately rather than waiting for background
turnover to catch up. (With `fnrt_opt_eff=0` no mass goes to storage.)
"""
function TrimFineRoot!(this::cnp_allom_prt_vartypes)
    nday_buffer = 0.0
    fnrt_opt_eff = 0.0   # transfer fraction to storage (0 -> all to turnover)

    use_unrestricted_contraction || return nothing

    ipft        = _iin(this, acnp_bc_in_id_pft)
    l2fr        = _rio(this, acnp_bc_inout_id_l2fr)
    dbh         = _rio(this, acnp_bc_inout_id_dbh)
    elongf_fnrt = _rin(this, acnp_bc_in_id_effnrt)
    canopy_trim = _rin(this, acnp_bc_in_id_ctrim)

    target_fnrt_c, _ = bfineroot(dbh, ipft, canopy_trim, l2fr, elongf_fnrt)

    fnrt_flux_c = max(0.0,
        this.variables[fnrt_c_id].val[1] *
        (1.0 - nday_buffer * (years_per_day / prt_params.root_long[ipft])) - target_fnrt_c)

    if fnrt_flux_c > nearzero
        turn_flux_c  = (1.0 - fnrt_opt_eff) * fnrt_flux_c
        store_flux_c = fnrt_opt_eff * fnrt_flux_c

        nc_fnrt = this.variables[fnrt_n_id].val[1] / this.variables[fnrt_c_id].val[1]
        pc_fnrt = this.variables[fnrt_p_id].val[1] / this.variables[fnrt_c_id].val[1]

        this.variables[fnrt_c_id].val[1]        -= fnrt_flux_c
        this.variables[fnrt_c_id].turnover[1]   += turn_flux_c
        this.variables[fnrt_c_id].net_alloc[1]  -= store_flux_c
        this.variables[store_c_id].val[1]       += store_flux_c
        this.variables[store_c_id].net_alloc[1] += store_flux_c

        this.variables[fnrt_n_id].val[1]        -= fnrt_flux_c * nc_fnrt
        this.variables[fnrt_n_id].turnover[1]   += turn_flux_c * nc_fnrt
        this.variables[fnrt_n_id].net_alloc[1]  -= store_flux_c * nc_fnrt
        this.variables[store_n_id].val[1]       += store_flux_c * nc_fnrt
        this.variables[store_n_id].net_alloc[1] += store_flux_c * nc_fnrt

        this.variables[fnrt_p_id].val[1]        -= fnrt_flux_c * pc_fnrt
        this.variables[fnrt_p_id].turnover[1]   += turn_flux_c * pc_fnrt
        this.variables[fnrt_p_id].net_alloc[1]  -= store_flux_c * pc_fnrt
        this.variables[store_p_id].val[1]       += store_flux_c * pc_fnrt
        this.variables[store_p_id].net_alloc[1] += store_flux_c * pc_fnrt
    end

    return nothing
end

# =====================================================================================
# CNPPrioritizedReplacement — Step 2 of the daily allocation
# =====================================================================================

"""
    CNPPrioritizedReplacement!(this, c_gain, n_gain, p_gain, target_c) -> (c_gain, n_gain, p_gain)

Prioritized replacement of tissue turnover (and payment of unpaid maintenance
respiration from storage), then bring pools up to their carbon allometric targets
in priority order, filling nutrients proportionally to deficit. Returns the
updated gain pools (caller reassigns). `target_c` is indexed by global organ id.
"""
function CNPPrioritizedReplacement!(this::cnp_allom_prt_vartypes,
                                    c_gain::Float64, n_gain::Float64, p_gain::Float64,
                                    target_c::AbstractVector)
    g = get_prt_global()

    leaf_status = _iin(this, acnp_bc_in_id_lstat)
    ipft        = _iin(this, acnp_bc_in_id_pft)

    curpri_org = fill(fates_unset_int, acnp_num_organs)
    deficit_c  = zeros(Float64, acnp_num_organs)
    deficit_n  = zeros(Float64, acnp_num_organs)
    deficit_p  = zeros(Float64, acnp_num_organs)

    n_max_priority = maximum(prt_params.organ_param_id)
    if n_max_priority > 10 || n_max_priority < 0
        fates_endrun("Unable to interpret prt_params.organ_param_id for cnp allocation; " *
                     "values should be non-zero and <10: $(prt_params.organ_param_id)")
    end

    # ---- Highest-priority pools (priority_code == 1) ---------------------------------
    i = 0
    for ii in 1:length(prt_params.organ_id)
        i_org = prt_params.organ_id[ii]

        # Don't allocate to leaves if "off"/"shedding", or to replace turnover if
        # not evergreen (prevents accidental re-flushing on drop day).
        if ((leaf_status == leaves_off || leaf_status == leaves_shedding) ||
            (prt_params.evergreen[ipft] != itrue)) && (i_org == leaf_organ)
            continue
        end

        priority_code = Int(prt_params.alloc_priority[ipft, ii])
        if priority_code == 1
            i += 1
            curpri_org[i] = i_org
            deficit_c[i]  = max(0.0, GetDeficit(this, carbon12_element, i_org, target_c[i_org]))
        end
    end
    n_curpri_org = i

    # ---- Replace maintenance turnover of the priority-1 pools ------------------------
    sum_c_demand = 0.0
    for i in 1:n_curpri_org
        i_org = curpri_org[i]
        i_var = sp_organ_map_get(g, i_org, carbon12_element)
        sum_c_demand += prt_params.leaf_stor_priority[ipft] * sum(this.variables[i_var].turnover)
    end

    sum_c_flux = max(0.0, min(sum_c_demand, this.variables[store_c_id].val[1] + c_gain))

    if sum_c_flux > nearzero
        for i in 1:n_curpri_org
            i_org = curpri_org[i]
            i_var = sp_organ_map_get(g, i_org, carbon12_element)
            c_flux = sum_c_flux * (prt_params.leaf_stor_priority[ipft] *
                                   sum(this.variables[i_var].turnover) / sum_c_demand)
            this.variables[i_var].val[1] += c_flux
            c_gain -= c_flux
        end
    end

    # ---- Nutrient demand for the priority-1 pools ------------------------------------
    for i in 1:n_curpri_org
        i_org = curpri_org[i]
        target_n = GetNutrientTarget(this, nitrogen_element, i_org, stoich_growth_min)
        deficit_n[i] = max(0.0, target_n - GetState(this, i_org, nitrogen_element, 1))
        target_p = GetNutrientTarget(this, phosphorus_element, i_org, stoich_growth_min)
        deficit_p[i] = max(0.0, target_p - GetState(this, i_org, phosphorus_element, 1))
    end

    n_gain = ProportionalNutrAllocation!(this, view(deficit_n, 1:n_curpri_org), n_gain,
                                         nitrogen_element, view(curpri_org, 1:n_curpri_org))
    p_gain = ProportionalNutrAllocation!(this, view(deficit_p, 1:n_curpri_org), p_gain,
                                         phosphorus_element, view(curpri_org, 1:n_curpri_org))

    # ---- IV. Reconcile the carbon balance with storage ------------------------------
    if c_gain < 0.0
        store_c_flux = -c_gain
        c_gain       = c_gain + store_c_flux
        this.variables[store_c_id].val[1] -= store_c_flux
    else
        store_below_target    = max(target_c[store_organ] - this.variables[store_c_id].val[1], 0.0)
        store_target_fraction = max(this.variables[store_c_id].val[1] / target_c[store_organ], 0.0)
        store_demand          = max(c_gain * (exp(-1.0 * store_target_fraction^4.0) - exp(-1.0)), 0.0)
        store_c_flux          = min(store_below_target, store_demand)
        c_gain                = c_gain - store_c_flux
        this.variables[store_c_id].val[1] += store_c_flux
    end

    # ---- Bring all pools up to allometric C targets in priority order ----------------
    for i_pri in 1:n_max_priority
        fill!(curpri_org, fates_unset_int)
        i = 0

        # Storage has a special hard-coded priority level of 2.
        if i_pri == 2
            curpri_org[1] = store_organ
            i = 1
        end

        for ii in 1:length(prt_params.organ_id)
            i_org = prt_params.organ_id[ii]
            priority_code = Int(prt_params.alloc_priority[ipft, ii])

            if (leaf_status == leaves_off || leaf_status == leaves_shedding) &&
               (i_org == leaf_organ)
                continue
            end

            if priority_code == i_pri
                i += 1
                curpri_org[i] = i_org
            end
        end
        n_curpri_org = i

        for i in 1:n_curpri_org
            i_org = curpri_org[i]
            deficit_c[i] = max(0.0, GetDeficit(this, carbon12_element, i_org, target_c[i_org]))
        end

        # Carbon up to target first (sets the carbon concentrations the nutrient
        # targets depend on).
        sum_c_demand = 0.0
        for i in 1:n_curpri_org
            sum_c_demand += deficit_c[i]
        end

        sum_c_flux = min(c_gain, sum_c_demand)
        if sum_c_flux > nearzero
            for i in 1:n_curpri_org
                i_org = curpri_org[i]
                c_flux = sum_c_flux * deficit_c[i] / sum_c_demand
                i_var = sp_organ_map_get(g, i_org, carbon12_element)
                this.variables[i_var].val[1] += c_flux
                deficit_c[i] = max(0.0, deficit_c[i] - c_flux)
                c_gain -= c_flux
            end
        end

        # Nutrient demand at this priority level.
        for i in 1:n_curpri_org
            i_org = curpri_org[i]
            target_n = GetNutrientTarget(this, nitrogen_element, i_org, stoich_growth_min)
            deficit_n[i] = max(0.0, target_n - GetState(this, i_org, nitrogen_element, 1))
            target_p = GetNutrientTarget(this, phosphorus_element, i_org, stoich_growth_min)
            deficit_p[i] = max(0.0, target_p - GetState(this, i_org, phosphorus_element, 1))
        end

        n_gain = ProportionalNutrAllocation!(this, view(deficit_n, 1:n_curpri_org), n_gain,
                                             nitrogen_element, view(curpri_org, 1:n_curpri_org))
        p_gain = ProportionalNutrAllocation!(this, view(deficit_p, 1:n_curpri_org), p_gain,
                                             phosphorus_element, view(curpri_org, 1:n_curpri_org))
    end

    return c_gain, n_gain, p_gain
end

# =====================================================================================
# CNPStatureGrowth — Step 3 of the daily allocation
# =====================================================================================

"""
    CNPStatureGrowth!(this, c_gain, n_gain, p_gain, target_c, target_dcdd) -> (c_gain, n_gain, p_gain)

Grow the plant's stature (allocate beyond current allometric targets), advancing
dbh through an Euler integration of `AllomCNPGrowthDeriv`, then proportionally
fill the resulting nutrient demand. The amount of carbon committed to growth is
limited by the equivalent-carbon estimate of N and P availability. Returns the
updated gain pools (caller reassigns).
"""
function CNPStatureGrowth!(this::cnp_allom_prt_vartypes,
                           c_gain::Float64, n_gain::Float64, p_gain::Float64,
                           target_c::AbstractVector, target_dcdd::AbstractVector)
    g = get_prt_global()

    max_substeps    = 300
    max_trunc_error = 1.0
    ODESolve        = 2   # 1=RKF45, 2=Euler

    leaf_status = Float64(_iin(this, acnp_bc_in_id_lstat))
    elongf_leaf = _rin(this, acnp_bc_in_id_efleaf)
    elongf_fnrt = _rin(this, acnp_bc_in_id_effnrt)
    elongf_stem = _rin(this, acnp_bc_in_id_efstem)
    dbh         = _rio(this, acnp_bc_inout_id_dbh)
    ipft        = _iin(this, acnp_bc_in_id_pft)
    crown_damage = _iin(this, acnp_bc_in_id_cdamage)
    canopy_trim = _rin(this, acnp_bc_in_id_ctrim)
    l2fr        = _rio(this, acnp_bc_inout_id_l2fr)

    limiter = 0
    if c_gain <= calloc_abs_error
        limiter = c_limited
        if (n_gain <= 0.1 * calloc_abs_error) || (p_gain <= 0.02 * calloc_abs_error)
            limiter = cnp_limited
        end
    else
        n_gain <= 0.1 * calloc_abs_error && (limiter = n_limited)
        p_gain <= 0.02 * calloc_abs_error && (limiter = p_limited)
    end
    limiter = 0
    _set_iout!(this, acnp_bc_out_id_limiter, limiter)

    # No point growing if any resource is tapped out, or leaves are off.
    if c_gain <= calloc_abs_error ||
       n_gain <= 0.1 * calloc_abs_error ||
       p_gain <= 0.02 * calloc_abs_error ||
       (leaf_status == leaves_off || leaf_status == leaves_shedding)
        return c_gain, n_gain, p_gain
    end

    intgr_params = fill(fates_unset_r8, num_intgr_parm)
    intgr_params[intgr_parm_ctrim]   = canopy_trim
    intgr_params[intgr_parm_pft]     = Float64(ipft)
    intgr_params[intgr_parm_l2fr]    = l2fr
    intgr_params[intgr_parm_cdamage] = Float64(crown_damage)
    intgr_params[intgr_parm_efleaf]  = elongf_leaf
    intgr_params[intgr_parm_effnrt]  = elongf_fnrt
    intgr_params[intgr_parm_efstem]  = elongf_stem

    state_mask  = fill(false, num_intgr_vars)
    mask_organs  = fill(fates_unset_int, acnp_num_organs)
    mask_gorgans = fill(fates_unset_int, acnp_num_organs)

    # Flag pools as growing (on allometry) or not (above target -> skip this step).
    ii = 0
    for i in 1:acnp_num_organs
        i_org = l2g_organ_list[i]
        cdeficit = GetDeficit(this, carbon12_element, i_org, target_c[i_org])

        if cdeficit > calloc_abs_error
            fates_endrun("A carbon pool reached the stature-growth step yet its " *
                         "deficit is too large to integrate. organ=$(i_org) " *
                         "c_gain=$(c_gain) leaf_status=$(leaf_status) " *
                         "cdeficit=$(cdeficit) target=$(target_c[i_org])")
        elseif (-cdeficit) > calloc_abs_error
            state_mask[i] = false   # above target -> will catch up next step
        else
            state_mask[i] = true
            if i_org != repro_organ
                ii += 1
                mask_organs[ii]  = i
                mask_gorgans[ii] = i_org
            end
        end
    end
    n_mask_organs = ii

    # Reproductive fraction of carbon (special-cased, not proportional).
    regen = ed_params().regeneration_model
    local repro_c_frac::Float64
    if regen == default_regeneration ||
       prt_params.allom_dbh_maxheight[ipft] < min_max_dbh_for_trees
        if dbh <= prt_params.dbh_repro_threshold[ipft]
            repro_c_frac = prt_params.seed_alloc[ipft]
        else
            repro_c_frac = prt_params.seed_alloc[ipft] + prt_params.seed_alloc_mature[ipft]
        end
    elseif (regen == TRS_regeneration || regen == TRS_no_seedling_dyn) &&
           prt_params.allom_dbh_maxheight[ipft] > min_max_dbh_for_trees
        ex = exp(prt_params.repro_alloc_b[ipft] + prt_params.repro_alloc_a[ipft] * dbh * mm_per_cm)
        repro_c_frac = prt_params.seed_alloc[ipft] * (ex / (1.0 + ex))
    else
        fates_endrun("unknown seed allocation and regeneration model: $(regen)")
    end

    if repro_c_frac > nearzero
        state_mask[repro_id] = true
        n_mask_organs += 1
        mask_organs[n_mask_organs]  = repro_id
        mask_gorgans[n_mask_organs] = repro_organ
    else
        state_mask[repro_id] = false
    end

    # NOTE (verbatim bug preservation): the Fortran computes total_dcostdd here
    # using `i_org = mask_gorgans(ii)` where `ii` is left over from the masking
    # loop (NOT the loop index `i`). This makes the value effectively unused
    # downstream (the integrator recomputes total_dcostdd internally), so we keep
    # the harmless computation but do not rely on it.
    total_dcostdd = 0.0
    for i in 1:n_mask_organs
        i_org = mask_gorgans[ii]
        i_org == fates_unset_int && continue
        total_dcostdd += target_dcdd[i_org]
    end

    # Decide how much carbon to commit to growth (C-only or NP-equivalent-limited).
    local c_gstature::Float64
    if grow_lim_type == grow_lim_conly
        c_gstature = c_gain
        limiter = 0
    else  # grow_lim_estNP
        avg_nc, avg_pc = EstimateGrowthNC(this, target_c, target_dcdd, state_mask)
        neq_cgain = n_gain / avg_nc
        peq_cgain = p_gain / avg_pc

        if c_gain < neq_cgain
            if c_gain < peq_cgain
                limiter = c_limited
                c_gstature = c_gain
            else
                limiter = p_limited
                c_gstature = peq_cgain
            end
        else
            if neq_cgain < peq_cgain
                limiter = n_limited
                c_gstature = neq_cgain
            else
                limiter = p_limited
                c_gstature = peq_cgain
            end
        end
    end
    _set_iout!(this, acnp_bc_out_id_limiter, limiter)

    if c_gstature > nearzero

        if ODESolve == 2
            this.ode_opt_step = c_gstature
        end

        ierr   = 1
        nsteps = 0
        totalC = c_gstature

        state_array     = zeros(Float64, num_intgr_vars)
        state_array_out = zeros(Float64, num_intgr_vars)

        for i in 1:acnp_num_organs
            i_org = l2g_organ_list[i]
            i_var = sp_organ_map_get(g, i_org, carbon12_element)
            state_array[i] = this.variables[i_var].val[1]
        end
        state_mask[dbh_id]  = true
        state_array[dbh_id] = dbh

        step_pass = false
        while ierr != 0
            deltaC = min(totalC, this.ode_opt_step)

            if ODESolve == 1
                this.ode_opt_step, step_pass =
                    RKF45(AllomCNPGrowthDeriv, state_array, state_mask, deltaC, totalC,
                          max_trunc_error, intgr_params, state_array_out)
            elseif ODESolve == 2
                Euler(AllomCNPGrowthDeriv, state_array, state_mask, deltaC, totalC,
                      intgr_params, state_array_out)

                # Check the solution is reasonably close to allometry; sum leaf bins.
                leafc_tp1 = state_array_out[leaf_id]
                i_var = sp_organ_map_get(g, leaf_organ, carbon12_element)
                nbins = g.state_descriptor[i_var].num_pos
                for i in 2:nbins
                    leafc_tp1 += this.variables[i_var].val[i]
                end

                step_pass = CheckIntegratedAllometries(
                    state_array_out[dbh_id], ipft, crown_damage, canopy_trim,
                    elongf_leaf, elongf_fnrt, elongf_stem, l2fr,
                    leafc_tp1, state_array_out[fnrt_id], state_array_out[sapw_id],
                    state_array_out[store_id], state_array_out[struct_id],
                    state_mask[leaf_id], state_mask[fnrt_id], state_mask[sapw_id],
                    state_mask[store_id], state_mask[struct_id], max_trunc_error)

                this.ode_opt_step = step_pass ? deltaC : 0.5 * deltaC
            else
                fates_endrun("An integrator was chosen that DNE: ODESolve=$(ODESolve)")
            end

            nsteps += 1

            if step_pass
                totalC = totalC - deltaC
                state_array .= state_array_out
            end

            if totalC < calloc_abs_error
                ierr = 0

                # Sum the total flux the integrator predicts, then proportionally
                # correct all pools so carbon is conserved exactly.
                sum_c_flux = 0.0
                for jj in 1:n_mask_organs
                    i     = mask_organs[jj]
                    i_org = mask_gorgans[jj]
                    i_var = sp_organ_map_get(g, i_org, carbon12_element)
                    sum_c_flux += (state_array[i] - this.variables[i_var].val[1])
                end

                c_flux_adj = c_gstature / sum_c_flux

                for jj in 1:n_mask_organs
                    i     = mask_organs[jj]
                    i_org = mask_gorgans[jj]
                    i_var = sp_organ_map_get(g, i_org, carbon12_element)
                    c_flux = (state_array[i] - this.variables[i_var].val[1]) * c_flux_adj
                    this.variables[i_var].val[1] += c_flux
                    c_gain -= c_flux
                end

                dbh = state_array[dbh_id]
                _set_rio!(this, acnp_bc_inout_id_dbh, dbh)
            else
                if nsteps > max_substeps
                    fates_endrun("CNP Plant Growth Integrator could not find a " *
                                 "solution in < $(max_substeps) tries. pft=$(ipft) " *
                                 "dbh=$(dbh) totalC=$(totalC) smallest deltaC=$(this.ode_opt_step)")
                end
            end
        end

        # Prioritize balanced-CNP nutrient transfer to the reproductive pool.
        if prioritize_repro_nutr_growth
            target_n = GetNutrientTarget(this, nitrogen_element, repro_organ, stoich_growth_min)
            dn = GetDeficit(this, nitrogen_element, repro_organ, target_n)
            n_flux = max(0.0, min(n_gain, dn))

            target_p = GetNutrientTarget(this, phosphorus_element, repro_organ, stoich_growth_min)
            dp = GetDeficit(this, phosphorus_element, repro_organ, target_p)
            p_flux = max(0.0, min(p_gain, dp))

            this.variables[repro_n_id].val[1] += n_flux
            this.variables[repro_p_id].val[1] += p_flux
            n_gain -= n_flux
            p_gain -= p_flux
        end

        # Proportional nutrient fluxes to each masked pool.
        deficit_n = zeros(Float64, acnp_num_organs)
        deficit_p = zeros(Float64, acnp_num_organs)
        for jj in 1:n_mask_organs
            i_org = mask_gorgans[jj]
            target_n = GetNutrientTarget(this, nitrogen_element, i_org, stoich_growth_min)
            target_p = GetNutrientTarget(this, phosphorus_element, i_org, stoich_growth_min)
            deficit_n[jj] = GetDeficit(this, nitrogen_element, i_org, target_n)
            deficit_p[jj] = GetDeficit(this, phosphorus_element, i_org, target_p)
        end

        n_gain = ProportionalNutrAllocation!(this, view(deficit_n, 1:n_mask_organs), n_gain,
                                             nitrogen_element, view(mask_gorgans, 1:n_mask_organs))
        p_gain = ProportionalNutrAllocation!(this, view(deficit_p, 1:n_mask_organs), p_gain,
                                             phosphorus_element, view(mask_gorgans, 1:n_mask_organs))
    end

    return c_gain, n_gain, p_gain
end

# =====================================================================================
# CNPAllocateRemainder — Step 4: allocate or efflux remaining resources
# =====================================================================================

"""
    CNPAllocateRemainder!(this, c_gain, n_gain, p_gain, target_c, target_dcdd)
        -> (c_gain, n_gain, p_gain, c_efflux, n_efflux, p_efflux)

After growth, bring nutrient pools toward their optimal targets, update the l2fr
controller, cram excess carbon into storage overflow (or burn/exude per
`store_c_overflow`), then efflux whatever cannot be allocated. Returns the
(zeroed) gains plus the three efflux quantities.
"""
function CNPAllocateRemainder!(this::cnp_allom_prt_vartypes,
                               c_gain::Float64, n_gain::Float64, p_gain::Float64,
                               target_c::AbstractVector, target_dcdd::AbstractVector)
    dbh         = _rio(this, acnp_bc_inout_id_dbh)
    canopy_trim = _rin(this, acnp_bc_in_id_ctrim)
    ipft        = _iin(this, acnp_bc_in_id_pft)
    crown_damage = _iin(this, acnp_bc_in_id_cdamage)
    resp_excess = _rio(this, acnp_bc_inout_id_resp_excess)

    deficit_n = zeros(Float64, acnp_num_organs)
    deficit_p = zeros(Float64, acnp_num_organs)

    # Bump nutrient pools toward optimal (allowing storage overflow).
    for i in 1:acnp_num_organs
        i_org = l2g_organ_list[i]
        target_n = GetNutrientTarget(this, nitrogen_element, i_org, stoich_growth_min)
        target_p = GetNutrientTarget(this, phosphorus_element, i_org, stoich_growth_min)

        if i_org == store_organ
            target_n = target_n * (1.0 + prt_params.store_ovrflw_frac[ipft])
            target_p = target_p * (1.0 + prt_params.store_ovrflw_frac[ipft])
        end

        deficit_n[i] = max(0.0, GetDeficit(this, nitrogen_element, i_org, target_n))
        deficit_p[i] = max(0.0, GetDeficit(this, phosphorus_element, i_org, target_p))
    end

    n_gain = ProportionalNutrAllocation!(this, view(deficit_n, 1:acnp_num_organs), n_gain,
                                         nitrogen_element, l2g_organ_list)
    p_gain = ProportionalNutrAllocation!(this, view(deficit_p, 1:acnp_num_organs), p_gain,
                                         phosphorus_element, l2g_organ_list)

    # Update the l2fr controller and fineroot target.
    CNPAdjustFRootTargets!(this, target_c, target_dcdd)

    # Cram remaining carbon into storage overflow (or burn/exude).
    if c_gain > calloc_abs_error
        if store_c_overflow == retain_c_store_overflow
            total_c_flux = c_gain
            this.variables[store_c_id].val[1] += total_c_flux
            c_gain -= total_c_flux
        elseif store_c_overflow == burn_c_store_overflow
            store_c_target, _ = bstore_allom(dbh, ipft, crown_damage, canopy_trim)
            store_c_target = store_c_target * (1.0 + prt_params.store_ovrflw_frac[ipft])
            total_c_flux = max(0.0, min(c_gain, store_c_target - this.variables[store_c_id].val[1]))
            this.variables[store_c_id].val[1] += total_c_flux
            c_gain -= total_c_flux

            resp_excess = resp_excess + c_gain
            c_gain = 0.0
        elseif store_c_overflow == exude_c_store_overflow
            store_c_target, _ = bstore_allom(dbh, ipft, crown_damage, canopy_trim)
            store_c_target = store_c_target * (1.0 + prt_params.store_ovrflw_frac[ipft])
            total_c_flux = max(0.0, min(c_gain, store_c_target - this.variables[store_c_id].val[1]))
            this.variables[store_c_id].val[1] += total_c_flux
            c_gain -= total_c_flux
        end
    end

    # Recover small negative gains from storage (numerical precision).
    if c_gain < -nearzero
        this.variables[store_c_id].val[1] += c_gain
        c_gain = 0.0
    end
    if n_gain < -nearzero
        this.variables[store_n_id].val[1] += n_gain
        n_gain = 0.0
    end
    if p_gain < -nearzero
        this.variables[store_p_id].val[1] += p_gain
        p_gain = 0.0
    end

    # Efflux the remainder (unless prescribed mode, in which case keep the gain to
    # report demand as the amount uptaken).
    npmode = ed_params().n_uptake_mode
    ppmode = ed_params().p_uptake_mode

    if npmode == prescribed_n_uptake
        n_efflux = 0.0
    else
        n_efflux = n_gain
        n_gain   = 0.0
    end

    if ppmode == prescribed_p_uptake
        p_efflux = 0.0
    else
        p_efflux = p_gain
        p_gain   = 0.0
    end

    c_efflux = c_gain
    c_gain   = 0.0

    _set_rio!(this, acnp_bc_inout_id_resp_excess, resp_excess)

    return c_gain, n_gain, p_gain, c_efflux, n_efflux, p_efflux
end

# =====================================================================================
# FastPRT (deferred-generic override) — stub (no sub-daily processes)
# =====================================================================================

"""
    FastPRT!(this::cnp_allom_prt_vartypes)

Fast (sub-daily) reactive transport. A stub for this hypothesis — there are no
fast-timestep processes.
"""
function FastPRT!(::cnp_allom_prt_vartypes)
    return nothing
end

# =====================================================================================
# DailyPRT (the deferred-generic override) — the main daily allocation driver
# =====================================================================================

"""
    DailyPRT!(this::cnp_allom_prt_vartypes, phase)

The CNP coupled C/N/P daily allocation. Phasing is only used to accommodate the
damage module, which is incompatible with CNP, so any phase but 1 returns
immediately. The sequence is:

  0. Transfer storage overflow (C above target, all stored N/P) into the daily
     gain pools.
  1/2. Prioritized replacement of turnover + bring pools to C allometric targets,
     filling nutrients proportionally ([`CNPPrioritizedReplacement!`](@ref)).
  3. Stature growth ([`CNPStatureGrowth!`](@ref)).
  4. Allocate / efflux the remainder ([`CNPAllocateRemainder!`](@ref)).

Mass conservation is checked after steps 2 and 3 (and, under `acnp_debug`, a final
full tally). Updates the per-organ daily `net_alloc` diagnostics and the
boundary-condition efflux/gain outputs, and forcefully trims fine-roots if a
smaller l2fr was generated ([`TrimFineRoot!`](@ref)).
"""
function DailyPRT!(this::cnp_allom_prt_vartypes, phase::Integer)
    g = get_prt_global()

    # Phasing only accommodates the damage module (incompatible with CNP).
    phase != 1 && return nothing

    npmode = ed_params().n_uptake_mode
    ppmode = ed_params().p_uptake_mode

    # In/out boundary conditions.
    dbh0   = _rio(this, acnp_bc_inout_id_dbh)
    l2fr   = _rio(this, acnp_bc_inout_id_l2fr)
    n_gain = _rio(this, acnp_bc_inout_id_netdn)
    p_gain = _rio(this, acnp_bc_inout_id_netdp)

    # Zero excess respiration (it is updated in step 4 if needed).
    resp_excess0 = 0.0
    _set_rio!(this, acnp_bc_inout_id_resp_excess, 0.0)

    c_gain      = _rin(this, acnp_bc_in_id_netdc)
    canopy_trim = _rin(this, acnp_bc_in_id_ctrim)
    ipft        = _iin(this, acnp_bc_in_id_pft)
    crown_damage = _iin(this, acnp_bc_in_id_cdamage)
    elongf_leaf  = _rin(this, acnp_bc_in_id_efleaf)
    elongf_fnrt  = _rin(this, acnp_bc_in_id_effnrt)
    elongf_stem  = _rin(this, acnp_bc_in_id_efstem)

    # Prescribed mode -> set gains very large (1 kg of pure nutrient).
    if npmode == prescribed_n_uptake
        n_gain = 1.0e3
    end
    if ppmode == prescribed_p_uptake
        p_gain = 1.0e3
    end

    n_gain0 = n_gain
    p_gain0 = p_gain
    c_gain0 = c_gain

    dbh = dbh0

    # ---- Carbon allometry targets at current stature --------------------------------
    target_c    = fill(fates_unset_r8, num_organ_types)
    target_dcdd = fill(fates_unset_r8, num_organ_types)

    sapw_area, target_c[sapw_organ], target_dcdd[sapw_organ] =
        bsap_allom(dbh, ipft, crown_damage, canopy_trim, elongf_stem)
    agw_c_target, agw_dcdd_target = bagw_allom(dbh, ipft, crown_damage, elongf_stem)
    bgw_c_target, bgw_dcdd_target = bbgw_allom(dbh, ipft, elongf_stem)
    target_c[struct_organ], target_dcdd[struct_organ] =
        bdead_allom(agw_c_target, bgw_c_target, target_c[sapw_organ], ipft;
                    dbagwdd=agw_dcdd_target, dbbgwdd=bgw_dcdd_target,
                    dbsapdd=target_dcdd[sapw_organ])
    target_c[leaf_organ], target_dcdd[leaf_organ] =
        bleaf(dbh, ipft, crown_damage, canopy_trim, elongf_leaf)
    target_c[fnrt_organ], target_dcdd[fnrt_organ] =
        bfineroot(dbh, ipft, canopy_trim, l2fr, elongf_fnrt)
    target_c[store_organ], target_dcdd[store_organ] =
        bstore_allom(dbh, ipft, crown_damage, canopy_trim)
    target_c[repro_organ]    = 0.0
    target_dcdd[repro_organ] = 0.0

    # ---- Remember original states for the final allocation tally ---------------------
    state_c0 = zeros(Float64, num_organ_types)
    state_n0 = zeros(Float64, num_organ_types)
    state_p0 = zeros(Float64, num_organ_types)
    for i in 1:acnp_num_organs
        i_org = l2g_organ_list[i]
        state_c0[i_org] = this.variables[sp_organ_map_get(g, i_org, carbon12_element)].val[1]
        state_n0[i_org] = this.variables[sp_organ_map_get(g, i_org, nitrogen_element)].val[1]
        state_p0[i_org] = this.variables[sp_organ_map_get(g, i_org, phosphorus_element)].val[1]
    end

    # Output-only efflux BCs.
    _set_rout!(this, acnp_bc_out_id_cefflux, 0.0)
    _set_rout!(this, acnp_bc_out_id_nefflux, 0.0)
    _set_rout!(this, acnp_bc_out_id_pefflux, 0.0)

    # ---- Step 0: transfer storage overflow into the daily gain pools -----------------
    store_flux = max(0.0, this.variables[store_c_id].val[1] - target_c[store_organ])
    c_gain += store_flux
    this.variables[store_c_id].val[1] -= store_flux

    n_gain += sum(this.variables[store_n_id].val)
    fill!(this.variables[store_n_id].val, 0.0)

    p_gain += sum(this.variables[store_p_id].val)
    fill!(this.variables[store_p_id].val, 0.0)

    # ---- Step 2: prioritized replacement --------------------------------------------
    c_gain, n_gain, p_gain = CNPPrioritizedReplacement!(this, c_gain, n_gain, p_gain, target_c)

    # Carbon balance check I.
    sum_c = 0.0
    for i in 1:acnp_num_organs
        i_org = l2g_organ_list[i]
        sum_c += this.variables[sp_organ_map_get(g, i_org, carbon12_element)].val[1]
    end
    if abs((c_gain0 - c_gain) - (sum_c - sum(state_c0))) > calloc_abs_error
        fates_endrun("Carbon not balancing I in CNP DailyPRT")
    end

    # ---- Step 3: stature growth -----------------------------------------------------
    c_gain, n_gain, p_gain = CNPStatureGrowth!(this, c_gain, n_gain, p_gain, target_c, target_dcdd)

    # Carbon balance check II.
    sum_c = 0.0
    for i in 1:acnp_num_organs
        i_org = l2g_organ_list[i]
        sum_c += this.variables[sp_organ_map_get(g, i_org, carbon12_element)].val[1]
    end
    if abs((c_gain0 - c_gain) - (sum_c - sum(state_c0))) > calloc_abs_error
        fates_endrun("Carbon not balancing II in CNP DailyPRT")
    end

    # ---- Step 4: allocate remainder / efflux ----------------------------------------
    c_gain, n_gain, p_gain, c_efflux, n_efflux, p_efflux =
        CNPAllocateRemainder!(this, c_gain, n_gain, p_gain, target_c, target_dcdd)
    _set_rout!(this, acnp_bc_out_id_cefflux, c_efflux)
    _set_rout!(this, acnp_bc_out_id_nefflux, n_efflux)
    _set_rout!(this, acnp_bc_out_id_pefflux, p_efflux)

    # ---- Sanity: all gain pools should be used up -----------------------------------
    if npmode != prescribed_n_uptake && abs(n_gain) > 0.1 * calloc_abs_error
        fates_endrun("Allocation should have used up all mass gain; n_gain=$(n_gain)")
    end
    if ppmode != prescribed_p_uptake && abs(p_gain) > 0.01 * calloc_abs_error
        fates_endrun("Allocation should have used up all mass gain; p_gain=$(p_gain)")
    end
    if abs(c_gain) > calloc_abs_error
        fates_endrun("Allocation should have used up all mass gain; c_gain=$(c_gain)")
    end

    # ---- Final tally + net_alloc diagnostics ----------------------------------------
    resp_excess = _rio(this, acnp_bc_inout_id_resp_excess)
    allocated_c = (resp_excess - resp_excess0) + c_efflux
    allocated_n = n_efflux
    allocated_p = p_efflux

    for i in 1:acnp_num_organs
        i_org = l2g_organ_list[i]

        i_var = sp_organ_map_get(g, i_org, carbon12_element)
        this.variables[i_var].net_alloc[1] += (this.variables[i_var].val[1] - state_c0[i_org])
        allocated_c += (this.variables[i_var].val[1] - state_c0[i_org])

        i_var = sp_organ_map_get(g, i_org, nitrogen_element)
        this.variables[i_var].net_alloc[1] += (this.variables[i_var].val[1] - state_n0[i_org])
        allocated_n += (this.variables[i_var].val[1] - state_n0[i_org])

        i_var = sp_organ_map_get(g, i_org, phosphorus_element)
        this.variables[i_var].net_alloc[1] += (this.variables[i_var].val[1] - state_p0[i_org])
        allocated_p += (this.variables[i_var].val[1] - state_p0[i_org])
    end

    if acnp_debug
        if abs(allocated_c - (c_gain0 - c_gain)) > calloc_abs_error ||
           abs(allocated_n - (n_gain0 - n_gain)) > calloc_abs_error ||
           abs(allocated_p - (p_gain0 - p_gain)) > calloc_abs_error
            fates_endrun("CNP allocation scheme did not balance mass. " *
                         "c_gain0=$(c_gain0) allocated_c=$(allocated_c) " *
                         "n_gain0=$(n_gain0) allocated_n=$(allocated_n) " *
                         "p_gain0=$(p_gain0) allocated_p=$(allocated_p)")
        end
    end

    # Report uptake (prescribed) or reset gains (predictive) for diagnostics.
    if npmode == prescribed_n_uptake
        n_gain = n_gain0 - n_gain
    else
        n_gain = n_gain0
    end
    if ppmode == prescribed_p_uptake
        p_gain = p_gain0 - p_gain
    else
        p_gain = p_gain0
    end
    _set_rio!(this, acnp_bc_inout_id_netdn, n_gain)
    _set_rio!(this, acnp_bc_inout_id_netdp, p_gain)

    # Forceful fine-root turnover if l2fr dropped.
    TrimFineRoot!(this)

    return nothing
end
