# PRTAllometricCarbonMod.jl
# Julia port of FATES src/fates/parteh/PRTAllometricCarbonMod.F90
#
# Plant Allocation and Reactive Transport Extensible Hypotheses (PARTEH):
# the CARBON-only, allometric-growth hypothesis. This is a concrete strategy that
# extends the abstract `AbstractPRTVartypes` base (from PRTGenericMod) and
# implements the deferred generic functions:
#
#   * DailyPRT!     => DailyPRTAllometricCarbon
#   * FastPRT!      => FastPRTAllometricCarbon (a no-op stub)
#
# The daily allocation distributes the net daily carbon balance among
# leaf / fine-root / sapwood / storage / structure / reproduction pools, driven
# by allometric targets and the allometry engine (FatesAllometryMod). When
# carbon remains after pools have been pushed to their allometric targets, growth
# proceeds by numerically integrating the pools along the dbh growth curve.
#
# Translation notes:
# - `fates_r8` -> Float64.
# - The Fortran allometry routines return their value/derivative via out-args; the
#   merged Julia versions return `(value, deriv)` tuples (bsap returns
#   `(area, value, deriv)`). Calls are adapted accordingly.
# - Boundary conditions are `prt_bctype` with `Ref`-wrapped `rval`/`ival` (the
#   `nothing`-flushable pointers); inout BCs (dbh, carbon_balance) are read/written
#   through their `Ref`.
# - `prt_global` (generic) and `prt_global_ac` (this hypothesis' descriptor)
#   mirror the Fortran module globals; `prt_global` is set to point at
#   `prt_global_ac` at init, exactly like the Fortran.
# - `regeneration_model` is read via `ed_params().regeneration_model`.
#
# Deps: PRTGenericMod, FatesAllometryMod, FatesIntegratorsMod, FatesConstantsMod,
# PRTParametersMod, EDParamsMod, FatesGlobals.

# -------------------------------------------------------------------------------------
# State variable object indices for this hypothesis.
# -------------------------------------------------------------------------------------
const ac_leaf_c_id   = 1   # Unique object index for leaf carbon
const ac_fnrt_c_id   = 2   # Unique object index for fine-root carbon
const ac_sapw_c_id   = 3   # Unique object index for sapwood carbon
const ac_store_c_id  = 4   # Unique object index for storage carbon
const ac_repro_c_id  = 5   # Unique object index for reproductive carbon
const ac_struct_c_id = 6   # Unique object index for structural carbon
const ac_num_vars    = 6   # THIS MUST MATCH THE LARGEST INDEX ABOVE

# For this hypothesis, we integrate dbh along with the other 6. Since this
# is a boundary condition, we do not add it to the state array, but we do want
# to include it with the integrator array.
const ac_dbh_id             = 7   # This is just used for the integrator
const ac_n_integration_vars = 7

# -------------------------------------------------------------------------------------
# Boundary Conditions (input/output)
# -------------------------------------------------------------------------------------
const ac_bc_inout_id_dbh   = 1    # Plant DBH
const ac_bc_inout_id_netdc = 2    # Index for the net daily C input BC
const ac_num_bc_inout      = 2    # Number of in & output boundary conditions

const ac_bc_in_id_pft     = 1     # Index for the PFT input BC
const ac_bc_in_id_ctrim   = 2     # Index for the canopy trim function
const ac_bc_in_id_lstat   = 3     # Leaf status (on, off, partial abscission)
const ac_bc_in_id_cdamage = 4     # Index for the crowndamage input BC
const ac_bc_in_id_efleaf  = 5     # Elongation factor (leaves)
const ac_bc_in_id_effnrt  = 6     # "Elongation factor" (fine roots)
const ac_bc_in_id_efstem  = 7     # "Elongation factor" (stem)
const ac_num_bc_in        = 7     # Number of input boundary conditions

# There are no purely output boundary conditions
const ac_num_bc_out       = 0     # Number of purely output boundary conditions

# Only 1 coordinate per variable for this hypothesis.
const ac_icd              = 1

# Maximum number of leaf age pools (scratch space).
const ac_max_nleafage     = 10

# -------------------------------------------------------------------------------------
# This is the core type that holds this specific plant reactive transport (PRT)
# module. It extends the generic prt_vartypes carrier; the concrete subtype carries
# the same SoA state fields and dispatches DailyPRT!/FastPRT! to the carbon methods.
# -------------------------------------------------------------------------------------
Base.@kwdef mutable struct callom_prt_vartypes <: AbstractPRTVartypes
    variables::Vector{prt_vartype} = prt_vartype[]   # state variables and fluxes
    bc_inout::Vector{prt_bctype}   = prt_bctype[]    # boundaries that may be changed
    bc_in::Vector{prt_bctype}      = prt_bctype[]    # protected
    bc_out::Vector{prt_bctype}     = prt_bctype[]    # overwritten
    ode_opt_step::Float64          = 0.0
end

# This is the instance of the mapping table and variable definitions. It is only
# allocated once per node and should be read-only everywhere except where it is
# populated in the init routine below. (Fortran: prt_global_ac)
const prt_global_ac = Ref{Union{Nothing,prt_global_type}}(nothing)

# =====================================================================================

"""
    InitPRTGlobalAllometricCarbon!()

Initialize and populate the global descriptor object for the carbon-only
allometric hypothesis: the names/symbols/units of each variable and the
organ/element mapping tables. Called once, very early in the call sequence,
before any plants are initialized. Also points the generic `prt_global` at this
hypothesis' descriptor `prt_global_ac`, mirroring the Fortran.
"""
function InitPRTGlobalAllometricCarbon!()

    g = prt_global_type()

    # The "state descriptor" object holds the names, symbols and units of each
    # variable; allocate one per variable.
    g.state_descriptor = [state_descriptor_type() for _ in 1:ac_num_vars]

    g.hyp_name = "Allometric Carbon Only"
    g.hyp_id   = prt_carbon_allom_hyp

    # Set mapping tables to zero
    ZeroGlobal!(g)

    # The number of leaf age classes is determined from the parameter file
    # (second dimension of leaf_long). Same value as in FatesInterfaceMod.F90.
    nleafage = size(prt_params.leaf_long, 2)

    if nleafage > ac_max_nleafage
        fates_endrun(
            "The allometric carbon PARTEH hypothesis sets a maximum number of " *
            "leaf age classes used for scratch space. The model wants to exceed " *
            "that ($(nleafage) > $(ac_max_nleafage)). Simply increase " *
            "ac_max_nleafage found in fates/PRTAllometricCarbonMod.jl")
    end

    # Register the variables. Each variable is associated with a global organ and
    # element identifier. Leaves are discretized further by age class.
    RegisterVarInGlobal!(g, ac_leaf_c_id,   "Leaf Carbon",          "leaf_c",   leaf_organ,   carbon12_element, nleafage)
    RegisterVarInGlobal!(g, ac_fnrt_c_id,   "Fine Root Carbon",     "fnrt_c",   fnrt_organ,   carbon12_element, ac_icd)
    RegisterVarInGlobal!(g, ac_sapw_c_id,   "Sapwood Carbon",       "sapw_c",   sapw_organ,   carbon12_element, ac_icd)
    RegisterVarInGlobal!(g, ac_store_c_id,  "Storage Carbon",       "store_c",  store_organ,  carbon12_element, ac_icd)
    RegisterVarInGlobal!(g, ac_struct_c_id, "Structural Carbon",    "struct_c", struct_organ, carbon12_element, ac_icd)
    RegisterVarInGlobal!(g, ac_repro_c_id,  "Reproductive Carbon",  "repro_c",  repro_organ,  carbon12_element, ac_icd)

    # Set the array sizes for input and output boundary conditions
    g.num_bc_in    = ac_num_bc_in
    g.num_bc_out   = ac_num_bc_out
    g.num_bc_inout = ac_num_bc_inout
    g.num_vars     = ac_num_vars

    # Store this hypothesis' descriptor and have the generic global point at it.
    prt_global_ac[] = g
    prt_global[]    = g

    return nothing
end

# =====================================================================================

"""
    DailyPRT!(this::callom_prt_vartypes, phase::Integer)

Main daily allocation routine for the carbon-only, allometric-growth hypothesis.

Allocation priorities (split across `phase` 1..3):
  1. Replace leaf/fine-root maintenance turnover (may draw from storage).
  2. Re-coup storage losses if carbon balance is negative, else top off storage.
  3. Replenish the rest of the leaf/fine-root deficit, then push all live pools
     and structure toward their allometric targets.
  4. If carbon remains, grow all pools (incl. structure and reproduction)
     concurrently by integrating along the dbh growth curve.

`dbh` and the daily `carbon_balance` are inout boundary conditions, read from and
written back through their `Ref`s.
"""
function DailyPRT!(this::callom_prt_vartypes, phase::Integer)

    g = get_prt_global()

    # Tunable parameters (Fortran module parameters).
    max_substeps    = 300       # Maximum allowable iterations
    max_trunc_error = 1.0       # Maximum allowable truncation error
    ODESolve        = 2         # 1=RKF45, 2=Euler
    iexp_leaf       = 1         # index of expanding (youngest) leaf age class

    # -----------------------------------------------------------------------------------
    # 0. Copy the input-only boundary conditions into readable local variables.
    # -----------------------------------------------------------------------------------
    ipft        = this.bc_in[ac_bc_in_id_pft].ival[]
    canopy_trim = this.bc_in[ac_bc_in_id_ctrim].rval[]
    leaf_status = this.bc_in[ac_bc_in_id_lstat].ival[]
    crowndamage = this.bc_in[ac_bc_in_id_cdamage].ival[]
    elongf_leaf = this.bc_in[ac_bc_in_id_efleaf].rval[]
    elongf_fnrt = this.bc_in[ac_bc_in_id_effnrt].rval[]
    elongf_stem = this.bc_in[ac_bc_in_id_efstem].rval[]

    # Set some logical flags to simplify "if" blocks
    is_hydecid_dormant =
        (prt_params.stress_decid[ipft] == ihard_stress_decid ||
         prt_params.stress_decid[ipft] == isemi_stress_decid) &&
        (leaf_status == leaves_off || leaf_status == leaves_shedding)
    is_deciduous =
        (prt_params.stress_decid[ipft] == ihard_stress_decid ||
         prt_params.stress_decid[ipft] == isemi_stress_decid) ||
        (prt_params.season_decid[ipft] == itrue)

    nleafage = g.state_descriptor[ac_leaf_c_id].num_pos  # Number of leaf age class

    # -----------------------------------------------------------------------------------
    # 1/2. Simpler local names for the state variables (associate in Fortran).
    #      The leaf variable is a vector (over age positions); others are scalars at icd.
    # -----------------------------------------------------------------------------------
    leaf_c   = this.variables[ac_leaf_c_id].val             # vector (1:nleafage)
    l2fr     = prt_params.allom_l2fr[ipft]

    # inout boundary-condition Refs
    dbh_ref = this.bc_inout[ac_bc_inout_id_dbh].rval
    cb_ref  = this.bc_inout[ac_bc_inout_id_netdc].rval
    dbh            = dbh_ref[]
    carbon_balance = cb_ref[]

    # scalar pool locals (mirror the Fortran scalar associations)
    fnrt_c   = this.variables[ac_fnrt_c_id].val[ac_icd]
    sapw_c   = this.variables[ac_sapw_c_id].val[ac_icd]
    store_c  = this.variables[ac_store_c_id].val[ac_icd]
    repro_c  = this.variables[ac_repro_c_id].val[ac_icd]
    struct_c = this.variables[ac_struct_c_id].val[ac_icd]

    # -----------------------------------------------------------------------------------
    # I. Remember the values of the state variables at the beginning of this routine.
    # -----------------------------------------------------------------------------------
    leaf_c0   = zeros(Float64, ac_max_nleafage)
    @views leaf_c0[1:nleafage] .= leaf_c[1:nleafage]
    fnrt_c0   = fnrt_c
    sapw_c0   = sapw_c
    store_c0  = store_c
    repro_c0  = repro_c
    struct_c0 = struct_c

    # -----------------------------------------------------------------------------------
    # II. Calculate target size of each biomass compartment for the given dbh.
    # -----------------------------------------------------------------------------------
    # Target sapwood biomass according to allometry and trimming [kgC]
    _sapw_area, target_sapw_c, _ = bsap_allom(dbh, ipft, crowndamage, canopy_trim, elongf_stem)
    # Target total above ground biomass in woody/fibrous tissues [kgC]
    target_agw_c, _              = bagw_allom(dbh, ipft, crowndamage, elongf_stem)
    # Target total below ground biomass in woody/fibrous tissues [kgC]
    target_bgw_c, _              = bbgw_allom(dbh, ipft, elongf_stem)
    # Target total dead (structural) biomass [kgC]
    target_struct_c, _           = bdead_allom(target_agw_c, target_bgw_c, target_sapw_c, ipft)
    # Target leaf biomass according to allometry and trimming
    target_leaf_c, _             = bleaf(dbh, ipft, crowndamage, canopy_trim, elongf_leaf)
    # Target fine-root biomass according to allometry and trimming [kgC]
    target_fnrt_c, _             = bfineroot(dbh, ipft, canopy_trim, l2fr, elongf_fnrt)
    # Target storage carbon [kgC]
    target_store_c, _            = bstore_allom(dbh, ipft, crowndamage, canopy_trim)

    # -----------------------------------------------------------------------------------
    # II 1/2. Update targets based on leaf elongation factor for deciduous dormancy.
    # -----------------------------------------------------------------------------------
    if is_hydecid_dormant
        target_leaf_c   = 0.0
        target_fnrt_c   = 0.0
        target_sapw_c   = 0.0
        target_struct_c = 0.0
        target_store_c  = target_store_c + max(0.0, carbon_balance)
    end

    # -----------------------------------------------------------------------------------
    # The following blocks allocate carbon.
    # -----------------------------------------------------------------------------------
    if phase == 1
        # -------------------------------------------------------------------------------
        # Phase 1: Replace losses, push pools towards targets
        # -------------------------------------------------------------------------------
        # III. Prioritize some amount of carbon to replace leaf/root turnover.
        if is_hydecid_dormant
            # Drought deciduous, dormant. Set demands to leaves and roots to zero.
            leaf_c_demand = 0.0
            fnrt_c_demand = 0.0
        elseif is_deciduous
            # Cold deciduous, or drought deciduous w/ leaves on. Maintain roots.
            leaf_c_demand = 0.0
            fnrt_c_demand = max(0.0,
                prt_params.leaf_stor_priority[ipft] * this.variables[ac_fnrt_c_id].turnover[ac_icd])
        else
            # Evergreen PFT. Try to meet demands for both leaves and fine roots.
            leaf_c_demand = max(0.0,
                prt_params.leaf_stor_priority[ipft] * sum(this.variables[ac_leaf_c_id].turnover))
            fnrt_c_demand = max(0.0,
                prt_params.leaf_stor_priority[ipft] * this.variables[ac_fnrt_c_id].turnover[ac_icd])
        end

        total_c_demand = leaf_c_demand + fnrt_c_demand

        if total_c_demand > nearzero
            # Pay this even if we don't have the carbon, but don't pay so much that
            # storage+carbon_balance can't cover it.
            allocation_factor = max(0.0, min(1.0, (store_c + carbon_balance) / total_c_demand))

            leaf_c_flux = leaf_c_demand * allocation_factor
            fnrt_c_flux = fnrt_c_demand * allocation_factor

            # Add carbon to the youngest age pool and fine roots
            leaf_c[iexp_leaf] = leaf_c[iexp_leaf] + leaf_c_flux
            fnrt_c            = fnrt_c + fnrt_c_flux

            # Remove fluxes from carbon balance (may become negative -> drawn from
            # storage in the next step).
            carbon_balance = carbon_balance - (leaf_c_flux + fnrt_c_flux)
        end

        # IV. If carbon balance is negative, re-coup losses from storage; if positive
        #     give some love to storage carbon.
        if carbon_balance < 0.0
            # store_c_flux is negative, so store_c is depleted
            store_c_flux   = carbon_balance
            carbon_balance = carbon_balance - store_c_flux
            store_c        = store_c + store_c_flux
        else
            # Accumulate some carbon in storage; if depleted, aim to increase storage
            # but not replenish completely so some carbon remains for growth.
            store_below_target    = max(0.0, target_store_c - store_c)
            store_target_fraction = max(0.0, store_c / target_store_c)

            store_c_flux = min(store_below_target, carbon_balance *
                max(exp(-1.0 * store_target_fraction^4.0) - exp(-1.0), 0.0))

            carbon_balance = carbon_balance - store_c_flux
            store_c        = store_c + store_c_flux
        end

    elseif phase == 2

        # V. If carbon is still available, prioritize allocation to replace the rest
        #    of the leaf/fineroot deficit. (carbon_balance is >= 0 beyond this point.)
        leaf_below_target  = max(0.0, target_leaf_c - sum(@view leaf_c[1:nleafage]))
        fnrt_below_target  = max(0.0, target_fnrt_c - fnrt_c)

        total_below_target = leaf_below_target + fnrt_below_target

        if (carbon_balance > nearzero) && (total_below_target > nearzero)
            allocation_factor = min(1.0, carbon_balance / total_below_target)

            leaf_c_flux = leaf_below_target * allocation_factor
            fnrt_c_flux = fnrt_below_target * allocation_factor

            leaf_c[iexp_leaf] = leaf_c[iexp_leaf] + leaf_c_flux
            fnrt_c            = fnrt_c + fnrt_c_flux

            carbon_balance = carbon_balance - (leaf_c_flux + fnrt_c_flux)
        end

        # VI. If carbon is still available, push all live pools back towards allometry
        #     (upwards only).
        if carbon_balance > nearzero
            leaf_below_target  = max(target_leaf_c - sum(@view leaf_c[1:nleafage]), 0.0)
            fnrt_below_target  = max(target_fnrt_c - fnrt_c, 0.0)
            sapw_below_target  = max(target_sapw_c - sapw_c, 0.0)
            store_below_target = max(target_store_c - store_c, 0.0)

            total_below_target = leaf_below_target + fnrt_below_target +
                                 sapw_below_target + store_below_target

            if total_below_target > nearzero
                allocation_factor = min(1.0, carbon_balance / total_below_target)

                leaf_c_flux  = leaf_below_target  * allocation_factor
                fnrt_c_flux  = fnrt_below_target  * allocation_factor
                sapw_c_flux  = sapw_below_target  * allocation_factor
                store_c_flux = store_below_target * allocation_factor

                leaf_c[iexp_leaf] = leaf_c[iexp_leaf] + leaf_c_flux
                fnrt_c            = fnrt_c + fnrt_c_flux
                sapw_c            = sapw_c + sapw_c_flux
                store_c           = store_c + store_c_flux

                carbon_balance = carbon_balance -
                    (leaf_c_flux + fnrt_c_flux + sapw_c_flux + store_c_flux)
            end
        end

        # VII. If carbon is still available, replenish the structural pool to allometry.
        if carbon_balance > nearzero
            struct_below_target = max(target_struct_c - struct_c, 0.0)
            if struct_below_target > 0.0
                struct_c_flux  = min(carbon_balance, struct_below_target)
                carbon_balance = carbon_balance - struct_c_flux
                struct_c       = struct_c + struct_c_flux
            end
        end

    elseif phase == 3
        # VII 1/2: Semi-deciduous plants with positive carbon balance but losing
        #          leaves should not invest in growth: stash the carbon to storage.
        if leaf_status == leaves_off || leaf_status == leaves_shedding
            store_c_flux   = carbon_balance
            carbon_balance = carbon_balance - store_c_flux
            store_c        = store_c + store_c_flux
        end

        if carbon_balance > calloc_abs_error

            # VIII. If carbon is yet still available, grow all pools at or below target
            #       (including structure and reproduction) by integrating along dbh.
            intgr_params = fill(un_initialized, ac_num_bc_in)
            intgr_params[ac_bc_in_id_ctrim]   = this.bc_in[ac_bc_in_id_ctrim].rval[]
            intgr_params[ac_bc_in_id_pft]     = Float64(this.bc_in[ac_bc_in_id_pft].ival[])
            intgr_params[ac_bc_in_id_cdamage] = Float64(this.bc_in[ac_bc_in_id_cdamage].ival[])
            intgr_params[ac_bc_in_id_lstat]   = Float64(this.bc_in[ac_bc_in_id_lstat].ival[])
            intgr_params[ac_bc_in_id_efleaf]  = this.bc_in[ac_bc_in_id_efleaf].rval[]
            intgr_params[ac_bc_in_id_effnrt]  = this.bc_in[ac_bc_in_id_effnrt].rval[]
            intgr_params[ac_bc_in_id_efstem]  = this.bc_in[ac_bc_in_id_efstem].rval[]

            # Check actual pools against targets; allow pools above target (grow=false
            # there) so the plant can grow into them.
            grow_leaf, grow_fnrt, grow_sapw, grow_store, grow_struct =
                _TargetAllometryCheck(sum(@view leaf_c0[1:nleafage]), fnrt_c0, sapw_c0, store_c0, struct_c0,
                    sum(@view leaf_c[1:nleafage]), fnrt_c, sapw_c, store_c, struct_c,
                    target_leaf_c, target_fnrt_c, target_sapw_c, target_store_c, target_struct_c,
                    carbon_balance, elongf_leaf, elongf_fnrt, elongf_stem, ipft, leaf_status)

            # Initialize the adaptive integrator arrays and flags.
            ierr   = 1
            totalC = carbon_balance
            nsteps = 0

            c_pool = zeros(Float64, ac_n_integration_vars)
            c_mask = fill(false, ac_n_integration_vars)

            c_pool[ac_leaf_c_id]   = sum(@view leaf_c[1:nleafage])
            c_pool[ac_fnrt_c_id]   = fnrt_c
            c_pool[ac_sapw_c_id]   = sapw_c
            c_pool[ac_store_c_id]  = store_c
            c_pool[ac_struct_c_id] = struct_c
            c_pool[ac_repro_c_id]  = repro_c
            c_pool[ac_dbh_id]      = dbh

            # Only grow leaves if in a leaf-on status; interrupt growth of all tissues
            # when drought-deciduous and dormant.
            if is_hydecid_dormant
                c_mask[ac_leaf_c_id]   = false
                c_mask[ac_fnrt_c_id]   = false
                c_mask[ac_sapw_c_id]   = false
                c_mask[ac_struct_c_id] = false
            else
                c_mask[ac_leaf_c_id]   = (leaf_status == leaves_on) ? grow_leaf : false
                c_mask[ac_fnrt_c_id]   = grow_fnrt
                c_mask[ac_sapw_c_id]   = grow_sapw
                c_mask[ac_struct_c_id] = grow_struct
            end
            c_mask[ac_store_c_id]  = grow_store
            c_mask[ac_repro_c_id]  = true   # Always calculate reproduction on growth
            c_mask[ac_dbh_id]      = true   # Always increment dbh on growth step

            # Euler: try to span the entire integration window in the first step.
            if ODESolve == 2
                this.ode_opt_step = totalC
            end

            c_pool_out = zeros(Float64, ac_n_integration_vars)

            step_pass = false
            while ierr != 0

                deltaC = min(totalC, this.ode_opt_step)

                if ODESolve == 1
                    this.ode_opt_step, step_pass =
                        RKF45(AllomCGrowthDeriv, c_pool, c_mask, deltaC, totalC,
                              max_trunc_error, intgr_params, c_pool_out)

                elseif ODESolve == 2
                    Euler(AllomCGrowthDeriv, c_pool, c_mask, deltaC, totalC,
                          intgr_params, c_pool_out)

                    # After the step, check how close the pools are to the allometric
                    # targets for the new dbh. If too far, halve the step and retry.
                    step_pass = CheckIntegratedAllometries(c_pool_out[ac_dbh_id], ipft,
                        crowndamage, canopy_trim, elongf_leaf, elongf_fnrt, elongf_stem, l2fr,
                        c_pool_out[ac_leaf_c_id], c_pool_out[ac_fnrt_c_id], c_pool_out[ac_sapw_c_id],
                        c_pool_out[ac_store_c_id], c_pool_out[ac_struct_c_id],
                        c_mask[ac_leaf_c_id], c_mask[ac_fnrt_c_id], c_mask[ac_sapw_c_id],
                        c_mask[ac_store_c_id], c_mask[ac_struct_c_id], max_trunc_error)
                    if step_pass
                        this.ode_opt_step = deltaC
                    else
                        this.ode_opt_step = 0.5 * deltaC
                    end
                else
                    fates_endrun("An integrator was chosen that does not exist. ODESolve = $(ODESolve)")
                end

                nsteps = nsteps + 1

                if step_pass  # If true, the step is accepted
                    totalC = totalC - deltaC
                    c_pool .= c_pool_out
                end

                if nsteps > max_substeps
                    fates_endrun(
                        "Plant Growth Integrator could not find a solution in less " *
                        "than $(max_substeps) tries. Aborting! leaf_status=$(leaf_status) " *
                        "carbon_balance=$(carbon_balance) deltaC=$(deltaC) totalC=$(totalC) " *
                        "crowndamage=$(crowndamage)")
                end

                # totalC should eventually be whittled down to near zero.
                if (totalC < calloc_abs_error) && step_pass
                    ierr          = 0
                    leaf_c_flux   = c_pool[ac_leaf_c_id]   - sum(@view leaf_c[1:nleafage])
                    fnrt_c_flux   = c_pool[ac_fnrt_c_id]   - fnrt_c
                    sapw_c_flux   = c_pool[ac_sapw_c_id]   - sapw_c
                    store_c_flux  = c_pool[ac_store_c_id]  - store_c
                    struct_c_flux = c_pool[ac_struct_c_id] - struct_c
                    repro_c_flux  = c_pool[ac_repro_c_id]  - repro_c

                    # Adjust flux partitions to match the remaining carbon balance.
                    flux_adj = carbon_balance / (leaf_c_flux + fnrt_c_flux + sapw_c_flux +
                                                 store_c_flux + struct_c_flux + repro_c_flux)

                    leaf_c_flux   = leaf_c_flux   * flux_adj
                    fnrt_c_flux   = fnrt_c_flux   * flux_adj
                    sapw_c_flux   = sapw_c_flux   * flux_adj
                    store_c_flux  = store_c_flux  * flux_adj
                    struct_c_flux = struct_c_flux * flux_adj
                    repro_c_flux  = repro_c_flux  * flux_adj

                    leaf_c[iexp_leaf] = leaf_c[iexp_leaf] + leaf_c_flux
                    fnrt_c            = fnrt_c + fnrt_c_flux
                    sapw_c            = sapw_c + sapw_c_flux
                    store_c           = store_c + store_c_flux
                    struct_c          = struct_c + struct_c_flux
                    repro_c           = repro_c + repro_c_flux

                    carbon_balance = carbon_balance -
                        (leaf_c_flux + fnrt_c_flux + sapw_c_flux +
                         store_c_flux + struct_c_flux + repro_c_flux)

                    dbh = c_pool[ac_dbh_id]

                    if abs(carbon_balance) > calloc_abs_error
                        fates_endrun(
                            "carbon conservation error while integrating pools along " *
                            "allometric curve. carbon_balance=$(carbon_balance) totalC=$(totalC)")
                    end
                end

            end  # while ierr != 0
        end  # if carbon_balance > calloc_abs_error
    end  # phase select

    # -----------------------------------------------------------------------------------
    # Write the scalar pool locals back into the state arrays, and (re)write dbh and
    # carbon_balance through their inout BC Refs.
    # -----------------------------------------------------------------------------------
    this.variables[ac_fnrt_c_id].val[ac_icd]   = fnrt_c
    this.variables[ac_sapw_c_id].val[ac_icd]   = sapw_c
    this.variables[ac_store_c_id].val[ac_icd]  = store_c
    this.variables[ac_repro_c_id].val[ac_icd]  = repro_c
    this.variables[ac_struct_c_id].val[ac_icd] = struct_c
    dbh_ref[] = dbh
    cb_ref[]  = carbon_balance

    # Track the net allocations and transport from this routine (AgeLeaves! handles
    # tracking allocation through aging).
    this.variables[ac_leaf_c_id].net_alloc[ac_icd] +=
        (leaf_c[ac_icd] - leaf_c0[ac_icd])
    this.variables[ac_fnrt_c_id].net_alloc[ac_icd] +=
        (fnrt_c - fnrt_c0)
    this.variables[ac_sapw_c_id].net_alloc[ac_icd] +=
        (sapw_c - sapw_c0)
    this.variables[ac_store_c_id].net_alloc[ac_icd] +=
        (store_c - store_c0)
    this.variables[ac_repro_c_id].net_alloc[ac_icd] +=
        (repro_c - repro_c0)
    this.variables[ac_struct_c_id].net_alloc[ac_icd] +=
        (struct_c - struct_c0)

    return nothing
end

# =====================================================================================

"""
    AllomCGrowthDeriv(c_pools, c_mask, cbalance, intgr_params) -> dCdx

Derivatives of the carbon pools with respect to the amount of carbon balance
(the independent variable). Based completely off allometry. `c_pools` carries
`dbh, leaf, fnrt, sapw, store, struct, repro` (indexed by the `ac_*_id` consts).
`c_mask` flags which pools are active. Returns a vector `dCdx` (same length as
`c_pools`) giving change in each pool per change in total allocatable carbon
[kgC/kgC].
"""
function AllomCGrowthDeriv(c_pools::AbstractVector, c_mask::AbstractVector{Bool},
                           cbalance::Real, intgr_params::AbstractVector)

    dbh    = c_pools[ac_dbh_id]
    # cleaf/cfnrt/... read as needed; the unused ones mirror the Fortran associate

    mask_leaf  = c_mask[ac_leaf_c_id]
    mask_fnrt  = c_mask[ac_fnrt_c_id]
    mask_sap   = c_mask[ac_sapw_c_id]
    mask_store = c_mask[ac_store_c_id]

    canopy_trim = intgr_params[ac_bc_in_id_ctrim]
    ipft        = Int(intgr_params[ac_bc_in_id_pft])
    elongf_leaf = intgr_params[ac_bc_in_id_efleaf]
    elongf_fnrt = intgr_params[ac_bc_in_id_effnrt]
    elongf_stem = intgr_params[ac_bc_in_id_efstem]
    crowndamage = Int(intgr_params[ac_bc_in_id_cdamage])
    l2fr        = prt_params.allom_l2fr[ipft]

    _ct_leaf,  ct_dleafdd = bleaf(dbh, ipft, crowndamage, canopy_trim, elongf_leaf)
    _ct_fnrt,  ct_dfnrtdd = bfineroot(dbh, ipft, canopy_trim, l2fr, elongf_fnrt)
    _sapw_area, _ct_sap, ct_dsapdd = bsap_allom(dbh, ipft, crowndamage, canopy_trim, elongf_stem)
    ct_agw,    ct_dagwdd  = bagw_allom(dbh, ipft, crowndamage, elongf_stem)
    ct_bgw,    ct_dbgwdd  = bbgw_allom(dbh, ipft, elongf_stem)
    _ct_dead,  ct_ddeaddd = bdead_allom(ct_agw, ct_bgw, _ct_sap, ipft;
                                        dbagwdd = ct_dagwdd, dbbgwdd = ct_dbgwdd, dbsapdd = ct_dsapdd)
    _ct_store, ct_dstoredd = bstore_allom(dbh, ipft, crowndamage, canopy_trim)

    # Reproductive allocation fraction. If the TRS is off, or the plant is a shrub
    # or grass (max-height dbh < min_max_dbh_for_trees), use FATES's default.
    reg_model = ed_params().regeneration_model
    if reg_model == default_regeneration ||
       prt_params.allom_dbh_maxheight[ipft] < min_max_dbh_for_trees

        if dbh <= prt_params.dbh_repro_threshold[ipft]   # cap on leaf biomass
            repro_fraction = prt_params.seed_alloc[ipft]
        else
            repro_fraction = prt_params.seed_alloc[ipft] + prt_params.seed_alloc_mature[ipft]
        end

    elseif (reg_model == TRS_regeneration || reg_model == TRS_no_seedling_dyn) &&
           prt_params.allom_dbh_maxheight[ipft] > min_max_dbh_for_trees

        repro_fraction = prt_params.seed_alloc[ipft] *
            (exp(prt_params.repro_alloc_b[ipft] + prt_params.repro_alloc_a[ipft] * dbh * mm_per_cm) /
             (1.0 + exp(prt_params.repro_alloc_b[ipft] + prt_params.repro_alloc_a[ipft] * dbh * mm_per_cm)))
    else
        fates_endrun(
            "unknown seed allocation and regeneration model, exiting. " *
            "regeneration_model: $(reg_model)")
    end

    dCdx = zeros(eltype(c_pools), length(c_pools))

    ct_dtotaldd = ct_ddeaddd
    if mask_leaf;  ct_dtotaldd = ct_dtotaldd + ct_dleafdd;  end
    if mask_fnrt;  ct_dtotaldd = ct_dtotaldd + ct_dfnrtdd;  end
    if mask_sap;   ct_dtotaldd = ct_dtotaldd + ct_dsapdd;   end
    if mask_store; ct_dtotaldd = ct_dtotaldd + ct_dstoredd; end

    # If all growth rates reach zero (asymptotic/capped allometry), give any
    # remaining carbon to reproduction. Fortran compares against `tiny()`, the
    # smallest positive normal Float64 -> `floatmin(Float64)`.
    if ct_dtotaldd <= floatmin(Float64)
        dCdx[ac_struct_c_id] = 0.0
        dCdx[ac_dbh_id]      = 0.0
        dCdx[ac_leaf_c_id]   = 0.0
        dCdx[ac_fnrt_c_id]   = 0.0
        dCdx[ac_sapw_c_id]   = 0.0
        dCdx[ac_store_c_id]  = 0.0
        dCdx[ac_repro_c_id]  = 1.0
    else
        dCdx[ac_struct_c_id] = (ct_ddeaddd / ct_dtotaldd) * (1.0 - repro_fraction)
        dCdx[ac_dbh_id]      = (1.0 / ct_dtotaldd) * (1.0 - repro_fraction)

        dCdx[ac_leaf_c_id]  = mask_leaf  ? (ct_dleafdd  / ct_dtotaldd) * (1.0 - repro_fraction) : 0.0
        dCdx[ac_fnrt_c_id]  = mask_fnrt  ? (ct_dfnrtdd  / ct_dtotaldd) * (1.0 - repro_fraction) : 0.0
        dCdx[ac_sapw_c_id]  = mask_sap   ? (ct_dsapdd   / ct_dtotaldd) * (1.0 - repro_fraction) : 0.0
        dCdx[ac_store_c_id] = mask_store ? (ct_dstoredd / ct_dtotaldd) * (1.0 - repro_fraction) : 0.0

        dCdx[ac_repro_c_id] = repro_fraction
    end

    return dCdx
end

# =====================================================================================

"""
    _TargetAllometryCheck(b0_*, ba_*, bt_*, carbon_balance, elongf_*, ipft, leaf_status)
        -> (grow_leaf, grow_fnrt, grow_sapw, grow_store, grow_struct)

Verify each pool is on (or below) its allometric target before the growth step.
`b0_*` are the initial values, `ba_*` the current (actual), `bt_*` the targets.
A pool is "fine" if it is not more than `calloc_abs_error` above its target. If
all pools look fine, returns growth flags (true where actual <= target within
tolerance). If any pool looks not-fine, errors with a detailed report.
"""
function _TargetAllometryCheck(b0_leaf, b0_fnrt, b0_sapw, b0_store, b0_struct,
                               ba_leaf, ba_fnrt, ba_sapw, ba_store, ba_struct,
                               bt_leaf, bt_fnrt, bt_sapw, bt_store, bt_struct,
                               carbon_balance, elongf_leaf, elongf_fnrt, elongf_stem,
                               ipft, leaf_status)

    # Test whether each pool looks reasonable (not above target beyond tolerance).
    fine_leaf   = (bt_leaf   - ba_leaf  ) <= calloc_abs_error
    fine_fnrt   = (bt_fnrt   - ba_fnrt  ) <= calloc_abs_error
    fine_sapw   = (bt_sapw   - ba_sapw  ) <= calloc_abs_error
    fine_store  = (bt_store  - ba_store ) <= calloc_abs_error
    fine_struct = (bt_struct - ba_struct) <= calloc_abs_error
    all_fine    = fine_leaf && fine_fnrt && fine_sapw && fine_store && fine_struct

    if all_fine
        grow_leaf   = (ba_leaf   - bt_leaf  ) <= calloc_abs_error
        grow_fnrt   = (ba_fnrt   - bt_fnrt  ) <= calloc_abs_error
        grow_sapw   = (ba_sapw   - bt_sapw  ) <= calloc_abs_error
        grow_store  = (ba_store  - bt_store ) <= calloc_abs_error
        grow_struct = (ba_struct - bt_struct) <= calloc_abs_error
        return grow_leaf, grow_fnrt, grow_sapw, grow_store, grow_struct
    else
        fates_endrun(
            "At least one tissue is not on-allometry at the growth step. " *
            "Biomass (initial|current|target|on-allometry): " *
            "Leaf [$(b0_leaf)|$(ba_leaf)|$(bt_leaf)|$(fine_leaf)] " *
            "Fnrt [$(b0_fnrt)|$(ba_fnrt)|$(bt_fnrt)|$(fine_fnrt)] " *
            "Sapw [$(b0_sapw)|$(ba_sapw)|$(bt_sapw)|$(fine_sapw)] " *
            "Store [$(b0_store)|$(ba_store)|$(bt_store)|$(fine_store)] " *
            "Struct [$(b0_struct)|$(ba_struct)|$(bt_struct)|$(fine_struct)] " *
            "PFT=$(ipft) leaf_status=$(leaf_status) " *
            "elongf_leaf=$(elongf_leaf) elongf_fnrt=$(elongf_fnrt) elongf_stem=$(elongf_stem) " *
            "carbon_balance=$(carbon_balance) calloc_abs_error=$(calloc_abs_error)")
    end
end

# =====================================================================================

"""
    FastPRT!(this::callom_prt_vartypes)

Fast (sub-daily) reactive transport for the carbon-only allometric hypothesis.
This is a stub: there are currently no fast-timestep processes in this hypothesis.
"""
function FastPRT!(::callom_prt_vartypes)
    # This routine does nothing (no fast-timestep processes in carbon-only RT).
    return nothing
end
