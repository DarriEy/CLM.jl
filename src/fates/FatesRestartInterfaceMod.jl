# FatesRestartInterfaceMod.jl
# Julia port of FATES src/fates/main/FatesRestartInterfaceMod.F90
#
# The `fates_restart_interface_type` owns the registry of restart variables
# (a `Vector{fates_restart_variable_type}`), the dimension/dim-kind bookkeeping
# (2 dims: cohort, column/site; 4 dim-kinds: cohort-int/r8, site-int/r8) and the
# load-bearing serialization that walks the site -> patch -> cohort demographic
# linked lists, PACKING the state into the flat restart vectors
# (`set_restart_vectors!`) and UNPACKING it back into freshly-rebuilt lists
# (`create_patchcohort_structure!` allocates the skeleton from the per-site /
# per-patch COUNTS, then `get_restart_vectors!` fills the fields and asserts the
# round-trip counts).
#
# Naming conventions (kept verbatim from Fortran):
#   <use_case>_<description>_<dimension>, e.g. ir_npatch_si = restart-var INDEX
#   into the registry for "number of patches per site". si/pa/co/ft/cl/... =
#   the site/patch/cohort/functional-type/canopy-layer dimensions.
#
# Index bookkeeping (the skeleton):
#   io_idx_si      one slot per site (restart_map.site_index[s])
#   io_idx_co_1st  per-patch base address; +fates_maxElementsPerPatch per patch
#   io_idx_co      per-cohort cursor; resets to io_idx_co_1st at each new patch,
#                  +1 per cohort within a patch (shortest -> taller walk).
# COUNTS driving the unpack: fates_PatchesPerSite (ir_npatch_si) and
# fates_CohortsPerPatch (ir_ncohort_pa).
#
# Scope of THIS port (Batch 18): the registry is reproduced in FULL and
# faithfully (every set_restart_var / RegisterCohortVector / DefineRMeanRestartVar
# / DefinePRTRestartVars call is here, gated by the same hlm flags). The pack /
# unpack routines port the DEMOGRAPHIC SKELETON + load-bearing fields exactly
# (patch area/age/land-use/nocomp/livegrass + cohort n/dbh/pft/height/coage/
# canopy_layer/crowndamage/status/trim/l2fr/... + the PatchesPerSite /
# CohortsPerPatch counts and their consistency checks). Batch-18-followup R1/R2/R4/R5
# now also pack/unpack the cohort diagnostic+mortality scalars (R1), patch
# fuel/scorch + ground-albedo + site hydraulics scalars (R2), the patch litter
# blocks (R4), and the site demographic/mortality/damage/flux-diagnostic arrays
# (R5), and update_3dpatch_radiation! re-derives the 3D albedo BC (R3).
# Batch-18-followup R6/R7/R8 complete the round-trip: R6 packs/unpacks the
# patch running means (tveg24 / tveg_lpa / tveg_longterm via the ported
# Get/SetRMeanRestartVar accessors) + the site disturbance-rate diagnostic; R7
# packs/unpacks the per-cohort PRT (PARTEH) pools (val + turnover/net_alloc/
# burned) — the restart cohort's `prt` is re-created by InitPRTObject! in
# create_patchcohort_structure!, then InitPRTBoundaryConditions /
# UpdateCohortBioPhysRates finalize it on unpack; R8 packs/unpacks the cohort
# hydraulics (th_ag/th_troot/th_aroot/errh2o) gated on hlm_use_planthydro — the
# restart cohort's `co_hydr` is re-allocated by InitHydrCohort and finalized by
# UpdatePlantPsiFTCFromTheta!. The full restart fill is now complete.
#
# Deps: FatesRestartVariableType.jl (fates_restart_variable_type, Init!/Flush!),
#       FatesIODimensionsMod (fates_io_dimension_type, fates_bounds_type,
#         dimname_cohort/column, Init!/SetThreadBounds!), FatesIOVariableKindMod
#       (fates_io_variable_kind_type, site/cohort kind constants, iotype_index,
#        Init!/set_active!), FatesGlobals (fates_log/fates_endrun),
#       FatesUtilsMod (check_hlm_list), FatesInterfaceTypesMod (hlm_* flags,
#        numpft, fates_maxElementsPerPatch/Site, hlm_name), EDTypesMod
#       (ed_site_type, area), FatesPatchMod (fates_patch_type, Create!),
#       FatesCohortMod (fates_cohort_type), EDCohortDynamicsMod (insert_cohort,
#        count_cohorts), EDPatchDynamicsMod (set_patchno),
#       FatesSizeAgeTypeIndicesMod (get_age_class_index), PRTGenericMod
#       (prt_global, prt_cnp_flex_allom_hyp, num_elements), EDParamsMod
#       (nlevleaf, regeneration_model), FatesHydraulicsMemMod (nshell,
#        n_hypool_ag/troot, nlevsoi_hyd_max), FatesConstantsMod (N_DIST_TYPES,
#        n_landuse_cats), FatesRadiationMemMod (num_swb).

# ---------------------------------------------------------------------------
# Module-level constants (mirror the Fortran parameters)
# ---------------------------------------------------------------------------
const fates_restart_num_dimensions = 2   # (cohort, column)
const fates_restart_num_dim_kinds  = 4   # (cohort-int, cohort-r8, site-int, site-r8)

const old_cohort = 0
const new_cohort = 1

const flushinvalid = -9999.0
const flushzero    = 0.0
const flushone     = 1.0

# Local debug flag
const restart_debug = false

# ---------------------------------------------------------------------------
# restart_map_type — maps FATES site/cohort indices to HIO array positions.
# Allocated by thread (clump). One entry per thread holds the per-site base
# addresses into the flat restart vectors.
# ---------------------------------------------------------------------------
Base.@kwdef mutable struct restart_map_type
    site_index::Vector{Int}    = Int[]  # maps site index -> HIO site position
    cohort1_index::Vector{Int} = Int[]  # maps site index -> HIO cohort 1st position
end

# ---------------------------------------------------------------------------
# fates_restart_interface_type — the registry + dim bookkeeping + index handles.
#
# The `ir_*` index handles (Fortran module-level integers) become a Dict so the
# many registration sites can stash/lookup the integer index of each registered
# restart variable by symbolic name without 200 individual fields. The
# load-bearing handles used by pack/unpack are also surfaced as direct fields
# for clarity and speed.
# ---------------------------------------------------------------------------
Base.@kwdef mutable struct fates_restart_interface_type
    rvars::Vector{fates_restart_variable_type} = fates_restart_variable_type[]
    num_restart_vars_::Int = 0

    # 4 dim-kinds (cohort-r8, site-r8, cohort-int, site-int) and 2 dim-bounds
    # (cohort, column). Allocated in Init!.
    dim_kinds::Vector{fates_io_variable_kind_type}  = fates_io_variable_kind_type[]
    dim_bounds::Vector{fates_io_dimension_type}     = fates_io_dimension_type[]

    restart_map::Vector{restart_map_type} = restart_map_type[]

    cohort_index_::Int = 0
    column_index_::Int = 0

    # ir_* index handles: name -> integer index into `rvars`. 0 == not used.
    ir::Dict{Symbol,Int} = Dict{Symbol,Int}()
    # base index of all PRT variables
    ir_prt_base::Int = 0
end

# small logging helper
_ri_log_println(args...) = (println(stderr, args...); nothing)

# Accessors for the ir_* handles, mirroring the named Fortran integers.
ir_get(this::fates_restart_interface_type, name::Symbol) = get(this.ir, name, 0)

# =====================================================================================
# Init / dim setup
# =====================================================================================

"""
    Init!(this, num_threads, fates_bounds)

Allocate the two IO dimension bounds (cohort, column) and the restart_map
(per-thread). Mirrors the Fortran `Init`.
"""
function Init!(this::fates_restart_interface_type, num_threads::Integer,
               fates_bounds::fates_bounds_type)
    this.dim_kinds  = [fates_io_variable_kind_type() for _ in 1:fates_restart_num_dim_kinds]
    this.dim_bounds = [fates_io_dimension_type() for _ in 1:fates_restart_num_dimensions]

    dim_count = 0

    dim_count += 1
    set_cohort_index!(this, dim_count)
    Init!(this.dim_bounds[dim_count], dimname_cohort, num_threads,
          fates_bounds.cohort_begin, fates_bounds.cohort_end)

    dim_count += 1
    set_column_index!(this, dim_count)
    Init!(this.dim_bounds[dim_count], dimname_column, num_threads,
          fates_bounds.column_begin, fates_bounds.column_end)

    # Allocate the mapping between FATES indices and the IO indices
    this.restart_map = [restart_map_type() for _ in 1:num_threads]
    return nothing
end

"""
    SetThreadBoundsEach!(this, thread_index, thread_bounds)

Set the per-thread (clump) bounds for the cohort and column dimensions.
"""
function SetThreadBoundsEach!(this::fates_restart_interface_type, thread_index::Integer,
                              thread_bounds::fates_bounds_type)
    idx = cohort_index(this)
    SetThreadBounds!(this.dim_bounds[idx], thread_index,
                     thread_bounds.cohort_begin, thread_bounds.cohort_end)

    idx = column_index(this)
    SetThreadBounds!(this.dim_bounds[idx], thread_index,
                     thread_bounds.column_begin, thread_bounds.column_end)
    return nothing
end

set_cohort_index!(this::fates_restart_interface_type, index::Integer) = (this.cohort_index_ = index; nothing)
cohort_index(this::fates_restart_interface_type) = this.cohort_index_
set_column_index!(this::fates_restart_interface_type, index::Integer) = (this.column_index_ = index; nothing)
column_index(this::fates_restart_interface_type) = this.column_index_

num_restart_vars(this::fates_restart_interface_type) = this.num_restart_vars_

# =====================================================================================

"""
    init_dim_kinds_maps!(this)

Initialize the 4 dim-kind descriptors: CO_R8, SI_R8, CO_INT, SI_INT (all 1D).
"""
function init_dim_kinds_maps!(this::fates_restart_interface_type)
    Init!(this.dim_kinds[1], cohort_r8, 1)
    Init!(this.dim_kinds[2], site_r8, 1)
    Init!(this.dim_kinds[3], cohort_int, 1)
    Init!(this.dim_kinds[4], site_int, 1)
    return nothing
end

"""
    set_dim_indices!(this, dk_name, idim, dim_index)

Wire a dim-kind's `idim`-th dimension to dimension `dim_index` and set its size.
"""
function set_dim_indices!(this::fates_restart_interface_type, dk_name::AbstractString,
                          idim::Integer, dim_index::Integer)
    ityp = iotype_index(strip(dk_name), fates_restart_num_dim_kinds, this.dim_kinds)

    if this.dim_kinds[ityp].ndims < idim
        _ri_log_println("Trying to define dimension size to a dim-type structure")
        _ri_log_println("but the dimension index does not exist")
        fates_endrun("set_dim_indices!: idim $idim > ndims for $dk_name")
    end

    if idim == 1
        this.dim_kinds[ityp].dim1_index = dim_index
    elseif idim == 2
        this.dim_kinds[ityp].dim2_index = dim_index
    end

    this.dim_kinds[ityp].dimsize[idim] =
        this.dim_bounds[dim_index].upper_bound - this.dim_bounds[dim_index].lower_bound + 1
    return nothing
end

"""
    assemble_restart_output_types!(this)

Build the dim-kind maps and wire each kind to its dimension (cohort kinds ->
cohort dim, site kinds -> column dim). Mirrors the Fortran routine.
"""
function assemble_restart_output_types!(this::fates_restart_interface_type)
    init_dim_kinds_maps!(this)

    set_dim_indices!(this, cohort_r8, 1, cohort_index(this))
    set_dim_indices!(this, cohort_int, 1, cohort_index(this))

    set_dim_indices!(this, site_r8, 1, column_index(this))
    set_dim_indices!(this, site_int, 1, column_index(this))
    return nothing
end

# =====================================================================================
# Variable registration
# =====================================================================================

"""
    set_restart_var!(this, vname, vtype, long_name, units, flushval, hlms, initialize) -> index

Register one restart variable. Returns the integer index of this variable (0 if
the HLM filter excludes it). When `initialize` is true the variable's payload is
allocated. `this.num_restart_vars_` is advanced as the counter. Mirrors the
Fortran `set_restart_var` (which threads `ivar`/`index` as inout args; here we
use the interface's own counter and return the index).
"""
function set_restart_var!(this::fates_restart_interface_type, vname::AbstractString,
                          vtype::AbstractString, long_name::AbstractString,
                          units::AbstractString, flushval::Real, hlms::AbstractString,
                          initialize::Bool)
    use_var = check_hlm_list(strip(hlms), strip(hlm_name[]))

    if use_var
        this.num_restart_vars_ += 1
        ivar = this.num_restart_vars_
        if initialize
            Init!(this.rvars[ivar], vname, units, long_name, vtype, flushval,
                  fates_restart_num_dim_kinds, this.dim_kinds, this.dim_bounds)
        end
        return ivar
    else
        return 0
    end
end

"""
    RegisterCohortVector!(this, symbol_base, vtype, long_name_base, units, veclength, flushval, hlms, initialize) -> index

Register a vector-valued cohort variable as `veclength` individually-named
scalar restart variables (`<base>_vec_NNN`). Returns the index of the FIRST
position (so the whole vector can be addressed as base, base+1, ...). Mirrors
the Fortran `RegisterCohortVector`.
"""
function RegisterCohortVector!(this::fates_restart_interface_type, symbol_base::AbstractString,
                               vtype::AbstractString, long_name_base::AbstractString,
                               units::AbstractString, veclength::Integer, flushval::Real,
                               hlms::AbstractString, initialize::Bool)
    index = this.num_restart_vars_ + 1
    for i_pos in 1:veclength
        pos_symbol = lpad(string(i_pos), 3, '0')
        symbol = strip(symbol_base) * "_vec_" * pos_symbol
        long_name = strip(long_name_base) * ", position:" * pos_symbol
        set_restart_var!(this, symbol, vtype, long_name, units, flushval, "CLM:ALM", initialize)
    end
    return index
end

"""
    DefineRMeanRestartVar!(this, vname, vtype, long_name, units, initialize) -> index

Register a running-mean variable as a triple (`_cmean`, `_lmean`, `_cindex`).
Returns the index of the first (`_cmean`). Mirrors the Fortran routine.
"""
function DefineRMeanRestartVar!(this::fates_restart_interface_type, vname::AbstractString,
                                vtype::AbstractString, long_name::AbstractString,
                                units::AbstractString, initialize::Bool)
    index = set_restart_var!(this, strip(vname) * "_cmean", vtype,
                             long_name * " current mean", units, flushzero, "CLM:ALM", initialize)
    set_restart_var!(this, strip(vname) * "_lmean", vtype,
                     long_name * " latest mean", units, flushzero, "CLM:ALM", initialize)
    set_restart_var!(this, strip(vname) * "_cindex", vtype,
                     long_name * " index", "index", flushzero, "CLM:ALM", initialize)
    return index
end

"""
    DefinePRTRestartVars!(this, initialize) -> nothing

Register every PARTEH (PRT) state + flux variable. For each PRT state descriptor
and each discrete position we register the instantaneous state (`_val_NNN`), the
turnover flux (`_turn_NNN`), the net-allocation flux (`_net_NNN`) and the burned
mass (`_burned_NNN`). Mirrors the Fortran routine. If `prt_global` is empty
(PARTEH not configured) this registers nothing.
"""
function DefinePRTRestartVars!(this::fates_restart_interface_type, initialize::Bool)
    pg = prt_global[]
    pg === nothing && return nothing
    for i_var in 1:pg.num_vars
        symbol_base = pg.state_descriptor[i_var].symbol
        name_base   = pg.state_descriptor[i_var].longname
        for i_pos in 1:pg.state_descriptor[i_var].num_pos
            pos_symbol = lpad(string(i_pos), 3, '0')
            for (suffix, descr) in (("_val_", ", state var, position:"),
                                    ("_turn_", ", turnover, position:"),
                                    ("_net_", ", net allocation/transp, position:"),
                                    ("_burned_", ", burned mass:"))
                symbol = strip(symbol_base) * suffix * pos_symbol
                long_name = strip(name_base) * descr * pos_symbol
                set_restart_var!(this, symbol, cohort_r8, long_name, "kg", flushzero,
                                 "CLM:ALM", initialize)
            end
        end
    end
    return nothing
end

# =====================================================================================
# define_restart_vars — the FULL registry. Run twice: first to COUNT, then to
# INITIALIZE (allocate payloads). The named ir_* handles are stored in this.ir.
# =====================================================================================

"""
    define_restart_vars!(this, initialize_variables)

Register every FATES restart variable, gated by the same hlm flags as the
Fortran. `initialize_variables=false` only counts; `=true` allocates payloads.
The integer index of every named variable is recorded in `this.ir[:name]`.
"""
function define_restart_vars!(this::fates_restart_interface_type, initialize_variables::Bool)
    this.num_restart_vars_ = 0  # reset counter (Fortran ivar=0)
    ir = this.ir

    S(name, vname, vtype, long, units, flush) =
        (ir[name] = set_restart_var!(this, vname, vtype, long, units, flush, "CLM:ALM", initialize_variables))
    RC(name, base, vtype, long, units, vlen, flush) =
        (ir[name] = RegisterCohortVector!(this, base, vtype, long, units, vlen, flush, "CLM:ALM", initialize_variables))
    RM(name, vname, vtype, long, units) =
        (ir[name] = DefineRMeanRestartVar!(this, vname, vtype, long, units, initialize_variables))

    # ---- Site level variables -------------------------------------------------------
    S(:ir_npatch_si, "fates_PatchesPerSite", site_int, "Total number of FATES patches per column", "none", flushinvalid)
    S(:ir_cd_status_si, "fates_cold_dec_status", site_int, "status flag for cold deciduous plants", "unitless", flushinvalid)
    S(:ir_nchill_days_si, "fates_chilling_days", site_int, "chilling day counter", "unitless", flushinvalid)
    S(:ir_ncold_days_si, "fates_cold_days", site_int, "cold day counter", "unitless", flushinvalid)
    S(:ir_cleafondate_si, "fates_cold_leafondate", site_int, "the model day of last cold leaf on", "absolute integer day", flushinvalid)
    S(:ir_cleafoffdate_si, "fates_cold_leafoffdate", site_int, "the model day last cold leaf off", "absolute integer day", flushinvalid)
    S(:ir_cndaysleafon_si, "fates_cold_ndaysleafon", site_int, "number of days since leaf on (cold deciduous)", "days", flushinvalid)
    S(:ir_cndaysleafoff_si, "fates_cold_ndaysleafoff", site_int, "number of days since leaf off (cold deciduous)", "days", flushinvalid)
    S(:ir_phenmodeldate_si, "fates_phen_model_date", site_int, "integer model day used for phen timing", "absolute integer day", flushinvalid)
    S(:ir_fireweather_index_si, "fates_acc_nesterov_id", site_r8, "a nesterov index accumulator", "unitless", flushzero)
    S(:ir_gdd_si, "fates_gdd_site", site_r8, "growing degree days at each site", "degC days", flushzero)
    S(:ir_min_allowed_landuse_fraction_si, "fates_min_allowed_landuse_fraction_site", site_r8, "minimum allowed land use fraction at each site", "fraction", flushzero)
    S(:ir_landuse_vector_gt_min_si, "fates_landuse_vector_gt_min_site", cohort_int, "minimum allowed land use fraction at each site", "logical", flushzero)
    S(:ir_area_bareground_si, "fates_area_bareground_site", site_r8, "minimum allowed land use fraction at each site", "fraction", flushzero)
    S(:ir_snow_depth_si, "fates_snow_depth_site", site_r8, "average snow depth", "m", flushzero)
    S(:ir_trunk_product_si, "fates_trunk_product_site", site_r8, "Accumulate trunk product flux at site", "kgC/m2", flushzero)
    S(:ir_landuse_config_si, "fates_landuse_config_site", site_int, "hlm_use_potentialveg status of run that created this restart file", "kgC/m2", flushzero)

    # ---- Variables stored within cohort vectors (some are patch/site collapsed) ------
    S(:ir_ncohort_pa, "fates_CohortsPerPatch", cohort_int, "the number of cohorts per patch", "unitless", flushinvalid)
    S(:ir_fcansno_pa, "fates_fcansno_pa", cohort_r8, "Fraction of canopy covered in snow", "unitless", flushinvalid)
    S(:ir_solar_zenith_flag_pa, "fates_solar_zenith_flag_pa", cohort_int, "switch specifying if zenith is positive", "unitless", flushinvalid)
    S(:ir_solar_zenith_angle_pa, "fates_solar_zenith_angle_pa", cohort_r8, "the angle of the solar zenith for each patch", "radians", flushinvalid)

    # ---- 1D cohort variables --------------------------------------------------------
    S(:ir_seed_prod_co, "fates_seed_prod", cohort_r8, "fates cohort - seed production", "kgC/plant", flushinvalid)
    S(:ir_canopy_layer_co, "fates_canopy_layer", cohort_int, "ed cohort - canopy_layer", "unitless", flushinvalid)
    S(:ir_canopy_layer_yesterday_co, "fates_canopy_layer_yesterday", cohort_r8, "ed cohort - canopy_layer_yesterday", "unitless", flushzero)
    S(:ir_crowndamage_co, "fates_crowndamage", cohort_int, "ed cohort - crowndamage class", "unitless", flushinvalid)
    S(:ir_canopy_trim_co, "fates_canopy_trim", cohort_r8, "ed cohort - canopy_trim", "fraction", flushzero)
    S(:ir_l2fr_co, "fates_l2fr", cohort_r8, "ed cohort - l2fr", "fraction", flushzero)

    if hlm_parteh_mode[] == prt_cnp_flex_allom_hyp
        S(:ir_cx_int_co, "fates_cx_int", cohort_r8, "ed cohort - emacx", "fraction", flushzero)
        S(:ir_emadcxdt_co, "fates_emadcxdt", cohort_r8, "ed cohort - emadcxdt", "fraction", flushzero)
        S(:ir_cx0_co, "fates_cx0", cohort_r8, "ed cohort - cx0", "fraction", flushzero)
        S(:ir_cnplimiter_co, "fates_cnplimiter", cohort_r8, "ed cohort - cnp limiter index", "index", flushzero)
        S(:ir_daily_nh4_uptake_co, "fates_daily_nh4_uptake", cohort_r8, "fates cohort- daily ammonium [NH4] uptake", "kg/plant/day", flushzero)
        S(:ir_daily_no3_uptake_co, "fates_daily_no3_uptake", cohort_r8, "fates cohort- daily ammonium [NO3] uptake", "kg/plant/day", flushzero)
        S(:ir_daily_n_fixation_co, "fates_daily_n_fixation", cohort_r8, "fates cohort- daily N symbiotic fixation", "kg/plant/day", flushzero)
        S(:ir_daily_p_uptake_co, "fates_daily_p_uptake", cohort_r8, "fates cohort- daily phosphorus uptake", "kg/plant/day", flushzero)
        S(:ir_daily_p_demand_co, "fates_daily_p_demand", cohort_r8, "fates cohort- daily phosphorus demand", "kgP/plant/day", flushzero)
        S(:ir_daily_n_demand_co, "fates_daily_n_demand", cohort_r8, "fates cohort- daily nitrogen demand", "kgN/plant/day", flushzero)
    end

    S(:ir_size_class_lasttimestep_co, "fates_size_class_lasttimestep", cohort_int, "ed cohort - size-class last timestep", "index", flushzero)
    S(:ir_dbh_co, "fates_dbh", cohort_r8, "ed cohort - diameter at breast height", "cm", flushzero)
    S(:ir_coage_co, "fates_coage", cohort_r8, "ed cohort - age in days", "days", flushzero)
    S(:ir_height_co, "fates_height", cohort_r8, "ed cohort - plant height", "m", flushzero)
    S(:ir_nplant_co, "fates_nplant", cohort_r8, "ed cohort - number of plants in the cohort", "/patch", flushzero)
    S(:ir_gpp_acc_co, "fates_gpp_acc", cohort_r8, "ed cohort - accumulated gpp over dynamics step", "kgC/indiv", flushzero)
    S(:ir_npp_acc_co, "fates_npp_acc", cohort_r8, "ed cohort - accumulated npp over dynamics step", "kgC/indiv", flushzero)
    S(:ir_resp_acc_co, "fates_resp_acc", cohort_r8, "ed cohort - accumulated respiration over dynamics step", "kgC/indiv", flushzero)
    S(:ir_gpp_acc_hold_co, "fates_gpp_acc_hold", cohort_r8, "ed cohort - current step gpp", "kgC/indiv/year", flushzero)
    S(:ir_npp_acc_hold_co, "fates_npp_acc_hold", cohort_r8, "ed cohort - current step npp", "kgC/indiv/year", flushzero)
    S(:ir_resp_acc_hold_co, "fates_resp_acc_hold", cohort_r8, "ed cohort - current step resp", "kgC/indiv/year", flushzero)
    S(:ir_resp_excess_co, "fates_resp_excess", cohort_r8, "ed cohort - maintenance respiration deficit", "kgC/indiv", flushzero)
    S(:ir_bmort_co, "fates_bmort", cohort_r8, "ed cohort - background mortality rate", "/year", flushzero)
    S(:ir_hmort_co, "fates_hmort", cohort_r8, "ed cohort - hydraulic mortality rate", "/year", flushzero)
    S(:ir_cmort_co, "fates_cmort", cohort_r8, "ed cohort - carbon starvation mortality rate", "/year", flushzero)
    S(:ir_frmort_co, "fates_frmort", cohort_r8, "ed cohort - freezing mortality rate", "/year", flushzero)
    S(:ir_smort_co, "fates_smort", cohort_r8, "ed cohort - senescence mortality rate", "/year", flushzero)
    S(:ir_asmort_co, "fates_asmort", cohort_r8, "ed cohort - age senescence mortality rate", "/year", flushzero)
    S(:ir_dgmort_co, "fates_dgmort", cohort_r8, "ed cohort - damage mortality rate", "/year", flushzero)
    S(:ir_lmort_direct_co, "fates_lmort_direct", cohort_r8, "ed cohort - directly logging mortality rate", "%/event", flushzero)
    S(:ir_lmort_collateral_co, "fates_lmort_collateral", cohort_r8, "ed cohort - collateral mortality rate", "%/event", flushzero)
    S(:ir_lmort_infra_co, "fates_lmort_in", cohort_r8, "ed cohort - mechanical mortality rate", "%/event", flushzero)
    S(:ir_ddbhdt_co, "fates_ddbhdt", cohort_r8, "ed cohort - differential: ddbh/dt", "cm/year", flushzero)
    S(:ir_resp_tstep_co, "fates_resp_tstep", cohort_r8, "ed cohort - autotrophic respiration over timestep", "kgC/indiv/timestep", flushzero)
    S(:ir_pft_co, "fates_pft", cohort_int, "ed cohort - plant functional type", "index", flushzero)
    S(:ir_status_co, "fates_status_coh", cohort_int, "ed cohort - plant phenology status", "unitless", flushzero)
    S(:ir_efleaf_co, "fates_efleaf_coh", cohort_r8, "ed cohort - leaf elongation factor", "unitless", flushzero)
    S(:ir_effnrt_co, "fates_effnrt_coh", cohort_r8, "ed cohort - fine-root elongation factor", "unitless", flushzero)
    S(:ir_efstem_co, "fates_efstem_coh", cohort_r8, "ed cohort - stem elongation factor", "unitless", flushzero)
    S(:ir_isnew_co, "fates_isnew", cohort_int, "ed cohort - binary flag specifying if a plant has experienced a full day cycle", "0/1", flushone)
    S(:ir_g_sb_laweight_co, "fates_gsblaweight", cohort_r8, "ed cohort - leaf-area weighted total stomatal+blayer conductance", "[m/s]*[m2]", flushzero)

    # ---- Mixed-dimension variables using the cohort vector --------------------------
    S(:ir_gnd_alb_dif_pasb, "fates_gnd_alb_dif", cohort_r8, "ground albedo of diffuse radiation vis and ir", "fraction", flushzero)
    S(:ir_gnd_alb_dir_pasb, "fates_gnd_alb_dir", cohort_r8, "ground albedo of direct radiation vis and ir", "fraction", flushzero)
    S(:ir_spread_si, "fates_spread", site_r8, "dynamic ratio of dbh to canopy area, by patch x canopy-layer", "cm/m2", flushzero)
    S(:ir_livegrass_pa, "fates_livegrass", cohort_r8, "total AGB from grass, by patch", "kgC/m2", flushzero)
    S(:ir_age_pa, "fates_age", cohort_r8, "age of the ED patch", "yr", flushzero)
    S(:ir_agesinceanthrodist_pa, "fates_age_since_anthro_dist", cohort_r8, "age of the ED patch since last anthropogenic disturbance", "yr", flushzero)
    S(:ir_patchdistturbcat_pa, "fates_patchdistturbcat", cohort_int, "Disturbance label of patch", "yr", flushzero)
    S(:ir_nocomp_pft_label_pa, "fates_nocomp_pft_label", cohort_int, "PFT label of patch in nocomp mode", "none", flushzero)
    S(:ir_area_pa, "fates_area", cohort_r8, "are of the ED patch", "m2", flushzero)

    if hlm_use_sp[] == ifalse
        S(:ir_scorch_ht_pa_pft, "fates_scorch_ht_pa_pft", cohort_r8, "scorch height", "m", flushzero)
        S(:ir_litter_moisture_pa_nfsc, "fates_litter_moisture_pa_nfsc", cohort_r8, "scorch height", "m", flushzero)
    end

    RC(:ir_year_net_up_co, "fates_year_net_up", cohort_r8, "yearly net uptake at leaf layers", "kg/m2/year", nlevleaf, flushzero)

    # ---- Patch-level litter pools (multi-element); site flux + mass-balance ----------
    if hlm_use_sp[] == ifalse
        RC(:ir_agcwd_litt, "fates_ag_cwd", cohort_r8, "above ground CWD", "kg/m2", num_elements[], flushzero)
        RC(:ir_bgcwd_litt, "fates_bg_cwd", cohort_r8, "below ground CWD", "kg/m2", num_elements[], flushzero)
        RC(:ir_leaf_litt, "fates_leaf_fines", cohort_r8, "above ground leaf litter", "kg/m2", num_elements[], flushzero)
        RC(:ir_fnrt_litt, "fates_fnrt_fines", cohort_r8, "fine root litter", "kg/m2", num_elements[], flushzero)
        RC(:ir_seed_litt, "fates_seed", cohort_r8, "seed bank (non-germinated)", "kg/m2", num_elements[], flushzero)
        RC(:ir_seedgerm_litt, "fates_seedgerm", cohort_r8, "seed bank (germinated)", "kg/m2", num_elements[], flushzero)
        RC(:ir_seed_decay_litt, "fates_seed_frag", cohort_r8, "seed bank fragmentation flux (non-germinated)", "kg/m2", num_elements[], flushzero)
        RC(:ir_seedgerm_decay_litt, "fates_seedgerm_frag", cohort_r8, "seed bank fragmentation flux (germinated)", "kg/m2", num_elements[], flushzero)
        RC(:ir_agcwd_frag_litt, "fates_ag_cwd_frag", cohort_r8, "above ground CWD frag flux", "kg/m2/day", num_elements[], flushzero)
        RC(:ir_bgcwd_frag_litt, "fates_bg_cwd_frag", cohort_r8, "below ground CWD frag flux", "kg/m2/day", num_elements[], flushzero)
        RC(:ir_lfines_frag_litt, "fates_lfines_frag", cohort_r8, "frag flux from leaf fines", "kg/m2/day", num_elements[], flushzero)
        RC(:ir_rfines_frag_litt, "fates_rfines_frag", cohort_r8, "frag flux from froot fines", "kg/m2/day", num_elements[], flushzero)
    end

    if hlm_use_sp[] == ifalse
        RC(:ir_cwdagin_flxdg, "fates_cwdagin", cohort_r8, "Input flux of AG CWD", "kg/ha", num_elements[], flushzero)
        RC(:ir_cwdbgin_flxdg, "fates_cwdbgin", cohort_r8, "Input flux of BG CWD", "kg/ha", num_elements[], flushzero)
        RC(:ir_leaflittin_flxdg, "fates_leaflittin", cohort_r8, "Input flux of leaf litter", "kg/ha", num_elements[], flushzero)
        RC(:ir_rootlittin_flxdg, "fates_rootlittin", cohort_r8, "Input flux of root litter", "kg/ha", num_elements[], flushzero)
        RC(:ir_oldstock_mbal, "fates_oldstock", site_r8, "Previous total mass of all fates state variables", "kg/ha", num_elements[], flushzero)
        RC(:ir_errfates_mbal, "fates_errfates", site_r8, "Previous total mass of error fates state variables", "kg/ha", num_elements[], flushzero)
        RC(:ir_liveveg_intflux_el, "fates_liveveg_intflux", site_r8, "total mass of live vegetation of each chemical species, integrated from fluxes", "kg/m2", num_elements[], flushzero)
        RC(:ir_liveveg_err_el, "fates_liveveg_err", site_r8, "total mass error of live vegetation of each chemical species, from integrated fluxes", "kg/m2", num_elements[], flushzero)
        RC(:ir_litter_intflux_el, "fates_litter_intflux", site_r8, "total mass of litter of each chemical species, integrated from fluxes", "kg/m2", num_elements[], flushzero)
        RC(:ir_litter_err_el, "fates_litter_err", site_r8, "total mass error of litter of each chemical species, from integrated fluxes", "kg/m2", num_elements[], flushzero)
    end

    RC(:ir_woodprod_harvest_mbal, "fates_woodprod_harv", cohort_r8, "Current wood product flux from harvest", "kg/m2/day", num_elements[], flushzero)
    RC(:ir_woodprod_landusechange_mbal, "fates_woodprod_luc", cohort_r8, "Current wood product flux from land use change", "kg/m2/day", num_elements[], flushzero)

    S(:ir_c_area_co, "fates_cohort_area", cohort_r8, "area of the fates cohort", "m2", flushzero)
    S(:ir_treelai_co, "fates_cohort_treelai", cohort_r8, "leaf area index of fates cohort", "m2/m2", flushzero)
    S(:ir_treesai_co, "fates_cohort_treesai", cohort_r8, "stem area index of fates cohort", "m2/m2", flushzero)

    if hlm_use_sp[] == itrue
        S(:ir_canopy_layer_tlai_pa, "fates_canopy_layer_tlai_pa", cohort_r8, "total patch level leaf area index of each fates canopy layer", "m2/m2", flushzero)
    end

    S(:ir_nclp_pa, "fates_nclp_pa", cohort_int, "total number of canopy layers", "-", flushzero)
    S(:ir_zstar_pa, "fates_zstar_pa", cohort_r8, "patch zstar", "-", flushzero)

    # ---- Plant hydraulics (if on) ---------------------------------------------------
    if hlm_use_planthydro[] == itrue
        if fates_maxElementsPerSite[] < (nshell * nlevsoi_hyd_max)
            _ri_log_println(" Fates plant hydraulics needs space to store site-level hydraulics info.")
            fates_endrun("define_restart_vars!: fates_maxElementsPerSite < nshell*nlevsoi_hyd_max")
        end
        RC(:ir_hydro_th_ag_covec, "fates_hydro_th_ag", cohort_r8, "water in aboveground compartments", "kg/plant", n_hypool_ag, flushzero)
        RC(:ir_hydro_th_troot, "fates_hydro_th_troot", cohort_r8, "water in transporting roots", "kg/plant", n_hypool_troot, flushzero)
        RC(:ir_hydro_th_aroot_covec, "fates_hydro_th_aroot", cohort_r8, "water in absorbing roots", "kg/plant", nlevsoi_hyd_max, flushzero)
        S(:ir_hydro_liqvol_shell_si, "fates_hydro_liqvol_shell", cohort_r8, "Volumetric water content of rhizosphere compartments (layerxshell)", "m3/m3", flushzero)
        S(:ir_hydro_recruit_si, "fates_hydro_recruit_h2o", site_r8, "Site level water mass used for new recruits", "kg", flushzero)
        S(:ir_hydro_dead_si, "fates_hydro_dead_h2o", site_r8, "Site level water bound in dead plants", "kg", flushzero)
        S(:ir_hydro_growturn_err_si, "fates_hydro_growturn_err", site_r8, "Site level error for hydraulics due to growth/turnover", "kg", flushzero)
        S(:ir_hydro_hydro_err_si, "fates_hydro_hydro_err", site_r8, "Site level error for hydrodynamics", "kg", flushzero)
        S(:ir_hydro_errh2o, "fates_errh2o", cohort_r8, "ed cohort - running plant h2o error for hydro", "kg/indiv", flushzero)
    end

    # ---- site x time level vars -----------------------------------------------------
    S(:ir_dd_status_sift, "fates_drought_dec_status", cohort_int, "status flag for drought deciduous plants", "unitless", flushinvalid)
    S(:ir_dleafondate_sift, "fates_drought_leafondate", cohort_int, "the day of year for drought based leaf-on", "day of year", flushinvalid)
    S(:ir_dleafoffdate_sift, "fates_drought_leafoffdate", cohort_int, "the day of year for drought based leaf-off", "day of year", flushinvalid)
    S(:ir_dndaysleafon_sift, "fates_drought_ndaysleafon", cohort_int, "number of days since leaf on (drought deciduous)", "days", flushinvalid)
    S(:ir_dndaysleafoff_sift, "fates_drought_ndaysleafoff", cohort_int, "number of days since leaf off (drought deciduous)", "days", flushinvalid)
    S(:ir_elong_factor_sift, "fates_elong_factor", cohort_r8, "leaf elongation factor (0 - completely abscissed; 1 - completely flushed)", "unitless", flushinvalid)
    S(:ir_recl2fr_sipfcl, "fates_recruit_l2fr", cohort_r8, "site-level mean recruit l2frs, by pft x canopy layer", "-", flushzero)
    S(:ir_liqvolmem_siwmft, "fates_liqvol_memory", cohort_r8, "last 10 days of volumetric soil water, by site x day-index", "m3/m3", flushzero)
    S(:ir_smpmem_siwmft, "fates_smp_memory", cohort_r8, "last 10 days of soil matric potential, by site x day-index", "mm", flushzero)
    S(:ir_vegtempmem_sitm, "fates_vegtemp_memory", cohort_r8, "last 10 days of 24-hour vegetation temperature, by site x day-index", "m3/m3", flushzero)
    S(:ir_recrate_sift, "fates_recrate", cohort_r8, "fates diagnostics on recruitment", "indiv/ha/day", flushzero)
    S(:ir_use_this_pft_sift, "fates_use_this_pft", cohort_int, "in fixed biogeog mode, is pft in gridcell?", "0/1", flushone)
    S(:ir_area_pft_sift, "fates_area_pft", cohort_r8, "in fixed biogeog mode, what is pft area in gridcell?", "0/1", flushzero)
    S(:ir_seed_in_sift, "fates_seed_in_site", cohort_r8, "Site-level seed mass input from neighboring gridcells per pft", "kg", flushzero)
    S(:ir_seed_out_sift, "fates_seed_out_site", cohort_r8, "Site-level seed mass output to neighboring gridcells per pft", "kg", flushzero)
    S(:ir_fmortrate_cano_siscpf, "fates_fmortrate_canopy", cohort_r8, "fates diagnostics on fire mortality canopy", "indiv/ha/year", flushzero)
    S(:ir_fmortrate_usto_siscpf, "fates_fmortrate_ustory", cohort_r8, "fates diagnostics on fire mortality ustory", "indiv/ha/year", flushzero)
    S(:ir_imortrate_siscpf, "fates_imortrate", cohort_r8, "fates diagnostics on impact mortality", "indiv/ha/year", flushzero)
    S(:ir_fmortrate_crown_siscpf, "fates_fmortrate_crown", cohort_r8, "fates diagnostics on crown fire mortality", "indiv/ha/year", flushzero)
    S(:ir_fmortrate_cambi_siscpf, "fates_fmortrate_cambi", cohort_r8, "fates diagnostics on fire cambial mortality", "indiv/ha/year", flushzero)
    S(:ir_termnindiv_cano_siscpf, "fates_termn_canopy", cohort_r8, "fates diagnostics on termin mortality canopy", "indiv/ha/day", flushzero)
    S(:ir_termnindiv_usto_siscpf, "fates_termn_ustory", cohort_r8, "fates diagnostics on term mortality ustory", "indiv/ha/day", flushzero)
    S(:ir_growflx_fusion_siscpf, "fates_growflx_fusion", cohort_r8, "fates diag: rate of indivs moving via fusion", "indiv/ha/day", flushzero)
    S(:ir_demorate_sisc, "fates_demorate", cohort_r8, "fates diagnoatic rate of indivs demoted", "indiv/ha/day", flushzero)
    S(:ir_promrate_sisc, "fates_promrate", cohort_r8, "fates diagnostic rate of indivs promoted", "indiv/ha/da", flushzero)
    S(:ir_imortcflux_sipft, "fates_imortcflux", cohort_r8, "biomass of indivs killed due to impact mort", "kgC/ha/day", flushzero)
    S(:ir_imortcarea_si, "fates_imortcarea", site_r8, "crownarea of indivs killed due to impact mort", "m2/ha/day", flushzero)
    S(:ir_fmortcflux_cano_sipft, "fates_fmortcflux_canopy", cohort_r8, "fates diagnostic biomass of canopy fire", "gC/m2/sec", flushzero)
    S(:ir_fmortcflux_usto_sipft, "fates_fmortcflux_ustory", cohort_r8, "fates diagnostic biomass of understory fire", "gC/m2/sec", flushzero)
    S(:ir_termcflux_cano_sipft, "fates_termcflux_canopy", cohort_r8, "fates diagnostic term carbon flux canopy", "", flushzero)
    S(:ir_termcflux_usto_sipft, "fates_termcflux_ustory", cohort_r8, "fates diagnostic term carbon flux understory", "", flushzero)
    S(:ir_abg_term_flux_siscpf, "fates_abg_term_flux", cohort_r8, "fates aboveground biomass loss from termination mortality", "", flushzero)
    S(:ir_abg_imort_flux_siscpf, "fates_abg_imort_flux", cohort_r8, "fates aboveground biomass loss from impact mortality", "", flushzero)
    S(:ir_abg_fmort_flux_siscpf, "fates_abg_fmort_flux", cohort_r8, "fates aboveground biomass loss from fire mortality", "", flushzero)
    S(:ir_democflux_si, "fates_democflux", site_r8, "fates diagnostic demotion carbon flux", "", flushzero)
    S(:ir_promcflux_si, "fates_promcflux", site_r8, "fates diagnostic promotion carbon flux ", "", flushzero)
    S(:ir_fmortcarea_cano_si, "fates_fmortcarea_canopy", site_r8, "fates diagnostic crownarea of canopy fire", "m2/sec", flushzero)
    S(:ir_fmortcarea_usto_si, "fates_fmortcarea_ustory", site_r8, "fates diagnostic crownarea of understory fire", "m2/sec", flushzero)
    S(:ir_termcarea_cano_si, "fates_termcarea_canopy", site_r8, "fates diagnostic term crownarea canopy", "", flushzero)
    S(:ir_termcarea_usto_si, "fates_termcarea_ustory", site_r8, "fates diagnostic term crownarea understory", "", flushzero)

    # ---- Damage variables -----------------------------------------------------------
    S(:ir_imortrate_sicdpf, "fates_imortrate_dam", cohort_r8, "fates diagnostics on impact mortality by damage class", "indiv/ha/year", flushzero)
    S(:ir_termnindiv_cano_sicdpf, "fates_termn_cano_dam", cohort_r8, "fates diagnostics on termination mortality by damage class -canopy", "indiv/ha/year", flushzero)
    S(:ir_termnindiv_usto_sicdpf, "fates_termn_usto_dam", cohort_r8, "fates diagnostics on termination mortality by damage class -understory", "indiv/ha/year", flushzero)
    S(:ir_fmortrate_cano_sicdpf, "fates_fmortrate_cano_dam", cohort_r8, "fates diagnostics on fire mortality by damage class", "indiv/ha/year", flushzero)
    S(:ir_fmortrate_usto_sicdpf, "fates_fmortrate_usto_dam", cohort_r8, "fates diagnostics on fire mortality by damage class", "indiv/ha/year", flushzero)
    S(:ir_imortcflux_sicdsc, "fates_imortcflux_dam", cohort_r8, "biomass of indivs killed due to impact mort by damage class", "kgC/ha/day", flushzero)
    S(:ir_termcflux_cano_sicdsc, "fates_termcflux_cano_dam", cohort_r8, "biomass of indivs killed due to termination  mort by damage class", "kgC/ha/day", flushzero)
    S(:ir_termcflux_usto_sicdsc, "fates_termcflux_usto_dam", cohort_r8, "biomass of indivs killed due to termination  mort by damage class", "kgC/ha/day", flushzero)
    S(:ir_fmortcflux_cano_sicdsc, "fates_fmortcflux_cano_dam", cohort_r8, "biomass of indivs killed due to fire mort by damage class", "kgC/ha/day", flushzero)
    S(:ir_fmortcflux_usto_sicdsc, "fates_fmortcflux_usto_dam", cohort_r8, "biomass of indivs killed due to fire mort by damage class", "kgC/ha/day", flushzero)
    S(:ir_crownarea_cano_si, "fates_crownarea_canopy_damage", site_r8, "fates area lost from damage each year", "m2/ha/year", flushzero)
    S(:ir_crownarea_usto_si, "fates_crownarea_understory_damage", site_r8, "fates area lost from damage each year", "m2/ha/year", flushzero)
    S(:ir_emanpp_si, "fates_emanpp", site_r8, "smoothed NPP (exp. moving avg) at the site level (for fixation)", "kg/m2/yr", flushzero)

    # ---- Running means --------------------------------------------------------------
    RM(:ir_tveg24_pa, "fates_tveg24patch", cohort_r8, "24-hour patch veg temp", "K")
    RM(:ir_disturbance_rates_siluludi, "fates_disturbance_rates", cohort_r8, "disturbance rates by donor land-use type, receiver land-use type, and disturbance type", "1/day")

    if ed_params().regeneration_model == TRS_regeneration
        RM(:ir_seedling_layer_par24_pa, "fates_seedling_layer_par24", cohort_r8, "24-hour seedling layer PAR", "W m2-1")
        RM(:ir_sdlng_emerg_smp_pa, "fates_sdlng_emerg_smp", cohort_r8, "seedling layer PAR on the seedling emergence timescale", "mm suction")
        RM(:ir_sdlng_mort_par_pa, "fates_sdlng_mort_par", cohort_r8, "seedling layer PAR on the seedling mortality timescale", "W m2-1")
        RM(:ir_sdlng2sap_par_pa, "fates_sdlng2sap_par", cohort_r8, "seedling layer PAR on the seedling to sapling transition timescale", "W m2-1")
        RM(:ir_sdlng_mdd_pa, "fates_sdlng_mdd", cohort_r8, "seedling moisture deficit days", "mm days")
    end

    RM(:ir_tveglpa_pa, "fates_tveglpapatch", cohort_r8, "running average (EMA) of patch veg temp for photo acclim", "K")
    RM(:ir_tveglongterm_pa, "fates_tveglongtermpatch", cohort_r8, "long-term (T_home) running average (EMA) of patch veg temp for photo acclim", "K")

    # ---- PRT states + fluxes --------------------------------------------------------
    this.ir_prt_base = this.num_restart_vars_
    DefinePRTRestartVars!(this, initialize_variables)

    return nothing
end

# =====================================================================================

"""
    initialize_restart_vars!(this)

Run `define_restart_vars!` twice — once to count (allocating the `rvars` vector),
once to initialize each variable's payload. Mirrors the Fortran.
"""
function initialize_restart_vars!(this::fates_restart_interface_type)
    define_restart_vars!(this, false)
    this.rvars = [fates_restart_variable_type() for _ in 1:num_restart_vars(this)]
    define_restart_vars!(this, true)
    return nothing
end

"""
    flush_rvars!(this, nc)

Flush every registered restart variable to its flushval (over thread nc).
"""
function flush_rvars!(this::fates_restart_interface_type, nc::Integer)
    for ivar in eachindex(this.rvars)
        Flush!(this.rvars[ivar], nc, this.dim_bounds, this.dim_kinds)
    end
    return nothing
end

# =====================================================================================
# Running-mean / cohort-vector pack-unpack helpers
# =====================================================================================

function GetCohortRealVector(this::fates_restart_interface_type, state_vector::AbstractVector{Float64},
                             len_state_vector::Integer, variable_index_base::Integer,
                             co_global_index::Integer)
    ir_pos_var = variable_index_base
    @inbounds for i_pos in 1:len_state_vector
        state_vector[i_pos] = this.rvars[ir_pos_var].r81d[co_global_index]
        ir_pos_var += 1
    end
    return nothing
end

function SetCohortRealVector(this::fates_restart_interface_type, state_vector::AbstractVector{Float64},
                             len_state_vector::Integer, variable_index_base::Integer,
                             co_global_index::Integer)
    ir_pos_var = variable_index_base
    @inbounds for i_pos in 1:len_state_vector
        this.rvars[ir_pos_var].r81d[co_global_index] = state_vector[i_pos]
        ir_pos_var += 1
    end
    return nothing
end

# Running-mean restart accessors. A running-mean variable is registered as a
# TRIPLE (`_cmean`, `_lmean`, `_cindex`) starting at `ir_var_index` (see
# `DefineRMeanRestartVar!`); these read/write the three payload slots at one
# position. Mirrors the Fortran `GetRMeanRestartVar` / `SetRMeanRestartVar`.
function GetRMeanRestartVar(this::fates_restart_interface_type, rmean_var::rmean_type,
                            ir_var_index::Integer, position_index::Integer)
    rmean_var.c_mean  = this.rvars[ir_var_index].r81d[position_index]
    rmean_var.l_mean  = this.rvars[ir_var_index + 1].r81d[position_index]
    rmean_var.c_index = round(Int, this.rvars[ir_var_index + 2].r81d[position_index])
    return nothing
end

function SetRMeanRestartVar(this::fates_restart_interface_type, rmean_var::rmean_type,
                            ir_var_index::Integer, position_index::Integer)
    this.rvars[ir_var_index].r81d[position_index]     = rmean_var.c_mean
    this.rvars[ir_var_index + 1].r81d[position_index] = rmean_var.l_mean
    this.rvars[ir_var_index + 2].r81d[position_index] = Float64(rmean_var.c_index)
    return nothing
end

# Real / int payload accessors (1-based io_idx already)
@inline _setr!(this, irname, idx, val) = (this.rvars[this.ir[irname]].r81d[idx] = Float64(val); nothing)
@inline _seti!(this, irname, idx, val) = (this.rvars[this.ir[irname]].int1d[idx] = round(Int, val); nothing)
@inline _getr(this, irname, idx) = this.rvars[this.ir[irname]].r81d[idx]
@inline _geti(this, irname, idx) = this.rvars[this.ir[irname]].int1d[idx]

# Vector-valued / per-element accessors: address rvar `this.ir[irname] + off`
# (the RegisterCohortVector / per-element handle base + offset; Fortran `ir_*+el`).
@inline _setr_off!(this, irname, off, idx, val) =
    (this.rvars[this.ir[irname] + off].r81d[idx] = Float64(val); nothing)
@inline _getr_off(this, irname, off, idx) = this.rvars[this.ir[irname] + off].r81d[idx]

# =====================================================================================
# R7: cohort PRT pool pack & unpack.
#
# The PRT (PARTEH) pools are registered (in `DefinePRTRestartVars!`) as a flat run
# of restart vars STARTING just after `this.ir_prt_base`, in the order
# (val, turnover, net_alloc, burned) per discrete position, per PRT variable —
# exactly the Fortran `ir_prt_var = ir_prt_base; do i_var; do i_pos; +1 each`.
# These helpers walk that same order so pack/unpack address the identical slots.
# =====================================================================================

"""
    _pack_cohort_prt!(this, ccohort, io_idx_co)

Pack one cohort's PRT pools (instantaneous `val` + the `turnover`/`net_alloc`/
`burned` fluxes, per discrete position, per PRT variable) into the flat restart
vectors at cohort slot `io_idx_co`. Mirrors the Fortran PRT block of
`set_restart_vectors`.
"""
function _pack_cohort_prt!(this::fates_restart_interface_type, ccohort, io_idx_co::Integer)
    pg = prt_global[]
    pg === nothing && return nothing
    ir_prt_var = this.ir_prt_base
    @inbounds for i_var in 1:pg.num_vars
        v = ccohort.prt.variables[i_var]
        for i_pos in 1:pg.state_descriptor[i_var].num_pos
            ir_prt_var += 1; this.rvars[ir_prt_var].r81d[io_idx_co] = v.val[i_pos]
            ir_prt_var += 1; this.rvars[ir_prt_var].r81d[io_idx_co] = v.turnover[i_pos]
            ir_prt_var += 1; this.rvars[ir_prt_var].r81d[io_idx_co] = v.net_alloc[i_pos]
            ir_prt_var += 1; this.rvars[ir_prt_var].r81d[io_idx_co] = v.burned[i_pos]
        end
    end
    return nothing
end

"""
    _unpack_cohort_prt!(this, ccohort, io_idx_co)

Read one cohort's PRT pools back into `ccohort.prt` (symmetric to
[`_pack_cohort_prt!`](@ref)). The cohort's `prt` object must already be allocated
(via `InitPRTObject!`) before this is called.
"""
function _unpack_cohort_prt!(this::fates_restart_interface_type, ccohort, io_idx_co::Integer)
    pg = prt_global[]
    pg === nothing && return nothing
    ir_prt_var = this.ir_prt_base
    @inbounds for i_var in 1:pg.num_vars
        v = ccohort.prt.variables[i_var]
        for i_pos in 1:pg.state_descriptor[i_var].num_pos
            ir_prt_var += 1; v.val[i_pos]       = this.rvars[ir_prt_var].r81d[io_idx_co]
            ir_prt_var += 1; v.turnover[i_pos]  = this.rvars[ir_prt_var].r81d[io_idx_co]
            ir_prt_var += 1; v.net_alloc[i_pos] = this.rvars[ir_prt_var].r81d[io_idx_co]
            ir_prt_var += 1; v.burned[i_pos]    = this.rvars[ir_prt_var].r81d[io_idx_co]
        end
    end
    return nothing
end

# =====================================================================================
# R5: site demographic / mortality / damage / flux-diag array pack & unpack.
#
# All of these site-level arrays are stored in the COHORT-vector address space,
# starting at the site's first-cohort base address (`io_idx_co_1st_init`, captured
# BEFORE the per-patch walk advances it), each with its OWN independent cursor.
# The pure-scalar site fields are stored at `io_idx_si`. Mirrors the Fortran
# `set_restart_vectors` / `get_restart_vectors` site loop exactly (the running-mean
# `disturbance_rates` block — R6 — is also packed/unpacked here).
# =====================================================================================

"""
    _pack_site_diagnostics!(this, site, io_idx_si, io_idx_base)

Pack the R5 site-level demographic / mortality / damage / flux-diagnostic arrays
for one site. `io_idx_base` is the site's first-cohort base address.
"""
function _pack_site_diagnostics!(this::fates_restart_interface_type, site,
                                 io_idx_si::Integer, io_idx_base::Integer)
    nsc  = nlevsclass[]
    npft = numpft[]

    # --- recruitment rate / use_this_pft (per pft, at base offset) ---
    for i_pft in 1:npft
        _setr!(this, :ir_recrate_sift, io_idx_base + i_pft - 1, site.recruitment_rate[i_pft])
        _seti!(this, :ir_use_this_pft_sift, io_idx_base + i_pft - 1, site.use_this_pft[i_pft])
    end

    # --- area_PFT (per pft x land-use, packed pft-major) ---
    for i_pft in 1:npft
        for i_landuse in 1:n_landuse_cats
            i_pflu = i_landuse + (i_pft - 1) * n_landuse_cats
            _setr!(this, :ir_area_pft_sift, io_idx_base + i_pflu - 1, site.area_PFT[i_pft, i_landuse])
        end
    end

    _setr!(this, :ir_min_allowed_landuse_fraction_si, io_idx_si, site.min_allowed_landuse_fraction)
    io_idx_si_lu = io_idx_base
    for i_landuse in 1:n_landuse_cats
        _seti!(this, :ir_landuse_vector_gt_min_si, io_idx_si_lu,
               site.landuse_vector_gt_min[i_landuse] ? itrue : ifalse)
        io_idx_si_lu += 1
    end
    _setr!(this, :ir_area_bareground_si, io_idx_si, site.area_bareground)

    # --- size-class x pft mortality / flux block ---
    io_idx_si_scpf = io_idx_base
    io_idx_si_scpf_term = io_idx_base
    for i_scls in 1:nsc
        for i_pft in 1:npft
            _setr!(this, :ir_fmortrate_cano_siscpf, io_idx_si_scpf, site.fmort_rate_canopy[i_scls, i_pft])
            _setr!(this, :ir_fmortrate_usto_siscpf, io_idx_si_scpf, site.fmort_rate_ustory[i_scls, i_pft])
            _setr!(this, :ir_imortrate_siscpf, io_idx_si_scpf, site.imort_rate[i_scls, i_pft])
            _setr!(this, :ir_fmortrate_crown_siscpf, io_idx_si_scpf, site.fmort_rate_crown[i_scls, i_pft])
            _setr!(this, :ir_fmortrate_cambi_siscpf, io_idx_si_scpf, site.fmort_rate_cambial[i_scls, i_pft])
            _setr!(this, :ir_growflx_fusion_siscpf, io_idx_si_scpf, site.growthflux_fusion[i_scls, i_pft])
            _setr!(this, :ir_abg_term_flux_siscpf, io_idx_si_scpf, site.term_abg_flux[i_scls, i_pft])
            _setr!(this, :ir_abg_imort_flux_siscpf, io_idx_si_scpf, site.imort_abg_flux[i_scls, i_pft])
            _setr!(this, :ir_abg_fmort_flux_siscpf, io_idx_si_scpf, site.fmort_abg_flux[i_scls, i_pft])
            io_idx_si_scpf += 1
            for i_term in 1:n_term_mort_types
                _setr!(this, :ir_termnindiv_cano_siscpf, io_idx_si_scpf_term, site.term_nindivs_canopy[i_term, i_scls, i_pft])
                _setr!(this, :ir_termnindiv_usto_siscpf, io_idx_si_scpf_term, site.term_nindivs_ustory[i_term, i_scls, i_pft])
                io_idx_si_scpf_term += 1
            end
        end
    end

    # --- per-pft carbon flux + drought phenology block ---
    io_idx_si_pft = io_idx_base
    io_idx_si_pft_term = io_idx_base
    for i_pft in 1:npft
        for i_term in 1:n_term_mort_types
            _setr!(this, :ir_termcflux_cano_sipft, io_idx_si_pft_term, site.term_carbonflux_canopy[i_term, i_pft])
            _setr!(this, :ir_termcflux_usto_sipft, io_idx_si_pft_term, site.term_carbonflux_ustory[i_term, i_pft])
            io_idx_si_pft_term += 1
        end
        _setr!(this, :ir_fmortcflux_cano_sipft, io_idx_si_pft, site.fmort_carbonflux_canopy[i_pft])
        _setr!(this, :ir_fmortcflux_usto_sipft, io_idx_si_pft, site.fmort_carbonflux_ustory[i_pft])
        _setr!(this, :ir_imortcflux_sipft, io_idx_si_pft, site.imort_carbonflux[i_pft])
        _seti!(this, :ir_dd_status_sift, io_idx_si_pft, site.dstatus[i_pft])
        _seti!(this, :ir_dleafondate_sift, io_idx_si_pft, site.dleafondate[i_pft])
        _seti!(this, :ir_dleafoffdate_sift, io_idx_si_pft, site.dleafoffdate[i_pft])
        _seti!(this, :ir_dndaysleafon_sift, io_idx_si_pft, site.dndaysleafon[i_pft])
        _seti!(this, :ir_dndaysleafoff_sift, io_idx_si_pft, site.dndaysleafoff[i_pft])
        _setr!(this, :ir_elong_factor_sift, io_idx_si_pft, site.elong_factor[i_pft])
        _setr!(this, :ir_seed_in_sift, io_idx_si_pft, site.seed_in[i_pft])
        _setr!(this, :ir_seed_out_sift, io_idx_si_pft, site.seed_out[i_pft])
        io_idx_si_pft += 1
    end

    # --- R6: site-level disturbance-rate diagnostic (lu_donor x lu_receiver x
    #     dist_type, packed donor-major). Registered as a running-mean triple but
    #     only the base (`_cmean`) slot carries the value (Fortran rio_ pointer). ---
    io_idx_si_luludi = io_idx_base
    for i_lu_donor in 1:n_landuse_cats
        for i_lu_receiver in 1:n_landuse_cats
            for i_dist in 1:N_DIST_TYPES
                _setr!(this, :ir_disturbance_rates_siluludi, io_idx_si_luludi,
                       site.disturbance_rates[i_dist, i_lu_donor, i_lu_receiver])
                io_idx_si_luludi += 1
            end
        end
    end

    # --- flux diagnostics + mass balance (per element), non-SP ---
    if hlm_use_sp[] == ifalse
        for el in 1:num_elements[]
            io_idx_si_cwd = io_idx_base
            io_idx_si_pft2 = io_idx_base
            elem = site.flux_diags.elem[el]
            mbal = site.mass_balance[el]
            ifb  = site.iflux_balance[el]
            for i_cwd in 1:ncwd
                _setr_off!(this, :ir_cwdagin_flxdg, el - 1, io_idx_si_cwd, elem.cwd_ag_input[i_cwd])
                _setr_off!(this, :ir_cwdbgin_flxdg, el - 1, io_idx_si_cwd, elem.cwd_bg_input[i_cwd])
                io_idx_si_cwd += 1
            end
            for i_pft in 1:npft
                _setr_off!(this, :ir_leaflittin_flxdg, el - 1, io_idx_si_pft2, elem.surf_fine_litter_input[i_pft])
                _setr_off!(this, :ir_rootlittin_flxdg, el - 1, io_idx_si_pft2, elem.root_litter_input[i_pft])
                _setr_off!(this, :ir_woodprod_harvest_mbal, el - 1, io_idx_si_pft2, mbal.wood_product_harvest[i_pft])
                _setr_off!(this, :ir_woodprod_landusechange_mbal, el - 1, io_idx_si_pft2, mbal.wood_product_landusechange[i_pft])
                io_idx_si_pft2 += 1
            end
            _setr_off!(this, :ir_oldstock_mbal, el - 1, io_idx_si, mbal.old_stock)
            _setr_off!(this, :ir_errfates_mbal, el - 1, io_idx_si, mbal.err_fates)
            _setr_off!(this, :ir_liveveg_intflux_el, el - 1, io_idx_si, ifb.iflux_liveveg)
            _setr_off!(this, :ir_liveveg_err_el, el - 1, io_idx_si, elem.err_liveveg)
            _setr_off!(this, :ir_litter_intflux_el, el - 1, io_idx_si, ifb.iflux_litter)
            _setr_off!(this, :ir_litter_err_el, el - 1, io_idx_si, elem.err_litter)
        end
    end

    _setr!(this, :ir_spread_si, io_idx_si, site.spread)

    # --- demotion / promotion rate (per size class) ---
    io_idx_si_sc = io_idx_base
    for i_scls in 1:nsc
        _setr!(this, :ir_demorate_sisc, io_idx_si_sc, site.demotion_rate[i_scls])
        _setr!(this, :ir_promrate_sisc, io_idx_si_sc, site.promotion_rate[i_scls])
        io_idx_si_sc += 1
    end

    _setr!(this, :ir_termcarea_cano_si, io_idx_si, site.term_crownarea_canopy)
    _setr!(this, :ir_termcarea_usto_si, io_idx_si, site.term_crownarea_ustory)
    _setr!(this, :ir_emanpp_si, io_idx_si, site.ema_npp)

    # --- damage cross-tab (gated) ---
    if hlm_use_tree_damage[] == itrue
        io_idx_si_cdpf = io_idx_base
        io_idx_si_cdsc = io_idx_base
        for i_scls in 1:nsc
            for i_cdam in 1:nlevdamage[]
                for i_pft in 1:npft
                    _setr!(this, :ir_imortrate_sicdpf, io_idx_si_cdpf, site.imort_rate_damage[i_cdam, i_scls, i_pft])
                    _setr!(this, :ir_termnindiv_cano_sicdpf, io_idx_si_cdpf, site.term_nindivs_canopy_damage[i_cdam, i_scls, i_pft])
                    _setr!(this, :ir_termnindiv_usto_sicdpf, io_idx_si_cdpf, site.term_nindivs_ustory_damage[i_cdam, i_scls, i_pft])
                    _setr!(this, :ir_imortcflux_sicdsc, io_idx_si_cdsc, site.imort_cflux_damage[i_cdam, i_scls])
                    _setr!(this, :ir_termcflux_cano_sicdsc, io_idx_si_cdsc, site.term_cflux_canopy_damage[i_cdam, i_scls])
                    _setr!(this, :ir_termcflux_usto_sicdsc, io_idx_si_cdsc, site.term_cflux_ustory_damage[i_cdam, i_scls])
                    _setr!(this, :ir_fmortrate_cano_sicdpf, io_idx_si_cdpf, site.fmort_rate_canopy_damage[i_cdam, i_scls, i_pft])
                    _setr!(this, :ir_fmortrate_usto_sicdpf, io_idx_si_cdpf, site.fmort_rate_ustory_damage[i_cdam, i_scls, i_pft])
                    _setr!(this, :ir_fmortcflux_cano_sicdsc, io_idx_si_cdsc, site.fmort_cflux_canopy_damage[i_cdam, i_scls])
                    _setr!(this, :ir_fmortcflux_usto_sicdsc, io_idx_si_cdsc, site.fmort_cflux_ustory_damage[i_cdam, i_scls])
                    io_idx_si_cdsc += 1
                    io_idx_si_cdpf += 1
                end
            end
        end
        _setr!(this, :ir_crownarea_cano_si, io_idx_si, site.crownarea_canopy_damage)
        _setr!(this, :ir_crownarea_usto_si, io_idx_si, site.crownarea_ustory_damage)
    end

    _setr!(this, :ir_democflux_si, io_idx_si, site.demotion_carbonflux)
    _setr!(this, :ir_promcflux_si, io_idx_si, site.promotion_carbonflux)
    _setr!(this, :ir_imortcarea_si, io_idx_si, site.imort_crownarea)
    _setr!(this, :ir_fmortcarea_cano_si, io_idx_si, site.fmort_crownarea_canopy)
    _setr!(this, :ir_fmortcarea_usto_si, io_idx_si, site.fmort_crownarea_ustory)

    # --- recruit l2fr (per pft x canopy layer) + soil-water / veg-temp memory ---
    io_idx_si_pfcl = io_idx_base
    for i in 1:nclmax
        for i_pft in 1:npft
            _setr!(this, :ir_recl2fr_sipfcl, io_idx_si_pfcl, site.rec_l2fr[i_pft, i])
            io_idx_si_pfcl += 1
        end
    end
    io_idx_si_wmem = io_idx_base
    for i in 1:numWaterMem
        for i_pft in 1:npft
            _setr!(this, :ir_liqvolmem_siwmft, io_idx_si_wmem, site.liqvol_memory[i, i_pft])
            _setr!(this, :ir_smpmem_siwmft, io_idx_si_wmem, site.smp_memory[i, i_pft])
            io_idx_si_wmem += 1
        end
    end
    io_idx_si_vtmem = io_idx_base
    for i in 1:num_vegtemp_mem
        _setr!(this, :ir_vegtempmem_sitm, io_idx_si_vtmem, site.vegtemp_memory[i])
        io_idx_si_vtmem += 1
    end

    # --- site-level hydraulics scalars (R2d, gated) ---
    if hlm_use_planthydro[] == itrue && site.si_hydr !== nothing
        sh = site.si_hydr
        _setr!(this, :ir_hydro_recruit_si, io_idx_si, sh.h2oveg_recruit)
        _setr!(this, :ir_hydro_dead_si, io_idx_si, sh.h2oveg_dead)
        _setr!(this, :ir_hydro_growturn_err_si, io_idx_si, sh.h2oveg_growturn_err)
        _setr!(this, :ir_hydro_hydro_err_si, io_idx_si, sh.h2oveg_hydro_err)
        io_idx_si_lyr_shell = io_idx_base
        for i in 1:sh.nlevrhiz
            for k in 1:nshell
                _setr!(this, :ir_hydro_liqvol_shell_si, io_idx_si_lyr_shell, sh.h2osoi_liqvol_shell[i, k])
                io_idx_si_lyr_shell += 1
            end
        end
    end
    return nothing
end

"""
    _unpack_site_diagnostics!(this, site, io_idx_si, io_idx_base)

Read back the R5 site-level arrays (symmetric to [`_pack_site_diagnostics!`](@ref)).
"""
function _unpack_site_diagnostics!(this::fates_restart_interface_type, site,
                                   io_idx_si::Integer, io_idx_base::Integer)
    nsc  = nlevsclass[]
    npft = numpft[]

    for i_pft in 1:npft
        site.recruitment_rate[i_pft] = _getr(this, :ir_recrate_sift, io_idx_base + i_pft - 1)
        site.use_this_pft[i_pft] = _geti(this, :ir_use_this_pft_sift, io_idx_base + i_pft - 1)
    end
    for i_pft in 1:npft
        for i_landuse in 1:n_landuse_cats
            i_pflu = i_landuse + (i_pft - 1) * n_landuse_cats
            site.area_PFT[i_pft, i_landuse] = _getr(this, :ir_area_pft_sift, io_idx_base + i_pflu - 1)
        end
    end
    site.min_allowed_landuse_fraction = _getr(this, :ir_min_allowed_landuse_fraction_si, io_idx_si)
    io_idx_si_lu = io_idx_base
    for i_landuse in 1:n_landuse_cats
        site.landuse_vector_gt_min[i_landuse] = (_geti(this, :ir_landuse_vector_gt_min_si, io_idx_si_lu) == itrue)
        io_idx_si_lu += 1
    end
    site.area_bareground = _getr(this, :ir_area_bareground_si, io_idx_si)

    io_idx_si_scpf = io_idx_base
    io_idx_si_scpf_term = io_idx_base
    for i_scls in 1:nsc
        for i_pft in 1:npft
            site.fmort_rate_canopy[i_scls, i_pft] = _getr(this, :ir_fmortrate_cano_siscpf, io_idx_si_scpf)
            site.fmort_rate_ustory[i_scls, i_pft] = _getr(this, :ir_fmortrate_usto_siscpf, io_idx_si_scpf)
            site.imort_rate[i_scls, i_pft] = _getr(this, :ir_imortrate_siscpf, io_idx_si_scpf)
            site.fmort_rate_crown[i_scls, i_pft] = _getr(this, :ir_fmortrate_crown_siscpf, io_idx_si_scpf)
            site.fmort_rate_cambial[i_scls, i_pft] = _getr(this, :ir_fmortrate_cambi_siscpf, io_idx_si_scpf)
            site.growthflux_fusion[i_scls, i_pft] = _getr(this, :ir_growflx_fusion_siscpf, io_idx_si_scpf)
            site.term_abg_flux[i_scls, i_pft] = _getr(this, :ir_abg_term_flux_siscpf, io_idx_si_scpf)
            site.imort_abg_flux[i_scls, i_pft] = _getr(this, :ir_abg_imort_flux_siscpf, io_idx_si_scpf)
            site.fmort_abg_flux[i_scls, i_pft] = _getr(this, :ir_abg_fmort_flux_siscpf, io_idx_si_scpf)
            io_idx_si_scpf += 1
            for i_term in 1:n_term_mort_types
                site.term_nindivs_canopy[i_term, i_scls, i_pft] = _getr(this, :ir_termnindiv_cano_siscpf, io_idx_si_scpf_term)
                site.term_nindivs_ustory[i_term, i_scls, i_pft] = _getr(this, :ir_termnindiv_usto_siscpf, io_idx_si_scpf_term)
                io_idx_si_scpf_term += 1
            end
        end
    end

    io_idx_si_pft = io_idx_base
    io_idx_si_pft_term = io_idx_base
    for i_pft in 1:npft
        for i_term in 1:n_term_mort_types
            site.term_carbonflux_canopy[i_term, i_pft] = _getr(this, :ir_termcflux_cano_sipft, io_idx_si_pft_term)
            site.term_carbonflux_ustory[i_term, i_pft] = _getr(this, :ir_termcflux_usto_sipft, io_idx_si_pft_term)
            io_idx_si_pft_term += 1
        end
        site.fmort_carbonflux_canopy[i_pft] = _getr(this, :ir_fmortcflux_cano_sipft, io_idx_si_pft)
        site.fmort_carbonflux_ustory[i_pft] = _getr(this, :ir_fmortcflux_usto_sipft, io_idx_si_pft)
        site.imort_carbonflux[i_pft] = _getr(this, :ir_imortcflux_sipft, io_idx_si_pft)
        site.dstatus[i_pft] = _geti(this, :ir_dd_status_sift, io_idx_si_pft)
        site.dleafondate[i_pft] = _geti(this, :ir_dleafondate_sift, io_idx_si_pft)
        site.dleafoffdate[i_pft] = _geti(this, :ir_dleafoffdate_sift, io_idx_si_pft)
        site.dndaysleafon[i_pft] = _geti(this, :ir_dndaysleafon_sift, io_idx_si_pft)
        site.dndaysleafoff[i_pft] = _geti(this, :ir_dndaysleafoff_sift, io_idx_si_pft)
        site.elong_factor[i_pft] = _getr(this, :ir_elong_factor_sift, io_idx_si_pft)
        site.seed_in[i_pft] = _getr(this, :ir_seed_in_sift, io_idx_si_pft)
        site.seed_out[i_pft] = _getr(this, :ir_seed_out_sift, io_idx_si_pft)
        io_idx_si_pft += 1
    end

    # --- R6: site-level disturbance-rate diagnostic (symmetric to pack) ---
    io_idx_si_luludi = io_idx_base
    for i_lu_donor in 1:n_landuse_cats
        for i_lu_receiver in 1:n_landuse_cats
            for i_dist in 1:N_DIST_TYPES
                site.disturbance_rates[i_dist, i_lu_donor, i_lu_receiver] =
                    _getr(this, :ir_disturbance_rates_siluludi, io_idx_si_luludi)
                io_idx_si_luludi += 1
            end
        end
    end

    if hlm_use_sp[] == ifalse
        for el in 1:num_elements[]
            io_idx_si_cwd = io_idx_base
            io_idx_si_pft2 = io_idx_base
            elem = site.flux_diags.elem[el]
            mbal = site.mass_balance[el]
            ifb  = site.iflux_balance[el]
            for i_cwd in 1:ncwd
                elem.cwd_ag_input[i_cwd] = _getr_off(this, :ir_cwdagin_flxdg, el - 1, io_idx_si_cwd)
                elem.cwd_bg_input[i_cwd] = _getr_off(this, :ir_cwdbgin_flxdg, el - 1, io_idx_si_cwd)
                io_idx_si_cwd += 1
            end
            for i_pft in 1:npft
                elem.surf_fine_litter_input[i_pft] = _getr_off(this, :ir_leaflittin_flxdg, el - 1, io_idx_si_pft2)
                elem.root_litter_input[i_pft] = _getr_off(this, :ir_rootlittin_flxdg, el - 1, io_idx_si_pft2)
                mbal.wood_product_harvest[i_pft] = _getr_off(this, :ir_woodprod_harvest_mbal, el - 1, io_idx_si_pft2)
                mbal.wood_product_landusechange[i_pft] = _getr_off(this, :ir_woodprod_landusechange_mbal, el - 1, io_idx_si_pft2)
                io_idx_si_pft2 += 1
            end
            mbal.old_stock = _getr_off(this, :ir_oldstock_mbal, el - 1, io_idx_si)
            mbal.err_fates = _getr_off(this, :ir_errfates_mbal, el - 1, io_idx_si)
            ifb.iflux_liveveg = _getr_off(this, :ir_liveveg_intflux_el, el - 1, io_idx_si)
            elem.err_liveveg = _getr_off(this, :ir_liveveg_err_el, el - 1, io_idx_si)
            ifb.iflux_litter = _getr_off(this, :ir_litter_intflux_el, el - 1, io_idx_si)
            elem.err_litter = _getr_off(this, :ir_litter_err_el, el - 1, io_idx_si)
        end
    end

    site.spread = _getr(this, :ir_spread_si, io_idx_si)

    io_idx_si_sc = io_idx_base
    for i_scls in 1:nsc
        site.demotion_rate[i_scls] = _getr(this, :ir_demorate_sisc, io_idx_si_sc)
        site.promotion_rate[i_scls] = _getr(this, :ir_promrate_sisc, io_idx_si_sc)
        io_idx_si_sc += 1
    end

    site.term_crownarea_canopy = _getr(this, :ir_termcarea_cano_si, io_idx_si)
    site.term_crownarea_ustory = _getr(this, :ir_termcarea_usto_si, io_idx_si)
    site.ema_npp = _getr(this, :ir_emanpp_si, io_idx_si)

    if hlm_use_tree_damage[] == itrue
        io_idx_si_cdpf = io_idx_base
        io_idx_si_cdsc = io_idx_base
        for i_scls in 1:nsc
            for i_cdam in 1:nlevdamage[]
                for i_pft in 1:npft
                    site.imort_rate_damage[i_cdam, i_scls, i_pft] = _getr(this, :ir_imortrate_sicdpf, io_idx_si_cdpf)
                    site.term_nindivs_canopy_damage[i_cdam, i_scls, i_pft] = _getr(this, :ir_termnindiv_cano_sicdpf, io_idx_si_cdpf)
                    site.term_nindivs_ustory_damage[i_cdam, i_scls, i_pft] = _getr(this, :ir_termnindiv_usto_sicdpf, io_idx_si_cdpf)
                    site.imort_cflux_damage[i_cdam, i_scls] = _getr(this, :ir_imortcflux_sicdsc, io_idx_si_cdsc)
                    site.term_cflux_canopy_damage[i_cdam, i_scls] = _getr(this, :ir_termcflux_cano_sicdsc, io_idx_si_cdsc)
                    site.term_cflux_ustory_damage[i_cdam, i_scls] = _getr(this, :ir_termcflux_usto_sicdsc, io_idx_si_cdsc)
                    site.fmort_rate_canopy_damage[i_cdam, i_scls, i_pft] = _getr(this, :ir_fmortrate_cano_sicdpf, io_idx_si_cdpf)
                    site.fmort_rate_ustory_damage[i_cdam, i_scls, i_pft] = _getr(this, :ir_fmortrate_usto_sicdpf, io_idx_si_cdpf)
                    site.fmort_cflux_canopy_damage[i_cdam, i_scls] = _getr(this, :ir_fmortcflux_cano_sicdsc, io_idx_si_cdsc)
                    site.fmort_cflux_ustory_damage[i_cdam, i_scls] = _getr(this, :ir_fmortcflux_usto_sicdsc, io_idx_si_cdsc)
                    io_idx_si_cdsc += 1
                    io_idx_si_cdpf += 1
                end
            end
        end
        site.crownarea_canopy_damage = _getr(this, :ir_crownarea_cano_si, io_idx_si)
        site.crownarea_ustory_damage = _getr(this, :ir_crownarea_usto_si, io_idx_si)
    end

    site.demotion_carbonflux = _getr(this, :ir_democflux_si, io_idx_si)
    site.promotion_carbonflux = _getr(this, :ir_promcflux_si, io_idx_si)
    site.imort_crownarea = _getr(this, :ir_imortcarea_si, io_idx_si)
    site.fmort_crownarea_canopy = _getr(this, :ir_fmortcarea_cano_si, io_idx_si)
    site.fmort_crownarea_ustory = _getr(this, :ir_fmortcarea_usto_si, io_idx_si)

    io_idx_si_pfcl = io_idx_base
    for i in 1:nclmax
        for i_pft in 1:npft
            site.rec_l2fr[i_pft, i] = _getr(this, :ir_recl2fr_sipfcl, io_idx_si_pfcl)
            io_idx_si_pfcl += 1
        end
    end
    io_idx_si_wmem = io_idx_base
    for i in 1:numWaterMem
        for i_pft in 1:npft
            site.liqvol_memory[i, i_pft] = _getr(this, :ir_liqvolmem_siwmft, io_idx_si_wmem)
            site.smp_memory[i, i_pft] = _getr(this, :ir_smpmem_siwmft, io_idx_si_wmem)
            io_idx_si_wmem += 1
        end
    end
    io_idx_si_vtmem = io_idx_base
    for i in 1:num_vegtemp_mem
        site.vegtemp_memory[i] = _getr(this, :ir_vegtempmem_sitm, io_idx_si_vtmem)
        io_idx_si_vtmem += 1
    end

    if hlm_use_planthydro[] == itrue && site.si_hydr !== nothing
        sh = site.si_hydr
        sh.h2oveg_recruit = _getr(this, :ir_hydro_recruit_si, io_idx_si)
        sh.h2oveg_dead = _getr(this, :ir_hydro_dead_si, io_idx_si)
        sh.h2oveg_growturn_err = _getr(this, :ir_hydro_growturn_err_si, io_idx_si)
        sh.h2oveg_hydro_err = _getr(this, :ir_hydro_hydro_err_si, io_idx_si)
        io_idx_si_lyr_shell = io_idx_base
        for i in 1:sh.nlevrhiz
            for k in 1:nshell
                sh.h2osoi_liqvol_shell[i, k] = _getr(this, :ir_hydro_liqvol_shell_si, io_idx_si_lyr_shell)
                io_idx_si_lyr_shell += 1
            end
        end
    end
    return nothing
end

# =====================================================================================
# R4: patch litter pack & unpack (per element, multi-dimensional nested loops).
# =====================================================================================

"""
    _pack_patch_litter!(this, site, cpatch, io_idx_co_1st)

Pack the patch litter pools (`cpatch.litter[el]`) for all elements, mirroring the
Fortran `el=0..num_elements-1` block (here `el=1..num_elements`, handle offset
`el-1`). Each element resets all per-element cursors to the patch base.
"""
function _pack_patch_litter!(this::fates_restart_interface_type, site, cpatch,
                             io_idx_co_1st::Integer)
    hlm_use_sp[] == ifalse || return nothing
    npft = numpft[]
    nlevsoil = site.nlevsoil
    for el in 1:num_elements[]
        io_idx_pa_pft  = io_idx_co_1st
        io_idx_pa_cwd  = io_idx_co_1st
        io_idx_pa_cwsl = io_idx_co_1st
        io_idx_pa_dcsl = io_idx_co_1st
        io_idx_pa_dc   = io_idx_co_1st
        litt = cpatch.litter[el]
        off = el - 1
        for i in 1:npft
            _setr_off!(this, :ir_seed_litt, off, io_idx_pa_pft, litt.seed[i])
            _setr_off!(this, :ir_seedgerm_litt, off, io_idx_pa_pft, litt.seed_germ[i])
            _setr_off!(this, :ir_seed_decay_litt, off, io_idx_pa_pft, litt.seed_decay[i])
            _setr_off!(this, :ir_seedgerm_decay_litt, off, io_idx_pa_pft, litt.seed_germ_decay[i])
            io_idx_pa_pft += 1
        end
        for i in 1:ndcmpy
            _setr_off!(this, :ir_leaf_litt, off, io_idx_pa_dc, litt.leaf_fines[i])
            _setr_off!(this, :ir_lfines_frag_litt, off, io_idx_pa_dc, litt.leaf_fines_frag[i])
            io_idx_pa_dc += 1
            for ilyr in 1:nlevsoil
                _setr_off!(this, :ir_fnrt_litt, off, io_idx_pa_dcsl, litt.root_fines[i, ilyr])
                _setr_off!(this, :ir_rfines_frag_litt, off, io_idx_pa_dcsl, litt.root_fines_frag[i, ilyr])
                io_idx_pa_dcsl += 1
            end
        end
        for i in 1:ncwd
            _setr_off!(this, :ir_agcwd_litt, off, io_idx_pa_cwd, litt.ag_cwd[i])
            _setr_off!(this, :ir_agcwd_frag_litt, off, io_idx_pa_cwd, litt.ag_cwd_frag[i])
            io_idx_pa_cwd += 1
            for ilyr in 1:nlevsoil
                _setr_off!(this, :ir_bgcwd_litt, off, io_idx_pa_cwsl, litt.bg_cwd[i, ilyr])
                _setr_off!(this, :ir_bgcwd_frag_litt, off, io_idx_pa_cwsl, litt.bg_cwd_frag[i, ilyr])
                io_idx_pa_cwsl += 1
            end
        end
    end
    return nothing
end

"""
    _unpack_patch_litter!(this, site, cpatch, io_idx_co_1st)

Read back the patch litter pools (symmetric to [`_pack_patch_litter!`](@ref)).
"""
function _unpack_patch_litter!(this::fates_restart_interface_type, site, cpatch,
                               io_idx_co_1st::Integer)
    hlm_use_sp[] == ifalse || return nothing
    npft = numpft[]
    nlevsoil = site.nlevsoil
    for el in 1:num_elements[]
        io_idx_pa_pft  = io_idx_co_1st
        io_idx_pa_cwd  = io_idx_co_1st
        io_idx_pa_cwsl = io_idx_co_1st
        io_idx_pa_dcsl = io_idx_co_1st
        io_idx_pa_dc   = io_idx_co_1st
        litt = cpatch.litter[el]
        off = el - 1
        for i in 1:npft
            litt.seed[i] = _getr_off(this, :ir_seed_litt, off, io_idx_pa_pft)
            litt.seed_germ[i] = _getr_off(this, :ir_seedgerm_litt, off, io_idx_pa_pft)
            litt.seed_decay[i] = _getr_off(this, :ir_seed_decay_litt, off, io_idx_pa_pft)
            litt.seed_germ_decay[i] = _getr_off(this, :ir_seedgerm_decay_litt, off, io_idx_pa_pft)
            io_idx_pa_pft += 1
        end
        for i in 1:ndcmpy
            litt.leaf_fines[i] = _getr_off(this, :ir_leaf_litt, off, io_idx_pa_dc)
            litt.leaf_fines_frag[i] = _getr_off(this, :ir_lfines_frag_litt, off, io_idx_pa_dc)
            io_idx_pa_dc += 1
            for ilyr in 1:nlevsoil
                litt.root_fines[i, ilyr] = _getr_off(this, :ir_fnrt_litt, off, io_idx_pa_dcsl)
                litt.root_fines_frag[i, ilyr] = _getr_off(this, :ir_rfines_frag_litt, off, io_idx_pa_dcsl)
                io_idx_pa_dcsl += 1
            end
        end
        for i in 1:ncwd
            litt.ag_cwd[i] = _getr_off(this, :ir_agcwd_litt, off, io_idx_pa_cwd)
            litt.ag_cwd_frag[i] = _getr_off(this, :ir_agcwd_frag_litt, off, io_idx_pa_cwd)
            io_idx_pa_cwd += 1
            for ilyr in 1:nlevsoil
                litt.bg_cwd[i, ilyr] = _getr_off(this, :ir_bgcwd_litt, off, io_idx_pa_cwsl)
                litt.bg_cwd_frag[i, ilyr] = _getr_off(this, :ir_bgcwd_frag_litt, off, io_idx_pa_cwsl)
                io_idx_pa_cwsl += 1
            end
        end
    end
    return nothing
end

# =====================================================================================
# set_restart_vectors! — PACK the linked-list demographic state into flat vectors.
# =====================================================================================

"""
    set_restart_vectors!(this, nc, nsites, sites)

Walk each site's oldest->youngest patch list, and within each patch the
shortest->taller cohort list, packing the DEMOGRAPHIC skeleton and load-bearing
fields into the registered restart vectors. Writes the per-site
`fates_PatchesPerSite` and per-patch `fates_CohortsPerPatch` COUNTS that drive
the unpack. Mirrors the Fortran `set_restart_vectors`.

Packs the demographic skeleton + load-bearing fields, plus the B18-followup
fill-groups: R1 (cohort diagnostic scalars + mortality fluxes), R2 (patch
fuel/scorch + radiation/albedo + site hydraulics scalars), R4 (patch litter
blocks), R5 (site demographic/mortality/damage/flux-diag arrays), R6 (patch
running means + site disturbance rates), R7 (cohort PRT pools) and R8 (cohort
hydraulics, gated on `hlm_use_planthydro`). The restart fill is complete.
"""
function set_restart_vectors!(this::fates_restart_interface_type, nc::Integer,
                              nsites::Integer, sites::AbstractVector{<:ed_site_type})
    flush_rvars!(this, nc)

    maxperpatch = fates_maxElementsPerPatch[]

    for s in 1:nsites
        site = sites[s]
        io_idx_si     = this.restart_map[nc].site_index[s]
        io_idx_co_1st = this.restart_map[nc].cohort1_index[s]

        # --- site-level scalars (after-patch in Fortran, order-independent here) ---
        _seti!(this, :ir_cd_status_si, io_idx_si, site.cstatus)
        _seti!(this, :ir_nchill_days_si, io_idx_si, site.nchilldays)
        _seti!(this, :ir_ncold_days_si, io_idx_si, site.ncolddays)
        _seti!(this, :ir_cleafondate_si, io_idx_si, site.cleafondate)
        _seti!(this, :ir_cleafoffdate_si, io_idx_si, site.cleafoffdate)
        _seti!(this, :ir_cndaysleafon_si, io_idx_si, site.cndaysleafon)
        _seti!(this, :ir_cndaysleafoff_si, io_idx_si, site.cndaysleafoff)
        _seti!(this, :ir_phenmodeldate_si, io_idx_si, site.phen_model_date)
        _setr!(this, :ir_gdd_si, io_idx_si, site.grow_deg_days)
        _setr!(this, :ir_snow_depth_si, io_idx_si, site.snow_depth)
        _setr!(this, :ir_trunk_product_si, io_idx_si, site.resources_management.trunk_product_site)
        _seti!(this, :ir_landuse_config_si, io_idx_si, hlm_use_potentialveg[])

        # capture the site's cohort base BEFORE the patch loop advances it (the
        # R5 site-level arrays all index from this initial base).
        io_idx_base = io_idx_co_1st

        # --- patch loop (oldest -> youngest) ---
        patchespersite = 0
        cpatch = site.oldest_patch
        while cpatch !== nothing
            patchespersite += 1
            io_idx_co = io_idx_co_1st

            # --- cohort loop (shortest -> taller) ---
            cohortsperpatch = 0
            ccohort = cpatch.shortest
            while ccohort !== nothing
                cohortsperpatch += 1

                _seti!(this, :ir_canopy_layer_co, io_idx_co, ccohort.canopy_layer)
                _setr!(this, :ir_canopy_layer_yesterday_co, io_idx_co, ccohort.canopy_layer_yesterday)
                _seti!(this, :ir_crowndamage_co, io_idx_co, ccohort.crowndamage)
                _setr!(this, :ir_canopy_trim_co, io_idx_co, ccohort.canopy_trim)
                _setr!(this, :ir_l2fr_co, io_idx_co, ccohort.l2fr)
                _setr!(this, :ir_seed_prod_co, io_idx_co, ccohort.seed_prod)
                _seti!(this, :ir_size_class_lasttimestep_co, io_idx_co, ccohort.size_class_lasttimestep)
                _setr!(this, :ir_dbh_co, io_idx_co, ccohort.dbh)
                _setr!(this, :ir_coage_co, io_idx_co, ccohort.coage)
                _setr!(this, :ir_height_co, io_idx_co, ccohort.height)
                _setr!(this, :ir_g_sb_laweight_co, io_idx_co, ccohort.g_sb_laweight)
                _setr!(this, :ir_nplant_co, io_idx_co, ccohort.n)
                _seti!(this, :ir_pft_co, io_idx_co, ccohort.pft)
                _seti!(this, :ir_status_co, io_idx_co, ccohort.status_coh)
                _setr!(this, :ir_efleaf_co, io_idx_co, ccohort.efleaf_coh)
                _setr!(this, :ir_effnrt_co, io_idx_co, ccohort.effnrt_coh)
                _setr!(this, :ir_efstem_co, io_idx_co, ccohort.efstem_coh)
                _seti!(this, :ir_isnew_co, io_idx_co, ccohort.isnew ? new_cohort : old_cohort)
                _setr!(this, :ir_c_area_co, io_idx_co, ccohort.c_area)
                _setr!(this, :ir_treelai_co, io_idx_co, ccohort.treelai)
                _setr!(this, :ir_treesai_co, io_idx_co, ccohort.treesai)

                # --- R1: cohort diagnostic scalars + mortality fluxes ---
                _setr!(this, :ir_gpp_acc_co, io_idx_co, ccohort.gpp_acc)
                _setr!(this, :ir_npp_acc_co, io_idx_co, ccohort.npp_acc)
                _setr!(this, :ir_resp_acc_co, io_idx_co, ccohort.resp_acc)
                _setr!(this, :ir_gpp_acc_hold_co, io_idx_co, ccohort.gpp_acc_hold)
                _setr!(this, :ir_npp_acc_hold_co, io_idx_co, ccohort.npp_acc_hold)
                _setr!(this, :ir_resp_acc_hold_co, io_idx_co, ccohort.resp_acc_hold)
                _setr!(this, :ir_resp_excess_co, io_idx_co, ccohort.resp_excess)
                _setr!(this, :ir_bmort_co, io_idx_co, ccohort.bmort)
                _setr!(this, :ir_hmort_co, io_idx_co, ccohort.hmort)
                _setr!(this, :ir_cmort_co, io_idx_co, ccohort.cmort)
                _setr!(this, :ir_smort_co, io_idx_co, ccohort.smort)
                _setr!(this, :ir_asmort_co, io_idx_co, ccohort.asmort)
                _setr!(this, :ir_dgmort_co, io_idx_co, ccohort.dgmort)
                _setr!(this, :ir_frmort_co, io_idx_co, ccohort.frmort)
                _setr!(this, :ir_lmort_direct_co, io_idx_co, ccohort.lmort_direct)
                _setr!(this, :ir_lmort_collateral_co, io_idx_co, ccohort.lmort_collateral)
                _setr!(this, :ir_lmort_infra_co, io_idx_co, ccohort.lmort_infra)
                _setr!(this, :ir_ddbhdt_co, io_idx_co, ccohort.ddbhdt)
                _setr!(this, :ir_resp_tstep_co, io_idx_co, ccohort.resp_tstep)

                # --- R7: PRT pool block (instantaneous state + 3 fluxes per pos) ---
                _pack_cohort_prt!(this, ccohort, io_idx_co)

                # --- R8: cohort hydraulics (gated; th_*/errh2o) ---
                if hlm_use_planthydro[] == itrue && ccohort.co_hydr !== nothing
                    ch = ccohort.co_hydr
                    SetCohortRealVector(this, ch.th_ag, n_hypool_ag,
                                        this.ir[:ir_hydro_th_ag_covec], io_idx_co)
                    SetCohortRealVector(this, ch.th_aroot, site.si_hydr.nlevrhiz,
                                        this.ir[:ir_hydro_th_aroot_covec], io_idx_co)
                    _setr!(this, :ir_hydro_th_troot, io_idx_co, ch.th_troot)
                    _setr!(this, :ir_hydro_errh2o, io_idx_co, ch.errh2o)
                end

                io_idx_co += 1
                ccohort = ccohort.taller
            end

            # --- patch scalars (at the patch base address) ---
            _setr!(this, :ir_livegrass_pa, io_idx_co_1st, cpatch.livegrass)
            _setr!(this, :ir_age_pa, io_idx_co_1st, cpatch.age)
            _seti!(this, :ir_patchdistturbcat_pa, io_idx_co_1st, cpatch.land_use_label)
            _setr!(this, :ir_agesinceanthrodist_pa, io_idx_co_1st, cpatch.age_since_anthro_disturbance)
            _seti!(this, :ir_nocomp_pft_label_pa, io_idx_co_1st, cpatch.nocomp_pft_label)
            _setr!(this, :ir_area_pa, io_idx_co_1st, cpatch.area)
            _seti!(this, :ir_ncohort_pa, io_idx_co_1st, cohortsperpatch)   # KEY COUNT
            _setr!(this, :ir_fcansno_pa, io_idx_co_1st, cpatch.fcansno)
            _seti!(this, :ir_solar_zenith_flag_pa, io_idx_co_1st, cpatch.solar_zenith_flag ? 1 : 0)
            _setr!(this, :ir_solar_zenith_angle_pa, io_idx_co_1st, cpatch.solar_zenith_angle)
            _seti!(this, :ir_nclp_pa, io_idx_co_1st, cpatch.ncl_p)
            _setr!(this, :ir_zstar_pa, io_idx_co_1st, cpatch.zstar)

            # --- R2: patch fuel/scorch + radiation/albedo (mixed-dim, at base) ---
            if hlm_use_sp[] == ifalse
                io_idx_pa_pft = io_idx_co_1st
                for i in 1:numpft[]
                    _setr!(this, :ir_scorch_ht_pa_pft, io_idx_pa_pft, cpatch.scorch_ht[i])
                    io_idx_pa_pft += 1
                end
                io_idx_pa_fc = io_idx_co_1st
                for i in 1:num_fuel_classes
                    moist = cpatch.fuel === nothing ? 0.0 : cpatch.fuel.effective_moisture[i]
                    _setr!(this, :ir_litter_moisture_pa_nfsc, io_idx_pa_fc, moist)
                    io_idx_pa_fc += 1
                end
            end
            io_idx_pa_ib = io_idx_co_1st
            for i in 1:num_swb
                _setr!(this, :ir_gnd_alb_dif_pasb, io_idx_pa_ib, cpatch.gnd_alb_dif[i])
                _setr!(this, :ir_gnd_alb_dir_pasb, io_idx_pa_ib, cpatch.gnd_alb_dir[i])
                io_idx_pa_ib += 1
            end

            # --- R4: patch litter blocks (per element nested loops) ---
            _pack_patch_litter!(this, site, cpatch, io_idx_co_1st)

            # --- R6: patch running means (24-hr veg temp + photo-acclim EMAs) ---
            cpatch.tveg24        === nothing || SetRMeanRestartVar(this, cpatch.tveg24,        this.ir[:ir_tveg24_pa],        io_idx_co_1st)
            cpatch.tveg_lpa      === nothing || SetRMeanRestartVar(this, cpatch.tveg_lpa,      this.ir[:ir_tveglpa_pa],       io_idx_co_1st)
            cpatch.tveg_longterm === nothing || SetRMeanRestartVar(this, cpatch.tveg_longterm, this.ir[:ir_tveglongterm_pa],  io_idx_co_1st)

            io_idx_co_1st += maxperpatch
            cpatch = cpatch.younger
        end

        # --- R5: site-level demographic / mortality / damage / flux arrays ---
        _pack_site_diagnostics!(this, site, io_idx_si, io_idx_base)

        _seti!(this, :ir_npatch_si, io_idx_si, patchespersite)   # KEY COUNT
    end
    return nothing
end

# =====================================================================================
# create_patchcohort_structure! — allocate the empty site -> patch -> cohort
# skeleton from the per-site / per-patch COUNTS.
# =====================================================================================

"""
    create_patchcohort_structure!(this, nc, nsites, sites; num_swb, current_tod)

Read only the `fates_PatchesPerSite` (ir_npatch_si) and `fates_CohortsPerPatch`
(ir_ncohort_pa) counts and allocate that many patches per site / cohorts per
patch, building the age-ordered patch list (oldest = idx 1) and the
height-ordered cohort list (first created = tallest, last = shortest). Mirrors
the Fortran `create_patchcohort_structure` skeleton-build. Each freshly-created
cohort gets its PRT (PARTEH) object re-created via `InitPRTObject!`, and — when
`hlm_use_planthydro` is on — its hydraulics object via `InitHydrCohort`, so the
R7/R8 pools have somewhere to land on unpack. (The HLM `init_site_vars` /
`zero_site` boundary-condition setup has no standalone equivalent yet — `# TODO:`
below.)
"""
function create_patchcohort_structure!(this::fates_restart_interface_type, nc::Integer,
                                       nsites::Integer, sites::AbstractVector{<:ed_site_type};
                                       n_swb::Integer = num_swb, current_tod::Integer = 0)
    maxperpatch = fates_maxElementsPerPatch[]
    regen = ed_params().regeneration_model

    for s in 1:nsites
        site = sites[s]
        io_idx_si     = this.restart_map[nc].site_index[s]
        io_idx_co_1st = this.restart_map[nc].cohort1_index[s]

        # TODO: init_site_vars(site, bc_in, bc_out); zero_site(site) — HLM BC setup.

        npatches = _geti(this, :ir_npatch_si, io_idx_si)
        if npatches < 0 || npatches > 10000
            fates_endrun("create_patchcohort_structure!: nonsensical npatches=$npatches")
        end

        site.youngest_patch = nothing
        site.oldest_patch   = nothing

        nlevsoil = site.nlevsoil > 0 ? site.nlevsoil : max(length(site.dz_soil), 1)

        for idx_pa in 1:npatches
            newp = fates_patch_type()
            Create!(newp, fates_unset_r8, fates_unset_r8, primaryland, fates_unset_int,
                    n_swb, numpft[], nlevsoil, current_tod, regen)
            newp.patchno = idx_pa

            # --- cohort list build: first created = tallest, last = shortest ---
            newp.tallest  = nothing
            newp.shortest = nothing
            prev_cohort   = nothing

            ncohorts = _geti(this, :ir_ncohort_pa, io_idx_co_1st)
            for fto in 1:ncohorts
                new_cohort_node = fates_cohort_type()
                # NanValues / ZeroValues equivalent: a fresh cohort defaults NaN/unset.
                if newp.tallest === nothing
                    newp.tallest = new_cohort_node     # 1st created = tallest
                end
                if prev_cohort !== nothing
                    new_cohort_node.taller = prev_cohort
                    prev_cohort.shorter    = new_cohort_node
                end
                newp.shortest = new_cohort_node        # every new cohort becomes shortest

                # Initialize the PARTEH (PRT) object so the unpack can fill its
                # pools (restart cohorts start with `prt === nothing`). Then, if
                # plant hydraulics is on, allocate the per-cohort hydro object.
                new_cohort_node.prt = InitPRTObject!()
                if hlm_use_planthydro[] == itrue
                    InitHydrCohort(site, new_cohort_node)
                end

                prev_cohort = new_cohort_node
            end

            # --- patch insertion: always appended at the youngest end (idx order) ---
            if idx_pa == 1
                site.youngest_patch = newp
                site.oldest_patch   = newp
                newp.younger = nothing
                newp.older   = nothing
            else
                newp.older   = site.youngest_patch
                newp.younger = nothing
                site.youngest_patch.younger = newp
                site.youngest_patch = newp
            end

            io_idx_co_1st += maxperpatch
        end
    end
    return nothing
end

# =====================================================================================
# get_restart_vectors! — UNPACK the flat vectors into the rebuilt skeleton.
# =====================================================================================

"""
    get_restart_vectors!(this, nc, nsites, sites)

Walk the (already-allocated) skeleton and fill patch + cohort fields from the
restart vectors, recomputing derived indices (age_class) and asserting the
round-trip counts (`cohortsperpatch == fates_CohortsPerPatch` per patch,
`patchespersite == fates_PatchesPerSite` per site). Mirrors the Fortran
`get_restart_vectors` for the demographic skeleton + load-bearing fields.

Reads back the demographic skeleton + the B18-followup R1/R2/R4/R5/R6/R7/R8
fill-groups symmetric to the pack side, including the PRT (R7) finalizers
`InitPRTBoundaryConditions` / `UpdateCohortBioPhysRates` and the hydraulics (R8)
finalizer `UpdatePlantPsiFTCFromTheta!`.
"""
function get_restart_vectors!(this::fates_restart_interface_type, nc::Integer,
                              nsites::Integer, sites::AbstractVector{<:ed_site_type})
    maxperpatch = fates_maxElementsPerPatch[]

    for s in 1:nsites
        site = sites[s]
        io_idx_si     = this.restart_map[nc].site_index[s]
        io_idx_co_1st = this.restart_map[nc].cohort1_index[s]

        # --- site-level scalars ---
        site.cstatus         = _geti(this, :ir_cd_status_si, io_idx_si)
        site.nchilldays      = _geti(this, :ir_nchill_days_si, io_idx_si)
        site.ncolddays       = _geti(this, :ir_ncold_days_si, io_idx_si)
        site.cleafondate     = _geti(this, :ir_cleafondate_si, io_idx_si)
        site.cleafoffdate    = _geti(this, :ir_cleafoffdate_si, io_idx_si)
        site.cndaysleafon    = _geti(this, :ir_cndaysleafon_si, io_idx_si)
        site.cndaysleafoff   = _geti(this, :ir_cndaysleafoff_si, io_idx_si)
        site.phen_model_date = _geti(this, :ir_phenmodeldate_si, io_idx_si)
        site.grow_deg_days   = _getr(this, :ir_gdd_si, io_idx_si)
        site.snow_depth      = _getr(this, :ir_snow_depth_si, io_idx_si)
        site.resources_management.trunk_product_site = _getr(this, :ir_trunk_product_si, io_idx_si)

        # capture the site's cohort base BEFORE the patch loop advances it.
        io_idx_base = io_idx_co_1st

        # --- patch loop (oldest -> youngest) ---
        patchespersite = 0
        cpatch = site.oldest_patch
        while cpatch !== nothing
            patchespersite += 1
            io_idx_co = io_idx_co_1st

            # --- cohort loop (shortest -> taller) ---
            cohortsperpatch = 0
            ccohort = cpatch.shortest
            while ccohort !== nothing
                cohortsperpatch += 1

                ccohort.canopy_layer           = _geti(this, :ir_canopy_layer_co, io_idx_co)
                ccohort.canopy_layer_yesterday = _getr(this, :ir_canopy_layer_yesterday_co, io_idx_co)
                ccohort.crowndamage            = _geti(this, :ir_crowndamage_co, io_idx_co)
                ccohort.canopy_trim            = _getr(this, :ir_canopy_trim_co, io_idx_co)
                ccohort.l2fr                   = _getr(this, :ir_l2fr_co, io_idx_co)
                ccohort.seed_prod              = _getr(this, :ir_seed_prod_co, io_idx_co)
                ccohort.size_class_lasttimestep = _geti(this, :ir_size_class_lasttimestep_co, io_idx_co)
                ccohort.dbh                    = _getr(this, :ir_dbh_co, io_idx_co)
                ccohort.coage                  = _getr(this, :ir_coage_co, io_idx_co)
                ccohort.g_sb_laweight          = _getr(this, :ir_g_sb_laweight_co, io_idx_co)
                ccohort.height                 = _getr(this, :ir_height_co, io_idx_co)
                ccohort.n                      = _getr(this, :ir_nplant_co, io_idx_co)
                ccohort.pft                    = _geti(this, :ir_pft_co, io_idx_co)
                ccohort.status_coh             = _geti(this, :ir_status_co, io_idx_co)
                ccohort.efleaf_coh             = _getr(this, :ir_efleaf_co, io_idx_co)
                ccohort.effnrt_coh             = _getr(this, :ir_effnrt_co, io_idx_co)
                ccohort.efstem_coh             = _getr(this, :ir_efstem_co, io_idx_co)
                ccohort.isnew                  = (_geti(this, :ir_isnew_co, io_idx_co) == new_cohort)
                ccohort.c_area                 = _getr(this, :ir_c_area_co, io_idx_co)
                ccohort.treelai                = _getr(this, :ir_treelai_co, io_idx_co)
                ccohort.treesai                = _getr(this, :ir_treesai_co, io_idx_co)

                # --- R1: cohort diagnostic scalars + mortality fluxes ---
                ccohort.gpp_acc          = _getr(this, :ir_gpp_acc_co, io_idx_co)
                ccohort.npp_acc          = _getr(this, :ir_npp_acc_co, io_idx_co)
                ccohort.resp_acc         = _getr(this, :ir_resp_acc_co, io_idx_co)
                ccohort.gpp_acc_hold     = _getr(this, :ir_gpp_acc_hold_co, io_idx_co)
                ccohort.npp_acc_hold     = _getr(this, :ir_npp_acc_hold_co, io_idx_co)
                ccohort.resp_acc_hold    = _getr(this, :ir_resp_acc_hold_co, io_idx_co)
                ccohort.resp_excess      = _getr(this, :ir_resp_excess_co, io_idx_co)
                ccohort.bmort            = _getr(this, :ir_bmort_co, io_idx_co)
                ccohort.hmort            = _getr(this, :ir_hmort_co, io_idx_co)
                ccohort.cmort            = _getr(this, :ir_cmort_co, io_idx_co)
                ccohort.smort            = _getr(this, :ir_smort_co, io_idx_co)
                ccohort.asmort           = _getr(this, :ir_asmort_co, io_idx_co)
                ccohort.dgmort           = _getr(this, :ir_dgmort_co, io_idx_co)
                ccohort.frmort           = _getr(this, :ir_frmort_co, io_idx_co)
                ccohort.lmort_direct     = _getr(this, :ir_lmort_direct_co, io_idx_co)
                ccohort.lmort_collateral = _getr(this, :ir_lmort_collateral_co, io_idx_co)
                ccohort.lmort_infra      = _getr(this, :ir_lmort_infra_co, io_idx_co)
                ccohort.ddbhdt           = _getr(this, :ir_ddbhdt_co, io_idx_co)
                ccohort.resp_tstep       = _getr(this, :ir_resp_tstep_co, io_idx_co)

                # --- R7: PRT pool block (the cohort.prt object was re-created by
                #     create_patchcohort_structure!'s InitPRTObject! call) ---
                _unpack_cohort_prt!(this, ccohort, io_idx_co)

                # PRT finalizers: re-point boundary conditions and refresh the
                # cohort biophysical rates now that pft + pools are filled.
                InitPRTBoundaryConditions(ccohort)
                UpdateCohortBioPhysRates(ccohort)

                # --- R8: cohort hydraulics (gated; th_*/errh2o + psi/ftc refresh) ---
                if hlm_use_planthydro[] == itrue && ccohort.co_hydr !== nothing
                    ch = ccohort.co_hydr
                    GetCohortRealVector(this, ch.th_ag, n_hypool_ag,
                                        this.ir[:ir_hydro_th_ag_covec], io_idx_co)
                    GetCohortRealVector(this, ch.th_aroot, site.si_hydr.nlevrhiz,
                                        this.ir[:ir_hydro_th_aroot_covec], io_idx_co)
                    ch.th_troot = _getr(this, :ir_hydro_th_troot, io_idx_co)
                    ch.errh2o   = _getr(this, :ir_hydro_errh2o, io_idx_co)
                    UpdatePlantPsiFTCFromTheta!(ccohort, site.si_hydr)
                end

                io_idx_co += 1
                ccohort = ccohort.taller
            end

            # --- consistency check (round-trip guarantee) ---
            ncohort_stored = _geti(this, :ir_ncohort_pa, io_idx_co_1st)
            if cohortsperpatch != ncohort_stored
                fates_endrun("get_restart_vectors!: cohortsperpatch ($cohortsperpatch) != stored ($ncohort_stored)")
            end

            # --- patch-level fields ---
            cpatch.livegrass        = _getr(this, :ir_livegrass_pa, io_idx_co_1st)
            cpatch.age              = _getr(this, :ir_age_pa, io_idx_co_1st)
            cpatch.land_use_label   = _geti(this, :ir_patchdistturbcat_pa, io_idx_co_1st)
            cpatch.age_since_anthro_disturbance = _getr(this, :ir_agesinceanthrodist_pa, io_idx_co_1st)
            cpatch.nocomp_pft_label = _geti(this, :ir_nocomp_pft_label_pa, io_idx_co_1st)
            cpatch.area             = _getr(this, :ir_area_pa, io_idx_co_1st)
            cpatch.age_class        = get_age_class_index(cpatch.age)   # recomputed, not stored
            cpatch.fcansno          = _getr(this, :ir_fcansno_pa, io_idx_co_1st)
            cpatch.solar_zenith_flag  = (_geti(this, :ir_solar_zenith_flag_pa, io_idx_co_1st) != 0)
            cpatch.solar_zenith_angle = _getr(this, :ir_solar_zenith_angle_pa, io_idx_co_1st)
            cpatch.ncl_p            = _geti(this, :ir_nclp_pa, io_idx_co_1st)
            cpatch.zstar            = _getr(this, :ir_zstar_pa, io_idx_co_1st)
            cpatch.countcohorts     = cohortsperpatch

            # --- R2: patch fuel/scorch + radiation/albedo ---
            if hlm_use_sp[] == ifalse
                io_idx_pa_pft = io_idx_co_1st
                for i in 1:numpft[]
                    cpatch.scorch_ht[i] = _getr(this, :ir_scorch_ht_pa_pft, io_idx_pa_pft)
                    io_idx_pa_pft += 1
                end
                if cpatch.fuel !== nothing
                    io_idx_pa_fc = io_idx_co_1st
                    for i in 1:num_fuel_classes
                        cpatch.fuel.effective_moisture[i] = _getr(this, :ir_litter_moisture_pa_nfsc, io_idx_pa_fc)
                        io_idx_pa_fc += 1
                    end
                end
            end
            io_idx_pa_ib = io_idx_co_1st
            for i in 1:num_swb
                cpatch.gnd_alb_dif[i] = _getr(this, :ir_gnd_alb_dif_pasb, io_idx_pa_ib)
                cpatch.gnd_alb_dir[i] = _getr(this, :ir_gnd_alb_dir_pasb, io_idx_pa_ib)
                io_idx_pa_ib += 1
            end

            # --- R4: patch litter blocks ---
            _unpack_patch_litter!(this, site, cpatch, io_idx_co_1st)

            # --- R6: patch running means (the rmean_type instances were created by
            #     Create! during create_patchcohort_structure!) ---
            cpatch.tveg24        === nothing || GetRMeanRestartVar(this, cpatch.tveg24,        this.ir[:ir_tveg24_pa],        io_idx_co_1st)
            cpatch.tveg_lpa      === nothing || GetRMeanRestartVar(this, cpatch.tveg_lpa,      this.ir[:ir_tveglpa_pa],       io_idx_co_1st)
            cpatch.tveg_longterm === nothing || GetRMeanRestartVar(this, cpatch.tveg_longterm, this.ir[:ir_tveglongterm_pa],  io_idx_co_1st)

            io_idx_co_1st += maxperpatch
            cpatch = cpatch.younger
        end

        # --- R5: site-level demographic / mortality / damage / flux arrays ---
        _unpack_site_diagnostics!(this, site, io_idx_si, io_idx_base)

        # --- consistency check (round-trip guarantee) ---
        npatch_stored = _geti(this, :ir_npatch_si, io_idx_si)
        if patchespersite != npatch_stored
            fates_endrun("get_restart_vectors!: patchespersite ($patchespersite) != stored ($npatch_stored)")
        end
    end
    return nothing
end

# =====================================================================================
# update_3dpatch_radiation — re-derives 3D albedo from packed restart radiation.
# =====================================================================================

"""
    update_3dpatch_radiation!(this, nsites, sites, bc_out)

R3: re-derive the per-patch 3D radiation / albedo boundary outputs after a
restart read, from the restart-stored ground albedos (`gnd_alb_dir/dif`) and
canopy structure. For each patch (ifp incremented over EVERY patch in the
oldest->youngest list, including bare-ground — mirroring the Fortran which does
NOT special-case bare ground for ifp), zero the patch rad-profile fields, and if
`solar_zenith_flag` is set: bare-ground patches (`maximum(nrad[1,:])==0`) copy
the ground albedos into `bc_out.albd/albi_parb`; vegetated patches run the
ported `PatchNormanRadiation` (norman) or `Solve!` (two-stream) solver. Mirrors
the Fortran `update_3dpatch_radiation`.
"""
function update_3dpatch_radiation!(this::fates_restart_interface_type, nsites::Integer,
                                   sites::AbstractVector{<:ed_site_type}, bc_out)
    radiation_model = ed_params().radiation_model

    for s in 1:nsites
        ifp = 0
        currentPatch = sites[s].oldest_patch
        while currentPatch !== nothing
            ifp += 1

            # --- zero the patch radiation-profile fields ---
            isempty(currentPatch.f_sun)      || fill!(currentPatch.f_sun, 0.0)
            isempty(currentPatch.fabd_sun_z) || fill!(currentPatch.fabd_sun_z, 0.0)
            isempty(currentPatch.fabd_sha_z) || fill!(currentPatch.fabd_sha_z, 0.0)
            isempty(currentPatch.fabi_sun_z) || fill!(currentPatch.fabi_sun_z, 0.0)
            isempty(currentPatch.fabi_sha_z) || fill!(currentPatch.fabi_sha_z, 0.0)
            isempty(currentPatch.fabd)       || fill!(currentPatch.fabd, 0.0)
            isempty(currentPatch.fabi)       || fill!(currentPatch.fabi, 0.0)
            isempty(currentPatch.nrmlzd_parprof_pft_dir_z) || fill!(currentPatch.nrmlzd_parprof_pft_dir_z, 0.0)
            isempty(currentPatch.nrmlzd_parprof_pft_dif_z) || fill!(currentPatch.nrmlzd_parprof_pft_dif_z, 0.0)
            isempty(currentPatch.rad_error)  || fill!(currentPatch.rad_error, 0.0)

            if currentPatch.solar_zenith_flag

                for ib in 1:num_swb
                    bc_out[s].albd_parb[ifp, ib] = 0.0
                    bc_out[s].albi_parb[ifp, ib] = 0.0
                    bc_out[s].fabi_parb[ifp, ib] = 0.0
                    bc_out[s].fabd_parb[ifp, ib] = 0.0
                    bc_out[s].ftdd_parb[ifp, ib] = 1.0
                    bc_out[s].ftid_parb[ifp, ib] = 1.0
                    bc_out[s].ftii_parb[ifp, ib] = 1.0
                end

                if maximum(@view currentPatch.nrad[1, :]) == 0
                    # No leaf layers — effectively bare ground.
                    for ib in 1:num_swb
                        bc_out[s].fabd_parb[ifp, ib] = 0.0
                        bc_out[s].fabi_parb[ifp, ib] = 0.0
                        bc_out[s].albd_parb[ifp, ib] = currentPatch.gnd_alb_dir[ib]
                        bc_out[s].albi_parb[ifp, ib] = currentPatch.gnd_alb_dif[ib]
                        bc_out[s].ftdd_parb[ifp, ib] = 1.0
                        bc_out[s].ftid_parb[ifp, ib] = 1.0
                        bc_out[s].ftii_parb[ifp, ib] = 1.0
                    end
                else
                    if radiation_model == norman_solver
                        albd = zeros(num_swb); albi = zeros(num_swb)
                        fabd = zeros(num_swb); fabi = zeros(num_swb)
                        ftdd = zeros(num_swb); ftid = zeros(num_swb)
                        ftii = zeros(num_swb)
                        PatchNormanRadiation(currentPatch, albd, albi, fabd, fabi, ftdd, ftid, ftii)
                        for ib in 1:num_swb
                            bc_out[s].albd_parb[ifp, ib] = albd[ib]
                            bc_out[s].albi_parb[ifp, ib] = albi[ib]
                            bc_out[s].fabd_parb[ifp, ib] = fabd[ib]
                            bc_out[s].fabi_parb[ifp, ib] = fabi[ib]
                            bc_out[s].ftdd_parb[ifp, ib] = ftdd[ib]
                            bc_out[s].ftid_parb[ifp, ib] = ftid[ib]
                            bc_out[s].ftii_parb[ifp, ib] = ftii[ib]
                        end
                    elseif radiation_model == twostr_solver
                        twostr = currentPatch.twostr
                        CanopyPrep!(twostr, currentPatch.fcansno)
                        ZenithPrep!(twostr, currentPatch.solar_zenith_angle)
                        for ib in 1:num_swb
                            twostr.band[ib].albedo_grnd_diff = currentPatch.gnd_alb_dif[ib]
                            twostr.band[ib].albedo_grnd_beam = currentPatch.gnd_alb_dir[ib]
                            (albedo_beam, albedo_diff, consv_err,
                             frac_abs_can_beam, frac_abs_can_diff,
                             frac_beam_grnd_beam, frac_diff_grnd_beam,
                             frac_diff_grnd_diff) =
                                Solve!(twostr, ib, normalized_upper_boundary, 1.0, 1.0,
                                       sites[s].taulambda_2str, sites[s].omega_2str, sites[s].ipiv_2str)
                            bc_out[s].albd_parb[ifp, ib] = albedo_beam
                            bc_out[s].albi_parb[ifp, ib] = albedo_diff
                            currentPatch.rad_error[ib]   = consv_err
                            bc_out[s].fabd_parb[ifp, ib] = frac_abs_can_beam
                            bc_out[s].fabi_parb[ifp, ib] = frac_abs_can_diff
                            bc_out[s].ftdd_parb[ifp, ib] = frac_beam_grnd_beam
                            bc_out[s].ftid_parb[ifp, ib] = frac_diff_grnd_beam
                            bc_out[s].ftii_parb[ifp, ib] = frac_diff_grnd_diff
                        end
                    end
                end
            end

            currentPatch = currentPatch.younger
        end
    end
    return nothing
end
