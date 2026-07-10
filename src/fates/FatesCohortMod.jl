# FatesCohortMod.jl
# Julia port of FATES src/fates/biogeochem/FatesCohortMod.F90
#
# The COHORT type — the fundamental FATES demographic unit. A cohort is a group
# of plants of the same PFT and (approximately) the same size, treated as one
# representative individual scaled by a number density `n`. `fates_cohort_type`
# integrates:
#   * the PRT (Plant Reactive Transport / PARTEH) allocation state object (`prt`)
#   * the per-cohort plant-hydraulics state object (`co_hydr`)
#   * vegetation structure (dbh, height, crown area/damage, PFT, canopy layer)
#   * carbon/nutrient fluxes, respiration components, mortality rates,
#     growth derivatives, fire effects, phenology status
#   * the doubly-linked-list `taller`/`shorter` pointers that order cohorts within
#     a patch by height.
#
# Translation notes (per project conventions):
#   * `fates_r8` -> Float64; cohort integers -> Int.
#   * The Fortran `class(prt_vartypes), pointer :: prt` (polymorphic) ->
#     `Union{AbstractPRTVartypes,Nothing}` (any concrete PRT hypothesis subtype).
#   * `type(ed_cohort_hydr_type), pointer :: co_hydr` ->
#     `Union{ed_cohort_hydr_type,Nothing}`.
#   * The linked-list `taller`/`shorter` pointers -> self-referential
#     `Union{fates_cohort_type,Nothing}` fields (legal for a Julia mutable struct).
#   * Type-bound procedures (Init/NanValues/ZeroValues/Create/Copy/FreeMemory/
#     CanUpperUnder/InitPRTBoundaryConditions/UpdateCohortBioPhysRates/Dump) ->
#     bang-functions dispatching on `fates_cohort_type`, preserving the Fortran
#     names. NaN -> `NaN`; `fates_unset_int` flushes follow Fortran exactly.
#   * The `nan` initialization of scalar fields is preserved as `NaN`.
#   * Module flag globals (hlm_parteh_mode/hlm_use_sp/hlm_use_planthydro/nleafage)
#     are `Ref{Int}` in FatesInterfaceTypesMod -> dereferenced with `[]`.
#   * Param singletons: `prt_params` (live instance), `edpftvarcon_inst()`,
#     `param_derived()` (function accessors).
#
# Deps: FatesConstantsMod (fates_unset_int, ifalse, itrue, nearzero,
# ican_upper, ican_ustory), EDParamsMod (nlevleaf), FatesGlobals (fates_endrun,
# fates_log), PRTGenericMod (AbstractPRTVartypes, organ/element ids,
# max_nleafage, CheckInitialConditions/CopyPRTVartypes!/GetState/
# DeallocatePRTVartypes!, prt_carbon_allom_hyp, prt_cnp_flex_allom_hyp),
# PRTParametersMod (prt_params), FatesParameterDerivedMod (param_derived),
# FatesHydraulicsMemMod (ed_cohort_hydr_type, CopyCohortHydraulics!,
# DeallocateHydrCohortArrays!), FatesInterfaceTypesMod (hlm_* flags, nleafage),
# EDPftvarcon (edpftvarcon_inst), FatesSizeAgeTypeIndicesMod
# (sizetype_class_index, coagetype_class_index), FatesAllometryMod (carea_allom,
# tree_lai, tree_sai), PRTAllometricCarbonMod / PRTAllometricCNPMod (BC ids).

# ===========================================================================
# FATES COHORT TYPE
# ===========================================================================

"""
    fates_cohort_type

The per-cohort state of a FATES demographic unit. Field names are preserved from
the Fortran `fates_cohort_type` for traceability. Scalar reals default to `NaN`
and integers to `fates_unset_int`, matching the Fortran NanValues initialization;
the `prt`/`co_hydr`/`taller`/`shorter` pointers default to `nothing`.

Construct an empty cohort and then call [`Create`](@ref) (or [`Init`](@ref) for a
bare allocation) to populate it.
"""
Base.@kwdef mutable struct fates_cohort_type
    # POINTERS --------------------------------------------------------------
    taller::Union{fates_cohort_type,Nothing}  = nothing  # pointer to next tallest cohort
    shorter::Union{fates_cohort_type,Nothing} = nothing  # pointer to next shorter cohort

    # Multi-species, multi-organ Plant Reactive Transport (PRT) -------------
    prt::Union{AbstractPRTVartypes,Nothing} = nothing
    l2fr::Float64 = NaN   # leaf to fineroot biomass ratio [kg root / kg leaf]

    # VEGETATION STRUCTURE --------------------------------------------------
    pft::Int                       = fates_unset_int  # pft index
    n::Float64                     = NaN              # number of individuals per area [/m2]
    dbh::Float64                   = NaN              # diameter at breast height [cm]
    coage::Float64                 = NaN              # age [years]
    height::Float64                = NaN              # height [m]
    indexnumber::Int               = fates_unset_int  # unique number for each cohort
    canopy_layer::Int              = fates_unset_int  # canopy status [1 = canopy, 2 = understory, ...]
    canopy_layer_yesterday::Float64 = NaN             # recent canopy status (real for fusion safety)
    crowndamage::Int               = fates_unset_int  # crown damage class [1 = undamaged, >1 = damaged]
    g_sb_laweight::Float64         = NaN              # total conductance weighted by leaf area [m/s]*[m2]
    canopy_trim::Float64           = NaN              # fraction of maximum leaf biomass targeted [0-1]
    leaf_cost::Float64             = NaN              # cost to maintain leaves [kgC/m2/year]
    excl_weight::Float64           = NaN              # proportion of cohort demoted each year
    prom_weight::Float64           = NaN              # proportion of cohort promoted each year
    nv::Int                        = fates_unset_int  # number of leaf layers
    status_coh::Int                = fates_unset_int  # growth status [2 = leaves on, 1 = leaves off]
    efleaf_coh::Float64            = NaN              # elongation factor for leaves [fraction]
    effnrt_coh::Float64            = NaN              # elongation factor for fine roots [fraction]
    efstem_coh::Float64            = NaN              # elongation factor for stem [fraction]
    c_area::Float64                = NaN              # areal extent of canopy [m2]
    treelai::Float64               = NaN              # individual LAI [m2 leaf / m2 crown]
    treesai::Float64               = NaN              # individual SAI [m2 stem / m2 crown]
    isnew::Bool                    = false            # flag for a new cohort (no npp/mortality yet)
    size_class::Int                = fates_unset_int  # diameter size bin index
    coage_class::Int               = fates_unset_int  # age bin index
    size_by_pft_class::Int         = fates_unset_int  # joint size-class x PFT index
    coage_by_pft_class::Int        = fates_unset_int  # joint cohort-age-class x PFT index
    size_class_lasttimestep::Int   = fates_unset_int  # size class at the last time step

    # CARBON AND NUTRIENT FLUXES --------------------------------------------
    gpp_tstep::Float64    = NaN   # Gross Primary Production [kgC/indiv/timestep]
    gpp_acc::Float64      = NaN
    gpp_acc_hold::Float64 = NaN

    npp_tstep::Float64    = NaN   # Net Primary Production
    npp_acc::Float64      = NaN
    npp_acc_hold::Float64 = NaN

    resp_tstep::Float64    = NaN  # Autotrophic respiration
    resp_acc::Float64      = NaN
    resp_acc_hold::Float64 = NaN

    c13disc_clm::Float64 = NaN    # carbon-13 discrimination, per indiv/timestep [ppm]
    c13disc_acc::Float64 = NaN    # carbon-13 discrimination, per indiv/day [ppm]

    # Canopy-top, 25degC, leaf-age-weighted biophysical rates
    vcmax25top::Float64 = NaN     # max carboxylation [umol CO2/m2/s]
    jmax25top::Float64  = NaN     # max electron transport [umol electrons/m2/s]
    tpu25top::Float64   = NaN     # triose phosphate utilization [umol CO2/m2/s]
    kp25top::Float64    = NaN     # initial slope of CO2 response (C4) at 25C

    ts_net_uptake::Vector{Float64}   = fill(NaN, nlevleaf)  # net uptake of leaf layers [kgC/m2/timestep]
    year_net_uptake::Vector{Float64} = fill(NaN, nlevleaf)  # net uptake of leaf layers [kgC/m2/year]

    # used for CNP
    cnp_limiter::Int    = fates_unset_int  # limiting element [0 none, 1 C, 2 N, 3 P]
    cx_int::Float64     = NaN  # time integral of log relative C storage over relative nutrient
    ema_dcxdt::Float64  = NaN  # derivative of log relative C storage over relative nutrient
    cx0::Float64        = NaN  # previous-step value of the above log ratio
    nc_repro::Float64   = NaN  # N:C ratio of a new recruit
    pc_repro::Float64   = NaN  # P:C ratio of a new recruit

    # Nutrient fluxes (if N, P, etc. are on)
    daily_nh4_uptake::Float64 = NaN  # daily uptake of mineralized ammonium [kgN/plant/day]
    daily_no3_uptake::Float64 = NaN  # daily uptake of mineralized nitrate  [kgN/plant/day]
    sym_nfix_daily::Float64   = NaN  # accumulated symbiotic N fixation [kgN/indiv/day]
    sym_nfix_tstep::Float64   = NaN  # symbiotic N fixation for the time-step [kgN/indiv/timestep]
    daily_n_gain::Float64     = NaN  # sum of fixation + uptake of mineralized NH4/NO3
    daily_p_gain::Float64     = NaN  # daily uptake of mineralized P [kgP/plant/day]
    daily_c_efflux::Float64   = NaN  # daily efflux of excess C from roots into labile pool [kgC/plant/day]
    daily_n_efflux::Float64   = NaN  # daily efflux of excess N [kgN/plant/day]
    daily_p_efflux::Float64   = NaN  # daily efflux of excess P [kgP/plant/day]
    daily_n_demand::Float64   = NaN  # daily N demanded by the plant [kgN/plant/day]
    daily_p_demand::Float64   = NaN  # daily P demanded by the plant [kgP/plant/day]
    seed_prod::Float64        = NaN  # diagnostic seed production rate [kgC/plant/day]

    twostr_col::Int = fates_unset_int  # column index in the two-stream solution this cohort is part of

    # RESPIRATION COMPONENTS ------------------------------------------------
    rdark::Float64            = NaN  # dark respiration [kgC/indiv/s]
    resp_g_tstep::Float64     = NaN  # growth respiration [kgC/indiv/timestep]
    resp_m::Float64           = NaN  # maintenance respiration [kgC/indiv/timestep]
    resp_m_unreduced::Float64 = NaN  # diagnostic unreduced maintenance respiration
    resp_excess::Float64      = NaN  # respiration of excess carbon [kgC/indiv/day]
    livestem_mr::Float64      = NaN  # aboveground live stem maintenance respiration [kgC/indiv/s]
    livecroot_mr::Float64     = NaN  # belowground live stem maintenance respiration [kgC/indiv/s]
    froot_mr::Float64         = NaN  # live fine root maintenance respiration [kgC/indiv/s]

    # DAMAGE ----------------------------------------------------------------
    branch_frac::Float64 = NaN  # fraction of aboveground woody biomass in branches [0-1]

    # MORTALITY -------------------------------------------------------------
    dmort::Float64  = NaN  # proportional mortality rate [/year]
    bmort::Float64  = NaN  # background mortality rate [indiv/year]
    cmort::Float64  = NaN  # carbon starvation mortality rate [indiv/year]
    hmort::Float64  = NaN  # hydraulic failure mortality rate [indiv/year]
    frmort::Float64 = NaN  # freezing mortality rate [indiv/year]
    smort::Float64  = NaN  # senescence mortality [indiv/year]
    asmort::Float64 = NaN  # age senescence mortality [indiv/year]
    dgmort::Float64 = NaN  # damage mortality [indiv/year]

    # Logging mortality rates
    lmort_direct::Float64     = NaN  # directly logged rate [fraction/logging activity]
    lmort_collateral::Float64 = NaN  # collaterally damaged rate
    lmort_infra::Float64      = NaN  # mechanically damaged rate
    l_degrad::Float64         = NaN  # rate of trees moved to degraded forest

    # GROWTH DERIVATIVES ----------------------------------------------------
    dndt::Float64     = NaN  # time derivative of cohort size [n/year]
    dhdt::Float64     = NaN  # time derivative of height [m/year]
    ddbhdt::Float64   = NaN  # time derivative of dbh [cm/year]
    dbdeaddt::Float64 = NaN  # time derivative of dead biomass [kgC/year]

    # FIRE ------------------------------------------------------------------
    fraction_crown_burned::Float64 = NaN  # proportion of crown affected by fire [0-1]
    cambial_mort::Float64          = NaN  # probability of death from cambial charring [0-1]
    crownfire_mort::Float64        = NaN  # probability of post-fire mortality from crown scorch [0-1]
    fire_mort::Float64             = NaN  # post-fire mortality from cambial+crown damage [0-1]

    # HYDRAULICS ------------------------------------------------------------
    co_hydr::Union{ed_cohort_hydr_type,Nothing} = nothing  # all cohort hydraulics data
end

# ===========================================================================
# Cohort routines (Fortran type-bound procedures)
# ===========================================================================

"""
    Init(this::fates_cohort_type, prt::AbstractPRTVartypes)

Create a new cohort and set default values for all variables. NaNs everything,
zeros what must be zeroed, attaches the (already-allocated and -initialized)
PARTEH `prt` object and verifies its initial conditions, then flags the cohort as
new. Mirrors the Fortran `Init`.
"""
function Init(this::fates_cohort_type, prt::AbstractPRTVartypes)
    NanValues(this)   # make everything in the cohort not-a-number
    ZeroValues(this)  # zero things that need to be zeroed

    # point to the PARTEH object
    this.prt = prt

    # The PARTEH cohort object should be allocated and already initialized.
    CheckInitialConditions(this.prt)

    # New cohorts have no mortality and have moved no carbon; flag accordingly.
    this.isnew = true

    return nothing
end

"""
    NanValues(this::fates_cohort_type)

Make all the cohort variables NaN / unset so they are not used before being
defined. Pointers are nulled (`nothing`). Mirrors the Fortran `NanValues`.
"""
function NanValues(this::fates_cohort_type)
    # set pointers to null
    this.taller  = nothing
    this.shorter = nothing
    this.prt     = nothing
    this.co_hydr = nothing

    # VEGETATION STRUCTURE
    this.l2fr                    = NaN
    this.pft                     = fates_unset_int
    this.n                       = NaN
    this.dbh                     = NaN
    this.coage                   = NaN
    this.height                  = NaN
    this.indexnumber             = fates_unset_int
    this.canopy_layer            = fates_unset_int
    this.canopy_layer_yesterday  = NaN
    this.crowndamage             = fates_unset_int
    this.g_sb_laweight           = NaN
    this.canopy_trim             = NaN
    this.leaf_cost               = NaN
    this.excl_weight             = NaN
    this.prom_weight             = NaN
    this.nv                      = fates_unset_int
    this.status_coh              = fates_unset_int
    this.efleaf_coh              = NaN
    this.effnrt_coh              = NaN
    this.efstem_coh              = NaN
    this.c_area                  = NaN
    this.treelai                 = NaN
    this.treesai                 = NaN
    this.isnew                   = false
    this.size_class              = fates_unset_int
    this.coage_class             = fates_unset_int
    this.size_by_pft_class       = fates_unset_int
    this.coage_by_pft_class      = fates_unset_int
    this.size_class_lasttimestep = fates_unset_int

    # CARBON AND NUTRIENT FLUXES
    this.gpp_tstep        = NaN
    this.gpp_acc          = NaN
    this.gpp_acc_hold     = NaN
    this.npp_tstep        = NaN
    this.npp_acc          = NaN
    this.npp_acc_hold     = NaN
    this.resp_tstep       = NaN
    this.resp_acc         = NaN
    this.resp_acc_hold    = NaN
    this.c13disc_clm      = NaN
    this.c13disc_acc      = NaN
    this.vcmax25top       = NaN
    this.jmax25top        = NaN
    this.tpu25top         = NaN
    this.kp25top          = NaN
    fill!(this.year_net_uptake, NaN)
    fill!(this.ts_net_uptake, NaN)
    this.cnp_limiter      = fates_unset_int
    this.cx_int           = NaN
    this.ema_dcxdt        = NaN
    this.cx0              = NaN
    this.nc_repro         = NaN
    this.pc_repro         = NaN
    this.daily_nh4_uptake = NaN
    this.daily_no3_uptake = NaN
    this.sym_nfix_daily   = NaN
    this.sym_nfix_tstep   = NaN
    this.daily_n_gain     = NaN
    this.daily_p_gain     = NaN
    this.daily_c_efflux   = NaN
    this.daily_n_efflux   = NaN
    this.daily_p_efflux   = NaN
    this.daily_n_demand   = NaN
    this.daily_p_demand   = NaN
    this.seed_prod        = NaN

    # RESPIRATION COMPONENTS
    this.rdark            = NaN
    this.resp_g_tstep     = NaN
    this.resp_m           = NaN
    this.resp_m_unreduced = NaN
    this.resp_excess      = NaN
    this.livestem_mr      = NaN
    this.livecroot_mr     = NaN
    this.froot_mr         = NaN

    # DAMAGE
    this.branch_frac      = NaN

    # MORTALITY
    this.dmort            = NaN
    this.bmort            = NaN
    this.cmort            = NaN
    this.frmort           = NaN
    this.smort            = NaN
    this.asmort           = NaN
    this.dgmort           = NaN
    this.lmort_direct     = NaN
    this.lmort_collateral = NaN
    this.lmort_infra      = NaN
    this.l_degrad         = NaN

    # GROWTH DERIVATIVES
    this.dndt             = NaN
    this.dhdt             = NaN
    this.ddbhdt           = NaN
    this.dbdeaddt         = NaN

    # FIRE
    this.fraction_crown_burned = NaN
    this.cambial_mort          = NaN
    this.crownfire_mort        = NaN
    this.fire_mort             = NaN

    return nothing
end

"""
    ZeroValues(this::fates_cohort_type)

Zero variables that need to be accounted for if this cohort is altered before they
are defined. Mirrors the Fortran `ZeroValues` exactly, including the special
initializations (`year_net_uptake = 999`, `daily_n/p_demand = -9`) and the
deliberately-NOT-zeroed accumulator-hold fields.
"""
function ZeroValues(this::fates_cohort_type)
    this.g_sb_laweight = 0.0

    this.leaf_cost   = 0.0
    this.excl_weight = 0.0
    this.prom_weight = 0.0
    this.nv          = 0
    this.status_coh  = 0
    this.efleaf_coh  = 0.0
    this.effnrt_coh  = 0.0
    this.efstem_coh  = 0.0

    this.treesai     = 0.0
    this.size_class  = 1
    this.coage_class = 1

    this.size_class_lasttimestep = 0
    this.gpp_tstep  = 0.0
    this.gpp_acc    = 0.0
    this.npp_tstep  = 0.0
    this.npp_acc    = 0.0
    this.resp_tstep = 0.0
    this.resp_acc   = 0.0

    # NOTE: gpp_acc_hold / npp_acc_hold / resp_acc_hold are deliberately left
    # unzeroed (kept NaN) to prevent uninitialized use, matching Fortran.

    this.c13disc_clm = 0.0
    this.c13disc_acc = 0.0

    fill!(this.ts_net_uptake, 0.0)
    fill!(this.year_net_uptake, 999.0)  # must be 999, or trimming of new cohorts breaks

    this.daily_nh4_uptake = 0.0
    this.daily_no3_uptake = 0.0

    # Fixation is integrated over the day; zero on creation and after allocation.
    this.sym_nfix_daily = 0.0
    this.daily_n_gain   = 0.0
    this.daily_p_gain   = 0.0

    # Daily nutrient fluxes are integrated; MUST be zeroed on creation AND after
    # allocation. They exist in carbon-only mode but are not used.
    this.daily_c_efflux = 0.0
    this.daily_n_efflux = 0.0
    this.daily_p_efflux = 0.0

    # initialize these as negative
    this.daily_n_demand   = -9.0
    this.daily_p_demand   = -9.0
    this.seed_prod        = 0.0
    this.rdark            = 0.0
    this.resp_g_tstep     = 0.0
    this.resp_m           = 0.0
    this.resp_m_unreduced = 0.0
    this.resp_excess      = 0.0
    this.livestem_mr      = 0.0
    this.livecroot_mr     = 0.0
    this.froot_mr         = 0.0

    this.dmort                 = 0.0
    this.lmort_direct          = 0.0
    this.lmort_collateral      = 0.0
    this.lmort_infra           = 0.0
    this.l_degrad              = 0.0
    this.fraction_crown_burned = 0.0
    this.cambial_mort          = 0.0
    this.crownfire_mort        = 0.0
    this.fire_mort             = 0.0

    return nothing
end

"""
    Create(this, prt, pft, nn, height, coage, dbh, status, ctrim, carea,
           clayer, crowndamage, spread, can_tlai, elongf_leaf, elongf_fnrt,
           elongf_stem)

Set up values for a newly created cohort: initialize it (attaching the PARTEH
`prt` object), copy in the prescribed structural/phenology state, set the
leaf-to-fineroot ratio from the PFT parameter, update the canopy-top biophysical
rates, classify the cohort into size/age bins, compute crown area, tree LAI and
tree SAI, and register the PARTEH boundary conditions. Mirrors the Fortran
`Create`. `can_tlai` is the patch-level total LAI of each canopy layer.
"""
function Create(this::fates_cohort_type, prt::AbstractPRTVartypes, pft::Integer,
                nn::Real, height::Real, coage::Real, dbh::Real, status::Integer,
                ctrim::Real, carea::Real, clayer::Integer, crowndamage::Integer,
                spread::Real, can_tlai::AbstractVector, elongf_leaf::Real,
                elongf_fnrt::Real, elongf_stem::Real)

    # initialize cohort
    Init(this, prt)

    # set values
    this.pft          = pft
    this.crowndamage  = crowndamage
    this.canopy_layer = clayer
    this.canopy_layer_yesterday = Float64(clayer)
    this.status_coh   = status
    this.n            = nn
    this.height       = height
    this.dbh          = dbh
    this.coage        = coage
    this.canopy_trim  = ctrim
    this.efleaf_coh   = elongf_leaf
    this.effnrt_coh   = elongf_fnrt
    this.efstem_coh   = elongf_stem

    # This routine may be called during restarts, when the actual cohort data is
    # unknown. Nominal values for size/number/type are passed in those cases.
    if this.dbh <= 0.0 || this.n == 0.0 || this.pft == 0
        fates_endrun("FATES: something is zero in cohort Create: " *
                     "dbh=$(this.dbh) n=$(this.n) pft=$(this.pft)")
    end

    # Initialize the leaf-to-fineroot biomass ratio. For C-only this stays
    # constant; for nutrient-enabled it is dynamic. New cohorts start at the
    # minimum (the allom_l2fr parameter).
    this.l2fr = prt_params.allom_l2fr[pft]

    if hlm_parteh_mode[] == prt_cnp_flex_allom_hyp
        this.cx_int      = 0.0  # Assume balanced N,P/C stores, log(1) = 0
        this.cx0         = 0.0  # Assume balanced N,P/C stores, log(1) = 0
        this.ema_dcxdt   = 0.0  # Assume unchanged dCX/dt
        this.cnp_limiter = 0    # Assume limitations are unknown
    end

    # Sets vcmax25top etc., which depend on the leaf age fractions (from PARTEH).
    UpdateCohortBioPhysRates(this)

    # calculate size classes
    this.size_class, this.size_by_pft_class = sizetype_class_index(this.dbh, this.pft)

    # If cohort age tracking is off we call this once so everything is in the
    # first bin (makes copy/terminate easier later).
    this.coage_class, this.coage_by_pft_class = coagetype_class_index(this.coage, this.pft)

    # assign or calculate canopy extent and depth
    if hlm_use_sp[] == ifalse
        _, this.c_area = carea_allom(this.dbh, this.n, spread, this.pft, this.crowndamage)
    else
        # set this from previously precision-controlled value in SP mode
        this.c_area = carea
    end

    # Query PARTEH for the leaf carbon [kg]
    leaf_c = GetState(this.prt, leaf_organ, carbon12_element)

    # calculate tree lai
    this.treelai = tree_lai(leaf_c, this.pft, this.c_area, this.n,
                            this.canopy_layer, can_tlai, this.vcmax25top)

    if hlm_use_sp[] == ifalse
        this.treesai = tree_sai(this.pft, this.dbh, this.crowndamage,
                                this.canopy_trim, this.efstem_coh, this.c_area, this.n,
                                this.canopy_layer, can_tlai, this.treelai,
                                this.vcmax25top, 2)
    end

    InitPRTBoundaryConditions(this)

    return nothing
end

"""
    Copy(this::fates_cohort_type, copyCohort::fates_cohort_type)

Copy all the variables in one cohort (`this`) into another (`copyCohort`). The
PRT state is copied via [`CopyPRTVartypes!`](@ref); the `taller`/`shorter`
linked-list pointers are nulled and `indexnumber` is unset (the copy is a fresh
node). The CNP-only fields (`cx_int`/`ema_dcxdt`/`cx0`) are copied only under the
flexible-CNP hypothesis, and the hydraulics state only when plant hydraulics is
on. Mirrors the Fortran `Copy`.
"""
function Copy(this::fates_cohort_type, copyCohort::fates_cohort_type)
    copyCohort.indexnumber = fates_unset_int

    # POINTERS
    copyCohort.taller  = nothing
    copyCohort.shorter = nothing

    # PRT
    CopyPRTVartypes!(copyCohort.prt, this.prt)
    copyCohort.l2fr = this.l2fr

    # VEGETATION STRUCTURE
    copyCohort.pft                     = this.pft
    copyCohort.n                       = this.n
    copyCohort.dbh                     = this.dbh
    copyCohort.coage                   = this.coage
    copyCohort.height                  = this.height
    copyCohort.canopy_layer            = this.canopy_layer
    copyCohort.canopy_layer_yesterday  = this.canopy_layer_yesterday
    copyCohort.crowndamage             = this.crowndamage
    copyCohort.g_sb_laweight           = this.g_sb_laweight
    copyCohort.canopy_trim             = this.canopy_trim
    copyCohort.leaf_cost               = this.leaf_cost
    copyCohort.excl_weight             = this.excl_weight
    copyCohort.prom_weight             = this.prom_weight
    copyCohort.nv                      = this.nv
    copyCohort.status_coh              = this.status_coh
    copyCohort.efleaf_coh              = this.efleaf_coh
    copyCohort.effnrt_coh              = this.effnrt_coh
    copyCohort.efstem_coh              = this.efstem_coh
    copyCohort.c_area                  = this.c_area
    copyCohort.treelai                 = this.treelai
    copyCohort.treesai                 = this.treesai
    copyCohort.isnew                   = this.isnew
    copyCohort.size_class              = this.size_class
    copyCohort.coage_class             = this.coage_class
    copyCohort.size_by_pft_class       = this.size_by_pft_class
    copyCohort.coage_by_pft_class      = this.coage_by_pft_class
    copyCohort.size_class_lasttimestep = this.size_class_lasttimestep

    # CARBON AND NUTRIENT FLUXES
    copyCohort.gpp_tstep      = this.gpp_tstep
    copyCohort.gpp_acc        = this.gpp_acc
    copyCohort.gpp_acc_hold   = this.gpp_acc_hold
    copyCohort.npp_tstep      = this.npp_tstep
    copyCohort.npp_acc        = this.npp_acc
    copyCohort.npp_acc_hold   = this.npp_acc_hold
    copyCohort.resp_tstep     = this.resp_tstep
    copyCohort.resp_acc       = this.resp_acc
    copyCohort.resp_acc_hold  = this.resp_acc_hold
    copyCohort.c13disc_clm    = this.c13disc_clm
    copyCohort.c13disc_acc    = this.c13disc_acc
    copyCohort.vcmax25top     = this.vcmax25top
    copyCohort.jmax25top      = this.jmax25top
    copyCohort.tpu25top       = this.tpu25top
    copyCohort.kp25top        = this.kp25top
    copyCohort.ts_net_uptake   .= this.ts_net_uptake
    copyCohort.year_net_uptake .= this.year_net_uptake
    copyCohort.cnp_limiter    = this.cnp_limiter

    if hlm_parteh_mode[] == prt_cnp_flex_allom_hyp
        copyCohort.cx_int    = this.cx_int
        copyCohort.ema_dcxdt = this.ema_dcxdt
        copyCohort.cx0       = this.cx0
    end

    copyCohort.nc_repro         = this.nc_repro
    copyCohort.daily_nh4_uptake = this.daily_nh4_uptake
    copyCohort.daily_no3_uptake = this.daily_no3_uptake
    copyCohort.sym_nfix_daily   = this.sym_nfix_daily
    copyCohort.sym_nfix_tstep   = this.sym_nfix_tstep
    copyCohort.daily_n_gain     = this.daily_n_gain
    copyCohort.daily_p_gain     = this.daily_p_gain
    copyCohort.daily_c_efflux   = this.daily_c_efflux
    copyCohort.daily_n_efflux   = this.daily_n_efflux
    copyCohort.daily_p_efflux   = this.daily_p_efflux
    copyCohort.daily_n_demand   = this.daily_n_demand
    copyCohort.daily_p_demand   = this.daily_p_demand
    copyCohort.seed_prod        = this.seed_prod

    # RESPIRATION COMPONENTS
    copyCohort.rdark            = this.rdark
    copyCohort.resp_g_tstep     = this.resp_g_tstep
    copyCohort.resp_m           = this.resp_m
    copyCohort.resp_m_unreduced = this.resp_m_unreduced
    copyCohort.resp_excess      = this.resp_excess
    copyCohort.livestem_mr      = this.livestem_mr
    copyCohort.livecroot_mr     = this.livecroot_mr
    copyCohort.froot_mr         = this.froot_mr

    # DAMAGE
    copyCohort.branch_frac = this.branch_frac

    # MORTALITY
    copyCohort.dmort            = this.dmort
    copyCohort.bmort            = this.bmort
    copyCohort.cmort            = this.cmort
    copyCohort.hmort            = this.hmort
    copyCohort.frmort           = this.frmort
    copyCohort.smort            = this.smort
    copyCohort.asmort           = this.asmort
    copyCohort.dgmort           = this.dgmort
    copyCohort.lmort_direct     = this.lmort_direct
    copyCohort.lmort_collateral = this.lmort_collateral
    copyCohort.lmort_infra      = this.lmort_infra
    copyCohort.l_degrad         = this.l_degrad

    # GROWTH DERIVATIVES
    copyCohort.dndt     = this.dndt
    copyCohort.dhdt     = this.dhdt
    copyCohort.ddbhdt   = this.ddbhdt
    copyCohort.dbdeaddt = this.dbdeaddt

    # FIRE
    copyCohort.fraction_crown_burned = this.fraction_crown_burned
    copyCohort.cambial_mort          = this.cambial_mort
    copyCohort.crownfire_mort        = this.crownfire_mort
    copyCohort.fire_mort             = this.fire_mort

    # HYDRAULICS
    if hlm_use_planthydro[] == itrue
        CopyCohortHydraulics!(copyCohort.co_hydr, this.co_hydr)
    end

    return nothing
end

"""
    FreeMemory(this::fates_cohort_type)

Release the dynamic memory and objects held within the cohort structure (the
hydraulics object's per-layer arrays and the PRT object's state) without
destroying the cohort itself. Julia is garbage-collected, so this drops the
references (sets `co_hydr`/`prt` to `nothing`) after delegating to the contained
objects' deallocators. Mirrors the Fortran `FreeMemory`.
"""
function FreeMemory(this::fates_cohort_type)
    # at this point, nothing should be pointing to current cohort
    if hlm_use_planthydro[] == itrue
        if this.co_hydr !== nothing
            DeallocateHydrCohortArrays!(this.co_hydr)
        end
        this.co_hydr = nothing
    end

    # deallocate the cohort's PRT structures
    if this.prt !== nothing
        DeallocatePRTVartypes!(this.prt)
    end
    this.prt = nothing

    return nothing
end

"""
    InitPRTBoundaryConditions(this::fates_cohort_type)

Set the boundary conditions that flow in and out of the PARTEH allocation
hypotheses. Each `RegisterBC*` call wires a `Ref` into the cohort's scalar state
so the allocation solver can read/update it. The set of registered conditions
depends on `hlm_parteh_mode`. Mirrors the Fortran `InitPRTBoundaryConditions`.

Because Julia scalars are passed by value, each registered boundary condition is
backed by a `Ref` to the corresponding cohort field; after a PARTEH solve the
caller must copy any updated inout/out `Ref` values back into the cohort. The
`Ref`s are stored on the PRT object's `bc_*` arrays.
"""
function InitPRTBoundaryConditions(this::fates_cohort_type)
    mode = hlm_parteh_mode[]

    if mode == prt_carbon_allom_hyp
        # Register boundary conditions for the Carbon-Only Allometric Hypothesis
        RegisterBCInout!(this.prt, ac_bc_inout_id_dbh;   bc_rval = Ref(this.dbh))
        RegisterBCInout!(this.prt, ac_bc_inout_id_netdc; bc_rval = Ref(this.npp_acc))
        RegisterBCIn!(this.prt, ac_bc_in_id_cdamage; bc_ival = Ref(this.crowndamage))
        RegisterBCIn!(this.prt, ac_bc_in_id_pft;     bc_ival = Ref(this.pft))
        RegisterBCIn!(this.prt, ac_bc_in_id_ctrim;   bc_rval = Ref(this.canopy_trim))
        RegisterBCIn!(this.prt, ac_bc_in_id_lstat;   bc_ival = Ref(this.status_coh))
        RegisterBCIn!(this.prt, ac_bc_in_id_efleaf;  bc_rval = Ref(this.efleaf_coh))
        RegisterBCIn!(this.prt, ac_bc_in_id_effnrt;  bc_rval = Ref(this.effnrt_coh))
        RegisterBCIn!(this.prt, ac_bc_in_id_efstem;  bc_rval = Ref(this.efstem_coh))

    elseif mode == prt_cnp_flex_allom_hyp
        # Register boundary conditions for the CNP Allometric Hypothesis
        RegisterBCIn!(this.prt, acnp_bc_in_id_pft;     bc_ival = Ref(this.pft))
        RegisterBCIn!(this.prt, acnp_bc_in_id_ctrim;   bc_rval = Ref(this.canopy_trim))
        RegisterBCIn!(this.prt, acnp_bc_in_id_lstat;   bc_ival = Ref(this.status_coh))
        RegisterBCIn!(this.prt, acnp_bc_in_id_efleaf;  bc_rval = Ref(this.efleaf_coh))
        RegisterBCIn!(this.prt, acnp_bc_in_id_effnrt;  bc_rval = Ref(this.effnrt_coh))
        RegisterBCIn!(this.prt, acnp_bc_in_id_efstem;  bc_rval = Ref(this.efstem_coh))
        RegisterBCIn!(this.prt, acnp_bc_in_id_netdc;   bc_rval = Ref(this.npp_acc))

        RegisterBCIn!(this.prt, acnp_bc_in_id_nc_repro; bc_rval = Ref(this.nc_repro))
        RegisterBCIn!(this.prt, acnp_bc_in_id_pc_repro; bc_rval = Ref(this.pc_repro))
        RegisterBCIn!(this.prt, acnp_bc_in_id_cdamage;  bc_ival = Ref(this.crowndamage))

        RegisterBCInout!(this.prt, acnp_bc_inout_id_dbh;         bc_rval = Ref(this.dbh))
        RegisterBCInout!(this.prt, acnp_bc_inout_id_resp_excess; bc_rval = Ref(this.resp_excess))
        RegisterBCInout!(this.prt, acnp_bc_inout_id_l2fr;        bc_rval = Ref(this.l2fr))
        RegisterBCInout!(this.prt, acnp_bc_inout_id_cx_int;      bc_rval = Ref(this.cx_int))
        RegisterBCInout!(this.prt, acnp_bc_inout_id_emadcxdt;    bc_rval = Ref(this.ema_dcxdt))
        RegisterBCInout!(this.prt, acnp_bc_inout_id_cx0;         bc_rval = Ref(this.cx0))

        RegisterBCInout!(this.prt, acnp_bc_inout_id_netdn; bc_rval = Ref(this.daily_n_gain))
        RegisterBCInout!(this.prt, acnp_bc_inout_id_netdp; bc_rval = Ref(this.daily_p_gain))

        RegisterBCOut!(this.prt, acnp_bc_out_id_cefflux; bc_rval = Ref(this.daily_c_efflux))
        RegisterBCOut!(this.prt, acnp_bc_out_id_nefflux; bc_rval = Ref(this.daily_n_efflux))
        RegisterBCOut!(this.prt, acnp_bc_out_id_pefflux; bc_rval = Ref(this.daily_p_efflux))
        RegisterBCOut!(this.prt, acnp_bc_out_id_limiter; bc_ival = Ref(this.cnp_limiter))

    else
        fates_endrun("You specified an unknown PRT module. Aborting.")
    end

    return nothing
end

"""
    sync_cohort_to_prt_bcs!(this)   /   sync_prt_bcs_to_cohort!(this)

Refresh the PARTEH boundary-condition `Ref`s from the live cohort fields before an
allocation/turnover solve, and copy the updated inout/out `Ref`s back afterward.
[`InitPRTBoundaryConditions`](@ref) registers every BC as a `Ref(field)` captured ONCE
at cohort creation — the Julia stand-in for the Fortran pointer that ALIASES the field.
Nothing else refreshes them, so without this the solve reads stale cold-start inputs
(most importantly `netdc`/`npp_acc` ≡ 0 → zero allocation → no growth) and its outputs
(grown `dbh`, effluxes) never reach the cohort. These two functions mirror the exact
BC↔field map in `InitPRTBoundaryConditions` for each PARTEH hypothesis: push covers all
`bc_in` + `bc_inout`; pull covers all `bc_inout` + `bc_out`.
"""
function sync_cohort_to_prt_bcs!(this::fates_cohort_type)
    prt = this.prt
    mode = hlm_parteh_mode[]
    if mode == prt_carbon_allom_hyp
        prt.bc_inout[ac_bc_inout_id_dbh].rval[]   = this.dbh
        prt.bc_inout[ac_bc_inout_id_netdc].rval[] = this.npp_acc
        prt.bc_in[ac_bc_in_id_cdamage].ival[] = this.crowndamage
        prt.bc_in[ac_bc_in_id_pft].ival[]     = this.pft
        prt.bc_in[ac_bc_in_id_ctrim].rval[]   = this.canopy_trim
        prt.bc_in[ac_bc_in_id_lstat].ival[]   = this.status_coh
        prt.bc_in[ac_bc_in_id_efleaf].rval[]  = this.efleaf_coh
        prt.bc_in[ac_bc_in_id_effnrt].rval[]  = this.effnrt_coh
        prt.bc_in[ac_bc_in_id_efstem].rval[]  = this.efstem_coh
    elseif mode == prt_cnp_flex_allom_hyp
        prt.bc_in[acnp_bc_in_id_pft].ival[]      = this.pft
        prt.bc_in[acnp_bc_in_id_ctrim].rval[]    = this.canopy_trim
        prt.bc_in[acnp_bc_in_id_lstat].ival[]    = this.status_coh
        prt.bc_in[acnp_bc_in_id_efleaf].rval[]   = this.efleaf_coh
        prt.bc_in[acnp_bc_in_id_effnrt].rval[]   = this.effnrt_coh
        prt.bc_in[acnp_bc_in_id_efstem].rval[]   = this.efstem_coh
        prt.bc_in[acnp_bc_in_id_netdc].rval[]    = this.npp_acc
        prt.bc_in[acnp_bc_in_id_nc_repro].rval[] = this.nc_repro
        prt.bc_in[acnp_bc_in_id_pc_repro].rval[] = this.pc_repro
        prt.bc_in[acnp_bc_in_id_cdamage].ival[]  = this.crowndamage
        prt.bc_inout[acnp_bc_inout_id_dbh].rval[]         = this.dbh
        prt.bc_inout[acnp_bc_inout_id_resp_excess].rval[] = this.resp_excess
        prt.bc_inout[acnp_bc_inout_id_l2fr].rval[]        = this.l2fr
        prt.bc_inout[acnp_bc_inout_id_cx_int].rval[]      = this.cx_int
        prt.bc_inout[acnp_bc_inout_id_emadcxdt].rval[]    = this.ema_dcxdt
        prt.bc_inout[acnp_bc_inout_id_cx0].rval[]         = this.cx0
        prt.bc_inout[acnp_bc_inout_id_netdn].rval[]       = this.daily_n_gain
        prt.bc_inout[acnp_bc_inout_id_netdp].rval[]       = this.daily_p_gain
    end
    return nothing
end

function sync_prt_bcs_to_cohort!(this::fates_cohort_type)
    prt = this.prt
    mode = hlm_parteh_mode[]
    if mode == prt_carbon_allom_hyp
        this.dbh     = prt.bc_inout[ac_bc_inout_id_dbh].rval[]
        this.npp_acc = prt.bc_inout[ac_bc_inout_id_netdc].rval[]
    elseif mode == prt_cnp_flex_allom_hyp
        this.dbh          = prt.bc_inout[acnp_bc_inout_id_dbh].rval[]
        this.resp_excess  = prt.bc_inout[acnp_bc_inout_id_resp_excess].rval[]
        this.l2fr         = prt.bc_inout[acnp_bc_inout_id_l2fr].rval[]
        this.cx_int       = prt.bc_inout[acnp_bc_inout_id_cx_int].rval[]
        this.ema_dcxdt    = prt.bc_inout[acnp_bc_inout_id_emadcxdt].rval[]
        this.cx0          = prt.bc_inout[acnp_bc_inout_id_cx0].rval[]
        this.daily_n_gain = prt.bc_inout[acnp_bc_inout_id_netdn].rval[]
        this.daily_p_gain = prt.bc_inout[acnp_bc_inout_id_netdp].rval[]
        this.daily_c_efflux = prt.bc_out[acnp_bc_out_id_cefflux].rval[]
        this.daily_n_efflux = prt.bc_out[acnp_bc_out_id_nefflux].rval[]
        this.daily_p_efflux = prt.bc_out[acnp_bc_out_id_pefflux].rval[]
        this.cnp_limiter    = prt.bc_out[acnp_bc_out_id_limiter].ival[]
    end
    return nothing
end

"""
    UpdateCohortBioPhysRates(this::fates_cohort_type)

Update the four canopy-top biophysical rates (`vcmax25top`, `jmax25top`,
`tpu25top`, `kp25top`) of a cohort based on the leaf-age mass fractions held in
its PARTEH object. With leaves present (and not in satellite-phenology mode) the
rates are the leaf-age-weighted PFT values; in SP mode they take the first leaf
age class; otherwise they are zeroed. Mirrors the Fortran
`UpdateCohortBioPhysRates`.
"""
function UpdateCohortBioPhysRates(this::fates_cohort_type)
    nla = nleafage[]
    frac_leaf_aclass = zeros(Float64, max_nleafage)

    # Fraction of leaves in each age class (same proportion across leaf layers).
    for iage in 1:nla
        frac_leaf_aclass[iage] = GetState(this.prt, leaf_organ, carbon12_element, iage)
    end

    ipft = this.pft
    p  = edpftvarcon_inst()
    pd = param_derived()

    if sum(@view frac_leaf_aclass[1:nla]) > nearzero && hlm_use_sp[] == ifalse
        s = sum(@view frac_leaf_aclass[1:nla])
        @views frac_leaf_aclass[1:nla] .= frac_leaf_aclass[1:nla] ./ s

        this.vcmax25top = sum(@views p.vcmax25top[ipft, 1:nla] .* frac_leaf_aclass[1:nla])
        this.jmax25top  = sum(@views pd.jmax25top[ipft, 1:nla] .* frac_leaf_aclass[1:nla])
        this.tpu25top   = sum(@views pd.tpu25top[ipft, 1:nla]  .* frac_leaf_aclass[1:nla])
        this.kp25top    = sum(@views pd.kp25top[ipft, 1:nla]   .* frac_leaf_aclass[1:nla])

    elseif hlm_use_sp[] == itrue
        this.vcmax25top = p.vcmax25top[ipft, 1]
        this.jmax25top  = pd.jmax25top[ipft, 1]
        this.tpu25top   = pd.tpu25top[ipft, 1]
        this.kp25top    = pd.kp25top[ipft, 1]

    else
        this.vcmax25top = 0.0
        this.jmax25top  = 0.0
        this.tpu25top   = 0.0
        this.kp25top    = 0.0
    end

    return nothing
end

"""
    CanUpperUnder(this::fates_cohort_type) -> Int

Determine whether a cohort's crown position is in the upper canopy
(`ican_upper`, when `canopy_layer == 1`) or the understory (`ican_ustory`). Used
only for diagnostics. Mirrors the Fortran `CanUpperUnder`.
"""
function CanUpperUnder(this::fates_cohort_type)
    return this.canopy_layer == 1 ? ican_upper : ican_ustory
end

"""
    Dump(this::fates_cohort_type)

Print out the attributes of a cohort (and its hydraulics object, if present) for
diagnostics. Mirrors the Fortran `Dump`.
"""
function Dump(this::fates_cohort_type)
    log = fates_log()
    println(log, "----------------------------------------")
    println(log, " Dumping Cohort Information             ")
    println(log, "----------------------------------------")
    println(log, "cohort%pft                    = ", this.pft)
    println(log, "cohort%n                      = ", this.n)
    println(log, "cohort%dbh                    = ", this.dbh)
    println(log, "cohort%height                 = ", this.height)
    println(log, "cohort%crowndamage            = ", this.crowndamage)
    println(log, "cohort%coage                  = ", this.coage)
    println(log, "cohort%l2fr                   = ", this.l2fr)
    if this.prt !== nothing
        println(log, "leaf carbon                   = ", GetState(this.prt, leaf_organ, carbon12_element))
        println(log, "fineroot carbon               = ", GetState(this.prt, fnrt_organ, carbon12_element))
        println(log, "sapwood carbon                = ", GetState(this.prt, sapw_organ, carbon12_element))
        println(log, "structural (dead) carbon      = ", GetState(this.prt, struct_organ, carbon12_element))
        println(log, "storage carbon                = ", GetState(this.prt, store_organ, carbon12_element))
        println(log, "reproductive carbon           = ", GetState(this.prt, repro_organ, carbon12_element))
    end
    println(log, "cohort%g_sb_laweight          = ", this.g_sb_laweight)
    println(log, "cohort%leaf_cost              = ", this.leaf_cost)
    println(log, "cohort%canopy_layer           = ", this.canopy_layer)
    println(log, "cohort%canopy_layer_yesterday = ", this.canopy_layer_yesterday)
    println(log, "cohort%nv                     = ", this.nv)
    println(log, "cohort%status_coh             = ", this.status_coh)
    println(log, "co%efleaf_coh                 = ", this.efleaf_coh)
    println(log, "co%effnrt_coh                 = ", this.effnrt_coh)
    println(log, "co%efstem_coh                 = ", this.efstem_coh)
    println(log, "cohort%canopy_trim            = ", this.canopy_trim)
    println(log, "cohort%excl_weight            = ", this.excl_weight)
    println(log, "cohort%prom_weight            = ", this.prom_weight)
    println(log, "cohort%size_class             = ", this.size_class)
    println(log, "cohort%size_by_pft_class      = ", this.size_by_pft_class)
    println(log, "cohort%coage_class            = ", this.coage_class)
    println(log, "cohort%coage_by_pft_class     = ", this.coage_by_pft_class)
    println(log, "cohort%gpp_acc_hold           = ", this.gpp_acc_hold)
    println(log, "cohort%gpp_acc                = ", this.gpp_acc)
    println(log, "cohort%gpp_tstep              = ", this.gpp_tstep)
    println(log, "cohort%npp_acc_hold           = ", this.npp_acc_hold)
    println(log, "cohort%npp_tstep              = ", this.npp_tstep)
    println(log, "cohort%npp_acc                = ", this.npp_acc)
    println(log, "cohort%resp_tstep             = ", this.resp_tstep)
    println(log, "cohort%resp_acc               = ", this.resp_acc)
    println(log, "cohort%resp_acc_hold          = ", this.resp_acc_hold)
    println(log, "cohort%rdark                  = ", this.rdark)
    println(log, "cohort%resp_m                 = ", this.resp_m)
    println(log, "cohort%resp_g_tstep           = ", this.resp_g_tstep)
    println(log, "cohort%livestem_mr            = ", this.livestem_mr)
    println(log, "cohort%livecroot_mr           = ", this.livecroot_mr)
    println(log, "cohort%froot_mr               = ", this.froot_mr)
    println(log, "cohort%treelai                = ", this.treelai)
    println(log, "cohort%treesai                = ", this.treesai)
    println(log, "cohort%c_area                 = ", this.c_area)
    println(log, "cohort%cmort                  = ", this.cmort)
    println(log, "cohort%bmort                  = ", this.bmort)
    println(log, "cohort%smort                  = ", this.smort)
    println(log, "cohort%asmort                 = ", this.asmort)
    println(log, "cohort%dgmort                 = ", this.dgmort)
    println(log, "cohort%hmort                  = ", this.hmort)
    println(log, "cohort%frmort                 = ", this.frmort)
    println(log, "cohort%lmort_direct           = ", this.lmort_direct)
    println(log, "cohort%lmort_collateral       = ", this.lmort_collateral)
    println(log, "cohort%lmort_infra            = ", this.lmort_infra)
    println(log, "cohort%isnew                  = ", this.isnew)
    println(log, "cohort%dndt                   = ", this.dndt)
    println(log, "cohort%dhdt                   = ", this.dhdt)
    println(log, "cohort%ddbhdt                 = ", this.ddbhdt)
    println(log, "cohort%dbdeaddt               = ", this.dbdeaddt)
    println(log, "cohort%fraction_crown_burned  = ", this.fraction_crown_burned)
    println(log, "cohort%fire_mort              = ", this.fire_mort)
    println(log, "cohort%crownfire_mort         = ", this.crownfire_mort)
    println(log, "cohort%cambial_mort           = ", this.cambial_mort)

    # The Fortran calls this%co_hydr%Dump(); the hydraulics memory module
    # (FatesHydraulicsMemMod) does not port a per-cohort Dump method, so we just
    # note its presence.
    this.co_hydr !== nothing && println(log, "cohort has an attached co_hydr object")

    println(log, "----------------------------------------")
    return nothing
end
