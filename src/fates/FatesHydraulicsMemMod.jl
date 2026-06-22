# FatesHydraulicsMemMod.jl
# Julia port of FATES src/fates/main/FatesHydraulicsMemMod.F90
#
# The plant-hydraulics MEMORY / state types. Two types are defined:
#   * ed_site_hydr_type  — site-level hydraulics state (rhizosphere shells, the
#     full soil->leaf connection/node map for the matrix solver, diagnostics, and
#     references to the per-soil-layer water-retention/conductance functions).
#   * ed_cohort_hydr_type — per-cohort water potentials, conductances, volumes,
#     and capacitances by plant compartment (leaf / stem / troot / aroot) and by
#     rhizosphere layer.
#
# Translation notes (per project conventions):
#   * fates_r8 -> Float64.
#   * Fortran allocatables -> SoA Julia arrays, with undef-0 defaults so a freshly
#     constructed struct is cheap; the type-bound `allocate`/`init` routines become
#     bang-functions that size the arrays from the compartment / shell dimensions.
#   * The compartment-index `const`s (n_hypool_leaf, ... nshell, the porous-media
#     indices) are Julia `const`. The Fortran `nan` initialization is preserved as
#     `NaN`; `fates_unset_r8` flushes follow the Fortran exactly.
#   * The Fortran `wrf_arr_type`/`wkf_arr_type` are single-pointer holder types that
#     only exist because Fortran cannot make an array of a polymorphic class.  Julia
#     arrays CAN hold abstract-typed elements directly, so `wrf_soil`/`wkf_soil`
#     become `Vector{WRFType}` / `Vector{WKFType}` (from FatesHydroWTFMod, Batch 1).
#
# Deps: FatesConstantsMod (fates_unset_r8), FatesGlobals (fates_endrun),
#       FatesHydroWTFMod (WRFType / WKFType).

# ---------------------------------------------------------------------------
# Solver options for hydraulics (Fortran integer parameters)
# ---------------------------------------------------------------------------
const hydr_solver_1DTaylor = 1
const hydr_solver_2DPicard = 2
const hydr_solver_2DNewton = 3

# Number of soil layers for indexing cohort fine-root quantities. NOTE (Fortran):
# do not set to 1 except for development/testing.
const nlevsoi_hyd_max = 40

# ---------------------------------------------------------------------------
# Number of distinct plant porous-media compartments (leaf, stem, troot, aroot)
# and the per-compartment node counts. n_hypool_troot/aroot CANNOT be changed.
# ---------------------------------------------------------------------------
const n_porous_media = 5
const n_plant_media  = 4
const n_hypool_leaf  = 1
const n_hypool_stem  = 1
const n_hypool_troot = 1   # CANNOT BE CHANGED
const n_hypool_aroot = 1   # THIS IS "PER-SOIL-LAYER"
const nshell         = 1

# Number of aboveground plant water-storage nodes
const n_hypool_ag = n_hypool_leaf + n_hypool_stem

# Total number of water-storage nodes (single-layer reference count) and the
# plant-only subset (excludes the rhizosphere shells).
const n_hypool_tot   = n_hypool_ag + n_hypool_troot + n_hypool_aroot + nshell
const n_hypool_plant = n_hypool_tot - nshell

# ---------------------------------------------------------------------------
# Porous-medium type indices (used to index over an arbitrary number of pools)
# ---------------------------------------------------------------------------
const stomata_p_media = 0
const leaf_p_media    = 1
const stem_p_media    = 2
const troot_p_media   = 3
const aroot_p_media   = 4
const rhiz_p_media    = 5

# P-V curve: total RWC at which elastic drainage begins (tfs)         [-]
const rwcft  = (1.0, 0.958, 0.958, 0.958)
# P-V curve: total RWC at which capillary reserves are exhausted (tfs) [-]
const rwccap = (1.0, 0.947, 0.947, 0.947)

# Mean fine-root radius expected in the bulk soil [m]
const fine_root_radius_const = 0.0001

# ===========================================================================
# Site-level plant-hydraulics state
# ===========================================================================

"""
    ed_site_hydr_type

Site-level plant-hydraulics state. Holds the rhizosphere-shell geometry and
hydraulic conductances, the absorbing-root lengths, water-balance diagnostics,
the full soil->leaf node/connection map for the matrix (2D) solver, and
references to the per-soil-layer water-retention (`wrf_soil`) and water-
conductance (`wkf_soil`) functions.

Mirrors the Fortran `ed_site_hydr_type`. Allocatable arrays default to undef-0;
call [`InitHydrSite!`](@ref) to size and initialize them.
"""
Base.@kwdef mutable struct ed_site_hydr_type
    # Plant hydraulics geometry --------------------------------------------
    nlevrhiz::Int = 0                                    # Number of rhizosphere levels (vertical layers)
    map_s2r::Vector{Int} = Int[]                         # soil -> rhizosphere level mapping
    map_r2s::Matrix{Int} = Matrix{Int}(undef, 0, 0)      # rhizosphere -> soil level (1: top, 2: bottom)
    zi_rhiz::Vector{Float64} = Float64[]                 # depth of bottom edge of each rhiz level [m]
    dz_rhiz::Vector{Float64} = Float64[]                 # width of each rhiz level [m]

    v_shell::Matrix{Float64}            = Matrix{Float64}(undef, 0, 0)  # rhiz shell volume (m3) over the site
    v_shell_init::Matrix{Float64}       = Matrix{Float64}(undef, 0, 0)  # previous shell volume (m3)
    r_node_shell::Matrix{Float64}       = Matrix{Float64}(undef, 0, 0)  # nodal radius of rhiz shell (m)
    r_node_shell_init::Matrix{Float64}  = Matrix{Float64}(undef, 0, 0)  # previous nodal radius (m)
    l_aroot_layer::Vector{Float64}      = Float64[]                     # total absorbing-root length by soil layer (m)
    l_aroot_layer_init::Vector{Float64} = Float64[]                     # previous absorbing-root length (m)
    kmax_upper_shell::Matrix{Float64}   = Matrix{Float64}(undef, 0, 0)  # max soil k to upper shell boundary (kg s-1 MPa-1)
    kmax_lower_shell::Matrix{Float64}   = Matrix{Float64}(undef, 0, 0)  # max soil k to lower shell boundary (kg s-1 MPa-1)
    r_out_shell::Matrix{Float64}        = Matrix{Float64}(undef, 0, 0)  # outer radius of rhiz shell (m)

    rs1::Vector{Float64} = Float64[]                     # mean fine-root radius (m) (currently constant)

    supsub_flag::Vector{Int} = Int[]                     # index of outermost shell at super-/sub-saturation
    h2osoi_liqvol_shell::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # volumetric water in shell (m3/m3)

    recruit_w_uptake::Vector{Float64} = Float64[]        # recruitment water uptake (kg H2O/m2/s)

    # Water-balance scalars -------------------------------------------------
    errh2o_hyd::Float64        = 0.0    # plant-hydraulics error summed to column level (mm)
    dwat_veg::Float64          = 0.0    # change in stored vegetation water, column level (kg)
    h2oveg::Float64            = 0.0    # stored water in vegetation (kg/m2)
    h2oveg_recruit::Float64    = 0.0    # stored water in recruits (kg/m2)
    h2oveg_dead::Float64       = 0.0    # stored water in dead vegetation (kg/m2)
    h2oveg_growturn_err::Float64 = 0.0  # error water pool for growth/turnover tissue-volume change (kg/m2)
    h2oveg_hydro_err::Float64  = 0.0    # error water pool for hydrodynamics (kg/m2)

    # Useful diagnostics ----------------------------------------------------
    sapflow_scpf::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # flow at base of tree (+up) [kg/ha/s], size x pft

    rootuptake_sl::Vector{Float64} = Float64[]           # root uptake per SOIL layer [kg/m2/s]
    rootl_sl::Vector{Float64}      = Float64[]           # absorbing-root length on the soil grid

    rootuptake0_scpf::Matrix{Float64}   = Matrix{Float64}(undef, 0, 0)  # 0-10 cm   [kg/ha/m/s]
    rootuptake10_scpf::Matrix{Float64}  = Matrix{Float64}(undef, 0, 0)  # 10-50 cm
    rootuptake50_scpf::Matrix{Float64}  = Matrix{Float64}(undef, 0, 0)  # 50-100 cm
    rootuptake100_scpf::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # 100+ cm

    # Water-retention / conductance functions per soil (rhiz) layer ---------
    # Fortran: class(wrf_arr_type)/class(wkf_arr_type) pointer arrays. In Julia an
    # abstract-element vector holds the polymorphic functions directly.
    wrf_soil::Vector{WRFType} = WRFType[]                # water-retention function per soil layer
    wkf_soil::Vector{WKFType} = WKFType[]                # water-conductivity (K) function per soil layer

    # Matrix-solver system definition (connection / type map, soil -> leaf) -
    num_connections::Int = 0
    num_nodes::Int = 0
    conn_up::Vector{Int}    = Int[]
    conn_dn::Vector{Int}    = Int[]
    pm_node::Vector{Int}    = Int[]
    node_layer::Vector{Int} = Int[]
    ipiv::Vector{Int}       = Int[]                      # unused, returned from DGESV

    residual::Vector{Float64}       = Float64[]
    ajac::Matrix{Float64}           = Matrix{Float64}(undef, 0, 0)  # Jacobian (N terms, N equations)
    th_node_init::Vector{Float64}   = Float64[]
    th_node_prev::Vector{Float64}   = Float64[]
    th_node::Vector{Float64}        = Float64[]          # relative water content (theta) of node [m3/m3]
    dth_node::Vector{Float64}       = Float64[]          # time derivative of node water content
    h_node::Vector{Float64}         = Float64[]
    v_node::Vector{Float64}         = Float64[]          # volume of the node [m3]
    z_node::Vector{Float64}         = Float64[]          # elevation potential of the node (datum 0 = surface)
    psi_node::Vector{Float64}       = Float64[]          # suction of the node [MPa]
    q_flux::Vector{Float64}         = Float64[]          # mass flux of pathways between nodes
    dftc_dpsi_node::Vector{Float64} = Float64[]          # d(fraction total conductivity)/d(suction)
    ftc_node::Vector{Float64}       = Float64[]          # fraction of total conductivity [-]
    kmax_up::Vector{Float64}        = Float64[]          # max conductivity, upstream side of compartment
    kmax_dn::Vector{Float64}        = Float64[]          # max conductivity, downstream side of compartment

    # Scratch arrays (fixed nlevsoi_hyd_max length, matching Fortran) -------
    cohort_recruit_water_layer::Vector{Float64} = zeros(Float64, nlevsoi_hyd_max)  # recruit water requirement
    recruit_water_avail_layer::Vector{Float64}  = zeros(Float64, nlevsoi_hyd_max)  # recruit water available (kg H2O/m2)
end

# ===========================================================================
# Cohort-level plant-hydraulics state
# ===========================================================================

"""
    ed_cohort_hydr_type

Per-cohort plant-hydraulics state: node heights, maximum hydraulic conductances,
compartment volumes/lengths, the state variable (relative water content `th`),
and diagnostics (water potential `psi`, fraction of total conductivity `ftc`,
btran, transpiration, iteration/error tracking) by plant compartment.

Aboveground (leaf+stem) compartments and the transporting-root use fixed-size
vectors/scalars; absorbing-root and rhizosphere-layer quantities are allocatable
over `nlevrhiz`. Call [`AllocateHydrCohortArrays!`](@ref) to size the
per-layer arrays.
"""
Base.@kwdef mutable struct ed_cohort_hydr_type
    # Node heights of compartments [m], referenced to soil surface ----------
    z_node_ag::Vector{Float64}  = zeros(Float64, n_hypool_ag)  # nodal height of stem & leaf compartments (positive)
    z_upper_ag::Vector{Float64} = zeros(Float64, n_hypool_ag)  # upper-boundary heights (positive)
    z_lower_ag::Vector{Float64} = zeros(Float64, n_hypool_ag)  # lower-boundary heights (positive)
    z_node_troot::Float64       = 0.0                          # height of transporting-root node

    # Maximum hydraulic conductances [kg H2O s-1 MPa-1] ---------------------
    kmax_petiole_to_leaf::Float64       = 0.0                          # petiole -> leaf (nominally very high)
    kmax_stem_upper::Vector{Float64}    = zeros(Float64, n_hypool_stem)  # upper stem compartments
    kmax_stem_lower::Vector{Float64}    = zeros(Float64, n_hypool_stem)  # lower stem compartments
    kmax_troot_upper::Float64           = 0.0                          # upper transporting root
    kmax_troot_lower::Vector{Float64}      = Float64[]  # troot portion joining each aroot compartment
    kmax_aroot_upper::Vector{Float64}      = Float64[]  # aroot xylem into the transporting root
    kmax_aroot_lower::Vector{Float64}      = Float64[]  # aroot edge -> node center (hybrid troot volume)
    kmax_aroot_radial_in::Vector{Float64}  = Float64[]  # radial, potential gradient into the root
    kmax_aroot_radial_out::Vector{Float64} = Float64[]  # radial, potential gradient out of the root

    # Compartment volumes and lengths ---------------------------------------
    v_ag_init::Vector{Float64} = zeros(Float64, n_hypool_ag)  # previous day's aboveground storage volume [m3]
    v_ag::Vector{Float64}      = zeros(Float64, n_hypool_ag)  # aboveground storage volume [m3]
    v_troot_init::Float64      = 0.0                          # previous day's belowground storage volume [m3]
    v_troot::Float64           = 0.0                          # belowground storage volume [m3]
    v_aroot_layer_init::Vector{Float64} = Float64[]  # previous day's absorbing-root volume by soil layer [m3]
    v_aroot_layer::Vector{Float64}      = Float64[]  # absorbing-root volume by soil layer [m3]
    l_aroot_layer::Vector{Float64}      = Float64[]  # absorbing-root length by soil layer [m]

    # State variable: relative water content by volume ("theta") ------------
    th_ag::Vector{Float64}   = zeros(Float64, n_hypool_ag)  # water in aboveground compartments [kgh2o/indiv]
    th_troot::Float64        = 0.0                          # water in belowground compartment [kgh2o/indiv]
    th_aroot::Vector{Float64} = Float64[]                   # water in absorbing roots [kgh2o/indiv]

    # Diagnostic: water potential [MPa] -------------------------------------
    psi_ag::Vector{Float64}   = zeros(Float64, n_hypool_ag)  # aboveground compartments
    psi_troot::Float64        = 0.0                          # transporting root
    psi_aroot::Vector{Float64} = Float64[]                   # absorbing roots

    # Diagnostic: fraction of total conductivity [-] ------------------------
    ftc_ag::Vector{Float64}   = zeros(Float64, n_hypool_ag)  # aboveground compartments
    ftc_troot::Float64        = 0.0                          # transporting root
    ftc_aroot::Vector{Float64} = Float64[]                   # absorbing roots

    btran::Float64 = 0.0   # leaf-water-potential limitation on gs [0-1]
    qtop::Float64  = 0.0   # mean transpiration flux rate [kg/cohort/s]

    # Error tracking / flagging ---------------------------------------------
    supsub_flag::Float64 = 0.0  # k index of last node to encounter super/sub saturation (+sup / -sub)
    iterh1::Float64      = 0.0   # max iterations to achieve tolerable water-balance error
    iterh2::Float64      = 0.0   # number of inner iterations
    iterlayer::Float64   = 0.0   # layer index with the highest iterations
    errh2o::Float64      = 0.0   # total water-balance error per unit crown area [kgh2o/m2]

    # Other -----------------------------------------------------------------
    is_newly_recruited::Bool = false  # whether the new cohort is newly recruited
end

# ===========================================================================
# Cohort-level routines
# ===========================================================================

"""
    AllocateHydrCohortArrays!(this::ed_cohort_hydr_type, nlevrhiz::Integer)

Allocate (size) the per-rhizosphere-layer arrays of a cohort hydraulics object,
zero-initialized. Mirrors the Fortran `AllocateHydrCohortArrays`.
"""
function AllocateHydrCohortArrays!(this::ed_cohort_hydr_type, nlevrhiz::Integer)
    this.kmax_troot_lower      = zeros(Float64, nlevrhiz)
    this.kmax_aroot_upper      = zeros(Float64, nlevrhiz)
    this.kmax_aroot_lower      = zeros(Float64, nlevrhiz)
    this.kmax_aroot_radial_in  = zeros(Float64, nlevrhiz)
    this.kmax_aroot_radial_out = zeros(Float64, nlevrhiz)
    this.v_aroot_layer_init    = zeros(Float64, nlevrhiz)
    this.v_aroot_layer         = zeros(Float64, nlevrhiz)
    this.l_aroot_layer         = zeros(Float64, nlevrhiz)
    this.th_aroot              = zeros(Float64, nlevrhiz)
    this.psi_aroot             = zeros(Float64, nlevrhiz)
    this.ftc_aroot             = zeros(Float64, nlevrhiz)
    return nothing
end

"""
    DeallocateHydrCohortArrays!(this::ed_cohort_hydr_type)

Release the per-rhizosphere-layer arrays of a cohort hydraulics object (resets
them to undef-0). Mirrors the Fortran `DeallocateHydrCohortArrays`.
"""
function DeallocateHydrCohortArrays!(this::ed_cohort_hydr_type)
    this.kmax_troot_lower      = Float64[]
    this.kmax_aroot_upper      = Float64[]
    this.kmax_aroot_lower      = Float64[]
    this.kmax_aroot_radial_in  = Float64[]
    this.kmax_aroot_radial_out = Float64[]
    this.v_aroot_layer_init    = Float64[]
    this.v_aroot_layer         = Float64[]
    this.l_aroot_layer         = Float64[]
    this.th_aroot              = Float64[]
    this.psi_aroot             = Float64[]
    this.ftc_aroot             = Float64[]
    return nothing
end

"""
    CopyCohortHydraulics!(ncohort_hydr, ocohort_hydr)

Copy the hydraulics state from `ocohort_hydr` into `ncohort_hydr`. Per-layer
arrays are copied element-wise (`.=` semantics via `copy`). Mirrors the Fortran
`CopyCohortHydraulics`.
"""
function CopyCohortHydraulics!(ncohort_hydr::ed_cohort_hydr_type,
                               ocohort_hydr::ed_cohort_hydr_type)
    # Node heights
    ncohort_hydr.z_node_ag    = copy(ocohort_hydr.z_node_ag)
    ncohort_hydr.z_upper_ag   = copy(ocohort_hydr.z_upper_ag)
    ncohort_hydr.z_lower_ag   = copy(ocohort_hydr.z_lower_ag)
    ncohort_hydr.z_node_troot = ocohort_hydr.z_node_troot

    # Compartment kmax's
    ncohort_hydr.kmax_petiole_to_leaf  = ocohort_hydr.kmax_petiole_to_leaf
    ncohort_hydr.kmax_stem_lower       = copy(ocohort_hydr.kmax_stem_lower)
    ncohort_hydr.kmax_stem_upper       = copy(ocohort_hydr.kmax_stem_upper)
    ncohort_hydr.kmax_troot_upper      = ocohort_hydr.kmax_troot_upper
    ncohort_hydr.kmax_troot_lower      = copy(ocohort_hydr.kmax_troot_lower)
    ncohort_hydr.kmax_aroot_upper      = copy(ocohort_hydr.kmax_aroot_upper)
    ncohort_hydr.kmax_aroot_lower      = copy(ocohort_hydr.kmax_aroot_lower)
    ncohort_hydr.kmax_aroot_radial_in  = copy(ocohort_hydr.kmax_aroot_radial_in)
    ncohort_hydr.kmax_aroot_radial_out = copy(ocohort_hydr.kmax_aroot_radial_out)

    # Compartment volumes
    ncohort_hydr.v_ag_init          = copy(ocohort_hydr.v_ag_init)
    ncohort_hydr.v_ag               = copy(ocohort_hydr.v_ag)
    ncohort_hydr.v_troot_init       = ocohort_hydr.v_troot_init
    ncohort_hydr.v_troot            = ocohort_hydr.v_troot
    ncohort_hydr.v_aroot_layer_init = copy(ocohort_hydr.v_aroot_layer_init)
    ncohort_hydr.v_aroot_layer      = copy(ocohort_hydr.v_aroot_layer)
    ncohort_hydr.l_aroot_layer      = copy(ocohort_hydr.l_aroot_layer)

    # State variables
    ncohort_hydr.th_ag     = copy(ocohort_hydr.th_ag)
    ncohort_hydr.th_troot  = ocohort_hydr.th_troot
    ncohort_hydr.th_aroot  = copy(ocohort_hydr.th_aroot)
    ncohort_hydr.psi_ag    = copy(ocohort_hydr.psi_ag)
    ncohort_hydr.psi_troot = ocohort_hydr.psi_troot
    ncohort_hydr.psi_aroot = copy(ocohort_hydr.psi_aroot)
    ncohort_hydr.ftc_ag    = copy(ocohort_hydr.ftc_ag)
    ncohort_hydr.ftc_troot = ocohort_hydr.ftc_troot
    ncohort_hydr.ftc_aroot = copy(ocohort_hydr.ftc_aroot)

    # Other
    ncohort_hydr.btran       = ocohort_hydr.btran
    ncohort_hydr.supsub_flag = ocohort_hydr.supsub_flag
    ncohort_hydr.iterh1      = ocohort_hydr.iterh1
    ncohort_hydr.iterh2      = ocohort_hydr.iterh2
    ncohort_hydr.iterlayer   = ocohort_hydr.iterlayer
    ncohort_hydr.errh2o      = ocohort_hydr.errh2o

    # BC plant hydraulics - flux terms
    ncohort_hydr.qtop = ocohort_hydr.qtop

    ncohort_hydr.is_newly_recruited = ocohort_hydr.is_newly_recruited
    return nothing
end

# ===========================================================================
# Site-level routines
# ===========================================================================

"""
    InitHydrSite!(this::ed_site_hydr_type, numpft, numlevsclass, hydr_solver_type, nlevsoil)

Allocate and initialize the site-level hydraulics arrays. `this.nlevrhiz` must
already be set. Allocatable rhizosphere arrays are sized from `nlevrhiz`/`nshell`
and initialized as in Fortran (`NaN`, `-999`, `fine_root_radius_const`, or 0).
The matrix-solver system (`conn_*`, `pm_node`, node arrays) is sized per the
solver type, then [`SetConnections!`](@ref) populates the static topology.

Mirrors the Fortran `InitHydrSite`.
"""
function InitHydrSite!(this::ed_site_hydr_type, numpft::Integer, numlevsclass::Integer,
                       hydr_solver_type::Integer, nlevsoil::Integer)
    nlevrhiz = this.nlevrhiz

    # In all cases, the 0 index of the layer bottom is a value of 0.
    this.zi_rhiz = fill(NaN, nlevrhiz)
    this.dz_rhiz = fill(NaN, nlevrhiz)
    this.map_s2r = fill(-999, nlevrhiz)
    this.map_r2s = fill(-999, nlevrhiz, 2)
    this.v_shell           = fill(NaN, nlevrhiz, nshell)
    this.v_shell_init      = fill(NaN, nlevrhiz, nshell)
    this.r_node_shell      = fill(NaN, nlevrhiz, nshell)
    this.r_node_shell_init = fill(NaN, nlevrhiz, nshell)
    this.r_out_shell       = fill(NaN, nlevrhiz, nshell)
    this.l_aroot_layer      = fill(NaN, nlevrhiz)
    this.l_aroot_layer_init = fill(NaN, nlevrhiz)
    this.kmax_upper_shell = fill(NaN, nlevrhiz, nshell)
    this.kmax_lower_shell = fill(NaN, nlevrhiz, nshell)
    this.supsub_flag = fill(-999, nlevrhiz)
    this.h2osoi_liqvol_shell = fill(NaN, nlevrhiz, nshell)
    this.rs1 = fill(fine_root_radius_const, nlevrhiz)
    this.recruit_w_uptake = fill(NaN, nlevrhiz)

    this.rootuptake_sl = fill(NaN, nlevsoil)
    this.rootl_sl      = zeros(Float64, nlevsoil)

    this.sapflow_scpf       = fill(NaN, numlevsclass, numpft)
    this.rootuptake0_scpf   = fill(NaN, numlevsclass, numpft)
    this.rootuptake10_scpf  = fill(NaN, numlevsclass, numpft)
    this.rootuptake50_scpf  = fill(NaN, numlevsclass, numpft)
    this.rootuptake100_scpf = fill(NaN, numlevsclass, numpft)

    this.errh2o_hyd     = NaN
    this.dwat_veg       = NaN
    this.h2oveg         = 0.0
    this.h2oveg_recruit = 0.0
    this.h2oveg_dead    = 0.0

    this.h2oveg_growturn_err = 0.0
    this.h2oveg_hydro_err    = 0.0

    # Separate water-transfer functions and parameters for each soil layer and
    # plant compartment type. Allocate the per-layer holder vectors (entries are
    # filled by the hydraulics-driver setup, downstream of this memory module).
    this.wrf_soil = Vector{WRFType}(undef, nlevrhiz)
    this.wkf_soil = Vector{WKFType}(undef, nlevrhiz)

    if hydr_solver_type == hydr_solver_2DNewton || hydr_solver_type == hydr_solver_2DPicard
        this.num_connections = n_hypool_leaf + n_hypool_stem + n_hypool_troot - 1 +
                               (n_hypool_aroot + nshell) * nlevrhiz
        this.num_nodes = n_hypool_leaf + n_hypool_stem + n_hypool_troot +
                         (n_hypool_aroot + nshell) * nlevrhiz

        # These are only in the Newton-matrix solve.
        this.conn_up        = zeros(Int, this.num_connections)
        this.conn_dn        = zeros(Int, this.num_connections)
        this.residual       = zeros(Float64, this.num_nodes)
        this.ajac           = zeros(Float64, this.num_nodes, this.num_nodes)
        this.th_node_init   = zeros(Float64, this.num_nodes)
        this.th_node_prev   = zeros(Float64, this.num_nodes)
        this.th_node        = zeros(Float64, this.num_nodes)
        this.dth_node       = zeros(Float64, this.num_nodes)
        this.h_node         = zeros(Float64, this.num_nodes)
        this.v_node         = zeros(Float64, this.num_nodes)
        this.z_node         = zeros(Float64, this.num_nodes)
        this.psi_node       = zeros(Float64, this.num_nodes)
        this.q_flux         = zeros(Float64, this.num_connections)
        this.dftc_dpsi_node = zeros(Float64, this.num_nodes)
        this.ftc_node       = zeros(Float64, this.num_nodes)
        this.pm_node        = zeros(Int, this.num_nodes)
        this.ipiv           = zeros(Int, this.num_nodes)
        this.node_layer     = zeros(Int, this.num_nodes)

        this.kmax_up = zeros(Float64, this.num_connections)
        this.kmax_dn = zeros(Float64, this.num_connections)
    else
        this.num_connections = n_hypool_leaf + n_hypool_stem +
                               n_hypool_troot + n_hypool_aroot + nshell - 1
        this.num_nodes = n_hypool_leaf + n_hypool_stem +
                         n_hypool_troot + n_hypool_aroot + nshell

        this.conn_up = zeros(Int, this.num_connections)
        this.conn_dn = zeros(Int, this.num_connections)
        this.pm_node = zeros(Int, this.num_nodes)
    end

    SetConnections!(this, hydr_solver_type)
    return nothing
end

"""
    FlushSiteScratch!(this::ed_site_hydr_type, hydr_solver_type::Integer)

Reset the matrix-solver scratch arrays to `fates_unset_r8`. Only acts for the
2D (Newton/Picard) solvers, matching the Fortran `FlushSiteScratch`.
"""
function FlushSiteScratch!(this::ed_site_hydr_type, hydr_solver_type::Integer)
    if hydr_solver_type == hydr_solver_2DNewton || hydr_solver_type == hydr_solver_2DPicard
        fill!(this.residual, fates_unset_r8)
        fill!(this.ajac, fates_unset_r8)
        fill!(this.th_node_init, fates_unset_r8)
        fill!(this.th_node_prev, fates_unset_r8)
        fill!(this.th_node, fates_unset_r8)
        fill!(this.dth_node, fates_unset_r8)
        fill!(this.h_node, fates_unset_r8)
        fill!(this.v_node, fates_unset_r8)
        fill!(this.z_node, fates_unset_r8)
        fill!(this.psi_node, fates_unset_r8)
        fill!(this.ftc_node, fates_unset_r8)
        fill!(this.dftc_dpsi_node, fates_unset_r8)
        # kmax_up / kmax_dn intentionally NOT flushed (commented out in Fortran).
        fill!(this.q_flux, fates_unset_r8)
    end
    return nothing
end

"""
    AggBCToRhiz(this::ed_site_hydr_type, var_in, j, weight) -> var_out

Aggregate a soil-layer property to the rhizosphere (root) layer `j` using a
harmonic mean (the Fortran default `mean_type = harmonic_mean`) over the soil
layers `map_r2s[j,1]:map_r2s[j,2]`. Mirrors the Fortran `AggBCToRhiz`.
"""
function AggBCToRhiz(this::ed_site_hydr_type, var_in::AbstractVector{<:Real},
                     j::Integer, weight::AbstractVector{<:Real})
    j_t = this.map_r2s[j, 1]
    j_b = this.map_r2s[j, 2]

    # Harmonic mean (mean_type == harmonic_mean).
    w = @view weight[j_t:j_b]
    v = @view var_in[j_t:j_b]
    return sum(w) / sum(w ./ v)
end

"""
    SetConnections!(this::ed_site_hydr_type, hydr_solver_type::Integer)

Populate the static node/connection topology (`conn_up`, `conn_dn`, `pm_node`,
and, for the 2D solvers, `node_layer`) describing the soil->leaf hydraulic
network. Mirrors the Fortran `SetConnections`.
"""
function SetConnections!(this::ed_site_hydr_type, hydr_solver_type::Integer)
    num_cnxs = 0
    num_nds  = 0
    for k in 1:n_hypool_leaf
        num_cnxs += 1
        num_nds  += 1
        this.conn_dn[num_cnxs] = k       # leaf is the dn, origin, bottom
        this.conn_up[num_cnxs] = k + 1
        this.pm_node[num_nds]  = leaf_p_media
    end
    for k in (n_hypool_leaf + 1):n_hypool_ag
        num_cnxs += 1
        num_nds  += 1
        this.conn_dn[num_cnxs] = k
        this.conn_up[num_cnxs] = k + 1
        this.pm_node[num_nds]  = stem_p_media
    end

    if hydr_solver_type == hydr_solver_2DNewton || hydr_solver_type == hydr_solver_2DPicard
        num_nds     = n_hypool_ag + n_hypool_troot
        node_tr_end = num_nds
        # nt_ab kept for traceability to Fortran (computed but unused there too).
        nt_ab       = n_hypool_ag + n_hypool_troot + n_hypool_aroot
        num_cnxs    = n_hypool_ag

        this.pm_node[num_nds] = troot_p_media
        this.node_layer[1:n_hypool_ag] .= 0
        this.node_layer[num_nds] = 1

        for j in 1:this.nlevrhiz
            for k in 1:(n_hypool_aroot + nshell)
                num_nds  += 1
                num_cnxs += 1
                this.node_layer[num_nds] = j
                if k == 1  # troot-aroot junction node
                    this.conn_dn[num_cnxs] = node_tr_end  # absorbing root
                    this.conn_up[num_cnxs] = num_nds
                    this.pm_node[num_nds]  = aroot_p_media
                else
                    this.conn_dn[num_cnxs] = num_nds - 1
                    this.conn_up[num_cnxs] = num_nds
                    this.pm_node[num_nds]  = rhiz_p_media
                end
            end
        end
    else
        this.pm_node[n_hypool_ag + 1] = troot_p_media
        this.pm_node[n_hypool_ag + 2] = aroot_p_media
        this.pm_node[(n_hypool_ag + 3):(n_hypool_ag + 2 + nshell)] .= rhiz_p_media
    end
    return nothing
end
