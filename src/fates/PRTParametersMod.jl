# PRTParametersMod.jl
# Julia port of FATES src/fates/parteh/PRTParametersMod.F90
#
# This module only holds the parameter definitions for PARTEH (Plant Allocation
# & Reactive Transport) and allometry. It does NOT hold any of the code used for
# initializing and filling that data — that is model dependent (i.e. FATES may
# have a different way than another TBM).
#
# Fortran `prt_param_type` -> Julia `prt_param_type` (mutable struct, SoA).
# `real(r8), allocatable :: x(:)`   -> Vector{Float64}
# `real(r8), allocatable :: x(:,:)` -> Matrix{Float64}
# `integer , allocatable :: x(:)`   -> Vector{Int}
# Allocatables start unallocated -> zero-length arrays; `allocate_prt_params!`
# sizes them, mirroring the Fortran allocate that the (model-dependent) reader
# would perform. Arrays are dimensioned (PFT), (PFT x organ), (PFT x ageclass),
# or (organ)/(all-organs) as documented per field.
# `fates_r8` -> Float64. Deps: FatesConstantsMod (fates_r8).

"""
    prt_param_type

Storage for the PARTEH and allometry parameters, indexed by PFT, organ,
element/age-class as documented per field. SoA layout (one array per
parameter). Construct empty (all zero-length arrays) then call
[`allocate_prt_params!`](@ref) with the dimension sizes to size every array.
"""
Base.@kwdef mutable struct prt_param_type
    # The following three PFT classes are mutually exclusive (dim: PFT).
    stress_decid::Vector{Int} = Int[]   # Is the plant stress deciduous?
                                        #   0 - No
                                        #   1 - Drought "hard" deciduous
                                        #   2 - Drought semi-deciduous
    season_decid::Vector{Int} = Int[]   # Is the plant seasonally deciduous (1=yes, 0=no)
    evergreen::Vector{Int}    = Int[]   # Is the plant an evergreen (1=yes, 0=no)

    # Drop fraction for tissues other than leaves (dim: PFT).
    phen_fnrt_drop_fraction::Vector{Float64} = Float64[]  # Abscission fraction of fine roots
    phen_stem_drop_fraction::Vector{Float64} = Float64[]  # Abscission fraction of stems
    phen_drought_threshold::Vector{Float64} = Float64[]   # Hard drought-deciduous abscission/flush threshold
    phen_moist_threshold::Vector{Float64}   = Float64[]   # Semi-deciduous complete-flush threshold
    phen_doff_time::Vector{Float64}         = Float64[]   # Min days leafless before flushing again

    # Growth and Turnover Parameters
    senleaf_long_fdrought::Vector{Float64} = Float64[]    # Leaf-longevity multiplier for senescent leaves (dim: PFT)
    leaf_long::Matrix{Float64}        = Matrix{Float64}(undef, 0, 0)  # Leaf turnover time (PFT x age-class) [yr]
    leaf_long_ustory::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # As above but for understory trees
    root_long::Vector{Float64}   = Float64[]   # Root turnover time (longevity) (dim: PFT) [yr]
    branch_long::Vector{Float64} = Float64[]   # Branchfall turnover on live trees (dim: PFT) [yr]
    turnover_nitr_retrans::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # N re-translocation fraction (PFT x organ)
    turnover_phos_retrans::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # P re-translocation fraction (PFT x organ)

    # Vertical N profile scalers (affect crown allometry and sla) (dim: PFT).
    leafn_vert_scaler_coeff1::Vector{Float64} = Float64[]  # Coefficient one for leaf-N decrease through canopy
    leafn_vert_scaler_coeff2::Vector{Float64} = Float64[]  # Coefficient two for leaf-N decrease through canopy

    grperc::Vector{Float64} = Float64[]   # Growth respiration per unit C gained (dim: PFT) [kg/kg]

    nitr_stoich_p1::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # N stoichiometry param 1 (PFT x organ)
    phos_stoich_p1::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # P stoichiometry param 1 (PFT x organ)

    nitr_store_ratio::Vector{Float64} = Float64[]  # Target stored-N per tissue-bound-N ratio (dim: PFT)
    phos_store_ratio::Vector{Float64} = Float64[]  # Target stored-P per tissue-bound-P ratio (dim: PFT)

    organ_id::Vector{Int} = Int[]   # Map of organ index in param file -> global PRTGenericMod organ list (dim: organ)
    alloc_priority::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # Allocation priority per organ (PFT x organ) [0-6]
    cushion::Vector{Float64} = Float64[]   # Labile C storage target as multiple of leaf pool (dim: PFT)
    leaf_stor_priority::Vector{Float64} = Float64[]  # Leaf turnover vs labile-C use prioritisation (dim: PFT)
    dbh_repro_threshold::Vector{Float64} = Float64[] # Diameter at which mature plants shift allocation (dim: PFT)
    seed_alloc_mature::Vector{Float64} = Float64[]   # Fraction of C balance to clonal reproduction (dim: PFT)
    seed_alloc::Vector{Float64} = Float64[]          # Fraction of C balance allocated to seeds (dim: PFT)
    repro_alloc_a::Vector{Float64} = Float64[]       # Sigmoidal shape param: dbh -> seed alloc fraction (dim: PFT)
    repro_alloc_b::Vector{Float64} = Float64[]       # Intercept param: dbh -> seed alloc fraction (dim: PFT)

    # Derived parameters
    organ_param_id::Vector{Int} = Int[]  # Sparse reverse lookup: dim by all PRT organs, -> param-file index or -1

    # Allometry Parameters --------------------------------------------------------
    # Root profile parameters (dim: PFT).
    fnrt_prof_mode::Vector{Float64} = Float64[]  # Fine root profile functional form
    fnrt_prof_a::Vector{Float64}    = Float64[]  # Fine root profile scaling parameter A
    fnrt_prof_b::Vector{Float64}    = Float64[]  # Fine root profile scaling parameter B

    c2b::Vector{Float64}          = Float64[]  # Carbon to biomass multiplier [kg/kgC] (dim: PFT)
    wood_density::Vector{Float64} = Float64[]  # Wood density [g cm^-3] (dim: PFT)
    woody::Vector{Int}            = Int[]      # Does the plant have wood? (1=yes, 0=no) (dim: PFT)
    slamax::Vector{Float64}       = Float64[]  # Maximum specific leaf area (at bottom) [m2/gC] (dim: PFT)
    slatop::Vector{Float64}       = Float64[]  # Specific leaf area at canopy top [m2/gC] (dim: PFT)
    allom_sai_scaler::Vector{Float64} = Float64[]      # (dim: PFT)
    allom_dbh_maxheight::Vector{Float64} = Float64[]   # dbh at which height growth ceases (dim: PFT)
    allom_hmode::Vector{Int}  = Int[]   # Height allometry function type (dim: PFT)
    allom_lmode::Vector{Int}  = Int[]   # Maximum leaf allometry function type (dim: PFT)
    allom_fmode::Vector{Int}  = Int[]   # Maximum root allometry function type (dim: PFT)
    allom_amode::Vector{Int}  = Int[]   # AGB allometry function type (dim: PFT)
    allom_cmode::Vector{Int}  = Int[]   # Coarse root allometry function type (dim: PFT)
    allom_smode::Vector{Int}  = Int[]   # Sapwood allometry function type (dim: PFT)
    allom_stmode::Vector{Int} = Int[]   # Storage allometry functional type (dim: PFT)
    allom_dmode::Vector{Int}  = Int[]   # Crown depth allometry function type (dim: PFT)
    allom_la_per_sa_int::Vector{Float64} = Float64[]  # Leaf area to sap area conversion, intercept [cm2/m2] (dim: PFT)
    allom_la_per_sa_slp::Vector{Float64} = Float64[]  # Leaf area to sap area conversion, slope [cm2/m2/cm] (dim: PFT)
    allom_l2fr::Vector{Float64}    = Float64[]  # Fine root biomass per leaf biomass ratio [kgC/kgC] (dim: PFT)
    allom_agb_frac::Vector{Float64} = Float64[] # Fraction of stem above ground [-] (dim: PFT)
    allom_d2h1::Vector{Float64} = Float64[]  # Param 1 for d2h allometry (intercept, "c") (dim: PFT)
    allom_d2h2::Vector{Float64} = Float64[]  # Param 2 for d2h allometry (slope, "m") (dim: PFT)
    allom_d2h3::Vector{Float64} = Float64[]  # Param 3 for d2h allometry (optional) (dim: PFT)
    allom_d2bl1::Vector{Float64} = Float64[] # Param 1 for d2bl allometry (intercept) (dim: PFT)
    allom_d2bl2::Vector{Float64} = Float64[] # Param 2 for d2bl allometry (slope) (dim: PFT)
    allom_d2bl3::Vector{Float64} = Float64[] # Param 3 for d2bl allometry (optional) (dim: PFT)
    allom_blca_expnt_diff::Vector{Float64} = Float64[]      # Leaf-biomass vs crown-area scaling exponent diff (dim: PFT)
    allom_d2ca_coefficient_max::Vector{Float64} = Float64[] # Upper (savanna) crown-area-to-dbh coeff (dim: PFT)
    allom_d2ca_coefficient_min::Vector{Float64} = Float64[] # Lower (closed-canopy) crown-area-to-dbh coeff (dim: PFT)
    allom_agb1::Vector{Float64} = Float64[]  # Param 1 for agb allometry (dim: PFT)
    allom_agb2::Vector{Float64} = Float64[]  # Param 2 for agb allometry (dim: PFT)
    allom_agb3::Vector{Float64} = Float64[]  # Param 3 for agb allometry (dim: PFT)
    allom_agb4::Vector{Float64} = Float64[]  # Param 4 for agb allometry (dim: PFT)
    allom_h2cd1::Vector{Float64} = Float64[] # Param 1 for crown depth allometry (dim: PFT)
    allom_h2cd2::Vector{Float64} = Float64[] # Exponent for crown depth allometry (dim: PFT)
    allom_zroot_max_dbh::Vector{Float64} = Float64[]  # dbh at which max rooting depth saturates [cm] (dim: PFT)
    allom_zroot_max_z::Vector{Float64}   = Float64[]  # Max rooting depth at allom_zroot_max_dbh [m] (dim: PFT)
    allom_zroot_min_dbh::Vector{Float64} = Float64[]  # dbh at which max recruit rooting depth defined [cm] (dim: PFT)
    allom_zroot_min_z::Vector{Float64}   = Float64[]  # Max rooting depth at allom_zroot_min_dbh [m] (dim: PFT)
    allom_zroot_k::Vector{Float64}       = Float64[]  # Scale coefficient of logistic rooting depth model (dim: PFT)

    # PID controller parameters (dim: PFT).
    pid_kp::Vector{Float64} = Float64[]  # Proportion constant in fine-root biomass PID controller
    pid_ki::Vector{Float64} = Float64[]  # Integral constant in fine-root biomass PID controller
    pid_kd::Vector{Float64} = Float64[]  # Derivative constant in fine-root biomass PID controller

    store_ovrflw_frac::Vector{Float64} = Float64[]  # Excess storage retention fraction beyond target (dim: PFT)

    nfix_mresp_scfrac::Vector{Float64} = Float64[]  # Maint-resp surcharge fraction to pay for N-fixation (dim: PFT)
end

"""
    allocate_prt_params!(p::prt_param_type, npft::Integer, norgan::Integer,
                         nleafage::Integer; nall_organs::Integer = norgan)

Allocate (size) all parameter arrays in `p`. Mirrors the Fortran allocate that
the (model-dependent) parameter reader would perform. Dimensions:
- PFT-indexed arrays  -> length `npft`
- (PFT x organ)       -> `npft × norgan`
- (PFT x age-class)   -> `npft × nleafage` (leaf longevity)
- organ-indexed       -> `organ_id` length `norgan`
- all-PRT-organs      -> `organ_param_id` length `nall_organs`

All elements are zero-initialised (`0` for integers, `0.0` for reals).
"""
function allocate_prt_params!(p::prt_param_type, npft::Integer, norgan::Integer,
                              nleafage::Integer; nall_organs::Integer = norgan)
    npft     >= 0 || error("allocate_prt_params!: npft must be >= 0")
    norgan   >= 0 || error("allocate_prt_params!: norgan must be >= 0")
    nleafage >= 0 || error("allocate_prt_params!: nleafage must be >= 0")

    zi = () -> zeros(Int, npft)
    zr = () -> zeros(Float64, npft)

    # Phenology / class flags (PFT)
    p.stress_decid = zi();  p.season_decid = zi();  p.evergreen = zi()
    p.phen_fnrt_drop_fraction = zr()
    p.phen_stem_drop_fraction = zr()
    p.phen_drought_threshold  = zr()
    p.phen_moist_threshold    = zr()
    p.phen_doff_time          = zr()

    # Growth and turnover
    p.senleaf_long_fdrought = zr()
    p.leaf_long        = zeros(Float64, npft, nleafage)
    p.leaf_long_ustory = zeros(Float64, npft, nleafage)
    p.root_long        = zr()
    p.branch_long      = zr()
    p.turnover_nitr_retrans = zeros(Float64, npft, norgan)
    p.turnover_phos_retrans = zeros(Float64, npft, norgan)

    p.leafn_vert_scaler_coeff1 = zr()
    p.leafn_vert_scaler_coeff2 = zr()

    p.grperc = zr()

    p.nitr_stoich_p1 = zeros(Float64, npft, norgan)
    p.phos_stoich_p1 = zeros(Float64, npft, norgan)

    p.nitr_store_ratio = zr()
    p.phos_store_ratio = zr()

    p.organ_id       = zeros(Int, norgan)
    p.alloc_priority = zeros(Float64, npft, norgan)
    p.cushion        = zr()
    p.leaf_stor_priority = zr()
    p.dbh_repro_threshold = zr()
    p.seed_alloc_mature = zr()
    p.seed_alloc        = zr()
    p.repro_alloc_a     = zr()
    p.repro_alloc_b     = zr()

    p.organ_param_id = zeros(Int, nall_organs)

    # Allometry
    p.fnrt_prof_mode = zr();  p.fnrt_prof_a = zr();  p.fnrt_prof_b = zr()
    p.c2b = zr();  p.wood_density = zr();  p.woody = zi()
    p.slamax = zr();  p.slatop = zr()
    p.allom_sai_scaler = zr()
    p.allom_dbh_maxheight = zr()
    p.allom_hmode = zi();  p.allom_lmode = zi();  p.allom_fmode = zi()
    p.allom_amode = zi();  p.allom_cmode = zi();  p.allom_smode = zi()
    p.allom_stmode = zi();  p.allom_dmode = zi()
    p.allom_la_per_sa_int = zr();  p.allom_la_per_sa_slp = zr()
    p.allom_l2fr = zr();  p.allom_agb_frac = zr()
    p.allom_d2h1 = zr();  p.allom_d2h2 = zr();  p.allom_d2h3 = zr()
    p.allom_d2bl1 = zr();  p.allom_d2bl2 = zr();  p.allom_d2bl3 = zr()
    p.allom_blca_expnt_diff = zr()
    p.allom_d2ca_coefficient_max = zr()
    p.allom_d2ca_coefficient_min = zr()
    p.allom_agb1 = zr();  p.allom_agb2 = zr();  p.allom_agb3 = zr();  p.allom_agb4 = zr()
    p.allom_h2cd1 = zr();  p.allom_h2cd2 = zr()
    p.allom_zroot_max_dbh = zr();  p.allom_zroot_max_z = zr()
    p.allom_zroot_min_dbh = zr();  p.allom_zroot_min_z = zr()
    p.allom_zroot_k = zr()

    # PID controller
    p.pid_kp = zr();  p.pid_ki = zr();  p.pid_kd = zr()

    p.store_ovrflw_frac = zr()
    p.nfix_mresp_scfrac = zr()

    return p
end

# Instantiation of the parameter object (Fortran: `prt_params`). A module-level
# singleton mirroring the Fortran `type(prt_param_type),public :: prt_params`.
const prt_params = prt_param_type()
