# ==========================================================================
# Ported from: src/main/pftconMod.F90
# PFT (Plant Functional Type) constants and data structure
#
# NOTE: Fortran uses 0-based PFT indexing (allocate(0:mxpft)).
# Julia uses 1-based arrays of size (MXPFT+1), where Julia index i
# corresponds to Fortran index (i-1).
# ==========================================================================

# --- Module-level PFT index constants ---
# These are mutable so they can be set during initialization (like Fortran).
# Stored as Fortran 0-based values for traceability. Add +1 when indexing Julia arrays.
noveg                 = 0   # not vegetated
ndllf_evr_tmp_tree    = 1   # Needleleaf evergreen temperate tree
ndllf_evr_brl_tree    = 2   # Needleleaf evergreen boreal tree
ndllf_dcd_brl_tree    = 3   # Needleleaf deciduous boreal tree
nbrdlf_evr_trp_tree   = 4   # Broadleaf evergreen tropical tree
nbrdlf_evr_tmp_tree   = 5   # Broadleaf evergreen temperate tree
nbrdlf_dcd_trp_tree   = 6   # Broadleaf deciduous tropical tree
nbrdlf_dcd_tmp_tree   = 7   # Broadleaf deciduous temperate tree
nbrdlf_dcd_brl_tree   = 8   # Broadleaf deciduous boreal tree
nbrdlf_evr_shrub      = 9   # Broadleaf evergreen shrub
nbrdlf_dcd_tmp_shrub  = 10  # Broadleaf deciduous temperate shrub
nbrdlf_dcd_brl_shrub  = 11  # Broadleaf deciduous boreal shrub
nc3_arctic_grass      = 12  # C3 arctic grass
nc3_nonarctic_grass   = 13  # C3 non-arctic grass
nc4_grass             = 14  # C4 grass
npcropmin             = 0   # first prognostic crop (set during init)
ntmp_corn             = 17  # temperate corn (rf)
nirrig_tmp_corn       = 18  # temperate corn (ir)
nswheat               = 19  # spring wheat (rf)
nirrig_swheat         = 20  # spring wheat (ir)
nwwheat               = 21  # winter wheat (rf)
nirrig_wwheat         = 22  # winter wheat (ir)
ntmp_soybean          = 23  # temperate soybean (rf)
nirrig_tmp_soybean    = 24  # temperate soybean (ir)
nbarley               = 25  # spring barley (rf)
nirrig_barley         = 26  # spring barley (ir)
nwbarley              = 27  # winter barley (rf)
nirrig_wbarley        = 28  # winter barley (ir)
nrye                  = 29  # spring rye (rf)
nirrig_rye            = 30  # spring rye (ir)
nwrye                 = 31  # winter rye (rf)
nirrig_wrye           = 32  # winter rye (ir)
ncassava              = 33
nirrig_cassava        = 34
ncitrus               = 35
nirrig_citrus         = 36
ncocoa                = 37
nirrig_cocoa          = 38
ncoffee               = 39
nirrig_coffee         = 40
ncotton               = 41
nirrig_cotton         = 42
ndatepalm             = 43
nirrig_datepalm       = 44
nfoddergrass          = 45
nirrig_foddergrass    = 46
ngrapes               = 47
nirrig_grapes         = 48
ngroundnuts           = 49
nirrig_groundnuts     = 50
nmillet               = 51
nirrig_millet         = 52
noilpalm              = 53
nirrig_oilpalm        = 54
npotatoes             = 55
nirrig_potatoes       = 56
npulses               = 57
nirrig_pulses         = 58
nrapeseed             = 59
nirrig_rapeseed       = 60
nrice                 = 61
nirrig_rice           = 62
nsorghum              = 63
nirrig_sorghum        = 64
nsugarbeet            = 65
nirrig_sugarbeet      = 66
nsugarcane            = 67
nirrig_sugarcane      = 68
nsunflower            = 69
nirrig_sunflower      = 70
nmiscanthus           = 71
nirrig_miscanthus     = 72
nswitchgrass          = 73
nirrig_switchgrass    = 74
ntrp_corn             = 75  # tropical corn (rf)
nirrig_trp_corn       = 76  # tropical corn (ir)
ntrp_soybean          = 77  # tropical soybean (rf)
nirrig_trp_soybean    = 78  # tropical soybean (ir)
npcropmax             = 0   # last prognostic crop (set during init)
nc3crop               = 15  # generic C3 crop (rf)
nc3irrig              = 16  # irrigated generic crop (ir)

# Number of CFTs known to model (set during init)
num_cfts_known_to_model = 0

# --- Module-level real parameters (Fortran `parameter`) ---
const REINICKERP       = 1.6       # parameter in allometric equation
const DWOOD_PARAM      = 2.5e5     # cn wood density (gC/m3); lpj:2.0e5
const ALLOM1           = 100.0     # parameters in allometric equations
const ALLOM2           = 40.0
const ALLOM3           = 0.5
const ALLOM1S          = 250.0     # modified for shrubs
const ALLOM2S          = 8.0
const ROOT_DENSITY_PARAM = 0.31e6  # root density (g biomass / m3 root)
const ROOT_RADIUS_PARAM  = 0.29e-3 # root radius (m)
const PFTNAME_LEN      = 40       # max length of pftname

# PFT names (Julia 1-based: index 1 = Fortran index 0)
const PFTNAMES = [
    "not_vegetated",                        #  0
    "needleleaf_evergreen_temperate_tree",  #  1
    "needleleaf_evergreen_boreal_tree",     #  2
    "needleleaf_deciduous_boreal_tree",     #  3
    "broadleaf_evergreen_tropical_tree",    #  4
    "broadleaf_evergreen_temperate_tree",   #  5
    "broadleaf_deciduous_tropical_tree",    #  6
    "broadleaf_deciduous_temperate_tree",   #  7
    "broadleaf_deciduous_boreal_tree",      #  8
    "broadleaf_evergreen_shrub",            #  9
    "broadleaf_deciduous_temperate_shrub",  # 10
    "broadleaf_deciduous_boreal_shrub",     # 11
    "c3_arctic_grass",                      # 12
    "c3_non-arctic_grass",                  # 13
    "c4_grass",                             # 14
    "c3_crop",                              # 15
    "c3_irrigated",                         # 16
    "temperate_corn",                       # 17
    "irrigated_temperate_corn",             # 18
    "spring_wheat",                         # 19
    "irrigated_spring_wheat",               # 20
    "winter_wheat",                         # 21
    "irrigated_winter_wheat",               # 22
    "temperate_soybean",                    # 23
    "irrigated_temperate_soybean",          # 24
    "barley",                               # 25
    "irrigated_barley",                     # 26
    "winter_barley",                        # 27
    "irrigated_winter_barley",              # 28
    "rye",                                  # 29
    "irrigated_rye",                        # 30
    "winter_rye",                           # 31
    "irrigated_winter_rye",                 # 32
    "cassava",                              # 33
    "irrigated_cassava",                    # 34
    "citrus",                               # 35
    "irrigated_citrus",                     # 36
    "cocoa",                                # 37
    "irrigated_cocoa",                      # 38
    "coffee",                               # 39
    "irrigated_coffee",                     # 40
    "cotton",                               # 41
    "irrigated_cotton",                     # 42
    "datepalm",                             # 43
    "irrigated_datepalm",                   # 44
    "foddergrass",                          # 45
    "irrigated_foddergrass",                # 46
    "grapes",                               # 47
    "irrigated_grapes",                     # 48
    "groundnuts",                           # 49
    "irrigated_groundnuts",                 # 50
    "millet",                               # 51
    "irrigated_millet",                     # 52
    "oilpalm",                              # 53
    "irrigated_oilpalm",                    # 54
    "potatoes",                             # 55
    "irrigated_potatoes",                   # 56
    "pulses",                               # 57
    "irrigated_pulses",                     # 58
    "rapeseed",                             # 59
    "irrigated_rapeseed",                   # 60
    "rice",                                 # 61
    "irrigated_rice",                       # 62
    "sorghum",                              # 63
    "irrigated_sorghum",                    # 64
    "sugarbeet",                            # 65
    "irrigated_sugarbeet",                  # 66
    "sugarcane",                            # 67
    "irrigated_sugarcane",                  # 68
    "sunflower",                            # 69
    "irrigated_sunflower",                  # 70
    "miscanthus",                           # 71
    "irrigated_miscanthus",                 # 72
    "switchgrass",                          # 73
    "irrigated_switchgrass",                # 74
    "tropical_corn",                        # 75
    "irrigated_tropical_corn",              # 76
    "tropical_soybean",                     # 77
    "irrigated_tropical_soybean",           # 78
]

# --- PFT constants data structure ---
# All 1D arrays are size (MXPFT+1) corresponding to Fortran 0:mxpft.
# Julia index i corresponds to Fortran PFT index (i-1).
Base.@kwdef mutable struct PftconType
    # Vegetation type flags
    noveg::Vector{Int}           = Int[]
    is_tree::Vector{Bool}        = Bool[]
    is_shrub::Vector{Bool}       = Bool[]
    is_grass::Vector{Bool}       = Bool[]

    # Optical/structural properties
    dleaf::Vector{Float64}       = Float64[]   # characteristic leaf dimension (m)
    c3psn::Vector{Float64}       = Float64[]   # photosynthetic pathway: 0=c4, 1=c3
    xl::Vector{Float64}          = Float64[]   # leaf/stem orientation index
    rhol::Matrix{Float64}        = Matrix{Float64}(undef, 0, 0)  # leaf reflectance (npft, numrad)
    rhos::Matrix{Float64}        = Matrix{Float64}(undef, 0, 0)  # stem reflectance (npft, numrad)
    taul::Matrix{Float64}        = Matrix{Float64}(undef, 0, 0)  # leaf transmittance (npft, numrad)
    taus::Matrix{Float64}        = Matrix{Float64}(undef, 0, 0)  # stem transmittance (npft, numrad)

    # Roughness parameters
    z0mr::Vector{Float64}        = Float64[]   # momentum roughness length ratio
    z0v_Cr::Vector{Float64}      = Float64[]   # Raupach92 drag coefficient
    z0v_Cs::Vector{Float64}      = Float64[]   # Raupach92 substrate drag coefficient
    z0v_c::Vector{Float64}       = Float64[]   # Raupach92 c parameter
    z0v_cw::Vector{Float64}      = Float64[]   # Raupach92 roughness sublayer coefficient
    z0v_LAIoff::Vector{Float64}  = Float64[]   # Raupach92 LAI offset
    z0v_LAImax::Vector{Float64}  = Float64[]   # Raupach92 LAI max
    displar::Vector{Float64}     = Float64[]   # displacement height ratio

    # Root parameters
    roota_par::Vector{Float64}   = Float64[]   # CLM rooting distribution parameter [1/m]
    rootb_par::Vector{Float64}   = Float64[]   # CLM rooting distribution parameter [1/m]

    # Crop/irrigation flags
    crop::Vector{Float64}        = Float64[]   # crop pft: 0=not crop, 1=crop
    irrigated::Vector{Float64}   = Float64[]   # irrigated: 0=not, 1=irrigated

    # Soil water potential
    smpso::Vector{Float64}       = Float64[]   # soil water potential at full stomatal opening (mm)
    smpsc::Vector{Float64}       = Float64[]   # soil water potential at full stomatal closure (mm)
    fnitr::Vector{Float64}       = Float64[]   # foliage nitrogen limitation factor

    # CN code
    dwood::Vector{Float64}       = Float64[]   # wood density (gC/m3)
    slatop::Vector{Float64}      = Float64[]   # SLA at top of canopy [m^2/gC]
    dsladlai::Vector{Float64}    = Float64[]   # dSLA/dLAI [m^2/gC]
    leafcn::Vector{Float64}      = Float64[]   # leaf C:N [gC/gN]
    biofuel_harvfrac::Vector{Float64} = Float64[]  # fraction harvested for biofuels
    repr_structure_harvfrac::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # (npft, repr_structure)
    flnr::Vector{Float64}        = Float64[]   # fraction of leaf N in Rubisco
    woody::Vector{Float64}       = Float64[]   # woody lifeform flag (0 or 1)
    lflitcn::Vector{Float64}     = Float64[]   # leaf litter C:N (gC/gN)
    frootcn::Vector{Float64}     = Float64[]   # fine root C:N (gC/gN)
    livewdcn::Vector{Float64}    = Float64[]   # live wood C:N (gC/gN)
    deadwdcn::Vector{Float64}    = Float64[]   # dead wood C:N (gC/gN)
    grperc::Vector{Float64}      = Float64[]   # growth respiration parameter
    grpnow::Vector{Float64}      = Float64[]   # growth respiration parameter
    rootprof_beta::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # rooting distribution (npft, nvariants)
    root_radius::Vector{Float64}   = Float64[] # root radius (m)
    root_density::Vector{Float64}  = Float64[] # root density (gC/m3)

    # Tree structure
    dbh::Vector{Float64}         = Float64[]   # diameter at breast height (m)
    fbw::Vector{Float64}         = Float64[]   # fraction of biomass that is water
    nstem::Vector{Float64}       = Float64[]   # stem density (#/m2)
    taper::Vector{Float64}       = Float64[]   # tapering ratio height:radius_breast_height
    rstem_per_dbh::Vector{Float64} = Float64[] # stem resistance per dbh (s/m/m)
    wood_density::Vector{Float64}  = Float64[] # wood density (kg/m3)
    crit_onset_gdd_sf::Vector{Float64} = Float64[]  # scale factor for crit_onset_gdd
    ndays_on::Vector{Float64}    = Float64[]   # number of days to complete leaf onset

    # Crop merge info
    mergetoclmpft::Vector{Int}   = Int[]
    is_pft_known_to_model::Vector{Bool} = Bool[]

    # Crop parameters
    graincn::Vector{Float64}     = Float64[]   # grain C:N (gC/gN)
    mxtmp::Vector{Float64}       = Float64[]   # parameter used in accFlds
    baset::Vector{Float64}       = Float64[]   # parameter used in accFlds
    declfact::Vector{Float64}    = Float64[]   # parameter used in CNAllocation
    bfact::Vector{Float64}       = Float64[]   # parameter used in CNAllocation
    aleaff::Vector{Float64}      = Float64[]   # parameter used in CNAllocation
    arootf::Vector{Float64}      = Float64[]   # parameter used in CNAllocation
    astemf::Vector{Float64}      = Float64[]   # parameter used in CNAllocation
    arooti::Vector{Float64}      = Float64[]   # parameter used in CNAllocation
    fleafi::Vector{Float64}      = Float64[]   # parameter used in CNAllocation
    allconsl::Vector{Float64}    = Float64[]   # parameter used in CNAllocation
    allconss::Vector{Float64}    = Float64[]   # parameter used in CNAllocation
    ztopmx::Vector{Float64}      = Float64[]   # parameter used in CNVegStructUpdate
    laimx::Vector{Float64}       = Float64[]   # parameter used in CNVegStructUpdate
    gddmin::Vector{Float64}      = Float64[]   # parameter used in CNPhenology
    hybgdd::Vector{Float64}      = Float64[]   # parameter used in CNPhenology
    lfemerg::Vector{Float64}     = Float64[]   # parameter used in CNPhenology
    grnfill::Vector{Float64}     = Float64[]   # parameter used in CNPhenology
    mxmat::Vector{Int}           = Int[]       # parameter used in CNPhenology
    mbbopt::Vector{Float64}      = Float64[]   # Ball-Berry slope
    medlynslope::Vector{Float64} = Float64[]   # Medlyn slope
    medlynintercept::Vector{Float64} = Float64[]  # Medlyn intercept
    mnNHplantdate::Vector{Int}   = Int[]       # min NH planting date (YYYYMMDD)
    mxNHplantdate::Vector{Int}   = Int[]       # max NH planting date (YYYYMMDD)
    mnSHplantdate::Vector{Int}   = Int[]       # min SH planting date (YYYYMMDD)
    mxSHplantdate::Vector{Int}   = Int[]       # max SH planting date (YYYYMMDD)
    planttemp::Vector{Float64}   = Float64[]   # planting temperature (K)
    minplanttemp::Vector{Float64} = Float64[]  # minimum planting temperature (K)

    # Allocation parameters
    froot_leaf::Vector{Float64}  = Float64[]   # new fine root C per new leaf C (gC/gC)
    stem_leaf::Vector{Float64}   = Float64[]   # new stem C per new leaf C (gC/gC)
    croot_stem::Vector{Float64}  = Float64[]   # new coarse root C per new stem C (gC/gC)
    flivewd::Vector{Float64}     = Float64[]   # fraction of new wood that is live
    fcur::Vector{Float64}        = Float64[]   # fraction of allocation to current growth
    fcurdv::Vector{Float64}      = Float64[]   # alternate fcur for cndv

    # Litter fractions
    lf_f::Matrix{Float64}        = Matrix{Float64}(undef, 0, 0)  # leaf litter fractions (npft, ndecomp_pools)
    fr_f::Matrix{Float64}        = Matrix{Float64}(undef, 0, 0)  # fine root litter fractions (npft, ndecomp_pools)
    lf_flab::Vector{Float64}     = Float64[]   # leaf litter labile fraction
    lf_fcel::Vector{Float64}     = Float64[]   # leaf litter cellulose fraction
    lf_flig::Vector{Float64}     = Float64[]   # leaf litter lignin fraction
    fr_flab::Vector{Float64}     = Float64[]   # fine root litter labile fraction
    fr_fcel::Vector{Float64}     = Float64[]   # fine root litter cellulose fraction
    fr_flig::Vector{Float64}     = Float64[]   # fine root litter lignin fraction

    # Phenology flags
    leaf_long::Vector{Float64}   = Float64[]   # leaf longevity (yrs)
    evergreen::Vector{Float64}   = Float64[]   # binary flag for evergreen (0 or 1)
    stress_decid::Vector{Float64} = Float64[]  # binary flag for stress-deciduous (0 or 1)
    season_decid::Vector{Float64} = Float64[]  # binary flag for seasonal-deciduous (0 or 1)
    season_decid_temperate::Vector{Float64} = Float64[]  # binary flag for seasonal-deciduous temperate

    # Product pools
    pconv::Vector{Float64}       = Float64[]   # proportion of deadstem to conversion flux
    pprod10::Vector{Float64}     = Float64[]   # proportion of deadstem to 10-yr product pool
    pprod100::Vector{Float64}    = Float64[]   # proportion of deadstem to 100-yr product pool
    pprodharv10::Vector{Float64} = Float64[]   # harvest mortality proportion to 10-yr pool

    # Fire parameters
    cc_leaf::Vector{Float64}     = Float64[]
    cc_lstem::Vector{Float64}    = Float64[]
    cc_dstem::Vector{Float64}    = Float64[]
    cc_other::Vector{Float64}    = Float64[]
    fm_leaf::Vector{Float64}     = Float64[]
    fm_lstem::Vector{Float64}    = Float64[]
    fm_dstem::Vector{Float64}    = Float64[]
    fm_other::Vector{Float64}    = Float64[]
    fm_root::Vector{Float64}     = Float64[]
    fm_lroot::Vector{Float64}    = Float64[]
    fm_droot::Vector{Float64}    = Float64[]
    fsr_pft::Vector{Float64}     = Float64[]
    fd_pft::Vector{Float64}      = Float64[]
    rswf_min::Vector{Float64}    = Float64[]
    rswf_max::Vector{Float64}    = Float64[]

    # Crop nitrogen parameters
    manunitro::Vector{Float64}   = Float64[]   # manure
    fleafcn::Vector{Float64}     = Float64[]   # C:N during grain fill; leaf
    ffrootcn::Vector{Float64}    = Float64[]   # C:N during grain fill; fine root
    fstemcn::Vector{Float64}     = Float64[]   # C:N during grain fill; stem

    # CLM5 flexible nitrogen
    i_vcad::Vector{Float64}      = Float64[]
    s_vcad::Vector{Float64}      = Float64[]
    i_flnr::Vector{Float64}     = Float64[]
    s_flnr::Vector{Float64}     = Float64[]

    # CNDV parameters (from LPJ)
    pftpar20::Vector{Float64}    = Float64[]   # tree max crown area (m2)
    pftpar28::Vector{Float64}    = Float64[]   # min coldest monthly mean temp
    pftpar29::Vector{Float64}    = Float64[]   # max coldest monthly mean temp
    pftpar30::Vector{Float64}    = Float64[]   # min growing degree days (>= 5 deg C)
    pftpar31::Vector{Float64}    = Float64[]   # upper limit warmest month temp (twmax)

    # FUN parameters
    a_fix::Vector{Float64}         = Float64[] # BNF parameter
    b_fix::Vector{Float64}         = Float64[] # BNF parameter
    c_fix::Vector{Float64}         = Float64[] # BNF parameter
    s_fix::Vector{Float64}         = Float64[] # BNF parameter
    akc_active::Vector{Float64}    = Float64[] # mycorrhizal uptake parameter
    akn_active::Vector{Float64}    = Float64[] # mycorrhizal uptake parameter
    ekc_active::Vector{Float64}    = Float64[] # mycorrhizal uptake parameter
    ekn_active::Vector{Float64}    = Float64[] # mycorrhizal uptake parameter
    kc_nonmyc::Vector{Float64}     = Float64[] # non-mycorrhizal uptake parameter
    kn_nonmyc::Vector{Float64}     = Float64[] # non-mycorrhizal uptake parameter
    kr_resorb::Vector{Float64}     = Float64[] # retranslocation parameter
    perecm::Vector{Float64}        = Float64[] # fraction of ECM-associated PFT
    fun_cn_flex_a::Vector{Float64} = Float64[] # FUN-flexcn parameter a
    fun_cn_flex_b::Vector{Float64} = Float64[] # FUN-flexcn parameter b
    fun_cn_flex_c::Vector{Float64} = Float64[] # FUN-flexcn parameter c
    FUN_fracfixers::Vector{Float64} = Float64[] # fraction of C for fixation
end

"""
    pftcon_allocate!(p::PftconType; ndecomp_pools=3, nvariants=NVARIANTS,
                     repr_structure_min=1, repr_structure_max=NREPR)

Allocate all arrays in pftcon_type. Arrays are size (MXPFT+1) for 1D
and (MXPFT+1, dim2) for 2D. Corresponds to Fortran InitAllocate.
"""
function pftcon_allocate!(p::PftconType;
                          ndecomp_pools::Int=3,
                          nvariants::Int=NVARIANTS,
                          repr_structure_min::Int=1,
                          repr_structure_max::Int=NREPR)
    n = MXPFT + 1  # size for 0:mxpft in Fortran

    p.noveg          = fill(typemax(Int), n)
    p.is_tree        = fill(false, n)
    p.is_shrub       = fill(false, n)
    p.is_grass       = fill(false, n)

    p.dleaf          = zeros(n)
    p.c3psn          = zeros(n)
    p.xl             = zeros(n)
    p.rhol           = zeros(n, NUMRAD)
    p.rhos           = zeros(n, NUMRAD)
    p.taul           = zeros(n, NUMRAD)
    p.taus           = zeros(n, NUMRAD)
    p.z0mr           = zeros(n)
    p.z0v_Cr         = zeros(n)
    p.z0v_Cs         = zeros(n)
    p.z0v_c          = zeros(n)
    p.z0v_cw         = zeros(n)
    p.z0v_LAIoff     = zeros(n)
    p.z0v_LAImax     = zeros(n)
    p.displar        = zeros(n)
    p.roota_par      = zeros(n)
    p.rootb_par      = zeros(n)
    p.crop           = zeros(n)
    p.mergetoclmpft  = zeros(Int, n)
    p.is_pft_known_to_model = fill(false, n)
    p.irrigated      = zeros(n)
    p.smpso          = zeros(n)
    p.smpsc          = zeros(n)
    p.fnitr          = zeros(n)
    p.slatop         = zeros(n)
    p.dsladlai       = zeros(n)
    p.leafcn         = zeros(n)
    p.biofuel_harvfrac = zeros(n)

    nrepr = max(0, repr_structure_max - repr_structure_min + 1)
    p.repr_structure_harvfrac = zeros(n, nrepr)

    p.flnr           = zeros(n)
    p.woody          = zeros(n)
    p.lflitcn        = zeros(n)
    p.frootcn        = zeros(n)
    p.livewdcn       = zeros(n)
    p.deadwdcn       = zeros(n)
    p.grperc         = zeros(n)
    p.grpnow         = zeros(n)
    p.rootprof_beta  = zeros(n, nvariants)
    p.graincn        = zeros(n)
    p.mxtmp          = zeros(n)
    p.baset          = zeros(n)
    p.declfact       = zeros(n)
    p.bfact          = zeros(n)
    p.aleaff         = zeros(n)
    p.arootf         = zeros(n)
    p.astemf         = zeros(n)
    p.arooti         = zeros(n)
    p.fleafi         = zeros(n)
    p.allconsl       = zeros(n)
    p.allconss        = zeros(n)
    p.ztopmx         = zeros(n)
    p.laimx          = zeros(n)
    p.gddmin         = zeros(n)
    p.hybgdd         = zeros(n)
    p.lfemerg        = zeros(n)
    p.grnfill        = zeros(n)
    p.mbbopt         = zeros(n)
    p.medlynslope    = zeros(n)
    p.medlynintercept = zeros(n)
    p.mxmat          = zeros(Int, n)
    p.mnNHplantdate  = zeros(Int, n)
    p.mxNHplantdate  = zeros(Int, n)
    p.mnSHplantdate  = zeros(Int, n)
    p.mxSHplantdate  = zeros(Int, n)
    p.planttemp      = zeros(n)
    p.minplanttemp   = zeros(n)
    p.froot_leaf     = zeros(n)
    p.stem_leaf      = zeros(n)
    p.croot_stem     = zeros(n)
    p.flivewd        = zeros(n)
    p.fcur           = zeros(n)
    p.fcurdv         = zeros(n)
    p.lf_f           = zeros(n, ndecomp_pools)
    p.fr_f           = zeros(n, ndecomp_pools)
    p.lf_flab        = zeros(n)
    p.lf_fcel        = zeros(n)
    p.lf_flig        = zeros(n)
    p.fr_flab        = zeros(n)
    p.fr_fcel        = zeros(n)
    p.fr_flig        = zeros(n)
    p.leaf_long      = zeros(n)
    p.evergreen      = zeros(n)
    p.stress_decid   = zeros(n)
    p.season_decid   = zeros(n)
    p.season_decid_temperate = zeros(n)
    p.dwood          = zeros(n)
    p.root_density   = zeros(n)
    p.root_radius    = zeros(n)
    p.pconv          = zeros(n)
    p.pprod10        = zeros(n)
    p.pprod100       = zeros(n)
    p.pprodharv10    = zeros(n)
    p.cc_leaf        = zeros(n)
    p.cc_lstem       = zeros(n)
    p.cc_dstem       = zeros(n)
    p.cc_other       = zeros(n)
    p.fm_leaf        = zeros(n)
    p.fm_lstem       = zeros(n)
    p.fm_dstem       = zeros(n)
    p.fm_other       = zeros(n)
    p.fm_root        = zeros(n)
    p.fm_lroot       = zeros(n)
    p.fm_droot       = zeros(n)
    p.fsr_pft        = zeros(n)
    p.fd_pft         = zeros(n)
    p.rswf_min       = zeros(n)
    p.rswf_max       = zeros(n)
    p.manunitro      = zeros(n)
    p.fleafcn        = zeros(n)
    p.ffrootcn       = zeros(n)
    p.fstemcn        = zeros(n)
    p.i_vcad         = zeros(n)
    p.s_vcad         = zeros(n)
    p.i_flnr         = zeros(n)
    p.s_flnr         = zeros(n)
    p.pftpar20       = zeros(n)
    p.pftpar28       = zeros(n)
    p.pftpar29       = zeros(n)
    p.pftpar30       = zeros(n)
    p.pftpar31       = zeros(n)
    p.a_fix          = zeros(n)
    p.b_fix          = zeros(n)
    p.c_fix          = zeros(n)
    p.s_fix          = zeros(n)
    p.akc_active     = zeros(n)
    p.akn_active     = zeros(n)
    p.ekc_active     = zeros(n)
    p.ekn_active     = zeros(n)
    p.kc_nonmyc      = zeros(n)
    p.kn_nonmyc      = zeros(n)
    p.kr_resorb      = zeros(n)
    p.perecm         = zeros(n)
    p.fun_cn_flex_a  = zeros(n)
    p.fun_cn_flex_b  = zeros(n)
    p.fun_cn_flex_c  = zeros(n)
    p.FUN_fracfixers = zeros(n)

    p.dbh            = zeros(n)
    p.fbw            = zeros(n)
    p.nstem          = zeros(n)
    p.taper          = zeros(n)
    p.rstem_per_dbh  = zeros(n)
    p.wood_density   = zeros(n)
    p.crit_onset_gdd_sf = zeros(n)
    p.ndays_on       = zeros(n)

    nothing
end

"""
    pftcon_init_for_testing!(p::PftconType; kwargs...)

Allocate arrays but don't read from file. Values can then be set by tests.
Corresponds to Fortran InitForTesting.
"""
function pftcon_init_for_testing!(p::PftconType;
                                  ndecomp_pools::Int=3,
                                  nvariants::Int=NVARIANTS,
                                  repr_structure_min::Int=1,
                                  repr_structure_max::Int=NREPR)
    pftcon_allocate!(p; ndecomp_pools=ndecomp_pools,
                     nvariants=nvariants,
                     repr_structure_min=repr_structure_min,
                     repr_structure_max=repr_structure_max)
end

"""
    set_is_pft_known_to_model!(p::PftconType)

Set is_pft_known_to_model based on mergetoclmpft.
Julia index = Fortran index + 1.
"""
function set_is_pft_known_to_model!(p::PftconType)
    p.is_pft_known_to_model .= false

    # Type 0 (Julia index 1) is always known
    p.is_pft_known_to_model[1] = true

    # For m = 1:mxpft (Fortran), Julia index = m+1
    for m in 1:MXPFT
        merge_type = p.mergetoclmpft[m + 1]  # Fortran value (0-based)
        p.is_pft_known_to_model[merge_type + 1] = true  # +1 for Julia indexing
    end
end

"""
    set_num_cfts_known_to_model!(p::PftconType; cft_lb, cft_ub)

Set the module-level variable num_cfts_known_to_model.
cft_lb and cft_ub are Fortran 0-based indices.
"""
function set_num_cfts_known_to_model!(p::PftconType;
                                      cft_lb::Int=varpar.cft_lb,
                                      cft_ub::Int=varpar.cft_ub)
    global num_cfts_known_to_model
    num_cfts_known_to_model = 0
    for m in cft_lb:cft_ub
        if p.is_pft_known_to_model[m + 1]  # +1 for Julia indexing
            num_cfts_known_to_model += 1
        end
    end
end

"""
    pftcon_clean!(p::PftconType)

Reset all arrays to empty. Corresponds to Fortran Clean/deallocate.
"""
function pftcon_clean!(p::PftconType)
    p.noveg          = Int[]
    p.is_tree        = Bool[]
    p.is_shrub       = Bool[]
    p.is_grass       = Bool[]
    p.dleaf          = Float64[]
    p.c3psn          = Float64[]
    p.xl             = Float64[]
    p.rhol           = Matrix{Float64}(undef, 0, 0)
    p.rhos           = Matrix{Float64}(undef, 0, 0)
    p.taul           = Matrix{Float64}(undef, 0, 0)
    p.taus           = Matrix{Float64}(undef, 0, 0)
    p.z0mr           = Float64[]
    p.z0v_Cr         = Float64[]
    p.z0v_Cs         = Float64[]
    p.z0v_c          = Float64[]
    p.z0v_cw         = Float64[]
    p.z0v_LAIoff     = Float64[]
    p.z0v_LAImax     = Float64[]
    p.displar        = Float64[]
    p.roota_par      = Float64[]
    p.rootb_par      = Float64[]
    p.crop           = Float64[]
    p.mergetoclmpft  = Int[]
    p.is_pft_known_to_model = Bool[]
    p.irrigated      = Float64[]
    p.smpso          = Float64[]
    p.smpsc          = Float64[]
    p.fnitr          = Float64[]
    p.slatop         = Float64[]
    p.dsladlai       = Float64[]
    p.leafcn         = Float64[]
    p.biofuel_harvfrac = Float64[]
    p.repr_structure_harvfrac = Matrix{Float64}(undef, 0, 0)
    p.flnr           = Float64[]
    p.woody          = Float64[]
    p.lflitcn        = Float64[]
    p.frootcn        = Float64[]
    p.livewdcn       = Float64[]
    p.deadwdcn       = Float64[]
    p.grperc         = Float64[]
    p.grpnow         = Float64[]
    p.rootprof_beta  = Matrix{Float64}(undef, 0, 0)
    p.graincn        = Float64[]
    p.mxtmp          = Float64[]
    p.baset          = Float64[]
    p.declfact       = Float64[]
    p.bfact          = Float64[]
    p.aleaff         = Float64[]
    p.arootf         = Float64[]
    p.astemf         = Float64[]
    p.arooti         = Float64[]
    p.fleafi         = Float64[]
    p.allconsl       = Float64[]
    p.allconss        = Float64[]
    p.ztopmx         = Float64[]
    p.laimx          = Float64[]
    p.gddmin         = Float64[]
    p.hybgdd         = Float64[]
    p.lfemerg        = Float64[]
    p.grnfill        = Float64[]
    p.mbbopt         = Float64[]
    p.medlynslope    = Float64[]
    p.medlynintercept = Float64[]
    p.mxmat          = Int[]
    p.mnNHplantdate  = Int[]
    p.mxNHplantdate  = Int[]
    p.mnSHplantdate  = Int[]
    p.mxSHplantdate  = Int[]
    p.planttemp      = Float64[]
    p.minplanttemp   = Float64[]
    p.froot_leaf     = Float64[]
    p.stem_leaf      = Float64[]
    p.croot_stem     = Float64[]
    p.flivewd        = Float64[]
    p.fcur           = Float64[]
    p.fcurdv         = Float64[]
    p.lf_f           = Matrix{Float64}(undef, 0, 0)
    p.fr_f           = Matrix{Float64}(undef, 0, 0)
    p.lf_flab        = Float64[]
    p.lf_fcel        = Float64[]
    p.lf_flig        = Float64[]
    p.fr_flab        = Float64[]
    p.fr_fcel        = Float64[]
    p.fr_flig        = Float64[]
    p.leaf_long      = Float64[]
    p.evergreen      = Float64[]
    p.stress_decid   = Float64[]
    p.season_decid   = Float64[]
    p.season_decid_temperate = Float64[]
    p.dwood          = Float64[]
    p.root_density   = Float64[]
    p.root_radius    = Float64[]
    p.pconv          = Float64[]
    p.pprod10        = Float64[]
    p.pprod100       = Float64[]
    p.pprodharv10    = Float64[]
    p.cc_leaf        = Float64[]
    p.cc_lstem       = Float64[]
    p.cc_dstem       = Float64[]
    p.cc_other       = Float64[]
    p.fm_leaf        = Float64[]
    p.fm_lstem       = Float64[]
    p.fm_dstem       = Float64[]
    p.fm_other       = Float64[]
    p.fm_root        = Float64[]
    p.fm_lroot       = Float64[]
    p.fm_droot       = Float64[]
    p.fsr_pft        = Float64[]
    p.fd_pft         = Float64[]
    p.rswf_min       = Float64[]
    p.rswf_max       = Float64[]
    p.manunitro      = Float64[]
    p.fleafcn        = Float64[]
    p.ffrootcn       = Float64[]
    p.fstemcn         = Float64[]
    p.i_vcad         = Float64[]
    p.s_vcad         = Float64[]
    p.i_flnr         = Float64[]
    p.s_flnr         = Float64[]
    p.pftpar20       = Float64[]
    p.pftpar28       = Float64[]
    p.pftpar29       = Float64[]
    p.pftpar30       = Float64[]
    p.pftpar31       = Float64[]
    p.a_fix          = Float64[]
    p.b_fix          = Float64[]
    p.c_fix          = Float64[]
    p.s_fix          = Float64[]
    p.akc_active     = Float64[]
    p.akn_active     = Float64[]
    p.ekc_active     = Float64[]
    p.ekn_active     = Float64[]
    p.kc_nonmyc      = Float64[]
    p.kn_nonmyc      = Float64[]
    p.kr_resorb      = Float64[]
    p.perecm         = Float64[]
    p.fun_cn_flex_a  = Float64[]
    p.fun_cn_flex_b  = Float64[]
    p.fun_cn_flex_c  = Float64[]
    p.FUN_fracfixers = Float64[]
    p.dbh            = Float64[]
    p.fbw            = Float64[]
    p.nstem          = Float64[]
    p.taper          = Float64[]
    p.rstem_per_dbh  = Float64[]
    p.wood_density   = Float64[]
    p.crit_onset_gdd_sf = Float64[]
    p.ndays_on       = Float64[]

    nothing
end

# Global pftcon instance
const pftcon = PftconType()
