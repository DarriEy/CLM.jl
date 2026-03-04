# ==========================================================================
# Ported from: src/biogeophys/PhotosynthesisMod.F90 (5210 lines)
# Leaf photosynthesis and stomatal conductance calculation as described by
# Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 and extended to
# a multi-layer canopy.
#
# Public subroutines:
#   photosynthesis!            — Leaf stomatal resistance and leaf photosynthesis
#   photosynthesis_total!      — Determine total photosynthesis
#   fractionation!             — C13 fractionation during photosynthesis
#   photosynthesis_hydrstress! — Hydraulic stress photosynthesis (PHS)
#   plc                        — Vulnerability curve value
#   d1plc                      — 1st derivative of vulnerability curve
#
# Private functions:
#   hybrid_solver!       — Hybrid solver for ci
#   ci_func!             — ci function evaluation
#   brent_solver!        — Brent solver for root finding
#   ft_photo             — Photosynthesis temperature response
#   fth_photo            — Photosynthesis temperature inhibition
#   fth25_photo          — Scaling factor for temperature inhibition
#   quadratic_solve      — Solve quadratic equation
#   hybrid_PHS!          — Hybrid solver for PHS ci
#   ci_func_PHS!         — ci function for PHS
#   brent_PHS!           — Brent solver for PHS
#   calcstress!          — Compute transpiration stress via plant hydraulics
#   getvegwp!            — Calculate vegetation water potential
#   getqflx!             — Calculate sunlit and shaded transpiration
#   spacF!               — Flux divergence across vegetation segments
#   spacA!               — Inverse Jacobian for vegetation water potential
# ==========================================================================

# --- Module-level constants ---

# Leaf respiration methods
const LEAFRESP_MTD_RYAN1991  = 1
const LEAFRESP_MTD_ATKIN2015 = 2

# PLC method type
const VEGETATION_WEIBULL = 0

# Segment indices (public for unit tests)
const SUN  = 1
const SHA  = 2
const XYL  = 3
const ROOT_SEG = 4   # renamed from 'root' to avoid clash
const VEG  = VEGETATION_WEIBULL
const SOIL_SEG = 1   # renamed from 'soil'

# Stomatal conductance method types
const STOMATALCOND_MTD_BB1987     = 1
const STOMATALCOND_MTD_MEDLYN2011 = 2

# Ball-Berry intercepts
const BBBOPT_C3 = 10000.0
const BBBOPT_C4 = 40000.0

# Medlyn parameters
const MEDLYN_RH_CAN_MAX  = 50.0
const MEDLYN_RH_CAN_FACT = 0.001

# Max CO2 partial pressure at leaf surface for PHS
const MAX_CS = 1.0e-06

# =====================================================================
# Photosynthesis parameters type
# =====================================================================

"""
    PhotoParamsData

Photosynthesis parameters (from `photo_params_type` in Fortran).
"""
Base.@kwdef mutable struct PhotoParamsData
    act25::Float64 = NaN         # Rubisco activity at 25 C (umol CO2/gRubisco/s)
    fnr::Float64 = NaN           # Mass ratio of total Rubisco molecular mass to nitrogen in Rubisco
    cp25_yr2000::Float64 = NaN   # CO2 compensation point at 25°C at present day O2 (mol/mol)
    kc25_coef::Float64 = NaN     # Michaelis-Menten const. at 25°C for CO2
    ko25_coef::Float64 = NaN     # Michaelis-Menten const. at 25°C for O2
    fnps::Float64 = NaN          # Fraction of light absorbed by non-photosynthetic pigment
    theta_psii::Float64 = NaN    # Empirical curvature parameter for electron transport rate
    theta_ip::Float64 = NaN      # Empirical curvature parameter for ap photosynthesis co-limitation
    vcmaxha::Float64 = NaN       # Activation energy for vcmax (J/mol)
    jmaxha::Float64 = NaN        # Activation energy for jmax (J/mol)
    tpuha::Float64 = NaN         # Activation energy for tpu (J/mol)
    lmrha::Float64 = NaN         # Activation energy for lmr (J/mol)
    kcha::Float64 = NaN          # Activation energy for kc (J/mol)
    koha::Float64 = NaN          # Activation energy for ko (J/mol)
    cpha::Float64 = NaN          # Activation energy for cp (J/mol)
    vcmaxhd::Float64 = NaN       # Deactivation energy for vcmax (J/mol)
    jmaxhd::Float64 = NaN        # Deactivation energy for jmax (J/mol)
    tpuhd::Float64 = NaN         # Deactivation energy for tpu (J/mol)
    lmrhd::Float64 = NaN         # Deactivation energy for lmr (J/mol)
    lmrse::Float64 = NaN         # Entropy term for lmr (J/mol/K)
    tpu25ratio::Float64 = NaN    # Ratio of tpu25top to vcmax25top
    kp25ratio::Float64 = NaN     # Ratio of kp25top to vcmax25top
    vcmaxse_sf::Float64 = NaN    # Scale factor for vcmaxse
    jmaxse_sf::Float64 = NaN     # Scale factor for jmaxse
    tpuse_sf::Float64 = NaN      # Scale factor for tpuse
    jmax25top_sf::Float64 = NaN  # Scale factor for jmax25top
    krmax::Vector{Float64} = Float64[]                  # (0:mxpft)
    kmax::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # (0:mxpft, nvegwcs)
    psi50::Matrix{Float64} = Matrix{Float64}(undef, 0, 0) # (0:mxpft, nvegwcs)
    ck::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)    # (0:mxpft, nvegwcs)
    lmr_intercept_atkin::Vector{Float64} = Float64[]    # (0:mxpft)
    theta_cj::Vector{Float64} = Float64[]               # (0:mxpft)
end

# Global params instance
const params_inst = PhotoParamsData()

"""
    photo_params_init!(pp::PhotoParamsData, mxpft::Int, nvegwcs::Int)

Allocate photosynthesis parameter arrays.
"""
function photo_params_init!(pp::PhotoParamsData, mxpft::Int=MXPFT, nvegwcs::Int=NVEGWCS)
    pp.krmax = fill(NaN, mxpft + 1)
    pp.theta_cj = fill(NaN, mxpft + 1)
    pp.kmax = fill(NaN, mxpft + 1, nvegwcs)
    pp.psi50 = fill(NaN, mxpft + 1, nvegwcs)
    pp.ck = fill(NaN, mxpft + 1, nvegwcs)
    pp.lmr_intercept_atkin = fill(NaN, mxpft + 1)
    return nothing
end

"""
    photo_params_clean!(pp::PhotoParamsData)

Deallocate photosynthesis parameter arrays.
"""
function photo_params_clean!(pp::PhotoParamsData)
    pp.krmax = Float64[]
    pp.theta_cj = Float64[]
    pp.kmax = Matrix{Float64}(undef, 0, 0)
    pp.psi50 = Matrix{Float64}(undef, 0, 0)
    pp.ck = Matrix{Float64}(undef, 0, 0)
    pp.lmr_intercept_atkin = Float64[]
    return nothing
end

# =====================================================================
# Photosynthesis data type
# =====================================================================

"""
    PhotosynthesisData

Photosynthesis state data (from `photosyns_type` in Fortran).
Holds all photosynthesis and stomatal conductance variables at patch level.
"""
Base.@kwdef mutable struct PhotosynthesisData
    # Logical/config
    c3flag_patch::Vector{Bool} = Bool[]

    # PHS-specific 3D arrays (np, 2, nlevcan)
    ac_phs_patch::Array{Float64,3} = Array{Float64}(undef, 0, 0, 0)
    aj_phs_patch::Array{Float64,3} = Array{Float64}(undef, 0, 0, 0)
    ap_phs_patch::Array{Float64,3} = Array{Float64}(undef, 0, 0, 0)
    ag_phs_patch::Array{Float64,3} = Array{Float64}(undef, 0, 0, 0)
    vcmax_z_phs_patch::Array{Float64,3} = Array{Float64}(undef, 0, 0, 0)
    tpu_z_phs_patch::Array{Float64,3} = Array{Float64}(undef, 0, 0, 0)
    kp_z_phs_patch::Array{Float64,3} = Array{Float64}(undef, 0, 0, 0)

    # Sunlit/shaded net photosynthesis (np, nlevcan)
    an_sun_patch::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    an_sha_patch::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)

    # Stomatal conductance sunlit/shaded (np, nlevcan)
    gs_mol_sun_patch::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    gs_mol_sha_patch::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    gs_mol_sun_ln_patch::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    gs_mol_sha_ln_patch::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)

    # Standard (non-PHS) 2D arrays (np, nlevcan)
    ac_patch::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    aj_patch::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    ap_patch::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    ag_patch::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    an_patch::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    vcmax_z_patch::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    tpu_z_patch::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    kp_z_patch::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    gs_mol_patch::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)

    # 1D patch-level
    cp_patch::Vector{Float64} = Float64[]
    kc_patch::Vector{Float64} = Float64[]
    ko_patch::Vector{Float64} = Float64[]
    qe_patch::Vector{Float64} = Float64[]
    bbb_patch::Vector{Float64} = Float64[]
    mbb_patch::Vector{Float64} = Float64[]
    gb_mol_patch::Vector{Float64} = Float64[]
    rh_leaf_patch::Vector{Float64} = Float64[]
    vpd_can_patch::Vector{Float64} = Float64[]

    # Photosynthesis outputs
    psnsun_patch::Vector{Float64} = Float64[]
    psnsha_patch::Vector{Float64} = Float64[]
    c13_psnsun_patch::Vector{Float64} = Float64[]
    c13_psnsha_patch::Vector{Float64} = Float64[]
    c14_psnsun_patch::Vector{Float64} = Float64[]
    c14_psnsha_patch::Vector{Float64} = Float64[]

    psnsun_z_patch::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    psnsha_z_patch::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    psnsun_wc_patch::Vector{Float64} = Float64[]
    psnsha_wc_patch::Vector{Float64} = Float64[]
    psnsun_wj_patch::Vector{Float64} = Float64[]
    psnsha_wj_patch::Vector{Float64} = Float64[]
    psnsun_wp_patch::Vector{Float64} = Float64[]
    psnsha_wp_patch::Vector{Float64} = Float64[]

    fpsn_patch::Vector{Float64} = Float64[]
    fpsn_wc_patch::Vector{Float64} = Float64[]
    fpsn_wj_patch::Vector{Float64} = Float64[]
    fpsn_wp_patch::Vector{Float64} = Float64[]

    lnca_patch::Vector{Float64} = Float64[]

    lmrsun_patch::Vector{Float64} = Float64[]
    lmrsha_patch::Vector{Float64} = Float64[]
    lmrsun_z_patch::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    lmrsha_z_patch::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)

    alphapsnsun_patch::Vector{Float64} = Float64[]
    alphapsnsha_patch::Vector{Float64} = Float64[]
    rc13_canair_patch::Vector{Float64} = Float64[]
    rc13_psnsun_patch::Vector{Float64} = Float64[]
    rc13_psnsha_patch::Vector{Float64} = Float64[]

    cisun_z_patch::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    cisha_z_patch::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)

    rssun_z_patch::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    rssha_z_patch::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    rssun_patch::Vector{Float64} = Float64[]
    rssha_patch::Vector{Float64} = Float64[]

    luvcmax25top_patch::Vector{Float64} = Float64[]
    lujmax25top_patch::Vector{Float64} = Float64[]
    lutpu25top_patch::Vector{Float64} = Float64[]

    # LUNA-specific
    vcmx25_z_patch::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    jmx25_z_patch::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    vcmx25_z_last_valid_patch::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    jmx25_z_last_valid_patch::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    pnlc_z_patch::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    enzs_z_patch::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    fpsn24_patch::Vector{Float64} = Float64[]

    # Configuration switches
    rootstem_acc::Bool = false
    light_inhibit::Bool = false
    leafresp_method::Int = LEAFRESP_MTD_RYAN1991
    stomatalcond_mtd::Int = STOMATALCOND_MTD_BB1987
    modifyphoto_and_lmr_forcrop::Bool = false
end

"""
    photosynthesis_data_init!(ps::PhotosynthesisData, np::Int; nlevcan::Int=NLEVCAN, use_luna::Bool=false)

Allocate and initialize all fields for `np` patches.
"""
function photosynthesis_data_init!(ps::PhotosynthesisData, np::Int;
                                   nlevcan::Int=NLEVCAN, use_luna::Bool=false)
    ps.c3flag_patch = fill(false, np)

    ps.ac_phs_patch = fill(NaN, np, 2, nlevcan)
    ps.aj_phs_patch = fill(NaN, np, 2, nlevcan)
    ps.ap_phs_patch = fill(NaN, np, 2, nlevcan)
    ps.ag_phs_patch = fill(NaN, np, 2, nlevcan)
    ps.vcmax_z_phs_patch = fill(NaN, np, 2, nlevcan)
    ps.tpu_z_phs_patch = fill(NaN, np, 2, nlevcan)
    ps.kp_z_phs_patch = fill(NaN, np, 2, nlevcan)

    ps.an_sun_patch = fill(NaN, np, nlevcan)
    ps.an_sha_patch = fill(NaN, np, nlevcan)

    ps.gs_mol_sun_patch = fill(NaN, np, nlevcan)
    ps.gs_mol_sha_patch = fill(NaN, np, nlevcan)
    ps.gs_mol_sun_ln_patch = fill(NaN, np, nlevcan)
    ps.gs_mol_sha_ln_patch = fill(NaN, np, nlevcan)

    ps.ac_patch = fill(NaN, np, nlevcan)
    ps.aj_patch = fill(NaN, np, nlevcan)
    ps.ap_patch = fill(NaN, np, nlevcan)
    ps.ag_patch = fill(NaN, np, nlevcan)
    ps.an_patch = fill(NaN, np, nlevcan)
    ps.vcmax_z_patch = fill(NaN, np, nlevcan)
    ps.tpu_z_patch = fill(NaN, np, nlevcan)
    ps.kp_z_patch = fill(NaN, np, nlevcan)
    ps.gs_mol_patch = fill(NaN, np, nlevcan)

    ps.cp_patch = fill(NaN, np)
    ps.kc_patch = fill(NaN, np)
    ps.ko_patch = fill(NaN, np)
    ps.qe_patch = fill(NaN, np)
    ps.bbb_patch = fill(NaN, np)
    ps.mbb_patch = fill(NaN, np)
    ps.gb_mol_patch = fill(NaN, np)
    ps.rh_leaf_patch = fill(NaN, np)
    ps.vpd_can_patch = fill(NaN, np)

    ps.psnsun_patch = fill(NaN, np)
    ps.psnsha_patch = fill(NaN, np)
    ps.c13_psnsun_patch = fill(NaN, np)
    ps.c13_psnsha_patch = fill(NaN, np)
    ps.c14_psnsun_patch = fill(NaN, np)
    ps.c14_psnsha_patch = fill(NaN, np)

    ps.psnsun_z_patch = fill(NaN, np, nlevcan)
    ps.psnsha_z_patch = fill(NaN, np, nlevcan)
    ps.psnsun_wc_patch = fill(NaN, np)
    ps.psnsha_wc_patch = fill(NaN, np)
    ps.psnsun_wj_patch = fill(NaN, np)
    ps.psnsha_wj_patch = fill(NaN, np)
    ps.psnsun_wp_patch = fill(NaN, np)
    ps.psnsha_wp_patch = fill(NaN, np)

    ps.fpsn_patch = fill(NaN, np)
    ps.fpsn_wc_patch = fill(NaN, np)
    ps.fpsn_wj_patch = fill(NaN, np)
    ps.fpsn_wp_patch = fill(NaN, np)

    ps.lnca_patch = fill(NaN, np)

    ps.lmrsun_z_patch = fill(NaN, np, nlevcan)
    ps.lmrsha_z_patch = fill(NaN, np, nlevcan)
    ps.lmrsun_patch = fill(NaN, np)
    ps.lmrsha_patch = fill(NaN, np)

    ps.alphapsnsun_patch = fill(NaN, np)
    ps.alphapsnsha_patch = fill(NaN, np)
    ps.rc13_canair_patch = fill(NaN, np)
    ps.rc13_psnsun_patch = fill(NaN, np)
    ps.rc13_psnsha_patch = fill(NaN, np)

    ps.cisun_z_patch = fill(NaN, np, nlevcan)
    ps.cisha_z_patch = fill(NaN, np, nlevcan)

    ps.rssun_z_patch = fill(2.0e4, np, nlevcan)   # large default (closed stomata)
    ps.rssha_z_patch = fill(2.0e4, np, nlevcan)
    ps.rssun_patch = fill(2.0e4, np)              # large default (closed stomata)
    ps.rssha_patch = fill(2.0e4, np)

    ps.luvcmax25top_patch = fill(NaN, np)
    ps.lujmax25top_patch = fill(NaN, np)
    ps.lutpu25top_patch = fill(NaN, np)

    if use_luna
        ps.vcmx25_z_patch = fill(30.0, np, nlevcan)
        ps.jmx25_z_patch = fill(60.0, np, nlevcan)
        ps.vcmx25_z_last_valid_patch = fill(30.0, np, nlevcan)
        ps.jmx25_z_last_valid_patch = fill(60.0, np, nlevcan)
        ps.pnlc_z_patch = fill(0.01, np, nlevcan)
        ps.fpsn24_patch = fill(NaN, np)
        ps.enzs_z_patch = fill(1.0, np, nlevcan)
    end

    return nothing
end

"""
    photosynthesis_data_clean!(ps::PhotosynthesisData)

Reset all fields to empty arrays.
"""
function photosynthesis_data_clean!(ps::PhotosynthesisData)
    ps.c3flag_patch = Bool[]
    ps.ac_phs_patch = Array{Float64}(undef, 0, 0, 0)
    ps.aj_phs_patch = Array{Float64}(undef, 0, 0, 0)
    ps.ap_phs_patch = Array{Float64}(undef, 0, 0, 0)
    ps.ag_phs_patch = Array{Float64}(undef, 0, 0, 0)
    ps.vcmax_z_phs_patch = Array{Float64}(undef, 0, 0, 0)
    ps.tpu_z_phs_patch = Array{Float64}(undef, 0, 0, 0)
    ps.kp_z_phs_patch = Array{Float64}(undef, 0, 0, 0)
    ps.an_sun_patch = Matrix{Float64}(undef, 0, 0)
    ps.an_sha_patch = Matrix{Float64}(undef, 0, 0)
    ps.gs_mol_sun_patch = Matrix{Float64}(undef, 0, 0)
    ps.gs_mol_sha_patch = Matrix{Float64}(undef, 0, 0)
    ps.gs_mol_sun_ln_patch = Matrix{Float64}(undef, 0, 0)
    ps.gs_mol_sha_ln_patch = Matrix{Float64}(undef, 0, 0)
    ps.ac_patch = Matrix{Float64}(undef, 0, 0)
    ps.aj_patch = Matrix{Float64}(undef, 0, 0)
    ps.ap_patch = Matrix{Float64}(undef, 0, 0)
    ps.ag_patch = Matrix{Float64}(undef, 0, 0)
    ps.an_patch = Matrix{Float64}(undef, 0, 0)
    ps.vcmax_z_patch = Matrix{Float64}(undef, 0, 0)
    ps.tpu_z_patch = Matrix{Float64}(undef, 0, 0)
    ps.kp_z_patch = Matrix{Float64}(undef, 0, 0)
    ps.gs_mol_patch = Matrix{Float64}(undef, 0, 0)
    ps.cp_patch = Float64[]
    ps.kc_patch = Float64[]
    ps.ko_patch = Float64[]
    ps.qe_patch = Float64[]
    ps.bbb_patch = Float64[]
    ps.mbb_patch = Float64[]
    ps.gb_mol_patch = Float64[]
    ps.rh_leaf_patch = Float64[]
    ps.vpd_can_patch = Float64[]
    ps.psnsun_patch = Float64[]
    ps.psnsha_patch = Float64[]
    ps.c13_psnsun_patch = Float64[]
    ps.c13_psnsha_patch = Float64[]
    ps.c14_psnsun_patch = Float64[]
    ps.c14_psnsha_patch = Float64[]
    ps.psnsun_z_patch = Matrix{Float64}(undef, 0, 0)
    ps.psnsha_z_patch = Matrix{Float64}(undef, 0, 0)
    ps.psnsun_wc_patch = Float64[]
    ps.psnsha_wc_patch = Float64[]
    ps.psnsun_wj_patch = Float64[]
    ps.psnsha_wj_patch = Float64[]
    ps.psnsun_wp_patch = Float64[]
    ps.psnsha_wp_patch = Float64[]
    ps.fpsn_patch = Float64[]
    ps.fpsn_wc_patch = Float64[]
    ps.fpsn_wj_patch = Float64[]
    ps.fpsn_wp_patch = Float64[]
    ps.lnca_patch = Float64[]
    ps.lmrsun_z_patch = Matrix{Float64}(undef, 0, 0)
    ps.lmrsha_z_patch = Matrix{Float64}(undef, 0, 0)
    ps.lmrsun_patch = Float64[]
    ps.lmrsha_patch = Float64[]
    ps.alphapsnsun_patch = Float64[]
    ps.alphapsnsha_patch = Float64[]
    ps.rc13_canair_patch = Float64[]
    ps.rc13_psnsun_patch = Float64[]
    ps.rc13_psnsha_patch = Float64[]
    ps.cisun_z_patch = Matrix{Float64}(undef, 0, 0)
    ps.cisha_z_patch = Matrix{Float64}(undef, 0, 0)
    ps.rssun_z_patch = Matrix{Float64}(undef, 0, 0)
    ps.rssha_z_patch = Matrix{Float64}(undef, 0, 0)
    ps.rssun_patch = Float64[]
    ps.rssha_patch = Float64[]
    ps.luvcmax25top_patch = Float64[]
    ps.lujmax25top_patch = Float64[]
    ps.lutpu25top_patch = Float64[]
    ps.vcmx25_z_patch = Matrix{Float64}(undef, 0, 0)
    ps.jmx25_z_patch = Matrix{Float64}(undef, 0, 0)
    ps.vcmx25_z_last_valid_patch = Matrix{Float64}(undef, 0, 0)
    ps.jmx25_z_last_valid_patch = Matrix{Float64}(undef, 0, 0)
    ps.pnlc_z_patch = Matrix{Float64}(undef, 0, 0)
    ps.enzs_z_patch = Matrix{Float64}(undef, 0, 0)
    ps.fpsn24_patch = Float64[]
    return nothing
end

# =====================================================================
# Helper functions: temperature response, quadratic solver
# =====================================================================

"""
    ft_photo(tl, ha)

Photosynthesis temperature response function.
"""
@inline function ft_photo(tl::Float64, ha::Float64)
    return exp(ha / (RGAS * 1.0e-3 * (TFRZ + 25.0)) * (1.0 - (TFRZ + 25.0) / tl))
end

"""
    fth_photo(tl, hd, se, scaleFactor)

Photosynthesis temperature inhibition function.
"""
@inline function fth_photo(tl::Float64, hd::Float64, se::Float64, scaleFactor::Float64)
    return scaleFactor / (1.0 + exp((-hd + se * tl) / (RGAS * 1.0e-3 * tl)))
end

"""
    fth25_photo(hd, se)

Scaling factor for photosynthesis temperature inhibition at 25°C.
"""
@inline function fth25_photo(hd::Float64, se::Float64)
    return 1.0 + exp((-hd + se * (TFRZ + 25.0)) / (RGAS * 1.0e-3 * (TFRZ + 25.0)))
end

"""
    quadratic_solve(a, b, c) -> (r1, r2)

Solve the quadratic equation `a*x^2 + b*x + c = 0` for real roots.
Returns the two roots `(r1, r2)` where `r1 >= r2`.
"""
function quadratic_solve(a::Float64, b::Float64, c::Float64)
    if a == 0.0
        if b == 0.0
            return (0.0, 0.0)
        else
            r1 = -c / b
            return (r1, r1)
        end
    end
    discriminant = b * b - 4.0 * a * c
    if discriminant < 0.0
        discriminant = 0.0
    end
    q = -0.5 * (b + sign(b) * sqrt(discriminant))
    if q == 0.0
        r1 = 0.0
        r2 = 0.0
    else
        r1 = q / a
        r2 = c / q
    end
    if r1 < r2
        r1, r2 = r2, r1
    end
    return (r1, r2)
end

"""
    plc(x, ivt, level, plc_method, params)

Return value of vulnerability curve at water potential `x`.
`ivt` is 1-based PFT index, `level` is segment index (SUN, SHA, XYL, ROOT_SEG).
"""
function plc(x::Float64, ivt::Int, level::Int, plc_method::Int, params::PhotoParamsData=params_inst)
    if plc_method == VEGETATION_WEIBULL
        val = 2.0^(-(x / params.psi50[ivt, level])^params.ck[ivt, level])
        if val < 0.005
            val = 0.0
        end
        return val
    else
        error("PLC: must choose valid plc_method")
    end
end

"""
    d1plc(x, ivt, level, plc_method, params)

Return 1st derivative of vulnerability curve at water potential `x`.
"""
function d1plc(x::Float64, ivt::Int, level::Int, plc_method::Int, params::PhotoParamsData=params_inst)
    if plc_method == VEGETATION_WEIBULL
        ck_val = params.ck[ivt, level]
        psi50_val = params.psi50[ivt, level]
        val = -ck_val * log(2.0) * (2.0^(-(x / psi50_val)^ck_val)) *
              ((x / psi50_val)^ck_val) / x
        return val
    else
        error("D1PLC: must choose valid plc_method")
    end
end

# =====================================================================
# ci_func! — evaluate f(ci) for the standard (non-PHS) method
# =====================================================================

"""
    ci_func!(ci, p, iv, forc_pbot_c, gb_mol, je, cair, oair, lmr_z, par_z,
             rh_can, photosyns)

Evaluate f(ci) = ci - (ca - (1.37rb+1.65rs))*patm*an.
Returns `(fval, gs_mol)`.
"""
function ci_func!(ci::Float64, p::Int, iv::Int, forc_pbot_c::Float64,
                  gb_mol::Float64, je::Float64, cair::Float64, oair::Float64,
                  lmr_z::Float64, par_z::Float64, rh_can::Float64,
                  ps::PhotosynthesisData;
                  c3psn_val::Float64=1.0,
                  medlynslope_val::Float64=6.0,
                  medlynintercept_val::Float64=100.0,
                  mbbopt_val::Float64=9.0)
    c3flag = ps.c3flag_patch[p]
    ac = ps.ac_patch
    aj = ps.aj_patch
    ap = ps.ap_patch
    ag = ps.ag_patch
    an = ps.an_patch
    vcmax_z = ps.vcmax_z_patch
    cp_p = ps.cp_patch[p]
    kc_p = ps.kc_patch[p]
    ko_p = ps.ko_patch[p]
    qe_p = ps.qe_patch[p]
    tpu_z = ps.tpu_z_patch
    kp_z = ps.kp_z_patch
    bbb_p = ps.bbb_patch[p]
    mbb_p = ps.mbb_patch[p]
    stomatalcond_mtd = ps.stomatalcond_mtd
    theta_cj_val = params_inst.theta_cj[1]  # default; overridden by ivt in real usage

    gs_mol = 0.0

    if c3flag
        ac[p, iv] = vcmax_z[p, iv] * max(ci - cp_p, 0.0) / (ci + kc_p * (1.0 + oair / ko_p))
        aj[p, iv] = je * max(ci - cp_p, 0.0) / (4.0 * ci + 8.0 * cp_p)
        ap[p, iv] = 3.0 * tpu_z[p, iv]
    else
        ac[p, iv] = vcmax_z[p, iv]
        aj[p, iv] = qe_p * par_z * 4.6
        ap[p, iv] = kp_z[p, iv] * max(ci, 0.0) / forc_pbot_c
    end

    # Gross photosynthesis: co-limit ac and aj, then co-limit ap
    aquad = theta_cj_val
    bquad = -(ac[p, iv] + aj[p, iv])
    cquad = ac[p, iv] * aj[p, iv]
    r1, r2 = quadratic_solve(aquad, bquad, cquad)
    ai = min(r1, r2)

    aquad = params_inst.theta_ip
    bquad = -(ai + ap[p, iv])
    cquad = ai * ap[p, iv]
    r1, r2 = quadratic_solve(aquad, bquad, cquad)
    ag[p, iv] = max(0.0, min(r1, r2))

    an[p, iv] = ag[p, iv] - lmr_z
    if an[p, iv] < 0.0
        return (0.0, gs_mol)
    end

    cs = cair - 1.4 / gb_mol * an[p, iv] * forc_pbot_c
    cs = max(cs, MAX_CS)

    if stomatalcond_mtd == STOMATALCOND_MTD_MEDLYN2011
        term = 1.6 * an[p, iv] / (cs / forc_pbot_c * 1.0e06)
        aquad = 1.0
        bquad = -(2.0 * (medlynintercept_val * 1.0e-06 + term) +
                  (medlynslope_val * term)^2 / (gb_mol * 1.0e-06 * rh_can))
        cquad = medlynintercept_val^2 * 1.0e-12 +
                (2.0 * medlynintercept_val * 1.0e-06 + term *
                 (1.0 - medlynslope_val^2 / rh_can)) * term
        r1, r2 = quadratic_solve(aquad, bquad, cquad)
        gs_mol = max(r1, r2) * 1.0e06
    elseif stomatalcond_mtd == STOMATALCOND_MTD_BB1987
        aquad = cs
        bquad = cs * (gb_mol - bbb_p) - mbb_p * an[p, iv] * forc_pbot_c
        cquad = -gb_mol * (cs * bbb_p + mbb_p * an[p, iv] * forc_pbot_c * rh_can)
        r1, r2 = quadratic_solve(aquad, bquad, cquad)
        gs_mol = max(r1, r2)
    end

    fval = ci - cair + an[p, iv] * forc_pbot_c * (1.4 * gs_mol + 1.6 * gb_mol) / (gb_mol * gs_mol)

    return (fval, gs_mol)
end

# =====================================================================
# hybrid_solver! — hybrid Newton-secant / Brent solver for ci
# =====================================================================

"""
    hybrid_solver!(x0, p, iv, forc_pbot_c, gb_mol, je, cair, oair,
                   lmr_z, par_z, rh_can, ps; kwargs...)

Hybrid solver for ci. Returns `(ci_solution, gs_mol, niter)`.
"""
function hybrid_solver!(x0::Float64, p::Int, iv::Int, forc_pbot_c::Float64,
                        gb_mol::Float64, je::Float64, cair::Float64, oair::Float64,
                        lmr_z::Float64, par_z::Float64, rh_can::Float64,
                        ps::PhotosynthesisData; kwargs...)
    eps_val = 1.0e-2
    eps1 = 1.0e-4
    itmax = 40

    f0, gs_mol = ci_func!(x0, p, iv, forc_pbot_c, gb_mol, je, cair, oair,
                           lmr_z, par_z, rh_can, ps; kwargs...)

    if f0 == 0.0
        return (x0, gs_mol, 0)
    end

    minx = x0
    minf = abs(f0)
    x1 = x0 * 0.99

    f1, gs_mol = ci_func!(x1, p, iv, forc_pbot_c, gb_mol, je, cair, oair,
                           lmr_z, par_z, rh_can, ps; kwargs...)

    if f1 == 0.0
        return (x1, gs_mol, 0)
    end
    if abs(f1) < minf
        minx = x1
        minf = abs(f1)
    end

    iter = 0
    while true
        iter += 1
        if (f1 - f0) == 0.0
            break
        end
        dx = -f1 * (x1 - x0) / (f1 - f0)
        x = x1 + dx
        tol = abs(x) * eps_val
        if abs(dx) < tol
            x0 = x
            break
        end
        x0 = x1
        f0 = f1
        x1 = x

        f1, gs_mol = ci_func!(x1, p, iv, forc_pbot_c, gb_mol, je, cair, oair,
                               lmr_z, par_z, rh_can, ps; kwargs...)

        if abs(f1) < minf
            minx = x1
            minf = abs(f1)
        end
        if abs(f1) <= eps1
            x0 = x1
            break
        end

        if f1 * f0 < 0.0
            x_brent = brent_solver!(x0, x1, f0, f1, tol, p, iv, forc_pbot_c,
                                     gb_mol, je, cair, oair, lmr_z, par_z,
                                     rh_can, ps; kwargs...)
            x0 = x_brent
            # get final gs_mol
            _, gs_mol = ci_func!(x0, p, iv, forc_pbot_c, gb_mol, je, cair, oair,
                                  lmr_z, par_z, rh_can, ps; kwargs...)
            break
        end

        if iter > itmax
            _, gs_mol = ci_func!(minx, p, iv, forc_pbot_c, gb_mol, je, cair, oair,
                                  lmr_z, par_z, rh_can, ps; kwargs...)
            x0 = minx
            break
        end
    end

    return (x0, gs_mol, iter)
end

# =====================================================================
# brent_solver! — Brent's method for root finding
# =====================================================================

"""
    brent_solver!(x1, x2, f1, f2, tol, p, iv, forc_pbot_c, gb_mol, je,
                  cair, oair, lmr_z, par_z, rh_can, ps; kwargs...)

Brent's method to find the root of ci_func between x1 and x2.
Returns the root `x`.
"""
function brent_solver!(x1::Float64, x2::Float64, f1::Float64, f2::Float64,
                       tol::Float64, p::Int, iv::Int, forc_pbot_c::Float64,
                       gb_mol::Float64, je::Float64, cair::Float64, oair::Float64,
                       lmr_z::Float64, par_z::Float64, rh_can::Float64,
                       ps::PhotosynthesisData; kwargs...)
    itmax = 20
    eps_val = 1.0e-2

    a = x1
    b = x2
    fa = f1
    fb = f2

    c = b
    fc = fb
    d = b - a
    e = d

    for iter in 1:itmax
        if (fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)
            c = a
            fc = fa
            d = b - a
            e = d
        end
        if abs(fc) < abs(fb)
            a = b
            b = c
            c = a
            fa = fb
            fb = fc
            fc = fa
        end
        tol1 = 2.0 * eps_val * abs(b) + 0.5 * tol
        xm = 0.5 * (c - b)

        if abs(xm) <= tol1 || fb == 0.0
            return b
        end

        if abs(e) >= tol1 && abs(fa) > abs(fb)
            s = fb / fa
            if a == c
                p_val = 2.0 * xm * s
                q_val = 1.0 - s
            else
                q_val = fa / fc
                r_val = fb / fc
                p_val = s * (2.0 * xm * q_val * (q_val - r_val) - (b - a) * (r_val - 1.0))
                q_val = (q_val - 1.0) * (r_val - 1.0) * (s - 1.0)
            end
            if p_val > 0.0
                q_val = -q_val
            end
            p_val = abs(p_val)
            if 2.0 * p_val < min(3.0 * xm * q_val - abs(tol1 * q_val), abs(e * q_val))
                e = d
                d = p_val / q_val
            else
                d = xm
                e = d
            end
        else
            d = xm
            e = d
        end
        a = b
        fa = fb
        if abs(d) > tol1
            b = b + d
        else
            b = b + copysign(tol1, xm)
        end

        fb, _ = ci_func!(b, p, iv, forc_pbot_c, gb_mol, je, cair, oair,
                          lmr_z, par_z, rh_can, ps; kwargs...)

        if fb == 0.0
            return b
        end
    end

    return b
end

# =====================================================================
# photosynthesis! — Main leaf photosynthesis and stomatal conductance
# =====================================================================

"""
    photosynthesis!(ps::PhotosynthesisData, ...)

Leaf photosynthesis and stomatal conductance calculation.
Ported from `subroutine Photosynthesis` in `PhotosynthesisMod.F90`.

This is a simplified interface for unit testing. In production CLM,
additional data types (atm2lnd, temperature, surfalb, solarabs, canopystate,
ozone) would be passed. Here we pass the needed arrays directly.
"""
function photosynthesis!(ps::PhotosynthesisData,
                         esat_tv::Vector{Float64},
                         eair::Vector{Float64},
                         oair::Vector{Float64},
                         cair::Vector{Float64},
                         rb::Vector{Float64},
                         btran::Vector{Float64},
                         dayl_factor::Vector{Float64},
                         leafn::Vector{Float64},
                         forc_pbot::Vector{Float64},
                         t_veg::Vector{Float64},
                         t10::Vector{Float64},
                         tgcm::Vector{Float64},
                         nrad::Vector{Int},
                         tlai_z::Matrix{Float64},
                         tlai::Vector{Float64},
                         par_z_in::Matrix{Float64},
                         lai_z_in::Matrix{Float64},
                         vcmaxcint::Vector{Float64},
                         o3coefv::Vector{Float64},
                         o3coefg::Vector{Float64},
                         c3psn_pft::Vector{Float64},
                         leafcn_pft::Vector{Float64},
                         flnr_pft::Vector{Float64},
                         fnitr_pft::Vector{Float64},
                         slatop_pft::Vector{Float64},
                         mbbopt_pft::Vector{Float64},
                         medlynintercept_pft::Vector{Float64},
                         medlynslope_pft::Vector{Float64},
                         ivt::Vector{Int},
                         col_of_patch::Vector{Int},
                         mask_patch::BitVector,
                         bounds_patch::UnitRange{Int},
                         phase::String;
                         nlevcan::Int=NLEVCAN,
                         use_cn::Bool=false,
                         use_luna::Bool=false,
                         use_c13::Bool=false,
                         leaf_mr_vcm::Float64=0.015,
                         crop_pft::Vector{Float64}=Float64[])
    # Aliases for output arrays based on phase
    if phase == "sun"
        ci_z = ps.cisun_z_patch
        rs = ps.rssun_patch
        rs_z = ps.rssun_z_patch
        lmr_out = ps.lmrsun_patch
        lmr_z = ps.lmrsun_z_patch
        psn = ps.psnsun_patch
        psn_z = ps.psnsun_z_patch
        psn_wc = ps.psnsun_wc_patch
        psn_wj = ps.psnsun_wj_patch
        psn_wp = ps.psnsun_wp_patch
    else
        ci_z = ps.cisha_z_patch
        rs = ps.rssha_patch
        rs_z = ps.rssha_z_patch
        lmr_out = ps.lmrsha_patch
        lmr_z = ps.lmrsha_z_patch
        psn = ps.psnsha_patch
        psn_z = ps.psnsha_z_patch
        psn_wc = ps.psnsha_wc_patch
        psn_wj = ps.psnsha_wj_patch
        psn_wp = ps.psnsha_wp_patch
    end

    ac = ps.ac_patch
    aj = ps.aj_patch
    ap = ps.ap_patch
    ag = ps.ag_patch
    an = ps.an_patch
    vcmax_z = ps.vcmax_z_patch
    tpu_z = ps.tpu_z_patch
    kp_z = ps.kp_z_patch
    gs_mol = ps.gs_mol_patch
    c3flag = ps.c3flag_patch
    cp = ps.cp_patch
    kc = ps.kc_patch
    ko = ps.ko_patch
    qe = ps.qe_patch
    bbb = ps.bbb_patch
    mbb = ps.mbb_patch
    gb_mol = ps.gb_mol_patch
    rh_leaf = ps.rh_leaf_patch
    vpd_can = ps.vpd_can_patch
    lnc = ps.lnca_patch
    stomatalcond_mtd = ps.stomatalcond_mtd
    leafresp_method = ps.leafresp_method
    light_inhibit = ps.light_inhibit

    lmrc = fth25_photo(params_inst.lmrhd, params_inst.lmrse)

    np = length(bounds_patch)
    jmax_z_local = zeros(np, nlevcan)
    bbbopt = zeros(np)
    kn = zeros(np)
    psn_wc_z = zeros(np, nlevcan)
    psn_wj_z = zeros(np, nlevcan)
    psn_wp_z = zeros(np, nlevcan)

    rsmax0 = 2.0e4

    # ---- Pass 1: Compute kinetics, N scaling, respiration ----
    for p in bounds_patch
        mask_patch[p] || continue
        c = col_of_patch[p]
        ivt_p = ivt[p]

        # C3 or C4
        if round(Int, c3psn_pft[ivt_p]) == 1
            c3flag[p] = true
        else
            c3flag[p] = false
        end

        # C3/C4 dependent parameters
        if c3flag[p]
            qe[p] = 0.0
            if stomatalcond_mtd == STOMATALCOND_MTD_BB1987
                bbbopt[p] = BBBOPT_C3
            end
        else
            qe[p] = 0.05
            if stomatalcond_mtd == STOMATALCOND_MTD_BB1987
                bbbopt[p] = BBBOPT_C4
            end
        end

        if stomatalcond_mtd == STOMATALCOND_MTD_BB1987
            bbb[p] = max(bbbopt[p] * btran[p], 1.0)
            mbb[p] = mbbopt_pft[ivt_p]
        end

        # Michaelis-Menten constants
        kc25 = params_inst.kc25_coef * forc_pbot[c]
        ko25 = params_inst.ko25_coef * forc_pbot[c]
        sco = 0.5 * 0.209 / params_inst.cp25_yr2000
        cp25 = 0.5 * oair[p] / sco

        kc[p] = kc25 * ft_photo(t_veg[p], params_inst.kcha)
        ko[p] = ko25 * ft_photo(t_veg[p], params_inst.koha)
        cp[p] = cp25 * ft_photo(t_veg[p], params_inst.cpha)
    end

    # ---- Pass 2: Nitrogen profile and vcmax ----
    for p in bounds_patch
        mask_patch[p] || continue
        ivt_p = ivt[p]

        lnc[p] = 1.0 / (slatop_pft[ivt_p] * leafcn_pft[ivt_p])

        vcmax25top = lnc[p] * flnr_pft[ivt_p] * params_inst.fnr * params_inst.act25 * dayl_factor[p]
        if !use_cn
            vcmax25top = vcmax25top * fnitr_pft[ivt_p]
        end

        jmax25top = ((2.59 - 0.035 * min(max((t10[p] - TFRZ), 11.0), 35.0)) * vcmax25top) *
                    params_inst.jmax25top_sf
        tpu25top = params_inst.tpu25ratio * vcmax25top
        kp25top = params_inst.kp25ratio * vcmax25top

        if dayl_factor[p] < 1.0e-12
            kn[p] = 0.0
        else
            kn[p] = exp(0.00963 * vcmax25top / dayl_factor[p] - 2.43)
        end

        if use_cn
            if leafresp_method == LEAFRESP_MTD_RYAN1991
                lmr25top = 2.525e-6 * (1.5^((25.0 - 20.0) / 10.0))
                lmr25top = lmr25top * lnc[p] / 12.0e-06
            else
                lmr25top = 0.0
            end
        else
            if c3flag[p]
                lmr25top = vcmax25top * leaf_mr_vcm
            else
                lmr25top = vcmax25top * 0.025
            end
        end

        # Loop through canopy layers
        laican = 0.0
        for iv in 1:nrad[p]
            if iv == 1
                laican = 0.5 * tlai_z[p, iv]
            else
                laican = laican + 0.5 * (tlai_z[p, iv-1] + tlai_z[p, iv])
            end

            if nlevcan == 1
                nscaler = vcmaxcint[p]
            else
                nscaler = exp(-kn[p] * laican)
            end

            # Maintenance respiration
            lmr25 = lmr25top * nscaler

            if c3flag[p]
                lmr_z[p, iv] = lmr25 * ft_photo(t_veg[p], params_inst.lmrha) *
                               fth_photo(t_veg[p], params_inst.lmrhd, params_inst.lmrse, lmrc)
            else
                lmr_z[p, iv] = lmr25 * 2.0^((t_veg[p] - (TFRZ + 25.0)) / 10.0)
                lmr_z[p, iv] = lmr_z[p, iv] / (1.0 + exp(1.3 * (t_veg[p] - (TFRZ + 55.0))))
            end

            if par_z_in[p, iv] <= 0.0  # night time
                vcmax_z[p, iv] = 0.0
                jmax_z_local[p, iv] = 0.0
                tpu_z[p, iv] = 0.0
                kp_z[p, iv] = 0.0
            else  # day time
                vcmax25 = vcmax25top * nscaler
                jmax25 = jmax25top * nscaler
                tpu25 = tpu25top * nscaler
                kp25 = kp25top * nscaler

                # Temperature adjustment
                vcmaxse = (668.39 - 1.07 * min(max((t10[p] - TFRZ), 11.0), 35.0)) * params_inst.vcmaxse_sf
                jmaxse = (659.70 - 0.75 * min(max((t10[p] - TFRZ), 11.0), 35.0)) * params_inst.jmaxse_sf
                tpuse = (668.39 - 1.07 * min(max((t10[p] - TFRZ), 11.0), 35.0)) * params_inst.tpuse_sf
                vcmaxc = fth25_photo(params_inst.vcmaxhd, vcmaxse)
                jmaxc = fth25_photo(params_inst.jmaxhd, jmaxse)
                tpuc = fth25_photo(params_inst.tpuhd, tpuse)

                vcmax_z[p, iv] = vcmax25 * ft_photo(t_veg[p], params_inst.vcmaxha) *
                                  fth_photo(t_veg[p], params_inst.vcmaxhd, vcmaxse, vcmaxc)
                jmax_z_local[p, iv] = jmax25 * ft_photo(t_veg[p], params_inst.jmaxha) *
                                       fth_photo(t_veg[p], params_inst.jmaxhd, jmaxse, jmaxc)
                tpu_z[p, iv] = tpu25 * ft_photo(t_veg[p], params_inst.tpuha) *
                                fth_photo(t_veg[p], params_inst.tpuhd, tpuse, tpuc)

                if !c3flag[p]
                    vcmax_z[p, iv] = vcmax25 * 2.0^((t_veg[p] - (TFRZ + 25.0)) / 10.0)
                    vcmax_z[p, iv] = vcmax_z[p, iv] / (1.0 + exp(0.2 * ((TFRZ + 15.0) - t_veg[p])))
                    vcmax_z[p, iv] = vcmax_z[p, iv] / (1.0 + exp(0.3 * (t_veg[p] - (TFRZ + 40.0))))
                end

                kp_z[p, iv] = kp25 * 2.0^((t_veg[p] - (TFRZ + 25.0)) / 10.0)
            end

            # Soil water stress
            vcmax_z[p, iv] = vcmax_z[p, iv] * btran[p]
            lmr_z[p, iv] = lmr_z[p, iv] * btran[p]

            # Light inhibition
            if light_inhibit && par_z_in[p, 1] > 0.0
                lmr_z[p, iv] = lmr_z[p, iv] * 0.67
            end
        end
    end

    # ---- Pass 3: Leaf-level photosynthesis and stomatal conductance ----
    for p in bounds_patch
        mask_patch[p] || continue
        c = col_of_patch[p]
        ivt_p = ivt[p]

        cf = forc_pbot[c] / (RGAS * 1.0e-3 * tgcm[p]) * 1.0e06
        gb = 1.0 / rb[p]
        gb_mol[p] = gb * cf

        for iv in 1:nrad[p]
            if par_z_in[p, iv] <= 0.0  # night time
                ac[p, iv] = 0.0
                aj[p, iv] = 0.0
                ap[p, iv] = 0.0
                ag[p, iv] = 0.0
                an[p, iv] = ag[p, iv] - lmr_z[p, iv]
                psn_z[p, iv] = 0.0
                psn_wc_z[p, iv] = 0.0
                psn_wj_z[p, iv] = 0.0
                psn_wp_z[p, iv] = 0.0
                if stomatalcond_mtd == STOMATALCOND_MTD_BB1987
                    rs_z[p, iv] = min(rsmax0, 1.0 / bbb[p] * cf)
                elseif stomatalcond_mtd == STOMATALCOND_MTD_MEDLYN2011
                    rs_z[p, iv] = min(rsmax0, 1.0 / medlynintercept_pft[ivt_p] * cf)
                end
                ci_z[p, iv] = 0.0
                rh_leaf[p] = 0.0
                gs_mol[p, iv] = stomatalcond_mtd == STOMATALCOND_MTD_BB1987 ? bbb[p] : medlynintercept_pft[ivt_p]
            else  # day time
                ceair = min(eair[p], esat_tv[p])
                if stomatalcond_mtd == STOMATALCOND_MTD_BB1987
                    rh_can = ceair / esat_tv[p]
                else
                    rh_can = max((esat_tv[p] - ceair), MEDLYN_RH_CAN_MAX) * MEDLYN_RH_CAN_FACT
                    vpd_can[p] = rh_can
                end

                # Electron transport rate
                qabs = 0.5 * (1.0 - params_inst.fnps) * par_z_in[p, iv] * 4.6
                aquad = params_inst.theta_psii
                bquad = -(qabs + jmax_z_local[p, iv])
                cquad = qabs * jmax_z_local[p, iv]
                r1, r2 = quadratic_solve(aquad, bquad, cquad)
                je = min(r1, r2)

                # Initial guess for ci
                if c3flag[p]
                    ci_z[p, iv] = 0.7 * cair[p]
                else
                    ci_z[p, iv] = 0.4 * cair[p]
                end

                # Use theta_cj for this PFT
                theta_cj_val = params_inst.theta_cj[ivt_p]

                # Solve for ci and gs
                ci_sol, gs_mol_val, _ = hybrid_solver!(ci_z[p, iv], p, iv,
                    forc_pbot[c], gb_mol[p], je, cair[p], oair[p],
                    lmr_z[p, iv], par_z_in[p, iv], rh_can, ps;
                    medlynslope_val=medlynslope_pft[ivt_p],
                    medlynintercept_val=medlynintercept_pft[ivt_p],
                    mbbopt_val=mbbopt_pft[ivt_p])

                gs_mol[p, iv] = gs_mol_val

                # Check for an < 0
                if an[p, iv] < 0.0
                    if stomatalcond_mtd == STOMATALCOND_MTD_BB1987
                        gs_mol[p, iv] = bbb[p]
                    else
                        gs_mol[p, iv] = medlynintercept_pft[ivt_p]
                    end
                end

                # Store sunlit/shaded gs
                if phase == "sun"
                    ps.gs_mol_sun_patch[p, iv] = gs_mol[p, iv]
                else
                    ps.gs_mol_sha_patch[p, iv] = gs_mol[p, iv]
                end

                # Final ci
                cs = cair[p] - 1.4 / gb_mol[p] * an[p, iv] * forc_pbot[c]
                cs = max(cs, MAX_CS)
                ci_z[p, iv] = cair[p] - an[p, iv] * forc_pbot[c] *
                    (1.4 * gs_mol[p, iv] + 1.6 * gb_mol[p]) / (gb_mol[p] * gs_mol[p, iv])
                ci_z[p, iv] = max(ci_z[p, iv], 1.0e-06)

                # Convert to resistance
                gs_val = gs_mol[p, iv] / cf
                rs_z[p, iv] = min(1.0 / gs_val, rsmax0)
                rs_z[p, iv] = rs_z[p, iv] / o3coefg[p]

                # Photosynthesis
                psn_z[p, iv] = ag[p, iv]
                psn_z[p, iv] = psn_z[p, iv] * o3coefv[p]

                psn_wc_z[p, iv] = 0.0
                psn_wj_z[p, iv] = 0.0
                psn_wp_z[p, iv] = 0.0

                if ac[p, iv] <= aj[p, iv] && ac[p, iv] <= ap[p, iv]
                    psn_wc_z[p, iv] = psn_z[p, iv]
                elseif aj[p, iv] < ac[p, iv] && aj[p, iv] <= ap[p, iv]
                    psn_wj_z[p, iv] = psn_z[p, iv]
                elseif ap[p, iv] < ac[p, iv] && ap[p, iv] < aj[p, iv]
                    psn_wp_z[p, iv] = psn_z[p, iv]
                end

                # Ball-Berry error check
                if stomatalcond_mtd == STOMATALCOND_MTD_BB1987
                    hs = (gb_mol[p] * ceair + gs_mol[p, iv] * esat_tv[p]) /
                         ((gb_mol[p] + gs_mol[p, iv]) * esat_tv[p])
                    rh_leaf[p] = hs
                end
            end
        end
    end

    # ---- Pass 4: Canopy integration ----
    for p in bounds_patch
        mask_patch[p] || continue

        psncan = 0.0
        psncan_wc = 0.0
        psncan_wj = 0.0
        psncan_wp = 0.0
        lmrcan = 0.0
        gscan = 0.0
        laican = 0.0

        for iv in 1:nrad[p]
            psncan += psn_z[p, iv] * lai_z_in[p, iv]
            psncan_wc += psn_wc_z[p, iv] * lai_z_in[p, iv]
            psncan_wj += psn_wj_z[p, iv] * lai_z_in[p, iv]
            psncan_wp += psn_wp_z[p, iv] * lai_z_in[p, iv]
            lmrcan += lmr_z[p, iv] * lai_z_in[p, iv]
            gscan += lai_z_in[p, iv] / (rb[p] + rs_z[p, iv])
            laican += lai_z_in[p, iv]
        end

        if laican > 0.0
            psn[p] = psncan / laican
            psn_wc[p] = psncan_wc / laican
            psn_wj[p] = psncan_wj / laican
            psn_wp[p] = psncan_wp / laican
            lmr_out[p] = lmrcan / laican
            rs[p] = laican / gscan - rb[p]
        else
            psn[p] = 0.0
            psn_wc[p] = 0.0
            psn_wj[p] = 0.0
            psn_wp[p] = 0.0
            lmr_out[p] = 0.0
            rs[p] = 0.0
        end
    end

    return nothing
end

# =====================================================================
# photosynthesis_total! — Determine total photosynthesis
# =====================================================================

"""
    photosynthesis_total!(ps::PhotosynthesisData, laisun, laisha, mask_patch, bounds_patch;
                          use_fates=false, use_cn=false, use_c13=false, use_c14=false)

Determine total photosynthesis by combining sunlit and shaded fluxes.
Ported from `subroutine PhotosynthesisTotal`.
"""
function photosynthesis_total!(ps::PhotosynthesisData,
                               laisun::Vector{Float64},
                               laisha::Vector{Float64},
                               mask_patch::BitVector,
                               bounds_patch::UnitRange{Int};
                               use_fates::Bool=false,
                               use_cn::Bool=false,
                               use_c13::Bool=false,
                               use_c14::Bool=false)
    for p in bounds_patch
        mask_patch[p] || continue

        if !use_fates
            ps.fpsn_patch[p] = ps.psnsun_patch[p] * laisun[p] + ps.psnsha_patch[p] * laisha[p]
            ps.fpsn_wc_patch[p] = ps.psnsun_wc_patch[p] * laisun[p] + ps.psnsha_wc_patch[p] * laisha[p]
            ps.fpsn_wj_patch[p] = ps.psnsun_wj_patch[p] * laisun[p] + ps.psnsha_wj_patch[p] * laisha[p]
            ps.fpsn_wp_patch[p] = ps.psnsun_wp_patch[p] * laisun[p] + ps.psnsha_wp_patch[p] * laisha[p]
        end
    end

    return nothing
end

# =====================================================================
# fractionation! — C13 fractionation during photosynthesis
# =====================================================================

"""
    fractionation!(ps::PhotosynthesisData, ...)

C13 fractionation during photosynthesis.
Ported from `subroutine Fractionation`.
"""
function fractionation!(ps::PhotosynthesisData,
                        forc_pbot::Vector{Float64},
                        forc_pco2::Vector{Float64},
                        par_z_in::Matrix{Float64},
                        nrad::Vector{Int},
                        c3psn_pft::Vector{Float64},
                        ivt::Vector{Int},
                        col_of_patch::Vector{Int},
                        gridcell_of_patch::Vector{Int},
                        mask_patch::BitVector,
                        bounds_patch::UnitRange{Int},
                        phase::String;
                        use_hydrstress::Bool=false)
    if phase == "sun"
        alphapsn = ps.alphapsnsun_patch
        if use_hydrstress
            gs_mol_ref = ps.gs_mol_sun_patch
            an_ref = ps.an_sun_patch
        else
            gs_mol_ref = ps.gs_mol_patch
            an_ref = ps.an_patch
        end
    else
        alphapsn = ps.alphapsnsha_patch
        if use_hydrstress
            gs_mol_ref = ps.gs_mol_sha_patch
            an_ref = ps.an_sha_patch
        else
            gs_mol_ref = ps.gs_mol_patch
            an_ref = ps.an_patch
        end
    end

    gb_mol = ps.gb_mol_patch

    for p in bounds_patch
        mask_patch[p] || continue
        c = col_of_patch[p]
        g = gridcell_of_patch[p]
        co2_p = forc_pco2[g]

        for iv in 1:nrad[p]
            if par_z_in[p, iv] <= 0.0
                alphapsn[p] = 1.0
            else
                ci = co2_p - (an_ref[p, iv] *
                    forc_pbot[c] *
                    (1.4 * gs_mol_ref[p, iv] + 1.6 * gb_mol[p]) /
                    (gb_mol[p] * gs_mol_ref[p, iv]))
                alphapsn[p] = 1.0 + (((c3psn_pft[ivt[p]] *
                    (4.4 + (22.6 * (ci / co2_p)))) +
                    ((1.0 - c3psn_pft[ivt[p]]) * 4.4)) / 1000.0)
            end
        end
    end

    return nothing
end

# =====================================================================
# TimeStepInit — Initialize at start of time step
# =====================================================================

"""
    photosynthesis_timestep_init!(ps::PhotosynthesisData, mask_nolake, bounds_patch;
                                  use_c13=false, use_c14=false)

Time step initialization. Zeros out photosynthesis fluxes for non-lake patches.
"""
function photosynthesis_timestep_init!(ps::PhotosynthesisData,
                                       mask_nolake::BitVector,
                                       bounds_patch::UnitRange{Int};
                                       use_c13::Bool=false,
                                       use_c14::Bool=false)
    for p in bounds_patch
        mask_nolake[p] || continue

        ps.psnsun_patch[p] = 0.0
        ps.psnsun_wc_patch[p] = 0.0
        ps.psnsun_wj_patch[p] = 0.0
        ps.psnsun_wp_patch[p] = 0.0

        ps.psnsha_patch[p] = 0.0
        ps.psnsha_wc_patch[p] = 0.0
        ps.psnsha_wj_patch[p] = 0.0
        ps.psnsha_wp_patch[p] = 0.0

        ps.fpsn_patch[p] = 0.0
        ps.fpsn_wc_patch[p] = 0.0
        ps.fpsn_wj_patch[p] = 0.0
        ps.fpsn_wp_patch[p] = 0.0

        if use_c13
            ps.alphapsnsun_patch[p] = 0.0
            ps.alphapsnsha_patch[p] = 0.0
            ps.c13_psnsun_patch[p] = 0.0
            ps.c13_psnsha_patch[p] = 0.0
        end
        if use_c14
            ps.c14_psnsun_patch[p] = 0.0
            ps.c14_psnsha_patch[p] = 0.0
        end
    end

    return nothing
end

# =====================================================================
# SetParamsForTesting — Set parameters for unit testing
# =====================================================================

"""
    set_params_for_testing!(ps::PhotosynthesisData)

Set parameters for unit testing (from `subroutine setParamsForTesting`).
"""
function set_params_for_testing!(ps::PhotosynthesisData)
    photo_params_init!(params_inst)
    params_inst.ck .= 3.95
    params_inst.psi50[1, :] .= -150000.0
    params_inst.psi50[2, :] .= -530000.0
    params_inst.psi50[3:12, :] .= -400000.0
    if size(params_inst.psi50, 1) >= 13
        params_inst.psi50[13:end, :] .= -340000.0
    end
    return nothing
end

# =====================================================================
# PHS (Plant Hydraulic Stress) Functions
# =====================================================================
# These implement the simultaneous sunlit/shaded solution for
# photosynthesis with plant hydraulic stress.
# Ported from subroutines in PhotosynthesisMod.F90.
# =====================================================================

# =====================================================================
# getqflx! — Calculate sunlit and shaded transpiration
# =====================================================================

"""
    getqflx!(p, c, gb_mol, gs_mol_sun, gs_mol_sha, qflx_sun, qflx_sha,
             qsatl, qaf, havegs, laisun, laisha, elai, esai,
             fdry, forc_rho, forc_pbot, tgcm)

Calculate sunlit and shaded transpiration using gb_mol and gs_mol.
When `havegs=true`, compute qflx from gs. When `havegs=false`, compute gs from qflx.
Returns `(qflx_sun, qflx_sha, gs_mol_sun, gs_mol_sha)`.
Ported from `subroutine getqflx`.
"""
function getqflx!(p::Int, c::Int, gb_mol::Float64,
                  gs_mol_sun::Float64, gs_mol_sha::Float64,
                  qflx_sun::Float64, qflx_sha::Float64,
                  qsatl::Float64, qaf::Float64, havegs::Bool,
                  laisun_p::Float64, laisha_p::Float64,
                  elai_p::Float64, esai_p::Float64,
                  fdry_p::Float64, forc_rho_c::Float64,
                  forc_pbot_c::Float64, tgcm_p::Float64)

    cf = forc_pbot_c / (RGAS * 1.0e-3 * tgcm_p) * 1.0e6
    wtl = (elai_p + esai_p) * gb_mol
    efpot = forc_rho_c * wtl * (qsatl - qaf)

    if havegs
        if efpot > 0.0 && elai_p > 0.0
            if gs_mol_sun > 0.0
                rppdry_sun = fdry_p / gb_mol * (laisun_p / (1.0 / gb_mol + 1.0 / gs_mol_sun)) / elai_p
                qflx_sun = efpot * rppdry_sun / cf
            else
                qflx_sun = 0.0
            end
            if gs_mol_sha > 0.0
                rppdry_sha = fdry_p / gb_mol * (laisha_p / (1.0 / gb_mol + 1.0 / gs_mol_sha)) / elai_p
                qflx_sha = efpot * rppdry_sha / cf
            else
                qflx_sha = 0.0
            end
        else
            qflx_sun = 0.0
            qflx_sha = 0.0
        end
    else
        if qflx_sun > 0.0
            gs_mol_sun = gb_mol * qflx_sun * cf * elai_p / (efpot * fdry_p * laisun_p - qflx_sun * cf * elai_p)
        else
            gs_mol_sun = 0.0
        end
        if qflx_sha > 0.0
            gs_mol_sha = gb_mol * qflx_sha * cf * elai_p / (efpot * fdry_p * laisha_p - qflx_sha * cf * elai_p)
        else
            gs_mol_sha = 0.0
        end
    end

    return (qflx_sun, qflx_sha, gs_mol_sun, gs_mol_sha)
end

# =====================================================================
# spacF! — Flux divergence across vegetation segments
# =====================================================================

"""
    spacF!(p, c, x, qflx_sun, qflx_sha, k_soil_root_p, smp_c, z_c,
           laisun_p, laisha_p, htop_p, tsai_p, ivt_p, nlevsoi)

Returns `f`, the flux divergence across each vegetation segment
calculated for vegwp as passed in via `x`.
Ported from `subroutine spacF`.
"""
function spacF!(p::Int, c::Int, x::Vector{Float64},
                qflx_sun::Float64, qflx_sha::Float64,
                k_soil_root_p::AbstractVector{Float64},
                smp_c::AbstractVector{Float64},
                z_c::AbstractVector{Float64},
                laisun_p::Float64, laisha_p::Float64,
                htop_p::Float64, tsai_p::Float64,
                ivt_p::Int, nlevsoi::Int)
    tol_lai = 0.001

    grav1 = htop_p * 1000.0
    grav2 = z_c[1:nlevsoi] .* 1000.0

    fsto1 = plc(x[SUN], ivt_p, SUN, VEG)
    fsto2 = plc(x[SHA], ivt_p, SHA, VEG)
    fx    = plc(x[XYL], ivt_p, XYL, VEG)
    fr    = plc(x[ROOT_SEG], ivt_p, ROOT_SEG, VEG)

    f = zeros(NVEGWCS)
    f[SUN] = qflx_sun * fsto1 - laisun_p * params_inst.kmax[ivt_p, SUN] * fx * (x[XYL] - x[SUN])
    f[SHA] = qflx_sha * fsto2 - laisha_p * params_inst.kmax[ivt_p, SHA] * fx * (x[XYL] - x[SHA])
    f[XYL] = laisun_p * params_inst.kmax[ivt_p, SUN] * fx * (x[XYL] - x[SUN]) +
             laisha_p * params_inst.kmax[ivt_p, SHA] * fx * (x[XYL] - x[SHA]) -
             tsai_p * params_inst.kmax[ivt_p, XYL] / htop_p * fr * (x[ROOT_SEG] - x[XYL] - grav1)
    f[ROOT_SEG] = tsai_p * params_inst.kmax[ivt_p, XYL] / htop_p * fr * (x[ROOT_SEG] - x[XYL] - grav1) +
                  sum(k_soil_root_p[1:nlevsoi] .* (x[ROOT_SEG] .+ grav2)) -
                  sum(k_soil_root_p[1:nlevsoi] .* smp_c[1:nlevsoi])

    if laisha_p < tol_lai
        temp = f[SUN]
        f[SUN] = f[SHA]
        f[SHA] = temp
    end

    return f
end

# =====================================================================
# spacA! — Inverse Jacobian for vegetation water potential
# =====================================================================

"""
    spacA!(p, c, x, qflx_sun, qflx_sha, k_soil_root_p,
           laisun_p, laisha_p, htop_p, tsai_p, ivt_p, nlevsoi)

Returns `(invA, flag)` where `invA` is the inverse matrix relating
delta(vegwp) to f: d(vegwp)=invA*f. `flag` is true if matrix is singular.
Ported from `subroutine spacA`.
"""
function spacA!(p::Int, c::Int, x::Vector{Float64},
                qflx_sun::Float64, qflx_sha::Float64,
                k_soil_root_p::AbstractVector{Float64},
                laisun_p::Float64, laisha_p::Float64,
                htop_p::Float64, tsai_p::Float64,
                ivt_p::Int, nlevsoi::Int)
    tol_lai = 0.001

    A = zeros(NVEGWCS, NVEGWCS)
    invA = zeros(NVEGWCS, NVEGWCS)

    grav1 = htop_p * 1000.0

    fsto1 = plc(x[SUN], ivt_p, SUN, VEG)
    fsto2 = plc(x[SHA], ivt_p, SHA, VEG)
    fx    = plc(x[XYL], ivt_p, XYL, VEG)
    fr    = plc(x[ROOT_SEG], ivt_p, ROOT_SEG, VEG)

    dfsto1 = d1plc(x[SUN], ivt_p, SUN, VEG)
    dfsto2 = d1plc(x[SHA], ivt_p, SHA, VEG)
    dfx    = d1plc(x[XYL], ivt_p, XYL, VEG)
    dfr    = d1plc(x[ROOT_SEG], ivt_p, ROOT_SEG, VEG)

    # Build A matrix: f = A * d(vegwp)
    A[1,1] = -laisun_p * params_inst.kmax[ivt_p, SUN] * fx - qflx_sun * dfsto1
    A[1,3] = laisun_p * params_inst.kmax[ivt_p, SUN] * dfx * (x[XYL] - x[SUN]) +
             laisun_p * params_inst.kmax[ivt_p, SUN] * fx
    A[2,2] = -laisha_p * params_inst.kmax[ivt_p, SHA] * fx - qflx_sha * dfsto2
    A[2,3] = laisha_p * params_inst.kmax[ivt_p, SHA] * dfx * (x[XYL] - x[SHA]) +
             laisha_p * params_inst.kmax[ivt_p, SHA] * fx
    A[3,1] = laisun_p * params_inst.kmax[ivt_p, SUN] * fx
    A[3,2] = laisha_p * params_inst.kmax[ivt_p, SHA] * fx
    A[3,3] = -laisun_p * params_inst.kmax[ivt_p, SUN] * dfx * (x[XYL] - x[SUN]) -
              laisun_p * params_inst.kmax[ivt_p, SUN] * fx -
              laisha_p * params_inst.kmax[ivt_p, SHA] * dfx * (x[XYL] - x[SHA]) -
              laisha_p * params_inst.kmax[ivt_p, SHA] * fx -
              tsai_p * params_inst.kmax[ivt_p, XYL] / htop_p * fr
    A[3,4] = tsai_p * params_inst.kmax[ivt_p, XYL] / htop_p * dfr * (x[ROOT_SEG] - x[XYL] - grav1) +
             tsai_p * params_inst.kmax[ivt_p, XYL] / htop_p * fr
    A[4,3] = tsai_p * params_inst.kmax[ivt_p, XYL] / htop_p * fr
    A[4,4] = -tsai_p * params_inst.kmax[ivt_p, XYL] / htop_p * fr -
              tsai_p * params_inst.kmax[ivt_p, XYL] / htop_p * dfr * (x[ROOT_SEG] - x[XYL] - grav1) -
              sum(k_soil_root_p[1:nlevsoi])

    # Matrix inversion
    if laisun_p > tol_lai && laisha_p > tol_lai
        # General 4x4 case
        determ = A[4,4]*A[2,2]*A[3,3]*A[1,1] - A[4,4]*A[2,2]*A[3,1]*A[1,3] -
                 A[4,4]*A[3,2]*A[2,3]*A[1,1] - A[4,3]*A[1,1]*A[2,2]*A[3,4]
        if abs(determ) <= 1.0e-50
            return (invA, true)
        end
        leading = 1.0 / determ

        invA[1,1] = leading*A[4,4]*A[2,2]*A[3,3] - leading*A[4,4]*A[3,2]*A[2,3] - leading*A[4,3]*A[2,2]*A[3,4]
        invA[2,1] = leading*A[2,3]*A[4,4]*A[3,1]
        invA[3,1] = -leading*A[4,4]*A[2,2]*A[3,1]
        invA[4,1] = leading*A[4,3]*A[2,2]*A[3,1]
        invA[1,2] = leading*A[1,3]*A[4,4]*A[3,2]
        invA[2,2] = leading*A[4,4]*A[3,3]*A[1,1] - leading*A[4,4]*A[3,1]*A[1,3] - leading*A[4,3]*A[1,1]*A[3,4]
        invA[3,2] = -leading*A[1,1]*A[4,4]*A[3,2]
        invA[4,2] = leading*A[4,3]*A[1,1]*A[3,2]
        invA[1,3] = -leading*A[1,3]*A[2,2]*A[4,4]
        invA[2,3] = -leading*A[2,3]*A[1,1]*A[4,4]
        invA[3,3] = leading*A[2,2]*A[1,1]*A[4,4]
        invA[4,3] = -leading*A[4,3]*A[1,1]*A[2,2]
        invA[1,4] = leading*A[1,3]*A[3,4]*A[2,2]
        invA[2,4] = leading*A[2,3]*A[3,4]*A[1,1]
        invA[3,4] = -leading*A[3,4]*A[1,1]*A[2,2]
        invA[4,4] = leading*A[2,2]*A[3,3]*A[1,1] - leading*A[2,2]*A[3,1]*A[1,3] - leading*A[3,2]*A[2,3]*A[1,1]
    else
        # 3x3 case when one of laisun/laisha is ~0
        if laisha_p <= tol_lai
            A[2,2] = A[1,1]
            A[3,2] = A[3,1]
            A[2,3] = A[1,3]
        end
        determ = A[2,2]*A[3,3]*A[4,4] - A[3,4]*A[2,2]*A[4,3] - A[2,3]*A[3,2]*A[4,4]
        if abs(determ) <= 1.0e-50
            return (invA, true)
        end

        invA[2,2] = A[3,3]*A[4,4] - A[3,4]*A[4,3]
        invA[2,3] = -A[2,3]*A[4,4]
        invA[2,4] = A[3,4]*A[2,3]
        invA[3,2] = -A[3,2]*A[4,4]
        invA[3,3] = A[2,2]*A[4,4]
        invA[3,4] = -A[3,4]*A[2,2]
        invA[4,2] = A[3,2]*A[4,3]
        invA[4,3] = -A[2,2]*A[4,3]
        invA[4,4] = A[2,2]*A[3,3] - A[2,3]*A[3,2]
        invA .= (1.0 / determ) .* invA
    end

    return (invA, false)
end

# =====================================================================
# getvegwp! — Calculate vegetation water potential
# =====================================================================

"""
    getvegwp!(p, c, gb_mol, gs_mol_sun, gs_mol_sha, qsatl, qaf,
              k_soil_root_p, smp_c, z_c, laisun_p, laisha_p,
              htop_p, tsai_p, ivt_p, nlevsoi,
              elai_p, esai_p, fdry_p, forc_rho_c, forc_pbot_c, tgcm_p)

Calculates transpiration and returns corresponding vegwp in x and soilflux.
Returns `(x, soilflux)`.
Ported from `subroutine getvegwp`.
"""
function getvegwp!(p::Int, c::Int, gb_mol::Float64,
                   gs_mol_sun::Float64, gs_mol_sha::Float64,
                   qsatl::Float64, qaf::Float64,
                   k_soil_root_p::AbstractVector{Float64},
                   smp_c::AbstractVector{Float64},
                   z_c::AbstractVector{Float64},
                   laisun_p::Float64, laisha_p::Float64,
                   htop_p::Float64, tsai_p::Float64,
                   ivt_p::Int, nlevsoi::Int,
                   elai_p::Float64, esai_p::Float64,
                   fdry_p::Float64, forc_rho_c::Float64,
                   forc_pbot_c::Float64, tgcm_p::Float64)
    x = zeros(NVEGWCS)

    grav1 = 1000.0 * htop_p
    grav2 = 1000.0 .* z_c[1:nlevsoi]

    # Compute transpiration demand
    havegs = true
    qflx_sun, qflx_sha, _, _ = getqflx!(p, c, gb_mol, gs_mol_sun, gs_mol_sha,
        0.0, 0.0, qsatl, qaf, havegs, laisun_p, laisha_p,
        elai_p, esai_p, fdry_p, forc_rho_c, forc_pbot_c, tgcm_p)

    # Calculate root water potential
    sum_ksr = sum(k_soil_root_p[1:nlevsoi])
    if abs(sum_ksr) == 0.0
        x[ROOT_SEG] = sum(smp_c[1:nlevsoi] .- grav2) / nlevsoi
    else
        x[ROOT_SEG] = (sum(k_soil_root_p[1:nlevsoi] .* (smp_c[1:nlevsoi] .- grav2)) - qflx_sun - qflx_sha) / sum_ksr
    end

    # Calculate xylem water potential
    fr = plc(x[ROOT_SEG], ivt_p, ROOT_SEG, VEG)
    if tsai_p > 0.0 && fr > 0.0
        x[XYL] = x[ROOT_SEG] - grav1 - (qflx_sun + qflx_sha) / (fr * params_inst.kmax[ivt_p, ROOT_SEG] / htop_p * tsai_p)
    else
        x[XYL] = x[ROOT_SEG] - grav1
    end

    # Calculate sun/sha leaf water potential
    fx = plc(x[XYL], ivt_p, XYL, VEG)
    if laisha_p > 0.0 && fx > 0.0
        x[SHA] = x[XYL] - qflx_sha / (fx * params_inst.kmax[ivt_p, XYL] * laisha_p)
    else
        x[SHA] = x[XYL]
    end
    if laisun_p > 0.0 && fx > 0.0
        x[SUN] = x[XYL] - qflx_sun / (fx * params_inst.kmax[ivt_p, XYL] * laisun_p)
    else
        x[SUN] = x[XYL]
    end

    # Calculate soil flux
    soilflux = 0.0
    for j in 1:nlevsoi
        soilflux += k_soil_root_p[j] * (smp_c[j] - x[ROOT_SEG] - grav2[j])
    end

    return (x, soilflux)
end

# =====================================================================
# calcstress! — Compute transpiration stress via plant hydraulics
# =====================================================================

"""
    calcstress!(p, c, x, gb_mol, gs_mol_sun, gs_mol_sha, qsatl, qaf,
                k_soil_root_p, smp_c, z_c, laisun_p, laisha_p,
                htop_p, tsai_p, ivt_p, nlevsoi,
                elai_p, esai_p, fdry_p, forc_rho_c, forc_pbot_c, tgcm_p,
                qflx_tran_veg_p, vegwp_pd_p, is_night, londeg_g, local_time)

Compute transpiration stress using a plant hydraulics approach.
Returns `(x, bsun, bsha)` where x is updated vegwp.
Ported from `subroutine calcstress`.
"""
function calcstress!(p::Int, c::Int, x::Vector{Float64},
                     gb_mol::Float64, gs_mol_sun_in::Float64, gs_mol_sha_in::Float64,
                     qsatl::Float64, qaf::Float64,
                     k_soil_root_p::AbstractVector{Float64},
                     smp_c::AbstractVector{Float64},
                     z_c::AbstractVector{Float64},
                     laisun_p::Float64, laisha_p::Float64,
                     htop_p::Float64, tsai_p::Float64,
                     ivt_p::Int, nlevsoi::Int,
                     elai_p::Float64, esai_p::Float64,
                     fdry_p::Float64, forc_rho_c::Float64,
                     forc_pbot_c::Float64, tgcm_p::Float64)
    itmax = 50
    tolf = 1.0e-6
    toldx = 1.0e-9
    tol_lai = 0.001

    # Night time flag: vegwp(sun) > 0 signals night
    night = false
    if x[SUN] > 0.0
        night = true
        x[SUN] = x[SHA]
    end

    gs0sun = gs_mol_sun_in
    gs0sha = gs_mol_sha_in

    # Compute transpiration demand
    havegs = true
    qflx_sun, qflx_sha, _, _ = getqflx!(p, c, gb_mol, gs0sun, gs0sha,
        0.0, 0.0, qsatl, qaf, havegs, laisun_p, laisha_p,
        elai_p, esai_p, fdry_p, forc_rho_c, forc_pbot_c, tgcm_p)

    flag = false

    if (laisun_p > tol_lai || laisha_p > tol_lai) && (qflx_sun > 0.0 || qflx_sha > 0.0)
        # Newton's method
        iter = 0
        while true
            iter += 1

            f = spacF!(p, c, x, qflx_sun, qflx_sha,
                       k_soil_root_p, smp_c, z_c,
                       laisun_p, laisha_p, htop_p, tsai_p, ivt_p, nlevsoi)

            if sqrt(sum(f .* f)) < tolf * (qflx_sun + qflx_sha)
                flag = false
                break
            end
            if iter > itmax
                flag = false
                break
            end

            invA, flag = spacA!(p, c, x, qflx_sun, qflx_sha,
                                k_soil_root_p, laisun_p, laisha_p,
                                htop_p, tsai_p, ivt_p, nlevsoi)

            if flag
                break
            end

            if laisun_p > tol_lai && laisha_p > tol_lai
                dx = invA * f
            else
                dx = zeros(NVEGWCS)
                dx[SUN] = 0.0
                dx[SHA:ROOT_SEG] = invA[SHA:ROOT_SEG, SHA:ROOT_SEG] * f[SHA:ROOT_SEG]
            end

            if maximum(abs.(dx)) > 50000.0
                dx .= 50000.0 .* dx ./ maximum(abs.(dx))
            end

            if laisun_p > tol_lai && laisha_p > tol_lai
                x .= x .+ dx
            elseif laisha_p > tol_lai
                x .= x .+ dx
                x[SUN] = x[XYL]  # psi_sun = psi_xyl because laisun==0
            else
                x[XYL:ROOT_SEG] .= x[XYL:ROOT_SEG] .+ dx[XYL:ROOT_SEG]
                x[SUN] = x[SUN] + dx[SHA]
                x[SHA] = x[XYL]  # psi_sha = psi_xyl because laisha==0
            end

            if sqrt(sum(dx .* dx)) < toldx
                break
            end

            # Force SPAC gradient toward atmosphere
            if x[XYL] > x[ROOT_SEG]
                x[XYL] = x[ROOT_SEG]
            end
            if x[SUN] > x[XYL]
                x[SUN] = x[XYL]
            end
            if x[SHA] > x[XYL]
                x[SHA] = x[XYL]
            end
        end
    else
        flag = true
    end

    bsun = 0.0
    bsha = 0.0

    if flag
        # Solve algebraically
        x_new, soilflux = getvegwp!(p, c, gb_mol, gs0sun, gs0sha, qsatl, qaf,
            k_soil_root_p, smp_c, z_c, laisun_p, laisha_p,
            htop_p, tsai_p, ivt_p, nlevsoi,
            elai_p, esai_p, fdry_p, forc_rho_c, forc_pbot_c, tgcm_p)
        x .= x_new
        bsun = plc(x[SUN], ivt_p, SUN, VEG)
        bsha = plc(x[SHA], ivt_p, SHA, VEG)
    else
        # Compute attenuated flux
        qsun = qflx_sun * plc(x[SUN], ivt_p, SUN, VEG)
        qsha = qflx_sha * plc(x[SHA], ivt_p, SHA, VEG)

        # Retrieve stressed stomatal conductance
        havegs = false
        _, _, gs0sun_stressed, gs0sha_stressed = getqflx!(p, c, gb_mol, gs0sun, gs0sha,
            qsun, qsha, qsatl, qaf, havegs, laisun_p, laisha_p,
            elai_p, esai_p, fdry_p, forc_rho_c, forc_pbot_c, tgcm_p)

        # Compute water stress
        if qflx_sun > 0.0
            bsun = gs0sun_stressed / gs_mol_sun_in
        else
            bsun = plc(x[SUN], ivt_p, SUN, VEG)
        end
        if qflx_sha > 0.0
            bsha = gs0sha_stressed / gs_mol_sha_in
        else
            bsha = plc(x[SHA], ivt_p, SHA, VEG)
        end
    end

    if bsun < 0.01
        bsun = 0.0
    end
    if bsha < 0.01
        bsha = 0.0
    end

    # Night time: set vegwp and transpiration
    if night
        gs0sun_night = bsun * gs_mol_sun_in
        gs0sha_night = bsha * gs_mol_sha_in
        x_new, soilflux = getvegwp!(p, c, gb_mol, gs0sun_night, gs0sha_night, qsatl, qaf,
            k_soil_root_p, smp_c, z_c, laisun_p, laisha_p,
            htop_p, tsai_p, ivt_p, nlevsoi,
            elai_p, esai_p, fdry_p, forc_rho_c, forc_pbot_c, tgcm_p)
        x .= x_new
    end

    return (x, bsun, bsha)
end

# =====================================================================
# ci_func_PHS! — ci function evaluation for PHS
# =====================================================================

"""
    ci_func_PHS!(x, cisun, cisha, p, iv, c, bsun, bsha, bflag,
                 gb_mol, gs0sun, gs0sha, jesun, jesha, cair, oair,
                 lmr_z_sun, lmr_z_sha, par_z_sun, par_z_sha, rh_can,
                 qsatl, qaf, ps, forc_pbot_c,
                 k_soil_root_p, smp_c, z_c,
                 laisun_p, laisha_p, htop_p, tsai_p, ivt_p, nlevsoi,
                 elai_p, esai_p, fdry_p, forc_rho_c, tgcm_p;
                 medlynslope_val, medlynintercept_val)

Evaluate f(ci) = ci - (ca - (1.37rb+1.65rs))*patm*an for sunlit and shaded leaves.
Returns `(fvalsun, fvalsha, gs_mol_sun, gs_mol_sha, bsun, bsha)`.
Ported from `subroutine ci_func_PHS`.
"""
function ci_func_PHS!(x::Vector{Float64}, cisun::Float64, cisha::Float64,
                      p::Int, iv::Int, c::Int,
                      bsun::Float64, bsha::Float64, bflag::Bool,
                      gb_mol::Float64, gs0sun::Float64, gs0sha::Float64,
                      gs_mol_sun::Float64, gs_mol_sha::Float64,
                      jesun::Float64, jesha::Float64,
                      cair::Float64, oair::Float64,
                      lmr_z_sun::Float64, lmr_z_sha::Float64,
                      par_z_sun::Float64, par_z_sha::Float64,
                      rh_can::Float64, qsatl::Float64, qaf::Float64,
                      ps::PhotosynthesisData, forc_pbot_c::Float64,
                      k_soil_root_p::AbstractVector{Float64},
                      smp_c::AbstractVector{Float64},
                      z_c::AbstractVector{Float64},
                      laisun_p::Float64, laisha_p::Float64,
                      htop_p::Float64, tsai_p::Float64,
                      ivt_p::Int, nlevsoi::Int,
                      elai_p::Float64, esai_p::Float64,
                      fdry_p::Float64, forc_rho_c::Float64,
                      tgcm_p::Float64;
                      medlynslope_val::Float64=6.0,
                      medlynintercept_val::Float64=100.0)
    c3flag = ps.c3flag_patch[p]
    ac = ps.ac_phs_patch
    aj = ps.aj_phs_patch
    ap = ps.ap_phs_patch
    ag = ps.ag_phs_patch
    vcmax_z = ps.vcmax_z_phs_patch
    tpu_z = ps.tpu_z_phs_patch
    kp_z = ps.kp_z_phs_patch
    an_sun = ps.an_sun_patch
    an_sha = ps.an_sha_patch
    cp_p = ps.cp_patch[p]
    kc_p = ps.kc_patch[p]
    ko_p = ps.ko_patch[p]
    qe_p = ps.qe_patch[p]
    bbb_p = ps.bbb_patch[p]
    mbb_p = ps.mbb_patch[p]
    stomatalcond_mtd = ps.stomatalcond_mtd

    fvalsun = 0.0
    fvalsha = 0.0

    if bflag
        x, bsun, bsha = calcstress!(p, c, x, gb_mol, gs0sun, gs0sha, qsatl, qaf,
            k_soil_root_p, smp_c, z_c, laisun_p, laisha_p,
            htop_p, tsai_p, ivt_p, nlevsoi,
            elai_p, esai_p, fdry_p, forc_rho_c, forc_pbot_c, tgcm_p)
    end

    if c3flag
        # C3: Rubisco-limited
        ac[p, SUN, iv] = bsun * vcmax_z[p, SUN, iv] * max(cisun - cp_p, 0.0) / (cisun + kc_p * (1.0 + oair / ko_p))
        ac[p, SHA, iv] = bsha * vcmax_z[p, SHA, iv] * max(cisha - cp_p, 0.0) / (cisha + kc_p * (1.0 + oair / ko_p))
        # C3: RuBP-limited
        aj[p, SUN, iv] = jesun * max(cisun - cp_p, 0.0) / (4.0 * cisun + 8.0 * cp_p)
        aj[p, SHA, iv] = jesha * max(cisha - cp_p, 0.0) / (4.0 * cisha + 8.0 * cp_p)
        # C3: Product-limited
        ap[p, SUN, iv] = 3.0 * tpu_z[p, SUN, iv]
        ap[p, SHA, iv] = 3.0 * tpu_z[p, SHA, iv]
    else
        # C4: Rubisco-limited
        ac[p, SUN, iv] = bsun * vcmax_z[p, SUN, iv]
        ac[p, SHA, iv] = bsha * vcmax_z[p, SHA, iv]
        # C4: RuBP-limited
        aj[p, SUN, iv] = qe_p * par_z_sun * 4.6
        aj[p, SHA, iv] = qe_p * par_z_sha * 4.6
        # C4: PEP carboxylase-limited
        ap[p, SUN, iv] = kp_z[p, SUN, iv] * max(cisun, 0.0) / forc_pbot_c
        ap[p, SHA, iv] = kp_z[p, SHA, iv] * max(cisha, 0.0) / forc_pbot_c
    end

    # Gross photosynthesis - Sunlit
    aquad = params_inst.theta_cj[ivt_p]
    bquad = -(ac[p, SUN, iv] + aj[p, SUN, iv])
    cquad = ac[p, SUN, iv] * aj[p, SUN, iv]
    r1, r2 = quadratic_solve(aquad, bquad, cquad)
    ai = min(r1, r2)

    aquad = params_inst.theta_ip
    bquad = -(ai + ap[p, SUN, iv])
    cquad = ai * ap[p, SUN, iv]
    r1, r2 = quadratic_solve(aquad, bquad, cquad)
    ag[p, SUN, iv] = max(0.0, min(r1, r2))

    # Gross photosynthesis - Shaded
    aquad = params_inst.theta_cj[ivt_p]
    bquad = -(ac[p, SHA, iv] + aj[p, SHA, iv])
    cquad = ac[p, SHA, iv] * aj[p, SHA, iv]
    r1, r2 = quadratic_solve(aquad, bquad, cquad)
    ai = min(r1, r2)

    aquad = params_inst.theta_ip
    bquad = -(ai + ap[p, SHA, iv])
    cquad = ai * ap[p, SHA, iv]
    r1, r2 = quadratic_solve(aquad, bquad, cquad)
    ag[p, SHA, iv] = max(0.0, min(r1, r2))

    # Net photosynthesis
    an_sun[p, iv] = ag[p, SUN, iv] - bsun * lmr_z_sun
    an_sha[p, iv] = ag[p, SHA, iv] - bsha * lmr_z_sha

    if an_sun[p, iv] < 0.0
        if stomatalcond_mtd == STOMATALCOND_MTD_MEDLYN2011
            gs_mol_sun = medlynintercept_val
        elseif stomatalcond_mtd == STOMATALCOND_MTD_BB1987
            gs_mol_sun = bbb_p
        end
        gs_mol_sun = max(bsun * gs_mol_sun, 1.0)
        fvalsun = 0.0
    end
    if an_sha[p, iv] < 0.0
        if stomatalcond_mtd == STOMATALCOND_MTD_MEDLYN2011
            gs_mol_sha = medlynintercept_val
        elseif stomatalcond_mtd == STOMATALCOND_MTD_BB1987
            gs_mol_sha = bbb_p
        end
        gs_mol_sha = max(bsha * gs_mol_sha, 1.0)
        fvalsha = 0.0
    end
    if an_sun[p, iv] < 0.0 && an_sha[p, iv] < 0.0
        return (fvalsun, fvalsha, gs_mol_sun, gs_mol_sha, bsun, bsha)
    end

    # Quadratic gs_mol calculation with an known
    # Sunlit
    cs_sun = 0.0
    if an_sun[p, iv] >= 0.0
        cs_sun = cair - 1.4 / gb_mol * an_sun[p, iv] * forc_pbot_c
        cs_sun = max(cs_sun, MAX_CS)
    end

    if stomatalcond_mtd == STOMATALCOND_MTD_MEDLYN2011
        if an_sun[p, iv] >= 0.0
            term = 1.6 * an_sun[p, iv] / (cs_sun / forc_pbot_c * 1.0e06)
            aquad = 1.0
            bquad = -(2.0 * (medlynintercept_val * 1.0e-06 + term) +
                      (medlynslope_val * term)^2 / (gb_mol * 1.0e-06 * rh_can))
            cquad = medlynintercept_val^2 * 1.0e-12 +
                    (2.0 * medlynintercept_val * 1.0e-06 + term *
                     (1.0 - medlynslope_val^2 / rh_can)) * term
            r1, r2 = quadratic_solve(aquad, bquad, cquad)
            gs_mol_sun = max(r1, r2) * 1.0e06
        end

        # Shaded
        if an_sha[p, iv] >= 0.0
            cs_sha = cair - 1.4 / gb_mol * an_sha[p, iv] * forc_pbot_c
            cs_sha = max(cs_sha, MAX_CS)
            term = 1.6 * an_sha[p, iv] / (cs_sha / forc_pbot_c * 1.0e06)
            aquad = 1.0
            bquad = -(2.0 * (medlynintercept_val * 1.0e-06 + term) +
                      (medlynslope_val * term)^2 / (gb_mol * 1.0e-06 * rh_can))
            cquad = medlynintercept_val^2 * 1.0e-12 +
                    (2.0 * medlynintercept_val * 1.0e-06 + term *
                     (1.0 - medlynslope_val^2 / rh_can)) * term
            r1, r2 = quadratic_solve(aquad, bquad, cquad)
            gs_mol_sha = max(r1, r2) * 1.0e06
        end
    elseif stomatalcond_mtd == STOMATALCOND_MTD_BB1987
        if an_sun[p, iv] >= 0.0
            aquad = cs_sun
            bquad = cs_sun * (gb_mol - max(bsun * bbb_p, 1.0)) - mbb_p * an_sun[p, iv] * forc_pbot_c
            cquad = -gb_mol * (cs_sun * max(bsun * bbb_p, 1.0) + mbb_p * an_sun[p, iv] * forc_pbot_c * rh_can)
            r1, r2 = quadratic_solve(aquad, bquad, cquad)
            gs_mol_sun = max(r1, r2)
        end

        # Shaded
        if an_sha[p, iv] >= 0.0
            cs_sha = cair - 1.4 / gb_mol * an_sha[p, iv] * forc_pbot_c
            cs_sha = max(cs_sha, MAX_CS)
            aquad = cs_sha
            bquad = cs_sha * (gb_mol - max(bsha * bbb_p, 1.0)) - mbb_p * an_sha[p, iv] * forc_pbot_c
            cquad = -gb_mol * (cs_sha * max(bsha * bbb_p, 1.0) + mbb_p * an_sha[p, iv] * forc_pbot_c * rh_can)
            r1, r2 = quadratic_solve(aquad, bquad, cquad)
            gs_mol_sha = max(r1, r2)
        end
    end

    # Derive new estimate for cisun, cisha
    if an_sun[p, iv] >= 0.0
        if gs_mol_sun > 0.0
            fvalsun = cisun - cair + an_sun[p, iv] * forc_pbot_c * (1.4 * gs_mol_sun + 1.6 * gb_mol) / (gb_mol * gs_mol_sun)
        else
            fvalsun = cisun - cair
        end
    end
    if an_sha[p, iv] >= 0.0
        if gs_mol_sha > 0.0
            fvalsha = cisha - cair + an_sha[p, iv] * forc_pbot_c * (1.4 * gs_mol_sha + 1.6 * gb_mol) / (gb_mol * gs_mol_sha)
        else
            fvalsha = cisha - cair
        end
    end

    return (fvalsun, fvalsha, gs_mol_sun, gs_mol_sha, bsun, bsha)
end

# =====================================================================
# brent_PHS! — Brent's method for PHS (simultaneous sun/shade)
# =====================================================================

"""
    brent_PHS!(x1sun, x2sun, f1sun, f2sun, x1sha, x2sha, f1sha, f2sha,
               tol, p, iv, c, gb_mol, jesun, jesha, cair, oair,
               lmr_z_sun, lmr_z_sha, par_z_sun, par_z_sha, rh_can,
               gs_mol_sun, gs_mol_sha, bsun, bsha, qsatl, qaf,
               ps, forc_pbot_c, k_soil_root_p, smp_c, z_c,
               laisun_p, laisha_p, htop_p, tsai_p, ivt_p, nlevsoi,
               elai_p, esai_p, fdry_p, forc_rho_c, tgcm_p; kwargs...)

Brent's method for simultaneous sun/shade ci root finding.
Returns `(xsun, xsha, gs_mol_sun, gs_mol_sha, bsun, bsha)`.
Ported from `subroutine brent_PHS`.
"""
function brent_PHS!(x1sun::Float64, x2sun::Float64, f1sun::Float64, f2sun::Float64,
                    x1sha::Float64, x2sha::Float64, f1sha::Float64, f2sha::Float64,
                    tol::Float64, p::Int, iv::Int, c::Int,
                    gb_mol::Float64, jesun::Float64, jesha::Float64,
                    cair::Float64, oair::Float64,
                    lmr_z_sun::Float64, lmr_z_sha::Float64,
                    par_z_sun::Float64, par_z_sha::Float64,
                    rh_can::Float64,
                    gs_mol_sun::Float64, gs_mol_sha::Float64,
                    bsun::Float64, bsha::Float64,
                    qsatl::Float64, qaf::Float64,
                    ps::PhotosynthesisData, forc_pbot_c::Float64,
                    k_soil_root_p::AbstractVector{Float64},
                    smp_c::AbstractVector{Float64},
                    z_c::AbstractVector{Float64},
                    laisun_p::Float64, laisha_p::Float64,
                    htop_p::Float64, tsai_p::Float64,
                    ivt_p::Int, nlevsoi::Int,
                    elai_p::Float64, esai_p::Float64,
                    fdry_p::Float64, forc_rho_c::Float64,
                    tgcm_p::Float64;
                    medlynslope_val::Float64=6.0,
                    medlynintercept_val::Float64=100.0)
    nphs = 2
    itmax = 20
    eps_val = 1.0e-4
    bflag_const = false

    a = [x1sun, x1sha]
    b = [x2sun, x2sha]
    fa = [f1sun, f1sha]
    fb = [f2sun, f2sha]

    c_arr = copy(b)
    fc = copy(fb)
    d = b .- a
    e = copy(d)

    x_dummy = zeros(NVEGWCS)

    iter = 0
    while iter < itmax
        iter += 1

        for phase in 1:nphs
            if (fb[phase] > 0.0 && fc[phase] > 0.0) || (fb[phase] < 0.0 && fc[phase] < 0.0)
                c_arr[phase] = a[phase]
                fc[phase] = fa[phase]
                d[phase] = b[phase] - a[phase]
                e[phase] = d[phase]
            end
            if abs(fc[phase]) < abs(fb[phase])
                a[phase] = b[phase]
                b[phase] = c_arr[phase]
                c_arr[phase] = a[phase]
                fa[phase] = fb[phase]
                fb[phase] = fc[phase]
                fc[phase] = fa[phase]
            end
        end

        tol1 = 2.0 .* eps_val .* abs.(b) .+ 0.5 * tol
        xm = 0.5 .* (c_arr .- b)

        if (abs(xm[SUN]) <= tol1[SUN] || fb[SUN] == 0.0) &&
           (abs(xm[SHA]) <= tol1[SHA] || fb[SHA] == 0.0)
            return (b[SUN], b[SHA], gs_mol_sun, gs_mol_sha, bsun, bsha)
        end

        p_arr = zeros(nphs)
        q_arr = zeros(nphs)
        r_arr = zeros(nphs)
        s_arr = zeros(nphs)

        for phase in 1:nphs
            if abs(e[phase]) >= tol1[phase] && abs(fa[phase]) > abs(fb[phase])
                s_arr[phase] = fb[phase] / fa[phase]
                if a[phase] == c_arr[phase]
                    p_arr[phase] = 2.0 * xm[phase] * s_arr[phase]
                    q_arr[phase] = 1.0 - s_arr[phase]
                else
                    q_arr[phase] = fa[phase] / fc[phase]
                    r_arr[phase] = fb[phase] / fc[phase]
                    p_arr[phase] = s_arr[phase] * (2.0 * xm[phase] * q_arr[phase] * (q_arr[phase] - r_arr[phase]) -
                                   (b[phase] - a[phase]) * (r_arr[phase] - 1.0))
                    q_arr[phase] = (q_arr[phase] - 1.0) * (r_arr[phase] - 1.0) * (s_arr[phase] - 1.0)
                end
                if p_arr[phase] > 0.0
                    q_arr[phase] = -q_arr[phase]
                end
                p_arr[phase] = abs(p_arr[phase])
                if 2.0 * p_arr[phase] < min(3.0 * xm[phase] * q_arr[phase] - abs(tol1[phase] * q_arr[phase]),
                                              abs(e[phase] * q_arr[phase]))
                    e[phase] = d[phase]
                    d[phase] = p_arr[phase] / q_arr[phase]
                else
                    d[phase] = xm[phase]
                    e[phase] = d[phase]
                end
            else
                d[phase] = xm[phase]
                e[phase] = d[phase]
            end
            a[phase] = b[phase]
            fa[phase] = fb[phase]
            if abs(d[phase]) > tol1[phase]
                b[phase] = b[phase] + d[phase]
            else
                b[phase] = b[phase] + copysign(tol1[phase], xm[phase])
            end
        end

        gs0sun_local = gs_mol_sun
        gs0sha_local = gs_mol_sha
        fb[SUN], fb[SHA], gs_mol_sun, gs_mol_sha, bsun, bsha = ci_func_PHS!(
            x_dummy, b[SUN], b[SHA], p, iv, c,
            bsun, bsha, bflag_const,
            gb_mol, gs0sun_local, gs0sha_local, gs_mol_sun, gs_mol_sha,
            jesun, jesha, cair, oair,
            lmr_z_sun, lmr_z_sha, par_z_sun, par_z_sha,
            rh_can, qsatl, qaf,
            ps, forc_pbot_c,
            k_soil_root_p, smp_c, z_c,
            laisun_p, laisha_p, htop_p, tsai_p, ivt_p, nlevsoi,
            elai_p, esai_p, fdry_p, forc_rho_c, tgcm_p;
            medlynslope_val=medlynslope_val, medlynintercept_val=medlynintercept_val)

        if fb[SUN] == 0.0 && fb[SHA] == 0.0
            break
        end
    end

    return (b[SUN], b[SHA], gs_mol_sun, gs_mol_sha, bsun, bsha)
end

# =====================================================================
# hybrid_PHS! — Hybrid solver for PHS (simultaneous sun/shade ci)
# =====================================================================

"""
    hybrid_PHS!(x0sun, x0sha, p, iv, c, gb_mol, jesun, jesha,
                cair, oair, lmr_z_sun, lmr_z_sha,
                par_z_sun, par_z_sha, rh_can, qsatl, qaf,
                ps, forc_pbot_c, vegwp_p, k_soil_root_p, smp_c, z_c,
                laisun_p, laisha_p, htop_p, tsai_p, ivt_p, nlevsoi,
                elai_p, esai_p, fdry_p, forc_rho_c, tgcm_p,
                qflx_tran_veg_ref; kwargs...)

Hybrid solver for simultaneous sun/shade ci with plant hydraulic stress.
Returns `(x0sun, x0sha, gs_mol_sun, gs_mol_sha, bsun, bsha, vegwp, iter1, iter2)`.
Ported from `subroutine hybrid_PHS`.
"""
function hybrid_PHS!(x0sun::Float64, x0sha::Float64,
                     p::Int, iv::Int, c::Int,
                     gb_mol::Float64,
                     jesun::Float64, jesha::Float64,
                     cair::Float64, oair::Float64,
                     lmr_z_sun::Float64, lmr_z_sha::Float64,
                     par_z_sun::Float64, par_z_sha::Float64,
                     rh_can::Float64, qsatl::Float64, qaf::Float64,
                     ps::PhotosynthesisData, forc_pbot_c::Float64,
                     vegwp_p::Vector{Float64},
                     k_soil_root_p::AbstractVector{Float64},
                     smp_c::AbstractVector{Float64},
                     z_c::AbstractVector{Float64},
                     laisun_p::Float64, laisha_p::Float64,
                     htop_p::Float64, tsai_p::Float64,
                     ivt_p::Int, nlevsoi::Int,
                     elai_p::Float64, esai_p::Float64,
                     fdry_p::Float64, forc_rho_c::Float64,
                     tgcm_p::Float64;
                     medlynslope_val::Float64=6.0,
                     medlynintercept_val::Float64=100.0)
    toldb = 1.0e-2
    eps_val = 1.0e-2
    eps1 = 1.0e-4
    itmax = 3

    x1sun = x0sun
    x1sha = x0sha
    bflag = false
    b0sun = -1.0
    b0sha = -1.0
    gs0sun = 0.0
    gs0sha = 0.0
    bsun = 1.0
    bsha = 1.0
    iter1 = 0
    gs_mol_sun = 0.0
    gs_mol_sha = 0.0

    while true  # outer loop for bsun/bsha
        x = copy(vegwp_p)
        iter1 += 1
        iter2 = 0
        x0sun = max(0.1, x1sun)
        x1sun = 0.99 * x1sun
        x0sha = max(0.1, x1sha)
        x1sha = 0.99 * x1sha
        tolsun = abs(x1sun) * eps_val
        tolsha = abs(x1sha) * eps_val

        # First ci_func_PHS call: updates bsun/bsha (except first iter)
        f0sun, f0sha, gs_mol_sun, gs_mol_sha, bsun, bsha = ci_func_PHS!(
            x, x0sun, x0sha, p, iv, c, bsun, bsha, bflag,
            gb_mol, gs0sun, gs0sha, gs_mol_sun, gs_mol_sha,
            jesun, jesha, cair, oair,
            lmr_z_sun, lmr_z_sha, par_z_sun, par_z_sha,
            rh_can, qsatl, qaf,
            ps, forc_pbot_c,
            k_soil_root_p, smp_c, z_c,
            laisun_p, laisha_p, htop_p, tsai_p, ivt_p, nlevsoi,
            elai_p, esai_p, fdry_p, forc_rho_c, tgcm_p;
            medlynslope_val=medlynslope_val, medlynintercept_val=medlynintercept_val)

        dbsun = b0sun - bsun
        dbsha = b0sha - bsha
        b0sun = bsun
        b0sha = bsha
        bflag = false

        # Second ci_func_PHS call: create second point
        f1sun, f1sha, gs_mol_sun, gs_mol_sha, bsun, bsha = ci_func_PHS!(
            x, x1sun, x1sha, p, iv, c, bsun, bsha, bflag,
            gb_mol, gs0sun, gs0sha, gs_mol_sun, gs_mol_sha,
            jesun, jesha, cair, oair,
            lmr_z_sun, lmr_z_sha, par_z_sun, par_z_sha,
            rh_can, qsatl, qaf,
            ps, forc_pbot_c,
            k_soil_root_p, smp_c, z_c,
            laisun_p, laisha_p, htop_p, tsai_p, ivt_p, nlevsoi,
            elai_p, esai_p, fdry_p, forc_rho_c, tgcm_p;
            medlynslope_val=medlynslope_val, medlynintercept_val=medlynintercept_val)

        minf = abs(f1sun + f1sha)
        minxsun = x1sun
        minxsha = x1sha

        while true  # inner loop for ci
            if abs(f0sun) < eps1 && abs(f0sha) < eps1
                x1sun = x0sun
                x1sha = x0sha
                break
            end
            if abs(f1sun) < eps1 && abs(f1sha) < eps1
                break
            end
            iter2 += 1

            if (f1sun - f0sun) == 0.0
                dxsun = 0.5 * (x1sun + x0sun) - x1sun
            else
                dxsun = -f1sun * (x1sun - x0sun) / (f1sun - f0sun)
            end
            if (f1sha - f0sha) == 0.0
                dxsha = 0.5 * (x1sha + x0sha) - x1sha
            else
                dxsha = -f1sha * (x1sha - x0sha) / (f1sha - f0sha)
            end
            x0sun = x1sun
            x1sun = x1sun + dxsun
            x0sha = x1sha
            x1sha = x1sha + dxsha

            f1sun, f1sha, gs_mol_sun, gs_mol_sha, bsun, bsha = ci_func_PHS!(
                x, x1sun, x1sha, p, iv, c, bsun, bsha, bflag,
                gb_mol, gs0sun, gs0sha, gs_mol_sun, gs_mol_sha,
                jesun, jesha, cair, oair,
                lmr_z_sun, lmr_z_sha, par_z_sun, par_z_sha,
                rh_can, qsatl, qaf,
                ps, forc_pbot_c,
                k_soil_root_p, smp_c, z_c,
                laisun_p, laisha_p, htop_p, tsai_p, ivt_p, nlevsoi,
                elai_p, esai_p, fdry_p, forc_rho_c, tgcm_p;
                medlynslope_val=medlynslope_val, medlynintercept_val=medlynintercept_val)

            if abs(dxsun) < tolsun && abs(dxsha) < tolsha
                x0sun = x1sun
                x0sha = x1sha
                break
            end

            if iter2 == 1
                minf = abs(f1sun + f1sha)
                minxsun = x1sun
                minxsha = x1sha
            else
                if abs(f1sun + f1sha) < minf
                    minf = abs(f1sun + f1sha)
                    minxsun = x1sun
                    minxsha = x1sha
                end
            end

            if abs(f1sun) < eps1 && abs(f1sha) < eps1
                break
            end

            if f1sun * f0sun < 0.0 && f1sha * f0sha < 0.0
                xsun, xsha, gs_mol_sun, gs_mol_sha, bsun, bsha = brent_PHS!(
                    x0sun, x1sun, f0sun, f1sun,
                    x0sha, x1sha, f0sha, f1sha,
                    tolsun, p, iv, c,
                    gb_mol, jesun, jesha, cair, oair,
                    lmr_z_sun, lmr_z_sha, par_z_sun, par_z_sha,
                    rh_can, gs_mol_sun, gs_mol_sha, bsun, bsha,
                    qsatl, qaf,
                    ps, forc_pbot_c,
                    k_soil_root_p, smp_c, z_c,
                    laisun_p, laisha_p, htop_p, tsai_p, ivt_p, nlevsoi,
                    elai_p, esai_p, fdry_p, forc_rho_c, tgcm_p;
                    medlynslope_val=medlynslope_val, medlynintercept_val=medlynintercept_val)
                x0sun = xsun
                x0sha = xsha
                break
            end

            if iter2 > itmax
                x1sun = minxsun
                x1sha = minxsha
                f1sun, f1sha, gs_mol_sun, gs_mol_sha, bsun, bsha = ci_func_PHS!(
                    x, x1sun, x1sha, p, iv, c, bsun, bsha, bflag,
                    gb_mol, gs0sun, gs0sha, gs_mol_sun, gs_mol_sha,
                    jesun, jesha, cair, oair,
                    lmr_z_sun, lmr_z_sha, par_z_sun, par_z_sha,
                    rh_can, qsatl, qaf,
                    ps, forc_pbot_c,
                    k_soil_root_p, smp_c, z_c,
                    laisun_p, laisha_p, htop_p, tsai_p, ivt_p, nlevsoi,
                    elai_p, esai_p, fdry_p, forc_rho_c, tgcm_p;
                    medlynslope_val=medlynslope_val, medlynintercept_val=medlynintercept_val)
                break
            end
        end  # inner loop

        # Update unstressed stomatal conductance
        if bsun > 0.01
            gs0sun = gs_mol_sun / bsun
        end
        if bsha > 0.01
            gs0sha = gs_mol_sha / bsha
        end
        bflag = true

        if abs(dbsun) < toldb && abs(dbsha) < toldb
            break
        end
        if iter1 > itmax
            break
        end
    end  # outer loop

    x0sun = x1sun
    x0sha = x1sha

    # Set vegwp for final solution
    vegwp_new, soilflux = getvegwp!(p, c, gb_mol, gs_mol_sun, gs_mol_sha, qsatl, qaf,
        k_soil_root_p, smp_c, z_c, laisun_p, laisha_p,
        htop_p, tsai_p, ivt_p, nlevsoi,
        elai_p, esai_p, fdry_p, forc_rho_c, forc_pbot_c, tgcm_p)
    vegwp_p .= vegwp_new

    if soilflux < 0.0
        soilflux = 0.0
    end

    return (x0sun, x0sha, gs_mol_sun, gs_mol_sha, bsun, bsha, vegwp_p, soilflux, iter1, iter2)
end

# =====================================================================
# photosynthesis_hydrstress! — Main PHS photosynthesis
# =====================================================================

"""
    photosynthesis_hydrstress!(ps::PhotosynthesisData, ...)

Leaf photosynthesis and stomatal conductance with plant hydraulic stress.
Sunlit and shaded photosynthesis and stomatal conductance are solved
simultaneously per Pierre Gentine/Daniel Kennedy PHS method.
Ported from `subroutine PhotosynthesisHydraulicStress`.
"""
function photosynthesis_hydrstress!(ps::PhotosynthesisData,
                                    esat_tv::Vector{Float64},
                                    eair::Vector{Float64},
                                    oair::Vector{Float64},
                                    cair::Vector{Float64},
                                    rb::Vector{Float64},
                                    btran::Vector{Float64},
                                    dayl_factor::Vector{Float64},
                                    leafn::Vector{Float64},
                                    qsatl::Vector{Float64},
                                    qaf::Vector{Float64},
                                    forc_pbot::Vector{Float64},
                                    forc_rho::Vector{Float64},
                                    t_veg::Vector{Float64},
                                    t10::Vector{Float64},
                                    tgcm::Vector{Float64},
                                    nrad::Vector{Int},
                                    tlai_z::Matrix{Float64},
                                    tlai::Vector{Float64},
                                    tsai::Vector{Float64},
                                    par_z_sun_in::Matrix{Float64},
                                    par_z_sha_in::Matrix{Float64},
                                    lai_z_sun_in::Matrix{Float64},
                                    lai_z_sha_in::Matrix{Float64},
                                    vcmaxcint_sun::Vector{Float64},
                                    vcmaxcint_sha::Vector{Float64},
                                    o3coefv_sun::Vector{Float64},
                                    o3coefg_sun::Vector{Float64},
                                    o3coefv_sha::Vector{Float64},
                                    o3coefg_sha::Vector{Float64},
                                    c3psn_pft::Vector{Float64},
                                    leafcn_pft::Vector{Float64},
                                    flnr_pft::Vector{Float64},
                                    fnitr_pft::Vector{Float64},
                                    slatop_pft::Vector{Float64},
                                    mbbopt_pft::Vector{Float64},
                                    medlynintercept_pft::Vector{Float64},
                                    medlynslope_pft::Vector{Float64},
                                    froot_leaf_pft::Vector{Float64},
                                    root_radius_pft::Vector{Float64},
                                    root_density_pft::Vector{Float64},
                                    crop_pft::Vector{Float64},
                                    ivt::Vector{Int},
                                    col_of_patch::Vector{Int},
                                    mask_patch::BitVector,
                                    bounds_patch::UnitRange{Int},
                                    froot_carbon::Vector{Float64},
                                    croot_carbon::Vector{Float64},
                                    k_soil_root::Matrix{Float64},
                                    root_conductance_out::Matrix{Float64},
                                    soil_conductance_out::Matrix{Float64},
                                    rootfr::Matrix{Float64},
                                    dz::Matrix{Float64},
                                    z_col::Matrix{Float64},
                                    hk_l::Matrix{Float64},
                                    hksat::Matrix{Float64},
                                    smp_l::Matrix{Float64},
                                    vegwp::Matrix{Float64},
                                    vegwp_ln::Matrix{Float64},
                                    laisun::Vector{Float64},
                                    laisha::Vector{Float64},
                                    elai::Vector{Float64},
                                    esai::Vector{Float64},
                                    htop::Vector{Float64},
                                    fdry::Vector{Float64},
                                    qflx_tran_veg::Vector{Float64};
                                    nlevcan::Int=NLEVCAN,
                                    nlevsoi::Int=varpar.nlevsoi,
                                    use_cn::Bool=false,
                                    use_luna::Bool=false,
                                    use_c13::Bool=false,
                                    leaf_mr_vcm::Float64=0.015,
                                    is_near_local_noon_fn::Function=(p) -> false)
    stomatalcond_mtd = ps.stomatalcond_mtd
    leafresp_method = ps.leafresp_method
    light_inhibit = ps.light_inhibit
    modifyphoto_and_lmr_forcrop = ps.modifyphoto_and_lmr_forcrop
    ac = ps.ac_phs_patch
    aj = ps.aj_phs_patch
    ap = ps.ap_phs_patch
    ag = ps.ag_phs_patch
    vcmax_z = ps.vcmax_z_phs_patch
    tpu_z = ps.tpu_z_phs_patch
    kp_z = ps.kp_z_phs_patch
    an_sun = ps.an_sun_patch
    an_sha = ps.an_sha_patch
    gs_mol_sun = ps.gs_mol_sun_patch
    gs_mol_sha = ps.gs_mol_sha_patch
    gs_mol_sun_ln = ps.gs_mol_sun_ln_patch
    gs_mol_sha_ln = ps.gs_mol_sha_ln_patch
    c3flag = ps.c3flag_patch
    cp = ps.cp_patch
    kc = ps.kc_patch
    ko = ps.ko_patch
    qe = ps.qe_patch
    bbb = ps.bbb_patch
    mbb = ps.mbb_patch
    gb_mol_arr = ps.gb_mol_patch
    rh_leaf = ps.rh_leaf_patch
    vpd_can = ps.vpd_can_patch
    lnc = ps.lnca_patch

    croot_lateral_length = 0.25
    c_to_b = 2.0

    lmrc = fth25_photo(params_inst.lmrhd, params_inst.lmrse)

    np = length(bounds_patch)
    jmax_z_local = zeros(np, 2, nlevcan)
    bbbopt = zeros(np)
    kn = zeros(np)
    psn_wc_z_sun = zeros(np, nlevcan)
    psn_wj_z_sun = zeros(np, nlevcan)
    psn_wp_z_sun = zeros(np, nlevcan)
    psn_wc_z_sha = zeros(np, nlevcan)
    psn_wj_z_sha = zeros(np, nlevcan)
    psn_wp_z_sha = zeros(np, nlevcan)
    rh_leaf_sun = zeros(np)
    rh_leaf_sha = zeros(np)
    bsun_arr = ones(np)
    bsha_arr = ones(np)

    # Aliases for output
    ci_z_sun = ps.cisun_z_patch
    ci_z_sha = ps.cisha_z_patch
    rs_sun = ps.rssun_patch
    rs_z_sun = ps.rssun_z_patch
    lmr_sun = ps.lmrsun_patch
    lmr_z_sun = ps.lmrsun_z_patch
    psn_sun = ps.psnsun_patch
    psn_z_sun = ps.psnsun_z_patch
    psn_wc_sun = ps.psnsun_wc_patch
    psn_wj_sun = ps.psnsun_wj_patch
    psn_wp_sun = ps.psnsun_wp_patch
    rs_sha = ps.rssha_patch
    rs_z_sha = ps.rssha_z_patch
    lmr_sha = ps.lmrsha_patch
    lmr_z_sha = ps.lmrsha_z_patch
    psn_sha = ps.psnsha_patch
    psn_z_sha = ps.psnsha_z_patch
    psn_wc_sha = ps.psnsha_wc_patch
    psn_wj_sha = ps.psnsha_wj_patch
    psn_wp_sha = ps.psnsha_wp_patch

    rsmax0 = 2.0e4

    # ---- Pass 1: Root-soil conductance ----
    for p in bounds_patch
        mask_patch[p] || continue
        c = col_of_patch[p]
        ivt_p = ivt[p]

        for j in 1:nlevsoi
            root_biomass_density = c_to_b * froot_carbon[p] * rootfr[p, j] / dz[c, j]
            root_biomass_density = max(c_to_b * 1.0, root_biomass_density)

            root_cross_sec_area = RPI * root_radius_pft[ivt_p]^2
            root_length_density = root_biomass_density / (root_density_pft[ivt_p] * root_cross_sec_area)

            rai_j = (tsai[p] + tlai[p]) * froot_leaf_pft[ivt_p] * rootfr[p, j]

            croot_average_length = croot_lateral_length
            r_soil = sqrt(1.0 / (RPI * root_length_density))
            soil_cond = min(hksat[c, j], hk_l[c, j]) / (1.0e3 * r_soil)

            fs_j = plc(smp_l[c, j], ivt_p, ROOT_SEG, VEG)
            root_cond = (fs_j * rai_j * params_inst.krmax[ivt_p]) / (croot_average_length + z_col[c, j])

            soil_cond = max(soil_cond, 1.0e-16)
            root_cond = max(root_cond, 1.0e-16)

            root_conductance_out[p, j] = root_cond
            soil_conductance_out[p, j] = soil_cond

            rs_resis = 1.0 / soil_cond + 1.0 / root_cond

            if rai_j * rootfr[p, j] > 0.0 && j > 1
                k_soil_root[p, j] = 1.0 / rs_resis
            else
                k_soil_root[p, j] = 0.0
            end
        end
    end

    # ---- Pass 2: Kinetics, N profile, vcmax, respiration ----
    for p in bounds_patch
        mask_patch[p] || continue
        c = col_of_patch[p]
        ivt_p = ivt[p]

        # C3/C4 flag
        if round(Int, c3psn_pft[ivt_p]) == 1
            c3flag[p] = true
        else
            c3flag[p] = false
        end

        if c3flag[p]
            qe[p] = 0.0
            if stomatalcond_mtd == STOMATALCOND_MTD_BB1987
                bbbopt[p] = BBBOPT_C3
            end
        else
            qe[p] = 0.05
            if stomatalcond_mtd == STOMATALCOND_MTD_BB1987
                bbbopt[p] = BBBOPT_C4
            end
        end

        if stomatalcond_mtd == STOMATALCOND_MTD_BB1987
            bbb[p] = bbbopt[p]
            mbb[p] = mbbopt_pft[ivt_p]
        end

        kc25 = params_inst.kc25_coef * forc_pbot[c]
        ko25 = params_inst.ko25_coef * forc_pbot[c]
        sco = 0.5 * 0.209 / params_inst.cp25_yr2000
        cp25 = 0.5 * oair[p] / sco

        kc[p] = kc25 * ft_photo(t_veg[p], params_inst.kcha)
        ko[p] = ko25 * ft_photo(t_veg[p], params_inst.koha)
        cp[p] = cp25 * ft_photo(t_veg[p], params_inst.cpha)

        # Nitrogen profile
        lnc[p] = 1.0 / (slatop_pft[ivt_p] * leafcn_pft[ivt_p])
        lnc[p] = min(lnc[p], 10.0)

        vcmax25top = lnc[p] * flnr_pft[ivt_p] * params_inst.fnr * params_inst.act25 * dayl_factor[p]
        if !use_cn
            vcmax25top = vcmax25top * fnitr_pft[ivt_p]
        end

        jmax25top = ((2.59 - 0.035 * min(max((t10[p] - TFRZ), 11.0), 35.0)) * vcmax25top) *
                    params_inst.jmax25top_sf
        tpu25top = params_inst.tpu25ratio * vcmax25top
        kp25top = params_inst.kp25ratio * vcmax25top

        ps.luvcmax25top_patch[p] = vcmax25top
        ps.lujmax25top_patch[p] = jmax25top
        ps.lutpu25top_patch[p] = tpu25top

        if dayl_factor[p] < 1.0e-12
            kn[p] = 0.0
        else
            kn[p] = exp(0.00963 * vcmax25top / dayl_factor[p] - 2.43)
        end

        if use_cn
            if leafresp_method == LEAFRESP_MTD_RYAN1991
                lmr25top = 2.525e-6 * (1.5^((25.0 - 20.0) / 10.0))
                lmr25top = lmr25top * lnc[p] / 12.0e-06
            else
                lmr25top = 0.0
            end
        else
            if c3flag[p]
                lmr25top = vcmax25top * leaf_mr_vcm
            else
                lmr25top = vcmax25top * 0.025
            end
        end

        # Canopy layer loop
        laican = 0.0
        for iv in 1:nrad[p]
            if iv == 1
                laican = 0.5 * tlai_z[p, iv]
            else
                laican = laican + 0.5 * (tlai_z[p, iv-1] + tlai_z[p, iv])
            end

            if nlevcan == 1
                nscaler_sun = vcmaxcint_sun[p]
                nscaler_sha = vcmaxcint_sha[p]
            else
                nscaler_sun = exp(-kn[p] * laican)
                nscaler_sha = exp(-kn[p] * laican)
            end

            # Maintenance respiration
            lmr25_sun = lmr25top * nscaler_sun
            lmr25_sha = lmr25top * nscaler_sha

            if c3flag[p]
                lmr_z_sun[p, iv] = lmr25_sun * ft_photo(t_veg[p], params_inst.lmrha) *
                                   fth_photo(t_veg[p], params_inst.lmrhd, params_inst.lmrse, lmrc)
                lmr_z_sha[p, iv] = lmr25_sha * ft_photo(t_veg[p], params_inst.lmrha) *
                                   fth_photo(t_veg[p], params_inst.lmrhd, params_inst.lmrse, lmrc)
            else
                lmr_z_sun[p, iv] = lmr25_sun * 2.0^((t_veg[p] - (TFRZ + 25.0)) / 10.0)
                lmr_z_sun[p, iv] = lmr_z_sun[p, iv] / (1.0 + exp(1.3 * (t_veg[p] - (TFRZ + 55.0))))
                lmr_z_sha[p, iv] = lmr25_sha * 2.0^((t_veg[p] - (TFRZ + 25.0)) / 10.0)
                lmr_z_sha[p, iv] = lmr_z_sha[p, iv] / (1.0 + exp(1.3 * (t_veg[p] - (TFRZ + 55.0))))
            end

            # Reduce lmr with low LAI
            lmr_z_sun[p, iv] *= min(0.2 * exp(3.218 * tlai_z[p, iv]), 1.0)
            lmr_z_sha[p, iv] *= min(0.2 * exp(3.218 * tlai_z[p, iv]), 1.0)

            if par_z_sun_in[p, iv] <= 0.0  # night time
                vcmax_z[p, SUN, iv] = 0.0
                jmax_z_local[p, SUN, iv] = 0.0
                tpu_z[p, SUN, iv] = 0.0
                kp_z[p, SUN, iv] = 0.0
                vcmax_z[p, SHA, iv] = 0.0
                jmax_z_local[p, SHA, iv] = 0.0
                tpu_z[p, SHA, iv] = 0.0
                kp_z[p, SHA, iv] = 0.0

                if use_c13
                    ps.alphapsnsun_patch[p] = 1.0
                    ps.alphapsnsha_patch[p] = 1.0
                end
            else  # day time
                vcmax25_sun = vcmax25top * nscaler_sun
                jmax25_sun = jmax25top * nscaler_sun
                tpu25_sun = tpu25top * nscaler_sun
                vcmax25_sha = vcmax25top * nscaler_sha
                jmax25_sha = jmax25top * nscaler_sha
                tpu25_sha = tpu25top * nscaler_sha
                kp25_sun = kp25top * nscaler_sun
                kp25_sha = kp25top * nscaler_sha

                vcmaxse = (668.39 - 1.07 * min(max((t10[p] - TFRZ), 11.0), 35.0)) * params_inst.vcmaxse_sf
                jmaxse = (659.70 - 0.75 * min(max((t10[p] - TFRZ), 11.0), 35.0)) * params_inst.jmaxse_sf
                tpuse = (668.39 - 1.07 * min(max((t10[p] - TFRZ), 11.0), 35.0)) * params_inst.tpuse_sf
                vcmaxc = fth25_photo(params_inst.vcmaxhd, vcmaxse)
                jmaxc = fth25_photo(params_inst.jmaxhd, jmaxse)
                tpuc = fth25_photo(params_inst.tpuhd, tpuse)

                vcmax_z[p, SUN, iv] = vcmax25_sun * ft_photo(t_veg[p], params_inst.vcmaxha) *
                                      fth_photo(t_veg[p], params_inst.vcmaxhd, vcmaxse, vcmaxc)
                jmax_z_local[p, SUN, iv] = jmax25_sun * ft_photo(t_veg[p], params_inst.jmaxha) *
                                            fth_photo(t_veg[p], params_inst.jmaxhd, jmaxse, jmaxc)
                tpu_z[p, SUN, iv] = tpu25_sun * ft_photo(t_veg[p], params_inst.tpuha) *
                                    fth_photo(t_veg[p], params_inst.tpuhd, tpuse, tpuc)

                vcmax_z[p, SHA, iv] = vcmax25_sha * ft_photo(t_veg[p], params_inst.vcmaxha) *
                                      fth_photo(t_veg[p], params_inst.vcmaxhd, vcmaxse, vcmaxc)
                jmax_z_local[p, SHA, iv] = jmax25_sha * ft_photo(t_veg[p], params_inst.jmaxha) *
                                            fth_photo(t_veg[p], params_inst.jmaxhd, jmaxse, jmaxc)
                tpu_z[p, SHA, iv] = tpu25_sha * ft_photo(t_veg[p], params_inst.tpuha) *
                                    fth_photo(t_veg[p], params_inst.tpuhd, tpuse, tpuc)

                if !c3flag[p]
                    vcmax_z[p, SUN, iv] = vcmax25_sun * 2.0^((t_veg[p] - (TFRZ + 25.0)) / 10.0)
                    vcmax_z[p, SUN, iv] /= (1.0 + exp(0.2 * ((TFRZ + 15.0) - t_veg[p])))
                    vcmax_z[p, SUN, iv] /= (1.0 + exp(0.3 * (t_veg[p] - (TFRZ + 40.0))))
                    vcmax_z[p, SHA, iv] = vcmax25_sha * 2.0^((t_veg[p] - (TFRZ + 25.0)) / 10.0)
                    vcmax_z[p, SHA, iv] /= (1.0 + exp(0.2 * ((TFRZ + 15.0) - t_veg[p])))
                    vcmax_z[p, SHA, iv] /= (1.0 + exp(0.3 * (t_veg[p] - (TFRZ + 40.0))))
                end

                kp_z[p, SUN, iv] = kp25_sun * 2.0^((t_veg[p] - (TFRZ + 25.0)) / 10.0)
                kp_z[p, SHA, iv] = kp25_sha * 2.0^((t_veg[p] - (TFRZ + 25.0)) / 10.0)
            end

            # Light inhibition
            if light_inhibit && par_z_sun_in[p, 1] > 0.0
                lmr_z_sun[p, iv] *= 0.67
            end
            if light_inhibit && par_z_sha_in[p, 1] > 0.0
                lmr_z_sha[p, iv] *= 0.67
            end
        end
    end

    # ---- Pass 3: Leaf-level photosynthesis ----
    for p in bounds_patch
        mask_patch[p] || continue
        c = col_of_patch[p]
        ivt_p = ivt[p]

        cf = forc_pbot[c] / (RGAS * 1.0e-3 * tgcm[p]) * 1.0e06
        gb = 1.0 / rb[p]
        gb_mol_arr[p] = gb * cf

        for iv in 1:nrad[p]
            if par_z_sun_in[p, iv] <= 0.0  # night time
                vegwp[p, SUN] = 1.0  # signal for night

                if stomatalcond_mtd == STOMATALCOND_MTD_BB1987
                    gsminsun = bbb[p]
                    gsminsha = bbb[p]
                elseif stomatalcond_mtd == STOMATALCOND_MTD_MEDLYN2011
                    gsminsun = medlynintercept_pft[ivt_p]
                    gsminsha = medlynintercept_pft[ivt_p]
                end

                vegwp_view = vegwp[p, :]
                _, bsun_val, bsha_val = calcstress!(p, c, vegwp_view,
                    gb_mol_arr[p], gsminsun, gsminsha, qsatl[p], qaf[p],
                    k_soil_root[p, :], smp_l[c, :], z_col[c, :],
                    laisun[p], laisha[p], htop[p], tsai[p], ivt_p, nlevsoi,
                    elai[p], esai[p], fdry[p], forc_rho[c], forc_pbot[c], tgcm[p])
                vegwp[p, :] .= vegwp_view
                bsun_arr[p] = bsun_val
                bsha_arr[p] = bsha_val

                ac[p, SUN, iv] = 0.0
                aj[p, SUN, iv] = 0.0
                ap[p, SUN, iv] = 0.0
                ag[p, SUN, iv] = 0.0
                if round(Int, crop_pft[ivt_p]) == 0 || !modifyphoto_and_lmr_forcrop
                    an_sun[p, iv] = ag[p, SUN, iv] - bsun_arr[p] * lmr_z_sun[p, iv]
                else
                    an_sun[p, iv] = ag[p, SUN, iv] - lmr_z_sun[p, iv]
                end
                psn_z_sun[p, iv] = 0.0
                psn_wc_z_sun[p, iv] = 0.0
                psn_wj_z_sun[p, iv] = 0.0
                psn_wp_z_sun[p, iv] = 0.0
                rs_z_sun[p, iv] = min(rsmax0, 1.0 / max(bsun_arr[p] * gsminsun, 1.0) * cf)
                ci_z_sun[p, iv] = 0.0
                rh_leaf_sun[p] = 0.0

                ac[p, SHA, iv] = 0.0
                aj[p, SHA, iv] = 0.0
                ap[p, SHA, iv] = 0.0
                ag[p, SHA, iv] = 0.0
                if round(Int, crop_pft[ivt_p]) == 0 || !modifyphoto_and_lmr_forcrop
                    an_sha[p, iv] = ag[p, SHA, iv] - bsha_arr[p] * lmr_z_sha[p, iv]
                else
                    an_sha[p, iv] = ag[p, SHA, iv] - lmr_z_sha[p, iv]
                end
                psn_z_sha[p, iv] = 0.0
                psn_wc_z_sha[p, iv] = 0.0
                psn_wj_z_sha[p, iv] = 0.0
                psn_wp_z_sha[p, iv] = 0.0
                rs_z_sha[p, iv] = min(rsmax0, 1.0 / max(bsha_arr[p] * gsminsha, 1.0) * cf)
                ci_z_sha[p, iv] = 0.0
                rh_leaf_sha[p] = 0.0

            else  # day time
                ceair = min(eair[p], esat_tv[p])
                if stomatalcond_mtd == STOMATALCOND_MTD_BB1987
                    rh_can = ceair / esat_tv[p]
                elseif stomatalcond_mtd == STOMATALCOND_MTD_MEDLYN2011
                    rh_can = max((esat_tv[p] - ceair), MEDLYN_RH_CAN_MAX) * MEDLYN_RH_CAN_FACT
                    vpd_can[p] = rh_can
                end

                # Electron transport - Sun
                qabs = 0.5 * (1.0 - params_inst.fnps) * par_z_sun_in[p, iv] * 4.6
                aquad = params_inst.theta_psii
                bquad = -(qabs + jmax_z_local[p, SUN, iv])
                cquad = qabs * jmax_z_local[p, SUN, iv]
                r1, r2 = quadratic_solve(aquad, bquad, cquad)
                je_sun = min(r1, r2)

                # Electron transport - Shade
                qabs = 0.5 * (1.0 - params_inst.fnps) * par_z_sha_in[p, iv] * 4.6
                aquad = params_inst.theta_psii
                bquad = -(qabs + jmax_z_local[p, SHA, iv])
                cquad = qabs * jmax_z_local[p, SHA, iv]
                r1, r2 = quadratic_solve(aquad, bquad, cquad)
                je_sha = min(r1, r2)

                # Initial ci guess
                if c3flag[p]
                    ci_z_sun[p, iv] = 0.7 * cair[p]
                    ci_z_sha[p, iv] = 0.7 * cair[p]
                else
                    ci_z_sun[p, iv] = 0.4 * cair[p]
                    ci_z_sha[p, iv] = 0.4 * cair[p]
                end

                # Solve for ci and gs via hybrid_PHS
                vegwp_view = vegwp[p, :]
                ci_z_sun[p, iv], ci_z_sha[p, iv], gs_mol_sun_val, gs_mol_sha_val,
                    bsun_val, bsha_val, _, soilflux, _, _ = hybrid_PHS!(
                    ci_z_sun[p, iv], ci_z_sha[p, iv], p, iv, c,
                    gb_mol_arr[p], je_sun, je_sha, cair[p], oair[p],
                    lmr_z_sun[p, iv], lmr_z_sha[p, iv],
                    par_z_sun_in[p, iv], par_z_sha_in[p, iv],
                    rh_can, qsatl[p], qaf[p],
                    ps, forc_pbot[c], vegwp_view,
                    k_soil_root[p, :], smp_l[c, :], z_col[c, :],
                    laisun[p], laisha[p], htop[p], tsai[p], ivt_p, nlevsoi,
                    elai[p], esai[p], fdry[p], forc_rho[c], tgcm[p];
                    medlynslope_val=medlynslope_pft[ivt_p],
                    medlynintercept_val=medlynintercept_pft[ivt_p])
                vegwp[p, :] .= vegwp_view
                gs_mol_sun[p, iv] = gs_mol_sun_val
                gs_mol_sha[p, iv] = gs_mol_sha_val
                bsun_arr[p] = bsun_val
                bsha_arr[p] = bsha_val
                qflx_tran_veg[p] = max(soilflux, 0.0)

                # Determine gs min/slope for error checking
                if stomatalcond_mtd == STOMATALCOND_MTD_MEDLYN2011
                    gsminsun = medlynintercept_pft[ivt_p]
                    gsminsha = medlynintercept_pft[ivt_p]
                elseif stomatalcond_mtd == STOMATALCOND_MTD_BB1987
                    gsminsun = bbb[p]
                    gsminsha = bbb[p]
                end

                # Check an < 0
                if an_sun[p, iv] < 0.0
                    gs_mol_sun[p, iv] = max(bsun_arr[p] * gsminsun, 1.0)
                end
                if an_sha[p, iv] < 0.0
                    gs_mol_sha[p, iv] = max(bsha_arr[p] * gsminsha, 1.0)
                end

                # Local noon gs
                if is_near_local_noon_fn(p)
                    gs_mol_sun_ln[p, iv] = gs_mol_sun[p, iv]
                    gs_mol_sha_ln[p, iv] = gs_mol_sha[p, iv]
                    vegwp_ln[p, :] .= vegwp[p, :]
                else
                    gs_mol_sun_ln[p, iv] = SPVAL
                    gs_mol_sha_ln[p, iv] = SPVAL
                    vegwp_ln[p, :] .= SPVAL
                end

                # Final cs and ci
                cs_sun = cair[p] - 1.4 / gb_mol_arr[p] * an_sun[p, iv] * forc_pbot[c]
                cs_sun = max(cs_sun, MAX_CS)
                ci_z_sun[p, iv] = cair[p] - an_sun[p, iv] * forc_pbot[c] *
                    (1.4 * gs_mol_sun[p, iv] + 1.6 * gb_mol_arr[p]) /
                    (gb_mol_arr[p] * gs_mol_sun[p, iv])
                ci_z_sun[p, iv] = max(ci_z_sun[p, iv], 1.0e-06)

                cs_sha = cair[p] - 1.4 / gb_mol_arr[p] * an_sha[p, iv] * forc_pbot[c]
                cs_sha = max(cs_sha, MAX_CS)
                ci_z_sha[p, iv] = cair[p] - an_sha[p, iv] * forc_pbot[c] *
                    (1.4 * gs_mol_sha[p, iv] + 1.6 * gb_mol_arr[p]) /
                    (gb_mol_arr[p] * gs_mol_sha[p, iv])
                ci_z_sha[p, iv] = max(ci_z_sha[p, iv], 1.0e-06)

                # Convert to resistance
                gs = gs_mol_sun[p, iv] / cf
                rs_z_sun[p, iv] = min(1.0 / gs, rsmax0)
                rs_z_sun[p, iv] = rs_z_sun[p, iv] / o3coefg_sun[p]
                gs = gs_mol_sha[p, iv] / cf
                rs_z_sha[p, iv] = min(1.0 / gs, rsmax0)
                rs_z_sha[p, iv] = rs_z_sha[p, iv] / o3coefg_sha[p]

                # Photosynthesis output
                psn_z_sun[p, iv] = ag[p, SUN, iv] * o3coefv_sun[p]
                psn_wc_z_sun[p, iv] = 0.0
                psn_wj_z_sun[p, iv] = 0.0
                psn_wp_z_sun[p, iv] = 0.0
                if ac[p, SUN, iv] <= aj[p, SUN, iv] && ac[p, SUN, iv] <= ap[p, SUN, iv]
                    psn_wc_z_sun[p, iv] = psn_z_sun[p, iv]
                elseif aj[p, SUN, iv] < ac[p, SUN, iv] && aj[p, SUN, iv] <= ap[p, SUN, iv]
                    psn_wj_z_sun[p, iv] = psn_z_sun[p, iv]
                elseif ap[p, SUN, iv] < ac[p, SUN, iv] && ap[p, SUN, iv] < aj[p, SUN, iv]
                    psn_wp_z_sun[p, iv] = psn_z_sun[p, iv]
                end

                psn_z_sha[p, iv] = ag[p, SHA, iv] * o3coefv_sha[p]
                psn_wc_z_sha[p, iv] = 0.0
                psn_wj_z_sha[p, iv] = 0.0
                psn_wp_z_sha[p, iv] = 0.0
                if ac[p, SHA, iv] <= aj[p, SHA, iv] && ac[p, SHA, iv] <= ap[p, SHA, iv]
                    psn_wc_z_sha[p, iv] = psn_z_sha[p, iv]
                elseif aj[p, SHA, iv] < ac[p, SHA, iv] && aj[p, SHA, iv] <= ap[p, SHA, iv]
                    psn_wj_z_sha[p, iv] = psn_z_sha[p, iv]
                elseif ap[p, SHA, iv] < ac[p, SHA, iv] && ap[p, SHA, iv] < aj[p, SHA, iv]
                    psn_wp_z_sha[p, iv] = psn_z_sha[p, iv]
                end

                # Relative humidity at leaf surface
                hs = (gb_mol_arr[p] * ceair + gs_mol_sun[p, iv] * esat_tv[p]) /
                     ((gb_mol_arr[p] + gs_mol_sun[p, iv]) * esat_tv[p])
                rh_leaf_sun[p] = hs
                hs = (gb_mol_arr[p] * ceair + gs_mol_sha[p, iv] * esat_tv[p]) /
                     ((gb_mol_arr[p] + gs_mol_sha[p, iv]) * esat_tv[p])
                rh_leaf_sha[p] = hs
            end
        end
    end

    # ---- Pass 4: Canopy integration ----
    for p in bounds_patch
        mask_patch[p] || continue

        # Sunlit canopy
        psncan_sun = 0.0; psncan_wc_sun = 0.0; psncan_wj_sun = 0.0; psncan_wp_sun = 0.0
        lmrcan_sun = 0.0; gscan_sun = 0.0; laican_sun = 0.0
        for iv in 1:nrad[p]
            psncan_sun += psn_z_sun[p, iv] * lai_z_sun_in[p, iv]
            psncan_wc_sun += psn_wc_z_sun[p, iv] * lai_z_sun_in[p, iv]
            psncan_wj_sun += psn_wj_z_sun[p, iv] * lai_z_sun_in[p, iv]
            psncan_wp_sun += psn_wp_z_sun[p, iv] * lai_z_sun_in[p, iv]
            if round(Int, crop_pft[ivt[p]]) == 0 && modifyphoto_and_lmr_forcrop
                lmrcan_sun += lmr_z_sun[p, iv] * lai_z_sun_in[p, iv] * bsun_arr[p]
            else
                lmrcan_sun += lmr_z_sun[p, iv] * lai_z_sun_in[p, iv]
            end
            gscan_sun += lai_z_sun_in[p, iv] / (rb[p] + rs_z_sun[p, iv])
            laican_sun += lai_z_sun_in[p, iv]
        end
        if laican_sun > 0.0
            psn_sun[p] = psncan_sun / laican_sun
            psn_wc_sun[p] = psncan_wc_sun / laican_sun
            psn_wj_sun[p] = psncan_wj_sun / laican_sun
            psn_wp_sun[p] = psncan_wp_sun / laican_sun
            lmr_sun[p] = lmrcan_sun / laican_sun
            rs_sun[p] = laican_sun / gscan_sun - rb[p]
        else
            psn_sun[p] = 0.0; psn_wc_sun[p] = 0.0; psn_wj_sun[p] = 0.0; psn_wp_sun[p] = 0.0
            lmr_sun[p] = 0.0; rs_sun[p] = 0.0
        end

        # Shaded canopy
        psncan_sha = 0.0; psncan_wc_sha = 0.0; psncan_wj_sha = 0.0; psncan_wp_sha = 0.0
        lmrcan_sha = 0.0; gscan_sha = 0.0; laican_sha = 0.0
        for iv in 1:nrad[p]
            psncan_sha += psn_z_sha[p, iv] * lai_z_sha_in[p, iv]
            psncan_wc_sha += psn_wc_z_sha[p, iv] * lai_z_sha_in[p, iv]
            psncan_wj_sha += psn_wj_z_sha[p, iv] * lai_z_sha_in[p, iv]
            psncan_wp_sha += psn_wp_z_sha[p, iv] * lai_z_sha_in[p, iv]
            if round(Int, crop_pft[ivt[p]]) == 0 && modifyphoto_and_lmr_forcrop
                lmrcan_sha += lmr_z_sha[p, iv] * lai_z_sha_in[p, iv] * bsha_arr[p]
            else
                lmrcan_sha += lmr_z_sha[p, iv] * lai_z_sha_in[p, iv]
            end
            gscan_sha += lai_z_sha_in[p, iv] / (rb[p] + rs_z_sha[p, iv])
            laican_sha += lai_z_sha_in[p, iv]
        end
        if laican_sha > 0.0
            psn_sha[p] = psncan_sha / laican_sha
            psn_wc_sha[p] = psncan_wc_sha / laican_sha
            psn_wj_sha[p] = psncan_wj_sha / laican_sha
            psn_wp_sha[p] = psncan_wp_sha / laican_sha
            lmr_sha[p] = lmrcan_sha / laican_sha
            rs_sha[p] = laican_sha / gscan_sha - rb[p]
        else
            psn_sha[p] = 0.0; psn_wc_sha[p] = 0.0; psn_wj_sha[p] = 0.0; psn_wp_sha[p] = 0.0
            lmr_sha[p] = 0.0; rs_sha[p] = 0.0
        end

        # btran from LAI-weighted bsun/bsha
        if laican_sha + laican_sun > 0.0
            btran[p] = bsun_arr[p] * (laican_sun / (laican_sun + laican_sha)) +
                       bsha_arr[p] * (laican_sha / (laican_sun + laican_sha))
        else
            btran[p] = bsun_arr[p]
        end
    end

    return nothing
end
