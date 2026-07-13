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
Base.@kwdef mutable struct PhotoParamsData{FT<:Real}
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
    krmax::Vector{FT} = Float64[]                  # (0:mxpft)
    kmax::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # (0:mxpft, nvegwcs)
    psi50::Matrix{FT} = Matrix{Float64}(undef, 0, 0) # (0:mxpft, nvegwcs)
    ck::Matrix{FT} = Matrix{Float64}(undef, 0, 0)    # (0:mxpft, nvegwcs)
    lmr_intercept_atkin::Vector{FT} = Float64[]    # (0:mxpft)
    theta_cj::Vector{FT} = Float64[]               # (0:mxpft)
end

# Global params instance
const params_inst = PhotoParamsData()

"""
    photo_params_init!(pp::PhotoParamsData, mxpft::Int, nvegwcs::Int)

Allocate photosynthesis parameter arrays.
"""
function photo_params_init!(pp::PhotoParamsData{FT}, mxpft::Int=MXPFT, nvegwcs::Int=NVEGWCS) where {FT}
    pp.krmax = fill(FT(NaN), mxpft + 1)
    pp.theta_cj = fill(FT(NaN), mxpft + 1)
    pp.kmax = fill(FT(NaN), mxpft + 1, nvegwcs)
    pp.psi50 = fill(FT(NaN), mxpft + 1, nvegwcs)
    pp.ck = fill(FT(NaN), mxpft + 1, nvegwcs)
    pp.lmr_intercept_atkin = fill(FT(NaN), mxpft + 1)
    return nothing
end

"""
    photo_params_clean!(pp::PhotoParamsData)

Deallocate photosynthesis parameter arrays.
"""
function photo_params_clean!(pp::PhotoParamsData{FT}) where {FT}
    pp.krmax = FT[]
    pp.theta_cj = FT[]
    pp.kmax = Matrix{FT}(undef, 0, 0)
    pp.psi50 = Matrix{FT}(undef, 0, 0)
    pp.ck = Matrix{FT}(undef, 0, 0)
    pp.lmr_intercept_atkin = FT[]
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
Base.@kwdef mutable struct PhotosynthesisData{FT<:Real,
                                  V<:AbstractVector{FT},
                                  M<:AbstractMatrix{FT},
                                  A3<:AbstractArray{FT,3},
                                  VB<:AbstractVector{Bool}}
    # Logical/config
    c3flag_patch::VB = Bool[]

    # PHS-specific 3D arrays (np, 2, nlevcan)
    ac_phs_patch::A3 = Array{Float64}(undef, 0, 0, 0)
    aj_phs_patch::A3 = Array{Float64}(undef, 0, 0, 0)
    ap_phs_patch::A3 = Array{Float64}(undef, 0, 0, 0)
    ag_phs_patch::A3 = Array{Float64}(undef, 0, 0, 0)
    vcmax_z_phs_patch::A3 = Array{Float64}(undef, 0, 0, 0)
    tpu_z_phs_patch::A3 = Array{Float64}(undef, 0, 0, 0)
    kp_z_phs_patch::A3 = Array{Float64}(undef, 0, 0, 0)

    # Sunlit/shaded net photosynthesis (np, nlevcan)
    an_sun_patch::M = Matrix{Float64}(undef, 0, 0)
    an_sha_patch::M = Matrix{Float64}(undef, 0, 0)

    # Stomatal conductance sunlit/shaded (np, nlevcan)
    gs_mol_sun_patch::M = Matrix{Float64}(undef, 0, 0)
    gs_mol_sha_patch::M = Matrix{Float64}(undef, 0, 0)
    gs_mol_sun_ln_patch::M = Matrix{Float64}(undef, 0, 0)
    gs_mol_sha_ln_patch::M = Matrix{Float64}(undef, 0, 0)

    # Standard (non-PHS) 2D arrays (np, nlevcan)
    ac_patch::M = Matrix{Float64}(undef, 0, 0)
    aj_patch::M = Matrix{Float64}(undef, 0, 0)
    ap_patch::M = Matrix{Float64}(undef, 0, 0)
    ag_patch::M = Matrix{Float64}(undef, 0, 0)
    an_patch::M = Matrix{Float64}(undef, 0, 0)
    vcmax_z_patch::M = Matrix{Float64}(undef, 0, 0)
    tpu_z_patch::M = Matrix{Float64}(undef, 0, 0)
    kp_z_patch::M = Matrix{Float64}(undef, 0, 0)
    gs_mol_patch::M = Matrix{Float64}(undef, 0, 0)

    # 1D patch-level
    cp_patch::V = Float64[]
    kc_patch::V = Float64[]
    ko_patch::V = Float64[]
    qe_patch::V = Float64[]
    bbb_patch::V = Float64[]
    mbb_patch::V = Float64[]
    gb_mol_patch::V = Float64[]
    rh_leaf_patch::V = Float64[]
    vpd_can_patch::V = Float64[]

    # Photosynthesis outputs
    psnsun_patch::V = Float64[]
    psnsha_patch::V = Float64[]
    c13_psnsun_patch::V = Float64[]
    c13_psnsha_patch::V = Float64[]
    c14_psnsun_patch::V = Float64[]
    c14_psnsha_patch::V = Float64[]

    psnsun_z_patch::M = Matrix{Float64}(undef, 0, 0)
    psnsha_z_patch::M = Matrix{Float64}(undef, 0, 0)
    psnsun_wc_patch::V = Float64[]
    psnsha_wc_patch::V = Float64[]
    psnsun_wj_patch::V = Float64[]
    psnsha_wj_patch::V = Float64[]
    psnsun_wp_patch::V = Float64[]
    psnsha_wp_patch::V = Float64[]

    fpsn_patch::V = Float64[]
    fpsn_wc_patch::V = Float64[]
    fpsn_wj_patch::V = Float64[]
    fpsn_wp_patch::V = Float64[]

    lnca_patch::V = Float64[]

    lmrsun_patch::V = Float64[]
    lmrsha_patch::V = Float64[]
    lmrsun_z_patch::M = Matrix{Float64}(undef, 0, 0)
    lmrsha_z_patch::M = Matrix{Float64}(undef, 0, 0)

    alphapsnsun_patch::V = Float64[]
    alphapsnsha_patch::V = Float64[]
    rc13_canair_patch::V = Float64[]
    rc13_psnsun_patch::V = Float64[]
    rc13_psnsha_patch::V = Float64[]

    cisun_z_patch::M = Matrix{Float64}(undef, 0, 0)
    cisha_z_patch::M = Matrix{Float64}(undef, 0, 0)

    rssun_z_patch::M = Matrix{Float64}(undef, 0, 0)
    rssha_z_patch::M = Matrix{Float64}(undef, 0, 0)
    rssun_patch::V = Float64[]
    rssha_patch::V = Float64[]

    luvcmax25top_patch::V = Float64[]
    lujmax25top_patch::V = Float64[]
    lutpu25top_patch::V = Float64[]

    # LUNA-specific
    vcmx25_z_patch::M = Matrix{Float64}(undef, 0, 0)
    jmx25_z_patch::M = Matrix{Float64}(undef, 0, 0)
    vcmx25_z_last_valid_patch::M = Matrix{Float64}(undef, 0, 0)
    jmx25_z_last_valid_patch::M = Matrix{Float64}(undef, 0, 0)
    pnlc_z_patch::M = Matrix{Float64}(undef, 0, 0)
    enzs_z_patch::M = Matrix{Float64}(undef, 0, 0)
    fpsn24_patch::V = Float64[]

    # Configuration switches
    rootstem_acc::Bool = false
    light_inhibit::Bool = false
    leafresp_method::Int = LEAFRESP_MTD_RYAN1991
    stomatalcond_mtd::Int = STOMATALCOND_MTD_MEDLYN2011
    modifyphoto_and_lmr_forcrop::Bool = false
end

PhotosynthesisData{FT}(; kwargs...) where {FT<:Real} =
    PhotosynthesisData{FT, Vector{FT}, Matrix{FT}, Array{FT,3}, Vector{Bool}}(; kwargs...)
Adapt.@adapt_structure PhotosynthesisData


"""
    photosynthesis_data_init!(ps, np::Int; nlevcan::Int=NLEVCAN, use_luna::Bool=false)

Allocate and initialize all fields for `np` patches.
"""
function photosynthesis_data_init!(ps::PhotosynthesisData{FT}, np::Int;
                                   nlevcan::Int=NLEVCAN, use_luna::Bool=false) where {FT}
    _nan = FT(NaN)

    ps.c3flag_patch = fill(false, np)

    ps.ac_phs_patch = fill(_nan, np, 2, nlevcan)
    ps.aj_phs_patch = fill(_nan, np, 2, nlevcan)
    ps.ap_phs_patch = fill(_nan, np, 2, nlevcan)
    ps.ag_phs_patch = fill(_nan, np, 2, nlevcan)
    ps.vcmax_z_phs_patch = fill(_nan, np, 2, nlevcan)
    ps.tpu_z_phs_patch = fill(_nan, np, 2, nlevcan)
    ps.kp_z_phs_patch = fill(_nan, np, 2, nlevcan)

    ps.an_sun_patch = fill(_nan, np, nlevcan)
    ps.an_sha_patch = fill(_nan, np, nlevcan)

    ps.gs_mol_sun_patch = fill(_nan, np, nlevcan)
    ps.gs_mol_sha_patch = fill(_nan, np, nlevcan)
    ps.gs_mol_sun_ln_patch = fill(_nan, np, nlevcan)
    ps.gs_mol_sha_ln_patch = fill(_nan, np, nlevcan)

    ps.ac_patch = fill(_nan, np, nlevcan)
    ps.aj_patch = fill(_nan, np, nlevcan)
    ps.ap_patch = fill(_nan, np, nlevcan)
    ps.ag_patch = fill(_nan, np, nlevcan)
    ps.an_patch = fill(_nan, np, nlevcan)
    ps.vcmax_z_patch = fill(_nan, np, nlevcan)
    ps.tpu_z_patch = fill(_nan, np, nlevcan)
    ps.kp_z_patch = fill(_nan, np, nlevcan)
    ps.gs_mol_patch = fill(_nan, np, nlevcan)

    ps.cp_patch = fill(_nan, np)
    ps.kc_patch = fill(_nan, np)
    ps.ko_patch = fill(_nan, np)
    ps.qe_patch = fill(_nan, np)
    ps.bbb_patch = fill(_nan, np)
    ps.mbb_patch = fill(_nan, np)
    ps.gb_mol_patch = fill(_nan, np)
    ps.rh_leaf_patch = fill(_nan, np)
    ps.vpd_can_patch = fill(_nan, np)

    ps.psnsun_patch = fill(_nan, np)
    ps.psnsha_patch = fill(_nan, np)
    ps.c13_psnsun_patch = fill(_nan, np)
    ps.c13_psnsha_patch = fill(_nan, np)
    ps.c14_psnsun_patch = fill(_nan, np)
    ps.c14_psnsha_patch = fill(_nan, np)

    ps.psnsun_z_patch = fill(_nan, np, nlevcan)
    ps.psnsha_z_patch = fill(_nan, np, nlevcan)
    ps.psnsun_wc_patch = fill(_nan, np)
    ps.psnsha_wc_patch = fill(_nan, np)
    ps.psnsun_wj_patch = fill(_nan, np)
    ps.psnsha_wj_patch = fill(_nan, np)
    ps.psnsun_wp_patch = fill(_nan, np)
    ps.psnsha_wp_patch = fill(_nan, np)

    ps.fpsn_patch = fill(_nan, np)
    ps.fpsn_wc_patch = fill(_nan, np)
    ps.fpsn_wj_patch = fill(_nan, np)
    ps.fpsn_wp_patch = fill(_nan, np)

    ps.lnca_patch = fill(_nan, np)

    ps.lmrsun_z_patch = fill(_nan, np, nlevcan)
    ps.lmrsha_z_patch = fill(_nan, np, nlevcan)
    ps.lmrsun_patch = fill(_nan, np)
    ps.lmrsha_patch = fill(_nan, np)

    ps.alphapsnsun_patch = fill(_nan, np)
    ps.alphapsnsha_patch = fill(_nan, np)
    ps.rc13_canair_patch = fill(_nan, np)
    ps.rc13_psnsun_patch = fill(_nan, np)
    ps.rc13_psnsha_patch = fill(_nan, np)

    ps.cisun_z_patch = fill(_nan, np, nlevcan)
    ps.cisha_z_patch = fill(_nan, np, nlevcan)

    ps.rssun_z_patch = fill(FT(2.0e4), np, nlevcan)   # large default (closed stomata)
    ps.rssha_z_patch = fill(FT(2.0e4), np, nlevcan)
    ps.rssun_patch = fill(FT(2.0e4), np)              # large default (closed stomata)
    ps.rssha_patch = fill(FT(2.0e4), np)

    ps.luvcmax25top_patch = fill(_nan, np)
    ps.lujmax25top_patch = fill(_nan, np)
    ps.lutpu25top_patch = fill(_nan, np)

    if use_luna
        ps.vcmx25_z_patch = fill(FT(30.0), np, nlevcan)
        ps.jmx25_z_patch = fill(FT(60.0), np, nlevcan)
        ps.vcmx25_z_last_valid_patch = fill(FT(30.0), np, nlevcan)
        ps.jmx25_z_last_valid_patch = fill(FT(60.0), np, nlevcan)
        ps.pnlc_z_patch = fill(FT(0.01), np, nlevcan)
        ps.fpsn24_patch = fill(_nan, np)
        ps.enzs_z_patch = fill(one(FT), np, nlevcan)
    end

    return nothing
end

"""
    photosyns_init_cold!(ps, patch_data, lun, bounds_patch)

Cold-start initialization for the photosynthesis state.

Ported from `photosyns_type%InitCold` in `PhotosynthesisMod.F90`. Fortran's
InitCold is small and does exactly two things:

  * `alphapsnsun/alphapsnsha` → `spval` (an isotope diagnostic that is only
    defined where the C13 discrimination solve runs; spval is the deliberate
    missing-value flag, NOT zero);
  * `psnsun/psnsha` (and the C13/C14 twins) → `0.0` on SPECIAL landunits
    (urban / glacier / lake / wetland). Those patches never enter the
    photosynthesis solve, so nothing else ever writes their photosynthesis
    rate — yet `psnsun/psnsha` are summed into the patch→column→gridcell GPP
    aggregation. Left at the allocator's NaN they poison it.

This routine had NO Julia port (unlike the ~20 `*_init_cold!` that were ported
and then never called); it is added here so the InitCold class is complete.
"""
function photosyns_init_cold!(ps::PhotosynthesisData,
                              patch_data::PatchData,
                              lun::LandunitData,
                              bounds_patch::UnitRange{Int})
    has_c13 = !isempty(ps.c13_psnsun_patch)
    has_c14 = !isempty(ps.c14_psnsun_patch)
    for p in bounds_patch
        l = patch_data.landunit[p]

        ps.alphapsnsun_patch[p] = SPVAL
        ps.alphapsnsha_patch[p] = SPVAL

        (l >= 1 && l <= length(lun.ifspecial) && lun.ifspecial[l]) || continue
        ps.psnsun_patch[p] = 0.0
        ps.psnsha_patch[p] = 0.0
        if has_c13
            ps.c13_psnsun_patch[p] = 0.0
            ps.c13_psnsha_patch[p] = 0.0
        end
        if has_c14
            ps.c14_psnsun_patch[p] = 0.0
            ps.c14_psnsha_patch[p] = 0.0
        end
    end
    return nothing
end

"""
    photosynthesis_data_clean!(ps::PhotosynthesisData)

Reset all fields to empty arrays.
"""
function photosynthesis_data_clean!(ps::PhotosynthesisData{FT}) where {FT}
    ps.c3flag_patch = Bool[]
    ps.ac_phs_patch = Array{Float64}(undef, 0, 0, 0)
    ps.aj_phs_patch = Array{Float64}(undef, 0, 0, 0)
    ps.ap_phs_patch = Array{Float64}(undef, 0, 0, 0)
    ps.ag_phs_patch = Array{Float64}(undef, 0, 0, 0)
    ps.vcmax_z_phs_patch = Array{Float64}(undef, 0, 0, 0)
    ps.tpu_z_phs_patch = Array{Float64}(undef, 0, 0, 0)
    ps.kp_z_phs_patch = Array{Float64}(undef, 0, 0, 0)
    ps.an_sun_patch = Matrix{FT}(undef, 0, 0)
    ps.an_sha_patch = Matrix{FT}(undef, 0, 0)
    ps.gs_mol_sun_patch = Matrix{FT}(undef, 0, 0)
    ps.gs_mol_sha_patch = Matrix{FT}(undef, 0, 0)
    ps.gs_mol_sun_ln_patch = Matrix{FT}(undef, 0, 0)
    ps.gs_mol_sha_ln_patch = Matrix{FT}(undef, 0, 0)
    ps.ac_patch = Matrix{FT}(undef, 0, 0)
    ps.aj_patch = Matrix{FT}(undef, 0, 0)
    ps.ap_patch = Matrix{FT}(undef, 0, 0)
    ps.ag_patch = Matrix{FT}(undef, 0, 0)
    ps.an_patch = Matrix{FT}(undef, 0, 0)
    ps.vcmax_z_patch = Matrix{FT}(undef, 0, 0)
    ps.tpu_z_patch = Matrix{FT}(undef, 0, 0)
    ps.kp_z_patch = Matrix{FT}(undef, 0, 0)
    ps.gs_mol_patch = Matrix{FT}(undef, 0, 0)
    ps.cp_patch = FT[]
    ps.kc_patch = FT[]
    ps.ko_patch = FT[]
    ps.qe_patch = FT[]
    ps.bbb_patch = FT[]
    ps.mbb_patch = FT[]
    ps.gb_mol_patch = FT[]
    ps.rh_leaf_patch = FT[]
    ps.vpd_can_patch = FT[]
    ps.psnsun_patch = FT[]
    ps.psnsha_patch = FT[]
    ps.c13_psnsun_patch = FT[]
    ps.c13_psnsha_patch = FT[]
    ps.c14_psnsun_patch = FT[]
    ps.c14_psnsha_patch = FT[]
    ps.psnsun_z_patch = Matrix{FT}(undef, 0, 0)
    ps.psnsha_z_patch = Matrix{FT}(undef, 0, 0)
    ps.psnsun_wc_patch = FT[]
    ps.psnsha_wc_patch = FT[]
    ps.psnsun_wj_patch = FT[]
    ps.psnsha_wj_patch = FT[]
    ps.psnsun_wp_patch = FT[]
    ps.psnsha_wp_patch = FT[]
    ps.fpsn_patch = FT[]
    ps.fpsn_wc_patch = FT[]
    ps.fpsn_wj_patch = FT[]
    ps.fpsn_wp_patch = FT[]
    ps.lnca_patch = FT[]
    ps.lmrsun_z_patch = Matrix{FT}(undef, 0, 0)
    ps.lmrsha_z_patch = Matrix{FT}(undef, 0, 0)
    ps.lmrsun_patch = FT[]
    ps.lmrsha_patch = FT[]
    ps.alphapsnsun_patch = FT[]
    ps.alphapsnsha_patch = FT[]
    ps.rc13_canair_patch = FT[]
    ps.rc13_psnsun_patch = FT[]
    ps.rc13_psnsha_patch = FT[]
    ps.cisun_z_patch = Matrix{FT}(undef, 0, 0)
    ps.cisha_z_patch = Matrix{FT}(undef, 0, 0)
    ps.rssun_z_patch = Matrix{FT}(undef, 0, 0)
    ps.rssha_z_patch = Matrix{FT}(undef, 0, 0)
    ps.rssun_patch = FT[]
    ps.rssha_patch = FT[]
    ps.luvcmax25top_patch = FT[]
    ps.lujmax25top_patch = FT[]
    ps.lutpu25top_patch = FT[]
    ps.vcmx25_z_patch = Matrix{FT}(undef, 0, 0)
    ps.jmx25_z_patch = Matrix{FT}(undef, 0, 0)
    ps.vcmx25_z_last_valid_patch = Matrix{FT}(undef, 0, 0)
    ps.jmx25_z_last_valid_patch = Matrix{FT}(undef, 0, 0)
    ps.pnlc_z_patch = Matrix{FT}(undef, 0, 0)
    ps.enzs_z_patch = Matrix{FT}(undef, 0, 0)
    ps.fpsn24_patch = FT[]
    return nothing
end

# =====================================================================
# Helper functions: temperature response, quadratic solver
# =====================================================================

# Explicit molar-units alias for the photosynthesis equations. RGAS is built from
# the same CESM shared constants and multiply order as CTSM.
const RGAS_PSN = RGAS

"""
    ft_photo(tl, ha)

Photosynthesis temperature response function.
"""
@inline function ft_photo(tl::Real, ha::Real)
    # eltype-generic (T from tl) so it lowers to valid Metal IR under Float32;
    # byte-identical to the Float64 literals on the CPU path.
    T = typeof(tl)
    return exp(T(ha) / (T(RGAS_PSN) * (T(TFRZ) + T(25.0))) * (one(T) - (T(TFRZ) + T(25.0)) / tl))
end

"""
    fth_photo(tl, hd, se, scaleFactor)

Photosynthesis temperature inhibition function.
"""
@inline function fth_photo(tl::Real, hd::Real, se::Real, scaleFactor::Real)
    T = typeof(tl)
    return T(scaleFactor) / (one(T) + exp((-T(hd) + T(se) * tl) / (T(RGAS_PSN) * tl)))
end

"""
    fth25_photo(hd, se)

Scaling factor for photosynthesis temperature inhibition at 25°C.
"""
@inline function fth25_photo(hd::Real, se::Real)
    T = typeof(se)
    return one(T) + exp((-T(hd) + T(se) * (T(TFRZ) + T(25.0))) / (T(RGAS_PSN) * (T(TFRZ) + T(25.0))))
end

"""
    quadratic_solve(a, b, c) -> (r1, r2)

Solve the quadratic equation `a*x^2 + b*x + c = 0` for real roots.
Returns the two roots `(r1, r2)` where `r1 >= r2`.
"""
function quadratic_solve(a::Real, b::Real, c::Real)
    # eltype-generic (promote a/b/c) so it lowers to valid Metal IR under Float32;
    # byte-identical to the Float64 literals on the CPU path.
    T = promote_type(typeof(a), typeof(b), typeof(c))
    if a == zero(a)
        if b == zero(b)
            return (zero(T), zero(T))
        else
            r1 = -c / b
            return (r1, r1)
        end
    end
    discriminant = b * b - T(4.0) * a * c
    if discriminant < zero(T)
        discriminant = zero(T)
    end
    q = -T(0.5) * (b + sign(b) * sqrt(discriminant))
    if q == zero(q)
        r1 = zero(T)
        r2 = zero(T)
    else
        r1 = q / a
        r2 = c / q
    end
    if r1 < r2
        r1, r2 = r2, r1
    end
    # CONTRACT: the roots are returned SORTED, r1 >= r2, ALWAYS.
    #
    # Callers must therefore use `r1` / `r2` directly to select the larger / smaller root. Do NOT
    # wrap them in smooth_max(r1, r2) / smooth_min(r1, r2): the selection is not a physical branch,
    # it is an already-resolved sort, so the smoothing recovers no derivative information — while
    # smooth_max(a,b) OVERSHOOTS max(a,b) by up to log(2)/k IN THE UNITS OF THE ROOTS. Those units
    # are not O(1): for the Medlyn stomatal-conductance quadratic the roots are a conductance in
    # mol H2O/m2/s (~0.001-0.5), so the generic k = 50 imposed a
    #     log(2)/50 = 0.0139 mol/m2/s = 13 900 umol/m2/s
    # FLOOR on gs_mol — 139x the medlynintercept of 100 umol/m2/s, i.e. stomata that can never
    # close. (Measured directly: one smoothed step moved gs_mol from 100 to 1.39e4.) The same
    # guards on the co-limitation quadratics fabricated ~0.008-0.024 umol/m2/s of photosynthesis on
    # patches whose true assimilation is exactly zero.
    return (r1, r2)
end

"""
    plc(x, ivt, level, plc_method, params)

Return value of vulnerability curve at water potential `x`.
`ivt` is 1-based PFT index, `level` is segment index (SUN, SHA, XYL, ROOT_SEG).
"""
function plc(x::Real, ivt::Int, level::Int, plc_method::Int, params::PhotoParamsData=params_inst)
    if plc_method == VEGETATION_WEIBULL
        # Water potential is always <= 0 physically; for x >= 0 (fully hydrated, no
        # cavitation) the curve is 1. Guarding here avoids (x/psi50)<0 raised to the
        # non-integer ck -> NaN, which only arises in the degenerate all-frozen PHS
        # solve (vegwp driven >= 0); Fortran never reaches that regime, so parity holds.
        x >= 0.0 && return 1.0
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
function d1plc(x::Real, ivt::Int, level::Int, plc_method::Int, params::PhotoParamsData=params_inst)
    if plc_method == VEGETATION_WEIBULL
        # plc is constant (1) for x >= 0, so its derivative is 0 there; this also
        # avoids the /x singularity and (x/psi50)<0 ^ ck NaN (see plc above).
        x >= 0.0 && return 0.0
        ck_val = params.ck[ivt, level]
        psi50_val = params.psi50[ivt, level]
        val = -ck_val * log(2.0) * (2.0^(-(x / psi50_val)^ck_val)) *
              ((x / psi50_val)^ck_val) / x
        return val
    else
        error("D1PLC: must choose valid plc_method")
    end
end

# GPU-callable positional cores for the vulnerability curve — the WEIBULL branch
# only (the always-valid device path; the host plc/d1plc keep the method dispatch +
# error() for back-compat/tests). psi50_v/ck_v are the per-(ivt,level) curve params
# the caller pulls from a threaded PsnPhsParams bundle. Byte-identical to plc/d1plc
# on Float64 (T(2.0)===2.0); eltype-generic so Float32/Metal carries no Float64.
@inline function _plc(x::Real, psi50_v::Real, ck_v::Real)
    T = typeof(x)
    x >= zero(T) && return one(T)   # no cavitation for x>=0; avoids (neg)^ck NaN (see plc)
    val = T(2.0)^(-(x / psi50_v)^ck_v)
    return val < T(0.005) ? zero(T) : val
end

@inline function _d1plc(x::Real, psi50_v::Real, ck_v::Real)
    T = typeof(x)
    x >= zero(T) && return zero(T)  # plc flat for x>=0; avoids /x singularity + NaN
    return -ck_v * log(T(2.0)) * (T(2.0)^(-(x / psi50_v)^ck_v)) *
           ((x / psi50_v)^ck_v) / x
end

# =====================================================================
# ci_func! — evaluate f(ci) for the standard (non-PHS) method
# =====================================================================

# Immutable device-view of the PhotosynthesisData fields the GPU kernels + solver
# cores touch. PhotosynthesisData is a MUTABLE struct → non-bitstype → cannot be
# passed by value to a Metal kernel; this immutable bundle (built per launch from
# ps via _psn_dv) can. Field NAMES mirror PhotosynthesisData exactly, so the
# kernel/core bodies (ps.ac_patch, …) are unchanged — `ps` just becomes a PsnDV
# on the device path (and the real mutable struct on the host/back-compat path,
# since the cores take `ps` untyped). {Vb=Bool vec, V=float vec, M=float matrix}.
# PHS fields get their OWN type params (M3 3D, M2 2D, Vp 1D): in the non-PHS AD path
# these are bundled but NOT written, so they stay Float64 while the active V/M fields
# become ForwardDiff.Dual — a single shared param can't unify (the canopy CfEn lesson).
Base.@kwdef struct PsnDV{Vb,V,M,M3,M2,Vp}
    c3flag_patch::Vb
    qe_patch::V; bbb_patch::V; mbb_patch::V; kc_patch::V; ko_patch::V; cp_patch::V
    lnca_patch::V; gb_mol_patch::V; rh_leaf_patch::V; vpd_can_patch::V
    psnsun_patch::V; psnsha_patch::V
    psnsun_wc_patch::V; psnsun_wj_patch::V; psnsun_wp_patch::V
    psnsha_wc_patch::V; psnsha_wj_patch::V; psnsha_wp_patch::V
    rssun_patch::V; rssha_patch::V; lmrsun_patch::V; lmrsha_patch::V
    vcmax_z_patch::M; tpu_z_patch::M; kp_z_patch::M
    ac_patch::M; aj_patch::M; ap_patch::M; ag_patch::M; an_patch::M
    gs_mol_patch::M; gs_mol_sun_patch::M; gs_mol_sha_patch::M
    lmrsun_z_patch::M; lmrsha_z_patch::M; psnsun_z_patch::M; psnsha_z_patch::M
    rssun_z_patch::M; rssha_z_patch::M; cisun_z_patch::M; cisha_z_patch::M
    # --- PHS (plant-hydraulic-stress) fields ---
    ac_phs_patch::M3; aj_phs_patch::M3; ap_phs_patch::M3; ag_phs_patch::M3
    vcmax_z_phs_patch::M3; tpu_z_phs_patch::M3; kp_z_phs_patch::M3
    vcmx25_z_patch::M; jmx25_z_patch::M   # LUNA-acclimated vcmax25/jmax25 (injected)
    an_sun_patch::M2; an_sha_patch::M2; gs_mol_sun_ln_patch::M2; gs_mol_sha_ln_patch::M2
    luvcmax25top_patch::Vp; lujmax25top_patch::Vp; lutpu25top_patch::Vp
    alphapsnsun_patch::Vp; alphapsnsha_patch::Vp
    stomatalcond_mtd::Int
end
Adapt.@adapt_structure PsnDV   # so KA adapts MtlArray fields -> MtlDeviceArray (bitstype) at launch

# Build a PsnDV from ps (array fields are shared refs — writes flow back to ps).
_psn_dv(ps) = PsnDV(; c3flag_patch = ps.c3flag_patch,
    qe_patch = ps.qe_patch, bbb_patch = ps.bbb_patch, mbb_patch = ps.mbb_patch,
    kc_patch = ps.kc_patch, ko_patch = ps.ko_patch, cp_patch = ps.cp_patch,
    lnca_patch = ps.lnca_patch, gb_mol_patch = ps.gb_mol_patch,
    rh_leaf_patch = ps.rh_leaf_patch, vpd_can_patch = ps.vpd_can_patch,
    psnsun_patch = ps.psnsun_patch, psnsha_patch = ps.psnsha_patch,
    psnsun_wc_patch = ps.psnsun_wc_patch, psnsun_wj_patch = ps.psnsun_wj_patch,
    psnsun_wp_patch = ps.psnsun_wp_patch, psnsha_wc_patch = ps.psnsha_wc_patch,
    psnsha_wj_patch = ps.psnsha_wj_patch, psnsha_wp_patch = ps.psnsha_wp_patch,
    rssun_patch = ps.rssun_patch, rssha_patch = ps.rssha_patch,
    lmrsun_patch = ps.lmrsun_patch, lmrsha_patch = ps.lmrsha_patch,
    vcmax_z_patch = ps.vcmax_z_patch, tpu_z_patch = ps.tpu_z_patch,
    kp_z_patch = ps.kp_z_patch, ac_patch = ps.ac_patch, aj_patch = ps.aj_patch,
    ap_patch = ps.ap_patch, ag_patch = ps.ag_patch, an_patch = ps.an_patch,
    gs_mol_patch = ps.gs_mol_patch, gs_mol_sun_patch = ps.gs_mol_sun_patch,
    gs_mol_sha_patch = ps.gs_mol_sha_patch, lmrsun_z_patch = ps.lmrsun_z_patch,
    lmrsha_z_patch = ps.lmrsha_z_patch, psnsun_z_patch = ps.psnsun_z_patch,
    psnsha_z_patch = ps.psnsha_z_patch, rssun_z_patch = ps.rssun_z_patch,
    rssha_z_patch = ps.rssha_z_patch, cisun_z_patch = ps.cisun_z_patch,
    cisha_z_patch = ps.cisha_z_patch,
    ac_phs_patch = ps.ac_phs_patch, aj_phs_patch = ps.aj_phs_patch,
    ap_phs_patch = ps.ap_phs_patch, ag_phs_patch = ps.ag_phs_patch,
    vcmax_z_phs_patch = ps.vcmax_z_phs_patch, tpu_z_phs_patch = ps.tpu_z_phs_patch,
    kp_z_phs_patch = ps.kp_z_phs_patch,
    vcmx25_z_patch = ps.vcmx25_z_patch, jmx25_z_patch = ps.jmx25_z_patch,
    an_sun_patch = ps.an_sun_patch,
    an_sha_patch = ps.an_sha_patch, gs_mol_sun_ln_patch = ps.gs_mol_sun_ln_patch,
    gs_mol_sha_ln_patch = ps.gs_mol_sha_ln_patch,
    luvcmax25top_patch = ps.luvcmax25top_patch, lujmax25top_patch = ps.lujmax25top_patch,
    lutpu25top_patch = ps.lutpu25top_patch,
    alphapsnsun_patch = ps.alphapsnsun_patch, alphapsnsha_patch = ps.alphapsnsha_patch,
    stomatalcond_mtd = ps.stomatalcond_mtd)

# Device-view bundle of the PHS vulnerability-curve / kinetics params. params_inst
# is a MUTABLE global PhotoParamsData (non-bitstype) → can't pass by value to a
# Metal kernel; this immutable @adapt_structure bundle of the needed param arrays
# can. Built device-resident at the working precision by _psn_phs_params.
Base.@kwdef struct PsnPhsParams{V,M,S}
    krmax::V; theta_cj::V
    kmax::M; psi50::M; ck::M
    theta_ip::S
end
Adapt.@adapt_structure PsnPhsParams

# Build a PsnPhsParams on the backend of `tmpl` (any state array) at its eltype,
# copying the needed params_inst fields. `T.(a)` first so copyto! into a Float32
# device array converts from the Float64 host params (copyto! across eltypes on a
# GPU array otherwise errors). CPU: T=Float64, T.(a) is a plain copy.
function _psn_phs_params(tmpl, params::PhotoParamsData=params_inst)
    T = eltype(tmpl)
    cpv(a) = copyto!(similar(tmpl, T, size(a)), T.(a))
    return PsnPhsParams(; krmax = cpv(params.krmax), theta_cj = cpv(params.theta_cj),
        kmax = cpv(params.kmax), psi50 = cpv(params.psi50), ck = cpv(params.ck),
        theta_ip = T(params.theta_ip))
end

"""
    ci_func!(ci, p, iv, forc_pbot_c, gb_mol, je, cair, oair, lmr_z, par_z,
             rh_can, photosyns)

Evaluate f(ci) = ci - (ca - (1.37rb+1.65rs))*patm*an.
Returns `(fval, gs_mol)`.
"""
function _ci_func_core!(ci::Real, p::Int, iv::Int, forc_pbot_c::Real,
                  gb_mol::Real, je::Real, cair::Real, oair::Real,
                  lmr_z::Real, par_z::Real, rh_can::Real,
                  ps,
                  c3psn_val::Real, medlynslope_val::Real,
                  medlynintercept_val::Real, mbbopt_val::Real,
                  theta_cj_val::Real, theta_ip_val::Real)
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
    # theta_cj_val / theta_ip_val are passed in (were params_inst.theta_cj[1] /
    # params_inst.theta_ip); threaded as args so this core is GPU-kernel callable.
    # eltype-generic (T from ci); byte-identical to Float64 literals on CPU.
    T = typeof(ci)

    gs_mol = zero(T)

    if c3flag
        # Guard against ko_p=0 at very cold temperatures (Arrhenius → 0)
        if ko_p > zero(T)
            ac[p, iv] = vcmax_z[p, iv] * smooth_max(ci - cp_p, zero(ci)) / (ci + kc_p * (one(T) + oair / ko_p))
        else
            ac[p, iv] = zero(T)
        end
        denom_j = T(4.0) * ci + T(8.0) * cp_p
        aj[p, iv] = denom_j > zero(T) ? je * smooth_max(ci - cp_p, zero(ci)) / denom_j : zero(T)
        ap[p, iv] = T(3.0) * tpu_z[p, iv]
    else
        ac[p, iv] = vcmax_z[p, iv]
        aj[p, iv] = qe_p * par_z * T(4.6)
        ap[p, iv] = kp_z[p, iv] * smooth_max(ci, zero(ci)) / forc_pbot_c
    end

    # Gross photosynthesis: co-limit ac and aj, then co-limit ap
    aquad = theta_cj_val
    bquad = -(ac[p, iv] + aj[p, iv])
    cquad = ac[p, iv] * aj[p, iv]
    r1, r2 = quadratic_solve(aquad, bquad, cquad)
    ai = r2

    aquad = theta_ip_val
    bquad = -(ai + ap[p, iv])
    cquad = ai * ap[p, iv]
    r1, r2 = quadratic_solve(aquad, bquad, cquad)
    ag[p, iv] = max(zero(r1), r2)

    an[p, iv] = ag[p, iv] - lmr_z
    if an[p, iv] < zero(T)
        return (zero(T), gs_mol)
    end

    cs = cair - T(1.4) / gb_mol * an[p, iv] * forc_pbot_c
    cs = smooth_max(cs, T(MAX_CS))

    if stomatalcond_mtd == STOMATALCOND_MTD_MEDLYN2011
        term = T(1.6) * an[p, iv] / (cs / forc_pbot_c * T(1.0e06))
        aquad = one(T)
        bquad = -(T(2.0) * (medlynintercept_val * T(1.0e-06) + term) +
                  (medlynslope_val * term)^2 / (gb_mol * T(1.0e-06) * rh_can))
        cquad = medlynintercept_val^2 * T(1.0e-12) +
                (T(2.0) * medlynintercept_val * T(1.0e-06) + term *
                 (one(T) - medlynslope_val^2 / rh_can)) * term
        r1, r2 = quadratic_solve(aquad, bquad, cquad)
        gs_mol = max(r1 * T(1.0e06), one(T))
    elseif stomatalcond_mtd == STOMATALCOND_MTD_BB1987
        aquad = cs
        bquad = cs * (gb_mol - bbb_p) - mbb_p * an[p, iv] * forc_pbot_c
        cquad = -gb_mol * (cs * bbb_p + mbb_p * an[p, iv] * forc_pbot_c * rh_can)
        r1, r2 = quadratic_solve(aquad, bquad, cquad)
        gs_mol = max(r1, bbb_p)
    end

    fval = ci - cair + an[p, iv] * forc_pbot_c * (T(1.4) * gs_mol + T(1.6) * gb_mol) / (gb_mol * gs_mol)

    return (fval, gs_mol)
end

# Back-compat / host wrapper: reads the module-global photosynthesis params
# (theta_cj[1], theta_ip) and forwards to the positional core. Callers on the host
# (tests, the non-kernel path) use this; GPU kernels call _ci_func_core! directly.
function ci_func!(ci::Real, p::Int, iv::Int, forc_pbot_c::Real,
                  gb_mol::Real, je::Real, cair::Real, oair::Real,
                  lmr_z::Real, par_z::Real, rh_can::Real,
                  ps::PhotosynthesisData;
                  c3psn_val::Real=1.0,
                  medlynslope_val::Real=6.0,
                  medlynintercept_val::Real=100.0,
                  mbbopt_val::Real=9.0)
    return _ci_func_core!(ci, p, iv, forc_pbot_c, gb_mol, je, cair, oair,
                          lmr_z, par_z, rh_can, ps,
                          c3psn_val, medlynslope_val, medlynintercept_val, mbbopt_val,
                          params_inst.theta_cj[1], params_inst.theta_ip)
end

# =====================================================================
# hybrid_solver! — hybrid Newton-secant / Brent solver for ci
# =====================================================================

"""
    hybrid_solver!(x0, p, iv, forc_pbot_c, gb_mol, je, cair, oair,
                   lmr_z, par_z, rh_can, ps; kwargs...)

Hybrid solver for ci. Returns `(ci_solution, gs_mol, niter)`.
"""
function _hybrid_solver_core!(x0::Real, p::Int, iv::Int, forc_pbot_c::Real,
                        gb_mol::Real, je::Real, cair::Real, oair::Real,
                        lmr_z::Real, par_z::Real, rh_can::Real,
                        ps,
                        c3psn_val::Real, medlynslope_val::Real,
                        medlynintercept_val::Real, mbbopt_val::Real,
                        theta_cj_val::Real, theta_ip_val::Real)
    T = typeof(x0)
    eps_val = T(1.0e-2)
    eps1 = T(1.0e-4)
    itmax = 40

    f0, gs_mol = _ci_func_core!(x0, p, iv, forc_pbot_c, gb_mol, je, cair, oair,
                           lmr_z, par_z, rh_can, ps,
                           c3psn_val, medlynslope_val, medlynintercept_val, mbbopt_val,
                           theta_cj_val, theta_ip_val)

    if f0 == zero(T)
        return (x0, gs_mol, 0)
    end

    minx = x0
    minf = abs(f0)
    x1 = x0 * T(0.99)

    f1, gs_mol = _ci_func_core!(x1, p, iv, forc_pbot_c, gb_mol, je, cair, oair,
                           lmr_z, par_z, rh_can, ps,
                           c3psn_val, medlynslope_val, medlynintercept_val, mbbopt_val,
                           theta_cj_val, theta_ip_val)

    if f1 == zero(T)
        return (x1, gs_mol, 0)
    end
    if abs(f1) < minf
        minx = x1
        minf = abs(f1)
    end

    iter = 0
    while true
        iter += 1
        if (f1 - f0) == zero(T)
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

        f1, gs_mol = _ci_func_core!(x1, p, iv, forc_pbot_c, gb_mol, je, cair, oair,
                               lmr_z, par_z, rh_can, ps,
                               c3psn_val, medlynslope_val, medlynintercept_val, mbbopt_val,
                               theta_cj_val, theta_ip_val)

        if abs(f1) < minf
            minx = x1
            minf = abs(f1)
        end
        if abs(f1) <= eps1
            x0 = x1
            break
        end

        if f1 * f0 < zero(T)
            x_brent = _brent_solver_core!(x0, x1, f0, f1, tol, p, iv, forc_pbot_c,
                                     gb_mol, je, cair, oair, lmr_z, par_z,
                                     rh_can, ps,
                                     c3psn_val, medlynslope_val, medlynintercept_val, mbbopt_val,
                                     theta_cj_val, theta_ip_val)
            x0 = x_brent
            # get final gs_mol
            _, gs_mol = _ci_func_core!(x0, p, iv, forc_pbot_c, gb_mol, je, cair, oair,
                                  lmr_z, par_z, rh_can, ps,
                                  c3psn_val, medlynslope_val, medlynintercept_val, mbbopt_val,
                                  theta_cj_val, theta_ip_val)
            break
        end

        if iter > itmax
            _, gs_mol = _ci_func_core!(minx, p, iv, forc_pbot_c, gb_mol, je, cair, oair,
                                  lmr_z, par_z, rh_can, ps,
                                  c3psn_val, medlynslope_val, medlynintercept_val, mbbopt_val,
                                  theta_cj_val, theta_ip_val)
            x0 = minx
            break
        end
    end

    return (x0, gs_mol, iter)
end

# Back-compat / host wrapper (reads module-global theta_cj[1]/theta_ip).
function hybrid_solver!(x0::Real, p::Int, iv::Int, forc_pbot_c::Real,
                        gb_mol::Real, je::Real, cair::Real, oair::Real,
                        lmr_z::Real, par_z::Real, rh_can::Real,
                        ps::PhotosynthesisData;
                        c3psn_val::Real=1.0, medlynslope_val::Real=6.0,
                        medlynintercept_val::Real=100.0, mbbopt_val::Real=9.0)
    return _hybrid_solver_core!(x0, p, iv, forc_pbot_c, gb_mol, je, cair, oair,
                          lmr_z, par_z, rh_can, ps,
                          c3psn_val, medlynslope_val, medlynintercept_val, mbbopt_val,
                          params_inst.theta_cj[1], params_inst.theta_ip)
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
function _brent_solver_core!(x1::Real, x2::Real, f1::Real, f2::Real,
                       tol::Real, p::Int, iv::Int, forc_pbot_c::Real,
                       gb_mol::Real, je::Real, cair::Real, oair::Real,
                       lmr_z::Real, par_z::Real, rh_can::Real,
                       ps,
                       c3psn_val::Real, medlynslope_val::Real,
                       medlynintercept_val::Real, mbbopt_val::Real,
                       theta_cj_val::Real, theta_ip_val::Real)
    T = typeof(x1)
    itmax = 20
    eps_val = T(1.0e-2)

    a = x1
    b = x2
    fa = f1
    fb = f2

    c = b
    fc = fb
    d = b - a
    e = d

    for iter in 1:itmax
        if (fb > zero(T) && fc > zero(T)) || (fb < zero(T) && fc < zero(T))
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
        tol1 = T(2.0) * eps_val * abs(b) + T(0.5) * tol
        xm = T(0.5) * (c - b)

        if abs(xm) <= tol1 || fb == zero(T)
            return b
        end

        if abs(e) >= tol1 && abs(fa) > abs(fb)
            s = fb / fa
            if a == c
                p_val = T(2.0) * xm * s
                q_val = one(T) - s
            else
                q_val = fa / fc
                r_val = fb / fc
                p_val = s * (T(2.0) * xm * q_val * (q_val - r_val) - (b - a) * (r_val - one(T)))
                q_val = (q_val - one(T)) * (r_val - one(T)) * (s - one(T))
            end
            if p_val > zero(T)
                q_val = -q_val
            end
            p_val = abs(p_val)
            if T(2.0) * p_val < min(T(3.0) * xm * q_val - abs(tol1 * q_val), abs(e * q_val))
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

        fb, _ = _ci_func_core!(b, p, iv, forc_pbot_c, gb_mol, je, cair, oair,
                          lmr_z, par_z, rh_can, ps,
                          c3psn_val, medlynslope_val, medlynintercept_val, mbbopt_val,
                          theta_cj_val, theta_ip_val)

        if fb == zero(T)
            return b
        end
    end

    return b
end

# Back-compat / host wrapper (reads module-global theta_cj[1]/theta_ip).
function brent_solver!(x1::Real, x2::Real, f1::Real, f2::Real,
                       tol::Real, p::Int, iv::Int, forc_pbot_c::Real,
                       gb_mol::Real, je::Real, cair::Real, oair::Real,
                       lmr_z::Real, par_z::Real, rh_can::Real,
                       ps::PhotosynthesisData;
                       c3psn_val::Real=1.0, medlynslope_val::Real=6.0,
                       medlynintercept_val::Real=100.0, mbbopt_val::Real=9.0)
    return _brent_solver_core!(x1, x2, f1, f2, tol, p, iv, forc_pbot_c,
                          gb_mol, je, cair, oair, lmr_z, par_z, rh_can, ps,
                          c3psn_val, medlynslope_val, medlynintercept_val, mbbopt_val,
                          params_inst.theta_cj[1], params_inst.theta_ip)
end

# =====================================================================
# Pass 3 kernel — leaf-level photosynthesis & stomatal conductance
#
# ONE THREAD PER PATCH: p = @index(Global) over ndrange = length(bounds_patch)
# (bounds_patch is 1:np; the per-(p,iv) state matrices are allocated np×nlevcan
# and indexed [p,iv]). The inner `for iv in 1:nrad[p]` loop stays INSIDE the
# kernel so each thread walks its own iv's serially — the per-patch outputs
# gb_mol[p], vpd_can[p], rh_leaf[p] are written within the iv loop, so flattening
# to (p,iv) would race; per-patch is race-free and matches serial semantics.
#
# This kernel calls the GPU-callable Ci-solver chain (_hybrid_solver_core! →
# _ci_func_core! → _brent_solver_core!) which read NO globals.
#
# ps (PhotosynthesisData) is passed so the solver cores can index its fields. The
# bare array args (ac, aj, ... ci_z, gs_mol, ...) are LOCAL aliases of ps fields /
# scratch matrices established in photosynthesis!; on the CPU backend they are the
# SAME memory as the corresponding ps fields the solver writes, so writes stay
# consistent. Written arrays are non-@Const; read-only are @Const.
# =====================================================================
# Scalar constants bundled into one isbits arg (Metal ~31 kernel-arg limit; the
# flat form passed ~50 args). The "written" arrays are ps fields, accessed via ps
# inside (phase ones selected by is_sun) rather than passed individually; only the
# local scratch psn_w*_z stay as args.
struct PsnCiScalars{S}
    fnps::S; theta_psii::S; theta_cj_1::S; theta_ip_val::S; medlyn_slope_override::S
    RGAS_::S; MAX_CS_::S; MEDLYN_RH_CAN_MAX_::S; MEDLYN_RH_CAN_FACT_::S; rsmax0::S
end

# AD-mode switch for the Pass-3 ci-solve. When true, `photosynth_ci_solve!` runs a
# plain host loop (over `_photosynth_ci_body!`) instead of launching the KA kernel —
# the KA kernel-launch path makes Enzyme reverse segfault on Julia 1.12. The reverse
# engine (compositional_reverse!) sets this true for the duration of the autodiff and
# restores it; normal forward/GPU runs leave it false (KA kernel). Byte-identical primal
# is guaranteed because both paths execute the same `_photosynth_ci_body!`.
const _PSN_CI_AD_HOSTLOOP = Ref(false)

# PARITY localization (Phase 3c): opt-in debug dump of converged PHS sunlit-leaf
# photosynthesis internals (vcmax_z → an → gs chain) for patches in PHS_PHOTO_DEBUG_PATCHES.
# Default OFF → forward/GPU/AD runs are byte-identical (the dump is a read-only loop
# over already-converged state in the host launcher, gated by this flag). Mirrors the
# Fortran FPHOTO/FVCM25 SourceMods dump in PhotosynthesisMod.F90.
const PHS_PHOTO_DEBUG          = Ref(false)
const PHS_PHOTO_DEBUG_PATH     = Ref("photo_internals_julia.txt")
const PHS_PHOTO_DEBUG_PATCHES  = Ref{Vector{Int}}([2, 3])

@kernel function _photosynth_ci_kernel!(ps, psn_wc_z, psn_wj_z, psn_wp_z,
        @Const(mask_patch), @Const(ivt),
        @Const(medlynslope_pft), @Const(medlynintercept_pft), @Const(mbbopt_pft),
        @Const(forc_pbot), @Const(tgcm), @Const(rb), @Const(par_z_in),
        @Const(jmax_z_local), @Const(eair), @Const(esat_tv), @Const(cair),
        @Const(oair), @Const(o3coefg), @Const(o3coefv), @Const(nrad),
        sc, is_sun::Bool, stomatalcond_mtd::Int,
        STOMATALCOND_MTD_BB1987_::Int, STOMATALCOND_MTD_MEDLYN2011_::Int)

    p = @index(Global)
    _photosynth_ci_body!(p, ps, psn_wc_z, psn_wj_z, psn_wp_z,
        mask_patch, ivt, medlynslope_pft, medlynintercept_pft, mbbopt_pft,
        forc_pbot, tgcm, rb, par_z_in, jmax_z_local, eair, esat_tv, cair,
        oair, o3coefg, o3coefv, nrad, sc, is_sun, stomatalcond_mtd,
        STOMATALCOND_MTD_BB1987_, STOMATALCOND_MTD_MEDLYN2011_)
end

# Per-patch Pass-3 body, factored out of `_photosynth_ci_kernel!` so the SAME code
# runs in two ways: (a) the KA kernel above calls it with `p = @index(Global)` (KA
# inlines it on the GPU device path → the Metal/CUDA kernel is functionally unchanged);
# (b) the host loop in `photosynth_ci_solve!` (taken only under `_PSN_CI_AD_HOSTLOOP[]`,
# i.e. Enzyme reverse on Julia 1.12) calls it directly, avoiding the KA-kernel-launch
# codegen that Enzyme reverse segfaults on under 1.12. Identical body ⇒ byte-identical
# primal between the two paths.
@inline function _photosynth_ci_body!(p, ps, psn_wc_z, psn_wj_z, psn_wp_z,
        mask_patch, ivt, medlynslope_pft, medlynintercept_pft, mbbopt_pft,
        forc_pbot, tgcm, rb, par_z_in, jmax_z_local, eair, esat_tv, cair,
        oair, o3coefg, o3coefv, nrad, sc, is_sun::Bool, stomatalcond_mtd::Int,
        STOMATALCOND_MTD_BB1987_::Int, STOMATALCOND_MTD_MEDLYN2011_::Int)
    @inbounds if mask_patch[p]
        ivt_p = ivt[p]
        T = eltype(forc_pbot)
        # Unpack the scalar bundle + alias ps fields (phase-selected via is_sun)
        # to locals so the body below is unchanged.
        fnps = sc.fnps; theta_psii = sc.theta_psii; theta_cj_1 = sc.theta_cj_1
        theta_ip_val = sc.theta_ip_val; medlyn_slope_override = sc.medlyn_slope_override
        RGAS_ = sc.RGAS_; MAX_CS_ = sc.MAX_CS_; MEDLYN_RH_CAN_MAX_ = sc.MEDLYN_RH_CAN_MAX_
        MEDLYN_RH_CAN_FACT_ = sc.MEDLYN_RH_CAN_FACT_; rsmax0 = sc.rsmax0
        ac = ps.ac_patch; aj = ps.aj_patch; ap = ps.ap_patch; ag = ps.ag_patch; an = ps.an_patch
        gs_mol = ps.gs_mol_patch; rh_leaf = ps.rh_leaf_patch; gb_mol = ps.gb_mol_patch
        vpd_can = ps.vpd_can_patch; bbb = ps.bbb_patch; c3flag = ps.c3flag_patch
        lmr_z = is_sun ? ps.lmrsun_z_patch : ps.lmrsha_z_patch
        psn_z = is_sun ? ps.psnsun_z_patch : ps.psnsha_z_patch
        rs_z  = is_sun ? ps.rssun_z_patch  : ps.rssha_z_patch
        ci_z  = is_sun ? ps.cisun_z_patch  : ps.cisha_z_patch

        # Medlyn slope: use override if set, else PFT default
        medlynslope_p = isnan(medlyn_slope_override) ? medlynslope_pft[ivt_p] : medlyn_slope_override

        cf = forc_pbot[p] / (RGAS_ * tgcm[p]) * T(1.0e06)
        gb = one(T) / rb[p]
        gb_mol[p] = gb * cf

        for iv in 1:nrad[p]
            if par_z_in[p, iv] <= zero(T)  # night time
                ac[p, iv] = zero(T)
                aj[p, iv] = zero(T)
                ap[p, iv] = zero(T)
                ag[p, iv] = zero(T)
                an[p, iv] = ag[p, iv] - lmr_z[p, iv]
                psn_z[p, iv] = zero(T)
                psn_wc_z[p, iv] = zero(T)
                psn_wj_z[p, iv] = zero(T)
                psn_wp_z[p, iv] = zero(T)
                if stomatalcond_mtd == STOMATALCOND_MTD_BB1987_
                    rs_z[p, iv] = smooth_min(rsmax0, one(T) / bbb[p] * cf)
                elseif stomatalcond_mtd == STOMATALCOND_MTD_MEDLYN2011_
                    rs_z[p, iv] = smooth_min(rsmax0, one(T) / medlynintercept_pft[ivt_p] * cf)
                end
                ci_z[p, iv] = zero(T)
                rh_leaf[p] = zero(T)
                gs_mol[p, iv] = stomatalcond_mtd == STOMATALCOND_MTD_BB1987_ ? bbb[p] : medlynintercept_pft[ivt_p]
            else  # day time
                ceair = smooth_min(eair[p], esat_tv[p])
                if stomatalcond_mtd == STOMATALCOND_MTD_BB1987_
                    rh_can = ceair / esat_tv[p]
                else
                    rh_can = smooth_max((esat_tv[p] - ceair), MEDLYN_RH_CAN_MAX_) * MEDLYN_RH_CAN_FACT_
                    vpd_can[p] = rh_can
                end

                # Electron transport rate
                qabs = T(0.5) * (one(T) - fnps) * par_z_in[p, iv] * T(4.6)
                aquad = theta_psii
                bquad = -(qabs + jmax_z_local[p, iv])
                cquad = qabs * jmax_z_local[p, iv]
                r1, r2 = quadratic_solve(aquad, bquad, cquad)
                je = r2

                # Initial guess for ci
                if c3flag[p]
                    ci_z[p, iv] = T(0.7) * cair[p]
                else
                    ci_z[p, iv] = T(0.4) * cair[p]
                end

                # Solve for ci and gs (positional GPU-callable core).
                # c3psn_val = 1.0 (wrapper default); theta_cj_1 = params_inst.theta_cj[1]
                # (NOT theta_cj[ivt_p] — that was dead code never threaded to the solver).
                ci_sol, gs_mol_val, _ = _hybrid_solver_core!(ci_z[p, iv], p, iv,
                    forc_pbot[p], gb_mol[p], je, cair[p], oair[p],
                    lmr_z[p, iv], par_z_in[p, iv], rh_can, ps,
                    one(T), medlynslope_p, medlynintercept_pft[ivt_p], mbbopt_pft[ivt_p],
                    theta_cj_1, theta_ip_val)

                gs_mol[p, iv] = gs_mol_val

                # Check for an < 0
                if an[p, iv] < zero(T)
                    if stomatalcond_mtd == STOMATALCOND_MTD_BB1987_
                        gs_mol[p, iv] = bbb[p]
                    else
                        gs_mol[p, iv] = medlynintercept_pft[ivt_p]
                    end
                end

                # Store sunlit/shaded gs
                if is_sun
                    ps.gs_mol_sun_patch[p, iv] = gs_mol[p, iv]
                else
                    ps.gs_mol_sha_patch[p, iv] = gs_mol[p, iv]
                end

                # Final ci
                cs = cair[p] - T(1.4) / gb_mol[p] * an[p, iv] * forc_pbot[p]
                cs = smooth_max(cs, MAX_CS_)
                ci_z[p, iv] = cair[p] - an[p, iv] * forc_pbot[p] *
                    (T(1.4) * gs_mol[p, iv] + T(1.6) * gb_mol[p]) / (gb_mol[p] * gs_mol[p, iv])
                ci_z[p, iv] = smooth_max(ci_z[p, iv], T(1.0e-06))

                # Convert to resistance (gs must be positive)
                gs_val = max(gs_mol[p, iv], one(T)) / cf
                rs_z[p, iv] = smooth_min(one(T) / gs_val, rsmax0)
                rs_z[p, iv] = rs_z[p, iv] / o3coefg[p]

                # Photosynthesis
                psn_z[p, iv] = ag[p, iv]
                psn_z[p, iv] = psn_z[p, iv] * o3coefv[p]

                psn_wc_z[p, iv] = zero(T)
                psn_wj_z[p, iv] = zero(T)
                psn_wp_z[p, iv] = zero(T)

                if ac[p, iv] <= aj[p, iv] && ac[p, iv] <= ap[p, iv]
                    psn_wc_z[p, iv] = psn_z[p, iv]
                elseif aj[p, iv] < ac[p, iv] && aj[p, iv] <= ap[p, iv]
                    psn_wj_z[p, iv] = psn_z[p, iv]
                elseif ap[p, iv] < ac[p, iv] && ap[p, iv] < aj[p, iv]
                    psn_wp_z[p, iv] = psn_z[p, iv]
                end

                # Ball-Berry error check
                if stomatalcond_mtd == STOMATALCOND_MTD_BB1987_
                    hs = (gb_mol[p] * ceair + gs_mol[p, iv] * esat_tv[p]) /
                         ((gb_mol[p] + gs_mol[p, iv]) * esat_tv[p])
                    rh_leaf[p] = hs
                end
            end
        end
    end
end

"""
    photosynth_ci_solve!(ps, ac, aj, ..., bounds_patch, ...)

Host launcher for the Pass-3 leaf photosynthesis / stomatal-conductance kernel
(`_photosynth_ci_kernel!`). GPU-hostile values (params_inst fields, the Medlyn
override, the `phase == "sun"` string compare, module-global constants) are
resolved to host scalars before the launch; the kernel reads no globals/Strings.
One thread per patch; the per-canopy-layer `iv` loop runs serially inside.
"""
function photosynth_ci_solve!(ps,
        ac, aj, ap, ag, an, psn_z, psn_wc_z, psn_wj_z, psn_wp_z,
        rs_z, ci_z, gs_mol, rh_leaf, gb_mol, vpd_can,
        mask_patch, ivt,
        medlynslope_pft, medlynintercept_pft, mbbopt_pft,
        forc_pbot, tgcm, rb, par_z_in, lmr_z, jmax_z_local, c3flag,
        eair, esat_tv, cair, oair, bbb, o3coefg, o3coefv, nrad,
        bounds_patch::UnitRange{Int}, phase::String,
        stomatalcond_mtd::Int, rsmax0,
        overrides::CalibrationOverrides)
    # Resolve GPU-hostile params to host scalars, eltype-converted + bundled into
    # one isbits arg (Metal ~31-arg limit). The ps-field arrays (ac/aj/.../ci_z/
    # rs_z/psn_z/gs_mol/…) are accessed via ps inside the kernel, not passed.
    # Struct-first kernel → manual backend + KA.synchronize.
    T = eltype(forc_pbot)
    sc = PsnCiScalars{T}(T(params_inst.fnps), T(params_inst.theta_psii),
        T(params_inst.theta_cj[1]), T(params_inst.theta_ip), T(overrides.medlyn_slope),
        T(RGAS_PSN), T(MAX_CS), T(MEDLYN_RH_CAN_MAX), T(MEDLYN_RH_CAN_FACT), T(rsmax0))
    is_sun = (phase == "sun")
    dv = _psn_dv(ps)
    if _PSN_CI_AD_HOSTLOOP[]
        # AD-mode path (set only by the reverse engine; see compositional_reverse!):
        # plain host loop over the SAME per-patch body the KA kernel runs. This avoids
        # the KernelAbstractions kernel-launch codegen that Enzyme reverse segfaults on
        # under Julia 1.12, while staying byte-identical to the kernel primal (shared
        # `_photosynth_ci_body!`). Normal/GPU runs never set the flag → KA kernel below.
        @inbounds for p in 1:length(bounds_patch)
            _photosynth_ci_body!(p, dv, psn_wc_z, psn_wj_z, psn_wp_z,
                mask_patch, ivt, medlynslope_pft, medlynintercept_pft, mbbopt_pft,
                forc_pbot, tgcm, rb, par_z_in, jmax_z_local, eair, esat_tv, cair, oair,
                o3coefg, o3coefv, nrad, sc, is_sun, stomatalcond_mtd,
                STOMATALCOND_MTD_BB1987, STOMATALCOND_MTD_MEDLYN2011)
        end
    else
        be = _kernel_backend(forc_pbot)
        _photosynth_ci_kernel!(be)(dv, psn_wc_z, psn_wj_z, psn_wp_z,
            mask_patch, ivt, medlynslope_pft, medlynintercept_pft, mbbopt_pft,
            forc_pbot, tgcm, rb, par_z_in, jmax_z_local, eair, esat_tv, cair, oair,
            o3coefg, o3coefv, nrad, sc, is_sun, stomatalcond_mtd,
            STOMATALCOND_MTD_BB1987, STOMATALCOND_MTD_MEDLYN2011;
            ndrange = length(bounds_patch))
        KA.synchronize(be)
    end
    return ps
end

# =====================================================================
# photosynthesis! per-patch passes (kernelized). `ps` (PhotosynthesisData,
# @adapt_structure'd) is passed whole as ONE arg and its fields are written
# directly inside the kernels; params_inst scalar fields are bundled into an
# isbits NamedTuple built on the host (module globals are not GPU-safe). All
# Float64 literals are eltype-converted (T = eltype of an input). Launched via
# manual backend + KA.synchronize (struct-first kernel). One thread per patch.
# =====================================================================

# Pass 1: C3/C4 flag, quantum efficiency, Ball-Berry intercept, Michaelis-Menten.
@kernel function _psn_pass1_kernel!(ps, @Const(mask_patch), @Const(ivt),
        @Const(c3psn_pft), @Const(mbbopt_pft), @Const(forc_pbot), @Const(oair),
        @Const(t_veg), @Const(btran), prm, stomatalcond_mtd::Int,
        bbbopt_c3, bbbopt_c4, stomatal_bb::Int)
    p = @index(Global)
    _psn_pass1_body!(p, ps, mask_patch, ivt, c3psn_pft, mbbopt_pft, forc_pbot, oair,
        t_veg, btran, prm, stomatalcond_mtd, bbbopt_c3, bbbopt_c4, stomatal_bb)
end

# Shared per-patch body for Pass 1 (KA kernel + AD host loop) — see _photosynth_ci_body!.
@inline function _psn_pass1_body!(p, ps, mask_patch, ivt, c3psn_pft, mbbopt_pft,
        forc_pbot, oair, t_veg, btran, prm, stomatalcond_mtd::Int,
        bbbopt_c3, bbbopt_c4, stomatal_bb::Int)
    @inbounds if mask_patch[p]
        T = eltype(forc_pbot)
        ivt_p = ivt[p]
        # c3psn is 0/1; round(Int,·)==1 ⟺ >0.5. round(Int,Float32) heap-allocates
        # on Metal (gpu_malloc), so use the comparison (byte-identical for 0/1).
        if c3psn_pft[ivt_p] > T(0.5)
            ps.c3flag_patch[p] = true
            ps.qe_patch[p] = zero(T)
            bbbopt_p = bbbopt_c3
        else
            ps.c3flag_patch[p] = false
            ps.qe_patch[p] = T(0.05)
            bbbopt_p = bbbopt_c4
        end
        if stomatalcond_mtd == stomatal_bb
            ps.bbb_patch[p] = smooth_max(bbbopt_p * btran[p], one(T))
            ps.mbb_patch[p] = mbbopt_pft[ivt_p]
        end
        kc25 = prm.kc25_coef * forc_pbot[p]
        ko25 = prm.ko25_coef * forc_pbot[p]
        sco  = T(0.5) * T(0.209) / prm.cp25_yr2000
        cp25 = T(0.5) * oair[p] / sco
        ps.kc_patch[p] = kc25 * ft_photo(t_veg[p], prm.kcha)
        ps.ko_patch[p] = ko25 * ft_photo(t_veg[p], prm.koha)
        ps.cp_patch[p] = cp25 * ft_photo(t_veg[p], prm.cpha)
    end
end

function psn_pass1_update!(ps, mask_patch, ivt, c3psn_pft, mbbopt_pft,
        forc_pbot, oair, t_veg, btran, stomatalcond_mtd::Int, bounds_patch)
    T = eltype(forc_pbot)
    prm = (kc25_coef = T(params_inst.kc25_coef), ko25_coef = T(params_inst.ko25_coef),
           cp25_yr2000 = T(params_inst.cp25_yr2000), kcha = T(params_inst.kcha),
           koha = T(params_inst.koha), cpha = T(params_inst.cpha))
    dv = _psn_dv(ps)
    if _PSN_CI_AD_HOSTLOOP[]
        @inbounds for p in 1:length(bounds_patch)
            _psn_pass1_body!(p, dv, mask_patch, ivt, c3psn_pft, mbbopt_pft,
                forc_pbot, oair, t_veg, btran, prm, stomatalcond_mtd,
                T(BBBOPT_C3), T(BBBOPT_C4), STOMATALCOND_MTD_BB1987)
        end
    else
        be = _kernel_backend(forc_pbot)
        _psn_pass1_kernel!(be)(dv, mask_patch, ivt, c3psn_pft, mbbopt_pft,
            forc_pbot, oair, t_veg, btran, prm, stomatalcond_mtd,
            T(BBBOPT_C3), T(BBBOPT_C4), STOMATALCOND_MTD_BB1987;
            ndrange = length(bounds_patch))
        KA.synchronize(be)
    end
    return nothing
end

# Pass 4: canopy integration over layers → per-patch psn/rs/lmr. Phase-dependent
# ps fields selected via is_sun (no String compare in-kernel). Internal iv loop
# runs serially per thread (per-patch reduction, no race). psn_wc_z/wj_z/wp_z are
# the local scratch matrices written by the Ci solve (Pass 3).
@kernel function _psn_pass4_kernel!(ps, @Const(mask_patch), @Const(nrad),
        @Const(lai_z_in), @Const(rb), @Const(psn_wc_z), @Const(psn_wj_z),
        @Const(psn_wp_z), is_sun::Bool)
    p = @index(Global)
    _psn_pass4_body!(p, ps, mask_patch, nrad, lai_z_in, rb,
        psn_wc_z, psn_wj_z, psn_wp_z, is_sun)
end

# Shared per-patch body for Pass 4 (KA kernel + AD host loop) — see _photosynth_ci_body!.
@inline function _psn_pass4_body!(p, ps, mask_patch, nrad, lai_z_in, rb,
        psn_wc_z, psn_wj_z, psn_wp_z, is_sun::Bool)
    @inbounds if mask_patch[p]
        T = eltype(rb)
        psn_z   = is_sun ? ps.psnsun_z_patch  : ps.psnsha_z_patch
        lmr_z   = is_sun ? ps.lmrsun_z_patch  : ps.lmrsha_z_patch
        rs_z    = is_sun ? ps.rssun_z_patch   : ps.rssha_z_patch
        psn_o   = is_sun ? ps.psnsun_patch    : ps.psnsha_patch
        psnwc_o = is_sun ? ps.psnsun_wc_patch : ps.psnsha_wc_patch
        psnwj_o = is_sun ? ps.psnsun_wj_patch : ps.psnsha_wj_patch
        psnwp_o = is_sun ? ps.psnsun_wp_patch : ps.psnsha_wp_patch
        lmr_o   = is_sun ? ps.lmrsun_patch    : ps.lmrsha_patch
        rs_o    = is_sun ? ps.rssun_patch     : ps.rssha_patch

        psncan = zero(T); psncan_wc = zero(T); psncan_wj = zero(T); psncan_wp = zero(T)
        lmrcan = zero(T); gscan = zero(T); laican = zero(T)
        for iv in 1:nrad[p]
            psncan    += psn_z[p, iv]    * lai_z_in[p, iv]
            psncan_wc += psn_wc_z[p, iv] * lai_z_in[p, iv]
            psncan_wj += psn_wj_z[p, iv] * lai_z_in[p, iv]
            psncan_wp += psn_wp_z[p, iv] * lai_z_in[p, iv]
            lmrcan    += lmr_z[p, iv]    * lai_z_in[p, iv]
            gscan     += lai_z_in[p, iv] / max(rb[p] + rs_z[p, iv], T(1.0e-6))
            laican    += lai_z_in[p, iv]
        end
        if laican > zero(T)
            psn_o[p]   = psncan / laican
            psnwc_o[p] = psncan_wc / laican
            psnwj_o[p] = psncan_wj / laican
            psnwp_o[p] = psncan_wp / laican
            lmr_o[p]   = lmrcan / laican
            rs_o[p]    = laican / gscan - rb[p]
        else
            psn_o[p] = zero(T); psnwc_o[p] = zero(T); psnwj_o[p] = zero(T)
            psnwp_o[p] = zero(T); lmr_o[p] = zero(T); rs_o[p] = zero(T)
        end
    end
end

function psn_pass4_update!(ps, mask_patch, nrad, lai_z_in, rb,
        psn_wc_z, psn_wj_z, psn_wp_z, is_sun::Bool, bounds_patch)
    dv = _psn_dv(ps)
    if _PSN_CI_AD_HOSTLOOP[]
        @inbounds for p in 1:length(bounds_patch)
            _psn_pass4_body!(p, dv, mask_patch, nrad, lai_z_in, rb,
                psn_wc_z, psn_wj_z, psn_wp_z, is_sun)
        end
    else
        be = _kernel_backend(rb)
        _psn_pass4_kernel!(be)(dv, mask_patch, nrad, lai_z_in, rb,
            psn_wc_z, psn_wj_z, psn_wp_z, is_sun; ndrange = length(bounds_patch))
        KA.synchronize(be)
    end
    return nothing
end

# Pass 2: nitrogen profile + per-canopy-layer vcmax/jmax/tpu/kp/lmr. Per-patch
# with an internal serial iv loop. params_inst scalars bundled in `prm`; overrides
# / leaf_mr_vcm / flags / lmrc resolved on the host. kn + jmax_z_local are local
# scratch matrices; lmr_z is the phase-selected ps field.
@kernel function _psn_pass2_kernel!(ps, kn, jmax_z_local, @Const(mask_patch),
        @Const(ivt), @Const(slatop_pft), @Const(leafcn_pft), @Const(flnr_pft),
        @Const(fnitr_pft), @Const(lmr_intercept_atkin), @Const(dayl_factor), @Const(t10), @Const(t_veg),
        @Const(btran), @Const(tlai_z), @Const(par_z_in), @Const(vcmaxcint),
        @Const(nrad), prm, vcmax25_scale, jmax25top_sf_val, leaf_mr_vcm,
        use_cn::Bool, light_inhibit::Bool, leafresp_method::Int, nlevcan::Int,
        is_sun::Bool, leafresp_ryan::Int)
    p = @index(Global)
    _psn_pass2_body!(p, ps, kn, jmax_z_local, mask_patch, ivt, slatop_pft, leafcn_pft,
        flnr_pft, fnitr_pft, lmr_intercept_atkin, dayl_factor, t10, t_veg, btran, tlai_z,
        par_z_in, vcmaxcint, nrad, prm, vcmax25_scale, jmax25top_sf_val, leaf_mr_vcm,
        use_cn, light_inhibit, leafresp_method, nlevcan, is_sun, leafresp_ryan)
end

# Shared per-patch body for Pass 2 (KA kernel + AD host loop) — see _photosynth_ci_body!.
@inline function _psn_pass2_body!(p, ps, kn, jmax_z_local, mask_patch, ivt, slatop_pft,
        leafcn_pft, flnr_pft, fnitr_pft, lmr_intercept_atkin, dayl_factor, t10, t_veg,
        btran, tlai_z, par_z_in, vcmaxcint, nrad, prm, vcmax25_scale, jmax25top_sf_val,
        leaf_mr_vcm, use_cn::Bool, light_inhibit::Bool, leafresp_method::Int, nlevcan::Int,
        is_sun::Bool, leafresp_ryan::Int)
    @inbounds if mask_patch[p]
        T = eltype(t_veg)
        ivt_p = ivt[p]
        lmr_z = is_sun ? ps.lmrsun_z_patch : ps.lmrsha_z_patch
        vcmax_z = ps.vcmax_z_patch; tpu_z = ps.tpu_z_patch; kp_z = ps.kp_z_patch
        c3flag_p = ps.c3flag_patch[p]

        lnc_p = one(T) / (slatop_pft[ivt_p] * leafcn_pft[ivt_p])
        ps.lnca_patch[p] = lnc_p
        vcmax25top = lnc_p * flnr_pft[ivt_p] * prm.fnr * prm.act25 * dayl_factor[p]
        if !use_cn
            vcmax25top = vcmax25top * fnitr_pft[ivt_p]
        end
        if !isnan(vcmax25_scale)
            vcmax25top = vcmax25top * vcmax25_scale
        end
        jmax25top = ((T(2.59) - T(0.035) * smooth_clamp(t10[p] - T(TFRZ), T(11.0), T(35.0))) *
                     vcmax25top) * jmax25top_sf_val
        tpu25top = prm.tpu25ratio * vcmax25top
        kp25top  = prm.kp25ratio * vcmax25top

        if dayl_factor[p] < T(1.0e-12)
            kn[p] = zero(T)
        else
            kn[p] = exp(T(0.00963) * vcmax25top / dayl_factor[p] - T(2.43))
        end

        if use_cn
            if leafresp_method == leafresp_ryan
                lmr25top = T(2.525e-6) * (T(1.5)^((T(25.0) - T(20.0)) / T(10.0)))
                lmr25top = lmr25top * lnc_p / T(12.0e-06)
            else
                # Atkin2015 (leafresp_method=2): PhotosynthesisMod.F90:1678
                # lmr25top = lmr_intercept_atkin(ivt) + lnc*0.2061 - 0.0402*(t10-tfrz), lnc>0.
                lmr25top = lnc_p > zero(T) ?
                    lmr_intercept_atkin[ivt_p] + lnc_p * T(0.2061) - T(0.0402) * (t10[p] - T(TFRZ)) :
                    zero(T)
            end
        else
            if c3flag_p
                lmr25top = vcmax25top * T(leaf_mr_vcm)
            else
                lmr25top = vcmax25top * T(0.025)
            end
        end

        laican = zero(T)
        for iv in 1:nrad[p]
            if iv == 1
                laican = T(0.5) * tlai_z[p, iv]
            else
                laican = laican + T(0.5) * (tlai_z[p, iv-1] + tlai_z[p, iv])
            end
            if nlevcan == 1
                nscaler = vcmaxcint[p]
            else
                nscaler = exp(-kn[p] * laican)
            end
            lmr25 = lmr25top * nscaler
            if c3flag_p
                lmr_z[p, iv] = lmr25 * ft_photo(t_veg[p], prm.lmrha) *
                               fth_photo(t_veg[p], prm.lmrhd, prm.lmrse, prm.lmrc)
            else
                lmr_z[p, iv] = lmr25 * T(2.0)^((t_veg[p] - (T(TFRZ) + T(25.0))) / T(10.0))
                lmr_z[p, iv] = lmr_z[p, iv] / (one(T) + exp(T(1.3) * (t_veg[p] - (T(TFRZ) + T(55.0)))))
            end
            if par_z_in[p, iv] <= zero(T)
                vcmax_z[p, iv] = zero(T); jmax_z_local[p, iv] = zero(T)
                tpu_z[p, iv] = zero(T); kp_z[p, iv] = zero(T)
            else
                vcmax25 = vcmax25top * nscaler
                jmax25  = jmax25top * nscaler
                tpu25   = tpu25top * nscaler
                kp25    = kp25top * nscaler
                vcmaxse = (T(668.39) - T(1.07) * smooth_clamp(t10[p] - T(TFRZ), T(11.0), T(35.0))) * prm.vcmaxse_sf
                jmaxse  = (T(659.70) - T(0.75) * smooth_clamp(t10[p] - T(TFRZ), T(11.0), T(35.0))) * prm.jmaxse_sf
                tpuse   = (T(668.39) - T(1.07) * smooth_clamp(t10[p] - T(TFRZ), T(11.0), T(35.0))) * prm.tpuse_sf
                vcmaxc = fth25_photo(prm.vcmaxhd, vcmaxse)
                jmaxc  = fth25_photo(prm.jmaxhd, jmaxse)
                tpuc   = fth25_photo(prm.tpuhd, tpuse)
                vcmax_z[p, iv] = vcmax25 * ft_photo(t_veg[p], prm.vcmaxha) *
                                 fth_photo(t_veg[p], prm.vcmaxhd, vcmaxse, vcmaxc)
                jmax_z_local[p, iv] = jmax25 * ft_photo(t_veg[p], prm.jmaxha) *
                                      fth_photo(t_veg[p], prm.jmaxhd, jmaxse, jmaxc)
                tpu_z[p, iv] = tpu25 * ft_photo(t_veg[p], prm.tpuha) *
                               fth_photo(t_veg[p], prm.tpuhd, tpuse, tpuc)
                if !c3flag_p
                    vcmax_z[p, iv] = vcmax25 * T(2.0)^((t_veg[p] - (T(TFRZ) + T(25.0))) / T(10.0))
                    vcmax_z[p, iv] = vcmax_z[p, iv] / (one(T) + exp(T(0.2) * ((T(TFRZ) + T(15.0)) - t_veg[p])))
                    vcmax_z[p, iv] = vcmax_z[p, iv] / (one(T) + exp(T(0.3) * (t_veg[p] - (T(TFRZ) + T(40.0)))))
                end
                kp_z[p, iv] = kp25 * T(2.0)^((t_veg[p] - (T(TFRZ) + T(25.0))) / T(10.0))
            end
            vcmax_z[p, iv] = vcmax_z[p, iv] * btran[p]
            lmr_z[p, iv]   = lmr_z[p, iv] * btran[p]
            if light_inhibit && par_z_in[p, 1] > zero(T)
                lmr_z[p, iv] = lmr_z[p, iv] * T(0.67)
            end
        end
    end
end

function psn_pass2_update!(ps, kn, jmax_z_local, mask_patch, ivt, slatop_pft,
        leafcn_pft, flnr_pft, fnitr_pft, dayl_factor, t10, t_veg, btran, tlai_z,
        par_z_in, vcmaxcint, nrad, use_cn::Bool, leaf_mr_vcm, nlevcan::Int,
        is_sun::Bool, overrides, bounds_patch)
    T = eltype(t_veg)
    lmrc = fth25_photo(params_inst.lmrhd, params_inst.lmrse)
    prm = (fnr = T(params_inst.fnr), act25 = T(params_inst.act25),
           tpu25ratio = T(params_inst.tpu25ratio), kp25ratio = T(params_inst.kp25ratio),
           lmrha = T(params_inst.lmrha), lmrhd = T(params_inst.lmrhd),
           lmrse = T(params_inst.lmrse), lmrc = T(lmrc),
           vcmaxse_sf = T(params_inst.vcmaxse_sf), jmaxse_sf = T(params_inst.jmaxse_sf),
           tpuse_sf = T(params_inst.tpuse_sf), vcmaxhd = T(params_inst.vcmaxhd),
           jmaxhd = T(params_inst.jmaxhd), tpuhd = T(params_inst.tpuhd),
           vcmaxha = T(params_inst.vcmaxha), jmaxha = T(params_inst.jmaxha),
           tpuha = T(params_inst.tpuha))
    jmax25top_sf_val = isnan(overrides.jmax25top_sf) ? params_inst.jmax25top_sf : overrides.jmax25top_sf
    # Atkin2015 leaf-resp base rate: per-PFT param, only indexed when
    # leafresp_method≠Ryan. Match slatop_pft's backend (host on CPU, device on GPU).
    _lmratk = T.(params_inst.lmr_intercept_atkin)
    lmr_intercept_atkin = slatop_pft isa Array ? _lmratk : convert(typeof(slatop_pft), _lmratk)
    dv = _psn_dv(ps)
    if _PSN_CI_AD_HOSTLOOP[]
        @inbounds for p in 1:length(bounds_patch)
            _psn_pass2_body!(p, dv, kn, jmax_z_local, mask_patch, ivt, slatop_pft,
                leafcn_pft, flnr_pft, fnitr_pft, lmr_intercept_atkin, dayl_factor, t10, t_veg, btran, tlai_z,
                par_z_in, vcmaxcint, nrad, prm, T(overrides.vcmax25_scale),
                T(jmax25top_sf_val), T(leaf_mr_vcm), use_cn, ps.light_inhibit,
                ps.leafresp_method, nlevcan, is_sun, LEAFRESP_MTD_RYAN1991)
        end
    else
        be = _kernel_backend(t_veg)
        _psn_pass2_kernel!(be)(dv, kn, jmax_z_local, mask_patch, ivt, slatop_pft,
            leafcn_pft, flnr_pft, fnitr_pft, lmr_intercept_atkin, dayl_factor, t10, t_veg, btran, tlai_z,
            par_z_in, vcmaxcint, nrad, prm, T(overrides.vcmax25_scale),
            T(jmax25top_sf_val), T(leaf_mr_vcm), use_cn, ps.light_inhibit,
            ps.leafresp_method, nlevcan, is_sun, LEAFRESP_MTD_RYAN1991;
            ndrange = length(bounds_patch))
        KA.synchronize(be)
    end
    return nothing
end

# =====================================================================
# PHS (plant-hydraulic-stress) kernelized passes 1/2/4. These mirror the
# single-phase Pass 1/2/4 kernels above, but the PHS path computes BOTH the
# sunlit AND shaded phase in a single body (writing the 3D ps.*_phs fields and
# the per-phase scratch). Pass 3 (the calcstress!/hybrid_PHS! Ci solve) stays a
# host loop, so the 6 psn_w*_z_{sun,sha} + bsun/bsha scratch arrays it writes are
# passed into the Pass 4 kernel as plain args (device-resident when Pass 3 is
# kernelized later). All Float64 literals → T(·); module globals → T(GLOBAL).
# =====================================================================

# PHS Pass 4: canopy integration (sun + sha) → per-patch psn/rs/lmr + btran blend.
# Internal iv loop runs serially per thread (per-patch reduction, no race).
@kernel function _psn_phs_pass4_kernel!(ps, @Const(mask_patch), @Const(nrad),
        @Const(ivt), @Const(crop_pft), @Const(lai_z_sun_in), @Const(lai_z_sha_in),
        @Const(rb), btran,
        @Const(psn_wc_z_sun), @Const(psn_wj_z_sun), @Const(psn_wp_z_sun),
        @Const(psn_wc_z_sha), @Const(psn_wj_z_sha), @Const(psn_wp_z_sha),
        @Const(bsun_arr), @Const(bsha_arr), modifyphoto_and_lmr_forcrop::Bool)
    p = @index(Global)
    _psn_phs_pass4_body!(p, ps, mask_patch, nrad, ivt, crop_pft, lai_z_sun_in,
        lai_z_sha_in, rb, btran, psn_wc_z_sun, psn_wj_z_sun, psn_wp_z_sun,
        psn_wc_z_sha, psn_wj_z_sha, psn_wp_z_sha, bsun_arr, bsha_arr,
        modifyphoto_and_lmr_forcrop)
end

# Shared per-patch body for PHS Pass 4 (KA kernel + AD host loop) — see _photosynth_ci_body!.
@inline function _psn_phs_pass4_body!(p, ps, mask_patch, nrad, ivt, crop_pft,
        lai_z_sun_in, lai_z_sha_in, rb, btran,
        psn_wc_z_sun, psn_wj_z_sun, psn_wp_z_sun,
        psn_wc_z_sha, psn_wj_z_sha, psn_wp_z_sha,
        bsun_arr, bsha_arr, modifyphoto_and_lmr_forcrop::Bool)
    @inbounds if mask_patch[p]
        T = eltype(rb)
        ivt_p = ivt[p]
        is_crop = crop_pft[ivt_p] > T(0.5)  # round(Int,·)==0 ⟺ <0.5 (0/1 valued)

        psn_z_sun  = ps.psnsun_z_patch; psn_z_sha  = ps.psnsha_z_patch
        lmr_z_sun  = ps.lmrsun_z_patch; lmr_z_sha  = ps.lmrsha_z_patch
        rs_z_sun   = ps.rssun_z_patch;  rs_z_sha   = ps.rssha_z_patch

        # Sunlit canopy
        psncan_sun = zero(T); psncan_wc_sun = zero(T); psncan_wj_sun = zero(T)
        psncan_wp_sun = zero(T); lmrcan_sun = zero(T); gscan_sun = zero(T)
        laican_sun = zero(T)
        for iv in 1:nrad[p]
            psncan_sun    += psn_z_sun[p, iv]    * lai_z_sun_in[p, iv]
            psncan_wc_sun += psn_wc_z_sun[p, iv] * lai_z_sun_in[p, iv]
            psncan_wj_sun += psn_wj_z_sun[p, iv] * lai_z_sun_in[p, iv]
            psncan_wp_sun += psn_wp_z_sun[p, iv] * lai_z_sun_in[p, iv]
            if !is_crop && modifyphoto_and_lmr_forcrop
                lmrcan_sun += lmr_z_sun[p, iv] * lai_z_sun_in[p, iv] * bsun_arr[p]
            else
                lmrcan_sun += lmr_z_sun[p, iv] * lai_z_sun_in[p, iv]
            end
            gscan_sun  += lai_z_sun_in[p, iv] / max(rb[p] + rs_z_sun[p, iv], T(1.0e-6))
            laican_sun += lai_z_sun_in[p, iv]
        end
        if laican_sun > zero(T)
            ps.psnsun_patch[p]    = psncan_sun / laican_sun
            ps.psnsun_wc_patch[p] = psncan_wc_sun / laican_sun
            ps.psnsun_wj_patch[p] = psncan_wj_sun / laican_sun
            ps.psnsun_wp_patch[p] = psncan_wp_sun / laican_sun
            ps.lmrsun_patch[p]    = lmrcan_sun / laican_sun
            ps.rssun_patch[p]     = max(laican_sun / gscan_sun - rb[p], zero(T))
        else
            ps.psnsun_patch[p] = zero(T); ps.psnsun_wc_patch[p] = zero(T)
            ps.psnsun_wj_patch[p] = zero(T); ps.psnsun_wp_patch[p] = zero(T)
            ps.lmrsun_patch[p] = zero(T); ps.rssun_patch[p] = zero(T)
        end

        # Shaded canopy
        psncan_sha = zero(T); psncan_wc_sha = zero(T); psncan_wj_sha = zero(T)
        psncan_wp_sha = zero(T); lmrcan_sha = zero(T); gscan_sha = zero(T)
        laican_sha = zero(T)
        for iv in 1:nrad[p]
            psncan_sha    += psn_z_sha[p, iv]    * lai_z_sha_in[p, iv]
            psncan_wc_sha += psn_wc_z_sha[p, iv] * lai_z_sha_in[p, iv]
            psncan_wj_sha += psn_wj_z_sha[p, iv] * lai_z_sha_in[p, iv]
            psncan_wp_sha += psn_wp_z_sha[p, iv] * lai_z_sha_in[p, iv]
            if !is_crop && modifyphoto_and_lmr_forcrop
                lmrcan_sha += lmr_z_sha[p, iv] * lai_z_sha_in[p, iv] * bsha_arr[p]
            else
                lmrcan_sha += lmr_z_sha[p, iv] * lai_z_sha_in[p, iv]
            end
            gscan_sha  += lai_z_sha_in[p, iv] / max(rb[p] + rs_z_sha[p, iv], T(1.0e-6))
            laican_sha += lai_z_sha_in[p, iv]
        end
        if laican_sha > zero(T)
            ps.psnsha_patch[p]    = psncan_sha / laican_sha
            ps.psnsha_wc_patch[p] = psncan_wc_sha / laican_sha
            ps.psnsha_wj_patch[p] = psncan_wj_sha / laican_sha
            ps.psnsha_wp_patch[p] = psncan_wp_sha / laican_sha
            ps.lmrsha_patch[p]    = lmrcan_sha / laican_sha
            ps.rssha_patch[p]     = max(laican_sha / gscan_sha - rb[p], zero(T))
        else
            ps.psnsha_patch[p] = zero(T); ps.psnsha_wc_patch[p] = zero(T)
            ps.psnsha_wj_patch[p] = zero(T); ps.psnsha_wp_patch[p] = zero(T)
            ps.lmrsha_patch[p] = zero(T); ps.rssha_patch[p] = zero(T)
        end

        # btran from LAI-weighted bsun/bsha (btran is mutated in place)
        if laican_sha + laican_sun > zero(T)
            btran[p] = bsun_arr[p] * (laican_sun / (laican_sun + laican_sha)) +
                       bsha_arr[p] * (laican_sha / (laican_sun + laican_sha))
        else
            btran[p] = bsun_arr[p]
        end
    end
end

function psn_phs_pass4_update!(ps, mask_patch, nrad, ivt, crop_pft,
        lai_z_sun_in, lai_z_sha_in, rb, btran,
        psn_wc_z_sun, psn_wj_z_sun, psn_wp_z_sun,
        psn_wc_z_sha, psn_wj_z_sha, psn_wp_z_sha,
        bsun_arr, bsha_arr, modifyphoto_and_lmr_forcrop::Bool, bounds_patch)
    dv = _psn_dv(ps)
    if _PSN_CI_AD_HOSTLOOP[]
        @inbounds for p in 1:length(bounds_patch)
            _psn_phs_pass4_body!(p, dv, mask_patch, nrad, ivt, crop_pft,
                lai_z_sun_in, lai_z_sha_in, rb, btran,
                psn_wc_z_sun, psn_wj_z_sun, psn_wp_z_sun,
                psn_wc_z_sha, psn_wj_z_sha, psn_wp_z_sha,
                bsun_arr, bsha_arr, modifyphoto_and_lmr_forcrop)
        end
    else
        be = _kernel_backend(rb)
        _psn_phs_pass4_kernel!(be)(dv, mask_patch, nrad, ivt, crop_pft,
            lai_z_sun_in, lai_z_sha_in, rb, btran,
            psn_wc_z_sun, psn_wj_z_sun, psn_wp_z_sun,
            psn_wc_z_sha, psn_wj_z_sha, psn_wp_z_sha,
            bsun_arr, bsha_arr, modifyphoto_and_lmr_forcrop;
            ndrange = length(bounds_patch))
        KA.synchronize(be)
    end
    return nothing
end

# PHS Pass 1: root-soil conductance k_soil_root (per-patch, internal serial j loop
# over soil levels). Vulnerability curve via _plc(smp, psi50, ck) from the threaded
# phs_params bundle (no error()/dispatch); krmax from phs_params too.
@kernel function _psn_phs_pass1_kernel!(@Const(mask_patch), @Const(col_of_patch),
        @Const(ivt), phs_params, @Const(froot_carbon), @Const(rootfr), @Const(dz),
        @Const(tsai), @Const(tlai), @Const(froot_leaf_pft), @Const(root_radius_pft),
        @Const(root_density_pft), @Const(hksat), @Const(hk_l), @Const(smp_l),
        @Const(z_col), k_soil_root, root_conductance_out, soil_conductance_out,
        c_to_b, croot_lateral_length, nlevsoi::Int, root_seg::Int)
    p = @index(Global)
    _psn_phs_pass1_body!(p, mask_patch, col_of_patch, ivt, phs_params, froot_carbon,
        rootfr, dz, tsai, tlai, froot_leaf_pft, root_radius_pft, root_density_pft,
        hksat, hk_l, smp_l, z_col, k_soil_root, root_conductance_out,
        soil_conductance_out, c_to_b, croot_lateral_length, nlevsoi, root_seg)
end

# Shared per-patch body for PHS Pass 1 (KA kernel + AD host loop) — see _photosynth_ci_body!.
@inline function _psn_phs_pass1_body!(p, mask_patch, col_of_patch, ivt, phs_params,
        froot_carbon, rootfr, dz, tsai, tlai, froot_leaf_pft, root_radius_pft,
        root_density_pft, hksat, hk_l, smp_l, z_col, k_soil_root, root_conductance_out,
        soil_conductance_out, c_to_b, croot_lateral_length, nlevsoi::Int, root_seg::Int)
    @inbounds if mask_patch[p]
        T = eltype(smp_l)
        c = col_of_patch[p]
        ivt_p = ivt[p]
        for j in 1:nlevsoi
            root_biomass_density = c_to_b * froot_carbon[p] * rootfr[p, j] / dz[c, j]
            root_biomass_density = max(c_to_b * one(T), root_biomass_density)

            root_cross_sec_area = T(RPI) * root_radius_pft[ivt_p]^2
            root_length_density = root_biomass_density / (root_density_pft[ivt_p] * root_cross_sec_area)

            rai_j = (tsai[p] + tlai[p]) * froot_leaf_pft[ivt_p] * rootfr[p, j]

            croot_average_length = croot_lateral_length
            r_soil = sqrt(one(T) / (T(RPI) * root_length_density))
            soil_cond = min(hksat[c, j], hk_l[c, j]) / (T(1.0e3) * r_soil)

            fs_j = _plc(smp_l[c, j], phs_params.psi50[ivt_p, root_seg],
                        phs_params.ck[ivt_p, root_seg])
            root_cond = (fs_j * rai_j * phs_params.krmax[ivt_p]) /
                        (croot_average_length + z_col[c, j])

            soil_cond = max(soil_cond, T(1.0e-16))   # hard: constant floor, zero derivative on the clamped branch. At k=50 this 1e-16 floor was silently raised to log(2)/50 = 0.0139 [s-1], which is ~13 orders above the floor and inflates even WET-soil conductance (~0.018) by ~38%.
            root_cond = max(root_cond, T(1.0e-16))   # hard: see soil_cond above.

            root_conductance_out[p, j] = root_cond
            soil_conductance_out[p, j] = soil_cond

            rs_resis = one(T) / soil_cond + one(T) / root_cond

            if rai_j * rootfr[p, j] > zero(T) && j > 1
                k_soil_root[p, j] = one(T) / rs_resis
            else
                k_soil_root[p, j] = zero(T)
            end
        end
    end
end

function psn_phs_pass1_update!(ps, mask_patch, col_of_patch, ivt, froot_carbon,
        rootfr, dz, tsai, tlai, froot_leaf_pft, root_radius_pft, root_density_pft,
        hksat, hk_l, smp_l, z_col, k_soil_root, root_conductance_out,
        soil_conductance_out, c_to_b, croot_lateral_length, nlevsoi::Int, bounds_patch)
    T = eltype(smp_l)
    phs_params = _psn_phs_params(smp_l)
    if _PSN_CI_AD_HOSTLOOP[]
        # AD-mode host loop (see compositional_reverse! / _photosynth_ci_body!):
        # byte-identical to the KA kernel (shared _psn_phs_pass1_body!), avoids the
        # KA-launch codegen Enzyme reverse segfaults on under Julia 1.12.
        @inbounds for p in 1:length(bounds_patch)
            _psn_phs_pass1_body!(p, mask_patch, col_of_patch, ivt, phs_params,
                froot_carbon, rootfr, dz, tsai, tlai, froot_leaf_pft, root_radius_pft,
                root_density_pft, hksat, hk_l, smp_l, z_col, k_soil_root,
                root_conductance_out, soil_conductance_out, T(c_to_b),
                T(croot_lateral_length), nlevsoi, ROOT_SEG)
        end
    else
        be = _kernel_backend(smp_l)
        # host-global pftcon params -> device backend/precision (no-op on host)
        _FTp = eltype(smp_l)
        _frl = _to_backend_like(smp_l, _FTp, froot_leaf_pft)
        _rr  = _to_backend_like(smp_l, _FTp, root_radius_pft)
        _rd  = _to_backend_like(smp_l, _FTp, root_density_pft)
        _psn_phs_pass1_kernel!(be)(mask_patch, col_of_patch, ivt, phs_params,
            froot_carbon, rootfr, dz, tsai, tlai, _frl, _rr,
            _rd, hksat, hk_l, smp_l, z_col, k_soil_root,
            root_conductance_out, soil_conductance_out, T(c_to_b),
            T(croot_lateral_length), nlevsoi, ROOT_SEG; ndrange = length(bounds_patch))
        KA.synchronize(be)
    end
    return nothing
end

# --- PHS Pass-2 device-view arg bundles (Metal 31-arg-buffer reduction) ---
# Loose @Const input arrays grouped into immutable @adapt_structure'd bundles so the
# kernel takes a handful of struct args (each = 1 indirect-arg-buffer). Fields are
# aliased to locals at the kernel top so the physics body is UNCHANGED.
# Index/Bool arrays (VB=mask Bool, VI=ivt/nrad Int). PFT params are Float64-pinned
# in the AD path → own struct/param (Vp) so they don't unify with Dual state.
Base.@kwdef struct _Psn2Idx{VB,VI}
    mask_patch::VB; ivt::VI; nrad::VI
end
Base.@kwdef struct _Psn2Pft{Vp}
    c3psn_pft::Vp; mbbopt_pft::Vp; slatop_pft::Vp; leafcn_pft::Vp
    flnr_pft::Vp; fnitr_pft::Vp; lmr_intercept_atkin::Vp
end
Base.@kwdef struct _Psn2In{V,M}     # read-only forcing vectors (V) + matrices (M)
    forc_pbot::V; oair::V; dayl_factor::V; t10::V; t_veg::V
    vcmaxcint_sun::V; vcmaxcint_sha::V
    tlai_z::M; par_z_sun_in::M; par_z_sha_in::M
end
Adapt.@adapt_structure _Psn2Idx
Adapt.@adapt_structure _Psn2Pft
Adapt.@adapt_structure _Psn2In
# Loose scalars bundled into one isbits struct (Float scalars share S; Bools/Ints own
# fields). Counts as 1 kernel arg instead of 15.
struct _Psn2Scalars{S}
    vcmax25_scale::S; jmax25top_sf_val::S; leaf_mr_vcm::S; bbbopt_c3::S; bbbopt_c4::S
    use_cn::Bool; use_c13::Bool; light_inhibit::Bool; use_luna::Bool
    leafresp_method::Int; nlevcan::Int; stomatalcond_mtd::Int
    stomatal_bb::Int; leafresp_ryan::Int; sun::Int; sha::Int
end

# PHS Pass 2: kinetics + N-profile + vcmax/jmax/tpu/kp/lmr, BOTH phases, writing
# the 3D ps.*_phs fields ([p,SUN,iv] & [p,SHA,iv]) in one iv pass. params_inst
# scalars bundled in `prm`; overrides / leaf_mr_vcm resolved on the host.
@kernel function _psn_phs_pass2_kernel!(ps, kn, jmax_z_local, prm, idx, pft, inp, sc2)
    p = @index(Global)
    _psn_phs_pass2_body!(p, ps, kn, jmax_z_local, prm, idx, pft, inp, sc2)
end

# Shared per-patch body for PHS Pass 2 (KA kernel + AD host loop) — see _photosynth_ci_body!.
@inline function _psn_phs_pass2_body!(p, ps, kn, jmax_z_local, prm, idx, pft, inp, sc2)
    # alias grouped fields to locals → physics body below is byte-identical
    mask_patch = idx.mask_patch; ivt = idx.ivt; nrad = idx.nrad
    c3psn_pft = pft.c3psn_pft; mbbopt_pft = pft.mbbopt_pft; slatop_pft = pft.slatop_pft
    leafcn_pft = pft.leafcn_pft; flnr_pft = pft.flnr_pft; fnitr_pft = pft.fnitr_pft
    lmr_intercept_atkin = pft.lmr_intercept_atkin
    forc_pbot = inp.forc_pbot; oair = inp.oair; dayl_factor = inp.dayl_factor
    t10 = inp.t10; t_veg = inp.t_veg
    vcmaxcint_sun = inp.vcmaxcint_sun; vcmaxcint_sha = inp.vcmaxcint_sha
    tlai_z = inp.tlai_z; par_z_sun_in = inp.par_z_sun_in; par_z_sha_in = inp.par_z_sha_in
    vcmax25_scale = sc2.vcmax25_scale; jmax25top_sf_val = sc2.jmax25top_sf_val
    leaf_mr_vcm = sc2.leaf_mr_vcm; bbbopt_c3 = sc2.bbbopt_c3; bbbopt_c4 = sc2.bbbopt_c4
    use_cn = sc2.use_cn; use_c13 = sc2.use_c13; light_inhibit = sc2.light_inhibit
    use_luna = sc2.use_luna
    leafresp_method = sc2.leafresp_method; nlevcan = sc2.nlevcan
    stomatalcond_mtd = sc2.stomatalcond_mtd; stomatal_bb = sc2.stomatal_bb
    leafresp_ryan = sc2.leafresp_ryan; sun = sc2.sun; sha = sc2.sha
    vcmx25_z = ps.vcmx25_z_patch; jmx25_z = ps.jmx25_z_patch
    @inbounds if mask_patch[p]
        T = eltype(t_veg)
        ivt_p = ivt[p]
        vcmax_z = ps.vcmax_z_phs_patch; tpu_z = ps.tpu_z_phs_patch
        kp_z = ps.kp_z_phs_patch
        lmr_z_sun = ps.lmrsun_z_patch; lmr_z_sha = ps.lmrsha_z_patch

        # C3/C4 flag (round(Int,·)==1 ⟺ >0.5; 0/1 valued — heap-free on Metal)
        if c3psn_pft[ivt_p] > T(0.5)
            ps.c3flag_patch[p] = true
        else
            ps.c3flag_patch[p] = false
        end
        c3flag_p = ps.c3flag_patch[p]

        if c3flag_p
            ps.qe_patch[p] = zero(T)
            bbbopt_p = bbbopt_c3
        else
            ps.qe_patch[p] = T(0.05)
            bbbopt_p = bbbopt_c4
        end

        if stomatalcond_mtd == stomatal_bb
            ps.bbb_patch[p] = bbbopt_p
            ps.mbb_patch[p] = mbbopt_pft[ivt_p]
        end

        kc25 = prm.kc25_coef * forc_pbot[p]
        ko25 = prm.ko25_coef * forc_pbot[p]
        sco  = T(0.5) * T(0.209) / prm.cp25_yr2000
        cp25 = T(0.5) * oair[p] / sco
        ps.kc_patch[p] = kc25 * ft_photo(t_veg[p], prm.kcha)
        ps.ko_patch[p] = ko25 * ft_photo(t_veg[p], prm.koha)
        ps.cp_patch[p] = cp25 * ft_photo(t_veg[p], prm.cpha)

        # Nitrogen profile
        lnc_p = one(T) / (slatop_pft[ivt_p] * leafcn_pft[ivt_p])
        lnc_p = smooth_min(lnc_p, T(10.0))
        ps.lnca_patch[p] = lnc_p

        vcmax25top = lnc_p * flnr_pft[ivt_p] * prm.fnr * prm.act25 * dayl_factor[p]
        if !use_cn
            vcmax25top = vcmax25top * fnitr_pft[ivt_p]
        end
        if !isnan(vcmax25_scale)
            vcmax25top = vcmax25top * vcmax25_scale
        end

        jmax25top = ((T(2.59) - T(0.035) * smooth_clamp(t10[p] - T(TFRZ), T(11.0), T(35.0))) *
                     vcmax25top) * jmax25top_sf_val
        tpu25top = prm.tpu25ratio * vcmax25top
        kp25top  = prm.kp25ratio * vcmax25top

        ps.luvcmax25top_patch[p] = vcmax25top
        ps.lujmax25top_patch[p] = jmax25top
        ps.lutpu25top_patch[p] = tpu25top

        if dayl_factor[p] < T(1.0e-12)
            kn[p] = zero(T)
        else
            kn[p] = exp(T(0.00963) * vcmax25top / dayl_factor[p] - T(2.43))
        end

        if use_cn
            if leafresp_method == leafresp_ryan
                lmr25top = T(2.525e-6) * (T(1.5)^((T(25.0) - T(20.0)) / T(10.0)))
                lmr25top = lmr25top * lnc_p / T(12.0e-06)
            else
                # Atkin2015 (leafresp_method=2): PhotosynthesisMod.F90:3289 (PHS path)
                lmr25top = lnc_p > zero(T) ?
                    lmr_intercept_atkin[ivt_p] + lnc_p * T(0.2061) - T(0.0402) * (t10[p] - T(TFRZ)) :
                    zero(T)
            end
        else
            if c3flag_p
                lmr25top = vcmax25top * T(leaf_mr_vcm)
            else
                lmr25top = vcmax25top * T(0.025)
            end
        end

        laican = zero(T)
        for iv in 1:nrad[p]
            if iv == 1
                laican = T(0.5) * tlai_z[p, iv]
            else
                laican = laican + T(0.5) * (tlai_z[p, iv-1] + tlai_z[p, iv])
            end

            if nlevcan == 1
                nscaler_sun = vcmaxcint_sun[p]
                nscaler_sha = vcmaxcint_sha[p]
            else
                nscaler_sun = exp(-kn[p] * laican)
                nscaler_sha = exp(-kn[p] * laican)
            end

            lmr25_sun = lmr25top * nscaler_sun
            lmr25_sha = lmr25top * nscaler_sha

            # LUNA leaf-respiration base rate: when LUNA is active without prognostic
            # CN, leaf maintenance respiration scales with the LUNA-acclimated vcmax25
            # (vcmx25_z) rather than the static lnc-based vcmax25top*nscaler. Matches
            # PhotosynthesisMod.F90:3335-3340 (use_luna .and. c3flag .and. crop==0 .and.
            # .not.use_cn). Mirrors the existing LUNA vcmax25 branch below (crop check
            # omitted, consistent with that branch).
            if use_luna && c3flag_p && !use_cn
                lmr25_sun = T(leaf_mr_vcm) * vcmx25_z[p, iv]
                lmr25_sha = T(leaf_mr_vcm) * vcmx25_z[p, iv]
            end

            if c3flag_p
                lmr_z_sun[p, iv] = lmr25_sun * ft_photo(t_veg[p], prm.lmrha) *
                                   fth_photo(t_veg[p], prm.lmrhd, prm.lmrse, prm.lmrc)
                lmr_z_sha[p, iv] = lmr25_sha * ft_photo(t_veg[p], prm.lmrha) *
                                   fth_photo(t_veg[p], prm.lmrhd, prm.lmrse, prm.lmrc)
            else
                lmr_z_sun[p, iv] = lmr25_sun * T(2.0)^((t_veg[p] - (T(TFRZ) + T(25.0))) / T(10.0))
                lmr_z_sun[p, iv] = lmr_z_sun[p, iv] / (one(T) + exp(T(1.3) * (t_veg[p] - (T(TFRZ) + T(55.0)))))
                lmr_z_sha[p, iv] = lmr25_sha * T(2.0)^((t_veg[p] - (T(TFRZ) + T(25.0))) / T(10.0))
                lmr_z_sha[p, iv] = lmr_z_sha[p, iv] / (one(T) + exp(T(1.3) * (t_veg[p] - (T(TFRZ) + T(55.0)))))
            end

            # Reduce lmr with low LAI
            lmr_z_sun[p, iv] *= smooth_min(T(0.2) * exp(T(3.218) * tlai_z[p, iv]), one(T))
            lmr_z_sha[p, iv] *= smooth_min(T(0.2) * exp(T(3.218) * tlai_z[p, iv]), one(T))

            if par_z_sun_in[p, iv] <= zero(T)  # night time
                vcmax_z[p, sun, iv] = zero(T); jmax_z_local[p, sun, iv] = zero(T)
                tpu_z[p, sun, iv] = zero(T); kp_z[p, sun, iv] = zero(T)
                vcmax_z[p, sha, iv] = zero(T); jmax_z_local[p, sha, iv] = zero(T)
                tpu_z[p, sha, iv] = zero(T); kp_z[p, sha, iv] = zero(T)

                if use_c13
                    ps.alphapsnsun_patch[p] = one(T)
                    ps.alphapsnsha_patch[p] = one(T)
                end
            else  # day time
                if use_luna && c3flag_p
                    # LUNA-acclimated vcmax25/jmax25 (injected from restart), overriding
                    # the static lnc-based vcmax25top. Sun uses the layer value directly;
                    # sha is scaled by the canopy sun/sha integration ratio. Matches
                    # PhotosynthesisMod.F90:3378 (use_luna .and. c3flag .and. crop==0).
                    vcmax25_sun = vcmx25_z[p, iv]
                    jmax25_sun = jmx25_z[p, iv]
                    tpu25_sun = prm.tpu25ratio * vcmax25_sun
                    vcmax25_sha = vcmax25_sun
                    jmax25_sha = jmax25_sun
                    tpu25_sha = prm.tpu25ratio * vcmax25_sha
                    if vcmaxcint_sun[p] > zero(T) && nlevcan == 1
                        _luna_r = vcmaxcint_sha[p] / vcmaxcint_sun[p]
                        vcmax25_sha = vcmax25_sun * _luna_r
                        jmax25_sha = jmax25_sun * _luna_r
                        tpu25_sha = tpu25_sun * _luna_r
                    end
                else
                    vcmax25_sun = vcmax25top * nscaler_sun
                    jmax25_sun = jmax25top * nscaler_sun
                    tpu25_sun = tpu25top * nscaler_sun
                    vcmax25_sha = vcmax25top * nscaler_sha
                    jmax25_sha = jmax25top * nscaler_sha
                    tpu25_sha = tpu25top * nscaler_sha
                end
                kp25_sun = kp25top * nscaler_sun
                kp25_sha = kp25top * nscaler_sha

                vcmaxse = (T(668.39) - T(1.07) * smooth_clamp(t10[p] - T(TFRZ), T(11.0), T(35.0))) * prm.vcmaxse_sf
                jmaxse  = (T(659.70) - T(0.75) * smooth_clamp(t10[p] - T(TFRZ), T(11.0), T(35.0))) * prm.jmaxse_sf
                tpuse   = (T(668.39) - T(1.07) * smooth_clamp(t10[p] - T(TFRZ), T(11.0), T(35.0))) * prm.tpuse_sf
                vcmaxc = fth25_photo(prm.vcmaxhd, vcmaxse)
                jmaxc  = fth25_photo(prm.jmaxhd, jmaxse)
                tpuc   = fth25_photo(prm.tpuhd, tpuse)

                vcmax_z[p, sun, iv] = vcmax25_sun * ft_photo(t_veg[p], prm.vcmaxha) *
                                      fth_photo(t_veg[p], prm.vcmaxhd, vcmaxse, vcmaxc)
                jmax_z_local[p, sun, iv] = jmax25_sun * ft_photo(t_veg[p], prm.jmaxha) *
                                            fth_photo(t_veg[p], prm.jmaxhd, jmaxse, jmaxc)
                tpu_z[p, sun, iv] = tpu25_sun * ft_photo(t_veg[p], prm.tpuha) *
                                    fth_photo(t_veg[p], prm.tpuhd, tpuse, tpuc)

                vcmax_z[p, sha, iv] = vcmax25_sha * ft_photo(t_veg[p], prm.vcmaxha) *
                                      fth_photo(t_veg[p], prm.vcmaxhd, vcmaxse, vcmaxc)
                jmax_z_local[p, sha, iv] = jmax25_sha * ft_photo(t_veg[p], prm.jmaxha) *
                                            fth_photo(t_veg[p], prm.jmaxhd, jmaxse, jmaxc)
                tpu_z[p, sha, iv] = tpu25_sha * ft_photo(t_veg[p], prm.tpuha) *
                                    fth_photo(t_veg[p], prm.tpuhd, tpuse, tpuc)

                if !c3flag_p
                    vcmax_z[p, sun, iv] = vcmax25_sun * T(2.0)^((t_veg[p] - (T(TFRZ) + T(25.0))) / T(10.0))
                    vcmax_z[p, sun, iv] /= (one(T) + exp(T(0.2) * ((T(TFRZ) + T(15.0)) - t_veg[p])))
                    vcmax_z[p, sun, iv] /= (one(T) + exp(T(0.3) * (t_veg[p] - (T(TFRZ) + T(40.0)))))
                    vcmax_z[p, sha, iv] = vcmax25_sha * T(2.0)^((t_veg[p] - (T(TFRZ) + T(25.0))) / T(10.0))
                    vcmax_z[p, sha, iv] /= (one(T) + exp(T(0.2) * ((T(TFRZ) + T(15.0)) - t_veg[p])))
                    vcmax_z[p, sha, iv] /= (one(T) + exp(T(0.3) * (t_veg[p] - (T(TFRZ) + T(40.0)))))
                end

                kp_z[p, sun, iv] = kp25_sun * T(2.0)^((t_veg[p] - (T(TFRZ) + T(25.0))) / T(10.0))
                kp_z[p, sha, iv] = kp25_sha * T(2.0)^((t_veg[p] - (T(TFRZ) + T(25.0))) / T(10.0))
            end

            # Light inhibition
            if light_inhibit && par_z_sun_in[p, 1] > zero(T)
                lmr_z_sun[p, iv] *= T(0.67)
            end
            if light_inhibit && par_z_sha_in[p, 1] > zero(T)
                lmr_z_sha[p, iv] *= T(0.67)
            end
        end
    end
end

function psn_phs_pass2_update!(ps, kn, jmax_z_local, mask_patch, ivt, c3psn_pft,
        mbbopt_pft, forc_pbot, oair, slatop_pft, leafcn_pft, flnr_pft, fnitr_pft,
        dayl_factor, t10, t_veg, tlai_z, par_z_sun_in, par_z_sha_in, vcmaxcint_sun,
        vcmaxcint_sha, nrad, use_cn::Bool, use_c13::Bool, leaf_mr_vcm, nlevcan::Int,
        stomatalcond_mtd::Int, light_inhibit::Bool, leafresp_method::Int,
        overrides, bounds_patch; use_luna::Bool=false)
    T = eltype(t_veg)
    lmrc = fth25_photo(params_inst.lmrhd, params_inst.lmrse)
    prm = (kc25_coef = T(params_inst.kc25_coef), ko25_coef = T(params_inst.ko25_coef),
           cp25_yr2000 = T(params_inst.cp25_yr2000), kcha = T(params_inst.kcha),
           koha = T(params_inst.koha), cpha = T(params_inst.cpha),
           fnr = T(params_inst.fnr), act25 = T(params_inst.act25),
           tpu25ratio = T(params_inst.tpu25ratio), kp25ratio = T(params_inst.kp25ratio),
           lmrha = T(params_inst.lmrha), lmrhd = T(params_inst.lmrhd),
           lmrse = T(params_inst.lmrse), lmrc = T(lmrc),
           vcmaxse_sf = T(params_inst.vcmaxse_sf), jmaxse_sf = T(params_inst.jmaxse_sf),
           tpuse_sf = T(params_inst.tpuse_sf), vcmaxhd = T(params_inst.vcmaxhd),
           jmaxhd = T(params_inst.jmaxhd), tpuhd = T(params_inst.tpuhd),
           vcmaxha = T(params_inst.vcmaxha), jmaxha = T(params_inst.jmaxha),
           tpuha = T(params_inst.tpuha))
    jmax25top_sf_val = isnan(overrides.jmax25top_sf) ? params_inst.jmax25top_sf : overrides.jmax25top_sf
    dv = _psn_dv(ps)
    # group loose @Const arrays into device-view bundles + scalars into one isbits
    # struct (Metal 31-arg-buffer reduction); same refs, so behaviour unchanged.
    idx = _Psn2Idx(; mask_patch = mask_patch, ivt = ivt, nrad = nrad)
    _lmratk = T.(params_inst.lmr_intercept_atkin)
    lmr_intercept_atkin = slatop_pft isa Array ? _lmratk : convert(typeof(slatop_pft), _lmratk)
    pft = _Psn2Pft(; c3psn_pft = c3psn_pft, mbbopt_pft = mbbopt_pft,
        slatop_pft = slatop_pft, leafcn_pft = leafcn_pft, flnr_pft = flnr_pft,
        fnitr_pft = fnitr_pft, lmr_intercept_atkin = lmr_intercept_atkin)
    inp = _Psn2In(; forc_pbot = forc_pbot, oair = oair, dayl_factor = dayl_factor,
        t10 = t10, t_veg = t_veg, vcmaxcint_sun = vcmaxcint_sun,
        vcmaxcint_sha = vcmaxcint_sha, tlai_z = tlai_z,
        par_z_sun_in = par_z_sun_in, par_z_sha_in = par_z_sha_in)
    sc2 = _Psn2Scalars{T}(T(overrides.vcmax25_scale), T(jmax25top_sf_val),
        T(leaf_mr_vcm), T(BBBOPT_C3), T(BBBOPT_C4),
        use_cn, use_c13, light_inhibit, use_luna, leafresp_method, nlevcan,
        stomatalcond_mtd, STOMATALCOND_MTD_BB1987, LEAFRESP_MTD_RYAN1991, SUN, SHA)
    if _PSN_CI_AD_HOSTLOOP[]
        @inbounds for p in 1:length(bounds_patch)
            _psn_phs_pass2_body!(p, dv, kn, jmax_z_local, prm, idx, pft, inp, sc2)
        end
    else
        be = _kernel_backend(t_veg)
        _psn_phs_pass2_kernel!(be)(dv, kn, jmax_z_local, prm, idx, pft, inp, sc2;
            ndrange = length(bounds_patch))
        KA.synchronize(be)
    end
    return nothing
end

# =====================================================================
# photosynthesis! — Main leaf photosynthesis and stomatal conductance
# =====================================================================

"""
    photosynthesis!(ps, ...)

Leaf photosynthesis and stomatal conductance calculation.
Ported from `subroutine Photosynthesis` in `PhotosynthesisMod.F90`.

This is a simplified interface for unit testing. In production CLM,
additional data types (atm2lnd, temperature, surfalb, solarabs, canopystate,
ozone) would be passed. Here we pass the needed arrays directly.
"""
function photosynthesis!(ps,
                         esat_tv::AbstractVector{<:Real},
                         eair::AbstractVector{<:Real},
                         oair::AbstractVector{<:Real},
                         cair::AbstractVector{<:Real},
                         rb::AbstractVector{<:Real},
                         btran::AbstractVector{<:Real},
                         dayl_factor::AbstractVector{<:Real},
                         leafn::AbstractVector{<:Real},
                         forc_pbot::AbstractVector{<:Real},
                         t_veg::AbstractVector{<:Real},
                         t10::AbstractVector{<:Real},
                         tgcm::AbstractVector{<:Real},
                         nrad::AbstractVector{<:Integer},
                         tlai_z::AbstractMatrix{<:Real},
                         tlai::AbstractVector{<:Real},
                         par_z_in::AbstractMatrix{<:Real},
                         lai_z_in::AbstractMatrix{<:Real},
                         vcmaxcint::AbstractVector{<:Real},
                         o3coefv::AbstractVector{<:Real},
                         o3coefg::AbstractVector{<:Real},
                         c3psn_pft::AbstractVector{<:Real},
                         leafcn_pft::AbstractVector{<:Real},
                         flnr_pft::AbstractVector{<:Real},
                         fnitr_pft::AbstractVector{<:Real},
                         slatop_pft::AbstractVector{<:Real},
                         mbbopt_pft::AbstractVector{<:Real},
                         medlynintercept_pft::AbstractVector{<:Real},
                         medlynslope_pft::AbstractVector{<:Real},
                         ivt::AbstractVector{<:Integer},
                         col_of_patch::AbstractVector{<:Integer},
                         mask_patch::AbstractVector{Bool},
                         bounds_patch::UnitRange{Int},
                         phase::String;
                         nlevcan::Int=NLEVCAN,
                         use_cn::Bool=false,
                         use_luna::Bool=false,
                         use_c13::Bool=false,
                         leaf_mr_vcm::Real=0.015,
                         crop_pft::AbstractVector{<:Real}=Float64[],
                         # Calibration overrides (NaN = use defaults)
                         overrides::CalibrationOverrides=CalibrationOverrides())
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

    np = length(bounds_patch)
    FT = eltype(t_veg)
    # Scratch via similar() of an input (device-resident on GPU) + fill!; the _z
    # matrices feed the kernelized passes. lmrc/bbbopt are no longer needed here —
    # the Pass 1/2 kernels resolve them internally.
    jmax_z_local = fill!(similar(t_veg, FT, np, nlevcan), zero(FT))
    kn           = fill!(similar(t_veg, FT, np), zero(FT))
    psn_wc_z     = fill!(similar(t_veg, FT, np, nlevcan), zero(FT))
    psn_wj_z     = fill!(similar(t_veg, FT, np, nlevcan), zero(FT))
    psn_wp_z     = fill!(similar(t_veg, FT, np, nlevcan), zero(FT))

    rsmax0 = 2.0e4

    # ---- Pass 1: Compute kinetics, N scaling, respiration (kernelized) ----
    # bbbopt is now a kernel-local scalar (was a per-patch scratch array used only
    # within this pass). Writes ps.c3flag/qe/bbb/mbb/kc/ko/cp (== local aliases).
    psn_pass1_update!(ps, mask_patch, ivt, c3psn_pft, mbbopt_pft,
        forc_pbot, oair, t_veg, btran, stomatalcond_mtd, bounds_patch)

    # ---- Pass 2: Nitrogen profile and vcmax (kernelized) ----
    psn_pass2_update!(ps, kn, jmax_z_local, mask_patch, ivt, slatop_pft,
        leafcn_pft, flnr_pft, fnitr_pft, dayl_factor, t10, t_veg, btran, tlai_z,
        par_z_in, vcmaxcint, nrad, use_cn, leaf_mr_vcm, nlevcan,
        phase == "sun", overrides, bounds_patch)

    # ---- Pass 3: Leaf-level photosynthesis and stomatal conductance ----
    # Kernelized: one thread per patch, with the per-canopy-layer iv loop running
    # serially inside each thread (per-patch outputs gb_mol/vpd_can/rh_leaf are
    # written within the iv loop, so a (p,iv) flatten would race). See
    # _photosynth_ci_kernel! / photosynth_ci_solve! above.
    photosynth_ci_solve!(ps,
        ac, aj, ap, ag, an, psn_z, psn_wc_z, psn_wj_z, psn_wp_z,
        rs_z, ci_z, gs_mol, rh_leaf, gb_mol, vpd_can,
        mask_patch, ivt,
        medlynslope_pft, medlynintercept_pft, mbbopt_pft,
        forc_pbot, tgcm, rb, par_z_in, lmr_z, jmax_z_local, c3flag,
        eair, esat_tv, cair, oair, bbb, o3coefg, o3coefv, nrad,
        bounds_patch, phase, stomatalcond_mtd, rsmax0, overrides)

    # ---- Pass 4: Canopy integration (kernelized) ----
    psn_pass4_update!(ps, mask_patch, nrad, lai_z_in, rb,
        psn_wc_z, psn_wj_z, psn_wp_z, phase == "sun", bounds_patch)

    return nothing
end

# =====================================================================
# photosynthesis_total! — Determine total photosynthesis
# =====================================================================

# Per-patch combine of sunlit + shaded fluxes (mask-based standalone form; the
# canopy_fluxes! call site uses the filterp-gather variant _cf_psn_total_kernel!).
@kernel function _psn_total_kernel!(fpsn, fpsn_wc, fpsn_wj, fpsn_wp,
        @Const(mask_patch), @Const(psnsun), @Const(psnsun_wc), @Const(psnsun_wj),
        @Const(psnsun_wp), @Const(psnsha), @Const(psnsha_wc), @Const(psnsha_wj),
        @Const(psnsha_wp), @Const(laisun), @Const(laisha))
    p = @index(Global)
    @inbounds if mask_patch[p]
        ls = laisun[p]; lh = laisha[p]
        fpsn[p]    = psnsun[p]    * ls + psnsha[p]    * lh
        fpsn_wc[p] = psnsun_wc[p] * ls + psnsha_wc[p] * lh
        fpsn_wj[p] = psnsun_wj[p] * ls + psnsha_wj[p] * lh
        fpsn_wp[p] = psnsun_wp[p] * ls + psnsha_wp[p] * lh
    end
end

"""
    photosynthesis_total!(ps, laisun, laisha, mask_patch, bounds_patch;
                          use_fates=false, use_cn=false, use_c13=false, use_c14=false)

Determine total photosynthesis by combining sunlit and shaded fluxes.
Ported from `subroutine PhotosynthesisTotal`.
"""
function photosynthesis_total!(ps,
                               laisun::AbstractVector{<:Real},
                               laisha::AbstractVector{<:Real},
                               mask_patch::AbstractVector{Bool},
                               bounds_patch::UnitRange{Int};
                               use_fates::Bool=false,
                               use_cn::Bool=false,
                               use_c13::Bool=false,
                               use_c14::Bool=false)
    if !use_fates
        _launch!(_psn_total_kernel!, ps.fpsn_patch, ps.fpsn_wc_patch, ps.fpsn_wj_patch,
            ps.fpsn_wp_patch, mask_patch, ps.psnsun_patch, ps.psnsun_wc_patch,
            ps.psnsun_wj_patch, ps.psnsun_wp_patch, ps.psnsha_patch, ps.psnsha_wc_patch,
            ps.psnsha_wj_patch, ps.psnsha_wp_patch, laisun, laisha;
            ndrange = length(bounds_patch))
    end

    return nothing
end

# =====================================================================
# fractionation! — C13 fractionation during photosynthesis
# =====================================================================

# Per-patch C13 fractionation. alphapsn/gs_mol_ref/an_ref are resolved on the host
# (phase + use_hydrstress) then passed loose; the per-canopy-layer iv loop runs
# serially in-thread (alphapsn[p] is an own-index write — last iv wins, no race).
# col_of_patch is unused in the body (Fortran kept `c` but never reads it).
@kernel function _psn_fractionation_kernel!(alphapsn, @Const(mask_patch), @Const(nrad),
        @Const(par_z_in), @Const(an_ref), @Const(gs_mol_ref), @Const(gb_mol),
        @Const(forc_pbot), @Const(forc_pco2), @Const(c3psn_pft), @Const(ivt),
        @Const(gridcell_of_patch))
    p = @index(Global)
    @inbounds if mask_patch[p]
        T = eltype(par_z_in)
        g = gridcell_of_patch[p]
        co2_p = forc_pco2[g]
        for iv in 1:nrad[p]
            if par_z_in[p, iv] <= zero(T)
                alphapsn[p] = one(T)
            else
                ci = co2_p - (an_ref[p, iv] * forc_pbot[p] *
                    (T(1.4) * gs_mol_ref[p, iv] + T(1.6) * gb_mol[p]) /
                    (gb_mol[p] * gs_mol_ref[p, iv]))
                alphapsn[p] = one(T) + (((c3psn_pft[ivt[p]] *
                    (T(4.4) + (T(22.6) * (ci / co2_p)))) +
                    ((one(T) - c3psn_pft[ivt[p]]) * T(4.4))) / T(1000.0))
            end
        end
    end
end

"""
    fractionation!(ps, ...)

C13 fractionation during photosynthesis.
Ported from `subroutine Fractionation`.
"""
function fractionation!(ps,
                        forc_pbot::AbstractVector{<:Real},
                        forc_pco2::AbstractVector{<:Real},
                        par_z_in::AbstractMatrix{<:Real},
                        nrad::AbstractVector{<:Integer},
                        c3psn_pft::AbstractVector{<:Real},
                        ivt::AbstractVector{<:Integer},
                        col_of_patch::AbstractVector{<:Integer},
                        gridcell_of_patch::AbstractVector{<:Integer},
                        mask_patch::AbstractVector{Bool},
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

    _launch!(_psn_fractionation_kernel!, alphapsn, mask_patch, nrad, par_z_in,
        an_ref, gs_mol_ref, gb_mol, forc_pbot, forc_pco2, c3psn_pft, ivt,
        gridcell_of_patch; ndrange = length(bounds_patch))

    return nothing
end

# =====================================================================
# TimeStepInit — Initialize at start of time step
# =====================================================================

# Per-patch zero of photosynthesis fluxes (mask-based standalone form; the
# canopy_fluxes! call site uses the filterp-gather variant _cf_psn_init_kernel!).
# use_c13/use_c14 are scalar Bool kernel args so the optional fields are written
# only when active (byte-identical to the host flag branches).
@kernel function _psn_timestep_init_kernel!(psnsun, psnsun_wc, psnsun_wj, psnsun_wp,
        psnsha, psnsha_wc, psnsha_wj, psnsha_wp, fpsn, fpsn_wc, fpsn_wj, fpsn_wp,
        alphapsnsun, alphapsnsha, c13_psnsun, c13_psnsha, c14_psnsun, c14_psnsha,
        @Const(mask_nolake), use_c13::Bool, use_c14::Bool)
    p = @index(Global)
    @inbounds if mask_nolake[p]
        z = zero(eltype(psnsun))
        psnsun[p] = z; psnsun_wc[p] = z; psnsun_wj[p] = z; psnsun_wp[p] = z
        psnsha[p] = z; psnsha_wc[p] = z; psnsha_wj[p] = z; psnsha_wp[p] = z
        fpsn[p]   = z; fpsn_wc[p]   = z; fpsn_wj[p]   = z; fpsn_wp[p]   = z
        if use_c13
            alphapsnsun[p] = z; alphapsnsha[p] = z
            c13_psnsun[p]  = z; c13_psnsha[p]  = z
        end
        if use_c14
            c14_psnsun[p] = z; c14_psnsha[p] = z
        end
    end
end

"""
    photosynthesis_timestep_init!(ps, mask_nolake, bounds_patch;
                                  use_c13=false, use_c14=false)

Time step initialization. Zeros out photosynthesis fluxes for non-lake patches.
"""
function photosynthesis_timestep_init!(ps,
                                       mask_nolake::AbstractVector{Bool},
                                       bounds_patch::UnitRange{Int};
                                       use_c13::Bool=false,
                                       use_c14::Bool=false)
    _launch!(_psn_timestep_init_kernel!, ps.psnsun_patch, ps.psnsun_wc_patch,
        ps.psnsun_wj_patch, ps.psnsun_wp_patch, ps.psnsha_patch, ps.psnsha_wc_patch,
        ps.psnsha_wj_patch, ps.psnsha_wp_patch, ps.fpsn_patch, ps.fpsn_wc_patch,
        ps.fpsn_wj_patch, ps.fpsn_wp_patch, ps.alphapsnsun_patch, ps.alphapsnsha_patch,
        ps.c13_psnsun_patch, ps.c13_psnsha_patch, ps.c14_psnsun_patch, ps.c14_psnsha_patch,
        mask_nolake, use_c13, use_c14; ndrange = length(bounds_patch))

    return nothing
end

# =====================================================================
# SetParamsForTesting — Set parameters for unit testing
# =====================================================================

"""
    set_params_for_testing!(ps::PhotosynthesisData)

Set parameters for unit testing (from `subroutine setParamsForTesting`).
"""
function set_params_for_testing!(ps::PhotosynthesisData{FT}) where {FT}
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
function getqflx!(p::Int, c::Int, gb_mol::Real,
                  gs_mol_sun::Real, gs_mol_sha::Real,
                  qflx_sun::Real, qflx_sha::Real,
                  qsatl::Real, qaf::Real, havegs::Bool,
                  laisun_p::Real, laisha_p::Real,
                  elai_p::Real, esai_p::Real,
                  fdry_p::Real, forc_rho_c::Real,
                  forc_pbot_c::Real, tgcm_p::Real)

    cf = forc_pbot_c / (RGAS * tgcm_p) * 1.0e6
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
# GPU-callable POSITIONAL device cores for the PHS Newton chain.
# Each mirrors a host helper above bit-for-bit on Float64, but:
#   * NO vectors / StaticArrays — the length-4 vegwp unknown is carried as 4
#     NAMED SCALARS (xsun=x[SUN], xsha=x[SHA], xxyl=x[XYL], xroot=x[ROOT_SEG]);
#   * NO params_inst / globals — psi50/ck/kmax/theta come from a threaded
#     PsnPhsParams bundle; vulnerability curve via _plc/_d1plc;
#   * NO error(); NO kwargs (positional); every Float64 literal T()-wrapped;
#   * length-nlevsoi reductions become explicit serial j-loops over the 2D
#     k_soil_root[p,·]/smp_l[c,·]/z_col[c,·] (passed whole + (p,c) indices, so
#     no device slices). The host's TWO separate sums in spacF / getvegwp are
#     kept as TWO separate accumulators (combining them would re-order rounding).
# These cores are what the per-patch Pass-3 kernel calls in-thread. The original
# host helpers above are retained for the standalone unit tests.
# =====================================================================

# _getqflx — pure scalar 4-tuple form of getqflx!. forc_pbot_c/tgcm_p/etc are
# scalars; T from gb_mol. Byte-identical to getqflx! on Float64.
@inline function _getqflx(gb_mol::Real,
                          gs_mol_sun::Real, gs_mol_sha::Real,
                          qflx_sun::Real, qflx_sha::Real,
                          qsatl::Real, qaf::Real, havegs::Bool,
                          laisun_p::Real, laisha_p::Real,
                          elai_p::Real, esai_p::Real,
                          fdry_p::Real, forc_rho_c::Real,
                          forc_pbot_c::Real, tgcm_p::Real, RGAS_::Real)
    T = typeof(gb_mol)
    cf = forc_pbot_c / (RGAS_ * tgcm_p) * T(1.0e6)
    wtl = (elai_p + esai_p) * gb_mol
    efpot = forc_rho_c * wtl * (qsatl - qaf)

    qsun = qflx_sun
    qsha = qflx_sha
    gsun = gs_mol_sun
    gsha = gs_mol_sha
    if havegs
        if efpot > zero(T) && elai_p > zero(T)
            if gs_mol_sun > zero(T)
                rppdry_sun = fdry_p / gb_mol * (laisun_p / (one(T) / gb_mol + one(T) / gs_mol_sun)) / elai_p
                qsun = efpot * rppdry_sun / cf
            else
                qsun = zero(T)
            end
            if gs_mol_sha > zero(T)
                rppdry_sha = fdry_p / gb_mol * (laisha_p / (one(T) / gb_mol + one(T) / gs_mol_sha)) / elai_p
                qsha = efpot * rppdry_sha / cf
            else
                qsha = zero(T)
            end
        else
            qsun = zero(T)
            qsha = zero(T)
        end
    else
        if qflx_sun > zero(T)
            gsun = gb_mol * qflx_sun * cf * elai_p / (efpot * fdry_p * laisun_p - qflx_sun * cf * elai_p)
        else
            gsun = zero(T)
        end
        if qflx_sha > zero(T)
            gsha = gb_mol * qflx_sha * cf * elai_p / (efpot * fdry_p * laisha_p - qflx_sha * cf * elai_p)
        else
            gsha = zero(T)
        end
    end
    return (qsun, qsha, gsun, gsha)
end

# _spacF — scalar-unrolled spacF!. Returns the 4 residual components as 4 scalars
# (fsun, fsha, fxyl, froot mapping f[SUN],f[SHA],f[XYL],f[ROOT_SEG]). The vegwp
# unknown is passed as 4 named scalars (xsun..xroot). kmax/psi50/ck via phs_params;
# the two nlevsoi sums kept as two serial accumulators (host order preserved).
@inline function _spacF(xsun::Real, xsha::Real, xxyl::Real, xroot::Real,
                        qflx_sun::Real, qflx_sha::Real,
                        k_soil_root, smp_l, z_col, p::Int, c::Int,
                        laisun_p::Real, laisha_p::Real,
                        htop_p::Real, tsai_p::Real, ivt_p::Int, nlevsoi::Int,
                        phs_params)
    T = typeof(xsun)
    tol_lai = T(0.001)

    psi50 = phs_params.psi50; ck = phs_params.ck; kmax = phs_params.kmax
    kmax_sun = kmax[ivt_p, SUN]; kmax_sha = kmax[ivt_p, SHA]; kmax_xyl = kmax[ivt_p, XYL]

    grav1 = htop_p * T(1000.0)

    fsto1 = _plc(xsun,  psi50[ivt_p, SUN],      ck[ivt_p, SUN])
    fsto2 = _plc(xsha,  psi50[ivt_p, SHA],      ck[ivt_p, SHA])
    fx    = _plc(xxyl,  psi50[ivt_p, XYL],      ck[ivt_p, XYL])
    fr    = _plc(xroot, psi50[ivt_p, ROOT_SEG], ck[ivt_p, ROOT_SEG])

    fsun = qflx_sun * fsto1 - laisun_p * kmax_sun * fx * (xxyl - xsun)
    fsha = qflx_sha * fsto2 - laisha_p * kmax_sha * fx * (xxyl - xsha)
    fxyl = laisun_p * kmax_sun * fx * (xxyl - xsun) +
           laisha_p * kmax_sha * fx * (xxyl - xsha) -
           tsai_p * kmax_xyl / htop_p * fr * (xroot - xxyl - grav1)
    # f[ROOT_SEG] = xyl-root flux + sum(ksr*(xroot+grav2)) - sum(ksr*smp). Keep the
    # two sums separate (host order); grav2[j] = z_col[c,j]*1000.
    s1 = zero(T)
    @inbounds for j in 1:nlevsoi
        s1 += k_soil_root[p, j] * (xroot + z_col[c, j] * T(1000.0))
    end
    s2 = zero(T)
    @inbounds for j in 1:nlevsoi
        s2 += k_soil_root[p, j] * smp_l[c, j]
    end
    froot = tsai_p * kmax_xyl / htop_p * fr * (xroot - xxyl - grav1) + s1 - s2

    if laisha_p < tol_lai
        tmp = fsun
        fsun = fsha
        fsha = tmp
    end
    return (fsun, fsha, fxyl, froot)
end

# _spacA — scalar-unrolled spacA!. Returns the (up to 16) invA entries as named
# scalars + a singular flag (the two host early-returns become flag=true with the
# already-zeroed invA). Unused A/invA entries are exactly zero as in the host
# (A/invA start zeroed). Closed-form inverse copied verbatim.
@inline function _spacA(xsun::Real, xsha::Real, xxyl::Real, xroot::Real,
                        qflx_sun::Real, qflx_sha::Real,
                        k_soil_root, p::Int,
                        laisun_p::Real, laisha_p::Real,
                        htop_p::Real, tsai_p::Real, ivt_p::Int, nlevsoi::Int,
                        phs_params)
    T = typeof(xsun)
    tol_lai = T(0.001)
    z = zero(T)
    # invA entries (row-major names iaRC); default zero like the host invA.
    ia11 = z; ia12 = z; ia13 = z; ia14 = z
    ia21 = z; ia22 = z; ia23 = z; ia24 = z
    ia31 = z; ia32 = z; ia33 = z; ia34 = z
    ia41 = z; ia42 = z; ia43 = z; ia44 = z

    psi50 = phs_params.psi50; ck = phs_params.ck; kmax = phs_params.kmax
    kmax_sun = kmax[ivt_p, SUN]; kmax_sha = kmax[ivt_p, SHA]; kmax_xyl = kmax[ivt_p, XYL]

    grav1 = htop_p * T(1000.0)

    fsto1 = _plc(xsun,  psi50[ivt_p, SUN],      ck[ivt_p, SUN])
    fsto2 = _plc(xsha,  psi50[ivt_p, SHA],      ck[ivt_p, SHA])
    fx    = _plc(xxyl,  psi50[ivt_p, XYL],      ck[ivt_p, XYL])
    fr    = _plc(xroot, psi50[ivt_p, ROOT_SEG], ck[ivt_p, ROOT_SEG])

    dfsto1 = _d1plc(xsun,  psi50[ivt_p, SUN],      ck[ivt_p, SUN])
    dfsto2 = _d1plc(xsha,  psi50[ivt_p, SHA],      ck[ivt_p, SHA])
    dfx    = _d1plc(xxyl,  psi50[ivt_p, XYL],      ck[ivt_p, XYL])
    dfr    = _d1plc(xroot, psi50[ivt_p, ROOT_SEG], ck[ivt_p, ROOT_SEG])

    # A matrix entries (named aRC); zero where the host leaves them zero.
    a11 = -laisun_p * kmax_sun * fx - qflx_sun * dfsto1
    a13 = laisun_p * kmax_sun * dfx * (xxyl - xsun) + laisun_p * kmax_sun * fx
    a22 = -laisha_p * kmax_sha * fx - qflx_sha * dfsto2
    a23 = laisha_p * kmax_sha * dfx * (xxyl - xsha) + laisha_p * kmax_sha * fx
    a31 = laisun_p * kmax_sun * fx
    a32 = laisha_p * kmax_sha * fx
    a33 = -laisun_p * kmax_sun * dfx * (xxyl - xsun) -
           laisun_p * kmax_sun * fx -
           laisha_p * kmax_sha * dfx * (xxyl - xsha) -
           laisha_p * kmax_sha * fx -
           tsai_p * kmax_xyl / htop_p * fr
    a34 = tsai_p * kmax_xyl / htop_p * dfr * (xroot - xxyl - grav1) +
          tsai_p * kmax_xyl / htop_p * fr
    a43 = tsai_p * kmax_xyl / htop_p * fr
    # sum(k_soil_root) serial accumulator
    sksr = zero(T)
    @inbounds for j in 1:nlevsoi
        sksr += k_soil_root[p, j]
    end
    a44 = -tsai_p * kmax_xyl / htop_p * fr -
           tsai_p * kmax_xyl / htop_p * dfr * (xroot - xxyl - grav1) -
           sksr

    if laisun_p > tol_lai && laisha_p > tol_lai
        # General 4x4 case
        determ = a44*a22*a33*a11 - a44*a22*a31*a13 -
                 a44*a32*a23*a11 - a43*a11*a22*a34
        if abs(determ) <= T(1.0e-50)
            return (ia11, ia12, ia13, ia14, ia21, ia22, ia23, ia24,
                    ia31, ia32, ia33, ia34, ia41, ia42, ia43, ia44, true)
        end
        leading = one(T) / determ
        ia11 = leading*a44*a22*a33 - leading*a44*a32*a23 - leading*a43*a22*a34
        ia21 = leading*a23*a44*a31
        ia31 = -leading*a44*a22*a31
        ia41 = leading*a43*a22*a31
        ia12 = leading*a13*a44*a32
        ia22 = leading*a44*a33*a11 - leading*a44*a31*a13 - leading*a43*a11*a34
        ia32 = -leading*a11*a44*a32
        ia42 = leading*a43*a11*a32
        ia13 = -leading*a13*a22*a44
        ia23 = -leading*a23*a11*a44
        ia33 = leading*a22*a11*a44
        ia43 = -leading*a43*a11*a22
        ia14 = leading*a13*a34*a22
        ia24 = leading*a23*a34*a11
        ia34 = -leading*a34*a11*a22
        ia44 = leading*a22*a33*a11 - leading*a22*a31*a13 - leading*a32*a23*a11
    else
        # 3x3 case when one of laisun/laisha is ~0. The host mutates A entries
        # when laisha<=tol; mirror with locals so the formulas read identically.
        a22l = a22; a32l = a32; a23l = a23
        if laisha_p <= tol_lai
            a22l = a11
            a32l = a31
            a23l = a13
        end
        determ = a22l*a33*a44 - a34*a22l*a43 - a23l*a32l*a44
        if abs(determ) <= T(1.0e-50)
            return (ia11, ia12, ia13, ia14, ia21, ia22, ia23, ia24,
                    ia31, ia32, ia33, ia34, ia41, ia42, ia43, ia44, true)
        end
        # build then scale by 1/determ (host: invA .= (1/determ).*invA)
        b22 = a33*a44 - a34*a43
        b23 = -a23l*a44
        b24 = a34*a23l
        b32 = -a32l*a44
        b33 = a22l*a44
        b34 = -a34*a22l
        b42 = a32l*a43
        b43 = -a22l*a43
        b44 = a22l*a33 - a23l*a32l
        invd = one(T) / determ
        ia22 = invd*b22; ia23 = invd*b23; ia24 = invd*b24
        ia32 = invd*b32; ia33 = invd*b33; ia34 = invd*b34
        ia42 = invd*b42; ia43 = invd*b43; ia44 = invd*b44
    end
    return (ia11, ia12, ia13, ia14, ia21, ia22, ia23, ia24,
            ia31, ia32, ia33, ia34, ia41, ia42, ia43, ia44, false)
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
function spacF!(p::Int, c::Int, x::AbstractVector{<:Real},
                qflx_sun::Real, qflx_sha::Real,
                k_soil_root_p::AbstractVector{Float64},
                smp_c::AbstractVector{Float64},
                z_c::AbstractVector{Float64},
                laisun_p::Real, laisha_p::Real,
                htop_p::Real, tsai_p::Real,
                ivt_p::Int, nlevsoi::Int)
    tol_lai = 0.001

    grav1 = htop_p * 1000.0
    grav2 = z_c[1:nlevsoi] .* 1000.0

    fsto1 = plc(x[SUN], ivt_p, SUN, VEG)
    fsto2 = plc(x[SHA], ivt_p, SHA, VEG)
    fx    = plc(x[XYL], ivt_p, XYL, VEG)
    fr    = plc(x[ROOT_SEG], ivt_p, ROOT_SEG, VEG)

    FT = eltype(x)
    f = zeros(FT, NVEGWCS)
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
function spacA!(p::Int, c::Int, x::AbstractVector{<:Real},
                qflx_sun::Real, qflx_sha::Real,
                k_soil_root_p::AbstractVector{Float64},
                laisun_p::Real, laisha_p::Real,
                htop_p::Real, tsai_p::Real,
                ivt_p::Int, nlevsoi::Int)
    tol_lai = 0.001

    FT = eltype(x)
    A = zeros(FT, NVEGWCS, NVEGWCS)
    invA = zeros(FT, NVEGWCS, NVEGWCS)

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
function getvegwp!(p::Int, c::Int, gb_mol::Real,
                   gs_mol_sun::Real, gs_mol_sha::Real,
                   qsatl::Real, qaf::Real,
                   k_soil_root_p::AbstractVector{Float64},
                   smp_c::AbstractVector{Float64},
                   z_c::AbstractVector{Float64},
                   laisun_p::Real, laisha_p::Real,
                   htop_p::Real, tsai_p::Real,
                   ivt_p::Int, nlevsoi::Int,
                   elai_p::Real, esai_p::Real,
                   fdry_p::Real, forc_rho_c::Real,
                   forc_pbot_c::Real, tgcm_p::Real)
    FT = typeof(gb_mol)
    x = zeros(FT, NVEGWCS)

    grav1 = 1000.0 * htop_p
    grav2 = 1000.0 .* z_c[1:nlevsoi]

    # Compute transpiration demand
    havegs = true
    qflx_sun, qflx_sha, _, _ = getqflx!(p, c, gb_mol, gs_mol_sun, gs_mol_sha,
        0.0, 0.0, qsatl, qaf, havegs, laisun_p, laisha_p,
        elai_p, esai_p, fdry_p, forc_rho_c, forc_pbot_c, tgcm_p)

    # Calculate root water potential. Threshold (not ==0): when ALL root layers are
    # frozen the smooth_max floor in pass-1 leaves k_soil_root ~1e-16 (not exactly 0,
    # an AD-safety artifact), so sum_ksr ~1e-17 would bypass a literal ==0 check and
    # the -qflx/sum_ksr term below would give a huge xroot → vegwp Newton NaN / bogus
    # vegwp. Below the threshold = effectively no root water access: use the k-weighted
    # mean of (smp-grav2) for xroot (the qflx→0 limit), which is bounded AND zero-
    # weights the bedrock layers (k=0, smp pinned at smpmin=-1e8) — an UNWEIGHTED mean
    # would be dragged to ~-4e7 by that dry bedrock. soilflux is the mass-balance
    # demand below (= Fortran else-branch's soilflux=qflx).
    sum_ksr = sum(k_soil_root_p[1:nlevsoi])
    no_root_access = abs(sum_ksr) <= oftype(sum_ksr, 1.0e-12)
    ksmp = sum(k_soil_root_p[1:nlevsoi] .* (smp_c[1:nlevsoi] .- grav2))
    if no_root_access
        x[ROOT_SEG] = abs(sum_ksr) > zero(sum_ksr) ? ksmp / sum_ksr : zero(eltype(x))
    else
        x[ROOT_SEG] = (ksmp - qflx_sun - qflx_sha) / sum_ksr
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

    # Calculate soil flux. With no root water access (frozen soil; sum_ksr below
    # threshold), soilflux = transpiration demand by mass balance (= what the
    # else-branch loop yields by cancellation). The threshold branch's unweighted-mean
    # xroot breaks that cancellation, leaving a spurious k·smp residual; use the
    # demand instead. See _getvegwp for the detailed rationale.
    soilflux = 0.0
    if no_root_access
        soilflux = qflx_sun + qflx_sha
    else
        for j in 1:nlevsoi
            soilflux += k_soil_root_p[j] * (smp_c[j] - x[ROOT_SEG] - grav2[j])
        end
    end

    return (x, soilflux)
end

# _getvegwp — scalar-unrolled getvegwp!. Returns the 4 vegwp scalars + soilflux:
# (xsun, xsha, xxyl, xroot, soilflux). Inlines _getqflx + _plc. The host's
# nlevsoi reductions become serial accumulators in HOST ORDER (sum_ksr, the
# combined sum(smp.-grav2), the combined sum(ksr.*(smp.-grav2)), and soilflux).
@inline function _getvegwp(gb_mol::Real, gs_mol_sun::Real, gs_mol_sha::Real,
                           qsatl::Real, qaf::Real,
                           k_soil_root, smp_l, z_col, p::Int, c::Int,
                           laisun_p::Real, laisha_p::Real,
                           htop_p::Real, tsai_p::Real, ivt_p::Int, nlevsoi::Int,
                           elai_p::Real, esai_p::Real,
                           fdry_p::Real, forc_rho_c::Real,
                           forc_pbot_c::Real, tgcm_p::Real,
                           phs_params, RGAS_::Real)
    T = typeof(gb_mol)
    psi50 = phs_params.psi50; ck = phs_params.ck; kmax = phs_params.kmax
    grav1 = T(1000.0) * htop_p

    # Compute transpiration demand
    qflx_sun, qflx_sha, _, _ = _getqflx(gb_mol, gs_mol_sun, gs_mol_sha,
        zero(T), zero(T), qsatl, qaf, true, laisun_p, laisha_p,
        elai_p, esai_p, fdry_p, forc_rho_c, forc_pbot_c, tgcm_p, RGAS_)

    # Root water potential. sum_ksr in host order.
    sum_ksr = zero(T)
    @inbounds for j in 1:nlevsoi
        sum_ksr += k_soil_root[p, j]
    end
    xroot = zero(T)
    # Threshold (not ==0): all-frozen soil leaves k_soil_root at the pass-1 smooth_max
    # floor (~1e-16, an AD artifact) so sum_ksr ~1e-17 ≠ 0; the full else-branch's
    # -qflx/sum_ksr term would then blow xroot → ±1e11 (Newton NaN / bogus vegwp).
    # Below threshold = effectively no root water access: use the k-weighted mean of
    # (smp-grav2) for xroot (the qflx→0 limit of the else-branch). This is bounded
    # (a weighted average, so within the smp range) AND zero-weights the bedrock
    # layers (k_soil_root=0 there, smp pinned at smpmin=-1e8); an UNWEIGHTED mean
    # would be dragged to ~-4e7 by that dry bedrock. soilflux is set to the
    # mass-balance demand below (= Fortran else-branch's soilflux=qflx).
    no_root_access = abs(sum_ksr) <= T(1.0e-12)
    s = zero(T)
    @inbounds for j in 1:nlevsoi
        s += k_soil_root[p, j] * (smp_l[c, j] - z_col[c, j] * T(1000.0))
    end
    if no_root_access
        xroot = abs(sum_ksr) > zero(T) ? s / sum_ksr : zero(T)
    else
        xroot = (s - qflx_sun - qflx_sha) / sum_ksr
    end

    # Xylem water potential
    fr = _plc(xroot, psi50[ivt_p, ROOT_SEG], ck[ivt_p, ROOT_SEG])
    if tsai_p > zero(T) && fr > zero(T)
        xxyl = xroot - grav1 - (qflx_sun + qflx_sha) / (fr * kmax[ivt_p, ROOT_SEG] / htop_p * tsai_p)
    else
        xxyl = xroot - grav1
    end

    # Sun/sha leaf water potential
    fx = _plc(xxyl, psi50[ivt_p, XYL], ck[ivt_p, XYL])
    if laisha_p > zero(T) && fx > zero(T)
        xsha = xxyl - qflx_sha / (fx * kmax[ivt_p, XYL] * laisha_p)
    else
        xsha = xxyl
    end
    if laisun_p > zero(T) && fx > zero(T)
        xsun = xxyl - qflx_sun / (fx * kmax[ivt_p, XYL] * laisun_p)
    else
        xsun = xxyl
    end

    # Soil flux (host order: grav2[j]=z_col[c,j]*1000). When there is effectively no
    # root water access (sum_ksr below threshold, e.g. frozen soil), the soil-to-root
    # flux equals the transpiration demand by mass balance — which is exactly what the
    # else-branch loop yields via cancellation (soilflux = Σk(smp-xroot-grav2) =
    # qflx_sun+qflx_sha). The threshold branch uses the UNWEIGHTED mean for xroot
    # (NaN-safety), which breaks that cancellation and would otherwise leave a spurious
    # residual ∝ k·smp (huge with very-negative frozen smp → bogus late-season
    # transpiration). Set soilflux to the mass-balance demand instead. Fortran's
    # getvegwp takes the else-branch (sum(k_soil_root)==0 only when ALL layers are 0),
    # so its soilflux is always = qflx; this matches that limit.
    soilflux = zero(T)
    if no_root_access
        soilflux = qflx_sun + qflx_sha
    else
        @inbounds for j in 1:nlevsoi
            soilflux += k_soil_root[p, j] * (smp_l[c, j] - xroot - z_col[c, j] * T(1000.0))
        end
    end

    return (xsun, xsha, xxyl, xroot, soilflux)
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
function calcstress!(p::Int, c::Int, x::AbstractVector{<:Real},
                     gb_mol::Real, gs_mol_sun_in::Real, gs_mol_sha_in::Real,
                     qsatl::Real, qaf::Real,
                     k_soil_root_p::AbstractVector{Float64},
                     smp_c::AbstractVector{Float64},
                     z_c::AbstractVector{Float64},
                     laisun_p::Real, laisha_p::Real,
                     htop_p::Real, tsai_p::Real,
                     ivt_p::Int, nlevsoi::Int,
                     elai_p::Real, esai_p::Real,
                     fdry_p::Real, forc_rho_c::Real,
                     forc_pbot_c::Real, tgcm_p::Real)
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

    FT = eltype(x)
    dx = zeros(FT, NVEGWCS)  # Pre-allocate Newton step vector

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
                # dx = invA * f (no-alloc manual multiply)
                for i in 1:NVEGWCS
                    dx[i] = 0.0
                    for j in 1:NVEGWCS
                        dx[i] += invA[i, j] * f[j]
                    end
                end
            else
                dx[SUN] = 0.0
                for i in SHA:ROOT_SEG
                    dx[i] = 0.0
                    for j in SHA:ROOT_SEG
                        dx[i] += invA[i, j] * f[j]
                    end
                end
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

# _calcstress — scalar-unrolled calcstress! Newton orchestrator. Takes/returns the
# vegwp unknown as 4 named scalars; returns (xsun, xsha, xxyl, xroot, bsun, bsha).
# The host `while true` Newton loop becomes a fixed-trip `for iter in 1:(itmax+1)`
# with a `done::Bool` that, once set, FREEZES the body (every later iteration is a
# no-op) → uniform trip count, NO divergent break (Metal-safe). `flag` is captured
# at the moment the loop finishes (it drives the algebraic-fallback branch below).
# The 4x4 / 3x3 solves are scalar-unrolled in the exact host index order.
@inline function _calcstress(xsun0::Real, xsha0::Real, xxyl0::Real, xroot0::Real,
                             gb_mol::Real, gs_mol_sun_in::Real, gs_mol_sha_in::Real,
                             qsatl::Real, qaf::Real,
                             k_soil_root, smp_l, z_col, p::Int, c::Int,
                             laisun_p::Real, laisha_p::Real,
                             htop_p::Real, tsai_p::Real, ivt_p::Int, nlevsoi::Int,
                             elai_p::Real, esai_p::Real,
                             fdry_p::Real, forc_rho_c::Real,
                             forc_pbot_c::Real, tgcm_p::Real,
                             phs_params, RGAS_::Real)
    T = typeof(xsun0)
    itmax = 50
    tolf = T(1.0e-6)
    toldx = T(1.0e-9)
    tol_lai = T(0.001)
    psi50 = phs_params.psi50; ck = phs_params.ck

    xsun = xsun0; xsha = xsha0; xxyl = xxyl0; xroot = xroot0

    # Night flag: vegwp(sun) > 0 signals night
    night = false
    if xsun > zero(T)
        night = true
        xsun = xsha
    end

    gs0sun = gs_mol_sun_in
    gs0sha = gs_mol_sha_in

    # Transpiration demand
    qflx_sun, qflx_sha, _, _ = _getqflx(gb_mol, gs0sun, gs0sha,
        zero(T), zero(T), qsatl, qaf, true, laisun_p, laisha_p,
        elai_p, esai_p, fdry_p, forc_rho_c, forc_pbot_c, tgcm_p, RGAS_)

    flag = false

    if (laisun_p > tol_lai || laisha_p > tol_lai) && (qflx_sun > zero(T) || qflx_sha > zero(T))
        # Newton's method (fixed-trip, frozen-once-done). The host increments iter
        # at the top then breaks on iter>itmax — so the body runs for iter=1..itmax+1
        # (the +1 pass only evaluates spacF + its convergence test before breaking).
        done = false
        @inbounds for iter in 1:(itmax + 1)
            if !done
                fsun, fsha, fxyl, froot = _spacF(xsun, xsha, xxyl, xroot,
                    qflx_sun, qflx_sha, k_soil_root, smp_l, z_col, p, c,
                    laisun_p, laisha_p, htop_p, tsai_p, ivt_p, nlevsoi, phs_params)

                normf = sqrt(fsun*fsun + fsha*fsha + fxyl*fxyl + froot*froot)
                if normf < tolf * (qflx_sun + qflx_sha)
                    flag = false
                    done = true
                elseif iter > itmax
                    flag = false
                    done = true
                else
                    ia11, ia12, ia13, ia14, ia21, ia22, ia23, ia24,
                    ia31, ia32, ia33, ia34, ia41, ia42, ia43, ia44, sing =
                        _spacA(xsun, xsha, xxyl, xroot, qflx_sun, qflx_sha,
                            k_soil_root, p, laisun_p, laisha_p, htop_p, tsai_p,
                            ivt_p, nlevsoi, phs_params)
                    flag = sing
                    if sing
                        done = true
                    else
                        # dx = invA * f. General 4x4 vs 3x3 (laisun==0 → dx[SUN]=0).
                        if laisun_p > tol_lai && laisha_p > tol_lai
                            dxsun  = ia11*fsun + ia12*fsha + ia13*fxyl + ia14*froot
                            dxsha  = ia21*fsun + ia22*fsha + ia23*fxyl + ia24*froot
                            dxxyl  = ia31*fsun + ia32*fsha + ia33*fxyl + ia34*froot
                            dxroot = ia41*fsun + ia42*fsha + ia43*fxyl + ia44*froot
                        else
                            dxsun  = zero(T)
                            dxsha  = ia22*fsha + ia23*fxyl + ia24*froot
                            dxxyl  = ia32*fsha + ia33*fxyl + ia34*froot
                            dxroot = ia42*fsha + ia43*fxyl + ia44*froot
                        end

                        # maximum(abs.(dx)) over all 4 (dx[SUN]=0 in 3x3 case).
                        mx = abs(dxsun)
                        adsha = abs(dxsha); if adsha > mx; mx = adsha; end
                        adxyl = abs(dxxyl); if adxyl > mx; mx = adxyl; end
                        adroot = abs(dxroot); if adroot > mx; mx = adroot; end
                        if mx > T(50000.0)
                            dxsun  = T(50000.0) * dxsun  / mx
                            dxsha  = T(50000.0) * dxsha  / mx
                            dxxyl  = T(50000.0) * dxxyl  / mx
                            dxroot = T(50000.0) * dxroot / mx
                        end

                        # Update x per lai case
                        if laisun_p > tol_lai && laisha_p > tol_lai
                            xsun  += dxsun
                            xsha  += dxsha
                            xxyl  += dxxyl
                            xroot += dxroot
                        elseif laisha_p > tol_lai
                            xsun  += dxsun
                            xsha  += dxsha
                            xxyl  += dxxyl
                            xroot += dxroot
                            xsun = xxyl  # psi_sun = psi_xyl (laisun==0)
                        else
                            xxyl  += dxxyl
                            xroot += dxroot
                            xsun  += dxsha   # x[SUN] += dx[SHA]
                            xsha = xxyl      # psi_sha = psi_xyl (laisha==0)
                        end

                        normdx = sqrt(dxsun*dxsun + dxsha*dxsha + dxxyl*dxxyl + dxroot*dxroot)
                        if normdx < toldx
                            done = true
                        else
                            # Force SPAC gradient toward atmosphere
                            if xxyl > xroot
                                xxyl = xroot
                            end
                            if xsun > xxyl
                                xsun = xxyl
                            end
                            if xsha > xxyl
                                xsha = xxyl
                            end
                        end
                    end
                end
            end
        end
    else
        flag = true
    end

    bsun = zero(T)
    bsha = zero(T)

    if flag
        # Solve algebraically
        xsun, xsha, xxyl, xroot, _ = _getvegwp(gb_mol, gs0sun, gs0sha, qsatl, qaf,
            k_soil_root, smp_l, z_col, p, c, laisun_p, laisha_p, htop_p, tsai_p,
            ivt_p, nlevsoi, elai_p, esai_p, fdry_p, forc_rho_c, forc_pbot_c, tgcm_p,
            phs_params, RGAS_)
        bsun = _plc(xsun, psi50[ivt_p, SUN], ck[ivt_p, SUN])
        bsha = _plc(xsha, psi50[ivt_p, SHA], ck[ivt_p, SHA])
    else
        qsun = qflx_sun * _plc(xsun, psi50[ivt_p, SUN], ck[ivt_p, SUN])
        qsha = qflx_sha * _plc(xsha, psi50[ivt_p, SHA], ck[ivt_p, SHA])

        _, _, gs0sun_stressed, gs0sha_stressed = _getqflx(gb_mol, gs0sun, gs0sha,
            qsun, qsha, qsatl, qaf, false, laisun_p, laisha_p,
            elai_p, esai_p, fdry_p, forc_rho_c, forc_pbot_c, tgcm_p, RGAS_)

        if qflx_sun > zero(T)
            bsun = gs0sun_stressed / gs_mol_sun_in
        else
            bsun = _plc(xsun, psi50[ivt_p, SUN], ck[ivt_p, SUN])
        end
        if qflx_sha > zero(T)
            bsha = gs0sha_stressed / gs_mol_sha_in
        else
            bsha = _plc(xsha, psi50[ivt_p, SHA], ck[ivt_p, SHA])
        end
    end

    if bsun < T(0.01)
        bsun = zero(T)
    end
    if bsha < T(0.01)
        bsha = zero(T)
    end

    soilflux = zero(T)
    if night
        gs0sun_night = bsun * gs_mol_sun_in
        gs0sha_night = bsha * gs_mol_sha_in
        xsun, xsha, xxyl, xroot, soilflux = _getvegwp(gb_mol, gs0sun_night, gs0sha_night,
            qsatl, qaf, k_soil_root, smp_l, z_col, p, c, laisun_p, laisha_p,
            htop_p, tsai_p, ivt_p, nlevsoi, elai_p, esai_p, fdry_p, forc_rho_c,
            forc_pbot_c, tgcm_p, phs_params, RGAS_)
    end

    # Return the night soil-water flux so the caller can set qflx_tran_veg at night.
    # Fortran sets qflx_tran_veg(p)=soilflux inside getvegwp (PhotosynthesisMod ~4039),
    # which the night calcstress->getvegwp path also hits; the Julia _getvegwp returns
    # soilflux for the caller, and the night psn pass would otherwise leave
    # qflx_tran_veg at its NaN init (transpiration NaN -> t_veg/t_grnd/h2osoi NaN).
    return (xsun, xsha, xxyl, xroot, bsun, bsha, soilflux)
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
function ci_func_PHS!(x::AbstractVector{<:Real}, cisun::Real, cisha::Real,
                      p::Int, iv::Int, c::Int,
                      bsun::Real, bsha::Real, bflag::Bool,
                      gb_mol::Real, gs0sun::Real, gs0sha::Real,
                      gs_mol_sun::Real, gs_mol_sha::Real,
                      jesun::Real, jesha::Real,
                      cair::Real, oair::Real,
                      lmr_z_sun::Real, lmr_z_sha::Real,
                      par_z_sun::Real, par_z_sha::Real,
                      rh_can::Real, qsatl::Real, qaf::Real,
                      ps, forc_pbot_c::Real,
                      k_soil_root_p::AbstractVector{Float64},
                      smp_c::AbstractVector{Float64},
                      z_c::AbstractVector{Float64},
                      laisun_p::Real, laisha_p::Real,
                      htop_p::Real, tsai_p::Real,
                      ivt_p::Int, nlevsoi::Int,
                      elai_p::Real, esai_p::Real,
                      fdry_p::Real, forc_rho_c::Real,
                      tgcm_p::Real;
                      medlynslope_val::Real=6.0,
                      medlynintercept_val::Real=100.0)
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
        ac[p, SUN, iv] = bsun * vcmax_z[p, SUN, iv] * smooth_max(cisun - cp_p, zero(cisun)) / (cisun + kc_p * (1.0 + oair / ko_p))
        ac[p, SHA, iv] = bsha * vcmax_z[p, SHA, iv] * smooth_max(cisha - cp_p, zero(cisha)) / (cisha + kc_p * (1.0 + oair / ko_p))
        # C3: RuBP-limited
        aj[p, SUN, iv] = jesun * smooth_max(cisun - cp_p, zero(cisun)) / (4.0 * cisun + 8.0 * cp_p)
        aj[p, SHA, iv] = jesha * smooth_max(cisha - cp_p, zero(cisha)) / (4.0 * cisha + 8.0 * cp_p)
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
        ap[p, SUN, iv] = kp_z[p, SUN, iv] * smooth_max(cisun, zero(cisun)) / forc_pbot_c
        ap[p, SHA, iv] = kp_z[p, SHA, iv] * smooth_max(cisha, zero(cisha)) / forc_pbot_c
    end

    # Gross photosynthesis - Sunlit
    aquad = params_inst.theta_cj[ivt_p]
    bquad = -(ac[p, SUN, iv] + aj[p, SUN, iv])
    cquad = ac[p, SUN, iv] * aj[p, SUN, iv]
    r1, r2 = quadratic_solve(aquad, bquad, cquad)
    ai = r2

    aquad = params_inst.theta_ip
    bquad = -(ai + ap[p, SUN, iv])
    cquad = ai * ap[p, SUN, iv]
    r1, r2 = quadratic_solve(aquad, bquad, cquad)
    ag[p, SUN, iv] = max(zero(r1), r2)

    # Gross photosynthesis - Shaded
    aquad = params_inst.theta_cj[ivt_p]
    bquad = -(ac[p, SHA, iv] + aj[p, SHA, iv])
    cquad = ac[p, SHA, iv] * aj[p, SHA, iv]
    r1, r2 = quadratic_solve(aquad, bquad, cquad)
    ai = r2

    aquad = params_inst.theta_ip
    bquad = -(ai + ap[p, SHA, iv])
    cquad = ai * ap[p, SHA, iv]
    r1, r2 = quadratic_solve(aquad, bquad, cquad)
    ag[p, SHA, iv] = max(zero(r1), r2)

    # Net photosynthesis
    an_sun[p, iv] = ag[p, SUN, iv] - bsun * lmr_z_sun
    an_sha[p, iv] = ag[p, SHA, iv] - bsha * lmr_z_sha

    if an_sun[p, iv] < 0.0
        if stomatalcond_mtd == STOMATALCOND_MTD_MEDLYN2011
            gs_mol_sun = medlynintercept_val
        elseif stomatalcond_mtd == STOMATALCOND_MTD_BB1987
            gs_mol_sun = bbb_p
        end
        gs_mol_sun = smooth_max(bsun * gs_mol_sun, 1.0)
        fvalsun = 0.0
    end
    if an_sha[p, iv] < 0.0
        if stomatalcond_mtd == STOMATALCOND_MTD_MEDLYN2011
            gs_mol_sha = medlynintercept_val
        elseif stomatalcond_mtd == STOMATALCOND_MTD_BB1987
            gs_mol_sha = bbb_p
        end
        gs_mol_sha = smooth_max(bsha * gs_mol_sha, 1.0)
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
        cs_sun = smooth_max(cs_sun, MAX_CS)
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
            gs_mol_sun = max(r1 * 1.0e06, 1.0)
        end

        # Shaded
        if an_sha[p, iv] >= 0.0
            cs_sha = cair - 1.4 / gb_mol * an_sha[p, iv] * forc_pbot_c
            cs_sha = smooth_max(cs_sha, MAX_CS)
            term = 1.6 * an_sha[p, iv] / (cs_sha / forc_pbot_c * 1.0e06)
            aquad = 1.0
            bquad = -(2.0 * (medlynintercept_val * 1.0e-06 + term) +
                      (medlynslope_val * term)^2 / (gb_mol * 1.0e-06 * rh_can))
            cquad = medlynintercept_val^2 * 1.0e-12 +
                    (2.0 * medlynintercept_val * 1.0e-06 + term *
                     (1.0 - medlynslope_val^2 / rh_can)) * term
            r1, r2 = quadratic_solve(aquad, bquad, cquad)
            gs_mol_sha = max(r1 * 1.0e06, 1.0)
        end
    elseif stomatalcond_mtd == STOMATALCOND_MTD_BB1987
        if an_sun[p, iv] >= 0.0
            aquad = cs_sun
            bsun_bbb = smooth_max(bsun * bbb_p, 1.0)
            bquad = cs_sun * (gb_mol - bsun_bbb) - mbb_p * an_sun[p, iv] * forc_pbot_c
            cquad = -gb_mol * (cs_sun * bsun_bbb + mbb_p * an_sun[p, iv] * forc_pbot_c * rh_can)
            r1, r2 = quadratic_solve(aquad, bquad, cquad)
            gs_mol_sun = max(r1, bbb_p)
        end

        # Shaded
        if an_sha[p, iv] >= 0.0
            cs_sha = cair - 1.4 / gb_mol * an_sha[p, iv] * forc_pbot_c
            cs_sha = smooth_max(cs_sha, MAX_CS)
            aquad = cs_sha
            bsha_bbb = smooth_max(bsha * bbb_p, 1.0)
            bquad = cs_sha * (gb_mol - bsha_bbb) - mbb_p * an_sha[p, iv] * forc_pbot_c
            cquad = -gb_mol * (cs_sha * bsha_bbb + mbb_p * an_sha[p, iv] * forc_pbot_c * rh_can)
            r1, r2 = quadratic_solve(aquad, bquad, cquad)
            gs_mol_sha = max(r1, bbb_p)
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

# _ci_func_PHS_core! — GPU-callable positional core of ci_func_PHS!. The vegwp
# unknown `x` is carried as 4 named scalars (xsun..xroot) and returned updated (it
# is only mutated when bflag routes to _calcstress). theta_cj_p/theta_ip_val/MAX_CS_
# threaded; ps fields (ac_phs/.../an_sun/an_sha, 3D + 2D) indexed in-place — on the
# device path `ps` is the PsnDV bundle. Returns
# (fvalsun, fvalsha, gs_mol_sun, gs_mol_sha, bsun, bsha, xsun, xsha, xxyl, xroot).
@inline function _ci_func_PHS_core!(xsun_w::Real, xsha_w::Real, xxyl_w::Real, xroot_w::Real,
        cisun::Real, cisha::Real, p::Int, iv::Int, c::Int,
        bsun::Real, bsha::Real, bflag::Bool,
        gb_mol::Real, gs0sun::Real, gs0sha::Real,
        gs_mol_sun::Real, gs_mol_sha::Real,
        jesun::Real, jesha::Real, cair::Real, oair::Real,
        lmr_z_sun::Real, lmr_z_sha::Real,
        par_z_sun::Real, par_z_sha::Real,
        rh_can::Real, qsatl::Real, qaf::Real,
        ps, forc_pbot_c::Real,
        k_soil_root, smp_l, z_col,
        laisun_p::Real, laisha_p::Real, htop_p::Real, tsai_p::Real,
        ivt_p::Int, nlevsoi::Int,
        elai_p::Real, esai_p::Real, fdry_p::Real, forc_rho_c::Real, tgcm_p::Real,
        medlynslope_val::Real, medlynintercept_val::Real,
        theta_cj_p::Real, theta_ip_val::Real, MAX_CS_::Real, RGAS_::Real,
        stomatalcond_mtd::Int, STOMATALCOND_MTD_BB1987_::Int,
        STOMATALCOND_MTD_MEDLYN2011_::Int, phs_params)
    T = typeof(cisun)
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

    fvalsun = zero(T)
    fvalsha = zero(T)

    xsun = xsun_w; xsha = xsha_w; xxyl = xxyl_w; xroot = xroot_w
    if bflag
        xsun, xsha, xxyl, xroot, bsun, bsha, _ = _calcstress(xsun, xsha, xxyl, xroot,
            gb_mol, gs0sun, gs0sha, qsatl, qaf,
            k_soil_root, smp_l, z_col, p, c, laisun_p, laisha_p, htop_p, tsai_p,
            ivt_p, nlevsoi, elai_p, esai_p, fdry_p, forc_rho_c, forc_pbot_c, tgcm_p,
            phs_params, RGAS_)
    end

    @inbounds if c3flag
        ac[p, SUN, iv] = bsun * vcmax_z[p, SUN, iv] * smooth_max(cisun - cp_p, zero(cisun)) / (cisun + kc_p * (one(T) + oair / ko_p))
        ac[p, SHA, iv] = bsha * vcmax_z[p, SHA, iv] * smooth_max(cisha - cp_p, zero(cisha)) / (cisha + kc_p * (one(T) + oair / ko_p))
        aj[p, SUN, iv] = jesun * smooth_max(cisun - cp_p, zero(cisun)) / (T(4.0) * cisun + T(8.0) * cp_p)
        aj[p, SHA, iv] = jesha * smooth_max(cisha - cp_p, zero(cisha)) / (T(4.0) * cisha + T(8.0) * cp_p)
        ap[p, SUN, iv] = T(3.0) * tpu_z[p, SUN, iv]
        ap[p, SHA, iv] = T(3.0) * tpu_z[p, SHA, iv]
    else
        ac[p, SUN, iv] = bsun * vcmax_z[p, SUN, iv]
        ac[p, SHA, iv] = bsha * vcmax_z[p, SHA, iv]
        aj[p, SUN, iv] = qe_p * par_z_sun * T(4.6)
        aj[p, SHA, iv] = qe_p * par_z_sha * T(4.6)
        ap[p, SUN, iv] = kp_z[p, SUN, iv] * smooth_max(cisun, zero(cisun)) / forc_pbot_c
        ap[p, SHA, iv] = kp_z[p, SHA, iv] * smooth_max(cisha, zero(cisha)) / forc_pbot_c
    end

    @inbounds begin
        # Gross photosynthesis - Sunlit
        aquad = theta_cj_p
        bquad = -(ac[p, SUN, iv] + aj[p, SUN, iv])
        cquad = ac[p, SUN, iv] * aj[p, SUN, iv]
        r1, r2 = quadratic_solve(aquad, bquad, cquad)
        ai = r2

        aquad = theta_ip_val
        bquad = -(ai + ap[p, SUN, iv])
        cquad = ai * ap[p, SUN, iv]
        r1, r2 = quadratic_solve(aquad, bquad, cquad)
        ag[p, SUN, iv] = max(zero(r1), r2)

        # Gross photosynthesis - Shaded
        aquad = theta_cj_p
        bquad = -(ac[p, SHA, iv] + aj[p, SHA, iv])
        cquad = ac[p, SHA, iv] * aj[p, SHA, iv]
        r1, r2 = quadratic_solve(aquad, bquad, cquad)
        ai = r2

        aquad = theta_ip_val
        bquad = -(ai + ap[p, SHA, iv])
        cquad = ai * ap[p, SHA, iv]
        r1, r2 = quadratic_solve(aquad, bquad, cquad)
        ag[p, SHA, iv] = max(zero(r1), r2)

        an_sun[p, iv] = ag[p, SUN, iv] - bsun * lmr_z_sun
        an_sha[p, iv] = ag[p, SHA, iv] - bsha * lmr_z_sha

        if an_sun[p, iv] < zero(T)
            if stomatalcond_mtd == STOMATALCOND_MTD_MEDLYN2011_
                gs_mol_sun = medlynintercept_val
            elseif stomatalcond_mtd == STOMATALCOND_MTD_BB1987_
                gs_mol_sun = bbb_p
            end
            gs_mol_sun = smooth_max(bsun * gs_mol_sun, one(T))
            fvalsun = zero(T)
        end
        if an_sha[p, iv] < zero(T)
            if stomatalcond_mtd == STOMATALCOND_MTD_MEDLYN2011_
                gs_mol_sha = medlynintercept_val
            elseif stomatalcond_mtd == STOMATALCOND_MTD_BB1987_
                gs_mol_sha = bbb_p
            end
            gs_mol_sha = smooth_max(bsha * gs_mol_sha, one(T))
            fvalsha = zero(T)
        end
        if an_sun[p, iv] < zero(T) && an_sha[p, iv] < zero(T)
            return (fvalsun, fvalsha, gs_mol_sun, gs_mol_sha, bsun, bsha, xsun, xsha, xxyl, xroot)
        end

        # Quadratic gs_mol calculation with an known
        cs_sun = zero(T)
        if an_sun[p, iv] >= zero(T)
            cs_sun = cair - T(1.4) / gb_mol * an_sun[p, iv] * forc_pbot_c
            cs_sun = smooth_max(cs_sun, MAX_CS_)
        end

        if stomatalcond_mtd == STOMATALCOND_MTD_MEDLYN2011_
            if an_sun[p, iv] >= zero(T)
                term = T(1.6) * an_sun[p, iv] / (cs_sun / forc_pbot_c * T(1.0e06))
                aquad = one(T)
                bquad = -(T(2.0) * (medlynintercept_val * T(1.0e-06) + term) +
                          (medlynslope_val * term)^2 / (gb_mol * T(1.0e-06) * rh_can))
                cquad = medlynintercept_val^2 * T(1.0e-12) +
                        (T(2.0) * medlynintercept_val * T(1.0e-06) + term *
                         (one(T) - medlynslope_val^2 / rh_can)) * term
                r1, r2 = quadratic_solve(aquad, bquad, cquad)
                gs_mol_sun = max(r1 * T(1.0e06), one(T))
            end
            if an_sha[p, iv] >= zero(T)
                cs_sha = cair - T(1.4) / gb_mol * an_sha[p, iv] * forc_pbot_c
                cs_sha = smooth_max(cs_sha, MAX_CS_)
                term = T(1.6) * an_sha[p, iv] / (cs_sha / forc_pbot_c * T(1.0e06))
                aquad = one(T)
                bquad = -(T(2.0) * (medlynintercept_val * T(1.0e-06) + term) +
                          (medlynslope_val * term)^2 / (gb_mol * T(1.0e-06) * rh_can))
                cquad = medlynintercept_val^2 * T(1.0e-12) +
                        (T(2.0) * medlynintercept_val * T(1.0e-06) + term *
                         (one(T) - medlynslope_val^2 / rh_can)) * term
                r1, r2 = quadratic_solve(aquad, bquad, cquad)
                gs_mol_sha = max(r1 * T(1.0e06), one(T))
            end
        elseif stomatalcond_mtd == STOMATALCOND_MTD_BB1987_
            if an_sun[p, iv] >= zero(T)
                aquad = cs_sun
                bsun_bbb = smooth_max(bsun * bbb_p, one(T))
                bquad = cs_sun * (gb_mol - bsun_bbb) - mbb_p * an_sun[p, iv] * forc_pbot_c
                cquad = -gb_mol * (cs_sun * bsun_bbb + mbb_p * an_sun[p, iv] * forc_pbot_c * rh_can)
                r1, r2 = quadratic_solve(aquad, bquad, cquad)
                gs_mol_sun = max(r1, bbb_p)
            end
            if an_sha[p, iv] >= zero(T)
                cs_sha = cair - T(1.4) / gb_mol * an_sha[p, iv] * forc_pbot_c
                cs_sha = smooth_max(cs_sha, MAX_CS_)
                aquad = cs_sha
                bsha_bbb = smooth_max(bsha * bbb_p, one(T))
                bquad = cs_sha * (gb_mol - bsha_bbb) - mbb_p * an_sha[p, iv] * forc_pbot_c
                cquad = -gb_mol * (cs_sha * bsha_bbb + mbb_p * an_sha[p, iv] * forc_pbot_c * rh_can)
                r1, r2 = quadratic_solve(aquad, bquad, cquad)
                gs_mol_sha = max(r1, bbb_p)
            end
        end

        # Derive new estimate for cisun, cisha
        if an_sun[p, iv] >= zero(T)
            if gs_mol_sun > zero(T)
                fvalsun = cisun - cair + an_sun[p, iv] * forc_pbot_c * (T(1.4) * gs_mol_sun + T(1.6) * gb_mol) / (gb_mol * gs_mol_sun)
            else
                fvalsun = cisun - cair
            end
        end
        if an_sha[p, iv] >= zero(T)
            if gs_mol_sha > zero(T)
                fvalsha = cisha - cair + an_sha[p, iv] * forc_pbot_c * (T(1.4) * gs_mol_sha + T(1.6) * gb_mol) / (gb_mol * gs_mol_sha)
            else
                fvalsha = cisha - cair
            end
        end
    end

    return (fvalsun, fvalsha, gs_mol_sun, gs_mol_sha, bsun, bsha, xsun, xsha, xxyl, xroot)
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
function brent_PHS!(x1sun::Real, x2sun::Real, f1sun::Real, f2sun::Real,
                    x1sha::Real, x2sha::Real, f1sha::Real, f2sha::Real,
                    tol::Real, p::Int, iv::Int, c::Int,
                    gb_mol::Real, jesun::Real, jesha::Real,
                    cair::Real, oair::Real,
                    lmr_z_sun::Real, lmr_z_sha::Real,
                    par_z_sun::Real, par_z_sha::Real,
                    rh_can::Real,
                    gs_mol_sun::Real, gs_mol_sha::Real,
                    bsun::Real, bsha::Real,
                    qsatl::Real, qaf::Real,
                    ps, forc_pbot_c::Real,
                    k_soil_root_p::AbstractVector{Float64},
                    smp_c::AbstractVector{Float64},
                    z_c::AbstractVector{Float64},
                    laisun_p::Real, laisha_p::Real,
                    htop_p::Real, tsai_p::Real,
                    ivt_p::Int, nlevsoi::Int,
                    elai_p::Real, esai_p::Real,
                    fdry_p::Real, forc_rho_c::Real,
                    tgcm_p::Real;
                    medlynslope_val::Real=6.0,
                    medlynintercept_val::Real=100.0)
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

    FT = typeof(x1sun)
    x_dummy = zeros(FT, NVEGWCS)

    # Pre-allocate Brent iteration workspace (reused each iteration)
    p_arr = zeros(FT, nphs)
    q_arr = zeros(FT, nphs)
    r_arr = zeros(FT, nphs)
    s_arr = zeros(FT, nphs)

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

        fill!(p_arr, 0.0)
        fill!(q_arr, 0.0)
        fill!(r_arr, 0.0)
        fill!(s_arr, 0.0)

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

# _brent_PHS_core! — GPU-callable positional core of brent_PHS!. The length-2
# [sun,sha] vectors are scalarized into _sun/_sha pairs; the `for phase in 1:nphs`
# loops are unrolled into explicit sun-then-sha blocks (host order). The host
# `while iter<itmax` with an early `return` becomes a fixed-trip `for iter in
# 1:itmax` + `done` flag that freezes the body (live b/gs/bsun/bsha are returned).
# bflag is constant false here → ci_func core never mutates vegwp (x passed through).
@inline function _brent_PHS_core!(
        x1sun::Real, x2sun::Real, f1sun::Real, f2sun::Real,
        x1sha::Real, x2sha::Real, f1sha::Real, f2sha::Real,
        tol::Real, p::Int, iv::Int, c::Int,
        gb_mol::Real, jesun::Real, jesha::Real, cair::Real, oair::Real,
        lmr_z_sun::Real, lmr_z_sha::Real, par_z_sun::Real, par_z_sha::Real,
        rh_can::Real, gs_mol_sun::Real, gs_mol_sha::Real, bsun::Real, bsha::Real,
        qsatl::Real, qaf::Real, ps, forc_pbot_c::Real,
        k_soil_root, smp_l, z_col,
        laisun_p::Real, laisha_p::Real, htop_p::Real, tsai_p::Real,
        ivt_p::Int, nlevsoi::Int,
        elai_p::Real, esai_p::Real, fdry_p::Real, forc_rho_c::Real, tgcm_p::Real,
        xsun_w::Real, xsha_w::Real, xxyl_w::Real, xroot_w::Real,
        medlynslope_val::Real, medlynintercept_val::Real,
        theta_cj_p::Real, theta_ip_val::Real, MAX_CS_::Real, RGAS_::Real,
        stomatalcond_mtd::Int, STOMATALCOND_MTD_BB1987_::Int,
        STOMATALCOND_MTD_MEDLYN2011_::Int, phs_params)
    T = typeof(x1sun)
    itmax = 20
    eps_val = T(1.0e-4)
    bflag_const = false

    a_sun = x1sun; a_sha = x1sha
    b_sun = x2sun; b_sha = x2sha
    fa_sun = f1sun; fa_sha = f1sha
    fb_sun = f2sun; fb_sha = f2sha

    c_sun = b_sun; c_sha = b_sha
    fc_sun = fb_sun; fc_sha = fb_sha
    d_sun = b_sun - a_sun; d_sha = b_sha - a_sha
    e_sun = d_sun; e_sha = d_sha

    done = false
    @inbounds for iter in 1:itmax
        if !done
            # phase = SUN
            if (fb_sun > zero(T) && fc_sun > zero(T)) || (fb_sun < zero(T) && fc_sun < zero(T))
                c_sun = a_sun; fc_sun = fa_sun; d_sun = b_sun - a_sun; e_sun = d_sun
            end
            if abs(fc_sun) < abs(fb_sun)
                a_sun = b_sun; b_sun = c_sun; c_sun = a_sun
                fa_sun = fb_sun; fb_sun = fc_sun; fc_sun = fa_sun
            end
            # phase = SHA
            if (fb_sha > zero(T) && fc_sha > zero(T)) || (fb_sha < zero(T) && fc_sha < zero(T))
                c_sha = a_sha; fc_sha = fa_sha; d_sha = b_sha - a_sha; e_sha = d_sha
            end
            if abs(fc_sha) < abs(fb_sha)
                a_sha = b_sha; b_sha = c_sha; c_sha = a_sha
                fa_sha = fb_sha; fb_sha = fc_sha; fc_sha = fa_sha
            end

            tol1_sun = T(2.0) * eps_val * abs(b_sun) + T(0.5) * tol
            tol1_sha = T(2.0) * eps_val * abs(b_sha) + T(0.5) * tol
            xm_sun = T(0.5) * (c_sun - b_sun)
            xm_sha = T(0.5) * (c_sha - b_sha)

            if (abs(xm_sun) <= tol1_sun || fb_sun == zero(T)) &&
               (abs(xm_sha) <= tol1_sha || fb_sha == zero(T))
                done = true
            else
                # Brent step — phase SUN
                if abs(e_sun) >= tol1_sun && abs(fa_sun) > abs(fb_sun)
                    s_sun = fb_sun / fa_sun
                    if a_sun == c_sun
                        p_sun = T(2.0) * xm_sun * s_sun
                        q_sun = one(T) - s_sun
                    else
                        q_sun = fa_sun / fc_sun
                        r_sun = fb_sun / fc_sun
                        p_sun = s_sun * (T(2.0) * xm_sun * q_sun * (q_sun - r_sun) -
                                (b_sun - a_sun) * (r_sun - one(T)))
                        q_sun = (q_sun - one(T)) * (r_sun - one(T)) * (s_sun - one(T))
                    end
                    if p_sun > zero(T)
                        q_sun = -q_sun
                    end
                    p_sun = abs(p_sun)
                    if T(2.0) * p_sun < min(T(3.0) * xm_sun * q_sun - abs(tol1_sun * q_sun),
                                            abs(e_sun * q_sun))
                        e_sun = d_sun
                        d_sun = p_sun / q_sun
                    else
                        d_sun = xm_sun
                        e_sun = d_sun
                    end
                else
                    d_sun = xm_sun
                    e_sun = d_sun
                end
                a_sun = b_sun
                fa_sun = fb_sun
                if abs(d_sun) > tol1_sun
                    b_sun = b_sun + d_sun
                else
                    b_sun = b_sun + copysign(tol1_sun, xm_sun)
                end

                # Brent step — phase SHA
                if abs(e_sha) >= tol1_sha && abs(fa_sha) > abs(fb_sha)
                    s_sha = fb_sha / fa_sha
                    if a_sha == c_sha
                        p_sha = T(2.0) * xm_sha * s_sha
                        q_sha = one(T) - s_sha
                    else
                        q_sha = fa_sha / fc_sha
                        r_sha = fb_sha / fc_sha
                        p_sha = s_sha * (T(2.0) * xm_sha * q_sha * (q_sha - r_sha) -
                                (b_sha - a_sha) * (r_sha - one(T)))
                        q_sha = (q_sha - one(T)) * (r_sha - one(T)) * (s_sha - one(T))
                    end
                    if p_sha > zero(T)
                        q_sha = -q_sha
                    end
                    p_sha = abs(p_sha)
                    if T(2.0) * p_sha < min(T(3.0) * xm_sha * q_sha - abs(tol1_sha * q_sha),
                                            abs(e_sha * q_sha))
                        e_sha = d_sha
                        d_sha = p_sha / q_sha
                    else
                        d_sha = xm_sha
                        e_sha = d_sha
                    end
                else
                    d_sha = xm_sha
                    e_sha = d_sha
                end
                a_sha = b_sha
                fa_sha = fb_sha
                if abs(d_sha) > tol1_sha
                    b_sha = b_sha + d_sha
                else
                    b_sha = b_sha + copysign(tol1_sha, xm_sha)
                end

                fb_sun, fb_sha, gs_mol_sun, gs_mol_sha, bsun, bsha, _, _, _, _ =
                    _ci_func_PHS_core!(xsun_w, xsha_w, xxyl_w, xroot_w,
                        b_sun, b_sha, p, iv, c, bsun, bsha, bflag_const,
                        gb_mol, gs_mol_sun, gs_mol_sha, gs_mol_sun, gs_mol_sha,
                        jesun, jesha, cair, oair,
                        lmr_z_sun, lmr_z_sha, par_z_sun, par_z_sha,
                        rh_can, qsatl, qaf, ps, forc_pbot_c,
                        k_soil_root, smp_l, z_col,
                        laisun_p, laisha_p, htop_p, tsai_p, ivt_p, nlevsoi,
                        elai_p, esai_p, fdry_p, forc_rho_c, tgcm_p,
                        medlynslope_val, medlynintercept_val,
                        theta_cj_p, theta_ip_val, MAX_CS_, RGAS_,
                        stomatalcond_mtd, STOMATALCOND_MTD_BB1987_,
                        STOMATALCOND_MTD_MEDLYN2011_, phs_params)

                if fb_sun == zero(T) && fb_sha == zero(T)
                    done = true
                end
            end
        end
    end

    return (b_sun, b_sha, gs_mol_sun, gs_mol_sha, bsun, bsha)
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
function hybrid_PHS!(x0sun::Real, x0sha::Real,
                     p::Int, iv::Int, c::Int,
                     gb_mol::Real,
                     jesun::Real, jesha::Real,
                     cair::Real, oair::Real,
                     lmr_z_sun::Real, lmr_z_sha::Real,
                     par_z_sun::Real, par_z_sha::Real,
                     rh_can::Real, qsatl::Real, qaf::Real,
                     ps, forc_pbot_c::Real,
                     vegwp_p::AbstractVector{<:Real},
                     k_soil_root_p::AbstractVector{Float64},
                     smp_c::AbstractVector{Float64},
                     z_c::AbstractVector{Float64},
                     laisun_p::Real, laisha_p::Real,
                     htop_p::Real, tsai_p::Real,
                     ivt_p::Int, nlevsoi::Int,
                     elai_p::Real, esai_p::Real,
                     fdry_p::Real, forc_rho_c::Real,
                     tgcm_p::Real;
                     medlynslope_val::Real=6.0,
                     medlynintercept_val::Real=100.0)
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
    iter2 = 0          # hoisted out of the while body so it's bound for the return tuple
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

# _hybrid_PHS_core! — GPU-callable positional core of hybrid_PHS!. The vegwp input
# (4 named scalars vegwp0_*) is the per-iteration seed `x` (host: x=copy(vegwp_p));
# its 4 components are reset each outer iteration and threaded through the ci_func
# core. Both host `while true` loops become fixed-trip `for`s + done/converged flags
# (uniform trip count, no divergent break). Final vegwp is recomputed via _getvegwp.
# Returns (x0sun, x0sha, gs_mol_sun, gs_mol_sha, bsun, bsha, vsun, vsha, vxyl, vroot,
# soilflux). All Float64 literals are T()-wrapped; medlyn/theta/MAX_CS threaded.
@inline function _hybrid_PHS_core!(
        x0sun::Real, x0sha::Real, p::Int, iv::Int, c::Int,
        gb_mol::Real, jesun::Real, jesha::Real, cair::Real, oair::Real,
        lmr_z_sun::Real, lmr_z_sha::Real, par_z_sun::Real, par_z_sha::Real,
        rh_can::Real, qsatl::Real, qaf::Real, ps, forc_pbot_c::Real,
        vegwp0_sun::Real, vegwp0_sha::Real, vegwp0_xyl::Real, vegwp0_root::Real,
        k_soil_root, smp_l, z_col,
        laisun_p::Real, laisha_p::Real, htop_p::Real, tsai_p::Real,
        ivt_p::Int, nlevsoi::Int,
        elai_p::Real, esai_p::Real, fdry_p::Real, forc_rho_c::Real, tgcm_p::Real,
        medlynslope_val::Real, medlynintercept_val::Real,
        theta_cj_p::Real, theta_ip_val::Real, MAX_CS_::Real, RGAS_::Real,
        stomatalcond_mtd::Int, STOMATALCOND_MTD_BB1987_::Int,
        STOMATALCOND_MTD_MEDLYN2011_::Int, phs_params)
    T = typeof(x0sun)
    toldb = T(1.0e-2)
    eps_val = T(1.0e-2)
    eps1 = T(1.0e-4)
    itmax = 3

    x1sun = x0sun
    x1sha = x0sha
    bflag = false
    b0sun = -one(T)
    b0sha = -one(T)
    gs0sun = zero(T)
    gs0sha = zero(T)
    bsun = one(T)
    bsha = one(T)
    gs_mol_sun = zero(T)
    gs_mol_sha = zero(T)
    dbsun = zero(T)
    dbsha = zero(T)

    outer_done = false
    @inbounds for io in 1:(itmax + 1)
        if !outer_done
            # x = copy(vegwp_p): reset the working vegwp scalars each outer iter
            xsun_w = vegwp0_sun; xsha_w = vegwp0_sha; xxyl_w = vegwp0_xyl; xroot_w = vegwp0_root

            x0sun = max(T(0.1), x1sun)
            x1sun = T(0.99) * x1sun
            x0sha = max(T(0.1), x1sha)
            x1sha = T(0.99) * x1sha
            tolsun = abs(x1sun) * eps_val
            tolsha = abs(x1sha) * eps_val

            # First ci_func_PHS call (updates bsun/bsha except first iter; may mutate x)
            f0sun, f0sha, gs_mol_sun, gs_mol_sha, bsun, bsha, xsun_w, xsha_w, xxyl_w, xroot_w =
                _ci_func_PHS_core!(xsun_w, xsha_w, xxyl_w, xroot_w,
                    x0sun, x0sha, p, iv, c, bsun, bsha, bflag,
                    gb_mol, gs0sun, gs0sha, gs_mol_sun, gs_mol_sha,
                    jesun, jesha, cair, oair, lmr_z_sun, lmr_z_sha,
                    par_z_sun, par_z_sha, rh_can, qsatl, qaf, ps, forc_pbot_c,
                    k_soil_root, smp_l, z_col, laisun_p, laisha_p, htop_p, tsai_p,
                    ivt_p, nlevsoi, elai_p, esai_p, fdry_p, forc_rho_c, tgcm_p,
                    medlynslope_val, medlynintercept_val,
                    theta_cj_p, theta_ip_val, MAX_CS_, RGAS_,
                    stomatalcond_mtd, STOMATALCOND_MTD_BB1987_,
                    STOMATALCOND_MTD_MEDLYN2011_, phs_params)

            dbsun = b0sun - bsun
            dbsha = b0sha - bsha
            b0sun = bsun
            b0sha = bsha
            bflag = false

            # Second ci_func_PHS call (second point)
            f1sun, f1sha, gs_mol_sun, gs_mol_sha, bsun, bsha, xsun_w, xsha_w, xxyl_w, xroot_w =
                _ci_func_PHS_core!(xsun_w, xsha_w, xxyl_w, xroot_w,
                    x1sun, x1sha, p, iv, c, bsun, bsha, bflag,
                    gb_mol, gs0sun, gs0sha, gs_mol_sun, gs_mol_sha,
                    jesun, jesha, cair, oair, lmr_z_sun, lmr_z_sha,
                    par_z_sun, par_z_sha, rh_can, qsatl, qaf, ps, forc_pbot_c,
                    k_soil_root, smp_l, z_col, laisun_p, laisha_p, htop_p, tsai_p,
                    ivt_p, nlevsoi, elai_p, esai_p, fdry_p, forc_rho_c, tgcm_p,
                    medlynslope_val, medlynintercept_val,
                    theta_cj_p, theta_ip_val, MAX_CS_, RGAS_,
                    stomatalcond_mtd, STOMATALCOND_MTD_BB1987_,
                    STOMATALCOND_MTD_MEDLYN2011_, phs_params)

            minf = abs(f1sun + f1sha)
            minxsun = x1sun
            minxsha = x1sha

            iter2 = 0
            inner_done = false
            for ii in 1:(itmax + 2)
                if !inner_done
                    if abs(f0sun) < eps1 && abs(f0sha) < eps1
                        x1sun = x0sun
                        x1sha = x0sha
                        inner_done = true
                    elseif abs(f1sun) < eps1 && abs(f1sha) < eps1
                        inner_done = true
                    else
                        iter2 += 1
                        if (f1sun - f0sun) == zero(T)
                            dxsun = T(0.5) * (x1sun + x0sun) - x1sun
                        else
                            dxsun = -f1sun * (x1sun - x0sun) / (f1sun - f0sun)
                        end
                        if (f1sha - f0sha) == zero(T)
                            dxsha = T(0.5) * (x1sha + x0sha) - x1sha
                        else
                            dxsha = -f1sha * (x1sha - x0sha) / (f1sha - f0sha)
                        end
                        x0sun = x1sun
                        x1sun = x1sun + dxsun
                        x0sha = x1sha
                        x1sha = x1sha + dxsha

                        f1sun, f1sha, gs_mol_sun, gs_mol_sha, bsun, bsha, xsun_w, xsha_w, xxyl_w, xroot_w =
                            _ci_func_PHS_core!(xsun_w, xsha_w, xxyl_w, xroot_w,
                                x1sun, x1sha, p, iv, c, bsun, bsha, bflag,
                                gb_mol, gs0sun, gs0sha, gs_mol_sun, gs_mol_sha,
                                jesun, jesha, cair, oair, lmr_z_sun, lmr_z_sha,
                                par_z_sun, par_z_sha, rh_can, qsatl, qaf, ps, forc_pbot_c,
                                k_soil_root, smp_l, z_col, laisun_p, laisha_p, htop_p, tsai_p,
                                ivt_p, nlevsoi, elai_p, esai_p, fdry_p, forc_rho_c, tgcm_p,
                                medlynslope_val, medlynintercept_val,
                                theta_cj_p, theta_ip_val, MAX_CS_, RGAS_,
                                stomatalcond_mtd, STOMATALCOND_MTD_BB1987_,
                                STOMATALCOND_MTD_MEDLYN2011_, phs_params)

                        if abs(dxsun) < tolsun && abs(dxsha) < tolsha
                            x0sun = x1sun
                            x0sha = x1sha
                            inner_done = true
                        else
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
                                inner_done = true
                            elseif f1sun * f0sun < zero(T) && f1sha * f0sha < zero(T)
                                xsun, xsha, gs_mol_sun, gs_mol_sha, bsun, bsha = _brent_PHS_core!(
                                    x0sun, x1sun, f0sun, f1sun, x0sha, x1sha, f0sha, f1sha,
                                    tolsun, p, iv, c, gb_mol, jesun, jesha, cair, oair,
                                    lmr_z_sun, lmr_z_sha, par_z_sun, par_z_sha, rh_can,
                                    gs_mol_sun, gs_mol_sha, bsun, bsha, qsatl, qaf, ps, forc_pbot_c,
                                    k_soil_root, smp_l, z_col, laisun_p, laisha_p, htop_p, tsai_p,
                                    ivt_p, nlevsoi, elai_p, esai_p, fdry_p, forc_rho_c, tgcm_p,
                                    xsun_w, xsha_w, xxyl_w, xroot_w,
                                    medlynslope_val, medlynintercept_val,
                                    theta_cj_p, theta_ip_val, MAX_CS_, RGAS_,
                                    stomatalcond_mtd, STOMATALCOND_MTD_BB1987_,
                                    STOMATALCOND_MTD_MEDLYN2011_, phs_params)
                                x0sun = xsun
                                x0sha = xsha
                                inner_done = true
                            elseif iter2 > itmax
                                x1sun = minxsun
                                x1sha = minxsha
                                f1sun, f1sha, gs_mol_sun, gs_mol_sha, bsun, bsha, xsun_w, xsha_w, xxyl_w, xroot_w =
                                    _ci_func_PHS_core!(xsun_w, xsha_w, xxyl_w, xroot_w,
                                        x1sun, x1sha, p, iv, c, bsun, bsha, bflag,
                                        gb_mol, gs0sun, gs0sha, gs_mol_sun, gs_mol_sha,
                                        jesun, jesha, cair, oair, lmr_z_sun, lmr_z_sha,
                                        par_z_sun, par_z_sha, rh_can, qsatl, qaf, ps, forc_pbot_c,
                                        k_soil_root, smp_l, z_col, laisun_p, laisha_p, htop_p, tsai_p,
                                        ivt_p, nlevsoi, elai_p, esai_p, fdry_p, forc_rho_c, tgcm_p,
                                        medlynslope_val, medlynintercept_val,
                                        theta_cj_p, theta_ip_val, MAX_CS_, RGAS_,
                                        stomatalcond_mtd, STOMATALCOND_MTD_BB1987_,
                                        STOMATALCOND_MTD_MEDLYN2011_, phs_params)
                                inner_done = true
                            end
                        end
                    end
                end
            end  # inner loop

            # Update unstressed stomatal conductance
            if bsun > T(0.01)
                gs0sun = gs_mol_sun / bsun
            end
            if bsha > T(0.01)
                gs0sha = gs_mol_sha / bsha
            end
            bflag = true

            if abs(dbsun) < toldb && abs(dbsha) < toldb
                outer_done = true
            elseif io > itmax
                outer_done = true
            end
        end
    end  # outer loop

    x0sun = x1sun
    x0sha = x1sha

    # Final vegwp (algebraic), recomputed from the converged gs_mol
    vsun, vsha, vxyl, vroot, soilflux = _getvegwp(gb_mol, gs_mol_sun, gs_mol_sha,
        qsatl, qaf, k_soil_root, smp_l, z_col, p, c, laisun_p, laisha_p, htop_p,
        tsai_p, ivt_p, nlevsoi, elai_p, esai_p, fdry_p, forc_rho_c, forc_pbot_c,
        tgcm_p, phs_params, RGAS_)

    if soilflux < zero(T)
        soilflux = zero(T)
    end

    return (x0sun, x0sha, gs_mol_sun, gs_mol_sha, bsun, bsha,
            vsun, vsha, vxyl, vroot, soilflux)
end

# =====================================================================
# PHS Pass 3 kernel — per-patch leaf photosynthesis + the simultaneous
# sun/shade plant-hydraulic-stress Newton solve, run ENTIRELY in-thread.
# ONE THREAD PER PATCH (p = @index(Global) over ndrange = length(bounds_patch));
# the inner `for iv in 1:nrad[p]` loop runs serially per thread (per-patch
# outputs vegwp[p,:]/qflx/rh_leaf/gb_mol are written within the iv loop, so
# flattening to (p,iv) would race — per-patch is race-free, matches serial order).
# Calls the GPU-callable PHS chain (_hybrid_PHS_core! → _ci_func_PHS_core! →
# _brent_PHS_core!/_calcstress → _spacF/_spacA/_getvegwp/_getqflx), which read
# NO globals. The vegwp length-4 unknown is carried as 4 named scalars; the
# nlevsoi soil stack (k_soil_root[p,·]/smp_l[c,·]/z_col[c,·]) is indexed in-kernel
# (full 2D arrays + (p,c), no device slices).
# =====================================================================
# Scalar bundle (Metal ~31-arg limit). theta_cj[ivt] comes from phs_params.theta_cj.
struct PsnPhsP3Scalars{S}
    fnps::S; theta_psii::S; theta_ip::S; medlyn_slope_override::S
    RGAS_::S; MAX_CS_::S; MEDLYN_RH_CAN_MAX_::S; MEDLYN_RH_CAN_FACT_::S
    rsmax0::S; SPVAL_::S
end

# --- PHS Pass-3 device-view arg bundles (Metal 31-arg-buffer reduction) ---
# The fused-Newton kernel has ~48 loose arrays; grouped into 4 @adapt_structure'd
# bundles (each = 1 indirect-arg-buffer, ≤17 array fields each). Fields aliased to
# locals at the kernel top so the physics body is byte-identical.
# Index/Bool arrays own params (VI=Int col/ivt/nrad, VB=mask+near-noon Bool); PFT
# params Float64-pinned in the AD path → own param (Vp).
Base.@kwdef struct _Psn3Idx{VB,VI,Vp}
    mask_patch::VB; is_near_local_noon::VB
    col_of_patch::VI; ivt::VI; nrad::VI
    medlynslope_pft::Vp; medlynintercept_pft::Vp; crop_pft::Vp
end
Base.@kwdef struct _Psn3InA{V}      # read-only per-patch/col forcing vectors (1/1)
    forc_pbot::V; tgcm::V; rb::V; eair::V; esat_tv::V; cair::V; oair::V; qsatl::V; qaf::V
    laisun::V; laisha::V; htop::V; tsai::V; elai::V; esai::V; fdry::V; forc_rho::V
end
Base.@kwdef struct _Psn3InB{V,M,M3}  # o3 vectors (V) + 2D matrices (M) + 3D jmax (M3)
    o3coefg_sun::V; o3coefg_sha::V; o3coefv_sun::V; o3coefv_sha::V
    par_z_sun_in::M; par_z_sha_in::M
    k_soil_root::M; smp_l::M; z_col::M
    jmax_z_local::M3
end
Base.@kwdef struct _Psn3Out{V,M}    # written outputs: vectors (V) + matrices (M)
    qflx_tran_veg::V; bsun_arr::V; bsha_arr::V; rh_leaf_sun::V; rh_leaf_sha::V
    vegwp::M; vegwp_ln::M
    psn_wc_z_sun::M; psn_wj_z_sun::M; psn_wp_z_sun::M
    psn_wc_z_sha::M; psn_wj_z_sha::M; psn_wp_z_sha::M
end
Adapt.@adapt_structure _Psn3Idx
Adapt.@adapt_structure _Psn3InA
Adapt.@adapt_structure _Psn3InB
Adapt.@adapt_structure _Psn3Out

@kernel function _psn_phs_pass3_kernel!(ps, phs_params, sc, idx, ina, inb, out,
        nlevsoi::Int, modifyphoto_and_lmr_forcrop::Bool, stomatalcond_mtd::Int,
        STOMATALCOND_MTD_BB1987_::Int, STOMATALCOND_MTD_MEDLYN2011_::Int)
    p = @index(Global)
    _psn_phs_pass3_body!(p, ps, phs_params, sc, idx, ina, inb, out, nlevsoi,
        modifyphoto_and_lmr_forcrop, stomatalcond_mtd,
        STOMATALCOND_MTD_BB1987_, STOMATALCOND_MTD_MEDLYN2011_)
end

# Shared per-patch body for PHS Pass 3 (KA kernel + AD host loop) — see _photosynth_ci_body!.
# This is the fused per-patch PHS Newton solve (hybrid_PHS → ci_func_PHS →
# brent_PHS/calcstress → spacF/spacA/getvegwp/getqflx via fixed-trip device cores).
@inline function _psn_phs_pass3_body!(p, ps, phs_params, sc, idx, ina, inb, out,
        nlevsoi::Int, modifyphoto_and_lmr_forcrop::Bool, stomatalcond_mtd::Int,
        STOMATALCOND_MTD_BB1987_::Int, STOMATALCOND_MTD_MEDLYN2011_::Int)
    # alias grouped fields to locals → physics body below is byte-identical
    mask_patch = idx.mask_patch; is_near_local_noon = idx.is_near_local_noon
    col_of_patch = idx.col_of_patch; ivt = idx.ivt; nrad = idx.nrad
    medlynslope_pft = idx.medlynslope_pft; medlynintercept_pft = idx.medlynintercept_pft
    crop_pft = idx.crop_pft
    forc_pbot = ina.forc_pbot; tgcm = ina.tgcm; rb = ina.rb; eair = ina.eair
    esat_tv = ina.esat_tv; cair = ina.cair; oair = ina.oair; qsatl = ina.qsatl
    qaf = ina.qaf; laisun = ina.laisun; laisha = ina.laisha; htop = ina.htop
    tsai = ina.tsai; elai = ina.elai; esai = ina.esai; fdry = ina.fdry
    forc_rho = ina.forc_rho
    o3coefg_sun = inb.o3coefg_sun; o3coefg_sha = inb.o3coefg_sha
    o3coefv_sun = inb.o3coefv_sun; o3coefv_sha = inb.o3coefv_sha
    par_z_sun_in = inb.par_z_sun_in; par_z_sha_in = inb.par_z_sha_in
    jmax_z_local = inb.jmax_z_local; k_soil_root = inb.k_soil_root
    smp_l = inb.smp_l; z_col = inb.z_col
    qflx_tran_veg = out.qflx_tran_veg; bsun_arr = out.bsun_arr; bsha_arr = out.bsha_arr
    rh_leaf_sun = out.rh_leaf_sun; rh_leaf_sha = out.rh_leaf_sha
    vegwp = out.vegwp; vegwp_ln = out.vegwp_ln
    psn_wc_z_sun = out.psn_wc_z_sun; psn_wj_z_sun = out.psn_wj_z_sun
    psn_wp_z_sun = out.psn_wp_z_sun; psn_wc_z_sha = out.psn_wc_z_sha
    psn_wj_z_sha = out.psn_wj_z_sha; psn_wp_z_sha = out.psn_wp_z_sha
    @inbounds if mask_patch[p]
        T = eltype(forc_pbot)
        c = col_of_patch[p]
        ivt_p = ivt[p]

        fnps = sc.fnps; theta_psii = sc.theta_psii; theta_ip_val = sc.theta_ip
        RGAS_ = sc.RGAS_; MAX_CS_ = sc.MAX_CS_; rsmax0 = sc.rsmax0; SPVAL_ = sc.SPVAL_
        MEDLYN_RH_CAN_MAX_ = sc.MEDLYN_RH_CAN_MAX_; MEDLYN_RH_CAN_FACT_ = sc.MEDLYN_RH_CAN_FACT_

        theta_cj_p = phs_params.theta_cj[ivt_p]
        # Medlyn slope: use override if set, else PFT default
        medlynslope_p = isnan(sc.medlyn_slope_override) ? medlynslope_pft[ivt_p] : sc.medlyn_slope_override
        medlynintercept_p = medlynintercept_pft[ivt_p]

        ac = ps.ac_phs_patch; aj = ps.aj_phs_patch; ap = ps.ap_phs_patch; ag = ps.ag_phs_patch
        an_sun = ps.an_sun_patch; an_sha = ps.an_sha_patch
        gs_mol_sun = ps.gs_mol_sun_patch; gs_mol_sha = ps.gs_mol_sha_patch
        gs_mol_sun_ln = ps.gs_mol_sun_ln_patch; gs_mol_sha_ln = ps.gs_mol_sha_ln_patch
        bbb = ps.bbb_patch; vpd_can = ps.vpd_can_patch; c3flag = ps.c3flag_patch
        gb_mol_arr = ps.gb_mol_patch
        ci_z_sun = ps.cisun_z_patch; ci_z_sha = ps.cisha_z_patch
        rs_z_sun = ps.rssun_z_patch; rs_z_sha = ps.rssha_z_patch
        lmr_z_sun = ps.lmrsun_z_patch; lmr_z_sha = ps.lmrsha_z_patch
        psn_z_sun = ps.psnsun_z_patch; psn_z_sha = ps.psnsha_z_patch

        is_crop = crop_pft[ivt_p] > T(0.5)  # round(Int,·)==0 ⟺ <0.5 (0/1 valued)

        cf = forc_pbot[p] / (RGAS_ * tgcm[p]) * T(1.0e06)
        gb = one(T) / rb[p]
        gb_mol_arr[p] = gb * cf

        for iv in 1:nrad[p]
            if par_z_sun_in[p, iv] <= zero(T)  # night time
                vegwp[p, SUN] = one(T)  # signal for night

                if stomatalcond_mtd == STOMATALCOND_MTD_BB1987_
                    gsminsun = bbb[p]
                    gsminsha = bbb[p]
                else
                    gsminsun = medlynintercept_p
                    gsminsha = medlynintercept_p
                end

                vsun, vsha, vxyl, vroot, bsun_val, bsha_val, soilflux_night = _calcstress(
                    vegwp[p, SUN], vegwp[p, SHA], vegwp[p, XYL], vegwp[p, ROOT_SEG],
                    gb_mol_arr[p], gsminsun, gsminsha, qsatl[p], qaf[p],
                    k_soil_root, smp_l, z_col, p, c,
                    laisun[p], laisha[p], htop[p], tsai[p], ivt_p, nlevsoi,
                    elai[p], esai[p], fdry[p], forc_rho[c], forc_pbot[p], tgcm[p],
                    phs_params, RGAS_)
                vegwp[p, SUN] = vsun; vegwp[p, SHA] = vsha
                vegwp[p, XYL] = vxyl; vegwp[p, ROOT_SEG] = vroot
                bsun_arr[p] = bsun_val
                bsha_arr[p] = bsha_val
                # Night transpiration: Fortran sets qflx_tran_veg(p)=soilflux inside the
                # night calcstress->getvegwp path; without this it stays at its NaN init.
                qflx_tran_veg[p] = max(soilflux_night, zero(T))   # HARD: transpiration is kg/m2/s. smooth_max(x,0) sits AT its kink here and returned log(2)/50 = 0.0139 mm/s = 1200 mm/day of water no store paid for.

                ac[p, SUN, iv] = zero(T)
                aj[p, SUN, iv] = zero(T)
                ap[p, SUN, iv] = zero(T)
                ag[p, SUN, iv] = zero(T)
                if !is_crop || !modifyphoto_and_lmr_forcrop
                    an_sun[p, iv] = ag[p, SUN, iv] - bsun_arr[p] * lmr_z_sun[p, iv]
                else
                    an_sun[p, iv] = ag[p, SUN, iv] - lmr_z_sun[p, iv]
                end
                psn_z_sun[p, iv] = zero(T)
                psn_wc_z_sun[p, iv] = zero(T)
                psn_wj_z_sun[p, iv] = zero(T)
                psn_wp_z_sun[p, iv] = zero(T)
                rs_z_sun[p, iv] = smooth_min(rsmax0, one(T) / smooth_max(bsun_arr[p] * gsminsun, one(T)) * cf)
                ci_z_sun[p, iv] = zero(T)
                rh_leaf_sun[p] = zero(T)

                ac[p, SHA, iv] = zero(T)
                aj[p, SHA, iv] = zero(T)
                ap[p, SHA, iv] = zero(T)
                ag[p, SHA, iv] = zero(T)
                if !is_crop || !modifyphoto_and_lmr_forcrop
                    an_sha[p, iv] = ag[p, SHA, iv] - bsha_arr[p] * lmr_z_sha[p, iv]
                else
                    an_sha[p, iv] = ag[p, SHA, iv] - lmr_z_sha[p, iv]
                end
                psn_z_sha[p, iv] = zero(T)
                psn_wc_z_sha[p, iv] = zero(T)
                psn_wj_z_sha[p, iv] = zero(T)
                psn_wp_z_sha[p, iv] = zero(T)
                rs_z_sha[p, iv] = smooth_min(rsmax0, one(T) / smooth_max(bsha_arr[p] * gsminsha, one(T)) * cf)
                ci_z_sha[p, iv] = zero(T)
                rh_leaf_sha[p] = zero(T)

            else  # day time
                ceair = smooth_min(eair[p], esat_tv[p])
                if stomatalcond_mtd == STOMATALCOND_MTD_BB1987_
                    rh_can = ceair / esat_tv[p]
                else
                    rh_can = smooth_max((esat_tv[p] - ceair), MEDLYN_RH_CAN_MAX_) * MEDLYN_RH_CAN_FACT_
                    vpd_can[p] = rh_can
                end

                # Electron transport - Sun
                qabs = T(0.5) * (one(T) - fnps) * par_z_sun_in[p, iv] * T(4.6)
                aquad = theta_psii
                bquad = -(qabs + jmax_z_local[p, SUN, iv])
                cquad = qabs * jmax_z_local[p, SUN, iv]
                r1, r2 = quadratic_solve(aquad, bquad, cquad)
                je_sun = r2

                # Electron transport - Shade
                qabs = T(0.5) * (one(T) - fnps) * par_z_sha_in[p, iv] * T(4.6)
                aquad = theta_psii
                bquad = -(qabs + jmax_z_local[p, SHA, iv])
                cquad = qabs * jmax_z_local[p, SHA, iv]
                r1, r2 = quadratic_solve(aquad, bquad, cquad)
                je_sha = r2

                # Initial ci guess
                if c3flag[p]
                    ci_z_sun[p, iv] = T(0.7) * cair[p]
                    ci_z_sha[p, iv] = T(0.7) * cair[p]
                else
                    ci_z_sun[p, iv] = T(0.4) * cair[p]
                    ci_z_sha[p, iv] = T(0.4) * cair[p]
                end

                # Solve for ci and gs via the in-thread hybrid_PHS chain
                cisun_sol, cisha_sol, gs_mol_sun_val, gs_mol_sha_val,
                    bsun_val, bsha_val, vsun, vsha, vxyl, vroot, soilflux =
                    _hybrid_PHS_core!(ci_z_sun[p, iv], ci_z_sha[p, iv], p, iv, c,
                        gb_mol_arr[p], je_sun, je_sha, cair[p], oair[p],
                        lmr_z_sun[p, iv], lmr_z_sha[p, iv],
                        par_z_sun_in[p, iv], par_z_sha_in[p, iv],
                        rh_can, qsatl[p], qaf[p], ps, forc_pbot[p],
                        vegwp[p, SUN], vegwp[p, SHA], vegwp[p, XYL], vegwp[p, ROOT_SEG],
                        k_soil_root, smp_l, z_col,
                        laisun[p], laisha[p], htop[p], tsai[p], ivt_p, nlevsoi,
                        elai[p], esai[p], fdry[p], forc_rho[c], tgcm[p],
                        medlynslope_p, medlynintercept_p,
                        theta_cj_p, theta_ip_val, MAX_CS_, RGAS_,
                        stomatalcond_mtd, STOMATALCOND_MTD_BB1987_,
                        STOMATALCOND_MTD_MEDLYN2011_, phs_params)
                ci_z_sun[p, iv] = cisun_sol
                ci_z_sha[p, iv] = cisha_sol
                vegwp[p, SUN] = vsun; vegwp[p, SHA] = vsha
                vegwp[p, XYL] = vxyl; vegwp[p, ROOT_SEG] = vroot
                gs_mol_sun[p, iv] = gs_mol_sun_val
                gs_mol_sha[p, iv] = gs_mol_sha_val
                bsun_arr[p] = bsun_val
                bsha_arr[p] = bsha_val
                qflx_tran_veg[p] = max(soilflux, zero(soilflux))   # HARD: see the night branch — mm/s axis, ReLU at its kink.

                # gs min for error checking
                if stomatalcond_mtd == STOMATALCOND_MTD_MEDLYN2011_
                    gsminsun = medlynintercept_p
                    gsminsha = medlynintercept_p
                else
                    gsminsun = bbb[p]
                    gsminsha = bbb[p]
                end

                if an_sun[p, iv] < zero(T)
                    gs_mol_sun[p, iv] = smooth_max(bsun_arr[p] * gsminsun, one(T))
                end
                if an_sha[p, iv] < zero(T)
                    gs_mol_sha[p, iv] = smooth_max(bsha_arr[p] * gsminsha, one(T))
                end

                # Local noon gs
                if is_near_local_noon[p]
                    gs_mol_sun_ln[p, iv] = gs_mol_sun[p, iv]
                    gs_mol_sha_ln[p, iv] = gs_mol_sha[p, iv]
                    vegwp_ln[p, SUN] = vegwp[p, SUN]; vegwp_ln[p, SHA] = vegwp[p, SHA]
                    vegwp_ln[p, XYL] = vegwp[p, XYL]; vegwp_ln[p, ROOT_SEG] = vegwp[p, ROOT_SEG]
                else
                    gs_mol_sun_ln[p, iv] = SPVAL_
                    gs_mol_sha_ln[p, iv] = SPVAL_
                    vegwp_ln[p, SUN] = SPVAL_; vegwp_ln[p, SHA] = SPVAL_
                    vegwp_ln[p, XYL] = SPVAL_; vegwp_ln[p, ROOT_SEG] = SPVAL_
                end

                # Final cs and ci
                cs_sun = cair[p] - T(1.4) / gb_mol_arr[p] * an_sun[p, iv] * forc_pbot[p]
                cs_sun = smooth_max(cs_sun, MAX_CS_)
                ci_z_sun[p, iv] = cair[p] - an_sun[p, iv] * forc_pbot[p] *
                    (T(1.4) * gs_mol_sun[p, iv] + T(1.6) * gb_mol_arr[p]) /
                    (gb_mol_arr[p] * gs_mol_sun[p, iv])
                ci_z_sun[p, iv] = smooth_max(ci_z_sun[p, iv], T(1.0e-06))

                cs_sha = cair[p] - T(1.4) / gb_mol_arr[p] * an_sha[p, iv] * forc_pbot[p]
                cs_sha = smooth_max(cs_sha, MAX_CS_)
                ci_z_sha[p, iv] = cair[p] - an_sha[p, iv] * forc_pbot[p] *
                    (T(1.4) * gs_mol_sha[p, iv] + T(1.6) * gb_mol_arr[p]) /
                    (gb_mol_arr[p] * gs_mol_sha[p, iv])
                ci_z_sha[p, iv] = smooth_max(ci_z_sha[p, iv], T(1.0e-06))

                # Convert to resistance
                gs = max(gs_mol_sun[p, iv], one(T)) / cf
                rs_z_sun[p, iv] = smooth_min(one(T) / gs, rsmax0)
                rs_z_sun[p, iv] = rs_z_sun[p, iv] / o3coefg_sun[p]
                gs = max(gs_mol_sha[p, iv], one(T)) / cf
                rs_z_sha[p, iv] = smooth_min(one(T) / gs, rsmax0)
                rs_z_sha[p, iv] = rs_z_sha[p, iv] / o3coefg_sha[p]

                # Photosynthesis output
                psn_z_sun[p, iv] = ag[p, SUN, iv] * o3coefv_sun[p]
                psn_wc_z_sun[p, iv] = zero(T)
                psn_wj_z_sun[p, iv] = zero(T)
                psn_wp_z_sun[p, iv] = zero(T)
                if ac[p, SUN, iv] <= aj[p, SUN, iv] && ac[p, SUN, iv] <= ap[p, SUN, iv]
                    psn_wc_z_sun[p, iv] = psn_z_sun[p, iv]
                elseif aj[p, SUN, iv] < ac[p, SUN, iv] && aj[p, SUN, iv] <= ap[p, SUN, iv]
                    psn_wj_z_sun[p, iv] = psn_z_sun[p, iv]
                elseif ap[p, SUN, iv] < ac[p, SUN, iv] && ap[p, SUN, iv] < aj[p, SUN, iv]
                    psn_wp_z_sun[p, iv] = psn_z_sun[p, iv]
                end

                psn_z_sha[p, iv] = ag[p, SHA, iv] * o3coefv_sha[p]
                psn_wc_z_sha[p, iv] = zero(T)
                psn_wj_z_sha[p, iv] = zero(T)
                psn_wp_z_sha[p, iv] = zero(T)
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
end

# Copy every array field of a host bundle back into the matching field of the
# device bundle (used by the PHS Pass-3 hybrid below). Read-only fields copy back
# unchanged (harmless); output fields carry the host-computed Pass-3 result onto
# the device. Bulk copyto! (no scalar indexing) so it is device-safe.
function _phs_pass3_copyback!(dev, host)
    for f in fieldnames(typeof(dev))
        d = getfield(dev, f); h = getfield(host, f)
        (d isa AbstractArray && h isa AbstractArray) && copyto!(d, h)
    end
    return nothing
end

# Host launcher for the PHS Pass-3 kernel. GPU-hostile values (params_inst fields,
# the Medlyn override, the is_near_local_noon function, module-global constants) are
# resolved to host scalars / a per-patch Bool vector before the launch.
function psn_phs_pass3_update!(ps, mask_patch, col_of_patch, ivt, nrad,
        medlynslope_pft, medlynintercept_pft, crop_pft,
        forc_pbot, tgcm, rb, eair, esat_tv, cair, oair, qsatl, qaf,
        laisun, laisha, htop, tsai, elai, esai, fdry, forc_rho,
        par_z_sun_in, par_z_sha_in, jmax_z_local,
        o3coefg_sun, o3coefg_sha, o3coefv_sun, o3coefv_sha,
        k_soil_root, smp_l, z_col,
        vegwp, vegwp_ln, qflx_tran_veg, bsun_arr, bsha_arr,
        rh_leaf_sun, rh_leaf_sha,
        psn_wc_z_sun, psn_wj_z_sun, psn_wp_z_sun,
        psn_wc_z_sha, psn_wj_z_sha, psn_wp_z_sha,
        nlevsoi::Int, modifyphoto_and_lmr_forcrop::Bool, stomatalcond_mtd::Int,
        rsmax0, overrides::CalibrationOverrides, is_near_local_noon_fn::Function,
        bounds_patch)
    T = eltype(forc_pbot)
    sc = PsnPhsP3Scalars{T}(T(params_inst.fnps), T(params_inst.theta_psii),
        T(params_inst.theta_ip), T(overrides.medlyn_slope),
        T(RGAS), T(MAX_CS), T(MEDLYN_RH_CAN_MAX), T(MEDLYN_RH_CAN_FACT),
        T(rsmax0), T(SPVAL))
    # is_near_local_noon resolved to a per-patch Bool vector on the host (the
    # closure is not GPU-callable). On the device path this is copied onto the
    # backend so the kernel reads a device Bool array.
    nln = [is_near_local_noon_fn(p) for p in eachindex(mask_patch)]
    is_near_local_noon = similar(mask_patch, Bool, length(nln))
    copyto!(is_near_local_noon, nln)
    phs_params = _psn_phs_params(forc_pbot)
    dv = _psn_dv(ps)
    # group loose @Const inputs + mutable outputs into device-view bundles (Metal
    # 31-arg-buffer reduction); same refs, so writes flow back and behaviour is
    # unchanged.
    idx = _Psn3Idx(; mask_patch = mask_patch, is_near_local_noon = is_near_local_noon,
        col_of_patch = col_of_patch, ivt = ivt, nrad = nrad,
        medlynslope_pft = medlynslope_pft, medlynintercept_pft = medlynintercept_pft,
        crop_pft = crop_pft)
    ina = _Psn3InA(; forc_pbot = forc_pbot, tgcm = tgcm, rb = rb, eair = eair,
        esat_tv = esat_tv, cair = cair, oair = oair, qsatl = qsatl, qaf = qaf,
        laisun = laisun, laisha = laisha, htop = htop, tsai = tsai, elai = elai,
        esai = esai, fdry = fdry, forc_rho = forc_rho)
    inb = _Psn3InB(; o3coefg_sun = o3coefg_sun, o3coefg_sha = o3coefg_sha,
        o3coefv_sun = o3coefv_sun, o3coefv_sha = o3coefv_sha,
        par_z_sun_in = par_z_sun_in, par_z_sha_in = par_z_sha_in,
        jmax_z_local = jmax_z_local, k_soil_root = k_soil_root, smp_l = smp_l,
        z_col = z_col)
    out = _Psn3Out(; qflx_tran_veg = qflx_tran_veg, bsun_arr = bsun_arr,
        bsha_arr = bsha_arr, rh_leaf_sun = rh_leaf_sun, rh_leaf_sha = rh_leaf_sha,
        vegwp = vegwp, vegwp_ln = vegwp_ln,
        psn_wc_z_sun = psn_wc_z_sun, psn_wj_z_sun = psn_wj_z_sun,
        psn_wp_z_sun = psn_wp_z_sun, psn_wc_z_sha = psn_wc_z_sha,
        psn_wj_z_sha = psn_wj_z_sha, psn_wp_z_sha = psn_wp_z_sha)
    be = _kernel_backend(forc_pbot)
    if _PSN_CI_AD_HOSTLOOP[]
        # AD-mode host loop (see compositional_reverse! / _photosynth_ci_body!): runs the
        # SAME fused per-patch Newton body the KA kernel runs, but as a plain loop so
        # Enzyme reverse on Julia 1.12 dodges the KA kernel-launch codegen. Byte-identical
        # primal (shared _psn_phs_pass3_body!).
        @inbounds for p in 1:length(bounds_patch)
            _psn_phs_pass3_body!(p, dv, phs_params, sc, idx, ina, inb, out,
                nlevsoi, modifyphoto_and_lmr_forcrop, stomatalcond_mtd,
                STOMATALCOND_MTD_BB1987, STOMATALCOND_MTD_MEDLYN2011)
        end
    elseif be isa KA.CPU
        _psn_phs_pass3_kernel!(be)(dv, phs_params, sc, idx, ina, inb, out,
            nlevsoi, modifyphoto_and_lmr_forcrop,
            stomatalcond_mtd, STOMATALCOND_MTD_BB1987, STOMATALCOND_MTD_MEDLYN2011;
            ndrange = length(bounds_patch))
        KA.synchronize(be)
    else
        # HYBRID FALLBACK: the fused PHS Pass-3 Newton kernel is CPU-bit-identical
        # but exceeds the Apple Metal shader compiler's capacity (the compiler
        # daemon crashes, XPC_ERROR_CONNECTION_INTERRUPTED — a compiler limit, not
        # a code bug; it would likely compile on CUDA). Run the SAME kernel on the
        # CPU over host copies of the device bundles, then copy the mutated arrays
        # back to the device. Passes 1/2/4 still run on the GPU; only this one
        # coupled Newton solve hops to the host. Float32 throughout, so the device
        # state after copyback matches a would-be on-device F32 result to last bits.
        dv_h  = Adapt.adapt(Array, dv)
        pp_h  = Adapt.adapt(Array, phs_params)
        idx_h = Adapt.adapt(Array, idx)
        ina_h = Adapt.adapt(Array, ina)
        inb_h = Adapt.adapt(Array, inb)
        out_h = Adapt.adapt(Array, out)
        _psn_phs_pass3_kernel!(KA.CPU())(dv_h, pp_h, sc, idx_h, ina_h, inb_h, out_h,
            nlevsoi, modifyphoto_and_lmr_forcrop,
            stomatalcond_mtd, STOMATALCOND_MTD_BB1987, STOMATALCOND_MTD_MEDLYN2011;
            ndrange = length(bounds_patch))
        KA.synchronize(KA.CPU())
        _phs_pass3_copyback!(dv, dv_h)
        _phs_pass3_copyback!(out, out_h)
    end
    # PARITY (Phase 3c): opt-in read-only dump of converged sunlit-leaf photosynthesis
    # internals, mirroring the Fortran FPHOTO/FVCM25 SourceMods lines. Default OFF →
    # byte-identical. Appends one block per pass-3 call (i.e. per canopy leaf-temp
    # Newton iteration) so the Julia trajectory lines up with the Fortran trajectory.
    if PHS_PHOTO_DEBUG[]
        _phs_photo_debug_dump(ps, bsun_arr, rh_leaf_sun, eair, esat_tv, cair, oair,
            par_z_sun_in, jmax_z_local, nrad, mask_patch)
    end
    return nothing
end

# Read-only dump helper for the PHS photosynthesis-internals parity localization.
# Reads only already-converged ps arrays + pass-3 inputs; recomputes je_sun the same
# way the body does. Writes nothing into model state.
function _phs_photo_debug_dump(ps, bsun_arr, rh_leaf_sun, eair, esat_tv, cair, oair,
        par_z_sun_in, jmax_z_local, nrad, mask_patch)
    fnps = params_inst.fnps; theta_psii = params_inst.theta_psii
    row(p, iv, key, vals...) = string("JPHOTO p=", p, " iv=", iv, " ", key, "= ",
                                      join(vals, " "))
    open(PHS_PHOTO_DEBUG_PATH[], "a") do io
        for p in PHS_PHOTO_DEBUG_PATCHES[]
            (1 <= p <= length(mask_patch) && mask_patch[p]) || continue
            for iv in 1:nrad[p]
                vc = ps.vcmax_z_phs_patch[p, SUN, iv]
                jm = jmax_z_local[p, SUN, iv]
                tp = ps.tpu_z_phs_patch[p, SUN, iv]
                lmr = ps.lmrsun_z_patch[p, iv]
                # je_sun (recomputed identically to the body)
                qabs = 0.5 * (1.0 - fnps) * par_z_sun_in[p, iv] * 4.6
                r1, r2 = quadratic_solve(theta_psii, -(qabs + jm), qabs * jm)
                je = min(r1, r2)
                println(io, row(p, iv, "vcmax_jmax_tpu_lmr_sun", vc, jm, tp, lmr))
                println(io, row(p, iv, "cp_kc_ko_jesun",
                               ps.cp_patch[p], ps.kc_patch[p], ps.ko_patch[p], je))
                println(io, row(p, iv, "ci_ac_aj_ap_ag_sun",
                               ps.cisun_z_patch[p, iv], ps.ac_phs_patch[p, SUN, iv],
                               ps.aj_phs_patch[p, SUN, iv], ps.ap_phs_patch[p, SUN, iv],
                               ps.ag_phs_patch[p, SUN, iv]))
                println(io, row(p, iv, "an_gsmol_bsun_rs_sun",
                               ps.an_sun_patch[p, iv], ps.gs_mol_sun_patch[p, iv],
                               bsun_arr[p], ps.rssun_z_patch[p, iv]))
                println(io, row(p, iv, "rhcan_rhleaf_cair_oair_esattv_eair",
                               ps.vpd_can_patch[p], rh_leaf_sun[p], cair[p], oair[p],
                               esat_tv[p], eair[p]))
                println(io, row(p, iv, "gbmol_par_z_sun",
                               ps.gb_mol_patch[p], par_z_sun_in[p, iv]))
                println(io, string("JVCM25 p=", p, " iv=", iv,
                               " luna_vcmx25z_jmx25z_vc25top_jm25top= ",
                               join((ps.vcmx25_z_patch[p, iv], ps.jmx25_z_patch[p, iv],
                                     ps.luvcmax25top_patch[p], ps.lujmax25top_patch[p]), " ")))
            end
        end
    end
    return nothing
end

# =====================================================================
# photosynthesis_hydrstress! — Main PHS photosynthesis
# =====================================================================

"""
    photosynthesis_hydrstress!(ps, ...)

Leaf photosynthesis and stomatal conductance with plant hydraulic stress.
Sunlit and shaded photosynthesis and stomatal conductance are solved
simultaneously per Pierre Gentine/Daniel Kennedy PHS method.
Ported from `subroutine PhotosynthesisHydraulicStress`.
"""
function photosynthesis_hydrstress!(ps,
                                    esat_tv::AbstractVector{<:Real},
                                    eair::AbstractVector{<:Real},
                                    oair::AbstractVector{<:Real},
                                    cair::AbstractVector{<:Real},
                                    rb::AbstractVector{<:Real},
                                    btran::AbstractVector{<:Real},
                                    dayl_factor::AbstractVector{<:Real},
                                    leafn::AbstractVector{<:Real},
                                    qsatl::AbstractVector{<:Real},
                                    qaf::AbstractVector{<:Real},
                                    forc_pbot::AbstractVector{<:Real},
                                    forc_rho::AbstractVector{<:Real},
                                    t_veg::AbstractVector{<:Real},
                                    t10::AbstractVector{<:Real},
                                    tgcm::AbstractVector{<:Real},
                                    nrad::AbstractVector{<:Integer},
                                    tlai_z::AbstractMatrix{<:Real},
                                    tlai::AbstractVector{<:Real},
                                    tsai::AbstractVector{<:Real},
                                    par_z_sun_in::AbstractMatrix{<:Real},
                                    par_z_sha_in::AbstractMatrix{<:Real},
                                    lai_z_sun_in::AbstractMatrix{<:Real},
                                    lai_z_sha_in::AbstractMatrix{<:Real},
                                    vcmaxcint_sun::AbstractVector{<:Real},
                                    vcmaxcint_sha::AbstractVector{<:Real},
                                    o3coefv_sun::AbstractVector{<:Real},
                                    o3coefg_sun::AbstractVector{<:Real},
                                    o3coefv_sha::AbstractVector{<:Real},
                                    o3coefg_sha::AbstractVector{<:Real},
                                    c3psn_pft::AbstractVector{<:Real},
                                    leafcn_pft::AbstractVector{<:Real},
                                    flnr_pft::AbstractVector{<:Real},
                                    fnitr_pft::AbstractVector{<:Real},
                                    slatop_pft::AbstractVector{<:Real},
                                    mbbopt_pft::AbstractVector{<:Real},
                                    medlynintercept_pft::AbstractVector{<:Real},
                                    medlynslope_pft::AbstractVector{<:Real},
                                    froot_leaf_pft::AbstractVector{<:Real},
                                    root_radius_pft::AbstractVector{<:Real},
                                    root_density_pft::AbstractVector{<:Real},
                                    crop_pft::AbstractVector{<:Real},
                                    ivt::AbstractVector{<:Integer},
                                    col_of_patch::AbstractVector{<:Integer},
                                    mask_patch::AbstractVector{Bool},
                                    bounds_patch::UnitRange{Int},
                                    froot_carbon::AbstractVector{<:Real},
                                    croot_carbon::AbstractVector{<:Real},
                                    k_soil_root::AbstractMatrix{<:Real},
                                    root_conductance_out::AbstractMatrix{<:Real},
                                    soil_conductance_out::AbstractMatrix{<:Real},
                                    rootfr::AbstractMatrix{<:Real},
                                    dz::AbstractMatrix{<:Real},
                                    z_col::AbstractMatrix{<:Real},
                                    hk_l::AbstractMatrix{<:Real},
                                    hksat::AbstractMatrix{<:Real},
                                    smp_l::AbstractMatrix{<:Real},
                                    vegwp::AbstractMatrix{<:Real},
                                    vegwp_ln::AbstractMatrix{<:Real},
                                    laisun::AbstractVector{<:Real},
                                    laisha::AbstractVector{<:Real},
                                    elai::AbstractVector{<:Real},
                                    esai::AbstractVector{<:Real},
                                    htop::AbstractVector{<:Real},
                                    fdry::AbstractVector{<:Real},
                                    qflx_tran_veg::AbstractVector{<:Real};
                                    nlevcan::Int=NLEVCAN,
                                    nlevsoi::Int=varpar.nlevsoi,
                                    use_cn::Bool=false,
                                    use_luna::Bool=false,
                                    use_c13::Bool=false,
                                    leaf_mr_vcm::Real=0.015,
                                    overrides::CalibrationOverrides=CalibrationOverrides(),
                                    is_near_local_noon_fn::Function=(p) -> false,
                                    bsun_patch=nothing,
                                    bsha_patch=nothing)
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
    FT = eltype(t_veg)
    # Device-resident scratch (fill!(similar(t_veg,…)) so the kernelized passes
    # write into backend arrays; on CPU this is a plain zeros/ones-filled Array).
    # bbbopt is now a kernel-LOCAL scalar in PHS Pass 2 (was a per-patch array).
    jmax_z_local = fill!(similar(t_veg, FT, np, 2, nlevcan), zero(FT))
    kn           = fill!(similar(t_veg, FT, np), zero(FT))
    psn_wc_z_sun = fill!(similar(t_veg, FT, np, nlevcan), zero(FT))
    psn_wj_z_sun = fill!(similar(t_veg, FT, np, nlevcan), zero(FT))
    psn_wp_z_sun = fill!(similar(t_veg, FT, np, nlevcan), zero(FT))
    psn_wc_z_sha = fill!(similar(t_veg, FT, np, nlevcan), zero(FT))
    psn_wj_z_sha = fill!(similar(t_veg, FT, np, nlevcan), zero(FT))
    psn_wp_z_sha = fill!(similar(t_veg, FT, np, nlevcan), zero(FT))
    rh_leaf_sun  = fill!(similar(t_veg, FT, np), zero(FT))   # written by the Pass 3 kernel
    rh_leaf_sha  = fill!(similar(t_veg, FT, np), zero(FT))
    # Use the caller's energyflux bsun/bsha_patch as the working array when supplied
    # (so the BSUN/BSHA history diagnostics get populated — they are NaN-initialized and
    # otherwise never written); else a local scratch. Init to 1.0 to match the scratch.
    bsun_arr = (bsun_patch === nothing) ? fill!(similar(t_veg, FT, np), one(FT)) :
               fill!(bsun_patch, one(eltype(bsun_patch)))
    bsha_arr = (bsha_patch === nothing) ? fill!(similar(t_veg, FT, np), one(FT)) :
               fill!(bsha_patch, one(eltype(bsha_patch)))

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

    # ---- Pass 1: Root-soil conductance (kernelized) ----
    # One thread per patch; the per-soil-level j loop runs serially inside each
    # thread (own-index k_soil_root[p,j] writes, no race). _plc(smp,psi50,ck) +
    # krmax come from the device-resident phs_params bundle.
    psn_phs_pass1_update!(ps, mask_patch, col_of_patch, ivt, froot_carbon,
        rootfr, dz, tsai, tlai, froot_leaf_pft, root_radius_pft, root_density_pft,
        hksat, hk_l, smp_l, z_col, k_soil_root, root_conductance_out,
        soil_conductance_out, c_to_b, croot_lateral_length, nlevsoi, bounds_patch)

    # ---- Pass 2: Kinetics, N profile, vcmax, respiration (kernelized) ----
    # One thread per patch; internal serial iv loop writes both the 3D ps.*_phs
    # fields ([p,SUN,iv] & [p,SHA,iv]) and per-phase lmr in one pass. bbbopt is a
    # kernel-local scalar; overrides/leaf_mr_vcm/jmax25top_sf resolved on the host.
    psn_phs_pass2_update!(ps, kn, jmax_z_local, mask_patch, ivt, c3psn_pft,
        mbbopt_pft, forc_pbot, oair, slatop_pft, leafcn_pft, flnr_pft, fnitr_pft,
        dayl_factor, t10, t_veg, tlai_z, par_z_sun_in, par_z_sha_in, vcmaxcint_sun,
        vcmaxcint_sha, nrad, use_cn, use_c13, leaf_mr_vcm, nlevcan,
        stomatalcond_mtd, light_inhibit, leafresp_method, overrides, bounds_patch;
        use_luna=use_luna)

    # ---- Pass 3: Leaf-level photosynthesis (kernelized) ----
    # One thread per patch; the entire per-patch PHS Newton solve (hybrid_PHS →
    # ci_func_PHS → brent_PHS/calcstress → spacF/spacA/getvegwp/getqflx) runs
    # in-thread via the positional device cores. Reads Pass 1/2 outputs
    # (k_soil_root, vcmax_z_phs, lmr_z, …) and writes what Pass 4 reads.
    psn_phs_pass3_update!(ps, mask_patch, col_of_patch, ivt, nrad,
        medlynslope_pft, medlynintercept_pft, crop_pft,
        forc_pbot, tgcm, rb, eair, esat_tv, cair, oair, qsatl, qaf,
        laisun, laisha, htop, tsai, elai, esai, fdry, forc_rho,
        par_z_sun_in, par_z_sha_in, jmax_z_local,
        o3coefg_sun, o3coefg_sha, o3coefv_sun, o3coefv_sha,
        k_soil_root, smp_l, z_col,
        vegwp, vegwp_ln, qflx_tran_veg, bsun_arr, bsha_arr,
        rh_leaf_sun, rh_leaf_sha,
        psn_wc_z_sun, psn_wj_z_sun, psn_wp_z_sun,
        psn_wc_z_sha, psn_wj_z_sha, psn_wp_z_sha,
        nlevsoi, modifyphoto_and_lmr_forcrop, stomatalcond_mtd,
        rsmax0, overrides, is_near_local_noon_fn, bounds_patch)

    # ---- Pass 4: Canopy integration (kernelized) ----
    # One thread per patch; internal serial iv loops reduce sun+sha per-canopy
    # fluxes → per-patch psn/rs/lmr, then blend btran from LAI-weighted bsun/bsha.
    psn_phs_pass4_update!(ps, mask_patch, nrad, ivt, crop_pft,
        lai_z_sun_in, lai_z_sha_in, rb, btran,
        psn_wc_z_sun, psn_wj_z_sun, psn_wp_z_sun,
        psn_wc_z_sha, psn_wj_z_sha, psn_wp_z_sha,
        bsun_arr, bsha_arr, modifyphoto_and_lmr_forcrop, bounds_patch)

    return nothing
end
