# TwoStreamMLPEMod.jl
# Julia port of FATES src/fates/radiation/TwoStreamMLPEMod.F90
#
# Two-stream radiation scattering of vegetation in "M"ultiple "L"ayers with
# "P"arallel "E"lements (MLPE).
#
# There may be numerous canopy layers. In each canopy layer, plant media for
# different functional types are grouped so that they inhabit their own
# exclusive footprint ("column") within the layer. Within each column there are
# further sub-layer discretizations organized by top-down integrated vegetation
# area index.
#
# There is a separate allocation and call sequence for each broad band; the
# two_stream_type is instantiated for each broad band.
#
# Assumptions: band index 1 = visible (vis)
#                         2 = near infrared (nir)
#                         3 = thermal (not used at the moment)
#
# Deps: FatesConstantsMod (fates_r8), FatesGlobals (fates_log/fates_endrun).
# Standalone solver — NOT added to CLMInstances or any dual-copied struct.

using LinearAlgebra: lu!, ldiv!

# ---------------------------------------------------------------------------
# Module-level parameter constants
# ---------------------------------------------------------------------------
const twostr_nearzero = 1.0e-20            # local nearzero (matches Fortran nearzero in this module)
const twostr_debug = true
const twostr_unset_r8 = 1.0e-36
const twostr_unset_int = -999
const twostr_vis = 1                       # named index of visible shortwave radiation
const twostr_nir = 2                       # named index for near infrared shortwave radiation

# Allowable error, as a fraction of total incident, for canopy radiation balance
const rel_err_thresh = 1.0e-6
const area_err_thresh = rel_err_thresh * 0.1

# Codes for how the upper boundary is specified
const normalized_upper_boundary = 1
const absolute_upper_boundary = 2

# Snow optical parameter constants for visible (index=1) and NIR (index=2)
const betad_snow = (0.5, 0.5)   # diffuse backscatter fraction    (CLM50 Tech Man)
const betab_snow = (0.5, 0.5)   # beam backscatter fraction       (CLM50 Tech Man)
const om_snow = (0.8, 0.4)      # scattering coefficient for snow (CLM50 Tech Man)

# Cap the maximum optical depth to keep exponents from blowing up
const kb_max = 30.0

# For air, use nominal values to prevent div0s; the key is that vai = 0
const k_air = 0.5
const om_air = 0.5
const beta_air = 0.5
const air_ft = 0

# ---------------------------------------------------------------------------
# rad_params_type — parameter constants, specific to plant material and band
# ---------------------------------------------------------------------------
Base.@kwdef mutable struct rad_params_type
    # From the parameter file
    rhol::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)            # leaf reflectance   (band x pft)
    rhos::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)            # stem reflectance   (band x pft)
    taul::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)            # leaf transmittance (band x pft)
    taus::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)            # stem transmittance (band x pft)
    xl::Vector{Float64} = Float64[]                                  # leaf/stem orientation (pft)
    clumping_index::Vector{Float64} = Float64[]                      # clumping index 0-1 (pft)

    # Derived parameters
    phi1::Vector{Float64} = Float64[]        # intermediate term for kd and kb
    phi2::Vector{Float64} = Float64[]        # intermediate term for kd and kb
    avmu::Vector{Float64} = Float64[]        # average inverse optical depth per unit leaf+stem area
    kd_leaf::Vector{Float64} = Float64[]     # mean optical depth per unit area leaves in diffuse
    kd_stem::Vector{Float64} = Float64[]     # mean optical depth per unit area stems in diffuse
    om_leaf::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)         # leaf scattering coeff (band x pft)
    om_stem::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)         # stem scattering coeff (band x pft)
end

# Module-level instance, mirroring Fortran `type(rad_params_type),public :: rad_params`
const rad_params = rad_params_type()

# ---------------------------------------------------------------------------
# scelg_type — scattering elements, geometry, wavelength-independent
# ---------------------------------------------------------------------------
Base.@kwdef mutable struct scelg_type
    pft::Int = twostr_unset_int    # pft index
    area::Float64 = twostr_unset_r8 # m2 col / m2 ground
    lai::Float64 = twostr_unset_r8  # m2 leaf area / m2 col
    sai::Float64 = twostr_unset_r8  # m2 stem area / m2 col
    Kb::Float64 = twostr_unset_r8       # optical depth of beam radiation
    Kb_leaf::Float64 = twostr_unset_r8  # optical depth of just leaves in beam radiation
    Kd::Float64 = twostr_unset_r8       # optical depth of diffuse radiation
    area_squeeze::Float64 = 1.0    # ratio of element area to area of its constituents
end

# ---------------------------------------------------------------------------
# scelb_type — scattering elements, wavelength (band) dependent
# ---------------------------------------------------------------------------
Base.@kwdef mutable struct scelb_type
    # Terms used in the final solution, also used for decomposing the solution
    Au::Float64 = twostr_unset_r8           # compound intercept term
    Ad::Float64 = twostr_unset_r8           # compound intercept term
    B1::Float64 = twostr_unset_r8           # compound term w/ lambdas (operates on e^{av})
    B2::Float64 = twostr_unset_r8           # compound term w/ lambdas (operates on e^{-av})
    lambda1_diff::Float64 = twostr_unset_r8 # compound term w/ B for diffuse forcing
    lambda2_diff::Float64 = twostr_unset_r8 # compound term w/ B for diffuse forcing
    lambda1_beam::Float64 = twostr_unset_r8 # compound term w/ B for beam forcing
    lambda2_beam::Float64 = twostr_unset_r8 # compound term w/ B for beam forcing

    a::Float64 = twostr_unset_r8      # complex term operating on veg area index
    om::Float64 = twostr_unset_r8     # scattering coefficient for media as a whole
    betad::Float64 = twostr_unset_r8  # backscatter fraction of diffuse radiation for media
    betab::Float64 = twostr_unset_r8  # backscatter fraction of beam radiation for media
    Rbeam0::Float64 = twostr_unset_r8 # normalized downwelling beam at top of element [-]
end

# ---------------------------------------------------------------------------
# band_type — per-band scattering coefficients
# ---------------------------------------------------------------------------
Base.@kwdef mutable struct band_type
    scelb::Matrix{scelb_type} = Matrix{scelb_type}(undef, 0, 0) # (layer, column), can be sparse
    ib::Int = twostr_unset_int          # band index, should be consistent with rad_params
    Rbeam_atm::Float64 = twostr_unset_r8        # downwelling beam from atmosphere [W/m2 ground]
    Rdiff_atm::Float64 = twostr_unset_r8        # downwelling diffuse from atmosphere [W/m2 ground]
    albedo_grnd_diff::Float64 = twostr_unset_r8 # ground albedo diffuse
    albedo_grnd_beam::Float64 = twostr_unset_r8 # ground albedo direct
end

# ---------------------------------------------------------------------------
# twostream_type — parent type holding almost everything in the solver
# ---------------------------------------------------------------------------
Base.@kwdef mutable struct twostream_type
    scelg::Matrix{scelg_type} = Matrix{scelg_type}(undef, 0, 0) # (layer, column), wavelength-indep info
    band::Vector{band_type} = band_type[]                        # per band (vis, nir, ...)

    n_bands::Int = 0                # number of bands
    n_lyr::Int = 0                  # number of (vertical) scattering element layers
    n_col::Vector{Int} = Int[]      # number of (horizontal) columns per layer
    n_scel::Int = 0                 # total number of scattering elements
    force_prep::Bool = false        # signal to update diffuse scattering coefficients
    frac_snow::Float64 = twostr_unset_r8      # current mean snow fraction of the canopy
    frac_snow_old::Float64 = twostr_unset_r8  # previous mean snow fraction of the canopy
    cosz::Float64 = twostr_unset_r8           # current cosine of the zenith angle
end

# ---------------------------------------------------------------------------
# TwoStreamLogInit / endrun
# ---------------------------------------------------------------------------
function TwoStreamLogInit(log_unit_in::Integer)
    # log unit is provided via FatesGlobals.fates_log(); kept for API parity.
    return nothing
end

# Module-local endrun mirrors the Fortran subroutine; routes to fates_endrun.
twostream_endrun(msg::AbstractString) = fates_endrun(msg)

# ===============================================================================================

"""
    AllocInitTwoStream!(this, band_indices, ncan, ncol)

Allocate and initialize a `twostream_type` for `ncan` canopy layers, `ncol`
columns, and the given band indices.
"""
function AllocInitTwoStream!(this::twostream_type, band_indices::AbstractVector{<:Integer},
                             ncan::Integer, ncol::Integer)
    nbands = length(band_indices)

    this.n_col = fill(twostr_unset_int, ncan)
    this.scelg = [scelg_type() for _ in 1:ncan, _ in 1:ncol]
    this.band = Vector{band_type}(undef, nbands)

    this.n_bands = nbands
    this.n_lyr = ncan
    this.frac_snow = twostr_unset_r8
    this.frac_snow_old = twostr_unset_r8

    for ib in 1:nbands
        b = band_type()
        b.scelb = [scelb_type() for _ in 1:ncan, _ in 1:ncol]
        b.albedo_grnd_diff = twostr_unset_r8
        b.albedo_grnd_beam = twostr_unset_r8
        b.ib = Int(band_indices[ib])
        this.band[ib] = b
    end

    return nothing
end

# ===============================================================================================

"""
    DeallocTwoStream!(this)

Reset the allocatable arrays of a `twostream_type` (Julia GC handles memory;
this clears the references for API parity).
"""
function DeallocTwoStream!(this::twostream_type)
    this.scelg = Matrix{scelg_type}(undef, 0, 0)
    this.n_col = Int[]
    this.band = band_type[]
    return nothing
end

# ===============================================================================================

"""
    AllocateRadParams(n_pft, n_bands)

Allocate the module-level `rad_params` arrays for the given number of PFTs and
broad bands.
"""
function AllocateRadParams(n_pft::Integer, n_bands::Integer)
    rad_params.rhol = Matrix{Float64}(undef, n_bands, n_pft)
    rad_params.rhos = Matrix{Float64}(undef, n_bands, n_pft)
    rad_params.taul = Matrix{Float64}(undef, n_bands, n_pft)
    rad_params.taus = Matrix{Float64}(undef, n_bands, n_pft)
    rad_params.xl = Vector{Float64}(undef, n_pft)
    rad_params.clumping_index = Vector{Float64}(undef, n_pft)

    rad_params.phi1 = Vector{Float64}(undef, n_pft)
    rad_params.phi2 = Vector{Float64}(undef, n_pft)
    rad_params.avmu = Vector{Float64}(undef, n_pft)
    rad_params.kd_leaf = Vector{Float64}(undef, n_pft)
    rad_params.kd_stem = Vector{Float64}(undef, n_pft)
    rad_params.om_leaf = Matrix{Float64}(undef, n_bands, n_pft)
    rad_params.om_stem = Matrix{Float64}(undef, n_bands, n_pft)
    return nothing
end

# ================================================================================================
# Radiation-intensity query routines
# ================================================================================================

"""
    GetRdDn(this, ican, icol, ib, vai) -> r_diff_dn

Downwelling diffuse radiation at integrated veg area index `vai` within element
(ican, icol) for band `ib`.

    Rdn = Ad e^(-Kb v) + λ1 B2 e^(av) + λ2 B1 e^(-av)
"""
function GetRdDn(this::twostream_type, ican::Integer, icol::Integer, ib::Integer, vai::Real)
    scelb = this.band[ib].scelb[ican, icol]
    scelg = this.scelg[ican, icol]

    r_diff_dn = this.band[ib].Rbeam_atm * (
        scelb.Ad * exp(-scelg.Kb * vai) +
        scelb.B2 * scelb.lambda1_beam * exp(scelb.a * vai) +
        scelb.B1 * scelb.lambda2_beam * exp(-scelb.a * vai)) +
        this.band[ib].Rdiff_atm * (
        scelb.B2 * scelb.lambda1_diff * exp(scelb.a * vai) +
        scelb.B1 * scelb.lambda2_diff * exp(-scelb.a * vai))

    if twostr_debug
        if isnan(r_diff_dn)
            twostream_endrun("GETRDN: NaN downwelling diffuse radiation in GetRdDn")
        end
    end

    return r_diff_dn
end

"""
    GetRdUp(this, ican, icol, ib, vai) -> r_diff_up

Upwelling diffuse radiation at integrated veg area index `vai` within element
(ican, icol) for band `ib`.

    Rup = Au e^(-Kb v) + λ1 B1 e^(av) + λ2 B2 e^(-av)
"""
function GetRdUp(this::twostream_type, ican::Integer, icol::Integer, ib::Integer, vai::Real)
    scelb = this.band[ib].scelb[ican, icol]
    scelg = this.scelg[ican, icol]

    r_diff_up = this.band[ib].Rbeam_atm * (
        scelb.Au * exp(-scelg.Kb * vai) +
        scelb.B1 * scelb.lambda1_beam * exp(scelb.a * vai) +
        scelb.B2 * scelb.lambda2_beam * exp(-scelb.a * vai)) +
        this.band[ib].Rdiff_atm * (
        scelb.B1 * scelb.lambda1_diff * exp(scelb.a * vai) +
        scelb.B2 * scelb.lambda2_diff * exp(-scelb.a * vai))

    return r_diff_up
end

"""
    GetRb(this, ican, icol, ib, vai) -> r_beam_dn

Downwelling direct-beam radiation at integrated veg area index `vai`.
"""
function GetRb(this::twostream_type, ican::Integer, icol::Integer, ib::Integer, vai::Real)
    return this.band[ib].Rbeam_atm *
           this.band[ib].scelb[ican, icol].Rbeam0 * exp(-this.scelg[ican, icol].Kb * vai)
end

"""
    GetAbsRad(this, ican, icol, ib, vai_top, vai_bot)
        -> (Rb_abs, Rd_abs, Rd_abs_leaf, Rb_abs_leaf, R_abs_stem, R_abs_snow, leaf_sun_frac)

Decompose radiation scattering and return absorbed radiation over the interval
[vai_top, vai_bot] within element (ican, icol). VAI is zero at the top of the
element and increases downwards.
"""
function GetAbsRad(this::twostream_type, ican::Integer, icol::Integer, ib::Integer,
                   vai_top::Real, vai_bot::Real)
    scelb = this.band[ib].scelb[ican, icol]
    scelg = this.scelg[ican, icol]
    ft = scelg.pft

    # If this is air, trivial solutions
    if ft == air_ft
        return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    end

    # The total vegetation area index of the element
    vai_max = scelg.lai + scelg.sai

    dvai = vai_bot - vai_top

    lai_top = vai_top * scelg.lai / (scelg.lai + scelg.sai)
    lai_bot = vai_bot * scelg.lai / (scelg.lai + scelg.sai)
    dlai = dvai * scelg.lai / (scelg.lai + scelg.sai)

    local leaf_sun_frac::Float64
    if dlai > twostr_nearzero
        leaf_sun_frac = max(0.001, min(0.999,
            scelb.Rbeam0 / (dlai * scelg.Kb_leaf / rad_params.clumping_index[ft]) *
            (exp(-scelg.Kb_leaf * lai_top) - exp(-scelg.Kb_leaf * lai_bot))))
    else
        leaf_sun_frac = 0.001
    end

    if twostr_debug
        if leaf_sun_frac > 1.0 || leaf_sun_frac < 0.0
            twostream_endrun("impossible leaf sun fraction in GetAbsRad")
        end
    end

    # Disentangle absorption between leaves and stems; weight by area*K*(1-om)
    frac_abs_snow = this.frac_snow * (1.0 - om_snow[ib]) / (1.0 - scelb.om)

    diff_wt_leaf = scelg.lai * (1.0 - rad_params.om_leaf[ib, ft]) * rad_params.kd_leaf[ft]
    diff_wt_stem = scelg.sai * (1.0 - rad_params.om_stem[ib, ft]) * rad_params.kd_stem[ft]

    beam_wt_leaf = scelg.lai * (1.0 - rad_params.om_leaf[ib, ft]) * scelg.Kb_leaf
    beam_wt_stem = scelg.sai * (1.0 - rad_params.om_stem[ib, ft]) * 1.0

    if twostr_debug
        if (vai_bot - vai_max) > rel_err_thresh
            twostream_endrun("GetAbsRad: requested VAI larger than total integrated VAI of element")
        end
        if (vai_bot - vai_top) < -rel_err_thresh
            twostream_endrun("GetAbsRad: lower-position VAI greater than upper-position VAI")
        end
    end

    # Absorbed radiation = energy balance (net) over the depth of interest
    Rb_net = GetRb(this, ican, icol, ib, vai_top) - GetRb(this, ican, icol, ib, vai_bot)

    Rd_net = (GetRdDn(this, ican, icol, ib, vai_top) - GetRdDn(this, ican, icol, ib, vai_bot)) +
             (GetRdUp(this, ican, icol, ib, vai_bot) - GetRdUp(this, ican, icol, ib, vai_top))

    # Net beam includes absorbed + re-scattered; re-scattered sources the diffuse balance
    Rb_abs = Rb_net * (1.0 - scelb.om)
    Rd_abs = Rd_net + Rb_net * scelb.om

    Rb_abs_leaf = (1.0 - frac_abs_snow) * Rb_abs * beam_wt_leaf / (beam_wt_leaf + beam_wt_stem)
    Rd_abs_leaf = (1.0 - frac_abs_snow) * Rd_abs * diff_wt_leaf / (diff_wt_leaf + diff_wt_stem)

    R_abs_snow = (Rb_abs + Rd_abs) * frac_abs_snow

    R_abs_stem = (1.0 - frac_abs_snow) *
                 (Rb_abs * beam_wt_stem / (beam_wt_leaf + beam_wt_stem) +
                  Rd_abs * diff_wt_stem / (diff_wt_leaf + diff_wt_stem))

    return (Rb_abs, Rd_abs, Rd_abs_leaf, Rb_abs_leaf, R_abs_stem, R_abs_snow, leaf_sun_frac)
end

# ================================================================================================

"""
    RadParamPrep()

Derive the band-independent and band-dependent rad_params terms from the input
optical properties. Called once during initialization.
"""
function RadParamPrep()
    numpft = size(rad_params.om_leaf, 2)
    nbands = size(rad_params.om_leaf, 1)

    for ft in 1:numpft
        if rad_params.xl[ft] < -0.4 || rad_params.xl[ft] > 0.6
            twostream_endrun("Leaf orientation factor (xl) should be between -0.4 and 0.6")
        end

        # Singularity protection: xl exactly 0 -> avmu blows up; nudge it.
        if abs(rad_params.xl[ft]) < 0.0001
            rad_params.xl[ft] = 0.0001
        end

        rad_params.phi1[ft] = 0.5 - 0.633 * rad_params.xl[ft] - 0.330 * rad_params.xl[ft] * rad_params.xl[ft]
        rad_params.phi2[ft] = 0.877 * (1.0 - 2.0 * rad_params.phi1[ft]) # 0 = horiz, 1 = vert leaves

        # Eq. 3.4 CLM50 Tech Man
        rad_params.avmu[ft] = (1.0 / rad_params.phi2[ft]) *
            (1.0 - (rad_params.phi1[ft] / rad_params.phi2[ft]) *
             log((rad_params.phi2[ft] + rad_params.phi1[ft]) / rad_params.phi1[ft]))

        for ib in 1:nbands
            rad_params.kd_leaf[ft] = rad_params.clumping_index[ft] / rad_params.avmu[ft]
            rad_params.kd_stem[ft] = 1.0

            rad_params.om_leaf[ib, ft] = rad_params.rhol[ib, ft] + rad_params.taul[ib, ft]
            rad_params.om_stem[ib, ft] = rad_params.rhos[ib, ft] + rad_params.taus[ib, ft]

            if rad_params.om_leaf[ib, ft] > 0.99
                twostream_endrun("RadParamPrep: extremely high leaf scattering coefficient (om = tau + rho > 0.99)")
            end
            if rad_params.om_stem[ib, ft] > 0.99
                twostream_endrun("RadParamPrep: extremely high stem scattering coefficient (om = tau + rho > 0.99)")
            end
        end
    end

    return nothing
end

# ================================================================================================

"""
    CanopyPrep!(this, frac_snow)

Pre-process things that change with canopy geometry or snow cover. Run only when
the canopy vegetation composition changes or the snow-cover fraction changes.
"""
function CanopyPrep!(this::twostream_type, frac_snow::Real)
    this.frac_snow = frac_snow

    if !this.force_prep
        if abs(this.frac_snow - this.frac_snow_old) < twostr_nearzero
            this.frac_snow_old = this.frac_snow
            return nothing
        end
    end

    this.frac_snow_old = this.frac_snow

    for ican in 1:this.n_lyr
        area_check = 0.0
        for icol in 1:this.n_col[ican]
            scelg = this.scelg[ican, icol]
            lai = scelg.lai
            sai = scelg.sai
            ft = scelg.pft

            vai = lai + sai

            # Mean element transmission coefficients w/o snow effects
            if ft == air_ft
                scelg.Kd = k_air
            else
                if twostr_debug && vai < twostr_nearzero
                    twostream_endrun("CanopyPrep: zero vai in non-air element")
                end
                scelg.Kd = (lai * rad_params.kd_leaf[ft] + sai * rad_params.kd_stem[ft]) / vai
            end

            area_check += scelg.area

            for ib in 1:this.n_bands
                scelb = this.band[ib].scelb[ican, icol]

                if ft == air_ft
                    scelb.om = om_air
                    scelb.betad = beta_air
                else
                    # Eq. 3.11 and 3.12 CLM5.0 Tech Man
                    om_veg = (lai * rad_params.om_leaf[ib, ft] + sai * rad_params.om_stem[ib, ft]) / vai

                    # Eq. 3.5 CLM5.0 Tech Man
                    scelb.om = this.frac_snow * om_snow[ib] + (1.0 - this.frac_snow) * om_veg

                    # Diffuse backscatter (from G. Bonan's code)
                    rho = (lai * rad_params.rhol[ib, ft] + sai * rad_params.rhos[ib, ft]) / vai
                    tau = (lai * rad_params.taul[ib, ft] + sai * rad_params.taus[ib, ft]) / vai

                    # Eq. 3.13 from CLM5.0 Tech Man
                    betad_veg = 0.5 / scelb.om *
                        (scelb.om + (rho - tau) * ((1.0 + rad_params.xl[ft]) / 2.0)^2)

                    # Eq. 3.6 from CLM5.0 Tech Man
                    betad_om = betad_veg * om_veg * (1.0 - this.frac_snow) +
                               om_snow[ib] * betad_snow[ib] * this.frac_snow

                    scelb.betad = betad_om / scelb.om

                    if twostr_debug && isnan(scelb.betad)
                        twostream_endrun("CanopyPrep: NaN diffuse backscatter (betad)")
                    end
                end

                a2 = scelg.Kd * scelg.Kd * (1.0 - scelb.om) *
                     (1.0 - scelb.om + 2.0 * scelb.om * scelb.betad)
                if a2 < 0.0
                    twostream_endrun("CanopyPrep: a^2 is less than zero")
                end

                # Avoid singularities (see Ad and Au, where a^2-Kb^2 is in the denominator)
                scelb.a = sqrt(a2)
            end
        end

        # Area-coverage check is disabled in the Fortran (the `if(.false.)`); keep it off.
    end

    return nothing
end

# ================================================================================================

"""
    ZenithPrep!(this, cosz_in)

Pre-process beam optical properties that change with the zenith angle. Must be
called AFTER `CanopyPrep!` (it relies on `scelb.om` and `scelb.a`).
"""
function ZenithPrep!(this::twostream_type, cosz_in::Real)
    if (cosz_in - 1.0) > twostr_nearzero
        twostream_endrun("ZenithPrep: the cosine of the zenith angle cannot exceed 1")
    elseif cosz_in < 0.0
        twostream_endrun("ZenithPrep: the cosine of the zenith angle should not be less than zero")
    end

    cosz = max(0.001, cosz_in)
    this.cosz = cosz

    Kb_stem_base = 1.0
    sing_tol = 0.01

    for ican in 1:this.n_lyr
        for icol in 1:this.n_col[ican]
            scelg = this.scelg[ican, icol]
            ft = scelg.pft

            local gdir::Float64, tmp0::Float64, tmp2::Float64

            if ft == air_ft
                # Simple provisions for a ghost element (air)
                scelg.Kb_leaf = k_air
                scelg.Kb = k_air
                gdir = 0.0
                tmp0 = 0.0
                tmp2 = 0.0
            else
                gdir = rad_params.phi1[ft] + rad_params.phi2[ft] * cosz

                Kb_stem = Kb_stem_base

                # How much direct light penetrates a single unit of lai?
                scelg.Kb_leaf = min(kb_max, rad_params.clumping_index[ft] * gdir / cosz)

                # Avoid singularities: ensure Kb =/ a. Identify the Kb_leaf that
                # gives a singularity and nudge it.
                if scelg.lai > twostr_nearzero
                    for ib in 1:this.n_bands
                        Kb_sing = (this.band[ib].scelb[ican, icol].a * (scelg.lai + scelg.sai) -
                                   scelg.sai * Kb_stem) / scelg.lai
                        if abs(scelg.Kb_leaf - Kb_sing) < sing_tol
                            scelg.Kb_leaf = Kb_sing + sing_tol
                        end
                    end
                else
                    for ib in 1:this.n_bands
                        Kb_sing = this.band[ib].scelb[ican, icol].a
                        if abs(Kb_stem - Kb_sing) < sing_tol
                            Kb_stem = Kb_sing + sing_tol
                        end
                    end
                end

                scelg.Kb = min(kb_max,
                    (scelg.lai * scelg.Kb_leaf + scelg.sai * Kb_stem) / (scelg.lai + scelg.sai))

                # Component terms for asu (single scattering albedo)
                tmp0 = gdir + rad_params.phi2[ft] * cosz
                tmp1 = rad_params.phi1[ft] * cosz
                tmp2 = 1.0 - tmp1 / tmp0 * log((tmp1 + tmp0) / tmp1)
            end

            for ib in 1:this.n_bands
                scelb = this.band[ib].scelb[ican, icol]

                if ft == air_ft
                    scelb.betab = beta_air
                else
                    # betab - upscatter parameter for direct beam radiation (G. Bonan)
                    # Eq. 3.16 CLM50 Tech Man
                    asu = 0.5 * gdir / tmp0 * tmp2

                    betab_veg = (1.0 + rad_params.avmu[ft] * scelg.Kb) /
                                (rad_params.avmu[ft] * scelg.Kb) * asu

                    om_veg = (scelg.lai * rad_params.om_leaf[ib, ft] +
                              scelg.sai * rad_params.om_stem[ib, ft]) / (scelg.lai + scelg.sai)

                    # Eq. 3.7 CLM50 Tech Man
                    betab_om = betab_veg * om_veg * (1.0 - this.frac_snow) +
                               om_snow[ib] * betab_snow[ib] * this.frac_snow

                    scelb.betab = betab_om / scelb.om

                    if twostr_debug && isnan(scelb.betab)
                        twostream_endrun("ZenithPrep: beam backscatter fraction is NaN")
                    end
                end
            end
        end
    end

    return nothing
end

# ================================================================================================

"""
    GetNSCel!(this)

Compute and store the total number of scattering elements.
"""
function GetNSCel!(this::twostream_type)
    this.n_scel = 0
    for ican in 1:this.n_lyr
        this.n_scel += this.n_col[ican]
    end
    return this.n_scel
end

# ===============================================================

"""
    Solve!(this, ib, upper_boundary_type, Rbeam_atm, Rdiff_atm, taulamb, omega, ipiv)
        -> (albedo_beam, albedo_diff, consv_err,
            frac_abs_can_beam, frac_abs_can_diff,
            frac_beam_grnd_beam, frac_diff_grnd_beam, frac_diff_grnd_diff)

Find the scattering coefficients for two-stream radiation in the canopy and the
canopy-top albedos, canopy-absorbed fractions, and ground-incident fractions for
band `ib`.

`taulamb`, `omega`, `ipiv` are caller-provided scratch arrays (sized at least
`2*n_scel`); they match the Fortran LAPACK-workspace arguments. `taulamb` is
overwritten with the linear-system solution; `ipiv` is unused in the Julia path
(retained for API parity).

The internal beam/diffuse scattering passes are solved with a normalized unit
incident in each stream, so the returned fractions and the internal
conservation check (`consv_err`) are only consistent when `Rbeam_atm` and
`Rdiff_atm` are passed as 1.0 (the normalized solution that ELM/CLM requests).
"""
function Solve!(this::twostream_type, ib::Integer,
                upper_boundary_type::Integer,
                Rbeam_atm::Real, Rdiff_atm::Real,
                taulamb::AbstractVector{<:Real},
                omega::AbstractMatrix{<:Real},
                ipiv::AbstractVector{<:Integer})

    # Outputs (initialized; the beam/diff isol passes each fill their half)
    albedo_beam = 0.0
    albedo_diff = 0.0
    frac_abs_can_beam = 0.0
    frac_abs_can_diff = 0.0
    frac_beam_grnd_beam = 0.0
    frac_diff_grnd_beam = 0.0
    frac_diff_grnd_diff = 0.0

    if (Rbeam_atm + Rdiff_atm) < twostr_nearzero
        twostream_endrun("Solve: no radiation; two stream should not have been called")
    end

    # ----------------------------------------------------------------------
    # Beam scattering: trivial solution that is a BC for the diffuse solver.
    # ----------------------------------------------------------------------
    Rbeam_top = 1.0
    for ican in 1:this.n_lyr
        Rbeam_bot = 0.0
        for icol in 1:this.n_col[ican]
            scelgp = this.scelg[ican, icol]
            scelbp = this.band[ib].scelb[ican, icol]
            scelbp.Rbeam0 = Rbeam_top
            Rbeam_bot += Rbeam_top * scelgp.area * exp(-scelgp.Kb * (scelgp.lai + scelgp.sai))
        end
        Rbeam_top = Rbeam_bot
    end

    # ----------------------------------------------------------------------
    # Element-level intermediate terms (B1, B2, Ad, Au)
    # ----------------------------------------------------------------------
    for ican in 1:this.n_lyr
        for icol in 1:this.n_col[ican]
            scelgp = this.scelg[ican, icol]
            scelbp = this.band[ib].scelb[ican, icol]

            b2 = -(scelgp.Kd * (1.0 - scelbp.om) * (1.0 - 2.0 * scelbp.betab) + scelgp.Kb) *
                 scelbp.om * scelgp.Kb * scelbp.Rbeam0

            b1 = -(scelgp.Kd * (1.0 - scelbp.om + 2.0 * scelbp.om * scelbp.betad) +
                   (1.0 - 2.0 * scelbp.betab) * scelgp.Kb) *
                 scelbp.om * scelgp.Kb * scelbp.Rbeam0

            nu_sqrd = (1.0 - scelbp.om) / (1.0 - scelbp.om + 2.0 * scelbp.om * scelbp.betad)

            if nu_sqrd < 0.0
                twostream_endrun("Solve: nu_sqrd is less than zero")
            end

            scelbp.B1 = 0.5 * (1.0 + sqrt(nu_sqrd))
            scelbp.B2 = 0.5 * (1.0 - sqrt(nu_sqrd))

            scelbp.Ad = -0.5 * (b1 + b2) / (scelbp.a * scelbp.a - scelgp.Kb * scelgp.Kb)
            scelbp.Au = -0.5 * (b1 - b2) / (scelbp.a * scelbp.a - scelgp.Kb * scelgp.Kb)
        end
    end

    # ----------------------------------------------------------------------
    # Set up and solve the linear systems: TAU = OMEGA*LAMBDA
    # ----------------------------------------------------------------------
    n_eq = 2 * this.n_scel

    # Two solutions: isol=1 beam only (for beam albedo), isol=2 diffuse only.
    for isol in 1:2

        if isol == 1
            this.band[ib].Rbeam_atm = 1.0
            this.band[ib].Rdiff_atm = 0.0
        else
            this.band[ib].Rbeam_atm = 0.0
            this.band[ib].Rdiff_atm = 1.0
        end

        # Zero the working slices
        @inbounds for j in 1:n_eq, i in 1:n_eq
            omega[i, j] = 0.0
        end
        @inbounds for i in 1:n_eq
            taulamb[i] = 0.0
        end

        # -------------------------------------------------------------------
        # I. Flux equations with the atmospheric boundary (upper canopy only)
        # -------------------------------------------------------------------
        qp = 0
        for icol in 1:this.n_col[1]
            scelbp = this.band[ib].scelb[1, icol]
            ilem = icol
            qp += 1
            k1 = 2 * (ilem - 1) + 1
            k2 = k1 + 1
            taulamb[qp] = this.band[ib].Rdiff_atm - this.band[ib].Rbeam_atm * scelbp.Ad
            omega[qp, k1] = scelbp.B2
            omega[qp, k2] = scelbp.B1
        end

        if this.n_lyr > 1
            # ----------------------------------------------------------------
            # II. Flux equations between canopy layers, DOWNWELLING
            # ----------------------------------------------------------------
            ilem_off = 0
            for ican in 2:this.n_lyr
                itop = ican - 1
                ibot = ican

                for jcol in 1:this.n_col[ibot]
                    qp += 1
                    ilem = ilem_off + this.n_col[itop] + jcol
                    k1 = 2 * (ilem - 1) + 1
                    k2 = k1 + 1

                    sb = this.band[ib].scelb[ibot, jcol]
                    taulamb[qp] = this.band[ib].Rbeam_atm * sb.Ad
                    omega[qp, k1] -= sb.B2
                    omega[qp, k2] -= sb.B1

                    for icol in 1:this.n_col[itop]
                        ilem2 = ilem_off + icol
                        k1 = 2 * (ilem2 - 1) + 1
                        k2 = k1 + 1
                        scelgp = this.scelg[itop, icol]
                        scelbp = this.band[ib].scelb[itop, icol]
                        vai = scelgp.lai + scelgp.sai

                        taulamb[qp] -= scelgp.area * this.band[ib].Rbeam_atm * scelbp.Ad * exp(-scelgp.Kb * vai)
                        omega[qp, k1] += scelgp.area * scelbp.B2 * exp(scelbp.a * vai)
                        omega[qp, k2] += scelgp.area * scelbp.B1 * exp(-scelbp.a * vai)
                    end
                end

                ilem_off += this.n_col[itop]
            end

            # ----------------------------------------------------------------
            # III. Flux equations between canopy layers, UPWELLING
            # ----------------------------------------------------------------
            ilem_off = 0
            for ican in 2:this.n_lyr
                itop = ican - 1
                ibot = ican

                for icol in 1:this.n_col[itop]
                    qp += 1

                    ilem = ilem_off + icol
                    k1 = 2 * (ilem - 1) + 1
                    k2 = k1 + 1
                    scelgp = this.scelg[itop, icol]
                    scelbp = this.band[ib].scelb[itop, icol]

                    vai = scelgp.lai + scelgp.sai
                    taulamb[qp] = this.band[ib].Rbeam_atm * scelbp.Au * exp(-scelgp.Kb * vai)
                    omega[qp, k1] -= scelbp.B1 * exp(scelbp.a * vai)
                    omega[qp, k2] -= scelbp.B2 * exp(-scelbp.a * vai)

                    for jcol in 1:this.n_col[ibot]
                        ilem2 = ilem_off + this.n_col[itop] + jcol
                        k1 = 2 * (ilem2 - 1) + 1
                        k2 = k1 + 1
                        scelgp = this.scelg[ibot, jcol]
                        scelbp = this.band[ib].scelb[ibot, jcol]

                        taulamb[qp] -= this.band[ib].Rbeam_atm * scelgp.area * scelbp.Au
                        omega[qp, k1] += scelgp.area * scelbp.B1
                        omega[qp, k2] += scelgp.area * scelbp.B2
                    end
                end

                ilem_off += this.n_col[itop]
            end
        end

        # -------------------------------------------------------------------
        # Flux balance between understory elements and the ground below them
        # -------------------------------------------------------------------
        ilem_off = 0
        for ican in 1:(this.n_lyr - 1)
            ilem_off += this.n_col[ican]
        end

        for jcol in 1:this.n_col[this.n_lyr]
            ilem = ilem_off + jcol
            qp += 1
            k1 = 2 * (ilem - 1) + 1
            k2 = k1 + 1

            scelgp = this.scelg[this.n_lyr, jcol]
            scelbp = this.band[ib].scelb[this.n_lyr, jcol]

            vai = scelgp.lai + scelgp.sai

            taulamb[qp] = this.band[ib].Rbeam_atm * (scelbp.Au * exp(-scelgp.Kb * vai) -
                this.band[ib].albedo_grnd_diff * scelbp.Ad * exp(-scelgp.Kb * vai) -
                this.band[ib].albedo_grnd_beam * scelbp.Rbeam0 * exp(-scelgp.Kb * vai))

            omega[qp, k1] -= scelbp.B1 * exp(scelbp.a * vai)
            omega[qp, k2] -= scelbp.B2 * exp(-scelbp.a * vai)

            omega[qp, k1] += this.band[ib].albedo_grnd_diff * scelbp.B2 * exp(scelbp.a * vai)
            omega[qp, k2] += this.band[ib].albedo_grnd_diff * scelbp.B1 * exp(-scelbp.a * vai)
        end

        # -------------------------------------------------------------------
        # Solve. dgesv overwrites TAU with LAMBDA. We use Julia's LU/ldiv.
        # -------------------------------------------------------------------
        A = @view omega[1:n_eq, 1:n_eq]
        rhs = @view taulamb[1:n_eq]

        if twostr_debug
            tau_temp = Vector{Float64}(rhs)
            omega_temp = Matrix{Float64}(A)
        end

        local fac
        try
            fac = lu!(Matrix{Float64}(A))
        catch
            twostream_endrun("Solve: could not find a solution via LU (singular matrix)")
        end
        sol = Vector{Float64}(rhs)
        ldiv!(fac, sol)
        @inbounds for i in 1:n_eq
            taulamb[i] = sol[i]
        end

        # Forward check on the solution error
        if twostr_debug
            for ilem in 1:n_eq
                acc = 0.0
                for j in 1:n_eq
                    acc += taulamb[j] * omega_temp[ilem, j]
                end
                temp_err = tau_temp[ilem] - acc
                if abs(temp_err) > rel_err_thresh
                    twostream_endrun("Solve: poor forward solution on two-stream solver")
                end
            end
        end

        # Save the solution terms
        ilem_off = 0
        if isol == 1  # Beam
            for ican in 1:this.n_lyr
                for icol in 1:this.n_col[ican]
                    ilem = ilem_off + icol
                    k1 = 2 * (ilem - 1) + 1
                    k2 = k1 + 1
                    scelbp = this.band[ib].scelb[ican, icol]
                    scelbp.lambda1_beam = taulamb[k1]
                    scelbp.lambda2_beam = taulamb[k2]
                    # diff terms multiplied by zero before use; set to zero to avoid NaN
                    scelbp.lambda1_diff = 0.0
                    scelbp.lambda2_diff = 0.0
                end
                ilem_off += this.n_col[ican]
            end
        else
            for ican in 1:this.n_lyr
                for icol in 1:this.n_col[ican]
                    ilem = ilem_off + icol
                    k1 = 2 * (ilem - 1) + 1
                    k2 = k1 + 1
                    scelbp = this.band[ib].scelb[ican, icol]
                    scelbp.lambda1_diff = taulamb[k1]
                    scelbp.lambda2_diff = taulamb[k2]
                end
                ilem_off += this.n_col[ican]
            end
        end

        # -------------------------------------------------------------------
        # Process total canopy absorbed radiation + ground-interface fluxes
        # -------------------------------------------------------------------
        if isol == 1  # Beam
            ican = 1
            albedo_beam = 0.0
            for icol in 1:this.n_col[ican]
                scelgp = this.scelg[ican, icol]
                albedo_beam += scelgp.area * GetRdUp(this, ican, icol, ib, 0.0)
            end

            frac_diff_grnd_beam = 0.0
            frac_beam_grnd_beam = 0.0
            ican = this.n_lyr
            for icol in 1:this.n_col[ican]
                scelgp = this.scelg[ican, icol]
                scelbp = this.band[ib].scelb[ican, icol]
                frac_diff_grnd_beam += scelgp.area * GetRdDn(this, ican, icol, ib, scelgp.lai + scelgp.sai)
                frac_beam_grnd_beam += scelgp.area * scelbp.Rbeam0 * exp(-scelgp.Kb * (scelgp.lai + scelgp.sai))
            end

            frac_abs_can_beam = 0.0
            for ican in 1:this.n_lyr
                for icol in 1:this.n_col[ican]
                    scelgp = this.scelg[ican, icol]
                    (rb_abs, rd_abs, _, _, _, _, _) =
                        GetAbsRad(this, ican, icol, ib, 0.0, scelgp.lai + scelgp.sai)
                    frac_abs_can_beam += scelgp.area * (rb_abs + rd_abs)
                end
            end
        else  # Diffuse
            albedo_diff = 0.0
            for icol in 1:this.n_col[1]
                scelgp = this.scelg[1, icol]
                albedo_diff += scelgp.area * GetRdUp(this, 1, icol, ib, 0.0)
            end

            frac_abs_can_diff = 0.0
            for ican in 1:this.n_lyr
                for icol in 1:this.n_col[ican]
                    scelgp = this.scelg[ican, icol]
                    (_, rd_abs, _, _, _, _, _) =
                        GetAbsRad(this, ican, icol, ib, 0.0, scelgp.lai + scelgp.sai)
                    frac_abs_can_diff += scelgp.area * rd_abs
                end
            end

            frac_diff_grnd_diff = 0.0
            ican = this.n_lyr
            for icol in 1:this.n_col[ican]
                scelgp = this.scelg[ican, icol]
                frac_diff_grnd_diff += scelgp.area * GetRdDn(this, ican, icol, ib, scelgp.lai + scelgp.sai)
            end
        end

    end # isol

    # ----------------------------------------------------------------------
    # Canopy radiation balance conservation check
    # Source = upwelling + canopy absorbed + ground absorbed
    # ----------------------------------------------------------------------
    consv_err = ((Rbeam_atm + Rdiff_atm) -
                 (albedo_diff + albedo_beam) -
                 (frac_abs_can_diff + frac_abs_can_beam) -
                 ((frac_diff_grnd_diff + frac_diff_grnd_beam) * (1.0 - this.band[ib].albedo_grnd_diff)) -
                 (frac_beam_grnd_beam * (1.0 - this.band[ib].albedo_grnd_beam))) / (Rbeam_atm + Rdiff_atm)

    consv_err = abs(consv_err)

    if consv_err > rel_err_thresh
        twostream_endrun("Total canopy flux balance not closing in TwoStreamMLPEMod:Solve")
    end

    # ----------------------------------------------------------------------
    # Re-cast albedos as a direct result of the components (albedo_corr=.true.)
    # ----------------------------------------------------------------------
    albedo_beam = Rbeam_atm - (frac_abs_can_beam +
        frac_diff_grnd_beam * (1.0 - this.band[ib].albedo_grnd_diff) +
        frac_beam_grnd_beam * (1.0 - this.band[ib].albedo_grnd_beam))

    albedo_diff = Rdiff_atm - (frac_abs_can_diff +
        frac_diff_grnd_diff * (1.0 - this.band[ib].albedo_grnd_diff))

    # ----------------------------------------------------------------------
    # Restore boundary conditions
    # ----------------------------------------------------------------------
    if upper_boundary_type == normalized_upper_boundary
        this.band[ib].Rbeam_atm = twostr_unset_r8
        this.band[ib].Rdiff_atm = twostr_unset_r8
    else
        this.band[ib].Rbeam_atm = Rbeam_atm
        this.band[ib].Rdiff_atm = Rdiff_atm
    end

    return (albedo_beam, albedo_diff, consv_err,
            frac_abs_can_beam, frac_abs_can_diff,
            frac_beam_grnd_beam, frac_diff_grnd_beam, frac_diff_grnd_diff)
end
