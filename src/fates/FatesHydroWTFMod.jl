# FatesHydroWTFMod.jl
# Julia port of FATES src/fates/biogeophys/FatesHydroWTFMod.F90
#
# Plant-hydraulics Water Transfer Functions (WTFs, a.k.a. pedotransfer functions).
# These come in two flavors:
#   * Water Retention Functions (WRF): map between volumetric water content (theta)
#     and matric potential (psi), plus the analytic derivative dpsi/dth.
#   * Water Conductance Functions (WKF): map matric potential (psi) to the fractional
#     loss of total conductance (FLC, here called "ftc"), plus dftc/dpsi.
#
# The Fortran uses an abstract `wrf_type`/`wkf_type` base + polymorphic `class(...)`
# methods. Here we use Julia abstract types (`WRFType`/`WKFType`) + concrete mutable
# structs + multiple dispatch. The Fortran method names are preserved:
#   th_from_psi, psi_from_th, dpsidth_from_th, set_wrf_param!, get_thsat
#   ftc_from_psi, dftcdpsi_from_psi, set_wkf_param!
#
# A WKF holds a reference to its matching WRF (Fortran `this%wrf` pointer) only so it
# can read `psi_min`. The analytic derivative formulas are preserved EXACTLY because
# they feed downstream AD / Newton solves.
#
# Deps: FatesConstantsMod (fates_r8, fates_unset_r8, nearzero), FatesGlobals (fates_endrun).

# ---------------------------------------------------------------------------
# Module-level parameters (Fortran `real(r8), parameter`)
# ---------------------------------------------------------------------------
const wtf_min_ftc = 0.0          # Minimum allowed fraction of total conductance
const min_ftc_scalar = 2.0       # Attenuation of the weighting used to impose ftc-min
                                 # at the minimum allowable psi. A value of two yields a
                                 # weighting factor of ~0.175 after 1 MPa, ~0.025 after 2 MPa.

# Bounds on saturated fraction, outside of which we use linear PV or stop flow.
# Saturated fraction is defined from the volumetric WC "th" and residual/saturation
# "th_res"/"th_sat": (th-th_res)/(th_sat-th_res).
const min_sf_interp = 0.01       # Linear interpolation below this saturated frac
const max_sf_interp = 0.998      # Linear interpolation above this saturated frac

const quad_a1 = 0.80             # Smoothing "A" term in the capillary-elastic region
const quad_a2 = 0.99             # Smoothing "A" term in the elastic-cavitation region

const min_psi_cch = -15.0        # Minimum suction (MPa) for Campbell/Clapp-Hornberger

# ---------------------------------------------------------------------------
# Abstract base types
# ---------------------------------------------------------------------------

"""
    WRFType

Abstract base for water-retention functions (theta <-> psi). Mirrors the Fortran
`wrf_type`. Concrete subtypes each carry the shared endpoint-interpolation fields:
`psi_max`, `psi_min`, `dpsidth_max`, `dpsidth_min`, `th_min`, `th_max`.
"""
abstract type WRFType end

"""
    WKFType

Abstract base for water-conductance functions (psi -> ftc). Mirrors the Fortran
`wkf_type`. Concrete subtypes carry a `wrf` reference to the matching WRFType so they
can read `psi_min`.
"""
abstract type WKFType end

# ===========================================================================
# Generic WRF helpers (Fortran `non_overridable` procedures), usable by all WRFs.
# These are linear extrapolations for solutions outside the expected range.
# ===========================================================================

"""
    set_min_max_from_satres!(this, th_res, th_sat)

Use `max_sf_interp`/`min_sf_interp` to define where the linear ranges start and stop.
Only valid for functions defined by a saturation and a residual value.
"""
function set_min_max_from_satres!(this::WRFType, th_res::Float64, th_sat::Float64)
    this.th_max      = max_sf_interp * (th_sat - th_res) + th_res
    this.th_min      = min_sf_interp * (th_sat - th_res) + th_res
    this.psi_max     = psi_from_th(this, this.th_max)
    this.dpsidth_max = dpsidth_from_th(this, this.th_max)
    this.psi_min     = psi_from_th(this, this.th_min)
    this.dpsidth_min = dpsidth_from_th(this, this.th_min)
    return nothing
end

"Calculate psi in the linear range below residual."
psi_linear_res(this::WRFType, th::Float64) =
    this.psi_min + this.dpsidth_min * (th - this.th_min)

"Calculate psi in the linear range above saturation."
psi_linear_sat(this::WRFType, th::Float64) =
    this.psi_max + this.dpsidth_max * (th - this.th_max)

"Calculate th from psi in the linear range above saturation."
th_linear_sat(this::WRFType, psi::Float64) =
    this.th_max + (psi - this.psi_max) / this.dpsidth_max

"Calculate th from psi in the linear range below residual."
th_linear_res(this::WRFType, psi::Float64) =
    this.th_min + (psi - this.psi_min) / this.dpsidth_min

get_thmin(this::WRFType) = this.th_min

# ===========================================================================

"""
    get_min_ftc_weight(psi_min, psi) -> (min_ftc_weight, dmin_ftc_weight_dpsi)

Weighting factor used to generate a smooth curve imposing that min_ftc happens at
min_psi. Usable with any WTF. Returns the weight and its derivative wrt psi.
"""
function get_min_ftc_weight(psi_min::Float64, psi::Float64)
    # If the difference between psi and psi_min exceeds 10 MPa, assume no effect (weight 0).
    arg = max(psi_min - psi, -10.0)
    min_ftc_weight = exp(min_ftc_scalar * arg)
    dmin_ftc_weight_dpsi = -min_ftc_scalar * exp(min_ftc_scalar * arg)

    if min_ftc_weight >= 1.0
        min_ftc_weight = 1.0
        dmin_ftc_weight_dpsi = 0.0
    elseif min_ftc_weight <= 0.0
        min_ftc_weight = 0.0
        dmin_ftc_weight_dpsi = 0.0
    end

    return min_ftc_weight, dmin_ftc_weight_dpsi
end

# ===========================================================================
# Base methods — these should never be actualized (the abstract type is never
# pointed to directly). They mirror the Fortran base routines that endrun.
# ===========================================================================

const _wtf_base_err = "The base water transfer function should never be actualized; " *
                      "check how the class pointer was setup"

set_wrf_param!(::WRFType, ::AbstractVector{<:Real}) = fates_endrun(_wtf_base_err)
get_thsat(::WRFType)                                = fates_endrun(_wtf_base_err)
th_from_psi(::WRFType, ::Real)                      = fates_endrun(_wtf_base_err)
psi_from_th(::WRFType, ::Real)                      = fates_endrun(_wtf_base_err)
dpsidth_from_th(::WRFType, ::Real)                  = fates_endrun(_wtf_base_err)

set_wkf_param!(::WKFType, ::AbstractVector{<:Real}) = fates_endrun(_wtf_base_err)
ftc_from_psi(::WKFType, ::Real)                     = fates_endrun(_wtf_base_err)
dftcdpsi_from_psi(::WKFType, ::Real)                = fates_endrun(_wtf_base_err)

# ===========================================================================
# Van Genuchten WTF
# ===========================================================================

"Van Genuchten water retention function."
Base.@kwdef mutable struct wrf_type_vg <: WRFType
    # Shared WRFType endpoint fields
    psi_max::Float64     = 0.0
    psi_min::Float64     = 0.0
    dpsidth_max::Float64 = 0.0
    dpsidth_min::Float64 = 0.0
    th_min::Float64      = 0.0
    th_max::Float64      = 0.0
    # VG params
    alpha::Float64  = 0.0   # Inverse air entry parameter [m3/MPa]
    n_vg::Float64   = 0.0   # pore size distribution parameter (psd)
    m_vg::Float64   = 0.0   # m in van Genuchten 1980
    th_sat::Float64 = 0.0   # Saturation volumetric water content [m3/m3]
    th_res::Float64 = 0.0   # Residual volumetric water content   [m3/m3]
end

"Van Genuchten water conductivity function."
Base.@kwdef mutable struct wkf_type_vg <: WKFType
    wrf::Union{WRFType,Nothing} = nothing  # matching WRF (for psi_min)
    alpha::Float64  = 0.0   # Inverse air entry parameter [m3/MPa]
    n_vg::Float64   = 0.0   # pore size distribution parameter
    m_vg::Float64   = 0.0   # m in van Genuchten 1980
    tort::Float64   = 0.0   # Tortuosity parameter (sometimes "l")
    th_sat::Float64 = 0.0   # Saturation volumetric water content [m3/m3]
    th_res::Float64 = 0.0   # Residual volumetric water content   [m3/m3]
end

function set_wrf_param!(this::wrf_type_vg, params_in::AbstractVector{<:Real})
    this.alpha  = params_in[1]
    this.n_vg   = params_in[2]
    this.m_vg   = params_in[3]
    this.th_sat = params_in[4]
    this.th_res = params_in[5]
    set_min_max_from_satres!(this, this.th_res, this.th_sat)
    return nothing
end

function set_wkf_param!(this::wkf_type_vg, params_in::AbstractVector{<:Real})
    this.alpha  = params_in[1]
    this.n_vg   = params_in[2]
    this.m_vg   = params_in[3]
    this.th_sat = params_in[4]
    this.th_res = params_in[5]
    this.tort   = params_in[6]
    return nothing
end

get_thsat(this::wrf_type_vg) = this.th_sat

function th_from_psi(this::wrf_type_vg, psi::Real)
    # Van Genuchten (1980): theta from matric potential.
    m = this.m_vg
    n = this.n_vg

    if psi > this.psi_max
        th = th_linear_sat(this, psi)            # Linear range for extreme values
    elseif psi < this.psi_min
        th = th_linear_res(this, psi)            # Linear range for extreme values
    else
        satfrac = (1.0 + (-this.alpha * psi)^n)^(-m)
        th = satfrac * (this.th_sat - this.th_res) + this.th_res
    end
    return th
end

function psi_from_th(this::wrf_type_vg, th::Real)
    # Van Genuchten (1980): matric potential from theta (inverted).
    if th > this.th_max
        psi = psi_linear_sat(this, th)
    elseif th < this.th_min
        psi = psi_linear_res(this, th)
    else
        m = this.m_vg
        n = this.n_vg
        satfrac = (th - this.th_res) / (this.th_sat - this.th_res)
        psi = -(1.0 / this.alpha) * (satfrac^(1.0 / (-m)) - 1.0)^(1 / n)
    end
    return psi
end

function dpsidth_from_th(this::wrf_type_vg, th::Real)
    a1 = 1.0 / this.alpha
    m1 = 1.0 / this.n_vg
    m2 = -1.0 / this.m_vg

    if th > this.th_max
        dpsidth = this.dpsidth_max
    elseif th < this.th_min
        dpsidth = this.dpsidth_min
    else
        satfrac = (th - this.th_res) / (this.th_sat - this.th_res)
        dsatfrac_dth = 1.0 / (this.th_sat - this.th_res)
        # psi = -a1 * (satfrac**m2 - 1)** m1
        # dpsidth = -m1*a1*(satfrac**m2-1)**(m1-1) * m2*satfrac**(m2-1)*dsatfracdth
        dpsidth = -m1 * a1 * (satfrac^m2 - 1.0)^(m1 - 1.0) *
                  m2 * satfrac^(m2 - 1.0) * dsatfrac_dth
    end
    return dpsidth
end

function ftc_from_psi(this::wkf_type_vg, psi::Real)
    psi_min = this.wrf.psi_min
    n = this.n_vg
    m = this.m_vg

    if psi < 0.0
        # VG 1980 assumes a positive pressure convention...
        psi_eff = -psi
        num = (1.0 - ((this.alpha * psi_eff)^n /
                      (1.0 + (this.alpha * psi_eff)^n))^m)^2.0
        den = (1.0 + (this.alpha * psi_eff)^n)^(this.tort * m)

        ftc = min(1.0, num / den)

        if ftc <= wtf_min_ftc
            ftc = wtf_min_ftc
        else
            # Add protections and ensure no conductance at incredibly low suction.
            min_ftc_weight, _ = get_min_ftc_weight(psi_min, Float64(psi))
            ftc = ftc * (1.0 - min_ftc_weight) + wtf_min_ftc * min_ftc_weight
        end
    else
        ftc = 1.0
    end
    return ftc
end

function dftcdpsi_from_psi(this::wkf_type_vg, psi::Real)
    # Derivative of fraction of total conductivity. Broken into terms; see tech note.
    psi_min = this.wrf.psi_min
    n = this.n_vg
    m = this.m_vg

    if psi >= 0.0
        dftcdpsi = 0.0
    else
        psi_eff = -psi  # switch VG 1980 convention

        ftc = ftc_from_psi(this, psi)

        if abs(ftc - wtf_min_ftc) < nearzero
            dftcdpsi = 0.0   # We cap ftc, so derivative is zero
        else
            t1  = (this.alpha * psi_eff)^(n * m)
            dt1 = this.alpha * (n * m) * (this.alpha * psi_eff)^(n * m - 1.0)

            t2  = (1.0 + (this.alpha * psi_eff)^n)^(-m)
            dt2 = (-m) *
                  (1.0 + (this.alpha * psi_eff)^n)^(-m - 1.0) *
                  n * (this.alpha * psi_eff)^(n - 1.0) * this.alpha

            t3  = (1.0 + (this.alpha * psi_eff)^n)^(this.tort * m)
            dt3 = this.tort * m *
                  (1.0 + (this.alpha * psi_eff)^n)^(this.tort * m - 1.0) *
                  n * (this.alpha * psi_eff)^(n - 1.0) * this.alpha

            dftcdpsi = 2.0 * (1.0 - t1 * t2) * (t1 * dt2 + t2 * dt1) / t3 -
                       t3^(-2.0) * dt3 * (1.0 - t1 * t2)^2.0

            min_ftc_weight, dmin_ftc_weight_dpsi = get_min_ftc_weight(psi_min, Float64(psi))

            # differentiate: ftc = ftc*(1 - w) + min_ftc*w
            dftcdpsi = dftcdpsi -
                       (dftcdpsi * min_ftc_weight + ftc * dmin_ftc_weight_dpsi) +
                       wtf_min_ftc * dmin_ftc_weight_dpsi
        end
    end
    return dftcdpsi
end

# ===========================================================================
# Campbell / Clapp-Hornberger (CCH)
# ===========================================================================

"Campbell/Clapp-Hornberger water retention function."
Base.@kwdef mutable struct wrf_type_cch <: WRFType
    psi_max::Float64     = 0.0
    psi_min::Float64     = 0.0
    dpsidth_max::Float64 = 0.0
    dpsidth_min::Float64 = 0.0
    th_min::Float64      = 0.0
    th_max::Float64      = 0.0
    th_sat::Float64  = 0.0   # Saturation volumetric water content [m3/m3]
    psi_sat::Float64 = 0.0   # Bubbling pressure (potential at saturation) [MPa]
    beta::Float64    = 0.0   # Clapp-Hornberger "beta" parameter [-]
end

"Campbell/Clapp-Hornberger water conductivity function."
Base.@kwdef mutable struct wkf_type_cch <: WKFType
    wrf::Union{WRFType,Nothing} = nothing
    th_sat::Float64  = 0.0
    psi_sat::Float64 = 0.0
    beta::Float64    = 0.0
end

function set_wrf_param!(this::wrf_type_cch, params_in::AbstractVector{<:Real})
    this.th_sat  = params_in[1]
    this.psi_sat = params_in[2]
    this.beta    = params_in[3]

    # Set DERIVED constants used for interpolating in extreme ranges.
    this.th_max  = max_sf_interp * this.th_sat
    this.psi_min = min_psi_cch
    # Temporary th_min (can't be uninitialized while computing psi_max).
    this.th_min  = 0.001

    this.psi_max     = psi_from_th(this, this.th_max)
    this.dpsidth_max = dpsidth_from_th(this, this.th_max)

    # Get the actual th_min equivalent to psi_min.
    this.th_min      = th_from_psi(this, min_psi_cch)
    this.dpsidth_min = dpsidth_from_th(this, this.th_min)
    return nothing
end

function set_wkf_param!(this::wkf_type_cch, params_in::AbstractVector{<:Real})
    this.th_sat  = params_in[1]
    this.psi_sat = params_in[2]
    this.beta    = params_in[3]
    return nothing
end

get_thsat(this::wrf_type_cch) = this.th_sat

function th_from_psi(this::wrf_type_cch, psi::Real)
    if psi > this.psi_max
        th = this.th_max + (psi - this.psi_max) / this.dpsidth_max  # Linear range
    else
        if psi < this.psi_min
            th = this.th_sat * (this.psi_min / this.psi_sat)^(-1.0 / this.beta)
        else
            th = this.th_sat * (psi / this.psi_sat)^(-1.0 / this.beta)
        end
    end
    return th
end

function psi_from_th(this::wrf_type_cch, th::Real)
    if th > this.th_max
        psi = this.psi_max + this.dpsidth_max * (th - max_sf_interp * this.th_sat)
    else
        if th < this.th_min
            psi = this.psi_sat * (this.th_min / this.th_sat)^(-this.beta)
        else
            psi = this.psi_sat * (th / this.th_sat)^(-this.beta)
        end
    end
    return psi
end

function dpsidth_from_th(this::wrf_type_cch, th::Real)
    if th > this.th_max
        dpsidth = this.dpsidth_max
    else
        if th < this.th_min
            # Edge case; th is capped at th_min normally.
            dpsidth = -this.beta * this.psi_sat / this.th_sat *
                      (this.th_min / this.th_sat)^(-this.beta - 1.0)
        else
            dpsidth = -this.beta * this.psi_sat / this.th_sat *
                      (th / this.th_sat)^(-this.beta - 1.0)
        end
    end
    return dpsidth
end

function ftc_from_psi(this::wkf_type_cch, psi::Real)
    psi_min = this.wrf.psi_min
    # ftc = (psi/psi_sat)^(-2-3/b)
    psi_eff = min(psi, this.psi_sat)
    ftc = (psi_eff / this.psi_sat)^(-2.0 - 3.0 / this.beta)

    if ftc <= wtf_min_ftc
        ftc = wtf_min_ftc
    else
        min_ftc_weight, _ = get_min_ftc_weight(psi_min, Float64(psi))
        ftc = ftc * (1.0 - min_ftc_weight) + wtf_min_ftc * min_ftc_weight
    end
    return ftc
end

function dftcdpsi_from_psi(this::wkf_type_cch, psi::Real)
    psi_min = this.wrf.psi_min
    # Capped FTC=1.0 at saturation => derivative zero there.
    if psi < this.psi_sat
        ftc = ftc_from_psi(this, psi)
        if abs(ftc - wtf_min_ftc) < nearzero
            dftcdpsi = 0.0
        else
            dftcdpsi = (-2.0 - 3.0 / this.beta) / this.psi_sat *
                       (psi / this.psi_sat)^(-3.0 - 3.0 / this.beta)

            min_ftc_weight, dmin_ftc_weight_dpsi = get_min_ftc_weight(psi_min, Float64(psi))

            dftcdpsi = dftcdpsi -
                       (dftcdpsi * min_ftc_weight + ftc * dmin_ftc_weight_dpsi) +
                       wtf_min_ftc * dmin_ftc_weight_dpsi
        end
    else
        dftcdpsi = 0.0
    end
    return dftcdpsi
end

# ===========================================================================
# Type1 smooth approximation of Campbell/Clapp-Hornberger
# Bisht et al. Geosci. Model Dev., 11, 4085-4102, 2018
# ===========================================================================

"Smooth Campbell/Clapp-Hornberger water retention function."
Base.@kwdef mutable struct wrf_type_smooth_cch <: WRFType
    psi_max::Float64     = 0.0
    psi_min::Float64     = 0.0
    dpsidth_max::Float64 = 0.0
    dpsidth_min::Float64 = 0.0
    th_min::Float64      = 0.0
    th_max::Float64      = 0.0
    th_sat::Float64  = 0.0
    psi_sat::Float64 = 0.0
    beta::Float64    = 0.0
    scch_pu::Float64 = 0.0   # breakpoint capillary pressure (lower smoothing limit) [MPa]
    scch_ps::Float64 = 0.0   # breakpoint capillary pressure (upper smoothing limit) [MPa]
    scch_b2::Float64 = 0.0   # quadratic-term coefficient in smoothing polynomial [-]
    scch_b3::Float64 = 0.0   # cubic-term coefficient in smoothing polynomial [-]
end

"Smooth Campbell/Clapp-Hornberger water conductivity function."
Base.@kwdef mutable struct wkf_type_smooth_cch <: WKFType
    wrf::Union{WRFType,Nothing} = nothing
    th_sat::Float64  = 0.0
    psi_sat::Float64 = 0.0
    beta::Float64    = 0.0
    scch_pu::Float64 = 0.0
    scch_ps::Float64 = 0.0
    scch_b2::Float64 = 0.0
    scch_b3::Float64 = 0.0
end

"""
    findGu_SBC_zeroCoeff(lambda, AA, gs) -> gu

Find `pu` (as a multiplier of `pc0`) that forces a coefficient of the smoothing cubic
to zero (b2=0 for AA=3, b3=0 for AA=2). Bracketed Newton-Raphson. Bisht et al. 2018.
"""
function findGu_SBC_zeroCoeff(lambda::Float64, AA::Integer, gs::Float64)
    relTol = 1.0e-12

    if lambda <= 0.0 || lambda >= 2.0 ||
       (AA != 2 && AA != 3) ||
       gs >= 1.0 || gs < 0.0
        fates_endrun("findGu_SBC_zeroCoeff: bad param")
    end

    # Initialize (solution if gs = 0).
    gu = (AA / (AA + lambda))^(-1.0 / lambda)

    if gs > 0.0
        guLeft = 1.0
        guRight = gu

        iter = 0
        deltaGu = 1.0e20  # something large
        while abs(deltaGu) >= relTol * abs(gu)
            iter += 1

            # Reset gu using bisection if necessary.
            if gu <= guLeft || gu >= guRight
                gu = guLeft + 0.5 * (guRight - guLeft)
            end

            # Find residual.
            guInv = 1.0 / gu
            guToMinusLam = gu^(-lambda)
            gsOnGu = gs * guInv
            resid = AA - guToMinusLam * (AA + lambda - lambda * gsOnGu)

            # Update bracket.
            if resid < 0.0
                guLeft = gu
            else
                guRight = gu
            end

            # Newton-Raphson step.
            dr_dGu = (1.0 + lambda) * (1.0 - gsOnGu) + (AA - 1)
            dr_dGu = lambda * guToMinusLam * guInv * dr_dGu
            deltaGu = resid / dr_dGu
            gu = gu - deltaGu

            if iter > 10000
                fates_endrun("findGu_SBC_zeroCoeff iteration not converging")
            end
        end
    end

    if gu != gu
        fates_endrun("gu = nan in findGu_SBC_zeroCoeff")
    end
    return gu
end

# Shared parameter setup for both smooth-CCH WRF and WKF. Returns the derived
# (scch_ps, scch_pu, scch_b2, scch_b3) tuple.
function _set_smooth_cch_coeffs(psi_sat::Float64, beta::Float64, styp::Integer)
    alpha  = -1.0 / psi_sat
    lambda = 1.0 / beta
    ps = -0.9 / alpha
    scch_ps = ps

    if styp == 1
        # Choose pu that forces scch_b2 = 0.
        pu = findGu_SBC_zeroCoeff(lambda, 3, -alpha * ps) / (-alpha)
        scch_pu = pu

        bcAtPu            = (-alpha * pu)^(-lambda)
        lambdaDeltaPuOnPu = lambda * (1.0 - ps / pu)
        oneOnDeltaPu      = 1.0 / (pu - ps)

        scch_b2 = 0.0
        scch_b3 = (2.0 - bcAtPu * (2.0 + lambdaDeltaPuOnPu)) *
                  oneOnDeltaPu * oneOnDeltaPu * oneOnDeltaPu
        if scch_b3 <= 0.0
            fates_endrun("set_wrf_param_smooth_cch b3 <=0")
        end
    else
        # Choose pu that forces scch_b3 = 0.
        pu = findGu_SBC_zeroCoeff(lambda, 2, -alpha * ps) / (-alpha)
        scch_pu = pu

        bcAtPu            = (-alpha * pu)^(-lambda)
        lambdaDeltaPuOnPu = lambda * (1.0 - ps / pu)
        oneOnDeltaPu      = 1.0 / (pu - ps)

        scch_b2 = -(3.0 - bcAtPu * (3.0 + lambdaDeltaPuOnPu)) * oneOnDeltaPu * oneOnDeltaPu
        if scch_b2 >= 0.0
            fates_endrun("set_wrf_param_smooth_cch b2 <= 0")
        end
        scch_b3 = 0.0
    end
    return scch_ps, scch_pu, scch_b2, scch_b3
end

function set_wrf_param!(this::wrf_type_smooth_cch, params_in::AbstractVector{<:Real})
    this.th_sat  = params_in[1]
    this.psi_sat = params_in[2]
    this.beta    = params_in[3]
    styp = Int(params_in[4])

    this.scch_ps, this.scch_pu, this.scch_b2, this.scch_b3 =
        _set_smooth_cch_coeffs(this.psi_sat, this.beta, styp)

    # Set DERIVED constants used for interpolating in extreme ranges.
    this.th_max      = max_sf_interp * this.th_sat
    tinyth           = eps(this.th_max)  # Fortran tiny(); a tiny offset below th_max
    this.psi_max     = psi_from_th(this, this.th_max - tinyth)
    this.dpsidth_max = dpsidth_from_th(this, this.th_max - tinyth)
    this.th_min      = fates_unset_r8
    this.psi_min     = fates_unset_r8
    this.dpsidth_min = fates_unset_r8
    return nothing
end

function set_wkf_param!(this::wkf_type_smooth_cch, params_in::AbstractVector{<:Real})
    this.th_sat  = params_in[1]
    this.psi_sat = params_in[2]
    this.beta    = params_in[3]
    styp = Int(params_in[4])

    this.scch_ps, this.scch_pu, this.scch_b2, this.scch_b3 =
        _set_smooth_cch_coeffs(this.psi_sat, this.beta, styp)
    return nothing
end

get_thsat(this::wrf_type_smooth_cch) = this.th_sat

function th_from_psi(this::wrf_type_smooth_cch, psi::Real)
    alpha  = -1.0 / this.psi_sat
    lambda = 1.0 / this.beta
    pc = psi

    if pc <= this.scch_pu
        # Unsaturated full Brooks-Corey regime (pc <= pu < 0).
        sat = (-alpha * pc)^(-lambda)
    elseif pc < this.scch_ps
        # Cubic smoothing regime (pu < pc < ps <= 0).
        deltaPc = pc - this.scch_ps
        sat = 1.0 + deltaPc * deltaPc * (this.scch_b2 + deltaPc * this.scch_b3)
    else
        # Saturated regime (pc >= ps).
        sat = 1.0
    end
    th = sat * this.th_sat
    return th
end

function psi_from_th(this::wrf_type_smooth_cch, th::Real)
    relTol = 1.0e-9
    alpha  = -1.0 / this.psi_sat
    lambda = 1.0 / this.beta

    sat = max(1.0e-6, th / this.th_sat)
    if sat < 1.0
        # Find pc satisfying the unmodified Brooks-Corey function.
        Se = sat
        pc = -(Se^(-1.0 / lambda)) / alpha
        if pc > this.scch_pu
            # Solution is in the cubic smoothing regime.
            if this.scch_b2 == 0.0
                # b3 > 0.
                pc = this.scch_ps - ((1.0 - Se) / this.scch_b3)^(1.0 / 3.0)
            elseif this.scch_b3 == 0.0
                # b2 < 0.
                pc = this.scch_ps - sqrt((Se - 1.0) / this.scch_b2)
            else
                # General cubic 1 + b2*x^2 + b3*x^3 = Se, x = pc - ps.
                # Bracketed Newton-Raphson on residual r = x^2*(b2 + b3*x) + (1 - Se).
                xL = this.scch_pu - this.scch_ps
                xR = 0.0
                xc = pc - this.scch_ps

                iter = 0
                dx = 1.0e20  # something large
                while abs(dx) >= -relTol * this.scch_pu
                    iter += 1

                    # Reset xc using bisection if necessary.
                    if xc <= xL || xc >= xR
                        xc = xL + 0.5 * (xR - xL)
                    end

                    # Find NR step.
                    dx = this.scch_b3 * xc
                    resid = xc * xc * (this.scch_b2 + dx) + 1.0 - Se
                    dx = resid / (xc * (2.0 * this.scch_b2 + 3.0 * dx))

                    # Update bracket.
                    if resid > 0.0
                        xR = xc
                    else
                        xL = xc
                    end

                    # Take the Newton-Raphson step.
                    xc = xc - dx

                    if iter > 10000
                        fates_endrun("psi_from_th_smooth_cch iteration not converging")
                    end
                end
                # Here xc = pc - ps.
                pc = xc + this.scch_ps
            end
        end
    else
        pc = 0.0
    end
    return pc
end

function dpsidth_from_th(this::wrf_type_smooth_cch, th::Real)
    sat_res = 0.0
    alpha   = -1.0 / this.psi_sat
    lambda  = 1.0 / this.beta

    pc = 1.0 * psi_from_th(this, th)
    if pc <= this.scch_pu
        # Unsaturated full Brooks-Corey regime (pc <= pu < 0).
        Se = (-alpha * pc)^(-lambda)
        dSe_dpc = -lambda * Se / pc
        dsat_dp = (1.0 - sat_res) * dSe_dpc
        dpsidth = 1.0 / (dsat_dp * this.th_sat)
    elseif pc < this.scch_ps
        # Cubic smoothing regime (pu < pc < ps <= 0).
        deltaPc = pc - this.scch_ps
        dSe_dpc = deltaPc * (2 * this.scch_b2 + 3 * deltaPc * this.scch_b3)
        dsat_dp = (1.0 - sat_res) * dSe_dpc
        dpsidth = 1.0 / (dsat_dp * this.th_sat)
    else
        # Saturated regime (pc >= ps).
        dpsidth = this.dpsidth_max
    end
    return dpsidth
end

function ftc_from_psi(this::wkf_type_smooth_cch, psi::Real)
    psi_min = this.wrf.psi_min
    pc = psi
    alpha  = -1.0 / this.psi_sat
    lambda = 1.0 / this.beta

    if pc <= this.scch_pu
        Se = (-alpha * pc)^(-lambda)
        kr = Se^(3.0 + 2.0 / lambda)
    elseif pc < this.scch_ps
        deltaPc = pc - this.scch_ps
        Se = 1.0 + deltaPc * deltaPc * (this.scch_b2 + deltaPc * this.scch_b3)
        kr = Se^(3.0 + 2.0 / lambda)
    else
        kr = 1.0
    end

    if kr <= wtf_min_ftc
        ftc = wtf_min_ftc
    else
        min_ftc_weight, _ = get_min_ftc_weight(psi_min, Float64(psi))
        ftc = kr * (1.0 - min_ftc_weight) + wtf_min_ftc * min_ftc_weight
    end
    return ftc
end

function dftcdpsi_from_psi(this::wkf_type_smooth_cch, psi::Real)
    psi_min = this.wrf.psi_min
    ftc = ftc_from_psi(this, psi)

    if abs(ftc - wtf_min_ftc) < nearzero
        dftcdpsi = 0.0
    else
        pc = psi
        alpha  = -1.0 / this.psi_sat
        lambda = 1.0 / this.beta

        if pc <= this.scch_pu
            Se = (-alpha * pc)^(-lambda)
            dSe_dpc = -lambda * Se / pc
            kr = Se^(3.0 + 2.0 / lambda)
            dkr_dSe = (3.0 + 2.0 / lambda) * kr / Se
            dkr_dp  = dkr_dSe * dSe_dpc
        elseif pc < this.scch_ps
            deltaPc = pc - this.scch_ps
            Se = 1.0 + deltaPc * deltaPc * (this.scch_b2 + deltaPc * this.scch_b3)
            dSe_dpc = deltaPc * (2 * this.scch_b2 + 3 * deltaPc * this.scch_b3)
            kr = Se^(2.5 + 2.0 / lambda)
            dkr_dSe = (2.5 + 2.0 / lambda) * kr / Se
            dkr_dp  = dkr_dSe * dSe_dpc
        else
            kr = 1.0
            dkr_dp = 0.0
        end
        dftcdpsi = dkr_dp

        min_ftc_weight, dmin_ftc_weight_dpsi = get_min_ftc_weight(psi_min, Float64(psi))

        dftcdpsi = dftcdpsi -
                   (dftcdpsi * min_ftc_weight + ftc * dmin_ftc_weight_dpsi) +
                   wtf_min_ftc * dmin_ftc_weight_dpsi
    end
    return dftcdpsi
end

# ===========================================================================
# TFS functions
# ===========================================================================

"TFS water retention function."
Base.@kwdef mutable struct wrf_type_tfs <: WRFType
    psi_max::Float64     = 0.0
    psi_min::Float64     = 0.0
    dpsidth_max::Float64 = 0.0
    dpsidth_min::Float64 = 0.0
    th_min::Float64      = 0.0
    th_max::Float64      = 0.0
    th_sat::Float64   = 0.0   # Saturation volumetric water content [m3/m3]
    th_res::Float64   = 0.0   # Residual volumetric water content   [m3/m3]
    pinot::Float64    = 0.0   # osmotic potential at full turgor [MPa]
    epsil::Float64    = 0.0   # bulk elastic modulus [MPa]
    rwc_ft::Float64   = 0.0   # RWC @ full turgor [-]
    cap_corr::Float64 = 0.0   # correction for nonzero psi0x
    cap_int::Float64  = 0.0   # intercept of capillary region
    cap_slp::Float64  = 0.0   # slope of capillary region
    pmedia::Int       = 0     # porous media index
end

"TFS water conductivity function."
Base.@kwdef mutable struct wkf_type_tfs <: WKFType
    wrf::Union{WRFType,Nothing} = nothing
    p50::Float64    = 0.0   # matric potential at 50% conductivity loss [MPa]
    avuln::Float64  = 0.0   # vulnerability curve parameter
    th_sat::Float64 = 0.0   # volumetric water content at saturation
end

function set_wkf_param!(this::wkf_type_tfs, params_in::AbstractVector{<:Real})
    this.p50   = params_in[1]
    this.avuln = params_in[2]
    return nothing
end

function set_wrf_param!(this::wrf_type_tfs, params_in::AbstractVector{<:Real})
    this.th_sat   = params_in[1]
    this.th_res   = params_in[2]
    this.pinot    = params_in[3]
    this.epsil    = params_in[4]
    this.rwc_ft   = params_in[5]
    this.cap_corr = params_in[6]
    this.cap_int  = params_in[7]
    this.cap_slp  = params_in[8]
    this.pmedia   = Int(params_in[9])
    set_min_max_from_satres!(this, this.th_res, this.th_sat)
    return nothing
end

get_thsat(this::wrf_type_tfs) = this.th_sat

# --- TFS PV-curve component functions (Fortran module subroutines) ---

"Solute water potential (negative) vs water content for the plant PV curve."
solutepsi(th, rwc_ft, th_sat, th_res, pinot) =
    pinot * (th_sat * rwc_ft - th_res) / (th - th_res)

"Derivative of solutepsi() wrt theta."
dsolutepsidth(th, th_sat, th_res, rwc_ft, pinot) =
    -1.0 * pinot * (th_sat * rwc_ft - th_res) * (th - th_res)^(-2.0)

"Pressure water potential (positive) vs water content for the plant PV curve."
pressurepsi(th, rwc_ft, th_sat, th_res, pinot, epsil) =
    epsil * (th - th_sat * rwc_ft) / (th_sat * rwc_ft - th_res) - pinot

"Derivative of pressurepsi() wrt theta."
dpressurepsidth(th_sat, th_res, rwc_ft, epsil) =
    epsil / (th_sat * rwc_ft - th_res)

"Water potential in the capillary region of the plant PV curve (sapwood only)."
capillarypsi(th, th_sat, cap_int, cap_slp) =
    cap_int + th * cap_slp / th_sat

"Derivative of capillarypsi() wrt theta."
dcapillarypsidth(cap_slp, th_sat) = cap_slp / th_sat

"""
    bisect_pv(this, lower, upper, psi) -> th

Bisection inverse of the plant PV curve (no analytical inverse exists due to the
quadratic smoothing). Returns th for the given psi.
"""
function bisect_pv(this::wrf_type_tfs, lower::Float64, upper::Float64, psi::Float64)
    xtol = 1.0e-16   # error tolerance for th [m3/m3]
    ytol = 1.0e-8    # error tolerance for psi [MPa]

    if psi > 0.0
        fates_endrun("Error: psi became positive during pv bisection; psi: " * string(psi))
    end

    y_lo = psi_from_th(this, lower)
    y_hi = psi_from_th(this, upper)

    f_lo = y_lo - psi
    f_hi = y_hi - psi
    chg  = upper - lower

    x_new = 0.5 * (lower + upper)
    nitr = 0
    while abs(chg) > xtol && nitr < 100
        x_new = 0.5 * (lower + upper)
        y_new = psi_from_th(this, x_new)
        f_new = y_new - psi
        if abs(f_new) <= ytol
            break
        end
        if (f_lo * f_new) < 0.0
            upper = x_new
        end
        if (f_hi * f_new) < 0.0
            lower = x_new
        end
        chg = upper - lower
        nitr += 1
    end

    if nitr == 100
        FatesWarn("Warning: number of iteration reaches 100 for bisect_pv")
    end
    return x_new
end

function th_from_psi(this::wrf_type_tfs, psi::Real)
    if psi > this.psi_max
        th = th_linear_sat(this, psi)             # Linear range for extreme values
    elseif psi < this.psi_min
        th = th_linear_res(this, psi)             # Linear range for extreme values
    else
        # Bisection search; define bounds.
        lower = this.th_min - 1.0e-9
        upper = this.th_max + 1.0e-9
        th = bisect_pv(this, lower, upper, Float64(psi))
        psi_check = psi_from_th(this, th)
        if psi_check > -1.0e-8
            fates_endrun("bisect_pv returned positive value for water potential?")
        end
    end
    return th
end

function psi_from_th(this::wrf_type_tfs, th::Real)
    if th > this.th_max
        psi = psi_linear_sat(this, th)
    elseif th < this.th_min
        psi = psi_linear_res(this, th)
    else
        th_corr = th * this.cap_corr

        # Two rounds of quadratic smoothing: (1) elastic+capillary, then (2) with cavitation.
        psi_sol   = solutepsi(th_corr, this.rwc_ft, this.th_sat, this.th_res, this.pinot)
        psi_press = pressurepsi(th_corr, this.rwc_ft, this.th_sat, this.th_res, this.pinot, this.epsil)
        psi_elastic = psi_sol + psi_press

        if this.pmedia == 1
            # leaves have no capillary region
            psi_capelast = psi_elastic
        elseif this.pmedia <= 4
            # sapwood has a capillary region
            psi_capillary = capillarypsi(th_corr, this.th_sat, this.cap_int, this.cap_slp)
            b = -1.0 * (psi_capillary + psi_elastic)
            c = psi_capillary * psi_elastic
            psi_capelast = (-b - sqrt(b * b - 4.0 * quad_a1 * c)) / (2.0 * quad_a1)
        else
            fates_endrun("TFS WRF was called for an inelligable porous media")
        end

        # Smooth the result of capillary-elastic with cavitation.
        psi_cavitation = psi_sol
        b = -1.0 * (psi_capelast + psi_cavitation)
        c = psi_capelast * psi_cavitation
        psi = (-b + sqrt(b * b - 4.0 * quad_a2 * c)) / (2.0 * quad_a2)
    end
    return psi
end

function dpsidth_from_th(this::wrf_type_tfs, th::Real)
    if th > this.th_max
        dpsidth = this.dpsidth_max
    elseif th < this.th_min
        dpsidth = this.dpsidth_min
    else
        th_corr = th * this.cap_corr

        psi_sol   = solutepsi(th_corr, this.rwc_ft, this.th_sat, this.th_res, this.pinot)
        psi_press = pressurepsi(th_corr, this.rwc_ft, this.th_sat, this.th_res, this.pinot, this.epsil)

        dsol_dth   = dsolutepsidth(th, this.th_sat, this.th_res, this.rwc_ft, this.pinot)
        dpress_dth = dpressurepsidth(this.th_sat, this.th_res, this.rwc_ft, this.epsil)

        delast_dth  = dsol_dth + dpress_dth
        psi_elastic = psi_sol + psi_press

        if this.pmedia == 1
            psi_capelast  = psi_elastic
            dcapelast_dth = delast_dth
        elseif this.pmedia <= 4
            psi_capillary = capillarypsi(th, this.th_sat, this.cap_int, this.cap_slp)
            b = -1.0 * (psi_capillary + psi_elastic)
            c = psi_capillary * psi_elastic
            psi_capelast = (-b - sqrt(b * b - 4.0 * quad_a1 * c)) / (2.0 * quad_a1)

            dcap_dth = dcapillarypsidth(this.cap_slp, this.th_sat)

            dbdth = -1.0 * (delast_dth + dcap_dth)
            dcdth = psi_elastic * dcap_dth + delast_dth * psi_capillary

            dcapelast_dth = 1.0 / (2.0 * quad_a1) *
                            (-dbdth - 0.5 * ((b * b - 4.0 * quad_a1 * c)^(-0.5)) *
                             (2.0 * b * dbdth - 4.0 * quad_a1 * dcdth))
        else
            fates_endrun("TFS WRF was called for an ineligible porous media")
        end

        # Smooth capillary-elastic with cavitation.
        psi_cavitation = psi_sol
        b = -1.0 * (psi_capelast + psi_cavitation)
        c = psi_capelast * psi_cavitation

        dcav_dth = dsol_dth
        dbdth = -1.0 * (dcapelast_dth + dcav_dth)
        dcdth = psi_capelast * dcav_dth + dcapelast_dth * psi_cavitation

        dpsidth = 1.0 / (2.0 * quad_a2) *
                  (-dbdth + 0.5 * ((b * b - 4.0 * quad_a2 * c)^(-0.5)) *
                   (2.0 * b * dbdth - 4.0 * quad_a2 * dcdth))
    end
    return dpsidth
end

function ftc_from_psi(this::wkf_type_tfs, psi::Real)
    psi_min = this.wrf.psi_min
    psi_eff = max(psi_min, min(-nearzero, psi))

    ftc = 1.0 / (1.0 + (psi_eff / this.p50)^this.avuln)

    if ftc <= wtf_min_ftc
        ftc = wtf_min_ftc
    else
        # Add protections; no conductance at incredibly low suction.
        min_ftc_weight, _ = get_min_ftc_weight(psi_min, psi_eff)
        ftc = ftc * (1.0 - min_ftc_weight) + wtf_min_ftc * min_ftc_weight
    end
    return ftc
end

function dftcdpsi_from_psi(this::wkf_type_tfs, psi::Real)
    psi_min = this.wrf.psi_min
    psi_eff = max(psi_min, min(-nearzero, psi))

    if psi_eff > 0.0
        dftcdpsi = 0.0
    else
        ftc = ftc_from_psi(this, psi_eff)

        if abs(ftc - wtf_min_ftc) < nearzero
            dftcdpsi = 0.0   # We cap ftc, so derivative is zero
        else
            fx  = 1.0 + (psi_eff / this.p50)^this.avuln
            dfx = this.avuln * (psi_eff / this.p50)^(this.avuln - 1.0) * (1.0 / this.p50)
            dftcdpsi = -fx^(-2.0) * dfx

            min_ftc_weight, dmin_ftc_weight_dpsi = get_min_ftc_weight(psi_min, psi_eff)

            dftcdpsi = dftcdpsi -
                       (dftcdpsi * min_ftc_weight + ftc * dmin_ftc_weight_dpsi) +
                       wtf_min_ftc * dmin_ftc_weight_dpsi
        end
    end
    return dftcdpsi
end
