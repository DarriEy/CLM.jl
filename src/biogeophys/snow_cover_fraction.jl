# ==========================================================================
# Ported from:
#   src/biogeophys/SnowCoverFractionBaseMod.F90
#   src/biogeophys/SnowCoverFractionSwensonLawrence2012Mod.F90
#   src/biogeophys/SnowCoverFractionNiuYang2007Mod.F90
#   src/biogeophys/SnowCoverFractionFactoryMod.F90
#
# Snow cover fraction parameterizations for CLM.
#
# Provides:
#   - Abstract base type (SnowCoverFractionBase) and two concrete subtypes:
#       SnowCoverFractionSwensonLawrence2012
#       SnowCoverFractionNiuYang2007
#   - update_snow_depth_and_frac!   — Update snow depth and frac_sno
#   - add_newsnow_to_intsnow!       — Update integrated snowfall
#   - frac_snow_during_melt         — Point function for melt-period frac_sno
#   - calc_frac_sno_eff!            — Calculate effective snow-covered fraction
#   - create_snow_cover_fraction    — Factory constructor
# ==========================================================================

# =========================================================================
# Abstract base type (SnowCoverFractionBaseMod.F90)
# =========================================================================

"""
    SnowCoverFractionBase

Abstract base type for snow cover fraction parameterizations.

Ported from `snow_cover_fraction_base_type` in `SnowCoverFractionBaseMod.F90`.
"""
abstract type SnowCoverFractionBase end

# =========================================================================
# SwensonLawrence2012 type (SnowCoverFractionSwensonLawrence2012Mod.F90)
# =========================================================================

"""
    SnowCoverFractionSwensonLawrence2012

Snow cover fraction parameterization from Swenson & Lawrence (2012).
Uses changes in SWE to determine fractional snow covered area.

Fields:
- `int_snow_max` : limit on integrated snowfall during melt (mm H2O)
- `accum_factor` : accumulation constant for fractional snow covered area (unitless)
- `n_melt`       : SCA shape parameter per column (unitless)

Ported from `snow_cover_fraction_swenson_lawrence_2012_type` in
`SnowCoverFractionSwensonLawrence2012Mod.F90`.
"""
Base.@kwdef mutable struct SnowCoverFractionSwensonLawrence2012 <: SnowCoverFractionBase
    int_snow_max::Float64 = 2000.0        # limit on integrated snowfall during melt [mm H2O]
    accum_factor::Float64 = 0.1           # accumulation constant for fsca [unitless]
    n_melt::Vector{Float64} = Float64[]   # SCA shape parameter per column [unitless]
end

# =========================================================================
# NiuYang2007 type (SnowCoverFractionNiuYang2007Mod.F90)
# =========================================================================

"""
    SnowCoverFractionNiuYang2007

Snow cover fraction parameterization from Niu & Yang (2007).
Computes fractional snow cover from snow depth, density and surface roughness.

Fields:
- `zlnd` : roughness length for soil (m)

Ported from `snow_cover_fraction_niu_yang_2007_type` in
`SnowCoverFractionNiuYang2007Mod.F90`.
"""
Base.@kwdef mutable struct SnowCoverFractionNiuYang2007 <: SnowCoverFractionBase
    zlnd::Float64 = 0.01   # Roughness length for soil (m)
end

# =========================================================================
# Initialization
# =========================================================================

"""
    snow_cover_fraction_init!(scf::SnowCoverFractionSwensonLawrence2012,
        ncols; col_lun_itype, col_gridcell, col_topo_std,
        allow_multiple_columns_grc,
        int_snow_max, accum_factor, n_melt_coef, n_melt_glcmec)

Initialize the SwensonLawrence2012 parameterization for `ncols` columns.
Computes per-column `n_melt` from topographic standard deviation.

Ported from `Init`, `ReadParams`, `ReadNamelist`, `SetDerivedParameters`
in `SnowCoverFractionSwensonLawrence2012Mod.F90`.
"""
function snow_cover_fraction_init!(
    scf::SnowCoverFractionSwensonLawrence2012,
    ncols::Int;
    col_lun_itype::Vector{Int},
    col_gridcell::Vector{Int},
    col_topo_std::Vector{<:Real},
    allow_multiple_columns_grc::Vector{Bool} = Bool[],
    int_snow_max::Real = 2000.0,
    accum_factor::Real = 0.1,
    n_melt_coef::Real = 200.0,
    n_melt_glcmec::Real = 10.0
)
    # Validate
    @assert int_snow_max > 0.0 "int_snow_max must be > 0, got $int_snow_max"
    @assert n_melt_glcmec > 0.0 "n_melt_glcmec must be > 0, got $n_melt_glcmec"

    scf.int_snow_max = int_snow_max
    scf.accum_factor = accum_factor
    scf.n_melt = Vector{Float64}(undef, ncols)

    has_glc_info = !isempty(allow_multiple_columns_grc)

    for c in 1:ncols
        if has_glc_info
            g = col_gridcell[c]
            if col_lun_itype[c] == ISTICE && g >= 1 && g <= length(allow_multiple_columns_grc) &&
               allow_multiple_columns_grc[g]
                # Ice columns with multiple elevation classes: use fixed n_melt
                # to avoid double-counting topographic variability
                scf.n_melt[c] = n_melt_glcmec
                continue
            end
        end
        scf.n_melt[c] = n_melt_coef / max(10.0, col_topo_std[c])
    end

    return nothing
end

"""
    snow_cover_fraction_init!(scf::SnowCoverFractionNiuYang2007;
        zlnd, use_subgrid_fluxes)

Initialize the NiuYang2007 parameterization.

Ported from `Init` and `ReadParams` in `SnowCoverFractionNiuYang2007Mod.F90`.
"""
function snow_cover_fraction_init!(
    scf::SnowCoverFractionNiuYang2007;
    zlnd::Real = 0.01,
    use_subgrid_fluxes::Bool = false
)
    if use_subgrid_fluxes
        error("NiuYang2007 snow cover fraction parameterization is incompatible " *
              "with use_subgrid_fluxes = true")
    end
    scf.zlnd = zlnd
    return nothing
end

# =========================================================================
# CalcFracSnoEff (SnowCoverFractionBaseMod.F90)
# =========================================================================

"""
    calc_frac_sno_eff!(frac_sno_eff, frac_sno, lun_itype_col, urbpoi,
        mask, bounds; use_subgrid_fluxes)

Calculate effective snow-covered fraction given frac_sno.

For urban or lake columns, or when use_subgrid_fluxes is false,
frac_sno_eff is forced to 0 or 1. Otherwise it equals frac_sno.

Ported from `CalcFracSnoEff` in `SnowCoverFractionBaseMod.F90`.
"""
# Per-column kernel: effective snow-covered fraction. `use_subgrid` carried as a
# Bool; ISTDLAK as an Int literal. Byte-identical on CPU Float64.
@kernel function _calc_frac_sno_eff_kernel!(frac_sno_eff, @Const(frac_sno),
        @Const(lun_itype_col), @Const(urbpoi), @Const(mask), use_subgrid::Bool,
        istdlak::Int, cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && mask[c]
        allow_fractional = !(urbpoi[c] || lun_itype_col[c] == istdlak || !use_subgrid)
        if allow_fractional
            frac_sno_eff[c] = frac_sno[c]
        else
            frac_sno_eff[c] = frac_sno[c] > zero(eltype(frac_sno)) ?
                one(eltype(frac_sno_eff)) : zero(eltype(frac_sno_eff))
        end
    end
end

function calc_frac_sno_eff!(
    frac_sno_eff::AbstractVector{<:Real},
    frac_sno::AbstractVector{<:Real},
    lun_itype_col::AbstractVector{<:Integer},
    urbpoi::AbstractVector{Bool},
    mask::AbstractVector{Bool},
    bounds::UnitRange{Int};
    use_subgrid_fluxes::Bool = true
)
    _launch!(_calc_frac_sno_eff_kernel!, frac_sno_eff, frac_sno, lun_itype_col,
             urbpoi, mask, use_subgrid_fluxes, ISTDLAK, first(bounds), last(bounds))
    return nothing
end

# =========================================================================
# FracSnowDuringMelt — SwensonLawrence2012 (science function)
# =========================================================================

"""
    frac_snow_during_melt(scf::SnowCoverFractionSwensonLawrence2012,
        c, h2osno_total, int_snow) -> Float64

Single-point function returning fractional snow cover during melt for the
SwensonLawrence2012 parameterization.

Ported from `FracSnowDuringMelt` in
`SnowCoverFractionSwensonLawrence2012Mod.F90`.
"""
function frac_snow_during_melt(
    scf::SnowCoverFractionSwensonLawrence2012,
    c::Int,
    h2osno_total::Real,
    int_snow::Real
)
    int_snow_limited = min(int_snow, scf.int_snow_max)
    smr = min(1.0, h2osno_total / int_snow_limited)
    frac_sno = 1.0 - (acos(min(1.0, 2.0 * smr - 1.0)) / RPI) ^ scf.n_melt[c]
    return frac_sno
end

"""
    frac_snow_during_melt(scf::SnowCoverFractionNiuYang2007,
        c, h2osno_total, int_snow) -> Float64

Single-point function for fractional snow cover during melt.
Not implemented for NiuYang2007; returns NaN.

Ported from `FracSnowDuringMelt` in `SnowCoverFractionNiuYang2007Mod.F90`.
"""
function frac_snow_during_melt(
    scf::SnowCoverFractionNiuYang2007,
    c::Int,
    h2osno_total::Real,
    int_snow::Real
)
    return NaN
end

# =========================================================================
# UpdateSnowDepthAndFrac — SwensonLawrence2012
# =========================================================================

"""
    update_snow_depth_and_frac!(scf::SnowCoverFractionSwensonLawrence2012,
        frac_sno, frac_sno_eff, snow_depth,
        lun_itype_col, urbpoi, h2osno_total, snowmelt, int_snow,
        newsnow, bifall, mask, bounds;
        use_subgrid_fluxes)

Update snow depth and snow fraction using the SwensonLawrence2012 parameterization.

FSCA is based on *changes* in SWE:
- Accumulation: frac_sno updated via tanh(accum_factor * newsnow)
- Melt: frac_sno via the hysteresis depletion curve

Ported from `UpdateSnowDepthAndFrac` in
`SnowCoverFractionSwensonLawrence2012Mod.F90`.
"""
# Per-column kernel: SwensonLawrence2012 frac_sno update. frac_snow_during_melt is
# folded inline (uses acos, RPI, n_melt[c]). Scalar params (accum_factor, int_snow_max,
# RPI) carried at the working precision; n_melt is a device-resident vector. Byte-
# identical on CPU Float64.
@kernel function _sl12_update_fracsno_kernel!(frac_sno, @Const(h2osno_total),
        @Const(snowmelt), @Const(int_snow), @Const(newsnow), @Const(n_melt),
        @Const(mask), accum_factor, int_snow_max, rpi, cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && mask[c]
        T = eltype(frac_sno)
        if h2osno_total[c] == zero(T)
            if newsnow[c] > zero(T)
                frac_sno[c] = tanh(accum_factor * newsnow[c])
            else
                frac_sno[c] = zero(T)
            end
        else  # h2osno_total > 0
            if snowmelt[c] > zero(T)
                int_snow_limited = min(int_snow[c], int_snow_max)
                smr = min(one(T), h2osno_total[c] / int_snow_limited)
                frac_sno[c] = one(T) -
                    (acos(min(one(T), T(2.0) * smr - one(T))) / rpi) ^ n_melt[c]
            end
            if newsnow[c] > zero(T)
                frac_sno[c] = frac_sno[c] +
                    tanh(accum_factor * newsnow[c]) * (one(T) - frac_sno[c])
            end
        end
    end
end

# Per-column kernel: SwensonLawrence2012 snow_depth update.
@kernel function _sl12_update_snowdepth_kernel!(snow_depth, @Const(h2osno_total),
        @Const(newsnow), @Const(bifall), @Const(frac_sno_eff), @Const(mask),
        cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && mask[c]
        T = eltype(snow_depth)
        if h2osno_total[c] > zero(T)
            if frac_sno_eff[c] > zero(T)
                snow_depth[c] = snow_depth[c] + newsnow[c] / (bifall[c] * frac_sno_eff[c])
            else
                snow_depth[c] = zero(T)
            end
        else  # h2osno_total == 0
            if newsnow[c] > zero(T)
                z_avg = newsnow[c] / bifall[c]
                snow_depth[c] = z_avg / frac_sno_eff[c]
            else
                snow_depth[c] = zero(T)
            end
        end
    end
end

function update_snow_depth_and_frac!(
    scf::SnowCoverFractionSwensonLawrence2012,
    # Outputs (modified in place)
    frac_sno::AbstractVector{<:Real},
    frac_sno_eff::AbstractVector{<:Real},
    snow_depth::AbstractVector{<:Real},
    # Inputs
    lun_itype_col::AbstractVector{<:Integer},
    urbpoi::AbstractVector{Bool},
    h2osno_total::AbstractVector{<:Real},
    snowmelt::AbstractVector{<:Real},
    int_snow::AbstractVector{<:Real},
    newsnow::AbstractVector{<:Real},
    bifall::AbstractVector{<:Real},
    mask::AbstractVector{Bool},
    bounds::UnitRange{Int};
    use_subgrid_fluxes::Bool = true,
    n_melt_dev::AbstractVector{<:Real} = scf.n_melt
)
    T = eltype(frac_sno)
    # Auto-initialize n_melt if not yet set (default: n_melt_coef/10 = 20).
    # Only applies to the host-resident default vector.
    ncols = length(frac_sno)
    if n_melt_dev === scf.n_melt && length(scf.n_melt) < ncols
        resize!(scf.n_melt, ncols)
        fill!(scf.n_melt, 200.0 / 10.0)
        n_melt_dev = scf.n_melt
    end

    # scf is a host-resident config struct; copy its n_melt onto the working
    # backend so the kernels don't get a host Array among device args. No-op on CPU.
    if !(frac_sno isa Array) && n_melt_dev isa Array
        nm = similar(frac_sno, T, length(n_melt_dev))
        copyto!(nm, T.(n_melt_dev))
        n_melt_dev = nm
    end

    # ---- Update frac_sno (folds frac_snow_during_melt inline) ----
    _launch!(_sl12_update_fracsno_kernel!, frac_sno, h2osno_total, snowmelt,
             int_snow, newsnow, n_melt_dev, mask,
             T(scf.accum_factor), T(scf.int_snow_max), T(RPI),
             first(bounds), last(bounds))

    # ---- Calculate frac_sno_eff ----
    calc_frac_sno_eff!(frac_sno_eff, frac_sno, lun_itype_col, urbpoi,
        mask, bounds; use_subgrid_fluxes=use_subgrid_fluxes)

    # ---- Update snow_depth ----
    _launch!(_sl12_update_snowdepth_kernel!, snow_depth, h2osno_total, newsnow,
             bifall, frac_sno_eff, mask, first(bounds), last(bounds))

    return nothing
end

# =========================================================================
# UpdateSnowDepthAndFrac — NiuYang2007
# =========================================================================

"""
    update_snow_depth_and_frac!(scf::SnowCoverFractionNiuYang2007,
        frac_sno, frac_sno_eff, snow_depth,
        lun_itype_col, urbpoi, h2osno_total, snowmelt, int_snow,
        newsnow, bifall, mask, bounds;
        use_subgrid_fluxes)

Update snow depth and snow fraction using the NiuYang2007 parameterization.

Snow cover fraction is computed from snow depth, snow density and surface
roughness length.

Ported from `UpdateSnowDepthAndFrac` in
`SnowCoverFractionNiuYang2007Mod.F90`.
"""
# Per-column kernel: NiuYang2007 snow_depth + frac_sno update. zlnd carried at the
# working precision; literals eltype-converted. Byte-identical on CPU Float64.
@kernel function _ny07_update_depthfrac_kernel!(snow_depth, frac_sno,
        @Const(h2osno_total), @Const(newsnow), @Const(bifall), @Const(mask),
        zlnd, cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && mask[c]
        T = eltype(snow_depth)
        if h2osno_total[c] == zero(T)
            snow_depth[c] = zero(T)
        end
        snow_depth[c] = snow_depth[c] + newsnow[c] / bifall[c]
        if snow_depth[c] > zero(T)
            frac_sno[c] = tanh(snow_depth[c] / (T(2.5) * zlnd *
                (min(T(800.0), (h2osno_total[c] + newsnow[c]) / snow_depth[c]) / T(100.0))^one(T)))
        else
            frac_sno[c] = zero(T)
        end
        if h2osno_total[c] > zero(T) && h2osno_total[c] < one(T)
            frac_sno[c] = min(frac_sno[c], h2osno_total[c])
        end
    end
end

function update_snow_depth_and_frac!(
    scf::SnowCoverFractionNiuYang2007,
    # Outputs (modified in place)
    frac_sno::AbstractVector{<:Real},
    frac_sno_eff::AbstractVector{<:Real},
    snow_depth::AbstractVector{<:Real},
    # Inputs
    lun_itype_col::AbstractVector{<:Integer},
    urbpoi::AbstractVector{Bool},
    h2osno_total::AbstractVector{<:Real},
    snowmelt::AbstractVector{<:Real},
    int_snow::AbstractVector{<:Real},
    newsnow::AbstractVector{<:Real},
    bifall::AbstractVector{<:Real},
    mask::AbstractVector{Bool},
    bounds::UnitRange{Int};
    use_subgrid_fluxes::Bool = false
)
    T = eltype(frac_sno)
    _launch!(_ny07_update_depthfrac_kernel!, snow_depth, frac_sno, h2osno_total,
             newsnow, bifall, mask, T(scf.zlnd), first(bounds), last(bounds))

    # Calculate frac_sno_eff
    calc_frac_sno_eff!(frac_sno_eff, frac_sno, lun_itype_col, urbpoi,
        mask, bounds; use_subgrid_fluxes=use_subgrid_fluxes)

    return nothing
end

# =========================================================================
# AddNewsnowToIntsnow — SwensonLawrence2012
# =========================================================================

"""
    add_newsnow_to_intsnow!(scf::SnowCoverFractionSwensonLawrence2012,
        int_snow, newsnow, h2osno_total, frac_sno,
        mask, bounds)

Add new snow to integrated snowfall for the SwensonLawrence2012 parameterization.
Resets int_snow after accumulation events to maintain consistency with
frac_sno and h2osno_total.

Ported from `AddNewsnowToIntsnow` in
`SnowCoverFractionSwensonLawrence2012Mod.F90`.
"""
# Per-column kernel: SwensonLawrence2012 int_snow accumulation. smooth_max on a
# non-AD eltype (Float32/Float64) is exact `max` (GPU-safe). n_melt is a device
# vector; RPI and the limit constants carried at the working precision. Byte-
# identical on CPU Float64.
@kernel function _sl12_add_intsnow_kernel!(int_snow, @Const(newsnow),
        @Const(h2osno_total), @Const(frac_sno), @Const(n_melt), @Const(mask),
        rpi, cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && mask[c]
        T = eltype(int_snow)
        if newsnow[c] > zero(T)
            fsno_base = smooth_max(one(T) - max(frac_sno[c], T(1.0e-6)), T(1.0e-10))
            temp_intsnow = (h2osno_total[c] + newsnow[c]) /
                (T(0.5) * (cos(rpi * fsno_base^(one(T) / n_melt[c])) + one(T)))
            int_snow[c] = min(T(1.0e8), temp_intsnow)
        end
        # NOTE(wjs, 2019-07-25): extra addition of new_snow kept for bit-for-bit
        # consistency with the Fortran (see clm45 branch history).
        int_snow[c] = int_snow[c] + newsnow[c]
    end
end

function add_newsnow_to_intsnow!(
    scf::SnowCoverFractionSwensonLawrence2012,
    int_snow::AbstractVector{<:Real},
    newsnow::AbstractVector{<:Real},
    h2osno_total::AbstractVector{<:Real},
    frac_sno::AbstractVector{<:Real},
    mask::AbstractVector{Bool},
    bounds::UnitRange{Int};
    n_melt_dev::AbstractVector{<:Real} = scf.n_melt
)
    T = eltype(int_snow)
    # Auto-initialize n_melt if not yet set (host-resident default vector only).
    ncols = length(int_snow)
    if n_melt_dev === scf.n_melt && length(scf.n_melt) < ncols
        resize!(scf.n_melt, ncols)
        fill!(scf.n_melt, 200.0 / 10.0)
        n_melt_dev = scf.n_melt
    end
    # scf is host-resident; move n_melt onto the working backend (no-op on CPU).
    if !(int_snow isa Array) && n_melt_dev isa Array
        nm = similar(int_snow, T, length(n_melt_dev))
        copyto!(nm, T.(n_melt_dev))
        n_melt_dev = nm
    end

    _launch!(_sl12_add_intsnow_kernel!, int_snow, newsnow, h2osno_total, frac_sno,
             n_melt_dev, mask, T(RPI), first(bounds), last(bounds))
    return nothing
end

# =========================================================================
# AddNewsnowToIntsnow — NiuYang2007
# =========================================================================

"""
    add_newsnow_to_intsnow!(scf::SnowCoverFractionNiuYang2007,
        int_snow, newsnow, h2osno_total, frac_sno,
        mask, bounds)

Add new snow to integrated snowfall for the NiuYang2007 parameterization.
Straightforward addition.

Ported from `AddNewsnowToIntsnow` in `SnowCoverFractionNiuYang2007Mod.F90`.
"""
# Per-column kernel: NiuYang2007 int_snow accumulation (straight add).
@kernel function _ny07_add_intsnow_kernel!(int_snow, @Const(newsnow), @Const(mask),
        cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && mask[c]
        int_snow[c] = int_snow[c] + newsnow[c]
    end
end

function add_newsnow_to_intsnow!(
    scf::SnowCoverFractionNiuYang2007,
    int_snow::AbstractVector{<:Real},
    newsnow::AbstractVector{<:Real},
    h2osno_total::AbstractVector{<:Real},
    frac_sno::AbstractVector{<:Real},
    mask::AbstractVector{Bool},
    bounds::UnitRange{Int};
    n_melt_dev::AbstractVector{<:Real} = newsnow
)
    _launch!(_ny07_add_intsnow_kernel!, int_snow, newsnow, mask,
             first(bounds), last(bounds))
    return nothing
end

# =========================================================================
# Factory (SnowCoverFractionFactoryMod.F90)
# =========================================================================

"""
    create_snow_cover_fraction(method::String) -> SnowCoverFractionBase

Factory function to create a snow cover fraction parameterization instance.
Returns an uninitialized instance; call `snow_cover_fraction_init!` to set
parameters.

Supported methods:
- `"SwensonLawrence2012"` (default in CLM5)
- `"NiuYang2007"`

Ported from `CreateAndInitSnowCoverFraction` in
`SnowCoverFractionFactoryMod.F90`.
"""
function create_snow_cover_fraction(method::String = "SwensonLawrence2012")
    if method == "SwensonLawrence2012"
        return SnowCoverFractionSwensonLawrence2012()
    elseif method == "NiuYang2007"
        return SnowCoverFractionNiuYang2007()
    else
        error("Unknown snow_cover_fraction_method: '$method'. " *
              "Supported: \"SwensonLawrence2012\", \"NiuYang2007\"")
    end
end
