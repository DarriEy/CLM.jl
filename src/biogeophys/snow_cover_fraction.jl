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
    col_topo_std::Vector{Float64},
    allow_multiple_columns_grc::Vector{Bool} = Bool[],
    int_snow_max::Float64 = 2000.0,
    accum_factor::Float64 = 0.1,
    n_melt_coef::Float64 = 200.0,
    n_melt_glcmec::Float64 = 10.0
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
    zlnd::Float64 = 0.01,
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
function calc_frac_sno_eff!(
    frac_sno_eff::Vector{Float64},
    frac_sno::Vector{Float64},
    lun_itype_col::Vector{Int},
    urbpoi::Vector{Bool},
    mask::BitVector,
    bounds::UnitRange{Int};
    use_subgrid_fluxes::Bool = true
)
    for c in bounds
        mask[c] || continue

        # Determine whether fractional frac_sno_eff is allowed
        allow_fractional = !(urbpoi[c] || lun_itype_col[c] == ISTDLAK || !use_subgrid_fluxes)

        if allow_fractional
            frac_sno_eff[c] = frac_sno[c]
        else
            frac_sno_eff[c] = frac_sno[c] > 0.0 ? 1.0 : 0.0
        end
    end

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
    h2osno_total::Float64,
    int_snow::Float64
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
    h2osno_total::Float64,
    int_snow::Float64
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
function update_snow_depth_and_frac!(
    scf::SnowCoverFractionSwensonLawrence2012,
    # Outputs (modified in place)
    frac_sno::Vector{Float64},
    frac_sno_eff::Vector{Float64},
    snow_depth::Vector{Float64},
    # Inputs
    lun_itype_col::Vector{Int},
    urbpoi::Vector{Bool},
    h2osno_total::Vector{Float64},
    snowmelt::Vector{Float64},
    int_snow::Vector{Float64},
    newsnow::Vector{Float64},
    bifall::Vector{Float64},
    mask::BitVector,
    bounds::UnitRange{Int};
    use_subgrid_fluxes::Bool = true
)
    # ---- Update frac_sno ----
    for c in bounds
        mask[c] || continue

        # FSCA parameterization based on *changes* in SWE
        if h2osno_total[c] == 0.0
            if newsnow[c] > 0.0
                frac_sno[c] = tanh(scf.accum_factor * newsnow[c])
            else
                # Reset frac_sno when no snow
                frac_sno[c] = 0.0
            end
        else  # h2osno_total > 0
            if snowmelt[c] > 0.0
                # Compute change from melt during previous time step
                frac_sno[c] = frac_snow_during_melt(scf,
                    c, h2osno_total[c], int_snow[c])
            end

            if newsnow[c] > 0.0
                # Update fsca by new snow event, add to previous fsca
                # Algebraically equivalent to Swenson & Lawrence 2012 eqn. 3,
                # but simpler and less prone to roundoff errors
                # (see https://github.com/ESCOMP/ctsm/issues/784)
                frac_sno[c] = frac_sno[c] + tanh(scf.accum_factor * newsnow[c]) * (1.0 - frac_sno[c])
            end
        end
    end

    # ---- Calculate frac_sno_eff ----
    calc_frac_sno_eff!(frac_sno_eff, frac_sno, lun_itype_col, urbpoi,
        mask, bounds; use_subgrid_fluxes=use_subgrid_fluxes)

    # ---- Update snow_depth ----
    for c in bounds
        mask[c] || continue

        if h2osno_total[c] > 0.0
            if frac_sno_eff[c] > 0.0
                snow_depth[c] = snow_depth[c] + newsnow[c] / (bifall[c] * frac_sno_eff[c])
            else
                snow_depth[c] = 0.0
            end
        else  # h2osno_total == 0
            if newsnow[c] > 0.0
                z_avg = newsnow[c] / bifall[c]
                snow_depth[c] = z_avg / frac_sno_eff[c]
            else
                snow_depth[c] = 0.0
            end
        end
    end

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
function update_snow_depth_and_frac!(
    scf::SnowCoverFractionNiuYang2007,
    # Outputs (modified in place)
    frac_sno::Vector{Float64},
    frac_sno_eff::Vector{Float64},
    snow_depth::Vector{Float64},
    # Inputs
    lun_itype_col::Vector{Int},
    urbpoi::Vector{Bool},
    h2osno_total::Vector{Float64},
    snowmelt::Vector{Float64},
    int_snow::Vector{Float64},
    newsnow::Vector{Float64},
    bifall::Vector{Float64},
    mask::BitVector,
    bounds::UnitRange{Int};
    use_subgrid_fluxes::Bool = false
)
    for c in bounds
        mask[c] || continue

        if h2osno_total[c] == 0.0
            # Reset snow_depth when no snow
            snow_depth[c] = 0.0
        end

        snow_depth[c] = snow_depth[c] + newsnow[c] / bifall[c]

        if snow_depth[c] > 0.0
            frac_sno[c] = tanh(snow_depth[c] / (2.5 * scf.zlnd *
                (min(800.0, (h2osno_total[c] + newsnow[c]) / snow_depth[c]) / 100.0)^1.0))
        else
            frac_sno[c] = 0.0
        end

        # Limit frac_sno when snow water is small but nonzero
        if h2osno_total[c] > 0.0 && h2osno_total[c] < 1.0
            frac_sno[c] = min(frac_sno[c], h2osno_total[c])
        end
    end

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
function add_newsnow_to_intsnow!(
    scf::SnowCoverFractionSwensonLawrence2012,
    int_snow::Vector{Float64},
    newsnow::Vector{Float64},
    h2osno_total::Vector{Float64},
    frac_sno::Vector{Float64},
    mask::BitVector,
    bounds::UnitRange{Int}
)
    for c in bounds
        mask[c] || continue

        if newsnow[c] > 0.0
            # Reset int_snow after accumulation events: make int_snow consistent
            # with new fsno and h2osno_total
            temp_intsnow = (h2osno_total[c] + newsnow[c]) /
                (0.5 * (cos(RPI * (1.0 - max(frac_sno[c], 1.0e-6))^(1.0 / scf.n_melt[c])) + 1.0))
            int_snow[c] = min(1.0e8, temp_intsnow)
        end

        # NOTE(wjs, 2019-07-25): Sean Swenson and Bill Sacks aren't sure whether this
        # extra addition of new_snow is correct: it seems to be double-adding newsnow,
        # but we're not positive that it's wrong. This seems to have been in place ever
        # since the clm45 branch came to the trunk.
        int_snow[c] = int_snow[c] + newsnow[c]
    end

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
function add_newsnow_to_intsnow!(
    scf::SnowCoverFractionNiuYang2007,
    int_snow::Vector{Float64},
    newsnow::Vector{Float64},
    h2osno_total::Vector{Float64},
    frac_sno::Vector{Float64},
    mask::BitVector,
    bounds::UnitRange{Int}
)
    for c in bounds
        mask[c] || continue
        int_snow[c] = int_snow[c] + newsnow[c]
    end

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
