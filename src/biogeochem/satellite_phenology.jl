# ==========================================================================
# Ported from: src/biogeochem/SatellitePhenologyMod.F90
# CLM Satellite Phenology model (SP) ecosystem dynamics (phenology, vegetation)
#
# Allows some subroutines to be used by the CLM Carbon Nitrogen model (CLMCN)
# so that DryDeposition code can get estimates of LAI differences between months.
#
# Public types:
#   SatellitePhenologyData — Interpolation state for monthly vegetation data
#
# Public functions:
#   satellite_phenology_init!    — Allocate and initialize interpolation arrays
#   satellite_phenology_clean!   — Deallocate interpolation arrays
#   satellite_phenology!         — Main SP dynamics: phenology, vegetation
#   interp_monthly_veg!          — Determine if new monthly data is needed
#   read_annual_vegetation!      — Read 12 months of LAI for dry deposition (stub)
#   read_monthly_vegetation!     — Read two consecutive months of veg data (stub)
# ==========================================================================

# ---------------------------------------------------------------------------
# Days per month (non-leap year)
# ---------------------------------------------------------------------------
const NDAYPM = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

# ---------------------------------------------------------------------------
# SatellitePhenologyData — Module-level interpolation state
# Ported from module variables in SatellitePhenologyMod.F90
# ---------------------------------------------------------------------------

"""
    SatellitePhenologyData

Holds the interpolation state for satellite phenology monthly vegetation data.
Contains the saved month index, time weights, and two-month interpolation arrays
for LAI, SAI, and vegetation heights.

Ported from module-level variables in `SatellitePhenologyMod.F90`.
"""
Base.@kwdef mutable struct SatellitePhenologyData{FT<:Real,
                                                  V<:AbstractVector{FT},
                                                  M<:AbstractMatrix{FT}}
    InterpMonths1::Int = -999                                          # saved month index
    timwt::V = zeros(2)                                                # time weights for month 1 and month 2
    mlai2t::M = Matrix{Float64}(undef, 0, 0)            # lai for interpolation (np × 2)
    msai2t::M = Matrix{Float64}(undef, 0, 0)            # sai for interpolation (np × 2)
    mhvt2t::M = Matrix{Float64}(undef, 0, 0)            # top vegetation height for interpolation (np × 2)
    mhvb2t::M = Matrix{Float64}(undef, 0, 0)            # bottom vegetation height for interpolation (np × 2)
end
SatellitePhenologyData{FT}(; kwargs...) where {FT<:Real} =
    SatellitePhenologyData{FT, Vector{FT}, Matrix{FT}}(; kwargs...)
Adapt.@adapt_structure SatellitePhenologyData

# ---------------------------------------------------------------------------
# satellite_phenology_init! — Allocate and initialize
# Ported from: SatellitePhenologyInit in SatellitePhenologyMod.F90
# ---------------------------------------------------------------------------

"""
    satellite_phenology_init!(sp::SatellitePhenologyData, np::Int)

Allocate and initialize all interpolation arrays for `np` patches.
Arrays are initialized to `NaN` (signaling NaN in Fortran).

Ported from `SatellitePhenologyInit` in `SatellitePhenologyMod.F90`.
"""
function satellite_phenology_init!(sp::SatellitePhenologyData{FT}, np::Int) where {FT}
    sp.InterpMonths1 = -999

    sp.mlai2t = fill(FT(NaN), np, 2)
    sp.msai2t = fill(FT(NaN), np, 2)
    sp.mhvt2t = fill(FT(NaN), np, 2)
    sp.mhvb2t = fill(FT(NaN), np, 2)

    sp.timwt = zeros(2)

    return nothing
end

# ---------------------------------------------------------------------------
# satellite_phenology_clean! — Deallocate
# ---------------------------------------------------------------------------

"""
    satellite_phenology_clean!(sp::SatellitePhenologyData)

Deallocate (reset to empty) all interpolation arrays.
"""
function satellite_phenology_clean!(sp::SatellitePhenologyData{FT}) where {FT}
    sp.InterpMonths1 = -999
    sp.timwt = zeros(2)
    sp.mlai2t = Matrix{FT}(undef, 0, 0)
    sp.msai2t = Matrix{FT}(undef, 0, 0)
    sp.mhvt2t = Matrix{FT}(undef, 0, 0)
    sp.mhvb2t = Matrix{FT}(undef, 0, 0)
    return nothing
end

# ---------------------------------------------------------------------------
# satellite_phenology! — Main SP ecosystem dynamics
# Ported from: SatellitePhenology in SatellitePhenologyMod.F90
# ---------------------------------------------------------------------------

"""
    satellite_phenology!(sp, canopystate, waterdiagbulk, patch, mask_patch, bounds_patch;
                          noveg=0, nbrdlf_dcd_brl_shrub=11,
                          use_lai_streams=false, use_fates_sp=false)

Ecosystem dynamics: phenology, vegetation.
Calculates leaf areas (tlai, elai), stem areas (tsai, esai) and height (htop).

Interpolates monthly vegetation data and adjusts LAI/SAI for snow burial.
Snow burial fraction for short vegetation accounts for a 20% bending factor,
as used in Lombardozzi et al. (2018) GRL 45(18), 9889-9897.

Ported from `SatellitePhenology` in `SatellitePhenologyMod.F90`.
"""
# Per-patch satellite-phenology LAI/SAI/height interpolation + snow burial.
# Fully per-patch independent (own-index reads/writes; column gather of
# frac_sno/snow_depth). Float64 literals -> T = eltype(tlai) for Metal.
@kernel function _satphen_kernel!(tlai, tsai, htop, hbot, elai, esai,
        frac_veg_nosno_alb,
        @Const(mask_patch), @Const(column), @Const(itype),
        @Const(frac_sno), @Const(snow_depth),
        @Const(timwt), @Const(mlai2t), @Const(msai2t), @Const(mhvt2t), @Const(mhvb2t),
        noveg::Int, nbrdlf_dcd_brl_shrub::Int, use_lai_streams::Bool, use_fates_sp::Bool)
    T = eltype(tlai)
    p = @index(Global)
    @inbounds if mask_patch[p]
        c = column[p]

        if !use_lai_streams
            tlai[p] = timwt[1] * mlai2t[p, 1] + timwt[2] * mlai2t[p, 2]
        end
        tsai[p] = timwt[1] * msai2t[p, 1] + timwt[2] * msai2t[p, 2]
        htop[p] = timwt[1] * mhvt2t[p, 1] + timwt[2] * mhvt2t[p, 2]
        hbot[p] = timwt[1] * mhvb2t[p, 1] + timwt[2] * mhvb2t[p, 2]

        if itype[p] > noveg && itype[p] <= nbrdlf_dcd_brl_shrub
            ol = min(max(snow_depth[c] - hbot[p], zero(T)), htop[p] - hbot[p])
            fb = one(T) - ol / max(T(1.0e-06), htop[p] - hbot[p])
        else
            fb = one(T) - (max(min(snow_depth[c], max(T(0.05), htop[p] * T(0.8))), zero(T)) /
                         max(T(0.05), htop[p] * T(0.8)))
        end

        if !use_fates_sp
            elai[p] = max(tlai[p] * (one(T) - frac_sno[c]) + tlai[p] * fb * frac_sno[c], zero(T))
            esai[p] = max(tsai[p] * (one(T) - frac_sno[c]) + tsai[p] * fb * frac_sno[c], zero(T))
            if elai[p] < T(0.05)
                elai[p] = zero(T)
            end
            if esai[p] < T(0.05)
                esai[p] = zero(T)
            end
            if (elai[p] + esai[p]) >= T(0.05)
                frac_veg_nosno_alb[p] = 1
            else
                frac_veg_nosno_alb[p] = 0
            end
        end
    end
end

function satellite_phenology!(sp::SatellitePhenologyData,
                               canopystate::CanopyStateData,
                               waterdiagbulk::WaterDiagnosticBulkData,
                               patch::PatchData,
                               mask_patch::AbstractVector{Bool},
                               bounds_patch::UnitRange{Int};
                               noveg::Int = 0,
                               nbrdlf_dcd_brl_shrub::Int = 11,
                               use_lai_streams::Bool = false,
                               use_fates_sp::Bool = false)
    # Note: lai_interp (use_lai_streams) is not ported; skip that call here.
    _launch!(_satphen_kernel!,
        canopystate.tlai_patch, canopystate.tsai_patch, canopystate.htop_patch,
        canopystate.hbot_patch, canopystate.elai_patch, canopystate.esai_patch,
        canopystate.frac_veg_nosno_alb_patch,
        mask_patch, patch.column, patch.itype,
        waterdiagbulk.frac_sno_col, waterdiagbulk.snow_depth_col,
        sp.timwt, sp.mlai2t, sp.msai2t, sp.mhvt2t, sp.mhvb2t,
        noveg, nbrdlf_dcd_brl_shrub, use_lai_streams, use_fates_sp)
    return nothing
end

# ---------------------------------------------------------------------------
# interp_monthly_veg! — Determine if new months need reading
# Ported from: interpMonthlyVeg in SatellitePhenologyMod.F90
# ---------------------------------------------------------------------------

"""
    interp_monthly_veg!(sp, canopystate, bounds_patch;
                         kmo, kda, do_read_monthly=false)

Determine time weights for monthly vegetation data interpolation.
Updates `sp.timwt` based on the current month and day. If the required
months change, signals that new data should be read.

In the Fortran source, this calls `readMonthlyVegetation` when the
interpolation months change. Here, the I/O is separated — callers
should check the return value or `sp.InterpMonths1` for changes.

Ported from `interpMonthlyVeg` in `SatellitePhenologyMod.F90`.

Returns `(months, needs_read)` where `months` is a 2-element tuple of the
months to interpolate, and `needs_read` is true if new data should be read.
"""
function interp_monthly_veg!(sp::SatellitePhenologyData;
                              kmo::Int,
                              kda::Int)
    t = (kda - 0.5) / NDAYPM[kmo]
    it1 = floor(Int, t + 0.5)
    it2 = it1 + 1
    months1 = kmo + it1 - 1
    months2 = kmo + it2 - 1
    if months1 < 1
        months1 = 12
    end
    if months2 > 12
        months2 = 1
    end
    # Bulk-write the two interpolation weights so this works whether sp.timwt is
    # a host Array or a device array (a scalar `sp.timwt[i] = ` would scalar-index
    # a GPU array). Byte-identical to the per-element writes on the CPU.
    Twt = eltype(sp.timwt)
    wt1 = Twt((it1 + 0.5) - t)
    copyto!(sp.timwt, Twt[wt1, one(Twt) - wt1])

    needs_read = false
    if sp.InterpMonths1 != months1
        sp.InterpMonths1 = months1
        needs_read = true
    end

    return (months1, months2), needs_read
end

# ---------------------------------------------------------------------------
# read_annual_vegetation! — Read 12 months of LAI for dry deposition
# Ported from: readAnnualVegetation in SatellitePhenologyMod.F90
# ---------------------------------------------------------------------------

"""
    read_annual_vegetation!(canopystate, patch, bounds_patch;
                             monthly_lai, noveg=0, maxveg=78)

Read 12 months of LAI data for dry deposition.

In the Fortran source, this reads from a NetCDF surface data file. Here the
monthly LAI data is provided directly as a 3D array `monthly_lai[g, l, k]`
where `g` is gridcell, `l` is PFT index (0:maxveg), and `k` is month (1:12).

Only vegetated patches get nonzero values; non-vegetated patches get 0.

Ported from `readAnnualVegetation` in `SatellitePhenologyMod.F90`.
"""
# Per-patch annual-LAI read (internal k=1:12 loop in-thread). No mask: all
# patches in [pmin,pmax]. monthly_lai is [g, l+1, k] (gather).
@kernel function _satphen_annlai_kernel!(annlai, @Const(gridcell), @Const(itype),
        @Const(monthly_lai), noveg::Int, pmin::Int, pmax::Int)
    T = eltype(annlai)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax
        g = gridcell[p]
        if itype[p] != noveg
            l = itype[p]
            for k in 1:12
                annlai[p, k] = monthly_lai[g, l + 1, k]
            end
        else
            for k in 1:12
                annlai[p, k] = zero(T)
            end
        end
    end
end

function read_annual_vegetation!(canopystate::CanopyStateData,
                                  patch::PatchData,
                                  bounds_patch::UnitRange{Int};
                                  monthly_lai::AbstractArray{<:Real,3},
                                  noveg::Int = 0,
                                  maxveg::Int = 78)
    _launch!(_satphen_annlai_kernel!, canopystate.annlai_patch,
        patch.gridcell, patch.itype, monthly_lai, noveg,
        first(bounds_patch), last(bounds_patch);
        ndrange = length(patch.itype))
    return nothing
end

# ---------------------------------------------------------------------------
# read_monthly_vegetation! — Read two consecutive months of veg data
# Ported from: readMonthlyVegetation in SatellitePhenologyMod.F90
# ---------------------------------------------------------------------------

"""
    read_monthly_vegetation!(sp, canopystate, patch, bounds_patch;
                              monthly_lai, monthly_sai, monthly_height_top,
                              monthly_height_bot, months,
                              noveg=0, maxveg=78)

Read monthly vegetation data for two consecutive months.

In the Fortran source, this reads from a NetCDF surface data file. Here the
monthly data is provided directly as 3D arrays `[g, l+1, month]` where `g`
is gridcell, `l` is PFT index (0:maxveg), and month is 1:12.

After reading, also computes `mlaidiff_patch` = mlai2t(:,1) - mlai2t(:,2).

Ported from `readMonthlyVegetation` in `SatellitePhenologyMod.F90`.
"""
# Per-patch two-month veg read + mlaidiff (internal k=1:2 in-thread). months
# passed as 2 Int scalars (m1,m2). monthly_* are [g, l+1, month] (gather).
@kernel function _satphen_monthly_kernel!(mlai2t, msai2t, mhvt2t, mhvb2t, mlaidiff,
        @Const(gridcell), @Const(itype),
        @Const(monthly_lai), @Const(monthly_sai), @Const(monthly_height_top),
        @Const(monthly_height_bot), m1::Int, m2::Int, noveg::Int, pmin::Int, pmax::Int)
    T = eltype(mlai2t)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax
        g = gridcell[p]
        for k in 1:2
            mk = k == 1 ? m1 : m2
            if itype[p] != noveg
                l = itype[p]
                mlai2t[p, k] = monthly_lai[g, l + 1, mk]
                msai2t[p, k] = monthly_sai[g, l + 1, mk]
                mhvt2t[p, k] = monthly_height_top[g, l + 1, mk]
                mhvb2t[p, k] = monthly_height_bot[g, l + 1, mk]
            else
                mlai2t[p, k] = zero(T)
                msai2t[p, k] = zero(T)
                mhvt2t[p, k] = zero(T)
                mhvb2t[p, k] = zero(T)
            end
        end
        mlaidiff[p] = mlai2t[p, 1] - mlai2t[p, 2]
    end
end

function read_monthly_vegetation!(sp::SatellitePhenologyData,
                                   canopystate::CanopyStateData,
                                   patch::PatchData,
                                   bounds_patch::UnitRange{Int};
                                   monthly_lai::AbstractArray{<:Real,3},
                                   monthly_sai::AbstractArray{<:Real,3},
                                   monthly_height_top::AbstractArray{<:Real,3},
                                   monthly_height_bot::AbstractArray{<:Real,3},
                                   months::Tuple{Int,Int},
                                   noveg::Int = 0,
                                   maxveg::Int = 78)
    _launch!(_satphen_monthly_kernel!, sp.mlai2t, sp.msai2t, sp.mhvt2t, sp.mhvb2t,
        canopystate.mlaidiff_patch, patch.gridcell, patch.itype,
        monthly_lai, monthly_sai, monthly_height_top, monthly_height_bot,
        months[1], months[2], noveg, first(bounds_patch), last(bounds_patch);
        ndrange = length(patch.itype))
    return nothing
end
