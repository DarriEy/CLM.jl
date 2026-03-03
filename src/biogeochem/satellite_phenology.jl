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
Base.@kwdef mutable struct SatellitePhenologyData
    InterpMonths1::Int = -999                                          # saved month index
    timwt::Vector{Float64} = zeros(2)                                  # time weights for month 1 and month 2
    mlai2t::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)            # lai for interpolation (np × 2)
    msai2t::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)            # sai for interpolation (np × 2)
    mhvt2t::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)            # top vegetation height for interpolation (np × 2)
    mhvb2t::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)            # bottom vegetation height for interpolation (np × 2)
end

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
function satellite_phenology_init!(sp::SatellitePhenologyData, np::Int)
    sp.InterpMonths1 = -999

    sp.mlai2t = fill(NaN, np, 2)
    sp.msai2t = fill(NaN, np, 2)
    sp.mhvt2t = fill(NaN, np, 2)
    sp.mhvb2t = fill(NaN, np, 2)

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
function satellite_phenology_clean!(sp::SatellitePhenologyData)
    sp.InterpMonths1 = -999
    sp.timwt = zeros(2)
    sp.mlai2t = Matrix{Float64}(undef, 0, 0)
    sp.msai2t = Matrix{Float64}(undef, 0, 0)
    sp.mhvt2t = Matrix{Float64}(undef, 0, 0)
    sp.mhvb2t = Matrix{Float64}(undef, 0, 0)
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
function satellite_phenology!(sp::SatellitePhenologyData,
                               canopystate::CanopyStateData,
                               waterdiagbulk::WaterDiagnosticBulkData,
                               patch::PatchData,
                               mask_patch::BitVector,
                               bounds_patch::UnitRange{Int};
                               noveg::Int = 0,
                               nbrdlf_dcd_brl_shrub::Int = 11,
                               use_lai_streams::Bool = false,
                               use_fates_sp::Bool = false)

    # Aliases matching Fortran associate block
    frac_sno           = waterdiagbulk.frac_sno_col
    snow_depth         = waterdiagbulk.snow_depth_col
    tlai               = canopystate.tlai_patch
    tsai               = canopystate.tsai_patch
    elai               = canopystate.elai_patch
    esai               = canopystate.esai_patch
    htop               = canopystate.htop_patch
    hbot               = canopystate.hbot_patch
    frac_veg_nosno_alb = canopystate.frac_veg_nosno_alb_patch

    timwt  = sp.timwt
    mlai2t = sp.mlai2t
    msai2t = sp.msai2t
    mhvt2t = sp.mhvt2t
    mhvb2t = sp.mhvb2t

    # Note: lai_interp (use_lai_streams) is not ported; skip that call here.
    # When lai_streams infrastructure is ported, add:
    #   if use_lai_streams
    #       lai_interp!(bounds_patch, canopystate)
    #   end

    for p in bounds_patch
        mask_patch[p] || continue
        c = patch.column[p]

        # Interpolate leaf area index, stem area index, and vegetation heights
        # between two monthly values
        if !use_lai_streams
            tlai[p] = timwt[1] * mlai2t[p, 1] + timwt[2] * mlai2t[p, 2]
        end

        tsai[p] = timwt[1] * msai2t[p, 1] + timwt[2] * msai2t[p, 2]
        htop[p] = timwt[1] * mhvt2t[p, 1] + timwt[2] * mhvt2t[p, 2]
        hbot[p] = timwt[1] * mhvb2t[p, 1] + timwt[2] * mhvb2t[p, 2]

        # Adjust lai and sai for burying by snow. If exposed lai and sai
        # are less than 0.05, set equal to zero to prevent numerical
        # problems associated with very small lai and sai.

        # Snow burial fraction for short vegetation (e.g. grasses, crops)
        # changes with vegetation height. Accounts for a 20% bending factor.
        # NOTE: This snow burial code is duplicated in CNVegStructUpdateMod.

        if patch.itype[p] > noveg && patch.itype[p] <= nbrdlf_dcd_brl_shrub
            ol = min(max(snow_depth[c] - hbot[p], 0.0), htop[p] - hbot[p])
            fb = 1.0 - ol / max(1.0e-06, htop[p] - hbot[p])
        else
            fb = 1.0 - (max(min(snow_depth[c], max(0.05, htop[p] * 0.8)), 0.0) /
                         max(0.05, htop[p] * 0.8))
        end

        # Area weight by snow covered fraction
        if !use_fates_sp
            # Do not set these in FATES_SP mode as they turn on the 'vegsol' filter
            # and also are duplicated by the FATES variables
            elai[p] = max(tlai[p] * (1.0 - frac_sno[c]) + tlai[p] * fb * frac_sno[c], 0.0)
            esai[p] = max(tsai[p] * (1.0 - frac_sno[c]) + tsai[p] * fb * frac_sno[c], 0.0)
            if elai[p] < 0.05
                elai[p] = 0.0
            end
            if esai[p] < 0.05
                esai[p] = 0.0
            end

            # Fraction of vegetation free of snow
            if (elai[p] + esai[p]) >= 0.05
                frac_veg_nosno_alb[p] = 1
            else
                frac_veg_nosno_alb[p] = 0
            end
        end
    end

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
    sp.timwt[1] = (it1 + 0.5) - t
    sp.timwt[2] = 1.0 - sp.timwt[1]

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
function read_annual_vegetation!(canopystate::CanopyStateData,
                                  patch::PatchData,
                                  bounds_patch::UnitRange{Int};
                                  monthly_lai::Array{Float64,3},
                                  noveg::Int = 0,
                                  maxveg::Int = 78)
    annlai = canopystate.annlai_patch

    for k in 1:12
        for p in bounds_patch
            g = patch.gridcell[p]
            if patch.itype[p] != noveg
                # monthly_lai is indexed [g, l+1, k] since Julia is 1-based
                # PFT l in Fortran (0:maxveg) maps to index l+1 in Julia
                l = patch.itype[p]
                annlai[p, k] = monthly_lai[g, l + 1, k]
            else
                annlai[p, k] = 0.0
            end
        end
    end

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
function read_monthly_vegetation!(sp::SatellitePhenologyData,
                                   canopystate::CanopyStateData,
                                   patch::PatchData,
                                   bounds_patch::UnitRange{Int};
                                   monthly_lai::Array{Float64,3},
                                   monthly_sai::Array{Float64,3},
                                   monthly_height_top::Array{Float64,3},
                                   monthly_height_bot::Array{Float64,3},
                                   months::Tuple{Int,Int},
                                   noveg::Int = 0,
                                   maxveg::Int = 78)

    mlai2t = sp.mlai2t
    msai2t = sp.msai2t
    mhvt2t = sp.mhvt2t
    mhvb2t = sp.mhvb2t

    for k in 1:2
        mk = k == 1 ? months[1] : months[2]

        for p in bounds_patch
            g = patch.gridcell[p]
            if patch.itype[p] != noveg
                l = patch.itype[p]
                mlai2t[p, k] = monthly_lai[g, l + 1, mk]
                msai2t[p, k] = monthly_sai[g, l + 1, mk]
                mhvt2t[p, k] = monthly_height_top[g, l + 1, mk]
                mhvb2t[p, k] = monthly_height_bot[g, l + 1, mk]
            else
                mlai2t[p, k] = 0.0
                msai2t[p, k] = 0.0
                mhvt2t[p, k] = 0.0
                mhvb2t[p, k] = 0.0
            end
        end
    end

    # Compute mlaidiff
    for p in bounds_patch
        canopystate.mlaidiff_patch[p] = mlai2t[p, 1] - mlai2t[p, 2]
    end

    return nothing
end
