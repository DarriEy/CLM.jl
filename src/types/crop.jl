# ==========================================================================
# Ported from: src/biogeochem/CropType.F90
# Crop type data structure, allocation, initialization, and accumulation
# ==========================================================================

# --- Crop phenology phase constants ---
const cphase_not_planted = 0.0
const cphase_planted     = 1.0
const cphase_leafemerge  = 2.0
const cphase_grainfill   = 3.0
const cphase_harvest     = 4.0

# --- Base temperature mapping modes ---
const BASET_MAP_CONSTANT = "constant"
const BASET_MAP_LATVARY  = "varytropicsbylat"

"""
    CropData

Crop model state variables. Holds phenology, sowing/harvest tracking,
growing degree-day accumulators, and related per-patch fields.

Ported from `crop_type` in `CropType.F90`.
"""
Base.@kwdef mutable struct CropData{FT<:Real,
                        V<:AbstractVector{FT},
                        M<:AbstractMatrix{FT},
                        VI<:AbstractVector{<:Integer},
                        MI<:AbstractMatrix{<:Integer},
                        VB<:AbstractVector{Bool}}
    # --- Scalar / config fields ---
    baset_mapping              ::String  = BASET_MAP_CONSTANT
    baset_latvary_intercept    ::FT = 12.0
    baset_latvary_slope        ::FT = 0.4

    # --- Patch-level 1D integer fields ---
    nyrs_crop_active_patch     ::VI     = Int[]           # years this crop patch has been active
    harvdate_patch             ::VI     = Int[]           # most recent harvest date (999 if currently/never planted)
    sowing_reason_patch        ::VI     = Int[]           # reason for most recent sowing
    sowing_count               ::VI     = Int[]           # number of sowing events this year
    harvest_count              ::VI     = Int[]           # number of harvest events this year

    # --- Patch-level 1D boolean fields ---
    croplive_patch             ::VB    = Bool[]          # true if planted, not harvested
    sown_in_this_window        ::VB    = Bool[]          # true if sown during current sowing window

    # --- Patch-level 1D real fields ---
    fertnitro_patch            ::V = Float64[]       # fertilizer nitrogen (gN/m2/yr)
    gddtsoi_patch              ::V = Float64[]       # GDD from planting (top two soil layers) (ddays)
    vf_patch                   ::V = Float64[]       # vernalization factor for cereal
    cphase_patch               ::V = Float64[]       # phenology phase (see cphase_* constants)
    latbaset_patch             ::V = Float64[]       # latitude-varying baset for hui (degree C)
    hui_patch                  ::V = Float64[]       # heat unit index (ddays)
    gddaccum_patch             ::V = Float64[]       # GDD from planting (air) (ddays)
    gdd20_baseline_patch       ::V = Float64[]       # GDD20 baseline (ddays)
    gdd20_season_start_patch   ::V = Float64[]       # GDD20 season start (day of year)
    gdd20_season_end_patch     ::V = Float64[]       # GDD20 season end (day of year)

    # --- Patch-level 2D integer fields (patch × mxsowings) ---
    rx_swindow_starts_thisyr_patch ::MI     = Matrix{Int}(undef, 0, 0)     # prescribed sowing window starts
    rx_swindow_ends_thisyr_patch   ::MI     = Matrix{Int}(undef, 0, 0)     # prescribed sowing window ends

    # --- Patch-level 2D real fields (patch × mxsowings) ---
    rx_cultivar_gdds_thisyr_patch  ::M = Matrix{Float64}(undef, 0, 0) # cultivar GDD targets
    sdates_thisyr_patch            ::M = Matrix{Float64}(undef, 0, 0) # actual sowing dates this year
    swindow_starts_thisyr_patch    ::M = Matrix{Float64}(undef, 0, 0) # sowing window starts this year
    swindow_ends_thisyr_patch      ::M = Matrix{Float64}(undef, 0, 0) # sowing window ends this year
    sowing_reason_thisyr_patch     ::M = Matrix{Float64}(undef, 0, 0) # reason for each sowing this year

    # --- Patch-level 2D real fields (patch × mxharvests) ---
    sdates_perharv_patch           ::M = Matrix{Float64}(undef, 0, 0) # sowing dates for harvested crops
    syears_perharv_patch           ::M = Matrix{Float64}(undef, 0, 0) # sowing years for harvested crops
    hdates_thisyr_patch            ::M = Matrix{Float64}(undef, 0, 0) # harvest dates this year
    gddaccum_thisyr_patch          ::M = Matrix{Float64}(undef, 0, 0) # accumulated GDD at harvest this year
    hui_thisyr_patch               ::M = Matrix{Float64}(undef, 0, 0) # accumulated HUI at harvest this year
    sowing_reason_perharv_patch    ::M = Matrix{Float64}(undef, 0, 0) # sowing reason for harvested crops
    harvest_reason_thisyr_patch    ::M = Matrix{Float64}(undef, 0, 0) # reason for each harvest this year
end

CropData{FT}(; kwargs...) where {FT<:Real} =
    CropData{FT, Vector{FT}, Matrix{FT}, Vector{Int}, Matrix{Int}, Vector{Bool}}(; kwargs...)
Adapt.@adapt_structure CropData


# ---------------------------------------------------------------------------
# Helper constructors (reuse if already defined)
# ---------------------------------------------------------------------------
if !@isdefined(nanvec)
    nanvec(n) = fill(NaN, n)
end
if !@isdefined(nanmat)
    nanmat(r, c) = fill(NaN, r, c)
end

"""
    crop_init!(cr, np; mxsowings=MXSOWINGS, mxharvests=MXHARVESTS)

Allocate all fields of `CropData` for `np` patches.
Corresponds to `InitAllocate` in the Fortran source.
"""
function crop_init!(cr::CropData, np::Int;
                    mxsowings::Int=MXSOWINGS,
                    mxharvests::Int=MXHARVESTS)

    # --- Patch-level 1D integer ---
    cr.nyrs_crop_active_patch  = zeros(Int, np)
    cr.harvdate_patch          = fill(typemax(Int), np)  # huge(1) in Fortran
    cr.sowing_reason_patch     = fill(-1, np)
    cr.sowing_count            = zeros(Int, np)
    cr.harvest_count           = zeros(Int, np)

    # --- Patch-level 1D boolean ---
    cr.croplive_patch          = falses(np)
    cr.sown_in_this_window     = falses(np)

    # --- Patch-level 1D real ---
    cr.fertnitro_patch         = fill(SPVAL, np)
    cr.hui_patch               = fill(SPVAL, np)
    cr.gddaccum_patch          = fill(SPVAL, np)
    cr.gddtsoi_patch           = fill(SPVAL, np)
    cr.vf_patch                = zeros(np)
    cr.cphase_patch            = fill(cphase_not_planted, np)
    cr.latbaset_patch          = fill(SPVAL, np)
    cr.gdd20_baseline_patch    = fill(SPVAL, np)
    cr.gdd20_season_start_patch = fill(SPVAL, np)
    cr.gdd20_season_end_patch  = fill(SPVAL, np)

    # --- Patch-level 2D integer (patch × mxsowings) ---
    cr.rx_swindow_starts_thisyr_patch = fill(-1, np, mxsowings)
    cr.rx_swindow_ends_thisyr_patch   = fill(-1, np, mxsowings)

    # --- Patch-level 2D real (patch × mxsowings) ---
    cr.rx_cultivar_gdds_thisyr_patch  = fill(SPVAL, np, mxsowings)
    cr.sdates_thisyr_patch            = fill(SPVAL, np, mxsowings)
    cr.swindow_starts_thisyr_patch    = fill(SPVAL, np, mxsowings)
    cr.swindow_ends_thisyr_patch      = fill(SPVAL, np, mxsowings)
    cr.sowing_reason_thisyr_patch     = fill(SPVAL, np, mxsowings)

    # --- Patch-level 2D real (patch × mxharvests) ---
    cr.sdates_perharv_patch           = fill(SPVAL, np, mxharvests)
    cr.syears_perharv_patch           = fill(SPVAL, np, mxharvests)
    cr.hdates_thisyr_patch            = fill(SPVAL, np, mxharvests)
    cr.gddaccum_thisyr_patch          = fill(SPVAL, np, mxharvests)
    cr.hui_thisyr_patch               = fill(SPVAL, np, mxharvests)
    cr.sowing_reason_perharv_patch    = fill(SPVAL, np, mxharvests)
    cr.harvest_reason_thisyr_patch    = fill(SPVAL, np, mxharvests)

    return nothing
end

"""
    crop_read_nml!(cr; baset_mapping, baset_latvary_intercept, baset_latvary_slope)

Set namelist parameters for the crop type.
Corresponds to `ReadNML` in the Fortran source.
"""
function crop_read_nml!(cr::CropData;
                        baset_mapping::String=BASET_MAP_CONSTANT,
                        baset_latvary_intercept::Real=12.0,
                        baset_latvary_slope::Real=0.4)

    if baset_mapping != BASET_MAP_CONSTANT && baset_mapping != BASET_MAP_LATVARY
        error("Bad value for baset_mapping: $baset_mapping")
    end

    cr.baset_mapping           = baset_mapping
    cr.baset_latvary_intercept = baset_latvary_intercept
    cr.baset_latvary_slope     = baset_latvary_slope

    return nothing
end

"""
    crop_init_history!(cr, bounds_patch)

Stub for history field registration (no-op in Julia port).
Corresponds to `InitHistory` in the Fortran source.
"""
function crop_init_history!(cr::CropData{FT}, bounds_patch::UnitRange{Int}) where {FT}
    return nothing
end

"""
    crop_init_cold!(cr, bounds_patch;
                    patch_landunit=nothing, lun_itype=nothing, istcrop=nothing,
                    patch_gridcell=nothing, patch_itype=nothing,
                    fert_cft=nothing, pftcon_baset=nothing, grc_latdeg=nothing)

Cold-start initialization of crop fields.
Corresponds to `InitCold` in the Fortran source.

For each crop patch (where `lun_itype[patch_landunit[p]] == istcrop`),
sets fertilizer nitrogen and optionally latitude-varying base temperature.
"""
function crop_init_cold!(cr::CropData, bounds_patch::UnitRange{Int};
                         patch_landunit::Union{AbstractVector{<:Integer},Nothing}=nothing,
                         lun_itype::Union{AbstractVector{<:Integer},Nothing}=nothing,
                         istcrop::Union{Int,Nothing}=nothing,
                         patch_gridcell::Union{AbstractVector{<:Integer},Nothing}=nothing,
                         patch_itype::Union{AbstractVector{<:Integer},Nothing}=nothing,
                         fert_cft::Union{Matrix{Float64},Nothing}=nothing,
                         pftcon_baset::Union{Vector{Float64},Nothing}=nothing,
                         grc_latdeg::Union{Vector{Float64},Nothing}=nothing)

    latvary_baset = cr.baset_mapping == BASET_MAP_LATVARY

    if !latvary_baset
        for p in bounds_patch
            cr.latbaset_patch[p] = NaN
        end
    end

    for p in bounds_patch
        cr.nyrs_crop_active_patch[p] = 0

        # If linkage data is provided, set crop-specific fields
        if patch_landunit !== nothing && lun_itype !== nothing && istcrop !== nothing
            l = patch_landunit[p]
            if lun_itype[l] == istcrop
                if fert_cft !== nothing && patch_gridcell !== nothing && patch_itype !== nothing
                    g   = patch_gridcell[p]
                    ivt = patch_itype[p] + 1  # 0-based Fortran → 1-based Julia
                    cr.fertnitro_patch[p] = fert_cft[g, ivt]
                end

                if latvary_baset && pftcon_baset !== nothing && grc_latdeg !== nothing &&
                   patch_gridcell !== nothing && patch_itype !== nothing
                    g   = patch_gridcell[p]
                    ivt = patch_itype[p] + 1  # 0-based Fortran → 1-based Julia
                    cr.latbaset_patch[p] = latbaset(pftcon_baset[ivt], grc_latdeg[g],
                                                    cr.baset_latvary_intercept,
                                                    cr.baset_latvary_slope)
                end
            end
        end
    end

    return nothing
end

"""
    crop_init_acc_buffer!(cr)

Stub for accumulation buffer initialization (no-op in Julia port).
Corresponds to `InitAccBuffer` in the Fortran source.
"""
function crop_init_acc_buffer!(cr::CropData{FT}) where {FT}
    return nothing
end

"""
    crop_init_acc_vars!(cr, bounds_patch)

Stub for accumulation variable extraction (no-op in Julia port).
Corresponds to `InitAccVars` in the Fortran source.
"""
function crop_init_acc_vars!(cr::CropData{FT}, bounds_patch::UnitRange{Int}) where {FT}
    return nothing
end

"""
    crop_restart!(cr, bounds_patch)

Stub for restart read/write (no-op in Julia port).
Corresponds to `Restart` in the Fortran source.
"""
function crop_restart!(cr::CropData{FT}, bounds_patch::UnitRange{Int}) where {FT}
    return nothing
end

"""
    crop_update_acc_vars!(cr, bounds_patch, t_ref2m_patch, t_soisno_col;
                          patch, col, pftcon, dtime, nlevsno=NLEVSNO)

Accumulate the three crop `runaccum` fields, every timestep:

- `hui_patch`     — heat unit index since planting (HUI)
- `gddaccum_patch`— growing degree-days since planting, from air temperature
- `gddtsoi_patch` — growing degree-days since planting, from the top two soil layers

For a live crop (`croplive_patch`) each accumulates a daily increment

    max(0, min(mxtmp[ivt], T - (TFRZ + baset))) * dtime/SECSPDAY

(with `T` = `t_ref2m` for HUI/GDDACCUM and the dz-weighted mean of the top two
soil layers for GDDTSOI), scaled by the vernalization factor `vf` for winter
wheat. When the crop is not live the accumulator is RESET to zero — Fortran's
per-patch `markreset_accum_field`, which (per `update_accum_field_runaccum`)
resets WITHOUT also accumulating on that step.

`baset` is `pftcon.baset[ivt]`, except under `baset_mapping ==
"varytropicsbylat"` for (spring wheat, sugarcane), where the latitude-varying
`latbaset_patch` is used instead.

BUG HISTORY: this was a no-op stub AND had no call site. `hui_patch`,
`gddaccum_patch` and `gddtsoi_patch` are allocated to SPVAL (1e36) and were never
written by anything, yet they ARE read — by `cn_crop_phenology!` (`hui >=
huigrain` gates grain fill), by crop allocation, and by `n_dynamics`' soybean
fixation. A crop run therefore saw HUI pinned at SPVAL, i.e. "infinitely mature"
from the first step.

Fortran reconciles the accumulator buffer against external writes to `hui_patch`
by CropPhenology (extract, subtract, re-update). Here `hui_patch` IS the
accumulator, so that reconciliation is a no-op by construction and is omitted.

Ported from `CropUpdateAccVars` in `CropType.F90`.
"""
function crop_update_acc_vars!(cr::CropData, bounds_patch::UnitRange{Int},
                                t_ref2m_patch::AbstractVector{<:Real},
                                t_soisno_col::AbstractMatrix{<:Real};
                                patch::PatchData,
                                col::ColumnData,
                                pftcon,
                                dtime::Real,
                                nlevsno::Int = varpar.nlevsno)
    isempty(bounds_patch) && return nothing
    length(cr.hui_patch) >= last(bounds_patch) || return nothing

    fracday = dtime / SECSPDAY
    latvary = (cr.baset_mapping == BASET_MAP_LATVARY)

    @inbounds for p in bounds_patch
        patch.active[p] || continue
        ivt = patch.itype[p]
        (ivt >= 1 && ivt <= length(pftcon.baset)) || continue

        if !cr.croplive_patch[p]
            # markreset_accum_field(..., p): reset, and do NOT accumulate this step.
            cr.hui_patch[p]      = 0.0
            cr.gddaccum_patch[p] = 0.0
            cr.gddtsoi_patch[p]  = 0.0
            continue
        end

        # Base temperature: latitude-varying for spring wheat / sugarcane under
        # the varytropicsbylat mapping, else the per-PFT constant.
        baset = if latvary && (ivt == nswheat || ivt == nirrig_swheat ||
                               ivt == nsugarcane || ivt == nirrig_sugarcane)
            cr.latbaset_patch[p]
        else
            pftcon.baset[ivt]
        end
        isfinite(baset) || continue
        mxtmp = pftcon.mxtmp[ivt]
        winter_wheat = (ivt == nwwheat || ivt == nirrig_wwheat)

        # --- HUI / GDDACCUM: air-temperature degree-days ---
        t2m = t_ref2m_patch[p]
        if isfinite(t2m)
            incr = max(0.0, min(mxtmp, t2m - (TFRZ + baset))) * fracday
            if winter_wheat
                incr *= cr.vf_patch[p]
            end
            cur_h = cr.hui_patch[p]
            cur_g = cr.gddaccum_patch[p]
            # runaccum clamps to [0, 99999]; a non-finite/SPVAL start re-seeds at 0.
            cr.hui_patch[p]      = min(max((isfinite(cur_h) && cur_h < SPVAL ? cur_h : 0.0) + incr, 0.0), 99999.0)
            cr.gddaccum_patch[p] = min(max((isfinite(cur_g) && cur_g < SPVAL ? cur_g : 0.0) + incr, 0.0), 99999.0)
        end

        # --- GDDTSOI: degree-days from the dz-weighted top TWO soil layers ---
        # (agroibis computes this to 0.05 m). t_soisno is snow-padded here.
        c = patch.column[p]
        dz1 = col.dz[c, nlevsno + 1]
        dz2 = col.dz[c, nlevsno + 2]
        ts1 = t_soisno_col[c, nlevsno + 1]
        ts2 = t_soisno_col[c, nlevsno + 2]
        if isfinite(ts1) && isfinite(ts2) && (dz1 + dz2) > 0
            tsoi = (ts1 * dz1 + ts2 * dz2) / (dz1 + dz2)
            incr = max(0.0, min(mxtmp, tsoi - (TFRZ + baset))) * fracday
            # UPSTREAM BUG, REPRODUCED DELIBERATELY: Fortran guards this branch with
            # `if (ivt == nwwheat .or. ivt == nwwheat)` — the second disjunct is a
            # typo for `nirrig_wwheat`. So GDDTSOI applies the vernalization factor
            # to RAINFED winter wheat only, while HUI/GDDACCUM above apply it to both
            # rainfed and irrigated. Matching CTSM bit-for-bit means keeping this.
            if ivt == nwwheat
                incr *= cr.vf_patch[p]
            end
            cur_t = cr.gddtsoi_patch[p]
            cr.gddtsoi_patch[p] = min(max((isfinite(cur_t) && cur_t < SPVAL ? cur_t : 0.0) + incr, 0.0), 99999.0)
        end
    end

    return nothing
end

"""
    crop_increment_year!(cr, mask_pcrop, bounds_patch;
                         kmo, kda, mcsec, is_first_step)

Increment the crop year counter for active crop patches at year boundaries.
Should be called every time step.

Corresponds to `CropIncrementYear` in the Fortran source.
"""
function crop_increment_year!(cr::CropData, mask_pcrop::BitVector,
                              bounds_patch::UnitRange{Int};
                              kmo::Int, kda::Int, mcsec::Int,
                              is_first_step::Bool)

    if kmo == 1 && kda == 1 && mcsec == 0 && !is_first_step
        for p in bounds_patch
            mask_pcrop[p] || continue
            cr.nyrs_crop_active_patch[p] += 1
        end
    end

    return nothing
end

"""
    latbaset(baset, latdeg, baset_latvary_intercept, baset_latvary_slope)

Compute latitude-varying base temperature for heat unit index calculation.

Ported from `latbaset` function in `CropType.F90`.
"""
function latbaset(baset::Real, latdeg::Real,
                  baset_latvary_intercept::Real,
                  baset_latvary_slope::Real)
    return baset + baset_latvary_intercept - min(baset_latvary_intercept,
                                                  baset_latvary_slope * abs(latdeg))
end

"""
    crop_check_dates!()

Stub for restart date checking (no-op in Julia port).
Corresponds to `checkDates` in the Fortran source.
"""
function crop_check_dates!()
    return nothing
end
