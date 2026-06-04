# ==========================================================================
# Ported from: src/biogeochem/CNNDynamicsMod.F90
# Mineral nitrogen dynamics: deposition, fixation, fertilization, soy fixation
# for coupled carbon-nitrogen code.
#
# Public functions:
#   n_dynamics_read_nml!     — Read/set N dynamics parameters
#   n_deposition!            — Update N deposition from atmospheric forcing
#   n_fixation!              — Update N fixation as f(NPP)
#   n_fert!                  — Update N fertilizer for crops (p2c)
#   n_soyfix!                — N fixation for soybeans
#   n_free_living_fixation!  — Free-living N fixation as f(ET)
# ==========================================================================

# ---------------------------------------------------------------------------
# NDynamicsParams — module-level parameters
# Ported from params_type in CNNDynamicsMod.F90
# ---------------------------------------------------------------------------

"""
    NDynamicsParams

Parameters for mineral nitrogen dynamics.
Holds intercept and slope for free-living N fixation with annual ET.

Ported from `params_type` in `CNNDynamicsMod.F90`.
"""
Base.@kwdef mutable struct NDynamicsParams
    freelivfix_intercept::Float64 = 0.0117   # intercept of line of free living fixation with annual ET
    freelivfix_slope_wET::Float64 = 0.0006   # slope of line of free living fixation with annual ET
end

# ---------------------------------------------------------------------------
# n_dynamics_read_nml! — read/set namelist parameters
# Ported from CNNDynamicsReadNML in CNNDynamicsMod.F90
# ---------------------------------------------------------------------------

"""
    n_dynamics_read_nml!(params; freelivfix_intercept, freelivfix_slope_wET)

Set N dynamics parameters from keyword arguments.
Corresponds to `CNNDynamicsReadNML` in the Fortran source.
"""
function n_dynamics_read_nml!(params::NDynamicsParams;
                               freelivfix_intercept::Real = 0.0117,
                               freelivfix_slope_wET::Real = 0.0006)
    params.freelivfix_intercept = freelivfix_intercept
    params.freelivfix_slope_wET = freelivfix_slope_wET
    return nothing
end

# ---------------------------------------------------------------------------
# n_deposition! — atmospheric N deposition
# Ported from CNNDeposition in CNNDynamicsMod.F90
# ---------------------------------------------------------------------------

# Per-column gather of atmospheric N deposition. cmin/cmax gate the active
# bounds range (one thread per column over 1:length(out)).
@kernel function _ndyn_deposition_kernel!(ndep_to_sminn_col, @Const(col_gridcell),
                                          @Const(forc_ndep), cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax
        ndep_to_sminn_col[c] = forc_ndep[col_gridcell[c]]
    end
end

ndyn_deposition!(ndep_to_sminn_col, col_gridcell, forc_ndep, cmin::Int, cmax::Int) =
    _launch!(_ndyn_deposition_kernel!, ndep_to_sminn_col, col_gridcell, forc_ndep, cmin, cmax)

"""
    n_deposition!(nf; forc_ndep, col_gridcell, bounds)

Update N deposition rate from atmospheric forcing. All atmospheric N
deposition goes to the soil mineral N pool.

Corresponds to `CNNDeposition` in the Fortran source.
"""
function n_deposition!(
    nf::SoilBiogeochemNitrogenFluxData;
    forc_ndep::AbstractVector{<:Real},
    col_gridcell::AbstractVector{<:Integer},
    bounds::UnitRange{Int})

    isempty(bounds) && return nothing
    ndyn_deposition!(nf.ndep_to_sminn_col, col_gridcell, forc_ndep,
                     first(bounds), last(bounds))

    return nothing
end

# ---------------------------------------------------------------------------
# n_free_living_fixation! — free-living N fixation as f(ET)
# Ported from CNFreeLivingFixation in CNNDynamicsMod.F90
# ---------------------------------------------------------------------------

# Per-column free-living N fixation from annual ET (masked). Pure arithmetic,
# fully independent per column. units: gN/m2/s.
@kernel function _ndyn_freelivfix_kernel!(ffix_to_sminn_col, @Const(mask), @Const(AnnET),
                                          freelivfix_slope, freelivfix_inter, secs_per_year)
    T = eltype(ffix_to_sminn_col)
    c = @index(Global)
    @inbounds if mask[c]
        ffix_to_sminn_col[c] =
            (freelivfix_slope * (max(zero(T), AnnET[c]) * secs_per_year) + freelivfix_inter) / secs_per_year
    end
end

ndyn_freelivfix!(ffix_to_sminn_col, mask, AnnET, freelivfix_slope, freelivfix_inter, secs_per_year) =
    _launch!(_ndyn_freelivfix_kernel!, ffix_to_sminn_col, mask, AnnET,
             eltype(ffix_to_sminn_col)(freelivfix_slope),
             eltype(ffix_to_sminn_col)(freelivfix_inter),
             eltype(ffix_to_sminn_col)(secs_per_year))

"""
    n_free_living_fixation!(nf, params; mask_soilc, bounds, AnnET, dayspyr)

Calculate free-living N fixation to soil mineral N as a function of
annual evapotranspiration (ET).

Corresponds to `CNFreeLivingFixation` in the Fortran source.
"""
function n_free_living_fixation!(
    nf::SoilBiogeochemNitrogenFluxData,
    params::NDynamicsParams;
    mask_soilc::AbstractVector{Bool},
    bounds::UnitRange{Int},
    AnnET::AbstractVector{<:Real},
    dayspyr::Real)

    secs_per_year = dayspyr * 24.0 * 3600.0

    freelivfix_slope = params.freelivfix_slope_wET
    freelivfix_inter = params.freelivfix_intercept

    ndyn_freelivfix!(nf.ffix_to_sminn_col, mask_soilc, AnnET,
                     freelivfix_slope, freelivfix_inter, secs_per_year)

    return nothing
end

# ---------------------------------------------------------------------------
# n_fixation! — symbiotic/asymbiotic N fixation as f(NPP)
# Ported from CNNFixation in CNNDynamicsMod.F90
# ---------------------------------------------------------------------------

# Per-column symbiotic N fixation using exponential relaxation with lagged NPP
# (masked, fully independent per column). col_is_fates[c] forces npp=0; SPVAL npp
# yields zero flux. secs_per_year = SECSPDAY * dayspyr.
@kernel function _ndyn_fixation_lag_kernel!(nfix_to_sminn_col, @Const(mask),
                                            @Const(col_lag_npp), @Const(col_is_fates),
                                            secs_per_year, spval)
    T = eltype(nfix_to_sminn_col)
    c = @index(Global)
    @inbounds if mask[c]
        npp = col_is_fates[c] ? zero(T) : col_lag_npp[c]
        if npp != spval
            t = (T(1.8) * (one(T) - exp(T(-0.003) * npp * secs_per_year))) / secs_per_year
            nfix_to_sminn_col[c] = max(zero(T), t)
        else
            nfix_to_sminn_col[c] = zero(T)
        end
    end
end

ndyn_fixation_lag!(nfix_to_sminn_col, mask, col_lag_npp, col_is_fates, secs_per_year, spval) =
    _launch!(_ndyn_fixation_lag_kernel!, nfix_to_sminn_col, mask, col_lag_npp,
             col_is_fates,
             eltype(nfix_to_sminn_col)(secs_per_year),
             eltype(nfix_to_sminn_col)(spval))

# Per-column symbiotic N fixation using annual-mean NPP (masked, independent).
@kernel function _ndyn_fixation_annsum_kernel!(nfix_to_sminn_col, @Const(mask),
                                               @Const(cannsum_npp), @Const(col_is_fates),
                                               secs_per_year)
    T = eltype(nfix_to_sminn_col)
    c = @index(Global)
    @inbounds if mask[c]
        npp = col_is_fates[c] ? zero(T) : cannsum_npp[c]
        t = (T(1.8) * (one(T) - exp(T(-0.003) * npp))) / secs_per_year
        nfix_to_sminn_col[c] = max(zero(T), t)
    end
end

ndyn_fixation_annsum!(nfix_to_sminn_col, mask, cannsum_npp, col_is_fates, secs_per_year) =
    _launch!(_ndyn_fixation_annsum_kernel!, nfix_to_sminn_col, mask, cannsum_npp,
             col_is_fates,
             eltype(nfix_to_sminn_col)(secs_per_year))

# Masked per-column zero (independent). Used to disable symbiotic fixation under
# FUN and to zero column-level accumulators before p2c aggregation.
@kernel function _ndyn_col_zero_kernel!(out, @Const(mask))
    T = eltype(out)
    c = @index(Global)
    @inbounds if mask[c]
        out[c] = zero(T)
    end
end

ndyn_col_zero!(out, mask) = _launch!(_ndyn_col_zero_kernel!, out, mask)

"""
    n_fixation!(nf, cf; mask_soilc, bounds, col_is_fates, dayspyr,
                nfix_timeconst, use_fun)

Calculate N fixation rate as a function of annual total NPP.
All N fixation goes to the soil mineral N pool.

When `nfix_timeconst ∈ (0, 500)`, uses exponential relaxation with
lagged NPP (`lag_npp_col`). Otherwise uses annual-mean NPP
(`annsum_npp_col`).

Corresponds to `CNNFixation` in the Fortran source.
"""
function n_fixation!(
    nf::SoilBiogeochemNitrogenFluxData,
    cf::CNVegCarbonFluxData;
    mask_soilc::AbstractVector{Bool},
    bounds::UnitRange{Int},
    col_is_fates::AbstractVector{Bool},
    dayspyr::Real,
    nfix_timeconst::Real,
    use_fun::Bool)

    cannsum_npp = cf.annsum_npp_col
    col_lag_npp = cf.lag_npp_col

    secs_per_year = SECSPDAY * dayspyr

    if nfix_timeconst > 0.0 && nfix_timeconst < 500.0
        # Use exponential relaxation with time constant nfix_timeconst for NPP - NFIX relation
        # (npp put in units of gC/m^2/year inside the kernel via secs_per_year)
        ndyn_fixation_lag!(nf.nfix_to_sminn_col, mask_soilc, col_lag_npp,
                           col_is_fates, secs_per_year, SPVAL)
    else
        # Use annual-mean values for NPP-NFIX relation
        ndyn_fixation_annsum!(nf.nfix_to_sminn_col, mask_soilc, cannsum_npp,
                              col_is_fates, secs_per_year)
    end

    # When FUN is active, disable symbiotic fixation
    if use_fun
        ndyn_col_zero!(nf.nfix_to_sminn_col, mask_soilc)
    end

    return nothing
end

# ---------------------------------------------------------------------------
# n_fert! — N fertilizer for crops (patch-to-column aggregation)
# Ported from CNNFert in CNNDynamicsMod.F90
# ---------------------------------------------------------------------------

"""
    n_fert!(soilbgc_nf, cnveg_nf; mask_soilc, bounds, patch, mask_soilp, bounds_p)

Aggregate patch-level fertilizer N flux to column level.
All fertilizer goes into the soil mineral N pool.

Corresponds to `CNNFert` in the Fortran source.
The Fortran `p2c` call is replaced by inline weighted averaging
using `patch.wtcol`.
"""
# Generic per-patch p2c weighted scatter (replaces Fortran p2c). Each active patch
# scatters patch_arr[p]*wtcol[p] into its column accumulator col_arr[column[p]].
# pmin/pmax gate the active bounds range (one thread per patch over 1:length(mask_p)).
# The accumulator MUST already be zeroed by the caller (ndyn_col_zero!). KA CPU
# iterates p ascending → byte-identical to the host += loop; on GPU the add is atomic.
@kernel function _ndyn_p2c_scatter_kernel!(col_arr, @Const(mask_p), @Const(pcol),
                                           @Const(patch_arr), @Const(wtcol),
                                           pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax && mask_p[p]
        _scatter_add!(col_arr, pcol[p], patch_arr[p] * wtcol[p])
    end
end

# ndrange is over PATCHES (length(mask_p)), not the column accumulator — the
# kernel index p must cover every patch, so pass ndrange explicitly.
ndyn_p2c_scatter!(col_arr, mask_p, pcol, patch_arr, wtcol, pmin::Int, pmax::Int) =
    _launch!(_ndyn_p2c_scatter_kernel!, col_arr, mask_p, pcol, patch_arr, wtcol, pmin, pmax;
             ndrange = length(mask_p))

function n_fert!(
    soilbgc_nf::SoilBiogeochemNitrogenFluxData,
    cnveg_nf::CNVegNitrogenFluxData;
    mask_soilc::AbstractVector{Bool},
    bounds::UnitRange{Int},
    patch::PatchData,
    mask_soilp::AbstractVector{Bool},
    bounds_p::UnitRange{Int})

    fert          = cnveg_nf.fert_patch
    fert_to_sminn = soilbgc_nf.fert_to_sminn_col

    # Zero out column-level accumulator (masked, independent per column)
    ndyn_col_zero!(fert_to_sminn, mask_soilc)

    # Patch-to-column weighted aggregation (replaces Fortran p2c call).
    # Per-patch scatter into the (already zeroed) column accumulator.
    isempty(bounds_p) || ndyn_p2c_scatter!(fert_to_sminn, mask_soilp, patch.column,
                                           fert, patch.wtcol, first(bounds_p), last(bounds_p))

    return nothing
end

# ---------------------------------------------------------------------------
# n_soyfix! — soybean N fixation
# Ported from CNSoyfix in CNNDynamicsMod.F90
# ---------------------------------------------------------------------------

# Per-patch soybean N fixation (EPICPHASE). Each patch writes only soyfixn[p];
# column-indexed fields (fpg, wf, sminn) are read-only gathers, so iterations are
# fully independent. Branch structure and thresholds mirror the Fortran exactly.
@kernel function _ndyn_soyfix_patch_kernel!(soyfixn, @Const(mask_p), @Const(pcol),
                                            @Const(itype), @Const(croplive),
                                            @Const(plant_ndemand), @Const(hui),
                                            @Const(gddmaturity), @Const(fpg),
                                            @Const(wf), @Const(sminn),
                                            ntmp_soybean::Int, nirrig_tmp_soybean::Int,
                                            ntrp_soybean::Int, nirrig_trp_soybean::Int)
    T = eltype(soyfixn)
    p = @index(Global)
    @inbounds if mask_p[p]
        c = pcol[p]

        sminnthreshold1   = T(30.0)
        sminnthreshold2   = T(10.0)
        GDDfracthreshold1 = T(0.15)
        GDDfracthreshold2 = T(0.30)
        GDDfracthreshold3 = T(0.55)
        GDDfracthreshold4 = T(0.75)

        if croplive[p] &&
           (itype[p] == ntmp_soybean ||
            itype[p] == nirrig_tmp_soybean ||
            itype[p] == ntrp_soybean ||
            itype[p] == nirrig_trp_soybean)

            if fpg[c] < one(T)
                soy_ndemand = plant_ndemand[p] - plant_ndemand[p] * fpg[c]

                # Soil water factor
                fxw = wf[c] / T(0.85)

                # Soil nitrogen factor
                fxn = zero(T)
                if sminn[c] > sminnthreshold1
                    fxn = zero(T)
                elseif sminn[c] > sminnthreshold2 && sminn[c] <= sminnthreshold1
                    fxn = T(1.5) - T(0.005) * (sminn[c] * T(10.0))
                elseif sminn[c] <= sminnthreshold2
                    fxn = one(T)
                end

                # Growth stage factor
                GDDfrac = hui[p] / gddmaturity[p]
                fxg = zero(T)
                if GDDfrac <= GDDfracthreshold1
                    fxg = zero(T)
                elseif GDDfrac > GDDfracthreshold1 && GDDfrac <= GDDfracthreshold2
                    fxg = T(6.67) * GDDfrac - one(T)
                elseif GDDfrac > GDDfracthreshold2 && GDDfrac <= GDDfracthreshold3
                    fxg = one(T)
                elseif GDDfrac > GDDfracthreshold3 && GDDfrac <= GDDfracthreshold4
                    fxg = T(3.75) - T(5.0) * GDDfrac
                else
                    fxg = zero(T)
                end

                # Calculate the nitrogen fixed by the soybean
                fxr = min(one(T), fxw, fxn) * fxg
                fxr = max(zero(T), fxr)
                s = fxr * soy_ndemand
                soyfixn[p] = min(s, soy_ndemand)
            else
                soyfixn[p] = zero(T)
            end
        else
            soyfixn[p] = zero(T)
        end
    end
end

function ndyn_soyfix_patch!(soyfixn, mask_p, pcol, itype, croplive, plant_ndemand,
                            hui, gddmaturity, fpg, wf, sminn,
                            ntmp_soybean, nirrig_tmp_soybean, ntrp_soybean, nirrig_trp_soybean)
    _launch!(_ndyn_soyfix_patch_kernel!, soyfixn, mask_p, pcol, itype, croplive,
             plant_ndemand, hui, gddmaturity, fpg, wf, sminn,
             ntmp_soybean, nirrig_tmp_soybean, ntrp_soybean, nirrig_trp_soybean)
end

"""
    n_soyfix!(soilbgc_nf, cnveg_nf, soilbgc_state, soilbgc_ns, cnveg_state, crop, wdiag;
              mask_soilc, bounds, mask_soilp, bounds_p, patch,
              ntmp_soybean, nirrig_tmp_soybean, ntrp_soybean, nirrig_trp_soybean)

Calculate N fixation for soybeans based on the EPICPHASE model
(Cabelguenne et al., Agricultural Systems 60: 175-196, 1999).
N fixation depends on soil moisture, plant growth phase, and
availability of nitrogen in the soil root zone.

Corresponds to `CNSoyfix` in the Fortran source.
"""
function n_soyfix!(
    soilbgc_nf::SoilBiogeochemNitrogenFluxData,
    cnveg_nf::CNVegNitrogenFluxData,
    soilbgc_state::SoilBiogeochemStateData,
    soilbgc_ns::SoilBiogeochemNitrogenStateData,
    cnveg_state::CNVegStateData,
    crop::CropData,
    wdiag::WaterDiagnosticBulkData;
    mask_soilc::AbstractVector{Bool},
    bounds::UnitRange{Int},
    mask_soilp::AbstractVector{Bool},
    bounds_p::UnitRange{Int},
    patch::PatchData,
    ntmp_soybean::Int,
    nirrig_tmp_soybean::Int,
    ntrp_soybean::Int,
    nirrig_trp_soybean::Int)

    # Aliases
    wf            = wdiag.wf_col
    hui           = crop.hui_patch
    croplive      = crop.croplive_patch
    gddmaturity   = cnveg_state.gddmaturity_patch
    plant_ndemand = cnveg_nf.plant_ndemand_patch
    soyfixn       = cnveg_nf.soyfixn_patch
    fpg           = soilbgc_state.fpg_col
    sminn         = soilbgc_ns.sminn_col
    soyfixn_to_sminn = soilbgc_nf.soyfixn_to_sminn_col

    # Per-patch soybean fixation (independent per patch; thresholds inside kernel)
    ndyn_soyfix_patch!(soyfixn, mask_soilp, patch.column, patch.itype, croplive,
                       plant_ndemand, hui, gddmaturity, fpg, wf, sminn,
                       ntmp_soybean, nirrig_tmp_soybean, ntrp_soybean, nirrig_trp_soybean)

    # Patch-to-column aggregation (replaces Fortran p2c call)
    ndyn_col_zero!(soyfixn_to_sminn, mask_soilc)

    # Per-patch scatter into the (already zeroed) column accumulator.
    isempty(bounds_p) || ndyn_p2c_scatter!(soyfixn_to_sminn, mask_soilp, patch.column,
                                           soyfixn, patch.wtcol, first(bounds_p), last(bounds_p))

    return nothing
end
