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

"""
    n_deposition!(nf; forc_ndep, col_gridcell, bounds)

Update N deposition rate from atmospheric forcing. All atmospheric N
deposition goes to the soil mineral N pool.

Corresponds to `CNNDeposition` in the Fortran source.
"""
function n_deposition!(
    nf::SoilBiogeochemNitrogenFluxData;
    forc_ndep::Vector{<:Real},
    col_gridcell::Vector{Int},
    bounds::UnitRange{Int})

    for c in bounds
        g = col_gridcell[c]
        nf.ndep_to_sminn_col[c] = forc_ndep[g]
    end

    return nothing
end

# ---------------------------------------------------------------------------
# n_free_living_fixation! — free-living N fixation as f(ET)
# Ported from CNFreeLivingFixation in CNNDynamicsMod.F90
# ---------------------------------------------------------------------------

"""
    n_free_living_fixation!(nf, params; mask_soilc, bounds, AnnET, dayspyr)

Calculate free-living N fixation to soil mineral N as a function of
annual evapotranspiration (ET).

Corresponds to `CNFreeLivingFixation` in the Fortran source.
"""
function n_free_living_fixation!(
    nf::SoilBiogeochemNitrogenFluxData,
    params::NDynamicsParams;
    mask_soilc::BitVector,
    bounds::UnitRange{Int},
    AnnET::Vector{<:Real},
    dayspyr::Real)

    secs_per_year = dayspyr * 24.0 * 3600.0

    freelivfix_slope = params.freelivfix_slope_wET
    freelivfix_inter = params.freelivfix_intercept

    for c in bounds
        mask_soilc[c] || continue
        # units: gN/m2/s
        nf.ffix_to_sminn_col[c] = (freelivfix_slope * (max(0.0, AnnET[c]) * secs_per_year) + freelivfix_inter) / secs_per_year
    end

    return nothing
end

# ---------------------------------------------------------------------------
# n_fixation! — symbiotic/asymbiotic N fixation as f(NPP)
# Ported from CNNFixation in CNNDynamicsMod.F90
# ---------------------------------------------------------------------------

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
    mask_soilc::BitVector,
    bounds::UnitRange{Int},
    col_is_fates::Vector{Bool},
    dayspyr::Real,
    nfix_timeconst::Real,
    use_fun::Bool)

    cannsum_npp = cf.annsum_npp_col
    col_lag_npp = cf.lag_npp_col

    if nfix_timeconst > 0.0 && nfix_timeconst < 500.0
        # Use exponential relaxation with time constant nfix_timeconst for NPP - NFIX relation
        for c in bounds
            mask_soilc[c] || continue

            if col_is_fates[c]
                # FATES N cycling not yet active; set npp to 0
                npp = 0.0
            else
                npp = col_lag_npp[c]
            end

            if npp != SPVAL
                # Need to put npp in units of gC/m^2/year here first
                t = (1.8 * (1.0 - exp(-0.003 * npp * (SECSPDAY * dayspyr)))) / (SECSPDAY * dayspyr)
                nf.nfix_to_sminn_col[c] = max(0.0, t)
            else
                nf.nfix_to_sminn_col[c] = 0.0
            end
        end
    else
        # Use annual-mean values for NPP-NFIX relation
        for c in bounds
            mask_soilc[c] || continue

            if col_is_fates[c]
                npp = 0.0
            else
                npp = cannsum_npp[c]
            end

            t = (1.8 * (1.0 - exp(-0.003 * npp))) / (SECSPDAY * dayspyr)
            nf.nfix_to_sminn_col[c] = max(0.0, t)
        end
    end

    # When FUN is active, disable symbiotic fixation
    if use_fun
        for c in bounds
            mask_soilc[c] || continue
            nf.nfix_to_sminn_col[c] = 0.0
        end
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
function n_fert!(
    soilbgc_nf::SoilBiogeochemNitrogenFluxData,
    cnveg_nf::CNVegNitrogenFluxData;
    mask_soilc::BitVector,
    bounds::UnitRange{Int},
    patch::PatchData,
    mask_soilp::BitVector,
    bounds_p::UnitRange{Int})

    fert          = cnveg_nf.fert_patch
    fert_to_sminn = soilbgc_nf.fert_to_sminn_col

    # Zero out column-level accumulator
    for c in bounds
        mask_soilc[c] || continue
        fert_to_sminn[c] = 0.0
    end

    # Patch-to-column weighted aggregation (replaces Fortran p2c call)
    for p in bounds_p
        mask_soilp[p] || continue
        c = patch.column[p]
        fert_to_sminn[c] += fert[p] * patch.wtcol[p]
    end

    return nothing
end

# ---------------------------------------------------------------------------
# n_soyfix! — soybean N fixation
# Ported from CNSoyfix in CNNDynamicsMod.F90
# ---------------------------------------------------------------------------

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
    mask_soilc::BitVector,
    bounds::UnitRange{Int},
    mask_soilp::BitVector,
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

    # Thresholds (same as Fortran)
    sminnthreshold1    = 30.0
    sminnthreshold2    = 10.0
    GDDfracthreshold1  = 0.15
    GDDfracthreshold2  = 0.30
    GDDfracthreshold3  = 0.55
    GDDfracthreshold4  = 0.75

    for p in bounds_p
        mask_soilp[p] || continue
        c = patch.column[p]

        # If soybean currently growing then calculate fixation
        if croplive[p] &&
            (patch.itype[p] == ntmp_soybean ||
             patch.itype[p] == nirrig_tmp_soybean ||
             patch.itype[p] == ntrp_soybean ||
             patch.itype[p] == nirrig_trp_soybean)

            # Difference between supply and demand
            if fpg[c] < 1.0
                soy_ndemand = 0.0
                soy_ndemand = plant_ndemand[p] - plant_ndemand[p] * fpg[c]

                # Soil water factor
                fxw = 0.0
                fxw = wf[c] / 0.85

                # Soil nitrogen factor (CHECK UNITS)
                if sminn[c] > sminnthreshold1
                    fxn = 0.0
                elseif sminn[c] > sminnthreshold2 && sminn[c] <= sminnthreshold1
                    fxn = 1.5 - 0.005 * (sminn[c] * 10.0)
                elseif sminn[c] <= sminnthreshold2
                    fxn = 1.0
                end

                # Growth stage factor
                GDDfrac = hui[p] / gddmaturity[p]

                if GDDfrac <= GDDfracthreshold1
                    fxg = 0.0
                elseif GDDfrac > GDDfracthreshold1 && GDDfrac <= GDDfracthreshold2
                    fxg = 6.67 * GDDfrac - 1.0
                elseif GDDfrac > GDDfracthreshold2 && GDDfrac <= GDDfracthreshold3
                    fxg = 1.0
                elseif GDDfrac > GDDfracthreshold3 && GDDfrac <= GDDfracthreshold4
                    fxg = 3.75 - 5.0 * GDDfrac
                else  # GDDfrac > GDDfracthreshold4
                    fxg = 0.0
                end

                # Calculate the nitrogen fixed by the soybean
                fxr = min(1.0, fxw, fxn) * fxg
                fxr = max(0.0, fxr)
                soyfixn[p] = fxr * soy_ndemand
                soyfixn[p] = min(soyfixn[p], soy_ndemand)

            else  # Nitrogen demand met, no fixation
                soyfixn[p] = 0.0
            end

        else  # Not live soybean, no fixation
            soyfixn[p] = 0.0
        end
    end

    # Patch-to-column aggregation (replaces Fortran p2c call)
    for c in bounds
        mask_soilc[c] || continue
        soyfixn_to_sminn[c] = 0.0
    end

    for p in bounds_p
        mask_soilp[p] || continue
        c = patch.column[p]
        soyfixn_to_sminn[c] += soyfixn[p] * patch.wtcol[p]
    end

    return nothing
end
