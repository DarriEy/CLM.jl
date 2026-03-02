# ==========================================================================
# Ported from: src/biogeophys/SoilMoistStressMod.F90
# Calculates soil moisture stress for plant GPP and transpiration.
# ==========================================================================

# --- Module-level control state ---

const MOIST_STRESS_CLM_DEFAULT = 0

Base.@kwdef mutable struct SoilMoistStressControl
    root_moist_stress_method::Int = MOIST_STRESS_CLM_DEFAULT
    perchroot::Bool = false       # true => btran is based only on unfrozen soil levels
    perchroot_alt::Bool = false   # true => btran is based on active layer (defined over two years)
end

const soil_moist_stress_ctrl = SoilMoistStressControl()

"""
    init_root_moist_stress!()

Specify the method to compute root soil moisture stress.
Ported from `init_root_moist_stress` in `SoilMoistStressMod.F90`.
"""
function init_root_moist_stress!()
    soil_moist_stress_ctrl.root_moist_stress_method = MOIST_STRESS_CLM_DEFAULT
    return nothing
end

"""
    set_perchroot_opt!(perchroot_global, perchroot_alt_global)

Set up local perchroot logical switches.
Ported from `set_perchroot_opt` in `SoilMoistStressMod.F90`.
"""
function set_perchroot_opt!(perchroot_global::Bool, perchroot_alt_global::Bool)
    soil_moist_stress_ctrl.perchroot = perchroot_global
    soil_moist_stress_ctrl.perchroot_alt = perchroot_alt_global
    return nothing
end

# --------------------------------------------------------------------------
# calc_effective_soilporosity!
# --------------------------------------------------------------------------

"""
    calc_effective_soilporosity!(watsat, h2osoi_ice, col_dz, eff_por,
                                 mask, bounds_col, nlevgrnd, nlevsno)

Compute the effective soil porosity.
Ported from `calc_effective_soilporosity` in `SoilMoistStressMod.F90`.

Arguments:
- `watsat`: soil porosity (nc × nlevgrnd), soil-only indexing
- `h2osoi_ice`: ice water content kg/m2 (nc × nlevsno+nlevgrnd), combined indexing
- `col_dz`: layer thickness (nc × nlevsno+nlevgrnd), combined indexing
- `eff_por`: effective porosity output (nc × nlevgrnd), soil-only indexing
- `mask`: BitVector mask for active columns
- `bounds_col`: column index range
- `nlevgrnd`: number of ground levels (ubj)
- `nlevsno`: number of snow levels (for index offset)
"""
function calc_effective_soilporosity!(watsat::Matrix{Float64},
                                      h2osoi_ice::Matrix{Float64},
                                      col_dz::Matrix{Float64},
                                      eff_por::Matrix{Float64},
                                      mask::BitVector,
                                      bounds_col::UnitRange{Int},
                                      nlevgrnd::Int,
                                      nlevsno::Int)
    joff = nlevsno  # offset for combined snow+soil indexing

    for j in 1:nlevgrnd
        for c in bounds_col
            mask[c] || continue
            # compute the volumetric ice content
            vol_ice = min(watsat[c, j], h2osoi_ice[c, j + joff] / (DENICE * col_dz[c, j + joff]))
            # compute the maximum soil space to fill liquid water and air
            eff_por[c, j] = watsat[c, j] - vol_ice
        end
    end
    return nothing
end

# --------------------------------------------------------------------------
# calc_effective_snowporosity!
# --------------------------------------------------------------------------

"""
    calc_effective_snowporosity!(h2osoi_ice, col_dz, jtop, eff_por,
                                 mask, bounds_col, lbj, nlevsno)

Compute the effective porosity of snow.
Ported from `calc_effective_snowporosity` in `SoilMoistStressMod.F90`.

Arguments:
- `h2osoi_ice`: ice water content kg/m2 (nc × nlevsno+nlevgrnd), combined indexing
- `col_dz`: layer thickness (nc × nlevsno+nlevgrnd), combined indexing
- `jtop`: top level for each column (Fortran snow layer index, e.g. -4 to 0)
- `eff_por`: effective porosity output (nc × nlevsno), snow-only indexing (1:nlevsno maps to lbj:0)
- `mask`: BitVector mask for active columns
- `bounds_col`: column index range
- `lbj`: lower bound snow layer index (Fortran, e.g. -nlevsno+1)
- `nlevsno`: number of snow levels
"""
function calc_effective_snowporosity!(h2osoi_ice::Matrix{Float64},
                                      col_dz::Matrix{Float64},
                                      jtop::Vector{Int},
                                      eff_por::Matrix{Float64},
                                      mask::BitVector,
                                      bounds_col::UnitRange{Int},
                                      lbj::Int,
                                      nlevsno::Int)
    joff = nlevsno  # offset: Fortran j → Julia index j + nlevsno

    # Loop over snow layers: Fortran j from lbj to 0
    for j in lbj:0
        jj = j + joff  # Julia 1-based index
        for c in bounds_col
            mask[c] || continue
            if j >= jtop[c]
                # compute the volumetric ice content
                vol_ice = min(1.0, h2osoi_ice[c, jj] / (DENICE * col_dz[c, jj]))
                # compute the maximum snow void space to fill liquid water and air
                eff_por[c, jj] = 1.0 - vol_ice
            end
        end
    end
    return nothing
end

# --------------------------------------------------------------------------
# calc_volumetric_h2oliq!
# --------------------------------------------------------------------------

"""
    calc_volumetric_h2oliq!(eff_porosity, h2osoi_liq, col_dz, jtop, vol_liq,
                             mask, bounds_col, lbj, ubj, nlevsno)

Compute the volumetric liquid water content.
Ported from `calc_volumetric_h2oliq` in `SoilMoistStressMod.F90`.

Arguments:
- `eff_porosity`: effective soil porosity (nc × dim), indexing depends on lbj
- `h2osoi_liq`: liquid water content kg/m2 (nc × nlevsno+nlevgrnd), combined indexing
- `col_dz`: layer thickness (nc × nlevsno+nlevgrnd), combined indexing
- `jtop`: top level for each column (Fortran snow layer index)
- `vol_liq`: volumetric liquid water content output (nc × dim)
- `mask`: BitVector mask for active columns
- `bounds_col`: column index range
- `lbj`: lower bound layer index (Fortran)
- `ubj`: upper bound layer index (Fortran)
- `nlevsno`: number of snow levels (for index offset)
"""
function calc_volumetric_h2oliq!(eff_porosity::Matrix{Float64},
                                  h2osoi_liq::Matrix{Float64},
                                  col_dz::Matrix{Float64},
                                  jtop::Vector{Int},
                                  vol_liq::Matrix{Float64},
                                  mask::BitVector,
                                  bounds_col::UnitRange{Int},
                                  lbj::Int,
                                  ubj::Int,
                                  nlevsno::Int)
    joff = nlevsno  # offset: Fortran j → Julia index j + nlevsno

    for j in lbj:ubj
        jj = j + joff  # Julia 1-based index into combined arrays
        for c in bounds_col
            mask[c] || continue
            if j >= jtop[c]
                # volume of liquid is no greater than effective void space
                vol_liq[c, jj] = min(eff_porosity[c, jj], h2osoi_liq[c, jj] / (col_dz[c, jj] * DENH2O))
            end
        end
    end
    return nothing
end

# --------------------------------------------------------------------------
# array_normalization! (from SimpleMathMod)
# --------------------------------------------------------------------------

"""
    array_normalization!(arr, mask, bounds, nlevgrnd)

Normalize each row of `arr` (patch × nlevgrnd) so that its sum equals 1.
If the sum is zero, the row is left as zeros.
Ported from `array_normalization` in `SimpleMathMod.F90`.
"""
function array_normalization!(arr::Matrix{Float64},
                               mask::BitVector,
                               bounds::UnitRange{Int},
                               nlevgrnd::Int)
    for p in bounds
        mask[p] || continue
        rootsum = 0.0
        for j in 1:nlevgrnd
            rootsum += arr[p, j]
        end
        if rootsum > 0.0
            for j in 1:nlevgrnd
                arr[p, j] /= rootsum
            end
        end
    end
    return nothing
end

# --------------------------------------------------------------------------
# normalize_unfrozen_rootfr!
# --------------------------------------------------------------------------

"""
    normalize_unfrozen_rootfr!(rootfr, t_soisno, altmax_lastyear_indx, altmax_indx,
                                patch_column, rootfr_unf,
                                mask_patch, bounds_patch, nlevgrnd, nlevsno)

Normalize root fraction for total unfrozen depth.
Ported from `normalize_unfrozen_rootfr` in `SoilMoistStressMod.F90`.
"""
function normalize_unfrozen_rootfr!(rootfr::Matrix{Float64},
                                     t_soisno::Matrix{Float64},
                                     altmax_lastyear_indx::Vector{Float64},
                                     altmax_indx::Vector{Float64},
                                     patch_column::Vector{Int},
                                     rootfr_unf::Matrix{Float64},
                                     mask_patch::BitVector,
                                     bounds_patch::UnitRange{Int},
                                     nlevgrnd::Int,
                                     nlevsno::Int)
    joff = nlevsno
    perchroot = soil_moist_stress_ctrl.perchroot
    perchroot_alt = soil_moist_stress_ctrl.perchroot_alt

    if perchroot || perchroot_alt
        if perchroot_alt
            # use total active layer (max thaw depth for current and prior year)
            for j in 1:nlevgrnd
                for p in bounds_patch
                    mask_patch[p] || continue
                    c = patch_column[p]
                    if j <= max(altmax_lastyear_indx[c], altmax_indx[c], 1.0)
                        rootfr_unf[p, j] = rootfr[p, j]
                    else
                        rootfr_unf[p, j] = 0.0
                    end
                end
            end
        else
            # use instantaneous temperature
            for j in 1:nlevgrnd
                for p in bounds_patch
                    mask_patch[p] || continue
                    c = patch_column[p]
                    if t_soisno[c, j + joff] >= TFRZ
                        rootfr_unf[p, j] = rootfr[p, j]
                    else
                        rootfr_unf[p, j] = 0.0
                    end
                end
            end
        end
    end

    # normalize the root fraction for each pft
    array_normalization!(rootfr_unf, mask_patch, bounds_patch, nlevgrnd)

    return nothing
end

# --------------------------------------------------------------------------
# soil_suction (inline replacement for SoilWaterRetentionCurveMod)
# --------------------------------------------------------------------------

"""
    soil_suction_clapp_hornberger(sucsat, s_node, bsw)

Compute soil matric potential using the Clapp-Hornberger parameterization.
Returns smp_node (matric potential in mm).
Ported from `soil_suction` in `SoilWaterRetentionCurveMod.F90`.
"""
function soil_suction_clapp_hornberger(sucsat::Float64, s_node::Float64, bsw::Float64)
    return -sucsat * s_node^(-bsw)
end

# --------------------------------------------------------------------------
# calc_root_moist_stress_clm45default!
# --------------------------------------------------------------------------

"""
    calc_root_moist_stress_clm45default!(...)

Compute root water stress using the default CLM4.5 approach.
Ported from `calc_root_moist_stress_clm45default` in `SoilMoistStressMod.F90`.
"""
function calc_root_moist_stress_clm45default!(rootfr_unf::Matrix{Float64},
                                               rootfr::Matrix{Float64},
                                               rootr::Matrix{Float64},
                                               btran::Vector{Float64},
                                               rresis::Matrix{Float64},
                                               smpso::Vector{Float64},
                                               smpsc::Vector{Float64},
                                               t_soisno::Matrix{Float64},
                                               watsat::Matrix{Float64},
                                               sucsat::Matrix{Float64},
                                               bsw::Matrix{Float64},
                                               eff_porosity::Matrix{Float64},
                                               h2osoi_liqvol::Matrix{Float64},
                                               patch_column::Vector{Int},
                                               patch_itype::Vector{Int},
                                               mask_patch::BitVector,
                                               bounds_patch::UnitRange{Int},
                                               nlevgrnd::Int,
                                               nlevsno::Int)
    btran0 = 0.0
    joff = nlevsno
    perchroot = soil_moist_stress_ctrl.perchroot
    perchroot_alt = soil_moist_stress_ctrl.perchroot_alt

    for j in 1:nlevgrnd
        for p in bounds_patch
            mask_patch[p] || continue
            c = patch_column[p]
            itype = patch_itype[p]

            # Root resistance factors
            if h2osoi_liqvol[c, j + joff] <= 0.0 || t_soisno[c, j + joff] <= TFRZ - 2.0
                rootr[p, j] = 0.0
            else
                s_node = max(h2osoi_liqvol[c, j + joff] / eff_porosity[c, j], 0.01)

                smp_node = soil_suction_clapp_hornberger(sucsat[c, j], s_node, bsw[c, j])
                smp_node = max(smpsc[itype], smp_node)

                rresis[p, j] = min((eff_porosity[c, j] / watsat[c, j]) *
                    (smp_node - smpsc[itype]) / (smpso[itype] - smpsc[itype]), 1.0)

                if !(perchroot || perchroot_alt)
                    rootr[p, j] = rootfr[p, j] * rresis[p, j]
                else
                    rootr[p, j] = rootfr_unf[p, j] * rresis[p, j]
                end

                btran[p] = btran[p] + max(rootr[p, j], 0.0)
            end
        end
    end

    # Normalize root resistances to get layer contribution to ET
    for j in 1:nlevgrnd
        for p in bounds_patch
            mask_patch[p] || continue
            if btran[p] > btran0
                rootr[p, j] = rootr[p, j] / btran[p]
            else
                rootr[p, j] = 0.0
            end
        end
    end

    return nothing
end

# --------------------------------------------------------------------------
# calc_root_moist_stress! (top-level dispatcher)
# --------------------------------------------------------------------------

"""
    calc_root_moist_stress!(...)

Compute root water stress using different approaches.
Ported from `calc_root_moist_stress` in `SoilMoistStressMod.F90`.
"""
function calc_root_moist_stress!(soilstate::SoilStateData,
                                  energyflux::EnergyFluxData,
                                  temperature::TemperatureData,
                                  waterstatebulk::WaterStateBulkData,
                                  waterdiagbulk::WaterDiagnosticBulkData,
                                  col::ColumnData,
                                  patchdata::PatchData,
                                  smpso::Vector{Float64},
                                  smpsc::Vector{Float64},
                                  altmax_lastyear_indx::Vector{Float64},
                                  altmax_indx::Vector{Float64},
                                  mask_patch::BitVector,
                                  bounds_patch::UnitRange{Int},
                                  nlevgrnd::Int,
                                  nlevsno::Int)

    np = length(bounds_patch)

    # Initialize rootfr_unf to zero
    rootfr_unf = zeros(Float64, length(mask_patch), nlevgrnd)

    # Initialize btran to zero for accumulation
    for p in bounds_patch
        mask_patch[p] || continue
        energyflux.btran_patch[p] = 0.0
    end

    # Define normalized rootfraction for unfrozen soil
    normalize_unfrozen_rootfr!(soilstate.rootfr_patch,
                                temperature.t_soisno_col,
                                altmax_lastyear_indx,
                                altmax_indx,
                                patchdata.column,
                                rootfr_unf,
                                mask_patch,
                                bounds_patch,
                                nlevgrnd,
                                nlevsno)

    method = soil_moist_stress_ctrl.root_moist_stress_method

    if method == MOIST_STRESS_CLM_DEFAULT
        calc_root_moist_stress_clm45default!(rootfr_unf,
                                              soilstate.rootfr_patch,
                                              soilstate.rootr_patch,
                                              energyflux.btran_patch,
                                              energyflux.rresis_patch,
                                              smpso,
                                              smpsc,
                                              temperature.t_soisno_col,
                                              soilstate.watsat_col,
                                              soilstate.sucsat_col,
                                              soilstate.bsw_col,
                                              soilstate.eff_porosity_col,
                                              waterdiagbulk.h2osoi_liqvol_col,
                                              patchdata.column,
                                              patchdata.itype,
                                              mask_patch,
                                              bounds_patch,
                                              nlevgrnd,
                                              nlevsno)
    else
        error("calc_root_moist_stress: a root moisture stress function must be specified!")
    end

    return nothing
end
