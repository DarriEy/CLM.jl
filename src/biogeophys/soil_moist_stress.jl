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
# --- Kernel: effective soil porosity over (column, soil layer) ---
@kernel function _smstress_effsoilpor_kernel!(eff_por, @Const(mask), @Const(watsat),
                                              @Const(h2osoi_ice), @Const(col_dz),
                                              joff::Int, denice)
    c, j = @index(Global, NTuple)
    @inbounds if mask[c]
        # compute the volumetric ice content
        vol_ice = min(watsat[c, j], h2osoi_ice[c, j + joff] / (denice * col_dz[c, j + joff]))
        # compute the maximum soil space to fill liquid water and air
        eff_por[c, j] = watsat[c, j] - vol_ice
    end
end

smstress_effsoilpor!(eff_por, mask, watsat, h2osoi_ice, col_dz, joff::Int, nlevgrnd::Int, denice) =
    _launch!(_smstress_effsoilpor_kernel!, eff_por, mask, watsat, h2osoi_ice, col_dz,
             joff, denice; ndrange = (length(mask), nlevgrnd))

function calc_effective_soilporosity!(watsat::Matrix{<:Real},
                                      h2osoi_ice::Matrix{<:Real},
                                      col_dz::Matrix{<:Real},
                                      eff_por::Matrix{<:Real},
                                      mask::BitVector,
                                      bounds_col::UnitRange{Int},
                                      nlevgrnd::Int,
                                      nlevsno::Int)
    joff = nlevsno  # offset for combined snow+soil indexing
    smstress_effsoilpor!(eff_por, mask, watsat, h2osoi_ice, col_dz, joff, nlevgrnd, DENICE)
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
# --- Kernel: effective snow porosity over (column, snow layer) ---
# jj is the Julia 1-based combined index (1..nlevsno); the Fortran snow index is
# j = jj - joff. Only layers at or below the column top (j >= jtop) are written.
@kernel function _smstress_effsnowpor_kernel!(eff_por, @Const(mask), @Const(jtop),
                                              @Const(h2osoi_ice), @Const(col_dz),
                                              joff::Int, denice)
    c, jj = @index(Global, NTuple)
    @inbounds if mask[c]
        j = jj - joff  # Fortran snow layer index (negative)
        if j >= jtop[c]
            # compute the volumetric ice content
            vol_ice = min(1.0, h2osoi_ice[c, jj] / (denice * col_dz[c, jj]))
            # compute the maximum snow void space to fill liquid water and air
            eff_por[c, jj] = 1.0 - vol_ice
        end
    end
end

smstress_effsnowpor!(eff_por, mask, jtop, h2osoi_ice, col_dz, joff::Int, nlevsno::Int, denice) =
    _launch!(_smstress_effsnowpor_kernel!, eff_por, mask, jtop, h2osoi_ice, col_dz,
             joff, denice; ndrange = (length(mask), nlevsno))

function calc_effective_snowporosity!(h2osoi_ice::Matrix{<:Real},
                                      col_dz::Matrix{<:Real},
                                      jtop::Vector{Int},
                                      eff_por::Matrix{<:Real},
                                      mask::BitVector,
                                      bounds_col::UnitRange{Int},
                                      lbj::Int,
                                      nlevsno::Int)
    joff = nlevsno  # offset: Fortran j → Julia index j + nlevsno
    # Fortran snow layers run lbj:0, i.e. Julia combined indices 1:nlevsno.
    smstress_effsnowpor!(eff_por, mask, jtop, h2osoi_ice, col_dz, joff, nlevsno, DENICE)
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
# --- Kernel: volumetric liquid water content over (column, layer) ---
# jrel is 1-based over the layer span lbj:ubj; the combined Julia index is
# jj = jrel + lbj - 1 + joff and the Fortran layer index is j = jj - joff.
# Only layers at or below the column top (j >= jtop) are written.
@kernel function _smstress_volh2oliq_kernel!(vol_liq, @Const(mask), @Const(jtop),
                                             @Const(eff_porosity), @Const(h2osoi_liq),
                                             @Const(col_dz), lbj::Int, joff::Int, denh2o)
    c, jrel = @index(Global, NTuple)
    @inbounds if mask[c]
        j = lbj + (jrel - 1)            # Fortran layer index
        jj = j + joff                   # Julia 1-based combined index
        if j >= jtop[c]
            # volume of liquid is no greater than effective void space
            vol_liq[c, jj] = min(eff_porosity[c, jj], h2osoi_liq[c, jj] / (col_dz[c, jj] * denh2o))
        end
    end
end

function smstress_volh2oliq!(vol_liq, mask, jtop, eff_porosity, h2osoi_liq, col_dz,
                             lbj::Int, ubj::Int, joff::Int, denh2o)
    nlay = ubj - lbj + 1
    _launch!(_smstress_volh2oliq_kernel!, vol_liq, mask, jtop, eff_porosity, h2osoi_liq,
             col_dz, lbj, joff, denh2o; ndrange = (length(mask), nlay))
end

function calc_volumetric_h2oliq!(eff_porosity::Matrix{<:Real},
                                  h2osoi_liq::Matrix{<:Real},
                                  col_dz::Matrix{<:Real},
                                  jtop::Vector{Int},
                                  vol_liq::Matrix{<:Real},
                                  mask::BitVector,
                                  bounds_col::UnitRange{Int},
                                  lbj::Int,
                                  ubj::Int,
                                  nlevsno::Int)
    joff = nlevsno  # offset: Fortran j → Julia index j + nlevsno
    smstress_volh2oliq!(vol_liq, mask, jtop, eff_porosity, h2osoi_liq, col_dz,
                        lbj, ubj, joff, DENH2O)
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
function array_normalization!(arr::Matrix{<:Real},
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
# --- Kernel: unfrozen rootfr from active layer depth, over (patch, layer) ---
@kernel function _smstress_rootfr_unf_alt_kernel!(rootfr_unf, @Const(mask_patch),
                                                  @Const(patch_column), @Const(rootfr),
                                                  @Const(altmax_lastyear_indx),
                                                  @Const(altmax_indx))
    p, j = @index(Global, NTuple)
    @inbounds if mask_patch[p]
        c = patch_column[p]
        if j <= max(altmax_lastyear_indx[c], altmax_indx[c], 1.0)
            rootfr_unf[p, j] = rootfr[p, j]
        else
            rootfr_unf[p, j] = 0.0
        end
    end
end

smstress_rootfr_unf_alt!(rootfr_unf, mask_patch, patch_column, rootfr,
                         altmax_lastyear_indx, altmax_indx, nlevgrnd::Int) =
    _launch!(_smstress_rootfr_unf_alt_kernel!, rootfr_unf, mask_patch, patch_column,
             rootfr, altmax_lastyear_indx, altmax_indx;
             ndrange = (length(mask_patch), nlevgrnd))

# --- Kernel: unfrozen rootfr from instantaneous temperature, over (patch, layer) ---
@kernel function _smstress_rootfr_unf_temp_kernel!(rootfr_unf, @Const(mask_patch),
                                                   @Const(patch_column), @Const(rootfr),
                                                   @Const(t_soisno), joff::Int, tfrz)
    p, j = @index(Global, NTuple)
    @inbounds if mask_patch[p]
        c = patch_column[p]
        if t_soisno[c, j + joff] >= tfrz
            rootfr_unf[p, j] = rootfr[p, j]
        else
            rootfr_unf[p, j] = 0.0
        end
    end
end

smstress_rootfr_unf_temp!(rootfr_unf, mask_patch, patch_column, rootfr, t_soisno,
                          joff::Int, nlevgrnd::Int, tfrz) =
    _launch!(_smstress_rootfr_unf_temp_kernel!, rootfr_unf, mask_patch, patch_column,
             rootfr, t_soisno, joff, tfrz; ndrange = (length(mask_patch), nlevgrnd))

function normalize_unfrozen_rootfr!(rootfr::Matrix{<:Real},
                                     t_soisno::Matrix{<:Real},
                                     altmax_lastyear_indx::Vector{<:Real},
                                     altmax_indx::Vector{<:Real},
                                     patch_column::Vector{Int},
                                     rootfr_unf::Matrix{<:Real},
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
            smstress_rootfr_unf_alt!(rootfr_unf, mask_patch, patch_column, rootfr,
                                     altmax_lastyear_indx, altmax_indx, nlevgrnd)
        else
            # use instantaneous temperature
            smstress_rootfr_unf_temp!(rootfr_unf, mask_patch, patch_column, rootfr,
                                      t_soisno, joff, nlevgrnd, TFRZ)
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
function soil_suction_clapp_hornberger(sucsat::Real, s_node::Real, bsw::Real)
    return -sucsat * s_node^(-bsw)
end

# --------------------------------------------------------------------------
# calc_root_moist_stress_clm45default!
# --------------------------------------------------------------------------

# --- Kernel: normalize root resistances by btran, over (patch, layer) ---
@kernel function _smstress_normalize_rootr_kernel!(rootr, @Const(mask_patch),
                                                   @Const(btran), btran0)
    p, j = @index(Global, NTuple)
    @inbounds if mask_patch[p]
        if btran[p] > btran0
            rootr[p, j] = rootr[p, j] / btran[p]
        else
            rootr[p, j] = 0.0
        end
    end
end

smstress_normalize_rootr!(rootr, mask_patch, btran, btran0, nlevgrnd::Int) =
    _launch!(_smstress_normalize_rootr_kernel!, rootr, mask_patch, btran, btran0;
             ndrange = (length(mask_patch), nlevgrnd))

"""
    calc_root_moist_stress_clm45default!(...)

Compute root water stress using the default CLM4.5 approach.
Ported from `calc_root_moist_stress_clm45default` in `SoilMoistStressMod.F90`.
"""
function calc_root_moist_stress_clm45default!(rootfr_unf::Matrix{<:Real},
                                               rootfr::Matrix{<:Real},
                                               rootr::Matrix{<:Real},
                                               btran::Vector{<:Real},
                                               rresis::Matrix{<:Real},
                                               smpso::Vector{<:Real},
                                               smpsc::Vector{<:Real},
                                               t_soisno::Matrix{<:Real},
                                               watsat::Matrix{<:Real},
                                               sucsat::Matrix{<:Real},
                                               bsw::Matrix{<:Real},
                                               eff_porosity::Matrix{<:Real},
                                               h2osoi_liqvol::Matrix{<:Real},
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
            itype = patch_itype[p] + 1  # 0-based Fortran PFT → 1-based Julia

            # Root resistance factors
            if h2osoi_liqvol[c, j + joff] <= 0.0 || t_soisno[c, j + joff] <= TFRZ - 2.0
                rootr[p, j] = 0.0
            else
                s_node = max(h2osoi_liqvol[c, j + joff] / eff_porosity[c, j], 0.01)
                s_node = min(s_node, 1.0)

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

    # Normalize root resistances to get layer contribution to ET.
    # btran[p] is fully accumulated above; here it is read-only, so each
    # (patch, layer) is independent.
    smstress_normalize_rootr!(rootr, mask_patch, btran, btran0, nlevgrnd)

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
                                  smpso::Vector{<:Real},
                                  smpsc::Vector{<:Real},
                                  altmax_lastyear_indx::Vector{<:Real},
                                  altmax_indx::Vector{<:Real},
                                  mask_patch::BitVector,
                                  bounds_patch::UnitRange{Int},
                                  nlevgrnd::Int,
                                  nlevsno::Int)

    np = length(bounds_patch)

    # Initialize rootfr_unf to zero
    FT = eltype(temperature.t_soisno_col)
    rootfr_unf = zeros(FT, length(mask_patch), nlevgrnd)

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
