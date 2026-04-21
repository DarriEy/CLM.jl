# ==========================================================================
# Carbon Isotope Tracer Transport (C13/C14)
#
# Ported from:
#   - src/biogeochem/CNC14DecayMod.F90        — C14 radioactive decay
#   - src/biogeophys/PhotosynthesisMod.F90     — C13/C14 fractionation during photosynthesis
#   - src/biogeochem/SpeciesIsotopeType.F90    — isotope ratio type
#   - src/biogeochem/SpeciesNonIsotopeType.F90 — C12 version
#
# Public functions:
#   c14_decay!                   — Apply radioactive decay to all C14 pools
#   c13_c14_photosynthesis!      — Compute C13/C14 ratios and photosynthesis fluxes
#   delta13C_to_ratio            — Convert delta-13C (per mil) to C13/(C12+C13) ratio
#   ratio_to_delta13C            — Convert C13/(C12+C13) ratio to delta-13C (per mil)
#   c14_bomb_factor              — Get C14/C ratio from atmospheric bomb spike by latitude
# ==========================================================================

# ---------------------------------------------------------------------------
# Physical constants for carbon isotopes
# ---------------------------------------------------------------------------

"""C14 half-life in years."""
const C14_HALF_LIFE_YEARS = 5730.0

"""PDB standard ratio for C13/C12 (Vienna PDB)."""
const C13_PDB_RATIO = 0.0112372

"""Standard atmospheric C14/C ratio (pre-bomb, ~1950)."""
const C14_ATM_RATIO_PREBOMB = 1.0e-12

# ---------------------------------------------------------------------------
# Utility: delta-13C <-> ratio conversions
# ---------------------------------------------------------------------------

"""
    delta13C_to_ratio(delta13C) -> Float64

Convert delta-13C (per mil, VPDB) to the mass ratio C13/(C12+C13).

The delta notation is: delta = (R_sample/R_standard - 1) * 1000
where R = C13/C12. We return C13/(C12+C13) = R/(1+R).
"""
function delta13C_to_ratio(delta13C::Real)
    R = C13_PDB_RATIO * (1.0 + delta13C / 1000.0)
    return R / (1.0 + R)
end

"""
    ratio_to_delta13C(ratio) -> Float64

Convert C13/(C12+C13) mass ratio to delta-13C (per mil, VPDB).
"""
function ratio_to_delta13C(ratio::Real)
    # ratio = R/(1+R), so R = ratio/(1-ratio)
    if ratio <= 0.0 || ratio >= 1.0
        return 0.0
    end
    R = ratio / (1.0 - ratio)
    return (R / C13_PDB_RATIO - 1.0) * 1000.0
end

# ---------------------------------------------------------------------------
# c14_bomb_factor — atmospheric C14/C ratio by latitude sector
# Ported from PhotosynthesisTotal in PhotosynthesisMod.F90
# ---------------------------------------------------------------------------

"""
    c14_bomb_factor(latdeg; use_c14_bombspike=false, atm_c14_filename="") -> Float64

Return the atmospheric C14/C ratio for the given latitude.
In the simple (non-timeseries) case, returns the pre-bomb ratio of 1.0
(i.e., C14 tracks C12 with ratio=1). In CLM, the bomb spike is read from
a time series file; here we provide the pre-bomb default.

Latitude sectors (Fortran convention):
- sector 1: lat >= 30N
- sector 2: -30 <= lat < 30
- sector 3: lat < -30S
"""
function c14_bomb_factor(latdeg::Real;
                         rc14_atm::Vector{<:Real} = [1.0, 1.0, 1.0])
    if latdeg >= 30.0
        return rc14_atm[1]
    elseif latdeg >= -30.0
        return rc14_atm[2]
    else
        return rc14_atm[3]
    end
end

# ---------------------------------------------------------------------------
# c13_c14_photosynthesis! — Compute C13/C14 photosynthesis fluxes
# Ported from PhotosynthesisTotal in PhotosynthesisMod.F90
# ---------------------------------------------------------------------------

"""
    c13_c14_photosynthesis!(ps, forc_pco2_grc, forc_pc13o2_grc,
                            gridcell_of_patch, latdeg_grc,
                            mask_patch, bounds_patch;
                            use_c13, use_c14, use_c13_timeseries,
                            rc13_atm, rc14_atm)

Compute C13 and C14 photosynthesis fluxes from total (C12) photosynthesis.

For C13:
  - rc13_canair = forc_pc13o2 / (forc_pco2 - forc_pc13o2)   [or from timeseries]
  - rc13_psnsun = rc13_canair / alphapsnsun                  [fractionation]
  - c13_psnsun  = psnsun * rc13_psnsun / (1 + rc13_psnsun)

For C14:
  - c14_psnsun = rc14_atm(sector) * psnsun
  - (C14 fractionation is doubled relative to C13 in the flux calc, not here)

Ported from lines 2088-2134 of `PhotosynthesisTotal` in `PhotosynthesisMod.F90`.
"""
function c13_c14_photosynthesis!(ps::PhotosynthesisData,
                                 forc_pco2_grc::Vector{<:Real},
                                 forc_pc13o2_grc::Vector{<:Real},
                                 gridcell_of_patch::Vector{Int},
                                 latdeg_grc::Vector{<:Real},
                                 mask_patch::BitVector,
                                 bounds_patch::UnitRange{Int};
                                 use_c13::Bool = false,
                                 use_c14::Bool = false,
                                 use_c13_timeseries::Bool = false,
                                 rc13_atm::Real = NaN,
                                 rc14_atm::Vector{<:Real} = [1.0, 1.0, 1.0])

    for p in bounds_patch
        mask_patch[p] || continue
        g = gridcell_of_patch[p]

        # --- C13 fractionation ---
        if use_c13
            if use_c13_timeseries
                ps.rc13_canair_patch[p] = rc13_atm
            else
                pco2 = forc_pco2_grc[g]
                pc13o2 = forc_pc13o2_grc[g]
                if pco2 > pc13o2 && pc13o2 > 0.0
                    ps.rc13_canair_patch[p] = pc13o2 / (pco2 - pc13o2)
                else
                    ps.rc13_canair_patch[p] = 0.0
                end
            end

            # Apply fractionation from Fractionation() subroutine
            alpha_sun = ps.alphapsnsun_patch[p]
            alpha_sha = ps.alphapsnsha_patch[p]

            if alpha_sun != 0.0
                ps.rc13_psnsun_patch[p] = ps.rc13_canair_patch[p] / alpha_sun
            else
                ps.rc13_psnsun_patch[p] = 0.0
            end

            if alpha_sha != 0.0
                ps.rc13_psnsha_patch[p] = ps.rc13_canair_patch[p] / alpha_sha
            else
                ps.rc13_psnsha_patch[p] = 0.0
            end

            # C13 photosynthesis = total_psn * rc13 / (1 + rc13)
            rc_sun = ps.rc13_psnsun_patch[p]
            rc_sha = ps.rc13_psnsha_patch[p]
            ps.c13_psnsun_patch[p] = ps.psnsun_patch[p] * (rc_sun / (1.0 + rc_sun))
            ps.c13_psnsha_patch[p] = ps.psnsha_patch[p] * (rc_sha / (1.0 + rc_sha))
        end

        # --- C14 photosynthesis ---
        if use_c14
            lat = latdeg_grc[g]
            rc14 = c14_bomb_factor(lat; rc14_atm=rc14_atm)
            ps.c14_psnsun_patch[p] = rc14 * ps.psnsun_patch[p]
            ps.c14_psnsha_patch[p] = rc14 * ps.psnsha_patch[p]
        end
    end

    return nothing
end

# ---------------------------------------------------------------------------
# c14_decay! — C14 radioactive decay
# Ported from: subroutine C14Decay in CNC14DecayMod.F90
# ---------------------------------------------------------------------------

"""
    c14_decay!(c14_cnveg_cs, c14_cnveg_cf, c14_soilbgc_cs, c14_soilbgc_cf;
               mask_soilc, mask_soilp, bounds_col, bounds_patch, bounds_gridcell,
               dt, days_per_year, nlevdecomp, ndecomp_pools,
               spinup_state, spinup_factor, col_gridcell, latdeg_grc,
               ivt, npcropmin, use_matrixcn, use_soil_matrixcn)

Apply radioactive decay to all C14 carbon pools.

The C14 half-life is 5730 years. The decay constant is:
    decay_const = ln(2) / (5730 * SECSPDAY * days_per_year)

Each pool is multiplied by (1 - decay_const * dt) per timestep.

For soil decomposition pools, the decay can be accelerated during
spinup by the same factor used for decomposition (spinup_factor).

Ported from `C14Decay` in `CNC14DecayMod.F90`.
"""
function c14_decay!(c14_cnveg_cs::CNVegCarbonStateData,
                    c14_cnveg_cf::CNVegCarbonFluxData,
                    c14_soilbgc_cs::SoilBiogeochemCarbonStateData,
                    c14_soilbgc_cf::SoilBiogeochemCarbonFluxData;
                    mask_soilc::BitVector,
                    mask_soilp::BitVector,
                    bounds_col::UnitRange{Int},
                    bounds_patch::UnitRange{Int},
                    bounds_gridcell::UnitRange{Int},
                    dt::Real,
                    days_per_year::Real = 365.0,
                    nlevdecomp::Int,
                    ndecomp_pools::Int,
                    spinup_state::Int = 0,
                    spinup_factor::Vector{<:Real} = ones(ndecomp_pools),
                    col_gridcell::Vector{Int} = Int[],
                    latdeg_grc::Vector{<:Real} = Float64[],
                    ivt::Vector{Int} = Int[],
                    npcropmin::Int = 17,
                    use_matrixcn::Bool = false,
                    use_soil_matrixcn::Bool = false)

    # Compute decay constant
    half_life = C14_HALF_LIFE_YEARS * SECSPDAY * days_per_year
    decay_const = -log(0.5) / half_life
    decay_factor = 1.0 - decay_const * dt

    # --- Gridcell-level: seedc ---
    for g in bounds_gridcell
        c14_cnveg_cs.seedc_grc[g] *= decay_factor
    end

    # --- Column-level: soil decomposition pools ---
    for l in 1:ndecomp_pools
        for j in 1:nlevdecomp
            for c in bounds_col
                mask_soilc[c] || continue

                if spinup_state >= 1
                    # Accelerate radioactive decay during spinup
                    sp_term = spinup_factor[l]
                    if abs(spinup_factor[l] - 1.0) > 1.0e-6
                        sp_term *= get_spinup_latitude_term(latdeg_grc[col_gridcell[c]])
                    end
                else
                    sp_term = 1.0
                end

                if !use_soil_matrixcn
                    c14_soilbgc_cs.decomp_cpools_vr_col[c, j, l] *= (1.0 - decay_const * sp_term * dt)
                else
                    idx = j + nlevdecomp * (l - 1)
                    c14_soilbgc_cf.matrix_decomp_fire_k_col[c, idx] -= sp_term * decay_const * dt
                end
            end
        end
    end

    # --- Patch-level: vegetation C pools ---
    for p in bounds_patch
        mask_soilp[p] || continue

        c14_cnveg_cs.cpool_patch[p]    *= decay_factor
        c14_cnveg_cs.xsmrpool_patch[p] *= decay_factor

        if !use_matrixcn
            c14_cnveg_cs.deadcrootc_patch[p]         *= decay_factor
            c14_cnveg_cs.deadcrootc_storage_patch[p]  *= decay_factor
            c14_cnveg_cs.deadcrootc_xfer_patch[p]     *= decay_factor
            c14_cnveg_cs.deadstemc_patch[p]           *= decay_factor
            c14_cnveg_cs.deadstemc_storage_patch[p]   *= decay_factor
            c14_cnveg_cs.deadstemc_xfer_patch[p]      *= decay_factor
            c14_cnveg_cs.frootc_patch[p]              *= decay_factor
            c14_cnveg_cs.frootc_storage_patch[p]      *= decay_factor
            c14_cnveg_cs.frootc_xfer_patch[p]         *= decay_factor
            c14_cnveg_cs.leafc_patch[p]               *= decay_factor
            c14_cnveg_cs.leafc_storage_patch[p]       *= decay_factor
            c14_cnveg_cs.leafc_xfer_patch[p]          *= decay_factor
            c14_cnveg_cs.livecrootc_patch[p]          *= decay_factor
            c14_cnveg_cs.livecrootc_storage_patch[p]  *= decay_factor
            c14_cnveg_cs.livecrootc_xfer_patch[p]     *= decay_factor
            c14_cnveg_cs.livestemc_patch[p]           *= decay_factor
            c14_cnveg_cs.livestemc_storage_patch[p]   *= decay_factor
            c14_cnveg_cs.livestemc_xfer_patch[p]      *= decay_factor
        else
            # Matrix mode: encode decay as additional fire transfer rate
            c14_cnveg_cf.matrix_fitransfer_patch[p, c14_cnveg_cf.ideadcroot_to_iout_fi[1]]   = decay_const
            c14_cnveg_cf.matrix_fitransfer_patch[p, c14_cnveg_cf.ideadcrootst_to_iout_fi[1]]  = decay_const
            c14_cnveg_cf.matrix_fitransfer_patch[p, c14_cnveg_cf.ideadcrootxf_to_iout_fi[1]]  = decay_const
            c14_cnveg_cf.matrix_fitransfer_patch[p, c14_cnveg_cf.ideadstem_to_iout_fi[1]]     = decay_const
            c14_cnveg_cf.matrix_fitransfer_patch[p, c14_cnveg_cf.ideadstemst_to_iout_fi[1]]   = decay_const
            c14_cnveg_cf.matrix_fitransfer_patch[p, c14_cnveg_cf.ideadstemxf_to_iout_fi[1]]   = decay_const
            c14_cnveg_cf.matrix_fitransfer_patch[p, c14_cnveg_cf.ifroot_to_iout_fi[1]]        = decay_const
            c14_cnveg_cf.matrix_fitransfer_patch[p, c14_cnveg_cf.ifrootst_to_iout_fi[1]]      = decay_const
            c14_cnveg_cf.matrix_fitransfer_patch[p, c14_cnveg_cf.ifrootxf_to_iout_fi[1]]      = decay_const
            c14_cnveg_cf.matrix_fitransfer_patch[p, c14_cnveg_cf.ileaf_to_iout_fi[1]]         = decay_const
            c14_cnveg_cf.matrix_fitransfer_patch[p, c14_cnveg_cf.ileafst_to_iout_fi[1]]       = decay_const
            c14_cnveg_cf.matrix_fitransfer_patch[p, c14_cnveg_cf.ileafxf_to_iout_fi[1]]       = decay_const
            c14_cnveg_cf.matrix_fitransfer_patch[p, c14_cnveg_cf.ilivecroot_to_iout_fi[1]]    = decay_const
            c14_cnveg_cf.matrix_fitransfer_patch[p, c14_cnveg_cf.ilivecrootst_to_iout_fi[1]]  = decay_const
            c14_cnveg_cf.matrix_fitransfer_patch[p, c14_cnveg_cf.ilivecrootxf_to_iout_fi[1]]  = decay_const
            c14_cnveg_cf.matrix_fitransfer_patch[p, c14_cnveg_cf.ilivestem_to_iout_fi[1]]     = decay_const
            c14_cnveg_cf.matrix_fitransfer_patch[p, c14_cnveg_cf.ilivestemst_to_iout_fi[1]]   = decay_const
            c14_cnveg_cf.matrix_fitransfer_patch[p, c14_cnveg_cf.ilivestemxf_to_iout_fi[1]]   = decay_const
        end

        c14_cnveg_cs.gresp_storage_patch[p] *= decay_factor
        c14_cnveg_cs.gresp_xfer_patch[p]    *= decay_factor
        c14_cnveg_cs.ctrunc_patch[p]        *= decay_factor

        # Crop seed deficit decay
        if length(ivt) >= p && ivt[p] >= npcropmin
            c14_cnveg_cs.cropseedc_deficit_patch[p] *= decay_factor
        end
    end

    return nothing
end
