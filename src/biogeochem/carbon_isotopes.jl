# ==========================================================================
# Carbon Isotope Tracer Transport (C13/C14)
#
# ┌────────────────────────────────────────────────────────────────────────┐
# │ STATUS: BEING WIRED LIVE (2026-07). Included in src/CLM.jl and tested    │
# │ (test_carbon_isotopes.jl in runtests.jl). C13/C14 tracers are gated by   │
# │ use_c13/use_c14 (default false → no default-path impact). The parallel   │
# │ C13/C14 carbon-state instances live on the CN vegetation facade          │
# │ (c13_/c14_cnveg_carbonstate_inst); c13_c14_photosynthesis! + c14_decay!  │
# │ thread into clm_drv! in isotope phase order, and are GPU-kernelized.     │
# └────────────────────────────────────────────────────────────────────────┘
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

"""Preindustrial atmospheric del13C (per mil). CTSM clm_varcon `preind_atm_del13c`."""
const PREIND_ATM_DEL13C = -6.0

"""Preindustrial atmospheric 13C/12C ratio. CTSM clm_varcon `preind_atm_ratio`."""
const PREIND_ATM_RATIO = C13_PDB_RATIO + (PREIND_ATM_DEL13C * C13_PDB_RATIO) / 1000.0

"""
Preindustrial atmospheric 13C/(12C+13C) mass ratio used to form the C13O2
partial pressure: `forc_pc13o2 = C13RATIO * forc_pco2`. CTSM clm_varcon
`c13ratio = preind_atm_ratio/(1+preind_atm_ratio)`. This yields
`rc13_canair = C13RATIO/(1-C13RATIO) = 0.0111718` (del13C_atm = -6 permil).
"""
const C13RATIO = PREIND_ATM_RATIO / (1.0 + PREIND_ATM_RATIO)

"""Standard atmospheric C14/C ratio (pre-bomb, ~1950). CTSM clm_varcon `c14ratio`."""
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
# Per-patch kernel: every write is own-index (ps.*_patch[p]). c14_bomb_factor is
# inlined (latitude→sector), rc14_atm passed as 3 scalars. Byte-identical on the
# KA CPU backend; device backend matches at working precision.
@kernel function _c13c14_photosynthesis_kernel!(
        rc13_canair, rc13_psnsun, rc13_psnsha, c13_psnsun, c13_psnsha,
        c14_psnsun, c14_psnsha,
        @Const(alphapsnsun), @Const(alphapsnsha), @Const(psnsun), @Const(psnsha),
        @Const(mask_patch), @Const(gridcell_of_patch),
        @Const(forc_pco2_grc), @Const(forc_pc13o2_grc), @Const(latdeg_grc),
        use_c13::Bool, use_c14::Bool, use_c13_timeseries::Bool,
        rc13_atm, rc14_1, rc14_2, rc14_3, begp::Int, endp::Int)
    p = @index(Global)
    T = eltype(rc13_canair)
    @inbounds if begp <= p <= endp && mask_patch[p]
        g = gridcell_of_patch[p]
        if use_c13
            if use_c13_timeseries
                rc13_canair[p] = rc13_atm
            else
                pco2 = forc_pco2_grc[g]; pc13o2 = forc_pc13o2_grc[g]
                rc13_canair[p] = (pco2 > pc13o2 && pc13o2 > zero(T)) ?
                                 pc13o2 / (pco2 - pc13o2) : zero(T)
            end
            alpha_sun = alphapsnsun[p]; alpha_sha = alphapsnsha[p]
            rc13_psnsun[p] = alpha_sun != zero(T) ? rc13_canair[p] / alpha_sun : zero(T)
            rc13_psnsha[p] = alpha_sha != zero(T) ? rc13_canair[p] / alpha_sha : zero(T)
            rc_sun = rc13_psnsun[p]; rc_sha = rc13_psnsha[p]
            c13_psnsun[p] = psnsun[p] * (rc_sun / (one(T) + rc_sun))
            c13_psnsha[p] = psnsha[p] * (rc_sha / (one(T) + rc_sha))
        end
        if use_c14
            lat = latdeg_grc[g]
            rc14 = lat >= T(30.0) ? rc14_1 : (lat >= T(-30.0) ? rc14_2 : rc14_3)
            c14_psnsun[p] = rc14 * psnsun[p]
            c14_psnsha[p] = rc14 * psnsha[p]
        end
    end
end

function c13_c14_photosynthesis!(ps::PhotosynthesisData,
                                 forc_pco2_grc::AbstractVector{<:Real},
                                 forc_pc13o2_grc::AbstractVector{<:Real},
                                 gridcell_of_patch::AbstractVector{<:Integer},
                                 latdeg_grc::AbstractVector{<:Real},
                                 mask_patch::AbstractVector{Bool},
                                 bounds_patch::UnitRange{Int};
                                 use_c13::Bool = false,
                                 use_c14::Bool = false,
                                 use_c13_timeseries::Bool = false,
                                 rc13_atm::Real = NaN,
                                 rc14_atm::AbstractVector{<:Real} = [1.0, 1.0, 1.0])

    isempty(bounds_patch) && return nothing
    T = eltype(ps.psnsun_patch)
    _launch!(_c13c14_photosynthesis_kernel!,
        ps.rc13_canair_patch, ps.rc13_psnsun_patch, ps.rc13_psnsha_patch,
        ps.c13_psnsun_patch, ps.c13_psnsha_patch, ps.c14_psnsun_patch, ps.c14_psnsha_patch,
        ps.alphapsnsun_patch, ps.alphapsnsha_patch, ps.psnsun_patch, ps.psnsha_patch,
        mask_patch, gridcell_of_patch, forc_pco2_grc, forc_pc13o2_grc, latdeg_grc,
        use_c13, use_c14, use_c13_timeseries,
        T(rc13_atm), T(rc14_atm[1]), T(rc14_atm[2]), T(rc14_atm[3]),
        first(bounds_patch), last(bounds_patch); ndrange = last(bounds_patch))
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
# --- Kernels for the (default, non-matrix) C14 decay path ------------------
# All writes are own-index; the soil kernel runs per-column with internal l/j
# loops (fac computed once per (c,l), applied over all j — byte-identical).
@kernel function _c14decay_seedc_kernel!(seedc_grc, decay_factor, begg::Int, endg::Int)
    g = @index(Global)
    @inbounds if begg <= g <= endg
        seedc_grc[g] = seedc_grc[g] * decay_factor
    end
end

@kernel function _c14decay_soil_kernel!(decomp_cpools_vr, @Const(mask_soilc),
        @Const(spinup_factor), @Const(spinup_lat_term),
        decay_const, dt, spinup_state::Int, nlevdecomp::Int, ndecomp_pools::Int,
        begc::Int, endc::Int)
    c = @index(Global)
    T = eltype(decomp_cpools_vr)
    @inbounds if begc <= c <= endc && mask_soilc[c]
        for l in 1:ndecomp_pools
            sp_term = one(T)
            if spinup_state >= 1
                sp_term = spinup_factor[l]
                if abs(spinup_factor[l] - one(T)) > T(1.0e-6)
                    sp_term = sp_term * spinup_lat_term[c]
                end
            end
            fac = one(T) - decay_const * sp_term * dt
            for j in 1:nlevdecomp
                decomp_cpools_vr[c, j, l] = decomp_cpools_vr[c, j, l] * fac
            end
        end
    end
end

# Split into two kernels over the same patches to stay under Metal's ~31
# top-level buffer limit (24 pool arrays + mask + ivt exceeds it as one kernel).
# All writes are own-index and independent, so two launches are equivalent.
@kernel function _c14decay_patch1_kernel!(
        cpool, xsmrpool,
        deadcrootc, deadcrootc_storage, deadcrootc_xfer,
        deadstemc, deadstemc_storage, deadstemc_xfer,
        frootc, frootc_storage, frootc_xfer, leafc,
        @Const(mask_soilp), decay_factor, begp::Int, endp::Int)
    p = @index(Global)
    @inbounds if begp <= p <= endp && mask_soilp[p]
        cpool[p]    = cpool[p]    * decay_factor
        xsmrpool[p] = xsmrpool[p] * decay_factor
        deadcrootc[p]         = deadcrootc[p]         * decay_factor
        deadcrootc_storage[p] = deadcrootc_storage[p] * decay_factor
        deadcrootc_xfer[p]    = deadcrootc_xfer[p]    * decay_factor
        deadstemc[p]          = deadstemc[p]          * decay_factor
        deadstemc_storage[p]  = deadstemc_storage[p]  * decay_factor
        deadstemc_xfer[p]     = deadstemc_xfer[p]     * decay_factor
        frootc[p]             = frootc[p]             * decay_factor
        frootc_storage[p]     = frootc_storage[p]     * decay_factor
        frootc_xfer[p]        = frootc_xfer[p]        * decay_factor
        leafc[p]              = leafc[p]              * decay_factor
    end
end

@kernel function _c14decay_patch2_kernel!(
        leafc_storage, leafc_xfer,
        livecrootc, livecrootc_storage, livecrootc_xfer,
        livestemc, livestemc_storage, livestemc_xfer,
        gresp_storage, gresp_xfer, ctrunc, cropseedc_deficit,
        @Const(mask_soilp), @Const(ivt),
        decay_factor, npcropmin::Int, ivt_len::Int, begp::Int, endp::Int)
    p = @index(Global)
    @inbounds if begp <= p <= endp && mask_soilp[p]
        leafc_storage[p]      = leafc_storage[p]      * decay_factor
        leafc_xfer[p]         = leafc_xfer[p]         * decay_factor
        livecrootc[p]         = livecrootc[p]         * decay_factor
        livecrootc_storage[p] = livecrootc_storage[p] * decay_factor
        livecrootc_xfer[p]    = livecrootc_xfer[p]    * decay_factor
        livestemc[p]          = livestemc[p]          * decay_factor
        livestemc_storage[p]  = livestemc_storage[p]  * decay_factor
        livestemc_xfer[p]     = livestemc_xfer[p]     * decay_factor
        gresp_storage[p]      = gresp_storage[p]      * decay_factor
        gresp_xfer[p]         = gresp_xfer[p]         * decay_factor
        ctrunc[p]             = ctrunc[p]             * decay_factor
        if p <= ivt_len && ivt[p] >= npcropmin
            cropseedc_deficit[p] = cropseedc_deficit[p] * decay_factor
        end
    end
end

function c14_decay!(c14_cnveg_cs::CNVegCarbonStateData,
                    c14_cnveg_cf::CNVegCarbonFluxData,
                    c14_soilbgc_cs::SoilBiogeochemCarbonStateData,
                    c14_soilbgc_cf::SoilBiogeochemCarbonFluxData;
                    mask_soilc::AbstractVector{Bool},
                    mask_soilp::AbstractVector{Bool},
                    bounds_col::UnitRange{Int},
                    bounds_patch::UnitRange{Int},
                    bounds_gridcell::UnitRange{Int},
                    dt::Real,
                    days_per_year::Real = 365.0,
                    nlevdecomp::Int,
                    ndecomp_pools::Int,
                    spinup_state::Int = 0,
                    spinup_factor::AbstractVector{<:Real} = ones(ndecomp_pools),
                    col_gridcell::AbstractVector{<:Integer} = Int[],
                    latdeg_grc::AbstractVector{<:Real} = Float64[],
                    ivt::AbstractVector{<:Integer} = Int[],
                    npcropmin::Int = 17,
                    use_matrixcn::Bool = false,
                    use_soil_matrixcn::Bool = false)

    # Compute decay constant
    half_life = C14_HALF_LIFE_YEARS * SECSPDAY * days_per_year
    decay_const = -log(0.5) / half_life
    decay_factor = 1.0 - decay_const * dt

    cs  = c14_cnveg_cs
    scs = c14_soilbgc_cs
    FT  = eltype(cs.cpool_patch)
    # backend-resident copies of host param arrays (identity on CPU)
    _onbk(ref, a, ET) = ref isa Array ? a : copyto!(similar(ref, ET, length(a)), a)

    # --- Gridcell-level: seedc ---
    if !isempty(bounds_gridcell)
        _launch!(_c14decay_seedc_kernel!, cs.seedc_grc, FT(decay_factor),
                 first(bounds_gridcell), last(bounds_gridcell);
                 ndrange = last(bounds_gridcell))
    end

    # --- Column-level: soil decomposition pools ---
    if use_soil_matrixcn
        # matrix path (host; exercised when matrix-CN is wired)
        for l in 1:ndecomp_pools, j in 1:nlevdecomp, c in bounds_col
            mask_soilc[c] || continue
            if spinup_state >= 1
                sp_term = spinup_factor[l]
                if abs(spinup_factor[l] - 1.0) > 1.0e-6
                    sp_term *= get_spinup_latitude_term(latdeg_grc[col_gridcell[c]])
                end
            else
                sp_term = 1.0
            end
            idx = j + nlevdecomp * (l - 1)
            c14_soilbgc_cf.matrix_decomp_fire_k_col[c, idx] -= sp_term * decay_const * dt
        end
    elseif !isempty(bounds_col)
        # spinup latitude term per column (device-safe; only read when spinup_state>=1)
        spinup_lat = fill!(similar(scs.decomp_cpools_vr_col, FT, last(bounds_col)), zero(FT))
        if spinup_state >= 1 && !isempty(col_gridcell) && !isempty(latdeg_grc)
            fire_spinup_latterm!(spinup_lat, mask_soilc,
                _onbk(scs.decomp_cpools_vr_col, col_gridcell, Int), latdeg_grc)
        end
        _launch!(_c14decay_soil_kernel!, scs.decomp_cpools_vr_col, mask_soilc,
                 _onbk(scs.decomp_cpools_vr_col, spinup_factor, FT), spinup_lat,
                 FT(decay_const), FT(dt), spinup_state, nlevdecomp, ndecomp_pools,
                 first(bounds_col), last(bounds_col); ndrange = last(bounds_col))
    end

    # --- Patch-level: vegetation C pools ---
    if use_matrixcn
        # matrix path (host; encode decay as a fire transfer rate)
        cf = c14_cnveg_cf
        for p in bounds_patch
            mask_soilp[p] || continue
            cs.cpool_patch[p]    *= decay_factor
            cs.xsmrpool_patch[p] *= decay_factor
            cf.matrix_fitransfer_patch[p, cf.ideadcroot_to_iout_fi[1]]   = decay_const
            cf.matrix_fitransfer_patch[p, cf.ideadcrootst_to_iout_fi[1]]  = decay_const
            cf.matrix_fitransfer_patch[p, cf.ideadcrootxf_to_iout_fi[1]]  = decay_const
            cf.matrix_fitransfer_patch[p, cf.ideadstem_to_iout_fi[1]]     = decay_const
            cf.matrix_fitransfer_patch[p, cf.ideadstemst_to_iout_fi[1]]   = decay_const
            cf.matrix_fitransfer_patch[p, cf.ideadstemxf_to_iout_fi[1]]   = decay_const
            cf.matrix_fitransfer_patch[p, cf.ifroot_to_iout_fi[1]]        = decay_const
            cf.matrix_fitransfer_patch[p, cf.ifrootst_to_iout_fi[1]]      = decay_const
            cf.matrix_fitransfer_patch[p, cf.ifrootxf_to_iout_fi[1]]      = decay_const
            cf.matrix_fitransfer_patch[p, cf.ileaf_to_iout_fi[1]]         = decay_const
            cf.matrix_fitransfer_patch[p, cf.ileafst_to_iout_fi[1]]       = decay_const
            cf.matrix_fitransfer_patch[p, cf.ileafxf_to_iout_fi[1]]       = decay_const
            cf.matrix_fitransfer_patch[p, cf.ilivecroot_to_iout_fi[1]]    = decay_const
            cf.matrix_fitransfer_patch[p, cf.ilivecrootst_to_iout_fi[1]]  = decay_const
            cf.matrix_fitransfer_patch[p, cf.ilivecrootxf_to_iout_fi[1]]  = decay_const
            cf.matrix_fitransfer_patch[p, cf.ilivestem_to_iout_fi[1]]     = decay_const
            cf.matrix_fitransfer_patch[p, cf.ilivestemst_to_iout_fi[1]]   = decay_const
            cf.matrix_fitransfer_patch[p, cf.ilivestemxf_to_iout_fi[1]]   = decay_const
            cs.gresp_storage_patch[p] *= decay_factor
            cs.gresp_xfer_patch[p]    *= decay_factor
            cs.ctrunc_patch[p]        *= decay_factor
            if length(ivt) >= p && ivt[p] >= npcropmin
                cs.cropseedc_deficit_patch[p] *= decay_factor
            end
        end
    elseif !isempty(bounds_patch)
        _launch!(_c14decay_patch1_kernel!,
            cs.cpool_patch, cs.xsmrpool_patch,
            cs.deadcrootc_patch, cs.deadcrootc_storage_patch, cs.deadcrootc_xfer_patch,
            cs.deadstemc_patch, cs.deadstemc_storage_patch, cs.deadstemc_xfer_patch,
            cs.frootc_patch, cs.frootc_storage_patch, cs.frootc_xfer_patch, cs.leafc_patch,
            mask_soilp, FT(decay_factor),
            first(bounds_patch), last(bounds_patch); ndrange = last(bounds_patch))
        _launch!(_c14decay_patch2_kernel!,
            cs.leafc_storage_patch, cs.leafc_xfer_patch,
            cs.livecrootc_patch, cs.livecrootc_storage_patch, cs.livecrootc_xfer_patch,
            cs.livestemc_patch, cs.livestemc_storage_patch, cs.livestemc_xfer_patch,
            cs.gresp_storage_patch, cs.gresp_xfer_patch, cs.ctrunc_patch,
            cs.cropseedc_deficit_patch,
            mask_soilp, _onbk(cs.cpool_patch, ivt, Int),
            FT(decay_factor), npcropmin, length(ivt),
            first(bounds_patch), last(bounds_patch); ndrange = last(bounds_patch))
    end

    return nothing
end
