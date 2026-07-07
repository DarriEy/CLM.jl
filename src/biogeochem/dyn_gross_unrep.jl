# ==========================================================================
# Ported from: src/dyn_subgrid/dynGrossUnrepMod.F90
#
# Handle reading of the gross unrepresented landcover change data
# (UNREPRESENTED_PFT_LULCC time series), as well as the disturbance-like C/N
# state updates that happen as a result. Structurally a close sibling of the
# gap-mortality routine (CNGapMortalityMod): a CN mortality-flux kernel plus a
# patch -> column aggregation.
#
# Public functions (preserve Fortran names):
#   dynGrossUnrep_init           — open the file, set up the dyn variable
#   dynGrossUnrep_interp!        — read the current-year gross-unrep fraction
#   cn_gross_unrep!              (= CNGrossUnrep)            — patch-level C/N fluxes
#   cn_gross_unrep_pft_to_column! (= CNGrossUnrepPftToColumn) — patch -> column
#
# Module-private Fortran singletons (dyngrossunrep_file, gru_inst,
# grossunrepfrac, do_grossunrep) are bundled into an explicit
# `DynGrossUnrepState` struct that the caller owns and passes in (NOT added to
# CLMInstances or any dual-copied struct).
#
# Conventions (CLAUDE.md): preserve Fortran names; mutable @kwdef structs;
# Float64; SoA arrays; BitVector masks; mask-based loops over a filtered range.
# ==========================================================================

# ---------------------------------------------------------------------------
# Module-private state singleton
# ---------------------------------------------------------------------------

"""
    DynGrossUnrepState

Bundles the four module-private Fortran globals of `dynGrossUnrepMod`:

- `dyngrossunrep_file` : the `DynFile` for the gross-unrep file.
- `gru_inst`           : the `DynVarTimeUninterp` reader for `UNREPRESENTED_PFT_LULCC`.
- `grossunrepfrac`     : gross unrepresented landcover change fraction
                         (gridcell, 0:natpft_size-1 → stored 1-based: gridcell, natpft_size).
- `do_grossunrep`      : whether we are in a period when gross unrep should be applied.

The caller owns one of these and passes it into `cn_gross_unrep!`; it is NOT
added to `CLMInstances` or any dual-copied struct.

Mirrors the module-private globals at the top of `dynGrossUnrepMod.F90`.
"""
Base.@kwdef mutable struct DynGrossUnrepState
    dyngrossunrep_file::Union{DynFile, Nothing} = nothing
    gru_inst::Union{DynVarTimeUninterp, Nothing} = nothing
    grossunrepfrac::AbstractMatrix{<:Real} = Matrix{Float64}(undef, 0, 0)  # (gridcell, natpft_size; Abstract so a device array can be stored)
    do_grossunrep::Bool = false
end

# Underlying gross-unrep variable name on file (Fortran: grossunrep_varname).
const GROSSUNREP_VARNAME = "UNREPRESENTED_PFT_LULCC"

# ---------------------------------------------------------------------------
# PFT constants needed by gross unrep (subset of pftcon)
# ---------------------------------------------------------------------------

"""
    PftConGrossUnrep

PFT-level parameters referenced by the gross-unrep routines. The subset of
`pftconMod` fields used by `dynGrossUnrepMod.F90`.

Fields:
- `pconv` : proportion of deadstem routed to the conversion (atmosphere) flux
            (Fortran `pftcon%pconv`, aliased there as `convfrac`).
- `lf_f`  : leaf litter fractions (pft, litr-pool).
- `fr_f`  : fine root litter fractions (pft, litr-pool).
"""
Base.@kwdef mutable struct PftConGrossUnrep{FT<:Real, V<:AbstractVector{FT}, M<:AbstractMatrix{FT}}
    pconv ::V = Float64[]                          # deadstem -> conversion flux fraction
    lf_f  ::M = Matrix{Float64}(undef, 0, 0)       # leaf litter fractions (pft, litr)
    fr_f  ::M = Matrix{Float64}(undef, 0, 0)       # fine root litter fractions (pft, litr)
end
PftConGrossUnrep{FT}(; kwargs...) where {FT<:Real} = PftConGrossUnrep{FT, Vector{FT}, Matrix{FT}}(; kwargs...)

# ---------------------------------------------------------------------------
# dynGrossUnrep_init — open the file + set up the dyn variable
# ---------------------------------------------------------------------------

"""
    dynGrossUnrep_init(begg, endg, natpft_size, grossunrep_filename;
                       current_year=nothing) -> DynGrossUnrepState

Initialize data structures for gross-unrep information. Opens the gross-unrep
file (via `DynFile`) and sets up the `UNREPRESENTED_PFT_LULCC`
`DynVarTimeUninterp` variable. Should be called once, during model
initialization.

`begg:endg` are proc-level gridcell bounds; `natpft_size` is the number of
natural PFTs. The Fortran allocates `grossunrepfrac(begg:endg, 0:natpft_size-1)`
— here it is a `(num_points, natpft_size)` matrix with the PFT axis stored
1-based (PFT 0 → column 1).

Year handling uses `YEAR_POSITION_START_OF_TIMESTEP`, matching the Fortran.

Ported from `dynGrossUnrep_init` in `dynGrossUnrepMod.F90`.
"""
function dynGrossUnrep_init(begg::Int, endg::Int, natpft_size::Int,
                            grossunrep_filename::String;
                            current_year::Union{Int,Nothing} = nothing)
    num_points = endg - begg + 1

    # Get the year from the START of the timestep for consistency with the other
    # dyn-file machinery (Fortran: YEAR_POSITION_START_OF_TIMESTEP).
    dyngrossunrep_file = dyn_file_open(grossunrep_filename,
                                       YEAR_POSITION_START_OF_TIMESTEP;
                                       current_year = current_year)

    # Get initial grossunrep data. The data are stored on file with shape
    # [num_points, natpft_size]; allow_nodata so a missing field zeroes out.
    gru_inst = dyn_var_time_uninterp(dyngrossunrep_file,
                                     GROSSUNREP_VARNAME,
                                     "grlnd",                # dim1name (spatial)
                                     1.0,                    # conversion_factor
                                     false,                  # do_check_sums_equal_1
                                     [num_points, natpft_size];
                                     allow_nodata = true)

    return DynGrossUnrepState(
        dyngrossunrep_file = dyngrossunrep_file,
        gru_inst           = gru_inst,
        grossunrepfrac     = zeros(Float64, num_points, natpft_size),
        do_grossunrep      = false,
    )
end

# ---------------------------------------------------------------------------
# dynGrossUnrep_interp! — read the current-year gross-unrep fraction
# ---------------------------------------------------------------------------

"""
    dynGrossUnrep_interp!(state::DynGrossUnrepState, year::Int)

Get gross-unrep data for the model time, when needed.

Gross-unrep data are stored as rates (not weights) and so time interpolation is
not necessary — the gross-unrep rate is held constant through the year.

`year` is the model year for the start of the timestep (the caller supplies the
year `get_prev_date` would have produced, per the file's year position).

After the call, `state.grossunrepfrac` holds the current-year fraction
(gridcell, PFT) and `state.do_grossunrep` says whether gross unrep is active.
Before the start of the time series gross unrep is OFF; past the end the last
year's rates are maintained (do_grossunrep stays true).

Ported from `dynGrossUnrep_interp` in `dynGrossUnrepMod.F90`.
"""
function dynGrossUnrep_interp!(state::DynGrossUnrepState, year::Int)
    set_current_year!(state.dyngrossunrep_file, year)

    # Zero out grossunrepfrac for this time step (Fortran: set to 0 then fill).
    fill!(state.grossunrepfrac, 0.0)

    if is_before_time_series(state.dyngrossunrep_file.time_info)
        # Turn off grossunrep before the start of the time series.
        state.do_grossunrep = false
    else
        # do_grossunrep stays true even past the end of the time series: the last
        # year's rates are maintained for all subsequent years.
        state.do_grossunrep = true
        data = get_current_data_2d(state.gru_inst)  # (num_points, natpft_size)
        state.grossunrepfrac .= data
    end

    return nothing
end

# --- Gross-unrep kernels. Like harvest, every output = input * (rate), but the
#     rate is m for most pools, m*cfrac for deadstem→atm, m*(1-cfrac) for wood
#     products. Compute the three rate arrays + the natural-pft mask once, then
#     reuse the generic masked-scale kernel (_dynh_scale_kernel!) per field. ---
@kernel function _gru_m_kernel!(gmask, m, mc, mpr, @Const(mask), @Const(ivt),
        @Const(gridcell), @Const(grossunrepfrac), @Const(convfrac),
        is_beg::Bool, do_gru::Bool, dtime, noveg::Int, nc4::Int)
    p = @index(Global)
    @inbounds begin
        T = eltype(m)
        isnat = mask[p] && ivt[p] > noveg && ivt[p] <= nc4
        gmask[p] = isnat
        mm = zero(T); cf = zero(T)
        if isnat
            cf = convfrac[ivt[p] + 1]
            if do_gru
                am = grossunrepfrac[gridcell[p], ivt[p] + 1]
                mm = is_beg ? am / dtime : zero(T)
            end
        end
        m[p] = mm; mc[p] = mm * cf; mpr[p] = mm * (one(T) - cf)
    end
end

# gross-unrep p2c scatter: only 9 patch fields feed it (< 31 buffers), so no
# packing — the per-patch kernel reads them directly and atomic-scatters.
@kernel function _gru_scatter_kernel!(@Const(mask), @Const(column), @Const(ivt),
        @Const(leaf_prof), @Const(froot_prof), @Const(croot_prof), @Const(lf_f), @Const(fr_f),
        @Const(wtcol), @Const(gru_leafc), @Const(gru_frootc), @Const(gru_livecrootc),
        @Const(gru_deadcrootc), @Const(gru_leafn), @Const(gru_frootn), @Const(gru_livecrootn),
        @Const(gru_deadcrootn), @Const(gru_retransn), litr_c, cwdc, litr_n, cwdn,
        nlevdecomp::Int, i_litr_min::Int, i_litr_max::Int, i_met_lit::Int)
    p = @index(Global)
    @inbounds if mask[p]
        c = column[p]; ivtp = ivt[p] + 1; wt = wtcol[p]
        for j in 1:nlevdecomp
            lp = leaf_prof[p,j]; fp = froot_prof[p,j]; cp = croot_prof[p,j]
            for i in i_litr_min:i_litr_max
                _scatter_add!(litr_c, c, j, i,
                    gru_leafc[p]*lf_f[ivtp,i]*wt*lp + gru_frootc[p]*fr_f[ivtp,i]*wt*fp)
                _scatter_add!(litr_n, c, j, i,
                    gru_leafn[p]*lf_f[ivtp,i]*wt*lp + gru_frootn[p]*fr_f[ivtp,i]*wt*fp)
            end
            _scatter_add!(cwdc, c, j, gru_livecrootc[p]*wt*cp)
            _scatter_add!(cwdc, c, j, gru_deadcrootc[p]*wt*cp)
            _scatter_add!(cwdn, c, j, gru_livecrootn[p]*wt*cp)
            _scatter_add!(cwdn, c, j, gru_deadcrootn[p]*wt*cp)
            _scatter_add!(litr_n, c, j, i_met_lit, gru_retransn[p]*wt*lp)
        end
    end
end

# ---------------------------------------------------------------------------
# cn_gross_unrep! — main kernel (= CNGrossUnrep)
# ---------------------------------------------------------------------------

"""
    cn_gross_unrep!(mask_soilp, state, pftcon, patch,
                    cnveg_cs, cnveg_ns, cnveg_cf, cnveg_nf, soilbiogeochem_state;
                    dt, noveg=0, nc4_grass=14, is_beg_curr_year=true,
                    nlevdecomp=1, i_litr_min=1, i_litr_max=3, i_met_lit=1)

Gross unrepresented landcover change "mortality" routine for the coupled
carbon-nitrogen code (CN). Computes the patch-level disturbance C/N fluxes from
each vegetation pool given the gross-unrep fraction, routing displayed leaf/root
pools to litter, wood to atmosphere and wood-product pools, and storage/transfer
pools to the atmosphere. Then aggregates to the column level via
`cn_gross_unrep_pft_to_column!`.

`mask_soilp` is the BitVector / Bool-vector soil-patch filter. `state` carries
`grossunrepfrac` and `do_grossunrep` from `dynGrossUnrep_interp!`. The natural
PFT range that receives gross unrep is `noveg < ivt <= nc4_grass`.

Only the non-matrix-CN path is ported (`use_matrixcn=.false.`).

Ported from `CNGrossUnrep` in `dynGrossUnrepMod.F90`.
"""
function cn_gross_unrep!(mask_soilp::AbstractVector{Bool},
                         state::DynGrossUnrepState,
                         pftcon::PftConGrossUnrep,
                         patch::PatchData,
                         cnveg_cs::CNVegCarbonStateData,
                         cnveg_ns::CNVegNitrogenStateData,
                         cnveg_cf::CNVegCarbonFluxData,
                         cnveg_nf::CNVegNitrogenFluxData,
                         soilbiogeochem_state::SoilBiogeochemStateData;
                         dt::Real = 1800.0,
                         noveg::Int = 0,
                         nc4_grass::Int = 14,
                         is_beg_curr_year::Bool = true,
                         nlevdecomp::Int = 1,
                         i_litr_min::Int = 1,
                         i_litr_max::Int = 3,
                         i_met_lit::Int = 1)

    dtime = dt

    # ---- aliases (Input: patch / pftcon / state) ----
    ivt            = patch.itype                  # pft vegetation type (0-based)
    gridcell       = patch.gridcell
    convfrac       = pftcon.pconv                 # Fortran convfrac => pftcon%pconv
    grossunrepfrac = state.grossunrepfrac
    do_grossunrep  = state.do_grossunrep

    # ---- aliases (Input: carbon state, displayed pools) ----
    leafc       = cnveg_cs.leafc_patch
    frootc      = cnveg_cs.frootc_patch
    livestemc   = cnveg_cs.livestemc_patch
    deadstemc   = cnveg_cs.deadstemc_patch
    livecrootc  = cnveg_cs.livecrootc_patch
    deadcrootc  = cnveg_cs.deadcrootc_patch
    xsmrpool    = cnveg_cs.xsmrpool_patch
    # storage pools
    leafc_storage      = cnveg_cs.leafc_storage_patch
    frootc_storage     = cnveg_cs.frootc_storage_patch
    livestemc_storage  = cnveg_cs.livestemc_storage_patch
    deadstemc_storage  = cnveg_cs.deadstemc_storage_patch
    livecrootc_storage = cnveg_cs.livecrootc_storage_patch
    deadcrootc_storage = cnveg_cs.deadcrootc_storage_patch
    gresp_storage      = cnveg_cs.gresp_storage_patch
    # transfer pools
    leafc_xfer      = cnveg_cs.leafc_xfer_patch
    frootc_xfer     = cnveg_cs.frootc_xfer_patch
    livestemc_xfer  = cnveg_cs.livestemc_xfer_patch
    deadstemc_xfer  = cnveg_cs.deadstemc_xfer_patch
    livecrootc_xfer = cnveg_cs.livecrootc_xfer_patch
    deadcrootc_xfer = cnveg_cs.deadcrootc_xfer_patch
    gresp_xfer      = cnveg_cs.gresp_xfer_patch

    # ---- aliases (Input: nitrogen state) ----
    leafn       = cnveg_ns.leafn_patch
    frootn      = cnveg_ns.frootn_patch
    livestemn   = cnveg_ns.livestemn_patch
    deadstemn   = cnveg_ns.deadstemn_patch
    livecrootn  = cnveg_ns.livecrootn_patch
    deadcrootn  = cnveg_ns.deadcrootn_patch
    retransn    = cnveg_ns.retransn_patch
    leafn_storage      = cnveg_ns.leafn_storage_patch
    frootn_storage     = cnveg_ns.frootn_storage_patch
    livestemn_storage  = cnveg_ns.livestemn_storage_patch
    deadstemn_storage  = cnveg_ns.deadstemn_storage_patch
    livecrootn_storage = cnveg_ns.livecrootn_storage_patch
    deadcrootn_storage = cnveg_ns.deadcrootn_storage_patch
    leafn_xfer      = cnveg_ns.leafn_xfer_patch
    frootn_xfer     = cnveg_ns.frootn_xfer_patch
    livestemn_xfer  = cnveg_ns.livestemn_xfer_patch
    deadstemn_xfer  = cnveg_ns.deadstemn_xfer_patch
    livecrootn_xfer = cnveg_ns.livecrootn_xfer_patch
    deadcrootn_xfer = cnveg_ns.deadcrootn_xfer_patch

    # ---- aliases (Output: carbon fluxes) ----
    gru_leafc_to_litter           = cnveg_cf.gru_leafc_to_litter_patch
    gru_frootc_to_litter          = cnveg_cf.gru_frootc_to_litter_patch
    gru_livestemc_to_atm          = cnveg_cf.gru_livestemc_to_atm_patch
    gru_deadstemc_to_atm          = cnveg_cf.gru_deadstemc_to_atm_patch
    gru_wood_productc_gain        = cnveg_cf.gru_wood_productc_gain_patch
    gru_livecrootc_to_litter      = cnveg_cf.gru_livecrootc_to_litter_patch
    gru_deadcrootc_to_litter      = cnveg_cf.gru_deadcrootc_to_litter_patch
    gru_xsmrpool_to_atm           = cnveg_cf.gru_xsmrpool_to_atm_patch
    gru_leafc_storage_to_atm      = cnveg_cf.gru_leafc_storage_to_atm_patch
    gru_frootc_storage_to_atm     = cnveg_cf.gru_frootc_storage_to_atm_patch
    gru_livestemc_storage_to_atm  = cnveg_cf.gru_livestemc_storage_to_atm_patch
    gru_deadstemc_storage_to_atm  = cnveg_cf.gru_deadstemc_storage_to_atm_patch
    gru_livecrootc_storage_to_atm = cnveg_cf.gru_livecrootc_storage_to_atm_patch
    gru_deadcrootc_storage_to_atm = cnveg_cf.gru_deadcrootc_storage_to_atm_patch
    gru_gresp_storage_to_atm      = cnveg_cf.gru_gresp_storage_to_atm_patch
    gru_leafc_xfer_to_atm         = cnveg_cf.gru_leafc_xfer_to_atm_patch
    gru_frootc_xfer_to_atm        = cnveg_cf.gru_frootc_xfer_to_atm_patch
    gru_livestemc_xfer_to_atm     = cnveg_cf.gru_livestemc_xfer_to_atm_patch
    gru_deadstemc_xfer_to_atm     = cnveg_cf.gru_deadstemc_xfer_to_atm_patch
    gru_livecrootc_xfer_to_atm    = cnveg_cf.gru_livecrootc_xfer_to_atm_patch
    gru_deadcrootc_xfer_to_atm    = cnveg_cf.gru_deadcrootc_xfer_to_atm_patch
    gru_gresp_xfer_to_atm         = cnveg_cf.gru_gresp_xfer_to_atm_patch

    # ---- aliases (Output: nitrogen fluxes) ----
    gru_leafn_to_litter           = cnveg_nf.gru_leafn_to_litter_patch
    gru_frootn_to_litter          = cnveg_nf.gru_frootn_to_litter_patch
    gru_livestemn_to_atm          = cnveg_nf.gru_livestemn_to_atm_patch
    gru_deadstemn_to_atm          = cnveg_nf.gru_deadstemn_to_atm_patch
    gru_wood_productn_gain        = cnveg_nf.gru_wood_productn_gain_patch
    gru_livecrootn_to_litter      = cnveg_nf.gru_livecrootn_to_litter_patch
    gru_deadcrootn_to_litter      = cnveg_nf.gru_deadcrootn_to_litter_patch
    gru_retransn_to_litter        = cnveg_nf.gru_retransn_to_litter_patch
    gru_leafn_storage_to_atm      = cnveg_nf.gru_leafn_storage_to_atm_patch
    gru_frootn_storage_to_atm     = cnveg_nf.gru_frootn_storage_to_atm_patch
    gru_livestemn_storage_to_atm  = cnveg_nf.gru_livestemn_storage_to_atm_patch
    gru_deadstemn_storage_to_atm  = cnveg_nf.gru_deadstemn_storage_to_atm_patch
    gru_livecrootn_storage_to_atm = cnveg_nf.gru_livecrootn_storage_to_atm_patch
    gru_deadcrootn_storage_to_atm = cnveg_nf.gru_deadcrootn_storage_to_atm_patch
    gru_leafn_xfer_to_atm         = cnveg_nf.gru_leafn_xfer_to_atm_patch
    gru_frootn_xfer_to_atm        = cnveg_nf.gru_frootn_xfer_to_atm_patch
    gru_livestemn_xfer_to_atm     = cnveg_nf.gru_livestemn_xfer_to_atm_patch
    gru_deadstemn_xfer_to_atm     = cnveg_nf.gru_deadstemn_xfer_to_atm_patch
    gru_livecrootn_xfer_to_atm    = cnveg_nf.gru_livecrootn_xfer_to_atm_patch
    gru_deadcrootn_xfer_to_atm    = cnveg_nf.gru_deadcrootn_xfer_to_atm_patch

    # compute the per-patch rate arrays (m, m*cfrac, m*(1-cfrac)) + natural-pft mask
    FT = eltype(leafc)
    np = length(mask_soilp)
    m   = fill!(similar(leafc, FT, np), zero(FT))
    mc  = fill!(similar(leafc, FT, np), zero(FT))
    mpr = fill!(similar(leafc, FT, np), zero(FT))
    gmask = fill!(similar(leafc, Bool, np), false)
    _launch!(_gru_m_kernel!, gmask, m, mc, mpr, mask_soilp, ivt, gridcell, grossunrepfrac,
        convfrac, is_beg_curr_year, do_grossunrep, FT(dtime), noveg, nc4_grass)
    # deadstem→atm scales by m*cfrac; wood-product gain by m*(1-cfrac); all else by m.
    for (o, inp) in ((gru_deadstemc_to_atm, deadstemc), (gru_deadstemn_to_atm, deadstemn))
        _launch!(_dynh_scale_kernel!, o, gmask, inp, mc)
    end
    for (o, inp) in ((gru_wood_productc_gain, deadstemc), (gru_wood_productn_gain, deadstemn))
        _launch!(_dynh_scale_kernel!, o, gmask, inp, mpr)
    end
    for (o, inp) in ((gru_leafc_to_litter, leafc), (gru_frootc_to_litter, frootc),
            (gru_livestemc_to_atm, livestemc), (gru_livecrootc_to_litter, livecrootc),
            (gru_deadcrootc_to_litter, deadcrootc), (gru_xsmrpool_to_atm, xsmrpool),
            (gru_leafc_storage_to_atm, leafc_storage), (gru_frootc_storage_to_atm, frootc_storage),
            (gru_livestemc_storage_to_atm, livestemc_storage),
            (gru_deadstemc_storage_to_atm, deadstemc_storage),
            (gru_livecrootc_storage_to_atm, livecrootc_storage),
            (gru_deadcrootc_storage_to_atm, deadcrootc_storage),
            (gru_gresp_storage_to_atm, gresp_storage), (gru_leafc_xfer_to_atm, leafc_xfer),
            (gru_frootc_xfer_to_atm, frootc_xfer), (gru_livestemc_xfer_to_atm, livestemc_xfer),
            (gru_deadstemc_xfer_to_atm, deadstemc_xfer),
            (gru_livecrootc_xfer_to_atm, livecrootc_xfer),
            (gru_deadcrootc_xfer_to_atm, deadcrootc_xfer), (gru_gresp_xfer_to_atm, gresp_xfer),
            (gru_leafn_to_litter, leafn), (gru_frootn_to_litter, frootn),
            (gru_livestemn_to_atm, livestemn), (gru_livecrootn_to_litter, livecrootn),
            (gru_deadcrootn_to_litter, deadcrootn), (gru_retransn_to_litter, retransn),
            (gru_leafn_storage_to_atm, leafn_storage), (gru_frootn_storage_to_atm, frootn_storage),
            (gru_livestemn_storage_to_atm, livestemn_storage),
            (gru_deadstemn_storage_to_atm, deadstemn_storage),
            (gru_livecrootn_storage_to_atm, livecrootn_storage),
            (gru_deadcrootn_storage_to_atm, deadcrootn_storage),
            (gru_leafn_xfer_to_atm, leafn_xfer), (gru_frootn_xfer_to_atm, frootn_xfer),
            (gru_livestemn_xfer_to_atm, livestemn_xfer),
            (gru_deadstemn_xfer_to_atm, deadstemn_xfer),
            (gru_livecrootn_xfer_to_atm, livecrootn_xfer),
            (gru_deadcrootn_xfer_to_atm, deadcrootn_xfer))
        _launch!(_dynh_scale_kernel!, o, gmask, inp, m)
    end

    # Gather all patch-level litterfall fluxes from grossunrep to the column for
    # litter C and N inputs.
    cn_gross_unrep_pft_to_column!(mask_soilp, pftcon, patch,
                                  cnveg_cf, cnveg_nf, soilbiogeochem_state;
                                  nlevdecomp = nlevdecomp,
                                  i_litr_min = i_litr_min,
                                  i_litr_max = i_litr_max,
                                  i_met_lit  = i_met_lit)

    return nothing
end

# ---------------------------------------------------------------------------
# cn_gross_unrep_pft_to_column! — patch -> column (= CNGrossUnrepPftToColumn)
# ---------------------------------------------------------------------------

"""
    cn_gross_unrep_pft_to_column!(mask_soilp, pftcon, patch,
                                  cnveg_cf, cnveg_nf, soilbiogeochem_state;
                                  nlevdecomp=1, i_litr_min=1, i_litr_max=3, i_met_lit=1)

Gather all patch-level gross-unrep litterfall fluxes to the column level and
assign them to the litter pools (leaf/root -> litter pools by litter fraction
and vertical profile), the coarse woody debris (CWD) pool (live/dead coarse
root), the metabolic litter pool (retranslocated N), and the wood-product gain
column accumulators.

The column accumulators are accumulated into (`+=`); the caller is responsible
for zeroing them.

Ported from `CNGrossUnrepPftToColumn` in `dynGrossUnrepMod.F90`.
"""
function cn_gross_unrep_pft_to_column!(mask_soilp::AbstractVector{Bool},
                                       pftcon::PftConGrossUnrep,
                                       patch::PatchData,
                                       cnveg_cf::CNVegCarbonFluxData,
                                       cnveg_nf::CNVegNitrogenFluxData,
                                       soilbiogeochem_state::SoilBiogeochemStateData;
                                       nlevdecomp::Int = 1,
                                       i_litr_min::Int = 1,
                                       i_litr_max::Int = 3,
                                       i_met_lit::Int = 1)

    # ---- aliases (Input) ----
    ivt   = patch.itype
    wtcol = patch.wtcol
    column = patch.column

    lf_f = pftcon.lf_f
    fr_f = pftcon.fr_f

    leaf_prof  = soilbiogeochem_state.leaf_prof_patch
    froot_prof = soilbiogeochem_state.froot_prof_patch
    croot_prof = soilbiogeochem_state.croot_prof_patch
    # stem_prof unused here (gross unrep routes wood C/N to CWD via croot_prof
    # and to product pools, matching the Fortran exactly).

    gru_leafc_to_litter      = cnveg_cf.gru_leafc_to_litter_patch
    gru_frootc_to_litter     = cnveg_cf.gru_frootc_to_litter_patch
    gru_livecrootc_to_litter = cnveg_cf.gru_livecrootc_to_litter_patch
    gru_deadcrootc_to_litter = cnveg_cf.gru_deadcrootc_to_litter_patch
    gru_wood_productc_gain   = cnveg_cf.gru_wood_productc_gain_patch

    gru_c_to_litr_c          = cnveg_cf.gru_c_to_litr_c_col        # (col, level, litr)
    gru_c_to_cwdc            = cnveg_cf.gru_c_to_cwdc_col          # (col, level)
    gru_wood_productc_gain_c = cnveg_cf.gru_wood_productc_gain_col # (col)

    gru_leafn_to_litter      = cnveg_nf.gru_leafn_to_litter_patch
    gru_frootn_to_litter     = cnveg_nf.gru_frootn_to_litter_patch
    gru_livecrootn_to_litter = cnveg_nf.gru_livecrootn_to_litter_patch
    gru_deadcrootn_to_litter = cnveg_nf.gru_deadcrootn_to_litter_patch
    gru_retransn_to_litter   = cnveg_nf.gru_retransn_to_litter_patch
    gru_wood_productn_gain   = cnveg_nf.gru_wood_productn_gain_patch

    gru_n_to_litr_n          = cnveg_nf.gru_n_to_litr_n_col        # (col, level, litr)
    gru_n_to_cwdn            = cnveg_nf.gru_n_to_cwdn_col          # (col, level)
    gru_wood_productn_gain_c = cnveg_nf.gru_wood_productn_gain_col # (col)

    _launch!(_gru_scatter_kernel!, mask_soilp, column, ivt, leaf_prof, froot_prof, croot_prof,
        lf_f, fr_f, wtcol, gru_leafc_to_litter, gru_frootc_to_litter, gru_livecrootc_to_litter,
        gru_deadcrootc_to_litter, gru_leafn_to_litter, gru_frootn_to_litter,
        gru_livecrootn_to_litter, gru_deadcrootn_to_litter, gru_retransn_to_litter,
        gru_c_to_litr_c, gru_c_to_cwdc, gru_n_to_litr_n, gru_n_to_cwdn,
        nlevdecomp, i_litr_min, i_litr_max, i_met_lit)

    # wood gross-unrep mortality C & N to the product pools (generic p2c scatter)
    np = length(mask_soilp)
    ndyn_p2c_scatter!(gru_wood_productc_gain_c, mask_soilp, column, gru_wood_productc_gain, wtcol, 1, np)
    ndyn_p2c_scatter!(gru_wood_productn_gain_c, mask_soilp, column, gru_wood_productn_gain, wtcol, 1, np)

    return nothing
end
