# ==========================================================================
# Ported from: src/biogeochem/CNFireNoFireMod.F90
# Fire dynamics module with fire explicitly turned off (the "nofire" method).
#
# Public functions:
#   need_lightning_and_popdens_nofire  — Returns false
#   cnfire_area_nofire!                — Zeros all burned-area state
#   cnfire_fluxes_nofire!              — No fire C/N fluxes (no-op)
# ==========================================================================

# ---------------------------------------------------------------------------
# need_lightning_and_popdens_nofire
# ---------------------------------------------------------------------------

"""
    need_lightning_and_popdens_nofire()

Returns `false` — the NoFire model needs no lightning / population density data.
Ported from `need_lightning_and_popdens` in `CNFireNoFireMod.F90`.
"""
need_lightning_and_popdens_nofire() = false

# ---------------------------------------------------------------------------
# cnfire_area_nofire! — zero all burned area
# ---------------------------------------------------------------------------

"""
    cnfire_area_nofire!(cnveg_state, cnveg_cs, patch, mask_soilc, bounds_c, bounds_p)

Compute column-level burned area for the NoFire method: a leafc patch-to-column
average (kept for diagnostic parity with the Fortran), then every fire-area
state variable is zeroed. With NoFire, tree carbon is still removed in land-use
change regions by the land-use code (not here).

Ported from `CNFireArea` in `CNFireNoFireMod.F90`.
"""
function cnfire_area_nofire!(
    cnveg_state::CNVegStateData,
    cnveg_cs::CNVegCarbonStateData,
    patch::PatchData,
    mask_soilc::AbstractVector{Bool},
    bounds_c::UnitRange{Int},
    bounds_p::UnitRange{Int},
)
    farea_burned = cnveg_state.farea_burned_col
    baf_crop     = cnveg_state.baf_crop_col
    baf_peatf    = cnveg_state.baf_peatf_col
    fbac         = cnveg_state.fbac_col
    fbac1        = cnveg_state.fbac1_col
    cropf_col    = cnveg_state.cropf_col
    lfc          = cnveg_state.lfc_col

    leafc        = cnveg_cs.leafc_patch
    leafc_col    = cnveg_cs.leafc_col

    # pft-to-column average of leafc (matches the Fortran p2c call)
    p2c!(leafc_col, leafc, patch, mask_soilc, bounds_c, bounds_p)

    for c in bounds_c
        mask_soilc[c] || continue
        farea_burned[c] = 0.0
        baf_crop[c]     = 0.0
        baf_peatf[c]    = 0.0
        fbac[c]         = 0.0
        fbac1[c]        = 0.0
        cropf_col[c]    = 0.0
        lfc[c]          = 0.0
    end

    return nothing
end

# ---------------------------------------------------------------------------
# cnfire_fluxes_nofire! — no fire fluxes
# ---------------------------------------------------------------------------

"""
    cnfire_fluxes_nofire!(...)

Fire effects routine for the NoFire method. The base `CNFireFluxes` is still
inherited in Fortran (NoFire only overrides `CNFireArea`, zeroing burned area),
so with `farea_burned == 0` every fire flux evaluates to zero. We delegate to
the shared `cnfire_fluxes!` with the standard combustion-completeness factors;
because burned area is zero the result is no fire fluxes.

Ported from the inherited `CNFireFluxes` behavior under `CNFireNoFireMod.F90`.
"""
function cnfire_fluxes_nofire!(
    mask_soilc::AbstractVector{Bool},
    mask_soilp::AbstractVector{Bool},
    bounds_c::UnitRange{Int},
    bounds_p::UnitRange{Int},
    cnfire_const::CNFireConstData,
    pftcon::PftConFireBase,
    patch::PatchData,
    col::ColumnData,
    grc::GridcellData,
    dgvs::DgvsFireData,
    cnveg_state::CNVegStateData,
    cnveg_cs::CNVegCarbonStateData,
    cnveg_cf::CNVegCarbonFluxData,
    cnveg_ns::CNVegNitrogenStateData,
    cnveg_nf::CNVegNitrogenFluxData,
    soilbgc_cf::SoilBiogeochemCarbonFluxData,
    decomp_cascade_con::DecompCascadeConData,
    leaf_prof_patch::AbstractMatrix{<:Real},
    froot_prof_patch::AbstractMatrix{<:Real},
    croot_prof_patch::AbstractMatrix{<:Real},
    stem_prof_patch::AbstractMatrix{<:Real},
    totsomc_col::AbstractVector{<:Real},
    decomp_cpools_vr_col::AbstractArray{<:Real,3},
    decomp_npools_vr_col::AbstractArray{<:Real,3},
    somc_fire_col::AbstractVector{<:Real};
    kwargs...
)
    return cnfire_fluxes!(
        mask_soilc, mask_soilp, bounds_c, bounds_p,
        cnfire_const, pftcon, patch, col, grc,
        dgvs, cnveg_state, cnveg_cs, cnveg_cf, cnveg_ns, cnveg_nf,
        soilbgc_cf, decomp_cascade_con,
        leaf_prof_patch, froot_prof_patch, croot_prof_patch, stem_prof_patch,
        totsomc_col, decomp_cpools_vr_col, decomp_npools_vr_col, somc_fire_col;
        kwargs...)
end
