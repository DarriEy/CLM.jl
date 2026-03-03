# ==========================================================================
# Ported from: src/biogeochem/CNCIsoFluxMod.F90
# Module for carbon isotopic flux variable update, non-mortality fluxes.
#
# Public functions:
#   c_iso_flux_calc!        — Core 1D isotopic flux calculation
#   c_iso_flux_calc_2d_flux! — 2D wrapper (flux is 2D, state is 1D)
#   c_iso_flux_calc_2d_both! — 2D wrapper (both flux and state are 2D)
#   c_iso_flux1!            — Non-mortality isotopic fluxes
#   c_iso_flux2!            — Gap mortality isotopic fluxes
#   c_iso_flux2h!           — Harvest mortality isotopic fluxes
#   c_iso_flux2g!           — Gross unrepresented landcover change isotopic fluxes
#   c_iso_flux3!            — Fire mortality isotopic fluxes
#
# Private functions:
#   cn_c_iso_litter_to_column!       — Phenology litterfall patch-to-column
#   cn_c_iso_gap_pft_to_column!      — Gap mortality patch-to-column
#   cn_c_iso_harvest_pft_to_column!  — Harvest mortality patch-to-column
#   cn_c_iso_gross_unrep_pft_to_column! — Gross unrep. LC change patch-to-column
# ==========================================================================

# ---------------------------------------------------------------------------
# c_iso_flux_calc! — Core 1D isotopic flux calculation
# Ported from: CIsoFluxCalc1d in CNCIsoFluxMod.F90
# ---------------------------------------------------------------------------

"""
    c_iso_flux_calc!(ciso_flux, ctot_flux, ciso_state, ctot_state,
                     mask, bounds, frax_c13, isotope)

Set the carbon isotopic flux for a 1D variable. For each index in `bounds`
where `mask` is true, compute:

    ciso_flux[i] = ctot_flux[i] * (ciso_state[i] / ctot_state[i]) * frax

where `frax` is the fractionation factor (doubled deviation for C14).

Ported from `CIsoFluxCalc1d` in `CNCIsoFluxMod.F90`.
"""
function c_iso_flux_calc!(ciso_flux::AbstractVector{Float64},
                          ctot_flux::AbstractVector{Float64},
                          ciso_state::AbstractVector{Float64},
                          ctot_state::AbstractVector{Float64},
                          mask::BitVector,
                          bounds::UnitRange{Int},
                          frax_c13::Float64,
                          isotope::String)
    # if C14, double the fractionation
    if isotope == "c14"
        frax = 1.0 + (1.0 - frax_c13) * 2.0
    elseif isotope == "c13"
        frax = frax_c13
    else
        error("c_iso_flux_calc!: isotope must be either \"c13\" or \"c14\", got \"$isotope\"")
    end

    for i in bounds
        mask[i] || continue
        if ctot_state[i] != 0.0 && ciso_state[i] != 0.0
            ciso_flux[i] = ctot_flux[i] * (ciso_state[i] / ctot_state[i]) * frax
        else
            ciso_flux[i] = 0.0
        end
    end
    return nothing
end

# ---------------------------------------------------------------------------
# c_iso_flux_calc_2d_flux! — 2D wrapper (flux is 2D, state is 1D)
# Ported from: CIsoFluxCalc2dFlux in CNCIsoFluxMod.F90
# ---------------------------------------------------------------------------

"""
    c_iso_flux_calc_2d_flux!(ciso_flux, ctot_flux, ciso_state, ctot_state,
                             mask, bounds, frax_c13, isotope)

Wrapper for `c_iso_flux_calc!` when the flux arrays are 2D (first dim × second dim)
but the state arrays are 1D. Loops over the second dimension.

Ported from `CIsoFluxCalc2dFlux` in `CNCIsoFluxMod.F90`.
"""
function c_iso_flux_calc_2d_flux!(ciso_flux::AbstractMatrix{Float64},
                                  ctot_flux::AbstractMatrix{Float64},
                                  ciso_state::AbstractVector{Float64},
                                  ctot_state::AbstractVector{Float64},
                                  mask::BitVector,
                                  bounds::UnitRange{Int},
                                  frax_c13::Float64,
                                  isotope::String)
    for k in axes(ciso_flux, 2)
        ciso_flux_1d = @view ciso_flux[:, k]
        ctot_flux_1d = @view ctot_flux[:, k]
        c_iso_flux_calc!(ciso_flux_1d, ctot_flux_1d,
                         ciso_state, ctot_state,
                         mask, bounds, frax_c13, isotope)
    end
    return nothing
end

# ---------------------------------------------------------------------------
# c_iso_flux_calc_2d_both! — 2D wrapper (both flux and state are 2D)
# Ported from: CIsoFluxCalc2dBoth in CNCIsoFluxMod.F90
# ---------------------------------------------------------------------------

"""
    c_iso_flux_calc_2d_both!(ciso_flux, ctot_flux, ciso_state, ctot_state,
                              mask, bounds, frax_c13, isotope)

Wrapper for `c_iso_flux_calc!` when both flux and state arrays are 2D.
Loops over the second dimension of the flux arrays.

Ported from `CIsoFluxCalc2dBoth` in `CNCIsoFluxMod.F90`.
"""
function c_iso_flux_calc_2d_both!(ciso_flux::AbstractMatrix{Float64},
                                   ctot_flux::AbstractMatrix{Float64},
                                   ciso_state::AbstractMatrix{Float64},
                                   ctot_state::AbstractMatrix{Float64},
                                   mask::BitVector,
                                   bounds::UnitRange{Int},
                                   frax_c13::Float64,
                                   isotope::String)
    for k in axes(ciso_flux, 2)
        ciso_flux_1d = @view ciso_flux[:, k]
        ctot_flux_1d = @view ctot_flux[:, k]
        ciso_state_1d = @view ciso_state[:, k]
        ctot_state_1d = @view ctot_state[:, k]
        c_iso_flux_calc!(ciso_flux_1d, ctot_flux_1d,
                         ciso_state_1d, ctot_state_1d,
                         mask, bounds, frax_c13, isotope)
    end
    return nothing
end

# ---------------------------------------------------------------------------
# c_iso_flux1! — Non-mortality isotopic fluxes
# Ported from: CIsoFlux1 in CNCIsoFluxMod.F90
# ---------------------------------------------------------------------------

"""
    c_iso_flux1!(soilbiogeochem_state, soilbiogeochem_cf, soilbiogeochem_cs,
                 cnveg_cf, cnveg_cs,
                 iso_soilbiogeochem_cf, iso_soilbiogeochem_cs,
                 iso_cnveg_cf, iso_cnveg_cs,
                 mask_soilc, mask_soilp, bounds_c, bounds_p,
                 cascade_donor_pool, nlevdecomp, ndecomp_cascade_transitions,
                 isotope; use_crop, use_grainproduct, nrepr, npcropmin,
                 repr_grain_min, repr_grain_max, repr_structure_min, repr_structure_max,
                 patch_column, patch_itype, patch_wtcol, lf_f, fr_f,
                 leaf_prof, froot_prof)

On the radiation time step, set the carbon isotopic flux variables
(except for gap-phase mortality and fire fluxes).

Ported from `CIsoFlux1` in `CNCIsoFluxMod.F90`.
"""
function c_iso_flux1!(soilbiogeochem_state::SoilBiogeochemStateData,
                      soilbiogeochem_cf::SoilBiogeochemCarbonFluxData,
                      soilbiogeochem_cs::SoilBiogeochemCarbonStateData,
                      cnveg_cf::CNVegCarbonFluxData,
                      cnveg_cs::CNVegCarbonStateData,
                      iso_soilbiogeochem_cf::SoilBiogeochemCarbonFluxData,
                      iso_soilbiogeochem_cs::SoilBiogeochemCarbonStateData,
                      iso_cnveg_cf::CNVegCarbonFluxData,
                      iso_cnveg_cs::CNVegCarbonStateData,
                      mask_soilc::BitVector,
                      mask_soilp::BitVector,
                      bounds_c::UnitRange{Int},
                      bounds_p::UnitRange{Int},
                      cascade_donor_pool::Vector{Int},
                      nlevdecomp::Int,
                      ndecomp_cascade_transitions::Int,
                      isotope::String;
                      use_crop::Bool=false,
                      use_grainproduct::Bool=false,
                      nrepr::Int=NREPR,
                      npcropmin::Int=NPCROPMIN,
                      repr_grain_min::Int=1,
                      repr_grain_max::Int=1,
                      repr_structure_min::Int=1,
                      repr_structure_max::Int=0,
                      patch_column::Vector{Int}=Int[],
                      patch_itype::Vector{Int}=Int[],
                      patch_wtcol::Vector{Float64}=Float64[],
                      lf_f::Matrix{Float64}=Matrix{Float64}(undef, 0, 0),
                      fr_f::Matrix{Float64}=Matrix{Float64}(undef, 0, 0),
                      leaf_prof::Matrix{Float64}=Matrix{Float64}(undef, 0, 0),
                      froot_prof::Matrix{Float64}=Matrix{Float64}(undef, 0, 0))

    # Helper for 1D patch-level iso flux calc
    calc1d! = (iso_f, tot_f, iso_s, tot_s) -> c_iso_flux_calc!(
        iso_f, tot_f, iso_s, tot_s, mask_soilp, bounds_p, 1.0, isotope)

    # Helper for 2D flux (flux 2D, state 1D)
    calc2df! = (iso_f, tot_f, iso_s, tot_s) -> c_iso_flux_calc_2d_flux!(
        iso_f, tot_f, iso_s, tot_s, mask_soilp, bounds_p, 1.0, isotope)

    # Helper for 2D both (flux 2D, state 2D)
    calc2db! = (iso_f, tot_f, iso_s, tot_s) -> c_iso_flux_calc_2d_both!(
        iso_f, tot_f, iso_s, tot_s, mask_soilp, bounds_p, 1.0, isotope)

    # --- patch-level non-mortality fluxes ---

    calc1d!(iso_cnveg_cf.leafc_xfer_to_leafc_patch, cnveg_cf.leafc_xfer_to_leafc_patch,
            iso_cnveg_cs.leafc_xfer_patch, cnveg_cs.leafc_xfer_patch)

    calc1d!(iso_cnveg_cf.frootc_xfer_to_frootc_patch, cnveg_cf.frootc_xfer_to_frootc_patch,
            iso_cnveg_cs.frootc_xfer_patch, cnveg_cs.frootc_xfer_patch)

    calc1d!(iso_cnveg_cf.livestemc_xfer_to_livestemc_patch, cnveg_cf.livestemc_xfer_to_livestemc_patch,
            iso_cnveg_cs.livestemc_xfer_patch, cnveg_cs.livestemc_xfer_patch)

    calc1d!(iso_cnveg_cf.deadstemc_xfer_to_deadstemc_patch, cnveg_cf.deadstemc_xfer_to_deadstemc_patch,
            iso_cnveg_cs.deadstemc_xfer_patch, cnveg_cs.deadstemc_xfer_patch)

    calc1d!(iso_cnveg_cf.livecrootc_xfer_to_livecrootc_patch, cnveg_cf.livecrootc_xfer_to_livecrootc_patch,
            iso_cnveg_cs.livecrootc_xfer_patch, cnveg_cs.livecrootc_xfer_patch)

    calc1d!(iso_cnveg_cf.deadcrootc_xfer_to_deadcrootc_patch, cnveg_cf.deadcrootc_xfer_to_deadcrootc_patch,
            iso_cnveg_cs.deadcrootc_xfer_patch, cnveg_cs.deadcrootc_xfer_patch)

    calc1d!(iso_cnveg_cf.leafc_to_litter_patch, cnveg_cf.leafc_to_litter_patch,
            iso_cnveg_cs.leafc_patch, cnveg_cs.leafc_patch)

    calc1d!(iso_cnveg_cf.frootc_to_litter_patch, cnveg_cf.frootc_to_litter_patch,
            iso_cnveg_cs.frootc_patch, cnveg_cs.frootc_patch)

    calc1d!(iso_cnveg_cf.livestemc_to_deadstemc_patch, cnveg_cf.livestemc_to_deadstemc_patch,
            iso_cnveg_cs.livestemc_patch, cnveg_cs.livestemc_patch)

    calc1d!(iso_cnveg_cf.livecrootc_to_deadcrootc_patch, cnveg_cf.livecrootc_to_deadcrootc_patch,
            iso_cnveg_cs.livecrootc_patch, cnveg_cs.livecrootc_patch)

    # Maintenance respiration from cpool
    calc1d!(iso_cnveg_cf.leaf_curmr_patch, cnveg_cf.leaf_curmr_patch,
            iso_cnveg_cs.cpool_patch, cnveg_cs.cpool_patch)
    calc1d!(iso_cnveg_cf.froot_curmr_patch, cnveg_cf.froot_curmr_patch,
            iso_cnveg_cs.cpool_patch, cnveg_cs.cpool_patch)
    calc1d!(iso_cnveg_cf.livestem_curmr_patch, cnveg_cf.livestem_curmr_patch,
            iso_cnveg_cs.cpool_patch, cnveg_cs.cpool_patch)
    calc1d!(iso_cnveg_cf.livecroot_curmr_patch, cnveg_cf.livecroot_curmr_patch,
            iso_cnveg_cs.cpool_patch, cnveg_cs.cpool_patch)

    # Excess maintenance respiration from totvegc
    calc1d!(iso_cnveg_cf.leaf_xsmr_patch, cnveg_cf.leaf_xsmr_patch,
            iso_cnveg_cs.totvegc_patch, cnveg_cs.totvegc_patch)
    calc1d!(iso_cnveg_cf.froot_xsmr_patch, cnveg_cf.froot_xsmr_patch,
            iso_cnveg_cs.totvegc_patch, cnveg_cs.totvegc_patch)
    calc1d!(iso_cnveg_cf.livestem_xsmr_patch, cnveg_cf.livestem_xsmr_patch,
            iso_cnveg_cs.totvegc_patch, cnveg_cs.totvegc_patch)
    calc1d!(iso_cnveg_cf.livecroot_xsmr_patch, cnveg_cf.livecroot_xsmr_patch,
            iso_cnveg_cs.totvegc_patch, cnveg_cs.totvegc_patch)

    # cpool allocations
    calc1d!(iso_cnveg_cf.cpool_to_xsmrpool_patch, cnveg_cf.cpool_to_xsmrpool_patch,
            iso_cnveg_cs.cpool_patch, cnveg_cs.cpool_patch)
    calc1d!(iso_cnveg_cf.cpool_to_leafc_patch, cnveg_cf.cpool_to_leafc_patch,
            iso_cnveg_cs.cpool_patch, cnveg_cs.cpool_patch)
    calc1d!(iso_cnveg_cf.cpool_to_leafc_storage_patch, cnveg_cf.cpool_to_leafc_storage_patch,
            iso_cnveg_cs.cpool_patch, cnveg_cs.cpool_patch)
    calc1d!(iso_cnveg_cf.cpool_to_frootc_patch, cnveg_cf.cpool_to_frootc_patch,
            iso_cnveg_cs.cpool_patch, cnveg_cs.cpool_patch)
    calc1d!(iso_cnveg_cf.cpool_to_frootc_storage_patch, cnveg_cf.cpool_to_frootc_storage_patch,
            iso_cnveg_cs.cpool_patch, cnveg_cs.cpool_patch)
    calc1d!(iso_cnveg_cf.cpool_to_livestemc_patch, cnveg_cf.cpool_to_livestemc_patch,
            iso_cnveg_cs.cpool_patch, cnveg_cs.cpool_patch)
    calc1d!(iso_cnveg_cf.cpool_to_livestemc_storage_patch, cnveg_cf.cpool_to_livestemc_storage_patch,
            iso_cnveg_cs.cpool_patch, cnveg_cs.cpool_patch)
    calc1d!(iso_cnveg_cf.cpool_to_deadstemc_patch, cnveg_cf.cpool_to_deadstemc_patch,
            iso_cnveg_cs.cpool_patch, cnveg_cs.cpool_patch)
    calc1d!(iso_cnveg_cf.cpool_to_deadstemc_storage_patch, cnveg_cf.cpool_to_deadstemc_storage_patch,
            iso_cnveg_cs.cpool_patch, cnveg_cs.cpool_patch)
    calc1d!(iso_cnveg_cf.cpool_to_livecrootc_patch, cnveg_cf.cpool_to_livecrootc_patch,
            iso_cnveg_cs.cpool_patch, cnveg_cs.cpool_patch)
    calc1d!(iso_cnveg_cf.cpool_to_livecrootc_storage_patch, cnveg_cf.cpool_to_livecrootc_storage_patch,
            iso_cnveg_cs.cpool_patch, cnveg_cs.cpool_patch)
    calc1d!(iso_cnveg_cf.cpool_to_deadcrootc_patch, cnveg_cf.cpool_to_deadcrootc_patch,
            iso_cnveg_cs.cpool_patch, cnveg_cs.cpool_patch)
    calc1d!(iso_cnveg_cf.cpool_to_deadcrootc_storage_patch, cnveg_cf.cpool_to_deadcrootc_storage_patch,
            iso_cnveg_cs.cpool_patch, cnveg_cs.cpool_patch)

    # Growth respiration from cpool
    calc1d!(iso_cnveg_cf.cpool_leaf_gr_patch, cnveg_cf.cpool_leaf_gr_patch,
            iso_cnveg_cs.cpool_patch, cnveg_cs.cpool_patch)
    calc1d!(iso_cnveg_cf.cpool_froot_gr_patch, cnveg_cf.cpool_froot_gr_patch,
            iso_cnveg_cs.cpool_patch, cnveg_cs.cpool_patch)
    calc1d!(iso_cnveg_cf.cpool_livestem_gr_patch, cnveg_cf.cpool_livestem_gr_patch,
            iso_cnveg_cs.cpool_patch, cnveg_cs.cpool_patch)
    calc1d!(iso_cnveg_cf.cpool_deadstem_gr_patch, cnveg_cf.cpool_deadstem_gr_patch,
            iso_cnveg_cs.cpool_patch, cnveg_cs.cpool_patch)
    calc1d!(iso_cnveg_cf.cpool_livecroot_gr_patch, cnveg_cf.cpool_livecroot_gr_patch,
            iso_cnveg_cs.cpool_patch, cnveg_cs.cpool_patch)
    calc1d!(iso_cnveg_cf.cpool_deadcroot_gr_patch, cnveg_cf.cpool_deadcroot_gr_patch,
            iso_cnveg_cs.cpool_patch, cnveg_cs.cpool_patch)

    # Storage growth respiration from cpool
    calc1d!(iso_cnveg_cf.cpool_leaf_storage_gr_patch, cnveg_cf.cpool_leaf_storage_gr_patch,
            iso_cnveg_cs.cpool_patch, cnveg_cs.cpool_patch)
    calc1d!(iso_cnveg_cf.cpool_froot_storage_gr_patch, cnveg_cf.cpool_froot_storage_gr_patch,
            iso_cnveg_cs.cpool_patch, cnveg_cs.cpool_patch)
    calc1d!(iso_cnveg_cf.cpool_livestem_storage_gr_patch, cnveg_cf.cpool_livestem_storage_gr_patch,
            iso_cnveg_cs.cpool_patch, cnveg_cs.cpool_patch)
    calc1d!(iso_cnveg_cf.cpool_deadstem_storage_gr_patch, cnveg_cf.cpool_deadstem_storage_gr_patch,
            iso_cnveg_cs.cpool_patch, cnveg_cs.cpool_patch)
    calc1d!(iso_cnveg_cf.cpool_livecroot_storage_gr_patch, cnveg_cf.cpool_livecroot_storage_gr_patch,
            iso_cnveg_cs.cpool_patch, cnveg_cs.cpool_patch)
    calc1d!(iso_cnveg_cf.cpool_deadcroot_storage_gr_patch, cnveg_cf.cpool_deadcroot_storage_gr_patch,
            iso_cnveg_cs.cpool_patch, cnveg_cs.cpool_patch)

    calc1d!(iso_cnveg_cf.cpool_to_gresp_storage_patch, cnveg_cf.cpool_to_gresp_storage_patch,
            iso_cnveg_cs.cpool_patch, cnveg_cs.cpool_patch)

    # Transfer growth respiration from gresp_xfer
    calc1d!(iso_cnveg_cf.transfer_leaf_gr_patch, cnveg_cf.transfer_leaf_gr_patch,
            iso_cnveg_cs.gresp_xfer_patch, cnveg_cs.gresp_xfer_patch)
    calc1d!(iso_cnveg_cf.transfer_froot_gr_patch, cnveg_cf.transfer_froot_gr_patch,
            iso_cnveg_cs.gresp_xfer_patch, cnveg_cs.gresp_xfer_patch)
    calc1d!(iso_cnveg_cf.transfer_livestem_gr_patch, cnveg_cf.transfer_livestem_gr_patch,
            iso_cnveg_cs.gresp_xfer_patch, cnveg_cs.gresp_xfer_patch)
    calc1d!(iso_cnveg_cf.transfer_deadstem_gr_patch, cnveg_cf.transfer_deadstem_gr_patch,
            iso_cnveg_cs.gresp_xfer_patch, cnveg_cs.gresp_xfer_patch)
    calc1d!(iso_cnveg_cf.transfer_livecroot_gr_patch, cnveg_cf.transfer_livecroot_gr_patch,
            iso_cnveg_cs.gresp_xfer_patch, cnveg_cs.gresp_xfer_patch)
    calc1d!(iso_cnveg_cf.transfer_deadcroot_gr_patch, cnveg_cf.transfer_deadcroot_gr_patch,
            iso_cnveg_cs.gresp_xfer_patch, cnveg_cs.gresp_xfer_patch)

    # Storage to transfer
    calc1d!(iso_cnveg_cf.leafc_storage_to_xfer_patch, cnveg_cf.leafc_storage_to_xfer_patch,
            iso_cnveg_cs.leafc_storage_patch, cnveg_cs.leafc_storage_patch)
    calc1d!(iso_cnveg_cf.frootc_storage_to_xfer_patch, cnveg_cf.frootc_storage_to_xfer_patch,
            iso_cnveg_cs.frootc_storage_patch, cnveg_cs.frootc_storage_patch)
    calc1d!(iso_cnveg_cf.livestemc_storage_to_xfer_patch, cnveg_cf.livestemc_storage_to_xfer_patch,
            iso_cnveg_cs.livestemc_storage_patch, cnveg_cs.livestemc_storage_patch)
    calc1d!(iso_cnveg_cf.deadstemc_storage_to_xfer_patch, cnveg_cf.deadstemc_storage_to_xfer_patch,
            iso_cnveg_cs.deadstemc_storage_patch, cnveg_cs.deadstemc_storage_patch)
    calc1d!(iso_cnveg_cf.livecrootc_storage_to_xfer_patch, cnveg_cf.livecrootc_storage_to_xfer_patch,
            iso_cnveg_cs.livecrootc_storage_patch, cnveg_cs.livecrootc_storage_patch)
    calc1d!(iso_cnveg_cf.deadcrootc_storage_to_xfer_patch, cnveg_cf.deadcrootc_storage_to_xfer_patch,
            iso_cnveg_cs.deadcrootc_storage_patch, cnveg_cs.deadcrootc_storage_patch)
    calc1d!(iso_cnveg_cf.gresp_storage_to_xfer_patch, cnveg_cf.gresp_storage_to_xfer_patch,
            iso_cnveg_cs.gresp_storage_patch, cnveg_cs.gresp_storage_patch)

    # Misc
    calc1d!(iso_cnveg_cf.soilc_change_patch, cnveg_cf.soilc_change_patch,
            iso_cnveg_cs.cpool_patch, cnveg_cs.cpool_patch)
    calc1d!(iso_cnveg_cf.cpool_to_resp_patch, cnveg_cf.cpool_to_resp_patch,
            iso_cnveg_cs.cpool_patch, cnveg_cs.cpool_patch)

    # --- Crop-specific fluxes ---
    if use_crop
        calc2db!(iso_cnveg_cf.reproductivec_xfer_to_reproductivec_patch,
                 cnveg_cf.reproductivec_xfer_to_reproductivec_patch,
                 iso_cnveg_cs.reproductivec_xfer_patch, cnveg_cs.reproductivec_xfer_patch)

        calc2db!(iso_cnveg_cf.repr_grainc_to_food_patch, cnveg_cf.repr_grainc_to_food_patch,
                 iso_cnveg_cs.reproductivec_patch, cnveg_cs.reproductivec_patch)

        calc1d!(iso_cnveg_cf.leafc_to_biofuelc_patch, cnveg_cf.leafc_to_biofuelc_patch,
                iso_cnveg_cs.leafc_patch, cnveg_cs.leafc_patch)

        calc1d!(iso_cnveg_cf.livestemc_to_biofuelc_patch, cnveg_cf.livestemc_to_biofuelc_patch,
                iso_cnveg_cs.livestemc_patch, cnveg_cs.livestemc_patch)

        calc1d!(iso_cnveg_cf.leafc_to_removedresiduec_patch, cnveg_cf.leafc_to_removedresiduec_patch,
                iso_cnveg_cs.leafc_patch, cnveg_cs.leafc_patch)

        calc1d!(iso_cnveg_cf.livestemc_to_removedresiduec_patch, cnveg_cf.livestemc_to_removedresiduec_patch,
                iso_cnveg_cs.livestemc_patch, cnveg_cs.livestemc_patch)

        calc2db!(iso_cnveg_cf.repr_grainc_to_seed_patch, cnveg_cf.repr_grainc_to_seed_patch,
                 iso_cnveg_cs.reproductivec_patch, cnveg_cs.reproductivec_patch)

        calc2db!(iso_cnveg_cf.repr_structurec_to_cropprod_patch, cnveg_cf.repr_structurec_to_cropprod_patch,
                 iso_cnveg_cs.reproductivec_patch, cnveg_cs.reproductivec_patch)

        calc2db!(iso_cnveg_cf.repr_structurec_to_litter_patch, cnveg_cf.repr_structurec_to_litter_patch,
                 iso_cnveg_cs.reproductivec_patch, cnveg_cs.reproductivec_patch)

        calc1d!(iso_cnveg_cf.crop_seedc_to_leaf_patch, cnveg_cf.crop_seedc_to_leaf_patch,
                iso_cnveg_cs.totvegc_patch, cnveg_cs.totvegc_patch)

        calc2df!(iso_cnveg_cf.reproductive_curmr_patch, cnveg_cf.reproductive_curmr_patch,
                 iso_cnveg_cs.cpool_patch, cnveg_cs.cpool_patch)

        calc2df!(iso_cnveg_cf.reproductive_xsmr_patch, cnveg_cf.reproductive_xsmr_patch,
                 iso_cnveg_cs.totvegc_patch, cnveg_cs.totvegc_patch)

        calc2df!(iso_cnveg_cf.cpool_reproductive_gr_patch, cnveg_cf.cpool_reproductive_gr_patch,
                 iso_cnveg_cs.cpool_patch, cnveg_cs.cpool_patch)

        calc2df!(iso_cnveg_cf.cpool_to_reproductivec_patch, cnveg_cf.cpool_to_reproductivec_patch,
                 iso_cnveg_cs.cpool_patch, cnveg_cs.cpool_patch)

        calc2df!(iso_cnveg_cf.cpool_to_reproductivec_storage_patch, cnveg_cf.cpool_to_reproductivec_storage_patch,
                 iso_cnveg_cs.cpool_patch, cnveg_cs.cpool_patch)

        calc2df!(iso_cnveg_cf.transfer_reproductive_gr_patch, cnveg_cf.transfer_reproductive_gr_patch,
                 iso_cnveg_cs.gresp_xfer_patch, cnveg_cs.gresp_xfer_patch)

        calc2df!(iso_cnveg_cf.cpool_reproductive_storage_gr_patch, cnveg_cf.cpool_reproductive_storage_gr_patch,
                 iso_cnveg_cs.cpool_patch, cnveg_cs.cpool_patch)

        calc2db!(iso_cnveg_cf.reproductivec_storage_to_xfer_patch, cnveg_cf.reproductivec_storage_to_xfer_patch,
                 iso_cnveg_cs.reproductivec_storage_patch, cnveg_cs.reproductivec_storage_patch)

        calc1d!(iso_cnveg_cf.livestemc_to_litter_patch, cnveg_cf.livestemc_to_litter_patch,
                iso_cnveg_cs.livestemc_patch, cnveg_cs.livestemc_patch)

        # crop_harvestc_to_cropprodc_patch
        for p in bounds_p
            mask_soilp[p] || continue
            iso_cnveg_cf.crop_harvestc_to_cropprodc_patch[p] =
                iso_cnveg_cf.leafc_to_biofuelc_patch[p] +
                iso_cnveg_cf.livestemc_to_biofuelc_patch[p] +
                iso_cnveg_cf.leafc_to_removedresiduec_patch[p] +
                iso_cnveg_cf.livestemc_to_removedresiduec_patch[p]
        end

        if use_grainproduct
            for k in repr_grain_min:repr_grain_max
                for p in bounds_p
                    mask_soilp[p] || continue
                    iso_cnveg_cf.crop_harvestc_to_cropprodc_patch[p] +=
                        iso_cnveg_cf.repr_grainc_to_food_patch[p, k]
                end
            end
        end

        for k in 1:nrepr
            for p in bounds_p
                mask_soilp[p] || continue
                iso_cnveg_cf.reproductive_mr_patch[p, k] =
                    iso_cnveg_cf.reproductive_xsmr_patch[p, k] +
                    iso_cnveg_cf.reproductive_curmr_patch[p, k]
            end
        end

        for k in repr_structure_min:repr_structure_max
            for p in bounds_p
                mask_soilp[p] || continue
                iso_cnveg_cf.crop_harvestc_to_cropprodc_patch[p] +=
                    iso_cnveg_cf.repr_structurec_to_cropprod_patch[p, k]
            end
        end
    end

    # --- Litter-to-column for isotopes ---
    cn_c_iso_litter_to_column!(iso_cnveg_cf, soilbiogeochem_state,
                               mask_soilp, bounds_p, nlevdecomp,
                               patch_column=patch_column, patch_itype=patch_itype,
                               patch_wtcol=patch_wtcol, lf_f=lf_f, fr_f=fr_f,
                               use_crop=use_crop, use_grainproduct=use_grainproduct,
                               npcropmin=npcropmin,
                               repr_grain_min=repr_grain_min, repr_grain_max=repr_grain_max,
                               repr_structure_min=repr_structure_min, repr_structure_max=repr_structure_max)

    # --- Column-level decomposition fluxes ---
    for c in bounds_c
        mask_soilc[c] || continue
        for j in 1:nlevdecomp
            for l in 1:ndecomp_cascade_transitions
                cdp = cascade_donor_pool[l]
                if soilbiogeochem_cs.decomp_cpools_vr_col[c, j, cdp] != 0.0
                    iso_soilbiogeochem_cf.decomp_cascade_hr_vr_col[c, j, l] =
                        soilbiogeochem_cf.decomp_cascade_hr_vr_col[c, j, l] *
                        (iso_soilbiogeochem_cs.decomp_cpools_vr_col[c, j, cdp] /
                         soilbiogeochem_cs.decomp_cpools_vr_col[c, j, cdp])
                else
                    iso_soilbiogeochem_cf.decomp_cascade_hr_vr_col[c, j, l] = 0.0
                end
            end
        end
    end

    for c in bounds_c
        mask_soilc[c] || continue
        for j in 1:nlevdecomp
            for l in 1:ndecomp_cascade_transitions
                cdp = cascade_donor_pool[l]
                if soilbiogeochem_cs.decomp_cpools_vr_col[c, j, cdp] != 0.0
                    iso_soilbiogeochem_cf.decomp_cascade_ctransfer_vr_col[c, j, l] =
                        soilbiogeochem_cf.decomp_cascade_ctransfer_vr_col[c, j, l] *
                        (iso_soilbiogeochem_cs.decomp_cpools_vr_col[c, j, cdp] /
                         soilbiogeochem_cs.decomp_cpools_vr_col[c, j, cdp])
                else
                    iso_soilbiogeochem_cf.decomp_cascade_ctransfer_vr_col[c, j, l] = 0.0
                end
            end
        end
    end

    return nothing
end

# ---------------------------------------------------------------------------
# c_iso_flux2! — Gap mortality isotopic fluxes
# Ported from: CIsoFlux2 in CNCIsoFluxMod.F90
# ---------------------------------------------------------------------------

"""
    c_iso_flux2!(soilbiogeochem_state, cnveg_cf, cnveg_cs,
                 iso_cnveg_cf, iso_cnveg_cs,
                 mask_soilp, bounds_p, isotope;
                 patch_column, patch_itype, patch_wtcol, lf_f, fr_f)

Set the carbon isotopic fluxes for gap mortality.

Ported from `CIsoFlux2` in `CNCIsoFluxMod.F90`.
"""
function c_iso_flux2!(soilbiogeochem_state::SoilBiogeochemStateData,
                      cnveg_cf::CNVegCarbonFluxData,
                      cnveg_cs::CNVegCarbonStateData,
                      iso_cnveg_cf::CNVegCarbonFluxData,
                      iso_cnveg_cs::CNVegCarbonStateData,
                      mask_soilp::BitVector,
                      bounds_p::UnitRange{Int},
                      nlevdecomp::Int,
                      isotope::String;
                      patch_column::Vector{Int}=Int[],
                      patch_itype::Vector{Int}=Int[],
                      patch_wtcol::Vector{Float64}=Float64[],
                      lf_f::Matrix{Float64}=Matrix{Float64}(undef, 0, 0),
                      fr_f::Matrix{Float64}=Matrix{Float64}(undef, 0, 0))

    calc1d! = (iso_f, tot_f, iso_s, tot_s) -> c_iso_flux_calc!(
        iso_f, tot_f, iso_s, tot_s, mask_soilp, bounds_p, 1.0, isotope)

    # patch-level gap mortality fluxes
    calc1d!(iso_cnveg_cf.m_leafc_to_litter_patch, cnveg_cf.m_leafc_to_litter_patch,
            iso_cnveg_cs.leafc_patch, cnveg_cs.leafc_patch)
    calc1d!(iso_cnveg_cf.m_leafc_storage_to_litter_patch, cnveg_cf.m_leafc_storage_to_litter_patch,
            iso_cnveg_cs.leafc_storage_patch, cnveg_cs.leafc_storage_patch)
    calc1d!(iso_cnveg_cf.m_leafc_xfer_to_litter_patch, cnveg_cf.m_leafc_xfer_to_litter_patch,
            iso_cnveg_cs.leafc_xfer_patch, cnveg_cs.leafc_xfer_patch)
    calc1d!(iso_cnveg_cf.m_frootc_to_litter_patch, cnveg_cf.m_frootc_to_litter_patch,
            iso_cnveg_cs.frootc_patch, cnveg_cs.frootc_patch)
    calc1d!(iso_cnveg_cf.m_frootc_storage_to_litter_patch, cnveg_cf.m_frootc_storage_to_litter_patch,
            iso_cnveg_cs.frootc_storage_patch, cnveg_cs.frootc_storage_patch)
    calc1d!(iso_cnveg_cf.m_frootc_xfer_to_litter_patch, cnveg_cf.m_frootc_xfer_to_litter_patch,
            iso_cnveg_cs.frootc_xfer_patch, cnveg_cs.frootc_xfer_patch)
    calc1d!(iso_cnveg_cf.m_livestemc_to_litter_patch, cnveg_cf.m_livestemc_to_litter_patch,
            iso_cnveg_cs.livestemc_patch, cnveg_cs.livestemc_patch)
    calc1d!(iso_cnveg_cf.m_livestemc_storage_to_litter_patch, cnveg_cf.m_livestemc_storage_to_litter_patch,
            iso_cnveg_cs.livestemc_storage_patch, cnveg_cs.livestemc_storage_patch)
    calc1d!(iso_cnveg_cf.m_livestemc_xfer_to_litter_patch, cnveg_cf.m_livestemc_xfer_to_litter_patch,
            iso_cnveg_cs.livestemc_xfer_patch, cnveg_cs.livestemc_xfer_patch)
    calc1d!(iso_cnveg_cf.m_deadstemc_to_litter_patch, cnveg_cf.m_deadstemc_to_litter_patch,
            iso_cnveg_cs.deadstemc_patch, cnveg_cs.deadstemc_patch)
    calc1d!(iso_cnveg_cf.m_deadstemc_storage_to_litter_patch, cnveg_cf.m_deadstemc_storage_to_litter_patch,
            iso_cnveg_cs.deadstemc_storage_patch, cnveg_cs.deadstemc_storage_patch)
    calc1d!(iso_cnveg_cf.m_deadstemc_xfer_to_litter_patch, cnveg_cf.m_deadstemc_xfer_to_litter_patch,
            iso_cnveg_cs.deadstemc_xfer_patch, cnveg_cs.deadstemc_xfer_patch)
    calc1d!(iso_cnveg_cf.m_livecrootc_to_litter_patch, cnveg_cf.m_livecrootc_to_litter_patch,
            iso_cnveg_cs.livecrootc_patch, cnveg_cs.livecrootc_patch)
    calc1d!(iso_cnveg_cf.m_livecrootc_storage_to_litter_patch, cnveg_cf.m_livecrootc_storage_to_litter_patch,
            iso_cnveg_cs.livecrootc_storage_patch, cnveg_cs.livecrootc_storage_patch)
    calc1d!(iso_cnveg_cf.m_livecrootc_xfer_to_litter_patch, cnveg_cf.m_livecrootc_xfer_to_litter_patch,
            iso_cnveg_cs.livecrootc_xfer_patch, cnveg_cs.livecrootc_xfer_patch)
    calc1d!(iso_cnveg_cf.m_deadcrootc_to_litter_patch, cnveg_cf.m_deadcrootc_to_litter_patch,
            iso_cnveg_cs.deadcrootc_patch, cnveg_cs.deadcrootc_patch)
    calc1d!(iso_cnveg_cf.m_deadcrootc_storage_to_litter_patch, cnveg_cf.m_deadcrootc_storage_to_litter_patch,
            iso_cnveg_cs.deadcrootc_storage_patch, cnveg_cs.deadcrootc_storage_patch)
    calc1d!(iso_cnveg_cf.m_deadcrootc_xfer_to_litter_patch, cnveg_cf.m_deadcrootc_xfer_to_litter_patch,
            iso_cnveg_cs.deadcrootc_xfer_patch, cnveg_cs.deadcrootc_xfer_patch)
    calc1d!(iso_cnveg_cf.m_gresp_storage_to_litter_patch, cnveg_cf.m_gresp_storage_to_litter_patch,
            iso_cnveg_cs.gresp_storage_patch, cnveg_cs.gresp_storage_patch)
    calc1d!(iso_cnveg_cf.m_gresp_xfer_to_litter_patch, cnveg_cf.m_gresp_xfer_to_litter_patch,
            iso_cnveg_cs.gresp_xfer_patch, cnveg_cs.gresp_xfer_patch)

    # shift patch-level gap mortality fluxes to column
    cn_c_iso_gap_pft_to_column!(iso_cnveg_cf, soilbiogeochem_state,
                                mask_soilp, bounds_p, nlevdecomp,
                                patch_column=patch_column, patch_itype=patch_itype,
                                patch_wtcol=patch_wtcol, lf_f=lf_f, fr_f=fr_f)

    return nothing
end

# ---------------------------------------------------------------------------
# c_iso_flux2h! — Harvest mortality isotopic fluxes
# Ported from: CIsoFlux2h in CNCIsoFluxMod.F90
# ---------------------------------------------------------------------------

"""
    c_iso_flux2h!(soilbiogeochem_state, cnveg_cf, cnveg_cs,
                  iso_cnveg_cf, iso_cnveg_cs,
                  mask_soilp, bounds_p, isotope; ...)

Set the carbon isotopic fluxes for harvest mortality.

Ported from `CIsoFlux2h` in `CNCIsoFluxMod.F90`.
"""
function c_iso_flux2h!(soilbiogeochem_state::SoilBiogeochemStateData,
                       cnveg_cf::CNVegCarbonFluxData,
                       cnveg_cs::CNVegCarbonStateData,
                       iso_cnveg_cf::CNVegCarbonFluxData,
                       iso_cnveg_cs::CNVegCarbonStateData,
                       mask_soilp::BitVector,
                       bounds_p::UnitRange{Int},
                       nlevdecomp::Int,
                       isotope::String;
                       patch_column::Vector{Int}=Int[],
                       patch_itype::Vector{Int}=Int[],
                       patch_wtcol::Vector{Float64}=Float64[],
                       lf_f::Matrix{Float64}=Matrix{Float64}(undef, 0, 0),
                       fr_f::Matrix{Float64}=Matrix{Float64}(undef, 0, 0))

    calc1d! = (iso_f, tot_f, iso_s, tot_s) -> c_iso_flux_calc!(
        iso_f, tot_f, iso_s, tot_s, mask_soilp, bounds_p, 1.0, isotope)

    calc1d!(iso_cnveg_cf.hrv_leafc_to_litter_patch, cnveg_cf.hrv_leafc_to_litter_patch,
            iso_cnveg_cs.leafc_patch, cnveg_cs.leafc_patch)
    calc1d!(iso_cnveg_cf.hrv_leafc_storage_to_litter_patch, cnveg_cf.hrv_leafc_storage_to_litter_patch,
            iso_cnveg_cs.leafc_storage_patch, cnveg_cs.leafc_storage_patch)
    calc1d!(iso_cnveg_cf.hrv_leafc_xfer_to_litter_patch, cnveg_cf.hrv_leafc_xfer_to_litter_patch,
            iso_cnveg_cs.leafc_xfer_patch, cnveg_cs.leafc_xfer_patch)
    calc1d!(iso_cnveg_cf.hrv_frootc_to_litter_patch, cnveg_cf.hrv_frootc_to_litter_patch,
            iso_cnveg_cs.frootc_patch, cnveg_cs.frootc_patch)
    calc1d!(iso_cnveg_cf.hrv_frootc_storage_to_litter_patch, cnveg_cf.hrv_frootc_storage_to_litter_patch,
            iso_cnveg_cs.frootc_storage_patch, cnveg_cs.frootc_storage_patch)
    calc1d!(iso_cnveg_cf.hrv_frootc_xfer_to_litter_patch, cnveg_cf.hrv_frootc_xfer_to_litter_patch,
            iso_cnveg_cs.frootc_xfer_patch, cnveg_cs.frootc_xfer_patch)
    calc1d!(iso_cnveg_cf.hrv_livestemc_to_litter_patch, cnveg_cf.hrv_livestemc_to_litter_patch,
            iso_cnveg_cs.livestemc_patch, cnveg_cs.livestemc_patch)
    calc1d!(iso_cnveg_cf.hrv_livestemc_storage_to_litter_patch, cnveg_cf.hrv_livestemc_storage_to_litter_patch,
            iso_cnveg_cs.livestemc_storage_patch, cnveg_cs.livestemc_storage_patch)
    calc1d!(iso_cnveg_cf.hrv_livestemc_xfer_to_litter_patch, cnveg_cf.hrv_livestemc_xfer_to_litter_patch,
            iso_cnveg_cs.livestemc_xfer_patch, cnveg_cs.livestemc_xfer_patch)
    calc1d!(iso_cnveg_cf.wood_harvestc_patch, cnveg_cf.wood_harvestc_patch,
            iso_cnveg_cs.deadstemc_patch, cnveg_cs.deadstemc_patch)
    calc1d!(iso_cnveg_cf.hrv_deadstemc_storage_to_litter_patch, cnveg_cf.hrv_deadstemc_storage_to_litter_patch,
            iso_cnveg_cs.deadstemc_storage_patch, cnveg_cs.deadstemc_storage_patch)
    calc1d!(iso_cnveg_cf.hrv_deadstemc_xfer_to_litter_patch, cnveg_cf.hrv_deadstemc_xfer_to_litter_patch,
            iso_cnveg_cs.deadstemc_xfer_patch, cnveg_cs.deadstemc_xfer_patch)
    calc1d!(iso_cnveg_cf.hrv_livecrootc_to_litter_patch, cnveg_cf.hrv_livecrootc_to_litter_patch,
            iso_cnveg_cs.livecrootc_patch, cnveg_cs.livecrootc_patch)
    calc1d!(iso_cnveg_cf.hrv_livecrootc_storage_to_litter_patch, cnveg_cf.hrv_livecrootc_storage_to_litter_patch,
            iso_cnveg_cs.livecrootc_storage_patch, cnveg_cs.livecrootc_storage_patch)
    calc1d!(iso_cnveg_cf.hrv_livecrootc_xfer_to_litter_patch, cnveg_cf.hrv_livecrootc_xfer_to_litter_patch,
            iso_cnveg_cs.livecrootc_xfer_patch, cnveg_cs.livecrootc_xfer_patch)
    calc1d!(iso_cnveg_cf.hrv_deadcrootc_to_litter_patch, cnveg_cf.hrv_deadcrootc_to_litter_patch,
            iso_cnveg_cs.deadcrootc_patch, cnveg_cs.deadcrootc_patch)
    calc1d!(iso_cnveg_cf.hrv_deadcrootc_storage_to_litter_patch, cnveg_cf.hrv_deadcrootc_storage_to_litter_patch,
            iso_cnveg_cs.deadcrootc_storage_patch, cnveg_cs.deadcrootc_storage_patch)
    calc1d!(iso_cnveg_cf.hrv_deadcrootc_xfer_to_litter_patch, cnveg_cf.hrv_deadcrootc_xfer_to_litter_patch,
            iso_cnveg_cs.deadcrootc_xfer_patch, cnveg_cs.deadcrootc_xfer_patch)
    calc1d!(iso_cnveg_cf.hrv_gresp_storage_to_litter_patch, cnveg_cf.hrv_gresp_storage_to_litter_patch,
            iso_cnveg_cs.gresp_storage_patch, cnveg_cs.gresp_storage_patch)
    calc1d!(iso_cnveg_cf.hrv_gresp_xfer_to_litter_patch, cnveg_cf.hrv_gresp_xfer_to_litter_patch,
            iso_cnveg_cs.gresp_xfer_patch, cnveg_cs.gresp_xfer_patch)
    calc1d!(iso_cnveg_cf.hrv_xsmrpool_to_atm_patch, cnveg_cf.hrv_xsmrpool_to_atm_patch,
            iso_cnveg_cs.totvegc_patch, cnveg_cs.totvegc_patch)

    # shift patch-level harvest mortality fluxes to column
    cn_c_iso_harvest_pft_to_column!(iso_cnveg_cf, soilbiogeochem_state,
                                    mask_soilp, bounds_p, nlevdecomp,
                                    patch_column=patch_column, patch_itype=patch_itype,
                                    patch_wtcol=patch_wtcol, lf_f=lf_f, fr_f=fr_f)

    return nothing
end

# ---------------------------------------------------------------------------
# c_iso_flux2g! — Gross unrepresented landcover change isotopic fluxes
# Ported from: CIsoFlux2g in CNCIsoFluxMod.F90
# ---------------------------------------------------------------------------

"""
    c_iso_flux2g!(soilbiogeochem_state, cnveg_cf, cnveg_cs,
                  iso_cnveg_cf, iso_cnveg_cs,
                  mask_soilp, bounds_p, isotope; ...)

Set the carbon isotopic fluxes for gross unrepresented landcover change mortality.

Ported from `CIsoFlux2g` in `CNCIsoFluxMod.F90`.
"""
function c_iso_flux2g!(soilbiogeochem_state::SoilBiogeochemStateData,
                       cnveg_cf::CNVegCarbonFluxData,
                       cnveg_cs::CNVegCarbonStateData,
                       iso_cnveg_cf::CNVegCarbonFluxData,
                       iso_cnveg_cs::CNVegCarbonStateData,
                       mask_soilp::BitVector,
                       bounds_p::UnitRange{Int},
                       nlevdecomp::Int,
                       isotope::String;
                       patch_column::Vector{Int}=Int[],
                       patch_itype::Vector{Int}=Int[],
                       patch_wtcol::Vector{Float64}=Float64[],
                       lf_f::Matrix{Float64}=Matrix{Float64}(undef, 0, 0),
                       fr_f::Matrix{Float64}=Matrix{Float64}(undef, 0, 0))

    calc1d! = (iso_f, tot_f, iso_s, tot_s) -> c_iso_flux_calc!(
        iso_f, tot_f, iso_s, tot_s, mask_soilp, bounds_p, 1.0, isotope)

    calc1d!(iso_cnveg_cf.gru_leafc_to_litter_patch, cnveg_cf.gru_leafc_to_litter_patch,
            iso_cnveg_cs.leafc_patch, cnveg_cs.leafc_patch)
    calc1d!(iso_cnveg_cf.gru_leafc_storage_to_atm_patch, cnveg_cf.gru_leafc_storage_to_atm_patch,
            iso_cnveg_cs.leafc_storage_patch, cnveg_cs.leafc_storage_patch)
    calc1d!(iso_cnveg_cf.gru_leafc_xfer_to_atm_patch, cnveg_cf.gru_leafc_xfer_to_atm_patch,
            iso_cnveg_cs.leafc_xfer_patch, cnveg_cs.leafc_xfer_patch)
    calc1d!(iso_cnveg_cf.gru_frootc_to_litter_patch, cnveg_cf.gru_frootc_to_litter_patch,
            iso_cnveg_cs.frootc_patch, cnveg_cs.frootc_patch)
    calc1d!(iso_cnveg_cf.gru_frootc_storage_to_atm_patch, cnveg_cf.gru_frootc_storage_to_atm_patch,
            iso_cnveg_cs.frootc_storage_patch, cnveg_cs.frootc_storage_patch)
    calc1d!(iso_cnveg_cf.gru_frootc_xfer_to_atm_patch, cnveg_cf.gru_frootc_xfer_to_atm_patch,
            iso_cnveg_cs.frootc_xfer_patch, cnveg_cs.frootc_xfer_patch)
    calc1d!(iso_cnveg_cf.gru_livestemc_to_atm_patch, cnveg_cf.gru_livestemc_to_atm_patch,
            iso_cnveg_cs.livestemc_patch, cnveg_cs.livestemc_patch)
    calc1d!(iso_cnveg_cf.gru_livestemc_storage_to_atm_patch, cnveg_cf.gru_livestemc_storage_to_atm_patch,
            iso_cnveg_cs.livestemc_storage_patch, cnveg_cs.livestemc_storage_patch)
    calc1d!(iso_cnveg_cf.gru_livestemc_xfer_to_atm_patch, cnveg_cf.gru_livestemc_xfer_to_atm_patch,
            iso_cnveg_cs.livestemc_xfer_patch, cnveg_cs.livestemc_xfer_patch)
    calc1d!(iso_cnveg_cf.gru_deadstemc_to_atm_patch, cnveg_cf.gru_deadstemc_to_atm_patch,
            iso_cnveg_cs.deadstemc_patch, cnveg_cs.deadstemc_patch)
    calc1d!(iso_cnveg_cf.gru_wood_productc_gain_patch, cnveg_cf.gru_wood_productc_gain_patch,
            iso_cnveg_cs.deadstemc_patch, cnveg_cs.deadstemc_patch)
    calc1d!(iso_cnveg_cf.gru_deadstemc_storage_to_atm_patch, cnveg_cf.gru_deadstemc_storage_to_atm_patch,
            iso_cnveg_cs.deadstemc_storage_patch, cnveg_cs.deadstemc_storage_patch)
    calc1d!(iso_cnveg_cf.gru_deadstemc_xfer_to_atm_patch, cnveg_cf.gru_deadstemc_xfer_to_atm_patch,
            iso_cnveg_cs.deadstemc_xfer_patch, cnveg_cs.deadstemc_xfer_patch)
    calc1d!(iso_cnveg_cf.gru_livecrootc_to_litter_patch, cnveg_cf.gru_livecrootc_to_litter_patch,
            iso_cnveg_cs.livecrootc_patch, cnveg_cs.livecrootc_patch)
    calc1d!(iso_cnveg_cf.gru_livecrootc_storage_to_atm_patch, cnveg_cf.gru_livecrootc_storage_to_atm_patch,
            iso_cnveg_cs.livecrootc_storage_patch, cnveg_cs.livecrootc_storage_patch)
    calc1d!(iso_cnveg_cf.gru_livecrootc_xfer_to_atm_patch, cnveg_cf.gru_livecrootc_xfer_to_atm_patch,
            iso_cnveg_cs.livecrootc_xfer_patch, cnveg_cs.livecrootc_xfer_patch)
    calc1d!(iso_cnveg_cf.gru_deadcrootc_to_litter_patch, cnveg_cf.gru_deadcrootc_to_litter_patch,
            iso_cnveg_cs.deadcrootc_patch, cnveg_cs.deadcrootc_patch)
    calc1d!(iso_cnveg_cf.gru_deadcrootc_storage_to_atm_patch, cnveg_cf.gru_deadcrootc_storage_to_atm_patch,
            iso_cnveg_cs.deadcrootc_storage_patch, cnveg_cs.deadcrootc_storage_patch)
    calc1d!(iso_cnveg_cf.gru_deadcrootc_xfer_to_atm_patch, cnveg_cf.gru_deadcrootc_xfer_to_atm_patch,
            iso_cnveg_cs.deadcrootc_xfer_patch, cnveg_cs.deadcrootc_xfer_patch)
    calc1d!(iso_cnveg_cf.gru_gresp_storage_to_atm_patch, cnveg_cf.gru_gresp_storage_to_atm_patch,
            iso_cnveg_cs.gresp_storage_patch, cnveg_cs.gresp_storage_patch)
    calc1d!(iso_cnveg_cf.gru_gresp_xfer_to_atm_patch, cnveg_cf.gru_gresp_xfer_to_atm_patch,
            iso_cnveg_cs.gresp_xfer_patch, cnveg_cs.gresp_xfer_patch)
    calc1d!(iso_cnveg_cf.gru_xsmrpool_to_atm_patch, cnveg_cf.gru_xsmrpool_to_atm_patch,
            iso_cnveg_cs.totvegc_patch, cnveg_cs.totvegc_patch)

    # shift patch-level gross unrep fluxes to column
    cn_c_iso_gross_unrep_pft_to_column!(iso_cnveg_cf, soilbiogeochem_state,
                                        mask_soilp, bounds_p, nlevdecomp,
                                        patch_column=patch_column, patch_itype=patch_itype,
                                        patch_wtcol=patch_wtcol, lf_f=lf_f, fr_f=fr_f)

    return nothing
end

# ---------------------------------------------------------------------------
# c_iso_flux3! — Fire mortality isotopic fluxes
# Ported from: CIsoFlux3 in CNCIsoFluxMod.F90
# ---------------------------------------------------------------------------

"""
    c_iso_flux3!(soilbiogeochem_state, soilbiogeochem_cs,
                 cnveg_cf, cnveg_cs,
                 iso_cnveg_cf, iso_cnveg_cs,
                 iso_soilbiogeochem_cs,
                 mask_soilp, bounds_p, nlevdecomp, ndecomp_pools,
                 i_met_lit, i_litr_max, isotope; ...)

Set the carbon isotopic fluxes for fire mortality.

Ported from `CIsoFlux3` in `CNCIsoFluxMod.F90`.
"""
function c_iso_flux3!(soilbiogeochem_state::SoilBiogeochemStateData,
                      soilbiogeochem_cs::SoilBiogeochemCarbonStateData,
                      cnveg_cf::CNVegCarbonFluxData,
                      cnveg_cs::CNVegCarbonStateData,
                      iso_cnveg_cf::CNVegCarbonFluxData,
                      iso_cnveg_cs::CNVegCarbonStateData,
                      iso_soilbiogeochem_cs::SoilBiogeochemCarbonStateData,
                      mask_soilp::BitVector,
                      bounds_p::UnitRange{Int},
                      nlevdecomp::Int,
                      ndecomp_pools::Int,
                      i_met_lit::Int,
                      i_litr_max::Int,
                      isotope::String;
                      patch_column::Vector{Int}=Int[],
                      patch_itype::Vector{Int}=Int[],
                      patch_wtcol::Vector{Float64}=Float64[],
                      lf_f::Matrix{Float64}=Matrix{Float64}(undef, 0, 0),
                      fr_f::Matrix{Float64}=Matrix{Float64}(undef, 0, 0),
                      leaf_prof::Matrix{Float64}=Matrix{Float64}(undef, 0, 0),
                      froot_prof::Matrix{Float64}=Matrix{Float64}(undef, 0, 0),
                      stem_prof::Matrix{Float64}=Matrix{Float64}(undef, 0, 0),
                      croot_prof::Matrix{Float64}=Matrix{Float64}(undef, 0, 0))

    calc1d! = (iso_f, tot_f, iso_s, tot_s) -> c_iso_flux_calc!(
        iso_f, tot_f, iso_s, tot_s, mask_soilp, bounds_p, 1.0, isotope)

    # patch-level fire mortality fluxes
    calc1d!(iso_cnveg_cf.m_leafc_to_fire_patch, cnveg_cf.m_leafc_to_fire_patch,
            iso_cnveg_cs.leafc_patch, cnveg_cs.leafc_patch)
    calc1d!(iso_cnveg_cf.m_leafc_storage_to_fire_patch, cnveg_cf.m_leafc_storage_to_fire_patch,
            iso_cnveg_cs.leafc_storage_patch, cnveg_cs.leafc_storage_patch)
    calc1d!(iso_cnveg_cf.m_leafc_xfer_to_fire_patch, cnveg_cf.m_leafc_xfer_to_fire_patch,
            iso_cnveg_cs.leafc_xfer_patch, cnveg_cs.leafc_xfer_patch)
    calc1d!(iso_cnveg_cf.m_frootc_to_fire_patch, cnveg_cf.m_frootc_to_fire_patch,
            iso_cnveg_cs.frootc_patch, cnveg_cs.frootc_patch)
    calc1d!(iso_cnveg_cf.m_frootc_storage_to_fire_patch, cnveg_cf.m_frootc_storage_to_fire_patch,
            iso_cnveg_cs.frootc_storage_patch, cnveg_cs.frootc_storage_patch)
    calc1d!(iso_cnveg_cf.m_frootc_xfer_to_fire_patch, cnveg_cf.m_frootc_xfer_to_fire_patch,
            iso_cnveg_cs.frootc_xfer_patch, cnveg_cs.frootc_xfer_patch)
    calc1d!(iso_cnveg_cf.m_livestemc_to_fire_patch, cnveg_cf.m_livestemc_to_fire_patch,
            iso_cnveg_cs.livestemc_patch, cnveg_cs.livestemc_patch)
    calc1d!(iso_cnveg_cf.m_livestemc_storage_to_fire_patch, cnveg_cf.m_livestemc_storage_to_fire_patch,
            iso_cnveg_cs.livestemc_storage_patch, cnveg_cs.livestemc_storage_patch)
    calc1d!(iso_cnveg_cf.m_livestemc_xfer_to_fire_patch, cnveg_cf.m_livestemc_xfer_to_fire_patch,
            iso_cnveg_cs.livestemc_xfer_patch, cnveg_cs.livestemc_xfer_patch)
    calc1d!(iso_cnveg_cf.m_deadstemc_to_fire_patch, cnveg_cf.m_deadstemc_to_fire_patch,
            iso_cnveg_cs.deadstemc_patch, cnveg_cs.deadstemc_patch)
    calc1d!(iso_cnveg_cf.m_deadstemc_to_litter_fire_patch, cnveg_cf.m_deadstemc_to_litter_fire_patch,
            iso_cnveg_cs.deadstemc_patch, cnveg_cs.deadstemc_patch)
    calc1d!(iso_cnveg_cf.m_deadstemc_storage_to_fire_patch, cnveg_cf.m_deadstemc_storage_to_fire_patch,
            iso_cnveg_cs.deadstemc_storage_patch, cnveg_cs.deadstemc_storage_patch)
    calc1d!(iso_cnveg_cf.m_deadstemc_xfer_to_fire_patch, cnveg_cf.m_deadstemc_xfer_to_fire_patch,
            iso_cnveg_cs.deadstemc_xfer_patch, cnveg_cs.deadstemc_xfer_patch)
    calc1d!(iso_cnveg_cf.m_livecrootc_to_fire_patch, cnveg_cf.m_livecrootc_to_fire_patch,
            iso_cnveg_cs.livecrootc_patch, cnveg_cs.livecrootc_patch)
    calc1d!(iso_cnveg_cf.m_livecrootc_storage_to_fire_patch, cnveg_cf.m_livecrootc_storage_to_fire_patch,
            iso_cnveg_cs.livecrootc_storage_patch, cnveg_cs.livecrootc_storage_patch)
    calc1d!(iso_cnveg_cf.m_livecrootc_xfer_to_fire_patch, cnveg_cf.m_livecrootc_xfer_to_fire_patch,
            iso_cnveg_cs.livecrootc_xfer_patch, cnveg_cs.livecrootc_xfer_patch)
    calc1d!(iso_cnveg_cf.m_deadcrootc_to_fire_patch, cnveg_cf.m_deadcrootc_to_fire_patch,
            iso_cnveg_cs.deadcrootc_patch, cnveg_cs.deadcrootc_patch)
    calc1d!(iso_cnveg_cf.m_deadcrootc_to_litter_fire_patch, cnveg_cf.m_deadcrootc_to_litter_fire_patch,
            iso_cnveg_cs.deadcrootc_patch, cnveg_cs.deadcrootc_patch)
    calc1d!(iso_cnveg_cf.m_deadcrootc_storage_to_fire_patch, cnveg_cf.m_deadcrootc_storage_to_fire_patch,
            iso_cnveg_cs.deadcrootc_storage_patch, cnveg_cs.deadcrootc_storage_patch)
    calc1d!(iso_cnveg_cf.m_deadcrootc_xfer_to_fire_patch, cnveg_cf.m_deadcrootc_xfer_to_fire_patch,
            iso_cnveg_cs.deadcrootc_xfer_patch, cnveg_cs.deadcrootc_xfer_patch)
    calc1d!(iso_cnveg_cf.m_gresp_storage_to_fire_patch, cnveg_cf.m_gresp_storage_to_fire_patch,
            iso_cnveg_cs.gresp_storage_patch, cnveg_cs.gresp_storage_patch)
    calc1d!(iso_cnveg_cf.m_gresp_xfer_to_fire_patch, cnveg_cf.m_gresp_xfer_to_fire_patch,
            iso_cnveg_cs.gresp_xfer_patch, cnveg_cs.gresp_xfer_patch)

    # fire litter fluxes
    calc1d!(iso_cnveg_cf.m_leafc_to_litter_fire_patch, cnveg_cf.m_leafc_to_litter_fire_patch,
            iso_cnveg_cs.leafc_patch, cnveg_cs.leafc_patch)
    calc1d!(iso_cnveg_cf.m_leafc_storage_to_litter_fire_patch, cnveg_cf.m_leafc_storage_to_litter_fire_patch,
            iso_cnveg_cs.leafc_storage_patch, cnveg_cs.leafc_storage_patch)
    calc1d!(iso_cnveg_cf.m_leafc_xfer_to_litter_fire_patch, cnveg_cf.m_leafc_xfer_to_litter_fire_patch,
            iso_cnveg_cs.leafc_xfer_patch, cnveg_cs.leafc_xfer_patch)
    calc1d!(iso_cnveg_cf.m_livestemc_to_litter_fire_patch, cnveg_cf.m_livestemc_to_litter_fire_patch,
            iso_cnveg_cs.livestemc_patch, cnveg_cs.livestemc_patch)
    calc1d!(iso_cnveg_cf.m_livestemc_storage_to_litter_fire_patch, cnveg_cf.m_livestemc_storage_to_litter_fire_patch,
            iso_cnveg_cs.livestemc_storage_patch, cnveg_cs.livestemc_storage_patch)
    calc1d!(iso_cnveg_cf.m_livestemc_xfer_to_litter_fire_patch, cnveg_cf.m_livestemc_xfer_to_litter_fire_patch,
            iso_cnveg_cs.livestemc_xfer_patch, cnveg_cs.livestemc_xfer_patch)
    calc1d!(iso_cnveg_cf.m_livestemc_to_deadstemc_fire_patch, cnveg_cf.m_livestemc_to_deadstemc_fire_patch,
            iso_cnveg_cs.livestemc_patch, cnveg_cs.livestemc_patch)
    calc1d!(iso_cnveg_cf.m_deadstemc_storage_to_litter_fire_patch, cnveg_cf.m_deadstemc_storage_to_litter_fire_patch,
            iso_cnveg_cs.deadstemc_storage_patch, cnveg_cs.deadstemc_storage_patch)
    calc1d!(iso_cnveg_cf.m_deadstemc_xfer_to_litter_fire_patch, cnveg_cf.m_deadstemc_xfer_to_litter_fire_patch,
            iso_cnveg_cs.deadstemc_xfer_patch, cnveg_cs.deadstemc_xfer_patch)
    calc1d!(iso_cnveg_cf.m_frootc_to_litter_fire_patch, cnveg_cf.m_frootc_to_litter_fire_patch,
            iso_cnveg_cs.frootc_patch, cnveg_cs.frootc_patch)
    calc1d!(iso_cnveg_cf.m_frootc_storage_to_litter_fire_patch, cnveg_cf.m_frootc_storage_to_litter_fire_patch,
            iso_cnveg_cs.frootc_storage_patch, cnveg_cs.frootc_storage_patch)
    calc1d!(iso_cnveg_cf.m_frootc_xfer_to_litter_fire_patch, cnveg_cf.m_frootc_xfer_to_litter_fire_patch,
            iso_cnveg_cs.frootc_xfer_patch, cnveg_cs.frootc_xfer_patch)
    calc1d!(iso_cnveg_cf.m_livecrootc_to_litter_fire_patch, cnveg_cf.m_livecrootc_to_litter_fire_patch,
            iso_cnveg_cs.livecrootc_patch, cnveg_cs.livecrootc_patch)
    calc1d!(iso_cnveg_cf.m_livecrootc_storage_to_litter_fire_patch, cnveg_cf.m_livecrootc_storage_to_litter_fire_patch,
            iso_cnveg_cs.livecrootc_storage_patch, cnveg_cs.livecrootc_storage_patch)
    calc1d!(iso_cnveg_cf.m_livecrootc_xfer_to_litter_fire_patch, cnveg_cf.m_livecrootc_xfer_to_litter_fire_patch,
            iso_cnveg_cs.livecrootc_xfer_patch, cnveg_cs.livecrootc_xfer_patch)
    calc1d!(iso_cnveg_cf.m_livecrootc_to_deadcrootc_fire_patch, cnveg_cf.m_livecrootc_to_deadcrootc_fire_patch,
            iso_cnveg_cs.livecrootc_patch, cnveg_cs.livecrootc_patch)
    calc1d!(iso_cnveg_cf.m_deadcrootc_storage_to_litter_fire_patch, cnveg_cf.m_deadcrootc_storage_to_litter_fire_patch,
            iso_cnveg_cs.deadcrootc_storage_patch, cnveg_cs.deadcrootc_storage_patch)
    calc1d!(iso_cnveg_cf.m_deadcrootc_xfer_to_litter_fire_patch, cnveg_cf.m_deadcrootc_xfer_to_litter_fire_patch,
            iso_cnveg_cs.deadcrootc_xfer_patch, cnveg_cs.deadcrootc_xfer_patch)
    calc1d!(iso_cnveg_cf.m_gresp_storage_to_litter_fire_patch, cnveg_cf.m_gresp_storage_to_litter_fire_patch,
            iso_cnveg_cs.gresp_storage_patch, cnveg_cs.gresp_storage_patch)
    calc1d!(iso_cnveg_cf.m_gresp_xfer_to_litter_fire_patch, cnveg_cf.m_gresp_xfer_to_litter_fire_patch,
            iso_cnveg_cs.gresp_xfer_patch, cnveg_cs.gresp_xfer_patch)

    # column-level fire mortality fluxes
    ivt = patch_itype
    wtcol = patch_wtcol

    for p in bounds_p
        mask_soilp[p] || continue
        cc = patch_column[p]
        for j in 1:nlevdecomp
            iso_cnveg_cf.fire_mortality_c_to_cwdc_col[cc, j] +=
                (iso_cnveg_cf.m_deadstemc_to_litter_fire_patch[p] +
                 iso_cnveg_cf.m_livestemc_to_litter_fire_patch[p]) *
                wtcol[p] * stem_prof[p, j]
            iso_cnveg_cf.fire_mortality_c_to_cwdc_col[cc, j] +=
                (iso_cnveg_cf.m_deadcrootc_to_litter_fire_patch[p] +
                 iso_cnveg_cf.m_livecrootc_to_litter_fire_patch[p]) *
                wtcol[p] * croot_prof[p, j]
        end

        for j in 1:nlevdecomp
            for l in 1:ndecomp_pools
                if soilbiogeochem_cs.decomp_cpools_vr_col[cc, j, l] != 0.0
                    iso_cnveg_cf.m_decomp_cpools_to_fire_vr_col[cc, j, l] =
                        cnveg_cf.m_decomp_cpools_to_fire_vr_col[cc, j, l] *
                        (iso_soilbiogeochem_cs.decomp_cpools_vr_col[cc, j, l] /
                         soilbiogeochem_cs.decomp_cpools_vr_col[cc, j, l])
                else
                    iso_cnveg_cf.m_decomp_cpools_to_fire_vr_col[cc, j, l] = 0.0
                end
            end
        end
    end

    for p in bounds_p
        mask_soilp[p] || continue
        cc = patch_column[p]
        for j in 1:nlevdecomp
            iso_cnveg_cf.m_c_to_litr_fire_col[cc, j, i_met_lit] +=
                ((iso_cnveg_cf.m_leafc_to_litter_fire_patch[p] * lf_f[ivt[p], i_met_lit] +
                  iso_cnveg_cf.m_leafc_storage_to_litter_fire_patch[p] +
                  iso_cnveg_cf.m_leafc_xfer_to_litter_fire_patch[p] +
                  iso_cnveg_cf.m_gresp_storage_to_litter_fire_patch[p] +
                  iso_cnveg_cf.m_gresp_xfer_to_litter_fire_patch[p]) * leaf_prof[p, j] +
                 (iso_cnveg_cf.m_frootc_to_litter_fire_patch[p] * fr_f[ivt[p], i_met_lit] +
                  iso_cnveg_cf.m_frootc_storage_to_litter_fire_patch[p] +
                  iso_cnveg_cf.m_frootc_xfer_to_litter_fire_patch[p]) * froot_prof[p, j] +
                 (iso_cnveg_cf.m_livestemc_storage_to_litter_fire_patch[p] +
                  iso_cnveg_cf.m_livestemc_xfer_to_litter_fire_patch[p] +
                  iso_cnveg_cf.m_deadstemc_storage_to_litter_fire_patch[p] +
                  iso_cnveg_cf.m_deadstemc_xfer_to_litter_fire_patch[p]) * stem_prof[p, j] +
                 (iso_cnveg_cf.m_livecrootc_storage_to_litter_fire_patch[p] +
                  iso_cnveg_cf.m_livecrootc_xfer_to_litter_fire_patch[p] +
                  iso_cnveg_cf.m_deadcrootc_storage_to_litter_fire_patch[p] +
                  iso_cnveg_cf.m_deadcrootc_xfer_to_litter_fire_patch[p]) * croot_prof[p, j]) * wtcol[p]

            for i in (i_met_lit+1):i_litr_max
                iso_cnveg_cf.m_c_to_litr_fire_col[cc, j, i] +=
                    (iso_cnveg_cf.m_leafc_to_litter_fire_patch[p] * lf_f[ivt[p], i] * leaf_prof[p, j] +
                     iso_cnveg_cf.m_frootc_to_litter_fire_patch[p] * fr_f[ivt[p], i] * froot_prof[p, j]) * wtcol[p]
            end
        end
    end

    return nothing
end

# ---------------------------------------------------------------------------
# cn_c_iso_litter_to_column! — Phenology litterfall patch-to-column
# Ported from: CNCIsoLitterToColumn in CNCIsoFluxMod.F90
# ---------------------------------------------------------------------------

function cn_c_iso_litter_to_column!(iso_cnveg_cf::CNVegCarbonFluxData,
                                    soilbiogeochem_state::SoilBiogeochemStateData,
                                    mask_soilp::BitVector,
                                    bounds_p::UnitRange{Int},
                                    nlevdecomp::Int;
                                    patch_column::Vector{Int}=Int[],
                                    patch_itype::Vector{Int}=Int[],
                                    patch_wtcol::Vector{Float64}=Float64[],
                                    lf_f::Matrix{Float64}=Matrix{Float64}(undef, 0, 0),
                                    fr_f::Matrix{Float64}=Matrix{Float64}(undef, 0, 0),
                                    use_crop::Bool=false,
                                    use_grainproduct::Bool=false,
                                    npcropmin::Int=NPCROPMIN,
                                    i_litr_min::Int=varpar.i_litr_min,
                                    i_litr_max::Int=varpar.i_litr_max,
                                    repr_grain_min::Int=1,
                                    repr_grain_max::Int=1,
                                    repr_structure_min::Int=1,
                                    repr_structure_max::Int=0)

    ivt = patch_itype
    wtcol = patch_wtcol
    leaf_prof = soilbiogeochem_state.leaf_prof_patch
    froot_prof = soilbiogeochem_state.froot_prof_patch

    for j in 1:nlevdecomp
        for p in bounds_p
            mask_soilp[p] || continue
            c = patch_column[p]

            for i in i_litr_min:i_litr_max
                iso_cnveg_cf.phenology_c_to_litr_c_col[c, j, i] +=
                    iso_cnveg_cf.leafc_to_litter_patch[p] * lf_f[ivt[p], i] * wtcol[p] * leaf_prof[p, j] +
                    iso_cnveg_cf.frootc_to_litter_patch[p] * fr_f[ivt[p], i] * wtcol[p] * froot_prof[p, j]
            end

            if use_crop && ivt[p] >= npcropmin
                # stem litter carbon fluxes
                for i in i_litr_min:i_litr_max
                    iso_cnveg_cf.phenology_c_to_litr_c_col[c, j, i] +=
                        iso_cnveg_cf.livestemc_to_litter_patch[p] * lf_f[ivt[p], i] * wtcol[p] * leaf_prof[p, j]
                end

                if !use_grainproduct
                    for i in i_litr_min:i_litr_max
                        for k in repr_grain_min:repr_grain_max
                            iso_cnveg_cf.phenology_c_to_litr_c_col[c, j, i] +=
                                iso_cnveg_cf.repr_grainc_to_food_patch[p, k] * lf_f[ivt[p], i] * wtcol[p] * leaf_prof[p, j]
                        end
                    end
                end

                # reproductive structure litter carbon fluxes
                for i in i_litr_min:i_litr_max
                    for k in repr_structure_min:repr_structure_max
                        iso_cnveg_cf.phenology_c_to_litr_c_col[c, j, i] +=
                            iso_cnveg_cf.repr_structurec_to_litter_patch[p, k] * lf_f[ivt[p], i] * wtcol[p] * leaf_prof[p, j]
                    end
                end
            end
        end
    end

    return nothing
end

# ---------------------------------------------------------------------------
# cn_c_iso_gap_pft_to_column! — Gap mortality patch-to-column
# Ported from: CNCIsoGapPftToColumn in CNCIsoFluxMod.F90
# ---------------------------------------------------------------------------

function cn_c_iso_gap_pft_to_column!(iso_cnveg_cf::CNVegCarbonFluxData,
                                     soilbiogeochem_state::SoilBiogeochemStateData,
                                     mask_soilp::BitVector,
                                     bounds_p::UnitRange{Int},
                                     nlevdecomp::Int;
                                     patch_column::Vector{Int}=Int[],
                                     patch_itype::Vector{Int}=Int[],
                                     patch_wtcol::Vector{Float64}=Float64[],
                                     lf_f::Matrix{Float64}=Matrix{Float64}(undef, 0, 0),
                                     fr_f::Matrix{Float64}=Matrix{Float64}(undef, 0, 0),
                                     i_litr_min::Int=varpar.i_litr_min,
                                     i_litr_max::Int=varpar.i_litr_max,
                                     i_met_lit::Int=varpar.i_met_lit)

    ivt = patch_itype
    wtcol = patch_wtcol
    leaf_prof = soilbiogeochem_state.leaf_prof_patch
    froot_prof = soilbiogeochem_state.froot_prof_patch
    croot_prof = soilbiogeochem_state.croot_prof_patch
    stem_prof = soilbiogeochem_state.stem_prof_patch

    for j in 1:nlevdecomp
        for p in bounds_p
            mask_soilp[p] || continue
            c = patch_column[p]

            for i in i_litr_min:i_litr_max
                # leaf gap mortality
                iso_cnveg_cf.gap_mortality_c_to_litr_c_col[c, j, i] +=
                    iso_cnveg_cf.m_leafc_to_litter_patch[p] * lf_f[ivt[p], i] * wtcol[p] * leaf_prof[p, j]
                # fine root gap mortality
                iso_cnveg_cf.gap_mortality_c_to_litr_c_col[c, j, i] +=
                    iso_cnveg_cf.m_frootc_to_litter_patch[p] * fr_f[ivt[p], i] * wtcol[p] * froot_prof[p, j]
            end

            # wood gap mortality to CWD
            iso_cnveg_cf.gap_mortality_c_to_cwdc_col[c, j] +=
                iso_cnveg_cf.m_livestemc_to_litter_patch[p] * wtcol[p] * stem_prof[p, j]
            iso_cnveg_cf.gap_mortality_c_to_cwdc_col[c, j] +=
                iso_cnveg_cf.m_deadstemc_to_litter_patch[p] * wtcol[p] * stem_prof[p, j]
            iso_cnveg_cf.gap_mortality_c_to_cwdc_col[c, j] +=
                iso_cnveg_cf.m_livecrootc_to_litter_patch[p] * wtcol[p] * croot_prof[p, j]
            iso_cnveg_cf.gap_mortality_c_to_cwdc_col[c, j] +=
                iso_cnveg_cf.m_deadcrootc_to_litter_patch[p] * wtcol[p] * croot_prof[p, j]

            # metabolic litter (storage + transfer gap mortality)
            iso_cnveg_cf.gap_mortality_c_to_litr_c_col[c, j, i_met_lit] +=
                iso_cnveg_cf.m_leafc_storage_to_litter_patch[p] * wtcol[p] * leaf_prof[p, j] +
                iso_cnveg_cf.m_frootc_storage_to_litter_patch[p] * wtcol[p] * froot_prof[p, j] +
                iso_cnveg_cf.m_livestemc_storage_to_litter_patch[p] * wtcol[p] * stem_prof[p, j] +
                iso_cnveg_cf.m_deadstemc_storage_to_litter_patch[p] * wtcol[p] * stem_prof[p, j] +
                iso_cnveg_cf.m_livecrootc_storage_to_litter_patch[p] * wtcol[p] * croot_prof[p, j] +
                iso_cnveg_cf.m_deadcrootc_storage_to_litter_patch[p] * wtcol[p] * croot_prof[p, j] +
                iso_cnveg_cf.m_gresp_storage_to_litter_patch[p] * wtcol[p] * leaf_prof[p, j] +
                iso_cnveg_cf.m_leafc_xfer_to_litter_patch[p] * wtcol[p] * leaf_prof[p, j] +
                iso_cnveg_cf.m_frootc_xfer_to_litter_patch[p] * wtcol[p] * froot_prof[p, j] +
                iso_cnveg_cf.m_livestemc_xfer_to_litter_patch[p] * wtcol[p] * stem_prof[p, j] +
                iso_cnveg_cf.m_deadstemc_xfer_to_litter_patch[p] * wtcol[p] * stem_prof[p, j] +
                iso_cnveg_cf.m_livecrootc_xfer_to_litter_patch[p] * wtcol[p] * croot_prof[p, j] +
                iso_cnveg_cf.m_deadcrootc_xfer_to_litter_patch[p] * wtcol[p] * croot_prof[p, j] +
                iso_cnveg_cf.m_gresp_xfer_to_litter_patch[p] * wtcol[p] * leaf_prof[p, j]
        end
    end

    return nothing
end

# ---------------------------------------------------------------------------
# cn_c_iso_harvest_pft_to_column! — Harvest mortality patch-to-column
# Ported from: CNCIsoHarvestPftToColumn in CNCIsoFluxMod.F90
# ---------------------------------------------------------------------------

function cn_c_iso_harvest_pft_to_column!(iso_cnveg_cf::CNVegCarbonFluxData,
                                         soilbiogeochem_state::SoilBiogeochemStateData,
                                         mask_soilp::BitVector,
                                         bounds_p::UnitRange{Int},
                                         nlevdecomp::Int;
                                         patch_column::Vector{Int}=Int[],
                                         patch_itype::Vector{Int}=Int[],
                                         patch_wtcol::Vector{Float64}=Float64[],
                                         lf_f::Matrix{Float64}=Matrix{Float64}(undef, 0, 0),
                                         fr_f::Matrix{Float64}=Matrix{Float64}(undef, 0, 0),
                                         i_litr_min::Int=varpar.i_litr_min,
                                         i_litr_max::Int=varpar.i_litr_max,
                                         i_met_lit::Int=varpar.i_met_lit)

    ivt = patch_itype
    wtcol = patch_wtcol
    leaf_prof = soilbiogeochem_state.leaf_prof_patch
    froot_prof = soilbiogeochem_state.froot_prof_patch
    croot_prof = soilbiogeochem_state.croot_prof_patch
    stem_prof = soilbiogeochem_state.stem_prof_patch

    for j in 1:nlevdecomp
        for p in bounds_p
            mask_soilp[p] || continue
            c = patch_column[p]

            for i in i_litr_min:i_litr_max
                iso_cnveg_cf.harvest_c_to_litr_c_col[c, j, i] +=
                    iso_cnveg_cf.hrv_leafc_to_litter_patch[p] * lf_f[ivt[p], i] * wtcol[p] * leaf_prof[p, j]
                iso_cnveg_cf.harvest_c_to_litr_c_col[c, j, i] +=
                    iso_cnveg_cf.hrv_frootc_to_litter_patch[p] * fr_f[ivt[p], i] * wtcol[p] * froot_prof[p, j]
            end

            # wood harvest to CWD
            iso_cnveg_cf.harvest_c_to_cwdc_col[c, j] +=
                iso_cnveg_cf.hrv_livestemc_to_litter_patch[p] * wtcol[p] * stem_prof[p, j]
            iso_cnveg_cf.harvest_c_to_cwdc_col[c, j] +=
                iso_cnveg_cf.hrv_livecrootc_to_litter_patch[p] * wtcol[p] * croot_prof[p, j]
            iso_cnveg_cf.harvest_c_to_cwdc_col[c, j] +=
                iso_cnveg_cf.hrv_deadcrootc_to_litter_patch[p] * wtcol[p] * croot_prof[p, j]

            # metabolic litter (storage + transfer harvest mortality)
            iso_cnveg_cf.harvest_c_to_litr_c_col[c, j, i_met_lit] +=
                iso_cnveg_cf.hrv_leafc_storage_to_litter_patch[p] * wtcol[p] * leaf_prof[p, j] +
                iso_cnveg_cf.hrv_frootc_storage_to_litter_patch[p] * wtcol[p] * froot_prof[p, j] +
                iso_cnveg_cf.hrv_livestemc_storage_to_litter_patch[p] * wtcol[p] * stem_prof[p, j] +
                iso_cnveg_cf.hrv_deadstemc_storage_to_litter_patch[p] * wtcol[p] * stem_prof[p, j] +
                iso_cnveg_cf.hrv_livecrootc_storage_to_litter_patch[p] * wtcol[p] * croot_prof[p, j] +
                iso_cnveg_cf.hrv_deadcrootc_storage_to_litter_patch[p] * wtcol[p] * croot_prof[p, j] +
                iso_cnveg_cf.hrv_gresp_storage_to_litter_patch[p] * wtcol[p] * leaf_prof[p, j] +
                iso_cnveg_cf.hrv_leafc_xfer_to_litter_patch[p] * wtcol[p] * leaf_prof[p, j] +
                iso_cnveg_cf.hrv_frootc_xfer_to_litter_patch[p] * wtcol[p] * froot_prof[p, j] +
                iso_cnveg_cf.hrv_livestemc_xfer_to_litter_patch[p] * wtcol[p] * stem_prof[p, j] +
                iso_cnveg_cf.hrv_deadstemc_xfer_to_litter_patch[p] * wtcol[p] * stem_prof[p, j] +
                iso_cnveg_cf.hrv_livecrootc_xfer_to_litter_patch[p] * wtcol[p] * croot_prof[p, j] +
                iso_cnveg_cf.hrv_deadcrootc_xfer_to_litter_patch[p] * wtcol[p] * croot_prof[p, j] +
                iso_cnveg_cf.hrv_gresp_xfer_to_litter_patch[p] * wtcol[p] * leaf_prof[p, j]
        end
    end

    # wood harvest to column
    for p in bounds_p
        mask_soilp[p] || continue
        c = patch_column[p]
        iso_cnveg_cf.wood_harvestc_col[c] +=
            iso_cnveg_cf.wood_harvestc_patch[p] * wtcol[p]
    end

    return nothing
end

# ---------------------------------------------------------------------------
# cn_c_iso_gross_unrep_pft_to_column! — Gross unrep. LC change patch-to-column
# Ported from: CNCIsoGrossUnrepPftToColumn in CNCIsoFluxMod.F90
# ---------------------------------------------------------------------------

function cn_c_iso_gross_unrep_pft_to_column!(iso_cnveg_cf::CNVegCarbonFluxData,
                                              soilbiogeochem_state::SoilBiogeochemStateData,
                                              mask_soilp::BitVector,
                                              bounds_p::UnitRange{Int},
                                              nlevdecomp::Int;
                                              patch_column::Vector{Int}=Int[],
                                              patch_itype::Vector{Int}=Int[],
                                              patch_wtcol::Vector{Float64}=Float64[],
                                              lf_f::Matrix{Float64}=Matrix{Float64}(undef, 0, 0),
                                              fr_f::Matrix{Float64}=Matrix{Float64}(undef, 0, 0),
                                              i_litr_min::Int=varpar.i_litr_min,
                                              i_litr_max::Int=varpar.i_litr_max)

    ivt = patch_itype
    wtcol = patch_wtcol
    leaf_prof = soilbiogeochem_state.leaf_prof_patch
    froot_prof = soilbiogeochem_state.froot_prof_patch
    croot_prof = soilbiogeochem_state.croot_prof_patch

    for j in 1:nlevdecomp
        for p in bounds_p
            mask_soilp[p] || continue
            c = patch_column[p]

            for i in i_litr_min:i_litr_max
                iso_cnveg_cf.gru_c_to_litr_c_col[c, j, i] +=
                    iso_cnveg_cf.gru_leafc_to_litter_patch[p] * lf_f[ivt[p], i] * wtcol[p] * leaf_prof[p, j] +
                    iso_cnveg_cf.gru_frootc_to_litter_patch[p] * fr_f[ivt[p], i] * wtcol[p] * froot_prof[p, j]
            end

            # coarse root to CWD
            iso_cnveg_cf.gru_c_to_cwdc_col[c, j] +=
                iso_cnveg_cf.gru_livecrootc_to_litter_patch[p] * wtcol[p] * croot_prof[p, j]
            iso_cnveg_cf.gru_c_to_cwdc_col[c, j] +=
                iso_cnveg_cf.gru_deadcrootc_to_litter_patch[p] * wtcol[p] * croot_prof[p, j]
        end
    end

    # wood product gain to column
    for p in bounds_p
        mask_soilp[p] || continue
        c = patch_column[p]
        iso_cnveg_cf.gru_wood_productc_gain_col[c] +=
            iso_cnveg_cf.gru_wood_productc_gain_patch[p] * wtcol[p]
    end

    return nothing
end
