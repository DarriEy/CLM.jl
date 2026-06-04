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
# GPU kernels (KernelAbstractions) for fully-independent per-element loops.
# ---------------------------------------------------------------------------

# Per-element isotopic flux: ciso_flux[i] = ctot_flux[i]*(ciso_state[i]/ctot_state[i])*frax,
# with the exact zero-state guard preserved. One thread per index i in cmin:cmax
# (gated by the mask). frax is a precomputed scalar (C13/C14 fractionation).
@kernel function _ciso_flux_calc_kernel!(ciso_flux, @Const(ctot_flux),
                                         @Const(ciso_state), @Const(ctot_state),
                                         @Const(mask), cmin::Int, cmax::Int, frax)
    i = @index(Global)
    @inbounds if cmin <= i <= cmax && mask[i]
        if ctot_state[i] != 0.0 && ciso_state[i] != 0.0
            ciso_flux[i] = ctot_flux[i] * (ciso_state[i] / ctot_state[i]) * frax
        else
            ciso_flux[i] = 0.0
        end
    end
end

ciso_flux_calc!(ciso_flux, ctot_flux, ciso_state, ctot_state, mask, cmin::Int, cmax::Int, frax) =
    _launch!(_ciso_flux_calc_kernel!, ciso_flux, ctot_flux, ciso_state, ctot_state,
             mask, cmin, cmax, frax; ndrange = length(ciso_flux))

# Column-level decomposition isotopic flux over (c, j, l):
#   out[c,j,l] = tot_flux[c,j,l] * (iso_cpool[c,j,cdp]/tot_cpool[c,j,cdp])
# guarded by tot_cpool[c,j,cdp] != 0, with cdp = cascade_donor_pool[l].
# One thread per (c,j,l); column gated by cmin:cmax and the mask.
@kernel function _ciso_decomp_cascade_kernel!(out, @Const(tot_flux),
                                              @Const(iso_cpool), @Const(tot_cpool),
                                              @Const(mask), @Const(cascade_donor_pool),
                                              cmin::Int, cmax::Int)
    c, j, l = @index(Global, NTuple)
    @inbounds if cmin <= c <= cmax && mask[c]
        cdp = cascade_donor_pool[l]
        if tot_cpool[c, j, cdp] != 0.0
            out[c, j, l] = tot_flux[c, j, l] * (iso_cpool[c, j, cdp] / tot_cpool[c, j, cdp])
        else
            out[c, j, l] = 0.0
        end
    end
end

function ciso_decomp_cascade!(out, tot_flux, iso_cpool, tot_cpool, mask,
                              cascade_donor_pool, cmin::Int, cmax::Int,
                              nlevdecomp::Int, ndecomp_cascade_transitions::Int)
    _launch!(_ciso_decomp_cascade_kernel!, out, tot_flux, iso_cpool, tot_cpool,
             mask, cascade_donor_pool, cmin, cmax;
             ndrange = (size(out, 1), nlevdecomp, ndecomp_cascade_transitions))
end

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
                          frax_c13::Real,
                          isotope::String)
    # if C14, double the fractionation
    if isotope == "c14"
        frax = 1.0 + (1.0 - frax_c13) * 2.0
    elseif isotope == "c13"
        frax = frax_c13
    else
        error("c_iso_flux_calc!: isotope must be either \"c13\" or \"c14\", got \"$isotope\"")
    end

    isempty(bounds) && return nothing
    ciso_flux_calc!(ciso_flux, ctot_flux, ciso_state, ctot_state,
                    mask, first(bounds), last(bounds), frax)
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
                                  frax_c13::Real,
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
                                   frax_c13::Real,
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
# Per-patch crop_harvestc_to_cropprodc init (own-index, no scatter):
#   out[p] = leafc_to_biofuelc[p] + livestemc_to_biofuelc[p]
#          + leafc_to_removedresiduec[p] + livestemc_to_removedresiduec[p]
@kernel function _ciso_crop_harvestc_init_kernel!(out,
                                                  @Const(leafc_to_biofuelc),
                                                  @Const(livestemc_to_biofuelc),
                                                  @Const(leafc_to_removedresiduec),
                                                  @Const(livestemc_to_removedresiduec),
                                                  @Const(mask_soilp), pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax && mask_soilp[p]
        out[p] = leafc_to_biofuelc[p] +
                 livestemc_to_biofuelc[p] +
                 leafc_to_removedresiduec[p] +
                 livestemc_to_removedresiduec[p]
    end
end

# Per-patch accumulation of repr_grainc_to_food over k (own-index += over k):
#   out[p] += sum_{k=kmin:kmax} repr_grainc_to_food[p,k]
@kernel function _ciso_crop_harvestc_add_grain_kernel!(out, @Const(repr_grainc_to_food),
                                                       @Const(mask_soilp),
                                                       kmin::Int, kmax::Int,
                                                       pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax && mask_soilp[p]
        for k in kmin:kmax
            out[p] += repr_grainc_to_food[p, k]
        end
    end
end

# Per-patch reproductive_mr[p,k] = reproductive_xsmr[p,k] + reproductive_curmr[p,k]
# over k (own-index per (p,k); k-loop stays inside the thread):
@kernel function _ciso_reproductive_mr_kernel!(reproductive_mr,
                                               @Const(reproductive_xsmr),
                                               @Const(reproductive_curmr),
                                               @Const(mask_soilp),
                                               kmin::Int, kmax::Int,
                                               pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax && mask_soilp[p]
        for k in kmin:kmax
            reproductive_mr[p, k] =
                reproductive_xsmr[p, k] + reproductive_curmr[p, k]
        end
    end
end

# Per-patch accumulation of repr_structurec_to_cropprod over k (own-index += over k):
#   out[p] += sum_{k=kmin:kmax} repr_structurec_to_cropprod[p,k]
@kernel function _ciso_crop_harvestc_add_struct_kernel!(out, @Const(repr_structurec_to_cropprod),
                                                        @Const(mask_soilp),
                                                        kmin::Int, kmax::Int,
                                                        pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax && mask_soilp[p]
        for k in kmin:kmax
            out[p] += repr_structurec_to_cropprod[p, k]
        end
    end
end

function c_iso_flux1!(soilbiogeochem_state::SoilBiogeochemStateData,
                      soilbiogeochem_cf::SoilBiogeochemCarbonFluxData,
                      soilbiogeochem_cs::SoilBiogeochemCarbonStateData,
                      cnveg_cf::CNVegCarbonFluxData,
                      cnveg_cs::CNVegCarbonStateData,
                      iso_soilbiogeochem_cf::SoilBiogeochemCarbonFluxData,
                      iso_soilbiogeochem_cs::SoilBiogeochemCarbonStateData,
                      iso_cnveg_cf::CNVegCarbonFluxData,
                      iso_cnveg_cs::CNVegCarbonStateData,
                      mask_soilc::AbstractVector{Bool},
                      mask_soilp::AbstractVector{Bool},
                      bounds_c::UnitRange{Int},
                      bounds_p::UnitRange{Int},
                      cascade_donor_pool::AbstractVector{<:Integer},
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
                      patch_wtcol::Vector{<:Real}=Float64[],
                      lf_f::Matrix{<:Real}=Matrix{Float64}(undef, 0, 0),
                      fr_f::Matrix{<:Real}=Matrix{Float64}(undef, 0, 0),
                      leaf_prof::Matrix{<:Real}=Matrix{Float64}(undef, 0, 0),
                      froot_prof::Matrix{<:Real}=Matrix{Float64}(undef, 0, 0))

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

        # crop_harvestc_to_cropprodc_patch (init)
        pmin_p = first(bounds_p)
        pmax_p = last(bounds_p)
        _launch!(_ciso_crop_harvestc_init_kernel!,
                 iso_cnveg_cf.crop_harvestc_to_cropprodc_patch,
                 iso_cnveg_cf.leafc_to_biofuelc_patch,
                 iso_cnveg_cf.livestemc_to_biofuelc_patch,
                 iso_cnveg_cf.leafc_to_removedresiduec_patch,
                 iso_cnveg_cf.livestemc_to_removedresiduec_patch,
                 mask_soilp, pmin_p, pmax_p;
                 ndrange = length(mask_soilp))

        if use_grainproduct
            _launch!(_ciso_crop_harvestc_add_grain_kernel!,
                     iso_cnveg_cf.crop_harvestc_to_cropprodc_patch,
                     iso_cnveg_cf.repr_grainc_to_food_patch,
                     mask_soilp, repr_grain_min, repr_grain_max, pmin_p, pmax_p;
                     ndrange = length(mask_soilp))
        end

        _launch!(_ciso_reproductive_mr_kernel!,
                 iso_cnveg_cf.reproductive_mr_patch,
                 iso_cnveg_cf.reproductive_xsmr_patch,
                 iso_cnveg_cf.reproductive_curmr_patch,
                 mask_soilp, 1, nrepr, pmin_p, pmax_p;
                 ndrange = length(mask_soilp))

        _launch!(_ciso_crop_harvestc_add_struct_kernel!,
                 iso_cnveg_cf.crop_harvestc_to_cropprodc_patch,
                 iso_cnveg_cf.repr_structurec_to_cropprod_patch,
                 mask_soilp, repr_structure_min, repr_structure_max, pmin_p, pmax_p;
                 ndrange = length(mask_soilp))
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
    if !isempty(bounds_c)
        cmin_c = first(bounds_c)
        cmax_c = last(bounds_c)

        ciso_decomp_cascade!(iso_soilbiogeochem_cf.decomp_cascade_hr_vr_col,
                             soilbiogeochem_cf.decomp_cascade_hr_vr_col,
                             iso_soilbiogeochem_cs.decomp_cpools_vr_col,
                             soilbiogeochem_cs.decomp_cpools_vr_col,
                             mask_soilc, cascade_donor_pool, cmin_c, cmax_c,
                             nlevdecomp, ndecomp_cascade_transitions)

        ciso_decomp_cascade!(iso_soilbiogeochem_cf.decomp_cascade_ctransfer_vr_col,
                             soilbiogeochem_cf.decomp_cascade_ctransfer_vr_col,
                             iso_soilbiogeochem_cs.decomp_cpools_vr_col,
                             soilbiogeochem_cs.decomp_cpools_vr_col,
                             mask_soilc, cascade_donor_pool, cmin_c, cmax_c,
                             nlevdecomp, ndecomp_cascade_transitions)
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
                      patch_wtcol::Vector{<:Real}=Float64[],
                      lf_f::Matrix{<:Real}=Matrix{Float64}(undef, 0, 0),
                      fr_f::Matrix{<:Real}=Matrix{Float64}(undef, 0, 0))

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
                       patch_wtcol::Vector{<:Real}=Float64[],
                       lf_f::Matrix{<:Real}=Matrix{Float64}(undef, 0, 0),
                       fr_f::Matrix{<:Real}=Matrix{Float64}(undef, 0, 0))

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
                       patch_wtcol::Vector{<:Real}=Float64[],
                       lf_f::Matrix{<:Real}=Matrix{Float64}(undef, 0, 0),
                       fr_f::Matrix{<:Real}=Matrix{Float64}(undef, 0, 0))

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

# Per-patch fire CWD scatter + decomp-cpool-to-fire vertical-resolved fraction.
# One thread per patch p (gated by mask); the internal j/l loops accumulate into
# this patch's column cc. The two CWD `+=` and the decomp `=` assignment are kept
# as separate scatter/assign ops matching the host's statement order, so the CPU
# float accumulation is byte-identical. (decomp `=` is column-only-valued, so the
# many-patch overwrite is benign.)
@kernel function _ciso_flux3_cwdfire_kernel!(
        fire_mortality_c_to_cwdc_col, m_decomp_cpools_to_fire_vr_col,
        @Const(m_deadstemc_to_litter_fire_patch), @Const(m_livestemc_to_litter_fire_patch),
        @Const(m_deadcrootc_to_litter_fire_patch), @Const(m_livecrootc_to_litter_fire_patch),
        @Const(stem_prof), @Const(croot_prof),
        @Const(tot_m_decomp_cpools_to_fire_vr_col),
        @Const(iso_decomp_cpools_vr_col), @Const(tot_decomp_cpools_vr_col),
        @Const(mask_soilp), @Const(patch_column), @Const(wtcol),
        nlevdecomp::Int, ndecomp_pools::Int, pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax && mask_soilp[p]
        cc = patch_column[p]
        for j in 1:nlevdecomp
            _scatter_add!(fire_mortality_c_to_cwdc_col, cc, j,
                (m_deadstemc_to_litter_fire_patch[p] +
                 m_livestemc_to_litter_fire_patch[p]) *
                wtcol[p] * stem_prof[p, j])
            _scatter_add!(fire_mortality_c_to_cwdc_col, cc, j,
                (m_deadcrootc_to_litter_fire_patch[p] +
                 m_livecrootc_to_litter_fire_patch[p]) *
                wtcol[p] * croot_prof[p, j])
        end

        for j in 1:nlevdecomp
            for l in 1:ndecomp_pools
                if tot_decomp_cpools_vr_col[cc, j, l] != 0.0
                    m_decomp_cpools_to_fire_vr_col[cc, j, l] =
                        tot_m_decomp_cpools_to_fire_vr_col[cc, j, l] *
                        (iso_decomp_cpools_vr_col[cc, j, l] /
                         tot_decomp_cpools_vr_col[cc, j, l])
                else
                    m_decomp_cpools_to_fire_vr_col[cc, j, l] = 0.0
                end
            end
        end
    end
end

# Per-patch fire litter scatter into the metabolic + non-metabolic litter pools.
# One thread per patch p; internal j (and the inner i) loops accumulate into this
# patch's column cc. Byte-identical accumulation order to the host `for p: for j`.
@kernel function _ciso_flux3_litrfire_kernel!(
        m_c_to_litr_fire_col,
        @Const(m_leafc_to_litter_fire_patch), @Const(m_leafc_storage_to_litter_fire_patch),
        @Const(m_leafc_xfer_to_litter_fire_patch), @Const(m_gresp_storage_to_litter_fire_patch),
        @Const(m_gresp_xfer_to_litter_fire_patch), @Const(m_frootc_to_litter_fire_patch),
        @Const(m_frootc_storage_to_litter_fire_patch), @Const(m_frootc_xfer_to_litter_fire_patch),
        @Const(m_livestemc_storage_to_litter_fire_patch), @Const(m_livestemc_xfer_to_litter_fire_patch),
        @Const(m_deadstemc_storage_to_litter_fire_patch), @Const(m_deadstemc_xfer_to_litter_fire_patch),
        @Const(m_livecrootc_storage_to_litter_fire_patch), @Const(m_livecrootc_xfer_to_litter_fire_patch),
        @Const(m_deadcrootc_storage_to_litter_fire_patch), @Const(m_deadcrootc_xfer_to_litter_fire_patch),
        @Const(lf_f), @Const(fr_f), @Const(leaf_prof), @Const(froot_prof),
        @Const(stem_prof), @Const(croot_prof),
        @Const(mask_soilp), @Const(ivt), @Const(patch_column), @Const(wtcol),
        nlevdecomp::Int, i_met_lit::Int, i_litr_max::Int, pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax && mask_soilp[p]
        cc = patch_column[p]
        for j in 1:nlevdecomp
            _scatter_add!(m_c_to_litr_fire_col, cc, j, i_met_lit,
                ((m_leafc_to_litter_fire_patch[p] * lf_f[ivt[p] + 1, i_met_lit] +
                  m_leafc_storage_to_litter_fire_patch[p] +
                  m_leafc_xfer_to_litter_fire_patch[p] +
                  m_gresp_storage_to_litter_fire_patch[p] +
                  m_gresp_xfer_to_litter_fire_patch[p]) * leaf_prof[p, j] +
                 (m_frootc_to_litter_fire_patch[p] * fr_f[ivt[p] + 1, i_met_lit] +
                  m_frootc_storage_to_litter_fire_patch[p] +
                  m_frootc_xfer_to_litter_fire_patch[p]) * froot_prof[p, j] +
                 (m_livestemc_storage_to_litter_fire_patch[p] +
                  m_livestemc_xfer_to_litter_fire_patch[p] +
                  m_deadstemc_storage_to_litter_fire_patch[p] +
                  m_deadstemc_xfer_to_litter_fire_patch[p]) * stem_prof[p, j] +
                 (m_livecrootc_storage_to_litter_fire_patch[p] +
                  m_livecrootc_xfer_to_litter_fire_patch[p] +
                  m_deadcrootc_storage_to_litter_fire_patch[p] +
                  m_deadcrootc_xfer_to_litter_fire_patch[p]) * croot_prof[p, j]) * wtcol[p])

            for i in (i_met_lit+1):i_litr_max
                _scatter_add!(m_c_to_litr_fire_col, cc, j, i,
                    (m_leafc_to_litter_fire_patch[p] * lf_f[ivt[p] + 1, i] * leaf_prof[p, j] +
                     m_frootc_to_litter_fire_patch[p] * fr_f[ivt[p] + 1, i] * froot_prof[p, j]) * wtcol[p])
            end
        end
    end
end

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
                      mask_soilp::AbstractVector{Bool},
                      bounds_p::UnitRange{Int},
                      nlevdecomp::Int,
                      ndecomp_pools::Int,
                      i_met_lit::Int,
                      i_litr_max::Int,
                      isotope::String;
                      patch_column::AbstractVector{<:Integer}=Int[],
                      patch_itype::AbstractVector{<:Integer}=Int[],
                      patch_wtcol::AbstractVector{<:Real}=Float64[],
                      lf_f::AbstractMatrix{<:Real}=Matrix{Float64}(undef, 0, 0),
                      fr_f::AbstractMatrix{<:Real}=Matrix{Float64}(undef, 0, 0),
                      leaf_prof::AbstractMatrix{<:Real}=Matrix{Float64}(undef, 0, 0),
                      froot_prof::AbstractMatrix{<:Real}=Matrix{Float64}(undef, 0, 0),
                      stem_prof::AbstractMatrix{<:Real}=Matrix{Float64}(undef, 0, 0),
                      croot_prof::AbstractMatrix{<:Real}=Matrix{Float64}(undef, 0, 0))

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

    if !isempty(bounds_p)
        _launch!(_ciso_flux3_cwdfire_kernel!,
            iso_cnveg_cf.fire_mortality_c_to_cwdc_col,
            iso_cnveg_cf.m_decomp_cpools_to_fire_vr_col,
            iso_cnveg_cf.m_deadstemc_to_litter_fire_patch,
            iso_cnveg_cf.m_livestemc_to_litter_fire_patch,
            iso_cnveg_cf.m_deadcrootc_to_litter_fire_patch,
            iso_cnveg_cf.m_livecrootc_to_litter_fire_patch,
            stem_prof, croot_prof,
            cnveg_cf.m_decomp_cpools_to_fire_vr_col,
            iso_soilbiogeochem_cs.decomp_cpools_vr_col,
            soilbiogeochem_cs.decomp_cpools_vr_col,
            mask_soilp, patch_column, wtcol,
            nlevdecomp, ndecomp_pools, first(bounds_p), last(bounds_p);
            ndrange = length(mask_soilp))

        _launch!(_ciso_flux3_litrfire_kernel!,
            iso_cnveg_cf.m_c_to_litr_fire_col,
            iso_cnveg_cf.m_leafc_to_litter_fire_patch,
            iso_cnveg_cf.m_leafc_storage_to_litter_fire_patch,
            iso_cnveg_cf.m_leafc_xfer_to_litter_fire_patch,
            iso_cnveg_cf.m_gresp_storage_to_litter_fire_patch,
            iso_cnveg_cf.m_gresp_xfer_to_litter_fire_patch,
            iso_cnveg_cf.m_frootc_to_litter_fire_patch,
            iso_cnveg_cf.m_frootc_storage_to_litter_fire_patch,
            iso_cnveg_cf.m_frootc_xfer_to_litter_fire_patch,
            iso_cnveg_cf.m_livestemc_storage_to_litter_fire_patch,
            iso_cnveg_cf.m_livestemc_xfer_to_litter_fire_patch,
            iso_cnveg_cf.m_deadstemc_storage_to_litter_fire_patch,
            iso_cnveg_cf.m_deadstemc_xfer_to_litter_fire_patch,
            iso_cnveg_cf.m_livecrootc_storage_to_litter_fire_patch,
            iso_cnveg_cf.m_livecrootc_xfer_to_litter_fire_patch,
            iso_cnveg_cf.m_deadcrootc_storage_to_litter_fire_patch,
            iso_cnveg_cf.m_deadcrootc_xfer_to_litter_fire_patch,
            lf_f, fr_f, leaf_prof, froot_prof, stem_prof, croot_prof,
            mask_soilp, ivt, patch_column, wtcol,
            nlevdecomp, i_met_lit, i_litr_max, first(bounds_p), last(bounds_p);
            ndrange = length(mask_soilp))
    end

    return nothing
end

# ---------------------------------------------------------------------------
# cn_c_iso_litter_to_column! — Phenology litterfall patch-to-column
# Ported from: CNCIsoLitterToColumn in CNCIsoFluxMod.F90
# ---------------------------------------------------------------------------

# Per-patch phenology litterfall scatter into column litter pools. One thread per
# patch p; the internal j (and inner i / reproductive-k) loops accumulate into this
# patch's column c. The crop / grainproduct config flags are Bool scalars and the
# per-patch `ivt[p] >= npcropmin` branch is kept in-kernel (the gating is interwoven
# with the per-patch type, so it cannot be hoisted to host). Accumulation order is
# byte-identical to the host `for j: for p`.
@kernel function _ciso_litter_p2c_kernel!(
        phenology_c_to_litr_c_col,
        @Const(leafc_to_litter_patch), @Const(frootc_to_litter_patch),
        @Const(livestemc_to_litter_patch),
        @Const(repr_grainc_to_food_patch), @Const(repr_structurec_to_litter_patch),
        @Const(lf_f), @Const(fr_f), @Const(leaf_prof), @Const(froot_prof),
        @Const(mask_soilp), @Const(ivt), @Const(patch_column), @Const(wtcol),
        nlevdecomp::Int, i_litr_min::Int, i_litr_max::Int,
        use_crop::Bool, use_grainproduct::Bool, npcropmin::Int,
        repr_grain_min::Int, repr_grain_max::Int,
        repr_structure_min::Int, repr_structure_max::Int, pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax && mask_soilp[p]
        c = patch_column[p]
        for j in 1:nlevdecomp
            for i in i_litr_min:i_litr_max
                _scatter_add!(phenology_c_to_litr_c_col, c, j, i,
                    leafc_to_litter_patch[p] * lf_f[ivt[p] + 1, i] * wtcol[p] * leaf_prof[p, j] +
                    frootc_to_litter_patch[p] * fr_f[ivt[p] + 1, i] * wtcol[p] * froot_prof[p, j])
            end

            if use_crop && ivt[p] >= npcropmin
                # stem litter carbon fluxes
                for i in i_litr_min:i_litr_max
                    _scatter_add!(phenology_c_to_litr_c_col, c, j, i,
                        livestemc_to_litter_patch[p] * lf_f[ivt[p] + 1, i] * wtcol[p] * leaf_prof[p, j])
                end

                if !use_grainproduct
                    for i in i_litr_min:i_litr_max
                        for k in repr_grain_min:repr_grain_max
                            _scatter_add!(phenology_c_to_litr_c_col, c, j, i,
                                repr_grainc_to_food_patch[p, k] * lf_f[ivt[p] + 1, i] * wtcol[p] * leaf_prof[p, j])
                        end
                    end
                end

                # reproductive structure litter carbon fluxes
                for i in i_litr_min:i_litr_max
                    for k in repr_structure_min:repr_structure_max
                        _scatter_add!(phenology_c_to_litr_c_col, c, j, i,
                            repr_structurec_to_litter_patch[p, k] * lf_f[ivt[p] + 1, i] * wtcol[p] * leaf_prof[p, j])
                    end
                end
            end
        end
    end
end

function cn_c_iso_litter_to_column!(iso_cnveg_cf::CNVegCarbonFluxData,
                                    soilbiogeochem_state::SoilBiogeochemStateData,
                                    mask_soilp::AbstractVector{Bool},
                                    bounds_p::UnitRange{Int},
                                    nlevdecomp::Int;
                                    patch_column::AbstractVector{<:Integer}=Int[],
                                    patch_itype::AbstractVector{<:Integer}=Int[],
                                    patch_wtcol::AbstractVector{<:Real}=Float64[],
                                    lf_f::AbstractMatrix{<:Real}=Matrix{Float64}(undef, 0, 0),
                                    fr_f::AbstractMatrix{<:Real}=Matrix{Float64}(undef, 0, 0),
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

    isempty(bounds_p) && return nothing
    _launch!(_ciso_litter_p2c_kernel!,
        iso_cnveg_cf.phenology_c_to_litr_c_col,
        iso_cnveg_cf.leafc_to_litter_patch,
        iso_cnveg_cf.frootc_to_litter_patch,
        iso_cnveg_cf.livestemc_to_litter_patch,
        iso_cnveg_cf.repr_grainc_to_food_patch,
        iso_cnveg_cf.repr_structurec_to_litter_patch,
        lf_f, fr_f, leaf_prof, froot_prof,
        mask_soilp, ivt, patch_column, wtcol,
        nlevdecomp, i_litr_min, i_litr_max,
        use_crop, use_grainproduct, npcropmin,
        repr_grain_min, repr_grain_max, repr_structure_min, repr_structure_max,
        first(bounds_p), last(bounds_p);
        ndrange = length(mask_soilp))

    return nothing
end

# ---------------------------------------------------------------------------
# cn_c_iso_gap_pft_to_column! — Gap mortality patch-to-column
# Ported from: CNCIsoGapPftToColumn in CNCIsoFluxMod.F90
# ---------------------------------------------------------------------------

# Per-patch gap-mortality scatter into column litter/CWD pools. One thread per
# patch p; internal j (and inner i) loops accumulate into this patch's column c.
# The four CWD `+=` adds are kept as four separate _scatter_add! calls in the host
# statement order, and the metabolic-litter add keeps the host's full RHS sum, so
# CPU float accumulation is byte-identical to the host `for j: for p`.
@kernel function _ciso_gap_p2c_kernel!(
        gap_mortality_c_to_litr_c_col, gap_mortality_c_to_cwdc_col,
        @Const(m_leafc_to_litter_patch), @Const(m_frootc_to_litter_patch),
        @Const(m_livestemc_to_litter_patch), @Const(m_deadstemc_to_litter_patch),
        @Const(m_livecrootc_to_litter_patch), @Const(m_deadcrootc_to_litter_patch),
        @Const(m_leafc_storage_to_litter_patch), @Const(m_frootc_storage_to_litter_patch),
        @Const(m_livestemc_storage_to_litter_patch), @Const(m_deadstemc_storage_to_litter_patch),
        @Const(m_livecrootc_storage_to_litter_patch), @Const(m_deadcrootc_storage_to_litter_patch),
        @Const(m_gresp_storage_to_litter_patch),
        @Const(m_leafc_xfer_to_litter_patch), @Const(m_frootc_xfer_to_litter_patch),
        @Const(m_livestemc_xfer_to_litter_patch), @Const(m_deadstemc_xfer_to_litter_patch),
        @Const(m_livecrootc_xfer_to_litter_patch), @Const(m_deadcrootc_xfer_to_litter_patch),
        @Const(m_gresp_xfer_to_litter_patch),
        @Const(lf_f), @Const(fr_f), @Const(leaf_prof), @Const(froot_prof),
        @Const(croot_prof), @Const(stem_prof),
        @Const(mask_soilp), @Const(ivt), @Const(patch_column), @Const(wtcol),
        nlevdecomp::Int, i_litr_min::Int, i_litr_max::Int, i_met_lit::Int,
        pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax && mask_soilp[p]
        c = patch_column[p]
        for j in 1:nlevdecomp
            for i in i_litr_min:i_litr_max
                # leaf gap mortality
                _scatter_add!(gap_mortality_c_to_litr_c_col, c, j, i,
                    m_leafc_to_litter_patch[p] * lf_f[ivt[p] + 1, i] * wtcol[p] * leaf_prof[p, j])
                # fine root gap mortality
                _scatter_add!(gap_mortality_c_to_litr_c_col, c, j, i,
                    m_frootc_to_litter_patch[p] * fr_f[ivt[p] + 1, i] * wtcol[p] * froot_prof[p, j])
            end

            # wood gap mortality to CWD
            _scatter_add!(gap_mortality_c_to_cwdc_col, c, j,
                m_livestemc_to_litter_patch[p] * wtcol[p] * stem_prof[p, j])
            _scatter_add!(gap_mortality_c_to_cwdc_col, c, j,
                m_deadstemc_to_litter_patch[p] * wtcol[p] * stem_prof[p, j])
            _scatter_add!(gap_mortality_c_to_cwdc_col, c, j,
                m_livecrootc_to_litter_patch[p] * wtcol[p] * croot_prof[p, j])
            _scatter_add!(gap_mortality_c_to_cwdc_col, c, j,
                m_deadcrootc_to_litter_patch[p] * wtcol[p] * croot_prof[p, j])

            # metabolic litter (storage + transfer gap mortality)
            _scatter_add!(gap_mortality_c_to_litr_c_col, c, j, i_met_lit,
                m_leafc_storage_to_litter_patch[p] * wtcol[p] * leaf_prof[p, j] +
                m_frootc_storage_to_litter_patch[p] * wtcol[p] * froot_prof[p, j] +
                m_livestemc_storage_to_litter_patch[p] * wtcol[p] * stem_prof[p, j] +
                m_deadstemc_storage_to_litter_patch[p] * wtcol[p] * stem_prof[p, j] +
                m_livecrootc_storage_to_litter_patch[p] * wtcol[p] * croot_prof[p, j] +
                m_deadcrootc_storage_to_litter_patch[p] * wtcol[p] * croot_prof[p, j] +
                m_gresp_storage_to_litter_patch[p] * wtcol[p] * leaf_prof[p, j] +
                m_leafc_xfer_to_litter_patch[p] * wtcol[p] * leaf_prof[p, j] +
                m_frootc_xfer_to_litter_patch[p] * wtcol[p] * froot_prof[p, j] +
                m_livestemc_xfer_to_litter_patch[p] * wtcol[p] * stem_prof[p, j] +
                m_deadstemc_xfer_to_litter_patch[p] * wtcol[p] * stem_prof[p, j] +
                m_livecrootc_xfer_to_litter_patch[p] * wtcol[p] * croot_prof[p, j] +
                m_deadcrootc_xfer_to_litter_patch[p] * wtcol[p] * croot_prof[p, j] +
                m_gresp_xfer_to_litter_patch[p] * wtcol[p] * leaf_prof[p, j])
        end
    end
end

function cn_c_iso_gap_pft_to_column!(iso_cnveg_cf::CNVegCarbonFluxData,
                                     soilbiogeochem_state::SoilBiogeochemStateData,
                                     mask_soilp::AbstractVector{Bool},
                                     bounds_p::UnitRange{Int},
                                     nlevdecomp::Int;
                                     patch_column::AbstractVector{<:Integer}=Int[],
                                     patch_itype::AbstractVector{<:Integer}=Int[],
                                     patch_wtcol::AbstractVector{<:Real}=Float64[],
                                     lf_f::AbstractMatrix{<:Real}=Matrix{Float64}(undef, 0, 0),
                                     fr_f::AbstractMatrix{<:Real}=Matrix{Float64}(undef, 0, 0),
                                     i_litr_min::Int=varpar.i_litr_min,
                                     i_litr_max::Int=varpar.i_litr_max,
                                     i_met_lit::Int=varpar.i_met_lit)

    ivt = patch_itype
    wtcol = patch_wtcol
    leaf_prof = soilbiogeochem_state.leaf_prof_patch
    froot_prof = soilbiogeochem_state.froot_prof_patch
    croot_prof = soilbiogeochem_state.croot_prof_patch
    stem_prof = soilbiogeochem_state.stem_prof_patch

    isempty(bounds_p) && return nothing
    _launch!(_ciso_gap_p2c_kernel!,
        iso_cnveg_cf.gap_mortality_c_to_litr_c_col,
        iso_cnveg_cf.gap_mortality_c_to_cwdc_col,
        iso_cnveg_cf.m_leafc_to_litter_patch,
        iso_cnveg_cf.m_frootc_to_litter_patch,
        iso_cnveg_cf.m_livestemc_to_litter_patch,
        iso_cnveg_cf.m_deadstemc_to_litter_patch,
        iso_cnveg_cf.m_livecrootc_to_litter_patch,
        iso_cnveg_cf.m_deadcrootc_to_litter_patch,
        iso_cnveg_cf.m_leafc_storage_to_litter_patch,
        iso_cnveg_cf.m_frootc_storage_to_litter_patch,
        iso_cnveg_cf.m_livestemc_storage_to_litter_patch,
        iso_cnveg_cf.m_deadstemc_storage_to_litter_patch,
        iso_cnveg_cf.m_livecrootc_storage_to_litter_patch,
        iso_cnveg_cf.m_deadcrootc_storage_to_litter_patch,
        iso_cnveg_cf.m_gresp_storage_to_litter_patch,
        iso_cnveg_cf.m_leafc_xfer_to_litter_patch,
        iso_cnveg_cf.m_frootc_xfer_to_litter_patch,
        iso_cnveg_cf.m_livestemc_xfer_to_litter_patch,
        iso_cnveg_cf.m_deadstemc_xfer_to_litter_patch,
        iso_cnveg_cf.m_livecrootc_xfer_to_litter_patch,
        iso_cnveg_cf.m_deadcrootc_xfer_to_litter_patch,
        iso_cnveg_cf.m_gresp_xfer_to_litter_patch,
        lf_f, fr_f, leaf_prof, froot_prof, croot_prof, stem_prof,
        mask_soilp, ivt, patch_column, wtcol,
        nlevdecomp, i_litr_min, i_litr_max, i_met_lit,
        first(bounds_p), last(bounds_p);
        ndrange = length(mask_soilp))

    return nothing
end

# ---------------------------------------------------------------------------
# cn_c_iso_harvest_pft_to_column! — Harvest mortality patch-to-column
# Ported from: CNCIsoHarvestPftToColumn in CNCIsoFluxMod.F90
# ---------------------------------------------------------------------------

# Per-patch harvest-mortality scatter into column litter/CWD pools. One thread per
# patch p; internal j (and inner i) loops accumulate into this patch's column c.
# CWD adds kept as separate _scatter_add! ops; metabolic-litter add keeps the full
# host RHS. Byte-identical CPU accumulation to the host `for j: for p`.
@kernel function _ciso_harvest_p2c_kernel!(
        harvest_c_to_litr_c_col, harvest_c_to_cwdc_col,
        @Const(hrv_leafc_to_litter_patch), @Const(hrv_frootc_to_litter_patch),
        @Const(hrv_livestemc_to_litter_patch), @Const(hrv_livecrootc_to_litter_patch),
        @Const(hrv_deadcrootc_to_litter_patch),
        @Const(hrv_leafc_storage_to_litter_patch), @Const(hrv_frootc_storage_to_litter_patch),
        @Const(hrv_livestemc_storage_to_litter_patch), @Const(hrv_deadstemc_storage_to_litter_patch),
        @Const(hrv_livecrootc_storage_to_litter_patch), @Const(hrv_deadcrootc_storage_to_litter_patch),
        @Const(hrv_gresp_storage_to_litter_patch),
        @Const(hrv_leafc_xfer_to_litter_patch), @Const(hrv_frootc_xfer_to_litter_patch),
        @Const(hrv_livestemc_xfer_to_litter_patch), @Const(hrv_deadstemc_xfer_to_litter_patch),
        @Const(hrv_livecrootc_xfer_to_litter_patch), @Const(hrv_deadcrootc_xfer_to_litter_patch),
        @Const(hrv_gresp_xfer_to_litter_patch),
        @Const(lf_f), @Const(fr_f), @Const(leaf_prof), @Const(froot_prof),
        @Const(croot_prof), @Const(stem_prof),
        @Const(mask_soilp), @Const(ivt), @Const(patch_column), @Const(wtcol),
        nlevdecomp::Int, i_litr_min::Int, i_litr_max::Int, i_met_lit::Int,
        pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax && mask_soilp[p]
        c = patch_column[p]
        for j in 1:nlevdecomp
            for i in i_litr_min:i_litr_max
                _scatter_add!(harvest_c_to_litr_c_col, c, j, i,
                    hrv_leafc_to_litter_patch[p] * lf_f[ivt[p] + 1, i] * wtcol[p] * leaf_prof[p, j])
                _scatter_add!(harvest_c_to_litr_c_col, c, j, i,
                    hrv_frootc_to_litter_patch[p] * fr_f[ivt[p] + 1, i] * wtcol[p] * froot_prof[p, j])
            end

            # wood harvest to CWD
            _scatter_add!(harvest_c_to_cwdc_col, c, j,
                hrv_livestemc_to_litter_patch[p] * wtcol[p] * stem_prof[p, j])
            _scatter_add!(harvest_c_to_cwdc_col, c, j,
                hrv_livecrootc_to_litter_patch[p] * wtcol[p] * croot_prof[p, j])
            _scatter_add!(harvest_c_to_cwdc_col, c, j,
                hrv_deadcrootc_to_litter_patch[p] * wtcol[p] * croot_prof[p, j])

            # metabolic litter (storage + transfer harvest mortality)
            _scatter_add!(harvest_c_to_litr_c_col, c, j, i_met_lit,
                hrv_leafc_storage_to_litter_patch[p] * wtcol[p] * leaf_prof[p, j] +
                hrv_frootc_storage_to_litter_patch[p] * wtcol[p] * froot_prof[p, j] +
                hrv_livestemc_storage_to_litter_patch[p] * wtcol[p] * stem_prof[p, j] +
                hrv_deadstemc_storage_to_litter_patch[p] * wtcol[p] * stem_prof[p, j] +
                hrv_livecrootc_storage_to_litter_patch[p] * wtcol[p] * croot_prof[p, j] +
                hrv_deadcrootc_storage_to_litter_patch[p] * wtcol[p] * croot_prof[p, j] +
                hrv_gresp_storage_to_litter_patch[p] * wtcol[p] * leaf_prof[p, j] +
                hrv_leafc_xfer_to_litter_patch[p] * wtcol[p] * leaf_prof[p, j] +
                hrv_frootc_xfer_to_litter_patch[p] * wtcol[p] * froot_prof[p, j] +
                hrv_livestemc_xfer_to_litter_patch[p] * wtcol[p] * stem_prof[p, j] +
                hrv_deadstemc_xfer_to_litter_patch[p] * wtcol[p] * stem_prof[p, j] +
                hrv_livecrootc_xfer_to_litter_patch[p] * wtcol[p] * croot_prof[p, j] +
                hrv_deadcrootc_xfer_to_litter_patch[p] * wtcol[p] * croot_prof[p, j] +
                hrv_gresp_xfer_to_litter_patch[p] * wtcol[p] * leaf_prof[p, j])
        end
    end
end

# Per-patch wood-harvest scatter into the column total (1D output).
@kernel function _ciso_harvest_wood_p2c_kernel!(
        wood_harvestc_col, @Const(wood_harvestc_patch),
        @Const(mask_soilp), @Const(patch_column), @Const(wtcol),
        pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax && mask_soilp[p]
        c = patch_column[p]
        _scatter_add!(wood_harvestc_col, c, wood_harvestc_patch[p] * wtcol[p])
    end
end

function cn_c_iso_harvest_pft_to_column!(iso_cnveg_cf::CNVegCarbonFluxData,
                                         soilbiogeochem_state::SoilBiogeochemStateData,
                                         mask_soilp::AbstractVector{Bool},
                                         bounds_p::UnitRange{Int},
                                         nlevdecomp::Int;
                                         patch_column::AbstractVector{<:Integer}=Int[],
                                         patch_itype::AbstractVector{<:Integer}=Int[],
                                         patch_wtcol::AbstractVector{<:Real}=Float64[],
                                         lf_f::AbstractMatrix{<:Real}=Matrix{Float64}(undef, 0, 0),
                                         fr_f::AbstractMatrix{<:Real}=Matrix{Float64}(undef, 0, 0),
                                         i_litr_min::Int=varpar.i_litr_min,
                                         i_litr_max::Int=varpar.i_litr_max,
                                         i_met_lit::Int=varpar.i_met_lit)

    ivt = patch_itype
    wtcol = patch_wtcol
    leaf_prof = soilbiogeochem_state.leaf_prof_patch
    froot_prof = soilbiogeochem_state.froot_prof_patch
    croot_prof = soilbiogeochem_state.croot_prof_patch
    stem_prof = soilbiogeochem_state.stem_prof_patch

    isempty(bounds_p) && return nothing
    _launch!(_ciso_harvest_p2c_kernel!,
        iso_cnveg_cf.harvest_c_to_litr_c_col,
        iso_cnveg_cf.harvest_c_to_cwdc_col,
        iso_cnveg_cf.hrv_leafc_to_litter_patch,
        iso_cnveg_cf.hrv_frootc_to_litter_patch,
        iso_cnveg_cf.hrv_livestemc_to_litter_patch,
        iso_cnveg_cf.hrv_livecrootc_to_litter_patch,
        iso_cnveg_cf.hrv_deadcrootc_to_litter_patch,
        iso_cnveg_cf.hrv_leafc_storage_to_litter_patch,
        iso_cnveg_cf.hrv_frootc_storage_to_litter_patch,
        iso_cnveg_cf.hrv_livestemc_storage_to_litter_patch,
        iso_cnveg_cf.hrv_deadstemc_storage_to_litter_patch,
        iso_cnveg_cf.hrv_livecrootc_storage_to_litter_patch,
        iso_cnveg_cf.hrv_deadcrootc_storage_to_litter_patch,
        iso_cnveg_cf.hrv_gresp_storage_to_litter_patch,
        iso_cnveg_cf.hrv_leafc_xfer_to_litter_patch,
        iso_cnveg_cf.hrv_frootc_xfer_to_litter_patch,
        iso_cnveg_cf.hrv_livestemc_xfer_to_litter_patch,
        iso_cnveg_cf.hrv_deadstemc_xfer_to_litter_patch,
        iso_cnveg_cf.hrv_livecrootc_xfer_to_litter_patch,
        iso_cnveg_cf.hrv_deadcrootc_xfer_to_litter_patch,
        iso_cnveg_cf.hrv_gresp_xfer_to_litter_patch,
        lf_f, fr_f, leaf_prof, froot_prof, croot_prof, stem_prof,
        mask_soilp, ivt, patch_column, wtcol,
        nlevdecomp, i_litr_min, i_litr_max, i_met_lit,
        first(bounds_p), last(bounds_p);
        ndrange = length(mask_soilp))

    # wood harvest to column
    _launch!(_ciso_harvest_wood_p2c_kernel!,
        iso_cnveg_cf.wood_harvestc_col,
        iso_cnveg_cf.wood_harvestc_patch,
        mask_soilp, patch_column, wtcol,
        first(bounds_p), last(bounds_p);
        ndrange = length(mask_soilp))

    return nothing
end

# ---------------------------------------------------------------------------
# cn_c_iso_gross_unrep_pft_to_column! — Gross unrep. LC change patch-to-column
# Ported from: CNCIsoGrossUnrepPftToColumn in CNCIsoFluxMod.F90
# ---------------------------------------------------------------------------

# Per-patch gross-unrepresented-LC-change scatter into column litter/CWD pools.
# One thread per patch p; internal j (and inner i) loops accumulate into column c.
# Byte-identical CPU accumulation to the host `for j: for p`.
@kernel function _ciso_gru_p2c_kernel!(
        gru_c_to_litr_c_col, gru_c_to_cwdc_col,
        @Const(gru_leafc_to_litter_patch), @Const(gru_frootc_to_litter_patch),
        @Const(gru_livecrootc_to_litter_patch), @Const(gru_deadcrootc_to_litter_patch),
        @Const(lf_f), @Const(fr_f), @Const(leaf_prof), @Const(froot_prof),
        @Const(croot_prof),
        @Const(mask_soilp), @Const(ivt), @Const(patch_column), @Const(wtcol),
        nlevdecomp::Int, i_litr_min::Int, i_litr_max::Int, pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax && mask_soilp[p]
        c = patch_column[p]
        for j in 1:nlevdecomp
            for i in i_litr_min:i_litr_max
                _scatter_add!(gru_c_to_litr_c_col, c, j, i,
                    gru_leafc_to_litter_patch[p] * lf_f[ivt[p] + 1, i] * wtcol[p] * leaf_prof[p, j] +
                    gru_frootc_to_litter_patch[p] * fr_f[ivt[p] + 1, i] * wtcol[p] * froot_prof[p, j])
            end

            # coarse root to CWD
            _scatter_add!(gru_c_to_cwdc_col, c, j,
                gru_livecrootc_to_litter_patch[p] * wtcol[p] * croot_prof[p, j])
            _scatter_add!(gru_c_to_cwdc_col, c, j,
                gru_deadcrootc_to_litter_patch[p] * wtcol[p] * croot_prof[p, j])
        end
    end
end

# Per-patch wood-product-gain scatter into the column total (1D output).
@kernel function _ciso_gru_wood_p2c_kernel!(
        gru_wood_productc_gain_col, @Const(gru_wood_productc_gain_patch),
        @Const(mask_soilp), @Const(patch_column), @Const(wtcol),
        pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax && mask_soilp[p]
        c = patch_column[p]
        _scatter_add!(gru_wood_productc_gain_col, c, gru_wood_productc_gain_patch[p] * wtcol[p])
    end
end

function cn_c_iso_gross_unrep_pft_to_column!(iso_cnveg_cf::CNVegCarbonFluxData,
                                              soilbiogeochem_state::SoilBiogeochemStateData,
                                              mask_soilp::AbstractVector{Bool},
                                              bounds_p::UnitRange{Int},
                                              nlevdecomp::Int;
                                              patch_column::AbstractVector{<:Integer}=Int[],
                                              patch_itype::AbstractVector{<:Integer}=Int[],
                                              patch_wtcol::AbstractVector{<:Real}=Float64[],
                                              lf_f::AbstractMatrix{<:Real}=Matrix{Float64}(undef, 0, 0),
                                              fr_f::AbstractMatrix{<:Real}=Matrix{Float64}(undef, 0, 0),
                                              i_litr_min::Int=varpar.i_litr_min,
                                              i_litr_max::Int=varpar.i_litr_max)

    ivt = patch_itype
    wtcol = patch_wtcol
    leaf_prof = soilbiogeochem_state.leaf_prof_patch
    froot_prof = soilbiogeochem_state.froot_prof_patch
    croot_prof = soilbiogeochem_state.croot_prof_patch

    isempty(bounds_p) && return nothing
    _launch!(_ciso_gru_p2c_kernel!,
        iso_cnveg_cf.gru_c_to_litr_c_col,
        iso_cnveg_cf.gru_c_to_cwdc_col,
        iso_cnveg_cf.gru_leafc_to_litter_patch,
        iso_cnveg_cf.gru_frootc_to_litter_patch,
        iso_cnveg_cf.gru_livecrootc_to_litter_patch,
        iso_cnveg_cf.gru_deadcrootc_to_litter_patch,
        lf_f, fr_f, leaf_prof, froot_prof, croot_prof,
        mask_soilp, ivt, patch_column, wtcol,
        nlevdecomp, i_litr_min, i_litr_max,
        first(bounds_p), last(bounds_p);
        ndrange = length(mask_soilp))

    # wood product gain to column
    _launch!(_ciso_gru_wood_p2c_kernel!,
        iso_cnveg_cf.gru_wood_productc_gain_col,
        iso_cnveg_cf.gru_wood_productc_gain_patch,
        mask_soilp, patch_column, wtcol,
        first(bounds_p), last(bounds_p);
        ndrange = length(mask_soilp))

    return nothing
end
