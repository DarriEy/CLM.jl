# ==========================================================================
# Ported from: src/biogeochem/CNAnnualUpdateMod.F90
# Module for updating annual summation variables
# ==========================================================================

"""
    cn_annual_update!(mask_bgc_soilc, mask_bgc_vegp, bounds_c, bounds_p,
                      col, patch, cnveg_state, cnveg_carbonflux;
                      dt, days_per_year, secspday)

On the radiation time step, update annual summation variables.

For each BGC soil column (not FATES), the annual-sum counter is incremented
by `dt`. When the counter reaches the number of seconds in a year
(`days_per_year * secspday`), the column is flagged as having reached its
end-of-year. For each BGC vegetation patch whose parent column has reached
end-of-year (and is not FATES), the temporary accumulator variables are
promoted to their annual counterparts and then reset to zero.

Finally, a patch-to-column averaging step computes column-level annual NPP
and annual average 2-m air temperature for the end-of-year columns.

Ported from `CNAnnualUpdate` in `CNAnnualUpdateMod.F90`.

# Arguments
- `mask_bgc_soilc::BitVector` : column mask for BGC soil columns
- `mask_bgc_vegp::BitVector`  : patch mask for BGC vegetation patches
- `bounds_c::UnitRange{Int}`  : column index range
- `bounds_p::UnitRange{Int}`  : patch index range
- `col::ColumnData`           : column-level data (is_fates, patchi, patchf)
- `patch::PatchData`          : patch-level data (column, active, wtcol)
- `cnveg_state::CNVegStateData`         : CN vegetation state (accumulators)
- `cnveg_carbonflux::CNVegCarbonFluxData` : CN vegetation carbon fluxes (NPP, litfall accumulators)
- `dt::Float64`               : radiation time step (seconds)
- `days_per_year::Float64`    : number of days in the current year
- `secspday::Float64`         : seconds per day (default `SECSPDAY`)
"""
# --- Per-column end-of-year kernel (loop 1) ---
# Writes end_of_year[c] and mutates annsum_counter_col[c] for BGC soil,
# non-FATES columns inside [c_lo, c_hi]. end_of_year must be pre-zeroed.
@kernel function _cn_annual_eoy_col_kernel!(end_of_year, annsum_counter_col,
                                            @Const(mask_bgc_soilc), @Const(is_fates),
                                            dt, secspyear, c_lo, c_hi)
    T = eltype(annsum_counter_col)
    c = @index(Global)
    @inbounds if c >= c_lo && c <= c_hi && mask_bgc_soilc[c] && !is_fates[c]
        annsum_counter_col[c] += dt
        if annsum_counter_col[c] >= secspyear
            end_of_year[c] = true
            annsum_counter_col[c] = zero(T)
        end
    end
end

# --- Per-patch annual update kernel (loop 2) ---
# For BGC veg patches inside [p_lo, p_hi] whose parent column is at end-of-year
# (and not FATES): promote tempsum/tempmax/tempavg accumulators to their annual
# counterparts (npp/litfall scaled by dt) and reset the temporaries to zero.
@kernel function _cn_annual_patch_kernel!(annsum_potential_gpp_patch,
                                          tempsum_potential_gpp_patch,
                                          annmax_retransn_patch,
                                          tempmax_retransn_patch,
                                          annavg_t2m_patch,
                                          tempavg_t2m_patch,
                                          annsum_npp_patch,
                                          tempsum_npp_patch,
                                          annsum_litfall_patch,
                                          tempsum_litfall_patch,
                                          @Const(mask_bgc_vegp), @Const(pcol),
                                          @Const(end_of_year), @Const(is_fates),
                                          dt, p_lo, p_hi)
    T = eltype(annsum_potential_gpp_patch)
    p = @index(Global)
    @inbounds if p >= p_lo && p <= p_hi && mask_bgc_vegp[p]
        c = pcol[p]
        if end_of_year[c] && !is_fates[c]
            # update annual plant ndemand accumulator
            annsum_potential_gpp_patch[p]  = tempsum_potential_gpp_patch[p]
            tempsum_potential_gpp_patch[p] = zero(T)

            # update annual total N retranslocation accumulator
            annmax_retransn_patch[p]  = tempmax_retransn_patch[p]
            tempmax_retransn_patch[p] = zero(T)

            # update annual average 2m air temperature accumulator
            annavg_t2m_patch[p]  = tempavg_t2m_patch[p]
            tempavg_t2m_patch[p] = zero(T)

            # update annual NPP accumulator, convert to annual total
            annsum_npp_patch[p]  = tempsum_npp_patch[p] * dt
            tempsum_npp_patch[p] = zero(T)

            # update annual litfall accumulator, convert to annual total
            annsum_litfall_patch[p]  = tempsum_litfall_patch[p] * dt
            tempsum_litfall_patch[p] = zero(T)
        end
    end
end

function cn_annual_update!(mask_bgc_soilc::AbstractVector{Bool},
                           mask_bgc_vegp::AbstractVector{Bool},
                           bounds_c::UnitRange{Int},
                           bounds_p::UnitRange{Int},
                           col::ColumnData,
                           patch::PatchData,
                           cnveg_state::CNVegStateData,
                           cnveg_carbonflux::CNVegCarbonFluxData;
                           dt::Real,
                           days_per_year::Real,
                           secspday::Real = SECSPDAY)

    secspyear = days_per_year * secspday

    # --- Build end_of_year flag per column ---
    # Vector{Bool} over all columns; only BGC soil, non-FATES columns inside
    # bounds_c can trigger end-of-year. Loop1 (per-column kernel) writes it;
    # loop2 (per-patch kernel) reads end_of_year[col[p]].
    end_of_year = fill(false, length(mask_bgc_soilc))

    # Convert Float64 scalar args to the working precision so no double is
    # materialized inside the kernels on a Float32-only backend (Metal).
    FTe = eltype(cnveg_state.annsum_counter_col)
    _launch!(_cn_annual_eoy_col_kernel!, end_of_year,
             cnveg_state.annsum_counter_col,
             mask_bgc_soilc, col.is_fates,
             FTe(dt), FTe(secspyear), first(bounds_c), last(bounds_c))

    # --- Patch-level annual update ---
    FTp = eltype(cnveg_carbonflux.annsum_npp_patch)
    _launch!(_cn_annual_patch_kernel!,
             cnveg_state.annsum_potential_gpp_patch,
             cnveg_state.tempsum_potential_gpp_patch,
             cnveg_state.annmax_retransn_patch,
             cnveg_state.tempmax_retransn_patch,
             cnveg_state.annavg_t2m_patch,
             cnveg_state.tempavg_t2m_patch,
             cnveg_carbonflux.annsum_npp_patch,
             cnveg_carbonflux.tempsum_npp_patch,
             cnveg_carbonflux.annsum_litfall_patch,
             cnveg_carbonflux.tempsum_litfall_patch,
             mask_bgc_vegp, patch.column, end_of_year, col.is_fates,
             FTp(dt), first(bounds_p), last(bounds_p))

    # any_vegp gates the p2c averaging: true iff some in-range, non-FATES BGC
    # veg patch sits on an end-of-year column. Resolved on the host.
    any_vegp = false
    for p in bounds_p
        mask_bgc_vegp[p] || continue
        c = patch.column[p]
        if end_of_year[c] && !col.is_fates[c]
            any_vegp = true
            break
        end
    end

    # --- Patch-to-column averaging for end-of-year columns ---
    if any_vegp
        # p2c_1d_filter! expects a BitVector mask; end_of_year is a Vector{Bool}
        # (kernel-writable). Convert once for the (host-side) averaging helper.
        eoy_mask = BitVector(end_of_year)
        p2c_1d_filter!(cnveg_carbonflux.annsum_npp_col,
                        cnveg_carbonflux.annsum_npp_patch,
                        eoy_mask, col, patch)

        p2c_1d_filter!(cnveg_state.annavg_t2m_col,
                        cnveg_state.annavg_t2m_patch,
                        eoy_mask, col, patch)
    end

    return nothing
end
