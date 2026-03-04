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
function cn_annual_update!(mask_bgc_soilc::BitVector,
                           mask_bgc_vegp::BitVector,
                           bounds_c::UnitRange{Int},
                           bounds_p::UnitRange{Int},
                           col::ColumnData,
                           patch::PatchData,
                           cnveg_state::CNVegStateData,
                           cnveg_carbonflux::CNVegCarbonFluxData;
                           dt::Float64,
                           days_per_year::Float64,
                           secspday::Float64 = SECSPDAY)

    secspyear = days_per_year * secspday

    # --- Build end_of_year flag per column ---
    # We use a BitVector over bounds_c; only BGC soil, non-FATES columns
    # can trigger end-of-year.
    end_of_year = falses(length(mask_bgc_soilc))

    for c in bounds_c
        mask_bgc_soilc[c] || continue
        col.is_fates[c] && continue

        cnveg_state.annsum_counter_col[c] += dt

        if cnveg_state.annsum_counter_col[c] >= secspyear
            end_of_year[c] = true
            cnveg_state.annsum_counter_col[c] = 0.0
        end
    end

    # --- Patch-level annual update ---
    any_vegp = false
    for p in bounds_p
        mask_bgc_vegp[p] || continue
        c = patch.column[p]

        if end_of_year[c] && !col.is_fates[c]
            any_vegp = true

            # update annual plant ndemand accumulator
            cnveg_state.annsum_potential_gpp_patch[p]  = cnveg_state.tempsum_potential_gpp_patch[p]
            cnveg_state.tempsum_potential_gpp_patch[p] = 0.0

            # update annual total N retranslocation accumulator
            cnveg_state.annmax_retransn_patch[p]  = cnveg_state.tempmax_retransn_patch[p]
            cnveg_state.tempmax_retransn_patch[p] = 0.0

            # update annual average 2m air temperature accumulator
            cnveg_state.annavg_t2m_patch[p]  = cnveg_state.tempavg_t2m_patch[p]
            cnveg_state.tempavg_t2m_patch[p] = 0.0

            # update annual NPP accumulator, convert to annual total
            cnveg_carbonflux.annsum_npp_patch[p] = cnveg_carbonflux.tempsum_npp_patch[p] * dt
            cnveg_carbonflux.tempsum_npp_patch[p] = 0.0

            # update annual litfall accumulator, convert to annual total
            cnveg_carbonflux.annsum_litfall_patch[p] = cnveg_carbonflux.tempsum_litfall_patch[p] * dt
            cnveg_carbonflux.tempsum_litfall_patch[p] = 0.0
        end
    end

    # --- Patch-to-column averaging for end-of-year columns ---
    if any_vegp
        p2c_1d_filter!(cnveg_carbonflux.annsum_npp_col,
                        cnveg_carbonflux.annsum_npp_patch,
                        end_of_year, col, patch)

        p2c_1d_filter!(cnveg_state.annavg_t2m_col,
                        cnveg_state.annavg_t2m_patch,
                        end_of_year, col, patch)
    end

    return nothing
end
