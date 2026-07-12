# ==========================================================================
# Ported from: src/main/lnd2atmMod.F90 (ice-runoff handling)
#
# Public functions:
#   handle_ice_runoff!           — total column ice runoff, and its optional
#                                  conversion to liquid runoff + a compensating
#                                  negative sensible heat flux
#   add_liq_from_ice_to_runoff!  — route the converted liquid into qflx_qrgwl /
#                                  qflx_runoff (so the water balance still closes)
#
# Only the ice-runoff part of `lnd2atm` is ported here: the rest of lnd2atm is
# gridcell averaging for the atmosphere/river coupler, which the standalone port
# does not need. The column-level `qflx_ice_runoff_col` computed here is the term
# `BalanceCheck` subtracts in the column water balance (BalanceCheckMod.F90:650),
# so it is required for the water-balance check to be meaningful.
# ==========================================================================

# ---- handle_ice_runoff! : per-column, fully independent (one thread per column) ----
# Mask is `col_active` (Fortran loops all columns and guards on col%active). Consts
# are eltype-converted so no Float64 reaches a Float32-only backend (Metal); on a
# Float64 CPU run every store is byte-identical to the plain loop.
@kernel function _l2a_handle_ice_runoff_kernel!(qflx_ice_runoff_col, qflx_liq_from_ice_col,
        eflx_sh_ice_to_liq_col,
        @Const(qflx_ice_runoff_snwcp_col), @Const(qflx_ice_runoff_xs_col),
        @Const(col_active), @Const(col_landunit), @Const(col_gridcell),
        @Const(lun_itype), @Const(ice_runoff_melted_grc),
        melt_non_icesheet_ice_runoff::Bool, istice::Int, hfus, cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && col_active[c]
        T = eltype(qflx_ice_runoff_col)

        qflx_ice_runoff   = qflx_ice_runoff_snwcp_col[c] + qflx_ice_runoff_xs_col[c]
        qflx_liq_from_ice = zero(T)
        eflx_sh_ice_to_liq = zero(T)

        if melt_non_icesheet_ice_runoff
            l = col_landunit[c]
            g = col_gridcell[c]
            # Non-icesheet landunits always convert; icesheet columns convert only
            # where the GLC behavior says their ice runoff is melted en route.
            do_conversion = (lun_itype[l] != istice) || ice_runoff_melted_grc[g]
            if do_conversion
                # Ice -> liquid absorbs energy, so it is a negative heat flux to atm.
                # qflx_ice_runoff_col is mm H2O/s == kg m-2 s-1, so multiply by hfus.
                eflx_sh_ice_to_liq = -qflx_ice_runoff * T(hfus)
                qflx_liq_from_ice  = qflx_ice_runoff
                qflx_ice_runoff    = zero(T)
            end
        end

        qflx_ice_runoff_col[c]    = qflx_ice_runoff
        qflx_liq_from_ice_col[c]  = qflx_liq_from_ice
        eflx_sh_ice_to_liq_col[c] = eflx_sh_ice_to_liq
    end
end

"""
    handle_ice_runoff!(l2a, waterflux, col_data, lun_data, bounds_c;
                       ice_runoff_melted_grc = nothing)

Compute the total column-level ice runoff

    qflx_ice_runoff_col = qflx_ice_runoff_snwcp_col + qflx_ice_runoff_xs_col

(solid runoff from snow capping + solid runoff from excess soil ice) and divide it
between (a) ice runoff and (b) liquid runoff with a compensating negative sensible
heat flux, controlled by `l2a.params.melt_non_icesheet_ice_runoff`.

Ice runoff is a crude parameterization of iceberg calving, which is only really
appropriate where an ice sheet terminates at the land-ocean boundary. When
`melt_non_icesheet_ice_runoff` is true, ice runoff from non-icesheet columns (and
from icesheet columns in gridcells flagged by `ice_runoff_melted_grc`) is instead
melted: it becomes `qflx_liq_from_ice_col` with `eflx_sh_ice_to_liq_col` absorbing
the latent heat of fusion. The Fortran namelist default is `.false.`.

`ice_runoff_melted_grc` mirrors `glc_behavior%ice_runoff_melted_grc`; glc2lnd
coupling is not ported, so it defaults to all-false (standalone CLM).

Results are written to `l2a.qflx_ice_runoff_col`, `l2a.qflx_liq_from_ice_col` and
`l2a.eflx_sh_ice_to_liq_col`.

Ported from `handle_ice_runoff` in `lnd2atmMod.F90`.
"""
function handle_ice_runoff!(
    l2a::Lnd2AtmData,
    waterflux::Union{WaterFluxData, WaterFluxBulkData},
    col_data::ColumnData,
    lun_data::LandunitData,
    bounds_c::UnitRange{Int};
    ice_runoff_melted_grc::Union{Nothing, AbstractVector{Bool}} = nothing
)
    wf = waterflux isa WaterFluxBulkData ? waterflux.wf : waterflux
    isempty(bounds_c) && return nothing

    # glc2lnd coupling is not ported -> no gridcell has its ice runoff melted en
    # route (same convention as glc_dyn_runoff_routing_grc in hydrology_drainage!).
    ng = maximum(Array(col_data.gridcell)[bounds_c])
    melted = ice_runoff_melted_grc === nothing ?
        fill!(similar(col_data.active, Bool, ng), false) : ice_runoff_melted_grc

    _launch!(_l2a_handle_ice_runoff_kernel!,
        l2a.qflx_ice_runoff_col, l2a.qflx_liq_from_ice_col, l2a.eflx_sh_ice_to_liq_col,
        wf.qflx_ice_runoff_snwcp_col, wf.qflx_ice_runoff_xs_col,
        col_data.active, col_data.landunit, col_data.gridcell, lun_data.itype,
        melted, l2a.params.melt_non_icesheet_ice_runoff, ISTICE,
        eltype(l2a.qflx_ice_runoff_col)(HFUS), first(bounds_c), last(bounds_c);
        ndrange = length(l2a.qflx_ice_runoff_col))

    return nothing
end

# ---- add_liq_from_ice_to_runoff! : per-column, fully independent ----
@kernel function _l2a_liq_from_ice_runoff_kernel!(qflx_qrgwl_col, qflx_runoff_col,
        @Const(qflx_liq_from_ice_col), @Const(col_active), cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && col_active[c]
        qflx_qrgwl_col[c]  += qflx_liq_from_ice_col[c]
        qflx_runoff_col[c] += qflx_liq_from_ice_col[c]
    end
end

"""
    add_liq_from_ice_to_runoff!(waterflux, l2a, col_data, bounds_c)

Route `qflx_liq_from_ice_col` (ice runoff that `handle_ice_runoff!` melted) into
`qflx_qrgwl_col`, and analogously into `qflx_runoff_col` (of which qflx_qrgwl is a
component). Without this the melted water would leave the column unaccounted for
in the water balance. A no-op when `melt_non_icesheet_ice_runoff` is false (the
default), since then `qflx_liq_from_ice_col` is identically zero.

Ported from the inline loop in `lnd2atm` in `lnd2atmMod.F90`.
"""
function add_liq_from_ice_to_runoff!(
    waterflux::Union{WaterFluxData, WaterFluxBulkData},
    l2a::Lnd2AtmData,
    col_data::ColumnData,
    bounds_c::UnitRange{Int}
)
    wf = waterflux isa WaterFluxBulkData ? waterflux.wf : waterflux
    isempty(bounds_c) && return nothing

    _launch!(_l2a_liq_from_ice_runoff_kernel!,
        wf.qflx_qrgwl_col, wf.qflx_runoff_col, l2a.qflx_liq_from_ice_col,
        col_data.active, first(bounds_c), last(bounds_c);
        ndrange = length(wf.qflx_qrgwl_col))

    return nothing
end
