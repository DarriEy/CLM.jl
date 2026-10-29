# ==========================================================================
# Ported from: src/main/lnd2atmMod.F90 (ice-runoff handling + the lnd->rof block)
#
# Public functions:
#   handle_ice_runoff!           — total column ice runoff, and its optional
#                                  conversion to liquid runoff + a compensating
#                                  negative sensible heat flux
#   add_liq_from_ice_to_runoff!  — route the converted liquid into qflx_qrgwl /
#                                  qflx_runoff (so the water balance still closes)
#   lnd2atm_rof!                 — the `! lnd -> rof` section of `lnd2atm`: every
#                                  column runoff term aggregated to the gridcell
#                                  fields the river model is forced with
#
# The patch->gridcell ENERGY/ALBEDO half of `lnd2atm` lives in
# `infrastructure/lnd2atm_mod.jl`. The column-level `qflx_ice_runoff_col` computed
# here is the term `BalanceCheck` subtracts in the column water balance
# (BalanceCheckMod.F90:650), so it is required for the water-balance check to be
# meaningful.
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

# =========================================================================
# lnd2atm_rof! — the `! lnd -> rof` block of `lnd2atm`
# =========================================================================

"""
    lnd2atm_rof!(l2a, waterflux, col, lun, grc, bounds_c, bounds_l, bounds_g;
                 use_hillslope_routing = false,
                 is_hillslope_column   = nothing)

Aggregate every column-level runoff term to the gridcell-level fields the river
model is forced with, and store them on `l2a`:

| gridcell field                  | source column field         | note |
|---------------------------------|-----------------------------|------|
| `qflx_rofliq_qsur_grc`          | `qflx_surf_col`             | surface runoff |
| `qflx_rofliq_qsub_grc`          | `qflx_drain_col`            | subsurface drainage |
| `qflx_rofliq_drain_perched_grc` | `qflx_drain_perched_col`    | perched water table |
| `qflx_rofliq_qgwl_grc`          | `qflx_qrgwl_col`            | minus `qflx_liq_dynbal_grc` |
| `qflx_rofliq_grc`               | `qflx_runoff_col`           | minus `qflx_liq_dynbal_grc` |
| `qflx_rofice_grc`               | `l2a.qflx_ice_runoff_col`   | minus `qflx_ice_dynbal_grc` |
| `qirrig_grc`                    | `qflx_sfc_irrig_col`        | irrigation withdrawal from rof |
| `qflx_evap_tot_grc`             | `qflx_evap_tot_col`         | total evapotranspiration |
| `qflx_rofliq_stream_grc`        | `volumetric_streamflow_lun` | hillslope routing only |

Every aggregation uses `c2l_scale_type='urbanf'` (walls x 3*canyon_hwr, roads x 3)
and `l2g_scale_type='unity'`, exactly as `lnd2atmMod.F90` and `BalanceCheckMod.F90`
do — those fluxes are per m2 of WALL / ROAD area, not of ground area, so `unity`
would make an urban gridcell fail to close even when every column closes to 1e-13.

Under `use_hillslope_routing`, hillslope columns send their surface / drainage /
perched runoff to the stream channel rather than straight to rof, so they are
EXCLUDED from the qsur/qsub/drain_perched gridcell averages to avoid
double-counting, and the stream discharge is summed (volume/time, NOT area
weighted) over the gridcell's active landunits into `qflx_rofliq_stream_grc`.
`is_hillslope_column` defaults to `col.is_hillslope_column`.

Must be called AFTER [`handle_ice_runoff!`](@ref) and
[`add_liq_from_ice_to_runoff!`](@ref) — CTSM orders it the same way, so that the
melted-ice liquid is already inside `qflx_qrgwl_col` / `qflx_runoff_col` when they
are aggregated.

Ported from the `! lnd -> rof` section of `lnd2atm` in `lnd2atmMod.F90`
(lines 340-470), writing the `waterlnd2atm_type` members of
`Waterlnd2atmType.F90`.
"""
function lnd2atm_rof!(
    l2a::Lnd2AtmData,
    waterflux::Union{WaterFluxData, WaterFluxBulkData},
    col::ColumnData,
    lun::LandunitData,
    grc::GridcellData,
    bounds_c::UnitRange{Int},
    bounds_l::UnitRange{Int},
    bounds_g::UnitRange{Int};
    use_hillslope_routing::Bool = false,
    is_hillslope_column::Union{Nothing, AbstractVector{Bool}} = nothing
)
    wf = waterflux isa WaterFluxBulkData ? waterflux.wf : waterflux
    isempty(bounds_g) && return nothing

    # --- hillslope stream channel (lnd2atmMod.F90:343-354) ---
    if use_hillslope_routing && !isempty(wf.volumetric_streamflow_lun)
        hillslope_streamflow_to_grc!(l2a.qflx_rofliq_stream_grc,
            wf.volumetric_streamflow_lun, lun, grc, bounds_l, bounds_g)
    else
        for g in bounds_g
            l2a.qflx_rofliq_stream_grc[g] = zero(eltype(l2a.qflx_rofliq_stream_grc))
        end
    end

    # --- surface / subsurface / perched runoff (lnd2atmMod.F90:356-415) ---
    if use_hillslope_routing
        # Hillslope runoff is sent to the stream, not to rof, and is already
        # accounted for in qflx_rofliq_stream_grc: zero those columns before
        # averaging. On the default path the mask is all-false, so the masked
        # copies equal the originals and this is byte-identical to the plain c2g.
        hill = is_hillslope_column === nothing ? col.is_hillslope_column : is_hillslope_column
        nothill = (!).(hill)
        c2g_urbanf!(l2a.qflx_rofliq_qsur_grc, wf.qflx_surf_col .* nothill,
                    col, lun, bounds_c, bounds_g)
        c2g_urbanf!(l2a.qflx_rofliq_qsub_grc, wf.qflx_drain_col .* nothill,
                    col, lun, bounds_c, bounds_g)
        c2g_urbanf!(l2a.qflx_rofliq_drain_perched_grc,
                    wf.qflx_drain_perched_col .* nothill, col, lun, bounds_c, bounds_g)
    else
        c2g_urbanf!(l2a.qflx_rofliq_qsur_grc, wf.qflx_surf_col,
                    col, lun, bounds_c, bounds_g)
        c2g_urbanf!(l2a.qflx_rofliq_qsub_grc, wf.qflx_drain_col,
                    col, lun, bounds_c, bounds_g)
        c2g_urbanf!(l2a.qflx_rofliq_drain_perched_grc, wf.qflx_drain_perched_col,
                    col, lun, bounds_c, bounds_g)
    end

    # --- qgwl / total liquid runoff, minus the dynbal correction ---
    # (lnd2atmMod.F90:437-453). qflx_liq_dynbal_grc is the water the dynamic
    # landunit-area adjustment took out of / put into storage; it is NOT real
    # runoff, so it is removed from what the river receives.
    #
    # NOTE on the isfinite guard: `waterflux_init!` allocates the two dynbal
    # gridcell fluxes as NaN, and nothing on the current CLM.jl path writes them
    # (the dynamic-landunit machinery keeps its own copy on `DynBalData`). CTSM
    # sets them every step, so there they are always finite. Subtracting an unset
    # NaN would silently poison the exported river forcing AND the gridcell
    # balance check, so a non-finite value is treated as "no dynbal correction",
    # which is what an unwritten field means.
    c2g_urbanf!(l2a.qflx_rofliq_qgwl_grc, wf.qflx_qrgwl_col, col, lun, bounds_c, bounds_g)
    c2g_urbanf!(l2a.qflx_rofliq_grc, wf.qflx_runoff_col, col, lun, bounds_c, bounds_g)
    if length(wf.qflx_liq_dynbal_grc) >= last(bounds_g)
        for g in bounds_g
            d = wf.qflx_liq_dynbal_grc[g]
            isfinite(d) || continue
            l2a.qflx_rofliq_qgwl_grc[g] -= d
            l2a.qflx_rofliq_grc[g]      -= d
        end
    end

    # --- irrigation withdrawal (lnd2atmMod.F90:455-458) ---
    c2g_urbanf!(l2a.qirrig_grc, wf.qflx_sfc_irrig_col, col, lun, bounds_c, bounds_g)

    # --- ice runoff, minus the dynbal correction (lnd2atmMod.F90:460-468) ---
    c2g_urbanf!(l2a.qflx_rofice_grc, l2a.qflx_ice_runoff_col, col, lun, bounds_c, bounds_g)
    if length(wf.qflx_ice_dynbal_grc) >= last(bounds_g)
        for g in bounds_g
            d = wf.qflx_ice_dynbal_grc[g]
            isfinite(d) && (l2a.qflx_rofice_grc[g] -= d)
        end
    end

    # --- total evapotranspiration (waterlnd2atm_type%qflx_evap_tot_grc) ---
    c2g_urbanf!(l2a.qflx_evap_tot_grc, wf.qflx_evap_tot_col, col, lun, bounds_c, bounds_g)

    return nothing
end
