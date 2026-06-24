# ==========================================================================
# Ported from: src/biogeochem/CNFireFactoryMod.F90
# Factory / dispatcher selecting the active CN fire method from the
# `cnfire_method` namelist string.
#
# The Fortran factory allocates a polymorphic fire_method_type whose virtual
# CNFireArea / CNFireFluxes / need_lightning_and_popdens dispatch to the chosen
# Li-year or NoFire implementation. Julia has no class hierarchy here, so we
# map the namelist string to a method symbol and dispatch with a small
# `cnfire_area!` / `cnfire_fluxes!` / `need_lightning_and_popdens` shim.
#
# Recognized cnfire_method strings (mirroring create_cnfire_method):
#   "nofire"                     -> :nofire
#   "li2014qianfrc"   (default)  -> :li2014
#   "li2016crufrc"               -> :li2016
#   "li2021gswpfrc"              -> :li2021
#   "li2024gswpfrc" / "li2024crujra" -> :li2024
# ==========================================================================

"""
    cnfire_method_symbol(name::AbstractString) -> Symbol

Map a `cnfire_method` namelist string to the internal method symbol.
Ported from the `select case (trim(fire_method))` in `create_cnfire_method`.
Throws on an unknown method (mirrors the Fortran `endrun`).
"""
function cnfire_method_symbol(name::AbstractString)
    s = lowercase(strip(name))
    if s == "nofire"
        return :nofire
    elseif s == "li2014qianfrc"
        return :li2014
    elseif s == "li2016crufrc"
        return :li2016
    elseif s == "li2021gswpfrc"
        return :li2021
    elseif s == "li2024gswpfrc" || s == "li2024crujra"
        return :li2024
    else
        error("create_cnfire_method: unknown cnfire_method: $(name)")
    end
end

"""
    need_lightning_and_popdens(method::Symbol) -> Bool

Whether the selected fire method needs lightning + population-density inputs.
All Li-family methods do; NoFire does not.
"""
function need_lightning_and_popdens(method::Symbol)
    if method === :nofire
        return need_lightning_and_popdens_nofire()
    elseif method === :li2014
        return need_lightning_and_popdens_li2014()
    elseif method === :li2016
        return need_lightning_and_popdens_li2016()
    elseif method === :li2021
        return need_lightning_and_popdens_li2021()
    elseif method === :li2024
        return need_lightning_and_popdens_li2024()
    else
        error("need_lightning_and_popdens: unknown method $(method)")
    end
end

"""
    cnfire_area!(method, fire_li2014, pftcon_li2014, fire_data, cnfire_const,
                 cnfire_params, pftcon, masks..., bounds..., patch, col, grc,
                 soilstate, h2osoi_vol_col, cnveg_state, cnveg_cs,
                 decomp_cascade_con, totlitc_col, decomp_cpools_vr_col,
                 t_soi17cm_col; kwargs...)

Dispatch column-level burned-area computation to the fire method selected by
`method` (a `Symbol` from `cnfire_method_symbol`). All Li methods share the same
positional/keyword interface as `cnfire_area_li2014!`; the dispatcher forwards
the common arguments and lets each method pick the keyword args it needs (extra
kwargs are ignored by methods that don't use them, e.g. `prec30_patch` for
Li2024, `rh30_patch` for Li2016/21/24). NoFire ignores the fire inputs and just
zeroes burned-area state.

Ported from the dispatch on `cnfire_method%CNFireArea(...)`.
"""
function cnfire_area!(
    method::Symbol,
    fire_li2014::CNFireLi2014Data,
    pftcon_li2014::PftConFireLi2014,
    fire_data::CNFireBaseData,
    cnfire_const::CNFireConstData,
    cnfire_params::CNFireParams,
    pftcon::PftConFireBase,
    mask_soilc::AbstractVector{Bool},
    mask_soilp::AbstractVector{Bool},
    mask_exposedveg::AbstractVector{Bool},
    mask_noexposedveg::AbstractVector{Bool},
    bounds_c::UnitRange{Int},
    bounds_p::UnitRange{Int},
    patch::PatchData,
    col::ColumnData,
    grc::GridcellData,
    soilstate::SoilStateData,
    h2osoi_vol_col::AbstractMatrix{<:Real},
    cnveg_state::CNVegStateData,
    cnveg_cs::CNVegCarbonStateData,
    decomp_cascade_con::DecompCascadeConData,
    totlitc_col::AbstractVector{<:Real},
    decomp_cpools_vr_col::AbstractArray{<:Real,3},
    t_soi17cm_col::AbstractVector{<:Real};
    kwargs...
)
    if method === :nofire
        return cnfire_area_nofire!(cnveg_state, cnveg_cs, patch,
                                   mask_soilc, bounds_c, bounds_p)
    end

    common = (fire_li2014, pftcon_li2014, fire_data, cnfire_const,
              cnfire_params, pftcon, mask_soilc, mask_soilp, mask_exposedveg,
              mask_noexposedveg, bounds_c, bounds_p, patch, col, grc,
              soilstate, h2osoi_vol_col, cnveg_state, cnveg_cs,
              decomp_cascade_con, totlitc_col, decomp_cpools_vr_col,
              t_soi17cm_col)

    if method === :li2014
        # li2014 is the default path and keeps a fixed keyword signature (no
        # `_ignored...` catch-all), so drop keys it doesn't accept (e.g. the
        # Li2024-only prec30_patch/croplive_patch) to preserve byte-identical
        # default behavior.
        li2014_keys = (:forc_rh_grc, :forc_wind_grc, :forc_t_col, :forc_rain_col,
                       :forc_snow_col, :prec60_patch, :prec10_patch, :fsat_col,
                       :wf_col, :wf2_col, :dt, :dayspyr, :kmo, :kda, :mcsec,
                       :nstep, :nlevgrnd, :nlevdecomp, :ndecomp_pools, :nc4_grass,
                       :nc3crop, :ndllf_evr_tmp_tree, :nbrdlf_evr_trp_tree,
                       :nbrdlf_dcd_trp_tree, :nbrdlf_evr_shrub, :noveg,
                       :transient_landcover, :soil_suction_fn)
        kw14 = NamedTuple(k => v for (k, v) in pairs(kwargs) if k in li2014_keys)
        return cnfire_area_li2014!(common...; kw14...)
    elseif method === :li2016
        return cnfire_area_li2016!(common...; kwargs...)
    elseif method === :li2021
        return cnfire_area_li2021!(common...; kwargs...)
    elseif method === :li2024
        return cnfire_area_li2024!(common...; kwargs...)
    else
        error("cnfire_area!: unknown method $(method)")
    end
end

"""
    cnfire_fluxes!(method, masks..., bounds..., cnfire_const, pftcon, ...; kwargs...)

Dispatch fire C/N flux computation by `method`. Every Li-family method inherits
the base `CNFireFluxes` with litter/CWD combustion-completeness factors 0.5/0.25,
so they all forward to `cnfire_fluxes_li2014!`. NoFire forwards too, but with
zero burned area the fluxes are all zero.

NOTE: this dispatcher is `cnfire_fluxes_dispatch!` to avoid clashing with the
base `cnfire_fluxes!` (the worker routine in `fire_base.jl`).

Ported from the dispatch on `cnfire_method%CNFireFluxes(...)`.
"""
function cnfire_fluxes_dispatch!(method::Symbol, args...; kwargs...)
    if method === :nofire
        return cnfire_fluxes_nofire!(args...; kwargs...)
    elseif method === :li2014
        return cnfire_fluxes_li2014!(args...; kwargs...)
    elseif method === :li2016
        return cnfire_fluxes_li2016!(args...; kwargs...)
    elseif method === :li2021
        return cnfire_fluxes_li2021!(args...; kwargs...)
    elseif method === :li2024
        return cnfire_fluxes_li2024!(args...; kwargs...)
    else
        error("cnfire_fluxes_dispatch!: unknown method $(method)")
    end
end
