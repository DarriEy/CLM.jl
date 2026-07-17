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

# ==========================================================================
# Fire-bundle initialization (Fortran: CNFireFactory's `create_cnfire_method`
# followed by `fire_method%FireInit(bounds, NLFilename)` = BaseFireInit ->
# InitAllocate + surfdataread; FireDataBaseType.F90 / CNFireBaseMod.F90).
#
# Lives here (not in fire_base.jl) because the signatures reference
# CNFireLi2014Data / PftConFireLi2014, which are defined in fire_li2014.jl —
# included after fire_base.jl.
# ==========================================================================

"""
    cnfire_init_allocate!(fire_data, fire_li2014, dgvs_fire, nc, np, ng)

Allocate the fire state arrays: `btran2_patch` (patch), `forc_hdm`/`forc_lnfm`
(gridcell), `gdp_lf_col`/`peatf_lf_col`/`abm_lf_col` (column), plus the DGVS
`nind_patch` the fire mortality decrements.

Ported from `BaseFireInit` (`InitAllocate`) in `CNFireBaseMod.F90` /
`FireDataBaseType.F90`. Fortran allocates to NaN; we allocate to zero and rely on
`cnfire_surfdata_read!` + `fire_stream_interp!` to fill the inputs, so a missing
input degrades to "no ignition" rather than NaN-poisoning the live (and fatal)
C/N balance check. `abm_lf_col` starts at 13 = "no crop-fire month" (the crop-fire
kernel fires only when `kmo == abm_lf[c]`, and month 13 never matches).

`dgvs_fire.nind_patch` is left alone when already non-empty, so a use_cndv run can
alias the live `DGVSData.nind_patch` array into it and see the fire mortality.
"""
function cnfire_init_allocate!(fire_data::CNFireBaseData,
                               fire_li2014::CNFireLi2014Data,
                               dgvs_fire::DgvsFireData,
                               nc::Int, np::Int, ng::Int)
    fire_data.btran2_patch   = zeros(np)
    fire_li2014.forc_hdm     = zeros(ng)
    fire_li2014.forc_lnfm    = zeros(ng)
    fire_li2014.gdp_lf_col   = zeros(nc)
    fire_li2014.peatf_lf_col = zeros(nc)
    fire_li2014.abm_lf_col   = fill(13, nc)
    isempty(dgvs_fire.nind_patch) && (dgvs_fire.nind_patch = zeros(np))
    return nothing
end

"""
    cnfire_surfdata_read!(fire_li2014, surf, col, bounds_col)

Scatter the surface-dataset fire inputs (`gdp`, `peatf`, `abm`) from gridcell to
column: `gdp_lf_col[c] = gdp[col.gridcell[c]]`, and likewise for peatf/abm.

Ported from `surfdataread` in `cpl/share_esmf/FireDataBaseType.F90`, including its
hard error when any of the three is missing from `fsurdat`. No unit conversion —
Fortran applies none.
"""
function cnfire_surfdata_read!(fire_li2014::CNFireLi2014Data,
                               surf, col::ColumnData,
                               bounds_col::UnitRange{Int})
    isempty(surf.gdp) &&
        error("cnfire_surfdata_read!: gdp NOT on surfdata file (required by the Li fire model)")
    isempty(surf.peatf) &&
        error("cnfire_surfdata_read!: peatf NOT on surfdata file (required by the Li fire model)")
    isempty(surf.abm) &&
        error("cnfire_surfdata_read!: abm NOT on surfdata file (required by the Li fire model)")

    for c in bounds_col
        g = col.gridcell[c]
        (1 <= g <= length(surf.gdp)) || continue
        fire_li2014.gdp_lf_col[c]   = surf.gdp[g]
        fire_li2014.peatf_lf_col[c] = surf.peatf[g]
        fire_li2014.abm_lf_col[c]   = surf.abm[g]
    end
    return nothing
end

"""
    cnfire_pftcon_init!(pftcon_fire, pftcon_fire_li2014, p)

Fill the fire PFT-parameter holders from the global `pftcon` `p` (populated by
`readParameters!`). These are exactly the `pftcon` members the Fortran fire
modules `associate` to; the holders exist only because the Julia fire kernels take
a narrow, GPU-adaptable parameter struct rather than the whole pftcon.

Arrays are ALIASED, not copied, so a later pftcon override is picked up.
"""
function cnfire_pftcon_init!(pftcon_fire::PftConFireBase,
                             pftcon_fire_li2014::PftConFireLi2014,
                             p)
    pftcon_fire.woody    = p.woody
    pftcon_fire.cc_leaf  = p.cc_leaf
    pftcon_fire.cc_lstem = p.cc_lstem
    pftcon_fire.cc_dstem = p.cc_dstem
    pftcon_fire.cc_other = p.cc_other
    pftcon_fire.fm_leaf  = p.fm_leaf
    pftcon_fire.fm_lstem = p.fm_lstem
    pftcon_fire.fm_other = p.fm_other
    pftcon_fire.fm_root  = p.fm_root
    pftcon_fire.fm_lroot = p.fm_lroot
    pftcon_fire.fm_droot = p.fm_droot
    pftcon_fire.lf_f     = p.lf_f
    pftcon_fire.fr_f     = p.fr_f
    pftcon_fire.smpso    = p.smpso
    pftcon_fire.smpsc    = p.smpsc
    pftcon_fire.rswf_min = p.rswf_min
    pftcon_fire.rswf_max = p.rswf_max

    pftcon_fire_li2014.fsr_pft = p.fsr_pft
    pftcon_fire_li2014.fd_pft  = p.fd_pft
    return nothing
end
