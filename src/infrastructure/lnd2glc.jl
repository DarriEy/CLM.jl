# ==========================================================================
# Ported from: src/main/lnd2glcMod.F90
# Land-to-ice-sheet (CLM -> CISM) coupling fields.
#
# Public functions:
#   lnd2glc_init!            — allocate the (gridcell x elevation-class) arrays
#   update_lnd2glc!          — fill tsrf / topo / qice from CLM's own state
#   bareland_normalization   — conservation factor for bare-land fluxes
#
# The fields sent from land to glc via the coupler are the 's2x' (sno-to-coupler)
# bundle. 'Sno' is a misnomer: the exchanged data describe the ICE beneath the
# snow, not the snow itself — CESM reserves 'ice' for sea ice.
#
# Index convention. Fortran declares the arrays as `(begg:endg, 0:maxpatch_glc)`:
# the second dimension is the glacier ELEVATION CLASS, with 0 meaning "the bare
# land (istsoil) portion of the gridcell" and 1..maxpatch_glc the glacier_mec
# classes. Julia cannot use a 0-based dimension, so the arrays here are
# `(maxpatch_glc + 1, ng)` — CLASS-MAJOR, so a whole gridcell's elevation profile
# is contiguous, which is the layout every coupler wants — and Fortran class `n`
# is stored at Julia row `n + 1`. `lnd2glc_class_index(n)` performs that map.
#
# CLM.jl DOES carry a genuine glacier elevation-class subgrid: `set_landunit_ice!`
# (init_gridcells.jl) creates one ISTICE column per non-zero `wt_glc_mec[g, m]`
# with `col.itype = ice_class_to_col_itype(m)`, and `topo_init_cold!` gives each
# of them its own `topo_col` from the surface dataset's `topo_glc_mec[g, m]`. So
# this is a real port over real per-elevation-class columns, not a stand-in.
# ==========================================================================

# --------------------------------------------------------------------------
# Data type
# --------------------------------------------------------------------------

"""
    Lnd2GlcData

Land-to-glacier (ice sheet) coupling fields, one value per gridcell per glacier
elevation class. Class 0 (Julia row 1) carries the bare-land information.

- `tsrf_grc` : surface temperature of the ice/bare land (K) — `t_soisno(c,1)`
- `topo_grc` : surface elevation (m) — `topo_col(c)`
- `qice_grc` : surface mass balance (mm H2O/s, i.e. kg m-2 s-1), positive =
  ice gain — `qflx_glcice_col(c)` times the bare-land normalization

Ported from `lnd2glc_type` in `lnd2glcMod.F90`.
"""
Base.@kwdef mutable struct Lnd2GlcData{FT<:Real, M<:AbstractMatrix{FT}}
    # (maxpatch_glc+1, ng) — see the index convention note at the top of the file.
    tsrf_grc::M = Matrix{Float64}(undef, 0, 0)
    topo_grc::M = Matrix{Float64}(undef, 0, 0)
    qice_grc::M = Matrix{Float64}(undef, 0, 0)
end

Lnd2GlcData{FT}(; kwargs...) where {FT<:Real} = Lnd2GlcData{FT, Matrix{FT}}(; kwargs...)
Adapt.@adapt_structure Lnd2GlcData

"""
    lnd2glc_class_index(n) -> Int

Map a Fortran glacier elevation class `n` (0 = bare land, 1..maxpatch_glc =
glacier_mec classes) to the row index of the `Lnd2GlcData` arrays.
"""
@inline lnd2glc_class_index(n::Integer) = Int(n) + 1

"""
    lnd2glc_init!(l2g, ng; maxpatch_glc = varpar.maxpatch_glc)

Allocate and zero the land->glc arrays as `(maxpatch_glc + 1, ng)`.

Ported from `InitAllocate` in `lnd2glcMod.F90`.
"""
function lnd2glc_init!(l2g::Lnd2GlcData{FT}, ng::Int;
                       maxpatch_glc::Int = varpar.maxpatch_glc) where {FT}
    nec1 = maxpatch_glc + 1
    l2g.tsrf_grc = zeros(FT, nec1, ng)
    l2g.topo_grc = zeros(FT, nec1, ng)
    l2g.qice_grc = zeros(FT, nec1, ng)
    return nothing
end

"""
    lnd2glc_clean!(l2g)

Reset all arrays to empty. Mirrors Fortran deallocation.
"""
function lnd2glc_clean!(l2g::Lnd2GlcData{FT}) where {FT}
    l2g.tsrf_grc = Matrix{FT}(undef, 0, 0)
    l2g.topo_grc = Matrix{FT}(undef, 0, 0)
    l2g.qice_grc = Matrix{FT}(undef, 0, 0)
    return nothing
end

# --------------------------------------------------------------------------
# bareland_normalization
# --------------------------------------------------------------------------

"""
    bareland_normalization(c, col, lun, grc) -> Float64

Normalization factor for fluxes coming from the bare-land portion of a gridcell.
Fluxes must be multiplied by this before being sent to the ice sheet.

The ice sheet sees only two covers: glaciated and bare. CLM subdivides the bare
part into several landunits but sends only the natural-vegetated landunit's value
(the flux is 0 in the others), so that value must be spread over the whole bare
fraction to conserve mass. For a gridcell that is 60% glacier_mec / 30% natural
veg / 10% lake, the ice sheet calls 40% of it "bare land", so a 1 m SMB computed
on the natural-veg landunit must be sent as `1 * 0.3/0.4 = 0.75` m.

If the gridcell is entirely glacier the factor is arbitrary and is set to 1 to
avoid dividing by zero.

Ported from `bareland_normalization` in `lnd2glcMod.F90`.
"""
function bareland_normalization(c::Int, col::ColumnData, lun::LandunitData,
                                grc::GridcellData)
    tol = 1.0e-13    # Fortran `tol` for subgrid weight equality
    g = col.gridcell[c]
    area_glacier = get_landunit_weight(grc, lun, g, ISTICE)
    if abs(area_glacier - 1.0) < tol
        # Whole gridcell is glacier: the normalization is arbitrary, do none.
        return 1.0
    else
        return col.wtgcell[c] / (1.0 - area_glacier)
    end
end

# --------------------------------------------------------------------------
# update_lnd2glc!
# --------------------------------------------------------------------------

"""
    update_lnd2glc!(l2g, temperature, waterflux, topo, col, lun, grc,
                    filter_do_smb_c, bounds_c, bounds_g; init = false,
                    use_hillslope = false)

Fill the land->glc coupling fields from CLM's own snow/energy/glacier state.

Defaults, applied to EVERY gridcell/class before the filter loop, and therefore
what is sent for points outside the `do_smb` filter:
- `qice = 0` — required for conservation (the ice sheet must not be handed mass
  CLM did not compute).
- `tsrf = TFRZ` — 0 degC; arbitrary but as reasonable as anything outside the filter.
- `topo = 0` — the historical CTSM choice.

Then, for every column `c` in the `do_smb` filter (`filter_do_smb_c` is a
`BitVector`/`Vector{Bool}` mask over columns, i.e. CLM.jl's `filt.do_smb_c`):

- ISTICE landunit → class `n = col_itype_to_ice_class(col.itype[c])`, flux
  normalization 1.
- ISTSOIL landunit → class `n = 0` (the bare-land slot), flux normalization
  [`bareland_normalization`](@ref).
- any other landunit → skipped: those types pass no information to the ice sheet.
  (For that to be acceptable CTSM requires virtual vegetated columns in any
  gridcell made solely of glacier plus another special landunit.)

`tsrf` (`t_soisno(c,1)`) and `topo` (`topo_col(c)`) are valid even during
initialization and are always set. `qice` (`qflx_glcice_col(c)`) is NOT valid
until the run loop, so with `init = true` it keeps its default of 0.

Assigning the same `(g, n)` slot twice is an error — it means the gridcell has
multiple columns in the istsoil landunit, which this routine cannot handle. CTSM
bypasses that check under `use_hillslope` (ESCOMP/CTSM#204); the same escape
hatch is available here.

Ported from `update_lnd2glc` in `lnd2glcMod.F90`.
"""
function update_lnd2glc!(
    l2g::Lnd2GlcData,
    temperature::TemperatureData,
    waterflux::Union{WaterFluxData, WaterFluxBulkData},
    topo::TopoData,
    col::ColumnData,
    lun::LandunitData,
    grc::GridcellData,
    filter_do_smb_c::AbstractVector{Bool},
    bounds_c::UnitRange{Int},
    bounds_g::UnitRange{Int};
    init::Bool = false,
    use_hillslope::Bool = false,
    nlevsno::Int = varpar.nlevsno
)
    wf = waterflux isa WaterFluxBulkData ? waterflux.wf : waterflux
    isempty(bounds_g) && return nothing
    FT = eltype(l2g.qice_grc)
    nec1 = size(l2g.qice_grc, 1)
    # Fortran `t_soisno(c,1)` is the TOP GROUND layer: its second dimension runs
    # -nlevsno+1 : nlevgrnd, so Fortran index j maps to Julia column j + nlevsno.
    jtop = nlevsno + 1

    # --- defaults (lnd2glcMod.F90:174-193) ---
    for g in bounds_g, n in 1:nec1
        l2g.qice_grc[n, g] = zero(FT)
        l2g.tsrf_grc[n, g] = FT(TFRZ)
        l2g.topo_grc[n, g] = zero(FT)
    end

    isempty(bounds_c) && return nothing

    # Guard against assigning the same (gridcell, class) slot twice.
    fields_assigned = falses(nec1, last(bounds_g))

    for c in bounds_c
        (c <= length(filter_do_smb_c) && filter_do_smb_c[c]) || continue
        l = col.landunit[c]
        g = col.gridcell[c]

        # Vertical index + flux normalization, by landunit type. Other landunit
        # types pass no information in the lnd2glc fields, so they are skipped.
        # (For that to be acceptable CTSM requires a virtual vegetated column in
        # any gridcell made solely of glacier plus another special landunit.)
        if lun.itype[l] == ISTICE
            n = col_itype_to_ice_class(Int(col.itype[c]))
            flux_normalization = 1.0
        elseif lun.itype[l] == ISTSOIL
            n = 0                                        # bare-land slot
            flux_normalization = bareland_normalization(c, col, lun, grc)
        else
            continue
        end

        idx = lnd2glc_class_index(n)
        (1 <= idx <= nec1) || error("update_lnd2glc!: elevation class $n out of range " *
                                    "for maxpatch_glc = $(nec1 - 1) (column $c)")

        if fields_assigned[idx, g] && !use_hillslope
            error("update_lnd2glc!: attempt to assign coupling fields twice for the " *
                  "same index (g = $g, n = $n). One possible cause is having multiple " *
                  "columns in the istsoil landunit, which this routine cannot handle.")
        end
        fields_assigned[idx, g] = true

        # t_soisno and topo_col are valid even in initialization, so tsrf and topo
        # are set regardless of `init`. qflx_glcice is not valid until the run
        # loop, so during initialization qice keeps the default set above.
        l2g.tsrf_grc[idx, g] = temperature.t_soisno_col[c, jtop]
        l2g.topo_grc[idx, g] = topo.topo_col[c]
        if !init
            qice = wf.qflx_glcice_col[c] * FT(flux_normalization)
            l2g.qice_grc[idx, g] = qice
            if abs(qice) > one(FT)
                @warn "update_lnd2glc!: qice out of bounds" gridcell=g class=n qice=qice
            end
        end
    end

    return nothing
end
