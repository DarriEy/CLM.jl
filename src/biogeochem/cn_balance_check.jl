# ==========================================================================
# Ported from: src/biogeochem/CNBalanceCheckMod.F90
# Carbon/nitrogen mass balance checking
# ==========================================================================

"""
    CNBalanceData

Data structure for carbon/nitrogen mass balance checking.
Stores beginning and end-of-timestep masses at column and gridcell levels,
plus warning/error thresholds.

Ported from `cn_balance_type` in `CNBalanceCheckMod.F90`.
"""
# Reparametrized {FT,V} + @adapt_structure so the begin/end mass arrays are
# device-movable (kernels write them on Metal). The 4 threshold scalars stay
# concrete ::Float64 — they are only read in the host-side warn/error scan,
# never inside a device kernel, so the arrays-only _F32 adaptor leaves them
# untouched (V<:AbstractVector{FT} is satisfied independent of the scalars).
Base.@kwdef mutable struct CNBalanceData{FT<:Real,
                              V<:AbstractVector{FT}}
    begcb_col ::V = Float64[]  # (gC/m2) column carbon mass, beginning of time step
    endcb_col ::V = Float64[]  # (gC/m2) column carbon mass, end of time step
    begnb_col ::V = Float64[]  # (gN/m2) column nitrogen mass, beginning of time step
    endnb_col ::V = Float64[]  # (gN/m2) column nitrogen mass, end of time step
    begcb_grc ::V = Float64[]  # (gC/m2) gridcell carbon mass, beginning of time step
    endcb_grc ::V = Float64[]  # (gC/m2) gridcell carbon mass, end of time step
    begnb_grc ::V = Float64[]  # (gN/m2) gridcell nitrogen mass, beginning of time step
    endnb_grc ::V = Float64[]  # (gN/m2) gridcell nitrogen mass, end of time step
    # Last-computed imbalances. Fortran keeps these as subroutine locals; retaining
    # them makes the check INSPECTABLE — the harnesses quantify the non-closure from
    # them, and the tests assert their VALUE rather than merely that a run did not
    # throw. (An `@test isfinite(x)` assertion is nearly worthless: PR #220 shipped a
    # gradient that was silently 7% wrong 98% of the time behind exactly that.)
    errcb_col ::V = Float64[]  # (gC/m2) column carbon imbalance, last step
    errnb_col ::V = Float64[]  # (gN/m2) column nitrogen imbalance, last step
    errcb_grc ::V = Float64[]  # (gC/m2) gridcell carbon imbalance, last step
    errnb_grc ::V = Float64[]  # (gN/m2) gridcell nitrogen imbalance, last step
    cwarning  ::Float64 = 1.0e-8   # (gC/m2) carbon balance warning threshold
    nwarning  ::Float64 = 1.0e-7   # (gN/m2) nitrogen balance warning threshold
    cerror    ::Float64 = 1.0e-7   # (gC/m2) carbon balance error threshold
    nerror    ::Float64 = 1.0e-3   # (gN/m2) nitrogen balance error threshold
end

CNBalanceData{FT}(; kwargs...) where {FT<:Real} =
    CNBalanceData{FT, Vector{FT}}(; kwargs...)
Adapt.@adapt_structure CNBalanceData

"""
    cn_balance_init!(bal::CNBalanceData, nc::Int, ng::Int)

Allocate and initialize a `CNBalanceData` instance for `nc` columns and `ng` gridcells.

Ported from `Init` and `InitAllocate` in `CNBalanceCheckMod.F90`.
"""
# --------------------------------------------------------------------------
# Global enable switch for the CN mass-conservation check.  DEFAULT: ON.
#
# The check is a genuine `error()`-on-failure guard — it is how CTSM catches a
# broken C or N budget, and a conservation check that cannot fail is worse than
# no check at all.
#
# It used to be default-OFF, and the stated reason was a supposed product-pool
# gap on the carbon side. That reason was wrong. The real reason the carbon half
# "could not run" was that EVERY term it reads was structurally dead:
# `gpp_col`/`er_col`/`hr_col`/`nbp_grc` were never aggregated (the column half of
# CNVegCarbonFluxType::Summary was a stub and the two SoilBiogeochem flux
# summaries were never called), and `totc_col` contained no vegetation carbon at
# all. With nothing to compare, the check was structurally incapable of failing.
#
# Once those were wired the check ran — and it FAILED, by 4.5e-2 gC/m2/step
# against a 1e-7 threshold. Chasing that down produced four independent, real
# carbon-destroying/creating bugs (dead `lf_f`/`fr_f` litter fractions; the
# phenology litter routing collapsed to one soil layer; the veg-C total missing
# cpool; and a PFT off-by-one that made every TREE non-woody in the C/N state
# update). See the PR for the full chain.
#
# Bow `use_cn`, summer + autumn windows and a 480-step free run now close to
#     |errcb_col| <= 1.0e-11   and   |errcb_grc| <= 1.0e-11
# against `cerror` = 1e-7 — four orders of margin. So the carbon check is ON and
# FATAL, exactly like the water balance check (#211).
#
# NITROGEN: `n_balance_check!` also runs and is fatal at Fortran's `nerror`
# (1e-3), which it passes (|errnb| ~ 3.1e-4). It does NOT yet clear the 1e-7
# `nwarning`: nitrogen is still being CREATED in the vegetation pools at
# ~3.3e-4 gN/m2/step (the veg N pools gain ~9.9e-4 while the soil supplies only
# ~6.7e-4). That is a REAL, open, documented non-closure — see
# docs/CN_BALANCE_STATUS.md. It is left as a WARNING (not silenced, not
# retuned) rather than an abort, because Fortran's own error threshold is 1e-3
# and we are inside it.
# --------------------------------------------------------------------------
const _CN_BALANCE_CHECK_ENABLED = Ref(true)

"""Return whether the CN mass-conservation check runs inside `clm_drv!`."""
cn_balance_check_enabled() = _CN_BALANCE_CHECK_ENABLED[]

"""
    cn_balance_check_enabled!(on::Bool)

Enable/disable the CN mass-conservation check inside `clm_drv!`. Default **`true`**
(fatal on failure) — see the note above. Turn it off only to *investigate* a
failure, never to hide one.
"""
cn_balance_check_enabled!(on::Bool) = (_CN_BALANCE_CHECK_ENABLED[] = on)

function cn_balance_init!(bal::CNBalanceData, nc::Int, ng::Int)
    bal.begcb_col = fill(NaN, nc)
    bal.endcb_col = fill(NaN, nc)
    bal.begnb_col = fill(NaN, nc)
    bal.endnb_col = fill(NaN, nc)
    bal.begcb_grc = fill(NaN, ng)
    bal.endcb_grc = fill(NaN, ng)
    bal.begnb_grc = fill(NaN, ng)
    bal.endnb_grc = fill(NaN, ng)
    bal.errcb_col = zeros(nc)
    bal.errnb_col = zeros(nc)
    bal.errcb_grc = zeros(ng)
    bal.errnb_grc = zeros(ng)

    bal.cwarning = 1.0e-8
    bal.nwarning = 1.0e-7
    bal.nerror   = 1.0e-3
    bal.cerror   = 1.0e-7
    return nothing
end

# --------------------------------------------------------------------------
# Column-to-gridcell aggregation (unity weighting)
# --------------------------------------------------------------------------

"""
    c2g_unity!(garr, carr, col_gridcell, col_wtgcell, bounds_c, bounds_g)

Column-to-gridcell area-weighted average with unity scaling.
Simple reimplementation of `c2g` from `subgridAveMod.F90` used by balance checks.
"""
@kernel function _c2g_zero_kernel!(garr, gmin::Int, gmax::Int)
    g = @index(Global)
    @inbounds if gmin <= g <= gmax
        garr[g] = zero(eltype(garr))
    end
end

# Per-column scatter into the owning gridcell. KA-CPU iterates threads in
# ascending column order, so `_scatter_add!` reproduces the host `for c` += order
# byte-for-byte; on the GPU it is an atomic add.
@kernel function _c2g_scatter_kernel!(garr, @Const(carr), @Const(col_gridcell),
        @Const(col_wtgcell), cmin::Int, cmax::Int, gmin::Int, gmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax
        g = col_gridcell[c]
        if gmin <= g <= gmax
            _scatter_add!(garr, g, carr[c] * col_wtgcell[c])
        end
    end
end

function c2g_unity!(
    garr::AbstractVector{<:Real},
    carr::AbstractVector{<:Real},
    col_gridcell::AbstractVector{<:Integer},
    col_wtgcell::AbstractVector{<:Real},
    bounds_c::UnitRange{Int},
    bounds_g::UnitRange{Int}
)
    isempty(bounds_g) && return nothing
    _launch!(_c2g_zero_kernel!, garr, first(bounds_g), last(bounds_g);
        ndrange = length(garr))
    isempty(bounds_c) && return nothing
    _launch!(_c2g_scatter_kernel!, garr, carr, col_gridcell, col_wtgcell,
        first(bounds_c), last(bounds_c), first(bounds_g), last(bounds_g);
        ndrange = length(carr))
    return nothing
end

# --------------------------------------------------------------------------
# c2g with Fortran's 'urbanf' column-to-landunit scaling
# --------------------------------------------------------------------------
#
# Fortran's WATER balance check aggregates every column term to the gridcell with
#     call c2g(..., c2l_scale_type='urbanf', l2g_scale_type='unity')
# (BalanceCheckMod.F90:270, :332, :713-724). 'urbanf' is the scaling for
# EXTENSIVE, per-m^2 urban quantities: an urban WALL flux is expressed per m^2 of
# VERTICAL WALL area, and a canyon-floor flux per m^2 of ROAD area, so both must
# be converted to per-m^2 of GROUND area before they can be area-averaged with the
# rest of the gridcell (subgridAveMod::set_c2l_scale):
#
#     sunwall / shadewall -> 3 * canyon_hwr
#     road (perv/imperv)  -> 3
#     roof                -> 1
#     everything non-urban-> 1
#
# The port was using `c2g_unity!` (scale == 1 everywhere) for the water balance —
# and its docstring even CLAIMED that was "equivalent to c2l_scale_type='urbanf'".
# It is not. On a gridcell with urban landunits the wall/road columns entered the
# gridcell sum unscaled, so `errh2o_grc` could not close even when EVERY COLUMN
# closed to machine precision (mexicocity: columns ~1e-13, gridcell -2.6e-4 mm).
# This was invisible until the water InitCold + urban-ponding fixes made those
# columns finite and closing.
#
# On a gridcell with NO urban landunit every scale is 1.0, so this is EXACTLY
# `c2g_unity!` there — the non-urban domains stay bit-identical.
#
# The scale is recomputed inside the kernel from fields that are already
# device-resident (col.itype / col.landunit / lun.urbpoi / lun.canyon_hwr), so
# there is no host scale vector to allocate or move to the GPU.
@inline function _c2l_urbanf_scale(::Type{T}, c, col_itype, col_landunit,
                                   lun_urbpoi, lun_canyon_hwr) where {T}
    l = col_landunit[c]
    (l >= 1 && l <= length(lun_urbpoi) && lun_urbpoi[l]) || return one(T)
    it = col_itype[c]
    if it == ICOL_SUNWALL || it == ICOL_SHADEWALL
        # canyon_hwr is populated by urbanparams_populate! on every real urban run
        # (clm_initialize! Step 13a). Minimal unit-test fixtures that flag a landunit
        # urban WITHOUT loading the morphology leave it unallocated; fall back to the
        # unity scale there rather than reading out of bounds.
        return l <= length(lun_canyon_hwr) ? T(3) * T(lun_canyon_hwr[l]) : one(T)
    elseif it == ICOL_ROAD_PERV || it == ICOL_ROAD_IMPERV
        return T(3)
    else                                    # ICOL_ROOF (and any other urban col)
        return one(T)
    end
end

@kernel function _c2g_scatter_urbanf_kernel!(garr, @Const(carr), @Const(col_gridcell),
        @Const(col_wtgcell), @Const(col_itype), @Const(col_landunit),
        @Const(lun_urbpoi), @Const(lun_canyon_hwr),
        cmin::Int, cmax::Int, gmin::Int, gmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax
        g = col_gridcell[c]
        if gmin <= g <= gmax
            T = eltype(garr)
            s = _c2l_urbanf_scale(T, c, col_itype, col_landunit,
                                  lun_urbpoi, lun_canyon_hwr)
            _scatter_add!(garr, g, carr[c] * s * col_wtgcell[c])
        end
    end
end

"""
    c2g_urbanf!(garr, carr, col, lun, bounds_c, bounds_g)

Column-to-gridcell area-weighted sum with Fortran's `c2l_scale_type='urbanf'`
scaling. Use this — not [`c2g_unity!`](@ref) — for every term of the WATER
balance, matching `BalanceCheckMod.F90`. Identical to `c2g_unity!` on gridcells
with no urban landunit.
"""
function c2g_urbanf!(
    garr::AbstractVector{<:Real},
    carr::AbstractVector{<:Real},
    col,
    lun,
    bounds_c::UnitRange{Int},
    bounds_g::UnitRange{Int}
)
    isempty(bounds_g) && return nothing
    _launch!(_c2g_zero_kernel!, garr, first(bounds_g), last(bounds_g);
        ndrange = length(garr))
    isempty(bounds_c) && return nothing
    _launch!(_c2g_scatter_urbanf_kernel!, garr, carr, col.gridcell, col.wtgcell,
        col.itype, col.landunit, lun.urbpoi, lun.canyon_hwr,
        first(bounds_c), last(bounds_c), first(bounds_g), last(bounds_g);
        ndrange = length(carr))
    return nothing
end

# --------------------------------------------------------------------------
# BeginCNGridcellBalance
# --------------------------------------------------------------------------

"""
    begin_cn_gridcell_balance!(bal, soilbgc_cstate, soilbgc_nstate,
        c_products, n_products, bounds_g;
        use_fates_bgc, hrv_xsmrpool_amount_left, gru_conv_cflux_amount_left,
        dwt_conv_cflux_amount_left)

Calculate beginning gridcell-level carbon/nitrogen balance for mass conservation check.

Should be called after CN state summaries have been computed and before
the dynamic landunit area updates.

Ported from `BeginCNGridcellBalance` in `CNBalanceCheckMod.F90`.

The dribbler amounts (`hrv_xsmrpool_amount_left`, `gru_conv_cflux_amount_left`,
`dwt_conv_cflux_amount_left`) must be pre-computed by the caller when
`use_fates_bgc` is false. When `use_fates_bgc` is true, they are ignored
and treated as zero.
"""
# One thread per gridcell. The dribbler arrays are only indexed on the
# !use_fates_bgc branch, and only when the dribblers actually exist: a config that
# never allocated them leaves them EMPTY, and reading them is an out-of-bounds read
# that `@inbounds` hides (it throws only under --check-bounds=yes). `has_dribbler`
# encodes `!isempty(hrv_xsmrpool_amount_left)` — the same guard the matching END
# kernel `_cbal_grc_kernel!` already applies. Without it, the begin/end balance pair
# was ASYMMETRIC as well as unsafe: the end kernel skipped the dribblers and the
# begin kernel read past the end of an empty vector.
@kernel function _cnbal_begin_grc_kernel!(begcb_grc, begnb_grc,
        @Const(totc), @Const(totn), @Const(c_tot_woodprod), @Const(c_cropprod1),
        @Const(n_tot_woodprod), @Const(n_cropprod1),
        @Const(hrv_left), @Const(gru_left), @Const(dwt_left), has_dribbler::Bool,
        use_fates_bgc::Bool, gmin::Int, gmax::Int)
    g = @index(Global)
    @inbounds if gmin <= g <= gmax
        begcb_grc[g] = totc[g] + c_tot_woodprod[g] + c_cropprod1[g]
        if !use_fates_bgc && has_dribbler
            begcb_grc[g] += hrv_left[g] + gru_left[g] + dwt_left[g]
        end
        begnb_grc[g] = totn[g] + n_tot_woodprod[g] + n_cropprod1[g]
    end
end

function begin_cn_gridcell_balance!(
    bal::CNBalanceData,
    soilbgc_cstate::SoilBiogeochemCarbonStateData,
    soilbgc_nstate::SoilBiogeochemNitrogenStateData,
    c_products::CNProductsData,
    n_products::CNProductsData,
    bounds_g::UnitRange{Int};
    use_fates_bgc::Bool=false,
    hrv_xsmrpool_amount_left::AbstractVector{<:Real}=Float64[],
    gru_conv_cflux_amount_left::AbstractVector{<:Real}=Float64[],
    dwt_conv_cflux_amount_left::AbstractVector{<:Real}=Float64[]
)
    totc = soilbgc_cstate.totc_grc
    totn = soilbgc_nstate.totn_grc
    c_cropprod1    = c_products.cropprod1_grc
    n_cropprod1    = n_products.cropprod1_grc
    c_tot_woodprod = c_products.tot_woodprod_grc
    n_tot_woodprod = n_products.tot_woodprod_grc

    isempty(bounds_g) && return nothing
    # Same dribbler-presence guard as the END kernel (`_cbal_grc_kernel!`).
    has_dribbler = !isempty(hrv_xsmrpool_amount_left)
    _launch!(_cnbal_begin_grc_kernel!, bal.begcb_grc, bal.begnb_grc,
        totc, totn, c_tot_woodprod, c_cropprod1, n_tot_woodprod, n_cropprod1,
        hrv_xsmrpool_amount_left, gru_conv_cflux_amount_left, dwt_conv_cflux_amount_left,
        has_dribbler, use_fates_bgc, first(bounds_g), last(bounds_g);
        ndrange = length(bal.begcb_grc))
    return nothing
end

# --------------------------------------------------------------------------
# BeginCNColumnBalance
# --------------------------------------------------------------------------

"""
    begin_cn_column_balance!(bal, soilbgc_cstate, soilbgc_nstate,
        mask_soil, bounds_c)

Calculate beginning column-level carbon/nitrogen balance for mass conservation check.

Should be called after CN state summaries have been recomputed for this time step
(after dynamic landunit area updates and associated filter updates).

Ported from `BeginCNColumnBalance` in `CNBalanceCheckMod.F90`.
"""
@kernel function _cnbal_begin_col_kernel!(begcb_col, begnb_col,
        @Const(mask_soil), @Const(totcolc), @Const(totcoln), cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && mask_soil[c]
        begcb_col[c] = totcolc[c]
        begnb_col[c] = totcoln[c]
    end
end

function begin_cn_column_balance!(
    bal::CNBalanceData,
    soilbgc_cstate::SoilBiogeochemCarbonStateData,
    soilbgc_nstate::SoilBiogeochemNitrogenStateData,
    mask_soil::AbstractVector{Bool},
    bounds_c::UnitRange{Int}
)
    totcolc = soilbgc_cstate.totc_col
    totcoln = soilbgc_nstate.totn_col

    isempty(bounds_c) && return nothing
    _launch!(_cnbal_begin_col_kernel!, bal.begcb_col, bal.begnb_col,
        mask_soil, totcolc, totcoln, first(bounds_c), last(bounds_c);
        ndrange = length(bal.begcb_col))
    return nothing
end

# --------------------------------------------------------------------------
# CBalanceCheck
# --------------------------------------------------------------------------

"""
    c_balance_check!(bal, soilbgc_cflux, soilbgc_cstate, cnveg_cflux,
        c_products, col_data, grc_data,
        mask_soil, bounds_c, bounds_g, dt;
        is_fates_col, use_fates_bgc,
        hrv_xsmrpool_amount_left, gru_conv_cflux_amount_left,
        dwt_conv_cflux_amount_left)

Perform carbon mass conservation check for column and gridcell.

Ported from `CBalanceCheck` in `CNBalanceCheckMod.F90`.

# Arguments
- `bal`: CNBalanceData with beginning-of-step masses already set
- `soilbgc_cflux`: SoilBiogeochemCarbonFluxData
- `soilbgc_cstate`: SoilBiogeochemCarbonStateData
- `cnveg_cflux`: CNVegCarbonFluxData
- `c_products`: CNProductsData for carbon
- `col_data`: ColumnData for subgrid mapping
- `grc_data`: GridcellData for lat/lon
- `mask_soil`: BitVector mask for soil columns
- `bounds_c`: column index range
- `bounds_g`: gridcell index range
- `dt`: timestep size (seconds)
- `is_fates_col`: per-column vector indicating FATES columns
- `use_fates_bgc`: global flag for FATES BGC
- `hrv_xsmrpool_amount_left`: dribbler amount (end-of-step, gridcell)
- `gru_conv_cflux_amount_left`: dribbler amount (end-of-step, gridcell)
- `dwt_conv_cflux_amount_left`: dribbler amount (end-of-step, gridcell)
"""
# One thread per column: writes col_endcb and col_errcb. has_fates encodes
# `!isempty(is_fates_col)` (when false the empty is_fates_col is never indexed).
# T(dt) keeps the timestep arithmetic in the device eltype (Metal has no Float64);
# on CPU T==Float64 so it is byte-identical to the original loop.
@kernel function _cbal_col_kernel!(col_endcb, col_errcb,
        @Const(mask_soil), @Const(is_fates_col), has_fates::Bool,
        @Const(totcolc), @Const(col_begcb), @Const(gpp), @Const(er),
        @Const(col_fire_closs), @Const(col_hrv_xsmrpool_to_atm),
        @Const(col_xsmrpool_to_atm), @Const(gru_conv_cflux), @Const(wood_harvestc),
        @Const(gru_wood_productc_gain), @Const(crop_harvestc_to_cropprodc),
        @Const(som_c_leached), @Const(fates_litter_flux), @Const(hr_col),
        dt, cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && mask_soil[c]
        T = eltype(col_errcb)
        col_endcb[c] = totcolc[c]
        if has_fates && is_fates_col[c]
            col_cinputs = fates_litter_flux[c]
            col_coutputs = hr_col[c]
        else
            col_cinputs = gpp[c]
            col_coutputs = er[c] + col_fire_closs[c] + col_hrv_xsmrpool_to_atm[c] +
                           col_xsmrpool_to_atm[c] + gru_conv_cflux[c]
            col_coutputs += wood_harvestc[c] +
                            gru_wood_productc_gain[c] +
                            crop_harvestc_to_cropprodc[c]
        end
        col_coutputs -= som_c_leached[c]
        col_errcb[c] = (col_cinputs - col_coutputs) * T(dt) -
                        (col_endcb[c] - col_begcb[c])
    end
end

# One thread per gridcell: writes grc_endcb and grc_errcb. has_dribbler encodes
# `!isempty(hrv_xsmrpool_amount_left)`.
@kernel function _cbal_grc_kernel!(grc_endcb, grc_errcb,
        @Const(grc_begcb), @Const(totgrcc), @Const(tot_woodprod_grc), @Const(cropprod1_grc),
        @Const(hrv_left), @Const(gru_left), @Const(dwt_left), has_dribbler::Bool,
        @Const(nbp_grc), @Const(dwt_seedc_to_leaf_grc), @Const(dwt_seedc_to_deadstem_grc),
        @Const(som_c_leached_grc), use_fates_bgc::Bool, dt, gmin::Int, gmax::Int)
    g = @index(Global)
    @inbounds if gmin <= g <= gmax
        T = eltype(grc_errcb)
        if !use_fates_bgc
            grc_endcb[g] = totgrcc[g] + tot_woodprod_grc[g] + cropprod1_grc[g]
            if has_dribbler
                grc_endcb[g] += hrv_left[g] + gru_left[g] + dwt_left[g]
            end
            grc_cinputs = nbp_grc[g] +
                          dwt_seedc_to_leaf_grc[g] + dwt_seedc_to_deadstem_grc[g]
            grc_coutputs = -som_c_leached_grc[g]
            grc_errcb[g] = (grc_cinputs - grc_coutputs) * T(dt) -
                            (grc_endcb[g] - grc_begcb[g])
        else
            grc_endcb[g] = grc_begcb[g]
            grc_errcb[g] = zero(T)
        end
    end
end

function c_balance_check!(
    bal::CNBalanceData,
    soilbgc_cflux::SoilBiogeochemCarbonFluxData,
    soilbgc_cstate::SoilBiogeochemCarbonStateData,
    cnveg_cflux::CNVegCarbonFluxData,
    c_products::CNProductsData,
    col_data::ColumnData,
    grc_data::GridcellData,
    mask_soil::AbstractVector{Bool},
    bounds_c::UnitRange{Int},
    bounds_g::UnitRange{Int},
    dt::Real;
    is_fates_col::AbstractVector{Bool}=Bool[],
    use_fates_bgc::Bool=false,
    hrv_xsmrpool_amount_left::AbstractVector{<:Real}=Float64[],
    gru_conv_cflux_amount_left::AbstractVector{<:Real}=Float64[],
    dwt_conv_cflux_amount_left::AbstractVector{<:Real}=Float64[]
)
    # Local aliases for column-level fields
    col_begcb = bal.begcb_col
    col_endcb = bal.endcb_col
    grc_begcb = bal.begcb_grc
    grc_endcb = bal.endcb_grc
    totcolc   = soilbgc_cstate.totc_col

    gpp                     = cnveg_cflux.gpp_col
    er                      = cnveg_cflux.er_col
    col_fire_closs          = cnveg_cflux.fire_closs_col
    col_hrv_xsmrpool_to_atm = cnveg_cflux.hrv_xsmrpool_to_atm_col
    col_xsmrpool_to_atm     = cnveg_cflux.xsmrpool_to_atm_col
    wood_harvestc           = cnveg_cflux.wood_harvestc_col
    gru_conv_cflux          = cnveg_cflux.gru_conv_cflux_col
    gru_wood_productc_gain  = cnveg_cflux.gru_wood_productc_gain_col
    crop_harvestc_to_cropprodc = cnveg_cflux.crop_harvestc_to_cropprodc_col

    som_c_leached    = soilbgc_cflux.som_c_leached_col
    hr_col           = soilbgc_cflux.hr_col
    fates_litter_flux = soilbgc_cflux.fates_litter_flux

    nbp_grc                   = cnveg_cflux.nbp_grc
    dwt_seedc_to_leaf_grc     = cnveg_cflux.dwt_seedc_to_leaf_grc
    dwt_seedc_to_deadstem_grc = cnveg_cflux.dwt_seedc_to_deadstem_grc

    cropprod1_grc    = c_products.cropprod1_grc
    tot_woodprod_grc = c_products.tot_woodprod_grc

    # Column-level error array (device-resident scratch via similar()).
    FT = eltype(bal.begcb_col)
    nc = length(bounds_c) > 0 ? last(bounds_c) : 0
    col_errcb = fill!(similar(bal.begcb_col, FT, nc), zero(FT))
    has_fates = !isempty(is_fates_col)

    # --- Column-level balance check (numeric pass on device/host) ---
    if !isempty(bounds_c)
        _launch!(_cbal_col_kernel!, col_endcb, col_errcb,
            mask_soil, is_fates_col, has_fates,
            totcolc, col_begcb, gpp, er, col_fire_closs, col_hrv_xsmrpool_to_atm,
            col_xsmrpool_to_atm, gru_conv_cflux, wood_harvestc, gru_wood_productc_gain,
            crop_harvestc_to_cropprodc, som_c_leached, fates_litter_flux, hr_col,
            FT(dt), first(bounds_c), last(bounds_c); ndrange = length(col_endcb))
    end

    # Publish the imbalance so callers/tests can assert its VALUE (see CNBalanceData).
    if length(bal.errcb_col) == length(col_errcb)
        copyto!(bal.errcb_col, col_errcb)
    end

    # --- Host-only warn/error scan (col_errcb is a plain Array only on the CPU;
    #     on the GPU we skip the String-building error path entirely) ---
    err_found = false
    err_index = 0
    if col_errcb isa Array
        for c in bounds_c
            mask_soil[c] || continue
            if abs(col_errcb[c]) > bal.cerror
                err_found = true
                err_index = c
            end
            if abs(col_errcb[c]) > bal.cwarning
                @warn "cbalance warning at c = $c" col_errcb=col_errcb[c] col_endcb=col_endcb[c]
            end
        end
    end

    if err_found
        c = err_index
        is_fates = !isempty(is_fates_col) && is_fates_col[c]
        g = col_data.gridcell[c]
        msg = string(
            "column cbalance error    = ", col_errcb[c], " c=", c, "\n",
            "is fates column?         = ", is_fates, "\n",
            "Latdeg,Londeg            = ", grc_data.latdeg[g], ",", grc_data.londeg[g], "\n",
            "begcb                    = ", col_begcb[c], "\n",
            "endcb                    = ", col_endcb[c], "\n",
            "delta store              = ", col_endcb[c] - col_begcb[c], "\n",
            "--- Inputs ---\n",
            is_fates ? "fates litter_flux        = $(fates_litter_flux[c]*dt)\n" :
                       "gpp                      = $(gpp[c]*dt)\n",
            "--- Outputs ---\n",
            !is_fates ? string(
                "er                       = ", er[c]*dt, "\n",
                "col_fire_closs           = ", col_fire_closs[c]*dt, "\n",
                "col_hrv_xsmrpool_to_atm  = ", col_hrv_xsmrpool_to_atm[c]*dt, "\n",
                "col_xsmrpool_to_atm      = ", col_xsmrpool_to_atm[c]*dt, "\n",
                "wood_harvestc            = ", wood_harvestc[c]*dt, "\n",
                "crop_harvestc_to_cropprodc = ", crop_harvestc_to_cropprodc[c]*dt, "\n"
            ) : string("hr                       = ", hr_col[c]*dt, "\n"),
            "-1*som_c_leached         = ", som_c_leached[c]*dt
        )
        error("CNBalanceCheck carbon error:\n$msg")
    end

    # --- Gridcell-level balance check ---

    # Column-to-gridcell aggregation (kernelized scatter; device-resident scratch)
    totgrcc = soilbgc_cstate.totc_grc
    c2g_unity!(totgrcc, totcolc, col_data.gridcell, col_data.wtgcell, bounds_c, bounds_g)

    ng = length(grc_data.latdeg)
    som_c_leached_grc = fill!(similar(bal.begcb_grc, FT, ng), zero(FT))
    c2g_unity!(som_c_leached_grc, som_c_leached, col_data.gridcell, col_data.wtgcell, bounds_c, bounds_g)

    grc_errcb = fill!(similar(bal.begcb_grc, FT, ng), zero(FT))
    has_dribbler = !isempty(hrv_xsmrpool_amount_left)

    if !isempty(bounds_g)
        _launch!(_cbal_grc_kernel!, grc_endcb, grc_errcb,
            grc_begcb, totgrcc, tot_woodprod_grc, cropprod1_grc,
            hrv_xsmrpool_amount_left, gru_conv_cflux_amount_left, dwt_conv_cflux_amount_left,
            has_dribbler, nbp_grc, dwt_seedc_to_leaf_grc, dwt_seedc_to_deadstem_grc,
            som_c_leached_grc, use_fates_bgc, FT(dt), first(bounds_g), last(bounds_g);
            ndrange = length(grc_endcb))
    end

    if length(bal.errcb_grc) == length(grc_errcb)
        copyto!(bal.errcb_grc, grc_errcb)
    end

    err_found = false
    err_index = 0
    if grc_errcb isa Array
        for g in bounds_g
            if abs(grc_errcb[g]) > bal.cerror
                err_found = true
                err_index = g
            end
            if abs(grc_errcb[g]) > bal.cwarning
                @warn "cbalance warning at g = $g" grc_errcb=grc_errcb[g] grc_endcb=grc_endcb[g]
            end
        end
    end

    if err_found
        g = err_index
        msg = string(
            "gridcell cbalance error  = ", grc_errcb[g], " g=", g, "\n",
            "latdeg, londeg           = ", grc_data.latdeg[g], ",", grc_data.londeg[g], "\n",
            "begcb                    = ", grc_begcb[g], "\n",
            "endcb                    = ", grc_endcb[g], "\n",
            "delta store              = ", grc_endcb[g] - grc_begcb[g], "\n",
            "--- Inputs ---\n",
            "nbp_grc                  = ", nbp_grc[g] * dt, "\n",
            "dwt_seedc_to_leaf_grc    = ", dwt_seedc_to_leaf_grc[g] * dt, "\n",
            "dwt_seedc_to_deadstem_grc = ", dwt_seedc_to_deadstem_grc[g] * dt, "\n",
            "--- Outputs ---\n",
            "-1*som_c_leached_grc     = ", som_c_leached_grc[g] * dt
        )
        error("CNBalanceCheck gridcell carbon error:\n$msg")
    end

    return nothing
end

# --------------------------------------------------------------------------
# NBalanceCheck
# --------------------------------------------------------------------------

"""
    n_balance_check!(bal, soilbgc_nflux, soilbgc_nstate, cnveg_nflux,
        n_products, col_data, grc_data,
        mask_soil, bounds_c, bounds_g, dt;
        is_fates_col, use_fates_bgc, use_nitrif_denitrif,
        use_crop, use_fun)

Perform nitrogen mass conservation check for column and gridcell.

Ported from `NBalanceCheck` in `CNBalanceCheckMod.F90`.

# Arguments
- `bal`: CNBalanceData with beginning-of-step masses already set
- `soilbgc_nflux`: SoilBiogeochemNitrogenFluxData
- `soilbgc_nstate`: SoilBiogeochemNitrogenStateData
- `cnveg_nflux`: CNVegNitrogenFluxData
- `n_products`: CNProductsData for nitrogen
- `col_data`: ColumnData for subgrid mapping
- `grc_data`: GridcellData for lat/lon
- `mask_soil`: BitVector mask for soil columns
- `bounds_c`: column index range
- `bounds_g`: gridcell index range
- `dt`: timestep size (seconds)
- `is_fates_col`: per-column flag for FATES columns
- `use_fates_bgc`: global FATES BGC flag
- `use_nitrif_denitrif`: flag for nitrification/denitrification
- `use_crop`: flag for crop model
- `use_fun`: flag for FUN model
"""
# Device-view bundle: the N column kernel touches 21 per-column input arrays,
# which would blow Metal's ~31 TOTAL-arg limit if passed loose. They are grouped
# into one immutable Adapt-able struct (one kernel arg) and aliased back to their
# Fortran-named locals at the kernel top so the body stays verbatim.
Base.@kwdef struct _NBalColIn{V}
    totcoln::V; ndep_to_sminn::V; nfix_to_sminn::V; ffix_to_sminn::V
    fert_to_sminn::V; soyfixn_to_sminn::V; supplement_to_sminn::V; denit::V
    sminn_leached::V; smin_no3_leached::V; smin_no3_runoff::V; f_n2o_nit::V
    som_n_leached::V; sminn_to_plant::V; fates_litter_flux::V; col_fire_nloss::V
    wood_harvestn::V; gru_conv_nflux::V; gru_wood_productn_gain::V
    crop_harvestn_to_cropprodn::V; col_begnb::V
end
Adapt.@adapt_structure _NBalColIn

# One thread per column. Per-column accumulation uses locals (`ninp`, `nout`,
# `noutp`) written once into the output arrays — byte-identical to the original
# += chain, and avoids the --check-bounds mis-lowering of repeated `out[c]+=`.
@kernel function _nbal_col_kernel!(col_endnb, col_ninputs, col_noutputs,
        col_ninputs_partial, col_noutputs_partial, col_errnb,
        nb::_NBalColIn, @Const(mask_soil), @Const(is_fates_col),
        has_fates::Bool, use_fun::Bool, use_crop::Bool, use_nitrif_denitrif::Bool,
        dt, cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && mask_soil[c]
        T = eltype(col_errnb)
        totcoln = nb.totcoln; ndep_to_sminn = nb.ndep_to_sminn; nfix_to_sminn = nb.nfix_to_sminn
        ffix_to_sminn = nb.ffix_to_sminn; fert_to_sminn = nb.fert_to_sminn
        soyfixn_to_sminn = nb.soyfixn_to_sminn; supplement_to_sminn = nb.supplement_to_sminn
        denit = nb.denit; sminn_leached = nb.sminn_leached; smin_no3_leached = nb.smin_no3_leached
        smin_no3_runoff = nb.smin_no3_runoff; f_n2o_nit = nb.f_n2o_nit; som_n_leached = nb.som_n_leached
        sminn_to_plant = nb.sminn_to_plant; fates_litter_flux = nb.fates_litter_flux
        col_fire_nloss = nb.col_fire_nloss; wood_harvestn = nb.wood_harvestn
        gru_conv_nflux = nb.gru_conv_nflux; gru_wood_productn_gain = nb.gru_wood_productn_gain
        crop_harvestn_to_cropprodn = nb.crop_harvestn_to_cropprodn; col_begnb = nb.col_begnb

        is_fates = has_fates && is_fates_col[c]
        col_endnb[c] = totcoln[c]

        ninp = ndep_to_sminn[c] + nfix_to_sminn[c] + supplement_to_sminn[c]
        if is_fates
            ninp += fates_litter_flux[c]
        end
        if use_fun
            ninp += ffix_to_sminn[c]
        end
        if use_crop
            ninp += fert_to_sminn[c] + soyfixn_to_sminn[c]
        end
        col_ninputs[c] = ninp
        col_ninputs_partial[c] = ninp

        nout = denit[c]
        if !is_fates
            nout += col_fire_nloss[c] + gru_conv_nflux[c]
            nout += wood_harvestn[c] + gru_wood_productn_gain[c] + crop_harvestn_to_cropprodn[c]
        else
            nout += sminn_to_plant[c]
        end
        if !use_nitrif_denitrif
            nout += sminn_leached[c]
        else
            nout += f_n2o_nit[c]
            nout += smin_no3_leached[c] + smin_no3_runoff[c]
        end
        nout -= som_n_leached[c]
        col_noutputs[c] = nout

        noutp = nout
        if !is_fates
            noutp -= wood_harvestn[c] + crop_harvestn_to_cropprodn[c]
        end
        col_noutputs_partial[c] = noutp

        col_errnb[c] = (ninp - nout) * T(dt) - (col_endnb[c] - col_begnb[c])
    end
end

# One thread per gridcell (only launched when !use_fates_bgc).
@kernel function _nbal_grc_kernel!(grc_endnb, grc_errnb,
        @Const(grc_begnb), @Const(totgrcn), @Const(tot_woodprod_grc), @Const(cropprod1_grc),
        @Const(grc_ninputs_partial), @Const(grc_noutputs_partial),
        @Const(dwt_seedn_to_leaf_grc), @Const(dwt_seedn_to_deadstem_grc),
        @Const(dwt_conv_nflux_grc), @Const(product_loss_grc), @Const(gru_wood_productn_gain_grc),
        dt, gmin::Int, gmax::Int)
    g = @index(Global)
    @inbounds if gmin <= g <= gmax
        T = eltype(grc_errnb)
        grc_endnb[g] = totgrcn[g] + tot_woodprod_grc[g] + cropprod1_grc[g]
        grc_ninputs = grc_ninputs_partial[g] +
                      dwt_seedn_to_leaf_grc[g] + dwt_seedn_to_deadstem_grc[g]
        grc_noutputs = grc_noutputs_partial[g] +
                       dwt_conv_nflux_grc[g] + product_loss_grc[g] -
                       gru_wood_productn_gain_grc[g]
        grc_errnb[g] = (grc_ninputs - grc_noutputs) * T(dt) -
                        (grc_endnb[g] - grc_begnb[g])
    end
end

function n_balance_check!(
    bal::CNBalanceData,
    soilbgc_nflux::SoilBiogeochemNitrogenFluxData,
    soilbgc_nstate::SoilBiogeochemNitrogenStateData,
    cnveg_nflux::CNVegNitrogenFluxData,
    n_products::CNProductsData,
    col_data::ColumnData,
    grc_data::GridcellData,
    mask_soil::AbstractVector{Bool},
    bounds_c::UnitRange{Int},
    bounds_g::UnitRange{Int},
    dt::Real;
    is_fates_col::AbstractVector{Bool}=Bool[],
    use_fates_bgc::Bool=false,
    use_nitrif_denitrif::Bool=true,
    use_crop::Bool=false,
    use_fun::Bool=false
)
    # Local aliases
    col_begnb = bal.begnb_col
    col_endnb = bal.endnb_col
    grc_begnb = bal.begnb_grc
    grc_endnb = bal.endnb_grc

    totcoln           = soilbgc_nstate.totn_col
    ndep_to_sminn     = soilbgc_nflux.ndep_to_sminn_col
    nfix_to_sminn     = soilbgc_nflux.nfix_to_sminn_col
    ffix_to_sminn     = soilbgc_nflux.ffix_to_sminn_col
    fert_to_sminn     = soilbgc_nflux.fert_to_sminn_col
    soyfixn_to_sminn  = soilbgc_nflux.soyfixn_to_sminn_col
    supplement_to_sminn = soilbgc_nflux.supplement_to_sminn_col
    denit             = soilbgc_nflux.denit_col
    sminn_leached     = soilbgc_nflux.sminn_leached_col
    smin_no3_leached  = soilbgc_nflux.smin_no3_leached_col
    smin_no3_runoff   = soilbgc_nflux.smin_no3_runoff_col
    f_n2o_nit         = soilbgc_nflux.f_n2o_nit_col
    som_n_leached     = soilbgc_nflux.som_n_leached_col
    sminn_to_plant    = soilbgc_nflux.sminn_to_plant_col
    fates_litter_flux = soilbgc_nflux.fates_litter_flux

    col_fire_nloss    = cnveg_nflux.fire_nloss_col
    wood_harvestn     = cnveg_nflux.wood_harvestn_col
    gru_conv_nflux    = cnveg_nflux.gru_conv_nflux_col
    gru_wood_productn_gain = cnveg_nflux.gru_wood_productn_gain_col
    crop_harvestn_to_cropprodn = cnveg_nflux.crop_harvestn_to_cropprodn_col

    dwt_seedn_to_leaf_grc     = cnveg_nflux.dwt_seedn_to_leaf_grc
    dwt_seedn_to_deadstem_grc = cnveg_nflux.dwt_seedn_to_deadstem_grc
    dwt_conv_nflux_grc        = cnveg_nflux.dwt_conv_nflux_grc
    gru_conv_nflux_grc        = cnveg_nflux.gru_conv_nflux_grc
    gru_wood_productn_gain_grc = cnveg_nflux.gru_wood_productn_gain_grc

    cropprod1_grc    = n_products.cropprod1_grc
    tot_woodprod_grc = n_products.tot_woodprod_grc
    product_loss_grc = n_products.product_loss_grc

    # Column-level scratch (device-resident via similar()).
    FT = eltype(bal.begnb_col)
    nc = length(bounds_c) > 0 ? last(bounds_c) : 0
    col_ninputs          = fill!(similar(bal.begnb_col, FT, nc), zero(FT))
    col_noutputs         = fill!(similar(bal.begnb_col, FT, nc), zero(FT))
    col_errnb            = fill!(similar(bal.begnb_col, FT, nc), zero(FT))
    col_ninputs_partial  = fill!(similar(bal.begnb_col, FT, nc), zero(FT))
    col_noutputs_partial = fill!(similar(bal.begnb_col, FT, nc), zero(FT))
    has_fates = !isempty(is_fates_col)

    nb_in = _NBalColIn(; totcoln, ndep_to_sminn, nfix_to_sminn, ffix_to_sminn,
        fert_to_sminn, soyfixn_to_sminn, supplement_to_sminn, denit, sminn_leached,
        smin_no3_leached, smin_no3_runoff, f_n2o_nit, som_n_leached, sminn_to_plant,
        fates_litter_flux, col_fire_nloss, wood_harvestn, gru_conv_nflux,
        gru_wood_productn_gain, crop_harvestn_to_cropprodn, col_begnb)

    if !isempty(bounds_c)
        backend = _kernel_backend(col_endnb)
        _nbal_col_kernel!(backend)(col_endnb, col_ninputs, col_noutputs,
            col_ninputs_partial, col_noutputs_partial, col_errnb,
            nb_in, mask_soil, is_fates_col,
            has_fates, use_fun, use_crop, use_nitrif_denitrif,
            FT(dt), first(bounds_c), last(bounds_c); ndrange = length(col_endnb))
        KA.synchronize(backend)
    end

    if length(bal.errnb_col) == length(col_errnb)
        copyto!(bal.errnb_col, col_errnb)
    end

    # Host-only warn/error scan (CPU path only; device just computes numbers).
    err_found = false
    err_index = 0
    if col_errnb isa Array
        for c in bounds_c
            mask_soil[c] || continue
            if abs(col_errnb[c]) > bal.nerror
                err_found = true
                err_index = c
            end
            if abs(col_errnb[c]) > bal.nwarning
                # KNOWN OPEN RESIDUAL — see docs/CN_BALANCE_STATUS.md. The N budget
                # clears Fortran's `nerror` (1e-3) but not its `nwarning` (1e-7):
                # nitrogen is still created in the VEGETATION pools at ~3.3e-4
                # gN/m2/step. Neither threshold is retuned and the check is not
                # silenced — `maxlog` only stops one documented, non-fatal residual
                # from emitting an identical line on all 480 steps of a free run.
                @warn "nbalance warning at c = $c (KNOWN open residual — see docs/CN_BALANCE_STATUS.md)" col_errnb=col_errnb[c] col_endnb=col_endnb[c] inputs_ffix=ffix_to_sminn[c]*dt inputs_nfix=nfix_to_sminn[c]*dt inputs_ndep=ndep_to_sminn[c]*dt outputs_lch=smin_no3_leached[c]*dt outputs_roff=smin_no3_runoff[c]*dt outputs_dnit=f_n2o_nit[c]*dt maxlog=3
            end
        end
    end

    if err_found
        c = err_index
        is_fates = !isempty(is_fates_col) && is_fates_col[c]
        g = col_data.gridcell[c]
        msg = string(
            "column nbalance error    = ", col_errnb[c], " c=", c, "\n",
            "Latdeg,Londeg            = ", grc_data.latdeg[g], ",", grc_data.londeg[g], "\n",
            "begnb                    = ", col_begnb[c], "\n",
            "endnb                    = ", col_endnb[c], "\n",
            "delta store              = ", col_endnb[c] - col_begnb[c], "\n",
            "input mass               = ", col_ninputs[c] * dt, "\n",
            "output mass              = ", col_noutputs[c] * dt, "\n",
            "net flux                 = ", (col_ninputs[c] - col_noutputs[c]) * dt, "\n",
            is_fates ?
                "inputs,ndep,nfix,suppn   = $(ndep_to_sminn[c]*dt),$(nfix_to_sminn[c]*dt),$(supplement_to_sminn[c]*dt)\n" :
                "inputs,ffix,nfix,ndep    = $(ffix_to_sminn[c]*dt),$(nfix_to_sminn[c]*dt),$(ndep_to_sminn[c]*dt)\n",
            is_fates ?
                "outputs,lch,roff,dnit,plnt = $(smin_no3_leached[c]*dt),$(smin_no3_runoff[c]*dt),$(f_n2o_nit[c]*dt),$(sminn_to_plant[c]*dt)\n" :
                "outputs,lch,roff,dnit    = $(smin_no3_leached[c]*dt),$(smin_no3_runoff[c]*dt),$(f_n2o_nit[c]*dt)\n"
        )
        error("CNBalanceCheck nitrogen error:\n$msg")
    end

    # --- Gridcell-level balance check (only when not using FATES BGC) ---
    if !use_fates_bgc
        totgrcn = soilbgc_nstate.totn_grc
        c2g_unity!(totgrcn, totcoln, col_data.gridcell, col_data.wtgcell, bounds_c, bounds_g)

        ng = length(grc_data.latdeg)
        grc_ninputs_partial  = fill!(similar(bal.begnb_grc, FT, ng), zero(FT))
        grc_noutputs_partial = fill!(similar(bal.begnb_grc, FT, ng), zero(FT))

        c2g_unity!(grc_ninputs_partial, col_ninputs_partial,
                   col_data.gridcell, col_data.wtgcell, bounds_c, bounds_g)
        c2g_unity!(grc_noutputs_partial, col_noutputs_partial,
                   col_data.gridcell, col_data.wtgcell, bounds_c, bounds_g)

        grc_errnb = fill!(similar(bal.begnb_grc, FT, ng), zero(FT))

        if !isempty(bounds_g)
            _launch!(_nbal_grc_kernel!, grc_endnb, grc_errnb,
                grc_begnb, totgrcn, tot_woodprod_grc, cropprod1_grc,
                grc_ninputs_partial, grc_noutputs_partial,
                dwt_seedn_to_leaf_grc, dwt_seedn_to_deadstem_grc,
                dwt_conv_nflux_grc, product_loss_grc, gru_wood_productn_gain_grc,
                FT(dt), first(bounds_g), last(bounds_g); ndrange = length(grc_endnb))
        end

        if length(bal.errnb_grc) == length(grc_errnb)
            copyto!(bal.errnb_grc, grc_errnb)
        end

        err_found = false
        err_index = 0
        if grc_errnb isa Array
            for g in bounds_g
                if abs(grc_errnb[g]) > bal.nerror
                    err_found = true
                    err_index = g
                end
                if abs(grc_errnb[g]) > bal.nwarning
                    # Same KNOWN open residual as the column warning above.
                    @warn "nbalance warning at g = $g (KNOWN open residual — see docs/CN_BALANCE_STATUS.md)" grc_errnb=grc_errnb[g] grc_endnb=grc_endnb[g] maxlog=3
                end
            end
        end

        if err_found
            g = err_index
            msg = string(
                "gridcell nbalance error  = ", grc_errnb[g], " g=", g, "\n",
                "latdeg, londeg           = ", grc_data.latdeg[g], ",", grc_data.londeg[g], "\n",
                "begnb                    = ", grc_begnb[g], "\n",
                "endnb                    = ", grc_endnb[g], "\n",
                "delta store              = ", grc_endnb[g] - grc_begnb[g], "\n",
                "input mass               = ", grc_ninputs_partial[g] * dt, "\n",
                "output mass              = ", grc_noutputs_partial[g] * dt, "\n",
                "net flux                 = ", (grc_ninputs_partial[g] - grc_noutputs_partial[g]) * dt, "\n",
                "--- Inputs ---\n",
                "grc_ninputs_partial      = ", grc_ninputs_partial[g] * dt, "\n",
                "dwt_seedn_to_leaf_grc    = ", dwt_seedn_to_leaf_grc[g] * dt, "\n",
                "dwt_seedn_to_deadstem_grc = ", dwt_seedn_to_deadstem_grc[g] * dt, "\n",
                "--- Outputs ---\n",
                "grc_noutputs_partial     = ", grc_noutputs_partial[g] * dt, "\n",
                "dwt_conv_nflux_grc       = ", dwt_conv_nflux_grc[g] * dt, "\n",
                "-gru_wood_productn_gain_grc = ", -gru_wood_productn_gain_grc[g] * dt, "\n",
                "product_loss_grc         = ", product_loss_grc[g] * dt
            )
            error("CNBalanceCheck gridcell nitrogen error:\n$msg")
        end
    end

    return nothing
end
