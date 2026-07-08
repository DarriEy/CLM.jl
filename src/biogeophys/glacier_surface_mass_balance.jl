# ==========================================================================
# Ported from: src/biogeophys/GlacierSurfaceMassBalanceMod.F90
# Computes fluxes that are specific to glaciers (surface mass balance).
#
# Public functions:
#   handle_ice_melt!               — Compute ice melt in glacier columns and
#                                    convert liquid meltwater back to ice
#   compute_surface_mass_balance!  — Compute glacier fluxes other than ice melt
#                                    (frz, net glcice, dyn water-flux balance term)
#   adjust_runoff_terms!           — Adjust liquid/ice runoff fluxes for glacier SMB
#
# These three routines correspond to the Fortran `glacier_smb_type` procedures
# HandleIceMelt, ComputeSurfaceMassBalance and AdjustRunoffTerms. They are kept
# separate so they can be sequenced in the driver loop based on what variables
# they depend on and affect (HandleIceMelt right after ice is melted, then SMB
# after qflx_snwcp_ice is computed, then AdjustRunoffTerms after qflx_qrgwl and
# qflx_ice_runoff_snwcp have their initial values).
#
# `glc_dyn_runoff_routing_grc` is the per-gridcell fraction coupled to a dynamic
# ice-sheet model. CLM.jl does not yet port glc2lnd coupling, so callers pass a
# zero vector by default, which reproduces the standalone-CLM (no dynamic ice
# sheet) behaviour exactly.
# ==========================================================================

# =========================================================================
# handle_ice_melt!
# =========================================================================

# One thread per column. Each thread owns column c: it zeros the melt flux and,
# for istice columns, accumulates meltwater into it over the ground layers (a
# per-column RMW done sequentially in-thread → race-free). Byte-identical to the
# original zero-loop + j-outer/c-inner loop: same adds in the same per-column order.
@kernel function _handle_ice_melt_kernel!(qflx_glcice_melt_col, h2osoi_ice_col,
        h2osoi_liq_col, @Const(mask_do_smb), @Const(col_landunit), @Const(lun_itype),
        lo::Int, hi::Int, dt, nlevsno::Int, nlevgrnd::Int, ISTICE_::Int)
    c = @index(Global)
    @inbounds if lo <= c <= hi && mask_do_smb[c]
        z = zero(eltype(qflx_glcice_melt_col))
        qflx_glcice_melt_col[c] = z
        if lun_itype[col_landunit[c]] == ISTICE_
            for j in 1:nlevgrnd
                jj = j + nlevsno
                liq = h2osoi_liq_col[c, jj]
                if liq > z   # ice layer with meltwater
                    qflx_glcice_melt_col[c] += liq / dt
                    h2osoi_ice_col[c, jj] += liq
                    h2osoi_liq_col[c, jj] = z
                end
            end
        end
    end
end

"""
    handle_ice_melt!(h2osoi_liq_col, h2osoi_ice_col, qflx_glcice_melt_col,
        col_landunit, lun_itype, mask_do_smb, bounds, dtime, nlevsno, nlevgrnd)

Compute ice melt in glacier (istice) columns, and convert liquid meltwater
back to ice.

For each `do_smb` column whose landunit is `ISTICE`, any liquid water in a
ground layer is treated as melt: it is accumulated into `qflx_glcice_melt_col`
(positive definite, mm H2O/s) and the layer is converted back to pure ice by
"borrowing" an equivalent ice mass from below the column (`h2osoi_ice += h2osoi_liq;
h2osoi_liq = 0`). This keeps the glacier surface ice-covered while letting the
meltwater run off (the borrowing is reconciled in the runoff/water-balance terms).

`qflx_glcice_melt_col` is first zeroed over the whole `do_smb` filter.

Index convention: Fortran ground layer `j` (1..nlevgrnd) maps to Julia matrix
column `j + nlevsno` in `h2osoi_*_col`.

Ported from `HandleIceMelt` in `GlacierSurfaceMassBalanceMod.F90`.
"""
function handle_ice_melt!(
    h2osoi_liq_col::AbstractMatrix{<:Real},
    h2osoi_ice_col::AbstractMatrix{<:Real},
    qflx_glcice_melt_col::AbstractVector{<:Real},
    col_landunit::AbstractVector{<:Integer},
    lun_itype::AbstractVector{<:Integer},
    mask_do_smb::AbstractVector{Bool},
    bounds::UnitRange{Int},
    dtime::Real,
    nlevsno::Int,
    nlevgrnd::Int,
)
    isempty(bounds) && return nothing
    FT = eltype(h2osoi_liq_col)
    # Index arrays move to the state backend preserving Int; dtime → array precision
    # so no Float64 scalar reaches a Metal kernel. No-ops on the host path.
    cl = _to_backend_like(qflx_glcice_melt_col, FT, col_landunit)
    li = _to_backend_like(qflx_glcice_melt_col, FT, lun_itype)
    _launch!(_handle_ice_melt_kernel!, qflx_glcice_melt_col, h2osoi_ice_col,
        h2osoi_liq_col, mask_do_smb, cl, li,
        first(bounds), last(bounds), FT(dtime), nlevsno, nlevgrnd, Int(ISTICE))
    return nothing
end

# =========================================================================
# compute_surface_mass_balance!
# =========================================================================

# One thread per column, two independent per-column guarded blocks: (a) zero the
# dyn-water-flux over mask_allc, then (b) the do_smb frz / net-glcice / dyn compute.
# The blocks are independent across columns, so per-column "zero then set" reproduces
# the original "zero-all-columns loop, then compute-all-columns loop" exactly.
@kernel function _smb_compute_kernel!(qflx_glcice_dyn_water_flux_col,
        qflx_glcice_col, qflx_glcice_frz_col, @Const(mask_allc), @Const(mask_do_smb),
        @Const(qflx_snwcp_ice_col), @Const(qflx_glcice_melt_col),
        @Const(snow_persistence_col), @Const(glc_dyn_runoff_routing_grc),
        @Const(col_landunit), @Const(col_gridcell), @Const(lun_itype),
        lo::Int, hi::Int, persistence_threshold, ISTICE_::Int)
    c = @index(Global)
    @inbounds if lo <= c <= hi
        z = zero(eltype(qflx_glcice_col))
        if mask_allc[c]
            qflx_glcice_dyn_water_flux_col[c] = z
        end
        if mask_do_smb[c]
            l = col_landunit[c]
            g = col_gridcell[c]
            if (snow_persistence_col[c] >= persistence_threshold) || (lun_itype[l] == ISTICE_)
                qflx_glcice_frz_col[c] = qflx_snwcp_ice_col[c]
            else
                qflx_glcice_frz_col[c] = z
            end
            qflx_glcice_col[c] = qflx_glcice_frz_col[c] - qflx_glcice_melt_col[c]
            qflx_glcice_dyn_water_flux_col[c] =
                glc_dyn_runoff_routing_grc[g] * (qflx_glcice_melt_col[c] - qflx_glcice_frz_col[c])
        end
    end
end

"""
    compute_surface_mass_balance!(qflx_glcice_col, qflx_glcice_frz_col,
        qflx_glcice_dyn_water_flux_col, qflx_snwcp_ice_col, qflx_glcice_melt_col,
        snow_persistence_col, glc_dyn_runoff_routing_grc,
        col_landunit, col_gridcell, lun_itype,
        mask_allc, mask_do_smb, bounds;
        glc_snow_persistence_max_days = 7300)

Compute glacier fluxes other than ice melt: ice growth (`qflx_glcice_frz`), net
glacial ice flux (`qflx_glcice = frz - melt`), and the water-balance correction
flux due to dynamic runoff routing (`qflx_glcice_dyn_water_flux`).

Surface mass balance (ice growth) is generated wherever snow has persisted long
enough (`snow_persistence >= glc_snow_persistence_max_days * SECSPDAY`, i.e.
glacial inception) or the column is already an `ISTICE` landunit; in those cases
the excess snow-capping flux is converted to ice growth, otherwise growth is 0.

`qflx_glcice_dyn_water_flux` is zeroed over the whole-column filter first (in
case the do_smb/glc_dyn membership changes mid-run), then set on do_smb columns.

Should be called after `handle_ice_melt!` and after `qflx_snwcp_ice_col` is computed.

Ported from `ComputeSurfaceMassBalance` in `GlacierSurfaceMassBalanceMod.F90`.
"""
function compute_surface_mass_balance!(
    qflx_glcice_col::AbstractVector{<:Real},
    qflx_glcice_frz_col::AbstractVector{<:Real},
    qflx_glcice_dyn_water_flux_col::AbstractVector{<:Real},
    qflx_snwcp_ice_col::AbstractVector{<:Real},
    qflx_glcice_melt_col::AbstractVector{<:Real},
    snow_persistence_col::AbstractVector{<:Real},
    glc_dyn_runoff_routing_grc::AbstractVector{<:Real},
    col_landunit::AbstractVector{<:Integer},
    col_gridcell::AbstractVector{<:Integer},
    lun_itype::AbstractVector{<:Integer},
    mask_allc::AbstractVector{Bool},
    mask_do_smb::AbstractVector{Bool},
    bounds::UnitRange{Int};
    glc_snow_persistence_max_days::Integer = 7300,
)
    isempty(bounds) && return nothing
    FT = eltype(qflx_glcice_col)
    # Convert max-days to working precision to avoid integer overflow (matches the
    # Fortran `real(glc_snow_persistence_max_days, r8)` cast before multiplying).
    persistence_threshold = FT(glc_snow_persistence_max_days) * FT(SECSPDAY)
    cl = _to_backend_like(qflx_glcice_col, FT, col_landunit)
    cg = _to_backend_like(qflx_glcice_col, FT, col_gridcell)
    li = _to_backend_like(qflx_glcice_col, FT, lun_itype)
    gd = _to_backend_like(qflx_glcice_col, FT, glc_dyn_runoff_routing_grc)
    _launch!(_smb_compute_kernel!, qflx_glcice_dyn_water_flux_col,
        qflx_glcice_col, qflx_glcice_frz_col, mask_allc, mask_do_smb,
        qflx_snwcp_ice_col, qflx_glcice_melt_col, snow_persistence_col, gd,
        cl, cg, li, first(bounds), last(bounds), persistence_threshold, Int(ISTICE))
    return nothing
end

# =========================================================================
# adjust_runoff_terms!
# =========================================================================

# One thread per do_smb column; every write is to the column's own index (no
# cross-column dependency), so the per-column RMWs are race-free.
@kernel function _adjust_runoff_kernel!(qflx_qrgwl, qflx_ice_runoff_snwcp,
        @Const(qflx_glcice_frz_col), @Const(qflx_glcice_melt_col),
        @Const(glc_dyn_runoff_routing_grc), @Const(col_gridcell),
        @Const(mask_do_smb), lo::Int, hi::Int)
    c = @index(Global)
    @inbounds if lo <= c <= hi && mask_do_smb[c]
        one_ = one(eltype(qflx_qrgwl))
        g = col_gridcell[c]
        # Ice melt is added to liquid runoff regardless of dynamic coupling.
        qflx_qrgwl[c] += qflx_glcice_melt_col[c]
        # Capped snow on the dynamically-coupled fraction is owned by the ice sheet.
        qflx_ice_runoff_snwcp[c] -= glc_dyn_runoff_routing_grc[g] * qflx_glcice_frz_col[c]
        # On the uncoupled fraction, remove one unit of ice runoff per unit of melt.
        qflx_ice_runoff_snwcp[c] -=
            (one_ - glc_dyn_runoff_routing_grc[g]) * qflx_glcice_melt_col[c]
    end
end

"""
    adjust_runoff_terms!(qflx_qrgwl, qflx_ice_runoff_snwcp,
        qflx_glcice_frz_col, qflx_glcice_melt_col, glc_dyn_runoff_routing_grc,
        col_gridcell, mask_do_smb, bounds)

Adjust liquid (`qflx_qrgwl`) and ice (`qflx_ice_runoff_snwcp`) runoff fluxes due
to glacier fluxes, over the `do_smb` filter.

- Ice melt is always added to liquid runoff (`qflx_qrgwl += qflx_glcice_melt`),
  whether or not coupled to a dynamic glacier model.
- For the fraction coupled to a dynamic glacier model, the capped snow is owned
  by the ice-sheet model, so it is removed from ice runoff
  (`-= glc_dyn_routing * qflx_glcice_frz`).
- For the uncoupled fraction, one unit of ice runoff is removed per unit of melt
  (`-= (1 - glc_dyn_routing) * qflx_glcice_melt`), correcting the
  accumulation/melt double-counting (conserves both mass and energy; can be
  locally negative but integrates to a reduction of the too-high ice runoff).

Should be called after `compute_surface_mass_balance!` and after `qflx_qrgwl`
and `qflx_ice_runoff_snwcp` have their initial values.

Ported from `AdjustRunoffTerms` in `GlacierSurfaceMassBalanceMod.F90`.
"""
function adjust_runoff_terms!(
    qflx_qrgwl::AbstractVector{<:Real},
    qflx_ice_runoff_snwcp::AbstractVector{<:Real},
    qflx_glcice_frz_col::AbstractVector{<:Real},
    qflx_glcice_melt_col::AbstractVector{<:Real},
    glc_dyn_runoff_routing_grc::AbstractVector{<:Real},
    col_gridcell::AbstractVector{<:Integer},
    mask_do_smb::AbstractVector{Bool},
    bounds::UnitRange{Int},
)
    isempty(bounds) && return nothing
    FT = eltype(qflx_qrgwl)
    cg = _to_backend_like(qflx_qrgwl, FT, col_gridcell)
    gd = _to_backend_like(qflx_qrgwl, FT, glc_dyn_runoff_routing_grc)
    _launch!(_adjust_runoff_kernel!, qflx_qrgwl, qflx_ice_runoff_snwcp,
        qflx_glcice_frz_col, qflx_glcice_melt_col, gd, cg, mask_do_smb,
        first(bounds), last(bounds))
    return nothing
end
