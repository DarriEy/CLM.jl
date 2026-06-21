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
    T = eltype(h2osoi_liq_col)
    dt = T(dtime)

    # Zero the melt flux over the whole do_smb filter first.
    for c in bounds
        mask_do_smb[c] || continue
        qflx_glcice_melt_col[c] = zero(T)
    end

    # Convert meltwater back to ice only for istice columns inside the do_smb filter.
    for j in 1:nlevgrnd
        jj = j + nlevsno
        for c in bounds
            mask_do_smb[c] || continue
            l = col_landunit[c]
            if lun_itype[l] == ISTICE
                liq = h2osoi_liq_col[c, jj]
                if liq > zero(T)   # ice layer with meltwater
                    qflx_glcice_melt_col[c] += liq / dt
                    # convert layer back to pure ice by borrowing ice from below
                    h2osoi_ice_col[c, jj] += liq
                    h2osoi_liq_col[c, jj] = zero(T)
                end
            end
        end
    end

    return nothing
end

# =========================================================================
# compute_surface_mass_balance!
# =========================================================================

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
    T = eltype(qflx_glcice_col)
    # Convert max-days to working precision to avoid integer overflow (matches the
    # Fortran `real(glc_snow_persistence_max_days, r8)` cast before multiplying).
    persistence_threshold = T(glc_snow_persistence_max_days) * T(SECSPDAY)

    # Zero the dynamic water-flux balance term over all columns (handles columns
    # leaving the do_smb / glc_dyn_runoff_routing masks mid-run).
    for c in bounds
        mask_allc[c] || continue
        qflx_glcice_dyn_water_flux_col[c] = zero(T)
    end

    for c in bounds
        mask_do_smb[c] || continue
        l = col_landunit[c]
        g = col_gridcell[c]

        if (snow_persistence_col[c] >= persistence_threshold) || (lun_itype[l] == ISTICE)
            qflx_glcice_frz_col[c] = qflx_snwcp_ice_col[c]
        else
            qflx_glcice_frz_col[c] = zero(T)
        end

        qflx_glcice_col[c] = qflx_glcice_frz_col[c] - qflx_glcice_melt_col[c]

        qflx_glcice_dyn_water_flux_col[c] =
            glc_dyn_runoff_routing_grc[g] * (qflx_glcice_melt_col[c] - qflx_glcice_frz_col[c])
    end

    return nothing
end

# =========================================================================
# adjust_runoff_terms!
# =========================================================================

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
    T = eltype(qflx_qrgwl)

    for c in bounds
        mask_do_smb[c] || continue
        g = col_gridcell[c]

        # Ice melt is added to liquid runoff regardless of dynamic coupling.
        qflx_qrgwl[c] += qflx_glcice_melt_col[c]

        # Capped snow on the dynamically-coupled fraction is owned by the ice sheet.
        qflx_ice_runoff_snwcp[c] -= glc_dyn_runoff_routing_grc[g] * qflx_glcice_frz_col[c]

        # On the uncoupled fraction, remove one unit of ice runoff per unit of melt.
        qflx_ice_runoff_snwcp[c] -=
            (one(T) - glc_dyn_runoff_routing_grc[g]) * qflx_glcice_melt_col[c]
    end

    return nothing
end
