# ==========================================================================
# Ported from: src/biogeochem/CNFireLi2014Mod.F90
# Fire dynamics module for Li et al. (2014) version.
# Based on Li et al. (2012a,b; 2013; 2014).
#
# Public types:
#   CNFireLi2014Data    — Li2014-specific fire data (forcing, GDP, peatland)
#   PftConFireLi2014    — PFT-level fire parameters specific to Li2014
#
# Public functions:
#   need_lightning_and_popdens_li2014  — Returns true
#   p2c!                               — Patch to column area-weighted average
#   cnfire_area_li2014!                — Compute column-level burned area
#   cnfire_fluxes_li2014!              — Fire C/N flux calculations (Li2014)
# ==========================================================================

# ---------------------------------------------------------------------------
# Li2014-specific data
# ---------------------------------------------------------------------------

"""
    CNFireLi2014Data

Li2014-specific fire data: population density, lightning frequency,
GDP, peatland fraction, and crop fire timing.
Ported from member data of `cnfire_li2014_type` in `CNFireLi2014Mod.F90`.
"""
Base.@kwdef mutable struct CNFireLi2014Data{FT<:Real,
                                            V<:AbstractVector{FT},
                                            VI<:AbstractVector{<:Integer}}
    forc_hdm     ::V = Float64[]  # human population density (#/km2) (gridcell)
    forc_lnfm    ::V = Float64[]  # lightning frequency (#/km2/hr) (gridcell)
    gdp_lf_col   ::V = Float64[]  # GDP data (k$/capita) (column)
    peatf_lf_col ::V = Float64[]  # peatland fraction data (0-1) (column)
    abm_lf_col   ::VI = Int[]      # prescribed crop fire month (1-12) (column)
end
CNFireLi2014Data{FT}(; kwargs...) where {FT<:Real} =
    CNFireLi2014Data{FT, Vector{FT}, Vector{Int}}(; kwargs...)
Adapt.@adapt_structure CNFireLi2014Data

# ---------------------------------------------------------------------------
# PFT-level fire parameters specific to Li2014
# ---------------------------------------------------------------------------

"""
    PftConFireLi2014

PFT-level fire parameters specific to the Li2014 fire module.
Contains fire spread rate and fire duration by PFT.
"""
Base.@kwdef mutable struct PftConFireLi2014{FT<:Real, V<:AbstractVector{FT}}
    fsr_pft ::V = Float64[]  # fire spread rate by PFT (km/hr)
    fd_pft  ::V = Float64[]  # fire duration by PFT (hr)
end
PftConFireLi2014{FT}(; kwargs...) where {FT<:Real} =
    PftConFireLi2014{FT, Vector{FT}}(; kwargs...)
Adapt.@adapt_structure PftConFireLi2014

# ---------------------------------------------------------------------------
# need_lightning_and_popdens_li2014
# ---------------------------------------------------------------------------

"""
    need_lightning_and_popdens_li2014()

Returns `true` — the Li2014 fire model requires lightning and population
density data.
Ported from `need_lightning_and_popdens` in `CNFireLi2014Mod.F90`.
"""
need_lightning_and_popdens_li2014() = true

# ---------------------------------------------------------------------------
# p2c! — patch to column area-weighted average
# ---------------------------------------------------------------------------

"""
    p2c!(col_out, patch_in, patch, mask_soilc, bounds_c, bounds_p)

Patch-to-column area-weighted average.
Computes: col_out[c] = Σ_p (patch_in[p] * patch.wtcol[p]) for all p on column c.
Ported from `p2c` in `subgridAveMod.F90`.
"""
# Zero col_out[c] for every masked column in [cmin, cmax].
@kernel function _fire_p2c_zero_kernel!(col_out, @Const(mask_soilc), cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && mask_soilc[c]
        col_out[c] = 0.0
    end
end

# Per-patch scatter: col_out[column[p]] += patch_in[p] * wtcol[p] for masked
# patches in [pmin, pmax] whose column c lies in [cmin, cmax] and is masked.
# _scatter_add! is atomic on GPU; on the KA CPU backend it iterates patches in
# ascending order so the accumulation is byte-identical to the host p-loop.
@kernel function _fire_p2c_scatter_kernel!(col_out, @Const(patch_in), @Const(column),
                                           @Const(wtcol), @Const(mask_soilc),
                                           cmin::Int, cmax::Int, pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax
        c = column[p]
        if cmin <= c <= cmax && mask_soilc[c]
            _scatter_add!(col_out, c, patch_in[p] * wtcol[p])
        end
    end
end

function p2c!(
    col_out::AbstractVector{<:Real},
    patch_in::AbstractVector{<:Real},
    patch_d::PatchData,
    mask_soilc::AbstractVector{Bool},
    bounds_c::UnitRange{Int},
    bounds_p::UnitRange{Int}
)
    _launch!(_fire_p2c_zero_kernel!, col_out, mask_soilc,
             first(bounds_c), last(bounds_c))
    _launch!(_fire_p2c_scatter_kernel!, col_out, patch_in, patch_d.column,
             patch_d.wtcol, mask_soilc,
             first(bounds_c), last(bounds_c), first(bounds_p), last(bounds_p);
             ndrange = length(patch_d.column))
    return nothing
end

# ---------------------------------------------------------------------------
# Kernelized per-column peatland fire (Li2014)
# Each column is fully independent: baf_peatf[c] depends only on column c's
# own inputs (its gridcell g). No accumulation / loop-carried deps.
# ---------------------------------------------------------------------------
@kernel function _fire_peatland_kernel!(
    baf_peatf,
    @Const(mask_soilc), @Const(gridcell), @Const(latdeg),
    @Const(prec60_col), @Const(peatf_lf), @Const(fsat),
    @Const(wf2), @Const(tsoi17),
    borealat, non_boreal_peatfire_c, boreal_peatfire_c,
    nonborpeat_fire_precip_denom, borpeat_fire_soilmoist_denom,
    secsphr, secspday, tfrz, rpi
)
    c = @index(Global)
    @inbounds if mask_soilc[c]
        g = gridcell[c]
        if latdeg[g] < borealat
            baf_peatf[c] = non_boreal_peatfire_c / secsphr *
                smooth_max(0.0, smooth_min(1.0,
                    (4.0 - prec60_col[c] * secspday / nonborpeat_fire_precip_denom) /
                    4.0))^2 * peatf_lf[c] * (1.0 - fsat[c])
        else
            baf_peatf[c] = boreal_peatfire_c / secsphr *
                exp(-rpi * (smooth_max(wf2[c], 0.0) / borpeat_fire_soilmoist_denom)) *
                smooth_max(0.0, smooth_min(1.0, (tsoi17[c] - tfrz) / 10.0)) * peatf_lf[c] *
                (1.0 - fsat[c])
        end
    end
end

function fire_peatland!(
    baf_peatf, mask_soilc, gridcell, latdeg, prec60_col, peatf_lf, fsat,
    wf2, tsoi17, borealat, non_boreal_peatfire_c, boreal_peatfire_c,
    nonborpeat_fire_precip_denom, borpeat_fire_soilmoist_denom,
    secsphr, secspday, tfrz, rpi
)
    _launch!(_fire_peatland_kernel!, baf_peatf,
        mask_soilc, gridcell, latdeg, prec60_col, peatf_lf, fsat,
        wf2, tsoi17, borealat, non_boreal_peatfire_c, boreal_peatfire_c,
        nonborpeat_fire_precip_denom, borpeat_fire_soilmoist_denom,
        secsphr, secspday, tfrz, rpi;
        ndrange = length(baf_peatf))
    return nothing
end

# ---------------------------------------------------------------------------
# cnfire_area_li2014! — Compute column-level burned area
# ---------------------------------------------------------------------------

# ===========================================================================
# KernelAbstractions kernels for the cnfire_area_li2014! driver sections.
# Each kernel replaces one host loop section; they are launched in order from
# the rewritten driver below. Patch→column accumulations use _scatter_add!
# (atomic on GPU, ascending-p plain += on the KA CPU backend → byte-identical).
# Accumulators are zeroed by a preceding zero-kernel, exactly as the host did.
# ===========================================================================

# --- Section: zero cropf_col / lfwt (per-column) ---------------------------
@kernel function _firea_zero_cropf_lfwt_kernel!(cropf_col, lfwt, @Const(mask_soilc))
    c = @index(Global)
    @inbounds if mask_soilc[c]
        cropf_col[c] = 0.0
        lfwt[c]      = 0.0
    end
end

# --- Section: accumulate cropf_col / lfwt (per-patch scatter) --------------
@kernel function _firea_cropf_lfwt_kernel!(cropf_col, lfwt, @Const(mask_soilp),
                                           @Const(column), @Const(itype), @Const(wtcol),
                                           nc4_grass::Int, ndllf_evr_tmp_tree::Int)
    p = @index(Global)
    @inbounds if mask_soilp[p]
        c = column[p]
        if itype[p] > nc4_grass
            _scatter_add!(cropf_col, c, wtcol[p])
        end
        if itype[p] >= ndllf_evr_tmp_tree && itype[p] <= nc4_grass
            _scatter_add!(lfwt, c, wtcol[p])
        end
    end
end

# --- Section: zero fuelc_crop (per-column) ---------------------------------
@kernel function _firea_zero_fuelc_crop_kernel!(fuelc_crop, @Const(mask_soilc))
    c = @index(Global)
    @inbounds if mask_soilc[c]
        fuelc_crop[c] = 0.0
    end
end

# --- Section: accumulate crop fuel (per-patch scatter, /cropf_col) ---------
@kernel function _firea_fuelc_crop_kernel!(fuelc_crop, @Const(mask_soilp),
                                           @Const(column), @Const(itype), @Const(wtcol),
                                           @Const(leafc), @Const(leafc_storage),
                                           @Const(leafc_xfer), @Const(leafc_col),
                                           @Const(totlitc), @Const(cropf_col),
                                           nc4_grass::Int)
    p = @index(Global)
    @inbounds if mask_soilp[p]
        c = column[p]
        if itype[p] > nc4_grass && wtcol[p] > 0.0 && leafc_col[c] > 0.0
            contrib = (leafc[p] + leafc_storage[p] + leafc_xfer[p]) *
                          wtcol[p] / cropf_col[c] +
                      totlitc[c] * leafc[p] / leafc_col[c] *
                          wtcol[p] / cropf_col[c]
            _scatter_add!(fuelc_crop, c, contrib)
        end
    end
end

# --- Section: initialize noncrop column variables (per-column) -------------
@kernel function _firea_init_noncrop_kernel!(fsr_col, fd_col, rootc_col, lgdp_col,
                                             lgdp1_col, lpop_col, btran_col, wtlf,
                                             trotr1_col, trotr2_col, dtrotr_col,
                                             @Const(mask_soilc),
                                             transient_landcover::Bool)
    c = @index(Global)
    @inbounds if mask_soilc[c]
        fsr_col[c]    = 0.0
        fd_col[c]     = 0.0
        rootc_col[c]  = 0.0
        lgdp_col[c]   = 0.0
        lgdp1_col[c]  = 0.0
        lpop_col[c]   = 0.0
        btran_col[c]  = 0.0
        wtlf[c]       = 0.0
        trotr1_col[c] = 0.0
        trotr2_col[c] = 0.0
        if transient_landcover
            dtrotr_col[c] = 0.0
        end
    end
end

# --- Section: accumulate btran_col / wtlf (per-patch scatter) --------------
@kernel function _firea_btran_wtlf_kernel!(btran_col, wtlf, @Const(mask_exposedveg),
                                           @Const(column), @Const(itype), @Const(wtcol),
                                           @Const(btran2), @Const(cropf_col),
                                           nc3crop::Int)
    p = @index(Global)
    @inbounds if mask_exposedveg[p]
        c = column[p]
        if itype[p] < nc3crop && cropf_col[c] < 1.0
            _scatter_add!(btran_col, c, btran2[p] * wtcol[p])
            _scatter_add!(wtlf, c, wtcol[p])
        end
    end
end

# --- Section: main noncrop per-patch column-var accumulation (scatter) -----
@kernel function _firea_noncrop_main_kernel!(trotr1_col, trotr2_col, dtrotr_col,
        rootc_col, fsr_col, lgdp_col, lgdp1_col, lpop_col, fd_col,
        @Const(mask_soilp), @Const(column), @Const(gridcell), @Const(itype),
        @Const(wtcol), @Const(cropf_col), @Const(lfwt), @Const(dwt_smoothed),
        @Const(frootc), @Const(frootc_storage), @Const(frootc_xfer),
        @Const(deadcrootc), @Const(deadcrootc_storage), @Const(deadcrootc_xfer),
        @Const(livecrootc), @Const(livecrootc_storage), @Const(livecrootc_xfer),
        @Const(fsr_pft), @Const(fd_pft), @Const(gdp_lf), @Const(forc_hdm),
        nc3crop::Int, nbrdlf_evr_trp_tree::Int, nbrdlf_dcd_trp_tree::Int,
        nbrdlf_evr_shrub::Int, noveg::Int, occur_hi_gdp_tree, rpi, secsphr,
        transient_landcover::Bool)
    p = @index(Global)
    @inbounds if mask_soilp[p]
        c = column[p]
        g = gridcell[c]
        if itype[p] < nc3crop && cropf_col[c] < 1.0
            if itype[p] == nbrdlf_evr_trp_tree && wtcol[p] > 0.0
                _scatter_add!(trotr1_col, c, wtcol[p])
            end
            if itype[p] == nbrdlf_dcd_trp_tree && wtcol[p] > 0.0
                _scatter_add!(trotr2_col, c, wtcol[p])
            end
            if transient_landcover
                if itype[p] == nbrdlf_evr_trp_tree || itype[p] == nbrdlf_dcd_trp_tree
                    if dwt_smoothed[p] < 0.0
                        _scatter_add!(dtrotr_col, c, -dwt_smoothed[p])
                    end
                end
            end
            _scatter_add!(rootc_col, c, (frootc[p] + frootc_storage[p] +
                frootc_xfer[p] + deadcrootc[p] +
                deadcrootc_storage[p] + deadcrootc_xfer[p] +
                livecrootc[p] + livecrootc_storage[p] +
                livecrootc_xfer[p]) * wtcol[p])
            _scatter_add!(fsr_col, c, fsr_pft[itype[p] + 1] * wtcol[p] / (1.0 - cropf_col[c]))
            if lfwt[c] != 0.0
                hdmlf = forc_hdm[g]
                if hdmlf > 0.1
                    if itype[p] != noveg
                        if itype[p] >= nbrdlf_evr_shrub
                            _scatter_add!(lgdp_col, c, (0.1 + 0.9 *
                                exp(-1.0 * rpi * (gdp_lf[c] / 8.0)^0.5)) * wtcol[p] /
                                (1.0 - cropf_col[c]))
                            _scatter_add!(lgdp1_col, c, (0.2 + 0.8 *
                                exp(-1.0 * rpi * (gdp_lf[c] / 7.0))) * wtcol[p] / lfwt[c])
                            _scatter_add!(lpop_col, c, (0.2 + 0.8 *
                                exp(-1.0 * rpi * (hdmlf / 450.0)^0.5)) * wtcol[p] / lfwt[c])
                        else
                            if gdp_lf[c] > 20.0
                                _scatter_add!(lgdp_col, c, occur_hi_gdp_tree *
                                    wtcol[p] / (1.0 - cropf_col[c]))
                            else
                                _scatter_add!(lgdp_col, c, wtcol[p] / (1.0 - cropf_col[c]))
                            end
                            if gdp_lf[c] > 20.0
                                _scatter_add!(lgdp1_col, c, 0.62 * wtcol[p] / lfwt[c])
                            else
                                if gdp_lf[c] > 8.0
                                    _scatter_add!(lgdp1_col, c, 0.83 * wtcol[p] / lfwt[c])
                                else
                                    _scatter_add!(lgdp1_col, c, wtcol[p] / lfwt[c])
                                end
                            end
                            _scatter_add!(lpop_col, c, (0.4 + 0.6 *
                                exp(-1.0 * rpi * (hdmlf / 125.0))) * wtcol[p] / lfwt[c])
                        end
                    end
                else
                    _scatter_add!(lgdp_col, c, wtcol[p] / (1.0 - cropf_col[c]))
                    _scatter_add!(lgdp1_col, c, wtcol[p] / lfwt[c])
                    _scatter_add!(lpop_col, c, wtcol[p] / lfwt[c])
                end
            end
            _scatter_add!(fd_col, c, fd_pft[itype[p]] * wtcol[p] * secsphr / (1.0 - cropf_col[c]))
        end
    end
end

# --- Section: transient-landcover lfc (per-column, date-gated) -------------
@kernel function _firea_lfc_kernel!(lfc, @Const(mask_soilc), @Const(dtrotr_col),
                                    date_is_jan1_0::Bool, date_is_jan1_dt::Bool,
                                    dayspyr, secspday, dt)
    c = @index(Global)
    @inbounds if mask_soilc[c]
        if dtrotr_col[c] > 0.0
            if date_is_jan1_0
                lfc[c] = 0.0
            end
            if date_is_jan1_dt
                lfc[c] = dtrotr_col[c] * dayspyr * secspday / dt
            end
        else
            lfc[c] = 0.0
        end
    end
end

# --- Section: zero baf_crop (per-column) -----------------------------------
@kernel function _firea_zero_baf_crop_kernel!(baf_crop, @Const(mask_soilc))
    c = @index(Global)
    @inbounds if mask_soilc[c]
        baf_crop[c] = 0.0
    end
end

# --- Section: burndate init (per-patch, date-gated) ------------------------
@kernel function _firea_burndate_init_kernel!(burndate, @Const(mask_soilp),
                                              date_is_jan1_0::Bool)
    p = @index(Global)
    @inbounds if mask_soilp[p]
        if date_is_jan1_0
            burndate[p] = 10000
        end
    end
end

# --- Section: crop burned-area fraction (per-patch scatter + burndate) -----
@kernel function _firea_crop_fire_kernel!(baf_crop, burndate, @Const(mask_soilp),
        @Const(column), @Const(gridcell), @Const(itype), @Const(wtcol),
        @Const(forc_t), @Const(abm_lf), @Const(forc_rain), @Const(forc_snow),
        @Const(fuelc_crop), @Const(gdp_lf), @Const(forc_hdm),
        nc4_grass::Int, kmo::Int, kda::Int, tfrz, rpi, lfuel, ufuel,
        cropfire_a1, secsphr)
    p = @index(Global)
    @inbounds if mask_soilp[p]
        c = column[p]
        g = gridcell[c]
        if forc_t[c] >= tfrz && itype[p] > nc4_grass &&
           kmo == abm_lf[c] && forc_rain[c] + forc_snow[c] == 0.0 &&
           burndate[p] >= 999 && wtcol[p] > 0.0
            hdmlf = forc_hdm[g]
            fhd = 0.04 + 0.96 * exp(-1.0 * rpi * (hdmlf / 350.0)^0.5)
            fgdp = 0.01 + 0.99 * exp(-1.0 * rpi * (gdp_lf[c] / 10.0))
            fb = smooth_max(0.0, smooth_min(1.0, (fuelc_crop[c] - lfuel) / (ufuel - lfuel)))
            _scatter_add!(baf_crop, c, cropfire_a1 / secsphr * fb * fhd * fgdp * wtcol[p])
            if fb * fhd * fgdp * wtcol[p] > 0.0
                burndate[p] = kda
            end
        end
    end
end

# --- Section: main per-column fire-spread loop -----------------------------
@kernel function _firea_fire_spread_kernel!(nfire, farea_burned, fuelc, fbac, fbac1,
        @Const(mask_soilc), @Const(gridcell), @Const(forc_hdm), @Const(cropf_col),
        @Const(trotr1_col), @Const(trotr2_col), @Const(baf_crop), @Const(baf_peatf),
        @Const(totlitc), @Const(totvegc), @Const(rootc_col), @Const(fuelc_crop),
        @Const(decomp_cpools_vr), @Const(dzsoi_decomp), @Const(wf), @Const(forc_rh),
        @Const(forc_t), @Const(forc_lnfm), @Const(lat), @Const(forc_wind),
        @Const(lgdp_col), @Const(lgdp1_col), @Const(lpop_col), @Const(fsr_col),
        @Const(fd_col), @Const(btran_col), @Const(wtlf), @Const(dtrotr_col),
        @Const(prec60_col), @Const(prec10_col), @Const(forc_rain), @Const(forc_snow),
        @Const(lfc),
        i_cwd::Int, nlevdecomp::Int, lfuel, ufuel, rh_low, rh_hgh, bt_min, bt_max,
        pot_hmn_ign_counts_alpha, ignition_efficiency, g0, cli_scale,
        defo_fire_precip_thresh_bet, defo_fire_precip_thresh_bdt,
        tfrz, rpi, secsphr, secspday, dayspyr, dt,
        date_is_jan1_0::Bool, transient_landcover::Bool)
    c = @index(Global)
    @inbounds if mask_soilc[c]
        g = gridcell[c]
        hdmlf = forc_hdm[g]
        nfire[c] = 0.0

        if cropf_col[c] < 1.0
            if trotr1_col[c] + trotr2_col[c] > 0.6
                farea_burned[c] = smooth_min(1.0, baf_crop[c] + baf_peatf[c])
            else
                fc = totlitc[c] + totvegc[c] - rootc_col[c] - fuelc_crop[c] * cropf_col[c]
                for j in 1:nlevdecomp
                    fc = fc + decomp_cpools_vr[c, j, i_cwd] * dzsoi_decomp[j]
                end
                fc = fc / (1.0 - cropf_col[c])
                fuelc[c] = fc
                fb       = smooth_max(0.0, smooth_min(1.0, (fc - lfuel) / (ufuel - lfuel)))
                m        = smooth_max(0.0, wf[c])
                fire_m   = exp(-rpi * (m / 0.69)^2) * (1.0 - smooth_max(0.0,
                    smooth_min(1.0, (forc_rh[g] - rh_low) / (rh_hgh - rh_low)))) *
                    smooth_min(1.0, exp(rpi * (forc_t[c] - tfrz) / 10.0))
                lh       = pot_hmn_ign_counts_alpha * 6.8 * hdmlf^0.43 / 30.0 / 24.0
                fs       = 1.0 - (0.01 + 0.98 * exp(-0.025 * hdmlf))
                ig       = (lh + forc_lnfm[g] /
                    (5.16 + 2.16 * cos(3.0 * lat[g])) *
                    ignition_efficiency) * (1.0 - fs) * (1.0 - cropf_col[c])
                nfire[c] = ig / secsphr * fb * fire_m * lgdp_col[c]
                Lb_lf    = 1.0 + 10.0 * (1.0 - exp(-0.06 * forc_wind[g]))
                if wtlf[c] > 0.0
                    spread_m = (1.0 - smooth_max(0.0, smooth_min(1.0,
                        (btran_col[c] / wtlf[c] - bt_min) /
                        (bt_max - bt_min)))) * (1.0 - smooth_max(0.0,
                        smooth_min(1.0, (forc_rh[g] - rh_low) / (rh_hgh - rh_low))))
                else
                    spread_m = 0.0
                end
                farea_burned[c] = smooth_min(1.0,
                    (g0 * spread_m * fsr_col[c] *
                    fd_col[c] / 1000.0)^2 * lgdp1_col[c] *
                    lpop_col[c] * nfire[c] * rpi * Lb_lf +
                    baf_crop[c] + baf_peatf[c])
            end

            if transient_landcover
                if trotr1_col[c] + trotr2_col[c] > 0.6
                    if date_is_jan1_0 || dtrotr_col[c] <= 0.0
                        fbac1[c]        = 0.0
                        farea_burned[c] = baf_crop[c] + baf_peatf[c]
                    else
                        cri = (defo_fire_precip_thresh_bet * trotr1_col[c] +
                               defo_fire_precip_thresh_bdt * trotr2_col[c]) /
                              (trotr1_col[c] + trotr2_col[c])
                        cli = (smooth_max(0.0, smooth_min(1.0, (cri - prec60_col[c] * secspday) / cri))^0.5) *
                              (smooth_max(0.0, smooth_min(1.0, (cri - prec10_col[c] * secspday) / cri))^0.5) *
                              smooth_max(0.0005, smooth_min(1.0, 19.0 * dtrotr_col[c] * dayspyr * secspday / dt - 0.001)) *
                              smooth_max(0.0, smooth_min(1.0, (0.25 - (forc_rain[c] + forc_snow[c]) * secsphr) / 0.25))
                        farea_burned[c] = cli * (cli_scale / secspday) + baf_crop[c] + baf_peatf[c]
                        fbac1[c] = smooth_max(0.0, cli * (cli_scale / secspday) - 2.0 * lfc[c] / dt)
                    end
                    fbac[c] = fbac1[c] + baf_crop[c] + baf_peatf[c]
                else
                    fbac[c] = farea_burned[c]
                end
            end
        else
            farea_burned[c] = smooth_min(1.0, baf_crop[c] + baf_peatf[c])
        end
    end
end

"""
    cnfire_area_li2014!(fire_li2014, pftcon_li2014, fire_data, cnfire_const,
        cnfire_params, pftcon, mask_soilc, mask_soilp, mask_exposedveg,
        mask_noexposedveg, bounds_c, bounds_p, patch, col, grc,
        soilstate, h2osoi_vol_col, cnveg_state, cnveg_cs,
        decomp_cascade_con, totlitc_col, decomp_cpools_vr_col,
        t_soi17cm_col; kwargs...)

Compute column-level burned area fraction for the Li et al. (2014)
fire model. Includes cropland fire, peatland fire, deforestation fire,
and natural fire components.

Ported from `CNFireArea` in `CNFireLi2014Mod.F90`.
"""
function cnfire_area_li2014!(
    # Li2014-specific data
    fire_li2014::CNFireLi2014Data,
    pftcon_li2014::PftConFireLi2014,
    # Fire base data & constants
    fire_data::CNFireBaseData,
    cnfire_const::CNFireConstData,
    cnfire_params::CNFireParams,
    pftcon::PftConFireBase,
    # Masks and bounds
    mask_soilc::BitVector,
    mask_soilp::BitVector,
    mask_exposedveg::BitVector,
    mask_noexposedveg::BitVector,
    bounds_c::UnitRange{Int},
    bounds_p::UnitRange{Int},
    # Core type data
    patch::PatchData,
    col::ColumnData,
    grc::GridcellData,
    # Soil state (for root wetness)
    soilstate::SoilStateData,
    h2osoi_vol_col::Matrix{<:Real},
    # CN Veg State and Carbon State
    cnveg_state::CNVegStateData,
    cnveg_cs::CNVegCarbonStateData,
    # Decomposition data
    decomp_cascade_con::DecompCascadeConData,
    # Input arrays
    totlitc_col::Vector{<:Real},
    decomp_cpools_vr_col::Array{<:Real,3},
    t_soi17cm_col::Vector{<:Real};
    # Atmospheric forcing (keyword args)
    forc_rh_grc::Vector{<:Real} = Float64[],
    forc_wind_grc::Vector{<:Real} = Float64[],
    forc_t_col::Vector{<:Real} = Float64[],
    forc_rain_col::Vector{<:Real} = Float64[],
    forc_snow_col::Vector{<:Real} = Float64[],
    prec60_patch::Vector{<:Real} = Float64[],
    prec10_patch::Vector{<:Real} = Float64[],
    # Saturated runoff
    fsat_col::Vector{<:Real} = Float64[],
    # Water diagnostic
    wf_col::Vector{<:Real} = Float64[],
    wf2_col::Vector{<:Real} = Float64[],
    # Time/date
    dt::Real = 1800.0,
    dayspyr::Real = 365.0,
    kmo::Int = 1,
    kda::Int = 1,
    mcsec::Int = 0,
    nstep::Int = 1,
    nlevgrnd::Int = 10,
    nlevdecomp::Int = 1,
    ndecomp_pools::Int = 7,
    # PFT indices
    nc4_grass::Int = 14,
    nc3crop::Int = 15,
    ndllf_evr_tmp_tree::Int = 1,
    nbrdlf_evr_trp_tree::Int = 4,
    nbrdlf_dcd_trp_tree::Int = 6,
    nbrdlf_evr_shrub::Int = 9,
    noveg::Int = 0,
    # Control flags
    transient_landcover::Bool = false,
    # Soil suction function
    soil_suction_fn::Function = default_soil_suction
)
    # --- Aliases (matching Fortran associate block) ---
    totlitc          = totlitc_col
    decomp_cpools_vr = decomp_cpools_vr_col
    tsoi17           = t_soi17cm_col
    lfuel            = cnfire_const.lfuel
    ufuel            = cnfire_const.ufuel
    rh_hgh           = cnfire_const.rh_hgh
    rh_low           = cnfire_const.rh_low
    bt_min           = cnfire_const.bt_min
    bt_max           = cnfire_const.bt_max
    cli_scale        = cnfire_const.cli_scale
    cropfire_a1      = cnfire_const.cropfire_a1
    non_boreal_peatfire_c    = cnfire_const.non_boreal_peatfire_c
    pot_hmn_ign_counts_alpha = cnfire_const.pot_hmn_ign_counts_alpha
    boreal_peatfire_c        = cnfire_const.boreal_peatfire_c
    defo_fire_precip_thresh_bet  = cnfire_const.defo_fire_precip_thresh_bet
    defo_fire_precip_thresh_bdt  = cnfire_const.defo_fire_precip_thresh_bdt
    borpeat_fire_soilmoist_denom = cnfire_const.borpeat_fire_soilmoist_denom
    nonborpeat_fire_precip_denom = cnfire_const.nonborpeat_fire_precip_denom

    fsr_pft          = pftcon_li2014.fsr_pft
    fd_pft           = pftcon_li2014.fd_pft

    btran2           = fire_data.btran2_patch
    gdp_lf           = fire_li2014.gdp_lf_col
    peatf_lf         = fire_li2014.peatf_lf_col
    abm_lf           = fire_li2014.abm_lf_col
    is_cwd           = decomp_cascade_con.is_cwd

    forc_rh          = forc_rh_grc
    forc_wind        = forc_wind_grc
    forc_t           = forc_t_col
    forc_rain        = forc_rain_col
    forc_snow        = forc_snow_col
    prec60           = prec60_patch
    prec10           = prec10_patch
    fsat             = fsat_col
    wf               = wf_col
    wf2              = wf2_col

    dwt_smoothed     = cnveg_state.dwt_smoothed_patch
    cropf_col        = cnveg_state.cropf_col
    baf_crop         = cnveg_state.baf_crop_col
    baf_peatf        = cnveg_state.baf_peatf_col
    burndate         = cnveg_state.burndate_patch
    fbac             = cnveg_state.fbac_col
    fbac1            = cnveg_state.fbac1_col
    farea_burned     = cnveg_state.farea_burned_col
    nfire            = cnveg_state.nfire_col
    fsr_col          = cnveg_state.fsr_col
    fd_col           = cnveg_state.fd_col
    lgdp_col         = cnveg_state.lgdp_col
    lgdp1_col        = cnveg_state.lgdp1_col
    lpop_col         = cnveg_state.lpop_col
    lfwt             = cnveg_state.lfwt_col
    trotr1_col       = cnveg_state.trotr1_col
    trotr2_col       = cnveg_state.trotr2_col
    dtrotr_col       = cnveg_state.dtrotr_col
    lfc              = cnveg_state.lfc_col
    wtlf             = cnveg_state.wtlf_col

    totvegc          = cnveg_cs.totvegc_col
    rootc_col        = cnveg_cs.rootc_col
    leafc_col        = cnveg_cs.leafc_col
    fuelc            = cnveg_cs.fuelc_col
    fuelc_crop       = cnveg_cs.fuelc_crop_col

    deadcrootc       = cnveg_cs.deadcrootc_patch
    deadcrootc_storage = cnveg_cs.deadcrootc_storage_patch
    deadcrootc_xfer  = cnveg_cs.deadcrootc_xfer_patch
    frootc           = cnveg_cs.frootc_patch
    frootc_storage   = cnveg_cs.frootc_storage_patch
    frootc_xfer      = cnveg_cs.frootc_xfer_patch
    livecrootc       = cnveg_cs.livecrootc_patch
    livecrootc_storage = cnveg_cs.livecrootc_storage_patch
    livecrootc_xfer  = cnveg_cs.livecrootc_xfer_patch
    leafc            = cnveg_cs.leafc_patch
    leafc_storage    = cnveg_cs.leafc_storage_patch
    leafc_xfer       = cnveg_cs.leafc_xfer_patch

    # --- Local arrays ---
    nc = length(bounds_c)
    FT = eltype(totvegc)
    btran_col = zeros(FT, last(bounds_c))

    # Temporary arrays for p2c results
    prec60_col = zeros(FT, last(bounds_c))
    prec10_col = zeros(FT, last(bounds_c))

    # --- Patch to column averaging ---
    p2c!(prec10_col, prec10, patch, mask_soilc, bounds_c, bounds_p)
    p2c!(prec60_col, prec60, patch, mask_soilc, bounds_c, bounds_p)
    p2c!(leafc_col,  leafc,  patch, mask_soilc, bounds_c, bounds_p)

    # --- On first time-step, set area burned to zero and exit ---
    if nstep == 0
        for c in bounds_c
            mask_soilc[c] || continue
            farea_burned[c] = 0.0
            baf_crop[c]     = 0.0
            baf_peatf[c]    = 0.0
            fbac[c]         = 0.0
            fbac1[c]        = 0.0
            cropf_col[c]    = 0.0
        end
        return nothing
    end

    # Host-resolved uniform date predicates (uniform across all columns/patches).
    date_is_jan1_0  = (kmo == 1 && kda == 1 && mcsec == 0)
    date_is_jan1_dt = (kmo == 1 && kda == 1 && mcsec == dt)

    # --- Calculate cropf_col and lfwt ---
    _launch!(_firea_zero_cropf_lfwt_kernel!, cropf_col, lfwt, mask_soilc)
    _launch!(_firea_cropf_lfwt_kernel!, cropf_col, lfwt, mask_soilp,
             patch.column, patch.itype, patch.wtcol,
             nc4_grass, ndllf_evr_tmp_tree; ndrange = length(mask_soilp))

    # --- Calculate crop fuel ---
    _launch!(_firea_zero_fuelc_crop_kernel!, fuelc_crop, mask_soilc)
    _launch!(_firea_fuelc_crop_kernel!, fuelc_crop, mask_soilp,
             patch.column, patch.itype, patch.wtcol,
             leafc, leafc_storage, leafc_xfer, leafc_col, totlitc, cropf_col,
             nc4_grass; ndrange = length(mask_soilp))

    # --- Initialize noncrop column variables ---
    _launch!(_firea_init_noncrop_kernel!, fsr_col, fd_col, rootc_col, lgdp_col,
             lgdp1_col, lpop_col, btran_col, wtlf, trotr1_col, trotr2_col,
             dtrotr_col, mask_soilc, transient_landcover)

    # --- Calculate root wetness (btran2) ---
    cnfire_calc_fire_root_wetness_li2014!(
        fire_data,
        mask_exposedveg,
        mask_noexposedveg,
        bounds_p,
        pftcon,
        patch,
        soilstate,
        h2osoi_vol_col,
        nlevgrnd;
        soil_suction_fn = soil_suction_fn
    )

    # --- Accumulate btran_col and wtlf from exposed vegetation ---
    _launch!(_firea_btran_wtlf_kernel!, btran_col, wtlf, mask_exposedveg,
             patch.column, patch.itype, patch.wtcol, btran2, cropf_col,
             nc3crop; ndrange = length(mask_exposedveg))

    # --- Main patch loop for noncrop column variables ---
    _launch!(_firea_noncrop_main_kernel!, trotr1_col, trotr2_col, dtrotr_col,
             rootc_col, fsr_col, lgdp_col, lgdp1_col, lpop_col, fd_col,
             mask_soilp, patch.column, col.gridcell, patch.itype, patch.wtcol,
             cropf_col, lfwt, dwt_smoothed,
             frootc, frootc_storage, frootc_xfer,
             deadcrootc, deadcrootc_storage, deadcrootc_xfer,
             livecrootc, livecrootc_storage, livecrootc_xfer,
             fsr_pft, fd_pft, gdp_lf, fire_li2014.forc_hdm,
             nc3crop, nbrdlf_evr_trp_tree, nbrdlf_dcd_trp_tree,
             nbrdlf_evr_shrub, noveg, cnfire_const.occur_hi_gdp_tree, RPI, SECSPHR,
             transient_landcover; ndrange = length(mask_soilp))

    # --- Estimate annual decreased fractional coverage of BET+BDT ---
    if transient_landcover
        _launch!(_firea_lfc_kernel!, lfc, mask_soilc, dtrotr_col,
                 date_is_jan1_0, date_is_jan1_dt, dayspyr, SECSPDAY, dt)
    end

    # --- Calculate burned area fraction in cropland ---
    _launch!(_firea_zero_baf_crop_kernel!, baf_crop, mask_soilc)

    _launch!(_firea_burndate_init_kernel!, burndate, mask_soilp,
             date_is_jan1_0; ndrange = length(mask_soilp))

    _launch!(_firea_crop_fire_kernel!, baf_crop, burndate, mask_soilp,
             patch.column, col.gridcell, patch.itype, patch.wtcol,
             forc_t, abm_lf, forc_rain, forc_snow, fuelc_crop, gdp_lf,
             fire_li2014.forc_hdm,
             nc4_grass, kmo, kda, TFRZ, RPI, lfuel, ufuel, cropfire_a1, SECSPHR;
             ndrange = length(mask_soilp))

    # --- Calculate peatland fire (kernelized; per-column independent) ---
    fire_peatland!(
        baf_peatf, mask_soilc, col.gridcell, grc.latdeg, prec60_col,
        peatf_lf, fsat, wf2, tsoi17,
        cnfire_const.borealat, non_boreal_peatfire_c, boreal_peatfire_c,
        nonborpeat_fire_precip_denom, borpeat_fire_soilmoist_denom,
        SECSPHR, SECSPDAY, TFRZ, RPI
    )

    # --- Find which pool is the CWD pool (tiny, host-resolved) ---
    i_cwd = 0
    for l in 1:ndecomp_pools
        if is_cwd[l]
            i_cwd = l
        end
    end

    # Resolve the global decomposition-thickness Ref to a concrete array (read on device).
    dzsoi_decomp_arr = dzsoi_decomp[]

    # --- Main column loop: fractional area affected by fire ---
    _launch!(_firea_fire_spread_kernel!, nfire, farea_burned, fuelc, fbac, fbac1,
             mask_soilc, col.gridcell, fire_li2014.forc_hdm, cropf_col,
             trotr1_col, trotr2_col, baf_crop, baf_peatf,
             totlitc, totvegc, rootc_col, fuelc_crop,
             decomp_cpools_vr, dzsoi_decomp_arr, wf, forc_rh,
             forc_t, fire_li2014.forc_lnfm, grc.lat, forc_wind,
             lgdp_col, lgdp1_col, lpop_col, fsr_col,
             fd_col, btran_col, wtlf, dtrotr_col,
             prec60_col, prec10_col, forc_rain, forc_snow, lfc,
             i_cwd, nlevdecomp, lfuel, ufuel, rh_low, rh_hgh, bt_min, bt_max,
             pot_hmn_ign_counts_alpha, cnfire_params.ignition_efficiency,
             cnfire_const.g0, cli_scale,
             defo_fire_precip_thresh_bet, defo_fire_precip_thresh_bdt,
             TFRZ, RPI, SECSPHR, SECSPDAY, dayspyr, dt,
             date_is_jan1_0, transient_landcover)

    return nothing
end

# ---------------------------------------------------------------------------
# cnfire_fluxes_li2014! — Fire C/N flux calculations (Li2014 version)
# ---------------------------------------------------------------------------

"""
    cnfire_fluxes_li2014!(mask_soilc, mask_soilp, bounds_c, bounds_p,
        cnfire_const, pftcon, patch, col, grc,
        dgvs, cnveg_state, cnveg_cs, cnveg_cf, cnveg_ns, cnveg_nf,
        soilbgc_cf, decomp_cascade_con,
        leaf_prof_patch, froot_prof_patch, croot_prof_patch, stem_prof_patch,
        totsomc_col, decomp_cpools_vr_col, decomp_npools_vr_col, somc_fire_col;
        kwargs...)

Fire effects routine for the Li et al. (2014) fire model.
Delegates to `cnfire_fluxes!` with Li2014-specific combustion completeness
factors: 0.5 for litter (changed from 0.4) and 0.25 for CWD (changed from 0.2)
per Li et al. (2014).

Ported from `CNFireFluxes` in `CNFireLi2014Mod.F90`.
"""
function cnfire_fluxes_li2014!(
    mask_soilc::BitVector,
    mask_soilp::BitVector,
    bounds_c::UnitRange{Int},
    bounds_p::UnitRange{Int},
    cnfire_const::CNFireConstData,
    pftcon::PftConFireBase,
    patch::PatchData,
    col::ColumnData,
    grc::GridcellData,
    dgvs::DgvsFireData,
    cnveg_state::CNVegStateData,
    cnveg_cs::CNVegCarbonStateData,
    cnveg_cf::CNVegCarbonFluxData,
    cnveg_ns::CNVegNitrogenStateData,
    cnveg_nf::CNVegNitrogenFluxData,
    soilbgc_cf::SoilBiogeochemCarbonFluxData,
    decomp_cascade_con::DecompCascadeConData,
    leaf_prof_patch::Matrix{<:Real},
    froot_prof_patch::Matrix{<:Real},
    croot_prof_patch::Matrix{<:Real},
    stem_prof_patch::Matrix{<:Real},
    totsomc_col::Vector{<:Real},
    decomp_cpools_vr_col::Array{<:Real,3},
    decomp_npools_vr_col::Array{<:Real,3},
    somc_fire_col::Vector{<:Real};
    dt::Real = 1800.0,
    dayspyr::Real = 365.0,
    nlevdecomp::Int = 1,
    ndecomp_pools::Int = 7,
    i_met_lit::Int = 1,
    i_litr_max::Int = 3,
    transient_landcover::Bool = false,
    use_cndv::Bool = false,
    use_matrixcn::Bool = false,
    spinup_factor_deadwood::Real = SPINUP_FACTOR_DEADWOOD_DEFAULT,
    secspday::Real = SECSPDAY,
    nc3crop::Int = 15,
    noveg::Int = 0,
    kmo::Int = 6,
    kda::Int = 15,
    mcsec::Int = 3600
)
    # Li2014 uses specific combustion completeness factors:
    # 0.5 for litter (changed from base 0.4) and 0.25 for CWD (changed from base 0.2)
    # per Li et al. (2014). Create a local copy with Li2014-specific values.
    li2014_const = CNFireConstData(
        borealat                     = cnfire_const.borealat,
        lfuel                        = cnfire_const.lfuel,
        ufuel                        = cnfire_const.ufuel,
        g0                           = cnfire_const.g0,
        rh_low                       = cnfire_const.rh_low,
        rh_hgh                       = cnfire_const.rh_hgh,
        bt_min                       = cnfire_const.bt_min,
        bt_max                       = cnfire_const.bt_max,
        cli_scale                    = cnfire_const.cli_scale,
        boreal_peatfire_c            = cnfire_const.boreal_peatfire_c,
        pot_hmn_ign_counts_alpha     = cnfire_const.pot_hmn_ign_counts_alpha,
        non_boreal_peatfire_c        = cnfire_const.non_boreal_peatfire_c,
        cropfire_a1                  = cnfire_const.cropfire_a1,
        occur_hi_gdp_tree            = cnfire_const.occur_hi_gdp_tree,
        cmb_cmplt_fact_litter        = 0.5,   # Li2014 value
        cmb_cmplt_fact_cwd           = 0.25,  # Li2014 value
        max_rh30_affecting_fuel      = cnfire_const.max_rh30_affecting_fuel,
        defo_fire_precip_thresh_bet  = cnfire_const.defo_fire_precip_thresh_bet,
        defo_fire_precip_thresh_bdt  = cnfire_const.defo_fire_precip_thresh_bdt,
        borpeat_fire_soilmoist_denom = cnfire_const.borpeat_fire_soilmoist_denom,
        nonborpeat_fire_precip_denom = cnfire_const.nonborpeat_fire_precip_denom,
    )

    return cnfire_fluxes!(
        mask_soilc, mask_soilp, bounds_c, bounds_p,
        li2014_const, pftcon, patch, col, grc,
        dgvs, cnveg_state, cnveg_cs, cnveg_cf, cnveg_ns, cnveg_nf,
        soilbgc_cf, decomp_cascade_con,
        leaf_prof_patch, froot_prof_patch, croot_prof_patch, stem_prof_patch,
        totsomc_col, decomp_cpools_vr_col, decomp_npools_vr_col, somc_fire_col;
        dt = dt,
        dayspyr = dayspyr,
        nlevdecomp = nlevdecomp,
        ndecomp_pools = ndecomp_pools,
        i_met_lit = i_met_lit,
        i_litr_max = i_litr_max,
        transient_landcover = transient_landcover,
        use_cndv = use_cndv,
        use_matrixcn = use_matrixcn,
        spinup_factor_deadwood = spinup_factor_deadwood,
        secspday = secspday,
        nc3crop = nc3crop,
        noveg = noveg,
        kmo = kmo,
        kda = kda,
        mcsec = mcsec,
    )
end
