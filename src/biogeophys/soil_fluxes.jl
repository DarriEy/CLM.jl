# ==========================================================================
# SoilFluxes — ported from SoilFluxesMod.F90
#
# Update surface fluxes based on the new ground temperature.
#
# GPU note: the whole driver runs on the device. Every per-(c)/per-(p) loop is a
# KernelAbstractions kernel; the patch->column errsoi average uses _scatter_add!
# (atomic on GPU, += on CPU/Dual). To stay under Metal's ~31 kernel-arg cap the
# many state arrays each loop touches are grouped into immutable device-view
# bundle structs (Adapt.@adapt_structure'd). The bundle FIELD NAMES are local
# aliases chosen so the kernel bodies read like the original scalar loops; on the
# host path the same bundles are built from CPU arrays, so CPU stays byte-identical.
# ==========================================================================

# --------------------------------------------------------------------------
# Kernel: temperature difference for flux corrections (per column, masked).
# Writes t_grnd0[c] and tinc[c]; both fully independent per column.
# --------------------------------------------------------------------------
@kernel function _soilflux_tinc_kernel!(tinc, @Const(mask), t_grnd0,
                                        @Const(snl), @Const(frac_sno_eff_col),
                                        @Const(frac_h2osfc_col), @Const(t_ssbef_col),
                                        @Const(t_h2osfc_bef_col), @Const(t_grnd_col),
                                        nlevsno::Int)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(tinc)
        if snl[c] < 0
            t_grnd0[c] = frac_sno_eff_col[c] *
                    t_ssbef_col[c, snl[c] + 1 + nlevsno] +
                (one(T) - frac_sno_eff_col[c] - frac_h2osfc_col[c]) *
                    t_ssbef_col[c, 1 + nlevsno] +
                frac_h2osfc_col[c] * t_h2osfc_bef_col[c]
        else
            t_grnd0[c] = (one(T) - frac_h2osfc_col[c]) *
                    t_ssbef_col[c, 1 + nlevsno] +
                frac_h2osfc_col[c] * t_h2osfc_bef_col[c]
        end
        tinc[c] = t_grnd_col[c] - t_grnd0[c]
    end
end

soilflux_tinc!(tinc, mask, t_grnd0, snl, frac_sno_eff_col, frac_h2osfc_col,
               t_ssbef_col, t_h2osfc_bef_col, t_grnd_col, nlevsno::Int) =
    _launch!(_soilflux_tinc_kernel!, tinc, mask, t_grnd0, snl, frac_sno_eff_col,
             frac_h2osfc_col, t_ssbef_col, t_h2osfc_bef_col, t_grnd_col, nlevsno)

# --------------------------------------------------------------------------
# Kernel: correct fluxes to present soil temperature (per patch, masked).
# Each patch updates only its own index; no cross-patch dependence.
# --------------------------------------------------------------------------
@kernel function _soilflux_correct_kernel!(eflx_sh_grnd_patch, @Const(mask),
                                           @Const(column), @Const(landunit),
                                           @Const(urbpoi), @Const(tinc),
                                           @Const(cgrnds_patch), @Const(cgrndl_patch),
                                           qflx_evap_soi_patch, qflx_ev_soil_patch,
                                           qflx_ev_h2osfc_patch, qflx_ev_snow_patch)
    p = @index(Global)
    @inbounds if mask[p]
        T = eltype(eflx_sh_grnd_patch)
        c = column[p]
        eflx_sh_grnd_patch[p] += tinc[c] * cgrnds_patch[p]
        qflx_evap_soi_patch[p] += tinc[c] * cgrndl_patch[p]

        l = landunit[p]
        if urbpoi[l]
            qflx_ev_soil_patch[p] = zero(T)
            qflx_ev_h2osfc_patch[p] = zero(T)
            qflx_ev_snow_patch[p] = qflx_evap_soi_patch[p]
        else
            qflx_ev_snow_patch[p] += tinc[c] * cgrndl_patch[p]
            qflx_ev_soil_patch[p] += tinc[c] * cgrndl_patch[p]
            qflx_ev_h2osfc_patch[p] += tinc[c] * cgrndl_patch[p]
        end
    end
end

soilflux_correct!(eflx_sh_grnd_patch, mask, column, landunit, urbpoi, tinc,
                  cgrnds_patch, cgrndl_patch, qflx_evap_soi_patch, qflx_ev_soil_patch,
                  qflx_ev_h2osfc_patch, qflx_ev_snow_patch) =
    _launch!(_soilflux_correct_kernel!, eflx_sh_grnd_patch, mask, column, landunit,
             urbpoi, tinc, cgrnds_patch, cgrndl_patch, qflx_evap_soi_patch,
             qflx_ev_soil_patch, qflx_ev_h2osfc_patch, qflx_ev_snow_patch)

# --------------------------------------------------------------------------
# Device-view bundles for the heavier per-patch loops. Each is an immutable,
# isbits-on-device struct (Adapt.@adapt_structure) holding only the arrays its
# loop body needs; field names mirror the original Julia variable paths so the
# kernel bodies read verbatim. Built per-launch from the live state structs, so
# all writes flow straight back into them (shared array refs).
# --------------------------------------------------------------------------

# Patch- and column-level water flux arrays the partition/limit/total loops touch.
Base.@kwdef struct _SFWaterDV{V}
    qflx_liqevap_from_top_layer_patch::V
    qflx_solidevap_from_top_layer_patch::V
    qflx_soliddew_to_top_layer_patch::V
    qflx_liqdew_to_top_layer_patch::V
    qflx_evap_soi_patch::V
    qflx_ev_snow_patch::V
    qflx_evap_tot_patch::V
    qflx_evap_veg_patch::V
    qflx_evap_can_patch::V
    qflx_tran_veg_patch::V
end
Adapt.@adapt_structure _SFWaterDV

_sf_water_dv(wfb) = _SFWaterDV(;
    qflx_liqevap_from_top_layer_patch  = wfb.wf.qflx_liqevap_from_top_layer_patch,
    qflx_solidevap_from_top_layer_patch = wfb.wf.qflx_solidevap_from_top_layer_patch,
    qflx_soliddew_to_top_layer_patch   = wfb.wf.qflx_soliddew_to_top_layer_patch,
    qflx_liqdew_to_top_layer_patch     = wfb.wf.qflx_liqdew_to_top_layer_patch,
    qflx_evap_soi_patch                = wfb.wf.qflx_evap_soi_patch,
    qflx_ev_snow_patch                 = wfb.qflx_ev_snow_patch,
    qflx_evap_tot_patch                = wfb.wf.qflx_evap_tot_patch,
    qflx_evap_veg_patch                = wfb.wf.qflx_evap_veg_patch,
    qflx_evap_can_patch                = wfb.wf.qflx_evap_can_patch,
    qflx_tran_veg_patch                = wfb.wf.qflx_tran_veg_patch)

# Energy flux arrays the total/ground-heat/lwrad loops touch.
Base.@kwdef struct _SFEnergyDV{V}
    eflx_sh_grnd_patch::V; htvp_col::V
    eflx_soil_grnd_patch::V; eflx_soil_grnd_r_patch::V; eflx_soil_grnd_u_patch::V
    dlrad_patch::V; eflx_lwrad_net_patch::V; eflx_lwrad_out_patch::V
    eflx_wasteheat_patch::V; eflx_heat_from_ac_patch::V; eflx_traffic_patch::V
    eflx_ventilation_patch::V
    eflx_sh_tot_patch::V; eflx_sh_veg_patch::V; eflx_sh_stem_patch::V
    eflx_lh_tot_patch::V; eflx_lh_tot_r_patch::V; eflx_sh_tot_r_patch::V
    eflx_lh_tot_u_patch::V; eflx_sh_tot_u_patch::V
    eflx_lh_vege_patch::V; eflx_lh_vegt_patch::V; eflx_lh_grnd_patch::V
    ulrad_patch::V; eflx_lwrad_net_r_patch::V; eflx_lwrad_out_r_patch::V
    eflx_lwrad_net_u_patch::V; eflx_lwrad_out_u_patch::V
end
Adapt.@adapt_structure _SFEnergyDV

_sf_energy_dv(ef) = _SFEnergyDV(;
    eflx_sh_grnd_patch = ef.eflx_sh_grnd_patch, htvp_col = ef.htvp_col,
    eflx_soil_grnd_patch = ef.eflx_soil_grnd_patch,
    eflx_soil_grnd_r_patch = ef.eflx_soil_grnd_r_patch,
    eflx_soil_grnd_u_patch = ef.eflx_soil_grnd_u_patch,
    dlrad_patch = ef.dlrad_patch, eflx_lwrad_net_patch = ef.eflx_lwrad_net_patch,
    eflx_lwrad_out_patch = ef.eflx_lwrad_out_patch,
    eflx_wasteheat_patch = ef.eflx_wasteheat_patch,
    eflx_heat_from_ac_patch = ef.eflx_heat_from_ac_patch,
    eflx_traffic_patch = ef.eflx_traffic_patch,
    eflx_ventilation_patch = ef.eflx_ventilation_patch,
    eflx_sh_tot_patch = ef.eflx_sh_tot_patch, eflx_sh_veg_patch = ef.eflx_sh_veg_patch,
    eflx_sh_stem_patch = ef.eflx_sh_stem_patch, eflx_lh_tot_patch = ef.eflx_lh_tot_patch,
    eflx_lh_tot_r_patch = ef.eflx_lh_tot_r_patch, eflx_sh_tot_r_patch = ef.eflx_sh_tot_r_patch,
    eflx_lh_tot_u_patch = ef.eflx_lh_tot_u_patch, eflx_sh_tot_u_patch = ef.eflx_sh_tot_u_patch,
    eflx_lh_vege_patch = ef.eflx_lh_vege_patch, eflx_lh_vegt_patch = ef.eflx_lh_vegt_patch,
    eflx_lh_grnd_patch = ef.eflx_lh_grnd_patch, ulrad_patch = ef.ulrad_patch,
    eflx_lwrad_net_r_patch = ef.eflx_lwrad_net_r_patch,
    eflx_lwrad_out_r_patch = ef.eflx_lwrad_out_r_patch,
    eflx_lwrad_net_u_patch = ef.eflx_lwrad_net_u_patch,
    eflx_lwrad_out_u_patch = ef.eflx_lwrad_out_u_patch)

# Column-level temperature / diagnostic / soil-state arrays the loops read.
Base.@kwdef struct _SFColDV{V,M,VI}
    frac_sno_eff_col::V; frac_h2osfc_col::V
    t_ssbef_col::M; t_h2osfc_bef_col::V; t_grnd_col::V; emg_col::V
    t_soisno_col::M; fact_col::M; t_skin_patch::V
    h2osoi_liq_col::M; h2osoi_ice_col::M
    sabg_soil_patch::V; sabg_snow_patch::V; sabg_patch::V
    frac_veg_nosno_patch::VI
end
Adapt.@adapt_structure _SFColDV

_sf_col_dv(temp, wdb, wsb, sa, cs) = _SFColDV(;
    frac_sno_eff_col = wdb.frac_sno_eff_col, frac_h2osfc_col = wdb.frac_h2osfc_col,
    t_ssbef_col = temp.t_ssbef_col, t_h2osfc_bef_col = temp.t_h2osfc_bef_col,
    t_grnd_col = temp.t_grnd_col, emg_col = temp.emg_col,
    t_soisno_col = temp.t_soisno_col, fact_col = temp.fact_col,
    t_skin_patch = temp.t_skin_patch,
    h2osoi_liq_col = wsb.ws.h2osoi_liq_col, h2osoi_ice_col = wsb.ws.h2osoi_ice_col,
    sabg_soil_patch = sa.sabg_soil_patch, sabg_snow_patch = sa.sabg_snow_patch,
    sabg_patch = sa.sabg_patch, frac_veg_nosno_patch = cs.frac_veg_nosno_patch)

# --------------------------------------------------------------------------
# Kernel: partition evaporation into liquid and solid (Loop 3, per patch).
# --------------------------------------------------------------------------
@kernel function _soilflux_partition_kernel!(@Const(_out), @Const(mask), @Const(column),
        @Const(landunit), @Const(snl), @Const(urbpoi), wfdv, cdv, nlevsno::Int, tfrz)
    p = @index(Global)
    @inbounds if mask[p]
        T = typeof(tfrz)
        c = column[p]
        l = landunit[p]
        j = snl[c] + 1            # Fortran-style top layer index
        j_jl = j + nlevsno        # Julia array index

        wfdv.qflx_liqevap_from_top_layer_patch[p] = zero(T)
        wfdv.qflx_solidevap_from_top_layer_patch[p] = zero(T)
        wfdv.qflx_soliddew_to_top_layer_patch[p] = zero(T)
        wfdv.qflx_liqdew_to_top_layer_patch[p] = zero(T)

        h2o_liq = cdv.h2osoi_liq_col[c, j_jl]
        h2o_ice = cdv.h2osoi_ice_col[c, j_jl]

        if !urbpoi[l]
            if wfdv.qflx_ev_snow_patch[p] >= zero(T)
                if (h2o_liq + h2o_ice) > zero(T)
                    wfdv.qflx_liqevap_from_top_layer_patch[p] = max(
                        wfdv.qflx_ev_snow_patch[p] * (h2o_liq / (h2o_liq + h2o_ice)), zero(T))
                else
                    wfdv.qflx_liqevap_from_top_layer_patch[p] = zero(T)
                end
                wfdv.qflx_solidevap_from_top_layer_patch[p] =
                    wfdv.qflx_ev_snow_patch[p] -
                    wfdv.qflx_liqevap_from_top_layer_patch[p]
            else
                if cdv.t_grnd_col[c] < tfrz
                    wfdv.qflx_soliddew_to_top_layer_patch[p] = abs(wfdv.qflx_ev_snow_patch[p])
                else
                    wfdv.qflx_liqdew_to_top_layer_patch[p] = abs(wfdv.qflx_ev_snow_patch[p])
                end
            end
        else  # Urban columns
            if wfdv.qflx_evap_soi_patch[p] >= zero(T)
                if (h2o_liq + h2o_ice) > zero(T)
                    wfdv.qflx_liqevap_from_top_layer_patch[p] = max(
                        wfdv.qflx_evap_soi_patch[p] * (h2o_liq / (h2o_liq + h2o_ice)), zero(T))
                else
                    wfdv.qflx_liqevap_from_top_layer_patch[p] = zero(T)
                end
                wfdv.qflx_solidevap_from_top_layer_patch[p] =
                    wfdv.qflx_evap_soi_patch[p] -
                    wfdv.qflx_liqevap_from_top_layer_patch[p]
            else
                if cdv.t_grnd_col[c] < tfrz
                    wfdv.qflx_soliddew_to_top_layer_patch[p] = abs(wfdv.qflx_evap_soi_patch[p])
                else
                    wfdv.qflx_liqdew_to_top_layer_patch[p] = abs(wfdv.qflx_evap_soi_patch[p])
                end
            end
        end
    end
end

# --------------------------------------------------------------------------
# Kernel: constrain evaporation from snow to <= available moisture (Loop 4).
# --------------------------------------------------------------------------
@kernel function _soilflux_evaplimit_kernel!(@Const(_out), @Const(mask), @Const(column),
        @Const(landunit), @Const(snl), @Const(itype), @Const(urbpoi),
        wfdv, edv, cdv, nlevsno::Int, dtime, icol_road_perv::Int)
    p = @index(Global)
    @inbounds if mask[p]
        T = typeof(dtime)
        c = column[p]
        j = snl[c] + 1
        j_jl = j + nlevsno

        h2o_liq = cdv.h2osoi_liq_col[c, j_jl]
        h2o_ice = cdv.h2osoi_ice_col[c, j_jl]

        # Snow layers; assumes for j < 1 that frac_sno_eff > 0
        if j < 1
            evaporation_limit = (h2o_ice + h2o_liq) / (cdv.frac_sno_eff_col[c] * dtime)
            if wfdv.qflx_ev_snow_patch[p] > evaporation_limit
                evaporation_demand = wfdv.qflx_ev_snow_patch[p]
                wfdv.qflx_ev_snow_patch[p] = evaporation_limit
                wfdv.qflx_evap_soi_patch[p] -=
                    cdv.frac_sno_eff_col[c] * (evaporation_demand - evaporation_limit)
                wfdv.qflx_liqevap_from_top_layer_patch[p] =
                    max(h2o_liq / (cdv.frac_sno_eff_col[c] * dtime), zero(T))
                wfdv.qflx_solidevap_from_top_layer_patch[p] =
                    max(h2o_ice / (cdv.frac_sno_eff_col[c] * dtime), zero(T))
                edv.eflx_sh_grnd_patch[p] +=
                    cdv.frac_sno_eff_col[c] * (evaporation_demand - evaporation_limit) * edv.htvp_col[c]
            end
        end

        # Top soil layer for urban columns (excluding pervious road)
        if urbpoi[landunit[p]] && (itype[c] != icol_road_perv) && (j == 1)
            evaporation_limit = (h2o_ice + h2o_liq) / dtime
            if wfdv.qflx_evap_soi_patch[p] > evaporation_limit
                evaporation_demand = wfdv.qflx_evap_soi_patch[p]
                wfdv.qflx_evap_soi_patch[p] = evaporation_limit
                wfdv.qflx_ev_snow_patch[p] = wfdv.qflx_evap_soi_patch[p]
                wfdv.qflx_liqevap_from_top_layer_patch[p] = max(h2o_liq / dtime, zero(T))
                wfdv.qflx_solidevap_from_top_layer_patch[p] = max(h2o_ice / dtime, zero(T))
                edv.eflx_sh_grnd_patch[p] +=
                    (evaporation_demand - evaporation_limit) * edv.htvp_col[c]
            end
        end

        # Limit only solid evaporation (sublimation) from top soil layer
        if j == 1 && cdv.frac_h2osfc_col[c] < one(T)
            evaporation_limit = h2o_ice / (dtime * (one(T) - cdv.frac_h2osfc_col[c]))
            if wfdv.qflx_solidevap_from_top_layer_patch[p] >= evaporation_limit
                evaporation_demand = wfdv.qflx_solidevap_from_top_layer_patch[p]
                wfdv.qflx_solidevap_from_top_layer_patch[p] = evaporation_limit
                wfdv.qflx_liqevap_from_top_layer_patch[p] +=
                    (evaporation_demand - evaporation_limit)
            end
        end
    end
end

# --------------------------------------------------------------------------
# Kernel: ground heat flux, total fluxes, history variables (Loop 5).
# eflx_lwrad_del[p] is an output (urban only); reused by Loop 7.
# --------------------------------------------------------------------------
@kernel function _soilflux_totals_kernel!(eflx_lwrad_del, @Const(mask),
        @Const(column), @Const(landunit), @Const(snl), @Const(lun_itype),
        @Const(urbpoi), @Const(forc_lwrad_col), @Const(tinc), @Const(t_grnd0),
        edv, wfdv, cdv, nlevsno::Int, sb, hvap,
        istsoil::Int, istcrop::Int)
    p = @index(Global)
    @inbounds if mask[p]
        T = eltype(eflx_lwrad_del)
        c = column[p]
        l = landunit[p]

        if !urbpoi[l]
            lw_grnd = cdv.frac_sno_eff_col[c] *
                    cdv.t_ssbef_col[c, snl[c] + 1 + nlevsno]^4 +
                (one(T) - cdv.frac_sno_eff_col[c] - cdv.frac_h2osfc_col[c]) *
                    cdv.t_ssbef_col[c, 1 + nlevsno]^4 +
                cdv.frac_h2osfc_col[c] * cdv.t_h2osfc_bef_col[c]^4

            edv.eflx_soil_grnd_patch[p] =
                ((one(T) - cdv.frac_sno_eff_col[c]) * cdv.sabg_soil_patch[p] +
                 cdv.frac_sno_eff_col[c] * cdv.sabg_snow_patch[p]) +
                edv.dlrad_patch[p] +
                (one(T) - cdv.frac_veg_nosno_patch[p]) * cdv.emg_col[c] * forc_lwrad_col[c] -
                cdv.emg_col[c] * sb * lw_grnd -
                cdv.emg_col[c] * sb * t_grnd0[c]^3 * (T(4) * tinc[c]) -
                (edv.eflx_sh_grnd_patch[p] + wfdv.qflx_evap_soi_patch[p] * edv.htvp_col[c])

            if lun_itype[l] == istsoil || lun_itype[l] == istcrop
                edv.eflx_soil_grnd_r_patch[p] = edv.eflx_soil_grnd_patch[p]
            end
        else
            eflx_lwrad_del[p] = T(4) * cdv.emg_col[c] * sb * t_grnd0[c]^3 * tinc[c]

            edv.eflx_soil_grnd_patch[p] = cdv.sabg_patch[p] + edv.dlrad_patch[p] -
                edv.eflx_lwrad_net_patch[p] - eflx_lwrad_del[p] -
                (edv.eflx_sh_grnd_patch[p] +
                 wfdv.qflx_evap_soi_patch[p] * edv.htvp_col[c] +
                 wfdv.qflx_tran_veg_patch[p] * hvap) +
                edv.eflx_wasteheat_patch[p] + edv.eflx_heat_from_ac_patch[p] +
                edv.eflx_traffic_patch[p] + edv.eflx_ventilation_patch[p]
            edv.eflx_soil_grnd_u_patch[p] = edv.eflx_soil_grnd_patch[p]
        end

        # Total fluxes (vegetation + ground)
        edv.eflx_sh_tot_patch[p] = edv.eflx_sh_veg_patch[p] + edv.eflx_sh_grnd_patch[p]
        if !urbpoi[l]
            edv.eflx_sh_tot_patch[p] += edv.eflx_sh_stem_patch[p]
        end
        wfdv.qflx_evap_tot_patch[p] =
            wfdv.qflx_evap_veg_patch[p] + wfdv.qflx_evap_soi_patch[p]

        edv.eflx_lh_tot_patch[p] =
            hvap * wfdv.qflx_evap_veg_patch[p] +
            edv.htvp_col[c] * wfdv.qflx_evap_soi_patch[p]
        if lun_itype[l] == istsoil || lun_itype[l] == istcrop
            edv.eflx_lh_tot_r_patch[p] = edv.eflx_lh_tot_patch[p]
            edv.eflx_sh_tot_r_patch[p] = edv.eflx_sh_tot_patch[p]
        elseif urbpoi[l]
            edv.eflx_lh_tot_u_patch[p] = edv.eflx_lh_tot_patch[p]
            edv.eflx_sh_tot_u_patch[p] = edv.eflx_sh_tot_patch[p]
        end

        # Variables needed by history tape
        wfdv.qflx_evap_can_patch[p] =
            wfdv.qflx_evap_veg_patch[p] - wfdv.qflx_tran_veg_patch[p]
        edv.eflx_lh_vege_patch[p] =
            (wfdv.qflx_evap_veg_patch[p] - wfdv.qflx_tran_veg_patch[p]) * hvap
        edv.eflx_lh_vegt_patch[p] = wfdv.qflx_tran_veg_patch[p] * hvap
        edv.eflx_lh_grnd_patch[p] = wfdv.qflx_evap_soi_patch[p] * edv.htvp_col[c]
    end
end

# --------------------------------------------------------------------------
# Kernel: base soil energy balance error (per patch, masked).
# Covers the per-patch terms only; the layer-summation contributions
# (loop-carried over j) are intentionally left as scalar loops.
# --------------------------------------------------------------------------
@kernel function _soilflux_errsoi_base_kernel!(errsoi_patch, @Const(mask),
                                               @Const(column), @Const(itype),
                                               @Const(eflx_soil_grnd_patch),
                                               @Const(xmf_col), @Const(xmf_h2osfc_col),
                                               @Const(frac_h2osfc_col),
                                               @Const(t_h2osfc_col),
                                               @Const(t_h2osfc_bef_col),
                                               @Const(c_h2osfc_col),
                                               @Const(eflx_h2osfc_to_snow_col),
                                               @Const(eflx_building_heat_errsoi_col),
                                               dtime,
                                               icol_sunwall::Int, icol_shadewall::Int,
                                               icol_roof::Int)
    p = @index(Global)
    @inbounds if mask[p]
        c = column[p]
        errsoi_patch[p] = eflx_soil_grnd_patch[p] -
            xmf_col[c] - xmf_h2osfc_col[c] -
            frac_h2osfc_col[c] *
            (t_h2osfc_col[c] - t_h2osfc_bef_col[c]) *
            (c_h2osfc_col[c] / dtime)
        errsoi_patch[p] += eflx_h2osfc_to_snow_col[c]

        if itype[c] == icol_sunwall || itype[c] == icol_shadewall || itype[c] == icol_roof
            errsoi_patch[p] += eflx_building_heat_errsoi_col[c]
        end
    end
end

soilflux_errsoi_base!(errsoi_patch, mask, column, itype, eflx_soil_grnd_patch,
                      xmf_col, xmf_h2osfc_col, frac_h2osfc_col, t_h2osfc_col,
                      t_h2osfc_bef_col, c_h2osfc_col, eflx_h2osfc_to_snow_col,
                      eflx_building_heat_errsoi_col, dtime,
                      icol_sunwall::Int, icol_shadewall::Int, icol_roof::Int) =
    _launch!(_soilflux_errsoi_base_kernel!, errsoi_patch, mask, column, itype,
             eflx_soil_grnd_patch, xmf_col, xmf_h2osfc_col, frac_h2osfc_col,
             t_h2osfc_col, t_h2osfc_bef_col, c_h2osfc_col, eflx_h2osfc_to_snow_col,
             eflx_building_heat_errsoi_col, dtime, icol_sunwall, icol_shadewall, icol_roof)

# --------------------------------------------------------------------------
# Kernel: errsoi layer-summation contributions (per patch; internal j loop so
# the per-patch accumulation has no cross-thread dependence). Handles both the
# non-urban (mask_nolakep, j up to nlevgrnd) and urban (mask_urbanp, j up to
# nlevurb) passes; `is_urban` selects which column types contribute.
# --------------------------------------------------------------------------
@kernel function _soilflux_errsoi_layers_kernel!(errsoi_patch, @Const(mask),
        @Const(column), @Const(itype), @Const(snl), @Const(frac_sno_eff_col),
        @Const(t_soisno_col), @Const(t_ssbef_col), @Const(fact_col),
        nlevsno::Int, jmax::Int, is_urban::Bool,
        icol_sunwall::Int, icol_shadewall::Int, icol_roof::Int)
    p = @index(Global)
    @inbounds if mask[p]
        c = column[p]
        is_wall_roof = (itype[c] == icol_sunwall || itype[c] == icol_shadewall ||
                        itype[c] == icol_roof)
        # non-urban pass: skip wall/roof columns; urban pass: only wall/roof columns
        if is_urban == is_wall_roof
            for j_f in (-nlevsno + 1):jmax
                j_jl = j_f + nlevsno
                if j_f >= snl[c] + 1 && j_f < 1
                    errsoi_patch[p] -= frac_sno_eff_col[c] *
                        (t_soisno_col[c, j_jl] - t_ssbef_col[c, j_jl]) / fact_col[c, j_jl]
                end
                if j_f >= 1
                    errsoi_patch[p] -=
                        (t_soisno_col[c, j_jl] - t_ssbef_col[c, j_jl]) / fact_col[c, j_jl]
                end
            end
        end
    end
end

# --------------------------------------------------------------------------
# Kernel: outgoing longwave radiation from vegetation + ground (Loop 7).
# --------------------------------------------------------------------------
@kernel function _soilflux_lwrad_kernel!(@Const(eflx_lwrad_del), @Const(mask),
        @Const(column), @Const(landunit), @Const(snl), @Const(lun_itype),
        @Const(urbpoi), @Const(forc_lwrad_col), @Const(tinc), @Const(t_grnd0),
        edv, cdv, nlevsno::Int, sb, istsoil::Int, istcrop::Int)
    p = @index(Global)
    @inbounds if mask[p]
        T = eltype(eflx_lwrad_del)
        c = column[p]
        l = landunit[p]

        if !urbpoi[l]
            lw_grnd = cdv.frac_sno_eff_col[c] *
                    cdv.t_ssbef_col[c, snl[c] + 1 + nlevsno]^4 +
                (one(T) - cdv.frac_sno_eff_col[c] - cdv.frac_h2osfc_col[c]) *
                    cdv.t_ssbef_col[c, 1 + nlevsno]^4 +
                cdv.frac_h2osfc_col[c] * cdv.t_h2osfc_bef_col[c]^4

            edv.eflx_lwrad_out_patch[p] = edv.ulrad_patch[p] +
                (one(T) - cdv.frac_veg_nosno_patch[p]) * (one(T) - cdv.emg_col[c]) * forc_lwrad_col[c] +
                (one(T) - cdv.frac_veg_nosno_patch[p]) * cdv.emg_col[c] * sb * lw_grnd +
                T(4) * cdv.emg_col[c] * sb * t_grnd0[c]^3 * tinc[c]

            if cdv.frac_veg_nosno_patch[p] == zero(T)
                cdv.t_skin_patch[p] = sqrt(sqrt(lw_grnd))
            end

            edv.eflx_lwrad_net_patch[p] = edv.eflx_lwrad_out_patch[p] - forc_lwrad_col[c]
            if lun_itype[l] == istsoil || lun_itype[l] == istcrop
                edv.eflx_lwrad_net_r_patch[p] = edv.eflx_lwrad_out_patch[p] - forc_lwrad_col[c]
                edv.eflx_lwrad_out_r_patch[p] = edv.eflx_lwrad_out_patch[p]
            end
        else
            edv.eflx_lwrad_out_patch[p] += eflx_lwrad_del[p]
            edv.eflx_lwrad_net_patch[p] += eflx_lwrad_del[p]
            edv.eflx_lwrad_net_u_patch[p] += eflx_lwrad_del[p]
            edv.eflx_lwrad_out_u_patch[p] = edv.eflx_lwrad_out_patch[p]
        end
    end
end

# --------------------------------------------------------------------------
# Kernels: patch->column errsoi average. Zero the accumulator (masked column),
# then scatter-add weighted patch contributions (atomic on GPU).
# --------------------------------------------------------------------------
@kernel function _soilflux_errsoi_zero_kernel!(errsoi_col, @Const(mask))
    c = @index(Global)
    @inbounds if mask[c]
        errsoi_col[c] = zero(eltype(errsoi_col))
    end
end

# Mirrors the original `for c (mask_c): for p in patchi[c]:patchf[c]` exactly: a
# patch p contributes to column c = column[p] iff that column is masked and p lies
# in the column's contiguous patch range (membership is what defines column[p]).
# No mask_nolakep filter — the scalar version summed every patch in the range.
@kernel function _soilflux_errsoi_p2c_kernel!(errsoi_col, @Const(mask_c),
        @Const(column), @Const(patchi), @Const(patchf),
        @Const(errsoi_patch), @Const(wtcol))
    p = @index(Global)
    @inbounds begin
        c = column[p]
        if c >= 1 && mask_c[c] && patchi[c] <= p && p <= patchf[c]
            _scatter_add!(errsoi_col, c, errsoi_patch[p] * wtcol[p])
        end
    end
end

# --------------------------------------------------------------------------
# Kernel: urban patch skin temperature from column top-layer soil temp.
# Per patch, masked; independent.
# --------------------------------------------------------------------------
@kernel function _soilflux_urban_tskin_kernel!(t_skin_patch, @Const(mask),
                                               @Const(column), @Const(snl),
                                               @Const(t_soisno_col), nlevsno::Int)
    p = @index(Global)
    @inbounds if mask[p]
        c = column[p]
        t_skin_patch[p] = t_soisno_col[c, snl[c] + 1 + nlevsno]
    end
end

soilflux_urban_tskin!(t_skin_patch, mask, column, snl, t_soisno_col, nlevsno::Int) =
    _launch!(_soilflux_urban_tskin_kernel!, t_skin_patch, mask, column, snl,
             t_soisno_col, nlevsno)

"""
    soil_fluxes!(...)

Update surface fluxes based on the new ground temperature.

Corrects sensible heat, latent heat, and soil heat fluxes for the change
in ground temperature between the previous and current time steps.
Partitions evaporation into liquid and solid components, constrains snow
evaporation to available moisture, computes total fluxes, outgoing
longwave radiation, and the soil energy balance error.

Ported from: SoilFluxesMod.F90 :: SoilFluxes
"""
function soil_fluxes!(
        # Data structures
        energyflux       ::EnergyFluxData,
        temperature      ::TemperatureData,
        canopystate      ::CanopyStateData,
        waterstatebulk   ::WaterStateBulkData,
        waterdiagbulk    ::WaterDiagnosticBulkData,
        waterfluxbulk    ::WaterFluxBulkData,
        solarabs         ::SolarAbsorbedData,
        patch_data       ::PatchData,
        col_data         ::ColumnData,
        lun_data         ::LandunitData,
        # Masks
        mask_nolakec     ::AbstractVector{Bool},
        mask_nolakep     ::AbstractVector{Bool},
        mask_urbanp      ::AbstractVector{Bool},
        # Bounds
        bounds_col       ::UnitRange{Int},
        bounds_patch     ::UnitRange{Int},
        # Atmospheric forcing (column-level)
        forc_lwrad_col   ::AbstractVector{<:Real},
        # Time step
        dtime            ::Real)

    nlevsno  = varpar.nlevsno
    nlevgrnd = varpar.nlevgrnd
    nlevurb  = varpar.nlevurb

    # Local work arrays (live on the same backend as the state arrays).
    endc = last(bounds_col)
    endp = last(bounds_patch)
    FT = eltype(temperature.t_grnd_col)
    proto_c = temperature.t_grnd_col           # backend/precision prototype (column)
    proto_p = energyflux.eflx_sh_grnd_patch    # backend/precision prototype (patch)
    tinc           = fill!(similar(proto_c, FT, endc), zero(FT))
    t_grnd0        = fill!(similar(proto_c, FT, endc), zero(FT))
    eflx_lwrad_del = fill!(similar(proto_p, FT, endp), zero(FT))

    # Scalar constants at working precision (no Float64 reaches a Float32 backend).
    dt    = convert(FT, dtime)
    tfrz  = convert(FT, TFRZ)
    sb    = convert(FT, SB)
    hvap  = convert(FT, HVAP)

    # Device-view bundles (shared array refs → writes flow back into the structs).
    wfdv = _sf_water_dv(waterfluxbulk)
    edv  = _sf_energy_dv(energyflux)
    cdv  = _sf_col_dv(temperature, waterdiagbulk, waterstatebulk, solarabs, canopystate)

    # =========================================================================
    # Loop 1: Calculate temperature difference for flux corrections (column)
    # =========================================================================

    soilflux_tinc!(tinc, mask_nolakec, t_grnd0, col_data.snl,
                   waterdiagbulk.frac_sno_eff_col, waterdiagbulk.frac_h2osfc_col,
                   temperature.t_ssbef_col, temperature.t_h2osfc_bef_col,
                   temperature.t_grnd_col, nlevsno)

    # =========================================================================
    # Loop 2: Correct fluxes to present soil temperature (patch)
    # =========================================================================

    soilflux_correct!(energyflux.eflx_sh_grnd_patch, mask_nolakep, patch_data.column,
                      patch_data.landunit, lun_data.urbpoi, tinc,
                      energyflux.cgrnds_patch, energyflux.cgrndl_patch,
                      waterfluxbulk.wf.qflx_evap_soi_patch, waterfluxbulk.qflx_ev_soil_patch,
                      waterfluxbulk.qflx_ev_h2osfc_patch, waterfluxbulk.qflx_ev_snow_patch)

    # =========================================================================
    # Loop 3: Partition evaporation into liquid and solid (patch)
    # =========================================================================

    _launch!(_soilflux_partition_kernel!, energyflux.eflx_sh_grnd_patch, mask_nolakep,
             patch_data.column, patch_data.landunit, col_data.snl, lun_data.urbpoi,
             wfdv, cdv, nlevsno, tfrz)

    # =========================================================================
    # Loop 4: Constrain evaporation from snow to be <= available moisture
    # =========================================================================

    _launch!(_soilflux_evaplimit_kernel!, energyflux.eflx_sh_grnd_patch, mask_nolakep,
             patch_data.column, patch_data.landunit, col_data.snl, col_data.itype,
             lun_data.urbpoi, wfdv, edv, cdv, nlevsno, dt, ICOL_ROAD_PERV)

    # =========================================================================
    # Loop 5: Ground heat flux, total fluxes, history variables (patch)
    # =========================================================================

    _launch!(_soilflux_totals_kernel!, eflx_lwrad_del, mask_nolakep,
             patch_data.column, patch_data.landunit, col_data.snl, lun_data.itype,
             lun_data.urbpoi, forc_lwrad_col, tinc, t_grnd0,
             edv, wfdv, cdv, nlevsno, sb, hvap, ISTSOIL, ISTCROP)

    # =========================================================================
    # Loop 6: Soil energy balance check (errsoi_patch)
    # =========================================================================

    soilflux_errsoi_base!(energyflux.errsoi_patch, mask_nolakep, patch_data.column,
                          col_data.itype, energyflux.eflx_soil_grnd_patch,
                          temperature.xmf_col, temperature.xmf_h2osfc_col,
                          waterdiagbulk.frac_h2osfc_col, temperature.t_h2osfc_col,
                          temperature.t_h2osfc_bef_col, temperature.c_h2osfc_col,
                          energyflux.eflx_h2osfc_to_snow_col,
                          energyflux.eflx_building_heat_errsoi_col, dt,
                          ICOL_SUNWALL, ICOL_SHADEWALL, ICOL_ROOF)

    # errsoi layer-summation: non-urban (j up to nlevgrnd) then urban (j up to nlevurb).
    _launch!(_soilflux_errsoi_layers_kernel!, energyflux.errsoi_patch, mask_nolakep,
             patch_data.column, col_data.itype, col_data.snl,
             waterdiagbulk.frac_sno_eff_col, temperature.t_soisno_col,
             temperature.t_ssbef_col, temperature.fact_col, nlevsno, nlevgrnd, false,
             ICOL_SUNWALL, ICOL_SHADEWALL, ICOL_ROOF)
    _launch!(_soilflux_errsoi_layers_kernel!, energyflux.errsoi_patch, mask_urbanp,
             patch_data.column, col_data.itype, col_data.snl,
             waterdiagbulk.frac_sno_eff_col, temperature.t_soisno_col,
             temperature.t_ssbef_col, temperature.fact_col, nlevsno, nlevurb, true,
             ICOL_SUNWALL, ICOL_SHADEWALL, ICOL_ROOF)

    # =========================================================================
    # Loop 7: Outgoing longwave radiation from vegetation + ground
    # =========================================================================

    _launch!(_soilflux_lwrad_kernel!, eflx_lwrad_del, mask_nolakep,
             patch_data.column, patch_data.landunit, col_data.snl, lun_data.itype,
             lun_data.urbpoi, forc_lwrad_col, tinc, t_grnd0,
             edv, cdv, nlevsno, sb, ISTSOIL, ISTCROP)

    # =========================================================================
    # Patch-to-column averaging for errsoi (replaces Fortran p2c call)
    # =========================================================================

    _launch!(_soilflux_errsoi_zero_kernel!, energyflux.errsoi_col, mask_nolakec)
    _launch!(_soilflux_errsoi_p2c_kernel!, energyflux.errsoi_col, mask_nolakec,
             patch_data.column, col_data.patchi, col_data.patchf,
             energyflux.errsoi_patch, patch_data.wtcol; ndrange = length(mask_nolakep))

    # =========================================================================
    # Assign column-level t_soisno(snl+1) to t_skin for urban patches
    # =========================================================================

    soilflux_urban_tskin!(temperature.t_skin_patch, mask_urbanp, patch_data.column,
                          col_data.snl, temperature.t_soisno_col, nlevsno)

    return nothing
end
