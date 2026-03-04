# ==========================================================================
# Ported from: src/main/lnd2atmMod.F90
# Land-to-atmosphere aggregation using subgrid averaging
#
# Public functions:
#   lnd2atm!         — Full land-to-atmosphere mapping
#   lnd2atm_minimal! — Albedo and radiative temperature only
# ==========================================================================

"""
    lnd2atm!(bounds, inst)

Map land model output fields from patch/column/landunit to gridcell level
for coupling to the atmosphere. Uses subgrid averaging (p2g, c2g).

Ported from `lnd2atm` in `lnd2atmMod.F90`.
"""
function lnd2atm!(bounds::BoundsType,
                  inst::CLMInstances)
    l2a = inst.lnd2atm
    ef  = inst.energyflux
    fv  = inst.frictionvel
    temp = inst.temperature
    wfb = inst.water.waterfluxbulk_inst
    wsb = inst.water.waterstatebulk_inst
    wdb = inst.water.waterdiagnosticbulk_inst
    sa  = inst.solarabs
    alb = inst.surfalb
    cs  = inst.canopystate
    col = inst.column
    lun = inst.landunit
    pch = inst.patch

    # --- Radiative temperature ---
    # t_rad = (eflx_lwrad_out / SB)^0.25
    for g in bounds.begg:bounds.endg
        lwout = 0.0
        count = 0
        for p in bounds.begp:bounds.endp
            if pch.active[p] && pch.gridcell[p] == g
                lwout += ef.eflx_lwrad_out_patch[p] * pch.wtgcell[p]
                count += 1
            end
        end
        if lwout > 0.0
            l2a.t_rad_grc[g] = (lwout / SB)^0.25
        else
            l2a.t_rad_grc[g] = temp.t_grnd_col[1]
        end
    end

    # --- 2m reference temperature (patch→gridcell) ---
    p2g_1d!(l2a.t_ref2m_grc, temp.t_ref2m_patch,
            bounds, "unity", "unity", "unity", pch, col, lun)

    # --- 10m wind speed (patch→gridcell) ---
    p2g_1d!(l2a.u_ref10m_grc, fv.u10_clm_patch,
            bounds, "unity", "unity", "unity", pch, col, lun)

    # --- Surface stress (patch→gridcell) ---
    p2g_1d!(l2a.taux_grc, ef.taux_patch,
            bounds, "unity", "unity", "unity", pch, col, lun)
    p2g_1d!(l2a.tauy_grc, ef.tauy_patch,
            bounds, "unity", "unity", "unity", pch, col, lun)

    # --- Roughness length (patch→gridcell) ---
    p2g_1d!(l2a.z0m_grc, fv.z0m_actual_patch,
            bounds, "unity", "unity", "unity", pch, col, lun)

    # --- Friction velocity and aerodynamic resistance ---
    p2g_1d!(l2a.fv_grc, fv.fv_patch,
            bounds, "unity", "unity", "unity", pch, col, lun)
    p2g_1d!(l2a.ram1_grc, fv.ram1_patch,
            bounds, "unity", "unity", "unity", pch, col, lun)

    # --- Energy fluxes ---
    # Sensible heat (patch→gridcell)
    p2g_1d!(l2a.eflx_sh_tot_grc, ef.eflx_sh_tot_patch,
            bounds, "unity", "unity", "unity", pch, col, lun)

    # Add precipitation conversion sensible heat (column→gridcell)
    if length(l2a.eflx_sh_precip_conversion_grc) > 0
        c2g_1d!(l2a.eflx_sh_precip_conversion_grc, ef.eflx_sh_precip_conversion_col,
                bounds, "unity", "unity", col, lun)
        for g in bounds.begg:bounds.endg
            l2a.eflx_sh_tot_grc[g] += l2a.eflx_sh_precip_conversion_grc[g]
        end
    end

    # Latent heat (patch→gridcell)
    p2g_1d!(l2a.eflx_lh_tot_grc, ef.eflx_lh_tot_patch,
            bounds, "unity", "unity", "unity", pch, col, lun)

    # Longwave outgoing (patch→gridcell)
    p2g_1d!(l2a.eflx_lwrad_out_grc, ef.eflx_lwrad_out_patch,
            bounds, "unity", "unity", "unity", pch, col, lun)

    # Absorbed solar (patch→gridcell)
    p2g_1d!(l2a.fsa_grc, sa.fsa_patch,
            bounds, "unity", "unity", "unity", pch, col, lun)

    # --- Albedos (patch→gridcell) ---
    for b in 1:NUMRAD
        # Direct beam
        garr_d = @view l2a.albd_grc[:, b]
        parr_d = @view alb.albd_patch[:, b]
        for g in bounds.begg:bounds.endg
            wsum = 0.0; val = 0.0
            for p in bounds.begp:bounds.endp
                if pch.active[p] && pch.gridcell[p] == g
                    w = pch.wtgcell[p]
                    val += parr_d[p] * w
                    wsum += w
                end
            end
            garr_d[g] = wsum > 0.0 ? val / wsum : 0.0
        end
        # Diffuse beam
        garr_i = @view l2a.albi_grc[:, b]
        parr_i = @view alb.albi_patch[:, b]
        for g in bounds.begg:bounds.endg
            wsum = 0.0; val = 0.0
            for p in bounds.begp:bounds.endp
                if pch.active[p] && pch.gridcell[p] == g
                    w = pch.wtgcell[p]
                    val += parr_i[p] * w
                    wsum += w
                end
            end
            garr_i[g] = wsum > 0.0 ? val / wsum : 0.0
        end
    end

    # --- Net carbon exchange ---
    # Simple: set to zero if not using CN
    for g in bounds.begg:bounds.endg
        l2a.net_carbon_exchange_grc[g] = 0.0
    end

    return nothing
end

"""
    lnd2atm_minimal!(bounds, inst)

Minimal land-to-atmosphere mapping: albedo and radiative temperature only.
Used for the first timestep before full physics runs.

Ported from `lnd2atm_minimal` in `lnd2atmMod.F90`.
"""
function lnd2atm_minimal!(bounds::BoundsType,
                           inst::CLMInstances)
    l2a = inst.lnd2atm
    ef  = inst.energyflux
    temp = inst.temperature
    alb = inst.surfalb
    col = inst.column
    lun = inst.landunit
    pch = inst.patch

    # Radiative temperature
    for g in bounds.begg:bounds.endg
        lwout = 0.0
        for p in bounds.begp:bounds.endp
            if pch.active[p] && pch.gridcell[p] == g
                lwout += ef.eflx_lwrad_out_patch[p] * pch.wtgcell[p]
            end
        end
        l2a.t_rad_grc[g] = lwout > 0.0 ? (lwout / SB)^0.25 : temp.t_grnd_col[1]
    end

    # Albedos
    for b in 1:NUMRAD
        garr_d = @view l2a.albd_grc[:, b]
        parr_d = @view alb.albd_patch[:, b]
        for g in bounds.begg:bounds.endg
            wsum = 0.0; val = 0.0
            for p in bounds.begp:bounds.endp
                if pch.active[p] && pch.gridcell[p] == g
                    w = pch.wtgcell[p]
                    val += parr_d[p] * w
                    wsum += w
                end
            end
            garr_d[g] = wsum > 0.0 ? val / wsum : 0.0
        end
    end

    return nothing
end
