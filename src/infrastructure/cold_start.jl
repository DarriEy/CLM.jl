# ==========================================================================
# Cold-start state initialization
# Ported from scattered *InitTimeConst subroutines in CLM Fortran
# ==========================================================================

"""
    cold_start_initialize!(inst, bounds, filt, surf)

Initialize all state variables for a cold start. Sets soil/snow temperatures,
moisture, soil properties from surface data, and satellite phenology.
"""
function cold_start_initialize!(inst::CLMInstances, bounds::BoundsType,
                                 filt::ClumpFilter, surf::SurfaceInputData;
                                 use_aquifer_layer::Bool = true)
    init_soil_properties!(inst, bounds, surf)
    init_temperatures!(inst, bounds)
    init_soil_moisture!(inst, bounds; use_aquifer_layer=use_aquifer_layer)
    init_snow_state!(inst, bounds)
    init_water_diagnostic_state!(inst, bounds)
    init_soil_state_defaults!(inst, bounds)
    init_eff_porosity!(inst, bounds)
    init_satellite_phenology!(inst, bounds, surf)
    init_misc_state!(inst, bounds)
    init_soil_hydrology_cold!(inst, bounds; use_aquifer_layer=use_aquifer_layer)
    init_lake_state!(inst, bounds)
    init_root_fractions!(inst, bounds)
    nothing
end

"""
    init_soil_properties!(inst, bounds, surf)

Compute soil physical properties (watsat, bsw, hksat, etc.) from sand/clay/organic.
Corresponds to SoilStateInitTimeConstMod.
"""
function init_soil_properties!(inst::CLMInstances, bounds::BoundsType,
                                surf::SurfaceInputData)
    col = inst.column
    ss = inst.soilstate
    nlevsoi = varpar.nlevsoi

    for c in bounds.begc:bounds.endc
        g = col.gridcell[c]
        gi = g - bounds.begg + 1

        for j in 1:nlevsoi
            # Surface data may only have 10 levels; use last available for deeper layers
            j_surf = min(j, size(surf.pct_sand, 2))
            sand = surf.pct_sand[gi, j_surf]
            clay = surf.pct_clay[gi, j_surf]
            om_frac = j_surf <= size(surf.organic, 2) ?
                min(surf.organic[gi, j_surf] / 130.0, 1.0) : 0.0

            # Mineral soil properties (Clapp-Hornberger)
            watsat_mineral = 0.489 - 0.00126 * sand
            bsw_mineral = 2.91 + 0.159 * clay
            sucsat_mineral = 10.0 * 10.0^(1.88 - 0.0131 * sand)
            xksat_mineral = 0.0070556 * 10.0^(-0.884 + 0.0153 * sand)

            # Organic soil properties
            watsat_organic = max(0.93 - 0.1 * (0.0 / 0.3), 0.83)
            bsw_organic = max(12.14 - 7.55 * (0.0 / 0.3), 2.7)
            sucsat_organic = min(10.3 - 0.098 * (0.0 / 0.3), 10.3)
            xksat_organic = max(0.28 - 0.2799 * (0.0 / 0.3), 0.0001)

            # Blend mineral + organic
            watsat_val = (1.0 - om_frac) * watsat_mineral + om_frac * watsat_organic
            bsw_val = (1.0 - om_frac) * bsw_mineral + om_frac * bsw_organic
            sucsat_val = (1.0 - om_frac) * sucsat_mineral + om_frac * sucsat_organic
            hksat_val = (1.0 - om_frac) * xksat_mineral + om_frac * xksat_organic

            ss.watsat_col[c, j] = watsat_val
            ss.bsw_col[c, j] = bsw_val
            ss.sucsat_col[c, j] = sucsat_val
            ss.hksat_col[c, j] = hksat_val

            # Thermal conductivity and heat capacity of dry soil
            ss.tkmg_col[c, j] = TKWAT^watsat_val * TKICE^(1.0 - watsat_val)
            ss.tksatu_col[c, j] = ss.tkmg_col[c, j] * TKICE^(watsat_val - watsat_val)
            ss.tkdry_col[c, j] = (0.135 * (1.0 - watsat_val) * 2700.0 + 64.7) /
                                  (2700.0 - 0.947 * (1.0 - watsat_val) * 2700.0)
            ss.csol_col[c, j] = (1.0 - watsat_val) * (2.0e6 * (1.0 - om_frac) +
                                 2.5e6 * om_frac)
            # Soil water retention (Clapp-Hornberger watdry, watopt, watfc)
            ss.watdry_col[c, j] = watsat_val * (316230.0 / sucsat_val)^(-1.0 / bsw_val)
            ss.watopt_col[c, j] = watsat_val * (158490.0 / sucsat_val)^(-1.0 / bsw_val)
            ss.watfc_col[c, j]  = watsat_val * (0.01 / hksat_val)^(1.0 / (2.0 * bsw_val + 3.0))
        end
    end

    nothing
end

"""
    init_temperatures!(inst, bounds)

Set initial temperature profiles for cold start.
"""
function init_temperatures!(inst::CLMInstances, bounds::BoundsType)
    temp = inst.temperature
    nlevsno = varpar.nlevsno
    nlevmaxurbgrnd = varpar.nlevmaxurbgrnd
    nlevlak = varpar.nlevlak

    # Estimate mean annual temperature from latitude for cold-start init.
    # At 274 K (previous default), deep soil at high latitudes never freezes,
    # causing excessive winter drainage. Use latitude-based estimate instead.
    lat = length(inst.gridcell.latdeg) > 0 ? inst.gridcell.latdeg[1] : 45.0
    T_INIT = max(250.0, 288.0 - 0.5 * abs(lat))  # rough MAT approximation

    for c in bounds.begc:bounds.endc
        temp.t_grnd_col[c] = T_INIT

        # Soil/snow temperature — set with depth-dependent offset
        # Deep soil converges to mean annual T; surface starts cooler in cold climates
        for j in 1:(nlevsno + nlevmaxurbgrnd)
            temp.t_soisno_col[c, j] = T_INIT
        end

        # Lake temperature
        for j in 1:nlevlak
            temp.t_lake_col[c, j] = T_INIT + 2.0
        end

        # Surface water temperature (must be valid even when frac_h2osfc=0
        # because 0.0 * NaN = NaN in IEEE 754)
        temp.t_h2osfc_col[c] = T_INIT
        temp.t_h2osfc_bef_col[c] = T_INIT
    end

    # Patch-level temperatures
    for p in bounds.begp:bounds.endp
        temp.t_veg_patch[p] = T_INIT
        temp.t_stem_patch[p] = T_INIT
        temp.t_skin_patch[p] = T_INIT
        temp.t_ref2m_patch[p] = T_INIT
    end

    nothing
end

"""
    init_soil_moisture!(inst, bounds)

Set initial soil moisture profiles for cold start.
"""
function init_soil_moisture!(inst::CLMInstances, bounds::BoundsType;
                              use_aquifer_layer::Bool = true)
    wsb = inst.water.waterstatebulk_inst
    ws = wsb.ws
    col = inst.column
    ss = inst.soilstate
    nlevsoi = varpar.nlevsoi
    nlevsno = varpar.nlevsno
    joff = nlevsno

    # Initialize soil moisture with depth-dependent profile.
    # Shallow layers get recharged by snowmelt every spring, so their init
    # value matters less. Deep layers (j>5) determine long-term btran and take
    # decades to equilibrate. Start them near residual to match Fortran equilibrium.
    VOL_INIT_SHALLOW = 0.08  # top 5 layers
    VOL_INIT_DEEP = 0.01     # layers 6+ (near residual water content)

    nlevtot = nlevsno + varpar.nlevmaxurbgrnd

    nlevgrnd = varpar.nlevgrnd

    for c in bounds.begc:bounds.endc
        ws.h2osno_no_layers_col[c] = 0.0
        ws.h2osfc_col[c] = 0.0

        # Initialize ALL layers (snow + soil + bedrock) to zero first
        for jj in 1:nlevtot
            ws.excess_ice_col[c, jj] = 0.0
            ws.h2osoi_liq_col[c, jj] = 0.0
            ws.h2osoi_ice_col[c, jj] = 0.0
        end

        # Set soil layers to initial volumetric water content
        for j in 1:nlevsoi
            vol_init = j <= 5 ? VOL_INIT_SHALLOW : VOL_INIT_DEEP
            vol_liq = min(vol_init, ss.watsat_col[c, j])
            ws.h2osoi_vol_col[c, j] = vol_liq
            ws.h2osoi_liq_col[c, j + joff] = vol_liq * col.dz[c, j + joff] * DENH2O
            ws.h2osoi_ice_col[c, j + joff] = 0.0
        end

        # Initialize h2osoi_vol for bedrock layers (nlevsoi+1:nlevgrnd) to zero
        for j in (nlevsoi + 1):nlevgrnd
            ws.h2osoi_vol_col[c, j] = 0.0
        end

        # Keep cold-start aquifer state consistent with WaterStateType%InitCold.
        # Fortran sets wa = aquifer_water_baseline (5000 mm) for both cases.
        ws.wa_col[c] = AQUIFER_WATER_BASELINE

        # Dynbal baselines (zero for cold start)
        ws.dynbal_baseline_liq_col[c] = 0.0
        ws.dynbal_baseline_ice_col[c] = 0.0
    end

    # Canopy water at patch level
    for p in bounds.begp:bounds.endp
        ws.liqcan_patch[p] = 0.0
        ws.snocan_patch[p] = 0.0
    end

    nothing
end

"""
    init_snow_state!(inst, bounds)

Initialize snow state for cold start (no snow).
"""
function init_snow_state!(inst::CLMInstances, bounds::BoundsType)
    col = inst.column

    for c in bounds.begc:bounds.endc
        col.snl[c] = 0
    end

    nothing
end

"""
    init_water_diagnostic_state!(inst, bounds)

Initialize water diagnostic state for cold start. Sets snow fractions,
surface water fractions, and related diagnostic fields to zero/defaults.
These must be valid before the first call to biogeophys_pre_flux_calcs!
and calculate_surface_humidity!.
"""
function init_water_diagnostic_state!(inst::CLMInstances, bounds::BoundsType)
    wdb = inst.water.waterdiagnosticbulk_inst
    nlevsno = varpar.nlevsno
    nlevgrnd = varpar.nlevgrnd

    for c in bounds.begc:bounds.endc
        # Snow cover fractions (no snow at cold start)
        wdb.frac_sno_col[c] = 0.0
        wdb.frac_sno_eff_col[c] = 0.0

        # Surface water fractions (no surface water at cold start)
        wdb.frac_h2osfc_col[c] = 0.0
        wdb.frac_h2osfc_nosnow_col[c] = 0.0

        # Snow diagnostics
        wdb.snow_depth_col[c] = 0.0
        wdb.h2osno_total_col[c] = 0.0
        wdb.snowice_col[c] = 0.0
        wdb.snowliq_col[c] = 0.0

        # Ground saturation humidity (will be recomputed by calculate_surface_humidity!)
        wdb.qg_col[c] = 0.0
        wdb.qg_snow_col[c] = 0.0
        wdb.qg_soil_col[c] = 0.0
        wdb.qg_h2osfc_col[c] = 0.0
        wdb.dqgdT_col[c] = 0.0

        # Snow grain radius defaults
        wdb.snw_rds_top_col[c] = 54.526
        wdb.h2osno_top_col[c] = 0.0
        wdb.sno_liq_top_col[c] = 0.0

        # Ice fraction old (needed by clm_drv_init!)
        for j in 1:nlevsno
            wdb.frac_iceold_col[c, j] = 0.0
        end

        # Snow layer fields
        for j in 1:(nlevsno + nlevgrnd)
            if j <= size(wdb.bw_col, 2)
                wdb.bw_col[c, j] = 0.0
            end
        end

        # Snow grain radius per layer (use fresh snow minimum = 54.526 microns)
        for j in 1:size(wdb.snw_rds_col, 2)
            wdb.snw_rds_col[c, j] = 54.526
        end
    end

    # Column-level plant stored water
    for c in bounds.begc:bounds.endc
        wdb.total_plant_stored_h2o_col[c] = 0.0
    end

    # Patch-level water diagnostics
    for p in bounds.begp:bounds.endp
        wdb.fwet_patch[p] = 0.0
        wdb.fdry_patch[p] = 1.0
        wdb.fcansno_patch[p] = 0.0
        wdb.h2ocan_patch[p] = 0.0
    end

    nothing
end

"""
    init_soil_state_defaults!(inst, bounds)

Initialize soil state fields that are not set by init_soil_properties!
but are needed by surface humidity and flux calculations.
"""
function init_soil_state_defaults!(inst::CLMInstances, bounds::BoundsType)
    ss = inst.soilstate

    col = inst.column
    surf = inst.surfdata

    for c in bounds.begc:bounds.endc
        # Minimum soil matric potential (mm) — very negative limit
        ss.smpmin_col[c] = -1.0e8

        # Soil evaporation beta factor (1.0 = no stress initially)
        ss.soilbeta_col[c] = 1.0

        # Soil evaporation alpha
        ss.soilalpha_col[c] = 1.0
        ss.soilalpha_u_col[c] = 1.0

        # Maximum saturated fraction (TOPMODEL wtfact) — from surfdata FMAX
        if surf !== nothing && length(surf.fmax) > 0
            g = col.gridcell[c]
            ss.wtfact_col[c] = surf.fmax[min(g, length(surf.fmax))]
        else
            ss.wtfact_col[c] = 0.38  # typical default
        end
    end

    nothing
end

"""
    init_satellite_phenology!(inst, bounds, surf)

Initialize satellite phenology (monthly LAI/SAI) from surface data.
"""
function init_satellite_phenology!(inst::CLMInstances, bounds::BoundsType,
                                    surf::SurfaceInputData)
    sp = inst.satellite_phenology
    cs = inst.canopystate
    pch = inst.patch

    for p in bounds.begp:bounds.endp
        g = pch.gridcell[p]
        gi = g - bounds.begg + 1
        ptype = pch.itype[p]

        # PFT type to natural PFT index: ptype=0 (bare ground) → index 1
        m = ptype - varpar.natpft_lb + 1

        if m >= 1 && m <= size(surf.monthly_lai, 2)
            lai_val = surf.monthly_lai[gi, m, 1]
            sai_val = surf.monthly_sai[gi, m, 1]
            htop_val = surf.monthly_htop[gi, m, 1]
            hbot_val = surf.monthly_hbot[gi, m, 1]

            # Set interpolation arrays (month 1 and 2)
            sp.mlai2t[p, 1] = lai_val
            sp.mlai2t[p, 2] = lai_val
            sp.msai2t[p, 1] = sai_val
            sp.msai2t[p, 2] = sai_val
            sp.mhvt2t[p, 1] = htop_val
            sp.mhvt2t[p, 2] = htop_val
            sp.mhvb2t[p, 1] = hbot_val
            sp.mhvb2t[p, 2] = hbot_val

            # Set current LAI/SAI on canopy state
            cs.tlai_patch[p] = lai_val
            cs.tsai_patch[p] = sai_val
            cs.htop_patch[p] = max(htop_val, 0.01)
            cs.hbot_patch[p] = max(hbot_val, 0.01)
            cs.elai_patch[p] = lai_val
            cs.esai_patch[p] = sai_val
        else
            sp.mlai2t[p, 1] = 0.0
            sp.mlai2t[p, 2] = 0.0
            sp.msai2t[p, 1] = 0.0
            sp.msai2t[p, 2] = 0.0
            sp.mhvt2t[p, 1] = 0.0
            sp.mhvt2t[p, 2] = 0.0
            sp.mhvb2t[p, 1] = 0.0
            sp.mhvb2t[p, 2] = 0.0
            cs.tlai_patch[p] = 0.0
            cs.tsai_patch[p] = 0.0
            cs.htop_patch[p] = 0.01
            cs.hbot_patch[p] = 0.01
            cs.elai_patch[p] = 0.0
            cs.esai_patch[p] = 0.0
        end

        # Set frac_veg_nosno_alb (matches satellite_phenology!)
        if (cs.elai_patch[p] + cs.esai_patch[p]) >= 0.05
            cs.frac_veg_nosno_alb_patch[p] = 1
        else
            cs.frac_veg_nosno_alb_patch[p] = 0
        end
    end

    nothing
end

"""
    init_misc_state!(inst, bounds)

Initialize miscellaneous state variables for cold start.
"""
function init_misc_state!(inst::CLMInstances, bounds::BoundsType)
    col = inst.column

    # Convective boundary layer height
    for c in bounds.begc:bounds.endc
        col.zii[c] = 1000.0
    end

    # Canopy state
    cs = inst.canopystate
    for p in bounds.begp:bounds.endp
        cs.laisun_patch[p] = 0.0
        cs.laisha_patch[p] = 0.0
    end

    # Wetness fractions on water diagnostics
    wdb = inst.water.waterdiagnosticbulk_inst
    for p in bounds.begp:bounds.endp
        wdb.fwet_patch[p] = 0.0
        wdb.fdry_patch[p] = 1.0
    end

    # Friction velocity defaults
    fv = inst.frictionvel
    for p in bounds.begp:bounds.endp
        fv.forc_hgt_u_patch[p] = 30.0
        fv.forc_hgt_t_patch[p] = 30.0
        fv.forc_hgt_q_patch[p] = 30.0
    end

    # Surface albedo cold-start defaults (needed before first surface_radiation! call)
    surfalb_init_cold!(inst.surfalb,
                       bounds.begc:bounds.endc,
                       bounds.begp:bounds.endp)

    # Emissivity cold-start
    temp = inst.temperature
    for c in bounds.begc:bounds.endc
        temp.emg_col[c] = 0.96
    end

    # Energy flux cold-start (zero-initialize critical flux fields)
    ef = inst.energyflux
    for p in bounds.begp:bounds.endp
        ef.eflx_sh_grnd_patch[p] = 0.0
        ef.eflx_sh_veg_patch[p] = 0.0
        ef.eflx_sh_tot_patch[p] = 0.0
        ef.eflx_lh_tot_patch[p] = 0.0
        ef.eflx_lh_grnd_patch[p] = 0.0
        ef.eflx_lh_vege_patch[p] = 0.0
        ef.eflx_lh_vegt_patch[p] = 0.0
        ef.eflx_soil_grnd_patch[p] = 0.0
        ef.eflx_lwrad_out_patch[p] = 0.0
        ef.eflx_lwrad_net_patch[p] = 0.0
        ef.eflx_gnet_patch[p] = 0.0
        ef.cgrnd_patch[p] = 0.0
        ef.cgrnds_patch[p] = 0.0
        ef.cgrndl_patch[p] = 0.0
        ef.dlrad_patch[p] = 0.0
        ef.ulrad_patch[p] = 0.0
        ef.dgnetdT_patch[p] = 0.0
        ef.netrad_patch[p] = 0.0
        ef.taux_patch[p] = 0.0
        ef.tauy_patch[p] = 0.0
        ef.btran_patch[p] = 0.0
        ef.btran_min_patch[p] = 1.0
        ef.btran_min_inst_patch[p] = 1.0
        ef.dhsdt_canopy_patch[p] = 0.0
        ef.errsoi_patch[p] = 0.0
        ef.errseb_patch[p] = 0.0
        ef.errsol_patch[p] = 0.0
        ef.errlon_patch[p] = 0.0
    end
    for c in bounds.begc:bounds.endc
        ef.eflx_bot_col[c] = 0.0
        ef.eflx_fgr12_col[c] = 0.0
        ef.eflx_snomelt_col[c] = 0.0
        ef.eflx_sh_precip_conversion_col[c] = 0.0
        ef.errsoi_col[c] = 0.0
        ef.errseb_col[c] = 0.0
        ef.errsol_col[c] = 0.0
        ef.errlon_col[c] = 0.0
    end
    for c in bounds.begc:bounds.endc
        for j in 1:varpar.nlevgrnd
            ef.eflx_fgr_col[c, j] = 0.0
        end
    end

    # Solar absorbed cold-start defaults
    sa = inst.solarabs
    for p in bounds.begp:bounds.endp
        sa.fsa_patch[p] = 0.0
        sa.sabv_patch[p] = 0.0
        sa.sabg_patch[p] = 0.0
        sa.sabg_soil_patch[p] = 0.0
        sa.sabg_snow_patch[p] = 0.0
        for j in 1:size(sa.sabg_lyr_patch, 2)
            sa.sabg_lyr_patch[p, j] = 0.0
        end
    end

    # Water state cold-start: canopy water
    ws_bulk = inst.water.waterstatebulk_inst.ws
    for p in bounds.begp:bounds.endp
        ws_bulk.snocan_patch[p] = 0.0
        ws_bulk.liqcan_patch[p] = 0.0
    end

    # Surface resistance cold-start
    if hasproperty(inst, :surfresis)
        sr = inst.surfresis
        for p in bounds.begp:bounds.endp
            if hasproperty(sr, :rhol_patch)
                sr.rhol_patch[p] = 0.0
            end
        end
    end

    nothing
end

# =========================================================================
# init_soil_hydrology_cold!
# =========================================================================

"""
    init_soil_hydrology_cold!(inst, bounds)

Initialize soil hydrology cold-start state (hkdepth, h2osfc_thresh, zwt,
zwt_perched, frost_table). Calls `soilhydrology_init_cold!` with column
metadata from the CLMInstances container.

Corresponds to SoilHydrologyInitTimeConstMod + SoilHydrologyType%InitCold.
"""
function init_soil_hydrology_cold!(inst::CLMInstances, bounds::BoundsType;
                                   use_aquifer_layer::Bool = true)
    col = inst.column
    lun = inst.landunit
    ws  = inst.water.waterstatebulk_inst.ws

    soilhydrology_init_cold!(
        inst.soilhydrology, bounds.begc:bounds.endc;
        use_aquifer_layer=use_aquifer_layer,
        zi_col       = col.zi,
        wa_col       = ws.wa_col,
        landunit_col = col.landunit,
        lakpoi       = BitVector(lun.lakpoi),
        urbpoi       = BitVector(lun.urbpoi),
        itype_col    = col.itype,
        nbedrock_col = col.nbedrock
    )
end

"""
    init_lake_state!(inst, bounds)

Initialize lake state variables for cold start. Uses `lakestate_init_cold!`
for ice fraction, savedtke1, and ust_lake, then sets remaining fields that
are otherwise left as NaN/SPVAL: betaprime, lakefetch, etal, lake_raw,
ks, ws, lake_icefracsurf, lake_icethick, lakeresist, ram1_lake.

Corresponds to LakeStateType%InitCold + additional defaults for fields
that Fortran initializes via restart reads.
"""
function init_lake_state!(inst::CLMInstances, bounds::BoundsType)
    ls  = inst.lakestate
    col = inst.column
    lun = inst.landunit
    pch = inst.patch

    bc_col   = bounds.begc:bounds.endc
    bc_patch = bounds.begp:bounds.endp

    # Build lake mask
    lake_mask = falses(bounds.endc)
    for c in bc_col
        l = col.landunit[c]
        if l >= 1 && l <= length(lun.lakpoi) && lun.lakpoi[l]
            lake_mask[c] = true
        end
    end

    # Call existing init_cold for icefrac, savedtke1, ust_lake
    lakestate_init_cold!(ls, bc_col; mask_lake=lake_mask)

    # Initialize remaining column-level fields for lake columns
    for c in bc_col
        lake_mask[c] || continue

        ls.betaprime_col[c]        = 1.0
        ls.lake_icefracsurf_col[c] = 0.0
        ls.lake_icethick_col[c]    = 0.0
        ls.lakeresist_col[c]       = 500.0
        ls.lakefetch_col[c]        = 25.0
        ls.etal_col[c]             = 1.0
        ls.lake_raw_col[c]         = 0.0
        ls.ks_col[c]               = 0.0
        ls.ws_col[c]               = 0.0
    end

    # Initialize patch-level fields for lake patches
    for p in bc_patch
        l = pch.landunit[p]
        if l >= 1 && l <= length(lun.lakpoi) && lun.lakpoi[l]
            ls.ram1_lake_patch[p] = 0.0
        end
    end

    # Initialize non-lake columns to safe defaults (zero instead of NaN)
    for c in bc_col
        lake_mask[c] && continue
        ls.betaprime_col[c]        = 0.0
        ls.lake_icefracsurf_col[c] = 0.0
        ls.lake_icethick_col[c]    = 0.0
        ls.lakeresist_col[c]       = 0.0
        ls.lakefetch_col[c]        = 0.0
        ls.etal_col[c]             = 0.0
        ls.lake_raw_col[c]         = 0.0
        ls.ks_col[c]               = 0.0
        ls.ws_col[c]               = 0.0
        for j in 1:size(ls.lake_icefrac_col, 2)
            ls.lake_icefrac_col[c, j] = 0.0
        end
        ls.savedtke1_col[c]        = 0.0
        ls.ust_lake_col[c]         = 0.0
    end
    for p in bc_patch
        l = pch.landunit[p]
        if !(l >= 1 && l <= length(lun.lakpoi) && lun.lakpoi[l])
            ls.ram1_lake_patch[p] = 0.0
        end
    end

    return nothing
end

# =========================================================================
# init_root_fractions!
# =========================================================================

"""
    init_root_fractions!(inst, bounds)

Initialize root fraction profiles (rootfr_patch for water, crootfr_patch for
carbon) using the Zeng (2001) method. Must be called after vertical structure
and soil properties are initialized.

Corresponds to SoilStateType%InitTimeConst → init_vegrootfr.
"""
function init_root_fractions!(inst::CLMInstances, bounds::BoundsType)
    ss  = inst.soilstate
    col = inst.column
    pch = inst.patch
    nlevsno  = varpar.nlevsno
    nlevsoi  = varpar.nlevsoi
    nlevgrnd = varpar.nlevgrnd
    joff     = nlevsno

    # Create soil-only copies of column vertical structure
    # col.zi layout: [snow interfaces..., surface(=0), soil interfaces...]
    #   col.zi[c, joff+1] = 0.0 (surface), col.zi[c, joff+1+j] = bottom of soil layer j
    # init_vegrootfr! expects: col_zi[c, lev] = bottom interface of soil layer lev
    col_zi_soil = Matrix(col.zi[:, (joff + 2):(joff + 1 + nlevgrnd)])
    col_z_soil  = Matrix(col.z[:,  (joff + 1):(joff + nlevgrnd)])
    col_dz_soil = Matrix(col.dz[:, (joff + 1):(joff + nlevgrnd)])

    patch_is_fates = fill(false, bounds.endp)
    bp = bounds.begp:bounds.endp

    # Water root fractions
    init_vegrootfr!(ss.rootfr_patch, col_zi_soil, col_z_soil, col_dz_soil,
                    col.nbedrock, pch.column, pch.itype, patch_is_fates,
                    pftcon, rooting_profile_config, bp, nlevsoi, nlevgrnd, "water")

    # Carbon root fractions
    init_vegrootfr!(ss.crootfr_patch, col_zi_soil, col_z_soil, col_dz_soil,
                    col.nbedrock, pch.column, pch.itype, patch_is_fates,
                    pftcon, rooting_profile_config, bp, nlevsoi, nlevgrnd, "carbon")

    return nothing
end

# =========================================================================
# init_eff_porosity!
# =========================================================================

"""
    init_eff_porosity!(inst, bounds)

Initialize effective porosity (watsat - vol_ice) for all soil layers.
Must be called after init_soil_properties! (sets watsat) and init_soil_moisture!
(sets h2osoi_ice). Needed before first call to calc_root_moist_stress!.
"""
function init_eff_porosity!(inst::CLMInstances, bounds::BoundsType)
    col = inst.column
    ss  = inst.soilstate
    wsb = inst.water.waterstatebulk_inst.ws
    nlevsoi  = varpar.nlevsoi
    nlevsno  = varpar.nlevsno
    joff     = nlevsno

    for j in 1:nlevsoi
        for c in bounds.begc:bounds.endc
            dz_val = col.dz[c, j + joff]
            if dz_val > 0.0
                vol_ice = min(ss.watsat_col[c, j],
                              wsb.h2osoi_ice_col[c, j + joff] / (dz_val * DENICE))
            else
                vol_ice = 0.0
            end
            ss.eff_porosity_col[c, j] = max(0.01, ss.watsat_col[c, j] - vol_ice)
        end
    end

    return nothing
end
