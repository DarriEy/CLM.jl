# ==========================================================================
# Ported from: src/biogeochem/CNVegStructUpdateMod.F90
# Vegetation structure updates (LAI, SAI, htop, hbot) from C state variables.
#
# On the radiation time step, uses C state variables and ecological process
# constants (epc) to diagnose vegetation structure (LAI, SAI, height).
#
# Public functions:
#   cn_veg_struct_update! -- Update vegetation structure from live C pools
# ==========================================================================

# ---------------------------------------------------------------------------
# Module-level constants
# ---------------------------------------------------------------------------

const DTSMONTH = 2592000.0            # seconds in a 30-day month (60*60*24*30)
const FRAC_SNO_THRESHOLD = 0.999      # frac_sno values > this treated as 1

# ---------------------------------------------------------------------------
# cn_veg_struct_update! -- Main vegetation structure update
# ---------------------------------------------------------------------------

"""
    cn_veg_struct_update!(mask_soilp, bounds, patch, canopystate,
        cnveg_carbonstate, waterdiagnosticbulk, frictionvel, cnveg_state,
        crop, pftcon_data;
        dt, use_cndv, use_biomass_heat_storage, spinup_factor_deadwood,
        noveg, nc3crop, nc3irrig, nbrdlf_evr_shrub, nbrdlf_dcd_brl_shrub,
        npcropmin, ntmp_corn, nirrig_tmp_corn, ntrp_corn, nirrig_trp_corn,
        nsugarcane, nirrig_sugarcane, nmiscanthus, nirrig_miscanthus,
        nswitchgrass, nirrig_switchgrass, c_to_b)

On the radiation time step, use C state variables and epc to diagnose
vegetation structure (LAI, SAI, height).

Updates `canopystate` fields: `tlai_patch`, `tsai_patch`, `htop_patch`,
`hbot_patch`, `elai_patch`, `esai_patch`, `frac_veg_nosno_alb_patch`,
`stem_biomass_patch`, `leaf_biomass_patch`.

Updates `cnveg_state` fields: `htmx_patch`, `peaklai_patch`.

Ported from `CNVegStructUpdate` in `CNVegStructUpdateMod.F90`.
"""
function cn_veg_struct_update!(mask_soilp::BitVector,
                                bounds::UnitRange{Int},
                                patch::PatchData,
                                canopystate::CanopyStateData,
                                cnveg_carbonstate::CNVegCarbonStateData,
                                waterdiagnosticbulk::WaterDiagnosticBulkData,
                                frictionvel::FrictionVelocityData,
                                cnveg_state::CNVegStateData,
                                crop::CropData,
                                pftcon_data::PftconType;
                                dt::Float64 = 1800.0,
                                use_cndv::Bool = false,
                                use_biomass_heat_storage::Bool = false,
                                spinup_factor_deadwood::Float64 = SPINUP_FACTOR_DEADWOOD_DEFAULT,
                                noveg::Int = CLM.noveg,
                                nc3crop::Int = CLM.nc3crop,
                                nc3irrig::Int = CLM.nc3irrig,
                                nbrdlf_evr_shrub::Int = CLM.nbrdlf_evr_shrub,
                                nbrdlf_dcd_brl_shrub::Int = CLM.nbrdlf_dcd_brl_shrub,
                                npcropmin::Int = CLM.npcropmin,
                                ntmp_corn::Int = CLM.ntmp_corn,
                                nirrig_tmp_corn::Int = CLM.nirrig_tmp_corn,
                                ntrp_corn::Int = CLM.ntrp_corn,
                                nirrig_trp_corn::Int = CLM.nirrig_trp_corn,
                                nsugarcane::Int = CLM.nsugarcane,
                                nirrig_sugarcane::Int = CLM.nirrig_sugarcane,
                                nmiscanthus::Int = CLM.nmiscanthus,
                                nirrig_miscanthus::Int = CLM.nirrig_miscanthus,
                                nswitchgrass::Int = CLM.nswitchgrass,
                                nirrig_switchgrass::Int = CLM.nirrig_switchgrass,
                                c_to_b::Float64 = C_TO_B)

    # --- Aliases (matching Fortran associate block) ---
    ivt            = patch.itype
    col            = patch.column

    woody          = pftcon_data.woody
    slatop         = pftcon_data.slatop
    dsladlai       = pftcon_data.dsladlai
    z0mr           = pftcon_data.z0mr
    displar        = pftcon_data.displar
    dwood          = pftcon_data.dwood
    ztopmx         = pftcon_data.ztopmx
    laimx          = pftcon_data.laimx
    nstem          = pftcon_data.nstem
    taper          = pftcon_data.taper
    fbw            = pftcon_data.fbw

    frac_sno       = waterdiagnosticbulk.frac_sno_col
    snow_depth     = waterdiagnosticbulk.snow_depth_col

    forc_hgt_u_patch = frictionvel.forc_hgt_u_patch

    leafc          = cnveg_carbonstate.leafc_patch
    deadstemc      = cnveg_carbonstate.deadstemc_patch
    livestemc      = cnveg_carbonstate.livestemc_patch

    farea_burned   = cnveg_state.farea_burned_col
    htmx           = cnveg_state.htmx_patch
    peaklai        = cnveg_state.peaklai_patch

    harvdate       = crop.harvdate_patch

    tlai           = canopystate.tlai_patch
    tsai           = canopystate.tsai_patch
    stem_biomass   = canopystate.stem_biomass_patch
    leaf_biomass   = canopystate.leaf_biomass_patch
    htop           = canopystate.htop_patch
    hbot           = canopystate.hbot_patch
    elai           = canopystate.elai_patch
    esai           = canopystate.esai_patch
    frac_veg_nosno_alb = canopystate.frac_veg_nosno_alb_patch

    # --- Patch loop ---
    for p in bounds
        mask_soilp[p] || continue

        c = col[p]

        if ivt[p] != noveg

            tlai_old = tlai[p]  # n-1 value
            tsai_old = tsai[p]  # n-1 value

            # Update the leaf area index based on leafC and SLA
            # Eq 3 from Thornton and Zimmerman, 2007, J Clim, 20, 3902-3923.
            if dsladlai[ivt[p]] > 0.0
                tlai[p] = (slatop[ivt[p]] * (exp(leafc[p] * dsladlai[ivt[p]]) - 1.0)) / dsladlai[ivt[p]]
            else
                tlai[p] = slatop[ivt[p]] * leafc[p]
            end
            tlai[p] = max(0.0, tlai[p])

            # Update the stem area index and height based on LAI, stem mass, and veg type.
            # tsai formula from Zeng et al. 2002, Journal of Climate, p1835
            # Assumes doalb time step .eq. CLM time step, SAI min and monthly decay factor
            # alpha are set by PFT, and alpha is scaled to CLM time step
            # tsai_min scaled by 0.5 to match MODIS satellite derived values
            if ivt[p] == nc3crop || ivt[p] == nc3irrig  # generic crops
                tsai_alpha = 1.0 - 1.0 * dt / DTSMONTH
                tsai_min = 0.1
            else
                tsai_alpha = 1.0 - 0.5 * dt / DTSMONTH
                tsai_min = 1.0
            end
            tsai_min = tsai_min * 0.5
            tsai[p] = max(tsai_alpha * tsai_old + max(tlai_old - tlai[p], 0.0), tsai_min)

            # Calculate vegetation physiological parameters used in biomass heat storage
            if use_biomass_heat_storage
                # Assumes fbw (fraction of biomass that is water) is the same for leaves and stems
                leaf_biomass[p] = max(0.0025, leafc[p]) *
                    c_to_b * 1.0e-3 / (1.0 - fbw[ivt[p]])
            end

            if woody[ivt[p]] == 1.0

                # Trees and shrubs: simple allometry with hard-wired stem taper
                # and nstem from PFT parameter file
                if use_cndv
                    # CNDV pathway (not commonly used)
                    # Placeholder: would use dgv_ecophyscon allom2/allom3 and nind/fpcgrid
                    # For now, use standard non-CNDV allometry
                    htop[p] = ((3.0 * deadstemc[p] * spinup_factor_deadwood * taper[ivt[p]] * taper[ivt[p]]) /
                        (pi * nstem[ivt[p]] * dwood[ivt[p]]))^(1.0 / 3.0)
                else
                    # correct height calculation if doing accelerated spinup
                    htop[p] = ((3.0 * deadstemc[p] * spinup_factor_deadwood * taper[ivt[p]] * taper[ivt[p]]) /
                        (pi * nstem[ivt[p]] * dwood[ivt[p]]))^(1.0 / 3.0)
                end

                if use_biomass_heat_storage
                    # Assumes fbw (fraction of biomass that is water) is the same for leaves and stems
                    stem_biomass[p] = (spinup_factor_deadwood * deadstemc[p] + livestemc[p]) *
                        c_to_b * 1.0e-3 / (1.0 - fbw[ivt[p]])
                end

                # Keep htop from getting too close to forcing height for windspeed
                # (Peter Thornton, 5/3/2004)
                htop[p] = min(htop[p], (forc_hgt_u_patch[p] / (displar[ivt[p]] + z0mr[ivt[p]])) - 3.0)

                # Constraint to keep htop from going to 0.0
                # (Peter Thornton, 8/11/2004)
                htop[p] = max(htop[p], 0.01)

                hbot[p] = max(0.0, min(3.0, htop[p] - 1.0))

            elseif ivt[p] >= npcropmin  # prognostic crops

                if tlai[p] >= laimx[ivt[p]]
                    peaklai[p] = 1  # used in CNAllocation
                end

                if ivt[p] == ntmp_corn || ivt[p] == nirrig_tmp_corn ||
                   ivt[p] == ntrp_corn || ivt[p] == nirrig_trp_corn ||
                   ivt[p] == nsugarcane || ivt[p] == nirrig_sugarcane ||
                   ivt[p] == nmiscanthus || ivt[p] == nirrig_miscanthus ||
                   ivt[p] == nswitchgrass || ivt[p] == nirrig_switchgrass
                    tsai[p] = 0.1 * tlai[p]
                else
                    tsai[p] = 0.2 * tlai[p]
                end

                # "stubble" after harvest
                if harvdate[p] < 999 && tlai[p] == 0.0
                    tsai[p] = 0.25 * (1.0 - farea_burned[c] * 0.90)
                    htmx[p] = 0.0
                    peaklai[p] = 0
                end

                # canopy top and bottom heights
                htop[p] = ztopmx[ivt[p]] * (min(tlai[p] / (laimx[ivt[p]] - 1.0), 1.0))^2
                htmx[p] = max(htmx[p], htop[p])
                htop[p] = max(0.05, max(htmx[p], htop[p]))
                hbot[p] = 0.02

            else  # generic crops and grasses

                # height for grasses depends only on LAI
                htop[p] = max(0.25, tlai[p] * 0.25)

                htop[p] = min(htop[p], (forc_hgt_u_patch[p] / (displar[ivt[p]] + z0mr[ivt[p]])) - 3.0)

                # Constraint to keep htop from going to 0.0
                htop[p] = max(htop[p], 0.01)

                hbot[p] = max(0.0, min(0.05, htop[p] - 0.20))
            end

        else
            # noveg
            tlai[p] = 0.0
            tsai[p] = 0.0
            htop[p] = 0.0
            hbot[p] = 0.0
        end

        # Adjust lai and sai for burying by snow.
        # Snow burial fraction for short vegetation (e.g. grasses, crops) changes
        # with vegetation height. Accounts for a 20% bending factor, as used in
        # Lombardozzi et al. (2018) GRL 45(18), 9889-9897.
        #
        # NOTE: The following snow burial code is duplicated in SatellitePhenologyMod.
        # Changes in one place should be accompanied by similar changes in the other.

        if ivt[p] > noveg && ivt[p] <= nbrdlf_dcd_brl_shrub
            ol = min(max(snow_depth[c] - hbot[p], 0.0), htop[p] - hbot[p])
            fb = 1.0 - ol / max(1.0e-06, htop[p] - hbot[p])
        else
            fb = 1.0 - (max(min(snow_depth[c], max(0.05, htop[p] * 0.8)), 0.0) /
                        (max(0.05, htop[p] * 0.8)))
            # depth of snow required for complete burial of grasses
        end

        if frac_sno[c] <= FRAC_SNO_THRESHOLD
            frac_sno_adjusted = frac_sno[c]
        else
            # avoid tiny but non-zero elai and esai that can cause radiation
            # and/or photosynthesis code to blow up
            frac_sno_adjusted = 1.0
        end

        elai[p] = max(tlai[p] * (1.0 - frac_sno_adjusted) + tlai[p] * fb * frac_sno_adjusted, 0.0)
        esai[p] = max(tsai[p] * (1.0 - frac_sno_adjusted) + tsai[p] * fb * frac_sno_adjusted, 0.0)

        # Fraction of vegetation free of snow
        if (elai[p] + esai[p]) > 0.0
            frac_veg_nosno_alb[p] = 1
        else
            frac_veg_nosno_alb[p] = 0
        end

    end  # patch loop

    return nothing
end
