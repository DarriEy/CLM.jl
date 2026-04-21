# ==========================================================================
# Ported from:
#   src/biogeochem/CNDVType.F90        (types + accumulator updates)
#   src/biogeochem/CNDVLightMod.F90    (light competition)
#   src/biogeochem/CNDVEstablishmentMod.F90  (establishment + mortality)
#   src/biogeochem/CNDVDriverMod.F90   (annual CNDV driver)
#   src/biogeochem/dynCNDVMod.F90      (weight interpolation)
#
# Dynamic Global Vegetation Model (CNDV) — the coupled CN-DGVM mode
# that computes annual establishment, light competition, and mortality,
# then interpolates patch weights to every time step.
#
# Public functions:
#   cndv_light!            — annual light competition
#   cndv_establishment!    — annual establishment + mortality
#   cndv_driver!           — annual CNDV driver (climate20 + light + estab)
#   dyn_cndv_init!         — initialize fpcgrid from wtcol at startup
#   dyn_cndv_interp!       — interpolate wtcol from annual fpcgrid
#   cndv_update_acc_vars!  — accumulate AGDDTW and AGDD every time step
# ==========================================================================

# Types DGVSData, DGVEcophysCon and their init functions are defined in
# src/types/dgvs.jl (dgvs_init!, dgvs_init_cold!, dgv_ecophyscon_init!)

# =========================================================================
#  Light competition — CNDVLightMod.F90
# =========================================================================

"""
    cndv_light!(dgvs, eco, deadstemc_patch, leafcmax_patch,
                pftcon, patch, mask_natvegp, bounds_patch;
                bounds_gridcell)

Calculate light competition and update FPC for the establishment routine.
Called once per year.

Ported from subroutine Light in CNDVLightMod.F90.

# Arguments
- `bounds_gridcell`: range of gridcell indices (default 1:1)
"""
function cndv_light!(dgvs::DGVSData,
                     eco::DGVEcophysCon,
                     deadstemc_patch::Vector{<:Real},
                     leafcmax_patch::Vector{<:Real},
                     pftcon_in::PftconType,
                     patch::PatchData,
                     mask_natvegp::BitVector,
                     bounds_patch::UnitRange{Int};
                     bounds_gridcell::UnitRange{Int} = 1:1)

    # --- Constants ---
    fpc_tree_max = 0.95
    taper = 200.0

    # --- Aliases ---
    ivt           = patch.itype
    crownarea_max = eco.crownarea_max
    reinickerp_a  = eco.reinickerp
    allom1_a      = eco.allom1
    dwood         = pftcon_in.dwood
    slatop        = pftcon_in.slatop
    dsladlai      = pftcon_in.dsladlai
    woody         = pftcon_in.woody
    is_tree       = pftcon_in.is_tree
    is_shrub      = pftcon_in.is_shrub

    crownarea = dgvs.crownarea_patch
    nind      = dgvs.nind_patch
    fpcgrid   = dgvs.fpcgrid_patch

    ng = length(bounds_gridcell)

    # --- Gridcell-level accumulators ---
    fpc_tree_total  = zeros(last(bounds_gridcell))
    fpc_inc_tree    = zeros(last(bounds_gridcell))
    fpc_grass_total = zeros(last(bounds_gridcell))
    fpc_shrub_total = zeros(last(bounds_gridcell))
    numtrees        = zeros(Int, last(bounds_gridcell))
    fpc_inc         = zeros(last(bounds_patch))

    # --- Pass 1: update LAI/FPC, accumulate gridcell totals ---
    for p in bounds_patch
        mask_natvegp[p] || continue
        g   = patch.gridcell[p]
        ivp = ivt[p] + 1  # +1 for Julia 1-based PFT indexing

        # Update crown area for woody PFTs
        if woody[ivp] == 1.0
            if fpcgrid[p] > 0.0 && nind[p] > 0.0
                stocking = nind[p] / fpcgrid[p]
                stemdiam = (24.0 * deadstemc_patch[p] /
                           (RPI * stocking * dwood[ivp] * taper))^(1.0/3.0)
            else
                stemdiam = 0.0
            end
            crownarea[p] = min(crownarea_max[ivp],
                               allom1_a[ivp] * stemdiam^reinickerp_a[ivp])
        end

        # Compute LAI per individual
        if crownarea[p] > 0.0 && nind[p] > 0.0
            lm_ind = leafcmax_patch[p] * fpcgrid[p] / nind[p]
            if dsladlai[ivp] > 0.0
                lai_ind = max(0.001,
                    (exp(lm_ind * dsladlai[ivp] + log(slatop[ivp])) -
                     slatop[ivp]) / dsladlai[ivp] / crownarea[p])
            else
                lai_ind = lm_ind * slatop[ivp] / crownarea[p]
            end
        else
            lai_ind = 0.0
        end

        fpc_ind = 1.0 - exp(-0.5 * lai_ind)
        fpcgrid_old = fpcgrid[p]
        fpcgrid[p] = crownarea[p] * nind[p] * fpc_ind
        fpc_inc[p] = max(0.0, fpcgrid[p] - fpcgrid_old)

        if woody[ivp] == 1.0
            if is_tree[ivp]
                numtrees[g] += 1
                fpc_tree_total[g] += fpcgrid[p]
                fpc_inc_tree[g] += fpc_inc[p]
            elseif is_shrub[ivp]
                fpc_shrub_total[g] += fpcgrid[p]
            end
        else  # grass
            fpc_grass_total[g] += fpcgrid[p]
        end
    end

    # --- Compute grass/shrub max allowed FPC ---
    fpc_grass_max = zeros(last(bounds_gridcell))
    fpc_shrub_max = zeros(last(bounds_gridcell))
    for g in bounds_gridcell
        fpc_grass_max[g] = 1.0 - min(fpc_tree_total[g], fpc_tree_max)
        fpc_shrub_max[g] = max(0.0, fpc_grass_max[g] - fpc_grass_total[g])
    end

    # --- Pass 2: light competition adjustments ---
    for p in bounds_patch
        mask_natvegp[p] || continue
        g   = patch.gridcell[p]
        ivp = ivt[p] + 1

        if woody[ivp] == 1.0 && is_tree[ivp]
            # Tree light competition
            if fpc_tree_total[g] > fpc_tree_max
                if fpc_inc_tree[g] > 0.0
                    excess = (fpc_tree_total[g] - fpc_tree_max) *
                             fpc_inc[p] / fpc_inc_tree[g]
                else
                    excess = (fpc_tree_total[g] - fpc_tree_max) /
                             Float64(numtrees[g])
                end

                if fpcgrid[p] > 0.0
                    nind_kill = nind[p] * excess / fpcgrid[p]
                    nind[p] = max(0.0, nind[p] - nind_kill)
                    fpcgrid[p] = max(0.0, fpcgrid[p] - excess)
                else
                    nind[p] = 0.0
                    fpcgrid[p] = 0.0
                end
            end

        elseif woody[ivp] == 0.0
            # Grass competition
            if fpc_grass_total[g] > fpc_grass_max[g]
                excess = (fpc_grass_total[g] - fpc_grass_max[g]) *
                         fpcgrid[p] / fpc_grass_total[g]
                fpcgrid[p] = max(0.0, fpcgrid[p] - excess)
            end

        elseif woody[ivp] == 1.0 && is_shrub[ivp]
            # Shrub competition
            if fpc_shrub_total[g] > fpc_shrub_max[g]
                excess = 1.0 - fpc_shrub_max[g] / fpc_shrub_total[g]

                if fpcgrid[p] > 0.0
                    nind_kill = nind[p] * excess / fpcgrid[p]
                    nind[p] = max(0.0, nind[p] - nind_kill)
                    fpcgrid[p] = max(0.0, fpcgrid[p] - excess)
                else
                    nind[p] = 0.0
                    fpcgrid[p] = 0.0
                end
            end
        end
    end

    return nothing
end

# =========================================================================
#  Establishment — CNDVEstablishmentMod.F90
# =========================================================================

"""
    cndv_establishment!(dgvs, eco, prec365_col, annsum_npp_patch,
                         annsum_litfall_patch, deadstemc_patch,
                         leafcmax_patch, pftcon, patch, lun,
                         mask_natvegp, bounds_patch;
                         bounds_gridcell)

Calculate establishment of new patches and annual mortality diagnostics.
Called once per year.

Ported from subroutine Establishment in CNDVEstablishmentMod.F90.
"""
function cndv_establishment!(dgvs::DGVSData,
                              eco::DGVEcophysCon,
                              prec365_col::Vector{<:Real},
                              annsum_npp_patch::Vector{<:Real},
                              annsum_litfall_patch::Vector{<:Real},
                              deadstemc_patch::Vector{<:Real},
                              leafcmax_patch::Vector{<:Real},
                              pftcon_in::PftconType,
                              patch::PatchData,
                              lun::LandunitData,
                              bounds_patch::UnitRange{Int};
                              bounds_gridcell::UnitRange{Int} = 1:1)

    # --- Constants ---
    ramp_agddtw   = 300.0
    nind_min       = 1.0e-10
    prec_min_estab = 100.0 / (365.0 * SECSPDAY)
    estab_max_val  = 0.24
    taper          = 200.0

    # --- Aliases ---
    ivt           = patch.itype
    slatop        = pftcon_in.slatop
    dsladlai      = pftcon_in.dsladlai
    dwood         = pftcon_in.dwood
    woody         = pftcon_in.woody

    crownarea_max = eco.crownarea_max
    twmax         = eco.twmax
    reinickerp_a  = eco.reinickerp
    allom1_a      = eco.allom1
    tcmax         = eco.tcmax
    tcmin         = eco.tcmin
    gddmin_eco    = eco.gddmin

    present    = dgvs.present_patch
    nind       = dgvs.nind_patch
    fpcgrid    = dgvs.fpcgrid_patch
    crownarea  = dgvs.crownarea_patch
    greffic    = dgvs.greffic_patch
    heatstress = dgvs.heatstress_patch
    agddtw     = dgvs.agddtw_patch
    agdd20     = dgvs.agdd20_patch
    tmomin20   = dgvs.tmomin20_patch
    pftmayexist = dgvs.pftmayexist_patch

    np = last(bounds_patch)

    # --- Gridcell-level accumulators ---
    FT = eltype(dgvs.nind_patch)
    ngrass        = zeros(Int, last(bounds_gridcell))
    npft_estab    = zeros(Int, last(bounds_gridcell))
    fpc_tree_total = zeros(FT, last(bounds_gridcell))
    fpc_total     = zeros(FT, last(bounds_gridcell))
    fpc_total_new = zeros(FT, last(bounds_gridcell))

    # --- Patch-level temporaries ---
    survive = fill(false, np)
    estab   = fill(false, np)
    dstemc  = zeros(FT, np)

    # ---------------------------------------------------------------
    # Initialize presence and local copies
    # ---------------------------------------------------------------
    for p in bounds_patch
        if nind[p] == 0.0
            present[p] = false
        end
        if !present[p]
            nind[p]    = 0.0
            fpcgrid[p] = 0.0
        end
        survive[p] = false
        estab[p]   = false
        dstemc[p]  = deadstemc_patch[p]
    end

    # ---------------------------------------------------------------
    # Bioclim: determine survival and establishment
    # ---------------------------------------------------------------
    for p in bounds_patch
        ivp = ivt[p] + 1  # Julia 1-based PFT index
        if tmomin20[p] >= tcmin[ivp] + TFRZ
            if tmomin20[p] <= tcmax[ivp] + TFRZ && agdd20[p] >= gddmin_eco[ivp]
                estab[p] = true
            end
            survive[p] = true
            # Seasonal deciduous patches that would have occurred in
            # regions without short winter day lengths
            if !pftmayexist[p]
                survive[p] = false
                estab[p]   = false
                pftmayexist[p] = true
            end
        end
    end

    # ---------------------------------------------------------------
    # Case 1: kill PFTs not adapted; Case 2: introduce new PFTs
    # ---------------------------------------------------------------
    for p in bounds_patch
        c   = patch.column[p]
        l   = patch.landunit[p]
        ivp = ivt[p] + 1

        # Case 1 — kill
        if present[p] && (!survive[p] || nind[p] < nind_min)
            present[p] = false
            fpcgrid[p] = 0.0
            nind[p]    = 0.0
        end

        # Case 2 — introduce
        if lun.itype[l] == ISTSOIL
            if !present[p] && prec365_col[c] >= prec_min_estab && estab[p]
                if twmax[ivp] > 999.0 || agddtw[p] == 0.0
                    present[p] = true
                    nind[p]    = 0.0
                    fpcgrid[p] = 0.000844
                    if woody[ivp] < 1.0
                        fpcgrid[p] = 0.05
                    end
                    leafcmax_patch[p] = 1.0
                    if dstemc[p] <= 0.0
                        dstemc[p] = 0.1
                    end
                end
            end
        end
    end

    # ---------------------------------------------------------------
    # Calculate total woody FPC and number of woody PFTs able to establish
    # ---------------------------------------------------------------
    for p in bounds_patch
        g   = patch.gridcell[p]
        ivp = ivt[p] + 1
        if present[p]
            if woody[ivp] == 1.0
                fpc_tree_total[g] += fpcgrid[p]
                if estab[p]
                    npft_estab[g] += 1
                end
            elseif woody[ivp] < 1.0 && ivt[p] > noveg
                ngrass[g] += 1
            end
        end
    end

    # ---------------------------------------------------------------
    # Sapling establishment for trees
    # ---------------------------------------------------------------
    for p in bounds_patch
        g   = patch.gridcell[p]
        ivp = ivt[p] + 1

        if present[p] && woody[ivp] == 1.0 && estab[p]
            estab_rate = estab_max_val *
                         (1.0 - exp(5.0 * (fpc_tree_total[g] - 1.0))) /
                         Float64(npft_estab[g])
            estab_grid = estab_rate * (1.0 - fpc_tree_total[g])
            nind[p] += estab_grid

            lm_ind = leafcmax_patch[p] * fpcgrid[p] / nind[p]
            if fpcgrid[p] > 0.0 && nind[p] > 0.0
                stocking = nind[p] / fpcgrid[p]
                stemdiam = (24.0 * dstemc[p] /
                           (RPI * stocking * dwood[ivp] * taper))^(1.0/3.0)
            else
                stemdiam = 0.0
            end
            crownarea[p] = min(crownarea_max[ivp],
                               allom1_a[ivp] * stemdiam^reinickerp_a[ivp])

            # Update LAI and FPC
            if crownarea[p] > 0.0
                if dsladlai[ivp] > 0.0
                    lai_ind = max(0.001,
                        (exp(lm_ind * dsladlai[ivp] + log(slatop[ivp])) -
                         slatop[ivp]) / dsladlai[ivp] / crownarea[p])
                else
                    lai_ind = lm_ind * slatop[ivp] / crownarea[p]
                end
            else
                lai_ind = 0.0
            end

            fpc_ind    = 1.0 - exp(-0.5 * lai_ind)
            fpcgrid[p] = crownarea[p] * nind[p] * fpc_ind
        end

        if present[p] && woody[ivp] == 1.0
            fpc_total_new[g] += fpcgrid[p]
        end
    end

    # ---------------------------------------------------------------
    # Adjustment: don't allow trees to exceed 95%
    # ---------------------------------------------------------------
    for p in bounds_patch
        g = patch.gridcell[p]
        ivp = ivt[p] + 1
        if fpc_total_new[g] > 0.95
            if woody[ivp] == 1.0 && present[p]
                nind[p]    = nind[p] * 0.95 / fpc_total_new[g]
                fpcgrid[p] = fpcgrid[p] * 0.95 / fpc_total_new[g]
            end
            fpc_total[g] = 0.95
        else
            fpc_total[g] = fpc_total_new[g]
        end
    end

    # ---------------------------------------------------------------
    # Grass establishment
    # ---------------------------------------------------------------
    for p in bounds_patch
        g   = patch.gridcell[p]
        ivp = ivt[p] + 1

        if present[p] && woody[ivp] < 1.0
            if leafcmax_patch[p] <= 0.0 || fpcgrid[p] <= 0.0
                present[p] = false
                nind[p]    = 0.0
            else
                nind[p]      = 1.0
                crownarea[p] = 1.0
                lm_ind = leafcmax_patch[p] * fpcgrid[p] / nind[p]
                if dsladlai[ivp] > 0.0
                    lai_ind = max(0.001,
                        (exp(lm_ind * dsladlai[ivp] + log(slatop[ivp])) -
                         slatop[ivp]) / dsladlai[ivp] / crownarea[p])
                else
                    lai_ind = lm_ind * slatop[ivp] / crownarea[p]
                end
                fpc_ind    = 1.0 - exp(-0.5 * lai_ind)
                fpcgrid[p] = crownarea[p] * nind[p] * fpc_ind
                fpc_total[g] += fpcgrid[p]
            end
        end
    end

    # ---------------------------------------------------------------
    # Adjustment: fpc_total > 1 due to grasses (nc3_arctic_grass and above)
    # ---------------------------------------------------------------
    for p in bounds_patch
        g   = patch.gridcell[p]
        ivp = ivt[p] + 1

        if fpc_total[g] > 1.0
            if ivt[p] >= nc3_arctic_grass && fpcgrid[p] > 0.0
                fpcgridtemp = fpcgrid[p]
                fpcgrid[p]  = max(0.0, fpcgrid[p] - (fpc_total[g] - 1.0))
                fpc_total[g] = fpc_total[g] - fpcgridtemp + fpcgrid[p]
            end
        end

        # Remove tiny fpcgrid amounts
        if fpcgrid[p] < 1.0e-15
            fpc_total[g] -= fpcgrid[p]
            fpcgrid[p]   = 0.0
            present[p]   = false
            nind[p]      = 0.0
        end

        # Set bare ground fpcgrid so everything adds to 1
        if fpc_total[g] < 1.0 && ivt[p] == noveg
            fpcgrid[p]   = 1.0 - fpc_total[g]
            fpc_total[g] += fpcgrid[p]
        end
    end

    # ---------------------------------------------------------------
    # Annual mortality diagnostics (heatstress, greffic)
    # Used hourly in gap_mortality
    # ---------------------------------------------------------------
    for p in bounds_patch
        g   = patch.gridcell[p]
        ivp = ivt[p] + 1

        if woody[ivp] == 1.0 && nind[p] > 0.0 &&
           leafcmax_patch[p] > 0.0 && fpcgrid[p] > 0.0

            if twmax[ivp] < 999.0
                heatstress[p] = max(0.0, min(1.0, agddtw[p] / ramp_agddtw))
            else
                heatstress[p] = 0.0
            end

            # Net individual living biomass increment
            bm_delta = max(0.0, annsum_npp_patch[p] - annsum_litfall_patch[p])
            lm_ind   = leafcmax_patch[p] * fpcgrid[p] / nind[p]

            # Growth efficiency (net biomass increment per unit leaf area)
            if dsladlai[ivp] > 0.0
                greffic[p] = bm_delta /
                    max(0.001,
                        (exp(lm_ind * dsladlai[ivp] + log(slatop[ivp])) -
                         slatop[ivp]) / dsladlai[ivp])
            else
                greffic[p] = bm_delta / (lm_ind * slatop[ivp])
            end
        else
            greffic[p]    = 0.0
            heatstress[p] = 0.0
        end
    end

    # ---------------------------------------------------------------
    # Error check: fpc sums to 1 per gridcell
    # ---------------------------------------------------------------
    for g in bounds_gridcell
        if abs(fpc_total[g] - 1.0) > 1.0e-6
            @warn "CNDV Establishment: fpc_total = $(fpc_total[g]) at gridcell $g"
        end
    end

    return nothing
end

# =========================================================================
#  CNDV Driver — CNDVDriverMod.F90
# =========================================================================

"""
    cndv_driver!(dgvs, eco, t_mo_min_patch, prec365_col,
                  annsum_npp_patch, annsum_litfall_patch,
                  deadstemc_patch, leafcmax_patch,
                  pftcon, patch, lun, mask_natvegp,
                  bounds_patch, kyr;
                  bounds_gridcell)

Drive the annual dynamic vegetation cycle:
1. Compute 20-yr running means of climate variables (climate20)
2. Rebuild naturally-vegetated patch filter from present_patch
3. Call `cndv_light!`
4. Call `cndv_establishment!`
5. Reset annual accumulators

Ported from subroutine CNDVDriver in CNDVDriverMod.F90.
"""
function cndv_driver!(dgvs::DGVSData,
                      eco::DGVEcophysCon,
                      t_mo_min_patch::Vector{<:Real},
                      prec365_col::Vector{<:Real},
                      annsum_npp_patch::Vector{<:Real},
                      annsum_litfall_patch::Vector{<:Real},
                      deadstemc_patch::Vector{<:Real},
                      leafcmax_patch::Vector{<:Real},
                      pftcon_in::PftconType,
                      patch::PatchData,
                      lun::LandunitData,
                      mask_natvegp::BitVector,
                      bounds_patch::UnitRange{Int},
                      kyr::Int;
                      bounds_gridcell::UnitRange{Int} = 1:1)

    # --- Aliases ---
    fpcgrid  = dgvs.fpcgrid_patch
    agdd20   = dgvs.agdd20_patch
    tmomin20 = dgvs.tmomin20_patch
    agdd     = dgvs.agdd_patch

    # ---------------------------------------------------------------
    # 1. Climate20: 20-yr running means
    # ---------------------------------------------------------------
    for p in bounds_patch
        if kyr == 2
            tmomin20[p] = t_mo_min_patch[p]
            agdd20[p]   = agdd[p]
        end
        tmomin20[p] = (19.0 * tmomin20[p] + t_mo_min_patch[p]) / 20.0
        agdd20[p]   = (19.0 * agdd20[p]   + agdd[p])           / 20.0
    end

    # ---------------------------------------------------------------
    # 2. Rebuild natveg filter from present_patch
    # ---------------------------------------------------------------
    for p in bounds_patch
        mask_natvegp[p] = dgvs.present_patch[p]
    end

    # ---------------------------------------------------------------
    # 3. Light competition
    # ---------------------------------------------------------------
    cndv_light!(dgvs, eco, deadstemc_patch, leafcmax_patch,
                pftcon_in, patch, mask_natvegp, bounds_patch;
                bounds_gridcell = bounds_gridcell)

    # ---------------------------------------------------------------
    # 4. Establishment (uses all patches, not just natveg filter)
    # ---------------------------------------------------------------
    cndv_establishment!(dgvs, eco, prec365_col,
                         annsum_npp_patch, annsum_litfall_patch,
                         deadstemc_patch, leafcmax_patch,
                         pftcon_in, patch, lun, bounds_patch;
                         bounds_gridcell = bounds_gridcell)

    # ---------------------------------------------------------------
    # 5. Reset annual variables
    # ---------------------------------------------------------------
    for p in bounds_patch
        leafcmax_patch[p]  = 0.0
        t_mo_min_patch[p]  = 1.0e36
    end

    return nothing
end

# =========================================================================
#  dynCNDV — dynCNDVMod.F90
# =========================================================================

"""
    dyn_cndv_init!(dgvs, patch, bounds_patch)

Initialize FPC grid from patch weights at startup.
Sets fpcgrid = wtcol, fpcgridold = wtcol for all patches.

Ported from subroutine dynCNDV_init in dynCNDVMod.F90.
"""
function dyn_cndv_init!(dgvs::DGVSData,
                         patch::PatchData,
                         bounds_patch::UnitRange{Int})
    for p in bounds_patch
        dgvs.fpcgrid_patch[p]    = patch.wtcol[p]
        dgvs.fpcgridold_patch[p] = patch.wtcol[p]
    end
    return nothing
end

"""
    dyn_cndv_interp!(dgvs, patch, lun, bounds_patch;
                      wt1, is_beg_curr_year)

Time-interpolate CNDV PFT weights from annual to time step.

# Arguments
- `wt1::Float64`: interpolation weight for old values (= 1 - yearfrac)
- `is_beg_curr_year::Bool`: true on Jan 1 first time step → update fpcgridold

Ported from subroutine dynCNDV_interp in dynCNDVMod.F90.
"""
function dyn_cndv_interp!(dgvs::DGVSData,
                           patch::PatchData,
                           lun::LandunitData,
                           bounds_patch::UnitRange{Int};
                           wt1::Real = 1.0,
                           is_beg_curr_year::Bool = false)
    fpcgrid    = dgvs.fpcgrid_patch
    fpcgridold = dgvs.fpcgridold_patch

    for p in bounds_patch
        l = patch.landunit[p]

        if lun.itype[l] == ISTSOIL && lun.wtgcell[l] > 0.0
            patch.wtcol[p] = fpcgrid[p] + wt1 * (fpcgridold[p] - fpcgrid[p])

            if is_beg_curr_year
                fpcgridold[p] = fpcgrid[p]
            end
        end
    end

    return nothing
end

# =========================================================================
#  Accumulator updates — CNDVType::UpdateAccVars
# =========================================================================

"""
    cndv_update_acc_vars!(dgvs, eco, t_a10_patch, t_ref2m_patch,
                           bounds_patch, dtime;
                           month, day, secs)

Update AGDDTW and AGDD accumulated fields every time step.
Reset annually at Jan 1.

# Arguments
- `t_a10_patch`: 10-day running mean of 2m temperature [K]
- `t_ref2m_patch`: 2m height surface air temperature [K]
- `dtime`: time step size [seconds]
- `month`, `day`, `secs`: current date components (for annual reset)

Ported from subroutine UpdateAccVars in CNDVType.F90.
"""
function cndv_update_acc_vars!(dgvs::DGVSData,
                                eco::DGVEcophysCon,
                                t_a10_patch::Vector{<:Real},
                                t_ref2m_patch::Vector{<:Real},
                                bounds_patch::UnitRange{Int},
                                dtime::Real;
                                month::Int = 1,
                                day::Int = 1,
                                secs::Int = 0)

    # Annual reset at first time step of Jan 1
    is_reset = (month == 1 && day == 1 && secs == Int(dtime))

    if is_reset
        for p in bounds_patch
            dgvs.agddtw_patch[p] = 0.0
            dgvs.agdd_patch[p]   = 0.0
        end
    end

    # ndllf_dcd_brl_tree = 3 (Fortran 0-based), Julia PFT index = 3+1 = 4
    twmax_brl = eco.twmax[ndllf_dcd_brl_tree + 1]

    for p in bounds_patch
        # AGDDTW: growing degree days above twmax of boreal deciduous tree
        dgvs.agddtw_patch[p] += max(0.0,
            (t_a10_patch[p] - TFRZ - twmax_brl) * dtime / SECSPDAY)

        # AGDD: growing degree days above 5 C
        dgvs.agdd_patch[p] += max(0.0,
            (t_ref2m_patch[p] - (TFRZ + 5.0)) * dtime / SECSPDAY)
    end

    return nothing
end
