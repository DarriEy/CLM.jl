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
# Device-view bundle of the per-PFT params cndv_light! reads. PftconType is a
# 148-field non-parametric struct (not device-movable) and DGVEcophysCon is
# mutable (not isbits), so the needed fields are extracted + moved to the backend
# (copyto!+similar) rather than passing those structs to the kernel.
Base.@kwdef struct _CndvLightPar{V,VB}
    dwood::V; slatop::V; dsladlai::V; woody::V
    is_tree::VB; is_shrub::VB
    crownarea_max::V; reinickerp::V; allom1::V
end
Adapt.@adapt_structure _CndvLightPar

# Pass 1: per-patch crownarea/LAI/FPC update + patch->gridcell scatter of the
# tree/shrub/grass FPC totals (KA-CPU ascending order → byte-identical to host).
@kernel function _cndv_light_pass1_kernel!(crownarea, fpcgrid, fpc_inc,
        @Const(nind), @Const(deadstemc), @Const(leafcmax),
        numtrees, fpc_tree_total, fpc_inc_tree, fpc_shrub_total, fpc_grass_total,
        par, @Const(mask), @Const(p_itype), @Const(p_gridcell),
        taperT, rpiT, pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax && mask[p]
        T = eltype(crownarea)
        g = p_gridcell[p]; ivp = p_itype[p] + 1
        woody = par.woody; is_tree = par.is_tree; is_shrub = par.is_shrub
        dwood = par.dwood; slatop = par.slatop; dsladlai = par.dsladlai
        crownarea_max = par.crownarea_max; reinickerp_a = par.reinickerp; allom1_a = par.allom1

        if woody[ivp] == one(T)
            if fpcgrid[p] > zero(T) && nind[p] > zero(T)
                stocking = nind[p] / fpcgrid[p]
                stemdiam = (T(24.0) * deadstemc[p] /
                           (rpiT * stocking * dwood[ivp] * taperT))^(one(T) / T(3.0))
            else
                stemdiam = zero(T)
            end
            crownarea[p] = min(crownarea_max[ivp], allom1_a[ivp] * stemdiam^reinickerp_a[ivp])
        end

        if crownarea[p] > zero(T) && nind[p] > zero(T)
            lm_ind = leafcmax[p] * fpcgrid[p] / nind[p]
            if dsladlai[ivp] > zero(T)
                lai_ind = max(T(0.001),
                    (exp(lm_ind * dsladlai[ivp] + log(slatop[ivp])) -
                     slatop[ivp]) / dsladlai[ivp] / crownarea[p])
            else
                lai_ind = lm_ind * slatop[ivp] / crownarea[p]
            end
        else
            lai_ind = zero(T)
        end

        fpc_ind = one(T) - exp(-T(0.5) * lai_ind)
        fpcgrid_old = fpcgrid[p]
        fpcgrid[p] = crownarea[p] * nind[p] * fpc_ind
        fpc_inc[p] = max(zero(T), fpcgrid[p] - fpcgrid_old)

        if woody[ivp] == one(T)
            if is_tree[ivp]
                _scatter_add!(numtrees, g, Int32(1))   # Int32: Metal has no 64-bit atomics
                _scatter_add!(fpc_tree_total, g, fpcgrid[p])
                _scatter_add!(fpc_inc_tree, g, fpc_inc[p])
            elseif is_shrub[ivp]
                _scatter_add!(fpc_shrub_total, g, fpcgrid[p])
            end
        else
            _scatter_add!(fpc_grass_total, g, fpcgrid[p])
        end
    end
end

# Per-gridcell grass/shrub max-allowed FPC.
@kernel function _cndv_light_max_kernel!(fpc_grass_max, fpc_shrub_max,
        @Const(fpc_tree_total), @Const(fpc_grass_total), fpc_tree_maxT, gmin::Int, gmax::Int)
    g = @index(Global)
    @inbounds if gmin <= g <= gmax
        T = eltype(fpc_grass_max)
        fpc_grass_max[g] = one(T) - min(fpc_tree_total[g], fpc_tree_maxT)
        fpc_shrub_max[g] = max(zero(T), fpc_grass_max[g] - fpc_grass_total[g])
    end
end

# Pass 2: per-patch light competition (reads the gridcell totals from pass 1).
@kernel function _cndv_light_pass2_kernel!(nind, fpcgrid,
        @Const(fpc_inc), @Const(fpc_tree_total), @Const(fpc_inc_tree), @Const(numtrees),
        @Const(fpc_grass_total), @Const(fpc_grass_max),
        @Const(fpc_shrub_total), @Const(fpc_shrub_max),
        par, @Const(mask), @Const(p_itype), @Const(p_gridcell),
        fpc_tree_maxT, pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax && mask[p]
        T = eltype(nind)
        g = p_gridcell[p]; ivp = p_itype[p] + 1
        woody = par.woody; is_tree = par.is_tree; is_shrub = par.is_shrub

        if woody[ivp] == one(T) && is_tree[ivp]
            if fpc_tree_total[g] > fpc_tree_maxT
                if fpc_inc_tree[g] > zero(T)
                    excess = (fpc_tree_total[g] - fpc_tree_maxT) * fpc_inc[p] / fpc_inc_tree[g]
                else
                    excess = (fpc_tree_total[g] - fpc_tree_maxT) / T(numtrees[g])
                end
                if fpcgrid[p] > zero(T)
                    nind_kill = nind[p] * excess / fpcgrid[p]
                    nind[p] = max(zero(T), nind[p] - nind_kill)
                    fpcgrid[p] = max(zero(T), fpcgrid[p] - excess)
                else
                    nind[p] = zero(T); fpcgrid[p] = zero(T)
                end
            end
        elseif woody[ivp] == zero(T)
            if fpc_grass_total[g] > fpc_grass_max[g]
                excess = (fpc_grass_total[g] - fpc_grass_max[g]) * fpcgrid[p] / fpc_grass_total[g]
                fpcgrid[p] = max(zero(T), fpcgrid[p] - excess)
            end
        elseif woody[ivp] == one(T) && is_shrub[ivp]
            if fpc_shrub_total[g] > fpc_shrub_max[g]
                excess = one(T) - fpc_shrub_max[g] / fpc_shrub_total[g]
                if fpcgrid[p] > zero(T)
                    nind_kill = nind[p] * excess / fpcgrid[p]
                    nind[p] = max(zero(T), nind[p] - nind_kill)
                    fpcgrid[p] = max(zero(T), fpcgrid[p] - excess)
                else
                    nind[p] = zero(T); fpcgrid[p] = zero(T)
                end
            end
        end
    end
end

function cndv_light!(dgvs::DGVSData,
                     eco::DGVEcophysCon,
                     deadstemc_patch::AbstractVector{<:Real},
                     leafcmax_patch::AbstractVector{<:Real},
                     pftcon_in::PftconType,
                     patch::PatchData,
                     mask_natvegp::AbstractVector{Bool},
                     bounds_patch::UnitRange{Int};
                     bounds_gridcell::UnitRange{Int} = 1:1)

    FT = eltype(dgvs.nind_patch)
    fpc_tree_maxT = FT(0.95)
    taperT = FT(200.0)
    rpiT   = FT(RPI)

    crownarea = dgvs.crownarea_patch
    nind      = dgvs.nind_patch
    fpcgrid   = dgvs.fpcgrid_patch

    isempty(bounds_patch) && return nothing
    ngl = last(bounds_gridcell)
    npl = last(bounds_patch)

    # device-resident gridcell accumulators + per-patch fpc_inc
    z(n) = fill!(similar(dgvs.nind_patch, FT, n), zero(FT))
    fpc_tree_total = z(ngl); fpc_inc_tree = z(ngl)
    fpc_grass_total = z(ngl); fpc_shrub_total = z(ngl)
    numtrees = fill!(similar(dgvs.nind_patch, Int32, ngl), Int32(0))
    fpc_inc = z(npl)
    fpc_grass_max = z(ngl); fpc_shrub_max = z(ngl)

    # extract the per-PFT params (PftconType/eco are not kernel-passable) and
    # place them on the backend of the dgvs state. _mdf down-converts the Float64
    # params to the working precision FT (Metal: no Float64); _mdb keeps Bool.
    _mdf(h) = copyto!(similar(dgvs.nind_patch, FT, size(h)...), FT.(h))
    _mdb(h) = copyto!(similar(dgvs.nind_patch, eltype(h), size(h)...), h)
    par = _CndvLightPar(;
        dwood = _mdf(pftcon_in.dwood), slatop = _mdf(pftcon_in.slatop),
        dsladlai = _mdf(pftcon_in.dsladlai), woody = _mdf(pftcon_in.woody),
        is_tree = _mdb(pftcon_in.is_tree), is_shrub = _mdb(pftcon_in.is_shrub),
        crownarea_max = _mdf(eco.crownarea_max), reinickerp = _mdf(eco.reinickerp),
        allom1 = _mdf(eco.allom1))

    backend = _kernel_backend(dgvs.nind_patch)
    _cndv_light_pass1_kernel!(backend)(crownarea, fpcgrid, fpc_inc,
        nind, deadstemc_patch, leafcmax_patch,
        numtrees, fpc_tree_total, fpc_inc_tree, fpc_shrub_total, fpc_grass_total,
        par, mask_natvegp, patch.itype, patch.gridcell,
        taperT, rpiT, first(bounds_patch), last(bounds_patch); ndrange = length(fpcgrid))
    KA.synchronize(backend)

    _cndv_light_max_kernel!(backend)(fpc_grass_max, fpc_shrub_max,
        fpc_tree_total, fpc_grass_total, fpc_tree_maxT,
        first(bounds_gridcell), last(bounds_gridcell); ndrange = ngl)
    KA.synchronize(backend)

    _cndv_light_pass2_kernel!(backend)(nind, fpcgrid,
        fpc_inc, fpc_tree_total, fpc_inc_tree, numtrees, fpc_grass_total, fpc_grass_max,
        fpc_shrub_total, fpc_shrub_max,
        par, mask_natvegp, patch.itype, patch.gridcell,
        fpc_tree_maxT, first(bounds_patch), last(bounds_patch); ndrange = length(fpcgrid))
    KA.synchronize(backend)

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
@kernel function _cndv_dyninit_kernel!(fpcgrid, fpcgridold, @Const(wtcol), pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax
        fpcgrid[p]    = wtcol[p]
        fpcgridold[p] = wtcol[p]
    end
end

function dyn_cndv_init!(dgvs::DGVSData,
                         patch::PatchData,
                         bounds_patch::UnitRange{Int})
    isempty(bounds_patch) && return nothing
    _launch!(_cndv_dyninit_kernel!, dgvs.fpcgrid_patch, dgvs.fpcgridold_patch,
        patch.wtcol, first(bounds_patch), last(bounds_patch);
        ndrange = length(dgvs.fpcgrid_patch))
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
@kernel function _cndv_interp_kernel!(wtcol, fpcgridold, @Const(fpcgrid),
        @Const(p_landunit), @Const(lun_itype), @Const(lun_wtgcell),
        wt1T, is_beg::Bool, pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax
        T = eltype(wtcol)
        l = p_landunit[p]
        if lun_itype[l] == ISTSOIL && lun_wtgcell[l] > zero(T)
            wtcol[p] = fpcgrid[p] + wt1T * (fpcgridold[p] - fpcgrid[p])
            if is_beg
                fpcgridold[p] = fpcgrid[p]
            end
        end
    end
end

function dyn_cndv_interp!(dgvs::DGVSData,
                           patch::PatchData,
                           lun::LandunitData,
                           bounds_patch::UnitRange{Int};
                           wt1::Real = 1.0,
                           is_beg_curr_year::Bool = false)
    isempty(bounds_patch) && return nothing
    FT = eltype(dgvs.fpcgrid_patch)
    _launch!(_cndv_interp_kernel!, patch.wtcol, dgvs.fpcgridold_patch, dgvs.fpcgrid_patch,
        patch.landunit, lun.itype, lun.wtgcell,
        FT(wt1), is_beg_curr_year, first(bounds_patch), last(bounds_patch);
        ndrange = length(patch.wtcol))
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
@kernel function _cndv_accreset_kernel!(agddtw, agdd, pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax
        T = eltype(agddtw)
        agddtw[p] = zero(T); agdd[p] = zero(T)
    end
end

@kernel function _cndv_accupdate_kernel!(agddtw, agdd, @Const(t_a10), @Const(t_ref2m),
        twmax_brlT, tfrzT, dtimeT, secspdayT, pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax
        T = eltype(agddtw)
        agddtw[p] = agddtw[p] + max(zero(T), (t_a10[p] - tfrzT - twmax_brlT) * dtimeT / secspdayT)
        agdd[p]   = agdd[p]   + max(zero(T), (t_ref2m[p] - (tfrzT + T(5.0))) * dtimeT / secspdayT)
    end
end

function cndv_update_acc_vars!(dgvs::DGVSData,
                                eco::DGVEcophysCon,
                                t_a10_patch::AbstractVector{<:Real},
                                t_ref2m_patch::AbstractVector{<:Real},
                                bounds_patch::UnitRange{Int},
                                dtime::Real;
                                month::Int = 1,
                                day::Int = 1,
                                secs::Int = 0)

    isempty(bounds_patch) && return nothing
    FT = eltype(dgvs.agdd_patch)

    # Annual reset at first time step of Jan 1
    is_reset = (month == 1 && day == 1 && secs == Int(dtime))

    backend = _kernel_backend(dgvs.agdd_patch)
    if is_reset
        _cndv_accreset_kernel!(backend)(dgvs.agddtw_patch, dgvs.agdd_patch,
            first(bounds_patch), last(bounds_patch); ndrange = length(dgvs.agdd_patch))
        KA.synchronize(backend)
    end

    # ndllf_dcd_brl_tree = 3 (Fortran 0-based), Julia PFT index = 3+1 = 4
    twmax_brl = eco.twmax[ndllf_dcd_brl_tree + 1]
    _cndv_accupdate_kernel!(backend)(dgvs.agddtw_patch, dgvs.agdd_patch,
        t_a10_patch, t_ref2m_patch, FT(twmax_brl), FT(TFRZ), FT(dtime), FT(SECSPDAY),
        first(bounds_patch), last(bounds_patch); ndrange = length(dgvs.agdd_patch))
    KA.synchronize(backend)

    return nothing
end
