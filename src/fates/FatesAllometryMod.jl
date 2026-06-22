# FatesAllometryMod.jl
# Julia port of FATES src/fates/biogeochem/FatesAllometryMod.F90 (3511 lines).
#
# A library of functions that calculate plant allometry and their derivatives.
# Most relationships are related to diameter [cm]; all derivatives w.r.t. change
# in diameter share those units.
#
# Allometric variables:
#   h      : height [m]
#   bagw   : above-ground woody biomass [kgC]
#   blmax  : maximum (pre-trim) leaf biomass [kgC]
#   bbgw   : below-ground woody biomass [kgC]
#   bfrmax : fine-root biomass on allometry [kgC]
#   bsap   : sapwood (above + below) biomass [kgC]
#   bdead  : structural (heartwood) biomass [kgC]
#
# Mode dispatch:  Fortran `select case` on a per-PFT allometry-mode integer ->
# Julia `if/elseif`. EVERY mode's value formula and analytic derivative are
# preserved verbatim (they feed Newton solves + downstream growth).
#
# Functions that returned a value + an `intent(out)` derivative -> Julia returns
# a tuple `(value, derivative)`. The lowest-level kernels (d2h_*, dh2bagw_*,
# d2blmax_*, etc.) always return both; the public wrappers do the same.
#
# Parameters are read from the merged `prt_params` singleton (indexed by PFT
# directly, e.g. `prt_params.allom_d2h1[ipft]`), `param_derived()`, and
# `ed_params()`. `GetCrownReduction(crowndamage)` returns the crown-loss
# fraction. `QuadraticRootsNSWC` / `QuadraticRootsSridharachary` live in the
# merged FatesUtilsMod (not used by this module, but available). `nlevleaf` is a
# const; `dinc_vai` / `dlower_vai` live on `ed_params()`.
#
# fates_r8 -> Float64. Standalone — nothing added to CLMInstances.

# =====================================================================================
# Parameter Checks / Integrated allometry consistency
# =====================================================================================

"""
    CheckIntegratedAllometries(dbh, ipft, crowndamage, canopy_trim,
        elongf_leaf, elongf_fnrt, elongf_stem, l2fr,
        bl, bfr, bsap, bstore, bdead,
        grow_leaf, grow_fr, grow_sap, grow_store, grow_dead, max_err) -> l_pass::Bool

Check the error on the carbon-allocation integration step: the integrated
quantities (`bl`, `bfr`, `bsap`, `bstore`, `bdead`) should closely match the
diagnosed (on-allometry) quantities at `dbh`. Returns `l_pass = true` if every
active pool is within `max_err`, else `false`. Mirrors Fortran
`CheckIntegratedAllometries` (returns the `l_pass` flag).
"""
function CheckIntegratedAllometries(dbh::Real, ipft::Integer, crowndamage::Integer,
        canopy_trim::Real, elongf_leaf::Real, elongf_fnrt::Real, elongf_stem::Real,
        l2fr::Real, bl::Real, bfr::Real, bsap::Real, bstore::Real, bdead::Real,
        grow_leaf::Bool, grow_fr::Bool, grow_sap::Bool, grow_store::Bool,
        grow_dead::Bool, max_err::Real)

    l_pass = true  # Default assumption is that step passed

    if grow_leaf
        bl_diag, _ = bleaf(dbh, ipft, crowndamage, canopy_trim, elongf_leaf)
        if abs(bl_diag - bl) > max_err
            l_pass = false
        end
    end

    if grow_fr
        bfr_diag, _ = bfineroot(dbh, ipft, canopy_trim, l2fr, elongf_fnrt)
        if abs(bfr_diag - bfr) > max_err
            l_pass = false
        end
    end

    if grow_sap
        _, bsap_diag, _ = bsap_allom(dbh, ipft, crowndamage, canopy_trim, elongf_stem)
        if abs(bsap_diag - bsap) > max_err
            l_pass = false
        end
    end

    if grow_store
        bstore_diag, _ = bstore_allom(dbh, ipft, crowndamage, canopy_trim)
        if abs(bstore_diag - bstore) > max_err
            l_pass = false
        end
    end

    if grow_dead
        _, bsap_diag, _ = bsap_allom(dbh, ipft, crowndamage, canopy_trim, elongf_stem)
        bagw_diag, _    = bagw_allom(dbh, ipft, crowndamage, elongf_stem)
        bbgw_diag, _    = bbgw_allom(dbh, ipft, elongf_stem)
        bdead_diag, _   = bdead_allom(bagw_diag, bbgw_diag, bsap_diag, ipft)
        if abs(bdead_diag - bdead) > max_err
            l_pass = false
        end
    end

    return l_pass
end

# =====================================================================================
# Generic height-to-diameter interface
# =====================================================================================

"""
    h2d_allom(h, ipft) -> (d, dddh)

Generic height-to-diameter wrapper. Dispatches on `prt_params.allom_hmode[ipft]`.
Returns diameter `d` [cm] and `dddh` [cm/m] (change in diameter per height).
"""
function h2d_allom(h::Real, ipft::Integer)
    p1 = prt_params.allom_d2h1[ipft]
    p2 = prt_params.allom_d2h2[ipft]
    p3 = prt_params.allom_d2h3[ipft]
    allom_hmode = prt_params.allom_hmode[ipft]

    if allom_hmode == 1       # O'Brien et al 1995, BCI
        return h2d_obrien(h, p1, p2)
    elseif allom_hmode == 2   # Poorter 2006
        return h2d_poorter2006(h, p1, p2, p3)
    elseif allom_hmode == 3   # 2 parameter power function
        return h2d_2pwr(h, p1, p2)
    elseif allom_hmode == 4   # Chave 2014
        return h2d_chave2014(h, p1, p2, p3)
    elseif allom_hmode == 5   # Martinez-Cano
        return h2d_martcano(h, p1, p2, p3)
    else
        fates_endrun("An undefined h2d allometry was specified: $(allom_hmode)")
    end
end

# =====================================================================================
# Generic diameter-to-height interface
# =====================================================================================

"""
    h_allom(d, ipft) -> (h, dhdd)

Generic diameter-to-height wrapper. Dispatches on `prt_params.allom_hmode[ipft]`.
Returns height `h` [m] and `dhdd` [m/cm].
"""
function h_allom(d::Real, ipft::Integer)
    dbh_maxh    = prt_params.allom_dbh_maxheight[ipft]
    p1          = prt_params.allom_d2h1[ipft]
    p2          = prt_params.allom_d2h2[ipft]
    p3          = prt_params.allom_d2h3[ipft]
    allom_hmode = prt_params.allom_hmode[ipft]

    if allom_hmode == 1       # obrien
        return d2h_obrien(d, p1, p2, dbh_maxh)
    elseif allom_hmode == 2   # poorter06
        return d2h_poorter2006(d, p1, p2, p3, dbh_maxh)
    elseif allom_hmode == 3   # 2 parameter power function h=a*d^b
        return d2h_2pwr(d, p1, p2, dbh_maxh)
    elseif allom_hmode == 4   # chave14
        return d2h_chave2014(d, p1, p2, p3, dbh_maxh)
    elseif allom_hmode == 5   # Martinez-Cano
        return d2h_martcano(d, p1, p2, p3, dbh_maxh)
    else
        fates_endrun("An undefined height allometry was specified: $(allom_hmode)")
    end
end

# =====================================================================================
# Generic AGB interface
# =====================================================================================

"""
    bagw_allom(d, ipft, crowndamage, elongf_stem) -> (bagw, dbagwdd)

Generic above-ground woody biomass wrapper [kgC]. Dispatches on
`prt_params.allom_amode[ipft]`, then applies crown-damage and stem-phenology
(`elongf_stem`) reductions. Returns `(bagw, dbagwdd)`.
"""
function bagw_allom(d::Real, ipft::Integer, crowndamage::Integer, elongf_stem::Real)
    p1           = prt_params.allom_agb1[ipft]
    p2           = prt_params.allom_agb2[ipft]
    p3           = prt_params.allom_agb3[ipft]
    p4           = prt_params.allom_agb4[ipft]
    wood_density = prt_params.wood_density[ipft]
    c2b          = prt_params.c2b[ipft]
    agb_frac     = prt_params.allom_agb_frac[ipft]
    allom_amode  = prt_params.allom_amode[ipft]

    branch_frac = param_derived().branch_frac[ipft]

    local bagw::Float64, dbagwdd::Float64
    if allom_amode == 1       # salda
        h, dhdd = h_allom(d, ipft)
        bagw, dbagwdd = dh2bagw_salda(d, h, dhdd, p1, p2, p3, p4, wood_density, c2b, agb_frac)
    elseif allom_amode == 2   # 2par_pwr (woodland dbh->drc switch handled upstream)
        bagw, dbagwdd = d2bagw_2pwr(d, p1, p2, c2b)
    elseif allom_amode == 3   # chave14
        h, dhdd = h_allom(d, ipft)
        bagw, dbagwdd = dh2bagw_chave2014(d, h, dhdd, p1, p2, wood_density, c2b)
    elseif allom_amode == 4   # 3par_pwr
        h, dhdd = h_allom(d, ipft)
        bagw, dbagwdd = dh2bagw_3pwr(d, h, dhdd, p1, p2, p3, wood_density, c2b)
    elseif allom_amode == 5   # 3par_pwr_grass
        h, dhdd = h_allom(d, ipft)
        bagw, dbagwdd = dh2bagw_3pwr_grass(d, h, dhdd, p1, p2, p3, c2b)
    else
        fates_endrun("An undefined AGB allometry was specified: $(allom_amode)")
    end

    # Potentially reduce AGB based on crown damage and/or phenology (elongf_stem).
    if crowndamage > 1
        crown_reduction = GetCrownReduction(crowndamage)
        bagw    = elongf_stem * (bagw - (bagw * branch_frac * crown_reduction))
        dbagwdd = elongf_stem * (dbagwdd - (dbagwdd * branch_frac * crown_reduction))
    else
        bagw    = elongf_stem * bagw
        dbagwdd = elongf_stem * dbagwdd
    end

    return bagw, dbagwdd
end

# =====================================================================================
# Generic diameter-to-maximum-leaf-biomass interface
# =====================================================================================

"""
    blmax_allom(d, ipft) -> (blmax, dblmaxdd)

Generic maximum (pre-trim, "on allometry") leaf biomass wrapper [kgC]. Dispatches
on `prt_params.allom_lmode[ipft]`. Returns `(blmax, dblmaxdd)`.
"""
function blmax_allom(d::Real, ipft::Integer)
    dbh_maxh    = prt_params.allom_dbh_maxheight[ipft]
    rho         = prt_params.wood_density[ipft]
    slatop      = prt_params.slatop[ipft]
    c2b         = prt_params.c2b[ipft]
    allom_lmode = prt_params.allom_lmode[ipft]
    p1          = prt_params.allom_d2bl1[ipft]
    p2          = prt_params.allom_d2bl2[ipft]
    p3          = prt_params.allom_d2bl3[ipft]

    if allom_lmode == 1       # salda
        return d2blmax_salda(d, p1, p2, p3, rho, dbh_maxh, c2b)
    elseif allom_lmode == 2   # 2par_pwr
        return d2blmax_2pwr(d, p1, p2, c2b)
    elseif allom_lmode == 3   # dh2blmax_2pwr
        return dh2blmax_2pwr(d, p1, p2, dbh_maxh, c2b)
    elseif allom_lmode == 4   # dh2blmax_3pwr
        h, dhdd = h_allom(d, ipft)
        return dh2blmax_3pwr(d, h, dhdd, p1, p2, p3, slatop, dbh_maxh, c2b)
    elseif allom_lmode == 5   # dh2blmax_3pwr_grass
        h, dhdd = h_allom(d, ipft)
        return dh2blmax_3pwr_grass(d, h, dhdd, p1, p2, p3, dbh_maxh, c2b)
    else
        fates_endrun("An undefined leaf allometry was specified: $(allom_lmode)")
    end
end

# =====================================================================================
# Generic crown area allometry wrapper
# =====================================================================================

"""
    carea_allom(dbh, nplant, site_spread, ipft, crowndamage; inverse=false) -> (dbh, c_area)

Generic crown area wrapper. Dispatches on `prt_params.allom_lmode[ipft]`. When
`inverse=false`, computes per-cohort `c_area` [m2] from `dbh`; when `inverse=true`,
computes `dbh` from `c_area`. Some modes cap the diameter at `allom_dbh_maxheight`;
in the inverse-capped case where the cap is hit, `dbh` is returned as
`fates_unset_r8` (signalling the caller that dbh<->area are no longer
proportional). Because `dbh` is `intent(inout)` in Fortran, this returns the
(possibly updated) `(dbh, c_area)` tuple.
"""
function carea_allom(dbh::Real, nplant::Real, site_spread::Real, ipft::Integer,
                     crowndamage::Integer; inverse::Bool=false)
    dbh_maxh    = prt_params.allom_dbh_maxheight[ipft]
    allom_lmode = prt_params.allom_lmode[ipft]
    d2bl_p2     = prt_params.allom_d2bl2[ipft]
    d2bl_ediff  = prt_params.allom_blca_expnt_diff[ipft]
    d2ca_min    = prt_params.allom_d2ca_coefficient_min[ipft]
    d2ca_max    = prt_params.allom_d2ca_coefficient_max[ipft]

    dbh    = Float64(dbh)
    c_area = 0.0          # only meaningful as input in the inverse case

    do_inverse = inverse
    # NOTE: in the inverse branch the Fortran divides the supplied per-cohort
    # c_area by nplant; the caller must therefore pass c_area through the
    # returned value chain. We expose c_area as a kw below for that use.
    return _carea_allom_impl(dbh, nplant, site_spread, ipft, crowndamage,
                             dbh_maxh, allom_lmode, d2bl_p2, d2bl_ediff,
                             d2ca_min, d2ca_max, do_inverse, c_area)
end

"""
    carea_allom(dbh, nplant, site_spread, ipft, crowndamage, c_area; inverse=false) -> (dbh, c_area)

Crown-area wrapper variant that also takes an input `c_area` (required for the
`inverse=true` path, where `dbh` is solved from `c_area`).
"""
function carea_allom(dbh::Real, nplant::Real, site_spread::Real, ipft::Integer,
                     crowndamage::Integer, c_area::Real; inverse::Bool=false)
    dbh_maxh    = prt_params.allom_dbh_maxheight[ipft]
    allom_lmode = prt_params.allom_lmode[ipft]
    d2bl_p2     = prt_params.allom_d2bl2[ipft]
    d2bl_ediff  = prt_params.allom_blca_expnt_diff[ipft]
    d2ca_min    = prt_params.allom_d2ca_coefficient_min[ipft]
    d2ca_max    = prt_params.allom_d2ca_coefficient_max[ipft]

    return _carea_allom_impl(Float64(dbh), nplant, site_spread, ipft, crowndamage,
                             dbh_maxh, allom_lmode, d2bl_p2, d2bl_ediff,
                             d2ca_min, d2ca_max, inverse, Float64(c_area))
end

function _carea_allom_impl(dbh::Float64, nplant::Real, site_spread::Real,
                           ipft::Integer, crowndamage::Integer, dbh_maxh::Real,
                           allom_lmode::Integer, d2bl_p2::Real, d2bl_ediff::Real,
                           d2ca_min::Real, d2ca_max::Real, do_inverse::Bool,
                           c_area::Float64)
    if do_inverse
        c_area = c_area / nplant
    end

    dbh_eff = dbh
    local capped_allom::Bool
    if allom_lmode == 1
        dbh_eff = min(dbh, dbh_maxh)
        dbh_eff, c_area = carea_2pwr(dbh_eff, site_spread, d2bl_p2, d2bl_ediff,
                                     d2ca_min, d2ca_max, crowndamage, c_area, do_inverse)
        capped_allom = true
    elseif allom_lmode == 2   # 2par_pwr
        dbh, c_area = carea_2pwr(dbh, site_spread, d2bl_p2, d2bl_ediff,
                                 d2ca_min, d2ca_max, crowndamage, c_area, do_inverse)
        capped_allom = false
    elseif allom_lmode == 3 || allom_lmode == 5
        dbh_eff = min(dbh, dbh_maxh)
        dbh_eff, c_area = carea_2pwr(dbh_eff, site_spread, d2bl_p2, d2bl_ediff,
                                     d2ca_min, d2ca_max, crowndamage, c_area, do_inverse)
        capped_allom = true
    elseif allom_lmode == 4
        dbh_eff = min(dbh, dbh_maxh)
        height, _ = h_allom(dbh, ipft)
        dbh_eff, height, c_area = carea_3pwr(dbh_eff, height, ipft, dbh_maxh,
                                             site_spread, d2bl_p2, d2bl_ediff,
                                             d2ca_min, d2ca_max, crowndamage,
                                             c_area, do_inverse)
        capped_allom = true
    else
        fates_endrun("An undefined leaf allometry was specified: $(allom_lmode)")
    end

    if capped_allom && do_inverse
        if dbh_eff < dbh_maxh
            dbh = dbh_eff
        else
            # We have hit the area cap; dbh<->area are no longer proportional.
            dbh = fates_unset_r8
        end
    end

    c_area = c_area * nplant

    return dbh, c_area
end

# =====================================================================================

"""
    bleaf(d, ipft, crowndamage, canopy_trim, elongf_leaf) -> (bl, dbldd)

Actual target leaf biomass [kgC] = `blmax * canopy_trim`, after crown-damage and
leaf-phenology (`elongf_leaf`) reductions. Returns `(bl, dbldd)`.
"""
function bleaf(d::Real, ipft::Integer, crowndamage::Integer, canopy_trim::Real,
               elongf_leaf::Real)
    blmax, dblmaxdd = blmax_allom(d, ipft)

    bl    = blmax * canopy_trim
    dbldd = dblmaxdd * canopy_trim

    if crowndamage > 1
        crown_reduction = GetCrownReduction(crowndamage)
        bl    = elongf_leaf * bl * (1.0 - crown_reduction)
        # NOTE: Fortran derivative here uses dblmaxdd*canopy_trim (== the
        # pre-damage dbldd), preserved verbatim.
        dbldd = elongf_leaf * dblmaxdd * canopy_trim * (1.0 - crown_reduction)
    else
        bl    = elongf_leaf * bl
        dbldd = elongf_leaf * dbldd
    end

    return bl, dbldd
end

# =====================================================================================

"""
    storage_fraction_of_target(c_store_target, c_store) -> frac

Ratio of the storage pool to its target, floored at 0 (can exceed 1). The target
is floored at `nearzero` to avoid division by zero.
"""
function storage_fraction_of_target(c_store_target::Real, c_store::Real)
    return max(0.0, c_store / max(c_store_target, nearzero))
end

# =====================================================================================

"""
    tree_lai(leaf_c, pft, c_area, nplant, cl, canopy_lai, vcmax25top) -> tree_lai

Individual-tree LAI from total leaf area and canopy area, accounting for an
exponential SLA profile through the canopy (with a `slamax` cap). `canopy_lai` is
the vector of total LAI per canopy layer; `cl` is this tree's canopy layer index.
"""
function tree_lai(leaf_c::Real, pft::Integer, c_area::Real, nplant::Real,
                  cl::Integer, canopy_lai::AbstractVector, vcmax25top::Real)

    if leaf_c < -1.1 * calloc_abs_error || pft == 0
        fates_endrun("negative leaf carbon in LAI calculation, or pft was zero: leaf_c=$(leaf_c), pft=$(pft)")
    end

    slat = g_per_kg * prt_params.slatop[pft]   # m2/g -> m2/kg
    leafc_per_unitarea = leaf_c / (c_area / nplant)  # KgC/m2

    if leafc_per_unitarea > 0.0
        if cl == 1
            canopy_lai_above = 0.0
        else
            canopy_lai_above = sum(@view canopy_lai[1:cl-1])
        end

        kn = decay_coeff_vcmax(vcmax25top,
                               prt_params.leafn_vert_scaler_coeff1[pft],
                               prt_params.leafn_vert_scaler_coeff2[pft])

        sla_max = g_per_kg * prt_params.slamax[pft]   # m2/gC -> m2/kgC
        leafc_slamax = (slat - sla_max * exp(-1.0 * kn * canopy_lai_above)) /
                       (-1.0 * kn * slat * sla_max)
        if leafc_slamax < 0.0
            leafc_slamax = 0.0
        end

        if leafc_per_unitarea <= leafc_slamax
            tlai = (log(exp(-1.0 * kn * canopy_lai_above) -
                        kn * slat * leafc_per_unitarea) +
                    (kn * canopy_lai_above)) / (-1.0 * kn)

            clim = (exp(-1.0 * kn * canopy_lai_above)) / (kn * slat)
            if leafc_per_unitarea >= clim
                fates_endrun("too much leafc_per_unitarea: $(leafc_per_unitarea) $(clim) $(pft) $(canopy_lai_above)")
            end
        else  # leafc_per_unitarea > leafc_slamax
            tlai = ((log(exp(-1.0 * kn * canopy_lai_above) -
                         kn * slat * leafc_slamax) +
                     (kn * canopy_lai_above)) / (-1.0 * kn)) +
                   (leafc_per_unitarea - leafc_slamax) * sla_max

            clim = (exp(-1.0 * kn * canopy_lai_above)) / (kn * slat)
            if leafc_slamax >= clim
                fates_endrun("too much leafc_slamax: $(leafc_per_unitarea) $(leafc_slamax) $(clim) $(pft) $(canopy_lai_above)")
            end
        end
    else
        tlai = 0.0
    end

    return tlai
end

# =====================================================================================

"""
    tree_sai(pft, dbh, crowndamage, canopy_trim, elongf_stem, c_area, nplant,
             cl, canopy_lai, treelai, vcmax25top, call_id) -> tree_sai

Individual-tree SAI = `elongf_stem * allom_sai_scaler * target_lai`, where
`target_lai` is computed from fully-flushed target leaf biomass. Aborts (matching
Fortran) if `treelai + tree_sai` exceeds the total vai-bin capacity `sum(dinc_vai)`.
"""
function tree_sai(pft::Integer, dbh::Real, crowndamage::Integer, canopy_trim::Real,
                  elongf_stem::Real, c_area::Real, nplant::Real, cl::Integer,
                  canopy_lai::AbstractVector, treelai::Real, vcmax25top::Real,
                  call_id::Integer)

    # Assume fully flushed leaves (elongf_leaf = 1) so SAI is leaf-phenology
    # independent; SAI can still be downscaled by stem phenology.
    target_bleaf, _ = bleaf(dbh, pft, crowndamage, canopy_trim, 1.0)

    target_lai = tree_lai(target_bleaf, pft, c_area, nplant, cl, canopy_lai, vcmax25top)

    sai = elongf_stem * prt_params.allom_sai_scaler[pft] * target_lai

    if (treelai + sai) > sum(ed_params().dinc_vai)
        fates_endrun("tree leaf+stem maxed out array size: lai=$(treelai) sai=$(sai) pft=$(pft) call_id=$(call_id) dbh=$(dbh)")
    end

    return sai
end

# =====================================================================================

"""
    leafc_from_treelai(treelai, pft, c_area, nplant, cl, vcmax25top) -> leafc

Inverse of [`tree_lai`](@ref): the leaf carbon [kg] needed to generate a given
`treelai`. Only valid in the top canopy layer (`cl == 1`).
"""
function leafc_from_treelai(treelai::Real, pft::Integer, c_area::Real, nplant::Real,
                            cl::Integer, vcmax25top::Real)
    if treelai < 0.0 || pft == 0
        fates_endrun("negative tree lai in leafc_from_treelai, or pft zero: treelai=$(treelai) pft=$(pft)")
    end
    if cl > 1
        fates_endrun("in sub-canopy layer in leafc_from_treelai (cl=$(cl), pft=$(pft)); not set up for lower canopy layers.")
    end

    slat    = g_per_kg * prt_params.slatop[pft]
    sla_max = g_per_kg * prt_params.slamax[pft]

    kn = decay_coeff_vcmax(vcmax25top,
                           prt_params.leafn_vert_scaler_coeff1[pft],
                           prt_params.leafn_vert_scaler_coeff2[pft])

    if treelai > 0.0
        leafc_slamax = max(0.0, (slat - sla_max) / (-1.0 * kn * slat * sla_max))

        tree_lai_at_slamax = (log(1.0 - kn * slat * leafc_slamax)) / (-1.0 * kn)

        if treelai < tree_lai_at_slamax
            leafc_per_unitarea = (1.0 - exp(treelai * (-1.0 * kn))) / (kn * slat)
        else  # we exceed the maximum sla
            leafc_linear_phase = (treelai - tree_lai_at_slamax) / sla_max
            leafc_per_unitarea = leafc_slamax + leafc_linear_phase
        end
        leafc = leafc_per_unitarea * (c_area / nplant)
    else
        leafc = 0.0
    end

    return leafc
end

# =====================================================================================
# Generic sapwood biomass interface
# =====================================================================================

"""
    bsap_allom(d, ipft, crowndamage, canopy_trim, elongf_stem) -> (sapw_area, bsap, dbsapdd)

Generic sapwood wrapper. Dispatches on `prt_params.allom_smode[ipft]`. Returns the
sapwood cross-section area `sapw_area` [m2], biomass `bsap` [kgC] and `dbsapdd`.
Mode 1 caps sapwood at 95% of total woody biomass.
"""
function bsap_allom(d::Real, ipft::Integer, crowndamage::Integer, canopy_trim::Real,
                    elongf_stem::Real)
    max_frac = 0.95   # cap on aboveground sapwood as fraction of woody/fibrous tissue

    agb_frac    = prt_params.allom_agb_frac[ipft]
    branch_frac = param_derived().branch_frac[ipft]
    allom_smode = prt_params.allom_smode[ipft]

    if allom_smode == 1
        # Linearly related to leaf area based on target leaf biomass and slatop.
        # Fully flushed leaves (elongf_leaf=1) -> sapwood independent of leaf phenology.
        h, dhdd = h_allom(d, ipft)
        bl, dbldd = bleaf(d, ipft, 1, canopy_trim, 1.0)
        sapw_area, bsap, dbsapdd = bsap_ltarg_slatop(d, h, dhdd, bl, dbldd, ipft)

        if crowndamage > 1
            crown_reduction = GetCrownReduction(crowndamage)
            bsap    = elongf_stem * (bsap - (bsap * agb_frac * branch_frac * crown_reduction))
            dbsapdd = elongf_stem * (dbsapdd - (dbsapdd * agb_frac * branch_frac * crown_reduction))
        else
            bsap    = elongf_stem * bsap
            dbsapdd = elongf_stem * dbsapdd
        end

        # Cap total woody biomass.
        bagw, dbagwdd = bagw_allom(d, ipft, crowndamage, elongf_stem)
        bbgw, dbbgwdd = bbgw_allom(d, ipft, elongf_stem)

        bsap_cap = max_frac * (bagw + bbgw)
        if bsap > bsap_cap
            bsap    = bsap_cap
            dbsapdd = max_frac * (dbagwdd + dbbgwdd)
        end

        return sapw_area, bsap, dbsapdd

    elseif allom_smode == 2
        # Grass-only: bsap = bagw + bbgw, no dead woody biomass, no crown damage.
        bagw, dbagwdd = bagw_allom(d, ipft, crowndamage, elongf_stem)
        bbgw, dbbgwdd = bbgw_allom(d, ipft, elongf_stem)
        bsap = bagw + bbgw

        bsap = elongf_stem * bsap
        # NOTE: Fortran sets dbsapdd = elongf_stem*(...) then immediately
        # overwrites it with the un-scaled (dbagwdd+dbbgwdd). Preserved verbatim:
        # the final dbsapdd does NOT include elongf_stem.
        dbsapdd = dbagwdd + dbbgwdd

        # sapw_area is not set for grass (dummy); return 0.0 to match the
        # uninitialized Fortran intent(out) being unused by grass callers.
        sapw_area = 0.0
        return sapw_area, bsap, dbsapdd
    else
        fates_endrun("An undefined sapwood allometry was specified: $(allom_smode)")
    end
end

# =====================================================================================
# Generic below-ground woody biomass interface
# =====================================================================================

"""
    bbgw_allom(d, ipft, elongf_stem) -> (bbgw, dbbgwdd)

Generic below-ground woody biomass wrapper [kgC]. Dispatches on
`prt_params.allom_cmode[ipft]`. Mode 1 = constant proportionality with `bagw`.
"""
function bbgw_allom(d::Real, ipft::Integer, elongf_stem::Real)
    allom_cmode = prt_params.allom_cmode[ipft]
    if allom_cmode == 1   # constant
        # bbgw not affected by damage -> use crowndamage=1; stem phenology enters
        # through bagw (which is already downscaled by elongf_stem).
        bagw, dbagwdd = bagw_allom(d, ipft, 1, elongf_stem)
        return bbgw_const(d, bagw, dbagwdd, ipft)
    else
        fates_endrun("An undefined coarse root allometry was specified: $(allom_cmode)")
    end
end

# =====================================================================================
# Fine root biomass allometry wrapper
# =====================================================================================

"""
    bfineroot(d, ipft, canopy_trim, l2fr, elongf_fnrt) -> (bfr, dbfrdd)

Generic actual target fine-root biomass [kgC]. Dispatches on
`prt_params.allom_fmode[ipft]` (1 = proportionality with trimmed bleaf, 2 = with
untrimmed bleaf), then applies the fine-root phenology factor `elongf_fnrt`.
"""
function bfineroot(d::Real, ipft::Integer, canopy_trim::Real, l2fr::Real,
                   elongf_fnrt::Real)
    allom_fmode = prt_params.allom_fmode[ipft]

    local bfr::Float64, dbfrdd::Float64
    if allom_fmode == 1   # constant proportionality with TRIMMED target bleaf
        blmax, dblmaxdd = blmax_allom(d, ipft)
        bfr    = blmax * l2fr * canopy_trim
        dbfrdd = dblmaxdd * l2fr * canopy_trim
    elseif allom_fmode == 2  # constant proportionality with UNTRIMMED target bleaf
        blmax, dblmaxdd = blmax_allom(d, ipft)
        bfr    = blmax * l2fr
        dbfrdd = dblmaxdd * l2fr
    else
        fates_endrun("An undefined fine root allometry was specified: $(allom_fmode)")
    end

    # Reduce fine-root biomass due to phenology.
    bfr    = elongf_fnrt * bfr
    dbfrdd = elongf_fnrt * dbfrdd

    return bfr, dbfrdd
end

# =====================================================================================
# Storage biomass interface
# =====================================================================================

"""
    bstore_allom(d, ipft, crowndamage, canopy_trim) -> (bstore, dbstoredd)

Generic allometric target storage carbon [kgC]. Dispatches on
`prt_params.allom_stmode[ipft]` (1 = cushion * trimmed bleaf, 2 = cushion *
untrimmed blmax).
"""
function bstore_allom(d::Real, ipft::Integer, crowndamage::Integer, canopy_trim::Real)
    allom_stmode = prt_params.allom_stmode[ipft]
    cushion      = prt_params.cushion[ipft]

    if allom_stmode == 1
        bl, dbldd = bleaf(d, ipft, crowndamage, canopy_trim, 1.0)
        return bstore_blcushion(d, bl, dbldd, cushion, ipft)
    elseif allom_stmode == 2
        blmax, dblmaxdd = blmax_allom(d, ipft)
        return bstore_blcushion(d, blmax, dblmaxdd, cushion, ipft)
    else
        fates_endrun("An undefined storage allometry was specified: $(allom_stmode)")
    end
end

# =====================================================================================
# Dead biomass interface
# =====================================================================================

"""
    bdead_allom(bagw, bbgw, bsap, ipft; dbagwdd, dbbgwdd, dbsapdd) -> (bdead, dbdeaddd)

Structural (heartwood) biomass [kgC] diagnosed by mass balance. Dispatches on
`prt_params.allom_amode[ipft]` (mode 1: `bagw/agb_frac`; modes 2-5:
`bagw + bbgw - bsap`). The derivative `dbdeaddd` is returned only if the relevant
input derivatives are supplied (else `NaN`, mirroring the Fortran absent-optional
behaviour — callers that need it pass the inputs).
"""
function bdead_allom(bagw::Real, bbgw::Real, bsap::Real, ipft::Integer;
                     dbagwdd::Union{Real,Nothing}=nothing,
                     dbbgwdd::Union{Real,Nothing}=nothing,
                     dbsapdd::Union{Real,Nothing}=nothing)
    agb_fraction = prt_params.allom_agb_frac[ipft]
    allom_amode  = prt_params.allom_amode[ipft]

    if allom_amode == 1
        # Saldarriaga: assume proportionality between bdead and bagw.
        bdead = bagw / agb_fraction
        dbdeaddd = (dbagwdd !== nothing) ? (dbagwdd / agb_fraction) : NaN
        return bdead, dbdeaddd
    elseif allom_amode in (2, 3, 4, 5)
        bdead = bagw + bbgw - bsap
        if dbagwdd !== nothing && dbbgwdd !== nothing && dbsapdd !== nothing
            dbdeaddd = dbagwdd + dbbgwdd - dbsapdd
        else
            dbdeaddd = NaN
        end
        return bdead, dbdeaddd
    else
        fates_endrun("An undefined AGB allometry was specified: $(allom_amode)")
    end
end

# =====================================================================================
# Specific bbgw relationships
# =====================================================================================

"""
    bbgw_const(d, bagw, dbagwdd, ipft) -> (bbgw, dbbgwdd)

Constant-proportionality below-ground woody biomass: `(1/agb_frac - 1) * bagw`.
"""
function bbgw_const(d::Real, bagw::Real, dbagwdd::Real, ipft::Integer)
    agb_fraction = prt_params.allom_agb_frac[ipft]
    bbgw    = (1.0 / agb_fraction - 1.0) * bagw
    dbbgwdd = (1.0 / agb_fraction - 1.0) * dbagwdd
    return bbgw, dbbgwdd
end

# =====================================================================================
# Specific d2bsap relationships
# =====================================================================================

"""
    bsap_ltarg_slatop(d, h, dhdd, bleaf, dbleafdd, ipft) -> (sapw_area, bsap, dbsapdd)

Sapwood carbon from its leaf-area-per-sapwood-area proportionality with the
plant's target leaf area (Calvo-Alvarado / Christofferson). Above + below ground.
Returns sapwood cross-section area [m2], biomass [kgC] and derivative.
"""
function bsap_ltarg_slatop(d::Real, h::Real, dhdd::Real, bleaf_in::Real,
                           dbleafdd::Real, ipft::Integer)
    la_per_sa_int = prt_params.allom_la_per_sa_int[ipft]
    la_per_sa_slp = prt_params.allom_la_per_sa_slp[ipft]
    slatop        = prt_params.slatop[ipft]
    wood_density  = prt_params.wood_density[ipft]
    c2b           = prt_params.c2b[ipft]
    agb_fraction  = prt_params.allom_agb_frac[ipft]

    # Sapwood biomass per linear height per kgC of leaf [m-1].
    hbl2bsap = slatop * g_per_kg * wood_density * kg_per_Megag / (c2b * cm2_per_m2)

    # Sapwood cross-section area [m2] (no c2b: wood_density already in biomass).
    la_per_sa = la_per_sa_int + h * la_per_sa_slp
    sapw_area = slatop * bleaf_in * g_per_kg / (la_per_sa * cm2_per_m2)

    # term1 combines two height-dependent functions.
    term1 = h / (la_per_sa_int + h * la_per_sa_slp)
    bsap  = (hbl2bsap / agb_fraction) * term1 * bleaf_in

    dterm1_dh = la_per_sa_int / (la_per_sa_int + la_per_sa_slp * h)^2.0
    dterm1_dd = dterm1_dh * dhdd
    dbsapdd   = hbl2bsap / agb_fraction * (bleaf_in * dterm1_dd + term1 * dbleafdd)

    return sapw_area, bsap, dbsapdd
end

# =====================================================================================
# Specific storage relationships
# =====================================================================================

"""
    bstore_blcushion(d, bl, dbldd, cushion, ipft) -> (bstore, dbstoredd)

Allometric target storage = `cushion * bl`.
"""
function bstore_blcushion(d::Real, bl::Real, dbldd::Real, cushion::Real, ipft::Integer)
    bstore    = bl * cushion
    dbstoredd = dbldd * cushion
    return bstore, dbstoredd
end

# =====================================================================================
# Specific d2blmax relationships
# =====================================================================================

"""
    d2blmax_salda(d, p1, p2, p3, rho, dbh_maxh, c2b) -> (blmax, dblmaxdd)

Saldarriaga leaf-biomass allometry, capped at `dbh_maxh`.
"""
function d2blmax_salda(d::Real, p1::Real, p2::Real, p3::Real, rho::Real,
                       dbh_maxh::Real, c2b::Real)
    if d < dbh_maxh
        blmax = p1 * d^p2 * rho^p3
    else
        blmax = p1 * dbh_maxh^p2 * rho^p3
    end

    if d < dbh_maxh
        dblmaxdd = p1 * p2 * d^(p2 - 1.0) * rho^p3
    else
        dblmaxdd = 0.0
    end

    return blmax, dblmaxdd
end

# =====================================================================================

"""
    d2blmax_2pwr(d, p1, p2, c2b) -> (blmax, dblmaxdd)

Two-parameter power-function leaf biomass: `(p1*d^p2)/c2b`.
"""
function d2blmax_2pwr(d::Real, p1::Real, p2::Real, c2b::Real)
    blmax    = (p1 * d^p2) / c2b
    dblmaxdd = p1 * p2 * d^(p2 - 1.0) / c2b
    return blmax, dblmaxdd
end

# =====================================================================================

"""
    dh2blmax_2pwr(d, p1, p2, dbh_maxh, c2b) -> (blmax, dblmaxdd)

Like [`d2blmax_2pwr`](@ref) but caps leaf biomass once the plant reaches `dbh_maxh`.
"""
function dh2blmax_2pwr(d::Real, p1::Real, p2::Real, dbh_maxh::Real, c2b::Real)
    blmax = (p1 * min(d, dbh_maxh)^p2) / c2b

    if d >= dbh_maxh
        dblmaxdd = 0.0
    else
        dblmaxdd = p1 * p2 * d^(p2 - 1.0) / c2b
    end

    return blmax, dblmaxdd
end

# =====================================================================================

"""
    dh2blmax_3pwr(d, h, dhdd, p1, p2, p3, slatop, dbh_maxh, c2b) -> (blmax, dblmaxdd)

Leaf biomass from `(d^2*h)` and SLA (Lescure/Longo-style), capped at `dbh_maxh`.
"""
function dh2blmax_3pwr(d::Real, h::Real, dhdd::Real, p1::Real, p2::Real, p3::Real,
                       slatop::Real, dbh_maxh::Real, c2b::Real)
    duse = min(d, dbh_maxh)

    blmax = p1 * (duse * duse * h)^p2 * slatop^p3 / c2b

    if d >= dbh_maxh
        dblmaxdd = 0.0
    else
        dblmaxdd = p2 * blmax * (2.0 / duse + dhdd / h)
    end

    return blmax, dblmaxdd
end

# =====================================================================================

"""
    dh2blmax_3pwr_grass(d, h, dhdd, p1, p2, p3, dbh_maxh, c2b) -> (blmax, dblmaxdd)

Grass leaf-biomass allometry (Gao et al. 2024): `p1*d^p2*h^p3/c2b`, capped at `dbh_maxh`.
"""
function dh2blmax_3pwr_grass(d::Real, h::Real, dhdd::Real, p1::Real, p2::Real,
                             p3::Real, dbh_maxh::Real, c2b::Real)
    duse = min(d, dbh_maxh)

    blmax = p1 * duse^p2 * h^p3 / c2b

    if d >= dbh_maxh
        dblmaxdd = 0.0
    else
        dblmaxdd = blmax * (p2 / duse + p3 * dhdd / h)
    end

    return blmax, dblmaxdd
end

# =====================================================================================
# Diameter to height (D2H) functions
# =====================================================================================

"""
    d2h_chave2014(d, p1, p2, p3, dbh_maxh) -> (h, dhdd)

Chave et al. 2014 height allometry (no hard height cap; `dbh_maxh` caps growth).
"""
function d2h_chave2014(d::Real, p1::Real, p2::Real, p3::Real, dbh_maxh::Real)
    p1e = p1   # eclim already removed
    if d >= dbh_maxh
        h = exp(p1e + p2 * log(dbh_maxh) + p3 * log(dbh_maxh)^2.0)
    else
        h = exp(p1e + p2 * log(d) + p3 * log(d)^2.0)
    end

    if d >= dbh_maxh
        dhdd = 0.0
    else
        dhdd = exp(p1e) * (2.0 * p3 * d^(p2 - 1.0 + p3 * log(d)) * log(d) +
                           p2 * d^(p2 - 1.0 + p3 * log(d)))
    end

    return h, dhdd
end

# =====================================================================================

"""
    d2h_poorter2006(d, p1, p2, p3, dbh_maxh) -> (h, dhdd)

Poorter et al. 2006 asymptotic (Weibull) height allometry. `p1` is the asymptote.
"""
function d2h_poorter2006(d::Real, p1::Real, p2::Real, p3::Real, dbh_maxh::Real)
    h = p1 * (1.0 - exp(p2 * min(d, dbh_maxh)^p3))

    if d >= dbh_maxh
        dhdd = 0.0
    else
        dhdd = -p1 * exp(p2 * d^p3) * p3 * p2 * d^(p3 - 1.0)
    end

    return h, dhdd
end

# =====================================================================================

"""
    d2h_2pwr(d, p1, p2, dbh_maxh) -> (h, dhdd)

Two-parameter power-function height: `h = p1*d^p2`, capped at `dbh_maxh`.
"""
function d2h_2pwr(d::Real, p1::Real, p2::Real, dbh_maxh::Real)
    h = p1 * min(d, dbh_maxh)^p2

    if d >= dbh_maxh
        dhdd = 0.0
    else
        dhdd = (p2 * p1) * d^(p2 - 1.0)
    end

    return h, dhdd
end

# =====================================================================================

"""
    d2h_obrien(d, p1, p2, dbh_maxh) -> (h, dhdd)

O'Brien et al. 1995 (BCI) log-log height allometry, capped at `dbh_maxh`.
"""
function d2h_obrien(d::Real, p1::Real, p2::Real, dbh_maxh::Real)
    h = 10.0^(log10(min(d, dbh_maxh)) * p1 + p2)

    if d >= dbh_maxh
        dhdd = 0.0
    else
        dhdd = p1 * 10.0^p2 * d^(p1 - 1.0)
    end

    return h, dhdd
end

# =====================================================================================

"""
    d2h_martcano(d, p1, p2, p3, dbh_maxh) -> (h, dhdd)

Martinez-Cano et al. 2016 three-parameter Michaelis-Menten height allometry,
capped at `dbh_maxh`.
"""
function d2h_martcano(d::Real, p1::Real, p2::Real, p3::Real, dbh_maxh::Real)
    h = (p1 * min(d, dbh_maxh)^p2) / (p3 + min(d, dbh_maxh)^p2)

    if d >= dbh_maxh
        dhdd = 0.0
    else
        dhdd = ((p2 * p1 * d^(p2 - 1.0)) * (p3 + d^p2) -
                (p2 * d^(p2 - 1.0)) * (p1 * d^p2)) /
               (p3 + d^p2)^2.0
    end

    return h, dhdd
end

# =====================================================================================
# Diameter to (above-ground) biomass
# =====================================================================================

"""
    dh2bagw_chave2014(d, h, dhdd, p1, p2, wood_density, c2b) -> (bagw, dbagwdd)

Chave et al. 2014 above-ground biomass from `(wood_density * d^2 * h)`.
"""
function dh2bagw_chave2014(d::Real, h::Real, dhdd::Real, p1::Real, p2::Real,
                           wood_density::Real, c2b::Real)
    bagw = (p1 * (wood_density * d^2.0 * h)^p2) / c2b

    dbagwdd1 = (p1 * wood_density^p2) / c2b
    dbagwdd2 = p2 * d^(2.0 * p2) * h^(p2 - 1.0) * dhdd
    dbagwdd3 = h^p2 * 2.0 * p2 * d^(2.0 * p2 - 1.0)
    dbagwdd  = dbagwdd1 * (dbagwdd2 + dbagwdd3)

    return bagw, dbagwdd
end

# =====================================================================================

"""
    dh2bagw_3pwr(d, h, dhdd, p1, p2, p3, wood_density, c2b) -> (bagw, dbagwdd)

Above-ground biomass from `(d^2*h)` with an independent wood-density exponent
(intermediate between Saldarriaga and Chave).
"""
function dh2bagw_3pwr(d::Real, h::Real, dhdd::Real, p1::Real, p2::Real, p3::Real,
                      wood_density::Real, c2b::Real)
    bagw = p1 * (d * d * h)^p2 * wood_density^p3 / c2b
    dbagwdd = p2 * bagw * (2.0 / d + dhdd / h)
    return bagw, dbagwdd
end

# =====================================================================================

"""
    dh2bagw_3pwr_grass(d, h, dhdd, p1, p2, p3, c2b) -> (bagw, dbagwdd)

Grass above-ground (non-leaf) biomass (Gao et al. 2024): `p1*d^p2*h^p3/c2b`.
"""
function dh2bagw_3pwr_grass(d::Real, h::Real, dhdd::Real, p1::Real, p2::Real,
                            p3::Real, c2b::Real)
    bagw = p1 * (d^p2) * (h^p3) / c2b
    dbagwdd = p2 * bagw / d + p3 * bagw * dhdd / h
    return bagw, dbagwdd
end

# =====================================================================================

"""
    d2bagw_2pwr(d, p1, p2, c2b) -> (bagw, dbagwdd)

Two-parameter power-function above-ground biomass: `(p1*d^p2)/c2b` (temperate/boreal).
"""
function d2bagw_2pwr(d::Real, p1::Real, p2::Real, c2b::Real)
    bagw    = (p1 * d^p2) / c2b
    dbagwdd = (p2 * p1 * d^(p2 - 1.0)) / c2b
    return bagw, dbagwdd
end

# =====================================================================================

"""
    dh2bagw_salda(d, h, dhdd, p1, p2, p3, p4, wood_density, c2b, allom_agb_frac) -> (bagw, dbagwdd)

Saldarriaga et al. 1988 stem biomass from height, dbh and wood density (the
above-ground fraction of total dead biomass).
"""
function dh2bagw_salda(d::Real, h::Real, dhdd::Real, p1::Real, p2::Real, p3::Real,
                       p4::Real, wood_density::Real, c2b::Real, allom_agb_frac::Real)
    bagw = allom_agb_frac * p1 * (h^p2) * (d^p3) * (wood_density^p4)

    term1 = allom_agb_frac * p1 * (wood_density^p4)
    term2 = (h^p2) * p3 * d^(p3 - 1.0)
    term3 = p2 * h^(p2 - 1.0) * (d^p3) * dhdd
    dbagwdd = term1 * (term2 + term3)

    return bagw, dbagwdd
end

# =====================================================================================
# Height-to-diameter conversions (effective diameter for un-capped plants)
# =====================================================================================

"""
    h2d_chave2014(h, p1, p2, p3) -> (de, ddedh)

Chave 2014 inverse: effective diameter `de` [cm] from height (via quadratic root
of the log relation) and `ddedh` [cm/m].
"""
function h2d_chave2014(h::Real, p1::Real, p2::Real, p3::Real)
    p1e   = p1   # eclim already removed
    eroot = (-p2 + sqrt(p2^2.0 + 4.0 * log(h * exp(-p1e)) * p3)) / (2.0 * p3)

    de = exp(eroot)

    dhpdd = exp(p1e) * (p3 * 2.0 * de^(p2 - 1.0) * log(de) * exp(p3 * log(de)^2) +
                        p2 * de^(p2 - 1.0) * exp(p3 * log(de)^2.0))
    ddedh = 1.0 / dhpdd

    return de, ddedh
end

# =====================================================================================

"""
    h2d_poorter2006(h, p1, p2, p3) -> (d, dddh)

Poorter 2006 inverse (sapling initialization): diameter from height and `dddh`.
"""
function h2d_poorter2006(h::Real, p1::Real, p2::Real, p3::Real)
    d = (log(1.0 - h / p1) / p2)^(1.0 / p3)
    dddh = -((log(1 - h / p1) / p2)^(1.0 / p3 - 1.0)) / (p2 * p3 * (p1 - h))
    return d, dddh
end

# =====================================================================================

"""
    h2d_2pwr(h, p1, p2) -> (d, dddh)

Two-parameter power-function inverse: `d = (h/p1)^(1/p2)`.
"""
function h2d_2pwr(h::Real, p1::Real, p2::Real)
    d = (h / p1)^(1.0 / p2)
    dddh = (1.0 / p2) * (1.0 / p1)^(1.0 / p2) * h^(1.0 / p2 - 1.0)
    return d, dddh
end

# =====================================================================================

"""
    h2d_obrien(h, p1, p2) -> (d, dddh)

O'Brien inverse: `d = 10^((log10(h)-p2)/p1)`.
"""
function h2d_obrien(h::Real, p1::Real, p2::Real)
    d = 10.0^((log10(h) - p2) / p1)
    dddh = 1.0 / (p1 * 10.0^p2 * d^(p1 - 1.0))
    return d, dddh
end

# =====================================================================================

"""
    h2d_martcano(h, p1, p2, p3) -> (d, dddh)

Martinez-Cano inverse: `d = ((h*p3)/(p1-h))^(1/p2)`.
"""
function h2d_martcano(h::Real, p1::Real, p2::Real, p3::Real)
    d = ((h * p3) / (p1 - h))^(1.0 / p2)

    dddh = (((1.0 / p2) * (h * p3)^(1.0 / p2 - 1.0)) * ((p1 - h)^(1.0 / p2)) -
            ((1.0 / p2) * (p1 - h)^(1.0 / p2 - 1.0)) * ((h * p3)^(1.0 / p2))) /
           ((p1 - h)^(1.0 / p2))^2.0

    return d, dddh
end

# =====================================================================================
# Crown depth
# =====================================================================================

"""
    CrownDepth(height, ipft) -> crown_depth

Depth of a plant's crown [m]. Dispatches on `prt_params.allom_dmode[ipft]`
(1 = linear fraction of height; 2 = Poorter-style power law, capped at height).
"""
function CrownDepth(height::Real, ipft::Integer)
    p1 = prt_params.allom_h2cd1[ipft]
    p2 = prt_params.allom_h2cd2[ipft]
    allom_dmode = prt_params.allom_dmode[ipft]

    if allom_dmode == 1   # linear with height
        crown_depth = p1 * height
    elseif allom_dmode == 2   # power law (Poorter), capped at height
        crown_depth = min(height, p1 * height^p2)
    else
        fates_endrun("Invalid crown depth mode for PFT $(ipft): allom_dmode=$(allom_dmode). Valid: 1 or 2.")
    end

    return crown_depth
end

# =====================================================================================
# Specific diameter-to-crown-area allometries
# =====================================================================================

"""
    carea_2pwr(dbh, spread, d2bl_p2, d2bl_ediff, d2ca_min, d2ca_max, crowndamage, c_area, inverse) -> (dbh, c_area)

Crown area for one plant [m2] as a power function of DBH and canopy spread. If
`inverse`, solves DBH from `c_area`. Returns the (possibly updated) `(dbh, c_area)`.
"""
function carea_2pwr(dbh::Real, spread::Real, d2bl_p2::Real, d2bl_ediff::Real,
                    d2ca_min::Real, d2ca_max::Real, crowndamage::Integer,
                    c_area::Real, inverse::Bool)
    dbh = Float64(dbh)
    c_area = Float64(c_area)

    # Default: same exponent as dbh->bleaf so canopy depth is growth-invariant,
    # adjustable via allom_blca_expnt_diff (default 0).
    crown_area_to_dbh_exponent = d2bl_p2 + d2bl_ediff

    spreadterm = spread * d2ca_max + (1.0 - spread) * d2ca_min

    if !inverse
        c_area = spreadterm * dbh^crown_area_to_dbh_exponent
        if crowndamage > 1
            crown_reduction = GetCrownReduction(crowndamage)
            c_area = c_area * (1.0 - crown_reduction)
        end
    else
        if crowndamage > 1
            crown_reduction = GetCrownReduction(crowndamage)
            c_area = c_area / (1.0 - crown_reduction)
        end
        dbh = (c_area / spreadterm)^(1.0 / crown_area_to_dbh_exponent)
    end

    return dbh, c_area
end

# =====================================================================================

"""
    carea_3pwr(dbh, height, ipft, dbh_maxh, spread, dh2bl_p2, dh2bl_ediff,
               dh2ca_min, dh2ca_max, crowndamage, c_area, inverse) -> (dbh, height, c_area)

Crown area for one plant [m2] as a power function of size `(dbh^2 * height)`. If
`inverse`, solves for DBH via [`size2dbh`](@ref). Returns `(dbh, height, c_area)`.
"""
function carea_3pwr(dbh::Real, height::Real, ipft::Integer, dbh_maxh::Real,
                    spread::Real, dh2bl_p2::Real, dh2bl_ediff::Real,
                    dh2ca_min::Real, dh2ca_max::Real, crowndamage::Integer,
                    c_area::Real, inverse::Bool)
    dbh    = Float64(dbh)
    height = Float64(height)
    c_area = Float64(c_area)

    dh2ca_p1 = spread * dh2ca_max + (1.0 - spread) * dh2ca_min
    dh2ca_p2 = dh2bl_p2 + dh2bl_ediff

    if !inverse
        sz = dbh * dbh * height
        c_area = dh2ca_p1 * sz^dh2ca_p2
        if crowndamage > 1
            crown_reduction = GetCrownReduction(crowndamage)
            c_area = c_area * (1.0 - crown_reduction)
        end
    else
        if crowndamage > 1
            crown_reduction = GetCrownReduction(crowndamage)
            c_area = c_area * (1.0 - crown_reduction)
        end
        sz = (c_area / dh2ca_p1)^(1.0 / dh2ca_p2)
        dbh = size2dbh(sz, ipft, dbh, dbh_maxh)
    end

    return dbh, height, c_area
end

# =====================================================================================
# Root profiles
# =====================================================================================

"""
    set_root_fraction(root_fraction, ft, zi; max_nlevroot=nothing) -> root_fraction

Fill `root_fraction` (modified in place and returned) with the normalized root
mass fraction per soil layer for PFT `ft`. Dispatches on
`prt_params.fnrt_prof_mode[ft]` (1 = Jackson beta, 2 = 1-param exponential,
3 = 2-param exponential). `zi` is the depth-of-layer-interface vector indexed 0:n
(pass a 0-based `OffsetArray`-like vector or a plain length n+1 vector treated as
`zi[0..n]`). A small correction is applied to the largest-fraction bin so the
profile sums to 1.
"""
function set_root_fraction(root_fraction::AbstractVector, ft::Integer,
                           zi::AbstractVector; max_nlevroot::Union{Integer,Nothing}=nothing)
    # zi is indexed 0:nlevroot. Accept either a 0-based offset vector or a plain
    # 1-based vector of length (nlevroot+1). Build a 0-based accessor.
    zi0 = _zi_zero_based(zi)
    nlev_full = length(zi0) - 1   # ubound(zi,1)

    if length(zi0) != (length(root_fraction) + 1)
        fates_endrun("layer interface array should be 1 larger than root fraction array")
    end

    # Zero everywhere (some layers may be inactive).
    fill!(root_fraction, 0.0)

    nlevroot = nlev_full
    if max_nlevroot !== nothing
        nlevroot = min(max_nlevroot, nlevroot)
    end

    jackson_beta_profile_type   = 1
    exponential_1p_profile_type = 2
    exponential_2p_profile_type = 3

    root_profile_type = round(Int, prt_params.fnrt_prof_mode[ft])

    rf_view = @view root_fraction[1:nlevroot]
    if root_profile_type == exponential_1p_profile_type
        exponential_1p_root_profile!(rf_view, zi0, nlevroot, prt_params.fnrt_prof_a[ft])
    elseif root_profile_type == jackson_beta_profile_type
        jackson_beta_root_profile!(rf_view, zi0, nlevroot, prt_params.fnrt_prof_a[ft])
    elseif root_profile_type == exponential_2p_profile_type
        exponential_2p_root_profile!(rf_view, zi0, nlevroot,
                                     prt_params.fnrt_prof_a[ft], prt_params.fnrt_prof_b[ft])
    else
        fates_endrun("An undefined root profile type was specified")
    end

    correction = 1.0 - sum(root_fraction)
    corr_id = argmax(root_fraction)
    root_fraction[corr_id] += correction

    return root_fraction
end

# Build a 0-based accessor closure-free wrapper. We represent zi as a function of
# the 0-based index by storing the underlying 1-based vector; helpers below index
# zi0[k] for k in 0:nlev meaning the (k+1)-th element.
struct _ZeroBasedVec{V<:AbstractVector}
    v::V
end
Base.getindex(z::_ZeroBasedVec, k::Integer) = z.v[k + 1]
Base.length(z::_ZeroBasedVec) = length(z.v)
_zi_zero_based(zi::AbstractVector) = _ZeroBasedVec(zi)

"""
    exponential_2p_root_profile!(root_fraction, zi0, nlevsoil, a, b)

Two-parameter exponential root profile (normalized). `zi0` is a 0-based accessor.
"""
function exponential_2p_root_profile!(root_fraction::AbstractVector, zi0, nlevsoil::Integer,
                                      a::Real, b::Real)
    sum_rootfr = 0.0
    for lev in 1:nlevsoil
        root_fraction[lev] = 0.5 * (exp(-a * zi0[lev-1]) +
                                    exp(-b * zi0[lev-1]) -
                                    exp(-a * zi0[lev]) -
                                    exp(-b * zi0[lev]))
        sum_rootfr += root_fraction[lev]
    end
    for lev in 1:nlevsoil
        root_fraction[lev] /= sum_rootfr
    end
    return root_fraction
end

"""
    exponential_1p_root_profile!(root_fraction, zi0, nlevsoil, a)

One-parameter exponential root profile (normalized).
"""
function exponential_1p_root_profile!(root_fraction::AbstractVector, zi0, nlevsoil::Integer,
                                      a::Real)
    sum_rootfr = 0.0
    for lev in 1:nlevsoil
        root_fraction[lev] = exp(-a * 0.5 * (zi0[lev] + zi0[lev-1]))
        sum_rootfr += root_fraction[lev]
    end
    for lev in 1:nlevsoil
        root_fraction[lev] /= sum_rootfr
    end
    return root_fraction
end

"""
    jackson_beta_root_profile!(root_fraction, zi0, nlevsoil, a)

Jackson et al. 1996 beta-distribution root profile (normalized).
"""
function jackson_beta_root_profile!(root_fraction::AbstractVector, zi0, nlevsoil::Integer,
                                    a::Real)
    sum_rootfr = 0.0
    for lev in 1:nlevsoil
        root_fraction[lev] = (a^(zi0[lev-1] * 100.0) - a^(zi0[lev] * 100.0))
        sum_rootfr += root_fraction[lev]
    end
    for lev in 1:nlevsoil
        root_fraction[lev] /= sum_rootfr
    end
    return root_fraction
end

# =====================================================================================

"""
    decay_coeff_vcmax(vcmax25top, slope_param, intercept_param) -> kn

Vertical canopy decay coefficient `kn = exp(slope_param*vcmax25top - intercept_param)`
(Lloyd et al. 2010), used to attenuate vcmax / SLA / leaf respiration with depth.
"""
function decay_coeff_vcmax(vcmax25top::Real, slope_param::Real, intercept_param::Real)
    return exp(slope_param * vcmax25top - intercept_param)
end

# =====================================================================================
# ForceDBH
# =====================================================================================

"""
    ForceDBH(ipft, crowndamage, canopy_trim, elongf_leaf, elongf_stem, d;
             bdead=nothing, bl=nothing) -> (d, h)

Solve for the diameter consistent with a known structural biomass (`bdead`, woody
PFTs) or leaf biomass (`bl`, grasses) via a damped Newton search, then return the
matching height. `d` is `intent(inout)` in Fortran (used as the initial guess);
returns the updated `(d, h)`.
"""
function ForceDBH(ipft::Integer, crowndamage::Integer, canopy_trim::Real,
                  elongf_leaf::Real, elongf_stem::Real, d::Real;
                  bdead::Union{Real,Nothing}=nothing,
                  bl::Union{Real,Nothing}=nothing)
    step_frac0  = 0.9
    max_counter = 200

    d = Float64(d)

    if prt_params.woody[ipft] == itrue
        if bdead === nothing
            fates_endrun("woody plants must use structure for dbh reset")
        end

        _, bt_sap, dbt_sap_dd = bsap_allom(d, ipft, crowndamage, canopy_trim, elongf_stem)
        bt_agw, dbt_agw_dd = bagw_allom(d, ipft, crowndamage, elongf_stem)
        bt_bgw, dbt_bgw_dd = bbgw_allom(d, ipft, elongf_stem)
        bt_dead, dbt_dead_dd = bdead_allom(bt_agw, bt_bgw, bt_sap, ipft;
                                           dbagwdd=dbt_agw_dd, dbbgwdd=dbt_bgw_dd,
                                           dbsapdd=dbt_sap_dd)

        counter = 0
        step_frac = step_frac0
        while (bdead - bt_dead) > calloc_abs_error && dbt_dead_dd > 0.0
            dd    = step_frac * (bdead - bt_dead) / dbt_dead_dd
            d_try = d + dd

            _, bt_sap, dbt_sap_dd = bsap_allom(d_try, ipft, crowndamage, canopy_trim, elongf_stem)
            bt_agw, dbt_agw_dd = bagw_allom(d_try, ipft, crowndamage, elongf_stem)
            bt_bgw, dbt_bgw_dd = bbgw_allom(d_try, ipft, elongf_stem)
            bt_dead_try, dbt_dead_dd_try = bdead_allom(bt_agw, bt_bgw, bt_sap, ipft;
                                                       dbagwdd=dbt_agw_dd,
                                                       dbbgwdd=dbt_bgw_dd,
                                                       dbsapdd=dbt_sap_dd)

            if bt_dead_try > (bdead + calloc_abs_error)
                step_frac = step_frac * 0.5
            else
                step_frac   = step_frac0
                d           = d_try
                bt_dead     = bt_dead_try
                dbt_dead_dd = dbt_dead_dd_try
            end
            counter += 1
            if counter > max_counter
                fates_endrun("Having trouble converging on dbh reset")
            end
        end
    else
        if bl === nothing
            fates_endrun("grasses must use leaf for dbh reset")
        end

        bt_leaf, dbt_leaf_dd = bleaf(d, ipft, crowndamage, canopy_trim, elongf_leaf)

        counter = 0
        step_frac = step_frac0
        while (bl - bt_leaf) > calloc_abs_error && dbt_leaf_dd > 0.0
            dd    = step_frac * (bl - bt_leaf) / dbt_leaf_dd
            d_try = d + dd

            # NOTE: Fortran passes elongf_stem (not elongf_leaf) here; verbatim.
            bt_leaf_try, dbt_leaf_dd_try = bleaf(d_try, ipft, crowndamage, canopy_trim, elongf_stem)

            if bt_leaf_try > (bl + calloc_abs_error)
                step_frac = step_frac * 0.5
            else
                step_frac   = step_frac0
                d           = d_try
                bt_leaf     = bt_leaf_try
                dbt_leaf_dd = dbt_leaf_dd_try
            end
            counter += 1
            if counter > max_counter
                fates_endrun("Having trouble converging on dbh reset")
            end
        end
    end

    h, _ = h_allom(d, ipft)
    return d, h
end

# =====================================================================================
# VegAreaLayer
# =====================================================================================

"""
    VegAreaLayer(tree_lai, tree_sai, tree_height, iv, nv, pft, snow_depth)
        -> (vai_top, vai_bot, elai_layer, esai_layer, tlai_layer, tsai_layer)

Exposed leaf/stem area (m2 leaf/stem per m2 ground, inside the crown) for veg
layer `iv` of `nv`. Returns the bin top/bottom VAI, exposed leaf/stem area indices,
and total leaf/stem area indices for the layer. Uses the constant-physical-depth
assumption (`layer_height_method = 1`) as in the Fortran default.
"""
function VegAreaLayer(tree_lai::Real, tree_sai::Real, tree_height::Real,
                      iv::Integer, nv::Integer, pft::Integer, snow_depth::Real)
    layer_height_const_depth = 1
    layer_height_method = layer_height_const_depth

    tree_vai = tree_lai + tree_sai

    crown_depth = CrownDepth(tree_height, pft)
    frac_crown_depth = crown_depth / tree_height

    dinc_vai   = ed_params().dinc_vai
    dlower_vai = ed_params().dlower_vai

    if tree_vai > 0.0
        if iv == 0
            vai_top = 0.0
            vai_bot = tree_vai
        else
            if iv > 1
                vai_top = dlower_vai[iv] - dinc_vai[iv]
            else
                vai_top = 0.0
            end

            if iv < nv
                vai_bot = dlower_vai[iv]
            else
                vai_bot = tree_vai
            end
        end

        if layer_height_method == layer_height_const_depth
            if iv == 0
                layer_top_height = tree_height
                layer_bot_height = tree_height * (1.0 - frac_crown_depth)
            else
                layer_top_height = tree_height * (1.0 - Float64(iv - 1) / Float64(nv) * frac_crown_depth)
                layer_bot_height = tree_height * (1.0 - Float64(iv) / Float64(nv) * frac_crown_depth)
            end
        else
            layer_top_height = tree_height * (1.0 - frac_crown_depth * vai_top / tree_vai)
            layer_bot_height = tree_height * (1.0 - frac_crown_depth * vai_bot / tree_vai)
        end

        fraction_exposed = 1.0 - max(0.0, (min(1.0, (snow_depth - layer_bot_height) /
                                                    (layer_top_height - layer_bot_height))))

        tlai = (vai_bot - vai_top) * tree_lai / tree_vai
        tsai = (vai_bot - vai_top) * tree_sai / tree_vai

        tlai_layer = tlai
        tsai_layer = tsai

        elai_layer = fraction_exposed * tlai
        esai_layer = fraction_exposed * tsai

        vai_bot = vai_top + fraction_exposed * (vai_bot - vai_top)
    else
        tlai_layer = 0.0
        tsai_layer = 0.0
        elai_layer = 0.0
        esai_layer = 0.0
        vai_bot = 0.0
        vai_top = 0.0
    end

    return vai_top, vai_bot, elai_layer, esai_layer, tlai_layer, tsai_layer
end

# =====================================================================================
# Cubic spline (utility; not name-spaced with allom_)
# =====================================================================================

"""
    cspline(x1, x2, y1, y2, dydx1, dydx2, x) -> (y, dydx)

Cubic Hermite spline interpolation between two endpoints with known slopes.
"""
function cspline(x1::Real, x2::Real, y1::Real, y2::Real, dydx1::Real, dydx2::Real, x::Real)
    t = (x - x1) / (x2 - x1)
    a = dydx1 * (x2 - x1) - (y2 - y1)
    b = -dydx2 * (x2 - x1) + (y2 - y1)

    y    = (1.0 - t) * y1 + t * y2 + t * (1.0 - t) * (a * (1.0 - t) + b * t)
    dydx = (y2 - y1) / (x2 - x1) + (1.0 - 2.0 * t) * (a * (1.0 - t) + b * t) / (x2 - x1) +
           t * (1.0 - t) * (b - a) / (x2 - x1)
    return y, dydx
end

# =====================================================================================
# size2dbh — find DBH from size (DBH^2 * Height) via Newton + Regula Falsi
# =====================================================================================

"""
    size2dbh(size, ipft, dbh, dbh_maxh) -> dbh

Solve for DBH given `size = DBH^2 * Height` (non-linear because height depends on
DBH). Uses Newton's method with a Regula-Falsi (Illinois) fallback. `dbh` is the
initial guess (`intent(inout)`); returns the converged DBH.
"""
function size2dbh(size::Real, ipft::Integer, dbh::Real, dbh_maxh::Real)
    toler      = 1.0e-12
    maxit_newt = 10
    maxit_rf   = 100

    dbh = Float64(dbh)

    # Maximum size beyond which height is constant -> analytic DBH.
    hgt, _ = h_allom(dbh_maxh, ipft)
    size_maxh = dbh_maxh * dbh_maxh * hgt
    if size >= size_maxh
        return sqrt(size / hgt)
    end

    # First guess: current DBH.
    adbh = dbh
    hgt, dhgtddbh = h_allom(adbh, ipft)
    afun  = adbh * adbh * hgt - size
    deriv = 2.0 * adbh * hgt + adbh * adbh * dhgtddbh

    zdbh = adbh
    zfun = afun

    # Newton's method loop.
    converged = false
    for _ in 1:maxit_newt
        if abs(deriv) < toler
            break
        end

        adbh = zdbh
        afun = zfun

        zdbh = adbh - afun / deriv
        hgt, dhgtddbh = h_allom(zdbh, ipft)
        zfun  = zdbh * zdbh * hgt - size
        deriv = 2.0 * zdbh * hgt + zdbh * zdbh * dhgtddbh

        converged = abs(adbh - zdbh) < toler * zdbh
        if converged
            return 0.5 * (adbh + zdbh)
        elseif abs(zfun) < nearzero
            return zdbh
        end
    end

    # Newton failed -> Regula Falsi. Need two guesses with opposite-sign evals.
    zside = false
    if afun * zfun <= -nearzero
        zside = true
    else
        if abs(zfun - afun) < 100.0 * toler * adbh
            delta = 100.0 * toler * adbh
        else
            delta = max(abs(afun * (zdbh - adbh) / (zfun - afun)), 100.0 * toler * adbh)
        end

        zdbh = adbh + delta
        zside = false
        for iti in 1:maxit_rf
            zdbh = adbh + Float64((-1)^iti * div(iti + 3, 2)) * delta
            hgt, _ = h_allom(zdbh, ipft)
            zfun  = zdbh * zdbh * hgt - size
            zside = afun * zfun < -nearzero
            if zside
                break
            end
        end

        if !zside
            fates_endrun("Second guess for Regula Falsi method not found (size=$(size), dbh=$(dbh)).")
        end
    end

    # Regula Falsi loop.
    converged = false
    rdbh = adbh
    for _ in 1:maxit_rf
        rdbh = (zfun * adbh - afun * zdbh) / (zfun - afun)

        converged = abs(rdbh - adbh) < toler * max(rdbh, adbh)
        if converged
            break
        end

        hgt, _ = h_allom(rdbh, ipft)
        rfun = rdbh * rdbh * hgt - size

        if abs(rfun) < nearzero
            converged = true
            break
        elseif rfun * afun <= -nearzero
            zdbh = rdbh
            zfun = rfun
            if zside
                afun = afun * 0.5
            end
            zside = true
        else
            adbh = rdbh
            afun = rfun
            if !zside
                zfun = zfun * 0.5
            end
            zside = false
        end
    end

    if converged
        return rdbh
    else
        fates_endrun("Size to DBH routine failed to converge! (size=$(size), dbh=$(dbh))")
    end
end
