# FatesTwoStreamUtilsMod.jl
# Julia port of FATES src/fates/radiation/FatesTwoStreamUtilsMod.F90
#
# Glue between the FATES canopy (patch/cohort canopy-layer x PFT structure) and
# the two-stream radiation solver (TwoStreamMLPEMod). These routines:
#   * build the two-stream scattering-element inputs from the patch canopy
#     profiles (FatesConstructRadElements),
#   * map absorbed radiation back onto cohorts/leaf-layers (FatesGetCohortAbsRad),
#   * compute the patch-average sunlit/shaded LAI (FatesPatchFSun),
#   * verify the per-patch radiation balance (CheckPatchRadiationBalance),
#   * transfer the FATES PFT optical parameters into rad_params (TransferRadParams).
#
# Deps: TwoStreamMLPEMod (twostream_type + AllocInitTwoStream!/DeallocTwoStream!/
# GetNSCel!/CanopyPrep!/ZenithPrep!/GetAbsRad/GetRdUp/RadParamPrep/AllocateRadParams,
# rad_params, air_ft, area_err_thresh, rel_err_thresh), fates_patch_type/
# fates_cohort_type, ed_site_type, EDPftvarcon_inst, FatesRadiationMemMod
# (num_swb, ivis, inir), FatesConstantsMod (nearzero), EDParamsMod (nclmax),
# FatesInterfaceTypesMod (numpft), FatesAllometryMod (VegAreaLayer).
#
# Standalone — NOT added to CLMInstances or any dual-copied struct.

# Local debug flag (mirrors Fortran `debug = .false.`)
const fates_twostreamutils_debug = false

# ===============================================================================================

"""
    FatesConstructRadElements!(site, fcansno_pa, coszen_pa)

Describe the two-stream scattering elements from cohort and patch data for every
patch at `site`. Each cohort becomes its own element ("column") in its canopy
layer; an extra air element is appended to any layer that is not fully occupied.
Element areas (and the LAI/SAI mass within them) are squeezed if a layer is
over-full. The per-patch `twostr` object is (re)allocated as needed and the
geometry- and zenith-dependent coefficients are pre-processed.

`fcansno_pa` and `coszen_pa` are per-patch (1-based ifp index) snow fraction and
cosine of the zenith angle. Also (re)allocates the site scratch arrays
`taulambda_2str`, `omega_2str`, `ipiv_2str` for the largest patch.
"""
function FatesConstructRadElements!(site::ed_site_type,
                                    fcansno_pa::AbstractVector{<:Real},
                                    coszen_pa::AbstractVector{<:Real})

    # DO NOT MAKE CANOPY_OPEN_FRAC > 0 UNTIL LAI COMPRESSION HAS BEEN THOUGHT
    # THROUGH (we can't decrease area without conserving total leaf+stem area).
    canopy_open_frac = 0.0

    n_col = zeros(Int, nclmax)  # Number of parallel column elements per layer

    max_elements = -1
    ifp = 0
    patch = site.oldest_patch
    while patch !== nothing
        ifp += 1
        twostr = patch.twostr

        # --------------------------------------------------------------------
        # Identify how many elements we need
        # --------------------------------------------------------------------
        fill!(n_col, 0)
        cohort = patch.tallest
        while cohort !== nothing
            ican = cohort.canopy_layer
            n_col[ican] += 1
            cohort = cohort.shorter
        end

        # Fraction of each layer occupied by cohorts. If the layer is not fully
        # occupied, an air element is needed for the non-occupied space.
        canopy_frac = zeros(Float64, 5)
        if patch.total_canopy_area > nearzero
            cohort = patch.tallest
            while cohort !== nothing
                ican = cohort.canopy_layer
                canopy_frac[ican] += cohort.c_area / patch.total_canopy_area
                cohort = cohort.shorter
            end
        end

        for ican in 1:patch.ncl_p
            if (1.0 - canopy_frac[ican]) > area_err_thresh
                n_col[ican] += 1
            end
        end

        # --------------------------------------------------------------------
        # Handle memory: (re)allocate the two-stream object if needed
        # --------------------------------------------------------------------
        maxcol = 0
        for ican in 1:patch.ncl_p
            if n_col[ican] > maxcol
                maxcol = n_col[ican]
            end
        end

        if isempty(twostr.scelg)
            AllocInitTwoStream!(twostr, [ivis, inir], patch.ncl_p, maxcol + 2)
        else
            if size(twostr.scelg, 2) < maxcol ||
               size(twostr.scelg, 2) > (maxcol + 4) ||
               size(twostr.scelg, 1) < patch.ncl_p

                DeallocTwoStream!(twostr)
                # A little more space than necessary to avoid re-allocating
                AllocInitTwoStream!(twostr, [ivis, inir], patch.ncl_p, maxcol + 2)
            end
        end

        # --------------------------------------------------------------------
        # Fill the elements with their basic data and reference cohort→element
        # --------------------------------------------------------------------
        fill!(n_col, 0)
        cohort = patch.tallest
        while cohort !== nothing
            ft = cohort.pft
            ican = cohort.canopy_layer

            patch.canopy_mask[ican, ft] = 1

            # Every cohort gets its own element right now
            n_col[ican] += 1

            # Passing layer index 0 returns the total plant LAIs and SAIs
            (vai_top, vai_bot, elai_cohort, esai_cohort, _, _) =
                VegAreaLayer(cohort.treelai, cohort.treesai, cohort.height,
                             0, cohort.nv, cohort.pft, site.snow_depth)

            scelg = twostr.scelg[ican, n_col[ican]]

            # If this layer is fully covered by snow, treat it as an air layer
            if (elai_cohort + esai_cohort) > nearzero
                scelg.pft = ft
            else
                scelg.pft = air_ft
            end

            scelg.area = cohort.c_area / patch.total_canopy_area
            scelg.lai  = elai_cohort
            scelg.sai  = esai_cohort

            # Cohort needs to know which column it is in
            cohort.twostr_col = n_col[ican]

            cohort = cohort.shorter
        end

        for ican in 1:patch.ncl_p

            # If the canopy is not full, add an air element
            if (1.0 - canopy_frac[ican]) > area_err_thresh
                n_col[ican] += 1
                scelg = twostr.scelg[ican, n_col[ican]]
                scelg.pft  = air_ft
                scelg.area = 1.0 - canopy_frac[ican]
                scelg.lai  = 0.0
                scelg.sai  = 0.0
            end

            # If the layer is overfull, remove some area from the first element
            # that is 10x larger than the threshold
            if (canopy_frac[ican] - 1.0) > area_err_thresh
                for icol in 1:n_col[ican]
                    scelg = twostr.scelg[ican, icol]
                    if scelg.area > 10.0 * (canopy_frac[ican] - 1.0)
                        area_ratio = (scelg.area + (1.0 - canopy_frac[ican])) / scelg.area
                        scelg.area = scelg.area * area_ratio
                        scelg.lai  = scelg.lai / area_ratio
                        scelg.sai  = scelg.sai / area_ratio
                        canopy_frac[ican] = 1.0
                        break
                    end
                end
            end
        end

        for ican in 1:patch.ncl_p
            twostr.n_col[ican] = n_col[ican]
        end

        # --------------------------------------------------------------------
        # Set up some non-element parameters
        # --------------------------------------------------------------------
        twostr.n_lyr = patch.ncl_p   # Number of layers

        GetNSCel!(twostr)            # Total number of elements

        max_elements = max(max_elements, twostr.n_scel)

        # Signal that geometry-dependent coefficients need updating
        twostr.force_prep = true

        CanopyPrep!(twostr, fcansno_pa[ifp])
        ZenithPrep!(twostr, coszen_pa[ifp])

        patch = patch.younger
    end

    # ------------------------------------------------------------------------
    # Re-evaluate the scratch space for solving two-stream radiation. The
    # scratch needs to be 2x the number of computational elements for the patch
    # with the most elements.
    # ------------------------------------------------------------------------
    allocate_scratch = false
    if !isempty(site.taulambda_2str) && max_elements > 0
        n_scr = length(site.taulambda_2str)
        if 2 * max_elements > n_scr
            allocate_scratch = true
        elseif 2 * max_elements < (n_scr - 24)
            allocate_scratch = true
        end
    else
        allocate_scratch = true
    end

    if allocate_scratch && max_elements > 0
        # Twice as many spaces as elements, plus extra to prevent reallocating
        n_scr = 2 * max_elements + 8
        site.taulambda_2str = Vector{Float64}(undef, n_scr)
        site.omega_2str = Matrix{Float64}(undef, n_scr, n_scr)
        site.ipiv_2str = Vector{Int}(undef, n_scr)
    end

    return nothing
end

# =============================================================================================

"""
    FatesPatchFSun(patch) -> (fsun, laisun, laisha)

Compute the patch-average sunlit fraction `fsun` and the patch-average sunlit /
shaded LAI by summing the area-weighted leaf sunlit fractions over all elements.
"""
function FatesPatchFSun(patch::fates_patch_type)
    laisun = 0.0
    laisha = 0.0

    twostr = patch.twostr

    for ican in 1:twostr.n_lyr
        for icol in 1:twostr.n_col[ican]
            scelg = twostr.scelg[ican, icol]

            (_, _, _, _, _, _, leaf_sun_frac) =
                GetAbsRad(twostr, ican, icol, ivis, 0.0, scelg.lai + scelg.sai)

            laisun += scelg.area * scelg.lai * leaf_sun_frac
            laisha += scelg.area * scelg.lai * (1.0 - leaf_sun_frac)
        end
    end

    if (laisun + laisha) > nearzero
        fsun = laisun / (laisun + laisha)
    else
        fsun = 0.5  # Nominal value, shouldn't affect results if no leaves/light
    end

    return (fsun, laisun, laisha)
end

# ============================================================================================

"""
    CheckPatchRadiationBalance(patch, snow_depth, ib, fabd, fabi)

Loop through the cohorts in the patch, sum the absorbed radiation per cohort/leaf-
layer, and compare it to the fraction the solver calculated. Ends the run if the
balance does not close within tolerance.
"""
function CheckPatchRadiationBalance(patch::fates_patch_type, snow_depth::Real, ib::Integer,
                                    fabd::Real, fabi::Real)
    twostr = patch.twostr

    check_fab = 0.0

    cohort = patch.tallest
    while cohort !== nothing

        cohort_vaitop = zeros(Float64, 50)
        cohort_vaibot = zeros(Float64, 50)
        cohort_layer_elai = zeros(Float64, 50)
        cohort_layer_esai = zeros(Float64, 50)

        for iv in 1:cohort.nv
            (vt, vb, el, es, _, _) =
                VegAreaLayer(cohort.treelai, cohort.treesai, cohort.height,
                             iv, cohort.nv, cohort.pft, snow_depth)
            cohort_vaitop[iv] = vt
            cohort_vaibot[iv] = vb
            cohort_layer_elai[iv] = el
            cohort_layer_esai[iv] = es
        end

        cohort_elai = sum(@view cohort_layer_elai[1:cohort.nv])
        cohort_esai = sum(@view cohort_layer_esai[1:cohort.nv])

        for iv in 1:cohort.nv
            (rb_abs, rd_abs, _, _, _) =
                FatesGetCohortAbsRad(patch, cohort, ib, cohort_vaitop[iv], cohort_vaibot[iv],
                                     cohort_elai, cohort_esai)

            check_fab += (rb_abs + rd_abs) * cohort.c_area / patch.total_canopy_area
        end

        cohort = cohort.shorter
    end

    in_fab = fabd * twostr.band[ib].Rbeam_atm + fabi * twostr.band[ib].Rdiff_atm

    if abs(check_fab - in_fab) > in_fab * 10.0 * rel_err_thresh
        fates_endrun("Absorbed radiation didnt balance after cohort sum")
    end

    return nothing
end

# =============================================================================================

"""
    FatesGetCohortAbsRad(patch, cohort, ib, vaitop, vaibot, cohort_elai, cohort_esai)
        -> (rb_abs, rd_abs, rb_abs_leaf, rd_abs_leaf, leaf_sun_frac)

Retrieve the absorbed radiation on leaves and stems, plus the leaf sunlit
fraction, over the VAI interval [vaitop, vaibot] of the cohort. VAI is exposed
leaf + stem area index. The cohort's VAI coordinate is mapped to the element's
own coordinate via the element-VAI / cohort-VAI ratio, the absorbed radiation is
obtained from the solver, then mapped back and disentangled into leaf vs total.
"""
function FatesGetCohortAbsRad(patch::fates_patch_type, cohort::fates_cohort_type, ib::Integer,
                              vaitop::Real, vaibot::Real, cohort_elai::Real, cohort_esai::Real)
    twostr = patch.twostr
    scelg = twostr.scelg[cohort.canopy_layer, cohort.twostr_col]
    scelb = twostr.band[ib].scelb[cohort.canopy_layer, cohort.twostr_col]

    if (cohort_elai + cohort_esai) < nearzero
        return (0.0, 0.0, 0.0, 0.0, 0.0)
    end

    evai_cvai = (scelg.lai + scelg.sai) / (cohort_elai + cohort_esai)  # element VAI / cohort VAI

    # Convert the vai coordinate from the cohort to the element
    vai_top_el = vaitop * evai_cvai
    vai_bot_el = vaibot * evai_cvai

    # Return absorbed radiation for the element over that band
    (Rb_abs_el, Rd_abs_el, rd_abs_leaf_el, rb_abs_leaf_el, r_abs_stem_el, r_abs_snow_el, leaf_sun_frac) =
        GetAbsRad(twostr, cohort.canopy_layer, cohort.twostr_col, ib, vai_top_el, vai_bot_el)

    # rd_abs_el and rb_abs_el both contain absorption by water; abs_leaf terms do not
    rd_abs = Rd_abs_el / evai_cvai
    rb_abs = Rb_abs_el / evai_cvai

    diff_wt_leaf = (1.0 - twostr.frac_snow) * cohort_elai *
                   (1.0 - rad_params.om_leaf[ib, cohort.pft]) * rad_params.kd_leaf[cohort.pft]
    diff_wt_elem = (cohort_elai + cohort_esai) * (1.0 - scelb.om) * scelg.Kd

    beam_wt_leaf = (1.0 - twostr.frac_snow) * cohort_elai *
                   (1.0 - rad_params.om_leaf[ib, cohort.pft]) * scelg.Kb_leaf
    beam_wt_elem = (cohort_elai + cohort_esai) * (1.0 - scelb.om) * scelg.Kb

    rd_abs_leaf = rd_abs * min(1.0, diff_wt_leaf / diff_wt_elem)
    rb_abs_leaf = rb_abs * min(1.0, beam_wt_leaf / beam_wt_elem)

    return (rb_abs, rd_abs, rb_abs_leaf, rd_abs_leaf, leaf_sun_frac)
end

# =============================================================================================

"""
    TransferRadParams()

Transfer the FATES PFT optical parameters from `EDPftvarcon_inst` into the
module-level `rad_params` of the two-stream solver and derive the dependent
terms. Note: EDPftvarcon stores rhol/rhos/taul/taus as (pft, band) whereas
rad_params stores (band, pft) — the transpose is applied here.
"""
function TransferRadParams()
    edp = edpftvarcon_inst()
    np = numpft[]

    AllocateRadParams(np, num_swb)

    for ft in 1:np
        for ib in 1:num_swb
            rad_params.rhol[ib, ft] = edp.rhol[ft, ib]
            rad_params.rhos[ib, ft] = edp.rhos[ft, ib]
            rad_params.taul[ib, ft] = edp.taul[ft, ib]
            rad_params.taus[ib, ft] = edp.taus[ft, ib]
        end
        rad_params.xl[ft] = edp.xl[ft]
        rad_params.clumping_index[ft] = edp.clumping_index[ft]
    end

    RadParamPrep()

    return nothing
end
