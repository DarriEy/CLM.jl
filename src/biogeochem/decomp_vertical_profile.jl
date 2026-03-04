# ==========================================================================
# Ported from: src/soilbiogeochem/SoilBiogeochemVerticalProfileMod.F90
# Calculate vertical profiles for distributing soil and litter C and N
#
# Public functions:
#   soil_bgc_vertical_profile!  -- main vertical profile computation
#
# Module-level parameter:
#   surfprof_exp  -- steepness of surface profile (1/e-folding depth) [1/m]
# ==========================================================================

# How steep the profile is for surface components (1/e_folding depth) (1/m)
const SURFPROF_EXP = 10.0

"""
    soil_bgc_vertical_profile!(
        soilbiogeochem_state::SoilBiogeochemStateData,
        active_layer::ActiveLayerData,
        soilstate::SoilStateData,
        col::ColumnData,
        patch::PatchData;
        mask_bgc_soilc::BitVector,
        mask_bgc_vegp::BitVector,
        nlevdecomp::Int,
        nlevdecomp_full::Int,
        dzsoi_decomp::Vector{Float64},
        zsoi::Vector{Float64},
        use_bedrock::Bool=false,
        zmin_bedrock::Float64=ZMIN_BEDROCK)

Calculate vertical profiles for distributing soil and litter C and N.

This routine computes profiles of leaves, stems, fine roots, coarse roots,
N fixation, and N deposition as a function of depth, normalised over the
active layer.  Where the soil is fully frozen or rootless, all inputs are
placed in the top decomposition layer.

The FATES branch uses the surface e-folding profile for both fixation and
deposition (partially because fixation may be free-living).

Ported from `SoilBiogeochemVerticalProfile` in
`SoilBiogeochemVerticalProfileMod.F90`.

# BUG note from Fortran
Because of this routine's placement in the driver sequence (called very
early, before weights and filters are updated), it computes values over
inactive as well as active points.
"""
function soil_bgc_vertical_profile!(
        soilbiogeochem_state::SoilBiogeochemStateData,
        active_layer::ActiveLayerData,
        soilstate::SoilStateData,
        col::ColumnData,
        patch::PatchData;
        mask_bgc_soilc::BitVector,
        mask_bgc_vegp::BitVector,
        nlevdecomp::Int,
        nlevdecomp_full::Int,
        dzsoi_decomp::Vector{Float64},
        zsoi::Vector{Float64},
        use_bedrock::Bool=false,
        zmin_bedrock::Float64=ZMIN_BEDROCK)

    # --- Aliases (matching Fortran associate block) ---
    altmax_lastyear_indx = active_layer.altmax_lastyear_indx_col
    crootfr              = soilstate.crootfr_patch

    nfixation_prof = soilbiogeochem_state.nfixation_prof_col
    ndep_prof      = soilbiogeochem_state.ndep_prof_col
    leaf_prof      = soilbiogeochem_state.leaf_prof_patch
    froot_prof     = soilbiogeochem_state.froot_prof_patch
    croot_prof     = soilbiogeochem_state.croot_prof_patch
    stem_prof      = soilbiogeochem_state.stem_prof_patch

    nc = length(mask_bgc_soilc)
    np = length(mask_bgc_vegp)

    # --- Local workspace ---
    cinput_rootfr     = zeros(Float64, np, nlevdecomp_full)
    col_cinput_rootfr = zeros(Float64, nc, nlevdecomp_full)

    delta = 1.0e-10

    # ------------------------------------------------------------------
    # 1. Surface profile: shallow exponential, optionally zeroed below
    #    bedrock when use_bedrock is true.
    # ------------------------------------------------------------------
    surface_prof = zeros(Float64, nlevdecomp)
    for j in 1:nlevdecomp
        surface_prof[j] = exp(-SURFPROF_EXP * zsoi[j]) / dzsoi_decomp[j]
        if use_bedrock
            if zsoi[j] > zmin_bedrock
                surface_prof[j] = 0.0
            end
        end
    end

    # ------------------------------------------------------------------
    # 2. Initialize all output profiles to zero.
    # ------------------------------------------------------------------
    for p in eachindex(mask_bgc_vegp)
        for j in 1:nlevdecomp_full
            leaf_prof[p, j]  = 0.0
            froot_prof[p, j] = 0.0
            croot_prof[p, j] = 0.0
            stem_prof[p, j]  = 0.0
        end
    end
    for c in eachindex(mask_bgc_soilc)
        for j in 1:nlevdecomp_full
            nfixation_prof[c, j] = 0.0
            ndep_prof[c, j]      = 0.0
        end
    end

    # ------------------------------------------------------------------
    # 3. Compute cinput_rootfr (root fraction normalised by layer
    #    thickness) for each vegetated patch.
    # ------------------------------------------------------------------
    for p in eachindex(mask_bgc_vegp)
        mask_bgc_vegp[p] || continue
        if patch.itype[p] != noveg
            for j in 1:nlevdecomp
                cinput_rootfr[p, j] = crootfr[p, j] / dzsoi_decomp[j]
            end
        else
            cinput_rootfr[p, 1] = 0.0
        end
    end

    # ------------------------------------------------------------------
    # 4. Integrate profiles over the active layer for each vegetated
    #    patch and compute normalised profiles.
    # ------------------------------------------------------------------
    for p in eachindex(mask_bgc_vegp)
        mask_bgc_vegp[p] || continue
        c = patch.column[p]

        # Integrate rootfr and surface_prof over active layer
        rootfr_tot       = 0.0
        surface_prof_tot = 0.0
        jmax = min(max(altmax_lastyear_indx[c], 1), nlevdecomp)
        for j in 1:jmax
            rootfr_tot       += cinput_rootfr[p, j] * dzsoi_decomp[j]
            surface_prof_tot += surface_prof[j]      * dzsoi_decomp[j]
        end

        if (altmax_lastyear_indx[c] > 0) && (rootfr_tot > 0.0) && (surface_prof_tot > 0.0)
            # Normal case: distribute over active layer
            for j in 1:jmax
                froot_prof[p, j] = cinput_rootfr[p, j] / rootfr_tot
                croot_prof[p, j] = cinput_rootfr[p, j] / rootfr_tot
                # Surface processes get the shallower profile
                leaf_prof[p, j]  = surface_prof[j] / surface_prof_tot
                stem_prof[p, j]  = surface_prof[j] / surface_prof_tot
            end
        else
            # Fully frozen or no roots: everything in top layer
            froot_prof[p, 1] = 1.0 / dzsoi_decomp[1]
            croot_prof[p, 1] = 1.0 / dzsoi_decomp[1]
            leaf_prof[p, 1]  = 1.0 / dzsoi_decomp[1]
            stem_prof[p, 1]  = 1.0 / dzsoi_decomp[1]
        end
    end

    # ------------------------------------------------------------------
    # 5. Aggregate root profile from patch to column using wtcol.
    # ------------------------------------------------------------------
    for p in eachindex(mask_bgc_vegp)
        mask_bgc_vegp[p] || continue
        c = patch.column[p]
        for j in 1:nlevdecomp
            col_cinput_rootfr[c, j] += cinput_rootfr[p, j] * patch.wtcol[p]
        end
    end

    # ------------------------------------------------------------------
    # 6. Column-native profiles: Ndep and Nfix.
    #    FATES columns use surface profile for both; non-FATES columns
    #    use col_cinput_rootfr for Nfix and surface_prof for Ndep.
    # ------------------------------------------------------------------
    for c in eachindex(mask_bgc_soilc)
        mask_bgc_soilc[c] || continue

        rootfr_tot       = 0.0
        surface_prof_tot = 0.0

        if col.is_fates[c]
            # FATES path
            jmax = min(max(altmax_lastyear_indx[c], 1), nlevdecomp)
            for j in 1:jmax
                surface_prof_tot += surface_prof[j] * dzsoi_decomp[j]
            end
            if (altmax_lastyear_indx[c] > 0) && (surface_prof_tot > 0.0)
                for j in 1:jmax
                    nfixation_prof[c, j] = surface_prof[j] / surface_prof_tot
                    ndep_prof[c, j]      = surface_prof[j] / surface_prof_tot
                end
            else
                nfixation_prof[c, 1] = 1.0 / dzsoi_decomp[1]
                ndep_prof[c, 1]      = 1.0 / dzsoi_decomp[1]
            end
        else
            # Non-FATES: redo column integration over active layer
            jmax = min(max(altmax_lastyear_indx[c], 1), nlevdecomp)
            for j in 1:jmax
                rootfr_tot       += col_cinput_rootfr[c, j] * dzsoi_decomp[j]
                surface_prof_tot += surface_prof[j]          * dzsoi_decomp[j]
            end
            if (altmax_lastyear_indx[c] > 0) && (rootfr_tot > 0.0) && (surface_prof_tot > 0.0)
                for j in 1:jmax
                    nfixation_prof[c, j] = col_cinput_rootfr[c, j] / rootfr_tot
                    ndep_prof[c, j]      = surface_prof[j]          / surface_prof_tot
                end
            else
                nfixation_prof[c, 1] = 1.0 / dzsoi_decomp[1]
                ndep_prof[c, 1]      = 1.0 / dzsoi_decomp[1]
            end
        end
    end

    # ------------------------------------------------------------------
    # 7. Integrity check: all column profiles must integrate to 1.
    # ------------------------------------------------------------------
    for c in eachindex(mask_bgc_soilc)
        mask_bgc_soilc[c] || continue
        ndep_prof_sum      = 0.0
        nfixation_prof_sum = 0.0
        for j in 1:nlevdecomp
            ndep_prof_sum      += ndep_prof[c, j]      * dzsoi_decomp[j]
            nfixation_prof_sum += nfixation_prof[c, j] * dzsoi_decomp[j]
        end
        if (abs(ndep_prof_sum - 1.0) > delta) || (abs(nfixation_prof_sum - 1.0) > delta)
            error("soil_bgc_vertical_profile!: column profile integrity check failed " *
                  "at c=$c  ndep_prof_sum=$ndep_prof_sum  nfixation_prof_sum=$nfixation_prof_sum  " *
                  "altmax_lastyear_indx=$(altmax_lastyear_indx[c])")
        end
    end

    # ------------------------------------------------------------------
    # 8. Integrity check: all patch profiles must integrate to 1.
    # ------------------------------------------------------------------
    for p in eachindex(mask_bgc_vegp)
        mask_bgc_vegp[p] || continue
        froot_prof_sum = 0.0
        croot_prof_sum = 0.0
        leaf_prof_sum  = 0.0
        stem_prof_sum  = 0.0
        for j in 1:nlevdecomp
            froot_prof_sum += froot_prof[p, j] * dzsoi_decomp[j]
            croot_prof_sum += croot_prof[p, j] * dzsoi_decomp[j]
            leaf_prof_sum  += leaf_prof[p, j]  * dzsoi_decomp[j]
            stem_prof_sum  += stem_prof[p, j]  * dzsoi_decomp[j]
        end
        if (abs(froot_prof_sum - 1.0) > delta) || (abs(croot_prof_sum - 1.0) > delta) ||
           (abs(leaf_prof_sum  - 1.0) > delta) || (abs(stem_prof_sum  - 1.0) > delta)
            error("soil_bgc_vertical_profile!: patch profile integrity check failed " *
                  "at p=$p  froot=$froot_prof_sum  croot=$croot_prof_sum  " *
                  "leaf=$leaf_prof_sum  stem=$stem_prof_sum")
        end
    end

    return nothing
end
