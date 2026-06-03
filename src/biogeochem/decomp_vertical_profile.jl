# ==========================================================================
# Ported from: src/soilbiogeochem/SoilBiogeochemVerticalProfileMod.F90
# Calculate vertical profiles for distributing soil and litter C and N
#
# Public functions:
#   soil_bgc_vertical_profile!  -- main vertical profile computation
#
# Module-level parameter:
#   surfprof_exp  -- steepness of surface profile (1/e-folding depth) [1/m]
#
# GPU note (Phase B kernelization):
#   The per-timestep compute body runs WHOLE-FUNCTION on Metal. Every
#   per-patch / per-column compute loop is a KernelAbstractions kernel
#   (one thread per patch / per column; internal j/k loops stay sequential
#   in-thread). The single patch->column aggregation is a true many-to-one
#   scatter (`col_cinput_rootfr[col,j] += cinput_rootfr[p,j]*wtcol[p]`,
#   col = patch_column[p]) done with `_scatter_add!` (atomic on GPU, plain
#   += on CPU/Dual) inside a per-PATCH kernel — the proven Phase-A pattern.
#   The two integrity checks reduce per-column / per-patch profile sums into
#   small device vectors via reduction kernels (own-index, race-free), then
#   the host inspects only those length-nc / length-np scalar sums to throw.
# ==========================================================================

# How steep the profile is for surface components (1/e_folding depth) (1/m)
const SURFPROF_EXP = 10.0

# --------------------------------------------------------------------------
# Step 1: surface profile (shallow exponential, optionally zeroed below
# bedrock). One thread per decomposition level j. `use_bedrock` is a Bool
# flag; `surfprof_exp` and `zmin_bedrock` are passed at the state precision.
# --------------------------------------------------------------------------
@kernel function _dvp_surface_prof_kernel!(surface_prof, @Const(zsoi), @Const(dzsoi_decomp),
        surfprof_exp, use_bedrock::Bool, zmin_bedrock)
    j = @index(Global)
    @inbounds begin
        T = eltype(surface_prof)
        sp = exp(-surfprof_exp * zsoi[j]) / dzsoi_decomp[j]
        if use_bedrock && zsoi[j] > zmin_bedrock
            sp = zero(T)
        end
        surface_prof[j] = sp
    end
end

# --------------------------------------------------------------------------
# Step 2a: zero all patch output profiles. One thread per patch; internal
# loop over nlevdecomp_full. (Unmasked — Fortran zeros all points.)
# --------------------------------------------------------------------------
@kernel function _dvp_init_patch_prof_kernel!(leaf_prof, froot_prof, croot_prof, stem_prof,
        cinput_rootfr, nlevdecomp_full::Int)
    p = @index(Global)
    @inbounds begin
        T = eltype(leaf_prof)
        for j in 1:nlevdecomp_full
            leaf_prof[p, j]  = zero(T)
            froot_prof[p, j] = zero(T)
            croot_prof[p, j] = zero(T)
            stem_prof[p, j]  = zero(T)
            cinput_rootfr[p, j] = zero(T)
        end
    end
end

# Step 2b: zero all column output profiles + the col aggregation workspace.
@kernel function _dvp_init_col_prof_kernel!(nfixation_prof, ndep_prof, col_cinput_rootfr,
        nlevdecomp_full::Int)
    c = @index(Global)
    @inbounds begin
        T = eltype(nfixation_prof)
        for j in 1:nlevdecomp_full
            nfixation_prof[c, j] = zero(T)
            ndep_prof[c, j]      = zero(T)
            col_cinput_rootfr[c, j] = zero(T)
        end
    end
end

# --------------------------------------------------------------------------
# Step 3+4: per-patch cinput_rootfr, then integrate over active layer and
# build the normalised patch profiles. One thread per patch (own-index),
# internal sequential j loops. noveg patches and frozen/rootless patches
# fall back to the top layer.
# --------------------------------------------------------------------------
@kernel function _dvp_patch_prof_kernel!(froot_prof, croot_prof, leaf_prof, stem_prof,
        cinput_rootfr, @Const(crootfr), @Const(surface_prof), @Const(dzsoi_decomp),
        @Const(mask_bgc_vegp), @Const(patch_itype), @Const(patch_column),
        @Const(altmax_lastyear_indx), noveg::Int, nlevdecomp::Int)
    p = @index(Global)
    @inbounds if mask_bgc_vegp[p]
        T = eltype(cinput_rootfr)
        # --- cinput_rootfr (root fraction normalised by layer thickness) ---
        if patch_itype[p] != noveg
            for j in 1:nlevdecomp
                cinput_rootfr[p, j] = crootfr[p, j] / dzsoi_decomp[j]
            end
        else
            cinput_rootfr[p, 1] = zero(T)
        end

        c = patch_column[p]

        # Integrate rootfr and surface_prof over the active layer.
        rootfr_tot       = zero(T)
        surface_prof_tot = zero(T)
        jmax = min(max(altmax_lastyear_indx[c], 1), nlevdecomp)
        for j in 1:jmax
            rootfr_tot       += cinput_rootfr[p, j] * dzsoi_decomp[j]
            surface_prof_tot += surface_prof[j]      * dzsoi_decomp[j]
        end

        if (altmax_lastyear_indx[c] > 0) && (rootfr_tot > zero(T)) && (surface_prof_tot > zero(T))
            for j in 1:jmax
                froot_prof[p, j] = cinput_rootfr[p, j] / rootfr_tot
                croot_prof[p, j] = cinput_rootfr[p, j] / rootfr_tot
                leaf_prof[p, j]  = surface_prof[j] / surface_prof_tot
                stem_prof[p, j]  = surface_prof[j] / surface_prof_tot
            end
        else
            froot_prof[p, 1] = one(T) / dzsoi_decomp[1]
            croot_prof[p, 1] = one(T) / dzsoi_decomp[1]
            leaf_prof[p, 1]  = one(T) / dzsoi_decomp[1]
            stem_prof[p, 1]  = one(T) / dzsoi_decomp[1]
        end
    end
end

# --------------------------------------------------------------------------
# Step 5: patch->column weighted aggregation of the root profile. TRUE
# many-to-one scatter (many patches -> one column): one thread per patch,
# `_scatter_add!` (atomic on GPU, += on CPU/Dual) into col_cinput_rootfr.
# col_cinput_rootfr MUST be zeroed first (done in Step 2b).
# --------------------------------------------------------------------------
@kernel function _dvp_scatter_col_rootfr_kernel!(col_cinput_rootfr, @Const(cinput_rootfr),
        @Const(mask_bgc_vegp), @Const(patch_column), @Const(patch_wtcol), nlevdecomp::Int)
    p = @index(Global)
    @inbounds if mask_bgc_vegp[p]
        c = patch_column[p]
        w = patch_wtcol[p]
        for j in 1:nlevdecomp
            _scatter_add!(col_cinput_rootfr, c, j, cinput_rootfr[p, j] * w)
        end
    end
end

# --------------------------------------------------------------------------
# Step 6: column-native Ndep / Nfix profiles. One thread per column; the
# FATES vs non-FATES branch is selected by the per-column is_fates flag.
# Internal sequential j loops.
# --------------------------------------------------------------------------
@kernel function _dvp_col_prof_kernel!(nfixation_prof, ndep_prof, @Const(col_cinput_rootfr),
        @Const(surface_prof), @Const(dzsoi_decomp), @Const(mask_bgc_soilc),
        @Const(col_is_fates), @Const(altmax_lastyear_indx), nlevdecomp::Int)
    c = @index(Global)
    @inbounds if mask_bgc_soilc[c]
        T = eltype(nfixation_prof)
        rootfr_tot       = zero(T)
        surface_prof_tot = zero(T)
        jmax = min(max(altmax_lastyear_indx[c], 1), nlevdecomp)

        if col_is_fates[c]
            for j in 1:jmax
                surface_prof_tot += surface_prof[j] * dzsoi_decomp[j]
            end
            if (altmax_lastyear_indx[c] > 0) && (surface_prof_tot > zero(T))
                for j in 1:jmax
                    nfixation_prof[c, j] = surface_prof[j] / surface_prof_tot
                    ndep_prof[c, j]      = surface_prof[j] / surface_prof_tot
                end
            else
                nfixation_prof[c, 1] = one(T) / dzsoi_decomp[1]
                ndep_prof[c, 1]      = one(T) / dzsoi_decomp[1]
            end
        else
            for j in 1:jmax
                rootfr_tot       += col_cinput_rootfr[c, j] * dzsoi_decomp[j]
                surface_prof_tot += surface_prof[j]          * dzsoi_decomp[j]
            end
            if (altmax_lastyear_indx[c] > 0) && (rootfr_tot > zero(T)) && (surface_prof_tot > zero(T))
                for j in 1:jmax
                    nfixation_prof[c, j] = col_cinput_rootfr[c, j] / rootfr_tot
                    ndep_prof[c, j]      = surface_prof[j]          / surface_prof_tot
                end
            else
                nfixation_prof[c, 1] = one(T) / dzsoi_decomp[1]
                ndep_prof[c, 1]      = one(T) / dzsoi_decomp[1]
            end
        end
    end
end

# --------------------------------------------------------------------------
# Step 7: per-column reduction of ndep / nfix profile sums (own-column, in
# ascending-j order). Writes into small device vectors the host inspects.
# --------------------------------------------------------------------------
@kernel function _dvp_col_sums_kernel!(ndep_sum, nfix_sum, @Const(ndep_prof),
        @Const(nfixation_prof), @Const(dzsoi_decomp), @Const(mask_bgc_soilc), nlevdecomp::Int)
    c = @index(Global)
    @inbounds if mask_bgc_soilc[c]
        T = eltype(ndep_sum)
        ds = zero(T)
        ns = zero(T)
        for j in 1:nlevdecomp
            ds += ndep_prof[c, j]      * dzsoi_decomp[j]
            ns += nfixation_prof[c, j] * dzsoi_decomp[j]
        end
        ndep_sum[c] = ds
        nfix_sum[c] = ns
    end
end

# --------------------------------------------------------------------------
# Step 8: per-patch reduction of froot/croot/leaf/stem profile sums.
# --------------------------------------------------------------------------
@kernel function _dvp_patch_sums_kernel!(froot_sum, croot_sum, leaf_sum, stem_sum,
        @Const(froot_prof), @Const(croot_prof), @Const(leaf_prof), @Const(stem_prof),
        @Const(dzsoi_decomp), @Const(mask_bgc_vegp), nlevdecomp::Int)
    p = @index(Global)
    @inbounds if mask_bgc_vegp[p]
        T = eltype(froot_sum)
        fs = zero(T); cs = zero(T); ls = zero(T); ss = zero(T)
        for j in 1:nlevdecomp
            fs += froot_prof[p, j] * dzsoi_decomp[j]
            cs += croot_prof[p, j] * dzsoi_decomp[j]
            ls += leaf_prof[p, j]  * dzsoi_decomp[j]
            ss += stem_prof[p, j]  * dzsoi_decomp[j]
        end
        froot_sum[p] = fs
        croot_sum[p] = cs
        leaf_sum[p]  = ls
        stem_sum[p]  = ss
    end
end

"""
    soil_bgc_vertical_profile!(
        soilbiogeochem_state::SoilBiogeochemStateData,
        active_layer::ActiveLayerData,
        soilstate::SoilStateData,
        col::ColumnData,
        patch::PatchData;
        mask_bgc_soilc,
        mask_bgc_vegp,
        nlevdecomp::Int,
        nlevdecomp_full::Int,
        dzsoi_decomp,
        zsoi,
        use_bedrock::Bool=false,
        zmin_bedrock=ZMIN_BEDROCK)

Calculate vertical profiles for distributing soil and litter C and N.

This routine computes profiles of leaves, stems, fine roots, coarse roots,
N fixation, and N deposition as a function of depth, normalised over the
active layer.  Where the soil is fully frozen or rootless, all inputs are
placed in the top decomposition layer.

The FATES branch uses the surface e-folding profile for both fixation and
deposition (partially because fixation may be free-living).

Ported from `SoilBiogeochemVerticalProfile` in
`SoilBiogeochemVerticalProfileMod.F90`.

The whole per-timestep compute body runs as KernelAbstractions kernels (CPU
and Apple Metal); see the module header for the kernel layout.

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
        mask_bgc_soilc::AbstractVector{Bool},
        mask_bgc_vegp::AbstractVector{Bool},
        nlevdecomp::Int,
        nlevdecomp_full::Int,
        dzsoi_decomp::AbstractVector{<:Real},
        zsoi::AbstractVector{<:Real},
        use_bedrock::Bool=false,
        zmin_bedrock::Real=ZMIN_BEDROCK)

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

    # --- Local workspace (allocated on the device of the state arrays) ---
    FT = eltype(crootfr)
    cinput_rootfr     = similar(crootfr, FT, np, nlevdecomp_full)
    col_cinput_rootfr = similar(nfixation_prof, FT, nc, nlevdecomp_full)
    surface_prof      = similar(crootfr, FT, nlevdecomp)

    # Integrity-check tolerance. Float64 keeps the exact Fortran value (1e-10);
    # the check is a guard only (it mutates no output), so on a reduced-precision
    # device path (Float32 on Metal) we widen it past that precision's round-off
    # to avoid spurious trips — the byte-identity guarantee is on the profiles.
    delta = FT === Float64 ? 1.0e-10 : 100 * eps(FT)

    # Scalar parameters passed into kernels must carry the state precision so no
    # Float64 literal/arg reaches a Float32-only backend (Metal).
    surfprof_exp_T = FT(SURFPROF_EXP)
    zmin_bedrock_T = FT(zmin_bedrock)

    # ------------------------------------------------------------------
    # 1. Surface profile (one thread per decomposition level).
    # ------------------------------------------------------------------
    _launch!(_dvp_surface_prof_kernel!, surface_prof, zsoi, dzsoi_decomp,
             surfprof_exp_T, use_bedrock, zmin_bedrock_T)

    # ------------------------------------------------------------------
    # 2. Zero all output profiles + the cinput / col-aggregation workspace.
    # ------------------------------------------------------------------
    _launch!(_dvp_init_patch_prof_kernel!, leaf_prof, froot_prof, croot_prof, stem_prof,
             cinput_rootfr, nlevdecomp_full; ndrange = np)
    _launch!(_dvp_init_col_prof_kernel!, nfixation_prof, ndep_prof, col_cinput_rootfr,
             nlevdecomp_full; ndrange = nc)

    # ------------------------------------------------------------------
    # 3+4. Per-patch cinput_rootfr and normalised patch profiles.
    # ------------------------------------------------------------------
    _launch!(_dvp_patch_prof_kernel!, froot_prof, croot_prof, leaf_prof, stem_prof,
             cinput_rootfr, crootfr, surface_prof, dzsoi_decomp,
             mask_bgc_vegp, patch.itype, patch.column,
             altmax_lastyear_indx, noveg, nlevdecomp; ndrange = np)

    # ------------------------------------------------------------------
    # 5. Aggregate root profile patch -> column (true scatter).
    # ------------------------------------------------------------------
    _launch!(_dvp_scatter_col_rootfr_kernel!, col_cinput_rootfr, cinput_rootfr,
             mask_bgc_vegp, patch.column, patch.wtcol, nlevdecomp; ndrange = np)

    # ------------------------------------------------------------------
    # 6. Column-native Ndep / Nfix profiles (FATES vs non-FATES).
    # ------------------------------------------------------------------
    _launch!(_dvp_col_prof_kernel!, nfixation_prof, ndep_prof, col_cinput_rootfr,
             surface_prof, dzsoi_decomp, mask_bgc_soilc,
             col.is_fates, altmax_lastyear_indx, nlevdecomp; ndrange = nc)

    # ------------------------------------------------------------------
    # 7. Integrity check: all column profiles must integrate to 1.
    #    Reduce per-column sums in a kernel, inspect the small sum vectors.
    # ------------------------------------------------------------------
    ndep_prof_sum_col = similar(nfixation_prof, FT, nc)
    nfix_prof_sum_col = similar(nfixation_prof, FT, nc)
    fill!(ndep_prof_sum_col, one(FT))
    fill!(nfix_prof_sum_col, one(FT))
    _launch!(_dvp_col_sums_kernel!, ndep_prof_sum_col, nfix_prof_sum_col, ndep_prof,
             nfixation_prof, dzsoi_decomp, mask_bgc_soilc, nlevdecomp; ndrange = nc)
    let ndsum = Array(ndep_prof_sum_col), nfsum = Array(nfix_prof_sum_col),
        amask = Array(mask_bgc_soilc), aalt = Array(altmax_lastyear_indx)
        for c in 1:nc
            amask[c] || continue
            if (abs(ndsum[c] - 1.0) > delta) || (abs(nfsum[c] - 1.0) > delta)
                error("soil_bgc_vertical_profile!: column profile integrity check failed " *
                      "at c=$c  ndep_prof_sum=$(ndsum[c])  nfixation_prof_sum=$(nfsum[c])  " *
                      "altmax_lastyear_indx=$(aalt[c])")
            end
        end
    end

    # ------------------------------------------------------------------
    # 8. Integrity check: all patch profiles must integrate to 1.
    # ------------------------------------------------------------------
    froot_prof_sum_p = similar(crootfr, FT, np)
    croot_prof_sum_p = similar(crootfr, FT, np)
    leaf_prof_sum_p  = similar(crootfr, FT, np)
    stem_prof_sum_p  = similar(crootfr, FT, np)
    fill!(froot_prof_sum_p, one(FT))
    fill!(croot_prof_sum_p, one(FT))
    fill!(leaf_prof_sum_p,  one(FT))
    fill!(stem_prof_sum_p,  one(FT))
    _launch!(_dvp_patch_sums_kernel!, froot_prof_sum_p, croot_prof_sum_p, leaf_prof_sum_p,
             stem_prof_sum_p, froot_prof, croot_prof, leaf_prof, stem_prof,
             dzsoi_decomp, mask_bgc_vegp, nlevdecomp; ndrange = np)
    let fsum = Array(froot_prof_sum_p), csum = Array(croot_prof_sum_p),
        lsum = Array(leaf_prof_sum_p), ssum = Array(stem_prof_sum_p),
        pmask = Array(mask_bgc_vegp)
        for p in 1:np
            pmask[p] || continue
            if (abs(fsum[p] - 1.0) > delta) || (abs(csum[p] - 1.0) > delta) ||
               (abs(lsum[p] - 1.0) > delta) || (abs(ssum[p] - 1.0) > delta)
                error("soil_bgc_vertical_profile!: patch profile integrity check failed " *
                      "at p=$p  froot=$(fsum[p])  croot=$(csum[p])  " *
                      "leaf=$(lsum[p])  stem=$(ssum[p])")
            end
        end
    end

    return nothing
end
