# FatesNormanRadMod.jl
# ============================================================================
# FATES Tier-F Batch 16 — port of FatesNormanRadMod.F90 (983 lines).
#
# The Norman (1979) multi-layer canopy radiative-transfer solver. For each
# patch it computes — per canopy layer (L) × PFT (ft) × leaf layer (iv) —
# the absorbed / reflected / transmitted direct + diffuse shortwave radiation,
# the patch albedo, and the sunlit/shaded leaf fractions over the canopy
# profile. The public entry point is `PatchNormanRadiation`.
#
# Translation notes / conventions (CLAUDE.md):
#   * Fortran variable names preserved for traceability.
#   * The patch radiation profile arrays (elai_profile, esai_profile,
#     canopy_area_profile, f_sun, fab{d,i}_{sun,sha}_z, nrmlzd_parprof_*) live
#     on `fates_patch_type` (FatesPatchMod.jl) and are indexed (L, ft, iv) —
#     the SAME order as Fortran (L, ft, iv). nrad / canopy_mask are (L, ft).
#   * Local work arrays are sized exactly as Fortran: (nclmax, maxpft, nlevleaf)
#     (and +num_swb where 4-D), so the `iv+1` accesses up to nrad+1 stay in
#     bounds (nrad <= nlevleaf).
#   * The `do radtype = 1, num_rad_stream_types` outer loop runs the whole
#     scattering solution twice: once for a unit of direct beam (radtype =
#     idirect) and once for a unit of diffuse (radtype = idiffuse), selected by
#     the binary switches forc_dir / forc_dif.
#
# Upstream-Fortran quirks preserved (see inline comments marked QUIRK):
#   * The albedo "error fix-up": the within-scheme conservation residual is
#     added back onto the output albedo for |error| in (1e-9, ∞). This is
#     deliberate in CTSM-FATES to prevent host-model crashes from small
#     occasional energy-balance errors; we reproduce it faithfully.
#   * There is a Fortran typo on the lai_change>0 down-flux branch:
#     `(ftweight(L,ft,1)-ftweight(L,ft,iv-1)/ftweight(L,ft,1))` — the division
#     binds before the subtraction (missing parens). Reproduced byte-for-byte.
#   * `rho_snow` / `tau_snow` / `albice` are module-level public canopy-snow
#     reflectance parameters, stored here (not on the parameter file) in common
#     with the ice albedo, exactly as upstream.
# ============================================================================

# Canopy-snow reflectance model parameters (waveband: 1=vis, 2=nir).
# QUIRK: stored here in common with the land-ice albedo rather than on the
# parameter file, mirroring upstream FATES.
const albice   = [0.80, 0.55]   # albedo land ice by waveband (1=vis, 2=nir)
const rho_snow = [0.80, 0.55]   # canopy snow reflectance by waveband
const tau_snow = [0.01, 0.01]   # canopy snow transmittance by waveband

"""
    PatchNormanRadiation(currentPatch,
                         albd_parb_out, albi_parb_out,
                         fabd_parb_out, fabi_parb_out,
                         ftdd_parb_out, ftid_parb_out, ftii_parb_out)

Perform the Norman radiation scattering for a single FATES patch. The seven
`*_parb_out` arguments are `Vector{Float64}` of length `num_swb` and are filled
in place (the Fortran `(ifp, ib)` slice for this patch):

  * `albd_parb_out` — direct-beam albedo
  * `albi_parb_out` — diffuse albedo
  * `fabd_parb_out` — fraction of direct radiation absorbed by the canopy
  * `fabi_parb_out` — fraction of diffuse radiation absorbed by the canopy
  * `ftdd_parb_out` — direct -> soil as direct
  * `ftid_parb_out` — direct -> soil as diffuse
  * `ftii_parb_out` — diffuse -> soil as diffuse

`currentPatch::fates_patch_type` is updated in place (f_sun, fab*_z,
nrmlzd_parprof_*, fabd/fabi, gnd-flux & sabs diagnostics, rad_error).
"""
function PatchNormanRadiation(currentPatch::fates_patch_type,
                              albd_parb_out::Vector{Float64},
                              albi_parb_out::Vector{Float64},
                              fabd_parb_out::Vector{Float64},
                              fabi_parb_out::Vector{Float64},
                              ftdd_parb_out::Vector{Float64},
                              ftid_parb_out::Vector{Float64},
                              ftii_parb_out::Vector{Float64})

    # ---- associate: live FATES PFT parameter table (EDPftvarcon_inst) --------
    pftcon = EDPftvarcon_inst[]
    rhol           = pftcon.rhol            # leaf reflectance   (ft, ib)
    rhos           = pftcon.rhos            # stem reflectance   (ft, ib)
    taul           = pftcon.taul            # leaf transmittance (ft, ib)
    taus           = pftcon.taus            # stem transmittance (ft, ib)
    xl             = pftcon.xl              # leaf/stem orientation index (ft)
    clumping_index = pftcon.clumping_index  # self-occlusion clumping (ft)

    npft = numpft[]

    tolerance = 1.0e-9

    # Binary switches used to turn radiation streams on/off (indexed by radtype)
    forc_dir = (1.0, 0.0)   # forc_dir(num_rad_stream_types)
    forc_dif = (0.0, 1.0)   # forc_dif(num_rad_stream_types)

    # ---- local work arrays (sized exactly as Fortran) ------------------------
    ftweight        = zeros(nclmax, maxpft, nlevleaf)
    k_dir           = zeros(maxpft)
    tr_dir_z        = zeros(nclmax, maxpft, nlevleaf)
    tr_dif_z        = zeros(nclmax, maxpft, nlevleaf)
    weighted_dir_tr = zeros(nclmax)
    weighted_fsun   = zeros(nclmax)
    weighted_dif_ratio = zeros(nclmax, num_swb)
    weighted_dif_down  = zeros(nclmax)
    weighted_dif_up    = zeros(nclmax)
    refl_dif    = zeros(nclmax, maxpft, nlevleaf, num_swb)
    tran_dif    = zeros(nclmax, maxpft, nlevleaf, num_swb)
    dif_ratio   = zeros(nclmax, maxpft, nlevleaf, num_swb)
    Dif_dn      = zeros(nclmax, maxpft, nlevleaf)
    Dif_up      = zeros(nclmax, maxpft, nlevleaf)
    lai_change  = zeros(nclmax, maxpft, nlevleaf)
    f_abs       = zeros(nclmax, maxpft, nlevleaf, num_swb)
    rho_layer   = zeros(nclmax, maxpft, nlevleaf, num_swb)
    tau_layer   = zeros(nclmax, maxpft, nlevleaf, num_swb)
    f_abs_leaf  = zeros(nclmax, maxpft, nlevleaf, num_swb)
    Abs_dir_z   = zeros(maxpft, nlevleaf)
    Abs_dif_z   = zeros(maxpft, nlevleaf)
    abs_rad     = zeros(num_swb)
    phi1b       = zeros(maxpft)
    phi2b       = zeros(maxpft)
    lai_reduction = zeros(nclmax)

    # ---- initialize the output arrays ----------------------------------------
    for ib in 1:num_swb
        albd_parb_out[ib] = 0.0
        albi_parb_out[ib] = 0.0
        fabd_parb_out[ib] = 0.0
        fabi_parb_out[ib] = 0.0
        ftdd_parb_out[ib] = 1.0
        ftid_parb_out[ib] = 1.0
        ftii_parb_out[ib] = 1.0
    end

    # =========================================================================
    # Layer-level optical properties: present PFT/canopy-layer combinations
    # =========================================================================
    for L in 1:currentPatch.ncl_p
        for ft in 1:npft
            currentPatch.canopy_mask[L, ft] = 0
            for iv in 1:currentPatch.nrad[L, ft]
                if currentPatch.canopy_area_profile[L, ft, iv] > 0.0
                    currentPatch.canopy_mask[L, ft] = 1

                    if currentPatch.elai_profile[L, ft, iv] +
                       currentPatch.esai_profile[L, ft, iv] > 0.0
                        frac_lai = currentPatch.elai_profile[L, ft, iv] /
                            (currentPatch.elai_profile[L, ft, iv] +
                             currentPatch.esai_profile[L, ft, iv])
                    else
                        frac_lai = 1.0
                    end
                    frac_sai = 1.0 - frac_lai

                    for ib in 1:num_swb  # vis, nir
                        rho_layer[L, ft, iv, ib] = frac_lai * rhol[ft, ib] +
                                                   frac_sai * rhos[ft, ib]
                        tau_layer[L, ft, iv, ib] = frac_lai * taul[ft, ib] +
                                                   frac_sai * taus[ft, ib]

                        # adjust reflectance / transmittance for canopy snow
                        rho_layer[L, ft, iv, ib] =
                            rho_layer[L, ft, iv, ib] * (1.0 - currentPatch.fcansno) +
                            rho_snow[ib] * currentPatch.fcansno
                        tau_layer[L, ft, iv, ib] =
                            tau_layer[L, ft, iv, ib] * (1.0 - currentPatch.fcansno) +
                            tau_snow[ib] * currentPatch.fcansno

                        # fraction of incoming light absorbed by leaves or stems
                        f_abs[L, ft, iv, ib] = 1.0 - tau_layer[L, ft, iv, ib] -
                                                     rho_layer[L, ft, iv, ib]

                        # fraction of the absorbed light absorbed by leaves
                        f_abs_leaf[L, ft, iv, ib] =
                            (1.0 - currentPatch.fcansno) * frac_lai *
                            (1.0 - rhol[ft, ib] - taul[ft, ib]) / f_abs[L, ft, iv, ib]
                    end  # ib
                end
            end  # iv
        end  # ft
    end  # L

    # =========================================================================
    # Direct-beam extinction coefficient k_dir (PFT specific)
    # =========================================================================
    cosz = max(0.001, currentPatch.solar_zenith_angle)  # copied from previous code
    for ft in 1:npft
        sb = (90.0 - (acos(cosz) * 180.0 / pi_const)) * (pi_const / 180.0)
        phi1b[ft] = 0.5 - 0.633 * xl[ft] - 0.330 * xl[ft] * xl[ft]
        phi2b[ft] = 0.877 * (1.0 - 2.0 * phi1b[ft])  # 0=horiz leaves, 1=vert
        gdir = phi1b[ft] + phi2b[ft] * sin(sb)
        # how much direct light penetrates a single unit of lai?
        k_dir[ft] = clumping_index[ft] * gdir / sin(sb)
    end  # ft

    # =========================================================================
    # Once for one unit of diffuse, once for one unit of direct radiation
    # =========================================================================
    for radtype in 1:num_rad_stream_types

        # ---- ftweight = canopy area profile (already area-corrected) ---------
        fill!(ftweight, 0.0)
        for L in 1:currentPatch.ncl_p
            for ft in 1:npft
                for iv in 1:currentPatch.nrad[L, ft]
                    ftweight[L, ft, iv] = currentPatch.canopy_area_profile[L, ft, iv]
                end
            end
        end

        # ---- direct/diffuse transmittance & sunlit fraction (top-down) -------
        for L in 1:currentPatch.ncl_p  # start at the top canopy layer (1 = top)

            weighted_dir_tr[L] = 0.0
            weighted_fsun[L]   = 0.0
            for ib in 1:num_swb
                weighted_dif_ratio[L, ib] = 0.0
            end

            for ft in 1:npft
                if currentPatch.canopy_mask[L, ft] == 1  # only if leaves present

                    # --- diffuse transmittance tr_dif_z (9-angle sky integral) ---
                    for iv in 1:currentPatch.nrad[L, ft]
                        tr_dif_z[L, ft, iv] = 0.0
                    end
                    for iv in 1:currentPatch.nrad[L, ft]
                        for j in 1:9
                            angle = (5.0 + (j - 1) * 10.0) * pi_const / 180.0
                            gdir = phi1b[ft] + phi2b[ft] * sin(angle)
                            tr_dif_z[L, ft, iv] = tr_dif_z[L, ft, iv] +
                                exp(-clumping_index[ft] * gdir / sin(angle) *
                                    (currentPatch.elai_profile[L, ft, iv] +
                                     currentPatch.esai_profile[L, ft, iv])) *
                                sin(angle) * cos(angle)
                        end
                        tr_dif_z[L, ft, iv] = tr_dif_z[L, ft, iv] * 2.0 *
                                              (10.0 * pi_const / 180.0)
                    end

                    # --- direct beam transmittance tr_dir_z (cumulative LAI) ---
                    if L == 1
                        tr_dir_z[L, ft, 1] = 1.0
                    else
                        tr_dir_z[L, ft, 1] = weighted_dir_tr[L - 1]
                    end
                    laisum = 0.0
                    for iv in 1:currentPatch.nrad[L, ft]
                        laisum += currentPatch.elai_profile[L, ft, iv] +
                                  currentPatch.esai_profile[L, ft, iv]
                        lai_change[L, ft, iv] = 0.0
                        if ftweight[L, ft, iv + 1] > 0.0 &&
                           ftweight[L, ft, iv + 1] < ftweight[L, ft, iv]
                            # partly-empty leaf layer: some flux goes straight through
                            lai_change[L, ft, iv] = ftweight[L, ft, iv] -
                                                    ftweight[L, ft, iv + 1]
                        end

                        # light coming straight through the canopy
                        if L == 1
                            tr_dir_z[L, ft, iv + 1] = exp(-k_dir[ft] * laisum) *
                                (ftweight[L, ft, iv] / ftweight[L, ft, 1])
                        else
                            tr_dir_z[L, ft, iv + 1] = weighted_dir_tr[L - 1] *
                                exp(-k_dir[ft] * laisum) *
                                (ftweight[L, ft, iv] / ftweight[L, ft, 1])
                        end

                        if iv == 1
                            # top layer
                            tr_dir_z[L, ft, iv + 1] = tr_dir_z[L, ft, iv + 1] +
                                tr_dir_z[L, ft, iv] *
                                ((ftweight[L, ft, 1] - ftweight[L, ft, iv]) /
                                 ftweight[L, ft, 1])
                        else
                            # lai_change(iv-1) affects the light onto iv+1
                            if lai_change[L, ft, iv - 1] > 0.0
                                tr_dir_z[L, ft, iv + 1] = tr_dir_z[L, ft, iv + 1] +
                                    tr_dir_z[L, ft, iv] * lai_change[L, ft, iv - 1] /
                                    ftweight[L, ft, 1]
                                tr_dir_z[L, ft, iv + 1] = tr_dir_z[L, ft, iv + 1] +
                                    tr_dir_z[L, ft, iv - 1] *
                                    (ftweight[L, ft, 1] - ftweight[L, ft, iv - 1]) /
                                    ftweight[L, ft, 1]
                            else
                                tr_dir_z[L, ft, iv + 1] = tr_dir_z[L, ft, iv + 1] +
                                    tr_dir_z[L, ft, iv] *
                                    ((ftweight[L, ft, 1] - ftweight[L, ft, iv]) /
                                     ftweight[L, ft, 1])
                            end
                        end
                    end

                    # weighted contribution from this PFT column
                    weighted_dir_tr[L] = weighted_dir_tr[L] +
                        tr_dir_z[L, ft, currentPatch.nrad[L, ft] + 1] *
                        ftweight[L, ft, 1]

                    # --- sunlit / shaded fraction of leaf layer ---
                    laisum = 0.0
                    for iv in 1:currentPatch.nrad[L, ft]
                        # cumulative lai at center of layer
                        if iv == 1
                            laisum = 0.5 * (currentPatch.elai_profile[L, ft, iv] +
                                            currentPatch.esai_profile[L, ft, iv])
                        else
                            laisum += currentPatch.elai_profile[L, ft, iv] +
                                      currentPatch.esai_profile[L, ft, iv]
                        end

                        if L == 1  # top canopy layer
                            currentPatch.f_sun[L, ft, iv] = exp(-k_dir[ft] * laisum) *
                                (ftweight[L, ft, iv] / ftweight[L, ft, 1])
                        else
                            currentPatch.f_sun[L, ft, iv] = weighted_fsun[L - 1] *
                                exp(-k_dir[ft] * laisum) *
                                (ftweight[L, ft, iv] / ftweight[L, ft, 1])
                        end

                        if iv > 1
                            if lai_change[L, ft, iv - 1] > 0.0
                                currentPatch.f_sun[L, ft, iv] =
                                    currentPatch.f_sun[L, ft, iv] +
                                    currentPatch.f_sun[L, ft, iv] *
                                    lai_change[L, ft, iv - 1] / ftweight[L, ft, 1]
                                currentPatch.f_sun[L, ft, iv] =
                                    currentPatch.f_sun[L, ft, iv] +
                                    currentPatch.f_sun[L, ft, iv - 1] *
                                    (ftweight[L, ft, 1] - ftweight[L, ft, iv - 1]) /
                                    ftweight[L, ft, 1]
                            else
                                currentPatch.f_sun[L, ft, iv] =
                                    currentPatch.f_sun[L, ft, iv] +
                                    currentPatch.f_sun[L, ft, iv - 1] *
                                    (ftweight[L, ft, 1] - ftweight[L, ft, iv]) /
                                    ftweight[L, ft, 1]
                            end
                        end
                    end  # iv

                    weighted_fsun[L] = weighted_fsun[L] +
                        currentPatch.f_sun[L, ft, currentPatch.nrad[L, ft]] *
                        ftweight[L, ft, 1]
                    # QUIRK (FTWA): first-layer ftweight is used as a proxy for the
                    # whole column — a slight source of error noted upstream.
                end  # present
            end  # ft
        end  # L

        # ---- dif_ratio: upward/forward diffuse flux ratio (bottom-up) --------
        for L in currentPatch.ncl_p:-1:1
            for ft in 1:npft
                if currentPatch.canopy_mask[L, ft] == 1

                    for ib in 1:num_swb  # vis, nir
                        # diffuse reflected / transmitted by a layer
                        for iv in 1:currentPatch.nrad[L, ft]
                            refl_dif[L, ft, iv, ib] =
                                (1.0 - tr_dif_z[L, ft, iv]) * rho_layer[L, ft, iv, ib]
                            tran_dif[L, ft, iv, ib] =
                                (1.0 - tr_dif_z[L, ft, iv]) * tau_layer[L, ft, iv, ib] +
                                tr_dif_z[L, ft, iv]
                        end

                        # soil diffuse reflectance (down -> up ratio)
                        iv = currentPatch.nrad[L, ft] + 1
                        if L == currentPatch.ncl_p  # nearest the soil
                            dif_ratio[L, ft, iv, ib] = currentPatch.gnd_alb_dif[ib]
                        else
                            dif_ratio[L, ft, iv, ib] = weighted_dif_ratio[L + 1, ib]
                        end

                        # canopy layers, upward from soil
                        for iv in currentPatch.nrad[L, ft]:-1:1
                            dif_ratio[L, ft, iv, ib] =
                                dif_ratio[L, ft, iv + 1, ib] *
                                tran_dif[L, ft, iv, ib] * tran_dif[L, ft, iv, ib] /
                                (1.0 - dif_ratio[L, ft, iv + 1, ib] *
                                       refl_dif[L, ft, iv, ib]) +
                                refl_dif[L, ft, iv, ib]
                            dif_ratio[L, ft, iv, ib] = dif_ratio[L, ft, iv, ib] *
                                ftweight[L, ft, iv] / ftweight[L, ft, 1]
                            dif_ratio[L, ft, iv, ib] = dif_ratio[L, ft, iv, ib] +
                                dif_ratio[L, ft, iv + 1, ib] *
                                (ftweight[L, ft, 1] - ftweight[L, ft, iv]) /
                                ftweight[L, ft, 1]
                        end
                        weighted_dif_ratio[L, ib] = weighted_dif_ratio[L, ib] +
                            dif_ratio[L, ft, 1, ib] * ftweight[L, ft, 1]
                    end  # num_swb
                end  # canopy_mask
            end  # ft
        end  # L

        # =====================================================================
        # Per-band flux solution
        # =====================================================================
        for ib in 1:num_swb

            currentPatch.rad_error[ib] = 0.0

            fill!(Dif_dn, 0.0)
            fill!(Dif_up, 0.0)

            # ---- first estimate: downward diffuse flux (top-down) ------------
            for L in 1:currentPatch.ncl_p
                weighted_dif_down[L] = 0.0
                for ft in 1:npft
                    if currentPatch.canopy_mask[L, ft] == 1
                        if L == 1
                            Dif_dn[L, ft, 1] = forc_dif[radtype]
                        else
                            Dif_dn[L, ft, 1] = weighted_dif_down[L - 1]
                        end
                        for iv in 1:currentPatch.nrad[L, ft]
                            denom = refl_dif[L, ft, iv, ib] * dif_ratio[L, ft, iv, ib]
                            denom = 1.0 - denom
                            Dif_dn[L, ft, iv + 1] = Dif_dn[L, ft, iv] *
                                tran_dif[L, ft, iv, ib] / denom *
                                ftweight[L, ft, iv] / ftweight[L, ft, 1]
                            if iv > 1
                                if lai_change[L, ft, iv - 1] > 0.0
                                    Dif_dn[L, ft, iv + 1] = Dif_dn[L, ft, iv + 1] +
                                        Dif_dn[L, ft, iv] *
                                        lai_change[L, ft, iv - 1] / ftweight[L, ft, 1]
                                    # QUIRK: upstream typo — missing parens so the
                                    # division binds before the subtraction:
                                    #   (ftw(1) - ftw(iv-1)/ftw(1))   [sic]
                                    Dif_dn[L, ft, iv + 1] = Dif_dn[L, ft, iv + 1] +
                                        Dif_dn[L, ft, iv - 1] *
                                        (ftweight[L, ft, 1] -
                                         ftweight[L, ft, iv - 1] / ftweight[L, ft, 1])
                                else
                                    Dif_dn[L, ft, iv + 1] = Dif_dn[L, ft, iv + 1] +
                                        Dif_dn[L, ft, iv] *
                                        (ftweight[L, ft, 1] - ftweight[L, ft, iv]) /
                                        ftweight[L, ft, 1]
                                end
                            else
                                Dif_dn[L, ft, iv + 1] = Dif_dn[L, ft, iv + 1] +
                                    Dif_dn[L, ft, iv] *
                                    (ftweight[L, ft, 1] - ftweight[L, ft, iv]) /
                                    ftweight[L, ft, 1]
                            end
                        end
                        weighted_dif_down[L] = weighted_dif_down[L] +
                            Dif_dn[L, ft, currentPatch.nrad[L, ft] + 1] *
                            ftweight[L, ft, 1]
                    end  # present
                end  # ft
                if L == currentPatch.ncl_p && currentPatch.ncl_p > 1
                    # incomplete understorey: add radiation through canopy gaps
                    weighted_dif_down[L] = weighted_dif_down[L] +
                        weighted_dif_down[L - 1] *
                        (1.0 - sum(@view ftweight[L, 1:npft, 1]))
                end
            end  # L

            # ---- first estimate: upward diffuse flux (bottom-up) -------------
            for L in currentPatch.ncl_p:-1:1
                weighted_dif_up[L] = 0.0
                for ft in 1:npft
                    if currentPatch.canopy_mask[L, ft] == 1
                        iv = currentPatch.nrad[L, ft] + 1
                        if L == currentPatch.ncl_p  # bottom layer
                            Dif_up[L, ft, iv] = currentPatch.gnd_alb_dif[ib] *
                                                Dif_dn[L, ft, iv]
                        else
                            Dif_up[L, ft, iv] = weighted_dif_up[L + 1]
                        end

                        for iv in currentPatch.nrad[L, ft]:-1:1
                            if lai_change[L, ft, iv] > 0.0
                                Dif_up[L, ft, iv] = dif_ratio[L, ft, iv, ib] *
                                    Dif_dn[L, ft, iv] *
                                    ftweight[L, ft, iv] / ftweight[L, ft, 1]
                                Dif_up[L, ft, iv] = Dif_up[L, ft, iv] +
                                    Dif_up[L, ft, iv + 1] * tran_dif[L, ft, iv, ib] *
                                    lai_change[L, ft, iv] / ftweight[L, ft, 1]
                                Dif_up[L, ft, iv] = Dif_up[L, ft, iv] +
                                    Dif_up[L, ft, iv + 1] *
                                    (ftweight[L, ft, 1] - ftweight[L, ft, iv]) /
                                    ftweight[L, ft, 1]
                            else
                                Dif_up[L, ft, iv] = dif_ratio[L, ft, iv, ib] *
                                    Dif_dn[L, ft, iv] * ftweight[L, ft, iv]
                                Dif_up[L, ft, iv] = Dif_up[L, ft, iv] +
                                    Dif_up[L, ft, iv + 1] * (1.0 - ftweight[L, ft, iv])
                            end
                        end

                        weighted_dif_up[L] = weighted_dif_up[L] +
                            Dif_up[L, ft, 1] * ftweight[L, ft, 1]
                    end  # present
                end  # ft
                if L == currentPatch.ncl_p && currentPatch.ncl_p > 1
                    # incomplete understorey: radiation up through canopy gaps
                    # diffuse -> diffuse
                    weighted_dif_up[L] = weighted_dif_up[L] +
                        (1.0 - sum(@view ftweight[L, 1:npft, 1])) *
                        weighted_dif_down[L - 1] * currentPatch.gnd_alb_dif[ib]
                    # direct -> diffuse
                    weighted_dif_up[L] = weighted_dif_up[L] +
                        forc_dir[radtype] * weighted_dir_tr[L - 1] *
                        (1.0 - sum(@view ftweight[L, 1:npft, 1])) *
                        currentPatch.gnd_alb_dir[ib]
                end
            end  # L

            # ---- iterative forward/upward diffuse incl. scattered direct -----
            irep = 1
            iter = 0
            while irep == 1 && iter < 50
                iter += 1
                irep = 0

                # downward (top-down)
                for L in 1:currentPatch.ncl_p
                    weighted_dif_down[L] = 0.0
                    for ft in 1:npft
                        if currentPatch.canopy_mask[L, ft] == 1
                            if L == 1
                                Dif_dn[L, ft, 1] = forc_dif[radtype]
                            else
                                Dif_dn[L, ft, 1] = weighted_dif_down[L - 1]
                            end
                            down_rad = 0.0
                            for iv in 1:currentPatch.nrad[L, ft]
                                down_rad = Dif_dn[L, ft, iv] * tran_dif[L, ft, iv, ib] +
                                    Dif_up[L, ft, iv + 1] * refl_dif[L, ft, iv, ib]
                                # plus direct beam intercepted & transmitted
                                down_rad = down_rad + forc_dir[radtype] *
                                    tr_dir_z[L, ft, iv] *
                                    (1.0 - exp(-k_dir[ft] *
                                        (currentPatch.elai_profile[L, ft, iv] +
                                         currentPatch.esai_profile[L, ft, iv]))) *
                                    tau_layer[L, ft, iv, ib]
                                # spread over the whole of incomplete layers
                                down_rad = down_rad *
                                    (ftweight[L, ft, iv] / ftweight[L, ft, 1])
                                if iv > 1
                                    if lai_change[L, ft, iv - 1] > 0.0
                                        down_rad = down_rad + Dif_dn[L, ft, iv] *
                                            lai_change[L, ft, iv - 1] / ftweight[L, ft, 1]
                                        down_rad = down_rad + Dif_dn[L, ft, iv - 1] *
                                            (ftweight[L, ft, 1] - ftweight[L, ft, iv - 1]) /
                                            ftweight[L, ft, 1]
                                    else
                                        down_rad = down_rad + Dif_dn[L, ft, iv] *
                                            (ftweight[L, ft, 1] - ftweight[L, ft, iv]) /
                                            ftweight[L, ft, 1]
                                    end
                                else
                                    down_rad = down_rad + Dif_dn[L, ft, iv] *
                                        (ftweight[L, ft, 1] - ftweight[L, ft, iv]) /
                                        ftweight[L, ft, 1]
                                end

                                if abs(down_rad - Dif_dn[L, ft, iv + 1]) > tolerance
                                    irep = 1
                                end
                                Dif_dn[L, ft, iv + 1] = down_rad
                            end  # iv

                            weighted_dif_down[L] = weighted_dif_down[L] +
                                Dif_dn[L, ft, currentPatch.nrad[L, ft] + 1] *
                                ftweight[L, ft, 1]
                        end  # present
                    end  # ft
                    if L == currentPatch.ncl_p && currentPatch.ncl_p > 1
                        weighted_dif_down[L] = weighted_dif_down[L] +
                            weighted_dif_down[L - 1] *
                            (1.0 - sum(@view ftweight[L, 1:npft, 1]))
                    end
                end  # L

                # upward (top-down — note Fortran goes 1:NCL_p here)
                for L in 1:currentPatch.ncl_p
                    weighted_dif_up[L] = 0.0
                    for ft in 1:npft
                        if currentPatch.canopy_mask[L, ft] == 1
                            iv = currentPatch.nrad[L, ft] + 1
                            if L == currentPatch.ncl_p  # reflect off soil
                                Dif_up[L, ft, iv] = Dif_dn[L, ft, iv] *
                                    currentPatch.gnd_alb_dif[ib] +
                                    forc_dir[radtype] * tr_dir_z[L, ft, iv] *
                                    currentPatch.gnd_alb_dir[ib]
                            else
                                Dif_up[L, ft, iv] = weighted_dif_up[L + 1]
                            end

                            for iv in currentPatch.nrad[L, ft]:-1:1
                                up_rad = Dif_dn[L, ft, iv] * refl_dif[L, ft, iv, ib]
                                up_rad = up_rad + forc_dir[radtype] *
                                    tr_dir_z[L, ft, iv] *
                                    (1.0 - exp(-k_dir[ft] *
                                        (currentPatch.elai_profile[L, ft, iv] +
                                         currentPatch.esai_profile[L, ft, iv]))) *
                                    rho_layer[L, ft, iv, ib]
                                up_rad = up_rad + Dif_up[L, ft, iv + 1] *
                                    tran_dif[L, ft, iv, ib]
                                up_rad = up_rad * ftweight[L, ft, iv] / ftweight[L, ft, 1]
                                up_rad = up_rad + Dif_up[L, ft, iv + 1] *
                                    (ftweight[L, ft, 1] - ftweight[L, ft, iv]) /
                                    ftweight[L, ft, 1]
                                # NB the lower-layer flux is homogenized, so we
                                # don't consider lai_change here.
                                if abs(up_rad - Dif_up[L, ft, iv]) > tolerance
                                    irep = 1
                                end
                                Dif_up[L, ft, iv] = up_rad
                            end  # iv

                            weighted_dif_up[L] = weighted_dif_up[L] +
                                Dif_up[L, ft, 1] * ftweight[L, ft, 1]
                        end  # present
                    end  # ft
                    if L == currentPatch.ncl_p && currentPatch.ncl_p > 1
                        weighted_dif_up[L] = weighted_dif_up[L] +
                            (1.0 - sum(@view ftweight[L, 1:npft, 1])) *
                            weighted_dif_down[L - 1] * currentPatch.gnd_alb_dif[ib]
                        weighted_dif_up[L] = weighted_dif_up[L] +
                            forc_dir[radtype] * weighted_dir_tr[L - 1] *
                            (1.0 - sum(@view ftweight[L, 1:npft, 1])) *
                            currentPatch.gnd_alb_dir[ib]
                    end
                end  # L
            end  # while iter

            # ---- absorbed flux densities -------------------------------------
            abs_rad[ib] = 0.0
            tr_soili = 0.0
            tr_soild = 0.0

            for L in 1:currentPatch.ncl_p
                fill!(Abs_dir_z, 0.0)
                fill!(Abs_dif_z, 0.0)
                for ft in 1:npft
                    if currentPatch.canopy_mask[L, ft] == 1
                        # absorbed direct & diffuse for leaf layers
                        for iv in 1:currentPatch.nrad[L, ft]
                            Abs_dir_z[ft, iv] = ftweight[L, ft, iv] *
                                forc_dir[radtype] * tr_dir_z[L, ft, iv] *
                                (1.0 - exp(-k_dir[ft] *
                                    (currentPatch.elai_profile[L, ft, iv] +
                                     currentPatch.esai_profile[L, ft, iv]))) *
                                f_abs[L, ft, iv, ib]
                            Abs_dif_z[ft, iv] = ftweight[L, ft, iv] *
                                ((Dif_dn[L, ft, iv] + Dif_up[L, ft, iv + 1]) *
                                 (1.0 - tr_dif_z[L, ft, iv]) * f_abs[L, ft, iv, ib])
                        end

                        # absorbed direct & diffuse for soil
                        if L == currentPatch.ncl_p
                            iv = currentPatch.nrad[L, ft] + 1
                            Abs_dif_z[ft, iv] = ftweight[L, ft, 1] *
                                Dif_dn[L, ft, iv] * (1.0 - currentPatch.gnd_alb_dif[ib])
                            Abs_dir_z[ft, iv] = ftweight[L, ft, 1] *
                                forc_dir[radtype] * tr_dir_z[L, ft, iv] *
                                (1.0 - currentPatch.gnd_alb_dir[ib])
                            tr_soild = tr_soild + ftweight[L, ft, 1] *
                                forc_dir[radtype] * tr_dir_z[L, ft, iv]
                            tr_soili = tr_soili + ftweight[L, ft, 1] *
                                Dif_dn[L, ft, iv]
                        end

                        # absorbed PAR (visible band only): shaded & sunlit
                        if ib == ivis
                            for iv in 1:currentPatch.nrad[L, ft]
                                if radtype == idirect
                                    currentPatch.fabd_sha_z[L, ft, iv] =
                                        Abs_dif_z[ft, iv] *
                                        (1.0 - currentPatch.f_sun[L, ft, iv]) *
                                        f_abs_leaf[L, ft, iv, ib]
                                    currentPatch.fabd_sun_z[L, ft, iv] =
                                        (Abs_dif_z[ft, iv] *
                                         currentPatch.f_sun[L, ft, iv] +
                                         Abs_dir_z[ft, iv]) * f_abs_leaf[L, ft, iv, ib]
                                else
                                    currentPatch.fabi_sha_z[L, ft, iv] =
                                        Abs_dif_z[ft, iv] *
                                        (1.0 - currentPatch.f_sun[L, ft, iv]) *
                                        f_abs_leaf[L, ft, iv, ib]
                                    currentPatch.fabi_sun_z[L, ft, iv] =
                                        Abs_dif_z[ft, iv] *
                                        currentPatch.f_sun[L, ft, iv] *
                                        f_abs_leaf[L, ft, iv, ib]
                                end
                            end
                        end  # ib == ivis

                        # --- sum fluxes ---
                        # solar radiation absorbed by ground
                        iv = currentPatch.nrad[L, ft] + 1
                        if L == currentPatch.ncl_p
                            abs_rad[ib] = abs_rad[ib] +
                                (Abs_dir_z[ft, iv] + Abs_dif_z[ft, iv])
                        end
                        # absorbed by vegetation
                        for iv in 1:currentPatch.nrad[L, ft]
                            if radtype == idirect
                                currentPatch.fabd[ib] = currentPatch.fabd[ib] +
                                    Abs_dir_z[ft, iv] + Abs_dif_z[ft, iv]
                            else
                                currentPatch.fabi[ib] = currentPatch.fabi[ib] +
                                    Abs_dif_z[ft, iv]
                            end
                        end

                        # albedo (top canopy layer)
                        if L == 1
                            if radtype == idirect
                                albd_parb_out[ib] = albd_parb_out[ib] +
                                    Dif_up[L, ft, 1] * ftweight[L, ft, 1]
                            else
                                albi_parb_out[ib] = albi_parb_out[ib] +
                                    Dif_up[L, ft, 1] * ftweight[L, ft, 1]
                            end
                        end

                        # normalized PAR profiles (visible band only)
                        if ib == ivis
                            for iv in 1:currentPatch.nrad[L, ft]
                                currentPatch.nrmlzd_parprof_pft_dir_z[radtype, L, ft, iv] =
                                    forc_dir[radtype] * tr_dir_z[L, ft, iv]
                                currentPatch.nrmlzd_parprof_pft_dif_z[radtype, L, ft, iv] =
                                    Dif_dn[L, ft, iv] + Dif_up[L, ft, iv]
                            end
                        end
                    end  # present
                end  # ft

                if radtype == idirect
                    fabd_parb_out[ib] = currentPatch.fabd[ib]
                else
                    fabi_parb_out[ib] = currentPatch.fabi[ib]
                end

                # radiation absorbed through unfilled part of lower canopy
                if currentPatch.ncl_p > 1 && L == currentPatch.ncl_p
                    abs_rad[ib] = abs_rad[ib] + weighted_dif_down[L - 1] *
                        (1.0 - sum(@view ftweight[L, 1:npft, 1])) *
                        (1.0 - currentPatch.gnd_alb_dif[ib])
                    abs_rad[ib] = abs_rad[ib] + forc_dir[radtype] *
                        weighted_dir_tr[L - 1] *
                        (1.0 - sum(@view ftweight[L, 1:npft, 1])) *
                        (1.0 - currentPatch.gnd_alb_dir[ib])
                    tr_soili = tr_soili + weighted_dif_down[L - 1] *
                        (1.0 - sum(@view ftweight[L, 1:npft, 1]))
                    tr_soild = tr_soild + forc_dir[radtype] * weighted_dir_tr[L - 1] *
                        (1.0 - sum(@view ftweight[L, 1:npft, 1]))
                end

                if radtype == idirect
                    currentPatch.tr_soil_dir[ib]     = tr_soild
                    currentPatch.tr_soil_dir_dif[ib] = tr_soili
                    currentPatch.sabs_dir[ib]        = abs_rad[ib]
                    ftdd_parb_out[ib] = tr_soild
                    ftid_parb_out[ib] = tr_soili
                else
                    currentPatch.tr_soil_dif[ib] = tr_soili
                    currentPatch.sabs_dif[ib]    = abs_rad[ib]
                    ftii_parb_out[ib] = tr_soili
                end
            end  # L

            # ---- conservation check & albedo fix-up --------------------------
            # (the ground-absorption debug error is computed in Fortran only when
            #  debug=.true.; we skip those prints.)

            if radtype == idirect
                error = (forc_dir[radtype] + forc_dif[radtype]) -
                    (fabd_parb_out[ib] + albd_parb_out[ib] + currentPatch.sabs_dir[ib])
            else
                error = (forc_dir[radtype] + forc_dif[radtype]) -
                    (fabi_parb_out[ib] + albi_parb_out[ib] + currentPatch.sabs_dif[ib])
            end

            # ignore the patch radiation error if veg-covered fraction is tiny
            if (currentPatch.total_canopy_area / currentPatch.area) > tolerance
                currentPatch.rad_error[ib] = currentPatch.rad_error[ib] + error *
                    currentPatch.total_canopy_area / currentPatch.area
            end

            fill!(lai_reduction, 0.0)
            for L in 1:currentPatch.ncl_p
                for ft in 1:npft
                    if currentPatch.canopy_mask[L, ft] == 1
                        for iv in 1:currentPatch.nrad[L, ft]
                            if lai_change[L, ft, iv] > 0.0
                                lai_reduction[L] = max(lai_reduction[L],
                                                       lai_change[L, ft, iv])
                            end
                        end
                    end
                end
            end

            # QUIRK: the within-ED radiation tolerance — add the residual back
            # onto the output albedo. Deliberate upstream to prevent host-model
            # crashes from small occasional energy-balance errors.
            if radtype == idirect
                if abs(error) > 1.0e-9 && abs(error) < 0.15
                    albd_parb_out[ib] = albd_parb_out[ib] + error
                end
                if abs(error) > 0.15
                    albd_parb_out[ib] = albd_parb_out[ib] + error
                end
            else
                if abs(error) > 1.0e-9 && abs(error) < 0.15
                    albi_parb_out[ib] = albi_parb_out[ib] + error
                end
                if abs(error) > 0.15
                    albi_parb_out[ib] = albi_parb_out[ib] + error
                end

                # QUIRK: this re-computation of `error` lives inside the diffuse
                # branch in the Fortran (after the diffuse fix-up) and feeds only
                # the debug print, so it has no effect on outputs — preserved for
                # fidelity.
                if radtype == idirect
                    error = (forc_dir[radtype] + forc_dif[radtype]) -
                        (fabd_parb_out[ib] + albd_parb_out[ib] +
                         currentPatch.sabs_dir[ib])
                else
                    error = (forc_dir[radtype] + forc_dif[radtype]) -
                        (fabi_parb_out[ib] + albi_parb_out[ib] +
                         currentPatch.sabs_dif[ib])
                end
            end
        end  # num_swb
    end  # radtype

    return nothing
end
