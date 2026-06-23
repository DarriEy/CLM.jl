# FatesRadiationDriveMod.jl
# ============================================================================
# FATES Tier-F Batch 17 — port of FatesRadiationDriveMod.F90 (504 lines).
#
# The host-facing canopy-radiation DRIVER. Two public entry points:
#
#   * `FatesNormalizedCanopyRadiation(nsites, sites, bc_in, bc_out)` — performs
#     the normalized (per-unit-downwelling) surface-albedo / canopy-scattering
#     solve for every vegetated patch at every site, selecting the Norman or the
#     two-stream solver, and packs albedo / absorbed-fraction / transmitted-to-
#     ground fractions into the `bc_out.*_parb` arrays the host land model reads
#     for the next time step.
#
#   * `FatesSunShadeFracs(nsites, sites, bc_in, bc_out)` — derives the patch
#     sunlit fraction (`fsun_pa`), sunlit/shaded LAI, and the actual (forced)
#     absorbed-PAR profiles through the canopy from the saved normalized
#     scattering profiles, for the photosynthesis driver.
#
# Translation notes / conventions (CLAUDE.md):
#   * Fortran variable names preserved (s/ifp/ib/cl/ft/iv/icol...).
#   * The Fortran loops the age-ordered patch linked list (oldest_patch ->
#     younger). `ifp` is the 1-based VEGETATED-patch index used to slice the
#     per-patch bc arrays; it is incremented ONLY for non-bareground patches
#     (matching upstream — bareground patches in SP/nocomp mode must NOT bump it
#     or the indexing breaks).
#   * `radiation_model`, `dinc_vai`, `dlower_vai` come from the live parameter
#     instance `ed_params()`; `numpft`, `hlm_hio_ignore_val` are module Refs.
#   * The per-patch bc_out arrays are (ifp, ib); Fortran passes the whole
#     `(ifp,:)` slice into PatchNormanRadiation — here we build length-num_swb
#     scratch band vectors and copy back, since the ported PatchNormanRadiation
#     takes `Vector{Float64}` band vectors (it fills them in place).
#
# Upstream-Fortran quirks preserved (see inline `# QUIRK`):
#   * `preserve_b4b = true` (module-level): keeps several night/no-canopy branches
#     turned OFF to stay bit-for-bit with the baseline two-stream tests. The
#     `.not.preserve_b4b` branches are reproduced but gated identically.
#   * The night branch (`solar_zenith_flag == false`) leaves bc_out untouched
#     while preserve_b4b is true (the trivial albedo=1 solution is the disabled
#     branch).
#   * `debug = false`: the two-stream radiation-balance check / Dump / absorbed-
#     fraction>1 endrun all live under `if(debug)` and are no-ops here.
# ============================================================================

# Local debug flag (mirrors Fortran `debug = .false.`)
const fates_radiationdrive_debug = false

# QUIRK: preserve_b4b will be removed upstream soon. Kept here to prevent round-
# off differences in the baseline two-stream tests (RGK 12-27-23).
const fates_radiationdrive_preserve_b4b = true

# ============================================================================

"""
    FatesNormalizedCanopyRadiation(nsites, sites, bc_in, bc_out)

Perform normalized (per unit downwelling radiative forcing) radiation scattering
of the vegetation canopy for every vegetated patch at each site. The call is
normalized because the host model wants an albedo for the next time step before
it knows the absolute beam/diffuse forcing; the normalized scattering and
absorption profiles are saved on the patch and scaled by the forcing later
(diagnostics, heating rates, absorbed leaf PAR for photosynthesis).

Fills, per vegetated patch `ifp` and shortwave band `ib`, the `bc_out` arrays:
`albd_parb`/`albi_parb` (direct/diffuse albedo), `fabd_parb`/`fabi_parb`
(fraction absorbed by canopy), `ftdd_parb`/`ftid_parb`/`ftii_parb` (down
direct->direct / direct->diffuse / diffuse->diffuse below canopy).
"""
function FatesNormalizedCanopyRadiation(nsites::Integer,
                                        sites::AbstractVector,
                                        bc_in::AbstractVector,
                                        bc_out::AbstractVector)

    radiation_model = ed_params().radiation_model

    # -------------------------------------------------------------------------
    # TODO (mv, 2014-10-29) the filter here is different than below — needed to
    # keep the VOC emissions bit-for-bit; FATES is still incompatible with the
    # VOC emission module (RGK 2016-08-06). Preserved as a comment for fidelity.
    # -------------------------------------------------------------------------

    for s in 1:nsites

        ifp = 0
        currentPatch = sites[s].oldest_patch
        while currentPatch !== nothing

            # Do not do albedo calculations for the bare-ground patch in SP mode,
            # and (more importantly) do not iterate ifp or it will mess up the
            # indexing wherein ifp=1 is the first vegetated patch.
            if currentPatch.nocomp_pft_label != nocomp_bareground

                ifp += 1

                # ---- Zero diagnostics on the patch --------------------------
                fill!(currentPatch.f_sun,      0.0)
                fill!(currentPatch.fabd_sun_z, 0.0)
                fill!(currentPatch.fabd_sha_z, 0.0)
                fill!(currentPatch.fabi_sun_z, 0.0)
                fill!(currentPatch.fabi_sha_z, 0.0)
                fill!(currentPatch.fabd,       0.0)
                fill!(currentPatch.fabi,       0.0)
                fill!(currentPatch.nrmlzd_parprof_pft_dir_z, 0.0)
                fill!(currentPatch.nrmlzd_parprof_pft_dif_z, 0.0)

                fill!(currentPatch.rad_error, hlm_hio_ignore_val[])

                currentPatch.solar_zenith_flag  = bc_in[s].filter_vegzen_pa[ifp]
                currentPatch.solar_zenith_angle = bc_in[s].coszen_pa[ifp]
                for ib in 1:num_swb
                    currentPatch.gnd_alb_dif[ib] = bc_in[s].albgr_dif_rb[ib]
                    currentPatch.gnd_alb_dir[ib] = bc_in[s].albgr_dir_rb[ib]
                end
                currentPatch.fcansno = bc_in[s].fcansno_pa[ifp]

                if radiation_model == twostr_solver
                    CanopyPrep!(currentPatch.twostr, bc_in[s].fcansno_pa[ifp])
                    ZenithPrep!(currentPatch.twostr, bc_in[s].coszen_pa[ifp])
                end

                if !currentPatch.solar_zenith_flag
                    # Sun below horizon, trivial solution.
                    # Note (RGK-MLO): Investigate twilight mechanics for non-zero
                    # diffuse radiation when cosz<=0.
                    # QUIRK: while preserve_b4b is true the trivial albedo=1
                    # solution is DISABLED to keep the baseline tests bit-for-bit;
                    # bc_out is simply left as the host initialized it.
                    if !fates_radiationdrive_preserve_b4b
                        for ib in 1:num_swb
                            bc_out[s].albd_parb[ifp, ib] = 1.0
                            bc_out[s].albi_parb[ifp, ib] = 1.0
                            bc_out[s].fabi_parb[ifp, ib] = 0.0
                            bc_out[s].fabd_parb[ifp, ib] = 0.0
                            bc_out[s].ftdd_parb[ifp, ib] = 0.0
                            bc_out[s].ftid_parb[ifp, ib] = 0.0
                            bc_out[s].ftii_parb[ifp, ib] = 0.0
                        end
                    end
                else

                    for ib in 1:num_swb
                        bc_out[s].albd_parb[ifp, ib] = 0.0
                        bc_out[s].albi_parb[ifp, ib] = 0.0
                        bc_out[s].fabi_parb[ifp, ib] = 0.0
                        bc_out[s].fabd_parb[ifp, ib] = 0.0
                        bc_out[s].ftdd_parb[ifp, ib] = 1.0
                        bc_out[s].ftid_parb[ifp, ib] = 1.0
                        bc_out[s].ftii_parb[ifp, ib] = 1.0
                    end

                    if maximum(@view currentPatch.nrad[1, :]) == 0
                        # No leaf layers in this patch — effectively bare ground.
                        for ib in 1:num_swb
                            bc_out[s].fabd_parb[ifp, ib] = 0.0
                            bc_out[s].fabi_parb[ifp, ib] = 0.0
                        end
                        fill!(currentPatch.rad_error, 0.0)

                        for ib in 1:num_swb
                            bc_out[s].albd_parb[ifp, ib] = bc_in[s].albgr_dir_rb[ib]
                            bc_out[s].albi_parb[ifp, ib] = bc_in[s].albgr_dif_rb[ib]
                            bc_out[s].ftdd_parb[ifp, ib] = 1.0
                            bc_out[s].ftid_parb[ifp, ib] = 0.0
                            bc_out[s].ftii_parb[ifp, ib] = 1.0
                        end

                    else

                        if radiation_model == norman_solver

                            # The ported PatchNormanRadiation fills band vectors
                            # (length num_swb) in place; copy back to the (ifp,:)
                            # slices afterward.
                            albd = zeros(num_swb); albi = zeros(num_swb)
                            fabd = zeros(num_swb); fabi = zeros(num_swb)
                            ftdd = zeros(num_swb); ftid = zeros(num_swb)
                            ftii = zeros(num_swb)

                            PatchNormanRadiation(currentPatch,
                                                 albd,   # Surface Albedo direct
                                                 albi,   # Surface Albedo (indirect) diffuse
                                                 fabd,   # Frac direct absorbed by canopy / unit incident
                                                 fabi,   # Frac diffuse absorbed by canopy / unit incident
                                                 ftdd,   # Down direct flux below canopy / unit direct at top
                                                 ftid,   # Down diffuse flux below canopy / unit direct at top
                                                 ftii)   # Down diffuse flux below canopy / unit diffuse at top

                            for ib in 1:num_swb
                                bc_out[s].albd_parb[ifp, ib] = albd[ib]
                                bc_out[s].albi_parb[ifp, ib] = albi[ib]
                                bc_out[s].fabd_parb[ifp, ib] = fabd[ib]
                                bc_out[s].fabi_parb[ifp, ib] = fabi[ib]
                                bc_out[s].ftdd_parb[ifp, ib] = ftdd[ib]
                                bc_out[s].ftid_parb[ifp, ib] = ftid[ib]
                                bc_out[s].ftii_parb[ifp, ib] = ftii[ib]
                            end

                        elseif radiation_model == twostr_solver

                            twostr = currentPatch.twostr

                            # (CanopyPrep!/ZenithPrep! already done above)
                            for ib in 1:num_swb

                                twostr.band[ib].albedo_grnd_diff = bc_in[s].albgr_dif_rb[ib]
                                twostr.band[ib].albedo_grnd_beam = bc_in[s].albgr_dir_rb[ib]

                                (albedo_beam, albedo_diff, consv_err,
                                 frac_abs_can_beam, frac_abs_can_diff,
                                 frac_beam_grnd_beam, frac_diff_grnd_beam,
                                 frac_diff_grnd_diff) =
                                    Solve!(twostr, ib,
                                           normalized_upper_boundary,
                                           1.0, 1.0,
                                           sites[s].taulambda_2str,   # inout (scratch)
                                           sites[s].omega_2str,       # inout (scratch)
                                           sites[s].ipiv_2str)        # inout (scratch)

                                bc_out[s].albd_parb[ifp, ib] = albedo_beam
                                bc_out[s].albi_parb[ifp, ib] = albedo_diff
                                currentPatch.rad_error[ib]   = consv_err
                                bc_out[s].fabd_parb[ifp, ib] = frac_abs_can_beam
                                bc_out[s].fabi_parb[ifp, ib] = frac_abs_can_diff
                                bc_out[s].ftdd_parb[ifp, ib] = frac_beam_grnd_beam
                                bc_out[s].ftid_parb[ifp, ib] = frac_diff_grnd_beam
                                bc_out[s].ftii_parb[ifp, ib] = frac_diff_grnd_diff

                                # QUIRK: the debug radiation-balance check
                                # (CheckPatchRadiationBalance + twostr%Dump +
                                # endrun on absorbed-fraction>1) lives under
                                # `if(debug)`; debug=false -> no-op here.
                                if fates_radiationdrive_debug
                                    twostr.band[ib].Rbeam_atm = 1.0
                                    twostr.band[ib].Rdiff_atm = 1.0
                                    CheckPatchRadiationBalance(currentPatch,
                                                               sites[s].snow_depth, ib,
                                                               bc_out[s].fabd_parb[ifp, ib],
                                                               bc_out[s].fabi_parb[ifp, ib])
                                    twostr.band[ib].Rbeam_atm = fates_unset_r8
                                    twostr.band[ib].Rdiff_atm = fates_unset_r8
                                end
                            end  # ib

                        end  # select radiation_model

                    end  # if_nrad

                end  # if_zenith_flag
            end  # if_notbareground

            currentPatch = currentPatch.younger
        end  # Loop linked-list patches
    end  # Loop sites

    return nothing
end

# ============================================================================

"""
    FatesSunShadeFracs(nsites, sites, bc_in, bc_out)

Compute the patch sunlit fraction (`bc_out.fsun_pa`), sunlit/shaded LAI
(`laisun_pa`/`laisha_pa`), and the absorbed-PAR profiles through the canopy
(`ed_parsun_z`/`ed_parsha_z`, `parprof_pft_dir_z`/`parprof_pft_dif_z`) by
scaling the saved normalized scattering profiles with the actual beam/diffuse
PAR forcing in `bc_in.solad_parb`/`bc_in.solai_parb`. Handles both the Norman
solver path (derive layer sunlit/shaded LAI from `f_sun`) and the two-stream
path (area-weight per-element absorbed radiation from the solver).
"""
function FatesSunShadeFracs(nsites::Integer,
                            sites::AbstractVector,
                            bc_in::AbstractVector,
                            bc_out::AbstractVector)

    radiation_model = ed_params().radiation_model
    dinc_vai   = ed_params().dinc_vai
    dlower_vai = ed_params().dlower_vai
    np = numpft[]

    for s in 1:nsites

        ifp = 0
        cpatch = sites[s].oldest_patch

        while cpatch !== nothing

            # Only for veg patches: do not do albedo calcs for the bare-ground
            # patch in SP mode, and (more importantly) do not iterate ifp or it
            # will mess up the indexing wherein ifp=1 is the first vegetated patch.
            if cpatch.nocomp_pft_label != nocomp_bareground

                ifp += 1

                # ---- Initialize diagnostics --------------------------------
                fill!(cpatch.ed_parsun_z, 0.0)
                fill!(cpatch.ed_parsha_z, 0.0)
                fill!(cpatch.ed_laisun_z, 0.0)
                fill!(cpatch.ed_laisha_z, 0.0)
                fill!(cpatch.parprof_pft_dir_z, 0.0)
                fill!(cpatch.parprof_pft_dif_z, 0.0)

                bc_out[s].fsun_pa[ifp] = 0.0

                # QUIRK: preserve_b4b — these two are disabled to stay bit-for-bit
                if !fates_radiationdrive_preserve_b4b
                    bc_out[s].laisun_pa[ifp] = 0.0
                    bc_out[s].laisha_pa[ifp] = calc_areaindex(cpatch, "elai")
                end

                sunlai = 0.0
                shalai = 0.0

                if radiation_model == norman_solver

                    # Loop over canopy layers to calculate laisun_z and laisha_z
                    # for each layer; derive canopy laisun/laisha/fsun from layer
                    # sums. (cpatch.f_sun is set by FatesNormalizedCanopyRadiation.)
                    for cl in 1:cpatch.ncl_p
                        for ft in 1:np
                            # QUIRK: preserve_b4b selects the layerwise path
                            if !fates_radiationdrive_preserve_b4b
                                for iv in 1:cpatch.nrad[cl, ft]
                                    sunlai += cpatch.elai_profile[cl, ft, iv] *
                                              cpatch.f_sun[cl, ft, iv]
                                    shalai += cpatch.elai_profile[cl, ft, iv]
                                end
                            else
                                for iv in 1:cpatch.nrad[cl, ft]
                                    cpatch.ed_laisun_z[cl, ft, iv] =
                                        cpatch.elai_profile[cl, ft, iv] *
                                        cpatch.f_sun[cl, ft, iv]
                                    cpatch.ed_laisha_z[cl, ft, iv] =
                                        cpatch.elai_profile[cl, ft, iv] *
                                        (1.0 - cpatch.f_sun[cl, ft, iv])
                                end
                                # needed for the VOC emissions, etc.
                                for iv in 1:cpatch.nrad[cl, ft]
                                    sunlai += cpatch.ed_laisun_z[cl, ft, iv]
                                    shalai += cpatch.ed_laisha_z[cl, ft, iv]
                                end
                            end
                        end
                    end

                    if !fates_radiationdrive_preserve_b4b
                        shalai = shalai - sunlai
                    end

                    if sunlai + shalai > 0.0
                        bc_out[s].fsun_pa[ifp] = sunlai / (sunlai + shalai)
                    else
                        bc_out[s].fsun_pa[ifp] = 0.0
                    end

                    if fates_radiationdrive_debug
                        if bc_out[s].fsun_pa[ifp] > 1.0
                            @warn "too much leaf area in profile" bc_out[s].fsun_pa[ifp] sunlai shalai
                        end
                    end

                    elai = calc_areaindex(cpatch, "elai")

                    bc_out[s].laisun_pa[ifp] = elai * bc_out[s].fsun_pa[ifp]
                    bc_out[s].laisha_pa[ifp] = elai * (1.0 - bc_out[s].fsun_pa[ifp])

                    # Absorbed PAR profile through canopy. If sun/shade big-leaf
                    # code, nrad=1 and fluxes from SurfaceAlbedo are canopy
                    # integrated so layer values equal big-leaf values.
                    for cl in 1:cpatch.ncl_p
                        for ft in 1:np
                            for iv in 1:cpatch.nrad[cl, ft]
                                cpatch.ed_parsun_z[cl, ft, iv] =
                                    bc_in[s].solad_parb[ifp, ipar] * cpatch.fabd_sun_z[cl, ft, iv] +
                                    bc_in[s].solai_parb[ifp, ipar] * cpatch.fabi_sun_z[cl, ft, iv]

                                cpatch.ed_parsha_z[cl, ft, iv] =
                                    bc_in[s].solad_parb[ifp, ipar] * cpatch.fabd_sha_z[cl, ft, iv] +
                                    bc_in[s].solai_parb[ifp, ipar] * cpatch.fabi_sha_z[cl, ft, iv]
                            end
                        end
                    end

                    # Convert normalized radiation error from fraction to W/m2
                    for ib in 1:num_swb
                        cpatch.rad_error[ib] = cpatch.rad_error[ib] *
                            (bc_in[s].solad_parb[ifp, ib] + bc_in[s].solai_parb[ifp, ib])
                    end

                    # Output the actual PAR profiles through canopy (diagnostics)
                    for cl in 1:cpatch.ncl_p
                        for ft in 1:np
                            for iv in 1:cpatch.nrad[cl, ft]
                                cpatch.parprof_pft_dir_z[cl, ft, iv] =
                                    bc_in[s].solad_parb[ifp, ipar] *
                                    cpatch.nrmlzd_parprof_pft_dir_z[idirect, cl, ft, iv] +
                                    bc_in[s].solai_parb[ifp, ipar] *
                                    cpatch.nrmlzd_parprof_pft_dir_z[idiffuse, cl, ft, iv]

                                cpatch.parprof_pft_dif_z[cl, ft, iv] =
                                    bc_in[s].solad_parb[ifp, ipar] *
                                    cpatch.nrmlzd_parprof_pft_dif_z[idirect, cl, ft, iv] +
                                    bc_in[s].solai_parb[ifp, ipar] *
                                    cpatch.nrmlzd_parprof_pft_dif_z[idiffuse, cl, ft, iv]
                            end
                        end
                    end

                else  # two-stream solver

                    # If there is no sun out, we have a trivial solution.
                    if cpatch.solar_zenith_flag

                        for ib in 1:num_swb
                            cpatch.twostr.band[ib].Rbeam_atm = bc_in[s].solad_parb[ifp, ib]
                            cpatch.twostr.band[ib].Rdiff_atm = bc_in[s].solai_parb[ifp, ib]
                        end

                        # Fraction of the canopy area associated with each pft and
                        # layer (used for weighting diagnostics).
                        area_vlpfcl = zeros(nlevleaf, maxpft, nclmax)
                        fill!(cpatch.f_sun, 0.0)

                        (fsun, laisun, laisha) = FatesPatchFSun(cpatch)
                        bc_out[s].fsun_pa[ifp]   = fsun
                        bc_out[s].laisun_pa[ifp] = laisun
                        bc_out[s].laisha_pa[ifp] = laisha

                        twostr = cpatch.twostr

                        for cl in 1:twostr.n_lyr
                            for icol in 1:twostr.n_col[cl]

                                scelg = twostr.scelg[cl, icol]
                                ft = scelg.pft
                                if ft > 0  # not air
                                    area_frac = scelg.area
                                    vai = scelg.sai + scelg.lai
                                    # nv = lowest VAI bin whose lower edge exceeds vai
                                    nv = _minloc_mask_gt(dlower_vai, vai)
                                    for iv in 1:nv

                                        vai_top = dlower_vai[iv] - dinc_vai[iv]
                                        vai_bot = min(dlower_vai[iv], scelg.sai + scelg.lai)

                                        cpatch.parprof_pft_dir_z[cl, ft, iv] =
                                            cpatch.parprof_pft_dir_z[cl, ft, iv] +
                                            area_frac * GetRb(twostr, cl, icol, ivis, vai_top)
                                        cpatch.parprof_pft_dif_z[cl, ft, iv] =
                                            cpatch.parprof_pft_dif_z[cl, ft, iv] +
                                            area_frac * GetRdDn(twostr, cl, icol, ivis, vai_top) +
                                            area_frac * GetRdUp(twostr, cl, icol, ivis, vai_top)

                                        (Rb_abs, Rd_abs, Rd_abs_leaf, Rb_abs_leaf,
                                         R_abs_stem, R_abs_snow, leaf_sun_frac) =
                                            GetAbsRad(twostr, cl, icol, ipar, vai_top, vai_bot)

                                        cpatch.f_sun[cl, ft, iv] =
                                            cpatch.f_sun[cl, ft, iv] + area_frac * leaf_sun_frac
                                        cpatch.ed_parsun_z[cl, ft, iv] =
                                            cpatch.ed_parsun_z[cl, ft, iv] +
                                            area_frac * (Rd_abs_leaf * leaf_sun_frac + Rb_abs_leaf)
                                        cpatch.ed_parsha_z[cl, ft, iv] =
                                            cpatch.ed_parsha_z[cl, ft, iv] +
                                            area_frac * Rd_abs_leaf * (1.0 - leaf_sun_frac)

                                        area_vlpfcl[iv, ft, cl] =
                                            area_vlpfcl[iv, ft, cl] + area_frac
                                    end  # iv
                                end  # not air
                            end  # icol

                            for ft in 1:np
                                for iv in 1:cpatch.nleaf[cl, ft]
                                    if area_vlpfcl[iv, ft, cl] < nearzero
                                        break  # exit do_iv
                                    end
                                    cpatch.parprof_pft_dir_z[cl, ft, iv] =
                                        cpatch.parprof_pft_dir_z[cl, ft, iv] / area_vlpfcl[iv, ft, cl]
                                    cpatch.parprof_pft_dif_z[cl, ft, iv] =
                                        cpatch.parprof_pft_dif_z[cl, ft, iv] / area_vlpfcl[iv, ft, cl]
                                    cpatch.f_sun[cl, ft, iv] =
                                        cpatch.f_sun[cl, ft, iv] / area_vlpfcl[iv, ft, cl]
                                    cpatch.ed_parsun_z[cl, ft, iv] =
                                        cpatch.ed_parsun_z[cl, ft, iv] / area_vlpfcl[iv, ft, cl]
                                    cpatch.ed_parsha_z[cl, ft, iv] =
                                        cpatch.ed_parsha_z[cl, ft, iv] / area_vlpfcl[iv, ft, cl]
                                end  # iv
                            end  # ft
                        end  # cl

                    end  # if_zenithflag

                end  # if_norm_twostr

            end  # if_notbareground

            cpatch = cpatch.younger
        end  # patches
    end  # sites

    return nothing
end

# ----------------------------------------------------------------------------
# Helper: Fortran `minloc(dlower_vai, DIM=1, MASK=(dlower_vai>vai))` — returns
# the index of the FIRST element (lowest index, since the array is monotone)
# satisfying dlower_vai > vai. Fortran minloc returns 0 when no element matches.
# ----------------------------------------------------------------------------
function _minloc_mask_gt(arr::AbstractVector{<:Real}, thresh::Real)
    best_idx = 0
    best_val = 0.0
    for i in eachindex(arr)
        if arr[i] > thresh
            if best_idx == 0 || arr[i] < best_val
                best_idx = i
                best_val = arr[i]
            end
        end
    end
    return best_idx
end
