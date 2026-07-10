# fates_driver_bc.jl
# FATES live-driver wiring — W3/W4 boundary-condition pack/unpack helpers.
#
# These bridge the SoA CLM column/patch state and the FATES site `bc_in`/`bc_out`
# vectors. The MVP maps a single FATES site (index `s`) 1:1 onto a single CLM
# column `c` and its (single) vegetated patch `p`, mirroring the Fortran host glue
# `f2hmap` (sites <-> istsoil columns) in clmfates_interfaceMod.F90, but for the
# carbon-only single-site cold-start built by `clm_fates_init!`.
#
# Each FATES `bc_*` per-patch array is indexed by `ifp` (1-based vegetated-patch
# index). For the NBG cold start there is exactly ONE vegetated patch, so `ifp = 1`.
#
# Pack helpers fill `inst.fates.bc_in[s]` from CLM forcing/state.
# Unpack helpers write `inst.fates.bc_out[s]` back into the CLM structs.
#
# Only the fields needed for the four W3/W4 hooks are touched — kept minimal and
# finite. All packers/unpackers are gated by the caller behind `config.use_fates`.
#
# MULTI-VEG-PATCH / MULTI-SITE GENERALIZATION
# -------------------------------------------
# A real FATES column can carry MORE THAN ONE vegetated patch (disturbance splits a
# patch into a disturbed + an undisturbed remnant), and a run can carry MORE THAN
# ONE FATES site/column. The Fortran host walks the site's age-ordered patch linked
# list (`oldest_patch -> younger`), and for each VEGETATED patch (i.e. skipping a
# `nocomp_bareground` patch, which does not advance the index) maps the FATES patch
# to its HLM patch index via
#
#     p = ifp + col%patchi(c)         (ifp = 1, 2, ..., site%youngest_patch%patchno)
#
# so the HLM bare-ground patch sits at `col%patchi(c)+0` and the vegetated patches
# at `+1, +2, ...`.  The per-patch FATES bc arrays (`*_pa`) are indexed by `ifp`.
#
# `fates_veg_patches(site, c, col)` reproduces exactly that walk: it returns the
# `(ifp, p)` pairs for the site's vegetated patches mapped onto column `c`'s HLM
# patch slots. The single-veg-patch MVP value `p = col.patchi[c]+1` is just the
# `ifp == 1` case of this walk; the carbon-only NBG cold start builds exactly one
# (vegetated, non-bareground) patch, so for that config the walk yields `(1, +1)`
# and the generalized code is byte-identical to the old fixed indexing.

"""
    fates_veg_patches(site, c, col) -> Vector{Tuple{Int,Int}}

Walk a FATES site's age-ordered patch linked list (`oldest_patch -> younger`) and
return the `(ifp, p)` pairs for its VEGETATED patches, where `ifp` is the 1-based
vegetated-patch index (used to slice the patch-dimensioned `bc_*_pa` arrays) and
`p = ifp + col.patchi[c]` is the corresponding HLM patch index in column `c`.

Mirrors the Fortran host `ifp` walk (clmfates_interfaceMod.F90): a
`nocomp_bareground` patch is skipped and does NOT advance `ifp` (so `ifp == 1` is
always the first vegetated patch). The HLM bare-ground patch lives at
`col.patchi[c]+0`; vegetated patches at `+1, +2, ...`.

Patches whose HLM index would exceed the column's reserved patch range
(`col.patchf[c]`) are dropped with a warning — the HLM patch allocation is static,
so a FATES column must reserve enough patch slots for its maximum patch count.
"""
function fates_veg_patches(site, c::Int, col)
    out = Tuple{Int,Int}[]
    pi = col.patchi[c]
    # Reserved HLM patch upper bound. Prefer patchf when it is a valid index
    # (>= patchi); otherwise fall back to patchi + npatches[c] - 1 (bare-ground at
    # +0 plus the reserved veg slots). Some helper-driven setups leave patchf unset
    # (ISPVAL) but do set npatches; this keeps the walk inside the HLM patch
    # allocation either way (matching the old single-patch p = patchi[c]+1 bound).
    pf = col.patchf[c]
    if pf < pi
        npc = col.npatches[c]
        pf = npc >= 1 ? pi + npc - 1 : pi   # +0 bare-ground only ⇒ no veg slot
    end
    ifp = 0
    cp = site.oldest_patch
    while cp !== nothing
        if cp.nocomp_pft_label != nocomp_bareground
            ifp += 1
            p = ifp + pi
            if p <= pf
                push!(out, (ifp, p))
            else
                @warn "FATES column $c has more vegetated patches than reserved HLM " *
                      "patch slots (patchi=$pi patchf=$pf); patch ifp=$ifp dropped. " *
                      "Increase the column's patch allocation."
            end
        end
        cp = cp.younger
    end
    return out
end

"""
    fates_set_filters!(inst; c, s)

Rebuild the HLM per-patch weights for FATES column `c` from the freshly-advanced
FATES patch areas, mirroring the Fortran host `wrap_update_hlmfates_dyn` +
`setFilters` sequence (clm_driver.F90:1151) that runs after `dynamics_driv`.

The downstream per-column averaging (`p2c_1d!`) weights a patch by
`pch.active[p] && pch.wtcol[p]`, so after a disturbance changes the FATES patch
population the HLM patch weights must be recomputed or the column aggregates would
still reflect the pre-disturbance single-patch layout. This sets, for column `c`:

  * the bare-ground HLM patch (`col.patchi[c]+0`): `is_bareground=true`,
    `wt_ed = max(0, 1 - sum(canopy_fraction_pa))`;
  * each vegetated HLM patch (`+ifp`): `is_veg=true`,
    `wt_ed = canopy_fraction_pa[ifp]`, and marks it `active`;
  * any reserved-but-now-unused veg slot beyond the live patch count: zeroed +
    `active=false` (a disturbance can MERGE patches → fewer live patches);
  * `pch.wtcol[p] = wt_ed[p]` for every patch in the column's range (the weight the
    column averaging reads).

This is the minimal `setFilters!`-equivalent: it propagates FATES patch-weight
changes into the weights the non-FATES column averaging consumes. A full CLM
BitVector-filter rebuild (reweight_wrapup) is NOT done — the FATES per-timestep
hooks index patches directly via `fates_veg_patches`, and `p2c_1d!` reads
`active`/`wtcol`, so the weights are the load-bearing quantity. Gated by the caller
behind `config.use_fates`.
"""
function fates_set_filters!(inst::CLMInstances; c::Int, s::Int)
    fates = inst.fates
    fates === nothing && return inst
    col = inst.column
    pch = inst.patch
    bc_out = fates.bc_out[s]
    site = fates.sites[s]

    pi = col.patchi[c]
    pf = col.patchf[c]

    # Number of vegetated patches in the site (youngest patch's patchno = the
    # total veg-patch count for the age-ordered list).
    npatch = site.youngest_patch === nothing ? 0 : site.youngest_patch.patchno

    # Reset the column's FATES weight flags over its whole patch range.
    for p in pi:pf
        if !isempty(pch.is_veg);        pch.is_veg[p] = false;        end
        if !isempty(pch.is_bareground); pch.is_bareground[p] = false; end
        if !isempty(pch.wt_ed);         pch.wt_ed[p] = 0.0;           end
    end

    # Bare-ground HLM patch (col.patchi[c]+0): weight = 1 - sum of canopy fractions.
    sum_canopy = 0.0
    for ifp in 1:npatch
        ifp <= length(bc_out.canopy_fraction_pa) || break
        sum_canopy += bc_out.canopy_fraction_pa[ifp]
    end
    bg = pi
    if !isempty(pch.is_bareground); pch.is_bareground[bg] = true; end
    wt_bg = max(0.0, 1.0 - sum_canopy)
    if !isempty(pch.wt_ed); pch.wt_ed[bg] = wt_bg; end
    pch.wtcol[bg] = wt_bg

    # Vegetated HLM patches (+ifp): weight = canopy_fraction_pa[ifp].
    for (ifp, p) in fates_veg_patches(site, c, col)
        wt = ifp <= length(bc_out.canopy_fraction_pa) ?
             bc_out.canopy_fraction_pa[ifp] : 0.0
        if !isempty(pch.is_veg); pch.is_veg[p] = true; end
        if !isempty(pch.wt_ed);  pch.wt_ed[p] = wt;    end
        pch.wtcol[p] = wt
        pch.active[p] = true
    end

    # Reserved-but-unused vegetated slots (a fuse/merge can leave fewer live
    # patches than were reserved): zero them and mark inactive so the column
    # averaging skips them.
    for p in (pi + npatch + 1):pf
        pch.wtcol[p] = 0.0
        pch.active[p] = false
    end

    return inst
end

"""
    _fates_first_fates_column(col) -> Int

Return the index of the first `col.is_fates` column, or 0 if there is none.
Convenience for the single-site MVP host glue (driver fire-weather derivation).
"""
function _fates_first_fates_column(col)
    for c in 1:length(col.is_fates)
        col.is_fates[c] && return c
    end
    return 0
end

"""
    fates_pack_bcin_radiation!(inst; s=1, c=1, p=1, coszen)

Fill the FATES radiation `bc_in` fields for site `s` from CLM column `c` /
vegetated patch `p` state:

  * `solad_parb` / `solai_parb` — direct/diffuse downwelling radiation per band
    (from atm2lnd `forc_solad_downscaled_col` / `forc_solai_grc`).
  * `coszen_pa` — cosine of the solar zenith angle (driver-supplied `coszen`).
  * `filter_vegzen_pa` — daylight flag (coszen > 0).
  * `albgr_dir_rb` / `albgr_dif_rb` — ground albedo per band (surfalb
    `albgrd_col` / `albgri_col`).
  * `fcansno_pa` — fraction of canopy covered in snow (0 stand-in; no CLM source
    on the carbon-only path).

The FATES radiation solver uses `num_swb` bands (vis, nir); CLM carries `NUMRAD`
(also 2) radiation bands in the same vis/nir order, so the band index maps directly.

The `ifp` keyword selects which vegetated-patch slot of the patch-dimensioned
`bc_in` arrays to fill (1-based, from the [`fates_veg_patches`](@ref) walk). The
ground albedo is site-broadband so it is (re)written each call. The forcing /
coszen are site/column-uniform, so each patch slot gets the same column value —
matching the Fortran host which broadcasts `forc_solad(g,:)` to every
`solad_parb(ifp,:)`. The driver hook loops `(ifp, p)` over the site's veg patches.
"""
function fates_pack_bcin_radiation!(inst::CLMInstances; s::Int = 1, c::Int = 1,
                                    p::Int = 1, ifp::Int = 1, coszen::Real)
    fates = inst.fates
    fates === nothing && return inst
    bc = fates.bc_in[s]

    a2l = inst.atm2lnd
    sa  = inst.surfalb

    g = inst.column.gridcell[c]

    # Radiation (patch, band) for this vegetated patch (ifp).
    for ib in 1:num_swb
        bc.solad_parb[ifp, ib] = a2l.forc_solad_downscaled_col[c, ib]
        bc.solai_parb[ifp, ib] = a2l.forc_solai_grc[g, ib]
    end

    bc.coszen_pa[ifp]        = coszen
    bc.filter_vegzen_pa[ifp] = coszen > 0.0
    bc.fcansno_pa[ifp]       = 0.0  # no snow-on-canopy source on the carbon-only path

    # Ground albedo (site broadband). Default to 0 if surfalb not yet computed
    # (the very first timestep before surface_albedo! runs).
    for ib in 1:num_swb
        ad = (!isempty(sa.albgrd_col) && isfinite(sa.albgrd_col[c, ib])) ?
             sa.albgrd_col[c, ib] : 0.0
        ai = (!isempty(sa.albgri_col) && isfinite(sa.albgri_col[c, ib])) ?
             sa.albgri_col[c, ib] : 0.0
        bc.albgr_dir_rb[ib] = ad
        bc.albgr_dif_rb[ib] = ai
    end

    return inst
end

"""
    fates_pack_bcin_btran!(inst; s=1, c=1, nlevsoil, smp_sl=nothing)

Fill the FATES BTRAN `bc_in` soil-hydrology fields for site `s` from CLM column
`c` state (reusing the existing CLM soil-suction `smp_l` path):

  * `tempk_sl` / `t_soisno_sl` — soil-layer temperature (K).
  * `h2o_liqvol_sl` — liquid volume in soil layer (m3/m3).
  * `watsat_sl` — saturated volumetric water (porosity).
  * `eff_porosity_sl` — effective porosity.
  * `smp_sl` — soil suction potential (mm, negative). Supplied via `smp_sl`
    (the host-computed CLM `smp_l_col` slice); if `nothing`, derived from the
    Clapp-Hornberger retention curve using watsat/sucsat/bsw + h2o_liqvol.

`max_rooting_depth_index_col` is left as set by `clm_fates_init!`.
"""
function fates_pack_bcin_btran!(inst::CLMInstances; s::Int = 1, c::Int = 1,
                                nlevsoil::Int,
                                smp_sl::Union{AbstractVector, Nothing} = nothing)
    fates = inst.fates
    fates === nothing && return inst
    bc = fates.bc_in[s]

    temp = inst.temperature
    wdb  = inst.water.waterdiagnosticbulk_inst
    ss   = inst.soilstate

    joff = varpar.nlevsno  # snow layers offset into t_soisno / h2osoi arrays

    for j in 1:nlevsoil
        tk = temp.t_soisno_col[c, joff + j]
        bc.tempk_sl[j]      = tk
        bc.t_soisno_sl[j]   = tk
        bc.h2o_liqvol_sl[j] = wdb.h2osoi_liqvol_col[c, joff + j]
        bc.watsat_sl[j]     = ss.watsat_col[c, j]
        bc.eff_porosity_sl[j] = ss.eff_porosity_col[c, j]
    end

    if smp_sl !== nothing
        for j in 1:nlevsoil
            bc.smp_sl[j] = smp_sl[j]
        end
    else
        # Derive soil matric suction from the Clapp-Hornberger retention curve,
        # matching the CLM soil-water-potential path: smp = -sucsat * s^(-bsw),
        # where s = liqvol/watsat clamped to [0.01, 1].
        for j in 1:nlevsoil
            wsat = ss.watsat_col[c, j]
            scl  = wsat > 0.0 ? clamp(wdb.h2osoi_liqvol_col[c, joff + j] / wsat, 0.01, 1.0) : 0.01
            sucsat = ss.sucsat_col[c, j]   # minimum soil suction (mm, positive)
            bsw    = ss.bsw_col[c, j]
            bc.smp_sl[j] = -sucsat * scl^(-bsw)
        end
    end

    return inst
end

"""
    fates_pack_bcin_photosynthesis!(inst; s=1, c=1, p=1, forc_pbot, forc_pco2,
                                    forc_po2, t_veg, tgcm, esat_tv, eair, rb,
                                    dayl_factor)

Fill the FATES photosynthesis `bc_in` fields for site `s` from CLM column `c` /
patch `p` canopy-solve locals:

  * `forc_pbot` — atmospheric pressure (Pa).
  * `cair_pa` / `oair_pa` — CO2 / O2 partial pressure (Pa).
  * `dayl_factor_pa` — daylength scaling factor (0-1).
  * `esat_tv_pa` / `eair_pa` — saturation vapor pressure at t_veg / canopy-air
    vapor pressure (Pa).
  * `rb_pa` — leaf boundary-layer resistance (s/m).
  * `t_veg_pa` — vegetation temperature (K).
  * `tgcm_pa` — air temperature at reference height (K).
  * `filter_photo_pa` — photosynthesis filter flag (2 = compute).

`forc_pco2`/`forc_po2` are the gridcell CO2/O2 partial pressures (Pa). The
soil-suction `smp_sl` and `coszen` were already packed by the btran/radiation
helpers; the radiation profiles come from `FatesSunShadeFracs`.
"""
function fates_pack_bcin_photosynthesis!(inst::CLMInstances; s::Int = 1, c::Int = 1,
                                         p::Int = 1, ifp::Int = 1,
                                         forc_pbot::Real, forc_pco2::Real,
                                         forc_po2::Real, t_veg::Real, tgcm::Real,
                                         esat_tv::Real, eair::Real, rb::Real,
                                         dayl_factor::Real)
    fates = inst.fates
    fates === nothing && return inst
    bc = fates.bc_in[s]

    bc.forc_pbot            = forc_pbot   # site scalar
    bc.cair_pa[ifp]         = forc_pco2
    bc.oair_pa[ifp]         = forc_po2
    bc.dayl_factor_pa[ifp]  = dayl_factor
    bc.esat_tv_pa[ifp]      = esat_tv
    bc.eair_pa[ifp]         = eair
    bc.rb_pa[ifp]           = rb
    bc.t_veg_pa[ifp]        = t_veg
    bc.tgcm_pa[ifp]         = tgcm
    bc.filter_photo_pa[ifp] = 2  # 2 = "compute photosynthesis" branch

    return inst
end

"""
    fates_pack_bcin_daily!(inst; s=1, c=1, nlevsoil, nlevdecomp, use_fates_bgc=false)

Fill the FATES *daily* demographic `bc_in` soil state for site `s` from CLM column
`c`, mirroring the carbon-only, no-fire, no-hydro, no-LUH part of the Fortran host
`dynamics_driv` (clmfates_interfaceMod.F90 Part I):

  * `tempk_sl` / `t_soisno_sl` — soil-layer temperature (K).
  * `h2o_liqvol_sl` — liquid volume in soil layer (m3/m3).
  * `watsat_sl` / `eff_porosity_sl` — porosity / effective porosity.
  * `smp_sl` — soil suction potential (mm) from the Clapp-Hornberger curve.
  * `w_scalar_sisl` / `t_scalar_sisl` — soil-decomposition moisture/temperature
    limitation scalars consumed by litter fragmentation (EDPhysiologyMod). Only
    packed when `use_fates_bgc` (Fortran gates these on the soilbiogeochem carbon
    flux being live); otherwise zeroed, matching the Fortran else-branch.
  * `max_rooting_depth_index_col` — left as set by `clm_fates_init!`.

Fire drivers (lightning/pop/RH24/wind24/precip24), SP-LAI streams, harvest/LUH
arrays are intentionally NOT packed — those configs are off on this MVP path.

This reuses the same soil-state slice as `fates_pack_bcin_btran!`; calling both in
one step is harmless (idempotent on the soil fields). The decomposition scalars
are the daily-step-specific additions.
"""
function fates_pack_bcin_daily!(inst::CLMInstances; s::Int = 1, c::Int = 1,
                                nlevsoil::Int, nlevdecomp::Int,
                                use_fates_bgc::Bool = false,
                                fire_weather::Union{NamedTuple,Nothing} = nothing)
    fates = inst.fates
    fates === nothing && return inst
    bc = fates.bc_in[s]

    temp = inst.temperature
    wdb  = inst.water.waterdiagnosticbulk_inst
    ss   = inst.soilstate

    joff = varpar.nlevsno  # snow layers offset into t_soisno / h2osoi arrays

    for j in 1:nlevsoil
        tk = temp.t_soisno_col[c, joff + j]
        bc.tempk_sl[j]        = tk
        bc.t_soisno_sl[j]     = tk
        bc.h2o_liqvol_sl[j]   = wdb.h2osoi_liqvol_col[c, joff + j]
        bc.watsat_sl[j]       = ss.watsat_col[c, j]
        bc.eff_porosity_sl[j] = ss.eff_porosity_col[c, j]

        wsat = ss.watsat_col[c, j]
        scl  = wsat > 0.0 ? clamp(wdb.h2osoi_liqvol_col[c, joff + j] / wsat, 0.01, 1.0) : 0.01
        bc.smp_sl[j] = -ss.sucsat_col[c, j] * scl^(-ss.bsw_col[c, j])
    end

    # Soil-decomposition limitation scalars: live only under use_fates_bgc, else 0
    # (matches the Fortran w_scalar/t_scalar gate in dynamics_driv). Clamp the loop
    # to the FATES bc array length (the cold start sizes these to FATES nlevdecomp,
    # which may differ from the CLM-side nlevdecomp passed by the driver).
    ndl = min(nlevdecomp, length(bc.w_scalar_sisl))
    if use_fates_bgc
        cf = inst.soilbiogeochem_carbonflux
        for j in 1:ndl
            wsv = (!isempty(cf.w_scalar_col) && isfinite(cf.w_scalar_col[c, j])) ?
                  cf.w_scalar_col[c, j] : 0.0
            tsv = (!isempty(cf.t_scalar_col) && isfinite(cf.t_scalar_col[c, j])) ?
                  cf.t_scalar_col[c, j] : 0.0
            bc.w_scalar_sisl[j] = wsv
            bc.t_scalar_sisl[j] = tsv
        end
    else
        for j in 1:ndl
            bc.w_scalar_sisl[j] = 0.0
            bc.t_scalar_sisl[j] = 0.0
        end
    end

    # SPITFIRE fire-weather drivers (per-patch 24-hour running means the fire model
    # reads in `update_fire_weather!` / `area_burnt_intensity!`). Packed only when
    # SPITFIRE is on (`fire_weather` supplied); the carbon-only default leaves these
    # zero (allocate_bcin already zeroed them), matching the no_fire path which never
    # reads them. The fire model indexes these by patchno (1-based vegetated patch),
    # so fill the whole allocated span with the representative site values; the
    # bareground patch (patchno 0) is skipped inside the fire model.
    if fire_weather !== nothing
        fw = fire_weather
        np_fire = length(bc.precip24_pa)
        @inbounds for ip in 1:np_fire
            bc.precip24_pa[ip]   = fw.precip24      # [mm/s]
            bc.relhumid24_pa[ip] = fw.relhumid24    # [%]
            bc.wind24_pa[ip]     = fw.wind24        # [m/s]
            bc.lightning24[ip]   = get(fw, :lightning24, 0.0)  # [#/km2/day]
            bc.pop_density[ip]   = get(fw, :pop_density, 0.0)  # [#/km2]
        end
    end

    return inst
end

"""
    fates_daily_dynamics_step!(inst; nlevsoil, nlevdecomp, use_fates_bgc=false)

Run one FATES *daily demographic step* for every FATES site attached to `inst`,
mirroring the Fortran host `dynamics_driv`:

  1. pack the daily `bc_in` (`fates_pack_bcin_daily!`) for each site's CLM column,
  2. `ed_ecosystem_dynamics(site, bc_in, bc_out)` — the demographic advance
     (phenology -> growth/allocation -> mortality -> recruitment -> cohort/patch
     dynamics, interleaved with `TotalBalanceCheck` mass-conservation audits),
  3. `ed_update_site(site, bc_in, bc_out, false)` — recompute site diagnostics +
     canopy structure (runs `TotalBalanceCheck` itself),
  4. `TotalBalanceCheck(site, -1)` — final mass-conservation audit (throws on
     imbalance > 1e-5),
  5. unpack the canopy structure (`fates_unpack_bcout_canopy_structure!`) back into
     CLM canopystate (elai/esai/htop/hbot/z0m/displa/dleaf) for EACH vegetated patch
     of the site (the `ifp` walk — a disturbance can split a patch into several),
  6. rebuild the HLM per-patch weights for the column
     (`fates_set_filters!`, mirroring the Fortran host `setFilters` after
     `dynamics_driv`) so the post-disturbance patch-weight changes propagate into
     the per-column averaging.

The site<->column map mirrors the W3/W4 hooks: FATES site `s` maps onto the `s`-th
`col.is_fates` column; each vegetated patch maps to HLM patch
`p = ifp + col.patchi[c]` (bare-ground at `+0`).

When SPITFIRE is active (`hlm_spitfire_mode > hlm_sf_nofire_def`), pass `fire_weather`
(a `NamedTuple` with `precip24` [mm/s], `relhumid24` [%], `wind24` [m/s], optionally
`lightning24`/`pop_density`) — the daily pack fills the per-patch fire-weather bc_in
the fire model reads. Omit it on the carbon-only / no_fire path (the fields stay zero).

Returns `inst`. Gated by the caller behind `config.use_fates`.
"""
function fates_daily_dynamics_step!(inst::CLMInstances; nlevsoil::Int,
                                    nlevdecomp::Int, use_fates_bgc::Bool = false,
                                    fire_weather::Union{NamedTuple,Nothing} = nothing)
    fates = inst.fates
    fates === nothing && return inst
    col = inst.column

    # 0. CNP: parse the nutrient-acquisition boundary conditions out to the cohorts
    #    BEFORE the demographic step (sets cohort daily_n/p_demand + daily_*_uptake
    #    from the prescribed-uptake params). Mirrors the Fortran host order
    #    (clmfates_interfaceMod: UnPackNutrientAquisitionBCs precedes the per-site
    #    ed_ecosystem_dynamics loop). A no-op (zeroes the uptake bc) in carbon-only
    #    mode, so it is safe to call unconditionally.
    UnPackNutrientAquisitionBCs(fates.sites, fates.bc_in)

    # Map each FATES site index s -> its CLM column c (the s-th is_fates column).
    site_col = zeros(Int, fates.nsites)

    s = 0
    for c in 1:length(col.is_fates)
        col.is_fates[c] || continue
        s += 1
        s <= fates.nsites || break
        site_col[s] = c

        # 1. pack the daily soil/decomp bc_in (+ SPITFIRE fire-weather when on).
        fates_pack_bcin_daily!(inst; s=s, c=c, nlevsoil=nlevsoil,
                               nlevdecomp=nlevdecomp, use_fates_bgc=use_fates_bgc,
                               fire_weather=fire_weather)

        site = fates.sites[s]
        bc_in  = fates.bc_in[s]
        bc_out = fates.bc_out[s]

        # 2-3. the demographic step + site diagnostics update.
        ed_ecosystem_dynamics(site, bc_in, bc_out)
        ed_update_site(site, bc_in, bc_out, false)

        # 4. final mass-conservation audit (throws on imbalance).
        TotalBalanceCheck(site, -1)
    end

    # 5. summarize the per-patch canopy totals (recomputes cohort crown areas
    #    `c_area` + patch `total_canopy_area` from allometry), then pack the FATES
    #    canopy state into bc_out for ALL sites — fills the per-vegetated-patch *_pa
    #    arrays (elai/esai/htop/.../canopy_fraction_pa) for every (possibly newly
    #    split) patch. canopy_summarization! is required before update_hlm_dynamics!:
    #    the latter derives canopy_fraction_pa/elai_pa from `total_canopy_area`, which
    #    is otherwise stale/NaN if no radiation pass has run this step (the real
    #    driver runs radiation first, but a daily step driven in isolation has not).
    #    Mirrors the Fortran host (canopy summarization precedes update_hlm_dynamics
    #    in wrap_update_hlmfates_dyn after dynamics). `fcolumn` (site->column map) is
    #    unused by the body, which indexes only by site.
    canopy_summarization!(fates.nsites, fates.sites, fates.bc_in)
    update_hlm_dynamics!(fates.nsites, fates.sites, site_col, fates.bc_out)

    # 6. per FATES column: unpack the canopy structure into CLM for each vegetated
    #    patch (the ifp walk — a disturbance can split one patch into several), then
    #    rebuild the HLM per-patch weights (setFilters-equivalent) from the new
    #    canopy_fraction_pa so the per-column averaging sees the right weights.
    for s in 1:fates.nsites
        c = site_col[s]
        c == 0 && continue
        for (ifp, p) in fates_veg_patches(fates.sites[s], c, col)
            fates_unpack_bcout_canopy_structure!(inst; s=s, c=c, p=p, ifp=ifp)
        end
        fates_set_filters!(inst; c=c, s=s)
    end

    # 7. fill the daily ("dynamics") history buffers from the freshly-advanced
    #    site/patch/cohort state. Mirrors the Fortran host's `update_history_dyn`
    #    call after `dynamics_driv`. Write-only diagnostics; gated on the history
    #    interface having been built (clm_fates_init!) — a no-op otherwise.
    if fates.hist !== nothing
        update_history_dyn!(fates.hist, 1, fates.nsites, fates.sites, fates.bc_in)
    end

    return inst
end

"""
    fates_hifrq_history_step!(inst; dt_tstep)

Fill the high-frequency ("hifrq", per-timestep) FATES history buffers from the
current per-step cohort/patch carbon-flux + radiation state. Mirrors the Fortran
host's `update_history_hifrq` call in the per-timestep path. Write-only
diagnostics; a no-op when the history interface has not been built. `dt_tstep`
is the model timestep length [s]. Gated by the caller behind `config.use_fates`.
"""
function fates_hifrq_history_step!(inst::CLMInstances; dt_tstep::Real)
    fates = inst.fates
    fates === nothing && return inst
    fates.hist === nothing && return inst
    update_history_hifrq!(fates.hist, 1, fates.nsites, fates.sites,
                          fates.bc_in, fates.bc_out, dt_tstep)
    # plant-hydraulics history (no-op unless hlm_use_planthydro is on + si_hydr built)
    update_history_hydraulics!(fates.hist, 1, fates.nsites, fates.sites,
                               fates.bc_in, dt_tstep)
    return inst
end

"""
    fates_pack_bcin_hydro!(inst; s=1, c=1, p=1, nlevsoil)

Fill the FATES plant-hydraulics `bc_in` soil-rhizosphere fields for site `s` from
CLM column `c` / vegetated patch `p`, packing exactly the inputs the
`hydraulics_drive` (FillDrainRhizShells + Hydraulics_BC) solve reads:

  * `watsat_sisl` / `eff_porosity_sl` — saturated / effective porosity.
  * `sucsat_sisl` — minimum soil suction (mm).
  * `bsw_sisl` — Clapp-Hornberger "b".
  * `hksat_sisl` — saturated hydraulic conductivity (mm/s).
  * `h2o_liq_sisl` — liquid water mass in the soil layer (kg/m2).
  * `smpmin_si` — host minimum soil matric potential (mm).
  * `qflx_transp_pa[ifp]` — the HLM canopy-solver transpiration demand for the
    vegetated patch (mm H2O/s, + into root). ifp = 1 for the single veg patch.

`zi_sisl`/`dz_sisl` are set once at cold-start and are not re-packed here.
Gated by the caller behind `config.use_fates && hlm_use_planthydro==itrue`.
"""
function fates_pack_bcin_hydro!(inst::CLMInstances; s::Int = 1, c::Int = 1,
                                p::Int = 1, nlevsoil::Int)
    fates = inst.fates
    fates === nothing && return inst
    bc = fates.bc_in[s]

    ss  = inst.soilstate
    wsb = inst.water.waterstatebulk_inst
    wfb = inst.water.waterfluxbulk_inst
    col = inst.column

    joff = varpar.nlevsno  # snow-layer offset into h2osoi_liq_col

    for j in 1:nlevsoil
        bc.watsat_sisl[j]     = ss.watsat_col[c, j]
        bc.eff_porosity_sl[j] = ss.eff_porosity_col[c, j]
        bc.sucsat_sisl[j]     = ss.sucsat_col[c, j]
        bc.bsw_sisl[j]        = ss.bsw_col[c, j]
        bc.hksat_sisl[j]      = ss.hksat_col[c, j]
        bc.h2o_liq_sisl[j]    = wsb.ws.h2osoi_liq_col[c, joff + j]
    end

    bc.smpmin_si = (!isempty(ss.smpmin_col) && isfinite(ss.smpmin_col[c])) ?
                   ss.smpmin_col[c] : -1.0e8

    # Patch transpiration demand per vegetated patch (ifp walk → HLM patch p).
    for (ifp, pp) in fates_veg_patches(fates.sites[s], c, col)
        bc.qflx_transp_pa[ifp] = wfb.wf.qflx_tran_veg_patch[pp]
    end

    return inst
end

"""
    fates_hydraulics_step!(inst; nlevsoil)

Run one FATES plant-hydraulics transpiration/uptake step for every FATES site
attached to `inst`, mirroring the Fortran host's `hydraulics_drive` call from the
canopy-flux/photosynthesis path:

  1. pack the soil-rhizosphere + transpiration `bc_in` (`fates_pack_bcin_hydro!`)
     for each site's CLM column / vegetated patch,
  2. `hydraulics_drive(nsites, sites, bc_in, bc_out, dtime)` — reconcile the
     rhizosphere shells with the host soil column then run `Hydraulics_BC`, which
     populates `co_hydr.ftc_ag/ftc_troot/ftc_aroot`, `co_hydr.btran`, and
     `psi_ag[1]` (leaf_psi) — the fields the mortality + photosynthesis branches
     consume — and fills `bc_out.qflx_soil2root_sisl` / `plant_stored_h2o_si`.

The site<->column<->veg-patch map mirrors the W3/W4 hooks. Returns `inst`. Gated
by the caller behind `config.use_fates && hlm_use_planthydro[]==itrue`; a no-op
when no site has a built `si_hydr`.
"""
function fates_hydraulics_step!(inst::CLMInstances; nlevsoil::Int, dtime::Real)
    fates = inst.fates
    fates === nothing && return inst
    hlm_use_planthydro[] == itrue || return inst
    col = inst.column

    s = 0
    for c in 1:length(col.is_fates)
        col.is_fates[c] || continue
        s += 1
        s <= fates.nsites || break
        # A site without a built hydraulics object (planthydro just toggled, or no
        # cold-start hydro init) has nothing to solve — skip it.
        fates.sites[s].si_hydr === nothing && continue
        # fates_pack_bcin_hydro! walks all the site's vegetated patches internally
        # (ifp → HLM patch p) for the per-patch transpiration demand.
        fates_pack_bcin_hydro!(inst; s=s, c=c, nlevsoil=nlevsoil)
    end

    # Only drive when at least one site is hydro-enabled.
    any(site -> site.si_hydr !== nothing, fates.sites) || return inst
    hydraulics_drive(fates.nsites, fates.sites, fates.bc_in, fates.bc_out, dtime)
    return inst
end

"""
    fates_unpack_bcout_sunfrac!(inst; s=1, c=1, p=1)

Write the FATES sun/shade `bc_out` back to CLM canopystate for column `c` /
patch `p`: `fsun_pa -> fsun_patch`, `laisun_pa -> laisun_patch`,
`laisha_pa -> laisha_patch`.
"""
function fates_unpack_bcout_sunfrac!(inst::CLMInstances; s::Int = 1, c::Int = 1,
                                     p::Int = 1, ifp::Int = 1)
    fates = inst.fates
    fates === nothing && return inst
    bc = fates.bc_out[s]
    cs = inst.canopystate

    cs.fsun_patch[p]   = bc.fsun_pa[ifp]
    cs.laisun_patch[p] = bc.laisun_pa[ifp]
    cs.laisha_patch[p] = bc.laisha_pa[ifp]

    return inst
end

"""
    fates_unpack_bcout_canopy_radiation!(inst; s=1, c=1, p=1)

Write the FATES normalized-canopy-radiation `bc_out` back to CLM surfalb for
patch `p`: `albd_parb -> albd_patch`, `albi_parb -> albi_patch`,
`fabd/fabi_parb -> fabd/fabi_patch`, `ftdd/ftid/ftii_parb -> ftdd/ftid/ftii_patch`.
"""
function fates_unpack_bcout_canopy_radiation!(inst::CLMInstances; s::Int = 1,
                                              c::Int = 1, p::Int = 1, ifp::Int = 1)
    fates = inst.fates
    fates === nothing && return inst
    bc = fates.bc_out[s]
    sa = inst.surfalb

    for ib in 1:num_swb
        sa.albd_patch[p, ib] = bc.albd_parb[ifp, ib]
        sa.albi_patch[p, ib] = bc.albi_parb[ifp, ib]
        sa.fabd_patch[p, ib] = bc.fabd_parb[ifp, ib]
        sa.fabi_patch[p, ib] = bc.fabi_parb[ifp, ib]
        sa.ftdd_patch[p, ib] = bc.ftdd_parb[ifp, ib]
        sa.ftid_patch[p, ib] = bc.ftid_parb[ifp, ib]
        sa.ftii_patch[p, ib] = bc.ftii_parb[ifp, ib]
    end

    return inst
end

"""
    fates_unpack_bcout_btran!(inst; s=1, c=1, p=1, nlevsoil)

Write the FATES BTRAN `bc_out` back to CLM: `btran_pa -> energyflux.btran_patch`,
`rootr_pasl -> soilstate.rootr_patch` (per soil layer).
"""
function fates_unpack_bcout_btran!(inst::CLMInstances; s::Int = 1, c::Int = 1,
                                   p::Int = 1, ifp::Int = 1, nlevsoil::Int)
    fates = inst.fates
    fates === nothing && return inst
    bc = fates.bc_out[s]
    ef = inst.energyflux
    ss = inst.soilstate

    ef.btran_patch[p] = bc.btran_pa[ifp]
    for j in 1:nlevsoil
        ss.rootr_patch[p, j] = bc.rootr_pasl[ifp, j]
    end

    return inst
end

"""
    fates_unpack_bcout_photosynthesis!(inst; s=1, c=1, p=1)

Write the FATES photosynthesis `bc_out` resistances back to CLM photosyns for
patch `p`: `rssun_pa -> rssun_patch`, `rssha_pa -> rssha_patch`.
"""
function fates_unpack_bcout_photosynthesis!(inst::CLMInstances; s::Int = 1,
                                            c::Int = 1, p::Int = 1, ifp::Int = 1)
    fates = inst.fates
    fates === nothing && return inst
    bc = fates.bc_out[s]
    ps = inst.photosyns

    ps.rssun_patch[p] = bc.rssun_pa[ifp]
    ps.rssha_patch[p] = bc.rssha_pa[ifp]

    return inst
end

"""
    fates_unpack_bcout_canopy_structure!(inst; s=1, c=1, p=1, ifp=1)

Write the FATES canopy-structure `bc_out` back to CLM canopystate for patch `p`
(vegetated-patch slot `ifp`): `elai/esai_pa`, `htop/hbot_pa`, `z0m/displa/dleaf_pa`.
The daily step runs `canopy_summarization!` + `update_hlm_dynamics!` before this, so
the `bc_out` `*_pa[ifp]` slots are finite.
"""
function fates_unpack_bcout_canopy_structure!(inst::CLMInstances; s::Int = 1,
                                              c::Int = 1, p::Int = 1, ifp::Int = 1)
    fates = inst.fates
    fates === nothing && return inst
    bc = fates.bc_out[s]
    cs = inst.canopystate

    cs.elai_patch[p]   = bc.elai_pa[ifp]
    cs.esai_patch[p]   = bc.esai_pa[ifp]
    cs.htop_patch[p]   = bc.htop_pa[ifp]
    cs.hbot_patch[p]   = bc.hbot_pa[ifp]
    cs.z0m_patch[p]    = bc.z0m_pa[ifp]
    cs.displa_patch[p] = bc.displa_pa[ifp]
    cs.dleaf_patch[p]  = bc.dleaf_pa[ifp]

    return inst
end
