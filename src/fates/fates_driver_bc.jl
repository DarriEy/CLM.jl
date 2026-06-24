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
"""
function fates_pack_bcin_radiation!(inst::CLMInstances; s::Int = 1, c::Int = 1,
                                    p::Int = 1, coszen::Real)
    fates = inst.fates
    fates === nothing && return inst
    bc = fates.bc_in[s]

    a2l = inst.atm2lnd
    sa  = inst.surfalb

    g = inst.column.gridcell[c]

    # Radiation (patch, band). ifp = 1 (single vegetated patch).
    for ib in 1:num_swb
        bc.solad_parb[1, ib] = a2l.forc_solad_downscaled_col[c, ib]
        bc.solai_parb[1, ib] = a2l.forc_solai_grc[g, ib]
    end

    bc.coszen_pa[1]        = coszen
    bc.filter_vegzen_pa[1] = coszen > 0.0
    bc.fcansno_pa[1]       = 0.0  # no snow-on-canopy source on the carbon-only path

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
        bc.h2o_liqvol_sl[j] = wdb.h2osoi_liqvol_col[c, j]
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
            scl  = wsat > 0.0 ? clamp(wdb.h2osoi_liqvol_col[c, j] / wsat, 0.01, 1.0) : 0.01
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
                                         p::Int = 1, forc_pbot::Real, forc_pco2::Real,
                                         forc_po2::Real, t_veg::Real, tgcm::Real,
                                         esat_tv::Real, eair::Real, rb::Real,
                                         dayl_factor::Real)
    fates = inst.fates
    fates === nothing && return inst
    bc = fates.bc_in[s]

    bc.forc_pbot          = forc_pbot
    bc.cair_pa[1]         = forc_pco2
    bc.oair_pa[1]         = forc_po2
    bc.dayl_factor_pa[1]  = dayl_factor
    bc.esat_tv_pa[1]      = esat_tv
    bc.eair_pa[1]         = eair
    bc.rb_pa[1]           = rb
    bc.t_veg_pa[1]        = t_veg
    bc.tgcm_pa[1]         = tgcm
    bc.filter_photo_pa[1] = 2  # 2 = "compute photosynthesis" branch

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
                                use_fates_bgc::Bool = false)
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
        bc.h2o_liqvol_sl[j]   = wdb.h2osoi_liqvol_col[c, j]
        bc.watsat_sl[j]       = ss.watsat_col[c, j]
        bc.eff_porosity_sl[j] = ss.eff_porosity_col[c, j]

        wsat = ss.watsat_col[c, j]
        scl  = wsat > 0.0 ? clamp(wdb.h2osoi_liqvol_col[c, j] / wsat, 0.01, 1.0) : 0.01
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
     CLM canopystate (elai/esai/htop/hbot/z0m/displa/dleaf) for each site's veg patch.

The site<->column<->veg-patch map mirrors the W3/W4 hooks: FATES site `s` maps onto
the `s`-th `col.is_fates` column, vegetated patch = `col.patchi[c] + 1`.

Returns `inst`. Gated by the caller behind `config.use_fates`.
"""
function fates_daily_dynamics_step!(inst::CLMInstances; nlevsoil::Int,
                                    nlevdecomp::Int, use_fates_bgc::Bool = false)
    fates = inst.fates
    fates === nothing && return inst
    col = inst.column

    s = 0
    for c in 1:length(col.is_fates)
        col.is_fates[c] || continue
        s += 1
        s <= fates.nsites || break
        p = col.patchi[c] + 1   # vegetated patch (bare-ground at +0)

        # 1. pack the daily soil/decomp bc_in.
        fates_pack_bcin_daily!(inst; s=s, c=c, nlevsoil=nlevsoil,
                               nlevdecomp=nlevdecomp, use_fates_bgc=use_fates_bgc)

        site = fates.sites[s]
        bc_in  = fates.bc_in[s]
        bc_out = fates.bc_out[s]

        # 2-3. the demographic step + site diagnostics update.
        ed_ecosystem_dynamics(site, bc_in, bc_out)
        ed_update_site(site, bc_in, bc_out, false)

        # 4. final mass-conservation audit (throws on imbalance).
        TotalBalanceCheck(site, -1)

        # 5. unpack canopy structure back into CLM.
        fates_unpack_bcout_canopy_structure!(inst; s=s, c=c, p=p)
    end

    # 6. fill the daily ("dynamics") history buffers from the freshly-advanced
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
    return inst
end

"""
    fates_unpack_bcout_sunfrac!(inst; s=1, c=1, p=1)

Write the FATES sun/shade `bc_out` back to CLM canopystate for column `c` /
patch `p`: `fsun_pa -> fsun_patch`, `laisun_pa -> laisun_patch`,
`laisha_pa -> laisha_patch`.
"""
function fates_unpack_bcout_sunfrac!(inst::CLMInstances; s::Int = 1, c::Int = 1, p::Int = 1)
    fates = inst.fates
    fates === nothing && return inst
    bc = fates.bc_out[s]
    cs = inst.canopystate

    cs.fsun_patch[p]   = bc.fsun_pa[1]
    cs.laisun_patch[p] = bc.laisun_pa[1]
    cs.laisha_patch[p] = bc.laisha_pa[1]

    return inst
end

"""
    fates_unpack_bcout_canopy_radiation!(inst; s=1, c=1, p=1)

Write the FATES normalized-canopy-radiation `bc_out` back to CLM surfalb for
patch `p`: `albd_parb -> albd_patch`, `albi_parb -> albi_patch`,
`fabd/fabi_parb -> fabd/fabi_patch`, `ftdd/ftid/ftii_parb -> ftdd/ftid/ftii_patch`.
"""
function fates_unpack_bcout_canopy_radiation!(inst::CLMInstances; s::Int = 1,
                                              c::Int = 1, p::Int = 1)
    fates = inst.fates
    fates === nothing && return inst
    bc = fates.bc_out[s]
    sa = inst.surfalb

    for ib in 1:num_swb
        sa.albd_patch[p, ib] = bc.albd_parb[1, ib]
        sa.albi_patch[p, ib] = bc.albi_parb[1, ib]
        sa.fabd_patch[p, ib] = bc.fabd_parb[1, ib]
        sa.fabi_patch[p, ib] = bc.fabi_parb[1, ib]
        sa.ftdd_patch[p, ib] = bc.ftdd_parb[1, ib]
        sa.ftid_patch[p, ib] = bc.ftid_parb[1, ib]
        sa.ftii_patch[p, ib] = bc.ftii_parb[1, ib]
    end

    return inst
end

"""
    fates_unpack_bcout_btran!(inst; s=1, c=1, p=1, nlevsoil)

Write the FATES BTRAN `bc_out` back to CLM: `btran_pa -> energyflux.btran_patch`,
`rootr_pasl -> soilstate.rootr_patch` (per soil layer).
"""
function fates_unpack_bcout_btran!(inst::CLMInstances; s::Int = 1, c::Int = 1,
                                   p::Int = 1, nlevsoil::Int)
    fates = inst.fates
    fates === nothing && return inst
    bc = fates.bc_out[s]
    ef = inst.energyflux
    ss = inst.soilstate

    ef.btran_patch[p] = bc.btran_pa[1]
    for j in 1:nlevsoil
        ss.rootr_patch[p, j] = bc.rootr_pasl[1, j]
    end

    return inst
end

"""
    fates_unpack_bcout_photosynthesis!(inst; s=1, c=1, p=1)

Write the FATES photosynthesis `bc_out` resistances back to CLM photosyns for
patch `p`: `rssun_pa -> rssun_patch`, `rssha_pa -> rssha_patch`.
"""
function fates_unpack_bcout_photosynthesis!(inst::CLMInstances; s::Int = 1,
                                            c::Int = 1, p::Int = 1)
    fates = inst.fates
    fates === nothing && return inst
    bc = fates.bc_out[s]
    ps = inst.photosyns

    ps.rssun_patch[p] = bc.rssun_pa[1]
    ps.rssha_patch[p] = bc.rssha_pa[1]

    return inst
end

"""
    fates_unpack_bcout_canopy_structure!(inst; s=1, c=1, p=1)

Write the FATES canopy-structure `bc_out` back to CLM canopystate for patch `p`:
`elai/esai_pa`, `htop/hbot_pa`, `z0m/displa/dleaf_pa`.
"""
function fates_unpack_bcout_canopy_structure!(inst::CLMInstances; s::Int = 1,
                                              c::Int = 1, p::Int = 1)
    fates = inst.fates
    fates === nothing && return inst
    bc = fates.bc_out[s]
    cs = inst.canopystate

    cs.elai_patch[p]   = bc.elai_pa[1]
    cs.esai_patch[p]   = bc.esai_pa[1]
    cs.htop_patch[p]   = bc.htop_pa[1]
    cs.hbot_patch[p]   = bc.hbot_pa[1]
    cs.z0m_patch[p]    = bc.z0m_pa[1]
    cs.displa_patch[p] = bc.displa_pa[1]
    cs.dleaf_patch[p]  = bc.dleaf_pa[1]

    return inst
end
