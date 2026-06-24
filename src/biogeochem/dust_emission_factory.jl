# ==========================================================================
# Ported from: src/biogeochem/DustEmisFactory.F90
# Dust-emission scheme selector + config.
#
# CTSM's create_dust_emissions() picks a dust_emis_base subclass from the
# `dust_emis_method` namelist item (shr_dust_emis_mod's is_dust_emis_zender /
# is_dust_emis_leung). Here the scheme is a string carried on a small
# `DustEmisConfig`; `dust_emission!` dispatches on it. Disabled by default so
# the default driver path is byte-identical (no emission call).
# ==========================================================================

# Recognized dust-emission method strings (match CTSM shr_dust_emis_mod).
const DUST_EMIS_ZENDER2003 = "Zender_2003"
const DUST_EMIS_LEUNG2023  = "Leung_2023"

"""
    DustEmisConfig

Configuration for the dust-emission scheme. `active=false` (the default) means
no dust-emission mobilization is run (the historical CLM.jl behavior — only dry
deposition). When `active=true`, `method` selects the scheme:
`"Zender_2003"` (default CLM5) or `"Leung_2023"`.

Mirrors the role of `dust_emis_method` + `DustEmisFactory.F90`.
"""
Base.@kwdef struct DustEmisConfig
    active::Bool = false
    method::String = DUST_EMIS_ZENDER2003
end

is_dust_emis_zender(cfg::DustEmisConfig) = cfg.method == DUST_EMIS_ZENDER2003
is_dust_emis_leung(cfg::DustEmisConfig)  = cfg.method == DUST_EMIS_LEUNG2023

"""
    dust_emission!(cfg::DustEmisConfig, dust, nolakep_mask, pch, lun, bounds_p, nl,
                   forc_rho, soilstate, canopystate, water, waterdiag, fv_inst)

Factory dispatcher (equivalent to `create_dust_emissions` + the virtual
`DustEmission` call). When `cfg.active` is false this is a no-op. Otherwise it
runs the scheme selected by `cfg.method`, filling `dust.flx_mss_vrt_dst_patch`.

Arguments are the high-level driver objects:
- `dust`        :: `DustEmisBaseData`
- `nolakep_mask`:: non-lake patch BitVector
- `pch`/`lun`   :: `PatchData` / `LandunitData` (g/l/c/p hierarchy + types)
- `bounds_p`    :: patch index range; `nl` :: number of landunits
- `forc_rho`    :: downscaled air density per column
- `soilstate`   :: `SoilStateData` (gwc_thr_col, mss_frc_cly_vld_col, watsat_col)
- `canopystate` :: `CanopyStateData` (tlai_patch, tsai_patch)
- `ws`          :: water-state sub-struct with `h2osoi_vol/liq/ice_col`
                   (i.e. `inst.water.waterstatebulk_inst.ws`)
- `wd`          :: water-diagnostic struct with `frac_sno_col`
                   (i.e. `inst.water.waterdiagnosticbulk_inst`)
- `fv_inst`     :: `FrictionVelocityData` (fv_patch, u10_patch, obu_patch)

Ported from `DustEmisFactory.F90`.
"""
function dust_emission!(cfg::DustEmisConfig, dust::DustEmisBaseData,
                        nolakep_mask::AbstractVector{Bool},
                        pch, lun, bounds_p::UnitRange{Int}, nl::Int,
                        forc_rho::AbstractVector{<:Real},
                        soilstate, canopystate, ws, wd, fv_inst)
    cfg.active || return nothing

    if is_dust_emis_zender(cfg)
        dust_emission_zender2003!(
            dust, nolakep_mask, pch.active, pch.column, pch.landunit,
            pch.wtlunit, lun.itype, bounds_p, nl,
            forc_rho, soilstate.gwc_thr_col, soilstate.mss_frc_cly_vld_col,
            soilstate.watsat_col, canopystate.tlai_patch, canopystate.tsai_patch,
            wd.frac_sno_col, ws.h2osoi_vol_col, ws.h2osoi_liq_col,
            ws.h2osoi_ice_col, fv_inst.fv_patch, fv_inst.u10_patch)
    elseif is_dust_emis_leung(cfg)
        dust_emission_leung2023!(
            dust, nolakep_mask, pch.active, pch.column, pch.landunit, pch.itype,
            pch.wtlunit, lun.itype, bounds_p, nl,
            forc_rho, soilstate.gwc_thr_col, soilstate.mss_frc_cly_vld_col,
            soilstate.watsat_col, canopystate.tlai_patch, canopystate.tsai_patch,
            wd.frac_sno_col, ws.h2osoi_vol_col, ws.h2osoi_liq_col,
            ws.h2osoi_ice_col, fv_inst.fv_patch, fv_inst.obu_patch)
    else
        error("Unrecognized dust_emis_method: $(cfg.method)")
    end

    return nothing
end
