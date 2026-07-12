# ==========================================================================
# Ported from: src/main/readParamsMod.F90
# Read parameters from clm5_params.nc — populates pftcon and scalar params
# ==========================================================================

"""
    readParameters!(paramfile::String)

Master parameter reader. Opens clm5_params.nc, reads PFT parameters into the
global `pftcon`, and reads scalar parameters used by various modules.
"""
function readParameters!(paramfile::String)
    ds = NCDataset(paramfile, "r")
    try
        # Read PFT parameters
        pftcon_read!(pftcon, ds)

        # Read initVertical scalar parameters
        readParams_initVertical!(ds)

        # Read saturated excess runoff parameters (SaturatedExcessRunoffMod)
        if haskey(ds, "fff")
            sat_excess_runoff_params.fff = ds["fff"][1]
        end

        # Read surface-water connectivity parameters (SurfaceWaterMod readParams).
        # clm5_params ships pc=0.4, mu=0.13889; the old hardcoded 0.5/1.0 throttled
        # h2osfc surface runoff on cold sites (melt water over-infiltrated).
        haskey(ds, "pc") && (SURFACE_WATER_PC[] = Float64(ds["pc"][1]))
        haskey(ds, "mu") && (SURFACE_WATER_MU[] = Float64(ds["mu"][1]))

        # Read photosynthesis parameters (PhotosynthesisMod)
        readParams_photosynthesis!(ds)

        # Read soil hydrology parameters (SoilHydrologyMod)
        haskey(ds, "n_baseflow") && (soilhydrology_params.n_baseflow = Float64(ds["n_baseflow"][1]))
        haskey(ds, "e_ice") && (soilhydrology_params.e_ice = Float64(ds["e_ice"][1]))
        haskey(ds, "perched_baseflow_scalar") && (soilhydrology_params.perched_baseflow_scalar = Float64(ds["perched_baseflow_scalar"][1]))
        haskey(ds, "aq_sp_yield_min") && (soilhydrology_params.aq_sp_yield_min = Float64(ds["aq_sp_yield_min"][1]))

        # Read snow hydrology parameters (SnowHydrologyMod)
        haskey(ds, "wimp") && (snowhydrology_params.wimp = Float64(ds["wimp"][1]))
        # Overburden-compaction viscosity params. ceta in particular DEFAULTS to 250
        # in the struct but the clm5 params file ships 450 — using the default makes
        # the Vionnet2012 viscosity eta = f1·f2·(bi/ceta)·… ~1.8× too large, i.e. the
        # snow pack under-compacts and snow_depth runs deep. Read them from the file.
        # upplim_destruct_metamorph: struct default is 100 but the clm5 file ships 175.
        # Destructive metamorphism is throttled by exp(-46e-3·(bi-upplim)) once bi>upplim,
        # so the too-low default cuts ddz1 ~10× for accumulation-season snow in the
        # 100–175 kg/m3 range → the pack under-densifies and SNOW_DEPTH runs deep.
        haskey(ds, "upplim_destruct_metamorph") && (snowhydrology_params.upplim_destruct_metamorph = Float64(ds["upplim_destruct_metamorph"][1]))
        haskey(ds, "ceta") && (snowhydrology_params.ceta = Float64(ds["ceta"][1]))
        haskey(ds, "eta0_vionnet") && (snowhydrology_params.eta0_vionnet = Float64(ds["eta0_vionnet"][1]))
        haskey(ds, "eta0_anderson") && (snowhydrology_params.eta0_anderson = Float64(ds["eta0_anderson"][1]))
        haskey(ds, "wind_snowcompact_fact") && (snowhydrology_params.wind_snowcompact_fact = Float64(ds["wind_snowcompact_fact"][1]))
        # van Kampenhout wind-drift compaction params (Fortran reads drift_gs/rho_max/
        # tau_ref via readNcdioScalar in SnowHydrology_readnl). drift_gs in particular
        # ships 0.00035 in the clm5 file but the struct default is 0.35 — a 1000× error.
        # It enters the drift mobility as −0.583·drift_gs, so the too-large default
        # suppresses wind-drift compaction; windy fresh snow (e.g. Baltimore) then stays
        # fluffy far too long → the pack sits below the 50 kg/m³ CombineSnowLayers
        # persistence threshold and snow_depth runs deep.
        haskey(ds, "drift_gs") && (snowhydrology_params.drift_gs = Float64(ds["drift_gs"][1]))
        haskey(ds, "rho_max") && (snowhydrology_params.rho_max = Float64(ds["rho_max"][1]))
        haskey(ds, "tau_ref") && (snowhydrology_params.tau_ref = Float64(ds["tau_ref"][1]))
        # SNOW_DENSITY_MAX / SNOW_DENSITY_MIN are CLM.jl calibration aliases for
        # rho_max / rho_min (calibration/param_injection.jl). They are read AFTER the
        # Fortran-named rho_max so the calibrated knob wins when a file carries both —
        # the precedence clm_run! had before this wiring was centralized here.
        haskey(ds, "SNOW_DENSITY_MAX") && (snowhydrology_params.rho_max = Float64(ds["SNOW_DENSITY_MAX"][1]))
        haskey(ds, "SNOW_DENSITY_MIN") && (snowhydrology_params.rho_min = Float64(ds["SNOW_DENSITY_MIN"][1]))
        haskey(ds, "ssi") && (snowhydrology_params.ssi = Float64(ds["ssi"][1]))
        # Aerosol scavenging factors (SnowHydrologyMod readParams) — consumed by
        # snow_hydrology's aerosol flushing.
        for f in (:scvng_fct_mlt_sf, :scvng_fct_mlt_bcphi, :scvng_fct_mlt_bcpho,
                  :scvng_fct_mlt_dst1, :scvng_fct_mlt_dst2, :scvng_fct_mlt_dst3,
                  :scvng_fct_mlt_dst4)
            k = String(f)
            haskey(ds, k) && setfield!(snowhydrology_params, f, Float64(ds[k][1]))
        end

        # SNICAR snow-aging / grain-radius parameters (SnowSnicarMod readParams).
        readParams_snicar!(ds)

        # NOTE: zsno/zlnd/zglc (FrictionVelocityMod) are per-instance roughness
        # lengths — clm_run! reads SNO_Z0MV into inst.frictionvel.zsno, which needs
        # `inst`. zlnd is also a WaterDiagnosticBulk param, wired here.
        haskey(ds, "zlnd") && (waterdiagbulk_params.zlnd = Float64(ds["zlnd"][1]))

        # Read interception parameters (CanopyHydrologyMod)
        haskey(ds, "interception_fraction") && (canopy_hydrology_params.interception_fraction = Float64(ds["interception_fraction"][1]))
        haskey(ds, "maximum_leaf_wetted_fraction") && (canopy_hydrology_params.maximum_leaf_wetted_fraction = Float64(ds["maximum_leaf_wetted_fraction"][1]))
        haskey(ds, "snowcan_unload_wind_fact") && (canopy_hydrology_params.snowcan_unload_wind_fact = Float64(ds["snowcan_unload_wind_fact"][1]))
        haskey(ds, "snowcan_unload_temp_fact") && (canopy_hydrology_params.snowcan_unload_temp_fact = Float64(ds["snowcan_unload_temp_fact"][1]))
        haskey(ds, "liq_canopy_storage_scalar") && (canopy_hydrology_params.liq_canopy_storage_scalar = Float64(ds["liq_canopy_storage_scalar"][1]))
        haskey(ds, "snow_canopy_storage_scalar") && (canopy_hydrology_params.snow_canopy_storage_scalar = Float64(ds["snow_canopy_storage_scalar"][1]))

        # Bare-ground flux params. Critically, wind_min (default 1.0 m/s) was not
        # being wired, so it stayed 0 — at low wind the reference wind was not
        # floored, collapsing the bare-ground friction velocity (and aerodynamic
        # conductance) and over-heating the summer surface.
        haskey(ds, "a_coef")   && (bareground_fluxes_params.a_coef   = Float64(ds["a_coef"][1]))
        haskey(ds, "a_exp")    && (bareground_fluxes_params.a_exp    = Float64(ds["a_exp"][1]))
        haskey(ds, "wind_min") && (bareground_fluxes_params.wind_min = Float64(ds["wind_min"][1]))

        # Canopy flux params — SEPARATE struct from the bare-ground one and likewise
        # not wired from file: a_coef/a_exp defaulted to 0.5/1.0 instead of the file's
        # 0.13/0.45, so the under-canopy soil-drag csoilb used a LINEAR (^1.0) Reynolds
        # term instead of ^0.45 → rah_below ~13x too high for sparse canopies (grass)
        # → the warm ground could not heat the canopy air → grass T_VEG too cold at
        # low sun. (csoilc/cv defaults already match the file but wire them too.)
        haskey(ds, "a_coef")   && (canopy_fluxes_params.a_coef = Float64(ds["a_coef"][1]))
        haskey(ds, "a_exp")    && (canopy_fluxes_params.a_exp  = Float64(ds["a_exp"][1]))
        haskey(ds, "csoilc")   && (canopy_fluxes_params.csoilc = Float64(ds["csoilc"][1]))
        haskey(ds, "cv")       && (canopy_fluxes_params.cv     = Float64(ds["cv"][1]))
        # Dry-litter layer (CanopyFluxesMod readParams) — the litter resistance term.
        haskey(ds, "lai_dl")   && (canopy_fluxes_params.lai_dl = Float64(ds["lai_dl"][1]))
        haskey(ds, "z_dl")     && (canopy_fluxes_params.z_dl   = Float64(ds["z_dl"][1]))

        # Dry-surface-layer soil-evaporation resistance (SurfaceResistanceMod readParams)
        haskey(ds, "d_max") && (surface_resistance_params.d_max = Float64(ds["d_max"][1]))
        haskey(ds, "frac_sat_soil_dsl_init") &&
            (surface_resistance_params.frac_sat_soil_dsl_init =
                 Float64(ds["frac_sat_soil_dsl_init"][1]))

        # Time-constant soil properties (SoilStateInitTimeConstMod readParams)
        readParams_soilstate!(ds)

        # CN phenology parameters (CNPhenologyMod readParams)
        readParams_CNPhenology!(cn_phenology_params, ds)

        # NOTE: fire params (prh30, ignition_efficiency; CNFireBaseMod readParams) are
        # NOT wired here because no module-level CNFireParams instance exists — the
        # fire methods take a caller-supplied `cnfire_params` and clm_driver! does not
        # instantiate one (fire is driven only from explicit harness/test calls).
        # Callers that do run fire should populate theirs with readParams_CNFire!.
        #
        # NOTE: sand_pf / clay_pf (SoilStateInitTimeConstMod) are PPE perturbation
        # knobs for cellsand/cellclay and are 0 (no-op) in the shipped params file;
        # the perturbation branch is deliberately not ported.

        # LUNA photosynthetic-N acclimation parameters (jmaxb0/jmaxb1/wc2wjb0 + scalars)
        luna_read_params!(luna_params_inst, ds)
    finally
        close(ds)
    end
    nothing
end

# --- Helper: read a 1D PFT variable from NetCDF ---

"""
    _read_pft_var!(ds, varname, dest; offset=0)

Read a 1D variable dimensioned on `pft` from the dataset and store it into
`dest`. NetCDF arrays are 0-based PFT (size npft_file), Julia arrays are
1-based (size MXPFT+1). `offset` shifts the start position in `dest`.
"""
function _read_pft_var!(ds::NCDataset, varname::String, dest::Vector{Float64};
                         offset::Int=0)
    if haskey(ds, varname)
        data = Array(ds[varname])
        npft_file = length(data)
        n = min(npft_file, length(dest) - offset)
        for i in 1:n
            val = data[i]
            dest[i + offset] = ismissing(val) ? 0.0 : Float64(val)
        end
    end
    nothing
end

function _read_pft_var_int!(ds::NCDataset, varname::String, dest::Vector{Int};
                             offset::Int=0)
    if haskey(ds, varname)
        data = Array(ds[varname])
        npft_file = length(data)
        n = min(npft_file, length(dest) - offset)
        for i in 1:n
            val = data[i]
            dest[i + offset] = ismissing(val) ? 0 : Int(round(val))
        end
    end
    nothing
end

"""
    _read_scalar(ds, varname, default) -> Float64

Read a scalar variable from the dataset with a fallback default.
"""
function _read_scalar(ds::NCDataset, varname::String, default::Float64)
    if haskey(ds, varname)
        v = Array(ds[varname])
        val = v isa AbstractArray ? v[1] : v
        return ismissing(val) ? default : Float64(val)
    end
    return default
end

"""
    pftcon_read!(p::PftconType, ds::NCDataset)

Populate PftconType arrays from an open clm5_params.nc dataset.
Corresponds to Fortran pftconMod::InitRead.
"""
function pftcon_read!(p::PftconType, ds::NCDataset)
    # Ensure pftcon is allocated
    if isempty(p.slatop)
        pftcon_allocate!(p)
    end

    # --- 1D float PFT variables ---
    # Mapping: NetCDF varname => pftcon field
    float_vars = [
        ("z0mr",           p.z0mr),
        ("displar",        p.displar),
        ("dleaf",          p.dleaf),
        ("c3psn",          p.c3psn),
        ("xl",             p.xl),
        ("roota_par",      p.roota_par),
        ("rootb_par",      p.rootb_par),
        ("slatop",         p.slatop),
        ("dsladlai",       p.dsladlai),
        ("leafcn",         p.leafcn),
        ("biofuel_harvfrac", p.biofuel_harvfrac),
        ("flnr",           p.flnr),
        ("smpso",          p.smpso),
        ("smpsc",          p.smpsc),
        ("fnitr",          p.fnitr),
        ("woody",          p.woody),
        ("lflitcn",        p.lflitcn),
        ("frootcn",        p.frootcn),
        ("livewdcn",       p.livewdcn),
        ("deadwdcn",       p.deadwdcn),
        ("grperc",         p.grperc),
        ("grpnow",         p.grpnow),
        ("froot_leaf",     p.froot_leaf),
        ("stem_leaf",      p.stem_leaf),
        ("croot_stem",     p.croot_stem),
        ("flivewd",        p.flivewd),
        ("fcur",           p.fcur),
        ("fcurdv",         p.fcurdv),
        ("lf_flab",        p.lf_flab),
        ("lf_fcel",        p.lf_fcel),
        ("lf_flig",        p.lf_flig),
        ("fr_flab",        p.fr_flab),
        ("fr_fcel",        p.fr_fcel),
        ("fr_flig",        p.fr_flig),
        ("leaf_long",      p.leaf_long),
        ("evergreen",      p.evergreen),
        ("stress_decid",   p.stress_decid),
        ("season_decid",   p.season_decid),
        ("season_decid_temperate", p.season_decid_temperate),
        ("crit_onset_gdd_sf", p.crit_onset_gdd_sf),
        ("ndays_on",       p.ndays_on),
        ("pftpar20",       p.pftpar20),
        ("pftpar28",       p.pftpar28),
        ("pftpar29",       p.pftpar29),
        ("pftpar30",       p.pftpar30),
        ("pftpar31",       p.pftpar31),
        ("a_fix",          p.a_fix),
        ("b_fix",          p.b_fix),
        ("c_fix",          p.c_fix),
        ("s_fix",          p.s_fix),
        ("akc_active",     p.akc_active),
        ("akn_active",     p.akn_active),
        ("ekc_active",     p.ekc_active),
        ("ekn_active",     p.ekn_active),
        ("kc_nonmyc",      p.kc_nonmyc),
        ("kn_nonmyc",      p.kn_nonmyc),
        ("kr_resorb",      p.kr_resorb),
        ("perecm",         p.perecm),
        ("fun_cn_flex_a",  p.fun_cn_flex_a),
        ("fun_cn_flex_b",  p.fun_cn_flex_b),
        ("fun_cn_flex_c",  p.fun_cn_flex_c),
        ("FUN_fracfixers", p.FUN_fracfixers),
        ("manunitro",      p.manunitro),
        ("fleafcn",        p.fleafcn),
        ("ffrootcn",       p.ffrootcn),
        ("fstemcn",        p.fstemcn),
        ("pconv",          p.pconv),
        ("pprod10",        p.pprod10),
        ("pprod100",       p.pprod100),
        ("pprodharv10",    p.pprodharv10),
        ("graincn",        p.graincn),
        ("mxtmp",          p.mxtmp),
        ("baset",          p.baset),
        ("declfact",       p.declfact),
        ("bfact",          p.bfact),
        ("aleaff",         p.aleaff),
        ("arootf",         p.arootf),
        ("astemf",         p.astemf),
        ("arooti",         p.arooti),
        ("fleafi",         p.fleafi),
        ("allconsl",       p.allconsl),
        ("allconss",       p.allconss),
        ("crop",           p.crop),
        ("irrigated",      p.irrigated),
        ("ztopmx",         p.ztopmx),
        ("laimx",          p.laimx),
        ("gddmin",         p.gddmin),
        ("hybgdd",         p.hybgdd),
        ("lfemerg",        p.lfemerg),
        ("grnfill",        p.grnfill),
        ("mbbopt",         p.mbbopt),
        ("medlynslope",    p.medlynslope),
        ("medlynintercept", p.medlynintercept),
        ("cc_leaf",        p.cc_leaf),
        ("cc_lstem",       p.cc_lstem),
        ("cc_dstem",       p.cc_dstem),
        ("cc_other",       p.cc_other),
        ("fm_leaf",        p.fm_leaf),
        ("fm_lstem",       p.fm_lstem),
        ("fm_dstem",       p.fm_dstem),
        ("fm_other",       p.fm_other),
        ("fm_root",        p.fm_root),
        ("fm_lroot",       p.fm_lroot),
        ("fm_droot",       p.fm_droot),
        ("fsr_pft",        p.fsr_pft),
        ("fd_pft",         p.fd_pft),
        ("rswf_min",       p.rswf_min),
        ("rswf_max",       p.rswf_max),
        ("nstem",          p.nstem),
        ("taper",          p.taper),
        # Raupach92 roughness parameters
        ("z0v_Cr",         p.z0v_Cr),
        ("z0v_Cs",         p.z0v_Cs),
        ("z0v_c",          p.z0v_c),
        ("z0v_cw",         p.z0v_cw),
        ("z0v_LAImax",     p.z0v_LAImax),
        ("z0v_LAIoff",     p.z0v_LAIoff),
        # Planting temperature
        ("planting_temp",       p.planttemp),
        ("min_planting_temp",   p.minplanttemp),
        # Flexible CN
        ("i_vcad",         p.i_vcad),
        ("s_vcad",         p.s_vcad),
        ("i_flnr",         p.i_flnr),
        ("s_flnr",         p.s_flnr),
    ]

    for (varname, dest) in float_vars
        _read_pft_var!(ds, varname, dest)
    end

    # --- 1D integer PFT variables ---
    _read_pft_var_int!(ds, "mergetoclmpft", p.mergetoclmpft)
    _read_pft_var_int!(ds, "mxmat", p.mxmat)
    _read_pft_var_int!(ds, "min_NH_planting_date", p.mnNHplantdate)
    _read_pft_var_int!(ds, "max_NH_planting_date", p.mxNHplantdate)
    _read_pft_var_int!(ds, "min_SH_planting_date", p.mnSHplantdate)
    _read_pft_var_int!(ds, "max_SH_planting_date", p.mxSHplantdate)

    # --- 2D optical properties (VIS=1, NIR=2) ---
    IVIS = 1; INIR = 2
    _read_pft_col!(ds, "rholvis", p.rhol, IVIS)
    _read_pft_col!(ds, "rholnir", p.rhol, INIR)
    _read_pft_col!(ds, "rhosvis", p.rhos, IVIS)
    _read_pft_col!(ds, "rhosnir", p.rhos, INIR)
    _read_pft_col!(ds, "taulvis", p.taul, IVIS)
    _read_pft_col!(ds, "taulnir", p.taul, INIR)
    _read_pft_col!(ds, "tausvis", p.taus, IVIS)
    _read_pft_col!(ds, "tausnir", p.taus, INIR)

    # --- 2D rootprof_beta [pft, variant] ---
    if haskey(ds, "rootprof_beta")
        data = Array(ds["rootprof_beta"])
        npft_file = size(data, 1)
        nvar = size(data, 2)
        for v in 1:min(nvar, size(p.rootprof_beta, 2))
            for i in 1:min(npft_file, size(p.rootprof_beta, 1))
                p.rootprof_beta[i, v] = Float64(data[i, v])
            end
        end
    end

    # --- Constants that don't come from file ---
    n = length(p.dwood)
    for i in 1:n
        p.dwood[i] = DWOOD_PARAM
        p.root_radius[i] = ROOT_RADIUS_PARAM
        p.root_density[i] = ROOT_DENSITY_PARAM
    end

    # --- Set vegetation type flags ---
    for i in 1:n
        p.is_tree[i]  = (p.woody[i] > 0.5) && (p.slatop[i] > 0.0)
        p.is_shrub[i] = false  # derived from type index later if needed
        p.is_grass[i] = (p.woody[i] < 0.5) && (p.crop[i] < 0.5) && (i > 1)
    end

    # Set mergetoclmpft-based flags
    set_is_pft_known_to_model!(p)

    nothing
end

"""
    _read_pft_col!(ds, varname, dest_matrix, col_idx)

Read a 1D PFT variable and store it into column `col_idx` of a 2D matrix.
"""
function _read_pft_col!(ds::NCDataset, varname::String, dest::Matrix{Float64}, col_idx::Int)
    if haskey(ds, varname)
        data = Array(ds[varname])
        npft_file = length(data)
        for i in 1:min(npft_file, size(dest, 1))
            val = data[i]
            dest[i, col_idx] = ismissing(val) ? 0.0 : Float64(val)
        end
    end
    nothing
end

# --- initVertical scalar parameters ---
# Module-level storage (replaces Fortran params_inst)
const _initvert_slopebeta  = Ref(0.0)
const _initvert_slopemax   = Ref(0.0)
const _initvert_zbedrock   = Ref(-1.0)   # negative means don't substitute
const _initvert_zbedrock_sf = Ref(1.0)

"""
    readParams_initVertical!(ds::NCDataset)

Read initVertical scalar parameters from the parameter file.
"""
function readParams_initVertical!(ds::NCDataset)
    _initvert_slopebeta[]  = _read_scalar(ds, "slopebeta", 0.0)
    _initvert_slopemax[]   = _read_scalar(ds, "slopemax", 0.0)
    _initvert_zbedrock[]   = _read_scalar(ds, "zbedrock", -1.0)
    _initvert_zbedrock_sf[] = _read_scalar(ds, "zbedrock_sf", 1.0)
    nothing
end

"""
    readParams_snicar!(ds::NCDataset)

Read the SNICAR snow-aging / grain-radius parameters (SnowSnicarMod readParams)
into the module-level `snicar_params`, and the shared `snw_rds_min` into
`waterdiagbulk_params` (WaterDiagnosticBulkType readParams reads the same
variable).

`xdrdt` is the arbitrary factor applied to the snow-aging rate `dr`
(SnowSnicarMod.F90:1659). `snw_aging_bst` is accepted as an alias because the
CLM.jl calibration layer (`calibration/param_injection.jl`) writes the knob
under that name.
"""
function readParams_snicar!(ds::NCDataset)
    # snw_aging_bst takes precedence: a calibrated file carries BOTH the stock
    # `xdrdt` and the injected `snw_aging_bst`, and the calibrated knob must win
    # (this is the precedence clm_run! used before the wiring moved here).
    if haskey(ds, "snw_aging_bst")
        snicar_params.xdrdt = Float64(ds["snw_aging_bst"][1])
    elseif haskey(ds, "xdrdt")
        snicar_params.xdrdt = Float64(ds["xdrdt"][1])
    end
    haskey(ds, "snw_rds_refrz") &&
        (snicar_params.snw_rds_refrz = Float64(ds["snw_rds_refrz"][1]))
    haskey(ds, "C2_liq_Brun89") &&
        (snicar_params.C2_liq_Brun89 = Float64(ds["C2_liq_Brun89"][1]))
    haskey(ds, "fresh_snw_rds_max") &&
        (snicar_params.fresh_snw_rds_max = Float64(ds["fresh_snw_rds_max"][1]))
    if haskey(ds, "snw_rds_min")
        v = Float64(ds["snw_rds_min"][1])
        snicar_params.snw_rds_min = v
        waterdiagbulk_params.snw_rds_min = v
    end
    nothing
end

"""
    readParams_soilstate!(ds::NCDataset)

Read the time-constant soil-property parameters (SoilStateInitTimeConstMod
readParams) into the module-level `soilstate_params`, consumed by
`init_soil_properties!`.
"""
function readParams_soilstate!(ds::NCDataset)
    for f in (:tkd_sand, :tkd_clay, :tkd_om, :tkm_om,
              :csol_sand, :csol_clay, :csol_om, :pd,
              :bsw_sf, :hksat_sf, :sucsat_sf, :watsat_sf, :om_frac_sf)
        k = String(f)
        haskey(ds, k) && setfield!(soilstate_params, f, Float64(ds[k][1]))
    end
    nothing
end

"""
    readParams_CNPhenology!(p::PhenologyParams, ds::NCDataset)

Read the CN phenology parameters (CNPhenologyMod readParams) into `p`.

Note `lwtop_ann` in the file maps to the `lwtop` field, matching Fortran
(`readNcdioScalar(ncid, 'lwtop_ann', subname, params_inst%lwtop)`).
"""
function readParams_CNPhenology!(p, ds::NCDataset)  # p::PhenologyParams (defined later in CLM.jl)
    for f in (:crit_dayl, :crit_dayl_at_high_lat, :crit_dayl_lat_slope,
              :ndays_off, :fstor2tran, :crit_onset_fdd, :crit_onset_swi,
              :soilpsi_on, :crit_offset_fdd, :crit_offset_swi, :soilpsi_off,
              :phenology_soil_depth, :snow5d_thresh_for_onset)
        k = String(f)
        haskey(ds, k) && setfield!(p, f, Float64(ds[k][1]))
    end
    haskey(ds, "lwtop_ann") && (p.lwtop = Float64(ds["lwtop_ann"][1]))
    nothing
end

"""
    readParams_CNFire!(p::CNFireParams, ds::NCDataset)

Read the fire parameters (CNFireBaseMod readParams) into `p`. Not called from
`readParameters!` — there is no module-level CNFireParams instance; harnesses
that run fire should call this on the instance they pass to `cn_driver_*`.
"""
function readParams_CNFire!(p, ds::NCDataset)  # p::CNFireParams (defined later in CLM.jl)
    haskey(ds, "prh30") && (p.prh30 = Float64(ds["prh30"][1]))
    haskey(ds, "ignition_efficiency") &&
        (p.ignition_efficiency = Float64(ds["ignition_efficiency"][1]))
    nothing
end

"""
    readParams_photosynthesis!(ds::NCDataset)

Read photosynthesis module parameters from the parameter file.
Initializes the global `params_inst` (PhotoParamsData).
"""
function readParams_photosynthesis!(ds::NCDataset)
    photo_params_init!(params_inst)

    # 1D PFT parameters
    _read_pft_var!(ds, "theta_cj", params_inst.theta_cj)
    _read_pft_var!(ds, "krmax", params_inst.krmax)
    _read_pft_var!(ds, "lmr_intercept_atkin", params_inst.lmr_intercept_atkin)

    # Scalar parameters (kinetics, activation/deactivation energies)
    params_inst.theta_ip     = _read_scalar(ds, "theta_ip", 0.95)
    params_inst.act25        = _read_scalar(ds, "act25", 72.0)
    params_inst.fnr          = _read_scalar(ds, "fnr", 7.16)
    params_inst.cp25_yr2000  = _read_scalar(ds, "cp25_yr2000", 42.75e-6)
    params_inst.kc25_coef    = _read_scalar(ds, "kc25_coef", 404.9e-6)
    params_inst.ko25_coef    = _read_scalar(ds, "ko25_coef", 278.4e-3)
    params_inst.fnps         = _read_scalar(ds, "fnps", 0.15)
    params_inst.theta_psii   = _read_scalar(ds, "theta_psii", 0.7)
    params_inst.vcmaxha      = _read_scalar(ds, "vcmaxha", 65330.0)
    params_inst.jmaxha       = _read_scalar(ds, "jmaxha", 43540.0)
    params_inst.tpuha        = _read_scalar(ds, "tpuha", 53100.0)
    params_inst.lmrha        = _read_scalar(ds, "lmrha", 46390.0)
    params_inst.kcha         = _read_scalar(ds, "kcha", 79430.0)
    params_inst.koha         = _read_scalar(ds, "koha", 36380.0)
    params_inst.cpha         = _read_scalar(ds, "cpha", 37830.0)
    params_inst.vcmaxhd      = _read_scalar(ds, "vcmaxhd", 149250.0)
    params_inst.jmaxhd       = _read_scalar(ds, "jmaxhd", 152040.0)
    params_inst.tpuhd        = _read_scalar(ds, "tpuhd", 150650.0)
    params_inst.lmrhd        = _read_scalar(ds, "lmrhd", 150650.0)
    params_inst.lmrse        = _read_scalar(ds, "lmrse", 490.0)
    params_inst.tpu25ratio   = _read_scalar(ds, "tpu25ratio", 0.167)
    params_inst.kp25ratio    = _read_scalar(ds, "kp25ratio", 20160.0)
    params_inst.vcmaxse_sf   = _read_scalar(ds, "vcmaxse_sf", 1.0)
    params_inst.jmaxse_sf    = _read_scalar(ds, "jmaxse_sf", 1.0)
    params_inst.tpuse_sf     = _read_scalar(ds, "tpuse_sf", 1.0)
    params_inst.jmax25top_sf = _read_scalar(ds, "jmax25top_sf", 1.0)

    # 2D parameters (pft × nvegwcs): kmax, psi50, ck
    if haskey(ds, "kmax")
        data = Array(ds["kmax"])
        data = replace(data, missing => NaN)
        n1 = min(size(data, 1), size(params_inst.kmax, 1))
        n2 = min(size(data, 2), size(params_inst.kmax, 2))
        params_inst.kmax[1:n1, 1:n2] .= Float64.(data[1:n1, 1:n2])
    end
    if haskey(ds, "psi50")
        data = Array(ds["psi50"])
        data = replace(data, missing => NaN)
        n1 = min(size(data, 1), size(params_inst.psi50, 1))
        n2 = min(size(data, 2), size(params_inst.psi50, 2))
        params_inst.psi50[1:n1, 1:n2] .= Float64.(data[1:n1, 1:n2])
    end
    if haskey(ds, "ck")
        data = Array(ds["ck"])
        data = replace(data, missing => NaN)
        n1 = min(size(data, 1), size(params_inst.ck, 1))
        n2 = min(size(data, 2), size(params_inst.ck, 2))
        params_inst.ck[1:n1, 1:n2] .= Float64.(data[1:n1, 1:n2])
    end

    nothing
end
