using Test
using NCDatasets
using CLM

# ==========================================================================
# Regression guard for the "param parsed but thrown away" bug class.
#
# Fortran reads ~110 scalars from the params file via readNcdioScalar. Several
# of them were, at various times, either never read by read_params.jl or read
# and discarded (`haskey(ds, "x") && nothing`), so the model silently ran on a
# hardcoded default — this caused a real 1000x error once (drift_gs) and a 2x
# error in the CN leaf-offset ramp (ndays_off: file 15 d, Julia default 30 d).
#
# This test writes a synthetic params file whose every scalar differs from the
# Julia struct default, runs readParameters!, and asserts the value actually
# landed in the struct the physics reads. Add any newly-wired scalar here.
# ==========================================================================

# varname => (getter, sentinel value distinct from the struct default)
const _PARAM_CASES = [
    # SnowHydrologyMod
    ("wimp",                      () -> CLM.snowhydrology_params.wimp,                       0.0123),
    ("ssi",                       () -> CLM.snowhydrology_params.ssi,                        0.0411),
    ("drift_gs",                  () -> CLM.snowhydrology_params.drift_gs,                   0.00042),
    ("rho_max",                   () -> CLM.snowhydrology_params.rho_max,                    357.0),
    ("tau_ref",                   () -> CLM.snowhydrology_params.tau_ref,                    172801.0),
    ("ceta",                      () -> CLM.snowhydrology_params.ceta,                       451.0),
    ("eta0_vionnet",              () -> CLM.snowhydrology_params.eta0_vionnet,               7.6e6),
    ("eta0_anderson",             () -> CLM.snowhydrology_params.eta0_anderson,              9.1e5),
    ("wind_snowcompact_fact",     () -> CLM.snowhydrology_params.wind_snowcompact_fact,      5.5),
    ("upplim_destruct_metamorph", () -> CLM.snowhydrology_params.upplim_destruct_metamorph,  176.0),
    ("scvng_fct_mlt_sf",          () -> CLM.snowhydrology_params.scvng_fct_mlt_sf,           0.51),
    ("scvng_fct_mlt_bcphi",       () -> CLM.snowhydrology_params.scvng_fct_mlt_bcphi,        0.21),
    ("scvng_fct_mlt_bcpho",       () -> CLM.snowhydrology_params.scvng_fct_mlt_bcpho,        0.031),
    ("scvng_fct_mlt_dst1",        () -> CLM.snowhydrology_params.scvng_fct_mlt_dst1,         0.021),
    ("scvng_fct_mlt_dst2",        () -> CLM.snowhydrology_params.scvng_fct_mlt_dst2,         0.022),
    ("scvng_fct_mlt_dst3",        () -> CLM.snowhydrology_params.scvng_fct_mlt_dst3,         0.011),
    ("scvng_fct_mlt_dst4",        () -> CLM.snowhydrology_params.scvng_fct_mlt_dst4,         0.012),
    # SnowSnicarMod — the snow-aging params (xdrdt is `snw_aging_bst` in the
    # CLM.jl calibration layer). xdrdt multiplies the grain-growth rate dr.
    ("xdrdt",                     () -> CLM.snicar_params.xdrdt,                             5.0),
    ("snw_rds_refrz",             () -> CLM.snicar_params.snw_rds_refrz,                     1500.0),
    ("C2_liq_Brun89",             () -> CLM.snicar_params.C2_liq_Brun89,                     4.3e-13),
    ("fresh_snw_rds_max",         () -> CLM.snicar_params.fresh_snw_rds_max,                 70.232),
    ("snw_rds_min",               () -> CLM.snicar_params.snw_rds_min,                       55.0),
    # SoilHydrologyMod
    ("n_baseflow",                () -> CLM.soilhydrology_params.n_baseflow,                 2.7),
    ("e_ice",                     () -> CLM.soilhydrology_params.e_ice,                      3.5),
    ("perched_baseflow_scalar",   () -> CLM.soilhydrology_params.perched_baseflow_scalar,    1.7e-6),
    ("aq_sp_yield_min",           () -> CLM.soilhydrology_params.aq_sp_yield_min,            0.021),
    # SaturatedExcessRunoffMod / SurfaceWaterMod
    ("fff",                       () -> CLM.sat_excess_runoff_params.fff,                    0.11),
    ("pc",                        () -> CLM.SURFACE_WATER_PC[],                              0.34),
    ("mu",                        () -> CLM.SURFACE_WATER_MU[],                              0.14),
    # CanopyHydrologyMod
    ("interception_fraction",         () -> CLM.canopy_hydrology_params.interception_fraction,         0.65),
    ("maximum_leaf_wetted_fraction",  () -> CLM.canopy_hydrology_params.maximum_leaf_wetted_fraction,  0.073),
    ("snowcan_unload_wind_fact",      () -> CLM.canopy_hydrology_params.snowcan_unload_wind_fact,      0.51),
    ("snowcan_unload_temp_fact",      () -> CLM.canopy_hydrology_params.snowcan_unload_temp_fact,      187001.0),
    ("liq_canopy_storage_scalar",     () -> CLM.canopy_hydrology_params.liq_canopy_storage_scalar,     0.11),
    ("snow_canopy_storage_scalar",    () -> CLM.canopy_hydrology_params.snow_canopy_storage_scalar,    6.1),
    # BareGroundFluxesMod / CanopyFluxesMod (a_coef/a_exp feed both structs)
    ("wind_min",                  () -> CLM.bareground_fluxes_params.wind_min,               1.1),
    ("csoilc",                    () -> CLM.canopy_fluxes_params.csoilc,                     0.0041),
    ("cv",                        () -> CLM.canopy_fluxes_params.cv,                         0.011),
    ("lai_dl",                    () -> CLM.canopy_fluxes_params.lai_dl,                     0.55),
    ("z_dl",                      () -> CLM.canopy_fluxes_params.z_dl,                       0.055),
    # SurfaceResistanceMod
    ("d_max",                     () -> CLM.surface_resistance_params.d_max,                 16.0),
    ("frac_sat_soil_dsl_init",    () -> CLM.surface_resistance_params.frac_sat_soil_dsl_init, 0.81),
    # SoilStateInitTimeConstMod
    ("tkd_sand",                  () -> CLM.soilstate_params.tkd_sand,                       8.9),
    ("tkd_clay",                  () -> CLM.soilstate_params.tkd_clay,                       2.93),
    ("tkd_om",                    () -> CLM.soilstate_params.tkd_om,                         0.051),
    ("tkm_om",                    () -> CLM.soilstate_params.tkm_om,                         0.26),
    ("csol_sand",                 () -> CLM.soilstate_params.csol_sand,                      2.129),
    ("csol_clay",                 () -> CLM.soilstate_params.csol_clay,                      2.386),
    ("csol_om",                   () -> CLM.soilstate_params.csol_om,                        2.6),
    ("pd",                        () -> CLM.soilstate_params.pd,                             2701.0),
    ("bsw_sf",                    () -> CLM.soilstate_params.bsw_sf,                         1.1),
    ("hksat_sf",                  () -> CLM.soilstate_params.hksat_sf,                       1.2),
    ("sucsat_sf",                 () -> CLM.soilstate_params.sucsat_sf,                      1.3),
    ("watsat_sf",                 () -> CLM.soilstate_params.watsat_sf,                      1.4),
    ("om_frac_sf",                () -> CLM.soilstate_params.om_frac_sf,                     1.5),
    # CNPhenologyMod (module-level instance consumed by cn_driver_no_leaching!)
    ("crit_dayl",                 () -> CLM.cn_phenology_params.crit_dayl,                   39300.0),
    ("crit_dayl_at_high_lat",     () -> CLM.cn_phenology_params.crit_dayl_at_high_lat,       54001.0),
    ("crit_dayl_lat_slope",       () -> CLM.cn_phenology_params.crit_dayl_lat_slope,         721.0),
    ("ndays_off",                 () -> CLM.cn_phenology_params.ndays_off,                   15.0),
    ("fstor2tran",                () -> CLM.cn_phenology_params.fstor2tran,                  0.51),
    ("crit_onset_fdd",            () -> CLM.cn_phenology_params.crit_onset_fdd,              16.0),
    ("crit_onset_swi",            () -> CLM.cn_phenology_params.crit_onset_swi,              17.0),
    ("soilpsi_on",                () -> CLM.cn_phenology_params.soilpsi_on,                  -0.61),
    ("crit_offset_fdd",           () -> CLM.cn_phenology_params.crit_offset_fdd,             18.0),
    ("crit_offset_swi",           () -> CLM.cn_phenology_params.crit_offset_swi,             19.0),
    ("soilpsi_off",               () -> CLM.cn_phenology_params.soilpsi_off,                 -0.81),
    ("lwtop_ann",                 () -> CLM.cn_phenology_params.lwtop,                       0.71),
    ("phenology_soil_depth",      () -> CLM.cn_phenology_params.phenology_soil_depth,        0.09),
    ("snow5d_thresh_for_onset",   () -> CLM.cn_phenology_params.snow5d_thresh_for_onset,     0.1),
    # WaterDiagnosticBulkType
    ("zlnd",                      () -> CLM.waterdiagbulk_params.zlnd,                       0.011),
    # PhotosynthesisMod (spot-check; the full set is read by readParams_photosynthesis!)
    ("act25",                     () -> CLM.params_inst.act25,                               60.0),
    ("fnr",                       () -> CLM.params_inst.fnr,                                 7.17),
    ("theta_ip",                  () -> CLM.params_inst.theta_ip,                            0.96),
    # initVerticalMod
    ("slopebeta",                 () -> CLM._initvert_slopebeta[],                           -3.1),
    ("slopemax",                  () -> CLM._initvert_slopemax[],                            0.41),
    ("zbedrock_sf",               () -> CLM._initvert_zbedrock_sf[],                         1.01),
]

"""Snapshot every global params struct readParameters! can mutate, so this test
does not leak calibrated values into the rest of the suite.

`deepcopy` is essential, not decorative: the PFT/photosynthesis params hold ARRAYS
that readParameters! fills IN PLACE, so a by-reference snapshot would restore
nothing. `pftcon` in particular is allocated by pftcon_read!; leaving it allocated
made a later driver test enter the irrigation branch (non-empty pftcon.irrigated)
that it is meant to skip."""
function _snapshot_globals()
    snap = Dict{Any,Any}()
    for obj in (CLM.snowhydrology_params, CLM.snicar_params, CLM.soilhydrology_params,
                CLM.sat_excess_runoff_params, CLM.canopy_hydrology_params,
                CLM.bareground_fluxes_params, CLM.canopy_fluxes_params,
                CLM.surface_resistance_params, CLM.soilstate_params,
                CLM.cn_phenology_params, CLM.waterdiagbulk_params,
                CLM.params_inst, CLM.luna_params_inst, CLM.pftcon)
        snap[obj] = Dict(f => deepcopy(getfield(obj, f)) for f in fieldnames(typeof(obj)))
    end
    snap[:refs] = (CLM.SURFACE_WATER_PC[], CLM.SURFACE_WATER_MU[],
                   CLM._initvert_slopebeta[], CLM._initvert_slopemax[],
                   CLM._initvert_zbedrock[], CLM._initvert_zbedrock_sf[])
    return snap
end

function _restore_globals!(snap)
    for (obj, fields) in snap
        obj === :refs && continue
        for (f, v) in fields
            setfield!(obj, f, v)
        end
    end
    (pc, mu, sb, sm, zb, zbsf) = snap[:refs]
    CLM.SURFACE_WATER_PC[] = pc;      CLM.SURFACE_WATER_MU[] = mu
    CLM._initvert_slopebeta[] = sb;   CLM._initvert_slopemax[] = sm
    CLM._initvert_zbedrock[] = zb;    CLM._initvert_zbedrock_sf[] = zbsf
    return nothing
end

"""Write a params file containing only the scalars under test."""
function _write_synthetic_params(path)
    NCDataset(path, "c") do ds
        for (name, _, val) in _PARAM_CASES
            defVar(ds, name, Float64(val), ())
        end
    end
    return path
end

@testset "read_params: scalars land in the structs the physics reads" begin
    snap = _snapshot_globals()
    try
        mktempdir() do dir
            pf = _write_synthetic_params(joinpath(dir, "synthetic_params.nc"))
            CLM.readParameters!(pf)
            for (name, getter, expected) in _PARAM_CASES
                @test getter() ≈ expected
            end
        end
    finally
        _restore_globals!(snap)
    end
end

@testset "read_params: snw_aging_bst alias feeds SNICAR xdrdt" begin
    # The CLM.jl calibration layer writes the snow-aging factor under the name
    # `snw_aging_bst` (calibration/param_injection.jl); Fortran's name is `xdrdt`
    # (SnowSnicarMod.F90:177). Both must reach snicar_params.xdrdt — this line used
    # to be `haskey(ds, "snw_aging_bst") && nothing`, discarding the value.
    snap = _snapshot_globals()
    try
        mktempdir() do dir
            pf = joinpath(dir, "alias_params.nc")
            NCDataset(pf, "c") do ds
                defVar(ds, "snw_aging_bst", 3.5, ())
            end
            CLM.readParameters!(pf)
            @test CLM.snicar_params.xdrdt ≈ 3.5
        end

        # A calibrated file carries BOTH the stock xdrdt and the injected
        # snw_aging_bst — the calibrated knob must win, or calibration is silently
        # ignored (the same failure mode as never reading it at all).
        mktempdir() do dir
            pf = joinpath(dir, "both_params.nc")
            NCDataset(pf, "c") do ds
                defVar(ds, "xdrdt", 1.0, ())
                defVar(ds, "snw_aging_bst", 4.5, ())
            end
            CLM.readParameters!(pf)
            @test CLM.snicar_params.xdrdt ≈ 4.5
        end

        # Same precedence for the SNOW_DENSITY_MAX alias over rho_max.
        mktempdir() do dir
            pf = joinpath(dir, "rho_params.nc")
            NCDataset(pf, "c") do ds
                defVar(ds, "rho_max", 350.0, ())
                defVar(ds, "SNOW_DENSITY_MAX", 275.0, ())
                defVar(ds, "SNOW_DENSITY_MIN", 60.0, ())
            end
            CLM.readParameters!(pf)
            @test CLM.snowhydrology_params.rho_max ≈ 275.0
            @test CLM.snowhydrology_params.rho_min ≈ 60.0
        end
    finally
        _restore_globals!(snap)
    end
end

@testset "read_params: every Fortran readNcdioScalar param is accounted for" begin
    # The complete set of scalars the Fortran reads from the params file
    # (`grep -rh readNcdioScalar src/**/*.F90` on CTSM). Every one must be either
    # (a) covered by _PARAM_CASES above — proven to reach the struct the physics
    # reads — or (b) in `known_unwired` with an explicit reason. This is what
    # turns the one-off sweep into a standing check: dropping a param out of the
    # wiring fails here.
    fortran_scalars = Set(split(
        "a_coef a_exp accum_factor act25 aq_sp_yield_min bsw_sf C2_liq_Brun89 ceta " *
        "clay_pf cp25_yr2000 cpha crit_dayl crit_dayl_at_high_lat crit_dayl_lat_slope " *
        "crit_offset_fdd crit_offset_swi crit_onset_fdd crit_onset_swi csoilc csol_clay " *
        "csol_om csol_sand cv d_max drift_gs e_ice enzyme_turnover_daily eta0_anderson " *
        "eta0_vionnet fff fnps fnr frac_sat_soil_dsl_init fresh_snw_rds_max fstor2tran " *
        "hksat_sf ignition_efficiency interception_fraction jmax25top_sf jmaxha jmaxhd " *
        "jmaxse_sf kc25_coef kcha ko25_coef koha kp25ratio lai_dl liq_canopy_storage_scalar " *
        "lmrha lmrhd lmrse luna_theta_cj lwtop_ann maximum_leaf_wetted_fraction minrelh mu " *
        "n_baseflow n_melt_coef ndays_off om_frac_sf pc pd perched_baseflow_scalar " *
        "phenology_soil_depth prh30 relhExp rho_max sand_pf scvng_fct_mlt_bcphi " *
        "scvng_fct_mlt_bcpho scvng_fct_mlt_dst1 scvng_fct_mlt_dst2 scvng_fct_mlt_dst3 " *
        "scvng_fct_mlt_dst4 scvng_fct_mlt_sf slopebeta slopemax snow_canopy_storage_scalar " *
        "snow5d_thresh_for_onset snowcan_unload_temp_fact snowcan_unload_wind_fact " *
        "snw_rds_min snw_rds_refrz soilpsi_off soilpsi_on ssi sucsat_sf tau_ref theta_ip " *
        "theta_psii tkd_clay tkd_om tkd_sand tkm_om tpu25ratio tpuha tpuhd tpuse_sf " *
        "upplim_destruct_metamorph vcmaxha vcmaxhd vcmaxse_sf watsat_sf wimp wind_min " *
        "wind_snowcompact_fact xdrdt z_dl zbedrock zbedrock_sf zglc zlnd zsno"))

    known_unwired = Set([
        # PPE perturbation knobs for cellsand/cellclay; 0 (no-op) in the shipped
        # file, perturbation branch deliberately not ported.
        "sand_pf", "clay_pf",
        # Per-instance roughness lengths (frictionvel) — wired in clm_run!, which
        # has `inst`. zlnd IS wired here (waterdiagbulk_params).
        "zsno", "zglc",
        # Per-column / per-instance snow-cover-fraction params — wired in clm_run!.
        "n_melt_coef", "accum_factor",
        # Fire params: no module-level CNFireParams instance (clm_driver! never runs
        # fire); harnesses populate theirs via readParams_CNFire!.
        "prh30", "ignition_efficiency",
        # zbedrock is read into _initvert_zbedrock (sentinel -1 = don't substitute).
        "zbedrock",
        # LUNA scalars, read by luna_read_params!.
        "enzyme_turnover_daily", "luna_theta_cj", "minrelh", "relhExp",
        # Photosynthesis scalars, read wholesale by readParams_photosynthesis!.
        "cp25_yr2000", "cpha", "fnps", "jmax25top_sf", "jmaxha", "jmaxhd",
        "jmaxse_sf", "kc25_coef", "kcha", "ko25_coef", "koha", "kp25ratio",
        "lmrha", "lmrhd", "lmrse", "theta_psii", "tpu25ratio", "tpuha", "tpuhd",
        "tpuse_sf", "vcmaxha", "vcmaxhd", "vcmaxse_sf",
        # a_coef/a_exp feed both flux structs; covered by the bareground/canopy asserts.
        "a_coef", "a_exp",
    ])
    covered = Set(name for (name, _, _) in _PARAM_CASES)
    unaccounted = setdiff(fortran_scalars, union(covered, known_unwired))
    @test isempty(unaccounted)
end
