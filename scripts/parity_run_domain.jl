#!/usr/bin/env julia
# ==========================================================================
# parity_run_domain.jl — domain-parameterized spun-up-vs-spun-up annual parity.
#
# Generalizes parity_run_spinup.jl (Bow) to any biome domain that has a Fortran
# CLM reference: initialize Julia from the Fortran spun-up restart, run the full
# reference year on observed forcing with the calibrated params, write daily
# history to paper/data/julia_clm_<domain>_<year>.nc. The plot script
# (paper/plot_parity_domain.py) compares it to that domain's Fortran h0.
#
#   DOMAIN=Stillwater julia +1.12 --project=. scripts/parity_run_domain.jl
#   DOMAIN=Aripuana   julia +1.12 --project=. scripts/parity_run_domain.jl
#   DOMAIN=Bow        julia +1.12 --project=. scripts/parity_run_domain.jl
#
# All three reference runs share the Bow lnd_in physics config (use_hydrstress +
# use_luna, Jackson1996 rooting, Vionnet2012 + wind-dependent snow density,
# baseflow_scalar/int_snow_max per lnd_in), so the only per-domain differences are
# the file paths, the year, dtime, and those two calibrated scalars.
# ==========================================================================

using NCDatasets, Dates, Printf, CLM
const run_clm! = getfield(CLM, Symbol("clm_run!"))
const DATA = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data"

# SNICAR optics/aging (shared); prefer the local copy, fall back to projects.
_snicar(f) = let d1 = "$DATA/installs/cesm-inputdata/lnd/clm2/snicardata",
                 d2 = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata"
    isfile(joinpath(d1, f)) ? joinpath(d1, f) : joinpath(d2, f)
end

# ---- Per-domain registry ----
const DOMAINS = Dict(
    # The CTSM reference is hourly: n8761 is 2003-01-01 01:00 immediately
    # following the prior 8760-step year. Using 1800 s here compares different
    # numerical protocols and creates artificial snow, canopy, and flux misses.
    "Bow" => (year = 2003, dtime = 3600, baseflow = 0.0022119554, int_snow = 3113.2227,
        caldir  = "$DATA/domain_Bow_at_Banff_lumped/optimization/CLM/dds_run_1/final_evaluation/settings/CLM/parameters",
        forcing = "$DATA/domain_Bow_at_Banff_lumped/data/forcing/CLM_input/clmforc.2003.nc",
        restart = "$DATA/clm_parity_run/Bow_at_Banff_lumped.clm2.r.2003-01-01-00000.nc",
        h0      = "$DATA/clm_parity_run/Bow_at_Banff_lumped.clm2.h0.2003-01-01-00000.nc"),
    "Stillwater" => (year = 2003, dtime = 3600, baseflow = 0.016035343, int_snow = 2000.0,
        caldir  = "$DATA/domain_Stillwater_Oklahoma/optimization/CLM/dds_clm_dds_calibration/final_evaluation/settings/CLM/parameters",
        forcing = "$DATA/domain_Stillwater_Oklahoma/data/forcing/CLM_input/clmforc.2003.nc",
        restart = "$DATA/domain_Stillwater_Oklahoma/optimization/CLM/dds_clm_dds_calibration/final_evaluation/Stillwater_Oklahoma.clm2.r.2003-01-01-00000.nc",
        h0      = "$DATA/domain_Stillwater_Oklahoma/optimization/CLM/dds_clm_dds_calibration/final_evaluation/Stillwater_Oklahoma.clm2.h0.2003-01-01-00000.nc"),
    "Aripuana" => (year = 2004, dtime = 3600, baseflow = 0.001, int_snow = 2000.0,
        caldir  = "$DATA/domain_Aripuana_Amazon/settings/CLM/parameters",
        forcing = "$DATA/domain_Aripuana_Amazon/data/forcing/CLM_input/clmforc.2004.nc",
        restart = "$DATA/domain_Aripuana_Amazon/simulations/clm_validation/CLM/Aripuana_Amazon.clm2.r.2004-01-01-00000.nc",
        h0      = "$DATA/domain_Aripuana_Amazon/simulations/clm_validation/CLM/Aripuana_Amazon.clm2.h0.2004-01-01-00000.nc"),
    # Boreal forest (Krycklan, Sweden). Fortran reference generated via SYMFLUENCE
    # run_model (spinup 2008-2012, then sim); parity year 2013 from the 2013-01-01
    # spun-up restart. The 2013 daily h0 records live in the 2012-12-30 file
    # (noleap interval-end stamping); the plot aligns by calendar date.
    "Krycklan" => (year = 2013, dtime = 3600, baseflow = 0.001, int_snow = 2000.0,
        caldir  = "$DATA/domain_Boreal_Krycklan_Sweden/settings/CLM/parameters",
        forcing = "$DATA/domain_Boreal_Krycklan_Sweden/data/forcing/CLM_input/clmforc.2013.nc",
        restart = "$DATA/domain_Boreal_Krycklan_Sweden/simulations/clm_boreal/CLM/Boreal_Krycklan_Sweden.clm2.r.2013-01-01-00000.nc",
        h0      = "$DATA/domain_Boreal_Krycklan_Sweden/simulations/clm_boreal/CLM/Boreal_Krycklan_Sweden.clm2.h0.2012-12-30-00000.nc"),
    # ---- Eval-prep biomes: Fortran references generated via SYMFLUENCE run_model
    # (lumped, spinup 2008-2012 → parity year 2013 in the .h0.2012-12-30 file;
    # Iceland spins 2015-2016 → year 2017 in .h0.2016-12-31). Same Bow lnd_in
    # physics (PHS+LUNA, non-aquifer, baseflow 0.001, int_snow 2000). Baltimore is
    # URBAN and Massa/Iceland are GLACIER — exotic-landunit robustness stress tests.
    "Tagus" => (year = 2013, dtime = 3600, baseflow = 0.001, int_snow = 2000.0,
        caldir  = "$DATA/domain_Mediterranean_Tagus_Spain/settings/CLM/parameters",
        forcing = "$DATA/domain_Mediterranean_Tagus_Spain/data/forcing/CLM_input/clmforc.2013.nc",
        restart = "$DATA/domain_Mediterranean_Tagus_Spain/simulations/clm_mediterranean/CLM/Mediterranean_Tagus_Spain.clm2.r.2013-01-01-00000.nc",
        h0      = "$DATA/domain_Mediterranean_Tagus_Spain/simulations/clm_mediterranean/CLM/Mediterranean_Tagus_Spain.clm2.h0.2012-12-30-00000.nc"),
    "Abisko" => (year = 2013, dtime = 3600, baseflow = 0.001, int_snow = 2000.0,
        caldir  = "$DATA/domain_Arctic_Abisko_Sweden/settings/CLM/parameters",
        forcing = "$DATA/domain_Arctic_Abisko_Sweden/data/forcing/CLM_input/clmforc.2013.nc",
        restart = "$DATA/domain_Arctic_Abisko_Sweden/simulations/clm_arctic/CLM/Arctic_Abisko_Sweden.clm2.r.2013-01-01-00000.nc",
        h0      = "$DATA/domain_Arctic_Abisko_Sweden/simulations/clm_arctic/CLM/Arctic_Abisko_Sweden.clm2.h0.2012-12-30-00000.nc"),
    "Massa" => (year = 2013, dtime = 3600, baseflow = 0.001, int_snow = 2000.0,
        caldir  = "$DATA/domain_Alps_Massa_Aletsch_CH/settings/CLM/parameters",
        forcing = "$DATA/domain_Alps_Massa_Aletsch_CH/data/forcing/CLM_input/clmforc.2013.nc",
        restart = "$DATA/domain_Alps_Massa_Aletsch_CH/simulations/clm_alpineglacier/CLM/Alps_Massa_Aletsch_CH.clm2.r.2013-01-01-00000.nc",
        h0      = "$DATA/domain_Alps_Massa_Aletsch_CH/simulations/clm_alpineglacier/CLM/Alps_Massa_Aletsch_CH.clm2.h0.2012-12-30-00000.nc"),
    "Baltimore" => (year = 2013, dtime = 3600, baseflow = 0.001, int_snow = 2000.0,
        caldir  = "$DATA/domain_Urban_DeadRun_Baltimore/settings/CLM/parameters",
        forcing = "$DATA/domain_Urban_DeadRun_Baltimore/data/forcing/CLM_input/clmforc.2013.nc",
        restart = "$DATA/domain_Urban_DeadRun_Baltimore/simulations/clm_urban/CLM/Urban_DeadRun_Baltimore.clm2.r.2013-01-01-00000.nc",
        h0      = "$DATA/domain_Urban_DeadRun_Baltimore/simulations/clm_urban/CLM/Urban_DeadRun_Baltimore.clm2.h0.2012-12-30-00000.nc"),
    # ---- New-biome expansion (2026-07-03): fresh SYMFLUENCE references generated
    # from scratch (define_domain → ERA5 forcing → surfdata → run_model via lldb
    # to bypass the ad-hoc-codesign SIGTRAP). Spinup 2014-2016 (HubbardBrook) or
    # Jul-Dec 2016 (others) → parity year 2017 in the .h0.2016-12-31 file.
    "HubbardBrook" => (year = 2017, dtime = 3600, baseflow = 0.001, int_snow = 2000.0,
        caldir  = "$DATA/domain_Temperate_HubbardBrook_USA/settings/CLM/parameters",
        forcing = "$DATA/domain_Temperate_HubbardBrook_USA/data/forcing/CLM_input/clmforc.2017.nc",
        restart = "$DATA/domain_Temperate_HubbardBrook_USA/simulations/clm_tempforest/CLM/Temperate_HubbardBrook_USA.clm2.r.2017-01-01-00000.nc",
        h0      = "$DATA/domain_Temperate_HubbardBrook_USA/simulations/clm_tempforest/CLM/Temperate_HubbardBrook_USA.clm2.h0.2016-12-31-00000.nc"),
    "Maritime" => (year = 2017, dtime = 3600, baseflow = 0.001, int_snow = 2000.0,
        caldir  = "$DATA/domain_Maritime_HJAndrews_USA/settings/CLM/parameters",
        forcing = "$DATA/domain_Maritime_HJAndrews_USA/data/forcing/CLM_input/clmforc.2017.nc",
        restart = "$DATA/domain_Maritime_HJAndrews_USA/simulations/clm_maritime/CLM/Maritime_HJAndrews_USA.clm2.r.2017-01-01-00000.nc",
        h0      = "$DATA/domain_Maritime_HJAndrews_USA/simulations/clm_maritime/CLM/Maritime_HJAndrews_USA.clm2.h0.2016-12-31-00000.nc"),
    # Hot desert (Walnut Gulch, Arizona) and tropical savanna (Donga, Benin):
    # Fortran references generated via SYMFLUENCE run_model, spinup 2014-2016 →
    # parity year 2017 in the .h0.2016-12-31 file. Same Bow lnd_in physics.
    "Desert" => (year = 2017, dtime = 3600, baseflow = 0.001, int_snow = 2000.0,
        caldir  = "$DATA/domain_Desert_WalnutGulch_USA/settings/CLM/parameters",
        forcing = "$DATA/domain_Desert_WalnutGulch_USA/data/forcing/CLM_input/clmforc.2017.nc",
        restart = "$DATA/domain_Desert_WalnutGulch_USA/simulations/clm_desert/CLM/Desert_WalnutGulch_USA.clm2.r.2017-01-01-00000.nc",
        h0      = "$DATA/domain_Desert_WalnutGulch_USA/simulations/clm_desert/CLM/Desert_WalnutGulch_USA.clm2.h0.2016-12-31-00000.nc"),
    "Savanna" => (year = 2017, dtime = 3600, baseflow = 0.001, int_snow = 2000.0,
        caldir  = "$DATA/domain_Savanna_Donga_Benin/settings/CLM/parameters",
        forcing = "$DATA/domain_Savanna_Donga_Benin/data/forcing/CLM_input/clmforc.2017.nc",
        restart = "$DATA/domain_Savanna_Donga_Benin/simulations/clm_savanna/CLM/Savanna_Donga_Benin.clm2.r.2017-01-01-00000.nc",
        h0      = "$DATA/domain_Savanna_Donga_Benin/simulations/clm_savanna/CLM/Savanna_Donga_Benin.clm2.h0.2016-12-31-00000.nc"),
    # Biome-gap expansion (2026-07-04): larch permafrost taiga (needleleaf-deciduous
    # boreal PFT), temperate peat bog (wetland hydrology), cold continental steppe.
    "Larch" => (year = 2017, dtime = 3600, baseflow = 0.001, int_snow = 2000.0,
        caldir  = "$DATA/domain_Yakutia_Larch_Russia/settings/CLM/parameters",
        forcing = "$DATA/domain_Yakutia_Larch_Russia/data/forcing/CLM_input/clmforc.2017.nc",
        restart = "$DATA/domain_Yakutia_Larch_Russia/simulations/clm_larch/CLM/Yakutia_Larch_Russia.clm2.r.2017-01-01-00000.nc",
        h0      = "$DATA/domain_Yakutia_Larch_Russia/simulations/clm_larch/CLM/Yakutia_Larch_Russia.clm2.h0.2016-12-31-00000.nc"),
    "Peatland" => (year = 2017, dtime = 3600, baseflow = 0.001, int_snow = 2000.0,
        caldir  = "$DATA/domain_Peatland_MerBleue_Canada/settings/CLM/parameters",
        forcing = "$DATA/domain_Peatland_MerBleue_Canada/data/forcing/CLM_input/clmforc.2017.nc",
        restart = "$DATA/domain_Peatland_MerBleue_Canada/simulations/clm_peatland/CLM/Peatland_MerBleue_Canada.clm2.r.2017-01-01-00000.nc",
        h0      = "$DATA/domain_Peatland_MerBleue_Canada/simulations/clm_peatland/CLM/Peatland_MerBleue_Canada.clm2.h0.2016-12-31-00000.nc"),
    "Steppe" => (year = 2017, dtime = 3600, baseflow = 0.001, int_snow = 2000.0,
        caldir  = "$DATA/domain_Steppe_Kherlen_Mongolia/settings/CLM/parameters",
        forcing = "$DATA/domain_Steppe_Kherlen_Mongolia/data/forcing/CLM_input/clmforc.2017.nc",
        restart = "$DATA/domain_Steppe_Kherlen_Mongolia/simulations/clm_steppe/CLM/Steppe_Kherlen_Mongolia.clm2.r.2017-01-01-00000.nc",
        h0      = "$DATA/domain_Steppe_Kherlen_Mongolia/simulations/clm_steppe/CLM/Steppe_Kherlen_Mongolia.clm2.h0.2016-12-31-00000.nc"),
    # PFT/climate + subsystem expansion (2026-07-04): temperate broadleaf-evergreen
    # (Eucalyptus), broadleaf-deciduous boreal (aspen), tropical montane páramo,
    # and cropland (crop PFTs).
    "Eucalyptus" => (year = 2017, dtime = 3600, baseflow = 0.001, int_snow = 2000.0,
        caldir  = "$DATA/domain_Eucalyptus_Tumbarumba_Australia/settings/CLM/parameters",
        forcing = "$DATA/domain_Eucalyptus_Tumbarumba_Australia/data/forcing/CLM_input/clmforc.2017.nc",
        restart = "$DATA/domain_Eucalyptus_Tumbarumba_Australia/simulations/clm_eucalyptus/CLM/Eucalyptus_Tumbarumba_Australia.clm2.r.2017-01-01-00000.nc",
        h0      = "$DATA/domain_Eucalyptus_Tumbarumba_Australia/simulations/clm_eucalyptus/CLM/Eucalyptus_Tumbarumba_Australia.clm2.h0.2016-12-31-00000.nc"),
    "Aspen" => (year = 2017, dtime = 3600, baseflow = 0.001, int_snow = 2000.0,
        caldir  = "$DATA/domain_Aspen_BOREAS_Canada/settings/CLM/parameters",
        forcing = "$DATA/domain_Aspen_BOREAS_Canada/data/forcing/CLM_input/clmforc.2017.nc",
        restart = "$DATA/domain_Aspen_BOREAS_Canada/simulations/clm_aspen/CLM/Aspen_BOREAS_Canada.clm2.r.2017-01-01-00000.nc",
        h0      = "$DATA/domain_Aspen_BOREAS_Canada/simulations/clm_aspen/CLM/Aspen_BOREAS_Canada.clm2.h0.2016-12-31-00000.nc"),
    "Paramo" => (year = 2017, dtime = 3600, baseflow = 0.001, int_snow = 2000.0,
        caldir  = "$DATA/domain_Paramo_Antisana_Ecuador/settings/CLM/parameters",
        forcing = "$DATA/domain_Paramo_Antisana_Ecuador/data/forcing/CLM_input/clmforc.2017.nc",
        restart = "$DATA/domain_Paramo_Antisana_Ecuador/simulations/clm_paramo/CLM/Paramo_Antisana_Ecuador.clm2.r.2017-01-01-00000.nc",
        h0      = "$DATA/domain_Paramo_Antisana_Ecuador/simulations/clm_paramo/CLM/Paramo_Antisana_Ecuador.clm2.h0.2016-12-31-00000.nc"),
    "Cropland" => (year = 2017, dtime = 3600, baseflow = 0.001, int_snow = 2000.0,
        caldir  = "$DATA/domain_Cropland_Mead_USA/settings/CLM/parameters",
        forcing = "$DATA/domain_Cropland_Mead_USA/data/forcing/CLM_input/clmforc.2017.nc",
        restart = "$DATA/domain_Cropland_Mead_USA/simulations/clm_crop/CLM/Cropland_Mead_USA.clm2.r.2017-01-01-00000.nc",
        h0      = "$DATA/domain_Cropland_Mead_USA/simulations/clm_crop/CLM/Cropland_Mead_USA.clm2.h0.2016-12-31-00000.nc"),
    "Iceland" => (year = 2017, dtime = 3600, baseflow = 0.001, int_snow = 2000.0,
        caldir  = "$DATA/domain_Iceland_Jokulsa_Fjollum/settings/CLM/parameters",
        forcing = "$DATA/domain_Iceland_Jokulsa_Fjollum/data/forcing/CLM_input/clmforc.2017.nc",
        restart = "$DATA/domain_Iceland_Jokulsa_Fjollum/simulations/clm_glacier/CLM/Iceland_Jokulsa_Fjollum.clm2.r.2017-01-01-00000.nc",
        h0      = "$DATA/domain_Iceland_Jokulsa_Fjollum/simulations/clm_glacier/CLM/Iceland_Jokulsa_Fjollum.clm2.h0.2016-12-31-00000.nc"),
)

const DOM = get(ENV, "DOMAIN", "Bow")
haskey(DOMAINS, DOM) || error("Unknown DOMAIN=$DOM; choose one of $(collect(keys(DOMAINS)))")
cfg = DOMAINS[DOM]
yr  = cfg.year

fsurdat   = joinpath(cfg.caldir, "surfdata_clm.nc")
paramfile = joinpath(cfg.caldir, "clm5_params.nc")
outdir    = abspath(joinpath(@__DIR__, "..", "paper", "data"))
fhistory  = joinpath(outdir, "julia_clm_$(lowercase(DOM))_phs_$(yr).nc")
isdir(outdir) || mkpath(outdir)

println("="^64)
println("  CLM.jl parity run — $DOM, $(yr) (PHS+LUNA, spun-up-vs-spun-up)")
println("="^64)
for (lbl, f) in (("surfdata", fsurdat), ("params", paramfile), ("forcing", cfg.forcing),
                 ("restart", cfg.restart), ("fortran h0", cfg.h0))
    println("  ", isfile(f) ? "OK  " : "MISS", "  $lbl: $f")
end
println()

# All reference domains use the Bow lnd_in config: Jackson1996 rooting + Vionnet2012
# wind-dependent snow density (matching their datm/lnd_in).
CLM.rooting_profile_config.rooting_profile_method_water  = CLM.JACKSON_1996_ROOT
CLM.rooting_profile_config.rooting_profile_method_carbon = CLM.JACKSON_1996_ROOT
CLM.snow_hydrology_set_control_for_testing!(;
    wind_dep_snow_density = true,
    overburden_compaction_method = CLM.OVERBURDEN_COMPACTION_VIONNET2012)

t0 = time()
inst = run_clm!(;
    fsurdat = fsurdat, paramfile = paramfile, fforcing = cfg.forcing,
    fhistory = fhistory,
    start_date = DateTime(yr, 1, 1), end_date = DateTime(yr + 1, 1, 1),
    dtime = cfg.dtime, use_cn = false, verbose = true,
    use_aquifer_layer = get(cfg, :aquifer, false), use_hydrstress = true, use_luna = true,
    h2osfcflag = 1,  # CLM5 default (surface water active); matters for wet sites
    baseflow_scalar = cfg.baseflow, int_snow_max = cfg.int_snow,
    ffortran_restart = cfg.restart,
    interp_forcing = true,
    fsnowoptics = isfile(_snicar("snicar_optics_5bnd_c013122.nc")) ? _snicar("snicar_optics_5bnd_c013122.nc") : "",
    fsnowaging  = isfile(_snicar("snicar_drdt_bst_fit_60_c070416.nc")) ? _snicar("snicar_drdt_bst_fit_60_c070416.nc") : "")
@printf("\n  Done: %s in %.1f s -> %s\n", DOM, time() - t0, fhistory)
