# =============================================================================
# CN/BGC ANNUAL Julia↔Fortran parity — Bow, use_cn=true, spun-up IC.
#
# The first END-TO-END annual CN parity axis (the existing CN checks are single-
# step pdump probes + a 28-step drift). Injects the Fortran BGC spun-up restart
# (2202-01-01, spinup_state=0) ONCE, then free-runs a full year and compares the
# Julia CN/BGC carbon-nitrogen state + fluxes against the Fortran .h0 (monthly
# means), alongside the coupled energy/water.
#
# Fortran reference (generated 2026-06-19, branch from bgc_spunup_final, forcing
# clmforc.2003 via year_align, daily-avg 365 records):
#   /Users/.../SYMFLUENCE_data/clm_cn_run/Bow_at_Banff_lumped.clm2.h0.2202-01-02-00000.nc
# h0 vars are (time,[lev,]lndgrid) → NCDatasets returns them reversed → f[1, day].
# The model year is 2202 but the forcing is clmforc.2003 (year_align), so we drive
# the Julia run with 2003 forcing while the clock runs 2202 (date mapped per step).
#
# Config matches the cn_run lnd_in: use_cn/use_fun/use_flexiblecn/use_nitrif_denitrif/
# use_luna/use_hydrstress = true, vcmax_opt=3, CENTURYKoven2013, nofire, dtime=3600.
#
# Usage: julia +1.12 --project=. scripts/fortran_parity_cn_annual.jl [NDAYS]
#        CN_PROBE=1 julia ... scripts/fortran_parity_cn_annual.jl 1   # step-1 NaN probe
#
# STATUS (2026-06-20): scaffold complete + the comparison framework works, but it
# surfaced a BLOCKER on the first end-to-end use_cn=true annual run: the soil
# hydrology NaNs at STEP 1 from the BGC restart. Localized (CN_PROBE):
#   * the soil water IS injected finite (h2osoi_liq soil[1:4]=[1.72,3.42,5.26,7.03],
#     wa/zwt finite); watsat/bsw/sucsat finite; vegwp seeded finite.
#   * after one clm_drv! step (night: coszen=0, sabv=0 — NOT the daytime canopy-
#     albedo NaN) h2osoi_liq, smp_l, t_veg, t_grnd all go NaN.
#   * NOT build_bow_inst's light setup (full clm_initialize! setup here, same NaN);
#     NOT the pedotransfer params; use_cn=false runs the same winter soil finite
#     (verify_vs_fortran) -> it is a use_cn-path soil-water NaN.
# The winter CN pools (SMINN/TLAI/LEAFC) still match Fortran because they are
# ~static in winter, so that is a weak check; the dynamic/energy parity is blocked.
# ROOT CAUSE (2026-06-20): the step-1 NaN is the PHS (use_hydrstress=TRUE) path on
# FROZEN winter soil from a bare .clm2.r restart — NOT use_cn. CONFIRMED: CN_NOPHS=1
# (use_hydrstress=false) runs FINITE. The .clm2.r restart has no derived hydrology
# fields, so they start NaN; the PHS plant-water solve in canopy_fluxes (which runs
# before HydrologyNoDrainage recomputes them) reads them at step 1. Seeding smp_l +
# vegwp here + the driver's compute_h2osoi_liqvol! are not enough — the residual NaN
# is inside the PHS solve itself (k_soil_root / vegwp Newton) on near-frozen soil
# (the 28-step CN drift only ran summer/unfrozen, so it never hit it). FIX (deep,
# deferred — touches the PHS solve, which the reverse-AD mirror also covers):
# initialize the PHS-read derived hydrology state from the restart + make the PHS
# k_soil_root/Newton robust to frozen (near-zero-liquid) layers. MEANWHILE CN_NOPHS=1
# gives a finite (PHS-off, so not exact-parity on water-stress-sensitive GPP/veg) run
# for sanity-checking the slowly-varying soil C/N pools.
#
# PROGRESS (2026-06-20): the PHS step-1 NaN is a MULTI-LAYER cascade — 4 layers fixed,
# still not finite. CORRECTED DIAGNOSIS: it is NOT "frozen" soil — the probe shows
# liqvol~0.086 (DRY, not frozen) at NIGHT (coszen=0). Dry soil → tiny hk_l (~3e-14) →
# k_soil_root floored to ~1e-16 → degenerate PHS solve driving vegwp to the extreme
# (~-4e7 mm). Layers fixed (all general robustness, suite-green):
#   (1) eff_porosity init from the restart (fortran_restart.jl) — was NaN at step 1.
#   (2)+(3) getvegwp!/_getvegwp sum_ksr guard ==0 → THRESHOLD <=1e-12 (the smooth_max
#       floor leaves k_soil_root~1e-16, so sum_ksr~1e-17≠0 divided → huge xroot).
#   (4) plc/_plc/d1plc/_d1plc now return 1/0 for x>=0 (photosynthesis.jl) — water
#       potential is always <=0 physically, and (x/psi50)<0 ^ (non-integer ck) = NaN;
#       the degenerate solve drives vegwp >=0 transiently. (Latent bug; Fortran never
#       hits x>=0 so parity holds.) vegwp now finite, but bsun/bsha STILL NaN.
# REMAINING (open frontier): bsun/bsha NaN in the NIGHT branch of calcstress! +
# smp_l(c,1)=NaN feeding it (a separate seed gap). ROOT ISSUE: the restart seed does
# not provide ALL the derived hydrology state the PHS-on path reads at step 1 (smp_l,
# the dry-soil k_soil_root regime), so the night PHS solve goes degenerate. PROPER FIX
# = either seed ALL PHS-read derived state consistently from the restart, OR add a
# clean dry/degenerate PHS bypass at the calcstress entry — both deeper than these
# point patches, and in the reverse-AD mirror's code (coordinate). PAUSED here.
# CN_NOPHS=1 = working PHS-off fallback (see CN_NOPHS findings below).
#
# RAW-POOL DIAGNOSTICS FIXED (2026-06-20): the "NaN C/N pools" was a pure DIAGNOSTIC
# artifact, not physics. The _col summary fields (totsomc_col/totvegc_col/totsomn_col)
# and gpp_patch are NOT refreshed in this CN driver path — their summary routines
# (soil_bgc_carbon_state_summary! / cnveg_carbon_flux_summary!) are unported
# placeholders (clm_driver.jl ~952), so they keep their restart-init NaN/0. This
# harness now aggregates the actual PROGNOSTIC pools (see RAW-POOL block below):
#   TOTSOMC/N = sum_j sum_{is_soil pools 4-6} decomp_*pools_vr[c,j,l]*dz_decomp[j]
#   TOTVEGC   = patch sum of all veg C pools (displayed+storage+xfer), wtgcell->grc
#   GPP       = psnsun_to_cpool + psnshade_to_cpool (real photosynthate flux). NOTE: as
#               of the summary-wiring fix, the model's own gpp_patch/npp_patch ARE now
#               populated (cnveg_carbon_flux_summary! wired into cn_driver_summarize_
#               fluxes!), so the GPP field reads gpp_patch directly (== the psn streams;
#               GPP_psn cross-checks) and NPP reads npp_patch (was stale 0 → now a real
#               seasonal value: -2.7e-7 winter resp. → +6.6e-6 Jun, vs Fortran 9.1e-6).
#               The soil-C/N _col STATE summaries are still stubs → TOTSOMC/TOTVEGC/
#               TOTSOMN still use raw-pool aggregation here (follow-up: wire those too).
# RESULT (3-day smoke vs Fortran h0, Jan): TOTSOMC rel 3.3e-8, TOTVEGC 5.9e-8,
# TOTECOSYSC 5.2e-7, TOTSOMN 1.7e-8, SMINN 1.4e-5, TLAI 8e-6, LEAFC 8.6e-4, TG 7.6e-5
# — i.e. the CN state was at PARITY all along; only the diagnostics were stale. The
# Jan match also confirms the forcing is shared (clmforc.2003), so the day-of-year
# comparison is valid. GPP in Jan is ~0 in both (winter); the seasonal leaf-out /
# summer-GPP trajectory is the real open test (full-year run).
#
# LEAF-OUT FIXED (2026-06-20): two unrelated bugs blocked the season_decid grass (p3,
# ivt=12, wt=0.35; p2 ivt=1 needleleaf-evergreen is correctly constant-LAI) from
# leafing out, so GPP/TLAI/TOTVEGC were flat all year. Found via the CN_PHENO probe
# (onset_gdd stuck at 0.0 even at peak summer):
#  (1) phenology.jl: the season_decid/stress_decid GDD-onset read soil temperature as
#      t_soisno_col[c, soil_layer] with soil_layer=1, but t_soisno_col is SNOW-PADDED
#      ([nlevsno + nlevgrnd]) so index 1 is the top snow slot (frozen year-round) →
#      onset_gdd never accumulated → onset never fired → dormant_flag stuck at 1. Fixed
#      by offsetting the t_soisno read by nlevsno (ts_layer = nlevsno + soil_layer) in
#      both kernels, keeping soilpsi_col (soil-only) at soil_layer. (soilpsi was right;
#      only t_soisno needed the offset — that mismatch was the giveaway.)
#  (2) cn_veg_struct_update! (leafc -> tlai/tsai/htop/elai) was DEFINED but never CALLED
#      in the driver path → even after onset fired and leafc grew, tlai stayed frozen at
#      the restart value, so the canopy never responded. NOW WIRED into the global driver
#      (clm_driver.jl, use_cn && doalb). Two fixes were needed to wire it cleanly:
#      (a) pass npcropmin=config.npcropmin — the CLM.npcropmin default is 0, so the crop
#          branch (ivt>=npcropmin) caught grasses (ivt=12) and gave tsai=0.2*tlai instead
#          of the 0.5 floor, breaking freewins surface-albedo parity; (b) skip patches
#          with non-finite leaf C in cn_veg_struct_update! (cold-start CN pools are NaN
#          off-Bow → would poison tlai; multisite robustness). No harness call needed now.
# RESULT (CN_NOPHS full year vs Fortran h0): TLAI flat 0.031 -> grows to 0.109 (Jul,
# vs Fortran 0.111); TOTVEGC Jul rel 1.0e-1 -> 2.2e-2; TOTECOSYSC 1.5e-3 -> 4.3e-4;
# SMINN 4.3e-2 -> 1.4e-2. The CN axis is now at good seasonal parity (PHS-off).
#  (3) GDD soil layer: phenology_soil_depth=0.08m maps (find_soil_layer_containing_depth)
#      to soil LAYER 3 (CLM5 interfaces 0.02/0.06/0.12 -> 0.08 falls in layer 3, node
#      ~0.09m), but cn_phenology_init!'s default stub returned layer 1 (~0.01m, too
#      shallow/responsive) -> onset too early, TLAI/GPP overshoot. Fixed by passing a
#      Fortran-matching find_soil_layer_fn from cn_driver.jl. TLAI parity improved ~10x:
#      Apr/May rel 1.1e-2/1.6e-2 -> 1.6e-3/1.7e-3 (Jun 1.8e-3, matches Fortran 0.111).
#      GPP overshoot ~1.3x -> ~0.85-0.92x (May 5.1e-6 vs 6.0e-6) — onset now tracks
#      Fortran's timing. Probe: CN_PHENO=1 CN_NOPHS=1 julia ... 365.
# =============================================================================
include(joinpath(@__DIR__, "fortran_parity_common.jl"))
using Statistics

const CN_DIR  = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_cn_run"
const CN_RST  = joinpath(CN_DIR, "Bow_at_Banff_lumped.clm2.r.2202-01-01-00000.nc")
const CN_H0   = joinpath(CN_DIR, "Bow_at_Banff_lumped.clm2.h0.2202-01-02-00000.nc")
const CN_FFORC = replace(FFORCING, "clmforc.2003.nc" => "clmforc.2003.nc")  # 2003 forcing
const FORCE_YEAR = 2003
const _read_fortran_restart! = getfield(CLM, Symbol("read_fortran_restart!"))

# model timestamp → forcing-file timestamp (year replaced; both 2202 & 2003 are 365-day)
_force_date(d) = DateTime(FORCE_YEAR, Dates.month(d), Dates.day(d), Dates.hour(d),
                          Dates.minute(d), Dates.second(d))
_fv(x) = ismissing(x) ? NaN : Float64(x)

function run_cn_annual(; ndays::Int = 365)
    isfile(CN_RST) || error("CN restart not found: $CN_RST")
    isfile(CN_H0)  || error("CN h0 reference not found: $CN_H0")
    # Run at the FORCING year (clock=2003 matches clmforc.2003); the CN restart
    # supplies the spun-up 2202 state; compare by day-of-year to the 2202 h0.
    # FULL clm_initialize! setup (mirrors clm_run.jl) — build_bow_inst is a light
    # injection setup that leaves smp_l/derived fields NaN, which the soil hydrology
    # then propagates into h2osoi/t_grnd over a long run from a bare restart.
    start_date = DateTime(FORCE_YEAR, 1, 1)
    CLM.rooting_profile_config.rooting_profile_method_water  = CLM.JACKSON_1996_ROOT
    CLM.rooting_profile_config.rooting_profile_method_carbon = CLM.JACKSON_1996_ROOT
    snowopt = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_optics_5bnd_c013122.nc"
    snowage = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_drdt_bst_fit_60_c070416.nc"
    (inst, bounds, filt, tm) = CLM.clm_initialize!(; fsurdat = FSURDAT, paramfile = FPARAM,
        start_date = start_date, dtime = 3600, use_cn = true, use_luna = true,
        use_bedrock = true, use_aquifer_layer = false, h2osfcflag = 0,
        fsnowoptics = isfile(snowopt) ? snowopt : "", fsnowaging = isfile(snowage) ? snowage : "",
        int_snow_max = 3113.2227)
    ng, nc, np = bounds.endg, bounds.endc, bounds.endp
    CLM.atm2lnd_read_namelist!(inst.atm2lnd; repartition_rain_snow = true, lapse_rate = 0.006,
        lapse_rate_longwave = 0.032, precip_repartition_nonglc_all_snow_t = 0.0,
        precip_repartition_nonglc_all_rain_t = 2.0, precip_repartition_glc_all_snow_t = -2.0,
        precip_repartition_glc_all_rain_t = 0.0)
    CLM.init_soil_hydrology_config(baseflow_scalar = 0.0022119554)
    let ds_p = NCDataset(FPARAM, "r")
        scf = inst.scf_method
        haskey(ds_p, "n_melt_coef") && (nmc = Float64(ds_p["n_melt_coef"][1]);
            for c in 1:nc; scf.n_melt[c] = nmc / max(10.0, inst.column.topo_std[c]); end)
        haskey(ds_p, "accum_factor") && (scf.accum_factor = Float64(ds_p["accum_factor"][1]))
        haskey(ds_p, "SNOW_DENSITY_MAX") && (CLM.snowhydrology_params.rho_max = Float64(ds_p["SNOW_DENSITY_MAX"][1]))
        haskey(ds_p, "SNOW_DENSITY_MIN") && (CLM.snowhydrology_params.rho_min = Float64(ds_p["SNOW_DENSITY_MIN"][1]))
        haskey(ds_p, "fresh_snw_rds_max") && (CLM.snowhydrology_params.snw_rds_min = Float64(ds_p["fresh_snw_rds_max"][1]))
        haskey(ds_p, "SNO_Z0MV") && (inst.frictionvel.zsno = Float64(ds_p["SNO_Z0MV"][1]))
        if haskey(ds_p, "snw_aging_bst"); CLM.snicar_params.xdrdt = Float64(ds_p["snw_aging_bst"][1])
        elseif haskey(ds_p, "xdrdt"); CLM.snicar_params.xdrdt = Float64(ds_p["xdrdt"][1]); end
        haskey(ds_p, "pc") && (pcv = Float64(ds_p["pc"][1]);
            pcv > 0 && for c in 1:nc; inst.soilhydrology.hkdepth_col[c] = 1.0 / pcv; end)
        close(ds_p)
    end
    # LUNA fields so the restart fills them; PHS needs finite vcmax/jmax seeds.
    if isempty(inst.photosyns.vcmx25_z_patch)
        inst.photosyns.vcmx25_z_patch = fill(30.0, np, CLM.NLEVCAN)
        inst.photosyns.jmx25_z_patch  = fill(60.0, np, CLM.NLEVCAN)
    end
    _read_fortran_restart!(CN_RST, inst, bounds)
    # Seed PHS leaf water potential from the restart (read_fortran_restart! leaves
    # vegwp NaN → the use_hydrstress Newton would propagate NaN into t_veg/t_grnd).
    let ds = NCDataset(CN_RST, "r")
        if haskey(ds, "vegwp")
            vw = ds["vegwp"][:, :]
            for pd in 1:size(vw, 2), seg in 1:size(vw, 1)
                seg <= size(inst.canopystate.vegwp_patch, 2) && pd <= np &&
                    (inst.canopystate.vegwp_patch[pd, seg] = Float64(vw[seg, pd]))
            end
        else
            fill!(inst.canopystate.vegwp_patch, -2.0e5)  # finite default [mm]
        end
        close(ds)
    end
    # Seed soil matric potential from the injected h2osoi_vol. The .clm2.r restart
    # has no smp_l (it's derived), so it starts NaN; with use_hydrstress the PHS
    # plant sink (canopy_fluxes) reads smp_l at step 1 BEFORE the hydrology recomputes
    # it → k_soil_root/qflx_rootsoi NaN → soil water NaN. (The 28-step CN drift
    # injected a pdump that included smp_l, so it never hit this.)
    CLM.update_smp_l!(inst.soilstate.smp_l_col, inst.water.waterstatebulk_inst.ws.h2osoi_vol_col,
        inst.soilstate.watsat_col, inst.soilstate.sucsat_col, inst.soilstate.bsw_col,
        inst.soilstate.smpmin_col, filt.nolakec, 1:nc, CLM.varpar.nlevgrnd)
    if get(ENV, "CN_PROBE", "") != ""
        cc = findfirst(filt.nolakec)
        _chk(nm, x) = @printf("    %-16s %s\n", nm, isfinite(x) ? @sprintf("%.4g", x) : "NaN")
        println("  [post-inject] soil col=$cc:")
        _chk("t_grnd", inst.temperature.t_grnd_col[cc])
        _chk("t_veg(p2)", inst.temperature.t_veg_patch[2])
        _chk("vegwp(p2,1)", inst.canopystate.vegwp_patch[2, 1])
        _chk("h2osoi_liq(c,1)", inst.water.waterstatebulk_inst.ws.h2osoi_liq_col[cc, 1])
        _chk("smp_l(c,1)", inst.soilstate.smp_l_col[cc, 1])
        _chk("albgrd(c,1)", inst.surfalb.albgrd_col[cc, 1])
        _chk("sabv(p2)", inst.solarabs.sabv_patch[2])
        ws = inst.water.waterstatebulk_inst.ws; nsno = CLM.varpar.nlevsno
        @printf("    h2osoi_liq SOIL[1:4] %s\n", string(round.(ws.h2osoi_liq_col[cc, (nsno+1):(nsno+4)], digits=2)))
        @printf("    h2osoi_ice SOIL[1:4] %s\n", string(round.(ws.h2osoi_ice_col[cc, (nsno+1):(nsno+4)], digits=2)))
        @printf("    wa=%s zwt=%s\n", string(ws.wa_col[cc]), string(inst.soilhydrology.zwt_col[cc]))
        @printf("    vegwp(p2,1:4)   %s\n", string(inst.canopystate.vegwp_patch[2, 1:4]))
        ss = inst.soilstate
        @printf("    watsat[1:4] %s\n", string(round.(ss.watsat_col[cc, 1:4], digits=3)))
        @printf("    bsw[1:4]    %s\n", string(round.(ss.bsw_col[cc, 1:4], digits=3)))
        @printf("    sucsat[1:4] %s\n", string(round.(ss.sucsat_col[cc, 1:4], digits=2)))
        @printf("    watsat[18:21] %s\n", string(round.(ss.watsat_col[cc, 18:21], digits=3)))
        @printf("    eff_poros[1:4] %s\n", string(round.(ss.eff_porosity_col[cc, 1:4], digits=3)))
        @printf("    rootfr(p2)[1:4] %s  k_soil_root(p2)[1:4] %s\n",
            string(round.(ss.rootfr_patch[2, 1:4], digits=3)), string(ss.k_soil_root_patch[2, 1:4]))
    end

    # Phenology-onset probe: which patches are deciduous + seeded onset state.
    if get(ENV, "CN_PHENO", "") != "" && hasproperty(inst.bgc_vegetation, :cnveg_state_inst)
        cvs0 = inst.bgc_vegetation.cnveg_state_inst
        pc = CLM.pftcon
        println("  [pheno-init] per-patch phenology type + seeded onset state:")
        for p in 1:np
            inst.patch.active[p] || continue
            it = inst.patch.itype[p]; ix = it + 1   # pft type → 1-based pftcon index
            eg = ix <= length(pc.evergreen) ? pc.evergreen[ix] : NaN
            sd = ix <= length(pc.season_decid) ? pc.season_decid[ix] : NaN
            st = ix <= length(pc.stress_decid) ? pc.stress_decid[ix] : NaN
            @printf("    p%-2d ivt=%-2d wt=%.3f evgrn=%g seas=%g stress=%g | dormant=%g gddflag=%g gdd=%.1f days_act=%.1f swi=%.1f onflag=%g\n",
                p, it, inst.patch.wtgcell[p], eg, sd, st,
                cvs0.dormant_flag_patch[p], cvs0.onset_gddflag_patch[p], cvs0.onset_gdd_patch[p],
                cvs0.days_active_patch[p], cvs0.onset_swi_patch[p], cvs0.onset_flag_patch[p])
        end
    end

    _phs = get(ENV, "CN_NOPHS", "") == ""   # CN_NOPHS=1 → diagnostic: disable PHS
    config  = CLM.CLMDriverConfig(use_cn = true, use_aquifer_layer = false,
                                  use_hydrstress = _phs, use_luna = true)
    filt_ia = CLM.clump_filter_inactive_and_active
    fr = CLM.ForcingReader(); CLM.forcing_reader_init!(fr, CN_FFORC)
    tf = replace(CN_FFORC, r"clmforc\.[^/]*\.nc$" => "topo_forcing.nc")
    if isfile(tf)
        dt = NCDataset(tf, "r")
        haskey(dt, "TOPO") && (ft = Float64(dt["TOPO"][1]);
            for g in 1:ng; inst.atm2lnd.forc_topo_grc[g] = ft; end;
            for c in 1:nc; inst.topo.topo_col[c] = ft; end)
        close(dt)
    end

    c_soil = findfirst(filt.nolakec)
    nl = CLM.varpar.nlevdecomp
    dzd = CLM.dzsoi_decomp[]
    scs() = inst.soilbiogeochem_carbonstate
    sns() = inst.soilbiogeochem_nitrogenstate
    ccs() = inst.bgc_vegetation.cnveg_carbonstate_inst
    ccf() = inst.bgc_vegetation.cnveg_carbonflux_inst
    # gridcell aggregate of a patch field (wtgcell-weighted sum over active patches)
    function p2g(arr)
        s = 0.0
        @inbounds for p in 1:np
            inst.patch.active[p] || continue
            v = arr[p]; isfinite(v) && (s += v * inst.patch.wtgcell[p])
        end
        s
    end
    sminn_tot() = (s = 0.0; @inbounds for j in 1:nl; s += sns().sminn_vr_col[c_soil, j] * dzd[j]; end; s)

    # --- RAW-POOL diagnostics --------------------------------------------------
    # The _col summary fields (totsomc_col/totvegc_col/totsomn_col) and gpp_patch
    # are NOT refreshed in this CN driver path: their summary routines
    # (soil_bgc_carbon_state_summary! / cnveg_carbon_flux_summary!) are unported
    # placeholders (clm_driver.jl ~952), so they keep their restart-init NaN/0.
    # Aggregate the actual prognostic pools instead, matching the summary formulas:
    #   TOTSOMC/N = sum_j sum_{is_soil pools} decomp_*pools_vr[c,j,l] * dz_decomp[j]
    #   TOTVEGC   = patch sum of all veg C pools (displayed + storage + xfer)
    #   GPP       = psnsun_to_cpool + psnshade_to_cpool (the real photosynthate flux)
    is_soil = (length(inst.decomp_cascade.is_soil) > 0) ? inst.decomp_cascade.is_soil : Bool[]
    function som_raw(vr)
        s = 0.0
        @inbounds for l in 1:size(vr, 3)
            (l <= length(is_soil) && is_soil[l]) || continue
            for j in 1:nl
                v = vr[c_soil, j, l]; isfinite(v) && (s += v * dzd[j])
            end
        end
        s
    end
    totsomc_raw() = som_raw(scs().decomp_cpools_vr_col)
    totsomn_raw() = som_raw(sns().decomp_npools_vr_col)
    vegc_stems = ("leafc","leafc_storage","leafc_xfer","frootc","frootc_storage",
        "frootc_xfer","livestemc","livestemc_storage","livestemc_xfer","deadstemc",
        "deadstemc_storage","deadstemc_xfer","livecrootc","livecrootc_storage",
        "livecrootc_xfer","deadcrootc","deadcrootc_storage","deadcrootc_xfer",
        "cpool","gresp_storage","gresp_xfer")
    vegc_arrs = [getfield(ccs(), Symbol(s * "_patch")) for s in vegc_stems
                 if hasproperty(ccs(), Symbol(s * "_patch"))]
    function totvegc_raw()
        s = 0.0
        @inbounds for p in 1:np
            inst.patch.active[p] || continue
            t = 0.0; for a in vegc_arrs; v = a[p]; isfinite(v) && (t += v); end
            s += t * inst.patch.wtgcell[p]
        end
        s
    end
    gpp_raw() = p2g(ccf().psnsun_to_cpool_patch) + p2g(ccf().psnshade_to_cpool_patch)
    function litcwd_raw()  # column-integrated litter + CWD C (the non-is_soil pools)
        vr = scs().decomp_cpools_vr_col
        s = 0.0
        @inbounds for l in 1:size(vr, 3)
            (l <= length(is_soil) && !is_soil[l]) || continue
            for j in 1:nl
                v = vr[c_soil, j, l]; isfinite(v) && (s += v * dzd[j])
            end
        end
        s
    end
    totecosysc_raw() = totsomc_raw() + totvegc_raw() + litcwd_raw()

    # Julia-side daily-mean accumulators: (h0 name, getter)
    fields = [
        ("TOTSOMC",    totsomc_raw),
        ("TOTVEGC",    totvegc_raw),
        ("TOTECOSYSC", totecosysc_raw),
        ("TOTSOMN",    totsomn_raw),
        ("SMINN",      sminn_tot),
        ("GPP",        () -> p2g(ccf().gpp_patch)),  # model's gpp_patch (cnveg flux summary now wired)
        ("GPP_psn",    gpp_raw),                      # cross-check: raw psn->cpool streams
        ("NPP",        () -> p2g(ccf().npp_patch)),   # npp_patch (summary now wired)
        ("TLAI",       () -> p2g(inst.canopystate.tlai_patch)),
        ("LEAFC",      () -> p2g(ccs().leafc_patch)),
        ("TG",         () -> inst.temperature.t_grnd_col[c_soil]),
        ("EFLX_LH_TOT",() -> p2g(inst.energyflux.eflx_lh_tot_patch)),
        ("FSH",        () -> p2g(inst.energyflux.eflx_sh_tot_patch)),
    ]
    daily = Dict(f[1] => Float64[] for f in fields)
    nf_total = 0

    @printf("CN annual: nc=%d np=%d soil_col=%d | injecting %s\n", nc, np, c_soil, basename(CN_RST))
    t0 = time()
    for d in 1:ndays
        acc = Dict(f[1] => 0.0 for f in fields); ns = 0
        for h in 1:24
            cur = start_date + Day(d - 1) + Hour(h - 1)
            calday = CLM.get_curr_calday(tm); (declin, _) = CLM.compute_orbital(calday)
            (declinm1, _) = CLM.compute_orbital(calday - 3600.0 / CLM.SECSPDAY)
            nextsw = calday + 3600.0 / CLM.SECSPDAY
            CLM.init_daylength!(inst.gridcell, declin, declinm1, CLM.ORB_OBLIQR_DEFAULT, 1:ng)
            CLM.advance_timestep!(tm)
            CLM.read_forcing_step!(fr, inst.atm2lnd, _force_date(cur), ng, nc)
            CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
            (yr, mon, dy, tod) = CLM.get_curr_date(tm)
            CLM.clm_drv!(config, inst, filt, filt_ia, bounds, true, nextsw, declin, declin,
                CLM.ORB_OBLIQR_DEFAULT, false, false, "", false;
                nstep = tm.nstep, is_first_step = false,
                is_beg_curr_day = CLM.is_beg_curr_day(tm), is_end_curr_day = CLM.is_end_curr_day(tm),
                is_beg_curr_year = CLM.is_beg_curr_year(tm), dtime = 3600.0, mon = mon, day = dy,
                photosyns = inst.photosyns)
            # cn_veg_struct_update! (leafc -> tlai/elai/htop) now runs INSIDE clm_drv!
            # (wired into the global driver), so no harness-side call is needed.
            if get(ENV, "CN_PROBE", "") != "" && d == 1 && h == 1
                ns2 = CLM.varpar.nlevsno
                ws = inst.water.waterstatebulk_inst.ws
                @printf("  [post-step1] t_grnd=%s t_veg(p2)=%s smp_l(c,1)=%s sabv(p2)=%s parsun(p2)=%s\n",
                    string(inst.temperature.t_grnd_col[c_soil]), string(inst.temperature.t_veg_patch[2]),
                    string(inst.soilstate.smp_l_col[c_soil, 1]), string(inst.solarabs.sabv_patch[2]),
                    string(inst.solarabs.parsun_z_patch[2, 1]))
                @printf("  [post-step1] h2osoi_liq soil[1:4]=%s albgrd(c,1)=%s coszen=%s\n",
                    string(round.(ws.h2osoi_liq_col[c_soil, (ns2+1):(ns2+4)], digits=2)),
                    string(inst.surfalb.albgrd_col[c_soil, 1]),
                    string(inst.surfalb.coszen_col[c_soil]))
                @printf("  [post-step1] vegwp(p2,1:4)=%s bsun=%s bsha=%s qtran=%s qevap=%s\n",
                    string(inst.canopystate.vegwp_patch[2, 1:4]),
                    string(inst.energyflux.bsun_patch[2]), string(inst.energyflux.bsha_patch[2]),
                    string(inst.water.waterfluxbulk_inst.wf.qflx_tran_veg_col[c_soil]),
                    string(inst.water.waterfluxbulk_inst.wf.qflx_evap_veg_patch[2]))
                @printf("  [post-step1] liqvol soil[1:3]=%s k_soil_root(p2)[1:3]=%s hk_l(c)[1:3]=%s\n",
                    string(round.(inst.water.waterdiagnosticbulk_inst.h2osoi_liqvol_col[c_soil, (ns2+1):(ns2+3)], digits=3)),
                    string(inst.soilstate.k_soil_root_patch[2, 1:3]),
                    string(inst.soilstate.hk_l_col[c_soil, 1:3]))
            end
            for (nm, g) in fields; acc[nm] += g(); end
            ns += 1
        end
        for (nm, _) in fields; push!(daily[nm], acc[nm] / ns); end
        if !isfinite(inst.temperature.t_grnd_col[c_soil]); nf_total += 1; end
        d % 30 == 0 && @printf("  day %3d  TOTSOMC=%.0f TOTVEGC=%.1f GPP=%.3e TLAI=%.3f Tg=%.1f\n",
                d, daily["TOTSOMC"][end], daily["TOTVEGC"][end], daily["GPP"][end],
                daily["TLAI"][end], daily["TG"][end])
        # phenology-onset evolution probe for the season_decid grass patch (p3)
        if get(ENV, "CN_PHENO", "") != "" && d % 15 == 0 &&
           hasproperty(inst.bgc_vegetation, :cnveg_state_inst)
            cvs = inst.bgc_vegetation.cnveg_state_inst
            pd = findfirst(p -> inst.patch.active[p] &&
                (inst.patch.itype[p]+1 <= length(CLM.pftcon.season_decid)) &&
                CLM.pftcon.season_decid[inst.patch.itype[p]+1] == 1.0, 1:np)
            if pd !== nothing
                tg = inst.gridcell.dayl
                @printf("    [pheno d%3d p%d] dormant=%g gddflag=%g gdd=%.1f days_act=%.1f swi=%.1f onflag=%g oncnt=%.0f offflag=%g leafc=%.2f tlai=%.3f dayl=%.0f\n",
                    d, pd, cvs.dormant_flag_patch[pd], cvs.onset_gddflag_patch[pd],
                    cvs.onset_gdd_patch[pd], cvs.days_active_patch[pd], cvs.onset_swi_patch[pd],
                    cvs.onset_flag_patch[pd], cvs.onset_counter_patch[pd],
                    cvs.offset_flag_patch[pd], ccs().leafc_patch[pd],
                    inst.canopystate.tlai_patch[pd],
                    (length(tg) > 0 ? tg[1] : NaN))
            end
        end
    end
    CLM.forcing_reader_close!(fr)
    @printf("  ran %d days in %.1fs\n\n", ndays, time() - t0)

    # ---- monthly comparison vs Fortran h0 ----
    fds = NCDataset(CN_H0, "r")
    mstart = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]
    mend   = [31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365]
    mons = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
    println("="^78)
    println("  CN ANNUAL MONTHLY MEANS — Julia vs Fortran (model 2202 / forcing 2003)")
    println("="^78)
    for (nm, _) in fields
        haskey(fds, nm) || continue
        jd = daily[nm]; nrec = min(length(jd), ndays)
        @printf("\n%-12s   %12s %12s %10s\n", nm, "Julia", "Fortran", "rel")
        for m in 1:12
            d1, d2 = mstart[m], min(mend[m], nrec)
            d1 > nrec && continue
            jm = mean(@view jd[d1:d2])
            fvals = [_fv(fds[nm][1, dd]) for dd in d1:d2]; fvals = filter(isfinite, fvals)
            fm = isempty(fvals) ? NaN : mean(fvals)
            rel = (isnan(fm) || isnan(jm)) ? NaN : abs(jm - fm) / (1 + abs(fm))
            @printf("  %-3s        %12.4g %12.4g %10.2e\n", mons[m], jm, fm, rel)
        end
    end
    close(fds)
    @printf("\nfiniteness: %s\n", nf_total == 0 ? "PASS (Tg finite all days)" : "FAIL ($nf_total bad days)")
    return inst, bounds
end

run_cn_annual(; ndays = length(ARGS) >= 1 && occursin(r"^\d+$", ARGS[1]) ? parse(Int, ARGS[1]) : 365)
