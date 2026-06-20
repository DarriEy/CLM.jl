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
# NEXT: chase the use_cn=true step-1 soil-hydrology NaN (likely the CN soil-water
# plant sink / smp_l coupling) — then this harness gives the first full CN parity.
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
    end

    config  = CLM.CLMDriverConfig(use_cn = true, use_aquifer_layer = false,
                                  use_hydrstress = true, use_luna = true)
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

    # Julia-side daily-mean accumulators: (h0 name, getter)
    fields = [
        ("TOTSOMC",    () -> scs().totsomc_col[c_soil]),
        ("TOTVEGC",    () -> ccs().totvegc_col[c_soil]),
        ("TOTECOSYSC", () -> scs().totecosysc_col[c_soil]),
        ("TOTSOMN",    () -> sns().totsomn_col[c_soil]),
        ("SMINN",      () -> sminn_tot()),
        ("GPP",        () -> p2g(ccf().gpp_patch)),
        ("NPP",        () -> p2g(ccf().npp_patch)),
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
            end
            for (nm, g) in fields; acc[nm] += g(); end
            ns += 1
        end
        for (nm, _) in fields; push!(daily[nm], acc[nm] / ns); end
        if !isfinite(inst.temperature.t_grnd_col[c_soil]); nf_total += 1; end
        d % 30 == 0 && @printf("  day %3d  TOTSOMC=%.0f TOTVEGC=%.1f GPP=%.3e TLAI=%.3f Tg=%.1f\n",
                d, daily["TOTSOMC"][end], daily["TOTVEGC"][end], daily["GPP"][end],
                daily["TLAI"][end], daily["TG"][end])
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
