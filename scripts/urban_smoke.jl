# =============================================================================
# Urban robustness smoke test — runs the CLM.jl driver on a surfdata that
# contains urban (isturb) land units from a cold start with a warm daytime
# forcing snapshot, and checks the driver runs without NaN/blowup on the urban
# column.
#
# Surfdata: the CTSM single-point Mexico City test point
# (surfdata_1x1_mexicocityMEX_hist_16pfts_CMIP6_2000_c231103.nc), which has
# PCT_URBAN > 0 across the three density classes (TBD/HD/MD = isturb 7/8/9) plus
# urban morphology (CANYON_HWR, HT_ROOF, roof/wall/road thermal & radiative
# fields). This exercises the urban land-unit path end-to-end, which the
# soil-only / glacier robustness tests never touch.
#
# The urban morphology is read by read_urban_input! (src/infrastructure/
# urban_input.jl) and pushed into the landunit + urbanparams state by
# urbanparams_populate!, both wired into clm_initialize!. Without that wiring the
# isturb columns run on fill(NaN) morphology and blow up.
#
# This is NOT a Fortran-parity test (no urban reference dumps yet) — it is a
# finiteness / robustness gate.
#
# Usage: julia +1.12 --project=. scripts/urban_smoke.jl
# =============================================================================
# NB: `Base.include(@__MODULE__, ...)`, not bare `include`. Several of these
# scripts are loaded by their tests into a fresh `Module(:X)`, which does NOT
# bind a bare `include` — that form throws UndefVarError there.
Base.include(@__MODULE__, joinpath(@__DIR__, "..", "test", "testdata.jl"))

using CLM, NCDatasets, Dates, Printf

const URB_FS = symfluence_path("installs", "clm", "python", "ctsm", "test", "testinputs",
                               "surfdata_1x1_mexicocityMEX_hist_16pfts_CMIP6_2000_c231103.nc")
const BOW_CAL = domain_params_dir("domain_Bow_at_Banff_lumped")
const URB_FP = joinpath(BOW_CAL, "clm5_params.nc")
const SNOWOPT = snicar_optics()
const SNOWAGE = snicar_aging()

# Warm, dry daytime forcing (Mexico City ~2240 m): mild air, no snow.
function set_urban_forcing!(inst, ng, nc)
    a = inst.atm2lnd
    T = 290.0; pbot = 78000.0; q = 0.006
    th = T * (1.0e5 / pbot)^0.286
    rho = pbot / (287.058 * T * (1.0 + 0.61 * q))
    vp = q * pbot / (0.622 + 0.378 * q)
    for g in 1:ng
        a.forc_t_not_downscaled_grc[g]    = T
        a.forc_th_not_downscaled_grc[g]   = th
        a.forc_pbot_not_downscaled_grc[g] = pbot
        a.forc_q_not_downscaled_grc[g]    = q
        a.forc_rho_not_downscaled_grc[g]  = rho
        a.forc_lwrad_not_downscaled_grc[g]= 320.0
        a.forc_rain_not_downscaled_grc[g] = 0.0
        a.forc_snow_not_downscaled_grc[g] = 0.0
        a.forc_u_grc[g] = 3.0; a.forc_v_grc[g] = 0.0
        # forc_wind_grc must be set (new-snow density / canyon wind need it; an
        # unset NaN wind -> NaN downstream). Real forcing always provides it.
        isempty(a.forc_wind_grc) || (a.forc_wind_grc[g] = sqrt(3.0^2 + 0.0^2))
        a.forc_hgt_grc[g] = 30.0
        a.forc_hgt_u_grc[g] = 30.0; a.forc_hgt_t_grc[g] = 30.0; a.forc_hgt_q_grc[g] = 30.0
        a.forc_vp_grc[g] = vp
        a.forc_pco2_grc[g] = 367.0e-6 * pbot
        a.forc_po2_grc[g]  = 0.209 * pbot
        a.forc_topo_grc[g] = 2240.0
        for b in 1:size(a.forc_solad_not_downscaled_grc, 2)
            a.forc_solad_not_downscaled_grc[g, b] = b == 1 ? 400.0 : 320.0
            a.forc_solai_grc[g, b]                = b == 1 ? 120.0 : 90.0
        end
    end
    for c in 1:nc; inst.topo.topo_col[c] = 2240.0; end
    return nothing
end

# true iff landunit type is one of the urban density classes (isturb 7/8/9).
_is_urban_lun(it) = it >= CLM.ISTURB_MIN && it <= CLM.ISTURB_MIN + CLM.NUMURBL - 1

function main(; nsteps::Int = 6)
    CLM.rooting_profile_config.rooting_profile_method_water  = CLM.JACKSON_1996_ROOT
    CLM.rooting_profile_config.rooting_profile_method_carbon = CLM.JACKSON_1996_ROOT
    # ~local solar noon at Mexico City's lon (~-99°): 12:00 local ≈ 18:00 UTC.
    start_date = DateTime(2000, 6, 15, 18)
    (inst, bounds, filt, tm) = CLM.clm_initialize!(; fsurdat=URB_FS, paramfile=URB_FP,
        start_date=start_date, dtime=3600, use_cn=false, use_luna=false,
        use_bedrock=true, use_aquifer_layer=false, h2osfcflag=0,
        fsnowoptics=SNOWOPT, fsnowaging=SNOWAGE, int_snow_max=2000.0)
    # Building-temperature method. Default = SIMPLE (CLM4.5, building_hac!).
    # The prognostic (CLM5.0 / Oleson 2015) building energy balance is now also
    # wired (building_temperature! in soil_temperature.jl) — set URB_PROG=1 to
    # exercise it end-to-end through the driver.
    CLM.urban_ctrl.building_temp_method = get(ENV, "URB_PROG", "") != "" ?
        CLM.BUILDING_TEMP_METHOD_PROG : CLM.BUILDING_TEMP_METHOD_SIMPLE

    ng, nc, np = bounds.endg, bounds.endc, bounds.endp
    @printf("Urban subgrid: ng=%d nc=%d np=%d  landunit_types=%s\n",
        ng, nc, np, string(Int.(inst.landunit.itype)))

    # locate the urban (isturb) columns + report their morphology was populated
    urb_cols = [c for c in 1:nc if _is_urban_lun(inst.landunit.itype[inst.column.landunit[c]])]
    @printf("Urban columns: %s\n", string(urb_cols))
    if !isempty(urb_cols) && get(ENV, "URB_DEBUG", "") != ""
        for c in urb_cols
            l = inst.column.landunit[c]
            @printf("  [urb-init c=%d l=%d itype=%d] canyon_hwr=%s ht_roof=%s wtroad_perv=%s\n",
                c, l, inst.landunit.itype[l],
                string(round(inst.landunit.canyon_hwr[l], digits=3)),
                string(round(inst.landunit.ht_roof[l], digits=2)),
                string(round(inst.landunit.wtroad_perv[l], digits=3)))
        end
    end

    config  = CLM.CLMDriverConfig(use_cn=false, use_aquifer_layer=false,
                                  use_hydrstress=false, use_luna=false)
    filt_ia = CLM.clump_filter_inactive_and_active

    nan_fields = String[]
    for i in 1:nsteps
        set_urban_forcing!(inst, ng, nc)
        calday = CLM.get_curr_calday(tm); (declin, _) = CLM.compute_orbital(calday)
        (declinm1, _) = CLM.compute_orbital(calday - 3600.0 / CLM.SECSPDAY)
        nextsw = calday + 3600.0 / CLM.SECSPDAY
        CLM.init_daylength!(inst.gridcell, declin, declinm1, CLM.ORB_OBLIQR_DEFAULT, 1:ng)
        CLM.advance_timestep!(tm)
        CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
        (yr, mon, d, tod) = CLM.get_curr_date(tm)
        CLM.clm_drv!(config, inst, filt, filt_ia, bounds, true, nextsw, declin, declin,
            CLM.ORB_OBLIQR_DEFAULT, false, false, "", false;
            nstep=tm.nstep, is_first_step=(i==1),
            is_beg_curr_day=CLM.is_beg_curr_day(tm), is_end_curr_day=CLM.is_end_curr_day(tm),
            is_beg_curr_year=CLM.is_beg_curr_year(tm), dtime=3600.0, mon=mon, day=d,
            photosyns=inst.photosyns)

        tg = inst.temperature.t_grnd_col
        ws = inst.water.waterstatebulk_inst.ws
        # audit the urban columns specifically (so a finite non-urban column can't
        # mask an urban NaN), plus the usual full-domain state arrays.
        audited = Dict("t_grnd"=>tg, "t_soisno"=>inst.temperature.t_soisno_col,
                       "h2osoi_liq"=>ws.h2osoi_liq_col, "h2osoi_ice"=>ws.h2osoi_ice_col,
                       "h2osno"=>ws.h2osno_no_layers_col,
                       "eflx_lwrad_out"=>inst.energyflux.eflx_lwrad_out_patch,
                       "eflx_sh_grnd"=>inst.energyflux.eflx_sh_grnd_patch)
        bad = [k for (k,v) in audited if any(!isfinite, v)]
        # explicit urban-column finiteness check on t_grnd / t_soisno
        urb_bad = [c for c in urb_cols if !isfinite(tg[c]) ||
                   any(!isfinite, @view inst.temperature.t_soisno_col[c, :])]
        rng(v) = (vv = filter(isfinite, v); isempty(vv) ? (NaN,NaN) : (minimum(vv), maximum(vv)))
        tgu = isempty(urb_cols) ? (NaN,NaN) : rng(@view tg[urb_cols])
        @printf("step %d: t_grnd_urban=[%.2f,%.2f] bad_fields=%s urb_bad_cols=%s\n",
            i, tgu[1], tgu[2], isempty(bad) ? "none" : string(bad),
            isempty(urb_bad) ? "none" : string(urb_bad))
        isempty(bad) || push!(nan_fields, "step$i:$(bad)")
        isempty(urb_bad) || push!(nan_fields, "step$i:urb_cols$(urb_bad)")
    end

    ok = isempty(nan_fields) && !isempty(urb_cols)
    if isempty(urb_cols)
        println("\n❌ URBAN SMOKE FAIL: no urban (isturb) columns in the subgrid")
    elseif isempty(nan_fields)
        println("\n✅ URBAN SMOKE PASS: driver ran $nsteps steps with finite urban columns $(urb_cols)")
    else
        println("\n❌ NaN appeared at: $(nan_fields)")
    end
    return (inst=inst, bounds=bounds, nbad=length(nan_fields), nurb=length(urb_cols))
end

# Returns true iff the urban domain ran finite for `nsteps` steps; missing if the
# urban surfdata / param file is absent (gated regression test skips).
function urban_smoke_ok(; nsteps::Int = 6)
    (isfile(URB_FS) && isfile(URB_FP)) ||
        (@info "urban smoke: inputs absent, skipping" URB_FS URB_FP; return missing)
    r = main(; nsteps=nsteps)
    return r.nbad == 0 && r.nurb > 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
