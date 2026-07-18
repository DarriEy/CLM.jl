# =============================================================================
# Snow / cold-soil robustness smoke test — runs the CLM.jl driver on the Bow
# natural-veg column (needleleaf-evergreen tree + grass = boreal-like) from a cold
# start with a cold, heavy-snow forcing, for enough steps to build a MULTI-LAYER
# snowpack. Exercises the snow-layer path (build/combine/divide snow layers, snow
# compaction, snow-layer phase change) and the frozen-soil path that the boreal and
# Arctic eval domains will hit — which the 6-step glacier_smoke (snl=0, <2mm snow)
# never reaches. Asserts no field goes NaN/Inf and that snow layers actually form.
#
# Usage: julia +1.12 --project=. scripts/snow_smoke.jl
# =============================================================================
include(joinpath(@__DIR__, "..", "test", "testdata.jl"))

using CLM, NCDatasets, Dates, Printf

const BOW_CAL = domain_params_dir("domain_Bow_at_Banff_lumped")
const SNOW_FS = joinpath(BOW_CAL, "surfdata_clm.nc")
const SNOW_FP = joinpath(BOW_CAL, "clm5_params.nc")
const SNOWOPT = snicar_optics()
const SNOWAGE = snicar_aging()

function set_cold_snow_forcing!(inst, ng, nc)
    a = inst.atm2lnd
    T = 258.0; pbot = 80000.0; q = 0.0008
    th = T * (1.0e5 / pbot)^0.286
    rho = pbot / (287.058 * T * (1.0 + 0.61 * q))
    vp = q * pbot / (0.622 + 0.378 * q)
    for g in 1:ng
        a.forc_t_not_downscaled_grc[g]    = T
        a.forc_th_not_downscaled_grc[g]   = th
        a.forc_pbot_not_downscaled_grc[g] = pbot
        a.forc_q_not_downscaled_grc[g]    = q
        a.forc_rho_not_downscaled_grc[g]  = rho
        a.forc_lwrad_not_downscaled_grc[g]= 210.0
        a.forc_rain_not_downscaled_grc[g] = 0.0
        a.forc_snow_not_downscaled_grc[g] = 3.0e-4   # heavy snowfall -> deep snowpack
        a.forc_u_grc[g] = 3.0; a.forc_v_grc[g] = 0.0
        isempty(a.forc_wind_grc) || (a.forc_wind_grc[g] = 3.0)
        a.forc_hgt_grc[g] = 30.0
        a.forc_hgt_u_grc[g] = 30.0; a.forc_hgt_t_grc[g] = 30.0; a.forc_hgt_q_grc[g] = 30.0
        a.forc_vp_grc[g] = vp
        a.forc_pco2_grc[g] = 367.0e-6 * pbot
        a.forc_po2_grc[g]  = 0.209 * pbot
        a.forc_topo_grc[g] = 1900.0
        for b in 1:size(a.forc_solad_not_downscaled_grc, 2)
            a.forc_solad_not_downscaled_grc[g, b] = b == 1 ? 120.0 : 100.0
            a.forc_solai_grc[g, b]                = b == 1 ? 50.0 : 40.0
        end
    end
    for c in 1:nc; inst.topo.topo_col[c] = 1900.0; end
    return nothing
end

function main(; nsteps::Int = 240)
    CLM.rooting_profile_config.rooting_profile_method_water  = CLM.JACKSON_1996_ROOT
    CLM.rooting_profile_config.rooting_profile_method_carbon = CLM.JACKSON_1996_ROOT
    start_date = DateTime(2006, 2, 15, 20)   # ~local noon at Bow's lon (winter)
    (inst, bounds, filt, tm) = CLM.clm_initialize!(; fsurdat=SNOW_FS, paramfile=SNOW_FP,
        start_date=start_date, dtime=3600, use_cn=false, use_luna=false,
        use_bedrock=true, use_aquifer_layer=false, h2osfcflag=0,
        fsnowoptics=SNOWOPT, fsnowaging=SNOWAGE, int_snow_max=2000.0)
    ng, nc, np = bounds.endg, bounds.endc, bounds.endp
    @printf("Snow subgrid: ng=%d nc=%d np=%d  landunit_types=%s\n",
        ng, nc, np, string(Int.(inst.landunit.itype)))
    config  = CLM.CLMDriverConfig(use_cn=false, use_aquifer_layer=false,
                                  use_hydrstress=false, use_luna=false)
    filt_ia = CLM.clump_filter_inactive_and_active
    c1 = findfirst(filt.nolakec)

    nan_fields = String[]; max_snl = 0
    for i in 1:nsteps
        set_cold_snow_forcing!(inst, ng, nc)
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
        audited = Dict("t_grnd"=>tg, "t_soisno"=>inst.temperature.t_soisno_col,
                       "h2osoi_liq"=>ws.h2osoi_liq_col, "h2osoi_ice"=>ws.h2osoi_ice_col,
                       "h2osno"=>ws.h2osno_no_layers_col,
                       "snow_depth"=>inst.water.waterdiagnosticbulk_inst.snow_depth_col,
                       "eflx_lwrad_out"=>inst.energyflux.eflx_lwrad_out_patch,
                       "eflx_sh_grnd"=>inst.energyflux.eflx_sh_grnd_patch)
        bad = [k for (k,v) in audited if any(!isfinite, v)]
        max_snl = max(max_snl, -minimum(inst.column.snl))   # snl<0 => snow layers
        isempty(bad) || (push!(nan_fields, "step$i:$(bad)"); break)
        i % 60 == 0 && @printf("step %3d: t_grnd=%.2f h2osno=%.1f snow_depth=%.3f snl=%d\n",
            i, tg[c1], ws.h2osno_no_layers_col[c1] + sum(@view ws.h2osoi_ice_col[c1, 1:CLM.varpar.nlevsno]),
            inst.water.waterdiagnosticbulk_inst.snow_depth_col[c1], inst.column.snl[c1])
    end
    if isempty(nan_fields)
        @printf("\n✅ SNOW SMOKE PASS: %d steps, all finite; max snow layers formed=%d\n", nsteps, max_snl)
    else
        @printf("\n❌ NaN appeared at: %s\n", string(nan_fields))
    end
    return (inst=inst, bounds=bounds, nbad=length(nan_fields), max_snl=max_snl)
end

# Gated regression entry: finite for nsteps AND snow layers actually formed.
function snow_smoke_ok(; nsteps::Int = 240)
    (isfile(SNOW_FS) && isfile(SNOW_FP)) ||
        (@info "snow smoke: inputs absent, skipping" SNOW_FS; return missing)
    r = main(; nsteps=nsteps)
    return r.nbad == 0 && r.max_snl >= 1
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
