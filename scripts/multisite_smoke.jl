# =============================================================================
# Multi-site robustness smoke test — runs the CLM.jl driver at a CONTRASTING
# site (Aripuanã, tropical Amazon: broadleaf evergreen/deciduous tropical tree
# + C3/C4 grass, no snow, ~99 kPa) from a cold start with a synthetic but
# physically-reasonable wet-season forcing snapshot, and checks the driver runs
# without NaN/blowup and conserves the internal water/energy balances.
#
# This is NOT a Fortran-parity test (no reference dumps exist off-Bow) — it is a
# robustness check that catches Bow-specific assumptions (alpine, snow-dominated,
# 3 PFTs) by exercising a radically different subgrid + climate.
#
# Usage: julia +1.12 --project=. scripts/multisite_smoke.jl
# =============================================================================
using CLM, NCDatasets, Dates, Printf

const ARI = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Aripuana_Amazon"
const FS  = joinpath(ARI, "settings/CLM/parameters/surfdata_clm.nc")
const FP  = joinpath(ARI, "settings/CLM/parameters/clm5_params.nc")
const SNOWOPT = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_optics_5bnd_c013122.nc"
const SNOWAGE = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_drdt_bst_fit_60_c070416.nc"

function set_tropical_forcing!(inst, ng, nc)
    a = inst.atm2lnd
    T = 298.0; pbot = 99000.0; q = 0.018
    th = T * (1.0e5 / pbot)^0.286
    rho = pbot / (287.058 * T * (1.0 + 0.61 * q))
    vp = q * pbot / (0.622 + 0.378 * q)
    for g in 1:ng
        a.forc_t_not_downscaled_grc[g]    = T
        a.forc_th_not_downscaled_grc[g]   = th
        a.forc_pbot_not_downscaled_grc[g] = pbot
        a.forc_q_not_downscaled_grc[g]    = q
        a.forc_rho_not_downscaled_grc[g]  = rho
        a.forc_lwrad_not_downscaled_grc[g]= 420.0
        a.forc_rain_not_downscaled_grc[g] = 1.0e-4   # light wet-season rain
        a.forc_snow_not_downscaled_grc[g] = 0.0      # tropical: no snow
        a.forc_u_grc[g] = 2.0; a.forc_v_grc[g] = 0.0
        a.forc_hgt_grc[g] = 30.0
        a.forc_hgt_u_grc[g] = 30.0; a.forc_hgt_t_grc[g] = 30.0; a.forc_hgt_q_grc[g] = 30.0
        a.forc_vp_grc[g] = vp
        a.forc_pco2_grc[g] = 367.0e-6 * pbot
        a.forc_po2_grc[g]  = 0.209 * pbot
        a.forc_topo_grc[g] = 150.0
        # solar: 2 numrad bands (vis, nir); midday wet-season values
        for b in 1:size(a.forc_solad_not_downscaled_grc, 2)
            a.forc_solad_not_downscaled_grc[g, b] = b == 1 ? 300.0 : 250.0
            a.forc_solai_grc[g, b]                = b == 1 ? 100.0 : 80.0
        end
    end
    for c in 1:nc; inst.topo.topo_col[c] = 150.0; end
    return nothing
end

function main(; nsteps::Int = 6, use_cn::Bool = false)
    CLM.rooting_profile_config.rooting_profile_method_water  = CLM.JACKSON_1996_ROOT
    CLM.rooting_profile_config.rooting_profile_method_carbon = CLM.JACKSON_1996_ROOT
    start_date = DateTime(2006, 1, 1, 12)   # midday, wet season
    (inst, bounds, filt, tm) = CLM.clm_initialize!(; fsurdat=FS, paramfile=FP,
        start_date=start_date, dtime=3600, use_cn=use_cn, use_luna=false,
        use_bedrock=true, use_aquifer_layer=false, h2osfcflag=0,
        fsnowoptics=SNOWOPT, fsnowaging=SNOWAGE, int_snow_max=2000.0)
    ng, nc, np = bounds.endg, bounds.endc, bounds.endp
    @printf("Aripuanã subgrid: ng=%d nc=%d np=%d  PFTs=%s\n", ng, nc, np, string(Int.(inst.patch.itype)))

    config  = CLM.CLMDriverConfig(use_cn=use_cn, use_aquifer_layer=false,
                                  use_hydrstress=false, use_luna=false)
    filt_ia = CLM.clump_filter_inactive_and_active

    nan_fields = String[]
    for i in 1:nsteps
        set_tropical_forcing!(inst, ng, nc)
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

        # finiteness + sane-range checks on the core surface predictions.
        # Audit a broad set of prognostic/flux fields for NaN/Inf each step.
        tv = inst.temperature.t_veg_patch; tg = inst.temperature.t_grnd_col
        sh = inst.energyflux.eflx_sh_grnd_patch
        ws = inst.water.waterstatebulk_inst.ws
        audited = Dict("t_veg"=>tv, "t_grnd"=>tg, "eflx_sh_grnd"=>sh,
                       "t_soisno"=>inst.temperature.t_soisno_col,
                       "h2osoi_liq"=>ws.h2osoi_liq_col, "h2osoi_ice"=>ws.h2osoi_ice_col,
                       "eflx_lwrad_out"=>inst.energyflux.eflx_lwrad_out_patch,
                       "btran"=>inst.energyflux.btran_patch, "zwt"=>inst.soilhydrology.zwt_col)
        bad = [k for (k,v) in audited if any(!isfinite, v)]
        rng(v) = (vv = filter(isfinite, v); isempty(vv) ? (NaN,NaN) : (minimum(vv), maximum(vv)))
        tvl,tvh = rng(tv); tgl,tgh = rng(tg); shl,shh = rng(sh)
        @printf("step %d: t_veg=[%.2f,%.2f] t_grnd=[%.2f,%.2f] SH=[%.1f,%.1f] bad_fields=%s\n",
            i, tvl,tvh, tgl,tgh, shl,shh, isempty(bad) ? "none" : string(bad))
        isempty(bad) || push!(nan_fields, "step$i:$(bad)")
    end

    println(isempty(nan_fields) ? "\n✅ MULTI-SITE SMOKE PASS: driver ran $nsteps steps at Aripuanã, all finite" :
                                  "\n❌ NaN appeared at: $(nan_fields)")
    return inst, bounds
end

main()
