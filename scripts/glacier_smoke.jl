# =============================================================================
# Glacier robustness smoke test — runs the CLM.jl driver on a pure-GLACIER
# (istice) column from a cold start with a cold, snowy alpine forcing snapshot,
# and checks the driver runs without NaN/blowup.
#
# The glacier surfdata (test_inputs/glacier/surfdata_glacier100.nc) is the Bow
# site with PCT_GLACIER=100 (all other landunit fractions 0) — see
# gen_glacier_surfdata.jl — so this isolates the istice land-unit path
# (generic soil/snow column with glacier albedo/runoff behavior) which the
# soil-only multisite robustness test never exercises.
#
# This is NOT a Fortran-parity test (no glacier reference dumps yet) — it is a
# finiteness/robustness gate: glacier columns hit driver branches (glacier
# filters, ice-melt handling, snow-on-ice) that have never been run end-to-end.
#
# Usage: julia +1.12 --project=. scripts/glacier_smoke.jl
#
# STATUS (2026-06-20): the istice cold-start path is a MULTI-LAYER NaN chase. Fixed so
# far: (1) cold-start albedo seed (cold_start.jl init_surface_albedo_cold!; was albedo=0
# -> over-absorption); (2) the harness clock is UTC -> use ~local noon (20:00 UTC at Bow's
# lon) + set forc_wind_grc (handle_new_snow's new-snow density needs it); (3) cv bedrock
# floor for ISTICE (soil_temperature.jl: deep bedrock layers had ice=liq=0 -> cv=0 ->
# fact=dtime/cv=Inf -> the band solve NaNs the whole glacier column). With (3) the glacier
# t_grnd is FINITE at step 1 (263.8 K). REMAINING (open): h2osno still NaN, traced into
# handle_new_snow! (snow accumulation) with all finite inputs (h2osno=0/int_snow=0/
# n_melt=0.4/wind=4/bifall finite) -> the NaN is in the snow accumulation / diagnostics
# (qflx_snow_grnd or the scf int_snow update) — a further cold-start snow layer. The
# REAL Iceland domain (Symfluence, building) with realistic forcing+spinup is the truer
# test of whether these synthetic-cold-start pathologies even manifest.
# =============================================================================
using CLM, NCDatasets, Dates, Printf

const GLAC_FS = joinpath(@__DIR__, "..", "test_inputs", "glacier", "surfdata_glacier100.nc")
const BOW_CAL = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/" *
                "domain_Bow_at_Banff_lumped/settings/CLM/parameters"
const GLAC_FP = joinpath(BOW_CAL, "clm5_params.nc")
const SNOWOPT = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_optics_5bnd_c013122.nc"
const SNOWAGE = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_drdt_bst_fit_60_c070416.nc"

# Cold, snowy alpine forcing (Bow ~1900 m): sub-freezing air, light snowfall.
function set_glacier_forcing!(inst, ng, nc)
    a = inst.atm2lnd
    T = 263.0; pbot = 80000.0; q = 0.001
    th = T * (1.0e5 / pbot)^0.286
    rho = pbot / (287.058 * T * (1.0 + 0.61 * q))
    vp = q * pbot / (0.622 + 0.378 * q)
    for g in 1:ng
        a.forc_t_not_downscaled_grc[g]    = T
        a.forc_th_not_downscaled_grc[g]   = th
        a.forc_pbot_not_downscaled_grc[g] = pbot
        a.forc_q_not_downscaled_grc[g]    = q
        a.forc_rho_not_downscaled_grc[g]  = rho
        a.forc_lwrad_not_downscaled_grc[g]= 220.0
        a.forc_rain_not_downscaled_grc[g] = 0.0
        a.forc_snow_not_downscaled_grc[g] = 5.0e-5   # light snowfall
        a.forc_u_grc[g] = 4.0; a.forc_v_grc[g] = 0.0
        # forc_wind_grc must be set (handle_new_snow!'s new-snow density needs it; an
        # unset NaN wind -> NaN snow density -> NaN h2osno). Real forcing provides it.
        isempty(a.forc_wind_grc) || (a.forc_wind_grc[g] = sqrt(4.0^2 + 0.0^2))
        a.forc_hgt_grc[g] = 30.0
        a.forc_hgt_u_grc[g] = 30.0; a.forc_hgt_t_grc[g] = 30.0; a.forc_hgt_q_grc[g] = 30.0
        a.forc_vp_grc[g] = vp
        a.forc_pco2_grc[g] = 367.0e-6 * pbot
        a.forc_po2_grc[g]  = 0.209 * pbot
        a.forc_topo_grc[g] = 1900.0
        for b in 1:size(a.forc_solad_not_downscaled_grc, 2)
            a.forc_solad_not_downscaled_grc[g, b] = b == 1 ? 250.0 : 200.0
            a.forc_solai_grc[g, b]                = b == 1 ? 90.0 : 70.0
        end
    end
    for c in 1:nc; inst.topo.topo_col[c] = 1900.0; end
    return nothing
end

function main(; nsteps::Int = 6)
    CLM.rooting_profile_config.rooting_profile_method_water  = CLM.JACKSON_1996_ROOT
    CLM.rooting_profile_config.rooting_profile_method_carbon = CLM.JACKSON_1996_ROOT
    # ~local solar noon at Bow's lon (~-115°): 12:00 local ≈ 20:00 UTC. The clock is
    # UTC, so 12:00 UTC would be ~04:00 local (night, coszen=0) — with a fixed solar
    # forcing that mismatch zeros the albedo and over-absorbs. June for a high sun angle.
    start_date = DateTime(2006, 6, 15, 20)
    (inst, bounds, filt, tm) = CLM.clm_initialize!(; fsurdat=GLAC_FS, paramfile=GLAC_FP,
        start_date=start_date, dtime=3600, use_cn=false, use_luna=false,
        use_bedrock=true, use_aquifer_layer=false, h2osfcflag=0,
        fsnowoptics=SNOWOPT, fsnowaging=SNOWAGE, int_snow_max=2000.0)
    ng, nc, np = bounds.endg, bounds.endc, bounds.endp
    @printf("Glacier subgrid: ng=%d nc=%d np=%d  landunit_types=%s\n",
        ng, nc, np, string(Int.(inst.landunit.itype)))

    # locate the glacier (istice) column + a patch on it
    gc = findfirst(c -> inst.landunit.itype[inst.column.landunit[c]] == 4, 1:nc)
    gp = gc === nothing ? 0 : findfirst(p -> inst.patch.column[p] == gc, 1:np)
    if gc !== nothing && get(ENV, "GLAC_DEBUG", "") != ""
        ss = inst.soilstate; ws = inst.water.waterstatebulk_inst.ws; ns = CLM.varpar.nlevsno
        @printf("  [glac-init c=%d p=%d] watsat[1:3]=%s sucsat[1:3]=%s bsw[1:3]=%s\n", gc, gp,
            string(round.(ss.watsat_col[gc,1:3],digits=3)), string(round.(ss.sucsat_col[gc,1:3],digits=1)),
            string(round.(ss.bsw_col[gc,1:3],digits=2)))
        @printf("  [glac-init] h2osno=%s snow_depth=%s snl=%s frac_sno=%s\n",
            string(ws.h2osno_no_layers_col[gc]),
            string(inst.water.waterdiagnosticbulk_inst.snow_depth_col[gc]),
            string(inst.column.snl[gc]),
            string(inst.water.waterdiagnosticbulk_inst.frac_sno_eff_col[gc]))
        @printf("  [glac-init] h2osoi_vol[1:3]=%s liq[1:3]=%s ice[1:3]=%s dz[1:3]=%s snl=%s\n",
            string(round.(ws.h2osoi_vol_col[gc,1:3],digits=3)),
            string(round.(ws.h2osoi_liq_col[gc,ns+1:ns+3],digits=2)),
            string(round.(ws.h2osoi_ice_col[gc,ns+1:ns+3],digits=1)),
            string(round.(inst.column.dz[gc,ns+1:ns+3],digits=3)), string(inst.column.snl[gc]))
    end

    config  = CLM.CLMDriverConfig(use_cn=false, use_aquifer_layer=false,
                                  use_hydrstress=false, use_luna=false)
    filt_ia = CLM.clump_filter_inactive_and_active

    nan_fields = String[]
    for i in 1:nsteps
        set_glacier_forcing!(inst, ng, nc)
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
                       "h2osno"=>ws.h2osno_no_layers_col, "snow_depth"=>inst.water.waterdiagnosticbulk_inst.snow_depth_col,
                       "eflx_lwrad_out"=>inst.energyflux.eflx_lwrad_out_patch,
                       "eflx_sh_grnd"=>inst.energyflux.eflx_sh_grnd_patch,
                       "qflx_snomelt"=>inst.water.waterfluxbulk_inst.wf.qflx_snomelt_col)
        bad = [k for (k,v) in audited if any(!isfinite, v)]
        rng(v) = (vv = filter(isfinite, v); isempty(vv) ? (NaN,NaN) : (minimum(vv), maximum(vv)))
        tgl,tgh = rng(tg); snl,snh = rng(ws.h2osno_no_layers_col)
        @printf("step %d: t_grnd=[%.2f,%.2f] h2osno=[%.2f,%.2f] bad_fields=%s\n",
            i, tgl,tgh, snl,snh, isempty(bad) ? "none" : string(bad))
        isempty(bad) || push!(nan_fields, "step$i:$(bad)")
        if i == 1
            ltypes = [inst.landunit.itype[inst.column.landunit[c]] for c in 1:nc]
            @printf("  [per-column] lun_type=%s  t_grnd_col=%s\n",
                string(Int.(ltypes)), string(tg[1:nc]))
            if gc !== nothing && get(ENV, "GLAC_DEBUG", "") != ""
                ns = CLM.varpar.nlevsno
                sa = inst.surfalb
                @printf("  [glac-step1 c=%d p=%d] t_grnd=%s tsoisno[1:3]=%s coszen=%s h2osno=%s\n",
                    gc, gp, string(tg[gc]),
                    string(inst.temperature.t_soisno_col[gc, ns+1:ns+3]),
                    string(round(sa.coszen_col[gc],digits=3)),
                    string(round(inst.water.waterstatebulk_inst.ws.h2osno_no_layers_col[gc],digits=3)))
                @printf("  [glac-step1] albgrd[1:2]=%s albgri[1:2]=%s albsod[1:2]=%s albsoi[1:2]=%s sabg=%s\n",
                    string(sa.albgrd_col[gc,1:2]), string(sa.albgri_col[gc,1:2]),
                    string(sa.albsod_col[gc,1:2]), string(sa.albsoi_col[gc,1:2]),
                    string(gp>0 ? inst.solarabs.sabg_patch[gp] : NaN))
            end
        end
    end

    println(isempty(nan_fields) ? "\n✅ GLACIER SMOKE PASS: driver ran $nsteps steps on a pure-glacier column, all finite" :
                                  "\n❌ NaN appeared at: $(nan_fields)")
    return (inst=inst, bounds=bounds, nbad=length(nan_fields))
end

# Returns true iff the glacier domain ran finite for `nsteps` steps; missing if
# the glacier surfdata / Bow param file is absent (gated regression test skips).
function glacier_smoke_ok(; nsteps::Int = 6)
    (isfile(GLAC_FS) && isfile(GLAC_FP)) ||
        (@info "glacier smoke: inputs absent, skipping" GLAC_FS GLAC_FP; return missing)
    r = main(; nsteps=nsteps)
    return r.nbad == 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
