# =============================================================================
# Glacier robustness smoke test — runs the CLM.jl driver on a pure-GLACIER
# (istice) column from a cold start with a cold, SNOWY alpine forcing, and
# checks the driver runs without NaN/blowup THROUGH SNOW-LAYER FORMATION.
#
# The glacier surfdata (test_inputs/glacier/surfdata_glacier100.nc) is the Bow
# site with PCT_GLACIER=100 (all other landunit fractions 0) — see
# gen_glacier_surfdata.jl — so this isolates the istice land-unit path
# (generic soil/snow column with glacier albedo/runoff behavior) which the
# soil-only multisite robustness test never exercises.
#
# This is NOT a Fortran-parity test (no glacier reference dumps yet) — it is a
# finiteness/robustness gate: glacier columns hit driver branches (glacier
# filters, ice-melt handling, snow-on-ice, SNICAR on a resolved snow layer)
# that have never been run end-to-end elsewhere.
#
# Usage: julia +1.12 --project=. scripts/glacier_smoke.jl
#
# ---------------------------------------------------------------------------
# THE GATE MUST FORM SNOW LAYERS. The original version of this harness used a
# LIGHT snowfall (5e-5 mm/s) over 6 steps, which never pushed snl below 0: the
# whole snow-layer code path (layer creation in handle_new_snow!, the snow band
# of the soil-temperature solve, SNICAR on a resolved snow layer, compaction /
# combine / divide) was structurally unreachable, so this test could not fail
# on it — and it didn't. A glacier column NaN'd on the step AFTER the first snow
# layer formed (the dead `aerosol_init_cold!` left the snow-layer aerosol masses
# at their allocation-time NaN → NaN mss_cnc → NaN SNICAR albedo/flx_abs → NaN
# sabg_lyr → NaN soil-temperature solve) and this harness stayed green.
#
# So the harness now:
#   (1) forces MODERATE snowfall (5e-4 mm/s) — a layer forms on step 1;
#   (2) runs long enough to cross layer creation, melt-back to the no-layer
#       state, and re-growth into a MULTI-layer snowpack;
#   (3) audits the per-ACTIVE-layer arrays (dz/z/t_soisno/h2osoi_*) as well as
#       the column scalars — the old audit never looked at dz/z at all;
#   (4) HARD-FAILS if snow layers never formed, so the gate cannot silently go
#       blind again.
# ---------------------------------------------------------------------------
#
# HISTORY (2026-06-20): the istice cold-start path first ran finite after a 4-layer fix
# (the soil-only multisite test never exercised istice):
#  (1) cold-start albedo seed (cold_start.jl init_surface_albedo_cold!; was albedo=0 ->
#      over-absorption sabg 488);
#  (2) harness: clock is UTC -> use ~local noon (20:00 UTC at Bow's lon) + set
#      forc_wind_grc (handle_new_snow's new-snow density needs it);
#  (3) cv bedrock floor for ISTICE (soil_temperature.jl: deep bedrock layers had
#      ice=liq=0 -> cv=0 -> fact=Inf -> band solve NaN'd the whole column);
#  (4) canopy-fall flux mask mismatch (canopy_hydrology.jl): qflx_snocanfall/liqcanfall/
#      snow_unload were set over mask_soilp but read over mask_nolakep in the fluxes-onto-
#      ground scatter; the glacier (istice) patch is nolakep-but-not-soilp -> its entries
#      stayed NaN -> poisoned qflx_snow_grnd -> int_snow/h2osno.
# (2026-07-12): + the dead-aerosol-InitCold fix above, which only became findable once
# this gate was reshaped to actually form snow layers.
# =============================================================================
include(joinpath(@__DIR__, "..", "test", "testdata.jl"))

using CLM, NCDatasets, Dates, Printf

const GLAC_FS = joinpath(@__DIR__, "..", "test_inputs", "glacier", "surfdata_glacier100.nc")
const BOW_CAL = domain_params_dir("domain_Bow_at_Banff_lumped")
const GLAC_FP = joinpath(BOW_CAL, "clm5_params.nc")
const SNOWOPT = snicar_optics()
const SNOWAGE = snicar_aging()

# Moderate snowfall [mm H2O/s]: 5e-4 * 3600 s = 1.8 kg/m2 per step, which crosses the
# snow-layer-creation threshold on step 1 and builds a multi-layer pack within a day.
# (The old 5e-5 never formed a layer — that is exactly what made this gate blind.)
const SNOWFALL_RATE = 5.0e-4

# Cold, snowy alpine forcing (Bow ~1900 m): sub-freezing air, moderate snowfall.
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
        a.forc_snow_not_downscaled_grc[g] = SNOWFALL_RATE
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

# Finiteness audit restricted to the layers that actually EXIST on each column: the
# port (like Fortran) leaves the unused snow-layer slots at their allocate-to-nan
# value, so a whole-row scan would false-positive on every column. Active layers run
# from nlevsno + snl + 1 (the top snow layer, or the top ground layer when snl == 0)
# to the bottom of the column.
function _active_layer_bad(name, arr, snl, nc, nlevsno)
    for c in 1:nc
        j0 = nlevsno + snl[c] + 1
        for j in j0:size(arr, 2)
            isfinite(arr[c, j]) || return name
        end
    end
    return nothing
end

function main(; nsteps::Int = 24)
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
    nlevsno = CLM.varpar.nlevsno
    @printf("Glacier subgrid: ng=%d nc=%d np=%d  landunit_types=%s\n",
        ng, nc, np, string(Int.(inst.landunit.itype)))

    # locate the glacier (istice) column
    gc = findfirst(c -> inst.landunit.itype[inst.column.landunit[c]] == 4, 1:nc)
    gc === nothing && error("glacier smoke: no istice column in the surfdata")

    config  = CLM.CLMDriverConfig(use_cn=false, use_aquifer_layer=false,
                                  use_hydrstress=false, use_luna=false)
    filt_ia = CLM.clump_filter_inactive_and_active

    nan_fields = String[]
    snl_seen = Int[]
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
        wd = inst.water.waterdiagnosticbulk_inst
        snl = inst.column.snl
        push!(snl_seen, snl[gc])

        # Column / patch scalars — every entry must be finite.
        audited = Dict("t_grnd"=>tg,
                       "h2osno"=>ws.h2osno_no_layers_col,
                       "snow_depth"=>wd.snow_depth_col,
                       "frac_sno_eff"=>wd.frac_sno_eff_col,
                       "eflx_lwrad_out"=>inst.energyflux.eflx_lwrad_out_patch,
                       "eflx_sh_grnd"=>inst.energyflux.eflx_sh_grnd_patch,
                       "qflx_snomelt"=>inst.water.waterfluxbulk_inst.wf.qflx_snomelt_col)
        bad = [k for (k, v) in audited if any(!isfinite, v)]
        # Per-layer arrays — audit the ACTIVE (snow + ground) layers, which is where the
        # snow-layer-formation NaN lived. The old harness never looked at dz/z at all.
        for (k, v) in ("t_soisno"=>inst.temperature.t_soisno_col,
                       "h2osoi_liq"=>ws.h2osoi_liq_col,
                       "h2osoi_ice"=>ws.h2osoi_ice_col,
                       "dz"=>inst.column.dz,
                       "z"=>inst.column.z)
            b = _active_layer_bad(k, v, snl, nc, nlevsno)
            b === nothing || push!(bad, b)
        end
        sort!(bad)

        rng(v) = (vv = filter(isfinite, v); isempty(vv) ? (NaN, NaN) : (minimum(vv), maximum(vv)))
        tgl, tgh = rng(tg)
        @printf("step %2d: snl[glc]=%2d snow_depth=%.4f t_grnd=[%.2f,%.2f] bad_fields=%s\n",
            i, snl[gc], wd.snow_depth_col[gc], tgl, tgh,
            isempty(bad) ? "none" : string(bad))
        isempty(bad) || push!(nan_fields, "step$i:$(bad)")
        if i == 1
            ltypes = [inst.landunit.itype[inst.column.landunit[c]] for c in 1:nc]
            @printf("  [per-column] lun_type=%s  t_grnd_col=%s\n",
                string(Int.(ltypes)), string(tg[1:nc]))
        end
    end

    # The gate is only meaningful if snow layers actually formed. If they never do, the
    # whole snow-layer path went unexercised and a green result proves nothing — so that
    # is reported (and returned) as a FAILURE, not a pass.
    layers_formed = any(<(0), snl_seen)
    min_snl = isempty(snl_seen) ? 0 : minimum(snl_seen)
    @printf("\nsnl history on the glacier column: %s  (deepest snowpack: %d layer(s))\n",
        string(snl_seen), -min_snl)
    if !layers_formed
        println("❌ GLACIER SMOKE BLIND: snow layers NEVER formed (snl stayed 0) — the " *
                "snow-layer path was not exercised, so this run proves nothing. " *
                "Raise SNOWFALL_RATE / nsteps.")
    elseif isempty(nan_fields)
        println("✅ GLACIER SMOKE PASS: driver ran $nsteps steps on a pure-glacier column, " *
                "snow layers formed (down to snl=$min_snl) and every field stayed finite " *
                "across the layer-formation transition")
    else
        println("❌ NaN appeared at: $(nan_fields)")
    end
    return (inst=inst, bounds=bounds, nbad=length(nan_fields),
            layers_formed=layers_formed, min_snl=min_snl, snl_seen=snl_seen)
end

# Returns true iff the glacier domain ran finite for `nsteps` steps AND actually formed
# snow layers (a run that never forms layers is a blind gate, and is reported as a
# FAILURE, not a pass). Returns `missing` if the glacier surfdata / Bow param file is
# absent (gated regression test skips).
function glacier_smoke_ok(; nsteps::Int = 24)
    r = glacier_smoke_result(; nsteps=nsteps)
    r === missing && return missing
    return r.nbad == 0 && r.layers_formed
end

# Same run, but returns the full result so a test can assert on the snl history too.
function glacier_smoke_result(; nsteps::Int = 24)
    (isfile(GLAC_FS) && isfile(GLAC_FP)) ||
        (@info "glacier smoke: inputs absent, skipping" GLAC_FS GLAC_FP; return missing)
    return main(; nsteps=nsteps)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
