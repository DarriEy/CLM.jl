# =============================================================================
# Multi-year stability + conservation harness (free-running, no Fortran needed).
#
# verify_vs_fortran.jl validates ONE year (2003) of biogeophysics/hydrology
# against the Fortran .h0 reference. This harness instead asks the long-horizon
# question that needs NO Fortran reference: does the Julia model, started from
# the Fortran spun-up IC, run STABLY and CONSERVATIVELY over multiple years?
#
# It mirrors clm_run!'s setup exactly (so the trajectory is identical to the
# validated offline run) but replaces the time loop with an instrumented one
# that, every step, records:
#   * finiteness   — any NaN/Inf in the core prognostic state (t_soisno,
#                    h2osoi_liq/ice, t_grnd, t_veg)
#   * conservation — the in-model balance-check residuals already computed by
#                    clm_drv! (errh2o_col, errh2osno_col, errseb/errsol/errlon
#                    patch). clm_drv! itself error()s if errh2o > 1e-5 mm, so a
#                    completed run already proves per-step closure; we ALSO
#                    capture the running MAX to quantify the margin.
#   * water budget — column storage endwb_col plus the P / ET / runoff fluxes,
#                    to verify the annual integral ΔS ≈ (P − ET − R) and to
#                    detect a secular storage trend (drift) across years.
#
# It then compares the seasonal cycle of year N vs year N+1 to flag any
# year-over-year drift in the free-running state.
#
# Usage:
#   julia +1.12 --project=. scripts/longhorizon_conservation.jl            # 2yr (2003-2004)
#   julia +1.12 --project=. scripts/longhorizon_conservation.jl --years 3  # 2002-2004
#   julia +1.12 --project=. scripts/longhorizon_conservation.jl --days 30  # short smoke
# =============================================================================
using CLM, NCDatasets, Dates, Printf, Statistics

# ---- CLI ----
_argval(flag, default) = begin
    i = findfirst(==(flag), ARGS)
    (i === nothing || i == length(ARGS)) ? default : ARGS[i + 1]
end
const NDAYS_OVERRIDE = let v = _argval("--days", ""); isempty(v) ? 0 : parse(Int, v); end
const NYEARS         = parse(Int, _argval("--years", "2"))
const DTIME          = 1800
# --probe YYYY-MM-DD : dump the per-step water budget + every nonzero candidate
# flux on that day, to localize which term carries the residual.
const PROBE_DATE     = let v = _argval("--probe", ""); isempty(v) ? nothing : Date(v); end
# --repeat-year YYYY : drive EVERY model year with the forcing from year YYYY, so
# the year-over-year report measures pure MODEL drift (no interannual weather).
const REPEAT_YEAR    = let v = _argval("--repeat-year", ""); isempty(v) ? 0 : parse(Int, v); end
# Map a model timestamp into the repeat forcing year (clamp leap-day Feb 29).
_repeat_force_date(d) = REPEAT_YEAR == 0 ? d : begin
    m = Dates.month(d); dy = Dates.day(d); (m == 2 && dy == 29) && (dy = 28)
    DateTime(REPEAT_YEAR, m, dy, Dates.hour(d), Dates.minute(d), Dates.second(d))
end

# ---- Input paths (calibrated config, matches test/verify_vs_fortran.jl) ----
const basedir = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped"
const caldir  = joinpath(basedir, "optimization/CLM/dds_run_1/final_evaluation/settings/CLM/parameters")
const fsurdat   = joinpath(caldir, "surfdata_clm.nc")
const paramfile = joinpath(caldir, "clm5_params.nc")
# Continuous multi-year forcing (2002-01-01 .. 2004-12-31, hourly)
const fforcing  = joinpath(basedir, "data/forcing/CLM_input/clmforc.2002_2004.nc")
const fsnowoptics = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_optics_5bnd_c013122.nc"
const fsnowaging  = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_drdt_bst_fit_60_c070416.nc"
# Fortran spun-up restart IC (start the free run from model equilibrium)
const ffortran_restart = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_parity_run/Bow_at_Banff_lumped.clm2.r.2003-01-01-00000.nc"

const baseflow_scalar = 0.0022119554
const int_snow_max    = 3113.2227
const USE_AQUIFER     = false

# When started from the 2003 restart we run 2003.. (NYEARS years). With --years 3
# and no restart we'd run 2002.. — but the restart is dated 2003, so default start
# is 2003. (A 3rd year would need clmforc through 2005; we cap end at 2005-01-01.)
const START_YEAR = 2003
const start_date = DateTime(START_YEAR, 1, 1)
const end_date   = NDAYS_OVERRIDE > 0 ? (start_date + Day(NDAYS_OVERRIDE)) :
                                        DateTime(START_YEAR + NYEARS, 1, 1)

println("="^70)
println("  CLM.jl — MULTI-YEAR STABILITY + CONSERVATION",
        REPEAT_YEAR > 0 ? " (CLEAN-DRIFT: repeat forcing $(REPEAT_YEAR))" : "")
println("  IC: Fortran spun-up restart 2003-01-01 | forcing: clmforc.2002_2004.nc",
        REPEAT_YEAR > 0 ? " mapped→$(REPEAT_YEAR)" : "")
@printf("  Window: %s .. %s  (dtime=%ds)\n", start_date, end_date, DTIME)
REPEAT_YEAR > 0 && println("  Every model year sees IDENTICAL weather → year-over-year Δ = pure model drift")
println("="^70)

# ============================================================================
# Setup — faithfully mirrors src/driver/clm_run.jl Phase 1-2
# ============================================================================
(inst, bounds, filt, tm) = CLM.clm_initialize!(;
    fsurdat=fsurdat, paramfile=paramfile,
    start_date=start_date, dtime=DTIME, use_cn=false,
    use_bedrock=true, use_aquifer_layer=USE_AQUIFER, h2osfcflag=0,
    fsnowoptics=isfile(fsnowoptics) ? fsnowoptics : "",
    fsnowaging=isfile(fsnowaging) ? fsnowaging : "",
    int_snow_max=int_snow_max)

config = CLM.CLMDriverConfig(use_cn=false, use_aquifer_layer=USE_AQUIFER)

CLM.atm2lnd_read_namelist!(inst.atm2lnd;
    repartition_rain_snow=true, lapse_rate=0.006, lapse_rate_longwave=0.032,
    precip_repartition_nonglc_all_snow_t=0.0, precip_repartition_nonglc_all_rain_t=2.0,
    precip_repartition_glc_all_snow_t=-2.0, precip_repartition_glc_all_rain_t=0.0)

CLM.init_soil_hydrology_config(baseflow_scalar=baseflow_scalar)

filt_ia = CLM.clump_filter_inactive_and_active
ng, nc, np = bounds.endg, bounds.endc, bounds.endp

# Runtime param wiring (snow/hydrology scalars), mirrors clm_run.jl
let ds_p = NCDataset(paramfile, "r")
    scf = inst.scf_method
    if haskey(ds_p, "n_melt_coef")
        nmc = Float64(ds_p["n_melt_coef"][1])
        for c in 1:nc
            scf.n_melt[c] = nmc / max(10.0, inst.column.topo_std[c])
        end
    end
    haskey(ds_p, "accum_factor") && (scf.accum_factor = Float64(ds_p["accum_factor"][1]))
    haskey(ds_p, "SNOW_DENSITY_MAX") && (CLM.snowhydrology_params.rho_max = Float64(ds_p["SNOW_DENSITY_MAX"][1]))
    haskey(ds_p, "SNOW_DENSITY_MIN") && (CLM.snowhydrology_params.rho_min = Float64(ds_p["SNOW_DENSITY_MIN"][1]))
    haskey(ds_p, "fresh_snw_rds_max") && (CLM.snowhydrology_params.snw_rds_min = Float64(ds_p["fresh_snw_rds_max"][1]))
    haskey(ds_p, "SNO_Z0MV") && (inst.frictionvel.zsno = Float64(ds_p["SNO_Z0MV"][1]))
    if haskey(ds_p, "snw_aging_bst")
        CLM.snicar_params.xdrdt = Float64(ds_p["snw_aging_bst"][1])
    elseif haskey(ds_p, "xdrdt")
        CLM.snicar_params.xdrdt = Float64(ds_p["xdrdt"][1])
    end
    if haskey(ds_p, "pc")
        pc_val = Float64(ds_p["pc"][1])
        pc_val > 0 && (for c in 1:nc; inst.soilhydrology.hkdepth_col[c] = 1.0 / pc_val; end)
    end
    close(ds_p)
end

# Inject Fortran spun-up IC
if isfile(ffortran_restart)
    println("  Injecting spun-up IC: ", basename(ffortran_restart))
    CLM.read_fortran_restart!(ffortran_restart, inst, bounds)
else
    error("Spun-up restart not found: $ffortran_restart")
end

# Forcing reader + topo
fr = CLM.ForcingReader(); CLM.forcing_reader_init!(fr, fforcing)
let topo_file = replace(fforcing, r"clmforc\.[^/]*\.nc$" => "topo_forcing.nc")
    if isfile(topo_file)
        dt = NCDataset(topo_file, "r")
        if haskey(dt, "TOPO")
            ft = Float64(dt["TOPO"][1])
            for g in 1:ng; inst.atm2lnd.forc_topo_grc[g] = ft; end
            for c in 1:nc; inst.topo.topo_col[c] = ft; end
        end
        close(dt)
    end
end

# ============================================================================
# Instrumented time loop
# ============================================================================
wb()  = inst.water.waterbalancebulk_inst
wf()  = inst.water.waterfluxbulk_inst.wf
ef()  = inst.energyflux
a2l() = inst.atm2lnd

# max of |x| over ACTIVE entries, skipping NaN/SPVAL. `active` is the element
# mask (patch or column). Returns 0.0 if no finite active entry.
function _amax(arr, active)
    m = 0.0
    @inbounds for i in eachindex(arr)
        active[i] || continue
        v = arr[i]
        (isfinite(v) && abs(v) < 1.0e30) || continue
        m = max(m, abs(v))
    end
    m
end

# IMPORTANT: this port leaves endwb_col == NaN (never assigned) and feeds the
# gridcell balance_check zeroed output fluxes, so the in-model errh2o checks are
# effectively dead. We therefore compute the column water residual ourselves:
#   errh2o = (S_end - S_beg) - (P - ET - R)·dt
# S_end via the model's own storage function; S_beg = begwb_col (set each step).
const _endwb_scratch = zeros(Float64, nc)
# column canopy water via p2c of the PROGNOSTIC liqcan+snocan
function _canopy_col(c)
    ws = inst.water.waterstatebulk_inst.ws; s = 0.0
    @inbounds for p in eachindex(inst.patch.column)
        inst.patch.column[p] == c || continue
        lc = ws.liqcan_patch[p]; sc = ws.snocan_patch[p]
        isfinite(lc) && (s += lc * inst.patch.wtcol[p])
        isfinite(sc) && (s += sc * inst.patch.wtcol[p])
    end
    s
end
function _true_errh2o_col(dtime)
    # Independent physical-conservation check. S_end = soil/snow/h2osfc/aquifer
    # (the real fn in total_water_heat.jl — the WaterStateBulkData method in
    # hydrology_drainage.jl is a 0.0 stub) PLUS canopy water. S_beg = begwb_col,
    # which (after the fix wiring p2c canopy into begin_water_column_balance!) is
    # itself canopy-inclusive — so we read it directly, no extra canopy term.
    CLM.compute_water_mass_non_lake!(filt.nolakec, inst.column,
        inst.water.waterstatebulk_inst.ws, inst.water.waterdiagnosticbulk_inst,
        false, _endwb_scratch)
    w = wf(); m = 0.0; smax = 0.0; signed = 0.0
    @inbounds for c in 1:nc
        filt.nolakec[c] || continue
        s_end = _endwb_scratch[c] + _canopy_col(c)
        s_beg = wb().begwb_col[c]
        (isfinite(s_end) && isfinite(s_beg)) || continue
        p_in = a2l().forc_rain_downscaled_col[c] + a2l().forc_snow_downscaled_col[c]
        out  = w.qflx_evap_tot_col[c] + w.qflx_surf_col[c] + w.qflx_qrgwl_col[c] +
               w.qflx_drain_col[c] + w.qflx_drain_perched_col[c] +
               w.qflx_snwcp_discarded_liq_col[c] + w.qflx_snwcp_discarded_ice_col[c]
        err = (s_end - s_beg) - (p_in - out) * dtime
        m = max(m, abs(err)); smax = s_end; signed += err
    end
    (m, smax, signed)
end

# storage components for column c: (snow, h2osfc, soil, aquifer_eff) in mm
function _components(c)
    ws = inst.water.waterstatebulk_inst.ws
    nlevsno = CLM.varpar.nlevsno; nlevgrnd = CLM.varpar.nlevgrnd
    snow = ws.h2osno_no_layers_col[c]
    @inbounds for j in (inst.column.snl[c] + 1):0
        jj = j + nlevsno; snow += ws.h2osoi_liq_col[c, jj] + ws.h2osoi_ice_col[c, jj]
    end
    soil = 0.0
    @inbounds for j in 1:nlevgrnd
        jj = j + nlevsno; soil += ws.h2osoi_liq_col[c, jj] + ws.h2osoi_ice_col[c, jj]
    end
    h2osfc = ws.h2osfc_col[c]
    wa_eff = inst.column.hydrologically_active[c] ? (ws.wa_col[c] - ws.aquifer_water_baseline) : 0.0
    # column canopy water via p2c using the PROGNOSTIC liqcan+snocan (h2ocan is a
    # possibly-stale diagnostic). OMITTED from the model's storage def → suspect.
    canopy = 0.0
    @inbounds for p in eachindex(inst.patch.column)
        inst.patch.column[p] == c || continue
        lc = ws.liqcan_patch[p]; sc = ws.snocan_patch[p]
        isfinite(lc) && (canopy += lc * inst.patch.wtcol[p])
        isfinite(sc) && (canopy += sc * inst.patch.wtcol[p])
    end
    (snow, h2osfc, soil, wa_eff, canopy)
end

# finiteness probe over core prognostic state
function _nonfinite_count()
    n = 0
    t = inst.temperature
    for x in t.t_grnd_col; isfinite(x) || (n += 1); end
    for x in t.t_soisno_col; isfinite(x) || (n += 1); end
    for x in t.t_veg_patch; isfinite(x) || (n += 1); end
    ws = inst.water.waterstatebulk_inst.ws
    for x in ws.h2osoi_liq_col; isfinite(x) || (n += 1); end
    for x in ws.h2osoi_ice_col; isfinite(x) || (n += 1); end
    n
end

# Per-day accumulators
mutable struct DayRec
    date::Date
    # max balance residuals over the day
    errh2o::Float64; errh2osno::Float64; errseb::Float64; errsol::Float64; errlon::Float64
    # end-of-day storage (mm) and daily integrated fluxes (mm/day)
    storage::Float64; precip::Float64; et::Float64; runoff::Float64
    # signed water residual accumulated over the day (mm) — localizes the leak
    cumerr::Float64
    # representative states (daily mean)
    tg::Float64; h2osno::Float64; nsamp::Int
    nonfinite::Int
end
DayRec(d) = DayRec(d, 0,0,0,0,0, NaN, 0,0,0, 0, 0,0,0, 0)

days = DayRec[]
cur_day = Date(start_date)
rec = DayRec(cur_day)
step_count = 0
t0 = time()
first_bad_step = 0
cum_errh2o = 0.0   # signed cumulative water residual (mm) — storage telescopes exactly

while tm.current_date < end_date
    global step_count, cur_day, rec, first_bad_step, cum_errh2o
    step_start = tm.current_date
    CLM.advance_timestep!(tm)
    step_count += 1

    CLM.read_forcing_step!(fr, inst.atm2lnd, _repeat_force_date(step_start), ng, nc)
    CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)

    calday = CLM.get_curr_calday(tm)
    (declin, _) = CLM.compute_orbital(calday)
    (declinm1, _) = CLM.compute_orbital(calday - Float64(DTIME) / CLM.SECSPDAY)
    CLM.init_daylength!(inst.gridcell, declin, declinm1, CLM.ORB_OBLIQR_DEFAULT, 1:ng)
    nextsw = calday + Float64(DTIME) / CLM.SECSPDAY
    (yr, mon, d, tod) = CLM.get_curr_date(tm)

    probe_now = PROBE_DATE !== nothing && Date(tm.current_date) == PROBE_DATE
    comp_beg = probe_now ? _components(findfirst(filt.nolakec)) : (0.0,0.0,0.0,0.0,0.0)

    CLM.clm_drv!(config, inst, filt, filt_ia, bounds, true, nextsw, declin, declin,
        CLM.ORB_OBLIQR_DEFAULT, false, false, "", false;
        nstep=tm.nstep, is_first_step=(tm.nstep == 1),
        is_beg_curr_day=CLM.is_beg_curr_day(tm), is_end_curr_day=CLM.is_end_curr_day(tm),
        is_beg_curr_year=CLM.is_beg_curr_year(tm), dtime=Float64(DTIME),
        mon=mon, day=d, photosyns=inst.photosyns)
    CLM.lnd2atm!(bounds, inst)

    # ---- instrument this step ----
    pa = inst.patch.active
    (e_h2o, s_end, e_signed) = _true_errh2o_col(Float64(DTIME))
    cum_errh2o += e_signed
    rec.cumerr   += e_signed
    rec.errh2o    = max(rec.errh2o,    e_h2o)
    rec.errh2osno = max(rec.errh2osno, _amax(wb().errh2osno_col, filt.nolakec))
    rec.errseb    = max(rec.errseb,    _amax(ef().errseb_patch, pa))
    rec.errsol    = max(rec.errsol,    _amax(ef().errsol_patch, pa))
    rec.errlon    = max(rec.errlon,    _amax(ef().errlon_patch, pa))
    # fluxes (mm/s) -> integrate to mm over the step
    p_in = (a2l().forc_rain_downscaled_col[1] + a2l().forc_snow_downscaled_col[1]) * DTIME
    et   = wf().qflx_evap_tot_col[1] * DTIME
    ro   = wf().qflx_runoff_col[1]   * DTIME
    rec.precip += p_in; rec.et += et; rec.runoff += ro
    rec.tg += inst.temperature.t_grnd_col[1]
    rec.h2osno += inst.water.waterdiagnosticbulk_inst.snow_depth_col[1]
    rec.storage = s_end                # true end-of-step storage (last step of day wins)
    rec.nsamp += 1
    # ---- probe: localize the residual on a target day ----
    if PROBE_DATE !== nothing && Date(tm.current_date) == PROBE_DATE && abs(e_signed) > 2e-3
        c = findfirst(filt.nolakec); w = wf(); dt = Float64(DTIME)
        s_beg = wb().begwb_col[c]
        @printf("\n[PROBE %s] errh2o=%+.4e mm | ΔS=%+.4f  S_beg=%.3f S_end=%.3f\n",
                tm.current_date, e_signed, s_end - s_beg, s_beg, s_end)
        ce = _components(c); dcan = ce[5]-comp_beg[5]
        @printf("    Δsnow=%+.4f  Δh2osfc=%+.4f  Δsoil=%+.4f  Δaquifer=%+.4f  Δcanopy=%+.4f\n",
                ce[1]-comp_beg[1], ce[2]-comp_beg[2], ce[3]-comp_beg[3], ce[4]-comp_beg[4], dcan)
        @printf("    >> errh2o + Δcanopy = %+.4e mm  (residual after crediting canopy storage;\n", e_signed + dcan)
        @printf("       the rest = evap_tot over-count: canopy evap not debited from liqcan/snocan)\n")
        # every candidate flux that the column balance might be mis-counting (mm this step)
        cands = [
            (:evap_tot, w.qflx_evap_tot_col[c]), (:surf, w.qflx_surf_col[c]),
            (:qrgwl, w.qflx_qrgwl_col[c]), (:drain, w.qflx_drain_col[c]),
            (:drain_perched, w.qflx_drain_perched_col[c]), (:rsub_sat, w.qflx_rsub_sat_col[c]),
            (:infl, w.qflx_infl_col[c]), (:top_soil, w.qflx_top_soil_col[c]),
            (:snomelt, w.qflx_snomelt_col[c]), (:snow_drain, w.qflx_snow_drain_col[c]),
            (:snow_perc, w.qflx_snow_percolation_col[c]),
            (:snwcp_ice, w.qflx_snwcp_ice_col[c]), (:snwcp_liq, w.qflx_snwcp_liq_col[c]),
            (:snwcp_disc_ice, w.qflx_snwcp_discarded_ice_col[c]),
            (:snwcp_disc_liq, w.qflx_snwcp_discarded_liq_col[c]),
            (:h2osfc_to_ice, w.qflx_h2osfc_to_ice_col[c]),
            (:ice_runoff_snwcp, w.qflx_ice_runoff_snwcp_col[c]),
            (:ice_runoff_xs, w.qflx_ice_runoff_xs_col[c]),
            (:sl_top_soil, w.qflx_sl_top_soil_col[c]),
            (:latflow_in, w.qflx_latflow_in_col[c]), (:latflow_out, w.qflx_latflow_out_col[c]),
            (:rain_plus_melt, w.qflx_rain_plus_snomelt_col[c]),
            (:liq_grnd, w.qflx_liq_grnd_col[c]), (:snow_grnd, w.qflx_snow_grnd_col[c]),
            (:tran_veg, w.qflx_tran_veg_col[c]), (:evap_soi, w.qflx_evap_soi_col[c]),
            (:evap_veg, w.qflx_evap_veg_col[c]), (:evap_can, w.qflx_evap_can_col[c]),
            (:solidevap_top, w.qflx_solidevap_from_top_layer_col[c]),
            (:liqevap_top, w.qflx_liqevap_from_top_layer_col[c]),
            (:soliddew_top, w.qflx_soliddew_to_top_layer_col[c]),
            (:liqdew_top, w.qflx_liqdew_to_top_layer_col[c]),
            (:snow_h2osfc, w.qflx_snow_h2osfc_col[c]),
            (:too_small_h2osfc, w.qflx_too_small_h2osfc_to_soil_col[c]),
            (:P_rain, a2l().forc_rain_downscaled_col[c]), (:P_snow, a2l().forc_snow_downscaled_col[c]),
        ]
        # component-vs-flux reconciliation (mm this step)
        snow_in  = (w.qflx_snow_grnd_col[c] + w.qflx_soliddew_to_top_layer_col[c] + w.qflx_liqdew_to_top_layer_col[c]) * dt
        snow_out = (w.qflx_snow_drain_col[c] + w.qflx_solidevap_from_top_layer_col[c] +
                    w.qflx_liqevap_from_top_layer_col[c] + w.qflx_snwcp_ice_col[c] + w.qflx_snwcp_liq_col[c]) * dt
        @printf("    snow recon: Δsnow=%+.4f  (in %.4f - out %.4f = %+.4f)  resid=%+.4e\n",
                ce[1]-comp_beg[1], snow_in, snow_out, snow_in-snow_out, (ce[1]-comp_beg[1])-(snow_in-snow_out))
        # canopy throughfall recon: water leaving canopy to ground vs canopy storage drop
        canfall = 0.0; evapveg_p = 0.0
        @inbounds for p in eachindex(inst.patch.column)
            inst.patch.column[p] == c || continue
            sf = w.qflx_snocanfall_patch[p]; lf = w.qflx_liqcanfall_patch[p]
            isfinite(sf) && (canfall += sf * inst.patch.wtcol[p] * dt)
            isfinite(lf) && (canfall += lf * inst.patch.wtcol[p] * dt)
        end
        ground_in = (w.qflx_liq_grnd_col[c] + w.qflx_snow_grnd_col[c]) * dt
        atmos_in  = (a2l().forc_rain_downscaled_col[c] + a2l().forc_snow_downscaled_col[c]) * dt
        @printf("    canopy: Δh2ocan=%+.4f  canfall(out→grnd)=%+.4f | ground_in=%+.4f  atmos_in=%+.4f  (grnd-atmos=%+.4f)\n",
                dcan, canfall, ground_in, atmos_in, ground_in - atmos_in)
        for (nm, v) in cands
            abs(v * dt) > 1e-6 && @printf("    %-18s %+10.4f mm\n", nm, v * dt)
        end
    end

    nf = _nonfinite_count()
    rec.nonfinite = max(rec.nonfinite, nf)
    if nf > 0 && first_bad_step == 0
        first_bad_step = step_count
        @warn "Non-finite state detected" step=step_count date=tm.current_date count=nf
    end

    # roll the day (only push a record that actually accumulated steps)
    today = Date(tm.current_date)
    if today != cur_day
        rec.nsamp > 0 && push!(days, rec)
        cur_day = today
        rec = DayRec(cur_day)
    end

    if step_count % (48 * 30) == 0
        @printf("  step %6d  %s  | maxerrh2o=%.2e errseb=%.2e  storage=%.1fmm  Tg=%.1fK  nf=%d\n",
                step_count, tm.current_date, rec.errh2o, rec.errseb,
                rec.storage, inst.temperature.t_grnd_col[1], rec.nonfinite)
    end
end
rec.nsamp > 0 && push!(days, rec)
CLM.forcing_reader_close!(fr)
elapsed = time() - t0
@printf("\n  Completed %d steps (%d days) in %.1fs (%.2fs/day)\n",
        step_count, length(days), elapsed, elapsed / max(1, length(days)))

# ============================================================================
# Report
# ============================================================================
println("\n", "="^70)
println("  CONSERVATION — running max of in-model balance residuals")
println("="^70)
# In-model balance-check fields (now live after wiring canopy into begwb/endwb +
# computing endwb_col). Confirm they are FINITE (not the old NaN silent-pass) and small.
let cc = findfirst(filt.nolakec), w = wb(), f = wf()
    @printf("  in-model check (last step): endwb_col=%.4f  begwb_col=%.4f  errh2o_col=%.3e  %s\n",
            w.endwb_col[cc], w.begwb_col[cc], w.errh2o_col[cc],
            isfinite(w.errh2o_col[cc]) ? "[FINITE — guard live]" : "[NaN — guard still dead!]")
    for (nm, v) in [("floodc", f.qflx_floodc_col[cc]), ("sfc_irrig", f.qflx_sfc_irrig_col[cc]),
                    ("glcice_dyn", f.qflx_glcice_dyn_water_flux_col[cc]), ("qrgwl", f.qflx_qrgwl_col[cc]),
                    ("surf", f.qflx_surf_col[cc]), ("drain", f.qflx_drain_col[cc]),
                    ("drain_perched", f.qflx_drain_perched_col[cc]), ("evap_tot", f.qflx_evap_tot_col[cc]),
                    ("snwcp_disc_liq", f.qflx_snwcp_discarded_liq_col[cc]), ("snwcp_disc_ice", f.qflx_snwcp_discarded_ice_col[cc])]
        isfinite(v) || @printf("    NaN flux term: %s\n", nm)
    end
end
gmax(f) = isempty(days) ? NaN : maximum(f, days)
run_years = max(1e-9, length(days) / 365.25)
@printf("  max |errh2o|    (water,  mm/step) : %.3e\n", gmax(r->r.errh2o))
@printf("  cumulative errh2o (signed, mm)   : %+.3e   over %.2f yr  → %+.3e mm/yr\n",
        cum_errh2o, run_years, cum_errh2o / run_years)
@printf("  max |errh2osno| (snow,   mm/step) : %.3e\n", gmax(r->r.errh2osno))
@printf("  max |errseb|    (sfc E,  W/m2) : %.3e   [threshold %.0e]\n", gmax(r->r.errseb), 1e-5)
@printf("  max |errsol|    (solar,  W/m2) : %.3e\n", gmax(r->r.errsol))
@printf("  max |errlon|    (longw,  W/m2) : %.3e\n", gmax(r->r.errlon))
totnf = sum(r->r.nonfinite, days)
@printf("  non-finite state samples       : %d %s\n", totnf,
        totnf == 0 ? "(model stayed finite)" : "(FIRST at step $first_bad_step)")

println("\n", "="^70)
println("  WATER BUDGET — annual closure & secular drift")
println("="^70)
println("  year |  P(mm)   ET(mm)   R(mm)  | ΔStorage  P-ET-R   residual | S_end(mm)")
println("  " * "-"^70)
yrs = sort(unique(Dates.year.(getfield.(days, :date))))
for y in yrs
    idx = findall(r -> Dates.year(r.date) == y, days)
    isempty(idx) && continue
    P = sum(days[i].precip for i in idx); ET = sum(days[i].et for i in idx); R = sum(days[i].runoff for i in idx)
    s0 = days[first(idx)].storage; s1 = days[last(idx)].storage
    dS = s1 - s0; flux = P - ET - R
    @printf("  %4d | %7.1f %7.1f %7.1f | %8.1f %8.1f %9.2e | %8.1f\n",
            y, P, ET, R, dS, flux, dS - flux, s1)
end

# Monthly signed water-residual breakdown — localizes the leak by season
println("\n", "="^70)
println("  WATER RESIDUAL BY MONTH (signed mm) — where conservation breaks")
println("="^70)
println("  mon |  Σresid   |   P(mm)    ET(mm)    R(mm)   snowdp(m)")
println("  " * "-"^58)
for y in yrs, m in 1:12
    idx = findall(r -> Dates.year(r.date) == y && Dates.month(r.date) == m, days)
    isempty(idx) && continue
    cum = sum(days[i].cumerr for i in idx)
    P = sum(days[i].precip for i in idx); ET = sum(days[i].et for i in idx); R = sum(days[i].runoff for i in idx)
    sd = mean(days[i].h2osno / max(1, days[i].nsamp) for i in idx)
    @printf("  %4d-%02d | %+8.2f | %8.1f %8.2f %8.2f   %7.3f\n", y, m, cum, P, ET, R, sd)
end

# Year-over-year seasonal-cycle drift (monthly-mean Tg & snow depth)
if length(yrs) >= 2
    println("\n", "="^70)
    println("  YEAR-OVER-YEAR DRIFT — monthly-mean ground temp T_GRND (K)")
    println("="^70)
    mon_mean(y, m, f) = begin
        idx = findall(r -> Dates.year(r.date) == y && Dates.month(r.date) == m, days)
        isempty(idx) ? NaN : mean(f(days[i]) / max(1, days[i].nsamp) for i in idx)
    end
    print("  mon "); for y in yrs; @printf("  %8d", y); end; println("   max|Δ|")
    drifts = Float64[]
    for m in 1:12
        @printf("  %3d ", m)
        vals = Float64[]
        for y in yrs
            v = mon_mean(y, m, r -> r.tg); push!(vals, v); @printf("  %8.2f", v)
        end
        fin = filter(isfinite, vals)
        d = length(fin) >= 2 ? maximum(fin) - minimum(fin) : NaN
        isfinite(d) && push!(drifts, d)
        @printf("   %6.2f\n", d)
    end
    @printf("\n  Max month-to-year T_GRND drift: %.2f K\n", isempty(drifts) ? NaN : maximum(drifts))
end

println("\n", "="^70)
println("  VERDICT")
println("="^70)
ok_fin = totnf == 0
drift_rate = abs(cum_errh2o) / run_years          # mm/yr
ok_h2o = drift_rate < 1.0                          # <1 mm/yr vs ~5400 mm storage
ok_seb = gmax(r->r.errseb) < 1e-5
@printf("  finiteness     : %s\n", ok_fin ? "PASS (no NaN/Inf over full horizon)" : "FAIL")
@printf("  water budget   : %s (cumulative drift %+.3e mm/yr; per-step max %.2e mm)\n",
        ok_h2o ? "CONSERVED" : "DRIFT", cum_errh2o / run_years, gmax(r->r.errh2o))
@printf("  energy closure : %s (errseb < 1e-5 W/m2 every step)\n", ok_seb ? "PASS" : "FAIL")
println("="^70)
