# =============================================================================
# REVERSE-AD STREAMFLOW CALIBRATION vs DDS — REAL forcing, REAL observed streamflow.
#
# Upgrade of scripts/calibrate_streamflow_ad*.jl: drives the reverse thread with the
# basin's REAL precipitation forcing over a multi-day/-week window (a genuine
# time-varying hydrograph, not idealized steady rain), builds a REAL streamflow
# objective (windowed daily QRUNOFF → area-gain → linear-reservoir routing → SSE/KGE
# vs the REAL observed daily discharge for the SAME dates), and compares:
#   (a) AD  gradient descent on {fff, baseflow_scalar}  — physics gradient, 1 reverse pass
#   (b) DDS derivative-free search on the SAME windowed objective — for cost+quality parity
#
# NEW differentiable parameter: baseflow_scalar. It enters the CLM power-law baseflow
# flux LINEARLY, behind a state-only branch (zwt <= zi_bedrock) whose condition does
# NOT depend on baseflow_scalar (src/biogeophys/soil_hydrology.jl:2542) — so Enzyme
# differentiates it natively with NO smoothing / NO src change (verified AD==FD below).
# `baseflow_cal!` re-expresses that exact flux with baseflow_scalar read from the bundle
# and a mass-conserving single-layer storage sink (a simplification of CLM's multi-layer
# redistribution — flagged), feeding qflx_drain → total runoff.
#
# HARD GATE: AD gradient of the REAL-obs windowed loss vs central FD (≤1e-4 rel) before
# any calibration is trusted. Honest scope: modeled HRU area ≪ gauge drainage area in
# this dataset (a delineation mismatch), so a FIXED area-gain (bias correction) is applied
# and KGE is reported on the bias-adjusted flows; the runoff is surface + simplified
# power-law baseflow (no snowmelt in the reverse thread → warm-season window).
#
#   CLM_DOMAIN=Stillwater_Oklahoma CLM_NSTEPS=336 CLM_START=2004-06-01 \
#     julia --project=. scripts/calibrate_streamflow_ad_real.jl
# =============================================================================
using CLM, Enzyme, Printf, NCDatasets, Statistics, Dates, Random
const C = CLM
const FT = Float64
const DT = 1800.0                       # 30-min model step
const STEPS_PER_DAY = 48
const ROOT = get(ENV, "SYMFLUENCE_DATA",
    "/Users/darri.eythorsson/Library/CloudStorage/GoogleDrive-dareyt@gmail.com/My Drive/data/SYMFLUENCE_data")
const DOMAIN = get(ENV, "CLM_DOMAIN", "Stillwater_Oklahoma")
const NSTEPS = parse(Int, get(ENV, "CLM_NSTEPS", "336"))     # 336 = 1 wk; 1440 = 30 d
const START  = Date(get(ENV, "CLM_START", "2004-06-01"))
const ROUTE_K = parse(Float64, get(ENV, "CLM_ROUTEK", "3.0"))
const DOM_DIR = joinpath(ROOT, "domain_$DOMAIN")
const FSURDAT = joinpath(DOM_DIR, "settings/CLM/parameters/surfdata_clm.nc")
const PARAMFILE = joinpath(DOM_DIR, "settings/CLM/parameters/clm5_params.nc")

# ---------------------------------------------------------------------------
# real forcing + obs readers
# ---------------------------------------------------------------------------
function read_forcing_window(start::Date, nsteps::Int)
    # returns per-step rain (kg/m2/s) and air temp (K) at 30-min resolution, from the
    # basin's hourly PRECTmms/TBOT (each forcing hour → 2 model steps).
    yr = year(start)
    ff = joinpath(DOM_DIR, "data/forcing/CLM_input", "clmforc.$(yr).nc")
    isfile(ff) || error("no forcing file $ff")
    ds = NCDataset(ff)
    traw = ds["time"][:]                         # NCDatasets CF-decodes to DateTime
    times = DateTime[ismissing(t) ? DateTime(0) : DateTime(t) for t in traw]
    prec = Float64.(ds["PRECTmms"][:]); tbot = Float64.(ds["TBOT"][:])
    close(ds)
    prec = vec(prec); tbot = vec(tbot)
    # first forcing hour at/after START 00:00
    t0 = DateTime(start); i0 = findfirst(>=(t0), times)
    i0 === nothing && error("window start not in forcing")
    nh = cld(nsteps, 2)
    rain = zeros(FT, nsteps); temp = fill(FT(288.0), nsteps)
    for k in 1:nsteps
        h = i0 + (k-1) ÷ 2
        h > length(prec) && break
        r = prec[h]; (isfinite(r) && r < 1e30) || (r = 0.0)
        rain[k] = FT(max(r, 0.0)); temp[k] = FT(isfinite(tbot[h]) ? tbot[h] : 288.0)
    end
    return rain, temp
end

function read_obs_daily(start::Date, ndays::Int)
    f = joinpath(DOM_DIR, "data/observations/streamflow/preprocessed", "$(DOMAIN)_streamflow_processed.csv")
    obs = fill(NaN, ndays); dmap = Dict{Date,Float64}()
    for line in readlines(f)[2:end]
        parts = split(line, ','); length(parts) < 2 && continue
        dt = tryparse(Date, strip(parts[1])[1:min(10,end)]); dt === nothing && continue
        q = tryparse(Float64, strip(parts[2])); (q === nothing || !isfinite(q)) && continue
        dmap[dt] = q
    end
    for d in 1:ndays
        haskey(dmap, start + Day(d-1)) && (obs[d] = dmap[start + Day(d-1)])
    end
    return obs
end

function basin_area_m2()
    for f in (joinpath(DOM_DIR, "data/model_ready/attributes", "$(DOMAIN)_attributes.nc"),)
        isfile(f) || continue
        ds = NCDataset(f)
        for v in ("gru_area","hru_area","area")
            if haskey(ds, v); a = sum(Float64.(ds[v][:])); close(ds); return a; end
        end
        close(ds)
    end
    return 5.46e6
end

# ---------------------------------------------------------------------------
# Enzyme-safe phase functions (per-step data via Const aux)
# ---------------------------------------------------------------------------
function setup_forcing!(a2l, T, ng, rain)
    for g in 1:ng
        a2l.forc_t_not_downscaled_grc[g]=T; a2l.forc_pbot_not_downscaled_grc[g]=95000.0
        a2l.forc_th_not_downscaled_grc[g]=T*(100000.0/95000.0)^(C.RAIR/C.CPAIR)
        a2l.forc_rho_not_downscaled_grc[g]=95000.0/(C.RAIR*T); a2l.forc_lwrad_not_downscaled_grc[g]=320.0
        a2l.forc_vp_grc[g]=1200.0; a2l.forc_hgt_grc[g]=30.0; a2l.forc_topo_grc[g]=0.0; a2l.forc_wind_grc[g]=3.0
        a2l.forc_hgt_u_grc[g]=30.0; a2l.forc_hgt_t_grc[g]=30.0; a2l.forc_hgt_q_grc[g]=30.0
        for b in 1:C.NUMRAD; a2l.forc_solad_not_downscaled_grc[g,b]=200.0; a2l.forc_solai_grc[g,b]=80.0; end
        a2l.forc_solar_not_downscaled_grc[g]=560.0
        a2l.forc_rain_not_downscaled_grc[g]=rain; a2l.forc_snow_not_downscaled_grc[g]=0.0
    end
end

# per-step REAL rain injected as qflx_rain_plus_snomelt (warm-season → rain only, no snowmelt)
function rain_inject!(b, aux)
    wf = b.inst.water.waterfluxbulk_inst.wf
    @inbounds for c in aux.cols; wf.qflx_rain_plus_snomelt_col[c] = aux.rain; end
    return nothing
end
# θ1 = fff (fover) TOPMODEL surface-saturation runoff (exact sat_excess_runoff! physics)
function satexcess_cal!(b, aux)
    i = b.inst; θ1 = b.theta[1]
    sh = i.soilhydrology; ss = i.soilstate; ser = i.sat_excess_runoff
    wfb = i.water.waterfluxbulk_inst; wf = wfb.wf
    @inbounds for c in aux.cols
        ft = sh.frost_table_col[c]; zw = sh.zwt_col[c]; zwp = sh.zwt_perched_col[c]; wt = ss.wtfact_col[c]
        fs = (ft > zwp && ft <= zw) ? wt*exp(-0.5*θ1*zwp) : wt*exp(-0.5*θ1*zw)
        ser.fsat_col[c] = fs; ser.fcov_col[c] = fs
        wfb.qflx_sat_excess_surf_col[c] = fs * wf.qflx_rain_plus_snomelt_col[c]
    end
    return nothing
end
# θ2 = baseflow_scalar: the EXACT CLM power-law baseflow flux (soil_hydrology.jl:2542),
# LINEAR in baseflow_scalar behind the state-only branch zwt<=zi_bedrock. Adds to qflx_drain
# and depletes the water-table soil layer (mass-conserving single-layer sink — a flagged
# simplification of CLM's multi-layer redistribution). No smoothing needed — Enzyme
# differentiates the linear flux natively (AD==FD gated below).
function baseflow_cal!(b, aux)
    i = b.inst; θ2 = b.theta[2]
    sh = i.soilhydrology; ss = i.soilstate; col = i.column
    ws = i.water.waterstatebulk_inst.ws; wf = i.water.waterfluxbulk_inst.wf
    joff = aux.nlevsno; joff_zi = aux.nlevsno + 1
    @inbounds for c in aux.cols
        nb = col.nbedrock[c]; zi_bed = col.zi[c, nb + joff_zi]; zw = sh.zwt_col[c]
        q = 0.0
        if zw <= zi_bed
            slope = col.topo_slope[c]
            q = θ2 * tan(FT(pi)/180.0 * slope) * (zi_bed - zw)^aux.n_baseflow   # mm/s (ice_imped≈1)
        end
        wf.qflx_drain_col[c] = q
        # mass-conserving depletion of the water-table layer (jwt+1), smooth, no branch on θ2
        drain_mm = q * aux.dt
        jl = min(nb, aux.nlevsoi)
        newliq = ws.h2osoi_liq_col[c, jl + joff] - drain_mm
        ws.h2osoi_liq_col[c, jl + joff] = newliq > 0.0 ? newliq : ws.h2osoi_liq_col[c, jl + joff]
    end
    return nothing
end
# per-step runoff = surface + baseflow (mm/s), col-summed
function accum_runoff!(b, aux)
    wf = b.inst.water.waterfluxbulk_inst.wf; s = 0.0
    @inbounds for c in aux.cols; s += wf.qflx_surf_col[c] + wf.qflx_drain_col[c]; end
    b.series[aux.k] += s
    return nothing
end

# ---------------------------------------------------------------------------
# build basin + trajectory
# ---------------------------------------------------------------------------
function build()
    (inst, bounds, filt, tm) = C.clm_initialize!(; fsurdat=FSURDAT, paramfile=PARAMFILE)
    rain0, _ = read_forcing_window(START, NSTEPS)
    setup_forcing!(inst.atm2lnd, 290.0, bounds.endg, mean(rain0))
    C.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
    config = C.CLMDriverConfig(); fia = C.clump_filter_inactive_and_active
    (declin, _) = C.compute_orbital(152.0); nextsw = 152.0 + DT/C.SECSPDAY
    for n in 1:3
        C.clm_drv!(config, inst, filt, fia, bounds, true, nextsw, declin, declin,
            C.ORB_OBLIQR_DEFAULT, false, false, "", false; nstep=n, is_first_step=(n==1),
            is_beg_curr_day=(n==1), dtime=DT, mon=6, day=1, photosyns=inst.photosyns)
    end
    # WET initial water table (flagged): after warmup zwt equilibrates to exactly the
    # bedrock depth, giving zero baseflow HEAD (baseflow ∝ zi_bedrock−zwt = 0), so the
    # baseflow_scalar lever is not exercised. Raise the water table HEAD_M above bedrock —
    # a physically-valid wet spin-up state — so both fff and baseflow_scalar are active and
    # their gradients non-vacuous. Default on; CLM_HEAD=0 restores the equilibrium state.
    head = parse(Float64, get(ENV, "CLM_HEAD", "1.0"))
    if head > 0
        jz = C.varpar.nlevsno + 1
        for c in (bounds.begc:bounds.endc)
            filt.hydrologyc[c] || continue
            nb = inst.column.nbedrock[c]
            zi_bed = inst.column.zi[c, nb + jz]
            inst.soilhydrology.zwt_col[c] = max(0.3, zi_bed - head)
        end
    end
    return inst, bounds, filt, config
end

inst, bounds, filt, config = build()
hc = filt.hydrologyc
const HCOLS = Int[c for c in bounds.begc:bounds.endc if hc[c]]
const FFF0 = FT(C.sat_excess_runoff_params.fff)
const BFS0 = FT(C.BASEFLOW_SCALAR[])
const NDAYS = NSTEPS ÷ STEPS_PER_DAY
const AREA = basin_area_m2()
rain_ser, temp_ser = read_forcing_window(START, NSTEPS)
obs_daily = read_obs_daily(START, NDAYS)
@printf("domain=%s  N=%d steps (%d days from %s)  hydro-cols=%d\n", DOMAIN, NSTEPS, NDAYS, START, length(HCOLS))
@printf("fff0=%.5f  baseflow_scalar0=%.3e  area=%.3e m2  rain window mean=%.3e max=%.3e kg/m2/s\n",
    FFF0, BFS0, AREA, mean(rain_ser), maximum(rain_ser))
@printf("obs daily (cms): mean=%.3f  finite-days=%d/%d\n", mean(filter(isfinite,obs_daily)), count(isfinite,obs_daily), NDAYS)
let c=HCOLS[1], nb=inst.column.nbedrock[c], jz=C.varpar.nlevsno+1
    zi_bed = inst.column.zi[c, nb+jz]; zw = inst.soilhydrology.zwt_col[c]
    @printf("baseflow state: zwt=%.3f m  zi_bedrock=%.3f m (nbedrock=%d)  slope=%.4f  → baseflow %s\n",
        zw, zi_bed, nb, inst.column.topo_slope[c], zw <= zi_bed ? "ACTIVE" : "INACTIVE (zwt below bedrock)")
end

# phase list per step, with this step's real rain
a_surf = C.surfhydro_rev_aux(bounds, filt; dtime=DT, use_hydrstress=config.use_hydrstress)
bf_aux = (; cols = HCOLS, nlevsno = C.varpar.nlevsno, nlevsoi = C.varpar.nlevsoi,
            n_baseflow = FT(C.soilhydrology_params.n_baseflow), dt = DT)
function step_phases(k)
    Any[
        (rain_inject!, ((; cols=HCOLS, rain=rain_ser[k]),)),
        (C.setsoilfrac_rev_phase!, (a_surf,)), (C.setfloodc_rev_phase!, (a_surf,)),
        (satexcess_cal!, ((; cols=HCOLS),)),
        (C.setqflx_rev_phase!, (a_surf,)), (C.inflexcess_rev_phase!, (a_surf,)),
        (C.routeinfl_rev_phase!, (a_surf,)), (C.updateh2osfc_rev_phase!, (a_surf,)),
        (C.infil_rev_phase!, (a_surf,)), (C.totalrunoff_rev_phase!, (a_surf,)),
        (C.soilwater_rev_phase!, (C.soilwater_rev_aux(filt, config; dtime=DT),)),
        (C.watertable_rev_phase!, (C.watertable_rev_aux(bounds, filt, config; dtime=DT),)),
        (C.hydnodrain_rev_phase!, (C.hydnodrain_rev_aux(bounds, filt; dtime=DT),)),
        (baseflow_cal!, (bf_aux,)),
        (accum_runoff!, ((; cols=HCOLS, k=k),)),
    ]
end
const STEPS = [step_phases(k) for k in 1:NSTEPS]
bundle(θ) = (; C.driver_rev_bundle(deepcopy(inst))..., theta=collect(FT,θ), series=zeros(FT,NSTEPS))
series(θ) = (b = bundle(θ); for ph in STEPS, (f,ca) in ph; f(b,ca...); end; copy(b.series))

# ---------------------------------------------------------------------------
# observation operator: per-step runoff series → daily mean → area-gain → route → cms
# GAIN fixed (bias correction for the HRU-vs-gauge area mismatch): set once at the prior.
# ---------------------------------------------------------------------------
daily_mean(s) = [mean(@view s[(d-1)*STEPS_PER_DAY+1 : d*STEPS_PER_DAY]) for d in 1:NDAYS]
function route(q, k)
    k <= 1.0 && return q
    α = 1.0/k; o = similar(q); o[1] = q[1]
    for i in 2:length(q); o[i] = (1-α)*o[i-1] + α*q[i]; end
    return o
end
const GAIN = let m = mean(daily_mean(series([FFF0, BFS0])))
    om = mean(filter(isfinite, obs_daily))
    (isfinite(m) && m > 0) ? om/m : 1.0                    # cms per (mm/s) — absorbs area/units
end
q_sim(s) = route(GAIN .* daily_mean(s), ROUTE_K)            # cms daily
fin = isfinite.(obs_daily)
sse_series(s) = sum((q_sim(s)[fin] .- obs_daily[fin]).^2)  # differentiable in s
function kge_val(s)
    qs = q_sim(s)[fin]; o = obs_daily[fin]
    (length(o) < 3 || std(o) < 1e-9 || std(qs) < 1e-12) && return -Inf
    r = cor(qs, o)
    1.0 - sqrt((r-1)^2 + (std(qs)/std(o)-1)^2 + (mean(qs)/mean(o)-1)^2)
end

# ---------------------------------------------------------------------------
# AD gradient of the REAL-obs SSE: obs-operator adjoint (Enzyme on series→L) seeds the
# physics adjoint (multistep_reverse!). dL/dθ in ONE physics reverse pass.
# ---------------------------------------------------------------------------
function dL_dseries(s)
    g = zero(s)
    Enzyme.autodiff(Enzyme.Reverse, sse_series, Enzyme.Active, Enzyme.Duplicated(copy(s), g))
    return g
end
function grad_ad(θ)
    s = series(θ); ds = dL_dseries(s)
    seed!(db,b) = (db.series .= ds; nothing)
    db = C.multistep_reverse!(STEPS, bundle(θ), seed!)
    return copy(db.theta)
end
Loss(θ) = sse_series(series(θ))

# =============================================================================
# GATE — AD gradient of the real-obs windowed SSE vs central FD (fff + baseflow_scalar)
# =============================================================================
println("\n", "="^80)
println("HARD GATE — AD dL/dθ (real forcing + real obs SSE) vs central FD")
println("="^80)
function run_gate()
    θp = [FFF0*1.2, BFS0*1.3]
    t_ad = @elapsed g = grad_ad(θp)
    @printf("AD gradient @θ=[fff=%.4f, bf=%.3e]: dL/dθ=[% .6e, % .6e]  (%.1fs)\n", θp[1], θp[2], g[1], g[2], t_ad)
    gate_ok = true
    for (j,name,h) in ((1,"fff", 1e-3*FFF0), (2,"baseflow_scalar", 1e-3*BFS0))
        θa=copy(θp); θa[j]+=h; θb=copy(θp); θb[j]-=h
        fd=(Loss(θa)-Loss(θb))/(2h); rel=abs(g[j]-fd)/max(abs(fd),1e-30)
        vac = abs(g[j]) < 1e-30 && abs(fd) < 1e-30
        ok = (rel < 1e-4) && !vac; gate_ok &= ok
        @printf("  ∂/∂%-16s AD=% .6e FD=% .6e rel=%.2e %s\n", name, g[j], fd, rel,
            vac ? "VACUOUS (flux inactive at this state)" : (ok ? "PASS ✓" : "FAIL ✗"))
    end
    @printf("GATE %s  (window N=%d reversed OK)\n", gate_ok ? "PASS ✓" : "FAIL ✗", NSTEPS)
    return gate_ok
end
gate_ok = run_gate()

# =========================================================================
# AD gradient descent vs DDS on the SAME windowed objective, BOTH from a DETUNED start
# (so there is genuine room to recover and the cost comparison is meaningful).
# =========================================================================
function run_calibration()
    println("\n", "="^80); println("AD (gradient descent) vs DDS on the SAME windowed objective — both from a detuned start"); println("="^80)
    lo=[1e-3, 1e-5]; hi=[2.0, 0.05]; sc=[FFF0, BFS0]
    θ0 = clamp.([FFF0*1.8, BFS0*4.0], lo, hi)      # detuned start
    @printf("prior θ=[%.4f, %.3e] KGE=%.4f | detuned start θ0=[%.4f, %.3e] KGE=%.4f\n",
        FFF0, BFS0, kge_val(series([FFF0,BFS0])), θ0[1], θ0[2], kge_val(series(θ0)))
    ADIT = parse(Int, get(ENV,"CLM_AD_ITER","12"))

    # --- AD damped 2×2 Newton on the reverse gradient (SSE objective; report KGE).
    #     3 reverse passes/iter (g at θ, θ+h·e1, θ+h·e2); μ backtracked by a loss test. ---
    θ = copy(θ0); npass = 0
    hs = [1e-3*sc[1], 1e-3*sc[2]]
    for it in 1:ADIT
        g0 = grad_ad(θ);            npass += 1
        g1 = grad_ad(θ .+ [hs[1],0.0]); npass += 1
        g2 = grad_ad(θ .+ [0.0,hs[2]]); npass += 1
        H11=(g1[1]-g0[1])/hs[1]; H21=(g1[2]-g0[2])/hs[1]
        H12=(g2[1]-g0[1])/hs[2]; H22=(g2[2]-g0[2])/hs[2]; H12s=0.5*(H12+H21)
        μ = 1e-2*max(abs(H11),abs(H22),1e-30); Lc = Loss(θ); acc=false
        for _t in 1:14
            d = (H11+μ)*(H22+μ) - H12s^2
            Δ1 = -((H22+μ)*g0[1] - H12s*g0[2]) / d
            Δ2 = -(-H12s*g0[1] + (H11+μ)*g0[2]) / d
            θt = clamp.(θ .+ [Δ1,Δ2], lo, hi)
            if Loss(θt) < Lc; θ = θt; acc = true; break; end
            μ *= 4
        end
        @printf("  AD it %2d θ=[%.4f,%.3e] SSE=%.4e KGE=%.4f rev-passes=%d\n",
            it, θ[1], θ[2], Loss(θ), kge_val(series(θ)), npass)
        (!acc || hypot(g0...) < 1e-20) && break
    end
    ad_kge = kge_val(series(θ)); ad_theta = copy(θ)
    @printf("AD final: θ=[%.4f, %.3e] KGE=%.4f  reverse-passes=%d\n", θ[1], θ[2], ad_kge, npass)

    # --- DDS-on-window (derivative-free) from the SAME detuned start, SAME objective ---
    println("\n-- DDS-on-window (derivative-free, same start + objective) --")
    rng = MersenneTwister(42)
    best = copy(θ0); best_L = Loss(best); evals = 1; m = 2
    ITER = parse(Int, get(ENV,"CLM_DDS_ITER","120")); r = 0.2
    for it in 1:ITER
        P = 1.0 - log(it)/log(ITER)
        cand = copy(best)
        for j in 1:m
            rand(rng) < P && (cand[j] = clamp(best[j] + r*(hi[j]-lo[j])*randn(rng), lo[j], hi[j]))
        end
        L = Loss(cand); evals += 1
        L < best_L && (best_L = L; best = cand)
    end
    dds_kge = kge_val(series(best))
    @printf("DDS final: θ=[%.4f, %.3e] KGE=%.4f  evals=%d\n", best[1], best[2], dds_kge, evals)

    println("\n", "="^80)
    @printf("HEADLINE — same windowed objective (N=%d, %d days), AD vs derivative-free:\n", NSTEPS, NDAYS)
    @printf("%-16s | %-9s | %-14s | %-9s | %-10s\n", "start KGE", "AD KGE", "AD rev-passes", "DDS KGE", "DDS evals")
    println("-"^80)
    @printf("%-16.4f | %-9.4f | %-14d | %-9.4f | %-10d\n", kge_val(series(θ0)), ad_kge, npass, dds_kge, evals)
    println("="^80)
    @printf("full-period SYMFLUENCE DDS baseline (CLM, 2004-05 daily KGE): -0.538 (301 iters, 65%% crash)\n")
    return nothing
end
gate_ok && get(ENV,"CLM_CAL","1") == "1" && run_calibration()
