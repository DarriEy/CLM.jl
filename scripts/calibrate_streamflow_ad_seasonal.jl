# =============================================================================
# SEASONAL SPLIT-SAMPLE streamflow calibration (Klemeš) via CHECKPOINTED reverse-AD.
#
# Extends calibrate_streamflow_ad_real.jl past the ~30-day flat-tape wall using
# multistep_reverse_binomial! (O(log N) coarse-checkpoint memory, O(N log N) recompute),
# so a ~90-day season reverses without OOM. Same real-forcing hydrology thread, same
# {fff, baseflow_scalar} params, same area-gain/routing/KGE observation operator as _real.jl.
#
# STAGE A (CLM_STAGE=A, default) — HARD gate BEFORE any long descent:
#   1. 30-day: flat multistep_reverse! == binomial == central FD  (binomial correctness +
#      identical-to-flat where both fit).
#   2. 90-day: binomial == central FD  (correctness at the new scale) + wall-time/gradient
#      + peak coarse-checkpoint count (memory-wall gone).
# STAGE B (CLM_STAGE=B, only after A is green) — the split-sample:
#   cal season (year Y) → checkpointed-reverse gradient descent of {fff, baseflow_scalar}
#   → eval on the SAME season in year Y+1 (independent). AD vs DDS; report cal/eval KGE.
#
#   CLM_STAGE=A CLM_NDAYS=90 CLM_CAL_START=2004-06-01 julia --project=. \
#       scripts/calibrate_streamflow_ad_seasonal.jl
# =============================================================================
using CLM, Enzyme, ForwardDiff, Printf, NCDatasets, Statistics, Dates, Random
const C = CLM
const FT = Float64
const DT = 1800.0
const STEPS_PER_DAY = 48
const ROOT = get(ENV, "SYMFLUENCE_DATA",
    "/Users/darri.eythorsson/Library/CloudStorage/GoogleDrive-dareyt@gmail.com/My Drive/data/SYMFLUENCE_data")
const DOMAIN = get(ENV, "CLM_DOMAIN", "Stillwater_Oklahoma")
const NDAYS = parse(Int, get(ENV, "CLM_NDAYS", "90"))
const NSTEPS = NDAYS * STEPS_PER_DAY
const ROUTE_K = parse(Float64, get(ENV, "CLM_ROUTEK", "3.0"))
const STAGE = get(ENV, "CLM_STAGE", "A")
const DOM_DIR = joinpath(ROOT, "domain_$DOMAIN")
const FSURDAT = joinpath(DOM_DIR, "settings/CLM/parameters/surfdata_clm.nc")
const PARAMFILE = joinpath(DOM_DIR, "settings/CLM/parameters/clm5_params.nc")

# ---------------------------------------------------------------------------
# readers (identical operator to _real.jl)
# ---------------------------------------------------------------------------
function read_forcing_window(start::Date, nsteps::Int)
    yr = year(start); ff = joinpath(DOM_DIR, "data/forcing/CLM_input", "clmforc.$(yr).nc")
    isfile(ff) || error("no forcing file $ff")
    ds = NCDataset(ff)
    traw = ds["time"][:]; times = DateTime[ismissing(t) ? DateTime(0) : DateTime(t) for t in traw]
    prec = vec(Float64.(ds["PRECTmms"][:])); close(ds)
    t0 = DateTime(start); i0 = findfirst(>=(t0), times); i0 === nothing && error("start not in forcing")
    rain = zeros(FT, nsteps)
    for k in 1:nsteps
        h = i0 + (k-1) ÷ 2; h > length(prec) && break
        r = prec[h]; rain[k] = FT((isfinite(r) && r < 1e30) ? max(r,0.0) : 0.0)
    end
    return rain
end
function read_obs_daily(start::Date, ndays::Int)
    f = joinpath(DOM_DIR, "data/observations/streamflow/preprocessed", "$(DOMAIN)_streamflow_processed.csv")
    obs = fill(NaN, ndays); ssum = Dict{Date,Float64}(); scnt = Dict{Date,Int}()
    for line in readlines(f)[2:end]                    # AVERAGE sub-daily readings → daily mean
        p = split(line, ','); length(p) < 2 && continue
        dt = tryparse(Date, strip(p[1])[1:min(10,end)]); dt === nothing && continue
        q = tryparse(Float64, strip(p[2])); (q === nothing || !isfinite(q)) && continue
        ssum[dt] = get(ssum, dt, 0.0) + q; scnt[dt] = get(scnt, dt, 0) + 1
    end
    for d in 1:ndays
        k = start+Day(d-1); haskey(scnt, k) && (obs[d] = ssum[k]/scnt[k])
    end
    return obs
end
function basin_area_m2()
    f = joinpath(DOM_DIR, "data/model_ready/attributes", "$(DOMAIN)_attributes.nc")
    isfile(f) || return 5.46e6
    ds = NCDataset(f)
    for v in ("gru_area","hru_area","area")
        haskey(ds, v) && (a = sum(Float64.(ds[v][:])); close(ds); return a)
    end
    close(ds); return 5.46e6
end

# ---------------------------------------------------------------------------
# phase functions (identical to _real.jl)
# ---------------------------------------------------------------------------
function setup_forcing!(a2l, T, ng, rain)
    for g in 1:ng
        a2l.forc_t_not_downscaled_grc[g]=T; a2l.forc_pbot_not_downscaled_grc[g]=95000.0
        a2l.forc_th_not_downscaled_grc[g]=T*(100000.0/95000.0)^(C.RAIR/C.CPAIR)
        a2l.forc_rho_not_downscaled_grc[g]=95000.0/(C.RAIR*T); a2l.forc_lwrad_not_downscaled_grc[g]=320.0
        a2l.forc_vp_grc[g]=1200.0; a2l.forc_hgt_grc[g]=30.0; a2l.forc_topo_grc[g]=0.0; a2l.forc_wind_grc[g]=3.0
        a2l.forc_hgt_u_grc[g]=30.0; a2l.forc_hgt_t_grc[g]=30.0; a2l.forc_hgt_q_grc[g]=30.0
        for b in 1:C.NUMRAD; a2l.forc_solad_not_downscaled_grc[g,b]=200.0; a2l.forc_solai_grc[g,b]=80.0; end
        a2l.forc_solar_not_downscaled_grc[g]=560.0; a2l.forc_rain_not_downscaled_grc[g]=rain; a2l.forc_snow_not_downscaled_grc[g]=0.0
    end
end
function rain_inject!(b, aux)
    wf = b.inst.water.waterfluxbulk_inst.wf
    @inbounds for c in aux.cols; wf.qflx_rain_plus_snomelt_col[c] = aux.rain; end
    return nothing
end
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
function baseflow_cal!(b, aux)
    i = b.inst; θ2 = b.theta[2]
    sh = i.soilhydrology; col = i.column
    ws = i.water.waterstatebulk_inst.ws; wf = i.water.waterfluxbulk_inst.wf
    jz = aux.nlevsno + 1
    @inbounds for c in aux.cols
        nb = col.nbedrock[c]; zi_bed = col.zi[c, nb + jz]; zw = sh.zwt_col[c]
        q = zw <= zi_bed ? θ2 * tan(FT(pi)/180.0 * col.topo_slope[c]) * (zi_bed - zw)^aux.n_baseflow : zero(zw)
        wf.qflx_drain_col[c] = q
        jl = min(nb, aux.nlevsoi); nl = ws.h2osoi_liq_col[c, jl + aux.nlevsno] - q*aux.dt
        ws.h2osoi_liq_col[c, jl + aux.nlevsno] = nl > 0 ? nl : ws.h2osoi_liq_col[c, jl + aux.nlevsno]
    end
    return nothing
end
function accum_runoff!(b, aux)
    wf = b.inst.water.waterfluxbulk_inst.wf; s = zero(eltype(b.series))
    @inbounds for c in aux.cols; s += wf.qflx_surf_col[c] + wf.qflx_drain_col[c]; end
    b.series[aux.k] += s
    return nothing
end

# ---------------------------------------------------------------------------
# build a season context: warmed inst + real forcing + obs + STEPS + operator closures.
# gain_fixed: if given, use that area-gain (cal→eval); else compute from this window's prior.
# ---------------------------------------------------------------------------
function make_season(start::Date; ndays::Int, gain_fixed::Union{Float64,Nothing}=nothing)
    nsteps = ndays * STEPS_PER_DAY
    (inst, bounds, filt, tm) = C.clm_initialize!(; fsurdat=FSURDAT, paramfile=PARAMFILE)
    rain_ser = read_forcing_window(start, nsteps)
    setup_forcing!(inst.atm2lnd, 290.0, bounds.endg, mean(rain_ser))
    C.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
    config = C.CLMDriverConfig(); fia = C.clump_filter_inactive_and_active
    doy = dayofyear(start); (declin,_) = C.compute_orbital(Float64(doy)); nextsw = doy + DT/C.SECSPDAY
    for n in 1:3
        C.clm_drv!(config, inst, filt, fia, bounds, true, nextsw, declin, declin,
            C.ORB_OBLIQR_DEFAULT, false, false, "", false; nstep=n, is_first_step=(n==1),
            is_beg_curr_day=(n==1), dtime=DT, mon=month(start), day=1, photosyns=inst.photosyns)
    end
    # Wet water-table HEAD above bedrock (flagged): head>0 activates baseflow. A LARGE head
    # makes runoff baseflow-dominated (smooth); a SMALL head lets rain-driven surface runoff
    # dominate the timing (better correlation with flashy obs). Tunable via CLM_HEAD.
    head = parse(Float64, get(ENV, "CLM_HEAD", "1.0"))
    if head > 0
        jz = C.varpar.nlevsno + 1
        for c in bounds.begc:bounds.endc
            filt.hydrologyc[c] || continue
            nb = inst.column.nbedrock[c]
            inst.soilhydrology.zwt_col[c] = max(0.3, inst.column.zi[c, nb+jz] - head)
        end
    end
    hcols = Int[c for c in bounds.begc:bounds.endc if filt.hydrologyc[c]]
    obs_daily = read_obs_daily(start, ndays)

    a_surf = C.surfhydro_rev_aux(bounds, filt; dtime=DT, use_hydrstress=config.use_hydrstress)
    bf_aux = (; cols=hcols, nlevsno=C.varpar.nlevsno, nlevsoi=C.varpar.nlevsoi,
                n_baseflow=FT(C.soilhydrology_params.n_baseflow), dt=DT)
    step_phases(k) = Any[
        (rain_inject!, ((; cols=hcols, rain=rain_ser[k]),)),
        (C.setsoilfrac_rev_phase!, (a_surf,)), (C.setfloodc_rev_phase!, (a_surf,)),
        (satexcess_cal!, ((; cols=hcols),)),
        (C.setqflx_rev_phase!, (a_surf,)), (C.inflexcess_rev_phase!, (a_surf,)),
        (C.routeinfl_rev_phase!, (a_surf,)), (C.updateh2osfc_rev_phase!, (a_surf,)),
        (C.infil_rev_phase!, (a_surf,)), (C.totalrunoff_rev_phase!, (a_surf,)),
        (C.soilwater_rev_phase!, (C.soilwater_rev_aux(filt, config; dtime=DT),)),
        (C.watertable_rev_phase!, (C.watertable_rev_aux(bounds, filt, config; dtime=DT),)),
        (C.hydnodrain_rev_phase!, (C.hydnodrain_rev_aux(bounds, filt; dtime=DT),)),
        (baseflow_cal!, (bf_aux,)),
        (accum_runoff!, ((; cols=hcols, k=k),))]
    steps = [step_phases(k) for k in 1:nsteps]
    bundle(θ) = (; C.driver_rev_bundle(deepcopy(inst))..., theta=collect(FT,θ), series=zeros(FT,nsteps))
    series(θ) = (b = bundle(θ); for ph in steps, (f,ca) in ph; f(b,ca...); end; copy(b.series))

    daily_mean(s) = [mean(@view s[(d-1)*STEPS_PER_DAY+1 : d*STEPS_PER_DAY]) for d in 1:ndays]
    function route(q, k)
        k <= 1.0 && return q
        α=1.0/k; o=similar(q); o[1]=q[1]; for i in 2:length(q); o[i]=(1-α)*o[i-1]+α*q[i]; end; o
    end
    gain = gain_fixed !== nothing ? gain_fixed : begin
        m = mean(daily_mean(series([FT(C.sat_excess_runoff_params.fff), FT(C.BASEFLOW_SCALAR[])])))
        om = mean(filter(isfinite, obs_daily))
        (isfinite(m) && m > 0 && isfinite(om)) ? om/m : 1.0
    end
    q_sim(s) = route(gain .* daily_mean(s), ROUTE_K)
    fin = isfinite.(obs_daily)
    sse_series(s) = sum((q_sim(s)[fin] .- obs_daily[fin]).^2)
    # KGE with an ARBITRARY gain g (so we can report both cal-fixed-gain (strict) and
    # own-gain (bias-corrected shape) eval-KGE).
    function kge_g(s, g)
        qs = route(g .* daily_mean(s), ROUTE_K)[fin]; o = obs_daily[fin]
        (length(o) < 3 || std(o) < 1e-9 || std(qs) < 1e-12) && return -Inf
        r = cor(qs, o); 1.0 - sqrt((r-1)^2 + (std(qs)/std(o)-1)^2 + (mean(qs)/mean(o)-1)^2)
    end
    kge_val(s) = kge_g(s, gain)
    # this window's OWN bias-correcting gain (mean obs / mean model at prior), for the
    # bias-corrected shape-generalization metric.
    own_gain(s) = (m=mean(daily_mean(s)); om=mean(o for o in obs_daily if isfinite(o)); (isfinite(m)&&m>0) ? om/m : gain)
    # EXACT adjoint of the LINEAR obs operator L = Σ_d fin (gain·route(daily(s)) − obs)².
    # (Hand-coded — the operator is linear, and Enzyme trips analyzing a closure that
    # captures a BitVector; ForwardDiff would need N partials. Validated vs central FD in
    # Stage A.) Chain: series → daily-mean → gain → linear-reservoir route → SSE.
    function dL_dseries(s)
        qs = q_sim(s)
        dqs = zeros(ndays)
        @inbounds for d in 1:ndays; fin[d] && (dqs[d] = 2*(qs[d]-obs_daily[d])); end
        α = ROUTE_K <= 1.0 ? 1.0 : 1.0/ROUTE_K
        dqg = zeros(ndays)
        if ROUTE_K <= 1.0
            dqg .= dqs
        else
            adj = copy(dqs)
            @inbounds for i in ndays:-1:2; dqg[i] += α*adj[i]; adj[i-1] += (1-α)*adj[i]; end
            dqg[1] += adj[1]
        end
        ds = zeros(length(s))
        @inbounds for d in 1:ndays
            g = gain * dqg[d] / STEPS_PER_DAY
            for k in (d-1)*STEPS_PER_DAY+1 : d*STEPS_PER_DAY; ds[k] = g; end
        end
        return ds
    end
    # gradient of the windowed SSE; method = :flat | :binomial
    function grad_ad(θ; method::Symbol=:binomial, peak::Union{Nothing,Base.RefValue{Int}}=nothing)
        s = series(θ); ds = dL_dseries(s)
        seed!(db,b) = (db.series .= ds; nothing)
        db = method === :flat ? C.multistep_reverse!(steps, bundle(θ), seed!) :
                                C.multistep_reverse_binomial!(steps, bundle(θ), seed!; peak_checkpoints=peak)
        return copy(db.theta)
    end
    Loss(θ) = sse_series(series(θ))
    return (; start, ndays, nsteps, hcols=length(hcols), obs=obs_daily, gain,
              series, grad_ad, Loss, kge_val, kge_g, own_gain, sse_series, nfin=count(fin))
end

# One throwaway init to populate the module param globals (fff/baseflow_scalar) so the
# constants below are the basin's real defaults, not the pre-init NaN.
let; C.clm_initialize!(; fsurdat=FSURDAT, paramfile=PARAMFILE); end
const FFF0 = FT(C.sat_excess_runoff_params.fff)
const BFS0 = FT(C.BASEFLOW_SCALAR[])
const CAL_START  = Date(get(ENV, "CLM_CAL_START",  "2004-06-01"))
const EVAL_START = Date(get(ENV, "CLM_EVAL_START", "2005-06-01"))
@printf("SEASONAL split-sample: domain=%s  NDAYS=%d (N=%d steps)  fff0=%.4f bf0=%.3e\n",
    DOMAIN, NDAYS, NSTEPS, FFF0, BFS0)

# =============================================================================
# STAGE A — checkpointed-gradient gate
# =============================================================================
function stage_a()
    println("\n", "="^82); println("STAGE A — checkpointed (binomial) reverse-AD gradient gate"); println("="^82)
    θp = [FFF0*1.2, BFS0*1.3]

    # (1) 30-day: flat == binomial == central FD
    println("[A1] 30-day: flat vs binomial vs central FD")
    c30 = make_season(CAL_START; ndays=30)
    gflat = c30.grad_ad(θp; method=:flat)
    gbin  = c30.grad_ad(θp; method=:binomial)
    @printf("     flat     dL/dθ = [% .6e, % .6e]\n", gflat[1], gflat[2])
    @printf("     binomial dL/dθ = [% .6e, % .6e]\n", gbin[1], gbin[2])
    ok_a1 = true
    for (j,name,h) in ((1,"fff",1e-3*FFF0),(2,"baseflow_scalar",1e-3*BFS0))
        θa=copy(θp); θa[j]+=h; θb=copy(θp); θb[j]-=h
        fd=(c30.Loss(θa)-c30.Loss(θb))/(2h)
        rel_bf=abs(gbin[j]-fd)/max(abs(fd),1e-30); rel_fb=abs(gbin[j]-gflat[j])/max(abs(gflat[j]),1e-30)
        ok = rel_bf<1e-5 && rel_fb<1e-8; ok_a1 &= ok
        @printf("     %-16s bin-vs-FD rel=%.2e  bin-vs-flat rel=%.2e  %s\n", name, rel_bf, rel_fb, ok ? "PASS ✓" : "FAIL ✗")
    end

    # (2) 90-day: binomial == central FD, wall-time, peak checkpoints
    @printf("\n[A2] %d-day: binomial vs central FD + wall-time + peak coarse-checkpoints\n", NDAYS)
    cN = make_season(CAL_START; ndays=NDAYS)
    peak = Ref(0)
    t_grad = @elapsed (gbinN = cN.grad_ad(θp; method=:binomial, peak=peak))
    @printf("     binomial dL/dθ = [% .6e, % .6e]   wall=%.1f s/gradient   peak coarse-checkpoints=%d (of %d steps)\n",
        gbinN[1], gbinN[2], t_grad, peak[], NSTEPS)
    ok_a2 = true
    for (j,name,h) in ((1,"fff",1e-3*FFF0),(2,"baseflow_scalar",1e-3*BFS0))
        θa=copy(θp); θa[j]+=h; θb=copy(θp); θb[j]-=h
        fd=(cN.Loss(θa)-cN.Loss(θb))/(2h); rel=abs(gbinN[j]-fd)/max(abs(fd),1e-30)
        vac = abs(gbinN[j])<1e-30 && abs(fd)<1e-30; ok = rel<1e-5 && !vac; ok_a2 &= ok
        @printf("     %-16s bin-vs-FD rel=%.2e  %s\n", name, rel, vac ? "VACUOUS" : (ok ? "PASS ✓" : "FAIL ✗"))
    end
    @printf("\nSTAGE A %s  (%d-day reversed OK, wall≈%.0fs/gradient)\n",
        (ok_a1 && ok_a2) ? "PASS ✓" : "FAIL ✗", NDAYS, t_grad)
    return ok_a1 && ok_a2, t_grad
end

# =============================================================================
# STAGE B — seasonal split-sample (cal year Y → eval year Y+1, same season)
# =============================================================================
function stage_b()
    println("\n", "="^82); println("STAGE B — seasonal split-sample (AD vs DDS)"); println("="^82)
    @printf("cal season  = %s .. %s\n", CAL_START, CAL_START + Day(NDAYS-1))
    @printf("eval season = %s .. %s (independent, held out)\n", EVAL_START, EVAL_START + Day(NDAYS-1))
    cal = make_season(CAL_START; ndays=NDAYS)
    eval = make_season(EVAL_START; ndays=NDAYS, gain_fixed=cal.gain)   # fixed area-gain from cal
    lo=[1e-3,1e-5]; hi=[2.0,0.05]; sc=[FFF0,BFS0]
    dscale = parse(Float64, get(ENV,"CLM_DETUNE","1.4"))              # detune magnitude (milder default)
    θprior=[FFF0,BFS0]; θ0=clamp.([FFF0*dscale, BFS0*dscale^2], lo, hi)
    # eval-KGE reported two ways: STRICT (cal-fixed gain — conservative, sees inter-annual mean
    # drift) and BIAS-CORRECTED (eval's own mean gain — isolates SHAPE/param generalization).
    eval_kge(θ) = (s=eval.series(θ); (eval.kge_g(s, cal.gain), eval.kge_g(s, eval.own_gain(s))))
    pk = eval_kge(θprior)
    @printf("prior cal-KGE=%.4f  eval-KGE strict=%.4f bc=%.4f | detuned start cal-KGE=%.4f\n",
        cal.kge_val(cal.series(θprior)), pk[1], pk[2], cal.kge_val(cal.series(θ0)))

    # --- AD: damped 2×2 Newton on the checkpointed cal gradient ---
    ADIT = parse(Int, get(ENV,"CLM_AD_ITER","5")); hs=[1e-3*sc[1],1e-3*sc[2]]
    θ=copy(θ0); npass=0
    for it in 1:ADIT
        g0=cal.grad_ad(θ); npass+=1
        g1=cal.grad_ad(θ.+[hs[1],0.0]); npass+=1
        g2=cal.grad_ad(θ.+[0.0,hs[2]]); npass+=1
        H11=(g1[1]-g0[1])/hs[1]; H21=(g1[2]-g0[2])/hs[1]
        H12=(g2[1]-g0[1])/hs[2]; H22=(g2[2]-g0[2])/hs[2]; H12s=0.5*(H12+H21)
        # LM with DIAGONAL damping (μ·|H_jj|): damps each direction by ITS OWN curvature so the
        # low-curvature (fff) direction is not over-damped by the high-curvature (baseflow) one.
        D1=abs(H11)+1e-30; D2=abs(H22)+1e-30; μ=1e-2; Lc=cal.Loss(θ); acc=false
        for _t in 1:20
            A11=H11+μ*D1; A22=H22+μ*D2; d=A11*A22-H12s^2
            Δ1=-(A22*g0[1]-H12s*g0[2])/d; Δ2=-(-H12s*g0[1]+A11*g0[2])/d
            θt=clamp.(θ.+[Δ1,Δ2],lo,hi); (cal.Loss(θt)<Lc) && (θ=θt; acc=true; break); μ*=3
        end
        @printf("  AD it %2d θ=[%.4f,%.3e] cal-SSE=%.4e cal-KGE=%.4f rev-passes=%d\n",
            it, θ[1], θ[2], cal.Loss(θ), cal.kge_val(cal.series(θ)), npass)
        acc || break
    end
    ad_cal=cal.kge_val(cal.series(θ)); (ad_es, ad_ebc)=eval_kge(θ); ad_theta=copy(θ)
    @printf("AD  calibrated θ=[%.4f,%.3e]  cal-KGE=%.4f  EVAL-KGE strict=%.4f bc=%.4f  reverse-passes=%d\n",
        θ[1], θ[2], ad_cal, ad_es, ad_ebc, npass)

    # --- DDS baseline on the SAME cal season, evaluated on the SAME eval season ---
    rng=MersenneTwister(7); best=copy(θ0); bestL=cal.Loss(best); evals=1
    ITER=parse(Int, get(ENV,"CLM_DDS_ITER","150")); r=0.2
    for it in 1:ITER
        P=1.0-log(it)/log(ITER); cand=copy(best)
        for j in 1:2; rand(rng)<P && (cand[j]=clamp(best[j]+r*(hi[j]-lo[j])*randn(rng),lo[j],hi[j])); end
        L=cal.Loss(cand); evals+=1; L<bestL && (bestL=L; best=cand)
    end
    dds_cal=cal.kge_val(cal.series(best)); (dds_es, dds_ebc)=eval_kge(best)
    @printf("DDS calibrated θ=[%.4f,%.3e]  cal-KGE=%.4f  EVAL-KGE strict=%.4f bc=%.4f  evals=%d\n",
        best[1], best[2], dds_cal, dds_es, dds_ebc, evals)

    println("\n", "="^82)
    @printf("SEASONAL SPLIT-SAMPLE (%s: cal %s, eval %s, %d-day seasons)\n", DOMAIN,
        Dates.format(CAL_START,"yyyy-mm-dd"), Dates.format(EVAL_START,"yyyy-mm-dd"), NDAYS)
    @printf("%-8s | %-8s | %-11s | %-11s | %-11s\n", "which", "cal-KGE", "eval-strict", "eval-bc", "model-evals")
    println("-"^82)
    @printf("%-8s | %-+8.4f | %-+11.4f | %-+11.4f | %d reverse-passes\n", "AD", ad_cal, ad_es, ad_ebc, npass)
    @printf("%-8s | %-+8.4f | %-+11.4f | %-+11.4f | %d forward-evals\n", "DDS", dds_cal, dds_es, dds_ebc, evals)
    println("="^82)
    @printf("headline: held-out eval-KGE (bias-corr, shape gen.) AD=%.4f vs DDS=%.4f | strict AD=%.4f vs DDS=%.4f\n",
        ad_ebc, dds_ebc, ad_es, dds_es)
end

# Fast SCREEN: prior cal/eval KGE + a coarse best-achievable KGE over a small fff×baseflow
# grid (cheap forward evals only) — to pick fittable basin×season combos before a full run.
function screen()
    println("\n", "="^82); println("SCREEN — prior + coarse best KGE (forward-only) for $DOMAIN"); println("="^82)
    cal = make_season(CAL_START; ndays=NDAYS)
    eval = make_season(EVAL_START; ndays=NDAYS, gain_fixed=cal.gain)
    @printf("cal %s (%dd) fin-obs=%d/%d  eval %s fin-obs=%d/%d  gain=%.3e  HRUarea=%.3e m2\n",
        Dates.format(CAL_START,"yyyy-mm-dd"), NDAYS, cal.nfin, NDAYS,
        Dates.format(EVAL_START,"yyyy-mm-dd"), eval.nfin, NDAYS, cal.gain, basin_area_m2())
    @printf("prior cal-KGE=%.4f  eval-KGE=%.4f\n",
        cal.kge_val(cal.series([FFF0,BFS0])), eval.kge_val(eval.series([FFF0,BFS0])))
    bestk=-Inf; bestθ=[FFF0,BFS0]
    for f in (0.01,0.05,0.1,0.3,0.5,1.0,1.5), bf in (1e-4,5e-4,1e-3,5e-3,2e-2)
        k = cal.kge_val(cal.series([f,bf])); (k>bestk) && (bestk=k; bestθ=[f,bf])
    end
    @printf("coarse-grid best cal-KGE=%.4f @θ=[%.3f,%.1e] → its eval-KGE=%.4f\n",
        bestk, bestθ[1], bestθ[2], eval.kge_val(eval.series(bestθ)))
end

# MULTISCREEN: sweep many (start, ndays) windows for this DOMAIN in ONE process. For each,
# report prior cal-KGE + coarse-grid best cal-KGE (forward-only, cheap) to locate windows the
# model can FIT. CLM_STARTS = comma-sep dates; CLM_WINS = comma-sep window lengths (days).
function multiscreen()
    println("\n", "="^90); println("MULTISCREEN — fittable-window survey for $DOMAIN  (prior & coarse-best cal-KGE)"); println("="^90)
    starts = [Date(strip(s)) for s in split(get(ENV,"CLM_STARTS","2004-04-01,2004-05-01,2004-06-01,2004-09-01,2005-04-01,2005-05-01,2005-06-01"), ",")]
    wins = [parse(Int,strip(w)) for w in split(get(ENV,"CLM_WINS","10,20"), ",")]
    @printf("%-12s %-4s | %-8s %-11s %-16s %-8s\n", "start","win","fin-obs","prior-KGE","best-KGE@θ","gain")
    println("-"^90)
    for st in starts, w in wins
        local c
        try; c = make_season(st; ndays=w); catch e; @printf("%-12s %-4d | ERR %s\n", st, w, sprint(showerror,e)[1:min(40,end)]); continue; end
        pk = c.kge_val(c.series([FFF0,BFS0]))
        bestk=-Inf; bestθ=[FFF0,BFS0]
        for f in (0.01,0.1,0.5,1.0), bf in (1e-4,1e-3,5e-3,2e-2)
            k=c.kge_val(c.series([f,bf])); (k>bestk) && (bestk=k; bestθ=[f,bf])
        end
        @printf("%-12s %-4d | %-8d %-+11.4f %-+8.4f@[%.2f,%.0e] %.2e\n",
            Dates.format(st,"yyyy-mm-dd"), w, c.nfin, pk, bestk, bestθ[1], bestθ[2], c.gain)
    end
end

if STAGE == "A"
    ok, _ = stage_a()
    println(ok ? "\n→ Stage A GREEN: proceed to Stage B (CLM_STAGE=B)." : "\n→ Stage A RED: do not run Stage B.")
elseif STAGE == "B"
    stage_b()
elseif STAGE == "SCREEN"
    screen()
elseif STAGE == "MULTISCREEN"
    multiscreen()
end
