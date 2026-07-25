# =============================================================================
# GATE 3 — reverse-AD streamflow calibration ACROSS BASINS vs the DDS baseline.
#
# Runs the windowed reverse-AD surface-runoff calibration (see
# scripts/calibrate_streamflow_ad.jl for the full method + Gate-1/2 detail) on
# several contrasting basins in ONE Julia process (compile cost amortized), and
# emits a results table:  basin | AD-vs-FD rel err | windowed loss ↓ | reverse-passes
# | DDS KGE (cited from SYMFLUENCE optimization/*.json if present).
#
# The differentiated parameters are the real TOPMODEL surface-runoff calibration
# params fff (fover) and fmax (max saturated fraction), injected into the production
# HydrologyNoDrainage reverse thread. HONEST SCOPE: the AD figure is the correctness
# (AD==FD) + windowed-loss reduction of the reverse gradient on each basin's real
# cold-start state under an idealized steady-rain window — NOT a full-period forced
# hydrograph KGE (subsurface-drainage/baseflow physics is a separate discrete routine
# not yet in the reverse tape, and the tractable AD window uses idealized forcing).
# The DDS KGE is the SYMFLUENCE full-period derivative-free baseline, cited for context.
#
#   SYMFLUENCE_DATA=<root> julia --project=. scripts/calibrate_streamflow_ad_multi.jl
# =============================================================================
using CLM, Enzyme, Printf, JSON, NCDatasets, Statistics
const C = CLM
const FT = Float64
const DT = 1800.0
const NSTEPS = parse(Int, get(ENV, "CLM_NSTEPS", "48"))
const RAIN_FLOOR = parse(Float64, get(ENV, "CLM_RAIN", "1e-6"))   # floor below all basin means → real rain distinguishes basins

# Per-basin steady rain input = mean of the basin's REAL precipitation forcing (PRECTmms,
# mm/s), floored so the runoff/gradient signal stays healthy. Gives genuine per-basin
# distinctness tied to real climate data rather than one synthetic rate for all basins.
function basin_rain(domain_dir)
    fdir = joinpath(domain_dir, "data/forcing/CLM_input")
    isdir(fdir) || return RAIN_FLOOR
    ncs = filter(f->endswith(f, ".nc"), readdir(fdir; join=true))
    isempty(ncs) && return RAIN_FLOOR
    try
        ds = NCDataset(first(ncs))
        v = haskey(ds, "PRECTmms") ? ds["PRECTmms"] : (haskey(ds, "PRECT") ? ds["PRECT"] : nothing)
        v === nothing && (close(ds); return RAIN_FLOOR)
        a = Float64.(v[:]); close(ds)
        m = mean(filter(x->isfinite(x) && x < 1e30, a))
        return max(isfinite(m) ? m : RAIN_FLOOR, RAIN_FLOOR)
    catch; return RAIN_FLOOR; end
end
const ROOT   = get(ENV, "SYMFLUENCE_DATA",
    "/Users/darri.eythorsson/Library/CloudStorage/GoogleDrive-dareyt@gmail.com/My Drive/data/SYMFLUENCE_data")

# ---- top-level Enzyme-safe phase functions (per-basin data via Const aux, never closures) ----
function setup_forcing!(a2l, T, ng, rain)
    for g in 1:ng
        a2l.forc_t_not_downscaled_grc[g]=T; a2l.forc_pbot_not_downscaled_grc[g]=85000.0
        a2l.forc_th_not_downscaled_grc[g]=T*(100000.0/85000.0)^(C.RAIR/C.CPAIR)
        a2l.forc_rho_not_downscaled_grc[g]=85000.0/(C.RAIR*T); a2l.forc_lwrad_not_downscaled_grc[g]=300.0
        a2l.forc_vp_grc[g]=800.0; a2l.forc_hgt_grc[g]=30.0; a2l.forc_topo_grc[g]=0.0; a2l.forc_wind_grc[g]=3.0
        a2l.forc_hgt_u_grc[g]=30.0; a2l.forc_hgt_t_grc[g]=30.0; a2l.forc_hgt_q_grc[g]=30.0
        for b in 1:C.NUMRAD; a2l.forc_solad_not_downscaled_grc[g,b]=200.0; a2l.forc_solai_grc[g,b]=80.0; end
        a2l.forc_solar_not_downscaled_grc[g]=560.0
        a2l.forc_rain_not_downscaled_grc[g]=rain; a2l.forc_snow_not_downscaled_grc[g]=0.0
    end
end

function rain_inject!(b, aux)
    wf = b.inst.water.waterfluxbulk_inst.wf
    @inbounds for c in aux.cols; wf.qflx_rain_plus_snomelt_col[c] = aux.rain; end
    return nothing
end
function satexcess_cal!(b, aux)
    i = b.inst; θ1 = b.theta[1]; θ2 = b.theta[2]
    sh = i.soilhydrology; ss = i.soilstate; ser = i.sat_excess_runoff
    wfb = i.water.waterfluxbulk_inst; wf = wfb.wf
    @inbounds for c in aux.cols
        ft = sh.frost_table_col[c]; zw = sh.zwt_col[c]; zwp = sh.zwt_perched_col[c]
        wt = θ2 * ss.wtfact_col[c]
        fs = (ft > zwp && ft <= zw) ? wt*exp(-0.5*θ1*zwp) : wt*exp(-0.5*θ1*zw)
        ser.fsat_col[c] = fs; ser.fcov_col[c] = fs
        wfb.qflx_sat_excess_surf_col[c] = fs * wf.qflx_rain_plus_snomelt_col[c]
    end
    return nothing
end
function accum_runoff!(b, aux)
    wf = b.inst.water.waterfluxbulk_inst.wf; s = 0.0
    @inbounds for c in aux.cols; s += wf.qflx_surf_col[c]; end
    b.series[aux.k] += s
    return nothing
end

# ---- build one basin's reverse-AD calibration context ----
function build_basin(fsurdat, paramfile, rain)
    (inst, bounds, filt, tm) = C.clm_initialize!(; fsurdat=fsurdat, paramfile=paramfile)
    setup_forcing!(inst.atm2lnd, 288.0, bounds.endg, rain)
    C.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
    config = C.CLMDriverConfig(); fia = C.clump_filter_inactive_and_active
    (declin, _) = C.compute_orbital(120.0); nextsw = 120.0 + DT/C.SECSPDAY
    for n in 1:3
        C.clm_drv!(config, inst, filt, fia, bounds, true, nextsw, declin, declin,
            C.ORB_OBLIQR_DEFAULT, false, false, "", false; nstep=n, is_first_step=(n==1),
            is_beg_curr_day=(n==1), dtime=DT, mon=5, day=1, photosyns=inst.photosyns)
    end
    hc = filt.hydrologyc
    hcols = Int[c for c in bounds.begc:bounds.endc if hc[c]]
    fff0 = FT(C.sat_excess_runoff_params.fff)
    a_surf = C.surfhydro_rev_aux(bounds, filt; dtime=DT, use_hydrstress=config.use_hydrstress)
    rax = (; cols = hcols, rain = FT(rain))
    sax = (; cols = hcols)
    step_phases(k) = Any[
        (rain_inject!, (rax,)), (C.setsoilfrac_rev_phase!, (a_surf,)),
        (C.setfloodc_rev_phase!, (a_surf,)), (satexcess_cal!, (sax,)),
        (C.setqflx_rev_phase!, (a_surf,)), (C.inflexcess_rev_phase!, (a_surf,)),
        (C.routeinfl_rev_phase!, (a_surf,)), (C.updateh2osfc_rev_phase!, (a_surf,)),
        (C.infil_rev_phase!, (a_surf,)), (C.totalrunoff_rev_phase!, (a_surf,)),
        (C.soilwater_rev_phase!, (C.soilwater_rev_aux(filt, config; dtime=DT),)),
        (C.watertable_rev_phase!, (C.watertable_rev_aux(bounds, filt, config; dtime=DT),)),
        (C.hydnodrain_rev_phase!, (C.hydnodrain_rev_aux(bounds, filt; dtime=DT),)),
        (accum_runoff!, ((; cols = hcols, k = k),)),
    ]
    steps = [step_phases(k) for k in 1:NSTEPS]
    bundle(θ) = (; C.driver_rev_bundle(deepcopy(inst))..., theta = collect(FT, θ), series = zeros(FT, NSTEPS))
    series(θ) = (b = bundle(θ); for ph in steps, (f,ca) in ph; f(b, ca...); end; copy(b.series))
    function grad_ad(θ, tgt)
        seed!(db, b) = (db.series .= 2 .* (b.series .- tgt); nothing)
        db = C.multistep_reverse!(steps, bundle(θ), seed!); copy(db.theta)
    end
    Loss(θ, tgt) = sum(abs2, series(θ) .- tgt)
    return (; fff0, series, grad_ad, Loss, nhc = length(hcols), rain = FT(rain))
end

# DDS baseline KGE from a SYMFLUENCE optimization best-params JSON (best_score = KGE), if any.
function dds_kge(domain_dir)
    best = nothing
    for (root, _, files) in walkdir(domain_dir)
        for f in files
            endswith(f, "best_params.json") || continue
            try
                j = JSON.parsefile(joinpath(root, f))
                if get(j, "metric", "") == "KGE" && haskey(j, "best_score")
                    s = Float64(j["best_score"])
                    (best === nothing || s > best) && (best = s)
                end
            catch; end
        end
    end
    return best
end

function run_basin(name, domain_dir)
    fsurdat = joinpath(domain_dir, "settings/CLM/parameters/surfdata_clm.nc")
    paramfile = joinpath(domain_dir, "settings/CLM/parameters/clm5_params.nc")
    (isfile(fsurdat) && isfile(paramfile)) || (return (; name, ok=false, msg="missing params"))
    rain = basin_rain(domain_dir)
    ctx = try build_basin(fsurdat, paramfile, rain) catch e
        return (; name, ok=false, msg=sprint(showerror, e)[1:min(80,end)]) end
  try
    θ0 = [ctx.fff0, 1.0]
    tgt = ctx.series(θ0)
    all(isfinite, tgt) || (return (; name, ok=false, msg="non-finite runoff series"))

    # Gate-1 AD-vs-FD at a perturbed θ
    θp = [ctx.fff0*1.3, 1.25]
    g = ctx.grad_ad(θp, tgt)
    relerrs = Float64[]
    for (j,h) in ((1,1e-3*ctx.fff0),(2,1e-3))
        θpp=copy(θp); θpp[j]+=h; θpm=copy(θp); θpm[j]-=h
        fd=(ctx.Loss(θpp,tgt)-ctx.Loss(θpm,tgt))/(2h)
        push!(relerrs, abs(g[j]-fd)/max(abs(fd),1e-30))
    end
    maxrel = maximum(relerrs)

    # Gate-2 well-posed 1-D fff calibration (fmax fixed) toward an fff-perturbed target.
    fmax_fix = 1.0; fff_obs = ctx.fff0*0.6; OBS = ctx.series([fff_obs, fmax_fix])
    fff = ctx.fff0; L0 = ctx.Loss([fff, fmax_fix], OBS); npass = 0
    g1(x) = (npass += 1; ctx.grad_ad([x, fmax_fix], OBS)[1])
    L1(x) = ctx.Loss([x, fmax_fix], OBS); h = 1e-3*ctx.fff0
    for _it in 1:12
        g0=g1(fff); gp=g1(fff+h); H=(gp-g0)/h; μ=1e-6*max(abs(H),1e-12); Lc=L1(fff); acc=false
        for _t in 1:15
            xt=clamp(fff-g0/(H+μ),1e-3,5.0); (L1(xt)<Lc) && (fff=xt; acc=true; break); μ=μ==0 ? 1e-9 : μ*4
        end
        (!acc || abs(g0)<1e-14) && break
    end
    Lf = L1(fff)
    return (; name, ok=true, nhc=ctx.nhc, maxrel, L0, Lf, reduction=L0/max(Lf,1e-300),
              npass, dds=dds_kge(domain_dir), fff_fit=fff, fff_obs=fff_obs, fff0=ctx.fff0, rain=ctx.rain)
  catch e
    return (; name, ok=false, msg="AD: "*sprint(showerror, e)[1:min(76,end)])
  end
end

# ---- basin list (env CLM_BASINS = comma-separated domain names, else default 3) ----
default_basins = ["Bow_at_Banff_lumped", "Stillwater_Oklahoma", "Iceland_Jokulsa_Fjollum",
                  "Alps_Massa_Aletsch_CH", "Boreal_Krycklan_Sweden"]
basins = haskey(ENV, "CLM_BASINS") ? split(ENV["CLM_BASINS"], ",") : default_basins

println("="^96)
println("GATE 3 — reverse-AD streamflow calibration across basins (N=$NSTEPS-step window)")
println("="^96)
results = Any[]
for b in basins
    dom = joinpath(ROOT, "domain_$(strip(b))")
    @printf("\n>>> %s\n", b)
    r = run_basin(strip(b), dom)
    push!(results, r)
    if r.ok
        @printf("    rain=%.2e mm/s (real forcing mean)  hydro-cols=%d  AD-vs-FD max rel=%.2e\n", r.rain, r.nhc, r.maxrel)
        @printf("    windowed loss %.3e→%.3e  reverse-passes=%d  DDS-KGE=%s\n",
            r.L0, r.Lf, r.npass, r.dds===nothing ? "n/a" : @sprintf("%.3f", r.dds))
    else
        @printf("    SKIPPED: %s\n", r.msg)
    end
end

println("\n", "="^96)
@printf("%-26s | %-12s | %-16s | %-8s | %-9s\n", "basin", "AD-vs-FD rel", "windowed loss ↓", "AD-passes", "DDS KGE")
println("-"^96)
for r in results
    if r.ok
        @printf("%-26s | %-12.1e | %-16s | %-8d | %-9s\n", r.name, r.maxrel,
            @sprintf("%.1e→%.1e", r.L0, r.Lf), r.npass, r.dds===nothing ? "n/a" : @sprintf("%.3f", r.dds))
    else
        @printf("%-26s | %-12s | %-16s | %-8s | %-9s\n", r.name, "—", r.msg[1:min(14,end)], "—", "—")
    end
end
println("="^96)
nok = count(r->r.ok && r.maxrel < 1e-4, results)
@printf("%d/%d basins: reverse-AD streamflow-loss gradient CORRECT (AD==FD ≤1e-4) on real basin state\n",
    nok, length(results))
