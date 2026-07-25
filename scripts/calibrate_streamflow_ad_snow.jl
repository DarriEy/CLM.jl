# =============================================================================
# SNOWMELT reverse-AD — Stage A feasibility + gate (task b).
#
# Goal: differentiate a windowed RUNOFF loss w.r.t. a real SNOW parameter, via the melt
# physics already inside soil_temperature! (the snow+soil phase change computes
# qflx_snomelt = smooth_max(0, wice0 − h2osoi_ice)/dtime — DIFFERENTIABLE), then connect
# qflx_snomelt → qflx_rain_plus_snomelt → the existing surface-runoff reverse thread.
#
# SNOW-PHASE CLASSIFICATION (src/driver/clm_driver.jl order):
#   1789 soil_temperature!      — snow+soil energy balance + PHASE CHANGE → qflx_snomelt
#                                 (smooth_max, soil_temperature.jl:2118)      [DIFFERENTIABLE]
#   1837 SnowWater (percolation)— bulk_flux_snow_percolation!/update_state    [differentiable]
#        → qflx_rain_plus_snomelt = qflx_liq_grnd + qflx_snomelt (snow_hydrology.jl:1125)
#   2143 snow_compaction!       — continuous densification                    [differentiable]
#   2159 combine_snow_layers!   — REMOVE thin layers (snl changes)      [DISCRETE-STRUCTURAL]
#   2170 divide_snow_layers!    — SPLIT thick layers (snl changes)      [DISCRETE-STRUCTURAL]
#   2191 hydrology_no_drainage! — diagnostics
# The discrete layer ops (combine/divide) run AFTER the melt→runoff path and change the
# state DIMENSION (snl), which reverse-AD cannot tape. FALLBACK (this script): a FIXED-LAYER
# window — differentiate melt→runoff, treat combine/divide as non-differentiated between-step
# transitions (omit them; valid when the layer count is steady in a mid-melt window).
#
# Injected snow param θ_snow: a multiplier on the absorbed shortwave in the snow layers
# (inst.solarabs.sabg_lyr_patch / sabg_snow_patch) — physically the snow shortwave
# absorption ≈ (1 − snow albedo), the dominant spring-melt energy driver.
#
#   CLM_SNOW_STAGE=probe julia --project=. scripts/calibrate_streamflow_ad_snow.jl
# =============================================================================
using CLM, Enzyme, Printf, Statistics, NCDatasets, Dates, Random, LinearAlgebra
const C = CLM
const FT = Float64
const DT = 1800.0
const STEPS_PER_DAY = 48
const FSURDAT = get(ENV, "CLM_FSURDAT",
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc")
const PARAMFILE = get(ENV, "CLM_PARAMFILE",
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc")
const STAGE = get(ENV, "CLM_SNOW_STAGE", "probe")

# cold/snow forcing to accumulate a pack; then a warmer high-solar melt forcing.
function set_forcing!(a2l, ng; T, solad, rain, snow)
    for g in 1:ng
        a2l.forc_t_not_downscaled_grc[g]=T; a2l.forc_pbot_not_downscaled_grc[g]=85000.0
        a2l.forc_th_not_downscaled_grc[g]=T*(100000.0/85000.0)^(C.RAIR/C.CPAIR)
        a2l.forc_rho_not_downscaled_grc[g]=85000.0/(C.RAIR*T); a2l.forc_lwrad_not_downscaled_grc[g]=250.0+ (T-260)*3
        a2l.forc_vp_grc[g]=300.0; a2l.forc_hgt_grc[g]=30.0; a2l.forc_topo_grc[g]=0.0; a2l.forc_wind_grc[g]=2.0
        a2l.forc_hgt_u_grc[g]=30.0; a2l.forc_hgt_t_grc[g]=30.0; a2l.forc_hgt_q_grc[g]=30.0
        for b in 1:C.NUMRAD; a2l.forc_solad_not_downscaled_grc[g,b]=solad; a2l.forc_solai_grc[g,b]=solad*0.4; end
        a2l.forc_solar_not_downscaled_grc[g]=solad*2.8
        a2l.forc_rain_not_downscaled_grc[g]=rain; a2l.forc_snow_not_downscaled_grc[g]=snow
    end
end

function build_snowpack(; nacc=48, nonset=2)
    (inst, bounds, filt, tm) = C.clm_initialize!(; fsurdat=FSURDAT, paramfile=PARAMFILE)
    config = C.CLMDriverConfig(); fia = C.clump_filter_inactive_and_active
    (declin,_) = C.compute_orbital(80.0); nextsw = 80.0 + DT/C.SECSPDAY
    runstep!(n; first, T, solad, rain, snow) = begin
        set_forcing!(inst.atm2lnd, bounds.endg; T, solad, rain, snow)
        C.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
        C.clm_drv!(config, inst, filt, fia, bounds, true, nextsw, declin, declin,
            C.ORB_OBLIQR_DEFAULT, false, false, "", false; nstep=n, is_first_step=first,
            is_beg_curr_day=first, dtime=DT, mon=3, day=20, photosyns=inst.photosyns)
    end
    # accumulate a NEAR-FREEZING snowpack (below 0°C so it stays frozen, but warm enough that a
    # melt-energy injection in the reverse window pushes the top layers over freezing). Low solar
    # during accumulation so the pack survives; melt is driven IN the window via sabg injection.
    for n in 1:nacc; runstep!(n; first=(n==1), T=270.5, solad=12.0, rain=0.0, snow=4.0e-4); end
    return inst, bounds, filt, config
end

inst, bounds, filt, config = build_snowpack()
hc = filt.hydrologyc
const HCOLS = Int[c for c in bounds.begc:bounds.endc if hc[c]]
c0 = HCOLS[1]
nsno = C.varpar.nlevsno
# shallow water table (1 m head above bedrock) so the baseflow lever is non-vacuous (flagged;
# same device as the seasonal script).
let jz=nsno+1
    for c in HCOLS; nb=inst.column.nbedrock[c]; inst.soilhydrology.zwt_col[c]=max(0.3, inst.column.zi[c,nb+jz]-1.0); end
end
snl = inst.column.snl[c0]
@printf("snowpack built: col %d  snl=%d (%d snow layers)  h2osno=%.3f mm\n",
    c0, snl, -snl, sum(inst.water.waterstatebulk_inst.ws.h2osoi_ice_col[c0, max(1,nsno+snl+1):nsno]))
@printf("t_soisno[snow layers] = %s\n",
    string(round.(inst.temperature.t_soisno_col[c0, max(1,nsno+snl+1):nsno], digits=2)))
let sl = inst.solarabs.sabg_lyr_patch
    fin = [v for v in sl if isfinite(v) && v != C.SPVAL]
    @printf("sabg_lyr_patch finite entries=%d  max=%.2f W/m2 | sabg_snow max=%.2f\n",
        length(fin), isempty(fin) ? 0.0 : maximum(fin), maximum(v for v in inst.solarabs.sabg_snow_patch if isfinite(v); init=0.0))
end

# ---- probe: run soil_temperature! once forward on the melt state; is qflx_snomelt>0? ----
# Inject a melt-energy field: set the active snow layers' absorbed shortwave to θ·BASE (the
# snow shortwave-absorption ≈ (1−albedo) melt-energy param). Only active snow layers (jtop..0),
# finite & non-SPVAL, are touched. This drives melt inside soil_temperature! in the window.
const SABG_BASE = 260.0   # W/m2 baseline snow-absorbed shortwave (θ_snow=1)
function inject_melt_energy!(s, θ)
    sl = s.solarabs.sabg_lyr_patch; ss = s.solarabs.sabg_snow_patch
    @inbounds for p in axes(sl,1)
        any_active = false
        for j in axes(sl,2)
            v = sl[p,j]
            if isfinite(v) && v != C.SPVAL
                sl[p,j] = θ*SABG_BASE; any_active = true
            end
        end
        (any_active && isfinite(ss[p])) && (ss[p] = θ*SABG_BASE)
    end
end

if STAGE == "probe"
    println("\n[PROBE] set sabg_lyr=θ·$(SABG_BASE) on snow layers, run soil_temperature! → qflx_snomelt")
    aux = C.soiltemp_rev_aux(bounds, filt; dtime=DT)
    function melt_at(θ)
        s = deepcopy(inst); inject_melt_energy!(s, θ)
        b = (; inst=s, scratch=C.cf_rev_scratch(FT, length(s.patch.column)))
        C.soiltemp_rev_phase!(b, aux)
        return s.water.waterfluxbulk_inst.wf.qflx_snomelt_col[c0]
    end
    q1 = melt_at(1.0); q15 = melt_at(1.5); q05 = melt_at(0.5)
    @printf("  qflx_snomelt: θ=0.5→%.6e  θ=1.0→%.6e  θ=1.5→%.6e\n", q05, q1, q15)
    @printf("  melt %s θ_snow (monotone-in-energy: %s)\n",
        (q15 != q1 || q05 != q1) ? "RESPONDS to" : "INSENSITIVE to", string(q05 <= q1 <= q15))
end

# =============================================================================
# STAGE A GATE — snow-inclusive reverse gradient of a windowed runoff loss vs central FD.
# θ = [θ_snow (melt energy), θ_fff (surface saturation), θ_baseflow]. Layer structure FIXED
# (combine/divide omitted — the blessed fallback). Melt → qflx_rain_plus_snomelt → surface runoff.
# =============================================================================
# active-snow-layer mask (fixed within the window — no combine/divide fires)
const SMASK = let sl=inst.solarabs.sabg_lyr_patch
    [isfinite(sl[p,j]) && sl[p,j] != C.SPVAL for p in axes(sl,1), j in axes(sl,2)]
end
function meltenergy_phase!(b, aux)
    θ = b.theta[1]; sl = b.inst.solarabs.sabg_lyr_patch; ss = b.inst.solarabs.sabg_snow_patch
    @inbounds for p in axes(sl,1)
        active=false
        for j in axes(sl,2)
            if aux.mask[p,j]; sl[p,j] = θ*aux.base; active=true; end
        end
        active && (ss[p] = θ*aux.base)
    end
    return nothing
end
function snomelt_connect!(b, aux)
    wf = b.inst.water.waterfluxbulk_inst.wf
    @inbounds for c in aux.cols; wf.qflx_rain_plus_snomelt_col[c] = wf.qflx_snomelt_col[c]; end
    return nothing
end
function satexcess_cal!(b, aux)
    i=b.inst; θ1=b.theta[2]; sh=i.soilhydrology; ss=i.soilstate; ser=i.sat_excess_runoff
    wfb=i.water.waterfluxbulk_inst; wf=wfb.wf
    @inbounds for c in aux.cols
        ft=sh.frost_table_col[c]; zw=sh.zwt_col[c]; zwp=sh.zwt_perched_col[c]; wt=ss.wtfact_col[c]
        fs=(ft>zwp && ft<=zw) ? wt*exp(-0.5*θ1*zwp) : wt*exp(-0.5*θ1*zw)
        ser.fsat_col[c]=fs; ser.fcov_col[c]=fs
        wfb.qflx_sat_excess_surf_col[c]=fs*wf.qflx_rain_plus_snomelt_col[c]
    end
    return nothing
end
function baseflow_cal!(b, aux)
    i=b.inst; θ2=b.theta[3]; sh=i.soilhydrology; col=i.column
    ws=i.water.waterstatebulk_inst.ws; wf=i.water.waterfluxbulk_inst.wf; jz=aux.nlevsno+1
    @inbounds for c in aux.cols
        nb=col.nbedrock[c]; zi_bed=col.zi[c,nb+jz]; zw=sh.zwt_col[c]
        q = zw<=zi_bed ? θ2*tan(FT(pi)/180.0*col.topo_slope[c])*(zi_bed-zw)^aux.n_baseflow : zero(zw)
        wf.qflx_drain_col[c]=q
        jl=min(nb,aux.nlevsoi); nl=ws.h2osoi_liq_col[c,jl+aux.nlevsno]-q*aux.dt
        ws.h2osoi_liq_col[c,jl+aux.nlevsno] = nl>0 ? nl : ws.h2osoi_liq_col[c,jl+aux.nlevsno]
    end
    return nothing
end
function accum_runoff!(b, aux)
    wf=b.inst.water.waterfluxbulk_inst.wf; s=zero(eltype(b.series))
    @inbounds for c in aux.cols; s += wf.qflx_surf_col[c]+wf.qflx_drain_col[c]; end
    b.series[aux.k]+=s; return nothing
end

if STAGE == "gate"
    NSTEPS = parse(Int, get(ENV,"CLM_NSTEPS","3"))
    a_surf = C.surfhydro_rev_aux(bounds, filt; dtime=DT, use_hydrstress=config.use_hydrstress)
    melt_aux = (; mask=SMASK, base=SABG_BASE)
    st_aux = C.soiltemp_rev_aux(bounds, filt; dtime=DT)
    sc_aux = (; cols=HCOLS)
    bf_aux = (; cols=HCOLS, nlevsno=nsno, nlevsoi=C.varpar.nlevsoi,
                n_baseflow=FT(C.soilhydrology_params.n_baseflow), dt=DT)
    step_phases(k) = Any[
        (meltenergy_phase!, (melt_aux,)),
        (C.soiltemp_rev_phase!, (st_aux,)),
        (snomelt_connect!, (sc_aux,)),
        (C.setsoilfrac_rev_phase!, (a_surf,)), (C.setfloodc_rev_phase!, (a_surf,)),
        (satexcess_cal!, (sc_aux,)),
        (C.setqflx_rev_phase!, (a_surf,)), (C.inflexcess_rev_phase!, (a_surf,)),
        (C.routeinfl_rev_phase!, (a_surf,)), (C.updateh2osfc_rev_phase!, (a_surf,)),
        (C.infil_rev_phase!, (a_surf,)), (C.totalrunoff_rev_phase!, (a_surf,)),
        (C.soilwater_rev_phase!, (C.soilwater_rev_aux(filt, config; dtime=DT),)),
        (C.watertable_rev_phase!, (C.watertable_rev_aux(bounds, filt, config; dtime=DT),)),
        (C.hydnodrain_rev_phase!, (C.hydnodrain_rev_aux(bounds, filt; dtime=DT),)),
        (baseflow_cal!, (bf_aux,)),
        (accum_runoff!, ((; cols=HCOLS, k=k),))]
    STEPS = [step_phases(k) for k in 1:NSTEPS]
    bundle(θ) = (; C.driver_rev_bundle(deepcopy(inst))..., theta=collect(FT,θ), series=zeros(FT,NSTEPS))
    series(θ) = (b=bundle(θ); for ph in STEPS, (f,ca) in ph; f(b,ca...); end; copy(b.series))
    Loss(θ, tgt) = sum(abs2, series(θ) .- tgt)
    function grad_ad(θ, tgt)
        seed!(db,b)=(db.series .= 2 .*(b.series .- tgt); nothing)
        db = C.multistep_reverse!(STEPS, bundle(θ), seed!); copy(db.theta)
    end

    println("\n", "="^80); println("STAGE A GATE — snow-melt reverse-AD gradient vs central FD"); println("="^80)
    θ0 = [1.0, FT(C.sat_excess_runoff_params.fff), FT(C.BASEFLOW_SCALAR[])]
    tgt = series(θ0)
    @printf("window N=%d  runoff series (mm/s)=%s\n", NSTEPS, string(round.(tgt, sigdigits=4)))
    θp = [1.2, θ0[2]*1.15, θ0[3]*1.3]
    t=@elapsed g = grad_ad(θp, tgt)
    @printf("AD dL/dθ @[snow=%.2f,fff=%.3f,bf=%.2e] = [% .5e, % .5e, % .5e] (%.1fs)\n",
        θp[1],θp[2],θp[3], g[1],g[2],g[3], t)
    names=["θ_snow (melt energy)","θ_fff (saturation)","θ_baseflow"]
    hs=[1e-3, 1e-3*θ0[2], 1e-3*θ0[3]]; passes=Bool[]
    for j in 1:3
        θa=copy(θp); θa[j]+=hs[j]; θb=copy(θp); θb[j]-=hs[j]
        fd=(Loss(θa,tgt)-Loss(θb,tgt))/(2hs[j]); rel=abs(g[j]-fd)/max(abs(fd),1e-30)
        vac=abs(g[j])<1e-30 && abs(fd)<1e-30; pass=(rel<1e-5)&&!vac; push!(passes,pass)
        @printf("  %-22s AD=% .5e FD=% .5e rel=%.2e %s\n", names[j], g[j], fd, rel,
            vac ? "VACUOUS" : (pass ? "PASS ✓" : "FAIL ✗"))
    end
    println("-"^80)
    ok = all(passes)
    @printf("STAGE A GATE %s — snow-melt-driven runoff gradient %s\n",
        ok ? "PASS ✓" : "FAIL ✗", ok ? "matches FD on θ_snow, θ_fff, θ_baseflow" : "see failures")
end

# =============================================================================
# STAGE B — REAL melt-season split-sample: melt energy driven by REAL solar forcing
# (sabg_lyr = θ_snow · FSDS[k] · ABSORP), spring window, {θ_snow, fff, baseflow_scalar}
# calibrated on the cal season, evaluated on the independent same-season-next-year window,
# AD (LM diagonal-damped, checkpointed) vs DDS. Fixed-layer window (combine/divide omitted).
# =============================================================================
const ROOT = get(ENV, "SYMFLUENCE_DATA",
    "/Users/darri.eythorsson/Library/CloudStorage/GoogleDrive-dareyt@gmail.com/My Drive/data/SYMFLUENCE_data")
const DOMAIN = get(ENV, "CLM_DOMAIN", "Boreal_Krycklan_Sweden")
const DOM_DIR = joinpath(ROOT, "domain_$DOMAIN")
const SNOW_FSURDAT = joinpath(DOM_DIR, "settings/CLM/parameters/surfdata_clm.nc")
const SNOW_PARAMFILE = joinpath(DOM_DIR, "settings/CLM/parameters/clm5_params.nc")
const ABSORP = parse(Float64, get(ENV, "CLM_ABSORP", "0.12"))   # net snow SW absorption→melt factor
const NACC = parse(Int, get(ENV, "CLM_NACC", "140"))            # accumulation steps (pack thickness)
const NDAYS = parse(Int, get(ENV, "CLM_NDAYS", "10"))
const ROUTE_K = parse(Float64, get(ENV, "CLM_ROUTEK", "3.0"))

read_fsds_window(fsurdir, start, nsteps) = begin
    yr=year(start); ff=joinpath(DOM_DIR,"data/forcing/CLM_input","clmforc.$(yr).nc")
    ds=NCDataset(ff); traw=ds["time"][:]
    times=DateTime[ismissing(t) ? DateTime(0) : DateTime(t) for t in traw]
    fs=vec(Float64.(ds["FSDS"][:])); close(ds)
    i0=findfirst(>=(DateTime(start)), times); i0===nothing && error("start not in forcing")
    out=zeros(FT,nsteps)
    for k in 1:nsteps; h=i0+(k-1)÷2; h>length(fs) && break; v=fs[h]; out[k]=FT((isfinite(v)&&v<1e30) ? max(v,0.0) : 0.0); end
    out
end
read_obs_daily_B(start, ndays) = begin
    f=joinpath(DOM_DIR,"data/observations/streamflow/preprocessed","$(DOMAIN)_streamflow_processed.csv")
    obs=fill(NaN,ndays); ssum=Dict{Date,Float64}(); scnt=Dict{Date,Int}()
    for line in readlines(f)[2:end]
        p=split(line,','); length(p)<2 && continue
        dt=tryparse(Date,strip(p[1])[1:min(10,end)]); dt===nothing && continue
        q=tryparse(Float64,strip(p[2])); (q===nothing||!isfinite(q)) && continue
        ssum[dt]=get(ssum,dt,0.0)+q; scnt[dt]=get(scnt,dt,0)+1
    end
    for d in 1:ndays; k=start+Day(d-1); haskey(scnt,k) && (obs[d]=ssum[k]/scnt[k]); end
    obs
end

function make_snow_season(start::Date; ndays::Int, gain_fixed=nothing)
    nsteps=ndays*STEPS_PER_DAY
    (I,B,F,cfg)=let
        (inst,bounds,filt,tm)=C.clm_initialize!(; fsurdat=SNOW_FSURDAT, paramfile=SNOW_PARAMFILE)
        config=C.CLMDriverConfig(); fia=C.clump_filter_inactive_and_active
        (declin,_)=C.compute_orbital(Float64(dayofyear(start))); nextsw=dayofyear(start)+DT/C.SECSPDAY
        runstep!(n;first,T,solad,snow)=begin
            set_forcing!(inst.atm2lnd,bounds.endg;T=T,solad=solad,rain=0.0,snow=snow)
            C.downscale_forcings!(bounds,inst.atm2lnd,inst.column,inst.landunit,inst.topo)
            C.clm_drv!(config,inst,filt,fia,bounds,true,nextsw,declin,declin,C.ORB_OBLIQR_DEFAULT,
                false,false,"",false; nstep=n,is_first_step=first,is_beg_curr_day=first,dtime=DT,
                mon=month(start),day=max(1,day(start)-5),photosyns=inst.photosyns)
        end
        for n in 1:NACC; runstep!(n;first=(n==1),T=270.8,solad=10.0,snow=8.0e-4); end  # near-freezing pack
        (inst,bounds,filt,config)
    end
    hcols=Int[c for c in B.begc:B.endc if F.hydrologyc[c]]
    let jz=C.varpar.nlevsno+1
        for c in hcols; nb=I.column.nbedrock[c]; I.soilhydrology.zwt_col[c]=max(0.3,I.column.zi[c,nb+jz]-1.0); end
    end
    smask=let sl=I.solarabs.sabg_lyr_patch
        [isfinite(sl[p,j]) && sl[p,j]!=C.SPVAL for p in axes(sl,1), j in axes(sl,2)]
    end
    fsds=read_fsds_window(SNOW_FSURDAT, start, nsteps); obs=read_obs_daily_B(start,ndays)
    a_surf=C.surfhydro_rev_aux(B,F; dtime=DT, use_hydrstress=cfg.use_hydrstress)
    st_aux=C.soiltemp_rev_aux(B,F; dtime=DT); sc_aux=(; cols=hcols)
    bf_aux=(; cols=hcols, nlevsno=C.varpar.nlevsno, nlevsoi=C.varpar.nlevsoi,
             n_baseflow=FT(C.soilhydrology_params.n_baseflow), dt=DT)
    step_phases(k)=Any[
        (meltenergy_phase!, ((; mask=smask, base=fsds[k]*ABSORP),)),
        (C.soiltemp_rev_phase!, (st_aux,)), (snomelt_connect!, (sc_aux,)),
        (C.setsoilfrac_rev_phase!,(a_surf,)),(C.setfloodc_rev_phase!,(a_surf,)),
        (satexcess_cal!,(sc_aux,)),
        (C.setqflx_rev_phase!,(a_surf,)),(C.inflexcess_rev_phase!,(a_surf,)),
        (C.routeinfl_rev_phase!,(a_surf,)),(C.updateh2osfc_rev_phase!,(a_surf,)),
        (C.infil_rev_phase!,(a_surf,)),(C.totalrunoff_rev_phase!,(a_surf,)),
        (C.soilwater_rev_phase!,(C.soilwater_rev_aux(F,cfg;dtime=DT),)),
        (C.watertable_rev_phase!,(C.watertable_rev_aux(B,F,cfg;dtime=DT),)),
        (C.hydnodrain_rev_phase!,(C.hydnodrain_rev_aux(B,F;dtime=DT),)),
        (baseflow_cal!,(bf_aux,)),(accum_runoff!,((;cols=hcols,k=k),))]
    steps=[step_phases(k) for k in 1:nsteps]
    bundle(θ)=(; C.driver_rev_bundle(deepcopy(I))..., theta=collect(FT,θ), series=zeros(FT,nsteps))
    series(θ)=(b=bundle(θ); for ph in steps,(f,ca) in ph; f(b,ca...); end; copy(b.series))
    # snl + min-ice diagnostic over the window (fixed-layer check)
    function layer_diag(θ)
        b=bundle(θ); snls=Int[]; minice=Inf
        for ph in steps; for (f,ca) in ph; f(b,ca...); end
            c=hcols[1]; push!(snls,b.inst.column.snl[c])
            ice=sum(b.inst.water.waterstatebulk_inst.ws.h2osoi_ice_col[c, max(1,C.varpar.nlevsno+b.inst.column.snl[c]+1):C.varpar.nlevsno])
            minice=min(minice,ice)
        end
        (snl_start=snls[1], snl_end=snls[end], snl_const=all(==(snls[1]),snls), minice=minice)
    end
    daily_mean(s)=[mean(@view s[(d-1)*STEPS_PER_DAY+1:d*STEPS_PER_DAY]) for d in 1:ndays]
    route(q,k)= k<=1.0 ? q : (α=1/k; o=similar(q); o[1]=q[1]; (for i in 2:length(q); o[i]=(1-α)*o[i-1]+α*q[i]; end); o)
    fin=isfinite.(obs)
    gain = gain_fixed!==nothing ? gain_fixed : begin
        m=mean(daily_mean(series([1.0,FT(C.sat_excess_runoff_params.fff),FT(C.BASEFLOW_SCALAR[])])))
        om=mean(o for o in obs if isfinite(o)); (isfinite(m)&&m>0&&isfinite(om)) ? om/m : 1.0
    end
    sse(s)=sum((route(gain.*daily_mean(s),ROUTE_K)[fin].-obs[fin]).^2)
    function kge_g(s,g)
        qs=route(g.*daily_mean(s),ROUTE_K)[fin]; o=obs[fin]
        (length(o)<3||std(o)<1e-9||std(qs)<1e-12) && return -Inf
        r=cor(qs,o); 1-sqrt((r-1)^2+(std(qs)/std(o)-1)^2+(mean(qs)/mean(o)-1)^2)
    end
    own_gain(s)=(m=mean(daily_mean(s)); om=mean(o for o in obs if isfinite(o)); (isfinite(m)&&m>0) ? om/m : gain)
    # exact linear obs-operator adjoint (as in the seasonal script)
    function dL(s)
        qs=route(gain.*daily_mean(s),ROUTE_K); dqs=zeros(ndays)
        for d in 1:ndays; fin[d] && (dqs[d]=2*(qs[d]-obs[d])); end
        α=ROUTE_K<=1 ? 1.0 : 1/ROUTE_K; dqg=zeros(ndays)
        if ROUTE_K<=1; dqg.=dqs else adj=copy(dqs); for i in ndays:-1:2; dqg[i]+=α*adj[i]; adj[i-1]+=(1-α)*adj[i]; end; dqg[1]+=adj[1] end
        ds=zeros(nsteps)
        for d in 1:ndays; g=gain*dqg[d]/STEPS_PER_DAY; for k in (d-1)*STEPS_PER_DAY+1:d*STEPS_PER_DAY; ds[k]=g; end; end
        ds
    end
    grad_ad(θ)=(s=series(θ); dd=dL(s); seed!(db,b)=(db.series.=dd;nothing); copy(C.multistep_reverse!(steps,bundle(θ),seed!).theta))
    Loss(θ)=sse(series(θ))
    return (; start,ndays,series,grad_ad,Loss,kge_g,own_gain,gain,layer_diag,nfin=count(fin))
end

if STAGE == "seasonB"
    CAL=Date(get(ENV,"CLM_CAL_START","2010-04-20")); EV=Date(get(ENV,"CLM_EVAL_START","2011-04-20"))
    fff0=FT(C.sat_excess_runoff_params.fff); bf0=FT(C.BASEFLOW_SCALAR[])
    println("\n", "="^84); @printf("STAGE B — REAL melt-season split-sample (%s)\n", DOMAIN); println("="^84)
    @printf("cal=%s  eval=%s  NDAYS=%d  ABSORP=%.2f  pack-accum=%d steps\n", CAL, EV, NDAYS, ABSORP, NACC)
    cal=make_snow_season(CAL; ndays=NDAYS)
    ld=cal.layer_diag([1.0,fff0,bf0])
    @printf("cal LAYER DIAG: snl %d→%d  constant=%s  min-ice=%.2f mm (fixed-layer %s)\n",
        ld.snl_start, ld.snl_end, ld.snl_const, ld.minice, (ld.snl_const && ld.minice>0.1) ? "VALID ✓" : "SHRINK NEEDED")
    tgt=cal.series([1.0,fff0,bf0])
    @printf("cal runoff series (mm/s, daily-ish): finite-obs=%d/%d  prior cal-KGE=%.4f\n",
        cal.nfin, NDAYS, cal.kge_g(cal.series([1.0,fff0,bf0]), cal.gain))

    # HARD GATE on the melt-season window
    θp=[1.2, fff0*1.15, bf0*1.3]; g=cal.grad_ad(θp)
    println("\n[GATE] melt-season windowed gradient vs central FD")
    hs=[1e-3,1e-3*fff0,1e-3*bf0]; names=["θ_snow","fff","baseflow"]; gok=true
    gp=Bool[]
    cfd(j,h)=(θa=copy(θp);θa[j]+=h;θb=copy(θp);θb[j]-=h;(cal.Loss(θa)-cal.Loss(θb))/(2h))
    for j in 1:3
        f1=cfd(j,hs[j]); f2=cfd(j,hs[j]/2); fd=(4*f2-f1)/3   # Richardson (O(h^4)) central FD
        rel=abs(g[j]-fd)/max(abs(fd),1e-30)
        vac=abs(g[j])<1e-30&&abs(fd)<1e-30; pass=(rel<1e-5)&&!vac; push!(gp,pass)
        @printf("  %-10s AD=% .6e FD(Rich)=% .6e rel=%.2e %s\n", names[j], g[j], fd, rel, vac ? "VACUOUS" : (pass ? "PASS" : "FAIL"))
    end
    gok=all(gp)
    @printf("[GATE] %s\n", gok ? "PASS ✓" : "FAIL ✗")

    if gok && get(ENV,"CLM_GATE_ONLY","0")=="0"
        evc=make_snow_season(EV; ndays=NDAYS, gain_fixed=cal.gain)
        lo=[1e-3,1e-3,1e-5]; hi=[5.0,2.0,0.05]; sc=[1.0,fff0,bf0]
        θ0=clamp.([1.0*1.4, fff0*1.4, bf0*1.4^2], lo, hi)
        eval_kge(θ)=(s=evc.series(θ); (evc.kge_g(s,cal.gain), evc.kge_g(s,evc.own_gain(s))))
        @printf("\nprior cal-KGE=%.4f  eval strict/bc=%.4f/%.4f | detuned start cal-KGE=%.4f\n",
            cal.kge_g(cal.series([1.0,fff0,bf0]),cal.gain), eval_kge([1.0,fff0,bf0])..., cal.kge_g(cal.series(θ0),cal.gain))
        # AD: LM diagonal-damped Newton  (wrapped in a let so loop assignments resolve)
        let
        θ=copy(θ0); npass=0; hstep=1e-3.*sc
        for it in 1:parse(Int,get(ENV,"CLM_AD_ITER","6"))
            g0=cal.grad_ad(θ); npass+=1
            gg=[cal.grad_ad(θ.+((1:3).==j).*hstep[j]) for j in 1:3]; npass+=3
            H=[ (gg[j][i]-g0[i])/hstep[j] for i in 1:3, j in 1:3 ]
            Hs=0.5.*(H+H'); D=abs.([Hs[i,i] for i in 1:3]).+1e-30; μ=1e-2; Lc=cal.Loss(θ); acc=false
            for _t in 1:20
                A=Hs+μ.*Diagonal(D); Δ=-(A\g0)
                θt=clamp.(θ.+Δ,lo,hi); (cal.Loss(θt)<Lc) && (θ=θt;acc=true;break); μ*=3
            end
            @printf("  AD it %2d θ=[snow=%.3f,fff=%.3f,bf=%.2e] cal-KGE=%.4f passes=%d\n",
                it,θ[1],θ[2],θ[3],cal.kge_g(cal.series(θ),cal.gain),npass)
            acc || break
        end
        ad_cal=cal.kge_g(cal.series(θ),cal.gain); (ad_es,ad_ebc)=eval_kge(θ)
        @printf("AD  θ=[%.3f,%.3f,%.2e] cal-KGE=%.4f EVAL strict/bc=%.4f/%.4f passes=%d\n",
            θ[1],θ[2],θ[3],ad_cal,ad_es,ad_ebc,npass)
        # DDS
        rng=MersenneTwister(11); best=copy(θ0); bestL=cal.Loss(best); evals=1
        for it in 1:parse(Int,get(ENV,"CLM_DDS_ITER","90"))
            P=1-log(it)/log(90); cand=copy(best)
            for j in 1:3; rand(rng)<P && (cand[j]=clamp(best[j]+0.2*(hi[j]-lo[j])*randn(rng),lo[j],hi[j])); end
            L=cal.Loss(cand); evals+=1; L<bestL && (bestL=L;best=cand)
        end
        dds_cal=cal.kge_g(cal.series(best),cal.gain); (dds_es,dds_ebc)=eval_kge(best)
        @printf("DDS θ=[%.3f,%.3f,%.2e] cal-KGE=%.4f EVAL strict/bc=%.4f/%.4f evals=%d\n",
            best[1],best[2],best[3],dds_cal,dds_es,dds_ebc,evals)
        println("\n", "="^84)
        @printf("MELT-SEASON SPLIT-SAMPLE (%s: cal %s → eval %s, %d-day)\n", DOMAIN,
            Dates.format(CAL,"yyyy-mm-dd"), Dates.format(EV,"yyyy-mm-dd"), NDAYS)
        @printf("%-6s | %-8s | %-11s | %-11s | %-10s\n","which","cal-KGE","eval-strict","eval-bc","evals")
        println("-"^84)
        @printf("%-6s | %-+8.4f | %-+11.4f | %-+11.4f | %d rev-passes\n","AD",ad_cal,ad_es,ad_ebc,npass)
        @printf("%-6s | %-+8.4f | %-+11.4f | %-+11.4f | %d fwd-evals\n","DDS",dds_cal,dds_es,dds_ebc,evals)
        println("="^84)
        end  # let
    end
end
