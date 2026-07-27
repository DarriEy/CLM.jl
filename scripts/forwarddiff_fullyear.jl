# =============================================================================
# FORWARD-MODE (ForwardDiff) full-year-scale streamflow gradient + snow split-sample.
#
# Builds on scripts/forwarddiff_feasibility.jl (which proved ForwardDiff propagates through
# the hydrology + soil_temperature! snow-melt physics, Dual-clean, ZERO src changes). Here:
#   STAGE combine  — the one untested spot: do combine_snow_layers!/divide_snow_layers!
#                    (discrete snow-layer merge/split, snl changes) run with DUAL state?
#   STAGE gate     — ForwardDiff gradient of a windowed melt-season streamflow loss vs central
#                    FD at increasing length (up to a year), + wall-time/gradient.
#   STAGE split    — ForwardDiff LM calibration of {θ_snow, fff, baseflow_scalar} on a melt
#                    season, evaluated on the independent same-season-next-year window, vs DDS.
#
#   STAGE sweep    — ADVERSARIAL: sweep θ_snow across combine/divide (snl) + imelt transitions,
#                    FwdDiff gradient vs central FD at each θ. Shows the gradient is only
#                    piecewise-valid: exact away from transitions, ~1e-3 imelt kink on melt
#                    steps, 11-85% at snl straddles → the SMOOTHNESS boundary (see
#                    docs/FORWARD_AD_SMOOTHNESS_BOUNDARY.md). Run with CLM_LAYERS=1.
#
#   CLM_FY_STAGE=combine|gate|split|sweep  julia --project=. scripts/forwarddiff_fullyear.jl
# =============================================================================
using CLM, ForwardDiff, Printf, Statistics, NCDatasets, Dates, Random, LinearAlgebra
const C = CLM
const FT = Float64
const DT = 1800.0
const STEPS_PER_DAY = 48
const ROOT = get(ENV, "SYMFLUENCE_DATA",
    "/Users/darri.eythorsson/Library/CloudStorage/GoogleDrive-dareyt@gmail.com/My Drive/data/SYMFLUENCE_data")
const DOMAIN = get(ENV, "CLM_DOMAIN", "Boreal_Krycklan_Sweden")
const DOM_DIR = joinpath(ROOT, "domain_$DOMAIN")
const FSURDAT = joinpath(DOM_DIR, "settings/CLM/parameters/surfdata_clm.nc")
const PARAMFILE = joinpath(DOM_DIR, "settings/CLM/parameters/clm5_params.nc")
const STAGE = get(ENV, "CLM_FY_STAGE", "combine")
# CLM_SMOOTH=always|never|auto — set SMOOTH_MODE so BOTH the ForwardDiff (Dual) and central-FD
# (Float64) paths evaluate the SAME function. Under the default :auto, Dual smooths but Float64
# is exact, so an FD-vs-AD gap conflates the smooth-vs-hard mismatch with real non-smoothness.
get(ENV,"CLM_SMOOTH","")!="" && (C.SMOOTH_MODE[] = Symbol(ENV["CLM_SMOOTH"]))
const ABSORP = parse(Float64, get(ENV, "CLM_ABSORP", "0.12"))
const NACC = parse(Int, get(ENV, "CLM_NACC", "140"))
const NDAYS = parse(Int, get(ENV, "CLM_NDAYS", "10"))
const ROUTE_K = parse(Float64, get(ENV, "CLM_ROUTEK", "3.0"))

# ---- type-aware Dual retype (from forwarddiff_feasibility.jl) ----
_isclm(x) = (t=typeof(x); isstructtype(t) && parentmodule(t)===CLM && !(x isa AbstractArray))
_hasfree(ty) = ty isa TypeVar || ty isa UnionAll || (ty isa DataType && any(_hasfree, ty.parameters))
_skip_retype(x) = (n=String(nameof(typeof(x))); n!="UrbanParamsData" &&
    (occursin("Config",n)||occursin("Constants",n)||occursin("Overrides",n)||occursin("Params",n)))
const RETYPE_FALLBACK = Set{Symbol}()
function retype_val(::Type{T}, v, dft) where {T}
    if v isa AbstractArray
        if eltype(v) <: AbstractFloat; return _hasfree(dft) ? T.(v) : v
        elseif !isempty(v) && _isclm(first(v)); return map(e->retype_struct(T,e), v)
        else; return v; end
    elseif v isa AbstractFloat; return _hasfree(dft) ? T(v) : v
    elseif _isclm(v); return retype_struct(T, v)
    else; return v; end
end
function retype_struct(::Type{T}, obj) where {T}
    _skip_retype(obj) && (push!(RETYPE_FALLBACK, nameof(typeof(obj))); return obj)
    DT2=typeof(obj); W=DT2.name.wrapper; body=Base.unwrap_unionall(W); dfts=fieldtypes(body); fns=fieldnames(DT2)
    vals=ntuple(i->retype_val(T, getfield(obj,fns[i]), dfts[i]), length(fns))
    any(i->typeof(vals[i])!==typeof(getfield(obj,fns[i])), eachindex(fns)) || return obj
    try; return W(; NamedTuple{fns}(vals)...); catch; push!(RETYPE_FALLBACK, nameof(typeof(obj))); return obj; end
end
retype(::Type{T}, x) where {T} = retype_struct(T, x)

# ---- forcing/obs readers ----
function read_series(var, start::Date, nsteps::Int)
    yr=year(start); ff=joinpath(DOM_DIR,"data/forcing/CLM_input","clmforc.$(yr).nc")
    ds=NCDataset(ff); traw=ds["time"][:]
    times=DateTime[ismissing(t) ? DateTime(0) : DateTime(t) for t in traw]
    v=vec(Float64.(ds[var][:])); close(ds)
    i0=findfirst(>=(DateTime(start)), times); i0===nothing && error("start not in forcing")
    out=zeros(FT,nsteps)
    for k in 1:nsteps; h=i0+(k-1)÷2; h>length(v) && break; x=v[h]; out[k]=FT((isfinite(x)&&x<1e30) ? max(x,0.0) : 0.0); end
    out
end
function read_obs_daily(start::Date, ndays::Int)
    f=joinpath(DOM_DIR,"data/observations/streamflow/preprocessed","$(DOMAIN)_streamflow_processed.csv")
    obs=fill(NaN,ndays); ss=Dict{Date,Float64}(); sc=Dict{Date,Int}()
    for line in readlines(f)[2:end]
        p=split(line,','); length(p)<2 && continue
        dt=tryparse(Date,strip(p[1])[1:min(10,end)]); dt===nothing && continue
        q=tryparse(Float64,strip(p[2])); (q===nothing||!isfinite(q)) && continue
        ss[dt]=get(ss,dt,0.0)+q; sc[dt]=get(sc,dt,0)+1
    end
    for d in 1:ndays; k=start+Day(d-1); haskey(sc,k) && (obs[d]=ss[k]/sc[k]); end
    obs
end

# ---- forcing helper + snowpack builder (Float64 base inst) ----
function set_forcing!(a2l, ng; T, solad, rain, snow)
    for g in 1:ng
        a2l.forc_t_not_downscaled_grc[g]=T; a2l.forc_pbot_not_downscaled_grc[g]=85000.0
        a2l.forc_th_not_downscaled_grc[g]=T*(100000.0/85000.0)^(C.RAIR/C.CPAIR)
        a2l.forc_rho_not_downscaled_grc[g]=85000.0/(C.RAIR*T); a2l.forc_lwrad_not_downscaled_grc[g]=250.0+(T-260)*3
        a2l.forc_vp_grc[g]=300.0; a2l.forc_hgt_grc[g]=30.0; a2l.forc_topo_grc[g]=0.0; a2l.forc_wind_grc[g]=2.0
        a2l.forc_hgt_u_grc[g]=30.0; a2l.forc_hgt_t_grc[g]=30.0; a2l.forc_hgt_q_grc[g]=30.0
        for b in 1:C.NUMRAD; a2l.forc_solad_not_downscaled_grc[g,b]=solad; a2l.forc_solai_grc[g,b]=solad*0.4; end
        a2l.forc_solar_not_downscaled_grc[g]=solad*2.8
        a2l.forc_rain_not_downscaled_grc[g]=rain; a2l.forc_snow_not_downscaled_grc[g]=snow
    end
end
function build_snowpack(start::Date; nacc=NACC)
    (inst,bounds,filt,tm)=C.clm_initialize!(; fsurdat=FSURDAT, paramfile=PARAMFILE)
    config=C.CLMDriverConfig(); fia=C.clump_filter_inactive_and_active
    (declin,_)=C.compute_orbital(Float64(dayofyear(start))); nextsw=dayofyear(start)+DT/C.SECSPDAY
    runstep!(n)=begin
        set_forcing!(inst.atm2lnd, bounds.endg; T=270.8, solad=12.0, rain=0.0, snow=8.0e-4)
        C.downscale_forcings!(bounds,inst.atm2lnd,inst.column,inst.landunit,inst.topo)
        C.clm_drv!(config,inst,filt,fia,bounds,true,nextsw,declin,declin,C.ORB_OBLIQR_DEFAULT,
            false,false,"",false; nstep=n,is_first_step=(n==1),is_beg_curr_day=(n==1),dtime=DT,
            mon=month(start),day=max(1,day(start)-5),photosyns=inst.photosyns)
    end
    for n in 1:nacc; runstep!(n); end
    jz=C.varpar.nlevsno+1
    for c in bounds.begc:bounds.endc
        filt.hydrologyc[c] || continue
        nb=inst.column.nbedrock[c]; inst.soilhydrology.zwt_col[c]=max(0.3, inst.column.zi[c,nb+jz]-1.0)
    end
    inst,bounds,filt,config
end

# ---- eltype-generic phases (read b.theta[1]=θ_snow, [2]=fff, [3]=baseflow) ----
function rain_inject!(b,aux); wf=b.inst.water.waterfluxbulk_inst.wf
    @inbounds for c in aux.cols; wf.qflx_rain_plus_snomelt_col[c]=aux.rain[c]; end; nothing; end
function meltenergy_phase!(b,aux); θ=b.theta[1]; sl=b.inst.solarabs.sabg_lyr_patch; ss=b.inst.solarabs.sabg_snow_patch
    @inbounds for p in axes(sl,1); act=false
        for j in axes(sl,2); if aux.mask[p,j]; sl[p,j]=θ*aux.energy; act=true; end; end
        act && (ss[p]=θ*aux.energy); end; nothing; end
function soiltemp_ft!(b,aux); T=eltype(b.theta); a0=aux.a0
    a=(; urbantv=fill(T(a0.urbantv[1]),length(a0.urbantv)), dtime=a0.dtime, bc_col=a0.bc_col,
        bc_lun=a0.bc_lun, bc_patch=a0.bc_patch, nolakec=a0.nolakec, nolakep=a0.nolakep, urbanl=a0.urbanl, urbanc=a0.urbanc)
    C.soiltemp_rev_phase!(b, a); nothing; end
function snomelt_connect!(b,aux); wf=b.inst.water.waterfluxbulk_inst.wf
    # melt-driven (matches the validated reverse-AD snow thread); rain kept off — on a frozen
    # spring soil, adding rain triggers a NaN in the infiltration-excess/h2osfc path.
    @inbounds for c in aux.cols; wf.qflx_rain_plus_snomelt_col[c]=wf.qflx_snomelt_col[c]; end; nothing; end
function satexcess_cal!(b,aux); i=b.inst; θ1=b.theta[2]; sh=i.soilhydrology; ss=i.soilstate; ser=i.sat_excess_runoff
    wfb=i.water.waterfluxbulk_inst; wf=wfb.wf
    @inbounds for c in aux.cols
        ft=sh.frost_table_col[c]; zw=sh.zwt_col[c]; zwp=sh.zwt_perched_col[c]; wt=ss.wtfact_col[c]
        fs=(ft>zwp && ft<=zw) ? wt*exp(-0.5*θ1*zwp) : wt*exp(-0.5*θ1*zw)
        ser.fsat_col[c]=fs; ser.fcov_col[c]=fs; wfb.qflx_sat_excess_surf_col[c]=fs*wf.qflx_rain_plus_snomelt_col[c]
    end; nothing; end
function baseflow_cal!(b,aux); i=b.inst; θ2=b.theta[3]; sh=i.soilhydrology; col=i.column
    ws=i.water.waterstatebulk_inst.ws; wf=i.water.waterfluxbulk_inst.wf; jz=aux.nlevsno+1
    @inbounds for c in aux.cols
        nb=col.nbedrock[c]; zi_bed=col.zi[c,nb+jz]; zw=sh.zwt_col[c]
        q= zw<=zi_bed ? θ2*tan(FT(pi)/180.0*col.topo_slope[c])*(zi_bed-zw)^aux.n_baseflow : zero(zw)
        wf.qflx_drain_col[c]=q
        jl=min(nb,aux.nlevsoi); nl=ws.h2osoi_liq_col[c,jl+aux.nlevsno]-q*aux.dt
        ws.h2osoi_liq_col[c,jl+aux.nlevsno]= nl>0 ? nl : ws.h2osoi_liq_col[c,jl+aux.nlevsno]
    end; nothing; end
function accum_runoff!(b,aux); wf=b.inst.water.waterfluxbulk_inst.wf; s=zero(eltype(b.series))
    @inbounds for c in aux.cols; s+=wf.qflx_surf_col[c]+wf.qflx_drain_col[c]; end
    b.series[aux.k]+=s; nothing; end

# ---- DISCRETE snow-layer management: combine/divide fire snl TRANSITIONS inside the
# differentiated loop (only wired in when CLM_LAYERS=1). Signatures match the STAGE combine
# probe above. These are the operations whose kink breaks gradient smoothness. ----
function combine_phase!(b,aux); i=b.inst; col=i.column; wsb=i.water.waterstatebulk_inst; ws=wsb.ws
    wd=i.water.waterdiagnosticbulk_inst; wf=i.water.waterfluxbulk_inst.wf; lun=i.landunit
    @inbounds for c in aux.bc; aux.snowc[c]= col.snl[c]<0; end
    C.combine_snow_layers!(col.snl, col.dz, col.zi, col.z, i.temperature.t_soisno_col,
        ws.h2osoi_ice_col, ws.h2osoi_liq_col, ws.h2osno_no_layers_col, wd.snow_depth_col,
        wd.frac_sno_col, wd.frac_sno_eff_col, wsb.int_snow_col, wd.snw_rds_col, i.aerosol,
        lun.itype, lun.urbpoi, col.landunit, aux.snowc, aux.bc, aux.nsno,
        wf.qflx_sl_top_soil_col, aux.dt); nothing; end
function divide_phase!(b,aux); i=b.inst; col=i.column; ws=i.water.waterstatebulk_inst.ws
    wd=i.water.waterdiagnosticbulk_inst
    @inbounds for c in aux.bc; aux.snowc[c]= col.snl[c]<0; end
    C.divide_snow_layers!(col.snl, col.dz, col.zi, col.z, i.temperature.t_soisno_col,
        ws.h2osoi_ice_col, ws.h2osoi_liq_col, wd.frac_sno_col, wd.snw_rds_col, i.aerosol, false,
        aux.snowc, aux.bc, aux.nsno); nothing; end

# ---- build a ForwardDiff season context ----
function make_ff_season(start::Date; ndays::Int, gain_fixed=nothing)
    nsteps=ndays*STEPS_PER_DAY
    inst,bounds,filt,config = build_snowpack(start)
    hcols=Int[c for c in bounds.begc:bounds.endc if filt.hydrologyc[c]]
    fsds=read_series("FSDS", start, nsteps); rain=read_series("PRECTmms", start, nsteps)
    obs=read_obs_daily(start, ndays)
    smask=let sl=inst.solarabs.sabg_lyr_patch
        [isfinite(sl[p,j]) && sl[p,j]!=C.SPVAL for p in axes(sl,1), j in axes(sl,2)]; end
    a_surf=C.surfhydro_rev_aux(bounds,filt; dtime=DT, use_hydrstress=config.use_hydrstress)
    st_a0=C.soiltemp_rev_aux(bounds,filt; dtime=DT)
    sw_aux=C.soilwater_rev_aux(filt,config; dtime=DT); wt_aux=C.watertable_rev_aux(bounds,filt,config; dtime=DT)
    hnd_aux=C.hydnodrain_rev_aux(bounds,filt; dtime=DT)
    bf_aux=(; cols=hcols, nlevsno=C.varpar.nlevsno, nlevsoi=C.varpar.nlevsoi,
             n_baseflow=FT(C.soilhydrology_params.n_baseflow), dt=DT)
    with_layers = get(ENV,"CLM_LAYERS","0")=="1"   # wire combine/divide (snl transitions) into the loop
    lay_aux=(; snowc=falses(length(inst.column.snl)), bc=bounds.begc:bounds.endc, nsno=C.varpar.nlevsno, dt=DT)
    step_phases(k)=Any[
        (meltenergy_phase!, ((; mask=smask, energy=fsds[k]*ABSORP),)),
        (soiltemp_ft!, ((; a0=st_a0),)),
        (snomelt_connect!, ((; cols=hcols, rainval=rain[k]),)),
        (C.setsoilfrac_rev_phase!,(a_surf,)),(C.setfloodc_rev_phase!,(a_surf,)),
        (satexcess_cal!,((; cols=hcols),)),
        (C.setqflx_rev_phase!,(a_surf,)),(C.inflexcess_rev_phase!,(a_surf,)),
        (C.routeinfl_rev_phase!,(a_surf,)),(C.updateh2osfc_rev_phase!,(a_surf,)),
        (C.infil_rev_phase!,(a_surf,)),(C.totalrunoff_rev_phase!,(a_surf,)),
        (C.soilwater_rev_phase!,(sw_aux,)),(C.watertable_rev_phase!,(wt_aux,)),
        (C.hydnodrain_rev_phase!,(hnd_aux,)),
        (baseflow_cal!,(bf_aux,)),(accum_runoff!,((; cols=hcols, k=k),)),
        (with_layers ? Any[(combine_phase!,(lay_aux,)),(divide_phase!,(lay_aux,))] : Any[])...]
    steps=[step_phases(k) for k in 1:nsteps]
    function series_T(θ::AbstractVector{Tt}) where {Tt}
        instT=retype(Tt, deepcopy(inst))    # fresh copy: never mutate the pristine base inst
                                            # (retype(Float64,·) returns its arg, so aliasing without this)
        b=(; inst=instT, scratch=C.cf_rev_scratch(Tt, length(instT.patch.column)),
             theta=collect(Tt,θ), series=zeros(Tt,nsteps))
        for ph in steps, (f,ca) in ph; f(b,ca...); end
        b.series
    end
    daily_mean(s)=[mean(@view s[(d-1)*STEPS_PER_DAY+1:d*STEPS_PER_DAY]) for d in 1:ndays]
    route(q,k)= k<=1.0 ? q : (α=1/k; o=similar(q); o[1]=q[1]; (for i in 2:length(q); o[i]=(1-α)*o[i-1]+α*q[i]; end); o)
    fin=isfinite.(obs)
    gain = gain_fixed!==nothing ? gain_fixed : begin
        m=mean(daily_mean(series_T([1.0,FT(C.sat_excess_runoff_params.fff),FT(C.BASEFLOW_SCALAR[])])))
        om=mean(o for o in obs if isfinite(o)); (isfinite(m)&&m>0&&isfinite(om)) ? om/m : 1.0; end
    sse(s)=sum((route(gain.*daily_mean(s),ROUTE_K)[fin].-obs[fin]).^2)
    loss(θ)=sse(series_T(θ))
    function kge_g(s,g)
        qs=route(g.*daily_mean(s),ROUTE_K)[fin]; o=obs[fin]
        (length(o)<3||std(o)<1e-9||std(qs)<1e-12) && return -Inf
        r=cor(qs,o); 1-sqrt((r-1)^2+(std(qs)/std(o)-1)^2+(mean(qs)/mean(o)-1)^2); end
    own_gain(s)=(m=mean(daily_mean(s)); om=mean(o for o in obs if isfinite(o)); (isfinite(m)&&m>0) ? om/m : gain)
    grad(θ)=ForwardDiff.gradient(loss, θ)
    mkbundle(θ)=(; inst=retype(FT, deepcopy(inst)), scratch=C.cf_rev_scratch(FT, length(inst.patch.column)),
                   theta=collect(FT,θ), series=zeros(FT,nsteps))
    # snl trajectory over the window (Float64 run) to detect discrete snl transitions firing
    function snl_traj(θ)
        b=mkbundle(θ); tr=Int[]
        for ph in steps; for (f,ca) in ph; f(b,ca...); end; push!(tr, b.inst.column.snl[hcols[1]]); end
        tr
    end
    (; inst,bounds,filt,config,hcols,nsteps,ndays,series_T,loss,grad,kge_g,own_gain,gain,fin,
       nfin=count(fin), start, steps, mkbundle, snl_traj)
end

# init once so the param globals (fff/baseflow_scalar) are the basin defaults, not pre-init NaN
let; C.clm_initialize!(; fsurdat=FSURDAT, paramfile=PARAMFILE); end
const FFF0=FT(C.sat_excess_runoff_params.fff); const BFS0=FT(C.BASEFLOW_SCALAR[])

# =============================================================================
if STAGE == "combine"
    println("="^80); println("STAGE combine — combine/divide snow layers with DUAL state (untested spot)"); println("="^80)
    inst,bounds,filt,config = build_snowpack(Date(2010,4,20))
    c0=[c for c in bounds.begc:bounds.endc if filt.hydrologyc[c]][1]
    @printf("snowpack: snl=%d\n", inst.column.snl[c0])
    D=ForwardDiff.Dual{:probe}(1.0,1.0,0.0); T=typeof(D)
    instD=retype(T, inst)
    wsb=instD.water.waterstatebulk_inst; ws=wsb.ws; wd=instD.water.waterdiagnosticbulk_inst; wf=instD.water.waterfluxbulk_inst.wf
    col=instD.column; lun=instD.landunit; nsno=C.varpar.nlevsno
    snowc=falses(length(col.snl)); for c in bounds.begc:bounds.endc; snowc[c]= col.snl[c]<0; end
    ok_c=false; ok_d=false
    try
        C.combine_snow_layers!(col.snl, col.dz, col.zi, col.z, instD.temperature.t_soisno_col,
            ws.h2osoi_ice_col, ws.h2osoi_liq_col, ws.h2osno_no_layers_col, wd.snow_depth_col,
            wd.frac_sno_col, wd.frac_sno_eff_col, wsb.int_snow_col, wd.snw_rds_col, instD.aerosol,
            lun.itype, lun.urbpoi, col.landunit, snowc, bounds.begc:bounds.endc, nsno,
            wf.qflx_sl_top_soil_col, DT)
        ok_c=true; @printf("combine_snow_layers! on Dual: OK  (snl[%d]=%d, eltype dz=%s)\n",
            c0, col.snl[c0], string(eltype(col.dz)))
    catch e; @printf("combine_snow_layers! Dual BREAK: %s\n", sprint(showerror,e)[1:min(160,end)]); end
    try
        C.divide_snow_layers!(col.snl, col.dz, col.zi, col.z, instD.temperature.t_soisno_col,
            ws.h2osoi_ice_col, ws.h2osoi_liq_col, wd.frac_sno_col, wd.snw_rds_col, instD.aerosol, false,
            snowc, bounds.begc:bounds.endc, nsno)
        ok_d=true; @printf("divide_snow_layers! on Dual: OK  (snl[%d]=%d)\n", c0, col.snl[c0])
    catch e; @printf("divide_snow_layers! Dual BREAK: %s\n", sprint(showerror,e)[1:min(140,end)]); end
    # --- ACCUMULATION side: handle_new_snow! (new snowfall → fresh layer) + snow_compaction! ---
    a2l=instD.atm2lnd; ok_h=false; ok_sc=false
    try
        for c in bounds.begc:bounds.endc; a2l.forc_snow_downscaled_col[c]=1e-3; a2l.forc_t_downscaled_col[c]=270.0; end
        snl0=col.snl[c0]
        C.handle_new_snow!(instD.temperature, wsb, wd, col, lun, filt.nolakec, bounds.begc:bounds.endc, DT, nsno;
            forc_t=a2l.forc_t_downscaled_col, forc_wind=a2l.forc_wind_grc,
            qflx_snow_grnd=wf.qflx_snow_grnd_col, qflx_snow_drain=wf.qflx_snow_drain_col,
            int_snow=wsb.int_snow_col, scf_method=instD.scf_method)
        ok_h=true; @printf("handle_new_snow! on Dual: OK  (snl %d→%d, eltype snowdepth=%s)\n",
            snl0, col.snl[c0], string(eltype(wd.snow_depth_col)))
    catch e; @printf("handle_new_snow! Dual BREAK: %s\n", sprint(showerror,e)[1:min(160,end)]); end
    try
        C.snow_compaction!(col.dz, DT, col.snl, instD.temperature.t_soisno_col, ws.h2osoi_ice_col, ws.h2osoi_liq_col,
            instD.temperature.imelt_col, wd.frac_sno_col, wd.frac_h2osfc_col, wd.swe_old_col, wsb.int_snow_col,
            wd.frac_iceold_col, a2l.forc_wind_grc, col.gridcell, col.landunit, lun.lakpoi, lun.urbpoi,
            snowc, bounds.begc:bounds.endc, nsno; use_subgrid_fluxes=false,
            n_melt=fill(20.0,length(col.snl)), int_snow_max=2000.0)
        ok_sc=true; @printf("snow_compaction! on Dual: OK  (eltype dz=%s)\n", string(eltype(col.dz)))
    catch e; @printf("snow_compaction! Dual BREAK: %s\n", sprint(showerror,e)[1:min(160,end)]); end
    @printf("STAGE combine: %s\n", (ok_c && ok_d && ok_h && ok_sc) ?
        "PASS — combine/divide + handle_new_snow + compaction all run with Duals (forward-mode)" : "see results above")
end

if STAGE == "gate"
    println("="^80); println("STAGE gate — ForwardDiff windowed melt gradient vs central FD + wall-time"); println("="^80)
    for nd in (parse(Int, get(ENV,"CLM_NDAYS","10")),)
        s=make_ff_season(Date(2010,4,20); ndays=nd)
        @printf("N=%d steps (%d days) hydro-cols=%d fin-obs=%d\n", s.nsteps, nd, length(s.hcols), s.nfin)
        θp=[1.2, FFF0*1.15, BFS0*1.3]
        # primal finiteness probe (Float64 series through the retype path)
        sp=s.series_T(Float64.(θp)); nf=count(isfinite,sp); fk=findfirst(!isfinite, sp)
        @printf("primal series: finite %d/%d  first-NaN step=%s  loss=%.4e\n", nf, length(sp), string(fk), s.loss(Float64.(θp)))
        if fk !== nothing   # walk step 1 phase-by-phase to find the culprit
            b=s.mkbundle(Float64.(θp)); c=s.hcols[1]
            for (pi,(f,ca)) in enumerate(s.steps[1])
                f(b, ca...)
                wfb=b.inst.water.waterfluxbulk_inst; wf=wfb.wf
                sh=b.inst.soilhydrology
                chk=(qse=wfb.qflx_sat_excess_surf_col[c], fsat=b.inst.sat_excess_runoff.fsat_col[c],
                     wtf=b.inst.soilstate.wtfact_col[c], zwt=sh.zwt_col[c], zwp=sh.zwt_perched_col[c],
                     ftbl=sh.frost_table_col[c], qrps=wf.qflx_rain_plus_snomelt_col[c])
                bad=[k for k in keys(chk) if !isfinite(chk[k])]
                !isempty(bad) && (@printf("  step1 phase %d (%s) → NaN in %s\n    values: %s\n", pi, string(nameof(f)), string(bad),
                    join([@sprintf("%s=%.4g",k,chk[k]) for k in keys(chk)], "  ")); break)
            end
        end
        t=@elapsed g=s.grad(θp)
        @printf("ForwardDiff dL/dθ=[snow=% .5e, fff=% .5e, bf=% .5e]  wall=%.1f s/gradient\n", g[1],g[2],g[3], t)
        names=["θ_snow","fff","baseflow"]; hs=[1e-3,1e-3*FFF0,1e-3*BFS0]; passes=Bool[]
        for j in 1:3
            θa=copy(θp); θa[j]+=hs[j]; θb=copy(θp); θb[j]-=hs[j]
            f1=(s.loss(θa)-s.loss(θb))/(2hs[j]); θa2=copy(θp); θa2[j]+=hs[j]/2; θb2=copy(θp); θb2[j]-=hs[j]/2
            f2=(s.loss(θa2)-s.loss(θb2))/(hs[j]); fd=(4*f2-f1)/3
            rel=abs(g[j]-fd)/max(abs(fd),1e-30); vac=abs(g[j])<1e-30&&abs(fd)<1e-30; pass=(rel<1e-5)&&!vac; push!(passes,pass)
            @printf("  %-9s FwdDiff=% .6e FD(Rich)=% .6e rel=%.2e %s\n", names[j], g[j], fd, rel, vac ? "VACUOUS" : (pass ? "PASS" : "CHECK"))
        end
        @printf("STAGE gate (N=%d): %s\n", s.nsteps, all(passes) ? "PASS ✓" : "see above")
    end
end

if STAGE == "split"
    println("="^80); println("STAGE split — ForwardDiff melt-season split-sample (AD vs DDS)"); println("="^80)
    CAL=Date(get(ENV,"CLM_CAL_START","2010-04-20")); EV=Date(get(ENV,"CLM_EVAL_START","2011-04-20"))
    @printf("cal=%s eval=%s NDAYS=%d\n", CAL, EV, NDAYS)
    cal=make_ff_season(CAL; ndays=NDAYS); evc=make_ff_season(EV; ndays=NDAYS, gain_fixed=cal.gain)
    lo=[1e-3,1e-3,1e-5]; hi=[5.0,2.0,0.05]; sc=[1.0,FFF0,BFS0]
    θ0=clamp.([1.4, FFF0*1.4, BFS0*1.96], lo, hi)
    eval_kge(θ)=(s=evc.series_T(θ); (evc.kge_g(s,cal.gain), evc.kge_g(s,evc.own_gain(s))))
    pk=eval_kge([1.0,FFF0,BFS0])
    @printf("prior cal-KGE=%.4f eval strict/bc=%.4f/%.4f | detuned start cal-KGE=%.4f\n",
        cal.kge_g(cal.series_T([1.0,FFF0,BFS0]),cal.gain), pk[1], pk[2], cal.kge_g(cal.series_T(θ0),cal.gain))
    # gate on cal window
    g=cal.grad([1.2,FFF0*1.15,BFS0*1.3]); @printf("[gate] ForwardDiff cal-grad=[%.3e,%.3e,%.3e]\n", g...)
    # AD LM diagonal-damped Newton
    let θ=copy(θ0); npass=0; hstep=1e-3.*sc
        for it in 1:parse(Int,get(ENV,"CLM_AD_ITER","4"))
            g0=cal.grad(θ); npass+=1
            gg=[cal.grad(θ.+((1:3).==j).*hstep[j]) for j in 1:3]; npass+=3
            H=[(gg[j][i]-g0[i])/hstep[j] for i in 1:3, j in 1:3]; Hs=0.5.*(H+H'); Dg=abs.([Hs[i,i] for i in 1:3]).+1e-30
            μ=1e-2; Lc=cal.loss(θ); acc=false
            for _t in 1:20; A=Hs+μ.*Diagonal(Dg); Δ=-(A\g0); θt=clamp.(θ.+Δ,lo,hi); (cal.loss(θt)<Lc)&&(θ=θt;acc=true;break); μ*=3; end
            @printf("  AD it %2d θ=[%.3f,%.3f,%.2e] cal-KGE=%.4f passes=%d\n", it,θ[1],θ[2],θ[3],cal.kge_g(cal.series_T(θ),cal.gain),npass)
            acc || break
        end
        ad_cal=cal.kge_g(cal.series_T(θ),cal.gain); (ad_es,ad_ebc)=eval_kge(θ)
        # DDS
        rng=MersenneTwister(11); best=copy(θ0); bestL=cal.loss(best); evals=1
        for it in 1:parse(Int,get(ENV,"CLM_DDS_ITER","50"))
            P=1-log(it)/log(50); cand=copy(best)
            for j in 1:3; rand(rng)<P && (cand[j]=clamp(best[j]+0.2*(hi[j]-lo[j])*randn(rng),lo[j],hi[j])); end
            L=cal.loss(cand); evals+=1; L<bestL && (bestL=L;best=cand)
        end
        dds_cal=cal.kge_g(cal.series_T(best),cal.gain); (dds_es,dds_ebc)=eval_kge(best)
        println("\n", "="^80)
        @printf("FORWARDDIFF MELT-SEASON SPLIT-SAMPLE (%s: cal %s → eval %s, %d-day)\n", DOMAIN,
            Dates.format(CAL,"yyyy-mm-dd"), Dates.format(EV,"yyyy-mm-dd"), NDAYS)
        @printf("%-6s | %-8s | %-11s | %-11s | %-10s\n","which","cal-KGE","eval-strict","eval-bc","evals")
        println("-"^80)
        @printf("%-6s | %-+8.4f | %-+11.4f | %-+11.4f | %d ForwardDiff-passes\n","AD",ad_cal,ad_es,ad_ebc,npass)
        @printf("%-6s | %-+8.4f | %-+11.4f | %-+11.4f | %d fwd-evals\n","DDS",dds_cal,dds_es,dds_ebc,evals)
        println("="^80)
    end
end

# =============================================================================
# STAGE sweep — ADVERSARIAL: does the ForwardDiff gradient survive discrete snl/imelt
# TRANSITIONS? Sweep θ_snow over a range that drives the pack through combine/divide (snl
# changes) and imelt freeze/thaw flips. At each θ: FwdDiff gradient vs CENTRAL FD (central FD
# re-runs the whole loop at θ±h so it is the ground truth for the ACTUAL kinked function).
# A spike in |FwdDiff-FD| where snl/imelt changes = the gradient is only piecewise-valid.
# Run with CLM_LAYERS=1 to include combine/divide (snl transitions); =0 for imelt-only.
# =============================================================================
if STAGE == "sweep"
    println("="^92); @printf("STAGE sweep — FwdDiff vs central FD across snl/imelt transitions (LAYERS=%s)\n", get(ENV,"CLM_LAYERS","0")); println("="^92)
    s=make_ff_season(Date(2010,4,20); ndays=parse(Int,get(ENV,"CLM_NDAYS","10")))
    fff0=FT(C.sat_excess_runoff_params.fff); bf0=FT(C.BASEFLOW_SCALAR[])
    lo=parse(Float64,get(ENV,"CLM_SNOW_LO","0.5")); hi=parse(Float64,get(ENV,"CLM_SNOW_HI","6.0"))
    nθ=parse(Int,get(ENV,"CLM_NSWEEP","28")); h=parse(Float64,get(ENV,"CLM_FDH","1e-4"))
    @printf("θ_snow sweep [%.2f,%.2f] × %d,  FD h=%.1e,  N=%d steps\n", lo, hi, nθ, h, s.nsteps)
    @printf("%-8s %-9s %-13s %-13s %-9s %-18s %s\n","θ_snow","loss","FwdDiff","centralFD","rel","snl(θ) uniq","snl-change@θ±h?")
    println("-"^92)
    worst=Ref((0.0, 0.0))   # Ref: top-level `for` is soft-scope, avoid the local-shadow gotcha
    for θs in range(lo, hi, length=nθ)
        θ=[θs, fff0, bf0]
        g=s.grad(θ)[1]
        fd=(s.loss([θs+h,fff0,bf0])-s.loss([θs-h,fff0,bf0]))/(2h)
        rel=abs(g-fd)/max(abs(fd),1e-30)
        tr=s.snl_traj(θ); uq=unique(tr)
        trp=s.snl_traj([θs+h,fff0,bf0]); trm=s.snl_traj([θs-h,fff0,bf0])
        changed = (tr != trp) || (tr != trm)   # snl trajectory differs across the FD bracket?
        @printf("%-8.3f %-9.3e % -12.5e % -12.5e %-9.2e %-18s %s\n", θs, s.loss(θ), g, fd, rel,
            string(uq), changed ? "YES <- STRADDLE" : "no")
        rel>worst[][2] && (worst[]=(θs, rel))
    end
    println("-"^92)
    @printf("worst FwdDiff-vs-FD rel over sweep = %.2e at θ_snow=%.3f\n", worst[][2], worst[][1])
    println("VERDICT: ", worst[][2]<0.05 ?
        "transitions gradient-BENIGN (FwdDiff==FD across snl/imelt changes → full-year is genuinely assembly)" :
        "FwdDiff DIVERGES from FD across transitions (gradient only piecewise-valid → needs smoothing; quantified above)")
end
