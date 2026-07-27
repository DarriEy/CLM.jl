# =============================================================================
# FEASIBILITY PROBE: does the SURFACE ENERGY BALANCE (canopy_fluxes! / friction velocity +
# turbulent fluxes + photosynthesis coupling) run Dual-clean under ForwardDiff?
#
# This settles the open question: a physical full-year hydrograph needs eflx_sh/eflx_lh from
# the surface energy balance (which soil_temperature! consumes for melt). Is that path
# Dual-clean, or is photosynthesis (Farquhar) / LUNA the whole-model forward-AD wall?
#
# Runs cf_rev_phases (the canopy energy-balance sub-phases) on a DUAL-typed inst, phase by
# phase, catching the FIRST break. Two configs: use_psn=false (pure energy balance: init →
# Monin-Obukhov friction → resistances → energy) and use_psn=true (+ photosynthesis/LUNA).
#
#   CLM_PSN=0|1  julia --project=. scripts/forwarddiff_canopy_probe.jl
# =============================================================================
using CLM, ForwardDiff, Printf
const C = CLM
const FT = Float64
const DT = 1800.0
const FSURDAT = get(ENV, "CLM_FSURDAT",
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc")
const PARAMFILE = get(ENV, "CLM_PARAMFILE",
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc")
const USE_PSN = get(ENV, "CLM_PSN", "0") == "1"

# ---- type-aware Dual retype (identical to forwarddiff_feasibility.jl / _fullyear.jl) ----
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

# ---- warm a real vegetated inst (summer, canopy active) ----
function setup_forcing!(a2l, T, ng)
    for g in 1:ng
        a2l.forc_t_not_downscaled_grc[g]=T; a2l.forc_pbot_not_downscaled_grc[g]=85000.0
        a2l.forc_th_not_downscaled_grc[g]=T*(100000.0/85000.0)^(C.RAIR/C.CPAIR)
        a2l.forc_rho_not_downscaled_grc[g]=85000.0/(C.RAIR*T); a2l.forc_lwrad_not_downscaled_grc[g]=340.0
        a2l.forc_vp_grc[g]=1500.0; a2l.forc_hgt_grc[g]=30.0; a2l.forc_topo_grc[g]=0.0; a2l.forc_wind_grc[g]=3.0
        a2l.forc_hgt_u_grc[g]=30.0; a2l.forc_hgt_t_grc[g]=30.0; a2l.forc_hgt_q_grc[g]=30.0
        for b in 1:C.NUMRAD; a2l.forc_solad_not_downscaled_grc[g,b]=300.0; a2l.forc_solai_grc[g,b]=120.0; end
        a2l.forc_solar_not_downscaled_grc[g]=840.0; a2l.forc_rain_not_downscaled_grc[g]=0.0
        a2l.forc_snow_not_downscaled_grc[g]=0.0
    end
end
function build()
    (inst,bounds,filt,tm)=C.clm_initialize!(; fsurdat=FSURDAT, paramfile=PARAMFILE)
    setup_forcing!(inst.atm2lnd, 298.0, bounds.endg)
    C.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
    config=C.CLMDriverConfig(); fia=C.clump_filter_inactive_and_active
    (declin,_)=C.compute_orbital(200.0); nextsw=200.0+DT/C.SECSPDAY
    C.interp_monthly_veg!(inst.satellite_phenology; kmo=7, kda=15)
    cs=inst.canopystate; wdb=inst.water.waterdiagnosticbulk_inst
    C.satellite_phenology!(inst.satellite_phenology, cs, wdb, inst.patch, filt.nolakep, bounds.begp:bounds.endp)
    for p in bounds.begp:bounds.endp; cs.frac_veg_nosno_patch[p]=cs.frac_veg_nosno_alb_patch[p]; end
    C.set_exposedvegp_filter!(filt, bounds, cs.frac_veg_nosno_patch)
    for n in 1:3
        C.clm_drv!(config,inst,filt,fia,bounds,true,nextsw,declin,declin,C.ORB_OBLIQR_DEFAULT,
            false,false,"",false; nstep=n,is_first_step=(n==1),is_beg_curr_day=(n==1),dtime=DT,
            mon=7,day=15,photosyns=inst.photosyns)
    end
    inst,bounds,filt,config
end

inst,bounds,filt,config = build()
@printf("PROBE: canopy surface-energy-balance under ForwardDiff  use_psn=%s\n", USE_PSN)
np = length(inst.patch.column)
D = ForwardDiff.Dual{:probe}(1.0,1.0,0.0); T = typeof(D)

# retype inst to Dual + build canopy aux FROM the Dual inst (so its Const copies are Dual)
instD = retype(T, inst)
@printf("retype OK: eltype(t_veg)=%s  structs left Float64: %d\n",
    string(eltype(instD.temperature.t_veg_patch)), length(RETYPE_FALLBACK))
# Dual-aware canopy aux (mirrors driver_reverse.jl canopy_rev_aux but with eltype T instead of
# the hardcoded Float64 — the aux BUILDER's Float64 pins are harness-side, not the physics).
function canopy_aux_T(inst, bounds, filt, ::Type{Tp}; use_psn, dtime) where {Tp}
    NP=bounds.endp; ev=filt.exposedvegp; fp=Int[p for p in bounds.begp:bounds.endp if ev[p]]; fn=length(fp)
    a2l=inst.atm2lnd; MP=C.MXPFT+1
    forc_q_col=fill!(similar(a2l.forc_pbot_downscaled_col), zero(Tp))
    C.compute_forc_q!(forc_q_col, inst.column.gridcell, a2l.forc_vp_grc, a2l.forc_pbot_downscaled_col)
    forc=(; lwrad=a2l.forc_lwrad_downscaled_col, q=forc_q_col, pbot=a2l.forc_pbot_downscaled_col,
           th=a2l.forc_th_downscaled_col, rho=a2l.forc_rho_downscaled_col, t=a2l.forc_t_downscaled_col,
           u_grc=a2l.forc_u_grc, v_grc=a2l.forc_v_grc, pco2=a2l.forc_pco2_grc, po2=a2l.forc_po2_grc,
           hgt_t=a2l.forc_hgt_t_grc, hgt_u=a2l.forc_hgt_u_grc, hgt_q=a2l.forc_hgt_q_grc,
           dayl=inst.gridcell.dayl, max_dayl=inst.gridcell.max_dayl,
           downreg=fill(Tp(1.0),NP), leafn=fill(Tp(1.0),NP))
    pft=(; dleaf=fill(Tp(0.04),MP), z0v_Cr=fill(Tp(0.35),MP), z0v_Cs=fill(Tp(0.003),MP),
          z0v_c=fill(Tp(0.25),MP), z0v_cw=fill(Tp(2.0),MP), z0v_LAImax=fill(Tp(8.0),MP), grnd_ch4=fill(Tp(0.0),NP))
    psn = if use_psn
        sa=inst.surfalb; so=inst.solarabs; cs2=inst.canopystate
        (; c3psn=Tp.(C.pftcon.c3psn), leafcn=Tp.(C.pftcon.leafcn), flnr=Tp.(C.pftcon.flnr),
           fnitr=Tp.(C.pftcon.fnitr), slatop=Tp.(C.pftcon.slatop), mbbopt=Tp.(C.pftcon.mbbopt),
           medlynintercept=fill(Tp(100.0),MP), medlynslope=fill(Tp(6.0),MP),
           nrad=copy(sa.nrad_patch), tlai_z=Tp.(sa.tlai_z_patch), parsun_z=Tp.(so.parsun_z_patch),
           parsha_z=Tp.(so.parsha_z_patch), laisun_z=Tp.(cs2.laisun_z_patch), laisha_z=Tp.(cs2.laisha_z_patch),
           vcmaxcint_sun=Tp.(sa.vcmaxcintsun_patch), vcmaxcint_sha=Tp.(sa.vcmaxcintsha_patch),
           o3coefv=fill(Tp(1.0),NP), o3coefg=fill(Tp(1.0),NP), t10=Tp.(inst.temperature.t_a10_patch))
    else
        (; c3psn=fill(Tp(1.0),MP), leafcn=fill(Tp(25.0),MP), flnr=fill(Tp(0.1),MP), fnitr=fill(Tp(1.0),MP),
           slatop=fill(Tp(0.01),MP), mbbopt=fill(Tp(9.0),MP), medlynintercept=fill(Tp(100.0),MP),
           medlynslope=fill(Tp(6.0),MP), nrad=fill(1,NP), tlai_z=fill(Tp(1.0),NP,C.NLEVCAN),
           parsun_z=fill(Tp(0.0),NP,C.NLEVCAN), parsha_z=fill(Tp(0.0),NP,C.NLEVCAN),
           laisun_z=fill(Tp(0.0),NP,C.NLEVCAN), laisha_z=fill(Tp(0.0),NP,C.NLEVCAN),
           vcmaxcint_sun=fill(Tp(1.0),NP), vcmaxcint_sha=fill(Tp(0.6),NP), o3coefv=fill(Tp(1.0),NP),
           o3coefg=fill(Tp(1.0),NP), t10=Tp.(inst.temperature.t_a10_patch))
    end
    return (; patch=inst.patch, col=inst.column, grid=inst.gridcell, forc=forc, pft=pft, psn=psn,
        filterp=fp, fn=fn, active=Bool[ev[p] for p in 1:NP], mask=ev, ivt=inst.patch.itype .+ 1,
        forc_pbot_patch=Tp[a2l.forc_pbot_downscaled_col[inst.patch.column[p]] for p in 1:NP],
        soilevap_beta=C.do_soilevap_beta(), soil_resis_sl14=C.do_soil_resistance_sl14(),
        nlevsno=C.varpar.nlevsno, dtime=Tp(dtime), use_psn=use_psn)
end
aux = try canopy_aux_T(instD, bounds, filt, T; use_psn=USE_PSN, dtime=DT)
    catch e; println("canopy_aux_T BREAK: ", sprint(showerror,e)[1:min(160,end)]); nothing end

if aux !== nothing
    b = (; inst=instD, scratch=C.cf_rev_scratch(T, np))
    cv = C.cf_rev_bundle(instD.canopystate, instD.energyflux, instD.frictionvel, instD.temperature,
        instD.solarabs, instD.soilstate, instD.water.waterfluxbulk_inst, instD.water.waterstatebulk_inst,
        instD.water.waterdiagnosticbulk_inst, instD.photosyns, b.scratch)
    phases = C.cf_rev_phases(aux, 3)
    @printf("running %d canopy sub-phases on Dual state ...\n", length(phases))
    firstbreak = nothing; bt = nothing
    for (pi,(f,ca)) in enumerate(phases)
        try
            f(cv, ca...)
        catch e
            firstbreak = (pi, nameof(f), sprint(showerror, e)); bt = stacktrace(catch_backtrace())
            break
        end
    end
    if firstbreak === nothing
        p = inst.patch.column === nothing ? 1 : 1
        ef = instD.energyflux
        @printf("ALL %d canopy phases ran Dual-clean ✓  (eflx_sh_grnd eltype=%s, t_veg eltype=%s)\n",
            length(phases), string(eltype(ef.eflx_sh_grnd_patch)), string(eltype(instD.temperature.t_veg_patch)))
        println("→ VERDICT: surface energy balance is Dual-clean under ForwardDiff (use_psn=$USE_PSN)")
        # GATE: d(sum eflx_sh_grnd)/d(wind-forcing scale θ) — ForwardDiff vs central FD.
        expp = [p for p in bounds.begp:bounds.endp if filt.exposedvegp[p]]
        bt = copy(inst.atm2lnd.forc_t_downscaled_col); bth = copy(inst.atm2lnd.forc_th_downscaled_col)
        function tveg_of(θ::AbstractVector{Q}) where {Q}   # melt-relevant: leaf/ground temperature
            iT = retype(Q, deepcopy(inst))
            iT.atm2lnd.forc_t_downscaled_col .= θ[1] .* bt; iT.atm2lnd.forc_th_downscaled_col .= θ[1] .* bth
            a = canopy_aux_T(iT, bounds, filt, Q; use_psn=USE_PSN, dtime=DT)
            bb = (; inst=iT, scratch=C.cf_rev_scratch(Q, np))
            vv = C.cf_rev_bundle(iT.canopystate, iT.energyflux, iT.frictionvel, iT.temperature,
                iT.solarabs, iT.soilstate, iT.water.waterfluxbulk_inst, iT.water.waterstatebulk_inst,
                iT.water.waterdiagnosticbulk_inst, iT.photosyns, bb.scratch)
            for (f2,ca2) in C.cf_rev_phases(a, 3); f2(vv, ca2...); end
            s = zero(Q); for p in expp; s += iT.temperature.t_veg_patch[p]; end; s
        end
        g = ForwardDiff.gradient(tveg_of, [1.0])[1]
        fd = (tveg_of([1.0+1e-5]) - tveg_of([1.0-1e-5]))/2e-5
        rel = abs(g-fd)/max(abs(fd),1e-30)
        @printf("[GATE] d(Σ t_veg)/d(air-temp scale): FwdDiff=% .6e  FD=% .6e  rel=%.2e  %s\n",
            g, fd, rel, rel<1e-5 ? "PASS ✓" : "CHECK")
        # ---- psn-DEGENERACY check: is the "Farquhar Dual-clean" claim NON-VACUOUS? ----
        # The t_veg gate flows through the ENERGY BALANCE; it does NOT prove the derivative goes
        # through the Farquhar solve. Gate a genuinely photosynthesis-dependent output (GPP =
        # Σ psnsun). If GPP≈0 in the probe state, d(GPP)/dθ = 0 → VACUOUS: psn runs without a
        # type-break but is un-exercised, so its Dual-cleanliness is UNCONFIRMED.
        if USE_PSN
            function gpp_of(θ::AbstractVector{Q}) where {Q}
                iT = retype(Q, deepcopy(inst))
                iT.atm2lnd.forc_t_downscaled_col .= θ[1] .* bt; iT.atm2lnd.forc_th_downscaled_col .= θ[1] .* bth
                a = canopy_aux_T(iT, bounds, filt, Q; use_psn=true, dtime=DT)
                bb = (; inst=iT, scratch=C.cf_rev_scratch(Q, np))
                vv = C.cf_rev_bundle(iT.canopystate, iT.energyflux, iT.frictionvel, iT.temperature,
                    iT.solarabs, iT.soilstate, iT.water.waterfluxbulk_inst, iT.water.waterstatebulk_inst,
                    iT.water.waterdiagnosticbulk_inst, iT.photosyns, bb.scratch)
                for (f2,ca2) in C.cf_rev_phases(a, 3); f2(vv, ca2...); end
                s = zero(Q); for p in expp; s += iT.photosyns.psnsun_patch[p]; end; s
            end
            gpp0 = gpp_of([1.0]); gg = ForwardDiff.gradient(gpp_of, [1.0])[1]
            @printf("[psn-degeneracy] GPP=Σpsnsun=% .4e  d(GPP)/d(air-temp)=% .4e  →  %s\n", gpp0, gg,
                abs(gg)<1e-30 ? "VACUOUS — psn un-exercised, Farquhar Dual-cleanliness UNCONFIRMED" :
                                "non-zero — Farquhar IS exercised (Dual-clean confirmed)")
        end
    else
        pi,fn,msg = firstbreak
        @printf("FIRST BREAK at canopy phase %d (%s):\n  %s\n", pi, string(fn), msg[1:min(260,end)])
        println("  --- stack (top frames) ---")
        for (i,fr) in enumerate(bt); i>12 && break; println("    ", fr); end
    end
end
