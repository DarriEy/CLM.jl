# =============================================================================
# STAGE A — ForwardDiff (forward-mode AD) FEASIBILITY GATE for a year-scale streamflow
# gradient. Reverse-AD (Enzyme + multistep_reverse!) is tape-bounded (~30-day window max,
# validated in calibrate_streamflow_ad_real.jl). Forward-mode dual numbers have NO tape
# (memory constant in time, cost ~(1+n_params) forward passes), so IF the physics is
# Dual-generic it can reach a full year. This probe tests that "IF".
#
# Route: the CLM state structs are parametric ({FT<:Real, V<:AbstractVector{FT}, …}) and
# Adapt-@adapt_structure'd (the Metal Float32 device path uses this). So a custom Adapt
# adaptor can re-type the whole inst tree from Float64 → ForwardDiff.Dual, INSIDE the
# differentiated function (so the Dual tag matches ForwardDiff's). We then run the SAME
# real-forcing hydrology thread and compare ForwardDiff vs central FD vs reverse-AD.
#
#   CLM_NSTEPS=48 julia --project=. scripts/forwarddiff_feasibility.jl
# =============================================================================
using CLM, ForwardDiff, Enzyme, Adapt, Printf, Statistics, Dates, NCDatasets
const C = CLM
const FT = Float64
const DT = 1800.0
const ROOT = get(ENV, "SYMFLUENCE_DATA",
    "/Users/darri.eythorsson/Library/CloudStorage/GoogleDrive-dareyt@gmail.com/My Drive/data/SYMFLUENCE_data")
const DOMAIN = get(ENV, "CLM_DOMAIN", "Stillwater_Oklahoma")
const NSTEPS = parse(Int, get(ENV, "CLM_NSTEPS", "48"))
const START  = Date(get(ENV, "CLM_START", "2004-06-01"))
const DOM_DIR = joinpath(ROOT, "domain_$DOMAIN")
const FSURDAT = joinpath(DOM_DIR, "settings/CLM/parameters/surfdata_clm.nc")
const PARAMFILE = joinpath(DOM_DIR, "settings/CLM/parameters/clm5_params.nc")

# ---- eltype adaptor: Float64 arrays → eltype T; everything else (Int/Bool arrays,
#      scalars, flags) passes through unchanged. ----
# TYPE-AWARE retype: the CLM state structs mix parametric (::FT, ::V, ::M) and CONCRETE
# (::Float64) float fields (e.g. TemperatureData.excess_ice_coldstart_depth::FT must become
# Dual, but CanopyStateData.leaf_mr_vcm::Float64 must STAY Float64). A value-based adaptor
# can't tell them apart, so we inspect each field's DECLARED type in the struct's UnionAll
# wrapper: a TypeVar (FT/V/M) → convert to T; a concrete Float64 → keep.
_isclm(x) = (t = typeof(x); isstructtype(t) && parentmodule(t) === CLM && !(x isa AbstractArray))
_hasfree(ty) = ty isa TypeVar || (ty isa UnionAll) ||
    (ty isa DataType && any(_hasfree, ty.parameters))
function retype_val(::Type{T}, v, dft) where {T}
    if v isa AbstractArray
        if eltype(v) <: AbstractFloat
            return _hasfree(dft) ? T.(v) : v
        elseif !isempty(v) && _isclm(first(v))
            return map(e -> retype_struct(T, e), v)      # array of CLM structs
        else
            return v                                      # Int/Bool/etc arrays
        end
    elseif v isa AbstractFloat
        return _hasfree(dft) ? T(v) : v
    elseif _isclm(v)
        return retype_struct(T, v)
    else
        return v
    end
end
const RETYPE_FALLBACK = Set{Symbol}()   # structs left Float64 (config/no-constructor)
# Config structs that are NOT differentiated state — leave Float64 (θ is seeded via the phase
# functions, not through these). CLMInstances fields are abstract-typed so the mix is accepted.
_skip_retype(x) = (n = String(nameof(typeof(x))); n != "UrbanParamsData" &&
                    (occursin("Config", n) || occursin("Constants", n) ||
                     occursin("Overrides", n) || occursin("Params", n)))
function retype_struct(::Type{T}, obj) where {T}
    _skip_retype(obj) && (push!(RETYPE_FALLBACK, nameof(typeof(obj))); return obj)
    DT = typeof(obj); W = DT.name.wrapper
    body = Base.unwrap_unionall(W); dfts = fieldtypes(body); fns = fieldnames(DT)
    vals = ntuple(i -> retype_val(T, getfield(obj, fns[i]), dfts[i]), length(fns))
    any(i -> typeof(vals[i]) !== typeof(getfield(obj, fns[i])), eachindex(fns)) || return obj
    try
        return W(; NamedTuple{fns}(vals)...)             # @kwdef KEYWORD constructor (wall-3 fix)
    catch e
        push!(RETYPE_FALLBACK, nameof(typeof(obj)))      # no usable constructor → leave Float64
        return obj
    end
end
retype(::Type{T}, x) where {T} = retype_struct(T, x)

# ---- real forcing window ----
function read_forcing_window(start::Date, nsteps::Int)
    yr = year(start); ff = joinpath(DOM_DIR, "data/forcing/CLM_input", "clmforc.$(yr).nc")
    ds = NCDataset(ff)
    traw = ds["time"][:]; times = DateTime[ismissing(t) ? DateTime(0) : DateTime(t) for t in traw]
    prec = vec(Float64.(ds["PRECTmms"][:])); close(ds)
    t0 = DateTime(start); i0 = findfirst(>=(t0), times); i0 === nothing && (i0 = 1)
    rain = zeros(FT, nsteps)
    for k in 1:nsteps
        h = i0 + (k-1) ÷ 2; h > length(prec) && break
        r = prec[h]; rain[k] = FT((isfinite(r) && r < 1e30) ? max(r,0.0) : 0.0)
    end
    return rain
end

# ---- phase fns (generic in eltype: read b.theta, write into b.inst arrays) ----
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

# ---- build Float64 inst + warmup + wet water table (as in the real script) ----
function build()
    (inst, bounds, filt, tm) = C.clm_initialize!(; fsurdat=FSURDAT, paramfile=PARAMFILE)
    rain0 = read_forcing_window(START, NSTEPS)
    for g in 1:bounds.endg
        a2l=inst.atm2lnd; T=290.0
        a2l.forc_t_not_downscaled_grc[g]=T; a2l.forc_pbot_not_downscaled_grc[g]=95000.0
        a2l.forc_th_not_downscaled_grc[g]=T*(100000.0/95000.0)^(C.RAIR/C.CPAIR)
        a2l.forc_rho_not_downscaled_grc[g]=95000.0/(C.RAIR*T); a2l.forc_lwrad_not_downscaled_grc[g]=320.0
        a2l.forc_vp_grc[g]=1200.0; a2l.forc_hgt_grc[g]=30.0; a2l.forc_topo_grc[g]=0.0; a2l.forc_wind_grc[g]=3.0
        a2l.forc_hgt_u_grc[g]=30.0; a2l.forc_hgt_t_grc[g]=30.0; a2l.forc_hgt_q_grc[g]=30.0
        for b in 1:C.NUMRAD; a2l.forc_solad_not_downscaled_grc[g,b]=200.0; a2l.forc_solai_grc[g,b]=80.0; end
        a2l.forc_solar_not_downscaled_grc[g]=560.0; a2l.forc_rain_not_downscaled_grc[g]=mean(rain0); a2l.forc_snow_not_downscaled_grc[g]=0.0
    end
    C.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
    config = C.CLMDriverConfig(); fia = C.clump_filter_inactive_and_active
    (declin,_) = C.compute_orbital(152.0); nextsw = 152.0 + DT/C.SECSPDAY
    for n in 1:3
        C.clm_drv!(config, inst, filt, fia, bounds, true, nextsw, declin, declin,
            C.ORB_OBLIQR_DEFAULT, false, false, "", false; nstep=n, is_first_step=(n==1),
            is_beg_curr_day=(n==1), dtime=DT, mon=6, day=1, photosyns=inst.photosyns)
    end
    jz = C.varpar.nlevsno + 1
    for c in bounds.begc:bounds.endc
        filt.hydrologyc[c] || continue
        nb = inst.column.nbedrock[c]
        inst.soilhydrology.zwt_col[c] = max(0.3, inst.column.zi[c, nb+jz] - 1.0)
    end
    return inst, bounds, filt, config
end

inst, bounds, filt, config = build()
hc = filt.hydrologyc
const HCOLS = Int[c for c in bounds.begc:bounds.endc if hc[c]]
const FFF0 = FT(C.sat_excess_runoff_params.fff)
const BFS0 = FT(C.BASEFLOW_SCALAR[])
rain_ser = read_forcing_window(START, NSTEPS)
@printf("STAGE A ForwardDiff feasibility: domain=%s N=%d steps  fff0=%.4f bf0=%.3e\n", DOMAIN, NSTEPS, FFF0, BFS0)

# soil_temperature! phase whose Float64 aux (urbantv) is rebuilt to the BUNDLE eltype, so
# the urban building-temp sub-path (_BTLunIn) doesn't mix Float64 aux with Dual state.
function soiltemp_ft!(b, bnds, flt)
    T = eltype(b.theta)
    a0 = C.soiltemp_rev_aux(bnds, flt; dtime=DT)
    a = (; urbantv = fill(T(a0.urbantv[1]), length(a0.urbantv)), dtime=a0.dtime,
           bc_col=a0.bc_col, bc_lun=a0.bc_lun, bc_patch=a0.bc_patch,
           nolakec=a0.nolakec, nolakep=a0.nolakep, urbanl=a0.urbanl, urbanc=a0.urbanc)
    C.soiltemp_rev_phase!(b, a)
end
a_surf = C.surfhydro_rev_aux(bounds, filt; dtime=DT, use_hydrstress=config.use_hydrstress)
bf_aux = (; cols=HCOLS, nlevsno=C.varpar.nlevsno, nlevsoi=C.varpar.nlevsoi,
            n_baseflow=FT(C.soilhydrology_params.n_baseflow), dt=DT)
function step_phases(k)
    Any[(rain_inject!, ((; cols=HCOLS, rain=max(rain_ser[k], 3e-4)),)),
        (soiltemp_ft!, (bounds, filt)),   # soil_temperature! — probing the urban sub-path wall
        (C.setsoilfrac_rev_phase!, (a_surf,)), (C.setfloodc_rev_phase!, (a_surf,)),
        (satexcess_cal!, ((; cols=HCOLS),)),
        (C.setqflx_rev_phase!, (a_surf,)), (C.inflexcess_rev_phase!, (a_surf,)),
        (C.routeinfl_rev_phase!, (a_surf,)), (C.updateh2osfc_rev_phase!, (a_surf,)),
        (C.infil_rev_phase!, (a_surf,)), (C.totalrunoff_rev_phase!, (a_surf,)),
        (C.soilwater_rev_phase!, (C.soilwater_rev_aux(filt, config; dtime=DT),)),
        (C.watertable_rev_phase!, (C.watertable_rev_aux(bounds, filt, config; dtime=DT),)),
        (C.hydnodrain_rev_phase!, (C.hydnodrain_rev_aux(bounds, filt; dtime=DT),)),
        (baseflow_cal!, (bf_aux,)),
        (accum_runoff!, ((; cols=HCOLS, k=k),))]
end
const STEPS = [step_phases(k) for k in 1:NSTEPS]

# ---- ForwardDiff path: retype base inst to eltype(θ) INSIDE the differentiated fn ----
function series_T(θ::AbstractVector{T}) where {T}
    instT = retype(T, inst)
    b = (; inst = instT, scratch = C.cf_rev_scratch(T, length(instT.patch.column)),
           theta = collect(T, θ), series = zeros(T, NSTEPS))
    for ph in STEPS, (f,ca) in ph; f(b, ca...); end
    return b.series
end
loss_fd(θ) = sum(abs2, series_T(θ))              # simple SSE-to-zero-target (gate is grad vs FD)

# reverse-AD gradient of the SAME loss/window (Float64 inst + Enzyme multistep_reverse!) — the
# cross-check that ForwardDiff (forward-mode, Dual) agrees with the validated reverse-mode path.
function grad_rev(θ)
    mkbundle(θ) = (; C.driver_rev_bundle(deepcopy(inst))..., theta=collect(FT,θ), series=zeros(FT,NSTEPS))
    b0 = mkbundle(θ); for ph in STEPS, (f,ca) in ph; f(b0,ca...); end
    dd = 2 .* b0.series                          # dL/dseries for L = sum(abs2, series)
    seed!(db,b) = (db.series .= dd; nothing)
    return copy(C.multistep_reverse!(STEPS, mkbundle(θ), seed!).theta)
end

# =============================================================================
# TEST 1 — does retype build a Dual inst + run ONE step? (isolate the break site)
# =============================================================================
println("\n[T1] retype inst → Dual eltype and run 1 step ...")
t1 = try
    D = ForwardDiff.Dual{:probe}(1.0, 1.0, 0.0)
    instD = retype(typeof(D), inst)
    @printf("     retype OK: eltype(zwt_col)=%s  eltype(qflx_surf)=%s\n",
        string(eltype(instD.soilhydrology.zwt_col)),
        string(eltype(instD.water.waterfluxbulk_inst.wf.qflx_surf_col)))
    @printf("     structs left Float64 (config/no-ctor): %s\n", string(sort(collect(RETYPE_FALLBACK))))
    b = (; inst=instD, scratch=C.cf_rev_scratch(typeof(D), length(instD.patch.column)),
           theta=[D, D], series=zeros(typeof(D), NSTEPS))
    for (f,ca) in STEPS[1]; f(b, ca...); end
    @printf("     1-step run OK: series[1]=%s\n", string(b.series[1]))
    true
catch e
    println("     BREAK: ", sprint(showerror, e))
    for (i,fr) in enumerate(stacktrace(catch_backtrace())); i>8 && break; println("       ", fr); end
    false
end

# =============================================================================
# TEST 2 — ForwardDiff gradient vs central FD vs reverse-AD (the correctness gate)
# =============================================================================
if t1
    println("\n[T2] ForwardDiff.gradient vs central FD ...")
    θ0 = [FFF0*1.2, BFS0*1.3]
    try
        t_fwd = @elapsed g_fd = ForwardDiff.gradient(loss_fd, θ0)
        @printf("     ForwardDiff grad = [% .6e, % .6e]  (%.1fs)\n", g_fd[1], g_fd[2], t_fwd)
        for (j,name,h) in ((1,"fff",1e-3*FFF0),(2,"baseflow_scalar",1e-3*BFS0))
            θa=copy(θ0); θa[j]+=h; θb=copy(θ0); θb[j]-=h
            cfd=(loss_fd(θa)-loss_fd(θb))/(2h); rel=abs(g_fd[j]-cfd)/max(abs(cfd),1e-30)
            @printf("     ∂/∂%-16s FwdDiff=% .6e centralFD=% .6e rel=%.2e %s\n",
                name, g_fd[j], cfd, rel, rel<1e-6 ? "PASS ✓" : "CHECK")
        end
        # cross-check vs reverse-AD (non-fatal: Enzyme reverse may not cover every phase Fwd does)
        try
            gr = grad_rev(θ0)
            for (j,name) in ((1,"fff"),(2,"baseflow_scalar"))
                rel=abs(g_fd[j]-gr[j])/max(abs(gr[j]),1e-30)
                @printf("     FwdDiff-vs-Reverse %-16s fwd=% .6e rev=% .6e rel=%.2e %s\n",
                    name, g_fd[j], gr[j], rel, rel<1e-6 ? "MATCH ✓" : "DIFF")
            end
        catch e
            @printf("     reverse cross-check n/a (reverse-mode): %s\n", sprint(showerror,e)[1:min(60,end)])
        end
    catch e
        println("     ForwardDiff.gradient BREAK: ", sprint(showerror, e))
        for (i,fr) in enumerate(stacktrace(catch_backtrace())); i>8 && break; println("       ", fr); end
    end
end
println("\nSTAGE A verdict printed above (T1 retype+run, T2 ForwardDiff==FD).")
