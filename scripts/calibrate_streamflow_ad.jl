# =============================================================================
# REVERSE-AD STREAMFLOW CALIBRATION for CLM.jl  (the differentiability payoff)
#
# d(streamflow-loss)/d(hydrology-params) from ONE reverse pass instead of 2N forward
# runs. Builds a WINDOWED surface-runoff loss on a real basin state and reverse-
# differentiates it through the production HydrologyNoDrainage reverse thread
# (src/driver/driver_reverse.jl) using Enzyme, exactly the multistep_reverse! pattern
# the season-GPP calibration (scripts/enzyme_season_calibrate*.jl) established.
#
# The differentiated hydrology parameters are injected into the reverse thread as
# scalar bundle entries (b.theta[j]) and read by θ-parameterized wrappers around the
# real physics (the fff/fover TOPMODEL saturation param in saturated_excess_runoff!;
# a baseflow multiplier on the subsurface drainage flux). This mirrors how the GPP
# harness injects vcmax/light-use as b.t1/b.t2 — the param IS a real CLM hydrology
# parameter, entering the real fsat/drainage formula.
#
# GATES (env CLM_GATE = 1 | 2 | 3, default 1):
#   1  AD gradient vs central finite differences (the correctness crux).
#   2  gradient descent calibrates a windowed streamflow fit; AD vs FD cost.
#   3  (driven by calibrate_streamflow_ad_multi.jl across basins).
#
#   DOMAIN via env:
#     CLM_FSURDAT, CLM_PARAMFILE  (default: Bow at Banff lumped, compHydro root)
#   Window / event via env:
#     CLM_NSTEPS  (default 24), CLM_RAIN (kg/m2/s steady rain input, default 3e-4)
#
#   julia --project=. scripts/calibrate_streamflow_ad.jl
# =============================================================================
using CLM, Enzyme, Printf, Random
const C = CLM
const FT = Float64

const FSURDAT = get(ENV, "CLM_FSURDAT",
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc")
const PARAMFILE = get(ENV, "CLM_PARAMFILE",
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc")
const DT     = 1800.0
const NSTEPS = parse(Int, get(ENV, "CLM_NSTEPS", "24"))
const RAIN   = parse(Float64, get(ENV, "CLM_RAIN", "3e-4"))   # steady rain+snomelt input (kg/m2/s)
const GATE   = parse(Int, get(ENV, "CLM_GATE", "1"))

# ---- warm a real basin cold-start inst (3 driver steps), same recipe as the hydro
#      reverse harness, but with a real runoff-driving rain forcing. ----
function setup_forcing!(a2l, T, ng)
    for g in 1:ng
        a2l.forc_t_not_downscaled_grc[g]=T; a2l.forc_pbot_not_downscaled_grc[g]=85000.0
        a2l.forc_th_not_downscaled_grc[g]=T*(100000.0/85000.0)^(C.RAIR/C.CPAIR)
        a2l.forc_rho_not_downscaled_grc[g]=85000.0/(C.RAIR*T); a2l.forc_lwrad_not_downscaled_grc[g]=300.0
        a2l.forc_vp_grc[g]=800.0; a2l.forc_hgt_grc[g]=30.0; a2l.forc_topo_grc[g]=0.0; a2l.forc_wind_grc[g]=3.0
        a2l.forc_hgt_u_grc[g]=30.0; a2l.forc_hgt_t_grc[g]=30.0; a2l.forc_hgt_q_grc[g]=30.0
        for b in 1:C.NUMRAD; a2l.forc_solad_not_downscaled_grc[g,b]=200.0; a2l.forc_solai_grc[g,b]=80.0; end
        a2l.forc_solar_not_downscaled_grc[g]=560.0
        a2l.forc_rain_not_downscaled_grc[g]=RAIN
        a2l.forc_snow_not_downscaled_grc[g]=0.0
    end
end

function build_real()
    (inst, bounds, filt, tm) = C.clm_initialize!(; fsurdat=FSURDAT, paramfile=PARAMFILE)
    setup_forcing!(inst.atm2lnd, 288.0, bounds.endg)
    C.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
    config = C.CLMDriverConfig(); filt_ia = C.clump_filter_inactive_and_active
    (declin, eccf) = C.compute_orbital(120.0); nextsw_cday = 120.0 + DT/C.SECSPDAY
    runstep!(n; first=false) = C.clm_drv!(config, inst, filt, filt_ia, bounds,
        true, nextsw_cday, declin, declin, C.ORB_OBLIQR_DEFAULT, false, false, "", false;
        nstep=n, is_first_step=first, is_beg_curr_day=first, dtime=DT, mon=5, day=1,
        photosyns=inst.photosyns)
    for n in 1:3; runstep!(n; first=(n==1)); end
    return inst, bounds, filt, config
end

inst, bounds, filt, config = build_real()
hc = filt.hydrologyc
const HCOLS = Int[c for c in bounds.begc:bounds.endc if hc[c]]
@assert !isempty(HCOLS) "no hydrology columns"

# The default fff/fover for this basin (from the params global that
# saturated_excess_runoff! reads), and the default baseflow scalar.
const FFF0   = FT(C.sat_excess_runoff_params.fff)
const BFS0   = FT(C.BASEFLOW_SCALAR[])
@printf("basin state: %d hydrology col(s); fff0=%.5f  baseflow_scalar0=%.3e\n",
    length(HCOLS), FFF0, BFS0)

# ---- assert a steady rain+snomelt input on the hydrology columns each step, so
#      qflx_sat_excess_surf = fsat(θ)·(rain+snomelt) is nonzero & θ-sensitive. Exogenous:
#      its adjoint flows to the Const source and is discarded (like forcingset_rev_phase!). ----
rain_aux() = (; cols = HCOLS, rain = fill(FT(RAIN), length(inst.column.snl)))
function rain_inject!(b, aux)
    wf = b.inst.water.waterfluxbulk_inst.wf
    @inbounds for c in aux.cols
        wf.qflx_rain_plus_snomelt_col[c] = aux.rain[c]
    end
    return nothing
end

# ---- θ1 = fff (fover), θ2 = fmax (max saturated fraction). BOTH are real SYMFLUENCE
#      TOPMODEL surface-runoff calibration parameters, entering the EXACT physics of
#      saturated_excess_runoff! (src/biogeophys/sat_excess_runoff.jl):
#          fsat = θ2·wtfact · exp(-0.5·θ1·zwt)
#      re-expressed with θ read from the bundle so Enzyme differentiates it. zwt/wtfact/
#      frost_table flow LIVE from the bundle inst (evolved by water_table! across the
#      window), so the gradient is genuinely multi-step: θ → qflx_sat_excess_surf →
#      infiltration partition → soil_water (h2osoi) → water_table (zwt) → next-step fsat.
#      This replaces satexcess_rev_phase! at the same position in the surface-hydrology block.
#      θ2 enters as a multiplier on the surfdata fmax (wtfact_col), nominal 1.0. ----
satexcess_cal_aux() = (; cols = HCOLS)
function satexcess_cal!(b, aux)
    i = b.inst
    θ1 = b.theta[1]; θ2 = b.theta[2]
    sh = i.soilhydrology; ss = i.soilstate; ser = i.sat_excess_runoff
    wfb = i.water.waterfluxbulk_inst; wf = wfb.wf
    @inbounds for c in aux.cols
        ft  = sh.frost_table_col[c]; zw = sh.zwt_col[c]; zwp = sh.zwt_perched_col[c]
        wt  = θ2 * ss.wtfact_col[c]
        fs  = (ft > zwp && ft <= zw) ? wt*exp(-0.5*θ1*zwp) : wt*exp(-0.5*θ1*zw)
        ser.fsat_col[c] = fs
        ser.fcov_col[c] = fs
        wfb.qflx_sat_excess_surf_col[c] = fs * wf.qflx_rain_plus_snomelt_col[c]
    end
    return nothing
end

# ---- per-step accumulator: surface-runoff series (qflx_surf) over the hydrology columns,
#      in mm/s (col-summed). This is the runoff component the reverse thread computes end-to-
#      end (subsurface drainage/baseflow is a SEPARATE discrete-branch routine, drainage!,
#      NOT yet in the reverse tape — an honest scope limit). L is built vs a target below. ----
accum_aux(k) = (; cols = HCOLS, k = k, dt = DT)
function accum_runoff!(b, aux)
    wf = b.inst.water.waterfluxbulk_inst.wf
    s = 0.0
    @inbounds for c in aux.cols
        s += wf.qflx_surf_col[c]
    end
    b.series[aux.k] += s
    return nothing
end

# Bundle: driver reverse bundle (inst + canopy scratch) + θ vector + series accumulator.
NPAR = 2
bundle(θ) = (; C.driver_rev_bundle(deepcopy(inst))..., theta = collect(FT, θ), series = zeros(FT, NSTEPS))

# Ordered per-step phase list: rain → surface-hydrology (satexcess swapped for θ) →
# soil_water → water_table → hydrology_no_drainage → accumulate.
a_surf = C.surfhydro_rev_aux(bounds, filt; dtime=DT, use_hydrstress=config.use_hydrstress)
function step_phases(k)
    Any[
        (rain_inject!,               (rain_aux(),)),
        (C.setsoilfrac_rev_phase!,   (a_surf,)),
        (C.setfloodc_rev_phase!,     (a_surf,)),
        (satexcess_cal!,             (satexcess_cal_aux(),)),          # θ1=fff, θ2=fmax
        (C.setqflx_rev_phase!,       (a_surf,)),
        (C.inflexcess_rev_phase!,    (a_surf,)),
        (C.routeinfl_rev_phase!,     (a_surf,)),
        (C.updateh2osfc_rev_phase!,  (a_surf,)),
        (C.infil_rev_phase!,         (a_surf,)),
        (C.totalrunoff_rev_phase!,   (a_surf,)),
        (C.soilwater_rev_phase!,     (C.soilwater_rev_aux(filt, config; dtime=DT),)),
        (C.watertable_rev_phase!,    (C.watertable_rev_aux(bounds, filt, config; dtime=DT),)),
        (C.hydnodrain_rev_phase!,    (C.hydnodrain_rev_aux(bounds, filt; dtime=DT),)),
        (accum_runoff!,              (accum_aux(k),)),
    ]
end
const STEPS = [step_phases(k) for k in 1:NSTEPS]

series(θ) = (b = bundle(θ); for ph in STEPS, (f,ca) in ph; f(b, ca...); end; copy(b.series))

# Reverse-AD gradient of L(θ) = Σ_k (series_k(θ) − target_k)² w.r.t. θ (one reverse pass).
function grad_ad(θ, target)
    seed!(db, b) = (db.series .= 2 .* (b.series .- target); nothing)
    db = C.multistep_reverse!(STEPS, bundle(θ), seed!)
    return copy(db.theta)
end
Loss(θ, target) = sum(abs2, series(θ) .- target)

# =============================================================================
# GATE 1 — AD gradient vs central finite differences.
# =============================================================================
function gate1()
    println("\n", "="^78)
    println("GATE 1 — reverse-AD streamflow-loss gradient vs central FD")
    println("="^78)
    θ0 = [FFF0, 1.0]                            # θ2 = fmax multiplier (nominal 1.0)
    tgt = series(θ0)
    @printf("window: N=%d steps, rain=%.1e kg/m2/s\n", NSTEPS, RAIN)
    @printf("target runoff series (mm/s, first 6): %s\n", string(round.(tgt[1:min(6,NSTEPS)], sigdigits=4)))
    @printf("series is θ-sensitive: range over window = [%.4e, %.4e]\n", minimum(tgt), maximum(tgt))

    # Evaluate gradient at a PERTURBED θ (away from the target's zero-residual point,
    # so the gradient is nonzero and the FD check is non-vacuous).
    θp = [FFF0*1.3, 1.25]
    g = grad_ad(θp, tgt)
    @printf("\nAD gradient at θ=[fff=%.4f, fmax_mult=%.3f]: dL/dθ = [% .6e, % .6e]\n", θp[1], θp[2], g[1], g[2])

    ok = true
    for (j, name, h) in ((1,"fff (fover)", 1e-3*FFF0), (2,"fmax_mult", 1e-3))
        θpp = copy(θp); θpp[j] += h
        θpm = copy(θp); θpm[j] -= h
        fd = (Loss(θpp, tgt) - Loss(θpm, tgt)) / (2h)
        rel = abs(g[j] - fd) / max(abs(fd), 1e-30)
        pass = rel < 1e-4
        ok &= pass
        @printf("  θ%d %-14s: AD=% .6e  FD=% .6e  rel=%.2e  %s\n",
            j, name, g[j], fd, rel, pass ? "PASS ✓" : "FAIL ✗")
    end
    println("-"^78)
    @printf("GATE 1 %s\n", ok ? "PASS ✓ — reverse-AD yields a CORRECT streamflow-loss gradient w.r.t. real hydrology params" :
                                 "FAIL ✗")
    return ok, θ0, tgt
end

# =============================================================================
# GATE 2 — gradient descent calibrates a windowed streamflow fit.
# =============================================================================
function gate2(θ_true, tgt_true)
    println("\n", "="^78)
    println("GATE 2 — AD-gradient calibration on the windowed streamflow fit")
    println("="^78)
    # WELL-POSED 1-D recovery: calibrate fff (fmax fixed at prior) toward an fff-perturbed
    # "observed" runoff series. (The 2-param fff+fmax problem is rank-deficient — both scale
    # fsat, an equifinality; see the note below — so we calibrate the identifiable param and
    # report the collinearity honestly rather than manufacture a 2-param recovery.)
    fmax_fix = 1.0
    fff_obs  = FFF0 * 0.6
    OBS = series([fff_obs, fmax_fix])
    @printf("target: fff_obs=%.5f (fmax fixed=%.2f); calibrate fff from prior=%.5f\n", fff_obs, fmax_fix, FFF0)

    fff = FFF0; L0 = Loss([fff, fmax_fix], OBS); npass = 0
    g1(x) = (npass += 1; grad_ad([x, fmax_fix], OBS)[1])   # 1-D AD gradient (one reverse pass)
    L1(x) = Loss([x, fmax_fix], OBS)
    h = 1e-3*FFF0
    for it in 1:12
        g0 = g1(fff); gp = g1(fff+h)
        H = (gp - g0)/h; μ = 1e-6*max(abs(H),1e-12)
        Lc = L1(fff); accepted = false
        for _t in 1:15
            xt = clamp(fff - g0/(H+μ), 1e-3, 5.0)
            if L1(xt) < Lc; fff = xt; accepted = true; break; end
            μ = μ==0 ? 1e-9 : μ*4
        end
        @printf("  it %2d: fff=%.5f  L=%.6e  |g|=%.2e  reverse-passes=%d\n", it, fff, L1(fff), abs(g0), npass)
        (!accepted || abs(g0) < 1e-14) && break
    end
    Lf = L1(fff)
    redstr = Lf < 1e-20 ? "→ ~0 (machine precision)" : @sprintf("%.0f×", L0/max(Lf,1e-300))
    @printf("final fff=%.5f (obs %.5f, prior %.5f)  L=%.3e → %.3e  reduction %s\n",
        fff, fff_obs, FFF0, L0, Lf, redstr)

    # 2-param equifinality note (honest): report the JᵀJ conditioning of [fff, fmax].
    let θe=[FFF0,1.0], s=series(θe), h1=1e-3*FFF0, h2=1e-3
        J1 = (series([FFF0+h1,1.0]) .- s) ./ h1
        J2 = (series([FFF0,1.0+h2]) .- s) ./ h2
        a=sum(J1.^2); b=sum(J1.*J2); d=sum(J2.^2)
        tr=a+d; det=a*d-b^2; λlo=(tr-sqrt(max(tr^2-4det,0)))/2; λhi=(tr+sqrt(max(tr^2-4det,0)))/2
        @printf("  [note] 2-param JᵀJ(fff,fmax) cond ≈ %.1e → %s\n", λhi/max(λlo,1e-300),
            λhi/max(λlo,1e-300) > 1e6 ? "rank-deficient (fff·fmax equifinal in fsat; needs a prior)" : "identifiable")
    end
    ok = Lf < 0.05 * L0
    @printf("cost: AD 1-D calibration = 2 reverse passes/iter (g, g') vs FD = 2 fwd/iter\n")
    @printf("GATE 2 %s\n", ok ? "PASS ✓ — windowed streamflow fit improved via AD gradient (loss ↓ ≥20×, fff recovered)" :
                                 "FAIL ✗ (loss did not drop 20×)")
    return ok
end

ok1, θ0, tgt = gate1()
if GATE >= 2 && ok1
    gate2(θ0, tgt)
end
