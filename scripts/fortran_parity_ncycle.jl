# ==========================================================================
# N-cycle Fortran parity: deposition + fixation + leaching.
#
# Scores the nitrogen INPUT and LOSS terms — the routines that were dead ports
# (defined, never called) until this PR — against an instrumented CTSM run at
# Bow-at-Banff.
#
# THE FORTRAN REFERENCE (see docs/N_CYCLE_PARITY.md)
#   run dir : SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_ndep{,_summer}
#   windows : autumn  nstep 1760377..1760560 (2202-10-29 →, dormant/frozen)
#             summer  nstep 1757873..1758060 (2202-07-16 →, ACTIVE DRAINAGE)
#   dtime   : 3600 s
#   ndep    : SYNTHETIC stream, spatially uniform + time-constant,
#             NDEP_month = 1.0e-8 g(N)/m2/s (= 0.3154 gN/m2/yr, in the observed
#             range for the Canadian Rockies).
#             The real CMIP6 fndep_* file is not on disk and could not be
#             re-fetched (the server cert does not cover the inputdata host;
#             disabling TLS verification was declined), so BOTH codes are driven
#             with this same synthetic stream. Because it is uniform and
#             constant, ESMF's bilinear regrid and its linear time interpolation
#             are both exact, so Fortran's forc_ndep_grc == 1.0e-8 analytically.
#             That isolates the N-cycle PHYSICS wiring from forcing-file
#             interpolation fidelity, which is what we are actually testing.
#   config  : use_cn, use_fun=.true., use_nitrif_denitrif=.true., use_crop=.false.
#             => Fortran runs CNNDeposition + CNFreeLivingFixation (ffix), NOT
#                CNNFixation, and puts ndep+ffix into smin_nh4 (the FUN branch of
#                SoilBiogeochemNStateUpdate1).
#
# BOUNDARIES
#   before_step                -> injected as the shared IC
#   after_hydrologydrainage    -> after EcosystemDynamicsPreDrainage (ndep+fix in,
#                                 decomposition applied; BEFORE leaching)
#   after_ecosysdyn_postdrain  -> after EcosystemDynamicsPostDrainage (leaching
#                                 applied). NEW boundary added for this work.
#
# Because Fortran dumps BOTH sides of the leaching call, the leaching flux itself
# has an exact per-level Fortran ground truth:
#     (leached + runoff)[j] == (no3_pre[j] - no3_post[j]) / dt
# which is what §"LEACHING" below scores Julia against.
#
# Usage:
#   julia --project=. scripts/fortran_parity_ncycle.jl [summer|autumn] [nprobe]
# ==========================================================================

include(joinpath(@__DIR__, "fortran_parity_common.jl"))
using Printf

const NDEP_VAL = 1.0e-8            # g(N)/m2/s, by construction of the stream
const FNDEP = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/ndepdata/" *
              "fndep_clm_UNIFORM1e-8_0.9x1.25_yr2000_CLMjl-parity.nc"
const DT = 3600.0

const WINDOWS = Dict(
    # name   => (dumpdir, first nstep, last nstep, model date of first nstep)
    "autumn" => ("/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_ndep",
                 1760377, 1760560, DateTime(2003, 10, 29, 0)),
    "summer" => ("/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_ndep_summer",
                 1757873, 1758060, DateTime(2003, 7, 16, 16)),
)

_dv(ds, n) = haskey(ds, n) ? vec(Float64.(Array(ds[n][:]))) : nothing

"""
Worst relative error over a vertical profile, scaled by the PROFILE's own
magnitude (max|F| over the levels), not by each element.

Per-element scaling is wrong here: the deep soil levels hold mineral N that is
zero (or ~1e-20) in Fortran, so an absolutely negligible Julia value there
produces a meaningless 1e+14 "relative error" that swamps the real signal. What
we care about is the error relative to how much N is in the column.
"""
function relerr(j, f)
    jj = vec(Float64.(j)); ff = vec(Float64.(f))
    n = min(length(jj), length(ff))
    scale = 0.0
    for i in 1:n
        isfinite(ff[i]) && (scale = max(scale, abs(ff[i])))
    end
    scale <= 0 && (scale = 1.0)
    worst = 0.0
    for i in 1:n
        a, b = jj[i], ff[i]
        (isfinite(a) && isfinite(b)) || continue
        worst = max(worst, abs(a - b) / scale)
    end
    return worst
end

function score_step(win::String, nstep::Int; with_ndep::Bool = true)
    (dumpdir, n0, _, date0) = WINDOWS[win]

    inst, _ = run_one_parity_step!(nstep;
        use_cn = true, dumpdir = dumpdir,
        step_date = date0 + Hour(nstep - n0),
        use_hydrstress = true,
        fndep = with_ndep ? FNDEP : "")

    sns = inst.soilbiogeochem_nitrogenstate
    snf = inst.soilbiogeochem_nitrogenflux
    nlev = CLM.varpar.nlevdecomp

    pre  = NCDataset(joinpath(dumpdir, "pdump_after_hydrologydrainage_n$(nstep).nc"), "r")
    post = NCDataset(joinpath(dumpdir, "pdump_after_ecosysdyn_postdrain_n$(nstep).nc"), "r")
    bs   = NCDataset(joinpath(dumpdir, "pdump_before_step_n$(nstep).nc"), "r")

    r = Dict{String, Any}()

    # ---------------- N INPUT fluxes (analytic Fortran ground truth) ----------
    r["ndep_J"] = snf.ndep_to_sminn_col[1]
    r["ndep_F"] = with_ndep ? NDEP_VAL : 0.0

    # Fortran ffix (CNFreeLivingFixation), evaluated from the SAME AnnET the dump
    # carries. Julia computes its own from the injected AnnET, so this is a real
    # check of the formula, not a tautology.
    annet = _dv(bs, "AnnET_VALUE")[1]
    spy = 86400.0 * 365.0
    r["ffix_J"] = snf.ffix_to_sminn_col[1]
    r["ffix_F"] = (0.0006 * (max(0.0, annet) * spy) + 0.0117) / spy
    r["AnnET"]  = annet

    # ---------------- mineral-N POOLS ----------------------------------------
    # Julia's clm_drv! runs PostDrainage inside the step, so the end-of-step inst
    # state is the POST-leaching state -> compare against `post`.
    for v in ("smin_nh4_vr", "smin_no3_vr", "sminn_vr")
        jf = v == "smin_nh4_vr" ? sns.smin_nh4_vr_col :
             v == "smin_no3_vr" ? sns.smin_no3_vr_col : sns.sminn_vr_col
        f = _dv(post, v)
        f === nothing || (r["post_" * v] = relerr(vec(jf[1, 1:nlev]), f[1:nlev]))
    end

    # ---------------- LEACHING (exact per-level Fortran ground truth) ---------
    no3_pre  = _dv(pre,  "smin_no3_vr")[1:nlev]
    no3_post = _dv(post, "smin_no3_vr")[1:nlev]
    leach_F  = (no3_pre .- no3_post) ./ DT                      # gN/m3/s
    leach_J  = vec(snf.smin_no3_leached_vr_col[1, 1:nlev]) .+
               vec(snf.smin_no3_runoff_vr_col[1, 1:nlev])

    r["leach_F_tot"] = sum(leach_F)
    r["leach_J_tot"] = sum(leach_J)
    r["leach_active"] = sum(abs, leach_F) > 0
    r["leach_rel"] = relerr(leach_J, leach_F)

    close(pre); close(post); close(bs)
    return r
end

# --------------------------------------------------------------------------
function main()
    win    = length(ARGS) >= 1 ? ARGS[1] : "summer"
    nprobe = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 5
    haskey(WINDOWS, win) || error("unknown window '$win' (summer|autumn)")
    (_, n0, n1, _) = WINDOWS[win]
    steps = round.(Int, range(n0 + 2, n1 - 2; length = nprobe))

    println("="^100)
    println("N-CYCLE FORTRAN PARITY — Bow at Banff — window=$win")
    println("ndep = $NDEP_VAL g(N)/m2/s (SYNTHETIC uniform stream, same file drives both codes)")
    println("="^100)

    @printf("\n%-9s %11s %11s %11s %11s | %10s %10s %10s | %10s\n",
            "nstep", "ndep J", "ndep F", "ffix J", "ffix F",
            "nh4 post", "no3 post", "sminn post", "leach rel")
    println("-"^100)

    wq = Dict("post_smin_nh4_vr"=>0.0, "post_smin_no3_vr"=>0.0, "post_sminn_vr"=>0.0)
    ndep_ok = true; ffix_ok = true
    leach_worst = 0.0; leach_n = 0

    for ns in steps
        r = score_step(win, ns)
        ndep_ok &= abs(r["ndep_J"] - r["ndep_F"]) <= 1e-12 * max(r["ndep_F"], 1e-30)
        ffix_ok &= abs(r["ffix_J"] - r["ffix_F"]) <= 1e-12 * max(r["ffix_F"], 1e-30)
        for k in keys(wq); wq[k] = max(wq[k], get(r, k, 0.0)); end
        if r["leach_active"]
            leach_worst = max(leach_worst, r["leach_rel"]); leach_n += 1
        end
        @printf("%-9d %11.4e %11.4e %11.4e %11.4e | %10.2e %10.2e %10.2e | %10.2e%s\n",
                ns, r["ndep_J"], r["ndep_F"], r["ffix_J"], r["ffix_F"],
                get(r, "post_smin_nh4_vr", NaN), get(r, "post_smin_no3_vr", NaN),
                get(r, "post_sminn_vr", NaN), r["leach_rel"],
                r["leach_active"] ? "" : "  (no drainage)")
    end

    println("-"^100)
    println("N-INPUT FLUXES vs Fortran (analytic ground truth):")
    println("  ndep_to_sminn : ", ndep_ok ? "EXACT (rel <= 1e-12)" : "*** MISMATCH ***")
    println("  ffix_to_sminn : ", ffix_ok ? "EXACT (rel <= 1e-12)" : "*** MISMATCH ***")
    println("\nN LOSS (leaching+runoff) vs Fortran (from the pre/post dump pair):")
    if leach_n == 0
        println("  NOT EXERCISED in this window — Fortran's own no3 pre->post change is")
        println("  identically zero (frozen/dry column, no drainage). Use the summer window.")
    else
        @printf("  worst |rel| over %d probes with active drainage: %.3e\n", leach_n, leach_worst)
    end
    println("\nMINERAL-N POOL worst |rel| (vs after_ecosysdyn_postdrain):")
    for k in sort(collect(keys(wq))); @printf("  %-20s %.3e\n", k, wq[k]); end

    # ------------- directional test: does wiring N move Julia TOWARD Fortran? --
    println("\n" * "="^100)
    println("DIRECTIONAL TEST — mineral-N error WITH the N inputs wired vs WITHOUT (= pre-PR behaviour)")
    println("="^100)
    @printf("%-9s %14s %14s %14s %14s %10s\n", "nstep",
            "sminn ON", "sminn OFF", "nh4 ON", "nh4 OFF", "improved?")
    println("-"^100)
    for ns in steps[1:min(3, length(steps))]
        on  = score_step(win, ns; with_ndep = true)
        off = score_step(win, ns; with_ndep = false)
        a = get(on, "post_sminn_vr", NaN);  b = get(off, "post_sminn_vr", NaN)
        c = get(on, "post_smin_nh4_vr", NaN); d = get(off, "post_smin_nh4_vr", NaN)
        @printf("%-9d %14.4e %14.4e %14.4e %14.4e %10s\n", ns, a, b, c, d,
                (a < b && c < d) ? "YES" : "no")
    end
end

main()
