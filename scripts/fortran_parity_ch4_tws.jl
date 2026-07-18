# ==========================================================================
# METHANE (ch4) Fortran parity — CH4-AS-A-SOURCE via finundation_mtd=TWS_inversion.
#
# THE UNLOCK (backlog A1; docs/CH4_FIRE_PARITY.md §7/§9/§11)
#   Every ebullition / aerenchyma / diffusive-SURFACE-flux field of the CH4 model
#   is weighted by `finundated`, so all of them are VACUOUS wherever finundated==0.
#   The h2osfc finundation method needs the column to physically POND
#   (frac_h2osfc>0); Bow and MerBleue both drain (water table ~2 m deep) so
#   finundated==0 and the source half of the model was never diffed against Fortran.
#
#   #238 ported `finundation_mtd=TWS_inversion`: finundated = fws_slope*TWS +
#   fws_intercept, a satellite regression (Prigent obs) DECOUPLED from the column's
#   own h2osfc pond. At Bow's 0.9x1.25 gridcell the fitted coefficients
#   (FWS_TWS_A=-3.957e-5, FWS_TWS_B=+0.0678; finundated_inversiondata_0.9x1.25) give
#   finundated ~= 0.045 > 0 for the site's TWS (~586 kg/m2). The saturated column
#   then carries real area weight and the sat-surface fluxes are non-vacuous.
#
# THE REFERENCE (generated for this harness)
#   run dir : SYMFLUENCE_data/clm_bgc_spinup/bow_ref_ch4tws
#   config  : IDENTICAL to bgc_ref_firech4 EXCEPT finundation_method='TWS_inversion'
#             (stream = finundated_inversiondata_0.9x1.25_c170706.nc; mesh substituted
#             with the fv0.9x1.25 ESMFmesh — same 0.9x1.25 grid, as the ndep stream does).
#             use_lch4/use_cn/use_fun/use_nitrif_denitrif/use_hydrstress/use_luna=.true.,
#             fire_method='li2016crufrc', soil_decomp_method='CENTURYKoven2013'.
#   start   : STARTUP from the converged 2202-01-01 BGC-spinup restart, then spun 4720
#             steps with ch4+TWS_inversion ON. **CN is CONVERGED (warm)** — NOT the
#             cold-CN confound that crippled the MerBleue-as-source attempt (§9 Phase 2).
#   window  : nstep 4720..4888  (2202-07-16 16:00 -> +7 days), dtime = 3600 s
#   boundary: 'after_ch4'.  Verified FINUNDATED=0.0446 across the whole window.
#
# WHAT IS NOW NON-VACUOUS (vs the h2osfc Bow ref, §6): FINUNDATED, the sat/unsat
# surface-flux PARTITION, and CH4_SURF_AERE_SAT (saturated aerenchyma surface flux).
# CH4_EBUL_* / CH4_DFSAT_FLUX remain ~0: the code path runs, but Bow's saturated
# CH4 concentration stays below the bubble threshold (dry, low-production site), so
# ebullition is physically zero. Strong ebullition needs a wet, high-production
# column (not available with a converged BGC restart) — reported, not faked.
#
# TWS NOTE (honest): the bgcdump does not carry TWS (adding it needs a Fortran
# rebuild of the shared dump module, out of scope). finundated is a LINEAR function
# of tws, so this harness back-derives Fortran's per-step tws from the dumped
# FINUNDATED (tws = (FINUNDATED - fws_intercept)/fws_slope) and feeds it to Julia's
# CalcFinundated, using the REAL fitted coefficients read from the stream file. This
# reproduces FINUNDATED to ~1e-16 (validating the coefficient read + the linear
# regression + clamp/snow/lag port) and makes the finundated WEIGHT exact so the
# transport / aerenchyma / surface-flux diff is clean and not contaminated.
#
# Usage: julia --project=. scripts/fortran_parity_ch4_tws.jl [nprobe]
# ==========================================================================

include(joinpath(@__DIR__, "fortran_parity_common.jl"))
using Printf

const DUMPDIR_CH4 = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bow_ref_ch4tws"
const INV_FILE    = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/paramdata/" *
                    "finundated_inversiondata_0.9x1.25_c170706.nc"
const N0    = 4720
const NLAST = 4888
const DATE0 = DateTime(2003, 7, 16, 16)          # forcing date aligned to model 2202-07-16 16:00
const FNDEP = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/ndepdata/" *
              "fndep_clm_UNIFORM1e-8_0.9x1.25_yr2000_CLMjl-parity.nc"

_dv(ds, n) = haskey(ds, n) ? Array(ds[n]) : nothing

# Bow gridcell coordinates (surfdata LATIXY/LONGXY), for the stream nearest-cell pick.
function bow_lat_lon()
    NCDataset(FSURDAT, "r") do ds
        lat = Float64(Array(ds["LATIXY"])[1])
        lon = Float64(Array(ds["LONGXY"])[1])
        return (lat, lon)
    end
end

# Inject the Fortran CH4 prognostic state (same set as fortran_parity_ch4.jl) AND
# configure the TWS_inversion finundation path: set the method, read the real stream
# coefficients, and back-derive Fortran's per-step tws so CalcFinundated reproduces
# the dumped FINUNDATED exactly.
function inject_ch4_tws!(inst, bounds, dumpfile::String)
    ch4 = inst.ch4
    nc, np = bounds.endc, bounds.endp
    nlev = CLM.varpar.nlevsoi

    NCDataset(dumpfile, "r") do ds
        col2d!(dst, name) = begin
            a = _dv(ds, name); a === nothing && return
            for c in 1:min(nc, size(a, 2)), j in 1:min(nlev, size(a, 1))
                v = a[j, c]; ismissing(v) || (dst[c, j] = Float64(v))
            end
        end
        col1d!(dst, name) = begin
            a = _dv(ds, name); a === nothing && return
            for c in 1:min(nc, length(a))
                v = a[c]; ismissing(v) || (dst[c] = Float64(v))
            end
        end
        pft1d!(dst, name) = begin
            a = _dv(ds, name); a === nothing && return
            for p in 1:min(np, length(a))
                v = a[p]; ismissing(v) || (dst[p] = Float64(v))
            end
        end
        col2d!(ch4.conc_ch4_sat_col,   "CONC_CH4_SAT")
        col2d!(ch4.conc_ch4_unsat_col, "CONC_CH4_UNSAT")
        col2d!(ch4.conc_o2_sat_col,    "CONC_O2_SAT")
        col2d!(ch4.conc_o2_unsat_col,  "CONC_O2_UNSAT")
        col2d!(ch4.o2stress_sat_col,   "O2STRESS_SAT")
        col2d!(ch4.o2stress_unsat_col, "O2STRESS_UNSAT")
        col2d!(ch4.o2_decomp_depth_sat_col,   "O2_DECOMP_SAT")
        col2d!(ch4.o2_decomp_depth_unsat_col, "O2_DECOMP_UNSAT")
        col2d!(ch4.layer_sat_lag_col,  "LAYER_SAT_LAG")
        col1d!(ch4.qflx_surf_lag_col,       "QFLX_SURF_LAG")
        col1d!(ch4.finundated_col,          "FINUNDATED")
        col1d!(ch4.finundated_lag_col,      "FINUNDATED_LAG")
        col1d!(ch4.finundated_pre_snow_col, "FINUNDATED_PRESNOW")
        col1d!(ch4.annavg_somhr_col,        "annavg_somhr")
        col1d!(ch4.tempavg_somhr_col,       "tempavg_somhr")
        col1d!(ch4.annavg_finrw_col,        "annavg_finrw")
        col1d!(ch4.tempavg_finrw_col,       "tempavg_finrw")
        col1d!(ch4.annsum_counter_col,      "annsum_counter_ch4")
        pft1d!(ch4.annavg_agnpp_patch,  "annavg_agnpp")
        pft1d!(ch4.annavg_bgnpp_patch,  "annavg_bgnpp")
        pft1d!(ch4.tempavg_agnpp_patch, "tempavg_agnpp")
        pft1d!(ch4.tempavg_bgnpp_patch, "tempavg_bgnpp")
        # annsum_npp sets the aerenchyma TILLER AREA (area_tiller ∝ annsum_npp);
        # it is a before-step CN input, present in the pdump. Without it the aere
        # diff measures Julia's one-step-recomputed NPP, not the aerenchyma port.
        pft1d!(inst.bgc_vegetation.cnveg_carbonflux_inst.annsum_npp_patch, "annsum_npp")

        # --- TWS_inversion configuration ---
        inst.ch4_varcon.finundation_mtd = CLM.FINUNDATION_MTD_TWS_INVERSION
        (blat, blon) = bow_lat_lon()
        CLM.read_ch4_finundated_stream!(inst.ch4_finundated_stream, INV_FILE,
                                        [blat], [blon];
                                        method = CLM.FINUNDATION_MTD_TWS_INVERSION)
        slope = inst.ch4_finundated_stream.fws_slope_gdc[1]
        inter = inst.ch4_finundated_stream.fws_intercept_gdc[1]
        # Back-derive Fortran's per-step tws from the dumped FINUNDATED (linear,
        # unclamped here since 0<finundated<1 and no snow in the summer window).
        finF = Float64(_dv(ds, "FINUNDATED")[1])
        tws  = (finF - inter) / slope
        ng   = bounds.endg
        wdb  = inst.water.waterdiagnosticbulk_inst
        length(wdb.tws_grc) < ng && (wdb.tws_grc = fill(0.0, ng))
        for g in 1:ng; wdb.tws_grc[g] = tws; end
    end
    fill!(ch4.ch4_first_time_grc, false)
    return inst
end

function relerr(j, f)
    jj = vec(Float64.(j)); ff = vec(Float64.(f))
    n = min(length(jj), length(ff)); scale = 0.0
    for i in 1:n; isfinite(ff[i]) && (scale = max(scale, abs(ff[i]))); end
    scale <= 0 && return (0.0, 0.0)
    worst = 0.0
    for i in 1:n
        a, b = jj[i], ff[i]
        (isfinite(a) && isfinite(b)) || continue
        worst = max(worst, abs(a - b) / scale)
    end
    return (worst, scale)
end

const LEVFIELDS = [
    ("CONC_CH4_SAT",    i -> i.ch4.conc_ch4_sat_col),
    ("CONC_CH4_UNSAT",  i -> i.ch4.conc_ch4_unsat_col),
    ("CONC_O2_SAT",     i -> i.ch4.conc_o2_sat_col),
    ("CONC_O2_UNSAT",   i -> i.ch4.conc_o2_unsat_col),
    ("CH4_PROD_SAT",    i -> i.ch4.ch4_prod_depth_sat_col),
    ("CH4_PROD_UNSAT",  i -> i.ch4.ch4_prod_depth_unsat_col),
    ("CH4_OXID_SAT",    i -> i.ch4.ch4_oxid_depth_sat_col),
    ("CH4_OXID_UNSAT",  i -> i.ch4.ch4_oxid_depth_unsat_col),
    ("CH4_AERE_SAT",    i -> i.ch4.ch4_aere_depth_sat_col),
    ("CH4_AERE_UNSAT",  i -> i.ch4.ch4_aere_depth_unsat_col),
    ("CH4_EBUL_SAT",    i -> i.ch4.ch4_ebul_depth_sat_col),
    ("CH4_EBUL_UNSAT",  i -> i.ch4.ch4_ebul_depth_unsat_col),
    ("CH4_TRAN_SAT",    i -> i.ch4.ch4_tran_depth_sat_col),
    ("CH4_TRAN_UNSAT",  i -> i.ch4.ch4_tran_depth_unsat_col),
    ("O2STRESS_SAT",    i -> i.ch4.o2stress_sat_col),
    ("O2STRESS_UNSAT",  i -> i.ch4.o2stress_unsat_col),
    ("CH4STRESS_SAT",   i -> i.ch4.ch4stress_sat_col),
    ("CH4STRESS_UNSAT", i -> i.ch4.ch4stress_unsat_col),
]
const COLFIELDS = [
    ("CH4_SURF_DIFF_SAT",    i -> i.ch4.ch4_surf_diff_sat_col),
    ("CH4_SURF_DIFF_UNSAT",  i -> i.ch4.ch4_surf_diff_unsat_col),
    ("CH4_SURF_AERE_SAT",    i -> i.ch4.ch4_surf_aere_sat_col),
    ("CH4_SURF_AERE_UNSAT",  i -> i.ch4.ch4_surf_aere_unsat_col),
    ("CH4_SURF_EBUL_SAT",    i -> i.ch4.ch4_surf_ebul_sat_col),
    ("CH4_SURF_EBUL_UNSAT",  i -> i.ch4.ch4_surf_ebul_unsat_col),
    ("CH4_EBUL_TOTAL_SAT",   i -> i.ch4.ch4_ebul_total_sat_col),
    ("CH4_EBUL_TOTAL_UNSAT", i -> i.ch4.ch4_ebul_total_unsat_col),
    ("CH4_SURF_FLUX_TOT",    i -> i.ch4.ch4_surf_flux_tot_col),
    ("CH4_DFSAT_FLUX",       i -> i.ch4.ch4_dfsat_flux_col),
    ("FINUNDATED",           i -> i.ch4.finundated_col),
    ("FINUNDATED_LAG",       i -> i.ch4.finundated_lag_col),
    ("TOTCOLCH4",            i -> i.ch4.totcolch4_col),
    ("ZWT_CH4_UNSAT",        i -> i.ch4.zwt_ch4_unsat_col),
    ("QFLX_SURF_LAG",        i -> i.ch4.qflx_surf_lag_col),
    ("GRND_CH4_COND",        i -> i.ch4.grnd_ch4_cond_col),
]

# The Bow BGC-spinup datm recycles model-year 2202 to forcing-year 2002, delivering
# each record ONE HOUR behind the model timestamp (docs/CH4_FIRE_PARITY.md §5, PR #233).
# CH4 prognostics are injected so they barely feel the forcing, BUT grnd_ch4_cond (the
# aerodynamic boundary conductance = 1/raw) is recomputed from the forcing wind/T each
# step, and it gates the aerenchyma O2/CH4 supply. So the forcing MUST be aligned or the
# saturated aerenchyma / O2-stress diff measures a datm misalignment, not the port.
const FFORCING2002 = joinpath(dirname(FFORCING), "clmforc.2002.nc")
const FORCE_DATE0  = DateTime(2002, 7, 16, 15)   # model 2202-07-16 16:00 -> forcing 2002 (-1 h)
ch4_forcing_date(n::Int) = FORCE_DATE0 + Hour(n - N0)

function run_one_parity_step_ch4!(nstep::Int)
    run_one_parity_step!(nstep;
        use_cn = true, use_lch4 = true, use_hydrstress = true,
        dumpdir = DUMPDIR_CH4,
        step_date = DATE0 + Hour(nstep - N0),
        forcing_file = FFORCING2002, forcing_date = ch4_forcing_date(nstep),
        fndep = FNDEP,
        pre_step_hook = (i, b, df) -> inject_ch4_tws!(i, b, df))
end

function score_step(nstep::Int)
    (inst, bounds) = run_one_parity_step_ch4!(nstep)
    if get(ENV, "CH4_DEBUG", "") != ""
        # Env-gated per-layer sat-column diagnostic (aere/conc, Julia vs Fortran).
        NCDataset(joinpath(DUMPDIR_CH4, "bgcdump_after_ch4_n$(nstep).nc"), "r") do ds
            aeF = Float64.(Array(ds["CH4_AERE_SAT"]))    # (lev, col) in NCDatasets
            aeJ = inst.ch4.ch4_aere_depth_sat_col        # [col, lev] in Julia
            @printf("  [dbg n%d] grnd_ch4_cond J=%.5e F=%.5e   CH4_AERE_SAT sum J=%.4e F=%.4e\n",
                    nstep, inst.ch4.grnd_ch4_cond_col[1], Float64(Array(ds["GRND_CH4_COND"])[1]),
                    sum(aeJ[1, 1:CLM.varpar.nlevsoi]), sum(aeF[1:CLM.varpar.nlevsoi, 1]))
        end
    end
    dump = joinpath(DUMPDIR_CH4, "bgcdump_after_ch4_n$(nstep).nc")
    nlev = CLM.varpar.nlevsoi
    rows = Tuple{String,Float64,Float64}[]
    NCDataset(dump, "r") do ds
        for (name, get) in LEVFIELDS
            f = _dv(ds, name); f === nothing && continue
            jf = get(inst)
            e, s = relerr(vec(jf[1, 1:nlev]), vec(Float64.(f[1:nlev, 1])))
            push!(rows, (name, e, s))
        end
        for (name, get) in COLFIELDS
            f = _dv(ds, name); f === nothing && continue
            jf = get(inst)
            e, s = relerr([jf[1]], [Float64(f[1])])
            push!(rows, (name, e, s))
        end
    end
    return rows
end

function main()
    nprobe = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 5
    steps  = round.(Int, range(N0 + 1, NLAST; length = nprobe))
    println("="^78)
    println("METHANE (ch4) FORTRAN PARITY — Bow TWS_inversion (finundated>0), boundary after_ch4")
    println("  reference : ", DUMPDIR_CH4)
    println("  window    : nstep ", N0, "..", NLAST, "  (probing ", nprobe, " steps)")
    println("  finundated: ~0.045 (TWS_inversion) — sat surface fluxes NON-vacuous")
    println("="^78)
    worst = Dict{String,Float64}(); scale = Dict{String,Float64}()
    for n in steps
        rows = score_step(n)
        for (name, e, s) in rows
            worst[name] = max(get(worst, name, 0.0), e)
            scale[name] = max(get(scale, name, 0.0), s)
        end
        @printf("  nstep %d scored (%d fields)\n", n, length(rows))
    end
    println()
    @printf("%-24s %14s %14s  %s\n", "FIELD", "worst rel.err", "|F| scale", "verdict")
    println("-"^78)
    nok = 0; ntot = 0
    for (name, _) in vcat(LEVFIELDS, COLFIELDS)
        haskey(worst, name) || continue
        e, s = worst[name], scale[name]; ntot += 1
        vac = s == 0.0; ok = e <= 1e-9; ok && (nok += 1)
        v = vac ? "— (F≡0, vacuous)" : (ok ? "OK" : "DIVERGES")
        @printf("%-24s %14.3e %14.3e  %s\n", name, e, s, v)
    end
    println("-"^78)
    @printf("  %d / %d fields within 1e-9\n", nok, ntot)
end

main()
