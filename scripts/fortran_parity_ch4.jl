# ==========================================================================
# METHANE (ch4) Fortran parity — the FIRST EVER diff of the CH4 subsystem.
#
# THE GAP
#   CH4 had never been compared to Fortran. Not because it was hard, but because
#   the reference run had `use_lch4 = .false.`, so the CLM5 restart carried ZERO
#   ch4 variables and there was nothing to diff against. The "GPU-validated"
#   claim in the README was internal-consistency only.
#
# THE REFERENCE (see docs/CH4_FIRE_PARITY.md)
#   run dir : SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_firech4
#   config  : use_lch4=.true., fire_method='li2016crufrc', finundation_method='h2osfc',
#             use_cn/use_fun/use_nitrif_denitrif/use_hydrstress/use_luna=.true.,
#             use_crop=.false., soil_decomp_method='CENTURYKoven2013'
#   start   : startup from the converged 2202-01-01 BGC-spinup restart (branch is
#             impossible: a branch run ABORTS on the ch4 restart fields that the
#             use_lch4=.false. spinup never wrote). 196 days of spin-up with ch4+fire
#             ON precede the window, so the ch4 state and the fire accumulators
#             (prec10/prec60/rh30) are physically converged, not cold.
#   window  : nstep 4720..4888  (2202-07-16 16:00 -> +7 days), dtime = 3600 s
#   boundary: 'after_ch4'  -- written by bgcdumpMod immediately after the ch4() call.
#
# NON-VACUITY (this test CAN fail)
#   Bow-at-Banff is a dry mountain site: frac_h2osfc ~ 0, so finundated ~ 0 and the
#   column is a net CH4 SINK, not a source. That is the physically right answer, and
#   it is NOT vacuous: BOTH the saturated and the unsaturated columns are integrated
#   every step regardless of finundated, so production (sat), oxidation (unsat),
#   aerenchyma, ebullition and the two tridiagonal transport solves are all exercised
#   with non-zero values. What is NOT exercised here is the high-inundation WETLAND
#   regime -- that needs a peatland/wetland site (see the doc).
#
# Usage: julia --project=. scripts/fortran_parity_ch4.jl [nprobe]
# ==========================================================================

include(joinpath(@__DIR__, "fortran_parity_common.jl"))
using Printf

const DUMPDIR_CH4 = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_firech4"
const N0   = 4720                                # first nstep of the window
const NLAST = 4888
const DATE0 = DateTime(2003, 7, 16, 16)          # forcing date aligned to model 2202-07-16 16:00
const FNDEP = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/ndepdata/" *
              "fndep_clm_UNIFORM1e-8_0.9x1.25_yr2000_CLMjl-parity.nc"

# NB: keep the shape. Both the restart-style pdump and the bgcdump store
# level-resolved column fields as (levels, column) when read into Julia
# (CDL declares them `VAR(column, lev*)`, which NCDatasets reverses).
_dv(ds, n) = haskey(ds, n) ? Array(ds[n]) : nothing

# --------------------------------------------------------------------------
# Inject the Fortran CH4 prognostic state.
#
# NOTE / KNOWN GAP: `src/infrastructure/fortran_restart.jl` (read_fortran_restart!)
# has NO ch4 handling at all -- CLM.jl cannot currently restart a methane run from a
# Fortran restart file. That is a real porting gap, reported in the scorecard. For the
# oracle we inject the ch4 state here in the harness so the single-step diff starts
# from EXACTLY Fortran's state.
# --------------------------------------------------------------------------
function inject_ch4_dump!(inst, bounds, dumpfile::String)
    ch4 = inst.ch4
    nc, np = bounds.endc, bounds.endp
    nlev = CLM.varpar.nlevsoi

    NCDataset(dumpfile, "r") do ds
        # dump is (levels, column); Julia state is [column, level].
        # The dump only carries the ACTIVE subgrid subset, so bound every loop by
        # the dump's own extent -- Julia's bounds.endc/endp include inactive units.
        col2d!(dst, name) = begin
            a = _dv(ds, name); a === nothing && return
            for c in 1:min(nc, size(a, 2)), j in 1:min(nlev, size(a, 1))
                v = a[j, c]
                ismissing(v) || (dst[c, j] = Float64(v))
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
        col2d!(ch4.o2_decomp_depth_sat_col,   "O2_DECOMP_DEPTH_SAT")
        col2d!(ch4.o2_decomp_depth_unsat_col, "O2_DECOMP_DEPTH_UNSAT")
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
    end
    # ch4_first_time_grc must be false: the Fortran state is warm, so the dfsat
    # first-step special case must NOT fire.
    fill!(ch4.ch4_first_time_grc, false)
    return inst
end

# --------------------------------------------------------------------------
# Profile-scaled worst relative error (same rationale as fortran_parity_ncycle.jl:
# per-element scaling is meaningless where the Fortran value is ~1e-30).
# --------------------------------------------------------------------------
function relerr(j, f)
    jj = vec(Float64.(j)); ff = vec(Float64.(f))
    n = min(length(jj), length(ff))
    scale = 0.0
    for i in 1:n
        isfinite(ff[i]) && (scale = max(scale, abs(ff[i])))
    end
    scale <= 0 && return (0.0, 0.0)          # Fortran identically zero
    worst = 0.0
    for i in 1:n
        a, b = jj[i], ff[i]
        (isfinite(a) && isfinite(b)) || continue
        worst = max(worst, abs(a - b) / scale)
    end
    return (worst, scale)
end

# field table: dump name => (julia getter, kind)
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
]

function score_step(nstep::Int)
    (inst, bounds) = run_one_parity_step_ch4!(nstep)

    dump = joinpath(DUMPDIR_CH4, "bgcdump_after_ch4_n$(nstep).nc")
    nlev = CLM.varpar.nlevsoi
    rows = Tuple{String,Float64,Float64}[]     # (name, relerr, fortran scale)

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

# a ch4-flavoured single-step runner: inject the biogeophysics state (the 16-field
# oracle registry) AND the ch4 prognostics, then run exactly one Julia step.
function run_one_parity_step_ch4!(nstep::Int)
    inst, bounds = run_one_parity_step!(nstep;
        use_cn = true, use_lch4 = true, use_hydrstress = true,
        dumpdir = DUMPDIR_CH4,
        step_date = DATE0 + Hour(nstep - N0),
        fndep = FNDEP,
        pre_step_hook = (i, b, df) -> inject_ch4_dump!(i, b, df))
    return (inst, bounds)
end

function main()
    nprobe = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 5
    steps  = round.(Int, range(N0 + 1, NLAST; length = nprobe))

    println("="^78)
    println("METHANE (ch4) FORTRAN PARITY — Bow-at-Banff, summer 2202, boundary after_ch4")
    println("  reference : ", DUMPDIR_CH4)
    println("  window    : nstep ", N0, "..", NLAST, "  (probing ", nprobe, " steps)")
    println("="^78)

    worst = Dict{String,Float64}()
    scale = Dict{String,Float64}()
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
        e, s = worst[name], scale[name]
        ntot += 1
        vac = s == 0.0
        ok  = e <= 1e-9
        ok && (nok += 1)
        v = vac ? "— (F≡0, vacuous)" : (ok ? "OK" : "DIVERGES")
        @printf("%-24s %14.3e %14.3e  %s\n", name, e, s, v)
    end
    println("-"^78)
    @printf("  %d / %d fields within 1e-9\n", nok, ntot)
end

main()
