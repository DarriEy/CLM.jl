# ==========================================================================
# METHANE (ch4) Fortran parity at MerBleue peatland — the SECOND CH4 reference
# site, and the first at an organic-soil peatland.
#
# WHY THIS SITE
#   Bow-at-Banff is a dry mountain column with finundated == 0, so the whole
#   ebullition / sat-unsat-exchange half of the CH4 model is vacuous there
#   (docs/CH4_FIRE_PARITY.md §6/§7). MerBleue is an ombrotrophic bog: the plan
#   was to generate the first `finundated>0` reference so CH4-as-a-SOURCE
#   (transport / aerenchyma / ebullition) could finally be diffed.
#
# GATE RESULT (Phase 1, this PR)
#   The MerBleue CLM configuration (surfdata PCT_NATVEG=100, FMAX=0.5,
#   baseflow_scalar=0.001, 20SL_8.5m + use_bedrock) equilibrates with the water
#   table ~1.9 m DEEP even in mid-July; the h2osfc surface-water store stays at
#   0, so `finundated == 0` across the whole window — same degenerate regime as
#   Bow. `CH4_EBUL_*` and `CH4_DFSAT_FLUX` are therefore identically zero
#   (vacuous). CH4-as-a-source is NOT unlocked here (see the PR / doc).
#
#   What this run DOES give: a second, independent, organic-soil column on which
#   the dominant remaining lead — the SATURATED-column collapse (inject Fortran's
#   CONC_CH4_SAT, one Julia step drives it to ~0) — can be checked for
#   site-independence. The sat column is integrated every step regardless of
#   finundated, and here CONC_CH4_SAT / CONC_O2_SAT / production / oxidation /
#   aerenchyma / transport are all non-zero in Fortran.
#
# REFERENCE
#   run dir : SYMFLUENCE_data/clm_bgc_spinup/merbleue_ref_ch4
#   config  : use_lch4=.true., fire_method='nofire', finundation_method='h2osfc',
#             use_cn/use_fun/use_nitrif_denitrif/use_hydrstress/use_luna=.true.,
#             use_crop=.false., soil_decomp_method='CENTURYKoven2013'
#   start   : startup + use_init_interp from the converged SP peatland restart
#             (2017-01-01). No BGC restart for MerBleue exists, so CN + CH4 are
#             init_interp cold-started, then spun 4720 steps (to mid-July 2017).
#             NOTE: cold-CN => lower soil C => the sat column is WEAKER than Bow's
#             converged-CN reference; this is a reproduction check, not a
#             stronger test.
#   window  : nstep 4720..4888  (2017-07-16 16:00 -> +7 days), dtime = 3600 s
#   boundary: 'after_ch4'
#
# Usage: julia --project=. scripts/fortran_parity_ch4_merbleue.jl [nprobe]
# ==========================================================================

include(joinpath(@__DIR__, "fortran_parity_common.jl"))
using Printf

# The single-step injection oracle injects a PARTIAL Fortran state (CH4 prognostics
# + the 16-field biogeophysics registry) into a cold-CN-started instance. That is
# NOT N/C mass-conservative by construction, so the CN veg balance check (fatal
# since #223) would abort on the injected transient. Disable it here — this harness
# scores CH4 transport, not CN conservation.
CLM.cn_balance_check_enabled!(false)

const MERB       = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Peatland_MerBleue_Canada"
const MERB_PAR   = joinpath(MERB, "settings/CLM/parameters")
const DUMPDIR_CH4 = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/merbleue_ref_ch4"
const M_FSURDAT  = joinpath(MERB_PAR, "surfdata_clm.nc")
const M_FPARAM   = joinpath(MERB_PAR, "clm5_params.nc")
const M_FORCING  = joinpath(MERB, "data/forcing/CLM_input/clmforc.2017.nc")
const N0    = 4720
const NLAST = 4888
# model year 2017 maps DIRECTLY to forcing year 2017 (year_align=2016, start
# 2017-01-01), no recycle. n4720 end-of-step = 2017-07-16 16:00.
const DATE0 = DateTime(2017, 7, 16, 16)
const FNDEP = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/ndepdata/" *
              "fndep_clm_ZERO_0.9x1.25_yr2000_CLMjl-parity.nc"
# MerBleue lnd_in calibration (differs from Bow)
const M_BASEFLOW = 0.001
const M_INT_SNOW = 2000.0

_dv(ds, n) = haskey(ds, n) ? Array(ds[n]) : nothing

function inject_ch4_dump!(inst, bounds, dumpfile::String)
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
        col1d!(ch4.finundated_pre_snow_col, "FINUNDATED_PRE_SNOW")
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
    ("CH4_EBUL_TOTAL_SAT",   i -> i.ch4.ch4_ebul_total_sat_col),
    ("CH4_SURF_FLUX_TOT",    i -> i.ch4.ch4_surf_flux_tot_col),
    ("CH4_DFSAT_FLUX",       i -> i.ch4.ch4_dfsat_flux_col),
    ("FINUNDATED",           i -> i.ch4.finundated_col),
    ("TOTCOLCH4",            i -> i.ch4.totcolch4_col),
    ("ZWT_CH4_UNSAT",        i -> i.ch4.zwt_ch4_unsat_col),
]

function run_one_parity_step_ch4!(nstep::Int)
    run_one_parity_step!(nstep;
        use_cn = true, use_lch4 = true, use_hydrstress = true,
        dumpdir = DUMPDIR_CH4,
        step_date = DATE0 + Hour(nstep - N0),
        forcing_file = M_FORCING,
        fsurdat = M_FSURDAT, paramfile = M_FPARAM,
        baseflow = M_BASEFLOW, int_snow = M_INT_SNOW,
        fndep = FNDEP,
        pre_step_hook = (i, b, df) -> inject_ch4_dump!(i, b, df))
end

function score_step(nstep::Int)
    (inst, bounds) = run_one_parity_step_ch4!(nstep)
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
    println("METHANE (ch4) FORTRAN PARITY — MerBleue peatland, summer 2017, boundary after_ch4")
    println("  reference : ", DUMPDIR_CH4)
    println("  window    : nstep ", N0, "..", NLAST, "  (probing ", nprobe, " steps)")
    println("  NOTE      : finundated==0 at this site — ebullition/dfsat are VACUOUS.")
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
