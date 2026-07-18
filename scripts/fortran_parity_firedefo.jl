# ==========================================================================
# FIRE (Li2016) TROPICAL DEFORESTATION branch — Fortran parity for `dtrotr`.
#
# THE GAP (backlog A2, last unexercised fire branch)
#   `dtrotr_col` / `lfc` / `lfc2` / `fbac1` — the tropical land-use ("slash and
#   burn") fire term — is `F≡0` vacuous at every migrated SYMFLUENCE domain for
#   TWO independent reasons:
#     (1) `dtrotr_col` accumulates `-dwt_smoothed(p)` ONLY inside
#         `if (transient_landcover)`, i.e. only when `run_has_transient_landcover()`
#         is true, which needs a `flanduse_timeseries` + `do_transient_pfts`.
#         No migrated domain has one ⇒ `dwt ≡ 0` ⇒ `dtrotr ≡ 0`.
#     (2) the branch is gated on `trotr1+trotr2 > 0.6` — the column must be >60%
#         tropical broadleaf tree. Bow is boreal (NET temperate + C3 arctic grass).
#
# THE UNBLOCK — a SYNTHETIC SHARED input, not a new site
#   Same technique as the synthetic lnfm/hdm streams (§2) and the synthetic
#   `peatf=0.5` surfdata (§12) of docs/CH4_FIRE_PARITY.md: give BOTH codes the
#   same fabricated input and diff the kernel. `scripts/make_firedefo_inputs.jl`
#   writes, into the reference case dir:
#     * `surfdata_defo.nc`           — Bow's surfdata with PCT_NAT_PFT re-composed
#                                      to 40% BET-tropical + 30% BDT-tropical
#                                      (trotr1+trotr2 = 0.70 > 0.6), every other
#                                      natpft given a small nonzero weight.
#     * `landuse_timeseries_defo.nc` — 2 slices; YEAR 2203 drops the tropical
#                                      fraction by 2.6 pp ⇒ dwt_smoothed < 0
#                                      ⇒ dtrotr > 0.
#   BOTH codes read BOTH files: CTSM via fsurdat/flanduse_timeseries, CLM.jl via
#   `fsurdat=` and the `flanduse=`/`flanduse_year=` hook in
#   fortran_parity_common.jl. So the diff measures the PORT, not the input.
#
# THE REFERENCE
#   run dir : SYMFLUENCE_data/clm_bgc_spinup/bow_ref_firedefo
#   config  : bow_ref_firepeat + do_transient_pfts + the two synthetic files +
#             finidat = the 2202-10-29 warm-CN restart & use_init_interp
#             (init_interp_method='general'). The run MUST CROSS Jan 1: CTSM
#             applies land-cover change on the first step of the year and then
#             DRIBBLES it over the year, so dwt_smoothed is nonzero only after a
#             year boundary inside the run.
#   window  : nstep 1536..1780, dtime = 3600 s, start 2202-10-29-00000
#             ⇒ nstep 1536 == 2203-01-01 00:00 exactly.
#   boundary: 'after_fire'
#
#   nstep 1536 is the `kmo==1 .and. kda==1 .and. mcsec==0` step: Fortran gives
#   DTROTR=0, TROTR1/2 = the PRE-transition 0.40/0.30, FAREA_BURNED=0. The probe
#   window deliberately starts at 1537 — the harness is a single-step oracle and
#   does NOT reproduce the Jan-1 land-cover-change event itself.
#
#   The window spans BOTH regimes of the branch, which is why it runs to 1780:
#     * nstep 1537..~1675 — `lfc > 0`: the conversion area has not yet all burned,
#       so `fbac1` clamps to 0 (fb*cli*cli_scale/secspday ~ 1e-7 can never exceed
#       2*lfc/dt ~ 1.4e-5) and the `lfc2` drawdown branch runs. LFC decays ~2.2e-4
#       per step from 0.0258.
#     * nstep ~1676..1780 — `lfc` exhausted: `lfc2` -> 0 and `fbac1`/`fbac` turn
#       ON (~1.2e-7), which is the branch that actually feeds the deforestation
#       CARBON flux (CNFireFluxes uses f = (fbac-baf_crop)/(1-cropf) when
#       transient_landcover). Probing only the first regime would leave
#       FBAC/FBAC1 F≡0 vacuous — a green scorecard blind to half the branch.
#
# THE INJECTION
#   `dwt_smoothed_patch` is a dribbler output, not a restart variable, so it is
#   injected from the dump's DWT_SMOOTHED. `lfc_col` is a restart variable but is
#   NOT in the pdump before_step set, so it is injected from the PREVIOUS step's
#   after_fire dump (lfc is updated in place by CNFireFluxes, so step n's input
#   is step n-1's output). Everything else — trotr1/trotr2 (from itype+wtcol),
#   dtrotr (from dwt_smoothed), cli, farea_burned, fbac/fbac1, lfc2 — is COMPUTED
#   by CLM.jl and diffed.
#
# Usage: julia +1.12 --project=. scripts/fortran_parity_firedefo.jl [nprobe]
# ==========================================================================

include(joinpath(@__DIR__, "fortran_parity_common.jl"))
using Printf

const DUMPDIR_DEFO = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bow_ref_firedefo"
const FSURDAT_DEFO = joinpath(DUMPDIR_DEFO, "surfdata_defo.nc")
const FLANDUSE     = joinpath(DUMPDIR_DEFO, "landuse_timeseries_defo.nc")

const N0    = 1536                       # 2203-01-01 00:00 (the jan1-0 step)
const NFIRST = 1537                      # first step with dtrotr > 0
const NLAST = 1780
const LANDUSE_YEAR = 2203                # post-transition flanduse slice

# --- clock alignment (this branch is DATE-GATED, so this must be exact) ------
# The Fortran dumps record `ymd`/`tod` from get_curr_date, which in CTSM is the
# date at the END of the dumped step: dump n1536 has tod=0, n1537 tod=3600,
# n1538 tod=7200 (verified in the files). CNFireLi2016 keys `lfc` off exactly
# that value — lfc is (re)set ONLY on the step with `mcsec == dt`, i.e. nstep
# 1537, and merely CARRIED (decremented by CNFireFluxes) on every later step.
#
# run_one_parity_step! takes the PRE-advance clock and calls advance_timestep!
# before reading get_curr_date, so `step_date` must be the step's START:
#     step_date(n) = end_of_step(n) - dt
# (fortran_parity_firepeat.jl passes the END-of-step time here; nothing in the
# peat branch is date-gated, so that one hour never showed up there. It does
# here: an hour late and `mcsec == dt` never fires, leaving lfc ≡ 0.)
defo_step_date(n::Int) = DateTime(2203, 1, 1, 0) + Hour(n - NFIRST)

# --- forcing alignment ------------------------------------------------------
# datm: year_align=2002, year_first=2002, year_last=2009 (8-year recycle).
# model year 2203 -> (2203-2002) mod 8 = 1 -> forcing year 2003.
# fortran_parity_firepeat.jl reads the forcing at the step START (its
# forcing_date is its step_date minus one hour, and its step_date is the
# end-of-step time). Since defo_step_date already IS the step start, the
# forcing time is the same wall clock, just in the mapped forcing year.
const FFORCING2003 = joinpath(dirname(FFORCING), "clmforc.2003.nc")
defo_forcing_date(n::Int) = DateTime(2003, 1, 1, 0) + Hour(n - NFIRST)

const FNDEP = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/ndepdata/" *
              "fndep_clm_UNIFORM1e-8_0.9x1.25_yr2000_CLMjl-parity.nc"
const FIRED = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/firedata"
const FLNFM = joinpath(FIRED, "clmforc.UNIFORM_lnfm_0.9x1.25_CLMjl-parity.nc")
const FHDM  = joinpath(FIRED, "clmforc.UNIFORM_hdm_0.9x1.25_CLMjl-parity.nc")
_dv(ds, n) = haskey(ds, n) ? Array(ds[n]) : nothing

# lifire_inparm namelist values of the reference run (identical to firepeat).
function apply_bow_lifire_inparm!(inst)
    k = inst.cnfire_const
    k.boreal_peatfire_c            = 0.09e-4
    k.non_boreal_peatfire_c        = 0.17e-3
    k.bt_min                       = 0.85
    k.bt_max                       = 0.98
    k.cli_scale                    = 0.033
    k.cmb_cmplt_fact_cwd           = 0.28
    k.cmb_cmplt_fact_litter        = 0.5
    k.cropfire_a1                  = 1.6e-4
    k.occur_hi_gdp_tree            = 0.33
    k.lfuel                        = 105.0
    k.ufuel                        = 1050.0
    k.pot_hmn_ign_counts_alpha     = 0.010
    k.rh_low                       = 30.0
    k.rh_hgh                       = 80.0
    k.max_rh30_affecting_fuel      = 90.0
    k.defo_fire_precip_thresh_bdt  = 1.8
    k.defo_fire_precip_thresh_bet  = 4.0
    k.borpeat_fire_soilmoist_denom = 0.3
    k.nonborpeat_fire_precip_denom = 1.0
    return inst
end

"""Inject the two deforestation-branch inputs CLM.jl cannot regenerate in a
single-step oracle: `dwt_smoothed_patch` (a dribbler output) and `lfc_col`
(carried across steps by CNFireFluxes, absent from the before_step pdump)."""
function inject_defo_state!(inst, nstep::Int)
    vs = inst.bgc_vegetation.cnveg_state_inst
    NCDataset(joinpath(DUMPDIR_DEFO, "bgcdump_after_fire_n$(nstep).nc"), "r") do ds
        dwt = _dv(ds, "DWT_SMOOTHED")
        if dwt !== nothing
            for p in 1:min(length(vs.dwt_smoothed_patch), length(dwt))
                vs.dwt_smoothed_patch[p] = Float64(dwt[p])
            end
        end
    end
    # lfc INPUT to step n == lfc OUTPUT of step n-1.
    prev = joinpath(DUMPDIR_DEFO, "bgcdump_after_fire_n$(nstep - 1).nc")
    if isfile(prev)
        NCDataset(prev, "r") do ds
            l = _dv(ds, "LFC")
            l === nothing || (vs.lfc_col[1] = Float64(l[1]))
        end
    end
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

const COLFIELDS = [
    ("TROTR1",       i -> i.bgc_vegetation.cnveg_state_inst.trotr1_col),
    ("TROTR2",       i -> i.bgc_vegetation.cnveg_state_inst.trotr2_col),
    ("DTROTR",       i -> i.bgc_vegetation.cnveg_state_inst.dtrotr_col),
    ("LFC",          i -> i.bgc_vegetation.cnveg_state_inst.lfc_col),
    ("LFC2",         i -> i.bgc_vegetation.cnveg_state_inst.lfc2_col),
    ("FBAC",         i -> i.bgc_vegetation.cnveg_state_inst.fbac_col),
    ("FBAC1",        i -> i.bgc_vegetation.cnveg_state_inst.fbac1_col),
    ("FAREA_BURNED", i -> i.bgc_vegetation.cnveg_state_inst.farea_burned_col),
    ("NFIRE",        i -> i.bgc_vegetation.cnveg_state_inst.nfire_col),
    ("BAF_CROP",     i -> i.bgc_vegetation.cnveg_state_inst.baf_crop_col),
    ("BAF_PEATF",    i -> i.bgc_vegetation.cnveg_state_inst.baf_peatf_col),
]

function score_step(nstep::Int)
    inst, _ = run_one_parity_step!(nstep;
        use_cn = true, use_lch4 = true, use_hydrstress = true,
        dumpdir = DUMPDIR_DEFO,
        fsurdat = FSURDAT_DEFO,
        flanduse = FLANDUSE, flanduse_year = LANDUSE_YEAR,
        step_date = defo_step_date(nstep),
        forcing_file = FFORCING2003, forcing_date = defo_forcing_date(nstep),
        fndep = FNDEP,
        cnfire_method = :li2016, flnfm = FLNFM, fhdm = FHDM,
        pre_step_hook = (i, b, df) -> begin
            apply_bow_lifire_inparm!(i)
            CLM.restore_atm2lnd_runmean_accum!(i.atm2lnd, df, b; dtime = 3600)
            inject_defo_state!(i, nstep)
        end)

    dump = joinpath(DUMPDIR_DEFO, "bgcdump_after_fire_n$(nstep).nc")
    rows = Tuple{String,Float64,Float64}[]
    NCDataset(dump, "r") do ds
        for (name, get) in COLFIELDS
            f = _dv(ds, name); f === nothing && continue
            jf = get(inst); isempty(jf) && continue
            e, s = relerr([jf[1]], [Float64(f[1])])
            push!(rows, (name, e, s))
        end
    end
    return rows
end

function main()
    nprobe = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 5
    steps  = round.(Int, range(NFIRST, NLAST; length = nprobe))
    println("="^80)
    println("FIRE (Li2016) TROPICAL DEFORESTATION BRANCH — Fortran parity (dtrotr), after_fire")
    println("  reference : ", DUMPDIR_DEFO)
    println("  window    : nstep ", NFIRST, "..", NLAST, "  (probing ", nprobe, " steps)")
    println("  synthetic : ", basename(FSURDAT_DEFO), " + ", basename(FLANDUSE), " (read by BOTH codes)")
    println("="^80)
    worst = Dict{String,Float64}(); scale = Dict{String,Float64}()
    for n in steps
        for (name, e, s) in score_step(n)
            worst[name] = max(get(worst, name, 0.0), e)
            scale[name] = max(get(scale, name, 0.0), s)
        end
        @printf("  nstep %d scored\n", n)
    end
    println()
    @printf("%-16s %14s %14s  %s\n", "FIELD", "worst rel.err", "|F| scale", "verdict")
    println("-"^70)
    nok = 0; ntot = 0
    for (name, _) in COLFIELDS
        haskey(worst, name) || continue
        e, s = worst[name], scale[name]
        vac = s == 0.0
        vac && (println(@sprintf("%-16s %14s %14.3e  %s", name, "—", s, "F≡0 (vacuous)")); continue)
        ntot += 1
        ok = e <= 1e-9; ok && (nok += 1)
        @printf("%-16s %14.3e %14.3e  %s\n", name, e, s, ok ? "OK" : "DIVERGES")
    end
    println("-"^70)
    @printf("  %d / %d NON-VACUOUS fields within 1e-9\n", nok, ntot)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
