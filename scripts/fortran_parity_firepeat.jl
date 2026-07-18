# ==========================================================================
# FIRE (Li2016) PEATLAND branch — Fortran parity for `baf_peatf`.
#
# THE GAP (backlog A2)
#   The peatland-fire term `baf_peatf` is identically zero at every migrated
#   SYMFLUENCE single-point site, because the surfdata `peatf` field is 0 in
#   ALL of them (the SYMFLUENCE surfdata pipeline never populates it). So the
#   `_fire_peatland_kernel!` branch had never been diffed against Fortran — it
#   was `F≡0` vacuous in the standing fire scorecard (`fortran_parity_fire.jl`).
#
# THE REFERENCE
#   run dir : SYMFLUENCE_data/clm_bgc_spinup/bow_ref_firepeat
#   config  : IDENTICAL to bow_ref_ch4tws (fire_method='li2016crufrc', warm-CN
#             from the converged 2202-01-01 BGC restart, use_lch4=.true.) EXCEPT
#             the surfdata `peatf` field is set to 0.5. Both CLM.jl and CTSM read
#             the SAME edited surfdata, so the diff measures the peatland-fire
#             KERNEL, not the surfdata value — exactly like the synthetic uniform
#             lnfm/hdm fire streams already used in fortran_parity_fire.jl.
#             `peatf` feeds ONLY the CNFire modules in CTSM (verified), so the
#             warm-CN state is otherwise unperturbed.
#   window  : nstep 4720..4888  (2202-07-16 16:00 -> +7 days), dtime = 3600 s
#   boundary: 'after_fire' (bgcdumpMod, immediately after EcosystemDynamicsPreDrainage)
#
#   Fortran BAF_PEATF across the window = 1.67382e-10 /s (nonzero, boreal branch:
#   Bow 51.4N > borealat 40, summer soil unfrozen).
#
# THE INJECTION
#   Julia's surfdata read gives peatf=0 (the un-edited Bow surfdata the harness
#   loads), so the pre-step hook sets cnfire_li2014.peatf_lf_col .= 0.5 to match
#   the Fortran reference. baf_peatf is a pure per-column function of injected
#   before-step state (prec60, wf2, tsoi17, fsat) + peatf, so the cold-CN
#   confound does NOT enter it.
#
# Usage: julia +1.12 --project=. scripts/fortran_parity_firepeat.jl [nprobe]
# ==========================================================================

include(joinpath(@__DIR__, "fortran_parity_common.jl"))
using Printf

const DUMPDIR_PEAT = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bow_ref_firepeat"
const N0    = 4720
const NLAST = 4888
const DATE0 = DateTime(2003, 7, 16, 16)
const PEATF_INJECT = 0.5

# --- forcing alignment (identical to fortran_parity_fire.jl §5) -------------
const FFORCING2002 = joinpath(dirname(FFORCING), "clmforc.2002.nc")
const FORCE_DATE0  = DateTime(2002, 7, 16, 15)
peat_forcing_date(n::Int) = FORCE_DATE0 + Hour(n - N0)

const FNDEP = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/ndepdata/" *
              "fndep_clm_UNIFORM1e-8_0.9x1.25_yr2000_CLMjl-parity.nc"
const FIRED = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/firedata"
const FLNFM = joinpath(FIRED, "clmforc.UNIFORM_lnfm_0.9x1.25_CLMjl-parity.nc")
const FHDM  = joinpath(FIRED, "clmforc.UNIFORM_hdm_0.9x1.25_CLMjl-parity.nc")
_dv(ds, n) = haskey(ds, n) ? Array(ds[n]) : nothing

# lifire_inparm namelist values of the reference run (matches Fortran &lifire_inparm)
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

inject_fire_accum!(inst, bounds, dumpfile::String) =
    (CLM.restore_atm2lnd_runmean_accum!(inst.atm2lnd, dumpfile, bounds; dtime = 3600); inst)

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
    ("BAF_PEATF",    i -> i.bgc_vegetation.cnveg_state_inst.baf_peatf_col),
    ("FAREA_BURNED", i -> i.bgc_vegetation.cnveg_state_inst.farea_burned_col),
    ("NFIRE",        i -> i.bgc_vegetation.cnveg_state_inst.nfire_col),
    ("FBAC",         i -> i.bgc_vegetation.cnveg_state_inst.fbac_col),
    ("SOMC_FIRE",    i -> i.soilbiogeochem_carbonflux.somc_fire_col),
]

function score_step(nstep::Int)
    inst, _ = run_one_parity_step!(nstep;
        use_cn = true, use_lch4 = true, use_hydrstress = true,
        dumpdir = DUMPDIR_PEAT,
        step_date = DATE0 + Hour(nstep - N0),
        forcing_file = FFORCING2002, forcing_date = peat_forcing_date(nstep),
        fndep = FNDEP,
        cnfire_method = :li2016, flnfm = FLNFM, fhdm = FHDM,
        pre_step_hook = (i, b, df) -> begin
            apply_bow_lifire_inparm!(i)
            inject_fire_accum!(i, b, df)
            fill!(i.cnfire_li2014.peatf_lf_col, PEATF_INJECT)   # match Fortran surfdata peatf
        end)

    dump = joinpath(DUMPDIR_PEAT, "bgcdump_after_fire_n$(nstep).nc")
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
    steps  = round.(Int, range(N0 + 1, NLAST; length = nprobe))
    println("="^80)
    println("FIRE (Li2016) PEATLAND-BRANCH FORTRAN PARITY — Bow, peatf=$PEATF_INJECT, after_fire")
    println("  reference : ", DUMPDIR_PEAT)
    println("  window    : nstep ", N0, "..", NLAST, "  (probing ", nprobe, " steps)")
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
        e, s = worst[name], scale[name]; ntot += 1
        vac = s == 0.0; ok = e <= 1e-9; ok && (nok += 1)
        v = vac ? "— (F≡0, vacuous)" : (ok ? "OK" : "DIVERGES")
        @printf("%-16s %14.3e %14.3e  %s\n", name, e, s, v)
    end
    println("-"^70)
    @printf("  %d / %d fields within 1e-9\n", nok, ntot)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
