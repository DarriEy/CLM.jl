# ==========================================================================
# FIRE (Li2016) Fortran parity — the FIRST EVER diff of the fire subsystem.
#
# THE GAP
#   Fire had never been compared to Fortran, for two independent reasons:
#     (a) the reference run had `fire_method = 'nofire'`, so Fortran never
#         ran fire either -- there was nothing to diff against; and
#     (b) on the Julia side the whole Li-family chain was a DEAD PORT: the
#         `_fire_active` gate in cn_driver.jl could never become true, because
#         `cnfire_method` had no route out of :nofire and the fire data bundle
#         was never constructed anywhere in src/.
#   Both are fixed now. This harness is the scorecard.
#
# THE REFERENCE (see docs/CH4_FIRE_PARITY.md)
#   run dir : SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_firech4
#   config  : fire_method='li2016crufrc', use_lch4=.true., use_cn/use_fun/
#             use_nitrif_denitrif/use_hydrstress/use_luna=.true., use_crop=.false.
#   window  : nstep 4720..4888  (2202-07-16 16:00 -> +7 days), dtime = 3600 s
#   boundary: 'after_fire' -- written by bgcdumpMod immediately after
#             EcosystemDynamicsPreDrainage. NB: CNFireArea + CNFireFluxes run
#             INSIDE CNDriverNoLeaching, i.e. inside PreDrainage -- NOT in
#             PostDrainage. That is the first point at which the fire fluxes exist.
#
# LIGHTNING / POPULATION DENSITY
#   The real lnfm/hdm inputdata files are absent from the local tree. BOTH codes are
#   driven by the SAME synthetic uniform-in-space, constant-in-time streams
#   (scripts/validation/make_fire_streams.jl):
#       lnfm = 4.57e-4 counts/km^2/hr,  hdm = 4.2 counts/km^2
#   Because the field is uniform+constant, ESMF bilinear regrid and linear time
#   interpolation are EXACT, so Fortran's forc_lnfm/forc_hdm are known analytically.
#   The diff therefore measures FIRE PHYSICS, not regrid fidelity.
#
# NON-VACUITY: Bow DOES burn under this forcing.  Fortran gives
#   FAREA_BURNED ~ 1.7e-8 s^-1, NFIRE ~ 9.9e-9, and non-zero per-patch fire C fluxes
#   (M_LEAFC_TO_FIRE, M_DEADSTEMC_TO_FIRE, ...). A 0-vs-0 pass is not possible here.
#
# Usage: julia --project=. scripts/fortran_parity_fire.jl [nprobe]
# ==========================================================================

include(joinpath(@__DIR__, "fortran_parity_common.jl"))
using Printf

const DUMPDIR_FIRE = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_firech4"
const N0    = 4720
const NLAST = 4888
const DATE0 = DateTime(2003, 7, 16, 16)
const FNDEP = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/ndepdata/" *
              "fndep_clm_UNIFORM1e-8_0.9x1.25_yr2000_CLMjl-parity.nc"
const FIRED = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/firedata"
const FLNFM = joinpath(FIRED, "clmforc.UNIFORM_lnfm_0.9x1.25_CLMjl-parity.nc")
const FHDM  = joinpath(FIRED, "clmforc.UNIFORM_hdm_0.9x1.25_CLMjl-parity.nc")

_dv(ds, n) = haskey(ds, n) ? Array(ds[n]) : nothing

"""
Apply the Bow reference run's `&lifire_inparm` namelist to the Julia fire constants.

This is CONFIGURATION, not tuning. The Bow case is built with
`-lnd_tuning_mode clm5_0_GSWP3v1`, whose `&lifire_inparm` overrides 10 of the
CTSM *module* defaults (which is what `CNFireConstData` in fire_base.jl carries).
Both codes must be driven by the SAME namelist or the comparison is meaningless —
exactly like `br_root`, `leafresp_method` and the rooting profile already handled
in `build_bow_inst`.

Verbatim from `bgc_ref_firech4/lnd_in`:

    boreal_peatfire_c        = 0.09d-4     non_boreal_peatfire_c    = 0.17d-3
    bt_min = 0.85  bt_max = 0.98           cli_scale                = 0.033
    cmb_cmplt_fact_cwd       = 0.28        cmb_cmplt_fact_litter    = 0.5
    cropfire_a1              = 1.6d-4      occur_hi_gdp_tree        = 0.33
    lfuel = 105.  ufuel = 1050.            pot_hmn_ign_counts_alpha = 0.010
    rh_low = 30.  rh_hgh = 80.             max_rh30_affecting_fuel  = 90.
    defo_fire_precip_thresh_bdt = 1.8      defo_fire_precip_thresh_bet = 4.0
    borpeat_fire_soilmoist_denom = 0.3     nonborpeat_fire_precip_denom = 1.0
"""
# ---------------------------------------------------------------------------
# Restore the Fortran runmean forcing ACCUMULATORS (rh30 / prec10 / prec30 /
# prec60) from the dump, via the ported `accumulRest` restore path
# (CLM.restore_atm2lnd_runmean_accum!, src/infrastructure/accumul_restart.jl).
#
# `read_fortran_restart!` maps ~170 prognostic restart variables but NONE of the
# `accumul` fields, so before this the Julia side started them cold. The restore
# copies each persisted `<NAME>_VALUE` into its Atm2LndData patch vector and reads
# `<NAME>_NSTEPS`. The single-step oracle reads the fire accumulators (during
# CNFireArea, inside PreDrainage) BEFORE the end-of-step accumulator update, so
# the VALUE restore is what the fire physics sees; `warm_nstep` is the counter a
# multi-step restart would run at so `accum_runmean` continues the mean.
#
# NB (see docs/CH4_FIRE_PARITY.md re-score): at Bow the fuel load keeps `afuel=0`,
# so `rh30` carries ZERO weight in `fire_m`. Restoring the accumulators is a real
# porting fix but does not move NFIRE here — the residual is the fire root-wetness
# `btran2`, not the rh30 accumulator.
# ---------------------------------------------------------------------------
function inject_fire_accum!(inst, bounds, dumpfile::String)
    CLM.restore_atm2lnd_runmean_accum!(inst.atm2lnd, dumpfile, bounds; dtime = 3600)
    return inst
end

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

function relerr(j, f)
    jj = vec(Float64.(j)); ff = vec(Float64.(f))
    n = min(length(jj), length(ff))
    scale = 0.0
    for i in 1:n
        isfinite(ff[i]) && (scale = max(scale, abs(ff[i])))
    end
    scale <= 0 && return (0.0, 0.0)
    worst = 0.0
    for i in 1:n
        a, b = jj[i], ff[i]
        (isfinite(a) && isfinite(b)) || continue
        worst = max(worst, abs(a - b) / scale)
    end
    return (worst, scale)
end

# ---- column-level fire diagnostics (CNVegStateData unless noted) ----
const COLFIELDS = [
    ("FAREA_BURNED", i -> i.bgc_vegetation.cnveg_state_inst.farea_burned_col),
    ("NFIRE",        i -> i.bgc_vegetation.cnveg_state_inst.nfire_col),
    ("BAF_CROP",     i -> i.bgc_vegetation.cnveg_state_inst.baf_crop_col),
    ("BAF_PEATF",    i -> i.bgc_vegetation.cnveg_state_inst.baf_peatf_col),
    ("LFC",          i -> i.bgc_vegetation.cnveg_state_inst.lfc_col),
    ("FBAC",         i -> i.bgc_vegetation.cnveg_state_inst.fbac_col),
    ("FBAC1",        i -> i.bgc_vegetation.cnveg_state_inst.fbac1_col),
    ("CROPF",        i -> i.bgc_vegetation.cnveg_state_inst.cropf_col),
    ("LGDP",         i -> i.bgc_vegetation.cnveg_state_inst.lgdp_col),
    ("LGDP1",        i -> i.bgc_vegetation.cnveg_state_inst.lgdp1_col),
    ("LPOP",         i -> i.bgc_vegetation.cnveg_state_inst.lpop_col),
    ("FSR_COL",      i -> i.bgc_vegetation.cnveg_state_inst.fsr_col),
    ("FD_COL",       i -> i.bgc_vegetation.cnveg_state_inst.fd_col),
    ("LFWT",         i -> i.bgc_vegetation.cnveg_state_inst.lfwt_col),
    ("TROTR1",       i -> i.bgc_vegetation.cnveg_state_inst.trotr1_col),
    ("TROTR2",       i -> i.bgc_vegetation.cnveg_state_inst.trotr2_col),
    ("DTROTR",       i -> i.bgc_vegetation.cnveg_state_inst.dtrotr_col),
    ("WTLF",         i -> i.bgc_vegetation.cnveg_state_inst.wtlf_col),
    ("FUELC",        i -> i.bgc_vegetation.cnveg_carbonstate_inst.fuelc_col),
    ("FUELC_CROP",   i -> i.bgc_vegetation.cnveg_carbonstate_inst.fuelc_crop_col),
    ("ROOTC_COL",    i -> i.bgc_vegetation.cnveg_carbonstate_inst.rootc_col),
    ("SOMC_FIRE",    i -> i.soilbiogeochem_carbonflux.somc_fire_col),
]

# ---- patch-level fire C/N fluxes ----
const PFTFIELDS = [
    ("M_LEAFC_TO_FIRE",        i -> i.bgc_vegetation.cnveg_carbonflux_inst.m_leafc_to_fire_patch),
    ("M_FROOTC_TO_FIRE",       i -> i.bgc_vegetation.cnveg_carbonflux_inst.m_frootc_to_fire_patch),
    ("M_LIVESTEMC_TO_FIRE",    i -> i.bgc_vegetation.cnveg_carbonflux_inst.m_livestemc_to_fire_patch),
    ("M_DEADSTEMC_TO_FIRE",    i -> i.bgc_vegetation.cnveg_carbonflux_inst.m_deadstemc_to_fire_patch),
    ("M_LEAFC_TO_LITTER_FIRE", i -> i.bgc_vegetation.cnveg_carbonflux_inst.m_leafc_to_litter_fire_patch),
]

# ---- level-resolved fire fluxes (nlevdecomp) ----
const LEVFIELDS = [
    ("FIRE_MORTALITY_C_TO_CWDC", i -> i.bgc_vegetation.cnveg_carbonflux_inst.fire_mortality_c_to_cwdc_col),
    ("FIRE_MORTALITY_N_TO_CWDN", i -> i.bgc_vegetation.cnveg_nitrogenflux_inst.fire_mortality_n_to_cwdn_col),
]

function score_step(nstep::Int)
    inst, _ = run_one_parity_step!(nstep;
        use_cn = true, use_lch4 = true, use_hydrstress = true,
        dumpdir = DUMPDIR_FIRE,
        step_date = DATE0 + Hour(nstep - N0),
        fndep = FNDEP,
        cnfire_method = :li2016, flnfm = FLNFM, fhdm = FHDM,
        pre_step_hook = (i, b, df) -> begin
            apply_bow_lifire_inparm!(i)
            inject_fire_accum!(i, b, df)
        end)

    dump = joinpath(DUMPDIR_FIRE, "bgcdump_after_fire_n$(nstep).nc")
    nlevd = CLM.varpar.nlevdecomp
    rows = Tuple{String,Float64,Float64}[]

    NCDataset(dump, "r") do ds
        for (name, get) in COLFIELDS
            f = _dv(ds, name); f === nothing && continue
            jf = get(inst); isempty(jf) && continue
            e, s = relerr([jf[1]], [Float64(f[1])])
            push!(rows, (name, e, s))
        end
        for (name, get) in PFTFIELDS
            f = _dv(ds, name); f === nothing && continue
            jf = get(inst); isempty(jf) && continue
            np = min(length(jf), length(f))
            e, s = relerr(vec(jf[1:np]), vec(Float64.(f[1:np])))
            push!(rows, (name, e, s))
        end
        for (name, get) in LEVFIELDS
            f = _dv(ds, name); f === nothing && continue
            jf = get(inst); isempty(jf) && continue
            nl = min(nlevd, size(f, 1))
            e, s = relerr(vec(jf[1, 1:nl]), vec(Float64.(f[1:nl, 1])))
            push!(rows, (name, e, s))
        end
    end
    return rows
end

function main()
    nprobe = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 5
    steps  = round.(Int, range(N0 + 1, NLAST; length = nprobe))

    println("="^80)
    println("FIRE (Li2016) FORTRAN PARITY — Bow-at-Banff, summer 2202, boundary after_fire")
    println("  reference : ", DUMPDIR_FIRE)
    println("  window    : nstep ", N0, "..", NLAST, "  (probing ", nprobe, " steps)")
    println("  forcing   : lnfm = 4.57e-4 counts/km2/hr, hdm = 4.2 counts/km2 (synthetic, uniform)")
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
    @printf("%-26s %14s %14s  %s\n", "FIELD", "worst rel.err", "|F| scale", "verdict")
    println("-"^80)
    nok = 0; ntot = 0
    for (name, _) in vcat(COLFIELDS, PFTFIELDS, LEVFIELDS)
        haskey(worst, name) || continue
        e, s = worst[name], scale[name]
        ntot += 1
        vac = s == 0.0
        ok  = e <= 1e-9
        ok && (nok += 1)
        v = vac ? "— (F≡0, vacuous)" : (ok ? "OK" : "DIVERGES")
        @printf("%-26s %14.3e %14.3e  %s\n", name, e, s, v)
    end
    println("-"^80)
    @printf("  %d / %d fields within 1e-9\n", nok, ntot)
end

main()
