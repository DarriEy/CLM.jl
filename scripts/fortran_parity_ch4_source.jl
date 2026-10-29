# ==========================================================================
# METHANE AS A SOURCE: multi-step Fortran parity with the CN confound removed.
#
# WHY THIS EXISTS (docs/CH4_FIRE_PARITY.md §16, HANDOFF §6 item 3)
#   fortran_parity_ch4_tws.jl is a SINGLE-STEP oracle: it re-injects the Fortran CH4
#   concentrations before every step. One hourly step barely moves a concentration,
#   so that oracle cannot tell whether the SOURCE term (production -> oxidation ->
#   aerenchyma / ebullition / diffusion -> surface flux) ACCUMULATES correctly. Its
#   one non-transport residual, CH4_PROD_SAT ~7.9%, was also never attributed. It was
#   blamed on "warm-start CN state", because production reads this step's somhr /
#   lithr / hr_vr and the Fortran dumps carry no HR.
#
# HOW THE CONFOUND IS REMOVED
#   1. CN STATE LOCKED TO FORTRAN EVERY STEP. Before each step the full Fortran restart
#      (pdump_before_step_n<N>: decomp C/N pools, veg C/N, soil T/water, ...) is injected
#      through CLM.read_fortran_restart!, the restart-injection route the CN/crop
#      harnesses use. Julia and Fortran decompose the IDENTICAL CN state each step.
#   2. CH4 STATE FREE-RUNNING. The CH4 prognostics (conc_ch4/o2 sat+unsat, o2stress,
#      layer_sat_lag, finundated_lag, annual averages) are injected ONCE, at the first
#      step, and then carried by Julia. Over the window, Julia's CH4 column integrates
#      its own production / oxidation / transport, so a wrong source term compounds and
#      shows up step by step.
#   3. HR IS ATTRIBUTED, NOT ASSUMED. In the saturated column every layer is below the
#      water table (jwt=0), so ch4Mod.F90 ch4_prod reduces to
#          CH4_PROD_SAT(j) = f_ch4_adj(j) * hr_vr(j) / catomw            (anoxia=.true.)
#          f_ch4_adj(j)    = min(f_ch4 * q10ch4^((T(j)-q10ch4base)/10)
#                                * [finundated_lag/finundated if lag<fin], 0.5)
#      With Fortran's own T_SOISNO (pdump_after_hydrologydrainage) and FINUNDATED[_LAG]
#      (bgcdump), this recovers Fortran's hr_vr per layer, which is then compared to Julia's
#      hr_vr. HR parity means the CN side is not the confound. Only the q10ch4 used
#      for the recovery is assumed: the paramfile's 1.33.
#
#   finundated is a forced input (TWS_inversion; tws back-derived from the dumped
#   FINUNDATED each step, as in fortran_parity_ch4_tws.jl).
#
# MODES
#   (default) locked : CN/biogeophysics re-injected every step, CH4 free-running.
#   free             : inject everything ONCE and free-run all of Julia (drift view).
#
# Usage: julia --project=. scripts/fortran_parity_ch4_source.jl [nsteps=168] [locked|free]
# ==========================================================================

include(joinpath(@__DIR__, "fortran_parity_common.jl"))
using Printf

# bow_ref_ch4o2 is the O2-demand ("bgcdump O2") rebuild of bow_ref_ch4tws: same case, same
# window, bit-identical in every field the two share, plus the O2-demand components
# (docs/CH4_FIRE_PARITY.md §17). Preferred when present; CH4_SOURCE_REFDIR overrides.
const _SD_SRC = get(ENV, "SYMFLUENCE_DATA", "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data")
const DUMPDIR_SRC = get(ENV, "CH4_SOURCE_REFDIR",
    isdir(joinpath(_SD_SRC, "clm_bgc_spinup", "bow_ref_ch4o2")) ?
        joinpath(_SD_SRC, "clm_bgc_spinup", "bow_ref_ch4o2") :
        joinpath(_SD_SRC, "clm_bgc_spinup", "bow_ref_ch4tws"))
const INV_FILE_SRC = "$(get(ENV, "CESM_INPUTDATA", "/Users/darri.eythorsson/projects/cesm-inputdata"))/lnd/clm2/paramdata/" *
                     "finundated_inversiondata_0.9x1.25_c170706.nc"
const FNDEP_SRC = "$(get(ENV, "CESM_INPUTDATA", "/Users/darri.eythorsson/projects/cesm-inputdata"))/lnd/clm2/ndepdata/" *
                  "fndep_clm_UNIFORM1e-8_0.9x1.25_yr2000_CLMjl-parity.nc"
const N0S    = 4720
const NLASTS = 4888
const DATE0S = DateTime(2003, 7, 16, 16)               # model clock (orbital) at N0
const FFORC_SRC   = joinpath(dirname(FFORCING), "clmforc.2002.nc")
const FORCE0_SRC  = DateTime(2002, 7, 16, 15)          # datm: forcing-year 2002, 1 h behind
const Q10CH4_F    = 1.33                               # clm5 paramfile q10ch4 (Fortran)

_dvs(ds, n) = haskey(ds, n) ? Array(ds[n]) : nothing

function _bow_latlon()
    NCDataset(FSURDAT, "r") do ds
        (Float64(Array(ds["LATIXY"])[1]), Float64(Array(ds["LONGXY"])[1]))
    end
end

# CH4 prognostic state from a Fortran restart-format dump (restart variable names).
function inject_ch4_restart!(inst, bounds, dumpfile::String)
    ch4 = inst.ch4
    nc, np = bounds.endc, bounds.endp
    nlev = CLM.varpar.nlevsoi
    NCDataset(dumpfile, "r") do ds
        col2d!(dst, name) = (a = _dvs(ds, name); a === nothing && return;
            for c in 1:min(nc, size(a, 2)), j in 1:min(nlev, size(a, 1))
                v = a[j, c]; ismissing(v) || (dst[c, j] = Float64(v))
            end)
        col1d!(dst, name) = (a = _dvs(ds, name); a === nothing && return;
            for c in 1:min(nc, length(a)); v = a[c]; ismissing(v) || (dst[c] = Float64(v)); end)
        pft1d!(dst, name) = (a = _dvs(ds, name); a === nothing && return;
            for p in 1:min(np, length(a)); v = a[p]; ismissing(v) || (dst[p] = Float64(v)); end)
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
    fill!(ch4.ch4_first_time_grc, false)
    return inst
end

# Per-step forced inputs: annsum_npp (aerenchyma tiller area, a CN restart field the
# reader does not map) and the TWS that reproduces the dumped FINUNDATED.
function force_step_inputs!(inst, bounds, beforefile::String, ch4file::String)
    NCDataset(beforefile, "r") do ds
        a = _dvs(ds, "annsum_npp")
        a === nothing || for p in 1:min(bounds.endp, length(a))
            ismissing(a[p]) || (inst.bgc_vegetation.cnveg_carbonflux_inst.annsum_npp_patch[p] = Float64(a[p]))
        end
    end
    st = inst.ch4_finundated_stream
    finF = NCDataset(ch4file, "r") do ds; Float64(_dvs(ds, "FINUNDATED")[1]); end
    tws  = (finF - st.fws_intercept_gdc[1]) / st.fws_slope_gdc[1]
    wdb  = inst.water.waterdiagnosticbulk_inst
    length(wdb.tws_grc) < bounds.endg && (wdb.tws_grc = fill(0.0, bounds.endg))
    fill!(wdb.tws_grc, tws)
    return nothing
end

function _seed_vegwp!(inst, dumpfile)
    NCDataset(dumpfile, "r") do ds
        haskey(ds, "vegwp") || return
        vw = ds["vegwp"][:, :]
        for pd in 1:size(vw, 2), seg in 1:4
            inst.canopystate.vegwp_patch[pd, seg] = Float64(vw[seg, pd])
        end
    end
end

# Fortran's hr_vr recovered from its own CH4_PROD_SAT (see header, point 3).
function recover_fortran_hr_vr(nstep::Int, params)
    nlev = CLM.varpar.nlevsoi
    prod, fin, lag = NCDataset(joinpath(DUMPDIR_SRC, "bgcdump_after_ch4_n$(nstep).nc"), "r") do ds
        (Float64.(Array(ds["CH4_PROD_SAT"])[1:nlev, 1]),
         Float64(Array(ds["FINUNDATED"])[1]), Float64(Array(ds["FINUNDATED_LAG"])[1]))
    end
    tsoi = NCDataset(joinpath(DUMPDIR_SRC, "pdump_after_hydrologydrainage_n$(nstep).nc"), "r") do ds
        a = Array(ds["T_SOISNO"])                       # (levtot, column)
        nsf = size(a, 1) - ds.dim["levgrnd"]            # the dump's own snow-layer count
        Float64.(a[nsf+1:nsf+nlev, 1])
    end
    hr = similar(prod)
    for j in 1:nlev
        f = params.f_ch4 * Q10CH4_F^((tsoi[j] - params.q10ch4base) / 10)
        lag < fin && (f *= lag / fin)
        f = min(f, 0.5)
        hr[j] = prod[j] * CLM.CATOMW / (f * params.cnscalefactor)
    end
    return hr
end

function _relerr(j, f)
    jj = vec(Float64.(j)); ff = vec(Float64.(f)); n = min(length(jj), length(ff))
    s = 0.0; for i in 1:n; isfinite(ff[i]) && (s = max(s, abs(ff[i]))); end
    s <= 0 && return (0.0, 0.0)
    w = 0.0
    for i in 1:n
        (isfinite(jj[i]) && isfinite(ff[i])) || continue
        w = max(w, abs(jj[i] - ff[i]) / s)
    end
    return (w, s)
end

const SRC_LEV = [
    ("CH4_PROD_SAT",   i -> i.ch4.ch4_prod_depth_sat_col),
    ("CH4_OXID_SAT",   i -> i.ch4.ch4_oxid_depth_sat_col),
    ("CH4_OXID_UNSAT", i -> i.ch4.ch4_oxid_depth_unsat_col),
    ("CH4_AERE_SAT",   i -> i.ch4.ch4_aere_depth_sat_col),
    ("CH4_AERE_UNSAT", i -> i.ch4.ch4_aere_depth_unsat_col),
    ("CH4_EBUL_SAT",   i -> i.ch4.ch4_ebul_depth_sat_col),
    ("CH4_TRAN_SAT",   i -> i.ch4.ch4_tran_depth_sat_col),
    ("CH4_TRAN_UNSAT", i -> i.ch4.ch4_tran_depth_unsat_col),
    ("CONC_CH4_SAT",   i -> i.ch4.conc_ch4_sat_col),
    ("CONC_CH4_UNSAT", i -> i.ch4.conc_ch4_unsat_col),
    ("CONC_O2_SAT",    i -> i.ch4.conc_o2_sat_col),
    ("CONC_O2_UNSAT",  i -> i.ch4.conc_o2_unsat_col),
    ("O2STRESS_SAT",   i -> i.ch4.o2stress_sat_col),
    ("CH4STRESS_SAT",  i -> i.ch4.ch4stress_sat_col),
    ("O2_DECOMP_SAT",  i -> i.ch4.o2_decomp_depth_sat_col),
    ("O2_DECOMP_UNSAT",i -> i.ch4.o2_decomp_depth_unsat_col),
]
# O2-demand components and their inputs. Only the "bgcdump O2-demand" rebuild
# (scripts/validation/fortran_sourcemods/bgcdumpMod.F90, O2 block) writes these;
# older references simply lack them and the rows are skipped.
function _o2_rr_demand(inst, nlev)
    nlevsno = CLM.varpar.nlevsno
    pch = inst.patch; ss = inst.soilstate
    rr = inst.bgc_vegetation.cnveg_carbonflux_inst.rr_patch
    out = zeros(nlev)
    for p in eachindex(pch.column)
        pch.column[p] == 1 || continue
        (pch.wtcol[p] > 0 && pch.itype[p] != CLM.noveg) || continue
        for j in 1:nlev
            out[j] += rr[p] * ss.crootfr_patch[p, j] * pch.wtcol[p]
        end
    end
    for j in 1:nlev
        out[j] /= CLM.CATOMW * inst.column.dz[1, nlevsno + j]
    end
    return out
end
const SRC_O2DEM = [
    ("O2_DEMAND_RR",   (i, n) -> _o2_rr_demand(i, n)),
    ("O2_DEMAND_NIT",  (i, n) -> i.soilbiogeochem_nitrogenflux.pot_f_nit_vr_col[1, 1:n] .* (2 / 14)),
    ("POT_F_NIT_VR",   (i, n) -> i.soilbiogeochem_nitrogenflux.pot_f_nit_vr_col[1, 1:n]),
    ("ANAEROBIC_FRAC", (i, n) -> i.soilbiogeochem_nitrogenflux.anaerobic_frac_col[1, 1:n]),
    ("K_NITR_VR",      (i, n) -> i.soilbiogeochem_nitrogenflux.k_nitr_vr_col[1, 1:n]),
    ("DIFFUS",         (i, n) -> i.soilbiogeochem_nitrogenflux.diffus_col[1, 1:n]),
    ("R_PSI",          (i, n) -> i.soilbiogeochem_nitrogenflux.r_psi_col[1, 1:n]),
    ("O_SCALAR",       (i, n) -> i.soilbiogeochem_carbonflux.o_scalar_col[1, 1:n]),
    ("HR_VR",          (i, n) -> i.soilbiogeochem_carbonflux.hr_vr_col[1, 1:n]),
    ("O2_OXID_SAT",    (i, n) -> i.ch4.o2_oxid_depth_sat_col[1, 1:n]),
    ("O2_OXID_UNSAT",  (i, n) -> i.ch4.o2_oxid_depth_unsat_col[1, 1:n]),
    ("O2_AERE_SAT",    (i, n) -> i.ch4.o2_aere_depth_sat_col[1, 1:n]),
    ("O2_AERE_UNSAT",  (i, n) -> i.ch4.o2_aere_depth_unsat_col[1, 1:n]),
    ("RR",             (i, n) -> i.bgc_vegetation.cnveg_carbonflux_inst.rr_patch),   # patch vector
    ("FROOT_MR",       (i, n) -> i.bgc_vegetation.cnveg_carbonflux_inst.froot_mr_patch),
    ("LIVECROOT_MR",   (i, n) -> i.bgc_vegetation.cnveg_carbonflux_inst.livecroot_mr_patch),
]
const SRC_COL = [
    ("CH4_SURF_DIFF_SAT",   i -> i.ch4.ch4_surf_diff_sat_col),
    ("CH4_SURF_DIFF_UNSAT", i -> i.ch4.ch4_surf_diff_unsat_col),
    ("CH4_SURF_AERE_SAT",   i -> i.ch4.ch4_surf_aere_sat_col),
    ("CH4_SURF_AERE_UNSAT", i -> i.ch4.ch4_surf_aere_unsat_col),
    ("CH4_SURF_EBUL_SAT",   i -> i.ch4.ch4_surf_ebul_sat_col),
    ("CH4_SURF_FLUX_TOT",   i -> i.ch4.ch4_surf_flux_tot_col),
    ("TOTCOLCH4",           i -> i.ch4.totcolch4_col),
    ("ZWT_CH4_UNSAT",       i -> i.ch4.zwt_ch4_unsat_col),
]

function run_source(; nsteps::Int = NLASTS - N0S, mode::Symbol = :locked, report_every::Int = 0)
    n1 = N0S + 1
    nsteps = min(nsteps, NLASTS - N0S)
    report_every = report_every > 0 ? report_every : (nsteps > 48 ? 12 : 1)
    (inst, bounds, filt, tm) = build_bow_inst(; dtime = 3600, start_date = DATE0S + Hour(n1 - N0S),
        use_cn = true, use_luna = true, use_lch4 = true, fndep = FNDEP_SRC)
    if isempty(inst.photosyns.vcmx25_z_patch)
        inst.photosyns.vcmx25_z_patch = fill(30.0, bounds.endp, CLM.NLEVCAN)
        inst.photosyns.jmx25_z_patch  = fill(60.0, bounds.endp, CLM.NLEVCAN)
    end
    # TWS_inversion finundation with the real stream coefficients at Bow's cell.
    inst.ch4_varcon.finundation_mtd = CLM.FINUNDATION_MTD_TWS_INVERSION
    (blat, blon) = _bow_latlon()
    CLM.read_ch4_finundated_stream!(inst.ch4_finundated_stream, INV_FILE_SRC, [blat], [blon];
                                    method = CLM.FINUNDATION_MTD_TWS_INVERSION)

    first = joinpath(DUMPDIR_SRC, "pdump_before_step_n$(n1).nc")
    inject_dump!(inst, bounds, first)
    inject_ch4_restart!(inst, bounds, first)
    _seed_vegwp!(inst, first)

    config = CLM.CLMDriverConfig(use_cn = true, use_aquifer_layer = false, use_hydrstress = true,
                                 use_luna = true, use_lch4 = true)
    filt_ia = CLM.clump_filter_inactive_and_active
    ng, nc = bounds.endg, bounds.endc
    fr = CLM.ForcingReader(); CLM.forcing_reader_init!(fr, FFORC_SRC)
    tf = replace(FFORC_SRC, r"clmforc\.[^/]*\.nc$" => "topo_forcing.nc")
    if isfile(tf)
        NCDataset(tf, "r") do dt
            if haskey(dt, "TOPO")
                ft = Float64(dt["TOPO"][1])
                for g in 1:ng; inst.atm2lnd.forc_topo_grc[g] = ft; end
                for c in 1:nc; inst.topo.topo_col[c] = ft; end
            end
        end
    end
    let cd0 = CLM.get_curr_calday(tm)
        (d0, _) = CLM.compute_orbital(cd0); (dm1, _) = CLM.compute_orbital(cd0 - 3600.0 / CLM.SECSPDAY)
        CLM.init_daylength!(inst.gridcell, d0, dm1, CLM.ORB_OBLIQR_DEFAULT, 1:ng)
    end

    nlev = CLM.varpar.nlevsoi
    worst = Dict{String,Float64}(); scale = Dict{String,Float64}()
    hrworst = 0.0; hrsum_rel = 0.0
    nonhr_worst = Ref(0.0)
    traj = Tuple{Int,Float64,Float64,Float64,Float64}[]   # n, TOTCOLCH4 J/F, FLUX J/F
    println("step  | hr_vr J/F  | CH4_PROD_SAT | CONC_CH4_SAT | TOTCOLCH4 J / F            | CH4_SURF_FLUX_TOT J / F")
    for i in 0:(nsteps - 1)
        n = n1 + i
        before = joinpath(DUMPDIR_SRC, "pdump_before_step_n$(n).nc")
        ch4f   = joinpath(DUMPDIR_SRC, "bgcdump_after_ch4_n$(n).nc")
        if mode === :locked && i > 0
            inject_dump!(inst, bounds, before)       # CN + biogeophysics locked to Fortran
            _seed_vegwp!(inst, before)
        end
        force_step_inputs!(inst, bounds, before, ch4f)

        CLM.advance_timestep!(tm)
        calday = CLM.get_curr_calday(tm); (declin, _) = CLM.compute_orbital(calday)
        CLM.read_forcing_step!(fr, inst.atm2lnd, FORCE0_SRC + Hour(n - N0S), ng, nc)
        CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
        (yr, mon, d, tod) = CLM.get_curr_date(tm)
        CLM.clm_drv!(config, inst, filt, filt_ia, bounds, true, calday, declin, declin,
            CLM.ORB_OBLIQR_DEFAULT, false, false, "", false;
            # nstep = the harness clock (1, 2, ...), as in run_one_parity_step!: the
            # balance-check hard-error gate (DAnstep > skip_steps) stays off, as it does
            # for every injection harness. The injected albedo/radiation state is a
            # restart snapshot, not a self-consistent Julia history.
            nstep = tm.nstep, is_first_step = false,
            is_beg_curr_day = CLM.is_beg_curr_day(tm), is_end_curr_day = CLM.is_end_curr_day(tm),
            is_beg_curr_year = CLM.is_beg_curr_year(tm), dtime = 3600.0, mon = mon, day = d,
            jday = dayofyear(tm.current_date), secs = tod, year = yr,
            photosyns = inst.photosyns)

        # --- HR attribution: Julia hr_vr vs Fortran hr_vr recovered from CH4_PROD_SAT ---
        hrF = recover_fortran_hr_vr(n, inst.ch4_params)
        hrJ = inst.soilbiogeochem_carbonflux.hr_vr_col[1, 1:nlev]
        (ehr, _) = _relerr(hrJ, hrF)
        hrworst = max(hrworst, ehr)
        hrsum_rel = max(hrsum_rel, abs(sum(hrJ) - sum(hrF)) / abs(sum(hrF)))

        rowerr = Dict{String,Float64}()
        # --- O2-demand attribution: o2_decomp = hr_vr/catomw + rr_vr/(catomw dz)
        #     + pot_f_nit_vr*2/14 (anoxia, o_scalar=1). Subtracting each side's own
        #     HR term leaves the root-respiration + nitrification O2 demand. ---
        NCDataset(ch4f, "r") do ds
            o2F = Float64.(ds["O2_DECOMP_SAT"][1:nlev, 1])
            o2J = inst.ch4.o2_decomp_depth_sat_col[1, 1:nlev]
            (e, _) = _relerr(o2J .- hrJ ./ CLM.CATOMW, o2F .- hrF ./ CLM.CATOMW)
            nonhr_worst[] = max(nonhr_worst[], e)
            if i == 0
                nit = inst.soilbiogeochem_nitrogenflux.pot_f_nit_vr_col[1, 1:nlev] .* (2 / 14)
                @printf("  O2 demand split (Julia, n%d, column max): HR %.3e  nitrif %.3e  roots %.3e mol/m3/s\n",
                        n, maximum(hrJ ./ CLM.CATOMW), maximum(nit),
                        maximum(o2J .- hrJ ./ CLM.CATOMW .- nit))
            end
        end
        NCDataset(ch4f, "r") do ds
            for (nm, g) in SRC_LEV
                f = _dvs(ds, nm); f === nothing && continue
                (e, s) = _relerr(g(inst)[1, 1:nlev], Float64.(f[1:nlev, 1]))
                worst[nm] = max(get(worst, nm, 0.0), e); scale[nm] = max(get(scale, nm, 0.0), s)
                rowerr[nm] = e
            end
            for (nm, g) in SRC_O2DEM
                f = _dvs(ds, nm); f === nothing && continue
                fv = ndims(f) == 1 ? Float64.(f) : Float64.(f[1:nlev, 1])
                jv = g(inst, nlev); m = min(length(jv), length(fv))
                (e, s) = _relerr(jv[1:m], fv[1:m])
                worst[nm] = max(get(worst, nm, 0.0), e); scale[nm] = max(get(scale, nm, 0.0), s)
                rowerr[nm] = e
            end
            for (nm, g) in SRC_COL
                f = _dvs(ds, nm); f === nothing && continue
                (e, s) = _relerr([g(inst)[1]], [Float64(f[1])])
                worst[nm] = max(get(worst, nm, 0.0), e); scale[nm] = max(get(scale, nm, 0.0), s)
                rowerr[nm] = e
            end
            push!(traj, (n, inst.ch4.totcolch4_col[1], Float64(ds["TOTCOLCH4"][1]),
                         inst.ch4.ch4_surf_flux_tot_col[1], Float64(ds["CH4_SURF_FLUX_TOT"][1])))
        end
        if (i % report_every == 0) || i == nsteps - 1
            t = traj[end]
            @printf("n%d | %9.2e | %12.2e | %12.2e | %.6e / %.6e | %.5e / %.5e\n",
                    n, ehr, rowerr["CH4_PROD_SAT"], rowerr["CONC_CH4_SAT"], t[2], t[3], t[4], t[5])
            if haskey(rowerr, "O2_DEMAND_RR")
                @printf("      O2: conc_sat %.2e  decomp_sat %.2e  rr %.2e  nit %.2e  anfr %.2e  o_scalar %.2e  oxid_sat %.2e  aere_sat %.2e  RR %.2e\n",
                        rowerr["CONC_O2_SAT"], rowerr["O2_DECOMP_SAT"], rowerr["O2_DEMAND_RR"],
                        rowerr["O2_DEMAND_NIT"], rowerr["ANAEROBIC_FRAC"], rowerr["O_SCALAR"],
                        rowerr["O2_OXID_SAT"], rowerr["O2_AERE_SAT"], rowerr["RR"])
            end
            flush(stdout)
        end
    end
    CLM.forcing_reader_close!(fr)

    println()
    @printf("hr_vr (Julia vs Fortran-recovered): worst profile rel.err %.3e, column-sum %.3e\n",
            hrworst, hrsum_rel)
    @printf("O2 demand minus HR (roots + nitrification), Julia vs Fortran: worst rel.err %.3e\n",
            nonhr_worst[])
    @printf("%-22s %14s %14s\n", "FIELD", "worst rel.err", "|F| scale")
    println("-"^54)
    for (nm, _) in vcat(SRC_LEV, SRC_O2DEM, SRC_COL)
        haskey(worst, nm) || continue
        @printf("%-22s %14.3e %14.3e%s\n", nm, worst[nm], scale[nm], scale[nm] == 0 ? "  (F≡0)" : "")
    end
    return (; worst, scale, hrworst, hrsum_rel, nonhr_worst = nonhr_worst[], traj)
end

if abspath(PROGRAM_FILE) == @__FILE__
    ns   = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : NLASTS - N0S
    mode = length(ARGS) >= 2 ? Symbol(ARGS[2]) : :locked
    println("="^78)
    println("CH4-AS-A-SOURCE — Bow TWS_inversion, CN locked to Fortran, CH4 free-running")
    println("  reference : ", DUMPDIR_SRC)
    println("  window    : n", N0S + 1, " .. n", N0S + ns, "  mode=", mode)
    println("="^78)
    run_source(; nsteps = ns, mode = mode)
end
