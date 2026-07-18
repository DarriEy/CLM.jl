# ==========================================================================
# CH4 SATURATED-COLUMN TRANSPORT ISOLATION TEST
#
# Purpose: prove whether the CONC_CH4_SAT single-step collapse (PR #235) is
# driven by the reaction-diffusion TRANSPORT solve (ch4_tran!, sat=1) or by an
# upstream source term.
#
# Method (transport is the ONLY variable):
#   1. Build a Bow (or MerBleue) instance; inject the Fortran before-step
#      biogeophysics state + the ch4 prognostics (conc_ch4_sat, conc_o2_sat,
#      o2stress_sat, o2_decomp_depth_sat, accumulators) — exactly Fortran's start.
#   2. Inject Fortran's per-layer POST-competition source/sink terms for the sat
#      column from the after_ch4 dump (CH4_PROD/OXID/AERE/EBUL_SAT), so the
#      production/oxidation confound is removed.
#   3. Zero jwt for the sat column (as the driver does) and call CLM.ch4_tran!
#      directly for sat=1 — nothing else runs.
#   4. Diff Julia's step-(n+1) conc_ch4_sat / conc_o2_sat against Fortran's
#      after_ch4 dump, per layer.
#
# If Julia collapses while Fortran stays non-zero -> transport solve is the bug.
#
# Usage: julia --project=. scripts/ch4_tran_sat_injection.jl [nstep]
# ==========================================================================

include(joinpath(@__DIR__, "fortran_parity_common.jl"))
using Printf

const DUMPDIR_CH4 = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_firech4"
const N0    = 4720
const DATE0 = DateTime(2003, 7, 16, 16)
const FNDEP = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/ndepdata/" *
              "fndep_clm_UNIFORM1e-8_0.9x1.25_yr2000_CLMjl-parity.nc"

_dv(ds, n) = haskey(ds, n) ? Array(ds[n]) : nothing

# --- ch4 prognostic injection (before-step state), same as fortran_parity_ch4.jl ---
function inject_ch4_before!(inst, bounds, dumpfile::String)
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
            for c in 1:min(nc, length(a)); v = a[c]; ismissing(v) || (dst[c] = Float64(v)); end
        end
        pft1d!(dst, name) = begin
            a = _dv(ds, name); a === nothing && return
            for p in 1:min(np, length(a)); v = a[p]; ismissing(v) || (dst[p] = Float64(v)); end
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
    fill!(ch4.ch4_first_time_grc, false)
    return inst
end

# --- inject Fortran's post-competition sat source/sink depth terms ---
function inject_sat_source!(inst, bounds, afterfile::String)
    ch4 = inst.ch4
    nc = bounds.endc
    nlev = CLM.varpar.nlevsoi
    NCDataset(afterfile, "r") do ds
        put!(dst, name) = begin
            a = _dv(ds, name); a === nothing && return
            for c in 1:min(nc, size(a, 2)), j in 1:min(nlev, size(a, 1))
                v = a[j, c]; ismissing(v) || (dst[c, j] = Float64(v))
            end
        end
        put!(ch4.ch4_prod_depth_sat_col, "CH4_PROD_SAT")
        put!(ch4.ch4_oxid_depth_sat_col, "CH4_OXID_SAT")
        put!(ch4.ch4_aere_depth_sat_col, "CH4_AERE_SAT")
        put!(ch4.ch4_ebul_depth_sat_col, "CH4_EBUL_SAT")
    end
    return inst
end

function run_isolation(nstep::Int)
    step_date = DATE0 + Hour(nstep - N0)
    (inst, bounds, filt, tm) = build_bow_inst(; dtime=3600, start_date=step_date,
        use_cn=true, use_luna=true, use_lch4=true,
        fndep=FNDEP, cnfire_method=:nofire)
    if isempty(inst.photosyns.vcmx25_z_patch)
        inst.photosyns.vcmx25_z_patch = fill(30.0, bounds.endp, CLM.NLEVCAN)
        inst.photosyns.jmx25_z_patch  = fill(60.0, bounds.endp, CLM.NLEVCAN)
    end
    beforefile = joinpath(DUMPDIR_CH4, "pdump_before_step_n$(nstep).nc")
    afterfile  = joinpath(DUMPDIR_CH4, "bgcdump_after_ch4_n$(nstep).nc")
    inject_dump!(inst, bounds, beforefile)          # 16-field biogeophys oracle
    inject_ch4_before!(inst, bounds, beforefile)     # ch4 prognostics (before-step)
    inject_sat_source!(inst, bounds, afterfile)      # Fortran exact sat source terms

    ch4 = inst.ch4
    col = inst.column
    ss  = inst.soilstate
    temp = inst.temperature
    wsb = inst.water.waterstatebulk_inst
    wdb = inst.water.waterdiagnosticbulk_inst
    nc = bounds.endc
    nlev = CLM.varpar.nlevsoi
    nlevsno = CLM.varpar.nlevsno
    _sl = (nlevsno + 1):(nlevsno + nlev)              # snow-offset soil slice
    _dl = 1:nlev                                       # soil slice
    _zl = (nlevsno + 1):(nlevsno + nlev + 1)           # zi slice

    # jwt = 0 for the saturated column (matches the driver's zero-kernel)
    jwt = zeros(Int, nc)

    # snapshot injected source terms to detect competition re-application
    prod0 = copy(ch4.ch4_prod_depth_sat_col[1, 1:nlev])
    oxid0 = copy(ch4.ch4_oxid_depth_sat_col[1, 1:nlev])
    aere0 = copy(ch4.ch4_aere_depth_sat_col[1, 1:nlev])
    ebul0 = copy(ch4.ch4_ebul_depth_sat_col[1, 1:nlev])
    conc0 = copy(ch4.conc_ch4_sat_col[1, 1:nlev])

    mask_soil = collect(Bool, filt.soilc)

    CLM.ch4_tran!(ch4, inst.ch4_params, inst.ch4_varcon, mask_soil, Array(col.gridcell),
        Matrix(ss.watsat_col[:, _dl]), Matrix(wsb.ws.h2osoi_vol_col[:, _dl]),
        Matrix(wsb.ws.h2osoi_liq_col[:, _sl]), Matrix(wsb.ws.h2osoi_ice_col[:, _sl]),
        wsb.ws.h2osfc_col,
        Matrix(ss.bsw_col[:, _dl]), Matrix(ss.cellorg_col[:, _dl]),
        Matrix(temp.t_soisno_col[:, _sl]), temp.t_grnd_col, temp.t_h2osfc_col,
        wdb.frac_h2osfc_col, wdb.snow_depth_col, Array(col.snl),
        Matrix(col.z[:, _sl]), Matrix(col.dz[:, _sl]), Matrix(col.zi[:, _zl]),
        jwt, 1, false, nlev, nlevsno, 3600.0, 130.0)

    # Fortran after-step reference
    concF = NCDataset(afterfile, "r") do ds
        a = _dv(ds, "CONC_CH4_SAT"); Float64.(a[1:nlev, 1])
    end
    o2F = NCDataset(afterfile, "r") do ds
        a = _dv(ds, "CONC_O2_SAT"); Float64.(a[1:nlev, 1])
    end

    concJ = ch4.conc_ch4_sat_col[1, 1:nlev]
    o2J   = ch4.conc_o2_sat_col[1, 1:nlev]

    println("="^78)
    println("CH4 sat-column TRANSPORT ISOLATION — nstep $nstep (Bow)")
    println("="^78)
    # competition no-op check
    dprod = maximum(abs.(ch4.ch4_prod_depth_sat_col[1,1:nlev] .- prod0))
    doxid = maximum(abs.(ch4.ch4_oxid_depth_sat_col[1,1:nlev] .- oxid0))
    daere = maximum(abs.(ch4.ch4_aere_depth_sat_col[1,1:nlev] .- aere0))
    debul = maximum(abs.(ch4.ch4_ebul_depth_sat_col[1,1:nlev] .- ebul0))
    @printf("competition delta on injected source (max abs): prod %.2e oxid %.2e aere %.2e ebul %.2e\n",
            dprod, doxid, daere, debul)
    println("(if ~0, the injected Fortran source passed through unchanged => transport isolated)\n")

    scale = maximum(abs.(concF))
    @printf("%4s %14s %14s %14s   %14s %14s\n", "lev", "conc_ch4 before", "conc_ch4 Julia", "conc_ch4 Fortran", "o2 Julia", "o2 Fortran")
    for j in 1:nlev
        @printf("%4d %14.6e %14.6e %14.6e   %14.6e %14.6e\n", j, conc0[j], concJ[j], concF[j], o2J[j], o2F[j])
    end
    worst = 0.0
    for j in 1:nlev
        (isfinite(concJ[j]) && isfinite(concF[j])) || continue
        worst = max(worst, abs(concJ[j] - concF[j]) / max(scale, eps()))
    end
    println("-"^78)
    @printf("CONC_CH4_SAT profile-scaled worst rel.err = %.4e   (Fortran |scale| = %.4e)\n", worst, scale)
    println("  ^ CLEAN transport test: all 4 CH4 source terms injected from Fortran.")
    o2scale = maximum(abs.(o2F))
    worst2 = 0.0
    for j in 1:nlev
        (isfinite(o2J[j]) && isfinite(o2F[j])) || continue
        worst2 = max(worst2, abs(o2J[j] - o2F[j]) / max(o2scale, eps()))
    end
    @printf("CONC_O2_SAT  profile-scaled worst rel.err = %.4e   (Fortran |scale| = %.4e)\n", worst2, o2scale)
    println("  ^ NOT a clean isolation: o2_oxid_depth / o2_aere_depth are not dumped as")
    println("    separate fields, so the O2 source term is incomplete here. Use the full-")
    println("    step harness (fortran_parity_ch4.jl) for the O2 verdict.")
    return worst
end

nstep = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 4721
run_isolation(nstep)
