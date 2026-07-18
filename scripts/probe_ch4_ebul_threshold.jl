# ==========================================================================
# CH4 EBULLITION threshold probe (docs/CH4_FIRE_PARITY.md §13 — the A1 residual).
#
# Ebullition fires only where the saturated dissolved-CH4 concentration crosses
# the bubble threshold `vgc > vgc_max*bubble_f = 0.0855`. Bow's PHYSICAL peak
# CONC_CH4_SAT (~0.029 mol/m^3) gives vgc ~= 0.045, so `_meth_ebul_kernel!` never
# fires at any on-disk site and the kernel's firing regime is numerically untested.
#
# This probe builds the Bow TWS_inversion single-step instance (finundated>0,
# warm CN), injects the Fortran before-step CH4 state, SCALES conc_ch4 by a factor
# (the "synthetic shared input" technique of the peatf §12 / lnfm §2 references),
# runs the step, and reports per-layer vgc + ch4_ebul. It confirms — with the
# MODEL's real watsat/pressure — how large a conc bump crosses the threshold, so
# the Fortran reference restart is scaled correctly before the box run.
#
# Usage: julia +1.12 --project=. scripts/probe_ch4_ebul_threshold.jl [factor] [nstep]
# ==========================================================================

include(joinpath(@__DIR__, "fortran_parity_common.jl"))
using Printf

const DUMPDIR_CH4 = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bow_ref_ch4tws"
const INV_FILE    = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/paramdata/" *
                    "finundated_inversiondata_0.9x1.25_c170706.nc"
const N0    = 4720
const DATE0 = DateTime(2003, 7, 16, 16)
const FNDEP = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/ndepdata/" *
              "fndep_clm_UNIFORM1e-8_0.9x1.25_yr2000_CLMjl-parity.nc"
const FFORCING2002 = joinpath(dirname(FFORCING), "clmforc.2002.nc")
const FORCE_DATE0  = DateTime(2002, 7, 16, 15)

_dv(ds, n) = haskey(ds, n) ? Array(ds[n]) : nothing

function bow_lat_lon()
    NCDataset(FSURDAT, "r") do ds
        (Float64(Array(ds["LATIXY"])[1]), Float64(Array(ds["LONGXY"])[1]))
    end
end

# Inject the Fortran CH4 before-step state + TWS_inversion config, then SCALE conc.
function inject_and_scale!(inst, bounds, dumpfile::String, scale::Float64)
    ch4 = inst.ch4
    nc, np = bounds.endc, bounds.endp
    nlev = CLM.varpar.nlevsoi
    NCDataset(dumpfile, "r") do ds
        col2d!(dst, name; f=1.0) = begin
            a = _dv(ds, name); a === nothing && return
            for c in 1:min(nc, size(a, 2)), j in 1:min(nlev, size(a, 1))
                v = a[j, c]; ismissing(v) || (dst[c, j] = Float64(v) * f)
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
        col2d!(ch4.conc_ch4_sat_col,   "CONC_CH4_SAT"; f=scale)
        col2d!(ch4.conc_ch4_unsat_col, "CONC_CH4_UNSAT"; f=scale)
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
        pft1d!(inst.bgc_vegetation.cnveg_carbonflux_inst.annsum_npp_patch, "annsum_npp")

        inst.ch4_varcon.finundation_mtd = CLM.FINUNDATION_MTD_TWS_INVERSION
        (blat, blon) = bow_lat_lon()
        CLM.read_ch4_finundated_stream!(inst.ch4_finundated_stream, INV_FILE,
                                        [blat], [blon]; method = CLM.FINUNDATION_MTD_TWS_INVERSION)
        slope = inst.ch4_finundated_stream.fws_slope_gdc[1]
        inter = inst.ch4_finundated_stream.fws_intercept_gdc[1]
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

function main()
    scale = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 5.0
    nstep = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : N0
    println("="^72)
    @printf("CH4 EBULLITION threshold probe — Bow TWS_inversion, nstep %d, conc x%.1f\n", nstep, scale)
    println("  threshold: vgc > vgc_max*bubble_f = 0.15*0.57 = 0.0855")
    println("="^72)
    (inst, bounds) = run_one_parity_step!(nstep;
        use_cn = true, use_lch4 = true, use_hydrstress = true,
        dumpdir = DUMPDIR_CH4,
        step_date = DATE0 + Hour(nstep - N0),
        forcing_file = FFORCING2002, forcing_date = FORCE_DATE0 + Hour(nstep - N0),
        fndep = FNDEP,
        pre_step_hook = (i, b, df) -> inject_and_scale!(i, b, df, scale))
    nlev  = CLM.varpar.nlevsoi
    ebulS = inst.ch4.ch4_ebul_depth_sat_col
    concS = inst.ch4.conc_ch4_sat_col
    ebulU = inst.ch4.ch4_ebul_depth_unsat_col
    @printf("%-4s %16s %16s\n", "lev", "conc_ch4_sat(post)", "ch4_ebul_sat")
    nfire = 0
    for j in 1:nlev
        e = ebulS[1, j]
        e > 0 && (nfire += 1)
        (e > 0 || j <= 8) && @printf("%-4d %16.6e %16.6e %s\n", j, concS[1, j], e, e > 0 ? "<== EBUL" : "")
    end
    @printf("\nSAT ebul total (sum layers) = %.6e   fired in %d layers\n",
            sum(ebulS[1, 1:nlev]), nfire)
    @printf("UNSAT ebul total            = %.6e\n", sum(ebulU[1, 1:nlev]))
    @printf("CH4_EBUL_TOTAL_SAT col      = %.6e\n", inst.ch4.ch4_ebul_total_sat_col[1])
    if nfire == 0
        println("\n[!] No ebullition — increase the scale factor.")
    else
        println("\n[OK] Ebullition kernel FIRES at this scale — safe to scale the Fortran restart identically.")
    end
end

main()
