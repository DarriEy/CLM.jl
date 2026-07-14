# ==========================================================================
# Surface-runoff (qflx_surf) Fortran parity — the whole surface-water chain.
#
# WHY THIS EXISTS
#   PR #221 validated `n_leaching!` and found the leaching OPERATOR faithful
#   (calibrate one scalar on soil level 1, predict levels 2..20 to <= 4.4e-07)
#   while the leaching TOTALS disagreed — because the WATER FLUX driving them,
#   `qflx_surf`, disagreed. This harness scores that water flux directly.
#
# THE FORTRAN REFERENCE
#   run dir : SYMFLUENCE_data/clm_ref_water_summer
#   file    : Bow_at_Banff_lumped.clm2.h1.2202-07-16-61200.nc
#   200 contiguous steps, hist_nhtfrq=1 + hist_avgflag_pertape='I', i.e. an
#   INSTANTANEOUS per-step dump of the surface-water chain -- not a time mean.
#   Record i  <->  nstep = 1757872 + i   (record 1 ends at 2202-07-16 17:00,
#   which is the step that STARTS at 16:00 = nstep 1757873, the summer window's
#   first step).
#
#   Eight of these fields are NOT stock CTSM history fields (QSATEXCS, QINFLEXC,
#   QINFLEXCS, QINSOIL, QINSOILL, QTS2SFC, QINH2OSFC, QSFCDRAIN); they were added
#   as hist_addfld1d registrations in a SourceMods copy of WaterFluxBulkType.F90
#   so the runoff partition could be bisected term-by-term. See
#   scripts/validation/fortran_pdump/README.md.
#
# BOUNDARY DISCIPLINE
#   The h1 record is written at the END of the step, so it is the SAME-step value
#   of everything computed during that step. The `before_step` pdump (previous
#   step's state) is used ONLY as the injected IC -- never as a comparison target.
#
# Usage:
#   julia --project=. scripts/fortran_parity_qflx_surf.jl [nprobe] [h2osfcflag]
# ==========================================================================

include(joinpath(@__DIR__, "fortran_parity_common.jl"))
using Printf

const DUMPDIR_SUMMER = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_ndep_summer"
const H1FILE = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_ref_water_summer/" *
               "Bow_at_Banff_lumped.clm2.h1.2202-07-16-61200.nc"
const FNDEP = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/ndepdata/" *
              "fndep_clm_UNIFORM1e-8_0.9x1.25_yr2000_CLMjl-parity.nc"

const N0   = 1757873                       # nstep of h1 record 1
# FORTRAN model date of N0. Fed through parity_forcing(), which applies the datm
# cycling rule to select data year 2002 (NOT 2003 — see fortran_parity_common.jl).
const DATE0 = DateTime(2202, 7, 16, 16)

# Fortran h1 record index for a given nstep
_rec(nstep::Int) = nstep - N0 + 1

"""Load the whole h1 file once: Dict(varname => Vector over records)."""
function load_h1(path::String)
    ds = NCDataset(path, "r")
    out = Dict{String, Vector{Float64}}()
    for v in ("QOVER", "QSATEXCS", "QINFLEXC", "QINFLEXCS", "QH2OSFC", "FH2OSFC",
              "H2OSFC", "FSAT", "FCOV", "QINSOIL", "QINSOILL", "QTOPSOIL", "QINFL",
              "QDRAI", "QRGWL", "QRUNOFF", "QTS2SFC", "QINH2OSFC", "QSFCDRAIN",
              "ZWT", "FROST_TABLE", "RAIN", "QSNOMELT", "FSNO")
        haskey(ds, v) || continue
        out[v] = vec(Float64.(Array(ds[v][:, :])))
    end
    close(ds)
    return out
end

# Julia counterpart of each Fortran history field (column 1).
function julia_chain(inst)
    wfb = inst.water.waterfluxbulk_inst
    wf  = wfb.wf
    wdb = inst.water.waterdiagnosticbulk_inst
    ws  = inst.water.waterstatebulk_inst.ws
    sh  = inst.soilhydrology
    ser = inst.sat_excess_runoff
    Dict{String, Float64}(
        "QOVER"     => wf.qflx_surf_col[1],
        "QSATEXCS"  => wfb.qflx_sat_excess_surf_col[1],
        "QINFLEXC"  => wfb.qflx_infl_excess_col[1],
        "QINFLEXCS" => wfb.qflx_infl_excess_surf_col[1],
        "QH2OSFC"   => wfb.qflx_h2osfc_surf_col[1],
        "FH2OSFC"   => wdb.frac_h2osfc_col[1],
        "H2OSFC"    => ws.h2osfc_col[1],
        "FSAT"      => ser.fsat_col[1],
        "FCOV"      => ser.fcov_col[1],
        "QINSOIL"   => wfb.qflx_in_soil_col[1],
        "QINSOILL"  => wfb.qflx_in_soil_limited_col[1],
        "QTOPSOIL"  => wf.qflx_top_soil_col[1],
        "QINFL"     => wf.qflx_infl_col[1],
        "QDRAI"     => wf.qflx_drain_col[1],
        "QRGWL"     => wf.qflx_qrgwl_col[1],
        "QTS2SFC"   => wfb.qflx_top_soil_to_h2osfc_col[1],
        "QINH2OSFC" => wfb.qflx_in_h2osfc_col[1],
        "QSFCDRAIN" => wfb.qflx_h2osfc_drain_col[1],
        "ZWT"       => sh.zwt_col[1],
        # qflx_rain_plus_snomelt has no stock history field; expose it for the
        # fsat-vs-forcing bisection (QSATEXCS = FSAT * this, exactly).
        "QRAINSNOM" => wf.qflx_rain_plus_snomelt_col[1],
        # Upstream of qflx_rain_plus_snomelt: the FORCING (forc_rain) and the
        # post-interception throughfall (qflx_liq_grnd) + snowmelt. Fortran:
        #   nosnowc: qflx_rain_plus_snomelt = qflx_liq_grnd + qflx_snomelt
        #   snowc  : qflx_rain_plus_snomelt = perc_bottom + (1-frac_sno_eff)*qflx_liq_grnd
        "RAIN"      => inst.atm2lnd.forc_rain_downscaled_col[1],
        "QSNOMELT"  => wf.qflx_snomelt_col[1],
        "FSNO"      => wdb.frac_sno_eff_col[1],
        "QLIQGRND"  => wf.qflx_liq_grnd_col[1],
        "SNL"       => Float64(inst.column.snl[1]),
    )
end

function main()
    nprobe = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 12
    h2oflag = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1

    F = load_h1(H1FILE)
    nrec = length(F["QOVER"])
    # probe within BOTH the h1 range and the pdump (before_step) range
    lo, hi = N0 + 1, min(N0 + nrec - 2, 1758060)
    steps = unique(round.(Int, range(lo, hi; length = nprobe)))

    println("="^118)
    println("qflx_surf FORTRAN PARITY — Bow at Banff — summer window — h2osfcflag=$h2oflag")
    println("Fortran ref: per-step INSTANTANEOUS history (clm2.h1), $nrec records")
    println("="^118)
    @printf("\n%-9s | %-11s %-11s %-9s | %-11s %-11s %-9s | %-11s %-11s\n",
            "nstep", "QOVER J", "QOVER F", "ratio", "FSAT J", "FSAT F", "ratio",
            "QRAINSNOM J", "QRAINSNOM F*")
    println("-"^118)

    rows = Tuple{Int, Dict{String,Float64}, Dict{String,Float64}}[]
    for ns in steps
        (forcing_file, step_date) = parity_forcing(DATE0 + Hour(ns - N0))
        inst, _ = run_one_parity_step!(ns; use_cn = true, dumpdir = DUMPDIR_SUMMER,
            step_date = step_date, forcing_file = forcing_file, use_hydrstress = true,
            fndep = FNDEP, h2osfcflag = h2oflag)
        J = julia_chain(inst)
        i = _rec(ns)
        Fs = Dict(k => v[i] for (k, v) in F)
        # Fortran qflx_rain_plus_snomelt, recovered exactly from QSATEXCS = FSAT*q
        Fs["QRAINSNOM"] = Fs["FSAT"] > 0 ? Fs["QSATEXCS"] / Fs["FSAT"] : NaN
        push!(rows, (ns, J, Fs))
        rat(a, b) = (abs(b) > 0 ? a / b : (abs(a) > 0 ? Inf : 1.0))
        @printf("%-9d | %11.4e %11.4e %9.3f | %11.4e %11.4e %9.5f | %11.4e %11.4e\n",
                ns, J["QOVER"], Fs["QOVER"], rat(J["QOVER"], Fs["QOVER"]),
                J["FSAT"], Fs["FSAT"], rat(J["FSAT"], Fs["FSAT"]),
                J["QRAINSNOM"], Fs["QRAINSNOM"])
    end

    # -------- per-field summary over the whole chain --------
    println("\n" * "="^118)
    println("WHOLE-CHAIN summary (max |J-F|, and max relative error scaled by the field's own range)")
    println("="^118)
    fields = ["FSAT", "FCOV", "RAIN", "QSNOMELT", "FSNO", "QRAINSNOM", "QSATEXCS",
              "QTOPSOIL", "QINSOIL", "QINFLEXC",
              "QINFLEXCS", "QH2OSFC", "FH2OSFC", "H2OSFC", "QINSOILL", "QINFL",
              "QTS2SFC", "QINH2OSFC", "QSFCDRAIN", "QOVER", "QDRAI", "QRGWL", "ZWT"]
    @printf("%-11s %13s %13s %13s %13s\n", "field", "max|abs|", "max|rel|", "mean J", "mean F")
    println("-"^118)
    for fl in fields
        (haskey(rows[1][2], fl) && haskey(rows[1][3], fl)) || continue
        scale = maximum(abs(r[3][fl]) for r in rows if isfinite(r[3][fl]))
        scale <= 0 && (scale = 1.0)
        mabs = 0.0
        for r in rows
            (isfinite(r[2][fl]) && isfinite(r[3][fl])) || continue
            mabs = max(mabs, abs(r[2][fl] - r[3][fl]))
        end
        mj = sum(r[2][fl] for r in rows) / length(rows)
        mf = sum(r[3][fl] for r in rows) / length(rows)
        flag = mabs / scale > 1e-6 ? "  <<< DIFF" : ""
        @printf("%-11s %13.4e %13.4e %13.4e %13.4e%s\n", fl, mabs, mabs / scale, mj, mf, flag)
    end

    # -------- the decisive identity --------
    println("\n" * "="^118)
    println("IDENTITY CHECK  QOVER == QSATEXCS + QINFLEXCS + QH2OSFC  (both codes), and")
    println("                QSATEXCS == FSAT * qflx_rain_plus_snomelt")
    println("="^118)
    for (ns, J, Fs) in rows[1:min(5, length(rows))]
        @printf("n=%d  F: QOVER=%.4e  = %.4e + %.4e + %.4e\n", ns, Fs["QOVER"],
                Fs["QSATEXCS"], Fs["QINFLEXCS"], Fs["QH2OSFC"])
        @printf("          J: QOVER=%.4e  = %.4e + %.4e + %.4e\n", J["QOVER"],
                J["QSATEXCS"], J["QINFLEXCS"], J["QH2OSFC"])
    end
end

main()
