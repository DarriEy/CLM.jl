# ==========================================================================
# fortran_parity_soilparam.jl — FORCING-FREE soil-parameter parity check.
#
# Hypothesis: CLM.jl computes watsat/bsw/sucsat from a pedotransfer function
# (cold_start.jl) instead of reading surfdata, so its soil hydraulic params
# may differ from Fortran's — which would corrupt the soil-water state (and
# hence BTRAN) independent of any forcing or time-stepping.
#
# Test: inject the Fortran `before_step` dump (shared water state), then for
# each soil layer compute Julia's matric potential from Julia's params and the
# SHARED water content, and compare to the Fortran `SMP` in the SAME dump.
# A mismatch isolates a pedotransfer/param difference — no forcing involved.
#
#   julia +1.12 --project=. scripts/fortran_parity_soilparam.jl [nstep]
# ==========================================================================

include(joinpath(@__DIR__, "fortran_parity_common.jl"))

const NSTEP = length(ARGS) >= 1 ? ARGS[1] : "8761"
const DUMP  = joinpath(DUMPDIR, "pdump_before_step_n$(NSTEP).nc")

const DENH2O = 1000.0
const DENICE = 917.0

println("="^72)
println("  FORCING-FREE soil-parameter parity check  (nstep=$NSTEP)")
println("="^72)

(inst, bounds, filt, tm) = build_bow_inst()
inject_dump!(inst, bounds, DUMP)

ds = NCDataset(DUMP, "r")
smp_f = _dumpvar(ds, "SMP")          # Fortran matric potential (mm), levgrnd
hk_f  = _dumpvar(ds, "HK")
watsat_f = _dumpvar(ds, "WATSAT_P")  # Fortran soil params (now dumped)
bsw_f    = _dumpvar(ds, "BSW_P")
sucsat_f = _dumpvar(ds, "SUCSAT_P")
close(ds)

# ---- Direct soil-param comparison (Julia pedotransfer vs Fortran) ----
if watsat_f !== nothing
    ss0 = inst.soilstate
    println("\n  Direct param comparison (Julia pedotransfer vs Fortran):")
    @printf("  %-4s | %9s %9s %8s | %9s %9s %8s | %9s %9s %8s\n",
            "lay","wsat_jl","wsat_ft","Δ","bsw_jl","bsw_ft","Δ","suc_jl","suc_ft","Δ")
    println("  " * "-"^96)
    for j in 1:CLM.varpar.nlevsoi
        @printf("  %-4d | %9.4f %9.4f %8.1e | %9.3f %9.3f %8.1e | %9.2f %9.2f %8.1e\n",
            j, ss0.watsat_col[1,j], watsat_f[j], abs(ss0.watsat_col[1,j]-watsat_f[j]),
            ss0.bsw_col[1,j], bsw_f[j], abs(ss0.bsw_col[1,j]-bsw_f[j]),
            ss0.sucsat_col[1,j], sucsat_f[j], abs(ss0.sucsat_col[1,j]-sucsat_f[j]))
    end
    println("  " * "-"^96)
end

ss   = inst.soilstate
ws   = inst.water.waterstatebulk_inst.ws
col  = inst.column
nlevsno = CLM.varpar.nlevsno
nlevsoi = CLM.varpar.nlevsoi
c = 1

@printf("\n  %-4s %9s %9s %9s | %11s %11s %9s | %9s\n",
        "lay", "watsat", "bsw", "sucsat", "smp_julia", "smp_fortran", "rel", "s_node")
println("  " * "-"^82)
maxrel = 0.0
for j in 1:nlevsoi
    lev = nlevsno + j
    dz  = col.dz[c, lev]
    liq = ws.h2osoi_liq_col[c, lev]
    ice = ws.h2osoi_ice_col[c, lev]
    wsat = ss.watsat_col[c, j]; suc = ss.sucsat_col[c, j]; b = ss.bsw_col[c, j]
    vol_liq = liq / (dz * DENH2O)
    eff_por = max(wsat - ice / (dz * DENICE), 0.01)
    s_node  = clamp(vol_liq / eff_por, 0.01, 1.0)
    smp_j   = CLM.soil_suction_clapp_hornberger(suc, s_node, b)   # mm
    smp_ft  = (smp_f === nothing || j > length(smp_f)) ? NaN : smp_f[j]
    rel = isnan(smp_ft) ? NaN : abs(smp_j - smp_ft) / (1.0 + max(abs(smp_j), abs(smp_ft)))
    isnan(rel) || (global maxrel = max(maxrel, rel))
    @printf("  %-4d %9.4f %9.3f %9.3f | %11.2f %11.2f %9.2e | %9.4f\n",
            j, wsat, b, suc, smp_j, smp_ft, rel, s_node)
end
println("  " * "-"^82)
@printf("  max|rel| matric-potential mismatch (Julia params vs Fortran SMP) = %.3e\n", maxrel)
println()
if maxrel > 1e-3
    println("  => SOIL PARAMS / matric-potential DIFFER. Pedotransfer (watsat/bsw/sucsat)")
    println("     or the suction formula is a root cause feeding the soil-water + BTRAN drift.")
else
    println("  => Soil params + matric potential MATCH at the shared state; the soil-water")
    println("     drift originates downstream (infiltration / Richards / snow), not params.")
end
