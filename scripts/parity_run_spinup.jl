#!/usr/bin/env julia
# ==========================================================================
# parity_run_spinup.jl — regenerate the Bow-at-Banff parity figure data.
#
# Runs the FAIR spun-up-vs-spun-up annual comparison (the same setup as
# test/verify_vs_fortran.jl --spinup): initialize Julia from the Fortran
# spun-up restart (clm2.r.2003-01-01) and run the full year 2003 on real
# observed forcing with the DDS-calibrated parameters, writing daily-average
# history. The companion Fortran reference is the CTSM h0.2003 history at the
# same site. Output (Julia daily NetCDF) lands at a FIXED path so the plotting
# script (paper/plot_parity_spinup.py) can consume it.
#
#   julia +1.12 --project=. scripts/parity_run_spinup.jl
#
# This is the honest comparison: the older README/paper figures cold-started
# Julia (soil too wet to equilibrate in one year with deep-drainage off → a
# +24–30% latent-heat/canopy-flux init artifact, NOT a translation error). The
# spun-up run removes that artifact (residual LH ≈ +1 W/m², T_GRND ≈ −0.2 K).
# ==========================================================================

using NCDatasets, Dates, Printf, CLM
const run_clm! = getfield(CLM, Symbol("clm_run!"))

const YEAR = 2003
const basedir = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped"
const caldir  = joinpath(basedir, "optimization/CLM/dds_run_1/final_evaluation/settings/CLM/parameters")

fsurdat  = joinpath(caldir, "surfdata_clm.nc")
paramfile = joinpath(caldir, "clm5_params.nc")
fforcing = joinpath(basedir, "data/forcing/CLM_input/clmforc.$(YEAR).nc")

# PHS+LUNA toggle (env CLM_PHS=1). The Fortran reference h0.2003 was generated
# with use_hydrstress=.true. + use_luna=.true. (lnd_in), so the apples-to-apples
# comparison is the PHS run; the non-PHS run is kept for the error-cancellation
# contrast. Distinct output paths so both survive.
const USE_PHS  = get(ENV, "CLM_PHS",  "0") == "1"
const USE_LUNA = get(ENV, "CLM_LUNA", USE_PHS ? "1" : "0") == "1"

# Fixed output path for the Julia daily history (consumed by the plot script).
outdir   = abspath(joinpath(@__DIR__, "..", "paper", "data"))
fhistory = joinpath(outdir,
    USE_PHS ? (USE_LUNA ? "julia_clm_bow_spinup_phs_$(YEAR).nc" :
                          "julia_clm_bow_spinup_phsnoluna_$(YEAR).nc") :
              "julia_clm_bow_spinup_$(YEAR).nc")

# Fortran spun-up restart IC + the Fortran reference history (h0).
ffortran_restart = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_parity_run/Bow_at_Banff_lumped.clm2.r.$(YEAR)-01-01-00000.nc"
f_h0 = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_parity_run/Bow_at_Banff_lumped.clm2.h0.$(YEAR)-01-01-00000.nc"

# SNICAR data files. The Fortran lnd_in references /Users/.../projects/cesm-inputdata
# (not present here); the same files live under SYMFLUENCE_data/installs/cesm-inputdata.
# Without them Julia silently falls back to the BATS snow albedo while Fortran ran full
# SNICAR — a scheme mismatch that melted Julia's spring snowpack ~a month early. Prefer
# the existing copy, fall back to the projects path.
const _snicar_dir1 = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/installs/cesm-inputdata/lnd/clm2/snicardata"
const _snicar_dir2 = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata"
_pick(f) = isfile(joinpath(_snicar_dir1, f)) ? joinpath(_snicar_dir1, f) : joinpath(_snicar_dir2, f)
fsnowoptics = _pick("snicar_optics_5bnd_c013122.nc")
fsnowaging  = _pick("snicar_drdt_bst_fit_60_c070416.nc")

# DDS-calibrated parameters (from the Bow user_nl_clm).
baseflow_scalar = 0.0022119554
int_snow_max    = 3113.2227

println("="^64)
println("  CLM.jl parity run — Bow at Banff, $(YEAR) (spun-up-vs-spun-up)")
println("  Mode: ", USE_PHS ? "PHS" : "non-PHS", " + ", USE_LUNA ? "LUNA" : "no-LUNA")
println("  Init: Fortran spun-up restart  | Forcing: clmforc.$(YEAR).nc (obs)")
println("  Julia history → $fhistory")
println("  Fortran ref   → $f_h0")
println("="^64)
for (lbl, f) in (("surfdata", fsurdat), ("params", paramfile), ("forcing", fforcing),
                 ("restart", ffortran_restart), ("fortran h0", f_h0))
    println("  ", isfile(f) ? "OK  " : "MISS", "  $lbl: $f")
end
isdir(outdir) || mkpath(outdir)
println()

t0 = time()
inst = run_clm!(;
    fsurdat=fsurdat, paramfile=paramfile, fforcing=fforcing,
    fhistory=fhistory,
    start_date=DateTime(YEAR, 1, 1),
    end_date=DateTime(YEAR + 1, 1, 1),
    dtime=1800, use_cn=false, verbose=true,
    use_aquifer_layer=false,
    use_hydrstress=USE_PHS, use_luna=USE_LUNA,
    baseflow_scalar=baseflow_scalar,
    int_snow_max=int_snow_max,
    ffortran_restart=ffortran_restart,
    fsnowoptics=isfile(fsnowoptics) ? fsnowoptics : "",
    fsnowaging=isfile(fsnowaging)   ? fsnowaging  : "")
elapsed = time() - t0
@printf("\n  Done: %d timesteps in %.1f s (%.2f s/day)\n", 365*48, elapsed, elapsed/365)

# Quick sanity readout: annual-mean Julia vs Fortran for the headline fluxes.
jds = NCDataset(fhistory); fds = NCDataset(f_h0)
jday(v) = [Float64(jds[v][1, d]) for d in 1:min(365, size(jds[v], 2)) if isfinite(Float64(jds[v][1, d])) && Float64(jds[v][1, d]) > -9990]
fday(v) = haskey(fds, v) ? [Float64(fds[v][1, d]) for d in 1:365 if !ismissing(fds[v][1, d])] : Float64[]
using Statistics
println("\n  Annual-mean check (Julia vs Fortran):")
for (jn, fn, u) in (("EFLX_LH_TOT","EFLX_LH_TOT","W/m2"), ("EFLX_SH_TOT","FSH","W/m2"),
                    ("TSA","TSA","K"), ("FSA","FSA","W/m2"), ("H2OSNO","H2OSNO","mm"),
                    ("FCTR","FCTR","W/m2"), ("FCEV","FCEV","W/m2"), ("FGEV","FGEV","W/m2"))
    j = jday(jn); f = fday(fn)
    (isempty(j) || isempty(f)) && continue
    @printf("    %-12s Julia %8.3f   Fortran %8.3f   Δ %+7.3f  (%s)\n",
            jn, mean(j), mean(f), mean(j)-mean(f), u)
end
close(jds); close(fds)
println("\n  ✓ wrote $fhistory")
