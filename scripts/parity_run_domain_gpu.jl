#!/usr/bin/env julia
# ==========================================================================
# parity_run_domain_gpu.jl — GPU (Metal Float32) annual parity run.
#
# Same spun-up-vs-spun-up annual run as parity_run_domain.jl, but the physics
# (clm_drv!) runs on a Metal device mirror while forcing-in / downscale / history
# stay host-side (the two-instance-mirror path added to CLM.clm_run! via the
# `device_adapt` kwarg). Writes paper/data/julia_clm_<domain>_metal_<year>.nc so
# the parity plotter can show a third series (GPU) vs Julia-CPU vs Fortran-CPU.
#
#   DOMAIN=Steppe julia +1.12 --project=scripts scripts/parity_run_domain_gpu.jl
#   GPU_NDAYS=7 DOMAIN=Steppe julia +1.12 --project=scripts …   # short shakeout
#
# Metal has no Float64, so this run is Float32 end-to-end — the third series shows
# the genuine Float32 drift over the year, not a bug.
# ==========================================================================

using Dates, Printf, CLM
import Metal
include(joinpath(@__DIR__, "gpu_adapt.jl"))   # mf(MtlArray, x): host tree -> device Float32
mf(x) = mf(Metal.MtlArray, x)

const run_clm! = getfield(CLM, Symbol("clm_run!"))
const DATA = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data"
_snicar(f) = let d1 = "$DATA/installs/cesm-inputdata/lnd/clm2/snicardata",
                 d2 = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata"
    isfile(joinpath(d1, f)) ? joinpath(d1, f) : joinpath(d2, f)
end

# Per-domain registry (mirror of parity_run_domain.jl; pilot carries Steppe —
# extend with more entries to scale the GPU comparison to all domains).
const DOMAINS = Dict(
    "Steppe" => (year = 2017, dtime = 3600, baseflow = 0.001, int_snow = 2000.0,
        caldir  = "$DATA/domain_Steppe_Kherlen_Mongolia/settings/CLM/parameters",
        forcing = "$DATA/domain_Steppe_Kherlen_Mongolia/data/forcing/CLM_input/clmforc.2017.nc",
        restart = "$DATA/domain_Steppe_Kherlen_Mongolia/simulations/clm_steppe/CLM/Steppe_Kherlen_Mongolia.clm2.r.2017-01-01-00000.nc",
        h0      = "$DATA/domain_Steppe_Kherlen_Mongolia/simulations/clm_steppe/CLM/Steppe_Kherlen_Mongolia.clm2.h0.2016-12-31-00000.nc"),
)

const DOM = get(ENV, "DOMAIN", "Steppe")
haskey(DOMAINS, DOM) || error("Unknown DOMAIN=$DOM; have $(collect(keys(DOMAINS)))")
cfg = DOMAINS[DOM]
yr  = cfg.year

if !Metal.functional()
    println("Metal not functional on this host — aborting GPU parity run."); exit(2)
end

fsurdat   = joinpath(cfg.caldir, "surfdata_clm.nc")
paramfile = joinpath(cfg.caldir, "clm5_params.nc")
outdir    = abspath(joinpath(@__DIR__, "..", "paper", "data"))
fhistory  = joinpath(outdir, "julia_clm_$(lowercase(DOM))_metal_$(yr).nc")
isdir(outdir) || mkpath(outdir)

# Optional short window for a device shakeout (GPU_NDAYS=7 → one week).
ndays = parse(Int, get(ENV, "GPU_NDAYS", "0"))
start_date = DateTime(yr, 1, 1)
end_date   = ndays > 0 ? start_date + Day(ndays) : DateTime(yr + 1, 1, 1)

println("="^64)
println("  CLM.jl GPU (Metal Float32) parity run — $DOM $(yr)  [$(ndays>0 ? "$ndays-day shakeout" : "full year")]")
println("="^64)
println("  output -> $fhistory")

# Same physics config as the CPU reference run.
CLM.rooting_profile_config.rooting_profile_method_water  = CLM.JACKSON_1996_ROOT
CLM.rooting_profile_config.rooting_profile_method_carbon = CLM.JACKSON_1996_ROOT
CLM.snow_hydrology_set_control_for_testing!(;
    wind_dep_snow_density = true,
    overburden_compaction_method = CLM.OVERBURDEN_COMPACTION_VIONNET2012)

t0 = time()
run_clm!(;
    fsurdat = fsurdat, paramfile = paramfile, fforcing = cfg.forcing,
    fhistory = fhistory,
    start_date = start_date, end_date = end_date,
    dtime = cfg.dtime, use_cn = false, verbose = true,
    use_aquifer_layer = false, use_hydrstress = true, use_luna = true,
    h2osfcflag = 1,
    baseflow_scalar = cfg.baseflow, int_snow_max = cfg.int_snow,
    ffortran_restart = cfg.restart,
    interp_forcing = true,
    forcing_phase_shift_s = parse(Int, get(ENV, "FORCING_SHIFT_S", "0")),
    fsnowoptics = isfile(_snicar("snicar_optics_5bnd_c013122.nc")) ? _snicar("snicar_optics_5bnd_c013122.nc") : "",
    fsnowaging  = isfile(_snicar("snicar_drdt_bst_fit_60_c070416.nc")) ? _snicar("snicar_drdt_bst_fit_60_c070416.nc") : "",
    device_adapt = mf)   # <-- the GPU switch: run clm_drv! on a Metal Float32 mirror
@printf("\n  Done (Metal): %s in %.1f s -> %s\n", DOM, time() - t0, fhistory)
