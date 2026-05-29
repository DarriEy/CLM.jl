#!/usr/bin/env julia
# ==========================================================================
# Export CLM test data for Haskell port verification
#
# Runs CLM.jl initialization + 1-year simulation at Bow-at-Banff and
# dumps all data as raw Float64 binary files with a JSON manifest.
#
# Output format:
#   - Raw little-endian Float64 arrays (column-major for matrices)
#   - JSON manifest with array names, shapes, and byte offsets
#
# Usage:
#   julia --project=. scripts/export_test_data.jl [output_dir]
# ==========================================================================

using Dates, JSON, Printf
using CLM

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

const DATA_DIR = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped"
const FSURDAT   = joinpath(DATA_DIR, "settings/CLM/parameters/surfdata_clm.nc")
const PARAMFILE = joinpath(DATA_DIR, "settings/CLM/parameters/clm5_params.nc")
const FORCING_DIR = joinpath(DATA_DIR, "data/forcing/CLM_input")
const FSNOWOPTICS = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_optics_5bnd_c013122.nc"
const FSNOWAGING  = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_drdt_bst_fit_60_c070416.nc"

const YEAR = 2002
const START_DATE = DateTime(YEAR, 1, 1)
const END_DATE   = DateTime(YEAR + 1, 1, 1)
const DTIME = 1800  # seconds
const BASEFLOW_SCALAR = 0.0022
const INT_SNOW_MAX = 3113.2

# Output directory
const OUTDIR = length(ARGS) >= 1 ? ARGS[1] :
    "/Users/darri.eythorsson/Library/CloudStorage/GoogleDrive-dareyt@gmail.com/My Drive/code/CLM-hs/test/data"

mkpath(OUTDIR)

# ---------------------------------------------------------------------------
# Helper: write raw Float64 array to binary file
# ---------------------------------------------------------------------------

function write_binary(filepath::String, data::AbstractArray{Float64})
    open(filepath, "w") do io
        write(io, vec(data))  # column-major order
    end
end

function write_binary(filepath::String, data::AbstractVector{<:Real})
    open(filepath, "w") do io
        write(io, Float64.(data))
    end
end

function write_binary(filepath::String, data::AbstractMatrix{<:Real})
    open(filepath, "w") do io
        write(io, Float64.(vec(data)))  # column-major
    end
end

function write_scalar(filepath::String, val::Real)
    open(filepath, "w") do io
        write(io, Float64(val))
    end
end

function write_int_binary(filepath::String, data::AbstractVector{<:Integer})
    open(filepath, "w") do io
        write(io, Int64.(data))
    end
end

# ---------------------------------------------------------------------------
# Phase 1: Initialize CLM
# ---------------------------------------------------------------------------

println("=== Phase 1: Initializing CLM ===")

(inst, bounds, filt, tm) = CLM.clm_initialize!(;
    fsurdat = FSURDAT,
    paramfile = PARAMFILE,
    start_date = START_DATE,
    dtime = DTIME,
    use_cn = false,
    use_aquifer_layer = false,
    use_hydrstress = false,
    use_luna = false,
    use_cndv = false,
    fsnowoptics = FSNOWOPTICS,
    fsnowaging = FSNOWAGING,
    int_snow_max = INT_SNOW_MAX
)

ng = bounds.endg - bounds.begg + 1
nc = bounds.endc - bounds.begc + 1
np = bounds.endp - bounds.begp + 1
nl = bounds.endl - bounds.begl + 1
nlevsoi = CLM.varpar.nlevsoi
nlevgrnd = CLM.varpar.nlevgrnd
nlevsno = CLM.varpar.nlevsno
nlevtot = nlevsno + nlevgrnd

println("  Grid: ng=$ng, nl=$nl, nc=$nc, np=$np")
println("  Levels: nlevsoi=$nlevsoi, nlevgrnd=$nlevgrnd, nlevsno=$nlevsno")

# ---------------------------------------------------------------------------
# Phase 2: Export surface data
# ---------------------------------------------------------------------------

println("=== Phase 2: Exporting surface data ===")

surf = inst.surfdata

surfdir = joinpath(OUTDIR, "surfdata")
mkpath(surfdir)

# Soil properties (per layer)
write_binary(joinpath(surfdir, "pct_sand.bin"), Float64.(surf.pct_sand))
write_binary(joinpath(surfdir, "pct_clay.bin"), Float64.(surf.pct_clay))
write_binary(joinpath(surfdir, "organic.bin"), Float64.(surf.organic))

# Soil color
write_int_binary(joinpath(surfdir, "soil_color.bin"), Int64.([surf.soil_color[1]]))

# Topography
write_binary(joinpath(surfdir, "fmax.bin"), Float64.([surf.fmax[1]]))
write_binary(joinpath(surfdir, "slope.bin"), Float64.([surf.slope[1]]))
write_binary(joinpath(surfdir, "std_elev.bin"), Float64.([surf.std_elev[1]]))

# Lake depth
lakedepth = hasproperty(surf, :lakedepth) ? Float64.([surf.lakedepth[1]]) : Float64.([0.0])
write_binary(joinpath(surfdir, "lakedepth.bin"), lakedepth)

# Zbedrock
zbedrock = hasproperty(surf, :zbedrock) ? Float64.([surf.zbedrock[1]]) : Float64.([8.5])
write_binary(joinpath(surfdir, "zbedrock.bin"), zbedrock)

# Landunit weights
wt_lunit = Float64.(surf.wt_lunit)
write_binary(joinpath(surfdir, "wt_lunit.bin"), wt_lunit)

# PFT weights
wt_nat_patch = Float64.(surf.wt_nat_patch)
write_binary(joinpath(surfdir, "wt_nat_patch.bin"), wt_nat_patch)

# Monthly phenology: LAI, SAI, HTOP, HBOT [ng, npft, 12]
write_binary(joinpath(surfdir, "monthly_lai.bin"), Float64.(surf.monthly_lai))
write_binary(joinpath(surfdir, "monthly_sai.bin"), Float64.(surf.monthly_sai))
write_binary(joinpath(surfdir, "monthly_htop.bin"), Float64.(surf.monthly_htop))
write_binary(joinpath(surfdir, "monthly_hbot.bin"), Float64.(surf.monthly_hbot))

# Latitude/longitude
write_scalar(joinpath(surfdir, "latitude.bin"), inst.gridcell.lat[1])
write_scalar(joinpath(surfdir, "longitude.bin"), inst.gridcell.lon[1])

println("  Surface data exported to $surfdir")

# ---------------------------------------------------------------------------
# Phase 3: Export parameters (PFT constants + scalar params)
# ---------------------------------------------------------------------------

println("=== Phase 3: Exporting parameters ===")

paramdir = joinpath(OUTDIR, "params")
mkpath(paramdir)

# PFT constants - all vectors
pftcon = CLM.pftcon
numpft = length(pftcon.z0mr)

pft_fields = [
    :z0mr, :displar, :dleaf, :c3psn, :vcmx25, :mp, :flnr, :slatop,
    :dsladlai, :leafcn, :fnitr, :woody, :lflitcn, :frootcn, :livewdcn,
    :deadwdcn, :froot_leaf, :stem_leaf, :croot_stem, :flivewd,
    :fcur, :lf_flab, :lf_fcel, :lf_flig, :fr_flab, :fr_fcel, :fr_flig,
    :leaf_long, :evergreen, :stress_decid, :season_decid,
    :roota_par, :rootb_par, :smpso, :smpsc,
    :xl, :rhol, :rhos, :taul, :taus,
    :qe25, :theta_cj, :bbbopt, :mbbopt, :nstor, :br_mr,
]

for f in pft_fields
    if hasproperty(pftcon, f)
        val = getproperty(pftcon, f)
        if val isa AbstractVector
            write_binary(joinpath(paramdir, "pftcon_$(f).bin"), Float64.(val))
        elseif val isa AbstractMatrix
            write_binary(joinpath(paramdir, "pftcon_$(f).bin"), Float64.(val))
        end
    end
end

# Photosynthesis scalar parameters (from CLM.params_inst)
pi = CLM.params_inst
for (name, val) in [
    ("theta_ip", pi.theta_ip), ("act25", pi.act25), ("fnr", pi.fnr),
    ("cp25_yr2000", pi.cp25_yr2000), ("kc25_coef", pi.kc25_coef),
    ("ko25_coef", pi.ko25_coef), ("fnps", pi.fnps), ("theta_psii", pi.theta_psii),
    ("vcmaxha", pi.vcmaxha), ("jmaxha", pi.jmaxha), ("tpuha", pi.tpuha),
    ("lmrha", pi.lmrha), ("kcha", pi.kcha), ("koha", pi.koha), ("cpha", pi.cpha),
    ("vcmaxhd", pi.vcmaxhd), ("jmaxhd", pi.jmaxhd), ("tpuhd", pi.tpuhd),
    ("lmrhd", pi.lmrhd), ("lmrse", pi.lmrse), ("tpu25ratio", pi.tpu25ratio),
    ("kp25ratio", pi.kp25ratio), ("vcmaxse_sf", pi.vcmaxse_sf),
    ("jmaxse_sf", pi.jmaxse_sf), ("tpuse_sf", pi.tpuse_sf),
    ("jmax25top_sf", pi.jmax25top_sf),
]
    write_scalar(joinpath(paramdir, "phot_$(name).bin"), val)
end

# Canopy fluxes params
cfp = CLM.canopy_fluxes_params
write_scalar(joinpath(paramdir, "canflux_csoilc.bin"), cfp.csoilc)
write_scalar(joinpath(paramdir, "canflux_cv.bin"), cfp.cv)

# Snow layer constants (stored as Ref{Vector})
write_binary(joinpath(paramdir, "snow_dzmin.bin"), Float64.(CLM.SNOW_DZMIN[]))
write_binary(joinpath(paramdir, "snow_dzmax_u.bin"), Float64.(CLM.SNOW_DZMAX_U[]))
write_binary(joinpath(paramdir, "snow_dzmax_l.bin"), Float64.(CLM.SNOW_DZMAX_L[]))

println("  Parameters exported to $paramdir")

# ---------------------------------------------------------------------------
# Phase 4: Export cold-start state
# ---------------------------------------------------------------------------

println("=== Phase 4: Exporting cold-start state ===")

statedir = joinpath(OUTDIR, "coldstart")
mkpath(statedir)

# Grid hierarchy
col = inst.column
pch = inst.patch
lun = inst.landunit
grc = inst.gridcell

# Column data
write_binary(joinpath(statedir, "col_z.bin"), Float64.(col.z))
write_binary(joinpath(statedir, "col_dz.bin"), Float64.(col.dz))
write_binary(joinpath(statedir, "col_zi.bin"), Float64.(col.zi))
write_int_binary(joinpath(statedir, "col_gridcell.bin"), Int64.(col.gridcell))
write_int_binary(joinpath(statedir, "col_landunit.bin"), Int64.(col.landunit))
write_int_binary(joinpath(statedir, "col_itype.bin"), Int64.(col.itype))
write_binary(joinpath(statedir, "col_wtgcell.bin"), Float64.(col.wtgcell))
write_binary(joinpath(statedir, "col_wtlunit.bin"), Float64.(col.wtlunit))

# Landunit data
write_int_binary(joinpath(statedir, "lun_itype.bin"), Int64.(lun.itype))
write_int_binary(joinpath(statedir, "lun_gridcell.bin"), Int64.(lun.gridcell))
write_binary(joinpath(statedir, "lun_wtgcell.bin"), Float64.(lun.wtgcell))
if hasproperty(lun, :urbpoi)
    write_binary(joinpath(statedir, "lun_urbpoi.bin"), Float64.(lun.urbpoi))
end
if hasproperty(lun, :lakpoi)
    write_binary(joinpath(statedir, "lun_lakpoi.bin"), Float64.(lun.lakpoi))
end
write_int_binary(joinpath(statedir, "lun_ncolumns.bin"), Int64.(lun.ncolumns))

# Patch data
write_int_binary(joinpath(statedir, "pch_column.bin"), Int64.(pch.column))
write_int_binary(joinpath(statedir, "pch_landunit.bin"), Int64.(pch.landunit))
write_int_binary(joinpath(statedir, "pch_gridcell.bin"), Int64.(pch.gridcell))
write_int_binary(joinpath(statedir, "pch_itype.bin"), Int64.(pch.itype))
write_binary(joinpath(statedir, "pch_wtgcell.bin"), Float64.(pch.wtgcell))
write_binary(joinpath(statedir, "pch_wtcol.bin"), Float64.(pch.wtcol))
write_binary(joinpath(statedir, "pch_wtlunit.bin"), Float64.(pch.wtlunit))
write_binary(joinpath(statedir, "pch_active.bin"), Float64.(pch.active))

# Temperature state
temp = inst.temperature
write_binary(joinpath(statedir, "t_soisno.bin"), Float64.(temp.t_soisno_col))
write_binary(joinpath(statedir, "t_grnd.bin"), Float64.(temp.t_grnd_col))
write_binary(joinpath(statedir, "t_h2osfc.bin"), Float64.(temp.t_h2osfc_col))
write_binary(joinpath(statedir, "t_ref2m.bin"), Float64.(temp.t_ref2m_patch))
write_binary(joinpath(statedir, "t_veg.bin"), Float64.(temp.t_veg_patch))

# Water state
wsb = inst.water.waterstatebulk_inst
ws = wsb.ws
write_binary(joinpath(statedir, "h2osoi_liq.bin"), Float64.(ws.h2osoi_liq_col))
write_binary(joinpath(statedir, "h2osoi_ice.bin"), Float64.(ws.h2osoi_ice_col))
write_binary(joinpath(statedir, "h2osoi_vol.bin"), Float64.(ws.h2osoi_vol_col))
write_binary(joinpath(statedir, "h2osno.bin"), Float64.(ws.h2osno_no_layers_col))
write_binary(joinpath(statedir, "h2osfc.bin"), Float64.(ws.h2osfc_col))
write_binary(joinpath(statedir, "h2ocan.bin"), Float64.(ws.liqcan_patch))

# Soil properties
ss = inst.soilstate
write_binary(joinpath(statedir, "watsat.bin"), Float64.(ss.watsat_col))
write_binary(joinpath(statedir, "bsw.bin"), Float64.(ss.bsw_col))
write_binary(joinpath(statedir, "sucsat.bin"), Float64.(ss.sucsat_col))
write_binary(joinpath(statedir, "hksat.bin"), Float64.(ss.hksat_col))
write_binary(joinpath(statedir, "tkmg.bin"), Float64.(ss.tkmg_col))
write_binary(joinpath(statedir, "tksatu.bin"), Float64.(ss.tksatu_col))
write_binary(joinpath(statedir, "tkdry.bin"), Float64.(ss.tkdry_col))
write_binary(joinpath(statedir, "csol.bin"), Float64.(ss.csol_col))
write_binary(joinpath(statedir, "watdry.bin"), Float64.(ss.watdry_col))
write_binary(joinpath(statedir, "watopt.bin"), Float64.(ss.watopt_col))
write_binary(joinpath(statedir, "watfc.bin"), Float64.(ss.watfc_col))

# Soil hydrology
sh = inst.soilhydrology
write_binary(joinpath(statedir, "zwt.bin"), Float64.(sh.zwt_col))

# Canopy state
cs = inst.canopystate
write_binary(joinpath(statedir, "elai.bin"), Float64.(cs.elai_patch))
write_binary(joinpath(statedir, "esai.bin"), Float64.(cs.esai_patch))
write_binary(joinpath(statedir, "tlai.bin"), Float64.(cs.tlai_patch))
write_binary(joinpath(statedir, "tsai.bin"), Float64.(cs.tsai_patch))
write_binary(joinpath(statedir, "htop.bin"), Float64.(cs.htop_patch))
write_binary(joinpath(statedir, "hbot.bin"), Float64.(cs.hbot_patch))

# Filter masks
write_binary(joinpath(statedir, "filt_soilc.bin"), Float64.(filt.soilc))
write_binary(joinpath(statedir, "filt_lakec.bin"), Float64.(filt.lakec))
write_binary(joinpath(statedir, "filt_nolakec.bin"), Float64.(filt.nolakec))
write_binary(joinpath(statedir, "filt_soilp.bin"), Float64.(filt.soilp))

# Dimensions
dims = Dict(
    "ng" => ng, "nl" => nl, "nc" => nc, "np" => np,
    "nlevsoi" => nlevsoi, "nlevgrnd" => nlevgrnd,
    "nlevsno" => nlevsno, "nlevtot" => nlevtot,
    "numpft" => numpft,
    "dtime" => DTIME, "year" => YEAR,
    "lat" => inst.gridcell.lat[1],
    "lon" => inst.gridcell.lon[1],
)

println("  Cold-start state exported to $statedir")

# ---------------------------------------------------------------------------
# Phase 5: Export forcing data for full year
# ---------------------------------------------------------------------------

println("=== Phase 5: Exporting forcing data ===")

forcdir = joinpath(OUTDIR, "forcing")
mkpath(forcdir)

# Read all forcing timesteps using Julia's forcing reader
fr = CLM.ForcingReader()
fforcing = joinpath(FORCING_DIR, "clmforc.$YEAR.nc")
CLM.forcing_reader_init!(fr, fforcing)

ntimes = fr.ntimes
println("  Forcing file: $fforcing ($ntimes timesteps)")

# Allocate arrays for all timesteps
forc_tbot   = zeros(ntimes)
forc_psrf   = zeros(ntimes)
forc_wind   = zeros(ntimes)
forc_flds   = zeros(ntimes)
forc_fsds   = zeros(ntimes)
forc_precip = zeros(ntimes)
forc_qbot   = zeros(ntimes)

# Use the already-initialized a2l from instances
a2l_tmp = inst.atm2lnd

for ti in 1:ntimes
    target_time = fr.times[ti]
    CLM.read_forcing_step!(fr, a2l_tmp, target_time, ng, nc)
    forc_tbot[ti]   = a2l_tmp.forc_t_not_downscaled_grc[1]
    forc_psrf[ti]   = a2l_tmp.forc_pbot_not_downscaled_grc[1]
    forc_wind[ti]   = a2l_tmp.forc_wind_grc[1]
    forc_flds[ti]   = a2l_tmp.forc_lwrad_not_downscaled_grc[1]
    forc_fsds[ti]   = a2l_tmp.forc_solad_not_downscaled_grc[1, 1] +
                      a2l_tmp.forc_solad_not_downscaled_grc[1, 2] +
                      a2l_tmp.forc_solai_grc[1, 1] +
                      a2l_tmp.forc_solai_grc[1, 2]
    forc_precip[ti] = a2l_tmp.forc_rain_not_downscaled_grc[1] +
                      a2l_tmp.forc_snow_not_downscaled_grc[1]
    forc_qbot[ti]   = a2l_tmp.forc_q_not_downscaled_grc[1]
end

CLM.forcing_reader_close!(fr)

write_binary(joinpath(forcdir, "tbot.bin"), forc_tbot)
write_binary(joinpath(forcdir, "psrf.bin"), forc_psrf)
write_binary(joinpath(forcdir, "wind.bin"), forc_wind)
write_binary(joinpath(forcdir, "flds.bin"), forc_flds)
write_binary(joinpath(forcdir, "fsds.bin"), forc_fsds)
write_binary(joinpath(forcdir, "precip.bin"), forc_precip)
write_binary(joinpath(forcdir, "qbot.bin"), forc_qbot)

dims["ntimes"] = ntimes

println("  Forcing exported: $ntimes timesteps × 7 fields")

# ---------------------------------------------------------------------------
# Phase 6: Run full year simulation, export daily-average reference output
# ---------------------------------------------------------------------------

println("=== Phase 6: Running full-year simulation ===")

inst_final = CLM.clm_run!(;
    fsurdat = FSURDAT,
    paramfile = PARAMFILE,
    fforcing = fforcing,
    fhistory = joinpath(OUTDIR, "julia_history.nc"),
    start_date = START_DATE,
    end_date = END_DATE,
    dtime = DTIME,
    use_cn = false,
    use_aquifer_layer = false,
    use_hydrstress = false,
    use_luna = false,
    use_cndv = false,
    irrigate = false,
    use_voc = false,
    baseflow_scalar = BASEFLOW_SCALAR,
    int_snow_max = INT_SNOW_MAX,
    fsnowoptics = FSNOWOPTICS,
    fsnowaging = FSNOWAGING,
    verbose = false,
)

# Extract daily averages from the history file we just wrote
using NCDatasets
ds = NCDataset(joinpath(OUTDIR, "julia_history.nc"), "r")

refdir = joinpath(OUTDIR, "reference")
mkpath(refdir)

# Write CSV for easy comparison
ref_vars = ["T_GRND", "TSA", "FSA", "EFLX_LH_TOT", "EFLX_SH_TOT",
            "H2OSNO", "QRUNOFF", "ZWT", "TLAI", "QFLX_EVAP_TOT",
            "SNOW_DEPTH", "FRAC_SNO"]

open(joinpath(OUTDIR, "julia_daily_avg.csv"), "w") do io
    # Header
    println(io, "day," * join(ref_vars, ","))

    ndays = length(ds["time"])
    for d in 1:ndays
        vals = String[]
        for v in ref_vars
            if haskey(ds, v)
                val = ds[v][1, d]
                push!(vals, @sprintf("%.8e", ismissing(val) ? -9999.0 : Float64(val)))
            else
                push!(vals, @sprintf("%.8e", -9999.0))
            end
        end
        println(io, "$d," * join(vals, ","))
    end
end

close(ds)

println("  Reference daily averages written to julia_daily_avg.csv")

# ---------------------------------------------------------------------------
# Phase 7: Write JSON manifest
# ---------------------------------------------------------------------------

println("=== Phase 7: Writing manifest ===")

manifest = Dict(
    "description" => "CLM.jl test data export for Haskell port verification",
    "site" => "Bow at Banff (single column)",
    "year" => YEAR,
    "created" => string(now()),
    "dimensions" => dims,
    "surfdata" => Dict(
        "pct_sand" => Dict("shape" => size(surf.pct_sand), "file" => "surfdata/pct_sand.bin"),
        "pct_clay" => Dict("shape" => size(surf.pct_clay), "file" => "surfdata/pct_clay.bin"),
        "organic" => Dict("shape" => size(surf.organic), "file" => "surfdata/organic.bin"),
        "monthly_lai" => Dict("shape" => size(surf.monthly_lai), "file" => "surfdata/monthly_lai.bin"),
        "monthly_sai" => Dict("shape" => size(surf.monthly_sai), "file" => "surfdata/monthly_sai.bin"),
        "monthly_htop" => Dict("shape" => size(surf.monthly_htop), "file" => "surfdata/monthly_htop.bin"),
        "monthly_hbot" => Dict("shape" => size(surf.monthly_hbot), "file" => "surfdata/monthly_hbot.bin"),
        "wt_lunit" => Dict("shape" => size(surf.wt_lunit), "file" => "surfdata/wt_lunit.bin"),
        "wt_nat_patch" => Dict("shape" => size(surf.wt_nat_patch), "file" => "surfdata/wt_nat_patch.bin"),
    ),
    "forcing" => Dict(
        "ntimes" => ntimes,
        "fields" => ["tbot", "psrf", "wind", "flds", "fsds", "precip", "qbot"],
        "units" => ["K", "Pa", "m/s", "W/m2", "W/m2", "mm/s", "kg/kg"],
    ),
    "coldstart" => Dict(
        "col_z" => Dict("shape" => size(col.z), "file" => "coldstart/col_z.bin"),
        "col_dz" => Dict("shape" => size(col.dz), "file" => "coldstart/col_dz.bin"),
        "col_zi" => Dict("shape" => size(col.zi), "file" => "coldstart/col_zi.bin"),
        "t_soisno" => Dict("shape" => size(temp.t_soisno_col), "file" => "coldstart/t_soisno.bin"),
        "h2osoi_liq" => Dict("shape" => size(ws.h2osoi_liq_col), "file" => "coldstart/h2osoi_liq.bin"),
        "h2osoi_ice" => Dict("shape" => size(ws.h2osoi_ice_col), "file" => "coldstart/h2osoi_ice.bin"),
        "watsat" => Dict("shape" => size(ss.watsat_col), "file" => "coldstart/watsat.bin"),
    ),
    "reference" => Dict(
        "variables" => ref_vars,
        "file" => "julia_daily_avg.csv",
    ),
)

open(joinpath(OUTDIR, "manifest.json"), "w") do io
    JSON.print(io, manifest, 2)
end

println("  Manifest written to manifest.json")
println("\n=== Export complete! ===")
println("  Output directory: $OUTDIR")
println("  Total files: $(length(readdir(OUTDIR, join=true)))")
