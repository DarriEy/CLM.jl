# fates_landuse_parity.jl — FATES land use (logging, LUH2 transitions / harvest,
# lupft, potential vegetation) vs Fortran CTSM.
#
# Steps the Julia port on the fates_parity_1pt site (Bow at Banff, 2002-07-01, 3600 s)
# in one land-use MODE and scores it against the hourly h0 of the same single-point
# Fortran case staged by stage_fates_landuse_run.sh (same modes, same edited FATES
# parameter file: logging event 2002-07-04, direct-logging dbhmin 0.1 cm; same
# single-point LUH2 files from make_luh_1pt.jl).
#
# The land-use diagnostics the port does not register as history are computed here
# from the FATES site state with the Fortran history formulas
# (FatesHistoryInterfaceMod update_history_dyn1):
#   FATES_NPATCHES_SECONDARY, FATES_SECONDARY_FOREST_FRACTION, FATES_PATCHAREA_LU_<k>,
#   FATES_DISTURBANCE_RATE_LOGGING (sum disturbance_rates(dtype_ilog,:,:) x 365),
#   FATES_HARVEST_WOODPROD_C_FLUX (sum wood_product_harvest / AREA),
#   FATES_WOOD_PRODUCT (trunk_product_site / AREA).
#
# usage:
#   julia --project=. scripts/validation/fates_fortran_parity/fates_landuse_parity.jl \
#         <mode> [fortran_h0.nc | fates_lu_<mode>_fortran.csv] [nsteps=240]
#   ... fates_landuse_parity.jl export <mode> <fortran_h0.nc>   # refresh the committed CSV
using CLM, Printf, Dates, NCDatasets
include(joinpath(@__DIR__, "..", "..", "fates_fortran_parity.jl"))   # FSURDAT, FPARAM, FORCDIR, START, DTIME

const LU_DATA = joinpath(@__DIR__, "..", "..", "..", "data", "fates", "luh_1pt")
const LUH_FILE = get(ENV, "LUH_FILE", joinpath(LU_DATA, "LUH2_1pt_Bow_hist_simyr1850-2015.nc"))
const LUPFT_FILE = get(ENV, "LUPFT_FILE", joinpath(LU_DATA, "fates_landuse_pft_map_1pt_Bow.nc"))
const LOG_EVENT_CODE = get(ENV, "LOG_EVENT_CODE", "20020704")
const LOG_DBHMIN = get(ENV, "LOG_DBHMIN", "0.1")

mode_flags(mode) = Dict(
    "logging"  => (; fates_harvest_mode = "event_code"),
    "luh"      => (; use_fates_luh = true, fluh_timeseries = LUH_FILE),
    "luh_area" => (; use_fates_luh = true, fluh_timeseries = LUH_FILE,
                   fates_harvest_mode = "luhdata_area"),
    "luh_mass" => (; use_fates_luh = true, fluh_timeseries = LUH_FILE,
                   fates_harvest_mode = "luhdata_mass"),
    "lupft"    => (; use_fates_lupft = true, fluh_timeseries = LUH_FILE,
                   flandusepftdat = LUPFT_FILE),
    "potveg"   => (; use_fates_potentialveg = true))[mode]

const SCALARS = ("FATES_GPP", "TLAI", "FSH", "FATES_NCOHORTS", "FATES_NPATCHES",
                 "FATES_VEGC", "FATES_LAI", "FATES_LITTER_IN",
                 "FATES_DISTURBANCE_RATE_TREEFALL", "FATES_NPATCHES_SECONDARY",
                 "FATES_SECONDARY_FOREST_FRACTION", "FATES_DISTURBANCE_RATE_LOGGING",
                 "FATES_HARVEST_WOODPROD_C_FLUX", "FATES_WOOD_PRODUCT",
                 "FATES_LUCHANGE_WOODPROD_C_FLUX")
# Sums over the land-use x land-use history matrices: the whole transition matrix and
# the off-diagonal (land-use change) part of the disturbance-rate matrix.
const MATRIX_SUMS = ("FATES_TRANSITION_MATRIX_LULU_SUM", "FATES_DISTURBANCE_RATE_LUC_SUM")
const NLU = 5
const ALL_FIELDS = (SCALARS..., MATRIX_SUMS...,
                    ("FATES_PATCHAREA_LU_" * string(k) for k in 1:NLU)...)

"The default FATES .cdl with the logging event / dbhmin edits of stage_fates_landuse_run.sh."
function landuse_cdl()
    out = joinpath(mktempdir(), "fates_params_landuse.cdl")
    open(out, "w") do io
        for l in eachline(CLM.FATES_PARAMS_DEFAULT_CDL)
            startswith(l, " fates_landuse_logging_event_code = ") &&
                (l = " fates_landuse_logging_event_code = $LOG_EVENT_CODE ;")
            startswith(l, " fates_landuse_logging_dbhmin = ") &&
                (l = " fates_landuse_logging_dbhmin = $LOG_DBHMIN ;")
            println(io, l)
        end
    end
    return out
end

function build(mode)
    fl = mode_flags(mode)
    inst, bounds, filt, tm = CLM.clm_initialize!(; fsurdat=FSURDAT, paramfile=FPARAM,
        use_fates=true, fl..., fates_paramfile=landuse_cdl(), fates_history_dimlevel=2,
        start_date=START, dtime=Int(DTIME), use_bedrock=false)
    for g in 1:bounds.endg
        isfinite(inst.gridcell.latdeg[g]) || (inst.gridcell.latdeg[g] = inst.gridcell.lat[g]*180/π)
        isfinite(inst.gridcell.londeg[g]) || (inst.gridcell.londeg[g] = inst.gridcell.lon[g]*180/π)
    end
    col, lun = inst.column, inst.landunit
    for c in 1:bounds.endc
        lun.itype[col.landunit[c]] == CLM.ISTSOIL && (col.is_fates[c] = true)
    end
    config = CLM.CLMDriverConfig(use_fates=true)
    return inst, bounds, filt, CLM.clump_filter_inactive_and_active, config, tm
end

_hv(h, ih) = h.hvars[h.ih[ih]]

function _p2c(inst, c, v)
    pch = inst.patch; col = inst.column
    s = 0.0; w = 0.0
    for p in col.patchi[c]:col.patchf[c]
        (pch.active[p] && pch.wtcol[p] > 0) || continue
        s += pch.wtcol[p] * v[p]; w += pch.wtcol[p]
    end
    return w > 0 ? s / w : NaN
end

function lu_record(inst, c)
    r = Dict{String,Float64}()
    h = inst.fates.hist
    r["FATES_GPP"] = _hv(h, :ih_gpp_si).r81d[1]
    r["TLAI"] = _p2c(inst, c, inst.canopystate.tlai_patch)
    r["FSH"] = _p2c(inst, c, inst.energyflux.eflx_sh_tot_patch)
    r["FATES_NCOHORTS"] = _hv(h, :ih_ncohorts_si).r81d[1]
    r["FATES_NPATCHES"] = _hv(h, :ih_npatches_si).r81d[1]
    r["FATES_VEGC"] = _hv(h, :ih_totvegc_si).r81d[1]
    r["FATES_LAI"] = _hv(h, :ih_lai_si).r81d[1]
    r["FATES_LITTER_IN"] = _hv(h, :ih_litter_in_si).r81d[1]
    r["FATES_DISTURBANCE_RATE_TREEFALL"] = _hv(h, :ih_fall_disturbance_rate_si).r81d[1]
    site = inst.fates.sites[1]
    ainv = 1.0 / CLM.area
    nsec = 0.0; fsec = 0.0; alu = zeros(NLU)
    p = site.oldest_patch
    while p !== nothing
        if p.land_use_label == CLM.secondaryland
            nsec += 1; fsec += p.area * ainv
        end
        1 <= p.land_use_label <= NLU && (alu[p.land_use_label] += p.area * ainv)
        p = p.younger
    end
    r["FATES_NPATCHES_SECONDARY"] = nsec
    r["FATES_SECONDARY_FOREST_FRACTION"] = fsec
    for k in 1:NLU; r["FATES_PATCHAREA_LU_" * string(k)] = alu[k]; end
    r["FATES_DISTURBANCE_RATE_LOGGING"] =
        sum(site.disturbance_rates[CLM.dtype_ilog, :, :]) * CLM.days_per_year
    mb = site.mass_balance[CLM.element_pos[CLM.carbon12_element]]
    r["FATES_HARVEST_WOODPROD_C_FLUX"] = ainv * sum(mb.wood_product_harvest)
    r["FATES_WOOD_PRODUCT"] = site.resources_management.trunk_product_site * ainv
    r["FATES_LUCHANGE_WOODPROD_C_FLUX"] = ainv * sum(mb.wood_product_landusechange)
    r["FATES_TRANSITION_MATRIX_LULU_SUM"] = sum(site.landuse_transition_matrix)
    luc = 0.0
    for i in 1:NLU, j in 1:NLU
        i == j && continue
        luc += sum(site.disturbance_rates[:, i, j]) * CLM.days_per_year
    end
    r["FATES_DISTURBANCE_RATE_LUC_SUM"] = luc
    return r
end

function run_julia(mode, nsteps::Int)
    inst, bounds, filt, filt_ia, config, _tm = build(mode)
    c = findfirst(inst.column.is_fates)
    ng, nc = bounds.endg, bounds.endc
    fr = CLM.ForcingReader()
    CLM.forcing_reader_init!(fr, joinpath(FORCDIR, "clmforc.$(year(START)).nc"))
    let tf = joinpath(FORCDIR, "topo_forcing.nc")
        if isfile(tf)
            NCDataset(tf, "r") do ds
                if haskey(ds, "TOPO")
                    ft = Float64(ds["TOPO"][1])
                    for g in 1:ng; inst.atm2lnd.forc_topo_grc[g] = ft; end
                    for cc in 1:nc; inst.topo.topo_col[cc] = ft; end
                end
            end
        end
    end
    recs = Vector{Dict{String,Float64}}()
    for i in 0:nsteps
        cur = START + Second(round(Int, max(i - 1, 0) * DTIME))
        tod = hour(cur)*3600 + minute(cur)*60 + second(cur)
        calday = dayofyear(cur) + tod / CLM.SECSPDAY
        (declin, _) = CLM.compute_orbital(calday)
        (declinm1, _) = CLM.compute_orbital(calday - DTIME / CLM.SECSPDAY)
        CLM.init_daylength!(inst.gridcell, declin, declinm1, CLM.ORB_OBLIQR_DEFAULT, 1:ng)
        nextsw = calday + DTIME / CLM.SECSPDAY
        CLM.read_forcing_step!(fr, inst.atm2lnd, cur, ng, nc)
        CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
        CLM.clm_drv!(config, inst, filt, filt_ia, bounds, true, nextsw, declin, declin,
            CLM.ORB_OBLIQR_DEFAULT, false, false, "", false;
            nstep=i, is_first_step=(i <= 1),
            is_beg_curr_day=(i > 0 && tod == 0), is_end_curr_day=false, is_beg_curr_year=false,
            dtime=DTIME, mon=month(cur), day=day(cur), secs=tod,
            jday=dayofyear(cur), photosyns=inst.photosyns)
        push!(recs, lu_record(inst, c))
    end
    CLM.forcing_reader_close!(fr)
    return recs
end

function read_fortran(path)
    if endswith(path, ".csv")
        lines = readlines(path)
        hdr = split(lines[1], ',')
        cols = [Float64[] for _ in hdr]
        for l in lines[2:end]
            isempty(l) && continue
            for (i, v) in enumerate(split(l, ','))
                push!(cols[i], parse(Float64, v))
            end
        end
        return Dict(String(h) => cols[i] for (i, h) in enumerate(hdr) if h != "nstep")
    end
    NCDataset(path, "r") do d
        out = Dict{String,Vector{Float64}}()
        for f in SCALARS
            haskey(d, f) && (out[f] = Float64.(vec(Array(d[f][:]))))
        end
        if haskey(d, "FATES_TRANSITION_MATRIX_LULU")
            a = Array(d["FATES_TRANSITION_MATRIX_LULU"][:, :, :])   # (lndgrid, lulu, time)
            out["FATES_TRANSITION_MATRIX_LULU_SUM"] = Float64.(vec(sum(a[1, :, :]; dims = 1)))
        end
        if haskey(d, "FATES_DISTURBANCE_RATE_MATRIX_LULU")
            a = Array(d["FATES_DISTURBANCE_RATE_MATRIX_LULU"][:, :, :])
            off = [k for k in 1:NLU^2 if (k - 1) % NLU != (k - 1) ÷ NLU]   # i + 5(j-1), i != j
            out["FATES_DISTURBANCE_RATE_LUC_SUM"] = Float64.(vec(sum(a[1, off, :]; dims = 1)))
        end
        if haskey(d, "FATES_PATCHAREA_LU")
            a = Array(d["FATES_PATCHAREA_LU"][:, :, :])     # (lndgrid, fates_levlanduse, time)
            for k in 1:NLU
                out["FATES_PATCHAREA_LU_" * string(k)] = Float64.(vec(a[1, k, :]))
            end
        end
        out
    end
end

"""
Per-field comparison over nstep >= 1: window means, end values, the mean relative
error |J-F| / max(|F|, 1e-3 peak|F|) and the relative error of the window mean.
Fields zero in both models are skipped.
"""
function score(jl, ft; io=stdout)
    n = min(length(jl), minimum(length(v) for v in values(ft)))
    @printf(io, "  %-34s %12s %12s %12s %12s %9s %9s\n", "field", "F mean", "J mean",
            "F end", "J end", "mean rel", "rel(mean)")
    out = Dict{String,NTuple{3,Float64}}()
    for f in ALL_FIELDS
        haskey(ft, f) || continue
        F = ft[f][2:n]; J = [jl[k][f] for k in 2:n]
        peak = max(maximum(abs, F), maximum(x -> isnan(x) ? 0.0 : abs(x), J))
        peak == 0 && continue
        floor_ = peak * 1e-3
        rel = [abs(J[k] - F[k]) / max(abs(F[k]), floor_) for k in eachindex(F)]
        Fm = sum(F) / length(F); Jm = sum(J) / length(J)
        relm = abs(Jm - Fm) / max(abs(Fm), floor_)
        mn = sum(rel) / length(rel)
        out[f] = (mn, relm, maximum(rel))
        @printf(io, "  %-34s %12.5e %12.5e %12.5e %12.5e %9.2e %9.2e\n", f, Fm, Jm,
                F[end], J[end], mn, relm)
    end
    return out
end

function export_fortran_csv(nc_path, csv_path)
    ft = read_fortran(nc_path)
    names = [f for f in ALL_FIELDS if haskey(ft, f)]
    n = length(ft[names[1]])
    open(csv_path, "w") do io
        println(io, join(["nstep"; names], ","))
        for k in 1:n
            println(io, join([string(k - 1); [repr(ft[f][k]) for f in names]], ","))
        end
    end
    return csv_path
end

fortran_csv(mode) = joinpath(@__DIR__, "fates_lu_$(mode)_fortran.csv")

function landuse_main(args)
    if !isempty(args) && args[1] == "export"
        export_fortran_csv(args[3], length(args) >= 4 ? args[4] : fortran_csv(args[2]))
        return 0
    end
    mode = args[1]
    ref = length(args) >= 2 ? args[2] : fortran_csv(mode)
    nsteps = length(args) >= 3 ? parse(Int, args[3]) : 240
    ft = read_fortran(ref)
    jl = run_julia(mode, nsteps)
    println("\n=== FATES land use `", mode, "`: Julia vs Fortran (", nsteps, " steps) ===")
    score(jl, ft)
    # the daily steps around the logging event (record k <-> nstep k-1)
    println("\n  daily-step trace (nstep: F / J)")
    for f in ("FATES_VEGC", "FATES_LITTER_IN", "FATES_SECONDARY_FOREST_FRACTION",
              "FATES_NCOHORTS", "FATES_NPATCHES", "FATES_GPP", "FATES_WOOD_PRODUCT",
              "FATES_PATCHAREA_LU_3")
        haskey(ft, f) || continue
        @printf("  %-32s", f)
        for n in (1, 24, 48, 72, 73, 96, 120, 168, 240)
            n + 1 <= min(length(jl), length(ft[f])) || continue
            @printf(" %d: %.5g/%.5g", n, ft[f][n + 1], jl[n + 1][f])
        end
        println()
    end
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(landuse_main(ARGS))
end
