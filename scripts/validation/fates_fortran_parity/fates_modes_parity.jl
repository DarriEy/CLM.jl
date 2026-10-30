# fates_modes_parity.jl — reduced-complexity FATES modes vs Fortran CTSM.
#
# Steps the Julia port on the fates_parity_1pt site (Bow at Banff, 2002-07-01, 3600 s)
# in one of three FATES modes and scores it against an hourly h0 history of the same
# single-point Fortran case staged by stage_fates_modes_run.sh:
#
#   sp       use_fates_sp + use_fates_nocomp + use_fates_fixed_biogeog  (no soil BGC)
#   nocomp   use_fates_nocomp + use_fates_fixed_biogeog                 (soil BGC on)
#   fixedbg  use_fates_fixed_biogeog                                    (soil BGC on)
#
# Scored: FATES_GPP, TLAI, ELAI, TSAI, HTOP, FSH, EFLX_LH_TOT, FATES_NCOHORTS,
# FATES_NPATCHES, FATES_VEGC, FATES_LAI, and per PFT FATES_NOCOMP_PATCHAREA_PF (nocomp
# modes), FATES_CANOPYCROWNAREA_PF, FATES_NPLANT_PF, FATES_LEAFC_PF, and the host
# fields of HIST_NAMES as written by the port's history writer (HIST_*).
#
# usage:
#   julia --project=. scripts/validation/fates_fortran_parity/fates_modes_parity.jl \
#         <mode> [fortran_h0.nc | fates_<mode>_fortran.csv] [nsteps=240]
#   (no reference => the committed fates_<mode>_fortran.csv)
#   ... fates_modes_parity.jl export <mode> <fortran_h0.nc>   # refresh the committed CSV
using CLM, Printf, Dates, NCDatasets
include(joinpath(@__DIR__, "..", "..", "fates_fortran_parity.jl"))   # FSURDAT, FPARAM, FORCDIR, START, DTIME

const MODE_FLAGS = Dict(
    "sp"      => (use_fates_sp = true,  use_fates_nocomp = true,  use_fates_fixed_biogeog = true),
    "nocomp"  => (use_fates_sp = false, use_fates_nocomp = true,  use_fates_fixed_biogeog = true),
    "fixedbg" => (use_fates_sp = false, use_fates_nocomp = false, use_fates_fixed_biogeog = true))

const SCALARS = ("FATES_GPP", "TLAI", "ELAI", "TSAI", "HTOP", "FSH", "EFLX_LH_TOT",
                 "FATES_NCOHORTS", "FATES_NPATCHES", "FATES_VEGC", "FATES_LAI")
const PF_FIELDS = ("FATES_NOCOMP_PATCHAREA_PF", "FATES_CANOPYCROWNAREA_PF",
                   "FATES_NPLANT_PF", "FATES_LEAFC_PF")
const NPFT = 14
# Host fields scored as the port's HISTORY WRITER aggregates them to the gridcell
# (CLM._aggregate_to_gridcell, the p2g of hist_dov2xy), not via a harness-side
# average: key "HIST_<h0 name>". BTRAN is the hourly instantaneous btran with the
# Fortran l2g_scale_type='veg' (EnergyFluxType.F90 hist_addfld BTRAN); the port's
# default "BTRAN" history field is the daily-min snapshot (BTRANMN).
const HIST_NAMES = ("BTRAN", "TV", "FSH", "EFLX_LH_TOT", "FSA", "FSR", "FIRA", "FCTR",
                    "FCEV", "FGEV", "FGR", "SABV", "SABG", "TSA", "Q2M", "TLAI", "FIRE")
const HIST_FIELDS = Tuple("HIST_" * n for n in HIST_NAMES)
const ALL_FIELDS = (SCALARS..., (f * "_" * string(k) for f in PF_FIELDS for k in 1:NPFT)...,
                    HIST_FIELDS...)

"History-writer field definitions for HIST_NAMES."
function hist_defs()
    d = Dict(f.name => f for f in CLM.default_hist_fields())
    d["BTRAN"] = CLM.HistFieldDef("BTRAN", "transpiration beta factor", "unitless", "patch",
        inst -> Float64.(inst.energyflux.btran_patch), :avg, "veg")
    return [d[n] for n in HIST_NAMES]
end
const HIST_DEFS = hist_defs()

function build_mode(mode)
    fl = MODE_FLAGS[mode]
    fates_case_snow_namelist!()   # the Fortran case lnd_in snow namelist (fates_fortran_parity.jl)
    inst, bounds, filt, tm = CLM.clm_initialize!(; fsurdat=FSURDAT, paramfile=FPARAM,
        fsnowoptics=FSNOWOPTICS, fsnowaging=FSNOWAGING,
        use_fates=true, fl..., fates_history_dimlevel=2, start_date=START, dtime=Int(DTIME), use_bedrock=false)
    for g in 1:bounds.endg
        isfinite(inst.gridcell.latdeg[g]) || (inst.gridcell.latdeg[g] = inst.gridcell.lat[g]*180/π)
        isfinite(inst.gridcell.londeg[g]) || (inst.gridcell.londeg[g] = inst.gridcell.lon[g]*180/π)
    end
    col, lun = inst.column, inst.landunit
    for c in 1:bounds.endc
        lun.itype[col.landunit[c]] == CLM.ISTSOIL && (col.is_fates[c] = true)
    end
    config = CLM.CLMDriverConfig(use_fates=true, use_fates_sp=fl.use_fates_sp)
    return inst, bounds, filt, CLM.clump_filter_inactive_and_active, config, tm
end

_hv(h, ih) = h.hvars[h.ih[ih]]
# 0 => not registered in this mode (e.g. the nocomp-only FATES_NOCOMP_PATCHAREA_PF)
_has(h, ih) = haskey(h.ih, ih) && h.ih[ih] > 0

"Column (= gridcell, single column) weighted average of a patch field over active patches."
function _p2c(inst, c, v)
    pch = inst.patch; col = inst.column
    s = 0.0; w = 0.0
    for p in col.patchi[c]:col.patchf[c]
        (pch.active[p] && pch.wtcol[p] > 0) || continue
        s += pch.wtcol[p] * v[p]; w += pch.wtcol[p]
    end
    return w > 0 ? s / w : NaN
end

function mode_record(inst, c)
    r = Dict{String,Float64}()
    h = inst.fates.hist
    cs = inst.canopystate; ef = inst.energyflux
    r["FATES_GPP"] = _hv(h, :ih_gpp_si).r81d[1]
    # Under SP, CTSM's TLAI/TSAI/HTOP history fields point at the *_hist copies (the
    # FATES totals; CanopyStateType.F90), tlai/tsai/htop being the prescribed inputs.
    sp = CLM.hlm_use_sp[] == CLM.itrue
    r["TLAI"] = _p2c(inst, c, sp ? cs.tlai_hist_patch : cs.tlai_patch)
    r["ELAI"] = _p2c(inst, c, cs.elai_patch)
    r["TSAI"] = _p2c(inst, c, sp ? cs.tsai_hist_patch : cs.tsai_patch)
    r["HTOP"] = _p2c(inst, c, sp ? cs.htop_hist_patch : cs.htop_patch)
    r["FSH"] = _p2c(inst, c, ef.eflx_sh_tot_patch)
    r["EFLX_LH_TOT"] = _p2c(inst, c, ef.eflx_lh_tot_patch)
    r["FATES_NCOHORTS"] = _hv(h, :ih_ncohorts_si).r81d[1]
    r["FATES_NPATCHES"] = _hv(h, :ih_npatches_si).r81d[1]
    r["FATES_VEGC"] = _hv(h, :ih_totvegc_si).r81d[1]
    r["FATES_LAI"] = _hv(h, :ih_lai_si).r81d[1]
    for (f, ih) in (("FATES_NOCOMP_PATCHAREA_PF", :ih_nocomp_pftpatchfraction_si_pft),
                    ("FATES_CANOPYCROWNAREA_PF", :ih_canopycrownarea_si_pft),
                    ("FATES_NPLANT_PF", :ih_nindivs_si_pft),
                    ("FATES_LEAFC_PF", :ih_leafbiomass_si_pft))
        a = _has(h, ih) ? _hv(h, ih).r82d : nothing
        for k in 1:NPFT
            r[f * "_" * string(k)] = a === nothing ? NaN : a[1, k]
        end
    end
    ng = maximum(inst.column.gridcell)
    for (key, def) in zip(HIST_FIELDS, HIST_DEFS)
        v = CLM._aggregate_to_gridcell(def, inst, ng)
        r[key] = (v === nothing || v[1] == -9999.0) ? NaN : v[1]
    end
    return r
end

"Step the port nsteps (plus the NUOPC nstep=0 init step), recording after every step."
function run_julia(mode, nsteps::Int; hook = nothing)
    inst, bounds, filt, filt_ia, config, _tm = build_mode(mode)
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
    for i in 0:nsteps   # identical clock convention to fates_fortran_parity.jl
        cur = START + Second(round(Int, max(i - 1, 0) * DTIME))
        tod = hour(cur)*3600 + minute(cur)*60 + second(cur)
        calday = dayofyear(cur) + tod / CLM.SECSPDAY
        (declin, _) = CLM.compute_orbital(calday)
        (declinm1, _) = CLM.compute_orbital(calday - DTIME / CLM.SECSPDAY)
        CLM.init_daylength!(inst.gridcell, declin, declinm1, CLM.ORB_OBLIQR_DEFAULT, 1:ng)
        nextsw = calday + DTIME / CLM.SECSPDAY
        CLM.read_forcing_step!(fr, inst.atm2lnd, cur, ng, nc)
        CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
        CLM.clm_drv!(config, inst, filt, filt_ia, bounds, CLM.nuopc_doalb(i), nextsw, declin, declin,
            CLM.ORB_OBLIQR_DEFAULT, false, false, "", false;
            nstep=i, is_first_step=(i <= 1),
            is_beg_curr_day=(i > 0 && tod == 0), is_end_curr_day=false, is_beg_curr_year=false,
            dtime=DTIME, mon=month(cur), day=day(cur), secs=tod,
            jday=dayofyear(cur), photosyns=inst.photosyns)
        hook === nothing || hook(inst, i)
        push!(recs, mode_record(inst, c))
    end
    CLM.forcing_reader_close!(fr)
    return recs
end

"Fortran series (record k <-> nstep k-1) from the hourly h0 (`.nc`) or the committed CSV."
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
        for (key, n) in zip(HIST_FIELDS, HIST_NAMES)
            haskey(d, n) && (out[key] = Float64.(vec(Array(d[n][:]))))
        end
        for f in PF_FIELDS
            haskey(d, f) || continue
            a = Array(d[f][:, :, :])            # (lndgrid, fates_levpft, time)
            for k in 1:NPFT
                out[f * "_" * string(k)] = Float64.(vec(a[1, k, :]))
            end
        end
        out
    end
end

"""
Per-field comparison over nstep >= 1: relative error |J-F| / max(|F|, floor) with a
floor of 1e-3 x the field's peak |F| (so night-time / absent-PFT zeros do not dominate),
plus the window means and the end values. PFT columns that are zero in both models
are skipped.
"""
function score(jl, ft; io=stdout)
    n = min(length(jl), minimum(length(v) for v in values(ft)))
    @printf(io, "  %-30s %12s %12s %12s %12s %9s %9s\n", "field", "F mean", "J mean",
            "F end", "J end", "mean rel", "rel(mean)")
    out = Dict{String,NTuple{3,Float64}}()
    for f in ALL_FIELDS
        haskey(ft, f) || continue
        F = ft[f][2:n]; J = [jl[k][f] for k in 2:n]
        all(isnan, J) && continue
        peak = max(maximum(abs, F), maximum(x -> isnan(x) ? 0.0 : abs(x), J))
        peak == 0 && continue
        floor_ = peak * 1e-3
        rel = [abs(J[k] - F[k]) / max(abs(F[k]), floor_) for k in eachindex(F)]
        Fm = sum(F) / length(F); Jm = sum(J) / length(J)
        relm = abs(Jm - Fm) / max(abs(Fm), floor_)
        mn = sum(rel) / length(rel)
        out[f] = (mn, relm, maximum(rel))
        @printf(io, "  %-30s %12.5e %12.5e %12.5e %12.5e %9.2e %9.2e\n", f, Fm, Jm,
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

fortran_csv(mode) = joinpath(@__DIR__, "fates_$(mode)_fortran.csv")

function modes_main(args)
    if !isempty(args) && args[1] == "export"
        export_fortran_csv(args[3], length(args) >= 4 ? args[4] : fortran_csv(args[2]))
        return 0
    end
    mode = args[1]
    haskey(MODE_FLAGS, mode) || error("mode must be one of $(keys(MODE_FLAGS))")
    ref = length(args) >= 2 ? args[2] : fortran_csv(mode)
    nsteps = length(args) >= 3 ? parse(Int, args[3]) : 240
    ft = read_fortran(ref)
    jl = run_julia(mode, nsteps)
    println("\n=== FATES mode `", mode, "`: Julia vs Fortran (", nsteps, " steps) ===")
    score(jl, ft)
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(modes_main(ARGS))
end
