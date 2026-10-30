# fates_planthydro_parity.jl — FATES plant hydraulics (use_fates_planthydro) vs Fortran CTSM.
#
# Steps the Julia port on the fates_parity_1pt site (Bow at Banff, 2002-07-01, 3600 s,
# FATES + FATES soil BGC + plant hydraulics with the default 1-D Taylor solver) and
# scores the hydraulics and host-water diagnostics step by step against an hourly h0
# history of the same single-point Fortran case (stage_planthydro_run.sh):
#
#   PSI_LEAF/STEM/TROOT/AROOT  n-weighted mean plant water potential      Pa
#   TH_LEAF/STEM               n-weighted mean leaf/stem water content    m3/m3
#   BTRAN_CO                   n-weighted mean cohort hydraulic btran     -
#   FATES_VEGH2O               plant-stored water                         kg/m2
#   FATES_ROOTUPTAKE           total root water uptake                    kg/m2/s
#   QROOTSINK_j                host soil sink by layer (= FATES uptake)   mm/s
#   QVEGT, BTRAN, EFLX_LH_TOT, FATES_GPP, SOILLIQ_j, H2OSOI_j
#
# The Fortran plant means are the FATES_*_SZPF history fields weighted by
# FATES_NPLANT_SZPF; the Julia ones are computed from the cohorts directly
# (Σ n·x / Σ n), the same quantity. The means are over the WOODY PFTs (1-11): for
# grasses (allom_smode = 2) ctsm5.3.012's FATES `bsap_allom` never assigns its
# intent(out) `sapw_area`, so the reference's grass hydraulics run on uninitialized
# stack memory (binary-layout dependent; see README §"PLANT HYDRAULICS"). The port
# uses upstream FATES's SapwoodAreaGrass. Grass PFTs are reported per PFT
# (PSI_LEAF_pft12..14) for transparency.
#
# usage:
#   julia --project=. scripts/validation/fates_fortran_parity/fates_planthydro_parity.jl \
#         [fortran_h0.nc | fates_planthydro_fortran.csv] [nsteps=240] [out.csv]
#   (no argument => the committed fates_planthydro_fortran.csv)
#   ... fates_planthydro_parity.jl export <fortran_h0.nc> [out.csv]   # refresh the extract
#   FATES_HYDRO_SOLVER=2 ... fates_planthydro_parity.jl fates_planthydro_picard_fortran.csv
#                                                        # 2-D Picard solver (README)
using CLM, Printf, Dates, NCDatasets
include(joinpath(@__DIR__, "..", "..", "fates_fortran_parity.jl"))   # FSURDAT, FPARAM, FORCDIR, START, DTIME

const NLAY = 6        # soil layers scored for the by-layer fields
const PLANT = ("PSI_LEAF", "PSI_STEM", "PSI_TROOT", "PSI_AROOT", "TH_LEAF", "TH_STEM", "BTRAN_CO")
const SZPF = Dict("PSI_LEAF" => "FATES_LEAF_H2OPOT_SZPF", "PSI_STEM" => "FATES_STEM_H2OPOT_SZPF",
                  "PSI_TROOT" => "FATES_TRANSROOT_H2OPOT_SZPF", "PSI_AROOT" => "FATES_ABSROOT_H2OPOT_SZPF",
                  "TH_LEAF" => "FATES_LEAF_H2O_SZPF", "TH_STEM" => "FATES_STEM_H2O_SZPF",
                  "BTRAN_CO" => "FATES_BTRAN_SZPF")
const SCALARS = ("QVEGT", "BTRAN", "EFLX_LH_TOT", "FATES_GPP", "FATES_VEGH2O", "FATES_ROOTUPTAKE",
                 "FATES_SAPFLOW", "FATES_NCOHORTS", "FATES_NPATCHES")
const LAYERED = ("QROOTSINK", "SOILLIQ", "H2OSOI")
const NPFT = 14
const NWOODY = 11   # PFTs 1-11 are woody in the default parameter file; 12-14 are grasses
# PFT-resolved leaf water potential / hydraulic btran (attribution of the means).
const PERPFT = Tuple(f * "_pft" * string(k) for f in ("PSI_LEAF", "BTRAN_CO") for k in 1:NPFT)
# PFT-resolved demography (attribution of post-daily-step divergence): density and
# leaf / storage carbon, = FATES_NPLANT_PF [m-2], FATES_LEAFC_PF / FATES_STOREC_PF [kgC m-2].
const DEMOG = Dict("NPLANT" => "FATES_NPLANT_PF", "LEAFC" => "FATES_LEAFC_PF",
                   "STOREC" => "FATES_STOREC_PF")
const PERPFT_DEMOG = Tuple(f * "_pft" * string(k) for f in ("NPLANT", "LEAFC", "STOREC") for k in 1:NPFT)
const ALL_FIELDS = (PLANT..., SCALARS...,
                    (f * "_" * string(j) for f in LAYERED for j in 1:NLAY)..., PERPFT...,
                    PERPFT_DEMOG...)

# FATES_HYDRO_SOLVER=1|2|3 selects the plant-hydraulics solver (unset: the parameter
# file's fates_hydro_solver, 1). Score a solver against a Fortran run made with the
# same value in its parameter file.
const HYDRO_SOLVER = haskey(ENV, "FATES_HYDRO_SOLVER") ? parse(Int, ENV["FATES_HYDRO_SOLVER"]) : nothing

function build_hydro()
    fates_case_snow_namelist!()   # the Fortran case lnd_in snow namelist (fates_fortran_parity.jl)
    inst, bounds, filt, tm = CLM.clm_initialize!(; fsurdat=FSURDAT, paramfile=FPARAM,
        fsnowoptics=FSNOWOPTICS, fsnowaging=FSNOWAGING,
        use_fates=true, use_fates_planthydro=true, fates_hydro_solver=HYDRO_SOLVER,
        start_date=START, dtime=Int(DTIME), use_bedrock=false)
    for g in 1:bounds.endg
        isfinite(inst.gridcell.latdeg[g]) || (inst.gridcell.latdeg[g] = inst.gridcell.lat[g]*180/π)
        isfinite(inst.gridcell.londeg[g]) || (inst.gridcell.londeg[g] = inst.gridcell.lon[g]*180/π)
    end
    for c in 1:bounds.endc
        inst.landunit.itype[inst.column.landunit[c]] == CLM.ISTSOIL && (inst.column.is_fates[c] = true)
    end
    config = CLM.CLMDriverConfig(use_fates=true)
    return inst, bounds, filt, CLM.clump_filter_inactive_and_active, config, tm
end

_hist(inst, ih::Symbol) = (h = inst.fates.hist; h.hvars[h.ih[ih]].r81d[1])

"Per-step Julia record for the (first) FATES column / site."
function hydro_record(inst, c)
    r = Dict{String,Float64}()
    site = inst.fates.sites[1]
    acc = Dict(k => 0.0 for k in PLANT); ntot = 0.0
    pl = zeros(NPFT); pb = zeros(NPFT); pn = zeros(NPFT)
    dn = zeros(NPFT); dl = zeros(NPFT); ds = zeros(NPFT)
    cp = site.oldest_patch
    while cp !== nothing
        cc = cp.tallest
        while cc !== nothing
            h = cc.co_hydr
            if h !== nothing && CLM.prt_params.woody[cc.pft] == 1   # woody mean (see header)
                n = cc.n
                ntot += n
                # Same definitions as update_history_hydraulics (FatesHistoryInterfaceMod
                # .F90:5930-5950): ag(1) = leaf, ag(2) = stem, v-weighted aroot mean.
                acc["PSI_LEAF"]  += n * h.psi_ag[1] * 1e6
                acc["PSI_STEM"]  += n * h.psi_ag[2] * 1e6
                acc["PSI_TROOT"] += n * h.psi_troot * 1e6
                acc["PSI_AROOT"] += n * 1e6 * sum(h.psi_aroot .* h.v_aroot_layer) / sum(h.v_aroot_layer)
                acc["TH_LEAF"]   += n * h.th_ag[1]
                acc["TH_STEM"]   += n * h.th_ag[2]
                acc["BTRAN_CO"]  += n * h.btran
            end
            if h !== nothing
                ft = cc.pft; n = cc.n
                pl[ft] += n * h.psi_ag[1] * 1e6; pb[ft] += n * h.btran; pn[ft] += n
            end
            ftd = cc.pft
            dn[ftd] += cc.n
            dl[ftd] += cc.n * CLM.GetState(cc.prt, CLM.leaf_organ, CLM.carbon12_element)
            ds[ftd] += cc.n * CLM.GetState(cc.prt, CLM.store_organ, CLM.carbon12_element)
            cc = cc.shorter
        end
        cp = cp.younger
    end
    for k in PLANT
        r[k] = ntot > 0 ? acc[k] / ntot : NaN
    end
    for k in 1:NPFT
        r["PSI_LEAF_pft$k"] = pn[k] > 0 ? pl[k] / pn[k] : NaN
        r["BTRAN_CO_pft$k"] = pn[k] > 0 ? pb[k] / pn[k] : NaN
    end
    for k in 1:NPFT
        r["NPLANT_pft$k"] = dn[k] * CLM.AREA_INV
        r["LEAFC_pft$k"]  = dl[k] * CLM.AREA_INV
        r["STOREC_pft$k"] = ds[k] * CLM.AREA_INV
    end
    wf = inst.water.waterfluxbulk_inst
    r["QVEGT"] = wf.wf.qflx_tran_veg_col[c]
    # BTRAN / LH: the history p2g, wtcol-weighted over EVERY active patch of the
    # column. BTRAN's l2g_scale_type='veg' (EnergyFluxType.F90:663) selects the
    # istsoil/istcrop LANDUNITS (subgridAveMod.F90:1077-1079), not FATES's patch%is_veg,
    # and p2c skips only inactive / zero-weight / spval patches (subgridAveMod.F90:203-204).
    # So the FATES bare-ground patch counts: its btran is 0, not spval, because
    # BareGroundFluxes zeroes btran on the no-exposed-veg filter
    # (BareGroundFluxesMod.F90:284), and wrap_btran writes only patchi+1..patchi+npatch
    # (clmfates_interfaceMod.F90:2509-2519). Host BTRAN is therefore the vegetated
    # patches' btran times their canopy fraction (sum of canopy_fraction_pa, ~0.62 here).
    pch = inst.patch; ef = inst.energyflux
    wb = 0.0; bt = 0.0; lh = 0.0; wl = 0.0
    for p in eachindex(pch.column)
        pch.column[p] == c || continue
        pch.active[p] || continue
        w = pch.wtcol[p]
        lh += w * ef.eflx_lh_tot_patch[p]; wl += w
        b = ef.btran_patch[p]
        if isfinite(b) && b != CLM.SPVAL
            bt += w * b; wb += w
        end
    end
    r["BTRAN"] = wb > 0 ? bt / wb : NaN
    r["EFLX_LH_TOT"] = wl > 0 ? lh / wl : NaN
    r["FATES_GPP"] = _hist(inst, :ih_gpp_si)
    r["FATES_VEGH2O"] = _hist(inst, :ih_h2oveg_si)
    r["FATES_ROOTUPTAKE"] = _hist(inst, :ih_rootuptake_si)
    r["FATES_SAPFLOW"] = _hist(inst, :ih_sapflow_si)
    let ncoh = 0, npat = 0, cp = inst.fates.sites[1].oldest_patch
        while cp !== nothing
            npat += 1
            cc = cp.tallest
            while cc !== nothing; ncoh += 1; cc = cc.shorter; end
            cp = cp.younger
        end
        r["FATES_NCOHORTS"] = ncoh; r["FATES_NPATCHES"] = npat
    end
    joff = CLM.varpar.nlevsno
    ws = inst.water.waterstatebulk_inst.ws
    for j in 1:NLAY
        r["QROOTSINK_$j"] = wf.qflx_rootsoi_col[c, j]
        r["SOILLIQ_$j"] = ws.h2osoi_liq_col[c, joff + j]
        r["H2OSOI_$j"] = ws.h2osoi_vol_col[c, j]
    end
    # Julia-only diagnostics for the water-balance check.
    r["PLANT_H2O_COL"] = inst.water.waterdiagnosticbulk_inst.total_plant_stored_h2o_col[c]
    r["ERRH2O_COL"] = inst.water.waterbalancebulk_inst.errh2o_col[c]
    return r
end

function run_julia(nsteps::Int; hook=nothing)
    inst, bounds, filt, filt_ia, config, _tm = build_hydro()
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
        push!(recs, hydro_record(inst, c))
        hook === nothing || hook(inst, i)
    end
    CLM.forcing_reader_close!(fr)
    return recs, inst
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
        f64(x) = Float64.(coalesce.(x, NaN))
        out = Dict{String,Vector{Float64}}()
        for f in SCALARS
            haskey(d, f) && (out[f] = f64(vec(Array(d[f][:]))))
        end
        for f in LAYERED
            haskey(d, f) || continue
            a = Array(d[f][:, :, :])
            for j in 1:NLAY
                out[f * "_" * string(j)] = f64(vec(a[1, j, :]))
            end
        end
        for (k, v) in DEMOG
            haskey(d, v) || continue
            x = f64(Array(d[v][:, :, :]))[1, :, :]                  # (pft, time)
            for ft in 1:NPFT
                out[k * "_pft" * string(ft)] = x[ft, :]
            end
        end
        np = f64(Array(d["FATES_NPLANT_SZPF"][:, :, :]))[1, :, :]   # (scpf, time)
        for (k, v) in SZPF
            x = f64(Array(d[v][:, :, :]))[1, :, :]
            nsc = size(x, 1) ÷ NPFT          # scpf = (pft-1)*nsizeclass + sizeclass
            wd = 1:NWOODY * nsc
            out[k] = [sum(x[wd, t] .* np[wd, t]) / sum(np[wd, t]) for t in axes(x, 2)]
            if k in ("PSI_LEAF", "BTRAN_CO")
                for ft in 1:NPFT
                    rg = (ft - 1) * nsc + 1:ft * nsc
                    out[k * "_pft" * string(ft)] =
                        [sum(x[rg, t] .* np[rg, t]) / sum(np[rg, t]) for t in axes(x, 2)]
                end
            end
        end
        out
    end
end

"""
    score(jl, ft)

Per-step relative error |J-F| / max(|F|, floor) over nstep >= 1, floor = 1e-6 x the
field's peak |F|. Also reports the Pearson correlation of the series.
"""
function score(jl, ft; io=stdout)
    n = min(length(jl), minimum(length(v) for v in values(ft)))
    @printf(io, "  %-17s %13s %13s %9s %9s %7s\n", "field", "Fortran(end)", "Julia(end)",
            "max rel", "mean rel", "corr")
    out = Dict{String,NTuple{3,Float64}}()
    for f in ALL_FIELDS
        haskey(ft, f) || continue
        F = ft[f][1:n]; J = [jl[k][f] for k in 1:n]
        idx = [k for k in 2:n if isfinite(F[k]) && isfinite(J[k])]
        isempty(idx) && continue
        peak = maximum(abs, F[idx])
        floor_ = max(peak * 1e-6, 1e-30)
        rel = [abs(J[k] - F[k]) / max(abs(F[k]), floor_) for k in idx]
        mx = maximum(rel); mn = sum(rel) / length(rel)
        Fi, Ji = F[idx], J[idx]
        sF, sJ = sum(abs2, Fi .- sum(Fi)/length(Fi)), sum(abs2, Ji .- sum(Ji)/length(Ji))
        cr = (sF > 0 && sJ > 0) ? sum((Fi .- sum(Fi)/length(Fi)) .* (Ji .- sum(Ji)/length(Ji))) / sqrt(sF*sJ) : NaN
        out[f] = (mx, mn, cr)
        @printf(io, "  %-17s %13.5e %13.5e %9.2e %9.2e %7.3f\n", f, F[n], J[n], mx, mn, cr)
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

function write_csv(path, jl, ft)
    jcols = sort(collect(keys(jl[1])))
    open(path, "w") do io
        println(io, join(["nstep"; [f * "_jl" for f in jcols]; [f * "_ft" for f in ALL_FIELDS]], ","))
        n = min(length(jl), minimum(length(v) for v in values(ft)))
        for k in 1:n
            println(io, join([string(k - 1); [string(jl[k][f]) for f in jcols];
                              [string(haskey(ft, f) ? ft[f][k] : NaN) for f in ALL_FIELDS]], ","))
        end
    end
end

const FORTRAN_CSV = joinpath(@__DIR__, "fates_planthydro_fortran.csv")

function planthydro_main(args)
    if !isempty(args) && args[1] == "export"
        export_fortran_csv(args[2], length(args) >= 3 ? args[3] : FORTRAN_CSV)
        return 0
    end
    ref = isempty(args) ? FORTRAN_CSV : args[1]
    nsteps = length(args) >= 2 ? parse(Int, args[2]) : 240
    ft = read_fortran(ref)
    jl, _ = run_julia(nsteps)
    println("\n=== FATES plant hydraulics: Julia vs Fortran (", nsteps, " steps) ===")
    score(jl, ft)
    errs = [abs(r["ERRH2O_COL"]) for r in jl[2:end] if isfinite(r["ERRH2O_COL"])]
    @printf("\n  column water balance |errh2o|: max %.3e mm, mean %.3e mm over %d steps\n",
            maximum(errs), sum(errs) / length(errs), length(errs))
    length(args) >= 3 && write_csv(args[3], jl, ft)
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(planthydro_main(ARGS))
end
