#!/usr/bin/env julia
# ==========================================================================
# fates_niche_modes.jl — run the Julia port of the `fates_parity_1pt` case (Bow at
# Banff, 2002-07-01, dtime 3600 s, all 14 PFTs, FATES soil BGC) with one FATES
# research mode switched on, writing one CSV row per model day of the diagnostics the
# matching Fortran daily history carries (scripts/validation/fates_fortran_parity/
# stage_niche_mode.py + fates_niche_compare.jl). See scripts/validation/fates_fortran_parity/README.md "RESEARCH MODES".
#
#   FATES_NICHE_DAYS=14 julia --project=. scripts/fates_niche_modes.jl <mode> <out.csv>
#
# <mode> ∈ base | tree_damage | st3 | prescribed_phys | coage | seeddisp | inventory
#   seeddisp  needs FATES_NICHE_PARAMS=<cdl> (a parameter file with the
#             fates_seed_dispersal_* parameters set; the default file has them unset)
#   inventory needs FATES_NICHE_INVENTORY=<control file>
# Optional FATES_NICHE_RESTART_DAY=<d>: at the end of model day d, pack the FATES
# state into the restart vectors, write them to <out>.fates_rst.nc, rebuild every
# site from that file (create_patchcohort_structure! + get_restart_vectors!) and
# continue: the restart round-trip (compare against the continuous run's CSV).
# Debug: FATES_NICHE_TRACE=d1,d2 prints every cohort after those days;
# FATES_NICHE_RSTDEBUG=1 lists the NaN-holding site/patch/cohort fields at the daily
# step entry of day FATES_NICHE_CENSUS_DAY (default: the day after the restart).
# ==========================================================================
Base.include(@__MODULE__, joinpath(@__DIR__, "fates_fortran_parity.jl"))

const NM_DAYS = parse(Int, get(ENV, "FATES_NICHE_DAYS", "14"))
const NM_RST  = parse(Int, get(ENV, "FATES_NICHE_RESTART_DAY", "0"))

function niche_kwargs(mode::AbstractString)
    pf = get(ENV, "FATES_NICHE_PARAMS", "")
    pkw = isempty(pf) ? (;) : (; fates_paramfile = pf)
    mode == "base"            && return (pkw, (;))
    mode == "tree_damage"     && return ((; pkw..., use_fates_tree_damage = true), (;))
    mode == "st3"             && return ((; pkw..., use_fates_ed_st3 = true), (;))
    mode == "prescribed_phys" && return ((; pkw..., use_fates_ed_prescribed_phys = true), (;))
    mode == "coage"           && return ((; pkw..., use_fates_cohort_age_tracking = true), (;))
    if mode == "seeddisp"
        cad = parse(Int, get(ENV, "FATES_NICHE_CADENCE", "1"))
        return ((; pkw..., fates_seeddisp_cadence = cad), (; fates_seeddisp_cadence = cad))
    end
    if mode == "inventory"
        ctl = ENV["FATES_NICHE_INVENTORY"]
        return ((; pkw..., use_fates_inventory_init = true,
                 fates_inventory_ctrl_filename = ctl), (;))
    end
    error("unknown mode $mode")
end

const NM_ORGANS = (:leaf, :fnrt, :sapw, :store, :struct, :repro)

"Per-day site diagnostics, in the units of the Fortran FATES history (per m2)."
function niche_row(site)
    np = _C.numpft[]
    A = _C.area                      # site area [m2]
    ndmg = _C.nlevdamage[]
    edp = _C.ed_params()
    caedges = edp.ED_val_history_coageclass_bin_edges
    nca = length(caedges)
    r = Dict{String,Float64}()
    ncoh = 0; npatch = 0
    pools = Dict(o => 0.0 for o in NM_ORGANS)
    nplant = zeros(np); nplant_cd = zeros(ndmg, np); nplant_ca = zeros(nca, np)
    seedbank = 0.0
    organ = Dict(:leaf => _C.leaf_organ, :fnrt => _C.fnrt_organ, :sapw => _C.sapw_organ,
                 :store => _C.store_organ, :struct => _C.struct_organ, :repro => _C.repro_organ)
    cp = site.oldest_patch
    while cp !== nothing
        npatch += 1
        lit = cp.litter[1]
        seedbank += (sum(lit.seed) + sum(lit.seed_germ)) * cp.area / A
        cc = cp.tallest
        while cc !== nothing
            ncoh += 1
            ft = cc.pft
            for o in NM_ORGANS
                pools[o] += _pool(cc, organ[o]) * cc.n / A
            end
            nplant[ft] += cc.n / A
            1 <= cc.crowndamage <= ndmg && (nplant_cd[cc.crowndamage, ft] += cc.n / A)
            if isfinite(cc.coage)
                ica = something(findlast(e -> cc.coage >= e, caedges), 1)
                nplant_ca[ica, ft] += cc.n / A
            end
            cc = cc.shorter
        end
        cp = cp.younger
    end
    r["ncoh"] = ncoh; r["npatch"] = npatch
    for o in NM_ORGANS; r["$(o)c"] = pools[o]; end
    r["vegc"] = sum(values(pools))
    r["seed_bank"] = seedbank
    for f in 1:np
        r["nplant_$f"] = nplant[f]
        r["seed_out_$f"] = isempty(site.seed_out) ? 0.0 : site.seed_out[f]
        for d in 1:ndmg; r["nplant_cd$(d)_$f"] = nplant_cd[d, f]; end
        for a in 1:nca;  r["nplant_ca$(a)_$f"] = nplant_ca[a, f]; end
    end
    return r
end

# --------------------------------------------------------------------------
# FATES restart round-trip: pack -> NetCDF file -> rebuild -> unpack.
# --------------------------------------------------------------------------
function fates_restart_roundtrip!(inst, path)
    _C.write_fates_restart(path, inst.fates)
    _C.read_fates_restart!(path, inst)
    return nothing
end

"debug: list fields holding NaN on the site / patches / cohorts at the daily-step entry"
function nan_census(site)
    out = String[]
    chk(tag, obj) = for fn in fieldnames(typeof(obj))
        isdefined(obj, fn) || continue
        v = getfield(obj, fn)
        bad = (v isa Float64 && isnan(v)) || (v isa AbstractArray{Float64} && any(isnan, v))
        bad && push!(out, "$tag.$fn")
    end
    chk("site", site)
    p = site.oldest_patch; ip = 0
    while p !== nothing
        ip += 1; chk("patch$ip", p)
        for (li, L) in enumerate(p.litter); chk("patch$ip.litter$li", L); end
        c = p.tallest
        while c !== nothing; chk("patch$ip.cohort", c); c = c.shorter; end
        p = p.younger
    end
    return sort!(unique(out))
end

function main_niche(mode::AbstractString, out::AbstractString)
    init_kw, cfg_kw = niche_kwargs(mode)
    inst, bounds, filt, filt_ia, config, _tm = build(; init_kw, cfg_kw)
    ng, nc = bounds.endg, bounds.endc
    let tf = joinpath(FORCDIR, "topo_forcing.nc")
        if isfile(tf)
            ds = NCDataset(tf, "r")
            if haskey(ds, "TOPO")
                ft = Float64(ds["TOPO"][1])
                for g in 1:ng; inst.atm2lnd.forc_topo_grc[g] = ft; end
                for c in 1:nc; inst.topo.topo_col[c] = ft; end
            end
            close(ds)
        end
    end
    io = open(out, "w")
    keys_ = nothing
    # daily means of host fields (the Fortran h0 FSH / TV / HR are daily averages):
    # gridcell-weighted over active patches; HR on the first FATES column
    hacc = Dict("fsh" => 0.0, "tv" => 0.0, "hr" => 0.0); hn = Ref(0)
    fc = findfirst(inst.column.is_fates)
    accum_host!() = begin
        pch = inst.patch; ws = 0.0; fsh = 0.0; tv = 0.0
        for p in eachindex(pch.wtgcell)
            (pch.active[p] && isfinite(pch.wtgcell[p]) && pch.wtgcell[p] > 0) || continue
            ws += pch.wtgcell[p]
            fsh += pch.wtgcell[p] * inst.energyflux.eflx_sh_tot_patch[p]
            tv  += pch.wtgcell[p] * inst.temperature.t_veg_patch[p]
        end
        hacc["fsh"] += fsh / ws; hacc["tv"] += tv / ws
        hacc["hr"]  += inst.soilbiogeochem_carbonflux.hr_col[fc]
        hn[] += 1
    end
    write_row(d) = begin
        r = niche_row(inst.fates.sites[1])
        for (k, v) in hacc
            r["host_$k"] = hn[] > 0 ? v / hn[] : NaN
            hacc[k] = 0.0
        end
        hn[] = 0
        if keys_ === nothing
            keys_ = sort!(collect(keys(r)))
            println(io, "day,", join(keys_, ","))
            r0 = r
        end
        println(io, d, ",", join((@sprintf("%.10g", r[k]) for k in keys_), ","))
        flush(io)
    end
    write_row(0)
    if get(ENV, "FATES_NICHE_RSTDEBUG", "") == "1"
        nd = Ref(0)
        _C.FATES_PARITY_HOOK[] = (s, bc, phase) -> begin
            phase == 10 || return nothing
            nd[] += 1
            nd[] == parse(Int, get(ENV, "FATES_NICHE_CENSUS_DAY", string(NM_RST + 1))) && println("  NANCENSUS day $(nd[]): ", join(nan_census(s), " "))
            return nothing
        end
    end
    fr = _C.ForcingReader()
    _C.forcing_reader_init!(fr, joinpath(FORCDIR, "clmforc.$(year(START)).nc"))
    spd = round(Int, 86400 / DTIME)
    nstep = 0; t0 = time()
    for d in 0:NM_DAYS-1
        for k in 0:spd-1
            for i in ((d == 0 && k == 0) ? (0:1) : (nstep+1:nstep+1))   # step-0 replay
                cur = START + Second(round(Int, (d * spd + k) * DTIME))
                tod = hour(cur)*3600 + minute(cur)*60 + second(cur)
                calday = dayofyear(cur) + tod / _C.SECSPDAY
                (declin, _e1)   = _C.compute_orbital(calday)
                (declinm1, _e2) = _C.compute_orbital(calday - DTIME / _C.SECSPDAY)
                _C.init_daylength!(inst.gridcell, declin, declinm1, _C.ORB_OBLIQR_DEFAULT, 1:ng)
                _C.read_forcing_step!(fr, inst.atm2lnd, cur, ng, nc)
                _C.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
                _C.clm_drv!(config, inst, filt, filt_ia, bounds, true, calday + DTIME/_C.SECSPDAY,
                    declin, declin, _C.ORB_OBLIQR_DEFAULT, false, false, "", false;
                    nstep=i, is_first_step=(i <= 1),
                    is_beg_curr_day=(i > 0 && tod == 0), is_end_curr_day=false,
                    is_beg_curr_year=false,
                    dtime=DTIME, mon=month(cur), day=day(cur), secs=tod,
                    jday=dayofyear(cur), year=year(cur), photosyns=inst.photosyns)
                nstep = i
                i >= 1 && accum_host!()
            end
        end
        write_row(d + 1)
        if string(d + 1) in split(get(ENV, "FATES_NICHE_TRACE", ""), ',')
            p = inst.fates.sites[1].oldest_patch; ip = 0
            while p !== nothing
                ip += 1; c = p.tallest
                while c !== nothing
                    @printf("  TRACE day=%d p=%d age=%.4f area=%.2f pft=%d n=%.5g dbh=%.5g h=%.4g layer=%d isnew=%d status=%d cd=%d coage=%.3g\n",
                            d + 1, ip, p.age, p.area, c.pft, c.n, c.dbh, c.height, c.canopy_layer,
                            c.isnew, c.status_coh, c.crowndamage, c.coage)
                    c = c.shorter
                end
                p = p.younger
            end
        end
        if d + 1 == NM_RST
            rst = replace(out, r"\.csv$" => "") * ".fates_rst.nc"
            fates_restart_roundtrip!(inst, rst)
            @printf("  [restart] wrote + re-read FATES restart %s after day %d\n", rst, d + 1)
        end
        (d + 1) % 5 == 0 && @printf("  day %d  %.0fs\n", d + 1, time() - t0)
    end
    _C.forcing_reader_close!(fr)
    close(io)
    println("  wrote $out")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main_niche(ARGS[1], ARGS[2]))
end
