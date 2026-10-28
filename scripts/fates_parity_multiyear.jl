#!/usr/bin/env julia
# ==========================================================================
# fates_parity_multiyear.jl — multi-year Julia run of the `fates_parity_1pt` case
# (Bow at Banff, 2002-07-01 start, dtime = 3600 s, use_bedrock = false, all 14 PFTs),
# writing one CSV row per model day of the demographic diagnostics that the matching
# Fortran CTSM-FATES daily history carries (FATES_NCOHORTS, FATES_RECRUITMENT_PF,
# FATES_NPLANT_PF, FATES_NPLANT_SZPF -> max occupied size class, FATES_VEGC,
# FATES_STRUCTC, FATES_NPP, ...). See docs/FATES_MULTIYEAR_ENERGY.md "Defect 2".
#
#   FATES_MY_DAYS=1460 FATES_MY_OUT=/path/out.csv \
#       julia --project=. scripts/fates_parity_multiyear.jl
#
# Reuses `build()` from scripts/fates_fortran_parity.jl (the flag-audited harness).
# The calendar is NOLEAP (Feb 29 skipped), matching the CESM run, and the forcing
# file rolls over each 1 January (clmforc.YYYY.nc), as DATM does.
# ==========================================================================
Base.include(@__MODULE__, joinpath(@__DIR__, "fates_fortran_parity.jl"))

const MY_DAYS = parse(Int, get(ENV, "FATES_MY_DAYS", "1460"))
const MY_OUT  = get(ENV, "FATES_MY_OUT", joinpath(@__DIR__, "fates_parity_multiyear_julia.csv"))
# optional per-cohort trace: comma-separated model days (1-based) to dump every cohort
const MY_TRACE = Set(parse.(Int, filter(!isempty, split(get(ENV, "FATES_MY_TRACE", ""), ','))))

"Noleap model day-start for elapsed day `d` (0-based) from START."
function noleap_day(d::Int)
    t = Date(START); k = 0
    while k < d
        t += Day(1)
        (month(t) == 2 && day(t) == 29) && (t += Day(1))
        k += 1
    end
    return t
end

function site_summary(site)
    np = _C.numpft[]
    ncoh = 0; npatch = 0; carbon = 0.0; structc = 0.0; leafc = 0.0
    maxdbh = 0.0; maxdbh_pft = 0; nplant = zeros(np); maxdbh_pf = zeros(np)
    npp = 0.0; gpp = 0.0
    cp = site.oldest_patch
    while cp !== nothing
        npatch += 1
        cc = cp.tallest
        while cc !== nothing
            ncoh += 1
            ft = cc.pft
            if cc.dbh > maxdbh; maxdbh = cc.dbh; maxdbh_pft = ft; end
            maxdbh_pf[ft] = max(maxdbh_pf[ft], cc.dbh)
            nplant[ft] += cc.n
            for org in (_C.leaf_organ, _C.fnrt_organ, _C.sapw_organ, _C.store_organ,
                        _C.struct_organ)
                carbon += _pool(cc, org) * cc.n
            end
            structc += _pool(cc, _C.struct_organ) * cc.n
            leafc   += _pool(cc, _C.leaf_organ) * cc.n
            isfinite(cc.npp_acc_hold) && (npp += cc.npp_acc_hold * cc.n)
            isfinite(cc.gpp_acc_hold) && (gpp += cc.gpp_acc_hold * cc.n)
            cc = cc.shorter
        end
        cp = cp.younger
    end
    return (; ncoh, npatch, carbon, structc, leafc, maxdbh, maxdbh_pft, nplant, maxdbh_pf, npp, gpp)
end

function cohort_trace(io, site, d)
    cp = site.oldest_patch; p = 0
    while cp !== nothing
        p += 1; c = 0; cc = cp.tallest
        _C.patch_pft_size_profile(cp)
        prof = join((@sprintf("%d/%d:%.4g", ft, z, cp.pft_agb_profile[ft, z])
                     for ft in axes(cp.pft_agb_profile, 1), z in axes(cp.pft_agb_profile, 2)
                     if cp.pft_agb_profile[ft, z] > 0), " ")
        @printf(io, "PATCH day=%d p=%d patchno=%d age=%.5g area=%.6g agbsum=%.5g dist=%s prof=[%s]\n",
                d, p, cp.patchno, cp.age, cp.area, sum(cp.pft_agb_profile),
                string(round.(cp.disturbance_rates, sigdigits=3)), prof)
        while cc !== nothing
            c += 1
            @printf(io, "TRACE day=%d p=%d c=%d pft=%d n=%.6g dbh=%.6g h=%.4g cl=%d status=%d leafc=%.5g fnrtc=%.5g sapwc=%.5g storec=%.5g structc=%.5g reproc=%.5g npp_hold=%.5g gpp_hold=%.5g ddbhdt=%.5g ctrim=%.4g isnew=%d\n",
                    d, p, c, cc.pft, cc.n, cc.dbh, cc.height, cc.canopy_layer, cc.status_coh,
                    _pool(cc, _C.leaf_organ), _pool(cc, _C.fnrt_organ), _pool(cc, _C.sapw_organ),
                    _pool(cc, _C.store_organ), _pool(cc, _C.struct_organ), _pool(cc, _C.repro_organ),
                    cc.npp_acc_hold, cc.gpp_acc_hold, cc.ddbhdt, cc.canopy_trim, cc.isnew ? 1 : 0)
            cc = cc.shorter
        end
        cp = cp.younger
    end
end

function main_multiyear()
    inst, bounds, filt, filt_ia, config, _tm = build()
    site = inst.fates.sites[1]
    np = _C.numpft[]
    ng, nc = bounds.endg, bounds.endc

    # recruitment captured from the site's per-PFT recruitment_rate accumulator between
    # the `integrate` (13) and `recruit` (14) phase boundaries of each daily step.
    rec_before = zeros(np); rec_today = zeros(np)
    _C.FATES_PARITY_HOOK[] = (s, bc, phase) -> begin
        s === site || return nothing
        if phase == 13
            rec_before .= s.recruitment_rate[1:np]
        elseif phase == 14
            rec_today .= s.recruitment_rate[1:np] .- rec_before
        end
        return nothing
    end

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

    io = open(MY_OUT, "w")
    tio = isempty(MY_TRACE) ? nothing : open(replace(MY_OUT, ".csv" => "_trace.txt"), "w")
    println(io, "day,date,ncoh,npatch,carbon,structc,leafc,maxdbh,maxdbh_pft,npp_d,gpp_d,recruit_n,",
            join(["nplant_$f" for f in 1:np], ","), ",", join(["recruit_$f" for f in 1:np], ","),
            ",", join(["maxdbh_$f" for f in 1:np], ","),
            ",cstatus,ncolddays,nchilldays,gdd,tveg24C,cndaysleafon,cndaysleafoff,doy,",
            join(["dstatus_$f" for f in 1:np], ","), ",", join(["elong_$f" for f in 1:np], ","))

    fr = _C.ForcingReader(); fyear = year(START)
    _C.forcing_reader_init!(fr, joinpath(FORCDIR, "clmforc.$(fyear).nc"))
    spd = round(Int, 86400 / DTIME)
    nstep = 0; t0 = time()
    for d in 0:MY_DAYS-1
        dstart = noleap_day(d)
        for k in 0:spd-1
            nsteps_this = (d == 0 && k == 0) ? (0:1) : (nstep+1:nstep+1)   # step-0 replay
            for i in nsteps_this
                cur = DateTime(dstart) + Second(round(Int, k * DTIME))
                if year(cur) != fyear
                    _C.forcing_reader_close!(fr); fyear = year(cur)
                    _C.forcing_reader_init!(fr, joinpath(FORCDIR, "clmforc.$(fyear).nc"))
                end
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
                    is_beg_curr_year=(i > 0 && tod == 0 && dayofyear(cur) == 1),
                    dtime=DTIME, mon=month(cur), day=day(cur), secs=tod,
                    jday=dayofyear(cur), year=year(cur), photosyns=inst.photosyns)
                nstep = i
            end
        end
        # the daily dynamics fire on the FIRST step of a day (tod == 0); the row for
        # model day d+1 is therefore written after day d's last step, i.e. it reflects
        # the dynamics that ran at the start of day d.
        sm = site_summary(site)
        @printf(io, "%d,%s,%d,%d,%.10g,%.10g,%.10g,%.8g,%d,%.8g,%.8g,%.8g,%s,%s,%s",
                d + 1, string(dstart), sm.ncoh, sm.npatch, sm.carbon, sm.structc, sm.leafc,
                sm.maxdbh, sm.maxdbh_pft, sm.npp * 1000, sm.gpp * 1000, sum(rec_today),
                join((@sprintf("%.8g", x) for x in sm.nplant), ","),
                join((@sprintf("%.8g", x) for x in rec_today), ","),
                join((@sprintf("%.6g", x) for x in sm.maxdbh_pf), ","))
        ncold = count(<(_C.ed_params().ED_val_phen_coldtemp), site.vegtemp_memory)
        @printf(io, ",%d,%d,%d,%.6g,%.5g,%d,%d,%d,%s,%s\n", site.cstatus, ncold, site.nchilldays,
                site.grow_deg_days, site.vegtemp_memory[1], site.cndaysleafon, site.cndaysleafoff,
                _C.hlm_day_of_year[], join(site.dstatus[1:np], ","),
                join((@sprintf("%.4g", x) for x in site.elong_factor[1:np]), ","))
        flush(io)
        (d + 1) in MY_TRACE && tio !== nothing && (cohort_trace(tio, site, d + 1); flush(tio))
        if (d + 1) % 30 == 0
            @printf("  day %5d %s  ncoh=%4d npatch=%d C=%9.2f maxdbh=%8.3f (pft %d)  %.0fs\n",
                    d + 1, string(dstart), sm.ncoh, sm.npatch, sm.carbon, sm.maxdbh,
                    sm.maxdbh_pft, time() - t0)
        end
        fill!(rec_today, 0.0)
    end
    _C.forcing_reader_close!(fr)
    _C.FATES_PARITY_HOOK[] = nothing
    close(io); tio === nothing || close(tio)
    println("  wrote $MY_OUT")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main_multiyear())
end
