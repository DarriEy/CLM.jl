# =============================================================================
# COLD-START PARITY PROOF — Aripuanã, Amazon (tropical, SP, use_cn=false).
#
# Unlike the injection harnesses, this does a TRUE cold start (clm_initialize! with
# NO Fortran-state injection) and FREE-RUNS, diffing against the native Fortran
# cold-start restarts (clm2.r every nstep). It proves Julia's COLD-START PATH —
# the initial-condition construction AND the first-steps evolution — matches
# Fortran, rather than only reproducing an injected single step.
#
# Fortran reference: the same clm_aripuana_parity cold-start run (startup, 24
# nsteps, restart every step). Julia cold-starts identically (startup, no finidat)
# and is diffed against clm2.r at base+N h after each step N.
#
# This harness opts into CLM.coldstart_match_fortran!(true) — the gated cold-start
# mode that initializes the soil to the EXACT CLM5/Fortran values (soil 272 K → frozen,
# h2osoi_vol = 0.15 for j<=nbedrock, frozen layers all-ice) instead of Julia's default
# latitude-based / 0.75*watsat init.
#
# RESULT (2026-06-18): with the matching cold-start, the true cold-start free-run
# tracks Fortran to ~1e-2 in temperatures (T_SOISNO/T_GRND/T_VEG) and the DEEP soil
# (layers 3-12 exact); global max|rel| fell from ~370 (default Julia cold start) to
# ~0.2-1.5. The remaining residual is H2OSOI_ICE/LIQ in the top 1-2 soil layers — the
# cold-start frozen-soil THAW transient (the surface layers thaw rapidly over the first
# steps; the exact liquid/ice split is the same melt-completeness sensitivity examined
# in fortran_parity_aripuana.jl). The deep frozen profile is byte-faithful.
#
# Usage: julia +1.12 --project=. scripts/fortran_parity_aripuana_coldstart.jl [N]
# =============================================================================
using CLM, NCDatasets, Dates, Printf

const ARI   = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Aripuana_Amazon"
const FS    = joinpath(ARI, "settings/CLM/parameters/surfdata_clm.nc")
const FP    = joinpath(ARI, "settings/CLM/parameters/clm5_params.nc")
const FFORC = joinpath(ARI, "data/forcing/CLM_input/clmforc.2002.nc")
const RDIR  = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_aripuana_parity"
const SNOWOPT = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_optics_5bnd_c013122.nc"
const SNOWAGE = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_drdt_bst_fit_60_c070416.nc"

restart_file(dt::DateTime) = joinpath(RDIR, @sprintf("Aripuana_Amazon.clm2.r.%04d-%02d-%02d-%05d.nc",
    year(dt), month(dt), day(dt), hour(dt)*3600 + minute(dt)*60 + second(dt)))

function registry(inst)
    ws = inst.water.waterstatebulk_inst.ws
    [("T_SOISNO",   :col2d, () -> inst.temperature.t_soisno_col),
     ("H2OSOI_LIQ", :col2d, () -> ws.h2osoi_liq_col),
     ("H2OSOI_ICE", :col2d, () -> ws.h2osoi_ice_col),
     ("T_GRND",     :col1d, () -> inst.temperature.t_grnd_col),
     ("ZWT",        :col1d, () -> inst.soilhydrology.zwt_col),
     ("T_VEG",      :patch, () -> inst.temperature.t_veg_patch)]
end

function diff_field(jl, ds, name, kind)
    haskey(ds, name) || return (NaN, 0)
    f = ds[name]; mx = 0.0; n = 0
    if kind == :col1d
        v = ismissing(f[1]) ? NaN : Float64(f[1]); isnan(v) && return (NaN, 0)
        return (abs(Float64(jl[1]) - v) / (1 + abs(v)), 1)
    elseif kind == :col2d
        d = f[:, 1]
        for j in 1:min(length(d), size(jl, 2))
            fv = ismissing(d[j]) ? NaN : Float64(d[j]); jv = Float64(jl[1, j])
            (isnan(fv) || isnan(jv)) && continue
            mx = max(mx, abs(jv - fv) / (1 + abs(fv))); n += 1
        end
    else
        d = f[:]
        for p in 1:min(length(d), length(jl))
            fv = ismissing(d[p]) ? NaN : Float64(d[p]); jv = Float64(jl[p])
            (isnan(fv) || isnan(jv)) && continue
            mx = max(mx, abs(jv - fv) / (1 + abs(fv))); n += 1
        end
    end
    return (mx, n)
end

function main(nsteps::Int)
    CLM.rooting_profile_config.rooting_profile_method_water  = CLM.JACKSON_1996_ROOT
    CLM.rooting_profile_config.rooting_profile_method_carbon = CLM.JACKSON_1996_ROOT
    CLM.snow_hydrology_set_control_for_testing!(wind_dep_snow_density=true)
    CLM.coldstart_match_fortran!(true)   # use the Fortran-exact cold-start IC (272 K / vol 0.15)
    base = DateTime(2002, 1, 1)
    # TRUE cold start: startup, no finidat, no injection.
    (inst, bounds, filt, tm) = CLM.clm_initialize!(; fsurdat=FS, paramfile=FP,
        start_date=base, dtime=3600, use_cn=false, use_luna=false, use_bedrock=true,
        use_aquifer_layer=false, h2osfcflag=0, fsnowoptics=SNOWOPT, fsnowaging=SNOWAGE,
        int_snow_max=2000.0)
    ng, nc, np = bounds.endg, bounds.endc, bounds.endp
    config  = CLM.CLMDriverConfig(use_cn=false, use_aquifer_layer=false,
                                  use_hydrstress=false, use_luna=false)
    filt_ia = CLM.clump_filter_inactive_and_active
    fr = CLM.ForcingReader(); CLM.forcing_reader_init!(fr, FFORC)
    tf = replace(FFORC, r"clmforc\.[^/]*\.nc$" => "topo_forcing.nc")
    if isfile(tf)
        dt = NCDataset(tf, "r")
        haskey(dt, "TOPO") && (ft = Float64(dt["TOPO"][1]);
            for g in 1:ng; inst.atm2lnd.forc_topo_grc[g] = ft; end;
            for c in 1:nc; inst.topo.topo_col[c] = ft; end)
        close(dt)
    end

    @printf("%-6s %12s | %s\n", "nstep", "global", "per-field max|rel|")
    for s in 1:nsteps
        calday = CLM.get_curr_calday(tm); (declin, _) = CLM.compute_orbital(calday)
        (declinm1, _) = CLM.compute_orbital(calday - 3600.0/CLM.SECSPDAY)
        nextsw = calday + 3600.0/CLM.SECSPDAY
        CLM.init_daylength!(inst.gridcell, declin, declinm1, CLM.ORB_OBLIQR_DEFAULT, 1:ng)
        CLM.advance_timestep!(tm)
        CLM.read_forcing_step!(fr, inst.atm2lnd, base + Hour(max(s - 1, 1)), ng, nc)
        for g in 1:ng
            inst.atm2lnd.forc_hgt_u_grc[g] = 30.0; inst.atm2lnd.forc_hgt_t_grc[g] = 30.0
            inst.atm2lnd.forc_hgt_q_grc[g] = 30.0
        end
        CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
        (yr, mon, d, tod) = CLM.get_curr_date(tm)
        CLM.clm_drv!(config, inst, filt, filt_ia, bounds, true, nextsw, declin, declin,
            CLM.ORB_OBLIQR_DEFAULT, false, false, "", false;
            nstep=tm.nstep, is_first_step=(s == 1),
            is_beg_curr_day=CLM.is_beg_curr_day(tm), is_end_curr_day=CLM.is_end_curr_day(tm),
            is_beg_curr_year=CLM.is_beg_curr_year(tm), dtime=3600.0, mon=mon, day=d,
            photosyns=inst.photosyns)

        efile = restart_file(base + Hour(s))
        if !isfile(efile); @printf("%-6d  (restart absent)\n", s); continue; end
        ds = NCDataset(efile, "r"); gmax = 0.0; per = Dict{String,Float64}()
        for (name, kind, getter) in registry(inst)
            rel, n = diff_field(getter(), ds, name, kind)
            n == 0 && continue
            per[name] = rel; gmax = max(gmax, rel)
        end
        close(ds)
        worst = sort(collect(per), by=x->-x[2])[1:min(3, length(per))]
        wstr = join([@sprintf("%s=%.1e", k, v) for (k, v) in worst], "  ")
        @printf("%-6d %12.3e | %s\n", s, gmax, wstr)
    end
    CLM.forcing_reader_close!(fr)
end

main(length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 24)
