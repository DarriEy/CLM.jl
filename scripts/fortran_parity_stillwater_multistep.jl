# =============================================================================
# MULTI-STEP off-Bow Julia↔Fortran parity — Stillwater, Oklahoma (SP, use_cn=false).
#
# Injects the Fortran before_step_n1 dump ONCE, then FREE-RUNS N consecutive
# clm_drv! steps (no re-injection) and diffs the prognostic state against the
# Fortran after_hydrologydrainage_nN dump after every step. Unlike the single-step
# harness (scripts/fortran_parity_stillwater.jl, which re-injects each step), this
# measures ACCUMULATED divergence — it catches drift/accumulation bugs the
# single-step check cannot, and validates that the forcing-record + veg-structure
# fixes hold over the full night→day→night cycle (nstep 1..24 = 2002-01-01).
#
# Forcing is read at the START of each step interval (the step yielding nstep s
# reads the record at base+(s-1)h, clamped to rec 1; the clmforc file starts 01:00).
#
# Usage: julia +1.12 --project=. scripts/fortran_parity_stillwater_multistep.jl [NSTEPS]
# =============================================================================
using CLM, NCDatasets, Dates, Printf

const STW   = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Stillwater_Oklahoma"
const FS    = joinpath(STW, "settings/CLM/parameters/surfdata_clm.nc")
const FP    = joinpath(STW, "settings/CLM/parameters/clm5_params.nc")
const FFORC = joinpath(STW, "data/forcing/CLM_input/clmforc.2002.nc")
const DUMPS = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_stillwater_parity"
const SNOWOPT = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_optics_5bnd_c013122.nc"
const SNOWAGE = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_drdt_bst_fit_60_c070416.nc"
const _read_fortran_restart! = getfield(CLM, Symbol("read_fortran_restart!"))

function registry(inst)
    ws = inst.water.waterstatebulk_inst.ws
    wd = inst.water.waterdiagnosticbulk_inst
    [("T_SOISNO",   :col2d, () -> inst.temperature.t_soisno_col),
     ("H2OSOI_LIQ", :col2d, () -> ws.h2osoi_liq_col),
     ("H2OSOI_ICE", :col2d, () -> ws.h2osoi_ice_col),
     ("T_GRND",     :col1d, () -> inst.temperature.t_grnd_col),
     ("ZWT",        :col1d, () -> inst.soilhydrology.zwt_col),
     ("H2OSFC",     :col1d, () -> ws.h2osfc_col),
     ("SNOW_DEPTH", :col1d, () -> wd.snow_depth_col),
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

function diff_against_dump(inst, nstep)
    f = joinpath(DUMPS, "pdump_after_hydrologydrainage_n$(nstep).nc")
    isfile(f) || return (NaN, Dict{String,Float64}())
    ds = NCDataset(f, "r"); gmax = 0.0; per = Dict{String,Float64}()
    for (name, kind, getter) in registry(inst)
        rel, n = diff_field(getter(), ds, name, kind)
        n == 0 && continue
        per[name] = rel; gmax = max(gmax, rel)
    end
    close(ds)
    return (gmax, per)
end

function main(; nsteps::Int = 24)
    CLM.rooting_profile_config.rooting_profile_method_water  = CLM.JACKSON_1996_ROOT
    CLM.rooting_profile_config.rooting_profile_method_carbon = CLM.JACKSON_1996_ROOT
    base = DateTime(2002, 1, 1)
    (inst, bounds, filt, tm) = CLM.clm_initialize!(; fsurdat=FS, paramfile=FP,
        start_date=base, dtime=3600, use_cn=false, use_luna=false, use_bedrock=true,
        use_aquifer_layer=false, h2osfcflag=0, fsnowoptics=SNOWOPT, fsnowaging=SNOWAGE,
        int_snow_max=2000.0)
    ng, nc, np = bounds.endg, bounds.endc, bounds.endp
    @printf("Stillwater inst: nc=%d np=%d  PFTs=%s\n", nc, np, string(Int.(inst.patch.itype)))

    # Inject the Fortran state at the START of step 1 (nstep-0 state), ONCE.
    _read_fortran_restart!(joinpath(DUMPS, "pdump_before_step_n1.nc"), inst, bounds;
                           inject_veg_struct=false)

    config  = CLM.CLMDriverConfig(use_cn=false, use_aquifer_layer=false,
                                  use_hydrstress=false, use_luna=false)
    filt_ia = CLM.clump_filter_inactive_and_active
    fr = CLM.ForcingReader(); CLM.forcing_reader_init!(fr, FFORC)
    # Topo (matches the single-step harness).
    tf = replace(FFORC, r"clmforc\.[^/]*\.nc$" => "topo_forcing.nc")
    if isfile(tf)
        dt = NCDataset(tf, "r")
        haskey(dt, "TOPO") && (ft = Float64(dt["TOPO"][1]);
            for g in 1:ng; inst.atm2lnd.forc_topo_grc[g] = ft; end;
            for c in 1:nc; inst.topo.topo_col[c] = ft; end)
        close(dt)
    end

    @printf("\n%-6s %-8s %12s | %s\n", "nstep", "tod", "global", "per-field max|rel|")
    for s in 1:nsteps
        calday = CLM.get_curr_calday(tm); (declin, _) = CLM.compute_orbital(calday)
        (declinm1, _) = CLM.compute_orbital(calday - 3600.0/CLM.SECSPDAY)
        nextsw = calday + 3600.0/CLM.SECSPDAY
        CLM.init_daylength!(inst.gridcell, declin, declinm1, CLM.ORB_OBLIQR_DEFAULT, 1:ng)
        CLM.advance_timestep!(tm)
        # Forcing at the START of this step's interval (rec s-1, clamped to rec 1).
        force_date = base + Hour(max(s - 1, 1))
        CLM.read_forcing_step!(fr, inst.atm2lnd, force_date, ng, nc)
        for g in 1:ng
            inst.atm2lnd.forc_hgt_u_grc[g] = 30.0
            inst.atm2lnd.forc_hgt_t_grc[g] = 30.0
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

        gmax, per = diff_against_dump(inst, s)
        worst = sort(collect(per), by=x->-x[2])[1:min(3, length(per))]
        wstr = join([@sprintf("%s=%.1e", k, v) for (k, v) in worst], "  ")
        @printf("%-6d %-8d %12.3e | %s\n", s, tod, gmax, wstr)
    end
    CLM.forcing_reader_close!(fr)
end

main(; nsteps = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 24)
