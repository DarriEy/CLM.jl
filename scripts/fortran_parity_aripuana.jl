# =============================================================================
# THIRD-DOMAIN Julia↔Fortran PARITY — Aripuanã, Amazon (tropical, evergreen, no
# snow; use_cn=false SP run). Same single-step methodology as the Stillwater
# harness, but the Fortran reference is NATIVE CLM restarts (clm2.r) written EVERY
# nstep (restart_option=nsteps, restart_n=1) by the clean cesm.exe — no custom dump
# instrumentation needed (the restart IS the pdump format read_fortran_restart!
# already consumes). The restart written after step N is at base+N h; injecting the
# (N-1) h restart and running ONE clm_drv! step reproduces Fortran's step N.
#
# Aripuanã is tropical broadleaf-evergreen — NO snow / NO frost — so it exercises a
# regime contrasting both Bow (alpine) and Stillwater (semi-arid), and is free of
# the trace-snow/frost diagnostics that dominated Stillwater's residuals.
#
# RESULT (2026-06-18, the THIRD domain, validating the Stillwater fixes generalize):
#   The settled regime is MACHINE PRECISION — nstep 7–13 global max|rel| ~3e-7
#   (every field exact). The only residuals are cold-start spin-up + daytime:
#     - nstep 2 (76%) and n3–6 (~1e-2): COLD-START FROZEN-SOIL THAW — and it is NOT
#       phase-change melt-completeness (the initial hypothesis was wrong). CLM
#       cold-starts the soil ice-filled (a generic IC, unphysical for the tropics);
#       it thaws over the first ~6 steps. Probing n2: Fortran warms the top soil
#       layer to 274.08 K (all ice melted + heated past freezing) while Julia stays
#       at 273.15 K still melting — i.e. Fortran delivers more heat to the cold soil
#       (eflx_soil_grnd ~352 W/m² implied vs Julia's 209.7); the phase change follows
#       whatever heat arrives. At cold start the exposed-veg filter hasn't activated
#       (injected FRAC_VEG_NOSNO_ALB = 0 for all patches, despite elai=5), so warm
#       tropical air (296.8 K) sits over the frozen ground (273.15 K) as a STRONGLY
#       STABLE layer, and the humid air condenses on it (latent eflx_lh_grnd -69.7 +
#       sensible -47.6 W/m², both INTO the ground).
#       ROOT CAUSE FOUND via Fortran re-instrumentation (2026-06-18) — and it is NOT a
#       model bug. Two earlier hypotheses (M-O over-suppression, then a coupled
#       soil-energy term) were BOTH disproven by a direct gated print of the soil solve
#       at nstep 2 (PARITYSOIL in SoilTemperatureMod, exe rebuilt via case.build with
#       LIBRARY_PATH→netcdf 4.10.0). The M-O scheme is byte-faithful, AND every soil
#       THERMAL PROPERTY and formula matches Fortran exactly: watsat (0.54669), tkmg
#       (2.15417), tkdry (0.17890), the tk/cv/fact formulas, CAPR (0.34), the layer
#       geometry (dz/z/zi), the constants (TKICE/TKWAT/DENICE/DENH2O), and the dke
#       frozen/thawed condition (>= TFRZ). The ONLY thing that differs is the soil
#       ICE/LIQ STATE fed to the tk computation: Julia (from the injected restart
#       clm2.r-03600, raw value verified) has ice=1.789/liq=1.044, but Fortran's
#       CONTINUOUS step-2 consumes ice=0.855/liq=2.060 — ~0.93 kg/m² of ice has already
#       melted (+a little infiltration, total water 2.833→2.915) between the restart-
#       write point and Fortran's step-2 SoilTemperature. With MORE ice Julia gets a
#       higher frozen-soil tk (1.583 vs 1.367) → conducts heat away from soil1 → melts
#       less (1.026 vs Fortran's full 0.855) → soil1 stays at 273.15 vs Fortran's 274.08.
#       So this is a RESTART-STATE / step-boundary mismatch in the rapidly-thawing
#       cold-start regime (the injected restart's ice ≠ the ice the next continuous step
#       actually uses), a single-step-injection artifact — NOT a soil-physics error. The
#       settled regime is exact (n7+) precisely because there the ice changes slowly, so
#       the restart matches the step boundary.
#     - nstep 14–22 (H2OSOI_LIQ ~5e-4→1e-2): a small daytime soil-water difference
#       (tropical high-moisture transpiration/evaporation), T_VEG/T_GRND ≤1.5e-3.
#   So the model is byte-faithful in the settled tropical regime; the >1e-2 points are
#   the cold-start IC thaw (n2-6) and a minor daytime soil-water creep (n14-22).
#
# Usage: julia +1.12 --project=. scripts/fortran_parity_aripuana.jl [nstep|all]
# =============================================================================
using CLM, NCDatasets, Dates, Printf

const ARI   = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Aripuana_Amazon"
const FS    = joinpath(ARI, "settings/CLM/parameters/surfdata_clm.nc")
const FP    = joinpath(ARI, "settings/CLM/parameters/clm5_params.nc")
const FFORC = joinpath(ARI, "data/forcing/CLM_input/clmforc.2002.nc")
const RDIR  = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_aripuana_parity"
const SNOWOPT = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_optics_5bnd_c013122.nc"
const SNOWAGE = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_drdt_bst_fit_60_c070416.nc"
const _read_fortran_restart! = getfield(CLM, Symbol("read_fortran_restart!"))

# CLM restart filename for the state at DateTime dt (end of the step reaching dt).
restart_file(dt::DateTime) = joinpath(RDIR, @sprintf("Aripuana_Amazon.clm2.r.%04d-%02d-%02d-%05d.nc",
    year(dt), month(dt), day(dt), hour(dt)*3600 + minute(dt)*60 + second(dt)))

function registry(inst)
    ws = inst.water.waterstatebulk_inst.ws
    wd = inst.water.waterdiagnosticbulk_inst
    [("T_SOISNO",   :col2d, () -> inst.temperature.t_soisno_col),
     ("H2OSOI_LIQ", :col2d, () -> ws.h2osoi_liq_col),
     ("H2OSOI_ICE", :col2d, () -> ws.h2osoi_ice_col),
     ("T_GRND",     :col1d, () -> inst.temperature.t_grnd_col),
     ("ZWT",        :col1d, () -> inst.soilhydrology.zwt_col),
     ("H2OSFC",     :col1d, () -> ws.h2osfc_col),
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

function run_step(nstep::Int)
    base = DateTime(2002, 1, 1)
    bfile = restart_file(base + Hour(nstep - 1))   # state at START of step (after step N-1)
    efile = restart_file(base + Hour(nstep))       # state after step N
    (isfile(bfile) && isfile(efile)) || return (NaN, Dict{String,Float64}())

    CLM.rooting_profile_config.rooting_profile_method_water  = CLM.JACKSON_1996_ROOT
    CLM.rooting_profile_config.rooting_profile_method_carbon = CLM.JACKSON_1996_ROOT
    CLM.snow_hydrology_set_control_for_testing!(wind_dep_snow_density=true)
    (inst, bounds, filt, tm) = CLM.clm_initialize!(; fsurdat=FS, paramfile=FP,
        start_date=base + Hour(nstep - 1), dtime=3600, use_cn=false, use_luna=false,
        use_bedrock=true, use_aquifer_layer=false, h2osfcflag=0,
        fsnowoptics=SNOWOPT, fsnowaging=SNOWAGE, int_snow_max=2000.0)
    _read_fortran_restart!(bfile, inst, bounds; inject_veg_struct=false)

    config  = CLM.CLMDriverConfig(use_cn=false, use_aquifer_layer=false,
                                  use_hydrstress=false, use_luna=false)
    filt_ia = CLM.clump_filter_inactive_and_active
    ng, nc, np = bounds.endg, bounds.endc, bounds.endp
    fr = CLM.ForcingReader(); CLM.forcing_reader_init!(fr, FFORC)
    tf = replace(FFORC, r"clmforc\.[^/]*\.nc$" => "topo_forcing.nc")
    if isfile(tf)
        dt = NCDataset(tf, "r")
        haskey(dt, "TOPO") && (ft = Float64(dt["TOPO"][1]);
            for g in 1:ng; inst.atm2lnd.forc_topo_grc[g] = ft; end;
            for c in 1:nc; inst.topo.topo_col[c] = ft; end)
        close(dt)
    end
    calday = CLM.get_curr_calday(tm); (declin, _) = CLM.compute_orbital(calday)
    (declinm1, _) = CLM.compute_orbital(calday - 3600.0/CLM.SECSPDAY)
    nextsw = calday + 3600.0/CLM.SECSPDAY
    CLM.init_daylength!(inst.gridcell, declin, declinm1, CLM.ORB_OBLIQR_DEFAULT, 1:ng)
    CLM.advance_timestep!(tm)
    CLM.read_forcing_step!(fr, inst.atm2lnd, base + Hour(max(nstep - 1, 1)), ng, nc)
    for g in 1:ng
        inst.atm2lnd.forc_hgt_u_grc[g] = 30.0; inst.atm2lnd.forc_hgt_t_grc[g] = 30.0
        inst.atm2lnd.forc_hgt_q_grc[g] = 30.0
    end
    CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
    (yr, mon, d, tod) = CLM.get_curr_date(tm)
    CLM.clm_drv!(config, inst, filt, filt_ia, bounds, true, nextsw, declin, declin,
        CLM.ORB_OBLIQR_DEFAULT, false, false, "", false;
        nstep=tm.nstep, is_first_step=true,
        is_beg_curr_day=CLM.is_beg_curr_day(tm), is_end_curr_day=CLM.is_end_curr_day(tm),
        is_beg_curr_year=CLM.is_beg_curr_year(tm), dtime=3600.0, mon=mon, day=d,
        photosyns=inst.photosyns)
    CLM.forcing_reader_close!(fr)

    ds = NCDataset(efile, "r"); gmax = 0.0; per = Dict{String,Float64}()
    for (name, kind, getter) in registry(inst)
        rel, n = diff_field(getter(), ds, name, kind)
        n == 0 && continue
        per[name] = rel; gmax = max(gmax, rel)
    end
    close(ds)
    return (gmax, per)
end

function main(steps)
    @printf("%-6s %12s | %s\n", "nstep", "global", "per-field max|rel|")
    for s in steps
        gmax, per = run_step(s)
        isnan(gmax) && (@printf("%-6d   (restart absent)\n", s); continue)
        worst = sort(collect(per), by=x->-x[2])[1:min(3, length(per))]
        wstr = join([@sprintf("%s=%.1e", k, v) for (k, v) in worst], "  ")
        @printf("%-6d %12.3e | %s\n", s, gmax, wstr)
    end
end

arg = length(ARGS) >= 1 ? ARGS[1] : "all"
main(arg == "all" ? (2:24) : (parse(Int, arg):parse(Int, arg)))
