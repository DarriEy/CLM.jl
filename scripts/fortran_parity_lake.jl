# =============================================================================
# LAKE config Julia↔Fortran parity — Bow PCT_LAKE=100, cold start (SP, use_cn=false).
#
# First validation of the config-gated LAKE landunit (deep-lake thermodynamics +
# ice). Both models COLD-START from the same surfdata_lake100.nc at 2003-01-01 and
# free-run 48 hourly steps; we diff the Julia lake-column/patch state against the
# Fortran .h0 history time series (gridcell = the lake value since PCT_LAKE=100).
#
# Fortran reference (generated 2026-06-19, cold start, "SUCCESSFUL TERMINATION"):
#   /Users/.../SYMFLUENCE_data/clm_lake_run/Bow_at_Banff_lumped.clm2.h0.2003-01-01-00000.nc
# h0 vars are stored (time, [levlak,] lndgrid) → NCDatasets returns them reversed,
# so a scalar field is f[1, s] and a per-level field is f[1, j, s] (s = step/record).
#
# Fields compared: TG, TLAKE (×10 levels), LAKEICEFRAC (×10) + LAKEICEFRAC_SURF,
# EFLX_LH_TOT, FSH (sensible), FIRE (LW out), FSA (absorbed solar), TSA, H2OSNO.
#
# Note on averaging: the h0 records are time:mean over each 1-step hour, and we read
# the Julia state at end-of-step — with 1 step/record these coincide for state vars;
# flux fields (LH/SH/FIRE/FSA) are the step's instantaneous value vs a 1-sample mean.
#
# RESULT (48-step cold-start, 2026-06-19): the lake THERMODYNAMICS track Fortran
# (TLAKE ~1%, H2OSNO tight, FSA==0 at the winter-night start), but the lake SURFACE
# TURBULENT-FLUX solve is broken: at step 1 Julia FSH=-1.0 / EFLX_LH=-0.14 W/m2 vs
# Fortran +84.9 / +53.8 (collapse to ~0 — implied aerodynamic resistance rah ~1e4
# s/m, ~100-1000x too high), so the lake surface radiatively decouples and TG
# over-cools to 250 K vs Fortran 271 K in ONE step (FIRE = sigma*Tg^4 follows:
# 220 vs 301). Same CLASS as the documented bare-ground Monin-Obukhov ustar
# collapse. It compounds: the cold surface over-freezes the lake (LAKEICE rel
# drifts to ~0.38 by step 48) and skews daytime albedo/FSA (~0.44). NEXT: chase the
# lake surface MO/aerodynamic-resistance in lake_fluxes! (ustar/rah; cf. ws_col,
# ks_col, and the stability iteration) — the thermodynamics downstream are fine.
#
# Usage: julia +1.12 --project=. scripts/fortran_parity_lake.jl [NSTEPS] [--all]
# =============================================================================
using CLM, NCDatasets, Dates, Printf

const SHOWALL = any(==("--all"), ARGS)

const LAKE   = "/Users/darri.eythorsson/Library/CloudStorage/GoogleDrive-dareyt@gmail.com/My Drive/code/clm_ports/CLM.jl/test_inputs/lake/surfdata_lake100.nc"
const FP     = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/optimization/CLM/dds_run_1/final_evaluation/settings/CLM/parameters/clm5_params.nc"
const FFORC  = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/data/forcing/CLM_input/clmforc.2003.nc"
const H0     = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_lake_run/Bow_at_Banff_lumped.clm2.h0.2003-01-01-00000.nc"
const SNOWOPT = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_optics_5bnd_c013122.nc"
const SNOWAGE = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_drdt_bst_fit_60_c070416.nc"

# NaN/missing-aware relative diff used throughout (rel to 1+|fortran|, like the
# other parity harnesses).
_rel(j, f) = (isnan(j) || isnan(f)) ? NaN : abs(j - f) / (1 + abs(f))
_fv(x) = ismissing(x) ? NaN : Float64(x)

function main(; nsteps::Int = 48)
    isfile(H0) || error("Fortran lake h0 reference not found: $H0")
    base = DateTime(2003, 1, 1)
    (inst, bounds, filt, tm) = CLM.clm_initialize!(; fsurdat=LAKE, paramfile=FP,
        start_date=base, dtime=3600, use_cn=false, use_luna=false, use_bedrock=true,
        use_aquifer_layer=false, h2osfcflag=0, fsnowoptics=SNOWOPT, fsnowaging=SNOWAGE,
        int_snow_max=2000.0)
    ng, nc, np = bounds.endg, bounds.endc, bounds.endp
    nlevlak = CLM.varpar.nlevlak

    # Locate the lake column and a patch on it (PCT_LAKE=100 → gridcell == lake).
    c_lake = findfirst(filt.lakec)
    c_lake === nothing && error("no lake column found (filt.lakec all false)")
    p_lake = findfirst(p -> inst.patch.column[p] == c_lake, 1:np)
    @printf("Lake inst: nc=%d np=%d nlevlak=%d | lake col=%d patch=%d\n",
            nc, np, nlevlak, c_lake, something(p_lake, -1))

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

    fds = NCDataset(H0, "r")
    # accessors: Julia value (this step) and the matching Fortran h0 record value
    temp = inst.temperature; ls = inst.lakestate; ef = inst.energyflux
    sa = inst.solarabs; ws = inst.water.waterstatebulk_inst.ws; wd = inst.water.waterdiagnosticbulk_inst
    # scalar fields: (name, julia getter, fortran h0 name)
    scalars = [
        ("TG",          () -> temp.t_grnd_col[c_lake],         "TG"),
        ("TSA",         () -> temp.t_ref2m_patch[p_lake],      "TSA"),
        ("LAKEICE_SRF", () -> ls.lake_icefracsurf_col[c_lake], "LAKEICEFRAC_SURF"),
        ("EFLX_LH",     () -> ef.eflx_lh_tot_patch[p_lake],    "EFLX_LH_TOT"),
        ("FSH",         () -> ef.eflx_sh_tot_patch[p_lake],    "FSH"),
        ("FIRE",        () -> ef.eflx_lwrad_out_patch[p_lake], "FIRE"),
        ("FSA",         () -> sa.fsa_patch[p_lake],            "FSA"),
        ("H2OSNO",      () -> CLM.history_h2osno_total_col(inst)[c_lake], "H2OSNO"),
    ]

    @printf("\n%-5s %-7s %11s | %s\n", "nstep", "tod", "global", "per-field max|rel|")
    gmax_run = 0.0
    for s in 1:nsteps
        calday = CLM.get_curr_calday(tm); (declin, _) = CLM.compute_orbital(calday)
        (declinm1, _) = CLM.compute_orbital(calday - 3600.0/CLM.SECSPDAY)
        nextsw = calday + 3600.0/CLM.SECSPDAY
        CLM.init_daylength!(inst.gridcell, declin, declinm1, CLM.ORB_OBLIQR_DEFAULT, 1:ng)
        CLM.advance_timestep!(tm)
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

        per = Dict{String,Float64}()
        # scalar fields
        for (nm, getj, fn) in scalars
            haskey(fds, fn) || continue
            per[nm] = _rel(Float64(getj()), _fv(fds[fn][1, s]))
        end
        # per-level lake fields: max|rel| over the 10 lake layers
        for (nm, jarr, fn) in [("TLAKE", temp.t_lake_col, "TLAKE"),
                               ("LAKEICE", ls.lake_icefrac_col, "LAKEICEFRAC")]
            haskey(fds, fn) || continue
            mx = 0.0
            for j in 1:nlevlak
                r = _rel(Float64(jarr[c_lake, j]), _fv(fds[fn][1, j, s]))
                isnan(r) || (mx = max(mx, r))
            end
            per[nm] = mx
        end
        gmax = maximum(v for v in values(per) if !isnan(v); init=0.0)
        gmax_run = max(gmax_run, gmax)
        if s == 1 && SHOWALL
            @printf("  patches: active=%s itype=%s column=%s\n",
                    string(inst.patch.active[1:np]), string(Int.(inst.patch.itype[1:np])),
                    string(Int.(inst.patch.column[1:np])))
            for (nm, getj, fn) in scalars
                haskey(fds, fn) && @printf("    %-12s julia=%+11.4f  fortran=%+11.4f\n",
                                           nm, Float64(getj()), _fv(fds[fn][1, s]))
            end
        end
        nshow = SHOWALL ? length(per) : min(3, length(per))
        worst = sort(collect(per), by=x->-(isnan(x[2]) ? -1.0 : x[2]))[1:nshow]
        wstr = join([@sprintf("%s=%.1e", k, v) for (k, v) in worst], "  ")
        @printf("%-5d %-7d %11.3e | %s\n", s, tod, gmax, wstr)
    end
    CLM.forcing_reader_close!(fr)
    close(fds)
    @printf("\nRun max |rel| over %d steps: %.3e\n", nsteps, gmax_run)
end

let na = filter(a -> occursin(r"^\d+$", a), ARGS)
    main(; nsteps = isempty(na) ? 48 : parse(Int, na[1]))
end
