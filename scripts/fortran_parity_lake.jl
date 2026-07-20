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
# so a scalar field is f[1, rec] and a per-level field is f[1, j, rec], where
# rec = s + LAKE_REC_SHIFT maps model step s onto its h0 record.
#
# Fields compared: TG, TLAKE (×10 levels), LAKEICEFRAC (×10) + LAKEICEFRAC_SURF,
# EFLX_LH_TOT, FSH (sensible), FIRE (LW out), FSA (absorbed solar), TSA, H2OSNO.
#
# Note on averaging: the h0 records are time:mean over each 1-step hour, and we read
# the Julia state at end-of-step — with 1 step/record these coincide for state vars;
# flux fields (LH/SH/FIRE/FSA) are the step's instantaneous value vs a 1-sample mean.
#
# RESULT (2026-07-20, see docs/LAKE_FLUX_RESIDUAL.md): the long-standing "surface
# turbulent-flux residual" was substantially a HARNESS bug — this script compared
# Julia step s against h0 record s, but record 1 is the nstep=0 cold-start dump
# (time_bounds [0,0]), so model step s is record s+1. Fixing the index (LAKE_REC_SHIFT,
# default 1) takes max|rel| from 1.298e+01 to 3.108e-01. Run 47 steps, not 48.
#
# What SURVIVES the fix (both lake-gated, both still open):
#   - TSA is a DEAD WRITE: t_ref2m_patch is written only by bareground/urban fluxes,
#     never on a lake patch, so it sits at the cold-start 283.000 K for all 48 steps.
#   - FSH/EFLX_LH ~6-12%: lake_fluxes.jl:336 implements ONE of the four stability
#     regimes in FrictionVelocityMod.F90:1010-1050 (the very-unstable branch, which is
#     the one a warm winter lake under cold air actually needs, is missing).
#
# NOTE: the old note below about t_grnd being "BOUNDED to ~[248,261] K, no rah reaches
# Fortran's 271 K" was an artifact of the off-by-one — it compared against record 1,
# the initial dump. TG now tracks to ~0.3 K. Do not re-chase it.
#
# FIXED 2026-06-19 (lake_fluxes.jl): (1) the MoninObukIni obu had the WRONG SIGN for
# unstable conditions (obu = -ustar^3*thv/(vkc*g*thvstar) returned the stable branch
# when the lake is warmer than the air) -> ustar collapsed; replaced with the CLM
# bulk-Richardson MoninObukIni (obu<0 + convective um for unstable). (2) wind_min
# 0.1->1.0 (Fortran param). (3) the final fluxes used rah_final=zldis/(vkc*ustar),
# dropping the log(z/z0)+stability profile; now use the iteration's converged
# rah/raw/ram. These make iteration 1 correctly unstable and improve early FSH/LH
# (rel 0.94->0.83). Also added the fetch-limited Charnock cur (CURM 0->0.1 + the
# fetch/depth exp term; the constant CUR0 left z0mg too small) — correct but inert.
#
# REMAINING (and it is NOT the aerodynamic resistance): instrumenting the surface-
# temp Newton (ax/bx, formula identical to Fortran) shows t_grnd = ax/bx as a
# function of rah is BOUNDED to ~[248, 261] K — NO rah reaches Fortran's 271 K. The
# only term that can hold t_grnd near the warm lake is the conduction
# tksur*tsur/dzsur; reaching 271 needs tksur ~= 15, but tksur = savedtke1 = tkwat =
# 0.57 (molecular) at cold start (which MATCHES Fortran's cold-start init). So the
# residual is the lake SURFACE-WATER THERMAL COUPLING — the eddy conductivity /
# cold-start t_lake profile, or LakeTemperature resetting t_grnd toward the warm
# column post-flux. NEXT: instrument the Fortran lake_fluxes/lake_temperature
# (rebuild w/ prints) to compare tksur / savedtke1 / t_lake at step 1 and whether
# t_grnd is updated after the flux solve.
#
# Usage: julia +1.12 --project=. scripts/fortran_parity_lake.jl [NSTEPS] [--all]
# =============================================================================
using CLM, NCDatasets, Dates, Printf

const SHOWALL = any(==("--all"), ARGS)
# Hour shift applied to the forcing lookup (diagnostic; 0 = historical behaviour).
const LAKE_FORC_OFFSET = parse(Int, get(ENV, "LAKE_FORC_OFFSET", "0"))
# LAKE_FORC_EXACT=1 drops the historical `max(s-1, 1)` clamp and uses the
# mechanistically-derived mapping step s -> forcing hour (s-1). The clamp makes
# steps 1 AND 2 read hour 1, i.e. it corrupts step 1 and duplicates a record.
const LAKE_FORC_EXACT = get(ENV, "LAKE_FORC_EXACT", "") == "1"
# Model-step -> h0-record mapping. The h0 `nstep` axis of the lake reference is
# [0, 1, 2, ...] and its `time_bounds` are [[0,0], [0,1h], [1h,2h], ...]: record 1
# is the nstep=0 COLD-START dump over a ZERO-WIDTH interval, not a step average.
# So model step s is h0 record s+1, and the historical `[1, s]` indexing compared
# every Julia step against the Fortran state ONE STEP BEHIND it.
const LAKE_REC_SHIFT = parse(Int, get(ENV, "LAKE_REC_SHIFT", "1"))

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

"""
    main(; nsteps=48) -> NamedTuple or `missing`

Runs the lake parity comparison and prints the per-step table. ALSO returns the
raw per-step Julia and Fortran values (`jv` / `fv`, `Dict{String,Vector{Float64}}`
keyed by the scalar field names below) so a TEST can assert VALUES rather than
re-print relative diffs — `test_lake_ref2m.jl` uses this. Returns `missing` when
the machine-local reference inputs are absent, so callers can skip.
"""
function main(; nsteps::Int = 48)
    if !isfile(H0) || !isfile(LAKE) || !isfile(FP) || !isfile(FFORC)
        return missing
    end
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
    # Raw per-step values, returned for value-level assertions in the test suite.
    jv = Dict{String,Vector{Float64}}(nm => Float64[] for (nm, _, _) in scalars)
    fvv = Dict{String,Vector{Float64}}(nm => Float64[] for (nm, _, _) in scalars)
    # Converged lake zeta per step, and the Monin-Obukhov stability REGIME it
    # lands in. CTSM's profile relations branch four ways on zeta; the port used
    # to implement only regime 2, so "which regime does this reference actually
    # exercise" was unanswerable. zeta is now carried out on frictionvel.zeta_patch.
    zetas = Float64[]
    regimes = Int[]
    # Lake eddy-diffusivity depth-decay rate ks [1/m] (LakeFluxesMod.F90:755) and the
    # roughness/resistance the surface flux hangs on. ks was hardwired to 0 (the
    # gridcell latitude was never wired into the lake kernel), which makes
    # exp(-ks*z) == 1 at every depth: the wind-driven eddy diffusivity then never
    # decays and the whole 50 m column is coupled to the skin. Carried out per step
    # so the suite can assert its VALUE, and the freeze it enables, per step.
    kss = Float64[]; z0mgs = Float64[]; rahs = Float64[]
    zetat_h = 0.465   # temperature/humidity profile transition (FrictionVelocityMod)
    regime_of(z) = z < -zetat_h ? 1 : z < 0.0 ? 2 : z <= 1.0 ? 3 : 4
    for s in 1:nsteps
        calday = CLM.get_curr_calday(tm); (declin, _) = CLM.compute_orbital(calday)
        (declinm1, _) = CLM.compute_orbital(calday - 3600.0/CLM.SECSPDAY)
        nextsw = calday + 3600.0/CLM.SECSPDAY
        CLM.init_daylength!(inst.gridcell, declin, declinm1, CLM.ORB_OBLIQR_DEFAULT, 1:ng)
        CLM.advance_timestep!(tm)
        # Forcing-hour mapping. The historical expression is `Hour(max(s-1, 1))`,
        # inherited from the SINGLE-step stillwater probe where max(.,1) only guards
        # N=1; in a LOOP it makes steps 1 and 2 BOTH read hour 1 and lags every later
        # step by one hour. LAKE_FORC_OFFSET lets us test the alignment empirically
        # (a harness forcing misalignment was the whole of the #233 fire "residual").
        force_date = LAKE_FORC_EXACT ? base + Hour(s - 1 + LAKE_FORC_OFFSET) :
                                       base + Hour(max(s - 1, 1) + LAKE_FORC_OFFSET)
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

        zc = Float64(inst.frictionvel.zeta_patch[p_lake])
        push!(zetas, zc); push!(regimes, regime_of(zc))
        push!(kss, Float64(inst.lakestate.ks_col[c_lake]))
        push!(z0mgs, Float64(inst.frictionvel.z0mg_patch[p_lake]))
        # rah inverted from the port's own flux definition (LakeFluxesMod.F90:526),
        # the same inversion LAKE_AERO_PROBE applies to the Fortran h0 fields.
        let rho = Float64(inst.atm2lnd.forc_rho_downscaled_col[c_lake]),
            thm = Float64(inst.atm2lnd.forc_t_downscaled_col[c_lake]) + 0.0098*30.0,
            fsh = Float64(inst.energyflux.eflx_sh_tot_patch[p_lake]),
            tg  = Float64(inst.temperature.t_grnd_col[c_lake])
            push!(rahs, abs(fsh) > 1e-8 ? rho*1004.64*(tg - thm)/fsh : NaN)
        end

        per = Dict{String,Float64}()
        # scalar fields
        for (nm, getj, fn) in scalars
            haskey(fds, fn) || continue
            jval = Float64(getj()); fval = _fv(fds[fn][1, s + LAKE_REC_SHIFT])
            push!(jv[nm], jval); push!(fvv[nm], fval)
            per[nm] = _rel(jval, fval)
        end
        # per-level lake fields: max|rel| over the 10 lake layers
        for (nm, jarr, fn) in [("TLAKE", temp.t_lake_col, "TLAKE"),
                               ("LAKEICE", ls.lake_icefrac_col, "LAKEICEFRAC")]
            haskey(fds, fn) || continue
            mx = 0.0
            for j in 1:nlevlak
                r = _rel(Float64(jarr[c_lake, j]), _fv(fds[fn][1, j, s + LAKE_REC_SHIFT]))
                isnan(r) || (mx = max(mx, r))
            end
            per[nm] = mx
        end
        # Diagnostic probes (env-gated). LAKE_ICE_PROBE: Julia-vs-Fortran top/bottom ice fraction
        # + t_grnd + t_lake[1] (localizes the freeze trajectory). LAKE_AERO_PROBE: ustar/z0mg/z0hg
        # + FSH/LH (localizes the surface-flux residual = the lake-MO aerodynamic resistance).
        if get(ENV, "LAKE_ICE_PROBE", "") == "1" && s <= 16
            # NOTE: the per-level Fortran fields take the SAME record map as the
            # scalars (s + LAKE_REC_SHIFT). An earlier revision of this probe read
            # them at [.., s], i.e. one record behind everything else.
            r = s + LAKE_REC_SHIFT
            @printf("  [ice s=%2d] ice1 J=%.5f F=%.5f | tg J=%.3f F=%.3f | tlk1 J=%.3f F=%.3f | tlk2 J=%.3f F=%.3f | dzlk1=%.4f dzlk2=%.4f tke1=%.4f gnet=%.2f\n",
                s, ls.lake_icefrac_col[c_lake,1], _fv(fds["LAKEICEFRAC"][1,1,r]),
                temp.t_grnd_col[c_lake], _fv(fds["TG"][1, r]),
                temp.t_lake_col[c_lake,1], _fv(fds["TLAKE"][1,1,r]),
                temp.t_lake_col[c_lake,2], _fv(fds["TLAKE"][1,2,r]),
                inst.column.dz_lake[c_lake,1], inst.column.dz_lake[c_lake,2],
                ls.savedtke1_col[c_lake], ef.eflx_gnet_patch[p_lake])
        end
        # LAKE_AERO_PROBE: the Monin-Obukhov state the surface fluxes hang on.
        # Fortran does not write ustar/obu/zeta/rah to h0, but FSH, TG, TBOT, QBOT,
        # PBOT, WIND and TAUX ARE there, and the flux definitions invert exactly:
        #   rah = rho*cp*(TG - thm)/FSH        (LakeFluxesMod.F90:526)
        #   ram = -rho*u/TAUX                  (LakeFluxesMod.F90:699)
        # so the Fortran aerodynamic resistances are RECOVERABLE from the reference
        # and can be compared against the port's directly. That is the only way to
        # tell WHICH of the four stability regimes the reference actually ran in.
        if get(ENV, "LAKE_AERO_PROBE", "") == "1" && s <= 16
            fv = inst.frictionvel
            rec = s + LAKE_REC_SHIFT
            rgas = 287.058; cp = 1004.64
            pbot = _fv(fds["PBOT"][1, rec]); tbot = _fv(fds["TBOT"][1, rec])
            qbot = _fv(fds["QBOT"][1, rec]); tgF = _fv(fds["TG"][1, rec])
            fshF = _fv(fds["FSH"][1, rec]); tauxF = _fv(fds["TAUX"][1, rec])
            windF = _fv(fds["WIND"][1, rec])
            rhoF = pbot / (rgas * tbot * (1 + 0.61 * qbot))
            thmF = tbot + 0.0098 * 30.0
            rahF = rhoF * cp * (tgF - thmF) / fshF
            ramF = -rhoF * windF / tauxF
            # same inversion on the Julia side, from the port's own state
            tgJ = temp.t_grnd_col[c_lake]; fshJ = ef.eflx_sh_tot_patch[p_lake]
            rhoJ = inst.atm2lnd.forc_rho_downscaled_col[c_lake]
            thmJ = inst.atm2lnd.forc_t_downscaled_col[c_lake] + 0.0098 * 30.0
            rahJ = rhoJ * cp * (tgJ - thmJ) / fshJ
            @printf("  [aero s=%2d] zeta=%9.3f obu=%9.3f ustar=%.4f z0mg=%.3e z0hg=%.3e | rah J=%7.1f F=%7.1f | ram F=%7.1f | TG J=%.2f F=%.2f\n",
                s, fv.zeta_patch[p_lake], fv.obu_patch[p_lake], fv.ustar_patch[p_lake],
                fv.z0mg_patch[p_lake], fv.z0hg_patch[p_lake], rahJ, rahF, ramF, tgJ, tgF)
        end
        # LAKE_DUMP: raw Julia-vs-Fortran values (not relative diffs) for the fields
        # in the open surface-flux residual, every step. Relative diffs hide a value
        # that is CONSTANT because it is never written (see TSA).
        if get(ENV, "LAKE_DUMP", "") == "1"
            @printf("  [dump s=%2d] TSA J=%8.3f F=%8.3f | FSH J=%8.3f F=%8.3f | LH J=%8.3f F=%8.3f | TG J=%8.3f F=%8.3f\n",
                s, temp.t_ref2m_patch[p_lake], _fv(fds["TSA"][1, s + LAKE_REC_SHIFT]),
                ef.eflx_sh_tot_patch[p_lake], _fv(fds["FSH"][1, s + LAKE_REC_SHIFT]),
                ef.eflx_lh_tot_patch[p_lake], _fv(fds["EFLX_LH_TOT"][1, s + LAKE_REC_SHIFT]),
                temp.t_grnd_col[c_lake], _fv(fds["TG"][1, s + LAKE_REC_SHIFT]))
        end
        gmax = maximum(v for v in values(per) if !isnan(v); init=0.0)
        gmax_run = max(gmax_run, gmax)
        if s == 1 && SHOWALL
            @printf("  patches: active=%s itype=%s column=%s\n",
                    string(inst.patch.active[1:np]), string(Int.(inst.patch.itype[1:np])),
                    string(Int.(inst.patch.column[1:np])))
            for (nm, getj, fn) in scalars
                haskey(fds, fn) && @printf("    %-12s julia=%+11.4f  fortran=%+11.4f\n",
                                           nm, Float64(getj()), _fv(fds[fn][1, s + LAKE_REC_SHIFT]))
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
    # Which stability regime did each step actually run in? A fix to the missing
    # regimes proves nothing if the reference never leaves the one already there.
    counts = [count(==(r), regimes) for r in 1:4]
    @printf("zeta range over run: [%.3f, %.3f] | regime steps: 1(very unstable)=%d 2(unstable)=%d 3(stable)=%d 4(very stable)=%d\n",
            minimum(zetas), maximum(zetas), counts[1], counts[2], counts[3], counts[4])
    return (; nsteps, jv, fv = fvv, gmax_run, zetas, regimes, kss, z0mgs, rahs)
end

# Auto-run only when executed as a script; `include`d from a test, callers drive
# `main` themselves and read its return value.
if abspath(PROGRAM_FILE) == @__FILE__
    let na = filter(a -> occursin(r"^\d+$", a), ARGS)
        main(; nsteps = isempty(na) ? 48 : parse(Int, na[1]))
    end
end
