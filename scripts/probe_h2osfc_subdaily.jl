# probe_h2osfc_subdaily.jl — dump per-timestep surface-water state for a domain.
#
# Confirms the FH2OSFC diurnal-cycle hypothesis: records h2osfc, frac_h2osfc and
# the surface-water in/out fluxes every timestep (not daily-averaged) for the
# soil column, so we can see whether Julia's within-day h2osfc is flat (steady
# ponding) vs cycles diurnally like Fortran.
#
#   DOMAIN=Krycklan julia +1.12 --project=. scripts/probe_h2osfc_subdaily.jl
#
# Output: scripts/../paper/data/probe_h2osfc_<domain>.csv
# Reuses the parity_run_domain.jl per-domain registry + config.

using NCDatasets, Dates, Printf, CLM
const run_clm! = getfield(CLM, Symbol("clm_run!"))
const DATA = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data"

_snicar(f) = let d1 = "$DATA/installs/cesm-inputdata/lnd/clm2/snicardata",
                 d2 = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata"
    isfile(joinpath(d1, f)) ? joinpath(d1, f) : joinpath(d2, f)
end

# Minimal per-domain registry (subset of parity_run_domain.jl fields we need).
const DOMAINS = Dict(
    "Krycklan" => (year = 2013, dtime = 3600, baseflow = 0.001, int_snow = 2000.0,
        caldir  = "$DATA/domain_Boreal_Krycklan_Sweden/settings/CLM/parameters",
        forcing = "$DATA/domain_Boreal_Krycklan_Sweden/data/forcing/CLM_input/clmforc.2013.nc",
        restart = "$DATA/domain_Boreal_Krycklan_Sweden/simulations/clm_boreal/CLM/Boreal_Krycklan_Sweden.clm2.r.2013-01-01-00000.nc"),
    "Abisko" => (year = 2013, dtime = 3600, baseflow = 0.001, int_snow = 2000.0,
        caldir  = "$DATA/domain_Arctic_Abisko_Sweden/settings/CLM/parameters",
        forcing = "$DATA/domain_Arctic_Abisko_Sweden/data/forcing/CLM_input/clmforc.2013.nc",
        restart = "$DATA/domain_Arctic_Abisko_Sweden/simulations/clm_arctic/CLM/Arctic_Abisko_Sweden.clm2.r.2013-01-01-00000.nc"),
    "Baltimore" => (year = 2013, dtime = 3600, baseflow = 0.001, int_snow = 2000.0,
        caldir  = "$DATA/domain_Urban_DeadRun_Baltimore/settings/CLM/parameters",
        forcing = "$DATA/domain_Urban_DeadRun_Baltimore/data/forcing/CLM_input/clmforc.2013.nc",
        restart = "$DATA/domain_Urban_DeadRun_Baltimore/simulations/clm_urban/CLM/Urban_DeadRun_Baltimore.clm2.r.2013-01-01-00000.nc"),
)

const DOM = get(ENV, "DOMAIN", "Krycklan")
cfg = DOMAINS[DOM]
yr  = cfg.year
fsurdat   = joinpath(cfg.caldir, "surfdata_clm.nc")
paramfile = joinpath(cfg.caldir, "clm5_params.nc")
outdir    = abspath(joinpath(@__DIR__, "..", "paper", "data"))
fhistory  = joinpath(outdir, "probe_hist_$(lowercase(DOM)).nc")

CLM.rooting_profile_config.rooting_profile_method_water  = CLM.JACKSON_1996_ROOT
CLM.rooting_profile_config.rooting_profile_method_carbon = CLM.JACKSON_1996_ROOT
CLM.snow_hydrology_set_control_for_testing!(;
    wind_dep_snow_density = true,
    overburden_compaction_method = CLM.OVERBURDEN_COMPACTION_VIONNET2012)

# Per-timestep records for the (first) soil column.
recs = NamedTuple[]
csoil = Ref(0)
pveg  = Ref(0)   # dominant vegetated patch (max wtgcell, itype != noveg)

function probe(inst, tm)
    # Locate a soil column once (istsoil landunit).
    if csoil[] == 0
        for c in eachindex(inst.column.snl)
            csoil[] = c; break
        end
    end
    c = csoil[]
    ws  = inst.water.waterstatebulk_inst.ws
    wd  = inst.water.waterdiagnosticbulk_inst
    wfb = inst.water.waterfluxbulk_inst
    wf  = wfb.wf
    getc(a) = (c <= length(a) ? Float64(a[c]) : NaN)
    # total snow water = unresolved + explicit layer ice+liq
    nsno = CLM.varpar.nlevsno
    h2osno_t = getc(ws.h2osno_no_layers_col)
    if c <= size(ws.h2osoi_ice_col, 1)
        for j in 1:nsno
            vi = ws.h2osoi_ice_col[c, j]; vl = ws.h2osoi_liq_col[c, j]
            isfinite(vi) && (h2osno_t += vi); isfinite(vl) && (h2osno_t += vl)
        end
    end
    snlc = c <= length(inst.column.snl) ? Int(inst.column.snl[c]) : 0
    # Locate dominant vegetated patch once (max wtgcell among itype != noveg=0).
    if pveg[] == 0
        want_it = parse(Int, get(ENV, "PVEG_ITYPE", "0"))  # 0 => max-wtgcell veg patch
        best = 0.0
        for pp in eachindex(inst.patch.itype)
            it = inst.patch.itype[pp]
            ok = want_it == 0 ? (it != 0) : (it == want_it)
            if ok && inst.patch.wtgcell[pp] > best
                best = inst.patch.wtgcell[pp]; pveg[] = pp
            end
        end
        pveg[] == 0 && (pveg[] = 1)
    end
    pv = pveg[]
    fv = inst.frictionvel; ps = inst.photosyns; cs = inst.canopystate; sa = inst.solarabs; alb = inst.surfalb
    getp(a) = (pv <= length(a) ? Float64(a[pv]) : NaN)
    getpz1(a) = (size(a,1) >= pv && size(a,2) >= 1 ? Float64(a[pv,1]) : NaN)  # canopy layer 1
    push!(recs, (
        date   = tm.current_date,
        h2osfc = getc(ws.h2osfc_col),
        frac   = getc(wd.frac_h2osfc_col),
        frac_nosnow = getc(wd.frac_h2osfc_nosnow_col),   # pre-snow-clamp frac
        frac_sno = getc(wd.frac_sno_col),                # snow cover fraction
        frac_sno_eff = getc(wd.frac_sno_eff_col),        # effective snow cover fraction
        h2osno = h2osno_t,                               # total SWE (mm)
        snl = Float64(snlc),                             # # snow layers (negative)
        snowdp = getc(wd.snow_depth_col),                # snow depth (m; >0 => snow present)
        q_in   = getc(wfb.qflx_top_soil_to_h2osfc_col),  # input to h2osfc (mm/s)
        q_drain= getc(wfb.qflx_h2osfc_drain_col),        # h2osfc -> soil (mm/s)
        q_spill= getc(wfb.qflx_h2osfc_surf_col),         # h2osfc -> runoff (mm/s)
        q_melt = getc(wf.qflx_snomelt_col),              # snowmelt (mm/s)
        q_infl = getc(wf.qflx_infl_col),                 # total infiltration (mm/s)
        # --- dominant veg-patch canopy internals (FCTR chase) ---
        uaf   = getp(fv.uaf_patch),
        qaf   = getp(fv.qaf_patch),
        rb    = getp(fv.rb1_patch),
        rssun = getp(ps.rssun_patch),
        rssha = getp(ps.rssha_patch),
        fdry  = getp(wd.fdry_patch),
        laisun= getp(cs.laisun_patch),
        laisha= getp(cs.laisha_patch),
        elai  = getp(cs.elai_patch),
        esai  = getp(cs.esai_patch),
        qtran = getp(wf.qflx_tran_veg_patch),
        parsun = getpz1(sa.parsun_z_patch),   # absorbed PAR per sunlit leaf area (W/m2)
        parsha = getpz1(sa.parsha_z_patch),   # absorbed PAR per shaded leaf area
        psnsun = getp(ps.psnsun_patch),
        psnsha = getp(ps.psnsha_patch),
        gsmolsun = getpz1(ps.gs_mol_sun_patch),
        gsmolsha = getpz1(ps.gs_mol_sha_patch),
        vcxcsun = getp(alb.vcmaxcintsun_patch),
        vcxcsha = getp(alb.vcmaxcintsha_patch),
        cisha = getpz1(ps.cisha_z_patch),    # shaded intercellular CO2 (Pa)
        ansha = getpz1(ps.an_sha_patch),     # shaded net assimilation (umol/m2/s)
        cisun = getpz1(ps.cisun_z_patch),
        ansun = getpz1(ps.an_sun_patch),
        lmrsha = getpz1(ps.lmrsha_z_patch),  # shaded leaf dark respiration
        lmrsun = getpz1(ps.lmrsun_z_patch),
        vcmx25z = getpz1(ps.vcmx25_z_patch),  # base LUNA vcmax25 (feeds lmr25)
        jmx25z = getpz1(ps.jmx25_z_patch),
        tlaiz1 = getpz1(alb.tlai_z_patch),    # per-layer-1 tlai increment (feeds lmr reduction!)
        nrad = getp(alb.nrad_patch),          # # radiative canopy layers
        tlaiz = getp(cs.tlai_patch),          # total LAI
    ))
end

println("Probe run: $DOM $yr (per-timestep surface water)")
t0 = time()
run_clm!(;
    fsurdat = fsurdat, paramfile = paramfile, fforcing = cfg.forcing,
    fhistory = fhistory,
    start_date = DateTime(yr, 1, 1), end_date = DateTime(yr + 1, 1, 1),
    dtime = cfg.dtime, use_cn = false, verbose = false,
    use_aquifer_layer = false, use_hydrstress = true, use_luna = true,
    h2osfcflag = 1,
    baseflow_scalar = cfg.baseflow, int_snow_max = cfg.int_snow,
    ffortran_restart = cfg.restart, interp_forcing = true,
    fsnowoptics = isfile(_snicar("snicar_optics_5bnd_c013122.nc")) ? _snicar("snicar_optics_5bnd_c013122.nc") : "",
    fsnowaging  = isfile(_snicar("snicar_drdt_bst_fit_60_c070416.nc")) ? _snicar("snicar_drdt_bst_fit_60_c070416.nc") : "",
    step_probe = probe)

csv = joinpath(outdir, "probe_h2osfc_$(lowercase(DOM)).csv")
open(csv, "w") do io
    println(io, "date,h2osfc,frac,frac_nosnow,frac_sno,frac_sno_eff,h2osno,snl,snowdp,q_in,q_drain,q_spill,q_melt,q_infl,uaf,qaf,rb,rssun,rssha,fdry,laisun,laisha,elai,esai,qtran,parsun,parsha,psnsun,psnsha,gsmolsun,gsmolsha,vcxcsun,vcxcsha,cisha,ansha,cisun,ansun,lmrsha,lmrsun,vcmx25z,jmx25z,tlaiz1,nrad,tlaiz")
    for r in recs
        @printf(io, "%s,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n",
            r.date, r.h2osfc, r.frac, r.frac_nosnow, r.frac_sno, r.frac_sno_eff, r.h2osno, r.snl, r.snowdp,
            r.q_in, r.q_drain, r.q_spill, r.q_melt, r.q_infl,
            r.uaf, r.qaf, r.rb, r.rssun, r.rssha, r.fdry, r.laisun, r.laisha, r.elai, r.esai, r.qtran,
            r.parsun, r.parsha, r.psnsun, r.psnsha, r.gsmolsun, r.gsmolsha, r.vcxcsun, r.vcxcsha,
            r.cisha, r.ansha, r.cisun, r.ansun, r.lmrsha, r.lmrsun, r.vcmx25z, r.jmx25z, r.tlaiz1, r.nrad, r.tlaiz)
    end
end
@printf("Done %s in %.1fs -> %s (%d steps)\n", DOM, time() - t0, csv, length(recs))
