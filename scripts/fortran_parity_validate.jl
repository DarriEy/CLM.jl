# ==========================================================================
# fortran_parity_validate.jl — single-step Julia ↔ Fortran parity from a
# SHARED initial condition.
#
# 1. Build a Bow single-column Julia CLM instance.
# 2. Inject the Fortran `before_step` dump as the shared IC.
# 3. Run clm_drv! for ONE timestep (forcing from clmforc.2003.nc at the
#    matching time — the same file Fortran's datm streams read).
# 4. Compare the resulting Julia state to the Fortran `after_hydrologydrainage`
#    dump (end-of-physics state for that step).
#
# This replaces the old free-running annual mean (which conflated per-routine
# translation error with year-long feedback accumulation) with a one-step,
# shared-IC comparison that localizes drift to a single timestep.
#
#   julia +1.12 --project=. scripts/fortran_parity_validate.jl [nstep]
# ==========================================================================

include(joinpath(@__DIR__, "fortran_parity_common.jl"))

const NSTEP        = length(ARGS) >= 1 ? ARGS[1] : "8761"
const DUMP_BEFORE  = joinpath(DUMPDIR, "pdump_before_step_n$(NSTEP).nc")
const DUMP_AFTER   = joinpath(DUMPDIR, "pdump_after_hydrologydrainage_n$(NSTEP).nc")
const DUMP_CANOPY  = joinpath(DUMPDIR, "pdump_after_canopyfluxes_n$(NSTEP).nc")

println("="^64)
println("  Julia ↔ Fortran parity — SINGLE-STEP, SHARED IC")
println("  nstep=$NSTEP   dtime=3600s   Bow-at-Banff (SP, use_cn=false)")
println("="^64)
for f in (DUMP_BEFORE, DUMP_AFTER); isfile(f) || error("dump missing: $f"); end

println("\n[1] Build Bow instance + inject shared IC (before_step) ...")
# Wall-clock date of this step's START: step 8761 begins 2003-01-01 00:00.
const STEP_DATE = DateTime(2003,1,1) + Hour(parse(Int, NSTEP) - 8761)
println("    step start date = ", STEP_DATE)
(inst, bounds, filt, tm) = build_bow_inst(; dtime=3600, start_date=STEP_DATE)
inject_dump!(inst, bounds, DUMP_BEFORE)

# ---- forcing + driver config (mirrors clm_run.jl loop body) ----
config   = CLM.CLMDriverConfig(use_cn=false, use_aquifer_layer=false)
filt_ia  = CLM.clump_filter_inactive_and_active
ng, nc, np = bounds.endg, bounds.endc, bounds.endp

fr = CLM.ForcingReader()
CLM.forcing_reader_init!(fr, FFORCING)

# topo for lapse correction (match clm_run.jl)
topo_file = replace(FFORCING, r"clmforc\.[^/]*\.nc$" => "topo_forcing.nc")
if isfile(topo_file)
    ds_topo = NCDataset(topo_file, "r")
    if haskey(ds_topo, "TOPO")
        ft = Float64(ds_topo["TOPO"][1])
        for g in 1:ng; inst.atm2lnd.forc_topo_grc[g] = ft; end
        for c in 1:nc; inst.topo.topo_col[c] = ft; end
    end
    close(ds_topo)
end

println("[2] Advance one step + read/downscale forcing ...")
step_start = tm.current_date   # forcing AND solar zenith both at the step START
# Orbital/cosz must use the SAME (step-start) time as the forcing, else at low
# sun a 1-hour offset swings absorbed solar and over-cools the ground.
calday = CLM.get_curr_calday(tm)
(declin, eccf) = CLM.compute_orbital(calday)
obliqr  = CLM.ORB_OBLIQR_DEFAULT
nextsw_cday = calday + 3600.0 / CLM.SECSPDAY

# Daylength is normally set in clm_initialize via init_daylength!; the real driver
# only RE-computes it when is_first_step==false. The teacher-forced harness resets
# nstep so is_first_step==true → update_daylength! would leave dayl=NaN, which makes
# dayl_factor → vcmax25 → photosynthesis → t_veg all NaN. Seed it from this step's
# declination (and the previous step's, for prev_dayl) to match real initialization.
(declinm1, _) = CLM.compute_orbital(calday - 3600.0 / CLM.SECSPDAY)
CLM.init_daylength!(inst.gridcell, declin, declinm1, obliqr, 1:bounds.endg)

CLM.advance_timestep!(tm)
CLM.read_forcing_step!(fr, inst.atm2lnd, step_start, ng, nc)
CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
(yr, mon, d, tod) = CLM.get_curr_date(tm)

println("    step_start = ", step_start, "  calday=", round(calday, digits=4))

println("[3] Run clm_drv! one step ...")
# Per-iteration leaf-temperature Newton dump (parity localization). Enable only
# when CANOPY_PERITER_DUMP env var is set, so the normal validate path is
# unchanged. Writes one line per (itlef, p) matching the instrumented Fortran.
if get(ENV, "CANOPY_PERITER_DUMP", "") != ""
    CLM.CANOPY_PERITER_DEBUG[] = true
    CLM.CANOPY_PERITER_NSTEP[] = parse(Int, NSTEP)
    CLM.CANOPY_PERITER_PATH[]  = joinpath(@__DIR__, "validation", "fortran_pdump",
                                          "periter_n$(NSTEP)_julia.txt")
    println("    per-iteration dump -> ", CLM.CANOPY_PERITER_PATH[])
end
CLM.clm_drv!(config, inst, filt, filt_ia, bounds,
             true, nextsw_cday, declin, declin, obliqr,
             false, false, "", false;
             nstep=tm.nstep, is_first_step=(tm.nstep==1),
             is_beg_curr_day=CLM.is_beg_curr_day(tm),
             is_end_curr_day=CLM.is_end_curr_day(tm),
             is_beg_curr_year=CLM.is_beg_curr_year(tm),
             dtime=3600.0, mon=mon, day=d, photosyns=inst.photosyns)
CLM.forcing_reader_close!(fr)

# Diagnostics most relevant to the summer soil-water/BTRAN divergence
println("\n  --- BTRAN / transpiration (summer water-stress) ---")
bt = inst.energyflux.btran_patch
qtran = inst.water.waterfluxbulk_inst.wf.qflx_tran_veg_patch
qet = inst.water.waterfluxbulk_inst.wf.qflx_evap_tot_patch
for p in 1:min(3, length(bt))
    @printf("    patch %d: btran=%.4f  qflx_tran_veg=%.4e  qflx_evap_tot=%.4e mm/s\n",
            p, bt[p], qtran[p], qet[p])
end
@printf("    (Fortran summer BTRANMN ~0.01; Julia high => soil-water/stress bug)\n")

println("\n  --- Surface energy budget (patch 1): Julia vs Fortran ---")
let dsE = NCDataset(DUMP_AFTER, "r")
    pairs = [
        ("SABG (abs solar grnd)", inst.solarabs.sabg_patch[1],            "SABG_P"),
        ("SABV (abs solar veg)",  inst.solarabs.sabv_patch[1],            "SABV_P"),
        ("eflx_sh_grnd",          inst.energyflux.eflx_sh_grnd_patch[1],  "EFLX_SHG_P"),
        ("eflx_lh_tot",           inst.energyflux.eflx_lh_tot_patch[1],   "EFLX_LH_P"),
        ("eflx_soil_grnd",        inst.energyflux.eflx_soil_grnd_patch[1],"EFLX_SOIG_P"),
        ("eflx_lwrad_net",        inst.energyflux.eflx_lwrad_net_patch[1], "EFLX_LWNET_P"),
        ("eflx_gnet",             inst.energyflux.eflx_gnet_patch[1],     "EFLX_GNET_P"),
        ("fv (ustar)",            inst.frictionvel.fv_patch[1],           "FV_P"),
        ("cgrnds (dSH/dTg)",      inst.energyflux.cgrnds_patch[1],        "CGRNDS_P"),
        ("cgrndl (dLH/dTg)",      inst.energyflux.cgrndl_patch[1],        "CGRNDL_P"),
    ]
    @printf("    %-22s %12s %12s %12s\n", "term [W/m2]", "Julia", "Fortran", "diff")
    for (nm, jv, fk) in pairs
        fv = haskey(dsE, fk) ? Float64(dsE[fk][1]) : NaN
        @printf("    %-22s %12.3f %12.3f %12.3f\n", nm, jv, fv, jv-fv)
    end
    close(dsE)
end

println("\n  --- Per-patch canopy energy/aero: Julia vs Fortran (after_canopyfluxes) ---")
let dsC = NCDataset(DUMP_CANOPY, "r")
    ef = inst.energyflux; fv = inst.frictionvel; te = inst.temperature
    sa = inst.solarabs;   cs = inst.canopystate; wf = inst.water.waterfluxbulk_inst.wf
    sb = inst.surfalb
    rows = [
        ("elai",            p->cs.elai_patch[p],            "elai"),
        ("frac_veg_nosno",  p->Float64(cs.frac_veg_nosno_patch[p]),     nothing),
        ("frac_vegnosno_alb",p->Float64(cs.frac_veg_nosno_alb_patch[p]),"FRAC_VEG_NOSNO_ALB"),
        ("fabd(band1)",     p->sb.fabd_patch[p,1],          nothing),
        ("fabi(band1)",     p->sb.fabi_patch[p,1],          nothing),
        ("T_VEG",           p->te.t_veg_patch[p],           "T_VEG"),
        ("taf (canopy air)",p->fv.taf_patch[p],             "taf"),
        ("SABV (solar veg)",p->sa.sabv_patch[p],            "SABV_P"),
        ("SABG (solar grnd)",p->sa.sabg_patch[p],           "SABG_P"),
        ("eflx_sh_veg",     p->ef.eflx_sh_veg_patch[p],     nothing),
        ("eflx_sh_grnd",    p->ef.eflx_sh_grnd_patch[p],    "EFLX_SHG_P"),
        ("eflx_lh_tot",     p->ef.eflx_lh_tot_patch[p],     "EFLX_LH_P"),
        ("eflx_gnet",       p->ef.eflx_gnet_patch[p],       "EFLX_GNET_P"),
        ("qflx_tran_veg",   p->wf.qflx_tran_veg_patch[p],   nothing),
        ("qflx_evap_veg",   p->wf.qflx_evap_veg_patch[p],   nothing),
        ("fv (ustar)",      p->fv.fv_patch[p],              "FV_P"),
        ("ram1",            p->fv.ram1_patch[p],            "RAM1_P"),
        ("obu",             p->fv.obu_patch[p],             "OBU"),
        ("forc_hgt_u",      p->fv.forc_hgt_u_patch[p],      "HGTU_P"),
        ("forc_hgt_t",      p->fv.forc_hgt_t_patch[p],      "HGTT_P"),
    ]
    npat = min(3, length(cs.elai_patch))
    @printf("    %-18s %30s   %30s\n", "", "Julia (p1/p2/p3)", "Fortran (p1/p2/p3)")
    for (nm, jf, fk) in rows
        jvals = [jf(p) for p in 1:npat]
        rdf(p) = (fk !== nothing && haskey(dsC, fk) && p <= length(dsC[fk]) && !ismissing(dsC[fk][p])) ?
                 Float64(dsC[fk][p]) : NaN
        fvals = [rdf(p) for p in 1:npat]
        js = join((@sprintf("%9.3f", v) for v in jvals), " ")
        fs = join((@sprintf("%9.3f", v) for v in fvals), " ")
        @printf("    %-18s %30s   %30s\n", nm, js, fs)
    end
    if hasproperty(fv, :num_iter_patch)
        @printf("    %-18s %30s   (itmax=%d)\n", "canopy num_iter",
                join((@sprintf("%9d", Int(fv.num_iter_patch[p])) for p in 1:npat), " "),
                CLM.canopy_fluxes_ctrl.itmax_canopy_fluxes)
    end
    close(dsC)
end

println("\n  --- Soil-water balance: per-layer ΔH2OSOI_LIQ (IC→after) Julia vs Fortran ---")
let dsB = NCDataset(DUMP_BEFORE, "r"), dsA = NCDataset(DUMP_AFTER, "r")
    ws = inst.water.waterstatebulk_inst.ws
    nlevsno = 12
    h0 = Float64.(dsB["H2OSOI_LIQ"][:, 1])   # injected IC (37)
    hF = Float64.(dsA["H2OSOI_LIQ"][:, 1])   # Fortran after (37)
    @printf("    %-5s %10s %10s %10s %10s %10s\n", "lay", "IC", "dJ", "dF", "dJ-dF", "J-F")
    for j in 1:18                              # soil layers 1..18 (model idx nlevsno+j)
        mi = nlevsno + j
        ic = h0[mi]; jf = ws.h2osoi_liq_col[1, mi]; ff = hF[mi]
        @printf("    %-5d %10.4f %10.4f %10.4f %10.3e %10.3e\n",
                j, ic, jf-ic, ff-ic, (jf-ic)-(ff-ic), jf-ff)
    end
    # Column water-table / surface-water state
    sh = inst.soilhydrology
    wtab = [
        ("ZWT",        sh.zwt_col[1],            "ZWT"),
        ("ZWT_PERCH",  sh.zwt_perched_col[1],    "ZWT_PERCH"),
        ("WA",         ws.wa_col[1],             "WA"),
        ("H2OSFC",     ws.h2osfc_col[1],         "H2OSFC"),
        ("FROST_TABLE",hasproperty(sh,:frost_table_col) ? sh.frost_table_col[1] : NaN, "FROST_TABLE"),
    ]
    println("    -- water-table / surface --")
    for (nm, jv, fk) in wtab
        ff = haskey(dsA, fk) ? Float64(dsA[fk][1]) : NaN
        ic = haskey(dsB, fk) ? Float64(dsB[fk][1]) : NaN
        @printf("    %-12s IC=%10.4f  J=%10.4f  F=%10.4f  J-F=%10.3e\n", nm, ic, jv, ff, jv-ff)
    end
    close(dsB); close(dsA)
end

println("\n[4] Compare Julia post-step state vs Fortran after_hydrologydrainage:")
results, gmax = compare_inst_to_dump(inst, DUMP_AFTER; label="1-step", tol=1e-9)

println()
println("="^64)
ndiff = count(r -> r[3] > 1e-6, results)
@printf("  global max|rel| after one step = %.3e\n", gmax)
println("  $ndiff / $(length(results)) fields diverge > 1e-6 from Fortran")
println("  (the first/largest divergence localizes the next module to inspect;")
println("   per-module boundary comparison is the next harness mode.)")
println("="^64)
