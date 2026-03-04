using Dates, Printf, CLM

basedir = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped"
init_fn = getfield(CLM, Symbol("clm_initialize!"))
(inst, bounds, filt, tm) = init_fn(;
    fsurdat=joinpath(basedir, "settings/CLM/parameters/surfdata_clm.nc"),
    paramfile=joinpath(basedir, "settings/CLM/parameters/clm5_params.nc"),
    start_date=DateTime(2002,7,1), dtime=1800)

nc = bounds.endc; np = bounds.endp
ss = inst.soilstate
nlevsno = CLM.varpar.nlevsno
nlevgrnd = CLM.varpar.nlevgrnd
nlevsoi = CLM.varpar.nlevsoi
joff = nlevsno

# Check eff_porosity after cold start
println("=== EFF_POROSITY after cold start (col 1, first 5 layers) ===")
for j in 1:5
    @printf("  Layer %d: eff_porosity=%.6f, watsat=%.6f\n",
            j, ss.eff_porosity_col[1, j], ss.watsat_col[1, j])
end

# Check h2osoi_liq and h2osoi_liqvol
wsb = inst.water.waterstatebulk_inst.ws
wdb = inst.water.waterdiagnosticbulk_inst
println("\n=== SOIL WATER after cold start (col 1, first 5 layers) ===")
for j in 1:5
    liq = wsb.h2osoi_liq_col[1, j + joff]
    ice = wsb.h2osoi_ice_col[1, j + joff]
    dz = inst.column.dz[1, j + joff]
    liqvol = liq / (dz * CLM.DENH2O)
    @printf("  Layer %d: h2osoi_liq=%.4f, h2osoi_ice=%.4f, dz=%.4f, liqvol=%.6f\n",
            j, liq, ice, dz, liqvol)
end

# Check soil temperatures
temp = inst.temperature
println("\n=== SOIL TEMPS after cold start (col 1, first 5 layers) ===")
for j in 1:5
    @printf("  Layer %d: t_soisno=%.2f K (TFRZ-2=%.2f)\n",
            j, temp.t_soisno_col[1, j + joff], CLM.TFRZ - 2.0)
end

# Check root fractions
println("\n=== ROOT FRACTIONS (patch 2 = PFT 1, first 5 layers) ===")
for j in 1:5
    @printf("  Layer %d: rootfr=%.6f\n", j, ss.rootfr_patch[2, j])
end

# Check pftcon smpso/smpsc for the relevant PFTs
println("\n=== PFT WATER STRESS PARAMS ===")
for p in 1:np
    itype = inst.patch.itype[p]
    idx = itype + 1  # 0-based → 1-based
    if idx >= 1 && idx <= length(CLM.pftcon.smpso)
        @printf("Patch %d: itype=%d, idx=%d, smpso=%.1f, smpsc=%.1f\n",
                p, itype, idx, CLM.pftcon.smpso[idx], CLM.pftcon.smpsc[idx])
    end
end

# Check filter
println("\n=== EXPOSED VEG FILTER ===")
for p in 1:np
    @printf("Patch %d: exposedvegp=%s, noexposedvegp=%s\n",
            p, filt.exposedvegp[p], filt.noexposedvegp[p])
end

# Run 48 steps (1 day in July)
fr = CLM.ForcingReader()
CLM.forcing_reader_init!(fr, joinpath(basedir, "data/forcing/CLM_input/clmforc.2002.nc"))
config = CLM.CLMDriverConfig()
filt_ia = CLM.clump_filter_inactive_and_active

for step in 1:48
    CLM.advance_timestep!(tm)
    CLM.read_forcing_step!(fr, inst.atm2lnd, tm.current_date, bounds.endg, nc)
    CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
    calday = CLM.get_curr_calday(tm)
    (declin, eccf) = CLM.compute_orbital(calday)
    (yr, mon, d, tod) = CLM.get_curr_date(tm)
    CLM.clm_drv!(config, inst, filt, filt_ia, bounds,
                 true, calday + 1800.0/86400, declin, declin, CLM.ORB_OBLIQR_DEFAULT,
                 false, false, "", false; nstep=tm.nstep,
                 is_first_step=(tm.nstep==1), is_beg_curr_day=CLM.is_beg_curr_day(tm),
                 is_beg_curr_year=CLM.is_beg_curr_year(tm), dtime=1800.0,
                 mon=mon, day=d, photosyns=inst.photosyns)
end

println("\n=== AFTER 48 STEPS (1 day in July) ===")
println("BTRAN values:")
for p in 1:np
    @printf("  Patch %d (PFT %d): btran=%.8f\n", p, inst.patch.itype[p], inst.energyflux.btran_patch[p])
end

println("\nSoil temperatures (col 1, first 5 layers):")
for j in 1:5
    @printf("  Layer %d: t_soisno=%.2f K\n", j, temp.t_soisno_col[1, j + joff])
end

println("\nh2osoi_liqvol (col 1, first 5 layers):")
for j in 1:5
    @printf("  Layer %d: h2osoi_liqvol=%.6f\n", j, wdb.h2osoi_liqvol_col[1, j + joff])
end

println("\neff_porosity (col 1, first 5 layers):")
for j in 1:5
    @printf("  Layer %d: eff_porosity=%.6f\n", j, ss.eff_porosity_col[1, j])
end

println("\nrootr (patch 2, first 5 layers):")
for j in 1:5
    @printf("  Layer %d: rootr=%.8f\n", j, ss.rootr_patch[2, j])
end

println("\n=== MANUAL BTRAN COMPUTATION (patch 2, layer 1) ===")
p = 2; c = inst.patch.column[p]; j = 1
itype = inst.patch.itype[p] + 1
liqvol = wdb.h2osoi_liqvol_col[c, j + joff]
ep = ss.eff_porosity_col[c, j]
ts = temp.t_soisno_col[c, j + joff]
watsat_val = ss.watsat_col[c, j]
sucsat_val = ss.sucsat_col[c, j]
bsw_val = ss.bsw_col[c, j]
smpso_val = CLM.pftcon.smpso[itype]
smpsc_val = CLM.pftcon.smpsc[itype]
rf = ss.rootfr_patch[p, j]

@printf("  c=%d, itype=%d, joff=%d\n", c, itype, joff)
@printf("  h2osoi_liqvol[%d,%d] = %.6f\n", c, j+joff, liqvol)
@printf("  eff_porosity[%d,%d] = %.6f\n", c, j, ep)
@printf("  t_soisno[%d,%d] = %.2f (TFRZ-2=%.2f)\n", c, j+joff, ts, CLM.TFRZ-2.0)
@printf("  watsat=%.6f, sucsat=%.2f, bsw=%.4f\n", watsat_val, sucsat_val, bsw_val)
@printf("  smpso=%.1f, smpsc=%.1f\n", smpso_val, smpsc_val)
@printf("  rootfr=%.6f\n", rf)

if liqvol > 0.0 && ts > CLM.TFRZ - 2.0
    s_node = max(liqvol / ep, 0.01)
    smp_node = -sucsat_val * s_node^(-bsw_val)
    smp_node_clamped = max(smpsc_val, smp_node)
    rresis_val = min((ep / watsat_val) * (smp_node_clamped - smpsc_val) / (smpso_val - smpsc_val), 1.0)
    rootr_val = rf * rresis_val
    @printf("  s_node=%.6f, smp_node=%.2f, smp_clamped=%.2f\n", s_node, smp_node, smp_node_clamped)
    @printf("  rresis=%.6f, rootr=%.8f\n", rresis_val, rootr_val)
else
    println("  SKIPPED: liqvol<=0 or frozen")
end

# Also check what the actual rresis_patch looks like
println("\nrresis (patch 2, first 5 layers):")
for jj in 1:5
    @printf("  Layer %d: rresis=%.8f\n", jj, inst.energyflux.rresis_patch[2, jj])
end

println("\n=== FILTER STATE AFTER 48 STEPS ===")
cs = inst.canopystate
for p in 1:np
    @printf("Patch %d: exposedvegp=%s, noexposedvegp=%s, nolakeurbanp=%s, frac_veg_nosno=%d, frac_veg_nosno_alb=%d\n",
            p, filt.exposedvegp[p], filt.noexposedvegp[p], filt.nolakeurbanp[p],
            cs.frac_veg_nosno_patch[p], cs.frac_veg_nosno_alb_patch[p])
end

println("\nDone.")
