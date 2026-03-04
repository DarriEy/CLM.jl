using NCDatasets, Dates, Printf, Statistics
using CLM

# Diagnostic: run a few summer timesteps and print radiation intermediates
const basedir = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped"
fsurdat  = joinpath(basedir, "settings/CLM/parameters/surfdata_clm.nc")
paramfile = joinpath(basedir, "settings/CLM/parameters/clm5_params.nc")
fforcing = joinpath(basedir, "data/forcing/CLM_input/clmforc.2002.nc")

# Initialize for July 1
(inst, bounds, filt, tm) = CLM.clm_initialize!(fsurdat=fsurdat, paramfile=paramfile,
    start_date=DateTime(2002, 7, 1), dtime=1800)

nc = bounds.endc
np = bounds.endp
ng = bounds.endg

# Print PFT info
println("=== PATCH INFO ===")
for p in 1:np
    ptype = inst.patch.itype[p]
    @printf("Patch %d: itype=%d, wtgcell=%.4f, landunit=%d\n",
            p, ptype, inst.patch.wtgcell[p], inst.patch.landunit[p])
    @printf("  elai=%.4f, esai=%.4f, tlai=%.4f, tsai=%.4f\n",
            inst.canopystate.elai_patch[p], inst.canopystate.esai_patch[p],
            inst.canopystate.tlai_patch[p], inst.canopystate.tsai_patch[p])
end

println("\n=== LANDUNIT INFO ===")
for l in 1:bounds.endl
    @printf("Landunit %d: itype=%d\n", l, inst.landunit.itype[l])
end

# Print PFT optical properties
println("\n=== PFT OPTICAL PROPERTIES ===")
for ptype in [0, 1, 2, 12, 13]
    idx = ptype + 1  # Julia 1-based for what SHOULD be accessed
    idx0 = ptype      # What is ACTUALLY accessed (0-based itype)
    if idx0 >= 1 && idx0 <= size(CLM.pftcon.rhol, 1)
        @printf("PFT %d (accessed as idx=%d): rhol_vis=%.4f, rhol_nir=%.4f, taul_vis=%.4f, taul_nir=%.4f\n",
                ptype, idx0,
                CLM.pftcon.rhol[idx0, 1], CLM.pftcon.rhol[idx0, 2],
                CLM.pftcon.taul[idx0, 1], CLM.pftcon.taul[idx0, 2])
    end
    if idx >= 1 && idx <= size(CLM.pftcon.rhol, 1)
        @printf("PFT %d (correct  idx=%d): rhol_vis=%.4f, rhol_nir=%.4f, taul_vis=%.4f, taul_nir=%.4f\n",
                ptype, idx,
                CLM.pftcon.rhol[idx, 1], CLM.pftcon.rhol[idx, 2],
                CLM.pftcon.taul[idx, 1], CLM.pftcon.taul[idx, 2])
    end
end

# Print xl values
println("\n=== XL (leaf orientation) ===")
for ptype in [0, 1, 2, 12, 13]
    idx0 = ptype  # what is actually accessed
    idx = ptype + 1  # what should be accessed
    if idx0 >= 1 && idx0 <= length(CLM.pftcon.xl) && idx <= length(CLM.pftcon.xl)
        @printf("PFT %d: xl[%d]=%.4f (accessed), xl[%d]=%.4f (correct)\n",
                ptype, idx0, CLM.pftcon.xl[idx0], idx, CLM.pftcon.xl[idx])
    end
end

# Check ISTSOIL constant
println("\n=== CONSTANTS ===")
@printf("ISTSOIL=%d, ISTCROP=%d\n", CLM.ISTSOIL, CLM.ISTCROP)
@printf("NUMRAD=%d, NLEVCAN=%d\n", CLM.NUMRAD, CLM.NLEVCAN)
println("\n=== GRIDCELL LAT/LON ===")
for g in 1:ng
    @printf("Grid %d: lat=%.6f, lon=%.6f\n", g, inst.gridcell.lat[g], inst.gridcell.lon[g])
end
# Also check if there's latdeg/londeg
if hasproperty(inst.gridcell, :latdeg)
    @printf("  latdeg=%s, londeg=%s\n", string(inst.gridcell.latdeg), string(inst.gridcell.londeg))
end

# Set initial temps to summer values
inst.temperature.t_grnd_col .= 280.0
inst.temperature.t_soisno_col .= 280.0
inst.temperature.t_veg_patch .= 280.0

# Open forcing file
fr = CLM.ForcingReader()
CLM.forcing_reader_init!(fr, fforcing)
dtime = 1800

# Run 24 half-hourly steps (12 hours from 6AM to 6PM July 1)
# July 1 starts at hour 4344 = 181*24
step_start = 181 * 24 + 1  # hour index in forcing file

println("\n=== RUNNING 3 STEPS ===")
for step in 1:3
    CLM.advance_timestep!(tm)
    calday = CLM.get_curr_calday(tm)
    (yr, mon, d, tod) = CLM.get_curr_date(tm)

    # Read forcing
    CLM.read_forcing_step!(fr, inst.atm2lnd, tm.current_date, ng, nc)
    CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)

    # Compute orbital params
    (declin, eccf) = CLM.compute_orbital(calday)
    declinp1 = declin
    obliqr = CLM.ORB_OBLIQR_DEFAULT
    nextsw_cday = calday + Float64(dtime) / CLM.SECSPDAY

    @printf("  calday=%.4f, nextsw_cday=%.4f, declinp1=%.4f rad (%.2f deg)\n",
            calday, nextsw_cday, declinp1, rad2deg(declinp1))
    # Manual coszen check
    lon_g = inst.gridcell.lon[1]; lat_g = inst.gridcell.lat[1]
    ha = 2.0 * π * mod(nextsw_cday, 1.0) + lon_g * π / 180.0 - π
    cz_manual = max(sin(lat_g * π / 180.0) * sin(declinp1) + cos(lat_g * π / 180.0) * cos(declinp1) * cos(ha), 0.0)
    @printf("  Manual coszen check: lon=%.4f lat=%.4f hour_angle=%.4f coszen=%.4f\n",
            lon_g, lat_g, ha, cz_manual)

    # Print forcing and pre-albedo state
    sol_ad = inst.atm2lnd.forc_solad_downscaled_col
    sol_ai = inst.atm2lnd.forc_solai_grc
    fsds = sum(sol_ad[1, :]) + sum(sol_ai[1, :])

    # Check surfalb state BEFORE driver call (from previous step's albedo calc)
    alb = inst.surfalb
    fabd_p1 = [alb.fabd_patch[p, 1] for p in 1:np]
    fabd_p2 = [alb.fabd_patch[p, 2] for p in 1:np]
    ftdd_p1 = [alb.ftdd_patch[p, 1] for p in 1:np]
    albgrd_c = [alb.albgrd_col[c, 1] for c in 1:nc]

    # Run full driver step
    filt_ia = CLM.clump_filter_inactive_and_active
    config = CLM.CLMDriverConfig()

    CLM.clm_drv!(config, inst, filt, filt_ia, bounds,
                 true, nextsw_cday, declinp1, declin, obliqr,
                 false, false, "", false;
                 nstep=tm.nstep,
                 is_first_step=(tm.nstep == 1),
                 is_beg_curr_day=CLM.is_beg_curr_day(tm),
                 is_beg_curr_year=CLM.is_beg_curr_year(tm),
                 dtime=Float64(dtime),
                 mon=mon,
                 day=d,
                 photosyns=inst.photosyns)

    # Print FSA result
    fsa_patches = [inst.solarabs.fsa_patch[p] for p in 1:np]
    sabv = [inst.solarabs.sabv_patch[p] for p in 1:np]
    sabg = [inst.solarabs.sabg_patch[p] for p in 1:np]

    # Post-driver albedo (computed for NEXT step)
    fabd_post = [alb.fabd_patch[p, 1] for p in 1:np]
    ftdd_post = [alb.ftdd_patch[p, 1] for p in 1:np]
    fabi_post = [alb.fabi_patch[p, 1] for p in 1:np]
    albgrd_post = [alb.albgrd_col[c, 1] for c in 1:nc]
    albgri_post = [alb.albgri_col[c, 1] for c in 1:nc]
    coszen_grc = alb.coszen_grc[1]

    @printf("Step %2d (tod=%5d): FSDS=%.1f | PRE: fabd_VIS=%s ftdd_VIS=%s albgrd_VIS=%s | FSA=%s sabv=%s sabg=%s | POST: coszen=%.4f fabd_VIS=%s ftdd_VIS=%s albgrd_VIS=%s\n",
            step, tod, fsds,
            join([@sprintf("%.3f", x) for x in fabd_p1], ","),
            join([@sprintf("%.3f", x) for x in ftdd_p1], ","),
            join([@sprintf("%.3f", x) for x in albgrd_c], ","),
            join([@sprintf("%.1f", x) for x in fsa_patches], ","),
            join([@sprintf("%.1f", x) for x in sabv], ","),
            join([@sprintf("%.1f", x) for x in sabg], ","),
            coszen_grc,
            join([@sprintf("%.3f", x) for x in fabd_post], ","),
            join([@sprintf("%.3f", x) for x in ftdd_post], ","),
            join([@sprintf("%.3f", x) for x in albgrd_post], ","))
end

println("\nDone.")
