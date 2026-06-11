# clmdrv_make_data.jl — cloned from test/test_clm_driver.jl make_driver_data.
# Builds a minimal Float64 driver state and sets the shared module globals
# (varpar / varcon / pftcon / SNOW_DZ* / urban_ctrl). Kept in sync with the test.

function make_driver_data(; ng=2, nl=3, nc=4, np=6)
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
    CLM.varcon_init!()

    nlevsno = CLM.varpar.nlevsno
    nlevgrnd = CLM.varpar.nlevgrnd
    nlevtot = nlevsno + nlevgrnd
    nlevcan = CLM.NLEVCAN

    inst = CLM.CLMInstances()
    CLM.clm_instInit!(inst; ng=ng, nl=nl, nc=nc, np=np,
                       nlevdecomp_full=CLM.varpar.nlevdecomp_full)

    for c in 1:nc
        inst.column.landunit[c] = 1
        inst.column.gridcell[c] = 1
        inst.column.snl[c] = 0
        inst.column.patchi[c] = 1
        inst.column.patchf[c] = 0
        inst.column.nbedrock[c] = CLM.varpar.nlevsoi
    end
    for l in 1:nl
        inst.landunit.itype[l] = CLM.ISTSOIL
        inst.landunit.urbpoi[l] = false
        inst.landunit.lakpoi[l] = false
        inst.landunit.active[l] = true
    end
    for p in 1:np
        inst.patch.active[p] = true
        inst.patch.landunit[p] = 1
        inst.patch.gridcell[p] = 1
        inst.patch.column[p] = min(p, nc)
        inst.patch.itype[p] = 1
    end

    ppercol = max(1, np ÷ nc)
    pidx = 1
    for c in 1:nc
        inst.column.patchi[c] = pidx
        inst.column.patchf[c] = min(pidx + ppercol - 1, np)
        npc = inst.column.patchf[c] - inst.column.patchi[c] + 1
        for p in inst.column.patchi[c]:inst.column.patchf[c]
            inst.patch.column[p] = c
            inst.patch.wtcol[p] = 1.0 / npc
        end
        pidx = inst.column.patchf[c] + 1
    end
    for p in pidx:np
        inst.patch.column[p] = nc
        inst.patch.wtcol[p] = 0.0
    end

    bounds = CLM.BoundsType(
        begg=1, endg=ng, begl=1, endl=nl, begc=1, endc=nc, begp=1, endp=np,
        begCohort=0, endCohort=0,
        level=CLM.BOUNDS_LEVEL_CLUMP, clump_index=1)

    filt = CLM.ClumpFilter()
    CLM.alloc_filters!(filt, nc, np, nl)
    filt.allc .= true
    filt.nolakec .= true
    filt.nolakep .= true
    filt.soilc .= true
    filt.soilp .= true
    filt.nolakeurbanp .= true
    filt.nourbanp .= true
    filt.nourbanc .= true
    filt.hydrologyc .= true
    filt.urbanc .= false
    filt.urbanl .= false
    filt.urbanp .= false
    filt.snowc .= false
    filt.nosnowc .= true
    filt.do_smb_c .= false
    filt.exposedvegp .= true
    filt.noexposedvegp .= false
    filt.lakec .= false
    filt.lakep .= false
    filt.lakesnowc .= false
    filt.lakenosnowc .= false
    filt.bgc_soilc .= false
    filt.bgc_vegp .= false
    filt.pcropp .= false
    filt.soilnopcropp .= true
    filt.actfirec .= false
    filt.actfirep .= false

    filt_ia = CLM.ClumpFilter()
    CLM.alloc_filters!(filt_ia, nc, np, nl)
    filt_ia.allc .= true
    filt_ia.nolakec .= true
    filt_ia.nolakep .= true
    filt_ia.nourbanc .= true
    filt_ia.nourbanp .= true

    CLM.urban_read_nml!(CLM.urban_ctrl)

    # Initialize the global photosynthesis params (the suite gets this as a
    # side effect of test_photosynthesis.jl running before the driver test).
    CLM.photo_params_init!(CLM.params_inst)

    dzmin = zeros(nlevsno)
    dzmax_u = zeros(nlevsno)
    dzmax_l = zeros(nlevsno)
    dzmin[1] = 0.010; dzmax_u[1] = 0.02; dzmax_l[1] = 0.03
    dzmin[2] = 0.015; dzmax_u[2] = 0.05; dzmax_l[2] = 0.07
    for j in 3:nlevsno
        dzmin[j] = dzmax_u[j-1] * 0.5
        dzmax_u[j] = 2.0 * dzmax_u[j-1] + 0.01
        dzmax_l[j] = dzmax_u[j] + dzmax_l[j-1]
        if j == nlevsno
            dzmax_u[j] = floatmax(Float64)
            dzmax_l[j] = floatmax(Float64)
        end
    end
    CLM.SNOW_DZMIN[] = dzmin
    CLM.SNOW_DZMAX_U[] = dzmax_u
    CLM.SNOW_DZMAX_L[] = dzmax_l

    config = CLM.CLMDriverConfig()

    npft = CLM.MXPFT + 1
    p = CLM.pftcon
    p.dleaf         = fill(0.04, npft)
    p.slatop        = fill(0.01, npft)
    p.leafcn        = fill(25.0, npft)
    p.flnr          = fill(0.1, npft)
    p.fnitr         = fill(0.1, npft)
    p.mbbopt        = fill(9.0, npft)
    p.c3psn         = fill(1.0, npft)
    p.woody         = fill(0.0, npft)
    p.smpso         = fill(-66000.0, npft)
    p.smpsc         = fill(-275000.0, npft)
    p.z0mr          = fill(0.055, npft)
    p.displar       = fill(0.67, npft)
    p.xl            = fill(0.1, npft)
    p.rhol          = fill(0.1, npft, 2)
    p.rhos          = fill(0.2, npft, 2)
    p.taul          = fill(0.05, npft, 2)
    p.taus          = fill(0.1, npft, 2)
    p.medlynintercept = fill(100.0, npft)
    p.medlynslope     = fill(6.0, npft)
    p.crop          = fill(0.0, npft)

    photosyns = CLM.PhotosynthesisData()
    CLM.photosynthesis_data_init!(photosyns, np)

    return inst, bounds, filt, filt_ia, config, photosyns
end

# --------------------------------------------------------------------------
# make_driver_data_physical — same topology/filters/pftcon/globals as
# make_driver_data(), but ADDITIONALLY seeds physically realistic cold-start
# initial conditions (atmospheric forcing, column geometry, soil texture +
# thermal/hydraulic properties, soil/snow temperature, soil water, canopy
# state) so that running clm_drv! on the CPU produces FINITE values across the
# bulk of the state instead of the ~112 finite entries the smoke fixture gives.
#
# Snow-free cold start (snl=0, no snowpack), above-freezing soil (~283 K),
# moderate moisture — chosen to keep every physics branch finite and active.
#
# make_driver_data() is intentionally NOT modified (it mirrors the test).
# --------------------------------------------------------------------------
function make_driver_data_physical(; ng=2, nl=3, nc=4, np=6)
    inst, bounds, filt, filt_ia, config, photosyns =
        make_driver_data(; ng=ng, nl=nl, nc=nc, np=np)

    vp        = CLM.varpar
    nlevsno   = vp.nlevsno
    nlevsoi   = vp.nlevsoi
    nlevgrnd  = vp.nlevgrnd
    nlevmax   = vp.nlevmaxurbgrnd
    joff      = nlevsno                       # snow-layer offset into t_soisno/dz/z/zi
    TFRZ      = CLM.TFRZ                       # 273.15 K

    # ---- 0. Gridcell geometry (lat/lon) -----------------------------------
    # grc.lat/lon default to NaN; the driver's coszen closure
    #   cosz = sin(lat)sin(decl) - cos(lat)cos(decl)cos(jday*2π + lon)
    # would then yield NaN, sending every column down the coszen<=0 ("night")
    # branch (fabd=fabi=0, the whole solar-radiation -> energy chain inactive).
    # Pick lat=0, lon=π so that with the harness DRV_ARGS (nextsw_cday=1.0,
    # declin=0): cosz = -cos(2π + π) = -cos(3π) = +1 → sun overhead, daytime.
    grc = inst.gridcell
    for g in 1:ng
        grc.lat[g]    = 0.0
        grc.lon[g]    = Float64(pi)
        grc.latdeg[g] = 0.0
        grc.londeg[g] = 180.0
    end

    # ---- 1. Atmospheric forcing (gridcell + downscaled-to-column) ----------
    a2l = inst.atm2lnd
    Tair = 283.0                              # 10 C air temperature
    pbot = 1.0e5                              # surface pressure (Pa)
    qair = 0.006                              # specific humidity (kg/kg)
    # vapor pressure from q: e = q*p / (0.622 + 0.378*q)
    eair = qair * pbot / (0.622 + 0.378 * qair)
    rho  = pbot / (287.04 * Tair * (1 + 0.61 * qair))   # moist air density
    for g in 1:ng
        a2l.forc_u_grc[g]    = 2.0
        a2l.forc_v_grc[g]    = 1.0
        a2l.forc_wind_grc[g] = sqrt(2.0^2 + 1.0^2)
        a2l.forc_hgt_grc[g]   = 30.0
        a2l.forc_topo_grc[g]  = 0.0
        a2l.forc_hgt_u_grc[g] = 30.0
        a2l.forc_hgt_t_grc[g] = 30.0
        a2l.forc_hgt_q_grc[g] = 30.0
        a2l.forc_vp_grc[g]    = eair
        a2l.forc_pco2_grc[g]  = 35.0          # CO2 partial pressure (Pa) ~ 350 ppm
        a2l.forc_po2_grc[g]   = 0.209 * pbot  # O2 partial pressure (Pa)
        a2l.forc_pch4_grc[g]  = 0.17          # CH4 partial pressure (Pa)
        a2l.forc_solar_not_downscaled_grc[g] = 400.0
        for b in 1:CLM.NUMRAD
            a2l.forc_solad_not_downscaled_grc[g, b] = 150.0   # direct beam (W/m2)
            a2l.forc_solai_grc[g, b]                = 50.0    # diffuse (W/m2)
        end
        a2l.forc_t_not_downscaled_grc[g]    = Tair
        a2l.forc_pbot_not_downscaled_grc[g] = pbot
        a2l.forc_th_not_downscaled_grc[g]   = Tair            # ~potential T at z=0
        a2l.forc_rho_not_downscaled_grc[g]  = rho
        a2l.forc_lwrad_not_downscaled_grc[g]= 300.0
        a2l.forc_rain_not_downscaled_grc[g] = 0.0
        a2l.forc_snow_not_downscaled_grc[g] = 0.0
        a2l.forc_q_not_downscaled_grc[g]    = qair
    end
    for c in 1:nc
        a2l.forc_t_downscaled_col[c]     = Tair
        a2l.forc_pbot_downscaled_col[c]  = pbot
        a2l.forc_th_downscaled_col[c]    = Tair
        a2l.forc_rho_downscaled_col[c]   = rho
        a2l.forc_lwrad_downscaled_col[c] = 300.0
        a2l.forc_solar_downscaled_col[c] = 400.0
        for b in 1:CLM.NUMRAD
            a2l.forc_solad_downscaled_col[c, b] = 150.0
        end
        a2l.forc_rain_downscaled_col[c]  = 0.0
        a2l.forc_snow_downscaled_col[c]  = 0.0
        a2l.forc_q_downscaled_col[c]     = qair
    end

    # ---- 2. Column geometry (dz/z/zi from varcon zsoi/dzsoi/zisoi) ---------
    # Snow levels (above joff) are unused with snl=0; zero their thicknesses.
    zsoi_g  = CLM.zsoi[]
    dzsoi_g = CLM.dzsoi[]
    zisoi_g = CLM.zisoi[]
    for c in 1:nc
        CLM._set_standard_soil!(inst.column, c, joff, nlevgrnd, nlevmax,
                                zsoi_g, dzsoi_g, zisoi_g)
        for j in 1:nlevsno
            inst.column.dz[c, j] = 0.0
            inst.column.z[c, j]  = 0.0
            inst.column.zi[c, j] = 0.0
        end
        inst.column.zi[c, nlevsno + 1] = 0.0   # snow/soil interface at surface
    end

    # ---- 3. Soil texture + hydraulic + thermal properties -----------------
    # soilstate_init! allocates the matrices the driver path reads; (re)allocate
    # any left at size (0,0) so the per-(c,j) fill below is safe.
    ss = inst.soilstate
    nsno_g = nlevsno + nlevmax
    if size(ss.watsat_col, 2)       < nlevmax;  ss.watsat_col       = fill(0.0, nc, nlevmax); end
    if size(ss.bsw_col, 2)          < nlevgrnd; ss.bsw_col          = fill(0.0, nc, nlevgrnd); end
    if size(ss.sucsat_col, 2)       < nlevgrnd; ss.sucsat_col       = fill(0.0, nc, nlevgrnd); end
    if size(ss.hksat_col, 2)        < nlevgrnd; ss.hksat_col        = fill(0.0, nc, nlevgrnd); end
    if size(ss.hksat_min_col, 2)    < nlevgrnd; ss.hksat_min_col    = fill(0.0, nc, nlevgrnd); end
    if size(ss.eff_porosity_col, 2) < nlevgrnd; ss.eff_porosity_col = fill(0.0, nc, nlevgrnd); end
    if size(ss.watfc_col, 2)        < nlevgrnd; ss.watfc_col        = fill(0.0, nc, nlevgrnd); end
    if size(ss.watdry_col, 2)       < nlevgrnd; ss.watdry_col       = fill(0.0, nc, nlevgrnd); end
    if size(ss.watopt_col, 2)       < nlevgrnd; ss.watopt_col       = fill(0.0, nc, nlevgrnd); end
    if size(ss.smp_l_col, 2)        < nlevgrnd; ss.smp_l_col        = fill(0.0, nc, nlevgrnd); end
    if size(ss.hk_l_col, 2)         < nlevgrnd; ss.hk_l_col         = fill(0.0, nc, nlevgrnd); end
    if size(ss.tkmg_col, 2)         < nlevgrnd; ss.tkmg_col         = fill(0.0, nc, nlevgrnd); end
    if size(ss.tkdry_col, 2)        < nlevgrnd; ss.tkdry_col        = fill(0.0, nc, nlevgrnd); end
    if size(ss.tksatu_col, 2)       < nlevgrnd; ss.tksatu_col       = fill(0.0, nc, nlevgrnd); end
    if size(ss.csol_col, 2)         < nlevgrnd; ss.csol_col         = fill(0.0, nc, nlevgrnd); end
    if size(ss.thk_col, 2)          < nsno_g;   ss.thk_col          = fill(0.0, nc, nsno_g);   end
    if size(ss.rootfr_col, 2)       < nlevgrnd; ss.rootfr_col       = fill(0.0, nc, nlevgrnd); end
    if size(ss.rootfr_patch, 2)     < nlevgrnd; ss.rootfr_patch     = fill(0.0, np, nlevgrnd); end
    if size(ss.rootr_patch, 2)      < nlevgrnd; ss.rootr_patch      = fill(0.0, np, nlevgrnd); end

    for c in 1:nc
        ss.smpmin_col[c] = -1.0e8
        for j in 1:nlevgrnd
            ss.watsat_col[c, j]       = 0.40
            ss.bsw_col[c, j]          = 5.0
            ss.sucsat_col[c, j]       = 100.0
            ss.hksat_col[c, j]        = 1.0e-4
            ss.hksat_min_col[c, j]    = 1.0e-4
            ss.eff_porosity_col[c, j] = 0.40
            ss.watfc_col[c, j]        = 0.25
            ss.watdry_col[c, j]       = 0.05
            ss.watopt_col[c, j]       = 0.35
            ss.smp_l_col[c, j]        = 0.0
            ss.hk_l_col[c, j]         = 0.0
            ss.tkmg_col[c, j]         = 2.0
            ss.tkdry_col[c, j]        = 0.25
            ss.tksatu_col[c, j]       = 1.5
            ss.csol_col[c, j]         = 2.0e6
        end
        for j in 1:nsno_g
            ss.thk_col[c, j] = 0.0
        end
    end
    # Root fractions (sum to 1 over the soil column) — uniform over nlevsoi.
    for c in 1:nc, j in 1:nlevgrnd
        ss.rootfr_col[c, j] = j <= nlevsoi ? 1.0 / nlevsoi : 0.0
    end
    for p in 1:np, j in 1:nlevgrnd
        ss.rootfr_patch[p, j] = j <= nlevsoi ? 1.0 / nlevsoi : 0.0
        ss.rootr_patch[p, j]  = j <= nlevsoi ? 1.0 / nlevsoi : 0.0
    end

    # ---- 4. Temperature state ---------------------------------------------
    temp = inst.temperature
    for c in 1:nc
        temp.t_grnd_col[c]   = 283.0
        temp.t_h2osfc_col[c] = 283.0
        temp.emg_col[c]      = 0.97
        temp.thv_col[c]      = 283.0 * (1 + 0.61 * qair)
        for j in 1:nlevsno
            temp.t_soisno_col[c, j] = TFRZ          # unused snow layers, keep finite
        end
        for j in 1:nlevmax
            temp.t_soisno_col[c, j + joff] = 283.0
        end
    end
    for p in 1:np
        temp.t_veg_patch[p] = 283.0
    end

    # ---- 5. Water state ----------------------------------------------------
    ws = inst.water.waterstatebulk_inst.ws
    DENH2O = CLM.DENH2O
    nlevtot = nlevsno + nlevmax
    for c in 1:nc
        ws.h2osfc_col[c]            = 0.0
        ws.h2osno_no_layers_col[c] = 0.0
        for j in 1:nlevtot
            ws.h2osoi_liq_col[c, j] = 0.0
            ws.h2osoi_ice_col[c, j] = 0.0
            ws.excess_ice_col[c, j] = 0.0
        end
        for j in 1:nlevmax
            ws.h2osoi_vol_col[c, j] = 0.0
        end
        # Soil layers (below joff): unfrozen liquid water at ~50% saturation.
        for j in 1:nlevgrnd
            jj  = j + joff
            vol = 0.20                       # vol. soil moisture (< watsat=0.40)
            ws.h2osoi_liq_col[c, jj] = vol * inst.column.dz[c, jj] * DENH2O
            ws.h2osoi_ice_col[c, jj] = 0.0
            ws.h2osoi_vol_col[c, j]  = vol
        end
    end

    # ---- 6b. Soil-color albedo constants (isoicol + albsat/albdry lookup) --
    # clm_initialize! normally fills inst.surfalb_con from surfdata.soil_color;
    # the fixture skips initialization, so the daytime soil-albedo kernel
    # (surface_albedo.jl:620) would index an empty isoicol. Populate it the same
    # way with a uniform mid soil-color class (valid in 1:20).
    soic2d = fill(5, ng)
    CLM.surface_albedo_init_time_const!(inst.surfalb_con, 20, soic2d,
                                        inst.column.gridcell[1:nc], 1:nc, 1:ng)

    # ---- 6. Canopy state ---------------------------------------------------
    cs = inst.canopystate
    for p in 1:np
        cs.tlai_patch[p] = 2.0
        cs.tsai_patch[p] = 0.5
        cs.elai_patch[p] = 2.0
        cs.esai_patch[p] = 0.5
        cs.htop_patch[p] = 1.0
        cs.hbot_patch[p] = 0.1
        cs.frac_veg_nosno_patch[p]     = 1
        cs.frac_veg_nosno_alb_patch[p] = 1
    end

    return inst, bounds, filt, filt_ia, config, photosyns
end
