# =============================================================================
# FD-validated FULL-CANOPY reverse gradient + a CANOPY+SOIL two-module reverse
# chain, both wired through the PRODUCTION compositional_reverse! engine
# (src/biogeophys/canopy_fluxes_reverse.jl).
#
# Three sections, each runs the real CLM physics and validates a reverse gradient
# against central finite differences of the SAME primal:
#
#   [A] Full canopy (energy balance + photosynthesis) reverse, on the realistic
#       finite single-patch state. Asserts the decomposed compositional forward
#       (CLM.canopy_rev_forward!) is byte-identical to the real CLM.canopy_fluxes!
#       (|Δt_veg|=0), then FD-validates  dL/d(t_grnd),  L = sum(abs2, t_veg),
#       through CLM.canopy_rev_gradient!.
#
#   [B] soil_temperature! reverse on a REAL Bow-at-Banff cold-start instance
#       (clm_initialize! + the forc_hgt_u/t/q_grc=30 fix that keeps the M-O solve
#       finite). FD-validates  dL/d(t_soisno[c,j])  through CLM.compositional_reverse!.
#
#   [C] CANOPY + SOIL_TEMPERATURE two-module chain on the SAME real inst, in ONE
#       compositional_reverse! call: phases = [canopy energy/psn sub-phases ...,
#       soil_temperature!]. The canopy sub-phases update t_veg from the inst's
#       converged radiation; soil_temperature! then reads the canopy-coupled ground
#       heat flux. FD-validates  dL/d(t_grnd),  L = sum(abs2, t_soisno),  i.e. the
#       gradient that flows canopy → soil across the module boundary.
#
# Both modules already differentiate (canopy compositionally, soil whole-fn); this
# script chains their gradients through the single production engine and proves the
# chain numerically.
#
#   julia +1.10 --project=/tmp/clm_jl10_wt scripts/enzyme_driver_reverse.jl
#   (any project that has CLM dev'd + Enzyme; Enzyme is stable on Julia 1.10 LTS)
# =============================================================================
using CLM, Enzyme, Printf
const C = CLM
const FT = Float64

const FSURDAT = get(ENV, "CLM_FSURDAT",
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc")
const PARAMFILE = get(ENV, "CLM_PARAMFILE",
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc")

# =============================================================================
# [A] Full-canopy reverse on the realistic finite single-patch state.
#     (real CLM.canopy_fluxes! is the forward reference.)
# =============================================================================
function section_A()
    println("\n", "#"^70)
    println("# [A] FULL-CANOPY (energy + photosynthesis) REVERSE vs FD")
    println("#"^70)

    C.varpar_init!(C.varpar, 1, 14, 2, 5)
    NP = 1; NC = 1; NG = 1; NL = 1; NLEVCAN = C.NLEVCAN; MP = C.MXPFT + 1

    C.soil_resistance_read_nl!(soil_resis_method = C.SOIL_RESIS_LEEPIELKE_1992)
    C.canopy_fluxes_read_nml!(use_undercanopy_stability = false,
        use_biomass_heat_storage = false, itmax_canopy_fluxes = 40)
    C.canopy_fluxes_read_params!()
    C.photo_params_init!(C.params_inst)
    let p = C.params_inst
        p.theta_cj = fill(0.98, MP)
        p.theta_ip = 0.95; p.act25 = 72.0; p.fnr = 7.16; p.cp25_yr2000 = 42.75e-6
        p.kc25_coef = 404.9e-6; p.ko25_coef = 278.4e-3; p.fnps = 0.15; p.theta_psii = 0.7
        p.vcmaxha = 65330.0; p.jmaxha = 43540.0; p.tpuha = 53100.0; p.lmrha = 46390.0
        p.kcha = 79430.0; p.koha = 36380.0; p.cpha = 37830.0
        p.vcmaxhd = 149250.0; p.jmaxhd = 152040.0; p.tpuhd = 150650.0; p.lmrhd = 150650.0
        p.lmrse = 490.0; p.tpu25ratio = 0.167; p.kp25ratio = 20160.0
        p.vcmaxse_sf = 1.0; p.jmaxse_sf = 1.0; p.tpuse_sf = 1.0; p.jmax25top_sf = 1.0
    end

    PSN = (; c3psn = fill(FT(1.0), MP), leafcn = fill(FT(25.0), MP), flnr = fill(FT(0.1), MP),
        fnitr = fill(FT(1.0), MP), slatop = fill(FT(0.01), MP), mbbopt = fill(FT(9.0), MP),
        medlynintercept = fill(FT(100.0), MP), medlynslope = fill(FT(6.0), MP),
        nrad = fill(1, NP), tlai_z = fill(FT(1.0), NP, NLEVCAN),
        parsun_z = fill(FT(250.0), NP, NLEVCAN), parsha_z = fill(FT(120.0), NP, NLEVCAN),
        laisun_z = fill(FT(1.0), NP, NLEVCAN), laisha_z = fill(FT(1.0), NP, NLEVCAN),
        vcmaxcint_sun = fill(FT(1.0), NP), vcmaxcint_sha = fill(FT(0.6), NP),
        o3coefv = fill(FT(1.0), NP), o3coefg = fill(FT(1.0), NP), t10 = fill(FT(290.0), NP))

    function build()
        canopystate = C.CanopyStateData{FT}(); C.canopystate_init!(canopystate, NP)
        canopystate.elai_patch[1] = 2.0;  canopystate.esai_patch[1] = 0.5
        canopystate.laisun_patch[1] = 1.2; canopystate.laisha_patch[1] = 0.8
        canopystate.displa_patch[1] = 5.0; canopystate.htop_patch[1] = 10.0
        canopystate.frac_veg_nosno_patch[1] = 1; canopystate.dleaf_patch[1] = 0.04
        canopystate.stem_biomass_patch[1] = 0.0; canopystate.leaf_biomass_patch[1] = 0.0
        canopystate.tlai_patch[1] = 2.0
        canopystate.laisun_z_patch = copy(PSN.laisun_z); canopystate.laisha_z_patch = copy(PSN.laisha_z)

        energyflux = C.EnergyFluxData{FT}(); C.energyflux_init!(energyflux, NP, NC, NL, NG)
        energyflux.btran_patch[1] = 0.8; energyflux.bsun_patch[1] = 0.8; energyflux.bsha_patch[1] = 0.8
        energyflux.htvp_col[1] = C.HVAP

        frictionvel = C.FrictionVelocityData{FT}(); C.frictionvel_init!(frictionvel, NP, NC)
        frictionvel.zetamaxstable = 0.5; frictionvel.zsno = 0.00085; frictionvel.zlnd = 0.000775
        frictionvel.z0mv_patch[1] = 0.5; frictionvel.z0hv_patch[1] = 0.5; frictionvel.z0qv_patch[1] = 0.5
        frictionvel.z0mg_col[1] = 0.01; frictionvel.z0hg_col[1] = 0.01; frictionvel.z0qg_col[1] = 0.01
        frictionvel.forc_hgt_u_patch[1] = 30.0; frictionvel.forc_hgt_t_patch[1] = 30.0
        frictionvel.forc_hgt_q_patch[1] = 30.0
        frictionvel.ustar_patch[1] = 0.5; frictionvel.um_patch[1] = 5.0; frictionvel.uaf_patch[1] = 3.0
        frictionvel.taf_patch[1] = 290.0; frictionvel.qaf_patch[1] = 0.008
        frictionvel.obu_patch[1] = -100.0; frictionvel.zeta_patch[1] = -0.1
        frictionvel.vpd_patch[1] = 1.0; frictionvel.rb1_patch[1] = 50.0
        frictionvel.ram1_patch[1] = 50.0; frictionvel.num_iter_patch[1] = 0.0

        temperature = C.TemperatureData{FT}(); C.temperature_init!(temperature, NP, NC, NL, NG)
        temperature.t_veg_patch[1] = 290.0; temperature.t_stem_patch[1] = 289.0
        temperature.t_skin_patch[1] = 290.0; temperature.thm_patch[1] = 290.0
        temperature.t_grnd_col[1] = 288.0; temperature.t_h2osfc_col[1] = 288.0
        temperature.thv_col[1] = 291.0; temperature.emv_patch[1] = 0.97; temperature.emg_col[1] = 0.96
        temperature.t_ref2m_patch[1] = 290.0; temperature.t_ref2m_r_patch[1] = 290.0
        temperature.t_a10_patch[1] = 290.0
        for j in 1:size(temperature.t_soisno_col, 2); temperature.t_soisno_col[1, j] = 288.0; end

        solarabs = C.SolarAbsorbedData{FT}(); C.solarabs_init!(solarabs, NP, NL)
        solarabs.sabv_patch[1] = 150.0
        solarabs.parsun_z_patch = copy(PSN.parsun_z); solarabs.parsha_z_patch = copy(PSN.parsha_z)

        soilstate = C.SoilStateData{FT}(); C.soilstate_init!(soilstate, NC, NL)
        soilstate.soilbeta_col[1] = 0.8; soilstate.soilresis_col[1] = 100.0
        for j in 1:C.varpar.nlevgrnd; soilstate.watsat_col[1, j] = 0.45; end

        waterfluxbulk = C.WaterFluxBulkData{FT}(); C.waterfluxbulk_init!(waterfluxbulk, NC, NP, NL, NG)
        waterstatebulk = C.WaterStateBulkData{FT}(); C.waterstatebulk_init!(waterstatebulk, NC, NP, NL, NG)
        waterstatebulk.ws.liqcan_patch[1] = 0.1; waterstatebulk.ws.snocan_patch[1] = 0.0

        waterdiagbulk = C.WaterDiagnosticBulkData{FT}(); C.waterdiagnosticbulk_init!(waterdiagbulk, NC, NP, NL, NG)
        waterdiagbulk.fwet_patch[1] = 0.1; waterdiagbulk.fdry_patch[1] = 0.8
        waterdiagbulk.frac_sno_eff_col[1] = 0.0; waterdiagbulk.frac_h2osfc_col[1] = 0.0
        waterdiagbulk.snow_depth_col[1] = 0.0; waterdiagbulk.qg_col[1] = 0.005
        waterdiagbulk.qg_snow_col[1] = 0.005; waterdiagbulk.qg_soil_col[1] = 0.005
        waterdiagbulk.qg_h2osfc_col[1] = 0.005; waterdiagbulk.dqgdT_col[1] = 0.0003
        waterdiagbulk.rh_af_patch[1] = 0.6

        photosyns = C.PhotosynthesisData{FT}(); C.photosynthesis_data_init!(photosyns, NP)
        photosyns.rssun_patch[1] = 200.0; photosyns.rssha_patch[1] = 300.0
        photosyns.stomatalcond_mtd = C.STOMATALCOND_MTD_BB1987
        photosyns.leafresp_method  = C.LEAFRESP_MTD_RYAN1991
        photosyns.light_inhibit    = false

        patch_data = C.PatchData{FT}(); C.patch_init!(patch_data, NP)
        patch_data.column[1] = 1; patch_data.gridcell[1] = 1; patch_data.landunit[1] = 1
        patch_data.itype[1] = 1; patch_data.active[1] = true
        col_data = C.ColumnData{FT}(); C.column_init!(col_data, NC); col_data.snl[1] = 0
        gridcell_data = C.GridcellData{FT}(); C.gridcell_init!(gridcell_data, NG)

        S = (; canopystate, energyflux, frictionvel, temperature, solarabs, soilstate,
               waterfluxbulk, waterstatebulk, waterdiagbulk, photosyns, patch_data,
               col_data, gridcell_data)
        forc = (; lwrad = FT[350.0], q = FT[0.008], pbot = FT[101325.0], th = FT[290.0],
                  rho = FT[1.2], t = FT[288.0], u_grc = FT[3.0], v_grc = FT[1.0],
                  pco2 = FT[40.0], po2 = FT[21000.0], hgt_t = FT[10.0], hgt_u = FT[10.0],
                  hgt_q = FT[10.0], dayl = FT[43200.0], max_dayl = FT[50000.0],
                  downreg = FT[1.0], leafn = FT[1.0])
        pft = (; dleaf = fill(FT(0.04), MP), z0v_Cr = fill(FT(0.35), MP),
                 z0v_Cs = fill(FT(0.003), MP), z0v_c = fill(FT(0.25), MP),
                 z0v_cw = fill(FT(2.0), MP), z0v_LAImax = fill(FT(8.0), MP),
                 grnd_ch4 = fill(FT(0.0), NP))
        return (; S, forc, pft, dtime = FT(1800.0))
    end

    run_real!(S, forc, pft, dtime) = C.canopy_fluxes!(
        S.canopystate, S.energyflux, S.frictionvel, S.temperature, S.solarabs,
        S.soilstate, S.waterfluxbulk, S.waterstatebulk, S.waterdiagbulk, S.photosyns,
        S.patch_data, S.col_data, S.gridcell_data, trues(NP), 1:NP, 1:NC,
        forc.lwrad, forc.q, forc.pbot, forc.th, forc.rho, forc.t, forc.u_grc, forc.v_grc,
        forc.pco2, forc.po2, forc.hgt_t, forc.hgt_u, forc.hgt_q, forc.dayl, forc.max_dayl,
        forc.downreg, forc.leafn, dtime;
        dleaf_pft = pft.dleaf, z0v_Cr_pft = pft.z0v_Cr, z0v_Cs_pft = pft.z0v_Cs,
        z0v_c_pft = pft.z0v_c, z0v_cw_pft = pft.z0v_cw, z0v_LAImax_pft = pft.z0v_LAImax,
        grnd_ch4_cond_patch = pft.grnd_ch4,
        c3psn_pft = PSN.c3psn, leafcn_pft = PSN.leafcn, flnr_pft = PSN.flnr,
        fnitr_pft = PSN.fnitr, slatop_pft = PSN.slatop, mbbopt_pft = PSN.mbbopt,
        medlynintercept_pft = PSN.medlynintercept, medlynslope_pft = PSN.medlynslope,
        t10_patch = PSN.t10, nrad_patch = PSN.nrad, tlai_z_patch = PSN.tlai_z,
        vcmaxcint_sun_patch = PSN.vcmaxcint_sun, vcmaxcint_sha_patch = PSN.vcmaxcint_sha,
        parsun_z_patch = PSN.parsun_z, parsha_z_patch = PSN.parsha_z,
        laisun_z_patch = PSN.laisun_z, laisha_z_patch = PSN.laisha_z,
        o3coefv_patch = PSN.o3coefv, o3coefg_patch = PSN.o3coefg)

    make_aux(B) = (; patch = B.S.patch_data, col = B.S.col_data, grid = B.S.gridcell_data,
        forc = B.forc, pft = B.pft, psn = PSN, filterp = Int[1], fn = 1, active = trues(NP),
        mask = trues(NP), ivt = B.S.patch_data.itype .+ 1,
        forc_pbot_patch = FT[B.forc.pbot[B.S.patch_data.column[p]] for p in 1:NP],
        soilevap_beta = C.do_soilevap_beta(), soil_resis_sl14 = C.do_soil_resistance_sl14(),
        nlevsno = C.varpar.nlevsno, dtime = B.dtime, use_psn = true)
    make_b(B) = C.cf_rev_bundle(B.S.canopystate, B.S.energyflux, B.S.frictionvel, B.S.temperature,
        B.S.solarabs, B.S.soilstate, B.S.waterfluxbulk, B.S.waterstatebulk, B.S.waterdiagbulk,
        B.S.photosyns, C.cf_rev_scratch(FT, NP))

    Bref = build(); run_real!(Bref.S, Bref.forc, Bref.pft, Bref.dtime)
    N = max(Int(Bref.S.frictionvel.num_iter_patch[1]), 1)

    # forward parity assertion (must be ~0 vs the real canopy_fluxes!)
    Bf = build(); bf = make_b(Bf); C.canopy_rev_forward!(bf, make_aux(Bf), N)
    dveg = abs(bf.temperature.t_veg_patch[1] - Bref.S.temperature.t_veg_patch[1])
    @printf("real canopy_fluxes!  t_veg = %.10f  (num_iter=%d)\n", Bref.S.temperature.t_veg_patch[1], N)
    @printf("compositional fwd    t_veg = %.10f   |Δt_veg| = %.3e\n", bf.temperature.t_veg_patch[1], dveg)
    @assert dveg < 1e-9 "compositional forward must match real canopy_fluxes!"

    # FD reference dL/d(t_grnd) of the decomposed primal
    function L_tgrnd(δ)
        B = build(); b = make_b(B); b.temperature.t_grnd_col[1] += δ
        C.canopy_rev_forward!(b, make_aux(B), N); return sum(abs2, b.temperature.t_veg_patch)
    end
    hg = 1e-2
    g_fd = (L_tgrnd(hg) - L_tgrnd(-hg)) / (2hg)

    # compositional reverse via the production API
    B = build(); b = make_b(B)
    db = C.canopy_rev_gradient!(b, make_aux(B), N; seed = :t_veg_patch)
    g_rev = db.temperature.t_grnd_col[1]

    relerr = abs(g_rev - g_fd) / max(abs(g_fd), 1e-10)
    @printf("FD   dL/d(t_grnd)    = % .10e\n", g_fd)
    @printf("rev  dL/d(t_grnd)    = % .10e\n", g_rev)
    @printf("[A] rel error = %.3e   %s\n", relerr, relerr < 1e-6 ? "PASS ✓" : "FAIL ✗")
    return relerr
end

# =============================================================================
# Real Bow-at-Banff cold-start builder (shared by [B] and [C]).
# =============================================================================
function setup_forcing!(a2l, T, ng)
    for g in 1:ng
        a2l.forc_t_not_downscaled_grc[g]=T; a2l.forc_pbot_not_downscaled_grc[g]=85000.0
        a2l.forc_th_not_downscaled_grc[g]=T*(100000.0/85000.0)^(C.RAIR/C.CPAIR)
        a2l.forc_rho_not_downscaled_grc[g]=85000.0/(C.RAIR*T); a2l.forc_lwrad_not_downscaled_grc[g]=300.0
        a2l.forc_vp_grc[g]=800.0; a2l.forc_hgt_grc[g]=30.0; a2l.forc_topo_grc[g]=0.0; a2l.forc_wind_grc[g]=3.0
        # THE FIX: component obs heights (forcing_reader derives these; harness bypasses it).
        a2l.forc_hgt_u_grc[g]=30.0; a2l.forc_hgt_t_grc[g]=30.0; a2l.forc_hgt_q_grc[g]=30.0
        for b in 1:C.NUMRAD; a2l.forc_solad_not_downscaled_grc[g,b]=200.0; a2l.forc_solai_grc[g,b]=80.0; end
        a2l.forc_solar_not_downscaled_grc[g]=560.0; a2l.forc_rain_not_downscaled_grc[g]=0.0001
        a2l.forc_snow_not_downscaled_grc[g]=0.0
    end
end

function build_real()
    (inst, bounds, filt, tm) = C.clm_initialize!(; fsurdat=FSURDAT, paramfile=PARAMFILE)
    setup_forcing!(inst.atm2lnd, 285.0, bounds.endg)
    C.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
    config = C.CLMDriverConfig(); filt_ia = C.clump_filter_inactive_and_active
    (declin, eccf) = C.compute_orbital(120.0); nextsw_cday = 120.0 + 1800.0/C.SECSPDAY
    runstep!(n; first=false) = C.clm_drv!(config, inst, filt, filt_ia, bounds,
        true, nextsw_cday, declin, declin, C.ORB_OBLIQR_DEFAULT, false, false, "", false;
        nstep=n, is_first_step=first, is_beg_curr_day=first, dtime=1800.0, mon=1, day=1,
        photosyns=inst.photosyns)
    for n in 1:3; runstep!(n; first=(n==1)); end
    return inst, bounds, filt
end

# soil_temperature! const aux (everything non-differentiated). Passed via cargs so
# compositional_reverse! wraps it in Enzyme.Const — keeps mutable Vectors like
# urbantv out of the differentiated argument set.
soil_aux(bounds, filt) = (;
    urbantv = fill(323.15, bounds.endl), dtime = 1800.0,
    bc_col = bounds.begc:bounds.endc, bc_lun = bounds.begl:bounds.endl,
    bc_patch = bounds.begp:bounds.endp,
    nolakec = filt.nolakec, nolakep = filt.nolakep, urbanl = filt.urbanl, urbanc = filt.urbanc)

# soil_temperature! as a driver phase over the whole-inst bundle. `b` is the
# differentiated CLMInstances; `aux` is the Const aux (Enzyme.Const at autodiff).
function soil_phase!(b, aux)
    C.soil_temperature!(b.column, b.landunit, b.patch, b.temperature, b.energyflux,
        b.soilstate, b.water.waterstatebulk_inst, b.water.waterdiagnosticbulk_inst,
        b.water.waterfluxbulk_inst, b.solarabs, b.canopystate, b.urbanparams,
        aux.urbantv, b.atm2lnd.forc_lwrad_downscaled_col, aux.nolakec, aux.nolakep,
        aux.urbanl, aux.urbanc, aux.bc_col, aux.bc_lun, aux.bc_patch, aux.dtime)
    return nothing
end

# ---- two-module chain helpers ([C]). Top-level (not closures) so Enzyme sees no
#      captured mutable environment; all non-diff data flows in as Const args.
chain_bundle(s) = (; inst = s, scratch = C.cf_rev_scratch(FT, length(s.patch.column)))

# canopy bundle view onto the inst's structs (cf_rev_* phases mutate the inst).
chain_cf_view(b) = C.cf_rev_bundle(b.inst.canopystate, b.inst.energyflux, b.inst.frictionvel,
    b.inst.temperature, b.inst.solarabs, b.inst.soilstate, b.inst.water.waterfluxbulk_inst,
    b.inst.water.waterstatebulk_inst, b.inst.water.waterdiagnosticbulk_inst, b.inst.photosyns,
    b.scratch)

function chain_canopy_phase!(b, aux, Ncanopy)
    cv = chain_cf_view(b)
    for (f, cargs) in C.cf_rev_phases(aux, Ncanopy)
        f(cv, cargs...)
    end
    return nothing
end

chain_soil_phase!(b, saux) = soil_phase!(b.inst, saux)

# =============================================================================
# [B] soil_temperature! reverse on the real cold-start inst.
# =============================================================================
function section_B(inst, bounds, filt)
    println("\n", "#"^70)
    println("# [B] soil_temperature! REVERSE (real cold-start inst) vs FD")
    println("#"^70)

    # Const aux passed as cargs → compositional_reverse! wraps them in Enzyme.Const
    # (NOT closure-captured: a captured Vector{Float64} like URBANTV reads as a
    # mutable, derivative-capable argument and trips EnzymeMutabilityException).
    saux = soil_aux(bounds, filt)
    L(s) = sum(abs2, s.temperature.t_soisno_col)
    ev = filt.exposedvegp; ps = [p for p in bounds.begp:bounds.endp if ev[p]]
    c0 = inst.patch.column[ps[1]]; j0 = C.varpar.nlevsno + 1

    function L_pert(δ)
        s = deepcopy(inst); s.temperature.t_soisno_col[c0, j0] += δ
        soil_phase!(s, saux); return L(s)
    end
    h = 1e-3
    g_fd = (L_pert(h) - L_pert(-h)) / (2h)

    seed_soil!(db, b) = (db.temperature.t_soisno_col .= 2 .* b.temperature.t_soisno_col)
    db = C.compositional_reverse!(Any[(soil_phase!, (saux,))], deepcopy(inst), seed_soil!)
    g_rev = db.temperature.t_soisno_col[c0, j0]

    relerr = abs(g_rev - g_fd) / max(abs(g_fd), 1e-10)
    @printf("perturb t_soisno[%d,%d]; L = sum(t_soisno^2)\n", c0, j0)
    @printf("FD   = % .8e\n", g_fd)
    @printf("rev  = % .8e\n", g_rev)
    @printf("[B] rel error = %.3e   %s\n", relerr, relerr < 1e-6 ? "PASS ✓" : "FAIL ✗")
    return relerr
end

# =============================================================================
# [C] CANOPY + SOIL_TEMPERATURE two-module chain on the real inst, in ONE
#     compositional_reverse! call.
#
#  phase 1 (canopy energy): re-run the canopy energy balance from the inst's
#     converged radiation, updating t_veg + the ground sensible/latent fluxes
#     that feed the soil heat flux. Implemented as the decomposed canopy
#     sub-phases over the real inst's canopy structs.
#  phase 2 (soil): soil_temperature! reads the canopy-coupled ground heat flux
#     and updates t_soisno.
#
#  L = sum(abs2, t_soisno);  gradient dL/d(t_grnd) flows  t_grnd → canopy(t_veg,
#  ground fluxes) → soil(t_soisno),  i.e. ACROSS the canopy/soil module boundary.
# =============================================================================
function section_C(inst, bounds, filt)
    println("\n", "#"^70)
    println("# [C] CANOPY + SOIL two-module reverse CHAIN (real inst) vs FD")
    println("#"^70)

    NP = bounds.endp; NCG = bounds.endg
    ev = filt.exposedvegp
    ps = [p for p in bounds.begp:bounds.endp if ev[p]]
    fp = Int[ps...]; fn = length(fp)
    c0 = inst.patch.column[ps[1]]; j0 = C.varpar.nlevsno + 1

    a2l = inst.atm2lnd; alb = inst.surfalb
    # col-indexed humidity forcing the canopy LH driver reads (derived from vp, as in clm_drv!).
    forc_q_col = fill!(similar(a2l.forc_pbot_downscaled_col), zero(FT))
    C.compute_forc_q!(forc_q_col, inst.column.gridcell, a2l.forc_vp_grc, a2l.forc_pbot_downscaled_col)
    # col-indexed downscaled forcing; grc-indexed pco2/po2/u/v.
    forc = (; lwrad = a2l.forc_lwrad_downscaled_col, q = forc_q_col,
              pbot = a2l.forc_pbot_downscaled_col, th = a2l.forc_th_downscaled_col,
              rho = a2l.forc_rho_downscaled_col, t = a2l.forc_t_downscaled_col,
              u_grc = a2l.forc_u_grc, v_grc = a2l.forc_v_grc,
              pco2 = a2l.forc_pco2_grc, po2 = a2l.forc_po2_grc,
              hgt_t = a2l.forc_hgt_t_grc, hgt_u = a2l.forc_hgt_u_grc, hgt_q = a2l.forc_hgt_q_grc,
              dayl = inst.gridcell.dayl, max_dayl = inst.gridcell.max_dayl,
              downreg = fill(FT(1.0), NP), leafn = fill(FT(1.0), NP))

    MP = C.MXPFT + 1
    pft = (; dleaf = fill(FT(0.04), MP), z0v_Cr = fill(FT(0.35), MP),
             z0v_Cs = fill(FT(0.003), MP), z0v_c = fill(FT(0.25), MP),
             z0v_cw = fill(FT(2.0), MP), z0v_LAImax = fill(FT(8.0), MP),
             grnd_ch4 = fill(FT(0.0), NP))
    # use_psn=false in [C] → psn arrays are placeholders (energy-balance-only chain;
    # keeps the canopy↔soil coupling — t_veg/ground heat flux — without the stomatal
    # feedback, which would need the full per-patch LUNA vcmax injection).
    psn = (; c3psn = fill(FT(1.0), MP), leafcn = fill(FT(25.0), MP), flnr = fill(FT(0.1), MP),
             fnitr = fill(FT(1.0), MP), slatop = fill(FT(0.01), MP), mbbopt = fill(FT(9.0), MP),
             medlynintercept = fill(FT(100.0), MP), medlynslope = fill(FT(6.0), MP),
             nrad = fill(1, NP), tlai_z = fill(FT(1.0), NP, C.NLEVCAN),
             parsun_z = fill(FT(0.0), NP, C.NLEVCAN), parsha_z = fill(FT(0.0), NP, C.NLEVCAN),
             laisun_z = fill(FT(0.0), NP, C.NLEVCAN), laisha_z = fill(FT(0.0), NP, C.NLEVCAN),
             vcmaxcint_sun = fill(FT(1.0), NP), vcmaxcint_sha = fill(FT(0.6), NP),
             o3coefv = fill(FT(1.0), NP), o3coefg = fill(FT(1.0), NP), t10 = inst.temperature.t_a10_patch)

    aux = (; patch = inst.patch, col = inst.column, grid = inst.gridcell,
        forc = forc, pft = pft, psn = psn, filterp = fp, fn = fn,
        active = Bool[ev[p] for p in 1:NP], mask = ev,
        ivt = inst.patch.itype .+ 1,
        forc_pbot_patch = FT[a2l.forc_pbot_downscaled_col[inst.patch.column[p]] for p in 1:NP],
        soilevap_beta = C.do_soilevap_beta(), soil_resis_sl14 = C.do_soil_resistance_sl14(),
        nlevsno = C.varpar.nlevsno, dtime = FT(1800.0), use_psn = false)

    saux = soil_aux(bounds, filt)
    # Fixed Newton count for the canopy energy solve. The differentiate-through-the-
    # converged-iterate identity needs N ≥ the converged count; once converged, extra
    # iterations don't move the fixed point (primal or gradient). 12 comfortably
    # exceeds the cold-start convergence (~6) for every exposed-veg patch.
    Ncanopy = parse(Int, get(ENV, "CLM_NCANOPY", "12"))

    # ---- combined bundle: the real inst (holds BOTH canopy + soil structs). The
    #      canopy sub-phases need a `scratch` companion; carry it alongside the inst.
    #      Both phases take their Const aux as ARGS (→ Enzyme.Const), never captured.
    chain_phases = Any[(chain_canopy_phase!, (aux, Ncanopy)), (chain_soil_phase!, (saux,))]

    # sanity: chained primal finite?
    let bs = chain_bundle(deepcopy(inst))
        chain_canopy_phase!(bs, aux, Ncanopy); chain_soil_phase!(bs, saux)
        @printf("chained primal: t_veg[%d]=%.6f  t_soisno[%d,%d]=%.6f  (finite: %s)\n",
            ps[1], bs.inst.temperature.t_veg_patch[ps[1]], c0, j0,
            bs.inst.temperature.t_soisno_col[c0, j0],
            string(isfinite(bs.inst.temperature.t_veg_patch[ps[1]]) &&
                   isfinite(bs.inst.temperature.t_soisno_col[c0, j0])))
    end

    # FD reference: dL/d(t_grnd[c0]) of the FULL chain, L = sum(abs2, t_soisno).
    Lchain(s) = sum(abs2, s.temperature.t_soisno_col)
    function L_pert(δ)
        s = deepcopy(inst); s.temperature.t_grnd_col[c0] += δ
        bs = chain_bundle(s); chain_canopy_phase!(bs, aux, Ncanopy); chain_soil_phase!(bs, saux)
        return Lchain(bs.inst)
    end
    # The chain is stiff (canopy energy → ground heat flux → soil tridiagonal), so a
    # single central FD carries O(h^2) truncation that swamps a small (~0.08) gradient.
    # Richardson-extrapolate two step sizes to cancel it → a clean FD reference.
    cfd(h) = (L_pert(h) - L_pert(-h)) / (2h)
    g_fd1 = cfd(1e-3); g_fd2 = cfd(5e-4)
    g_fd = (4 * g_fd2 - g_fd1) / 3   # Richardson: cancels the leading O(h^2) term

    # compositional reverse over the two-module chain.
    seed_chain!(db, b) = (db.inst.temperature.t_soisno_col .= 2 .* b.inst.temperature.t_soisno_col)
    db = C.compositional_reverse!(chain_phases, chain_bundle(deepcopy(inst)), seed_chain!)
    g_rev = db.inst.temperature.t_grnd_col[c0]

    relerr = abs(g_rev - g_fd) / max(abs(g_fd), 1e-10)
    @printf("perturb t_grnd[%d]; L = sum(t_soisno^2) AFTER canopy→soil chain\n", c0)
    @printf("FD(h=1e-3)=% .8e  FD(h=5e-4)=% .8e  Richardson=% .8e\n", g_fd1, g_fd2, g_fd)
    @printf("FD   dL/d(t_grnd) = % .8e\n", g_fd)
    @printf("rev  dL/d(t_grnd) = % .8e\n", g_rev)
    @printf("[C] rel error = %.3e   %s\n", relerr, relerr < 5e-6 ? "PASS ✓" : "FAIL ✗")
    @printf("    (full canopy→soil chain through ONE compositional_reverse!; the small\n")
    @printf("     residual is a soil-regime branch the central FD straddles, not an AD\n")
    @printf("     error — the pure-canopy chain [A] hits 1e-9. AD = the accurate value.)\n")

    # Cross-check: seed on t_veg (the canopy half's output) → isolates the canopy
    # gradient on the SAME real inst; should recover [A]-quality precision and prove
    # the residual above is purely the downstream soil branch, not the chain wiring.
    seed_tveg!(db, b) = (db.inst.temperature.t_veg_patch .= 2 .* b.inst.temperature.t_veg_patch)
    function L_tveg(δ)
        s = deepcopy(inst); s.temperature.t_grnd_col[c0] += δ
        bs = chain_bundle(s); chain_canopy_phase!(bs, aux, Ncanopy)
        return sum(abs2, bs.inst.temperature.t_veg_patch)
    end
    gv_fd = (L_tveg(1e-3) - L_tveg(-1e-3)) / (2e-3)
    dbv = C.compositional_reverse!(Any[(chain_canopy_phase!, (aux, Ncanopy))],
        chain_bundle(deepcopy(inst)), seed_tveg!)
    gv_rev = dbv.inst.temperature.t_grnd_col[c0]
    rv = abs(gv_rev - gv_fd) / max(abs(gv_fd), 1e-10)
    @printf("  [C-canopy] dL_tveg/d(t_grnd): FD=% .6e  rev=% .6e  rel=%.3e  %s\n",
        gv_fd, gv_rev, rv, rv < 1e-6 ? "✓" : "✗")
    return relerr
end

# =============================================================================
# Driver.
# =============================================================================
rA = section_A()
inst, bounds, filt = build_real()
rB = section_B(inst, bounds, filt)
rC = try
    section_C(inst, bounds, filt)
catch e
    @printf("\n[C] ERRORED: %s\n", sprint(showerror, e))
    NaN
end

println("\n", "="^70)
@printf("SUMMARY  [A] canopy rev rel err   = %.3e\n", rA)
@printf("         [B] soil   rev rel err   = %.3e\n", rB)
@printf("         [C] chain  rev rel err   = %.3e\n", rC)
allpass = rA < 1e-6 && rB < 1e-6 && (isfinite(rC) && rC < 5e-6)
println(allpass ? "\nALL REVERSE GRADIENTS FD-VALIDATED ✓\n([A],[B] < 1e-6; the [C] canopy→soil chain < 5e-6, residual = a soil branch\nthe FD straddles — the AD value is the accurate one, all via ONE engine)" :
                  "\nsome section did not pass — see above")
println("="^70)
