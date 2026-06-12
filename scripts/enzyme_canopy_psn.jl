# Whole-canopy compositional reverse-AD (+photosynthesis) — validation of the
# PRODUCTIONIZED src API (src/biogeophys/canopy_fluxes_reverse.jl).
#
# Builds the finite single-patch state, runs real canopy_fluxes!(+psn) as the
# reference, then exercises CLM.canopy_rev_forward! (primal parity) and
# CLM.canopy_rev_gradient! (compositional reverse), validating dL/d(t_grnd) vs
# central finite differences. The decomposition/checkpointing now lives in src;
# this script is just the Enzyme-on-1.10 validation harness + finite-state builder.
#
#   julia +1.10 --project=/tmp/clm_jl10_env scripts/enzyme_canopy_psn.jl
#
using CLM, Printf
const C = CLM
const FT = Float64
C.varpar_init!(C.varpar, 1, 14, 2, 5)
const NP = 1; const NC = 1; const NG = 1; const NL = 1
const NLEVCAN = C.NLEVCAN; const MP = C.MXPFT + 1

C.soil_resistance_read_nl!(soil_resis_method = C.SOIL_RESIS_LEEPIELKE_1992)
C.canopy_fluxes_read_nml!(use_undercanopy_stability = false,
    use_biomass_heat_storage = false, itmax_canopy_fluxes = 40)
C.canopy_fluxes_read_params!()

# photosynthesis scalar params (else photosynthesis! is NaN)
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

const PSN = (; c3psn = fill(FT(1.0), MP), leafcn = fill(FT(25.0), MP), flnr = fill(FT(0.1), MP),
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

function run_real!(S, forc, pft, dtime)
    C.canopy_fluxes!(S.canopystate, S.energyflux, S.frictionvel, S.temperature, S.solarabs,
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
end

# aux + bundle for the src compositional-reverse API.
function make_aux(B)
    (; patch = B.S.patch_data, col = B.S.col_data, grid = B.S.gridcell_data,
       forc = B.forc, pft = B.pft, psn = PSN, filterp = Int[1], fn = 1, active = trues(NP),
       mask = trues(NP), ivt = B.S.patch_data.itype .+ 1,
       forc_pbot_patch = FT[B.forc.pbot[B.S.patch_data.column[p]] for p in 1:NP],
       soilevap_beta = C.do_soilevap_beta(), soil_resis_sl14 = C.do_soil_resistance_sl14(),
       nlevsno = C.varpar.nlevsno, dtime = B.dtime, use_psn = true)
end
make_b(B) = C.cf_rev_bundle(B.S.canopystate, B.S.energyflux, B.S.frictionvel, B.S.temperature,
    B.S.solarabs, B.S.soilstate, B.S.waterfluxbulk, B.S.waterstatebulk, B.S.waterdiagbulk,
    B.S.photosyns, C.cf_rev_scratch(FT, NP))

# ---- reference + decomposed-forward parity (via src) ----
Bref = build(); run_real!(Bref.S, Bref.forc, Bref.pft, Bref.dtime)
const N = max(Int(Bref.S.frictionvel.num_iter_patch[1]), 1)
@printf("REAL canopy_fluxes!(+psn):  t_veg=%.10f  rssun=%.6f  num_iter=%d\n",
    Bref.S.temperature.t_veg_patch[1], Bref.S.photosyns.rssun_patch[1], N)
let B = build(); b = make_b(B); C.canopy_rev_forward!(b, make_aux(B), N)
    @printf("src canopy_rev_forward! :    t_veg=%.10f  rssun=%.6f  |Δt_veg|=%.3e\n",
        b.temperature.t_veg_patch[1], b.photosyns.rssun_patch[1],
        abs(b.temperature.t_veg_patch[1] - Bref.S.temperature.t_veg_patch[1]))
end

# ---- FD reference dL/d(t_grnd) ----
function L_tgrnd(δ)
    B = build(); b = make_b(B); b.temperature.t_grnd_col[1] += δ
    C.canopy_rev_forward!(b, make_aux(B), N); return sum(abs2, b.temperature.t_veg_patch)
end
hg = 1e-2
g_fd = (L_tgrnd(hg) - L_tgrnd(-hg)) / (2hg)
@printf("\nFD   dL/d(t_grnd)     = % .10e\n", g_fd)

# ---- compositional reverse (via src) ----
let B = build(); b = make_b(B)
    db = C.canopy_rev_gradient!(b, make_aux(B), N; seed = :t_veg_patch)
    g_comp = db.temperature.t_grnd_col[1]
    @printf("comp dL/d(t_grnd)     = % .10e\n", g_comp)
    println("\n", "="^64)
    relerr = abs(g_comp - g_fd) / max(abs(g_fd), 1e-10)
    @printf("abs error = %.3e   rel error = %.3e\n", abs(g_comp - g_fd), relerr)
    if isfinite(g_comp) && isfinite(g_fd) && relerr < 1e-4
        println("\nWHOLE-CANOPY (+PHOTOSYNTHESIS) COMPOSITIONAL REVERSE-AD VALIDATED ✓ (src API)")
        println("init + $N iters × (friction,resist,psn_sun,psn_sha,energy); gradient incl.")
        println("the stomatal feedback matches finite differences.")
    else
        println("\nMISMATCH ✗ — investigate")
    end
    println("="^64)
end
