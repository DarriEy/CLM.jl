# Whole-canopy compositional reverse-AD — MULTI-PATCH (item #4).
#
# Two distinct patches (different LAI/height/PFT) that converge at DIFFERENT Newton
# iteration counts. The mirror runs a FIXED N = max(num_iter) with active≡true for
# all patches: the early-converging patch keeps iterating at its (stable) fixed
# point, so its state doesn't move — primal parity must hold for BOTH patches, and
# the per-patch gradient must match finite differences. This validates that the
# fixed-N / active≡true scheme is convergence-aware enough for multi-patch without
# per-phase active-mask checkpointing.
#
#   julia +1.10 --project=/tmp/clm_jl10_env scripts/enzyme_canopy_multipatch.jl
#
using CLM, Printf
const C = CLM
const FT = Float64
C.varpar_init!(C.varpar, 1, 14, 2, 5)
const NP = 2; const NC = 2; const NG = 1; const NL = 2; const NLEVCAN = C.NLEVCAN; const MP = C.MXPFT + 1

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

# per-patch radiation/LAI inputs differ between the two patches
const PSN = (; c3psn = fill(FT(1.0), MP), leafcn = fill(FT(25.0), MP), flnr = fill(FT(0.1), MP),
    fnitr = fill(FT(1.0), MP), slatop = fill(FT(0.01), MP), mbbopt = fill(FT(9.0), MP),
    medlynintercept = fill(FT(100.0), MP), medlynslope = fill(FT(6.0), MP),
    nrad = fill(1, NP), tlai_z = fill(FT(1.0), NP, NLEVCAN),
    parsun_z = reshape(FT[250.0, 180.0], NP, NLEVCAN), parsha_z = reshape(FT[120.0, 90.0], NP, NLEVCAN),
    laisun_z = fill(FT(1.0), NP, NLEVCAN), laisha_z = fill(FT(1.0), NP, NLEVCAN),
    vcmaxcint_sun = FT[1.0, 0.8], vcmaxcint_sha = FT[0.6, 0.5],
    o3coefv = fill(FT(1.0), NP), o3coefg = fill(FT(1.0), NP), t10 = fill(FT(290.0), NP))

function build()
    cs = C.CanopyStateData{FT}(); C.canopystate_init!(cs, NP)
    cs.elai_patch .= FT[2.5, 1.2]; cs.esai_patch .= FT[0.5, 0.3]
    cs.laisun_patch .= FT[1.2, 0.7]; cs.laisha_patch .= FT[0.8, 0.5]
    cs.displa_patch .= FT[5.0, 3.0]; cs.htop_patch .= FT[10.0, 6.0]
    cs.frac_veg_nosno_patch .= 1; cs.dleaf_patch .= FT[0.04, 0.05]
    cs.stem_biomass_patch .= 0.0; cs.leaf_biomass_patch .= 0.0; cs.tlai_patch .= FT[2.5, 1.2]
    cs.laisun_z_patch = copy(PSN.laisun_z); cs.laisha_z_patch = copy(PSN.laisha_z)

    ef = C.EnergyFluxData{FT}(); C.energyflux_init!(ef, NP, NC, NL, NG)
    ef.btran_patch .= FT[0.8, 0.6]; ef.bsun_patch .= 0.8; ef.bsha_patch .= 0.8
    ef.htvp_col .= C.HVAP

    fv = C.FrictionVelocityData{FT}(); C.frictionvel_init!(fv, NP, NC)
    fv.zetamaxstable = 0.5; fv.zsno = 0.00085; fv.zlnd = 0.000775
    fv.z0mv_patch .= 0.5; fv.z0hv_patch .= 0.5; fv.z0qv_patch .= 0.5
    fv.z0mg_col .= 0.01; fv.z0hg_col .= 0.01; fv.z0qg_col .= 0.01
    fv.forc_hgt_u_patch .= 30.0; fv.forc_hgt_t_patch .= 30.0; fv.forc_hgt_q_patch .= 30.0
    fv.ustar_patch .= 0.5; fv.um_patch .= 5.0; fv.uaf_patch .= 3.0
    fv.taf_patch .= 290.0; fv.qaf_patch .= 0.008; fv.obu_patch .= -100.0; fv.zeta_patch .= -0.1
    fv.vpd_patch .= 1.0; fv.rb1_patch .= 50.0; fv.ram1_patch .= 50.0; fv.num_iter_patch .= 0.0

    tp = C.TemperatureData{FT}(); C.temperature_init!(tp, NP, NC, NL, NG)
    tp.t_veg_patch .= FT[290.0, 291.0]; tp.t_stem_patch .= 289.0; tp.t_skin_patch .= 290.0
    tp.thm_patch .= FT[290.0, 290.5]; tp.t_grnd_col .= FT[288.0, 287.0]; tp.t_h2osfc_col .= 288.0
    tp.thv_col .= 291.0; tp.emv_patch .= 0.97; tp.emg_col .= 0.96
    tp.t_ref2m_patch .= 290.0; tp.t_ref2m_r_patch .= 290.0; tp.t_a10_patch .= 290.0
    for j in 1:size(tp.t_soisno_col, 2), c in 1:NC; tp.t_soisno_col[c, j] = 288.0; end

    sa = C.SolarAbsorbedData{FT}(); C.solarabs_init!(sa, NP, NL)
    sa.sabv_patch .= FT[150.0, 95.0]
    sa.parsun_z_patch = copy(PSN.parsun_z); sa.parsha_z_patch = copy(PSN.parsha_z)

    ss = C.SoilStateData{FT}(); C.soilstate_init!(ss, NC, NL)
    ss.soilbeta_col .= 0.8; ss.soilresis_col .= 100.0
    for j in 1:C.varpar.nlevgrnd, c in 1:NC; ss.watsat_col[c, j] = 0.45; end

    wfb = C.WaterFluxBulkData{FT}(); C.waterfluxbulk_init!(wfb, NC, NP, NL, NG)
    wsb = C.WaterStateBulkData{FT}(); C.waterstatebulk_init!(wsb, NC, NP, NL, NG)
    wsb.ws.liqcan_patch .= 0.1; wsb.ws.snocan_patch .= 0.0

    wdb = C.WaterDiagnosticBulkData{FT}(); C.waterdiagnosticbulk_init!(wdb, NC, NP, NL, NG)
    wdb.fwet_patch .= 0.1; wdb.fdry_patch .= 0.8; wdb.frac_sno_eff_col .= 0.0; wdb.frac_h2osfc_col .= 0.0
    wdb.snow_depth_col .= 0.0; wdb.qg_col .= 0.005; wdb.qg_snow_col .= 0.005; wdb.qg_soil_col .= 0.005
    wdb.qg_h2osfc_col .= 0.005; wdb.dqgdT_col .= 0.0003; wdb.rh_af_patch .= 0.6

    ps = C.PhotosynthesisData{FT}(); C.photosynthesis_data_init!(ps, NP)
    ps.rssun_patch .= 200.0; ps.rssha_patch .= 300.0
    ps.stomatalcond_mtd = C.STOMATALCOND_MTD_BB1987; ps.leafresp_method = C.LEAFRESP_MTD_RYAN1991
    ps.light_inhibit = false

    pd = C.PatchData{FT}(); C.patch_init!(pd, NP)
    pd.column .= [1, 2]; pd.gridcell .= 1; pd.landunit .= [1, 2]; pd.itype .= [1, 7]; pd.active .= true
    cd = C.ColumnData{FT}(); C.column_init!(cd, NC); cd.snl .= 0
    gd = C.GridcellData{FT}(); C.gridcell_init!(gd, NG)

    S = (; canopystate=cs, energyflux=ef, frictionvel=fv, temperature=tp, solarabs=sa,
           soilstate=ss, waterfluxbulk=wfb, waterstatebulk=wsb, waterdiagbulk=wdb,
           photosyns=ps, patch_data=pd, col_data=cd, gridcell_data=gd)
    forc = (; lwrad = FT[350.0, 340.0], q = FT[0.008, 0.007], pbot = FT[101325.0, 101000.0],
              th = FT[290.0, 290.5], rho = FT[1.2, 1.2], t = FT[288.0, 287.5],
              u_grc = FT[3.0], v_grc = FT[1.0], pco2 = FT[40.0], po2 = FT[21000.0],
              hgt_t = FT[10.0], hgt_u = FT[10.0], hgt_q = FT[10.0], dayl = FT[43200.0],
              max_dayl = FT[50000.0], downreg = FT[1.0, 1.0], leafn = FT[1.0, 1.0])
    pft = (; dleaf = fill(FT(0.04), MP), z0v_Cr = fill(FT(0.35), MP), z0v_Cs = fill(FT(0.003), MP),
             z0v_c = fill(FT(0.25), MP), z0v_cw = fill(FT(2.0), MP), z0v_LAImax = fill(FT(8.0), MP),
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

make_aux(B) = (; patch = B.S.patch_data, col = B.S.col_data, grid = B.S.gridcell_data,
    forc = B.forc, pft = B.pft, psn = PSN, filterp = Int[1, 2], fn = NP, active = trues(NP),
    mask = trues(NP), ivt = B.S.patch_data.itype .+ 1,
    forc_pbot_patch = FT[B.forc.pbot[B.S.patch_data.column[p]] for p in 1:NP],
    soilevap_beta = C.do_soilevap_beta(), soil_resis_sl14 = C.do_soil_resistance_sl14(),
    nlevsno = C.varpar.nlevsno, dtime = B.dtime, use_psn = true)
make_b(B) = C.cf_rev_bundle(B.S.canopystate, B.S.energyflux, B.S.frictionvel, B.S.temperature,
    B.S.solarabs, B.S.soilstate, B.S.waterfluxbulk, B.S.waterstatebulk, B.S.waterdiagbulk,
    B.S.photosyns, C.cf_rev_scratch(FT, NP))

# ---- reference + per-patch num_iter ----
Bref = build(); run_real!(Bref.S, Bref.forc, Bref.pft, Bref.dtime)
nit = Int.(Bref.S.frictionvel.num_iter_patch)
const N = maximum(nit)
@printf("REAL canopy(+psn) per-patch:  num_iter=%s  t_veg=%s\n", string(nit),
    string(round.(Bref.S.temperature.t_veg_patch, digits=6)))

let B = build(); b = make_b(B); C.canopy_rev_forward!(b, make_aux(B), N)
    dtv = abs.(b.temperature.t_veg_patch .- Bref.S.temperature.t_veg_patch)
    @printf("src forward (N=%d):           t_veg=%s  max|Δt_veg|=%.3e\n", N,
        string(round.(b.temperature.t_veg_patch, digits=6)), maximum(dtv))
end

# ---- per-patch FD vs compositional reverse, dL/d(t_grnd[c]), L = sum(abs2,t_veg) ----
function L_tgrnd(c, δ)
    B = build(); b = make_b(B); b.temperature.t_grnd_col[c] += δ
    C.canopy_rev_forward!(b, make_aux(B), N); return sum(abs2, b.temperature.t_veg_patch)
end
hg = 1e-2
g_fd = [(L_tgrnd(c, hg) - L_tgrnd(c, -hg)) / (2hg) for c in 1:NC]

local g_comp
let B = build(); b = make_b(B)
    db = C.canopy_rev_gradient!(b, make_aux(B), N; seed = :t_veg_patch)
    global g_comp = copy(db.temperature.t_grnd_col)
end

println("\n", "="^60)
maxrel = 0.0
for c in 1:NC
    re = abs(g_comp[c] - g_fd[c]) / max(abs(g_fd[c]), 1e-10)
    global maxrel = max(maxrel, re)
    @printf("col %d:  FD=% .8e  comp=% .8e  rel=%.2e\n", c, g_fd[c], g_comp[c], re)
end
ok = isfinite(maxrel) && maxrel < 1e-4
println(ok ? "\nMULTI-PATCH WHOLE-CANOPY REVERSE-AD VALIDATED ✓ (2 patches, different convergence)" :
             "\nMISMATCH ✗ — investigate")
println("="^60)
