# Whole-canopy compositional reverse-AD — Stage 1: decomposed-mirror PRIMAL PARITY.
#
# A standalone "mirror" of canopy_fluxes_core!'s ENERGY-BALANCE path (no
# photosynthesis), built by calling the SAME init kernels + Newton-loop sub-phase
# functions in the SAME order. The mirror is the basis for the compositional reverse
# pass (Stage 2); here we only prove the mirror's PRIMAL matches the real
# canopy_fluxes! on the finite single-patch state — so the eventual reverse gradient
# is the REAL canopy gradient, not an approximation.
#
# Simplifications (all gradient-preserving for L = sum(abs2, t_veg)):
#   - post-iteration diagnostics skipped (they don't feed t_veg)
#   - photosynthesis / biomass-heat / daylength kernels skipped (energy-balance path)
#   - fixed N Newton iterations, active≡true (no convergence early-exit); once the
#     solve has converged, extra iterations don't move the fixed point.
#
#   julia +1.10 --project=/tmp/clm_jl10_env scripts/enzyme_canopy_mirror.jl
#
using CLM, Printf
const C = CLM
const FT = Float64
C.varpar_init!(C.varpar, 1, 14, 2, 5)
const NP = 1; const NC = 1; const NG = 1; const NL = 1

C.soil_resistance_read_nl!(soil_resis_method = C.SOIL_RESIS_LEEPIELKE_1992)
C.canopy_fluxes_read_nml!(use_undercanopy_stability = false,
    use_biomass_heat_storage = false, itmax_canopy_fluxes = 40)
C.canopy_fluxes_read_params!()

# ---- finite single-patch state (from gpu_validate_canopy_e2e.jl build()) ----
function build()
    canopystate = C.CanopyStateData{FT}(); C.canopystate_init!(canopystate, NP)
    canopystate.elai_patch[1] = 2.0;  canopystate.esai_patch[1] = 0.5
    canopystate.laisun_patch[1] = 1.2; canopystate.laisha_patch[1] = 0.8
    canopystate.displa_patch[1] = 5.0; canopystate.htop_patch[1] = 10.0
    canopystate.frac_veg_nosno_patch[1] = 1; canopystate.dleaf_patch[1] = 0.04
    canopystate.stem_biomass_patch[1] = 0.0; canopystate.leaf_biomass_patch[1] = 0.0

    energyflux = C.EnergyFluxData{FT}(); C.energyflux_init!(energyflux, NP, NC, NL, NG)
    energyflux.btran_patch[1] = 0.5; energyflux.bsun_patch[1] = 0.5; energyflux.bsha_patch[1] = 0.5
    energyflux.htvp_col[1] = C.HVAP
    energyflux.cgrnds_patch[1] = 0.0; energyflux.cgrndl_patch[1] = 0.0; energyflux.cgrnd_patch[1] = 0.0
    energyflux.dhsdt_canopy_patch[1] = 0.0; energyflux.eflx_sh_stem_patch[1] = 0.0

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
    for j in 1:size(temperature.t_soisno_col, 2); temperature.t_soisno_col[1, j] = 288.0; end

    solarabs = C.SolarAbsorbedData{FT}(); C.solarabs_init!(solarabs, NP, NL)
    solarabs.sabv_patch[1] = 150.0

    soilstate = C.SoilStateData{FT}(); C.soilstate_init!(soilstate, NC, NL)
    soilstate.soilbeta_col[1] = 0.8; soilstate.soilresis_col[1] = 100.0
    for j in 1:C.varpar.nlevgrnd; soilstate.watsat_col[1, j] = 0.45; end

    waterfluxbulk = C.WaterFluxBulkData{FT}(); C.waterfluxbulk_init!(waterfluxbulk, NC, NP, NL, NG)
    waterfluxbulk.wf.qflx_tran_veg_patch[1] = 0.0; waterfluxbulk.wf.qflx_evap_veg_patch[1] = 0.0
    waterfluxbulk.wf.qflx_evap_soi_patch[1] = 0.0

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
    mp = C.MXPFT + 1
    pft = (; dleaf = fill(FT(0.04), mp), z0v_Cr = fill(FT(0.35), mp),
             z0v_Cs = fill(FT(0.003), mp), z0v_c = fill(FT(0.25), mp),
             z0v_cw = fill(FT(2.0), mp), z0v_LAImax = fill(FT(8.0), mp),
             grnd_ch4 = fill(FT(0.0), NP))
    return (; S, forc, pft, dtime = FT(1800.0))
end

# real canopy_fluxes! (the reference) on a copy
run_real!(S, forc, pft, dtime) = C.canopy_fluxes!(
    S.canopystate, S.energyflux, S.frictionvel, S.temperature, S.solarabs,
    S.soilstate, S.waterfluxbulk, S.waterstatebulk, S.waterdiagbulk, S.photosyns,
    S.patch_data, S.col_data, S.gridcell_data, trues(NP), 1:NP, 1:NC,
    forc.lwrad, forc.q, forc.pbot, forc.th, forc.rho, forc.t,
    forc.u_grc, forc.v_grc, forc.pco2, forc.po2, forc.hgt_t, forc.hgt_u, forc.hgt_q,
    forc.dayl, forc.max_dayl, forc.downreg, forc.leafn, dtime;
    dleaf_pft = pft.dleaf, z0v_Cr_pft = pft.z0v_Cr, z0v_Cs_pft = pft.z0v_Cs,
    z0v_c_pft = pft.z0v_c, z0v_cw_pft = pft.z0v_cw, z0v_LAImax_pft = pft.z0v_LAImax,
    grnd_ch4_cond_patch = pft.grnd_ch4)

# ---------------------------------------------------------------------------
# The decomposed MIRROR (energy-balance path, fixed N iterations, active≡true).
# Scratch held in a flat NamedTuple so Stage 2 can checkpoint/shadow it.
# ---------------------------------------------------------------------------
mkscr() = (; zldis=zeros(FT,NP), dth=zeros(FT,NP), dthv=zeros(FT,NP), dqh=zeros(FT,NP),
    ur=zeros(FT,NP), temp1=zeros(FT,NP), temp12m=zeros(FT,NP), temp2=zeros(FT,NP),
    temp22m=zeros(FT,NP), rb=zeros(FT,NP), rah=zeros(FT,NP,2), raw=zeros(FT,NP,2),
    wtg=zeros(FT,NP), wta0=zeros(FT,NP), wtl0=zeros(FT,NP), wtstem0=zeros(FT,NP),
    wtal=zeros(FT,NP), wtga=zeros(FT,NP), wtgq=zeros(FT,NP), wtaq0=zeros(FT,NP),
    wtlq0=zeros(FT,NP), wtalq=zeros(FT,NP), el=zeros(FT,NP), qsatl=zeros(FT,NP),
    qsatldT=zeros(FT,NP), air=zeros(FT,NP), bir=zeros(FT,NP), cir=zeros(FT,NP),
    del_arr=zeros(FT,NP), del2=zeros(FT,NP), dele=zeros(FT,NP), delq=zeros(FT,NP),
    det_arr=zeros(FT,NP), efeb=zeros(FT,NP), efe=zeros(FT,NP), err_arr=zeros(FT,NP),
    obuold=zeros(FT,NP), tlbef=zeros(FT,NP), tl_ini=zeros(FT,NP), ts_ini=zeros(FT,NP),
    co2_arr=zeros(FT,NP), o2_arr=zeros(FT,NP), svpts=zeros(FT,NP), eah=zeros(FT,NP),
    dt_veg=zeros(FT,NP), fm=zeros(FT,NP), nmozsgn=zeros(Int,NP),
    cp_leaf=zeros(FT,NP), rstem=zeros(FT,NP), frac_rad_abs_by_stem=zeros(FT,NP),
    sa_stem=zeros(FT,NP), sa_leaf=zeros(FT,NP), sa_internal=zeros(FT,NP), uuc=zeros(FT,NP),
    lw_stem=zeros(FT,NP), lw_leaf=zeros(FT,NP))

const PARAMS = C.canopy_fluxes_params
const SOILEVAP_BETA = C.do_soilevap_beta()
const SOIL_RESIS_SL14 = C.do_soil_resistance_sl14()
const NLEVSNO = C.varpar.nlevsno

# init phase: the 5 energy-balance-relevant init kernels.
function canopy_init!(S, forc, pft, sc, filterp, fn)
    cs = S.canopystate; fv = S.frictionvel; tp = S.temperature; pd = S.patch_data
    C._launch!(C._cf_init_zero_kernel!, sc.del_arr, sc.efeb, sc.wtlq0, sc.wtalq, sc.wtgq,
        sc.wtaq0, sc.obuold, S.energyflux.dhsdt_canopy_patch, S.energyflux.eflx_sh_stem_patch,
        filterp; ndrange = fn)
    C._launch!(C._cf_no_biomass_kernel!, sc.frac_rad_abs_by_stem, sc.sa_leaf, sc.sa_stem,
        sc.sa_internal, sc.cp_leaf, zeros(FT, NP), sc.rstem, filterp, cs.elai_patch, cs.esai_patch;
        ndrange = fn)
    fv.rb1_patch[1:NP] .= zero(FT)
    C._launch!(C._cf_z0_kernel!, cs.displa_patch, fv.z0mv_patch, fv.z0hv_patch, fv.z0qv_patch,
        fv.forc_hgt_u_patch, fv.forc_hgt_t_patch, fv.forc_hgt_q_patch, filterp, pd.column,
        pd.gridcell, pd.itype, cs.elai_patch, cs.esai_patch, cs.htop_patch, fv.z0mg_col,
        forc.hgt_u, forc.hgt_t, forc.hgt_q, pft.z0v_LAImax, pft.z0v_Cs, pft.z0v_Cr,
        pft.z0v_c, pft.z0v_cw, 1; ndrange = fn)
    lwq_out = C.CfLwqOut(; air=sc.air, bir=sc.bir, cir=sc.cir, qsatl=sc.qsatl, el=sc.el,
        qsatldT=sc.qsatldT, co2_arr=sc.co2_arr, o2_arr=sc.o2_arr, taf=fv.taf_patch,
        qaf=fv.qaf_patch, ur=sc.ur, dth=sc.dth, dqh=sc.dqh, delq=sc.delq, dthv=sc.dthv, zldis=sc.zldis)
    lwq_p = C.CfLwqP(; emv=tp.emv_patch, t_veg=tp.t_veg_patch, thm=tp.thm_patch,
        forc_hgt_u=fv.forc_hgt_u_patch, displa=cs.displa_patch)
    lwq_c = C.CfLwqC(; emg=tp.emg_col, forc_lwrad=forc.lwrad, forc_pbot=forc.pbot,
        t_grnd=tp.t_grnd_col, forc_q=forc.q, qg=S.waterdiagbulk.qg_col, forc_th=forc.th)
    lwq_g = C.CfLwqG(; forc_pco2=forc.pco2, forc_po2=forc.po2, forc_u=forc.u_grc, forc_v=forc.v_grc)
    let be = C._kernel_backend(sc.air)
        C._cf_longwave_qsat_kernel!(be)(lwq_out, lwq_p, lwq_c, lwq_g, sc.nmozsgn,
            filterp, pd.column, pd.gridcell, convert(FT, PARAMS.wind_min); ndrange = fn)
        C.KA.synchronize(be)
    end
    C.cf_moninobukini_update!(fv, tp, pd, filterp, fn, sc.ur, sc.dthv, sc.zldis, sc.tl_ini, sc.ts_ini)
    return nothing
end

# one Newton step (energy-balance: friction -> resist -> energy). active≡true.
function newton_step!(S, forc, pft, sc, filterp, fn, active, iter)
    fv = S.frictionvel
    C.friction_velocity!(fv, fn, filterp[1:fn],
        S.canopystate.displa_patch, fv.z0mv_patch, fv.z0hv_patch, fv.z0qv_patch,
        fv.obu_patch, iter + 1, sc.ur, fv.um_patch, fv.ustar_patch,
        sc.temp1, sc.temp2, sc.temp12m, sc.temp22m, sc.fm; active = active)
    C.cf_resist_update!(fv, S.canopystate, S.temperature, S.waterdiagbulk, S.patch_data,
        filterp, fn, active, sc.temp1, sc.temp2, sc.tlbef, sc.del2, sc.del_arr, sc.rah,
        sc.raw, sc.uuc, sc.rb, sc.svpts, sc.eah, sc.el, pft.dleaf, pft.grnd_ch4, forc.pbot,
        PARAMS.csoilc, false, false, false, PARAMS)
    C.cf_energy_update!(S.canopystate, S.energyflux, fv, S.temperature, S.solarabs,
        S.soilstate, S.waterfluxbulk, S.waterstatebulk, S.waterdiagbulk, S.photosyns,
        S.patch_data, S.col_data, filterp, fn, active, sc.rah, sc.raw, sc.rb, sc.rstem,
        sc.sa_leaf, sc.sa_stem, sc.sa_internal, sc.frac_rad_abs_by_stem, sc.air, sc.bir,
        sc.cir, sc.cp_leaf, sc.tl_ini, sc.tlbef, sc.zldis, sc.temp1, sc.temp2, sc.ur,
        sc.efeb, sc.wtg, sc.wtl0, sc.wta0, sc.wtstem0, sc.wtga, sc.wtal, sc.lw_stem,
        sc.lw_leaf, sc.wtgq, sc.wtlq0, sc.wtaq0, sc.wtalq, sc.efe, sc.dt_veg, sc.del_arr,
        sc.err_arr, sc.qsatl, sc.el, sc.qsatldT, sc.dth, sc.dqh, sc.delq, sc.obuold,
        sc.nmozsgn, forc.q, forc.th, forc.pbot, forc.rho,
        SOILEVAP_BETA, SOIL_RESIS_SL14, false, false, NLEVSNO, S.__dtime, PARAMS)
    return nothing
end

# full mirror forward: init + N Newton iterations (active≡true throughout).
function mirror_forward!(B, N)
    S = B.S; forc = B.forc; pft = B.pft
    Snd = merge(S, (; __dtime = B.dtime))   # thread dtime through newton_step!
    sc = mkscr()
    filterp = Int[1]; fn = 1; active = trues(NP)
    canopy_init!(Snd, forc, pft, sc, filterp, fn)
    for it in 0:(N-1)
        newton_step!(Snd, forc, pft, sc, filterp, fn, active, it)
    end
    return nothing
end

# ---- run reference + mirror, compare ----
Bref = build()
run_real!(Bref.S, Bref.forc, Bref.pft, Bref.dtime)
nref = Int(Bref.S.frictionvel.num_iter_patch[1])
@printf("REAL canopy_fluxes!:  t_veg=%.10f  ustar=%.10f  ram1=%.10f  num_iter=%d\n",
    Bref.S.temperature.t_veg_patch[1], Bref.S.frictionvel.ustar_patch[1],
    Bref.S.frictionvel.ram1_patch[1], nref)

for N in (nref, nref + 5, max(nref, 40))
    B = build()
    mirror_forward!(B, N)
    dtveg = abs(B.S.temperature.t_veg_patch[1] - Bref.S.temperature.t_veg_patch[1])
    @printf("MIRROR N=%-3d:         t_veg=%.10f  ustar=%.10f  ram1=%.10f  |Δt_veg|=%.3e\n",
        N, B.S.temperature.t_veg_patch[1], B.S.frictionvel.ustar_patch[1],
        B.S.frictionvel.ram1_patch[1], dtveg)
end

# ===========================================================================
# Stage 2: COMPOSITIONAL REVERSE over init + N Newton iterations, vs FD.
#   input  = temperature.t_grnd_col[1]  (boundary temp, read-only in canopy)
#   output L = sum(abs2, t_veg)
# ===========================================================================
using Enzyme
const N = nref          # fixed iteration count (== the converged count)

# --- split the state into the differentiated bundle + captured Const aux ---
B0 = build()
const PATCH = B0.S.patch_data
const COL   = B0.S.col_data
const GRID  = B0.S.gridcell_data
const FORC  = B0.forc
const PFT   = B0.pft
const FILTERP = Int[1]
const FN = 1
const ACTIVE = trues(NP)
const DTIME = B0.dtime

make_bundle(S) = (; S.canopystate, S.energyflux, S.frictionvel, S.temperature,
    S.solarabs, S.soilstate, S.waterfluxbulk, S.waterstatebulk, S.waterdiagbulk,
    S.photosyns, scratch = mkscr())

# Reconstruct the S-view (differentiated structs from bundle + Const indices/dtime).
Sview(b) = (; b.canopystate, b.energyflux, b.frictionvel, b.temperature, b.solarabs,
    b.soilstate, b.waterfluxbulk, b.waterstatebulk, b.waterdiagbulk, b.photosyns,
    patch_data = PATCH, col_data = COL, gridcell_data = GRID, __dtime = DTIME)

# --- phases as concrete callables (each → its own Enzyme thunk) ---
ph_init!(b)     = (canopy_init!(Sview(b), FORC, PFT, b.scratch, FILTERP, FN); nothing)
ph_resist!(b)   = (fv = b.frictionvel;
    C.cf_resist_update!(fv, b.canopystate, b.temperature, b.waterdiagbulk, PATCH,
        FILTERP, FN, ACTIVE, b.scratch.temp1, b.scratch.temp2, b.scratch.tlbef,
        b.scratch.del2, b.scratch.del_arr, b.scratch.rah, b.scratch.raw, b.scratch.uuc,
        b.scratch.rb, b.scratch.svpts, b.scratch.eah, b.scratch.el, PFT.dleaf,
        PFT.grnd_ch4, FORC.pbot, PARAMS.csoilc, false, false, false, PARAMS); nothing)
ph_energy!(b)   = (newton_energy!(Sview(b), FORC, b.scratch); nothing)
ph_friction!(b, it) = (fv = b.frictionvel;
    C.friction_velocity!(fv, FN, FILTERP[1:FN], b.canopystate.displa_patch,
        fv.z0mv_patch, fv.z0hv_patch, fv.z0qv_patch, fv.obu_patch, it + 1,
        b.scratch.ur, fv.um_patch, fv.ustar_patch, b.scratch.temp1, b.scratch.temp2,
        b.scratch.temp12m, b.scratch.temp22m, b.scratch.fm; active = ACTIVE); nothing)

# energy sub-phase split out (so it can be a standalone thunk)
function newton_energy!(S, forc, sc)
    C.cf_energy_update!(S.canopystate, S.energyflux, S.frictionvel, S.temperature,
        S.solarabs, S.soilstate, S.waterfluxbulk, S.waterstatebulk, S.waterdiagbulk,
        S.photosyns, S.patch_data, S.col_data, FILTERP, FN, ACTIVE, sc.rah, sc.raw,
        sc.rb, sc.rstem, sc.sa_leaf, sc.sa_stem, sc.sa_internal, sc.frac_rad_abs_by_stem,
        sc.air, sc.bir, sc.cir, sc.cp_leaf, sc.tl_ini, sc.tlbef, sc.zldis, sc.temp1,
        sc.temp2, sc.ur, sc.efeb, sc.wtg, sc.wtl0, sc.wta0, sc.wtstem0, sc.wtga, sc.wtal,
        sc.lw_stem, sc.lw_leaf, sc.wtgq, sc.wtlq0, sc.wtaq0, sc.wtalq, sc.efe, sc.dt_veg,
        sc.del_arr, sc.err_arr, sc.qsatl, sc.el, sc.qsatldT, sc.dth, sc.dqh, sc.delq,
        sc.obuold, sc.nmozsgn, forc.q, forc.th, forc.pbot, forc.rho, SOILEVAP_BETA,
        SOIL_RESIS_SL14, false, false, NLEVSNO, DTIME, PARAMS)
    return nothing
end

# --- build the forward phase list (concrete callables) ---
phases = Any[ph_init!]
for it in 0:(N-1)
    push!(phases, let i = it; b -> ph_friction!(b, i); end)
    push!(phases, ph_resist!)
    push!(phases, ph_energy!)
end

# --- FD reference: dL/d(t_grnd) of the mirror primal ---
function L_tgrnd(δ)
    B = build(); B.S.temperature.t_grnd_col[1] += δ; mirror_forward!(B, N)
    return sum(abs2, B.S.temperature.t_veg_patch)
end
hg = 1e-2
g_fd = (L_tgrnd(hg) - L_tgrnd(-hg)) / (2hg)
@printf("\nFD   dL/d(t_grnd)     = % .10e\n", g_fd)

# --- compositional reverse ---
Enzyme.API.strictAliasing!(false)
revmode = Enzyme.set_runtime_activity(Enzyme.Reverse)

b = make_bundle(build().S)
checkpoints = Any[]
for ph in phases
    push!(checkpoints, deepcopy(b))
    ph(b)
end
tveg_final = b.temperature.t_veg_patch[1]
@printf("mirror(bundle) t_veg = %.10f\n", tveg_final)

db = Enzyme.make_zero(b)
db.temperature.t_veg_patch .= 2.0 .* b.temperature.t_veg_patch
for k in length(phases):-1:1
    bk = deepcopy(checkpoints[k])
    Enzyme.autodiff(revmode, phases[k], Enzyme.Const, Enzyme.Duplicated(bk, db))
end
g_comp = db.temperature.t_grnd_col[1]
@printf("comp dL/d(t_grnd)     = % .10e\n", g_comp)

println("\n", "="^64)
relerr = abs(g_comp - g_fd) / max(abs(g_fd), 1e-10)
@printf("abs error = %.3e   rel error = %.3e\n", abs(g_comp - g_fd), relerr)
if isfinite(g_comp) && isfinite(g_fd) && relerr < 1e-4
    println("\nWHOLE-CANOPY COMPOSITIONAL REVERSE-AD VALIDATED ✓")
    println("init + $N Newton iterations, $(length(phases)) sub-phase Enzyme calls,")
    println("gradient through the full energy-balance solve matches finite differences.")
else
    println("\nMISMATCH ✗ — investigate")
end
println("="^64)
