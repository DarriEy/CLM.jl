# ==========================================================================
# Ported from: src/fates/biogeophys/FatesPlantRespPhotosynthMod.F90 (2518 lines)
# FATES Batch 16 (Tier F).
#
# Calculates the plant respiration and photosynthetic fluxes for the FATES
# model. This is the FATES leaf-scale Farquhar/Collatz photosynthesis (C3/C4),
# stomatal conductance (Ball-Berry / Medlyn), leaf maintenance respiration
# (Ryan 1991 / Atkin 2017), and the scaling of leaf-layer fluxes up to the
# cohort. This code is similar to, and was originally based off of, the
# (non-FATES) `PhotosynthesisMod` photosynthesis subroutine -- the
# temperature-response functions (ft1_f/fth_f/fth25_f), the iterative ci solve,
# and the Ball-Berry / Medlyn quadratic conductance forms mirror that port.
#
# Public entry:
#   FatesPlantRespPhotosynthDrive(nsites, sites, bc_in, bc_out, dtime)
#
# Procedures (all ported):
#   FatesPlantRespPhotosynthDrive
#   RootLayerNFixation
#   LeafLayerPhotosynthesis
#   LeafHumidityStomaResis
#   ScaleLeafLayerFluxToCohort
#   ft1_f / fth_f / fth25_f               (temperature response functions)
#   UpdateCanopyNCanNRadPresent
#   GetCanopyGasParameters
#   LeafLayerMaintenanceRespiration_Ryan_1991
#   LeafLayerMaintenanceRespiration_Atkin_etal_2017
#   LeafLayerBiophysicalRates
#   lowstorage_maintresp_reduction
#
# Upstream-Fortran quirks preserved (see inline comments):
#   * `preserve_b4b` round-off-preservation branch for two-stream baselines.
#   * The ci iteration exits after at most 5 iterations (niter==5), NOT 10 as
#     the comment in the Fortran misleadingly states.
#   * C4 RuBP-limited aj uses `quant_eff(c3c4_path_index)` -- index 0 (=0.05),
#     because for C4 the path index is 0.
#   * Stem-only (no leaf) layers fall back to a fraction of cuticular
#     conductance (`stem_cuticle_loss_frac`).
#   * Ryan-1991 lmr25top uses the Q10=1.5 base-rate-at-20C -> 25C conversion.
# ==========================================================================

# --- Module-level constants (Fortran module parameters) ---

# maximum stomatal resistance [s/m]
const _frpp_rsmax0 = 2.0e8

# Ratio of H2O/CO2 gas diffusion in stomatal airspace (approximate)
const h2o_co2_stoma_diffuse_ratio = 1.6
# Ratio of H2O/CO2 gas diffusion in the leaf boundary layer (approximate)
const h2o_co2_bl_diffuse_ratio = 1.4

# Constants used to define C3 versus C4 photosynth pathways
const c3_path_index = 1
const c4_path_index = 0

# Constants used to define conductance models
const medlyn_model = 2
const ballberry_model = 1

# Net vs gross assimilation conductance model
const net_assim_model = 1
const gross_assim_model = 2

const _frpp_preserve_b4b = true
const _frpp_debug = false

# =====================================================================
# Temperature-response functions
# =====================================================================

"""
    ft1_f(tl, ha) -> ans

Photosynthesis temperature response (Arrhenius form normalized to 25 C).
`tl` leaf temperature [K], `ha` activation energy [J/mol]. Uses
`rgas_J_K_kmol * 1e-3` (= J/K/mol) verbatim from the Fortran.
"""
function ft1_f(tl::Real, ha::Real)
    return exp(ha / (rgas_J_K_kmol * 1.0e-3 * (t_water_freeze_k_1atm + 25.0)) *
               (1.0 - (t_water_freeze_k_1atm + 25.0) / tl))
end

"""
    fth_f(tl, hd, se, scaleFactor) -> ans

Photosynthesis high-temperature inhibition (Leuning 2002). `hd` deactivation
energy [J/mol], `se` entropy term [J/mol/K], `scaleFactor` normalizes to 1.0 at
25 C.
"""
function fth_f(tl::Real, hd::Real, se::Real, scaleFactor::Real)
    return scaleFactor / (1.0 + exp((-hd + se * tl) / (rgas_J_K_kmol * 1.0e-3 * tl)))
end

"""
    fth25_f(hd, se) -> ans

Scaling factor for the photosynthesis temperature inhibition (so `fth_f` == 1.0
at 25 C).
"""
function fth25_f(hd::Real, se::Real)
    return 1.0 + exp((-hd + se * (t_water_freeze_k_1atm + 25.0)) /
                     (rgas_J_K_kmol * 1.0e-3 * (t_water_freeze_k_1atm + 25.0)))
end

# =====================================================================
# RootLayerNFixation
# =====================================================================

"""
    RootLayerNFixation(t_soil, ft, dtime, fnrt_mr_layer) -> (fnrt_mr_nfix_layer, nfix_layer)

Symbiotic N fixation (Houlton et al 2008 / Fisher et al 2010). `t_soil` [K],
`fnrt_mr_layer` non-fixation fine-root MR [kgC/s]. Returns the added MR surcharge
to pay for fixation [kgC] and the N fixed in this layer [kgN].
"""
function RootLayerNFixation(t_soil::Real, ft::Integer, dtime::Real, fnrt_mr_layer::Real)
    # N fixation parameters from Houlton et al (2008) and Fisher et al (2010)
    s_fix = -6.25   # s parameter from FUN model (Fisher et al 2010)
    a_fix = -3.62   # a parameter from Houlton et al. 2010
    b_fix = 0.27    # b parameter from Houlton et al. 2010
    c_fix = 25.15   # c parameter from Houlton et al. 2010

    # Amount of C spent (as part of MR respiration) on symbiotic fixation [kgC/s]
    fnrt_mr_nfix_layer = fnrt_mr_layer * prt_params.nfix_mresp_scfrac[ft]

    # Unit carbon cost for nitrogen fixation, temperature dependent [kgC/kgN]
    c_cost_nfix = s_fix * (exp(a_fix + b_fix * (t_soil - t_water_freeze_k_1atm) *
                  (1.0 - 0.5 * (t_soil - t_water_freeze_k_1atm) / c_fix)) - 2.0)

    # Time-integrated C spent on fixation in this layer [kgC/plant/layer/tstep]
    c_spent_nfix = fnrt_mr_nfix_layer * dtime

    # Amount of nitrogen fixed [kgN/plant/layer/tstep]
    nfix_layer = c_spent_nfix / c_cost_nfix

    return fnrt_mr_nfix_layer, nfix_layer
end

# =====================================================================
# GetCanopyGasParameters
# =====================================================================

"""
    GetCanopyGasParameters(can_press, can_o2_partialpress, veg_tempk, air_tempk,
                           air_vpress, veg_esat, rb)
        -> (mm_kco2, mm_ko2, co2_cpoint, cf, gb_mol, ceair)

Michaelis-Menten parameters for CO2/O2, the CO2 compensation point, the molar
<-> velocity conversion factor `cf` [umol/m3], boundary-layer conductance and
constrained vapor pressure.
"""
function GetCanopyGasParameters(can_press::Real, can_o2_partialpress::Real,
                                veg_tempk::Real, air_tempk::Real, air_vpress::Real,
                                veg_esat::Real, rb::Real)
    # Intensive values (per mol of air), Bernacchi et al (2001)
    mm_kc25_umol_per_mol = 404.9
    mm_ko25_mmol_per_mol = 278.4
    co2_cpoint_umol_per_mol = 42.75

    # Activation energies, Bernacchi et al (2001/2003)
    kcha = 79430.0   # activation energy for kc (J/mol)
    koha = 36380.0   # activation energy for ko (J/mol)
    cpha = 37830.0   # activation energy for cp (J/mol)

    kc25 = (mm_kc25_umol_per_mol / umol_per_mol) * can_press
    ko25 = (mm_ko25_mmol_per_mol / mmol_per_mol) * can_press
    sco  = 0.5 * 0.209 / (co2_cpoint_umol_per_mol / umol_per_mol)
    cp25 = 0.5 * can_o2_partialpress / sco

    if veg_tempk > 150.0 && veg_tempk < 350.0
        mm_kco2    = kc25 * ft1_f(veg_tempk, kcha)
        mm_ko2     = ko25 * ft1_f(veg_tempk, koha)
        co2_cpoint = cp25 * ft1_f(veg_tempk, cpha)
    else
        mm_kco2    = 1.0
        mm_ko2     = 1.0
        co2_cpoint = 1.0
    end

    # cf is the conversion factor between molar and velocity forms [umol/m3]
    cf = can_press / (rgas_J_K_kmol * air_tempk) * umol_per_kmol
    gb_mol = (1.0 / rb) * cf

    # Constrain eair >= 0.05*esat_tv (hs > 0) and eair <= veg_esat (hs <= 1)
    ceair = min(max(air_vpress, 0.05 * veg_esat), veg_esat)

    return mm_kco2, mm_ko2, co2_cpoint, cf, gb_mol, ceair
end

# =====================================================================
# LeafLayerMaintenanceRespiration_Ryan_1991
# =====================================================================

"""
    LeafLayerMaintenanceRespiration_Ryan_1991(lnc_top, nscaler, ft, veg_tempk) -> lmr

Leaf (dark) maintenance respiration per Ryan (1991) base rate, with C3/C4
temperature response. `lnc_top` leaf N per area at canopy top [gN/m2]. Returns
lmr [umol CO2/m2/s].
"""
function LeafLayerMaintenanceRespiration_Ryan_1991(lnc_top::Real, nscaler::Real,
                                                   ft::Integer, veg_tempk::Real)
    # Parameters
    lmrha = 46390.0     # activation energy for lmr (J/mol)
    lmrhd = 150650.0    # deactivation energy for lmr (J/mol)
    lmrse = 490.0       # entropy term for lmr (J/mol/K)
    lmrc  = 1.15912391  # high-temperature-inhibition scaling factor (25 C = 1.0)

    edpft = EDPftvarcon_inst[]

    # Base rate is at 20C -> adjust to 25C using CN Q10 = 1.5
    lmr25top = edpft.maintresp_leaf_ryan1991_baserate[ft] * (1.5^((25.0 - 20.0) / 10.0))
    lmr25top = lmr25top * lnc_top / (umolC_to_kgC * g_per_kg)

    lmr25 = lmr25top * nscaler

    c3c4_path_index = round(Int, edpft.c3psn[ft])

    if c3c4_path_index == c3_path_index
        # temperature sensitivity of C3 plants
        lmr = lmr25 * ft1_f(veg_tempk, lmrha) * fth_f(veg_tempk, lmrhd, lmrse, lmrc)
    else
        # temperature sensitivity of C4 plants
        lmr = lmr25 * 2.0^((veg_tempk - (t_water_freeze_k_1atm + 25.0)) / 10.0)
        lmr = lmr / (1.0 + exp(1.3 * (veg_tempk - (t_water_freeze_k_1atm + 55.0))))
    end

    return lmr
end

# =====================================================================
# LeafLayerMaintenanceRespiration_Atkin_etal_2017
# =====================================================================

"""
    LeafLayerMaintenanceRespiration_Atkin_etal_2017(lnc_top, cumulative_lai,
        vcmax25top, ft, veg_tempk, tgrowth) -> lmr

Leaf dark respiration per Atkin et al (2017). `tgrowth` is the lagged veg
temperature averaged over the acclimation timescale [K].
"""
function LeafLayerMaintenanceRespiration_Atkin_etal_2017(lnc_top::Real,
        cumulative_lai::Real, vcmax25top::Real, ft::Integer, veg_tempk::Real,
        tgrowth::Real)
    edpft = EDPftvarcon_inst[]

    # base respiration rate, PFT-dependent (umol CO2/m2/s)
    r_0 = edpft.maintresp_leaf_atkin2017_baserate[ft]

    kn = decay_coeff_vcmax(vcmax25top,
                           edpft.maintresp_leaf_vert_scaler_coeff1[ft],
                           edpft.maintresp_leaf_vert_scaler_coeff2[ft])

    rdark_scaler = exp(-kn * cumulative_lai)

    r_t_ref = max(0.0, rdark_scaler * (r_0 + lmr_r_1 * lnc_top +
                  lmr_r_2 * max(0.0, (tgrowth - t_water_freeze_k_1atm))))

    # Quirk: when r_t_ref == 0 the Fortran emits a "Rdark capped at 0" warning.

    lmr = r_t_ref * exp(lmr_b * (veg_tempk - t_water_freeze_k_1atm - lmr_TrefC) +
              lmr_c * ((veg_tempk - t_water_freeze_k_1atm)^2 - lmr_TrefC^2))

    return lmr
end

# =====================================================================
# LeafLayerBiophysicalRates
# =====================================================================

"""
    LeafLayerBiophysicalRates(parsun_per_la, ft, vcmax25top_ft, jmax25top_ft,
        co2_rcurve_islope25top_ft, nscaler, veg_tempk, dayl_factor, t_growth,
        t_home, btran) -> (vcmax, jmax, co2_rcurve_islope)

Localized (PFT- and leaf-layer-specific) maximum carboxylation rate (vcmax),
maximum electron transport rate (jmax) and the C4 initial slope of the CO2
response curve, factoring temperature, daylength, leaf-N scaling and btran.
"""
function LeafLayerBiophysicalRates(parsun_per_la::Real, ft::Integer,
        vcmax25top_ft::Real, jmax25top_ft::Real, co2_rcurve_islope25top_ft::Real,
        nscaler::Real, veg_tempk::Real, dayl_factor::Real, t_growth::Real,
        t_home::Real, btran::Real)
    edpft = EDPftvarcon_inst[]
    p = ed_params()

    jvr = NaN  # ratio of Jmax25/Vcmax25 (only used by the Kumarathunge model)

    if p.photo_tempsens_model == photosynth_acclim_model_none
        vcmaxha = edpft.vcmaxha[ft]
        jmaxha  = edpft.jmaxha[ft]
        vcmaxhd = edpft.vcmaxhd[ft]
        jmaxhd  = edpft.jmaxhd[ft]
        vcmaxse = edpft.vcmaxse[ft]
        jmaxse  = edpft.jmaxse[ft]
    elseif p.photo_tempsens_model == photosynth_acclim_model_kumarathunge_etal_2019
        t_growth_celsius = t_growth - t_water_freeze_k_1atm
        t_home_celsius   = t_home - t_water_freeze_k_1atm
        vcmaxha = (42.6 + (1.14 * t_growth_celsius)) * 1.0e3  # J/mol
        jmaxha  = 40.71 * 1.0e3                               # J/mol
        vcmaxhd = 200.0 * 1.0e3                               # J/mol
        jmaxhd  = 200.0 * 1.0e3                               # J/mol
        vcmaxse = (645.13 - (0.38 * t_growth_celsius))
        jmaxse  = 658.77 - (0.84 * t_home_celsius) - 0.52 * (t_growth_celsius - t_home_celsius)
        jvr     = 2.56 - (0.0375 * t_home_celsius) - (0.0202 * (t_growth_celsius - t_home_celsius))
    else
        fates_endrun("error, incorrect leaf photosynthesis temperature acclimation model specified")
    end

    vcmaxc = fth25_f(vcmaxhd, vcmaxse)
    jmaxc  = fth25_f(jmaxhd, jmaxse)

    if parsun_per_la <= 0.0
        vcmax = 0.0
        jmax  = 0.0
        co2_rcurve_islope = 0.0
    else  # day time
        # daylength factor local
        dayl_factor_local = (p.dayl_switch == itrue) ? dayl_factor : 1.0

        vcmax25 = vcmax25top_ft * nscaler * dayl_factor_local
        if p.photo_tempsens_model == photosynth_acclim_model_none
            jmax25 = jmax25top_ft * nscaler * dayl_factor_local
        else  # kumarathunge
            jmax25 = vcmax25 * jvr
        end

        co2_rcurve_islope25 = co2_rcurve_islope25top_ft * nscaler

        c3c4_path_index = round(Int, edpft.c3psn[ft])

        if c3c4_path_index == c3_path_index
            vcmax = vcmax25 * ft1_f(veg_tempk, vcmaxha) * fth_f(veg_tempk, vcmaxhd, vcmaxse, vcmaxc)
        else
            vcmax = vcmax25 * 2.0^((veg_tempk - (t_water_freeze_k_1atm + 25.0)) / 10.0)
            vcmax = vcmax / (1.0 + exp(0.2 * ((t_water_freeze_k_1atm + 15.0) - veg_tempk)))
            vcmax = vcmax / (1.0 + exp(0.3 * (veg_tempk - (t_water_freeze_k_1atm + 40.0))))
        end

        jmax = jmax25 * ft1_f(veg_tempk, jmaxha) * fth_f(veg_tempk, jmaxhd, jmaxse, jmaxc)

        # q10 response of product-limited psn
        co2_rcurve_islope = co2_rcurve_islope25 * 2.0^((veg_tempk - (t_water_freeze_k_1atm + 25.0)) / 10.0)
    end

    # Adjust for water limitations
    vcmax = vcmax * btran

    return vcmax, jmax, co2_rcurve_islope
end

# =====================================================================
# lowstorage_maintresp_reduction
# =====================================================================

"""
    lowstorage_maintresp_reduction(frac, pft) -> maintresp_reduction_factor

Reduces maintenance respiration when the storage pool is low (curve in [0,1]).
`frac` is the storage:target-leaf-biomass ratio.
"""
function lowstorage_maintresp_reduction(frac::Real, pft::Integer)
    edpft = EDPftvarcon_inst[]
    if frac < 1.0
        if abs(edpft.maintresp_reduction_curvature[pft] - 1.0) > nearzero
            maintresp_reduction_factor = (1.0 - edpft.maintresp_reduction_intercept[pft]) +
                edpft.maintresp_reduction_intercept[pft] *
                (1.0 - edpft.maintresp_reduction_curvature[pft]^frac) /
                (1.0 - edpft.maintresp_reduction_curvature[pft])
        else  # avoid nan answer for linear case
            maintresp_reduction_factor = (1.0 - edpft.maintresp_reduction_intercept[pft]) +
                edpft.maintresp_reduction_intercept[pft] * frac
        end
    else
        maintresp_reduction_factor = 1.0
    end
    return maintresp_reduction_factor
end

# =====================================================================
# LeafHumidityStomaResis (only meaningful with plant hydraulics)
# =====================================================================

"""
    LeafHumidityStomaResis(leaf_psi, veg_tempk, ceair, can_press, veg_esat, rb,
                           gstoma, ft) -> rstoma_out

Inner-leaf humidity as a function of mesophyll water potential (Vesala et al
2017). Returns total stomatal+boundary-layer resistance [s/m].
"""
function LeafHumidityStomaResis(leaf_psi::Real, veg_tempk::Real, ceair::Real,
                                can_press::Real, veg_esat::Real, rb::Real,
                                gstoma::Real, ft::Integer)
    edpft = EDPftvarcon_inst[]

    # Note: to disable this control, set k_lwp to zero -> lwp_star = 1
    k_lwp = edpft.hydr_k_lwp[ft]
    if leaf_psi < 0.0
        lwp_star = exp(k_lwp * leaf_psi * molar_mass_water / (rgas_J_K_mol * veg_tempk))
    else
        lwp_star = 1.0
    end

    qs   = molar_mass_ratio_vapdry * ceair / (can_press - (1.0 - molar_mass_ratio_vapdry) * ceair)
    qsat = molar_mass_ratio_vapdry * veg_esat / (can_press - (1.0 - molar_mass_ratio_vapdry) * veg_esat)
    qsat_adj = qsat * lwp_star

    if qsat_adj < qs
        # inner leaf vapor pressure <= leaf surface: shut transpiration off
        rstoma_out = _frpp_rsmax0
    else
        rstoma_out = (qsat - qs) * (1.0 / gstoma + rb) / (qsat_adj - qs) - rb
    end

    if rstoma_out < nearzero
        fates_endrun("LeafHumidityStomaResis: negative stomatal resistance")
    end

    return rstoma_out
end

# =====================================================================
# LeafLayerPhotosynthesis
# =====================================================================

"""
    LeafLayerPhotosynthesis(f_sun_lsl, parsun_lsl, parsha_lsl, laisun_lsl,
        laisha_lsl, canopy_area_lsl, ft, vcmax, jmax, co2_rcurve_islope,
        veg_tempk, veg_esat, can_press, can_co2_ppress, can_o2_ppress, btran,
        stomatal_intercept_btran, cf, gb_mol, ceair, mm_kco2, mm_ko2, co2_cpoint,
        lmr, leaf_psi, rb) -> (psn_out, rstoma_out, anet_av_out, c13disc_z)

Leaf-sublayer photosynthesis + stomatal conductance. Iteratively solves the
intercellular CO2 (ci) consistent with co-limited (Rubisco / RuBP / export)
Farquhar-Collatz assimilation and a Ball-Berry or Medlyn stomatal model.
Returns gross psn [umolC/m2/s], stomatal resistance [s/m], net assimilation
averaged over sun/shade, and the C13 discrimination.
"""
function LeafLayerPhotosynthesis(f_sun_lsl::Real, parsun_lsl::Real, parsha_lsl::Real,
        laisun_lsl::Real, laisha_lsl::Real, canopy_area_lsl::Real, ft::Integer,
        vcmax::Real, jmax::Real, co2_rcurve_islope::Real, veg_tempk::Real,
        veg_esat::Real, can_press::Real, can_co2_ppress::Real, can_o2_ppress::Real,
        btran::Real, stomatal_intercept_btran::Real, cf::Real, gb_mol::Real,
        ceair::Real, mm_kco2::Real, mm_ko2::Real, co2_cpoint::Real, lmr::Real,
        leaf_psi::Real, rb::Real)

    p = ed_params()
    edpft = EDPftvarcon_inst[]
    bb_slope     = edpft.bb_slope
    medlyn_slope = edpft.medlyn_slope
    stomatal_intercept = edpft.stomatal_intercept

    theta_cj_c3 = p.theta_cj_c3
    theta_cj_c4 = p.theta_cj_c4

    # Parameters
    fnps = 0.15            # fraction of light absorbed by non-photosynthetic pigments
    photon_to_e = 0.5      # two photons -> one electron in PSII
    wm2_to_umolm2s = 4.6   # W/m2 -> umol photons m-2 s-1
    stem_cuticle_loss_frac = 0.1
    theta_psii = 0.7
    init_a2l_co2_c3 = 0.7
    init_a2l_co2_c4 = 0.4
    quant_eff = (0.05, 0.0)  # C4 (index 0) then C3 (index 1); access via path index + 1
    theta_ip = 0.999
    min_la_to_solve = 1.0e-10

    # photosynthetic pathway: 0. = c4, 1. = c3
    c3c4_path_index = round(Int, edpft.c3psn[ft])

    init_co2_inter_c = (c3c4_path_index == c3_path_index) ?
                       init_a2l_co2_c3 * can_co2_ppress : init_a2l_co2_c4 * can_co2_ppress

    # outputs
    psn_out = 0.0
    rstoma_out = 0.0
    anet_av_out = 0.0
    c13disc_z = 0.0

    if parsun_lsl <= 0.0   # night time
        anet_av_out = -lmr
        psn_out = 0.0
        # cuticular conductance already bounded by max resistance
        rstoma_out = cf / stomatal_intercept_btran
        c13disc_z = 0.0
        return psn_out, rstoma_out, anet_av_out, c13disc_z
    end

    # day time
    if laisun_lsl + laisha_lsl > 0.0
        psn_out = 0.0
        rstoma_out = 0.0
        anet_av_out = 0.0
        gstoma = 0.0

        # locals carried across sunsha (overwritten each iter, used post-loop)
        gs_mol = 0.0
        co2_inter_c = init_co2_inter_c
        anet = 0.0

        for sunsha in 1:2
            # Electron transport rate, convert PAR W/m2 -> umol photons/m2/s
            if sunsha == 1  # sunlit
                if (laisun_lsl * canopy_area_lsl) > min_la_to_solve
                    qabs = parsun_lsl / (laisun_lsl * canopy_area_lsl)
                    qabs = qabs * photon_to_e * (1.0 - fnps) * wm2_to_umolm2s
                else
                    qabs = 0.0
                end
            else  # shaded
                if (parsha_lsl > nearzero) && (laisha_lsl * canopy_area_lsl) > min_la_to_solve
                    qabs = parsha_lsl / (laisha_lsl * canopy_area_lsl)
                    qabs = qabs * photon_to_e * (1.0 - fnps) * wm2_to_umolm2s
                else
                    qabs = 0.0
                end
            end

            aquad = theta_psii
            bquad = -(qabs + jmax)
            cquad = qabs * jmax
            r1, r2 = QuadraticRootsSridharachary(aquad, bquad, cquad)
            je = min(r1, r2)

            # Initialize intercellular co2
            co2_inter_c = init_co2_inter_c

            niter = 0
            loop_continue = true
            while loop_continue
                niter += 1
                co2_inter_c_old = co2_inter_c

                if c3c4_path_index == c3_path_index
                    # C3: Rubisco-limited
                    ac = vcmax * max(co2_inter_c - co2_cpoint, 0.0) /
                         (co2_inter_c + mm_kco2 * (1.0 + can_o2_ppress / mm_ko2))
                    # C3: RuBP-limited
                    aj = je * max(co2_inter_c - co2_cpoint, 0.0) /
                         (4.0 * co2_inter_c + 8.0 * co2_cpoint)
                    # Co-limit ac and aj
                    aquad = theta_cj_c3
                    bquad = -(ac + aj)
                    cquad = ac * aj
                    r1, r2 = QuadraticRootsSridharachary(aquad, bquad, cquad)
                    agross = min(r1, r2)
                else
                    # C4: Rubisco-limited
                    ac = vcmax
                    # C4: RuBP-limited
                    if sunsha == 1  # sunlit
                        if (laisun_lsl * canopy_area_lsl) > 1.0e-10
                            aj = quant_eff[c3c4_path_index + 1] * parsun_lsl * wm2_to_umolm2s
                            aj = aj / (laisun_lsl * canopy_area_lsl)
                        else
                            aj = 0.0
                        end
                    else
                        # Quirk: the C4 SHADED RuBP term has NO div-by-zero
                        # guard (unlike the sunlit branch above) -- preserved
                        # verbatim. A leaf layer with laisha_lsl == 0 would yield
                        # NaN; in practice both sun and shade leaf area are > 0.
                        aj = quant_eff[c3c4_path_index + 1] * parsha_lsl * wm2_to_umolm2s
                        aj = aj / (laisha_lsl * canopy_area_lsl)
                    end
                    # C4: PEP carboxylase-limited (CO2-limited)
                    ap = co2_rcurve_islope * max(co2_inter_c, 0.0) / can_press

                    aquad = theta_cj_c4
                    bquad = -(ac + aj)
                    cquad = ac * aj
                    r1, r2 = QuadraticRootsSridharachary(aquad, bquad, cquad)
                    ai = min(r1, r2)

                    aquad = theta_ip
                    bquad = -(ai + ap)
                    cquad = ai * ap
                    r1, r2 = QuadraticRootsSridharachary(aquad, bquad, cquad)
                    agross = min(r1, r2)
                end

                anet = agross - lmr

                if p.stomatal_assim_model == gross_assim_model
                    if p.stomatal_model == medlyn_model
                        fates_endrun("Gross Assimilation conductance is incompatible with the Medlyn model")
                    end
                    a_gs = agross
                else
                    if anet < 0.0
                        loop_continue = false
                    end
                    a_gs = anet
                end

                # With an <= 0, gs_mol = stomatal_intercept_btran
                leaf_co2_ppress = can_co2_ppress - h2o_co2_bl_diffuse_ratio / gb_mol * a_gs * can_press
                leaf_co2_ppress = max(leaf_co2_ppress, 1.0e-06)

                if p.stomatal_model == medlyn_model
                    # Medlyn et al. (2011), numerics adapted from CLM5
                    vpd = max((veg_esat - ceair), 50.0) * 0.001  # KPa
                    term = h2o_co2_stoma_diffuse_ratio * anet / (leaf_co2_ppress / can_press)
                    aquad = 1.0
                    bquad = -(2.0 * (stomatal_intercept_btran + term) +
                              (medlyn_slope[ft] * term)^2 / (gb_mol * vpd))
                    cquad = stomatal_intercept_btran * stomatal_intercept_btran +
                            (2.0 * stomatal_intercept_btran + term *
                             (1.0 - medlyn_slope[ft] * medlyn_slope[ft] / vpd)) * term
                    r1, r2 = QuadraticRootsSridharachary(aquad, bquad, cquad)
                    gs_mol = max(r1, r2)
                elseif p.stomatal_model == ballberry_model
                    # Ball et al. (1987)
                    aquad = leaf_co2_ppress
                    bquad = leaf_co2_ppress * (gb_mol - stomatal_intercept_btran) -
                            bb_slope[ft] * a_gs * can_press
                    cquad = -gb_mol * (leaf_co2_ppress * stomatal_intercept_btran +
                            bb_slope[ft] * anet * can_press * ceair / veg_esat)
                    r1, r2 = QuadraticRootsSridharachary(aquad, bquad, cquad)
                    gs_mol = max(r1, r2)
                end

                # New estimate for co2_inter_c
                co2_inter_c = can_co2_ppress - anet * can_press *
                    (h2o_co2_bl_diffuse_ratio * gs_mol + h2o_co2_stoma_diffuse_ratio * gb_mol) /
                    (gb_mol * gs_mol)

                # Quirk: exit after convergence OR niter==5 (comment says 10).
                if (abs(co2_inter_c - co2_inter_c_old) / can_press * 1.0e06 <= 2.0e-06) || niter == 5
                    loop_continue = false
                end
            end  # iter_loop

            # End ci iteration. an < 0 -> gs_mol = bbb
            if anet < 0.0
                gs_mol = stomatal_intercept_btran
            end

            # Final estimates
            leaf_co2_ppress = can_co2_ppress - h2o_co2_bl_diffuse_ratio / gb_mol * anet * can_press
            leaf_co2_ppress = max(leaf_co2_ppress, 1.0e-06)
            co2_inter_c = can_co2_ppress - anet * can_press *
                (h2o_co2_bl_diffuse_ratio * gs_mol + h2o_co2_stoma_diffuse_ratio * gb_mol) /
                (gb_mol * gs_mol)

            # gs_mol (umol/m2/s) -> gs (m/s)
            gs = gs_mol / cf

            # C13 discrimination (Ubierna & Farquhar 2014, simplified). b=27, alpha_s=4.4
            c13disc_z = 4.4 + (27.0 - 4.4) *
                        min(can_co2_ppress, max(co2_inter_c, 0.0)) / can_co2_ppress

            # agross = anet + lmr (anet holds the final agross - lmr); identical b4b
            agross_local = anet + lmr
            if sunsha == 1  # sunlit
                psn_out     += agross_local * f_sun_lsl
                anet_av_out += anet * f_sun_lsl
                gstoma      += 1.0 / (min(1.0 / gs, _frpp_rsmax0)) * f_sun_lsl
            else
                psn_out     += agross_local * (1.0 - f_sun_lsl)
                anet_av_out += anet * (1.0 - f_sun_lsl)
                gstoma      += 1.0 / (min(1.0 / gs, _frpp_rsmax0)) * (1.0 - f_sun_lsl)
            end

            if gs_mol < 0.0
                fates_endrun("Negative stomatal conductance")
            end
        end  # sunsha loop

        # Stomatal resistance of the leaf-layer
        if hlm_use_planthydro[] == itrue && edpft.hydr_k_lwp[ft] > nearzero
            rstoma_out = LeafHumidityStomaResis(leaf_psi, veg_tempk, ceair, can_press,
                                                veg_esat, rb, gstoma, ft)
        else
            rstoma_out = 1.0 / gstoma
        end

    else
        # No leaf area: layer present only because of stems.
        psn_out = 0.0
        anet_av_out = 0.0
        rstoma_out = min(_frpp_rsmax0, cf / (stem_cuticle_loss_frac * stomatal_intercept[ft]))
        c13disc_z = 0.0
    end

    return psn_out, rstoma_out, anet_av_out, c13disc_z
end

# =====================================================================
# ScaleLeafLayerFluxToCohort
# =====================================================================

"""
    ScaleLeafLayerFluxToCohort(nv, psn_llz, lmr_llz, rs_llz, elai_llz,
        c13disc_llz, c_area, nplant, rb, maintresp_reduction_factor)
        -> (g_sb_laweight, gpp, rdark, c13disc_clm, cohort_eleaf_area)

Integrate leaf-layer carbon fluxes over the cohort's leaf layers. `gpp`/`rdark`
returned as [kgC/plant/s] (rate; the caller multiplies by dt). `g_sb_laweight`
is NOT normalized -- units [m/s]*[m2 effective leaf].
"""
function ScaleLeafLayerFluxToCohort(nv::Integer, psn_llz::AbstractVector,
        lmr_llz::AbstractVector, rs_llz::AbstractVector, elai_llz::AbstractVector,
        c13disc_llz::AbstractVector, c_area::Real, nplant::Real, rb::Real,
        maintresp_reduction_factor::Real)
    cohort_eleaf_area = 0.0
    g_sb_laweight = 0.0
    gpp = 0.0
    rdark = 0.0
    c13disc_clm = 0.0

    for il in 1:nv
        # Cohort's effective leaf area in this layer [m2]
        cohort_layer_eleaf_area = elai_llz[il] * c_area
        cohort_eleaf_area += cohort_layer_eleaf_area

        # Leaf conductance (stomatal + boundary layer) weighted by leaf area
        g_sb_laweight += 1.0 / (rs_llz[il] + rb) * cohort_layer_eleaf_area

        # GPP [umolC/m2leaf/s] * [m2 leaf] -> [umolC/s]
        gpp += psn_llz[il] * cohort_layer_eleaf_area

        # Dark respiration
        rdark += lmr_llz[il] * cohort_layer_eleaf_area
    end

    if nv > 1
        # weighted-mean d13c flux over leaf layers (1:nv-1, per Fortran)
        sum_weight = sum(psn_llz[1:nv-1] .* elai_llz[1:nv-1])
        if sum_weight == 0.0
            c13disc_clm = 0.0
        else
            c13disc_clm = sum(c13disc_llz[1:nv-1] .* psn_llz[1:nv-1] .* elai_llz[1:nv-1]) / sum_weight
        end
    end

    # Convert dark respiration and GPP from [umol/s] to [kgC/plant/s];
    # apply maintenance respiration reduction to dark resp.
    rdark = rdark * umolC_to_kgC * maintresp_reduction_factor / nplant
    gpp   = gpp * umolC_to_kgC / nplant

    return g_sb_laweight, gpp, rdark, c13disc_clm, cohort_eleaf_area
end

# =====================================================================
# UpdateCanopyNCanNRadPresent
# =====================================================================

"""
    UpdateCanopyNCanNRadPresent(currentPatch)

Fills `currentPatch.nleaf`, `currentPatch.nrad` (== nleaf) and
`currentPatch.canopy_mask` from the cohort list and canopy-area profile.
"""
function UpdateCanopyNCanNRadPresent(currentPatch)
    fill!(currentPatch.nleaf, 0)

    currentCohort = currentPatch.tallest
    while currentCohort !== nothing
        currentPatch.nleaf[currentCohort.canopy_layer, currentCohort.pft] =
            max(currentPatch.nleaf[currentCohort.canopy_layer, currentCohort.pft],
                currentCohort.nv)
        currentCohort = currentCohort.shorter
    end

    # NRAD = NCAN
    currentPatch.nrad .= currentPatch.nleaf

    for cl in 1:currentPatch.ncl_p
        for ft in 1:numpft[]
            currentPatch.canopy_mask[cl, ft] = 0
            for iv in 1:currentPatch.nrad[cl, ft]
                if currentPatch.canopy_area_profile[cl, ft, iv] > 0.0
                    currentPatch.canopy_mask[cl, ft] = 1
                end
            end
        end
    end

    return nothing
end

# =====================================================================
# FatesPlantRespPhotosynthDrive  (public entry)
# =====================================================================

"""
    FatesPlantRespPhotosynthDrive(nsites, sites, bc_in, bc_out, dtime)

Leaf photosynthesis, stomatal conductance and whole-plant maintenance/growth
respiration over all sites/patches/cohorts/leaf-layers. Mutates cohort flux
accumulators and `bc_out` canopy resistances, mirroring the Fortran driver.

NOTE: the two-stream radiation path and plant-hydraulics path are gated off by
default (`radiation_model == norman_solver`, `hlm_use_planthydro == 0`). The
two-stream branch is ported (calls `FatesGetCohortAbsRad`), as is the
plant-hydraulics conductance branch.
"""
function FatesPlantRespPhotosynthDrive(nsites::Integer, sites, bc_in, bc_out, dtime::Real)
    p = ed_params()
    edpft = EDPftvarcon_inst[]

    slatop  = prt_params.slatop
    woody   = prt_params.woody
    stomatal_intercept = edpft.stomatal_intercept

    dinc_vai   = p.dinc_vai
    dlower_vai = p.dlower_vai
    radiation_model = p.radiation_model

    # scratch arrays (leaf layer, pft, canopy layer)
    lmr_z       = zeros(Float64, nlevleaf, maxpft, nclmax)
    rs_z        = zeros(Float64, nlevleaf, maxpft, nclmax)
    anet_av_z   = zeros(Float64, nlevleaf, maxpft, nclmax)
    psn_z       = zeros(Float64, nlevleaf, maxpft, nclmax)
    rate_mask_z = falses(nlevleaf, maxpft, nclmax)
    # c13disc_z is indexed (cl, ft, iv) in the Fortran
    c13disc_z   = zeros(Float64, nclmax, maxpft, nlevleaf)

    for s in 1:nsites
        # Pre-process PFT-dependent (but env-independent) root fractions
        rootfr_ft = zeros(Float64, numpft[], bc_in[s].nlevsoil)
        for ft in 1:numpft[]
            set_root_fraction(view(rootfr_ft, ft, :), ft, bc_in[s].zi_sisl;
                              max_nlevroot = bc_in[s].max_rooting_depth_index_col)
        end

        ifp = 0
        currentPatch = sites[s].oldest_patch
        while currentPatch !== nothing
            if currentPatch.nocomp_pft_label != nocomp_bareground
                ifp += 1
                NCL_p = currentPatch.ncl_p

                # Part I. Zero output boundary conditions
                bc_out[s].rssun_pa[ifp] = 0.0
                bc_out[s].rssha_pa[ifp] = 0.0
                psn_z .= 0.0
                g_sb_leaves = 0.0
                patch_la = 0.0

                # Part II. Patch photosynthesis filter
                if bc_in[s].filter_photo_pa[ifp] == 2
                    # Part III.
                    UpdateCanopyNCanNRadPresent(currentPatch)

                    # Part IV. Environmentally derived gas parameters
                    mm_kco2, mm_ko2, co2_cpoint, cf, gb_mol, ceair =
                        GetCanopyGasParameters(bc_in[s].forc_pbot, bc_in[s].oair_pa[ifp],
                            bc_in[s].t_veg_pa[ifp], bc_in[s].tgcm_pa[ifp],
                            bc_in[s].eair_pa[ifp], bc_in[s].esat_tv_pa[ifp],
                            bc_in[s].rb_pa[ifp])

                    # Part VI. Loop over leaf layers via cohorts
                    @views rate_mask_z[:, 1:numpft[], :] .= false

                    if currentPatch.countcohorts > 0
                        currentCohort = currentPatch.tallest
                        while currentCohort !== nothing
                            ft = currentCohort.pft
                            cl = currentCohort.canopy_layer

                            # Cohort-specific elai profile + top/bottom VAI edges
                            cohort_vaitop      = zeros(Float64, 75)
                            cohort_vaibot      = zeros(Float64, 75)
                            cohort_layer_elai  = zeros(Float64, 75)
                            cohort_layer_esai  = zeros(Float64, 75)
                            cohort_elai = 0.0
                            cohort_esai = 0.0
                            if currentCohort.treesai > 0.0
                                for iv in 1:currentCohort.nv
                                    vt, vb, ela, esa, _, _ = VegAreaLayer(
                                        currentCohort.treelai, currentCohort.treesai,
                                        currentCohort.height, iv, currentCohort.nv,
                                        currentCohort.pft, sites[s].snow_depth)
                                    cohort_vaitop[iv]     = vt
                                    cohort_vaibot[iv]     = vb
                                    cohort_layer_elai[iv] = ela
                                    cohort_layer_esai[iv] = esa
                                end
                                cohort_elai = sum(@view cohort_layer_elai[1:currentCohort.nv])
                                cohort_esai = sum(@view cohort_layer_esai[1:currentCohort.nv])
                            end

                            # Maintenance respiration storage reduction
                            store_c_target, _ = bleaf(currentCohort.dbh, currentCohort.pft,
                                currentCohort.crowndamage, currentCohort.canopy_trim, 1.0)
                            frac = storage_fraction_of_target(store_c_target,
                                GetState(currentCohort.prt, store_organ, carbon12_element))
                            maintresp_reduction_factor =
                                lowstorage_maintresp_reduction(frac, currentCohort.pft)

                            # Are there leaves of this pft in this layer?
                            if currentPatch.canopy_mask[cl, ft] == 1
                                for iv in 1:currentCohort.nv
                                    # Recalc leaf biophysical rates if needed
                                    if (!rate_mask_z[iv, ft, cl]) ||
                                       (hlm_use_planthydro[] == itrue) ||
                                       (radiation_model == twostr_solver) ||
                                       (nleafage[] > 1) ||
                                       (hlm_parteh_mode[] != prt_carbon_allom_hyp)

                                        if hlm_use_planthydro[] == itrue
                                            stomatal_intercept_btran = max(cf / _frpp_rsmax0,
                                                stomatal_intercept[ft] * currentCohort.co_hydr.btran)
                                            btran_eff = currentCohort.co_hydr.btran
                                            leaf_inc = dinc_vai[iv] *
                                                currentCohort.treelai / (currentCohort.treelai + currentCohort.treesai)
                                            lai_canopy_above = sum(@view currentPatch.canopy_layer_tlai[1:cl-1])
                                            lai_layers_above = (dlower_vai[iv] - dinc_vai[iv]) *
                                                currentCohort.treelai / (currentCohort.treelai + currentCohort.treesai)
                                            lai_current = min(leaf_inc, currentCohort.treelai - lai_layers_above)
                                            cumulative_lai = lai_canopy_above + lai_layers_above + 0.5 * lai_current
                                            leaf_psi = currentCohort.co_hydr.psi_ag[1]
                                        else
                                            stomatal_intercept_btran = max(cf / _frpp_rsmax0,
                                                stomatal_intercept[ft] * currentPatch.btran_ft[ft])
                                            btran_eff = currentPatch.btran_ft[ft]
                                            cumulative_lai = sum(@view currentPatch.canopy_layer_tlai[1:cl-1]) +
                                                sum(@view currentPatch.tlai_profile[cl, ft, 1:iv-1]) +
                                                0.5 * currentPatch.tlai_profile[cl, ft, iv]
                                            leaf_psi = fates_unset_r8
                                        end

                                        if do_fates_salinity
                                            btran_eff = btran_eff * currentPatch.bstress_sal_ft[ft]
                                        end

                                        # Leaf nitrogen decay coefficient (Lloyd et al 2010)
                                        kn = decay_coeff_vcmax(currentCohort.vcmax25top,
                                            prt_params.leafn_vert_scaler_coeff1[ft],
                                            prt_params.leafn_vert_scaler_coeff2[ft])
                                        nscaler = exp(-kn * cumulative_lai)

                                        # Leaf N content per area at canopy top [gN/m2]
                                        if hlm_parteh_mode[] == prt_carbon_allom_hyp
                                            lnc_top = prt_params.nitr_stoich_p1[ft,
                                                prt_params.organ_param_id[leaf_organ]] / slatop[ft]
                                        elseif hlm_parteh_mode[] == prt_cnp_flex_allom_hyp
                                            leaf_c = GetState(currentCohort.prt, leaf_organ, carbon12_element)
                                            if (leaf_c * slatop[ft]) > nearzero
                                                leaf_n = GetState(currentCohort.prt, leaf_organ, nitrogen_element)
                                                lnc_top = leaf_n / (slatop[ft] * leaf_c)
                                            else
                                                lnc_top = prt_params.nitr_stoich_p1[ft,
                                                    prt_params.organ_param_id[leaf_organ]] / slatop[ft]
                                            end
                                        else
                                            lnc_top = prt_params.nitr_stoich_p1[ft,
                                                prt_params.organ_param_id[leaf_organ]] / slatop[ft]
                                        end

                                        # Part VII: dark respiration (leaf maintenance)
                                        if p.maintresp_leaf_model == lmrmodel_ryan_1991
                                            lmr_z[iv, ft, cl] = LeafLayerMaintenanceRespiration_Ryan_1991(
                                                lnc_top, nscaler, ft, bc_in[s].t_veg_pa[ifp])
                                        elseif p.maintresp_leaf_model == lmrmodel_atkin_etal_2017
                                            lmr_z[iv, ft, cl] = LeafLayerMaintenanceRespiration_Atkin_etal_2017(
                                                lnc_top, cumulative_lai, currentCohort.vcmax25top, ft,
                                                bc_in[s].t_veg_pa[ifp], GetMean(currentPatch.tveg_lpa))
                                        else
                                            fates_endrun("error, incorrect leaf respiration model specified")
                                        end

                                        # PAR per sunlit/shaded leaf area
                                        if radiation_model == norman_solver
                                            laisun = currentPatch.ed_laisun_z[cl, ft, iv]
                                            laisha = currentPatch.ed_laisha_z[cl, ft, iv]
                                            par_per_sunla = currentPatch.ed_parsun_z[cl, ft, iv]
                                            par_per_shala = currentPatch.ed_parsha_z[cl, ft, iv]
                                            canopy_area   = currentPatch.canopy_area_profile[cl, ft, iv]
                                            fsun = currentPatch.f_sun[cl, ft, iv]
                                        else  # Two-stream
                                            if cohort_layer_elai[iv] > nearzero && currentPatch.solar_zenith_flag
                                                rb_abs, rd_abs, rb_abs_leaf, rd_abs_leaf, fsun =
                                                    FatesGetCohortAbsRad(currentPatch, currentCohort, ipar,
                                                        cohort_vaitop[iv], cohort_vaibot[iv], cohort_elai, cohort_esai)
                                                laisun = fsun * cohort_layer_elai[iv]
                                                laisha = (1.0 - fsun) * cohort_layer_elai[iv]
                                                if fsun > nearzero
                                                    par_per_sunla = rd_abs_leaf * fsun + rb_abs_leaf
                                                else
                                                    par_per_sunla = 0.0
                                                end
                                                par_per_shala = rd_abs_leaf * (1.0 - fsun)
                                                canopy_area = 1.0
                                            else
                                                par_per_sunla = 0.0
                                                par_per_shala = 0.0
                                                laisun = 0.5 * cohort_layer_elai[iv]
                                                laisha = 0.5 * cohort_layer_elai[iv]
                                                canopy_area = 1.0
                                                fsun = 0.5  # avoid div0
                                            end
                                        end

                                        # Part VIII: localized biophysical rates
                                        vcmax_z, jmax_z, kp_z = LeafLayerBiophysicalRates(
                                            par_per_sunla, ft, currentCohort.vcmax25top,
                                            currentCohort.jmax25top, currentCohort.kp25top,
                                            nscaler, bc_in[s].t_veg_pa[ifp], bc_in[s].dayl_factor_pa[ifp],
                                            GetMean(currentPatch.tveg_lpa), GetMean(currentPatch.tveg_longterm),
                                            btran_eff)

                                        # Part IX: actual photosynthesis
                                        psn_z[iv, ft, cl], rs_z[iv, ft, cl], anet_av_z[iv, ft, cl],
                                        c13disc_z[cl, ft, iv] = LeafLayerPhotosynthesis(
                                            fsun, par_per_sunla, par_per_shala, laisun, laisha,
                                            canopy_area, ft, vcmax_z, jmax_z, kp_z,
                                            bc_in[s].t_veg_pa[ifp], bc_in[s].esat_tv_pa[ifp],
                                            bc_in[s].forc_pbot, bc_in[s].cair_pa[ifp], bc_in[s].oair_pa[ifp],
                                            btran_eff, stomatal_intercept_btran, cf, gb_mol, ceair,
                                            mm_kco2, mm_ko2, co2_cpoint, lmr_z[iv, ft, cl],
                                            leaf_psi, bc_in[s].rb_pa[ifp])

                                        rate_mask_z[iv, ft, cl] = true
                                    end
                                end  # leaf_layer_loop

                                # Zero cohort flux accumulators
                                currentCohort.npp_tstep  = 0.0
                                currentCohort.resp_tstep = 0.0
                                currentCohort.gpp_tstep  = 0.0
                                currentCohort.rdark      = 0.0
                                currentCohort.resp_m     = 0.0
                                currentCohort.ts_net_uptake .= 0.0
                                currentCohort.c13disc_clm = 0.0

                                nv = currentCohort.nv

                                # Transfer leaf-layer fluxes to cohort
                                if radiation_model == norman_solver
                                    elai_llz = view(currentPatch.elai_profile, cl, ft, 1:nv)
                                else
                                    elai_llz = view(cohort_layer_elai, 1:nv)
                                end
                                g_sb_laweight, gpp, rdark, c13disc_clm, cohort_eleaf_area =
                                    ScaleLeafLayerFluxToCohort(nv,
                                        view(psn_z, 1:nv, ft, cl), view(lmr_z, 1:nv, ft, cl),
                                        view(rs_z, 1:nv, ft, cl), elai_llz,
                                        view(c13disc_z, cl, ft, 1:nv), currentCohort.c_area,
                                        currentCohort.n, bc_in[s].rb_pa[ifp], maintresp_reduction_factor)
                                currentCohort.g_sb_laweight = g_sb_laweight
                                currentCohort.gpp_tstep = gpp
                                currentCohort.rdark = rdark
                                currentCohort.c13disc_clm = c13disc_clm

                                # Net uptake transferred directly
                                @views currentCohort.ts_net_uptake[1:nv] .=
                                    anet_av_z[1:nv, ft, cl] .* umolC_to_kgC
                            else
                                cohort_eleaf_area = 0.0
                                currentCohort.gpp_tstep = 0.0
                                currentCohort.rdark = 0.0
                                currentCohort.g_sb_laweight = 0.0
                                currentCohort.ts_net_uptake .= 0.0
                            end  # canopy_mask_if

                            # Part VIII: sapwood + fine root maintenance respiration
                            sapw_c = GetState(currentCohort.prt, sapw_organ, carbon12_element)
                            fnrt_c = GetState(currentCohort.prt, fnrt_organ, carbon12_element)

                            if hlm_use_tree_damage[] == itrue
                                crown_reduction = GetCrownReduction(currentCohort.crowndamage)
                            else
                                crown_reduction = 0.0
                            end

                            agb_frac = prt_params.allom_agb_frac[currentCohort.pft]
                            branch_frac = param_derived().branch_frac[currentCohort.pft]
                            sapw_c_undamaged = sapw_c / (1.0 - (agb_frac * branch_frac * crown_reduction))
                            sapw_c_bgw = sapw_c_undamaged * (1.0 - agb_frac)
                            sapw_c_agw = sapw_c - sapw_c_bgw

                            if hlm_parteh_mode[] == prt_carbon_allom_hyp
                                live_stem_n = sapw_c_agw * prt_params.nitr_stoich_p1[ft,
                                    prt_params.organ_param_id[sapw_organ]]
                                live_croot_n = sapw_c_bgw * prt_params.nitr_stoich_p1[ft,
                                    prt_params.organ_param_id[sapw_organ]]
                                fnrt_n = fnrt_c * prt_params.nitr_stoich_p1[ft,
                                    prt_params.organ_param_id[fnrt_organ]]
                            elseif hlm_parteh_mode[] == prt_cnp_flex_allom_hyp
                                live_stem_n = prt_params.allom_agb_frac[currentCohort.pft] *
                                    GetState(currentCohort.prt, sapw_organ, nitrogen_element)
                                live_croot_n = (1.0 - prt_params.allom_agb_frac[currentCohort.pft]) *
                                    GetState(currentCohort.prt, sapw_organ, nitrogen_element)
                                fnrt_n = GetState(currentCohort.prt, fnrt_organ, nitrogen_element)
                                if hlm_use_tree_damage[] == itrue
                                    sapw_n = GetState(currentCohort.prt, sapw_organ, nitrogen_element)
                                    sapw_n_undamaged = sapw_n / (1.0 - (agb_frac * branch_frac * crown_reduction))
                                    sapw_n_bgw = sapw_n_undamaged * (1.0 - agb_frac)
                                    sapw_n_agw = sapw_n - sapw_n_bgw
                                    live_croot_n = sapw_n_bgw
                                    live_stem_n = sapw_n_agw
                                end
                            else
                                live_stem_n = 0.0
                                live_croot_n = 0.0
                                fnrt_n = 0.0
                            end

                            # Live stem MR (kgC/plant/s) (above-ground sapwood)
                            if round(Int, woody[ft]) == itrue
                                tcwood = p.q10_mr^((bc_in[s].t_veg_pa[ifp] - t_water_freeze_k_1atm - 20.0) / 10.0)
                                currentCohort.livestem_mr = live_stem_n * p.maintresp_nonleaf_baserate *
                                    tcwood * maintresp_reduction_factor
                            else
                                currentCohort.livestem_mr = 0.0
                            end

                            # Fine root MR + N fixation
                            currentCohort.froot_mr = 0.0
                            currentCohort.sym_nfix_tstep = 0.0
                            for j in 1:bc_in[s].nlevsoil
                                tcsoi = p.q10_mr^((bc_in[s].t_soisno_sl[j] - t_water_freeze_k_1atm - 20.0) / 10.0)
                                fnrt_mr_layer = fnrt_n * p.maintresp_nonleaf_baserate * tcsoi *
                                    rootfr_ft[ft, j] * maintresp_reduction_factor
                                fnrt_mr_nfix_layer, nfix_layer = RootLayerNFixation(
                                    bc_in[s].t_soisno_sl[j], ft, dtime, fnrt_mr_layer)
                                currentCohort.froot_mr += fnrt_mr_nfix_layer + fnrt_mr_layer
                                currentCohort.sym_nfix_tstep += nfix_layer
                            end

                            # Coarse root MR (below-ground sapwood)
                            if round(Int, woody[ft]) == itrue
                                currentCohort.livecroot_mr = 0.0
                                for j in 1:bc_in[s].nlevsoil
                                    tcsoi = p.q10_mr^((bc_in[s].t_soisno_sl[j] - t_water_freeze_k_1atm - 20.0) / 10.0)
                                    currentCohort.livecroot_mr += live_croot_n * p.maintresp_nonleaf_baserate *
                                        tcsoi * rootfr_ft[ft, j] * maintresp_reduction_factor
                                end
                            else
                                currentCohort.livecroot_mr = 0.0
                            end

                            # Part IX: unit conversions + net/sum fluxes
                            currentCohort.resp_m = currentCohort.livestem_mr +
                                currentCohort.livecroot_mr + currentCohort.froot_mr
                            currentCohort.resp_m += currentCohort.rdark
                            currentCohort.resp_m_unreduced = currentCohort.resp_m / maintresp_reduction_factor

                            # kgC/indiv/s -> kgC/indiv/timestep
                            currentCohort.resp_m        *= dtime
                            currentCohort.gpp_tstep     *= dtime
                            currentCohort.ts_net_uptake .*= dtime

                            currentCohort.resp_g_tstep = prt_params.grperc[ft] *
                                max(0.0, currentCohort.gpp_tstep - currentCohort.resp_m)
                            currentCohort.resp_tstep = currentCohort.resp_m + currentCohort.resp_g_tstep
                            currentCohort.npp_tstep  = currentCohort.gpp_tstep - currentCohort.resp_tstep

                            g_sb_leaves += currentCohort.g_sb_laweight
                            patch_la += cohort_eleaf_area

                            currentCohort = currentCohort.shorter
                        end  # do_cohort_drive
                    end  # if_any_cohorts

                    # Normalize canopy total conductance by effective LAI
                    if _frpp_preserve_b4b
                        patch_la = patch_la / currentPatch.total_canopy_area
                    end

                    if patch_la > floatmin(Float64)
                        if _frpp_preserve_b4b
                            elai = calc_areaindex(currentPatch, "elai")
                            g_sb_leaves = g_sb_leaves / (elai * currentPatch.total_canopy_area)
                        else
                            g_sb_leaves = g_sb_leaves / max(0.1 * currentPatch.total_canopy_area, patch_la)
                        end

                        if g_sb_leaves > (1.0 / _frpp_rsmax0)
                            r_sb_leaves = 1.0 / g_sb_leaves
                            if r_sb_leaves < bc_in[s].rb_pa[ifp]
                                fates_endrun("Combined canopy resistance smaller than its boundary layer component")
                            end
                            r_stomata = r_sb_leaves - bc_in[s].rb_pa[ifp]
                        else
                            r_stomata = _frpp_rsmax0
                        end

                        bc_out[s].rssun_pa[ifp] = r_stomata
                        bc_out[s].rssha_pa[ifp] = r_stomata
                        currentPatch.c_stomata = cf / r_stomata
                    else
                        bc_out[s].rssun_pa[ifp] = _frpp_rsmax0
                        bc_out[s].rssha_pa[ifp] = _frpp_rsmax0
                        currentPatch.c_stomata = cf / _frpp_rsmax0
                    end

                    currentPatch.c_lblayer = cf / bc_in[s].rb_pa[ifp]
                end  # if_filter2
            end  # if_notbare

            currentPatch = currentPatch.younger
        end  # patch loop
    end  # site loop

    return nothing
end
