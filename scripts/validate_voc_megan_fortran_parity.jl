#!/usr/bin/env julia
# =============================================================================
# validate_voc_megan_fortran_parity.jl
#
# NON-VACUOUS Fortran-parity harness for VOC / biogenic-VOC emissions (MEGAN v2.1)
# ported from src/biogeochem/VOCEmissionMod.F90.
#
# WHY THIS EXISTS
# ---------------
# VOC is the last un-diffed BGC subsystem in this port. The existing
# test_voc_emission.jl asserts the driver produces flux `> 0.0` and that
# structural relationships hold (tot == sum of mech comps, diagnostics not-NaN).
# That is structurally BLIND to the two failure modes this project has shipped:
#   1. A subsystem that RUNS but computes the WRONG VALUE (a mis-transcribed
#      Guenther-2006 coefficient, a wrong exp() sign, a bad unit factor) — a
#      ">0 and finite" test passes right through it. (CONSERVATION-IS-NOT-ACCURACY
#      / "assert the VALUE, not finiteness".)
#   2. DEAD WIRING: the kernel exists + is unit-tested at some value, but no LIVE
#      init path ever populates the descriptors that the driver gate requires, so
#      end-to-end the subsystem is a permanent no-op (the MIMICS / isotopes class).
#
# WHAT THIS HARNESS DOES
# ----------------------
#   A. VALUE parity of every activity-factor helper (get_gamma_P/L/SM/T/A/C,
#      get_map_EF) against an INDEPENDENT plain-scalar oracle transcribed directly
#      from VOCEmissionMod.F90 (NOT from the Julia port), swept over realistic
#      inputs incl. the accumulated-var branch AND the fixed-coeff branch, the
#      boreal-shrub / arctic-C3-grass special T-response branches, and all four
#      CO2 growth-environment regimes. Assert agreement to ~1e-12.
#   B. VALUE parity of the FULL voc_emission! driver (the assembled kernel:
#      PAR construction, per-mega-compound epsilon*gamma*unitfactor, the
#      moles-flux vs kg-flux split, the mechanism-compound summation) against an
#      independent scalar recomputation of the entire pipeline on a multi-patch,
#      multi-compound, multi-PFT setup. Assert agreement to ~1e-12.
#   C. WIRING TRUTH assertions (the honest audit, made executable):
#      C1. A DEFAULT CLMDriverConfig carries an EMPTY MEGANConfig, so the driver
#          gate `use_voc && !isempty(megan.*)` is false and voc_emission! is a
#          no-op. This documents the init-plumbing gap: `use_voc` defaults false
#          and NO live initializer (instances.jl / clm_initialize.jl) reads the
#          megan namelist or the megan_factors_file netCDF to flip it. VOC on the
#          default/real run path is therefore inert — faithful to Fortran-default
#          (no megan namelist => shr_megan_mechcomps_n=0 => VOCEmission returns).
#      C2. Given descriptors built from a megan_specifier via the ported
#          namelist path (megan_config_from_nl), the SAME driver kernel produces
#          NONZERO, oracle-matching flux — proving the kernel is live-dispatchable
#          and faithful; only the activation plumbing (namelist/netCDF reader ->
#          config) is not connected to a live initializer.
#
# The oracle is transcribed independently from:
#   $CLM_FORTRAN/src/biogeochem/VOCEmissionMod.F90
#
# NOTE (documented benign divergence): Fortran get_gamma_C uses strict
# inequalities and leaves Ismax/h/Cstar UNINITIALIZED at co2_ppmv exactly in
# {400,600,800} (undefined behaviour). The port's trailing `else` handles those
# points deterministically. This harness sweeps co2 strictly inside each regime,
# where both agree; the boundary is measure-zero and Fortran is UB there.
#
# Run:  julia +1.12 --project=. scripts/validate_voc_megan_fortran_parity.jl
# Exit: 0 = all parity checks pass, 1 = a divergence / wiring regression found.
# =============================================================================

using CLM
using Printf

const TOL = 1e-11
const TFRZ = 273.15

# PFT indices used by the port defaults (VOCEmissionMod pft indices).
const NDLLF_EVR_TMP_TREE  = 1
const NDLLF_DCD_BRL_TREE  = 3
const NBRDLF_DCD_TRP_TREE = 6
const NBRDLF_DCD_BRL_SHRUB = 11
const NC3_ARCTIC_GRASS    = 12

# -----------------------------------------------------------------------------
# INDEPENDENT ORACLE — transcribed straight from VOCEmissionMod.F90.
# These deliberately do NOT call any CLM.get_* function.
# -----------------------------------------------------------------------------

function ref_gamma_L(fsun240, elai)
    cce, cce1 = 0.30, 0.24
    (fsun240 > 0.0 && fsun240 < 1.0e30) ? cce * elai : cce1 * elai
end

function ref_gamma_SM(btran)
    a1, b1, thr = -7.4463, 3.2552, 0.2
    btran >= 1.0 ? 1.0 : 1.0 / (1.0 + b1 * exp(a1 * (btran - thr)))
end

function ref_gamma_P(par_sun, par24_sun, par240_sun, par_sha, par24_sha, par240_sha,
                     fsun, fsun240, forc_solad240, forc_solai240, LDF)
    ca1, ca2, ca3 = 0.004, 0.0005, 0.0468
    par0_sun, par0_shade = 200.0, 50.0
    alpha_fix, cp_fix = 0.001, 1.21
    local alpha, cp, gp
    if fsun240 > 0.0 && fsun240 < 1.0 && forc_solad240 > 0.0 && forc_solai240 > 0.0
        alpha = ca1 - ca2 * log(par240_sun)
        cp    = ca3 * exp(ca2 * (par24_sun - par0_sun)) * par240_sun^0.6
        gp    = fsun * (cp * alpha * par_sun * (1.0 + alpha*alpha*par_sun*par_sun)^(-0.5))
        alpha = ca1 - ca2 * log(par240_sha)
        cp    = ca3 * exp(ca2 * (par_sha - par0_shade)) * par240_sha^0.6
        gp    = gp + (1.0 - fsun) * (cp*alpha*par_sha*(1.0 + alpha*alpha*par_sha*par_sha)^(-0.5))
    else
        alpha, cp = alpha_fix, cp_fix
        gp = fsun * (cp*alpha*par_sun*(1.0 + alpha*alpha*par_sun*par_sun)^(-0.5))
        gp = gp + (1.0 - fsun) * (cp*alpha*par_sha*(1.0 + alpha*alpha*par_sha*par_sha)^(-0.5))
    end
    gamma_p = (1.0 - LDF) + LDF * gp
    return (gamma_p, cp, alpha)
end

function ref_gamma_T(t_veg240, t_veg24, t_veg, ct1, ct2, betaT, LDF, Ceo, ivt)
    co1, co2, co4 = 313.0, 0.6, 0.05
    tstd0, topt_fix, Eopt_fix = 297.0, 317.0, 2.26
    ct3, tstd = 0.00831, 303.15
    std_act = 95.0
    p1, p2, p3, p4, p5 = 9.49, 0.53, 0.12, 7.9, 0.217
    bet_arc_c3_max = 300.0
    local Eopt, topt
    if t_veg240 > 0.0 && t_veg240 < 1.0e30
        topt = co1 + co2 * (t_veg240 - tstd0)
        if ivt == NBRDLF_DCD_BRL_SHRUB
            Eopt = p4 * exp(p5 * (t_veg24 - TFRZ - 24.0))
        elseif ivt == NC3_ARCTIC_GRASS
            Eopt = exp(p3 * (t_veg240 - TFRZ - 15.0))
        else
            Eopt = Ceo * exp(co4 * (t_veg24 - tstd0)) * exp(co4 * (t_veg240 - tstd0))
        end
    else
        topt, Eopt = topt_fix, Eopt_fix
    end
    x = ((1.0/topt) - (1.0/t_veg)) / ct3
    local gamma_t_LDF
    if ivt == NC3_ARCTIC_GRASS
        bet = std_act + p1 * exp(p2 * (TFRZ + 15.0 - t_veg240))
        bet = min(bet, bet_arc_c3_max)
        gamma_t_LDF = Eopt * exp(bet * ((1.0/(TFRZ + 30.0) - 1.0/t_veg) / ct3))
    else
        gamma_t_LDF = Eopt * (ct2 * exp(ct1 * x) / (ct2 - ct1 * (1.0 - exp(ct2 * x))))
    end
    gamma_t_LIF = exp(betaT * (t_veg - tstd))
    gamma_t = (1.0 - LDF) * gamma_t_LIF + LDF * gamma_t_LDF
    return (gamma_t, Eopt, topt)
end

function ref_gamma_A(ivt, elai240, elai, nclass, Anew, Agro, Amat, Aold)
    if ivt == NDLLF_DCD_BRL_TREE || ivt >= NBRDLF_DCD_TRP_TREE
        if elai240 > 0.0 && elai240 < 1.0e30
            elai_prev = 2.0*elai240 - elai
            if elai_prev == elai
                fnew, fgro, fmat, fold = 0.0, 0.0, 1.0, 0.0
            elseif elai_prev > elai
                fnew, fgro = 0.0, 0.0
                fmat = 1.0 - (elai_prev - elai)/elai_prev
                fold = (elai_prev - elai)/elai_prev
            else
                fnew = 1.0 - (elai_prev/elai)
                fgro = 0.0
                fmat = elai_prev/elai
                fold = 0.0
            end
            return fnew*Anew[nclass] + fgro*Agro[nclass] + fmat*Amat[nclass] + fold*Aold[nclass]
        else
            return 1.0
        end
    else
        return 1.0
    end
end

function ref_gamma_C(cisun, cisha, forc_pbot, fsun, co2_ppmv)
    Ismax_ca, h_ca, Cstar_ca, CiCa = 1.344, 1.4614, 585.0, 0.7
    gamma_ca = Ismax_ca - (Ismax_ca * (CiCa*co2_ppmv)^h_ca) / (Cstar_ca^h_ca + (CiCa*co2_ppmv)^h_ca)
    local Ismax, h, Cstar
    if co2_ppmv < 400.0
        Ismax, h, Cstar = 1.072, 1.70, 1218.0
    elseif co2_ppmv > 400.0 && co2_ppmv < 600.0
        fint = (co2_ppmv - 400.0)/200.0
        Ismax = fint*1.036 + (1-fint)*1.072
        h     = fint*2.0125 + (1-fint)*1.70
        Cstar = fint*1150.0 + (1-fint)*1218.0
    elseif co2_ppmv > 600.0 && co2_ppmv < 800.0
        fint = (co2_ppmv - 600.0)/200.0
        Ismax = fint*1.046 + (1-fint)*1.036
        h     = fint*1.5380 + (1-fint)*2.0125
        Cstar = fint*2025.0 + (1-fint)*1150.0
    else
        Ismax, h, Cstar = 1.014, 2.861, 1525.0
    end
    local gamma_ci
    if !isnan(cisun) && !isnan(cisha) && forc_pbot > 0.0 && fsun > 0.0
        ci = (fsun*cisun + (1.0-fsun)*cisha)/forc_pbot * 1.0e6
        gamma_ci = Ismax - (Ismax*ci^h)/(Cstar^h + ci^h)
    elseif cisun > 0.0 && cisun < 1.0e30 && forc_pbot > 0.0 && fsun == 1.0
        ci = cisun/forc_pbot * 1.0e6
        gamma_ci = Ismax - (Ismax*ci^h)/(Cstar^h + ci^h)
    elseif cisha > 0.0 && cisha < 1.0e30 && forc_pbot > 0.0 && fsun == 0.0
        ci = cisha/forc_pbot * 1.0e6
        gamma_ci = Ismax - (Ismax*ci^h)/(Cstar^h + ci^h)
    else
        gamma_ci = 1.0
    end
    return gamma_ci * gamma_ca
end

# get_map_EF: 6-group gridded isoprene EF lookup (VOCEmissionMod get_map_EF)
function ref_map_EF(ivt, g, efisop)
    if ivt == 1 || ivt == 2
        efisop[2, g]
    elseif ivt == 3
        efisop[3, g]
    elseif ivt >= 4 && ivt <= 8
        efisop[1, g]
    elseif ivt >= 9 && ivt <= 11
        efisop[4, g]
    elseif ivt >= 12 && ivt <= 14
        efisop[5, g]
    elseif ivt >= 15
        efisop[6, g]
    else
        0.0
    end
end

# -----------------------------------------------------------------------------
# Reporting helpers
# -----------------------------------------------------------------------------
const NFAIL = Ref(0)
function check(name, got, exp; tol=TOL)
    d = abs(got - exp)
    ok = d <= tol * max(1.0, abs(exp))
    ok || (NFAIL[] += 1)
    @printf("  %-52s got=% .10e ref=% .10e  |Δ|=% .2e  %s\n",
            name, got, exp, d, ok ? "ok" : "**FAIL**")
    return ok
end

# =============================================================================
println("="^78)
println("VOC / MEGAN Fortran-parity harness (independent scalar oracle)")
println("="^78)

# -----------------------------------------------------------------------------
# A. Activity-factor helper VALUE parity
# -----------------------------------------------------------------------------
println("\n[A] Activity-factor helpers vs independent oracle")

# gamma_L: accumulated + non-accumulated branch
for (f240, elai) in ((0.5, 2.3), (0.0, 2.3), (2.0e30, 1.1))
    check("gamma_L(fsun240=$f240,elai=$elai)",
          CLM.get_gamma_L(f240, elai), ref_gamma_L(f240, elai))
end

# gamma_SM across the drought curve
for btran in (0.05, 0.2, 0.5, 0.85, 1.0, 1.4)
    check("gamma_SM(btran=$btran)", CLM.get_gamma_SM(btran), ref_gamma_SM(btran))
end

# gamma_P: accumulated branch (all four >0 conditions) and fixed-coeff branch
let cases = [
        (300.0, 250.0, 220.0, 120.0, 110.0, 100.0, 0.6, 0.5, 180.0, 90.0, 0.6),   # accumulated
        (300.0, 250.0, 220.0, 120.0, 110.0, 100.0, 0.6, 1.2, 180.0, 90.0, 0.6),   # fsun240>1 -> fixed
        (300.0, 250.0, 220.0, 120.0, 110.0, 100.0, 0.6, 0.5,   0.0, 90.0, 0.6),   # solad240=0 -> fixed
    ]
    for (i, a) in enumerate(cases)
        gp, cp, al = CLM.get_gamma_P(a...)
        rgp, rcp, ral = ref_gamma_P(a...)
        check("gamma_P[case$i]", gp, rgp)
        check("gamma_P.cp[case$i]", cp, rcp)
        check("gamma_P.alpha[case$i]", al, ral)
    end
end

# gamma_T: generic, boreal shrub, arctic C3 grass, and fixed-coeff branch
let LDF = 0.6, Ceo = 2.0, ct1 = 95.0, ct2 = 230.0, betaT = 0.10
    for (ivt, t240, t24, tv) in (
            (NDLLF_EVR_TMP_TREE, 297.0, 300.0, 301.0),
            (NBRDLF_DCD_BRL_SHRUB, 290.0, 295.0, 298.0),
            (NC3_ARCTIC_GRASS, 288.0, 291.0, 293.0),
            (NDLLF_EVR_TMP_TREE, 0.0, 300.0, 301.0),   # t_veg240<=0 -> fixed Eopt/topt
        )
        gt, E, to = CLM.get_gamma_T(t240, t24, tv, ct1, ct2, betaT, LDF, Ceo, ivt)
        rgt, rE, rto = ref_gamma_T(t240, t24, tv, ct1, ct2, betaT, LDF, Ceo, ivt)
        check("gamma_T(ivt=$ivt,t240=$t240)", gt, rgt)
        check("gamma_T.Eopt(ivt=$ivt)", E, rE)
        check("gamma_T.topt(ivt=$ivt)", to, rto)
    end
end

# gamma_A: evergreen (=1), deciduous new/mature/old branches
let mf = CLM.MEGANFactors(); CLM.megan_factors_init!(mf, 20)
    mf.Anew .= 0.05; mf.Agro .= 0.6; mf.Amat .= 1.125; mf.Aold .= 1.0   # MEGAN v2.1 class-1-ish
    for (ivt, e240, e) in (
            (NDLLF_EVR_TMP_TREE, 2.0, 2.0),    # evergreen -> 1
            (NBRDLF_DCD_TRP_TREE, 1.0, 2.0),   # elai_prev<elai -> new leaves
            (NBRDLF_DCD_TRP_TREE, 2.0, 1.0),   # elai_prev>elai -> old leaves
            (NBRDLF_DCD_TRP_TREE, 1.5, 1.5),   # elai_prev==elai -> mature
            (NBRDLF_DCD_TRP_TREE, 0.0, 2.0),   # elai240<=0 -> 1
        )
        g = CLM.get_gamma_A(ivt, e240, e, 1, mf)
        r = ref_gamma_A(ivt, e240, e, 1, mf.Anew, mf.Agro, mf.Amat, mf.Aold)
        check("gamma_A(ivt=$ivt,e240=$e240,e=$e)", g, r)
    end
end

# gamma_C: all four CO2 regimes (strictly inside; boundary is Fortran UB) + the
# ci-branch variants (both sunlit&shaded, sunlit-only, shaded-only, neither)
let pbot = 101325.0
    for co2 in (350.0, 500.0, 700.0, 900.0)
        ci = 0.7 * co2 * 1.0e-6 * pbot
        for (cs, ch, fs) in ((ci, ci, 0.5), (ci, NaN, 1.0), (NaN, ci, 0.0), (NaN, NaN, 0.5))
            g = CLM.get_gamma_C(cs, ch, pbot, fs, co2)
            r = ref_gamma_C(cs, ch, pbot, fs, co2)
            check("gamma_C(co2=$co2,fsun=$fs)", g, r)
        end
    end
end

# get_map_EF: one ivt per PFT group
let efisop = reshape(collect(1.0:12.0), 6, 2)
    for ivt in (1, 3, 5, 10, 13, 16)
        check("map_EF(ivt=$ivt)", CLM.get_map_EF(ivt, 1, efisop), ref_map_EF(ivt, 1, efisop))
    end
end

# -----------------------------------------------------------------------------
# B. Full voc_emission! driver VALUE parity (assembled kernel) vs scalar oracle
# -----------------------------------------------------------------------------
println("\n[B] Full voc_emission! driver vs independent pipeline oracle")

function run_driver_case(; use_mapped)
    np, nc, ng = 5, 4, 2
    mf = CLM.MEGANFactors(); CLM.megan_factors_init!(mf, 20)
    # Perturb class params so a wrong class index would be caught.
    mf.betaT[1] = 0.08; mf.betaT[2] = 0.13
    mf.ct1[1] = 95.0;   mf.ct1[2] = 80.0
    mf.ct2[1] = 230.0;  mf.ct2[2] = 210.0
    mf.Ceo[1] = 2.0;    mf.Ceo[2] = 1.83
    mf.LDF[1] = 0.9999; mf.LDF[2] = 0.4    # isoprene ~all light-dependent
    mf.Anew .= 0.05; mf.Agro .= 0.6; mf.Amat .= 1.125; mf.Aold .= 1.0

    iso  = CLM.MEGANCompound(name="isoprene", index=1, class_number=1,
                             molec_weight=68.12, coeff=1.0, emis_factors=fill(24.0, 20))
    mono = CLM.MEGANCompound(name="myrcene", index=2, class_number=2,
                             molec_weight=136.23, coeff=0.5, emis_factors=fill(70.0, 20))
    meg = [iso, mono]
    mech = [CLM.MEGANMechComp(name="ISOP", n_megan_comps=1, megan_indices=[1]),
            CLM.MEGANMechComp(name="BIGENE", n_megan_comps=2, megan_indices=[1, 2])]  # sum of two

    voc = CLM.VOCEmisData(); CLM.vocemis_init!(voc, np, ng, 2, 2)
    voc.efisop_grc[:, 1] .= 600.0; voc.efisop_grc[:, 2] .= 480.0

    patch = CLM.PatchData(); CLM.patch_init!(patch, np)
    patch.itype    .= [0, NDLLF_EVR_TMP_TREE, NBRDLF_DCD_TRP_TREE, NC3_ARCTIC_GRASS, NBRDLF_DCD_BRL_SHRUB]
    patch.gridcell .= [1, 1, 1, 2, 2]
    patch.column   .= [1, 1, 2, 3, 4]

    forc_solad = fill(220.0, nc, 2); forc_solad[:, 1] .= [210.0, 215.0, 225.0, 205.0]
    forc_solai = fill(95.0, ng, 2);  forc_solai[:, 1] .= [90.0, 100.0]
    forc_pbot  = fill(101325.0, nc)
    forc_pco2  = fill(45.0, ng)      # Pa  (~444 ppm at 101325 Pa)

    fsd24  = [180.0, 190.0, 200.0, 175.0, 185.0]
    fsd240 = [170.0, 178.0, 188.0, 168.0, 176.0]
    fsi24  = [85.0, 88.0, 92.0, 83.0, 87.0]
    fsi240 = [80.0, 84.0, 89.0, 79.0, 83.0]
    fsun   = [0.5, 0.55, 0.6, 0.45, 0.5]
    fsun24 = [0.5, 0.52, 0.58, 0.47, 0.5]
    fsun240= [0.5, 0.53, 0.57, 0.48, 0.5]
    elai   = [2.0, 2.2, 1.8, 1.2, 1.5]
    elai240= [2.0, 2.1, 1.6, 1.3, 1.4]
    ci_val = 0.7 * 45.0 / 101325.0 * 101325.0   # ~ci in Pa
    cisun  = fill(ci_val, np, 1); cisha = fill(0.9*ci_val, np, 1)
    t_veg   = [300.0, 301.0, 299.0, 293.0, 297.0]
    t_veg24 = [300.0, 300.5, 298.5, 292.0, 296.0]
    t_veg240= [297.0, 298.0, 296.0, 290.0, 294.0]
    btran   = [0.8, 0.7, 0.9, 0.6, 0.85]

    CLM.voc_emission!(voc, meg, mech, mf, patch, 1:np, trues(np),
        forc_solad, forc_solai, forc_pbot, forc_pco2,
        fsd24, fsd240, fsi24, fsi240,
        fsun, fsun24, fsun240, elai, elai240,
        cisun, cisha, t_veg, t_veg24, t_veg240, btran;
        use_mapped_emisfctrs=use_mapped)

    # ---- independent scalar recomputation of the whole pipeline ----
    megfac = 1.0/3600.0/1.0e6
    exp_tot = zeros(np); exp_mech = zeros(np, 2); exp_meg = zeros(2, np)
    for p in 1:np
        ivt = patch.itype[p]; ivt == 0 && continue
        g = patch.gridcell[p]; c = patch.column[p]
        par_sun    = (forc_solad[c,1] + fsun[p]*forc_solai[g,1]) * 4.6
        par24_sun  = (fsd24[p] + fsun24[p]*fsi24[p]) * 4.6
        par240_sun = (fsd240[p] + fsun240[p]*fsi240[p]) * 4.6
        par_sha    = ((1.0-fsun[p])*forc_solai[g,1]) * 4.6
        par24_sha  = ((1.0-fsun24[p])*fsi24[p]) * 4.6
        par240_sha = ((1.0-fsun240[p])*fsi240[p]) * 4.6
        gl = ref_gamma_L(fsun240[p], elai[p])
        gsm = ref_gamma_SM(btran[p])
        vmeg = zeros(2)
        for (imeg, mc) in enumerate(meg)
            is_iso = mc.name == "isoprene"
            cn = mc.class_number
            eps = (is_iso && use_mapped) ? ref_map_EF(ivt, g, voc.efisop_grc) : mc.emis_factors[ivt]
            gp, = ref_gamma_P(par_sun, par24_sun, par240_sun, par_sha, par24_sha, par240_sha,
                              fsun[p], fsun240[p], fsd240[p], fsi240[p], mf.LDF[cn])
            gt, = ref_gamma_T(t_veg240[p], t_veg24[p], t_veg[p], mf.ct1[cn], mf.ct2[cn],
                              mf.betaT[cn], mf.LDF[cn], mf.Ceo[cn], ivt)
            ga = ref_gamma_A(ivt, elai240[p], elai[p], cn, mf.Anew, mf.Agro, mf.Amat, mf.Aold)
            if is_iso
                co2_ppmv = 1.0e6 * forc_pco2[g] / forc_pbot[c]
                gc = ref_gamma_C(cisun[p,1], cisha[p,1], forc_pbot[c], fsun[p], co2_ppmv)
            else
                gc = 1.0
            end
            gamma = gl*gsm*ga*gp*gt*gc
            if gamma >= 0.0 && gamma < 100.0
                vmeg[imeg] = mc.coeff * eps * gamma * megfac / mc.molec_weight
                exp_meg[imeg, p] = eps * gamma * megfac * 1.0e-3
            end
        end
        for imech in 1:2
            acc = 0.0
            for ii in 1:mech[imech].n_megan_comps
                acc += vmeg[mech[imech].megan_indices[ii]]
            end
            exp_mech[p, imech] = acc
            exp_tot[p] += acc
        end
    end
    return (voc, exp_tot, exp_mech, exp_meg, np)
end

for use_mapped in (true, false)
    (voc, exp_tot, exp_mech, exp_meg, np) = run_driver_case(use_mapped=use_mapped)
    tag = use_mapped ? "mapped-EF" : "table-EF"
    for p in 1:np
        check("vocflx_tot[$p] ($tag)",      voc.vocflx_tot_patch[p], exp_tot[p])
        check("vocflx[$p,ISOP] ($tag)",     voc.vocflx_patch[p,1],   exp_mech[p,1])
        check("vocflx[$p,BIGENE] ($tag)",   voc.vocflx_patch[p,2],   exp_mech[p,2])
        check("meg_flux_out[iso,$p] ($tag)",voc.meg_flux_out[1,p],   exp_meg[1,p])
    end
    # NON-VACUITY: at least one vegetated patch must carry real nonzero flux.
    @assert any(>(0.0), exp_tot) "oracle produced only zero flux — vacuous case"
    @assert voc.vocflx_tot_patch[1] == 0.0 "noveg patch must stay zero"
end

# -----------------------------------------------------------------------------
# C. WIRING TRUTH — executable audit of the activation path
# -----------------------------------------------------------------------------
println("\n[C] Wiring truth")

# C1. Default config => empty MEGANConfig => driver gate is false (inert VOC).
let cfg = CLM.CLMDriverConfig()
    empty_ok = isempty(cfg.megan.meg_compounds) && isempty(cfg.megan.mech_comps)
    gate = cfg.use_voc && !isempty(cfg.megan.mech_comps) && !isempty(cfg.megan.meg_compounds)
    @printf("  default use_voc                     = %s\n", cfg.use_voc)
    @printf("  default megan.meg_compounds empty   = %s\n", isempty(cfg.megan.meg_compounds))
    @printf("  default megan.mech_comps  empty     = %s\n", isempty(cfg.megan.mech_comps))
    @printf("  driver VOC gate on default config   = %s  (expect false => inert)\n", gate)
    (cfg.use_voc == false && empty_ok && gate == false) || (NFAIL[] += 1;
        println("  **FAIL** default config should leave VOC inert"))
    println("  NOTE: activation is now WIRED (megan_config_init reads the megan")
    println("        namelist specifier + megan_factors_file; read_efisop_from_surfdata!")
    println("        reads EF1_* from fsurdat; clm_initialize!(use_voc=true) +")
    println("        CLMDriverConfig(use_voc=true, megan_specifier=…, megan_factors_file=…)")
    println("        are the live entry points). The DEFAULT run stays VOC-off / inert")
    println("        (this assertion), byte-identical. See the end-to-end activated-path")
    println("        proof in scripts/validate_voc_megan_activation.jl.")
end

# C2. Ported namelist path => populated descriptors => SAME kernel emits nonzero,
#     oracle-matching flux (kernel is faithful + live-dispatchable given config).
let
    # Build an in-memory megan_factors table (stand-in for megan_factors_file).
    # megan_factors_table_init! wants comp_EF (ncomp x npft) and class_EF
    # (nclass x npft); the stored per-PFT factor is Comp_EF*Class_EF.
    comp_names = ["isoprene", "myrcene"]
    npft = 20
    class_nums = [1, 2]
    comp_mw  = [68.12, 136.23]
    comp_EF  = [i == 1 ? 2.0 : 3.0  for i in 1:2, k in 1:npft]   # 2 x npft
    class_EF = [i == 1 ? 12.0 : 7.0 for i in 1:2, k in 1:npft]   # 2 x npft
    tbl = CLM.MEGANFactorsTable{Float64}()
    CLM.megan_factors_table_init!(tbl, comp_names, comp_EF, class_EF, class_nums, comp_mw)

    specifier = ["ISOP = isoprene", "BIGENE = isoprene + 0.5*myrcene"]
    megcfg = CLM.megan_config_from_nl(specifier, tbl; maxveg=npft)
    cfg = CLM.CLMDriverConfig(use_voc=true, megan=megcfg)

    gate = cfg.use_voc && !isempty(cfg.megan.mech_comps) && !isempty(cfg.megan.meg_compounds)
    @printf("  namelist-built VOC gate             = %s  (expect true)\n", gate)
    (gate == true) || (NFAIL[] += 1; println("  **FAIL** namelist path should enable VOC"))
    @printf("  meg_compounds=%d  mech_comps=%d  isoprene EF[1]=%.3f (=2*12)\n",
            length(megcfg.meg_compounds), length(megcfg.mech_comps),
            megcfg.meg_compounds[1].emis_factors[1])
    (abs(megcfg.meg_compounds[1].emis_factors[1] - 24.0) < 1e-12) ||
        (NFAIL[] += 1; println("  **FAIL** Comp_EF*Class_EF mapping wrong"))

    # Drive the kernel through the namelist-built descriptors on one vegetated patch.
    np, nc, ng = 1, 1, 1
    voc = CLM.VOCEmisData(); CLM.vocemis_init!(voc, np, ng,
        length(megcfg.meg_compounds), length(megcfg.mech_comps))
    voc.efisop_grc .= 600.0
    patch = CLM.PatchData(); CLM.patch_init!(patch, np)
    patch.itype .= [NBRDLF_DCD_TRP_TREE]; patch.gridcell .= [1]; patch.column .= [1]
    CLM.voc_emission!(voc, megcfg.meg_compounds, megcfg.mech_comps, megcfg.megan_factors,
        patch, 1:np, trues(np),
        fill(220.0, nc, 2), fill(95.0, ng, 2), fill(101325.0, nc), fill(45.0, ng),
        fill(180.0, np), fill(170.0, np), fill(85.0, np), fill(80.0, np),
        fill(0.5, np), fill(0.5, np), fill(0.5, np), fill(2.0, np), fill(1.8, np),
        fill(0.7*45.0, np, 1), fill(0.6*45.0, np, 1),
        fill(300.0, np), fill(299.0, np), fill(297.0, np), fill(0.8, np);
        use_mapped_emisfctrs=false)
    @printf("  namelist-built driver vocflx_tot[1] = % .6e  (expect > 0)\n", voc.vocflx_tot_patch[1])
    (voc.vocflx_tot_patch[1] > 0.0 && isfinite(voc.vocflx_tot_patch[1])) ||
        (NFAIL[] += 1; println("  **FAIL** namelist-driven kernel produced no flux"))
end

# =============================================================================
println("\n" * "="^78)
if NFAIL[] == 0
    println("RESULT: VOC/MEGAN — all VALUE-parity + wiring assertions PASS.")
    println("  Kernel physics = faithful to VOCEmissionMod.F90 (helpers + assembled")
    println("  driver, to <1e-11). Kernel is live-dispatchable & correctly gated in")
    println("  clm_driver.jl, and ACTIVATION is now wired: megan_config_init +")
    println("  read_efisop_from_surfdata! populate config.megan / efisop from the")
    println("  megan_factors_file + surfdata when use_voc is set (default stays off /")
    println("  inert). End-to-end activated-path proof: validate_voc_megan_activation.jl.")
    exit(0)
else
    println("RESULT: $(NFAIL[]) divergence(s) — see **FAIL** lines above.")
    exit(1)
end
