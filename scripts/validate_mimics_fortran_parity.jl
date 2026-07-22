# ==========================================================================
# validate_mimics_fortran_parity.jl
#
# NUMERIC validation of the MIMICS (Wieder 2015) soil-decomposition cascade now
# that it is WIRED onto the live path (decomp_method==2). Replaces the old
# dead-wiring tripwire (scripts/mimics_wiring_audit.jl on the unmerged branch
# bgc/mimics-fortran-parity, which ASSERTED that MIMICS was dispatched from
# nowhere). MIMICS had never been numerically validated — its unit tests are
# qualitative/monotonicity only, with no Fortran oracle.
#
# Ground truth = an INDEPENDENT plain-scalar re-implementation of the Fortran
# `decomp_rates_mimics` kinetics (SoilBiogeochemDecompCascadeMIMICSMod.F90, the
# multi-level nlevdecomp>1 / non-FATES / non-anoxia path), fed byte-identical
# inputs and the SAME clm5_params `mimics_*` values readParams reads on the
# `soil_decomp_method='MIMICSWieder2015'` Fortran path (see
# mimics_default_read_params!). This is the CNDV/VOC-harness technique: a
# value-level oracle, not a self-consistency check. We assert the VALUES of
# decomp_k (all 8 pools), pathfrac (all 15 transitions), rf, the microbial C:N
# (cn_col), and the w/o scalars — not mere finiteness ([[conservation-is-not-accuracy]]).
#
# Also asserts the ACTIVATION wiring: CLMDriverConfig(decomp_method=2) derives the
# MIMICS 8-pool / 15-transition layout (i_litr_max=2, i_cwd=8), and the CN driver
# config carries decomp_method=2 through the facade sync.
#
#   julia +1.12 --project=. scripts/validate_mimics_fortran_parity.jl
#
# NOTE ON A FULL cesm.exe MIMICS RUN: the prebuilt Fortran exe on this machine is
# the SP compset (I2000Clm50SpRs, use_cn=.false.); MIMICS requires use_cn (soil
# BGC). No prebuilt Clm50Bgc use_cn exe exists — a real Fortran MIMICS reference
# needs a fresh I2000Clm50BgcRs build with soil_decomp_method='MIMICSWieder2015',
# not tractable this session. So this scalar oracle IS the value-level ground truth
# (the fallback the task sanctions), with the EXACT clm50_params mimics_* values.
#
# ✔ LIVE-PATH GAP NOW CLOSED: the fmet input `litr_lig_c_to_n_col`
# (ligninNratioAvg) is computed every step on the live path by
# `soil_bgc_carbon_flux_lignin_n_ratio!` (src/types/soil_bgc_carbon_flux.jl),
# called from cn_driver_summarize_fluxes! gated on decomp_method==MIMICS_DECOMP —
# a faithful port of the mimics_decomp block of
# SoilBiogeochemCarbonFluxType.F90::Summary (F90:961-1016). An activated MIMICS run
# now sees the REAL lignin:N ratio, not the former constant fmet=0.6375. That
# computation is oracle-checked in test/test_mimics_wiring.jl. This harness still
# exercises decomp_rates_mimics! itself across a RANGE of ligninNratioAvg (so it is
# NOT a vacuous 0.0 test). A fresh MIMICSWieder2015 Fortran build for a full
# value-level parity of the summary itself is still not tractable this session;
# the scalar oracle in test_mimics_wiring.jl is the sanctioned value-level ground
# truth for that piece ([[conservation-is-not-accuracy]]).
# ==========================================================================

using CLM
using Printf
using Random

const C = CLM

# ---- constants (must match CLM.varcon so the oracle is apples-to-apples) ----
const SECSPHR   = 3600.0
const SECSPDAY  = 86400.0
const TFRZ      = 273.15
const G_TO_MG   = 1.0e3
const CM3_TO_M3 = 1.0e-6
const PCT_TO_FRAC = 1.0e-2

# ---------------------------------------------------------------------------
# Independent scalar oracle of decomp_rates_mimics (multi-level, non-FATES,
# non-anoxia). Returns NamedTuple of expected output arrays.
# ---------------------------------------------------------------------------
function mimics_oracle(p::C.DecompMIMICSParams; t_soisno, soilpsi, cpools_vr, dz,
                       cellclay, ligninNratioAvg, annsum_npp_col,
                       minpsi, maxpsi, tau_cwd, days_per_year,
                       i_cop_mic, i_oli_mic)
    nc, nlev = size(t_soisno)
    # pool indices (fixed non-FATES layout)
    i_met=1; i_str=2; i_avl=3; i_chem=4; i_phys=5; i_cop=6; i_oli=7; i_cwd=8
    # transition indices
    i_l1m1=1; i_l1m2=2; i_l2m1=3; i_l2m2=4; i_s1m1=5; i_s1m2=6; i_s2s1=7
    i_s3s1=8; i_m1s1=9; i_m1s2=10; i_m1s3=11; i_m2s1=12; i_m2s2=13; i_m2s3=14; i_cwdl2=15

    decomp_k = fill(NaN, nc, nlev, 8)
    pathfrac = fill(NaN, nc, nlev, 15)
    rf       = fill(NaN, nc, nlev, 15)
    w_scalar = fill(NaN, nc, nlev)
    o_scalar = fill(NaN, nc, nlev)
    cn_cop   = fill(NaN, nc)
    cn_oli   = fill(NaN, nc)

    # unpack params
    vmod=p.mimics_vmod; vint=p.mimics_vint; vslope=p.mimics_vslope
    kmod=p.mimics_kmod; kint=p.mimics_kint; kslope=p.mimics_kslope
    mge=p.mimics_mge
    fmet_p1,fmet_p2,fmet_p3,fmet_p4 = p.mimics_fmet[1:4]
    fcr1,fcr2 = p.mimics_fchem_r[1:2]; fck1,fck2 = p.mimics_fchem_k[1:2]
    fpr1,fpr2 = p.mimics_fphys_r[1:2]; fpk1,fpk2 = p.mimics_fphys_k[1:2]
    ds1,ds2 = p.mimics_desorp[1:2]; ps1,ps2 = p.mimics_p_scalar[1:2]
    tr1,tr2 = p.mimics_tau_r[1:2]; tk1,tk2 = p.mimics_tau_k[1:2]
    tmn,tmx,tmf = p.mimics_tau_mod_min, p.mimics_tau_mod_max, p.mimics_tau_mod_factor
    ko_r,ko_k = p.mimics_ko_r, p.mimics_ko_k
    densdep = p.mimics_densdep; desQ10=p.mimics_desorpQ10; tref=p.mimics_t_soi_ref
    cn_mod=p.mimics_cn_mod_num; cnr=p.mimics_cn_r; cnk=p.mimics_cn_k
    nue = p.mimics_nue_into_mic  # unused here (rf only)

    k_frag = 1.0/(SECSPDAY*days_per_year*tau_cwd)

    # rf (time-independent)
    rf_l1m1 = 1-mge[1]; rf_l2m1 = 1-mge[2]; rf_s1m1 = 1-mge[3]
    rf_l1m2 = 1-mge[4]; rf_l2m2 = 1-mge[5]; rf_s1m2 = 1-mge[6]

    for c in 1:nc
        fmet = fmet_p1*(fmet_p2 - fmet_p3*min(fmet_p4, ligninNratioAvg[c]))
        tau_mod = min(tmx, max(tmn, sqrt(tmf*max(0.0, annsum_npp_col[c]))))
        tau_m1 = tr1*exp(tr2*fmet)*tau_mod/SECSPHR
        tau_m2 = tk1*exp(tk2*fmet)*tau_mod/SECSPHR
        cn_cop[c] = cnr*sqrt(cn_mod/fmet)
        cn_oli[c] = cnk*sqrt(cn_mod/fmet)
        fchem_m1 = min(1.0,max(0.0, fcr1*exp(fcr2*fmet)))
        fchem_m2 = min(1.0,max(0.0, fck1*exp(fck2*fmet)))

        for j in 1:nlev
            # w_scalar / o_scalar
            psi = min(soilpsi[c,j], maxpsi)
            w = psi > minpsi ? log(minpsi/psi)/log(minpsi/maxpsi) : 0.0
            w_scalar[c,j] = w
            o_scalar[c,j] = 1.0
            wdo = w * 1.0 * 1.0   # * depth_scalar(=1) * o_scalar(=1)

            clay = PCT_TO_FRAC*min(100.0, cellclay[c,j])
            desorp = ds1*exp(ds2*clay)
            fphys_m1 = min(1.0, fpr1*exp(fpr2*clay))
            fphys_m2 = min(1.0, fpk1*exp(fpk2*clay))
            p_scalar = 1.0/(ps1*exp(ps2*sqrt(clay)))

            tdeg = t_soisno[c,j]-TFRZ
            vmax_l1m1=exp(vslope[1]*tdeg+vint[1])*vmod[1]/SECSPHR
            vmax_l2m1=exp(vslope[2]*tdeg+vint[2])*vmod[2]/SECSPHR
            vmax_s1m1=exp(vslope[3]*tdeg+vint[3])*vmod[3]/SECSPHR
            vmax_l1m2=exp(vslope[4]*tdeg+vint[4])*vmod[4]/SECSPHR
            vmax_l2m2=exp(vslope[5]*tdeg+vint[5])*vmod[5]/SECSPHR
            vmax_s1m2=exp(vslope[6]*tdeg+vint[6])*vmod[6]/SECSPHR
            km_l1m1=exp(kslope[1]*tdeg+kint[1])*kmod[1]
            km_l2m1=exp(kslope[2]*tdeg+kint[2])*kmod[2]
            km_s1m1=exp(kslope[3]*tdeg+kint[3])*kmod[3]*p_scalar
            km_l1m2=exp(kslope[4]*tdeg+kint[4])*kmod[4]
            km_l2m2=exp(kslope[5]*tdeg+kint[5])*kmod[5]
            km_s1m2=exp(kslope[6]*tdeg+kint[6])*kmod[6]*p_scalar

            desorption = (desorp/SECSPHR)*desQ10*exp((tdeg-tref)/10.0)
            m1c = (cpools_vr[c,j,i_cop_mic]/dz[c,j])*G_TO_MG*CM3_TO_M3
            m2c = (cpools_vr[c,j,i_oli_mic]/dz[c,j])*G_TO_MG*CM3_TO_M3

            # met_lit
            t1=vmax_l1m1*m1c/(km_l1m1+m1c); t2=vmax_l1m2*m2c/(km_l1m2+m2c)
            decomp_k[c,j,i_met]=(t1+t2)*wdo
            if (t1+t2)!=0.0
                pathfrac[c,j,i_l1m1]=t1/(t1+t2); pathfrac[c,j,i_l1m2]=t2/(t1+t2)
            else; pathfrac[c,j,i_l1m1]=0.0; pathfrac[c,j,i_l1m2]=0.0; end
            # str_lit
            t1=vmax_l2m1*m1c/(km_l2m1+m1c); t2=vmax_l2m2*m2c/(km_l2m2+m2c)
            decomp_k[c,j,i_str]=(t1+t2)*wdo
            if (t1+t2)!=0.0
                pathfrac[c,j,i_l2m1]=t1/(t1+t2); pathfrac[c,j,i_l2m2]=t2/(t1+t2)
            else; pathfrac[c,j,i_l2m1]=0.0; pathfrac[c,j,i_l2m2]=0.0; end
            # avl_som
            t1=vmax_s1m1*m1c/(km_s1m1+m1c); t2=vmax_s1m2*m2c/(km_s1m2+m2c)
            decomp_k[c,j,i_avl]=(t1+t2)*wdo
            if (t1+t2)!=0.0
                pathfrac[c,j,i_s1m1]=t1/(t1+t2); pathfrac[c,j,i_s1m2]=t2/(t1+t2)
            else; pathfrac[c,j,i_s1m1]=0.0; pathfrac[c,j,i_s1m2]=0.0; end
            # phys_som
            decomp_k[c,j,i_phys]=desorption*1.0
            # chem_som (OXIDAT)
            t1=vmax_l2m1*m1c/(ko_r*km_l2m1+m1c); t2=vmax_l2m2*m2c/(ko_k*km_l2m2+m2c)
            decomp_k[c,j,i_chem]=(t1+t2)*wdo
            # cop_mic
            decomp_k[c,j,i_cop]=tau_m1*m1c^(densdep-1.0)*wdo
            favl=min(1.0,max(0.0,1.0-fphys_m1-fchem_m1))
            pathfrac[c,j,i_m1s1]=favl; pathfrac[c,j,i_m1s2]=fchem_m1
            # oli_mic
            decomp_k[c,j,i_oli]=tau_m2*m2c^(densdep-1.0)*wdo
            favl2=min(1.0,max(0.0,1.0-fphys_m2-fchem_m2))
            pathfrac[c,j,i_m2s1]=favl2; pathfrac[c,j,i_m2s2]=fchem_m2
            # cwd
            decomp_k[c,j,i_cwd]=k_frag*wdo

            # fixed pathfracs
            pathfrac[c,j,i_s2s1]=1.0; pathfrac[c,j,i_s3s1]=1.0
            pathfrac[c,j,i_m1s3]=fphys_m1; pathfrac[c,j,i_m2s3]=fphys_m2
            pathfrac[c,j,i_cwdl2]=1.0
            # rf
            rf[c,j,i_l1m1]=rf_l1m1; rf[c,j,i_l1m2]=rf_l1m2
            rf[c,j,i_l2m1]=rf_l2m1; rf[c,j,i_l2m2]=rf_l2m2
            rf[c,j,i_s1m1]=rf_s1m1; rf[c,j,i_s1m2]=rf_s1m2
            for k in (i_s2s1,i_s3s1,i_m1s1,i_m1s2,i_m1s3,i_m2s1,i_m2s2,i_m2s3)
                rf[c,j,k]=0.0
            end
            rf[c,j,i_cwdl2]=0.0  # rf_cwdl2 default 0
        end
    end
    return (; decomp_k, pathfrac, rf, w_scalar, o_scalar, cn_cop, cn_oli)
end

# ---- max relative diff, NaN-aware + finite assertion on the reference ----
function maxreldiff(port, ref)
    m = 0.0
    for i in eachindex(port, ref)
        pr = port[i]; rf = ref[i]
        if isnan(rf)
            continue  # oracle didn't set this slot
        end
        @assert isfinite(rf) "oracle produced non-finite reference at $i"
        d = abs(pr - rf)/(1.0 + max(abs(pr), abs(rf)))
        m = max(m, d)
    end
    return m
end

function run_oracle_parity()
    Random.seed!(20260721)
    nc = 4; nlev = 5
    # MIMICS state + params + cascade via the real init (validates unpacking too)
    p = C.DecompMIMICSParams(); C.mimics_default_read_params!(p)
    st = C.DecompMIMICSState()
    casc = C.DecompCascadeConData(); casc.initial_stock_soildepth = 0.0
    cnp = C.CNSharedParamsData()
    C.cn_shared_params_read!(cnp; q10_mr=1.5, minpsi_hr=-2.0, maxpsi_hr=-0.002,
        rf_cwdl2=0.0, tau_cwd=3.3333333, cwd_flig=0.24, decomp_depth_efolding=10.0,
        froz_q10=1.5, mino2lim=0.2, organic_max=130.0)
    # varied clay across columns so the texture params differ
    cellclay = [10.0 + 15.0*(c-1) + 2.0*(j-1) for c in 1:nc, j in 1:nlev]
    nue = C.init_decompcascade_mimics!(st, casc, p, cnp; cellclay=cellclay,
        bounds=1:nc, nlevdecomp=nlev, ndecomp_pools_max=8,
        ndecomp_cascade_transitions_max=15, spinup_state=0, use_fates=false)

    # carbon flux instance (outputs live here)
    cf = C.SoilBiogeochemCarbonFluxData()
    C.soil_bgc_carbon_flux_init!(cf, nc, nlev, 8, 15; nlevdecomp=nlev)

    # inputs — realistic ranges, varied per (c,j)
    t_soisno = [275.0 + 8.0*sin(0.7c + 0.3j) for c in 1:nc, j in 1:nlev]
    soilpsi  = [-0.02 - 0.5*rand() for c in 1:nc, j in 1:nlev]
    dz       = [0.02 + 0.05*(j-1) for c in 1:nc, j in 1:nlev]
    cpools_vr = zeros(nc, nlev, 8)
    for c in 1:nc, j in 1:nlev, k in 1:8
        cpools_vr[c,j,k] = 5.0 + 40.0*rand()   # nonzero microbe pools ⇒ nonzero M-M
    end
    ligninNratioAvg = [15.0 + 20.0*rand() for _ in 1:nc]
    annsum_npp_col  = [50.0 + 300.0*rand() for _ in 1:nc]

    # --- port ---
    C.decomp_rates_mimics!(cf, st, p, cnp, casc;
        mask_bgc_soilc = trues(nc),
        bounds = 1:nc, nlevdecomp = nlev,
        t_soisno = t_soisno, soilpsi = soilpsi,
        decomp_cpools_vr = cpools_vr, col_dz = dz,
        ligninNratioAvg = ligninNratioAvg, annsum_npp_col = annsum_npp_col,
        days_per_year = 365.0, dt = 1800.0,
        spinup_state = 0, use_lch4 = false, anoxia = false, use_fates = false)

    # --- oracle ---
    ref = mimics_oracle(p; t_soisno, soilpsi, cpools_vr, dz, cellclay,
        ligninNratioAvg, annsum_npp_col,
        minpsi = cnp.minpsi, maxpsi = cnp.maxpsi, tau_cwd = cnp.tau_cwd,
        days_per_year = 365.0, i_cop_mic = st.i_cop_mic, i_oli_mic = st.i_oli_mic)

    # --- compare ---
    tol = 1e-12
    results = Pair{String,Float64}[]
    push!(results, "decomp_k (8 pools)"   => maxreldiff(Array(cf.decomp_k_col), ref.decomp_k))
    push!(results, "pathfrac (15 trans)"  => maxreldiff(Array(cf.pathfrac_decomp_cascade_col), ref.pathfrac))
    push!(results, "rf (15 trans)"        => maxreldiff(Array(cf.rf_decomp_cascade_col), ref.rf))
    push!(results, "w_scalar"             => maxreldiff(Array(cf.w_scalar_col), ref.w_scalar))
    push!(results, "o_scalar"             => maxreldiff(Array(cf.o_scalar_col), ref.o_scalar))
    push!(results, "cn_col[cop_mic]"      => maxreldiff(Array(cf.cn_col)[:, st.i_cop_mic], ref.cn_cop))
    push!(results, "cn_col[oli_mic]"      => maxreldiff(Array(cf.cn_col)[:, st.i_oli_mic], ref.cn_oli))

    println("── MIMICS decomp_rates_mimics! vs independent Fortran scalar oracle ──")
    ok = true
    for (name, d) in results
        pass = d <= tol
        ok &= pass
        @printf("  %-22s  max|Δrel| = %.3e   %s\n", name, d, pass ? "PASS" : "FAIL")
    end
    # non-vacuous: assert the port actually produced nonzero, finite decomp_k
    dk = Array(cf.decomp_k_col)
    @assert all(isfinite, dk) "decomp_k has non-finite entries"
    @assert count(x -> x > 0, dk) > 0.5*length(dk) "decomp_k mostly zero — vacuous"
    println("  (non-vacuous: ", count(x->x>0, dk), "/", length(dk), " decomp_k entries > 0)")
    return ok
end

function run_activation_wiring()
    println("── MIMICS activation wiring (config derivation + facade sync) ──")
    ok = true
    cfg = C.CLMDriverConfig(use_cn=true, decomp_method=2)
    checks = [
        ("decomp_method",              cfg.decomp_method, 2),
        ("ndecomp_pools",              cfg.ndecomp_pools, 8),
        ("ndecomp_cascade_transitions",cfg.ndecomp_cascade_transitions, 15),
        ("i_litr_min",                 cfg.i_litr_min, 1),
        ("i_litr_max",                 cfg.i_litr_max, 2),
        ("i_cwd",                      cfg.i_cwd, 8),
    ]
    for (name, got, want) in checks
        pass = got == want; ok &= pass
        @printf("  %-30s = %-3d (want %d)  %s\n", name, got, want, pass ? "PASS" : "FAIL")
    end
    # CENTURY default untouched
    cfgc = C.CLMDriverConfig(use_cn=true)
    cpass = cfgc.decomp_method==1 && cfgc.ndecomp_pools==7 &&
            cfgc.ndecomp_cascade_transitions==10 && cfgc.i_litr_max==3 && cfgc.i_cwd==7
    ok &= cpass
    @printf("  %-30s %s\n", "CENTURY default (1/7/10/3/7)", cpass ? "PASS" : "FAIL")

    # facade sync propagates decomp_method into the CN driver config
    veg = C.CNVegetationData()
    veg.config.use_cn = true
    veg.config.decomp_method = 2
    C._sync_driver_config!(veg)
    spass = veg.driver_config.decomp_method == 2
    ok &= spass
    @printf("  %-30s %s\n", "facade _sync decomp_method=2", spass ? "PASS" : "FAIL")
    return ok
end

ok1 = run_oracle_parity()
ok2 = run_activation_wiring()
println()
println("✔ LIVE-PATH GAP CLOSED: litr_lig_c_to_n_col (ligninNratioAvg) is now computed")
println("  every step on the live path (soil_bgc_carbon_flux_lignin_n_ratio!, gated on")
println("  decomp_method==MIMICS) ⇒ an activated MIMICS run uses the REAL lignin:N ratio,")
println("  no longer the constant fmet=0.6375. That computation is oracle-checked in")
println("  test/test_mimics_wiring.jl; decomp_rates_mimics! itself is proved above across")
println("  varied ligninNratioAvg. See header note.")
println()
if ok1 && ok2
    println("ALL MIMICS PARITY + WIRING CHECKS PASSED")
    exit(0)
else
    println("MIMICS VALIDATION FAILED")
    exit(1)
end
