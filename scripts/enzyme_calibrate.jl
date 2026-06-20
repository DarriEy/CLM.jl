# =============================================================================
# #1 — GRADIENT-BASED PARAMETER CALIBRATION via the REVERSE chain. Backprop a model-output
# loss to PHYSICAL PARAMETERS through CLM.compositional_reverse! and recover their known values
# by gradient descent — the payoff of the reverse-AD work: d(Loss)/d(θ) for ALL parameters from
# ONE backward pass (vs ForwardDiff's one pass per parameter).
#
# Target physics: decomp_rate_constants_bgc! → the soil-C decomposition response scalars. We
# calibrate TWO real environmental-sensitivity knobs SIMULTANEOUSLY (one reverse pass per step):
#   Q10     (CNSharedParamsData) — temperature sensitivity: t_scalar ∝ Q10^((Tsoi-Tref)/10)
#   minpsi  (CNSharedParamsData) — dry-end moisture limit:  w_scalar = log(minpsi/ψ)/log(minpsi/maxpsi)
# cn_params is in the DIFFERENTIATED bundle, so one reverse pass yields db.cn_params.Q10 and
# db.cn_params.minpsi at once. Synthetic-observation recovery: make (t_scalar, w_scalar) from a
# "true" θ, then gradient-descend from a wrong θ₀ back to it.
#
# TWO gotchas this script resolves (both surfaced while making Q10 reverse-calibrate):
#  (1) t_scalar is normalized to be Q10-INDEPENDENT at Tref=15°C (normalize_q10_to_century_tfunc),
#      so d(t_scalar)/d(Q10)=0 there exactly — evaluate at Tsoi≠Tref (here +20°C) where Q10 bites.
#  (2) ALL decomposition taus must be set: an unset tau_s3_bgc → k_s3=1/(…·0)=Inf in decomp_k,
#      and Enzyme's 0·Inf (zero seed × non-finite output) POISONS the seeded shadow to NaN.
#      With finite decomp_k the Q10 (scalar-through-KA-kernel) reverse matches FD exactly.
#   julia +1.10 --project=/tmp/clm_jl10_<id> scripts/enzyme_calibrate.jl
# =============================================================================
using CLM, Enzyme, Printf
const C = CLM
C.varpar_init!(C.varpar, 1, 14, 2, 5)

const NC=1; const NLD=1; const NDP=7; const NCT=10; const NLEV=max(NLD,5)
const BGC_STATE = C.DecompBGCState()
const CASCADE   = C.DecompCascadeConData()
# All taus set → finite decomp_k (else k_s3=Inf poisons the reverse shadow via 0·Inf).
const PARAMS    = C.DecompBGCParams(tau_l1_bgc=1.0/18.5, tau_l2_l3_bgc=1.0/4.9, tau_s1_bgc=1.0/7.3,
    tau_s2_bgc=1.0/0.2, tau_s3_bgc=1.0/0.0045, bgc_initial_Cstocks=fill(200.0,NDP))
const SOILPSI   = fill(-1.0, NC, NLEV)
const ZSOI      = [0.01,0.04,0.09,0.16,0.26,0.40,0.58,0.80,1.06,1.36]
const COLDZ     = fill(0.1, NC, NLEV)
const TSOISNO   = fill(C.TFRZ + 20.0, NC, NLEV)        # +20°C ≠ Tref(15) → Q10 is active
const MAXPSI    = -0.1
const AUX = (; mask=BitVector([true]), bounds=1:NC, nlevdecomp=NLD)

make_cn_params(Q10, minpsi) = C.CNSharedParamsData(Q10=Q10, minpsi=minpsi, maxpsi=MAXPSI,
    rf_cwdl2=0.0, tau_cwd=10.0, cwd_flig=0.24, froz_q10=1.5, decomp_depth_efolding=0.5, mino2lim=0.0)
fresh_cf() = (cf=C.SoilBiogeochemCarbonFluxData(); C.soil_bgc_carbon_flux_init!(cf, NC, NLEV, NDP, NCT); cf)
make_bundle(Q10, minpsi) = (; cf=fresh_cf(), cn_params=make_cn_params(Q10, minpsi))

phase!(b, a) = (C.decomp_rate_constants_bgc!(b.cf, BGC_STATE, PARAMS, b.cn_params, CASCADE;
    mask_bgc_soilc=a.mask, bounds=a.bounds, nlevdecomp=a.nlevdecomp, t_soisno=TSOISNO,
    soilpsi=SOILPSI, days_per_year=365.0, dt=1800.0, zsoi_vals=ZSOI, col_dz=COLDZ); nothing)

scalars(Q10, minpsi) = (b=make_bundle(Q10, minpsi); phase!(b, AUX);
    (b.cf.t_scalar_col[1,1], b.cf.w_scalar_col[1,1]))
const Q10_TRUE, MINPSI_TRUE = 1.50, -10.0
const TS_OBS, WS_OBS = scalars(Q10_TRUE, MINPSI_TRUE)
@printf("true θ:  Q10=%.4f  minpsi=%.4f   →   t_scalar=%.6f  w_scalar=%.6f (synthetic obs)\n",
    Q10_TRUE, MINPSI_TRUE, TS_OBS, WS_OBS)

# Loss = (t_scalar/obs−1)² + (w_scalar/obs−1)². Reverse → gradient wrt BOTH params in one pass.
function loss_and_grad(Q10, minpsi)
    bf = make_bundle(Q10, minpsi); phase!(bf, AUX)
    ts, ws = bf.cf.t_scalar_col[1,1], bf.cf.w_scalar_col[1,1]
    Lval = (ts/TS_OBS - 1)^2 + (ws/WS_OBS - 1)^2
    seed!(db, b) = (db.cf.t_scalar_col[1,1] = 2*(b.cf.t_scalar_col[1,1]/TS_OBS - 1)/TS_OBS;
                    db.cf.w_scalar_col[1,1] = 2*(b.cf.w_scalar_col[1,1]/WS_OBS - 1)/WS_OBS)
    db = C.compositional_reverse!(Any[(phase!,(AUX,))], make_bundle(Q10, minpsi), seed!)
    return Lval, db.cn_params.Q10, db.cn_params.minpsi
end
loss(Q10, minpsi) = (b=make_bundle(Q10,minpsi); phase!(b,AUX);
    (b.cf.t_scalar_col[1,1]/TS_OBS-1)^2 + (b.cf.w_scalar_col[1,1]/WS_OBS-1)^2)

# 1) FD-verify BOTH reverse param-gradients at the (wrong) start point.
let q0=2.30, m0=-4.0
    _, gQ, gM = loss_and_grad(q0, m0)
    fQ = (loss(q0+1e-5,m0)-loss(q0-1e-5,m0))/2e-5
    fM = (loss(q0,m0+1e-4)-loss(q0,m0-1e-4))/2e-4
    @printf("grad check @θ₀:  dL/dQ10  rev=% .5e FD=% .5e (%s)   dL/dminpsi  rev=% .5e FD=% .5e (%s)\n",
        gQ, fQ, abs(gQ-fQ)/max(abs(fQ),1e-12)<1e-4 ? "✓" : "✗",
        gM, fM, abs(gM-fM)/max(abs(fM),1e-12)<1e-4 ? "✓" : "✗")
end

# 2) Recover BOTH params by gradient descent: Q10 in log space (>0), minpsi in log(−minpsi).
println("="^80)
function calibrate()
    lq = log(2.30); lm = log(4.0); local L      # θ₀: Q10=2.30, minpsi=-4.0 (both wrong)
    @printf("start θ₀: Q10=%.4f  minpsi=%.4f\n", exp(lq), -exp(lm))
    for it in 1:220
        Q10 = exp(lq); minpsi = -exp(lm)
        L, gQ, gM = loss_and_grad(Q10, minpsi)
        lr = 4.0/sqrt(it)
        lq -= lr * (Q10 * gQ)                    # dL/d(log Q10)     = Q10·dL/dQ10
        lm -= lr * (minpsi * gM)                 # dL/d(log(−minpsi))= (−minpsi)·dL/dminpsi
        (it==1 || it%40==0) && @printf("  it %3d  L=%.3e  Q10=%.4f (→%.2f)  minpsi=%.4f (→%.1f)\n",
            it, L, exp(lq), Q10_TRUE, -exp(lm), MINPSI_TRUE)
    end
    return exp(lq), -exp(lm)
end
Q10f, minf = calibrate()
eQ = abs(Q10f-Q10_TRUE)/abs(Q10_TRUE); eM = abs(minf-MINPSI_TRUE)/abs(MINPSI_TRUE)
println("="^80)
@printf("RECOVERED: Q10=%.5f (true %.4f, err %.1e)   minpsi=%.5f (true %.1f, err %.1e)   %s\n",
    Q10f, Q10_TRUE, eQ, minf, MINPSI_TRUE, eM,
    (eQ < 1e-3 && eM < 1e-3) ? "PASS ✓ (2-param reverse-AD calibration recovered both)" : "FAIL ✗")
println("="^80)
