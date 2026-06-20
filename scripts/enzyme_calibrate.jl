# =============================================================================
# #1 — GRADIENT-BASED PARAMETER CALIBRATION via the REVERSE chain. Backprop a model-output
# loss to a PHYSICAL PARAMETER through CLM.compositional_reverse! and recover its known value
# by gradient descent — the payoff of the reverse-AD work: d(Loss)/d(θ) from one backward pass.
#
# Target physics: decomp_rate_constants_bgc! → the soil-moisture decomposition scalar
#   w_scalar = log(minpsi/ψ) / log(minpsi/maxpsi)
# Calibrated knob: minpsi (CNSharedParamsData) — the dry-end limit of the decomposition
# moisture response (controls how fast soil C turns over under drought). cn_params is in the
# DIFFERENTIATED bundle, so the reverse yields db.cn_params.minpsi directly. The reverse
# gradient is first checked against central finite differences, then used to recover minpsi
# from a synthetic w_scalar observation.
#   julia +1.10 --project=/tmp/clm_jl10_<id> scripts/enzyme_calibrate.jl
# =============================================================================
using CLM, Enzyme, Printf
const C = CLM
C.varpar_init!(C.varpar, 1, 14, 2, 5)

const NC=1; const NLD=1; const NDP=7; const NCT=10; const NLEV=max(NLD,5)
const BGC_STATE = C.DecompBGCState()
const CASCADE   = C.DecompCascadeConData()
const PARAMS    = C.DecompBGCParams(bgc_initial_Cstocks=fill(200.0,NDP))
const SOILPSI   = fill(-1.0, NC, NLEV)
const ZSOI      = [0.01,0.04,0.09,0.16,0.26,0.40,0.58,0.80,1.06,1.36]
const COLDZ     = fill(0.1, NC, NLEV)
const TSOISNO   = fill(C.TFRZ + 15.0, NC, NLEV)
const MAXPSI    = -0.1                                   # held fixed at its true value
const AUX = (; mask=BitVector([true]), bounds=1:NC, nlevdecomp=NLD)

make_cn_params(minpsi) = C.CNSharedParamsData(Q10=1.5, minpsi=minpsi, maxpsi=MAXPSI,
    rf_cwdl2=0.0, tau_cwd=10.0, cwd_flig=0.24, froz_q10=1.5, decomp_depth_efolding=0.5, mino2lim=0.0)
fresh_cf() = (cf=C.SoilBiogeochemCarbonFluxData(); C.soil_bgc_carbon_flux_init!(cf, NC, NLEV, NDP, NCT); cf)
make_bundle(minpsi) = (; cf=fresh_cf(), cn_params=make_cn_params(minpsi))

phase!(b, a) = (C.decomp_rate_constants_bgc!(b.cf, BGC_STATE, PARAMS, b.cn_params, CASCADE;
    mask_bgc_soilc=a.mask, bounds=a.bounds, nlevdecomp=a.nlevdecomp, t_soisno=TSOISNO,
    soilpsi=SOILPSI, days_per_year=365.0, dt=1800.0, zsoi_vals=ZSOI, col_dz=COLDZ); nothing)

wscalar(minpsi) = (b=make_bundle(minpsi); phase!(b, AUX); b.cf.w_scalar_col[1,1])
const MINPSI_TRUE = -10.0
const W_OBS = wscalar(MINPSI_TRUE)
@printf("true θ:  minpsi=%.4f   →   w_scalar(ψ=-1) = %.6f (synthetic observation)\n", MINPSI_TRUE, W_OBS)

# Loss = (w_scalar/W_OBS − 1)². Reverse → d(Loss)/d(minpsi).
function loss_and_grad(minpsi)
    bf = make_bundle(minpsi); phase!(bf, AUX)
    Lval = (bf.cf.w_scalar_col[1,1]/W_OBS - 1)^2
    seed!(db, b) = (db.cf.w_scalar_col[1,1] = 2*(b.cf.w_scalar_col[1,1]/W_OBS - 1)/W_OBS)
    db = C.compositional_reverse!(Any[(phase!,(AUX,))], make_bundle(minpsi), seed!)
    return Lval, db.cn_params.minpsi
end
loss(minpsi) = (b=make_bundle(minpsi); phase!(b,AUX); (b.cf.w_scalar_col[1,1]/W_OBS - 1)^2)

# 1) Verify the reverse gradient against central FD at the (wrong) start point.
let m0 = -4.0
    _, grev = loss_and_grad(m0)
    gfd = (loss(m0+1e-4) - loss(m0-1e-4)) / 2e-4
    @printf("grad check @minpsi=%.1f:  reverse=% .6e  FD=% .6e  rel=%.2e  %s\n",
        m0, grev, gfd, abs(grev-gfd)/max(abs(gfd),1e-14),
        abs(grev-gfd)/max(abs(gfd),1e-14) < 1e-5 ? "MATCH ✓" : "MISMATCH ✗")
end

# 2) Recover minpsi by gradient descent (in log(−minpsi); minpsi<0). θ₀ deliberately wrong.
println("="^74)
function calibrate()
    l = log(4.0); local L
    @printf("start θ₀: minpsi=%.4f\n", -exp(l))
    for it in 1:220
        minpsi = -exp(l)
        L, g = loss_and_grad(minpsi)
        lr = 4.0 / sqrt(it)                    # decaying step (the landscape is steep far out)
        l -= lr * (minpsi * g)                 # dL/d(log(−minpsi)) = (−minpsi)·dL/dminpsi
        (it==1 || it%40==0) && @printf("  it %3d  L=%.3e  minpsi=%.5f (→ %.1f)\n", it, L, -exp(l), MINPSI_TRUE)
    end
    return -exp(l)
end
mf = calibrate()
err = abs(mf - MINPSI_TRUE)/abs(MINPSI_TRUE)
println("="^74)
@printf("RECOVERED: minpsi=%.6f  (true %.4f, rel err %.2e)   %s\n",
    mf, MINPSI_TRUE, err, err < 1e-3 ? "PASS ✓ (reverse-AD calibration recovered the parameter)" : "FAIL ✗")
println("="^74)
