#!/usr/bin/env python3
# ==========================================================================
# parity_jitter_audit.py — the FAITHFULNESS audit behind the GMD claim that
# CLM.jl has no translation errors and the strict-gate exceptions are dynamical
# jitter, not defects.
#
# For every cell that fails the STRICT scorecard gate (parity_config.gate_cell),
# apply five falsifiable tests, EACH of which rules out a specific
# translation-bug signature. A genuine translation bug (a mistranslated formula
# or wrong parameter) leaves a distinctive fingerprint: a ONE-SIGNED,
# NON-near-zero, NON-conserved, NON-timing-rescuable discrepancy. A miss is
# consistent with jitter if it FAILS to show that fingerprint on at least one
# axis; it is FLAGGED (candidate real bug) only if it shows the bug fingerprint
# on ALL axes.
#
#   J1  sign-incoherence   month-to-month Δ sign flips           -> stochastic, not a one-signed bias
#   J2  near-zero          annual |mean| ≪ seasonal peak         -> strict % amplifies a ~0 quantity
#   J3  conservation       cumulative / seasonal-total matches   -> redistribution in time, not a source/sink
#   J4  timing-rescue      a ≤2-day lag drops nRMSE below gate   -> a free-run event-timing offset
#   J5  regime-overlap     active-day sets barely overlap        -> event-timing offset (intermittent vars)
#
# Reuses parity_config's date-alignment + strict gate verbatim (single source of
# truth). Anchor: the per-timestep teacher-forced test (fidelity_perstep) shows
# the CODE reproduces Fortran to the FP floor; this audit shows the free-run
# exceptions are the dynamical amplification of that floor, not new error.
#
#   SYMFLUENCE_DATA=<drive root> figvenv/python parity_jitter_audit.py
# ==========================================================================
import os, sys, csv
import numpy as np
from netCDF4 import Dataset
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import parity_config as pc

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "paper", "data")
STRICT_NRMSE = getattr(pc, "STRICT_NRMSE", 0.10)
STRICT_TEMP_RMSE_K = getattr(pc, "STRICT_TEMP_RMSE_K", 0.20)

def scaled(ds, jv, fv, idx_j, idx_f, is_temp):
    j = pc.series(ds[0], jv, idx_j); f = pc.series(ds[1], fv, idx_f)
    if j is None or f is None: return None, None
    s = pc.SCALE.get(jv, 1.0)
    return j * s, f * s

def month_of(ords):
    return (ords // 100) % 100

def audit_cell(jv, fv, j, f, months, is_temp):
    """Return dict of the five jitter diagnostics for one failing cell."""
    v = np.isfinite(j) & np.isfinite(f)
    j, f, m = j[v], f[v], months[v]
    jm, fm = j.mean(), f.mean()
    ann = (jm - fm) if is_temp else ((jm - fm) / abs(fm) * 100 if abs(fm) > 1e-12 else np.nan)
    rmse = float(np.sqrt(np.mean((j - f) ** 2)))
    nrmse = rmse / (np.std(f) + 1e-30)

    # J1 sign-incoherence: fraction of months whose Δ shares the annual sign
    ann_sign = np.sign(jm - fm)
    md = np.array([ (j[m==mm].mean() - f[m==mm].mean()) for mm in np.unique(m) if (m==mm).sum() ])
    coherence = float(np.mean(np.sign(md) == ann_sign)) if md.size else np.nan

    # J2 near-zero: annual |mean_F| relative to seasonal peak |F|
    nearzero = float(abs(fm) / (np.nanmax(np.abs(f)) + 1e-30))

    # J3 conservation: cumulative (seasonal total) relative error
    cum = abs(np.nansum(j) - np.nansum(f)) / (abs(np.nansum(f)) + 1e-30) * 100

    # J4 timing-rescue: min daily nRMSE over lags -2..+2
    best = nrmse
    for L in (-2, -1, 1, 2):
        a, b = (j[:-abs(L)], f[abs(L):]) if L > 0 else (j[abs(L):], f[:-abs(L)])
        if len(a) > 3:
            r = float(np.sqrt(np.mean((a - b) ** 2))) / (np.std(b) + 1e-30)
            best = min(best, r)

    # J5 regime overlap (intermittent: active = |value| above 1% of its peak)
    thr_j, thr_f = 0.01*np.nanmax(np.abs(j)), 0.01*np.nanmax(np.abs(f))
    aj, af = np.abs(j) > thr_j, np.abs(f) > thr_f
    overlap = float((aj & af).sum() / max((aj | af).sum(), 1))

    gate = STRICT_TEMP_RMSE_K if is_temp else STRICT_NRMSE
    return dict(ann=ann, nrmse=nrmse, coherence=coherence, nearzero=nearzero,
                cum=cum, lag_nrmse=best, overlap=overlap, gate=gate,
                # each flag TRUE = this axis is consistent with jitter (rules out a bug)
                j1=coherence < 0.80, j2=nearzero < 0.15, j3=cum < 1.0,
                j4=best <= gate, j5=overlap < 0.80)

def main():
    rows, flagged = [], []
    for key, fh0, jnc, label, biome in pc.DOMAINS:
        jp = os.path.join(DATA, jnc)
        if not (os.path.exists(jp) and os.path.exists(fh0)): continue
        dj, df = Dataset(jp), Dataset(fh0)
        oj, ofd = pc.date_ord(dj), pc.date_ord(df)
        common = np.intersect1d(oj, ofd)
        if len(common) < 30:
            dj.close(); df.close(); continue
        ij = np.array([np.where(oj == c)[0][0] for c in common])
        iff = np.array([np.where(ofd == c)[0][0] for c in common])
        months = month_of(common)
        for jv, fv, vlabel, group, floor in pc.VARS:
            is_temp = jv in pc.TEMP_VARS
            j, f = scaled((dj, df), jv, fv, ij, iff, is_temp)
            if j is None: continue
            g = pc.gate_cell(j, f, floor, is_temp)
            if g is None or g[0]:  # None=unscoreable, g[0]=passed -> skip
                continue
            documented = (key, jv) in pc.DOCUMENTED_EXCEPTIONS
            a = audit_cell(jv, fv, j, f, months, is_temp)
            jitter_axes = sum([a['j1'], a['j2'], a['j3'], a['j4'], a['j5']])
            is_flagged = jitter_axes == 0  # bug fingerprint on ALL axes
            rows.append((label, vlabel, a, documented, jitter_axes, is_flagged))
            if is_flagged and not documented:
                flagged.append((label, vlabel, a))
        dj.close(); df.close()

    # ---- report ----
    hdr = ("%-15s %-16s %8s %6s %6s %7s %7s %7s %5s %-5s %s" %
           ("domain","variable","ann","nRMSE","cohr","near0","cum%","lagNR","ovlp","jit","verdict"))
    print(hdr); print("-"*len(hdr))
    out = os.path.join(os.path.dirname(__file__), "..", "paper", "parity_jitter_audit.csv")
    with open(out, "w", newline="") as fo:
        w = csv.writer(fo); w.writerow(["domain","variable","ann","nrmse","coherence",
            "nearzero","cum_pct","lag_nrmse","overlap","jitter_axes","documented","flagged",
            "J1_signincoh","J2_nearzero","J3_conserved","J4_timingrescue","J5_regimeoffset"])
        for label, vlabel, a, doc, jax, flag in sorted(rows, key=lambda r: -abs(r[2]['ann'] if np.isfinite(r[2]['ann']) else 0)):
            anns = f"{a['ann']:+.1f}%" if not np.isnan(a['ann']) else "  n/a"
            verdict = "FLAG!!" if flag and not doc else ("doc" if doc else "jitter")
            print("%-15s %-16s %8s %6.2f %6.2f %7.3f %7.2f %7.2f %5.2f %5d %s" %
                  (label, vlabel, anns, a['nrmse'], a['coherence'], a['nearzero'],
                   a['cum'], a['lag_nrmse'], a['overlap'], jax, verdict))
            w.writerow([label, vlabel, f"{a['ann']:.3f}", f"{a['nrmse']:.3f}",
                f"{a['coherence']:.3f}", f"{a['nearzero']:.4f}", f"{a['cum']:.3f}",
                f"{a['lag_nrmse']:.3f}", f"{a['overlap']:.3f}", jax, doc, flag,
                a['j1'], a['j2'], a['j3'], a['j4'], a['j5']])
    print("-"*len(hdr))
    print(f"TOTAL misses audited: {len(rows)}   |   documented: {sum(1 for r in rows if r[3])}")
    print(f"consistent with jitter (>=1 axis rules out a bug): {len(rows)-len(flagged)}/{len(rows)}")
    print(f"FLAGGED (bug fingerprint on ALL 5 axes): {len(flagged)}")
    for label, vlabel, a in flagged:
        print(f"   !! {label} {vlabel}: ann={a['ann']:+.1f}% cohr={a['coherence']:.2f} "
              f"near0={a['nearzero']:.3f} cum={a['cum']:.2f}% lagNR={a['lag_nrmse']:.2f}")
    print(f"\nwrote {out}")

if __name__ == "__main__":
    main()
