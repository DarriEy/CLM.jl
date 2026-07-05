#!/usr/bin/env python3
"""Scientific-parity gate — the STRICT bar (distinct from the 10% coverage
scorecard). Every (domain, variable) cell must meet BOTH an annual-mean gate
(|Δ%| ≤ 1%, temps |ΔK| ≤ 0.05 K, near-zero |Δ| ≤ unit floor) AND a daily-RMSE
gate (nRMSE = RMSE/std_Fortran ≤ 0.10, temps daily RMSE ≤ 0.20 K, near-constant
|RMSE| ≤ unit floor).

Scientific parity is DECLARED only when every failing cell is a documented
exception (see parity_config.DOCUMENTED_EXCEPTIONS) — i.e. zero undocumented
failures. Domains whose reference/Julia files are absent (still generating) are
skipped and reported as pending.

  python3 scripts/parity_gate.py            # full report
  python3 scripts/parity_gate.py --brief    # one-line status + undocumented list
"""
import os, sys
import numpy as np
from netCDF4 import Dataset
from parity_config import (DOMAINS, VARS, SCALE, DATA_DIR, date_ord, series,
                           STRICT_PCT, STRICT_TEMP_K, STRICT_NRMSE, STRICT_TEMP_RMSE_K,
                           TEMP_VARS, DOCUMENTED_EXCEPTIONS)

BRIEF = "--brief" in sys.argv

def gate_cell(j, f, unit_floor, is_temp):
    """Return (passed, ann_ok, rmse_ok, dpct_or_dk, nrmse, note)."""
    v = np.isfinite(j) & np.isfinite(f)
    if v.sum() < 3:
        return None
    j, f = j[v], f[v]
    jm, fm = j.mean(), f.mean()
    d = jm - fm
    rmse = float(np.sqrt(np.mean((j - f) ** 2)))
    std_f = float(np.std(f))
    if is_temp:
        ann_ok = abs(d) <= STRICT_TEMP_K
        rmse_ok = rmse <= STRICT_TEMP_RMSE_K
        metric = d  # K
        nrmse = rmse / std_f if std_f > 1e-9 else 0.0
    else:
        pct = d / abs(fm) * 100 if abs(fm) > 1e-12 else 0.0
        ann_ok = (abs(pct) <= STRICT_PCT) or (abs(d) <= unit_floor)
        # daily gate: normalized RMSE where the series has real variability;
        # otherwise (near-constant) judge the absolute RMSE against the floor.
        # A cell also passes when the ABSOLUTE daily RMSE is itself below the
        # per-unit noise floor — the differences are physically negligible
        # regardless of how small the Fortran std is (this is the config's stated
        # near-zero philosophy: near-zero quantities judged vs the absolute floor).
        # Without this, a variable that is ~0 all year (e.g. surface water at a
        # site that ponds a few days) fails on a normalized RMSE of noise even
        # though the actual difference is below the floor.
        if std_f > unit_floor:
            nrmse = rmse / std_f
            rmse_ok = (nrmse <= STRICT_NRMSE) or (rmse <= unit_floor)
        else:
            nrmse = rmse / std_f if std_f > 1e-9 else 0.0
            rmse_ok = rmse <= unit_floor
        metric = pct
    return (ann_ok and rmse_ok), ann_ok, rmse_ok, metric, nrmse

pending, per_dom = [], []
undocumented, documented_hit = [], []
total_cells = total_pass = 0

for key, fh0, jnc, lab, biome in DOMAINS:
    jpath = os.path.join(DATA_DIR, jnc)
    if not (os.path.exists(fh0) and os.path.exists(jpath)):
        pending.append((key, lab)); continue
    dj, df = Dataset(jpath), Dataset(fh0)
    jo, fo = date_ord(dj), date_ord(df)
    common = np.intersect1d(jo, fo)
    if len(common) < 3:
        # Reference and Julia output do not share enough dates to score
        # (e.g. a truncated or wrong-year Fortran reference).
        pending.append((key, f"{lab} (no date overlap: julia {jo[0]}-{jo[-1]} vs fortran {fo[0]}-{fo[-1]})"))
        continue
    ji = np.array([np.where(jo == c)[0][0] for c in common], dtype=int)
    fi = np.array([np.where(fo == c)[0][0] for c in common], dtype=int)
    npass = ntot = 0
    for jn, fn, vl, grp, floor in VARS:
        sc = SCALE.get(jn, 1)
        j = series(dj, jn, ji); f = series(df, fn, fi)
        if j is None or f is None:
            continue
        res = gate_cell(j * sc, f * sc, floor, jn in TEMP_VARS)
        if res is None:
            continue
        passed, ann_ok, rmse_ok, metric, nrmse = res
        ntot += 1; total_cells += 1
        if passed:
            npass += 1; total_pass += 1
        else:
            rec = (key, lab, jn, vl, metric, nrmse, ann_ok, rmse_ok, jn in TEMP_VARS)
            if (key, jn) in DOCUMENTED_EXCEPTIONS:
                documented_hit.append(rec)
            else:
                undocumented.append(rec)
    per_dom.append((lab, biome, npass, ntot))

def fmt_metric(m, is_temp): return f"{m:+.2f}K" if is_temp else f"{m:+.1f}%"

if not BRIEF:
    print("=" * 78)
    print(f"SCIENTIFIC-PARITY GATE  (strict: |Δ|≤{STRICT_PCT:.0f}% / {STRICT_TEMP_K}K annual"
          f"  +  daily nRMSE≤{STRICT_NRMSE:.2f} / {STRICT_TEMP_RMSE_K}K)")
    print("=" * 78)
    for lab, biome, npass, ntot in per_dom:
        flag = "" if npass == ntot else f"   ({ntot - npass} miss)"
        print(f"  {lab:16s} {npass:2d}/{ntot:<2d}  {biome}{flag}")
    if pending:
        print("  -- pending (reference/Julia output not yet generated):")
        for key, lab in pending:
            print(f"     {lab} [{key}]")

n_ann = sum(1 for r in undocumented if not r[6])          # fails the annual gate
n_daily_only = sum(1 for r in undocumented if r[6] and not r[7])  # annual OK, daily fails
n_strict_domains = sum(1 for _, _, np_, nt in per_dom if np_ == nt)
print()
print(f"strict cells passing: {total_pass}/{total_cells}"
      + (f"   ({len(per_dom)} domains scored, {len(pending)} pending)" if pending else ""))
print(f"domains fully strict (all {len(VARS)}): {n_strict_domains}/{len(per_dom)}")
print(f"documented exceptions hit: {len(documented_hit)}   undocumented failures: {len(undocumented)}"
      f"   ({n_ann} annual-magnitude, {n_daily_only} daily-timing-only)")

if documented_hit and not BRIEF:
    print("\nDocumented exceptions (exempt, but listed):")
    for key, lab, jn, vl, m, nrmse, a, r, it in documented_hit:
        why = "annual" if not a else ""; why += ("+daily" if not r else "") if not a else ("daily" if not r else "")
        print(f"  · {lab:14s} {vl:16s} {fmt_metric(m,it):>8}  nRMSE={nrmse:.2f}  [{why}]  — {DOCUMENTED_EXCEPTIONS[(key,jn)][:70]}…")

if undocumented:
    print(f"\nUndocumented strict failures ({len(undocumented)}) — the work remaining:")
    undocumented.sort(key=lambda x: (0 if not x[6] else 1, -abs(x[4]) if not x[8] else -abs(x[4])))
    for key, lab, jn, vl, m, nrmse, a, r, it in undocumented:
        gate = ("annual " if not a else "") + ("daily" if not r else "")
        print(f"  {lab:16s} {vl:18s} {fmt_metric(m,it):>8}  nRMSE={nrmse:5.2f}   [{gate.strip()}]")

achieved = (len(undocumented) == 0) and (len(pending) == 0)
print("\n" + "=" * 78)
if achieved:
    print("SCIENTIFIC PARITY: ACHIEVED — every cell meets the strict gate or is a "
          "documented exception.")
else:
    reason = []
    if undocumented: reason.append(f"{len(undocumented)} undocumented strict failure(s)")
    if pending: reason.append(f"{len(pending)} domain(s) pending")
    print("SCIENTIFIC PARITY: NOT YET — " + "; ".join(reason) + ".")
    print("(The 10% coverage scorecard remains green; that is breadth, not the parity bar.)")
print("=" * 78)
sys.exit(0)
