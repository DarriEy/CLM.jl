"""Before/after strict-gate comparison of two Julia output directories.

usage: python3 scripts/parity_gate_diff.py SCRIPTS_DIR BEFORE_DIR AFTER_DIR [DOMAIN...]

Scores each (domain, var) cell of both runs with parity_config.gate_cell (the same
cell rule as parity_gate.py) and reports per-site pass counts plus every cell whose
metric moved. Both runs are scored by the SAME field names: there is no BTRAN ->
BTRANMN alias, because since campaign #3 the Julia files carry BTRANMN (the
Fortran daily minimum of hourly-mean btran) and BTRAN (the hourly mean) as separate
fields, and scoring BTRAN against Fortran's BTRANMN compares different statistics.
"""
import os, sys
import numpy as np
sys.path.insert(0, sys.argv[1])
from netCDF4 import Dataset
from parity_config import DOMAINS, VARS, SCALE, date_ord, series, gate_cell, TEMP_VARS

bdir, adir = sys.argv[2], sys.argv[3]
want = set(sys.argv[4:])
ALIAS_BEFORE = {}

def score(path, fh0, alias):
    dj, df = Dataset(path), Dataset(fh0)
    jo, fo = date_ord(dj), date_ord(df)
    common = np.intersect1d(jo, fo)
    ji = np.array([np.where(jo == c)[0][0] for c in common], dtype=int)
    fi = np.array([np.where(fo == c)[0][0] for c in common], dtype=int)
    out = {}
    for jn, fn, vl, grp, floor in VARS:
        sc = SCALE.get(jn, 1)
        j = series(dj, alias.get(jn, jn), ji); f = series(df, fn, fi)
        if j is None or f is None:
            continue
        res = gate_cell(j * sc, f * sc, floor, jn in TEMP_VARS)
        if res is None:
            continue
        out[jn] = res
    return out

tot_b = tot_a = tot_n = 0
for key, fh0, jnc, lab, biome in DOMAINS:
    if want and key not in want:
        continue
    bp, ap = os.path.join(bdir, jnc), os.path.join(adir, jnc)
    if not (os.path.exists(bp) and os.path.exists(ap) and os.path.exists(fh0)):
        print(f"{key}: missing files"); continue
    B = score(bp, fh0, ALIAS_BEFORE); A = score(ap, fh0, {})
    names = [v[0] for v in VARS if v[0] in B and v[0] in A]
    nb = sum(B[n][0] for n in names); na = sum(A[n][0] for n in names)
    tot_b += nb; tot_a += na; tot_n += len(names)
    print(f"\n== {key}: before {nb}/{len(names)}  after {na}/{len(names)}")
    for n in names:
        b, a = B[n], A[n]
        # metric = annual delta (% or K); nrmse = daily normalized RMSE
        bm, bn, am, an = b[3], b[4], a[3], a[4]
        if abs(am - bm) < 1e-9 and abs(an - bn) < 1e-9:
            continue
        t = "K" if n in TEMP_VARS else "%"
        dirn = "TOWARD" if (abs(am) <= abs(bm) and an <= bn + 1e-12) else (
               "AWAY" if (abs(am) > abs(bm) and an > bn) else "mixed")
        flag = ""
        if b[0] != a[0]:
            flag = "  PASS->FAIL" if b[0] else "  FAIL->PASS"
        print(f"  {n:16s} ann {bm:+9.3f}{t} -> {am:+9.3f}{t}   nRMSE {bn:8.4f} -> {an:8.4f}  {dirn}{flag}")
print(f"\nTOTAL before {tot_b}/{tot_n}  after {tot_a}/{tot_n}")
