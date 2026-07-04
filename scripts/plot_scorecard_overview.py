#!/usr/bin/env python3
"""Multi-biome parity scorecard overview heatmap for the README.

Rows = biomes, columns = the 51 validated output variables (grouped by category).
Cell colour = signed relative error Julia-vs-Fortran (annual mean), diverging around
0. Near-zero quantities (surface water on dry/glacier sites) are judged by an
absolute noise floor and shown as agreement. Renders paper/parity_scorecard.png.

Paths point at the local SYMFLUENCE_data references; regenerate after a sweep.
"""
import os, numpy as np
from datetime import datetime
from netCDF4 import Dataset
import cftime
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Rectangle

# COVERAGE tier (|Δ|<10% / 0.5K or unit floor). The domain registry + variable
# set + per-unit floors live in the shared parity_config; the strict scientific-
# parity gate uses the same definitions via scripts/parity_gate.py.
from parity_config import DOMAINS, VARS, SCALE, DATA_DIR, date_ord, series

# Only score domains whose Fortran reference AND Julia output both exist; newly
# added biomes whose references are still generating are skipped until ready.
DOMAINS = [d for d in DOMAINS
           if os.path.exists(d[1]) and os.path.exists(os.path.join(DATA_DIR, d[2]))]

grid = np.full((len(DOMAINS), len(VARS)), np.nan)
passed = np.zeros((len(DOMAINS), len(VARS)), bool)
for i, (key, fh0, jnc, lab, biome) in enumerate(DOMAINS):
    dj = Dataset(os.path.join(DATA_DIR, jnc)); df = Dataset(fh0)
    jo, fo = date_ord(dj), date_ord(df)
    common = np.intersect1d(jo, fo)
    ji = np.array([np.where(jo==c)[0][0] for c in common])
    fi = np.array([np.where(fo==c)[0][0] for c in common])
    for k, (jn, fn, vl, grp, floor) in enumerate(VARS):
        sc = SCALE.get(jn, 1)
        j = series(dj, jn, ji); f = series(df, fn, fi)
        if j is None or f is None: continue
        jm, fm = np.nanmean(j)*sc, np.nanmean(f)*sc
        d = jm - fm
        pct = d/abs(fm)*100 if abs(fm) > 1e-9 else 0.0
        grid[i, k] = pct
        ok = (abs(d) < 0.5) if grp == "State" and "T" == vl.split()[-1][:1] and floor == 0.5 and vl not in ("Sat fraction",) else False
        # temperature vars: |Δ|<0.5K; else |%|<10 OR |Δ|<unit floor
        if floor == 0.5 and ("T" in vl):
            ok = abs(d) < 0.5
        else:
            ok = (abs(pct) < 10) or (abs(d) < floor)
        passed[i, k] = ok

npass = passed.sum(1)
print("Scorecard (within 10% / 0.5K):")
for i,(key,_,_,lab,biome) in enumerate(DOMAINS):
    print(f"  {lab:16s} {npass[i]:2d}/{len(VARS)}  {biome}")
print(f"TOTAL: {passed.sum()}/{passed.size}")

# ---- Plot ----
# For cells that pass only via the absolute noise floor (a physically near-zero
# quantity, e.g. surface water on a dry prairie), the relative % is meaningless —
# neutralise the displayed value so it isn't a misleading saturated colour.
disp = grid.copy()
for i in range(len(DOMAINS)):
    for k in range(len(VARS)):
        if passed[i, k] and np.isfinite(grid[i, k]) and abs(grid[i, k]) >= 10:
            disp[i, k] = 0.0

n_b, n_v = len(DOMAINS), len(VARS)
fig, ax = plt.subplots(figsize=(19, 6.6))
clip = np.clip(disp, -10, 10)
norm = TwoSlopeNorm(vmin=-10, vcenter=0, vmax=10)
im = ax.imshow(clip, aspect="auto", cmap="RdBu_r", norm=norm)
# mark the (few) genuine out-of-tolerance cells with an ×
for i in range(n_b):
    for k in range(n_v):
        if not passed[i, k] and np.isfinite(grid[i, k]):
            ax.text(k, i, "×", ha="center", va="center", color="black", fontsize=11, fontweight="bold")

ax.set_yticks(range(n_b))
ax.set_yticklabels([f"{lab}\n{biome}" for _,_,_,lab,biome in DOMAINS], fontsize=8.5)
ax.set_xticks(range(n_v))
ax.set_xticklabels([vl for _,_,vl,_,_ in VARS], rotation=90, fontsize=7.2)
# group separators + headers
groups = [g for _,_,_,g,_ in VARS]
bounds_g = [k for k in range(1, n_v) if groups[k] != groups[k-1]]
for b in bounds_g:
    ax.axvline(b-0.5, color="white", lw=2.5)
start = 0
for gi, g in enumerate(["Energy","Water","Snow","State","Carbon"]):
    idxs = [k for k in range(n_v) if groups[k]==g]
    if not idxs: continue
    ax.text(np.mean(idxs), -0.85, g, ha="center", va="bottom", fontsize=10, fontweight="bold")
for i in range(n_b+1): ax.axhline(i-0.5, color="white", lw=0.6)

cb = fig.colorbar(im, ax=ax, pad=0.008, fraction=0.018, extend="both")
cb.set_label("Julia − Fortran, annual mean (% of Fortran; clipped ±10%)", fontsize=9)
total = passed.sum()
ax.set_title(f"CLM.jl vs Fortran CLM5 — {n_b} biomes × {n_v} variables:  "
             f"{total}/{passed.size} within tolerance (|Δ| < 10% or < 0.5 K)  ·  × = out of tolerance",
             fontsize=13, fontweight="bold", pad=26)
plt.tight_layout()
out = os.path.join(os.path.dirname(__file__), "..", "scripts", "parity_scorecard.png")
plt.savefig(out, dpi=140, bbox_inches="tight")
print("saved", out)
