#!/usr/bin/env python3
"""Curated daily-timeseries 'hero' figure for the README — CLM.jl vs Fortran CLM5.

Six representative variables for one domain with a strong seasonal cycle (default
Bow at Banff, alpine: full snow-accumulation/melt cycle, large seasonal swings in
radiation, turbulent fluxes, temperature, and photosynthesis). Fortran and CLM.jl
daily series are overlaid so the reader can see the agreement directly, complementing
the multi-biome heatmap. Renders scripts/parity_timeseries.png.

  DOMAIN=Bow python3 scripts/plot_timeseries_hero.py

Paths point at the local SYMFLUENCE_data references + paper/data Julia outputs;
regenerate after a parity sweep.
"""
import os, numpy as np, cftime
from datetime import datetime
from netCDF4 import Dataset
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

_D = os.environ.get("SYMFLUENCE_DATA", "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data")
DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "paper", "data")
plt.rcParams.update({'font.size': 9, 'font.family': 'serif'})

# (fortran h0, julia nc, title). Mirror plot_parity_full.py's registry.
DOMAINS = {
    "Bow": (f"{_D}/clm_parity_run/Bow_at_Banff_lumped.clm2.h0.2003-01-01-00000.nc",
            "julia_clm_bow_phs_2003.nc", "Bow at Banff, Canada — alpine (2003)"),
    "Krycklan": (f"{_D}/domain_Boreal_Krycklan_Sweden/simulations/clm_boreal/CLM/Boreal_Krycklan_Sweden.clm2.h0.2012-12-30-00000.nc",
            "julia_clm_krycklan_phs_2013.nc", "Krycklan, Sweden — boreal forest (2013)"),
    "HubbardBrook": (f"{_D}/domain_Temperate_HubbardBrook_USA/simulations/clm_tempforest/CLM/Temperate_HubbardBrook_USA.clm2.h0.2016-12-31-00000.nc",
            "julia_clm_hubbardbrook_phs_2017.nc", "Hubbard Brook, USA — temperate forest (2017)"),
}
DOMAIN = os.environ.get("DOMAIN", "Bow")
F_H0, F_JU, TITLE = DOMAINS[DOMAIN]

# (julia, fortran, label, unit, scale)
PANELS = [
    ("FSA", "FSA", "Absorbed shortwave", "W m⁻²", 1),
    ("EFLX_SH_TOT", "FSH", "Sensible heat", "W m⁻²", 1),
    ("EFLX_LH_TOT", "EFLX_LH_TOT", "Latent heat", "W m⁻²", 1),
    ("H2OSNO", "H2OSNO", "Snow water equivalent", "mm", 1),
    ("TSA", "TSA", "2-m air temperature", "K", 1),
    ("FPSN", "FPSN", "Gross photosynthesis", "µmol m⁻² s⁻¹", 1),
]

ds_f = Dataset(F_H0, "r")
ds_j = Dataset(os.path.join(DATA_DIR, F_JU), "r")

def date_ord(ds):
    t = ds["time"]; dts = cftime.num2date(t[:], t.units, getattr(t, "calendar", "noleap"))
    return np.array([d.year * 10000 + d.month * 100 + d.day for d in dts])

jo, fo = date_ord(ds_j), date_ord(ds_f)
common = np.intersect1d(jo, fo)
ji = np.array([np.where(jo == c)[0][0] for c in common])
fi = np.array([np.where(fo == c)[0][0] for c in common])
days = np.array([datetime(c // 10000, (c // 100) % 100, c % 100) for c in common])

def load(ds, var, idx):
    if var not in ds.variables:
        return None
    a = np.ma.filled(ds[var][:].astype("f8"), np.nan)
    if a.ndim > 1:
        a = a.reshape(a.shape[0], -1)[:, 0]
    a = a.flatten()[idx]
    return np.where(np.isfinite(a) & (a > -9990), a, np.nan)

fig, axes = plt.subplots(2, 3, figsize=(13.5, 6.2), squeeze=False)
for i, (jn, fn, lab, unit, sc) in enumerate(PANELS):
    ax = axes[i // 3, i % 3]
    j = load(ds_j, jn, ji); f = load(ds_f, fn, fi)
    if j is None or f is None:
        ax.text(0.5, 0.5, f"{lab}\nmissing", transform=ax.transAxes, ha="center"); ax.set_axis_off(); continue
    j, f = j * sc, f * sc
    valid = np.isfinite(f) & np.isfinite(j)
    r = np.corrcoef(f[valid], j[valid])[0, 1] if valid.sum() > 2 else np.nan
    fo_m, jo_m = np.nanmean(f), np.nanmean(j)
    d = jo_m - fo_m
    dpct = d / abs(fo_m) * 100 if abs(fo_m) > 1e-8 else np.nan
    ax.plot(days, f, "k-", lw=1.6, label="Fortran CLM5", alpha=0.85)
    ax.plot(days, j, "#2171b5", lw=1.3, ls="--", label="CLM.jl", alpha=0.95)
    ax.fill_between(days, f, j, color="#2171b5", alpha=0.15)
    st = f"r = {r:.3f}   Δ = {d:+.2f} ({dpct:+.1f}%)" if not np.isnan(dpct) else f"r = {r:.3f}   Δ = {d:+.3f}"
    ax.set_title(f"{lab}  ({unit})", fontsize=10, fontweight="bold")
    ax.text(0.03, 0.96, st, transform=ax.transAxes, fontsize=8, va="top",
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.8", alpha=0.9))
    if i == 0:
        ax.legend(fontsize=8.5, loc="upper right", framealpha=0.9)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
    ax.tick_params(labelsize=8)
    ax.grid(True, alpha=0.25, lw=0.5)

fig.suptitle(f"CLM.jl vs Fortran CLM5 (PHS + LUNA), daily means — {TITLE}",
             fontsize=13, fontweight="bold")
plt.tight_layout(rect=[0, 0, 1, 0.96])
out = os.path.join(os.path.dirname(__file__), "parity_timeseries.png")
plt.savefig(out, dpi=140, bbox_inches="tight")
print("saved", out)
ds_f.close(); ds_j.close()
