#!/usr/bin/env python3
"""
Compare Julia CLM.jl vs Fortran CLM output — daily time series for key variables.

Usage:
    python3 scripts/plot_comparison.py [--year 2003] [--julia FILE] [--save]
    python3 scripts/plot_comparison.py --year 2002 --cycle 3        # spinup mode
    python3 scripts/plot_comparison.py --year 2003 --julia /tmp/clm_julia_restart_2003.nc --save
"""

import argparse
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timedelta
from pathlib import Path

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.titleweight": "bold",
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7.5,
    "figure.dpi": 120,
})

# ── Paths ──────────────────────────────────────────────────────────────────
FORTRAN_DIR = Path(
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/"
    "domain_Bow_at_Banff_lumped/optimization/CLM/dds_run_1/final_evaluation"
)
JULIA_DIR = Path("/tmp/clm_spinup")

# (display_name, julia_name, fortran_name, unit, scale_julia, scale_fortran)
VARIABLES = [
    ("2-m Air Temperature",      "TSA",          "TSA",          "K",      1.0,     1.0),
    ("Ground Temperature",       "T_GRND",       "TG",           "K",      1.0,     1.0),
    ("Absorbed Solar Rad.",      "FSA",          "FSA",          "W/m\u00b2",  1.0, 1.0),
    ("Latent Heat Flux",         "EFLX_LH_TOT", "EFLX_LH_TOT", "W/m\u00b2",  1.0, 1.0),
    ("Sensible Heat Flux",       "EFLX_SH_TOT", "FSH",          "W/m\u00b2",  1.0, 1.0),
    ("Canopy Evaporation",       "FCEV",         "FCEV",         "W/m\u00b2",  1.0, 1.0),
    ("Canopy Transpiration",     "FCTR",         "FCTR",         "W/m\u00b2",  1.0, 1.0),
    ("Ground Evaporation",       "FGEV",         "FGEV",         "W/m\u00b2",  1.0, 1.0),
    ("Snow Water Equivalent",    "H2OSNO",       "H2OSNO",       "mm",     1.0,     1.0),
    ("Snow Depth",               "SNOW_DEPTH",   "SNOW_DEPTH",   "m",      1.0,     1.0),
    ("Total Runoff",             "QRUNOFF",      "QRUNOFF",      "mm/d",   86400.0, 86400.0),
    ("Surface Runoff",           "QOVER",        "QOVER",        "mm/d",   86400.0, 86400.0),
]

COL_FORTRAN = "#2166ac"
COL_JULIA   = "#d6604d"


def load_var(ds, name, scale=1.0):
    v = ds[name][:]
    v = np.squeeze(v)
    if hasattr(v, "mask"):
        v = np.where(v.mask, np.nan, v.data)
    return v.astype(float) * scale


def smooth(x, window=7):
    kernel = np.ones(window) / window
    return np.convolve(x, kernel, mode="same")


def main():
    parser = argparse.ArgumentParser(description="CLM.jl vs Fortran CLM comparison plots")
    parser.add_argument("--year", type=int, default=2003)
    parser.add_argument("--cycle", type=int, default=None)
    parser.add_argument("--julia", type=str, default=None)
    parser.add_argument("--fortran", type=str, default=None)
    parser.add_argument("--save", action="store_true")
    parser.add_argument("--smooth", type=int, default=7)
    args = parser.parse_args()

    year, win = args.year, args.smooth

    # ── Resolve files ──────────────────────────────────────────────────────
    fortran_file = Path(args.fortran) if args.fortran else \
        FORTRAN_DIR / f"Bow_at_Banff_lumped.clm2.h0.{year}-01-01-00000.nc"

    if args.julia:
        julia_file = Path(args.julia)
        mode_label = "restart init" if "restart" in str(julia_file).lower() else "custom"
    elif args.cycle is not None:
        julia_file = JULIA_DIR / f"clm_hist_c{args.cycle:02d}_{year}.nc"
        mode_label = f"spinup cycle {args.cycle}"
    else:
        restart_file = Path(f"/tmp/clm_julia_restart_{year}.nc")
        if restart_file.exists():
            julia_file, mode_label = restart_file, "Fortran restart init"
        else:
            julia_file = JULIA_DIR / f"clm_hist_c03_{year}.nc"
            mode_label = "spinup cycle 3"

    for f in [fortran_file, julia_file]:
        if not f.exists():
            raise FileNotFoundError(f)

    ds_f = nc.Dataset(str(fortran_file))
    ds_j = nc.Dataset(str(julia_file))
    ndays = min(ds_j.dimensions["time"].size, ds_f.dimensions["time"].size)
    dates = [datetime(year, 1, 1) + timedelta(days=i) for i in range(ndays)]

    # ── Figure ─────────────────────────────────────────────────────────────
    nvar = len(VARIABLES)
    ncols, nrows = 3, (nvar + 2) // 3
    fig, axes = plt.subplots(nrows, ncols, figsize=(15, 2.8 * nrows),
                             constrained_layout=True)
    axes = axes.flatten()

    stats_text = []

    for i, (label, jname, fname, unit, sj, sf) in enumerate(VARIABLES):
        ax = axes[i]

        try:
            j_data = load_var(ds_j, jname, sj)
            f_data = load_var(ds_f, fname, sf)
        except (KeyError, IndexError):
            ax.set_visible(False)
            continue

        n = min(len(j_data), len(f_data), ndays)
        j_data, f_data, d = j_data[:n], f_data[:n], dates[:n]

        j_s = smooth(j_data, win) if win > 1 else j_data
        f_s = smooth(f_data, win) if win > 1 else f_data

        ax.plot(d, f_s, color=COL_FORTRAN, lw=1.2, label="Fortran")
        ax.plot(d, j_s, color=COL_JULIA, lw=1.2, label="Julia", ls="--")
        ax.fill_between(d, f_s, j_s, alpha=0.08, color=COL_JULIA)

        # Stats
        mask = np.isfinite(j_data) & np.isfinite(f_data)
        if mask.sum() > 0:
            bias = np.mean(j_data[mask] - f_data[mask])
            r = np.corrcoef(j_data[mask], f_data[mask])[0, 1]
            f_mean = np.mean(np.abs(f_data[mask]))
            pct = (bias / f_mean * 100) if f_mean > 1e-10 else 0.0
            stats_text.append(
                f"{label:28s}  bias={bias:+8.3f} {unit:6s} ({pct:+5.1f}%)  r={r:.4f}")

        ax.set_title(label)
        ax.set_ylabel(unit)
        ax.grid(True, alpha=0.15, lw=0.5)
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
        ax.xaxis.set_minor_locator(mdates.MonthLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))

    # Single shared legend at bottom of first row
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2,
               frameon=True, framealpha=0.9, edgecolor="0.8",
               bbox_to_anchor=(0.5, 1.02))

    for i in range(nvar, len(axes)):
        axes[i].set_visible(False)

    fig.suptitle(
        f"CLM.jl  vs  Fortran CLM5  \u2014  Bow at Banff {year}   "
        f"({mode_label}, {win}-day avg)",
        fontsize=12, fontweight="bold", y=1.06)

    ds_f.close()
    ds_j.close()

    # ── Stats table ────────────────────────────────────────────────────────
    print(f"\n{'='*78}")
    print(f"  Julia ({mode_label}) vs Fortran \u2014 {year}")
    print(f"{'='*78}")
    for s in stats_text:
        print(f"  {s}")
    print(f"{'='*78}\n")

    if args.save:
        outpath = Path(__file__).parent / "comparison.png"
        fig.savefig(str(outpath), dpi=200, bbox_inches="tight")
        print(f"Saved to {outpath}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
