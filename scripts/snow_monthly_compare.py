#!/usr/bin/env python3
"""Monthly Julia-vs-Fortran snow comparison for a domain.
Usage: python3 snow_monthly_compare.py <domain_lower> <julia.nc> <fortran_h0.nc>
Aligns by calendar month. Reports H2OSNO, SNOW_DEPTH, FSNO(FRAC_SNO), QSNOMELT.
Derives bulk snow density = H2OSNO / SNOW_DEPTH to separate accumulation vs compaction.
"""
import sys, numpy as np, netCDF4 as nc, cftime

dom, jf, ff = sys.argv[1], sys.argv[2], sys.argv[3]

VARS = ['H2OSNO', 'SNOW_DEPTH', 'FSNO', 'QSNOMELT']
# Fortran names -> the same; Julia may use FRAC_SNO for FSNO
ALIASES = {'FSNO': ['FSNO', 'FRAC_SNO'], 'QSNOMELT': ['QSNOMELT', 'QFLX_SNOMELT', 'QSNOMELT_ICE'],
           'SNOW_DEPTH': ['SNOW_DEPTH', 'SNOWDP'], 'H2OSNO': ['H2OSNO']}

def load(path):
    d = nc.Dataset(path)
    t = d.variables['time']
    dates = cftime.num2date(t[:], t.units, calendar=getattr(t, 'calendar', 'noleap'))
    months = np.array([dt.month for dt in dates])
    out = {}
    for v in VARS:
        arr = None
        for name in ALIASES[v]:
            if name in d.variables:
                arr = np.ma.filled(d.variables[name][:].astype(float), np.nan)
                break
        if arr is None:
            out[v] = None; continue
        # collapse to a single column timeseries
        while arr.ndim > 1:
            arr = np.nanmean(arr, axis=tuple(range(1, arr.ndim)))
        out[v] = arr
    return months, out

jm, jd = load(jf)
fm, fd = load(ff)

print(f"\n{'='*78}\n  SNOW MONTHLY PARITY — {dom.upper()}  (Julia vs Fortran h0)\n{'='*78}")
for v in VARS:
    if jd[v] is None or fd[v] is None:
        print(f"\n{v}: MISSING (julia={jd[v] is not None}, fortran={fd[v] is not None})")
        continue
    print(f"\n{v}")
    print(f"  {'Mon':>3} {'Julia':>11} {'Fortran':>11} {'Diff':>11} {'%err':>8}")
    for mo in range(1, 13):
        jv = np.nanmean(jd[v][jm == mo]) if (jm == mo).any() else np.nan
        fv = np.nanmean(fd[v][fm == mo]) if (fm == mo).any() else np.nan
        diff = jv - fv
        pct = 100 * diff / fv if (fv not in (0, np.nan) and abs(fv) > 1e-9) else np.nan
        print(f"  {mo:>3} {jv:>11.4g} {fv:>11.4g} {diff:>11.4g} {pct:>8.1f}")
    # annual mean
    jam, fam = np.nanmean(jd[v]), np.nanmean(fd[v])
    print(f"  ANN {jam:>11.4g} {fam:>11.4g} {jam-fam:>11.4g} {100*(jam-fam)/fam if abs(fam)>1e-9 else float('nan'):>8.1f}")

# bulk snow density diagnostic
print(f"\nBULK SNOW DENSITY  (H2OSNO/SNOW_DEPTH, kg/m3) — depth-vs-SWE compaction signature")
print(f"  {'Mon':>3} {'Julia':>11} {'Fortran':>11} {'Diff':>11}")
jdens = np.where(jd['SNOW_DEPTH'] > 0.01, jd['H2OSNO'] / np.maximum(jd['SNOW_DEPTH'], 1e-6), np.nan)
fdens = np.where(fd['SNOW_DEPTH'] > 0.01, fd['H2OSNO'] / np.maximum(fd['SNOW_DEPTH'], 1e-6), np.nan)
for mo in range(1, 13):
    jv = np.nanmean(jdens[jm == mo]) if (jm == mo).any() else np.nan
    fv = np.nanmean(fdens[fm == mo]) if (fm == mo).any() else np.nan
    print(f"  {mo:>3} {jv:>11.5g} {fv:>11.5g} {jv-fv:>11.5g}")
