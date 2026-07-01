#!/usr/bin/env python3
"""30-var annual parity table: Julia vs Fortran h0.
Usage: python3 parity_table.py <label> <julia.nc> <fortran_h0.nc>
Reports annual-mean Julia, Fortran, abs diff, %err and a within-tolerance count.
Tolerance: |%err| <= 5% for non-tiny vars, or |diff| below a per-var abs floor.
"""
import sys, numpy as np, netCDF4 as nc

label, jf, ff = sys.argv[1], sys.argv[2], sys.argv[3]

# (name, alias_list, abs_floor) — abs_floor: if |Fortran annual| < floor, judge by abs diff
VARS = [
    ('FSH', ['FSH'], 2.0), ('EFLX_LH_TOT', ['EFLX_LH_TOT'], 2.0),
    ('FGEV', ['FGEV'], 1.5), ('FCEV', ['FCEV'], 1.5), ('FCTR', ['FCTR'], 1.5),
    ('FSA', ['FSA'], 3.0), ('FIRA', ['FIRA'], 3.0), ('FGR', ['FGR'], 2.0),
    ('FPSN', ['FPSN'], 0.3), ('TV', ['TV'], 0.5), ('TG', ['TG'], 0.5),
    ('TSA', ['TSA'], 0.5), ('TSOI', ['TSOI'], 0.5),
    ('H2OSNO', ['H2OSNO'], 5.0), ('SNOW_DEPTH', ['SNOW_DEPTH', 'SNOWDP'], 0.03),
    ('FSNO', ['FSNO', 'FRAC_SNO'], 0.03), ('QSNOMELT', ['QSNOMELT'], 2e-6),
    ('QRUNOFF', ['QRUNOFF'], 2e-6), ('QOVER', ['QOVER'], 2e-6),
    ('QDRAI', ['QDRAI'], 2e-6), ('QINFL', ['QINFL'], 2e-6),
    ('ZWT', ['ZWT'], 0.1), ('SOILWATER_10CM', ['SOILWATER_10CM'], 2.0),
    ('TOTSOILLIQ', ['TOTSOILLIQ'], 10.0), ('TWS', ['TWS'], 10.0),
    ('RAIN', ['RAIN'], 1e-6), ('SNOW', ['SNOW'], 1e-6),
    ('TLAI', ['TLAI'], 0.1), ('QFLX_EVAP_TOT', ['QFLX_EVAP_TOT'], 2e-6),
    ('QSOIL', ['QSOIL'], 2e-6),
]

def ann(path, aliases):
    d = nc.Dataset(path)
    for nm in aliases:
        if nm in d.variables:
            a = np.ma.filled(d.variables[nm][:].astype(float), np.nan)
            while a.ndim > 1:
                a = np.nanmean(a, axis=tuple(range(1, a.ndim)))
            return np.nanmean(a)
    return None

print(f"\n{'='*72}\n  PARITY TABLE — {label}\n{'='*72}")
print(f"  {'VAR':>16} {'Julia':>12} {'Fortran':>12} {'%err':>8} {'ok':>4}")
nok = ntot = 0
for name, al, floor in VARS:
    jv = ann(jf, al); fv = ann(ff, al)
    if jv is None or fv is None:
        print(f"  {name:>16} {'--':>12} {'--':>12}   MISSING")
        continue
    ntot += 1
    diff = jv - fv
    if abs(fv) > 1e-12:
        pct = 100 * diff / fv
    else:
        pct = float('nan')
    ok = (abs(diff) <= floor) or (not np.isnan(pct) and abs(pct) <= 5.0)
    nok += ok
    print(f"  {name:>16} {jv:>12.4g} {fv:>12.4g} {pct:>8.1f} {'OK' if ok else 'XX':>4}")
print(f"\n  within tolerance: {nok}/{ntot}")
