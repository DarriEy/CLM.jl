#!/usr/bin/env python3
"""Annual parity table: Julia vs Fortran h0 (30 core + extended diagnostics).
Usage: python3 parity_table.py <label> <julia.nc> <fortran_h0.nc>
Reports annual-mean Julia, Fortran, abs diff, %err and a within-tolerance count.
Tolerance: |%err| <= 5% for non-tiny vars, or |diff| below a per-var abs floor.
"""
import sys, numpy as np, netCDF4 as nc

label, jf, ff = sys.argv[1], sys.argv[2], sys.argv[3]

# (name, alias_list, abs_floor) — abs_floor: if |Fortran annual| < floor, judge by abs diff
CORE_VARS = [
    ('FSH', ['FSH'], 2.0), ('EFLX_LH_TOT', ['EFLX_LH_TOT'], 2.0),
    ('FGEV', ['FGEV'], 1.5), ('FCEV', ['FCEV'], 1.5), ('FCTR', ['FCTR'], 1.5),
    ('FSA', ['FSA'], 3.0), ('FIRA', ['FIRA'], 3.0), ('FGR', ['FGR'], 2.0),
    ('FPSN', ['FPSN'], 0.3), ('TV', ['TV'], 0.5), ('TG', ['TG', 'T_GRND'], 0.5),
    ('TSA', ['TSA'], 0.5), ('TSOI', ['TSOI'], 0.5),
    ('H2OSNO', ['H2OSNO'], 5.0), ('SNOW_DEPTH', ['SNOW_DEPTH', 'SNOWDP'], 0.03),
    ('FSNO', ['FSNO', 'FRAC_SNO'], 0.03), ('QSNOMELT', ['QSNOMELT', 'QFLX_SNOMELT'], 2e-6),
    ('QRUNOFF', ['QRUNOFF'], 2e-6), ('QOVER', ['QOVER'], 2e-6),
    ('QDRAI', ['QDRAI', 'QFLX_DRAIN'], 2e-6), ('QINFL', ['QINFL', 'QFLX_INFL'], 2e-6),
    ('ZWT', ['ZWT'], 0.1), ('SOILWATER_10CM', ['SOILWATER_10CM'], 2.0),
    ('TOTSOILLIQ', ['TOTSOILLIQ'], 10.0),
    # NOTE: gridcell TWS = endwb landunit-area-weighted + river storage; that
    # landunit weighting is not reproduced by the lumped single-column parity
    # harness (Julia would report the soil column's endwb un-weighted, ~14x the
    # glacier-diluted Fortran gridcell value at Bow). Total-water conservation is
    # covered by the T2 conservation harness instead.
    ('RAIN', ['RAIN'], 1e-6), ('SNOW', ['SNOW'], 1e-6),
    ('TLAI', ['TLAI'], 0.1), ('QFLX_EVAP_TOT', ['QFLX_EVAP_TOT'], 2e-6),
    ('QSOIL', ['QSOIL'], 2e-6),
]

# Extended diagnostic set — four families that pin the currently-open residuals:
#   radiation band-split (canopy-albedo bias), surface water / perched table
#   (cold-site drainage/snow residuals), canopy structure, near-surface turbulence.
EXT_VARS = [
    # radiation band-split (W/m2)
    ('FSR', ['FSR'], 3.0),
    ('FSRVD', ['FSRVD'], 1.5), ('FSRND', ['FSRND'], 1.5),
    ('FSRVI', ['FSRVI'], 1.5), ('FSRNI', ['FSRNI'], 1.5),
    ('FSDSVD', ['FSDSVD'], 3.0), ('FSDSND', ['FSDSND'], 3.0),
    ('FSDSVI', ['FSDSVI'], 3.0), ('FSDSNI', ['FSDSNI'], 3.0),
    ('SABV', ['SABV'], 3.0), ('SABG', ['SABG'], 3.0),
    ('FIRE', ['FIRE'], 3.0),
    # surface water / perched water table
    ('H2OSFC', ['H2OSFC'], 2.0), ('FH2OSFC', ['FH2OSFC'], 0.02),
    ('QH2OSFC', ['QH2OSFC'], 2e-6), ('ZWT_PERCH', ['ZWT_PERCH'], 0.2),
    # near-surface turbulence / momentum. TAUX/TAUY wind-stress components now match:
    # the CLMNCEP datm splits the scalar forcing wind equally (Sa_u=Sa_v=wind/sqrt2,
    # datm_datamode_clmncep_mod.F90:435), and forcing_reader.jl now does the same
    # (was u=wind,v=0 → all stress on taux ×sqrt2, tauy==0). Wind SPEED is unchanged
    # so no physics moved; only the diagnostic u/v split. taux/tauy now agree to <1%.
    ('U10', ['U10'], 0.3), ('Q2M', ['Q2M'], 3e-4), ('RH2M', ['RH2M'], 2.0),
    ('TAUX', ['TAUX'], 0.01), ('TAUY', ['TAUY'], 0.01),
    ('FSAT', ['FSAT'], 0.02),
    # canopy structure + transpiration/canopy-evap fluxes
    ('TSAI', ['TSAI'], 0.05), ('ELAI', ['ELAI'], 0.1),
    ('LAISUN', ['LAISUN'], 0.1), ('LAISHA', ['LAISHA'], 0.1),
    ('BTRANMN', ['BTRANMN', 'BTRAN'], 0.02),
    ('QVEGT', ['QVEGT'], 2e-6), ('QVEGE', ['QVEGE'], 2e-6),
    # soil thermal / ice. Julia writes a column-summed SOILICE; the matching
    # Fortran field is TOTSOILICE (Fortran's own SOILICE is per-layer). List
    # TOTSOILICE first so Fortran uses it and Julia falls through to SOILICE.
    ('TSOI_10CM', ['TSOI_10CM'], 0.5), ('SOILICE', ['TOTSOILICE', 'SOILICE'], 5.0),
]

VARS = CORE_VARS + EXT_VARS

def ann(path, aliases):
    d = nc.Dataset(path)
    for nm in aliases:
        if nm in d.variables:
            a = np.ma.filled(d.variables[nm][:].astype(float), np.nan)
            while a.ndim > 1:
                a = np.nanmean(a, axis=tuple(range(1, a.ndim)))
            return np.nanmean(a)
    return None

core_names = {v[0] for v in CORE_VARS}
print(f"\n{'='*72}\n  PARITY TABLE — {label}\n{'='*72}")
print(f"  {'VAR':>16} {'Julia':>12} {'Fortran':>12} {'%err':>8} {'ok':>4}")
nok = ntot = 0
core_ok = core_tot = ext_ok = ext_tot = 0
for name, al, floor in VARS:
    if name == EXT_VARS[0][0]:
        print(f"  {'-'*16} extended diagnostics {'-'*8}")
    jv = ann(jf, al); fv = ann(ff, al)
    if jv is None or fv is None:
        print(f"  {name:>16} {'--':>12} {'--':>12}   MISSING")
        continue
    ntot += 1
    is_core = name in core_names
    if is_core: core_tot += 1
    else: ext_tot += 1
    diff = jv - fv
    if abs(fv) > 1e-12:
        pct = 100 * diff / fv
    else:
        pct = float('nan')
    ok = (abs(diff) <= floor) or (not np.isnan(pct) and abs(pct) <= 5.0)
    nok += ok
    if ok:
        if is_core: core_ok += 1
        else: ext_ok += 1
    print(f"  {name:>16} {jv:>12.4g} {fv:>12.4g} {pct:>8.1f} {'OK' if ok else 'XX':>4}")
print(f"\n  core within tolerance:     {core_ok}/{core_tot}")
print(f"  extended within tolerance: {ext_ok}/{ext_tot}")
print(f"  TOTAL within tolerance:    {nok}/{ntot}")
