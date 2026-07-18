#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Acquire single-point CLM DATM forcing (``clmforc.YYYY.nc``) from the public
ARCO-ERA5 zarr store on GCS, for the crop-CFT parity point.

Why this exists
---------------
The crop-CFT single-point surfdata (``surfdata_cropCFT_USplains_1pt.nc``,
lon -96.71 / lat 44.80) has no matching datm forcing on the box; every existing
SYMFLUENCE domain is elsewhere.  Running the full SYMFLUENCE acquisition
workflow would require defining a whole domain (DEM download, delineation,
basin-averaging) for what is, at a single point, a direct ERA5 grid-cell
extraction.  This script does exactly that extraction and writes the same
NetCDF layout SYMFLUENCE's ``CLMForcingProcessor`` emits, so the generated
files drop straight into a ``datm.streams.xml``.

The ERA5 -> CLM conversions are the standard ones:

===========  ==========================================  ===============
CLM var      ERA5 source                                 conversion
===========  ==========================================  ===============
TBOT         2m_temperature                              none (K)
PSRF         surface_pressure                            none (Pa)
QBOT         2m_dewpoint_temperature + surface_pressure  Tetens -> kg/kg
WIND         10m_{u,v}_component_of_wind                 sqrt(u^2+v^2)
FSDS         surface_solar_radiation_downwards           J/m2/h -> /3600
FLDS         surface_thermal_radiation_downwards         J/m2/h -> /3600
PRECTmms     total_precipitation                         m/h -> *1000/3600
===========  ==========================================  ===============

Validation mode
---------------
``--validate-merbleue`` re-derives forcing at the MerBleue point/period and
compares it against that domain's *existing* SYMFLUENCE-produced
``clmforc.2016.nc``.  This checks the conversion chain (especially the
accumulated-flux /3600 and the dewpoint->q inversion) against a known-good
reference produced by an independent code path, rather than trusting the units
documented in the store.  Per the repeated "it was the harness, not the port"
lesson, the forcing generator is validated *before* it is used as a parity
input.

Usage
-----
    python3 scripts/acquire_crop_forcing.py --validate-merbleue
    python3 scripts/acquire_crop_forcing.py --lat 44.80 --lon -96.71 \
        --years 2016 2017 --outdir <dir>
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

ZARR_STORE = "gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3"

ERA5_VARS = [
    "2m_temperature",
    "2m_dewpoint_temperature",
    "10m_u_component_of_wind",
    "10m_v_component_of_wind",
    "surface_pressure",
    "total_precipitation",
    "surface_solar_radiation_downwards",
    "surface_thermal_radiation_downwards",
]

# CLM DATM variable -> (units, long_name)
CLM_VARS = {
    "TBOT": ("K", "TBOT"),
    "QBOT": ("kg/kg", "QBOT"),
    "WIND": ("m/s", "WIND"),
    "FSDS": ("W/m2", "FSDS"),
    "FLDS": ("W/m2", "FLDS"),
    "PSRF": ("Pa", "PSRF"),
    "PRECTmms": ("mm/s", "PRECTmms"),
}

FILL = 1.0e36


def open_arco() -> xr.Dataset:
    """Open the public ARCO-ERA5 store anonymously, narrowed to our variables."""
    import gcsfs

    gcs = gcsfs.GCSFileSystem(token="anon")  # public, anonymous
    ds = xr.open_zarr(gcs.get_mapper(ZARR_STORE), consolidated=True, chunks={})
    missing = [v for v in ERA5_VARS if v not in ds.data_vars]
    if missing:
        raise RuntimeError(f"ARCO store is missing required variables: {missing}")
    return ds[ERA5_VARS]


def dewpoint_to_q(d2m: np.ndarray, sp: np.ndarray) -> np.ndarray:
    """Specific humidity (kg/kg) from dewpoint (K) and surface pressure (Pa).

    Saturation vapour pressure at the dewpoint *is* the actual vapour pressure.
    Uses the Tetens/Bolton formula over water, matching the ECMWF convention
    for 2m dewpoint.
    """
    tdc = d2m - 273.15
    e = 611.2 * np.exp(17.67 * tdc / (tdc + 243.5))  # Pa
    return 0.622 * e / (sp - 0.378 * e)


def extract_point(ds: xr.Dataset, lat: float, lon: float, start: str, end: str) -> xr.Dataset:
    """Select the nearest ERA5 cell and time window, and materialise it."""
    lon360 = lon % 360.0
    sel = ds.sel(latitude=lat, longitude=lon360, method="nearest")
    sel = sel.sel(time=slice(start, end))
    print(
        f"  ERA5 cell: lat={float(sel.latitude):.3f} lon={float(sel.longitude):.3f} "
        f"({lon360:.3f} requested), {sel.sizes['time']} hours",
        flush=True,
    )
    return sel.compute()


def to_clm(sel: xr.Dataset) -> dict[str, np.ndarray]:
    """Apply the ERA5 -> CLM variable/unit conversions."""
    sp = sel["surface_pressure"].values
    out = {
        "TBOT": sel["2m_temperature"].values,
        "PSRF": sp,
        "QBOT": dewpoint_to_q(sel["2m_dewpoint_temperature"].values, sp),
        "WIND": np.sqrt(
            sel["10m_u_component_of_wind"].values ** 2
            + sel["10m_v_component_of_wind"].values ** 2
        ),
        # Accumulated over the preceding hour: J/m2 -> W/m2, m -> mm/s.
        "FSDS": sel["surface_solar_radiation_downwards"].values / 3600.0,
        "FLDS": sel["surface_thermal_radiation_downwards"].values / 3600.0,
        "PRECTmms": sel["total_precipitation"].values * 1000.0 / 3600.0,
    }
    # ERA5 accumulations can carry tiny negative round-off; CLM rejects negative
    # fluxes.  Clamp only the physically non-negative ones.
    for k in ("FSDS", "FLDS", "PRECTmms"):
        out[k] = np.maximum(out[k], 0.0)
    return out


def write_clmforc(
    path: Path, times: pd.DatetimeIndex, data: dict[str, np.ndarray], lat: float, lon: float
) -> None:
    """Write one ``clmforc.YYYY.nc`` in the exact SYMFLUENCE/CLM DATM layout."""
    # DATM reads time as "hours since 1900-01-01" on a standard calendar, with
    # the point dimensions named LATIXY/LONGXY (size 1 each).
    hours = (times - pd.Timestamp("1900-01-01")) / pd.Timedelta("1h")
    ds = xr.Dataset(
        data_vars={
            name: (
                ("time", "LATIXY", "LONGXY"),
                data[name].reshape(-1, 1, 1).astype("float64"),
                {"units": units, "long_name": long_name, "_FillValue": FILL},
            )
            for name, (units, long_name) in CLM_VARS.items()
        },
        coords={
            "time": ("time", np.asarray(hours, dtype="float64")),
            "LATIXY": ("LATIXY", np.array([lat], dtype="float64")),
            "LONGXY": ("LONGXY", np.array([lon % 360.0], dtype="float64")),
        },
    )
    ds["time"].attrs = {"units": "hours since 1900-01-01", "calendar": "standard"}
    ds["LATIXY_data"] = (
        ("LATIXY", "LONGXY"),
        np.array([[lat]], dtype="float64"),
        {"units": "degrees_north", "_FillValue": FILL},
    )
    ds["LONGXY_data"] = (
        ("LATIXY", "LONGXY"),
        np.array([[lon % 360.0]], dtype="float64"),
        {"units": "degrees_east", "_FillValue": FILL},
    )
    enc = {v: {"_FillValue": FILL} for v in ds.data_vars}
    enc["time"] = {"_FillValue": None, "units": "hours since 1900-01-01", "calendar": "standard"}
    for c in ("LATIXY", "LONGXY"):
        enc[c] = {"_FillValue": None}
    ds.to_netcdf(path, encoding=enc)
    print(f"  wrote {path} ({len(times)} times)", flush=True)


def validate_merbleue() -> int:
    """Re-derive MerBleue 2016 forcing and compare to the existing clmforc file.

    The existing file was produced by SYMFLUENCE's own ERA5 -> basin-average ->
    CLM chain, so agreement confirms this script's conversions.  Basin-averaging
    over a small lumped catchment vs. a single nearest ERA5 cell means exact
    equality is not expected; we check correlation and bias per variable.
    """
    ref_path = Path(
        "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Peatland_MerBleue_Canada"
        "/data/forcing/CLM_input/clmforc.2016.nc"
    )
    if not ref_path.exists():
        print(f"reference not found: {ref_path}", file=sys.stderr)
        return 1
    ref = xr.open_dataset(ref_path)
    ref_times = xr.decode_cf(ref)["time"].values

    lat, lon = 45.409, -75.519  # MerBleue pour point
    print(f"validating against {ref_path.name} at {lat}/{lon}")
    ds = open_arco()
    sel = extract_point(
        ds, lat, lon, str(pd.Timestamp(ref_times[0])), str(pd.Timestamp(ref_times[-1]))
    )
    mine = to_clm(sel)

    n = min(len(sel["time"]), ref.sizes["time"])
    print(f"\n{'var':10s} {'r':>8s} {'bias':>12s} {'mine_mean':>12s} {'ref_mean':>12s}")
    ok = True
    for name in CLM_VARS:
        a = np.asarray(mine[name][:n], dtype=float)
        b = np.asarray(ref[name].values[:n]).reshape(-1)
        good = np.isfinite(a) & np.isfinite(b) & (b < FILL / 10)
        a, b = a[good], b[good]
        r = float(np.corrcoef(a, b)[0, 1]) if a.std() > 0 and b.std() > 0 else float("nan")
        bias = float(a.mean() - b.mean())
        print(f"{name:10s} {r:8.4f} {bias:12.5g} {a.mean():12.5g} {b.mean():12.5g}")
        # Correlation is the structural check; a lumped basin average vs a
        # single cell legitimately differs in magnitude.
        if not (np.isnan(r) or r > 0.90):
            print(f"  ^^ {name} correlation below 0.90 -- conversion suspect")
            ok = False
    print("\nVALIDATION", "PASSED" if ok else "FAILED")
    return 0 if ok else 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--validate-merbleue", action="store_true")
    ap.add_argument("--lat", type=float, default=44.80)
    ap.add_argument("--lon", type=float, default=-96.71)
    ap.add_argument("--years", type=int, nargs="+", default=[2016, 2017])
    ap.add_argument("--outdir", type=Path)
    args = ap.parse_args()

    if args.validate_merbleue:
        return validate_merbleue()

    if args.outdir is None:
        ap.error("--outdir is required unless --validate-merbleue")
    args.outdir.mkdir(parents=True, exist_ok=True)

    ds = open_arco()
    for year in args.years:
        print(f"year {year}:", flush=True)
        # One extra hour, as the SYMFLUENCE files carry 8783/8784 steps
        # (hour 0 of Jan 1 through hour 23 of Dec 31 plus the wrap point).
        sel = extract_point(ds, args.lat, args.lon, f"{year}-01-01", f"{year}-12-31T23:00")
        times = pd.DatetimeIndex(sel["time"].values)
        write_clmforc(
            args.outdir / f"clmforc.{year}.nc", times, to_clm(sel), args.lat, args.lon
        )
    print("done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
