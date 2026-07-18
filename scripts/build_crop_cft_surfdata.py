#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""
Build a single-point crop-CFT (cft=64, "78-pft") CLM surface dataset for crop
parity validation, by extracting one crop-bearing gridcell from the global
unstructured ne3np4 78-pft surfdata that ships in the local CESM inputdata tree.

This unblocks the docs/CROP_PARITY.md "no crop-resolved surfdata exists on disk"
blocker: a real crop-CFT surfdata *does* exist locally (the ne3np4 global 78-pft
file); it is merely unstructured (gridcell-dimensioned), so the stock CTSM
`subset_data` tool (which expects a structured lsmlat/lsmlon global surfdata)
cannot slice it. This script does the equivalent gridcell->(lsmlat=1,lsmlon=1)
reshape that `subset_data.create_surfdata_at_point` performs, directly on the
unstructured file, and emits a matching single-point ESMF mesh.

The chosen cell (default g=352) is a US Great Plains corn/soybean/spring-wheat
point (lon 263.29 / -96.71, lat 44.80, PCT_CROP=45%) — the managed CFTs that
A3 (crop pools + crop N) and A4 (irrigation) exercise. Per CROP_PARITY.md the
site need not be Mead: single-step parity validates the crop *code path*.

Requires: xarray, numpy (present in the box python3).

Usage:
    python3 scripts/build_crop_cft_surfdata.py [--gridcell N] [--outdir DIR]
"""
from __future__ import annotations
import argparse
import os

import numpy as np
import xarray as xr

# Global unstructured 78-pft (crop-resolved) surfdata in the local inputdata tree.
DEFAULT_SRC = (
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/installs/cesm-inputdata/"
    "lnd/clm2/surfdata_esmf/ctsm5.4.0/surfdata_ne3np4_hist_2000_78pfts_c251022.nc"
)
DEFAULT_OUTDIR = (
    "/Users/darri.eythorsson/Library/CloudStorage/GoogleDrive-dareyt@gmail.com/"
    "My Drive/data/SYMFLUENCE_data/crop_cft_surfdata"
)


def build(src: str, outdir: str, gridcell_1based: int) -> None:
    os.makedirs(outdir, exist_ok=True)
    g0 = gridcell_1based - 1
    ds = xr.open_dataset(src, decode_times=False)
    lon = float(ds.LONGXY[g0])
    lat = float(ds.LATIXY[g0])

    # gridcell -> (lsmlat=1, lsmlon=1) reshape (subset_data.create_surfdata_at_point).
    sub = ds.isel(gridcell=g0)
    sub = sub.expand_dims({"lsmlat": 1, "lsmlon": 1})
    sub = sub.transpose(..., "lsmlat", "lsmlon")  # CLM expects lat/lon trailing
    if "PFTDATA_MASK" not in sub:
        sub["PFTDATA_MASK"] = xr.DataArray(
            np.ones((1, 1), dtype=np.int32), dims=["lsmlat", "lsmlon"])
    if "LANDFRAC_PFT" not in sub:
        sub["LANDFRAC_PFT"] = xr.DataArray(
            np.ones((1, 1)), dims=["lsmlat", "lsmlon"])
    sub.attrs["subset_note"] = (
        f"single gridcell g={gridcell_1based} extracted from "
        f"{os.path.basename(src)} (lon={lon:.3f} lat={lat:.3f}) "
        "for CLM.jl crop parity")

    surf_out = os.path.join(outdir, "surfdata_cropCFT_USplains_1pt.nc")
    sub.to_netcdf(surf_out, format="NETCDF4_CLASSIC")

    # Matching single-point ESMF mesh (mirrors SYMFLUENCE CLMDomainGenerator.
    # generate_esmf_mesh): a half-degree box, one quad element centred on the
    # surfdata gridcell. centerCoords MUST equal LONGXY/LATIXY.
    d = 0.5
    nlon = [lon - d, lon + d, lon + d, lon - d]
    nlat = [lat - d, lat - d, lat + d, lat + d]
    r_earth = 6.371e6
    area_rad = (1.0 * 1e6) / (r_earth ** 2)  # placeholder ~1 km^2 element area
    mesh = xr.Dataset(
        {
            "nodeCoords": xr.DataArray(
                np.array([[nlon[i], nlat[i]] for i in range(4)]),
                dims=["nodeCount", "coordDim"], attrs={"units": "degrees"}),
            "elementConn": xr.DataArray(
                np.array([[1, 2, 3, 4]], dtype=np.int32),
                dims=["elementCount", "maxNodePElement"],
                attrs={"long_name": "Node indices per element", "start_index": 1}),
            "numElementConn": xr.DataArray(
                np.array([4], dtype=np.int32), dims=["elementCount"]),
            "centerCoords": xr.DataArray(
                np.array([[lon, lat]]), dims=["elementCount", "coordDim"],
                attrs={"units": "degrees"}),
            "elementArea": xr.DataArray(
                np.array([area_rad]), dims=["elementCount"],
                attrs={"units": "radians^2"}),
            "elementMask": xr.DataArray(
                np.array([1], dtype=np.int32), dims=["elementCount"]),
            "origGridDims": xr.DataArray(
                np.array([1, 1], dtype=np.int32), dims=["origGridRank"]),
        },
        attrs={"gridType": "unstructured", "version": "0.9",
               "title": "ESMF mesh for Cropland_USplains crop point"})
    mesh_out = os.path.join(outdir, "esmf_mesh_croppt.nc")
    mesh.to_netcdf(mesh_out, format="NETCDF4")

    # Verify.
    v = xr.open_dataset(surf_out, decode_times=False)
    dims = {k: v.sizes[k] for k in ("lsmlat", "lsmlon", "natpft", "cft", "lsmpft")}
    print("WROTE", surf_out)
    print("WROTE", mesh_out)
    print("dims:", dims)
    print("PCT_CROP", float(v.PCT_CROP.values.squeeze()),
          "PCT_NATVEG", float(v.PCT_NATVEG.values.squeeze()),
          "PCT_CFT_sum", float(v.PCT_CFT.values.sum()))
    print("LONGXY", float(v.LONGXY.values.squeeze()),
          "LATIXY", float(v.LATIXY.values.squeeze()))
    assert dims["cft"] == 64 and dims["natpft"] == 15, "not a crop-CFT layout"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--src", default=DEFAULT_SRC)
    ap.add_argument("--outdir", default=DEFAULT_OUTDIR)
    ap.add_argument("--gridcell", type=int, default=352,
                    help="1-based gridcell index in the ne3np4 surfdata")
    args = ap.parse_args()
    build(args.src, args.outdir, args.gridcell)


if __name__ == "__main__":
    main()
