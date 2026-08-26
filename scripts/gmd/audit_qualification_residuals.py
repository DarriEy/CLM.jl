#!/usr/bin/env python3
"""Reproduce diagnostics for the three residual cells in the GMD qualification.

The script reads the newly generated Julia outputs from CLM_QUALIFICATION_DIR and the
Fortran references under SYMFLUENCE_DATA. It prints JSON so the result can be archived
without scraping prose. No historical paper output is read.
"""
import json
import os
import sys

import numpy as np
from netCDF4 import Dataset

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import parity_config as pc


ROOT = os.environ["SYMFLUENCE_DATA"]
OUT = os.environ.get("CLM_QUALIFICATION_DIR", "/private/tmp/gmd-qualification")


def aligned(jpath, fpath, names):
    dj, df = Dataset(jpath), Dataset(fpath)
    jo, fo = pc.date_ord(dj), pc.date_ord(df)
    common = np.intersect1d(jo, fo)
    ji = np.array([np.where(jo == day)[0][0] for day in common], dtype=int)
    fi = np.array([np.where(fo == day)[0][0] for day in common], dtype=int)
    values = {}
    for out_name, jname, fname in names:
        values[out_name] = (pc.series(dj, jname, ji), pc.series(df, fname, fi))
    dj.close()
    df.close()
    return common, values


def finite_pair(pair):
    j, f = pair
    ok = np.isfinite(j) & np.isfinite(f)
    return j[ok], f[ok]


def pct(j, f):
    return 100.0 * (float(np.mean(j)) - float(np.mean(f))) / abs(float(np.mean(f)))


def nrmse(j, f):
    return float(np.sqrt(np.mean((j - f) ** 2)) / max(np.std(f), 1e-30))


def snow_audit(site, jpath, fpath):
    days, v = aligned(jpath, fpath, [
        ("swe", "H2OSNO", "H2OSNO"),
        ("depth", "SNOW_DEPTH", "SNOW_DEPTH"),
        ("grid_depth", "SNOWDP", "SNOWDP"),
        ("cover", "FRAC_SNO", "FSNO"),
        ("melt", "QFLX_SNOMELT", "QSNOMELT"),
    ])
    jswe, fswe = v["swe"]
    jdepth, fdepth = v["depth"]
    jcover, fcover = v["cover"]
    valid = np.isfinite(jswe) & np.isfinite(fswe) & np.isfinite(jdepth) & np.isfinite(fdepth)
    # The scientific floor used by the parity gate is 0.02 mm SWE. Require both
    # trajectories above it to remove trace frost from the physical-pack comparison.
    physical = valid & (jswe > 0.02) & (fswe > 0.02)
    density = physical & (jdepth > 1e-6) & (fdepth > 1e-6)
    jden, fden = jswe[density] / jdepth[density], fswe[density] / fdepth[density]
    jd, fd = jdepth[physical], fdepth[physical]
    jc, fc = finite_pair(v["cover"])
    jm, fm = finite_pair(v["melt"])
    jgrid, fgrid = finite_pair(v["grid_depth"])
    all_jd, all_fd = jdepth[valid], fdepth[valid]
    jeff, feff = (jdepth * jcover)[valid], (fdepth * fcover)[valid]
    depth_diff = np.where(valid, np.abs(jdepth - fdepth), -np.inf)
    worst = np.argsort(depth_diff)[::-1][:10]
    return {
        "site": site,
        "annual_depth_pct": pct(all_jd, all_fd),
        "annual_swe_pct": pct(jswe[valid], fswe[valid]),
        "physical_pack_days": int(np.sum(physical)),
        "physical_pack_depth_pct": pct(jd, fd),
        "all_day_depth_nrmse": nrmse(all_jd, all_fd),
        "effective_depth_pct": pct(jeff, feff),
        "effective_depth_nrmse": nrmse(jeff, feff),
        "reported_grid_depth_pct": pct(jgrid, fgrid),
        "reported_grid_depth_nrmse": nrmse(jgrid, fgrid),
        "physical_pack_depth_nrmse": nrmse(jd, fd),
        "physical_pack_density_pct": pct(jden, fden),
        "snow_cover_mean_difference": float(np.mean(jc) - np.mean(fc)),
        "cumulative_melt_pct": 100.0 * (float(np.sum(jm)) - float(np.sum(fm))) /
            max(abs(float(np.sum(fm))), 1e-30),
        "trace_only_days_julia": int(np.sum(valid & (jswe > 0) & (jswe <= 0.02))),
        "trace_only_days_fortran": int(np.sum(valid & (fswe > 0) & (fswe <= 0.02))),
        "largest_depth_difference_days": [
            {"date": int(days[i]), "julia_depth": float(jdepth[i]),
             "fortran_depth": float(fdepth[i]), "julia_swe": float(jswe[i]),
             "fortran_swe": float(fswe[i]), "julia_cover": float(jcover[i]),
             "fortran_cover": float(fcover[i]),
             "difference": float(jdepth[i] - fdepth[i])}
            for i in worst
        ],
    }


def bow_btran_audit(jpath, fpath):
    days, v = aligned(jpath, fpath, [
        ("btran", "BTRAN", "BTRANMN"),
        ("qvegt", "QVEGT", "QVEGT"),
        ("fctr", "FCTR", "FCTR"),
        ("latent", "EFLX_LH_TOT", "EFLX_LH_TOT"),
        ("soil_liq", "TOTSOILLIQ", "TOTSOILLIQ"),
        ("swe", "H2OSNO", "H2OSNO"),
    ])
    jb, fb = v["btran"]
    ok = np.isfinite(jb) & np.isfinite(fb)
    diff = np.abs(jb - fb)
    order = np.argsort(np.where(ok, diff, -np.inf))[::-1][:10]
    consumers = {}
    for key in ("qvegt", "fctr", "latent", "soil_liq", "swe"):
        j, f = finite_pair(v[key])
        consumers[key + "_annual_pct"] = pct(j, f)
    return {
        "annual_pct": pct(jb[ok], fb[ok]),
        "largest_absolute_difference_days": [
            {"date": int(days[i]), "julia": float(jb[i]), "fortran": float(fb[i]),
             "difference": float(jb[i] - fb[i])} for i in order
        ],
        **consumers,
    }


def main():
    bow_j = os.path.join(OUT, "Bow", "julia_clm_bow_phs_2003.nc")
    stw_j = os.path.join(OUT, "Stillwater", "julia_clm_stillwater_phs_2003.nc")
    bow_f = os.path.join(ROOT, "clm_parity_run",
                         "Bow_at_Banff_lumped.clm2.h0.2003-01-01-00000.nc")
    stw_f = os.path.join(ROOT, "domain_Stillwater_Oklahoma", "optimization", "CLM",
                         "dds_clm_dds_calibration", "final_evaluation",
                         "Stillwater_Oklahoma.clm2.h0.2003-01-01-00000.nc")
    result = {
        "bow_btran": bow_btran_audit(bow_j, bow_f),
        "snow": [snow_audit("Bow", bow_j, bow_f),
                 snow_audit("Stillwater", stw_j, stw_f)],
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
