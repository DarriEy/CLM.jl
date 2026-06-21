#!/usr/bin/env python3
"""
csfs_camels_to_obs.py — extract a single gauge's daily discharge from a locally
downloaded CSFS/CAMELS dataset and write it as a Symfluence-style preprocessed
streamflow CSV (datetime,discharge_cms) into a domain's observations dir, so
scripts/run_clm_streamflow.jl can score against it.

The CSFS connectors auto-download these datasets via ensure_dataset but choke on
non-UTF-8 station names; we read the (already-downloaded) files directly here.

Datasets handled:
  camels_ch  — timeseries/observation_based/CAMELS_CH_obs_based_<id>.csv
               cols: date, discharge_vol(m3/s)  (daily, 1981-2020)
  camels_se  — catchment time series/catchment_id_<id>_<NAME>.csv
  cedex      — CEDEX anuario afliq (indroea;fecha;...;caudal)  [Spain/Tagus]

Usage:
  python3 csfs_camels_to_obs.py --dataset camels_ch --gauge 2161 \
      --domain Alps_Massa_Aletsch_CH
"""
import argparse, csv, glob, os, sys

CSFS_DATASETS = "/Users/darri.eythorsson/compHydro/CSFS/data/datasets"
SYMFLUENCE_DATA = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data"


def _read_camels_ch(gauge):
    f = os.path.join(CSFS_DATASETS, "camels_ch", "camels_ch", "timeseries",
                     "observation_based", f"CAMELS_CH_obs_based_{gauge}.csv")
    if not os.path.isfile(f):
        sys.exit(f"camels_ch obs file not found: {f}")
    rows = []
    with open(f, encoding="latin-1") as fh:
        r = csv.DictReader(fh)
        qcol = next((c for c in r.fieldnames if c.startswith("discharge_vol")), None)
        for row in r:
            d, q = row.get("date"), row.get(qcol)
            if not d or q in (None, "", "NaN"):
                continue
            try:
                rows.append((d, float(q)))
            except ValueError:
                pass
    return rows


def _read_camels_se(gauge):
    base = os.path.join(CSFS_DATASETS, "camels_se")
    hits = glob.glob(os.path.join(base, "**", f"catchment_id_{gauge}_*.csv"), recursive=True)
    if not hits:
        sys.exit(f"camels_se timeseries for id {gauge} not found under {base}")
    rows = []
    with open(hits[0], encoding="latin-1") as fh:
        r = csv.DictReader(fh)
        # CAMELS-SE Qobs column is typically 'Qobs_m3s' / 'discharge'; date col 'date'
        qcol = next((c for c in r.fieldnames if "q" in c.lower() and ("m3" in c.lower() or "obs" in c.lower() or "disch" in c.lower())), None)
        dcol = next((c for c in r.fieldnames if c.lower() in ("date", "datetime", "time")), r.fieldnames[0])
        for row in r:
            d, q = row.get(dcol), row.get(qcol)
            if not d or q in (None, "", "NaN", "-9999"):
                continue
            try:
                rows.append((d[:10], float(q)))
            except ValueError:
                pass
    return rows


def _read_cedex(gauge):
    # CEDEX afliq: indroea;fecha;altura;caudal  (daily caudal m3/s)
    hits = (glob.glob(os.path.join(CSFS_DATASETS, "*cedex*", "**", "afliq*.csv"), recursive=True)
            or glob.glob(os.path.join(CSFS_DATASETS, "spain*", "**", "afliq*.csv"), recursive=True))
    if not hits:
        sys.exit("CEDEX afliq.csv not found locally — fetch via spain_cedex first")
    rows = []
    with open(hits[0], encoding="latin-1") as fh:
        r = csv.DictReader(fh, delimiter=";")
        for row in r:
            if str(row.get("indroea", "")).strip() != str(gauge):
                continue
            d, q = row.get("fecha"), row.get("caudal")
            if not d or q in (None, "", "NaN"):
                continue
            try:
                rows.append((d[:10].replace("/", "-"), float(q)))
            except ValueError:
                pass
    return rows


READERS = {"camels_ch": _read_camels_ch, "camels_se": _read_camels_se, "cedex": _read_cedex}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True, choices=list(READERS))
    ap.add_argument("--gauge", required=True)
    ap.add_argument("--domain", required=True, help="Symfluence domain name (without domain_ prefix)")
    a = ap.parse_args()

    rows = READERS[a.dataset](a.gauge)
    rows = [(d, q) for d, q in rows if q == q and q >= 0]  # drop NaN/neg
    if not rows:
        sys.exit(f"no valid discharge rows for {a.dataset} gauge {a.gauge}")
    rows.sort()
    outdir = os.path.join(SYMFLUENCE_DATA, f"domain_{a.domain}",
                          "data", "observations", "streamflow", "preprocessed")
    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, f"{a.domain}_streamflow_processed.csv")
    with open(out, "w") as fh:
        fh.write("datetime,discharge_cms\n")
        for d, q in rows:
            fh.write(f"{d} 00:00:00,{q:.4f}\n")
    print(f"wrote {len(rows)} rows ({rows[0][0]}..{rows[-1][0]}) -> {out}")


if __name__ == "__main__":
    main()
