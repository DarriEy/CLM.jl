#!/usr/bin/env python3
"""
csfs_connector_to_obs.py — fetch a gauge's daily discharge via a live CSFS
connector (one that serves history, e.g. spain_cedex, sweden_smhi) and write a
Symfluence-style preprocessed streamflow CSV (datetime,discharge_cms) into a
domain's observations dir, so scripts/run_clm_streamflow.jl can score against it.

Complements csfs_camels_to_obs.py (which reads frozen CAMELS dataset files);
this one drives the connector API for connectors that stream/serve a date range.

Usage:
  python3 csfs_connector_to_obs.py --provider spain_cedex --station 3070 \
      --domain Mediterranean_Tagus_Spain --start 2007-01-01 --end 2019-12-31
"""
import argparse, asyncio, os, sys
from datetime import datetime, timezone

SYMFLUENCE_DATA = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data"


async def fetch(provider, station, start, end):
    from csfs.core.registry import get_connector, discover
    discover()
    async with get_connector(provider)() as conn:
        sts = await conn.fetch_stations()
        hits = [s for s in sts if str(station) in str(s.native_id)]
        nid = hits[0].native_id if hits else f"{provider}:{station}"
        name = hits[0].name if hits else station
        chunk = await conn.fetch_observations(nid, start, end)
        rows = [(o.timestamp, o.discharge_m3s) for o in chunk.observations
                if o.discharge_m3s is not None and o.discharge_m3s == o.discharge_m3s
                and o.discharge_m3s >= 0]
        return name, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--provider", required=True)
    ap.add_argument("--station", required=True)
    ap.add_argument("--domain", required=True)
    ap.add_argument("--start", default="2007-01-01")
    ap.add_argument("--end", default="2019-12-31")
    a = ap.parse_args()
    s = datetime.fromisoformat(a.start).replace(tzinfo=timezone.utc)
    e = datetime.fromisoformat(a.end).replace(tzinfo=timezone.utc)

    name, rows = asyncio.run(fetch(a.provider, a.station, s, e))
    if not rows:
        sys.exit(f"no discharge rows for {a.provider} station {a.station}")
    # collapse to daily means (defensive — most of these are already daily)
    daily = {}
    for ts, q in rows:
        daily.setdefault(str(ts)[:10], []).append(q)
    out_rows = sorted((d, sum(v) / len(v)) for d, v in daily.items())

    outdir = os.path.join(SYMFLUENCE_DATA, f"domain_{a.domain}",
                          "data", "observations", "streamflow", "preprocessed")
    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, f"{a.domain}_streamflow_processed.csv")
    with open(out, "w") as fh:
        fh.write("datetime,discharge_cms\n")
        for d, q in out_rows:
            fh.write(f"{d} 00:00:00,{q:.4f}\n")
    print(f"[{a.provider}:{a.station} {name}] wrote {len(out_rows)} rows "
          f"({out_rows[0][0]}..{out_rows[-1][0]}) -> {out}")


if __name__ == "__main__":
    main()
