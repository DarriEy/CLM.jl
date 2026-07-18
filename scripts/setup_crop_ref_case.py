#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Build the CTSM CN-crop reference case for CLM.jl crop/irrigation parity
(backlog A3 crop pools + crop N, A4 irrigation, A2 fire ``baf_crop``).

This is the *reproducible* generator for the Fortran ground-truth case that
``docs/CROP_PARITY.md`` describes.  Run it on the box with the instrumented
CTSM build; it writes a self-contained run directory (namelists + forcing +
run script) that can be executed directly, in the style of the existing
``clm_bgc_spinup/bow_ref_*`` references.

Why a self-contained run dir and not a CIME case
------------------------------------------------
The working CN references on this box (``bow_ref_ch4ebul``, ``bgc_ref_summer``,
``merbleue_ref_ch4``) are all run this way: SYMFLUENCE writes ``lnd_in`` /
``drv_in`` / ``datm_in`` / ``nuopc.*`` directly and invokes the shared
``cesm.exe`` in place, bypassing ``bldnml`` and ``case.submit`` entirely.  This
script mirrors that pattern so the crop case is consistent with every other
reference the parity harnesses consume.

NO REBUILD IS REQUIRED — verified, see docs/CROP_PARITY.md
----------------------------------------------------------
The shared ``cesm.exe`` (configured ``CLM_BLDNML_OPTS="-bgc sp"``) already
contains the entire crop code path: ``CropType`` (15 symbols),
``irrigationMod`` (33), ``CNPhenologyMod::CropPhenology``/``CropPhase``, and
``CNNDynamicsMod::CNNFert``/``CNSoyfix``.  ``-bgc`` only steers ``bldnml``'s
*namelist defaults*; CTSM has no CPP gate on crop (``use_crop`` is branched on
at runtime, e.g. ``clm_varpar.F90:329``).  So flipping the namelist flags is
sufficient.

Usage
-----
    python3 scripts/setup_crop_ref_case.py            # build the case
    python3 scripts/setup_crop_ref_case.py --years 6  # longer spin-up
    bash <casedir>/run_crop_ref.sh                    # then run it
"""
from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

SYM = Path("/Users/darri.eythorsson/compHydro/SYMFLUENCE_data")
DRIVE = Path(
    "/Users/darri.eythorsson/Library/CloudStorage/GoogleDrive-dareyt@gmail.com/My Drive"
)

# The crop-CFT single-point assets built by scripts/build_crop_cft_surfdata.py (PR #249).
#
# They live on the Drive, whose path contains a SPACE ("My Drive"). `nuopc.runconfig`
# fields are UNQUOTED (`mesh_lnd = /path/...`), so a space silently truncates the
# mesh path and the run would fall back to a different grid -- the exact harness
# mismatch class behind #233/#240/#243/#248/#251. So the assets are staged into a
# space-free directory on the Fortran box, and the case only ever references that.
CROP_ASSETS_SRC = DRIVE / "data/SYMFLUENCE_data/crop_cft_surfdata"
CROP_ASSETS = SYM / "crop_cft_surfdata"
FSURDAT = CROP_ASSETS / "surfdata_cropCFT_USplains_1pt.nc"
MESH = CROP_ASSETS / "esmf_mesh_croppt.nc"

# Donor CN case: a known-good, currently-passing CN reference on this box.
DONOR = SYM / "clm_bgc_spinup/bow_ref_ch4ebul"
# Its paramfile carries pft=79 and every crop parameter (gddmin/hybgdd/mxmat/
# baset/fertnitro/...), so it covers the crop CFT layout unchanged.
PARAMFILE = (
    SYM
    / "domain_Bow_at_Banff_lumped/optimization/CLM/dds_run_1/final_evaluation"
    / "settings/CLM/parameters/clm5_params.nc"
)
DONOR_FORCING = SYM / "domain_Peatland_MerBleue_Canada/data/forcing/CLM_input"

CASE = SYM / "clm_bgc_spinup/crop_ref_usplains"

CROP_LAT, CROP_LON = 44.80, -96.71
FORCING_YEARS = (2016, 2017)

EXE = SYM / "installs/clm/cases/symfluence_build/bld/cesm.exe"

# --- lnd_in: the flags that define a CN-crop cold start -----------------------
# Every one of these is set EXPLICITLY. An unset flag is not a matched flag --
# that was the #251 bug (use_bedrock defaulted differently than assumed), so the
# case pins each side of the Julia/Fortran comparison rather than inheriting it.
LND_FLAGS = {
    # the crop path itself
    "use_crop": ".true.",
    "use_cn": ".true.",
    "irrigate": ".true.",
    "use_fertilizer": ".true.",
    "use_grainproduct": ".true.",
    "create_crop_landunit": ".true.",
    # cold start -- no finidat; the crop lifecycle must establish itself
    "finidat": "''",
    # CN configuration, matched to the working Bow CN references
    "use_nitrif_denitrif": ".true.",
    "use_fun": ".true.",
    "soil_decomp_method": "'CENTURYKoven2013'",
    # keep methane OFF so the diff isolates the crop path
    "use_lch4": ".false.",
    "use_fates": ".false.",
    "use_luna": ".true.",
    "use_hydrstress": ".true.",
    # crop-CFT assets
    "fsurdat": f"'{FSURDAT}'",
    "paramfile": f"'{PARAMFILE}'",
    # annual history; the parity diff reads the pdump/BGCDUMP, not h0
    "hist_nhtfrq": "-8760",
    "hist_mfilt": "20",
}

# Crop / irrigation / crop-fire diagnostics on the h0 stream, so the spin-up can
# be checked for "did a crop actually plant, grow and harvest" without parsing
# the dumps.
HIST_FIELDS = [
    "GDD0", "GDD8", "GDD10", "GDD020", "GDD820", "GDD1020",
    "GDDPLANT", "GDDHARV", "HUI", "CPHASE", "CROPPROD1C", "GRAINC",
    "GRAINN", "GRAINC_TO_FOOD", "QIRRIG", "TOTVEGC", "TOTVEGN",
    "NFERTILIZATION", "LEAFC", "TLAI", "BAF_CROP", "NFIRE",
]


def sh(cmd: list[str]) -> None:
    subprocess.run(cmd, check=True)


def patch_namelist(path: Path, flags: dict[str, str]) -> None:
    """Set each key in a Fortran namelist, appending to group 1 if absent."""
    text = path.read_text()
    for key, val in flags.items():
        pat = re.compile(rf"^\s*{re.escape(key)}\s*=.*$", re.MULTILINE)
        if pat.search(text):
            text = pat.sub(f" {key} = {val}", text, count=1)
        else:
            # append inside the first namelist group
            text = re.sub(r"\n/", f"\n {key} = {val}\n/", text, count=1)
    path.write_text(text)


def replace_paths(path: Path, old_sub: str, new: str) -> int:
    """Repoint every path containing ``old_sub`` at ``new`` (mesh redirection).

    The delimiter class must exclude ``<`` and ``>``: in ``datm.streams.xml``
    the path sits inside ``<meshfile>...</meshfile>`` with no surrounding
    whitespace or quotes, so a class of only ``[^\\s'"]`` swallows the tags
    themselves and silently deletes the element -- which CTSM reports much
    later as the opaque "ERROR: mesh file name must be provided".
    """
    text = path.read_text()
    text, n = re.subn(rf"[^\s'\"<>]*{re.escape(old_sub)}[^\s'\"<>]*", new, text)
    path.write_text(text)
    return n


def verify_case(case: Path) -> int:
    """Assert the generated case is structurally sound. Returns 0 if OK."""
    errs: list[str] = []

    lnd = (case / "lnd_in").read_text()
    for key, want in LND_FLAGS.items():
        m = re.search(rf"^\s*{re.escape(key)}\s*=\s*(.+?)\s*$", lnd, re.MULTILINE)
        if m is None:
            errs.append(f"lnd_in: {key} absent (an unset flag is not a matched flag)")
        elif m.group(1) != want:
            errs.append(f"lnd_in: {key} = {m.group(1)!r}, expected {want!r}")

    streams = (case / "datm.streams.xml").read_text()
    n_streams = len(re.findall(r"<stream_info ", streams))
    meshes = re.findall(r"<meshfile>([^<]*)</meshfile>", streams)
    if len(meshes) != n_streams:
        errs.append(
            f"datm.streams.xml: {len(meshes)} <meshfile> for {n_streams} streams "
            "(every stream needs one; CTSM reports this as 'mesh file name must be provided')"
        )
    for mf in meshes:
        if mf != str(MESH):
            errs.append(f"datm.streams.xml: meshfile {mf!r} is not the crop mesh")

    forcing_files = re.findall(r"<file>([^<]*clmforc[^<]*)</file>", streams)
    if not forcing_files:
        errs.append("datm.streams.xml: no clmforc data files")
    for f in set(forcing_files):
        if not Path(f).exists():
            errs.append(f"datm.streams.xml: forcing file missing on disk: {f}")
    # The forcing years must be exactly the cycled window -- a year/hour
    # mismatch between harness and reference was the #233 root cause.
    years = sorted({int(re.search(r"clmforc\.(\d{4})\.nc", f).group(1)) for f in forcing_files})
    if years != sorted(FORCING_YEARS):
        errs.append(f"datm.streams.xml: forcing years {years} != {sorted(FORCING_YEARS)}")
    for tag, want in (("year_first", FORCING_YEARS[0]), ("year_last", FORCING_YEARS[-1])):
        got = {int(v) for v in re.findall(rf"<{tag}>(\d+)</{tag}>", streams)}
        if got != {want}:
            errs.append(f"datm.streams.xml: {tag} = {got}, expected {{{want}}}")

    for f in re.findall(r"<file>([^<]*topo_forcing[^<]*)</file>", streams):
        if "MerBleue" not in f:
            errs.append(
                f"datm.streams.xml: topo stream {f!r} is not the met donor's -- "
                "the lapse-rate reference elevation would be wrong"
            )

    rc_text = (case / "nuopc.runconfig").read_text()
    for key, want in (("mesh_lnd", str(MESH)), ("mesh_atm", str(MESH)), ("mesh_mask", str(MESH))):
        m = re.search(rf"^\s*{key}\s*=\s*(.+?)\s*$", rc_text, re.MULTILINE)
        if m is None or m.group(1) != want:
            errs.append(f"nuopc.runconfig: {key} = {m and m.group(1)!r}, expected {want!r}")
    if re.search(r"^\s*restart_ymd\s*=\s*-999\s*$", rc_text, re.MULTILINE) is None:
        errs.append("nuopc.runconfig: restart_ymd sentinel -999 was clobbered")

    for f in ("fd.yaml", "lnd_in", "drv_in", "datm_in", "nuopc.runconfig",
              "nuopc.runseq", "drv_flds_in", "datm.streams.xml"):
        if not (case / f).exists():
            errs.append(f"missing case file: {f}")

    if errs:
        print("\n>> CASE VERIFICATION FAILED:", file=sys.stderr)
        for e in errs:
            print(f"   - {e}", file=sys.stderr)
        return 1
    print(">> case verified: flags, meshes, forcing years, topo donor, run control")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--years",
        type=int,
        default=4,
        help="model years of cold-start spin-up (crop must plant/grow/harvest)",
    )
    ap.add_argument("--case", type=Path, default=CASE)
    ap.add_argument(
        "--real-era5",
        action="store_true",
        help="pull site-native ERA5 instead of relabeling the donor forcing "
        "(slow: ARCO chunks one GLOBAL field per hour+variable)",
    )
    args = ap.parse_args()
    case = args.case

    # Stage the crop assets off the space-containing Drive path (see CROP_ASSETS).
    CROP_ASSETS.mkdir(parents=True, exist_ok=True)
    for name in (FSURDAT.name, MESH.name):
        src = CROP_ASSETS_SRC / name
        if not src.exists():
            print(f"ERROR: missing crop asset: {src}", file=sys.stderr)
            return 1
        shutil.copy(src, CROP_ASSETS / name)
    print(f">> staged crop assets -> {CROP_ASSETS} (space-free path)")

    for p in (FSURDAT, MESH, DONOR, PARAMFILE, EXE):
        if not p.exists():
            print(f"ERROR: missing required input: {p}", file=sys.stderr)
            return 1
    assert " " not in str(MESH), "mesh path must not contain spaces (nuopc.runconfig is unquoted)"

    case.mkdir(parents=True, exist_ok=True)
    print(f">> case: {case}")

    # 1. namelist scaffold from the known-good CN donor -----------------------
    # fd.yaml is the NUOPC field dictionary -- omitting it aborts in
    # ESMF_IO_YAMLRead before the model ever starts, with a bare
    # "Caught exception reading content from file".
    for f in (
        "lnd_in", "drv_in", "datm_in", "datm.streams.xml",
        "nuopc.runconfig", "nuopc.runseq", "drv_flds_in", "fd.yaml",
    ):
        shutil.copy(DONOR / f, case / f)
    print(f">> scaffolded namelists from {DONOR.name}")

    # 2. forcing at the crop point --------------------------------------------
    forcing = case / "forcing"
    here = Path(__file__).parent
    if args.real_era5:
        sh([sys.executable, str(here / "acquire_crop_forcing.py"),
            "--lat", str(CROP_LAT), "--lon", str(CROP_LON),
            "--years", *[str(y) for y in FORCING_YEARS],
            "--outdir", str(forcing)])
    else:
        sh([sys.executable, str(here / "acquire_crop_forcing.py"),
            "--relabel", *[str(DONOR_FORCING / f"clmforc.{y}.nc") for y in FORCING_YEARS],
            "--outdir", str(forcing),
            "--lat", str(CROP_LAT), "--lon", str(CROP_LON)])
    print(f">> forcing in {forcing}")

    # 3. crop surfdata + the single-point crop mesh everywhere ----------------
    # Mesh appears in FOUR places; missing any one silently leaves the run on
    # the donor's grid, which is exactly the class of harness mismatch that
    # produced #233/#240/#243/#248/#251.
    for f in ("nuopc.runconfig", "datm_in", "datm.streams.xml"):
        replace_paths(case / f, "esmf_mesh.nc", str(MESH))
    patch_namelist(case / "lnd_in", LND_FLAGS)
    patch_namelist(
        case / "lnd_in",
        {"hist_fincl1": ",".join(f"'{v}'" for v in HIST_FIELDS)},
    )

    # 4. point the datm streams at our forcing, on the right YEARS ------------
    # taxmode=cycle repeats [year_first, year_last] for the whole spin-up.
    s = (case / "datm.streams.xml").read_text()
    s = re.sub(r"<file>[^<]*clmforc[^<]*</file>", "", s)
    s = re.sub(
        r"<datafiles>\s*",
        "<datafiles>\n"
        + "\n".join(f"      <file>{forcing}/clmforc.{y}.nc</file>" for y in FORCING_YEARS)
        + "\n   ",
        s,
    )
    for tag, val in (
        ("year_first", FORCING_YEARS[0]),
        ("year_last", FORCING_YEARS[-1]),
        ("year_align", FORCING_YEARS[0]),
    ):
        s = re.sub(rf"<{tag}>\d+</{tag}>", f"<{tag}>{val}</{tag}>", s)
    # The topo stream supplies the elevation the atm forcing is lapse-rate
    # downscaled FROM. Left at the donor CN case's value it would be Bow's
    # ~1400 m against a ~500 m plains surfdata -- a multi-K temperature shift,
    # straight into GDD accumulation and the crop planting decision. It must
    # match whichever site the met actually came from.
    s = re.sub(
        r"<file>[^<]*topo_forcing\.nc</file>",
        f"<file>{DONOR_FORCING / 'topo_forcing.nc'}</file>",
        s,
    )
    (case / "datm.streams.xml").write_text(s)

    # 5. run control: cold start at the first forcing year --------------------
    # Each key is anchored to start-of-line + whitespace. Without the anchor,
    # `start_ymd` also matches `restart_ymd` and silently clobbers its -999
    # sentinel -- caught by verifying the generated namelist rather than
    # trusting the edit.
    rc_keys = {
        "case_name": "Crop_USplains",
        "start_type": "startup",
        "start_ymd": f"{FORCING_YEARS[0]}0101",
        "stop_option": "nyears",
        "stop_n": str(args.years),
        "restart_option": "nyears",
        "restart_n": "1",
        "restart_ymd": "-999",
    }
    rc = (case / "nuopc.runconfig").read_text()
    for key, val in rc_keys.items():
        rc, n = re.subn(
            rf"^(\s*){re.escape(key)}\s*=.*$", rf"\g<1>{key} = {val}", rc, flags=re.MULTILINE
        )
        if n != 1:
            print(f"ERROR: nuopc.runconfig key {key!r} matched {n} times (expected 1)",
                  file=sys.stderr)
            return 1
    (case / "nuopc.runconfig").write_text(rc)

    # 6. the run script -------------------------------------------------------
    # macOS-26 notes (memory: macos26-lldb-xzone-brk-artifact): run cesm.exe
    # PLAINLY -- the xzone allocator's inline `brk #0x1` guards halt lldb but
    # run fine unwrapped. MallocNanoZone=0 and the HWLOC_COMPONENTS=-opencl
    # workaround (Metal SIGILL at MPI startup) are both required.
    run_sh = case / "run_crop_ref.sh"
    run_sh.write_text(f"""#!/bin/bash
# CN-CROP reference spin-up -- CLM.jl backlog A3 (crop pools + crop N),
# A4 (irrigation), A2 (fire baf_crop). Generated by scripts/setup_crop_ref_case.py.
#
# Cold start ({args.years} model years) on the crop-CFT single-point surfdata, cycling
# {FORCING_YEARS[0]}-{FORCING_YEARS[-1]} forcing, so a crop plants, accumulates GDD/HUI, grows and
# harvests. NO REBUILD NEEDED: the shared cesm.exe already contains CropType,
# irrigationMod, CropPhenology and CNNFert/CNSoyfix (see docs/CROP_PARITY.md).
cd "{case}"
export NETCDF_PATH=/usr/local
export ESMFMKFILE=/Users/darri.eythorsson/.local/esmf/lib/libO/Darwin.gfortranclang.64.mpiuni.default/esmf.mk
export NETCDF_C_PATH=/opt/homebrew/Cellar/netcdf/4.10.0
export NETCDF_FORTRAN_PATH=/opt/homebrew/Cellar/netcdf-fortran/4.6.2
# Metal SIGILL at MPI startup comes in via libhwloc->OpenCL; disable that component.
export HWLOC_COMPONENTS=-opencl
export MallocNanoZone=0
# Dump windows OFF for the spin-up (they would emit a file per step for years).
# Set these to the growing-season nstep bracket for the parity diff run.
export PDUMP_NSTEP_LO=${{PDUMP_NSTEP_LO:-999999999}}
export PDUMP_NSTEP_HI=${{PDUMP_NSTEP_HI:-0}}
export BGCDUMP_NSTEP_LO=${{BGCDUMP_NSTEP_LO:-999999999}}
export BGCDUMP_NSTEP_HI=${{BGCDUMP_NSTEP_HI:-0}}
mkdir -p init_generated_files
EXE="{EXE}"
# macOS-26 xzone brk guard fires nondeterministically on a plain run; retry.
ok=0
for attempt in 1 2 3 4 5; do
  echo "=== attempt $attempt ===" >> cesm_run.log
  "$EXE" >> cesm_run.log 2>&1
  rc=$?
  echo "EXIT_CODE=$rc (attempt $attempt)" >> cesm_run.log
  if [ $rc -eq 0 ]; then ok=1; break; fi
done
echo "FINAL_OK=$ok"
""")
    run_sh.chmod(0o755)

    # 7. VERIFY the generated case ------------------------------------------
    # The generator is not trusted; every field it wrote is read back. Three
    # silent malformations have already been caught this way (space-truncated
    # mesh path, start_ymd clobbering restart_ymd, a wrong-site topo stream),
    # plus an XML regex that deleted the <meshfile> elements outright. A
    # malformed case does not fail loudly -- it either aborts with an opaque
    # ESMF error or, worse, runs to completion on the wrong inputs.
    if (rc := verify_case(case)) != 0:
        return rc

    print(f">> wrote {run_sh}")
    print(f">> flags: " + " ".join(f"{k}={v}" for k, v in LND_FLAGS.items()
                                   if k.startswith(("use_", "irrigate", "create_"))))
    print(f">> run it:  bash {run_sh}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
