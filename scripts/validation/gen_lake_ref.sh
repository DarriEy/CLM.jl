#!/usr/bin/env bash
# =============================================================================
# gen_lake_ref.sh — regenerate the Fortran CTSM lake-landunit h0 reference that
# scripts/fortran_parity_lake.jl diffs the Julia deep-lake port against.
#
# The original `clm_lake_run` h0 domain was removed when local SYMFLUENCE_data
# shrank (2026-07). This reconstructs it from surviving assets — NO SourceMods
# and NO rebuild are needed (it is a pure namelist reconfiguration of the
# existing instrumented cesm.exe); it is a cold-start SP run on a PCT_LAKE=100
# surfdata, which the lake physics init/step exercise independent of use_cn.
#
# Verified 2026-07-18 (rc=0, SUCCESSFUL TERMINATION on attempt 1): the h0 carries
# TG, TLAKE(×10 levlak), LAKEICEFRAC(×10), LAKEICEFRAC_SURF, EFLX_LH_TOT, FSH,
# FIRE, FSA, TSA, H2OSNO over 48 hourly records. Result vs Julia: TLAKE/LAKEICE
# thermodynamics track (~1%); the residual is the lake SURFACE TURBULENT-FLUX
# solve (FSH/EFLX_LH ~10%, TSA ~7-11%) — the documented open residual (the FSA
# ~13× blow-ups are relative-diff artifacts on the near-zero dawn/dusk solar).
#
# Base staging (a 2003 single-point SP run dir): clm_parity_run — this survives
# only on the Drive migration copy
#   /Users/.../GoogleDrive-*/My Drive/data/SYMFLUENCE_data/clm_parity_run
# (rsync it back to compHydro/SYMFLUENCE_data/ if you also need the SP :parity
# oracle dumps locally).
# =============================================================================
set -euo pipefail

CASE="${CASE:-/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/installs/clm/cases/symfluence_build}"
BASE="${BASE:-/Users/darri.eythorsson/Library/CloudStorage/GoogleDrive-dareyt@gmail.com/My Drive/data/SYMFLUENCE_data/clm_parity_run}"
LK="${LK:-/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_lake_run}"
# CTSM-format PCT_LAKE=100 surfdata lives in the Julia repo test inputs:
LAKESURF="${LAKESURF:-$(cd "$(dirname "$0")/../.." && pwd)/test_inputs/lake/surfdata_lake100.nc}"

# --- stage the run dir from the 2003 SP base ---
rm -rf "$LK"; mkdir -p "$LK/timing/checkpoints" "$LK/init_generated_files"
for f in drv_in datm_in datm.streams.xml drv_flds_in nuopc.runconfig nuopc.runseq lnd_in fd.yaml \
         Bow_at_Banff_lumped.cpl.r.2003-01-01-00000.nc Bow_at_Banff_lumped.datm.r.2003-01-01-00000.nc; do
  cp -f "$BASE/$f" "$LK/"
done
cp -f "$LAKESURF" "$LK/surfdata_lake100.nc"

# --- lnd_in: lake surfdata + cold start + hourly lake history (48 samples) ---
python3 - "$LK" <<'PY'
import sys, re
LK = sys.argv[1]
p = LK + "/lnd_in"; s = open(p).read()
s = re.sub(r"fsurdat\s*=\s*'[^']*'", "fsurdat = '%s/surfdata_lake100.nc'" % LK, s)
s = re.sub(r"nrevsn\s*=\s*'[^']*'",  "nrevsn = ''", s)      # startup: no branch restart
s = re.sub(r"finidat\s*=\s*'[^']*'", "finidat = ''", s)     # cold start
s = re.sub(r"hist_fincl1\s*=.*",
           "hist_fincl1 = 'TG','TLAKE','LAKEICEFRAC','LAKEICEFRAC_SURF','EFLX_LH_TOT','FSH','FIRE','FSA','TSA','H2OSNO'", s)
s = re.sub(r"hist_nhtfrq\s*=.*", "hist_nhtfrq = 1", s)      # hourly
s = re.sub(r"hist_mfilt\s*=.*",  "hist_mfilt = 48", s)      # 48 records / h0 file
open(p, "w").write(s)
p = LK + "/nuopc.runconfig"; s = open(p).read()
s = re.sub(r"start_type\s*=\s*\w+",     "start_type = startup", s)
s = re.sub(r"start_ymd\s*=\s*\d+",      "start_ymd = 20030101", s)
s = re.sub(r"stop_option\s*=\s*\w+",    "stop_option = nsteps", s)
s = re.sub(r"stop_n\s*=\s*-?\d+",       "stop_n = 48", s)
s = re.sub(r"restart_option\s*=\s*\w+", "restart_option = nsteps", s)
s = re.sub(r"restart_n\s*=\s*-?\d+",    "restart_n = 48", s)
open(p, "w").write(s)
PY

# --- codesign a local copy of the exe (macOS 26 ad-hoc-sign SIGTRAP fix) + run ---
cp -f "$CASE/bld/cesm.exe" "$LK/cesm_lake.exe"; codesign -f -s - "$LK/cesm_lake.exe"
cd "$LK"; source "$CASE/.env_mach_specific.sh"
export GFORTRAN_ERROR_BACKTRACE=1 HWLOC_COMPONENTS=-opencl MallocNanoZone=0
for a in $(seq 1 25); do        # retry: xzone brk SIGTRAP is flaky at init
  printf './Bow_at_Banff_lumped.cpl.r.2003-01-01-00000.nc'  > rpointer.cpl
  printf './Bow_at_Banff_lumped.datm.r.2003-01-01-00000.nc' > rpointer.atm
  "$LK/cesm_lake.exe" > "$LK/cesm_lake.log" 2>&1 && { echo "SUCCESS attempt $a"; break; }
  echo "attempt $a failed (rc=$?), retrying"
done
ls "$LK"/*.clm2.h0.2003-01-01-00000.nc && echo ">> lake h0 reference regenerated"
