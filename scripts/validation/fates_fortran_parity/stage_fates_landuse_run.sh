#!/usr/bin/env bash
# =============================================================================
# stage_fates_landuse_run.sh — Fortran references for FATES land use: logging,
# LUH2 land-use transitions, LUH2 harvest, land-use x PFT (lupft), potential veg.
#
# Same single-point case as stage_fates_modes_run.sh (Bow at Banff, archived
# fates_parity_1pt namelists, uninstrumented cesm.exe, 2002-07-01, dtime 3600 s,
# hourly h0), plus the namelist overrides of one MODE:
#
#   MODE=logging   fates_harvest_mode='event_code' (FATES logging, no LUH data)
#   MODE=luh       use_fates_luh (LUH2 states + transitions), no harvest
#   MODE=luh_area  use_fates_luh + fates_harvest_mode='luhdata_area'
#   MODE=luh_mass  use_fates_luh + fates_harvest_mode='luhdata_mass'
#   MODE=lupft     use_fates_lupft (=> luh + nocomp + fixed biogeography, flandusepftdat)
#   MODE=potveg    use_fates_potentialveg (=> luh, no LUH2 file)
#
# The LUH2 files are the single-point subsets made by make_luh_1pt.jl (nearest 4x5
# cell of the CESM inputdata LUH2 timeseries / landuse x pft map):
#   LUH_FILE    (default data/fates/luh_1pt/LUH2_1pt_Bow_hist_simyr1850-2015.nc)
#   LUPFT_FILE  (default data/fates/luh_1pt/fates_landuse_pft_map_1pt_Bow.nc)
#
# The FATES parameter file is the default .cdl with the logging event moved into the
# 10-day window and the direct-logging size threshold lowered so the cold-start
# saplings are harvested (the default dbhmin = 50 cm logs nothing in 10 days):
#   LOG_EVENT_CODE (default 20020704: a one-off YYYYMMDD event on day 4)
#   LOG_DBHMIN     (default 0.1 cm)
# The same edited .cdl is written to $RUN/fates_params_landuse.cdl for the port.
#
#   run dir : $RUN (default $SYMFLUENCE_DATA/runs/fates_lu_<MODE>_1pt)
#
# then score the port:
#   julia --project=. scripts/validation/fates_fortran_parity/fates_landuse_parity.jl \
#         <MODE> $RUN/fates_lu_<MODE>_1pt.clm2.h0.2002-07-01-00000.nc 240
# =============================================================================
set -euo pipefail
MODE=${MODE:?set MODE=logging|luh|luh_area|luh_mass|lupft|potveg}
D=${SYMFLUENCE_DATA:?set SYMFLUENCE_DATA}
IN=${CESM_INPUTDATA:-$D/installs/cesm-inputdata}
EXE=${EXE:-$D/installs/clm/bin/cesm.exe}
SITE=${SITE:-$D/domain_Bow_at_Banff_lumped_era5/settings/CLM}
ARCH=${ARCH:-$D/installs/clm/case-archive/fates_parity_1pt/CaseDocs}
CASE_NAME=fates_lu_${MODE}_1pt
RUN=${RUN:-$D/runs/$CASE_NAME}
STEPS=${STEPS:-240}
START_YMD=${START_YMD:-20020701}
CDL=${CDL:-$D/installs/clm/ctsm5.3.012/src/fates/parameter_files/fates_params_default.cdl}
LUH_FILE=${LUH_FILE:-$(cd "$(dirname "$0")/../../.." && pwd)/data/fates/luh_1pt/LUH2_1pt_Bow_hist_simyr1850-2015.nc}
LUPFT_FILE=${LUPFT_FILE:-$(cd "$(dirname "$0")/../../.." && pwd)/data/fates/luh_1pt/fates_landuse_pft_map_1pt_Bow.nc}
LOG_EVENT_CODE=${LOG_EVENT_CODE:-20020704}
LOG_DBHMIN=${LOG_DBHMIN:-0.1}
OLD_D=/Users/darri.eythorsson/compHydro/SYMFLUENCE_data
OLD_IN=/Users/darri.eythorsson/projects/cesm-inputdata

case "$MODE" in
  logging)  FLAGS="fates_harvest_mode='event_code'" ;;
  luh)      FLAGS="use_fates_luh=.true. fluh_timeseries='$LUH_FILE'" ;;
  luh_area) FLAGS="use_fates_luh=.true. fluh_timeseries='$LUH_FILE' fates_harvest_mode='luhdata_area'" ;;
  luh_mass) FLAGS="use_fates_luh=.true. fluh_timeseries='$LUH_FILE' fates_harvest_mode='luhdata_mass'" ;;
  lupft)    FLAGS="use_fates_lupft=.true. use_fates_luh=.true. use_fates_nocomp=.true. use_fates_fixed_biogeog=.true. fluh_timeseries='$LUH_FILE' flandusepftdat='$LUPFT_FILE'" ;;
  potveg)   FLAGS="use_fates_potentialveg=.true. use_fates_luh=.true." ;;
  *) echo "unknown MODE=$MODE"; exit 1 ;;
esac

rm -rf "$RUN"; mkdir -p "$RUN/timing/checkpoints"
for f in nuopc.runseq fd.yaml datm_in datm.streams.xml drv_in; do
  sed -e "s#$OLD_D#$D#g" -e "s#$OLD_IN#$IN#g" "$SITE/$f" > "$RUN/$f"
done
cp "$ARCH/drv_flds_in" "$RUN/"
sed -e "s#^ fates_landuse_logging_event_code = .*# fates_landuse_logging_event_code = $LOG_EVENT_CODE ;#" \
    -e "s#^ fates_landuse_logging_dbhmin = .*# fates_landuse_logging_dbhmin = $LOG_DBHMIN ;#" \
    "$CDL" > "$RUN/fates_params_landuse.cdl"
grep -q "fates_landuse_logging_event_code = $LOG_EVENT_CODE ;" "$RUN/fates_params_landuse.cdl"
grep -q "fates_landuse_logging_dbhmin = $LOG_DBHMIN ;" "$RUN/fates_params_landuse.cdl"
ncgen -o "$RUN/fates_params_landuse.nc" "$RUN/fates_params_landuse.cdl"

python3 - "$SITE" "$ARCH" "$RUN" "$D" "$IN" "$START_YMD" "$STEPS" "$OLD_D" "$OLD_IN" \
          "$CASE_NAME" $FLAGS <<'PY'
import sys, re
SITE, ARCH, RUN, D, IN, start_ymd, steps, OLD_D, OLD_IN, case = sys.argv[1:11]
flags = dict(a.split('=', 1) for a in sys.argv[11:])
fix = lambda line: line.replace(OLD_D, D).replace(OLD_IN, IN)

repl = {'case_name': case, 'case_desc': 'FATES land-use parity',
        'start_ymd': start_ymd, 'stop_option': 'nsteps', 'stop_n': steps,
        'restart_option': 'never', 'diro': RUN}
out = []
for line in open(f'{SITE}/nuopc.runconfig'):
    line = fix(line)
    m = re.match(r'^(\s*)(\w+)( = )(.*)$', line.rstrip('\n'))
    out.append(f"{m.group(1)}{m.group(2)}{m.group(3)}{repl[m.group(2)]}\n"
               if m and m.group(2) in repl else line)
open(f'{RUN}/nuopc.runconfig', 'w').writelines(out)

P = f'{SITE}/parameters'
over = {'fsurdat': f"'{P}/surfdata_clm.nc'", 'paramfile': f"'{P}/clm5_params.nc'",
        'fates_paramfile': f"'{RUN}/fates_params_landuse.nc'", 'use_lch4': '.false.',
        **flags}
FIELDS = ['TLAI', 'FSH', 'EFLX_LH_TOT', 'FATES_GPP', 'FATES_LAI', 'FATES_NCOHORTS',
          'FATES_NPATCHES', 'FATES_VEGC', 'FATES_NPATCHES_SECONDARY',
          'FATES_SECONDARY_FOREST_FRACTION', 'FATES_SECONDARY_FOREST_VEGC',
          'FATES_DISTURBANCE_RATE_LOGGING', 'FATES_DISTURBANCE_RATE_TREEFALL',
          'FATES_HARVEST_WOODPROD_C_FLUX', 'FATES_LUCHANGE_WOODPROD_C_FLUX',
          'FATES_WOOD_PRODUCT', 'FATES_HARVEST_DEBT', 'FATES_LITTER_IN',
          'FATES_PATCHAREA_LU', 'FATES_TRANSITION_MATRIX_LULU',
          'FATES_DISTURBANCE_RATE_MATRIX_LULU', 'FATES_NPLANT_PF', 'FATES_LEAFC_PF']
seen, out = set(), []
for line in open(f'{ARCH}/lnd_in'):
    line = fix(line)
    m = re.match(r'^(\s*)(\w+)(\s*=\s*)(.*)$', line.rstrip('\n'))
    if m and m.group(2) in over:
        k = m.group(2); seen.add(k)
        out.append(f"{m.group(1)}{k}{m.group(3)}{over[k]}\n")
    else:
        out.append(line)
    if line.strip() == '&clm_inparm':
        out.append(" hist_nhtfrq = 1\n hist_mfilt = 1000\n hist_fincl1 = " +
                   ",".join(f"'{f}'" for f in FIELDS) + "\n")
        # file names absent from the archived namelist (empty there) go in here
        for k in ('fluh_timeseries', 'flandusepftdat'):
            if k in over:
                out.append(f" {k} = {over[k]}\n"); seen.add(k)
open(f'{RUN}/lnd_in', 'w').writelines(out)
assert not (set(over) - seen), set(over) - seen
print('[ok] namelists staged in', RUN, flags)
PY

# Same launch environment and retry loop as setup_case.sh (macOS 26 xzone/ESMF trap).
cd "$RUN"
for attempt in $(seq 1 30); do
  rm -f "$CASE_NAME".clm2.h0.*.nc
  if MallocNanoZone=0 HWLOC_COMPONENTS=-opencl "$EXE" > cesm.stdout 2>&1 \
     && ls "$CASE_NAME".clm2.h0.*.nc > /dev/null 2>&1; then
    echo "[ok] run succeeded on attempt $attempt"; break
  fi
  if grep -q "ENDRUN" lnd.log 2>/dev/null; then grep -B8 -m1 ENDRUN lnd.log; exit 1; fi
  echo "[..] attempt $attempt did not complete; retrying"
done
ls "$RUN"/"$CASE_NAME".clm2.h0.*.nc
