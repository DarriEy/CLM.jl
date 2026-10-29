#!/usr/bin/env bash
# =============================================================================
# stage_spitfire_run.sh — Fortran reference for FATES SPITFIRE (fates_spitfire_mode=1).
#
# Stages and runs the archived fates_parity_1pt case (Bow at Banff single point,
# compset I2000Clm50FatesRs => use_fates, use_fates_bgc) with an UNINSTRUMENTED
# cesm.exe, fates_spitfire_mode = 1 (scalar_lightning: a constant ignition rate,
# no lightning / population-density streams), and a DAILY-MEAN h0 history of the
# fire fields plus the host fire-weather drivers. No rebuild is needed.
#
#   case config : $SYMFLUENCE_DATA/installs/clm/case-archive/fates_parity_1pt/CaseDocs/lnd_in
#   site        : $SYMFLUENCE_DATA/domain_Bow_at_Banff_lumped_era5/settings/CLM (read only)
#   run dir     : $RUN (default $SYMFLUENCE_DATA/runs/fates_fire_1pt)
#
# Same overrides as setup_case.sh (site fsurdat/paramfile, the FATES default param
# file ncgen'd from the CTSM source .cdl — byte-identical to data/fates/ — and
# use_lch4 off) plus fates_spitfire_mode = 1, start 2002-07-01, dtime 3600 s, STEPS
# steps (default 17520 = 730 days).
#
# then score the port:
#   julia --project=. scripts/validation/fates_fortran_parity/fates_spitfire_parity.jl \
#         export $RUN/fates_fire_1pt.clm2.h0.2002-07-01-00000.nc
#   julia --project=. scripts/validation/fates_fortran_parity/fates_spitfire_parity.jl
# =============================================================================
set -euo pipefail
D=${SYMFLUENCE_DATA:?set SYMFLUENCE_DATA}
IN=${CESM_INPUTDATA:-$D/installs/cesm-inputdata}
EXE=${EXE:-$D/installs/clm/bin/cesm.exe}
SITE=${SITE:-$D/domain_Bow_at_Banff_lumped_era5/settings/CLM}
ARCH=${ARCH:-$D/installs/clm/case-archive/fates_parity_1pt/CaseDocs}
RUN=${RUN:-$D/runs/fates_fire_1pt}
STEPS=${STEPS:-17520}
# SPITFIRE_MODE=2..5 adds the lightning stream and population density (read for every mode > 1, used from mode 3): the
# synthetic spatially-uniform CLM.jl parity streams in $IN/lnd/clm2/firedata
# (scripts/validation/make_fire_streams.jl), on the 0.9x1.25 ESMF mesh.
SPITFIRE_MODE=${SPITFIRE_MODE:-1}
export SPITFIRE_MODE
START_YMD=${START_YMD:-20020701}
CDL=${CDL:-$D/installs/clm/ctsm5.3.012/src/fates/parameter_files/fates_params_default.cdl}
# Paths baked into the archived namelists (the machine the case was built on).
OLD_D=/Users/darri.eythorsson/compHydro/SYMFLUENCE_data
OLD_IN=/Users/darri.eythorsson/projects/cesm-inputdata

rm -rf "$RUN"; mkdir -p "$RUN/timing/checkpoints"
for f in nuopc.runseq fd.yaml datm_in datm.streams.xml drv_in; do
  sed -e "s#$OLD_D#$D#g" -e "s#$OLD_IN#$IN#g" "$SITE/$f" > "$RUN/$f"
done
cp "$ARCH/drv_flds_in" "$RUN/"
ncgen -o "$RUN/fates_params_default.nc" "$CDL"

python3 - "$SITE" "$ARCH" "$RUN" "$D" "$IN" "$START_YMD" "$STEPS" "$OLD_D" "$OLD_IN" <<'PY'
import sys, re, os
SITE, ARCH, RUN, D, IN, start_ymd, steps, OLD_D, OLD_IN = sys.argv[1:10]
fix = lambda line: line.replace(OLD_D, D).replace(OLD_IN, IN)

repl = {'case_name': 'fates_fire_1pt', 'case_desc': 'FATES SPITFIRE parity',
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
        'fates_paramfile': f"'{RUN}/fates_params_default.nc'", 'use_lch4': '.false.',
        # scalar_lightning: constant ignitions (fates_fire_nignitions), no streams.
        'fates_spitfire_mode': os.environ['SPITFIRE_MODE']}
mode = int(os.environ['SPITFIRE_MODE'])
FD = f'{IN}/lnd/clm2/firedata'
MESH = f"'{IN}/share/meshes/fv0.9x1.25_141008_polemod_ESMFmesh.nc'"
if mode >= 2:   # FATESFireData reads BOTH streams for every mode > 1 (hdm_init)
    over['stream_fldfilename_popdens'] = f"'{FD}/clmforc.UNIFORM_hdm_0.9x1.25_CLMjl-parity.nc'"
    over['stream_meshfile_popdens'] = MESH
FIELDS = ['FATES_NESTEROV_INDEX', 'FATES_IGNITIONS', 'FATES_FDI', 'FATES_ROS',
          'FATES_EFFECT_WSPEED', 'FATES_FUELCONSUMED', 'FATES_FIRE_INTENSITY',
          'FATES_FIRE_INTENSITY_BURNFRAC', 'FATES_BURNFRAC', 'FATES_FUEL_AMOUNT',
          'FATES_FUEL_EFF_MOIST', 'FATES_FUEL_MEF', 'FATES_FUEL_BULKD', 'FATES_FUEL_SAV',
          'FATES_FIRE_CLOSS', 'FATES_DISTURBANCE_RATE_FIRE', 'FATES_MORTALITY_FIRE_SZ',
          'FATES_FUEL_MOISTURE_FC', 'FATES_FUEL_AMOUNT_FC', 'FATES_TVEG24', 'FATES_VEGC',
          'FATES_NPATCHES', 'FATES_GPP', 'FATES_LITTER_AG_CWD_EL', 'FATES_LITTER_AG_FINE_EL',
          'WIND', 'RH', 'RAIN', 'SNOW', 'TBOT', 'TOTSOMC', 'TOTLITC', 'HR']
seen, out = set(), []
for line in open(f'{ARCH}/lnd_in'):
    line = fix(line)
    m = re.match(r'^(\s*)(\w+)(\s*=\s*)(.*)$', line.rstrip('\n'))
    if m and m.group(2) in over:
        k = m.group(2); seen.add(k)
        out.append(f"{m.group(1)}{k}{m.group(3)}{over[k]}\n")
    else:
        out.append(line)
    if line.strip() == '&light_streams' and mode >= 2:
        out.append(f" lightngmapalgo = 'bilinear'\n stream_fldfilename_lightng = "
                   f"'{FD}/clmforc.UNIFORM_lnfm_0.9x1.25_CLMjl-parity.nc'\n"
                   f" stream_meshfile_lightng = {MESH}\n"
                   " stream_year_first_lightng = 2000\n stream_year_last_lightng = 2000\n")
    if line.strip() == '&clm_inparm':
        out.append(" hist_nhtfrq = -24\n hist_mfilt = 100000\n hist_fincl1 = " +
                   ",".join(f"'{f}'" for f in FIELDS) + "\n")
open(f'{RUN}/lnd_in', 'w').writelines(out)
assert not (set(over) - seen), set(over) - seen
print('[ok] namelists staged in', RUN)
PY

# Same launch environment and retry loop as setup_case.sh (macOS 26 xzone/ESMF trap).
cd "$RUN"
for attempt in $(seq 1 30); do
  rm -f fates_fire_1pt.clm2.h0.*.nc
  if MallocNanoZone=0 HWLOC_COMPONENTS=-opencl "$EXE" > cesm.stdout 2>&1 \
     && ls fates_fire_1pt.clm2.h0.*.nc > /dev/null 2>&1; then
    echo "[ok] run succeeded on attempt $attempt"; break
  fi
  echo "[..] attempt $attempt did not complete; retrying"
done
ls "$RUN"/fates_fire_1pt.clm2.h0.*.nc
