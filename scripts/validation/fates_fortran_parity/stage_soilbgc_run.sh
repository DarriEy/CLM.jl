#!/usr/bin/env bash
# =============================================================================
# stage_soilbgc_run.sh — Fortran reference for the FATES -> CLM soil-BGC coupling.
#
# Stages and runs the archived fates_parity_1pt case (Bow at Banff single point,
# compset I2000Clm50FatesRs => use_fates, use_fates_bgc, CENTURYKoven2013,
# suplnitro='ALL', use_nitrif_denitrif) with an UNINSTRUMENTED cesm.exe and an
# HOURLY h0 history of the soil-BGC fields. No rebuild is needed: FATES is always
# compiled into CTSM, and the soil-BGC fields are standard history output.
#
#   case config : $SYMFLUENCE_DATA/installs/clm/case-archive/fates_parity_1pt/CaseDocs/lnd_in
#   site        : $SYMFLUENCE_DATA/domain_Bow_at_Banff_lumped_era5/settings/CLM (read only)
#   run dir     : $RUN (default $SYMFLUENCE_DATA/runs/fates_bgc_1pt)
#
# Same overrides as setup_case.sh (site fsurdat/paramfile, the FATES default param
# file ncgen'd from the CTSM source .cdl — byte-identical to data/fates/ — and
# use_lch4 off), start 2002-07-01, dtime 3600 s, STEPS steps (default 240 = 10 days).
#
# then score the port:
#   julia --project=. scripts/validation/fates_fortran_parity/fates_soilbgc_parity.jl \
#         $RUN/fates_bgc_1pt.clm2.h0.2002-07-01-00000.nc 240
# =============================================================================
set -euo pipefail
D=${SYMFLUENCE_DATA:?set SYMFLUENCE_DATA}
IN=${CESM_INPUTDATA:-$D/installs/cesm-inputdata}
EXE=${EXE:-$D/installs/clm/bin/cesm.exe}
SITE=${SITE:-$D/domain_Bow_at_Banff_lumped_era5/settings/CLM}
ARCH=${ARCH:-$D/installs/clm/case-archive/fates_parity_1pt/CaseDocs}
RUN=${RUN:-$D/runs/fates_bgc_1pt}
STEPS=${STEPS:-240}
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
import sys, re
SITE, ARCH, RUN, D, IN, start_ymd, steps, OLD_D, OLD_IN = sys.argv[1:10]
fix = lambda line: line.replace(OLD_D, D).replace(OLD_IN, IN)

repl = {'case_name': 'fates_bgc_1pt', 'case_desc': 'FATES soil BGC parity',
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
        'fates_paramfile': f"'{RUN}/fates_params_default.nc'", 'use_lch4': '.false.'}
FIELDS = ['HR', 'LITTERC_HR', 'SOILC_HR', 'TOTLITC', 'TOTSOMC', 'LIT_MET_C', 'LIT_CEL_C',
          'LIT_LIG_C', 'SOM_ACT_C', 'SOM_SLO_C', 'SOM_PAS_C', 'HR_vr', 'LIT_MET_C_vr',
          'LIT_CEL_C_vr', 'LIT_LIG_C_vr', 'SOM_ACT_C_vr', 'T_SCALAR', 'W_SCALAR',
          'O_SCALAR', 'FATES_LITTER_OUT', 'FATES_LITTER_IN', 'FATES_HET_RESP', 'FATES_NEP',
          'FPI', 'TSOI', 'H2OSOI', 'SOILC_vr', 'FATES_GPP', 'FATES_LITTER_AG_FINE_EL',
          'FATES_LITTER_BG_FINE_EL']
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
open(f'{RUN}/lnd_in', 'w').writelines(out)
assert not (set(over) - seen), set(over) - seen
print('[ok] namelists staged in', RUN)
PY

# Same launch environment and retry loop as setup_case.sh (macOS 26 xzone/ESMF trap).
cd "$RUN"
for attempt in $(seq 1 30); do
  rm -f fates_bgc_1pt.clm2.h0.*.nc
  if MallocNanoZone=0 HWLOC_COMPONENTS=-opencl "$EXE" > cesm.stdout 2>&1 \
     && ls fates_bgc_1pt.clm2.h0.*.nc > /dev/null 2>&1; then
    echo "[ok] run succeeded on attempt $attempt"; break
  fi
  echo "[..] attempt $attempt did not complete; retrying"
done
ls "$RUN"/fates_bgc_1pt.clm2.h0.*.nc
