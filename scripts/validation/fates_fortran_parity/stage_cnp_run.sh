#!/usr/bin/env bash
# =============================================================================
# stage_cnp_run.sh — Fortran reference for FATES CNP (fates_parteh_mode = 2).
#
# Stages and runs the archived fates_parity_1pt case (Bow at Banff single point,
# compset I2000Clm50FatesRs => use_fates, use_fates_bgc, CENTURYKoven2013,
# suplnitro='ALL', use_nitrif_denitrif) with an UNINSTRUMENTED cesm.exe and an
# HOURLY h0 history of the soil-BGC fields. No rebuild is needed: FATES is always
# compiled into CTSM, and the soil-BGC fields are standard history output.
#
#   case config : $SYMFLUENCE_DATA/installs/clm/case-archive/fates_parity_1pt/CaseDocs/lnd_in
#   site        : $SYMFLUENCE_DATA/domain_Bow_at_Banff_lumped_era5/settings/CLM (read only)
#   run dir     : $RUN (default $SYMFLUENCE_DATA/runs/fates_cnp_1pt)
#
# Same overrides as stage_soilbgc_run.sh, plus the CTSM FatesColdPRT2 testmod CNP
# settings: fates_parteh_mode = 2 and fates_cnp_prescribed_nuptake/puptake = 1 for
# all PFTs (CTSM only runs CNP with prescribed uptake: FatesInterfaceMod check_allset
# endruns otherwise, and the ctsm5.3 host never fills bc_in%plant_*_uptake_flux).
# Start 2002-07-01, dtime 3600 s, STEPS steps (default 720 = 30 days).
#
# then score the port:
#   julia --project=. scripts/validation/fates_fortran_parity/fates_cnp_parity.jl \
#         $RUN/fates_cnp_1pt.clm2.h0.2002-07-01-00000.nc 720
# =============================================================================
set -euo pipefail
D=${SYMFLUENCE_DATA:?set SYMFLUENCE_DATA}
IN=${CESM_INPUTDATA:-$D/installs/cesm-inputdata}
EXE=${EXE:-$D/installs/clm/bin/cesm.exe}
SITE=${SITE:-$D/domain_Bow_at_Banff_lumped_era5/settings/CLM}
ARCH=${ARCH:-$D/installs/clm/case-archive/fates_parity_1pt/CaseDocs}
RUN=${RUN:-$D/runs/fates_cnp_1pt}
STEPS=${STEPS:-720}
# PARTEH=1 stages the same case carbon-only (control run; N/P fields dropped).
PARTEH=${PARTEH:-2}
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
# FatesColdPRT2: prescribed N and P uptake = 1.0 for every PFT.
sed -E -e 's/^( fates_cnp_prescribed_nuptake = ).*;/\11, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ;/' \
       -e 's/^( fates_cnp_prescribed_puptake = ).*;/\11, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ;/' \
       "$CDL" > "$RUN/fates_params_prt2.cdl"
grep -q '^ fates_cnp_prescribed_nuptake = 1, 1' "$RUN/fates_params_prt2.cdl"
grep -q '^ fates_cnp_prescribed_puptake = 1, 1' "$RUN/fates_params_prt2.cdl"
ncgen -o "$RUN/fates_params_prt2.nc" "$RUN/fates_params_prt2.cdl"

python3 - "$SITE" "$ARCH" "$RUN" "$D" "$IN" "$START_YMD" "$STEPS" "$OLD_D" "$OLD_IN" "$PARTEH" <<'PY'
import sys, re
SITE, ARCH, RUN, D, IN, start_ymd, steps, OLD_D, OLD_IN, parteh = sys.argv[1:11]
fix = lambda line: line.replace(OLD_D, D).replace(OLD_IN, IN)

repl = {'case_name': 'fates_cnp_1pt', 'case_desc': 'FATES CNP parity',
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
        'fates_paramfile': f"'{RUN}/fates_params_prt2.nc'", 'use_lch4': '.false.',
        'fates_parteh_mode': parteh}
FIELDS = ['FATES_VEGN', 'FATES_STOREN', 'FATES_LEAFN', 'FATES_FROOTN', 'FATES_SAPWOODN',
          'FATES_REPRON', 'FATES_STOREN_TF', 'FATES_VEGP', 'FATES_STOREP', 'FATES_LEAFP',
          'FATES_FROOTP', 'FATES_SAPWOODP', 'FATES_REPROP', 'FATES_STOREP_TF',
          'FATES_NH4UPTAKE', 'FATES_NO3UPTAKE', 'FATES_NDEMAND', 'FATES_NEFFLUX',
          'FATES_NFIX_SYM', 'FATES_PUPTAKE', 'FATES_PDEMAND', 'FATES_PEFFLUX',
          'FATES_EXCESS_RESP', 'FATES_VEGC', 'FATES_LEAFC', 'FATES_STOREC', 'FATES_FROOTC',
          'FATES_SAPWOODC', 'FATES_REPROC', 'FATES_NPP', 'FATES_GPP', 'FATES_NCOHORTS',
          'FATES_L2FR', 'FATES_LAI', 'SMINN', 'SMIN_NH4', 'SMIN_NO3', 'HR', 'TOTLITC',
          'TOTSOMC', 'FATES_LITTER_IN', 'FATES_LITTER_OUT', 'FATES_VEGN_SZPF',
          'FATES_STOREN_SZPF', 'FATES_VEGP_SZPF', 'FATES_STOREN_TF_CANOPY_SZPF',
          'FATES_NDEMAND_SZPF', 'FATES_LEAFC_PF', 'FATES_NPLANT_PF', 'FATES_MORTALITY_PF',
          'FATES_STOREC_PF', 'FATES_VEGC_PF', 'FATES_STOREC_TF', 'FATES_LEAFN_SZPF', 'FATES_COLD_STATUS',
          'FATES_NCOLDDAYS', 'FATES_GDD', 'FATES_DAYSINCE_COLDLEAFON', 'FATES_DAYSINCE_COLDLEAFOFF',
          'FATES_NCHILLDAYS', 'FATES_TVEG24']
if parteh == '1':   # the N/P FATES fields are only registered in CNP mode
    NP_ELEM = re.compile(r'^FATES_(VEG|STORE|LEAF|FROOT|SAPWOOD|REPRO)[NP](_|$)|'
                         r'^FATES_(NH4UPTAKE|NO3UPTAKE|NDEMAND|NEFFLUX|NFIX_SYM|'
                         r'PUPTAKE|PDEMAND|PEFFLUX)')
    FIELDS = [f for f in FIELDS if not NP_ELEM.match(f)]
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
  rm -f fates_cnp_1pt.clm2.h0.*.nc
  if MallocNanoZone=0 HWLOC_COMPONENTS=-opencl "$EXE" > cesm.stdout 2>&1 \
     && ls fates_cnp_1pt.clm2.h0.*.nc > /dev/null 2>&1; then
    echo "[ok] run succeeded on attempt $attempt"; break
  fi
  echo "[..] attempt $attempt did not complete; retrying"
done
ls "$RUN"/fates_cnp_1pt.clm2.h0.*.nc
