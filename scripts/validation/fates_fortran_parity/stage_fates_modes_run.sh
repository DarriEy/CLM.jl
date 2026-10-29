#!/usr/bin/env bash
# =============================================================================
# stage_fates_modes_run.sh — Fortran references for the reduced-complexity FATES
# modes (use_fates_sp / use_fates_nocomp / use_fates_fixed_biogeog).
#
# Same single-point case as stage_soilbgc_run.sh (Bow at Banff, archived
# fates_parity_1pt namelists, uninstrumented cesm.exe, 2002-07-01, dtime 3600 s,
# hourly h0), plus the namelist overrides of one MODE:
#
#   MODE=sp        use_fates_sp + use_fates_nocomp + use_fates_fixed_biogeog
#                  (CTSM requires all three; use_fates_bgc is then .false.)
#   MODE=nocomp    use_fates_nocomp + use_fates_fixed_biogeog
#   MODE=fixedbg   use_fates_fixed_biogeog only
#
# The fixed-biogeography PFT areas come from the site surfdata PCT_NAT_PFT; SP
# reads its prescribed LAI/SAI/height from the same file's MONTHLY_* fields.
#
#   run dir : $RUN (default $SYMFLUENCE_DATA/runs/fates_<MODE>_1pt)
#
# then score the port:
#   julia --project=. scripts/validation/fates_fortran_parity/fates_modes_parity.jl \
#         <MODE> $RUN/fates_<MODE>_1pt.clm2.h0.2002-07-01-00000.nc 240
# =============================================================================
set -euo pipefail
MODE=${MODE:?set MODE=sp|nocomp|fixedbg}
D=${SYMFLUENCE_DATA:?set SYMFLUENCE_DATA}
IN=${CESM_INPUTDATA:-$D/installs/cesm-inputdata}
EXE=${EXE:-$D/installs/clm/bin/cesm.exe}
SITE=${SITE:-$D/domain_Bow_at_Banff_lumped_era5/settings/CLM}
ARCH=${ARCH:-$D/installs/clm/case-archive/fates_parity_1pt/CaseDocs}
CASE_NAME=fates_${MODE}_1pt
RUN=${RUN:-$D/runs/$CASE_NAME}
STEPS=${STEPS:-240}
START_YMD=${START_YMD:-20020701}
CDL=${CDL:-$D/installs/clm/ctsm5.3.012/src/fates/parameter_files/fates_params_default.cdl}
OLD_D=/Users/darri.eythorsson/compHydro/SYMFLUENCE_data
OLD_IN=/Users/darri.eythorsson/projects/cesm-inputdata

case "$MODE" in
  # SP turns CTSM's soil BGC off (use_fates_bgc = .not. use_fates_sp), and
  # controlMod then requires soil_decomp_method = 'None'.
  sp)      FLAGS="use_fates_sp=.true. use_fates_nocomp=.true. use_fates_fixed_biogeog=.true. soil_decomp_method='None'" ;;
  nocomp)  FLAGS="use_fates_sp=.false. use_fates_nocomp=.true. use_fates_fixed_biogeog=.true." ;;
  fixedbg) FLAGS="use_fates_sp=.false. use_fates_nocomp=.false. use_fates_fixed_biogeog=.true." ;;
  *) echo "unknown MODE=$MODE"; exit 1 ;;
esac

rm -rf "$RUN"; mkdir -p "$RUN/timing/checkpoints"
for f in nuopc.runseq fd.yaml datm_in datm.streams.xml drv_in; do
  sed -e "s#$OLD_D#$D#g" -e "s#$OLD_IN#$IN#g" "$SITE/$f" > "$RUN/$f"
done
cp "$ARCH/drv_flds_in" "$RUN/"
ncgen -o "$RUN/fates_params_default.nc" "$CDL"

python3 - "$SITE" "$ARCH" "$RUN" "$D" "$IN" "$START_YMD" "$STEPS" "$OLD_D" "$OLD_IN" \
          "$CASE_NAME" $FLAGS <<'PY'
import sys, re
SITE, ARCH, RUN, D, IN, start_ymd, steps, OLD_D, OLD_IN, case = sys.argv[1:11]
flags = dict(a.split('=', 1) for a in sys.argv[11:])
fix = lambda line: line.replace(OLD_D, D).replace(OLD_IN, IN)

repl = {'case_name': case, 'case_desc': 'FATES reduced-complexity mode parity',
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
        **flags}
# (GPP is not an HLM field under FATES; FATES_GPP is. The FATES_NOCOMP_*_PF fields
# are only registered when nocomp is on.)
FIELDS = ['TLAI', 'TSAI', 'ELAI', 'ESAI', 'HTOP', 'FSH', 'EFLX_LH_TOT', 'FSA', 'TV',
          'FATES_GPP', 'FATES_LAI', 'FATES_NCOHORTS', 'FATES_NPATCHES', 'FATES_VEGC',
          'FATES_LEAFC', 'FATES_GPP_PF', 'FATES_VEGC_PF', 'FATES_LEAFC_PF',
          'FATES_NPLANT_PF', 'FATES_CANOPYCROWNAREA_PF', 'FATES_CROWNAREA_PF']
if flags.get('use_fates_nocomp') == '.true.':
    FIELDS += ['FATES_NOCOMP_PATCHAREA_PF', 'FATES_NOCOMP_NPATCHES_PF']
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
  # A namelist/model ENDRUN is deterministic -- only the malloc trap is worth a retry.
  if grep -q "ENDRUN" lnd.log 2>/dev/null; then grep -B8 -m1 ENDRUN lnd.log; exit 1; fi
  echo "[..] attempt $attempt did not complete; retrying"
done
ls "$RUN"/"$CASE_NAME".clm2.h0.*.nc
