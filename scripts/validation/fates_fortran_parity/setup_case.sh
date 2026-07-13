#!/usr/bin/env bash
# =============================================================================
# setup_case.sh — build an instrumented, FATES-enabled Fortran CTSM executable and
# run it SINGLE-POINT to emit a PHASE-TAGGED, TIME-STEPPED schema-v2 pdump: the
# ground truth for the Julia CLM.jl FATES port's DYNAMICS (not just its cold start).
#
# ── WHAT CHANGED vs the cold-start-only revision (PR #213) ───────────────────
# PR #213 proved FATES *initialization*: 27/27 fields at nstep=0. But the FATES
# parameter file is byte-identical between port and install, so those cohorts are
# parameter-identical BY CONSTRUCTION — it tested nothing about photosynthesis,
# allocation, growth, mortality, recruitment or disturbance.
#
# This revision runs the SAME single-point case FORWARD and dumps at every FATES
# phase boundary (see FatesParityDumpMod.F90 for the phase table):
#     phase  1  fast        every timestep, after the photosynthesis/flux solve
#     phase 10..18          the daily step, sub-phase by sub-phase
#                           (dyn_in / phenology / distrates / integrate / recruit /
#                            cohortfuse / spawn / patchfuse / updatesite)
# so a Julia-vs-Fortran divergence can be ATTRIBUTED to the phase that produced it.
#
# ── WHY SINGLE-POINT ─────────────────────────────────────────────────────────
# FATES does not need a global grid. A global GSWP3 compset pulls in hundreds of GB
# of DATM stream files that are not staged here, and DATM aborts in
# shr_stream_init_from_xml -> pio_openfile during ATM InitializeRealize, BEFORE the
# LND component is realized — so FATES is never reached. This project already has a
# PROVEN single-point Fortran CLM path: CIME is used ONLY to build cesm.exe, and the
# model is run directly in a hand-staged run dir whose namelists point at local site
# forcing + a 1x1 ESMF mesh. We reuse exactly that path and add FATES.
#
#   build : compset I2000Clm50FatesRs = 2000_DATM%GSWP3v1_CLM50%FATES_SICE_SOCN_
#           SROF_SGLC_SWAV  (FATES + SROF: components are CPL/ATM/LND only).
#           --res is build-only; the runtime grid comes from the namelists.
#   run   : namelists copied from the Bow single-point settings dir, with a
#           FATES-enabled lnd_in. DATM reads local clmforc.YYYY.nc; nx=ny=1.
#           dtime = lnd_cpl_dt = 3600 s -> 24 steps/day.
#
# ── RUN WINDOW ───────────────────────────────────────────────────────────────
# Default start 2002-07-01 (NOT the January of the cold-start case): Bow at Banff is
# a 51N montane site, so a January window leaves every cohort dormant with GPP == 0 —
# which would "validate" the fast phase against a column of zeros. Starting in
# mid-summer exercises real photosynthesis, real daily allocation and real growth.
# STEPS=120 -> 5 days -> 5 firings of the daily dynamics (nstep = 24,48,72,96,120).
#
# ── OUTPUT ───────────────────────────────────────────────────────────────────
#   $RUN/fates_pdump_fortran.txt          (all phases, all steps, one file)
# score the Julia port against it:
#   julia +1.12 --project=. scripts/fates_fortran_parity.jl compare <dump>
# =============================================================================
set -euo pipefail

# ---- paths (edit for your environment) --------------------------------------
CLMROOT=${CLMROOT:-/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/installs/clm}
INPUTDATA=${INPUTDATA:-/Users/darri.eythorsson/projects/cesm-inputdata}
CASE=${CASE:-$CLMROOT/cases/fates_parity_1pt}
RUN=${RUN:-/Users/darri.eythorsson/projects/scratch/fates_parity_1pt/run1pt}
# Proven single-point settings dir (namelists + 1x1 mesh + surfdata + params):
SITE=${SITE:-/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped_era5/settings/CLM}
COMPSET=${COMPSET:-I2000Clm50FatesRs}   # FATES + SROF (no MOSART)
RES=${RES:-f09_g17}                     # build-only; runtime grid is the 1x1 mesh
START_YMD=${START_YMD:-20020701}        # summer -> the fast phase is not all zeros
STEPS=${STEPS:-120}                     # 3600 s steps -> 120 = 5 days
SKIP_BUILD=${SKIP_BUILD:-0}             # 1 -> reuse the existing cesm.exe
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# CIME needs a python with a working distutils AND a stdlib-compatible ElementTree.
# (System python3.14 has no distutils; miniforge base's lxml shadows ElementTree.)
CIME_PY=${CIME_PY:-/Users/darri.eythorsson/.local/share/mamba/envs/ewc-summa/bin}
export PATH="$CIME_PY:$PATH"

# ---- 1. FATES parameter file: ncgen from the shared, byte-identical .cdl -----
# data/fates/fates_params_default.cdl (port) == src/fates/parameter_files/...cdl
# (install) — verified byte-identical, so cold-start cohorts are parameter-identical
# by construction. 14 PFTs -> 14 cold-start cohorts.
CDL="$CLMROOT/src/fates/parameter_files/fates_params_default.cdl"
FATES_PARAM="$INPUTDATA/lnd/clm2/paramdata/fates_params_default.nc"
mkdir -p "$(dirname "$FATES_PARAM")"
ncgen -o "$FATES_PARAM" "$CDL"
echo "[ok] FATES params -> $FATES_PARAM"

# ---- 2. create + configure the case (build only) ----------------------------
if [ "$SKIP_BUILD" = "0" ]; then
  rm -rf "$CASE"
  # `u` answers create_newcase's "EXEROOT/RUNDIR already exists" prompts with
  # "use existing" -> an INCREMENTAL rebuild of the object tree (only the changed
  # SourceMods and their dependents recompile). Non-interactive: without the piped
  # answers create_newcase dies with EOFError on a second run.
  printf 'u\nu\n' | "$CLMROOT/cime/scripts/create_newcase" --case "$CASE" \
      --compset "$COMPSET" --res "$RES" --run-unsupported \
      --machine homebrew --compiler gnu
  cd "$CASE"
  ./xmlchange STOP_N=1,STOP_OPTION=nsteps,RUN_STARTDATE=2000-01-01
  ./xmlchange NTASKS=1,NTASKS_ESP=1     # ESMF is built mpiuni => one PET
  ./case.setup

  # ---- 3. install the phase-tagged pdump instrumentation ---------------------
  # FATES compiles on the CLM Filepath, so SourceMods/src.clm overrides FATES files
  # too (EDMainMod.F90) as well as CLM's own (clm_driver.F90, clmfates_interfaceMod).
  python3 "$HERE/instrument.py" "$CLMROOT" "$CASE/SourceMods/src.clm"

  # ---- 4. build --------------------------------------------------------------
  # The homebrew mach config hardcodes a netcdf Cellar version; if that version is
  # gone the final link fails with "ld: library 'netcdf' not found". Point the stale
  # versioned dir at the installed one (link-time only).
  NCDIR=/opt/homebrew/Cellar/netcdf
  [ -e "$NCDIR/4.9.3_2" ] || ln -sfn "$(brew --prefix netcdf)" "$NCDIR/4.9.3_2"
  ./case.build
else
  cd "$CASE"
fi
EXE="$(./xmlquery -value EXEROOT)/cesm.exe"
nm "$EXE" | grep -q fatesparitydumpmod_MOD_fates_parity_dump \
    && echo "[ok] instrumentation linked into $EXE"
CIME_RUNDIR="$(./xmlquery -value RUNDIR)"

# ---- 5. stage the SINGLE-POINT run directory --------------------------------
# Driver/DATM config: the proven Bow single-point namelists (1x1 ESMF mesh,
# datamode=CLMNCEP, local clmforc.YYYY.nc). LND config: this case's CIME-generated
# FATES lnd_in with the data-file paths repointed at the site files.
rm -rf "$RUN"; mkdir -p "$RUN/timing/checkpoints"
cp "$SITE"/{nuopc.runseq,fd.yaml,datm_in,datm.streams.xml,drv_in} "$RUN/"
cp "$CIME_RUNDIR/drv_flds_in" "$RUN/"
echo "$CASE" > "$RUN/CASEROOT"

python3 - "$SITE/nuopc.runconfig" "$RUN/nuopc.runconfig" "$START_YMD" "$STEPS" <<'PY'
import sys, re
src, dst, start_ymd, steps = sys.argv[1:5]
repl = {'case_name': 'fates_parity_1pt', 'case_desc': 'FATES time-stepped parity',
        'start_ymd': start_ymd, 'stop_option': 'nsteps', 'stop_n': steps,
        'restart_option': 'never'}
out = []
for line in open(src):
    m = re.match(r'^(\s*)(\w+)( = )(.*)$', line.rstrip('\n'))
    out.append(f"{m.group(1)}{m.group(2)}{m.group(3)}{repl[m.group(2)]}\n"
               if m and m.group(2) in repl else line)
open(dst, 'w').writelines(out)
print(f"[ok] single-point nuopc.runconfig (start={start_ymd}, {steps} steps) -> {dst}")
PY

python3 - "$CIME_RUNDIR/lnd_in" "$RUN/lnd_in" "$SITE/parameters" "$FATES_PARAM" <<'PY'
import sys, re
gen, dst, P, fates_param = sys.argv[1:5]
# Repoint the CIME-generated FATES namelist at the single-point data files.
# Everything else is left exactly as `CLM build-namelist -bgc fates` produced it.
over = {
    'fsurdat':         f"'{P}/surfdata_clm.nc'",
    'paramfile':       f"'{P}/clm5_params.nc'",
    'fates_paramfile': f"'{fates_param}'",
    # use_lch4 pulls in finundated_inversiondata streams that are not staged, and
    # FATES does not need methane here. (The popdens/lightning fire streams and the
    # crop-calendar mesh are also absent but are never READ: fates_spitfire_mode=0
    # selects FATESFireNoData, whose need_lightning_and_popdens() is .false., and
    # cropcals_rx / use_crop are .false.)
    'use_lch4':        '.false.',
}
seen, out = set(), []
for line in open(gen):
    m = re.match(r'^(\s*)(\w+)(\s*=\s*)(.*)$', line.rstrip('\n'))
    if m and m.group(2) in over:
        k = m.group(2); seen.add(k)
        out.append(f"{m.group(1)}{k}{m.group(3)}{over[k]}\n")
    else:
        out.append(line)
open(dst, 'w').writelines(out)
assert not (set(over) - seen), set(over) - seen
print("[ok] FATES single-point lnd_in ->", dst)
PY

# ---- 6. run -----------------------------------------------------------------
# HWLOC_COMPONENTS=-opencl : REQUIRED. hwloc's OpenCL probe walks into the Metal
#   driver at MPI startup and dies with SIGILL inside AGXMetalG16X on macOS 26;
#   without this the exe never reaches CLM at all.
# MallocNanoZone=0 : ESMF otherwise trips false heap-corruption reports.
# Run cesm.exe DIRECTLY, never under lldb: macOS 26's xzone allocator uses inline
#   `brk #0x1` guards that halt a debugger spuriously.
# RETRY LOOP: the optimized build SIGTRAPs inside libsystem_malloc on roughly half of
#   all launches. The trap lands in a DIFFERENT allocation each time (netcdf NC3_open,
#   PatchType::Init, ...), which is the known ESMF/xzone false positive, not a CLM bug:
#   every successful run emits a BYTE-IDENTICAL dump.
cd "$RUN"
DUMP=fates_pdump_fortran.txt
for attempt in $(seq 1 30); do
    rm -f "$DUMP"
    if MallocNanoZone=0 HWLOC_COMPONENTS=-opencl "$EXE" > cesm.stdout 2>&1 \
       && [ -f "$DUMP" ]; then
        echo "[ok] run succeeded on attempt $attempt"
        break
    fi
    echo "[..] attempt $attempt did not complete (xzone malloc trap?); retrying"
done
[ -f "$DUMP" ] || { echo "[FAIL] no dump after 30 attempts"; tail -30 cesm.stdout; exit 1; }

echo
echo "=== Fortran FATES time-stepped ground truth ================================="
python3 - "$RUN/$DUMP" <<'PY'
import sys, collections
c = collections.Counter()
for line in open(sys.argv[1]):
    if line.startswith('#'):
        continue
    f = dict(t.split('=', 1) for t in line.split()[1:] if '=' in t)
    c[(int(f.get('t', -1)), int(f.get('ph', -1)))] += 1
steps = sorted({t for t, _ in c})
phases = sorted({p for _, p in c})
print(f"  records : {sum(c.values())}")
print(f"  steps   : {len(steps)}  (nstep {steps[0]}..{steps[-1]})")
print(f"  phases  : {phases}")
PY
echo
echo "dump: $RUN/$DUMP"
echo "score the Julia port against it:"
echo "  julia +1.12 --project=. scripts/fates_fortran_parity.jl compare $RUN/$DUMP"
