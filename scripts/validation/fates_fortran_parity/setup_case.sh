#!/usr/bin/env bash
# =============================================================================
# setup_case.sh — build an instrumented, FATES-enabled Fortran CTSM executable
# and run it SINGLE-POINT to emit a cold-start (nstep=0) schema-v1 pdump: the
# bit-parity ground truth for the Julia CLM.jl FATES port.
#
# THIS SCRIPT RUNS END-TO-END AND PRODUCES THE DUMP. (The previous revision stood
# up a GLOBAL ne3 GSWP3 case that aborted in DATM on 743 absent forcing files;
# that blocker was self-inflicted — see "WHY SINGLE-POINT" below.)
#
# ── WHY SINGLE-POINT ─────────────────────────────────────────────────────────
# FATES does not need a global grid. Driving the case from a global GSWP3 compset
# pulls in hundreds of GB of DATM stream files that are not staged here, and DATM
# aborts in shr_stream_init_from_xml -> pio_openfile during ATM InitializeRealize,
# BEFORE the LND component is ever realized — so FATES cold start is never reached.
#
# This project already has a PROVEN single-point Fortran CLM path (the one behind
# every pdump parity comparison to date): CIME is used ONLY to build cesm.exe, and
# the model is then run directly in a hand-staged run directory whose namelists
# point at local site forcing + a 1x1 ESMF mesh. We reuse exactly that path and
# just add FATES. Result: zero missing input files.
#
#   build : compset I2000Clm50FatesRs = 2000_DATM%GSWP3v1_CLM50%FATES_SICE_SOCN_
#           SROF_SGLC_SWAV — the FATES analogue of the working I2000Clm50SpRs
#           build (SROF, i.e. no MOSART: components are CPL/ATM/LND only).
#           --res is build-only; the runtime grid comes from the namelists.
#   run   : namelists copied from the Bow single-point settings dir, with a
#           FATES-enabled lnd_in. DATM reads local clmforc.YYYY.nc; nx=ny=1.
#
# ── OUTPUT ───────────────────────────────────────────────────────────────────
#   $RUN/fates_pdump_fortran_n0.txt   (1 site x 1 patch x 14 cohorts, one per PFT)
# then score the Julia port against it:
#   julia +1.12 --project=scripts scripts/fates_fortran_parity.jl compare <dump>
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
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# CIME needs a python with a working distutils AND a stdlib-compatible ElementTree.
# (System python3.14 has no distutils; miniforge base's lxml shadows ElementTree.)
CIME_PY=${CIME_PY:-/Users/darri.eythorsson/.local/share/mamba/envs/ewc-summa/bin}
export PATH="$CIME_PY:$PATH"

# ---- 1. FATES parameter file: ncgen from the shared, byte-identical .cdl -----
# data/fates/fates_params_default.cdl (port) == src/fates/parameter_files/...cdl
# (install) — verified byte-identical, so cold-start cohorts are parameter-identical
# by construction and need NO tuning to compare. 14 PFTs -> 14 cold-start cohorts.
CDL="$CLMROOT/src/fates/parameter_files/fates_params_default.cdl"
FATES_PARAM="$INPUTDATA/lnd/clm2/paramdata/fates_params_default.nc"
mkdir -p "$(dirname "$FATES_PARAM")"
ncgen -o "$FATES_PARAM" "$CDL"
echo "[ok] FATES params -> $FATES_PARAM"

# ---- 2. create + configure the case (build only) ----------------------------
rm -rf "$CASE"
"$CLMROOT/cime/scripts/create_newcase" --case "$CASE" --compset "$COMPSET" \
    --res "$RES" --run-unsupported --machine homebrew --compiler gnu
cd "$CASE"
./xmlchange STOP_N=1,STOP_OPTION=nsteps,RUN_STARTDATE=2000-01-01
./xmlchange NTASKS=1,NTASKS_ESP=1     # ESMF is built mpiuni => one PET
./case.setup

# ---- 3. install the FATES cold-start pdump instrumentation ------------------
# FATES compiles on the CLM Filepath, so the FATES SourceMod lives in src.clm.
SM="$CASE/SourceMods/src.clm"
cp "$HERE/FatesParityDumpMod.F90" "$SM/FatesParityDumpMod.F90"
python3 - "$CLMROOT/src/utils/clmfates_interfaceMod.F90" "$SM/clmfates_interfaceMod.F90" <<'PY'
import sys
src, dst = sys.argv[1], sys.argv[2]
out, used, called = [], False, False
for line in open(src).read().splitlines(keepends=True):
    if (not used) and 'use EDInitMod' in line and 'set_site_properties' in line:
        out.append(line)
        out.append("   use FatesParityDumpMod    , only : fates_parity_dump   ! FATES parity instrumentation\n")
        used = True; continue
    if (not called) and "t_stopf('fates_initcoldstart')" in line:
        out += [
          "     ! ---- FATES cold-start parity instrumentation (schema-v1 pdump) ------------\n",
          "     do nc = 1, nclumps\n",
          "        if ( this%fates(nc)%nsites>0 ) then\n",
          "           call fates_parity_dump(this%fates(nc)%sites, this%fates(nc)%nsites, &\n",
          "                get_nstep(), 'fates_pdump_fortran_n0.txt')\n",
          "        end if\n",
          "     end do\n",
          "     ! --------------------------------------------------------------------------\n",
        ]
        called = True
    out.append(line)
open(dst, 'w').writelines(out)
assert used and called, (used, called)
print("[ok] instrumented clmfates_interfaceMod.F90 (use + call injected)")
PY

# ---- 4. build ---------------------------------------------------------------
# The homebrew mach config hardcodes a netcdf Cellar version; if that version is
# gone the final link fails with "ld: library 'netcdf' not found". Point the stale
# versioned dir at the installed one (link-time only; the runtime dylib resolves
# via its install-name regardless).
NCDIR=/opt/homebrew/Cellar/netcdf
[ -e "$NCDIR/4.9.3_2" ] || ln -sfn "$(brew --prefix netcdf)" "$NCDIR/4.9.3_2"
./case.build
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

python3 - "$SITE/nuopc.runconfig" "$RUN/nuopc.runconfig" <<'PY'
import sys, re
src, dst = sys.argv[1], sys.argv[2]
# 1 step suffices: the dump fires in init_coldstart, before any timestep runs.
repl = {'case_name': 'fates_parity_1pt', 'case_desc': 'FATES cold-start parity',
        'stop_option': 'nsteps', 'stop_n': '1', 'restart_option': 'never'}
out = []
for line in open(src):
    m = re.match(r'^(\s*)(\w+)( = )(.*)$', line.rstrip('\n'))
    out.append(f"{m.group(1)}{m.group(2)}{m.group(3)}{repl[m.group(2)]}\n"
               if m and m.group(2) in repl else line)
open(dst, 'w').writelines(out)
print("[ok] single-point nuopc.runconfig ->", dst)
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
    # FATES cold start does not need methane. (The popdens/lightning fire streams
    # and the crop-calendar mesh are also absent but are never READ: fates_spitfire_mode=0
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
# RETRY LOOP: the optimized build still SIGTRAPs inside libsystem_malloc on roughly
#   half of all launches. The trap lands in a DIFFERENT allocation each time (netcdf
#   NC3_open, PatchType::Init, ...), which is the ESMF/xzone false-positive, not a
#   CLM bug — every successful run emits a BYTE-IDENTICAL dump (verified 6/6), and a
#   DEBUG (-O0 + bounds-checked) build runs clean and produces the same dump.
cd "$RUN"
for attempt in $(seq 1 20); do
    rm -f fates_pdump_fortran_n0.txt
    if MallocNanoZone=0 HWLOC_COMPONENTS=-opencl "$EXE" > cesm.stdout 2>&1 \
       && [ -f fates_pdump_fortran_n0.txt ]; then
        echo "[ok] run succeeded on attempt $attempt"
        break
    fi
    echo "[..] attempt $attempt hit the xzone malloc trap; retrying"
done
[ -f fates_pdump_fortran_n0.txt ] || { echo "[FAIL] no dump after 20 attempts"; exit 1; }

echo
echo "=== Fortran FATES cold-start ground truth ==================================="
head -4 "$RUN/fates_pdump_fortran_n0.txt"
echo "  ..."
echo
echo "dump: $RUN/fates_pdump_fortran_n0.txt"
echo "score the Julia port against it:"
echo "  julia +1.12 --project=scripts scripts/fates_fortran_parity.jl compare \\"
echo "      $RUN/fates_pdump_fortran_n0.txt"
