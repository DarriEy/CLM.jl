#!/usr/bin/env bash
# =============================================================================
# gen_fortran_refs.sh — generate Fortran CTSM reference dumps for T1 parity.
#
# T1 is the only validation tier that needs external ground truth: a per-step
# snapshot of Fortran CLM state for a given config, which the harness diffs the
# Julia port against (oracle :parity, via compare_inst_to_dump). This script is
# the REPRODUCIBLE generator for those dumps — run it on the box that has the
# instrumented CTSM build, NOT in CI (the Julia harness only consumes the output).
#
# It is deliberately out-of-band: each config = a Fortran build+run+dump,
# minutes-to-hours of serial work (DESIGN §4). Spend it only where T2/T3 can't
# substitute — i.e. new physics paths whose absolute correctness matters and
# isn't implied by an invariant (crop, lch4/methane, a non-soil landunit).
#
# Prereq: an instrumented cesm.exe whose SourceMods restFileMod emits per-boundary
# `pdump_<boundary>_n<nstep>.nc` snapshots (the existing Bow dumps in $DUMPDIR were
# made this way). The instrumentation + the validated short-run recipe are banked
# in the project memory `fortran-clm-short-run-recipe`.
#
# Usage:
#   scripts/validation/gen_fortran_refs.sh <config> <nstep_start> <nstep_count>
#   e.g.  scripts/validation/gen_fortran_refs.sh lch4 13458 13   # methane path, 13 steps
# =============================================================================
set -euo pipefail

CONFIG="${1:?config name, e.g. crop|lch4|lake|glacier}"
NSTEP_START="${2:?first nstep to dump}"
NSTEP_COUNT="${3:-1}"

# --- paths (the instrumented Fortran build + the case) ---------------------------
# Discovered live on the build box (2026-06-25): a complete, BUILT, instrumented case.
#   exe:        $CASE_DIR/bld/cesm.exe  (arm64, pdump SourceMods baked in)
#   SourceMods: $CASE_DIR/SourceMods/src.clm/{restFileMod,SoilTemperatureMod,...}.F90
#   existing SP dumps: $DUMP_OUT/pdump_*_n{11881..11900,13458..13470}.nc
CLM_ROOT="${CLM_ROOT:-/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/installs/clm}"
CASE_DIR="${CASE_DIR:-$CLM_ROOT/cases/symfluence_build}"
DUMP_OUT="${DUMP_OUT:-/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_parity_run}"
# ⚠️ SAFETY: the existing validated SP dumps in $DUMP_OUT are the live :parity oracle
# reference — do NOT let a new run overwrite them. Each config writes its run output to
# an ISOLATED RUNDIR and collects only into a per-config REF_DIR; the parity oracle is
# pointed at REF_DIR explicitly, never at a clobbered $DUMP_OUT.
RUNDIR="${RUNDIR:-$DUMP_OUT/gen_$CONFIG/run}"
REF_DIR="$DUMP_OUT/refs_$CONFIG"

# --- per-config namelist overrides (user_nl_clm) ---------------------------------
# Each new config flips exactly the flags whose Fortran-vs-Julia parity we want to
# anchor. Extend this case as new reference paths are added.
declare -a NL_OVERRIDES
case "$CONFIG" in
  crop)    NL_OVERRIDES=( "use_crop=.true." "use_cn=.true." ) ;;
  lch4)    NL_OVERRIDES=( "use_lch4=.true." "use_cn=.true." ) ;;
  lake)    NL_OVERRIDES=( )  ;;   # driven by a deep-lake surfdata, not a flag
  glacier) NL_OVERRIDES=( )  ;;   # driven by an istice (PCT_GLACIER=100) surfdata
  cn)      NL_OVERRIDES=( "use_cn=.true." ) ;;
  sp)      NL_OVERRIDES=( ) ;;
  *) echo "unknown config '$CONFIG' — add its user_nl_clm overrides here" >&2; exit 2 ;;
esac

echo ">> generating Fortran refs: config=$CONFIG steps=$NSTEP_START..$((NSTEP_START+NSTEP_COUNT-1))"
echo ">> case=$CASE_DIR  ->  refs=$REF_DIR"

if [ ! -d "$CASE_DIR" ]; then
  cat >&2 <<EOF
ERROR: case dir not found: $CASE_DIR
This script must run on the machine with the instrumented CTSM build. Create the
case from the Bow parity case (the one that produced the existing \$DUMPDIR dumps),
apply the SourceMods restFileMod dump instrumentation, then re-run. See the banked
recipe 'fortran-clm-short-run-recipe' for the exact case.build / case.submit steps.
EOF
  exit 1
fi

# --- apply namelist overrides ----------------------------------------------------
NL="$CASE_DIR/user_nl_clm"
cp -f "$NL" "$NL.bak.$CONFIG" 2>/dev/null || true
for kv in "${NL_OVERRIDES[@]:-}"; do
  [ -z "$kv" ] && continue
  key="${kv%%=*}"
  grep -q "^$key" "$NL" && sed -i.tmp "s|^$key.*|$kv|" "$NL" || echo "$kv" >> "$NL"
done

# --- run the short window that brackets the requested dump steps ------------------
# The SourceMods dump fires on nstep ∈ [NSTEP_START, NSTEP_START+NSTEP_COUNT). We run
# from the spun-up restart so the bracket lands on the right calendar steps (see the
# recipe for stop_n/continue-mode caveats). RUNDIR is ISOLATED so we never clobber the
# existing validated SP dumps in $DUMP_OUT.
mkdir -p "$RUNDIR"
( cd "$CASE_DIR"
  ./xmlchange RUNDIR="$RUNDIR"
  ./xmlchange STOP_OPTION=nsteps,STOP_N=$((NSTEP_COUNT+1))
  # if a SourceMods .F90 changed, rebuild first: ./case.build
  ./case.submit --no-batch )

# --- collect the dumps into a per-config ref dir (from the isolated RUNDIR) -------
mkdir -p "$REF_DIR"
moved=0
for n in $(seq "$NSTEP_START" $((NSTEP_START+NSTEP_COUNT-1))); do
  for b in before_step after_canopyfluxes after_soiltemperature \
           after_soilfluxes after_hydrologynodrainage after_hydrologydrainage; do
    f="$RUNDIR/pdump_${b}_n${n}.nc"
    [ -f "$f" ] && { cp -f "$f" "$REF_DIR/"; moved=$((moved+1)); }
  done
done
echo ">> collected $moved dump files into $REF_DIR"
echo ">> point the harness :parity oracle at these by extending PARITY_NSTEP / a"
echo "   per-config dump dir in scripts/validation/validate.jl (REF_DIR=$REF_DIR)."
