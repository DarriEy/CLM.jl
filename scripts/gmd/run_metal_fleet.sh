#!/bin/zsh
# ==========================================================================
# run_metal_fleet.sh — run the full Metal device-parity fleet for the GMD
# accelerator-qualification gate.
#
# Runs every scripts/gpu_validate_*.jl and scripts/gpu_ad_*.jl harness under
# `julia --project=scripts` with the Metal backend and records, per harness:
# name, exit status, wall-clock seconds, and the log path. Classification is
# STRUCTURAL (process exit status only) — never grep over log text; "0 failed"
# in a log has previously been misread as a failure by text classifiers.
#
# The whole fleet always runs; a dead or skipped harness must be visible as a
# row, not silently absent (the #258→#309 changed-files-sweep lesson).
#
# macOS quirks handled here (docs/HANDOFF_2026-07-24.md §7):
#  - no GNU timeout: per-harness background watchdog kill;
#  - spurious xzone SIGTRAP at startup: one retry with MallocNanoZone=0.
#
# Usage:
#   scripts/gmd/run_metal_fleet.sh <results_dir> [harness_glob_filter]
#
# Environment consumed (all optional):
#   SYMFLUENCE_DATA, CESM_INPUTDATA, CLM_FSURDAT, CLM_PARAMFILE — forwarded to
#   the data-resolving harnesses.
#   FLEET_TIMEOUT_S — per-harness wall limit, default 1800.
# Emits <results_dir>/fleet.jsonl (one JSON object per harness) and
# <results_dir>/logs/<name>.log.
# ==========================================================================
set -u
REPO_ROOT="${0:A:h:h:h}"
RESULTS_DIR="${1:?usage: run_metal_fleet.sh <results_dir> [filter]}"
FILTER="${2:-}"
TIMEOUT_S="${FLEET_TIMEOUT_S:-1800}"
mkdir -p "$RESULTS_DIR/logs"
JSONL="$RESULTS_DIR/fleet.jsonl"
: > "$JSONL"

export HWLOC_COMPONENTS=-opencl

cd "$REPO_ROOT"
harnesses=(scripts/gpu_validate_*.jl(N) scripts/gpu_ad_*.jl(N))
if [[ -n "$FILTER" ]]; then
  harnesses=(${(M)harnesses:#*${FILTER}*})
fi
echo "fleet: ${#harnesses[@]} harnesses, timeout ${TIMEOUT_S}s each"

run_one() {
  local script="$1" log="$2" env_prefix="$3"
  # Watchdog: kill the harness if it exceeds TIMEOUT_S (no GNU timeout on macOS).
  ( eval "$env_prefix" julia --project=scripts "$script" >"$log" 2>&1 ) &
  local pid=$!
  ( sleep "$TIMEOUT_S" && kill -9 "$pid" 2>/dev/null ) &
  local watchdog=$!
  wait "$pid" 2>/dev/null
  local rc=$?
  kill "$watchdog" 2>/dev/null; wait "$watchdog" 2>/dev/null
  return $rc
}

for script in "${harnesses[@]}"; do
  name="${${script:t}%.jl}"
  log="$RESULTS_DIR/logs/${name}.log"
  start=$(date +%s)
  run_one "$script" "$log" ""
  rc=$?
  # Retry once with MallocNanoZone=0 on the startup SIGTRAP quirk (exit 133).
  retried=false
  if [[ $rc -eq 133 ]]; then
    retried=true
    run_one "$script" "$log" "MallocNanoZone=0"
    rc=$?
  fi
  end=$(date +%s)
  timed_out=false
  [[ $rc -eq 137 ]] && timed_out=true
  printf '{"harness":"%s","exit_status":%d,"wall_s":%d,"timed_out":%s,"retried_nanozone":%s,"log":"logs/%s.log"}\n' \
    "$name" "$rc" $((end - start)) "$timed_out" "$retried" "$name" >> "$JSONL"
  echo "[$(date +%H:%M:%S)] $name exit=$rc wall=$((end - start))s"
done

pass=$(grep -c '"exit_status":0,' "$JSONL")
echo "fleet complete: $pass/${#harnesses[@]} exit-zero; results in $JSONL"
