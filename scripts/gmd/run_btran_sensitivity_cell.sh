#!/usr/bin/env bash
set -euo pipefail

# Execute one already-configured CTSM/CLM.jl Bow sensitivity cell on Linux.
# Case creation remains separate until the definitive CIME machine definition and
# licensed input bundle are frozen. Every path is explicit to prevent developer-path
# fallback from becoming paper evidence.

: "${CTSM_CASE:?set CTSM_CASE to the configured CTSM case directory}"
: "${CTSM_ROOT:?set CTSM_ROOT to the exact CTSM source checkout}"
: "${CTSM_RUN_DIR:?set CTSM_RUN_DIR to its isolated target run directory}"
: "${SYMFLUENCE_DATA:?set SYMFLUENCE_DATA to the staged input archive}"
: "${BTRAN_CELL_ID:?set BTRAN_CELL_ID (for example linux-gnu-release)}"

CLM_JL_ROOT=${CLM_JL_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}
BTRAN_RESULTS_DIR=${BTRAN_RESULTS_DIR:-$CLM_JL_ROOT/repro/results/btran-sensitivity}
BTRAN_NSTEP=${BTRAN_NSTEP:-2901}
CTSM_LAUNCH=${CTSM_LAUNCH:-"$CTSM_CASE/bld/cesm.exe"}

if [[ $(uname -s) != Linux ]]; then
  echo "error: definitive sensitivity cells must run on Linux" >&2
  exit 2
fi
for path in "$CTSM_CASE/case.build" "$CTSM_RUN_DIR" "$CLM_JL_ROOT/Project.toml"; do
  [[ -e $path ]] || { echo "error: missing required path: $path" >&2; exit 2; }
done
if [[ -n $(git -C "$CLM_JL_ROOT" status --porcelain=v1 --untracked-files=all) ]]; then
  echo "error: CLM.jl worktree must be clean" >&2
  exit 2
fi

mkdir -p "$BTRAN_RESULTS_DIR/$BTRAN_CELL_ID"
cell_dir=$BTRAN_RESULTS_DIR/$BTRAN_CELL_ID

(cd "$CTSM_CASE" && ./case.build) 2>&1 | tee "$cell_dir/ctsm-build.log"
(cd "$CTSM_RUN_DIR" && bash -lc "$CTSM_LAUNCH") 2>&1 | tee "$cell_dir/ctsm-run.log"

comparison_log=$cell_dir/comparison.log
env SYMFLUENCE_DATA="$SYMFLUENCE_DATA" julia --project="$CLM_JL_ROOT" \
  "$CLM_JL_ROOT/scripts/validation/compare_btran_shared_state.jl" \
  "$BTRAN_NSTEP" "$CTSM_RUN_DIR" 2>&1 | tee "$comparison_log"

julia_trace=$CTSM_RUN_DIR/btran_canopy_periter_julia_n${BTRAN_NSTEP}.txt
fortran_trace=$CTSM_RUN_DIR/btran_canopy_periter_fortran_n${BTRAN_NSTEP}.txt
python3 "$CLM_JL_ROOT/scripts/gmd/summarize_btran_oracle.py" \
  --cell-id "$BTRAN_CELL_ID" \
  --comparison-log "$comparison_log" \
  --julia-trace "$julia_trace" \
  --fortran-trace "$fortran_trace" \
  --ctsm-root "$CTSM_ROOT" \
  --ctsm-executable "$CTSM_LAUNCH" \
  --case-dir "$CTSM_CASE" \
  --output "$cell_dir/result.json"
julia --project="$CLM_JL_ROOT" "$CLM_JL_ROOT/scripts/gmd/record_environment.jl" \
  > "$cell_dir/environment.json"

echo "wrote $cell_dir/result.json"
