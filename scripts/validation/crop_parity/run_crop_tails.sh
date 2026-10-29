#!/bin/bash
# Run the instrumented crop-parity cesm.exe in a staged run dir (see setup_crop_tails_run.sh).
#   RUN=<run dir>  CROPDUMP_WINDOWS="lo:hi,lo:hi"  run_crop_tails.sh
# mpiuni build: run the exe directly (no mpirun, not under lldb). MallocNanoZone=0 +
# retry is the macOS-26 xzone SIGTRAP workaround from scripts/validation/fortran_pdump.
set -uo pipefail
SD=${SYMFLUENCE_DATA:?set SYMFLUENCE_DATA}
RUN=${RUN:?set RUN}
EXE=${EXE:-$SD/installs/clm/cases/crop_parity_macos/bld/cesm.exe}
export MallocNanoZone=0
export HWLOC_COMPONENTS=-opencl
export CROPDUMP_WINDOWS=${CROPDUMP_WINDOWS:-}
cd "$RUN"
ok=0
for attempt in 1 2 3; do
    echo "=== attempt $attempt CROPDUMP_WINDOWS=$CROPDUMP_WINDOWS ===" >> cesm_run.log
    "$EXE" >> cesm_run.log 2>&1
    rc=$?
    echo "EXIT_CODE=$rc (attempt $attempt)" >> cesm_run.log
    if [ $rc -eq 0 ]; then ok=1; break; fi
done
echo "FINAL_OK=$ok"
