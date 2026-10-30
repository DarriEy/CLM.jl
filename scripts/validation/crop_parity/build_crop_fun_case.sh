#!/bin/bash
# Build the instrumented FUN crop-parity CTSM case (soybean fixation + FUN N fluxes).
#
#   CASE = the CIME case (default: $SYMFLUENCE_DATA/installs/clm/cases/crop_fun_macos)
#
# A separate case (own bld/cesm.exe) so the shared crop_parity_macos binary is
# untouched. It carries the crop_parity SourceMods (CNDriverMod hooks + cropdumpMod)
# PLUS a SoilBiogeochemCompetitionMod SourceMod that dumps every CNFUN input and
# output around both CNFUN call sites ('before_cnfun' / 'after_cnfun').
#
# One-time case creation (clone of the crop_parity build), then this script:
#   PATH=/usr/bin:/bin:$CLM/toolchain/bin $CLM/ctsm5.3.012/cime/scripts/create_clone \
#       --clone $CLM/cases/crop_parity_macos --case $CASE --cime-output-root $CASE/output
#   cd $CASE && ./xmlchange EXEROOT=$CASE/bld,RUNDIR=$CASE/run && ./case.setup
# (CREATE=1 does both.)  Environment traps: see build_crop_case.sh.
set -euo pipefail
SD=${SYMFLUENCE_DATA:?set SYMFLUENCE_DATA}
CLM=$SD/installs/clm
CASE=${CASE:-$CLM/cases/crop_fun_macos}
HERE=$(cd "$(dirname "$0")" && pwd)
export PATH=$CLM/cases/bow_btran_macos/.runtime-bin:/usr/bin:/bin:$CLM/toolchain/bin:/opt/homebrew/bin:/usr/sbin:/sbin

if [ "${CREATE:-0}" = 1 ] && [ ! -d "$CASE" ]; then
    "$CLM/ctsm5.3.012/cime/scripts/create_clone" --clone "$CLM/cases/crop_parity_macos" \
        --case "$CASE" --cime-output-root "$CASE/output"
    (cd "$CASE" && ./xmlchange EXEROOT="$CASE/bld",RUNDIR="$CASE/run" && ./case.setup)
fi

python3 "$HERE/apply_cndriver_hooks.py" "$CLM/ctsm5.3.012/src/biogeochem/CNDriverMod.F90" \
    "$CASE/SourceMods/src.clm/CNDriverMod.F90"
python3 "$HERE/apply_competition_fun_hooks.py" \
    "$CLM/ctsm5.3.012/src/soilbiogeochem/SoilBiogeochemCompetitionMod.F90" \
    "$CASE/SourceMods/src.clm/SoilBiogeochemCompetitionMod.F90"
cp "$HERE/cropdumpMod.F90" "$CASE/SourceMods/src.clm/"
cp "$HERE/gnu_homebrew.cmake" "$CASE/cmake_macros/gnu_homebrew.cmake"
rm -f "$CASE"/bld/gnu/mpich/nodebug/nothreads/nuopc/nuopc/esmf/clm/obj/{Srcfiles,Depends}

export ESMFMKFILE=$CLM/toolchain/lib/esmf.mk
export NETCDF=$CLM/toolchain NETCDF_PATH=$CLM/toolchain PNETCDF_PATH=$CLM/toolchain
cd "$CASE"
./case.build
set +o pipefail
if strings bld/cesm.exe | grep before_cnfun >/dev/null; then echo "fun hooks present in bld/cesm.exe"; else echo "ERROR: fun hooks missing (stale exe?)"; exit 1; fi
