#!/bin/bash
# Build the instrumented crop-parity CTSM case (macOS replacement box).
#
#   CASE   = the CIME case (default: $SYMFLUENCE_DATA/installs/clm/cases/crop_parity_macos)
#
# One-time case creation (clone of the working SP build; CTSM compiles the CN/crop
# code unconditionally, so use_cn/use_crop are runtime namelist flags):
#   PATH=/usr/bin:/bin:$CLM/toolchain/bin $CLM/ctsm5.3.012/cime/scripts/create_clone \
#       --clone $CLM/cases/bow_btran_macos --case $CASE --cime-output-root $CASE/output
#   cd $CASE && ./xmlchange EXEROOT=$CASE/bld,RUNDIR=$CASE/run && ./case.setup
#
# Environment traps on this box (each cost a failed build):
#   * CIME needs python <= 3.11 (distutils): /usr/bin/python3 (3.9) works, the
#     homebrew 3.15 does not.
#   * build-namelist needs XML::LibXML: /usr/bin/perl has it, the toolchain perl does not.
#     => /usr/bin MUST precede the toolchain bin on PATH.
#   * ESMFMKFILE is not set by env_mach_specific.xml; point it at the toolchain.
set -euo pipefail
SD=${SYMFLUENCE_DATA:?set SYMFLUENCE_DATA}
CLM=$SD/installs/clm
CASE=${CASE:-$CLM/cases/crop_parity_macos}
HERE=$(cd "$(dirname "$0")" && pwd)

python3 "$HERE/apply_cndriver_hooks.py" "$CLM/ctsm5.3.012/src/biogeochem/CNDriverMod.F90" \
    "$CASE/SourceMods/src.clm/CNDriverMod.F90"
cp "$HERE/cropdumpMod.F90" "$CASE/SourceMods/src.clm/"
# create_clone does not carry the case's patched cmake macro: without
# -DNETCDF_ENABLE_LEGACY_MACROS the bundled PIO fails on `_FillValue`.
cp "$HERE/gnu_homebrew.cmake" "$CASE/cmake_macros/gnu_homebrew.cmake"
# a NEW SourceMods module needs the dependency lists regenerated
rm -f "$CASE"/bld/gnu/mpich/nodebug/nothreads/nuopc/nuopc/esmf/clm/obj/{Srcfiles,Depends}

export ESMFMKFILE=$CLM/toolchain/lib/esmf.mk
export NETCDF=$CLM/toolchain NETCDF_PATH=$CLM/toolchain PNETCDF_PATH=$CLM/toolchain   # Macros: -Wl,-rpath $(NETCDF)/lib
export PATH=$CLM/cases/bow_btran_macos/.runtime-bin:/usr/bin:/bin:$CLM/toolchain/bin:/opt/homebrew/bin:/usr/sbin:/sbin
cd "$CASE"
./case.build
# A 2 s "successful" clm build means nothing recompiled; prove the hooks are in:
set +o pipefail
if strings bld/cesm.exe | grep after_cropalloc >/dev/null; then echo "cropdump hooks present in bld/cesm.exe"; else echo "ERROR: hooks missing (stale exe?)"; exit 1; fi
