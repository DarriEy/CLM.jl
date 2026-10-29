#!/bin/bash
# Regenerate the bow_ref_ch4o2 reference: the bow_ref_ch4tws case (Bow, use_lch4,
# TWS_inversion, CN, nitrif_denitrif) rebuilt with the O2-demand bgcdump block
# (scripts/validation/fortran_sourcemods/bgcdumpMod.F90). docs/CH4_FIRE_PARITY.md §17.
#
#   SYMFLUENCE_DATA=... STEP=case|build|run|stage  build_and_run_ch4o2.sh
#
# case : clone an isolated CIME case (never touches the shared cesm.exe):
#          $CLM/cases/ch4o2_macos   (clone of crop_parity_macos, which builds on this box)
#        SourceMods = case-archive/symfluence_build/SourceMods/src.clm (the instrumentation
#        that produced bow_ref_ch4tws: clm_driver hooks, restFile_write_dump, public ch4Mod)
#        with bgcdumpMod.F90 replaced by the repo copy.
# build: ./case.build (~10 min from scratch).
# run  : stage a relocated copy of bow_ref_ch4tws's namelists in $RUN, run n4720..n4888.
# stage: copy the dumps to $SYMFLUENCE_DATA/clm_bgc_spinup/bow_ref_ch4o2.
#
# Verified 2026-10-29: every field shared with bow_ref_ch4tws is bit-identical.
set -euo pipefail
SD=${SYMFLUENCE_DATA:?set SYMFLUENCE_DATA}
CLM=$SD/installs/clm
CASE=${CASE:-$CLM/cases/ch4o2_macos}
RUN=${RUN:-$SD/runs/ch4o2_ref}
HERE=$(cd "$(dirname "$0")/.." && pwd)
ID=$SD/installs/cesm-inputdata
TPL=$SD/clm_bgc_spinup/bow_ref_ch4tws
export ESMFMKFILE=$CLM/toolchain/lib/esmf.mk
export NETCDF=$CLM/toolchain NETCDF_PATH=$CLM/toolchain PNETCDF_PATH=$CLM/toolchain
# CIME needs python <= 3.11 and /usr/bin/perl (XML::LibXML) ahead of the toolchain.
export PATH=$CLM/cases/bow_btran_macos/.runtime-bin:/usr/bin:/bin:$CLM/toolchain/bin:/opt/homebrew/bin:/usr/sbin:/sbin

case ${STEP:?set STEP=case|build|run|stage} in
case)
    "$CLM/ctsm5.3.012/cime/scripts/create_clone" --clone "$CLM/cases/crop_parity_macos" \
        --case "$CASE" --cime-output-root "$CASE/output"
    cd "$CASE" && ./xmlchange EXEROOT="$CASE/bld",RUNDIR="$CASE/run" && ./case.setup
    SM=$CASE/SourceMods/src.clm; A=$CLM/case-archive/symfluence_build/SourceMods/src.clm
    rm -f "$SM"/CNDriverMod.F90 "$SM"/cropdumpMod.F90
    for f in CanopyFluxesMod ch4Mod clm_driver LunaMod pdumpMod PhotosynthesisMod restFileMod \
             SnowCoverFractionSwensonLawrence2012Mod WaterFluxBulkType; do cp "$A/$f.F90" "$SM/"; done
    cp "$HERE/fortran_sourcemods/bgcdumpMod.F90" "$SM/"
    cp "$HERE/crop_parity/gnu_homebrew.cmake" "$CASE/cmake_macros/gnu_homebrew.cmake" ;;
build)
    rm -f "$CASE"/bld/gnu/mpich/nodebug/nothreads/nuopc/nuopc/esmf/clm/obj/{Srcfiles,Depends}
    cd "$CASE" && ./case.build
    strings bld/cesm.exe | grep -q O2_DEMAND_RR && echo "O2-demand bgcdump present" ;;
run)
    mkdir -p "$RUN/init_generated_files"
    for f in lnd_in datm_in drv_in drv_flds_in nuopc.runconfig nuopc.runseq fd.yaml datm.streams.xml; do
        cp "$TPL/$f" "$RUN/$f"; done
    cd "$RUN"
    OLD_SD=/Users/darri.eythorsson/compHydro/SYMFLUENCE_data
    OLD_ID=/Users/darri.eythorsson/projects/cesm-inputdata
    for f in lnd_in datm_in drv_in drv_flds_in nuopc.runconfig datm.streams.xml; do
        sed -i '' -e "s|$OLD_SD/installs/cesm-inputdata|$ID|g" -e "s|$OLD_ID|$ID|g" -e "s|$OLD_SD|$SD|g" "$f"; done
    sed -i '' -e "s|^     stop_n = .*|     stop_n = 4890|" -e "s|^     restart_n = .*|     restart_n = 99999|" nuopc.runconfig
    export MallocNanoZone=0 HWLOC_COMPONENTS=-opencl
    export PDUMP_NSTEP_LO=4720 PDUMP_NSTEP_HI=4888 BGCDUMP_NSTEP_LO=4720 BGCDUMP_NSTEP_HI=4888
    for attempt in 1 2 3; do                          # macOS-26 xzone SIGTRAP retry
        "$CASE/bld/cesm.exe" >> cesm_run.log 2>&1 && { echo "run OK"; break; }
    done ;;
stage)
    D=$SD/clm_bgc_spinup/bow_ref_ch4o2; mkdir -p "$D"
    cp "$RUN"/pdump_before_step_n*.nc "$RUN"/pdump_after_hydrologydrainage_n*.nc \
       "$RUN"/bgcdump_after_ch4_n*.nc "$RUN"/{lnd_in,datm_in,drv_in,drv_flds_in,nuopc.runconfig,nuopc.runseq,datm.streams.xml} "$D/" ;;
esac
