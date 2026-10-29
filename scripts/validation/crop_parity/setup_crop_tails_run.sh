#!/bin/bash
# Stage a self-contained CTSM run directory for the crop-tail parity reference
# (grain fill + grain fluxes, soybean N fixation, vernalization).
#
# Derived from the 2026-07-19 cold non-FUN crop case
# ($SYMFLUENCE_DATA/clm_bgc_spinup/crop_ref_nofun_cold, docs/CROP_LIFECYCLE_MAP.md §3e),
# with paths relocated to this box and three changes:
#   * fsurdat = the USplains crop surfdata with the foddergrass CFT (pft 45, which
#     CTSM merges into spring wheat) re-assigned to WINTER WHEAT (pft 21) -- the
#     ne3np4 78-pft surfdata has PCT_CFT == 0 for every winter cereal (21/22/27/28/
#     31/32) in every gridcell, so no real winter-wheat cell exists to extract.
#     Built by scripts/validation/crop_parity/make_winterwheat_surfdata.jl.
#   * use_fun = .false. (CNSoyfix only runs without FUN), use_fertilizer = .true.
#   * cold start (default) STOP_N STOP_OPTION from 2016-01-01 (forcing 2016-2017 cycled).
#
# Usage: RUN=<dir under $SYMFLUENCE_DATA/runs> [STOP_OPTION=nyears STOP_N=3 START_YMD=.. START_TOD=..
#        RESTART_OPTION=.. RESTART_N=.. FINIDAT=..] setup_crop_tails_run.sh
set -euo pipefail
SD=${SYMFLUENCE_DATA:?set SYMFLUENCE_DATA}
ID=${CESM_INPUTDATA:-$SD/installs/cesm-inputdata}
RUN=${RUN:-$SD/runs/crop_tails_ref}
STOP_OPTION=${STOP_OPTION:-nyears}
STOP_N=${STOP_N:-3}
START_YMD=${START_YMD:-20160101}
START_TOD=${START_TOD:-14400}
RESTART_OPTION=${RESTART_OPTION:-nyears}
RESTART_N=${RESTART_N:-1}
FINIDAT=${FINIDAT:-}          # startup from a (possibly edited) restart -- the PR #312 recipe
FLEXIBLECN=${FLEXIBLECN:-.true.}   # .false. for the CNSoyfix reference (see setup_soyfix_continue.sh)
SURF=${SURF:-$SD/runs/crop_tails_inputs/surfdata_cropCFT_USplains_1pt_wwheat.nc}
TPL=$SD/clm_bgc_spinup/crop_ref_nofun_cold

mkdir -p "$RUN/init_generated_files"
for f in lnd_in datm_in drv_in drv_flds_in nuopc.runconfig nuopc.runseq fd.yaml datm.streams.xml; do
    cp "$TPL/$f" "$RUN/$f"
done
cd "$RUN"
OLD_SD=/Users/darri.eythorsson/compHydro/SYMFLUENCE_data
OLD_ID=/Users/darri.eythorsson/projects/cesm-inputdata
for f in lnd_in datm_in drv_in drv_flds_in nuopc.runconfig datm.streams.xml; do
    sed -i '' -e "s|$OLD_SD/installs/cesm-inputdata|$ID|g" -e "s|$OLD_ID|$ID|g" -e "s|$OLD_SD|$SD|g" "$f"
done
sed -i '' -e "s|^ fsurdat = .*| fsurdat = '$SURF'|" \
          -e "s|^ use_fun = .*| use_fun = .false.|" \
          -e "s|^ use_fertilizer = .*| use_fertilizer = .true.|" \
          -e "s|^ nrevsn = .*||" \
          -e "s|^ finidat = .*| finidat = '$FINIDAT'|" \
          -e "s|^ use_flexiblecn = .*| use_flexiblecn = $FLEXIBLECN|" lnd_in
# daily patch-level history to locate phase events (tape 2)
perl -i -pe 's/^ hist_mfilt = .*/ hist_mfilt = 20, 1000/; s/^ hist_nhtfrq = .*/ hist_nhtfrq = -8760, -24\n hist_dov2xy = .true., .false.\n hist_fincl2 = \x27CPHASE\x27,\x27HUI\x27,\x27SOYFIXN\x27,\x27FPG\x27,\x27TLAI\x27,\x27PLANT_NDEMAND\x27,\x27GDDHARV\x27,\x27LEAFC\x27/' lnd_in
sed -i '' -e "s|^     stop_n = .*|     stop_n = $STOP_N|" -e "s|^     stop_option = .*|     stop_option = $STOP_OPTION|" \
          -e "s|^     start_ymd = .*|     start_ymd = $START_YMD|" -e "s|^     start_tod = .*|     start_tod = $START_TOD|" \
          -e "s|^     restart_n = .*|     restart_n = $RESTART_N|" -e "s|^     restart_option = .*|     restart_option = $RESTART_OPTION|" \
          -e "s|^     case_name = .*|     case_name = Crop_tails|" nuopc.runconfig
echo "staged $RUN (fsurdat=$SURF, start=$START_YMD/$START_TOD, stop=$STOP_N $STOP_OPTION, finidat='$FINIDAT')"
grep -nE "fsurdat|use_fun|use_fertilizer|use_crop|finidat|hist_" lnd_in
