#!/bin/bash
# CNSoyfix reference: CONTINUE a crop-tails run from an edited (N-limited) restart.
#
#   SPIN=<run dir that wrote the restart>  STAMP=2017-07-24-14400
#   EDITED=<restart from make_nlimited_restart.jl>  RUN=<new run dir>  NDAYS=6
#
# A `continue` run (not startup+finidat) keeps every accumulator and the step
# counter, so the only difference from the unperturbed trajectory is the injected
# soil-N state. (A startup+finidat run resets plant N demand to ~0 and it takes
# days to recover; observed 2026-10-28.) Needs the matching .rh0/.rh1, cpl and datm
# restarts, which the spin-up run writes alongside the .r file.
#
# FLEXIBLECN=.false. (default here) swaps NutrientCompetitionFlexibleCNMod for the
# CLM4.5-default competition. In ctsm5.3.012 the FlexibleCN non-FUN plant N demand
# is UNCONDITIONALLY the Michaelis-Menten uptake form
#     Vmax_N * frootc * sminn/(sminn+1) * t_scalar * nscalar
# (mm_nuptake_opt is not consulted on that branch), so demand collapses with the
# mineral pool, demand*dt < supply structurally, and fpg < 1 -- the CNSoyfix gate
# -- is essentially unreachable: zeroing the soil N drove plant_ndemand to 0.0 on
# the next step and fpg stayed 1 (observed 2026-10-28, with and without
# mm_nuptake_opt). The CLM4.5 default demand is availc*n_allometry/c_allometry.
set -euo pipefail
SPIN=${SPIN:?}; STAMP=${STAMP:?}; EDITED=${EDITED:?}; RUN=${RUN:?}; NDAYS=${NDAYS:-6}
FLEXIBLECN=${FLEXIBLECN:-.false.}
mkdir -p "$RUN/init_generated_files"
for f in lnd_in datm_in drv_in drv_flds_in nuopc.runconfig nuopc.runseq fd.yaml datm.streams.xml; do
    cp "$SPIN/$f" "$RUN/$f"
done
cp "$SPIN"/Crop_tails.clm2.rh*."$STAMP".nc "$SPIN"/Crop_tails.cpl.r."$STAMP".nc "$SPIN"/Crop_tails.datm.r."$STAMP".nc "$RUN/"
cp "$SPIN"/Crop_tails.clm2.h[0-9].*.nc "$RUN/"   # continue appends to the open history tapes
cp "$EDITED" "$RUN/Crop_tails.clm2.r.$STAMP.nc"
echo "./Crop_tails.clm2.r.$STAMP.nc" > "$RUN/rpointer.lnd"
echo "Crop_tails.cpl.r.$STAMP.nc" > "$RUN/rpointer.cpl"
echo "Crop_tails.datm.r.$STAMP.nc" > "$RUN/rpointer.atm"
cd "$RUN"
sed -i '' -e "s|^     start_type = .*|     start_type = continue|" \
          -e "s|^     stop_option = .*|     stop_option = ndays|" -e "s|^     stop_n = .*|     stop_n = $NDAYS|" \
          -e "s|^     restart_option = .*|     restart_option = never|" nuopc.runconfig
sed -i '' -e "s|^ use_flexiblecn = .*| use_flexiblecn = $FLEXIBLECN|" lnd_in
grep -n "use_flexiblecn" lnd_in
echo "staged continue run $RUN from $EDITED"
