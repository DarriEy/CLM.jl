#!/usr/bin/env bash
# =============================================================================
# regen_parity_references.sh — regenerate the Fortran ground-truth dumps that
# the scripts/fortran_parity_*.jl harnesses SKIP when their reference is not
# staged on disk. Banked 2026-07-24 after the shared SYMFLUENCE_data rundirs
# (clm_parity_run, clm_cn_coldstart, clm_aripuana_parity, clm_ref_water_summer)
# were pruned, leaving only the surviving restarts.
#
# KEY RECIPE FACTS (all validated 2026-07-24):
#   * The prebuilt exe already carries the pdump instrumentation
#     (PDUMP_NSTEP_LO/HI env, restFile_write_dump) AND the WaterFluxBulkType
#     runoff-chain hist registrations — NO REBUILD NEEDED.
#     CASE=.../installs/clm/cases/symfluence_build/bld/cesm.exe
#   * The macOS-26 startup SIGTRAP (exit 133, no CLM output) is FLAKY: set
#     MallocNanoZone=0 and simply RETRY the exe in a loop until it clears
#     startup (~1 in 3 attempts succeeds). HWLOC_COMPONENTS=-opencl too.
#   * A CTSM `continue`/`branch` run needs the .rh0 restart-history file AND the
#     partial history tape it references. To AVOID that, use start_type=startup
#     with finidat = the surviving clm2.r restart and use_init_interp=.false.:
#     it reads prognostic state only, resets nstep to 0, starts fresh history,
#     and reproduces the physics step from the injected IC (symmetric to what
#     the Julia harness does). Deterministic: a startup+finidat run from the
#     2202-01-01 BGC restart reproduced the ORIGINAL summer trajectory to
#     T_GRND 4e-4 K over 4721 steps (validated vs the surviving before_step dump).
#   * pdump `before_step_nN` = state at top of step; pair it with the SAME step's
#     `after_hydrologydrainage_nN`. Relabel by DATE to the nstep the harness
#     expects (harness nstep 8761 <-> 2003-01-01 00:00; +1 h per +1 nstep).
#   * The current in-exe pdumpMod writes SMP/HK/KSOILROOT but NOT the newer
#     WATSAT_P/BSW_P/SUCSAT_P param fields — so fortran_parity_soilparam.jl's
#     direct-param block is skipped; regenerating those needs a rebuild.
#
# This script is documentation-as-code: paths are the 2026-07 box's; adjust and
# run the blocks you need. It is intentionally not a turnkey driver.
# =============================================================================
set -u
SD=/Users/darri.eythorsson/compHydro/SYMFLUENCE_data
CASE="$SD/installs/clm/cases/symfluence_build"
EXE="$CASE/bld/cesm.exe"
FE="$SD/domain_Bow_at_Banff_lumped/optimization/CLM/dds_run_1/final_evaluation"

run_with_retry() {  # $1 = rundir ; env PDUMP_NSTEP_LO/HI already exported
  cd "$1"
  source "$CASE/.env_mach_specific.sh"
  export MallocNanoZone=0 HWLOC_COMPONENTS=-opencl
  for try in $(seq 1 40); do
    "$EXE" > cesm.run.log 2>&1
    [ "$(wc -l < cesm.run.log)" -gt 25 ] && { echo "cleared startup on try $try"; return 0; }
  done
  echo "ALL 40 attempts SIGTRAP'd at startup"; return 1
}

# --- Bow SP family: smoke / soilparam / snow / validate (DUMPDIR=clm_parity_run)
#     startup+finidat from the surviving 2003-01-01 SP restart; use n0 pair,
#     stage as pdump_{before_step,after_hydrologydrainage}_n8761.nc.
bow_sp() {
  R=/private/tmp/claude-501/sp_parity_run; rm -rf "$R"; mkdir -p "$R"; cd "$R"
  for f in lnd_in datm_in drv_in drv_flds_in nuopc.runconfig nuopc.runseq datm.streams.xml fd.yaml; do cp "$FE/$f" .; done
  sed -i '' "s|finidat = ''|finidat = '$FE/Bow_at_Banff_lumped.clm2.r.2003-01-01-00000.nc'|" lnd_in
  sed -i '' 's/start_type = startup/start_type = startup/; s/start_ymd = 20020101/start_ymd = 20030101/' nuopc.runconfig
  sed -i '' 's/stop_option = nyears/stop_option = nsteps/; s/stop_n = 8/stop_n = 6/' nuopc.runconfig
  PDUMP_NSTEP_LO=0 PDUMP_NSTEP_HI=4 run_with_retry "$R" || return 1
  D="$SD/clm_parity_run"; mkdir -p "$D"
  cp pdump_before_step_n0.nc "$D/pdump_before_step_n8761.nc"
  cp pdump_after_hydrologydrainage_n0.nc "$D/pdump_after_hydrologydrainage_n8761.nc"
}

# --- CN cold start: cn_coldstart (clm_cn_coldstart/pdump_before_step_n0.nc)
#     true cold start (finidat=''), use_cn=true, start 2002-01-01, nstep 0.
cn_coldstart() {
  R=/private/tmp/claude-501/cn_coldstart_run; rm -rf "$R"; mkdir -p "$R"; cd "$R"
  for f in datm_in drv_in drv_flds_in nuopc.runconfig nuopc.runseq datm.streams.xml fd.yaml; do cp /private/tmp/claude-501/sp_parity_run/$f .; done
  cp "$SD/clm_bgc_spinup/bgc_ref_ndep_summer/lnd_in.used" lnd_in
  sed -i '' "s|finidat = .*|finidat = ' '|" lnd_in
  sed -i '' 's/start_ymd = 20030101/start_ymd = 20020101/' nuopc.runconfig
  PDUMP_NSTEP_LO=0 PDUMP_NSTEP_HI=2 run_with_retry "$R" || return 1
  D="$SD/clm_cn_coldstart"; mkdir -p "$D"; cp pdump_before_step_n0.nc "$D/pdump_before_step_n0.nc"
}

# --- Aripuana cold start: aripuana_coldstart (clm_aripuana_parity/*.clm2.r.*)
#     SP cold start from the Drive Aripuana namelists, restart EVERY step x24.
#     Also recreate the domain_Aripuana_Amazon tree the harness hard-codes as
#     symlinks -> scratch aripuana_data (params) + forcing.
# --- qflx_surf: h1 instantaneous flux tape at the 2202-07-16 summer window.
#     startup+finidat from the 2202-01-01 BGC restart, +~4925 steps, tape 2
#     (hist_nhtfrq=-8760,1 ; avgflag 'A','I' ; fincl2 = the runoff chain), then
#     slice records for nstep 4721..4920 into
#     clm_ref_water_summer/...h1.2202-07-16-61200.nc (record 1 = ref nstep 1757873).
#     The before_step injection dumps already survive in bgc_ref_ndep_summer/.

echo "Recipe reference only — call bow_sp / cn_coldstart etc. as needed."
