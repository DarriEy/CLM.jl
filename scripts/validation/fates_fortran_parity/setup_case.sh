#!/usr/bin/env bash
# =============================================================================
# setup_case.sh — stand up a FATES-enabled Fortran CTSM case that emits a
# cold-start (nstep=0) schema-v1 pdump, as the bit-parity ground truth for the
# Julia CLM.jl FATES port (scripts/fates_fortran_parity.jl).
#
# This script performs every step that has been VERIFIED to work against the
# local CTSM install (ctsm5.3.012), plus installs the Fortran instrumentation
# SourceMod. It stops short of `case.build`/`case.submit`, which are gated on
# atmospheric forcing that is NOT staged in this environment (see the BLOCKER
# section at the bottom and README.md).
#
# Verified working here: create_newcase, case.setup, user_nl_* config,
# preview_namelists (all namelists generate), ncgen of the FATES param file,
# and SourceMods instrumentation install. The single hard wall is DATM forcing.
# =============================================================================
set -euo pipefail

# ---- paths (edit for your environment) --------------------------------------
CLMROOT=${CLMROOT:-/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/installs/clm}
INPUTDATA=${INPUTDATA:-/Users/darri.eythorsson/projects/cesm-inputdata}
CASE=${CASE:-$CLMROOT/cases/fates_parity_ne3}
COMPSET=${COMPSET:-I2000Clm60Fates}     # 2000_DATM%GSWP3v1_CLM60%FATES_...
RES=${RES:-ne3_ne3_mg37}                 # lnd=ne3np4 — the only staged surfdata/mesh
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---- 1. FATES parameter file: ncgen from the shared, byte-identical .cdl -----
# data/fates/fates_params_default.cdl (port) == src/fates/parameter_files/...cdl
# (install), verified byte-identical, so cold-start cohorts need NO param tuning.
CDL="$CLMROOT/src/fates/parameter_files/fates_params_default.cdl"
FATES_PARAM="$INPUTDATA/lnd/clm2/paramdata/fates_params_default.nc"
ncgen -o "$FATES_PARAM" "$CDL"
echo "[ok] generated FATES param nc -> $FATES_PARAM"

# ---- 2. create the case -----------------------------------------------------
rm -rf "$CASE"
"$CLMROOT/cime/scripts/create_newcase" --case "$CASE" --compset "$COMPSET" \
    --res "$RES" --run-unsupported --machine homebrew --compiler gnu
cd "$CASE"

# ---- 3. cold start, 1 step, single serial task (ESMF is built mpiuni) -------
./xmlchange STOP_N=1,STOP_OPTION=nsteps,RUN_STARTDATE=2000-01-01
./xmlchange NTASKS=1,NTASKS_ESP=1     # mpiuni ESMF => run as one PET (NTASKS>1 aborts)
./case.setup

# ---- 4. namelist: point FATES param + surfdata at staged files --------------
cat >> user_nl_clm <<EOF
fates_paramfile = '\$DIN_LOC_ROOT/lnd/clm2/paramdata/fates_params_default.nc'
fsurdat = '\$DIN_LOC_ROOT/lnd/clm2/surfdata_esmf/ctsm5.4.0/surfdata_ne3np4_hist_2000_78pfts_c251022.nc'
EOF
./preview_namelists     # verified: all component namelists generate cleanly

# ---- 5. install the FATES cold-start pdump instrumentation ------------------
SM="$CASE/SourceMods/src.clm"
cp "$HERE/FatesParityDumpMod.F90" "$SM/FatesParityDumpMod.F90"
python3 - "$CLMROOT/src/utils/clmfates_interfaceMod.F90" "$SM/clmfates_interfaceMod.F90" <<'PY'
import sys
src, dst = sys.argv[1], sys.argv[2]
lines = open(src).read().splitlines(keepends=True)
out, used, called = [], False, False
for line in lines:
    if (not used) and 'use EDInitMod' in line and 'set_site_properties' in line:
        out.append(line)
        out.append("   use FatesParityDumpMod    , only : fates_parity_dump   ! FATES parity instrumentation\n")
        used = True; continue
    if (not called) and "t_stopf('fates_initcoldstart')" in line:
        out += [
          "     ! ---- FATES cold-start parity instrumentation (schema-v1 pdump) ------------\n",
          "     do nc = 1, nclumps\n",
          "        if ( this%fates(nc)%nsites>0 ) then\n",
          "           call fates_parity_dump(this%fates(nc)%sites, this%fates(nc)%nsites, &\n",
          "                get_nstep(), 'fates_pdump_fortran_n0.txt')\n",
          "        end if\n",
          "     end do\n",
          "     ! --------------------------------------------------------------------------\n",
        ]
        called = True
    out.append(line)
open(dst, 'w').writelines(out)
assert used and called, (used, called)
print("[ok] instrumented clmfates_interfaceMod.F90 (use + call injected)")
PY

echo
echo "=== case created + configured + instrumented ==================================="
echo "BUILD NOTE: the homebrew mach config hardcodes a netcdf Cellar version. If the"
echo "  final link fails with \"ld: library 'netcdf' not found\", point it at the"
echo "  installed one, e.g.:  ln -sfn \"\$(brew --prefix netcdf)\" /opt/homebrew/Cellar/netcdf/<old-ver>"
echo "  (runtime netcdf loads fine via the dylib install-name; the symlink is link-time only.)"
echo
echo "NEXT (gated on DATM forcing — see README.md BLOCKER):"
echo "  ./check_input_data --download   # needs the GSWP3 stream files + meshes"
echo "  ./case.build"
echo "  ./case.submit                   # run cesm.exe DIRECTLY, NOT under lldb (macOS-26)"
echo "  julia +1.12 --project=scripts scripts/fates_fortran_parity.jl compare \\"
echo "      \$(./xmlquery -value RUNDIR)/fates_pdump_fortran_n0.txt"
