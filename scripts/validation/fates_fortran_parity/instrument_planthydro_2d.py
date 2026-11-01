"""instrument_planthydro_2d.py — per-iteration traces of the FATES 2-D hydraulics solvers.

Applied on top of instrument_planthydro.py's output; writes to fates_log (lnd.log):
  NWT ft dbh nsteps nwtn_iter nsd residual_amax sum(residual) dtime tm      (MatSolve2D)
  PIC ft dbh icnv picd_iter dtime tm wb_error max|dth| max|th-th_prev| cfl_max (PicardSolve2D)
the same quantities as the port's CLM._HYDRO2D_TRACE hook. Does not change the run.

  python3 instrument_planthydro.py $SRC/FatesPlantHydraulicsMod.F90 $OUT/a.F90
  python3 instrument_planthydro_2d.py $OUT/a.F90 $OUT/src/FatesPlantHydraulicsMod.F90
  julia build_patched_cesm.jl $OUT/exe SourceMods_planthydro/FatesAllometryMod.F90 \\
        $OUT/src/FatesPlantHydraulicsMod.F90
  EXE=$OUT/exe/cesm.exe CDL=<params with fates_hydro_solver = 2|3> bash stage_planthydro_run.sh
usage: python3 instrument_planthydro_2d.py <src.F90> <dst.F90>"""
import sys
src, dst = sys.argv[1:3]
t = open(src).read()

def rep(anchor, new, count=1):
    global t
    assert t.count(anchor) == count, (anchor[:80], t.count(anchor))
    t = t.replace(anchor, new)

# Newton: after the residual scan + iter cap check, log one line per iteration
rep("""         if ( nwtn_iter > max_newton_iter) then
            icnv = icnv_fail_round
            if(debug)then""",
"""         write(fates_log(),'(a,i4,es25.16e3,3i6,4es25.16e3)') 'NWT ', ft, cohort%dbh, nsteps, nwtn_iter, nsd, &
              residual_amax, sum(residual(:)), dtime, tm
         if ( nwtn_iter > max_newton_iter) then
            icnv = icnv_fail_round
            if(debug)then""")

# Picard: one line per iteration after the water-balance error
rep("""          wb_error = qtop*dtime - (sum( th_node_prev(:)*v_node(:) ) - sum( th_node(:)*v_node(:) ))*denh2o
""",
"""          wb_error = qtop*dtime - (sum( th_node_prev(:)*v_node(:) ) - sum( th_node(:)*v_node(:) ))*denh2o
          write(fates_log(),'(a,i4,es25.16e3,2i6,6es25.16e3)') 'PIC ', ft, cohort%dbh, icnv, picd_iter, &
               dtime, tm, wb_error, maxval(abs(residual(:))), maxval(abs(th_node(:) - th_prev(:))), cfl_max
""")
open(dst, 'w').write(t)
print('ok', dst)
