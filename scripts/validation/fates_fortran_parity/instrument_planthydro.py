"""instrument_planthydro.py — add a per-call / per-cohort dump to FATES `hydraulics_bc`.

Writes `hydro_pdump.txt` in the run directory, with three record kinds per call:
  SPRE   call nlevrhiz  h2oveg, h2o_liq_sisl(1:nlevsoil), shells after FillDrainRhizShells,
                        v_shell, l_aroot_layer, qflx_transp_pa(:)
  COH    call ifp pft nlevrhiz  n, dbh, q_top [kg/s/plant], g_sb_laweight, gscan_patch,
                        total_canopy_area, pre-solve th_ag(1:2)/th_troot/th_aroot, post-solve
                        ditto, psi_ag(1), btran, iterh1, iterh2, wb_err, dwat, sapflow,
                        rootuptake, l_aroot_layer, volumes, kmax_*, z_node_*
  SPOST  call nlevrhiz  h2oveg, shells after the update, qflx_soil2root_sisl, kmax_*_shell
The dump does not change the run (the relinked cesm.exe reproduces the reference h0
bit for bit, compare_h0_bitwise.jl). Build and run, with the SapwoodAreaGrass SourceMod:

  python3 instrument_planthydro.py $CTSM/src/fates/biogeophys/FatesPlantHydraulicsMod.F90 $OUT/FatesPlantHydraulicsMod.F90
  julia build_patched_cesm.jl $OUT/exe SourceMods_planthydro/FatesAllometryMod.F90 $OUT/FatesPlantHydraulicsMod.F90
  EXE=$OUT/exe/cesm.exe RUN=$OUT/run bash stage_planthydro_run.sh

then replay the port against it with fates_planthydro_replay.jl.
usage: python3 instrument_planthydro.py <src.F90> <dst.F90>"""
import sys
src, dst = sys.argv[1:3]
t = open(src).read()

def rep(anchor, new, count=1):
    global t
    assert t.count(anchor) == count, (anchor, t.count(anchor))
    t = t.replace(anchor, new)

# locals
rep("  real(r8) :: z_fr                ! Maximum fine root depth, used in disaggregation\n",
    "  real(r8) :: z_fr                ! Maximum fine root depth, used in disaggregation\n"
    "  integer, save :: pd_call = 0, pd_u = -1\n"
    "  real(r8) :: pd_pre(64)\n"
    "  integer :: pd_nl\n")
# open + counter
rep("  call RecruitWUptake(nsites,sites,bc_in,dtime,recruitflag)\n",
    "  pd_call = pd_call + 1\n"
    "  if (pd_call == 1) open(newunit=pd_u, file='hydro_pdump.txt', status='replace', action='write')\n"
    "  call RecruitWUptake(nsites,sites,bc_in,dtime,recruitflag)\n")
# site-level pre (shells after FillDrainRhizShells, host liquid)
rep("     prev_h2oveg    = csite_hydr%h2oveg\n",
    "     prev_h2oveg    = csite_hydr%h2oveg\n"
    "     pd_nl = csite_hydr%nlevrhiz\n"
    "     write(pd_u,'(a,2i8,9999es25.16e3)') 'SPRE ', pd_call, pd_nl, csite_hydr%h2oveg, &\n"
    "          bc_in(s)%h2o_liq_sisl(1:bc_in(s)%nlevsoil), csite_hydr%h2osoi_liqvol_shell(1:pd_nl,:), &\n"
    "          csite_hydr%v_shell(1:pd_nl,:), csite_hydr%l_aroot_layer(1:pd_nl), bc_in(s)%qflx_transp_pa(:)\n")
# cohort pre
rep("              ccohort_hydr%qtop = qflx_tran_veg_indiv*dtime\n",
    "              ccohort_hydr%qtop = qflx_tran_veg_indiv*dtime\n"
    "              pd_pre(1:2) = ccohort_hydr%th_ag(1:2); pd_pre(3) = ccohort_hydr%th_troot\n"
    "              pd_pre(4:3+pd_nl) = ccohort_hydr%th_aroot(1:pd_nl)\n")
# cohort post
rep("              ccohort_hydr%btran = wkf_plant(stomata_p_media,ft)%p%ftc_from_psi(ccohort_hydr%psi_ag(1))\n\n"
    "              ccohort => ccohort%shorter\n           enddo co_loop1",
    "              ccohort_hydr%btran = wkf_plant(stomata_p_media,ft)%p%ftc_from_psi(ccohort_hydr%psi_ag(1))\n"
    "              write(pd_u,'(a,4i8,9999es25.16e3)') 'COH ', pd_call, ifp, ft, pd_nl, ccohort%n, ccohort%dbh, &\n"
    "                   qflx_tran_veg_indiv, ccohort%g_sb_laweight, gscan_patch, cpatch%total_canopy_area, &\n"
    "                   pd_pre(1:3+pd_nl), ccohort_hydr%th_ag(1:2), ccohort_hydr%th_troot, &\n"
    "                   ccohort_hydr%th_aroot(1:pd_nl), ccohort_hydr%psi_ag(1), ccohort_hydr%btran, &\n"
    "                   ccohort_hydr%iterh1, ccohort_hydr%iterh2, wb_err_plant, dwat_plant, sapflow, &\n"
    "                   rootuptake(1:pd_nl), ccohort_hydr%l_aroot_layer(1:pd_nl), ccohort_hydr%v_ag(1:2), &\n"
    "                   ccohort_hydr%v_troot, ccohort_hydr%v_aroot_layer(1:pd_nl), &\n"
    "                   ccohort_hydr%kmax_petiole_to_leaf, ccohort_hydr%kmax_stem_upper(1), &\n"
    "                   ccohort_hydr%kmax_stem_lower(1), ccohort_hydr%kmax_troot_upper, &\n"
    "                   ccohort_hydr%kmax_troot_lower(1:pd_nl), ccohort_hydr%kmax_aroot_upper(1:pd_nl), &\n"
    "                   ccohort_hydr%kmax_aroot_lower(1:pd_nl), ccohort_hydr%kmax_aroot_radial_in(1:pd_nl), &\n"
    "                   ccohort_hydr%kmax_aroot_radial_out(1:pd_nl), ccohort_hydr%z_node_ag(1:2), &\n"
    "                   ccohort_hydr%z_node_troot\n"
    "              ccohort => ccohort%shorter\n           enddo co_loop1")
# site post (shells after update + soil2root)
rep("     csite_hydr%h2oveg_hydro_err = csite_hydr%h2oveg_hydro_err + csite_hydr%errh2o_hyd\n",
    "     csite_hydr%h2oveg_hydro_err = csite_hydr%h2oveg_hydro_err + csite_hydr%errh2o_hyd\n"
    "     write(pd_u,'(a,2i8,9999es25.16e3)') 'SPOST ', pd_call, pd_nl, csite_hydr%h2oveg, &\n"
    "          csite_hydr%h2osoi_liqvol_shell(1:pd_nl,:), bc_out(s)%qflx_soil2root_sisl(1:bc_in(s)%nlevsoil), &\n"
    "          csite_hydr%kmax_upper_shell(1:pd_nl,:), csite_hydr%kmax_lower_shell(1:pd_nl,:)\n"
    "     flush(pd_u)\n")
open(dst, 'w').write(t)
print('ok', dst)
