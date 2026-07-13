module FatesParityDumpMod

  ! ======================================================================================
  ! FATES parity instrumentation (Fortran ground-truth side) — SCHEMA v2, PHASE-TAGGED.
  !
  ! v1 (PR #213) emitted a single cold-start (nstep=0) snapshot. That validated FATES
  ! *initialization* only: the FATES parameter file is byte-identical between the port
  ! and the install, so the nstep=0 cohorts are parameter-identical BY CONSTRUCTION.
  ! It said nothing about the FATES *dynamics* — photosynthesis, allocation, growth,
  ! mortality, recruitment, disturbance, canopy promotion/demotion — which are the bulk
  ! of the port.
  !
  ! v2 emits a TIME SERIES of PHASE-TAGGED snapshots so a divergence can be attributed
  ! to the FATES phase that produced it. Every record carries `t=<nstep>` and
  ! `ph=<phase>` identity keys in addition to the s/p/c ordinals, and all records are
  ! APPENDED to one file.
  !
  ! PHASES (identical codes on the Julia side; see scripts/fates_fortran_parity.jl):
  !    0  coldstart   EDInitMod::init_cohorts (nstep=0)
  !    1  fast        every timestep, AFTER the photosynthesis/flux solve and the
  !                   running-mean + hifrq-history updates, BEFORE the daily dynamics.
  !                   (Called from clm_driver's use_fates block, exactly where the
  !                   `if (is_beg_curr_day())` dynamics branch begins.)
  !   10  dyn_in      ed_ecosystem_dynamics entry (the daily step's INPUT state)
  !   11  phenology   after phenology()
  !   12  distrates   after fire_model() + disturbance_rates()
  !   13  integrate   after ed_integrate_state_variables()  <- growth / allocation /
  !                   maintenance turnover / PRT / mortality derivative
  !   14  recruit     after the per-patch recruitment() loop
  !   15  cohortfuse  after sort_cohorts / terminate_cohorts / fuse_cohorts
  !   16  spawn       after spawn_patches()   (disturbance -> new patch)
  !   17  patchfuse   after fuse_patches + terminate_patches (= end of
  !                   ed_ecosystem_dynamics)
  !   18  updatesite  after ed_update_site() (canopy_spread, canopy_structure =
  !                   promotion/demotion, trim_canopy) — the end of the daily step
  !
  ! Records (one keyed record per line):
  !   META   t= ph= nsites= ncoh= npatch=
  !   SITE   t= ph= s= carbon= maxdbh= ncoh= npatch= pbot= tsoil1= h2ovol1= effporo1=
  !   PATCH  t= ph= s= p= patchno= area= age= ncoh= btran_ft= tcanarea= nocomp=
  !          tveg= tgcm= esat_tv= eair= oair= cair= rb= daylfac= solad1= solai1= fphoto=
  !   COHORT t= ph= s= p= c= pft= n= dbh= height= coage= canopy_layer= status= carea=
  !          treelai= treesai= canopy_trim= l2fr= nv= isnew=
  !          leafc= fnrtc= sapwc= storec= structc= reproc=
  !          gpp= npp= resp= gpp_acc= npp_acc= resp_acc= gpp_hold= npp_hold= resp_hold=
  !          vcmax= jmax= tpu= rdark= resp_m= resp_g=
  !          ddbhdt= dhdt= dndt= bmort= cmort= hmort= frmort= smort= asmort= dmort=
  !          seed_prod=
  !
  ! Traversal order is IDENTICAL to the Julia side:
  !   patches  site%oldest_patch -> %younger   (numbered p=1..; this is also the FATES
  !                                             `ifp` index that addresses bc_in arrays)
  !   cohorts  patch%tallest     -> %shorter   (numbered c=1..)
  !
  ! The PATCH record carries the HLM->FATES boundary conditions (t_veg, eair, cair, rb,
  ! daylength factor, incident solar) and the SITE record the soil BCs. That is the
  ! discriminator that makes divergence ATTRIBUTABLE: if bc_in matches to round-off but
  ! the FATES output does not, the divergence is INSIDE FATES; if bc_in already differs,
  ! it came from the host land model, not from the FATES port.
  !
  ! Reals are written with 17 significant digits (es-format) so the dump is
  ! bit-parity-capable (the Julia side uses %.17g).
  ! ======================================================================================

  use FatesConstantsMod,     only : r8 => fates_r8
  use EDTypesMod,            only : ed_site_type
  use FatesPatchMod,         only : fates_patch_type
  use FatesCohortMod,        only : fates_cohort_type
  use FatesInterfaceTypesMod, only : bc_in_type
  use PRTGenericMod,         only : leaf_organ, fnrt_organ, sapw_organ
  use PRTGenericMod,         only : store_organ, struct_organ, repro_organ
  use PRTGenericMod,         only : carbon12_element

  implicit none
  private

  ! ---- module state, set by the host instrumentation --------------------------------
  integer, public :: parity_nstep = 0        ! current HLM timestep (set in clm_driver)
  integer, public :: parity_site  = 1        ! current FATES site ordinal (set in dynamics_driv)
  logical, public :: parity_opened = .false. ! first write truncates, later ones append
  ! Armed only around the per-site ed_ecosystem_dynamics/ed_update_site calls in
  ! dynamics_driv, so the EDMainMod sub-phase hooks do NOT fire during the cold-start
  ! and restart-read calls to ed_update_site (which the Julia harness does not mirror).
  logical, public :: parity_daily = .false.

  character(len=*), parameter, public :: parity_path = 'fates_pdump_fortran.txt'

  public :: fates_parity_dump        ! dump a whole site array (coldstart / fast phases)
  public :: fates_parity_dump_site   ! dump ONE site (the ed_ecosystem_dynamics phases)

contains

  ! ------------------------------------------------------------------------------------
  function rstr(x) result(s)
    real(r8), intent(in) :: x
    character(len=32)    :: s
    write(s,'(es25.17e3)') x
    s = adjustl(s)
  end function rstr

  function istr(i) result(s)
    integer, intent(in) :: i
    character(len=16)   :: s
    write(s,'(i0)') i
    s = adjustl(s)
  end function istr

  ! ------------------------------------------------------------------------------------
  ! Open the dump file: truncate on the very first write of the run, append thereafter.
  ! ------------------------------------------------------------------------------------
  subroutine open_dump(iu, ios)
    integer, intent(out) :: iu, ios
    if (.not. parity_opened) then
       open(newunit=iu, file=parity_path, status='replace', action='write', iostat=ios)
       if (ios == 0) then
          write(iu,'(a)') '# FATES parity dump - schema v2 - generator=fortran'
          parity_opened = .true.
       end if
    else
       open(newunit=iu, file=parity_path, status='old', action='write', &
            position='append', iostat=ios)
    end if
  end subroutine open_dump

  ! ------------------------------------------------------------------------------------
  ! Dump every site in `sites(1:nsites)` at (nstep, phase).
  ! ------------------------------------------------------------------------------------
  subroutine fates_parity_dump(sites, bc_in, nsites, nstep, phase)

    type(ed_site_type), intent(inout) :: sites(:)
    type(bc_in_type),   intent(in)    :: bc_in(:)
    integer,            intent(in)    :: nsites, nstep, phase

    type(fates_patch_type),  pointer :: cpatch
    type(fates_cohort_type), pointer :: ccohort
    integer :: iu, ios, s, ncoh_tot, npatch_tot

    call open_dump(iu, ios)
    if (ios /= 0) return

    ncoh_tot = 0; npatch_tot = 0
    do s = 1, nsites
       cpatch => sites(s)%oldest_patch
       do while (associated(cpatch))
          npatch_tot = npatch_tot + 1
          ccohort => cpatch%tallest
          do while (associated(ccohort))
             ncoh_tot = ncoh_tot + 1
             ccohort => ccohort%shorter
          end do
          cpatch => cpatch%younger
       end do
    end do
    write(iu,'(a)') 'META t='//trim(istr(nstep))//' ph='//trim(istr(phase)) &
         //' nsites='//trim(istr(nsites))//' ncoh='//trim(istr(ncoh_tot)) &
         //' npatch='//trim(istr(npatch_tot))

    do s = 1, nsites
       call write_site(iu, sites(s), bc_in(s), s, nstep, phase)
    end do
    close(iu)

  end subroutine fates_parity_dump

  ! ------------------------------------------------------------------------------------
  ! Dump ONE site (used from inside ed_ecosystem_dynamics / ed_update_site, which are
  ! called per-site and do not see the site array). No META record: the comparator
  ! reconciles structure from the SITE record's ncoh/npatch.
  ! ------------------------------------------------------------------------------------
  subroutine fates_parity_dump_site(csite, bc_in, s, nstep, phase)

    type(ed_site_type), intent(inout) :: csite
    type(bc_in_type),   intent(in)    :: bc_in
    integer,            intent(in)    :: s, nstep, phase

    integer :: iu, ios

    if (.not. parity_daily) return
    call open_dump(iu, ios)
    if (ios /= 0) return
    call write_site(iu, csite, bc_in, s, nstep, phase)
    close(iu)

  end subroutine fates_parity_dump_site

  ! ------------------------------------------------------------------------------------
  subroutine write_site(iu, csite, bc_in, s, nstep, phase)

    integer,            intent(in)    :: iu, s, nstep, phase
    type(ed_site_type), intent(inout) :: csite
    type(bc_in_type),   intent(in)    :: bc_in

    type(fates_patch_type),  pointer :: cpatch
    type(fates_cohort_type), pointer :: ccohort
    integer  :: p, c, ncoh_p, ncoh_s, npatch_s
    real(r8) :: site_c, site_maxdbh, btr
    real(r8) :: leafc, fnrtc, sapwc, storec, structc, reproc
    real(r8) :: pbot, tsoil1, h2ovol1, effporo1
    character(len=48) :: hdr

    ! ---- site totals -------------------------------------------------------------
    site_c = 0.0_r8; site_maxdbh = 0.0_r8; ncoh_s = 0; npatch_s = 0
    cpatch => csite%oldest_patch
    do while (associated(cpatch))
       npatch_s = npatch_s + 1
       ccohort => cpatch%tallest
       do while (associated(ccohort))
          ncoh_s = ncoh_s + 1
          if (ccohort%dbh > site_maxdbh) site_maxdbh = ccohort%dbh
          site_c = site_c + pool_sum(ccohort)*ccohort%n
          ccohort => ccohort%shorter
       end do
       cpatch => cpatch%younger
    end do

    ! ---- soil boundary conditions (guard un-allocated / zero-length arrays) -------
    pbot = bc_in%forc_pbot
    tsoil1   = -9999.0_r8
    h2ovol1  = -9999.0_r8
    effporo1 = -9999.0_r8
    if (allocated(bc_in%t_soisno_sl))    then
       if (size(bc_in%t_soisno_sl)    >= 1) tsoil1   = bc_in%t_soisno_sl(1)
    end if
    if (allocated(bc_in%h2o_liqvol_sl))  then
       if (size(bc_in%h2o_liqvol_sl)  >= 1) h2ovol1  = bc_in%h2o_liqvol_sl(1)
    end if
    if (allocated(bc_in%eff_porosity_sl)) then
       if (size(bc_in%eff_porosity_sl) >= 1) effporo1 = bc_in%eff_porosity_sl(1)
    end if

    hdr = 't='//trim(istr(nstep))//' ph='//trim(istr(phase))//' s='//trim(istr(s))

    write(iu,'(a)') 'SITE '//trim(hdr) &
         //' carbon='//trim(rstr(site_c)) &
         //' maxdbh='//trim(rstr(site_maxdbh)) &
         //' ncoh='//trim(istr(ncoh_s)) &
         //' npatch='//trim(istr(npatch_s)) &
         //' pbot='//trim(rstr(pbot)) &
         //' tsoil1='//trim(rstr(tsoil1)) &
         //' h2ovol1='//trim(rstr(h2ovol1)) &
         //' effporo1='//trim(rstr(effporo1))

    ! ---- patches + cohorts -------------------------------------------------------
    p = 0
    cpatch => csite%oldest_patch
    do while (associated(cpatch))
       p = p + 1
       btr = maxval(cpatch%btran_ft)
       ncoh_p = 0
       ccohort => cpatch%tallest
       do while (associated(ccohort))
          ncoh_p = ncoh_p + 1
          ccohort => ccohort%shorter
       end do

       write(iu,'(a)') 'PATCH '//trim(hdr)//' p='//trim(istr(p)) &
            //' patchno='//trim(istr(cpatch%patchno)) &
            //' area='//trim(rstr(cpatch%area)) &
            //' age='//trim(rstr(cpatch%age)) &
            //' ncoh='//trim(istr(ncoh_p)) &
            //' btran_ft='//trim(rstr(btr)) &
            //' tcanarea='//trim(rstr(cpatch%total_canopy_area)) &
            //' nocomp='//trim(istr(cpatch%nocomp_pft_label)) &
            //' tveg='//trim(rstr(bcp(bc_in%t_veg_pa, p))) &
            //' tgcm='//trim(rstr(bcp(bc_in%tgcm_pa, p))) &
            //' esat_tv='//trim(rstr(bcp(bc_in%esat_tv_pa, p))) &
            //' eair='//trim(rstr(bcp(bc_in%eair_pa, p))) &
            //' oair='//trim(rstr(bcp(bc_in%oair_pa, p))) &
            //' cair='//trim(rstr(bcp(bc_in%cair_pa, p))) &
            //' rb='//trim(rstr(bcp(bc_in%rb_pa, p))) &
            //' daylfac='//trim(rstr(bcp(bc_in%dayl_factor_pa, p))) &
            //' solad1='//trim(rstr(bcm(bc_in%solad_parb, p))) &
            //' solai1='//trim(rstr(bcm(bc_in%solai_parb, p)))

       c = 0
       ccohort => cpatch%tallest
       do while (associated(ccohort))
          c = c + 1
          leafc   = ccohort%prt%GetState(leaf_organ,   carbon12_element)
          fnrtc   = ccohort%prt%GetState(fnrt_organ,   carbon12_element)
          sapwc   = ccohort%prt%GetState(sapw_organ,   carbon12_element)
          storec  = ccohort%prt%GetState(store_organ,  carbon12_element)
          structc = ccohort%prt%GetState(struct_organ, carbon12_element)
          reproc  = ccohort%prt%GetState(repro_organ,  carbon12_element)

          write(iu,'(a)') 'COHORT '//trim(hdr)//' p='//trim(istr(p)) &
               //' c='//trim(istr(c)) &
               //' pft='//trim(istr(ccohort%pft)) &
               //' n='//trim(rstr(ccohort%n)) &
               //' dbh='//trim(rstr(ccohort%dbh)) &
               //' height='//trim(rstr(ccohort%height)) &
               //' coage='//trim(rstr(ccohort%coage)) &
               //' canopy_layer='//trim(istr(ccohort%canopy_layer)) &
               //' status='//trim(rstr(real(ccohort%status_coh,r8))) &
               //' carea='//trim(rstr(ccohort%c_area)) &
               //' treelai='//trim(rstr(ccohort%treelai)) &
               //' treesai='//trim(rstr(ccohort%treesai)) &
               //' canopy_trim='//trim(rstr(ccohort%canopy_trim)) &
               //' l2fr='//trim(rstr(ccohort%l2fr)) &
               //' nv='//trim(istr(ccohort%nv)) &
               //' isnew='//trim(istr(merge(1,0,ccohort%isnew))) &
               //' leafc='//trim(rstr(leafc)) &
               //' fnrtc='//trim(rstr(fnrtc)) &
               //' sapwc='//trim(rstr(sapwc)) &
               //' storec='//trim(rstr(storec)) &
               //' structc='//trim(rstr(structc)) &
               //' reproc='//trim(rstr(reproc)) &
               //' gpp='//trim(rstr(ccohort%gpp_tstep)) &
               //' npp='//trim(rstr(ccohort%npp_tstep)) &
               //' resp='//trim(rstr(ccohort%resp_tstep)) &
               //' gpp_acc='//trim(rstr(ccohort%gpp_acc)) &
               //' npp_acc='//trim(rstr(ccohort%npp_acc)) &
               //' resp_acc='//trim(rstr(ccohort%resp_acc)) &
               //' gpp_hold='//trim(rstr(ccohort%gpp_acc_hold)) &
               //' npp_hold='//trim(rstr(ccohort%npp_acc_hold)) &
               //' resp_hold='//trim(rstr(ccohort%resp_acc_hold)) &
               //' vcmax='//trim(rstr(ccohort%vcmax25top)) &
               //' jmax='//trim(rstr(ccohort%jmax25top)) &
               //' tpu='//trim(rstr(ccohort%tpu25top)) &
               //' rdark='//trim(rstr(ccohort%rdark)) &
               //' resp_m='//trim(rstr(ccohort%resp_m)) &
               //' resp_g='//trim(rstr(ccohort%resp_g_tstep)) &
               //' ddbhdt='//trim(rstr(ccohort%ddbhdt)) &
               //' dhdt='//trim(rstr(ccohort%dhdt)) &
               //' dndt='//trim(rstr(ccohort%dndt)) &
               //' bmort='//trim(rstr(ccohort%bmort)) &
               //' cmort='//trim(rstr(ccohort%cmort)) &
               //' hmort='//trim(rstr(ccohort%hmort)) &
               //' frmort='//trim(rstr(ccohort%frmort)) &
               //' smort='//trim(rstr(ccohort%smort)) &
               //' asmort='//trim(rstr(ccohort%asmort)) &
               //' dmort='//trim(rstr(ccohort%dmort)) &
               //' seed_prod='//trim(rstr(ccohort%seed_prod))
          ccohort => ccohort%shorter
       end do
       cpatch => cpatch%younger
    end do

  end subroutine write_site

  ! ------------------------------------------------------------------------------------
  ! Safe accessors for the per-patch bc_in arrays: FATES sizes them to maxpatch, but a
  ! phase may run before they are packed (or the array may be unallocated in a mode that
  ! does not use it). Out-of-range / unallocated -> sentinel.
  ! ------------------------------------------------------------------------------------
  function bcp(a, i) result(v)
    real(r8), allocatable, intent(in) :: a(:)
    integer,               intent(in) :: i
    real(r8) :: v
    v = -9999.0_r8
    if (allocated(a)) then
       if (i >= lbound(a,1) .and. i <= ubound(a,1)) v = a(i)
    end if
  end function bcp

  function bcm(a, i) result(v)
    real(r8), allocatable, intent(in) :: a(:,:)
    integer,               intent(in) :: i
    real(r8) :: v
    v = -9999.0_r8
    if (allocated(a)) then
       if (i >= lbound(a,1) .and. i <= ubound(a,1) .and. size(a,2) >= 1) v = a(i,1)
    end if
  end function bcm

  ! ------------------------------------------------------------------------------------
  function pool_sum(ccohort) result(tot)
    type(fates_cohort_type), intent(inout) :: ccohort
    real(r8) :: tot
    tot = ccohort%prt%GetState(leaf_organ,   carbon12_element) &
        + ccohort%prt%GetState(fnrt_organ,   carbon12_element) &
        + ccohort%prt%GetState(sapw_organ,   carbon12_element) &
        + ccohort%prt%GetState(store_organ,  carbon12_element) &
        + ccohort%prt%GetState(struct_organ, carbon12_element)
  end function pool_sum

end module FatesParityDumpMod
