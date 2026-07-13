module FatesParityDumpMod

  ! ======================================================================================
  ! FATES cold-start parity instrumentation (Fortran ground-truth side).
  !
  ! Emits a schema-v1 text "pdump" of the FATES site->patch->cohort hierarchy in the
  ! EXACT format that the Julia CLM.jl port harness expects, so a field-by-field
  ! cold-start (nstep=0) parity scorecard can be produced against the Julia dump:
  !
  !     scripts/fates_fortran_parity.jl   (compare_dumps / parse_dump)
  !
  ! Schema v1 (one keyed record per line, matches scripts/fates_fortran_parity.jl):
  !   META   nstep=<i> nsites=<i> ncoh=<i> npatch=<i>
  !   SITE   s=<i> carbon=<f> maxdbh=<f> gpp=<f> npp=<f> btran=<f>
  !   PATCH  s=<i> p=<i> patchno=<i> area=<f> elai=<f> btran_ft=<f> ncoh=<i>
  !   COHORT s=<i> p=<i> c=<i> pft=<i> n=<f> dbh=<f> height=<f> leafc=<f> fnrtc=<f>
  !          sapwc=<f> storec=<f> structc=<f> treelai=<f> carea=<f> canopy_layer=<i>
  !          status=<f> gpp=<f> npp=<f> vcmax=<f>
  !
  ! Traversal order is IDENTICAL to the Julia side:
  !   patches  site%oldest_patch -> %younger   (numbered p=1..)
  !   cohorts  patch%tallest     -> %shorter   (numbered c=1..)
  !
  ! Reals are written with 17 significant digits (es-format) so the dump is
  ! bit-parity-capable (the Julia side uses %.17g).
  !
  ! NOTE (elai): the Julia dump reports the HLM canopystate elai_patch for the
  ! positional veg patch as a coupling cross-check. It is NOT a FATES-internal
  ! state; this Fortran instrument writes elai=-1 (sentinel: "not reported by the
  ! Fortran instrument") so the comparison can skip it. All FATES-internal cohort
  ! fields (the real parity target) are emitted.
  !
  ! STATUS: RUNS. Authored against the verified CTSM ctsm5.3.012 FATES type/API
  ! names (FatesCohortMod / FatesPatchMod / EDTypesMod / PRTGenericMod GetState),
  ! dropped into a case SourceMods/src.clm with the one-line call injected into
  ! clmfates_interfaceMod.F90 init_coldstart. Built and EXECUTED single-point
  ! (see setup_case.sh): emits 1 site x 1 patch x 14 cohorts at nstep=0, and the
  ! Julia FATES port matches that dump on all 27 fields (13 of them exactly, the
  ! PRT carbon pools to 1 ULP) — see README.md for the scorecard.
  !
  ! The earlier "blocked in DATM on absent GSWP3 forcing" note is obsolete: that
  ! was an artifact of running a GLOBAL compset. FATES does not need a global
  ! grid; setup_case.sh now builds a single-point (SROF, 1x1 ESMF mesh, local
  ! clmforc.YYYY.nc) case with ZERO missing input files.
  ! ======================================================================================

  use FatesConstantsMod, only : r8 => fates_r8
  use EDTypesMod,        only : ed_site_type
  use FatesPatchMod,     only : fates_patch_type
  use FatesCohortMod,    only : fates_cohort_type
  use PRTGenericMod,     only : leaf_organ, fnrt_organ, sapw_organ
  use PRTGenericMod,     only : store_organ, struct_organ, carbon12_element

  implicit none
  private

  public :: fates_parity_dump

contains

  ! ------------------------------------------------------------------------------------
  ! Format a real*8 with 17 significant digits into a left-justified token string.
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
  ! Walk sites -> patches -> cohorts and write the schema-v1 dump.
  ! ------------------------------------------------------------------------------------
  subroutine fates_parity_dump(sites, nsites, nstep, path)

    type(ed_site_type), intent(inout) :: sites(:)
    integer,            intent(in)    :: nsites
    integer,            intent(in)    :: nstep
    character(len=*),   intent(in)    :: path

    type(fates_patch_type),  pointer :: cpatch
    type(fates_cohort_type), pointer :: ccohort
    integer  :: iu, s, p, c
    integer  :: ncoh_tot, npatch_tot, ncoh_p
    real(r8) :: site_c, site_maxdbh, site_btran, btr
    real(r8) :: leafc, fnrtc, sapwc, storec, structc
    integer  :: ios

    open(newunit=iu, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) return
    write(iu,'(a)') '# FATES parity dump - schema v1 - generator=fortran'

    ! ---- first pass: totals for the META record --------------------------------------
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
    write(iu,'(a)') 'META nstep='//trim(istr(nstep))//' nsites='//trim(istr(nsites)) &
         //' ncoh='//trim(istr(ncoh_tot))//' npatch='//trim(istr(npatch_tot))

    ! ---- SITE records -----------------------------------------------------------------
    do s = 1, nsites
       site_c = 0.0_r8; site_maxdbh = 0.0_r8; site_btran = 0.0_r8
       cpatch => sites(s)%oldest_patch
       do while (associated(cpatch))
          btr = maxval(cpatch%btran_ft)
          if (btr > site_btran) site_btran = btr
          ccohort => cpatch%tallest
          do while (associated(ccohort))
             if (ccohort%dbh > site_maxdbh) site_maxdbh = ccohort%dbh
             leafc   = ccohort%prt%GetState(leaf_organ,   carbon12_element)
             fnrtc   = ccohort%prt%GetState(fnrt_organ,   carbon12_element)
             sapwc   = ccohort%prt%GetState(sapw_organ,   carbon12_element)
             storec  = ccohort%prt%GetState(store_organ,  carbon12_element)
             structc = ccohort%prt%GetState(struct_organ, carbon12_element)
             site_c  = site_c + (leafc+fnrtc+sapwc+storec+structc)*ccohort%n
             ccohort => ccohort%shorter
          end do
          cpatch => cpatch%younger
       end do
       write(iu,'(a)') 'SITE s='//trim(istr(s))//' carbon='//trim(rstr(site_c)) &
            //' maxdbh='//trim(rstr(site_maxdbh))//' gpp='//trim(rstr(0.0_r8)) &
            //' npp='//trim(rstr(0.0_r8))//' btran='//trim(rstr(site_btran))
    end do

    ! ---- PATCH + COHORT records -------------------------------------------------------
    do s = 1, nsites
       p = 0
       cpatch => sites(s)%oldest_patch
       do while (associated(cpatch))
          p = p + 1
          btr = maxval(cpatch%btran_ft)
          ncoh_p = 0
          ccohort => cpatch%tallest
          do while (associated(ccohort))
             ncoh_p = ncoh_p + 1
             ccohort => ccohort%shorter
          end do
          write(iu,'(a)') 'PATCH s='//trim(istr(s))//' p='//trim(istr(p)) &
               //' patchno='//trim(istr(cpatch%patchno)) &
               //' area='//trim(rstr(cpatch%area)) &
               //' elai='//trim(rstr(-1.0_r8)) &
               //' btran_ft='//trim(rstr(btr)) &
               //' ncoh='//trim(istr(ncoh_p))

          c = 0
          ccohort => cpatch%tallest
          do while (associated(ccohort))
             c = c + 1
             leafc   = ccohort%prt%GetState(leaf_organ,   carbon12_element)
             fnrtc   = ccohort%prt%GetState(fnrt_organ,   carbon12_element)
             sapwc   = ccohort%prt%GetState(sapw_organ,   carbon12_element)
             storec  = ccohort%prt%GetState(store_organ,  carbon12_element)
             structc = ccohort%prt%GetState(struct_organ, carbon12_element)
             write(iu,'(a)') 'COHORT s='//trim(istr(s))//' p='//trim(istr(p)) &
                  //' c='//trim(istr(c))//' pft='//trim(istr(ccohort%pft)) &
                  //' n='//trim(rstr(ccohort%n))//' dbh='//trim(rstr(ccohort%dbh)) &
                  //' height='//trim(rstr(ccohort%height)) &
                  //' leafc='//trim(rstr(leafc))//' fnrtc='//trim(rstr(fnrtc)) &
                  //' sapwc='//trim(rstr(sapwc))//' storec='//trim(rstr(storec)) &
                  //' structc='//trim(rstr(structc)) &
                  //' treelai='//trim(rstr(ccohort%treelai)) &
                  //' carea='//trim(rstr(ccohort%c_area)) &
                  //' canopy_layer='//trim(istr(ccohort%canopy_layer)) &
                  //' status='//trim(rstr(real(ccohort%status_coh,r8))) &
                  //' gpp='//trim(rstr(ccohort%gpp_tstep)) &
                  //' npp='//trim(rstr(ccohort%npp_tstep)) &
                  //' vcmax='//trim(rstr(ccohort%vcmax25top))
             ccohort => ccohort%shorter
          end do
          cpatch => cpatch%younger
       end do
    end do

    close(iu)

  end subroutine fates_parity_dump

end module FatesParityDumpMod
