module pdumpMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Parity-dump instrumentation for the CLM.jl validation harness.
  !
  ! Reconstructs the per-boundary "pdump_<boundary>_n<nstep>.nc" snapshots that
  ! the CLM.jl :parity oracle (scripts/fortran_parity_common.jl
  ! compare_inst_to_dump / run_one_parity_step!) reads from the parity run dir.
  !
  ! The oracle only reads the 16 prognostic fields in its _parity_registry:
  !   T_SOISNO, H2OSOI_LIQ, H2OSOI_ICE  (column, levtot)
  !   T_GRND, WA, H2OSFC, ZWT, ZWT_PERCH, INT_SNOW, SNOW_DEPTH, frac_sno (column)
  !   T_VEG, elai, tlai, LIQCAN, SNOCAN  (pft)
  ! plus the timemgr nstep/ymd/tod metadata.  This module emits exactly that
  ! oracle-readable subset with the SAME dim names and the SAME levtot layout
  ! as the legacy full-restart pdump files (snow layers -nlevsno+1..0 mapped to
  ! levtot 1..nlevsno, then ground 1..nlevgrnd mapped to nlevsno+1..nlevtot;
  ! the full array is written verbatim, inactive snow slots carry their stored
  ! value, typically 0 -- matching the legacy restart-style write).
  !
  ! It is NOT a byte-for-byte replica of the 319-variable restart-style pdump;
  ! it is the minimal snapshot the harness consumes, and works for ANY config
  ! (SP, use_cn, use_hydrstress, ...) because all 16 fields exist in every
  ! configuration.
  !
  ! Modeled on the proven dump_cv_sidecar in SoilTemperatureMod.F90.
  !
  ! Call hooks (see clm_driver.F90 diff in the workstream report):
  !   - before_step           : top of the main clump physics loop
  !   - after_hydrologydrainage: immediately after the HydrologyDrainage call
  ! (those are the two boundaries the oracle currently uses; the other four
  !  boundary labels can be added with the same one-line call.)
  !-----------------------------------------------------------------------

  use shr_kind_mod    , only : r8 => shr_kind_r8
  use decompMod       , only : bounds_type
  use clm_varctl      , only : iulog
  use spmdMod         , only : masterproc

  implicit none
  private

  ! nstep window so a production run does not flood the run directory.
  ! Defaults bracket the legacy Bow parity target (n13461); override at will.
  integer, parameter :: pdump_nstep_lo = 13455
  integer, parameter :: pdump_nstep_hi = 13470

  public :: pdump_write

contains

  !-----------------------------------------------------------------------
  subroutine pdump_write(bounds, label)
    !
    ! !DESCRIPTION:
    ! Write the oracle-readable parity snapshot for the current timestep and
    ! the given boundary label to pdump_<label>_n<nstep>.nc in the run dir.
    !
    ! !USES:
    use netcdf
    use clm_varpar      , only : nlevsno, nlevgrnd
    use clm_time_manager, only : get_nstep, get_curr_date
    use clm_instMod     , only : temperature_inst, canopystate_inst, &
                                 soilhydrology_inst, water_inst
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    character(len=*) , intent(in) :: label   ! e.g. 'before_step', 'after_hydrologydrainage'
    !
    ! !LOCAL VARIABLES:
    real(r8), parameter :: fillval = 1.e+36_r8
    integer  :: nstep, yr, mon, day, tod, ymd
    integer  :: ncol, npft, nlevtot
    integer  :: c, p, j, lt, ci, pi
    integer  :: ncid, ier
    integer  :: d_col, d_pft, d_lt
    integer  :: v_tsoi, v_hliq, v_hice
    integer  :: v_tgrnd, v_wa, v_h2osfc, v_zwt, v_zwtp, v_intsnow, v_snowdp, v_fsno
    integer  :: v_tveg, v_elai, v_tlai, v_liqcan, v_snocan
    integer  :: v_nstep, v_ymd, v_tod
    character(len=256) :: fname
    real(r8), allocatable :: tsoi(:,:), hliq(:,:), hice(:,:)   ! (levtot, col)
    real(r8), allocatable :: c1(:)                              ! (col)
    real(r8), allocatable :: p1(:)                              ! (pft)
    !-----------------------------------------------------------------------

    nstep = get_nstep()
    if (nstep < pdump_nstep_lo .or. nstep > pdump_nstep_hi) return
    if (.not. masterproc) return

    call get_curr_date(yr, mon, day, tod)
    ymd = yr*10000 + mon*100 + day

    ncol    = bounds%endc - bounds%begc + 1
    npft    = bounds%endp - bounds%begp + 1
    nlevtot = nlevsno + nlevgrnd

    allocate(tsoi(nlevtot, ncol), hliq(nlevtot, ncol), hice(nlevtot, ncol))
    allocate(c1(ncol), p1(npft))

    ! ---- column level-resolved fields (levtot = j + nlevsno) -------------
    tsoi(:,:) = fillval; hliq(:,:) = fillval; hice(:,:) = fillval
    do c = bounds%begc, bounds%endc
       ci = c - bounds%begc + 1
       do j = -nlevsno+1, nlevgrnd
          lt = j + nlevsno
          tsoi(lt, ci) = temperature_inst%t_soisno_col(c, j)
          hliq(lt, ci) = water_inst%waterstatebulk_inst%h2osoi_liq_col(c, j)
          hice(lt, ci) = water_inst%waterstatebulk_inst%h2osoi_ice_col(c, j)
       end do
    end do

    ! ---- file -------------------------------------------------------------
    write(fname,'(a,a,a,i0,a)') 'pdump_', trim(label), '_n', nstep, '.nc'

    ier = nf90_create(trim(fname), NF90_CLOBBER, ncid)
    if (ier /= NF90_NOERR) then
       write(iulog,*) 'pdump_write: nf90_create failed: ', trim(nf90_strerror(ier))
       deallocate(tsoi, hliq, hice, c1, p1)
       return
    end if

    ier = nf90_def_dim(ncid, 'column', ncol,    d_col)
    ier = nf90_def_dim(ncid, 'pft',    npft,    d_pft)
    ier = nf90_def_dim(ncid, 'levtot', nlevtot, d_lt)

    ! column, levtot -> Fortran (levtot, col) [first dimid fastest]
    ier = nf90_def_var(ncid, 'T_SOISNO',   NF90_DOUBLE, (/d_lt, d_col/), v_tsoi)
    ier = nf90_def_var(ncid, 'H2OSOI_LIQ', NF90_DOUBLE, (/d_lt, d_col/), v_hliq)
    ier = nf90_def_var(ncid, 'H2OSOI_ICE', NF90_DOUBLE, (/d_lt, d_col/), v_hice)
    ier = put_fill(ncid, v_tsoi, fillval)
    ier = put_fill(ncid, v_hliq, fillval)
    ier = put_fill(ncid, v_hice, fillval)

    ier = nf90_def_var(ncid, 'T_GRND',     NF90_DOUBLE, (/d_col/), v_tgrnd)
    ier = nf90_def_var(ncid, 'WA',         NF90_DOUBLE, (/d_col/), v_wa)
    ier = nf90_def_var(ncid, 'H2OSFC',     NF90_DOUBLE, (/d_col/), v_h2osfc)
    ier = nf90_def_var(ncid, 'ZWT',        NF90_DOUBLE, (/d_col/), v_zwt)
    ier = nf90_def_var(ncid, 'ZWT_PERCH',  NF90_DOUBLE, (/d_col/), v_zwtp)
    ier = nf90_def_var(ncid, 'INT_SNOW',   NF90_DOUBLE, (/d_col/), v_intsnow)
    ier = nf90_def_var(ncid, 'SNOW_DEPTH', NF90_DOUBLE, (/d_col/), v_snowdp)
    ier = nf90_def_var(ncid, 'frac_sno',   NF90_DOUBLE, (/d_col/), v_fsno)

    ier = nf90_def_var(ncid, 'T_VEG',      NF90_DOUBLE, (/d_pft/), v_tveg)
    ier = nf90_def_var(ncid, 'elai',       NF90_DOUBLE, (/d_pft/), v_elai)
    ier = nf90_def_var(ncid, 'tlai',       NF90_DOUBLE, (/d_pft/), v_tlai)
    ier = nf90_def_var(ncid, 'LIQCAN',     NF90_DOUBLE, (/d_pft/), v_liqcan)
    ier = nf90_def_var(ncid, 'SNOCAN',     NF90_DOUBLE, (/d_pft/), v_snocan)

    ier = nf90_def_var(ncid, 'nstep',                NF90_INT, v_nstep)
    ier = nf90_def_var(ncid, 'timemgr_rst_curr_ymd', NF90_INT, v_ymd)
    ier = nf90_def_var(ncid, 'timemgr_rst_curr_tod', NF90_INT, v_tod)

    ier = nf90_enddef(ncid)

    ! ---- write ------------------------------------------------------------
    ier = nf90_put_var(ncid, v_tsoi, tsoi)
    ier = nf90_put_var(ncid, v_hliq, hliq)
    ier = nf90_put_var(ncid, v_hice, hice)

    call col1d(bounds, water_inst%waterstatebulk_inst%wa_col,              c1); ier = nf90_put_var(ncid, v_wa,      c1)
    call col1d(bounds, water_inst%waterstatebulk_inst%h2osfc_col,          c1); ier = nf90_put_var(ncid, v_h2osfc,  c1)
    call col1d(bounds, water_inst%waterstatebulk_inst%int_snow_col,        c1); ier = nf90_put_var(ncid, v_intsnow, c1)
    call col1d(bounds, temperature_inst%t_grnd_col,                        c1); ier = nf90_put_var(ncid, v_tgrnd,   c1)
    call col1d(bounds, soilhydrology_inst%zwt_col,                         c1); ier = nf90_put_var(ncid, v_zwt,     c1)
    call col1d(bounds, soilhydrology_inst%zwt_perched_col,                 c1); ier = nf90_put_var(ncid, v_zwtp,    c1)
    call col1d(bounds, water_inst%waterdiagnosticbulk_inst%snow_depth_col, c1); ier = nf90_put_var(ncid, v_snowdp,  c1)
    call col1d(bounds, water_inst%waterdiagnosticbulk_inst%frac_sno_col,   c1); ier = nf90_put_var(ncid, v_fsno,    c1)

    call pft1d(bounds, temperature_inst%t_veg_patch,                p1); ier = nf90_put_var(ncid, v_tveg,   p1)
    call pft1d(bounds, canopystate_inst%elai_patch,                 p1); ier = nf90_put_var(ncid, v_elai,   p1)
    call pft1d(bounds, canopystate_inst%tlai_patch,                 p1); ier = nf90_put_var(ncid, v_tlai,   p1)
    call pft1d(bounds, water_inst%waterstatebulk_inst%liqcan_patch, p1); ier = nf90_put_var(ncid, v_liqcan, p1)
    call pft1d(bounds, water_inst%waterstatebulk_inst%snocan_patch, p1); ier = nf90_put_var(ncid, v_snocan, p1)

    ier = nf90_put_var(ncid, v_nstep, nstep)
    ier = nf90_put_var(ncid, v_ymd,   ymd)
    ier = nf90_put_var(ncid, v_tod,   tod)

    ier = nf90_close(ncid)
    if (ier /= NF90_NOERR) then
       write(iulog,*) 'pdump_write: nf90_close failed: ', trim(nf90_strerror(ier))
    else
       write(iulog,*) 'pdump_write: wrote ', trim(fname), ' nstep=', nstep, &
            ' ymd=', ymd, ' tod=', tod
    end if

    deallocate(tsoi, hliq, hice, c1, p1)

  end subroutine pdump_write

  !-----------------------------------------------------------------------
  integer function put_fill(ncid, varid, fillval) result(ier)
    use netcdf
    integer , intent(in) :: ncid, varid
    real(r8), intent(in) :: fillval
    ier = nf90_put_att(ncid, varid, '_FillValue',    fillval)
    ier = nf90_put_att(ncid, varid, 'missing_value', fillval)
  end function put_fill

  subroutine col1d(bounds, src, dst)
    type(bounds_type), intent(in)  :: bounds
    real(r8)         , intent(in)  :: src(bounds%begc:)
    real(r8)         , intent(out) :: dst(:)
    integer :: c
    do c = bounds%begc, bounds%endc
       dst(c - bounds%begc + 1) = src(c)
    end do
  end subroutine col1d

  subroutine pft1d(bounds, src, dst)
    type(bounds_type), intent(in)  :: bounds
    real(r8)         , intent(in)  :: src(bounds%begp:)
    real(r8)         , intent(out) :: dst(:)
    integer :: p
    do p = bounds%begp, bounds%endp
       dst(p - bounds%begp + 1) = src(p)
    end do
  end subroutine pft1d

end module pdumpMod
