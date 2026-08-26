module btranOracleMod

  ! Read-only shared-state oracle for the Bow BTRAN/LUNA discrepancy.
  ! Add `use btranOracleMod, only : btran_oracle_write` to clm_driver and call
  ! `btran_oracle_write(bounds, label)` beside each pdump_write hook.  The
  ! compile-time window is intentionally small; change only these two bounds
  ! when moving from the six-step neutrality smoke to the target trajectory.

  use shr_kind_mod, only : r8 => shr_kind_r8
  use decompMod   , only : bounds_type
  use clm_varctl  , only : iulog
  use spmdMod     , only : masterproc

  implicit none
  private

  integer, parameter :: oracle_nstep_lo = 2832
  integer, parameter :: oracle_nstep_hi = 3024

  public :: btran_oracle_write

contains

  subroutine btran_oracle_write(bounds, label)
    use netcdf
    use clm_varpar      , only : nlevgrnd, nlevsoi, nlevcan, nvegwcs
    use clm_time_manager, only : get_nstep, get_curr_date
    use clm_instMod     , only : energyflux_inst, canopystate_inst, soilstate_inst, &
                                 atm2lnd_inst, solarabs_inst, temperature_inst, &
                                 photosyns_inst, frictionvel_inst, water_inst, surfalb_inst

    type(bounds_type), intent(in) :: bounds
    character(len=*), intent(in) :: label

    integer :: nstep, yr, mon, day, tod, ymd, ncol, npft
    integer :: ncid, ier, dcol, dpft, dground, dsoil, dcan, dvegwcs
    integer :: v_nstep, v_ymd, v_tod, v_btran, v_vegwp, v_smp, v_hk, v_rootcond
    integer :: v_pbot240, v_pco2240, v_po2240, v_elai240
    integer :: v_par240d, v_par240x, v_ta10, v_tv10d, v_tv10n, v_rh10, v_rb10
    integer :: v_vcmx25, v_jmx25, v_pnlc, v_enzs
    integer :: v_kroot, v_soilcond, v_rootfr, v_bsw, v_sucsat
    integer :: v_laisun, v_laisha, v_elai, v_esai, v_tsai, v_fdry
    integer :: v_rho, v_pbot, v_thm, v_qtran, v_bsun, v_bsha
    integer :: v_nrad, v_tlaiz, v_fsunz, v_fabdsun, v_fabdsha, v_fabisun, v_fabisha
    integer :: v_vcintsun, v_vcintsha
    character(len=256) :: fname
    real(r8), allocatable :: p1(:), pveg(:,:), pcana(:,:), cground(:,:), psoil(:,:)
    integer, allocatable :: pi1(:)

    nstep = get_nstep()
    if (nstep < oracle_nstep_lo .or. nstep > oracle_nstep_hi) return
    if (.not. masterproc) return

    call get_curr_date(yr, mon, day, tod)
    ymd = yr*10000 + mon*100 + day
    ncol = bounds%endc - bounds%begc + 1
    npft = bounds%endp - bounds%begp + 1
    allocate(p1(npft), pveg(nvegwcs,npft), pcana(nlevcan,npft))
    allocate(cground(nlevgrnd,ncol), psoil(nlevsoi,npft))
    allocate(pi1(npft))

    write(fname,'(a,a,a,i0,a)') 'btran_oracle_', trim(label), '_n', nstep, '.nc'
    ier = nf90_create(trim(fname), NF90_CLOBBER, ncid)
    if (ier /= NF90_NOERR) then
       write(iulog,*) 'btran_oracle_write: create failed: ', trim(nf90_strerror(ier))
       call cleanup()
       return
    end if

    ier = nf90_def_dim(ncid, 'column', ncol, dcol)
    ier = nf90_def_dim(ncid, 'pft', npft, dpft)
    ier = nf90_def_dim(ncid, 'levgrnd', nlevgrnd, dground)
    ier = nf90_def_dim(ncid, 'levsoi', nlevsoi, dsoil)
    ier = nf90_def_dim(ncid, 'levcan', nlevcan, dcan)
    ier = nf90_def_dim(ncid, 'nvegwcs', nvegwcs, dvegwcs)

    ier = nf90_def_var(ncid, 'nstep', NF90_INT, v_nstep)
    ier = nf90_def_var(ncid, 'timemgr_rst_curr_ymd', NF90_INT, v_ymd)
    ier = nf90_def_var(ncid, 'timemgr_rst_curr_tod', NF90_INT, v_tod)
    ier = nf90_def_var(ncid, 'BTRAN', NF90_DOUBLE, (/dpft/), v_btran)
    ier = nf90_def_var(ncid, 'VEGWP', NF90_DOUBLE, (/dvegwcs,dpft/), v_vegwp)
    ier = nf90_def_var(ncid, 'SMP', NF90_DOUBLE, (/dground,dcol/), v_smp)
    ier = nf90_def_var(ncid, 'HK', NF90_DOUBLE, (/dground,dcol/), v_hk)
    ier = nf90_def_var(ncid, 'ROOT_CONDUCTANCE', NF90_DOUBLE, (/dsoil,dpft/), v_rootcond)
    ier = nf90_def_var(ncid, 'PBOT240', NF90_DOUBLE, (/dpft/), v_pbot240)
    ier = nf90_def_var(ncid, 'PCO2_240', NF90_DOUBLE, (/dpft/), v_pco2240)
    ier = nf90_def_var(ncid, 'PO2_240', NF90_DOUBLE, (/dpft/), v_po2240)
    ier = nf90_def_var(ncid, 'ELAI240', NF90_DOUBLE, (/dpft/), v_elai240)
    ier = nf90_def_var(ncid, 'PAR240D_Z', NF90_DOUBLE, (/dcan,dpft/), v_par240d)
    ier = nf90_def_var(ncid, 'PAR240X_Z', NF90_DOUBLE, (/dcan,dpft/), v_par240x)
    ier = nf90_def_var(ncid, 'T_A10', NF90_DOUBLE, (/dpft/), v_ta10)
    ier = nf90_def_var(ncid, 'T_VEG10_DAY', NF90_DOUBLE, (/dpft/), v_tv10d)
    ier = nf90_def_var(ncid, 'T_VEG10_NIGHT', NF90_DOUBLE, (/dpft/), v_tv10n)
    ier = nf90_def_var(ncid, 'RH10_AF', NF90_DOUBLE, (/dpft/), v_rh10)
    ier = nf90_def_var(ncid, 'RB10', NF90_DOUBLE, (/dpft/), v_rb10)
    ier = nf90_def_var(ncid, 'VCMX25_Z', NF90_DOUBLE, (/dcan,dpft/), v_vcmx25)
    ier = nf90_def_var(ncid, 'JMX25_Z', NF90_DOUBLE, (/dcan,dpft/), v_jmx25)
    ier = nf90_def_var(ncid, 'PNLC_Z', NF90_DOUBLE, (/dcan,dpft/), v_pnlc)
    ier = nf90_def_var(ncid, 'ENZS_Z', NF90_DOUBLE, (/dcan,dpft/), v_enzs)
    ier = nf90_def_var(ncid, 'K_SOIL_ROOT', NF90_DOUBLE, (/dsoil,dpft/), v_kroot)
    ier = nf90_def_var(ncid, 'SOIL_CONDUCTANCE', NF90_DOUBLE, (/dsoil,dpft/), v_soilcond)
    ier = nf90_def_var(ncid, 'ROOTFR', NF90_DOUBLE, (/dground,dpft/), v_rootfr)
    ier = nf90_def_var(ncid, 'BSW', NF90_DOUBLE, (/dground,dcol/), v_bsw)
    ier = nf90_def_var(ncid, 'SUCSAT', NF90_DOUBLE, (/dground,dcol/), v_sucsat)
    ier = nf90_def_var(ncid, 'LAISUN', NF90_DOUBLE, (/dpft/), v_laisun)
    ier = nf90_def_var(ncid, 'LAISHA', NF90_DOUBLE, (/dpft/), v_laisha)
    ier = nf90_def_var(ncid, 'ELAI', NF90_DOUBLE, (/dpft/), v_elai)
    ier = nf90_def_var(ncid, 'ESAI', NF90_DOUBLE, (/dpft/), v_esai)
    ier = nf90_def_var(ncid, 'TSAI', NF90_DOUBLE, (/dpft/), v_tsai)
    ier = nf90_def_var(ncid, 'FDRY', NF90_DOUBLE, (/dpft/), v_fdry)
    ier = nf90_def_var(ncid, 'FORC_RHO', NF90_DOUBLE, (/dcol/), v_rho)
    ier = nf90_def_var(ncid, 'FORC_PBOT', NF90_DOUBLE, (/dcol/), v_pbot)
    ier = nf90_def_var(ncid, 'THM', NF90_DOUBLE, (/dpft/), v_thm)
    ier = nf90_def_var(ncid, 'QFLX_TRAN_VEG', NF90_DOUBLE, (/dpft/), v_qtran)
    ier = nf90_def_var(ncid, 'BSUN', NF90_DOUBLE, (/dpft/), v_bsun)
    ier = nf90_def_var(ncid, 'BSHA', NF90_DOUBLE, (/dpft/), v_bsha)
    ier = nf90_def_var(ncid, 'NRAD', NF90_INT, (/dpft/), v_nrad)
    ier = nf90_def_var(ncid, 'TLAI_Z', NF90_DOUBLE, (/dcan,dpft/), v_tlaiz)
    ier = nf90_def_var(ncid, 'FSUN_Z', NF90_DOUBLE, (/dcan,dpft/), v_fsunz)
    ier = nf90_def_var(ncid, 'FABD_SUN_Z', NF90_DOUBLE, (/dcan,dpft/), v_fabdsun)
    ier = nf90_def_var(ncid, 'FABD_SHA_Z', NF90_DOUBLE, (/dcan,dpft/), v_fabdsha)
    ier = nf90_def_var(ncid, 'FABI_SUN_Z', NF90_DOUBLE, (/dcan,dpft/), v_fabisun)
    ier = nf90_def_var(ncid, 'FABI_SHA_Z', NF90_DOUBLE, (/dcan,dpft/), v_fabisha)
    ier = nf90_def_var(ncid, 'VCMAXCINTSUN', NF90_DOUBLE, (/dpft/), v_vcintsun)
    ier = nf90_def_var(ncid, 'VCMAXCINTSHA', NF90_DOUBLE, (/dpft/), v_vcintsha)
    ier = nf90_enddef(ncid)

    ier = nf90_put_var(ncid, v_nstep, nstep)
    ier = nf90_put_var(ncid, v_ymd, ymd)
    ier = nf90_put_var(ncid, v_tod, tod)
    call pft1d(bounds, energyflux_inst%btran_patch, p1); ier = nf90_put_var(ncid, v_btran, p1)
    call pft2d(bounds, nvegwcs, canopystate_inst%vegwp_patch, pveg); ier = nf90_put_var(ncid, v_vegwp, pveg)
    call col2d(bounds, nlevgrnd, soilstate_inst%smp_l_col, cground); ier = nf90_put_var(ncid, v_smp, cground)
    call col2d(bounds, nlevgrnd, soilstate_inst%hk_l_col, cground); ier = nf90_put_var(ncid, v_hk, cground)
    call pft2d(bounds, nlevsoi, soilstate_inst%root_conductance_patch, psoil); ier = nf90_put_var(ncid, v_rootcond, psoil)
    call pft1d(bounds, atm2lnd_inst%forc_pbot240_downscaled_patch, p1); ier = nf90_put_var(ncid, v_pbot240, p1)
    call pft1d(bounds, atm2lnd_inst%forc_pco2_240_patch, p1); ier = nf90_put_var(ncid, v_pco2240, p1)
    call pft1d(bounds, atm2lnd_inst%forc_po2_240_patch, p1); ier = nf90_put_var(ncid, v_po2240, p1)
    call pft1d(bounds, canopystate_inst%elai240_patch, p1); ier = nf90_put_var(ncid, v_elai240, p1)
    call pft2d(bounds, nlevcan, solarabs_inst%par240d_z_patch, pcana); ier = nf90_put_var(ncid, v_par240d, pcana)
    call pft2d(bounds, nlevcan, solarabs_inst%par240x_z_patch, pcana); ier = nf90_put_var(ncid, v_par240x, pcana)
    call pft1d(bounds, temperature_inst%t_a10_patch, p1); ier = nf90_put_var(ncid, v_ta10, p1)
    call pft1d(bounds, temperature_inst%t_veg10_day_patch, p1); ier = nf90_put_var(ncid, v_tv10d, p1)
    call pft1d(bounds, temperature_inst%t_veg10_night_patch, p1); ier = nf90_put_var(ncid, v_tv10n, p1)
    call pft1d(bounds, water_inst%waterdiagnosticbulk_inst%rh10_af_patch, p1); ier = nf90_put_var(ncid, v_rh10, p1)
    call pft1d(bounds, frictionvel_inst%rb10_patch, p1); ier = nf90_put_var(ncid, v_rb10, p1)
    call pft2d(bounds, nlevcan, photosyns_inst%vcmx25_z_patch, pcana); ier = nf90_put_var(ncid, v_vcmx25, pcana)
    call pft2d(bounds, nlevcan, photosyns_inst%jmx25_z_patch, pcana); ier = nf90_put_var(ncid, v_jmx25, pcana)
    call pft2d(bounds, nlevcan, photosyns_inst%pnlc_z_patch, pcana); ier = nf90_put_var(ncid, v_pnlc, pcana)
    call pft2d(bounds, nlevcan, photosyns_inst%enzs_z_patch, pcana); ier = nf90_put_var(ncid, v_enzs, pcana)
    call pft2d(bounds, nlevsoi, soilstate_inst%k_soil_root_patch, psoil); ier = nf90_put_var(ncid, v_kroot, psoil)
    call pft2d(bounds, nlevsoi, soilstate_inst%soil_conductance_patch, psoil); ier = nf90_put_var(ncid, v_soilcond, psoil)
    deallocate(psoil); allocate(psoil(nlevgrnd,npft))
    call pft2d(bounds, nlevgrnd, soilstate_inst%rootfr_patch, psoil); ier = nf90_put_var(ncid, v_rootfr, psoil)
    call col2d(bounds, nlevgrnd, soilstate_inst%bsw_col, cground); ier = nf90_put_var(ncid, v_bsw, cground)
    call col2d(bounds, nlevgrnd, soilstate_inst%sucsat_col, cground); ier = nf90_put_var(ncid, v_sucsat, cground)
    call pft1d(bounds, canopystate_inst%laisun_patch, p1); ier = nf90_put_var(ncid, v_laisun, p1)
    call pft1d(bounds, canopystate_inst%laisha_patch, p1); ier = nf90_put_var(ncid, v_laisha, p1)
    call pft1d(bounds, canopystate_inst%elai_patch, p1); ier = nf90_put_var(ncid, v_elai, p1)
    call pft1d(bounds, canopystate_inst%esai_patch, p1); ier = nf90_put_var(ncid, v_esai, p1)
    call pft1d(bounds, canopystate_inst%tsai_patch, p1); ier = nf90_put_var(ncid, v_tsai, p1)
    call pft1d(bounds, water_inst%waterdiagnosticbulk_inst%fdry_patch, p1); ier = nf90_put_var(ncid, v_fdry, p1)
    call col1d(bounds, atm2lnd_inst%forc_rho_downscaled_col, cground(1,:)); ier = nf90_put_var(ncid, v_rho, cground(1,:))
    call col1d(bounds, atm2lnd_inst%forc_pbot_downscaled_col, cground(1,:)); ier = nf90_put_var(ncid, v_pbot, cground(1,:))
    call pft1d(bounds, temperature_inst%thm_patch, p1); ier = nf90_put_var(ncid, v_thm, p1)
    call pft1d(bounds, water_inst%waterfluxbulk_inst%qflx_tran_veg_patch, p1); ier = nf90_put_var(ncid, v_qtran, p1)
    call pft1d(bounds, energyflux_inst%bsun_patch, p1); ier = nf90_put_var(ncid, v_bsun, p1)
    call pft1d(bounds, energyflux_inst%bsha_patch, p1); ier = nf90_put_var(ncid, v_bsha, p1)
    call pftint1d(bounds, surfalb_inst%nrad_patch, pi1); ier = nf90_put_var(ncid, v_nrad, pi1)
    call pft2d(bounds, nlevcan, surfalb_inst%tlai_z_patch, pcana); ier = nf90_put_var(ncid, v_tlaiz, pcana)
    call pft2d(bounds, nlevcan, surfalb_inst%fsun_z_patch, pcana); ier = nf90_put_var(ncid, v_fsunz, pcana)
    call pft2d(bounds, nlevcan, surfalb_inst%fabd_sun_z_patch, pcana); ier = nf90_put_var(ncid, v_fabdsun, pcana)
    call pft2d(bounds, nlevcan, surfalb_inst%fabd_sha_z_patch, pcana); ier = nf90_put_var(ncid, v_fabdsha, pcana)
    call pft2d(bounds, nlevcan, surfalb_inst%fabi_sun_z_patch, pcana); ier = nf90_put_var(ncid, v_fabisun, pcana)
    call pft2d(bounds, nlevcan, surfalb_inst%fabi_sha_z_patch, pcana); ier = nf90_put_var(ncid, v_fabisha, pcana)
    call pft1d(bounds, surfalb_inst%vcmaxcintsun_patch, p1); ier = nf90_put_var(ncid, v_vcintsun, p1)
    call pft1d(bounds, surfalb_inst%vcmaxcintsha_patch, p1); ier = nf90_put_var(ncid, v_vcintsha, p1)

    ier = nf90_close(ncid)
    if (ier == NF90_NOERR) then
       write(iulog,*) 'btran_oracle_write: wrote ', trim(fname), ' nstep=', nstep
    else
       write(iulog,*) 'btran_oracle_write: close failed: ', trim(nf90_strerror(ier))
    end if
    call cleanup()

  contains

    subroutine cleanup()
      deallocate(p1, pveg, pcana, cground, psoil, pi1)
    end subroutine cleanup

  end subroutine btran_oracle_write

  subroutine pft1d(bounds, src, dst)
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(in) :: src(bounds%begp:)
    real(r8), intent(out) :: dst(:)
    integer :: p
    do p = bounds%begp, bounds%endp
       dst(p-bounds%begp+1) = src(p)
    end do
  end subroutine pft1d

  subroutine pftint1d(bounds, src, dst)
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: src(bounds%begp:)
    integer, intent(out) :: dst(:)
    integer :: p
    do p = bounds%begp, bounds%endp
       dst(p-bounds%begp+1) = src(p)
    end do
  end subroutine pftint1d

  subroutine pft2d(bounds, nlev, src, dst)
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: nlev
    real(r8), intent(in) :: src(bounds%begp:,1:)
    real(r8), intent(out) :: dst(:,:)
    integer :: p, j
    do p = bounds%begp, bounds%endp
       do j = 1, nlev
          dst(j,p-bounds%begp+1) = src(p,j)
       end do
    end do
  end subroutine pft2d

  subroutine col2d(bounds, nlev, src, dst)
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: nlev
    real(r8), intent(in) :: src(bounds%begc:,1:)
    real(r8), intent(out) :: dst(:,:)
    integer :: c, j
    do c = bounds%begc, bounds%endc
       do j = 1, nlev
          dst(j,c-bounds%begc+1) = src(c,j)
       end do
    end do
  end subroutine col2d

  subroutine col1d(bounds, src, dst)
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(in) :: src(bounds%begc:)
    real(r8), intent(out) :: dst(:)
    integer :: c
    do c = bounds%begc, bounds%endc
       dst(c-bounds%begc+1) = src(c)
    end do
  end subroutine col1d

end module btranOracleMod
