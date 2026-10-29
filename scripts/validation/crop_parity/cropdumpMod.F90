module cropdumpMod

  !-----------------------------------------------------------------------
  ! CLM.jl crop-tail parity instrumentation (grain fill, grain fluxes,
  ! vernalization, soybean N fixation).
  !
  ! cropdump_write(bounds, label) writes one NetCDF file per call,
  !     cropdump_<label>_n<nstep>.nc
  ! carrying every patch/column field the three crop tails read or write.
  ! The SAME field set is written at every label; the Julia harness
  ! (scripts/fortran_parity_crop_tails.jl) decides which label supplies the
  ! inputs and which the outputs of each routine under test.
  !
  ! Labels placed by the companion CNDriverMod.F90 SourceMod:
  !   'after_soyfix'    immediately after CNSoyfix (before CNPhenology phase 1):
  !                     CNSoyfix's inputs (fpg, plant_ndemand, sminn, wf, hui,
  !                     gddmaturity) are exactly what it read, and the
  !                     vernalization state (vf/cumvd/hdidx) is the PRE-call value.
  !   'after_cropalloc' immediately after calc_crop_allocation_fractions +
  !                     calc_allometry: hui is the value allocation read, aleaf/
  !                     astem/aroot/arepr are its outputs, vf/cumvd/hdidx are the
  !                     POST-vernalization values.
  !   'after_cnalloc'   immediately after calc_plant_nutrient_competition:
  !                     cpool_to_reproductivec / npool_to_reproductiven are live.
  !   'end_noleach'     end of CNDriverNoLeaching: harvest fluxes
  !                     (repr_grainc_to_food/seed) and updated pools.
  !
  ! nstep windows come from the environment:
  !   CROPDUMP_WINDOWS="lo1:hi1,lo2:hi2,..."   (up to 16 windows, inclusive)
  ! With no variable set nothing is written.
  !
  ! Single-point, masterproc-only (mpiuni); patch/column indices are written
  ! bounds-relative and 1-based.
  !-----------------------------------------------------------------------

  use shr_kind_mod    , only : r8 => shr_kind_r8
  use decompMod       , only : bounds_type
  use clm_varctl      , only : iulog
  use spmdMod         , only : masterproc

  implicit none
  private

  integer, parameter :: maxwin = 16
  integer, save :: nwin = 0
  integer, save :: wlo(maxwin), whi(maxwin)
  logical, save :: initialized = .false.

  public :: cropdump_write
  public :: cropdump_soyfix_sweep

contains

  subroutine init_windows()
    character(len=1024) :: s
    integer :: st, i, j, k, ios
    character(len=64) :: tok
    nwin = 0
    call get_environment_variable('CROPDUMP_WINDOWS', s, status=st)
    if (st == 0 .and. len_trim(s) > 0) then
       i = 1
       do while (i <= len_trim(s) .and. nwin < maxwin)
          j = index(s(i:), ',')
          if (j == 0) then
             tok = s(i:len_trim(s)); i = len_trim(s) + 1
          else
             tok = s(i:i+j-2); i = i + j
          end if
          k = index(tok, ':')
          if (k > 0) then
             nwin = nwin + 1
             read(tok(1:k-1), *, iostat=ios) wlo(nwin)
             if (ios /= 0) error stop 'cropdump: bad CROPDUMP_WINDOWS lo'
             read(tok(k+1:), *, iostat=ios) whi(nwin)
             if (ios /= 0) error stop 'cropdump: bad CROPDUMP_WINDOWS hi'
          end if
       end do
    end if
    initialized = .true.
    if (masterproc) then
       write(iulog,*) 'cropdump_write: ', nwin, ' nstep windows'
       do k = 1, nwin
          write(iulog,'(a,i0,a,i0)') '   cropdump window ', wlo(k), ' .. ', whi(k)
       end do
    end if
  end subroutine init_windows

  !-----------------------------------------------------------------------
  subroutine cropdump_write(bounds, label, cs, ccs, ccf, cns, cnf, crop_inst, &
       soilbiogeochem_state_inst, soilbiogeochem_nitrogenstate_inst,        &
       soilbiogeochem_nitrogenflux_inst, wd, temperature_inst, canopystate_inst)
    ! Instances are passed in (not taken from clm_instMod): clm_instMod ->
    ! CNVegetationFacade -> CNDriverMod -> cropdumpMod would be a module cycle.
    use netcdf
    use CropReprPoolsMod, only : nrepr
    use clm_time_manager, only : get_nstep, get_curr_date
    use PatchType       , only : patch
    use CNVegStateType               , only : cnveg_state_type
    use CNVegCarbonStateType         , only : cnveg_carbonstate_type
    use CNVegCarbonFluxType          , only : cnveg_carbonflux_type
    use CNVegNitrogenStateType       , only : cnveg_nitrogenstate_type
    use CNVegNitrogenFluxType        , only : cnveg_nitrogenflux_type
    use CropType                     , only : crop_type
    use SoilBiogeochemStateType      , only : soilbiogeochem_state_type
    use SoilBiogeochemNitrogenStateType, only : soilbiogeochem_nitrogenstate_type
    use SoilBiogeochemNitrogenFluxType , only : soilbiogeochem_nitrogenflux_type
    use WaterDiagnosticBulkType      , only : waterdiagnosticbulk_type
    use TemperatureType              , only : temperature_type
    use CanopyStateType              , only : canopystate_type
    type(bounds_type), intent(in) :: bounds
    character(len=*) , intent(in) :: label
    type(cnveg_state_type)                 , intent(in) :: cs
    type(cnveg_carbonstate_type)           , intent(in) :: ccs
    type(cnveg_carbonflux_type)            , intent(in) :: ccf
    type(cnveg_nitrogenstate_type)         , intent(in) :: cns
    type(cnveg_nitrogenflux_type)          , intent(in) :: cnf
    type(crop_type)                        , intent(in) :: crop_inst
    type(soilbiogeochem_state_type)        , intent(in) :: soilbiogeochem_state_inst
    type(soilbiogeochem_nitrogenstate_type), intent(in) :: soilbiogeochem_nitrogenstate_inst
    type(soilbiogeochem_nitrogenflux_type) , intent(in) :: soilbiogeochem_nitrogenflux_inst
    type(waterdiagnosticbulk_type)         , intent(in) :: wd
    type(temperature_type)                 , intent(in) :: temperature_inst
    type(canopystate_type)                 , intent(in) :: canopystate_inst

    integer :: nstep, yr, mon, day, tod, k, ncid, ier, d_col, d_pft, v, p, c, np, nc
    logical :: inwin
    character(len=256) :: fname
    real(r8), allocatable :: pb(:), cb(:)
    integer , allocatable :: ib(:)

    if (.not. initialized) call init_windows()
    if (nwin == 0) return
    nstep = get_nstep()
    inwin = .false.
    do k = 1, nwin
       if (nstep >= wlo(k) .and. nstep <= whi(k)) inwin = .true.
    end do
    if (.not. inwin) return
    if (.not. masterproc) return

    call get_curr_date(yr, mon, day, tod)
    np = bounds%endp - bounds%begp + 1
    nc = bounds%endc - bounds%begc + 1
    allocate(pb(np), cb(nc), ib(np))

    write(fname,'(a,a,a,i0,a)') 'cropdump_', trim(label), '_n', nstep, '.nc'
    ier = nf90_create(trim(fname), NF90_CLOBBER, ncid)
    if (ier /= NF90_NOERR) then
       write(iulog,*) 'cropdump_write: nf90_create failed: ', trim(nf90_strerror(ier))
       return
    end if
    ier = nf90_def_dim(ncid, 'column', nc, d_col)
    ier = nf90_def_dim(ncid, 'pft',    np, d_pft)
    ier = nf90_put_att(ncid, NF90_GLOBAL, 'label', trim(label))
    ier = nf90_put_att(ncid, NF90_GLOBAL, 'nrepr', nrepr)
    ier = nf90_enddef(ncid)

    call wi0('nstep', nstep); call wi0('ymd', yr*10000+mon*100+day); call wi0('tod', tod)


    ! ---------------- patch integers ----------------
    do p = bounds%begp, bounds%endp; ib(p-bounds%begp+1) = patch%itype(p); end do
    call wip('ITYPE')
    do p = bounds%begp, bounds%endp; ib(p-bounds%begp+1) = patch%column(p) - bounds%begc + 1; end do
    call wip('PCOLUMN')
    do p = bounds%begp, bounds%endp
       ib(p-bounds%begp+1) = merge(1, 0, crop_inst%croplive_patch(p))
    end do
    call wip('CROPLIVE')
    do p = bounds%begp, bounds%endp; ib(p-bounds%begp+1) = cs%peaklai_patch(p); end do
    call wip('PEAKLAI')
    do p = bounds%begp, bounds%endp; ib(p-bounds%begp+1) = cs%idop_patch(p); end do
    call wip('IDOP')
    do p = bounds%begp, bounds%endp; ib(p-bounds%begp+1) = crop_inst%harvdate_patch(p); end do
    call wip('HARVDATE')
    do p = bounds%begp, bounds%endp; ib(p-bounds%begp+1) = crop_inst%harvest_count(p); end do
    call wip('HARVEST_COUNT')

    ! ---------------- patch reals ----------------
    call wp('WTCOL',        patch%wtcol)
    call wp('HUI',          crop_inst%hui_patch)
    call wp('GDDACCUM',     crop_inst%gddaccum_patch)
    call wp('GDDTSOI',      crop_inst%gddtsoi_patch)
    call wp('VF',           crop_inst%vf_patch)
    call wp('CPHASE',       crop_inst%cphase_patch)
    call wp('GDDMATURITY',  cs%gddmaturity_patch)
    call wp('HUILEAF',      cs%huileaf_patch)
    call wp('HUIGRAIN',     cs%huigrain_patch)
    call wp('CUMVD',        cs%cumvd_patch)
    call wp('HDIDX',        cs%hdidx_patch)
    call wp('ALEAF',        cs%aleaf_patch)
    call wp('ASTEM',        cs%astem_patch)
    call wp('AROOT',        cs%aroot_patch)
    call wp('ALEAFI',       cs%aleafi_patch)
    call wp('ASTEMI',       cs%astemi_patch)
    do p = bounds%begp, bounds%endp; pb(p-bounds%begp+1) = cs%arepr_patch(p,nrepr); end do
    call wpb('AREPR')
    call wp('C_ALLOMETRY',  cs%c_allometry_patch)
    call wp('N_ALLOMETRY',  cs%n_allometry_patch)
    call wp('ONSET_FLAG',   cs%onset_flag_patch)
    call wp('ONSET_COUNTER',cs%onset_counter_patch)
    call wp('OFFSET_FLAG',  cs%offset_flag_patch)
    call wp('OFFSET_COUNTER',cs%offset_counter_patch)
    call wp('BGLFR',        cs%bglfr_patch)
    call wp('BGTR',         cs%bgtr_patch)
    call wp('AVAILC',       ccf%availc_patch)
    call wp('PLANT_CALLOC', ccf%plant_calloc_patch)
    call wp('CPOOL_TO_LEAFC',     ccf%cpool_to_leafc_patch)
    call wp('CPOOL_TO_LIVESTEMC', ccf%cpool_to_livestemc_patch)
    call wp('CPOOL_TO_FROOTC',    ccf%cpool_to_frootc_patch)
    do p = bounds%begp, bounds%endp; pb(p-bounds%begp+1) = ccf%cpool_to_reproductivec_patch(p,nrepr); end do
    call wpb('CPOOL_TO_REPRC')
    do p = bounds%begp, bounds%endp; pb(p-bounds%begp+1) = ccf%cpool_to_reproductivec_storage_patch(p,nrepr); end do
    call wpb('CPOOL_TO_REPRC_STOR')
    do p = bounds%begp, bounds%endp; pb(p-bounds%begp+1) = ccf%repr_grainc_to_food_patch(p,nrepr); end do
    call wpb('REPR_GRAINC_TO_FOOD')
    do p = bounds%begp, bounds%endp; pb(p-bounds%begp+1) = ccf%repr_grainc_to_seed_patch(p,nrepr); end do
    call wpb('REPR_GRAINC_TO_SEED')
    do p = bounds%begp, bounds%endp; pb(p-bounds%begp+1) = ccs%reproductivec_patch(p,nrepr); end do
    call wpb('REPRC')
    call wp('CPOOL',        ccs%cpool_patch)
    call wp('LEAFC',        ccs%leafc_patch)
    call wp('LIVESTEMC',    ccs%livestemc_patch)
    call wp('FROOTC',       ccs%frootc_patch)
    call wp('CROPSEEDC_DEFICIT', ccs%cropseedc_deficit_patch)
    call wp('PLANT_NDEMAND',cnf%plant_ndemand_patch)
    call wp('PLANT_NALLOC', cnf%plant_nalloc_patch)
    call wp('SOYFIXN',      cnf%soyfixn_patch)
    call wp('FERT',         cnf%fert_patch)
    call wp('FERT_COUNTER', cnf%fert_counter_patch)
    do p = bounds%begp, bounds%endp; pb(p-bounds%begp+1) = cnf%npool_to_reproductiven_patch(p,nrepr); end do
    call wpb('NPOOL_TO_REPRN')
    do p = bounds%begp, bounds%endp; pb(p-bounds%begp+1) = cnf%npool_to_reproductiven_storage_patch(p,nrepr); end do
    call wpb('NPOOL_TO_REPRN_STOR')
    do p = bounds%begp, bounds%endp; pb(p-bounds%begp+1) = cnf%repr_grainn_to_food_patch(p,nrepr); end do
    call wpb('REPR_GRAINN_TO_FOOD')
    do p = bounds%begp, bounds%endp; pb(p-bounds%begp+1) = cns%reproductiven_patch(p,nrepr); end do
    call wpb('REPRN')
    call wp('NPOOL',        cns%npool_patch)
    call wp('RETRANSN',     cns%retransn_patch)
    call wp('LEAFN',        cns%leafn_patch)
    call wp('CROPSEEDN_DEFICIT', cns%cropseedn_deficit_patch)
    ! harvest (final-offset-step) fluxes of CNOffsetLitterfall and their pool inputs
    call wp('LIVESTEMN',    cns%livestemn_patch)
    call wp('FROOTN',       cns%frootn_patch)
    do p = bounds%begp, bounds%endp; pb(p-bounds%begp+1) = cnf%repr_grainn_to_seed_patch(p,nrepr); end do
    call wpb('REPR_GRAINN_TO_SEED')
    call wp('LEAFC_TO_LITTER',      ccf%leafc_to_litter_patch)
    call wp('FROOTC_TO_LITTER',     ccf%frootc_to_litter_patch)
    call wp('LIVESTEMC_TO_LITTER',  ccf%livestemc_to_litter_patch)
    call wp('LEAFC_TO_BIOFUELC',    ccf%leafc_to_biofuelc_patch)
    call wp('LIVESTEMC_TO_BIOFUELC',ccf%livestemc_to_biofuelc_patch)
    call wp('LEAFN_TO_LITTER',      cnf%leafn_to_litter_patch)
    call wp('LEAFN_TO_RETRANSN',    cnf%leafn_to_retransn_patch)
    call wp('FROOTN_TO_LITTER',     cnf%frootn_to_litter_patch)
    call wp('LIVESTEMN_TO_LITTER',  cnf%livestemn_to_litter_patch)
    call wp('T_REF2M',      temperature_inst%t_ref2m_patch)
    call wp('T_REF2M_MIN',  temperature_inst%t_ref2m_min_patch)
    call wp('T_REF2M_MAX',  temperature_inst%t_ref2m_max_patch)
    call wp('TLAI',         canopystate_inst%tlai_patch)

    ! ---------------- column reals ----------------
    call wc('FPG',              soilbiogeochem_state_inst%fpg_col)
    call wc('SMINN',            soilbiogeochem_nitrogenstate_inst%sminn_col)
    call wc('WF',               wd%wf_col)
    call wc('SNOW_DEPTH',       wd%snow_depth_col)
    call wc('SOYFIXN_TO_SMINN', soilbiogeochem_nitrogenflux_inst%soyfixn_to_sminn_col)
    call wc('FERT_TO_SMINN',    soilbiogeochem_nitrogenflux_inst%fert_to_sminn_col)


    ier = nf90_close(ncid)
    deallocate(pb, cb, ib)
    write(iulog,'(a,a,a,i0)') 'cropdump_write: wrote ', trim(fname), ' at nstep = ', nstep

  contains

    subroutine wi0(name, val)
      character(len=*), intent(in) :: name
      integer, intent(in) :: val
      ier = nf90_redef(ncid)
      ier = nf90_def_var(ncid, name, NF90_INT, v)
      ier = nf90_enddef(ncid)
      ier = nf90_put_var(ncid, v, val)
    end subroutine wi0

    subroutine wip(name)
      character(len=*), intent(in) :: name
      ier = nf90_redef(ncid)
      ier = nf90_def_var(ncid, name, NF90_INT, (/d_pft/), v)
      ier = nf90_enddef(ncid)
      ier = nf90_put_var(ncid, v, ib)
    end subroutine wip

    subroutine wpb(name)
      character(len=*), intent(in) :: name
      ier = nf90_redef(ncid)
      ier = nf90_def_var(ncid, name, NF90_DOUBLE, (/d_pft/), v)
      ier = nf90_enddef(ncid)
      ier = nf90_put_var(ncid, v, pb)
    end subroutine wpb

    subroutine wp(name, src)
      character(len=*), intent(in) :: name
      real(r8), intent(in) :: src(bounds%begp:)
      integer :: pp
      do pp = bounds%begp, bounds%endp
         pb(pp-bounds%begp+1) = src(pp)
      end do
      call wpb(name)
    end subroutine wp

    subroutine wc(name, src)
      character(len=*), intent(in) :: name
      real(r8), intent(in) :: src(bounds%begc:)
      integer :: cc
      do cc = bounds%begc, bounds%endc
         cb(cc-bounds%begc+1) = src(cc)
      end do
      ier = nf90_redef(ncid)
      ier = nf90_def_var(ncid, name, NF90_DOUBLE, (/d_col/), v)
      ier = nf90_enddef(ncid)
      ier = nf90_put_var(ncid, v, cb)
    end subroutine wc

  end subroutine cropdump_write

  !-----------------------------------------------------------------------
  subroutine cropdump_soyfix_sweep(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
       wd, crop_inst, cs, cnf, soilbiogeochem_state_inst, soilbiogeochem_nitrogenstate_inst, &
       soilbiogeochem_nitrogenflux_inst)
    !
    ! A Fortran-side unit harness for CNSoyfix, run INSIDE the coupled model.
    !
    ! The coupled crop column never reaches the CNSoyfix gate fpg < 1 (FlexibleCN's
    ! non-FUN plant N demand is the Michaelis-Menten form, which collapses with the
    ! mineral pool -- see docs/CROP_LIFECYCLE_MAP.md §3i). So, at each step in the
    ! dump window and only when CROPDUMP_SOYFIX_SWEEP=1, this routine:
    !   1. saves fpg_col, sminn_col, wf_col, hui_patch, soyfixn_patch, soyfixn_to_sminn_col;
    !   2. for each of NCASE injected (fpg, sminn, wf, hui/gddmaturity) combinations,
    !      overwrites those inputs on every column/patch and calls the REAL CNSoyfix;
    !   3. writes inputs + outputs of every case to cropdump_soyfix_sweep_n<nstep>.nc;
    !   4. restores every saved array, so the model trajectory is untouched.
    ! plant_ndemand, croplive, itype, gddmaturity and the patch->column weights stay
    ! the live model values. The case grid crosses every branch of fxn (sminn > 30,
    ! 10 < sminn <= 30, <= 10), fxg (GDDfrac <= .15, (.15,.30], (.30,.55], (.55,.75],
    ! > .75), the fpg gate (fpg < 1 vs >= 1) and the fxw cap (wf/0.85 > 1).
    !
    use netcdf
    use clm_time_manager, only : get_nstep
    use PatchType       , only : patch
    use CNNDynamicsMod  , only : CNSoyfix
    use CNVegStateType               , only : cnveg_state_type
    use CNVegNitrogenFluxType        , only : cnveg_nitrogenflux_type
    use CropType                     , only : crop_type
    use SoilBiogeochemStateType      , only : soilbiogeochem_state_type
    use SoilBiogeochemNitrogenStateType, only : soilbiogeochem_nitrogenstate_type
    use SoilBiogeochemNitrogenFluxType , only : soilbiogeochem_nitrogenflux_type
    use WaterDiagnosticBulkType      , only : waterdiagnosticbulk_type
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_soilc, filter_soilc(:), num_soilp, filter_soilp(:)
    ! intent(in) derived types: only pointer ASSOCIATION is protected; the targets
    ! are overwritten for the sweep and restored before returning.
    type(waterdiagnosticbulk_type)         , intent(in)    :: wd
    type(crop_type)                        , intent(in)    :: crop_inst
    type(cnveg_state_type)                 , intent(in)    :: cs
    type(cnveg_nitrogenflux_type)          , intent(inout) :: cnf
    type(soilbiogeochem_state_type)        , intent(in)    :: soilbiogeochem_state_inst
    type(soilbiogeochem_nitrogenstate_type), intent(in)    :: soilbiogeochem_nitrogenstate_inst
    type(soilbiogeochem_nitrogenflux_type) , intent(inout) :: soilbiogeochem_nitrogenflux_inst

    integer, parameter :: nfpg = 4, nsmn = 3, nwf = 2, ngf = 6
    real(r8), parameter :: fpgv(nfpg) = (/ 0.0_r8, 0.35_r8, 0.8_r8, 1.0_r8 /)
    real(r8), parameter :: smnv(nsmn) = (/ 5.0_r8, 18.0_r8, 40.0_r8 /)
    real(r8), parameter :: wfv(nwf)   = (/ 0.3_r8, 0.95_r8 /)
    real(r8), parameter :: gfv(ngf)   = (/ 0.10_r8, 0.22_r8, 0.42_r8, 0.62_r8, 0.85_r8, -1._r8 /)  ! -1 = live hui
    integer, parameter :: ncase = nfpg*nsmn*nwf*ngf
    character(len=32) :: env_value
    integer :: st, nstep, k, i1, i2, i3, i4, p, c, np, nc, ncid, ier, d_c, d_p, d_k, v
    real(r8), allocatable :: s_fpg(:), s_smn(:), s_wf(:), s_hui(:), s_sfx(:), s_sfx2s(:)
    real(r8), allocatable :: o_fpg(:,:), o_smn(:,:), o_wf(:,:), o_hui(:,:), o_sfx(:,:), o_sfx2s(:,:)
    real(r8), allocatable :: pb(:)
    integer , allocatable :: ib(:)
    character(len=256) :: fname
    logical :: inwin

    call get_environment_variable('CROPDUMP_SOYFIX_SWEEP', env_value, status=st)
    if (st /= 0) return
    if (trim(env_value) /= '1') return
    if (.not. initialized) call init_windows()
    nstep = get_nstep()
    inwin = .false.
    do k = 1, nwin
       if (nstep >= wlo(k) .and. nstep <= whi(k)) inwin = .true.
    end do
    if (.not. inwin) return
    if (.not. masterproc) return

    np = bounds%endp - bounds%begp + 1
    nc = bounds%endc - bounds%begc + 1
    allocate(s_fpg(nc), s_smn(nc), s_wf(nc), s_sfx2s(nc), s_hui(np), s_sfx(np), pb(np), ib(np))
    allocate(o_fpg(nc,ncase), o_smn(nc,ncase), o_wf(nc,ncase), o_sfx2s(nc,ncase))
    allocate(o_hui(np,ncase), o_sfx(np,ncase))

    associate(fpg => soilbiogeochem_state_inst%fpg_col, sminn => soilbiogeochem_nitrogenstate_inst%sminn_col, &
              wf => wd%wf_col, hui => crop_inst%hui_patch, sfx => cnf%soyfixn_patch, &
              sfx2s => soilbiogeochem_nitrogenflux_inst%soyfixn_to_sminn_col, gddm => cs%gddmaturity_patch)
    do c = bounds%begc, bounds%endc
       s_fpg(c-bounds%begc+1) = fpg(c); s_smn(c-bounds%begc+1) = sminn(c)
       s_wf(c-bounds%begc+1) = wf(c);   s_sfx2s(c-bounds%begc+1) = sfx2s(c)
    end do
    do p = bounds%begp, bounds%endp
       s_hui(p-bounds%begp+1) = hui(p); s_sfx(p-bounds%begp+1) = sfx(p)
    end do

    k = 0
    do i1 = 1, nfpg
    do i2 = 1, nsmn
    do i3 = 1, nwf
    do i4 = 1, ngf
       k = k + 1
       do c = bounds%begc, bounds%endc
          fpg(c) = fpgv(i1); sminn(c) = smnv(i2); wf(c) = wfv(i3)
       end do
       do p = bounds%begp, bounds%endp
          if (gfv(i4) > 0._r8 .and. gddm(p) > 0._r8 .and. gddm(p) < 1.e30_r8) then
             hui(p) = gfv(i4) * gddm(p)
          else
             hui(p) = s_hui(p-bounds%begp+1)
          end if
       end do
       call CNSoyfix(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
            wd, crop_inst, cs, cnf, soilbiogeochem_state_inst, &
            soilbiogeochem_nitrogenstate_inst, soilbiogeochem_nitrogenflux_inst)
       do c = bounds%begc, bounds%endc
          o_fpg(c-bounds%begc+1,k) = fpg(c); o_smn(c-bounds%begc+1,k) = sminn(c)
          o_wf(c-bounds%begc+1,k) = wf(c);   o_sfx2s(c-bounds%begc+1,k) = sfx2s(c)
       end do
       do p = bounds%begp, bounds%endp
          o_hui(p-bounds%begp+1,k) = hui(p); o_sfx(p-bounds%begp+1,k) = sfx(p)
       end do
    end do
    end do
    end do
    end do

    ! restore -- the sweep must not perturb the trajectory
    do c = bounds%begc, bounds%endc
       fpg(c) = s_fpg(c-bounds%begc+1); sminn(c) = s_smn(c-bounds%begc+1)
       wf(c) = s_wf(c-bounds%begc+1);   sfx2s(c) = s_sfx2s(c-bounds%begc+1)
    end do
    do p = bounds%begp, bounds%endp
       hui(p) = s_hui(p-bounds%begp+1); sfx(p) = s_sfx(p-bounds%begp+1)
    end do

    write(fname,'(a,i0,a)') 'cropdump_soyfix_sweep_n', nstep, '.nc'
    ier = nf90_create(trim(fname), NF90_CLOBBER, ncid)
    ier = nf90_def_dim(ncid, 'column', nc, d_c)
    ier = nf90_def_dim(ncid, 'pft', np, d_p)
    ier = nf90_def_dim(ncid, 'case', ncase, d_k)
    ier = nf90_enddef(ncid)
    call w2('FPG', o_fpg, d_c); call w2('SMINN', o_smn, d_c); call w2('WF', o_wf, d_c)
    call w2('SOYFIXN_TO_SMINN', o_sfx2s, d_c); call w2('HUI', o_hui, d_p); call w2('SOYFIXN', o_sfx, d_p)
    do p = bounds%begp, bounds%endp; pb(p-bounds%begp+1) = cnf%plant_ndemand_patch(p); end do
    call w1('PLANT_NDEMAND')
    do p = bounds%begp, bounds%endp; pb(p-bounds%begp+1) = gddm(p); end do
    call w1('GDDMATURITY')
    do p = bounds%begp, bounds%endp; pb(p-bounds%begp+1) = patch%wtcol(p); end do
    call w1('WTCOL')
    do p = bounds%begp, bounds%endp; ib(p-bounds%begp+1) = patch%itype(p); end do
    call wi1('ITYPE')
    do p = bounds%begp, bounds%endp; ib(p-bounds%begp+1) = patch%column(p) - bounds%begc + 1; end do
    call wi1('PCOLUMN')
    do p = bounds%begp, bounds%endp; ib(p-bounds%begp+1) = merge(1, 0, crop_inst%croplive_patch(p)); end do
    call wi1('CROPLIVE')
    ier = nf90_redef(ncid); ier = nf90_def_var(ncid, 'nstep', NF90_INT, v); ier = nf90_enddef(ncid)
    ier = nf90_put_var(ncid, v, nstep)
    ier = nf90_close(ncid)
    end associate
    deallocate(s_fpg, s_smn, s_wf, s_sfx2s, s_hui, s_sfx, pb, ib, o_fpg, o_smn, o_wf, o_sfx2s, o_hui, o_sfx)
    write(iulog,'(a,a)') 'cropdump_soyfix_sweep: wrote ', trim(fname)

  contains
    subroutine w2(name, a, d1)
      character(len=*), intent(in) :: name
      real(r8), intent(in) :: a(:,:)
      integer, intent(in) :: d1
      ier = nf90_redef(ncid)
      ier = nf90_def_var(ncid, name, NF90_DOUBLE, (/d1, d_k/), v)
      ier = nf90_enddef(ncid)
      ier = nf90_put_var(ncid, v, a)
    end subroutine w2
    subroutine w1(name)
      character(len=*), intent(in) :: name
      ier = nf90_redef(ncid)
      ier = nf90_def_var(ncid, name, NF90_DOUBLE, (/d_p/), v)
      ier = nf90_enddef(ncid)
      ier = nf90_put_var(ncid, v, pb)
    end subroutine w1
    subroutine wi1(name)
      character(len=*), intent(in) :: name
      ier = nf90_redef(ncid)
      ier = nf90_def_var(ncid, name, NF90_INT, (/d_p/), v)
      ier = nf90_enddef(ncid)
      ier = nf90_put_var(ncid, v, ib)
    end subroutine wi1
  end subroutine cropdump_soyfix_sweep

end module cropdumpMod
