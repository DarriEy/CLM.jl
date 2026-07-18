module bgcdumpMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! CLM.jl parity instrumentation for the METHANE (ch4) and FIRE subsystems.
  !
  ! Neither subsystem has EVER been diffed against Fortran, because the CLM5
  ! restart file does not persist their diagnostics:
  !   * ch4  : conc_ch4/conc_o2 ARE restart vars, but the prod/oxid/aere/ebul
  !            /surface-flux DIAGNOSTICS are not; and in the previously-used
  !            reference config use_lch4 = .false., so the restart carried
  !            ZERO ch4 variables at all.
  !   * fire : only burndate/lfc are restart vars.  farea_burned, nfire, the
  !            baf_* terms and the entire fire C/N flux set are diagnostics.
  !
  ! This module writes those diagnostics directly, at two new boundaries:
  !
  !   'after_fire' -- immediately after EcosystemDynamicsPreDrainage, which is
  !                   where CNFireArea + CNFireFluxes actually run (they live
  !                   inside CNDriverNoLeaching, NOT in PostDrainage).
  !   'after_ch4'  -- immediately after the ch4() call in clm_driver.
  !
  ! Both labels write the SAME full field set; the Julia harness selects which
  ! fields are meaningful at which boundary.  (The fire flux arrays are not
  ! zeroed until the next step's CNZeroFluxes, so they are still live at
  ! 'after_ch4' too; the ch4 fields at 'after_fire' are the PREVIOUS step's and
  ! must not be compared -- the harness does not read them there.)
  !
  ! Gated to an nstep window via BGCDUMP_NSTEP_LO / BGCDUMP_NSTEP_HI so a long
  ! spin-up run does not flood the run directory.
  !
  ! Modeled on the proven pdumpMod in this same SourceMods directory.
  !
  ! NOTE: requires the companion SourceMods copy of ch4Mod.F90, in which the
  ! ch4_type component pointers are `public` rather than `private` (a pure
  ! visibility change; no physics is touched).
  !-----------------------------------------------------------------------

  use shr_kind_mod    , only : r8 => shr_kind_r8
  use decompMod       , only : bounds_type
  use clm_varctl      , only : iulog
  use spmdMod         , only : masterproc

  implicit none
  private

  integer, save :: bgcdump_nstep_lo = 999999999
  integer, save :: bgcdump_nstep_hi = -1
  logical, save :: bgcdump_window_initialized = .false.

  public :: bgcdump_write

contains

  !-----------------------------------------------------------------------
  subroutine bgcdump_write(bounds, label)
    !
    ! !USES:
    use netcdf
    use clm_varpar      , only : nlevsoi, nlevdecomp_full, ndecomp_pools
    use clm_time_manager, only : get_nstep, get_curr_date
    use clm_instMod     , only : ch4_inst, bgc_vegetation_inst, &
                                 soilbiogeochem_carbonflux_inst, &
                                 crop_inst, soilbiogeochem_nitrogenflux_inst
    use clm_varctl      , only : use_crop
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    character(len=*) , intent(in) :: label
    !
    ! !LOCAL VARIABLES:
    real(r8), parameter :: fillval = 1.e+36_r8
    integer  :: nstep, yr, mon, day, tod, ymd
    integer  :: env_status, read_status
    integer  :: ncol, npft
    integer  :: ncid, ier
    integer  :: d_col, d_pft, d_soi, d_dcmp, d_pool
    integer  :: v_nstep, v_ymd, v_tod
    character(len=256) :: fname
    character(len=32)  :: env_value

    ! ---- ch4 2D (levsoi, col)
    integer :: v_cch4s, v_cch4u, v_co2s, v_co2u
    integer :: v_prods, v_produ, v_oxids, v_oxidu
    integer :: v_aeres, v_aereu, v_ebuls, v_ebulu
    integer :: v_trans, v_tranu
    integer :: v_o2ss, v_o2su, v_c4ss, v_c4su
    integer :: v_o2ds, v_o2du
    ! ---- ch4 1D (col)
    integer :: v_sdifs, v_sdifu, v_saers, v_saeru, v_sebus, v_sebuu
    integer :: v_ebts, v_ebtu, v_sflux, v_dfsat
    integer :: v_finu, v_finlag, v_finpre, v_tcch4, v_tcch4b
    integer :: v_zwtu, v_lsl, v_qsl, v_sif, v_gcond
    ! ---- fire 1D (col)
    integer :: v_fab, v_nfire, v_bafc, v_bafp, v_lfc, v_fbac, v_fbac1
    integer :: v_lfc2, v_dwts   ! deforestation-fire (dtrotr) branch
    integer :: v_fuelc, v_fuelcc, v_cropf, v_lgdp, v_lgdp1, v_lpop
    integer :: v_fsr, v_fd, v_rootc, v_lfwt, v_tr1, v_tr2, v_dtr, v_wtlf
    integer :: v_fclossc, v_fnlossc, v_somcfire
    ! ---- fire 2D
    integer :: v_fmc2cwdc, v_fmn2cwdn, v_mdc2fire
    ! ---- fire patch
    integer :: v_burndate, v_fclossp, v_fnlossp
    integer :: v_mlfc, v_mfrc, v_mlsc, v_mdsc, v_mlfc2l
    ! ---- crop patch (CropType) + crop N (backlog A3)
    ! Written only when use_crop; on a non-crop run crop_inst's arrays are
    ! unallocated, so every access below must stay behind the guard.
    integer :: v_croplive, v_harvdate, v_cphase, v_hui, v_gddaccum
    integer :: v_gddtsoi, v_vf, v_latbaset, v_fertnitro, v_nyrscrop
    integer :: v_fert, v_soyfixn, v_fert2sminn, v_soyfix2sminn
    integer :: pp

    real(r8), allocatable :: s2(:,:)      ! (nlevsoi, ncol)
    real(r8), allocatable :: d2(:,:)      ! (nlevdecomp_full, ncol)
    real(r8), allocatable :: pl2(:,:)     ! (ndecomp_pools, ncol)
    real(r8), allocatable :: c1(:)        ! (ncol)
    real(r8), allocatable :: p1(:)        ! (npft)
    integer , allocatable :: pi1(:)       ! (npft)
    !-----------------------------------------------------------------------

    if (.not. bgcdump_window_initialized) then
       call get_environment_variable('BGCDUMP_NSTEP_LO', env_value, status=env_status)
       if (env_status == 0 .and. len_trim(env_value) > 0) then
          read(env_value, *, iostat=read_status) bgcdump_nstep_lo
          if (read_status /= 0) error stop 'invalid BGCDUMP_NSTEP_LO'
       end if
       call get_environment_variable('BGCDUMP_NSTEP_HI', env_value, status=env_status)
       if (env_status == 0 .and. len_trim(env_value) > 0) then
          read(env_value, *, iostat=read_status) bgcdump_nstep_hi
          if (read_status /= 0) error stop 'invalid BGCDUMP_NSTEP_HI'
       end if
       bgcdump_window_initialized = .true.
       if (masterproc) then
          write(iulog,'(a,i0,a,i0)') 'bgcdump_write: nstep window ', &
               bgcdump_nstep_lo, ' .. ', bgcdump_nstep_hi
       end if
    end if

    nstep = get_nstep()
    if (nstep < bgcdump_nstep_lo .or. nstep > bgcdump_nstep_hi) return
    if (.not. masterproc) return

    call get_curr_date(yr, mon, day, tod)
    ymd = yr*10000 + mon*100 + day

    ncol = bounds%endc - bounds%begc + 1
    npft = bounds%endp - bounds%begp + 1

    allocate(s2(nlevsoi, ncol), d2(nlevdecomp_full, ncol), pl2(ndecomp_pools, ncol))
    allocate(c1(ncol), p1(npft), pi1(npft))

    write(fname,'(a,a,a,i0,a)') 'bgcdump_', trim(label), '_n', nstep, '.nc'

    ier = nf90_create(trim(fname), NF90_CLOBBER, ncid)
    if (ier /= NF90_NOERR) then
       write(iulog,*) 'bgcdump_write: nf90_create failed: ', trim(nf90_strerror(ier))
       deallocate(s2, d2, pl2, c1, p1, pi1)
       return
    end if

    ier = nf90_def_dim(ncid, 'column',  ncol,            d_col)
    ier = nf90_def_dim(ncid, 'pft',     npft,            d_pft)
    ier = nf90_def_dim(ncid, 'levsoi',  nlevsoi,         d_soi)
    ier = nf90_def_dim(ncid, 'levdcmp', nlevdecomp_full, d_dcmp)
    ier = nf90_def_dim(ncid, 'dpool',   ndecomp_pools,   d_pool)

    ! ================= ch4, level-resolved =================
    ier = nf90_def_var(ncid,'CONC_CH4_SAT',      NF90_DOUBLE,(/d_soi,d_col/), v_cch4s)
    ier = nf90_def_var(ncid,'CONC_CH4_UNSAT',    NF90_DOUBLE,(/d_soi,d_col/), v_cch4u)
    ier = nf90_def_var(ncid,'CONC_O2_SAT',       NF90_DOUBLE,(/d_soi,d_col/), v_co2s)
    ier = nf90_def_var(ncid,'CONC_O2_UNSAT',     NF90_DOUBLE,(/d_soi,d_col/), v_co2u)
    ier = nf90_def_var(ncid,'CH4_PROD_SAT',      NF90_DOUBLE,(/d_soi,d_col/), v_prods)
    ier = nf90_def_var(ncid,'CH4_PROD_UNSAT',    NF90_DOUBLE,(/d_soi,d_col/), v_produ)
    ier = nf90_def_var(ncid,'CH4_OXID_SAT',      NF90_DOUBLE,(/d_soi,d_col/), v_oxids)
    ier = nf90_def_var(ncid,'CH4_OXID_UNSAT',    NF90_DOUBLE,(/d_soi,d_col/), v_oxidu)
    ier = nf90_def_var(ncid,'CH4_AERE_SAT',      NF90_DOUBLE,(/d_soi,d_col/), v_aeres)
    ier = nf90_def_var(ncid,'CH4_AERE_UNSAT',    NF90_DOUBLE,(/d_soi,d_col/), v_aereu)
    ier = nf90_def_var(ncid,'CH4_EBUL_SAT',      NF90_DOUBLE,(/d_soi,d_col/), v_ebuls)
    ier = nf90_def_var(ncid,'CH4_EBUL_UNSAT',    NF90_DOUBLE,(/d_soi,d_col/), v_ebulu)
    ier = nf90_def_var(ncid,'CH4_TRAN_SAT',      NF90_DOUBLE,(/d_soi,d_col/), v_trans)
    ier = nf90_def_var(ncid,'CH4_TRAN_UNSAT',    NF90_DOUBLE,(/d_soi,d_col/), v_tranu)
    ier = nf90_def_var(ncid,'O2STRESS_SAT',      NF90_DOUBLE,(/d_soi,d_col/), v_o2ss)
    ier = nf90_def_var(ncid,'O2STRESS_UNSAT',    NF90_DOUBLE,(/d_soi,d_col/), v_o2su)
    ier = nf90_def_var(ncid,'CH4STRESS_SAT',     NF90_DOUBLE,(/d_soi,d_col/), v_c4ss)
    ier = nf90_def_var(ncid,'CH4STRESS_UNSAT',   NF90_DOUBLE,(/d_soi,d_col/), v_c4su)
    ier = nf90_def_var(ncid,'O2_DECOMP_SAT',     NF90_DOUBLE,(/d_soi,d_col/), v_o2ds)
    ier = nf90_def_var(ncid,'O2_DECOMP_UNSAT',   NF90_DOUBLE,(/d_soi,d_col/), v_o2du)

    ! ================= ch4, column =================
    ier = nf90_def_var(ncid,'CH4_SURF_DIFF_SAT',   NF90_DOUBLE,(/d_col/), v_sdifs)
    ier = nf90_def_var(ncid,'CH4_SURF_DIFF_UNSAT', NF90_DOUBLE,(/d_col/), v_sdifu)
    ier = nf90_def_var(ncid,'CH4_SURF_AERE_SAT',   NF90_DOUBLE,(/d_col/), v_saers)
    ier = nf90_def_var(ncid,'CH4_SURF_AERE_UNSAT', NF90_DOUBLE,(/d_col/), v_saeru)
    ier = nf90_def_var(ncid,'CH4_SURF_EBUL_SAT',   NF90_DOUBLE,(/d_col/), v_sebus)
    ier = nf90_def_var(ncid,'CH4_SURF_EBUL_UNSAT', NF90_DOUBLE,(/d_col/), v_sebuu)
    ier = nf90_def_var(ncid,'CH4_EBUL_TOTAL_SAT',  NF90_DOUBLE,(/d_col/), v_ebts)
    ier = nf90_def_var(ncid,'CH4_EBUL_TOTAL_UNSAT',NF90_DOUBLE,(/d_col/), v_ebtu)
    ier = nf90_def_var(ncid,'CH4_SURF_FLUX_TOT',   NF90_DOUBLE,(/d_col/), v_sflux)
    ier = nf90_def_var(ncid,'CH4_DFSAT_FLUX',      NF90_DOUBLE,(/d_col/), v_dfsat)
    ier = nf90_def_var(ncid,'FINUNDATED',          NF90_DOUBLE,(/d_col/), v_finu)
    ier = nf90_def_var(ncid,'FINUNDATED_LAG',      NF90_DOUBLE,(/d_col/), v_finlag)
    ier = nf90_def_var(ncid,'FINUNDATED_PRE_SNOW', NF90_DOUBLE,(/d_col/), v_finpre)
    ier = nf90_def_var(ncid,'TOTCOLCH4',           NF90_DOUBLE,(/d_col/), v_tcch4)
    ier = nf90_def_var(ncid,'TOTCOLCH4_BEF',       NF90_DOUBLE,(/d_col/), v_tcch4b)
    ier = nf90_def_var(ncid,'ZWT_CH4_UNSAT',       NF90_DOUBLE,(/d_col/), v_zwtu)
    ier = nf90_def_var(ncid,'LAYER_SAT_LAG',       NF90_DOUBLE,(/d_soi,d_col/), v_lsl)
    ier = nf90_def_var(ncid,'QFLX_SURF_LAG',       NF90_DOUBLE,(/d_col/), v_qsl)
    ier = nf90_def_var(ncid,'SIF',                 NF90_DOUBLE,(/d_col/), v_sif)
    ier = nf90_def_var(ncid,'GRND_CH4_COND',       NF90_DOUBLE,(/d_col/), v_gcond)

    ! ================= fire, column =================
    ier = nf90_def_var(ncid,'FAREA_BURNED', NF90_DOUBLE,(/d_col/), v_fab)
    ier = nf90_def_var(ncid,'NFIRE',        NF90_DOUBLE,(/d_col/), v_nfire)
    ier = nf90_def_var(ncid,'BAF_CROP',     NF90_DOUBLE,(/d_col/), v_bafc)
    ier = nf90_def_var(ncid,'BAF_PEATF',    NF90_DOUBLE,(/d_col/), v_bafp)
    ier = nf90_def_var(ncid,'LFC',          NF90_DOUBLE,(/d_col/), v_lfc)
    ier = nf90_def_var(ncid,'FBAC',         NF90_DOUBLE,(/d_col/), v_fbac)
    ier = nf90_def_var(ncid,'FBAC1',        NF90_DOUBLE,(/d_col/), v_fbac1)
    ier = nf90_def_var(ncid,'LFC2',         NF90_DOUBLE,(/d_col/), v_lfc2)
    ier = nf90_def_var(ncid,'FUELC',        NF90_DOUBLE,(/d_col/), v_fuelc)
    ier = nf90_def_var(ncid,'FUELC_CROP',   NF90_DOUBLE,(/d_col/), v_fuelcc)
    ier = nf90_def_var(ncid,'CROPF',        NF90_DOUBLE,(/d_col/), v_cropf)
    ier = nf90_def_var(ncid,'LGDP',         NF90_DOUBLE,(/d_col/), v_lgdp)
    ier = nf90_def_var(ncid,'LGDP1',        NF90_DOUBLE,(/d_col/), v_lgdp1)
    ier = nf90_def_var(ncid,'LPOP',         NF90_DOUBLE,(/d_col/), v_lpop)
    ier = nf90_def_var(ncid,'FSR_COL',      NF90_DOUBLE,(/d_col/), v_fsr)
    ier = nf90_def_var(ncid,'FD_COL',       NF90_DOUBLE,(/d_col/), v_fd)
    ier = nf90_def_var(ncid,'ROOTC_COL',    NF90_DOUBLE,(/d_col/), v_rootc)
    ier = nf90_def_var(ncid,'LFWT',         NF90_DOUBLE,(/d_col/), v_lfwt)
    ier = nf90_def_var(ncid,'TROTR1',       NF90_DOUBLE,(/d_col/), v_tr1)
    ier = nf90_def_var(ncid,'TROTR2',       NF90_DOUBLE,(/d_col/), v_tr2)
    ier = nf90_def_var(ncid,'DTROTR',       NF90_DOUBLE,(/d_col/), v_dtr)
    ier = nf90_def_var(ncid,'WTLF',         NF90_DOUBLE,(/d_col/), v_wtlf)
    ier = nf90_def_var(ncid,'FIRE_CLOSS_COL', NF90_DOUBLE,(/d_col/), v_fclossc)
    ier = nf90_def_var(ncid,'FIRE_NLOSS_COL', NF90_DOUBLE,(/d_col/), v_fnlossc)
    ier = nf90_def_var(ncid,'SOMC_FIRE',      NF90_DOUBLE,(/d_col/), v_somcfire)

    ier = nf90_def_var(ncid,'FIRE_MORTALITY_C_TO_CWDC', NF90_DOUBLE,(/d_dcmp,d_col/), v_fmc2cwdc)
    ier = nf90_def_var(ncid,'FIRE_MORTALITY_N_TO_CWDN', NF90_DOUBLE,(/d_dcmp,d_col/), v_fmn2cwdn)
    ier = nf90_def_var(ncid,'M_DECOMP_CPOOLS_TO_FIRE',  NF90_DOUBLE,(/d_pool,d_col/), v_mdc2fire)

    ! ================= fire, patch =================
    ier = nf90_def_var(ncid,'BURNDATE',            NF90_INT,   (/d_pft/), v_burndate)
    ier = nf90_def_var(ncid,'FIRE_CLOSS_PATCH',    NF90_DOUBLE,(/d_pft/), v_fclossp)
    ier = nf90_def_var(ncid,'FIRE_NLOSS_PATCH',    NF90_DOUBLE,(/d_pft/), v_fnlossp)
    ier = nf90_def_var(ncid,'M_LEAFC_TO_FIRE',     NF90_DOUBLE,(/d_pft/), v_mlfc)
    ier = nf90_def_var(ncid,'M_FROOTC_TO_FIRE',    NF90_DOUBLE,(/d_pft/), v_mfrc)
    ier = nf90_def_var(ncid,'M_LIVESTEMC_TO_FIRE', NF90_DOUBLE,(/d_pft/), v_mlsc)
    ier = nf90_def_var(ncid,'M_DEADSTEMC_TO_FIRE', NF90_DOUBLE,(/d_pft/), v_mdsc)
    ier = nf90_def_var(ncid,'M_LEAFC_TO_LITTER_FIRE', NF90_DOUBLE,(/d_pft/), v_mlfc2l)
    ier = nf90_def_var(ncid,'DWT_SMOOTHED',        NF90_DOUBLE,(/d_pft/), v_dwts)

    ! ================= crop phenology + crop N (backlog A3/A4) =================
    ! The crop lifecycle state CLM.jl must reproduce: the planting flag and date,
    ! the phase, and the degree-day/heat-unit accumulators that drive them. Plus
    ! the two crop N fluxes (CNNFert's p2c target and CNSoyfix's patch output),
    ! which is what the still-unwired n_fert!/n_soyfix! need to be diffed against.
    if (use_crop) then
       ier = nf90_def_var(ncid,'CROPLIVE',   NF90_INT,   (/d_pft/), v_croplive)
       ier = nf90_def_var(ncid,'HARVDATE',   NF90_INT,   (/d_pft/), v_harvdate)
       ier = nf90_def_var(ncid,'NYRS_CROP_ACTIVE', NF90_INT, (/d_pft/), v_nyrscrop)
       ier = nf90_def_var(ncid,'CPHASE',     NF90_DOUBLE,(/d_pft/), v_cphase)
       ier = nf90_def_var(ncid,'HUI',        NF90_DOUBLE,(/d_pft/), v_hui)
       ier = nf90_def_var(ncid,'GDDACCUM',   NF90_DOUBLE,(/d_pft/), v_gddaccum)
       ier = nf90_def_var(ncid,'GDDTSOI',    NF90_DOUBLE,(/d_pft/), v_gddtsoi)
       ier = nf90_def_var(ncid,'VF',         NF90_DOUBLE,(/d_pft/), v_vf)
       ier = nf90_def_var(ncid,'LATBASET',   NF90_DOUBLE,(/d_pft/), v_latbaset)
       ier = nf90_def_var(ncid,'FERTNITRO',  NF90_DOUBLE,(/d_pft/), v_fertnitro)
       ier = nf90_def_var(ncid,'FERT',       NF90_DOUBLE,(/d_pft/), v_fert)
       ier = nf90_def_var(ncid,'SOYFIXN',    NF90_DOUBLE,(/d_pft/), v_soyfixn)
       ier = nf90_def_var(ncid,'FERT_TO_SMINN',    NF90_DOUBLE,(/d_col/), v_fert2sminn)
       ier = nf90_def_var(ncid,'SOYFIXN_TO_SMINN', NF90_DOUBLE,(/d_col/), v_soyfix2sminn)
    end if

    ier = nf90_def_var(ncid,'nstep', NF90_INT, v_nstep)
    ier = nf90_def_var(ncid,'ymd',   NF90_INT, v_ymd)
    ier = nf90_def_var(ncid,'tod',   NF90_INT, v_tod)

    ier = nf90_enddef(ncid)

    ! ================= write ch4 level-resolved =================
    call put2d(ncid, v_cch4s, bounds, ch4_inst%conc_ch4_sat_col,        nlevsoi, s2)
    call put2d(ncid, v_cch4u, bounds, ch4_inst%conc_ch4_unsat_col,      nlevsoi, s2)
    call put2d(ncid, v_co2s,  bounds, ch4_inst%conc_o2_sat_col,         nlevsoi, s2)
    call put2d(ncid, v_co2u,  bounds, ch4_inst%conc_o2_unsat_col,       nlevsoi, s2)
    call put2d(ncid, v_prods, bounds, ch4_inst%ch4_prod_depth_sat_col,  nlevsoi, s2)
    call put2d(ncid, v_produ, bounds, ch4_inst%ch4_prod_depth_unsat_col,nlevsoi, s2)
    call put2d(ncid, v_oxids, bounds, ch4_inst%ch4_oxid_depth_sat_col,  nlevsoi, s2)
    call put2d(ncid, v_oxidu, bounds, ch4_inst%ch4_oxid_depth_unsat_col,nlevsoi, s2)
    call put2d(ncid, v_aeres, bounds, ch4_inst%ch4_aere_depth_sat_col,  nlevsoi, s2)
    call put2d(ncid, v_aereu, bounds, ch4_inst%ch4_aere_depth_unsat_col,nlevsoi, s2)
    call put2d(ncid, v_ebuls, bounds, ch4_inst%ch4_ebul_depth_sat_col,  nlevsoi, s2)
    call put2d(ncid, v_ebulu, bounds, ch4_inst%ch4_ebul_depth_unsat_col,nlevsoi, s2)
    call put2d(ncid, v_trans, bounds, ch4_inst%ch4_tran_depth_sat_col,  nlevsoi, s2)
    call put2d(ncid, v_tranu, bounds, ch4_inst%ch4_tran_depth_unsat_col,nlevsoi, s2)
    call put2d(ncid, v_o2ss,  bounds, ch4_inst%o2stress_sat_col,        nlevsoi, s2)
    call put2d(ncid, v_o2su,  bounds, ch4_inst%o2stress_unsat_col,      nlevsoi, s2)
    call put2d(ncid, v_c4ss,  bounds, ch4_inst%ch4stress_sat_col,       nlevsoi, s2)
    call put2d(ncid, v_c4su,  bounds, ch4_inst%ch4stress_unsat_col,     nlevsoi, s2)
    call put2d(ncid, v_o2ds,  bounds, ch4_inst%o2_decomp_depth_sat_col, nlevsoi, s2)
    call put2d(ncid, v_o2du,  bounds, ch4_inst%o2_decomp_depth_unsat_col,nlevsoi,s2)

    ! ================= write ch4 column =================
    call put1d(ncid, v_sdifs,  bounds, ch4_inst%ch4_surf_diff_sat_col,   c1)
    call put1d(ncid, v_sdifu,  bounds, ch4_inst%ch4_surf_diff_unsat_col, c1)
    call put1d(ncid, v_saers,  bounds, ch4_inst%ch4_surf_aere_sat_col,   c1)
    call put1d(ncid, v_saeru,  bounds, ch4_inst%ch4_surf_aere_unsat_col, c1)
    call put1d(ncid, v_sebus,  bounds, ch4_inst%ch4_surf_ebul_sat_col,   c1)
    call put1d(ncid, v_sebuu,  bounds, ch4_inst%ch4_surf_ebul_unsat_col, c1)
    call put1d(ncid, v_ebts,   bounds, ch4_inst%ch4_ebul_total_sat_col,  c1)
    call put1d(ncid, v_ebtu,   bounds, ch4_inst%ch4_ebul_total_unsat_col,c1)
    call put1d(ncid, v_sflux,  bounds, ch4_inst%ch4_surf_flux_tot_col,   c1)
    call put1d(ncid, v_dfsat,  bounds, ch4_inst%ch4_dfsat_flux_col,      c1)
    call put1d(ncid, v_finu,   bounds, ch4_inst%finundated_col,          c1)
    call put1d(ncid, v_finlag, bounds, ch4_inst%finundated_lag_col,      c1)
    call put1d(ncid, v_finpre, bounds, ch4_inst%finundated_pre_snow_col, c1)
    call put1d(ncid, v_tcch4,  bounds, ch4_inst%totcolch4_col,           c1)
    call put1d(ncid, v_tcch4b, bounds, ch4_inst%totcolch4_bef_col,       c1)
    call put1d(ncid, v_zwtu,   bounds, ch4_inst%zwt_ch4_unsat_col,       c1)
    call put2d(ncid, v_lsl, bounds, ch4_inst%layer_sat_lag_col, nlevsoi, s2)
    call put1d(ncid, v_qsl,    bounds, ch4_inst%qflx_surf_lag_col,       c1)
    call put1d(ncid, v_sif,    bounds, ch4_inst%sif_col,                 c1)
    call put1d(ncid, v_gcond,  bounds, ch4_inst%grnd_ch4_cond_col,       c1)

    ! ================= write fire column =================
    ! NB: fuelc/fuelc_crop/rootc live on cnveg_carbonSTATE_inst, not cnveg_state_inst.
    associate( cvs => bgc_vegetation_inst%cnveg_state_inst,        &
               cvc => bgc_vegetation_inst%cnveg_carbonflux_inst,   &
               cvn => bgc_vegetation_inst%cnveg_nitrogenflux_inst, &
               cvcs => bgc_vegetation_inst%cnveg_carbonstate_inst )

    call put1d(ncid, v_fab,    bounds, cvs%farea_burned_col, c1)
    call put1d(ncid, v_nfire,  bounds, cvs%nfire_col,        c1)
    call put1d(ncid, v_bafc,   bounds, cvs%baf_crop_col,     c1)
    call put1d(ncid, v_bafp,   bounds, cvs%baf_peatf_col,    c1)
    call put1d(ncid, v_lfc,    bounds, cvs%lfc_col,          c1)
    call put1d(ncid, v_fbac,   bounds, cvs%fbac_col,         c1)
    call put1d(ncid, v_fbac1,  bounds, cvs%fbac1_col,        c1)
    call put1d(ncid, v_lfc2,   bounds, cvs%lfc2_col,         c1)
    call put1d(ncid, v_fuelc,  bounds, cvcs%fuelc_col,        c1)
    call put1d(ncid, v_fuelcc, bounds, cvcs%fuelc_crop_col,   c1)
    call put1d(ncid, v_cropf,  bounds, cvs%cropf_col,        c1)
    call put1d(ncid, v_lgdp,   bounds, cvs%lgdp_col,         c1)
    call put1d(ncid, v_lgdp1,  bounds, cvs%lgdp1_col,        c1)
    call put1d(ncid, v_lpop,   bounds, cvs%lpop_col,         c1)
    call put1d(ncid, v_fsr,    bounds, cvs%fsr_col,          c1)
    call put1d(ncid, v_fd,     bounds, cvs%fd_col,           c1)
    call put1d(ncid, v_rootc,  bounds, cvcs%rootc_col,        c1)
    call put1d(ncid, v_lfwt,   bounds, cvs%lfwt_col,         c1)
    call put1d(ncid, v_tr1,    bounds, cvs%trotr1_col,       c1)
    call put1d(ncid, v_tr2,    bounds, cvs%trotr2_col,       c1)
    call put1d(ncid, v_dtr,    bounds, cvs%dtrotr_col,       c1)
    call put1d(ncid, v_wtlf,   bounds, cvs%wtlf_col,         c1)

    call put1d(ncid, v_fclossc,  bounds, cvc%fire_closs_col, c1)
    call put1d(ncid, v_fnlossc,  bounds, cvn%fire_nloss_col, c1)
    call put1d(ncid, v_somcfire, bounds, soilbiogeochem_carbonflux_inst%somc_fire_col, c1)

    call put2d(ncid, v_fmc2cwdc, bounds, cvc%fire_mortality_c_to_cwdc_col, nlevdecomp_full, d2)
    call put2d(ncid, v_fmn2cwdn, bounds, cvn%fire_mortality_n_to_cwdn_col, nlevdecomp_full, d2)
    call put2d(ncid, v_mdc2fire, bounds, cvc%m_decomp_cpools_to_fire_col,  ndecomp_pools,   pl2)

    ! ================= write fire patch =================
    call puti1p(ncid, v_burndate, bounds, cvs%burndate_patch, pi1)
    call put1p(ncid, v_fclossp, bounds, cvc%fire_closs_patch,            p1)
    call put1p(ncid, v_fnlossp, bounds, cvn%fire_nloss_patch,            p1)
    call put1p(ncid, v_mlfc,    bounds, cvc%m_leafc_to_fire_patch,       p1)
    call put1p(ncid, v_mfrc,    bounds, cvc%m_frootc_to_fire_patch,      p1)
    call put1p(ncid, v_mlsc,    bounds, cvc%m_livestemc_to_fire_patch,   p1)
    call put1p(ncid, v_mdsc,    bounds, cvc%m_deadstemc_to_fire_patch,   p1)
    call put1p(ncid, v_mlfc2l,  bounds, cvc%m_leafc_to_litter_fire_patch,p1)
    call put1p(ncid, v_dwts,    bounds, cvs%dwt_smoothed_patch,          p1)

    ! ================= write crop phenology + crop N =================
    if (use_crop) then
       ! croplive is LOGICAL; NetCDF has no logical type, so map to 0/1 by hand
       ! rather than via puti1p (which expects an integer array).
       do pp = bounds%begp, bounds%endp
          pi1(pp - bounds%begp + 1) = merge(1, 0, crop_inst%croplive_patch(pp))
       end do
       ier = nf90_put_var(ncid, v_croplive, pi1)

       call puti1p(ncid, v_harvdate,  bounds, crop_inst%harvdate_patch,        pi1)
       call puti1p(ncid, v_nyrscrop,  bounds, crop_inst%nyrs_crop_active_patch,pi1)
       call put1p(ncid, v_cphase,     bounds, crop_inst%cphase_patch,     p1)
       call put1p(ncid, v_hui,        bounds, crop_inst%hui_patch,        p1)
       call put1p(ncid, v_gddaccum,   bounds, crop_inst%gddaccum_patch,   p1)
       call put1p(ncid, v_gddtsoi,    bounds, crop_inst%gddtsoi_patch,    p1)
       call put1p(ncid, v_vf,         bounds, crop_inst%vf_patch,         p1)
       call put1p(ncid, v_latbaset,   bounds, crop_inst%latbaset_patch,   p1)
       call put1p(ncid, v_fertnitro,  bounds, crop_inst%fertnitro_patch,  p1)
       ! CNPhenology's fertilizer output (the producer CLM.jl does not port)
       ! and CNSoyfix's patch output, plus their two p2c column targets.
       call put1p(ncid, v_fert,       bounds, cvn%fert_patch,             p1)
       call put1p(ncid, v_soyfixn,    bounds, cvn%soyfixn_patch,          p1)
       call put1d(ncid, v_fert2sminn,   bounds, &
            soilbiogeochem_nitrogenflux_inst%fert_to_sminn_col,    c1)
       call put1d(ncid, v_soyfix2sminn, bounds, &
            soilbiogeochem_nitrogenflux_inst%soyfixn_to_sminn_col, c1)
    end if

    end associate

    ier = nf90_put_var(ncid, v_nstep, nstep)
    ier = nf90_put_var(ncid, v_ymd,   ymd)
    ier = nf90_put_var(ncid, v_tod,   tod)

    ier = nf90_close(ncid)

    deallocate(s2, d2, pl2, c1, p1, pi1)

    if (masterproc) then
       write(iulog,'(a,a,a,i0)') 'bgcdump_write: wrote ', trim(fname), ' at nstep = ', nstep
    end if

  end subroutine bgcdump_write

  !-----------------------------------------------------------------------
  ! helpers: gather bounds-offset arrays into 1-based buffers and write.
  ! SPVAL/NaN are written verbatim -- the Julia harness must see them.
  !-----------------------------------------------------------------------
  subroutine put1d(ncid, varid, bounds, src, buf)
    use netcdf
    integer          , intent(in)    :: ncid, varid
    type(bounds_type), intent(in)    :: bounds
    real(r8)         , intent(in)    :: src(bounds%begc:)
    real(r8)         , intent(inout) :: buf(:)
    integer :: c, ier
    do c = bounds%begc, bounds%endc
       buf(c - bounds%begc + 1) = src(c)
    end do
    ier = nf90_put_var(ncid, varid, buf)
  end subroutine put1d

  subroutine put1p(ncid, varid, bounds, src, buf)
    use netcdf
    integer          , intent(in)    :: ncid, varid
    type(bounds_type), intent(in)    :: bounds
    real(r8)         , intent(in)    :: src(bounds%begp:)
    real(r8)         , intent(inout) :: buf(:)
    integer :: p, ier
    do p = bounds%begp, bounds%endp
       buf(p - bounds%begp + 1) = src(p)
    end do
    ier = nf90_put_var(ncid, varid, buf)
  end subroutine put1p

  subroutine puti1p(ncid, varid, bounds, src, buf)
    use netcdf
    integer          , intent(in)    :: ncid, varid
    type(bounds_type), intent(in)    :: bounds
    integer          , intent(in)    :: src(bounds%begp:)
    integer          , intent(inout) :: buf(:)
    integer :: p, ier
    do p = bounds%begp, bounds%endp
       buf(p - bounds%begp + 1) = src(p)
    end do
    ier = nf90_put_var(ncid, varid, buf)
  end subroutine puti1p

  subroutine put2d(ncid, varid, bounds, src, nlev, buf)
    use netcdf
    integer          , intent(in)    :: ncid, varid, nlev
    type(bounds_type), intent(in)    :: bounds
    real(r8)         , intent(in)    :: src(bounds%begc:, 1:)
    real(r8)         , intent(inout) :: buf(:,:)
    integer :: c, j, ier
    do c = bounds%begc, bounds%endc
       do j = 1, nlev
          buf(j, c - bounds%begc + 1) = src(c, j)
       end do
    end do
    ier = nf90_put_var(ncid, varid, buf)
  end subroutine put2d

end module bgcdumpMod
