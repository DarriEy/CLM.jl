module SoilBudgetDumpMod

  ! ======================================================================================
  ! CLM host top-soil-layer WATER BUDGET dump — companion to FatesParityDumpMod, for the
  ! PR #247 D1 "top layer over-dries" divergence (h2ovol1 J 0.11-0.14 vs F 0.32-0.34).
  !
  ! FatesParityDumpMod dumps only the top-layer water STATE that FATES sees (h2ovol1).
  ! This module dumps the per-column WATER BUDGET that SETS that state — the per-layer
  ! h2osoi, and every flux into/out of the soil column: infiltration, surface runoff,
  ! subsurface drainage (+ perched), soil evaporation, transpiration and per-layer root
  ! extraction, the SL14 dry-surface-layer resistance/thickness, and the water table.
  ! Emitted once per soil column per timestep, appended to soil_budget_fortran.txt.
  !
  ! Called from clm_driver.F90 at the FATES "fast" phase, which is AFTER HydrologyDrainage
  ! for the timestep, so every qflx_* is the value that acted this step.
  !
  !   SOILBUD t=<nstep> c=<col> hvol1 hvol2 hvol3 hliq1 hice1 tsoi1 watsat1 effporo1 dz1
  !           soilresis dsl zwt q_infl q_surf q_drain q_drainp q_ev_soil q_liqevtop
  !           q_tran_veg q_root1 q_root2 q_root3
  ! Reals are 17 significant digits (es-format), matching the Julia probe.
  ! ======================================================================================

  use shr_kind_mod          , only : r8 => shr_kind_r8
  use decompMod             , only : bounds_type
  use ColumnType            , only : col
  use WaterStateBulkType    , only : waterstatebulk_type
  use WaterFluxBulkType     , only : waterfluxbulk_type
  use SoilStateType         , only : soilstate_type
  use SoilHydrologyType     , only : soilhydrology_type
  use TemperatureType       , only : temperature_type

  implicit none
  private

  logical, public :: soilbud_opened = .false.
  character(len=*), parameter, public :: soilbud_path = 'soil_budget_fortran.txt'

  public :: soil_budget_dump

contains

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

  subroutine soil_budget_dump(bounds, num_soilc, filter_soilc, nstep, &
       waterstatebulk_inst, waterfluxbulk_inst, soilstate_inst, &
       soilhydrology_inst, temperature_inst)

    type(bounds_type)         , intent(in) :: bounds
    integer                   , intent(in) :: num_soilc
    integer                   , intent(in) :: filter_soilc(:)
    integer                   , intent(in) :: nstep
    type(waterstatebulk_type) , intent(in) :: waterstatebulk_inst
    type(waterfluxbulk_type)  , intent(in) :: waterfluxbulk_inst
    type(soilstate_type)      , intent(in) :: soilstate_inst
    type(soilhydrology_type)  , intent(in) :: soilhydrology_inst
    type(temperature_type)    , intent(in) :: temperature_inst

    integer :: iu, ios, fc, c

    associate( &
         h2osoi_vol => waterstatebulk_inst%h2osoi_vol_col , &
         h2osoi_liq => waterstatebulk_inst%h2osoi_liq_col , &
         h2osoi_ice => waterstatebulk_inst%h2osoi_ice_col , &
         t_soisno   => temperature_inst%t_soisno_col      , &
         watsat     => soilstate_inst%watsat_col          , &
         effporo    => soilstate_inst%eff_porosity_col    , &
         soilresis  => soilstate_inst%soilresis_col       , &
         dsl        => soilstate_inst%dsl_col             , &
         zwt        => soilhydrology_inst%zwt_col         , &
         q_infl     => waterfluxbulk_inst%qflx_infl_col   , &
         q_surf     => waterfluxbulk_inst%qflx_surf_col   , &
         q_drain    => waterfluxbulk_inst%qflx_drain_col  , &
         q_drainp   => waterfluxbulk_inst%qflx_drain_perched_col , &
         q_ev_soil  => waterfluxbulk_inst%qflx_ev_soil_col , &
         q_liqevtop => waterfluxbulk_inst%qflx_liqevap_from_top_layer_col , &
         q_tran_veg => waterfluxbulk_inst%qflx_tran_veg_col , &
         q_rootsoi  => waterfluxbulk_inst%qflx_rootsoi_col , &
         dz         => col%dz )

    if (.not. soilbud_opened) then
       open(newunit=iu, file=soilbud_path, status='replace', action='write', iostat=ios)
       if (ios /= 0) return
       write(iu,'(a)') '# CLM soil budget dump - generator=fortran'
       soilbud_opened = .true.
    else
       open(newunit=iu, file=soilbud_path, status='old', action='write', &
            position='append', iostat=ios)
       if (ios /= 0) return
    end if

    do fc = 1, num_soilc
       c = filter_soilc(fc)
       write(iu,'(a)') 'SOILBUD t='//trim(istr(nstep))//' c='//trim(istr(c)) &
            //' hvol1='//trim(rstr(h2osoi_vol(c,1))) &
            //' hvol2='//trim(rstr(h2osoi_vol(c,2))) &
            //' hvol3='//trim(rstr(h2osoi_vol(c,3))) &
            //' hliq1='//trim(rstr(h2osoi_liq(c,1))) &
            //' hice1='//trim(rstr(h2osoi_ice(c,1))) &
            //' tsoi1='//trim(rstr(t_soisno(c,1))) &
            //' watsat1='//trim(rstr(watsat(c,1))) &
            //' effporo1='//trim(rstr(effporo(c,1))) &
            //' dz1='//trim(rstr(dz(c,1))) &
            //' soilresis='//trim(rstr(soilresis(c))) &
            //' dsl='//trim(rstr(dsl(c))) &
            //' zwt='//trim(rstr(zwt(c))) &
            //' q_infl='//trim(rstr(q_infl(c))) &
            //' q_surf='//trim(rstr(q_surf(c))) &
            //' q_drain='//trim(rstr(q_drain(c))) &
            //' q_drainp='//trim(rstr(q_drainp(c))) &
            //' q_ev_soil='//trim(rstr(q_ev_soil(c))) &
            //' q_liqevtop='//trim(rstr(q_liqevtop(c))) &
            //' q_tran_veg='//trim(rstr(q_tran_veg(c))) &
            //' q_root1='//trim(rstr(q_rootsoi(c,1))) &
            //' q_root2='//trim(rstr(q_rootsoi(c,2))) &
            //' q_root3='//trim(rstr(q_rootsoi(c,3)))
    end do
    close(iu)

    end associate

  end subroutine soil_budget_dump

end module SoilBudgetDumpMod
