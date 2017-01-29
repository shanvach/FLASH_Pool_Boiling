!!****if* source/Simulation/SimulationMain/RadShock/RadShock2dFull/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(IN) :: blockID) 
!!                       
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.
!! 
!! ARGUMENTS
!!
!!  blockID -        the number of the block to initialize
!!  
!!
!!
!!***
subroutine Simulation_initBlock(blockId)
  use Simulation_data
  use Eos_interface, ONLY : Eos, Eos_wrapped
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
       Grid_getCellCoords, Grid_putPointData
  use Driver_interface, ONLY: Driver_abortFlash
  use PhysicalConstants_interface, ONLY: PhysicalConstants_get
  use RadTrans_interface, ONLY: RadTrans_mgdEFromT

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"   

  ! compute the maximum length of a vector in each coordinate direction 
  ! (including guardcells)

  integer, intent(in) :: blockId
  
  integer :: i, j, k, n
  integer :: blkLimits(2, MDIM)
  integer :: blkLimitsGC(2, MDIM)
  integer :: axis(MDIM)
  real, allocatable :: xcent(:), ycent(:)

  ! get the coordinate information for the current block from the database
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)

  ! get the coordinate information for the current block from the database
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  allocate(xcent(blkLimitsGC(HIGH, IAXIS)))
  call Grid_getCellCoords(IAXIS, blockId, CENTER, .true., &
       xcent, blkLimitsGC(HIGH, IAXIS))

  allocate(ycent(blkLimitsGC(HIGH, JAXIS)))
  call Grid_getCellCoords(JAXIS, blockId, CENTER, .true., &
       ycent, blkLimitsGC(HIGH, JAXIS))
  !------------------------------------------------------------------------------

  ! Loop over cells and set the initial state
  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

           axis(IAXIS) = i
           axis(JAXIS) = j
           axis(KAXIS) = k

           call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, 0.0)
           
           if(  xcent(i) <= sim_tubeRadius .and. &
                ycent(j) >  sim_vacThickness + sim_slabThickness) then 
              ! ********************
              ! *   XENON REGION   *
              ! ********************
              call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, sim_rhoXe)
              call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, axis, sim_teleXe)
              call Grid_putPointData(blockId, CENTER, TION_VAR, EXTERIOR, axis, sim_tionXe)
              call Grid_putPointData(blockId, CENTER, TELE_VAR, EXTERIOR, axis, sim_teleXe)
              call Grid_putPointData(blockId, CENTER, TRAD_VAR, EXTERIOR, axis, sim_tradXe)
              call setMaterial(blockId, axis, XE_SPEC)
              call RadTrans_mgdEFromT(blockId, axis, sim_tradXe)
           elseif( xcent(i) <= sim_tubeRadius + sim_tubeThickness .and. &
                xcent(i) >  sim_tubeRadius .and. &
                ycent(j) >  sim_vacThickness + sim_slabThickness + &
                sim_goldThickness + sim_windowThickness) then

              ! **********************
              ! *   PLASTIC REGION   *
              ! **********************
              call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, sim_rhoCh)
              call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, axis, sim_teleCh)
              call Grid_putPointData(blockId, CENTER, TION_VAR, EXTERIOR, axis, sim_tionCh)
              call Grid_putPointData(blockId, CENTER, TELE_VAR, EXTERIOR, axis, sim_teleCh)
              call Grid_putPointData(blockId, CENTER, TRAD_VAR, EXTERIOR, axis, sim_tradCh)
              call setMaterial(blockId, axis, POLI_SPEC)
              call RadTrans_mgdEFromT(blockId, axis, sim_tradCh)

           elseif(xcent(i) > sim_tubeRadius .and. &
                xcent(i) <= sim_goldRadius .and. &
                ycent(j) > sim_vacThickness + sim_slabThickness .and. &
                ycent(j) <= sim_vacThickness + sim_slabThickness + sim_goldThickness) then

              ! *******************
              ! *   GOLD REGION   *
              ! *******************
              call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, sim_rhoAu)
              call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, axis, sim_teleAu)
              call Grid_putPointData(blockId, CENTER, TION_VAR, EXTERIOR, axis, sim_tionAu)
              call Grid_putPointData(blockId, CENTER, TELE_VAR, EXTERIOR, axis, sim_teleAu)
              call Grid_putPointData(blockId, CENTER, TRAD_VAR, EXTERIOR, axis, sim_tradAu)
              call setMaterial(blockId, axis, GOLD_SPEC)
              call RadTrans_mgdEFromT(blockId, axis, sim_tradAu)

           elseif( ( &
                xcent(i) > sim_tubeRadius .and. &
                xcent(i) <= sim_acrylicRadius .and. &
                ycent(j) > sim_vacThickness + sim_slabThickness + sim_goldThickness .and. &
                ycent(j) <= sim_vacThickness + sim_slabThickness + &
                sim_goldThickness + sim_windowThickness ) .or. ( &
                
                xcent(i) > sim_tubeRadius + sim_tubeThickness .and. &
                xcent(i) <= sim_acrylicRadius .and. &
                ycent(j) > sim_vacThickness + sim_slabThickness + sim_goldThickness + &
                sim_windowThickness .and. &
                ycent(j) <= sim_vacThickness + sim_slabThickness + sim_goldThickness + &
                sim_acrylicThickness ) ) then

              ! **********************
              ! *   ACRYLIC REGION   *
              ! **********************
              call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, sim_rhoAc)
              call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, axis, sim_teleAc)
              call Grid_putPointData(blockId, CENTER, TION_VAR, EXTERIOR, axis, sim_tionAc)
              call Grid_putPointData(blockId, CENTER, TELE_VAR, EXTERIOR, axis, sim_teleAc)
              call Grid_putPointData(blockId, CENTER, TRAD_VAR, EXTERIOR, axis, sim_tradAc)
              call setMaterial(blockId, axis, ACRY_SPEC)
              call RadTrans_mgdEFromT(blockId, axis, sim_tradAc)


           elseif( xcent(i) <= sim_slabRadius .and. &
              ycent(j) <= sim_vacThickness + sim_slabThickness .and. &
              ycent(j) >  sim_vacThickness) then
              ! ************************
              ! *   BERYLLIUM REGION   *
              ! ************************
              call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, sim_rhoBe)
              call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, axis, sim_teleBe)
              call Grid_putPointData(blockId, CENTER, TION_VAR, EXTERIOR, axis, sim_tionBe)
              call Grid_putPointData(blockId, CENTER, TELE_VAR, EXTERIOR, axis, sim_teleBe)
              call Grid_putPointData(blockId, CENTER, TRAD_VAR, EXTERIOR, axis, sim_tradBe)
              call setMaterial(blockId, axis, BE_SPEC)
              call RadTrans_mgdEFromT(blockId, axis, sim_tradBe)
           else
              ! *********************
              ! *   VACUUM REGION   *
              ! *********************
              call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, sim_rhoVa)
              call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, axis, sim_teleVa)
              call Grid_putPointData(blockId, CENTER, TION_VAR, EXTERIOR, axis, sim_tionVa)
              call Grid_putPointData(blockId, CENTER, TELE_VAR, EXTERIOR, axis, sim_teleVa)
              call Grid_putPointData(blockId, CENTER, TRAD_VAR, EXTERIOR, axis, sim_tradVa)
              call setMaterial(blockId, axis, VACU_SPEC)
              call RadTrans_mgdEFromT(blockId, axis, sim_tradVa)
           end if

        enddo
     enddo
  enddo

  deallocate(xcent)
  deallocate(ycent)

  return

contains

  subroutine setMaterial(blockID, axis, spec)
    implicit none

    ! The species to fill the cell with:
    integer, intent(in) :: spec 
    integer, intent(in) :: blockID
    integer, intent(in) :: axis(MDIM)

    ! It is not good to set the mass fraction of any species to
    ! zero. The parameter SMALL_FREQ is the floor value set on the
    ! mass fraction of any species.

    integer :: n

    if(NSPECIES == 0) return

    do n = SPECIES_BEGIN, SPECIES_END
       if(n == spec) then
          call Grid_putPointData(blockID, CENTER, n, EXTERIOR, axis, 1.0-real(NSPECIES-1)*sim_smallX)
       else
          call Grid_putPointData(blockID, CENTER, n, EXTERIOR, axis, sim_smallX)
       end if
    end do

  end subroutine setMaterial

end subroutine Simulation_initBlock
