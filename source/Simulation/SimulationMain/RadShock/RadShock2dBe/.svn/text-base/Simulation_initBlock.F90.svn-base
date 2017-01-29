!!****if* source/Simulation/SimulationMain/RadShock/RadShock2dBe/Simulation_initBlock
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


  real :: dens_bar

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
                ycent(j) >  sim_vacThickness + sim_slabThickness) then

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


!            if (xcent(i) <= 312.5e-4 .and. &
!                ycent(j) <= 100.0e-4 .and. &
!                ycent(j) >  100.0e-4 - sim_gradSize) then
!               ! ***************************
!               ! *     GRADIENT REGION     *
!               ! ***************************
!               dens_bar = sim_rhoVa + (ycent(j) - (100.0e-4-sim_gradSize))*(sim_rhoBe - sim_rhoVa) / sim_gradSize
!               call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, dens_bar)
!               call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, axis, sim_teleBe)
! #ifdef TION_VAR
!               call Grid_putPointData(blockId, CENTER, TION_VAR, EXTERIOR, axis, sim_tionBe)
! #endif
! #ifdef TELE_VAR
!               call Grid_putPointData(blockId, CENTER, TELE_VAR, EXTERIOR, axis, sim_teleBe)
! #endif
!               call setMaterial(blockId, axis, BE_SPEC)
!            endif

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
    real, parameter :: SMALL_FRAC = 1.e-10

    integer :: n

    if(NSPECIES == 0) return

    do n = SPECIES_BEGIN, SPECIES_END
       if(n == spec) then
          call Grid_putPointData(blockID, CENTER, n, EXTERIOR, axis, 1.0-real(NSPECIES-1)*SMALL_FRAC)
       else
          call Grid_putPointData(blockID, CENTER, n, EXTERIOR, axis, SMALL_FRAC)
       end if
    end do

  end subroutine setMaterial

end subroutine Simulation_initBlock
