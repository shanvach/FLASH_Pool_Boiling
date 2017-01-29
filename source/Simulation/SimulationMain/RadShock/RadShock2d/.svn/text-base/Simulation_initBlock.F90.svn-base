!!****if* source/Simulation/SimulationMain/RadShock/RadShock2d/Simulation_initBlock
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

           call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, sim_vely)

           ! The -1 factor in the xcent index is used to put a buffer
           ! region 1 cell thick
           if(xcent(i-sim_nbuffer) <= 287.5e-04) then 
              ! ********************
              ! *   XENON REGION   *
              ! ********************
              call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, sim_rhoXe)
              call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, axis, sim_teleXe)
#ifdef TION_VAR
              call Grid_putPointData(blockId, CENTER, TION_VAR, EXTERIOR, axis, sim_tionXe)
#endif
#ifdef TELE_VAR
              call Grid_putPointData(blockId, CENTER, TELE_VAR, EXTERIOR, axis, sim_teleXe)
#endif
              call RadTrans_mgdEFromT(blockId, axis, sim_tradXe)
              call Grid_putPointData(blockId, CENTER, TRAD_VAR, EXTERIOR, axis, sim_tradXe)
              call setMaterial(blockId, axis, XE_SPEC)

           elseif(xcent(i) <= 312.5e-04) then

              ! **********************
              ! *   PLASTIC REGION   *
              ! **********************
              call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, sim_rhoCh)
              call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, axis, sim_teleCh)
#ifdef TION_VAR
              call Grid_putPointData(blockId, CENTER, TION_VAR, EXTERIOR, axis, sim_tionCh)
#endif
#ifdef TELE_VAR
              call Grid_putPointData(blockId, CENTER, TELE_VAR, EXTERIOR, axis, sim_teleCh)
#endif
              call RadTrans_mgdEFromT(blockId, axis, sim_tradCh)
              call Grid_putPointData(blockId, CENTER, TRAD_VAR, EXTERIOR, axis, sim_tradCh)
              call setMaterial(blockId, axis, POLI_SPEC)

           else

              ! *********************
              ! *   VACUUM REGION   *
              ! *********************
              call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, sim_rhoVa)
              call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, axis, sim_teleVa)
#ifdef TION_VAR
              call Grid_putPointData(blockId, CENTER, TION_VAR, EXTERIOR, axis, sim_tionVa)
#endif
#ifdef TELE_VAR
              call Grid_putPointData(blockId, CENTER, TELE_VAR, EXTERIOR, axis, sim_teleVa)
#endif
              call RadTrans_mgdEFromT(blockId, axis, sim_tradVa)
              call Grid_putPointData(blockId, CENTER, TRAD_VAR, EXTERIOR, axis, sim_tradVa)
              call setMaterial(blockId, axis, VACU_SPEC)
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
    real, parameter :: SMALL_FRAC = 1.e-10

    integer :: n

    if(NSPECIES == 0) return

    do n = SPECIES_BEGIN, SPECIES_END
       if(n == spec) then
          call Grid_putPointData(blockID, CENTER, n, EXTERIOR, axis, 1.0-(NSPECIES-1)*SMALL_FRAC)
       else
          call Grid_putPointData(blockID, CENTER, n, EXTERIOR, axis, SMALL_FRAC)
       end if
    end do

  end subroutine setMaterial

end subroutine Simulation_initBlock
