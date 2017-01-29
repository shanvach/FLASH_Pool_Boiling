!!****if* source/Simulation/SimulationMain/SuOlson/Simulation_initBlock
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
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
       Grid_getCellCoords, Grid_putPointData, Grid_getPointData
  use Driver_interface, ONLY: Driver_abortFlash
  use RadTrans_interface, ONLY: RadTrans_mgdEFromT
  use Eos_interface, ONLY : Eos_wrapped

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"    

  ! compute the maximum length of a vector in each coordinate direction 
  ! (including guardcells)

  integer, intent(in) :: blockId
  
  integer :: i, j, k
  integer :: blkLimits(2, MDIM)
  integer :: blkLimitsGC(2, MDIM)
  integer :: axis(MDIM)
  real, allocatable :: xcent(:), ycent(:), zcent(:)

  real :: rho   ! Density
  real :: tion  ! Ion temperature
  real :: tele  ! Electron temperature 
  real :: trad  ! Radiation temperature
  real :: tradActual
  real :: velx

  
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  
  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)                     
           
           rho  = 1.0
           tele = 1.0E-10
           tion = 1.0E-10
           trad = 1.0E-10
           velx = 0.0
           
           axis(IAXIS) = i
           axis(JAXIS) = j
           axis(KAXIS) = k
           
           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rho)
           call Grid_putPointData(blockId, CENTER, TION_VAR, EXTERIOR, axis, tion)
           call Grid_putPointData(blockId, CENTER, TELE_VAR, EXTERIOR, axis, tele)
           call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, velx)           
           
           call RadTrans_mgdEFromT(blockId, axis, trad, tradActual)
           
           call Grid_putPointData(blockId, CENTER, TRAD_VAR, EXTERIOR, axis, tradActual)           
           
        enddo
     enddo
  enddo
    
  call Eos_wrapped(MODE_DENS_TEMP_GATHER, blkLimits, blockID)
  
  
  return

end subroutine Simulation_initBlock
