!!****if* source/Simulation/SimulationMain/INavierStokes/2D/LidDrivenCavity/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(in) :: blockID) 
!!                       
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.
!!
!!  Reference:
!!
!! 
!! ARGUMENTS
!!
!!  blockID -          the number of the block to update
!!  myPE   -           my processor number
!!
!! 
!!
!!***

subroutine Simulation_initBlock(blockId)
  use Grid_interface, only  :  Grid_getBlkPtr, Grid_releaseBlkPtr
  use Driver_data, only : dr_meshMe

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent(in) :: blockID
 
  real, pointer, dimension(:,:,:,:) :: solnData

  call Grid_getBlkPtr(blockID, solnData, CENTER)
  solnData(NSRC_VAR,:,:,:) = 1.0
  solnData(NSRC_VAR,:,:,:) = real(dr_meshMe)

  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  return
 
end subroutine Simulation_initBlock
