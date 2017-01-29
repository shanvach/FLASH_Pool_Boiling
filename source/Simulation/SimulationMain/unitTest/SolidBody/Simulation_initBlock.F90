!!****if* source/Simulation/SimulationMain/unitTest/SolidBody/Simulation_initBlock
!!
!! NAME
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!  Simulation_initBlock(integer(in) :: blockID)
!!                         
!!
!! DESCRIPTION   
!!
!! ARGUMENTS
!!      blockID:     integer(in)      the current block number to be filled
!!      
!!
!! PARAMETERS
!!
!!***


 subroutine Simulation_initBlock(blockID)

!============================================================================
  use Simulation_data, ONLY : sim_meshMe
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr
  implicit none 

#include "constants.h"
#include "Flash.h"


  integer, intent(IN) :: blockID
  real, pointer, dimension(:,:,:,:) :: solnData

  call Grid_getBlkPtr(blockID, solnData)
  solnData = sim_meshMe
  call Grid_releaseBlkPtr(blockID,solnData)

end subroutine Simulation_initBlock
